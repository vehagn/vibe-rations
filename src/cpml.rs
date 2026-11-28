//! Convolutional perfectly matched layers (Komatitsch & Martin 2007, Geophysics 72(5)).
//!
//! Inside a layer, each derivative ∂ normal to the boundary becomes ∂/κ + ψ, where the memory
//! variable follows ψ ← b·ψ + a·∂ every step. Memory is kept only for the nodes in the layers
//! (the "slab"): the first `lo` and the last `hi` nodes along the axis, in that order.

use crate::grid::Grid;
use crate::simd::Arithmetic;

/// How the layers damp.
#[derive(Clone, Copy, Debug)]
pub struct PmlSettings {
    /// Exponent of the damping profile.
    pub power: f64,
    /// Theoretical reflection coefficient; by default from the layer's width.
    pub reflection: Option<f64>,
    /// κ at the outer edge (1 means no coordinate stretching).
    pub kappa_max: f64,
    /// α at the inner edge, in 1/s; π times the dominant frequency is usual.
    pub alpha_max: f64,
    /// The fastest P-wave speed in the model, which sets the damping.
    pub vp_max: f64,
    /// Plain damping of every field inside the side layers when there is a free surface, as a
    /// fraction of their damping profile. CPMLs alone can grow without bound where surface
    /// waves in a medium with a high Poisson's ratio run into them along the free surface
    /// (Festa et al. 2005); damping that doesn't depend on direction, as in the multiaxial PML
    /// (Meza-Fajardo & Papageorgiou 2008), keeps them stable.
    pub sponge: f64,
}

/// The absorbing layers' thickness at the start and end of each axis of `grid`, for layers
/// `width` nodes thick: none on top with a free surface, and none across a 2D grid, which
/// has no sides along y.
pub fn layer_widths(grid: &Grid, width: usize, free_surface: bool) -> [(usize, usize); 3] {
    std::array::from_fn(|a| {
        if !grid.spans(a) {
            (0, 0)
        } else if a == 2 && free_surface {
            (0, width)
        } else {
            (width, width)
        }
    })
}

/// The reflection coefficient Collino & Tsogka (2001) suggest for a layer `width` nodes
/// thick: 1e-2 for 5 nodes, 1e-3 for 10, 1e-4 for 20.
pub fn default_reflection(width: usize) -> f64 {
    let w = width.max(2) as f64;
    10f64.powf(-((w.log10() - 1.0) / 2f64.log10() + 3.0))
}

/// Coefficients for each node in a slab.
#[derive(Clone, Debug, Default)]
pub struct Profile {
    pub a: Vec<f32>,
    pub b: Vec<f32>,
    pub inv_kappa: Vec<f32>,
    /// The damping d itself, in 1/s.
    pub d: Vec<f32>,
    /// Depth into the layer, from 0 at its inner edge to 1 at its outer edge.
    pub x: Vec<f32>,
}

/// Where along the staggered grid a derivative sits.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Stagger {
    /// At the nodes, i.
    Node,
    /// Half a node further along the axis, i + 1/2.
    Half,
}

/// The layers at both ends of one axis.
#[derive(Clone, Debug)]
pub struct CpmlAxis {
    n: usize,
    lo: usize,
    hi: usize,
    node: Profile,
    half: Profile,
}

impl CpmlAxis {
    /// Layers `lo` and `hi` nodes thick at the start and end of an axis of `n` nodes spaced `h`.
    ///
    /// The inner edge of each layer sits half a node outside the last node that is not in it,
    /// so that both the node and the half-node profiles are mirror images at the two ends.
    pub fn new(n: usize, lo: usize, hi: usize, h: f64, dt: f64, settings: &PmlSettings) -> Self {
        assert!(lo + hi <= n, "layers wider than the axis");
        let mut node = Profile::default();
        let mut half = Profile::default();
        for slab in 0..lo + hi {
            let (width, at_node, at_half) = if slab < lo {
                let p = slab as f64;
                (lo, lo as f64 - 0.5 - p, lo as f64 - 1.0 - p)
            } else {
                let p = (n - hi + slab - lo) as f64;
                let edge = (n - hi) as f64;
                (hi, p - edge + 0.5, p - edge + 1.0)
            };
            let r = settings
                .reflection
                .unwrap_or_else(|| default_reflection(width));
            let thickness = width as f64 * h;
            let d0 = -(settings.power + 1.0) * settings.vp_max * r.ln() / (2.0 * thickness);
            push(&mut node, at_node / width as f64, d0, dt, settings);
            push(&mut half, at_half / width as f64, d0, dt, settings);
        }
        Self {
            n,
            lo,
            hi,
            node,
            half,
        }
    }

    /// Nodes in the slab.
    pub fn width(&self) -> usize {
        self.lo + self.hi
    }

    /// The nodes along the axis that are not in a layer.
    pub fn interior(&self) -> std::ops::Range<usize> {
        self.lo..self.n - self.hi
    }

    /// The nodes along the axis that are in the layers: those at the start, then the end.
    pub fn ranges(&self) -> [std::ops::Range<usize>; 2] {
        [0..self.lo, self.n - self.hi..self.n]
    }

    /// exp(-sponge · d · x³ · dt) at each node along the axis, x being the depth into the
    /// layer from 0 to 1: 1 outside the layers.
    pub fn sponge_factors(&self, sponge: f64, dt: f64) -> Vec<f32> {
        (0..self.n)
            .map(|p| match self.slab(p) {
                Some(s) => {
                    // d x³ rises later than d itself, which reflects less.
                    let x = self.node.x[s] as f64;
                    (-sponge * self.node.d[s] as f64 * x.powi(3) * dt).exp() as f32
                }
                None => 1.0,
            })
            .collect()
    }

    /// The slab index of node `p` along the axis, or None if it is not in a layer.
    #[inline(always)]
    pub fn slab(&self, p: usize) -> Option<usize> {
        if p < self.lo {
            Some(p)
        } else if p >= self.n - self.hi {
            Some(self.lo + p - (self.n - self.hi))
        } else {
            None
        }
    }

    #[inline(always)]
    pub fn profile(&self, at: Stagger) -> &Profile {
        match at {
            Stagger::Node => &self.node,
            Stagger::Half => &self.half,
        }
    }

    /// The slab indices of the nodes in `nodes` that are in a layer: contiguous, since nodes
    /// that reach into both layers hold every node between them.
    pub fn slabs(&self, nodes: std::ops::Range<usize>) -> std::ops::Range<usize> {
        let start = self.n - self.hi;
        let low = nodes.start.min(self.lo)..nodes.end.min(self.lo);
        let high = self.lo + nodes.start.max(start) - start..self.lo + nodes.end.max(start) - start;
        match (low.is_empty(), high.is_empty()) {
            (true, true) => 0..0,
            (false, true) => low,
            (true, false) => high,
            (false, false) => low.start..high.end,
        }
    }

    /// Damps derivatives taken along this axis at the consecutive nodes `nodes` of a row: `d`
    /// holds them, and `psi` the memory of those in a layer ([`CpmlAxis::slabs`]). Always
    /// inlined, like the rest of the update loops (crate::simd).
    #[inline(always)]
    pub fn apply_along<A: Arithmetic>(
        &self,
        nodes: std::ops::Range<usize>,
        d: &mut [f32],
        psi: &mut [f32],
        at: Stagger,
    ) {
        let p = self.profile(at);
        let (lo, start) = (self.lo, self.n - self.hi);
        let first = self.slabs(nodes.clone()).start;
        // Nodes, as offsets into `d`, and their slab indices, in each layer.
        let low = (nodes.start.min(lo), nodes.end.min(lo));
        let high = (nodes.start.max(start), nodes.end.max(start));
        for ((a, b), slab) in [(low, low.0), (high, lo + high.0 - start)] {
            if b > a {
                let s = slab..slab + (b - a);
                update::<A>(
                    &mut d[a - nodes.start..b - nodes.start],
                    &mut psi[s.start - first..s.end - first],
                    &p.a[s.clone()],
                    &p.b[s.clone()],
                    &p.inv_kappa[s],
                );
            }
        }
    }

    /// Damps a row of derivatives taken along this axis while the row runs across it, so the
    /// whole row sits at slab index `slab`.
    #[inline(always)]
    pub fn apply_across<A: Arithmetic>(
        &self,
        slab: usize,
        d: &mut [f32],
        psi: &mut [f32],
        at: Stagger,
    ) {
        let p = self.profile(at);
        let (a, b, k) = (p.a[slab], p.b[slab], p.inv_kappa[slab]);
        for (d, psi) in d.iter_mut().zip(psi) {
            *psi = A::product_plus(b, *psi, a * *d);
            *d = A::product_plus(*d, k, *psi);
        }
    }
}

#[inline(always)]
fn update<A: Arithmetic>(d: &mut [f32], psi: &mut [f32], a: &[f32], b: &[f32], inv_kappa: &[f32]) {
    for ((((d, psi), a), b), k) in d.iter_mut().zip(psi).zip(a).zip(b).zip(inv_kappa) {
        *psi = A::product_plus(*b, *psi, a * *d);
        *d = A::product_plus(*d, *k, *psi);
    }
}

/// Appends the coefficients at normalised depth `x` into the layer (0 at its inner edge).
fn push(p: &mut Profile, x: f64, d0: f64, dt: f64, s: &PmlSettings) {
    let x = x.clamp(0.0, 1.0);
    let xp = x.powf(s.power);
    let d = d0 * xp;
    let kappa = 1.0 + (s.kappa_max - 1.0) * xp;
    let alpha = s.alpha_max * (1.0 - x);
    let b = (-(d / kappa + alpha) * dt).exp();
    let a = if d > 0.0 {
        d * (b - 1.0) / (kappa * (d + kappa * alpha))
    } else {
        0.0
    };
    p.a.push(a as f32);
    p.b.push(b as f32);
    p.inv_kappa.push((1.0 / kappa) as f32);
    p.d.push(d as f32);
    p.x.push(x as f32);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simd::Separate;

    fn settings() -> PmlSettings {
        PmlSettings {
            power: 2.0,
            reflection: None,
            kappa_max: 1.0,
            alpha_max: std::f64::consts::PI * 10.0,
            vp_max: 3000.0,
            sponge: 0.0,
        }
    }

    #[test]
    fn reflection_follows_collino_and_tsogka() {
        assert!((default_reflection(10) - 1e-3).abs() < 1e-12);
        assert!((default_reflection(20) - 1e-4).abs() < 1e-12);
    }

    #[test]
    fn maps_nodes_to_the_slab() {
        let c = CpmlAxis::new(50, 10, 8, 10.0, 1e-3, &settings());
        assert_eq!(c.width(), 18);
        assert_eq!(c.slab(0), Some(0));
        assert_eq!(c.slab(9), Some(9));
        assert_eq!(c.slab(10), None);
        assert_eq!(c.slab(41), None);
        assert_eq!(c.slab(42), Some(10));
        assert_eq!(c.slab(49), Some(17));
    }

    #[test]
    fn profiles_mirror_at_both_ends() {
        let c = CpmlAxis::new(40, 10, 10, 10.0, 1e-3, &settings());
        let (node, half) = (c.profile(Stagger::Node), c.profile(Stagger::Half));
        for s in 0..10 {
            // Node i mirrors node n-1-i; half-node i+1/2 mirrors half-node n-2-i.
            assert_eq!(node.a[s], node.a[19 - s]);
            if s < 9 {
                assert_eq!(half.a[s], half.a[18 - s]);
            }
        }
        // Damping grows outwards, and the half-node just inside the layer is undamped.
        assert!(node.a[0] < node.a[9] && node.a[9] < 0.0);
        assert_eq!(half.a[9], 0.0);
    }

    #[test]
    fn an_axis_without_layers_does_nothing() {
        let c = CpmlAxis::new(10, 0, 0, 1.0, 1e-3, &settings());
        assert_eq!(c.width(), 0);
        assert_eq!(c.slab(0), None);
        let mut d = [1.0f32; 10];
        assert_eq!(c.slabs(0..10), 0..0);
        c.apply_along::<Separate>(0..10, &mut d, &mut [], Stagger::Node);
        assert_eq!(d, [1.0; 10]);
    }

    /// Damping a row in parts, each with its part of the memory, is damping it whole, bit
    /// for bit, wherever the parts split it.
    #[test]
    fn damps_a_row_in_parts_as_a_whole() {
        let c = CpmlAxis::new(40, 6, 5, 10.0, 1e-3, &settings());
        assert_eq!(c.slabs(0..40), 0..11);
        assert_eq!(c.slabs(3..20), 3..6);
        assert_eq!(c.slabs(6..35), 0..0);
        assert_eq!(c.slabs(30..38), 6..9);
        assert_eq!(c.slabs(4..37), 4..8);
        let row: Vec<f32> = (0..40).map(|i| (i as f32 * 0.37).sin()).collect();
        let memory: Vec<f32> = (0..11).map(|s| (s as f32 * 1.3).cos()).collect();
        let (mut whole, mut psi) = (row.clone(), memory.clone());
        c.apply_along::<Separate>(0..40, &mut whole, &mut psi, Stagger::Half);
        for cuts in [
            vec![0, 40],
            vec![0, 3, 7, 33, 36, 40],
            vec![0, 6, 34, 40],
            vec![0, 1, 39, 40],
        ] {
            let (mut parts, mut psi_parts) = (row.clone(), memory.clone());
            for w in cuts.windows(2) {
                let (nodes, slabs) = (w[0]..w[1], c.slabs(w[0]..w[1]));
                c.apply_along::<Separate>(
                    nodes.clone(),
                    &mut parts[nodes],
                    &mut psi_parts[slabs],
                    Stagger::Half,
                );
            }
            assert_eq!(parts, whole, "{cuts:?}");
            assert_eq!(psi_parts, psi, "{cuts:?}");
        }
    }

    #[test]
    fn sponge_damps_only_inside_the_layers_and_more_outwards() {
        let c = CpmlAxis::new(40, 10, 10, 10.0, 1e-3, &settings());
        let f = c.sponge_factors(0.2, 1e-3);
        assert!(f[10..30].iter().all(|&v| v == 1.0));
        // It rises gently: just inside, the factor rounds to 1 in f32.
        assert!(f[0] < f[5] && f[5] < f[9] && f[9] <= 1.0);
        assert_eq!(f[0], f[39]);
        assert!(c.sponge_factors(0.0, 1e-3).iter().all(|&v| v == 1.0));
    }

    /// For a constant derivative the memory variable settles at a d / (1 - b), so the
    /// damped derivative tends to d (1/κ + a / (1 - b)).
    #[test]
    fn memory_variables_follow_their_recurrence() {
        let c = CpmlAxis::new(40, 10, 10, 10.0, 1e-3, &settings());
        let p = c.profile(Stagger::Node);
        let (a, b, k) = (p.a[3], p.b[3], p.inv_kappa[3]);
        let mut psi = [0.0f32; 2];
        let mut d = [1.0f32; 2];
        c.apply_across::<Separate>(3, &mut d, &mut psi, Stagger::Node);
        assert!((psi[0] - a).abs() < 1e-7 && (d[0] - (k + a)).abs() < 1e-6);
        for _ in 0..20_000 {
            d = [1.0; 2];
            c.apply_across::<Separate>(3, &mut d, &mut psi, Stagger::Node);
        }
        let settled = a / (1.0 - b);
        assert!(
            (psi[1] - settled).abs() < 1e-4 * settled.abs(),
            "{} vs {settled}",
            psi[1]
        );
        assert!((d[1] - (k + settled)).abs() < 1e-4);
    }
}
