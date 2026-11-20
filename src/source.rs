//! Point sources: source time functions and how sources enter the staggered grid.

use std::f64::consts::PI;

use crate::solver::{Field, Fields, PointOp, Solver};
use crate::tdt::Warp;

/// A source time function.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Wavelet {
    /// (1 - 2a) e^-a with a = (π f0 (t - t0))²: zero mean, peak 1 at t0, peak frequency f0.
    Ricker { f0: f64, t0: f64 },
    /// √π f0 e^-a with a = (π f0 (t - t0))²: unit area, centred on t0.
    Gaussian { f0: f64, t0: f64 },
}

impl Wavelet {
    pub fn value(&self, t: f64) -> f64 {
        match *self {
            Wavelet::Ricker { f0, t0 } => {
                let a = (PI * f0 * (t - t0)).powi(2);
                (1.0 - 2.0 * a) * (-a).exp()
            }
            Wavelet::Gaussian { f0, t0 } => {
                let a = (PI * f0 * (t - t0)).powi(2);
                PI.sqrt() * f0 * (-a).exp()
            }
        }
    }

    /// The time derivative.
    pub fn derivative(&self, t: f64) -> f64 {
        match *self {
            Wavelet::Ricker { f0, t0 } => {
                let c = (PI * f0).powi(2);
                let tau = t - t0;
                let a = c * tau * tau;
                2.0 * c * tau * (-a).exp() * (2.0 * a - 3.0)
            }
            Wavelet::Gaussian { f0, t0 } => {
                let c = (PI * f0).powi(2);
                let tau = t - t0;
                -2.0 * c * tau * PI.sqrt() * f0 * (-c * tau * tau).exp()
            }
        }
    }

    pub fn dominant_frequency(&self) -> f64 {
        match *self {
            Wavelet::Ricker { f0, .. } | Wavelet::Gaussian { f0, .. } => f0,
        }
    }

    /// The frequency above which the spectrum is negligible (a few percent of its peak or less),
    /// which sets the shortest wavelength the grid must resolve.
    pub fn max_frequency(&self) -> f64 {
        2.5 * self.dominant_frequency()
    }
}

/// The moment tensor (xx, yy, zz, xy, xz, yz) of a shear dislocation with the given strike,
/// dip and rake in degrees and scalar moment `m0`, in Aki & Richards' convention: x north,
/// y east, z down.
pub fn double_couple(strike: f64, dip: f64, rake: f64, m0: f64) -> [f64; 6] {
    let (phi, delta, lambda) = (strike.to_radians(), dip.to_radians(), rake.to_radians());
    let (sd, cd) = delta.sin_cos();
    let (s2d, c2d) = (2.0 * delta).sin_cos();
    let (sl, cl) = lambda.sin_cos();
    let (sp, cp) = phi.sin_cos();
    let (s2p, c2p) = (2.0 * phi).sin_cos();
    [
        -m0 * (sd * cl * s2p + s2d * sl * sp * sp),
        m0 * (sd * cl * s2p - s2d * sl * cp * cp),
        m0 * s2d * sl,
        m0 * (sd * cl * c2p + 0.5 * s2d * sl * s2p),
        -m0 * (cd * cl * cp + c2d * sl * sp),
        -m0 * (cd * cl * sp - c2d * sl * cp),
    ]
}

/// What a source does.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Mechanism {
    /// A moment tensor (xx, yy, zz, xy, xz, yz), scaling the wavelet as moment rate.
    Moment([f64; 6]),
    /// A force (x, y, z) in N, scaling the wavelet.
    Force([f64; 3]),
}

/// A point source on a node, ready to add to the wavefield.
#[derive(Clone, Debug)]
pub struct Source {
    pub node: [usize; 3],
    pub mechanism: Mechanism,
    pub wavelet: Wavelet,
    /// How it enters the wavefield, per unit of the wavelet.
    pub op: PointOp,
    /// The wavelet's values to inject at each step, pre-warped to remove leapfrog's time error
    /// ([`Source::remove_time_dispersion`]); without them, the wavelet itself.
    pub series: Option<Vec<f64>>,
}

impl Source {
    /// A moment tensor enters as σij -= dt·Ṁij/V and a force as vi += dt·Fi/(ρV), spread onto
    /// the staggered positions around the node with the receivers' interpolation weights
    /// ([`Solver::stagger_weights`]), so that a source is the transpose of a receiver at the
    /// same node. Shear components spread along both of their axes.
    pub fn new(solver: &Solver, node: [usize; 3], mechanism: Mechanism, wavelet: Wavelet) -> Self {
        let g = solver.grid;
        let volume = g.cell_volume();
        let idx = |n: [usize; 3]| g.idx(n[0], n[1], n[2]);
        // The positions around the node along each (axis, mirror parity), with their weights.
        let spread = |axes: &[(usize, f32)]| {
            let mut out = vec![(node, 1.0f32)];
            for &(axis, parity) in axes {
                out = out
                    .into_iter()
                    .flat_map(|(n, w)| {
                        solver
                            .stagger_weights(axis, n[axis], parity)
                            .into_iter()
                            .map(move |(m, v)| {
                                let mut n = n;
                                n[axis] = m;
                                (n, w * v)
                            })
                    })
                    .collect();
            }
            out
        };
        let mut terms = Vec::new();
        match mechanism {
            Mechanism::Moment(m) => {
                let scale = |v: f64| (-solver.dt * v / volume) as f32;
                for (field, v) in [(Field::Sxx, m[0]), (Field::Syy, m[1]), (Field::Szz, m[2])] {
                    if v != 0.0 {
                        terms.push((field, idx(node), scale(v)));
                    }
                }
                // σxz and σyz are odd about the free surface.
                for (field, v, axes) in [
                    (Field::Sxy, m[3], [(0, 1.0), (1, 1.0)]),
                    (Field::Sxz, m[4], [(0, 1.0), (2, -1.0)]),
                    (Field::Syz, m[5], [(1, 1.0), (2, -1.0)]),
                ] {
                    if v != 0.0 {
                        let s = scale(v);
                        terms.extend(
                            spread(&axes)
                                .into_iter()
                                .map(|(n, w)| (field, idx(n), s * w)),
                        );
                    }
                }
            }
            Mechanism::Force(f) => {
                for (axis, field) in [(0, Field::Vx), (1, Field::Vy), (2, Field::Vz)] {
                    if f[axis] != 0.0 {
                        let s = (f[axis] / volume) as f32;
                        let b = solver.material.buoyancy(axis);
                        terms.extend(
                            spread(&[(axis, 1.0)])
                                .into_iter()
                                .map(|(n, w)| (field, idx(n), s * w * b[idx(n)])),
                        );
                    }
                }
            }
        }
        // On the free surface a value stands for half a cell, so the same source moves it twice
        // as much.
        let plane = g.plane_len();
        for (field, index, w) in &mut terms {
            *w /= solver.cell_share(*field, *index / plane - g.halo);
        }
        Self {
            node,
            mechanism,
            wavelet,
            op: PointOp { terms },
            series: None,
        }
    }

    /// When, in a step that starts at t, the source enters: a force with the velocities at t,
    /// a moment tensor with the stresses at t + dt/2.
    pub fn offset(&self, dt: f64) -> f64 {
        if self.is_force() { 0.0 } else { 0.5 * dt }
    }

    /// Pre-warps the wavelet for a run of `steps` steps of `dt` (crate::tdt), keeping
    /// frequencies up to `f_max`, so that un-warping the traces afterwards removes leapfrog's
    /// time error.
    pub fn remove_time_dispersion(&mut self, dt: f64, steps: usize, f_max: f64) {
        let t0 = self.offset(dt);
        let samples: Vec<f64> = (0..steps)
            .map(|n| self.wavelet.value(t0 + n as f64 * dt))
            .collect();
        self.series = Some(Warp::to_leapfrog(dt, t0, dt, steps, f_max).apply(&samples));
    }

    /// Whether the source acts on velocity (a force) rather than stress (a moment tensor).
    pub fn is_force(&self) -> bool {
        matches!(self.mechanism, Mechanism::Force(_))
    }

    /// Adds the source during step `n`, at time `t` (see [`Source::offset`]).
    pub fn inject(&self, fields: &mut Fields, n: usize, t: f64) {
        let value = match &self.series {
            Some(series) => series.get(n).copied().unwrap_or(0.0),
            None => self.wavelet.value(t),
        };
        self.op.inject(fields, value as f32);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn integrate(w: &Wavelet, f: impl Fn(&Wavelet, f64) -> f64) -> f64 {
        let dt = 1e-5;
        (0..100_000).map(|n| f(w, n as f64 * dt) * dt).sum()
    }

    #[test]
    fn ricker_peaks_at_one_and_has_zero_mean() {
        let w = Wavelet::Ricker { f0: 10.0, t0: 0.12 };
        assert!((w.value(0.12) - 1.0).abs() < 1e-12);
        assert!(integrate(&w, Wavelet::value).abs() < 1e-6);
    }

    #[test]
    fn gaussian_has_unit_area() {
        let w = Wavelet::Gaussian { f0: 10.0, t0: 0.12 };
        assert!((integrate(&w, Wavelet::value) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn derivatives_match_finite_differences() {
        for w in [
            Wavelet::Ricker { f0: 7.0, t0: 0.2 },
            Wavelet::Gaussian { f0: 7.0, t0: 0.2 },
        ] {
            for t in [0.1, 0.17, 0.2, 0.26] {
                let e = 1e-6;
                let fd = (w.value(t + e) - w.value(t - e)) / (2.0 * e);
                assert!(
                    (fd - w.derivative(t)).abs() < 1e-5 * (1.0 + fd.abs()),
                    "{w:?} at {t}"
                );
            }
        }
    }

    #[test]
    fn double_couples_are_traceless() {
        let m = double_couple(30.0, 60.0, 90.0, 1.0);
        assert!((m[0] + m[1] + m[2]).abs() < 1e-12);
        // A vertical strike-slip fault striking north is pure xy.
        let m = double_couple(0.0, 90.0, 0.0, 2.0);
        let expect = [0.0, 0.0, 0.0, 2.0, 0.0, 0.0];
        assert!(
            m.iter().zip(expect).all(|(a, b)| (a - b).abs() < 1e-12),
            "{m:?}"
        );
    }

    use crate::cpml::PmlSettings;
    use crate::grid::Grid;
    use crate::model::Model;
    use crate::receiver::Component;

    fn solver(free_surface: bool) -> Solver {
        let g = Grid::new([12, 12, 12], [10.0; 3]);
        let pml = PmlSettings {
            power: 2.0,
            reflection: None,
            kappa_max: 1.0,
            alpha_max: 30.0,
            vp_max: 3000.0,
            sponge: 0.0,
        };
        Solver::new(
            Model::homogeneous(g, 3000.0, 1700.0, 2500.0),
            1e-3,
            3,
            free_surface,
            &pml,
        )
    }

    fn total(s: &Source) -> f32 {
        s.op.terms.iter().map(|t| t.2).sum()
    }

    fn ricker() -> Wavelet {
        Wavelet::Ricker { f0: 10.0, t0: 0.12 }
    }

    #[test]
    fn an_explosion_loads_the_normal_stresses() {
        let s = solver(false);
        let src = Source::new(
            &s,
            [6, 6, 6],
            Mechanism::Moment([1e9, 1e9, 1e9, 0.0, 0.0, 0.0]),
            ricker(),
        );
        let w = (-1e-3 * 1e9 / 1000.0) as f32;
        let p = s.grid.idx(6, 6, 6);
        assert_eq!(
            src.op.terms,
            vec![(Field::Sxx, p, w), (Field::Syy, p, w), (Field::Szz, p, w)]
        );
    }

    #[test]
    fn a_shear_moment_spreads_over_sixteen_positions() {
        let s = solver(false);
        let src = Source::new(
            &s,
            [6, 6, 6],
            Mechanism::Moment([0.0, 0.0, 0.0, 1e9, 0.0, 0.0]),
            ricker(),
        );
        assert_eq!(src.op.terms.len(), 16);
        assert!(src.op.terms.iter().all(|t| t.0 == Field::Sxy));
        assert!((total(&src) - (-1e-3 * 1e9 / 1000.0) as f32).abs() < 1e-3);
    }

    /// In a homogeneous model a force is the receiver at its node, scaled by dt F/(ρV).
    #[test]
    fn a_force_is_a_scaled_receiver() {
        let s = solver(true);
        for node in [[6, 6, 6], [6, 6, 1]] {
            let src = Source::new(&s, node, Mechanism::Force([0.0, 0.0, 2e9]), ricker());
            let rec = s.receiver_op(Component::Vz, node);
            let scale = (1e-3 / 2500.0 * 2e9 / 1000.0) as f32;
            assert_eq!(src.op.terms.len(), rec.terms.len());
            for (a, b) in src.op.terms.iter().zip(&rec.terms) {
                assert_eq!((a.0, a.1), (b.0, b.1));
                assert!((a.2 - scale * b.2).abs() < 1e-6 * scale, "{a:?} vs {b:?}");
            }
        }
    }

    /// On the free surface vx stands for half a cell, so a horizontal force there moves it
    /// twice as much; vz sits half a node below, a whole cell, and doesn't.
    #[test]
    fn a_force_on_the_free_surface_acts_on_half_a_cell() {
        let (free, open) = (solver(true), solver(false));
        let x = Mechanism::Force([1e9, 0.0, 0.0]);
        let on_free = Source::new(&free, [6, 6, 0], x, ricker());
        let on_open = Source::new(&open, [6, 6, 0], x, ricker());
        for (a, b) in on_free.op.terms.iter().zip(&on_open.op.terms) {
            assert!((a.2 - 2.0 * b.2).abs() < 1e-6 * b.2.abs());
        }
        let z = Source::new(
            &free,
            [6, 6, 0],
            Mechanism::Force([0.0, 0.0, 1e9]),
            ricker(),
        );
        assert!((total(&z) - (1e-3 / 2500.0 * 1e9 / 1000.0) as f32).abs() < 1e-3);
    }

    #[test]
    fn vertical_shear_sources_vanish_on_the_free_surface() {
        let s = solver(true);
        let m = Mechanism::Moment([0.0, 0.0, 0.0, 0.0, 1e9, 1e9]);
        assert!(Source::new(&s, [6, 6, 0], m, ricker()).op.terms.is_empty());
        assert!(!Source::new(&s, [6, 6, 1], m, ricker()).op.terms.is_empty());
    }
}
