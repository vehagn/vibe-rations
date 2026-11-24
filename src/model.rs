//! Elastic models on the grid's nodes: P- and S-wave speed and density.

use std::path::Path;

use std::f64::consts::PI;
use std::sync::OnceLock;

use anyhow::{Context, Result, bail, ensure};

use crate::config::{Interfaces, LayerConfig, ModelConfig, ModelFiles, Parameter};
use crate::grid::Grid;
use crate::npy;

/// How far either side of an interface its band-limited step reaches, in nodes. Measured with
/// plane waves at an interface (`just interfaces`): ±2 doubles the error of ±4, ±8 barely
/// improves on it.
pub const BAND_HALF_WIDTH: f64 = 4.0;

/// The most band-limiting may stiffen a modulus beyond the stiffer layer's, as a factor. The
/// band-limited step overshoots by 7.3%, which across a contrast of r in compliance takes the
/// stiffer side's compliance c to c(1 - 0.073(r - 1)): tiny or negative for r near 14.7 (a soft
/// sediment on rock, water on rock), and stiff enough to cut the time step long before. An
/// interface whose contrast would stiffen beyond this, r above 3.7, is cell-averaged instead
/// (`band_limits`), which never overshoots and still puts it where it lies. Measured with
/// plane waves (`just interfaces`): at a contrast of 48 in shear modulus, cell averages
/// reflected S within 2.4–9.1% where a sinc step clamped to keep the moduli in range was
/// off by 6.1–11.9%; at 2.8 in P modulus, the sinc's 1.2% (P) and 12.6% (S) beat cell
/// averages' 2.9–4.2% and 18–22%.
const MAX_STIFFENING: f64 = 1.25;

/// The unit step at u = 0 (u in nodes), band-limited to the grid's Nyquist wavenumber: the
/// integral of a sinc tapered by a Hann window over ±[`BAND_HALF_WIDTH`] nodes, normalised.
/// 0 below -w and 1 above w, ½ + (½ - f(-u)) between, and it overshoots by 7.3% a node either
/// side of the interface. Band-limiting the model to the grid this way is Mittet's (2017)
/// remedy for interfaces between nodes; the taper and its width are ours.
pub fn band_limited_step(u: f64) -> f64 {
    let w = BAND_HALF_WIDTH;
    if u <= -w {
        return 0.0;
    }
    if u >= w {
        return 1.0;
    }
    static TOTAL: OnceLock<f64> = OnceLock::new();
    let total = *TOTAL.get_or_init(|| 2.0 * kernel_integral(w));
    0.5 + u.signum() * kernel_integral(u.abs()) / total
}

/// ∫ from 0 to `a` of the tapered sinc, by Simpson's rule on 64 intervals or more per node.
fn kernel_integral(a: f64) -> f64 {
    let w = BAND_HALF_WIDTH;
    let kernel = |x: f64| {
        let sinc = if x == 0.0 {
            1.0
        } else {
            (PI * x).sin() / (PI * x)
        };
        sinc * (0.5 + 0.5 * (PI * x / w).cos())
    };
    let n = 2 * ((32.0 * a).ceil() as usize).max(1);
    let h = a / n as f64;
    let inner: f64 = (1..n)
        .map(|m| kernel(m as f64 * h) * if m % 2 == 1 { 4.0 } else { 2.0 })
        .sum();
    (kernel(0.0) + kernel(a) + inner) * h / 3.0
}

/// Whether the interface between `above` and `below` takes the band-limited step: unless
/// the largest ratio of the compliances it averages (P, and S between two solids) would let
/// the step's overshoot stiffen a modulus by more than [`MAX_STIFFENING`].
fn band_limits(above: &LayerConfig, below: &LayerConfig) -> bool {
    let ratio = |a: f64, b: f64| a.max(b) / a.min(b);
    let p = |l: &LayerConfig| l.rho * l.vp * l.vp;
    let s = |l: &LayerConfig| l.rho * l.vs * l.vs;
    let mut r = ratio(p(above), p(below));
    if above.vs > 0.0 && below.vs > 0.0 {
        r = r.max(ratio(s(above), s(below)));
    }
    let overshoot = band_limited_step(1.0) - 1.0;
    1.0 - overshoot * (r - 1.0) >= 1.0 / MAX_STIFFENING
}

/// [vp, vs, ρ] at each depth of nodes, k = 0 to nz - 1, of the layers `cfg` describes.
///
/// Band-limited, each layer weighs the difference of the steps at its top and its bottom
/// ([`band_limited_step`], or the cell's share where the contrast is too large for its
/// overshoot: [`MAX_STIFFENING`]), so the weights sum to one; density is their weighted mean
/// and so is each compliance, the average that puts an interface where it lies (moduli or
/// their logarithms leave it off by up to 0.13 node). A fluid stays where node sampling puts it: a node is fluid if the layer
/// it lies in is. A solid node takes the S compliance of the solids alone, since a fluid's is
/// infinite, and near a fluid at most half its P modulus as shear modulus, so that its bulk
/// modulus stays positive (λ ≥ 0) where the fluid's P compliance softens it.
fn layer_profile(grid: Grid, cfg: &ModelConfig) -> Result<Vec<[f64; 3]>> {
    let mut layers = cfg.layers.clone();
    layers.sort_by(|a, b| a.top.total_cmp(&b.top));
    ensure!(
        layers[0].top <= 0.0,
        "the shallowest layer must start at depth 0 or above"
    );
    let h = grid.h[2];
    // The deepest layer whose top is at or above depth z; a node exactly on an interface
    // takes the layer below it.
    let lies_in = |z: f64| {
        layers
            .iter()
            .rposition(|l| l.top <= z + 1e-6 * h)
            .expect("the first layer starts at the surface")
    };
    let depths = (0..grid.n[2]).map(|k| k as f64 * h);
    if cfg.interfaces == Interfaces::Nodes {
        return Ok(depths
            .map(|z| {
                let l = &layers[lies_in(z)];
                [l.vp, l.vs, l.rho]
            })
            .collect());
    }
    let sinc: Vec<bool> = layers
        .windows(2)
        .map(|pair| band_limits(&pair[0], &pair[1]))
        .collect();
    let profile = depths
        .map(|z| {
            // The step at layer j's top, j ≥ 1; the first layer reaches up forever.
            let step = |j: usize| -> f64 {
                if j == 0 {
                    return 1.0;
                }
                if j == layers.len() {
                    return 0.0;
                }
                let u = (z - layers[j].top) / h;
                if sinc[j - 1] {
                    band_limited_step(u)
                } else {
                    // The cell's share below the interface: a node on it is half each.
                    (u + 0.5).clamp(0.0, 1.0)
                }
            };
            let weights: Vec<f64> = (0..layers.len()).map(|j| step(j) - step(j + 1)).collect();
            let sum = |f: &dyn Fn(&LayerConfig) -> f64, solid_only: bool| -> f64 {
                layers
                    .iter()
                    .zip(&weights)
                    .filter(|(l, _)| !solid_only || l.vs > 0.0)
                    .map(|(l, w)| w * f(l))
                    .sum()
            };
            let rho = sum(&|l| l.rho, false);
            let p_modulus = 1.0 / sum(&|l| 1.0 / (l.rho * l.vp * l.vp), false);
            let vs = if layers[lies_in(z)].vs == 0.0 {
                0.0
            } else {
                let solid = sum(&|_| 1.0, true);
                let mut mu = solid / sum(&|l| 1.0 / (l.rho * l.vs * l.vs), true);
                if solid != 1.0 && layers.iter().any(|l| l.vs == 0.0) {
                    mu = mu.min(0.5 * p_modulus);
                }
                (mu / rho).sqrt()
            };
            [(p_modulus / rho).sqrt(), vs, rho]
        })
        .collect();
    Ok(profile)
}

/// Where a configuration's `[model.files]` come from: the disk, relative to a directory, or
/// bytes held in memory under the names the configuration gives them, as in a browser, which
/// has no files.
#[derive(Clone, Copy, Debug)]
pub enum Files<'a> {
    Dir(&'a Path),
    Memory(&'a [(String, Vec<u8>)]),
}

impl Files<'_> {
    /// The array `name` holds, and its shape.
    fn read(&self, name: &Path) -> Result<(Vec<usize>, Vec<f32>)> {
        match self {
            Files::Dir(base) => {
                let path = base.join(name);
                npy::read(&path).with_context(|| format!("reading {}", path.display()))
            }
            Files::Memory(files) => {
                let (_, bytes) = files
                    .iter()
                    .find(|(n, _)| Path::new(n) == name)
                    .with_context(|| format!("{} wasn't given", name.display()))?;
                npy::parse(bytes).with_context(|| format!("reading {}", name.display()))
            }
        }
    }
}

/// Whether an array of `shape` is a section, (nz, nx), the same at every y, rather than the
/// whole grid, (nz, ny, nx); an error if it is neither.
fn is_section(shape: &[usize], grid: Grid, name: &Path) -> Result<bool> {
    let [nx, _, nz] = grid.n;
    if shape == grid.shape() {
        Ok(false)
    } else if shape == [nz, nx] {
        Ok(true)
    } else {
        bail!(
            "{} has shape {shape:?}; the grid needs (nz, ny, nx) = {:?}, or (nz, nx) = {:?} \
             for a section the same at every y",
            name.display(),
            grid.shape(),
            [nz, nx]
        )
    }
}

/// The files a configured model is read from, or None if it is layered: it takes one or the
/// other.
fn model_files(cfg: &ModelConfig) -> Result<Option<&ModelFiles>> {
    match (&cfg.files, cfg.layers.is_empty()) {
        (Some(_), false) => bail!("the model takes either layers or files, not both"),
        (None, true) => bail!("the model needs [[model.layer]]s or [model.files]"),
        (files, _) => Ok(files.as_ref()),
    }
}

/// Isotropic elastic properties at every interior node, x fastest.
#[derive(Clone, Debug)]
pub struct Model {
    pub grid: Grid,
    /// P-wave speed, m/s.
    pub vp: Vec<f32>,
    /// S-wave speed, m/s; 0 in a fluid.
    pub vs: Vec<f32>,
    /// Density, kg/m³.
    pub rho: Vec<f32>,
}

impl Model {
    pub fn homogeneous(grid: Grid, vp: f64, vs: f64, rho: f64) -> Self {
        let n = grid.interior_len();
        Self {
            grid,
            vp: vec![vp as f32; n],
            vs: vec![vs as f32; n],
            rho: vec![rho as f32; n],
        }
    }

    /// Builds the model a configuration describes; `base` is the directory `.npy` paths are
    /// relative to.
    pub fn from_config(grid: Grid, cfg: &ModelConfig, base: &Path) -> Result<Self> {
        Self::from_files(grid, cfg, Files::Dir(base))
    }

    /// Builds the model a configuration describes, its files from `files`. A file is either
    /// the whole grid, (nz, ny, nx), or a section, (nz, nx), repeated along y: a 2.5D model.
    pub fn from_files(grid: Grid, cfg: &ModelConfig, files: Files) -> Result<Self> {
        let mut model = match model_files(cfg)? {
            None => Self::layered(grid, cfg)?,
            Some(names) => {
                let load = |p: &Path| -> Result<Vec<f32>> {
                    let (shape, data) = files.read(p)?;
                    if !is_section(&shape, grid, p)? {
                        return Ok(data);
                    }
                    let mut all = Vec::with_capacity(grid.interior_len());
                    for row in data.chunks_exact(grid.n[0]) {
                        for _ in 0..grid.n[1] {
                            all.extend_from_slice(row);
                        }
                    }
                    Ok(all)
                };
                Self {
                    grid,
                    vp: load(&names.vp)?,
                    vs: load(&names.vs)?,
                    rho: load(&names.rho)?,
                }
            }
        };
        for s in &cfg.spheres {
            ensure!(s.radius > 0.0, "a sphere's radius must be positive");
            model.fill_sphere(s.center, s.radius, [s.vp, s.vs, s.rho]);
        }
        model.validate()?;
        Ok(model)
    }

    /// The fastest P speed and the slowest wave (S in solids, P in fluids) of the model `cfg`
    /// describes, found without building it, so that a model too large to hold can still be
    /// planned: from the layers that reach a node and the spheres that contain one, or from
    /// the files, read through a piece at a time. A layer that spheres cover entirely still
    /// counts, so in that rare case the answer is conservative.
    pub fn extremes(grid: Grid, cfg: &ModelConfig, base: &Path) -> Result<[f64; 2]> {
        Self::extremes_in(grid, cfg, Files::Dir(base))
    }

    /// [`Model::extremes`] with the files from `files`: files on disk are streamed, files in
    /// memory read whole.
    pub fn extremes_in(grid: Grid, cfg: &ModelConfig, files: Files) -> Result<[f64; 2]> {
        let slow = |vp: f64, vs: f64| if vs > 0.0 { vs } else { vp };
        let (mut fast, mut slowest) = (f64::NEG_INFINITY, f64::INFINITY);
        let mut take = |vp: f64, vs: f64| {
            fast = fast.max(vp);
            slowest = slowest.min(slow(vp, vs));
        };
        match model_files(cfg)? {
            None => {
                for [vp, vs, _] in layer_profile(grid, cfg)? {
                    // As the model stores them.
                    take(vp as f32 as f64, vs as f32 as f64);
                }
            }
            // A section's values are the grid's, so its extremes are too.
            Some(names) => match files {
                Files::Dir(base) => {
                    let open = |p: &Path| -> Result<npy::Stream> {
                        let path = base.join(p);
                        let s = npy::stream(&path)
                            .with_context(|| format!("reading {}", path.display()))?;
                        is_section(&s.shape, grid, &path)?;
                        Ok(s)
                    };
                    let (mut vp, mut vs) = (open(&names.vp)?, open(&names.vs)?);
                    ensure!(
                        vp.shape == vs.shape,
                        "{} and {} differ in shape",
                        names.vp.display(),
                        names.vs.display()
                    );
                    let (mut a, mut b) = (Vec::new(), Vec::new());
                    while vp.next_chunk(&mut a, 1 << 16)? {
                        vs.next_chunk(&mut b, 1 << 16)?;
                        for (&p, &s) in a.iter().zip(&b) {
                            take(p as f64, s as f64);
                        }
                    }
                }
                Files::Memory(_) => {
                    let read = |p: &Path| -> Result<Vec<f32>> {
                        let (shape, data) = files.read(p)?;
                        is_section(&shape, grid, p)?;
                        Ok(data)
                    };
                    let (vp, vs) = (read(&names.vp)?, read(&names.vs)?);
                    ensure!(
                        vp.len() == vs.len(),
                        "{} and {} differ in shape",
                        names.vp.display(),
                        names.vs.display()
                    );
                    for (&p, &s) in vp.iter().zip(&vs) {
                        take(p as f64, s as f64);
                    }
                }
            },
        }
        for s in &cfg.spheres {
            // The node nearest the centre: the distance is a sum over the axes, so the nearest
            // index along each one.
            let r2: f64 = (0..3)
                .map(|a| {
                    let i = (s.center[a] / grid.h[a])
                        .round()
                        .clamp(0.0, (grid.n[a] - 1) as f64);
                    (i * grid.h[a] - s.center[a]).powi(2)
                })
                .sum();
            if r2 <= s.radius * s.radius {
                take(s.vp, s.vs);
            }
        }
        Ok([fast, slowest])
    }

    fn layered(grid: Grid, cfg: &ModelConfig) -> Result<Self> {
        let mut model = Self::homogeneous(grid, 0.0, 0.0, 0.0);
        let plane = grid.n[0] * grid.n[1];
        for (k, [vp, vs, rho]) in layer_profile(grid, cfg)?.into_iter().enumerate() {
            let range = k * plane..(k + 1) * plane;
            model.vp[range.clone()].fill(vp as f32);
            model.vs[range.clone()].fill(vs as f32);
            model.rho[range].fill(rho as f32);
        }
        Ok(model)
    }

    /// Sets [vp, vs, rho] at every node within `radius` of `center`.
    pub fn fill_sphere(&mut self, center: [f64; 3], radius: f64, props: [f64; 3]) {
        let g = self.grid;
        for k in 0..g.n[2] {
            for j in 0..g.n[1] {
                for i in 0..g.n[0] {
                    let p = g.position([i, j, k]);
                    let r2: f64 = (0..3).map(|a| (p[a] - center[a]).powi(2)).sum();
                    if r2 <= radius * radius {
                        let n = g.interior_idx(i, j, k);
                        self.vp[n] = props[0] as f32;
                        self.vs[n] = props[1] as f32;
                        self.rho[n] = props[2] as f32;
                    }
                }
            }
        }
    }

    /// Checks that every node is physical: positive density and P speed, a non-negative S
    /// speed, and a positive bulk modulus (vp² > 4/3 vs²).
    pub fn validate(&self) -> Result<()> {
        let n = self.grid.interior_len();
        ensure!(
            self.vp.len() == n && self.vs.len() == n && self.rho.len() == n,
            "model arrays don't match the grid"
        );
        for idx in 0..n {
            let (vp, vs, rho) = (self.vp[idx], self.vs[idx], self.rho[idx]);
            let ok = vp.is_finite()
                && vs.is_finite()
                && rho.is_finite()
                && vp > 0.0
                && vs >= 0.0
                && rho > 0.0
                && 3.0 * vp * vp > 4.0 * vs * vs;
            if !ok {
                let i = idx % self.grid.n[0];
                let j = idx / self.grid.n[0] % self.grid.n[1];
                let k = idx / (self.grid.n[0] * self.grid.n[1]);
                anyhow::bail!(
                    "unphysical model at node ({i}, {j}, {k}): vp {vp}, vs {vs}, rho {rho}"
                );
            }
        }
        Ok(())
    }

    /// The values of `p` at every node.
    pub fn values(&self, p: Parameter) -> &[f32] {
        match p {
            Parameter::Vp => &self.vp,
            Parameter::Vs => &self.vs,
            Parameter::Rho => &self.rho,
        }
    }

    /// The values of `p` at every node, to change.
    pub fn values_mut(&mut self, p: Parameter) -> &mut [f32] {
        match p {
            Parameter::Vp => &mut self.vp,
            Parameter::Vs => &mut self.vs,
            Parameter::Rho => &mut self.rho,
        }
    }

    /// Bytes its arrays hold.
    pub fn memory_bytes(&self) -> usize {
        4 * (self.vp.len() + self.vs.len() + self.rho.len())
    }

    /// Writes `vp.npy`, `vs.npy` and `rho.npy` into `dir`, shaped (nz, ny, nx), ready for
    /// `[model.files]`.
    pub fn write(&self, dir: &Path) -> Result<()> {
        for p in Parameter::ALL {
            npy::write(
                &dir.join(format!("{}.npy", p.name())),
                &self.grid.shape(),
                self.values(p),
            )?;
        }
        Ok(())
    }

    pub fn vp_range(&self) -> [f64; 2] {
        range(&self.vp)
    }

    pub fn vs_range(&self) -> [f64; 2] {
        range(&self.vs)
    }

    pub fn rho_range(&self) -> [f64; 2] {
        range(&self.rho)
    }

    /// The slowest wave anywhere: S in solids, P in fluids. It sets the shortest wavelength.
    pub fn slowest_wave(&self) -> f64 {
        self.vp
            .iter()
            .zip(&self.vs)
            .map(|(&vp, &vs)| if vs > 0.0 { vs } else { vp })
            .fold(f32::INFINITY, f32::min) as f64
    }
}

fn range(v: &[f32]) -> [f64; 2] {
    let (lo, hi) = v
        .iter()
        .fold((f32::INFINITY, f32::NEG_INFINITY), |(lo, hi), &x| {
            (lo.min(x), hi.max(x))
        });
    [lo as f64, hi as f64]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::{Interfaces, LayerConfig, SphereConfig};

    fn layer(top: f64, vp: f64) -> LayerConfig {
        LayerConfig {
            top,
            vp,
            vs: vp / 2.0,
            rho: 2000.0,
        }
    }

    #[test]
    fn layers_cover_the_nodes_below_their_tops() {
        let g = Grid::new([4, 4, 10], [10.0; 3]);
        let cfg = ModelConfig {
            // Out of order on purpose.
            layers: vec![layer(30.0, 3000.0), layer(0.0, 2000.0)],
            files: None,
            spheres: vec![],
            interfaces: Interfaces::Nodes,
        };
        let m = Model::from_config(g, &cfg, Path::new(".")).unwrap();
        assert_eq!(m.vp[g.interior_idx(0, 0, 2)], 2000.0);
        assert_eq!(m.vp[g.interior_idx(0, 0, 3)], 3000.0); // on the interface: the layer below
        assert_eq!(m.vp[g.interior_idx(3, 3, 9)], 3000.0);
        assert_eq!(m.vp_range(), [2000.0, 3000.0]);
        assert_eq!(m.slowest_wave(), 1000.0);
    }

    #[test]
    fn the_band_limited_step_is_odd_and_overshoots_by_seven_percent() {
        assert_eq!(band_limited_step(0.0), 0.5);
        assert_eq!(band_limited_step(-BAND_HALF_WIDTH), 0.0);
        assert_eq!(band_limited_step(BAND_HALF_WIDTH), 1.0);
        for u in [0.3, 0.5, 1.0, 1.7, 2.5, 3.9] {
            let (a, b) = (band_limited_step(u), band_limited_step(-u));
            assert!((a + b - 1.0).abs() < 1e-12, "{u}: {a} + {b}");
        }
        let peak = (0..=800)
            .map(|m| band_limited_step(m as f64 / 200.0))
            .fold(0.0, f64::max);
        assert!((peak - 1.0733).abs() < 1e-4, "{peak}");
        // Smooth enough that the step doesn't jump where the quadrature gains an interval.
        let (a, b) = (band_limited_step(0.999_999), band_limited_step(1.000_001));
        assert!((a - b).abs() < 1e-9, "{a} {b}");
    }

    fn band_limited(layers: Vec<LayerConfig>, nz: usize) -> (Grid, ModelConfig, Model) {
        let g = Grid::new([2, 2, nz], [10.0; 3]);
        let cfg = ModelConfig {
            layers,
            files: None,
            spheres: vec![],
            interfaces: Interfaces::BandLimited,
        };
        let m = Model::from_config(g, &cfg, Path::new(".")).unwrap();
        // What a run is planned with is what it gets.
        let extremes = Model::extremes(g, &cfg, Path::new(".")).unwrap();
        assert_eq!(extremes, [m.vp_range()[1], m.slowest_wave()]);
        (g, cfg, m)
    }

    fn at(m: &Model, k: usize) -> [f64; 3] {
        let n = m.grid.interior_idx(0, 0, k);
        [m.vp[n], m.vs[n], m.rho[n]].map(f64::from)
    }

    fn compliances([vp, vs, rho]: [f64; 3]) -> [f64; 2] {
        [1.0 / (rho * vp * vp), 1.0 / (rho * vs * vs)]
    }

    #[test]
    fn band_limited_layers_follow_their_interfaces_between_the_nodes() {
        let (upper, lower) = (layer(0.0, 2000.0), layer(100.0, 3000.0));
        let (_, _, m) = band_limited(vec![upper.clone(), lower.clone()], 21);
        let props = |l: &LayerConfig| [l.vp, l.vs, l.rho];
        // Beyond four nodes either side, each layer as it is.
        for k in [0, 5, 15, 20] {
            let want = props(if k < 10 { &upper } else { &lower });
            for (a, b) in at(&m, k).iter().zip(want) {
                assert!((a - b).abs() < 1e-3 * b, "node {k}: {a} vs {b}");
            }
        }
        // On the interface, the mean density and the mean compliances.
        let [c, d] = compliances(at(&m, 10));
        let [c1, d1] = compliances(props(&upper));
        let [c2, d2] = compliances(props(&lower));
        assert!((c / (0.5 * (c1 + c2)) - 1.0).abs() < 1e-6);
        assert!((d / (0.5 * (d1 + d2)) - 1.0).abs() < 1e-6);
        assert!((at(&m, 10)[2] - 2000.0).abs() < 1e-3);
        // Moved by a node, the profile moves by a node; moved by part of one, it changes.
        let moved = |top: f64| {
            let (_, _, m) = band_limited(vec![upper.clone(), layer(top, 3000.0)], 21);
            m
        };
        let (a, b) = (moved(100.0), moved(110.0));
        for k in 1..21 {
            assert_eq!(at(&a, k - 1), at(&b, k));
        }
        assert_ne!(at(&moved(103.0), 10), at(&a, 10));
    }

    #[test]
    fn band_limited_layers_keep_a_fluid_where_the_nodes_put_it() {
        let water = LayerConfig {
            top: 0.0,
            vp: 1500.0,
            vs: 0.0,
            rho: 1000.0,
        };
        // 2.7 times water's P modulus, so the interface takes the sinc step.
        let sediment = LayerConfig {
            top: 50.0,
            vp: 1800.0,
            vs: 1000.0,
            rho: 1900.0,
        };
        let (_, _, m) = band_limited(vec![water, sediment], 15);
        for k in 0..5 {
            assert_eq!(at(&m, k)[1], 0.0, "node {k} should be water");
        }
        // The node on the seafloor is sediment, its shear modulus held to half its P modulus;
        // further down the sediment's own.
        let [vp, vs, _] = at(&m, 5);
        assert!(
            vs > 0.0 && vs * vs <= 0.5 * vp * vp * (1.0 + 1e-6),
            "{vp} {vs}"
        );
        assert!((at(&m, 9)[1] - 1000.0).abs() < 1e-3 * 1000.0);
        // Water's P compliance reaches into the sediment and the sediment's into the water.
        assert!(at(&m, 5)[0] < 1800.0 && at(&m, 4)[0] > 1500.0);
    }

    #[test]
    fn band_limiting_stiffens_no_modulus_much_at_a_large_contrast() {
        // A soft sediment on rock: 48 times the shear modulus. The sinc step's overshoot would
        // make the rock's shear compliance negative next to the interface, so it is
        // cell-averaged: nothing stiffer than the rock, and the node 3 m above the interface
        // takes 20% rock.
        let sediment = LayerConfig {
            top: 0.0,
            vp: 1800.0,
            vs: 400.0,
            rho: 1900.0,
        };
        let rock = LayerConfig {
            top: 103.0,
            vp: 4500.0,
            vs: 2600.0,
            rho: 2600.0,
        };
        let (_, _, m) = band_limited(vec![sediment.clone(), rock.clone()], 21);
        let [c_rock, d_rock] = compliances([rock.vp, rock.vs, rock.rho]);
        for k in 0..21 {
            let [c, d] = compliances(at(&m, k));
            assert!(
                c >= c_rock * (1.0 - 1e-6) && d >= d_rock * (1.0 - 1e-6),
                "node {k}"
            );
        }
        let [_, d_sediment] = compliances([sediment.vp, sediment.vs, sediment.rho]);
        let [_, d] = compliances(at(&m, 10));
        assert!((d / (0.8 * d_sediment + 0.2 * d_rock) - 1.0).abs() < 1e-5);
        assert_eq!(at(&m, 9)[1], sediment.vs);
        assert_eq!(at(&m, 11)[1], rock.vs);
    }

    #[test]
    fn spheres_overwrite_the_layers() {
        let g = Grid::new([11, 11, 11], [1.0; 3]);
        let cfg = ModelConfig {
            layers: vec![layer(0.0, 2000.0)],
            files: None,
            spheres: vec![SphereConfig {
                center: [5.0, 5.0, 5.0],
                radius: 2.0,
                vp: 1500.0,
                vs: 0.0,
                rho: 1000.0,
            }],
            interfaces: Interfaces::Nodes,
        };
        let m = Model::from_config(g, &cfg, Path::new(".")).unwrap();
        assert_eq!(m.vs[g.interior_idx(5, 5, 7)], 0.0);
        assert_eq!(m.vs[g.interior_idx(5, 5, 8)], 1000.0);
        assert_eq!(m.slowest_wave(), 1000.0);
    }

    #[test]
    fn rejects_unphysical_speeds() {
        let g = Grid::new([8, 8, 8], [1.0; 3]);
        assert!(
            Model::homogeneous(g, 1000.0, 900.0, 2000.0)
                .validate()
                .is_err()
        );
        assert!(
            Model::homogeneous(g, 1000.0, 500.0, 0.0)
                .validate()
                .is_err()
        );
        assert!(
            Model::homogeneous(g, 1500.0, 0.0, 1000.0)
                .validate()
                .is_ok()
        );
    }

    use crate::config::ModelFiles;

    #[test]
    fn loads_npy_files_relative_to_the_configuration() {
        let dir = std::env::temp_dir().join(format!("vibe-rations-model-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let g = Grid::new([4, 3, 2], [10.0; 3]);
        let ramp = |base: f32| (0..24).map(|n| base + n as f32).collect::<Vec<_>>();
        for (name, base) in [("vp", 3000.0), ("vs", 1500.0), ("rho", 2500.0)] {
            npy::write(&dir.join(format!("{name}.npy")), &[2, 3, 4], &ramp(base)).unwrap();
        }
        let cfg = ModelConfig {
            layers: vec![],
            files: Some(ModelFiles {
                vp: "vp.npy".into(),
                vs: "vs.npy".into(),
                rho: "rho.npy".into(),
            }),
            spheres: vec![],
            interfaces: Interfaces::Nodes,
        };
        let m = Model::from_config(g, &cfg, &dir).unwrap();
        // x fastest: node (3, 2, 1) is the last value.
        assert_eq!(m.vp[g.interior_idx(3, 2, 1)], 3023.0);
        assert_eq!(m.rho[g.interior_idx(1, 0, 0)], 2501.0);
        let wrong = Grid::new([3, 4, 2], [10.0; 3]);
        let e = format!("{:#}", Model::from_config(wrong, &cfg, &dir).unwrap_err());
        assert!(e.contains("shape"), "{e}");
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn repeats_a_section_along_y_from_disk_or_memory() {
        let dir = std::env::temp_dir().join(format!("vibe-rations-section-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let g = Grid::new([4, 3, 2], [10.0; 3]);
        // (nz, nx) = (2, 4): a value per x and z.
        let section = |base: f32| (0..8).map(|n| base + n as f32).collect::<Vec<_>>();
        let mut held = Vec::new();
        for (name, base) in [("vp", 3000.0), ("vs", 1500.0), ("rho", 2500.0)] {
            let path = dir.join(format!("{name}.npy"));
            npy::write(&path, &[2, 4], &section(base)).unwrap();
            held.push((format!("{name}.npy"), std::fs::read(&path).unwrap()));
        }
        let cfg = ModelConfig {
            layers: vec![],
            files: Some(ModelFiles {
                vp: "vp.npy".into(),
                vs: "vs.npy".into(),
                rho: "rho.npy".into(),
            }),
            spheres: vec![],
            interfaces: Interfaces::Nodes,
        };
        let from_disk = Model::from_config(g, &cfg, &dir).unwrap();
        let from_memory = Model::from_files(g, &cfg, Files::Memory(&held)).unwrap();
        for m in [&from_disk, &from_memory] {
            for j in 0..3 {
                assert_eq!(m.vp[g.interior_idx(3, j, 1)], 3007.0);
                assert_eq!(m.vs[g.interior_idx(1, j, 0)], 1501.0);
            }
        }
        assert_eq!(from_disk.rho, from_memory.rho);
        let extremes = [
            Model::extremes(g, &cfg, &dir).unwrap(),
            Model::extremes_in(g, &cfg, Files::Memory(&held)).unwrap(),
        ];
        assert_eq!(extremes, [[3007.0, 1500.0]; 2]);
        // A file the configuration names but memory doesn't hold.
        let e = format!(
            "{:#}",
            Model::from_files(g, &cfg, Files::Memory(&held[..2])).unwrap_err()
        );
        assert!(e.contains("rho.npy wasn't given"), "{e}");
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn needs_either_layers_or_files() {
        let g = Grid::new([8, 8, 8], [1.0; 3]);
        let files = Some(ModelFiles {
            vp: "a".into(),
            vs: "b".into(),
            rho: "c".into(),
        });
        let none = ModelConfig {
            layers: vec![],
            files: None,
            spheres: vec![],
            interfaces: Interfaces::Nodes,
        };
        assert!(Model::from_config(g, &none, Path::new(".")).is_err());
        let both = ModelConfig {
            layers: vec![layer(0.0, 2000.0)],
            files,
            spheres: vec![],
            interfaces: Interfaces::Nodes,
        };
        let e = format!(
            "{:#}",
            Model::from_config(g, &both, Path::new(".")).unwrap_err()
        );
        assert!(e.contains("not both"), "{e}");
        let deep = ModelConfig {
            layers: vec![layer(10.0, 2000.0)],
            files: None,
            spheres: vec![],
            interfaces: Interfaces::Nodes,
        };
        assert!(Model::from_config(g, &deep, Path::new(".")).is_err());
    }
}
