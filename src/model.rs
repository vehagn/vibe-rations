//! Elastic models on the grid's nodes: P- and S-wave speed and density.

use std::path::Path;

use anyhow::{Context, Result, bail, ensure};

use crate::config::ModelConfig;
use crate::grid::Grid;
use crate::npy;

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
    let [nx, ny, nz] = grid.n;
    if shape == [nz, ny, nx] {
        Ok(false)
    } else if shape == [nz, nx] {
        Ok(true)
    } else {
        bail!(
            "{} has shape {shape:?}; the grid needs (nz, ny, nx) = {:?}, or (nz, nx) = {:?} \
             for a section the same at every y",
            name.display(),
            [nz, ny, nx],
            [nz, nx]
        )
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
        let mut model = match (&cfg.files, cfg.layers.is_empty()) {
            (Some(_), false) => bail!("the model takes either layers or files, not both"),
            (None, true) => bail!("the model needs [[model.layer]]s or [model.files]"),
            (None, false) => Self::layered(grid, cfg)?,
            (Some(names), true) => {
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
        match (&cfg.files, cfg.layers.is_empty()) {
            (Some(_), false) => bail!("the model takes either layers or files, not both"),
            (None, true) => bail!("the model needs [[model.layer]]s or [model.files]"),
            (None, false) => {
                let mut layers = cfg.layers.clone();
                layers.sort_by(|a, b| a.top.total_cmp(&b.top));
                ensure!(
                    layers[0].top <= 0.0,
                    "the shallowest layer must start at depth 0 or above"
                );
                for k in 0..grid.n[2] {
                    let z = k as f64 * grid.h[2];
                    let layer = layers
                        .iter()
                        .rev()
                        .find(|l| l.top <= z + 1e-6 * grid.h[2])
                        .expect("the first layer starts at the surface");
                    take(layer.vp, layer.vs);
                }
            }
            // A section's values are the grid's, so its extremes are too.
            (Some(names), true) => match files {
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
        let mut layers = cfg.layers.clone();
        layers.sort_by(|a, b| a.top.total_cmp(&b.top));
        ensure!(
            layers[0].top <= 0.0,
            "the shallowest layer must start at depth 0 or above"
        );
        let mut model = Self::homogeneous(grid, 0.0, 0.0, 0.0);
        for k in 0..grid.n[2] {
            let z = k as f64 * grid.h[2];
            // The deepest layer whose top is at or above this node; a node exactly on an
            // interface takes the layer below it.
            let layer = layers
                .iter()
                .rev()
                .find(|l| l.top <= z + 1e-6 * grid.h[2])
                .expect("the first layer starts at the surface");
            let plane = grid.n[0] * grid.n[1];
            let range = k * plane..(k + 1) * plane;
            model.vp[range.clone()].fill(layer.vp as f32);
            model.vs[range.clone()].fill(layer.vs as f32);
            model.rho[range].fill(layer.rho as f32);
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
    use crate::config::{LayerConfig, SphereConfig};

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
        };
        let m = Model::from_config(g, &cfg, Path::new(".")).unwrap();
        assert_eq!(m.vp[g.interior_idx(0, 0, 2)], 2000.0);
        assert_eq!(m.vp[g.interior_idx(0, 0, 3)], 3000.0); // on the interface: the layer below
        assert_eq!(m.vp[g.interior_idx(3, 3, 9)], 3000.0);
        assert_eq!(m.vp_range(), [2000.0, 3000.0]);
        assert_eq!(m.slowest_wave(), 1000.0);
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
        };
        assert!(Model::from_config(g, &none, Path::new(".")).is_err());
        let both = ModelConfig {
            layers: vec![layer(0.0, 2000.0)],
            files,
            spheres: vec![],
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
        };
        assert!(Model::from_config(g, &deep, Path::new(".")).is_err());
    }
}
