//! Elastic models on the grid's nodes: P- and S-wave speed and density.

use std::path::Path;

use anyhow::{Context, Result, ensure};

use crate::config::ModelConfig;
use crate::grid::Grid;
use crate::npy;

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
        let mut model = match (&cfg.files, cfg.layers.is_empty()) {
            (Some(_), false) => anyhow::bail!("the model takes either layers or files, not both"),
            (None, true) => anyhow::bail!("the model needs [[model.layer]]s or [model.files]"),
            (None, false) => Self::layered(grid, cfg)?,
            (Some(files), true) => {
                let load = |p: &Path| -> Result<Vec<f32>> {
                    let path = base.join(p);
                    let (shape, data) =
                        npy::read(&path).with_context(|| format!("reading {}", path.display()))?;
                    let want = [grid.n[2], grid.n[1], grid.n[0]];
                    ensure!(
                        shape == want,
                        "{} has shape {shape:?}; the grid needs (nz, ny, nx) = {want:?}",
                        path.display()
                    );
                    Ok(data)
                };
                Self {
                    grid,
                    vp: load(&files.vp)?,
                    vs: load(&files.vs)?,
                    rho: load(&files.rho)?,
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
