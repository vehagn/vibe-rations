//! The memory a run needs, worked out from its configuration without building it, so that a
//! model too large for the machine (or the browser) is caught before anything is allocated.
//!
//! Every large array is counted exactly, from the same quantities the run uses: the wavefield
//! and material on the padded grid, the model on the interior, the CPML memory variables, the
//! sources' pre-warped time series, the traces, and the copy a volume snapshot makes while it is
//! written. What is left (the program, the stencil and receiver operators, per-thread scratch
//! rows) is a few MiB.

use std::fmt::Write as _;
use std::path::Path;

use anyhow::Result;
use serde::Serialize;

use crate::config::{Config, PlaneKind};
use crate::fd::Stencil;
use crate::grid::Grid;
use crate::model::Model;
use crate::simulation::timing;

/// One part of a run's memory.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct Part {
    pub what: &'static str,
    pub bytes: u64,
    /// Held only for a moment (while a snapshot is written), not for the whole run.
    pub transient: bool,
}

/// A run's memory: its parts, what it holds throughout, and its peak.
#[derive(Clone, Debug, Default, PartialEq, Serialize)]
pub struct Memory {
    pub parts: Vec<Part>,
    /// Held for the whole run.
    pub resident: u64,
    /// The most at any one time: the resident memory and the largest transient part.
    pub peak: u64,
}

/// What decides a run's memory.
#[derive(Clone, Debug)]
pub struct Plan {
    /// The grid, with its halo.
    pub grid: Grid,
    /// Absorbing-layer nodes along each axis, both ends together.
    pub layer_nodes: [usize; 3],
    /// Whether the side layers also damp plainly (`pml_sponge` along a free surface).
    pub sponge: bool,
    pub steps: usize,
    /// Sources holding a pre-warped time series, one value per step.
    pub warped_sources: usize,
    /// Receivers times recorded components.
    pub traces: usize,
    pub samples: usize,
    pub volume_snapshots: bool,
}

impl Plan {
    /// The plan of a run of `steps` steps on `grid` (with its halo) as `cfg` sets it up.
    pub fn new(cfg: &Config, grid: Grid, steps: usize) -> Result<Self> {
        let b = &cfg.boundary;
        let w = b.pml_width;
        let receivers = cfg.receivers.positions()?.len();
        Ok(Self {
            grid,
            layer_nodes: [2 * w, 2 * w, if b.free_surface { w } else { 2 * w }],
            sponge: b.free_surface && b.pml_sponge > 0.0 && w > 0,
            steps,
            warped_sources: if cfg.scheme.remove_time_dispersion {
                cfg.sources.len()
            } else {
                0
            },
            traces: receivers * cfg.receivers.components.len(),
            samples: steps.div_ceil(cfg.output.record_every.max(1)),
            volume_snapshots: cfg.snapshots.iter().any(|s| s.plane == PlaneKind::Volume),
        })
    }

    pub fn memory(&self) -> Memory {
        let g = &self.grid;
        let [nx, ny, nz] = g.n.map(|n| n as u64);
        let padded = g.padded_len() as u64;
        let interior = g.interior_len() as u64;
        let [w0, w1, w2] = self.layer_nodes.map(|w| w as u64);
        // Three memory variables per slab node, for velocities and for stresses.
        let psi = 2 * 3 * (nz * ny * w0 + nz * w1 * nx + w2 * ny * nx);
        let sponge = if self.sponge { nx + ny + nz } else { 0 };
        let mut parts = vec![
            part("wavefield: 9 arrays on the padded grid", 4 * 9 * padded),
            part("material: 8 arrays on the padded grid", 4 * 8 * padded),
            part("model: vp, vs and rho", 4 * 3 * interior),
            part("absorbing layers' memory", 4 * (psi + sponge)),
            part(
                "sources' pre-warped time series",
                8 * (self.warped_sources * self.steps) as u64,
            ),
            part("traces", 4 * (self.traces * self.samples) as u64),
        ];
        if self.volume_snapshots {
            // The nodes (three indices each), their values and their shading.
            parts.push(Part {
                what: "a volume snapshot, while it is written",
                bytes: (24 + 4 + 4) * interior,
                transient: true,
            });
        }
        parts.retain(|p| p.bytes > 0);
        let resident = parts.iter().filter(|p| !p.transient).map(|p| p.bytes).sum();
        let peak = resident
            + parts
                .iter()
                .filter(|p| p.transient)
                .map(|p| p.bytes)
                .max()
                .unwrap_or(0);
        Memory {
            parts,
            resident,
            peak,
        }
    }
}

fn part(what: &'static str, bytes: u64) -> Part {
    Part {
        what,
        bytes,
        transient: false,
    }
}

/// The memory the run `cfg` describes would need, without building it: the model's extremes
/// (which set the time step, and so the number of steps and samples) come from its layers and
/// spheres, or from its files read through a piece at a time. `base` is the directory relative
/// paths start from.
pub fn estimate(cfg: &Config, base: &Path) -> Result<Memory> {
    let taylor = Stencil::taylor(cfg.scheme.order)?;
    let grid = cfg.grid.grid()?.with_halo(taylor.half_length());
    let [vp_max, slowest] = Model::extremes(grid, &cfg.model, base)?;
    let t = timing(cfg, &grid, vp_max, slowest)?;
    Ok(Plan::new(cfg, grid, t.steps)?.memory())
}

/// Memory the system can still give this process (Linux): /proc/meminfo's MemAvailable, or
/// less if a cgroup (a container's memory limit) allows less.
pub fn available() -> Option<u64> {
    let system = (|| {
        let text = std::fs::read_to_string("/proc/meminfo").ok()?;
        let line = text.lines().find(|l| l.starts_with("MemAvailable:"))?;
        let kib: u64 = line.split_whitespace().nth(1)?.parse().ok()?;
        Some(kib * 1024)
    })();
    let cgroup = (|| {
        let read = |f: &str| std::fs::read_to_string(format!("/sys/fs/cgroup/{f}")).ok();
        let max: u64 = read("memory.max")?.trim().parse().ok()?;
        let used: u64 = read("memory.current")?.trim().parse().ok()?;
        Some(max.saturating_sub(used))
    })();
    match (system, cgroup) {
        (Some(a), Some(b)) => Some(a.min(b)),
        (a, b) => a.or(b),
    }
}

/// Bytes for people: B, KiB, MiB or GiB.
pub fn human(bytes: u64) -> String {
    const UNITS: [&str; 4] = ["B", "KiB", "MiB", "GiB"];
    let mut v = bytes as f64;
    let mut unit = 0;
    while v >= 1024.0 && unit < UNITS.len() - 1 {
        v /= 1024.0;
        unit += 1;
    }
    match unit {
        0 => format!("{bytes} B"),
        _ if v < 10.0 => format!("{v:.2} {}", UNITS[unit]),
        _ => format!("{v:.0} {}", UNITS[unit]),
    }
}

impl Memory {
    /// A table for people, one part a line, then the totals.
    pub fn table(&self) -> String {
        let mut s = String::new();
        for p in &self.parts {
            let _ = writeln!(
                s,
                "  {:>10}  {}{}",
                human(p.bytes),
                p.what,
                if p.transient { " (only then)" } else { "" }
            );
        }
        let _ = writeln!(s, "  {:>10}  held throughout the run", human(self.resident));
        let _ = write!(s, "  {:>10}  at most, plus a few MiB", human(self.peak));
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulation::Simulation;

    const BASE: &str = r#"
        [grid]
        points = [40, 36, 30]
        spacing = 10.0
        [time]
        duration = 0.2
        [boundary]
        pml_width = 8
        [[model.layer]]
        top = 0.0
        vp = 2000.0
        vs = 1000.0
        rho = 2000.0
        [[model.layer]]
        top = 150.0
        vp = 3500.0
        vs = 2000.0
        rho = 2400.0
        [[model.sphere]]
        center = [200.0, 180.0, 100.0]
        radius = 30.0
        vp = 1500.0
        vs = 0.0
        rho = 1000.0
        [[source]]
        position = [200.0, 180.0, 60.0]
        explosion = 1e12
        wavelet = { kind = "ricker", frequency = 10.0 }
        [[source]]
        position = [100.0, 100.0, 60.0]
        force = [0.0, 0.0, 1e9]
        wavelet = { kind = "ricker", frequency = 8.0 }
        [receivers]
        components = ["vx", "vz", "p"]
        [[receivers.line]]
        from = [100.0, 180.0, 0.0]
        to = [300.0, 180.0, 0.0]
        count = 11
        [output]
        record_every = 2
    "#;

    fn sim(text: &str) -> Simulation {
        Simulation::new(&Config::from_toml(text).unwrap(), Path::new(".")).unwrap()
    }

    /// The estimate, made without building anything, counts exactly what a run allocates.
    #[test]
    fn counts_what_a_run_allocates() {
        for extra in [
            "",
            "[scheme]\nremove_time_dispersion = false\norder = 4",
            "[[snapshot]]\nfield = \"vz\"\nplane = \"volume\"\ninterval = 0.1\nformats = [\"npy\"]",
        ] {
            let text = format!("{BASE}\n{extra}");
            let cfg = Config::from_toml(&text).unwrap();
            let m = estimate(&cfg, Path::new(".")).unwrap();
            let s = sim(&text);
            let series: usize = s
                .sources
                .iter()
                .map(|x| x.series.as_ref().map_or(0, Vec::len))
                .sum();
            let traces: usize = s.receivers.traces.iter().map(Vec::len).sum();
            let allocated = (s.solver.memory_bytes() + 8 * series + 4 * traces) as u64;
            assert_eq!(m.resident, allocated, "{extra}");
            assert_eq!(s.info.memory, m, "{extra}");
            if extra.contains("volume") {
                assert_eq!(
                    m.peak,
                    m.resident + 32 * s.solver.grid.interior_len() as u64
                );
            } else {
                assert_eq!(m.peak, m.resident);
            }
        }
    }

    #[test]
    fn grows_with_the_grid_and_reads_for_people() {
        let small = estimate(&Config::from_toml(BASE).unwrap(), Path::new(".")).unwrap();
        let big = BASE.replace("points = [40, 36, 30]", "points = [400, 360, 300]");
        let big = estimate(&Config::from_toml(&big).unwrap(), Path::new(".")).unwrap();
        // About a thousand times the nodes, though only the arrays on the grid scale so.
        assert!(big.resident > 500 * small.resident);
        assert_eq!(human(512), "512 B");
        assert_eq!(human(78 << 20), "78 MiB");
        assert_eq!(human(5 * (1 << 30) + (1 << 29)), "5.50 GiB");
        assert!(small.table().contains("wavefield"));
    }

    /// Files are read through for their extremes, which match the built model's.
    #[test]
    fn reads_model_files_in_pieces() {
        let dir = std::env::temp_dir().join(format!("vibe-rations-memory-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let (nx, ny, nz) = (20, 18, 16);
        let n = nx * ny * nz;
        let vp: Vec<f32> = (0..n).map(|i| 2000.0 + (i % 97) as f32).collect();
        let vs: Vec<f32> = (0..n)
            .map(|i| {
                if i % 50 == 0 {
                    0.0
                } else {
                    900.0 + (i % 31) as f32
                }
            })
            .collect();
        let rho = vec![2000.0f32; n];
        for (name, data) in [("vp", &vp), ("vs", &vs), ("rho", &rho)] {
            crate::npy::write(&dir.join(format!("{name}.npy")), &[nz, ny, nx], data).unwrap();
        }
        let cfg = Config::from_toml(&format!(
            r#"
            [grid]
            points = [{nx}, {ny}, {nz}]
            spacing = 10.0
            [time]
            duration = 0.05
            [boundary]
            pml_width = 4
            [model.files]
            vp = "vp.npy"
            vs = "vs.npy"
            rho = "rho.npy"
            [[source]]
            position = [100.0, 90.0, 40.0]
            explosion = 1.0
            wavelet = {{ kind = "ricker", frequency = 10.0 }}
            "#
        ))
        .unwrap();
        let grid = cfg.grid.grid().unwrap().with_halo(4);
        let model = Model::from_config(grid, &cfg.model, &dir).unwrap();
        let extremes = Model::extremes(grid, &cfg.model, &dir).unwrap();
        assert_eq!(extremes, [model.vp_range()[1], model.slowest_wave()]);
        let s = Simulation::new(&cfg, &dir).unwrap();
        assert_eq!(estimate(&cfg, &dir).unwrap(), s.info.memory);
        std::fs::remove_dir_all(dir).unwrap();
    }
}
