//! The memory a run needs, worked out from its configuration without building it, so that a
//! model too large for the machine (or the browser) is caught before anything is allocated.
//!
//! Every large array is counted exactly, from the same quantities the run uses: the wavefield
//! and material on the padded grid, the model on the interior, the CPML memory variables, the
//! sources' pre-warped time series, the traces, and the copy a volume snapshot makes while it is
//! written. A gradient run ([`crate::gradient::gradient`]) holds much more: a second solver for
//! the adjoint, the residuals, checkpoints of the forward wavefield, one replayed segment of its
//! velocities, the correlation sums and the forward and adjoint strains. What is left (the
//! program, the stencil and receiver operators, per-thread scratch rows) is a few MiB.

use std::fmt::Write as _;
use std::path::Path;

use anyhow::{Context, Result};
use serde::Serialize;

use crate::config::{Config, PlaneKind};
use crate::fd::Stencil;
use crate::gradient::checkpoint_interval;
use crate::grid::Grid;
use crate::model::Model;
use crate::receiver::Component;
use crate::simulation::timing;

/// One part of a run's memory.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct Part {
    pub what: String,
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

/// What sets an inversion's memory besides its gradient runs.
#[derive(Clone, Copy, Debug)]
pub struct InversionSize {
    pub shots: usize,
    /// The inverted parameters at their free nodes.
    pub unknowns: usize,
    /// Pairs L-BFGS keeps.
    pub memory: usize,
    pub precondition: bool,
    /// Whether a stage low-pass filters the traces, holding a filtered copy of the observed.
    pub filtered: bool,
    /// Whether it estimates the receivers' side of the illumination, held between stages.
    pub receiver_side: bool,
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
    /// Receivers recording pressure: a gradient run keeps a second copy of their residuals.
    pub pressure_traces: usize,
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
            pressure_traces: if cfg.receivers.components.contains(&Component::P) {
                receivers
            } else {
                0
            },
            samples: steps.div_ceil(cfg.output.record_every.max(1)),
            volume_snapshots: cfg.snapshots.iter().any(|s| s.plane == PlaneKind::Volume),
        })
    }

    /// The memory of a forward run.
    pub fn memory(&self) -> Memory {
        let mut parts = self.forward_parts();
        if self.volume_snapshots {
            // The nodes (three indices each), their values and their shading.
            parts.push(Part {
                what: "a volume snapshot, while it is written".into(),
                bytes: (24 + 4 + 4) * self.grid.interior_len() as u64,
                transient: true,
            });
        }
        Memory::of(parts)
    }

    /// The memory of a gradient run ([`crate::gradient::gradient`]) with a checkpoint every
    /// `checkpoint_every` steps (by default the interval that about minimises it), and of the
    /// observed traces the `gradient` command reads. It peaks while the forward run is
    /// replayed: on top of the forward run (which writes no snapshots then), the adjoint's
    /// solver and residuals, every checkpoint, one segment's velocities, the correlation sums,
    /// the illumination and the two strains. All but the sums and the illumination are freed
    /// before the gradient is put together.
    pub fn gradient_memory(&self, checkpoint_every: Option<usize>) -> Memory {
        let observed = part("observed traces", self.trace_bytes());
        Memory::of(self.gradient_parts(checkpoint_every, observed))
    }

    /// The memory of an inversion ([`crate::inversion::Inversion`]) whose shots are each like
    /// this plan's run (one source): one shot's gradient run at a time, every shot's observed
    /// traces, and the inversion's own arrays, with L-BFGS's memory full. It peaks in a gradient
    /// run after a step: the model as it stands, the start and the model stepped to, the
    /// gradient and illumination summed over the shots, and per unknown the L-BFGS pairs, the
    /// unknowns, their gradient, the step, the preconditioner and the node's index.
    pub fn inversion_memory(&self, inversion: &InversionSize) -> Memory {
        let interior = self.grid.interior_len() as u64;
        let n = inversion.unknowns as u64;
        let shots = inversion.shots;
        let copies = 1 + u64::from(inversion.filtered);
        let observed = part(
            if inversion.filtered {
                format!("observed traces of {shots} shot(s), as recorded and filtered")
            } else {
                format!("observed traces of {shots} shot(s)")
            },
            copies * shots as u64 * self.trace_bytes(),
        );
        let mut parts = self.gradient_parts(None, observed);
        parts.extend([
            part(
                "models: as it stands, the start, and the one stepped to",
                3 * 4 * 3 * interior,
            ),
            part(
                "the gradient and illumination summed over the shots",
                8 * 4 * interior,
            ),
            part(
                format!(
                    "L-BFGS: {} pairs of changes of {n} unknowns",
                    inversion.memory
                ),
                16 * inversion.memory as u64 * n,
            ),
            part(
                "the unknowns, their gradient, the step and the preconditioner",
                8 * n * (3 + u64::from(inversion.precondition)),
            ),
            part(
                "the free nodes' indices",
                std::mem::size_of::<usize>() as u64 * n,
            ),
        ]);
        if inversion.receiver_side {
            // Before the free nodes' indices, which a test expects last.
            let at = parts.len() - 1;
            parts.insert(
                at,
                part("the receivers' side of the illumination", 8 * interior),
            );
        }
        Memory::of(parts)
    }

    /// A gradient run's parts, with `observed` traces.
    fn gradient_parts(&self, checkpoint_every: Option<usize>, observed: Part) -> Vec<Part> {
        let padded = self.grid.padded_len() as u64;
        let interior = self.grid.interior_len() as u64;
        let every = checkpoint_interval(self.steps, checkpoint_every);
        let checkpoints = self.steps.div_ceil(every);
        // The segment holds the velocities before its first step and after each of the others.
        let segment = if self.steps == 0 { 0 } else { every + 1 };
        let mut parts = self.forward_parts();
        parts.extend([
            part("the adjoint's solver, the same again", self.solver_bytes()),
            observed,
            part("residuals", self.trace_bytes()),
            part(
                "pressure residuals, as they enter the adjoint strain",
                4 * (self.pressure_traces * self.samples) as u64,
            ),
            part(
                format!(
                    "checkpoints: {checkpoints}, every {every} steps, of the wavefield and \
                     layers' memory"
                ),
                checkpoints as u64 * 4 * (9 * padded + self.psi()),
            ),
            part(
                format!("a replayed segment: {segment} copies of the velocities"),
                segment as u64 * 4 * 3 * padded,
            ),
            part(
                "correlation sums: 8 per padded node, in double precision",
                8 * 8 * padded,
            ),
            part(
                "illumination: Σ|v|² per node, in double precision",
                8 * interior,
            ),
            part(
                "forward and adjoint strains: 12 arrays on the padded grid",
                4 * 12 * padded,
            ),
        ]);
        parts
    }

    fn trace_bytes(&self) -> u64 {
        4 * (self.traces * self.samples) as u64
    }

    /// The parts of a forward run held throughout: the solver's, the sources' and the traces.
    fn forward_parts(&self) -> Vec<Part> {
        let g = &self.grid;
        let padded = g.padded_len() as u64;
        let interior = g.interior_len() as u64;
        vec![
            part("wavefield: 9 arrays on the padded grid", 4 * 9 * padded),
            part("material: 8 arrays on the padded grid", 4 * 8 * padded),
            part("model: vp, vs and rho", 4 * 3 * interior),
            part(
                "absorbing layers' memory",
                4 * (self.psi() + self.sponge_factors()),
            ),
            part(
                "sources' pre-warped time series",
                8 * (self.warped_sources * self.steps) as u64,
            ),
            part("traces", self.trace_bytes()),
        ]
    }

    /// Bytes of a solver (`Solver::memory_bytes`).
    fn solver_bytes(&self) -> u64 {
        let padded = self.grid.padded_len() as u64;
        let interior = self.grid.interior_len() as u64;
        4 * ((9 + 8) * padded + 3 * interior + self.psi() + self.sponge_factors())
    }

    /// The CPML memory variables: three per slab node, for velocities and for stresses.
    fn psi(&self) -> u64 {
        let [nx, ny, nz] = self.grid.n.map(|n| n as u64);
        let [w0, w1, w2] = self.layer_nodes.map(|w| w as u64);
        2 * 3 * (nz * ny * w0 + nz * w1 * nx + w2 * ny * nx)
    }

    /// The sponge's damping factors, one per node along each axis.
    fn sponge_factors(&self) -> u64 {
        let [nx, ny, nz] = self.grid.n.map(|n| n as u64);
        if self.sponge { nx + ny + nz } else { 0 }
    }
}

fn part(what: impl Into<String>, bytes: u64) -> Part {
    Part {
        what: what.into(),
        bytes,
        transient: false,
    }
}

/// The memory the run `cfg` describes would need, without building it: the model's extremes
/// (which set the time step, and so the number of steps and samples) come from its layers and
/// spheres, or from its files read through a piece at a time. `base` is the directory relative
/// paths start from.
pub fn estimate(cfg: &Config, base: &Path) -> Result<Memory> {
    Ok(plan(cfg, base)?.memory())
}

/// The memory a gradient run of `cfg` would need ([`Plan::gradient_memory`]), worked out like
/// [`estimate`].
pub fn estimate_gradient(
    cfg: &Config,
    base: &Path,
    checkpoint_every: Option<usize>,
) -> Result<Memory> {
    Ok(plan(cfg, base)?.gradient_memory(checkpoint_every))
}

/// The memory an inversion of `cfg` (`invert`) would need ([`Plan::inversion_memory`]),
/// worked out like [`estimate`]. Each source is a shot, run with the time step all of them
/// set. Unknowns are counted at every node outside the layers and below `fixed_above`; vs isn't
/// inverted where it starts at zero, so with fluids in the model this is a little more than
/// the inversion needs.
pub fn estimate_inversion(cfg: &Config, base: &Path) -> Result<Memory> {
    let settings = cfg
        .inversion
        .as_ref()
        .context("the configuration has no [inversion] section")?;
    let (grid, steps) = timed(cfg, base)?;
    let size = InversionSize {
        shots: cfg.sources.len(),
        unknowns: settings.parameters.len() * free_nodes(cfg, &grid, settings.fixed_above),
        memory: settings.memory,
        precondition: settings.precondition || !settings.scale.is_unit(),
        filtered: settings.stages()?.iter().any(|s| s.max_frequency.is_some()),
        receiver_side: settings.precondition && settings.receiver_runs > 0,
    };
    Ok(Plan::new(&cfg.shot(0)?, grid, steps)?.inversion_memory(&size))
}

fn plan(cfg: &Config, base: &Path) -> Result<Plan> {
    let (grid, steps) = timed(cfg, base)?;
    Plan::new(cfg, grid, steps)
}

/// The grid with its halo, and the steps of the run.
fn timed(cfg: &Config, base: &Path) -> Result<(Grid, usize)> {
    let taylor = Stencil::taylor(cfg.scheme.order)?;
    let grid = cfg.grid.grid()?.with_halo(taylor.half_length());
    let [vp_max, slowest] = Model::extremes(grid, &cfg.model, base)?;
    Ok((grid, timing(cfg, &grid, vp_max, slowest)?.steps))
}

/// Nodes an inversion may change (`Inversion::new`): outside the absorbing layers and not
/// above `fixed_above`.
fn free_nodes(cfg: &Config, grid: &Grid, fixed_above: Option<f64>) -> usize {
    let b = &cfg.boundary;
    let w = b.pml_width;
    let [nx, ny, nz] = grid.n;
    let top = if b.free_surface { 0 } else { w };
    let fixed = fixed_above.unwrap_or(f64::NEG_INFINITY);
    let depths = (top..nz.saturating_sub(w))
        .filter(|&k| k as f64 * grid.h[2] >= fixed)
        .count();
    nx.saturating_sub(2 * w) * ny.saturating_sub(2 * w) * depths
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
    /// The totals of `parts`, without those that hold nothing.
    fn of(mut parts: Vec<Part>) -> Self {
        parts.retain(|p| p.bytes > 0);
        let resident = parts.iter().filter(|p| !p.transient).map(|p| p.bytes).sum();
        let peak = resident
            + parts
                .iter()
                .filter(|p| p.transient)
                .map(|p| p.bytes)
                .max()
                .unwrap_or(0);
        Self {
            parts,
            resident,
            peak,
        }
    }

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
        if self.peak > self.resident {
            let _ = writeln!(s, "  {:>10}  held throughout the run", human(self.resident));
        }
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

    /// A gradient run's plan counts what `gradient` holds at its peak, with the observed
    /// traces: with and without pressure, the transform, a free surface (and its sponge), and
    /// with the default interval and one that leaves a short last segment.
    #[test]
    fn counts_what_a_gradient_run_allocates() {
        let free_surface = BASE.replace("pml_width = 8", "pml_width = 8\nfree_surface = true");
        let velocity_only = BASE.replace(r#"["vx", "vz", "p"]"#, r#"["vx", "vz"]"#);
        for (text, every) in [
            (BASE.to_string(), None),
            (
                format!("{BASE}\n[scheme]\nremove_time_dispersion = false\norder = 4"),
                None,
            ),
            (free_surface, Some(7)),
            (velocity_only, Some(1000)),
        ] {
            let cfg = Config::from_toml(&text).unwrap();
            let m = estimate_gradient(&cfg, Path::new("."), every).unwrap();
            let mut s = sim(&text);
            let observed = s.receivers.traces.clone();
            let weights = vec![1.0; observed.len()];
            let run = crate::gradient::gradient(
                &mut s,
                |r| crate::gradient::l2_misfit(&r.traces, &observed, &weights),
                every,
            )
            .unwrap();
            let held = run.memory_bytes + 4 * observed.iter().map(Vec::len).sum::<usize>();
            assert_eq!(m.peak, held as u64, "{every:?}\n{}", m.table());
            assert_eq!(m.peak, m.resident);
            let forward = estimate(&cfg, Path::new(".")).unwrap();
            assert!(m.peak > 3 * forward.peak);
        }
    }

    /// An inversion's plan counts what it holds at its peak, a gradient run after a step with
    /// L-BFGS's memory full, from its configuration alone: two shots recording velocity and
    /// pressure, the top held fixed, the gradient preconditioned.
    #[test]
    fn counts_what_an_inversion_allocates() {
        use crate::inversion::{Inversion, Observed};
        let model = |vs: f64, stages: &str| {
            format!(
                r#"
                [grid]
                points = [30, 28, 24]
                spacing = 10.0
                [time]
                duration = 0.2
                [boundary]
                pml_width = 6
                [output]
                record_every = 2
                [[model.layer]]
                top = 0.0
                vp = 2200.0
                vs = 1200.0
                rho = 2000.0
                [[model.sphere]]
                center = [150.0, 140.0, 120.0]
                radius = 30.0
                vp = 2200.0
                vs = {vs}
                rho = 2000.0
                [[source]]
                position = [90.0, 140.0, 70.0]
                force = [0.0, 0.0, 1e9]
                wavelet = {{ kind = "ricker", frequency = 14.0 }}
                [[source]]
                position = [210.0, 120.0, 70.0]
                explosion = 1e12
                wavelet = {{ kind = "ricker", frequency = 14.0 }}
                [receivers]
                components = ["vx", "vz", "p"]
                [[receivers.line]]
                from = [80.0, 140.0, 100.0]
                to = [220.0, 140.0, 100.0]
                count = 8
                [inversion]
                observed = "unused"
                memory = 1
                fixed_above = 85.0
                {stages}
                "#
            )
        };
        let truth = Config::from_toml(&model(1350.0, "")).unwrap();
        let observed: Vec<Observed> = (0..2)
            .map(|s| {
                let mut sim = Simulation::new(&truth.shot(s).unwrap(), Path::new(".")).unwrap();
                sim.run(None, false).unwrap();
                sim.receivers.traces.into()
            })
            .collect();
        // As recorded, with a stage that filters (and so keeps filtered copies), and with the
        // receivers' side of the illumination and a scale.
        let filtering = "[[inversion.stage]]\nmax_frequency = 10.0\niterations = 2";
        let receivers = "receiver_runs = 1\n[inversion.scale]\nvp = 0.5";
        for stages in ["", filtering, receivers] {
            let cfg = Config::from_toml(&model(1200.0, stages)).unwrap();
            let m = estimate_inversion(&cfg, Path::new(".")).unwrap();
            let settings = cfg.inversion.clone().unwrap();
            let start =
                Model::from_config(Simulation::grid(&cfg).unwrap(), &cfg.model, Path::new("."))
                    .unwrap();
            let mut inversion = Inversion::new(&cfg, start, observed.clone(), &settings).unwrap();
            assert_eq!(
                m.parts.last().unwrap().bytes as usize,
                8 * inversion.unknowns(),
                "the free nodes"
            );
            let first = inversion.iterate().unwrap().unwrap();
            let second = inversion.iterate().unwrap().unwrap();
            // The second gradient run had L-BFGS's one pair to hold.
            assert!(
                first.memory == 1 && !second.restarted,
                "{stages}: {first:?} {second:?}"
            );
            assert_eq!(
                m.peak,
                inversion.memory_bytes as u64,
                "{stages}\n{}",
                m.table()
            );
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
