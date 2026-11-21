//! A run: the solver with its sources, receivers and snapshots, set up from a configuration.

use std::f64::consts::PI;
use std::fmt::Write as _;
use std::fs;
use std::path::Path;
use std::time::{Duration, Instant};

use anyhow::{Context, Result, bail, ensure};
use rayon::prelude::*;
use serde::Serialize;

use crate::config::{Coefficients, Config, Format, PlaneKind, default_courant};
use crate::cpml::{PmlSettings, Stagger};
use crate::fd::Stencil;
use crate::grid::Grid;
use crate::image;
use crate::memory::{Memory, Plan, human};
use crate::model::{Files, Model};
use crate::npy;
use crate::receiver::{Component, Receivers};
use crate::snapshot::Snapshot;
use crate::solver::{Band, Solver};
use crate::source::{Mechanism, Source, Wavelet};
use crate::tdt::Warp;

/// Fewer grid points than this per shortest wavelength and numerical dispersion shows: a few
/// percent of misfit after ten wavelengths at the default time step (README, "Accuracy").
pub fn min_points_per_wavelength(order: usize, coefficients: Coefficients) -> f64 {
    match (order, coefficients) {
        (2, _) => 10.0,
        (4, _) => 5.0,
        (6, _) => 4.5,
        (_, Coefficients::Taylor) => 4.0,
        (8, Coefficients::Optimised) => 3.5,
        (_, Coefficients::Optimised) => 3.0,
    }
}

/// The largest stable time step: the Courant–Friedrichs–Lewy limit of the staggered scheme,
/// dt ≤ 1 / (vp_max Σ|cₙ| √(1/dx² + 1/dy² + 1/dz²)).
pub fn stable_dt(grid: &Grid, vp_max: f64, stencil: &Stencil) -> f64 {
    let s: f64 = grid.h.iter().map(|h| 1.0 / (h * h)).sum();
    1.0 / (vp_max * stencil.sum_abs() * s.sqrt())
}

pub struct Simulation {
    pub solver: Solver,
    pub sources: Vec<Source>,
    pub receivers: Receivers,
    pub snapshots: Vec<Snapshot>,
    /// Steps to take.
    pub steps: usize,
    /// Record the receivers every this many steps.
    pub record_every: usize,
    pub info: RunInfo,
    /// Whether the traces have been un-warped ([`Simulation::finish`]).
    finished: bool,
}

/// The largest phase error per node optimised coefficients may have at the edge of their
/// band (`Stencil::optimised`). 1e-3 let waves inside the band, P waves in particular, do
/// worse than with Taylor's coefficients (research log).
pub const PHASE_TOLERANCE: f64 = 1e-4;

/// The time-dispersion transforms keep frequencies up to this multiple of the sources' highest
/// (2.5 f0): a Ricker wavelet's spectrum there is 1e-9 of its peak.
const TRANSFORM_BAND: f64 = 2.0;

/// What a run is set up to do; written to `run.json` alongside the results.
#[derive(Clone, Debug, Serialize)]
pub struct RunInfo {
    pub points: [usize; 3],
    pub spacing: [f64; 3],
    pub extent: [f64; 3],
    pub vp_range: [f64; 2],
    pub vs_range: [f64; 2],
    pub rho_range: [f64; 2],
    pub dt: f64,
    pub dt_stable: f64,
    pub steps: usize,
    pub duration: f64,
    pub max_frequency: f64,
    pub slowest_wave: f64,
    pub points_per_wavelength: f64,
    pub pml_width: usize,
    pub free_surface: bool,
    /// Order of the spatial derivatives.
    pub order: usize,
    /// Whether leapfrog's time error is removed from the traces (crate::tdt).
    pub remove_time_dispersion: bool,
    /// How the derivative coefficients were chosen, and for optimised ones the largest k·h
    /// they were designed for.
    pub coefficients: Coefficients,
    pub band: Option<f64>,
    /// The most memory the run holds at once (`memory.peak`).
    pub memory_bytes: usize,
    /// What it holds, part by part (crate::memory).
    pub memory: Memory,
    pub sources: Vec<SourceInfo>,
    pub receivers: usize,
    pub components: Vec<String>,
    /// Time of the first receiver sample; the velocities lag the stresses by dt/2.
    pub first_sample_time: f64,
    pub sample_interval: f64,
    pub samples: usize,
    pub snapshots: usize,
    pub warnings: Vec<String>,
}

#[derive(Clone, Debug, Serialize)]
pub struct SourceInfo {
    /// Where the source acts, after snapping to a node.
    pub position: [f64; 3],
    pub mechanism: String,
    pub dominant_frequency: f64,
}

/// What follows from a configuration and its model's extremes before anything is allocated:
/// the stencil, the time step and the number of steps. [`Simulation::new`] and
/// [`crate::memory::estimate`] share it.
pub(crate) struct Timing {
    pub stencil: Stencil,
    pub dt: f64,
    pub dt_stable: f64,
    pub steps: usize,
    pub ppw: f64,
    pub slowest: f64,
    pub f_max: f64,
    pub f_dominant: f64,
    pub wavelets: Vec<Wavelet>,
}

/// The run's timing on `grid` (with its halo) for a model whose fastest P speed is `vp_max`
/// and slowest wave `slowest`.
pub(crate) fn timing(cfg: &Config, grid: &Grid, vp_max: f64, slowest: f64) -> Result<Timing> {
    ensure!(
        !cfg.sources.is_empty(),
        "the configuration has no [[source]]"
    );
    let wavelets = cfg
        .sources
        .iter()
        .map(|s| s.wavelet.wavelet())
        .collect::<Result<Vec<_>>>()?;
    let f_dominant = wavelets
        .iter()
        .map(|w| w.dominant_frequency())
        .fold(0.0, f64::max);
    let f_max = wavelets
        .iter()
        .map(|w| w.max_frequency())
        .fold(0.0, f64::max);

    let h_max = grid.h.iter().copied().fold(0.0, f64::max);
    let ppw = slowest / f_max / h_max;
    let stencil = match cfg.scheme.coefficients {
        Coefficients::Taylor => Stencil::taylor(cfg.scheme.order)?,
        // The grid can't hold waves shorter than two nodes anyway.
        Coefficients::Optimised => {
            if let Some(p) = cfg.scheme.optimised_for {
                ensure!(
                    p > 2.0,
                    "optimised_for must be more than 2 points per wavelength, not {p}"
                );
            }
            let design = cfg.scheme.optimised_for.unwrap_or(ppw);
            Stencil::optimised(
                cfg.scheme.order,
                (2.0 * PI / design).min(0.9 * PI),
                PHASE_TOLERANCE,
            )?
        }
    };

    let dt_stable = stable_dt(grid, vp_max, &stencil);
    let dt = match cfg.time.dt {
        Some(dt) => {
            ensure!(dt > 0.0, "dt must be positive");
            ensure!(
                dt <= dt_stable,
                "dt = {dt:e} s is unstable; the limit for this grid and model is {dt_stable:e} s"
            );
            dt
        }
        None => {
            let c = cfg.time.courant.unwrap_or_else(|| {
                default_courant(cfg.scheme.order, cfg.scheme.remove_time_dispersion)
            });
            ensure!(c > 0.0 && c <= 1.0, "courant must be in (0, 1]");
            c * dt_stable
        }
    };
    ensure!(cfg.time.duration > 0.0, "duration must be positive");
    let steps = (cfg.time.duration / dt).ceil() as usize;
    Ok(Timing {
        stencil,
        dt,
        dt_stable,
        steps,
        ppw,
        slowest,
        f_max,
        f_dominant,
        wavelets,
    })
}

impl Simulation {
    /// Sets up a run; `base` is the directory relative paths in the configuration start from.
    pub fn new(cfg: &Config, base: &Path) -> Result<Self> {
        Self::with_files(cfg, Files::Dir(base))
    }

    /// Sets up a run with the model's files from `files`, e.g. held in memory by a browser.
    pub fn with_files(cfg: &Config, files: Files) -> Result<Self> {
        // from_files validates the model.
        let model = Model::from_files(Self::grid(cfg)?, &cfg.model, files)?;
        Self::set_up(cfg, model)
    }

    /// The grid a configuration describes, with the halo its stencil needs: what a model
    /// for [`Simulation::with_model`] must be on.
    pub fn grid(cfg: &Config) -> Result<Grid> {
        let stencil = Stencil::taylor(cfg.scheme.order)?;
        Ok(cfg.grid.grid()?.with_halo(stencil.half_length()))
    }

    /// Sets up a run of `model` instead of the configuration's own, e.g. a model being
    /// inverted for. Everything else, the time step included, comes from the configuration
    /// and the model as for [`Simulation::new`].
    pub fn with_model(cfg: &Config, model: Model) -> Result<Self> {
        let grid = Self::grid(cfg)?;
        ensure!(
            model.grid == grid,
            "the model's grid ({:?} nodes spaced {:?} m, halo {}) isn't the configuration's \
             ({:?} spaced {:?}, halo {})",
            model.grid.n,
            model.grid.h,
            model.grid.halo,
            grid.n,
            grid.h,
            grid.halo
        );
        model.validate()?;
        let mut sim = Self::set_up(cfg, model)?;
        if cfg.scheme.coefficients == Coefficients::Optimised && cfg.scheme.optimised_for.is_none()
        {
            sim.info.warnings.push(
                "the optimised coefficients follow this model's slowest wave, so each model gets \
                 its own operator; pin [scheme] optimised_for to keep one across an inversion \
                 or a finite-difference check"
                    .into(),
            );
        }
        Ok(sim)
    }

    /// The rest of [`Simulation::new`], for a validated model on the configuration's grid.
    fn set_up(cfg: &Config, model: Model) -> Result<Self> {
        let grid = model.grid;
        let mut warnings = Vec::new();
        let [_, vp_max] = model.vp_range();
        let Timing {
            stencil,
            dt,
            dt_stable,
            steps,
            ppw,
            slowest,
            f_max,
            f_dominant,
            wavelets,
        } = timing(cfg, &grid, vp_max, model.slowest_wave())?;

        let min_ppw = min_points_per_wavelength(cfg.scheme.order, cfg.scheme.coefficients);
        if ppw < min_ppw {
            warnings.push(format!(
                "only {ppw:.1} grid points per shortest wavelength ({slowest} m/s at {f_max:.1} Hz); \
                 at order {} {min_ppw} or more keeps numerical dispersion small",
                cfg.scheme.order
            ));
        }

        let b = &cfg.boundary;
        let w = b.pml_width;
        for (a, axis) in ["x", "y", "z"].iter().enumerate() {
            let layers = if a == 2 && b.free_surface { w } else { 2 * w };
            ensure!(
                grid.n[a] > layers + 4,
                "{} nodes along {axis} leave no room inside absorbing layers {w} nodes thick",
                grid.n[a]
            );
        }
        if w > 0 && w < 8 {
            warnings.push(format!(
                "absorbing layers {w} nodes thick will reflect noticeably; 10 or more is better"
            ));
        }
        ensure!(b.pml_sponge >= 0.0, "pml_sponge can't be negative");
        let pml_vp = b.pml_vp_max.unwrap_or(vp_max);
        ensure!(pml_vp > 0.0, "pml_vp_max must be positive");
        let pml = PmlSettings {
            power: b.pml_power,
            reflection: b.pml_reflection,
            kappa_max: b.pml_kappa_max,
            alpha_max: PI * b.pml_alpha_frequency.unwrap_or(f_dominant),
            vp_max: pml_vp,
            sponge: b.pml_sponge,
        };
        let solver = Solver::new(model, dt, w, b.free_surface, &pml).with_stencil(stencil);
        let in_pml = |node: [usize; 3]| {
            node.iter()
                .zip(&solver.cpml)
                .any(|(&n, c)| c.slab(n).is_some())
        };

        let mut sources = Vec::new();
        let mut source_info = Vec::new();
        for (s, wavelet) in cfg.sources.iter().zip(wavelets) {
            let node = grid
                .nearest_node(s.position)
                .with_context(|| format!("source at {:?} is off the grid", s.position))?;
            if in_pml(node) {
                warnings.push(format!(
                    "source at {:?} is inside an absorbing layer",
                    s.position
                ));
            }
            let mechanism = s.mechanism()?;
            source_info.push(SourceInfo {
                position: grid.position(node),
                mechanism: match mechanism {
                    Mechanism::Moment(m) => format!("moment tensor {m:?}"),
                    Mechanism::Force(f) => format!("force {f:?}"),
                },
                dominant_frequency: wavelet.dominant_frequency(),
            });
            let mut source = Source::new(&solver, node, mechanism, wavelet);
            if cfg.scheme.remove_time_dispersion {
                source.remove_time_dispersion(dt, steps, TRANSFORM_BAND * f_max);
            }
            sources.push(source);
        }

        let record_every = cfg.output.record_every;
        ensure!(record_every >= 1, "record_every must be at least 1");
        let samples = steps.div_ceil(record_every);
        let mut nodes = Vec::new();
        for p in cfg.receivers.positions()? {
            let node = grid
                .nearest_node(p)
                .with_context(|| format!("receiver at {p:?} is off the grid"))?;
            nodes.push(node);
        }
        let in_layers = nodes.iter().filter(|&&n| in_pml(n)).count();
        if in_layers > 0 {
            warnings.push(format!(
                "{in_layers} receiver(s) inside the absorbing layers will record damped waves"
            ));
        }
        ensure!(
            !cfg.receivers.components.is_empty() || nodes.is_empty(),
            "receivers need at least one component"
        );
        // A second copy records nothing new, and the adjoint would back-propagate both.
        for (i, c) in cfg.receivers.components.iter().enumerate() {
            ensure!(
                !cfg.receivers.components[..i].contains(c),
                "the receivers list {} twice",
                c.name()
            );
        }
        let receivers = Receivers::new(&solver, nodes, cfg.receivers.components.clone(), samples);

        let mut snapshots = Vec::new();
        for s in &cfg.snapshots {
            ensure!(s.interval > 0.0, "a snapshot interval must be positive");
            let index = match s.plane {
                PlaneKind::Xy => Some(2),
                PlaneKind::Xz => Some(1),
                PlaneKind::Yz => Some(0),
                PlaneKind::Volume => None,
            }
            .map(|a| {
                let r = (s.at / grid.h[a]).round();
                ensure!(
                    (0.0..grid.n[a] as f64).contains(&r),
                    "snapshot plane at {} m is off the grid",
                    s.at
                );
                Ok(r as usize)
            })
            .transpose()?
            .unwrap_or(0);
            ensure!(
                !s.formats.is_empty(),
                "a snapshot needs at least one format"
            );
            snapshots.push(Snapshot {
                component: s.field,
                plane: s.plane,
                index,
                every: ((s.interval / dt).round() as usize).max(1),
                png: s.formats.contains(&Format::Png),
                npy: s.formats.contains(&Format::Npy),
                clip: s.clip.map(|c| c as f32),
            });
        }

        let memory = Plan::new(cfg, grid, steps)?.memory();
        let m = &solver.model;
        let info = RunInfo {
            points: grid.n,
            spacing: grid.h,
            extent: grid.extent(),
            vp_range: m.vp_range(),
            vs_range: m.vs_range(),
            rho_range: m.rho_range(),
            dt,
            dt_stable,
            steps,
            duration: steps as f64 * dt,
            max_frequency: f_max,
            slowest_wave: slowest,
            points_per_wavelength: ppw,
            pml_width: w,
            free_surface: b.free_surface,
            order: cfg.scheme.order,
            remove_time_dispersion: cfg.scheme.remove_time_dispersion,
            coefficients: cfg.scheme.coefficients,
            band: solver.stencil.band(),
            memory_bytes: memory.peak as usize,
            memory,
            sources: source_info,
            receivers: receivers.nodes.len(),
            components: receivers
                .components
                .iter()
                .map(|c| c.name().to_string())
                .collect(),
            first_sample_time: 0.5 * dt,
            sample_interval: record_every as f64 * dt,
            samples,
            snapshots: snapshots.len(),
            warnings,
        };
        Ok(Self {
            solver,
            sources,
            receivers,
            snapshots,
            steps,
            record_every,
            info,
            finished: false,
        })
    }

    /// Receivers outside the absorbing layers that still read a velocity the layers damp. A
    /// velocity is interpolated from staggered values either side of its node, so on the first
    /// node outside a layer it reads one inside. That hardly matters to the traces, but the
    /// adjoint injects the residual there, where the scheme is not self-adjoint, so adjoints
    /// and gradients are no longer exact: a node further in, they are.
    pub fn receivers_next_to_layers(&self) -> usize {
        let s = &self.solver;
        let in_layers =
            |node: &[usize; 3]| node.iter().zip(&s.cpml).any(|(&n, c)| c.slab(n).is_some());
        let damped = |axis: usize, m: usize| {
            let c = &s.cpml[axis];
            c.slab(m)
                .is_some_and(|i| c.profile(Stagger::Half).d[i] > 0.0)
        };
        let r = &self.receivers;
        r.nodes
            .iter()
            .filter(|node| !in_layers(node))
            .filter(|node| {
                r.components.iter().any(|c| {
                    let axis = match c {
                        Component::Vx => 0,
                        Component::Vy => 1,
                        Component::Vz => 2,
                        Component::P => return false,
                    };
                    s.stagger_weights(axis, node[axis], 1.0)
                        .iter()
                        .any(|&(m, _)| damped(axis, m))
                })
            })
            .count()
    }

    /// The band a step back reads ([`Band`]): around the nodes outside the absorbing layers,
    /// as deep as the stencil reaches.
    pub fn band(&self) -> Band {
        let s = &self.solver;
        Band::new(
            &s.grid,
            &s.interior(),
            s.free_surface,
            s.stencil.half_length(),
        )
    }

    /// Takes the last step back ([`Solver::step_back_with`]), withdrawing what the sources
    /// added during it; `before` is `band`'s record of the state before that step.
    pub fn step_back(&mut self, band: &Band, before: &[f32]) {
        let sources = &self.sources;
        let n = self.solver.step.saturating_sub(1);
        self.solver.step_back_with(
            band,
            before,
            |f, t| {
                sources
                    .iter()
                    .filter(|s| s.is_force())
                    .for_each(|s| s.withdraw(f, n, t))
            },
            |f, t| {
                sources
                    .iter()
                    .filter(|s| !s.is_force())
                    .for_each(|s| s.withdraw(f, n, t))
            },
        );
    }

    /// Advances the wavefield one step, with the sources, without recording anything.
    pub fn step(&mut self) {
        let sources = &self.sources;
        let n = self.solver.step;
        self.solver.step_with(
            |f, t| {
                sources
                    .iter()
                    .filter(|s| s.is_force())
                    .for_each(|s| s.inject(f, n, t))
            },
            |f, t| {
                sources
                    .iter()
                    .filter(|s| !s.is_force())
                    .for_each(|s| s.inject(f, n, t))
            },
        );
    }

    /// Advances one step and records the receivers if this step is one to record.
    pub fn advance(&mut self) {
        let n = self.solver.step;
        self.step();
        if n.is_multiple_of(self.record_every) {
            self.receivers.record(&self.solver);
        }
    }

    /// Advances one step like [`Simulation::advance`], then writes the snapshots that fall
    /// due into `dir`, which must exist.
    pub fn advance_writing(&mut self, dir: &Path) -> Result<()> {
        self.advance();
        for s in &self.snapshots {
            if s.due(self.solver.step) {
                s.write(dir, &self.solver)?;
            }
        }
        Ok(())
    }

    /// Steps still to take.
    pub fn remaining(&self) -> usize {
        self.steps.saturating_sub(self.solver.step)
    }

    /// The transform that un-warps the traces of `component` (crate::tdt), when the run
    /// removes time dispersion. The misfit's derivative with respect to the recorded samples is
    /// its transpose applied to the derivative with respect to the un-warped ones.
    pub fn trace_warp(&self, component: Component) -> Option<Warp> {
        if !self.info.remove_time_dispersion {
            return None;
        }
        let dt = self.solver.dt;
        // A sample follows its step: velocities from half a step before its end, stresses
        // from its end.
        let first = match component {
            Component::P => dt,
            _ => 0.5 * dt,
        };
        let f_max = TRANSFORM_BAND * self.info.max_frequency;
        Some(Warp::from_leapfrog(
            dt,
            first,
            self.record_every as f64 * dt,
            self.receivers.samples,
            f_max,
        ))
    }

    /// After the last step, un-warps the traces if the run removes time dispersion. Only the
    /// first call after the last step does anything; [`Simulation::run`] makes it.
    pub fn finish(&mut self) {
        if self.finished || self.remaining() > 0 {
            return;
        }
        self.finished = true;
        let samples = self.receivers.samples;
        let components = self.receivers.components.clone();
        for (c, data) in components.into_iter().zip(0..) {
            let Some(warp) = self.trace_warp(c) else {
                continue;
            };
            self.receivers.traces[data]
                .par_chunks_mut(samples)
                .for_each(|trace| {
                    let x: Vec<f64> = trace.iter().map(|&v| v as f64).collect();
                    for (t, y) in trace.iter_mut().zip(warp.apply(&x)) {
                        *t = y as f32;
                    }
                });
        }
    }

    /// Bytes the run holds: its solver's, the sources' pre-warped series and the traces
    /// (the arrays [`crate::memory::Plan`] counts beforehand).
    pub fn memory_bytes(&self) -> usize {
        let series: usize = self
            .sources
            .iter()
            .map(|s| s.series.as_ref().map_or(0, Vec::len))
            .sum();
        let traces: usize = self.receivers.traces.iter().map(Vec::len).sum();
        self.solver.memory_bytes() + 8 * series + 4 * traces
    }

    /// Whether the traces have been un-warped ([`Simulation::finish`] with
    /// `remove_time_dispersion`), so that they differ from the samples the solver recorded.
    pub fn traces_unwarped(&self) -> bool {
        self.finished && self.info.remove_time_dispersion
    }

    /// The derivative of a misfit with respect to the samples the solver recorded, given its
    /// derivative `residuals` with respect to the traces as they are now (in the layout of
    /// [`Receivers::traces`]). Once the traces are un-warped that is the un-warping's
    /// transpose applied to each trace (`Warp::apply_transpose`); before, the residuals
    /// themselves. The adjoint run ([`crate::adjoint::Adjoint`]) takes the former.
    pub fn recorded_residuals(&self, residuals: Vec<Vec<f32>>) -> Vec<Vec<f32>> {
        if !self.traces_unwarped() {
            return residuals;
        }
        let samples = self.receivers.samples;
        self.receivers
            .components
            .iter()
            .zip(residuals)
            .map(|(&c, mut r)| {
                if let Some(warp) = self.trace_warp(c) {
                    r.par_chunks_mut(samples).for_each(|trace| {
                        let y: Vec<f64> = trace.iter().map(|&v| v as f64).collect();
                        for (t, x) in trace.iter_mut().zip(warp.apply_transpose(&y)) {
                            *t = x as f32;
                        }
                    });
                }
                r
            })
            .collect()
    }

    /// Runs every step, writing snapshots into `out` if given, and reporting progress on
    /// stderr if `progress`.
    pub fn run(&mut self, out: Option<&Path>, progress: bool) -> Result<()> {
        let snap_dir = out.map(|d| d.join("snapshots"));
        if let Some(dir) = &snap_dir
            && !self.snapshots.is_empty()
        {
            fs::create_dir_all(dir).with_context(|| format!("creating {}", dir.display()))?;
        }
        // The clock only for the reports: WebAssembly has none, and the browser and the run
        // server's inversions run quietly.
        let start = progress.then(Instant::now);
        let mut last_report = start;
        let cells = self.solver.grid.interior_len() as f64;
        for n in 0..self.steps {
            match &snap_dir {
                Some(dir) => self.advance_writing(dir)?,
                None => self.advance(),
            }
            let done = n + 1;
            let now = start.map(|_| Instant::now());
            let report = match (now, last_report) {
                (Some(now), Some(last)) => {
                    now - last >= Duration::from_secs(2) || done == self.steps
                }
                _ => false,
            };
            if report || done % 100 == 0 || done == self.steps {
                let peak = self.solver.peak_velocity();
                if !peak.is_finite() {
                    bail!(
                        "the wavefield blew up by step {done}; check the model, or lower dt or courant"
                    );
                }
                if let (true, Some(now), Some(start)) = (report, now, start) {
                    last_report = Some(now);
                    let elapsed = (now - start).as_secs_f64();
                    let eta = elapsed / done as f64 * (self.steps - done) as f64;
                    eprintln!(
                        "step {done:>6}/{}  t {:.4} s  peak |v| {peak:.3e} m/s  {:.0} Mcells/s  {} left",
                        self.steps,
                        done as f64 * self.solver.dt,
                        cells * done as f64 / elapsed / 1e6,
                        human_time(eta),
                    );
                }
            }
        }
        self.finish();
        Ok(())
    }

    /// Times of the receiver samples.
    pub fn sample_times(&self) -> Vec<f64> {
        (0..self.receivers.samples)
            .map(|m| (m * self.record_every) as f64 * self.solver.dt + 0.5 * self.solver.dt)
            .collect()
    }

    /// Writes the traces (`<component>.npy`, receiver-major), their times (`time.npy`), the
    /// receivers' positions (`receivers.csv`), gathers (`<component>_gather.png`) and the
    /// run's settings (`run.json`) into `dir`.
    pub fn write_outputs(&self, dir: &Path) -> Result<()> {
        fs::create_dir_all(dir).with_context(|| format!("creating {}", dir.display()))?;
        let r = &self.receivers;
        if !r.nodes.is_empty() {
            for (c, data) in r.components.iter().zip(&r.traces) {
                npy::write(
                    &dir.join(format!("{}.npy", c.name())),
                    &[r.nodes.len(), r.samples],
                    data,
                )?;
                if r.nodes.len() >= 2 {
                    image::write_gather(
                        &dir.join(format!("{}_gather.png", c.name())),
                        r.nodes.len(),
                        r.samples,
                        data,
                    )?;
                }
            }
            let times: Vec<f32> = self.sample_times().iter().map(|&t| t as f32).collect();
            npy::write(&dir.join("time.npy"), &[times.len()], &times)?;
            let mut csv = String::from("receiver,x,y,z\n");
            for (n, &node) in r.nodes.iter().enumerate() {
                let [x, y, z] = self.solver.grid.position(node);
                writeln!(csv, "{n},{x},{y},{z}")?;
            }
            fs::write(dir.join("receivers.csv"), csv)?;
        }
        fs::write(
            dir.join("run.json"),
            serde_json::to_string_pretty(&self.info)? + "\n",
        )?;
        Ok(())
    }
}

impl RunInfo {
    /// A few lines describing the run, for people.
    pub fn summary(&self) -> String {
        let [nx, ny, nz] = self.points;
        let [dx, dy, dz] = self.spacing;
        let [ex, ey, ez] = self.extent;
        let mut s = String::new();
        let _ = writeln!(
            s,
            "grid      {nx} × {ny} × {nz} nodes, spacing {dx} × {dy} × {dz} m, extent {ex} × {ey} × {ez} m"
        );
        let _ = writeln!(
            s,
            "model     vp {:.0}–{:.0} m/s, vs {:.0}–{:.0} m/s, rho {:.0}–{:.0} kg/m³",
            self.vp_range[0],
            self.vp_range[1],
            self.vs_range[0],
            self.vs_range[1],
            self.rho_range[0],
            self.rho_range[1]
        );
        let _ = writeln!(
            s,
            "time      dt {:.4} ms ({:.0}% of the {:.4} ms limit), {} steps to {:.3} s",
            self.dt * 1e3,
            100.0 * self.dt / self.dt_stable,
            self.dt_stable * 1e3,
            self.steps,
            self.duration
        );
        let _ = writeln!(
            s,
            "accuracy  {:.1} points per shortest wavelength ({:.0} m/s at {:.1} Hz), order-{} stencils{}{}",
            self.points_per_wavelength,
            self.slowest_wave,
            self.max_frequency,
            self.order,
            match self.band {
                Some(b) => format!(" optimised up to k·h = {b:.2}"),
                None => String::new(),
            },
            if self.remove_time_dispersion {
                ", time dispersion removed from the traces"
            } else {
                ""
            }
        );
        let top = if self.free_surface {
            "free surface on top"
        } else {
            "absorbing on every side"
        };
        let _ = writeln!(
            s,
            "boundary  absorbing layers {} nodes thick, {top}",
            self.pml_width
        );
        let _ = writeln!(
            s,
            "output    {} source(s), {} receiver(s) × {} samples of {}, {} snapshot series",
            self.sources.len(),
            self.receivers,
            self.samples,
            self.components.join(", "),
            self.snapshots
        );
        let _ = write!(
            s,
            "memory    {} at most; `vibe-rations memory` itemises it",
            human(self.memory.peak)
        );
        s
    }
}

/// A duration for people: 12s, 3m05s, 1h02m.
pub fn human_time(seconds: f64) -> String {
    let s = seconds.round() as u64;
    if s >= 3600 {
        format!("{}h{:02}m", s / 3600, s % 3600 / 60)
    } else if s >= 60 {
        format!("{}m{:02}s", s / 60, s % 60)
    } else {
        format!("{s}s")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stability_limit_matches_the_uniform_grid_formula() {
        let g = Grid::new([10; 3], [10.0; 3]);
        let expect = 10.0 / (3000.0 * (9.0 / 8.0 + 1.0 / 24.0) * 3f64.sqrt());
        let fourth = Stencil::taylor(4).unwrap();
        assert!((stable_dt(&g, 3000.0, &fourth) - expect).abs() < 1e-15);
        // Longer stencils have larger coefficient sums, so shorter steps.
        let eighth = Stencil::taylor(8).unwrap();
        assert!(stable_dt(&g, 3000.0, &eighth) < stable_dt(&g, 3000.0, &fourth));
    }

    #[test]
    fn formats_durations() {
        assert_eq!(human_time(42.4), "42s");
        assert_eq!(human_time(125.0), "2m05s");
        assert_eq!(human_time(3725.0), "1h02m");
    }

    use crate::config::Config;
    use crate::solver::Band;

    const BASE: &str = r#"
        [grid]
        points = [30, 30, 30]
        spacing = 10.0
        [time]
        duration = 0.05
        [boundary]
        pml_width = 8
        [[model.layer]]
        top = 0.0
        vp = 2000.0
        vs = 1000.0
        rho = 2000.0
        [[source]]
        position = [150.0, 150.0, 100.0]
        explosion = 1e12
        wavelet = { kind = "ricker", frequency = 5.0 }
        [receivers]
        [[receivers.point]]
        position = [150.0, 150.0, 0.0]
    "#;

    fn setup(text: &str) -> Result<Simulation> {
        Simulation::new(&Config::from_toml(text)?, Path::new("."))
    }

    fn error(text: &str) -> String {
        format!("{:#}", setup(text).err().expect("an error"))
    }

    /// The largest difference between `a` and `b` at the nodes `interior` holds, in each
    /// array relative to the largest value it has in `scale`.
    fn interior_error(
        g: &Grid,
        interior: &[std::ops::Range<usize>; 3],
        a: &crate::solver::Fields,
        b: &crate::solver::Fields,
        scale: &[f32; 9],
    ) -> f64 {
        use crate::solver::Field::*;
        let mut worst = 0.0f64;
        for (c, f) in [Vx, Vy, Vz, Sxx, Syy, Szz, Sxy, Sxz, Syz]
            .into_iter()
            .enumerate()
        {
            let (a, b) = (a.get(f), b.get(f));
            for k in interior[2].clone() {
                for j in interior[1].clone() {
                    for i in interior[0].clone() {
                        let p = g.idx(i, j, k);
                        worst = worst.max(((a[p] - b[p]).abs() / scale[c]) as f64);
                    }
                }
            }
        }
        worst
    }

    /// Stepping back from the final state, with the band saved on the way forward, retraces
    /// the run outside the layers: sources, free surface and all, to round-off.
    #[test]
    fn steps_back_inside_the_band() {
        let text = BASE
            .replace("pml_width = 8", "pml_width = 5")
            .replace("duration = 0.05", "duration = 0.15")
            .replace(
                "[receivers]",
                "[[source]]\nposition = [120.0, 160.0, 50.0]\nforce = [1e9, 0.0, 5e8]\n\
                 wavelet = { kind = \"ricker\", frequency = 8.0 }\n[receivers]",
            );
        for (label, text) in [
            ("free surface", text.clone()),
            (
                "layers all round",
                text.replace("pml_width = 5", "pml_width = 5\nfree_surface = false"),
            ),
        ] {
            let mut sim = setup(&text).unwrap();
            let band = sim.band();
            let interior = sim.solver.interior();
            let free = sim.solver.free_surface;
            let size = Band::size(&interior, free, sim.solver.stencil.half_length());
            assert_eq!(
                (band.len(), band.memory_bytes()),
                (size.values, size.runs_bytes())
            );
            let len = band.len();
            let steps = sim.steps;
            let mut records = vec![0.0; steps * len];
            let mut states = Vec::with_capacity(steps);
            for n in 0..steps {
                band.save(&sim.solver.fields, &mut records[n * len..(n + 1) * len]);
                states.push(sim.solver.fields.clone());
                sim.step();
            }
            let mut scale = [0.0f32; 9];
            for f in &states {
                use crate::solver::Field::*;
                for (c, a) in [Vx, Vy, Vz, Sxx, Syy, Szz, Sxy, Sxz, Syz]
                    .into_iter()
                    .enumerate()
                {
                    scale[c] = f.get(a).iter().fold(scale[c], |m, v| m.max(v.abs()));
                }
            }
            let mut worst = 0.0f64;
            for n in (0..steps).rev() {
                sim.step_back(&band, &records[n * len..(n + 1) * len]);
                assert_eq!(sim.solver.step, n);
                let e = interior_error(
                    &sim.solver.grid,
                    &interior,
                    &sim.solver.fields,
                    &states[n],
                    &scale,
                );
                worst = worst.max(e);
            }
            println!("{label}: {steps} steps back, at most {worst:.1e} of each array's peak");
            assert!(worst < 1e-4, "{label}: {worst}");
        }
    }

    #[test]
    fn steps_at_a_fraction_of_the_stability_limit() {
        // Eighth order by default, with the time error removed, at 0.85 of the limit.
        let sim = setup(BASE).unwrap();
        let i = &sim.info;
        assert_eq!(i.order, 8);
        assert!(i.remove_time_dispersion);
        assert!((i.dt - 0.85 * i.dt_stable).abs() < 1e-15);
        assert_eq!(sim.steps, (0.05 / i.dt).ceil() as usize);
        assert_eq!(i.samples, sim.steps);
        assert!(i.warnings.is_empty(), "{:?}", i.warnings);
        // Opting out, eighth order needs 0.3; fourth order keeps 0.85, where its time and
        // space errors partly cancel.
        let kept = setup(&format!("{BASE}\n[scheme]\nremove_time_dispersion = false")).unwrap();
        assert!((kept.info.dt - 0.3 * kept.info.dt_stable).abs() < 1e-15);
        let fourth = setup(&format!(
            "{BASE}\n[scheme]\norder = 4\nremove_time_dispersion = false"
        ))
        .unwrap();
        assert!((fourth.info.dt - 0.85 * fourth.info.dt_stable).abs() < 1e-15);
        // A grid fine enough for eighth order is too coarse for second.
        let ppw = i.points_per_wavelength;
        assert!(
            ppw >= min_points_per_wavelength(8, Coefficients::Taylor)
                && ppw < min_points_per_wavelength(2, Coefficients::Taylor)
        );
        let second = setup(&format!("{BASE}\n[scheme]\norder = 2")).unwrap();
        assert_eq!(second.info.warnings.len(), 1, "{:?}", second.info.warnings);
        let t = sim.sample_times();
        assert!((t[0] - 0.5 * i.dt).abs() < 1e-15);
        assert!((t[1] - t[0] - i.dt).abs() < 1e-15);
    }

    /// Optimised coefficients, the default, are designed for the run's shortest wavelength
    /// unless told otherwise, and their larger Σ|cₙ| shortens the time step.
    #[test]
    fn optimises_coefficients_for_the_grid() {
        let taylor = setup(&format!("{BASE}\n[scheme]\ncoefficients = \"taylor\"")).unwrap();
        assert_eq!(taylor.info.band, None);
        let sim = setup(BASE).unwrap();
        let band = sim.info.band.expect("a band");
        assert!(band <= 2.0 * PI / sim.info.points_per_wavelength + 1e-12);
        assert_eq!(sim.info.coefficients, Coefficients::Optimised);
        assert!(sim.info.dt_stable < taylor.info.dt_stable);
        assert_ne!(sim.solver.stencil, taylor.solver.stencil);
        assert!(sim.info.summary().contains("optimised up to"));
        let fixed = setup(&format!("{BASE}\n[scheme]\noptimised_for = 6.0")).unwrap();
        assert!((fixed.info.band.unwrap() - 2.0 * PI / 6.0).abs() < 1e-12);
        assert!(setup(&format!("{BASE}\n[scheme]\noptimised_for = 1.5")).is_err());
    }

    #[test]
    fn takes_a_stable_dt_and_rejects_an_unstable_one() {
        let sim = setup(&BASE.replace("duration = 0.05", "duration = 0.05\ndt = 1e-3")).unwrap();
        assert_eq!(sim.info.dt, 1e-3);
        let e = error(&BASE.replace("duration = 0.05", "duration = 0.05\ndt = 0.01"));
        assert!(e.contains("unstable"), "{e}");
        let e = error(&BASE.replace("duration = 0.05", "duration = 0.05\ncourant = 1.5"));
        assert!(e.contains("courant"), "{e}");
    }

    #[test]
    fn rejects_what_does_not_fit() {
        let e = error(&BASE.replace("[150.0, 150.0, 100.0]", "[150.0, 150.0, 400.0]"));
        assert!(e.contains("source") && e.contains("off the grid"), "{e}");
        let e = error(&BASE.replace("[150.0, 150.0, 0.0]", "[150.0, -20.0, 0.0]"));
        assert!(e.contains("receiver") && e.contains("off the grid"), "{e}");
        let e = error(&BASE.replace("pml_width = 8", "pml_width = 13"));
        assert!(e.contains("no room"), "{e}");
        let start = BASE.find("[[source]]").unwrap();
        let end = BASE.find("[receivers]").unwrap();
        let e = error(&format!("{}{}", &BASE[..start], &BASE[end..]));
        assert!(e.contains("no [[source]]"), "{e}");
        let e = error(&format!("{BASE}\n[output]\nrecord_every = 0"));
        assert!(e.contains("record_every"), "{e}");
        let e = error(&BASE.replace(
            "[receivers]",
            "[receivers]\ncomponents = [\"vz\", \"p\", \"p\"]",
        ));
        assert!(e.contains("list p twice"), "{e}");
    }

    #[test]
    fn runs_a_model_given_in_memory_on_the_configured_grid() {
        let cfg = Config::from_toml(BASE).unwrap();
        let grid = Simulation::grid(&cfg).unwrap();
        assert_eq!(grid.halo, 4);
        let model = Model::homogeneous(grid, 2500.0, 1200.0, 2100.0);
        let sim = Simulation::with_model(&cfg, model).unwrap();
        assert_eq!(sim.info.vp_range, [2500.0, 2500.0]);
        let spacing = Model::homogeneous(
            Grid::new([30, 30, 30], [10.0, 10.0, 5.0]).with_halo(4),
            2500.0,
            1200.0,
            2100.0,
        );
        let e = format!("{:#}", Simulation::with_model(&cfg, spacing).err().unwrap());
        assert!(e.contains("[10.0, 10.0, 5.0]"), "{e}");
        let wrong = Model::homogeneous(
            Grid::new([30, 30, 31], [10.0; 3]).with_halo(4),
            2500.0,
            1200.0,
            2100.0,
        );
        let e = format!("{:#}", Simulation::with_model(&cfg, wrong).err().unwrap());
        assert!(e.contains("grid"), "{e}");
        let mut bad = Model::homogeneous(grid, 2500.0, 1200.0, 2100.0);
        bad.vp[0] = -1.0;
        assert!(Simulation::with_model(&cfg, bad).is_err());
    }

    #[test]
    fn warns_about_coarse_grids_thin_layers_and_receivers_in_layers() {
        let text = BASE
            .replace("frequency = 5.0", "frequency = 20.0")
            .replace("pml_width = 8", "pml_width = 5")
            .replace("[150.0, 150.0, 0.0]", "[20.0, 150.0, 0.0]");
        let w = setup(&text).unwrap().info.warnings;
        assert_eq!(w.len(), 3, "{w:?}");
        assert!(w[0].contains("points per shortest wavelength"), "{w:?}");
        assert!(w[1].contains("5 nodes thick"), "{w:?}");
        assert!(w[2].contains("receiver(s) inside"), "{w:?}");
    }

    /// A velocity is interpolated from as many staggered values as the stencil spans: at the
    /// default order 8 four either side, so with 8-node layers (nodes 0–7, the last one's
    /// half-node value undamped) a horizontal velocity reads a damped value up to node 10 and
    /// not from node 11; pressure, and vertical velocity at the free surface, never do.
    #[test]
    fn finds_velocity_receivers_next_to_the_layers() {
        let next = |at: &str, components: &str| {
            let text = BASE.replace("[150.0, 150.0, 0.0]", at).replace(
                "[receivers]",
                &format!("[receivers]\ncomponents = {components}"),
            );
            let sim = setup(&text).unwrap();
            assert!(sim.info.warnings.is_empty(), "{:?}", sim.info.warnings);
            sim.receivers_next_to_layers()
        };
        assert_eq!(next("[80.0, 150.0, 0.0]", r#"["vx"]"#), 1);
        assert_eq!(next("[100.0, 150.0, 0.0]", r#"["vx"]"#), 1);
        assert_eq!(next("[110.0, 150.0, 0.0]", r#"["vx"]"#), 0);
        assert_eq!(next("[80.0, 150.0, 0.0]", r#"["vz", "p"]"#), 0);
        // The bottom layer starts at node 22, its first half-node value damped: vz reads it
        // from node 19 (indices 15–22).
        assert_eq!(next("[150.0, 150.0, 190.0]", r#"["vz"]"#), 1);
        assert_eq!(next("[150.0, 150.0, 180.0]", r#"["vz"]"#), 0);
    }

    /// The layers' damping scales with the speed it is set for, the model's fastest unless
    /// pml_vp_max says otherwise.
    #[test]
    fn pml_vp_max_sets_the_damping() {
        let damping = |text: &str| {
            let sim = setup(text).unwrap();
            sim.solver.cpml[0].profile(crate::cpml::Stagger::Node).d[0]
        };
        let model = damping(BASE);
        let twice = damping(&BASE.replace("pml_width = 8", "pml_width = 8\npml_vp_max = 4000.0"));
        assert!((twice / model - 2.0).abs() < 1e-6, "{twice} vs {model}");
        let e = error(&BASE.replace("pml_width = 8", "pml_width = 8\npml_vp_max = 0.0"));
        assert!(e.contains("pml_vp_max"), "{e}");
    }

    #[test]
    fn records_every_nth_step_and_snaps_snapshot_planes() {
        let text = format!(
            "{BASE}\n[[snapshot]]\nfield = \"vz\"\nplane = \"xz\"\nat = 154.0\ninterval = 0.01\n\
             [output]\nrecord_every = 3"
        );
        let sim = setup(&text).unwrap();
        assert_eq!(sim.receivers.samples, sim.steps.div_ceil(3));
        let t = sim.sample_times();
        assert!((t[1] - t[0] - 3.0 * sim.info.dt).abs() < 1e-15);
        assert_eq!(sim.snapshots[0].index, 15);
        assert_eq!(
            sim.snapshots[0].every,
            (0.01 / sim.info.dt).round() as usize
        );
        let e = error(&text.replace("at = 154.0", "at = 900.0"));
        assert!(e.contains("snapshot plane"), "{e}");
    }

    #[test]
    fn writes_traces_with_their_shapes() {
        let mut sim = setup(BASE).unwrap();
        sim.run(None, false).unwrap();
        assert_eq!(sim.receivers.recorded(), sim.receivers.samples);
        let dir = std::env::temp_dir().join(format!("vibe-rations-sim-{}", std::process::id()));
        sim.write_outputs(&dir).unwrap();
        let (shape, _) = npy::read(&dir.join("vz.npy")).unwrap();
        assert_eq!(shape, vec![1, sim.receivers.samples]);
        let (shape, _) = npy::read(&dir.join("time.npy")).unwrap();
        assert_eq!(shape, vec![sim.receivers.samples]);
        // One receiver makes no gather.
        assert!(!dir.join("vz_gather.png").exists());
        assert!(sim.info.summary().contains("free surface on top"));
        fs::remove_dir_all(dir).unwrap();
    }
}
