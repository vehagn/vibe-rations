//! Full waveform inversion: L-BFGS (Nocedal 1980; Liu & Nocedal 1989) on the gradients of
//! the `gradient` module, as Fichtner (2011, ch. 7) describes it for seismic models.
//!
//! The unknowns are relative changes, x = m/m₀ - 1 for each inverted parameter at each free
//! node, so that vp, vs and ρ share a scale; the gradient with respect to x is m₀ ∂χ/∂m. Free
//! nodes are those outside the absorbing layers, where the gradient is exact, and below
//! `fixed_above`; vs stays zero where it starts at zero, so fluids stay fluid.
//!
//! The time step, the layers' damping and the optimised coefficients' band are fixed at the
//! starting model's (`[time] dt`, `[boundary] pml_vp_max`, `[scheme] optimised_for`), so the
//! misfit is a smooth function of the model and the gradient is exactly its derivative. vp is
//! kept below what that time step allows with that operator.
//!
//! Each iteration: the L-BFGS direction (the two-loop recursion), scaled so that no value
//! changes by more than `initial_step` when there are no model and gradient changes to learn
//! the scale from, and by more than `max_step` ever; trial models projected onto the bounds;
//! Armijo backtracking with misfit-only forward runs, shortening the step by quadratic
//! interpolation; then the gradient at the accepted model, whose change is remembered if it
//! has positive curvature (sᵀy > 0).

use std::collections::VecDeque;
use std::path::Path;
use std::sync::Arc;

use anyhow::{Context, Result, bail, ensure};
use serde::Serialize;

use crate::config::{Coefficients, Config, InversionConfig, Parameter, StageConfig, shot_dir};
use crate::gradient::gradient;
use crate::misfit::{Measure, Misfit};
use crate::model::Model;
use crate::npy;
use crate::simulation::{Simulation, stable_dt};

/// The limited-memory BFGS approximation of the inverse Hessian: the last few changes of the
/// unknowns (s) and of the gradient (y).
#[derive(Clone, Debug, Default)]
pub struct Lbfgs {
    memory: usize,
    pairs: VecDeque<(Vec<f64>, Vec<f64>)>,
}

impl Lbfgs {
    pub fn new(memory: usize) -> Self {
        Self {
            memory,
            pairs: VecDeque::with_capacity(memory),
        }
    }

    /// Pairs remembered.
    pub fn len(&self) -> usize {
        self.pairs.len()
    }

    pub fn is_empty(&self) -> bool {
        self.pairs.is_empty()
    }

    pub fn clear(&mut self) {
        self.pairs.clear();
    }

    /// Bytes its pairs hold.
    pub fn memory_bytes(&self) -> usize {
        8 * self
            .pairs
            .iter()
            .map(|(s, y)| s.len() + y.len())
            .sum::<usize>()
    }

    /// -H g by the two-loop recursion (Nocedal & Wright 2006, algorithm 7.4), with the initial
    /// inverse Hessian γP: P the diagonal `preconditioner` (or the identity), γ = sᵀy / yᵀPy of
    /// the newest pair, or 1 if there are none.
    pub fn direction(&self, g: &[f64], preconditioner: Option<&[f64]>) -> Vec<f64> {
        let mut q = g.to_vec();
        let mut alphas = Vec::with_capacity(self.pairs.len());
        for (s, y) in self.pairs.iter().rev() {
            let rho = 1.0 / dot(y, s);
            let a = rho * dot(s, &q);
            axpy(&mut q, -a, y);
            alphas.push(a);
        }
        let apply = |v: &mut [f64]| {
            if let Some(p) = preconditioner {
                v.iter_mut().zip(p).for_each(|(v, p)| *v *= p);
            }
        };
        if let Some((s, y)) = self.pairs.back() {
            let mut py = y.clone();
            apply(&mut py);
            let gamma = dot(s, y) / dot(y, &py);
            q.iter_mut().for_each(|v| *v *= gamma);
        }
        apply(&mut q);
        for ((s, y), a) in self.pairs.iter().zip(alphas.iter().rev()) {
            let rho = 1.0 / dot(y, s);
            let b = rho * dot(y, &q);
            axpy(&mut q, a - b, s);
        }
        q.iter_mut().for_each(|v| *v = -*v);
        q
    }

    /// Remembers the change `s` of the unknowns and `y` of the gradient, forgetting the oldest
    /// pair beyond the memory, unless the curvature sᵀy isn't clearly positive (then the BFGS
    /// update would not stay positive definite); says whether it did.
    pub fn update(&mut self, s: Vec<f64>, y: Vec<f64>) -> bool {
        let sy = dot(&s, &y);
        // Written so that a NaN curvature counts as not positive.
        let positive = sy > 1e-12 * norm(&s) * norm(&y);
        if !positive || self.memory == 0 {
            return false;
        }
        if self.pairs.len() == self.memory {
            self.pairs.pop_front();
        }
        self.pairs.push_back((s, y));
        true
    }
}

/// How far a step may go.
#[derive(Clone, Copy, Debug)]
pub struct StepLimits {
    /// The largest change of any unknown when L-BFGS has nothing to scale its direction by.
    pub initial: f64,
    /// The largest change of any unknown, ever.
    pub max: f64,
}

/// A step L-BFGS took.
#[derive(Clone, Debug)]
pub struct Step {
    pub x: Vec<f64>,
    pub value: f64,
    /// The fraction of the (scaled) direction taken.
    pub length: f64,
    /// Values computed to find it.
    pub trials: usize,
    /// Whether the memory was dropped because its direction didn't go downhill.
    pub restarted: bool,
}

/// Looks for a step from `x` (value `f`, gradient `g`) along the L-BFGS direction (with the
/// diagonal `preconditioner`), scaled to `limits` and projected by `project`, shortened until
/// the value (from `value`) falls by at least 10⁻⁴ of what the gradient predicts for the
/// projected step (Armijo).
///
/// Unknowns held at a bound that the gradient pushes against are left out of the direction,
/// so that the search runs along the bound's face instead of into it (projected gradient). A
/// trial `value` can't compute, or whose value isn't finite, counts as too long a step. If the
/// L-BFGS direction finds no step in eight tries, the memory is dropped and steepest descent
/// tried; None if that fails too. An error only if every trial was one.
#[allow(clippy::too_many_arguments)]
pub fn line_search(
    lbfgs: &mut Lbfgs,
    x: &[f64],
    f: f64,
    g: &[f64],
    preconditioner: Option<&[f64]>,
    limits: StepLimits,
    project: impl Fn(&mut [f64]),
    mut value: impl FnMut(&[f64]) -> Result<f64>,
) -> Result<Option<Step>> {
    // Unknowns a small step down the gradient leaves where they are: pinned to a bound.
    let largest_g = g.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    if largest_g == 0.0 {
        return Ok(None);
    }
    let mut probe: Vec<f64> = x
        .iter()
        .zip(g)
        .map(|(x, g)| x - 1e-6 * g / largest_g)
        .collect();
    project(&mut probe);
    // (A zero gradient also leaves an unknown in place, but no bound holds it there.)
    let pinned: Vec<bool> = probe
        .iter()
        .zip(x)
        .zip(g)
        .map(|((p, x), g)| p == x && *g != 0.0)
        .collect();
    let reduced: Vec<f64> = g
        .iter()
        .zip(&pinned)
        .map(|(&g, &p)| if p { 0.0 } else { g })
        .collect();

    let mut restarted = false;
    let mut last_error = None;
    let mut errors = 0;
    let mut trials = 0;
    loop {
        let mut d = lbfgs.direction(&reduced, preconditioner);
        if dot(&reduced, &d) >= 0.0 && !lbfgs.is_empty() {
            lbfgs.clear();
            restarted = true;
            d = lbfgs.direction(&reduced, preconditioner);
        }
        d.iter_mut()
            .zip(&pinned)
            .filter(|(_, p)| **p)
            .for_each(|(d, _)| *d = 0.0);
        let largest = d.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        if largest == 0.0 || dot(&reduced, &d) >= 0.0 {
            return Ok(None);
        }
        let cap = if lbfgs.is_empty() {
            limits.initial
        } else {
            limits.max
        };
        if largest > cap || lbfgs.is_empty() {
            let scale = cap / largest;
            d.iter_mut().for_each(|v| *v *= scale);
        }
        const ARMIJO: f64 = 1e-4;
        let mut alpha = 1.0;
        for _ in 0..8 {
            trials += 1;
            let mut t: Vec<f64> = x.iter().zip(&d).map(|(a, b)| a + alpha * b).collect();
            project(&mut t);
            let s: Vec<f64> = t.iter().zip(x).map(|(a, b)| a - b).collect();
            // What the gradient predicts for the step actually taken, after projection.
            let predicted = dot(g, &s);
            let ft = match value(&t) {
                Ok(v) if v.is_finite() => v,
                Ok(_) => {
                    alpha *= 0.5;
                    continue;
                }
                Err(e) => {
                    errors += 1;
                    last_error = Some(e);
                    alpha *= 0.5;
                    continue;
                }
            };
            if predicted < 0.0 && ft <= f + ARMIJO * predicted {
                return Ok(Some(Step {
                    x: t,
                    value: ft,
                    length: alpha,
                    trials,
                    restarted,
                }));
            }
            // The minimum of the parabola through f, its slope and ft, kept within [0.1, 0.5]
            // of the step that failed.
            let slope = predicted / alpha;
            let curvature = ft - f - slope * alpha;
            let next = if curvature > 0.0 && slope < 0.0 {
                -slope * alpha * alpha / (2.0 * curvature)
            } else {
                0.5 * alpha
            };
            alpha = next.clamp(0.1 * alpha, 0.5 * alpha);
        }
        if lbfgs.is_empty() {
            break;
        }
        lbfgs.clear();
        restarted = true;
    }
    match last_error {
        Some(e) if errors == trials => Err(e),
        _ => Ok(None),
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

fn axpy(y: &mut [f64], a: f64, x: &[f64]) {
    y.iter_mut().zip(x).for_each(|(y, x)| *y += a * x);
}

/// One shot of the survey: its configuration (one source, the inversion's time step and
/// layers) and observed traces, with the misfit the current stage measures them by.
#[derive(Clone, Debug)]
pub struct Shot {
    pub config: Config,
    pub observed: Arc<Vec<Vec<f32>>>,
    pub misfit: Misfit,
    /// The receivers' components, samples per trace and sample interval.
    components: Vec<crate::receiver::Component>,
    samples: usize,
    interval: f64,
}

/// The misfit, its gradient with respect to the unknowns, and the forward wavefields' energy
/// at every interior node, all summed over the shots.
pub struct Evaluation {
    pub misfit: f64,
    pub gradient: Vec<f64>,
    pub illumination: Vec<f64>,
    /// The most memory held at once while it was computed: the inversion's own arrays, the
    /// sums over the shots and one shot's gradient run ([`crate::gradient::Run::memory_bytes`]).
    pub memory_bytes: usize,
}

/// What an iteration did.
#[derive(Clone, Debug, Serialize)]
pub struct Iteration {
    pub iteration: usize,
    /// The stage it belongs to, counted from 0.
    pub stage: usize,
    pub misfit: f64,
    /// The fraction of the scaled L-BFGS direction taken.
    pub step: f64,
    /// Forward runs spent finding the step.
    pub trials: usize,
    /// The largest relative change of any value.
    pub largest_change: f64,
    /// Pairs L-BFGS remembers after it.
    pub memory: usize,
    pub restarted: bool,
}

/// An inversion in progress.
pub struct Inversion {
    pub shots: Vec<Shot>,
    /// The model as it now stands.
    pub model: Model,
    start: Model,
    /// The inverted parameters and, for each, the nodes (interior indices) it may change at.
    free: Vec<(Parameter, Vec<usize>)>,
    bounds: Vec<[f64; 2]>,
    limits: StepLimits,
    /// The stages, each with its misfit, and the one under way.
    stages: Vec<StageConfig>,
    stage: usize,
    /// The sources' dominant frequency, Hz: the phase misfit's scale on unfiltered traces.
    source_frequency: f64,
    /// The misfit at the start of each stage begun, as that stage measures it.
    pub stage_misfits: Vec<f64>,
    settings: InversionConfig,
    /// The water level of the illumination, if the inversion divides by it.
    water_level: Option<f64>,
    /// Each free parameter's factor on the preconditioner (`[inversion] scale`).
    scale: Vec<f64>,
    /// Runs that estimate the receivers' side of the illumination at each stage's start, and
    /// that estimate (Σₜ |v|² at every interior node), if there are any.
    receiver_runs: usize,
    receiver_side: Option<Vec<f64>>,
    preconditioner: Option<Vec<f64>>,
    lbfgs: Lbfgs,
    x: Vec<f64>,
    misfit: f64,
    g: Vec<f64>,
    /// Iterations done.
    pub iteration: usize,
    pub history: Vec<Iteration>,
    /// The most memory held at once so far, counted from the arrays
    /// ([`crate::memory::Plan::inversion_memory`] works it out beforehand).
    pub memory_bytes: usize,
}

impl Inversion {
    /// Sets up an inversion from a configuration with an `[inversion]` section, its own
    /// model to start from and the observed traces under `base`, as `invert` does.
    pub fn from_config(cfg: &Config, base: &Path) -> Result<Self> {
        let settings = cfg
            .inversion
            .clone()
            .context("the configuration has no [inversion] section")?;
        let model = Model::from_config(Simulation::grid(cfg)?, &cfg.model, base)?;
        let dir = base.join(&settings.observed);
        let observed = (0..cfg.sources.len())
            .map(|s| {
                // run --shots writes shot-NNN/ even for one source; a plain run doesn't.
                let shot = dir.join(shot_dir(s));
                let shot = if shot.is_dir() || cfg.sources.len() > 1 {
                    shot
                } else {
                    dir.clone()
                };
                Observed::load(cfg, &shot)
            })
            .collect::<Result<Vec<_>>>()?;
        Self::new(cfg, model, observed, &settings)
    }

    /// Sets up an inversion of `start` against `observed` traces, one set per source of `cfg`
    /// (each source a shot).
    pub fn new(
        cfg: &Config,
        start: Model,
        observed: Vec<Observed>,
        settings: &InversionConfig,
    ) -> Result<Self> {
        ensure!(
            observed.len() == cfg.sources.len(),
            "{} sets of observed traces for {} sources",
            observed.len(),
            cfg.sources.len()
        );
        ensure!(!settings.parameters.is_empty(), "nothing to invert for");
        for (i, p) in settings.parameters.iter().enumerate() {
            ensure!(
                !settings.parameters[..i].contains(p),
                "{} is listed twice",
                p.name()
            );
        }
        ensure!(
            settings.initial_step > 0.0 && settings.max_step >= settings.initial_step,
            "the steps must be positive, with max_step at least initial_step"
        );

        // Hold the operator, the time step and the layers at the starting model's. Optimised
        // coefficients follow the model's slowest wave unless pinned; pin them where the start
        // puts them (the band is capped at 0.9π, which 2/0.9 points per wavelength gives).
        let mut fixed = cfg.clone();
        if cfg.scheme.coefficients == Coefficients::Optimised && cfg.scheme.optimised_for.is_none()
        {
            let unpinned = Simulation::with_model(cfg, start.clone())?;
            fixed.scheme.optimised_for = Some(unpinned.info.points_per_wavelength.max(2.0 / 0.9));
        }
        let first = Simulation::with_model(&fixed, start.clone())?;
        let dt = first.solver.dt;
        fixed.time.dt = Some(dt);
        fixed.time.courant = None;
        fixed.boundary.pml_vp_max = Some(
            cfg.boundary
                .pml_vp_max
                .unwrap_or_else(|| start.vp_range()[1]),
        );
        let stages = settings.stages()?;
        let source_frequency = cfg
            .sources
            .iter()
            .map(|s| Ok(s.wavelet.wavelet()?.dominant_frequency()))
            .collect::<Result<Vec<f64>>>()?
            .into_iter()
            .fold(0.0, f64::max);
        let measure = Measure::of(&stages[0], settings, source_frequency);
        let mut shots = Vec::with_capacity(cfg.sources.len());
        for (s, observed) in observed.into_iter().enumerate() {
            let config = fixed.shot(s)?;
            let sim = Simulation::with_model(&config, start.clone())?;
            let r = &sim.receivers;
            observed.check(&sim).with_context(|| format!("shot {s}"))?;
            let observed = Arc::new(observed.traces);
            let (components, samples) = (r.components.clone(), r.samples);
            let interval = sim.info.sample_interval;
            let misfit = Misfit::new(measure, &observed, &components, samples, interval)
                .with_context(|| format!("shot {s}"))?;
            shots.push(Shot {
                config,
                observed,
                misfit,
                components,
                samples,
                interval,
            });
        }
        // Every stage's filter must suit the traces before the first one starts.
        for stage in &stages[1..] {
            let shot = &shots[0];
            Misfit::new(
                Measure::of(stage, settings, source_frequency),
                &shot.observed,
                &shot.components,
                shot.samples,
                shot.interval,
            )?;
        }

        // The nodes each parameter may change at.
        let g = start.grid;
        let [nx, ny, nz] = g.n;
        let layers = first.solver.cpml.each_ref().map(|c| c.clone());
        let fixed_depth = settings.fixed_above.unwrap_or(f64::NEG_INFINITY);
        let mut free = Vec::new();
        for &p in &settings.parameters {
            let mut nodes = Vec::new();
            for k in 0..nz {
                if (k as f64) * g.h[2] < fixed_depth {
                    continue;
                }
                for j in 0..ny {
                    for i in 0..nx {
                        let in_layers = [i, j, k]
                            .iter()
                            .zip(&layers)
                            .any(|(&n, c)| c.slab(n).is_some());
                        let n = g.interior_idx(i, j, k);
                        if !in_layers && !(p == Parameter::Vs && start.vs[n] == 0.0) {
                            nodes.push(n);
                        }
                    }
                }
            }
            ensure!(!nodes.is_empty(), "no node is free to change {}", p.name());
            free.push((p, nodes));
        }

        // Bounds; vp also below what the time step allows with this operator (with a margin
        // for rounding).
        let vp_limit = 0.999 * stable_dt(&g, 1.0, &first.solver.stencil) / dt;
        // Its memory would only add to the first gradient run's.
        drop(first);
        let around = |range: [f64; 2]| [0.5 * range[0], 2.0 * range[1]];
        let smallest_vs = start
            .vs
            .iter()
            .filter(|&&v| v > 0.0)
            .fold(f64::INFINITY, |m, &v| m.min(v as f64));
        let mut bounds = Vec::with_capacity(free.len());
        for (p, nodes) in &free {
            let [lo, hi] = match p {
                Parameter::Vp => settings
                    .vp_range
                    .unwrap_or_else(|| around(start.vp_range())),
                Parameter::Vs => settings
                    .vs_range
                    .unwrap_or_else(|| [0.5 * smallest_vs, 2.0 * start.vs_range()[1]]),
                Parameter::Rho => settings
                    .rho_range
                    .unwrap_or_else(|| around(start.rho_range())),
            };
            ensure!(
                lo >= 0.0 && lo < hi,
                "the {} range [{lo}, {hi}] is empty",
                p.name()
            );
            let hi = if *p == Parameter::Vp {
                ensure!(
                    lo < vp_limit,
                    "vp_range starts at {lo} m/s, above the {vp_limit:.0} m/s the time step allows"
                );
                hi.min(vp_limit)
            } else {
                hi
            };
            let values = match p {
                Parameter::Vp => &start.vp,
                Parameter::Vs => &start.vs,
                Parameter::Rho => &start.rho,
            };
            // A little slack: the default vp limit sits 0.1% below the stability limit.
            if let Some(&n) = nodes.iter().find(|&&n| {
                let v = values[n] as f64;
                v < lo * (1.0 - 1e-6) || v > hi * (1.0 + 2e-3)
            }) {
                bail!(
                    "the starting model's {} at node {n} ({}) is outside [{lo}, {hi}]",
                    p.name(),
                    values[n]
                );
            }
            bounds.push([lo, hi]);
        }
        // λ ≥ 0 (vs ≤ vp/√2) is a bound too: the start must keep it where either is free.
        for (p, nodes) in &free {
            if matches!(p, Parameter::Vp | Parameter::Vs)
                && let Some(&n) = nodes.iter().find(|&&n| {
                    2.0 * (start.vs[n] as f64).powi(2) > (start.vp[n] as f64).powi(2) * (1.0 + 1e-6)
                })
            {
                bail!(
                    "the starting model has vs > vp/√2 (λ < 0) at node {n}, which the \
                     inversion's bounds exclude"
                );
            }
        }

        let n: usize = free.iter().map(|(_, nodes)| nodes.len()).sum();
        let mut inversion = Self {
            shots,
            model: start.clone(),
            start,
            free,
            bounds,
            limits: StepLimits {
                initial: settings.initial_step,
                max: settings.max_step,
            },
            stages,
            stage: 0,
            source_frequency,
            stage_misfits: Vec::new(),
            settings: settings.clone(),
            water_level: settings.precondition.then_some(settings.water_level),
            scale: settings
                .parameters
                .iter()
                .map(|&p| settings.scale.of(p))
                .collect(),
            receiver_runs: settings.receiver_runs,
            receiver_side: None,
            preconditioner: None,
            lbfgs: Lbfgs::new(settings.memory),
            x: vec![0.0; n],
            misfit: 0.0,
            g: vec![0.0; n],
            iteration: 0,
            history: Vec::new(),
            memory_bytes: 0,
        };
        ensure!(
            settings.water_level > 0.0,
            "the water level must be positive"
        );
        inversion.illuminate_from_receivers(0)?;
        let e = inversion.evaluate(&inversion.model)?;
        inversion.memory_bytes = inversion.memory_bytes.max(e.memory_bytes);
        inversion.misfit = e.misfit;
        inversion.g = e.gradient;
        inversion.preconditioner = inversion.preconditioner_for(&e.illumination);
        inversion.stage_misfits.push(e.misfit);
        Ok(inversion)
    }

    /// The stages, each with its misfit.
    pub fn stages(&self) -> &[StageConfig] {
        &self.stages
    }

    /// The stage under way, counted from 0.
    pub fn stage(&self) -> usize {
        self.stage
    }

    /// Starts stage `s` from the model as it stands: each shot's misfit as the stage measures
    /// it, and L-BFGS afresh, since the pairs it learnt describe another misfit.
    pub fn begin_stage(&mut self, s: usize) -> Result<()> {
        ensure!(
            s < self.stages.len(),
            "there is no stage {s}; the inversion has {}",
            self.stages.len()
        );
        let measure = Measure::of(&self.stages[s], &self.settings, self.source_frequency);
        for shot in &mut self.shots {
            shot.misfit = Misfit::new(
                measure,
                &shot.observed,
                &shot.components,
                shot.samples,
                shot.interval,
            )?;
        }
        self.stage = s;
        self.lbfgs.clear();
        self.illuminate_from_receivers(s as u64)?;
        let e = self.evaluate(&self.model)?;
        self.memory_bytes = self.memory_bytes.max(e.memory_bytes);
        ensure!(
            e.misfit.is_finite(),
            "the misfit isn't finite at the start of stage {s}"
        );
        self.misfit = e.misfit;
        self.g = e.gradient;
        self.preconditioner = self.preconditioner_for(&e.illumination);
        self.stage_misfits.push(e.misfit);
        Ok(())
    }

    /// The misfit of the model as it now stands.
    pub fn misfit(&self) -> f64 {
        self.misfit
    }

    /// Unknowns: the inverted parameters at their free nodes.
    pub fn unknowns(&self) -> usize {
        self.x.len()
    }

    /// The model the unknowns `x` describe, within the bounds.
    pub fn model_at(&self, x: &[f64]) -> Model {
        let mut m = self.model.clone();
        let mut offset = 0;
        for (p, nodes) in &self.free {
            let (now, start) = match p {
                Parameter::Vp => (&mut m.vp, &self.start.vp),
                Parameter::Vs => (&mut m.vs, &self.start.vs),
                Parameter::Rho => (&mut m.rho, &self.start.rho),
            };
            for (i, &n) in nodes.iter().enumerate() {
                now[n] = (start[n] as f64 * (1.0 + x[offset + i])) as f32;
            }
            offset += nodes.len();
        }
        m
    }

    /// Projects unknowns onto the bounds: each value into its range, and λ ≥ 0. vp goes first,
    /// whatever order the parameters are listed in, so vs is capped at the projected vp/√2;
    /// where vs isn't inverted, vp is kept at or above √2 vs instead.
    fn project(&self, x: &mut [f64]) {
        let mut offsets = Vec::with_capacity(self.free.len());
        let mut offset = 0;
        for (_, nodes) in &self.free {
            offsets.push(offset);
            offset += nodes.len();
        }
        let find = |p: Parameter| self.free.iter().position(|(q, _)| *q == p);
        let (vp, vs) = (find(Parameter::Vp), find(Parameter::Vs));
        let order = [Parameter::Vp, Parameter::Rho, Parameter::Vs]
            .into_iter()
            .filter_map(find);
        for f in order {
            let (p, nodes) = &self.free[f];
            let [lo, hi] = self.bounds[f];
            let start = match p {
                Parameter::Vp => &self.start.vp,
                Parameter::Vs => &self.start.vs,
                Parameter::Rho => &self.start.rho,
            };
            for (i, &n) in nodes.iter().enumerate() {
                let m0 = start[n] as f64;
                let mut v = (m0 * (1.0 + x[offsets[f] + i])).clamp(lo, hi);
                match p {
                    Parameter::Vs => {
                        let vp_n = match vp
                            .and_then(|v| self.free[v].1.binary_search(&n).ok().map(|i| (v, i)))
                        {
                            Some((v, i)) => self.start.vp[n] as f64 * (1.0 + x[offsets[v] + i]),
                            None => self.model.vp[n] as f64,
                        };
                        v = v.min(vp_n / 2f64.sqrt());
                    }
                    Parameter::Vp => {
                        let vs_free = vs.is_some_and(|v| self.free[v].1.binary_search(&n).is_ok());
                        if !vs_free {
                            let floor = 2f64.sqrt() * self.model.vs[n] as f64;
                            v = v.max(floor.min(m0).min(hi));
                        }
                    }
                    Parameter::Rho => {}
                }
                x[offsets[f] + i] = v / m0 - 1.0;
            }
        }
    }

    /// The misfit of `model`, summed over the shots: a forward run each.
    pub fn misfit_of(&self, model: &Model) -> Result<f64> {
        let mut total = 0.0;
        for shot in &self.shots {
            let mut sim = Simulation::with_model(&shot.config, model.clone())?;
            sim.run(None, false)?;
            total += shot.misfit.value(&sim.receivers.traces)?;
        }
        Ok(total)
    }

    /// The misfit of `model`, its gradient with respect to the unknowns and the forward
    /// wavefields' energy at every node, summed over the shots.
    pub fn evaluate(&self, model: &Model) -> Result<Evaluation> {
        self.evaluate_holding(model, 0)
    }

    /// [`Inversion::evaluate`], while the caller holds `extra` bytes besides the inversion.
    fn evaluate_holding(&self, model: &Model, extra: usize) -> Result<Evaluation> {
        let n = model.grid.interior_len();
        let mut total = 0.0;
        let mut sum = [vec![0.0; n], vec![0.0; n], vec![0.0; n]];
        let mut illumination = vec![0.0; n];
        let held = self.memory_held() + extra + 8 * 4 * n;
        let mut memory_bytes = held;
        for shot in &self.shots {
            let mut sim = Simulation::with_model(&shot.config, model.clone())?;
            let run = gradient(&mut sim, |r| shot.misfit.evaluate(&r.traces), None)?;
            total += run.misfit;
            memory_bytes = memory_bytes.max(held + run.memory_bytes);
            illumination
                .iter_mut()
                .zip(&run.illumination)
                .for_each(|(a, b)| *a += b);
            for (s, g) in sum.iter_mut().zip(run.gradient.velocities(model)) {
                s.iter_mut().zip(g).for_each(|(s, g)| *s += g);
            }
        }
        let [vp, vs, rho] = sum;
        let mut g = Vec::with_capacity(self.x.len());
        for (p, nodes) in &self.free {
            let (grad, start) = match p {
                Parameter::Vp => (&vp, &self.start.vp),
                Parameter::Vs => (&vs, &self.start.vs),
                Parameter::Rho => (&rho, &self.start.rho),
            };
            g.extend(nodes.iter().map(|&n| grad[n] * start[n] as f64));
        }
        Ok(Evaluation {
            misfit: total,
            gradient: g,
            illumination,
            memory_bytes,
        })
    }

    /// Bytes the inversion's own arrays hold: the observed traces, the model as it stands and
    /// the start, the free nodes, the unknowns, their gradient, the preconditioner and L-BFGS's
    /// pairs.
    fn memory_held(&self) -> usize {
        let observed: usize = self
            .shots
            .iter()
            .flat_map(|s| s.observed.iter())
            .map(Vec::len)
            .sum();
        let misfits: usize = self.shots.iter().map(|s| s.misfit.memory_bytes()).sum();
        let receiver_side = 8 * self.receiver_side.as_ref().map_or(0, Vec::len);
        let model = |m: &Model| 4 * (m.vp.len() + m.vs.len() + m.rho.len());
        let free: usize = self.free.iter().map(|(_, nodes)| nodes.len()).sum();
        4 * observed
            + misfits
            + receiver_side
            + model(&self.model)
            + model(&self.start)
            + std::mem::size_of::<usize>() * free
            + 8 * (self.x.len() + self.g.len() + self.preconditioner.as_ref().map_or(0, Vec::len))
            + self.lbfgs.memory_bytes()
    }

    /// For each unknown, its parameter's factor w divided by (I / I_max + water level), I the
    /// energy at its node, if the inversion divides by it; w alone if not, and none if every w
    /// is 1 too.
    ///
    /// L-BFGS with an initial matrix γ D² P is L-BFGS with γ P on the unknowns scaled by D⁻¹,
    /// D = diag(√w): the pairs (s, y) and the products sᵀy, the Armijo test and the step
    /// limits are unchanged, so a factor is a change of variables for its parameter, and the
    /// factors can differ between stages without forgetting the pairs.
    fn preconditioner_for(&self, illumination: &[f64]) -> Option<Vec<f64>> {
        let unit = self.scale.iter().all(|&w| w == 1.0);
        let scaled = |f: &dyn Fn(usize) -> f64| -> Vec<f64> {
            self.free
                .iter()
                .zip(&self.scale)
                .flat_map(|((_, nodes), &w)| nodes.iter().map(move |&n| w * f(n)))
                .collect()
        };
        // The source side times, if it was estimated, the receivers' side.
        let both = |n: usize| illumination[n] * self.receiver_side.as_ref().map_or(1.0, |r| r[n]);
        let max = self
            .free
            .iter()
            .flat_map(|(_, nodes)| nodes.iter().map(|&n| both(n)))
            .fold(0.0, f64::max);
        match self.water_level {
            Some(level) if max > 0.0 => Some(scaled(&|n| 1.0 / (both(n) / max + level))),
            _ if !unit => Some(scaled(&|_| 1.0)),
            _ => None,
        }
    }

    /// Estimates the receivers' side of the illumination in the model as it stands, if the
    /// inversion asks for it, with signs drawn from `seed`.
    fn illuminate_from_receivers(&mut self, seed: u64) -> Result<()> {
        if self.receiver_runs == 0 || self.water_level.is_none() {
            return Ok(());
        }
        let mut rng = seed;
        let mut sign = move |_: usize, _: usize, _: usize| {
            if splitmix64(&mut rng) >> 63 == 0 {
                1.0
            } else {
                -1.0
            }
        };
        let (energy, peak) = self.receiver_energy(&self.model, self.receiver_runs, &mut sign)?;
        self.memory_bytes = self.memory_bytes.max(peak);
        self.receiver_side = Some(energy);
        Ok(())
    }

    /// Σₜ |v|² at every interior node of waves sent out from every receiver of the first shot
    /// at once, averaged over `runs` runs, and the most memory held meanwhile. Each receiver
    /// component sends the first source's wavelet times √(its weight in the misfit, relative to
    /// the largest, averaged over the shots), times `sign(run, receiver, component)`,
    /// injected as the adjoint run injects residuals (`Solver::adjoint_op`), so velocity and
    /// pressure receivers count as they do in the gradient.
    ///
    /// With independent random signs its expectation is Σₐ wₐ Σₜ |Gₐ|², Gₐ the waves from
    /// receiver component a alone: the receivers' side of the Gauss–Newton Hessian's diagonal,
    /// which the source side's energy multiplies (Shin et al. 2001; Plessix & Mulder 2004).
    /// Firing them all at once with random signs, as encoded simultaneous sources (Krebs et
    /// al. 2009), costs a forward run per run instead of one per receiver component; the cross
    /// terms only average out, so a few runs give a rough estimate, enough for a
    /// preconditioner.
    fn receiver_energy(
        &self,
        model: &Model,
        runs: usize,
        sign: &mut dyn FnMut(usize, usize, usize) -> f64,
    ) -> Result<(Vec<f64>, usize)> {
        use crate::receiver::Component;
        let shot = &self.shots[0];
        let wavelet = shot.config.sources[0].wavelet.wavelet()?;
        let mut weights = vec![0.0; shot.components.len()];
        for s in &self.shots {
            let w = s.misfit.component_weights();
            let max = w.iter().fold(0.0f64, |m, &v| m.max(v));
            if max > 0.0 {
                weights.iter_mut().zip(&w).for_each(|(a, b)| *a += b / max);
            }
        }
        let n = model.grid.interior_len();
        let mut energy = vec![0.0f64; n];
        let held = self.memory_held() + 8 * n;
        let mut peak = held;
        for run in 0..runs.max(1) {
            let mut sim = Simulation::with_model(&shot.config, model.clone())?;
            peak = peak.max(held + sim.memory_bytes());
            let mut forces = Vec::new();
            let mut stresses = Vec::new();
            for (c, (&component, &w)) in shot.components.iter().zip(&weights).enumerate() {
                if w <= 0.0 {
                    continue;
                }
                for (r, &node) in sim.receivers.nodes.iter().enumerate() {
                    let op = sim.solver.adjoint_op(component, node);
                    let amplitude = sign(run, r, c) * (w / self.shots.len() as f64).sqrt();
                    match component {
                        Component::P => stresses.push((op, amplitude)),
                        _ => forces.push((op, amplitude)),
                    }
                }
            }
            for _ in 0..sim.steps {
                let send = |ops: &[(crate::solver::PointOp, f64)],
                            fields: &mut crate::solver::Fields,
                            t: f64| {
                    let w = wavelet.value(t);
                    for (op, a) in ops {
                        op.inject(fields, (a * w) as f32);
                    }
                };
                sim.solver
                    .step_with(|f, t| send(&forces, f, t), |f, t| send(&stresses, f, t));
                add_energy(&sim.solver, &mut energy);
            }
        }
        energy.iter_mut().for_each(|e| *e /= runs.max(1) as f64);
        Ok((energy, peak))
    }

    /// One iteration, or None if no step along the L-BFGS direction (or, after a restart,
    /// down the gradient) lowers the misfit.
    pub fn iterate(&mut self) -> Result<Option<Iteration>> {
        let mut lbfgs = std::mem::take(&mut self.lbfgs);
        let step = {
            let this = &*self;
            line_search(
                &mut lbfgs,
                &this.x,
                this.misfit,
                &this.g,
                this.preconditioner.as_deref(),
                this.limits,
                |t| this.project(t),
                |t| this.misfit_of(&this.model_at(t)),
            )
        };
        self.lbfgs = lbfgs;
        let Some(step) = step? else {
            return Ok(None);
        };
        let model = self.model_at(&step.x);
        // Besides the inversion: the step and the model it leads to.
        let extra = 8 * step.x.len() + 4 * (model.vp.len() + model.vs.len() + model.rho.len());
        let Evaluation {
            misfit,
            gradient: g_new,
            illumination,
            memory_bytes,
        } = self.evaluate_holding(&model, extra)?;
        self.memory_bytes = self.memory_bytes.max(memory_bytes);
        if !misfit.is_finite() {
            bail!(
                "the misfit isn't finite after iteration {}",
                self.iteration + 1
            );
        }
        let s: Vec<f64> = step.x.iter().zip(&self.x).map(|(a, b)| a - b).collect();
        let y: Vec<f64> = g_new.iter().zip(&self.g).map(|(a, b)| a - b).collect();
        let largest_change = step
            .x
            .iter()
            .zip(&self.x)
            .map(|(a, b)| ((1.0 + a) / (1.0 + b) - 1.0).abs())
            .fold(0.0, f64::max);
        self.lbfgs.update(s, y);
        self.preconditioner = self.preconditioner_for(&illumination);
        self.x = step.x;
        self.g = g_new;
        self.misfit = misfit;
        self.model = model;
        self.iteration += 1;
        let it = Iteration {
            iteration: self.iteration,
            stage: self.stage,
            misfit,
            step: step.length,
            trials: step.trials,
            largest_change,
            memory: self.lbfgs.len(),
            restarted: step.restarted,
        };
        self.history.push(it.clone());
        Ok(Some(it))
    }
}

/// Adds |v|² at every interior node, each component taken where it is stored.
fn add_energy(solver: &crate::solver::Solver, energy: &mut [f64]) {
    use rayon::prelude::*;
    let g = solver.grid;
    let [nx, ny, _] = g.n;
    let v = [&solver.fields.vx, &solver.fields.vy, &solver.fields.vz];
    energy
        .par_chunks_mut(nx * ny)
        .enumerate()
        .for_each(|(k, plane)| {
            for j in 0..ny {
                let row = g.idx(0, j, k);
                for i in 0..nx {
                    let p = row + i;
                    plane[j * nx + i] += v.iter().map(|c| (c[p] as f64).powi(2)).sum::<f64>();
                }
            }
        });
}

/// The next number of a SplitMix64 sequence (Steele, Lea & Flood 2014): reproducible signs
/// without a dependency.
fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// One shot's observed traces, in the layout of [`crate::receiver::Receivers::traces`], with
/// what is known of how they were recorded.
#[derive(Clone, Debug, Default)]
pub struct Observed {
    pub traces: Vec<Vec<f32>>,
    /// Each component's array shape, (receivers, samples), if it came from a file.
    pub shapes: Option<Vec<Vec<usize>>>,
    /// The sample times, if known (`time.npy`).
    pub times: Option<Vec<f32>>,
}

impl From<Vec<Vec<f32>>> for Observed {
    fn from(traces: Vec<Vec<f32>>) -> Self {
        Self {
            traces,
            ..Self::default()
        }
    }
}

impl Observed {
    /// The traces `run` wrote into `dir`: `<component>.npy` for every component the
    /// configuration records, and `time.npy` if it is there.
    pub fn load(cfg: &Config, dir: &Path) -> Result<Self> {
        let mut traces = Vec::new();
        let mut shapes = Vec::new();
        for c in &cfg.receivers.components {
            let path = dir.join(format!("{}.npy", c.name()));
            let (shape, data) =
                npy::read(&path).with_context(|| format!("reading {}", path.display()))?;
            traces.push(data);
            shapes.push(shape);
        }
        let time = dir.join("time.npy");
        let times = if time.is_file() {
            Some(npy::read(&time)?.1)
        } else {
            None
        };
        Ok(Self {
            traces,
            shapes: Some(shapes),
            times,
        })
    }

    /// Fails unless these traces fit `sim`'s receivers: every component, each (receivers,
    /// samples), sampled at the same times.
    pub fn check(&self, sim: &Simulation) -> Result<()> {
        let r = &sim.receivers;
        let want = [r.nodes.len(), r.samples];
        let hint = "observed and synthetic traces must share their sampling: give both runs \
                    the same [time] dt and record_every";
        ensure!(
            self.traces.len() == r.components.len(),
            "{} observed components for the {} recorded",
            self.traces.len(),
            r.components.len()
        );
        for (c, t) in r.components.iter().zip(&self.traces) {
            ensure!(
                t.len() == want[0] * want[1],
                "{} holds {} values for {} receivers × {} samples; {hint}",
                c.name(),
                t.len(),
                want[0],
                want[1]
            );
        }
        if let Some(shapes) = &self.shapes {
            for (c, shape) in r.components.iter().zip(shapes) {
                ensure!(
                    shape[..] == want[..],
                    "{}.npy has shape {shape:?}, not (receivers, samples) = {want:?}",
                    c.name()
                );
            }
        }
        if let Some(times) = &self.times {
            let synthetic = sim.sample_times();
            let off = times.len() != synthetic.len()
                || times.iter().zip(&synthetic).any(|(&a, &b)| {
                    (a as f64 - b).abs() > 1e-6 * sim.info.duration.max(1.0) + 1e-5 * b.abs()
                });
            ensure!(
                !off,
                "the observed traces are sampled at other times; {hint}"
            );
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Minimises with L-BFGS and the line search as Inversion does, from `x`, returning the
    /// end point and the iterations it took.
    fn minimise(
        mut x: Vec<f64>,
        f: impl Fn(&[f64]) -> (f64, Vec<f64>),
        project: impl Fn(&mut [f64]),
        limits: StepLimits,
        tolerance: f64,
    ) -> (Vec<f64>, usize) {
        let mut lbfgs = Lbfgs::new(5);
        let (mut value, mut g) = f(&x);
        for it in 0..500 {
            if norm(&g) < tolerance {
                return (x, it);
            }
            let Some(step) = line_search(&mut lbfgs, &x, value, &g, None, limits, &project, |t| {
                Ok(f(t).0)
            })
            .unwrap() else {
                return (x, it);
            };
            let (v, g_new) = f(&step.x);
            assert_eq!(v, step.value);
            let s = step.x.iter().zip(&x).map(|(a, b)| a - b).collect();
            let y = g_new.iter().zip(&g).map(|(a, b)| a - b).collect();
            lbfgs.update(s, y);
            (x, value, g) = (step.x, v, g_new);
        }
        (x, 500)
    }

    const FREE: StepLimits = StepLimits {
        initial: 0.5,
        max: 10.0,
    };

    /// On a quadratic, L-BFGS with enough memory finds the minimum in few iterations.
    #[test]
    fn minimises_a_quadratic() {
        // f = ½ xᵀ A x - bᵀ x with A = diag(1..5) plus a coupling, minimum A⁻¹ b.
        let a = |x: &[f64]| -> Vec<f64> {
            (0..5)
                .map(|i| {
                    (i + 1) as f64 * x[i]
                        + 0.3 * if i > 0 { x[i - 1] } else { 0.0 }
                        + 0.3 * if i < 4 { x[i + 1] } else { 0.0 }
                })
                .collect()
        };
        let b = [1.0, -2.0, 0.5, 3.0, 1.0];
        let f = |x: &[f64]| {
            let ax = a(x);
            (
                0.5 * dot(x, &ax) - dot(&b, x),
                ax.iter().zip(&b).map(|(p, q)| p - q).collect(),
            )
        };
        let (x, its) = minimise(vec![0.0; 5], f, |_| {}, FREE, 1e-10);
        let residual: Vec<f64> = a(&x).iter().zip(&b).map(|(p, q)| p - q).collect();
        assert!(norm(&residual) < 1e-10, "{residual:?}");
        assert!(its < 20, "{its} iterations");
    }

    /// The Rosenbrock valley defeats steepest descent; L-BFGS gets to (1, 1).
    #[test]
    fn minimises_rosenbrock() {
        let f = |x: &[f64]| {
            let (a, b) = (x[0], x[1]);
            (
                (1.0 - a).powi(2) + 100.0 * (b - a * a).powi(2),
                vec![
                    -2.0 * (1.0 - a) - 400.0 * a * (b - a * a),
                    200.0 * (b - a * a),
                ],
            )
        };
        let (x, its) = minimise(vec![-1.2, 1.0], f, |_| {}, FREE, 1e-8);
        assert!(
            (x[0] - 1.0).abs() < 1e-6 && (x[1] - 1.0).abs() < 1e-6,
            "{x:?} after {its}"
        );
        assert!(its < 100, "{its} iterations");
    }

    /// With the minimum outside a box, the projected iterations end on the box's face.
    #[test]
    fn stays_within_bounds() {
        let f = |x: &[f64]| {
            let target = [2.0, -3.0];
            (
                (x[0] - target[0]).powi(2) + (x[1] - target[1]).powi(2),
                vec![2.0 * (x[0] - target[0]), 2.0 * (x[1] - target[1])],
            )
        };
        let project = |x: &mut [f64]| x.iter_mut().for_each(|v| *v = v.clamp(-1.0, 1.0));
        let (x, _) = minimise(vec![0.0, 0.0], f, project, FREE, 1e-12);
        assert!(
            (x[0] - 1.0).abs() < 1e-9 && (x[1] + 1.0).abs() < 1e-9,
            "{x:?}"
        );
    }

    /// The first step, with nothing to scale it by, changes no unknown by more than the
    /// initial limit; later ones by no more than the maximum.
    #[test]
    fn steps_within_their_limits() {
        let mut lbfgs = Lbfgs::new(3);
        let x = [0.0, 0.0];
        let g = [4.0, -1.0];
        let limits = StepLimits {
            initial: 0.02,
            max: 0.1,
        };
        let step = line_search(
            &mut lbfgs,
            &x,
            1.0,
            &g,
            None,
            limits,
            |_| {},
            |t| Ok(1.0 + dot(&g, t)),
        )
        .unwrap()
        .unwrap();
        assert!((step.x[0] + 0.02).abs() < 1e-15 && (step.x[1] - 0.005).abs() < 1e-15);
        // A pair that says the curvature is tiny asks for a huge step: capped at max.
        assert!(lbfgs.update(vec![1.0, 0.0], vec![1e-6, 0.0]));
        let step = line_search(
            &mut lbfgs,
            &x,
            1.0,
            &g,
            None,
            limits,
            |_| {},
            |t| Ok(1.0 + dot(&g, t)),
        )
        .unwrap()
        .unwrap();
        let largest = step.x.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        assert!((largest - 0.1).abs() < 1e-12, "{:?}", step.x);
    }

    #[test]
    fn keeps_only_pairs_with_positive_curvature() {
        let mut lbfgs = Lbfgs::new(2);
        assert!(!lbfgs.update(vec![1.0, 0.0], vec![-1.0, 0.0]));
        assert!(!lbfgs.update(vec![1.0, 0.0], vec![0.0, 1.0]));
        for i in 0..3 {
            assert!(lbfgs.update(vec![1.0, i as f64], vec![2.0, i as f64]));
        }
        assert_eq!(lbfgs.len(), 2);
        // Without pairs the direction is steepest descent.
        assert_eq!(Lbfgs::new(2).direction(&[1.0, -2.0], None), vec![-1.0, 2.0]);
        // And with a preconditioner, the preconditioned steepest descent.
        assert_eq!(
            Lbfgs::new(2).direction(&[1.0, -2.0], Some(&[0.5, 2.0])),
            vec![-0.5, 4.0]
        );
    }

    /// The oldest pair goes first; NaN or non-positive curvature, or no memory, keeps none;
    /// with one pair and a preconditioner the direction is the two-loop recursion by hand.
    #[test]
    fn remembers_the_newest_pairs_and_scales_by_the_preconditioner() {
        let pairs = [
            (vec![1.0, 0.0], vec![2.0, 0.5]),
            (vec![0.0, 1.0], vec![0.3, 3.0]),
            (vec![1.0, 1.0], vec![1.5, 1.0]),
        ];
        let mut three = Lbfgs::new(2);
        for (s, y) in pairs.clone() {
            assert!(three.update(s, y));
        }
        let mut last_two = Lbfgs::new(2);
        for (s, y) in pairs[1..].iter().cloned() {
            last_two.update(s, y);
        }
        let g = [0.7, -1.1];
        assert_eq!(three.direction(&g, None), last_two.direction(&g, None));
        assert!(!three.update(vec![f64::NAN, 1.0], vec![1.0, 1.0]));
        assert!(!Lbfgs::new(0).update(vec![1.0], vec![1.0]));

        let (s, y) = (vec![1.0, 2.0], vec![3.0, 1.0]);
        let p = [0.5, 4.0];
        let mut one = Lbfgs::new(3);
        one.update(s.clone(), y.clone());
        let rho = 1.0 / dot(&y, &s);
        let a = rho * dot(&s, &g);
        let q: Vec<f64> = g.iter().zip(&y).map(|(g, y)| g - a * y).collect();
        let py: Vec<f64> = y.iter().zip(&p).map(|(y, p)| y * p).collect();
        let gamma = dot(&s, &y) / dot(&y, &py);
        let r: Vec<f64> = q.iter().zip(&p).map(|(q, p)| gamma * p * q).collect();
        let b = rho * dot(&y, &r);
        let want: Vec<f64> = r.iter().zip(&s).map(|(r, s)| -(r + s * (a - b))).collect();
        let got = one.direction(&g, Some(&p));
        assert!(
            got.iter().zip(&want).all(|(a, b)| (a - b).abs() < 1e-12),
            "{got:?}"
        );
    }

    /// A factor w on the preconditioner is a change of variables x = D z, D = diag(√w): the
    /// direction with γ D² P and pairs (s, y) is D times the direction for z, with pairs
    /// (D⁻¹ s, D y), gradient D g and γ P.
    #[test]
    fn preconditioner_factors_change_variables() {
        let pairs = [
            (vec![0.3, -0.2, 0.5], vec![1.1, -0.4, 0.9]),
            (vec![-0.1, 0.4, 0.2], vec![0.2, 1.3, 0.1]),
        ];
        let g = [0.7, -1.1, 0.4];
        let p = [0.5, 4.0, 1.5];
        let w = [1.0, 0.09, 2.5];
        let d: Vec<f64> = w.iter().map(|w: &f64| w.sqrt()).collect();
        let (mut x, mut z) = (Lbfgs::new(5), Lbfgs::new(5));
        for (s, y) in &pairs {
            x.update(s.clone(), y.clone());
            z.update(
                s.iter().zip(&d).map(|(s, d)| s / d).collect(),
                y.iter().zip(&d).map(|(y, d)| y * d).collect(),
            );
        }
        let wp: Vec<f64> = w.iter().zip(&p).map(|(w, p)| w * p).collect();
        let dg: Vec<f64> = g.iter().zip(&d).map(|(g, d)| g * d).collect();
        let direct = x.direction(&g, Some(&wp));
        let changed = z.direction(&dg, Some(&p));
        for ((a, b), d) in direct.iter().zip(&changed).zip(&d) {
            assert!((a - b * d).abs() < 1e-14, "{direct:?} {changed:?}");
        }
    }

    /// Trials that can't be computed count as too long; only if every one fails is that an
    /// error.
    #[test]
    fn shortens_steps_it_cannot_compute() {
        let mut lbfgs = Lbfgs::new(2);
        let g = [1.0];
        let step = line_search(
            &mut lbfgs,
            &[0.0],
            1.0,
            &g,
            None,
            FREE,
            |_| {},
            |t| {
                if t[0] < -0.2 {
                    bail!("unphysical")
                } else {
                    Ok(1.0 + t[0])
                }
            },
        )
        .unwrap()
        .unwrap();
        assert!(
            step.x[0] >= -0.2 && step.x[0] < 0.0 && step.trials > 1,
            "{step:?}"
        );
        let e = line_search(
            &mut lbfgs,
            &[0.0],
            1.0,
            &g,
            None,
            FREE,
            |_| {},
            |_| bail!("no disk"),
        );
        assert!(e.is_err());
        let nan = line_search(
            &mut lbfgs,
            &[0.0],
            1.0,
            &g,
            None,
            FREE,
            |_| {},
            |_| Ok(f64::NAN),
        );
        assert!(nan.unwrap().is_none());
    }

    /// With the minimum on a face of the box, the search slides along the face to it rather
    /// than stopping where the gradient first points out of the box.
    #[test]
    fn slides_along_a_bound() {
        let f = |x: &[f64]| {
            let target = [2.0, 0.3];
            (
                (x[0] - target[0]).powi(2) + 4.0 * (x[1] - target[1]).powi(2),
                vec![2.0 * (x[0] - target[0]), 8.0 * (x[1] - target[1])],
            )
        };
        let project = |x: &mut [f64]| x.iter_mut().for_each(|v| *v = v.clamp(-1.0, 1.0));
        let (x, _) = minimise(vec![0.9, -0.8], f, project, FREE, 1e-12);
        assert!(
            (x[0] - 1.0).abs() < 1e-12 && (x[1] - 0.3).abs() < 1e-6,
            "{x:?}"
        );
    }

    /// A direction from a memory that misleads: the search drops it and goes downhill.
    #[test]
    fn falls_back_to_steepest_descent() {
        // The value rises steeply as soon as x₀ moves, which the gradient (0 along x₀) can't
        // tell; the remembered pair mixes x₀ into the direction, so only steepest descent,
        // along x₁ alone, finds a lower value.
        let value = |t: &[f64]| 1.0 + 1e3 * t[0].abs() + t[1];
        let mut lbfgs = Lbfgs::new(2);
        assert!(lbfgs.update(vec![1.0, 1.0], vec![1.0, 2.0]));
        let d = lbfgs.direction(&[0.0, 1.0], None);
        assert!(d[0] != 0.0 && d[1] < 0.0, "{d:?}");
        let step = line_search(
            &mut lbfgs,
            &[0.0, 0.0],
            1.0,
            &[0.0, 1.0],
            None,
            FREE,
            |_| {},
            |t| Ok(value(t)),
        )
        .unwrap()
        .unwrap();
        assert!(step.value < 1.0 && step.x[0] == 0.0, "{step:?}");
        assert!(step.restarted && lbfgs.is_empty() && step.trials > 8);
    }

    /// A direction that doesn't go downhill and a value that never falls: no step.
    #[test]
    fn gives_up_when_nothing_lowers_the_value() {
        let mut lbfgs = Lbfgs::new(2);
        let step = line_search(
            &mut lbfgs,
            &[0.0],
            1.0,
            &[1.0],
            None,
            FREE,
            |_| {},
            |_| Ok(2.0),
        );
        assert!(step.unwrap().is_none());
        assert!(
            line_search(
                &mut lbfgs,
                &[0.0],
                1.0,
                &[0.0],
                None,
                FREE,
                |_| {},
                |_| Ok(0.0)
            )
            .unwrap()
            .is_none()
        );
    }

    use crate::config::{Config, MisfitKind};
    use crate::model::Model;

    /// A small survey: one explosion, two receivers on the free surface.
    fn survey(layers: &str, extra: &str) -> Config {
        Config::from_toml(&format!(
            r#"
            [grid]
            points = [16, 16, 14]
            spacing = 10.0
            [time]
            duration = 0.03
            dt = 1e-3
            [boundary]
            pml_width = 4
            {layers}
            [[source]]
            position = [80.0, 80.0, 60.0]
            explosion = 1e12
            wavelet = {{ kind = "ricker", frequency = 25.0 }}
            [receivers]
            components = ["vx", "vz"]
            [[receivers.point]]
            position = [60.0, 80.0, 0.0]
            [[receivers.point]]
            position = [100.0, 80.0, 0.0]
            [inversion]
            observed = "observed"
            {extra}
            "#
        ))
        .unwrap()
    }

    const SOLID: &str = "[[model.layer]]\ntop = 0.0\nvp = 2000.0\nvs = 1300.0\nrho = 2000.0";

    fn start(cfg: &Config) -> Model {
        Model::from_config(Simulation::grid(cfg).unwrap(), &cfg.model, Path::new(".")).unwrap()
    }

    /// The traces the configuration's own model records, shaped and timed as run writes them.
    fn recorded(cfg: &Config) -> Observed {
        let mut sim = Simulation::new(cfg, Path::new(".")).unwrap();
        sim.run(None, false).unwrap();
        let r = &sim.receivers;
        Observed {
            shapes: Some(vec![vec![r.nodes.len(), r.samples]; r.components.len()]),
            times: Some(sim.sample_times().iter().map(|&t| t as f32).collect()),
            traces: sim.receivers.traces,
        }
    }

    fn setup(cfg: &Config) -> Result<Inversion> {
        let settings = cfg.inversion.clone().unwrap();
        Inversion::new(cfg, start(cfg), vec![recorded(cfg)], &settings)
    }

    /// vs stays zero in a water layer; the time step and the layers are fixed at the start's;
    /// the illumination preconditions every unknown.
    #[test]
    fn keeps_fluids_fluid_and_the_time_step_and_layers_fixed() {
        let water = "[[model.layer]]\ntop = 0.0\nvp = 1500.0\nvs = 0.0\nrho = 1000.0\n\
                     [[model.layer]]\ntop = 30.0\nvp = 2000.0\nvs = 1300.0\nrho = 2000.0";
        let cfg = survey(water, r#"parameters = ["vs", "vp"]"#);
        let inv = setup(&cfg).unwrap();
        let interior = 8 * 8 * 10; // 16 - 2·4 by 16 - 2·4 by 14 - 4 nodes
        let fluid = 8 * 8 * 3; // z = 0, 10, 20 m
        assert_eq!(inv.unknowns(), (interior - fluid) + interior);
        let shot = &inv.shots[0].config;
        assert_eq!(shot.time.dt, Some(1e-3));
        assert_eq!(shot.boundary.pml_vp_max, Some(2000.0));
        let p = inv.preconditioner.as_ref().unwrap();
        assert_eq!(p.len(), inv.unknowns());
        assert!(p.iter().all(|&v| v.is_finite() && v > 0.0));
        // Pushed hard, vs in the fluid stays 0 and the model stays valid.
        let mut x = vec![0.5; inv.unknowns()];
        inv.project(&mut x);
        let m = inv.model_at(&x);
        assert!(m.vs[..16 * 16 * 3].iter().all(|&v| v == 0.0));
        m.validate().unwrap();
    }

    /// Stages run in turn: each measures the misfit its own way, on its own filtered traces,
    /// and starts L-BFGS afresh; iterations say which stage they belong to.
    #[test]
    fn stages_change_the_misfit_and_restart_lbfgs() {
        let truth = survey(
            "[[model.layer]]\ntop = 0.0\nvp = 2200.0\nvs = 1400.0\nrho = 2000.0",
            "",
        );
        let cfg = survey(
            SOLID,
            "parameters = [\"vs\"]\n\
             [[inversion.stage]]\nmax_frequency = 30.0\niterations = 1\nmisfit = \"traveltime\"\n\
             [[inversion.stage]]\niterations = 1",
        );
        let settings = cfg.inversion.clone().unwrap();
        let mut inv = Inversion::new(&cfg, start(&cfg), vec![recorded(&truth)], &settings).unwrap();
        assert_eq!(inv.stages().len(), 2);
        assert_eq!(inv.stages()[1].misfit, Some(MisfitKind::L2));
        // The first stage measures time shifts of filtered traces.
        let shot = &inv.shots[0];
        let traveltime = Misfit::new(
            Measure::of(&inv.stages()[0], &settings, inv.source_frequency),
            &shot.observed,
            &shot.components,
            shot.samples,
            shot.interval,
        )
        .unwrap();
        let mut sim = Simulation::with_model(&shot.config, inv.model.clone()).unwrap();
        sim.run(None, false).unwrap();
        assert_eq!(
            inv.misfit(),
            traveltime.value(&sim.receivers.traces).unwrap()
        );
        let it = inv.iterate().unwrap().expect("a step");
        assert_eq!(it.stage, 0);
        let n = inv.unknowns();
        assert!(inv.lbfgs.update(vec![1.0; n], vec![1.0; n]));

        inv.begin_stage(1).unwrap();
        assert_eq!((inv.stage(), inv.lbfgs.len()), (1, 0));
        let l2 = Misfit::new(
            Measure::l2(true),
            &inv.shots[0].observed,
            &inv.shots[0].components,
            inv.shots[0].samples,
            inv.shots[0].interval,
        )
        .unwrap();
        assert_eq!(inv.misfit_of(&inv.model).unwrap(), inv.misfit());
        let mut sim = Simulation::with_model(&inv.shots[0].config, inv.model.clone()).unwrap();
        sim.run(None, false).unwrap();
        assert_eq!(inv.misfit(), l2.value(&sim.receivers.traces).unwrap());
        assert_eq!(inv.stage_misfits.len(), 2);
        let it = inv.iterate().unwrap().expect("a step");
        assert_eq!(it.stage, 1);
        assert!(inv.begin_stage(2).is_err());

        // Iterations go with the stages when there are any.
        let both = survey(SOLID, "iterations = 3\n[[inversion.stage]]\niterations = 1");
        assert!(setup(&both).is_err());
        // A corner above the traces' Nyquist frequency is refused before anything runs.
        let high = survey(
            SOLID,
            "[[inversion.stage]]\niterations = 1\n[[inversion.stage]]\nmax_frequency = 900.0\niterations = 1",
        );
        let e = setup(&high).err().expect("refused").to_string();
        assert!(e.contains("Nyquist"), "{e}");
    }

    /// `[inversion.scale]` multiplies each parameter's part of the preconditioner, and alone
    /// makes one; with every factor 1 nothing changes.
    #[test]
    fn scale_multiplies_each_parameters_preconditioner() {
        let plain = setup(&survey(SOLID, r#"parameters = ["vs", "rho"]"#)).unwrap();
        let scaled = setup(&survey(
            SOLID,
            "parameters = [\"vs\", \"rho\"]\n[inversion.scale]\nrho = 0.25\nvp = 7.0",
        ))
        .unwrap();
        let (p, q) = (
            plain.preconditioner.as_ref().unwrap(),
            scaled.preconditioner.as_ref().unwrap(),
        );
        let half = plain.unknowns() / 2;
        assert_eq!(p[..half], q[..half]);
        assert!(
            p[half..]
                .iter()
                .zip(&q[half..])
                .all(|(a, b)| *b == 0.25 * a)
        );
        let alone = setup(&survey(
            SOLID,
            "parameters = [\"vs\", \"rho\"]\nprecondition = false\n[inversion.scale]\nrho = 0.25",
        ))
        .unwrap();
        let r = alone.preconditioner.as_ref().unwrap();
        assert!(r[..half].iter().all(|&v| v == 1.0) && r[half..].iter().all(|&v| v == 0.25));
        let none = setup(&survey(SOLID, "precondition = false")).unwrap();
        assert!(none.preconditioner.is_none());
    }

    /// Sending from every receiver component at once with the signs of a Hadamard matrix's
    /// rows, as many runs as components, gives on average exactly the energy of each
    /// component sent alone, summed: the cross terms cancel. Random signs only cancel them on
    /// average; the inversion then divides by both sides of the illumination.
    #[test]
    fn the_receivers_side_adds_each_receivers_energy() {
        let cfg = survey(SOLID, "receiver_runs = 3\nparameters = [\"vs\"]");
        let inv = setup(&cfg).unwrap();
        // Two receivers, vx and vz: four sources in all, (receiver, component) → 2r + c.
        let hadamard = [
            [1.0, 1.0, 1.0, 1.0],
            [1.0, -1.0, 1.0, -1.0],
            [1.0, 1.0, -1.0, -1.0],
            [1.0, -1.0, -1.0, 1.0],
        ];
        let (coded, _) = inv
            .receiver_energy(&inv.model, 4, &mut |k, r, c| hadamard[k][2 * r + c])
            .unwrap();
        let (alone, _) = inv
            .receiver_energy(&inv.model, 4, &mut |k, r, c| {
                if 2 * r + c == k { 1.0 } else { 0.0 }
            })
            .unwrap();
        let max = alone.iter().fold(0.0f64, |m, &v| m.max(v));
        assert!(max > 0.0);
        // Four single runs averaged: a quarter of the sum.
        let worst = coded
            .iter()
            .zip(&alone)
            .fold(0.0f64, |m, (a, b)| m.max((a - 4.0 * b).abs()));
        assert!(worst < 1e-5 * 4.0 * max, "{worst} of {}", 4.0 * max);
        // The inversion's own estimate multiplies the source side.
        let plain = setup(&survey(SOLID, "parameters = [\"vs\"]")).unwrap();
        let (p, q) = (
            plain.preconditioner.as_ref().unwrap(),
            inv.preconditioner.as_ref().unwrap(),
        );
        assert!(inv.receiver_side.is_some() && p != q);
        assert!(q.iter().all(|&v| v.is_finite() && v > 0.0));
    }

    /// λ ≥ 0 survives projection in every case: vp inverted alone (vp held above √2 vs), and
    /// vs listed before vp (vs capped at the projected vp).
    #[test]
    fn projection_keeps_lambda_non_negative() {
        let inv = setup(&survey(SOLID, r#"parameters = ["vp"]"#)).unwrap();
        let mut x = vec![-0.4; inv.unknowns()];
        inv.project(&mut x);
        let m = inv.model_at(&x);
        m.validate().unwrap();
        let n = inv.free[0].1[0];
        assert!(
            (m.vp[n] as f64 - 2f64.sqrt() * 1300.0).abs() < 1e-2,
            "{}",
            m.vp[n]
        );

        let inv = setup(&survey(
            SOLID,
            "parameters = [\"vs\", \"vp\"]\nvp_range = [1500.0, 2100.0]",
        ))
        .unwrap();
        let half = inv.unknowns() / 2;
        let mut x: Vec<f64> = (0..inv.unknowns())
            .map(|i| if i < half { 0.3 } else { 0.15 })
            .collect();
        inv.project(&mut x);
        let m = inv.model_at(&x);
        for &n in &inv.free[0].1 {
            assert!(m.vp[n] <= 2100.0 + 1e-3);
            assert!(2.0 * (m.vs[n] as f64).powi(2) <= (m.vp[n] as f64).powi(2) * (1.0 + 1e-6));
        }
    }

    /// By default vs may fall to half the smallest non-zero vs, not to 0.
    #[test]
    fn vs_stays_above_half_its_smallest() {
        let inv = setup(&survey(SOLID, r#"parameters = ["vs"]"#)).unwrap();
        let mut x = vec![-0.9; inv.unknowns()];
        inv.project(&mut x);
        let m = inv.model_at(&x);
        let n = inv.free[0].1[0];
        assert!((m.vs[n] - 650.0).abs() < 1e-3, "{}", m.vs[n]);
    }

    /// Ranges the time step can't run, empty ranges and a start outside them are refused
    /// before any wave is computed.
    #[test]
    fn refuses_impossible_bounds() {
        let e = |extra: &str| format!("{:#}", setup(&survey(SOLID, extra)).err().unwrap());
        let msg = e("vp_range = [6000.0, 9000.0]");
        assert!(msg.contains("time step allows"), "{msg}");
        let msg = e("vs_range = [1000.0, 900.0]");
        assert!(msg.contains("empty"), "{msg}");
        let msg = e("rho_range = [2500.0, 3000.0]\nparameters = [\"rho\"]");
        assert!(msg.contains("outside"), "{msg}");
        let lambda = "[[model.layer]]\ntop = 0.0\nvp = 2000.0\nvs = 1500.0\nrho = 2000.0";
        let msg = format!("{:#}", setup(&survey(lambda, "")).err().unwrap());
        assert!(msg.contains("λ < 0"), "{msg}");
    }

    /// Observed traces must have the receivers' shape and sampling; a one-source survey reads
    /// shot-000/ if run --shots wrote it, the directory itself if not.
    #[test]
    fn reads_and_checks_the_observed_traces() {
        let cfg = survey(SOLID, "iterations = 1");
        let good = recorded(&cfg);
        let settings = cfg.inversion.clone().unwrap();
        let try_with = |o: Observed| {
            Inversion::new(&cfg, start(&cfg), vec![o], &settings)
                .err()
                .map(|e| format!("{e:#}"))
        };
        assert_eq!(try_with(good.clone()), None);
        let samples = good.traces[0].len() / 2;
        let transposed = Observed {
            shapes: Some(vec![vec![samples, 2]; 2]),
            ..good.clone()
        };
        assert!(try_with(transposed).unwrap().contains("shape"));
        let shifted = Observed {
            times: good
                .times
                .as_ref()
                .map(|t| t.iter().map(|&t| 1.5 * t).collect()),
            ..good.clone()
        };
        assert!(
            try_with(shifted)
                .unwrap()
                .contains("sampled at other times")
        );
        let short = Observed::from(vec![good.traces[0][1..].to_vec(), good.traces[1].clone()]);
        assert!(try_with(short).unwrap().contains("[time] dt"));

        // Laid out as files.
        let dir =
            std::env::temp_dir().join(format!("vibe-rations-observed-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let write = |into: &Path| {
            std::fs::create_dir_all(into).unwrap();
            for (c, t) in cfg.receivers.components.iter().zip(&good.traces) {
                npy::write(&into.join(format!("{}.npy", c.name())), &[2, samples], t).unwrap();
            }
            let times = good.times.clone().unwrap();
            npy::write(&into.join("time.npy"), &[times.len()], &times).unwrap();
        };
        let toml_dir = |cfg: &Config| {
            let mut c = cfg.clone();
            c.inversion.as_mut().unwrap().observed = dir.join("observed");
            c
        };
        write(&dir.join("observed"));
        Inversion::from_config(&toml_dir(&cfg), &dir).unwrap();
        std::fs::remove_dir_all(dir.join("observed")).unwrap();
        write(&dir.join("observed").join(shot_dir(0)));
        Inversion::from_config(&toml_dir(&cfg), &dir).unwrap();
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
