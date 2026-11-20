//! Full waveform inversion: L-BFGS (Nocedal 1980; Liu & Nocedal 1989) on the gradients of
//! the `gradient` module, as Fichtner (2011, ch. 7) describes it for seismic models.
//!
//! The unknowns are relative changes, x = m/m₀ - 1 for each inverted parameter at each free
//! node, so that vp, vs and ρ share a scale; the gradient with respect to x is m₀ ∂χ/∂m. Free
//! nodes are those outside the absorbing layers, where the gradient is exact, and below
//! `fixed_above`; vs stays zero where it starts at zero, so fluids stay fluid.
//!
//! The time step and the layers' damping are fixed at the starting model's (`[time] dt`,
//! `[boundary] pml_vp_max`), so the misfit is a smooth function of the model and the gradient
//! is exactly its derivative. vp is kept below what that time step allows.
//!
//! Each iteration: the L-BFGS direction (the two-loop recursion), scaled so that no value
//! changes by more than `initial_step` when there are no model and gradient changes to learn
//! the scale from, and by more than `max_step` ever; trial models projected onto the bounds;
//! Armijo backtracking with misfit-only forward runs, shortening the step by quadratic
//! interpolation; then the gradient at the accepted model, whose change is remembered if it
//! has positive curvature (sᵀy > 0).

use std::collections::VecDeque;
use std::path::Path;

use anyhow::{Context, Result, bail, ensure};
use serde::Serialize;

use crate::config::{Config, InversionConfig, Parameter, shot_dir};
use crate::fd::Stencil;
use crate::gradient::{balancing_weights, gradient, l2_misfit};
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
/// diagonal `preconditioner`): scaled to `limits`, projected by `project`, shortened until the
/// value (from `value`) falls by at least 10⁻⁴ of what the gradient predicts (Armijo). None if
/// eight tries don't find one.
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
    let mut restarted = false;
    let mut d = lbfgs.direction(g, preconditioner);
    if dot(g, &d) >= 0.0 {
        lbfgs.clear();
        restarted = true;
        d = lbfgs.direction(g, preconditioner);
    }
    let largest = d.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    if largest == 0.0 {
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
    for trial in 1..=8 {
        let mut t: Vec<f64> = x.iter().zip(&d).map(|(a, b)| a + alpha * b).collect();
        project(&mut t);
        let s: Vec<f64> = t.iter().zip(x).map(|(a, b)| a - b).collect();
        // What the gradient predicts for the step actually taken, after projection.
        let predicted = dot(g, &s);
        let ft = value(&t)?;
        if predicted < 0.0 && ft <= f + ARMIJO * predicted {
            return Ok(Some(Step {
                x: t,
                value: ft,
                length: alpha,
                trials: trial,
                restarted,
            }));
        }
        // The minimum of the parabola through f, its slope and ft, kept within [0.1, 0.5] of
        // the step that failed.
        let slope = predicted / alpha;
        let curvature = ft - f - slope * alpha;
        let next = if curvature > 0.0 && slope < 0.0 {
            -slope * alpha * alpha / (2.0 * curvature)
        } else {
            0.5 * alpha
        };
        alpha = next.clamp(0.1 * alpha, 0.5 * alpha);
    }
    Ok(None)
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
/// layers) and observed traces, with the weight of each component in the misfit.
#[derive(Clone, Debug)]
pub struct Shot {
    pub config: Config,
    pub observed: Vec<Vec<f32>>,
    pub weights: Vec<f64>,
}

/// The misfit, its gradient with respect to the unknowns, and the forward wavefields' energy
/// at every interior node, all summed over the shots.
pub struct Evaluation {
    pub misfit: f64,
    pub gradient: Vec<f64>,
    pub illumination: Vec<f64>,
}

/// What an iteration did.
#[derive(Clone, Debug, Serialize)]
pub struct Iteration {
    pub iteration: usize,
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
    /// The water level of the preconditioner, if there is one.
    water_level: Option<f64>,
    preconditioner: Option<Vec<f64>>,
    lbfgs: Lbfgs,
    x: Vec<f64>,
    misfit: f64,
    g: Vec<f64>,
    /// Iterations done.
    pub iteration: usize,
    pub history: Vec<Iteration>,
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
                let shot = if cfg.sources.len() == 1 {
                    dir.clone()
                } else {
                    dir.join(shot_dir(s))
                };
                load_traces(cfg, &shot)
            })
            .collect::<Result<Vec<_>>>()?;
        Self::new(cfg, model, observed, &settings)
    }

    /// Sets up an inversion of `start` against `observed` traces, one set per source of `cfg`
    /// (each source a shot), in the layout of the receivers' traces.
    pub fn new(
        cfg: &Config,
        start: Model,
        observed: Vec<Vec<Vec<f32>>>,
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

        // Hold the time step and the layers at the starting model's.
        let first = Simulation::with_model(cfg, start.clone())?;
        let dt = first.solver.dt;
        let mut fixed = cfg.clone();
        fixed.time.dt = Some(dt);
        fixed.time.courant = None;
        fixed.boundary.pml_vp_max = Some(
            cfg.boundary
                .pml_vp_max
                .unwrap_or_else(|| start.vp_range()[1]),
        );
        let mut shots = Vec::with_capacity(cfg.sources.len());
        for (s, observed) in observed.into_iter().enumerate() {
            let config = fixed.shot(s)?;
            let sim = Simulation::with_model(&config, start.clone())?;
            let r = &sim.receivers;
            ensure!(
                observed.len() == r.components.len()
                    && observed
                        .iter()
                        .all(|o| o.len() == r.nodes.len() * r.samples),
                "shot {s}: the observed traces don't fit {} receivers × {} samples of {} \
                 component(s); observed and synthetic traces must share their sampling",
                r.nodes.len(),
                r.samples,
                r.components.len()
            );
            let weights = if settings.balance {
                balancing_weights(&observed)
            } else {
                vec![1.0; observed.len()]
            };
            shots.push(Shot {
                config,
                observed,
                weights,
            });
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

        // Bounds; vp also below what the time step allows (with a margin for rounding).
        let stencil = Stencil::taylor(cfg.scheme.order)?;
        let vp_limit = 0.999 * stable_dt(&g, 1.0, &stencil) / dt;
        let around = |range: [f64; 2]| [0.5 * range[0], 2.0 * range[1]];
        let bounds = settings
            .parameters
            .iter()
            .map(|p| {
                let [lo, hi] = match p {
                    Parameter::Vp => settings.vp_range.unwrap_or_else(|| {
                        let [lo, hi] = around(start.vp_range());
                        [lo, hi.min(vp_limit)]
                    }),
                    Parameter::Vs => settings.vs_range.unwrap_or_else(|| {
                        let [_, hi] = around(start.vs_range());
                        [0.0, hi]
                    }),
                    Parameter::Rho => settings
                        .rho_range
                        .unwrap_or_else(|| around(start.rho_range())),
                };
                ensure!(lo >= 0.0 && lo < hi, "the {} range is empty", p.name());
                Ok([
                    lo,
                    if *p == Parameter::Vp {
                        hi.min(vp_limit)
                    } else {
                        hi
                    },
                ])
            })
            .collect::<Result<Vec<_>>>()?;

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
            water_level: settings.precondition.then_some(settings.water_level),
            preconditioner: None,
            lbfgs: Lbfgs::new(settings.memory),
            x: vec![0.0; n],
            misfit: 0.0,
            g: vec![0.0; n],
            iteration: 0,
            history: Vec::new(),
        };
        ensure!(
            settings.water_level > 0.0,
            "the water level must be positive"
        );
        let e = inversion.evaluate(&inversion.model)?;
        inversion.misfit = e.misfit;
        inversion.g = e.gradient;
        inversion.preconditioner = inversion.preconditioner_for(&e.illumination);
        Ok(inversion)
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

    /// Projects unknowns onto the bounds: each value into its range, then vs below vp/√2.
    fn project(&self, x: &mut [f64]) {
        let mut offset = 0;
        let vp_at = |n: usize, x: &[f64]| -> f64 {
            // vp as x has it, if vp is inverted, else as it stands.
            let mut o = 0;
            for (p, nodes) in &self.free {
                if *p == Parameter::Vp
                    && let Ok(i) = nodes.binary_search(&n)
                {
                    return self.start.vp[n] as f64 * (1.0 + x[o + i]);
                }
                o += nodes.len();
            }
            self.model.vp[n] as f64
        };
        for ((p, nodes), &[lo, hi]) in self.free.iter().zip(&self.bounds) {
            let start = match p {
                Parameter::Vp => &self.start.vp,
                Parameter::Vs => &self.start.vs,
                Parameter::Rho => &self.start.rho,
            };
            for (i, &n) in nodes.iter().enumerate() {
                let m0 = start[n] as f64;
                let mut v = (m0 * (1.0 + x[offset + i])).clamp(lo, hi);
                if *p == Parameter::Vs {
                    v = v.min(vp_at(n, x) / 2f64.sqrt());
                }
                x[offset + i] = v / m0 - 1.0;
            }
            offset += nodes.len();
        }
    }

    /// The misfit of `model`, summed over the shots: a forward run each.
    pub fn misfit_of(&self, model: &Model) -> Result<f64> {
        let mut total = 0.0;
        for shot in &self.shots {
            let mut sim = Simulation::with_model(&shot.config, model.clone())?;
            sim.run(None, false)?;
            total += l2_misfit(&sim.receivers.traces, &shot.observed, &shot.weights)?.0;
        }
        Ok(total)
    }

    /// The misfit of `model`, its gradient with respect to the unknowns and the forward
    /// wavefields' energy at every node, summed over the shots.
    pub fn evaluate(&self, model: &Model) -> Result<Evaluation> {
        let n = model.grid.interior_len();
        let mut total = 0.0;
        let mut sum = [vec![0.0; n], vec![0.0; n], vec![0.0; n]];
        let mut illumination = vec![0.0; n];
        for shot in &self.shots {
            let mut sim = Simulation::with_model(&shot.config, model.clone())?;
            let run = gradient(
                &mut sim,
                |r| l2_misfit(&r.traces, &shot.observed, &shot.weights),
                None,
            )?;
            total += run.misfit;
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
        })
    }

    /// 1 / (I / I_max + water level) for each unknown, I the energy at its node, if the
    /// inversion preconditions.
    fn preconditioner_for(&self, illumination: &[f64]) -> Option<Vec<f64>> {
        let level = self.water_level?;
        let max = self
            .free
            .iter()
            .flat_map(|(_, nodes)| nodes.iter().map(|&n| illumination[n]))
            .fold(0.0, f64::max);
        if max <= 0.0 {
            return None;
        }
        Some(
            self.free
                .iter()
                .flat_map(|(_, nodes)| nodes.iter().map(|&n| 1.0 / (illumination[n] / max + level)))
                .collect(),
        )
    }

    /// One iteration, or None if no step along the L-BFGS direction (or, after a restart,
    /// down the gradient) lowers the misfit.
    pub fn iterate(&mut self) -> Result<Option<Iteration>> {
        let x = self.x.clone();
        let g = self.g.clone();
        let mut lbfgs = std::mem::take(&mut self.lbfgs);
        let step = {
            let this = &*self;
            line_search(
                &mut lbfgs,
                &x,
                this.misfit,
                &g,
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
        let Evaluation {
            misfit,
            gradient: g_new,
            illumination,
        } = self.evaluate(&model)?;
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

/// The observed traces in `dir`, `<component>.npy` for every component the configuration
/// records, each (receivers, samples).
pub fn load_traces(cfg: &Config, dir: &Path) -> Result<Vec<Vec<f32>>> {
    cfg.receivers
        .components
        .iter()
        .map(|c| {
            let path = dir.join(format!("{}.npy", c.name()));
            let (_, data) =
                npy::read(&path).with_context(|| format!("reading {}", path.display()))?;
            Ok(data)
        })
        .collect()
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
}
