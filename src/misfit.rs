//! Misfits between synthetic and observed traces, each with its derivative with respect to
//! every synthetic sample: the adjoint source [`crate::gradient::gradient`] takes.
//!
//! Every misfit here is a function of the sampled traces with its exact derivative, so the
//! gradient stays that of the discrete run and finite differences check it. The traces can be
//! low-pass filtered first, both synthetic and observed, for frequency continuation (Bunks et
//! al. 1995; Fichtner et al. 2013): the misfit is then χ(F s, F d) and its derivative Fᵀ ∂χ.

use std::sync::Arc;

use anyhow::{Context, Result, ensure};
use rayon::prelude::*;

use crate::config::{Config, InversionConfig, MisfitKind, StageConfig};
use crate::gradient::{balancing_weights, l2_misfit};
use crate::receiver::Component;

/// How a stage measures the misfit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Measure {
    pub kind: MisfitKind,
    /// The low-pass filter's corner, Hz, if the traces are filtered.
    pub max_frequency: Option<f64>,
    pub filter_order: usize,
    /// The largest time shift the `traveltime` misfit looks for, s; any by default.
    pub max_shift: Option<f64>,
    /// For `l2`: weigh each component by 1/(its observed energy, after filtering).
    pub balance: bool,
    /// For `phase`: the data's dominant frequency, Hz, which sets the time–frequency cells.
    pub dominant_frequency: Option<f64>,
}

impl Measure {
    /// What `stage` of an inversion with `settings` measures, of data from sources of
    /// dominant frequency `source_frequency`: the stage's corner frequency, if it filters, is
    /// then the data's dominant frequency.
    pub fn of(stage: &StageConfig, settings: &InversionConfig, source_frequency: f64) -> Self {
        Self {
            kind: stage.misfit.unwrap_or(settings.misfit),
            max_frequency: stage.max_frequency,
            filter_order: settings.filter_order,
            max_shift: settings.max_shift,
            balance: settings.balance,
            dominant_frequency: Some(stage.max_frequency.unwrap_or(source_frequency)),
        }
    }

    /// The highest dominant frequency of `cfg`'s sources' wavelets: the phase misfit's scale
    /// on traces as recorded.
    pub fn source_frequency(cfg: &Config) -> Result<f64> {
        Ok(cfg
            .sources
            .iter()
            .map(|s| Ok(s.wavelet.wavelet()?.dominant_frequency()))
            .collect::<Result<Vec<f64>>>()?
            .into_iter()
            .fold(0.0, f64::max))
    }

    /// The plain L2 misfit of the traces as recorded, with or without balancing.
    pub fn l2(balance: bool) -> Self {
        Self {
            kind: MisfitKind::L2,
            max_frequency: None,
            filter_order: 4,
            max_shift: None,
            balance,
            dominant_frequency: None,
        }
    }
}

/// A digital Butterworth low-pass (Butterworth 1930) applied forwards, then backwards, with
/// zero initial conditions: zero phase, and the magnitude response squared, so the corner
/// frequency is where the amplitude halves and the roll-off is twice the order's.
///
/// Each pass is a causal recursive filter started at rest, which on N samples is the lower
/// triangular Toeplitz matrix H of its impulse response. Reversing time turns H into
/// R H R = Hᵀ, so the forward–backward filter F = R H R H = Hᵀ H is symmetric: it is its own
/// transpose, and the adjoint source of χ(F s, F d) is F applied to ∂χ/∂(F s).
///
/// The forward pass rings on after the last sample, so each trace is extended with zeros for
/// as long as that ringing lasts, filtered, and cut back: Pᵀ F P, still symmetric, and exact
/// for a trace at rest by its end. Without the zeros the last ~3/f_c of every trace came out
/// distorted. A trace at rest at its start needs no zeros there: the forward pass runs first.
/// Steady-state initial conditions or mirrored padding, as some libraries use, would make the
/// filter unsymmetric.
#[derive(Clone, Debug)]
pub struct Lowpass {
    /// Second-order sections (b0, b1, b2, a1, a2), with a first-order one as (b0, b1, 0, a1, 0)
    /// for odd orders.
    sections: Vec<[f64; 5]>,
    /// Zeros appended to each trace before filtering.
    pad: usize,
}

impl Lowpass {
    /// A Butterworth low-pass of `order` with its corner (each pass's -3 dB point) at `corner`
    /// Hz, for samples `interval` seconds apart: the analog prototype's poles, prewarped and
    /// mapped by the bilinear transform.
    pub fn butterworth(corner: f64, interval: f64, order: usize) -> Result<Self> {
        ensure!(interval > 0.0, "the sample interval must be positive");
        ensure!(
            (1..=16).contains(&order),
            "the filter order must be 1 to 16"
        );
        let nyquist = 0.5 / interval;
        // Near the Nyquist frequency the filter barely filters, and its poles approach the
        // unit circle, so their ringing (and the zeros below) would grow without bound.
        ensure!(
            corner > 0.0 && corner <= 0.8 * nyquist,
            "the corner frequency {corner} Hz must lie between 0 and 0.8 of the traces' \
             Nyquist frequency, {nyquist} Hz"
        );
        let k = (std::f64::consts::PI * corner * interval).tan();
        let mut sections = Vec::with_capacity(order.div_ceil(2));
        for i in 0..order / 2 {
            // s² + 2ζ s + 1 for each pair of conjugate poles, ζ = sin(π(2i + 1)/(2n)).
            let a = 2.0 * (std::f64::consts::PI * (2 * i + 1) as f64 / (2 * order) as f64).sin();
            let norm = 1.0 / (1.0 + a * k + k * k);
            let b0 = k * k * norm;
            sections.push([
                b0,
                2.0 * b0,
                b0,
                2.0 * (k * k - 1.0) * norm,
                (1.0 - a * k + k * k) * norm,
            ]);
        }
        if order % 2 == 1 {
            let norm = 1.0 / (1.0 + k);
            sections.push([k * norm, k * norm, 0.0, (k - 1.0) * norm, 0.0]);
        }
        // The impulse response decays as rⁿ, r the largest pole radius of the digital filter
        // (z² + a1 z + a2 = 0 for a section, z = -a1 for the first-order one): append zeros
        // until rⁿ is below 10⁻⁶, with a margin of 10%. Far below the Nyquist frequency that is
        // the analog prototype's exp(-2π f_c sin(π/2n) t); near it, the bilinear transform's
        // warping slows the decay.
        let radius = sections
            .iter()
            .map(|&[_, _, b2, a1, a2]| {
                if b2 == 0.0 && a2 == 0.0 {
                    a1.abs()
                } else if a1 * a1 < 4.0 * a2 {
                    a2.sqrt()
                } else {
                    let root = (a1 * a1 - 4.0 * a2).sqrt();
                    (0.5 * (-a1 + root)).abs().max((0.5 * (-a1 - root)).abs())
                }
            })
            .fold(0.0f64, f64::max);
        let pad = (1.1 * 1e6f64.ln() / -radius.ln()).ceil() as usize;
        Ok(Self { sections, pad })
    }

    /// Filters `x` in place, forwards then backwards, as if zeros followed it.
    pub fn apply(&self, x: &mut [f64]) {
        let mut y = Vec::with_capacity(x.len() + self.pad);
        y.extend_from_slice(x);
        y.resize(x.len() + self.pad, 0.0);
        self.causal(&mut y);
        y.reverse();
        self.causal(&mut y);
        y.reverse();
        x.copy_from_slice(&y[..x.len()]);
    }

    /// One causal pass from rest, section by section (transposed direct form II).
    fn causal(&self, x: &mut [f64]) {
        for &[b0, b1, b2, a1, a2] in &self.sections {
            let (mut z1, mut z2) = (0.0, 0.0);
            for v in x.iter_mut() {
                let input = *v;
                let out = b0 * input + z1;
                z1 = b1 * input - a1 * out + z2;
                z2 = b2 * input - a2 * out;
                *v = out;
            }
        }
    }

    /// Filters every trace of a component's array, receiver-major with `samples` per trace.
    pub fn apply_traces(&self, data: &[f32], samples: usize) -> Vec<f64> {
        let mut out: Vec<f64> = data.iter().map(|&v| v as f64).collect();
        if samples > 0 {
            out.chunks_mut(samples).for_each(|t| self.apply(t));
        }
        out
    }
}

/// One shot's misfit as a [`Measure`] takes it: the observed traces, filtered if it filters,
/// and the weights it gives them. Traces are laid out as [`crate::receiver::Receivers::traces`]
/// keeps them: one array per component, receiver-major, `samples` per trace.
#[derive(Clone, Debug)]
pub struct Misfit {
    measure: Measure,
    filter: Option<Lowpass>,
    samples: usize,
    interval: f64,
    /// The observed traces, filtered if the measure filters; shared with the caller if not.
    observed: Arc<Vec<Vec<f32>>>,
    /// `l2`: each component's weight. Otherwise each trace's (component-major), 0 for one
    /// too quiet to measure.
    weights: Vec<f64>,
    /// `traveltime`: the largest lag sought, in samples.
    max_lag: usize,
    /// `phase`: the time–frequency cells.
    gabor: Option<Gabor>,
}

/// Velocity traces quieter than this fraction of the loudest velocity trace's energy (and
/// pressure traces than the loudest pressure trace's) aren't measured by the trace-by-trace
/// misfits: a shift or a phase of noise would only mislead.
const QUIET: f64 = 1e-4;

impl Misfit {
    /// Prepares `measure` of synthetic traces against `observed` ones, sampled every
    /// `interval` seconds, `samples` per trace, of `components`.
    pub fn new(
        measure: Measure,
        observed: &Arc<Vec<Vec<f32>>>,
        components: &[Component],
        samples: usize,
        interval: f64,
    ) -> Result<Self> {
        ensure!(
            observed.len() == components.len(),
            "{} observed components for {} recorded",
            observed.len(),
            components.len()
        );
        ensure!(samples > 0, "the traces have no samples");
        for (c, o) in components.iter().zip(observed.iter()) {
            ensure!(
                o.len().is_multiple_of(samples),
                "{} holds {} values, not a whole number of {samples}-sample traces",
                c.name(),
                o.len()
            );
        }
        let filter = measure
            .max_frequency
            .map(|f| Lowpass::butterworth(f, interval, measure.filter_order))
            .transpose()?;
        let observed = match &filter {
            Some(f) => Arc::new(
                observed
                    .iter()
                    .map(|o| {
                        f.apply_traces(o, samples)
                            .into_iter()
                            .map(|v| v as f32)
                            .collect()
                    })
                    .collect(),
            ),
            None => Arc::clone(observed),
        };
        let weights = match measure.kind {
            MisfitKind::L2 if measure.balance => balancing_weights(&observed, components),
            MisfitKind::L2 => vec![1.0; observed.len()],
            MisfitKind::Traveltime | MisfitKind::Phase => {
                trace_weights(&observed, components, samples)
            }
        };
        if let Some(shift) = measure.max_shift {
            ensure!(shift > 0.0, "max_shift must be positive");
        }
        let max_lag = if measure.kind == MisfitKind::Traveltime {
            ensure!(
                samples >= 3,
                "the traveltime misfit needs at least three samples per trace"
            );
            match measure.max_shift {
                Some(shift) => ((shift / interval).round() as usize).clamp(1, samples - 1),
                None => samples - 1,
            }
        } else {
            0
        };
        let gabor = match measure.kind {
            MisfitKind::Phase => {
                let f0 = measure
                    .dominant_frequency
                    .context("the phase misfit needs the data's dominant frequency")?;
                Some(Gabor::new(
                    f0,
                    measure.max_frequency.unwrap_or(2.5 * f0),
                    interval,
                )?)
            }
            _ => None,
        };
        Ok(Self {
            measure,
            filter,
            samples,
            interval,
            observed,
            weights,
            max_lag,
            gabor,
        })
    }

    pub fn measure(&self) -> &Measure {
        &self.measure
    }

    /// The misfit of `synthetic` traces.
    pub fn value(&self, synthetic: &[Vec<f32>]) -> Result<f64> {
        Ok(self.compute(synthetic, false)?.0)
    }

    /// The misfit of `synthetic` traces and its derivative with respect to every sample: the
    /// adjoint source.
    pub fn evaluate(&self, synthetic: &[Vec<f32>]) -> Result<(f64, Vec<Vec<f32>>)> {
        self.compute(synthetic, true)
    }

    /// How strongly the misfit weighs each trace (per component, per receiver) in the
    /// Gauss–Newton Hessian: the square of its adjoint source's size per unit change of the
    /// trace. The L2 misfit's is its weight. The time shift's derivative is about the
    /// observed slope over its energy, ∂T/∂s ≈ ḋ/Σḋ², so its weight is 1/Σḋ²; the phase
    /// misfit compares unit phasors, about ‖δs‖²/E per trace (E the observed energy), divided
    /// by the component's number of measured traces as its normalisation does. Quiet traces,
    /// which the trace-by-trace misfits leave out, weigh 0. Built from the observed traces only.
    pub fn receiver_weights(&self) -> Vec<Vec<f64>> {
        let n = self.samples;
        let mut first = 0;
        self.observed
            .iter()
            .enumerate()
            .map(|(c, o)| {
                let traces = o.len() / n;
                let active = match self.measure.kind {
                    MisfitKind::L2 => 1,
                    _ => self.weights[first..first + traces]
                        .iter()
                        .filter(|&&w| w > 0.0)
                        .count()
                        .max(1),
                };
                let w: Vec<f64> = o
                    .chunks(n)
                    .enumerate()
                    .map(|(r, d)| match self.measure.kind {
                        MisfitKind::L2 => self.weights[c],
                        kind => {
                            let w = self.weights[first + r];
                            let energy: f64 = if kind == MisfitKind::Traveltime {
                                d.windows(3)
                                    .map(|x| (0.5 * (x[2] as f64 - x[0] as f64)).powi(2))
                                    .sum()
                            } else {
                                d.iter().map(|&x| (x as f64).powi(2)).sum::<f64>() * active as f64
                            };
                            if w > 0.0 && energy > 0.0 {
                                w / energy
                            } else {
                                0.0
                            }
                        }
                    })
                    .collect();
                if self.measure.kind != MisfitKind::L2 {
                    first += traces;
                }
                w
            })
            .collect()
    }

    /// Bytes it holds besides the observed traces it shares with the caller: the filtered
    /// observed traces if it filters, the weights, and the phase misfit's kernels.
    pub fn memory_bytes(&self) -> usize {
        let filtered = if self.filter.is_some() {
            4 * self.observed.iter().map(Vec::len).sum::<usize>()
        } else {
            0
        };
        filtered + 8 * self.weights.len() + self.gabor.as_ref().map_or(0, Gabor::memory_bytes)
    }

    /// [`Misfit::memory_bytes`] of a misfit `measure` would set up for one shot of
    /// `receivers` receivers recording `components` components, `samples` samples per trace
    /// `interval` seconds apart, worked out beforehand.
    pub fn held_bytes(
        measure: &Measure,
        receivers: usize,
        components: usize,
        samples: usize,
        interval: f64,
    ) -> Result<usize> {
        let traces = receivers * components;
        let filtered = if measure.max_frequency.is_some() {
            4 * traces * samples
        } else {
            0
        };
        let weights = 8 * match measure.kind {
            MisfitKind::L2 => components,
            _ => traces,
        };
        let kernels = match measure.kind {
            MisfitKind::Phase => {
                let f0 = measure
                    .dominant_frequency
                    .context("the phase misfit needs the data's dominant frequency")?;
                Gabor::new(f0, measure.max_frequency.unwrap_or(2.5 * f0), interval)?.memory_bytes()
            }
            _ => 0,
        };
        Ok(filtered + weights + kernels)
    }

    fn compute(&self, synthetic: &[Vec<f32>], derivative: bool) -> Result<(f64, Vec<Vec<f32>>)> {
        ensure!(
            synthetic.len() == self.observed.len(),
            "{} synthetic components for {} observed",
            synthetic.len(),
            self.observed.len()
        );
        for (s, o) in synthetic.iter().zip(self.observed.iter()) {
            ensure!(
                s.len() == o.len(),
                "{} synthetic samples for {} observed",
                s.len(),
                o.len()
            );
        }
        if self.measure.kind == MisfitKind::L2 && self.filter.is_none() {
            return l2_misfit(synthetic, &self.observed, &self.weights);
        }
        let n = self.samples;
        let mut total = 0.0f64;
        let mut residuals = Vec::with_capacity(synthetic.len());
        let mut first_trace = 0;
        for (c, (s, o)) in synthetic.iter().zip(self.observed.iter()).enumerate() {
            let traces = s.len() / n;
            // Each trace's misfit and derivative, and for `phase` the sum of its squared cell
            // weights, by which the component's misfit is then divided.
            let per_trace: Vec<(f64, Vec<f64>, f64)> = s
                .par_chunks(n)
                .zip(o.par_chunks(n))
                .enumerate()
                .map(|(r, (s, o))| {
                    let w = match self.measure.kind {
                        MisfitKind::L2 => self.weights[c],
                        _ => self.weights[first_trace + r],
                    };
                    let mut s: Vec<f64> = s.iter().map(|&v| v as f64).collect();
                    if w == 0.0 {
                        return (0.0, vec![0.0; if derivative { n } else { 0 }], 0.0);
                    }
                    if let Some(f) = &self.filter {
                        f.apply(&mut s);
                    }
                    let (m, mut g, cells) = match (self.measure.kind, &self.gabor) {
                        (MisfitKind::Traveltime, _) => {
                            let (m, g) =
                                traveltime_trace(&s, o, self.max_lag, self.interval, derivative);
                            (m, g, 1.0)
                        }
                        (MisfitKind::Phase, Some(gabor)) => phase_trace(&s, o, gabor, derivative),
                        _ => {
                            let (m, g) = l2_trace(&s, o, derivative);
                            (m, g, 1.0)
                        }
                    };
                    if derivative {
                        g.iter_mut().for_each(|v| *v *= w);
                        if let Some(f) = &self.filter {
                            f.apply(&mut g);
                        }
                    }
                    (w * m, g, w * cells)
                })
                .collect();
            first_trace += traces;
            let scale = if self.measure.kind == MisfitKind::Phase {
                let cells: f64 = per_trace.iter().map(|t| t.2).sum();
                if cells > 0.0 { 1.0 / cells } else { 0.0 }
            } else {
                1.0
            };
            total += scale * per_trace.iter().map(|t| t.0).sum::<f64>();
            residuals.push(if derivative {
                per_trace
                    .into_iter()
                    .flat_map(|t| t.1)
                    .map(|v| (scale * v) as f32)
                    .collect()
            } else {
                Vec::new()
            });
        }
        ensure!(total.is_finite(), "the misfit isn't finite");
        Ok((total, residuals))
    }
}

/// Each trace's weight in the trace-by-trace misfits, component-major: 1, or 0 for a trace
/// quieter than [`QUIET`] times the loudest of its kind (velocity or pressure).
fn trace_weights(observed: &[Vec<f32>], components: &[Component], samples: usize) -> Vec<f64> {
    let energy: Vec<Vec<f64>> = observed
        .iter()
        .map(|o| {
            o.chunks(samples)
                .map(|t| t.iter().map(|&v| (v as f64).powi(2)).sum())
                .collect()
        })
        .collect();
    let pressure = |c: &Component| matches!(c, Component::P);
    let loudest = |p: bool| {
        components
            .iter()
            .zip(&energy)
            .filter(|(c, _)| pressure(c) == p)
            .flat_map(|(_, e)| e)
            .fold(0.0f64, |m, &e| m.max(e))
    };
    let [velocity, pressure_max] = [loudest(false), loudest(true)];
    components
        .iter()
        .zip(&energy)
        .flat_map(|(c, e)| {
            let floor = QUIET * if pressure(c) { pressure_max } else { velocity };
            e.iter()
                .map(move |&e| if e > 0.0 && e >= floor { 1.0 } else { 0.0 })
        })
        .collect()
}

/// ½ Σ (s - d)² over a (filtered) trace, and its derivative s - d.
fn l2_trace(s: &[f64], d: &[f32], derivative: bool) -> (f64, Vec<f64>) {
    let r: Vec<f64> = s.iter().zip(d).map(|(&a, &b)| a - b as f64).collect();
    let m = 0.5 * r.iter().map(|x| x * x).sum::<f64>();
    (m, if derivative { r } else { Vec::new() })
}

/// The cross-correlation time shift T between a synthetic trace s and an observed one d, and
/// ½ T² with its derivative with respect to s (Luo & Schuster 1991, whose Δτ is -T).
///
/// C(k) = Σₙ s[n + k] d[n] peaks at the lag k* (|k*| ≤ `max_lag`) that best moves d onto s, so
/// T > 0 when the synthetic arrives late. A parabola through C(k* - 1), C(k*) and C(k* + 1)
/// refines it: T = (k* + δ) Δt, δ = (a - c) / (2(a - 2b + c)) with (a, b, c) those three
/// values, so |δ| ≤ ½ wherever k* is inside the lags sought. T is continuous where k* jumps
/// (both parabolas then give the same half-sample shift), and differentiable between:
/// ∂C(j)/∂s[m] = d[m - j], and by the quotient rule ∂δ/∂a = (c - b)/D², ∂δ/∂b = (a - c)/D²,
/// ∂δ/∂c = (b - a)/D², D = a - 2b + c. The exact derivative is built from the observed trace
/// (its slope at the shift found); the familiar adjoint source made of the synthetic's own
/// slope (Tromp et al. 2005) is its limit when d is a shifted copy of s.
///
/// At the largest lag sought δ is held within ±½ (with no derivative while held), so T
/// saturates at ±(`max_lag` + ½) Δt instead of extrapolating. A synthetic that correlates
/// nowhere positively with d (a silent one, say) gets that saturated shift too, so a model
/// can't lower the misfit by silencing a trace.
pub fn traveltime_trace(
    s: &[f64],
    d: &[f32],
    max_lag: usize,
    interval: f64,
    derivative: bool,
) -> (f64, Vec<f64>) {
    let n = s.len();
    let max_lag = max_lag.clamp(1, n.saturating_sub(1).max(1)) as isize;
    let c = |k: isize| -> f64 {
        let (from, to) = (0.max(-k), (n as isize).min(n as isize - k));
        (from..to)
            .map(|i| s[(i + k) as usize] * d[i as usize] as f64)
            .sum()
    };
    let mut best = (0, f64::NEG_INFINITY);
    for k in -max_lag..=max_lag {
        let v = c(k);
        if v > best.1 {
            best = (k, v);
        }
    }
    let held = |t: f64| {
        (
            0.5 * t * t,
            if derivative { vec![0.0; n] } else { Vec::new() },
        )
    };
    let (k, b) = best;
    if b <= 0.0 {
        return held((max_lag as f64 + 0.5) * interval);
    }
    let (a, cc) = (c(k - 1), c(k + 1));
    let den = a - 2.0 * b + cc;
    let raw = if den < 0.0 {
        (a - cc) / (2.0 * den)
    } else {
        f64::NAN
    };
    // Only at the largest lags can the parabola open upwards or its top lie beyond ±½.
    if raw.is_nan() || raw.abs() > 0.5 {
        let delta = if k > 0 { 0.5 } else { -0.5 };
        return held((k as f64 + delta) * interval);
    }
    let t = (k as f64 + raw) * interval;
    let m = 0.5 * t * t;
    if !derivative {
        return (m, Vec::new());
    }
    let mut g = vec![0.0; n];
    let d2 = den * den;
    // ∂(½T²)/∂s = T Δt ∂δ/∂s, each lag's share through ∂C(j)/∂s[i] = d[i - j].
    for (j, dd) in [
        (k - 1, (cc - b) / d2),
        (k, (a - cc) / d2),
        (k + 1, (b - a) / d2),
    ] {
        let f = t * interval * dd;
        for (i, gi) in g.iter_mut().enumerate() {
            let src = i as isize - j;
            if (0..n as isize).contains(&src) {
                *gi += f * d[src as usize] as f64;
            }
        }
    }
    (m, g)
}

/// The phase misfit's regularisation: z is compared as z/√(|z|² + ε²), ε this fraction of the
/// observed trace's largest |z|.
const PHASE_EPSILON: f64 = 0.05;
/// Cells where the observed trace's |z| is below this fraction of its largest carry no weight.
const PHASE_CUT: f64 = 0.05;

/// The cells of a time–frequency representation that the phase misfit compares: a Gabor
/// transform (Fichtner et al. 2008) with Gaussian windows of width σ = 1/f₀, of unit norm,
/// every σ/2 and truncated at ±3σ, at frequencies every 1/(4πσ) Hz from the first above zero
/// to 1.2 times the top of the band, f_hi (below the Nyquist frequency, whose cells are real).
/// Each frequency's weight tapers off below f₀/4 and above f_hi.
#[derive(Clone, Debug)]
struct Gabor {
    /// Samples between window centres, and each window's half-width in samples.
    spacing: usize,
    half: usize,
    /// For each frequency, its taper and its kernel (Δ/√2π) g(mΔ) e^(-iωmΔ),
    /// m = -half..=half, as (re, im).
    frequencies: Vec<(f64, Vec<(f64, f64)>)>,
}

impl Gabor {
    fn new(f0: f64, f_hi: f64, interval: f64) -> Result<Self> {
        use std::f64::consts::PI;
        ensure!(
            f0 > 0.0 && f_hi > 0.0,
            "the phase misfit needs positive frequencies"
        );
        let sigma = 1.0 / f0;
        let spacing = ((0.5 * sigma / interval).round() as usize).max(1);
        let half = (3.0 * sigma / interval).ceil() as usize;
        let df = 1.0 / (4.0 * PI * sigma);
        let nyquist = 0.5 / interval;
        let norm = interval / (2.0 * PI).sqrt() * (PI * sigma * sigma).powf(-0.25);
        let mut frequencies = Vec::new();
        for j in 1.. {
            let f = j as f64 * df;
            if f > 1.2 * f_hi || f >= nyquist * (1.0 - 1e-9) {
                break;
            }
            let low = 1.0 - (-(4.0 * f / f0).powi(2)).exp();
            let high = if f <= f_hi {
                1.0
            } else {
                (-10.0 * (f / f_hi - 1.0)).exp()
            };
            let w = 2.0 * PI * f;
            let kernel = (0..=2 * half)
                .map(|i| {
                    let t = (i as f64 - half as f64) * interval;
                    let g = norm * (-0.5 * (t / sigma).powi(2)).exp();
                    let (sin, cos) = (w * t).sin_cos();
                    (g * cos, -g * sin)
                })
                .collect();
            frequencies.push((low * high, kernel));
        }
        ensure!(
            !frequencies.is_empty(),
            "the phase misfit's band holds no frequency below the traces' Nyquist frequency"
        );
        Ok(Self {
            spacing,
            half,
            frequencies,
        })
    }

    /// Bytes its kernels hold.
    fn memory_bytes(&self) -> usize {
        self.frequencies.iter().map(|(_, k)| 16 * k.len()).sum()
    }

    /// For each window centre of a trace of `n` samples, the samples it reaches, and the
    /// kernel index of the first.
    fn windows(&self, n: usize) -> impl Iterator<Item = (std::ops::Range<usize>, usize)> + '_ {
        (0..n).step_by(self.spacing).map(move |c| {
            let from = c.saturating_sub(self.half);
            let to = (c + self.half + 1).min(n);
            (from..to, from + self.half - c)
        })
    }

    /// z for every cell, window by window, frequency by frequency.
    fn transform(&self, x: &[f64]) -> Vec<(f64, f64)> {
        let mut z = Vec::with_capacity(x.len().div_ceil(self.spacing) * self.frequencies.len());
        for (range, first) in self.windows(x.len()) {
            for (_, kernel) in &self.frequencies {
                let (mut re, mut im) = (0.0, 0.0);
                for (&v, &(kr, ki)) in x[range.clone()].iter().zip(&kernel[first..]) {
                    re += v * kr;
                    im += v * ki;
                }
                z.push((re, im));
            }
        }
        z
    }

    /// out[n] += Re Σ a G over every cell: the transform's transpose applied to the complex
    /// coefficients a, one per cell.
    fn add_transpose(&self, a: &[(f64, f64)], out: &mut [f64]) {
        let mut cells = a.iter();
        for (range, first) in self.windows(out.len()) {
            for (_, kernel) in &self.frequencies {
                let &(ar, ai) = cells.next().expect("a coefficient for every cell");
                if ar == 0.0 && ai == 0.0 {
                    continue;
                }
                for (o, &(kr, ki)) in out[range.clone()].iter_mut().zip(&kernel[first..]) {
                    *o += ar * kr - ai * ki;
                }
            }
        }
    }
}

/// The time–frequency phase misfit of a synthetic trace s against an observed one d, before
/// its component's normalisation: ½ Σ W² |u - u₀|² over the cells, its derivative, and Σ W².
///
/// u = z/√(|z|² + ε²) is the synthetic's Gabor coefficient z brought to unit length where it
/// is large, and u₀ likewise the observed one's. Where both are, |u - u₀|² = 2 - 2 cos Δφ ≈
/// Δφ², Fichtner et al.'s (2008) phase difference squared; ε keeps the misfit and its
/// derivative bounded where the synthetic lacks an arrival the observed trace has (|z| → 0),
/// and there is no branch cut at Δφ = ±π. W = taper(f) |z₀|/max|z₀| where that ratio is at
/// least [`PHASE_CUT`], else 0: built from the observed trace alone, so the misfit is one
/// smooth function of s. With q = |z|² + ε² and e = u - u₀, dχ = Re(A dz) with
/// A = W² q^(-½) (ē - Re(ē z) z̄/q), and dz = Σ G ds.
fn phase_trace(s: &[f64], d: &[f32], gabor: &Gabor, derivative: bool) -> (f64, Vec<f64>, f64) {
    let d: Vec<f64> = d.iter().map(|&v| v as f64).collect();
    let z0 = gabor.transform(&d);
    let largest = z0.iter().map(|&(r, i)| r.hypot(i)).fold(0.0f64, f64::max);
    let mut g = if derivative {
        vec![0.0; s.len()]
    } else {
        Vec::new()
    };
    if largest == 0.0 {
        return (0.0, g, 0.0);
    }
    let z = gabor.transform(s);
    let eps2 = (PHASE_EPSILON * largest).powi(2);
    let per_window = gabor.frequencies.len();
    let mut misfit = 0.0;
    let mut weights = 0.0;
    let mut a = Vec::with_capacity(if derivative { z.len() } else { 0 });
    for (i, (&(zr, zi), &(or, oi))) in z.iter().zip(&z0).enumerate() {
        let r = or.hypot(oi) / largest;
        let w2 = if r >= PHASE_CUT {
            (gabor.frequencies[i % per_window].0 * r).powi(2)
        } else {
            0.0
        };
        if w2 == 0.0 {
            if derivative {
                a.push((0.0, 0.0));
            }
            continue;
        }
        weights += w2;
        let q = zr * zr + zi * zi + eps2;
        let q0 = or * or + oi * oi + eps2;
        let (er, ei) = (
            zr / q.sqrt() - or / q0.sqrt(),
            zi / q.sqrt() - oi / q0.sqrt(),
        );
        misfit += 0.5 * w2 * (er * er + ei * ei);
        if derivative {
            // ē = er - i ei and z̄ = zr - i zi; Re(ē z) = er zr + ei zi.
            let p = (er * zr + ei * zi) / q;
            let f = w2 / q.sqrt();
            a.push((f * (er - p * zr), f * (p * zi - ei)));
        }
    }
    if derivative {
        gabor.add_transpose(&a, &mut g);
    }
    (misfit, g, weights)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// The amplitude the forward–backward filter passes at `f` Hz, measured on a long sine.
    fn gain(filter: &Lowpass, f: f64, interval: f64) -> f64 {
        let n = 20_000;
        let mut x: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * f * i as f64 * interval).sin())
            .collect();
        filter.apply(&mut x);
        // Projected on the sine and cosine away from both ends, where the filter has settled,
        // over a whole number of periods.
        let w = 2.0 * PI * f * interval;
        let period = 1.0 / (f * interval);
        let len = ((n / 2) as f64 / period).floor() * period;
        let (a, b) = (n / 4..n / 4 + len as usize).fold((0.0, 0.0), |(s, c), i| {
            (
                s + x[i] * (w * i as f64).sin(),
                c + x[i] * (w * i as f64).cos(),
            )
        });
        2.0 * (a * a + b * b).sqrt() / len
    }

    /// Each pass has the bilinear-transformed Butterworth response,
    /// |H|² = 1 / (1 + (tan(πfΔ) / tan(πf_cΔ))²ⁿ), and the two passes square it: half the
    /// amplitude at the corner.
    #[test]
    fn has_the_butterworth_response_squared() {
        let interval = 1e-3;
        let corner = 20.0;
        for order in [1, 2, 3, 4, 6] {
            let f = Lowpass::butterworth(corner, interval, order).unwrap();
            let tan = |hz: f64| (PI * hz * interval).tan();
            for hz in [2.0, 10.0, 20.0, 30.0, 80.0] {
                let expect = 1.0 / (1.0 + (tan(hz) / tan(corner)).powi(2 * order as i32));
                let got = gain(&f, hz, interval);
                assert!(
                    (got - expect).abs() < 1e-3 * expect.max(1e-3),
                    "order {order}, {hz} Hz: {got} vs {expect}"
                );
            }
        }
    }

    #[test]
    fn keeps_a_low_frequency_pulse_in_place() {
        // Zero phase: a symmetric pulse well inside the band keeps its peak where it was.
        let interval = 1e-3;
        let f = Lowpass::butterworth(40.0, interval, 4).unwrap();
        let mut x: Vec<f64> = (0..1000)
            .map(|i| {
                let t = i as f64 * interval - 0.5;
                let a = (PI * 5.0 * t).powi(2);
                (1.0 - 2.0 * a) * (-a).exp()
            })
            .collect();
        let before = x.clone();
        f.apply(&mut x);
        let peak = |v: &[f64]| (0..v.len()).max_by(|&a, &b| v[a].total_cmp(&v[b])).unwrap();
        assert_eq!(peak(&x), peak(&before));
        let err = x
            .iter()
            .zip(&before)
            .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()));
        assert!(err < 1e-3, "{err}");
    }

    #[test]
    fn is_its_own_transpose() {
        // ⟨F x, y⟩ = ⟨x, F y⟩ for any x, y: the forward–backward filter is symmetric.
        let n = 700;
        let x: Vec<f64> = (0..n).map(|i| ((i * 7919) % 101) as f64 - 50.0).collect();
        let y: Vec<f64> = (0..n).map(|i| ((i * 104729) % 97) as f64 - 48.0).collect();
        for (corner, order) in [(5.0, 4), (60.0, 3), (1.0, 8)] {
            let f = Lowpass::butterworth(corner, 2e-3, order).unwrap();
            let (mut fx, mut fy) = (x.clone(), y.clone());
            f.apply(&mut fx);
            f.apply(&mut fy);
            let a: f64 = fx.iter().zip(&y).map(|(a, b)| a * b).sum();
            let b: f64 = x.iter().zip(&fy).map(|(a, b)| a * b).sum();
            assert!((a - b).abs() <= 1e-12 * a.abs().max(b.abs()), "{a} {b}");
        }
    }

    fn ricker(f0: f64, t0: f64, interval: f64, n: usize) -> Vec<f64> {
        (0..n)
            .map(|i| {
                let a = (PI * f0 * (i as f64 * interval - t0)).powi(2);
                (1.0 - 2.0 * a) * (-a).exp()
            })
            .collect()
    }

    fn to_f32(x: &[f64]) -> Vec<f32> {
        x.iter().map(|&v| v as f32).collect()
    }

    #[test]
    fn measures_a_known_shift() {
        let dt = 1e-3;
        let d = to_f32(&ricker(15.0, 0.2, dt, 500));
        // The synthetic arrives 3.4 samples late; the parabola finds it to a few hundredths.
        let s = ricker(15.0, 0.2 + 3.4 * dt, dt, 500);
        let (m, _) = traveltime_trace(&s, &d, 50, dt, false);
        let t = (2.0 * m).sqrt();
        assert!((t / dt - 3.4).abs() < 0.03, "{}", t / dt);
        let (_, g) = traveltime_trace(&s, &d, 50, dt, true);
        // s(t + ε) ≈ s + εṡ arrives earlier, which lowers T and the misfit: g · ṡ < 0.
        let sdot: f64 = (1..499).map(|i| g[i] * (s[i + 1] - s[i - 1])).sum();
        assert!(sdot < 0.0, "{sdot}");
        // Early by 7.8 samples.
        let s = ricker(15.0, 0.2 - 7.8 * dt, dt, 500);
        let (m, _) = traveltime_trace(&s, &d, 50, dt, false);
        assert!(
            ((2.0 * m).sqrt() / dt - 7.8).abs() < 0.1,
            "{}",
            (2.0 * m).sqrt() / dt
        );
        // A silent synthetic gets the largest shift, with no derivative.
        let (m, g) = traveltime_trace(&vec![0.0; 500], &d, 50, dt, true);
        assert!(((2.0 * m).sqrt() / dt - 50.5).abs() < 1e-9, "{m}");
        assert!(g.iter().all(|&v| v == 0.0));
    }

    /// Beyond the largest lag sought, T saturates at it plus half a sample, continuously.
    #[test]
    fn the_time_shift_saturates_at_the_largest_lag() {
        let dt = 1e-3;
        let d = to_f32(&ricker(10.0, 0.3, dt, 600));
        let shift = |x: f64| {
            let s = ricker(10.0, 0.3 + x * dt, dt, 600);
            (2.0 * traveltime_trace(&s, &d, 8, dt, false).0).sqrt() / dt
        };
        let mut last = 0.0;
        for i in 0..=60 {
            let x = 6.0 + 0.1 * i as f64;
            let t = shift(x);
            assert!(t <= 8.5 + 1e-9 && t >= last - 1e-9, "{x}: {t} after {last}");
            last = t;
        }
        assert!((last - 8.5).abs() < 1e-9, "{last}");
    }

    /// Central differences in f64 along a few directions, between the lag jumps.
    #[test]
    fn the_traveltime_derivative_matches_finite_differences() {
        let dt = 2e-3;
        let n = 300;
        let d = to_f32(&ricker(10.0, 0.25, dt, n));
        let s: Vec<f64> = ricker(10.0, 0.25 + 5.3 * dt, dt, n)
            .iter()
            .zip(ricker(25.0, 0.4, dt, n))
            .map(|(a, b)| a + 0.2 * b)
            .collect();
        let (_, g) = traveltime_trace(&s, &d, 40, dt, true);
        for dir in 0..4 {
            let u: Vec<f64> = (0..n)
                .map(|i| ((i * (dir + 3) * 7919) % 113) as f64 / 56.0 - 1.0)
                .collect();
            let at = |h: f64| {
                let x: Vec<f64> = s.iter().zip(&u).map(|(a, b)| a + h * b).collect();
                traveltime_trace(&x, &d, 40, dt, false).0
            };
            let h = 1e-6;
            let fd = (at(h) - at(-h)) / (2.0 * h);
            let an: f64 = g.iter().zip(&u).map(|(a, b)| a * b).sum();
            assert!(
                (fd - an).abs() <= 1e-6 * an.abs().max(1e-12),
                "{dir}: {fd} vs {an}"
            );
        }
    }

    /// Through a [`Misfit`]: the filter, the weights and the f32 traces, for every kind, with
    /// and without filtering.
    #[test]
    fn misfits_differentiate_their_traces() {
        let dt = 2e-3;
        let n = 250;
        let components = [Component::Vx, Component::P];
        let trace = |f0: f64, t0: f64, a: f64| {
            ricker(f0, t0, dt, n)
                .iter()
                .map(|v| a * v)
                .collect::<Vec<_>>()
        };
        // Two receivers per component; pressure a million times velocity.
        let obs = Arc::new(vec![
            to_f32(&[trace(12.0, 0.2, 1e-3), trace(12.0, 0.26, 8e-4)].concat()),
            to_f32(&[trace(12.0, 0.21, 1e3), trace(12.0, 0.27, 9e2)].concat()),
        ]);
        let syn: Vec<Vec<f32>> = vec![
            to_f32(&[trace(11.0, 0.205, 1.1e-3), trace(12.5, 0.252, 7e-4)].concat()),
            to_f32(&[trace(12.0, 0.214, 1.2e3), trace(11.5, 0.266, 9.5e2)].concat()),
        ];
        for kind in [MisfitKind::L2, MisfitKind::Traveltime, MisfitKind::Phase] {
            for max_frequency in [None, Some(15.0)] {
                let measure = Measure {
                    kind,
                    max_frequency,
                    filter_order: 4,
                    max_shift: Some(0.03),
                    balance: true,
                    dominant_frequency: Some(max_frequency.unwrap_or(12.0)),
                };
                let m = Misfit::new(measure, &obs, &components, n, dt).unwrap();
                let (value, g) = m.evaluate(&syn).unwrap();
                assert_eq!(value, m.value(&syn).unwrap());
                // Down the gradient, scaled per component to its traces' size.
                let u: Vec<Vec<f32>> = g
                    .iter()
                    .zip(&syn)
                    .map(|(g, s)| {
                        let gmax = g.iter().fold(0f32, |a, v| a.max(v.abs()));
                        let smax = s.iter().fold(0f32, |a, v| a.max(v.abs()));
                        g.iter()
                            .map(|v| -v / gmax.max(f32::MIN_POSITIVE) * smax)
                            .collect()
                    })
                    .collect();
                let predicted: f64 = g
                    .iter()
                    .zip(&u)
                    .flat_map(|(g, u)| g.iter().zip(u).map(|(a, b)| *a as f64 * *b as f64))
                    .sum();
                let at = |h: f32| {
                    let x: Vec<Vec<f32>> = syn
                        .iter()
                        .zip(&u)
                        .map(|(s, u)| s.iter().zip(u).map(|(a, b)| a + h * b).collect())
                        .collect();
                    m.value(&x).unwrap()
                };
                let h = 1e-3;
                let fd = (at(h) - at(-h)) / (2.0 * h as f64);
                assert!(
                    (fd - predicted).abs() < 2e-3 * predicted.abs(),
                    "{kind:?}, {max_frequency:?}: {fd} vs {predicted}"
                );
            }
        }
    }

    /// ⟨G x, a⟩ = ⟨x, Gᵀ a⟩ (real parts), and the misfit and its derivative in f64.
    #[test]
    fn the_phase_misfit_differentiates_its_cells() {
        let dt = 1e-3;
        let n = 400;
        let gabor = Gabor::new(12.0, 30.0, dt).unwrap();
        let x: Vec<f64> = ricker(12.0, 0.15, dt, n)
            .iter()
            .zip(ricker(20.0, 0.25, dt, n))
            .map(|(a, b)| a - 0.5 * b)
            .collect();
        let z = gabor.transform(&x);
        let a: Vec<(f64, f64)> = (0..z.len())
            .map(|i| ((i % 7) as f64 - 3.0, (i % 5) as f64 - 2.0))
            .collect();
        let mut back = vec![0.0; n];
        gabor.add_transpose(&a, &mut back);
        let lhs: f64 = z.iter().zip(&a).map(|(z, a)| z.0 * a.0 - z.1 * a.1).sum();
        let rhs: f64 = x.iter().zip(&back).map(|(x, b)| x * b).sum();
        assert!((lhs - rhs).abs() < 1e-12 * lhs.abs(), "{lhs} {rhs}");

        // Identical traces have no misfit, and amplitude hardly matters: ε only weighs in
        // where |z| is near the cut.
        let d = to_f32(&x);
        let same: Vec<f64> = d.iter().map(|&v| v as f64).collect();
        assert_eq!(phase_trace(&same, &d, &gabor, false).0, 0.0);
        let louder: Vec<f64> = same.iter().map(|v| 3.0 * v).collect();
        let (m, _, w) = phase_trace(&louder, &d, &gabor, false);
        // Measured: 1.1e-4, against 0.29 for the shifted synthetic below.
        assert!(m / w < 3e-4, "{}", m / w);

        // A shifted synthetic, with an arrival the observed trace lacks.
        let s: Vec<f64> = ricker(12.0, 0.158, dt, n)
            .iter()
            .zip(ricker(20.0, 0.262, dt, n))
            .zip(ricker(15.0, 0.33, dt, n))
            .map(|((a, b), c)| a - 0.45 * b + 0.3 * c)
            .collect();
        let (m, g, w) = phase_trace(&s, &d, &gabor, true);
        assert!(m / w > 0.1, "{}", m / w);
        for dir in 0..3 {
            let u: Vec<f64> = (0..n)
                .map(|i| ((i * (dir + 2) * 7919) % 89) as f64 / 44.0 - 1.0)
                .collect();
            let at = |h: f64| {
                let y: Vec<f64> = s.iter().zip(&u).map(|(a, b)| a + h * b).collect();
                phase_trace(&y, &d, &gabor, false).0
            };
            let h = 1e-5;
            let fd = (at(h) - at(-h)) / (2.0 * h);
            let an: f64 = g.iter().zip(&u).map(|(a, b)| a * b).sum();
            assert!((fd - an).abs() < 1e-6 * an.abs(), "{dir}: {fd} vs {an}");
        }
    }

    #[test]
    fn without_a_filter_l2_is_the_plain_misfit() {
        let obs = Arc::new(vec![vec![1.0f32, 2.0, 0.5, 0.0]]);
        let syn = vec![vec![0.5f32, 2.5, 0.0, 1.0]];
        let m = Misfit::new(Measure::l2(true), &obs, &[Component::Vz], 2, 1e-3).unwrap();
        let weights = balancing_weights(&obs, &[Component::Vz]);
        assert_eq!(
            m.evaluate(&syn).unwrap(),
            l2_misfit(&syn, &obs, &weights).unwrap()
        );
        // It shares the observed traces, holding only its weight.
        assert_eq!(m.memory_bytes(), 8);
    }

    #[test]
    fn quiet_traces_are_left_out() {
        use Component::{P, Vx, Vy};
        let w = trace_weights(
            &[
                vec![1.0, 0.0, 1e-3, 0.0],
                vec![1e-3, 0.0, 0.0, 0.0],
                vec![1e-3, 0.0, 0.0, 0.0],
            ],
            &[Vx, Vy, P],
            2,
        );
        // vy's trace is quiet next to vx's loudest; pressure is measured against pressure.
        assert_eq!(w, vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0]);
    }

    /// The numbers in the research log's entry on these misfits: the parabola's bias, how far
    /// Tromp et al.'s (2005) adjoint source is from the exact one, the filter's symmetry and
    /// the phase misfit's cells. `cargo test --release --lib misfit_study -- --ignored
    /// --nocapture`.
    #[test]
    #[ignore]
    fn misfit_study() {
        // Parabola bias: true fractional shifts over a sample, at P samples per peak period.
        for p in [5.0, 10.0, 20.0, 50.0] {
            let f0 = 10.0;
            let dt = 1.0 / (p * f0);
            let n = (1.0 / dt) as usize;
            let d = to_f32(&ricker(f0, 0.5, dt, n));
            let mut worst = 0.0f64;
            for i in 0..=40 {
                let x = 3.0 + i as f64 / 40.0 - 0.5;
                let s = ricker(f0, 0.5 + x * dt, dt, n);
                let t = (2.0 * traveltime_trace(&s, &d, 20, dt, false).0).sqrt() / dt;
                worst = worst.max((t - x).abs());
            }
            eprintln!(
                "bias P={p}: {worst:.2e} samples, {:.2e} of a period",
                worst / p
            );
        }
        // Exact derivative against Tromp et al.'s ṡ / ∫ s s̈ (in samples, central differences).
        let dt = 1e-3;
        let n = 600;
        for (label, s, d) in [
            (
                "shifted copy, 3 samples",
                ricker(12.0, 0.303, dt, n),
                ricker(12.0, 0.3, dt, n),
            ),
            (
                "shifted copy, 3.3 samples",
                ricker(12.0, 0.3033, dt, n),
                ricker(12.0, 0.3, dt, n),
            ),
            (
                "different shapes",
                ricker(12.0, 0.303, dt, n)
                    .iter()
                    .zip(ricker(20.0, 0.33, dt, n))
                    .map(|(a, b)| a + 0.3 * b)
                    .collect(),
                ricker(13.0, 0.3, dt, n),
            ),
        ] {
            let d = to_f32(&d);
            let (m, g) = traveltime_trace(&s, &d, 40, dt, true);
            let t = (2.0 * m).sqrt();
            let sdot: Vec<f64> = (0..n)
                .map(|i| {
                    if i == 0 || i == n - 1 {
                        0.0
                    } else {
                        (s[i + 1] - s[i - 1]) / (2.0 * dt)
                    }
                })
                .collect();
            let sdd: f64 = (1..n - 1)
                .map(|i| s[i] * (s[i + 1] - 2.0 * s[i] + s[i - 1]) / (dt * dt))
                .sum::<f64>()
                * dt;
            let ttl: Vec<f64> = sdot.iter().map(|v| t * v / sdd * dt).collect();
            let dot: f64 = g.iter().zip(&ttl).map(|(a, b)| a * b).sum();
            let na = g.iter().map(|v| v * v).sum::<f64>().sqrt();
            let nb = ttl.iter().map(|v| v * v).sum::<f64>().sqrt();
            let diff = g
                .iter()
                .zip(&ttl)
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>()
                .sqrt();
            eprintln!(
                "{label}: T = {:.3} samples; Tromp vs exact: cosine {:.4}, relative difference {:.2e}",
                t / dt,
                dot / (na * nb),
                diff / na
            );
        }
        // The filter's symmetry, as measured.
        let x: Vec<f64> = (0..700).map(|i| ((i * 7919) % 101) as f64 - 50.0).collect();
        let y: Vec<f64> = (0..700)
            .map(|i| ((i * 104729) % 97) as f64 - 48.0)
            .collect();
        let f = Lowpass::butterworth(5.0, 2e-3, 4).unwrap();
        let (mut fx, mut fy) = (x.clone(), y.clone());
        f.apply(&mut fx);
        f.apply(&mut fy);
        let a: f64 = fx.iter().zip(&y).map(|(a, b)| a * b).sum();
        let b: f64 = x.iter().zip(&fy).map(|(a, b)| a * b).sum();
        eprintln!("filter symmetry: {:.1e}", (a - b).abs() / a.abs());
        // The record's end: a 12 Hz pulse ending 0.15 s (0.5/f_c) and 1 s (3/f_c) before it,
        // filtered at 3 Hz, against the same pulse in a record ten times as long; without the
        // zeros (the passes alone) and with them.
        let f = Lowpass::butterworth(3.0, 1e-3, 4).unwrap();
        for (label, len) in [("0.5/fc", 600), ("3/fc", 1450)] {
            let short = ricker(12.0, 0.4, 1e-3, len);
            let mut long = short.clone();
            long.resize(10 * len, 0.0);
            f.apply(&mut long);
            let peak = long.iter().fold(0.0f64, |m, v| m.max(v.abs()));
            let mut bare = short.clone();
            f.causal(&mut bare);
            bare.reverse();
            f.causal(&mut bare);
            bare.reverse();
            let mut padded = short.clone();
            f.apply(&mut padded);
            let err = |x: &[f64]| {
                x.iter()
                    .zip(&long)
                    .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()))
                    / peak
            };
            eprintln!(
                "end of record, pulse ending {label} before it: {:.1e} without zeros, {:.1e} with {}",
                err(&bare),
                err(&padded),
                f.pad
            );
        }
        // Phase: cells and the transpose.
        let gabor = Gabor::new(12.0, 30.0, 1e-3).unwrap();
        eprintln!(
            "phase cells: spacing {} samples, half-width {}, {} frequencies",
            gabor.spacing,
            gabor.half,
            gabor.frequencies.len()
        );
    }

    /// A trace at rest by its end comes out as if the record were longer, and one cut short
    /// by the record's end is still filtered symmetrically.
    #[test]
    fn filters_up_to_the_last_sample() {
        let dt = 1e-3;
        let f = Lowpass::butterworth(3.0, dt, 4).unwrap();
        // The pulse ends 0.15 s before the record does, half a period of the corner.
        let short = ricker(12.0, 0.4, dt, 600);
        let mut long = short.clone();
        long.resize(6000, 0.0);
        let mut a = short.clone();
        f.apply(&mut a);
        f.apply(&mut long);
        let peak = long.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        let err = a
            .iter()
            .zip(&long)
            .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()));
        assert!(err < 1e-6 * peak, "{err} of {peak}");
        // Symmetric still.
        let x: Vec<f64> = (0..300).map(|i| ((i * 7919) % 101) as f64 - 50.0).collect();
        let y: Vec<f64> = (0..300)
            .map(|i| ((i * 104729) % 97) as f64 - 48.0)
            .collect();
        let (mut fx, mut fy) = (x.clone(), y.clone());
        f.apply(&mut fx);
        f.apply(&mut fy);
        let l: f64 = fx.iter().zip(&y).map(|(a, b)| a * b).sum();
        let r: f64 = x.iter().zip(&fy).map(|(a, b)| a * b).sum();
        assert!((l - r).abs() <= 1e-12 * l.abs().max(r.abs()), "{l} {r}");
    }

    /// Near the Nyquist frequency the bilinear transform slows the poles' decay; the zeros
    /// follow the digital poles, so the end of the record is still exact there.
    #[test]
    fn pads_enough_near_the_nyquist_frequency() {
        let dt = 1e-3;
        for (corner, order) in [(380.0, 4), (300.0, 8), (5.0, 3)] {
            let f = Lowpass::butterworth(corner, dt, order).unwrap();
            let short: Vec<f64> = ricker(0.6 * corner, 0.3, dt, 400);
            let mut long = short.clone();
            long.resize(40 * short.len(), 0.0);
            let mut a = short.clone();
            f.apply(&mut a);
            f.apply(&mut long);
            let peak = long.iter().fold(0.0f64, |m, v| m.max(v.abs()));
            let err = a
                .iter()
                .zip(&long)
                .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()));
            assert!(err < 1e-6 * peak, "{corner} Hz: {err} of {peak}");
        }
    }

    /// The receivers' side weighs traces as the misfits' adjoint sources do: L2 by its
    /// weights, the time shifts and phases by each trace's observed slope or energy, so a
    /// pressure trace (a million times a velocity trace) weighs about 10⁻¹² of it, and a quiet
    /// trace nothing.
    #[test]
    fn receiver_weights_follow_the_adjoint_sources() {
        let dt = 2e-3;
        let n = 250;
        let trace = |a: f64| {
            ricker(12.0, 0.2, dt, n)
                .iter()
                .map(|v| a * v)
                .collect::<Vec<_>>()
        };
        let obs = Arc::new(vec![
            to_f32(&[trace(1e-3), trace(1e-9)].concat()),
            to_f32(&[trace(1e3), trace(2e3)].concat()),
        ]);
        let components = [Component::Vx, Component::P];
        for kind in [MisfitKind::L2, MisfitKind::Traveltime, MisfitKind::Phase] {
            let measure = Measure {
                kind,
                max_frequency: None,
                filter_order: 4,
                max_shift: None,
                balance: true,
                dominant_frequency: Some(12.0),
            };
            let w = Misfit::new(measure, &obs, &components, n, dt)
                .unwrap()
                .receiver_weights();
            match kind {
                MisfitKind::L2 => {
                    assert_eq!(w[0], vec![w[0][0]; 2]);
                    assert_eq!(w[1], vec![w[1][0]; 2]);
                }
                _ => {
                    // The quiet vx trace weighs nothing; twice the pressure, a quarter.
                    assert_eq!(w[0][1], 0.0);
                    assert!((w[1][1] / w[1][0] - 0.25).abs() < 1e-6, "{kind:?} {w:?}");
                }
            }
            let ratio = w[1][0] / w[0][0];
            assert!(ratio > 1e-13 && ratio < 1e-11, "{kind:?}: {ratio}");
        }
    }

    /// One-sample traces are fine for every misfit but the time shifts, which say so.
    #[test]
    fn short_traces_are_refused_only_where_they_must_be() {
        let obs = Arc::new(vec![vec![1.0f32, 2.0]]);
        for kind in [MisfitKind::L2, MisfitKind::Traveltime] {
            let measure = Measure {
                kind,
                max_frequency: None,
                filter_order: 4,
                max_shift: Some(0.01),
                balance: true,
                dominant_frequency: Some(10.0),
            };
            let m = Misfit::new(measure, &obs, &[Component::Vz], 1, 1e-3);
            assert_eq!(m.is_ok(), kind == MisfitKind::L2, "{kind:?}");
        }
    }

    /// What a misfit holds is what [`Misfit::held_bytes`] works out beforehand.
    #[test]
    fn counts_what_it_holds() {
        let samples = 300;
        let obs = Arc::new(vec![vec![0.5f32; 4 * samples], vec![1.0f32; 4 * samples]]);
        let components = [Component::Vx, Component::P];
        for kind in [MisfitKind::L2, MisfitKind::Traveltime, MisfitKind::Phase] {
            for max_frequency in [None, Some(20.0)] {
                let measure = Measure {
                    kind,
                    max_frequency,
                    filter_order: 4,
                    max_shift: None,
                    balance: true,
                    dominant_frequency: Some(12.0),
                };
                let m = Misfit::new(measure, &obs, &components, samples, 1e-3).unwrap();
                let planned = Misfit::held_bytes(&measure, 4, 2, samples, 1e-3).unwrap();
                assert_eq!(m.memory_bytes(), planned, "{kind:?} {max_frequency:?}");
            }
        }
    }

    #[test]
    fn refuses_corners_outside_the_band() {
        assert!(Lowpass::butterworth(0.0, 1e-3, 4).is_err());
        assert!(Lowpass::butterworth(500.0, 1e-3, 4).is_err());
        assert!(Lowpass::butterworth(420.0, 1e-3, 4).is_err());
        assert!(Lowpass::butterworth(10.0, 1e-3, 0).is_err());
    }
}
