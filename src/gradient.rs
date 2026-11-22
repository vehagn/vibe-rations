//! Gradients of a misfit with respect to the model, from the forward and adjoint wavefields
//! (Tarantola 1984; Tromp, Tape & Liu 2005; Fichtner et al. 2006; Fichtner 2011, ch. 9).
//!
//! The gradient is that of the discrete forward run, not of the wave equation: it comes from
//! transposing the solver's own step (see the `adjoint` module), so outside the absorbing
//! layers it agrees with finite differences of the solver to its round-off. With the forward
//! equations written as v⁽ⁿ⁺½⁾ - v⁽ⁿ⁻½⁾ = B (D σ⁽ⁿ⁾ + f⁽ⁿ⁾) and σ⁽ⁿ⁺¹⁾ - σ⁽ⁿ⁾ = M e⁽ⁿ⁺½⁾ + g⁽ⁿ⁾,
//! where e is the strain rate the stress update uses ([`Solver::strain_rate`]), their
//! Lagrange multipliers are λ̂v = W B⁻¹ ṽ and λ̂σ = -W ε̃, with W the cell weights (volume times
//! the share of the cell, half on the free surface), ṽ the adjoint velocities and ε̃ the
//! adjoint strain: the adjoint strain rates summed backwards in time, plus the pressure
//! residuals, so that the adjoint stresses are σ̃ = M ε̃. Then
//!
//! ∂χ/∂B = Σₙ W ṽ⁽ⁿ⁺½⁾ (v⁽ⁿ⁺½⁾ - v⁽ⁿ⁻½⁾) / B²  at each velocity position, and
//! ∂χ/∂M = -Σₙ W ε̃⁽ⁿ⁺¹⁾ ⊗ e⁽ⁿ⁺½⁾  at each stress position,
//!
//! the zero-lag correlations of the classic density and stiffness kernels. The chain rule
//! through the material averaging (arithmetic for density, harmonic for μ) and the
//! free surface's plane-stress moduli takes them to ρ, λ and μ at the nodes.
//!
//! The adjoint strain is summed rather than recovered from the adjoint stress as M⁻¹ σ̃: M is
//! singular in fluids, and inverting it would cancel digits where μ ≪ λ.
//!
//! The forward wavefield is needed backwards in time. Here it is recomputed from checkpoints
//! of the solver's state, which replays it exactly; reconstructing it from saved boundary
//! values instead (Dussaud et al. 2008) is cheaper and adds an error of its own.
//!
//! Inside the absorbing layers the strain rates leave out the CPML, and neither the layers
//! nor the sponge are self-adjoint, so the gradient there is only approximate. Just outside
//! them it is off by about 1e-6 with a free surface, whose sponge damps the velocities after
//! the stress update has used them, and by more if a velocity receiver sits on the first node
//! outside a layer ([`Simulation::receivers_next_to_layers`]).

use anyhow::{Result, bail, ensure};
use rayon::prelude::*;
use serde::Serialize;

use crate::adjoint::Adjoint;
use crate::grid::Grid;
use crate::model::Model;
use crate::receiver::{Component, Receivers};
use crate::simulation::Simulation;
use crate::solver::{Band, Field, Solver, State};

/// ½ Σ_c w_c Σ (d - d_obs)² over every recorded sample, and its derivative with respect to
/// each sample, w_c (d - d_obs): the adjoint source. Traces in the layout of
/// [`Receivers::traces`], one weight per component. Pressure is about 10⁶ times velocity, so
/// weight the components (e.g. with [`balancing_weights`]) when both are recorded.
pub fn l2_misfit(
    synthetic: &[Vec<f32>],
    observed: &[Vec<f32>],
    weights: &[f64],
) -> Result<(f64, Vec<Vec<f32>>)> {
    ensure!(
        synthetic.len() == observed.len() && synthetic.len() == weights.len(),
        "{} synthetic components, {} observed and {} weights",
        synthetic.len(),
        observed.len(),
        weights.len()
    );
    let mut misfit = 0.0;
    let mut residuals = Vec::with_capacity(synthetic.len());
    for ((s, o), &w) in synthetic.iter().zip(observed).zip(weights) {
        ensure!(
            s.len() == o.len(),
            "{} observed samples for {} synthetic ones",
            o.len(),
            s.len()
        );
        let r: Vec<f64> = s
            .iter()
            .zip(o)
            .map(|(&a, &b)| a as f64 - b as f64)
            .collect();
        misfit += 0.5 * w * r.iter().map(|x| x * x).sum::<f64>();
        residuals.push(r.iter().map(|&x| (w * x) as f32).collect());
    }
    Ok((misfit, residuals))
}

/// Weights 1/‖d_obs,c‖² that make each recorded component's misfit relative to its own
/// energy, so that pressure and velocity count alike. A velocity component that observed less
/// than 10⁻⁴ of the loudest velocity component's energy (vy on a symmetry plane, say) is
/// weighted as if it had that much, so its noise can't outweigh the others' signal; one that
/// observed nothing at all gets 0. `components` names the traces, as the receivers do.
pub fn balancing_weights(observed: &[Vec<f32>], components: &[Component]) -> Vec<f64> {
    let energy: Vec<f64> = observed
        .iter()
        .map(|o| o.iter().map(|&x| (x as f64).powi(2)).sum())
        .collect();
    let velocity = |c: &Component| !matches!(c, Component::P);
    let loudest = components
        .iter()
        .zip(&energy)
        .filter(|(c, _)| velocity(c))
        .fold(0.0f64, |m, (_, &e)| m.max(e));
    components
        .iter()
        .zip(&energy)
        .map(|(c, &e)| {
            let floor = if velocity(c) { 1e-4 * loudest } else { 0.0 };
            if e > 0.0 { 1.0 / e.max(floor) } else { 0.0 }
        })
        .collect()
}

/// The derivative of a misfit with respect to the model at every interior node (x fastest,
/// like [`Model`]), each with the other two parameters of its set held fixed.
#[derive(Clone, Debug)]
pub struct Gradient {
    /// ∂χ/∂ρ with λ and μ fixed, in units of χ per kg/m³.
    pub rho: Vec<f64>,
    /// ∂χ/∂λ, per Pa.
    pub lambda: Vec<f64>,
    /// ∂χ/∂μ, per Pa. Where μ = 0 the shear moduli's harmonic mean has a kink (it is 0 for
    /// any μ ≤ 0 and rises as 4μ above), so at a fluid node next to solid ones this holds
    /// only the part through the normal stresses; ∂χ/∂vs, where μ = ρ vs² smooths the kink
    /// out, is exact.
    pub mu: Vec<f64>,
}

impl Gradient {
    /// The same gradient for the parameters the model is given in: ∂χ/∂vp and ∂χ/∂vs with
    /// the other speed and ρ fixed, and ∂χ/∂ρ with both speeds fixed. With λ = ρ(vp² - 2vs²)
    /// and μ = ρ vs²: ∂χ/∂vp = 2ρ vp ∂χ/∂λ, ∂χ/∂vs = 2ρ vs (∂χ/∂μ - 2 ∂χ/∂λ) and
    /// ∂χ/∂ρ|vp,vs = ∂χ/∂ρ|λ,μ + (vp² - 2vs²) ∂χ/∂λ + vs² ∂χ/∂μ.
    pub fn velocities(&self, model: &Model) -> [Vec<f64>; 3] {
        let n = self.rho.len();
        let (mut vp, mut vs, mut rho) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
        for i in 0..n {
            let (a, b, r) = (model.vp[i] as f64, model.vs[i] as f64, model.rho[i] as f64);
            let (gl, gm) = (self.lambda[i], self.mu[i]);
            vp[i] = 2.0 * r * a * gl;
            vs[i] = 2.0 * r * b * (gm - 2.0 * gl);
            rho[i] = self.rho[i] + (a * a - 2.0 * b * b) * gl + b * b * gm;
        }
        [vp, vs, rho]
    }
}

/// What a gradient run found.
pub struct Run {
    pub gradient: Gradient,
    /// The misfit the adjoint source came from, as the closure returned it.
    pub misfit: f64,
    /// Steps between checkpoints (with [`Strategy::Reconstruct`], the segments' length).
    pub checkpoint_every: usize,
    /// The most memory the run held at once, counted from its arrays: the forward run, the
    /// adjoint run with its residuals, the checkpoints, one replayed segment's velocities, the
    /// correlation sums, the illumination and the strains.
    /// [`crate::memory::Plan::gradient_memory`] works it out
    /// beforehand (with the observed traces, which are the caller's).
    pub memory_bytes: usize,
    /// Σₙ |v⁽ⁿ⁺½⁾|² of the forward run at every interior node (x fastest): how much of the
    /// source's wavefield reaches it, the diagonal of a source-side pseudo-Hessian (Shin et
    /// al. 2001) to precondition an inversion with.
    pub illumination: Vec<f64>,
    /// How the forward wavefield was had back, with the interval or segment it used.
    pub strategy: Strategy,
}

/// How a gradient run gets the forward wavefield back, step by step in reverse, to correlate
/// with the adjoint's.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Strategy {
    /// Replay it a segment at a time from checkpoints `every` steps apart (by default
    /// [`checkpoint_interval`]'s), keeping a segment's velocities: exact, bit for bit.
    Checkpoints { every: Option<usize> },
    /// Rebuild it backwards ([`Solver::step_back_with`]) from the band around the absorbing
    /// layers ([`Band`]), saved for segments of `segment` steps: the last on the way forward,
    /// the others while replaying them from a checkpoint at their start (by default the length
    /// that needs least memory, [`reconstruction_segment`]). A segment as long as the run
    /// needs neither checkpoints nor replays. Exact outside the layers up to single-precision
    /// round-off, which accumulates over a segment; inside them it means nothing.
    Reconstruct { segment: Option<usize> },
}

impl Default for Strategy {
    fn default() -> Self {
        Self::Checkpoints { every: None }
    }
}

/// Steps per segment when the forward wavefield is rebuilt: `requested`, or by default the
/// length k that needs least memory for ⌈N/k⌉ - 1 checkpoints of `checkpoint_bytes` each
/// and k band records of `record_bytes`, about √(N C/B). Saving the band for the whole run
/// (k = N) spares the replays, but needs more memory unless the band is small next to a
/// checkpoint, as on large grids.
pub fn reconstruction_segment(
    steps: usize,
    checkpoint_bytes: usize,
    record_bytes: usize,
    requested: Option<usize>,
) -> usize {
    let steps = steps.max(1);
    match requested {
        Some(k) => k.clamp(1, steps),
        None => (1..=steps)
            .min_by_key(|&k| (steps.div_ceil(k) - 1) * checkpoint_bytes + k * record_bytes)
            .expect("at least one step"),
    }
}

/// Steps between checkpoints: `requested`, or by default ⌈√(3N)⌉, which about minimises the
/// memory of ⌈N/k⌉ checkpoints of nine arrays and a replayed segment of k + 1 velocity triples.
pub fn checkpoint_interval(steps: usize, requested: Option<usize>) -> usize {
    requested
        .unwrap_or_else(|| ((3 * steps) as f64).sqrt().ceil() as usize)
        .clamp(1, steps.max(1))
}

/// The correlations summed over the run, on the padded grid, per position: the density
/// terms Σ ṽ Δv for vx, vy, vz, then the stiffness terms for l2m and lam (on the free surface
/// the plane-stress a and b instead) and for mxy, mxz, myz, already weighted by -W.
type Sums = Vec<[f64; 8]>;

/// Runs `sim` (fresh, at step 0) forward, recording its receivers and keeping a checkpoint
/// every `checkpoint_every` steps (by default [`checkpoint_interval`]'s); hands the
/// receivers to `adjoint_source`, which returns the misfit and its derivative with respect to
/// every sample (e.g. [`l2_misfit`]); then runs the adjoint backwards, replaying the forward
/// run a segment at a time from the checkpoints, and correlates the two. The receivers'
/// traces are as `run` leaves them: un-warped when the run removes time dispersion.
///
/// The time step and the absorbing layers stay as `sim` set them up, so the gradient is that
/// of the forward run with them fixed. With dt from a Courant fraction, both depend on the
/// fastest vp in the model; perturbing it changes them too.
///
/// Afterwards `sim` holds the recorded traces, but its solver is wherever the last replay
/// left it (the end of the first segment), not at the end of the run.
pub fn gradient(
    sim: &mut Simulation,
    adjoint_source: impl FnOnce(&Receivers) -> Result<(f64, Vec<Vec<f32>>)>,
    checkpoint_every: Option<usize>,
) -> Result<Run> {
    gradient_with(
        sim,
        adjoint_source,
        Strategy::Checkpoints {
            every: checkpoint_every,
        },
    )
}

/// [`gradient`], getting the forward wavefield back as `strategy` says. With
/// [`Strategy::Reconstruct`] the gradient outside the layers is the replay's up to round-off
/// (the band reaches the positions just outside the region that feed its first nodes), and
/// inside them it means nothing; afterwards the solver holds the rebuilt state before the
/// first step, and its CPML memory is stale.
pub fn gradient_with(
    sim: &mut Simulation,
    adjoint_source: impl FnOnce(&Receivers) -> Result<(f64, Vec<Vec<f32>>)>,
    strategy: Strategy,
) -> Result<Run> {
    ensure!(sim.solver.step == 0, "the forward run must start at rest");
    let steps = sim.steps;
    let band = matches!(strategy, Strategy::Reconstruct { .. }).then(|| sim.band());
    let (every, strategy) = match strategy {
        Strategy::Checkpoints { every } => {
            let every = checkpoint_interval(steps, every);
            (every, Strategy::Checkpoints { every: Some(every) })
        }
        Strategy::Reconstruct { segment } => {
            let record = 4 * band.as_ref().map_or(0, Band::len);
            let k = reconstruction_segment(steps, sim.solver.state_bytes(), record, segment);
            (k, Strategy::Reconstruct { segment: Some(k) })
        }
    };
    let rebuild = matches!(strategy, Strategy::Reconstruct { .. });
    let segments = steps.div_ceil(every);
    // The last segment rebuilt needs no checkpoint: its band is saved on the way forward.
    let last = if rebuild {
        segments.saturating_sub(1) * every
    } else {
        steps
    };
    let mut records = band
        .as_ref()
        .map_or(Vec::new(), |b| vec![0.0f32; every * b.len()]);

    // Forward, keeping the state at the start of every segment.
    let mut checkpoints: Vec<State> = Vec::with_capacity(segments);
    for n in 0..steps {
        if n.is_multiple_of(every) && n < last {
            checkpoints.push(sim.solver.state());
        }
        if let Some(b) = &band
            && n >= last
        {
            let r = (n - last) * b.len();
            b.save(&sim.solver.fields, &mut records[r..r + b.len()]);
        }
        sim.advance();
        if (n.is_multiple_of(100) || n + 1 == steps) && !sim.solver.peak_velocity().is_finite() {
            bail!("the forward wavefield blew up by step {}", n + 1);
        }
    }
    // The traces as `run` writes them, and as observed data from it are: un-warped when the
    // run removes time dispersion (crate::tdt). The misfit is of those, so its derivative goes
    // back through the un-warping's transpose to the samples the solver recorded, which is
    // what the adjoint run injects. The gradient is then that of the whole discrete run,
    // the transforms included.
    sim.finish();
    let (misfit, residuals) = adjoint_source(&sim.receivers)?;
    ensure!(
        misfit.is_finite() && residuals.iter().flatten().all(|r| r.is_finite()),
        "the misfit or its derivative isn't finite; check the observed data"
    );
    let residuals = sim.recorded_residuals(residuals);
    let mut c = Correlator::new(sim, residuals, steps)?;
    let len = sim.solver.grid.padded_len();
    // Besides the checkpoints and what holds the forward field: both runs, the sums, the
    // illumination and the two strains.
    let fixed = sim.memory_bytes() + c.memory_bytes();
    let held_checkpoints = checkpoints.iter().map(State::memory_bytes).sum::<usize>();
    let mut peak_bytes = 0;

    if let Some(band) = &band {
        let blen = band.len();
        // The velocities after the step being taken back, for the density term.
        let mut after: [Vec<f32>; 3] = std::array::from_fn(|_| vec![0.0f32; len]);
        peak_bytes = fixed
            + held_checkpoints
            + 4 * records.len()
            + band.memory_bytes()
            + 4 * after.iter().map(Vec::len).sum::<usize>();
        for s in (0..segments).rev() {
            let start = s * every;
            let end = (start + every).min(steps);
            if start < last {
                // Replay the segment for its band; it ends where the next one started.
                sim.solver.restore(&checkpoints[s]);
                for n in start..end {
                    let r = (n - start) * blen;
                    band.save(&sim.solver.fields, &mut records[r..r + blen]);
                    sim.step();
                }
            }
            for n in (start..end).rev() {
                let f = &sim.solver.fields;
                for (a, v) in after.iter_mut().zip([&f.vx, &f.vy, &f.vz]) {
                    a.copy_from_slice(v);
                }
                let r = (n - start) * blen;
                sim.step_back(band, &records[r..r + blen]);
                let f = &sim.solver.fields;
                c.step(
                    &sim.solver,
                    n,
                    [&f.vx, &f.vy, &f.vz],
                    [&after[0], &after[1], &after[2]],
                );
            }
            c.check()?;
        }
    } else {
        let mut segment: Vec<[Vec<f32>; 3]> = Vec::with_capacity(every + 1);
        for (s, checkpoint) in checkpoints.iter().enumerate().rev() {
            let start = s * every;
            let end = (start + every).min(steps);
            // The velocities before the segment's first step, then after each of its steps.
            sim.solver.restore(checkpoint);
            segment.clear();
            segment.push(velocities(&sim.solver));
            for _ in start..end {
                sim.step();
                segment.push(velocities(&sim.solver));
            }
            let held = fixed + held_checkpoints + 4 * 3 * len * segment.len();
            peak_bytes = peak_bytes.max(held);

            for n in (start..end).rev() {
                let (before, after) = (&segment[n - start], &segment[n - start + 1]);
                c.step(
                    &sim.solver,
                    n,
                    [&before[0], &before[1], &before[2]],
                    [&after[0], &after[1], &after[2]],
                );
            }
            c.check()?;
        }
    }

    // The gradient's own arrays mustn't come on top of the rest.
    drop((checkpoints, records));
    let (sums, illumination) = c.finish();
    Ok(Run {
        gradient: to_nodes(&sim.solver, &sums),
        misfit,
        checkpoint_every: every,
        memory_bytes: peak_bytes,
        illumination,
        strategy,
    })
}

/// The adjoint run and what it is correlated into, a forward step at a time in reverse.
struct Correlator {
    adjoint: Adjoint,
    pressure: PressureResiduals,
    sums: Sums,
    illumination: Vec<f64>,
    strain: [Vec<f32>; 6],
    adjoint_strain: [Vec<f32>; 6],
}

impl Correlator {
    fn new(sim: &Simulation, residuals: Vec<Vec<f32>>, steps: usize) -> Result<Self> {
        let pressure = PressureResiduals::new(sim, &residuals);
        let adjoint = Adjoint::new(sim, residuals)?;
        let g = sim.solver.grid;
        let len = g.padded_len();
        let mut adjoint_strain = std::array::from_fn::<_, 6, _>(|_| vec![0.0f32; len]);
        pressure.add(&mut adjoint_strain, steps);
        Ok(Self {
            adjoint,
            pressure,
            sums: vec![[0.0; 8]; len],
            illumination: vec![0.0f64; g.interior_len()],
            strain: std::array::from_fn::<_, 6, _>(|_| vec![0.0f32; len]),
            adjoint_strain,
        })
    }

    fn memory_bytes(&self) -> usize {
        self.adjoint.memory_bytes()
            + self.pressure.memory_bytes()
            + std::mem::size_of_val(self.sums.as_slice())
            + std::mem::size_of_val(self.illumination.as_slice())
            + 4 * self
                .strain
                .iter()
                .chain(&self.adjoint_strain)
                .map(Vec::len)
                .sum::<usize>()
    }

    /// Correlates forward step `n`, whose velocities were `before` and `after` it, with the
    /// adjoint step that pairs with it, which it takes.
    fn step(&mut self, forward: &Solver, n: usize, before: [&[f32]; 3], after: [&[f32]; 3]) {
        let g = forward.grid;
        // λ̂σ⁽ⁿ⁺¹⁾ = -W ε̃⁽ⁿ⁺¹⁾ with the forward strain rate at n + ½.
        forward.strain_rate(after, &mut self.strain);
        correlate_stiffness(forward, &mut self.sums, &self.adjoint_strain, &self.strain);
        // ṽ⁽ⁿ⁺½⁾ with the forward velocity change over step n.
        self.adjoint.step();
        assert_eq!(
            self.adjoint.forward_step(),
            Some(n),
            "adjoint and forward out of step"
        );
        correlate_density(&self.adjoint.solver, &mut self.sums, before, after);
        illuminate(&g, &mut self.illumination, after);
        // ε̃⁽ⁿ⁾ = ε̃⁽ⁿ⁺¹⁾ + e(ṽ⁽ⁿ⁺½⁾) - W⁻¹ r_σ⁽ⁿ⁾.
        let f = &self.adjoint.solver.fields;
        self.adjoint
            .solver
            .strain_rate([&f.vx, &f.vy, &f.vz], &mut self.strain);
        self.adjoint_strain
            .par_iter_mut()
            .zip(&self.strain)
            .for_each(|(a, e)| a.iter_mut().zip(e).for_each(|(a, e)| *a += e));
        self.pressure.add(&mut self.adjoint_strain, n);
    }

    fn check(&self) -> Result<()> {
        ensure!(
            self.adjoint.solver.peak_velocity().is_finite(),
            "the adjoint wavefield blew up"
        );
        Ok(())
    }

    /// The sums and the illumination, freeing the rest.
    fn finish(self) -> (Sums, Vec<f64>) {
        (self.sums, self.illumination)
    }
}

fn velocities(s: &Solver) -> [Vec<f32>; 3] {
    [
        s.fields.vx.clone(),
        s.fields.vy.clone(),
        s.fields.vz.clone(),
    ]
}

/// The pressure residuals, as they enter the adjoint strain.
struct PressureResiduals {
    /// For each receiver, the terms it reads pressure from: (normal component 0–2, padded
    /// index, its interpolation weight / (3W) there). On a node, one per component.
    terms: Vec<Vec<(usize, usize, f32)>>,
    /// The pressure residuals, (receivers, samples), if pressure was recorded.
    values: Option<Vec<f32>>,
    samples: usize,
    record_every: usize,
}

impl PressureResiduals {
    fn new(sim: &Simulation, residuals: &[Vec<f32>]) -> Self {
        let r = &sim.receivers;
        let s = &sim.solver;
        let g = s.grid;
        let volume = g.cell_volume() as f32;
        let plane = g.plane_len();
        let terms = r
            .points
            .iter()
            .map(|&p| {
                [Field::Sxx, Field::Syy, Field::Szz]
                    .into_iter()
                    .enumerate()
                    .flat_map(|(c, f)| {
                        s.field_weights(f, p).into_iter().map(move |(i, w)| {
                            let k = i / plane - g.halo;
                            (
                                c,
                                i,
                                w * (1.0 / (3.0 * volume * s.cell_share(Field::Sxx, k))),
                            )
                        })
                    })
                    .collect()
            })
            .collect();
        let values = r
            .components
            .iter()
            .position(|&c| c == Component::P)
            .and_then(|c| residuals.get(c).cloned());
        Self {
            terms,
            values,
            samples: r.samples,
            record_every: sim.record_every,
        }
    }

    fn memory_bytes(&self) -> usize {
        4 * self.values.as_ref().map_or(0, Vec::len)
    }

    /// Adds -W⁻¹ r_σ⁽ⁿ⁾, the residual of the pressure read from σ⁽ⁿ⁾ (recorded after forward
    /// step n - 1), to the normal components of the adjoint strain.
    fn add(&self, adjoint_strain: &mut [Vec<f32>; 6], n: usize) {
        let Some(values) = &self.values else {
            return;
        };
        let Some(step) = n.checked_sub(1) else {
            return;
        };
        if !step.is_multiple_of(self.record_every) {
            return;
        }
        let sample = step / self.record_every;
        if sample >= self.samples {
            return;
        }
        for (r, terms) in self.terms.iter().enumerate() {
            let delta = values[r * self.samples + sample];
            for &(c, p, w) in terms {
                adjoint_strain[c][p] += w * delta;
            }
        }
    }
}

/// The cell weight W of `field` on interior plane `k`.
fn weight(s: &Solver, field: Field, k: usize) -> f64 {
    s.grid.cell_volume() * s.cell_share(field, k) as f64
}

/// Adds -W ε̃ ⊗ e for the stiffness at every stress position.
fn correlate_stiffness(
    s: &Solver,
    sums: &mut Sums,
    adjoint_strain: &[Vec<f32>; 6],
    strain: &[Vec<f32>; 6],
) {
    let g = s.grid;
    let [nx, ny, nz] = g.n;
    let len = g.plane_len();
    let halo = g.halo;
    let surface = s.free_surface;
    sums[halo * len..(halo + nz) * len]
        .par_chunks_mut(len)
        .enumerate()
        .for_each(|(k, plane)| {
            let w_normal = weight(s, Field::Sxx, k);
            let w_xy = weight(s, Field::Sxy, k);
            let w_xz = weight(s, Field::Sxz, k);
            let w_yz = weight(s, Field::Syz, k);
            let on_surface = surface && k == 0;
            for j in 0..ny {
                let row = g.idx(0, j, k);
                let local = row - (halo + k) * len;
                for i in 0..nx {
                    let p = row + i;
                    let a = |c: usize| adjoint_strain[c][p] as f64;
                    let e = |c: usize| strain[c][p] as f64;
                    let (xx, yy, zz) = (-w_normal * a(0), -w_normal * a(1), -w_normal * a(2));
                    let (exx, eyy, ezz) = (e(0), e(1), e(2));
                    let sum = &mut plane[local + i];
                    if on_surface {
                        // σxx += a ∂x vx + b ∂y vy, σyy += b ∂x vx + a ∂y vy.
                        sum[3] += xx * exx + yy * eyy;
                        sum[4] += xx * eyy + yy * exx;
                    } else {
                        sum[3] += xx * exx + yy * eyy + zz * ezz;
                        sum[4] += xx * (eyy + ezz) + yy * (exx + ezz) + zz * (exx + eyy);
                    }
                    sum[5] += -w_xy * a(3) * e(3);
                    sum[6] += -w_xz * a(4) * e(4);
                    sum[7] += -w_yz * a(5) * e(5);
                }
            }
        });
}

/// Adds ṽ (v⁽ⁿ⁺½⁾ - v⁽ⁿ⁻½⁾), with ṽ the adjoint velocities `adjoint` holds now, for every
/// velocity position. The weight and buoyancy come in when the sums go to the nodes.
fn correlate_density(adjoint: &Solver, sums: &mut Sums, before: [&[f32]; 3], after: [&[f32]; 3]) {
    let g = adjoint.grid;
    let [nx, ny, nz] = g.n;
    let len = g.plane_len();
    let halo = g.halo;
    let f = &adjoint.fields;
    let back = [&f.vx, &f.vy, &f.vz];
    sums[halo * len..(halo + nz) * len]
        .par_chunks_mut(len)
        .enumerate()
        .for_each(|(k, plane)| {
            for j in 0..ny {
                let row = g.idx(0, j, k);
                let local = row - (halo + k) * len;
                for i in 0..nx {
                    let p = row + i;
                    let sum = &mut plane[local + i];
                    for c in 0..3 {
                        sum[c] += back[c][p] as f64 * (after[c][p] - before[c][p]) as f64;
                    }
                }
            }
        });
}

/// Adds |v|² at every interior node, each component taken where it is stored.
fn illuminate(g: &Grid, illumination: &mut [f64], v: [&[f32]; 3]) {
    let [nx, ny, _] = g.n;
    illumination
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

/// Takes the sums from the staggered positions to ρ, λ and μ at the nodes.
fn to_nodes(s: &Solver, sums: &Sums) -> Gradient {
    let g = s.grid;
    let [nx, ny, nz] = g.n;
    let dt = s.dt;
    let model = &s.model;
    let n = g.interior_len();
    let mut out = Gradient {
        rho: vec![0.0; n],
        lambda: vec![0.0; n],
        mu: vec![0.0; n],
    };
    let node =
        |i: usize, j: usize, k: usize| g.interior_idx(i.min(nx - 1), j.min(ny - 1), k.min(nz - 1));
    let mu_at = |q: usize| model.rho[q] as f64 * (model.vs[q] as f64).powi(2);
    let m = &s.material;
    let fields = [Field::Vx, Field::Vy, Field::Vz];
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let p = g.idx(i, j, k);
                let here = node(i, j, k);
                let sum = &sums[p];
                // Density: B = 2dt/(ρ₁ + ρ₂), ∂B/∂ρ = -B²/(2dt), so ∂χ/∂ρ = -W Σ ṽ Δv / (2dt)
                // for each of the two nodes.
                for (c, &field) in fields.iter().enumerate() {
                    let d = -weight(s, field, k) * sum[c] / (2.0 * dt);
                    let mut next = [i, j, k];
                    next[c] += 1;
                    out.rho[here] += d;
                    out.rho[node(next[0], next[1], next[2])] += d;
                }
                // l2m = dt(λ + 2μ), lam = dt λ; on the free surface through a and b.
                let (g_l2m, g_lam) = if s.free_surface && k == 0 {
                    let (l2, l) = (m.l2m[p] as f64, m.lam[p] as f64);
                    let r = if l2 > 0.0 { l / l2 } else { 0.0 };
                    let (ga, gb) = (sum[3], sum[4]);
                    (
                        ga * (1.0 + r * r) + gb * r * r,
                        ga * (-2.0 * r) + gb * (1.0 - 2.0 * r),
                    )
                } else {
                    (sum[3], sum[4])
                };
                out.lambda[here] += dt * (g_l2m + g_lam);
                out.mu[here] += 2.0 * dt * g_l2m;
                // Shear: dt × the harmonic mean H of μ at four nodes, ∂H/∂μ_q = H²/(4 μ_q²).
                for (c, axes) in [(5, [0, 1]), (6, [0, 2]), (7, [1, 2])] {
                    let corners = [[0, 0], [1, 0], [0, 1], [1, 1]].map(|[a, b]| {
                        let mut q = [i, j, k];
                        q[axes[0]] += a;
                        q[axes[1]] += b;
                        node(q[0], q[1], q[2])
                    });
                    let mus = corners.map(mu_at);
                    if mus.iter().any(|&x| x <= 0.0) {
                        continue;
                    }
                    let h = 4.0 / mus.iter().map(|x| 1.0 / x).sum::<f64>();
                    for (&q, &mu) in corners.iter().zip(&mus) {
                        out.mu[q] += dt * sum[c] * h * h / (4.0 * mu * mu);
                    }
                }
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::Grid;

    #[test]
    fn the_l2_misfit_and_its_residuals() {
        let synthetic = vec![vec![1.0, 2.0], vec![0.5, 0.0]];
        let observed = vec![vec![1.0, 1.0], vec![0.0, 0.0]];
        let (m, r) = l2_misfit(&synthetic, &observed, &[1.0, 1.0]).unwrap();
        assert!((m - 0.5 * (1.0 + 0.25)).abs() < 1e-12);
        assert_eq!(r, vec![vec![0.0, 1.0], vec![0.5, 0.0]]);
        let (m, r) = l2_misfit(&synthetic, &observed, &[2.0, 0.0]).unwrap();
        assert!((m - 1.0).abs() < 1e-12);
        assert_eq!(r, vec![vec![0.0, 2.0], vec![0.0, 0.0]]);
        assert!(l2_misfit(&synthetic, &observed[..1], &[1.0]).is_err());
        assert!(l2_misfit(&synthetic, &observed, &[1.0]).is_err());
        assert!(l2_misfit(&synthetic, &[vec![1.0], vec![0.0, 0.0]], &[1.0, 1.0]).is_err());
    }

    #[test]
    fn balancing_weights_divide_by_the_observed_energy() {
        use Component::{P, Vx, Vy, Vz};
        let w = balancing_weights(&[vec![3.0, 4.0], vec![0.0; 2]], &[Vx, Vz]);
        assert_eq!(w, vec![1.0 / 25.0, 0.0]);
        // A nearly silent velocity component is weighed as if it had 10⁻⁴ of the loudest's
        // energy; pressure, in other units, isn't compared with velocity.
        let w = balancing_weights(&[vec![100.0], vec![0.1], vec![0.5]], &[Vx, Vy, P]);
        assert_eq!(w, vec![1e-4, 1.0, 4.0]);
        let _ = Vz;
    }

    /// The replay is exact, so how often the run keeps checkpoints doesn't change the
    /// gradient by a bit: every step, intervals that don't divide the run, the default, and
    /// one longer than the run.
    #[test]
    fn checkpoint_intervals_give_identical_gradients() {
        let (steps, reference) = small_gradient(Strategy::Checkpoints { every: Some(1) });
        assert!(steps > 20 && steps % 9 != 0, "{steps} steps");
        assert!(reference.misfit > 0.0);
        assert!(reference.gradient.mu.iter().any(|&g| g != 0.0));
        for every in [Some(9), None, Some(steps), Some(1 << 40)] {
            let (_, r) = small_gradient(Strategy::Checkpoints { every });
            assert!(r.checkpoint_every <= steps, "{every:?}");
            assert_eq!(r.misfit, reference.misfit, "{every:?}");
            assert_eq!(r.gradient.rho, reference.gradient.rho, "{every:?}");
            assert_eq!(r.gradient.lambda, reference.gradient.lambda, "{every:?}");
            assert_eq!(r.gradient.mu, reference.gradient.mu, "{every:?}");
        }
    }

    /// Rebuilt backwards from the band, the forward wavefield gives the replay's gradient at
    /// every node outside the layers to round-off, however long the segments: with one step
    /// each, most of them replayed; seven; the default; and the whole run, none replayed.
    #[test]
    fn rebuilt_gradients_match_the_replay() {
        let (steps, reference) = small_gradient(Strategy::Checkpoints { every: None });
        let sim = Simulation::new(&small_config(), std::path::Path::new(".")).unwrap();
        let [x, y, z] = sim.solver.interior();
        let g = sim.solver.grid;
        let nodes: Vec<usize> = z
            .flat_map(|k| {
                let x = x.clone();
                y.clone()
                    .flat_map(move |j| x.clone().map(move |i| g.interior_idx(i, j, k)))
            })
            .collect();
        for segment in [Some(1), Some(7), None, Some(steps)] {
            let (_, r) = small_gradient(Strategy::Reconstruct { segment });
            assert_eq!(r.misfit, reference.misfit, "{segment:?}");
            assert!(matches!(
                r.strategy,
                Strategy::Reconstruct { segment: Some(_) }
            ));
            for (name, a, b) in [
                ("rho", &r.gradient.rho, &reference.gradient.rho),
                ("lambda", &r.gradient.lambda, &reference.gradient.lambda),
                ("mu", &r.gradient.mu, &reference.gradient.mu),
            ] {
                let scale = nodes.iter().fold(0.0f64, |m, &n| m.max(b[n].abs()));
                let worst = nodes
                    .iter()
                    .fold(0.0f64, |m, &n| m.max((a[n] - b[n]).abs() / scale));
                println!("rebuilt, segment {segment:?}, {name}: off by {worst:.1e} of the largest");
                assert!(scale > 0.0 && worst < 1e-4, "{segment:?}, {name}: {worst}");
            }
        }
    }

    fn small_config() -> crate::config::Config {
        crate::config::Config::from_toml(
            r#"
            [grid]
            points = [20, 18, 16]
            spacing = 10.0
            [time]
            duration = 0.03
            [boundary]
            pml_width = 4
            [[model.layer]]
            top = 0.0
            vp = 2000.0
            vs = 1000.0
            rho = 2000.0
            [[model.layer]]
            top = 70.0
            vp = 2500.0
            vs = 1400.0
            rho = 2200.0
            [[source]]
            position = [100.0, 90.0, 50.0]
            force = [0.0, 0.0, 1e9]
            wavelet = { kind = "ricker", frequency = 40.0 }
            [receivers]
            components = ["vx", "vz", "p"]
            [[receivers.point]]
            position = [70.0, 90.0, 0.0]
            [[receivers.point]]
            position = [120.0, 80.0, 100.0]
            [output]
            record_every = 3
            "#,
        )
        .unwrap()
    }

    /// A gradient of a small run with two layers, a force and three components, against
    /// silent observed traces.
    fn small_gradient(strategy: Strategy) -> (usize, Run) {
        let mut sim = Simulation::new(&small_config(), std::path::Path::new(".")).unwrap();
        let steps = sim.steps;
        let zeros: Vec<Vec<f32>> = sim
            .receivers
            .traces
            .iter()
            .map(|t| vec![0.0; t.len()])
            .collect();
        let weights = [1.0, 1.0, 1e-18];
        let r = gradient_with(
            &mut sim,
            |r| l2_misfit(&r.traces, &zeros, &weights),
            strategy,
        )
        .unwrap();
        (steps, r)
    }

    #[test]
    fn refuses_a_misfit_that_is_not_finite() {
        use crate::config::Config;
        let cfg = Config::from_toml(
            r#"
            [grid]
            points = [16, 16, 16]
            spacing = 10.0
            [time]
            duration = 0.01
            [boundary]
            pml_width = 4
            [[model.layer]]
            top = 0.0
            vp = 2000.0
            vs = 1000.0
            rho = 2000.0
            [[source]]
            position = [80.0, 80.0, 50.0]
            explosion = 1e12
            wavelet = { kind = "ricker", frequency = 40.0 }
            [[receivers.point]]
            position = [70.0, 80.0, 0.0]
            "#,
        )
        .unwrap();
        let mut sim = Simulation::new(&cfg, std::path::Path::new(".")).unwrap();
        let e = gradient(
            &mut sim,
            |r| {
                let mut nan = r.traces.clone();
                nan[0][0] = f32::NAN;
                l2_misfit(&r.traces, &nan, &[1.0; 3])
            },
            None,
        )
        .err()
        .unwrap();
        assert!(e.to_string().contains("isn't finite"), "{e}");
    }

    /// The velocity parametrisation is the chain rule: a linear functional of (ρ, λ, μ)
    /// differentiated through λ = ρ(vp² - 2vs²), μ = ρ vs².
    #[test]
    fn converts_to_the_velocity_parametrisation() {
        let g = Grid::new([1, 1, 1], [1.0; 3]);
        let model = Model::homogeneous(g, 3000.0, 1500.0, 2200.0);
        let grad = Gradient {
            rho: vec![0.7],
            lambda: vec![-2e-10],
            mu: vec![5e-10],
        };
        let f = |vp: f64, vs: f64, rho: f64| {
            0.7 * rho - 2e-10 * rho * (vp * vp - 2.0 * vs * vs) + 5e-10 * rho * vs * vs
        };
        let [gvp, gvs, grho] = grad.velocities(&model);
        let (vp, vs, rho) = (3000.0, 1500.0, 2200.0);
        let d = |a: f64, b: f64, step: f64| (a - b) / (2.0 * step);
        let close = |a: f64, b: f64| (a - b).abs() < 1e-6 * b.abs();
        assert!(close(
            gvp[0],
            d(f(vp + 1.0, vs, rho), f(vp - 1.0, vs, rho), 1.0)
        ));
        assert!(close(
            gvs[0],
            d(f(vp, vs + 1.0, rho), f(vp, vs - 1.0, rho), 1.0)
        ));
        assert!(close(
            grho[0],
            d(f(vp, vs, rho + 1.0), f(vp, vs, rho - 1.0), 1.0)
        ));
    }
}
