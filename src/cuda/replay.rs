//! A gradient run on a device: [`crate::gradient::gradient`]'s checkpointed replay, with
//! every step, strain rate, correlation and sum on the device, in the CPU's order and
//! arithmetic, so the gradient, the illumination and the traces are the CPU's bit for bit.
//!
//! What crosses: the checkpoints, down on the way forward and up again for each segment's
//! replay (the host keeps them, since the device has less memory); the traces once, for the
//! misfit; the residuals once, as the adjoint's sources; and the sums at the end. A segment's
//! velocities, both runs and the sums stay on the device.

use anyhow::{Result, bail, ensure};

use super::{Args, Device, GpuSolver, Kernel, LINE, Shared, bytes, bytes_mut};
use crate::adjoint::Adjoint;
use crate::gradient::{PressureResiduals, Run, Strategy, Sums, checkpoint_interval, to_nodes};
use crate::receiver::Receivers;
use crate::simulation::Simulation;
use crate::solver::PointOp;

/// [`crate::gradient::gradient`] on `device`: runs `sim` (fresh, at step 0) forward with a
/// checkpoint every `checkpoint_every` steps (by default [`checkpoint_interval`]'s), hands the
/// receivers to `adjoint_source`, then runs the adjoint back, replaying the forward run a
/// segment at a time, and correlates the two. The same results bit for bit, whatever either
/// interval, unless the device fuses products into sums ([`Device::name`] says).
///
/// Afterwards `sim` holds the recorded traces and the forward run's last state.
pub fn gradient<D: Device>(
    sim: &mut Simulation,
    adjoint_source: impl FnOnce(&Receivers) -> Result<(f64, Vec<Vec<f32>>)>,
    checkpoint_every: Option<usize>,
    device: D,
) -> Result<Run> {
    ensure!(sim.solver.step == 0, "the forward run must start at rest");
    let steps = sim.steps;
    let every = checkpoint_interval(steps, checkpoint_every);
    let device = Shared::new(device);
    let g = sim.solver.grid;
    let len = g.padded_len();

    // Forward, keeping the state at the start of every segment.
    let mut forward = GpuSolver::new(device.clone(), &sim.solver)?;
    forward.attach(&sim.sources, steps, &sim.receivers)?;
    let mut checkpoints = Vec::with_capacity(steps.div_ceil(every));
    for n in 0..steps {
        if n.is_multiple_of(every) {
            let mut state = sim.solver.blank_state();
            forward.download_state(&mut state)?;
            checkpoints.push(state);
        }
        forward.step()?;
        if n.is_multiple_of(sim.record_every) {
            forward.record(&mut sim.receivers)?;
        }
        if (n.is_multiple_of(100) || n + 1 == steps) && !forward.peak_velocity()?.is_finite() {
            bail!("the forward wavefield blew up by step {}", n + 1);
        }
    }
    forward.fetch(&mut sim.receivers)?;
    forward.download(&mut sim.solver)?;
    // As on the CPU: the misfit of the traces as `run` writes them, its derivative taken back
    // to the samples the solver recorded.
    sim.finish();
    let (misfit, residuals) = adjoint_source(&sim.receivers)?;
    ensure!(
        misfit.is_finite() && residuals.iter().flatten().all(|r| r.is_finite()),
        "the misfit or its derivative isn't finite; check the observed data"
    );
    let residuals = sim.recorded_residuals(residuals);
    let pressure = PressureResiduals::new(sim, &residuals);

    // The adjoint run, from the state Adjoint::new starts it in.
    let mut back = {
        let adjoint = Adjoint::new(sim, residuals)?;
        let mut back = GpuSolver::new(device.clone(), &adjoint.solver)?;
        let ([velocity, stress], values) = adjoint.injections();
        let columns = velocity.len() + stress.len();
        let velocity: Vec<(&PointOp, usize)> = velocity.into_iter().zip(0..).collect();
        let stress: Vec<(&PointOp, usize)> = stress.into_iter().zip(velocity.len()..).collect();
        back.attach_ops([&velocity, &stress], columns, &values, steps, &[], 0)?;
        back.args.skip = 1;
        back
    };

    let mut c = Correlations::new(device, &back, &pressure, every, len, g.interior_len())?;
    c.add_pressure(&pressure, steps)?;
    for (s, checkpoint) in checkpoints.iter().enumerate().rev() {
        let start = s * every;
        let end = (start + every).min(steps);
        // The velocities before the segment's first step, then after each of its steps.
        forward.upload_state(checkpoint)?;
        c.keep(&forward, 0)?;
        for i in 0..end - start {
            forward.step()?;
            c.keep(&forward, i + 1)?;
        }
        for n in (start..end).rev() {
            let (before, after) = (n - start, n - start + 1);
            // λ̂σ⁽ⁿ⁺¹⁾ = -W ε̃⁽ⁿ⁺¹⁾ with the forward strain rate at n + ½.
            c.strain_rate(&forward, after)?;
            back.step()?;
            ensure!(
                steps.checked_sub(back.step) == Some(n),
                "adjoint and forward out of step"
            );
            c.correlate(&back, before, after)?;
            // ε̃⁽ⁿ⁾ = ε̃⁽ⁿ⁺¹⁾ + e(ṽ⁽ⁿ⁺½⁾) - W⁻¹ r_σ⁽ⁿ⁾.
            c.add_adjoint_strain_rate(&back)?;
            c.add_pressure(&pressure, n)?;
        }
        if !back.peak_velocity()?.is_finite() {
            bail!("the adjoint wavefield blew up");
        }
    }

    let (sums, illumination) = c.finish()?;
    // On the host: the run, the checkpoints and the sums.
    let held: usize = checkpoints
        .iter()
        .map(crate::solver::State::memory_bytes)
        .sum();
    let memory_bytes = sim.memory_bytes() + held + 64 * len + 8 * illumination.len();
    Ok(Run {
        gradient: to_nodes(&sim.solver, &sums),
        misfit,
        checkpoint_every: every,
        memory_bytes,
        illumination,
        strategy: Strategy::Checkpoints { every: Some(every) },
    })
}

/// The arrays a gradient run correlates into, on the device, and a segment's velocities.
struct Correlations<D: Device> {
    device: Shared<D>,
    /// What the kernels read: the adjoint run's arguments with these arrays.
    args: Args,
    /// The forward velocities before the segment's first step and after each of its steps.
    segment: Vec<[u64; 3]>,
    /// The forward strain rate and the adjoint strain.
    strain: [u64; 6],
    adjoint_strain: [u64; 6],
    /// The pressure residuals' terms and values (by sample), if pressure was recorded.
    pressure: Option<Args>,
    len: usize,
    interior: usize,
}

impl<D: Device> Correlations<D> {
    fn new(
        mut device: Shared<D>,
        back: &GpuSolver<Shared<D>>,
        pressure: &PressureResiduals,
        every: usize,
        len: usize,
        interior: usize,
    ) -> Result<Self> {
        let need = 4 * (12 + 3 * (every + 1)) * len + 64 * len + 8 * interior;
        if let Some(free) = device.free_bytes()? {
            ensure!(
                need <= free,
                "the gradient's correlations need {:.2} GiB more on the GPU, which has {:.2} \
                 GiB free; a shorter checkpoint interval needs less",
                need as f64 / (1u64 << 30) as f64,
                free as f64 / (1u64 << 30) as f64
            );
        }
        let zeros = |device: &mut Shared<D>, bytes: usize| -> Result<u64> {
            let at = device.alloc(bytes)?;
            device.upload(at, &vec![0u8; bytes])?;
            Ok(at)
        };
        let mut args = back.args;
        let strain = std::array::from_fn(|_| 0);
        let mut c = Self {
            args,
            segment: Vec::with_capacity(every + 1),
            strain,
            adjoint_strain: strain,
            pressure: None,
            len,
            interior,
            device,
        };
        for a in &mut c.strain {
            *a = c.device.alloc(4 * len)?;
        }
        for a in &mut c.adjoint_strain {
            *a = zeros(&mut c.device, 4 * len)?;
        }
        for _ in 0..=every {
            let slot = [(); 3].map(|_| c.device.alloc(4 * len));
            let [x, y, z] = slot;
            c.segment.push([x?, y?, z?]);
        }
        args.sums = zeros(&mut c.device, 64 * len)?;
        args.lit = zeros(&mut c.device, 8 * interior)?;
        args.volume = back.grid.cell_volume();
        args.ea = c.adjoint_strain;
        c.args = args;

        let (terms, values, samples) = pressure.parts();
        if let Some(values) = values {
            let receivers = terms.len();
            let (mut field, mut index, mut weight, mut source) = (vec![], vec![], vec![], vec![]);
            for (r, terms) in terms.iter().enumerate() {
                for &(component, p, w) in terms {
                    field.push(component as i32);
                    index.push(p as i32);
                    weight.push(w);
                    source.push(r as i32);
                }
            }
            // By sample, a row of the receivers' values each.
            let mut by_sample = vec![0.0f32; values.len()];
            for r in 0..receivers {
                for s in 0..samples {
                    by_sample[s * receivers + r] = values[r * samples + s];
                }
            }
            ensure!(
                field.len() < i32::MAX as usize && by_sample.len() < i32::MAX as usize,
                "too many pressure terms"
            );
            let mut put = |data: &[u8]| -> Result<u64> {
                let at = c.device.alloc(data.len())?;
                c.device.upload(at, data)?;
                Ok(at)
            };
            c.pressure = Some(Args {
                term_field: put(bytes(&field))?,
                term_index: put(bytes(&index))?,
                term_weight: put(bytes(&weight))?,
                term_source: put(bytes(&source))?,
                values: put(bytes(&by_sample))?,
                sources: receivers as i32,
                from: 0,
                to: field.len() as i32,
                which: 1,
                skip: 0,
                e: c.adjoint_strain,
                ..c.args
            });
        }
        Ok(c)
    }

    /// Copies the forward run's velocities into the segment's slot `i`.
    fn keep(&mut self, forward: &GpuSolver<Shared<D>>, i: usize) -> Result<()> {
        let args = Args {
            v: self.segment[i],
            ..forward.args
        };
        let grid = [(self.len as u32).div_ceil(LINE), 1, 1];
        self.device.launch(Kernel::Copy, grid, [LINE, 1, 1], &args)
    }

    /// The forward strain rate of the velocities in slot `after`.
    fn strain_rate(&mut self, forward: &GpuSolver<Shared<D>>, after: usize) -> Result<()> {
        let args = Args {
            v: self.segment[after],
            e: self.strain,
            which: 0,
            ..self.args
        };
        self.sweep(forward, Kernel::Strain, &args)
    }

    /// Adds the correlations of the forward step whose velocities are in slots `before` and
    /// `after` with the adjoint step `back` has just taken.
    fn correlate(
        &mut self,
        back: &GpuSolver<Shared<D>>,
        before: usize,
        after: usize,
    ) -> Result<()> {
        let args = Args {
            f: back.args.f,
            e: self.strain,
            before: self.segment[before],
            after: self.segment[after],
            ..self.args
        };
        self.sweep(back, Kernel::Correlate, &args)
    }

    /// Adds the strain rate of the adjoint velocities to the adjoint strain.
    fn add_adjoint_strain_rate(&mut self, back: &GpuSolver<Shared<D>>) -> Result<()> {
        let [vx, vy, vz, ..] = back.args.f;
        let args = Args {
            v: [vx, vy, vz],
            e: self.adjoint_strain,
            which: 1,
            ..self.args
        };
        self.sweep(back, Kernel::Strain, &args)
    }

    /// [`PressureResiduals::add`] for `n`.
    fn add_pressure(&mut self, pressure: &PressureResiduals, n: usize) -> Result<()> {
        let (Some(args), Some(sample)) = (&self.pressure, pressure.sample_for(n)) else {
            return Ok(());
        };
        let args = Args {
            row: sample as i32,
            ..*args
        };
        self.device.launch(Kernel::Inject, [1; 3], [1; 3], &args)
    }

    fn sweep(&mut self, on: &GpuSolver<Shared<D>>, kernel: Kernel, args: &Args) -> Result<()> {
        let grid = std::array::from_fn(|a| (on.grid.n[a] as u32).div_ceil(on.block[a]));
        self.device.launch(kernel, grid, on.block, args)
    }

    /// The sums and the illumination.
    fn finish(mut self) -> Result<(Sums, Vec<f64>)> {
        let mut sums: Sums = vec![[0.0; 8]; self.len];
        self.device.download(self.args.sums, bytes_mut(&mut sums))?;
        let mut illumination = vec![0.0f64; self.interior];
        self.device
            .download(self.args.lit, bytes_mut(&mut illumination))?;
        Ok((sums, illumination))
    }
}
