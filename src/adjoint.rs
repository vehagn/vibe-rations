//! The adjoint wavefield: residuals at the receivers, propagated backwards in time
//! (Tarantola 1984; Tromp, Tape & Liu 2005; Fichtner et al. 2006).
//!
//! The adjoint run is the forward solver itself, stepped from rest in reversed order, with each
//! residual entering through the transpose of the receiver that recorded it
//! ([`Solver::adjoint_op`]). That is the exact transpose of the discrete forward scheme, not
//! only of the wave equation, wherever the scheme is self-adjoint: everywhere outside the
//! absorbing layers. Why, in brief:
//!
//! One forward step n is v⁽ⁿ⁺½⁾ = v⁽ⁿ⁻½⁾ + B D σ⁽ⁿ⁾ + B f⁽ⁿ⁾, then σ⁽ⁿ⁺¹⁾ = σ⁽ⁿ⁾ + M D′ v⁽ⁿ⁺½⁾ + g⁽ⁿ⁾,
//! with B = dt/ρ, M = dt·C, and D′ = -W_σ⁻¹ Dᵀ W_v for the cell weights W (volume times the
//! cell's share, half on the free surface): the staggered derivatives are negative adjoints of
//! each other, which is also what makes the scheme conserve energy. Scaling each row by W B⁻¹
//! or W M⁻¹ makes the stepping symmetric except for the time differences, and transposing it
//! turns them round: the adjoint variables λ obey
//!
//! λv⁽ⁿ⁺½⁾ = λv⁽ⁿ⁺³ᐟ²⁾ - B D λσ⁽ⁿ⁺¹⁾ + B W_v⁻¹ rᵥ⁽ⁿ⁾,
//! λσ⁽ⁿ⁾ = λσ⁽ⁿ⁺¹⁾ - M D′ λv⁽ⁿ⁺½⁾ + M W_σ⁻¹ r_σ⁽ⁿ⁾,
//!
//! where r is the derivative of the misfit with respect to what the receivers read, Rᵀ times
//! the residuals. With ṽ = λv and σ̃ = -λσ this is the forward step again, taken backwards:
//! adjoint step m = N - 1 - n reads σ̃ at n + 1, writes ṽ at n + ½ and σ̃ at n, and injects a
//! velocity residual recorded after forward step n as a unit force (B W_v⁻¹ Rᵀ), and a pressure
//! residual, which the forward run read from σ⁽ⁿ⁺¹⁾, one step earlier as an isotropic stress
//! source (-M W_σ⁻¹ Rᵀ). A pressure sample from the last forward step is the adjoint run's
//! initial state.
//!
//! So after adjoint step m, the adjoint velocities pair with the forward velocities after
//! forward step n = N - 1 - m, at (n + ½)·dt, and the adjoint stresses (with their sign
//! flipped) with the forward stresses before it, at n·dt. The derivative of the misfit with
//! respect to a force f⁽ⁿ⁾ at a node is the adjoint velocity read there by the receiver's
//! operator after adjoint step m.
//!
//! The absorbing layers are not self-adjoint in themselves: the CPML's memory variables and the
//! sponge transpose into something else. Seen from outside them, they nearly are: with the
//! source and receivers clear of the layers, the dot-product test in `tests/physics.rs` holds
//! to single-precision round-off (below 1e-7), as it does with rigid walls. A receiver inside
//! a layer breaks it by up to 1.4e-3, and so does, less, a velocity receiver on the first node
//! outside one, whose interpolation reaches into the layer
//! ([`Simulation::receivers_next_to_layers`]; docs/research-log.md).

use anyhow::{Result, bail, ensure};

use crate::receiver::Component;
use crate::simulation::Simulation;
use crate::solver::{Fields, PointOp, Solver};

/// A run of the adjoint wavefield for a forward [`Simulation`]'s receivers.
pub struct Adjoint {
    /// The solver, holding the adjoint wavefield.
    pub solver: Solver,
    /// Steps of the forward run, which the adjoint run retraces in reverse.
    pub steps: usize,
    record_every: usize,
    samples: usize,
    components: Vec<Component>,
    /// How each residual enters: component-major, like the receivers' traces.
    ops: Vec<Vec<PointOp>>,
    residuals: Vec<Vec<f32>>,
}

impl Adjoint {
    /// Bytes its solver and residuals hold.
    pub fn memory_bytes(&self) -> usize {
        self.solver.memory_bytes() + 4 * self.residuals.iter().map(Vec::len).sum::<usize>()
    }

    /// Sets up the adjoint run for `sim`'s receivers, at rest, with `residuals` in the layout
    /// of [`crate::receiver::Receivers::traces`]: one array per recorded component,
    /// (receivers, samples), each the derivative of the misfit with respect to that sample as
    /// the solver recorded it. For a misfit of un-warped traces (time dispersion removed),
    /// pass them through [`Simulation::recorded_residuals`] first.
    pub fn new(sim: &Simulation, residuals: Vec<Vec<f32>>) -> Result<Self> {
        let r = &sim.receivers;
        ensure!(
            residuals.len() == r.components.len(),
            "{} residual arrays for {} recorded components",
            residuals.len(),
            r.components.len()
        );
        for (c, res) in r.components.iter().zip(&residuals) {
            ensure!(
                res.len() == r.nodes.len() * r.samples,
                "the {} residuals hold {} values; {} receivers × {} samples need {}",
                c.name(),
                res.len(),
                r.nodes.len(),
                r.samples,
                r.nodes.len() * r.samples
            );
        }
        let solver = sim.solver.at_rest();
        let ops = r
            .components
            .iter()
            .map(|&c| r.nodes.iter().map(|&n| solver.adjoint_op(c, n)).collect())
            .collect();
        let mut adjoint = Self {
            solver,
            steps: sim.steps,
            record_every: sim.record_every,
            samples: r.samples,
            components: r.components.clone(),
            ops,
            residuals,
        };
        // Stresses recorded after the last forward step have nothing after them to pass
        // through: they are where the adjoint run starts.
        if let Some(last) = adjoint.steps.checked_sub(1) {
            let Self {
                solver,
                record_every,
                samples,
                components,
                ops,
                residuals,
                ..
            } = &mut adjoint;
            let at = Samples {
                record_every: *record_every,
                samples: *samples,
                components,
                ops,
                residuals,
            };
            at.inject(&mut solver.fields, last, true);
            if solver.free_surface {
                solver.mirror_stress();
            }
        }
        Ok(adjoint)
    }

    /// The forward step the wavefield now pairs with, if a step has been taken: its velocities
    /// are the adjoint of the forward velocities after forward step n, at (n + ½)·dt, and its
    /// stresses, with their sign flipped, the adjoint of the forward stresses at n·dt.
    pub fn forward_step(&self) -> Option<usize> {
        let m = self.solver.step.checked_sub(1)?;
        self.steps.checked_sub(m + 1)
    }

    /// Steps still to take.
    pub fn remaining(&self) -> usize {
        self.steps.saturating_sub(self.solver.step)
    }

    /// Takes the next step back.
    pub fn step(&mut self) {
        assert!(self.remaining() > 0, "the adjoint run is over");
        let n = self.steps - 1 - self.solver.step;
        let Self {
            solver,
            record_every,
            samples,
            components,
            ops,
            residuals,
            ..
        } = self;
        let at = Samples {
            record_every: *record_every,
            samples: *samples,
            components,
            ops,
            residuals,
        };
        solver.step_with(
            |f, _| at.inject(f, n, false),
            |f, _| {
                if let Some(earlier) = n.checked_sub(1) {
                    at.inject(f, earlier, true);
                }
            },
        );
    }

    /// Takes every remaining step, calling `after` with the forward step the wavefield then
    /// pairs with ([`Adjoint::forward_step`]) and the solver.
    pub fn run(&mut self, mut after: impl FnMut(usize, &Solver)) -> Result<()> {
        while self.remaining() > 0 {
            self.step();
            let n = self.forward_step().expect("a step was taken");
            if self.solver.step.is_multiple_of(100) || n == 0 {
                let peak = self.solver.peak_velocity();
                if !peak.is_finite() {
                    bail!("the adjoint wavefield blew up by step {}", self.solver.step);
                }
            }
            after(n, &self.solver);
        }
        Ok(())
    }
}

/// The residuals, borrowed apart from the solver so they can be injected during its step.
struct Samples<'a> {
    record_every: usize,
    samples: usize,
    components: &'a [Component],
    ops: &'a [Vec<PointOp>],
    residuals: &'a [Vec<f32>],
}

impl Samples<'_> {
    /// Adds the residuals of the samples recorded after forward step `n`, if it recorded any:
    /// those of the stress components (`stress`) or of the velocity components.
    fn inject(&self, fields: &mut Fields, n: usize, stress: bool) {
        if !n.is_multiple_of(self.record_every) {
            return;
        }
        let s = n / self.record_every;
        if s >= self.samples {
            return;
        }
        for ((&c, ops), res) in self.components.iter().zip(self.ops).zip(self.residuals) {
            if (c == Component::P) != stress {
                continue;
            }
            for (r, op) in ops.iter().enumerate() {
                let value = res[r * self.samples + s];
                if value != 0.0 {
                    op.inject(fields, value);
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    use crate::config::Config;
    use crate::solver::Field;
    use crate::source::{Mechanism, Source, Wavelet};

    const BASE: &str = r#"
        [grid]
        points = [24, 22, 20]
        spacing = 10.0
        [time]
        duration = 0.02
        [boundary]
        pml_width = 0
        [[model.layer]]
        top = 0.0
        vp = 2000.0
        vs = 1000.0
        rho = 2000.0
        [[source]]
        position = [100.0, 100.0, 100.0]
        force = [0.0, 0.0, 1.0]
        wavelet = { kind = "ricker", frequency = 20.0 }
        [receivers]
        components = ["vx", "vz", "p"]
        [[receivers.point]]
        position = [120.0, 110.0, 90.0]
        [[receivers.point]]
        position = [80.0, 100.0, 0.0]
    "#;

    fn setup(text: &str) -> Simulation {
        Simulation::new(&Config::from_toml(text).unwrap(), Path::new(".")).unwrap()
    }

    fn zeros(sim: &Simulation) -> Vec<Vec<f32>> {
        let r = &sim.receivers;
        vec![vec![0.0; r.nodes.len() * r.samples]; r.components.len()]
    }

    #[test]
    fn rejects_residuals_that_do_not_fit_the_receivers() {
        let sim = setup(BASE);
        let mut res = zeros(&sim);
        res.pop();
        assert!(Adjoint::new(&sim, res).is_err());
        let mut res = zeros(&sim);
        res[1].push(0.0);
        let e = format!("{:#}", Adjoint::new(&sim, res).err().unwrap());
        assert!(e.contains("vz residuals"), "{e}");
    }

    #[test]
    fn a_velocity_residual_is_a_unit_force() {
        let sim = setup(BASE);
        let wavelet = Wavelet::Ricker { f0: 1.0, t0: 0.0 };
        for (c, f) in [
            (Component::Vx, [1.0, 0.0, 0.0]),
            (Component::Vz, [0.0, 0.0, 1.0]),
        ] {
            for node in [[12, 11, 9], [8, 10, 0]] {
                let force = Source::new(&sim.solver, node, Mechanism::Force(f), wavelet).op;
                let op = sim.solver.adjoint_op(c, node);
                assert_eq!(op.terms.len(), force.terms.len());
                for (a, b) in op.terms.iter().zip(&force.terms) {
                    assert_eq!((a.0, a.1), (b.0, b.1));
                    assert!((a.2 - b.2).abs() <= 1e-6 * b.2.abs(), "{c:?} at {node:?}");
                }
            }
        }
    }

    /// Below the surface a pressure residual loads the normal stresses with dt·K/V; on it, σxx
    /// and σyy with the plane-stress modulus over half a cell, and σzz, which stays zero, not.
    #[test]
    fn a_pressure_residual_is_an_isotropic_stress_source() {
        let sim = setup(BASE);
        let s = &sim.solver;
        let (dt, volume) = (s.dt, 1000.0);
        let (l2m, lam) = (8e9, 4e9);
        let close = |a: f32, b: f64| ((a as f64 - b) / b).abs() < 1e-6;
        let op = s.adjoint_op(Component::P, [12, 11, 9]);
        let k = lam + 2.0 * 2e9 / 3.0;
        assert_eq!(op.terms.len(), 3);
        assert!(
            op.terms.iter().all(|t| close(t.2, dt * k / volume)),
            "{op:?}"
        );
        let op = s.adjoint_op(Component::P, [8, 10, 0]);
        let q = lam * lam / l2m;
        let horizontal = (l2m - q) + (lam - q);
        assert_eq!(
            op.terms.iter().map(|t| t.0).collect::<Vec<_>>(),
            vec![Field::Sxx, Field::Syy]
        );
        assert!(
            op.terms
                .iter()
                .all(|t| close(t.2, dt * horizontal / (3.0 * 0.5 * volume))),
            "{op:?}"
        );
    }

    /// Adjoint steps taken when the residual of sample `s` of component `c` has entered, 0 if
    /// it is the initial state.
    fn entered(sim: &Simulation, c: usize, s: usize) -> usize {
        let mut res = zeros(sim);
        res[c][s] = 1.0;
        let mut adjoint = Adjoint::new(sim, res).unwrap();
        let quiet = |f: &Fields| {
            [Field::Vx, Field::Vz, Field::Sxx]
                .iter()
                .all(|&k| f.get(k).iter().all(|&x| x == 0.0))
        };
        while quiet(&adjoint.solver.fields) {
            adjoint.step();
        }
        adjoint.solver.step
    }

    /// A velocity sample from forward step n enters during adjoint step N - 1 - n; a pressure
    /// sample, read from the stresses after step n, one adjoint step earlier, so one from the
    /// last forward step is the initial state.
    #[test]
    fn residuals_enter_in_reverse_order() {
        for every in [1, 3] {
            let sim = setup(&format!("{BASE}\n[output]\nrecord_every = {every}"));
            let steps = sim.steps;
            for s in [0, 1, sim.receivers.samples - 1] {
                let n = every * s;
                // Components: vx, vz, p.
                assert_eq!(
                    entered(&sim, 0, s),
                    steps - n,
                    "vx sample {s}, every {every}"
                );
                assert_eq!(
                    entered(&sim, 2, s),
                    steps - 1 - n,
                    "p sample {s}, every {every}"
                );
            }
        }
    }

    #[test]
    fn pairs_each_step_with_the_forward_step_it_retraces() {
        let sim = setup(BASE);
        let mut adjoint = Adjoint::new(&sim, zeros(&sim)).unwrap();
        assert_eq!(adjoint.forward_step(), None);
        let mut seen = Vec::new();
        adjoint.run(|n, _| seen.push(n)).unwrap();
        let want: Vec<usize> = (0..sim.steps).rev().collect();
        assert_eq!(seen, want);
        assert_eq!(adjoint.remaining(), 0);
    }
}
