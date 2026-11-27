//! The Hessian's products: the change of the gradient for a change of the model, H δm, as the
//! discrete run computes the gradient (the tangent-linear model of
//! [`crate::gradient::gradient`], forward over reverse; in the continuous setting, Fichtner &
//! Trampert 2011's second-order adjoints).
//!
//! The gradient is g = T S(u, ũ): sums S over the run, bilinear in the forward wavefield u and
//! the adjoint wavefield ũ, taken to the nodes by T, whose factors depend on the model only
//! through the free surface's plane-stress moduli and the shear moduli's harmonic means. So
//!
//! H δm = T S(δu, ũ) + T S(u, δũ) + δT S(u, ũ),
//!
//! with δu the Born run ([`crate::born`]) and δũ the adjoint run's own change: the solver's
//! step on δṽ, δσ̃, driven by
//!
//! - the change of the residuals, the misfit's second derivative in the samples applied to
//!   the Born traces (for L2, `Misfit::gauss_newton_source`), entering as residuals do; and
//! - the adjoint run scattered by the change of the material, as the Born run scatters the
//!   forward one: (δB/B)(ṽ⁽ⁿ⁺½⁾ - ṽ⁽ⁿ⁺³ᐟ²⁾) on the velocities and δM (ε̃⁽ⁿ⁾ - ε̃⁽ⁿ⁺¹⁾) on the
//!   stresses, ε̃ the adjoint strain, and δM ε̃ᴺ in the stresses it starts from. A residual
//!   enters through the material too (B for a velocity, M for pressure); these differences
//!   include what the residuals added, so they take that in.
//!
//! Only the first drives δũ at zero residual, where ũ = 0: what is left is the Gauss–Newton
//! product, JᵀWJ δm. The rest, the second-order part, grows with the residual.
//!
//! The forward and Born runs are replayed together from checkpoints, so the product is exact
//! up to single-precision round-off outside the absorbing layers, where the gradient is, and
//! the gradient that comes with it is [`crate::gradient::gradient`]'s, bit for bit.

use anyhow::{Result, bail, ensure};

use crate::adjoint::Adjoint;
use crate::born::{Change, Scattering};
use crate::gradient::{
    Gradient, PressureResiduals, Sums, checkpoint_interval, correlate, to_nodes, velocities,
};
use crate::model::Model;
use crate::receiver::Receivers;
use crate::simulation::Simulation;
use crate::solver::{Solver, State};

/// What a Hessian product found.
pub struct Product {
    /// The gradient, as [`crate::gradient::gradient`] gives it with the same checkpoints.
    pub gradient: Gradient,
    /// H δm: the change of the gradient, in ρ, λ and μ, for the change of the model.
    pub product: Gradient,
    /// The misfit the adjoint source came from, as the closure returned it.
    pub misfit: f64,
    /// The most memory held at once, counted from its arrays: the forward, Born and two
    /// adjoint runs, both checkpoints, a replayed segment of both, the sums and the strains.
    pub memory_bytes: usize,
}

impl Product {
    /// The product in vp, vs and ρ, for the change of them `speeds` (in that order) that the
    /// change in ρ, λ and μ was made from ([`Change::from_velocities`]) at the nodes of
    /// `model`: the product taken through [`Gradient::velocities`], plus that conversion's own
    /// change applied to the gradient, since its factors depend on the model.
    pub fn velocities(&self, model: &Model, speeds: [&[f64]; 3]) -> [Vec<f64>; 3] {
        let [mut vp, mut vs, mut rho] = self.product.velocities(model);
        let [dvp, dvs, drho] = speeds;
        let g = &self.gradient;
        for i in 0..vp.len() {
            let (a, b, r) = (model.vp[i] as f64, model.vs[i] as f64, model.rho[i] as f64);
            let (gl, gm) = (g.lambda[i], g.mu[i]);
            vp[i] += 2.0 * (drho[i] * a + r * dvp[i]) * gl;
            vs[i] += 2.0 * (drho[i] * b + r * dvs[i]) * (gm - 2.0 * gl);
            rho[i] += (2.0 * a * dvp[i] - 4.0 * b * dvs[i]) * gl + 2.0 * b * dvs[i] * gm;
        }
        [vp, vs, rho]
    }
}

/// Runs `sim` (fresh, at step 0) forward with a Born run beside it for the model change
/// `change`, keeping a checkpoint of both every `checkpoint_every` steps (by default
/// [`checkpoint_interval`]'s); hands the receivers to `adjoint_source`, which returns the
/// misfit and its derivative with respect to every sample, as for
/// [`crate::gradient::gradient`], and the Born traces to `curvature`, which returns the
/// misfit's second derivative in the samples applied to them; then runs the adjoint and its
/// change backwards, replaying both forward runs a segment at a time, and correlates the
/// four. Both closures see traces as `run` leaves them, un-warped when the run removes time
/// dispersion.
///
/// As for the gradient, the time step and the absorbing layers stay as `sim` set them up.
/// Afterwards `sim` holds the recorded traces, its solver wherever the last replay left it.
pub fn hessian_product(
    sim: &mut Simulation,
    change: &Change,
    adjoint_source: impl FnOnce(&Receivers) -> Result<(f64, Vec<Vec<f32>>)>,
    curvature: impl FnOnce(&[Vec<f32>]) -> Result<Vec<Vec<f32>>>,
    checkpoint_every: Option<usize>,
) -> Result<Product> {
    let steps = sim.steps;
    let every = checkpoint_interval(steps, checkpoint_every);
    let mut scattering = Scattering::new(sim, change)?;

    // Forward, keeping both runs' states at the start of every segment.
    let mut checkpoints: Vec<[State; 2]> = Vec::with_capacity(steps.div_ceil(every));
    for n in 0..steps {
        if n.is_multiple_of(every) {
            checkpoints.push([sim.solver.state(), scattering.solver.state()]);
        }
        scattering.advance(sim, true);
        if (n.is_multiple_of(100) || n + 1 == steps)
            && !(sim.solver.peak_velocity().is_finite()
                && scattering.solver.peak_velocity().is_finite())
        {
            bail!("the forward or Born wavefield blew up by step {}", n + 1);
        }
    }
    sim.finish();
    let (misfit, residuals) = adjoint_source(&sim.receivers)?;
    let changes = curvature(&scattering.traces(sim)?)?;
    ensure!(
        misfit.is_finite()
            && residuals
                .iter()
                .chain(&changes)
                .flatten()
                .all(|r| r.is_finite()),
        "the misfit or its derivatives aren't finite; check the observed data"
    );
    let residuals = sim.recorded_residuals(residuals);
    let changes = sim.recorded_residuals(changes);
    let pressure = PressureResiduals::new(sim, &residuals);
    let changed_pressure = PressureResiduals::new(sim, &changes);
    let mut adjoint = Adjoint::new(sim, residuals)?;
    let mut changed = Adjoint::new(sim, changes)?;

    let g = sim.solver.grid;
    let len = g.padded_len();
    let six = || std::array::from_fn::<_, 6, _>(|_| vec![0.0f32; len]);
    // The adjoint strains ε̃ and δε̃, from ε̃ᴺ and δε̃ᴺ; the changed run starts with δM ε̃ᴺ
    // besides its own residuals.
    let (mut adjoint_strain, mut changed_adjoint_strain) = (six(), six());
    pressure.add(&mut adjoint_strain, steps);
    changed_pressure.add(&mut changed_adjoint_strain, steps);
    {
        let (forward, material) = (&sim.solver, &scattering.material);
        changed.add_to_start(|f| material.add_stress(forward, f, &adjoint_strain));
    }
    // The forward and Born strain rates, and ε̃⁽ⁿ⁾ - ε̃⁽ⁿ⁺¹⁾.
    let (mut strain, mut born_strain, mut adjoint_step) = (six(), six(), six());
    // The adjoint velocities before the step under way.
    let mut ahead: [Vec<f32>; 3] = std::array::from_fn(|_| vec![0.0f32; len]);
    let mut sums: Sums = vec![[0.0; 8]; len];
    let mut changed_sums: Sums = vec![[0.0; 8]; len];
    let mut unused = vec![0.0f64; g.interior_len()];
    let fixed = sim.memory_bytes()
        + scattering.memory_bytes()
        + adjoint.memory_bytes()
        + changed.memory_bytes()
        + pressure.memory_bytes()
        + changed_pressure.memory_bytes()
        + 2 * std::mem::size_of_val(sums.as_slice())
        + std::mem::size_of_val(unused.as_slice())
        + 4 * len * (5 * 6 + 3);
    let held: usize = checkpoints.iter().flatten().map(State::memory_bytes).sum();
    let mut peak_bytes = 0;

    let mut segment: Vec<[[Vec<f32>; 3]; 2]> = Vec::with_capacity(every + 1);
    for (s, [forward_state, born_state]) in checkpoints.iter().enumerate().rev() {
        let start = s * every;
        let end = (start + every).min(steps);
        // Both runs' velocities before the segment's first step, then after each of its steps.
        sim.solver.restore(forward_state);
        scattering.solver.restore(born_state);
        segment.clear();
        segment.push([velocities(&sim.solver), velocities(&scattering.solver)]);
        for _ in start..end {
            scattering.advance(sim, false);
            segment.push([velocities(&sim.solver), velocities(&scattering.solver)]);
        }
        peak_bytes = peak_bytes.max(fixed + held + 4 * 6 * len * segment.len());

        let (forward, material) = (&sim.solver, &scattering.material);
        for n in (start..end).rev() {
            let [before, after] =
                [n - start, n - start + 1].map(|i| segment[i][0].each_ref().map(Vec::as_slice));
            let [born_before, born_after] =
                [n - start, n - start + 1].map(|i| segment[i][1].each_ref().map(Vec::as_slice));
            forward.strain_rate(after, &mut strain);
            forward.strain_rate(born_after, &mut born_strain);
            for (a, v) in ahead.iter_mut().zip(adjoint.solver.fields.velocities()) {
                a.copy_from_slice(v);
            }
            adjoint.step();
            assert_eq!(
                adjoint.forward_step(),
                Some(n),
                "adjoint and forward out of step"
            );
            // e(ṽ⁽ⁿ⁺½⁾) - W⁻¹ r_σ⁽ⁿ⁾: what M multiplied in the step's stresses, residuals
            // included, as in the adjoint strain's update below.
            let back = adjoint.solver.fields.velocities();
            forward.strain_rate(back, &mut adjoint_step);
            pressure.add(&mut adjoint_step, n);
            changed.step_adding(
                |f| material.add_velocity(f, &ahead, back),
                |f| material.add_stress(forward, f, &adjoint_step),
            );
            // S(u, ũ), as the gradient takes it; S(u, δũ) + S(δu, ũ).
            let adjoint_strains = [&adjoint_strain, &strain];
            correlate(
                forward,
                &adjoint.solver,
                &mut sums,
                &mut unused,
                adjoint_strains,
                before,
                after,
            );
            correlate(
                forward,
                &changed.solver,
                &mut changed_sums,
                &mut unused,
                [&changed_adjoint_strain, &strain],
                before,
                after,
            );
            correlate(
                forward,
                &adjoint.solver,
                &mut changed_sums,
                &mut unused,
                [&adjoint_strain, &born_strain],
                born_before,
                born_after,
            );
            // ε̃⁽ⁿ⁾ and δε̃⁽ⁿ⁾, as the gradient sums ε̃.
            let a = &adjoint.solver;
            a.add_strain_rate(a.fields.velocities(), &mut adjoint_strain);
            pressure.add(&mut adjoint_strain, n);
            let c = &changed.solver;
            c.add_strain_rate(c.fields.velocities(), &mut changed_adjoint_strain);
            changed_pressure.add(&mut changed_adjoint_strain, n);
        }
        ensure!(
            adjoint.solver.peak_velocity().is_finite()
                && changed.solver.peak_velocity().is_finite(),
            "an adjoint wavefield blew up"
        );
    }

    drop((checkpoints, segment));
    let gradient = to_nodes(&sim.solver, &sums);
    let mut product = to_nodes(&sim.solver, &changed_sums);
    add_factors_change(&sim.solver, &sums, change, &mut product);
    Ok(Product {
        gradient,
        product,
        misfit,
        memory_bytes: peak_bytes,
    })
}

/// Adds δT S: the change of `gradient::to_nodes`'s factors for the model change `change`,
/// applied to the sums `sums`. Only two kinds depend on the model: on the free surface, the
/// plane-stress moduli's, a = l2m - lam²/l2m and b = lam - lam²/l2m, through r = lam/l2m; and
/// the shear moduli's harmonic means', H²/(4μ_q²) with H = 4/Σ 1/μ_q.
fn add_factors_change(s: &Solver, sums: &Sums, change: &Change, out: &mut Gradient) {
    let g = s.grid;
    let [nx, ny, nz] = g.n;
    let dt = s.dt;
    let model = &s.model;
    let m = &s.material;
    let node =
        |i: usize, j: usize, k: usize| g.interior_idx(i.min(nx - 1), j.min(ny - 1), k.min(nz - 1));
    let mu_at = |q: usize| model.rho[q] as f64 * (model.vs[q] as f64).powi(2);
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let p = g.idx(i, j, k);
                let here = node(i, j, k);
                let sum = &sums[p];
                let (l2, l) = (m.l2m[p] as f64, m.lam[p] as f64);
                if s.free_surface && k == 0 && l2 > 0.0 {
                    // ∂χ/∂l2m = ga (1 + r²) + gb r² and ∂χ/∂lam = -2r ga + gb (1 - 2r).
                    let dl2 = dt * (change.lambda[here] + 2.0 * change.mu[here]);
                    let dl = dt * change.lambda[here];
                    let r = l / l2;
                    let dr = (dl * l2 - l * dl2) / (l2 * l2);
                    let (ga, gb) = (sum[3], sum[4]);
                    let d_l2m = 2.0 * r * (ga + gb) * dr;
                    let d_lam = -2.0 * (ga + gb) * dr;
                    out.lambda[here] += dt * (d_l2m + d_lam);
                    out.mu[here] += 2.0 * dt * d_l2m;
                }
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
                    let dh: f64 = corners
                        .iter()
                        .zip(&mus)
                        .map(|(&q, &mu)| h * h / (4.0 * mu * mu) * change.mu[q])
                        .sum();
                    for (&q, &mu) in corners.iter().zip(&mus) {
                        let d =
                            h * dh / (2.0 * mu * mu) - h * h * change.mu[q] / (2.0 * mu.powi(3));
                        out.mu[q] += dt * sum[c] * d;
                    }
                }
            }
        }
    }
}
