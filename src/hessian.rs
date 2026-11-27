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
//!
//! # The matrix between a few nodes
//!
//! With the run written as A(m) u = s, its adjoint as Aᵀ ũ = Rᵀ r and δuᵦ = -A⁻¹ (∂ᵦA) u the
//! Born run of a unit change of the parameter b, the Hessian's entries are
//!
//! Hₐᵦ = (R δuₐ)ᵀ W (R δuᵦ) - ũᵀ (∂ₐA) δuᵦ - ũᵀ (∂ᵦA) δuₐ - ũᵀ (∂ₐ∂ᵦA) u:
//!
//! the Born traces' inner product (Gauss–Newton's JᵀWJ); B, with Bₐᵦ the gradient's
//! correlation at a taken with the Born wavefield of b in the forward one's place, which is
//! T S(δu, ũ) above; its transpose, which is the share of T S(u, δũ) that the adjoint run
//! scattered by the change drives (ũᵀ (∂ᵦA) A⁻¹ (∂ₐA) u read from its other end); and δT S.
//! So the entries between a set of nodes need no adjoint run of their own per column
//! ([`hessian_block`]): one forward and one adjoint run, their values kept at the positions
//! around the nodes for every step, then per column one Born run driven from the kept
//! forward values and correlated there with the kept adjoint ones. The second-order part,
//! B + Bᵀ + δT S, then comes at its own size instead of as the difference of the product
//! and JᵀWJ, which loses as many digits as it is smaller.

use anyhow::{Result, bail, ensure};

use rayon::prelude::*;

use crate::adjoint::Adjoint;
use crate::born::{Change, MaterialChange, Scattering, add_stress_at, add_velocity_at, unwarped};
use crate::gradient::{
    Gradient, PressureResiduals, Sums, checkpoint_interval, correlate, correlate_at, to_nodes,
    velocities,
};
use crate::grid::Grid;
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
    /// The product's two wavefield terms, in ρ, λ and μ: T S(δu, ũ), the Born run against
    /// the adjoint run, and T S(u, δũ), the forward run against the adjoint run's change
    /// (Fichtner & Trampert 2011's first two Hessian kernels; the second holds the
    /// Gauss–Newton product, the part of δũ the Born traces drive). The product is their sum
    /// and δT S(u, ũ).
    pub terms: [Gradient; 2],
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
    let mut born_sums: Sums = vec![[0.0; 8]; len];
    let mut unused = vec![0.0f64; g.interior_len()];
    let fixed = sim.memory_bytes()
        + scattering.memory_bytes()
        + adjoint.memory_bytes()
        + changed.memory_bytes()
        + pressure.memory_bytes()
        + changed_pressure.memory_bytes()
        + 3 * std::mem::size_of_val(sums.as_slice())
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
            // S(u, ũ), as the gradient takes it; S(u, δũ); S(δu, ũ).
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
                &mut born_sums,
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
    let terms = [
        to_nodes(&sim.solver, &born_sums),
        to_nodes(&sim.solver, &changed_sums),
    ];
    let [scattered, changed] = &terms;
    let add = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x + y).collect();
    let mut product = Gradient {
        rho: add(&changed.rho, &scattered.rho),
        lambda: add(&changed.lambda, &scattered.lambda),
        mu: add(&changed.mu, &scattered.mu),
    };
    add_factors_change(&sim.solver, &sums, change, &mut product);
    Ok(Product {
        gradient,
        product,
        terms,
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

/// The Hessian between a set of nodes, entry by entry, in ρ, λ and μ ([`hessian_block`]).
///
/// Each matrix is 3n × 3n for n nodes, row-major, a parameter's nodes together: the entry
/// of parameter p (ρ, λ, μ: 0, 1, 2) at the a-th node against parameter q at the b-th is at
/// `(p * n + a) * 3 * n + q * n + b`.
pub struct Block {
    /// The nodes, interior indices as [`Model`] counts them, in the order given.
    pub nodes: Vec<usize>,
    /// The gradient, as [`crate::gradient::gradient`] gives it with the same checkpoints.
    pub gradient: Gradient,
    /// The misfit the adjoint source came from, as the closure returned it.
    pub misfit: f64,
    /// J restricted to the nodes: for each column of the matrices, the Born traces of a unit
    /// change of its parameter at its node, in the layout of [`Receivers::traces`] and
    /// un-warped like the forward run's.
    pub jacobian: Vec<Vec<Vec<f32>>>,
    /// JᵀWJ: the inner products of the Born traces, through the misfit's second derivative
    /// in the samples.
    pub gauss_newton: Vec<f64>,
    /// B: a column is T S(δu, ũ) at the nodes, δu the Born run of a unit change of the
    /// column's parameter at its node. Not symmetric: the Hessian holds it and its transpose.
    pub scattered: Vec<f64>,
    /// δT S(u, ũ): the change of [`crate::gradient::to_nodes`]'s factors applied to the
    /// gradient's sums, between nodes that share a shear position or lie on a free surface.
    pub factors: Vec<f64>,
    /// The most memory held at once, counted from its arrays: the gradient run's, or the
    /// Born runs under way with every column's traces.
    pub memory_bytes: usize,
}

impl Block {
    /// The matrices' side: three parameters at each node.
    pub fn side(&self) -> usize {
        3 * self.nodes.len()
    }

    /// The second-order part of the Hessian between the nodes, B + Bᵀ + δT S: what the
    /// adjoint wavefield, and so the residual, adds to Gauss–Newton's JᵀWJ.
    pub fn second_order(&self) -> Vec<f64> {
        let side = self.side();
        (0..side * side)
            .map(|i| {
                let (a, b) = (i / side, i % side);
                self.scattered[i] + self.scattered[b * side + a] + self.factors[i]
            })
            .collect()
    }

    /// The Hessian between the nodes, in ρ, λ and μ.
    pub fn hessian(&self) -> Vec<f64> {
        let mut h = self.second_order();
        h.iter_mut()
            .zip(&self.gauss_newton)
            .for_each(|(h, g)| *h += g);
        h
    }

    /// `matrix`, between the nodes in ρ, λ and μ, taken to vp, vs and ρ (0, 1, 2, as
    /// [`Gradient::velocities`] orders them) at the nodes of `model`: CᵀMC, with C the
    /// derivative of ρ, λ = ρ(vp² - 2vs²) and μ = ρ vs² in vp, vs and ρ at each node.
    pub fn in_velocities(&self, model: &Model, matrix: &[f64]) -> Vec<f64> {
        let n = self.nodes.len();
        let side = 3 * n;
        let jacobians: Vec<[[f64; 3]; 3]> = self
            .nodes
            .iter()
            .map(|&i| {
                let (a, b, r) = (model.vp[i] as f64, model.vs[i] as f64, model.rho[i] as f64);
                [
                    [0.0, 0.0, 1.0],
                    [2.0 * r * a, -4.0 * r * b, a * a - 2.0 * b * b],
                    [0.0, 2.0 * r * b, b * b],
                ]
            })
            .collect();
        let mut out = vec![0.0; side * side];
        for (i, out) in out.iter_mut().enumerate() {
            let (row, column) = (i / side, i % side);
            let (v, a, w, b) = (row / n, row % n, column / n, column % n);
            for p in 0..3 {
                for q in 0..3 {
                    *out += jacobians[a][p][v]
                        * matrix[(p * n + a) * side + q * n + b]
                        * jacobians[b][q][w];
                }
            }
        }
        out
    }

    /// [`Block::jacobian`] for unit changes of vp, vs and ρ (the columns in that order, a
    /// parameter's nodes together) at the nodes of `model`: J C.
    pub fn jacobian_in_velocities(&self, model: &Model) -> Vec<Vec<Vec<f32>>> {
        let n = self.nodes.len();
        (0..3 * n)
            .map(|column| {
                let (v, b) = (column / n, column % n);
                let i = self.nodes[b];
                let (vp, vs, rho) = (model.vp[i] as f64, model.vs[i] as f64, model.rho[i] as f64);
                // ∂(ρ, λ, μ)/∂(vp, vs, ρ)'s column.
                let factors = [
                    [0.0, 2.0 * rho * vp, 0.0],
                    [0.0, -4.0 * rho * vs, 2.0 * rho * vs],
                    [1.0, vp * vp - 2.0 * vs * vs, vs * vs],
                ][v];
                let [of_rho, of_lambda, of_mu] =
                    [0, 1, 2].map(|p: usize| &self.jacobian[p * n + b]);
                of_rho
                    .iter()
                    .zip(of_lambda)
                    .zip(of_mu)
                    .map(|((r, l), m)| {
                        r.iter()
                            .zip(l)
                            .zip(m)
                            .map(|((&r, &l), &m)| {
                                (factors[0] * r as f64
                                    + factors[1] * l as f64
                                    + factors[2] * m as f64) as f32
                            })
                            .collect()
                    })
                    .collect()
            })
            .collect()
    }

    /// The Hessian between the nodes in vp, vs and ρ: [`Block::in_velocities`] of
    /// [`Block::hessian`], plus the conversion's own second derivatives applied to the
    /// gradient, three by three at each node, as [`Product::velocities`] adds them.
    pub fn velocities(&self, model: &Model) -> Vec<f64> {
        let n = self.nodes.len();
        let side = 3 * n;
        let mut h = self.in_velocities(model, &self.hessian());
        for (a, &i) in self.nodes.iter().enumerate() {
            let (vp, vs, rho) = (model.vp[i] as f64, model.vs[i] as f64, model.rho[i] as f64);
            let (gl, gm) = (self.gradient.lambda[i], self.gradient.mu[i]);
            let shear = gm - 2.0 * gl;
            let mut add = |v: usize, w: usize, x: f64| h[(v * n + a) * side + w * n + a] += x;
            add(0, 0, 2.0 * rho * gl);
            add(0, 2, 2.0 * vp * gl);
            add(2, 0, 2.0 * vp * gl);
            add(1, 1, 2.0 * rho * shear);
            add(1, 2, 2.0 * vs * shear);
            add(2, 1, 2.0 * vs * shear);
        }
        h
    }
}

/// The staggered positions a change at `nodes` moves, and whose sums reach them through
/// [`to_nodes`]: each node's own and the one before it along every axis and pair of axes
/// (padded indices, sorted).
fn positions_around(g: &Grid, nodes: &[usize]) -> Vec<usize> {
    let [nx, ny, _] = g.n;
    let mut points: Vec<usize> = nodes
        .iter()
        .flat_map(|&m| {
            let at = [m % nx, (m / nx) % ny, m / (nx * ny)];
            (0..8usize).filter_map(move |c| {
                let i = at[0].checked_sub(c & 1)?;
                let j = at[1].checked_sub((c >> 1) & 1)?;
                let k = at[2].checked_sub(c >> 2)?;
                Some(g.idx(i, j, k))
            })
        })
        .collect();
    points.sort_unstable();
    points.dedup();
    points
}

/// What one Born run of [`hessian_block`] leaves.
struct Column {
    /// The size of the change it ran with; the rest is per unit of it.
    scale: f64,
    traces: Vec<Vec<f32>>,
    /// T S(δu, ũ) and δT S(u, ũ) at the nodes: ρ, λ, μ, a parameter's nodes together.
    scattered: Vec<f64>,
    factors: Vec<f64>,
}

/// The Hessian's entries between the nodes `nodes` (interior indices, each once), in ρ, λ and
/// μ, from one forward and one adjoint run of `sim` (fresh, at step 0) and one Born run per
/// node and parameter (the module's text): the gradient run of
/// [`crate::gradient::gradient`] with checkpoints every `checkpoint_every` steps, keeping
/// both wavefields' values at the positions around the nodes for every step (36 bytes per
/// position and step), then the Born runs, `at_once` of them side by side (1 unless the grid
/// is too small to keep the threads busy; each holds a solver of its own). `adjoint_source`
/// and `curvature` are [`hessian_product`]'s; `curvature` is called once per column.
///
/// Exact as the products are, where the gradient is: nodes inside the absorbing layers, or
/// near enough for their positions to lie there, get entries that mean nothing.
pub fn hessian_block(
    sim: &mut Simulation,
    nodes: &[usize],
    adjoint_source: impl FnOnce(&Receivers) -> Result<(f64, Vec<Vec<f32>>)>,
    curvature: impl Fn(&[Vec<f32>]) -> Result<Vec<Vec<f32>>>,
    checkpoint_every: Option<usize>,
    at_once: usize,
) -> Result<Block> {
    let steps = sim.steps;
    let every = checkpoint_interval(steps, checkpoint_every);
    let g = sim.solver.grid;
    let len = g.padded_len();
    let interior = g.interior_len();
    ensure!(
        nodes.iter().all(|&m| m < interior),
        "a node is outside the grid's {interior}"
    );
    {
        let mut sorted = nodes.to_vec();
        sorted.sort_unstable();
        sorted.dedup();
        ensure!(sorted.len() == nodes.len(), "a node is given twice");
    }
    let points = positions_around(&g, nodes);
    let np = points.len();

    // The gradient run, as `gradient_with` replays it from checkpoints.
    let mut checkpoints: Vec<State> = Vec::with_capacity(steps.div_ceil(every));
    for n in 0..steps {
        if n.is_multiple_of(every) {
            checkpoints.push(sim.solver.state());
        }
        sim.advance();
        if (n.is_multiple_of(100) || n + 1 == steps) && !sim.solver.peak_velocity().is_finite() {
            bail!("the forward wavefield blew up by step {}", n + 1);
        }
    }
    sim.finish();
    let (misfit, residuals) = adjoint_source(&sim.receivers)?;
    ensure!(
        misfit.is_finite() && residuals.iter().flatten().all(|r| r.is_finite()),
        "the misfit or its derivative isn't finite; check the observed data"
    );
    let residuals = sim.recorded_residuals(residuals);
    let pressure = PressureResiduals::new(sim, &residuals);
    let mut adjoint = Adjoint::new(sim, residuals)?;
    let six = || std::array::from_fn::<_, 6, _>(|_| vec![0.0f32; len]);
    let (mut strain, mut adjoint_strain) = (six(), six());
    pressure.add(&mut adjoint_strain, steps);
    let mut sums: Sums = vec![[0.0; 8]; len];
    let mut unused = vec![0.0f64; interior];
    // Per step and position: the forward velocities' change and strain rates, and the
    // adjoint velocities and strain, as `correlate` reads them.
    let mut forward_kept = vec![[0.0f32; 9]; steps * np];
    let mut adjoint_kept = vec![[0.0f32; 9]; steps * np];
    let kept = std::mem::size_of_val(forward_kept.as_slice())
        + std::mem::size_of_val(adjoint_kept.as_slice());
    let held: usize = checkpoints.iter().map(State::memory_bytes).sum();
    let fixed = sim.memory_bytes()
        + adjoint.memory_bytes()
        + pressure.memory_bytes()
        + std::mem::size_of_val(sums.as_slice())
        + std::mem::size_of_val(unused.as_slice())
        + 4 * len * 12
        + kept;
    let mut peak_bytes = 0;
    let mut segment: Vec<[Vec<f32>; 3]> = Vec::with_capacity(every + 1);
    for (s, state) in checkpoints.iter().enumerate().rev() {
        let start = s * every;
        let end = (start + every).min(steps);
        sim.solver.restore(state);
        segment.clear();
        segment.push(velocities(&sim.solver));
        for _ in start..end {
            sim.step();
            segment.push(velocities(&sim.solver));
        }
        peak_bytes = peak_bytes.max(fixed + held + 4 * 3 * len * segment.len());
        let forward = &sim.solver;
        for n in (start..end).rev() {
            let [before, after] =
                [n - start, n - start + 1].map(|i| segment[i].each_ref().map(Vec::as_slice));
            forward.strain_rate(after, &mut strain);
            adjoint.step();
            assert_eq!(
                adjoint.forward_step(),
                Some(n),
                "adjoint and forward out of step"
            );
            let back = adjoint.solver.fields.velocities();
            let rows = forward_kept[n * np..(n + 1) * np]
                .iter_mut()
                .zip(&mut adjoint_kept[n * np..(n + 1) * np]);
            for (&p, (f, a)) in points.iter().zip(rows) {
                for c in 0..3 {
                    f[c] = after[c][p] - before[c][p];
                    a[c] = back[c][p];
                }
                for c in 0..6 {
                    f[3 + c] = strain[c][p];
                    a[3 + c] = adjoint_strain[c][p];
                }
            }
            correlate(
                forward,
                &adjoint.solver,
                &mut sums,
                &mut unused,
                [&adjoint_strain, &strain],
                before,
                after,
            );
            let a = &adjoint.solver;
            a.add_strain_rate(a.fields.velocities(), &mut adjoint_strain);
            pressure.add(&mut adjoint_strain, n);
        }
        ensure!(
            adjoint.solver.peak_velocity().is_finite(),
            "the adjoint wavefield blew up"
        );
    }
    drop((
        checkpoints,
        segment,
        adjoint,
        pressure,
        strain,
        adjoint_strain,
        unused,
    ));
    let gradient = to_nodes(&sim.solver, &sums);

    // The Born runs: per column the solver's step on δv, δσ with the kept forward values as
    // its sources, correlated with the kept adjoint values where they were kept.
    let n = nodes.len();
    let side = 3 * n;
    let sim = &*sim;
    let forward = &sim.solver;
    let plane = g.plane_len();
    let run = |column: usize| -> Result<Column> {
        let (parameter, node) = (column / n, nodes[column % n]);
        let model = &forward.model;
        // A change the size of the node's own values, for a Born wavefield about the forward
        // one's size in single precision.
        let (rho, vp) = (model.rho[node] as f64, model.vp[node] as f64);
        let scale = if parameter == 0 { rho } else { rho * vp * vp };
        let mut change = Change {
            rho: vec![0.0; interior],
            lambda: vec![0.0; interior],
            mu: vec![0.0; interior],
        };
        [&mut change.rho, &mut change.lambda, &mut change.mu][parameter][node] = scale;
        let Some(material) = MaterialChange::new(forward, &change).at(&points) else {
            bail!("a node's change moves the material beyond the positions kept");
        };
        let mut solver = forward.at_rest();
        let r = &sim.receivers;
        let mut receivers =
            Receivers::new(&solver, r.points.clone(), r.components.clone(), r.samples);
        let mut born_sums: Sums = vec![[0.0; 8]; len];
        let mut before = vec![[0.0f32; 3]; np];
        for step in 0..steps {
            let kept = &forward_kept[step * np..(step + 1) * np];
            let v = solver.fields.velocities();
            for (b, &p) in before.iter_mut().zip(&points) {
                *b = [v[0][p], v[1][p], v[2][p]];
            }
            solver.step_with(
                |f, _| add_velocity_at(f, &points, &material, kept),
                |f, _| add_stress_at(forward, f, &points, &material, kept),
            );
            if step.is_multiple_of(sim.record_every) {
                receivers.record(&solver);
            }
            let v = solver.fields.velocities();
            let adjoint = &adjoint_kept[step * np..(step + 1) * np];
            for ((&p, a), b) in points.iter().zip(adjoint).zip(&before) {
                let k = p / plane - g.halo;
                correlate_at(
                    forward,
                    k,
                    &mut born_sums[p],
                    a,
                    forward.strain_rate_at(v, p, k),
                    std::array::from_fn(|c| v[c][p] - b[c]),
                );
            }
        }
        ensure!(
            solver.peak_velocity().is_finite(),
            "a Born wavefield blew up"
        );
        let scattered = to_nodes(forward, &born_sums);
        let mut factors = Gradient {
            rho: vec![0.0; interior],
            lambda: vec![0.0; interior],
            mu: vec![0.0; interior],
        };
        add_factors_change(forward, &sums, &change, &mut factors);
        let at_nodes = |g: &Gradient| -> Vec<f64> {
            [&g.rho, &g.lambda, &g.mu]
                .iter()
                .flat_map(|values| nodes.iter().map(|&m| values[m] / scale))
                .collect()
        };
        Ok(Column {
            scale,
            traces: unwarped(sim, &receivers),
            scattered: at_nodes(&scattered),
            factors: at_nodes(&factors),
        })
    };
    let columns: Vec<usize> = (0..side).collect();
    let at_once = at_once.clamp(1, side.max(1));
    let columns: Vec<Column> = if at_once == 1 {
        columns.iter().map(|&c| run(c)).collect::<Result<_>>()?
    } else {
        columns
            .par_chunks(side.div_ceil(at_once))
            .map(|chunk| chunk.iter().map(|&c| run(c)).collect::<Result<Vec<_>>>())
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .flatten()
            .collect()
    };

    let mut block = Block {
        nodes: nodes.to_vec(),
        gradient,
        misfit,
        jacobian: Vec::with_capacity(side),
        gauss_newton: vec![0.0; side * side],
        scattered: vec![0.0; side * side],
        factors: vec![0.0; side * side],
        memory_bytes: 0,
    };
    for (b, column) in columns.iter().enumerate() {
        let curved = curvature(&column.traces)?;
        let through: Vec<f64> = columns
            .par_iter()
            .map(|other| {
                let dot: f64 = other
                    .traces
                    .iter()
                    .flatten()
                    .zip(curved.iter().flatten())
                    .map(|(&x, &y)| x as f64 * y as f64)
                    .sum();
                dot / (other.scale * column.scale)
            })
            .collect();
        let rows = through.iter().zip(&column.scattered).zip(&column.factors);
        for (a, ((&through, &scattered), &factors)) in rows.enumerate() {
            block.gauss_newton[a * side + b] = through;
            block.scattered[a * side + b] = scattered;
            block.factors[a * side + b] = factors;
        }
    }
    let traces: usize = columns
        .iter()
        .flat_map(|c| &c.traces)
        .map(|t| 4 * t.len())
        .sum();
    for column in columns {
        let per_unit = (1.0 / column.scale) as f32;
        let mut traces = column.traces;
        traces.iter_mut().flatten().for_each(|t| *t *= per_unit);
        block.jacobian.push(traces);
    }
    // A Born run under way: its solver, the material's change and the sums, the change and
    // two kernels at the nodes.
    let born = forward.memory_bytes() + 4 * len * 8 + 64 * len + 8 * interior * 9;
    let second = sim.memory_bytes()
        + std::mem::size_of_val(sums.as_slice())
        + kept
        + at_once * born
        + 2 * traces
        + 3 * 8 * side * side;
    block.memory_bytes = peak_bytes.max(second);
    Ok(block)
}
