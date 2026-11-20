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
//! nor the sponge are self-adjoint, so the gradient there is only approximate.

use anyhow::{Result, bail, ensure};
use rayon::prelude::*;

use crate::adjoint::Adjoint;
use crate::model::Model;
use crate::receiver::{Component, Receivers};
use crate::simulation::Simulation;
use crate::solver::{Field, Solver, State};

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

/// Weights 1/‖d_obs,c‖² that make each component's misfit relative to its own energy, so
/// that pressure and velocity count alike. A component that observed nothing gets 0.
pub fn balancing_weights(observed: &[Vec<f32>]) -> Vec<f64> {
    observed
        .iter()
        .map(|o| {
            let energy: f64 = o.iter().map(|&x| (x as f64).powi(2)).sum();
            if energy > 0.0 { 1.0 / energy } else { 0.0 }
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
    /// Steps between checkpoints, and the bytes the checkpoints and one replayed segment took.
    pub checkpoint_every: usize,
    pub checkpoint_bytes: usize,
}

/// The correlations summed over the run, on the padded grid, per position: the density
/// terms Σ ṽ Δv for vx, vy, vz, then the stiffness terms for l2m and lam (on the free surface
/// the plane-stress a and b instead) and for mxy, mxz, myz, already weighted by -W.
type Sums = Vec<[f64; 8]>;

/// Runs `sim` (fresh, at step 0) forward, recording its receivers and keeping a checkpoint
/// every `checkpoint_every` steps (by default √(3N), which about minimises memory); hands the
/// receivers to `adjoint_source`, which returns the misfit and its derivative with respect to
/// every sample (e.g. [`l2_misfit`]); then runs the adjoint backwards, replaying the forward
/// run a segment at a time from the checkpoints, and correlates the two.
///
/// The time step and the absorbing layers stay as `sim` set them up, so the gradient is that
/// of the forward run with them fixed. With dt from a Courant fraction, both depend on the
/// fastest vp in the model; perturbing it changes them too.
pub fn gradient(
    sim: &mut Simulation,
    adjoint_source: impl FnOnce(&Receivers) -> Result<(f64, Vec<Vec<f32>>)>,
    checkpoint_every: Option<usize>,
) -> Result<Run> {
    ensure!(sim.solver.step == 0, "the forward run must start at rest");
    let steps = sim.steps;
    let every = checkpoint_every
        .unwrap_or_else(|| ((3 * steps) as f64).sqrt().ceil() as usize)
        .max(1);

    // Forward, keeping the state at the start of every segment.
    let mut checkpoints: Vec<State> = Vec::with_capacity(steps.div_ceil(every));
    for n in 0..steps {
        if n.is_multiple_of(every) {
            checkpoints.push(sim.solver.state());
        }
        sim.advance();
        if n.is_multiple_of(100) && !sim.solver.peak_velocity().is_finite() {
            bail!("the forward wavefield blew up by step {}", n + 1);
        }
    }
    let (misfit, residuals) = adjoint_source(&sim.receivers)?;
    let pressure = PressureResiduals::new(sim, &residuals);
    let mut adjoint = Adjoint::new(sim, residuals)?;

    let g = sim.solver.grid;
    let len = g.padded_len();
    let mut sums: Sums = vec![[0.0; 8]; len];
    let mut strain = std::array::from_fn::<_, 6, _>(|_| vec![0.0f32; len]);
    let mut adjoint_strain = std::array::from_fn::<_, 6, _>(|_| vec![0.0f32; len]);
    pressure.add(&mut adjoint_strain, steps);

    let mut segment: Vec<[Vec<f32>; 3]> = Vec::with_capacity(every);
    let mut peak_bytes = 0;
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
        let held = checkpoints.iter().map(State::memory_bytes).sum::<usize>()
            + 4 * 3 * len * segment.len();
        peak_bytes = peak_bytes.max(held);

        for n in (start..end).rev() {
            let (before, after) = (&segment[n - start], &segment[n - start + 1]);
            // λ̂σ⁽ⁿ⁺¹⁾ = -W ε̃⁽ⁿ⁺¹⁾ with the forward strain rate at n + ½.
            sim.solver
                .strain_rate([&after[0], &after[1], &after[2]], &mut strain);
            correlate_stiffness(&sim.solver, &mut sums, &adjoint_strain, &strain);
            // ṽ⁽ⁿ⁺½⁾ with the forward velocity change over step n.
            adjoint.step();
            debug_assert_eq!(adjoint.forward_step(), Some(n));
            correlate_density(&adjoint.solver, &mut sums, before, after);
            // ε̃⁽ⁿ⁾ = ε̃⁽ⁿ⁺¹⁾ + e(ṽ⁽ⁿ⁺½⁾) - W⁻¹ r_σ⁽ⁿ⁾.
            let f = &adjoint.solver.fields;
            adjoint
                .solver
                .strain_rate([&f.vx, &f.vy, &f.vz], &mut strain);
            adjoint_strain
                .par_iter_mut()
                .zip(&strain)
                .for_each(|(a, e)| a.iter_mut().zip(e).for_each(|(a, e)| *a += e));
            pressure.add(&mut adjoint_strain, n);
        }
        if !adjoint.solver.peak_velocity().is_finite() {
            bail!("the adjoint wavefield blew up");
        }
    }

    Ok(Run {
        gradient: to_nodes(&sim.solver, &sums),
        misfit,
        checkpoint_every: every,
        checkpoint_bytes: peak_bytes,
    })
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
    /// Padded index of each pressure receiver's node and the weight 1/(3W) there.
    nodes: Vec<(usize, f32)>,
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
        let nodes = r
            .nodes
            .iter()
            .map(|&[i, j, k]| {
                (
                    g.idx(i, j, k),
                    1.0 / (3.0 * volume * s.cell_share(Field::Sxx, k)),
                )
            })
            .collect();
        let values = r
            .components
            .iter()
            .position(|&c| c == Component::P)
            .and_then(|c| residuals.get(c).cloned());
        Self {
            nodes,
            values,
            samples: r.samples,
            record_every: sim.record_every,
        }
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
        for (r, &(p, w)) in self.nodes.iter().enumerate() {
            let delta = values[r * self.samples + sample];
            for e in &mut adjoint_strain[..3] {
                e[p] += w * delta;
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
fn correlate_density(
    adjoint: &Solver,
    sums: &mut Sums,
    before: &[Vec<f32>; 3],
    after: &[Vec<f32>; 3],
) {
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
        let w = balancing_weights(&[vec![3.0, 4.0], vec![0.0; 2]]);
        assert_eq!(w, vec![1.0 / 25.0, 0.0]);
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
