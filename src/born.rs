//! The Born run: the change of a run's traces for a small change of its model, J δm, from a
//! second solver stepped beside the forward one (the tangent-linear model of the discrete run;
//! Tarantola 1984's Born approximation, Fichtner & Trampert 2011's scattered field δ₁u).
//!
//! With the forward equations written as in the `gradient` module, v⁽ⁿ⁺½⁾ - v⁽ⁿ⁻½⁾ =
//! B (D σ⁽ⁿ⁾ + f⁽ⁿ⁾) and σ⁽ⁿ⁺¹⁾ - σ⁽ⁿ⁾ = M e⁽ⁿ⁺½⁾ + g⁽ⁿ⁾, a change δB, δM of the staggered
//! material changes the fields by δv, δσ with
//!
//! δv⁽ⁿ⁺½⁾ - δv⁽ⁿ⁻½⁾ = B D δσ⁽ⁿ⁾ + (δB / B) (v⁽ⁿ⁺½⁾ - v⁽ⁿ⁻½⁾),
//! δσ⁽ⁿ⁺¹⁾ - δσ⁽ⁿ⁾ = M e(δv⁽ⁿ⁺½⁾) + δM e⁽ⁿ⁺½⁾:
//!
//! the solver's own step on δv and δσ, with the forward run's velocity change and strain rate
//! as sources, which its layers, sponge and free surface then treat as they treat the
//! forward run's. Receivers read δv and δσ as they read the forward fields, and the
//! time-dispersion transform, linear, un-warps them alike. The sources are the gradient's
//! terms (`gradient::correlate`) untransposed, through the same averaging and the same
//! plane-stress moduli on the free surface, so the Born run and the gradient are each
//! other's transposes: ⟨J δm, r⟩ = ⟨δm, Jᵀ r⟩ to single-precision round-off, wherever the
//! gradient is exact (outside the absorbing layers).
//!
//! With the gradient it gives the Gauss–Newton Hessian's products, JᵀWJ δm: the Born traces,
//! weighed by the misfit's Hessian in the samples (`Misfit::gauss_newton_source`), are the
//! adjoint source of a gradient run.

use anyhow::{Result, ensure};
use rayon::prelude::*;

use crate::model::Model;
use crate::receiver::Receivers;
use crate::simulation::Simulation;
use crate::solver::{Fields, Solver};

/// A change of the model at every interior node (x fastest, like [`Model`]), in the
/// parameters the solver's material is made of.
#[derive(Clone, Debug)]
pub struct Change {
    /// δρ, kg/m³.
    pub rho: Vec<f64>,
    /// δλ, Pa.
    pub lambda: Vec<f64>,
    /// δμ, Pa.
    pub mu: Vec<f64>,
}

impl Change {
    /// The change of ρ, λ and μ for changes of vp, vs and ρ (`speeds`, in that order) at the
    /// nodes of `model`, to first order: with λ = ρ(vp² - 2vs²) and μ = ρ vs², δλ = (vp² -
    /// 2vs²) δρ + 2ρ vp δvp - 4ρ vs δvs and δμ = vs² δρ + 2ρ vs δvs. The transpose of
    /// [`crate::gradient::Gradient::velocities`].
    pub fn from_velocities(model: &Model, speeds: [&[f64]; 3]) -> Self {
        let [dvp, dvs, drho] = speeds;
        let n = model.vp.len();
        let (mut rho, mut lambda, mut mu) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
        for i in 0..n {
            let (a, b, r) = (model.vp[i] as f64, model.vs[i] as f64, model.rho[i] as f64);
            rho[i] = drho[i];
            lambda[i] =
                (a * a - 2.0 * b * b) * drho[i] + 2.0 * r * a * dvp[i] - 4.0 * r * b * dvs[i];
            mu[i] = b * b * drho[i] + 2.0 * r * b * dvs[i];
        }
        Self { rho, lambda, mu }
    }

    /// Σ δm · g over the nodes, for a gradient `g` in the same parameters: the change of the
    /// misfit the gradient predicts.
    pub fn dot(&self, g: &crate::gradient::Gradient) -> f64 {
        let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
        dot(&self.rho, &g.rho) + dot(&self.lambda, &g.lambda) + dot(&self.mu, &g.mu)
    }
}

/// The change of the solver's staggered material ([`crate::solver::Material`]), padded like
/// it: δB/B at the velocity positions, dt δ(λ + 2μ) and dt δλ at the nodes (on the free
/// surface the plane-stress moduli's changes instead), and dt δμ̄ at the shear positions.
pub(crate) struct MaterialChange {
    /// δB/B at the vx, vy and vz positions.
    b: [Vec<f32>; 3],
    l2m: Vec<f32>,
    lam: Vec<f32>,
    /// At the σxy, σxz and σyz positions.
    shear: [Vec<f32>; 3],
}

impl MaterialChange {
    /// Forward differentiation of `Material::new` (the free surface's moduli as the stress
    /// update forms them), node by node as `gradient::to_nodes` transposes it.
    pub(crate) fn new(s: &Solver, change: &Change) -> Self {
        let g = s.grid;
        let [nx, ny, nz] = g.n;
        let len = g.padded_len();
        let dt = s.dt;
        let model = &s.model;
        let m = &s.material;
        let mut out = Self {
            b: std::array::from_fn(|_| vec![0.0; len]),
            l2m: vec![0.0; len],
            lam: vec![0.0; len],
            shear: std::array::from_fn(|_| vec![0.0; len]),
        };
        let node = |i: usize, j: usize, k: usize| {
            g.interior_idx(i.min(nx - 1), j.min(ny - 1), k.min(nz - 1))
        };
        let mu_at = |q: usize| model.rho[q] as f64 * (model.vs[q] as f64).powi(2);
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let p = g.idx(i, j, k);
                    let here = node(i, j, k);
                    // B = 2dt/(ρ₁ + ρ₂): δB/B = -(δρ₁ + δρ₂)/(ρ₁ + ρ₂).
                    for c in 0..3 {
                        let mut next = [i, j, k];
                        next[c] += 1;
                        let there = node(next[0], next[1], next[2]);
                        let sum = model.rho[here] as f64 + model.rho[there] as f64;
                        out.b[c][p] = (-(change.rho[here] + change.rho[there]) / sum) as f32;
                    }
                    let (dl2m, dlam) = (
                        dt * (change.lambda[here] + 2.0 * change.mu[here]),
                        dt * change.lambda[here],
                    );
                    if s.free_surface && k == 0 {
                        // a = l2m - lam²/l2m and b = lam - lam²/l2m.
                        let (l2, l) = (m.l2m[p] as f64, m.lam[p] as f64);
                        let r = if l2 > 0.0 { l / l2 } else { 0.0 };
                        out.l2m[p] = ((1.0 + r * r) * dl2m - 2.0 * r * dlam) as f32;
                        out.lam[p] = (r * r * dl2m + (1.0 - 2.0 * r) * dlam) as f32;
                    } else {
                        out.l2m[p] = dl2m as f32;
                        out.lam[p] = dlam as f32;
                    }
                    // dt × the harmonic mean H of μ at four nodes, ∂H/∂μ_q = H²/(4 μ_q²).
                    for (c, axes) in [[0, 1], [0, 2], [1, 2]].into_iter().enumerate() {
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
                        let d: f64 = corners
                            .iter()
                            .zip(&mus)
                            .map(|(&q, &mu)| h * h / (4.0 * mu * mu) * change.mu[q])
                            .sum();
                        out.shear[c][p] = (dt * d) as f32;
                    }
                }
            }
        }
        out
    }

    /// Adds (δB/B) (v⁽ⁿ⁺½⁾ - v⁽ⁿ⁻½⁾) to the Born run's velocities.
    pub(crate) fn add_velocity(&self, f: &mut Fields, before: &[Vec<f32>; 3], after: [&[f32]; 3]) {
        let [vx, vy, vz] = [&mut f.vx, &mut f.vy, &mut f.vz];
        for (c, v) in [vx, vy, vz].into_iter().enumerate() {
            v.par_iter_mut()
                .zip(&self.b[c])
                .zip(after[c].par_iter().zip(&before[c]))
                .for_each(|((v, &b), (&a, &p))| *v += b * (a - p));
        }
    }

    /// Adds δM e to the Born run's stresses, `e` the forward strain rates
    /// ([`Solver::strain_rate`]), or any strain at the stresses' positions. On a free surface
    /// the plane-stress moduli act on εxx and εyy alone, as the stress update's do, and σzz
    /// stays zero.
    pub(crate) fn add_stress(&self, s: &Solver, f: &mut Fields, e: &[Vec<f32>; 6]) {
        let [exx, eyy, ezz, exy, exz, eyz] = e;
        let g = s.grid;
        // The surface plane's padded indices, if there is a free surface.
        let plane = g.plane_len();
        let surface = s.free_surface.then(|| g.halo * plane..(g.halo + 1) * plane);
        f.sxx
            .par_iter_mut()
            .zip(f.syy.par_iter_mut())
            .zip(f.szz.par_iter_mut())
            .enumerate()
            .for_each(|(p, ((xx, yy), zz))| {
                let (l2, l) = (self.l2m[p], self.lam[p]);
                if surface.as_ref().is_some_and(|r| r.contains(&p)) {
                    *xx += l2 * exx[p] + l * eyy[p];
                    *yy += l2 * eyy[p] + l * exx[p];
                } else {
                    *xx += l2 * exx[p] + l * (eyy[p] + ezz[p]);
                    *yy += l2 * eyy[p] + l * (exx[p] + ezz[p]);
                    *zz += l2 * ezz[p] + l * (exx[p] + eyy[p]);
                }
            });

        for (sigma, (m, e)) in [&mut f.sxy, &mut f.sxz, &mut f.syz]
            .into_iter()
            .zip(self.shear.iter().zip([exy, exz, eyz]))
        {
            sigma
                .par_iter_mut()
                .zip(m)
                .zip(e)
                .for_each(|((v, &m), &e)| *v += m * e);
        }
        if s.free_surface {
            let g = s.grid;
            for j in 0..g.n[1] {
                let row = g.idx(0, j, 0);
                f.szz[row..row + g.n[0]].fill(0.0);
            }
        }
    }

    /// The change at the positions `points` (padded indices) alone, for a change of a few
    /// nodes: per position δB/B for vx, vy and vz, the two normal moduli's changes and the
    /// three shear moduli's. `None` if it is not zero everywhere else.
    pub(crate) fn at(&self, points: &[usize]) -> Option<Vec<[f32; 8]>> {
        let arrays = [
            &self.b[0],
            &self.b[1],
            &self.b[2],
            &self.l2m,
            &self.lam,
            &self.shear[0],
            &self.shear[1],
            &self.shear[2],
        ];
        let count = |values: &mut dyn Iterator<Item = f32>| values.filter(|&v| v != 0.0).count();
        let all: usize = arrays.iter().map(|a| count(&mut a.iter().copied())).sum();
        let here: Vec<[f32; 8]> = points.iter().map(|&p| arrays.map(|a| a[p])).collect();
        (count(&mut here.iter().flatten().copied()) == all).then_some(here)
    }

    pub(crate) fn memory_bytes(&self) -> usize {
        4 * self
            .b
            .iter()
            .chain(&self.shear)
            .chain([&self.l2m, &self.lam])
            .map(Vec::len)
            .sum::<usize>()
    }
}

/// [`MaterialChange::add_velocity`] for a change that is zero outside `points`
/// ([`MaterialChange::at`]), `forward` holding, per position, the forward velocities' change
/// over the step first: the same additions, so the same Born run bit for bit.
pub(crate) fn add_velocity_at(
    f: &mut Fields,
    points: &[usize],
    change: &[[f32; 8]],
    forward: &[[f32; 9]],
) {
    for ((&p, m), u) in points.iter().zip(change).zip(forward) {
        f.vx[p] += m[0] * u[0];
        f.vy[p] += m[1] * u[1];
        f.vz[p] += m[2] * u[2];
    }
}

/// [`MaterialChange::add_stress`] likewise, `forward` holding the forward strain rates after
/// the velocities' change.
pub(crate) fn add_stress_at(
    s: &Solver,
    f: &mut Fields,
    points: &[usize],
    change: &[[f32; 8]],
    forward: &[[f32; 9]],
) {
    let g = s.grid;
    let plane = g.plane_len();
    let surface = s.free_surface.then(|| g.halo * plane..(g.halo + 1) * plane);
    for ((&p, m), u) in points.iter().zip(change).zip(forward) {
        let (l2, l) = (m[3], m[4]);
        let [exx, eyy, ezz, exy, exz, eyz] = [u[3], u[4], u[5], u[6], u[7], u[8]];
        if surface.as_ref().is_some_and(|r| r.contains(&p)) {
            f.sxx[p] += l2 * exx + l * eyy;
            f.syy[p] += l2 * eyy + l * exx;
        } else {
            f.sxx[p] += l2 * exx + l * (eyy + ezz);
            f.syy[p] += l2 * eyy + l * (exx + ezz);
            f.szz[p] += l2 * ezz + l * (exx + eyy);
        }
        f.sxy[p] += m[5] * exy;
        f.sxz[p] += m[6] * exz;
        f.syz[p] += m[7] * eyz;
    }
    if s.free_surface {
        for j in 0..g.n[1] {
            let row = g.idx(0, j, 0);
            f.szz[row..row + g.n[0]].fill(0.0);
        }
    }
}

/// A Born run's traces as `receivers` recorded them, un-warped as `sim`'s are once it has
/// finished.
pub(crate) fn unwarped(sim: &Simulation, receivers: &Receivers) -> Vec<Vec<f32>> {
    let samples = receivers.samples;
    let mut traces = receivers.traces.clone();
    for (&c, data) in receivers.components.iter().zip(&mut traces) {
        if let Some(warp) = sim.trace_warp(c).filter(|_| sim.traces_unwarped()) {
            data.par_chunks_mut(samples).for_each(|trace| {
                let x: Vec<f64> = trace.iter().map(|&v| v as f64).collect();
                for (t, y) in trace.iter_mut().zip(warp.apply(&x)) {
                    *t = y as f32;
                }
            });
        }
    }
    traces
}

/// The Born run beside a forward run: the change of the material, the solver holding δv and
/// δσ, and what it records.
pub(crate) struct Scattering {
    pub(crate) material: MaterialChange,
    pub(crate) solver: Solver,
    pub(crate) receivers: Receivers,
    /// The forward velocities before the step under way, and its strain rates after.
    before: [Vec<f32>; 3],
    strain: [Vec<f32>; 6],
}

impl Scattering {
    /// At rest beside `sim`, which must be too, for the model change `change`.
    pub(crate) fn new(sim: &Simulation, change: &Change) -> Result<Self> {
        ensure!(sim.solver.step == 0, "the forward run must start at rest");
        let n = sim.solver.grid.interior_len();
        ensure!(
            change.rho.len() == n && change.lambda.len() == n && change.mu.len() == n,
            "the change has {} values, the grid {n} nodes",
            change.rho.len()
        );
        let solver = sim.solver.at_rest();
        let r = &sim.receivers;
        let receivers = Receivers::new(&solver, r.points.clone(), r.components.clone(), r.samples);
        let len = sim.solver.grid.padded_len();
        Ok(Self {
            material: MaterialChange::new(&sim.solver, change),
            solver,
            receivers,
            before: std::array::from_fn(|_| vec![0.0; len]),
            strain: std::array::from_fn(|_| vec![0.0; len]),
        })
    }

    /// Advances the forward run a step and the Born run with it, recording both if `record`
    /// and the step is one to record ([`Simulation::advance`]); without, as
    /// [`Simulation::step`] replays.
    pub(crate) fn advance(&mut self, sim: &mut Simulation, record: bool) {
        let n = sim.solver.step;
        for (b, v) in self.before.iter_mut().zip(sim.solver.fields.velocities()) {
            b.copy_from_slice(v);
        }
        if record {
            sim.advance();
        } else {
            sim.step();
        }
        let after = sim.solver.fields.velocities();
        sim.solver.strain_rate(after, &mut self.strain);
        let (forward, material, before, strain) =
            (&sim.solver, &self.material, &self.before, &self.strain);
        self.solver.step_with(
            |f, _| material.add_velocity(f, before, after),
            |f, _| material.add_stress(forward, f, strain),
        );
        if record && n.is_multiple_of(sim.record_every) {
            self.receivers.record(&self.solver);
        }
    }

    /// The Born traces, un-warped as `sim`'s are once it has finished.
    pub(crate) fn traces(&mut self, sim: &Simulation) -> Result<Vec<Vec<f32>>> {
        ensure!(
            self.solver.peak_velocity().is_finite(),
            "the Born wavefield blew up"
        );
        Ok(unwarped(sim, &self.receivers))
    }

    /// Bytes it holds.
    pub(crate) fn memory_bytes(&self) -> usize {
        let traces: usize = self.receivers.traces.iter().map(Vec::len).sum();
        let scratch: usize = self.before.iter().chain(&self.strain).map(Vec::len).sum();
        self.solver.memory_bytes() + self.material.memory_bytes() + 4 * (scratch + traces)
    }
}

/// What a Born run found.
pub struct Born {
    /// The change of the traces, in the layout of [`Receivers::traces`], un-warped like the
    /// forward run's.
    pub traces: Vec<Vec<f32>>,
    /// The most memory held at once besides the forward run's own.
    pub memory_bytes: usize,
}

/// Runs `sim` from rest with a Born run beside it: the change of its traces for the model
/// change `change`, to first order. `sim` ends as [`Simulation::run`] leaves it, its traces
/// finished.
pub fn born(sim: &mut Simulation, change: &Change) -> Result<Born> {
    let mut scattering = Scattering::new(sim, change)?;
    for _ in 0..sim.steps {
        scattering.advance(sim, true);
    }
    sim.finish();
    Ok(Born {
        traces: scattering.traces(sim)?,
        memory_bytes: scattering.memory_bytes(),
    })
}
