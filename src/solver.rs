//! The time stepper: staggered-grid velocity–stress finite differences, fourth order in space
//! by default (any even order up to 16 with [`Solver::with_stencil`]), second order in time.
//!
//! Where each quantity lives, for the cell whose padded index is that of node (i, j, k):
//!
//! | quantity           | position                  |
//! |--------------------|---------------------------|
//! | σxx, σyy, σzz, λ, μ | (i, j, k)                 |
//! | vx                 | (i + ½, j, k)             |
//! | vy                 | (i, j + ½, k)             |
//! | vz                 | (i, j, k + ½)             |
//! | σxy                | (i + ½, j + ½, k)         |
//! | σxz                | (i + ½, j, k + ½)         |
//! | σyz                | (i, j + ½, k + ½)         |
//!
//! Each step advances velocity from t - dt/2 to t + dt/2 using stress at t, then stress from
//! t to t + dt using the new velocity. The work is split across threads into blocks of rows,
//! each worked through a row (constant j and k) at a time, plane by plane upwards: the nine
//! derivatives a row needs go into scratch rows, the CPML corrects those in a layer, and the
//! fields are updated from them. Every node's arithmetic is the same however the work is
//! split, so the results are too, bit for bit.
//!
//! The free surface sits on the σzz plane k = 0, with stress imaging (Levander 1988; Graves
//! 1996): σzz is zero there, σxx and σyy are updated as if it were (∂z vz = -λ/(λ+2μ)
//! (∂x vx + ∂y vy)), and σzz, σxz, σyz above the surface mirror those below with the opposite
//! sign. Velocities above the surface mirror those below with the same sign. That choice makes
//! the discrete z-derivatives of stress and of velocity negative adjoints of each other next to
//! the surface, as they are everywhere else, so the scheme keeps conserving energy there and
//! stays stable; every stencil stays fourth order.
//!
//! Receivers and sources go through [`PointOp`]s: receivers interpolate staggered quantities
//! onto nodes, and sources spread onto the staggered positions with the same weights, so
//! injecting at a point is exactly the transpose of recording there (what back-propagating
//! adjoint sources needs).

#![allow(clippy::needless_range_loop)]

use std::marker::PhantomData;
use std::ops::Range;

use rayon::prelude::*;

use crate::cpml::{CpmlAxis, PmlSettings, Stagger};
use crate::fd::Stencil;
use crate::grid::{AsPoint, Grid};
use crate::model::Model;
use crate::receiver::Component;
use crate::simd::Kernels;

/// How many values below a free surface estimate a field's slope there, one-sidedly
/// ([`Solver::continue_with_slope`]): fifth order. Four points misread a vertically incident P
/// wave, whose even profile curves sharply, by 2.3% at kh = 0.6; six by 0.06%.
const SLOPE_POINTS: usize = 6;

/// The nine wavefield arrays, padded.
#[derive(Clone, Debug)]
pub struct Fields {
    pub vx: Vec<f32>,
    pub vy: Vec<f32>,
    pub vz: Vec<f32>,
    pub sxx: Vec<f32>,
    pub syy: Vec<f32>,
    pub szz: Vec<f32>,
    pub sxy: Vec<f32>,
    pub sxz: Vec<f32>,
    pub syz: Vec<f32>,
}

/// One of the wavefield arrays.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Field {
    Vx,
    Vy,
    Vz,
    Sxx,
    Syy,
    Szz,
    Sxy,
    Sxz,
    Syz,
}

impl Field {
    /// Every field, in the order of [`Fields::arrays`].
    pub const ALL: [Field; 9] = [
        Field::Vx,
        Field::Vy,
        Field::Vz,
        Field::Sxx,
        Field::Syy,
        Field::Szz,
        Field::Sxy,
        Field::Sxz,
        Field::Syz,
    ];

    /// Whether it is a particle velocity rather than a stress.
    pub fn is_velocity(self) -> bool {
        matches!(self, Field::Vx | Field::Vy | Field::Vz)
    }

    /// Where the field sits along each axis, in nodes past each node (0 or ½), and its sign
    /// mirrored across the free surface: even for the velocities and σxx, σyy and σxy, odd for
    /// σzz, σxz and σyz.
    pub fn lattice(self) -> ([f64; 3], f32) {
        match self {
            Field::Vx => ([0.5, 0.0, 0.0], 1.0),
            Field::Vy => ([0.0, 0.5, 0.0], 1.0),
            Field::Vz => ([0.0, 0.0, 0.5], 1.0),
            Field::Sxx | Field::Syy => ([0.0; 3], 1.0),
            Field::Szz => ([0.0; 3], -1.0),
            Field::Sxy => ([0.5, 0.5, 0.0], 1.0),
            Field::Sxz => ([0.5, 0.0, 0.5], -1.0),
            Field::Syz => ([0.0, 0.5, 0.5], -1.0),
        }
    }
}

impl Fields {
    pub fn zeros(len: usize) -> Self {
        Self {
            vx: vec![0.0; len],
            vy: vec![0.0; len],
            vz: vec![0.0; len],
            sxx: vec![0.0; len],
            syy: vec![0.0; len],
            szz: vec![0.0; len],
            sxy: vec![0.0; len],
            sxz: vec![0.0; len],
            syz: vec![0.0; len],
        }
    }

    pub fn get(&self, f: Field) -> &[f32] {
        match f {
            Field::Vx => &self.vx,
            Field::Vy => &self.vy,
            Field::Vz => &self.vz,
            Field::Sxx => &self.sxx,
            Field::Syy => &self.syy,
            Field::Szz => &self.szz,
            Field::Sxy => &self.sxy,
            Field::Sxz => &self.sxz,
            Field::Syz => &self.syz,
        }
    }

    pub fn get_mut(&mut self, f: Field) -> &mut [f32] {
        match f {
            Field::Vx => &mut self.vx,
            Field::Vy => &mut self.vy,
            Field::Vz => &mut self.vz,
            Field::Sxx => &mut self.sxx,
            Field::Syy => &mut self.syy,
            Field::Szz => &mut self.szz,
            Field::Sxy => &mut self.sxy,
            Field::Sxz => &mut self.sxz,
            Field::Syz => &mut self.syz,
        }
    }

    /// Every array, in the order of [`Field::ALL`].
    pub fn arrays(&self) -> [&Vec<f32>; 9] {
        [
            &self.vx, &self.vy, &self.vz, &self.sxx, &self.syy, &self.szz, &self.sxy, &self.sxz,
            &self.syz,
        ]
    }

    /// Every array, mutably, in the order of [`Field::ALL`].
    pub fn arrays_mut(&mut self) -> [&mut Vec<f32>; 9] {
        [
            &mut self.vx,
            &mut self.vy,
            &mut self.vz,
            &mut self.sxx,
            &mut self.syy,
            &mut self.szz,
            &mut self.sxy,
            &mut self.sxz,
            &mut self.syz,
        ]
    }

    /// The particle velocities vx, vy and vz.
    pub fn velocities(&self) -> [&[f32]; 3] {
        [&self.vx, &self.vy, &self.vz]
    }

    /// Bytes its arrays hold.
    pub fn memory_bytes(&self) -> usize {
        4 * self.arrays().iter().map(|v| v.len()).sum::<usize>()
    }
}

/// Material properties where the update equations need them, multiplied by dt and padded.
#[derive(Clone, Debug)]
pub struct Material {
    /// dt/ρ at the vx, vy and vz positions (arithmetic mean of ρ).
    pub bx: Vec<f32>,
    pub by: Vec<f32>,
    pub bz: Vec<f32>,
    /// dt·λ and dt·(λ + 2μ) at the nodes.
    pub lam: Vec<f32>,
    pub l2m: Vec<f32>,
    /// dt·μ at the σxy, σxz and σyz positions (harmonic mean of μ).
    pub mxy: Vec<f32>,
    pub mxz: Vec<f32>,
    pub myz: Vec<f32>,
}

impl Material {
    pub fn new(model: &Model, dt: f64) -> Self {
        let g = model.grid;
        let [nx, ny, nz] = g.n;
        let len = g.padded_len();
        let mut m = Self {
            bx: vec![0.0; len],
            by: vec![0.0; len],
            bz: vec![0.0; len],
            lam: vec![0.0; len],
            l2m: vec![0.0; len],
            mxy: vec![0.0; len],
            mxz: vec![0.0; len],
            myz: vec![0.0; len],
        };
        let mu: Vec<f64> = model
            .vs
            .iter()
            .zip(&model.rho)
            .map(|(&vs, &rho)| rho as f64 * vs as f64 * vs as f64)
            .collect();
        // Neighbours past the last node take its values.
        let at = |v: &[f64], i: usize, j: usize, k: usize| {
            v[g.interior_idx(i.min(nx - 1), j.min(ny - 1), k.min(nz - 1))]
        };
        let rho: Vec<f64> = model.rho.iter().map(|&r| r as f64).collect();
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let p = g.idx(i, j, k);
                    let n = g.interior_idx(i, j, k);
                    let vp = model.vp[n] as f64;
                    let l2m = rho[n] * vp * vp;
                    m.l2m[p] = (dt * l2m) as f32;
                    m.lam[p] = (dt * (l2m - 2.0 * mu[n])) as f32;
                    let r = rho[n];
                    m.bx[p] = (dt * 2.0 / (r + at(&rho, i + 1, j, k))) as f32;
                    m.by[p] = (dt * 2.0 / (r + at(&rho, i, j + 1, k))) as f32;
                    m.bz[p] = (dt * 2.0 / (r + at(&rho, i, j, k + 1))) as f32;
                    let shear = |a: [usize; 3], b: [usize; 3]| {
                        let mus = [
                            at(&mu, i, j, k),
                            at(&mu, i + a[0], j + a[1], k + a[2]),
                            at(&mu, i + b[0], j + b[1], k + b[2]),
                            at(&mu, i + a[0] + b[0], j + a[1] + b[1], k + a[2] + b[2]),
                        ];
                        (dt * harmonic_mean(mus)) as f32
                    };
                    m.mxy[p] = shear([1, 0, 0], [0, 1, 0]);
                    m.mxz[p] = shear([1, 0, 0], [0, 0, 1]);
                    m.myz[p] = shear([0, 1, 0], [0, 0, 1]);
                }
            }
        }
        m
    }

    /// dt/ρ at the position of velocity component `axis` (0 = x).
    pub fn buoyancy(&self, axis: usize) -> &[f32] {
        [&self.bx, &self.by, &self.bz][axis]
    }
}

/// Zero if any value is zero, so shear stress vanishes next to a fluid.
fn harmonic_mean(v: [f64; 4]) -> f64 {
    if v.iter().any(|&x| x <= 0.0) {
        0.0
    } else {
        4.0 / v.iter().map(|x| 1.0 / x).sum::<f64>()
    }
}

/// Everything that changes as the solver steps: the fields, the CPML memory and the step
/// count. Restoring one replays the run from there bit for bit, since the updates don't
/// depend on how the work is split between threads.
#[derive(Clone, Debug)]
pub struct State {
    fields: Fields,
    psi_v: [Vec<f32>; 3],
    psi_s: [Vec<f32>; 3],
    step: usize,
}

impl State {
    /// Bytes it holds.
    pub fn memory_bytes(&self) -> usize {
        let psi: usize = self.psi_v.iter().chain(&self.psi_s).map(Vec::len).sum();
        self.fields.memory_bytes() + 4 * psi
    }
}

/// The wavefield and everything needed to advance it.
pub struct Solver {
    pub grid: Grid,
    pub dt: f64,
    pub free_surface: bool,
    pub model: Model,
    pub material: Material,
    pub fields: Fields,
    /// The layers along x, y and z.
    pub cpml: [CpmlAxis; 3],
    /// CPML memory per axis, for the derivatives of stress (velocity update) and of
    /// velocity (stress update). See [`Outputs::row`] for the layout.
    psi_v: [Vec<f32>; 3],
    psi_s: [Vec<f32>; 3],
    /// The spatial derivative operator.
    pub stencil: Stencil,
    /// Factors along x, y and z for the plain damping inside the side layers, if any.
    sponge: Option<[Vec<f32>; 3]>,
    /// Whether point operators continue evenly mirrored fields above a free surface with their
    /// slope ([`Solver::continue_with_slope`]); off, they read them as mirrored.
    pub surface_slope: bool,
    /// Steps taken.
    pub step: usize,
    /// The compiled update loops it runs, which this CPU can ([`Solver::set_kernels`]).
    kernels: Kernels,
    /// How the updates split the nodes between tasks ([`blocks`]).
    blocks: Vec<Block>,
}

impl Solver {
    /// A solver at rest. Layers `pml_width` nodes thick line every side but the top when
    /// there is a free surface.
    pub fn new(
        model: Model,
        dt: f64,
        pml_width: usize,
        free_surface: bool,
        pml: &PmlSettings,
    ) -> Self {
        let g = model.grid;
        let [nx, ny, nz] = g.n;
        let cpml: [CpmlAxis; 3] = std::array::from_fn(|a| {
            let lo = if a == 2 && free_surface { 0 } else { pml_width };
            CpmlAxis::new(g.n[a], lo, pml_width, g.h[a], dt, pml)
        });
        // Only the side layers need it, and only along a free surface.
        let sponge = (free_surface && pml.sponge > 0.0 && pml_width > 0).then(|| {
            [
                cpml[0].sponge_factors(pml.sponge, dt),
                cpml[1].sponge_factors(pml.sponge, dt),
                vec![1.0; nz],
            ]
        });
        let sizes = [
            nz * ny * 3 * cpml[0].width(),
            nz * cpml[1].width() * 3 * nx,
            cpml[2].width() * ny * 3 * nx,
        ];
        Self {
            grid: g,
            dt,
            free_surface,
            material: Material::new(&model, dt),
            model,
            fields: Fields::zeros(g.padded_len()),
            cpml,
            psi_v: sizes.map(|s| vec![0.0; s]),
            psi_s: sizes.map(|s| vec![0.0; s]),
            sponge,
            stencil: Stencil::taylor(4).expect("order 4 exists"),
            surface_slope: false,
            step: 0,
            kernels: Kernels::chosen(),
            blocks: blocks(&g, block_size(&g)),
        }
    }

    /// The compiled update loops it runs (crate::simd).
    pub fn kernels(&self) -> Kernels {
        self.kernels
    }

    /// Runs the update loops compiled as `kernels` say where this CPU can, else the plain
    /// ones; returns those it will run. Both give the same results bit for bit.
    pub fn set_kernels(&mut self, kernels: Kernels) -> Kernels {
        // Only a CPU that has them may run the AVX2 loops: `sweep` relies on it.
        self.kernels = kernels.or_plain();
        self.kernels
    }

    /// Reads and spreads fields mirrored evenly above a free surface continued with their
    /// slope there ([`Solver::continue_with_slope`]). Set it before making sources and
    /// receivers: their weights come from it.
    pub fn with_surface_slope(mut self, on: bool) -> Self {
        self.surface_slope = on;
        self
    }

    /// Uses `stencil` for the spatial derivatives. The grid's halo must be at least its
    /// half-length; choose the time step from its [`Stencil::sum_abs`].
    pub fn with_stencil(mut self, stencil: Stencil) -> Self {
        assert!(
            stencil.half_length() <= self.grid.halo,
            "a stencil of order {} needs a halo of {}, not {}",
            stencil.order(),
            stencil.half_length(),
            self.grid.halo
        );
        self.stencil = stencil;
        self
    }

    /// The same solver at rest: the model, time step, layers and stencil, with every field and
    /// every CPML memory variable zero and no steps taken.
    pub fn at_rest(&self) -> Self {
        Self {
            grid: self.grid,
            dt: self.dt,
            free_surface: self.free_surface,
            model: self.model.clone(),
            material: self.material.clone(),
            fields: Fields::zeros(self.grid.padded_len()),
            cpml: self.cpml.clone(),
            psi_v: self.psi_v.each_ref().map(|p| vec![0.0; p.len()]),
            psi_s: self.psi_s.each_ref().map(|p| vec![0.0; p.len()]),
            stencil: self.stencil.clone(),
            sponge: self.sponge.clone(),
            surface_slope: self.surface_slope,
            step: 0,
            kernels: self.kernels,
            blocks: self.blocks.clone(),
        }
    }

    /// What it would take to resume from here.
    pub fn state(&self) -> State {
        State {
            fields: self.fields.clone(),
            psi_v: self.psi_v.clone(),
            psi_s: self.psi_s.clone(),
            step: self.step,
        }
    }

    /// Bytes a [`State`] of it holds, without making one.
    pub fn state_bytes(&self) -> usize {
        let psi: usize = self.psi_v.iter().chain(&self.psi_s).map(Vec::len).sum();
        4 * (9 * self.grid.padded_len() + psi)
    }

    /// Goes back (or forward) to a state this solver, or one [`Solver::at_rest`] from it,
    /// was in.
    pub fn restore(&mut self, state: &State) {
        self.fields.clone_from(&state.fields);
        self.psi_v.clone_from(&state.psi_v);
        self.psi_s.clone_from(&state.psi_s);
        self.step = state.step;
    }

    /// Bytes the solver holds.
    pub fn memory_bytes(&self) -> usize {
        let m = &self.material;
        let material = [&m.bx, &m.by, &m.bz, &m.lam, &m.l2m, &m.mxy, &m.mxz, &m.myz];
        let psi = self.psi_v.iter().chain(&self.psi_s);
        let sponge = self.sponge.iter().flatten();
        let values: usize = material
            .into_iter()
            .chain(psi)
            .chain(sponge)
            .map(Vec::len)
            .sum();
        self.fields.memory_bytes() + self.model.memory_bytes() + 4 * values
    }

    /// Whether the node nearest `point` (a node, or grid coordinates) is in an absorbing layer.
    pub fn in_layers(&self, point: impl AsPoint) -> bool {
        point
            .point()
            .iter()
            .zip(&self.cpml)
            .any(|(&x, c)| c.slab(x.round() as usize).is_some())
    }

    /// Time of the velocities now held: steps × dt - dt/2.
    pub fn velocity_time(&self) -> f64 {
        (self.step as f64 - 0.5) * self.dt
    }

    /// Advances one step. `force` and `moment` add the sources' contributions: forces to
    /// velocity at time step·dt, moment tensors to stress at (step + ½)·dt.
    pub fn step_with(
        &mut self,
        force: impl FnOnce(&mut Fields, f64),
        moment: impl FnOnce(&mut Fields, f64),
    ) {
        let t = self.step as f64 * self.dt;
        self.update_velocity(1.0, true);
        force(&mut self.fields, t);
        if self.free_surface {
            self.mirror_velocity();
        }
        self.update_stress(1.0, true);
        moment(&mut self.fields, t + 0.5 * self.dt);
        self.damp_layers();
        if self.free_surface {
            self.mirror_stress();
        }
        self.step += 1;
    }

    /// Advances one step without sources.
    pub fn step(&mut self) {
        self.step_with(|_, _| {}, |_, _| {});
    }

    /// Takes back the last step inside the layers' [`Band`], from `before`, the band's values
    /// before that step ([`Band::save`]): the reverse of [`Solver::step_with`], whose `force`
    /// and `moment` must withdraw what the step's sources added. Outside the band the result
    /// means nothing; inside the region it encloses, it is the state before the step up to
    /// single-precision round-off, so the forward wavefield can be rebuilt backwards from the
    /// final state (Dussaud et al. 2008). The CPML memory is left as it was, so restore a
    /// [`State`] before stepping forward again.
    ///
    /// Every operation of the step is undone in reverse order with the same arithmetic: the
    /// increments are computed from the same values as forward and subtracted, the sponge and
    /// the CPML act only in the layers, σzz on the free surface was zero before the step as
    /// after it, and the mirrors are recomputed.
    pub fn step_back_with(
        &mut self,
        band: &Band,
        before: &[f32],
        force: impl FnOnce(&mut Fields, f64),
        moment: impl FnOnce(&mut Fields, f64),
    ) {
        assert!(self.step > 0, "no step to take back");
        let n = self.step - 1;
        let t = n as f64 * self.dt;
        moment(&mut self.fields, t + 0.5 * self.dt);
        self.update_stress(-1.0, false);
        band.restore_stress(&mut self.fields, before);
        if self.free_surface {
            self.mirror_stress();
        }
        force(&mut self.fields, t);
        self.update_velocity(-1.0, false);
        band.restore_velocity(&mut self.fields, before);
        if self.free_surface {
            self.mirror_velocity();
        }
        self.step = n;
    }

    /// Advances the stresses alone, as the second half of [`Solver::step_with`] does, with
    /// `moment` adding to them, without counting a step: the first half step of a run that
    /// starts with stress and velocity at the same time (the Kirchhoff–Helmholtz rebuild's,
    /// [`crate::surface`]).
    pub fn advance_stress_with(&mut self, moment: impl FnOnce(&mut Fields)) {
        self.update_stress(1.0, true);
        moment(&mut self.fields);
        if self.free_surface {
            self.mirror_stress();
        }
    }

    /// Clears the absorbing layers' memory, as for a run starting from the fields it holds.
    pub fn clear_layers(&mut self) {
        for p in self.psi_v.iter_mut().chain(&mut self.psi_s) {
            p.fill(0.0);
        }
    }

    /// The nodes outside the absorbing layers along x, y and z, where the update is the same
    /// at every node and can be taken back.
    pub fn interior(&self) -> [std::ops::Range<usize>; 3] {
        self.cpml.each_ref().map(CpmlAxis::interior)
    }

    /// Adds `sign` times the velocity increment; with `layers` false, without the CPML.
    fn update_velocity(&mut self, sign: f32, layers: bool) {
        let g = self.grid;
        let Fields {
            vx,
            vy,
            vz,
            sxx,
            syy,
            szz,
            sxy,
            sxz,
            syz,
        } = &mut self.fields;
        let update = VelocityUpdate {
            g,
            sign,
            ih: g.h.map(|h| (1.0 / h) as f32),
            stress: [&sxx[..], &syy[..], &szz[..], &sxy[..], &sxz[..], &syz[..]],
            m: &self.material,
            st: &self.stencil,
            cpml: &self.cpml,
        };
        let mut none = [Vec::new(), Vec::new(), Vec::new()];
        let psi = if layers { &mut self.psi_v } else { &mut none };
        let outputs = Outputs::new([vx, vy, vz], psi);
        sweep(self.kernels, &update, &g, &self.cpml, &self.blocks, outputs);
    }

    /// Adds `sign` times the stress increment; with `layers` false, without the CPML.
    fn update_stress(&mut self, sign: f32, layers: bool) {
        let g = self.grid;
        let Fields {
            vx,
            vy,
            vz,
            sxx,
            syy,
            szz,
            sxy,
            sxz,
            syz,
        } = &mut self.fields;
        let update = StressUpdate {
            g,
            sign,
            ih: g.h.map(|h| (1.0 / h) as f32),
            free_surface: self.free_surface,
            v: [&vx[..], &vy[..], &vz[..]],
            m: &self.material,
            st: &self.stencil,
            cpml: &self.cpml,
        };
        let mut none = [Vec::new(), Vec::new(), Vec::new()];
        let psi = if layers { &mut self.psi_s } else { &mut none };
        let outputs = Outputs::new([sxx, syy, szz, sxy, sxz, syz], psi);
        sweep(self.kernels, &update, &g, &self.cpml, &self.blocks, outputs);
    }

    /// The strain rates the stress update applies the stiffness to, from the velocities
    /// `v` (vx, vy, vz, padded and mirrored as the solver leaves them), into `out` (padded,
    /// interior only): ∂x vx, ∂y vy, ∂z vz, ∂y vx + ∂x vy, ∂z vx + ∂x vz and ∂z vy + ∂y vz,
    /// each at the position of the stress it changes. On the free surface ∂z vz is left at
    /// zero, since the update replaces it there. Without the CPML: inside the absorbing
    /// layers these are not what the update uses.
    pub fn strain_rate(&self, v: [&[f32]; 3], out: &mut [Vec<f32>; 6]) {
        self.strain_rate_into(v, out, false);
    }

    /// Adds [`Solver::strain_rate`] to `sum`: the same values added as if made first, in one
    /// pass instead of two.
    pub fn add_strain_rate(&self, v: [&[f32]; 3], sum: &mut [Vec<f32>; 6]) {
        self.strain_rate_into(v, sum, true);
    }

    fn strain_rate_into(&self, v: [&[f32]; 3], out: &mut [Vec<f32>; 6], add: bool) {
        let g = self.grid;
        let update = StrainRate {
            g,
            ih: g.h.map(|h| (1.0 / h) as f32),
            free_surface: self.free_surface,
            v,
            st: &self.stencil,
            add,
        };
        let [exx, eyy, ezz, exy, exz, eyz] = out;
        let mut none = [Vec::new(), Vec::new(), Vec::new()];
        let outputs = Outputs::new([exx, eyy, ezz, exy, exz, eyz], &mut none);
        sweep(self.kernels, &update, &g, &self.cpml, &self.blocks, outputs);
    }

    /// Damps every field inside the side layers by the sponge factors.
    fn damp_layers(&mut self) {
        let Some(factors) = &self.sponge else {
            return;
        };
        let g = self.grid;
        let halo = g.halo;
        let [nx, ny, nz] = g.n;
        let (len, nxt) = (g.plane_len(), g.padded()[0]);
        let [fx, fy, fz] = factors;
        let x_ranges = self.cpml[0].ranges();
        for v in self.fields.arrays_mut() {
            v[g.halo * len..(g.halo + nz) * len]
                .par_chunks_mut(len)
                .enumerate()
                .for_each(|(kz, plane)| {
                    for j in 0..ny {
                        let o = (j + halo) * nxt + halo;
                        let row = &mut plane[o..o + nx];
                        let yz = fy[j] * fz[kz];
                        if yz < 1.0 {
                            for (v, &x) in row.iter_mut().zip(fx) {
                                *v *= yz * x;
                            }
                        } else {
                            for r in x_ranges.clone() {
                                for i in r {
                                    row[i] *= fx[i];
                                }
                            }
                        }
                    }
                });
        }
    }

    /// The share of a full cell that a value on interior plane `k` stands for: half on the
    /// free surface for the quantities that sit on it (all but vz, σxz and σyz), whose other
    /// half is their mirror image above it. A source there acts on half the mass, so it is
    /// divided by this.
    pub fn cell_share(&self, field: Field, k: usize) -> f32 {
        let on_surface = !matches!(field, Field::Vz | Field::Sxz | Field::Syz);
        if self.free_surface && k == 0 && on_surface {
            0.5
        } else {
            1.0
        }
    }

    /// Mirrors vx, vy and vz across the free surface with the same sign.
    fn mirror_velocity(&mut self) {
        let len = self.grid.plane_len();
        let s = self.grid.halo; // the surface plane, and how many planes the stencils reach above
        let f = &mut self.fields;
        for m in 1..=s {
            // vx and vy sit on their plane: s - m mirrors s + m.
            mirror_plane(&mut f.vx, len, s - m, s + m, 1.0);
            mirror_plane(&mut f.vy, len, s - m, s + m, 1.0);
            // vz sits half a node below its plane: s - m (at s - m + ½) mirrors s + m - 1.
            mirror_plane(&mut f.vz, len, s - m, s + m - 1, 1.0);
        }
    }

    /// Mirrors σzz, σxz and σyz across the free surface with the opposite sign.
    pub(crate) fn mirror_stress(&mut self) {
        let len = self.grid.plane_len();
        let s = self.grid.halo;
        let f = &mut self.fields;
        for m in 1..=s {
            mirror_plane(&mut f.szz, len, s - m, s + m, -1.0);
            mirror_plane(&mut f.sxz, len, s - m, s + m - 1, -1.0);
            mirror_plane(&mut f.syz, len, s - m, s + m - 1, -1.0);
        }
    }

    /// Weights on the interior indices along `axis` that interpolate a quantity stored half a
    /// node forward of each node (at m + ½ for index m) onto node `at`: the stencil's own
    /// interpolation from the 2L values around the node, L its half-length (for Taylor
    /// coefficients Lagrange's, the cubic (-1, 9, 9, -1)/16 at order 4), Lagrange from fewer
    /// near the ends, and a linear extrapolation on the first node. Along z with a free surface, values above it are its mirror images, with
    /// sign `parity` (+1 for vz, -1 for σxz and σyz), so the interpolation reaches up to the
    /// surface.
    pub fn stagger_weights(&self, axis: usize, at: usize, parity: f32) -> Vec<(usize, f32)> {
        let n = self.grid.n[axis];
        let mirror = axis == 2 && self.free_surface;
        // Indices at - m to at + m - 1 must exist, or be mirror images.
        let room = if mirror { usize::MAX } else { at };
        let m = self.stencil.half_length().min(n - at).min(room);
        if m == 0 {
            return self.fold(axis, [(0, 1.5), (1, -0.5)], 0.5, parity);
        }
        let w = if m == self.stencil.half_length() {
            self.stencil.interpolation().to_vec()
        } else {
            crate::fd::midpoint_weights(m)
        };
        let (at, m) = (at as isize, m as isize);
        // Symmetric about the node, and rounded to single precision before folding, so the
        // slope's share comes from the weights as they are applied.
        let weights = (-m..m).map(|o| {
            let w = w[if o < 0 { -o - 1 } else { o } as usize];
            (at + o, w as f32 as f64)
        });
        self.fold(axis, weights, 0.5, parity)
    }

    /// Interpolation weights on indices m along `axis` of a quantity stored at m + `shift`,
    /// as weights on the interior indices: one above a free surface (m < 0) folds onto its
    /// mirror image with sign `parity` (a node's index -m mirrors m; a value half a node down,
    /// at index m, position m + ½, mirrors index -1 - m), and evenly mirrored fields are
    /// continued with their slope ([`Solver::continue_with_slope`]). Zero weights are left out.
    fn fold(
        &self,
        axis: usize,
        weights: impl IntoIterator<Item = (isize, f64)>,
        shift: f64,
        parity: f32,
    ) -> Vec<(usize, f32)> {
        let mut out = Vec::new();
        let mut folded = 0.0;
        for (m, w) in weights {
            let (m, w) = if m < 0 {
                folded += w * -(m as f64 + shift);
                (if shift == 0.0 { -m } else { -1 - m }, parity * w as f32)
            } else {
                (m, w as f32)
            };
            accumulate(&mut out, m as usize, w);
        }
        if parity > 0.0 {
            self.continue_with_slope(&mut out, axis, folded, shift);
        }
        out.retain(|&(_, w)| w != 0.0);
        out
    }

    /// Values folded onto their mirror images across a free surface, with even parity, read
    /// the field as if it had no slope there; but velocities and the stresses that don't
    /// vanish on the surface do have one (σzz = 0 gives ∂z vz = -λ/(λ+2μ) (∂x vx + ∂y vy)),
    /// so evenly mirrored they have a kink at the surface, and an interpolation across it
    /// reads them as if from the mean depth of its folded weights below: vz on the surface
    /// from 0.34 h below at order 8, a Rayleigh wave 11% too strongly per unit kh at each end
    /// (the research log's 2026-11-23 entry). The true continuation is v(-ζ) = v(ζ) - 2ζ ∂z v(0)
    /// to second order, so each weight w folded from -ζ adds -2 w ζ times the slope at the
    /// surface, `folded` (Σ w ζ, in nodes) in all, the slope estimated one-sidedly from the
    /// field's first values below (at m + `shift`): linear in the fields and independent of
    /// the model, so a source, the transpose, stays exact, and so does the gradient. Only with
    /// [`Solver::surface_slope`]: the estimate reaches six nodes down, so an interface there
    /// spoils it, and it misreads waves whose vertical wavenumber passes about 1.1 per node.
    fn continue_with_slope(
        &self,
        out: &mut Vec<(usize, f32)>,
        axis: usize,
        folded: f64,
        shift: f64,
    ) {
        if !self.surface_slope || folded == 0.0 || axis != 2 {
            return;
        }
        let points: Vec<f64> = (0..SLOPE_POINTS.min(self.grid.n[2]))
            .map(|m| m as f64 + shift)
            .collect();
        for (m, d) in crate::fd::lagrange_slope_at_zero(&points)
            .into_iter()
            .enumerate()
        {
            accumulate(out, m, (-2.0 * folded * d) as f32);
        }
    }

    /// Weights on the interior indices along `axis` that interpolate a quantity stored at
    /// m + `shift` (0 or ½) for index m onto grid coordinate `x`: the value itself where `x`
    /// is one of its positions, [`Solver::stagger_weights`] onto a node, and otherwise the
    /// stencil's interpolation from the 2L values around `x` ([`Stencil::interpolation_at`]):
    /// Lagrange's from 2L values moved inwards near an end, and along z with a free surface
    /// the mirror images above it, with sign `parity`.
    pub fn axis_weights(&self, axis: usize, x: f64, shift: f64, parity: f32) -> Vec<(usize, f32)> {
        let n = self.grid.n[axis] as isize;
        let u = x - shift;
        if u >= 0.0 && u.fract() == 0.0 {
            return vec![(u as usize, 1.0)];
        }
        if shift == 0.5 && x.fract() == 0.0 {
            return self.stagger_weights(axis, x as usize, parity);
        }
        let mirror = axis == 2 && self.free_surface;
        let l = self.stencil.half_length() as isize;
        let m0 = u.floor() as isize;
        let (mut lo, mut hi) = (m0 - l + 1, m0 + l);
        let weights: Vec<f64> = if (mirror || lo >= 0) && hi < n {
            self.stencil.interpolation_at(u - m0 as f64)
        } else {
            if hi >= n {
                hi = n - 1;
                lo = hi - 2 * l + 1;
            }
            if lo < 0 && !mirror {
                lo = 0;
                hi = (2 * l - 1).min(n - 1);
            }
            let points: Vec<f64> = (lo..=hi).map(|m| m as f64 + shift - x).collect();
            crate::fd::lagrange_at_zero(&points)
        };
        self.fold(axis, (lo..=hi).zip(weights), shift, parity)
    }

    /// Padded indices and weights that interpolate `field` onto `point`: the product of
    /// [`Solver::axis_weights`] along the three axes, x outermost.
    pub fn field_weights(&self, field: Field, point: impl AsPoint) -> Vec<(usize, f32)> {
        let p = point.point();
        let (shift, parity) = field.lattice();
        let [ax, ay, az]: [Vec<(usize, f32)>; 3] = std::array::from_fn(|a| {
            self.axis_weights(a, p[a], shift[a], if a == 2 { parity } else { 1.0 })
        });
        let g = self.grid;
        let mut out = Vec::with_capacity(ax.len() * ay.len() * az.len());
        for &(i, wx) in &ax {
            for &(j, wy) in &ay {
                for &(k, wz) in &az {
                    out.push((g.idx(i, j, k), wx * wy * wz));
                }
            }
        }
        out
    }

    /// The linear functional that reads component `c` at `point` (a node, or any point in grid
    /// coordinates).
    pub fn receiver_op(&self, c: Component, point: impl AsPoint) -> PointOp {
        let p = point.point();
        // Pressure is -(σxx + σyy + σzz)/3.
        let scale = if c == Component::P { -1.0 / 3.0 } else { 1.0 };
        PointOp {
            terms: c
                .fields()
                .iter()
                .flat_map(|&f| {
                    self.field_weights(f, p)
                        .into_iter()
                        .map(move |(i, w)| (f, i, w * scale))
                })
                .collect(),
        }
    }

    /// How a residual in component `c` at `node` enters the adjoint wavefield: the transpose of
    /// [`Solver::receiver_op`], times the material the update applies there (dt/ρ for a
    /// velocity, dt·C for a stress), divided by the cell's volume and share as a source is.
    ///
    /// A velocity residual is then a unit force: the same operator as a force source of 1 N.
    /// A pressure residual is an isotropic stress source dt·K/V per unit, K = λ + 2μ/3; on the
    /// free surface, where σzz = 0, it loads σxx and σyy with the plane-stress modulus instead.
    /// See the `adjoint` module for why.
    pub fn adjoint_op(&self, c: Component, point: impl AsPoint) -> PointOp {
        let g = self.grid;
        let volume = g.cell_volume() as f32;
        let plane = g.plane_len();
        let m = &self.material;
        let point = point.point();
        let mut op = self.receiver_op(c, point);
        match c {
            Component::Vx | Component::Vy | Component::Vz => {
                let b = m.buoyancy(match c {
                    Component::Vx => 0,
                    Component::Vy => 1,
                    _ => 2,
                });
                for (f, i, w) in &mut op.terms {
                    *w *= b[*i] / (volume * self.cell_share(*f, *i / plane - g.halo));
                }
            }
            Component::P => {
                // The weights the receiver reads σxx, σyy and σzz with at each node: the same
                // but where σzz's odd mirror image folds in near the free surface.
                let mut nodes: Vec<(usize, [f32; 3])> = Vec::new();
                for (c, &field) in Component::P.fields().iter().enumerate() {
                    for (p, w) in self.field_weights(field, point) {
                        match nodes.iter_mut().find(|(q, _)| *q == p) {
                            Some((_, ws)) => ws[c] = w,
                            None => {
                                let mut ws = [0.0; 3];
                                ws[c] = w;
                                nodes.push((p, ws));
                            }
                        }
                    }
                }
                op.terms.clear();
                for (p, w) in nodes {
                    let k = p / plane - g.halo;
                    let (l, l2) = (m.lam[p], m.l2m[p]);
                    let surface = self.free_surface && k == 0;
                    // The update's stiffness applied to the weights, as in update_stress: on
                    // the free surface plane stress (σzz stays zero). Equal weights (on a node,
                    // or away from the surface) take the row sums.
                    let q = if surface && l2 > 0.0 { l * l / l2 } else { 0.0 };
                    let loads: [f32; 3] = if w[0] == w[1] && w[1] == w[2] {
                        let (horizontal, vertical) = if surface {
                            ((l2 - q) + (l - q), 0.0)
                        } else {
                            (l2 + 2.0 * l, l2 + 2.0 * l)
                        };
                        [w[0] * horizontal, w[0] * horizontal, w[0] * vertical]
                    } else if surface {
                        [
                            (l2 - q) * w[0] + (l - q) * w[1],
                            (l - q) * w[0] + (l2 - q) * w[1],
                            0.0,
                        ]
                    } else {
                        [
                            l2 * w[0] + l * w[1] + l * w[2],
                            l * w[0] + l2 * w[1] + l * w[2],
                            l * w[0] + l * w[1] + l2 * w[2],
                        ]
                    };
                    let scale = 1.0 / (3.0 * volume * self.cell_share(Field::Sxx, k));
                    for (&field, v) in Component::P.fields().iter().zip(loads) {
                        if v != 0.0 {
                            op.terms.push((field, p, v * scale));
                        }
                    }
                }
            }
        }
        op
    }

    /// A component's value at a point (a node, or grid coordinates) now.
    pub fn sample(&self, c: Component, point: impl AsPoint) -> f32 {
        self.receiver_op(c, point).sample(&self.fields)
    }

    /// The largest particle-velocity component anywhere, or infinity if any is not finite.
    pub fn peak_velocity(&self) -> f32 {
        self.fields
            .velocities()
            .iter()
            .map(|v| {
                v.par_iter()
                    .map(|x| {
                        if x.is_finite() {
                            x.abs()
                        } else {
                            f32::INFINITY
                        }
                    })
                    .reduce(|| 0.0, f32::max)
            })
            .fold(0.0, f32::max)
    }
}

/// A linear functional on the wavefield, Σ w · field[index]: how a receiver reads a point. Its
/// transpose, field[index] += w · value, is how a source acts there.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PointOp {
    /// Field, padded index and weight of each term.
    pub terms: Vec<(Field, usize, f32)>,
}

impl PointOp {
    pub fn sample(&self, fields: &Fields) -> f32 {
        self.terms
            .iter()
            .map(|&(f, i, w)| w * fields.get(f)[i])
            .sum()
    }

    pub fn inject(&self, fields: &mut Fields, value: f32) {
        for &(f, i, w) in &self.terms {
            fields.get_mut(f)[i] += w * value;
        }
    }
}

/// What a step back ([`Solver::step_back_with`]) reads from outside the region it rebuilds,
/// the nodes outside the absorbing layers ([`Solver::interior`]): every value within the
/// stencil's reach, its half-length L, of that region along an axis. Velocities are kept on
/// the whole shell, edges included, since the gradient's strain rates just outside the region
/// read them along its faces; each stress only on the faces across which the velocity update
/// differentiates it (σxx on the x faces, σxy on the x and y faces, and so on). Above a free
/// surface there is nothing to keep: the mirrors rebuild it.
#[derive(Clone, Debug)]
pub struct Band {
    /// Runs of padded indices, (start, length): the velocities' shell, shared by vx, vy and vz.
    shell: Vec<(usize, usize)>,
    /// The runs of σxx, σyy, σzz, σxy, σxz and σyz.
    faces: [Vec<(usize, usize)>; 6],
    /// Values of each array in a record: the shell, then each stress's faces.
    shell_len: usize,
    face_len: [usize; 6],
}

/// The values and runs of a [`Band`], counted without building it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BandSize {
    /// Values in one record, one state's band.
    pub values: usize,
    /// Runs of contiguous indices it is read and written through.
    pub runs: usize,
}

impl BandSize {
    /// Bytes the runs take.
    pub fn runs_bytes(&self) -> usize {
        self.runs * std::mem::size_of::<(usize, usize)>()
    }
}

impl Band {
    /// The band around `interior` on `grid`, `reach` nodes deep (the stencil's half-length).
    pub fn new(
        grid: &Grid,
        interior: &[std::ops::Range<usize>; 3],
        free_surface: bool,
        reach: usize,
    ) -> Self {
        assert!(reach <= grid.halo, "the band reaches past the halo");
        let [_, sy, sz] = grid.strides();
        let h = grid.halo as isize;
        let at = |i: isize, j: isize, k: isize| {
            ((i + h) as usize) + ((j + h) as usize) * sy + ((k + h) as usize) * sz
        };
        let l = reach as isize;
        let [x, y, z] = interior
            .each_ref()
            .map(|r| (r.start as isize, r.end as isize));
        // Above a free surface the mirrors stand in for the band.
        let up = if free_surface { 0 } else { l };
        let runs = |out: &mut Vec<(usize, usize)>, j: isize, k: isize, from: isize, to: isize| {
            if to > from {
                out.push((at(from, j, k), (to - from) as usize));
            }
        };
        let inside = |v: isize, r: (isize, isize)| v >= r.0 && v < r.1;
        let mut shell = Vec::new();
        for k in z.0 - up..z.1 + l {
            for j in y.0 - l..y.1 + l {
                if inside(j, y) && inside(k, z) {
                    runs(&mut shell, j, k, x.0 - l, x.0);
                    runs(&mut shell, j, k, x.1, x.1 + l);
                } else {
                    runs(&mut shell, j, k, x.0 - l, x.1 + l);
                }
            }
        }
        let (mut fx, mut fy, mut fz) = (Vec::new(), Vec::new(), Vec::new());
        for k in z.0..z.1 {
            for j in y.0..y.1 {
                runs(&mut fx, j, k, x.0 - l, x.0);
                runs(&mut fx, j, k, x.1, x.1 + l);
            }
            for j in (y.0 - l..y.0).chain(y.1..y.1 + l) {
                runs(&mut fy, j, k, x.0, x.1);
            }
        }
        for k in (z.0 - up..z.0).chain(z.1..z.1 + l) {
            for j in y.0..y.1 {
                runs(&mut fz, j, k, x.0, x.1);
            }
        }
        let both = |a: &Vec<(usize, usize)>, b: &Vec<(usize, usize)>| [a.as_slice(), b].concat();
        let faces = [
            fx.clone(),
            fy.clone(),
            fz.clone(),
            both(&fx, &fy),
            both(&fx, &fz),
            both(&fy, &fz),
        ];
        let total = |r: &[(usize, usize)]| r.iter().map(|&(_, n)| n).sum::<usize>();
        Self {
            shell_len: total(&shell),
            face_len: faces.each_ref().map(|f| total(f)),
            shell,
            faces,
        }
    }

    /// [`Band::new`]'s size, counted.
    pub fn size(
        interior: &[std::ops::Range<usize>; 3],
        free_surface: bool,
        reach: usize,
    ) -> BandSize {
        let [x, y, z] = interior.each_ref().map(|r| r.len());
        let l = reach;
        let up = if free_surface { 0 } else { l };
        let some = |n: usize, runs: usize| if n > 0 { runs } else { 0 };
        let shell = (x + 2 * l) * (y + 2 * l) * (z + up + l) - x * y * z;
        let shell_runs = if l == 0 {
            0
        } else {
            (y + 2 * l) * (z + up + l) - y * z + some(l, 2 * y * z)
        };
        let (fx, fy, fz) = (2 * l * y * z, x * 2 * l * z, x * y * (up + l));
        let runs = some(l, 2 * y * z) + some(x, 2 * l * z) + some(x, y * (up + l));
        BandSize {
            values: 3 * shell + 3 * (fx + fy + fz),
            runs: shell_runs + 3 * runs,
        }
    }

    /// Values in one record.
    pub fn len(&self) -> usize {
        3 * self.shell_len + self.face_len.iter().sum::<usize>()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Bytes its runs take.
    pub fn memory_bytes(&self) -> usize {
        let runs = self.shell.len() + self.faces.iter().map(Vec::len).sum::<usize>();
        runs * std::mem::size_of::<(usize, usize)>()
    }

    /// Copies the band of `fields` into `out`, one record.
    pub fn save(&self, fields: &Fields, out: &mut [f32]) {
        assert_eq!(out.len(), self.len());
        let mut o = 0;
        for (f, runs) in self.arrays() {
            let v = fields.get(f);
            for &(start, n) in runs {
                out[o..o + n].copy_from_slice(&v[start..start + n]);
                o += n;
            }
        }
    }

    /// Writes the velocities of `record` back into `fields`.
    pub fn restore_velocity(&self, fields: &mut Fields, record: &[f32]) {
        self.restore(fields, record, Field::is_velocity);
    }

    /// Writes the stresses of `record` back into `fields`.
    pub fn restore_stress(&self, fields: &mut Fields, record: &[f32]) {
        self.restore(fields, record, |f| !f.is_velocity());
    }

    fn restore(&self, fields: &mut Fields, record: &[f32], which: impl Fn(Field) -> bool) {
        assert_eq!(record.len(), self.len());
        let mut o = 0;
        for (f, runs) in self.arrays() {
            let len: usize = runs.iter().map(|&(_, n)| n).sum();
            if which(f) {
                let v = fields.get_mut(f);
                let mut p = o;
                for &(start, n) in runs {
                    v[start..start + n].copy_from_slice(&record[p..p + n]);
                    p += n;
                }
            }
            o += len;
        }
    }

    /// Each array and its runs, in the order of a record ([`Field::ALL`]'s): the velocities
    /// on the shell, then each stress on its faces.
    fn arrays(&self) -> impl Iterator<Item = (Field, &[(usize, usize)])> {
        let shell = std::iter::repeat_n(self.shell.as_slice(), 3);
        Field::ALL
            .into_iter()
            .zip(shell.chain(self.faces.iter().map(Vec::as_slice)))
    }
}

/// Adds weight `w` to index `m`'s entry in `out`, or makes one.
fn accumulate(out: &mut Vec<(usize, f32)>, m: usize, w: f32) {
    match out.iter_mut().find(|(i, _)| *i == m) {
        Some((_, sum)) => *sum += w,
        None => out.push((m, w)),
    }
}

/// Sets plane `dst` of a padded array to `sign` times plane `src`.
fn mirror_plane(v: &mut [f32], len: usize, dst: usize, src: usize, sign: f32) {
    debug_assert!(dst < src);
    let (above, below) = v.split_at_mut(src * len);
    for (d, s) in above[dst * len..(dst + 1) * len]
        .iter_mut()
        .zip(&below[..len])
    {
        *d = sign * *s;
    }
}

/// out += sign · b · (x + y + z). A sign of ±1 changes nothing but the sign, so taking a
/// step back subtracts exactly what the step added.
#[inline(always)]
fn add_sum(out: &mut [f32], sign: f32, b: &[f32], x: &[f32], y: &[f32], z: &[f32]) {
    let n = out.len();
    let (b, x, y, z) = (&b[..n], &x[..n], &y[..n], &z[..n]);
    for i in 0..n {
        out[i] += sign * (b[i] * (x[i] + y[i] + z[i]));
    }
}

/// out += sign · m · (x + y)
#[inline(always)]
fn add_sum2(out: &mut [f32], sign: f32, m: &[f32], x: &[f32], y: &[f32]) {
    let n = out.len();
    let (m, x, y) = (&m[..n], &x[..n], &y[..n]);
    for i in 0..n {
        out[i] += sign * (m[i] * (x[i] + y[i]));
    }
}

/// Scratch rows for one thread: the nine derivatives of the row being updated.
struct Rows {
    d: [Vec<f32>; 9],
}

impl Rows {
    fn new(n: usize) -> Self {
        Self {
            d: std::array::from_fn(|_| vec![0.0; n]),
        }
    }

    /// The nine rows, `n` values each.
    #[inline(always)]
    fn take(&mut self, n: usize) -> [&mut [f32]; 9] {
        self.d.each_mut().map(|d| &mut d[..n])
    }
}

/// The nodes one task of a sweep updates: nodes `i` of rows `j` of interior planes `k`,
/// plane by plane upwards, so that what the z-derivatives read of one plane is still in the
/// cache for the next ([`block_size`]).
#[derive(Clone, Debug, PartialEq, Eq)]
struct Block {
    i: Range<usize>,
    j: Range<usize>,
    k: Range<usize>,
}

/// Rows in a block's band: see [`block_size`].
const BAND_ROWS: usize = 4;

/// Nodes along x, y and z of the blocks the updates split the interior of `g` into: whole
/// rows, bands of [`BAND_ROWS`] rows, and whole columns of planes.
///
/// A block is worked through plane by plane upwards, so of the 2L planes of three fields that
/// the z-derivatives read for a plane, all but one are read again for the next a moment later
/// and come from the cache instead of memory. Updated a whole plane at a time, as before, a
/// large grid's planes had left the cache by then. Narrow bands keep what a block reads per
/// plane small; the y-derivatives read 2L rows more than the band's own, so very narrow ones
/// read more again. Measured on gauss (Xeon E5-1650 v3, 15 MiB L3, 12 threads; the research
/// log's 2026-11-25 entry), bands of 2–4 whole rows were fastest on every grid tried, and
/// parts of rows slower: shorter rows defeat the hardware prefetchers.
fn block_size(g: &Grid) -> [usize; 3] {
    [g.n[0], BAND_ROWS.min(g.n[1]), g.n[2]]
}

/// Splits `0..n` into the fewest nearly equal parts of at most `size`.
fn split(n: usize, size: usize) -> impl Iterator<Item = Range<usize>> + Clone {
    let parts = n.div_ceil(size.max(1)).max(1);
    (0..parts).map(move |p| p * n / parts..(p + 1) * n / parts)
}

/// The blocks of at most `size` nodes along each axis that tile the interior of `g`, each
/// node in exactly one: what makes [`Outputs::row`] sound. Within a column of blocks the
/// planes come in order, so a thread that takes consecutive blocks marches up the column.
fn blocks(g: &Grid, size: [usize; 3]) -> Vec<Block> {
    let [nx, ny, nz] = g.n;
    let mut out = Vec::new();
    for i in split(nx, size[0]) {
        for j in split(ny, size[1]) {
            for k in split(nz, size[2]) {
                out.push(Block {
                    i: i.clone(),
                    j: j.clone(),
                    k,
                });
            }
        }
    }
    out
}

/// A mutable array whose disjoint parts the tasks of a sweep write at once. Splitting it
/// ahead of time, as `chunks_mut` does whole planes, would take a slice per row of every
/// block and every array, each update; handing the rows out as they are reached takes nothing.
struct Shared<'a> {
    ptr: *mut f32,
    len: usize,
    array: PhantomData<&'a mut [f32]>,
}

// SAFETY: it is a `&mut [f32]` whose parts go to different threads, each part to one only
// (`Shared::part`'s contract).
unsafe impl Send for Shared<'_> {}
unsafe impl Sync for Shared<'_> {}

impl<'a> Shared<'a> {
    fn new(v: &'a mut [f32]) -> Self {
        Self {
            ptr: v.as_mut_ptr(),
            len: v.len(),
            array: PhantomData,
        }
    }

    /// `len` values from `start`.
    ///
    /// # Safety
    ///
    /// Nothing else may read or write them while the slice lives.
    #[inline(always)]
    unsafe fn part(&self, start: usize, len: usize) -> &'a mut [f32] {
        assert!(
            start <= self.len && len <= self.len - start,
            "outside the array"
        );
        // SAFETY: inside the array, which outlives 'a, and the caller's alone.
        unsafe { std::slice::from_raw_parts_mut(self.ptr.add(start), len) }
    }
}

/// What a sweep writes: its output arrays, padded, and the CPML memory along x, y and z
/// where it keeps any (see [`Outputs::row`] for the layout).
struct Outputs<'a, const N: usize> {
    out: [Shared<'a>; N],
    psi: [Option<Shared<'a>>; 3],
}

/// What a row of a block writes: its nodes in each output array, and the CPML memory for
/// them, each of three terms: along x if they reach into the x layers, along y if the row
/// is in a y layer and along z if its plane is in a z layer, these two with the slab index.
struct RowOut<'a, const N: usize> {
    out: [&'a mut [f32]; N],
    x: Option<[&'a mut [f32]; 3]>,
    y: Option<(usize, [&'a mut [f32]; 3])>,
    z: Option<(usize, [&'a mut [f32]; 3])>,
}

impl<'a, const N: usize> Outputs<'a, N> {
    /// `out` the padded output arrays; `psi` the CPML memory, empty where the sweep keeps
    /// none.
    fn new(out: [&'a mut Vec<f32>; N], psi: &'a mut [Vec<f32>; 3]) -> Self {
        Self {
            out: out.map(|v| Shared::new(v)),
            psi: psi
                .each_mut()
                .map(|v| (!v.is_empty()).then(|| Shared::new(v))),
        }
    }

    /// Nodes `i` of row `j` of interior plane `k`. The memory along x holds, per plane and
    /// row, each term's slab ([3][slab]); along y, per plane and slab row, each term's row
    /// ([3][nx]); along z, per slab plane and row, the same.
    ///
    /// # Safety
    ///
    /// No other task may hold any of the same nodes' values at the same time.
    #[inline(always)]
    unsafe fn row(
        &self,
        g: &Grid,
        cpml: &[CpmlAxis; 3],
        i: Range<usize>,
        j: usize,
        k: usize,
    ) -> RowOut<'a, N> {
        let [nx, ny, _] = g.n;
        let n = i.len();
        let start = g.idx(i.start, j, k);
        // SAFETY (all four): these nodes' values, the caller's alone. A node's memory along x
        // is at its slab index, along y and z at its x index, so other nodes' are elsewhere.
        let out = self.out.each_ref().map(|a| unsafe { a.part(start, n) });
        let terms = |a: &Shared<'a>, base: usize, term: usize, from: usize, len: usize| {
            std::array::from_fn(|t| unsafe { a.part(base + t * term + from, len) })
        };
        let x = self.psi[0].as_ref().and_then(|a| {
            let (w, slabs) = (cpml[0].width(), cpml[0].slabs(i.clone()));
            (!slabs.is_empty()).then(|| terms(a, (k * ny + j) * 3 * w, w, slabs.start, slabs.len()))
        });
        let y = self.psi[1].as_ref().zip(cpml[1].slab(j)).map(|(a, s)| {
            let base = (k * cpml[1].width() + s) * 3 * nx;
            (s, terms(a, base, nx, i.start, n))
        });
        let z = self.psi[2]
            .as_ref()
            .zip(cpml[2].slab(k))
            .map(|(a, s)| (s, terms(a, (s * ny + j) * 3 * nx, nx, i.start, n)));
        RowOut { out, x, y, z }
    }
}

/// One update's arithmetic on nodes `i` of row `j` of interior plane `k`, with the scratch
/// rows `d`, writing `row`. [`sweep`] compiles it once for each instruction set
/// (crate::simd), inlining it, so it and everything it calls in its loops must be
/// `#[inline(always)]`: a call left out of line would run the plain loops.
trait Update<const N: usize>: Sync {
    fn row(&self, d: [&mut [f32]; 9], i: Range<usize>, j: usize, k: usize, row: RowOut<'_, N>);
}

/// Runs `update` on every node of `blocks`, which tile the interior ([`blocks`]), in
/// parallel and with the loops `kernels` names, writing `outputs`.
fn sweep<const N: usize>(
    kernels: Kernels,
    update: &impl Update<N>,
    g: &Grid,
    cpml: &[CpmlAxis; 3],
    blocks: &[Block],
    outputs: Outputs<'_, N>,
) {
    let width = blocks.iter().map(|b| b.i.len()).max().unwrap_or(0);
    blocks.par_iter().for_each_init(
        || Rows::new(width),
        |rows, block| {
            #[cfg(target_arch = "x86_64")]
            if kernels == Kernels::Avx2 {
                // SAFETY: a solver holds Kernels::Avx2 only where the CPU has AVX2
                // (Solver::set_kernels, Kernels::chosen).
                return unsafe { block_avx2(update, g, cpml, &outputs, rows, block) };
            }
            #[cfg(not(target_arch = "x86_64"))]
            debug_assert_eq!(kernels, Kernels::Plain);
            block_plain(update, g, cpml, &outputs, rows, block);
        },
    );
}

/// A block's rows, plane by plane upwards.
#[inline(always)]
fn block<const N: usize>(
    update: &impl Update<N>,
    g: &Grid,
    cpml: &[CpmlAxis; 3],
    outputs: &Outputs<'_, N>,
    rows: &mut Rows,
    block: &Block,
) {
    for k in block.k.clone() {
        for j in block.j.clone() {
            // SAFETY: the blocks tile the interior, so no other task holds these nodes, and
            // this one holds a row at a time.
            let row = unsafe { outputs.row(g, cpml, block.i.clone(), j, k) };
            update.row(rows.take(block.i.len()), block.i.clone(), j, k, row);
        }
    }
}

fn block_plain<const N: usize>(
    update: &impl Update<N>,
    g: &Grid,
    cpml: &[CpmlAxis; 3],
    outputs: &Outputs<'_, N>,
    rows: &mut Rows,
    b: &Block,
) {
    block(update, g, cpml, outputs, rows, b);
}

/// [`block_plain`] compiled with AVX2, without FMA: the same operations, eight values at a
/// time.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
fn block_avx2<const N: usize>(
    update: &impl Update<N>,
    g: &Grid,
    cpml: &[CpmlAxis; 3],
    outputs: &Outputs<'_, N>,
    rows: &mut Rows,
    b: &Block,
) {
    block(update, g, cpml, outputs, rows, b);
}

/// What the velocity update reads.
struct VelocityUpdate<'a> {
    g: Grid,
    sign: f32,
    /// 1/h along x, y and z.
    ih: [f32; 3],
    /// σxx, σyy, σzz, σxy, σxz, σyz.
    stress: [&'a [f32]; 6],
    m: &'a Material,
    st: &'a Stencil,
    cpml: &'a [CpmlAxis; 3],
}

impl Update<3> for VelocityUpdate<'_> {
    #[inline(always)]
    fn row(&self, d: [&mut [f32]; 9], i: Range<usize>, j: usize, kz: usize, row: RowOut<'_, 3>) {
        let Self {
            g,
            sign,
            ih,
            stress: [sxx, syy, szz, sxy, sxz, syz],
            m,
            st,
            cpml: [cx, cy, cz],
        } = self;
        let [_, sy, sz] = g.strides();
        let [d0, d1, d2, d3, d4, d5, d6, d7, d8] = d;
        let r = g.idx(i.start, j, kz);
        let n = i.len();
        st.forward(d0, sxx, r, 1, ih[0]); // ∂x σxx at i+½
        st.backward(d1, sxy, r, sy, ih[1]); // ∂y σxy at j
        st.backward(d2, sxz, r, sz, ih[2]); // ∂z σxz at k
        st.backward(d3, sxy, r, 1, ih[0]); // ∂x σxy at i
        st.forward(d4, syy, r, sy, ih[1]); // ∂y σyy at j+½
        st.backward(d5, syz, r, sz, ih[2]); // ∂z σyz at k
        st.backward(d6, sxz, r, 1, ih[0]); // ∂x σxz at i
        st.backward(d7, syz, r, sy, ih[1]); // ∂y σyz at j
        st.forward(d8, szz, r, sz, ih[2]); // ∂z σzz at k+½

        if let Some([p0, p1, p2]) = row.x {
            cx.apply_along(i.clone(), d0, p0, Stagger::Half);
            cx.apply_along(i.clone(), d3, p1, Stagger::Node);
            cx.apply_along(i.clone(), d6, p2, Stagger::Node);
        }
        if let Some((s, [p0, p1, p2])) = row.y {
            cy.apply_across(s, d1, p0, Stagger::Node);
            cy.apply_across(s, d4, p1, Stagger::Half);
            cy.apply_across(s, d7, p2, Stagger::Node);
        }
        if let Some((s, [p0, p1, p2])) = row.z {
            cz.apply_across(s, d2, p0, Stagger::Node);
            cz.apply_across(s, d5, p1, Stagger::Node);
            cz.apply_across(s, d8, p2, Stagger::Half);
        }

        let [vx, vy, vz] = row.out;
        add_sum(vx, *sign, &m.bx[r..r + n], d0, d1, d2);
        add_sum(vy, *sign, &m.by[r..r + n], d3, d4, d5);
        add_sum(vz, *sign, &m.bz[r..r + n], d6, d7, d8);
    }
}

/// What the stress update reads.
struct StressUpdate<'a> {
    g: Grid,
    sign: f32,
    ih: [f32; 3],
    free_surface: bool,
    /// vx, vy, vz.
    v: [&'a [f32]; 3],
    m: &'a Material,
    st: &'a Stencil,
    cpml: &'a [CpmlAxis; 3],
}

impl Update<6> for StressUpdate<'_> {
    #[inline(always)]
    fn row(&self, d: [&mut [f32]; 9], i: Range<usize>, j: usize, kz: usize, row: RowOut<'_, 6>) {
        let Self {
            g,
            sign,
            ih,
            free_surface,
            v: [vx, vy, vz],
            m,
            st,
            cpml: [cx, cy, cz],
        } = self;
        let sign = *sign;
        let [_, sy, sz] = g.strides();
        let surface = *free_surface && kz == 0;
        let [d0, d1, d2, d3, d4, d5, d6, d7, d8] = d;
        let r = g.idx(i.start, j, kz);
        let n = i.len();
        st.backward(d0, vx, r, 1, ih[0]); // ∂x vx at i
        st.backward(d1, vy, r, sy, ih[1]); // ∂y vy at j
        if surface {
            d2.fill(0.0); // ∂z vz: replaced by the condition σzz = 0
        } else {
            st.backward(d2, vz, r, sz, ih[2]); // ∂z vz at k
        }
        st.forward(d3, vx, r, sy, ih[1]); // ∂y vx at j+½
        st.forward(d4, vy, r, 1, ih[0]); // ∂x vy at i+½
        st.forward(d5, vx, r, sz, ih[2]); // ∂z vx at k+½
        st.forward(d7, vy, r, sz, ih[2]); // ∂z vy at k+½
        st.forward(d6, vz, r, 1, ih[0]); // ∂x vz at i+½
        st.forward(d8, vz, r, sy, ih[1]); // ∂y vz at j+½

        if let Some([p0, p1, p2]) = row.x {
            cx.apply_along(i.clone(), d0, p0, Stagger::Node);
            cx.apply_along(i.clone(), d4, p1, Stagger::Half);
            cx.apply_along(i.clone(), d6, p2, Stagger::Half);
        }
        if let Some((s, [p0, p1, p2])) = row.y {
            cy.apply_across(s, d1, p0, Stagger::Node);
            cy.apply_across(s, d3, p1, Stagger::Half);
            cy.apply_across(s, d8, p2, Stagger::Half);
        }
        if let Some((s, [p0, p1, p2])) = row.z {
            cz.apply_across(s, d2, p0, Stagger::Node);
            cz.apply_across(s, d5, p1, Stagger::Half);
            cz.apply_across(s, d7, p2, Stagger::Half);
        }

        let [sxx, syy, szz, sxy, sxz, syz] = row.out;
        let lam = &m.lam[r..r + n];
        let l2m = &m.l2m[r..r + n];
        let (sxx, syy, szz) = (&mut sxx[..n], &mut syy[..n], &mut szz[..n]);
        let (d0, d1, d2) = (&d0[..n], &d1[..n], &d2[..n]);
        if surface {
            for i in 0..n {
                // With σzz = 0, ∂z vz = -λ/(λ+2μ) (∂x vx + ∂y vy).
                let (l, l2) = (lam[i], l2m[i]);
                let q = if l2 > 0.0 { l * l / l2 } else { 0.0 };
                let (a, b) = (l2 - q, l - q);
                sxx[i] += sign * (a * d0[i] + b * d1[i]);
                syy[i] += sign * (b * d0[i] + a * d1[i]);
                szz[i] = 0.0;
            }
        } else {
            for i in 0..n {
                let (l, l2) = (lam[i], l2m[i]);
                sxx[i] += sign * (l2 * d0[i] + l * (d1[i] + d2[i]));
                syy[i] += sign * (l2 * d1[i] + l * (d0[i] + d2[i]));
                szz[i] += sign * (l2 * d2[i] + l * (d0[i] + d1[i]));
            }
        }
        add_sum2(sxy, sign, &m.mxy[r..r + n], d3, d4);
        add_sum2(sxz, sign, &m.mxz[r..r + n], d5, d6);
        add_sum2(syz, sign, &m.myz[r..r + n], d7, d8);
    }
}

/// What the strain rates read.
struct StrainRate<'a> {
    g: Grid,
    ih: [f32; 3],
    free_surface: bool,
    /// vx, vy, vz.
    v: [&'a [f32]; 3],
    st: &'a Stencil,
    /// Whether to add them to the output instead of writing them.
    add: bool,
}

impl Update<6> for StrainRate<'_> {
    #[inline(always)]
    fn row(&self, d: [&mut [f32]; 9], i: Range<usize>, j: usize, kz: usize, row: RowOut<'_, 6>) {
        let Self {
            g,
            ih,
            free_surface,
            v: [vx, vy, vz],
            st,
            add,
        } = self;
        let [_, sy, sz] = g.strides();
        let [d0, d1, d2, d3, d4, d5, d6, d7, d8] = d;
        let r = g.idx(i.start, j, kz);
        let n = i.len();
        st.backward(d0, vx, r, 1, ih[0]);
        st.backward(d1, vy, r, sy, ih[1]);
        if *free_surface && kz == 0 {
            d2.fill(0.0);
        } else {
            st.backward(d2, vz, r, sz, ih[2]);
        }
        st.forward(d3, vx, r, sy, ih[1]);
        st.forward(d4, vy, r, 1, ih[0]);
        st.forward(d5, vx, r, sz, ih[2]);
        st.forward(d6, vz, r, 1, ih[0]);
        st.forward(d7, vy, r, sz, ih[2]);
        st.forward(d8, vz, r, sy, ih[1]);
        let [exx, eyy, ezz, exy, exz, eyz] = row.out;
        let (exy, exz, eyz) = (&mut exy[..n], &mut exz[..n], &mut eyz[..n]);
        if *add {
            // Each sum of two derivatives is rounded before it is added, as if written first.
            for i in 0..n {
                exx[i] += d0[i];
                eyy[i] += d1[i];
                ezz[i] += d2[i];
                exy[i] += d3[i] + d4[i];
                exz[i] += d5[i] + d6[i];
                eyz[i] += d7[i] + d8[i];
            }
        } else {
            exx.copy_from_slice(d0);
            eyy.copy_from_slice(d1);
            ezz.copy_from_slice(d2);
            for i in 0..n {
                exy[i] = d3[i] + d4[i];
                exz[i] = d5[i] + d6[i];
                eyz[i] = d7[i] + d8[i];
            }
        }
    }
}

/// What the unit tests across the crate make solvers with.
#[cfg(test)]
pub(crate) mod testing {
    use super::*;

    /// Layers with a quadratic profile set for 3000 m/s, and no sponge.
    pub fn pml() -> PmlSettings {
        PmlSettings {
            power: 2.0,
            reflection: None,
            kappa_max: 1.0,
            alpha_max: 30.0,
            vp_max: 3000.0,
            sponge: 0.0,
        }
    }

    /// A solver at rest in a homogeneous solid (vp 3000 m/s, vs 1700 m/s, ρ 2500 kg/m³) on
    /// `grid`, stepping 1 ms, with layers `pml_width` nodes thick ([`pml`]).
    pub fn solver(grid: Grid, pml_width: usize, free_surface: bool) -> Solver {
        let model = Model::homogeneous(grid, 3000.0, 1700.0, 2500.0);
        Solver::new(model, 1e-3, pml_width, free_surface, &pml())
    }
}

#[cfg(test)]
mod tests {
    use super::testing::{pml, solver};
    use super::*;

    #[test]
    fn harmonic_mean_vanishes_next_to_a_fluid() {
        assert_eq!(harmonic_mean([1.0, 1.0, 1.0, 0.0]), 0.0);
        assert!((harmonic_mean([2.0, 2.0, 2.0, 2.0]) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn a_wavefield_at_rest_stays_at_rest() {
        let g = Grid::new([16, 14, 12], [10.0; 3]);
        let mut s = solver(g, 4, true);
        for _ in 0..5 {
            s.step();
        }
        assert_eq!(s.peak_velocity(), 0.0);
    }

    /// The free surface's mirror images make σzz odd about the surface plane, and σxz odd
    /// about the plane half a node above the first σxz row.
    #[test]
    fn mirrors_stress_across_the_free_surface() {
        let g = Grid::new([10, 10, 10], [10.0; 3]);
        let mut s = solver(g, 3, true);
        let p = g.plane_len();
        for (n, v) in s.fields.szz.iter_mut().enumerate() {
            *v = n as f32;
        }
        s.fields.sxz.clone_from(&s.fields.szz);
        s.mirror_stress();
        let at = |k: usize| (k * p)..(k * p + p);
        let f = &s.fields;
        for (a, b) in f.szz[at(1)].iter().zip(&f.szz[at(3)]) {
            assert_eq!(*a, -*b);
        }
        for (a, b) in f.szz[at(0)].iter().zip(&f.szz[at(4)]) {
            assert_eq!(*a, -*b);
        }
        for (a, b) in f.sxz[at(1)].iter().zip(&f.sxz[at(2)]) {
            assert_eq!(*a, -*b);
        }
        for (a, b) in f.sxz[at(0)].iter().zip(&f.sxz[at(3)]) {
            assert_eq!(*a, -*b);
        }
    }

    /// Velocities mirror with the same sign: vx about the surface plane, vz about the plane
    /// half a node above its first row.
    #[test]
    fn mirrors_velocity_across_the_free_surface() {
        let g = Grid::new([10, 10, 10], [10.0; 3]);
        let mut s = solver(g, 3, true);
        let p = g.plane_len();
        for (n, v) in s.fields.vx.iter_mut().enumerate() {
            *v = n as f32;
        }
        s.fields.vz.clone_from(&s.fields.vx);
        s.mirror_velocity();
        let at = |k: usize| (k * p)..(k * p + p);
        let f = &s.fields;
        assert_eq!(f.vx[at(1)], f.vx[at(3)]);
        assert_eq!(f.vx[at(0)], f.vx[at(4)]);
        assert_eq!(f.vz[at(1)], f.vz[at(2)]);
        assert_eq!(f.vz[at(0)], f.vz[at(3)]);
    }

    #[test]
    fn interpolation_weights_sum_to_one_and_fold_at_the_free_surface() {
        let g = Grid::new([10, 10, 10], [10.0; 3]);
        let (open, free) = (solver(g, 3, false), solver(g, 3, true));
        let sum = |w: &[(usize, f32)]| w.iter().map(|&(_, w)| w).sum::<f32>();
        for at in 0..10 {
            for s in [&open, &free] {
                assert!((sum(&s.stagger_weights(2, at, 1.0)) - 1.0).abs() < 1e-6);
            }
        }
        assert_eq!(open.stagger_weights(0, 5, 1.0).len(), 4);
        assert_eq!(open.stagger_weights(0, 1, 1.0), vec![(0, 0.5), (1, 0.5)]);
        assert_eq!(open.stagger_weights(0, 0, 1.0), vec![(0, 1.5), (1, -0.5)]);
        // vz is even about the surface: the cubic folds onto the two values below it.
        assert_eq!(
            free.stagger_weights(2, 0, 1.0),
            vec![(1, -0.125), (0, 1.125)]
        );
        // σxz is odd: on the surface it vanishes.
        assert!(free.stagger_weights(2, 0, -1.0).is_empty());
        // Away from the surface, mirroring changes nothing.
        assert_eq!(
            free.stagger_weights(2, 5, -1.0),
            open.stagger_weights(2, 5, 1.0)
        );
    }

    /// With `surface_slope`, fields mirrored evenly are continued above a free surface with
    /// their slope, so the interpolation reads anything with a slope there, here
    /// a + b z + c z², exactly, on the surface and between it and the stencil's reach, on half
    /// nodes and on nodes, at orders 4 and 8. Without it vz on the surface reads b as if from
    /// 0.34 h below at order 8.
    #[test]
    fn interpolation_continues_even_fields_with_their_slope_across_the_free_surface() {
        let f = |z: f64| 2.0 + 3.0 * z - 0.7 * z * z;
        for order in [4, 8] {
            let g = Grid::new([12, 12, 12], [10.0; 3]).with_halo(order / 2);
            let s = solver(g, 0, true)
                .with_stencil(Stencil::taylor(order).unwrap())
                .with_surface_slope(true);
            for (shift, x) in [
                (0.5, 0.0),
                (0.5, 0.3),
                (0.5, 1.0),
                (0.5, 2.7),
                (0.0, 0.3),
                (0.0, 1.5),
            ] {
                let w = s.axis_weights(2, x, shift, 1.0);
                let read: f64 = w.iter().map(|&(m, w)| w as f64 * f(m as f64 + shift)).sum();
                assert!(
                    (read - f(x)).abs() < 1e-5,
                    "order {order}, shift {shift}, z {x}: {read} for {}",
                    f(x)
                );
            }
            // Odd fields vanish on the surface: an odd function is folded exactly as it is.
            let odd = |z: f64| 3.0 * z;
            let w = s.axis_weights(2, 0.3, 0.5, -1.0);
            let read: f64 = w.iter().map(|&(m, w)| w as f64 * odd(m as f64 + 0.5)).sum();
            assert!((read - odd(0.3)).abs() < 1e-5, "{read}");
        }
    }

    /// The interpolation is as long as the stencil, and exact for polynomials up to one
    /// degree less than its length.
    #[test]
    fn interpolation_matches_the_stencil() {
        let g = Grid::new([20, 20, 20], [10.0; 3]).with_halo(4);
        let s = solver(g, 3, false).with_stencil(Stencil::taylor(8).unwrap());
        let eighth: Vec<(usize, f32)> = [-5.0, 49.0, -245.0, 1225.0, 1225.0, -245.0, 49.0, -5.0]
            .iter()
            .enumerate()
            .map(|(k, w)| (6 + k, w / 2048.0))
            .collect();
        assert_eq!(s.stagger_weights(0, 10, 1.0), eighth);
        // A degree-7 polynomial at the staggered positions m + ½, read at node 10.
        let p = |x: f64| (x - 3.0).powi(7) / 1e5 - x * x;
        let read: f64 = s
            .stagger_weights(0, 10, 1.0)
            .iter()
            .map(|&(m, w)| w as f64 * p(m as f64 + 0.5))
            .sum();
        assert!((read - p(10.0)).abs() < 1e-4 * p(10.0).abs(), "{read}");
        // Shorter where fewer values exist: six points two nodes from the far end.
        assert_eq!(s.stagger_weights(0, 17, 1.0).len(), 6);
        assert_eq!(s.stagger_weights(0, 2, 1.0).len(), 4);
    }

    /// The blocks tile the interior, and the rows they hand out cover every interior value
    /// of the output arrays and every value of the CPML memory exactly once, whatever their
    /// size: what makes `Outputs::row` sound.
    #[test]
    fn blocks_hand_out_every_value_once() {
        let g = Grid::new([19, 11, 13], [1.0; 3]);
        let [nx, ny, nz] = g.n;
        let cpml: [CpmlAxis; 3] = std::array::from_fn(|a| {
            CpmlAxis::new(g.n[a], if a == 2 { 0 } else { 3 }, 3, 1.0, 1e-3, &pml())
        });
        for size in [
            [19, 11, 1],
            [5, 4, 3],
            [8, 11, 13],
            [1, 1, 1],
            [100, 100, 100],
        ] {
            let mut a = vec![0.0f32; g.padded_len()];
            let mut b = a.clone();
            let mut psi = [
                vec![0.0; nz * ny * 3 * 6],
                vec![0.0; nz * 6 * 3 * nx],
                vec![0.0; 3 * ny * 3 * nx],
            ];
            let outputs = Outputs::new([&mut a, &mut b], &mut psi);
            for block in blocks(&g, size) {
                let sizes = [&block.i, &block.j, &block.k].map(Range::len);
                assert!(
                    (0..3).all(|a| sizes[a] >= 1 && sizes[a] <= size[a]),
                    "{block:?}"
                );
                for k in block.k.clone() {
                    for j in block.j.clone() {
                        // SAFETY: one row at a time, in one thread.
                        let row = unsafe { outputs.row(&g, &cpml, block.i.clone(), j, k) };
                        let y = row.y.into_iter().flat_map(|(_, t)| t);
                        let z = row.z.into_iter().flat_map(|(_, t)| t);
                        let parts = row.out.into_iter().chain(row.x.into_iter().flatten());
                        for part in parts.chain(y).chain(z) {
                            part.iter_mut().for_each(|v| *v += 1.0);
                        }
                    }
                }
            }
            for v in [&a, &b] {
                let interior = |p: usize| {
                    let (i, j, k) = (
                        p % g.padded()[0],
                        p / g.padded()[0] % g.padded()[1],
                        p / g.plane_len(),
                    );
                    [i, j, k]
                        .iter()
                        .zip(g.n)
                        .all(|(&x, n)| x >= g.halo && x < g.halo + n)
                };
                assert!((0..v.len()).all(|p| v[p] == if interior(p) { 1.0 } else { 0.0 }));
            }
            assert!(psi.iter().flatten().all(|&v| v == 1.0), "{size:?}");
        }
    }

    fn rel(a: f32, b: f64) -> bool {
        ((a as f64 - b) / b).abs() < 1e-6
    }

    /// Density is averaged arithmetically and μ harmonically onto the staggered positions.
    #[test]
    fn material_takes_means_at_staggered_positions() {
        let g = Grid::new([8, 8, 8], [10.0; 3]);
        let mut model = Model::homogeneous(g, 2000.0, 1000.0, 2000.0);
        for k in 4..8 {
            for j in 0..8 {
                for i in 0..8 {
                    let n = g.interior_idx(i, j, k);
                    model.vp[n] = 3000.0;
                    model.vs[n] = 1500.0;
                    model.rho[n] = 3000.0;
                }
            }
        }
        let dt = 1e-3;
        let m = Material::new(&model, dt);
        let (mu1, mu2) = (2000.0 * 1e6, 3000.0 * 1500.0 * 1500.0);
        let p = g.idx(3, 3, 3); // the last node above the interface
        assert!(rel(m.bx[p], dt / 2000.0));
        assert!(rel(m.bz[p], dt * 2.0 / 5000.0));
        assert!(rel(m.mxy[p], dt * mu1));
        assert!(rel(m.mxz[p], dt * 2.0 * mu1 * mu2 / (mu1 + mu2)));
        assert!(rel(m.l2m[p], dt * 2000.0 * 4e6));
        assert!(rel(m.lam[p], dt * (2000.0 * 4e6 - 2.0 * mu1)));
    }

    /// Injecting through an operator is its transpose: ⟨inject(a), g⟩ = a · sample(g).
    #[test]
    fn point_ops_inject_the_transpose_of_what_they_sample() {
        let g = Grid::new([10, 10, 10], [10.0; 3]);
        let s = solver(g, 3, true);
        let mut probe = Fields::zeros(g.padded_len());
        for (n, v) in probe.vz.iter_mut().enumerate() {
            *v = ((n * 7919) % 1000) as f32 / 1000.0 - 0.5;
        }
        probe.vx.clone_from(&probe.vz);
        probe.sxx.clone_from(&probe.vz);
        for (c, node) in [
            (Component::Vz, [5, 5, 0]),
            (Component::Vz, [5, 5, 1]),
            (Component::Vx, [0, 5, 5]),
            (Component::P, [5, 5, 5]),
        ] {
            let op = s.receiver_op(c, node);
            let mut f = Fields::zeros(g.padded_len());
            op.inject(&mut f, 2.0);
            let dot: f32 = [Field::Vx, Field::Vz, Field::Sxx, Field::Syy, Field::Szz]
                .iter()
                .map(|&k| {
                    f.get(k)
                        .iter()
                        .zip(probe.get(k))
                        .map(|(a, b)| a * b)
                        .sum::<f32>()
                })
                .sum();
            assert!(
                (dot - 2.0 * op.sample(&probe)).abs() < 1e-5,
                "{c:?} at {node:?}"
            );
        }
    }

    #[test]
    fn values_on_the_free_surface_stand_for_half_a_cell() {
        let g = Grid::new([10, 10, 10], [10.0; 3]);
        let (free, open) = (solver(g, 3, true), solver(g, 3, false));
        for f in [
            Field::Vx,
            Field::Vy,
            Field::Sxx,
            Field::Syy,
            Field::Szz,
            Field::Sxy,
        ] {
            assert_eq!(free.cell_share(f, 0), 0.5, "{f:?}");
            assert_eq!(free.cell_share(f, 1), 1.0, "{f:?}");
            assert_eq!(open.cell_share(f, 0), 1.0, "{f:?}");
        }
        for f in [Field::Vz, Field::Sxz, Field::Syz] {
            assert_eq!(free.cell_share(f, 0), 1.0, "{f:?}");
        }
    }

    /// Restoring a state replays the run from there bit for bit, threads and all.
    #[test]
    fn restoring_a_state_replays_the_run_exactly() {
        let g = Grid::new([20, 18, 16], [10.0; 3]).with_halo(4);
        let model = Model::homogeneous(g, 3000.0, 1700.0, 2500.0);
        let pml = PmlSettings {
            sponge: 0.2,
            ..pml()
        };
        let mut s =
            Solver::new(model, 5e-4, 4, true, &pml).with_stencil(Stencil::taylor(8).unwrap());
        let kick = |f: &mut Fields, t: f64| f.vz[g.idx(10, 9, 3)] += (t * 400.0).sin() as f32;
        for _ in 0..30 {
            s.step_with(kick, |_, _| {});
        }
        let saved = s.state();
        for _ in 0..25 {
            s.step_with(kick, |_, _| {});
        }
        let first = s.state();
        s.restore(&saved);
        assert_eq!(s.step, 30);
        for _ in 0..25 {
            s.step_with(kick, |_, _| {});
        }
        assert_eq!(s.fields.vz, first.fields.vz);
        assert_eq!(s.fields.sxy, first.fields.sxy);
        assert_eq!(s.psi_s, first.psi_s);
        assert!(s.peak_velocity() > 0.0);
        assert_eq!(s.at_rest().peak_velocity(), 0.0);
    }

    /// A run with a free surface over two media, layers with their sponge, a force and a
    /// moment between nodes, with the update loops `kernels` names and the nodes split into
    /// blocks of `size` (by default as the solver splits them): every value of its state, its
    /// strain rates and its state after stepping back again, as bits.
    fn bits_of_a_run(kernels: Kernels, size: Option<[usize; 3]>) -> Vec<Vec<u32>> {
        let g = Grid::new([21, 19, 17], [10.0; 3]).with_halo(4);
        let mut model = Model::homogeneous(g, 3000.0, 1700.0, 2500.0);
        for n in 9 * 21 * 19..g.interior_len() {
            (model.vp[n], model.vs[n], model.rho[n]) = (3600.0, 2000.0, 2600.0);
        }
        let pml = PmlSettings {
            sponge: 0.2,
            ..pml()
        };
        let steps = 40;
        let mut s =
            Solver::new(model, 5e-4, 5, true, &pml).with_stencil(Stencil::taylor(8).unwrap());
        assert_eq!(s.set_kernels(kernels), kernels);
        if let Some(size) = size {
            s.blocks = blocks(&g, size);
        }
        let force = s.receiver_op(Component::Vz, [10.3, 9.0, 2.6]);
        let moment = s.receiver_op(Component::P, [9.0, 8.5, 0.0]);
        let wavelet = |t: f64| (t * 400.0).sin() as f32 * 1e3;
        let band = Band::new(&g, &s.interior(), true, s.stencil.half_length());
        let len = band.len();
        let mut records = vec![0.0; steps * len];
        for n in 0..steps {
            band.save(&s.fields, &mut records[n * len..(n + 1) * len]);
            s.step_with(
                |f, t| force.inject(f, wavelet(t)),
                |f, t| moment.inject(f, wavelet(t)),
            );
        }
        let state = s.state();
        assert!(state.fields.vz.iter().any(|&v| v != 0.0));
        let mut strain: [Vec<f32>; 6] = std::array::from_fn(|_| vec![0.0; g.padded_len()]);
        s.strain_rate(s.fields.velocities(), &mut strain);
        for n in (0..steps).rev() {
            s.step_back_with(
                &band,
                &records[n * len..(n + 1) * len],
                |f, t| force.inject(f, -wavelet(t)),
                |f, t| moment.inject(f, -wavelet(t)),
            );
        }
        let arrays = state.fields.arrays().into_iter().chain(&state.psi_v);
        let arrays = arrays.chain(&state.psi_s).chain(&strain);
        arrays
            .chain(s.fields.arrays())
            .map(|a| a.iter().map(|x| x.to_bits()).collect())
            .collect()
    }

    /// The AVX2 loops do the plain loops' arithmetic eight values at a time (crate::simd), so
    /// they step, take strain rates and step back to the same values bit for bit; rows of 21
    /// nodes leave a remainder after the vectors. Without AVX2 there is nothing to compare,
    /// and the plain loops run alone.
    #[test]
    fn the_avx2_loops_step_as_the_plain_ones_bit_for_bit() {
        let plain = bits_of_a_run(Kernels::Plain, None);
        if Kernels::Avx2.available() {
            let avx2 = bits_of_a_run(Kernels::Avx2, None);
            assert_eq!(plain.len(), avx2.len());
            for (a, (p, v)) in plain.iter().zip(&avx2).enumerate() {
                assert!(p == v, "array {a} differs");
            }
        }
    }

    /// However the nodes are split into blocks, every node's arithmetic is the same, so the
    /// results are too, bit for bit: a block per plane (as the solver split them before), per
    /// node, uneven ones, and one for the whole grid.
    #[test]
    fn any_split_into_blocks_steps_the_same_bit_for_bit() {
        let kernels = Kernels::chosen();
        let reference = bits_of_a_run(kernels, Some([21, 19, 1]));
        for size in [[1, 1, 1], [5, 4, 3], [8, 19, 17], [21, 2, 6], [64, 64, 64]] {
            let run = bits_of_a_run(kernels, Some(size));
            for (a, (r, v)) in reference.iter().zip(&run).enumerate() {
                assert!(r == v, "{size:?}: array {a} differs");
            }
        }
    }

    /// Adding the strain rates to a sum in the sweep gives what writing them and adding them
    /// gave, bit for bit, free surface included.
    #[test]
    fn adds_strain_rates_as_if_written_first() {
        let g = Grid::new([13, 11, 9], [10.0, 7.0, 5.0]).with_halo(4);
        let model = Model::homogeneous(g, 3000.0, 1700.0, 2500.0);
        let s = Solver::new(model, 1e-4, 0, true, &pml()).with_stencil(Stencil::taylor(8).unwrap());
        let wave = |seed: f32| -> Vec<f32> {
            (0..g.padded_len())
                .map(|p| (p as f32 * 0.731 + seed).sin() * 1e-3)
                .collect()
        };
        let v = [wave(0.1), wave(1.7), wave(2.9)];
        let v = [&v[0][..], &v[1][..], &v[2][..]];
        let start: [Vec<f32>; 6] = std::array::from_fn(|c| wave(5.0 + c as f32));
        let mut strain: [Vec<f32>; 6] = std::array::from_fn(|_| vec![0.0; g.padded_len()]);
        s.strain_rate(v, &mut strain);
        let mut added = start.clone();
        s.add_strain_rate(v, &mut added);
        for c in 0..6 {
            for k in 0..9 {
                for j in 0..11 {
                    for i in 0..13 {
                        let p = g.idx(i, j, k);
                        let want = start[c][p] + strain[c][p];
                        assert_eq!(added[c][p].to_bits(), want.to_bits(), "{c} at {i} {j} {k}");
                    }
                }
            }
        }
    }

    /// A linear velocity field has the same strain rate everywhere the stencil stays inside
    /// the grid; on the free surface ∂z vz is left out.
    #[test]
    fn strain_rates_of_a_linear_field() {
        let g = Grid::new([16, 14, 12], [10.0, 5.0, 4.0]).with_halo(4);
        let model = Model::homogeneous(g, 3000.0, 1700.0, 2500.0);
        let s = Solver::new(model, 1e-4, 0, true, &pml()).with_stencil(Stencil::taylor(8).unwrap());
        let [px, py, pz] = g.padded();
        // Gradient of each component along x, y, z (per metre), on the staggered positions.
        let grad = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let shift = [[0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5]];
        let mut v: [Vec<f32>; 3] = std::array::from_fn(|_| vec![0.0; g.padded_len()]);
        for (c, field) in v.iter_mut().enumerate() {
            for kk in 0..pz {
                for jj in 0..py {
                    for ii in 0..px {
                        let pos = [ii, jj, kk].map(|p| p as f64 - g.halo as f64);
                        let x: f64 = (0..3)
                            .map(|a| grad[c][a] * (pos[a] + shift[c][a]) * g.h[a])
                            .sum();
                        field[(kk * py + jj) * px + ii] = x as f32;
                    }
                }
            }
        }
        let mut out: [Vec<f32>; 6] = std::array::from_fn(|_| vec![0.0; g.padded_len()]);
        s.strain_rate([&v[0], &v[1], &v[2]], &mut out);
        let want = [1.0, 5.0, 9.0, 2.0 + 4.0, 3.0 + 7.0, 6.0 + 8.0];
        for k in [0, 5] {
            let p = g.idx(8, 7, k);
            for (c, &w) in want.iter().enumerate() {
                let w = if c == 2 && k == 0 { 0.0 } else { w };
                assert!(
                    (out[c][p] - w).abs() < 1e-3,
                    "component {c} at k = {k}: {}",
                    out[c][p]
                );
            }
        }
        // The halo stays untouched.
        assert_eq!(out[0][0], 0.0);
    }

    #[test]
    fn pressure_reads_the_mean_normal_stress() {
        let g = Grid::new([10, 10, 10], [10.0; 3]);
        let s = solver(g, 3, false);
        let op = s.receiver_op(Component::P, [4, 5, 6]);
        let p = g.idx(4, 5, 6);
        assert_eq!(
            op.terms,
            vec![
                (Field::Sxx, p, -1.0 / 3.0),
                (Field::Syy, p, -1.0 / 3.0),
                (Field::Szz, p, -1.0 / 3.0)
            ]
        );
    }
}
