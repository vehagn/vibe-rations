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
//! t to t + dt using the new velocity. The work is split across threads by z-plane, and each
//! plane is processed a row (constant j) at a time: the nine derivatives a row needs go into
//! scratch rows, the CPML corrects those in a layer, and the fields are updated from them.
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

use rayon::prelude::*;

use crate::cpml::{CpmlAxis, PmlSettings, Stagger};
use crate::fd::Stencil;
use crate::grid::Grid;
use crate::model::Model;

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
        let f = &self.fields;
        let fields: usize = [
            &f.vx, &f.vy, &f.vz, &f.sxx, &f.syy, &f.szz, &f.sxy, &f.sxz, &f.syz,
        ]
        .iter()
        .map(|v| v.len())
        .sum();
        let psi: usize = self.psi_v.iter().chain(&self.psi_s).map(Vec::len).sum();
        4 * (fields + psi)
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
    /// velocity (stress update). See [`planes`] for the layout.
    psi_v: [Vec<f32>; 3],
    psi_s: [Vec<f32>; 3],
    /// The spatial derivative operator.
    pub stencil: Stencil,
    /// Factors along x, y and z for the plain damping inside the side layers, if any.
    sponge: Option<[Vec<f32>; 3]>,
    /// Steps taken.
    pub step: usize,
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
            step: 0,
        }
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
            step: 0,
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
        let f = &self.fields;
        let m = &self.material;
        let d = &self.model;
        let arrays = [
            &f.vx, &f.vy, &f.vz, &f.sxx, &f.syy, &f.szz, &f.sxy, &f.sxz, &f.syz, &m.bx, &m.by,
            &m.bz, &m.lam, &m.l2m, &m.mxy, &m.mxz, &m.myz, &d.vp, &d.vs, &d.rho,
        ];
        let psi = self.psi_v.iter().chain(&self.psi_s);
        let sponge = self.sponge.iter().flatten();
        4 * arrays
            .into_iter()
            .chain(psi)
            .chain(sponge)
            .map(Vec::len)
            .sum::<usize>()
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
        self.update_velocity();
        force(&mut self.fields, t);
        if self.free_surface {
            self.mirror_velocity();
        }
        self.update_stress();
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

    fn update_velocity(&mut self) {
        let g = self.grid;
        let [nx, ny, _] = g.n;
        let [_, sy, sz] = g.strides();
        let nxt = g.padded()[0];
        let ih = g.h.map(|h| (1.0 / h) as f32);
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
        let (sxx, syy, szz, sxy, sxz, syz) =
            (&sxx[..], &syy[..], &szz[..], &sxy[..], &sxz[..], &syz[..]);
        let m = &self.material;
        let st = &self.stencil;
        let halo = g.halo;
        let [cx, cy, cz] = &self.cpml;
        let wx = cx.width();
        planes(&g, [vx, vy, vz], &mut self.psi_v, &self.cpml)
            .into_par_iter()
            .for_each_init(
                || Rows::new(nx),
                |rows, plane| {
                    let Plane {
                        kz,
                        out: [vx, vy, vz],
                        mut psi_x,
                        mut psi_y,
                        mut psi_z,
                    } = plane;
                    let [d0, d1, d2, d3, d4, d5, d6, d7, d8] = &mut rows.d;
                    for j in 0..ny {
                        let r = g.idx(0, j, kz);
                        st.forward(d0, sxx, r, 1, ih[0]); // ∂x σxx at i+½
                        st.backward(d1, sxy, r, sy, ih[1]); // ∂y σxy at j
                        st.backward(d2, sxz, r, sz, ih[2]); // ∂z σxz at k
                        st.backward(d3, sxy, r, 1, ih[0]); // ∂x σxy at i
                        st.forward(d4, syy, r, sy, ih[1]); // ∂y σyy at j+½
                        st.backward(d5, syz, r, sz, ih[2]); // ∂z σyz at k
                        st.backward(d6, sxz, r, 1, ih[0]); // ∂x σxz at i
                        st.backward(d7, syz, r, sy, ih[1]); // ∂y σyz at j
                        st.forward(d8, szz, r, sz, ih[2]); // ∂z σzz at k+½

                        if let Some(p) = psi_x.as_deref_mut() {
                            let [p0, p1, p2] = split3(&mut p[3 * wx * j..3 * wx * (j + 1)]);
                            cx.apply_along(d0, p0, Stagger::Half);
                            cx.apply_along(d3, p1, Stagger::Node);
                            cx.apply_along(d6, p2, Stagger::Node);
                        }
                        if let (Some(p), Some(s)) = (psi_y.as_deref_mut(), cy.slab(j)) {
                            let [p0, p1, p2] = split3(&mut p[3 * nx * s..3 * nx * (s + 1)]);
                            cy.apply_across(s, d1, p0, Stagger::Node);
                            cy.apply_across(s, d4, p1, Stagger::Half);
                            cy.apply_across(s, d7, p2, Stagger::Node);
                        }
                        if let (Some(p), Some(s)) = (psi_z.as_deref_mut(), cz.slab(kz)) {
                            let [p0, p1, p2] = split3(&mut p[3 * nx * j..3 * nx * (j + 1)]);
                            cz.apply_across(s, d2, p0, Stagger::Node);
                            cz.apply_across(s, d5, p1, Stagger::Node);
                            cz.apply_across(s, d8, p2, Stagger::Half);
                        }

                        let o = (j + halo) * nxt + halo;
                        add_sum(&mut vx[o..o + nx], &m.bx[r..r + nx], d0, d1, d2);
                        add_sum(&mut vy[o..o + nx], &m.by[r..r + nx], d3, d4, d5);
                        add_sum(&mut vz[o..o + nx], &m.bz[r..r + nx], d6, d7, d8);
                    }
                },
            );
    }

    fn update_stress(&mut self) {
        let g = self.grid;
        let [nx, ny, _] = g.n;
        let [_, sy, sz] = g.strides();
        let nxt = g.padded()[0];
        let ih = g.h.map(|h| (1.0 / h) as f32);
        let free_surface = self.free_surface;
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
        let (vx, vy, vz) = (&vx[..], &vy[..], &vz[..]);
        let m = &self.material;
        let st = &self.stencil;
        let halo = g.halo;
        let [cx, cy, cz] = &self.cpml;
        let wx = cx.width();
        planes(
            &g,
            [sxx, syy, szz, sxy, sxz, syz],
            &mut self.psi_s,
            &self.cpml,
        )
        .into_par_iter()
        .for_each_init(
            || Rows::new(nx),
            |rows, plane| {
                let Plane {
                    kz,
                    out: [sxx, syy, szz, sxy, sxz, syz],
                    mut psi_x,
                    mut psi_y,
                    mut psi_z,
                } = plane;
                let surface = free_surface && kz == 0;
                let [d0, d1, d2, d3, d4, d5, d6, d7, d8] = &mut rows.d;
                for j in 0..ny {
                    let r = g.idx(0, j, kz);
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

                    if let Some(p) = psi_x.as_deref_mut() {
                        let [p0, p1, p2] = split3(&mut p[3 * wx * j..3 * wx * (j + 1)]);
                        cx.apply_along(d0, p0, Stagger::Node);
                        cx.apply_along(d4, p1, Stagger::Half);
                        cx.apply_along(d6, p2, Stagger::Half);
                    }
                    if let (Some(p), Some(s)) = (psi_y.as_deref_mut(), cy.slab(j)) {
                        let [p0, p1, p2] = split3(&mut p[3 * nx * s..3 * nx * (s + 1)]);
                        cy.apply_across(s, d1, p0, Stagger::Node);
                        cy.apply_across(s, d3, p1, Stagger::Half);
                        cy.apply_across(s, d8, p2, Stagger::Half);
                    }
                    if let (Some(p), Some(s)) = (psi_z.as_deref_mut(), cz.slab(kz)) {
                        let [p0, p1, p2] = split3(&mut p[3 * nx * j..3 * nx * (j + 1)]);
                        cz.apply_across(s, d2, p0, Stagger::Node);
                        cz.apply_across(s, d5, p1, Stagger::Half);
                        cz.apply_across(s, d7, p2, Stagger::Half);
                    }

                    let o = (j + halo) * nxt + halo;
                    let lam = &m.lam[r..r + nx];
                    let l2m = &m.l2m[r..r + nx];
                    let (sxx, syy, szz) = (
                        &mut sxx[o..o + nx],
                        &mut syy[o..o + nx],
                        &mut szz[o..o + nx],
                    );
                    if surface {
                        for i in 0..nx {
                            // With σzz = 0, ∂z vz = -λ/(λ+2μ) (∂x vx + ∂y vy).
                            let (l, l2) = (lam[i], l2m[i]);
                            let q = if l2 > 0.0 { l * l / l2 } else { 0.0 };
                            let (a, b) = (l2 - q, l - q);
                            sxx[i] += a * d0[i] + b * d1[i];
                            syy[i] += b * d0[i] + a * d1[i];
                            szz[i] = 0.0;
                        }
                    } else {
                        for i in 0..nx {
                            let (l, l2) = (lam[i], l2m[i]);
                            sxx[i] += l2 * d0[i] + l * (d1[i] + d2[i]);
                            syy[i] += l2 * d1[i] + l * (d0[i] + d2[i]);
                            szz[i] += l2 * d2[i] + l * (d0[i] + d1[i]);
                        }
                    }
                    add_sum2(&mut sxy[o..o + nx], &m.mxy[r..r + nx], d3, d4);
                    add_sum2(&mut sxz[o..o + nx], &m.mxz[r..r + nx], d5, d6);
                    add_sum2(&mut syz[o..o + nx], &m.myz[r..r + nx], d7, d8);
                }
            },
        );
    }

    /// The strain rates the stress update applies the stiffness to, from the velocities
    /// `v` (vx, vy, vz, padded and mirrored as the solver leaves them), into `out` (padded,
    /// interior only): ∂x vx, ∂y vy, ∂z vz, ∂y vx + ∂x vy, ∂z vx + ∂x vz and ∂z vy + ∂y vz,
    /// each at the position of the stress it changes. On the free surface ∂z vz is left at
    /// zero, since the update replaces it there. Without the CPML: inside the absorbing
    /// layers these are not what the update uses.
    pub fn strain_rate(&self, v: [&[f32]; 3], out: &mut [Vec<f32>; 6]) {
        let g = self.grid;
        let [nx, ny, _] = g.n;
        let [_, sy, sz] = g.strides();
        let nxt = g.padded()[0];
        let ih = g.h.map(|h| (1.0 / h) as f32);
        let [vx, vy, vz] = v;
        let st = &self.stencil;
        let halo = g.halo;
        let free_surface = self.free_surface;
        let [exx, eyy, ezz, exy, exz, eyz] = out;
        let mut none = [Vec::new(), Vec::new(), Vec::new()];
        planes(&g, [exx, eyy, ezz, exy, exz, eyz], &mut none, &self.cpml)
            .into_par_iter()
            .for_each_init(
                || Rows::new(nx),
                |rows, plane| {
                    let kz = plane.kz;
                    let [exx, eyy, ezz, exy, exz, eyz] = plane.out;
                    let [d0, d1, d2, d3, d4, d5, d6, d7, d8] = &mut rows.d;
                    for j in 0..ny {
                        let r = g.idx(0, j, kz);
                        let o = (j + halo) * nxt + halo;
                        st.backward(d0, vx, r, 1, ih[0]);
                        st.backward(d1, vy, r, sy, ih[1]);
                        if free_surface && kz == 0 {
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
                        exx[o..o + nx].copy_from_slice(d0);
                        eyy[o..o + nx].copy_from_slice(d1);
                        ezz[o..o + nx].copy_from_slice(d2);
                        for i in 0..nx {
                            exy[o + i] = d3[i] + d4[i];
                            exz[o + i] = d5[i] + d6[i];
                            eyz[o + i] = d7[i] + d8[i];
                        }
                    }
                },
            );
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
        let f = &mut self.fields;
        for v in [
            &mut f.vx, &mut f.vy, &mut f.vz, &mut f.sxx, &mut f.syy, &mut f.szz, &mut f.sxy,
            &mut f.sxz, &mut f.syz,
        ] {
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
        let raw: Vec<(isize, f32)> = if m == 0 {
            vec![(0, 1.5), (1, -0.5)]
        } else {
            let w: Vec<f32> = if m == self.stencil.half_length() {
                self.stencil
                    .interpolation()
                    .iter()
                    .map(|&w| w as f32)
                    .collect()
            } else {
                crate::fd::midpoint_weights(m)
                    .iter()
                    .map(|&w| w as f32)
                    .collect()
            };
            let m = m as isize;
            (-m..m)
                .map(|o| (o, w[if o < 0 { -o - 1 } else { o } as usize]))
                .collect()
        };
        let at = at as isize;
        let mut out: Vec<(usize, f32)> = Vec::with_capacity(raw.len());
        for (offset, w) in raw {
            let m = at + offset;
            // Index -1 - m mirrors index m across the surface (positions -m - ½ and m + ½).
            let (m, w) = if m < 0 { (-1 - m, parity * w) } else { (m, w) };
            match out.iter_mut().find(|(i, _)| *i == m as usize) {
                Some((_, sum)) => *sum += w,
                None => out.push((m as usize, w)),
            }
        }
        out.retain(|&(_, w)| w != 0.0);
        out
    }

    /// The linear functional that reads component `c` at `node`.
    pub fn receiver_op(&self, c: crate::receiver::Component, node: [usize; 3]) -> PointOp {
        use crate::receiver::Component;
        let g = self.grid;
        let staggered = |axis: usize, field: Field| PointOp {
            terms: self
                .stagger_weights(axis, node[axis], 1.0)
                .into_iter()
                .map(|(m, w)| {
                    let mut n = node;
                    n[axis] = m;
                    (field, g.idx(n[0], n[1], n[2]), w)
                })
                .collect(),
        };
        match c {
            Component::Vx => staggered(0, Field::Vx),
            Component::Vy => staggered(1, Field::Vy),
            Component::Vz => staggered(2, Field::Vz),
            Component::P => {
                let p = g.idx(node[0], node[1], node[2]);
                PointOp {
                    terms: [Field::Sxx, Field::Syy, Field::Szz]
                        .map(|f| (f, p, -1.0 / 3.0))
                        .to_vec(),
                }
            }
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
    pub fn adjoint_op(&self, c: crate::receiver::Component, node: [usize; 3]) -> PointOp {
        use crate::receiver::Component;
        let g = self.grid;
        let volume = g.cell_volume() as f32;
        let plane = g.plane_len();
        let m = &self.material;
        let mut op = self.receiver_op(c, node);
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
                let p = g.idx(node[0], node[1], node[2]);
                let (l, l2) = (m.lam[p], m.l2m[p]);
                // The row sums of the update's stiffness: what the stress update does with an
                // isotropic strain, as in update_stress.
                let (horizontal, vertical) = if self.free_surface && node[2] == 0 {
                    let q = if l2 > 0.0 { l * l / l2 } else { 0.0 };
                    ((l2 - q) + (l - q), 0.0)
                } else {
                    (l2 + 2.0 * l, l2 + 2.0 * l)
                };
                let scale = 1.0 / (3.0 * volume * self.cell_share(Field::Sxx, node[2]));
                op.terms = [
                    (Field::Sxx, p, horizontal * scale),
                    (Field::Syy, p, horizontal * scale),
                    (Field::Szz, p, vertical * scale),
                ]
                .into_iter()
                .filter(|t| t.2 != 0.0)
                .collect();
            }
        }
        op
    }

    /// A component's value at a node now.
    pub fn sample(&self, c: crate::receiver::Component, node: [usize; 3]) -> f32 {
        self.receiver_op(c, node).sample(&self.fields)
    }

    /// The largest particle-velocity component anywhere, or infinity if any is not finite.
    pub fn peak_velocity(&self) -> f32 {
        let f = &self.fields;
        [&f.vx, &f.vy, &f.vz]
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

/// out += b · (x + y + z)
#[inline]
fn add_sum(out: &mut [f32], b: &[f32], x: &[f32], y: &[f32], z: &[f32]) {
    let n = out.len();
    let (b, x, y, z) = (&b[..n], &x[..n], &y[..n], &z[..n]);
    for i in 0..n {
        out[i] += b[i] * (x[i] + y[i] + z[i]);
    }
}

/// out += m · (x + y)
#[inline]
fn add_sum2(out: &mut [f32], m: &[f32], x: &[f32], y: &[f32]) {
    let n = out.len();
    let (m, x, y) = (&m[..n], &x[..n], &y[..n]);
    for i in 0..n {
        out[i] += m[i] * (x[i] + y[i]);
    }
}

fn split3(v: &mut [f32]) -> [&mut [f32]; 3] {
    let n = v.len() / 3;
    let (a, rest) = v.split_at_mut(n);
    let (b, c) = rest.split_at_mut(n);
    [a, b, c]
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
}

/// What one thread may write while it updates interior z-plane `kz`.
struct Plane<'a, const N: usize> {
    kz: usize,
    /// The plane of each output field, padded.
    out: [&'a mut [f32]; N],
    /// CPML memory for the plane: along x, rows j of [3 terms][slab]; along y, slab rows of
    /// [3 terms][nx]; along z (planes in the slab only), rows j of [3 terms][nx].
    psi_x: Option<&'a mut [f32]>,
    psi_y: Option<&'a mut [f32]>,
    psi_z: Option<&'a mut [f32]>,
}

/// Splits the output fields and CPML memory into one disjoint [`Plane`] per interior z-plane.
fn planes<'a, const N: usize>(
    g: &Grid,
    out: [&'a mut Vec<f32>; N],
    psi: &'a mut [Vec<f32>; 3],
    cpml: &[CpmlAxis; 3],
) -> Vec<Plane<'a, N>> {
    let [nx, ny, nz] = g.n;
    let len = g.plane_len();
    let halo = g.halo;
    let mut out = out.map(|v| v[halo * len..(halo + nz) * len].chunks_mut(len));
    let [px, py, pz] = psi;
    let mut px = (!px.is_empty()).then(move || px.chunks_mut(ny * 3 * cpml[0].width()));
    let mut py = (!py.is_empty()).then(move || py.chunks_mut(cpml[1].width() * 3 * nx));
    let mut pz = (!pz.is_empty()).then(move || pz.chunks_mut(ny * 3 * nx));
    (0..nz)
        .map(|kz| Plane {
            kz,
            out: std::array::from_fn(|f| out[f].next().expect("one chunk per plane")),
            psi_x: px.as_mut().and_then(Iterator::next),
            psi_y: py.as_mut().and_then(Iterator::next),
            psi_z: cpml[2]
                .slab(kz)
                .and_then(|_| pz.as_mut().and_then(Iterator::next)),
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pml() -> PmlSettings {
        PmlSettings {
            power: 2.0,
            reflection: None,
            kappa_max: 1.0,
            alpha_max: 30.0,
            vp_max: 3000.0,
            sponge: 0.0,
        }
    }

    #[test]
    fn harmonic_mean_vanishes_next_to_a_fluid() {
        assert_eq!(harmonic_mean([1.0, 1.0, 1.0, 0.0]), 0.0);
        assert!((harmonic_mean([2.0, 2.0, 2.0, 2.0]) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn a_wavefield_at_rest_stays_at_rest() {
        let g = Grid::new([16, 14, 12], [10.0; 3]);
        let mut s = Solver::new(
            Model::homogeneous(g, 3000.0, 1700.0, 2500.0),
            1e-3,
            4,
            true,
            &pml(),
        );
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
        let mut s = Solver::new(
            Model::homogeneous(g, 3000.0, 1700.0, 2500.0),
            1e-3,
            3,
            true,
            &pml(),
        );
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
        let mut s = Solver::new(
            Model::homogeneous(g, 3000.0, 1700.0, 2500.0),
            1e-3,
            3,
            true,
            &pml(),
        );
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
        let model = Model::homogeneous(g, 3000.0, 1700.0, 2500.0);
        let open = Solver::new(model.clone(), 1e-3, 3, false, &pml());
        let free = Solver::new(model, 1e-3, 3, true, &pml());
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

    /// The interpolation is as long as the stencil, and exact for polynomials up to one
    /// degree less than its length.
    #[test]
    fn interpolation_matches_the_stencil() {
        let g = Grid::new([20, 20, 20], [10.0; 3]).with_halo(4);
        let s = Solver::new(
            Model::homogeneous(g, 3000.0, 1700.0, 2500.0),
            1e-3,
            3,
            false,
            &pml(),
        )
        .with_stencil(Stencil::taylor(8).unwrap());
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

    #[test]
    fn splits_every_plane_once() {
        let g = Grid::new([9, 8, 12], [1.0; 3]);
        let cpml: [CpmlAxis; 3] = std::array::from_fn(|a| {
            CpmlAxis::new(g.n[a], if a == 2 { 0 } else { 3 }, 3, 1.0, 1e-3, &pml())
        });
        let mut a = vec![0.0; g.padded_len()];
        let mut psi = [
            vec![0.0; 12 * 8 * 3 * 6],
            vec![0.0; 12 * 6 * 3 * 9],
            vec![0.0; 3 * 8 * 3 * 9],
        ];
        let p = planes(&g, [&mut a], &mut psi, &cpml);
        assert_eq!(p.len(), 12);
        assert!(p.iter().all(|p| p.psi_x.is_some() && p.psi_y.is_some()));
        let with_z: Vec<usize> = p
            .iter()
            .filter(|p| p.psi_z.is_some())
            .map(|p| p.kz)
            .collect();
        assert_eq!(with_z, vec![9, 10, 11]);
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
        let s = Solver::new(
            Model::homogeneous(g, 3000.0, 1700.0, 2500.0),
            1e-3,
            3,
            true,
            &pml(),
        );
        let mut probe = Fields::zeros(g.padded_len());
        for (n, v) in probe.vz.iter_mut().enumerate() {
            *v = ((n * 7919) % 1000) as f32 / 1000.0 - 0.5;
        }
        probe.vx.clone_from(&probe.vz);
        probe.sxx.clone_from(&probe.vz);
        for (c, node) in [
            (crate::receiver::Component::Vz, [5, 5, 0]),
            (crate::receiver::Component::Vz, [5, 5, 1]),
            (crate::receiver::Component::Vx, [0, 5, 5]),
            (crate::receiver::Component::P, [5, 5, 5]),
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
        let model = Model::homogeneous(g, 3000.0, 1700.0, 2500.0);
        let free = Solver::new(model.clone(), 1e-3, 3, true, &pml());
        let open = Solver::new(model, 1e-3, 3, false, &pml());
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
        let s = Solver::new(
            Model::homogeneous(g, 3000.0, 1700.0, 2500.0),
            1e-3,
            3,
            false,
            &pml(),
        );
        let op = s.receiver_op(crate::receiver::Component::P, [4, 5, 6]);
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
