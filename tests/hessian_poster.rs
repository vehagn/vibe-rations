//! The Hessian matrix along a column under the source, parameter block by parameter block:
//! the object of Hagen, Arntsen & Raknes's poster (SEG 2018), built from this solver's exact
//! Hessian instead of differenced wavefields, and measured against closed forms.
//!
//! `cargo test --release --test hessian_poster -- --ignored --nocapture --exact <name>`
//!
//! The poster's simple model, as its figure files and the code behind them give it: a
//! homogeneous full space (ρ 1000 kg/m³, vp 1500 m/s, vs 1000 m/s; absorbing layers all
//! round, no free surface), a box 1000 m wide and deep, a pressure source with a 10 Hz
//! Ricker wavelet 30 m below the box's top at its middle, pressure receivers at every node
//! of that depth across the box, the L2 misfit, a record of one second, and observed data
//! from a model with a smooth disc at the box's centre (radius 200 m, +100 in each of ρ, vp
//! and vs). Its Hessian matrix has one column per grid cell of the vertical line under the
//! source and per parameter: the change of the gradient on that line when that cell's ρ, vp
//! or vs changes.
//!
//! Two versions. In 2D (`poster_2d…`, a grid one node across) it is the poster's own
//! survey, a line source, like for like with the code behind it. In 3D the disc is a cylinder
//! across the plane, as a 2D model means it, the source a point and the receivers a line in
//! its plane, which is what the point scatterer's closed forms are for. The matrix is
//! [`vibe_rations::hessian::hessian_block`]'s: one forward run, one adjoint run and one Born
//! run per column give JᵀJ (the Born traces' inner products), B (each Born run against the
//! adjoint run) and the local term, so H = JᵀJ + B + Bᵀ + local has its second-order part at
//! its own scale, not as a difference of numbers thousands of times larger.
//!
//! - `coarse`, `poster`, `poster_long`, `poster_2d`, `poster_2d_long`: the matrices into
//!   `out/hessian-poster/<case>/block`, with the closed forms: the start's traces against an
//!   explosion's pressure (in 2D a line of them), and in 3D every Born trace against a point
//!   scatterer's (`common::explosion`) and JᵀJ against the same from those.
//! - `poster_2d_floor`: the 2D survey over a reflecting half-space 20 m under the box,
//!   recorded until its reflection is back. The disc's edge is smooth (30 m against a P
//!   wavelength of 150 m), so its own reflection is weak and the plain survey's residual
//!   small; a wave returned from below crosses the disc twice and its residual, a third of
//!   a period's advance, is 10⁴ times larger in energy. The gradient on the column then has
//!   one sign at every depth, and so has the second-order part's diagonal: the look of the
//!   poster's own figure, which the plain survey's exact matrix hasn't.
//! - `coarse_products`, `poster_products`, `poster_2d_products`: six of those columns
//!   against whole products ([`hessian_product`]: a Born and two adjoint runs each), and the
//!   kernels in the plane.
//! - `coarse_differences`, `poster_2d_differences`: the same six against differences of the
//!   gradient.
//! - `coarse_along`, `poster_along`, `poster_2d_along`: the misfit along directions on the
//!   column read from files: the Hessian's curvature and a step's worth, measured without
//!   the Hessian.
//! - `coarse_without_residual`: with the start's own data the second-order part is zero and
//!   a product is JᵀJ's column.
//! - `coarse_forward`, `poster_forward`: the start's traces alone, for choosing the margins.
//! - `trials_…`: the 2D matrix on a column of fixed depths with one thing changed at a time
//!   (the spacing, the stencil and time step, the layers and the record, the wavelet's
//!   frequency, the walls' place, a small inclusion in place of the disc), into
//!   `out/hessian-poster/trials/<name>`: which of its features are the survey's and which
//!   the grid's or the boundary's.
//! - `trials_differenced`: Born traces and columns from differences of whole runs with the
//!   steps the code behind the poster took, against the exact ones.

mod common;

use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

use rayon::prelude::*;
use vibe_rations::born::Change;
use vibe_rations::config::{Config, Parameter};
use vibe_rations::gradient::gradient;
use vibe_rations::hessian::{Block, hessian_block, hessian_product};
use vibe_rations::inversion::{Inversion, Observed};
use vibe_rations::model::Model;
use vibe_rations::npy;
use vibe_rations::simulation::Simulation;

use common::explosion::{Explosion, evaluate, velocity_scattering};
use common::line::along_line;
use common::stokes::Medium;
use common::{norm, sub};

/// The box's width and depth, m.
const BOX: f64 = 1000.0;
/// The poster's start: vp, vs (m/s) and ρ (kg/m³).
const START: [f64; 3] = [1500.0, 1000.0, 1000.0];
/// The anomaly the observed data hold: its radius, the width of its edge (the standard
/// deviation of the Gaussian a sharp edge is smoothed by) and its size in each parameter.
const RADIUS: f64 = 200.0;
const EDGE: f64 = 30.0;
const ANOMALY: f64 = 100.0;
/// How far under the box a case's floor is, m.
const FLOOR_BELOW: f64 = 20.0;
/// The explosion's moment rate, N·m/s.
const RATE: f64 = 1e9;
/// The blocks' order, as on the poster.
const ORDER: [Parameter; 3] = [Parameter::Rho, Parameter::Vp, Parameter::Vs];

#[derive(Clone, Copy, Debug)]
struct Case {
    name: &'static str,
    spacing: f64,
    frequency: f64,
    delay: f64,
    /// The source's and the receivers' depth below the box's top.
    depth: f64,
    /// The column's first depth; it runs to the box's bottom, node by node.
    column_top: f64,
    /// The first depth whose columns are compared with the closed forms: nearer the source
    /// and the receivers a point scatterer's field isn't one the grid can hold.
    exact_from: f64,
    /// The absorbing layers' width, nodes.
    pml: usize,
    /// Nodes between the box and the absorbing layers, on its four sides in the plane.
    around: usize,
    /// Nodes outside the layers on each side of the plane.
    across: usize,
    duration: f64,
    /// A 2D run: the grid one node across, the source a line, as the poster's own code had
    /// it. Otherwise 3D with a point source, which the closed forms are for.
    planar: bool,
    /// A reflector under the box, 20 m below its bottom: vp, vs and ρ of a half-space
    /// there, in the start and in the observed data's model alike.
    floor: Option<[f64; 3]>,
    /// How far to the side of the column the source is, m.
    aside: f64,
    /// The column's last depth and its step, m (0: every node), and the receivers' step
    /// (0: every node): set in metres, grids of different spacings share them.
    column_bottom: f64,
    column_step: f64,
    receiver_step: f64,
    /// What the observed data's model adds to the start.
    anomaly: Anomaly,
    /// The start: vp, vs (m/s) and ρ (kg/m³).
    start: [f64; 3],
    /// The stencil's order, and the time step as a fraction of the stability limit.
    order: usize,
    courant: f64,
    /// The time stepping of the code behind the poster instead: steps of this length with
    /// Taylor's coefficients and leapfrog's time error left in the traces.
    plain_steps: Option<f64>,
}

/// The anomaly in the observed data's model, +[`ANOMALY`] in each of ρ, vp and vs.
#[derive(Clone, Copy, Debug)]
enum Anomaly {
    /// The poster's disc at the box's centre.
    Disc,
    /// A Gaussian of this standard deviation at this depth under the source: a small
    /// inclusion (the EGU talk's was 30 m square: the area of one of 12 m), smooth so that
    /// every grid samples the same one.
    Blob { depth: f64, sigma: f64 },
}

/// The poster's model at half its frequency on twice its spacing, with the same number of
/// nodes per wavelength: the same study at a thirtieth of the cost.
const COARSE: Case = Case {
    name: "coarse",
    spacing: 20.0,
    frequency: 5.0,
    delay: 0.3,
    depth: 40.0,
    column_top: 80.0,
    exact_from: 80.0,
    pml: 12,
    around: 4,
    across: 6,
    duration: 2.0,
    planar: false,
    floor: None,
    aside: 0.0,
    column_bottom: BOX,
    column_step: 0.0,
    receiver_step: 0.0,
    anomaly: Anomaly::Disc,
    start: START,
    order: 8,
    courant: 0.85,
    plain_steps: None,
};

/// The poster's spacing, wavelet, depths and record length (1000 steps of 1 ms).
const POSTER: Case = Case {
    name: "poster",
    spacing: 10.0,
    frequency: 10.0,
    delay: 0.15,
    depth: 30.0,
    column_top: 0.0,
    exact_from: 50.0,
    pml: 12,
    around: 4,
    across: 6,
    duration: 1.0,
    planar: false,
    floor: None,
    aside: 0.0,
    column_bottom: BOX,
    column_step: 0.0,
    receiver_step: 0.0,
    anomaly: Anomaly::Disc,
    start: START,
    order: 8,
    courant: 0.85,
    plain_steps: None,
};

/// The poster's survey recorded for long enough to hold the reflections of the whole
/// column: the poster's own second ends before those from below about 600 m arrive.
const POSTER_LONG: Case = Case {
    name: "poster-long",
    duration: 1.65,
    ..POSTER
};

/// The poster's survey in 2D, as its own code ran it: a line source, and the disc a disc.
const POSTER_2D: Case = Case {
    name: "poster-2d",
    planar: true,
    ..POSTER
};

const POSTER_2D_LONG: Case = Case {
    name: "poster-2d-long",
    planar: true,
    ..POSTER_LONG
};

/// The same over a reflecting floor, recorded until its reflection has come back: the wave
/// then crosses the disc twice, and the residual holds what the disc did to it on the way.
const POSTER_2D_FLOOR: Case = Case {
    name: "poster-2d-floor",
    floor: Some([3000.0, 2000.0, 2000.0]),
    ..POSTER_2D_LONG
};

/// The poster's survey as the code behind its figure ran it. Its absorbing layers, ten cells
/// under a stencil sixteen long, had two cells at work, at under 2% of their damping: the
/// model was a box with walls two cells outside it, and the record, 2000 steps of 1 ms, long
/// enough for the wave to cross it three times. The source was one cell to the side of the
/// column.
const POSTER_2D_BOX: Case = Case {
    name: "poster-2d-box",
    pml: 0,
    around: 2,
    duration: 2.0,
    aside: 10.0,
    ..POSTER_2D
};

impl Case {
    /// Where the box's point `x` across and `z` down is on the grid, in the source's plane.
    fn at(&self, x: f64, z: f64) -> [f64; 3] {
        let h = self.spacing;
        let edge = (self.pml + self.around) as f64 * h;
        [edge + x, self.row() as f64 * h, edge + z]
    }

    /// The plane of the source: its node along y.
    fn row(&self) -> usize {
        if self.planar {
            0
        } else {
            self.pml + self.across
        }
    }

    fn points(&self) -> [usize; 3] {
        let n = (BOX / self.spacing).round() as usize + 1 + 2 * (self.pml + self.around);
        [n, 2 * self.row() + 1, n]
    }

    fn source(&self) -> [f64; 3] {
        self.at(BOX / 2.0 + self.aside, self.depth)
    }

    fn receivers(&self) -> Vec<[f64; 3]> {
        let step = self.step(self.receiver_step);
        let n = (BOX / step).round() as usize;
        (0..=n)
            .map(|i| self.at(i as f64 * step, self.depth))
            .collect()
    }

    /// A step in metres, or the grid's for 0.
    fn step(&self, metres: f64) -> f64 {
        if metres > 0.0 { metres } else { self.spacing }
    }

    /// The column's depths below the box's top.
    fn depths(&self) -> Vec<f64> {
        let step = self.step(self.column_step);
        let (first, last) = (self.column_top / step, self.column_bottom / step);
        (first.round() as usize..=last.round() as usize)
            .map(|k| k as f64 * step)
            .collect()
    }

    fn explosion(&self) -> Explosion {
        Explosion {
            rate: RATE,
            f0: self.frequency,
            t0: self.delay,
        }
    }

    fn survey(&self) -> String {
        let [nx, ny, nz] = self.points();
        let (h, duration) = (self.spacing, self.duration);
        let [vp, vs, rho] = self.start;
        // The observed data's model is faster than the start, whose time step the inversion
        // keeps: the usual fraction of the stability limit there.
        let fastest = self.floor.map_or(vp, |f| f[0]).max(vp + ANOMALY);
        // The layers' own Courant number follows the model's fastest node.
        let courant = self.courant * self.floor.map_or(vp, |f| f[0]).max(vp) / fastest;
        let stepping = match self.plain_steps {
            Some(dt) => format!("dt = {dt}"),
            None => format!("courant = {courant}"),
        };
        let scheme = match self.plain_steps {
            Some(_) => "coefficients = \"taylor\"\nremove_time_dispersion = false",
            None => "",
        };
        let order = self.order;
        let floor = self.floor.map_or(String::new(), |[vp, vs, rho]| {
            let top = self.at(0.0, BOX)[2] + FLOOR_BELOW;
            format!("[[model.layer]]\ntop = {top}\nvp = {vp}\nvs = {vs}\nrho = {rho}")
        });
        let source = self.source();
        let receivers = self.receivers();
        let (from, to, count) = (
            receivers[0],
            receivers[receivers.len() - 1],
            receivers.len(),
        );
        let (f0, t0, pml) = (self.frequency, self.delay, self.pml);
        format!(
            r#"
            [grid]
            points = [{nx}, {ny}, {nz}]
            spacing = {h}
            [scheme]
            order = {order}
            {scheme}
            [time]
            duration = {duration}
            {stepping}
            [boundary]
            pml_width = {pml}
            free_surface = false
            pml_vp_max = {fastest}
            [[model.layer]]
            top = 0.0
            vp = {vp}
            vs = {vs}
            rho = {rho}
            {floor}
            [[source]]
            position = {source:?}
            explosion = {RATE:e}
            wavelet = {{ kind = "ricker", frequency = {f0}, delay = {t0} }}
            [receivers]
            components = ["p"]
            [[receivers.line]]
            from = {from:?}
            to = {to:?}
            count = {count}
            [inversion]
            observed = "unused"
            parameters = ["vp", "vs", "rho"]
            iterations = 1
            balance = false
            reconstruct = false
            precondition = false
            "#
        )
    }

    fn dir(&self) -> PathBuf {
        Path::new("out/hessian-poster").join(self.name)
    }
}

/// ½ erfc(u/√2): a unit step down at 0 smoothed by a Gaussian of unit standard deviation
/// (Abramowitz & Stegun 7.1.26, to 1.5e-7).
fn smooth_step(u: f64) -> f64 {
    let x = u.abs() / std::f64::consts::SQRT_2;
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    let half_erfc = 0.5 * poly * (-x * x).exp();
    if u >= 0.0 { half_erfc } else { 1.0 - half_erfc }
}

/// The start with the anomaly added, the model the observed data come from: the poster's
/// disc as a 2D model means it, a cylinder across the plane, which the layers there continue.
fn truth(case: &Case, start: &Model) -> Model {
    let mut model = start.clone();
    let g = model.grid;
    let centre = match case.anomaly {
        Anomaly::Disc => case.at(BOX / 2.0, BOX / 2.0),
        Anomaly::Blob { depth, .. } => case.at(BOX / 2.0, depth),
    };
    let [nx, ny, nz] = g.n;
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let mut d = sub(g.position([i, j, k]), centre);
                d[1] = 0.0;
                let shape = match case.anomaly {
                    Anomaly::Disc => smooth_step((norm(d) - RADIUS) / EDGE),
                    Anomaly::Blob { sigma, .. } => (-0.5 * (norm(d) / sigma).powi(2)).exp(),
                };
                let a = (ANOMALY * shape) as f32;
                let n = g.interior_idx(i, j, k);
                model.vp[n] += a;
                model.vs[n] += a;
                model.rho[n] += a;
            }
        }
    }
    model
}

struct Survey {
    case: Case,
    /// The inversion of the observed data from the start: one shot, whose configuration
    /// (time step, layers and coefficients pinned) every run here uses.
    inv: Inversion,
    /// The column's nodes (interior indices), top to bottom.
    column: Vec<usize>,
    /// The pressure samples' times.
    times: Vec<f64>,
}

/// The survey of `case`, its observed data from the model with the anomaly, or with
/// `residual` false from the start itself, so that the residuals vanish to the last bit.
fn setup(case: Case, residual: bool) -> Survey {
    let cfg = Config::from_toml(&case.survey()).unwrap();
    let base = Path::new(".");
    let start = Model::from_config(Simulation::grid(&cfg).unwrap(), &cfg.model, base).unwrap();
    let observed_in = if residual {
        truth(&case, &start)
    } else {
        start.clone()
    };
    let settings = cfg.inversion.clone().unwrap();
    let record = |shot: &Config, model: &Model| -> Observed {
        let mut sim = Simulation::with_model(shot, model.clone()).unwrap();
        sim.run(None, false).unwrap();
        sim.receivers.traces.into()
    };
    // Twice: the first inversion only pins the shot's configuration, which then records the
    // observed traces.
    let first = vec![record(&cfg.shot(0).unwrap(), &start)];
    let pinned = Inversion::new(&cfg, start.clone(), first, &settings).unwrap();
    let observed = vec![record(&pinned.shots[0].config, &observed_in)];
    let inv = Inversion::new(&cfg, start, observed, &settings).unwrap();

    let g = inv.model.grid;
    let nearest = |p: [f64; 3]| {
        let c = g.coordinates(p).expect("on the grid");
        g.interior_idx(
            c[0].round() as usize,
            c[1].round() as usize,
            c[2].round() as usize,
        )
    };
    let column = case
        .depths()
        .iter()
        .map(|&z| nearest(case.at(BOX / 2.0, z)))
        .collect();
    let sim = Simulation::with_model(&inv.shots[0].config, inv.model.clone()).unwrap();
    // A stress sample follows its step: half a step after the velocities'.
    let times = sim
        .sample_times()
        .iter()
        .map(|t| t + 0.5 * sim.solver.dt)
        .collect();
    Survey {
        case,
        inv,
        column,
        times,
    }
}

fn node_position(model: &Model, n: usize) -> [f64; 3] {
    let g = model.grid;
    let [nx, ny, _] = g.n;
    g.position([n % nx, (n / nx) % ny, n / (nx * ny)])
}

fn write_f64(path: &Path, values: &[f64]) {
    let bytes: Vec<u8> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
    // Written beside and renamed, so a stopped run leaves no half-written column.
    let partial = path.with_extension("partial");
    fs::write(&partial, bytes).unwrap();
    fs::rename(&partial, path).unwrap();
}

fn read_f64(path: &Path) -> Vec<f64> {
    fs::read(path)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
        .as_chunks::<8>()
        .0
        .iter()
        .map(|&b| f64::from_le_bytes(b))
        .collect()
}

/// What one product leaves: the gradient's change on the column for vp, vs and ρ (in
/// [`Parameter`]'s order, a column of the Hessian), then the same of its two wavefield terms
/// alone (the Born run against the adjoint, the forward run against the adjoint's change);
/// the Born traces; and the product.
struct Column {
    kernels: Vec<f64>,
    born: Vec<f32>,
    product: vibe_rations::hessian::Product,
    speeds: [Vec<f64>; 3],
}

/// H applied to a unit change of `p` at the node `node`: 1 m/s or 1 kg/m³.
fn product(s: &Survey, p: Parameter, node: usize) -> Column {
    let shot = &s.inv.shots[0];
    let model = &s.inv.model;
    let n = model.grid.interior_len();
    let mut speeds = [vec![0.0; n], vec![0.0; n], vec![0.0; n]];
    speeds[p as usize][node] = 1.0;
    let by_ref = speeds.each_ref().map(Vec::as_slice);
    let change = Change::from_velocities(model, by_ref);
    let mut sim = Simulation::with_model(&shot.config, model.clone()).unwrap();
    let mut born = Vec::new();
    let product = hessian_product(
        &mut sim,
        &change,
        |r| shot.misfit.evaluate(&r.traces),
        |d| {
            born = d[0].clone();
            shot.misfit.gauss_newton_source(d)
        },
        None,
    )
    .unwrap();
    let [first, second] = &product.terms;
    let kernels = [
        product.velocities(model, by_ref),
        first.velocities(model),
        second.velocities(model),
    ]
    .iter()
    .flatten()
    .flat_map(|k| s.column.iter().map(|&m| k[m]))
    .collect();
    Column {
        kernels,
        born,
        product,
        speeds,
    }
}

/// The plane of the source, (nz, nx), of three fields over the interior nodes, as one array.
fn plane(model: &Model, fields: [&[f64]; 3], j: usize) -> (Vec<usize>, Vec<f32>) {
    let [nx, _, nz] = model.grid.n;
    let mut out = Vec::with_capacity(3 * nx * nz);
    for f in fields {
        for k in 0..nz {
            for i in 0..nx {
                out.push(f[model.grid.interior_idx(i, j, k)] as f32);
            }
        }
    }
    (vec![3, nz, nx], out)
}

/// ‖M - Mᵀ‖ / ‖M‖ (Frobenius) of a square matrix of side `n`.
fn asymmetry(m: &[f64], n: usize) -> f64 {
    let mut num = 0.0;
    for a in 0..n {
        for b in 0..n {
            num += (m[a * n + b] - m[b * n + a]).powi(2);
        }
    }
    (num / m.iter().map(|x| x * x).sum::<f64>()).sqrt()
}

/// Σ over the samples of every pair of trace sets' products: the L2 misfit's JᵀJ.
fn gram<T: Copy + Into<f64> + Sync>(traces: &[Vec<T>]) -> Vec<f64> {
    let side = traces.len();
    (0..side * side)
        .into_par_iter()
        .map(|i| {
            let (a, b) = (&traces[i / side], &traces[i % side]);
            a.iter()
                .zip(b)
                .map(|(&x, &y)| x.into() * y.into())
                .sum::<f64>()
        })
        .collect()
}

/// A column of the Hessian from central differences of the gradient for steps of ±`step`
/// and ±`step`/2 of `p` at the node, extrapolated (Richardson), on the column's nodes.
fn difference_column(s: &Survey, p: Parameter, node: usize, step: f64) -> Vec<f64> {
    let shot = &s.inv.shots[0];
    let at = |e: f64| -> Vec<f64> {
        let mut model = s.inv.model.clone();
        model.values_mut(p)[node] += e as f32;
        let mut sim = Simulation::with_model(&shot.config, model.clone()).unwrap();
        let run = gradient(&mut sim, |r| shot.misfit.evaluate(&r.traces), None).unwrap();
        run.gradient
            .velocities(&model)
            .iter()
            .flat_map(|k| s.column.iter().map(|&m| k[m]))
            .collect()
    };
    let central = |e: f64| -> Vec<f64> {
        let (plus, minus) = (at(e), at(-e));
        plus.iter()
            .zip(&minus)
            .map(|(a, b)| (a - b) / (2.0 * e))
            .collect()
    };
    let (whole, half) = (central(step), central(0.5 * step));
    whole
        .iter()
        .zip(&half)
        .map(|(w, h)| (4.0 * h - w) / 3.0)
        .collect()
}

fn relative(a: &[f64], b: &[f64]) -> f64 {
    let num: f64 = a.iter().zip(b).map(|(x, y)| (x - y).powi(2)).sum();
    (num / b.iter().map(|y| y * y).sum::<f64>()).sqrt()
}

/// The traces of the survey's shot in `model`.
fn record(s: &Survey, model: &Model) -> Vec<f32> {
    let mut sim = Simulation::with_model(&s.inv.shots[0].config, model.clone()).unwrap();
    sim.run(None, false).unwrap();
    sim.receivers.traces.swap_remove(0)
}

/// The start's traces against the explosion's pressure: each receiver's misfit, and the
/// traces.
fn forward_misfits(s: &Survey) -> (Vec<f64>, Vec<f32>) {
    let case = &s.case;
    let traces = record(s, &s.inv.model);
    let medium = medium(case);
    let explosion = case.explosion();
    let samples = s.times.len();
    let misfits = case
        .receivers()
        .iter()
        .zip(traces.chunks(samples))
        .map(|(&x, trace)| {
            let r = norm(sub(x, case.source()));
            if r == 0.0 {
                return f64::NAN;
            }
            let exact: Vec<f64> = s
                .times
                .iter()
                .map(|&t| {
                    if case.planar {
                        // A line of explosions: the point's pressure summed along it.
                        along_line(medium.vp, r, t, |y| {
                            explosion.pressure(&medium, r.hypot(y), t)
                        })
                    } else {
                        explosion.pressure(&medium, r, t)
                    }
                })
                .collect();
            let got: Vec<f64> = trace.iter().map(|&v| v as f64).collect();
            relative(&got, &exact)
        })
        .collect();
    (misfits, traces)
}

fn medium(case: &Case) -> Medium {
    let [vp, vs, rho] = case.start;
    Medium { vp, vs, rho }
}

fn print_forward(misfits: &[f64], h: f64) {
    println!("the start's pressure against the explosion's, by the receiver's offset:");
    let middle = misfits.len() / 2;
    for i in (0..=middle).rev().step_by((middle / 6).max(1)) {
        // A receiver on the source has no closed form.
        if misfits[i].is_finite() {
            let offset = (middle - i) as f64 * h;
            println!("  {offset:>6.0} m  {:>8.4}%", 100.0 * misfits[i]);
        }
    }
    let worst = misfits
        .iter()
        .filter(|m| m.is_finite())
        .fold(0.0f64, |a, &m| a.max(m));
    println!("  worst {:.4}%", 100.0 * worst);
}

/// The Born traces for a unit change of vp, vs and ρ ([`Parameter`]'s order, a parameter's
/// nodes together, as the block's columns are) at each of the column's nodes, from the point
/// scatterer's closed form; zeros above the case's `exact_from`.
fn exact_traces(s: &Survey) -> Vec<Vec<f64>> {
    let case = &s.case;
    let (medium, explosion) = (medium(case), case.explosion());
    let (source, receivers) = (case.source(), case.receivers());
    let samples = s.times.len();
    // The wavelet and its derivatives are below 1e-9 of their peaks this far from its centre.
    let reach = 1.6 / case.frequency;
    let depths = case.depths();
    let per_node: Vec<[Vec<f64>; 3]> = s
        .column
        .par_iter()
        .zip(&depths)
        .map(|(&node, &depth)| {
            let y = node_position(&s.inv.model, node);
            let mut out = [(); 3].map(|_| vec![0.0; receivers.len() * samples]);
            if depth < case.exact_from {
                return out;
            }
            for (r, &x) in receivers.iter().enumerate() {
                let terms = velocity_scattering(&medium, source, x, y, case.spacing);
                for (q, terms) in terms.iter().enumerate() {
                    let (first, last) = terms.iter().fold((f64::MAX, f64::MIN), |(a, b), t| {
                        (a.min(t.delay), b.max(t.delay))
                    });
                    for (m, &t) in s.times.iter().enumerate() {
                        if t - case.delay < first - reach || t - case.delay > last + reach {
                            continue;
                        }
                        out[q][r * samples + m] = evaluate(terms, &explosion, t);
                    }
                }
            }
            out
        })
        .collect();
    let mut traces = Vec::with_capacity(3 * s.column.len());
    for p in 0..3 {
        for node in &per_node {
            traces.push(node[p].clone());
        }
    }
    traces
}

/// The matrix between the column's nodes from one forward run, one adjoint run and a Born
/// run per column ([`hessian_block`]), `at_once` Born runs side by side.
fn block(s: &Survey, at_once: usize) -> Block {
    let shot = &s.inv.shots[0];
    let mut sim = Simulation::with_model(&shot.config, s.inv.model.clone()).unwrap();
    hessian_block(
        &mut sim,
        &s.column,
        |r| shot.misfit.evaluate(&r.traces),
        |d| shot.misfit.gauss_newton_source(d),
        None,
        at_once,
    )
    .unwrap()
}

/// A block's matrices in vp, vs and ρ ([`Parameter`]'s order, a parameter's nodes together;
/// row: the gradient's parameter and node, column: the changed one's): the Hessian, JᵀJ, B
/// (each Born run against the adjoint run) and the local term, the Hessian less JᵀJ, B and
/// Bᵀ.
struct BlockMatrices {
    hessian: Vec<f64>,
    gauss_newton: Vec<f64>,
    scattered: Vec<f64>,
    local: Vec<f64>,
}

impl BlockMatrices {
    fn of(s: &Survey, b: &Block) -> Self {
        let model = &s.inv.model;
        let side = b.side();
        let hessian = b.velocities(model);
        let gauss_newton = b.in_velocities(model, &b.gauss_newton);
        let scattered = b.in_velocities(model, &b.scattered);
        let local = (0..side * side)
            .map(|i| {
                let t = (i % side) * side + i / side;
                hessian[i] - gauss_newton[i] - scattered[i] - scattered[t]
            })
            .collect();
        Self {
            hessian,
            gauss_newton,
            scattered,
            local,
        }
    }

    const FILES: [&'static str; 4] = [
        "hessian.f64",
        "gauss_newton.f64",
        "scattered.f64",
        "local.f64",
    ];

    fn write(&self, dir: &Path) {
        let parts = [
            &self.hessian,
            &self.gauss_newton,
            &self.scattered,
            &self.local,
        ];
        for (name, m) in Self::FILES.iter().zip(parts) {
            write_f64(&dir.join(name), m);
        }
    }

    /// The matrices a run of the case left in `dir`.
    fn read(dir: &Path) -> Self {
        let [hessian, gauss_newton, scattered, local] =
            Self::FILES.map(|name| read_f64(&dir.join(name)));
        Self {
            hessian,
            gauss_newton,
            scattered,
            local,
        }
    }
}

/// The Born traces of the block's columns in vp, vs and ρ, in double precision.
fn born_traces(s: &Survey, b: &Block) -> Vec<Vec<f64>> {
    b.jacobian_in_velocities(&s.inv.model)
        .iter()
        .map(|components| components[0].iter().map(|&v| v as f64).collect())
        .collect()
}

/// A trace set per row, in single precision.
fn write_traces(path: &Path, traces: &[Vec<f64>]) {
    let flat: Vec<f32> = traces.iter().flatten().map(|&v| v as f32).collect();
    npy::write(path, &[traces.len(), traces[0].len()], &flat).unwrap();
}

/// The case's block with its checks against the closed forms, written into `block` under
/// the case's directory (raw little-endian doubles unless .npy): the matrices
/// ([`BlockMatrices`]), the gradient on the column and in the plane, the Born traces, the
/// closed form's (`jacobian_exact.npy`) and its JᵀJ, each Born trace's misfit against the
/// closed form, and the start's and the observed traces.
fn run(case: Case, at_once: usize) {
    let s = setup(case, true);
    let dir = case.dir().join("block");
    fs::create_dir_all(&dir).unwrap();
    let [nx, ny, nz] = case.points();
    let n = s.column.len();
    let side = 3 * n;
    println!(
        "{}: {nx} × {ny} × {nz} nodes spaced {} m, {} steps of {:.3} ms, {} receivers, {side} columns",
        case.name,
        case.spacing,
        s.times.len(),
        1e3 * (s.times[1] - s.times[0]),
        case.receivers().len(),
    );
    let started = Instant::now();
    let b = block(&s, at_once);
    println!(
        "the block between {n} nodes in {:.0} s, {:.2} GiB at most; misfit {:.6e}",
        started.elapsed().as_secs_f64(),
        b.memory_bytes as f64 / (1u64 << 30) as f64,
        b.misfit
    );
    let m = BlockMatrices::of(&s, &b);
    m.write(&dir);
    let model = &s.inv.model;
    let g = b.gradient.velocities(model);
    let on_column: Vec<f64> = g
        .iter()
        .flat_map(|g| s.column.iter().map(|&m| g[m]))
        .collect();
    write_f64(&dir.join("gradient.f64"), &on_column);
    write_f64(&dir.join("depths.f64"), &case.depths());
    write_f64(&dir.join("times.f64"), &s.times);
    let (shape, data) = plane(model, [&g[0][..], &g[1][..], &g[2][..]], case.row());
    npy::write(&dir.join("gradient.npy"), &shape, &data).unwrap();
    let traces = born_traces(&s, &b);
    write_traces(&dir.join("jacobian.npy"), &traces);

    let norm = |m: &[f64]| m.iter().map(|x| x * x).sum::<f64>().sqrt();
    let second: Vec<f64> = m
        .hessian
        .iter()
        .zip(&m.gauss_newton)
        .map(|(h, g)| h - g)
        .collect();
    println!(
        "‖JᵀJ‖ {:.3e}, ‖B‖ {:.3e}, ‖local‖ {:.3e}; the second-order part is {:.3e} of the Hessian; B's asymmetry {:.2}",
        norm(&m.gauss_newton),
        norm(&m.scattered),
        norm(&m.local),
        norm(&second) / norm(&m.hessian),
        asymmetry(&m.scattered, side)
    );

    // The closed forms.
    println!();
    let (forward, traces_start) = forward_misfits(&s);
    if case.floor.is_none() && case.pml > 0 {
        write_f64(&dir.join("forward_misfit.f64"), &forward);
        print_forward(&forward, case.step(case.receiver_step));
    } else {
        println!("over a floor or between walls the explosion's closed form doesn't hold");
    }
    let gather = [case.receivers().len(), s.times.len()];
    npy::write(&dir.join("start.npy"), &gather, &traces_start).unwrap();
    let observed = record(&s, &truth(&case, model));
    npy::write(&dir.join("observed.npy"), &gather, &observed).unwrap();
    if case.planar {
        // The point scatterer's closed form is for a point source.
        return;
    }
    let t = Instant::now();
    let exact = exact_traces(&s);
    write_traces(&dir.join("jacobian_exact.npy"), &exact);
    let depths = case.depths();
    let from = depths
        .iter()
        .position(|&z| z >= case.exact_from)
        .expect("a column below `exact_from`");
    let misfits: Vec<f64> = traces
        .iter()
        .zip(&exact)
        .enumerate()
        .map(|(column, (got, exact))| {
            if column % n < from {
                f64::NAN
            } else {
                relative(got, exact)
            }
        })
        .collect();
    write_f64(&dir.join("born_misfit.f64"), &misfits);
    println!(
        "\nBorn traces against a point scatterer's closed form ({:.0} s), ‖difference‖ / ‖exact‖, from {:.0} m down:",
        t.elapsed().as_secs_f64(),
        depths[from]
    );
    println!("{:>8} {:>10} {:>10} {:>10}", "depth", "δvp", "δvs", "δρ");
    for c in (from..n).step_by(((n - from) / 12).max(1)) {
        println!(
            "{:>6.0} m {:>9.3}% {:>9.3}% {:>9.3}%",
            depths[c],
            100.0 * misfits[c],
            100.0 * misfits[n + c],
            100.0 * misfits[2 * n + c]
        );
    }
    for p in Parameter::ALL {
        let mut all = misfits[p as usize * n + from..(p as usize + 1) * n].to_vec();
        all.sort_by(f64::total_cmp);
        println!(
            "  δ{}: median {:.3}%, worst {:.3}%",
            p.name(),
            100.0 * all[all.len() / 2],
            100.0 * all[all.len() - 1]
        );
    }
    let exact_gn = gram(&exact);
    write_f64(&dir.join("gauss_newton_exact.f64"), &exact_gn);
    println!("JᵀJ against the closed form's, block by block (‖difference‖ / ‖exact‖):");
    println!("{:>8} {:>10} {:>10} {:>10}", "", "δvp", "δvs", "δρ");
    let (mut all_num, mut all_den) = (0.0, 0.0);
    for q in Parameter::ALL {
        let mut line = format!("{:>8}", format!("H^{}", q.name()));
        for p in 0..3 {
            let (mut num, mut den) = (0.0, 0.0);
            for i in from..n {
                for k in from..n {
                    let at = (q as usize * n + i) * side + p * n + k;
                    num += (m.gauss_newton[at] - exact_gn[at]).powi(2);
                    den += exact_gn[at].powi(2);
                }
            }
            line += &format!(" {:>9.3}%", 100.0 * (num / den).sqrt());
            all_num += num;
            all_den += den;
        }
        println!("{line}");
    }
    println!(
        "  all nine blocks: {:.3}%",
        100.0 * (all_num / all_den).sqrt()
    );
}

/// The trials' survey: the poster's in 2D with working layers and the long record, its
/// column every 20 m from 100 to 900 m and its receivers every 10 m whatever the grid.
const TRIAL: Case = Case {
    name: "trials/h10",
    column_top: 100.0,
    column_bottom: 900.0,
    column_step: 20.0,
    receiver_step: 10.0,
    ..POSTER_2D_LONG
};

/// The same over the floor, and between the walls of the code behind the poster.
const TRIAL_FLOOR: Case = Case {
    name: "trials/floor-h10",
    floor: POSTER_2D_FLOOR.floor,
    ..TRIAL
};
const TRIAL_BOX: Case = Case {
    name: "trials/box",
    pml: 0,
    around: 2,
    duration: 2.0,
    aside: 10.0,
    ..TRIAL
};

/// A small inclusion 500 m under the source in place of the disc, and the column every
/// 10 m for 700 m about it: the EGU talk's homogeneous case in the poster's box.
const TRIAL_INCLUSION: Case = Case {
    name: "trials/inclusion8-h10",
    anomaly: Anomaly::Blob {
        depth: 500.0,
        sigma: 12.0,
    },
    column_top: 150.0,
    column_bottom: 850.0,
    column_step: 10.0,
    frequency: 8.0,
    delay: 0.1875,
    ..TRIAL
};

/// `base` on a grid of spacing `h`, its margins and layers as wide in metres: at least
/// twelve nodes of layer, and four between the box and the layers.
fn spaced(base: Case, name: &'static str, h: f64) -> Case {
    let nodes = |n: usize| (n as f64 * base.spacing / h).round() as usize;
    let walls = base.pml == 0;
    Case {
        name,
        spacing: h,
        pml: if walls { 0 } else { nodes(base.pml).max(12) },
        around: if walls {
            nodes(base.around)
        } else {
            nodes(base.around).max(4)
        },
        ..base
    }
}

/// `base` with a Ricker wavelet of `f0`, delayed by a period and a half.
fn at_frequency(base: Case, name: &'static str, f0: f64) -> Case {
    Case {
        name,
        frequency: f0,
        delay: 1.5 / f0,
        ..base
    }
}

/// The time stepping and stencil of the code behind the poster: sixteen points, Taylor's
/// coefficients, steps of 1 ms, the time error kept.
fn as_2018(base: Case, name: &'static str) -> Case {
    Case {
        name,
        order: 16,
        plain_steps: Some(1e-3),
        ..base
    }
}

fn run_each(cases: &[Case]) {
    for &case in cases {
        run(case, SIDE_BY_SIDE);
        println!();
    }
}

#[test]
#[ignore = "minutes"]
fn trials_spacing() {
    run_each(&[
        spaced(TRIAL, "trials/h20", 20.0),
        TRIAL,
        spaced(TRIAL, "trials/h5", 5.0),
    ]);
}

#[test]
#[ignore = "minutes"]
fn trials_scheme() {
    let order = |name, order| Case {
        name,
        order,
        ..TRIAL
    };
    run_each(&[
        order("trials/order4", 4),
        order("trials/order12", 12),
        order("trials/order16", 16),
        Case {
            name: "trials/courant04",
            courant: 0.4,
            ..TRIAL
        },
        as_2018(TRIAL, "trials/as2018"),
    ]);
}

#[test]
#[ignore = "minutes"]
fn trials_layers() {
    run_each(&[
        Case {
            name: "trials/layers24",
            pml: 24,
            ..TRIAL
        },
        Case {
            name: "trials/far",
            around: 16,
            ..TRIAL
        },
        Case {
            name: "trials/record1",
            duration: 1.0,
            ..TRIAL
        },
        Case {
            name: "trials/record3",
            duration: 3.0,
            ..TRIAL
        },
    ]);
}

#[test]
#[ignore = "minutes"]
fn trials_frequency() {
    let h5 = spaced(TRIAL, "trials/h5", 5.0);
    run_each(&[
        at_frequency(TRIAL, "trials/f5-h10", 5.0),
        at_frequency(h5, "trials/f5-h5", 5.0),
        at_frequency(TRIAL, "trials/f20-h10", 20.0),
        at_frequency(h5, "trials/f20-h5", 20.0),
        at_frequency(TRIAL, "trials/f32-h10", 32.0),
        at_frequency(h5, "trials/f32-h5", 32.0),
        as_2018(at_frequency(TRIAL, "", 32.0), "trials/f32-as2018"),
    ]);
}

#[test]
#[ignore = "minutes"]
fn trials_floor() {
    run_each(&[TRIAL_FLOOR, spaced(TRIAL_FLOOR, "trials/floor-h5", 5.0)]);
}

#[test]
#[ignore = "minutes"]
fn trials_box() {
    let walls = |name, around| Case {
        name,
        around,
        ..TRIAL_BOX
    };
    let record = |name, duration| Case {
        name,
        duration,
        ..TRIAL_BOX
    };
    run_each(&[
        TRIAL_BOX,
        walls("trials/box-walls1", 1),
        walls("trials/box-walls3", 3),
        walls("trials/box-walls4", 4),
        record("trials/box-record19", 1.9),
        record("trials/box-record21", 2.1),
        spaced(TRIAL_BOX, "trials/box-h5", 5.0),
        as_2018(TRIAL_BOX, "trials/box-as2018"),
    ]);
}

#[test]
#[ignore = "minutes"]
fn trials_inclusion() {
    let h5 = spaced(TRIAL_INCLUSION, "trials/inclusion8-h5", 5.0);
    run_each(&[
        TRIAL_INCLUSION,
        h5,
        as_2018(TRIAL_INCLUSION, "trials/inclusion8-as2018"),
        at_frequency(TRIAL_INCLUSION, "trials/inclusion32-h10", 32.0),
        at_frequency(h5, "trials/inclusion32-h5", 32.0),
        as_2018(
            at_frequency(TRIAL_INCLUSION, "", 32.0),
            "trials/inclusion32-as2018",
        ),
    ]);
}

/// The inclusion in the EGU talk's medium (vp 2000 m/s, vs 1000 m/s, ρ 1500 kg/m³), at its
/// two frequencies.
const TRIAL_TALK: Case = Case {
    name: "trials/talk8-h10",
    start: [2000.0, 1000.0, 1500.0],
    ..TRIAL_INCLUSION
};

#[test]
#[ignore = "minutes"]
fn trials_talk() {
    let h5 = spaced(TRIAL_TALK, "trials/talk8-h5", 5.0);
    run_each(&[
        TRIAL_TALK,
        h5,
        as_2018(TRIAL_TALK, "trials/talk8-as2018"),
        at_frequency(TRIAL_TALK, "trials/talk32-h10", 32.0),
        at_frequency(h5, "trials/talk32-h5", 32.0),
        as_2018(at_frequency(TRIAL_TALK, "", 32.0), "trials/talk32-as2018"),
    ]);
}

#[test]
#[ignore = "half an hour"]
fn trials_talk_finest() {
    run(
        at_frequency(spaced(TRIAL_TALK, "", 2.5), "trials/talk32-h2.5", 32.0),
        SIDE_BY_SIDE,
    );
}

#[test]
#[ignore = "half an hour"]
fn trials_finest() {
    let h = 2.5;
    run_each(&[
        spaced(TRIAL, "trials/h2.5", h),
        at_frequency(spaced(TRIAL, "", h), "trials/f32-h2.5", 32.0),
        at_frequency(
            spaced(TRIAL_INCLUSION, "", h),
            "trials/inclusion32-h2.5",
            32.0,
        ),
    ]);
}

/// What differencing whole runs, as the code behind the poster did, costs: the Born traces
/// from one-sided differences of the traces for a change of `step` in one cell, and the
/// Hessian's column from one-sided differences of the gradient, against the exact ones.
fn differenced(case: Case) {
    let s = setup(case, true);
    let n = s.column.len();
    let shot = &s.inv.shots[0];
    let base = record(&s, &s.inv.model);
    let gradient_of = |model: &Model| -> Vec<f64> {
        let mut sim = Simulation::with_model(&shot.config, model.clone()).unwrap();
        let run = gradient(&mut sim, |r| shot.misfit.evaluate(&r.traces), None).unwrap();
        run.gradient
            .velocities(model)
            .iter()
            .flat_map(|k| s.column.iter().map(|&m| k[m]))
            .collect()
    };
    let gradient_base = gradient_of(&s.inv.model);
    println!(
        "{}: one-sided differences of whole runs against the exact Born traces and columns",
        case.name
    );
    println!("{:>16} {:>8} {:>12} {:>12}", "", "step", "traces", "column");
    for (p, c) in sampled(n) {
        let node = s.column[c];
        let exact = product(&s, p, node);
        let column = &exact.kernels[..3 * n];
        let born: Vec<f64> = exact.born.iter().map(|&v| v as f64).collect();
        for step in [100.0, 50.0, 5.0, 0.5, 0.05] {
            let mut model = s.inv.model.clone();
            model.values_mut(p)[node] += step as f32;
            let traces: Vec<f64> = record(&s, &model)
                .iter()
                .zip(&base)
                .map(|(&a, &b)| (a as f64 - b as f64) / step)
                .collect();
            let g: Vec<f64> = gradient_of(&model)
                .iter()
                .zip(&gradient_base)
                .map(|(a, b)| (a - b) / step)
                .collect();
            println!(
                "{:>16} {step:>8} {:>12.3e} {:>12.3e}",
                format!("δ{} at {:.0} m", p.name(), case.depths()[c]),
                relative(&traces, &born),
                relative(&g, column)
            );
        }
    }
}

#[test]
#[ignore = "minutes"]
fn trials_differenced() {
    differenced(TRIAL);
    println!();
    differenced(TRIAL_BOX);
}

#[test]
#[ignore = "about 8 minutes: `cargo test --release --test hessian_poster -- --ignored --nocapture --exact coarse`"]
fn coarse() {
    run(COARSE, 1);
}

#[test]
#[ignore = "about an hour"]
fn poster() {
    run(POSTER, 1);
}

#[test]
#[ignore = "hours"]
fn poster_long() {
    run(POSTER_LONG, 1);
}

#[test]
#[ignore = "a minute"]
fn poster_2d() {
    run(POSTER_2D, SIDE_BY_SIDE);
}

#[test]
#[ignore = "minutes"]
fn poster_2d_long() {
    run(POSTER_2D_LONG, SIDE_BY_SIDE);
}

#[test]
#[ignore = "minutes"]
fn poster_2d_floor() {
    run(POSTER_2D_FLOOR, SIDE_BY_SIDE);
}

#[test]
#[ignore = "a minute"]
fn poster_2d_box() {
    run(POSTER_2D_BOX, SIDE_BY_SIDE);
}

#[test]
#[ignore = "a minute, after `poster_2d_box`, whose block it reads"]
fn poster_2d_box_products() {
    products(POSTER_2D_BOX);
}

#[test]
#[ignore = "a minute"]
fn poster_2d_box_towards() {
    towards(POSTER_2D_BOX);
}

/// Born runs side by side on a 2D grid, which is too small for one run to use every thread.
const SIDE_BY_SIDE: usize = 8;

/// The columns the checks below take: each parameter's at a quarter and at three quarters
/// of the column.
fn sampled(n: usize) -> Vec<(Parameter, usize)> {
    ORDER
        .iter()
        .flat_map(|&p| [n / 4, 3 * n / 4].map(|c| (p, c)))
        .collect()
}

/// The block `run` left against whole products ([`hessian_product`]: a Born and two adjoint
/// runs each) for six of its columns: the Hessian's column, its Born-against-adjoint term
/// (B's column) and the rest (JᵀJ's column and Bᵀ's). And the same Hessian kernels in the
/// plane, with the Gauss–Newton product's (a gradient run driven by the Born traces), for
/// pictures.
fn products(case: Case) {
    let s = setup(case, true);
    let dir = case.dir().join("block");
    let m = BlockMatrices::read(&dir);
    let n = s.column.len();
    let side = 3 * n;
    let model = &s.inv.model;
    let shot = &s.inv.shots[0];
    let j = case.row();
    println!("the block against whole products, ‖difference‖ / ‖the product's column‖:");
    for (p, c) in sampled(n) {
        let t = Instant::now();
        let col = product(&s, p, s.column[c]);
        let column = p as usize * n + c;
        let apart = |part: usize, ours: &dyn Fn(usize) -> f64| {
            let theirs = &col.kernels[part * side..(part + 1) * side];
            let got: Vec<f64> = (0..side).map(ours).collect();
            relative(&got, theirs)
        };
        println!(
            "  δ{} at {:>5.0} m: the Hessian {:.1e}, the Born runs against the adjoint run {:.1e}, the forward run against the adjoint run's change {:.1e} ({:.0} s)",
            p.name(),
            case.depths()[c],
            apart(0, &|row| m.hessian[row * side + column]),
            apart(1, &|row| m.scattered[row * side + column]),
            apart(2, &|row| m.gauss_newton[row * side + column]
                + m.scattered[column * side + row]),
            t.elapsed().as_secs_f64()
        );
        // In the plane: the product, its two wavefield terms and Gauss–Newton's.
        let by_ref = col.speeds.each_ref().map(Vec::as_slice);
        let [first, second] = &col.product.terms;
        let source = shot
            .misfit
            .gauss_newton_source(std::slice::from_ref(&col.born))
            .unwrap();
        let mut sim = Simulation::with_model(&shot.config, model.clone()).unwrap();
        let gauss_newton = gradient(&mut sim, |_| Ok((0.0, source)), None)
            .unwrap()
            .gradient;
        for (name, k) in [
            ("product", col.product.velocities(model, by_ref)),
            ("born-adjoint", first.velocities(model)),
            ("forward-changed", second.velocities(model)),
            ("gauss-newton", gauss_newton.velocities(model)),
        ] {
            let (shape, data) = plane(model, [&k[0][..], &k[1][..], &k[2][..]], j);
            let file = format!("kernel-{name}-{}-{c:03}.npy", p.name());
            npy::write(&dir.join(file), &shape, &data).unwrap();
        }
    }
}

#[test]
#[ignore = "minutes, after `coarse`, whose block it reads"]
fn coarse_products() {
    products(COARSE);
}

#[test]
#[ignore = "about half an hour, after `poster`, whose block it reads"]
fn poster_products() {
    products(POSTER);
}

#[test]
#[ignore = "a minute, after `poster_2d`, whose block it reads"]
fn poster_2d_products() {
    products(POSTER_2D);
}

#[test]
#[ignore = "two minutes, after `poster_2d_floor`, whose block it reads"]
fn poster_2d_floor_products() {
    products(POSTER_2D_FLOOR);
}

/// The block `run` left against differences of the gradient, for six of its columns.
fn differences(case: Case) {
    let s = setup(case, true);
    let m = BlockMatrices::read(&case.dir().join("block"));
    let n = s.column.len();
    let side = 3 * n;
    println!(
        "the block's columns against extrapolated differences of the gradient (largest difference / the column's largest entry):"
    );
    for (p, c) in sampled(n) {
        let t = Instant::now();
        let step = 0.02 * START[p as usize];
        let differences = difference_column(&s, p, s.column[c], step);
        let column: Vec<f64> = (0..side)
            .map(|row| m.hessian[row * side + p as usize * n + c])
            .collect();
        let largest = column.iter().fold(0.0f64, |a, x| a.max(x.abs()));
        let worst = column
            .iter()
            .zip(&differences)
            .fold(0.0f64, |a, (k, d)| a.max((k - d).abs()));
        println!(
            "  δ{} at {:>5.0} m: {:.2e} ({:.0} s)",
            p.name(),
            case.depths()[c],
            worst / largest,
            t.elapsed().as_secs_f64()
        );
    }
}

#[test]
#[ignore = "minutes, after `coarse`, whose block it reads"]
fn coarse_differences() {
    differences(COARSE);
}

#[test]
#[ignore = "a minute, after `poster_2d`, whose block it reads"]
fn poster_2d_differences() {
    differences(POSTER_2D);
}

#[test]
#[ignore = "two minutes, after `poster_2d_floor`, whose block it reads"]
fn poster_2d_floor_differences() {
    differences(POSTER_2D_FLOOR);
}

/// With the start's own data the residuals vanish to the last bit: the block's second-order
/// part is zero, and a whole product is JᵀJ's column, the same numbers by two routes (the
/// second adjoint run, and the Born traces' inner products).
#[test]
#[ignore = "about 10 minutes"]
fn coarse_without_residual() {
    let s = setup(COARSE, false);
    let b = block(&s, 1);
    let m = BlockMatrices::of(&s, &b);
    let n = s.column.len();
    let side = 3 * n;
    assert!(m.scattered.iter().chain(&m.local).all(|&x| x == 0.0));
    for (p, c) in sampled(n) {
        let col = product(&s, p, s.column[c]);
        let column: Vec<f64> = (0..side)
            .map(|row| m.gauss_newton[row * side + p as usize * n + c])
            .collect();
        println!(
            "δ{} at {:.0} m: the product and JᵀJ's column differ by {:.2e} of the column",
            p.name(),
            s.case.depths()[c],
            relative(&col.kernels[..side], &column)
        );
    }
}

/// The misfit along directions on the column, measured: for every `<name>.f64` in the case's
/// `directions` (vp, vs and ρ's changes at the column's nodes, a parameter's nodes together,
/// as the block's columns are), the misfit of the start changed by [`ALONG`]'s multiples of
/// it, written to `<name>.misfit`. Central second differences of those are the curvature
/// vᵀHv without the Hessian, and positive multiples a line search.
fn along(case: Case) {
    let s = setup(case, true);
    let shot = &s.inv.shots[0];
    let n = s.column.len();
    let dir = case.dir().join("directions");
    let mut files: Vec<PathBuf> = fs::read_dir(&dir)
        .unwrap_or_else(|e| panic!("{}: {e}", dir.display()))
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.extension().is_some_and(|e| e == "f64"))
        .collect();
    files.sort();
    for file in files {
        let direction = read_f64(&file);
        assert_eq!(direction.len(), 3 * n, "{}", file.display());
        let misfits: Vec<f64> = ALONG
            .iter()
            .map(|&t| {
                let mut model = s.inv.model.clone();
                for p in Parameter::ALL {
                    let values = model.values_mut(p);
                    for (c, &node) in s.column.iter().enumerate() {
                        values[node] += (t * direction[p as usize * n + c]) as f32;
                    }
                }
                let traces = vec![record(&s, &model)];
                shot.misfit.evaluate(&traces).unwrap().0
            })
            .collect();
        println!(
            "{}: {}",
            file.file_stem().unwrap().to_string_lossy(),
            misfits
                .iter()
                .map(|m| format!("{m:.9e}"))
                .collect::<Vec<_>>()
                .join(" ")
        );
        write_f64(&file.with_extension("misfit"), &misfits);
    }
}

/// The multiples of a direction [`along`] measures the misfit at.
const ALONG: [f64; 7] = [-1.0, -0.5, 0.0, 0.25, 0.5, 1.0, 2.0];

#[test]
#[ignore = "seconds per direction, after `coarse`"]
fn coarse_along() {
    along(COARSE);
}

#[test]
#[ignore = "seconds per direction, after `poster`"]
fn poster_along() {
    along(POSTER);
}

#[test]
#[ignore = "seconds, after `poster_2d`"]
fn poster_2d_along() {
    along(POSTER_2D);
}

/// The multiples of the whole anomaly [`towards`] measures the misfit at: 1 is the model
/// the observed data came from.
const TOWARDS: [f64; 17] = [
    -0.5, -0.375, -0.25, -0.125, 0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0, 1.125,
    1.25, 1.375, 1.5,
];

/// Along the whole anomaly, not one cell of it: the direction from the start to the model
/// the observed data came from, in all three parameters and in each alone. For each, the
/// slope gᵀv, the two curvatures vᵀHv (a whole product) and vᵀJᵀJv (its Born traces), and
/// the misfit measured at [`TOWARDS`]'s multiples, into `towards/<name>.f64`: which
/// quadratic model, Newton's or Gauss–Newton's, follows the misfit to its minimum.
fn towards(case: Case) {
    let s = setup(case, true);
    let dir = case.dir().join("towards");
    fs::create_dir_all(&dir).unwrap();
    let shot = &s.inv.shots[0];
    let model = &s.inv.model;
    let target = truth(&case, model);
    let mut sim = Simulation::with_model(&shot.config, model.clone()).unwrap();
    let g = gradient(&mut sim, |r| shot.misfit.evaluate(&r.traces), None)
        .unwrap()
        .gradient
        .velocities(model);
    let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
    println!(
        "{}: along the anomaly, the slope, the curvatures and the misfit",
        case.name
    );
    for (name, which) in [
        ("all", [true; 3]),
        ("vp", [true, false, false]),
        ("vs", [false, true, false]),
        ("rho", [false, false, true]),
    ] {
        let speeds: [Vec<f64>; 3] = Parameter::ALL.map(|p| {
            let (from, to) = (model.values(p), target.values(p));
            let on = if which[p as usize] { 1.0 } else { 0.0 };
            from.iter()
                .zip(to)
                .map(|(&a, &b)| on * (b - a) as f64)
                .collect()
        });
        let by_ref = speeds.each_ref().map(Vec::as_slice);
        let change = Change::from_velocities(model, by_ref);
        let mut sim = Simulation::with_model(&shot.config, model.clone()).unwrap();
        let mut gauss_newton = 0.0;
        let product = hessian_product(
            &mut sim,
            &change,
            |r| shot.misfit.evaluate(&r.traces),
            |d| {
                let source = shot.misfit.gauss_newton_source(d)?;
                gauss_newton = source[0]
                    .iter()
                    .zip(&d[0])
                    .map(|(&w, &d)| w as f64 * d as f64)
                    .sum();
                Ok(source)
            },
            None,
        )
        .unwrap();
        let hv = product.velocities(model, by_ref);
        let slope: f64 = (0..3).map(|p| dot(&g[p], &speeds[p])).sum();
        let hessian: f64 = (0..3).map(|p| dot(&hv[p], &speeds[p])).sum();
        let misfits: Vec<f64> = TOWARDS
            .iter()
            .map(|&t| {
                let mut changed = model.clone();
                for p in Parameter::ALL {
                    let values = changed.values_mut(p);
                    for (v, d) in values.iter_mut().zip(&speeds[p as usize]) {
                        *v += (t * d) as f32;
                    }
                }
                shot.misfit.value(&[record(&s, &changed)]).unwrap()
            })
            .collect();
        println!(
            "  {name:>3}: gᵀv {slope:+.4e}, vᵀHv {hessian:+.4e}, vᵀJᵀJv {gauss_newton:+.4e}; the misfit {}",
            misfits
                .iter()
                .map(|m| format!("{m:.4e}"))
                .collect::<Vec<_>>()
                .join(" ")
        );
        let mut out = vec![slope, hessian, gauss_newton];
        out.extend(TOWARDS);
        out.extend(misfits);
        write_f64(&dir.join(format!("{name}.f64")), &out);
    }
}

#[test]
#[ignore = "a minute"]
fn poster_2d_towards() {
    towards(POSTER_2D_LONG);
}

#[test]
#[ignore = "two minutes"]
fn poster_2d_floor_towards() {
    towards(POSTER_2D_FLOOR);
}

#[test]
#[ignore = "seconds per direction, after `poster_2d_long`"]
fn poster_2d_long_along() {
    along(POSTER_2D_LONG);
}

#[test]
#[ignore = "seconds per direction, after `poster_2d_floor`"]
fn poster_2d_floor_along() {
    along(POSTER_2D_FLOOR);
}

/// The margins' choice: the start's traces against the closed form for narrower and wider
/// surveys across the plane and thinner layers.
#[test]
#[ignore = "minutes"]
fn coarse_forward() {
    for (pml, across) in [(12, 15), (12, 6), (10, 6), (8, 6), (12, 3)] {
        println!("layers {pml} nodes wide, {across} nodes on each side of the plane:");
        let case = Case {
            pml,
            across,
            ..COARSE
        };
        print_forward(&forward_misfits(&setup(case, false)).0, case.spacing);
    }
}

#[test]
#[ignore = "minutes"]
fn poster_forward() {
    print_forward(&forward_misfits(&setup(POSTER, false)).0, POSTER.spacing);
}
