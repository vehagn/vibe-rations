//! The inversion against an analytic solution of the inverse problem: `cargo test --release
//! --test analytic_inversion -- --ignored --nocapture --test-threads 1` (about 30 minutes on
//! six cores; a test's name after `analytic_inversion` runs that one alone).
//!
//! Full waveform inversion is nonlinear and has no closed form, but its linearisation about a
//! homogeneous full space has: Born's approximation with Stokes's solution as the Green's
//! function (`common::stokes`) gives the sensitivity of every recorded sample to every unknown,
//! the Jacobian J, and the Gauss–Newton step solves min ‖J δm - δd‖² (Tarantola 1984; Pratt,
//! Shin & Hicks 1998). Here the unknowns are the relative changes of vs at the nodes of a coarse
//! lattice, trilinear hats 45 m apart (125 of them), around a weak Gaussian (σ = 45 m, a quarter
//! of a dominant S wavelength) of +3% or +10% in vs, which that basis represents exactly. Three
//! point forces 300 m from it and eleven three-component receivers around it record 0.75 s.
//!
//! J is built column by column: a node's change of μ and λ is a point moment source whose
//! moment is the forward strain there (times δc), and a moment source's field is the spatial
//! derivative of Stokes's solution for its time function, which is tabulated. Each column is
//! summed over the nodes and shear positions as the staggered grid samples the model (as
//! `Problem::grid_integral` does; `J_grid`) or over the nodes alone (`J_cont`, the continuous
//! Born integral by the nodes' quadrature). `the_sensitivities_are_the_exact_kernels` checks
//! Jᵀr against the adjoint kernels of `common::stokes`, derived the other way round.
//!
//! - `born_and_gauss_newton`: Born data against the solver's data differences (the truth's
//!   traces less the start's, at +ε and -ε to split the nonlinearity from the solver's error),
//!   and the Gauss–Newton solutions for Born and for the solver's data, by damping.
//! - `lbfgs_against_gauss_newton`: the solver's gradient against Jᵀδd; the first L-BFGS step
//!   against it; L-BFGS on the same 125 unknowns (the library's `Lbfgs` and `line_search`), and
//!   the `Inversion` as it runs (every node free, preconditioned by the illumination), each
//!   iteration against the truth and the Gauss–Newton solution.

mod common;

use std::time::Instant;

use rayon::prelude::*;

use common::stokes::{Medium, Problem, Ricker, Table, stokes, velocity_kernels};
use common::*;
use vibe_rations::config::{Config, Parameter};
use vibe_rations::inversion::{Inversion, Lbfgs, StepLimits, line_search};
use vibe_rations::model::Model;
use vibe_rations::simulation::Simulation;

const H: f64 = 15.0;
const PML: usize = 12;
/// From the target's centre to the absorbing layers.
const HALF: f64 = 390.0;
const F0: f64 = 10.0;
const T0: f64 = 0.12;
const FORCE: f64 = 1e9;
const DURATION: f64 = 0.75;
/// The sources and the receivers on the axes, from the centre.
const FAR: f64 = 300.0;
/// The receivers at a cube's corners, from the centre along each axis.
const CORNER: f64 = 195.0;
/// The coarse lattice: `NC` nodes `COARSE` apart along each axis, centred on the target.
const COARSE: f64 = 45.0;
const NC: usize = 5;
const K: usize = NC * NC * NC;
const SIGMA: f64 = 45.0;
/// The target box: nodes up to `REACH` from the centre along each axis, all those where a
/// coarse hat is non-zero (they end 135 m out).
const REACH: i64 = 8;
const BN: usize = 2 * REACH as usize + 1;
/// The forward strain is tabulated this often for the moment sources' time functions.
const TABLE_STEP: f64 = 2.5e-4;
/// Damping, as a fraction of JᵀJ's largest eigenvalue.
const DAMPING: [f64; 6] = [0.0, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2];

fn medium() -> Medium {
    Medium {
        vp: VP,
        vs: VS,
        rho: RHO,
    }
}

/// The target's centre: a node.
fn centre() -> [f64; 3] {
    [PML as f64 * H + HALF; 3]
}

fn shift(x: [f64; 3], o: [f64; 3]) -> [f64; 3] {
    [x[0] + o[0], x[1] + o[1], x[2] + o[2]]
}

/// Each source's position and force axis: on the -x, -y and -z axes, each force across its
/// line to the centre.
fn sources() -> Vec<([f64; 3], usize)> {
    let c = centre();
    vec![
        (shift(c, [-FAR, 0.0, 0.0]), 1),
        (shift(c, [0.0, -FAR, 0.0]), 2),
        (shift(c, [0.0, 0.0, -FAR]), 0),
    ]
}

/// Opposite each source on the axes, and at the corners of a cube.
fn receivers() -> Vec<[f64; 3]> {
    let c = centre();
    let mut r = vec![
        shift(c, [FAR, 0.0, 0.0]),
        shift(c, [0.0, FAR, 0.0]),
        shift(c, [0.0, 0.0, FAR]),
    ];
    for sz in [-1.0, 1.0] {
        for sy in [-1.0, 1.0] {
            for sx in [-1.0, 1.0] {
                r.push(shift(c, [sx * CORNER, sy * CORNER, sz * CORNER]));
            }
        }
    }
    r
}

/// The survey's configuration; `pinned` gives dt and the band of the optimised coefficients.
fn text(pinned: Option<(f64, f64)>, iterations: usize) -> String {
    let n = (2.0 * HALF / H).round() as usize + 1 + 2 * PML;
    let (dt, band) = match pinned {
        Some((dt, ppw)) => (format!("dt = {dt}"), format!("optimised_for = {ppw}")),
        None => (String::new(), String::new()),
    };
    let mut s = format!(
        r#"
        [grid]
        points = [{n}, {n}, {n}]
        spacing = {H}
        [time]
        duration = {DURATION}
        {dt}
        [boundary]
        pml_width = {PML}
        free_surface = false
        pml_vp_max = {VP}
        [scheme]
        {band}
        [[model.layer]]
        top = 0.0
        vp = {VP}
        vs = {VS}
        rho = {RHO}
        [inversion]
        observed = "unused"
        parameters = ["vs"]
        balance = false
        iterations = {iterations}
        "#
    );
    for (x, axis) in sources() {
        let mut f = [0.0; 3];
        f[axis] = FORCE;
        s += &format!(
            "[[source]]\nposition = {x:?}\nforce = {f:?}\n\
             wavelet = {{ kind = \"ricker\", frequency = {F0} }}\n"
        );
    }
    s += "[receivers]\ncomponents = [\"vx\", \"vy\", \"vz\"]\n";
    for x in receivers() {
        s += &format!("[[receivers.point]]\nposition = {x:?}\n");
    }
    s
}

/// The survey with its operator pinned at the homogeneous start's (dt, the layers' P speed,
/// the optimised coefficients' band), as an inversion pins it, so the truth's traces are
/// computed by the same operator.
struct Survey {
    cfg: Config,
    start: Model,
    times: Vec<f64>,
}

fn survey(iterations: usize) -> Survey {
    let first = setup(&text(None, iterations));
    let pinned = (
        first.solver.dt,
        first.info.points_per_wavelength.max(2.0 / 0.9),
    );
    let times = first.sample_times();
    let cfg = Config::from_toml(&text(Some(pinned), iterations)).expect("a configuration");
    let start = Model::homogeneous(Simulation::grid(&cfg).unwrap(), VP, VS, RHO);
    Survey { cfg, start, times }
}

impl Survey {
    /// Each shot's traces of `model`, as `Receivers::traces` holds them.
    fn record(&self, model: &Model) -> Vec<Vec<Vec<f32>>> {
        (0..self.cfg.sources.len())
            .map(|s| {
                let mut sim =
                    Simulation::with_model(&self.cfg.shot(s).unwrap(), model.clone()).unwrap();
                sim.run(None, false).expect("a stable run");
                sim.receivers.traces
            })
            .collect()
    }

    fn rows_per_shot(&self) -> usize {
        3 * receivers().len() * self.times.len()
    }

    /// The interior index of the box node at `b` (offsets in nodes from the centre).
    fn node(&self, b: [i64; 3]) -> usize {
        let c = (centre()[0] / H).round() as i64;
        let [i, j, k] = b.map(|o| (c + o) as usize);
        self.start.grid.interior_idx(i, j, k)
    }
}

/// Every shot's traces end to end, in f64.
fn flatten(traces: &[Vec<Vec<f32>>]) -> Vec<f64> {
    traces
        .iter()
        .flatten()
        .flatten()
        .map(|&v| v as f64)
        .collect()
}

/// The box's nodes, offsets in nodes from the centre, x fastest.
fn box_nodes() -> Vec<[i64; 3]> {
    let r = -REACH..=REACH;
    r.clone()
        .flat_map(|k| {
            let r = r.clone();
            r.clone()
                .flat_map(move |j| r.clone().map(move |i| [i, j, k]))
        })
        .collect()
}

fn box_index(b: [i64; 3]) -> Option<usize> {
    let a = b.map(|o| o + REACH);
    let n = BN as i64;
    a.iter()
        .all(|&x| (0..n).contains(&x))
        .then(|| ((a[2] * n + a[1]) * n + a[0]) as usize)
}

/// The coarse hats' values at the node `b` (offsets in nodes from the centre): (hat, value)
/// where non-zero.
fn hats(b: [i64; 3]) -> Vec<(usize, f64)> {
    let along = |o: i64, c: usize| {
        let d = o as f64 * H - (c as f64 - (NC / 2) as f64) * COARSE;
        (1.0 - d.abs() / COARSE).max(0.0)
    };
    let mut out = Vec::new();
    for k in 0..K {
        let c = [k % NC, (k / NC) % NC, k / (NC * NC)];
        let w = along(b[0], c[0]) * along(b[1], c[1]) * along(b[2], c[2]);
        if w > 0.0 {
            out.push((k, w));
        }
    }
    out
}

/// The truth's coefficients: the Gaussian at the coarse nodes, peak `eps`.
fn truth_coefficients(eps: f64) -> Vec<f64> {
    (0..K)
        .map(|k| {
            let c = [k % NC, (k / NC) % NC, k / (NC * NC)];
            let r2: f64 = c
                .iter()
                .map(|&c| ((c as f64 - (NC / 2) as f64) * COARSE).powi(2))
                .sum();
            eps * (-0.5 * r2 / (SIGMA * SIGMA)).exp()
        })
        .collect()
}

/// The relative change of vs at the box's nodes that coefficients `c` describe.
fn nodal(c: &[f64]) -> Vec<f64> {
    box_nodes()
        .into_iter()
        .map(|b| hats(b).iter().map(|&(k, w)| w * c[k]).sum())
        .collect()
}

/// The start with vs changed by `c`'s hats.
fn model_of(s: &Survey, c: &[f64]) -> Model {
    let mut m = s.start.clone();
    for (b, d) in box_nodes().into_iter().zip(nodal(c)) {
        let n = s.node(b);
        m.vs[n] = (VS * (1.0 + d)) as f32;
    }
    m
}

/// Where the solver puts a model's terms: a node (λ, μ's normal part) or a shear position
/// (μ's shear part for axes i, j).
#[derive(Clone, Copy, Debug)]
enum Kind {
    Node,
    Edge(usize, usize),
}

/// A place a term of the model enters, with the basis functions' weights there: on the grid
/// (`grid`: a node's own values; an edge's, the mean of its four nodes', as the harmonic mean
/// linearises) and as point values at the nodes (`point`, for the continuous Born integral),
/// and the box nodes its share of a nodal gradient goes to.
struct Place {
    x: [f64; 3],
    kind: Kind,
    grid: Vec<(usize, f64)>,
    point: Vec<(usize, f64)>,
    nodes: Vec<(usize, f64)>,
}

/// The places around `centre` whose basis weights (`basis` at a node, by its offset in nodes)
/// may be non-zero, those within `reach` nodes, together with every shear position next to
/// one of those nodes.
fn places(
    centre: [f64; 3],
    reach: i64,
    basis: &dyn Fn([i64; 3]) -> Vec<(usize, f64)>,
) -> Vec<Place> {
    let at = |b: [i64; 3]| shift(centre, b.map(|o| o as f64 * H));
    let inside = |b: [i64; 3]| b.iter().all(|o| o.abs() <= reach);
    let mut out = Vec::new();
    for k in -reach..=reach {
        for j in -reach..=reach {
            for i in -reach..=reach {
                let b = [i, j, k];
                let w = basis(b);
                out.push(Place {
                    x: at(b),
                    kind: Kind::Node,
                    grid: w.clone(),
                    point: w,
                    nodes: box_index(b).map(|q| (q, 1.0)).into_iter().collect(),
                });
            }
        }
    }
    for (p, q) in [(0, 1), (0, 2), (1, 2)] {
        let r = 3 - p - q;
        for a in -reach - 1..=reach {
            for c in -reach - 1..=reach {
                for e in -reach..=reach {
                    let mut b = [0; 3];
                    b[p] = a;
                    b[q] = c;
                    b[r] = e;
                    let mut grid: Vec<(usize, f64)> = Vec::new();
                    let mut nodes = Vec::new();
                    for (dp, dq) in [(0, 0), (1, 0), (0, 1), (1, 1)] {
                        let mut n = b;
                        n[p] += dp;
                        n[q] += dq;
                        if !inside(n) {
                            continue;
                        }
                        for (k, w) in basis(n) {
                            match grid.iter_mut().find(|(j, _)| *j == k) {
                                Some((_, v)) => *v += 0.25 * w,
                                None => grid.push((k, 0.25 * w)),
                            }
                        }
                        if let Some(q) = box_index(n) {
                            nodes.push((q, 0.25));
                        }
                    }
                    let mut x = at(b);
                    x[p] += 0.5 * H;
                    x[q] += 0.5 * H;
                    out.push(Place {
                        x,
                        kind: Kind::Edge(p, q),
                        grid,
                        point: Vec::new(),
                        nodes,
                    });
                }
            }
        }
    }
    out
}

/// The six components of the symmetric moment, in this order.
const PAIRS: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];

fn pair(i: usize, j: usize) -> usize {
    PAIRS
        .iter()
        .position(|&p| p == (i.min(j), i.max(j)))
        .unwrap()
}

/// The moment per unit δμ of a change of vs alone (δλ = -2δμ) at `x`, from the source's
/// displacement u: Sᵢⱼ = ∂ᵢuⱼ + ∂ⱼuᵢ - 2δᵢⱼ ∇·u, tabulated from 0 to the record's end with its
/// first three time derivatives and its integral, each as a [`Table`] (whose pulse gives the
/// velocity a force with that time function causes).
fn moments(x: [f64; 3], source: [f64; 3], force: usize) -> Vec<Table> {
    let m = medium();
    let n = (DURATION / TABLE_STEP).ceil() as usize;
    let step = DURATION / n as f64;
    let d = sub(x, source);
    let mut v = vec![vec![[0.0; 5]; n + 1]; 6];
    for t in 0..=n {
        let tau = t as f64 * step;
        for level in 0..4 {
            let ricker = Ricker {
                f0: F0,
                t0: T0,
                amplitude: FORCE,
                level,
            };
            let (_, g) = stokes(&m, force, d, tau, &ricker);
            let div = g[0][0] + g[1][1] + g[2][2];
            for (vc, &(i, j)) in v.iter_mut().zip(&PAIRS) {
                let iso = if i == j { 2.0 * div } else { 0.0 };
                vc[t][level as usize] = g[i][j] + g[j][i] - iso;
            }
        }
    }
    v.into_iter()
        .map(|mut v| {
            // ∫₀ˢ, as `Problem::adjoint_pulse` integrates.
            for i in 1..=n {
                let (a, b) = (v[i - 1], v[i]);
                v[i][4] = a[4] + step / 2.0 * (a[0] + b[0]) + step * step / 12.0 * (a[1] - b[1]);
            }
            Table { step, v }
        })
        .collect()
}

/// What a unit relative change of vs at `place` does to the records of a source at `source`
/// (force along `force`): Born's scattered velocity, in the layout of the traces
/// ((component · receivers + receiver) · samples + sample), as the grid takes the place's
/// terms and, for a node, as a point scatterer takes all of them. The body force of the
/// scatterer is ∂ⱼ(δcᵢⱼₖₗ ∂ₗuₖ), so the scattered displacement is h³ δμ Σᵢⱼ ∂ⱼGₙᵢ * Sᵢⱼ, ∂ⱼ
/// along the offset from scatterer to receiver.
fn response(
    place: &Place,
    source: [f64; 3],
    force: usize,
    receivers: &[[f64; 3]],
    times: &[f64],
) -> (Vec<f64>, Vec<f64>) {
    let m = medium();
    let tables = moments(place.x, source, force);
    let scale = 2.0 * RHO * VS * VS * H * H * H;
    let (nr, ns) = (receivers.len(), times.len());
    let mut grid = vec![0.0; 3 * nr * ns];
    let mut point = if place.point.is_empty() {
        Vec::new()
    } else {
        vec![0.0; 3 * nr * ns]
    };
    for (r, &xr) in receivers.iter().enumerate() {
        let d = sub(xr, place.x);
        for (s, &t) in times.iter().enumerate() {
            let mut g = [0.0; 3];
            let mut p = [0.0; 3];
            match place.kind {
                Kind::Node => {
                    for i in 0..3 {
                        let (_, gr) = stokes(&m, i, d, t, &tables[pair(i, i)]);
                        for n in 0..3 {
                            g[n] += gr[n][i];
                        }
                    }
                    if !point.is_empty() {
                        p = g;
                        for i in 0..3 {
                            for j in 0..3 {
                                if i != j {
                                    let (_, gr) = stokes(&m, i, d, t, &tables[pair(i, j)]);
                                    for n in 0..3 {
                                        p[n] += gr[n][j];
                                    }
                                }
                            }
                        }
                    }
                }
                Kind::Edge(i, j) => {
                    let table = &tables[pair(i, j)];
                    let (_, gi) = stokes(&m, i, d, t, table);
                    let (_, gj) = stokes(&m, j, d, t, table);
                    for n in 0..3 {
                        g[n] += gi[n][j] + gj[n][i];
                    }
                }
            }
            for n in 0..3 {
                let row = (n * nr + r) * ns + s;
                grid[row] = scale * g[n];
                if !point.is_empty() {
                    point[row] = scale * p[n];
                }
            }
        }
    }
    (grid, point)
}

/// The Jacobian for the coarse hats, column-major (`j[k * rows + row]`), as the grid samples
/// the model and as point values at the nodes, and Jᵀ`data` at each box node, as the grid
/// samples it.
struct Jacobian {
    rows: usize,
    grid: Vec<f64>,
    point: Vec<f64>,
    nodal: Vec<f64>,
}

fn jacobian(s: &Survey, data: Option<&[f64]>) -> Jacobian {
    let rows_shot = s.rows_per_shot();
    let rows = rows_shot * sources().len();
    let mut out = Jacobian {
        rows,
        grid: vec![0.0; K * rows],
        point: vec![0.0; K * rows],
        nodal: vec![0.0; BN * BN * BN],
    };
    let all = places(centre(), REACH, &hats);
    let threads = rayon::current_num_threads();
    let rec = receivers();
    for (shot, (xs, force)) in sources().into_iter().enumerate() {
        let off = shot * rows_shot;
        let d = data.map(|d| &d[off..off + rows_shot]);
        let parts: Vec<(Vec<f64>, Vec<f64>, Vec<f64>)> = (0..threads)
            .into_par_iter()
            .map(|t| {
                let mut grid = vec![0.0; K * rows_shot];
                let mut point = vec![0.0; K * rows_shot];
                let mut nodal = vec![0.0; BN * BN * BN];
                for place in all.iter().skip(t).step_by(threads) {
                    if place.grid.is_empty() && place.nodes.is_empty() {
                        continue;
                    }
                    let (g, p) = response(place, xs, force, &rec, &s.times);
                    for &(k, w) in &place.grid {
                        let col = &mut grid[k * rows_shot..(k + 1) * rows_shot];
                        col.iter_mut().zip(&g).for_each(|(a, b)| *a += w * b);
                    }
                    for &(k, w) in &place.point {
                        let col = &mut point[k * rows_shot..(k + 1) * rows_shot];
                        col.iter_mut().zip(&p).for_each(|(a, b)| *a += w * b);
                    }
                    if let Some(d) = d {
                        let dot: f64 = g.iter().zip(d).map(|(a, b)| a * b).sum();
                        for &(q, w) in &place.nodes {
                            nodal[q] += w * dot;
                        }
                    }
                }
                (grid, point, nodal)
            })
            .collect();
        for (grid, point, nodal) in parts {
            for k in 0..K {
                let to = k * rows + off;
                let from = k * rows_shot;
                for r in 0..rows_shot {
                    out.grid[to + r] += grid[from + r];
                    out.point[to + r] += point[from + r];
                }
            }
            out.nodal.iter_mut().zip(nodal).for_each(|(a, b)| *a += b);
        }
    }
    out
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn norm2(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

fn cosine(a: &[f64], b: &[f64]) -> f64 {
    dot(a, b) / (norm2(a) * norm2(b))
}

/// ‖a - b‖ / ‖b‖.
fn relative(a: &[f64], b: &[f64]) -> f64 {
    let d: Vec<f64> = a.iter().zip(b).map(|(x, y)| x - y).collect();
    norm2(&d) / norm2(b)
}

fn apply(j: &[f64], rows: usize, c: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0; rows];
    for (k, &ck) in c.iter().enumerate() {
        out.iter_mut()
            .zip(&j[k * rows..(k + 1) * rows])
            .for_each(|(o, v)| *o += ck * v);
    }
    out
}

fn transpose_apply(j: &[f64], rows: usize, d: &[f64]) -> Vec<f64> {
    (0..K)
        .into_par_iter()
        .map(|k| dot(&j[k * rows..(k + 1) * rows], d))
        .collect()
}

/// JᵀJ's eigenvalues, largest first, and eigenvectors (`v[i]`), by cyclic Jacobi rotations in
/// f64: small enough not to need a library.
struct Eigen {
    values: Vec<f64>,
    vectors: Vec<Vec<f64>>,
}

fn normal_eigen(j: &[f64], rows: usize) -> Eigen {
    let n = K;
    let mut a: Vec<f64> = (0..n * n)
        .into_par_iter()
        .map(|e| {
            let (p, q) = (e / n, e % n);
            dot(&j[p * rows..(p + 1) * rows], &j[q * rows..(q + 1) * rows])
        })
        .collect();
    let original = a.clone();
    let mut v = vec![0.0; n * n];
    (0..n).for_each(|i| v[i * n + i] = 1.0);
    let scale: f64 = a.iter().map(|x| x * x).sum::<f64>();
    for _ in 0..100 {
        let off: f64 = (0..n)
            .flat_map(|p| (0..n).filter(move |&q| q != p).map(move |q| (p, q)))
            .map(|(p, q)| a[p * n + q].powi(2))
            .sum();
        if off < 1e-30 * scale {
            break;
        }
        for p in 0..n {
            for q in p + 1..n {
                let apq = a[p * n + q];
                if apq == 0.0 {
                    continue;
                }
                let theta = (a[q * n + q] - a[p * n + p]) / (2.0 * apq);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let (kp, kq) = (a[k * n + p], a[k * n + q]);
                    a[k * n + p] = c * kp - s * kq;
                    a[k * n + q] = s * kp + c * kq;
                }
                for k in 0..n {
                    let (pk, qk) = (a[p * n + k], a[q * n + k]);
                    a[p * n + k] = c * pk - s * qk;
                    a[q * n + k] = s * pk + c * qk;
                }
                for k in 0..n {
                    let (kp, kq) = (v[k * n + p], v[k * n + q]);
                    v[k * n + p] = c * kp - s * kq;
                    v[k * n + q] = s * kp + c * kq;
                }
            }
        }
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&x, &y| a[y * n + y].total_cmp(&a[x * n + x]));
    let values: Vec<f64> = order.iter().map(|&i| a[i * n + i]).collect();
    let vectors: Vec<Vec<f64>> = order
        .iter()
        .map(|&i| (0..n).map(|k| v[k * n + i]).collect())
        .collect();
    // A v = λ v, to the rotations' round-off.
    let worst = values
        .iter()
        .zip(&vectors)
        .map(|(&l, x)| {
            (0..n)
                .map(|p| (dot(&original[p * n..(p + 1) * n], x) - l * x[p]).abs())
                .fold(0.0, f64::max)
        })
        .fold(0.0, f64::max);
    assert!(worst < 1e-10 * values[0], "eigenvectors off by {worst:e}");
    Eigen { values, vectors }
}

impl Eigen {
    /// The damped least-squares solution (JᵀJ + μ λ_max I)⁻¹ Jᵀd, from b = Jᵀd.
    fn solve(&self, b: &[f64], mu: f64) -> Vec<f64> {
        let alpha = mu * self.values[0];
        let mut c = vec![0.0; K];
        for (&l, v) in self.values.iter().zip(&self.vectors) {
            let w = dot(v, b) / (l + alpha);
            if w.is_finite() {
                c.iter_mut().zip(v).for_each(|(c, v)| *c += w * v);
            }
        }
        c
    }
}

/// The error of a relative change against the truth's, over the same nodes, relative to the
/// start's (the truth's own size), and the cosine between them.
fn score(got: &[f64], truth: &[f64]) -> (f64, f64) {
    (relative(got, truth), cosine(got, truth))
}

/// ½‖d - d_obs‖² of the solver's traces of `model`.
fn misfit_of(s: &Survey, model: &Model, observed: &[f64]) -> f64 {
    let d = flatten(&s.record(model));
    0.5 * d
        .iter()
        .zip(observed)
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
}

fn pc(x: f64) -> String {
    format!("{:.2}%", 100.0 * x)
}

/// Born data against the solver's at two strengths, and Gauss–Newton solutions by damping.
#[test]
#[ignore = "a study, not a check: run with `cargo test --release --test analytic_inversion -- --ignored`"]
fn born_and_gauss_newton() {
    let t = Instant::now();
    let s = survey(0);
    let base = flatten(&s.record(&s.start));
    println!(
        "{} samples of {} s per trace, {} rows, {K} unknowns; start recorded in {:.0} s",
        s.times.len(),
        s.times[1] - s.times[0],
        base.len(),
        t.elapsed().as_secs_f64()
    );
    let t = Instant::now();
    let j = jacobian(&s, None);
    println!("J in {:.0} s", t.elapsed().as_secs_f64());
    let t = Instant::now();
    let eigen = normal_eigen(&j.grid, j.rows);
    let eigen_point = normal_eigen(&j.point, j.rows);
    let l = &eigen.values;
    println!(
        "JᵀJ (grid): eigenvalues {:.3e} to {:.3e}, condition {:.2e}; quantiles of λ/λ_max: \
         10% {:.1e}, 50% {:.1e}, 90% {:.1e} ({:.0} s)",
        l[0],
        l[K - 1],
        l[0] / l[K - 1],
        l[K / 10] / l[0],
        l[K / 2] / l[0],
        l[9 * K / 10] / l[0],
        t.elapsed().as_secs_f64()
    );
    for eps in [0.03, 0.1] {
        let truth = truth_coefficients(eps);
        let truth_nodal = nodal(&truth);
        let diff = |c: &[f64]| -> Vec<f64> {
            flatten(&s.record(&model_of(&s, c)))
                .iter()
                .zip(&base)
                .map(|(a, b)| a - b)
                .collect()
        };
        let plus = diff(&truth);
        let minus = diff(&truth.iter().map(|c| -c).collect::<Vec<_>>());
        let born = apply(&j.grid, j.rows, &truth);
        let born_point = apply(&j.point, j.rows, &truth);
        let linear: Vec<f64> = plus
            .iter()
            .zip(&minus)
            .map(|(a, b)| 0.5 * (a - b))
            .collect();
        let quadratic: Vec<f64> = plus
            .iter()
            .zip(&minus)
            .map(|(a, b)| 0.5 * (a + b))
            .collect();
        let peak = base.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        let dpeak = plus.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        println!("\n+{:.0}% in vs at the peak", 100.0 * eps);
        println!(
            "  solver's data difference: peak {} of the traces' peak, norm {} of theirs",
            pc(dpeak / peak),
            pc(norm2(&plus) / norm2(&base))
        );
        println!(
            "  Born (grid) against the solver's difference: {}, cosine {:.5}",
            pc(relative(&born, &plus)),
            cosine(&born, &plus)
        );
        println!(
            "  the solver's part odd in ε (linear, + O(ε³)) against Born: {}; even (the \
             nonlinearity's second order): {} of the difference",
            pc(relative(&linear, &born)),
            pc(norm2(&quadratic) / norm2(&plus))
        );
        println!(
            "  Born as point values at the nodes against Born on the grid: {}",
            pc(relative(&born_point, &born))
        );
        let observed: Vec<f64> = base.iter().zip(&plus).map(|(a, b)| a + b).collect();
        let chi0 = 0.5 * norm2(&plus).powi(2);
        println!(
            "  {:<34} {:>8} | {:>7} {:>8} | {:>9} {:>9}",
            "Gauss–Newton: data, J, damping", "", "error", "cosine", "lin. χ/χ₀", "χ/χ₀"
        );
        // The odd part has no second-order nonlinearity, so it shows what the solver's own
        // error does to the solution.
        let cases: [(&str, &[f64], &Eigen, &[f64]); 4] = [
            ("Born", &born, &eigen, &j.grid),
            ("solver", &plus, &eigen, &j.grid),
            ("solver, odd in ε", &linear, &eigen, &j.grid),
            ("solver, J point", &plus, &eigen_point, &j.point),
        ];
        for (name, data, eig, jm) in cases {
            let b = transpose_apply(jm, j.rows, data);
            for mu in DAMPING {
                let c = eig.solve(&b, mu);
                let (e, cos) = score(&nodal(&c), &truth_nodal);
                let pred = apply(jm, j.rows, &c);
                let lin = relative(&pred, data);
                // The solver's misfit of the solution, against the solver's data.
                let misfit = if name == "Born" || name.contains("odd") {
                    String::from("-")
                } else {
                    pc(misfit_of(&s, &model_of(&s, &c), &observed) / chi0)
                };
                println!(
                    "  {:<34} {:>8.0e} | {:>7} {:>+8.5} | {:>9} {:>9}",
                    name,
                    mu,
                    pc(e),
                    cos,
                    pc(lin * lin),
                    misfit
                );
            }
        }
    }
}

/// The residual a 2% faster medium leaves at one receiver, against Jᵀ of it: Jᵀr Δt must be
/// the exact kernels integrated against the hats, as the grid samples them
/// (`Problem::grid_integral`) and as point values (`Problem::kernels`). Both sides are the
/// same correlations of Stokes's fields, one scattered forward from the hat to the receiver,
/// the other sent back from the receiver.
#[test]
fn the_sensitivities_are_the_exact_kernels() {
    let m = medium();
    let source = [0.0, 0.0, 0.0];
    let receiver = [420.0, 0.0, 0.0];
    let p = Problem {
        medium: m,
        observed: Medium {
            vp: 1.02 * VP,
            vs: 1.02 * VS,
            rho: RHO,
        },
        source,
        force: 1,
        amplitude: FORCE,
        f0: F0,
        t0: T0,
        receiver,
        component: 1,
        duration: DURATION,
    };
    let dt = 1e-3;
    let times: Vec<f64> = (0..750).map(|i| (i as f64 + 0.5) * dt).collect();
    let residual: Vec<f64> = times.iter().map(|&t| p.residual(0, t)).collect();
    // One hat 30 m wide each way, at the midpoint and off the line.
    let reach = 2;
    let hat = |b: [i64; 3]| {
        let w: f64 = b
            .iter()
            .map(|&o| (1.0 - (o as f64).abs() / reach as f64).max(0.0))
            .product();
        if w > 0.0 { vec![(0, w)] } else { Vec::new() }
    };
    for c in [[210.0, 0.0, 0.0], [195.0, 60.0, -45.0]] {
        let all = places(c, reach, &hat);
        let (grid, point) = all
            .par_iter()
            .map(|place| {
                let wg: f64 = place.grid.iter().map(|x| x.1).sum();
                let wp: f64 = place.point.iter().map(|x| x.1).sum();
                if wg == 0.0 && wp == 0.0 {
                    return (0.0, 0.0);
                }
                let (g, pt) = response(place, source, 1, &[receiver], &times);
                // The scattered vy.
                let pick = |v: &[f64]| v[times.len()..2 * times.len()].to_vec();
                let g = dot(&pick(&g), &residual) * wg;
                let pt = if wp > 0.0 {
                    dot(&pick(&pt), &residual) * wp
                } else {
                    0.0
                };
                (g, pt)
            })
            .reduce(|| (0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1));
        let (grid, point) = (grid * dt, point * dt);
        let nodes: Vec<[f64; 3]> = all
            .iter()
            .filter(|x| matches!(x.kind, Kind::Node))
            .map(|x| x.x)
            .collect();
        let weight = |x: [f64; 3]| {
            let b = sub(x, c).map(|v| (v / H).round() as i64);
            hat(b).first().map_or(0.0, |w| w.1)
        };
        let exact_grid = VS * velocity_kernels(&m, p.grid_integral(&nodes, H, weight, 1e-3))[1];
        let kernels = p.kernels(&nodes, 1e-3);
        let exact_point: f64 = nodes
            .iter()
            .zip(kernels)
            .map(|(&x, k)| weight(x) * velocity_kernels(&m, k)[1])
            .sum::<f64>()
            * VS
            * H.powi(3);
        let (eg, ep) = (grid / exact_grid - 1.0, point / exact_point - 1.0);
        println!(
            "hat at {c:?}: Jᵀr Δt {grid:.6e} against {exact_grid:.6e} on the grid ({eg:+.1e}), \
             {point:.6e} against {exact_point:.6e} as point values ({ep:+.1e})"
        );
        assert!(eg.abs() < 1e-8 && ep.abs() < 1e-8, "{eg:e} {ep:e}");
    }
}

/// L-BFGS against Gauss–Newton at +3%: the gradient, the first step, L-BFGS on the coarse
/// unknowns and the inversion as it runs.
#[test]
#[ignore = "a study, not a check: run with `cargo test --release --test analytic_inversion -- --ignored`"]
fn lbfgs_against_gauss_newton() {
    const EPS: f64 = 0.03;
    const ITERATIONS: usize = 8;
    let s = survey(ITERATIONS);
    let truth = truth_coefficients(EPS);
    let truth_nodal = nodal(&truth);
    let start_traces = s.record(&s.start);
    let observed_traces = s.record(&model_of(&s, &truth));
    let base = flatten(&start_traces);
    let observed = flatten(&observed_traces);
    let delta: Vec<f64> = observed.iter().zip(&base).map(|(a, b)| a - b).collect();
    let chi0 = 0.5 * norm2(&delta).powi(2);
    let t = Instant::now();
    let j = jacobian(&s, Some(&delta));
    println!("J in {:.0} s", t.elapsed().as_secs_f64());
    let eigen = normal_eigen(&j.grid, j.rows);
    let b = transpose_apply(&j.grid, j.rows, &delta);
    println!("Gauss–Newton on the solver's data (J on the grid), by damping:");
    let mut references = Vec::new();
    for mu in DAMPING {
        let c = eigen.solve(&b, mu);
        let (e, cos) = score(&nodal(&c), &truth_nodal);
        println!("  {mu:>6.0e}: error {}, cosine {cos:+.5}", pc(e));
        references.push(c);
    }

    let settings = s.cfg.inversion.clone().unwrap();
    let t = Instant::now();
    let mut inv = Inversion::new(
        &s.cfg,
        s.start.clone(),
        observed_traces.into_iter().map(Into::into).collect(),
        &settings,
    )
    .unwrap();
    println!(
        "\nThe inversion set up in {:.0} s; its misfit {:.6e}, ½‖δd‖² {chi0:.6e}",
        t.elapsed().as_secs_f64(),
        inv.misfit()
    );
    // Free nodes: every node outside the layers, z slowest, as `Inversion::new` lists them.
    let ni = s.start.grid.n[0] - 2 * PML;
    assert_eq!(inv.unknowns(), ni * ni * ni);
    let free = |n: usize| {
        let g = s.start.grid;
        let (i, j, k) = (n % g.n[0], (n / g.n[0]) % g.n[1], n / (g.n[0] * g.n[1]));
        ((k - PML) * ni + (j - PML)) * ni + (i - PML)
    };
    let boxed: Vec<usize> = box_nodes().into_iter().map(|b| s.node(b)).collect();
    let g0 = inv.gradient_of(Parameter::Vs).unwrap();
    let g_box: Vec<f64> = boxed.iter().map(|&n| -g0[n] as f64).collect();
    println!(
        "The solver's gradient against Jᵀδd at the box's nodes: {}, cosine {:+.5}",
        pc(relative(&g_box, &j.nodal)),
        cosine(&g_box, &j.nodal)
    );
    let coarse = |g: &[f64]| -> Vec<f64> {
        let mut c = vec![0.0; K];
        for (b, &n) in box_nodes().into_iter().zip(&boxed) {
            for (k, w) in hats(b) {
                c[k] += w * g[free(n)];
            }
        }
        c
    };
    let e0 = inv.evaluate(&s.start).unwrap();
    let gc = coarse(&e0.gradient);
    let minus_gc: Vec<f64> = gc.iter().map(|g| -g).collect();
    println!(
        "  on the hats: {}, cosine {:+.5}",
        pc(relative(&minus_gc, &b)),
        cosine(&minus_gc, &b)
    );
    // The preconditioner the inversion divides by: the source side's illumination, its
    // largest over the free nodes, with the water level.
    let imax = (0..s.start.grid.interior_len())
        .filter(|&n| {
            let g = s.start.grid;
            let i = [n % g.n[0], (n / g.n[0]) % g.n[1], n / (g.n[0] * g.n[1])];
            i.iter().all(|&a| a >= PML && a < g.n[0] - PML)
        })
        .map(|n| e0.illumination[n])
        .fold(0.0, f64::max);
    let precondition = |n: usize| 1.0 / (e0.illumination[n] / imax + settings.water_level);
    let preconditioned: Vec<f64> = boxed
        .iter()
        .zip(&j.nodal)
        .map(|(&n, a)| precondition(n) * a)
        .collect();
    let p_box: Vec<f64> = boxed.iter().map(|&n| precondition(n)).collect();
    println!(
        "  the preconditioner over the box: {:.1} to {:.1} (1/water level = {:.0})",
        p_box.iter().fold(f64::INFINITY, |m, &v| m.min(v)),
        p_box.iter().fold(0.0f64, |m, &v| m.max(v)),
        1.0 / settings.water_level
    );
    let share = |v: &dyn Fn(usize) -> f64| -> f64 {
        let all: f64 = (0..s.start.grid.interior_len()).map(|n| v(n).powi(2)).sum();
        let inside: f64 = boxed.iter().map(|&n| v(n).powi(2)).sum();
        (inside / all).sqrt()
    };
    println!(
        "  the gradient's norm in the box: {} of its whole; preconditioned, {}",
        pc(share(&|n| g0[n] as f64)),
        pc(share(&|n| g0[n] as f64 * precondition(n)))
    );

    // L-BFGS on the hats' coefficients: the library's direction and line search, the solver's
    // gradient taken onto the hats, no preconditioner.
    println!(
        "\nL-BFGS on the {K} hats: iteration, misfit left, error, cosine, error against \
         Gauss–Newton (1e-6, 1e-4)"
    );
    let t = Instant::now();
    let mut lbfgs = Lbfgs::new(settings.memory);
    let limits = StepLimits {
        initial: settings.initial_step,
        max: settings.max_step,
    };
    let mut c = vec![0.0; K];
    let (mut f, mut g) = (e0.misfit, gc.clone());
    let against = |c: &[f64], r: &[f64]| relative(&nodal(c), &nodal(r));
    for it in 1..=ITERATIONS {
        let step = line_search(
            &mut lbfgs,
            &c,
            f,
            &g,
            None,
            limits,
            |_| {},
            |x| inv.misfit_of(&model_of(&s, x)),
        )
        .unwrap();
        let Some(step) = step else {
            println!("  no step lowers the misfit");
            break;
        };
        let e = inv.evaluate(&model_of(&s, &step.x)).unwrap();
        let g_new = coarse(&e.gradient);
        let sv: Vec<f64> = step.x.iter().zip(&c).map(|(a, b)| a - b).collect();
        let yv: Vec<f64> = g_new.iter().zip(&g).map(|(a, b)| a - b).collect();
        lbfgs.update(sv, yv);
        c = step.x;
        f = e.misfit;
        g = g_new;
        let (err, cos) = score(&nodal(&c), &truth_nodal);
        println!(
            "  {it}: {} | {} {cos:+.5} | {} {} ({} trials)",
            pc(f / chi0),
            pc(err),
            pc(against(&c, &references[2])),
            pc(against(&c, &references[3])),
            step.trials
        );
    }
    println!("  {:.0} s", t.elapsed().as_secs_f64());

    // The inversion as it runs: every node free, preconditioned.
    println!(
        "\nThe inversion (every node free, preconditioned): iteration, misfit left, error and \
         cosine in the box, error over the whole grid, error against Gauss–Newton (1e-4)"
    );
    let t = Instant::now();
    let gn = nodal(&references[3]);
    let whole = |m: &Model| -> f64 {
        let truth_model = model_of(&s, &truth);
        let d: f64 =
            m.vs.iter()
                .zip(&truth_model.vs)
                .map(|(a, b)| ((a - b) as f64 / VS).powi(2))
                .sum();
        d.sqrt() / norm2(&truth_nodal)
    };
    for it in 1..=ITERATIONS {
        let Some(step) = inv.iterate().unwrap() else {
            println!("  no step lowers the misfit");
            break;
        };
        let change: Vec<f64> = boxed
            .iter()
            .map(|&n| inv.model.vs[n] as f64 / VS - 1.0)
            .collect();
        let (err, cos) = score(&change, &truth_nodal);
        println!(
            "  {it}: {} | {} {cos:+.5} | {} | {} (step {:.3}, {} trials)",
            pc(inv.misfit() / chi0),
            pc(err),
            pc(whole(&inv.model)),
            pc(relative(&change, &gn)),
            step.step,
            step.trials
        );
        if it == 1 {
            let all: Vec<f64> = inv.model.vs.iter().map(|&v| v as f64 / VS - 1.0).collect();
            let largest = (0..all.len())
                .max_by(|&a, &b| all[a].abs().total_cmp(&all[b].abs()))
                .unwrap();
            let x = s.start.grid.position([
                largest % s.start.grid.n[0],
                (largest / s.start.grid.n[0]) % s.start.grid.n[1],
                largest / (s.start.grid.n[0] * s.start.grid.n[1]),
            ]);
            let nearest = sources()
                .iter()
                .map(|s| s.0)
                .chain(receivers())
                .map(|y| norm(sub(x, y)))
                .fold(f64::INFINITY, f64::min);
            println!(
                "     its first step against Jᵀδd in the box: cosine {:+.5}; against the \
                 preconditioned Jᵀδd {:+.5}; its norm in the box {} of its whole; its \
                 largest change {} at {:.0} m from a source or receiver",
                cosine(&change, &j.nodal),
                cosine(&change, &preconditioned),
                pc(share(&|n| all[n])),
                pc(all[largest]),
                nearest
            );
        }
    }
    println!("  {:.0} s", t.elapsed().as_secs_f64());
}

/// The solver's Gauss–Newton products against Born's in the full space: JᵀJ c with J built
/// from Stokes's solution as the grid samples the model, against the inversion's
/// `gauss_newton_product` of the hats' nodal values, taken back to the hats (Pᵀ H P c, P the
/// hats' interpolation), for the truth's coefficients, a single hat at the centre and two
/// random sets. The misfit is unweighted (`balance = false`), as J is; the forward wavefield
/// is replayed from checkpoints, so the products are exact for the discrete run and what
/// differs is the run's error against Stokes's. About 6 minutes.
///
/// Measured: against JᵀJc with J as the grid samples the model, 1.6e-5 to 2.8e-5 (the quadratic
/// forms to 1e-6 to 2.5e-5); against the continuous Born integral's (J at the nodes'
/// quadrature), 1.3–3.3%: the grid's sampling of the model, not the wavefields.
#[test]
#[ignore]
fn gauss_newton_products_against_born() {
    let s = survey(1);
    let t = Instant::now();
    let j = jacobian(&s, None);
    println!("J in {:.0} s", t.elapsed().as_secs_f64());
    let mut settings = s.cfg.inversion.clone().unwrap();
    settings.reconstruct = false;
    let observed = s.record(&s.start);
    let inv = Inversion::new(
        &s.cfg,
        s.start.clone(),
        observed.into_iter().map(Into::into).collect(),
        &settings,
    )
    .unwrap();
    // Free nodes: every node outside the layers, z slowest, as `Inversion::new` lists them.
    let ni = s.start.grid.n[0] - 2 * PML;
    assert_eq!(inv.unknowns(), ni * ni * ni);
    let free = |n: usize| {
        let g = s.start.grid;
        let (i, j, k) = (n % g.n[0], (n / g.n[0]) % g.n[1], n / (g.n[0] * g.n[1]));
        ((k - PML) * ni + (j - PML)) * ni + (i - PML)
    };
    let boxed: Vec<usize> = box_nodes().into_iter().map(|b| s.node(b)).collect();
    let mut seed = 0x5EED_u64;
    let mut random = || -> Vec<f64> {
        (0..K)
            .map(|_| {
                seed ^= seed << 13;
                seed ^= seed >> 7;
                seed ^= seed << 17;
                (seed >> 11) as f64 / (1u64 << 53) as f64 - 0.5
            })
            .collect()
    };
    let mut centre_hat = vec![0.0; K];
    centre_hat[K / 2] = 1.0;
    let cases = [
        ("the truth's Gaussian", truth_coefficients(1.0)),
        ("the centre's hat", centre_hat),
        ("random", random()),
        ("random", random()),
    ];
    for (name, c) in cases {
        let t = Instant::now();
        let mut v = vec![0.0; inv.unknowns()];
        for (&n, d) in boxed.iter().zip(nodal(&c)) {
            v[free(n)] = d;
        }
        let hv = inv.gauss_newton_product(&s.start, &v).unwrap();
        let mut solver = vec![0.0; K];
        for (b, &n) in box_nodes().into_iter().zip(&boxed) {
            for (k, w) in hats(b) {
                solver[k] += w * hv[free(n)];
            }
        }
        let born = transpose_apply(&j.grid, j.rows, &apply(&j.grid, j.rows, &c));
        let point = transpose_apply(&j.point, j.rows, &apply(&j.point, j.rows, &c));
        println!(
            "{name}: PᵀHPc against JᵀJc {:.1e} (cᵀPᵀHPc {:.7e} against {:.7e}); against the \
             continuous Born integral's {:.1e} ({:.0} s)",
            relative(&solver, &born),
            dot(&c, &solver),
            dot(&c, &born),
            relative(&solver, &point),
            t.elapsed().as_secs_f64()
        );
    }
}
