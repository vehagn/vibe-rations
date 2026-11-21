//! How accurate the solver's gradients are against the gradient of the wave equation itself:
//! `just gradient-accuracy`.
//!
//! The finite-difference checks (`tests/physics.rs`) show the gradient is that of the
//! discrete run; this measures how far that is from the exact one. In a homogeneous full
//! space a point force's wavefield is Stokes's solution, so the kernels of the L2 misfit of a
//! velocity trace are closed-form correlations (`common::stokes`). The observed trace comes
//! from a medium 2% faster, so the residual is the kind an inversion sees: mostly a shift of
//! a tenth of a period.
//!
//! - S: a force along y, its vy recorded 720 m along x (S broadside: the vs kernel);
//! - P: a force along x, its vx recorded on the same line (P along its axis: the vp kernel).
//!
//! 720 m is 4.2 dominant and 10.4 shortest S wavelengths (Ricker, 10 Hz), 2.4 and 6 for P.
//! The solver's gradient at a node, times the sample interval over the cell volume, is
//! compared with the kernel there, on the horizontal plane through source and receiver: in
//! the first Fresnel zone and on the whole plane, both outside a dominant S wavelength's
//! half around the source and the receiver, where the kernels are singular. And the
//! derivative along a smooth change (a Gaussian, σ = 60 m, midway), the number a line search
//! uses, against its exact value integrated on a finer grid.

mod common;

use std::time::Instant;

use common::stokes::{Medium, Problem, Pulse, Ricker, stokes, velocity_kernels};
use common::*;
use vibe_rations::gradient::{gradient, l2_misfit};

/// Source to receiver, and the margins beyond: every grid spacing here divides them.
const L: f64 = 720.0;
const PAD: f64 = 360.0;
const HALF: f64 = 360.0;
const PML: usize = 12;
const F0: f64 = 10.0;
const T0: f64 = 0.12;
const FORCE: f64 = 1e9;
const DURATION: f64 = 0.75;
/// The shortest S wavelength, at 2.5 f0, by which points per wavelength are counted.
const SHORTEST: f64 = VS / (2.5 * F0);
/// The dominant S wavelength.
const DOMINANT: f64 = VS / F0;
const BLOB: f64 = 60.0;

#[derive(Clone, Copy, Debug, PartialEq)]
enum Wave {
    S,
    P,
}

impl Wave {
    fn axis(self) -> usize {
        match self {
            Wave::S => 1,
            Wave::P => 0,
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct Scheme {
    order: usize,
    courant: f64,
    transform: bool,
    coefficients: &'static str,
}

const DEFAULT: Scheme = Scheme {
    order: 8,
    courant: 0.85,
    transform: true,
    coefficients: "optimised",
};

fn medium() -> Medium {
    Medium {
        vp: VP,
        vs: VS,
        rho: RHO,
    }
}

fn source(h: f64) -> [f64; 3] {
    let o = PML as f64 * h;
    [o + PAD, o + HALF, o + HALF]
}

fn problem(wave: Wave, h: f64) -> Problem {
    let s = source(h);
    Problem {
        medium: medium(),
        observed: Medium {
            vp: 1.02 * VP,
            vs: 1.02 * VS,
            rho: RHO,
        },
        source: s,
        force: wave.axis(),
        amplitude: FORCE,
        f0: F0,
        t0: T0,
        receiver: [s[0] + L, s[1], s[2]],
        component: wave.axis(),
        duration: DURATION,
    }
}

fn config(wave: Wave, h: f64, s: Scheme) -> String {
    let nx = ((2.0 * PAD + L) / h).round() as usize + 1 + 2 * PML;
    let ny = (2.0 * HALF / h).round() as usize + 1 + 2 * PML;
    let p = problem(wave, h);
    let force = match wave {
        Wave::S => format!("[0.0, {FORCE}, 0.0]"),
        Wave::P => format!("[{FORCE}, 0.0, 0.0]"),
    };
    let component = ["vx", "vy", "vz"][wave.axis()];
    format!(
        r#"
        [grid]
        points = [{nx}, {ny}, {ny}]
        spacing = {h}
        [time]
        duration = {DURATION}
        courant = {courant}
        [boundary]
        pml_width = {PML}
        free_surface = false
        [scheme]
        order = {order}
        remove_time_dispersion = {transform}
        coefficients = "{coefficients}"
        [[model.layer]]
        top = 0.0
        vp = {VP}
        vs = {VS}
        rho = {RHO}
        [[source]]
        position = {source:?}
        force = {force}
        wavelet = {{ kind = "ricker", frequency = {F0} }}
        [receivers]
        components = ["{component}"]
        [[receivers.point]]
        position = {receiver:?}
        "#,
        courant = s.courant,
        order = s.order,
        transform = s.transform,
        coefficients = s.coefficients,
        source = p.source,
        receiver = p.receiver,
    )
}

/// One gradient run: the solver's gradient for vp, vs, ρ on the plane through source and
/// receiver, as kernels (times the sample interval over the cell volume), and the derivative
/// along the Gaussian for each.
struct Outcome {
    points: Vec<[f64; 3]>,
    kernels: Vec<[f64; 3]>,
    blob: [f64; 3],
    /// The nodes where the Gaussian counts, for the exact derivative as the grid samples it.
    blob_nodes: Vec<[f64; 3]>,
    /// The recorded trace's misfit against Stokes's solution.
    trace_error: f64,
    /// The recorded residual's against the exact residual.
    residual_error: f64,
    seconds: f64,
}

/// The smooth change: a Gaussian of width `BLOB`, cut off 4 widths out (at e⁻⁸), which keeps
/// it clear of the source and the receiver.
fn blob(x: [f64; 3], centre: [f64; 3]) -> f64 {
    let d = sub(x, centre);
    let r2 = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]) / (BLOB * BLOB);
    if r2 > 16.0 { 0.0 } else { (-0.5 * r2).exp() }
}

fn midpoint(p: &Problem) -> [f64; 3] {
    std::array::from_fn(|a| 0.5 * (p.source[a] + p.receiver[a]))
}

/// `exact_residual`: the adjoint source is the exact residual rather than the recorded
/// trace's, to separate the error the fields carry into the kernel from the residual's.
fn run(wave: Wave, h: f64, scheme: Scheme, exact_residual: bool) -> Outcome {
    let p = problem(wave, h);
    let mut sim = setup(&config(wave, h, scheme));
    let times = sim.sample_times();
    let interval = times[1] - times[0];
    let observed: Vec<f32> = times.iter().map(|&t| p.observed(t) as f32).collect();
    let residual: Vec<f32> = times.iter().map(|&t| p.residual(0, t) as f32).collect();
    let start = Instant::now();
    let result = gradient(
        &mut sim,
        |r| {
            if exact_residual {
                let misfit = 0.5 * residual.iter().map(|&x| (x as f64).powi(2)).sum::<f64>();
                Ok((misfit, vec![residual.clone()]))
            } else {
                l2_misfit(&r.traces, std::slice::from_ref(&observed), &[1.0])
            }
        },
        None,
    )
    .expect("a gradient");
    let seconds = start.elapsed().as_secs_f64();
    let recorded = trace(&sim, sim.receivers.components[0], 0);
    let exact: Vec<f64> = times.iter().map(|&t| p.synthetic(t)).collect();
    let trace_error = misfit(&recorded, &exact);
    let recorded_residual: Vec<f64> = recorded
        .iter()
        .zip(&observed)
        .map(|(&s, &d)| s - d as f64)
        .collect();
    let exact_residual: Vec<f64> = residual.iter().map(|&r| r as f64).collect();
    let residual_error = misfit(&recorded_residual, &exact_residual);

    let grid = sim.solver.grid;
    let interior = sim.solver.interior();
    let g = &result.gradient;
    let m = medium();
    let scale = interval / grid.cell_volume();
    let at = |i: usize, j: usize, k: usize| {
        let n = grid.interior_idx(i, j, k);
        velocity_kernels(&m, [g.rho[n], g.lambda[n], g.mu[n]]).map(|x| x * scale)
    };
    let kc = (p.source[2] / h).round() as usize;
    let mut points = Vec::new();
    let mut kernels = Vec::new();
    for j in interior[1].clone() {
        for i in interior[0].clone() {
            points.push(grid.position([i, j, kc]));
            kernels.push(at(i, j, kc));
        }
    }
    let centre = midpoint(&p);
    let mut along = [0.0; 3];
    let mut blob_nodes = Vec::new();
    for k in interior[2].clone() {
        for j in interior[1].clone() {
            for i in interior[0].clone() {
                let w = blob(grid.position([i, j, k]), centre);
                if w > 0.0 {
                    blob_nodes.push(grid.position([i, j, k]));
                    let x = at(i, j, k);
                    // Relative changes: vp(1 + δ), as the inversion's unknowns are.
                    let rel = [m.vp, m.vs, m.rho];
                    for a in 0..3 {
                        along[a] += x[a] * rel[a] * w * grid.cell_volume();
                    }
                }
            }
        }
    }
    Outcome {
        points,
        kernels,
        blob: along,
        blob_nodes,
        trace_error,
        residual_error,
        seconds,
    }
}

/// The exact derivative along the Gaussian, ∫ K m δ dV, on a grid `q` apart.
fn exact_blob(p: &Problem, q: f64) -> [f64; 3] {
    let centre = midpoint(p);
    let reach = (4.0 * BLOB / q).ceil() as i64;
    let mut points = Vec::new();
    for k in -reach..=reach {
        for j in -reach..=reach {
            for i in -reach..=reach {
                let x = [
                    centre[0] + i as f64 * q,
                    centre[1] + j as f64 * q,
                    centre[2] + k as f64 * q,
                ];
                points.push(x);
            }
        }
    }
    let m = medium();
    let kernels = p.kernels(&points, 1e-3);
    let rel = [m.vp, m.vs, m.rho];
    let mut along = [0.0; 3];
    for (x, k) in points.iter().zip(kernels) {
        let w = blob(*x, centre) * q * q * q;
        let v = velocity_kernels(&m, k);
        for a in 0..3 {
            along[a] += v[a] * rel[a] * w;
        }
    }
    along
}

/// Which points of the plane are compared.
#[derive(Clone, Copy)]
enum Region {
    /// The first Fresnel zone: a detour of less than half a dominant S wavelength.
    Fresnel,
    Plane,
}

fn selected(p: &Problem, x: [f64; 3], region: Region) -> bool {
    let (a, b) = (norm(sub(x, p.source)), norm(sub(x, p.receiver)));
    let clear = a.min(b) >= 0.5 * DOMINANT;
    clear
        && match region {
            Region::Fresnel => a + b - L < 0.5 * DOMINANT,
            Region::Plane => true,
        }
}

/// For each of vp, vs, ρ: ‖g - k‖/‖k‖ and the cosine between them.
fn compare(got: &[[f64; 3]], exact: &[[f64; 3]]) -> [(f64, f64); 3] {
    std::array::from_fn(|a| {
        let (mut gk, mut gg, mut kk, mut dd) = (0.0, 0.0, 0.0, 0.0);
        for (g, k) in got.iter().zip(exact) {
            gk += g[a] * k[a];
            gg += g[a] * g[a];
            kk += k[a] * k[a];
            dd += (g[a] - k[a]).powi(2);
        }
        ((dd / kk).sqrt(), gk / (gg * kk).sqrt())
    })
}

/// For each of vp, vs, ρ, (‖g - k‖/‖k‖, cosine): the solver's gradient against the exact
/// kernels at the nodes (`total`) and against the exact kernels as the grid samples the model
/// (`fields`: the error its wavefields cause), and the latter against the former (`sampling`).
struct Comparison {
    total: [(f64, f64); 3],
    fields: [(f64, f64); 3],
    sampling: [(f64, f64); 3],
}

struct Measured {
    fresnel: Comparison,
    plane: Comparison,
    /// The derivative along the Gaussian against the exact one (`blob`), against the exact
    /// one as the grid samples the model (`blob_fields`), and the latter against the former.
    blob: [f64; 3],
    blob_fields: [f64; 3],
    blob_sampling: [f64; 3],
    trace_error: f64,
    residual_error: f64,
    seconds: f64,
}

fn measure(
    wave: Wave,
    h: f64,
    scheme: Scheme,
    exact_residual: bool,
    blob_ref: [f64; 3],
) -> Measured {
    let p = problem(wave, h);
    let o = run(wave, h, scheme, exact_residual);
    let m = medium();
    let out = |region| {
        let (pts, got): (Vec<[f64; 3]>, Vec<[f64; 3]>) = o
            .points
            .iter()
            .zip(&o.kernels)
            .filter(|(x, _)| selected(&p, **x, region))
            .map(|(x, k)| (*x, *k))
            .unzip();
        let velocities = |k: Vec<[f64; 3]>| -> Vec<[f64; 3]> {
            k.into_iter().map(|k| velocity_kernels(&m, k)).collect()
        };
        let exact = velocities(p.kernels(&pts, 1e-3));
        let grid = velocities(p.grid_kernels(&pts, h, 1e-3));
        Comparison {
            total: compare(&got, &exact),
            fields: compare(&got, &grid),
            sampling: compare(&grid, &exact),
        }
    };
    let fresnel = out(Region::Fresnel);
    let plane = out(Region::Plane);
    let centre = midpoint(&p);
    let rel = [m.vp, m.vs, m.rho];
    let sampled = velocity_kernels(
        &m,
        p.grid_integral(&o.blob_nodes, h, |x| blob(x, centre), 1e-3),
    );
    let sampled: [f64; 3] = std::array::from_fn(|a| sampled[a] * rel[a]);
    Measured {
        fresnel,
        plane,
        blob: std::array::from_fn(|a| o.blob[a] / blob_ref[a] - 1.0),
        blob_fields: std::array::from_fn(|a| o.blob[a] / sampled[a] - 1.0),
        blob_sampling: std::array::from_fn(|a| sampled[a] / blob_ref[a] - 1.0),
        trace_error: o.trace_error,
        residual_error: o.residual_error,
        seconds: o.seconds,
    }
}

fn main_parameter(wave: Wave) -> (usize, &'static str) {
    match wave {
        Wave::S => (1, "vs"),
        Wave::P => (0, "vp"),
    }
}

fn header(wave: Wave) {
    let (_, name) = main_parameter(wave);
    println!(
        "{:<30} {:>6} {:>6} | {:>24} {:>6} {:>6} | {:>14} {:>6} | {:>23} {:>7} {:>7} | {:>13} | {:>4}",
        "",
        "trace",
        "resid.",
        format!("Fresnel zone, {name}: all"),
        "fields",
        "grid",
        format!("plane, {name}: all"),
        "fields",
        format!("Gaussian, {name}: all"),
        "fields",
        "grid",
        "ρ: all fields",
        "s"
    );
}

fn row(label: &str, m: &Measured, wave: Wave) {
    let (main, _) = main_parameter(wave);
    let pc = |x: f64| format!("{:>5.2}%", 100.0 * x);
    let signed = |x: f64| format!("{:>+6.2}%", 100.0 * x);
    let (e, c) = m.fresnel.total[main];
    println!(
        "{label:<30} {} {} | {:>15} {:>+8.5} {} {} | {:>5} {:>+8.5} {} | {:>23} {} {} | {} {} | {:>4.0}",
        pc(m.trace_error),
        pc(m.residual_error),
        pc(e),
        c,
        pc(m.fresnel.fields[main].0),
        pc(m.fresnel.sampling[main].0),
        pc(m.plane.total[main].0),
        m.plane.total[main].1,
        pc(m.plane.fields[main].0),
        signed(m.blob[main]),
        signed(m.blob_fields[main]),
        signed(m.blob_sampling[main]),
        pc(m.fresnel.total[2].0),
        pc(m.fresnel.fields[2].0),
        m.seconds
    );
}

const SPACINGS: [f64; 7] = [20.0, 18.0, 15.0, 12.0, 10.0, 9.0, 8.0];

fn reference(wave: Wave) -> [f64; 3] {
    exact_blob(&problem(wave, 10.0), 10.0)
}

/// The default scheme (order 8, optimised coefficients, C = 0.85, the time error removed) on
/// grids from 3.5 to 8.7 points per shortest S wavelength.
#[test]
#[ignore = "a study, not a check: run with `just gradient-accuracy`"]
fn by_grid_spacing() {
    for wave in [Wave::S, Wave::P] {
        let reference = reference(wave);
        println!("\n{wave:?}, the default scheme, by grid spacing");
        header(wave);
        for h in SPACINGS {
            let m = measure(wave, h, DEFAULT, false, reference);
            row(&format!("h = {h} m: {:.2} ppw", SHORTEST / h), &m, wave);
        }
    }
}

/// Stencil orders and coefficients at 4.6 points per shortest S wavelength.
#[test]
#[ignore = "a study, not a check: run with `just gradient-accuracy`"]
fn by_order() {
    let wave = Wave::S;
    let reference = reference(wave);
    let h = 15.0;
    println!(
        "\n{wave:?}, h = {h} m ({:.2} ppw), by order and coefficients",
        SHORTEST / h
    );
    header(wave);
    for coefficients in ["taylor", "optimised"] {
        for order in [4, 8, 12, 16] {
            let scheme = Scheme {
                order,
                coefficients,
                ..DEFAULT
            };
            let m = measure(wave, h, scheme, false, reference);
            row(&format!("order {order}, {coefficients}"), &m, wave);
        }
    }
}

/// The time step, with leapfrog's time error kept and removed (order 8, optimised).
#[test]
#[ignore = "a study, not a check: run with `just gradient-accuracy`"]
fn by_time_step() {
    let wave = Wave::S;
    let reference = reference(wave);
    let h = 15.0;
    println!(
        "\n{wave:?}, h = {h} m ({:.2} ppw), by time step",
        SHORTEST / h
    );
    header(wave);
    for transform in [false, true] {
        for courant in [0.85, 0.5, 0.3] {
            let scheme = Scheme {
                courant,
                transform,
                ..DEFAULT
            };
            let m = measure(wave, h, scheme, false, reference);
            let removed = if transform { "removed" } else { "kept" };
            row(&format!("C = {courant}, time error {removed}"), &m, wave);
        }
    }
}

/// The adjoint source from the recorded trace, as an inversion has it, against the exact
/// residual: how much of the error comes in through the residual.
#[test]
#[ignore = "a study, not a check: run with `just gradient-accuracy`"]
fn by_residual() {
    let wave = Wave::S;
    let reference = reference(wave);
    println!("\n{wave:?}, the default scheme, the residual recorded or exact");
    header(wave);
    for h in [20.0, 15.0, 10.0] {
        for exact in [false, true] {
            let m = measure(wave, h, DEFAULT, exact, reference);
            let which = if exact { "exact" } else { "recorded" };
            row(&format!("h = {h} m, {which}"), &m, wave);
        }
    }
}

/// The baseline `just gradient-accuracy` measures, at its coarsest grid (3.5 points per
/// shortest S wavelength, the default scheme): the error the wavefields put into the
/// gradient, against the exact kernels as the grid samples the model. Measured: S, 0.26% in
/// the Fresnel zone and +0.19% along the Gaussian; P, 0.12% and -0.10%.
#[test]
fn the_gradient_matches_the_exact_kernels() {
    for (wave, zone, along) in [(Wave::S, 0.0035, 0.0025), (Wave::P, 0.0016, 0.0013)] {
        let (main, name) = main_parameter(wave);
        let m = measure(wave, 20.0, DEFAULT, false, reference(wave));
        let (e, c) = m.fresnel.fields[main];
        println!(
            "{wave:?}: {name} off by {:.2}% in the Fresnel zone (cosine {c:.5}), {:+.2}% along the Gaussian",
            100.0 * e,
            100.0 * m.blob_fields[main]
        );
        assert!(e < zone, "{wave:?}: {e}");
        assert!(
            m.blob_fields[main].abs() < along,
            "{wave:?}: {}",
            m.blob_fields[main]
        );
    }
}

/// The analytic pieces against each other: Stokes's velocity against the physics tests',
/// its gradient against finite differences of its displacement, the tabulated adjoint pulse
/// against the residual it tabulates.
#[test]
fn the_reference_is_consistent() {
    let m = medium();
    let w = vibe_rations::source::Wavelet::Ricker { f0: F0, t0: T0 };
    let velocity = Ricker {
        f0: F0,
        t0: T0,
        amplitude: FORCE,
        level: 1,
    };
    for d in [[300.0, 40.0, -70.0], [-90.0, 500.0, 120.0]] {
        for t in [0.15, 0.2, 0.3, 0.4] {
            let a = stokes(&m, 0, d, t, &velocity).0;
            let b = force_velocity(&w, FORCE, d, t);
            let scale = norm(b).max(1e-12);
            assert!(
                norm(sub(a, b)) < 1e-5 * scale + 1e-15,
                "{d:?} {t}: {a:?} {b:?}"
            );
        }
    }
    let displacement = Ricker {
        level: 0,
        ..velocity
    };
    let e = 1e-3;
    for j in 0..3 {
        for d in [[200.0, 50.0, 30.0], [-40.0, 260.0, -90.0]] {
            for t in [0.17, 0.25] {
                let (_, grad) = stokes(&m, j, d, t, &displacement);
                let scale = grad.iter().flatten().fold(0.0f64, |a, &x| a.max(x.abs()));
                for l in 0..3 {
                    let mut dp = d;
                    let mut dm = d;
                    dp[l] += e;
                    dm[l] -= e;
                    let up = stokes(&m, j, dp, t, &displacement).0;
                    let um = stokes(&m, j, dm, t, &displacement).0;
                    for i in 0..3 {
                        let fd = (up[i] - um[i]) / (2.0 * e);
                        assert!(
                            (fd - grad[i][l]).abs() < 1e-6 * scale,
                            "j {j} d {d:?} t {t} ∂{l}u{i}: {fd} against {}",
                            grad[i][l]
                        );
                    }
                }
            }
        }
    }
    let p = problem(Wave::S, 10.0);
    let table = p.adjoint_pulse(1e-5);
    let peak = |k| (0..750).fold(0.0f64, |a, i| a.max(p.residual(k, i as f64 * 1e-3).abs()));
    let (peak0, peak1) = (peak(0), peak(1));
    for s in [0.1234567, 0.25, 0.3333, 0.41] {
        let [g, _, g1, _] = table.eval(s);
        assert!(
            (g + p.residual(1, DURATION - s)).abs() < 1e-7 * peak1,
            "h'({s})"
        );
        assert!(
            (g1 - p.residual(0, DURATION - s)).abs() < 1e-7 * peak0,
            "h({s})"
        );
    }
}

/// The record ends after the residual: what the adjoint source starts with is negligible.
#[test]
fn the_record_holds_the_residual() {
    for wave in [Wave::S, Wave::P] {
        let p = problem(wave, 10.0);
        let peak = (0..750).fold(0.0f64, |a, i| a.max(p.residual(0, i as f64 * 1e-3).abs()));
        let end = (0..10).fold(0.0f64, |a, i| {
            a.max(p.residual(0, DURATION - i as f64 * 1e-3).abs())
        });
        assert!(end < 1e-6 * peak, "{wave:?}: {end:e} of {peak:e}");
    }
}
