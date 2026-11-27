//! How accurate each stencil order and time step is over FWI-like distances: `just accuracy`.
//!
//! Not part of the regular suite (it takes a few minutes); it prints tables of misfits
//! against the closed-form solutions at coarse grids where numerical dispersion shows: for
//! every combination of stencil order and Courant fraction with leapfrog's time error kept,
//! then with it removed (the default), and the free-surface test on finer grids.
//!
//! - P: an explosion, 4 grid points per shortest P wavelength;
//! - S: a point force seen broadside (S in the far field), 3.5, 3.9 and 4.6 points per
//!   shortest S wavelength.
//!
//! Receivers sit 300 to 1800 m from the source: 1.5–9 dominant P wavelengths, 1.7–10 S
//! wavelengths. The README's "Accuracy" section summarises the results.

mod common;

use std::time::Instant;

use common::*;
use vibe_rations::receiver::Component;
use vibe_rations::simulation::Simulation;
use vibe_rations::source::Wavelet;
use vibe_rations::tdt::Warp;

const DISTANCES: [f64; 4] = [300.0, 900.0, 1500.0, 1800.0];
const ORDERS: [usize; 4] = [4, 8, 12, 16];
const COURANTS: [f64; 4] = [0.85, 0.5, 0.4, 0.3];

enum Wave {
    P,
    S,
}

/// A line of receivers along x from a source, in a domain long in x and narrow across.
fn config(wave: &Wave, h: f64, order: usize, courant: f64) -> String {
    config_with(wave, h, order, courant, false)
}

fn config_with(wave: &Wave, h: f64, order: usize, courant: f64, transform: bool) -> String {
    config_full(wave, h, order, courant, transform, "taylor")
}

fn config_full(
    wave: &Wave,
    h: f64,
    order: usize,
    courant: f64,
    transform: bool,
    coefficients: &str,
) -> String {
    let pml = 12;
    let pad = (pml + 4) as f64 * h;
    let x0 = pad;
    let nx = ((x0 + DISTANCES[3] + pad) / h).round() as usize + 1;
    let n_across = 2 * pml + 17;
    let c = (n_across / 2) as f64 * h;
    let receivers: String = DISTANCES
        .iter()
        .map(|d| format!("[[receivers.point]]\nposition = [{}, {c}, {c}]\n", x0 + d))
        .collect();
    let (source, wavelet, duration, component) = match wave {
        Wave::P => ("explosion = 1e13", 15.0, 0.8, "vx"),
        // A force along y, seen along x: shear in the far field.
        Wave::S => ("force = [0.0, 1e9, 0.0]", 10.0, 1.3, "vy"),
    };
    format!(
        r#"
        [grid]
        points = [{nx}, {n_across}, {n_across}]
        spacing = {h}
        [time]
        duration = {duration}
        courant = {courant}
        [boundary]
        pml_width = {pml}
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
        position = [{x0}, {c}, {c}]
        {source}
        wavelet = {{ kind = "ricker", frequency = {wavelet} }}
        [receivers]
        components = ["{component}"]
        {receivers}
        "#
    )
}

/// The misfit at each distance, and the speed in Mcell updates per second.
fn measure(wave: &Wave, h: f64, order: usize, courant: f64) -> (Vec<f64>, f64) {
    let (sim, speed) = simulate_timed(&config(wave, h, order, courant));
    (misfits(wave, &sim, false), speed)
}

fn simulate_timed(toml: &str) -> (Simulation, f64) {
    let mut sim: Simulation = setup(toml);
    let start = Instant::now();
    sim.run(None, false).expect("a stable run");
    let seconds = start.elapsed().as_secs_f64();
    let speed = sim.solver.grid.interior_len() as f64 * sim.steps as f64 / seconds / 1e6;
    (sim, speed)
}

/// The misfit at each distance; `band_limit` first limits the traces to the band the
/// time-dispersion transform keeps, to compare like with like.
fn misfits(wave: &Wave, sim: &Simulation, band_limit: bool) -> Vec<f64> {
    let times = sim.sample_times();
    let dt = sim.solver.dt;
    let limit = Warp::band_limit(0.5 * dt, dt, times.len(), 2.0 * sim.info.max_frequency);
    DISTANCES
        .iter()
        .enumerate()
        .map(|(r, &d)| {
            let exact: Vec<f64> = match wave {
                Wave::P => {
                    let w = Wavelet::Ricker {
                        f0: 15.0,
                        t0: 1.2 / 15.0,
                    };
                    times
                        .iter()
                        .map(|&t| explosion_velocity(&w, 1e13, d, t))
                        .collect()
                }
                Wave::S => {
                    let w = Wavelet::Ricker { f0: 10.0, t0: 0.12 };
                    // A force along y seen along x is a force along x seen along y.
                    times
                        .iter()
                        .map(|&t| force_velocity(&w, 1e9, [0.0, d, 0.0], t)[0])
                        .collect()
                }
            };
            let c = match wave {
                Wave::P => Component::Vx,
                Wave::S => Component::Vy,
            };
            let mut got = trace(sim, c, r);
            if band_limit {
                got = limit.apply(&got);
            }
            misfit(&got, &exact)
        })
        .collect()
}

#[test]
#[ignore = "a study, not a check: run with `just accuracy`"]
fn stencil_order_and_time_step() {
    for (wave, h, label) in [
        (Wave::P, 20.0, "P, 4.0 points per wavelength"),
        (Wave::S, 20.0, "S, 3.5 points per wavelength"),
        (Wave::S, 15.0, "S, 4.6 points per wavelength"),
    ] {
        println!("\n{label}: L2 misfit (%) at {DISTANCES:?} m");
        println!("order  courant  {:>32}  Mcell/s", "");
        for order in ORDERS {
            for courant in COURANTS {
                let (m, speed) = measure(&wave, h, order, courant);
                let cells: Vec<String> = m.iter().map(|x| format!("{:7.2}", 100.0 * x)).collect();
                println!("{order:5}  {courant:7}  {}  {speed:7.0}", cells.join(" "));
            }
        }
    }
}

/// How fine a grid the fourth-order scheme needs to match the eighth-order one at Courant 0.3
/// on a 15 m grid (S misfit about 1.3% after ten wavelengths).
#[test]
#[ignore = "a study, not a check: run with `just accuracy`"]
fn fourth_order_on_finer_grids() {
    println!("\nS, order 4 on finer grids: L2 misfit (%) at {DISTANCES:?} m");
    println!("    h  ppw  courant  {:>32}  Mcell/s", "");
    for h in [10.0, 7.5, 5.0] {
        for courant in [0.85, 0.5] {
            let (m, speed) = measure(&Wave::S, h, 4, courant);
            let cells: Vec<String> = m.iter().map(|x| format!("{:7.2}", 100.0 * x)).collect();
            let ppw = VS / 25.0 / h;
            println!(
                "{h:5} {ppw:4.1}  {courant:7}  {}  {speed:7.0}",
                cells.join(" ")
            );
        }
    }
}

/// Removing leapfrog's time error (`remove_time_dispersion`): the transformed runs at large
/// time steps against today's default, C = 0.3, with its traces as they are and limited to the
/// band the transform keeps.
#[test]
#[ignore = "a study, not a check: run with `just accuracy`"]
fn removing_time_dispersion() {
    let row = |label: &str, m: &[f64], speed: f64| {
        let cells: Vec<String> = m.iter().map(|x| format!("{:7.2}", 100.0 * x)).collect();
        println!("{label:24} {}  {speed:7.0}", cells.join(" "));
    };
    for (wave, h, label) in [
        (Wave::P, 20.0, "P, 4.0 points per wavelength"),
        (Wave::S, 20.0, "S, 3.5 points per wavelength"),
        // The spacing divides the distances, so the receivers sit on nodes.
        (Wave::S, 300.0 / 17.0, "S, 3.9 points per wavelength"),
        (Wave::S, 15.0, "S, 4.6 points per wavelength"),
    ] {
        println!("\n{label}, time dispersion removed: L2 misfit (%) at {DISTANCES:?} m");
        println!("{:24} {:>32}  Mcell/s", "order, courant", "");
        for order in ORDERS {
            let (plain, speed) = simulate_timed(&config(&wave, h, order, 0.3));
            row(
                &format!("{order:2}, 0.3"),
                &misfits(&wave, &plain, false),
                speed,
            );
            row(
                &format!("{order:2}, 0.3, band-limited"),
                &misfits(&wave, &plain, true),
                speed,
            );
            for courant in [0.85, 0.95] {
                let (sim, speed) = simulate_timed(&config_with(&wave, h, order, courant, true));
                row(
                    &format!("{order:2}, {courant}, transformed"),
                    &misfits(&wave, &sim, false),
                    speed,
                );
            }
        }
    }
}

/// The free-surface test (`the_free_surface_doubles_a_vertically_incident_p_wave`) on finer
/// grids: whether its misfit is the scheme's, which shrinks with the spacing, or the
/// reference's (doubling is exact only for plane waves), which doesn't.
#[test]
#[ignore = "a study, not a check: run with `just accuracy`"]
fn free_surface_on_finer_grids() {
    let depth = 300.0;
    let w = Wavelet::Ricker {
        f0: 15.0,
        t0: 1.2 / 15.0,
    };
    println!("\nP doubling at the free surface, 300 m above an explosion: L2 misfit (%)");
    println!("    h  order  transform  misfit  Mcell/s");
    for h in [10.0, 7.5, 5.0] {
        for (order, transform) in [(8, true), (8, false), (16, true)] {
            let n = |metres: f64| (metres / h).round() as usize + 1;
            let (sim, speed) = simulate_timed(&format!(
                r#"
                [grid]
                points = [{}, {}, {}]
                spacing = {h}
                [time]
                duration = 0.3
                [boundary]
                pml_width = {}
                [scheme]
                order = {order}
                remove_time_dispersion = {transform}
                coefficients = "taylor"
                [[model.layer]]
                top = 0.0
                vp = {VP}
                vs = {VS}
                rho = {RHO}
                [[source]]
                position = [400.0, 400.0, {depth}]
                explosion = 1e13
                wavelet = {{ kind = "ricker", frequency = 15.0 }}
                [receivers]
                components = ["vz"]
                [[receivers.point]]
                position = [400.0, 400.0, 0.0]
                "#,
                n(800.0),
                n(800.0),
                n(600.0),
                (120.0 / h).round() as usize,
            ));
            let exact: Vec<f64> = sim
                .sample_times()
                .iter()
                .map(|&t| -2.0 * explosion_velocity(&w, 1e13, depth, t))
                .collect();
            let e = misfit(&trace(&sim, Component::Vz, 0), &exact);
            println!(
                "{h:5}  {order:5}  {transform:9}  {:6.2}  {speed:7.0}",
                100.0 * e
            );
        }
    }
}

/// Lamb's problem (`common::lamb::survey`) on finer grids: how the free surface's error in
/// the Rayleigh wave falls with the spacing (order 8, 10 m to 2.5 m), and whether the
/// stencil's order or the time-dispersion transform matter to it (at 10 and 5 m). About two
/// and a half minutes, most of it the 2.5 m grid's 11 million nodes.
#[test]
#[ignore = "a study, not a check: run with `just accuracy`"]
fn lambs_problem_on_finer_grids() {
    use common::lamb::survey;
    println!(
        "\nLamb's problem, vz on the surface against Pekeris: misfit (%), amplitude ⟨v, w⟩/⟨w, w⟩, \
         the largest swing's lateness (samples) and size"
    );
    println!("    h  scheme                    r (m)  misfit  amplitude  shift   peak  Mcell/s");
    for h in [10.0, 7.5, 5.0, 3.75, 2.5] {
        let mut variants = vec![("order 8", "", "")];
        if h == 10.0 || h == 5.0 {
            variants.push(("order 16", "order = 16", ""));
            variants.push((
                "order 8, time error kept",
                "remove_time_dispersion = false",
                "courant = 0.3",
            ));
        }
        for (label, scheme, time) in variants {
            let (sim, speed) = simulate_timed(&survey::config(h, scheme, time));
            for [r, e, amplitude, shift, peak] in survey::compare(&sim) {
                println!(
                    "{h:5}  {label:24} {r:6.1}  {:6.2}  {amplitude:9.4}  {shift:5.2}  {peak:.4}  {speed:7.0}",
                    100.0 * e
                );
            }
        }
    }
}

/// Optimised coefficients (`coefficients = "optimised"`) against Taylor's, both with the time
/// error removed at C = 0.85, down to three points per shortest wavelength. The spacings divide
/// the distances, so the receivers sit on nodes.
#[test]
#[ignore = "a study, not a check: run with `just accuracy`"]
fn optimised_coefficients() {
    for (wave, h, label) in [
        (Wave::P, 300.0 / 11.0, "P, 2.9 points per wavelength"),
        (Wave::P, 300.0 / 13.0, "P, 3.5 points per wavelength"),
        (Wave::P, 20.0, "P, 4.0 points per wavelength"),
        (Wave::S, 300.0 / 13.0, "S, 3.0 points per wavelength"),
        (Wave::S, 20.0, "S, 3.5 points per wavelength"),
        (Wave::S, 300.0 / 17.0, "S, 3.9 points per wavelength"),
        (Wave::S, 15.0, "S, 4.6 points per wavelength"),
    ] {
        println!("\n{label}, time dispersion removed: L2 misfit (%) at {DISTANCES:?} m");
        println!("{:24} {:>32}  Mcell/s", "order, coefficients", "");
        for order in [8, 12, 16] {
            for coefficients in ["taylor", "optimised"] {
                let (sim, speed) =
                    simulate_timed(&config_full(&wave, h, order, 0.85, true, coefficients));
                let cells: Vec<String> = misfits(&wave, &sim, false)
                    .iter()
                    .map(|x| format!("{:7.2}", 100.0 * x))
                    .collect();
                println!(
                    "{:24} {}  {speed:7.0}",
                    format!("{order:2}, {coefficients}"),
                    cells.join(" ")
                );
            }
        }
    }
}

/// Does the free surface's amplitude excess on Lamb's problem come from the kink of the
/// evenly mirrored vz? Above the surface vz is mirrored evenly, so continued across it vz has
/// a kink, and interpolating vz onto the surface reads it as if from c·h below (c, the mean
/// depth of the receiver's weights over h, is 0.34 at order 8). A force on the surface, the
/// receiver's transpose, acts from the same depth. A Rayleigh wave's vz grows with depth near
/// the surface (∂z w / w = 0.227 k for ν = ¼), so each end on the surface over-reads by about
/// 0.227 k c h: at a fixed h, in proportion to frequency. vz itself lives at h/2, so a source
/// and receivers there interpolate nothing along z and should show none of it.
///
/// At h = 5 m, at several frequencies: (a) source and receivers on the surface, as in the
/// Lamb test; (b) the source on the surface, the receivers at h/2; (c) both at h/2. The exact
/// answers at h/2 are Pekeris's surface traces taken down by the Rayleigh wave's depth
/// dependence for ν = ¼, once per end below the surface: exact for the Rayleigh wave, which
/// carries the swing, an approximation for the much weaker body waves. Each end's own excess
/// comes from pairs of runs that differ in it alone, receiver by receiver, so whatever else
/// the free surface gets wrong cancels: the receiver's, (a) - (b); the source's, (b) - (c).
/// (a′) and (b′) are (a) and (b) with `[boundary] surface_slope`, which continues the mirrored
/// fields with their slope instead: the same ends again, primed.
#[test]
#[ignore = "a study, not a check: run with `just accuracy`"]
fn the_free_surface_excess_by_frequency_and_depth() {
    use common::lamb::{Lamb, gamma, survey};
    use common::stokes::Ricker;
    use std::f64::consts::PI;
    use vibe_rations::solver::Field;

    let h = 5.0;
    let vs = survey::vs();
    let c_r = vs / gamma();
    // The Rayleigh wave's vertical motion at depth z over that at the surface, for k z = x:
    // (q e^(-q x) - a e^(-s x)) / (q - a), q and s its P and S decay rates over k; its slope
    // at the surface, over k.
    let q = (1.0 - (c_r / survey::VP).powi(2)).sqrt();
    let s = (1.0 - (c_r / vs).powi(2)).sqrt();
    let a = 2.0 * q / (1.0 + s * s);
    let rayleigh = |x: f64| (q * (-q * x).exp() - a * (-s * x).exp()) / (q - a);
    let slope = (a * s - q * q) / (q - a);
    // A trace taken down to each of `depths` in turn, frequency by frequency (a DFT, zero
    // padded to twice the length).
    let deepen = |x: &[f64], dt: f64, depths: &[f64]| -> Vec<f64> {
        let (n, m) = (x.len(), 2 * x.len());
        let twiddle: Vec<(f64, f64)> = (0..m)
            .map(|j| {
                let angle = -2.0 * PI * j as f64 / m as f64;
                (angle.cos(), angle.sin())
            })
            .collect();
        let spectrum: Vec<(f64, f64)> = (0..m)
            .map(|k| {
                let (mut re, mut im) = (0.0, 0.0);
                for (t, &v) in x.iter().enumerate() {
                    let (c, s) = twiddle[(k * t) % m];
                    re += v * c;
                    im += v * s;
                }
                let f = k.min(m - k) as f64 / (m as f64 * dt);
                let gain: f64 = depths
                    .iter()
                    .map(|&z| rayleigh(2.0 * PI * f / c_r * z))
                    .product();
                (re * gain, im * gain)
            })
            .collect();
        (0..n)
            .map(|t| {
                let sum: f64 = spectrum
                    .iter()
                    .enumerate()
                    .map(|(k, &(re, im))| {
                        let (c, s) = twiddle[(k * t) % m];
                        re * c + im * s
                    })
                    .sum();
                sum / m as f64
            })
            .collect()
    };
    let peak = |v: &[f64]| {
        let k = (1..v.len() - 1)
            .max_by(|&a, &b| v[a].abs().total_cmp(&v[b].abs()))
            .unwrap();
        let (a, b, c) = (v[k - 1].abs(), v[k].abs(), v[k + 1].abs());
        let d = 0.5 * (a - c) / (a - 2.0 * b + c);
        b - 0.25 * (a - c) * d
    };
    let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
    let range = |v: &[f64]| {
        let lo = v.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = v.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        format!("{:+6.2} to {:+6.2}%", 100.0 * lo, 100.0 * hi)
    };

    println!(
        "\nThe free surface's excess on Lamb's problem at h = {h} m, the four receivers' range: \
         the largest swing's size and the amplitude ⟨v, e⟩/⟨e, e⟩ against the exact answer e, \
         and the misfit"
    );
    println!(
        "  f0   kh_R  case                         swing              amplitude          misfit"
    );
    let mut ends = Vec::new();
    for f0 in [7.5, 10.0, 15.0, 20.0, 30.0] {
        let lamb = Lamb {
            vs,
            rho: survey::RHO,
            force: 1e9,
            wavelet: Ricker {
                f0,
                t0: 1.2 / f0,
                amplitude: 1.0,
                level: 0,
            },
        };
        let kh = 2.0 * PI * f0 / c_r * h;
        // The Rayleigh wave past the farthest receiver, and the wavelet's length.
        let duration = 330.0 / c_r + 2.4 / f0 + 0.03;
        // [swing, amplitude] excess per receiver, per case; the surface receiver's c.
        let mut cases: Vec<Vec<[f64; 2]>> = Vec::new();
        let mut c = 0.0;
        for (label, source_z, receiver_z, slope) in [
            ("(a) both on the surface", 0.0, 0.0, false),
            ("(b) receivers at h/2", 0.0, 0.5 * h, false),
            ("(c) both at h/2", 0.5 * h, 0.5 * h, false),
            ("(a′) (a), surface slope", 0.0, 0.0, true),
            ("(b′) (b), surface slope", 0.0, 0.5 * h, true),
        ] {
            let n = |metres: f64| (metres / h).round() as usize + 1;
            let at = |p: [f64; 3], z: f64| [p[0], p[1], z];
            let receivers: String = survey::RECEIVERS
                .iter()
                .map(|&p| format!("[[receivers.point]]\nposition = {:?}\n", at(p, receiver_z)))
                .collect();
            let toml = format!(
                r#"
                [grid]
                points = [{}, {}, {}]
                spacing = {h}
                [time]
                duration = {duration}
                [boundary]
                pml_width = {}
                free_surface = true
                surface_slope = {slope}
                [[model.layer]]
                top = 0.0
                vp = {}
                vs = {vs}
                rho = {}
                [[source]]
                position = {:?}
                force = [0.0, 0.0, 1e9]
                wavelet = {{ kind = "ricker", frequency = {f0} }}
                [receivers]
                components = ["vz"]
                {receivers}
                "#,
                n(750.0),
                n(630.0),
                n(360.0),
                (120.0 / h).round() as usize,
                survey::VP,
                survey::RHO,
                at(survey::SOURCE, source_z),
            );
            let (sim, _) = simulate_timed(&toml);
            if receiver_z == 0.0 && !slope {
                // The mean depth, over h, of the vz values a surface receiver reads (vz lives
                // at (k + ½) h; values above the surface fold onto those below).
                let g = sim.solver.grid;
                c = sim.receivers.ops[0][0]
                    .terms
                    .iter()
                    .filter(|t| t.0 == Field::Vz)
                    .map(|&(_, i, w)| {
                        let k = (i / g.plane_len()) as f64 - g.halo as f64;
                        w as f64 * (k + 0.5)
                    })
                    .sum::<f64>();
            }
            let times = sim.sample_times();
            let dt = times[1] - times[0];
            let mut excess = Vec::new();
            let mut misfits = Vec::new();
            for (i, p) in survey::RECEIVERS.iter().enumerate() {
                let r = norm(sub(*p, survey::SOURCE));
                let surface: Vec<f64> = times.iter().map(|&t| lamb.velocity(r, t)).collect();
                let exact = deepen(&surface, dt, &[source_z, receiver_z]);
                let v = trace(&sim, Component::Vz, i);
                excess.push([
                    peak(&v) / peak(&exact) - 1.0,
                    dot(&v, &exact) / dot(&exact, &exact) - 1.0,
                ]);
                misfits.push(misfit(&v, &exact));
            }
            let column = |j: usize| excess.iter().map(|e| e[j]).collect::<Vec<_>>();
            println!(
                "{f0:5} {kh:5.3}  {label:26}  {}  {}  {:5.1}–{:4.1}%",
                range(&column(0)),
                range(&column(1)),
                100.0 * misfits.iter().copied().fold(f64::INFINITY, f64::min),
                100.0 * misfits.iter().copied().fold(f64::NEG_INFINITY, f64::max),
            );
            cases.push(excess);
        }
        ends.push((f0, kh, c, cases));
    }
    println!(
        "\nEach end's own excess, receiver by receiver: the receiver's (a) - (b), the source's \
         (b) - (c); per unit kh_R; and 0.227 c kh_R at f0, the first-order estimate"
    );
    println!(
        "  f0   kh_R   end        swing             per kh   amplitude          per kh   estimate"
    );
    for (f0, kh, c, cases) in ends {
        for (label, x, y) in [
            ("receiver", 0, 1),
            ("source", 1, 2),
            ("receiver′", 3, 4),
            ("source′", 4, 2),
        ] {
            let diff = |j: usize| {
                cases[x]
                    .iter()
                    .zip(&cases[y])
                    .map(|(p, q)| p[j] - q[j])
                    .collect::<Vec<_>>()
            };
            let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
            let (swing, amplitude) = (diff(0), diff(1));
            println!(
                "{f0:5} {kh:5.3}  {label:9} {}  {:5.1}   {}  {:5.1}   {:+5.2}%",
                range(&swing),
                100.0 * mean(&swing) / kh,
                range(&amplitude),
                100.0 * mean(&amplitude) / kh,
                100.0 * slope * c * kh,
            );
        }
    }
}

/// The Born run and the Hessian against Lamb's problem (`just accuracy
/// hessian_against_lambs_problem`). For uniform changes of the half-space, of its speeds (vp
/// and vs together, so that it stays a Poisson solid) and of its density, Pekeris's closed
/// form gives the traces' derivatives and the L2 misfit's 2 × 2 Hessian,
/// χ″_ij = w ⟨∂_i d, ∂_j d⟩ + w ⟨d - d_obs, ∂_ij d⟩, here against observed data from Pekeris
/// for a solid 2% faster and 1% denser (w = 1/‖d_obs‖²). The solver's Born run and Hessian
/// products (every node changed, the layers' too) are set beside it, and beside the solver's
/// own derivatives along the same directions: central differences of its traces for the
/// speeds (extrapolated, Δt, the layers and the operator held), and the exact scaling for
/// density (d ∝ 1/ρ at fixed speeds, in the scheme as in the solid). So the error the free
/// surface puts into the traces can be followed into J and H, at 10 m and 5 m, and at 5 m
/// with `surface_slope`; and again with the solver's own data for the faster, denser solid,
/// which leaves its error out of the residual. About 20 minutes.
///
/// Measured (the research log's "The Hessian", 2026-11-27): the speeds' derivative is 1.1–1.3
/// times as far from Pekeris's as the traces are (17.7–34.7% at 10 m, 6.8–11.5% at 5 m,
/// 5.4–10.6% with the slope); the Born run and the products are the solver's own to 0.8–2.1e-3
/// and 2–6e-4, the layers' share (density's Born run 1.1–1.5e-2 of its smaller derivative).
/// The Gauss–Newton entries are off by about the traces' amplitude error squared (+7.8 to
/// +9.4% at 5 m, -2.9 to -0.5% with the slope); the second-order ones, which correlate the
/// residual, by 8–70% against Pekeris's data, whose residual holds the solver's error, and
/// by 0.6–15% against its own.
#[test]
#[ignore]
fn hessian_against_lambs_problem() {
    use common::lamb::{Lamb, survey};
    use vibe_rations::born::{Change, born};
    use vibe_rations::config::Config;
    use vibe_rations::hessian::hessian_product;
    use vibe_rations::model::Model;

    const FASTER: f64 = 0.02;
    const DENSER: f64 = 0.01;
    let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
    let rel = |a: f64, b: f64| (a - b) / b.abs();
    for (h, slope) in [(10.0, false), (5.0, false), (5.0, true)] {
        let started = Instant::now();
        let label = format!("{h} m{}", if slope { ", surface slope" } else { "" });
        let text = survey::config(h, "", "");
        let text = if slope {
            text.replace(
                "free_surface = true",
                "free_surface = true\nsurface_slope = true",
            )
        } else {
            text
        };
        let mut cfg = Config::from_toml(&text).unwrap();
        let base = std::path::Path::new(".");
        let model = Model::from_config(Simulation::grid(&cfg).unwrap(), &cfg.model, base).unwrap();
        let sim = Simulation::with_model(&cfg, model.clone()).unwrap();
        cfg.time.dt = Some(sim.solver.dt);
        cfg.time.courant = None;
        cfg.boundary.pml_vp_max = Some(model.vp_range()[1]);
        cfg.scheme.optimised_for = Some(sim.info.points_per_wavelength.max(2.0 / 0.9));
        let times = sim.sample_times();
        drop(sim);

        // Pekeris's traces for the speeds and density scaled by 1 + a and 1 + b, receiver by
        // receiver, at the solver's sample times.
        let exact = |a: f64, b: f64| -> Vec<f64> {
            let lamb = Lamb {
                vs: survey::vs() * (1.0 + a),
                rho: survey::RHO * (1.0 + b),
                ..survey::lamb()
            };
            survey::RECEIVERS
                .iter()
                .flat_map(|p| {
                    let r = norm(sub(*p, survey::SOURCE));
                    times
                        .iter()
                        .map(move |&t| lamb.velocity(r, t))
                        .collect::<Vec<_>>()
                })
                .collect()
        };
        let observed = exact(FASTER, DENSER);
        let w = 1.0 / dot(&observed, &observed);
        let observed32: Vec<f32> = observed.iter().map(|&v| v as f32).collect();
        // Five-point differences in a, from traces at a = ±s, ±2s.
        let differences = |at: &dyn Fn(f64) -> Vec<f64>, s: f64| {
            let (m2, m1, z, p1, p2) = (at(-2.0 * s), at(-s), at(0.0), at(s), at(2.0 * s));
            let first: Vec<f64> = (0..z.len())
                .map(|i| (8.0 * (p1[i] - m1[i]) - (p2[i] - m2[i])) / (12.0 * s))
                .collect();
            let second: Vec<f64> = (0..z.len())
                .map(|i| {
                    (-p2[i] + 16.0 * p1[i] - 30.0 * z[i] + 16.0 * m1[i] - m2[i]) / (12.0 * s * s)
                })
                .collect();
            (z, first, second)
        };
        // The 2 × 2 Hessian, speeds then density, from d, ∂_a d and ∂_aa d, with density's
        // derivatives from d ∝ 1/ρ: ∂_b d = -d, ∂_bb d = 2d, ∂_ab d = -∂_a d.
        let hessian = |d: &[f64], da: &[f64], daa: &[f64], observed: &[f64]| {
            let r: Vec<f64> = d.iter().zip(observed).map(|(d, o)| d - o).collect();
            let db: Vec<f64> = d.iter().map(|v| -v).collect();
            let dab: Vec<f64> = da.iter().map(|v| -v).collect();
            let dbb: Vec<f64> = d.iter().map(|v| 2.0 * v).collect();
            let gn = [
                [w * dot(da, da), w * dot(da, &db)],
                [w * dot(&db, da), w * dot(&db, &db)],
            ];
            let second = [
                [w * dot(&r, daa), w * dot(&r, &dab)],
                [w * dot(&r, &dab), w * dot(&r, &dbb)],
            ];
            (gn, second)
        };
        let (d, da, daa) = differences(&|a| exact(a, 0.0), 1e-3);
        let (gn_exact, second_exact) = hessian(&d, &da, &daa, &observed);

        // The solver's traces for the speeds scaled by 1 + a everywhere.
        let scaled = |a: f64| {
            let mut m = model.clone();
            m.vp.iter_mut()
                .for_each(|v| *v = (*v as f64 * (1.0 + a)) as f32);
            m.vs.iter_mut()
                .for_each(|v| *v = (*v as f64 * (1.0 + a)) as f32);
            m
        };
        let traces = |m: Model| -> Vec<f64> {
            let mut sim = Simulation::with_model(&cfg, m).unwrap();
            sim.run(None, false).unwrap();
            sim.receivers.traces[0].iter().map(|&v| v as f64).collect()
        };
        let (s, sa, saa) = differences(&|a| traces(scaled(a)), 1e-2);
        let (gn_solver, second_solver) = hessian(&s, &sa, &saa, &observed);
        // The same with the solver's own data for the faster, denser solid: its error left
        // out of the residual.
        let own = {
            let mut m = scaled(FASTER);
            m.rho
                .iter_mut()
                .for_each(|v| *v = (*v as f64 * (1.0 + DENSER)) as f32);
            traces(m)
        };
        let (_, second_own) = hessian(&s, &sa, &saa, &own);

        // Born runs and Hessian products for the two directions.
        let n = model.vp.len();
        let to_f64 = |v: &[f32]| v.iter().map(|&x| x as f64).collect::<Vec<f64>>();
        let speeds_up = [to_f64(&model.vp), to_f64(&model.vs), vec![0.0; n]];
        let denser = [vec![0.0; n], vec![0.0; n], to_f64(&model.rho)];
        let directions = [speeds_up, denser];
        let mut jacobian = Vec::new();
        let mut products = Vec::new();
        for speeds in &directions {
            let speeds = speeds.each_ref().map(Vec::as_slice);
            let change = Change::from_velocities(&model, speeds);
            let mut sim = Simulation::with_model(&cfg, model.clone()).unwrap();
            jacobian.push(to_f64(&born(&mut sim, &change).unwrap().traces[0]));
            let mut sim = Simulation::with_model(&cfg, model.clone()).unwrap();
            let product = hessian_product(
                &mut sim,
                &change,
                |r| {
                    let res = r.traces[0]
                        .iter()
                        .zip(&observed32)
                        .map(|(&s, &o)| (w * (s as f64 - o as f64)) as f32)
                        .collect();
                    Ok((0.0, vec![res]))
                },
                |d| Ok(vec![d[0].iter().map(|&x| (w * x as f64) as f32).collect()]),
                None,
            )
            .unwrap();
            products.push(product.velocities(&model, speeds));
        }
        let mut full = [[0.0; 2]; 2];
        for (i, u) in directions.iter().enumerate() {
            for (j, p) in products.iter().enumerate() {
                full[i][j] = (0..3).map(|c| dot(&u[c], &p[c])).sum();
            }
        }
        let gn_born = [
            [
                w * dot(&jacobian[0], &jacobian[0]),
                w * dot(&jacobian[0], &jacobian[1]),
            ],
            [
                w * dot(&jacobian[1], &jacobian[0]),
                w * dot(&jacobian[1], &jacobian[1]),
            ],
        ];

        println!(
            "Lamb's problem, {label} ({:.0} s):",
            started.elapsed().as_secs_f64()
        );
        let samples = times.len();
        for (k, p) in survey::RECEIVERS.iter().enumerate() {
            let at = k * samples..(k + 1) * samples;
            let r = norm(sub(*p, survey::SOURCE));
            println!(
                "  {r:5.1} m: traces against Pekeris {:.2}%; ∂d/∂(speeds): solver's against \
                 Pekeris's {:.2}%, Born against the solver's {:.1e}",
                100.0 * misfit(&s[at.clone()], &d[at.clone()]),
                100.0 * misfit(&sa[at.clone()], &da[at.clone()]),
                misfit(&jacobian[0][at.clone()], &sa[at.clone()]),
            );
        }
        println!(
            "  ∂d/∂(density): Born against -d {:.1e}; the residual's size: Pekeris's {:.2}%, \
             the solver's against Pekeris's data {:.2}%, against its own {:.2}%",
            misfit(&jacobian[1], &s.iter().map(|v| -v).collect::<Vec<_>>()),
            100.0 * misfit(&d, &observed),
            100.0 * misfit(&s, &observed),
            100.0 * misfit(&s, &own),
        );
        let names = ["speeds, speeds", "speeds, density", "density, density"];
        for (name, (i, j)) in names.iter().zip([(0, 0), (0, 1), (1, 1)]) {
            let (exact, solver) = (
                gn_exact[i][j] + second_exact[i][j],
                gn_solver[i][j] + second_solver[i][j],
            );
            println!(
                "  H({name}): Pekeris {exact:.5e} (Gauss–Newton {:.5e}, second order {:+.5e}); \
                 the solver's {:+.2}% off it (GN {:+.2}%, second order {:+.2}%); the product \
                 {:.1e} from the solver's, Born's GN {:.1e}, the off-diagonal products agree to \
                 {:.1e}; with the solver's own data {:+.2}% (second order {:+.2}%)",
                gn_exact[i][j],
                second_exact[i][j],
                100.0 * rel(solver, exact),
                100.0 * rel(gn_solver[i][j], gn_exact[i][j]),
                100.0 * rel(second_solver[i][j], second_exact[i][j]),
                rel(full[i][j], solver),
                rel(gn_born[i][j], gn_solver[i][j]),
                rel(full[0][1], full[1][0]),
                100.0 * rel(gn_solver[i][j] + second_own[i][j], exact),
                100.0 * rel(second_own[i][j], second_exact[i][j]),
            );
        }
    }
}
