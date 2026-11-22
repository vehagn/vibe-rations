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
