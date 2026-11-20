//! How accurate each stencil order and time step is over FWI-like distances: `just accuracy`.
//!
//! Not part of the regular suite (it takes a minute or two); it prints a table of misfits
//! against the closed-form solutions for every combination of stencil order and Courant
//! fraction, at coarse grids where numerical dispersion shows:
//!
//! - P: an explosion, 4 grid points per shortest P wavelength;
//! - S: a point force seen broadside (S in the far field), 3.5 and 4.6 points per shortest S
//!   wavelength.
//!
//! Receivers sit 300 to 1800 m from the source: 1.5–9 dominant P wavelengths, 1.7–10 S
//! wavelengths. The README's "Accuracy" section summarises the results.

mod common;

use std::time::Instant;

use common::*;
use vibe_rations::receiver::Component;
use vibe_rations::simulation::Simulation;
use vibe_rations::source::Wavelet;

const DISTANCES: [f64; 4] = [300.0, 900.0, 1500.0, 1800.0];
const ORDERS: [usize; 4] = [4, 8, 12, 16];
const COURANTS: [f64; 4] = [0.85, 0.5, 0.4, 0.3];

enum Wave {
    P,
    S,
}

/// A line of receivers along x from a source, in a domain long in x and narrow across.
fn config(wave: &Wave, h: f64, order: usize, courant: f64) -> String {
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
    let mut sim: Simulation = setup(&config(wave, h, order, courant));
    let start = Instant::now();
    sim.run(None, false).expect("a stable run");
    let seconds = start.elapsed().as_secs_f64();
    let speed = sim.solver.grid.interior_len() as f64 * sim.steps as f64 / seconds / 1e6;
    let times = sim.sample_times();
    let misfits = DISTANCES
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
            misfit(&trace(&sim, c, r), &exact)
        })
        .collect();
    (misfits, speed)
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
