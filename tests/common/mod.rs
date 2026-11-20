//! What the whole-simulation tests share: setup helpers and the closed-form solutions for a
//! homogeneous full space (Aki & Richards 2002, §4.2–4.3).

#![allow(dead_code)]

use std::f64::consts::PI;
use std::path::Path;

use vibe_rations::config::Config;
use vibe_rations::receiver::Component;
use vibe_rations::simulation::Simulation;
use vibe_rations::source::Wavelet;

pub const VP: f64 = 3000.0;
pub const VS: f64 = 1732.0;
pub const RHO: f64 = 2500.0;

pub fn setup(toml: &str) -> Simulation {
    let cfg = Config::from_toml(toml).expect("a valid configuration");
    Simulation::new(&cfg, Path::new(".")).expect("a valid setup")
}

pub fn simulate(toml: &str) -> Simulation {
    let mut sim = setup(toml);
    sim.run(None, false).expect("a stable run");
    sim
}

/// A homogeneous model with the given grid, boundary and sources/receivers sections.
pub fn homogeneous(
    points: [usize; 3],
    pml: usize,
    free_surface: bool,
    duration: f64,
    rest: &str,
) -> String {
    format!(
        r#"
        [grid]
        points = {points:?}
        spacing = 10.0
        [time]
        duration = {duration}
        [boundary]
        pml_width = {pml}
        free_surface = {free_surface}
        [[model.layer]]
        top = 0.0
        vp = {VP}
        vs = {VS}
        rho = {RHO}
        {rest}
        "#
    )
}

/// ‖a - b‖ / ‖b‖ over a whole trace.
pub fn misfit(a: &[f64], b: &[f64]) -> f64 {
    let num: f64 = a.iter().zip(b).map(|(x, y)| (x - y).powi(2)).sum();
    let den: f64 = b.iter().map(|y| y * y).sum();
    (num / den).sqrt()
}

pub fn trace(sim: &Simulation, c: Component, r: usize) -> Vec<f64> {
    sim.receivers
        .trace(c, r)
        .expect("recorded")
        .iter()
        .map(|&v| v as f64)
        .collect()
}

pub fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub fn norm(a: [f64; 3]) -> f64 {
    a.iter().map(|x| x * x).sum::<f64>().sqrt()
}

/// Radial particle velocity at distance r from an explosion whose moment rate is m0·w(t):
/// v = m0/(4πρα²) [w(τ)/r² + w'(τ)/(α r)], τ = t - r/α.
pub fn explosion_velocity(w: &Wavelet, m0: f64, r: f64, t: f64) -> f64 {
    let tau = t - r / VP;
    m0 / (4.0 * PI * RHO * VP * VP) * (w.value(tau) / (r * r) + w.derivative(tau) / (VP * r))
}

/// Particle velocity at offset `d` from a point force f·w(t) along x (Stokes' solution,
/// differentiated in time).
pub fn force_velocity(w: &Wavelet, f: f64, d: [f64; 3], t: f64) -> [f64; 3] {
    let r = norm(d);
    let g = d.map(|x| x / r);
    // ∫ τ w'(t - τ) dτ from r/α to r/β, by the trapezium rule.
    let (a, b) = (r / VP, r / VS);
    let n = 4000;
    let h = (b - a) / n as f64;
    let near: f64 = (0..=n)
        .map(|m| {
            let tau = a + m as f64 * h;
            let weight = if m == 0 || m == n { 0.5 } else { 1.0 };
            weight * tau * w.derivative(t - tau)
        })
        .sum::<f64>()
        * h;
    std::array::from_fn(|i| {
        let delta = if i == 0 { 1.0 } else { 0.0 };
        let gg = g[i] * g[0];
        f / (4.0 * PI * RHO)
            * ((3.0 * gg - delta) / r.powi(3) * near
                + gg / (VP * VP * r) * w.derivative(t - r / VP)
                - (gg - delta) / (VS * VS * r) * w.derivative(t - r / VS))
    })
}
