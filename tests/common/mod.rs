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

/// How far a rebuilt wavefield drifts from the run it retraces: runs `sim` forward, saving
/// its band before every step and its whole state every `every` steps, then steps back to the
/// start. For each saved state, (its step, the largest difference at the nodes outside the
/// layers in any of the nine arrays, relative to that array's largest value there over the
/// run, and the state's own largest value there on the same scale).
pub fn rebuild_drift(sim: Simulation, every: usize) -> Vec<(usize, f64, f64)> {
    rebuild_drift_where(sim, |n| n % every == 0)
}

/// [`rebuild_drift`], saving the states before the steps `keep` picks.
pub fn rebuild_drift_where(
    mut sim: Simulation,
    keep: impl Fn(usize) -> bool,
) -> Vec<(usize, f64, f64)> {
    use vibe_rations::solver::Fields;
    let band = sim.band();
    let len = band.len();
    let steps = sim.steps;
    let [x, y, z] = sim.solver.interior();
    let g = sim.solver.grid;
    let nodes: Vec<usize> = z
        .flat_map(|k| {
            let x = x.clone();
            y.clone()
                .flat_map(move |j| x.clone().map(move |i| g.idx(i, j, k)))
        })
        .collect();
    let arrays = |f: &Fields| -> [Vec<f32>; 9] {
        [
            &f.vx, &f.vy, &f.vz, &f.sxx, &f.syy, &f.szz, &f.sxy, &f.sxz, &f.syz,
        ]
        .map(|a| nodes.iter().map(|&p| a[p]).collect())
    };
    let mut records = vec![0.0f32; steps * len];
    let mut saved: Vec<(usize, Fields)> = Vec::new();
    let mut scale = [0.0f32; 9];
    for n in 0..=steps {
        for (s, a) in scale.iter_mut().zip(arrays(&sim.solver.fields)) {
            *s = a.iter().fold(*s, |m, v| m.max(v.abs()));
        }
        if n == steps {
            break;
        }
        if keep(n) {
            saved.push((n, sim.solver.fields.clone()));
        }
        band.save(&sim.solver.fields, &mut records[n * len..(n + 1) * len]);
        sim.step();
    }
    let mut drift = Vec::with_capacity(saved.len());
    for n in (0..steps).rev() {
        sim.step_back(&band, &records[n * len..(n + 1) * len]);
        if let Some((_, f)) = saved.iter().find(|(m, _)| *m == n) {
            let (now, then) = (arrays(&sim.solver.fields), arrays(f));
            let worst = (0..9)
                .filter(|&c| scale[c] > 0.0)
                .map(|c| {
                    now[c]
                        .iter()
                        .zip(&then[c])
                        .fold(0.0f64, |m, (a, b)| m.max(((a - b).abs() / scale[c]) as f64))
                })
                .fold(0.0, f64::max);
            let level = (0..9)
                .filter(|&c| scale[c] > 0.0)
                .map(|c| {
                    then[c]
                        .iter()
                        .fold(0.0f64, |m, v| m.max((v.abs() / scale[c]) as f64))
                })
                .fold(0.0, f64::max);
            drift.push((n, worst, level));
        }
    }
    drift.reverse();
    drift
}
