//! How far the forward wavefield rebuilt backwards from the band around the layers drifts
//! from the run it retraces, in the field and in the gradient: `just drift`.
//!
//! Not part of the regular suite: a study for the research log ("Rebuilding the forward
//! wavefield backwards"). Each table gives the largest difference at the nodes outside the
//! layers, in any of the nine arrays, relative to that array's largest value there over the
//! run, after stepping back the given number of steps from the end.

mod common;

use std::time::Instant;

use common::*;
use vibe_rations::gradient::{Strategy, gradient_with, l2_misfit};

/// A two-layer model with a sphere, a force and an explosion, recorded by a line of
/// receivers; `boundary` goes into [boundary], `top` is the model's first layer.
fn model(order: usize, courant: f64, duration: f64, boundary: &str, top: &str) -> String {
    format!(
        r#"
        [grid]
        points = [48, 40, 36]
        spacing = 10.0
        [time]
        duration = {duration}
        courant = {courant}
        [boundary]
        {boundary}
        [scheme]
        order = {order}
        optimised_for = 4.0
        [output]
        record_every = 2
        {top}
        [[model.layer]]
        top = 60.0
        vp = 2200.0
        vs = 1100.0
        rho = 2000.0
        [[model.layer]]
        top = 155.0
        vp = 3000.0
        vs = 1700.0
        rho = 2300.0
        [[model.sphere]]
        center = [240.0, 200.0, 120.0]
        radius = 35.0
        vp = 2400.0
        vs = 1200.0
        rho = 2100.0
        [[source]]
        position = [150.0, 160.0, 90.0]
        force = [3e8, 0.0, 1e9]
        wavelet = {{ kind = "ricker", frequency = 14.0 }}
        [[source]]
        position = [330.0, 230.0, 210.0]
        explosion = 2e12
        wavelet = {{ kind = "ricker", frequency = 14.0, delay = 0.1 }}
        [receivers]
        components = ["vx", "vy", "vz", "p"]
        [[receivers.line]]
        from = [120.0, 200.0, 100.0]
        to = [360.0, 200.0, 100.0]
        count = 8
        "#
    )
}

const SOLID_TOP: &str = "[[model.layer]]\ntop = 0.0\nvp = 2200.0\nvs = 1100.0\nrho = 2000.0";
const WATER_TOP: &str = "[[model.layer]]\ntop = 0.0\nvp = 1500.0\nvs = 0.0\nrho = 1000.0";
const CPML: &str = "pml_width = 8";

/// The drift after stepping back each of `back` steps from the end (the closest saved state).
fn drift_at(drift: &[(usize, f64, f64)], steps: usize, back: &[usize]) -> Vec<String> {
    back.iter()
        .map(|&b| {
            if b > steps {
                return format!("{:>9}", "");
            }
            let target = steps - b;
            let (_, e, _) = drift
                .iter()
                .min_by_key(|(n, _, _)| n.abs_diff(target))
                .expect("saved states");
            format!("{e:9.1e}")
        })
        .collect()
}

/// The drift as the rebuilt field steps back further and further: order 8 at C = 0.85 along
/// a free surface, about 3000 steps, with one step back from the exact final state for scale.
#[test]
#[ignore = "a study, not a check: run with `just drift`"]
fn drift_with_run_length() {
    let text = model(8, 0.85, 3.6, CPML, SOLID_TOP);
    let steps = setup(&text).steps;
    let drift = rebuild_drift(setup(&text), 10);
    let back = [1, 10, 30, 100, 300, 1000, 2000, steps];
    // One step back from the exact state at the end: one step's round-off.
    let one = rebuild_drift_where(setup(&text), |n| n + 1 == steps)[0].1;
    println!("\nrun length: order 8, C 0.85, free surface, {steps} steps");
    println!("one step back from the exact state: {one:.1e}");
    println!("steps back {}", back.map(|b| format!("{b:>9}")).join(""));
    println!("drift      {}", drift_at(&drift, steps, &back).join(""));
    // Where it grows: near the start, while the field is loud.
    println!("\nstep         drift    field's own level");
    for n in [1000, 500, 300, 200, 150, 100, 50, 20, 0] {
        if let Some((_, e, level)) = drift.iter().find(|(m, _, _)| *m == n) {
            println!("{n:>6}  {e:12.1e}  {level:12.1e}");
        }
    }
}

/// The drift after stepping all the way back over 1.2 s, by stencil order at C = 0.85 and by
/// Courant fraction at order 8.
#[test]
#[ignore = "a study, not a check: run with `just drift`"]
fn drift_with_order_and_time_step() {
    let back = [100, 300, 1000, 3000];
    println!("\norder and time step: 1.2 s, free surface");
    println!(
        "{:16} {:>6} {}",
        "",
        "steps",
        back.map(|b| format!("{b:>9}")).join("")
    );
    for (order, courant) in [
        (4, 0.85),
        (8, 0.85),
        (12, 0.85),
        (16, 0.85),
        (8, 0.5),
        (8, 0.3),
    ] {
        let text = model(order, courant, 1.2, CPML, SOLID_TOP);
        let steps = setup(&text).steps;
        let drift = rebuild_drift(setup(&text), 10);
        let worst = drift.iter().fold(0.0f64, |m, &(_, e, _)| m.max(e));
        println!(
            "{:16} {steps:>6} {}   all the way: {worst:.1e}",
            format!("order {order:2}, C {courant}"),
            drift_at(&drift, steps, &back).join("")
        );
    }
}

/// The drift with each kind of boundary, order 8 at C = 0.85 over 1.2 s: in exact arithmetic
/// the rebuilt field is exact whatever the boundary, so any difference here is how round-off
/// grows, not the boundary's own error.
#[test]
#[ignore = "a study, not a check: run with `just drift`"]
fn drift_with_boundaries() {
    let back = [100, 300, 1000];
    println!("\nboundaries: order 8, C 0.85, 1.2 s");
    for (label, boundary, top) in [
        ("layers, free surface", CPML, SOLID_TOP),
        (
            "layers all round",
            "pml_width = 8\nfree_surface = false",
            SOLID_TOP,
        ),
        (
            "rigid walls",
            "pml_width = 0\nfree_surface = false",
            SOLID_TOP,
        ),
        ("rigid walls, free surface", "pml_width = 0", SOLID_TOP),
        ("layers, water at the surface", CPML, WATER_TOP),
    ] {
        let text = model(8, 0.85, 1.2, boundary, top);
        let steps = setup(&text).steps;
        let drift = rebuild_drift(setup(&text), 10);
        let worst = drift.iter().fold(0.0f64, |m, &(_, e, _)| m.max(e));
        println!(
            "{label:30} {}   all the way ({steps}): {worst:.1e}",
            drift_at(&drift, steps, &back).join("")
        );
    }
}

/// The gradient from the rebuilt field against the replay's, at the nodes outside the layers,
/// for segments of increasing length up to the whole run, and the time and memory of each
/// against checkpoints.
#[test]
#[ignore = "a study, not a check: run with `just drift`"]
fn drift_in_the_gradient() {
    for duration in [1.2, 3.6] {
        let text = model(8, 0.85, duration, CPML, SOLID_TOP);
        let observed = simulate(&text.replace("vp = 2400.0", "vp = 2300.0"))
            .receivers
            .traces;
        let weights = vec![1.0; observed.len()];
        let run = |strategy| {
            let mut sim = setup(&text);
            let start = Instant::now();
            let r = gradient_with(
                &mut sim,
                |r| l2_misfit(&r.traces, &observed, &weights),
                strategy,
            )
            .unwrap();
            (
                sim.solver.interior(),
                sim.solver.grid,
                r,
                start.elapsed().as_secs_f64(),
            )
        };
        let (interior, g, reference, t0) = run(Strategy::default());
        let steps = setup(&text).steps;
        let nodes: Vec<usize> = interior[2]
            .clone()
            .flat_map(|k| {
                let x = interior[0].clone();
                interior[1]
                    .clone()
                    .flat_map(move |j| x.clone().map(move |i| g.interior_idx(i, j, k)))
            })
            .collect();
        println!("\ngradient: order 8, C 0.85, free surface, {steps} steps");
        println!("{:28} {:>10} {:>10} {:>8}", "", "off by", "memory", "time");
        println!(
            "{:28} {:>10} {:>10} {t0:7.1}s",
            format!("checkpoints every {}", reference.checkpoint_every),
            "",
            vibe_rations::memory::human(reference.memory_bytes as u64),
        );
        for segment in [Some(10), Some(100), None, Some(1000), Some(usize::MAX)] {
            let (_, _, r, t) = run(Strategy::Reconstruct { segment });
            let mut worst = 0.0f64;
            for (a, b) in [
                (&r.gradient.rho, &reference.gradient.rho),
                (&r.gradient.lambda, &reference.gradient.lambda),
                (&r.gradient.mu, &reference.gradient.mu),
            ] {
                let scale = nodes.iter().fold(0.0f64, |m, &n| m.max(b[n].abs()));
                worst = nodes
                    .iter()
                    .fold(worst, |m, &n| m.max((a[n] - b[n]).abs() / scale));
            }
            println!(
                "{:28} {worst:>10.1e} {:>10} {t:7.1}s",
                format!(
                    "rebuilt, segments of {}{}",
                    r.checkpoint_every,
                    if segment.is_none() { " (default)" } else { "" }
                ),
                vibe_rations::memory::human(r.memory_bytes as u64),
            );
        }
    }
}
