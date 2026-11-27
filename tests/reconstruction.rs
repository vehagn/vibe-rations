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

// The forward wavefield rebuilt inside a closed surface one node wide from the velocities and
// tractions recorded on it (Aaker et al. 2020; `vibe_rations::surface`), against the band's
// exact step back: the research log's "Rebuilding inside a surface" entry.

/// A larger model than [`model`], so that a surface up to eight nodes in from the layers
/// still holds the sources nine nodes clear of it: two layers, a sphere, a force and an
/// explosion at 14 Hz, a line of receivers. `boundary` goes into [boundary], `top` is the
/// model's first layer.
fn surface_model(order: usize, courant: f64, duration: f64, boundary: &str, top: &str) -> String {
    format!(
        r#"
        [grid]
        points = [64, 56, 48]
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
        top = 205.0
        vp = 3000.0
        vs = 1700.0
        rho = 2300.0
        [[model.sphere]]
        center = [320.0, 280.0, 160.0]
        radius = 35.0
        vp = 2400.0
        vs = 1200.0
        rho = 2100.0
        [[source]]
        position = [260.0, 260.0, 140.0]
        force = [3e8, 0.0, 1e9]
        wavelet = {{ kind = "ricker", frequency = 14.0 }}
        [[source]]
        position = [380.0, 300.0, 220.0]
        explosion = 2e12
        wavelet = {{ kind = "ricker", frequency = 14.0, delay = 0.1 }}
        [receivers]
        components = ["vx", "vy", "vz", "p"]
        [[receivers.line]]
        from = [200.0, 280.0, 100.0]
        to = [440.0, 280.0, 100.0]
        count = 8
        "#
    )
}

/// The band's drift and the surface's after stepping back each of `back` steps, by the
/// closest saved state.
fn surface_row(label: &str, text: &str, back: &[usize], every: usize, inset: Option<usize>) {
    let steps = setup(text).steps;
    let band = rebuild_drift(setup(text), every);
    let (surface, len) = surface_drift_where(setup(text), inset, 0, |n| n % every == 0);
    let worst = |d: &[(usize, f64, f64)]| d.iter().fold(0.0f64, |m, &(_, e, _)| m.max(e));
    println!(
        "{label:34} band    {}   all the way ({steps}): {:.1e}",
        drift_at(&band, steps, back).join(""),
        worst(&band)
    );
    println!(
        "{:34} surface {}   all the way: {:.1e}  ({len} values a step)",
        "",
        drift_at(&surface, steps, back).join(""),
        worst(&surface)
    );
}

/// How the surface's error grows as the rebuild steps back, against the band's: order 8 at
/// C = 0.85 along a free surface, about 3000 steps from the end in one go.
#[test]
#[ignore = "a study, not a check: run with `just drift`"]
fn surface_drift_with_run_length() {
    let text = surface_model(8, 0.85, 3.6, CPML, SOLID_TOP);
    let steps = setup(&text).steps;
    let back = [1, 10, 100, 300, 1000, 2000, steps];
    println!("\nsurface, run length: order 8, C 0.85, free surface, {steps} steps");
    println!(
        "{:42} {}",
        "steps back",
        back.map(|b| format!("{b:>9}")).join("")
    );
    surface_row("", &text, &back, 50, None);
    let (drift, _) = surface_drift_where(setup(&text), None, 0, |n| n % 50 == 0);
    println!("\nstep         surface's error    field's own level");
    for n in [3000, 2000, 1000, 500, 300, 200, 150, 100, 50, 0] {
        if let Some((_, e, level)) = drift.iter().find(|(m, _, _)| *m == n) {
            println!("{n:>6}  {e:16.1e}  {level:16.1e}");
        }
    }
    // Close to the face the error is larger: by how far in the nodes compared lie.
    println!("\nby the nodes' distance inside the surface (all the way back)");
    for margin in [1, 2, 4, 6, 9] {
        let (d, _) = surface_drift_where(setup(&text), None, margin, |n| n % 50 == 0);
        let worst = d.iter().fold(0.0f64, |m, &(_, e, _)| m.max(e));
        println!("  {margin} or more nodes in: {worst:.1e}");
    }
}

/// The surface's error by stencil order at C = 0.85 over 1.2 s, the surface L nodes in from
/// the layers, and at order 8 by how far in it sits.
#[test]
#[ignore = "a study, not a check: run with `just drift`"]
fn surface_drift_with_order_and_inset() {
    let back = [100, 300, 1000];
    println!("\nsurface, order and inset: 1.2 s, free surface");
    println!(
        "{:42} {}",
        "steps back",
        back.map(|b| format!("{b:>9}")).join("")
    );
    for order in [4, 8, 12, 16] {
        let text = surface_model(order, 0.85, 1.2, CPML, SOLID_TOP);
        surface_row(&format!("order {order:2}, inset L"), &text, &back, 10, None);
    }
    let text = surface_model(8, 0.85, 1.2, CPML, SOLID_TOP);
    for inset in [6, 8] {
        let (d, _) = surface_drift_where(setup(&text), Some(inset), 0, |n| n % 10 == 0);
        let steps = setup(&text).steps;
        println!(
            "{:34} surface {}   all the way: {:.1e}",
            format!("order  8, inset {inset}"),
            drift_at(&d, steps, &back).join(""),
            d.iter().fold(0.0f64, |m, &(_, e, _)| m.max(e))
        );
    }
}

/// The surface's error by how well the waves crossing it are sampled: the wavelets' peak
/// frequency lowered from 14 Hz (3.1 points per S wavelength at 2.5 f₀ in the 1100 m/s rock)
/// to 10 and 7 Hz (4.4 and 6.3), order 8 at C = 0.85 over 1.2 s along a free surface.
#[test]
#[ignore = "a study, not a check: run with `just drift`"]
fn surface_drift_with_frequency() {
    let back = [100, 300, 1000];
    println!("\nsurface, frequency: order 8, C 0.85, 1.2 s, free surface");
    println!(
        "{:42} {}",
        "steps back",
        back.map(|b| format!("{b:>9}")).join("")
    );
    for f in [14.0, 10.0, 7.0] {
        let text = surface_model(8, 0.85, 1.2, CPML, SOLID_TOP)
            .replace(
                "frequency = 14.0, delay = 0.1",
                &format!("frequency = {f:.1}, delay = {:.3}", 1.4 / f),
            )
            .replace("frequency = 14.0", &format!("frequency = {f:.1}"));
        let ppw = 1100.0 / (2.5 * f) / 10.0;
        surface_row(&format!("{f:4.1} Hz, {ppw:.1} ppw"), &text, &back, 10, None);
    }
}

/// The surface's error with each kind of boundary, order 8 at C = 0.85 over 1.2 s; and with
/// water above a seabed at, or near, the surface's top face (layers all round), which Aaker
/// et al. warn the interpolators can't take.
#[test]
#[ignore = "a study, not a check: run with `just drift`"]
fn surface_drift_with_boundaries() {
    let back = [100, 300, 1000];
    println!("\nsurface, boundaries: order 8, C 0.85, 1.2 s");
    println!(
        "{:42} {}",
        "steps back",
        back.map(|b| format!("{b:>9}")).join("")
    );
    // With layers all round the surface has a top face, 12 nodes down (120 m): the force
    // goes down to 260 m, 14 nodes clear of it.
    let all_round = "pml_width = 8\nfree_surface = false";
    let deeper = |text: String| text.replace("[260.0, 260.0, 140.0]", "[260.0, 260.0, 260.0]");
    for (label, text) in [
        (
            "layers, free surface",
            surface_model(8, 0.85, 1.2, CPML, SOLID_TOP),
        ),
        (
            "layers all round",
            deeper(surface_model(8, 0.85, 1.2, all_round, SOLID_TOP)),
        ),
        (
            "rigid walls, free surface",
            surface_model(8, 0.85, 1.2, "pml_width = 0", SOLID_TOP),
        ),
        (
            "layers, 60 m of water at the surface",
            surface_model(8, 0.85, 1.2, CPML, WATER_TOP),
        ),
    ] {
        surface_row(label, &text, &back, 10, None);
    }
    // Water above a seabed, layers all round: the seabed 6 nodes above the top face, on it,
    // 3 nodes below it and 8 below.
    for seabed in [60.0, 120.0, 150.0, 200.0] {
        let text = deeper(surface_model(8, 0.85, 1.2, all_round, WATER_TOP))
            .replace("top = 60.0", &format!("top = {seabed:.1}"));
        surface_row(
            &format!("all round, seabed at {seabed} m"),
            &text,
            &back,
            10,
            None,
        );
    }
}

/// The gradient rebuilt inside the surface and from the band against the replay's, at the
/// nodes where the surface's is meant to hold (L + 1 in from it), for segments of default
/// length and the whole run, with the time and memory of each: on the solid model over 1.2
/// and 3.6 s, and with 60 m of water at the top, inside the box (which is open there).
#[test]
#[ignore = "a study, not a check: run with `just drift`"]
fn surface_in_the_gradient() {
    use vibe_rations::surface::Surface;
    for (duration, top, label) in [
        (1.2, SOLID_TOP, "free surface"),
        (3.6, SOLID_TOP, "free surface"),
        (1.2, WATER_TOP, "60 m of water under the free surface"),
    ] {
        let text = surface_model(8, 0.85, duration, CPML, top);
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
            (r, start.elapsed().as_secs_f64())
        };
        let (reference, t0) = run(Strategy::default());
        let sim = setup(&text);
        let steps = sim.steps;
        let g = sim.solver.grid;
        let reach = sim.solver.stencil.half_length();
        let region = Surface::new(&sim.solver, None)
            .unwrap()
            .bounds
            .region(reach);
        let nodes: Vec<usize> = region[2]
            .clone()
            .flat_map(|k| {
                let x = region[0].clone();
                region[1]
                    .clone()
                    .flat_map(move |j| x.clone().map(move |i| g.interior_idx(i, j, k)))
            })
            .collect();
        let [vp0, vs0, rho0] = reference.gradient.velocities(&sim.solver.model);
        println!(
            "\ngradient: order 8, C 0.85, {label}, {steps} steps, compared at {} nodes",
            nodes.len()
        );
        println!(
            "{:36} {:>9} {:>9} {:>9} {:>10} {:>8}",
            "", "vp", "vs", "rho", "memory", "time"
        );
        println!(
            "{:36} {:>9} {:>9} {:>9} {:>10} {t0:7.1}s",
            format!("checkpoints every {}", reference.checkpoint_every),
            "",
            "",
            "",
            vibe_rations::memory::human(reference.memory_bytes as u64),
        );
        for (label, strategy) in [
            (
                "band, default segments",
                Strategy::Reconstruct { segment: None },
            ),
            (
                "band, the whole run",
                Strategy::Reconstruct {
                    segment: Some(usize::MAX),
                },
            ),
            (
                "surface, default segments",
                Strategy::Surface {
                    segment: None,
                    inset: None,
                },
            ),
            (
                "surface, the whole run",
                Strategy::Surface {
                    segment: Some(usize::MAX),
                    inset: None,
                },
            ),
        ] {
            let (r, t) = run(strategy);
            let [vp, vs, rho] = r.gradient.velocities(&sim.solver.model);
            let off = |a: &[f64], b: &[f64]| {
                let scale = nodes.iter().fold(0.0f64, |m, &n| m.max(b[n].abs()));
                nodes
                    .iter()
                    .fold(0.0f64, |m, &n| m.max((a[n] - b[n]).abs() / scale))
            };
            println!(
                "{:36} {:>9.1e} {:>9.1e} {:>9.1e} {:>10} {t:7.1}s",
                format!("{label} ({})", r.checkpoint_every),
                off(&vp, &vp0),
                off(&vs, &vs0),
                off(&rho, &rho0),
                vibe_rations::memory::human(r.memory_bytes as u64),
            );
        }
    }
}
