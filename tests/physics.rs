//! Whole simulations checked against closed-form solutions and against each other.
//!
//! The homogeneous full-space solutions are from Aki & Richards (2002), §4.2–4.3. Each test
//! keeps receivers clear of the absorbing layers and stops before their (small) reflections
//! could matter, unless reflections are what it measures.

mod common;

use common::*;
use vibe_rations::receiver::Component;
use vibe_rations::source::Wavelet;

/// At the default: 0.005–0.007%; with Taylor's coefficients 0.002–0.006%, as they are near
/// exact for waves this long (8 points per shortest P wavelength). Opting out of the
/// time-dispersion transform, at C = 0.3: 0.07–0.19%.
#[test]
fn an_explosion_matches_the_full_space_solution() {
    let source = [400.0, 400.0, 400.0];
    let receivers = [
        [550.0, 400.0, 400.0],
        [400.0, 400.0, 650.0],
        [500.0, 500.0, 500.0],
        [400.0, 300.0, 400.0],
    ];
    let points: String = receivers
        .iter()
        .map(|p| format!("[[receivers.point]]\nposition = {p:?}\n"))
        .collect();
    let sim = simulate(&homogeneous(
        [81, 81, 81],
        12,
        false,
        0.3,
        &format!(
            r#"
            [[source]]
            position = {source:?}
            explosion = 1e13
            wavelet = {{ kind = "ricker", frequency = 15.0 }}
            [receivers]
            components = ["vx", "vy", "vz"]
            {points}
            "#
        ),
    ));
    let w = Wavelet::Ricker {
        f0: 15.0,
        t0: 1.2 / 15.0,
    };
    let times = sim.sample_times();
    for (r, &p) in receivers.iter().enumerate() {
        let d = sub(p, source);
        let dist = norm(d);
        let [vx, vy, vz] = [Component::Vx, Component::Vy, Component::Vz].map(|c| trace(&sim, c, r));
        let radial: Vec<f64> = (0..times.len())
            .map(|m| (vx[m] * d[0] + vy[m] * d[1] + vz[m] * d[2]) / dist)
            .collect();
        let exact: Vec<f64> = times
            .iter()
            .map(|&t| explosion_velocity(&w, 1e13, dist, t))
            .collect();
        let e = misfit(&radial, &exact);
        println!("explosion, r = {dist:.1} m: misfit {:.3}%", 100.0 * e);
        assert!(e < 2e-4, "misfit {e} at {p:?}");
    }
}

/// At the default, with leapfrog's time error removed and C = 0.85: 0.001–0.005%.
#[test]
fn a_point_force_matches_stokes_solution() {
    for (what, e) in point_force_misfits("") {
        assert!(e < 2e-4, "misfit {e} in {what}");
    }
}

/// With Taylor's coefficients instead of optimised ones: 0.001–0.008%.
#[test]
fn a_point_force_matches_stokes_solution_with_taylor_coefficients() {
    for (what, e) in point_force_misfits("[scheme]\ncoefficients = \"taylor\"") {
        assert!(e < 2e-4, "misfit {e} in {what}");
    }
}

/// Opting out of the time-dispersion transform (which then steps at C = 0.3): 0.03–0.09%. On
/// this grid (6.9 points per S wavelength) the time error dominates even at 0.3.
#[test]
fn a_point_force_matches_stokes_solution_without_the_transform() {
    for (what, e) in point_force_misfits("[scheme]\nremove_time_dispersion = false") {
        assert!(e < 2e-3, "misfit {e} in {what}");
    }
}

/// The misfits of a point force against Stokes' solution along, across and diagonal to it,
/// with `scheme` added to the configuration.
fn point_force_misfits(scheme: &str) -> Vec<(String, f64)> {
    let source = [400.0, 400.0, 400.0];
    // Along the force (P and near field), across it (S and near field), and in between.
    let receivers = [
        [600.0, 400.0, 400.0],
        [400.0, 600.0, 400.0],
        [540.0, 540.0, 400.0],
    ];
    let points: String = receivers
        .iter()
        .map(|p| format!("[[receivers.point]]\nposition = {p:?}\n"))
        .collect();
    let sim = simulate(&homogeneous(
        [81, 81, 81],
        12,
        false,
        0.4,
        &format!(
            r#"
            [[source]]
            position = {source:?}
            force = [1e9, 0.0, 0.0]
            wavelet = {{ kind = "ricker", frequency = 10.0 }}
            [receivers]
            components = ["vx", "vy", "vz"]
            {points}
            {scheme}
            "#
        ),
    ));
    let w = Wavelet::Ricker { f0: 10.0, t0: 0.12 };
    // Forces enter with velocity, at whole steps; the velocity samples are half a step later.
    let times = sim.sample_times();
    let mut out = Vec::new();
    for (r, &p) in receivers.iter().enumerate() {
        let d = sub(p, source);
        let exact: Vec<[f64; 3]> = times
            .iter()
            .map(|&t| force_velocity(&w, 1e9, d, t))
            .collect();
        for (i, c) in [Component::Vx, Component::Vy].into_iter().enumerate() {
            let want: Vec<f64> = exact.iter().map(|v| v[i]).collect();
            if want.iter().all(|v| v.abs() < 1e-12) {
                continue;
            }
            let e = misfit(&trace(&sim, c, r), &want);
            let what = format!("{} at offset {d:?}", c.name());
            println!(
                "force{}, {what}: misfit {:.3}%",
                match scheme {
                    "" => "",
                    s if s.contains("taylor") => ", Taylor's coefficients",
                    _ => ", without the transform",
                },
                100.0 * e
            );
            out.push((what, e));
        }
    }
    out
}

/// At the epicentre a P wave meets the free surface head-on, and the surface moves twice as
/// much as the same wave would in a full space.
///
/// The misfit is 2.52% at the default (2.53% with Taylor's coefficients, 2.40% opting out of the time-dispersion transform, at
/// C = 0.3, where the time error partly cancels the free surface's). Most of it is the
/// reference: doubling is exact only for a plane wave. On finer grids the misfit tends to
/// about 2.2% (`free_surface_on_finer_grids` in `just accuracy`), so the scheme's own share at
/// this spacing is about 0.36% (0.23% opting out).
#[test]
fn the_free_surface_doubles_a_vertically_incident_p_wave() {
    let depth = 300.0;
    let sim = simulate(&homogeneous(
        [81, 81, 61],
        12,
        true,
        0.3,
        &format!(
            r#"
            [[source]]
            position = [400.0, 400.0, {depth}]
            explosion = 1e13
            wavelet = {{ kind = "ricker", frequency = 15.0 }}
            [receivers]
            components = ["vz"]
            [[receivers.point]]
            position = [400.0, 400.0, 0.0]
            "#
        ),
    ));
    let w = Wavelet::Ricker {
        f0: 15.0,
        t0: 1.2 / 15.0,
    };
    let exact: Vec<f64> = sim
        .sample_times()
        .iter()
        // Up is -z.
        .map(|&t| -2.0 * explosion_velocity(&w, 1e13, depth, t))
        .collect();
    let e = misfit(&trace(&sim, Component::Vz, 0), &exact);
    println!("free surface: misfit {:.2}%", 100.0 * e);
    assert!(e < 0.026, "misfit {e}");
}

/// The same receivers, the same distance from the source, in a domain small enough that the
/// absorbing layers are close and in one large enough that they are out of reach. With a free
/// surface, the source stays at the same depth and only the sides and bottom move away.
fn layer_reflection(free_surface: bool) -> f64 {
    // The source and receivers sit 300 m deep with a free surface; without one, mid-domain.
    let run = |n: usize, nz: usize, centre: f64| {
        let c = centre;
        let z = if free_surface { 300.0 } else { c };
        simulate(&homogeneous(
            [n, n, nz],
            10,
            free_surface,
            0.3,
            &format!(
                r#"
                [[source]]
                position = [{c}, {c}, {z}]
                explosion = 1e13
                wavelet = {{ kind = "ricker", frequency = 15.0 }}
                [receivers]
                components = ["vx", "vz"]
                [[receivers.point]]
                position = [{}, {c}, {z}]
                [[receivers.point]]
                position = [{}, {}, {}]
                "#,
                c + 120.0,
                c + 80.0,
                c + 80.0,
                z + 80.0,
            ),
        ))
    };
    let small = run(61, 61, 300.0);
    let large = run(121, if free_surface { 101 } else { 121 }, 600.0);
    let peak = |t: &[f64]| t.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    let mut worst_overall = 0.0f64;
    for r in 0..2 {
        let components = [Component::Vx, Component::Vz];
        // Measured against the strongest component: on the x axis vz is zero by symmetry.
        let strongest = components
            .iter()
            .map(|&c| peak(&trace(&large, c, r)))
            .fold(0.0, f64::max);
        for c in components {
            let (s, l) = (trace(&small, c, r), trace(&large, c, r));
            let worst = s
                .iter()
                .zip(&l)
                .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()))
                / strongest;
            println!(
                "absorbing layers (free surface {free_surface}), receiver {r} {}: worst difference {:.4}% of peak",
                c.name(),
                100.0 * worst
            );
            worst_overall = worst_overall.max(worst);
        }
    }
    worst_overall
}

/// Measured: 1.7e-4 of the peak (about -75 dB) at the default, order 8 at C = 0.85 with the
/// time-dispersion transform; 3.9e-4 at order 8 and C = 0.3 without it, 1.5e-4 at order 4 and
/// C = 0.85 without it. Longer stencils reach further across the layers' inner edge.
#[test]
fn absorbing_layers_reflect_little() {
    let worst = layer_reflection(false);
    assert!(worst < 2e-4, "{worst}");
}

/// Along a free surface the side layers also damp plainly, to stay stable, which costs some
/// reflection. Measured: 1.9e-3 at the default; 2.1e-3 at order 8 and C = 0.3 without the
/// transform, 1.8e-3 at order 4.
#[test]
fn absorbing_layers_reflect_little_along_a_free_surface() {
    let worst = layer_reflection(true);
    assert!(worst < 2.2e-3, "{worst}");
}

/// Swapping a point force and a receiver gives the same trace, with the components swapped
/// too: G_ij(B, A) = G_ji(A, B). This holds for the discrete scheme only if the sources are
/// the transpose of the receivers and the free surface keeps the scheme self-adjoint, which
/// is what back-propagating residuals in FWI relies on. A layered model with a sphere in it
/// and a receiver on the free surface make it a hard case.
#[test]
fn swapping_source_and_receiver_gives_the_same_trace() {
    // Longer stencils reach further above the free surface; the mirroring must stay
    // self-adjoint for every order.
    for order in [4, 8] {
        reciprocity(order);
    }
}

fn reciprocity(order: usize) {
    let model = format!(
        r#"
        [grid]
        points = [60, 50, 40]
        spacing = 10.0
        [time]
        duration = 0.5
        [boundary]
        pml_width = 10
        [scheme]
        order = {order}
        [[model.layer]]
        top = 0.0
        vp = 2200.0
        vs = 1100.0
        rho = 2000.0
        [[model.layer]]
        top = 105.0
        vp = 3500.0
        vs = 2000.0
        rho = 2400.0
        [[model.sphere]]
        center = [300.0, 250.0, 150.0]
        radius = 40.0
        vp = 1800.0
        vs = 900.0
        rho = 1900.0
    "#
    );
    let run = |at: [f64; 3], force: [f64; 3], to: [f64; 3]| {
        simulate(&format!(
            r#"{model}
            [[source]]
            position = {at:?}
            force = {force:?}
            wavelet = {{ kind = "ricker", frequency = 12.0 }}
            [receivers]
            components = ["vx", "vy", "vz"]
            [[receivers.point]]
            position = {to:?}
            "#
        ))
    };
    let a = [180.0, 220.0, 0.0]; // on the free surface
    let b = [420.0, 260.0, 230.0];
    let unit =
        |axis: usize| std::array::from_fn::<f64, 3, _>(|i| if i == axis { 1e9 } else { 0.0 });
    let components = [Component::Vx, Component::Vy, Component::Vz];
    for (i, j) in [(2, 0), (0, 1), (2, 2)] {
        let ab = trace(&run(a, unit(i), b), components[j], 0);
        let ba = trace(&run(b, unit(j), a), components[i], 0);
        let peak = ab.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        let worst = ab
            .iter()
            .zip(&ba)
            .fold(0.0f64, |m, (x, y)| m.max((x - y).abs()));
        println!(
            "reciprocity, order {order}, force {} at A to {} at B: worst difference {:.2e} of peak",
            ["x", "y", "z"][i],
            components[j].name(),
            worst / peak
        );
        assert!(
            peak > 0.0 && worst < 1e-4 * peak,
            "{i}->{j}: {worst} vs peak {peak}"
        );
    }
}

/// A long run with a free surface, a water layer, soft sediments and the absorbing layers
/// meeting the surface at the corners: the energy must leave and stay gone.
#[test]
fn a_long_run_stays_stable_and_drains() {
    for order in [4, 8] {
        long_run(order);
    }
}

fn long_run(order: usize) {
    let mut sim = setup(&format!(
        r#"
        [grid]
        points = [48, 48, 40]
        spacing = 10.0
        [time]
        duration = 0.4
        [boundary]
        pml_width = 10
        [scheme]
        order = {order}
        [[model.layer]]
        top = 0.0
        vp = 1500.0
        vs = 0.0
        rho = 1000.0
        [[model.layer]]
        top = 20.0
        vp = 2000.0
        vs = 600.0
        rho = 1900.0
        [[model.layer]]
        top = 100.0
        vp = 4000.0
        vs = 2300.0
        rho = 2600.0
        [[source]]
        position = [240.0, 240.0, 200.0]
        moment = {{ xx = 1e13, yy = -5e12, xy = 3e12, xz = 4e12 }}
        wavelet = {{ kind = "ricker", frequency = 6.0 }}
        "#
    ));
    let mut peak = 0.0f32;
    for _ in 0..sim.steps {
        sim.step();
        peak = peak.max(sim.solver.peak_velocity());
    }
    assert!(peak > 0.0 && peak.is_finite());
    // Nine seconds, whatever the time step.
    let steps = (9.0 / sim.solver.dt) as usize;
    for n in 0..steps {
        sim.step();
        if n % 500 == 0 {
            let now = sim.solver.peak_velocity();
            assert!(
                now.is_finite() && now <= peak,
                "step {n}: {now} after a peak of {peak}"
            );
        }
    }
    let end = sim.solver.peak_velocity();
    println!(
        "long run, order {order}: peak {peak:.3e}, after {} s {end:.3e}",
        sim.solver.step as f64 * sim.solver.dt
    );
    assert!(end < 1e-4 * peak, "{end} left of a peak of {peak}");
}
