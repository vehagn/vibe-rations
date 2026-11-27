//! Whole simulations checked against closed-form solutions and against each other.
//!
//! The homogeneous full-space solutions are from Aki & Richards (2002), §4.2–4.3. Each test
//! keeps receivers clear of the absorbing layers and stops before their (small) reflections
//! could matter, unless reflections are what it measures.

mod common;

use common::*;
use vibe_rations::gradient::Strategy;
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

/// The same between nodes: the force at (40.37, 39.62, 40.89) nodes and the receivers the same
/// offsets from it, so that every source and receiver interpolates along all three axes, at
/// fractions of a node from 0.37 to 0.89. Measured: 0.001–0.005%, as on the nodes.
#[test]
fn a_point_force_between_nodes_matches_stokes_solution() {
    for (what, e) in point_force_misfits_at([403.7, 396.2, 408.9], "") {
        assert!(e < 2e-4, "misfit {e} in {what}");
    }
}

/// The misfits of a point force against Stokes' solution along, across and diagonal to it,
/// with `scheme` added to the configuration.
fn point_force_misfits(scheme: &str) -> Vec<(String, f64)> {
    point_force_misfits_at([400.0, 400.0, 400.0], scheme)
}

/// [`point_force_misfits`] with the force at `source`.
fn point_force_misfits_at(source: [f64; 3], scheme: &str) -> Vec<(String, f64)> {
    // Along the force (P and near field), across it (S and near field), and in between.
    let receivers = [[200.0, 0.0, 0.0], [0.0, 200.0, 0.0], [140.0, 140.0, 0.0]]
        .map(|d: [f64; 3]| std::array::from_fn::<f64, 3, _>(|a| source[a] + d[a]));
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
                    "" if source[0] != 400.0 => " between nodes",
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

/// Pekeris's closed form for Lamb's problem (`common::lamb`) is what it should be: zero at
/// the P arrival, where the surface starts upward; continuous at the S arrival; Boussinesq's
/// static displacement, F(1 - ν)/(2πμr), after the Rayleigh wave; and its velocity
/// converged in the number of intervals.
#[test]
fn lambs_reference_is_consistent() {
    use common::lamb::{gamma, survey};
    let lamb = survey::lamb();
    let (r, vs) = (240.0, survey::vs());
    let at = |tau: f64| lamb.step(r, tau * r / vs);
    let w = lamb.static_displacement(r);
    let mu = survey::RHO * vs * vs;
    assert!((w / (0.75 * 1e9 / (2.0 * std::f64::consts::PI * mu * r)) - 1.0).abs() < 1e-12);
    assert!(
        at(1.0 / 3f64.sqrt() + 1e-12).abs() < 1e-5 * w,
        "{}",
        at(0.5774)
    );
    assert!(at(0.6) < 0.0, "the surface starts upward");
    assert!((at(1.0 - 1e-12) - at(1.0 + 1e-12)).abs() < 1e-6 * w);
    assert_eq!(at(gamma() + 1e-12), w);
    let times: Vec<f64> = (0..420).map(|k| k as f64 * 1e-3).collect();
    let [coarse, fine] = [2000, 4000].map(|n| {
        times
            .iter()
            .map(|&t| lamb.velocity_with(r, t, n))
            .collect::<Vec<_>>()
    });
    let peak = fine.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    let worst = coarse
        .iter()
        .zip(&fine)
        .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()));
    println!(
        "lamb reference: 2000 against 4000 intervals, {:.1e} of the peak",
        worst / peak
    );
    assert!(worst < 1e-9 * peak, "{}", worst / peak);
}

/// Lamb's problem: a vertical force on the free surface of a Poisson solid, and the vertical
/// velocity on the surface against Pekeris's closed form (`common::lamb::survey`) at 150,
/// 240 and 330 m along x and 254.6 m along the diagonal. The Rayleigh wave dominates it: this
/// is the free surface's own accuracy, which the P-doubling test can't separate from its
/// reference.
///
/// Measured at the default, h = 10 m (10.6 points per Rayleigh wavelength at f0, 4.6 per
/// shortest S wavelength): 15.5%, 21.4%, 28.2% and 22.5%, growing with distance. The Rayleigh
/// wave runs about 1% fast (its peak 1.0 to 1.9 samples early), and its largest swing is 5%
/// too large. Order 16 gives the same, so the free surface, not the stencil, sets it; it falls
/// to 3.3–3.8% at 2.5 m (`lambs_problem_on_finer_grids` in `just accuracy`). The limits sit
/// just above the measurements, so a free surface that gets worse fails.
#[test]
fn a_force_on_the_free_surface_matches_lambs_problem() {
    use common::lamb::survey;
    let sim = simulate(&survey::config(10.0, "", ""));
    for ((r, e), limit) in survey::misfits(&sim)
        .into_iter()
        .zip([0.16, 0.22, 0.29, 0.23])
    {
        println!("lamb, vz at {r:.1} m: misfit {:.2}%", 100.0 * e);
        assert!(e < limit, "{r} m: {e}");
    }
}

/// With `[boundary] surface_slope`, sources and receivers near the free surface read the
/// evenly mirrored fields continued with their slope instead of across the mirror's kink
/// (`Solver::continue_with_slope`): sources stay the transpose of receivers, between nodes and
/// just below the surface too. Measured at order 8: reciprocity up to 1.1e-6 of the peak (7.3e-7
/// without: the slope's weights add single-precision round-off), the dot-product test 1.6e-8
/// on the nodes and 9.6e-8 between them.
#[test]
fn the_surface_slope_keeps_sources_the_transpose_of_receivers() {
    reciprocity(8, true);
    for shift in [[0.0; 3], [3.7, -2.4, 4.1]] {
        let e = dot_product(8, "rigid", true, "white", 1, shift, true);
        assert!(e < 1.5e-7, "shifted {shift:?}: {e}");
    }
}

/// Lamb's problem at 5 m with `[boundary] surface_slope`: the Rayleigh wave's largest swing,
/// 5.2–5.4% too large read across the mirror's kink, is within half a per cent, and the
/// misfits fall from 6.4, 7.8, 9.5 and (diagonal) 8.1% to 4.3, 6.2, 8.2 and 6.5%. What is
/// left is mostly the Rayleigh wave's speed, about 0.3% fast here (the research log's
/// 2026-11-23 entries).
#[test]
fn the_surface_slope_takes_the_kink_out_of_lambs_problem() {
    use common::lamb::survey;
    let sim = simulate(&survey::config(5.0, "", "").replace(
        "free_surface = true",
        "free_surface = true\nsurface_slope = true",
    ));
    for ([r, e, _, _, swing], limit) in survey::compare(&sim)
        .into_iter()
        .zip([0.044, 0.063, 0.084, 0.067])
    {
        println!(
            "lamb with the surface slope, h = 5 m, vz at {r:.1} m: misfit {:.2}%, largest swing {:+.2}%",
            100.0 * e,
            100.0 * (swing - 1.0)
        );
        assert!(e < limit, "{r} m: {e}");
        assert!((swing - 1.0).abs() < 0.005, "{r} m: swing {swing}");
    }
}

/// Plane P and S waves at normal incidence on a planar interface between two layers, against
/// the exact reflection and transmission (`common::plane`; R = -0.304 for P, -0.238 for S).
/// This checks the material averaging where a layer meets another: density at the
/// velocities, the moduli at the stresses, the shear modulus harmonically.
///
/// With the interface midway between two planes of nodes, where sampling the layers onto
/// the nodes is exact, at the default, at 20 (P) and 11.6 (S) points per wavelength at f0:
/// the direct waves 0.010% and 0.026%, reflections 0.72% and 7.5%, transmissions 0.078% and
/// 0.44%. These errors fall as h² and don't depend on the stencil's order (4 to 16): the
/// interface, not the stencil, sets them (`just interfaces`). With a layer's top on a plane
/// of nodes, which takes the layer below, the interface acts half a node higher: the
/// reflection comes h/vp early (35% and 57% misfit); the test pins that too, so a change to
/// how layers meet the nodes shows up here.
#[test]
fn a_planar_interface_reflects_and_transmits_plane_waves() {
    use common::plane::{PlaneWave, Wave};
    // (wave, direct, reflected, transmitted): misfits just above those measured.
    for (wave, limits) in [
        (Wave::P, [1.2e-4, 7.5e-3, 8.5e-4]),
        (Wave::S, [3e-4, 7.7e-2, 4.7e-3]),
    ] {
        let p = PlaneWave::interface(wave, 1, 0.5);
        let fits = p.fits(&p.run(None));
        for ((name, fit), limit) in ["direct", "reflected", "transmitted"]
            .iter()
            .zip(fits)
            .zip(limits)
        {
            println!(
                "interface, {wave:?} {name}: misfit {:.3}%, amplitude {:.4}, shift {:.3} h/v",
                100.0 * fit.misfit,
                fit.ratio,
                fit.shift / p.node_time()
            );
            assert!(fit.misfit < limit, "{wave:?} {name}: misfit {}", fit.misfit);
            assert!(
                fit.shift.abs() < 0.02 * p.node_time(),
                "{wave:?} {name}: {fit:?}"
            );
        }
        let p = PlaneWave::interface(wave, 1, 0.0);
        let [_, reflected, _] = p.fits(&p.run(None));
        let shift = reflected.shift / p.node_time();
        println!("interface on a node, {wave:?} reflected: shift {shift:.3} h/v");
        assert!((shift + 1.0).abs() < 0.02, "{wave:?}: shift {shift}");
    }
}

/// The same with `interfaces = "band-limited"`: the interface on a node and a quarter of the
/// way to the next act where they lie, with the same errors. Measured: P reflected 1.19% and
/// 1.17%, transmitted 0.125%; S 12.6% and 12.8%, 0.75%; shifts under 0.01 h/v. Node-sampled,
/// the same interfaces reflect 35% and 18% (P), 57% and 30% (S) off, h/v early and h/(2v)
/// late.
#[test]
fn band_limited_interfaces_act_where_they_lie() {
    use common::plane::{PlaneWave, Wave};
    use vibe_rations::config::Interfaces;
    for (wave, limits) in [
        (Wave::P, [1.2e-4, 1.25e-2, 1.35e-3]),
        (Wave::S, [3e-4, 0.13, 8e-3]),
    ] {
        for delta in [0.0, 0.25] {
            let mut p = PlaneWave::interface(wave, 1, delta);
            p.interfaces = Interfaces::BandLimited;
            let fits = p.fits(&p.run(None));
            for ((name, fit), limit) in ["direct", "reflected", "transmitted"]
                .iter()
                .zip(fits)
                .zip(limits)
            {
                println!(
                    "band-limited interface {delta} node below one, {wave:?} {name}: misfit \
                     {:.3}%, amplitude {:.4}, shift {:.3} h/v",
                    100.0 * fit.misfit,
                    fit.ratio,
                    fit.shift / p.node_time()
                );
                assert!(fit.misfit < limit, "{wave:?} {name}: misfit {}", fit.misfit);
                assert!(
                    fit.shift.abs() < 0.02 * p.node_time(),
                    "{wave:?} {name}: {fit:?}"
                );
            }
        }
    }
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
        reciprocity(order, false);
    }
}

fn reciprocity(order: usize, surface_slope: bool) {
    let model = format!(
        r#"
        [grid]
        points = [60, 50, 40]
        spacing = 10.0
        [time]
        duration = 0.5
        [boundary]
        pml_width = 10
        surface_slope = {surface_slope}
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
    let unit =
        |axis: usize| std::array::from_fn::<f64, 3, _>(|i| if i == axis { 1e9 } else { 0.0 });
    let components = [Component::Vx, Component::Vy, Component::Vz];
    // On the nodes, A on the free surface; between nodes, A on the surface and then just
    // below it, between the surface and the first node, where both interpolate across the
    // mirror images.
    for (label, a, b) in [
        ("", [180.0, 220.0, 0.0], [420.0, 260.0, 230.0]),
        (
            ", between nodes",
            [183.7, 216.1, 0.0],
            [416.3, 263.8, 227.4],
        ),
        (
            ", between nodes, A below the surface",
            [176.2, 223.9, 4.1],
            [416.3, 263.8, 227.4],
        ),
    ] {
        for (i, j) in [(2, 0), (0, 1), (2, 2)] {
            let ab = trace(&run(a, unit(i), b), components[j], 0);
            let ba = trace(&run(b, unit(j), a), components[i], 0);
            let peak = ab.iter().fold(0.0f64, |m, v| m.max(v.abs()));
            let worst = ab
                .iter()
                .zip(&ba)
                .fold(0.0f64, |m, (x, y)| m.max((x - y).abs()));
            let slope = if surface_slope { ", surface slope" } else { "" };
            println!(
                "reciprocity, order {order}{slope}{label}, force {} at A to {} at B: worst difference {:.2e} of peak",
                ["x", "y", "z"][i],
                components[j].name(),
                worst / peak
            );
            assert!(
                peak > 0.0 && worst < 1e-4 * peak,
                "{label} {i}->{j}: {worst} vs peak {peak}"
            );
        }
    }
}

/// A long run with a free surface, a water layer, soft sediments and the absorbing layers
/// meeting the surface at the corners: the energy must leave and stay gone.
///
/// Measured, what is left after 9.4 s relative to the peak, orders 4 / 8: 2.5e-6 / 9.6e-7 with
/// band-limited layers, 2.7e-6 / 1.1e-6 node-sampled (until 2026-11-22).
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

/// The adjoint run is the transpose of the forward run. For a force history w at A, recorded
/// as d = F w at receivers in every component, and residuals d̃ at those receivers,
/// ⟨F w, d̃⟩ = ⟨w, Fᵀ d̃⟩, where Fᵀ d̃ is the adjoint velocity read at A after each adjoint step.
/// The scheme is self-adjoint everywhere but in the absorbing layers, so this holds to
/// round-off for receivers outside them, whatever the waves do inside. "Outside" includes
/// their interpolation: at order 8 a velocity reads four staggered values either side, and a
/// receiver whose values reach into a layer (one 120 m from the side) gave 2.4e-7 here, and
/// 2e-7 to 1.2e-6 over other seeds; `dot_product` now refuses such receivers.
#[test]
fn the_adjoint_run_is_the_transpose_of_the_forward_run() {
    for order in [4, 8] {
        for (walls, free_surface, signal, every) in [
            ("rigid", true, "white", 1),
            ("cpml", true, "band-limited", 2),
            ("cpml", false, "white", 1),
        ] {
            let e = dot_product(order, walls, free_surface, signal, every, [0.0; 3], false);
            // Measured 2.3e-8 to 8.3e-8: single precision over about 1300 steps. Leaving the
            // adjoint's initial pressure state unmirrored gives 2.4e-7 to 2.8e-6.
            assert!(e < 1.5e-7, "order {order}, {walls}, {signal}: {e}");
        }
    }
}

/// The same with the source and every receiver between nodes, the shallow ones 4.1 m below
/// the free surface, between it and the first node: the operators that interpolate along all
/// three axes and across the mirror images, pressure's included, are transposes too.
#[test]
fn the_adjoint_run_is_the_transpose_of_the_forward_run_between_nodes() {
    for (walls, signal, every) in [("rigid", "white", 1), ("cpml", "band-limited", 2)] {
        let e = dot_product(8, walls, true, signal, every, [3.7, -2.4, 4.1], false);
        assert!(e < 1.5e-7, "{walls}, {signal}: {e}");
    }
}

/// Uniform values in [0, 1) from a fixed seed (xorshift64*).
fn noise(n: usize, seed: u64) -> Vec<f64> {
    let mut x = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
    (0..n)
        .map(|_| {
            x ^= x >> 12;
            x ^= x << 25;
            x ^= x >> 27;
            (x.wrapping_mul(0x2545_F491_4F6C_DD1D) >> 11) as f64 / (1u64 << 53) as f64
        })
        .collect()
}

/// Zero-mean noise smoothed by a Ricker wavelet of `f0` Hz sampled every `dt`, scaled to a
/// peak of 1.
fn band_limited(n: usize, seed: u64, f0: f64, dt: f64) -> Vec<f64> {
    let w = Wavelet::Ricker { f0, t0: 0.0 };
    let half = (1.5 / f0 / dt).ceil() as isize;
    let kernel: Vec<f64> = (-half..=half).map(|k| w.value(k as f64 * dt)).collect();
    let x: Vec<f64> = noise(n, seed).iter().map(|v| 2.0 * v - 1.0).collect();
    let y: Vec<f64> = (0..n as isize)
        .map(|i| {
            (-half..=half)
                .filter(|k| (0..n as isize).contains(&(i - k)))
                .map(|k| x[(i - k) as usize] * kernel[(k + half) as usize])
                .sum()
        })
        .collect();
    let peak = y.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    y.iter().map(|v| v / peak).collect()
}

/// |⟨F w, d̃⟩ - ⟨w, Fᵀ d̃⟩| / ⟨F w, d̃⟩ in the reciprocity test's model, with every component
/// recorded at five receivers, one on the free surface and one just below it if there is one. The residuals are the
/// data times random positive weights, so that the terms of ⟨F w, d̃⟩ can't cancel.
/// `shift` (m) moves the source and the receivers off the nodes.
fn dot_product(
    order: usize,
    walls: &str,
    free_surface: bool,
    signal: &str,
    every: usize,
    shift: [f64; 3],
    surface_slope: bool,
) -> f64 {
    let pml = if walls == "rigid" { 0 } else { 10 };
    // Without a free surface the top 100 m are an absorbing layer, and the shallow receivers
    // sit where their interpolation (four nodes up at order 8) stays out of it.
    let shallow = if free_surface { 0.0 } else { 130.0 };
    let receivers: String = [
        [180.0, 220.0, shallow],
        [420.0, 260.0, 230.0],
        [300.0, 150.0, shallow + 60.0],
        [250.0, 300.0, 140.0],
        // One node below the free surface, where the stresses' mirror images reach; 140 m
        // from the side, so that vy's interpolation (four nodes at order 8) stays out of the
        // layer.
        [400.0, 140.0, shallow + 10.0],
    ]
    .iter()
    .map(|p| {
        let p: [f64; 3] = std::array::from_fn(|a| p[a] + shift[a]);
        format!("[[receivers.point]]\nposition = {p:?}\n")
    })
    .collect();
    let source: [f64; 3] = std::array::from_fn(|a| [330.0, 240.0, 200.0][a] + shift[a]);
    let mut sim = setup(&format!(
        r#"
        [grid]
        points = [60, 50, 40]
        spacing = 10.0
        [time]
        duration = 0.5
        [boundary]
        pml_width = {pml}
        free_surface = {free_surface}
        surface_slope = {surface_slope}
        [scheme]
        order = {order}
        [output]
        record_every = {every}
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
        [[source]]
        position = {source:?}
        force = [0.6, 0.0, 0.8]
        wavelet = {{ kind = "ricker", frequency = 12.0 }}
        [receivers]
        components = ["vx", "vy", "vz", "p"]
        {receivers}
        "#
    ));
    // A receiver reading damped values would make the test measure the layers instead.
    assert_eq!(
        sim.receivers_next_to_layers(),
        0,
        "receivers next to the layers"
    );
    let (steps, dt) = (sim.steps, sim.solver.dt);
    let w: Vec<f64> = match signal {
        "white" => noise(steps, 1).iter().map(|v| 2.0 * v - 1.0).collect(),
        _ => band_limited(steps, 1, 12.0, dt),
    };
    // Forward: the force history w through the source's operator, recorded as configured.
    let op = sim.sources[0].op.clone();
    for (n, &value) in w.iter().enumerate() {
        sim.solver
            .step_with(|f, _| op.inject(f, value as f32), |_, _| {});
        if n.is_multiple_of(every) {
            sim.receivers.record(&sim.solver);
        }
    }
    let len = sim.receivers.points.len() * sim.receivers.samples;
    let residuals: Vec<Vec<f64>> = sim
        .receivers
        .traces
        .iter()
        .zip(10..)
        .map(|(d, seed)| {
            let weights = match signal {
                "white" => noise(len, seed),
                _ => band_limited(len, seed, 12.0 / every as f64, dt * every as f64)
                    .iter()
                    .map(|v| 1.0 + 0.5 * v)
                    .collect(),
            };
            // Pressure is ~10⁶ times velocity: each component gets the same say.
            let norm: f64 = d.iter().map(|&x| (x as f64).powi(2)).sum();
            assert!(norm > 0.0, "every component sees the source");
            d.iter()
                .zip(weights)
                .map(|(&d, u)| d as f64 * u / norm)
                .collect()
        })
        .collect();
    let forward: f64 = sim
        .receivers
        .traces
        .iter()
        .zip(&residuals)
        .flat_map(|(d, r)| d.iter().zip(r).map(|(&a, &b)| a as f64 * b))
        .sum();
    // Adjoint: the residuals back from the receivers, read at A along the force.
    let node = sim.sources[0].point;
    let mut adjoint = vibe_rations::adjoint::Adjoint::new(
        &sim,
        residuals
            .iter()
            .map(|r| r.iter().map(|&x| x as f32).collect())
            .collect(),
    )
    .expect("residuals that fit");
    let mut back = vec![0.0f64; steps];
    adjoint
        .run(|n, s| {
            back[n] = 0.6 * s.sample(Component::Vx, node) as f64
                + 0.8 * s.sample(Component::Vz, node) as f64;
        })
        .expect("a stable adjoint run");
    let backward: f64 = w.iter().zip(&back).map(|(a, b)| a * b).sum();
    let e = (forward - backward).abs() / forward;
    println!(
        "dot product, order {order}, {walls} walls, free surface {free_surface}, {signal}, \
         every {every}, shifted {shift:?}: ⟨Fw, d̃⟩ {forward:.6e}, ⟨w, Fᵀd̃⟩ {backward:.6e}, \
         difference {e:.1e}"
    );
    e
}

/// The gradient is that of the solver's forward run: perturbing vp, vs or ρ at every node
/// outside the absorbing layers by random amounts down the gradient changes the misfit as the
/// gradient says, to within the finite difference's own error, with the time step and the
/// layers' damping held at the model's. The model has a free surface, two band-limited layers
/// and a sphere; a force and an explosion; receivers on the surface and at depth recording
/// every component, every second step.
///
/// Measured, extrapolated (orders 4 / 8), at the defaults (optimised coefficients pinned,
/// C = 0.85, traces un-warped, band-limited layers): vp 8.5e-6 / 2.3e-6, vs 6.9e-5 / 7.4e-5,
/// ρ 2.5e-6 / 1.2e-5 (node-sampled layers: vp 1.0e-5 / 6.3e-6, vs 7.1e-5 / 7.2e-5, ρ 4.8e-6 /
/// 1.3e-6). The plain central differences are off by 1e-4 to 1.4e-1 at these steps, falling
/// 4× per halving: their own h² error. For vp the extrapolated error falls as h⁴ from 1e-4 at
/// h = 3.2e-2 to a floor of about 3e-6 from 8e-3 down. Until 2026-11-22 the check perturbed
/// random signs, and vp only below the fastest speed, which set dt and the layers: the whole
/// lower layer stayed fixed. Then (node-sampled): vp 2.8e-4 / 3.9e-6, vs 2.1e-5 / 1.0e-5,
/// ρ 4.4e-5 / 1.6e-5, order 4's vp near the limit because its change cancelled to 4.5e-5.
/// Held dt and layers but random signs left vp's change cancelling tenfold and the check on
/// round-off, 2.5e-5 to 5e-4 by step without a trend (it failed at 3.0e-4). Before the
/// receiver line moved in from 100 m (at order 8 the first receiver read into a layer), with
/// Taylor coefficients at C = 0.3: vp 1.7e-4 / 6.5e-5, vs 1.6e-6 / 5.8e-6, ρ 1.5e-5 / 1.2e-5.
#[test]
fn the_gradient_matches_finite_differences() {
    for order in [4, 8] {
        gradient_check(
            &format!("order {order}"),
            8,
            &|sphere| gradient_model(order, sphere),
            Strategy::default(),
        );
    }
}

/// The same between nodes: every source and receiver moved by (-3.7, -2.4, 4.1) m, so that the
/// surface line, which records pressure too, sits 4.1 m down, between the surface and the
/// first node. This covers what the gradient takes from residuals read between nodes,
/// pressure's included (`gradient::PressureResiduals`), which the dot-product test doesn't.
#[test]
fn the_gradient_matches_finite_differences_between_nodes() {
    gradient_check(
        "between nodes, order 8",
        8,
        &|sphere| shifted(&gradient_model(8, sphere), [-3.7, -2.4, 4.1]),
        Strategy::default(),
    );
}

/// `config` with every `position = [x, y, z]` moved by `shift` (m).
fn shifted(config: &str, shift: [f64; 3]) -> String {
    config
        .lines()
        .map(
            |line| match line.trim_start().strip_prefix("position = [") {
                Some(rest) => {
                    let indent = &line[..line.len() - line.trim_start().len()];
                    let p: Vec<f64> = rest
                        .trim_end_matches(']')
                        .split(',')
                        .map(|v| v.trim().parse().unwrap())
                        .collect();
                    let q: [f64; 3] = std::array::from_fn(|a| p[a] + shift[a]);
                    format!("{indent}position = {q:?}")
                }
                None => line.to_string(),
            },
        )
        .collect::<Vec<_>>()
        .join("\n")
}

/// The same where the first check doesn't reach: rigid walls, so every node is perturbed,
/// the last ones (whose missing neighbours the material averaging clamps) included; uneven
/// spacing; a water layer under the free surface, with fluid nodes on it and against the
/// solid; a double couple; and every third step recorded.
///
/// Measured, extrapolated, with Taylor coefficients at C = 0.85, band-limited layers, down the
/// gradient: vp 5.4e-6, vs 9.1e-6, ρ 5.1e-6. Along random signs, node-sampled, vp below the
/// fastest only (until 2026-11-22): vp 6.4e-5, vs 1.2e-5, ρ 1.6e-4 (at C = 0.3: vp 4.6e-5,
/// vs 8.9e-7, ρ 2.2e-4); ρ's plain differences were then 37% off at h = 1.6e-2, and with
/// optimised coefficients 61%, extrapolating to 6.0e-4, a finite-difference limit (perturbing
/// the water layer alone, ρ agreed to 2.7e-5). This check is about the geometry, and the first
/// one covers optimised coefficients.
#[test]
fn the_gradient_matches_finite_differences_with_walls_and_water() {
    gradient_check("walls and water", 0, &walls_and_water, Strategy::default());
}

/// The same checks with the forward wavefield rebuilt backwards from the band around the
/// layers instead of replayed (`Strategy::Reconstruct`, default segments).
///
/// Measured, extrapolated (replayed in brackets): order 4 vp 8.6e-6 (8.5e-6), vs 6.9e-5
/// (6.9e-5), ρ 2.6e-6 (2.5e-6); order 8 vp 2.3e-6, vs 7.4e-5, ρ 1.2e-5 (all as replayed);
/// walls and water vp 5.4e-6, vs 9.1e-6, ρ 5.1e-6 (as replayed). The rebuilt gradient moves the
/// predicted change by 4e-7 at most, well inside what the finite differences resolve. Along
/// random signs, node-sampled (until 2026-11-22): order 4 vp 2.6e-4 (2.8e-4), vs 2.1e-5, ρ
/// 4.4e-5; order 8 vp 1.4e-6 (3.9e-6), vs 1.0e-5, ρ 1.6e-5; walls and water vp 6.3e-5
/// (6.4e-5), vs 1.2e-5, ρ 1.6e-4.
#[test]
fn the_gradient_matches_finite_differences_when_rebuilt() {
    let rebuilt = Strategy::Reconstruct { segment: None };
    for order in [4, 8] {
        gradient_check(
            &format!("rebuilt, order {order}"),
            8,
            &|sphere| gradient_model(order, sphere),
            rebuilt,
        );
    }
    gradient_check("rebuilt, walls and water", 0, &walls_and_water, rebuilt);
}

/// Rebuilt backwards from the band around the layers, the forward wavefield retraces the run
/// at every node outside them: at orders 4 and 8, along a free surface and with layers all
/// round, and with water under the surface and rigid walls (no layers, so the band is the
/// halo, always zero). The error is single-precision round-off, which accumulates as the
/// rebuilt field steps back.
///
/// Measured, largest over the run relative to each array's peak, with band-limited layers:
/// order 4, 253 steps, 8.0e-7 with a free surface and 9.0e-7 with layers all round; order 8,
/// 285 steps, 8.0e-7 and 1.0e-6; walls and water, 216 steps, 1.8e-6. Node-sampled (until
/// 2026-11-22): order 4, 251 steps, 7.5e-7 and 6.3e-7; order 8, 283 steps, 6.2e-7 and 1.2e-6;
/// walls and water, 213 steps, 9.6e-7. The water's overshoot beside the rock stiffens a few
/// nodes, and the drift grows with the loudest stretch of the field. The tolerance, 5e-6, was
/// four times the largest node-sampled value.
#[test]
fn the_rebuilt_wavefield_follows_the_run() {
    for order in [4, 8] {
        let model = gradient_model(order, [2300.0, 1150.0, 2050.0]);
        for (label, text) in [
            ("free surface", model.clone()),
            (
                "layers all round",
                model.replace("pml_width = 8", "pml_width = 8\nfree_surface = false"),
            ),
        ] {
            let drift = rebuild_drift(setup(&text), 10);
            let worst = drift.iter().fold(0.0f64, |m, &(_, e, _)| m.max(e));
            println!(
                "rebuilt wavefield, order {order}, {label}: {} steps, at most {worst:.1e}",
                setup(&text).steps
            );
            assert!(worst < 5e-6, "order {order}, {label}: {worst}");
        }
    }
    let water = walls_and_water([2300.0, 1150.0, 2050.0]);
    let drift = rebuild_drift(setup(&water), 10);
    let worst = drift.iter().fold(0.0f64, |m, &(_, e, _)| m.max(e));
    println!(
        "rebuilt wavefield, walls and water: {} steps, at most {worst:.1e}",
        setup(&water).steps
    );
    assert!(worst < 5e-6, "walls and water: {worst}");
}

/// The gradient from the rebuilt wavefield against the one from the replay, at every node
/// outside the layers: the default segments and the whole run as one, on the models of the
/// finite-difference checks.
///
/// Measured, largest difference relative to each parameter's largest value (default segments,
/// then the whole run), with band-limited layers: order 4, 3.9e-7 (85 steps) and 5.0e-7 (253);
/// order 8, 3.6e-7 (57) and 5.4e-7 (285); walls and water, 3.3e-7 (31) and 6.2e-7 (216).
/// Node-sampled (until 2026-11-22): 2.9e-7 (84) and 5.3e-7 (251); 3.2e-7 (57) and 5.5e-7 (283);
/// 3.9e-7 (31) and 7.1e-7 (213). The tolerance is four times the largest.
#[test]
fn rebuilt_gradients_match_the_replay() {
    use vibe_rations::gradient::{gradient_with, l2_misfit};
    let order4 = |sphere| gradient_model(4, sphere);
    let order8 = |sphere| gradient_model(8, sphere);
    type Config<'a> = &'a dyn Fn([f64; 3]) -> String;
    let cases: [(&str, Config); 3] = [
        ("order 4", &order4),
        ("order 8", &order8),
        ("walls and water", &walls_and_water),
    ];
    for (label, config) in cases {
        let text = config([2050.0, 1000.0, 2000.0]);
        let observed = simulate(&config([1900.0, 950.0, 1950.0])).receivers.traces;
        let weights = vec![1.0; observed.len()];
        let run = |strategy| {
            let mut sim = setup(&text);
            let r = gradient_with(
                &mut sim,
                |r| l2_misfit(&r.traces, &observed, &weights),
                strategy,
            )
            .unwrap();
            (sim.solver.interior(), sim.solver.grid, r)
        };
        let (interior, g, reference) = run(Strategy::default());
        let nodes: Vec<usize> = interior[2]
            .clone()
            .flat_map(|k| {
                let x = interior[0].clone();
                interior[1]
                    .clone()
                    .flat_map(move |j| x.clone().map(move |i| g.interior_idx(i, j, k)))
            })
            .collect();
        for segment in [None, Some(usize::MAX)] {
            let (_, _, r) = run(Strategy::Reconstruct { segment });
            let mut worst = 0.0f64;
            for (a, b) in [
                (&r.gradient.rho, &reference.gradient.rho),
                (&r.gradient.lambda, &reference.gradient.lambda),
                (&r.gradient.mu, &reference.gradient.mu),
            ] {
                let scale = nodes.iter().fold(0.0f64, |m, &n| m.max(b[n].abs()));
                assert!(scale > 0.0, "{label}: no gradient");
                let e = nodes
                    .iter()
                    .fold(0.0f64, |m, &n| m.max((a[n] - b[n]).abs() / scale));
                worst = worst.max(e);
            }
            println!(
                "rebuilt gradient, {label}, segments of {}: off by {worst:.1e}",
                r.checkpoint_every
            );
            assert!(worst < 3e-6, "{label}, {segment:?}: {worst}");
        }
    }
}

/// Its sources and receivers sit on the nodes they snapped to before positions were exact.
fn walls_and_water(sphere: [f64; 3]) -> String {
    let [vp, vs, rho] = sphere;
    format!(
        r#"
        [grid]
        points = [34, 30, 26]
        spacing = [11.0, 12.0, 10.5]
        [time]
        duration = 0.3
        [boundary]
        pml_width = 0
        [scheme]
        coefficients = "taylor"
        [output]
        record_every = 3
        [[model.layer]]
        top = 0.0
        vp = 1500.0
        vs = 0.0
        rho = 1000.0
        [[model.layer]]
        top = 21.0
        vp = 2200.0
        vs = 1100.0
        rho = 2000.0
        [[model.layer]]
        top = 150.0
        vp = 3000.0
        vs = 1700.0
        rho = 2300.0
        [[model.sphere]]
        center = [180.0, 170.0, 100.0]
        radius = 30.0
        vp = {vp}
        vs = {vs}
        rho = {rho}
        [[source]]
        position = [121.0, 204.0, 157.5]
        double_couple = {{ strike = 30.0, dip = 60.0, rake = 90.0, m0 = 1e12 }}
        wavelet = {{ kind = "ricker", frequency = 10.0 }}
        [[source]]
        position = [253.0, 120.0, 84.0]
        force = [0.0, 5e8, 1e9]
        wavelet = {{ kind = "ricker", frequency = 10.0, delay = 0.15 }}
        [receivers]
        components = ["vx", "vy", "vz", "p"]
        # The nodes a line from 44 to 330 m snapped to before positions were exact.
        [[receivers.point]]
        position = [44.0, 180.0, 0.0]
        [[receivers.point]]
        position = [77.0, 180.0, 0.0]
        [[receivers.point]]
        position = [121.0, 180.0, 0.0]
        [[receivers.point]]
        position = [154.0, 180.0, 0.0]
        [[receivers.point]]
        position = [187.0, 180.0, 0.0]
        [[receivers.point]]
        position = [220.0, 180.0, 0.0]
        [[receivers.point]]
        position = [264.0, 180.0, 0.0]
        [[receivers.point]]
        position = [297.0, 180.0, 0.0]
        [[receivers.point]]
        position = [330.0, 180.0, 0.0]
        [[receivers.point]]
        position = [198.0, 96.0, 220.5]
        "#
    )
}

fn gradient_model(order: usize, sphere: [f64; 3]) -> String {
    let [vp, vs, rho] = sphere;
    format!(
        r#"
        [grid]
        points = [48, 40, 36]
        spacing = 10.0
        [time]
        duration = 0.35
        [boundary]
        pml_width = 8
        [scheme]
        order = {order}
        # Perturbed models must keep the operator: the band would follow their slowest wave.
        optimised_for = 4.0
        [output]
        record_every = 2
        [[model.layer]]
        top = 0.0
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
        vp = {vp}
        vs = {vs}
        rho = {rho}
        [[source]]
        position = [150.0, 160.0, 60.0]
        force = [3e8, 0.0, 1e9]
        wavelet = {{ kind = "ricker", frequency = 14.0 }}
        [[source]]
        position = [330.0, 230.0, 210.0]
        explosion = 2e12
        wavelet = {{ kind = "ricker", frequency = 14.0, delay = 0.1 }}
        [receivers]
        components = ["vx", "vy", "vz", "p"]
        # From 120 m in: at order 8 a velocity's interpolation reaches four nodes, and a
        # receiver reading into a layer makes the gradient there approximate.
        # The nodes a line of 8 from 120 to 360 m snapped to before positions were exact.
        [[receivers.point]]
        position = [120.0, 200.0, 0.0]
        [[receivers.point]]
        position = [150.0, 200.0, 0.0]
        [[receivers.point]]
        position = [190.0, 200.0, 0.0]
        [[receivers.point]]
        position = [220.0, 200.0, 0.0]
        [[receivers.point]]
        position = [260.0, 200.0, 0.0]
        [[receivers.point]]
        position = [290.0, 200.0, 0.0]
        [[receivers.point]]
        position = [330.0, 200.0, 0.0]
        [[receivers.point]]
        position = [360.0, 200.0, 0.0]
        [[receivers.point]]
        position = [240.0, 120.0, 180.0]
        [[receivers.point]]
        position = [300.0, 280.0, 90.0]
        "#
    )
}

/// Checks the gradient of the balanced L2 misfit against finite differences for the model
/// `config` makes with a guessed sphere, against data from another sphere. Every node outside
/// layers `pml` nodes thick is perturbed.
fn gradient_check(
    label: &str,
    pml: usize,
    config: &dyn Fn([f64; 3]) -> String,
    strategy: vibe_rations::gradient::Strategy,
) {
    use vibe_rations::config::Config;
    use vibe_rations::gradient::{balancing_weights, gradient_with, l2_misfit};
    use vibe_rations::model::Model;
    use vibe_rations::simulation::Simulation;

    let base = std::path::Path::new(".");
    let mut cfg = Config::from_toml(&config([2050.0, 1000.0, 2000.0])).unwrap();
    let model = Model::from_config(Simulation::grid(&cfg).unwrap(), &cfg.model, base).unwrap();
    // The time step and the layers' damping held at the model's, as an inversion holds them:
    // with dt from a Courant fraction both follow the fastest node, which a perturbation may
    // move (band-limited layers put it on an overshoot beside a slightly slower layer).
    let dt = Simulation::with_model(&cfg, model.clone())
        .unwrap()
        .solver
        .dt;
    cfg.time.dt = Some(dt);
    cfg.time.courant = None;
    cfg.boundary.pml_vp_max.get_or_insert(model.vp_range()[1]);
    let mut truth = Config::from_toml(&config([1900.0, 950.0, 1950.0])).unwrap();
    (truth.time, truth.boundary) = (cfg.time.clone(), cfg.boundary.clone());
    let mut observed = Simulation::new(&truth, base).unwrap();
    observed.run(None, false).unwrap();
    let observed = observed.receivers;
    let weights = balancing_weights(&observed.traces, &observed.components);
    let observed = observed.traces;
    let misfit = |traces: &[Vec<f32>]| l2_misfit(traces, &observed, &weights);
    let misfit_of = |model: Model| {
        let mut sim = Simulation::with_model(&cfg, model).unwrap();
        sim.run(None, false).unwrap();
        misfit(&sim.receivers.traces).unwrap().0
    };

    let mut sim = Simulation::with_model(&cfg, model.clone()).unwrap();
    assert_eq!(
        sim.receivers_next_to_layers(),
        0,
        "receivers next to the layers"
    );
    let run = gradient_with(&mut sim, |r| misfit(&r.traces), strategy).unwrap();
    let [g_vp, g_vs, g_rho] = run.gradient.velocities(&model);

    // Every node outside the layers, the fastest included: dt and the layers are held. Down
    // the gradient, by a random amount at each node: along random signs the predicted change
    // cancels (vp's tenfold here, once the fast layer is perturbed too) and the finite
    // differences' round-off, about 1e-7 of a unit of misfit, swamps it.
    let g = model.grid;
    let [nx, ny, nz] = g.n;
    let inside = |i: usize, j: usize, k: usize| {
        (pml..nx - pml).contains(&i) && (pml..ny - pml).contains(&j) && k < nz - pml
    };
    for (name, grad, seed) in [("vp", &g_vp, 3), ("vs", &g_vs, 4), ("rho", &g_rho, 5)] {
        let u: Vec<f64> = noise(g.interior_len(), seed)
            .iter()
            .zip(grad.iter())
            .map(|(v, d)| -v * d.signum())
            .collect();
        let values = |m: &Model| match name {
            "vp" => m.vp.clone(),
            "vs" => m.vs.clone(),
            _ => m.rho.clone(),
        };
        let base_values = values(&model);
        let mut mask = vec![0.0f64; g.interior_len()];
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let n = g.interior_idx(i, j, k);
                    if inside(i, j, k) {
                        mask[n] = u[n];
                    }
                }
            }
        }
        // The change the gradient predicts for m → m (1 + h u), per unit h.
        let predicted: f64 = (0..mask.len())
            .map(|n| grad[n] * base_values[n] as f64 * mask[n])
            .sum();
        let perturbed = |h: f64| {
            let mut m = model.clone();
            let target = match name {
                "vp" => &mut m.vp,
                "vs" => &mut m.vs,
                _ => &mut m.rho,
            };
            for (v, &w) in target.iter_mut().zip(&mask) {
                *v = (*v as f64 * (1.0 + h * w)) as f32;
            }
            m
        };
        // Central differences at h and h/2, then Richardson's extrapolation to cancel their
        // h² error: smaller steps drown in single-precision round-off first.
        let h = 1.6e-2;
        let [coarse, fine] =
            [h, h / 2.0].map(|h| (misfit_of(perturbed(h)) - misfit_of(perturbed(-h))) / (2.0 * h));
        let extrapolated = (4.0 * fine - coarse) / 3.0;
        let e = |fd: f64| (fd - predicted).abs() / predicted.abs();
        println!(
            "gradient, {label}, {name}: predicted {predicted:.6e}; finite differences \
             off by {:.1e} (h {h:.1e}) and {:.1e} (h/2), extrapolated by {:.1e}",
            e(coarse),
            e(fine),
            e(extrapolated)
        );
        assert!(
            e(extrapolated) < 1.5e-4,
            "{name}, {label}: {}",
            e(extrapolated)
        );
    }
}

/// The Born run is the forward run's derivative and the gradient its transpose: for a random
/// change of vp, of vs or of ρ at every node outside the layers, its traces match the traces'
/// central differences (extrapolated in the step, with Δt and the layers held), and
/// ⟨J δm, r⟩ = ⟨δm, Jᵀ r⟩ with r the Born traces themselves times random positive weights, so
/// that the terms can't cancel. In the gradient checks' models: a free surface, two
/// band-limited layers and a sphere, a force and an explosion, at orders 4 and 8; and rigid
/// walls, uneven spacing, water under the surface and a double couple.
///
/// Measured (vp, vs, ρ), against the extrapolated differences, the largest component's error
/// relative to its size: order 4 3.5e-5, 1.1e-5, 8.1e-5; order 8 4.3e-5, 1.1e-5, 8.6e-5; walls
/// and water 6.3e-5, 3.0e-5, 1.6e-4 (ρ's plain differences stall at about 1e-4: round-off).
/// The transposes: order 4 1.4e-8, 5.9e-9, 4.2e-9; order 8 1.9e-9, 5.1e-9, 2.3e-8; walls and
/// water 1.1e-7, 1.6e-7, 6.1e-8.
#[test]
fn the_born_run_is_the_derivative_and_the_gradient_its_transpose() {
    for order in [4, 8] {
        born_check(&format!("order {order}"), 8, &|sphere| {
            gradient_model(order, sphere)
        });
    }
    born_check("walls and water", 0, &walls_and_water);
}

/// [`the_born_run_is_the_derivative_and_the_gradient_its_transpose`] in the model `config`
/// makes, every node outside layers `pml` nodes thick changed.
fn born_check(label: &str, pml: usize, config: &dyn Fn([f64; 3]) -> String) {
    use vibe_rations::born::{Change, born};
    use vibe_rations::config::Config;
    use vibe_rations::gradient::gradient_with;
    use vibe_rations::model::Model;
    use vibe_rations::simulation::Simulation;

    let base = std::path::Path::new(".");
    let mut cfg = Config::from_toml(&config([2050.0, 1000.0, 2000.0])).unwrap();
    let model = Model::from_config(Simulation::grid(&cfg).unwrap(), &cfg.model, base).unwrap();
    let dt = Simulation::with_model(&cfg, model.clone())
        .unwrap()
        .solver
        .dt;
    cfg.time.dt = Some(dt);
    cfg.time.courant = None;
    cfg.boundary.pml_vp_max.get_or_insert(model.vp_range()[1]);
    let traces = |m: Model| {
        let mut sim = Simulation::with_model(&cfg, m).unwrap();
        sim.run(None, false).unwrap();
        sim.receivers.traces
    };
    let g = model.grid;
    let [nx, ny, nz] = g.n;
    let mut mask = vec![0.0f64; g.interior_len()];
    for k in 0..nz - pml {
        for j in pml..ny - pml {
            for i in pml..nx - pml {
                mask[g.interior_idx(i, j, k)] = 1.0;
            }
        }
    }
    for (p, seed) in [(0, 13), (1, 14), (2, 15)] {
        let name = ["vp", "vs", "rho"][p];
        let u: Vec<f64> = noise(g.interior_len(), seed)
            .iter()
            .zip(&mask)
            .map(|(v, m)| v * m)
            .collect();
        let values = |m: &Model| [&m.vp, &m.vs, &m.rho][p].clone();
        let start = values(&model);
        let mut speeds = [vec![0.0; u.len()], vec![0.0; u.len()], vec![0.0; u.len()]];
        speeds[p] = u.iter().zip(&start).map(|(u, &m)| u * m as f64).collect();
        let change = Change::from_velocities(&model, speeds.each_ref().map(Vec::as_slice));
        let mut sim = Simulation::with_model(&cfg, model.clone()).unwrap();
        let jv = born(&mut sim, &change).unwrap().traces;

        let perturbed = |h: f64| {
            let mut m = model.clone();
            let target = [&mut m.vp, &mut m.vs, &mut m.rho]
                .into_iter()
                .nth(p)
                .unwrap();
            for (v, &w) in target.iter_mut().zip(&u) {
                *v = (*v as f64 * (1.0 + h * w)) as f32;
            }
            m
        };
        let difference = |h: f64| -> Vec<Vec<f64>> {
            let (a, b) = (traces(perturbed(h)), traces(perturbed(-h)));
            a.iter()
                .zip(&b)
                .map(|(a, b)| {
                    a.iter()
                        .zip(b)
                        .map(|(a, b)| (*a as f64 - *b as f64) / (2.0 * h))
                        .collect()
                })
                .collect()
        };
        let h = 1.6e-2;
        let (coarse, fine) = (difference(h), difference(h / 2.0));
        // Each component's error relative to its own size: pressure is 10⁶ times velocity.
        let error = |fd: &[Vec<f64>]| {
            fd.iter()
                .zip(&jv)
                .map(|(fd, jv)| {
                    let diff: f64 = fd
                        .iter()
                        .zip(jv)
                        .map(|(a, &b)| (a - b as f64).powi(2))
                        .sum();
                    let size: f64 = jv.iter().map(|&b| (b as f64).powi(2)).sum();
                    (diff / size).sqrt()
                })
                .fold(0.0f64, f64::max)
        };
        let extrapolated: Vec<Vec<f64>> = coarse
            .iter()
            .zip(&fine)
            .map(|(c, f)| c.iter().zip(f).map(|(c, f)| (4.0 * f - c) / 3.0).collect())
            .collect();

        // r: the Born traces with random positive weights, each component scaled to unit size.
        let weights = noise(jv.iter().map(Vec::len).sum(), seed + 100);
        let mut w = weights.iter();
        let r: Vec<Vec<f32>> = jv
            .iter()
            .map(|t| {
                let size: f64 = t.iter().map(|&v| (v as f64).powi(2)).sum();
                t.iter()
                    .map(|&v| (v as f64 * w.next().unwrap() / size) as f32)
                    .collect()
            })
            .collect();
        let forward: f64 = jv
            .iter()
            .flatten()
            .zip(r.iter().flatten())
            .map(|(&a, &b)| a as f64 * b as f64)
            .sum();
        let mut sim = Simulation::with_model(&cfg, model.clone()).unwrap();
        let run = gradient_with(&mut sim, |_| Ok((0.0, r.clone())), Strategy::default()).unwrap();
        let backward = change.dot(&run.gradient);
        let transpose = (forward - backward).abs() / forward.abs();
        println!(
            "Born, {label}, {name}: against central differences {:.1e} (h {h:.1e}), {:.1e} \
             (h/2), extrapolated {:.1e}; ⟨Jδm, r⟩ {forward:.6e}, ⟨δm, Jᵀr⟩ {backward:.6e}, \
             difference {transpose:.1e}",
            error(&coarse),
            error(&fine),
            error(&extrapolated)
        );
        assert!(
            error(&extrapolated) < 5e-4,
            "{label}, {name}: {}",
            error(&extrapolated)
        );
        assert!(transpose < 1e-6, "{label}, {name}: {transpose}");
    }
}

/// The Hessian's product is the gradient's derivative, second-order terms and all: for a
/// random change of vp, of vs or of ρ at every node outside the layers, at a model whose data
/// misfit the observed (the gradient checks' guessed sphere against another), H δm matches
/// the central differences of the balanced L2 misfit's gradient (in vp, vs and ρ, extrapolated
/// in the step, with Δt and the layers held), and differs from the Gauss–Newton product
/// JᵀWJ δm, which leaves out the terms through the adjoint wavefield.
///
/// Measured, for changes of (vp, vs, ρ), the largest of the three parameters' errors relative
/// to its size: order 8 1.6e-5, 3.6e-6, 1.8e-5, with the Gauss–Newton product 4.4e-2, 1.7e-2,
/// 1.5e-2 apart; walls and water 9.5e-6, 8.2e-6, 6.3e-5, and 1.6e-2, 2.5e-2, 2.5e-2 apart. So
/// the second-order terms are right to a few 10⁻³ of their size or better.
#[test]
fn the_hessian_product_is_the_gradients_derivative() {
    hessian_check("order 8", 8, &|sphere| gradient_model(8, sphere));
    hessian_check("walls and water", 0, &walls_and_water);
}

/// [`the_hessian_product_is_the_gradients_derivative`] in the model `config` makes, every node
/// outside layers `pml` nodes thick changed.
fn hessian_check(label: &str, pml: usize, config: &dyn Fn([f64; 3]) -> String) {
    use vibe_rations::born::Change;
    use vibe_rations::config::Config;
    use vibe_rations::gradient::{balancing_weights, gradient_with, l2_misfit};
    use vibe_rations::hessian::hessian_product;
    use vibe_rations::model::Model;
    use vibe_rations::simulation::Simulation;

    let base = std::path::Path::new(".");
    let mut cfg = Config::from_toml(&config([2050.0, 1000.0, 2000.0])).unwrap();
    let model = Model::from_config(Simulation::grid(&cfg).unwrap(), &cfg.model, base).unwrap();
    let dt = Simulation::with_model(&cfg, model.clone())
        .unwrap()
        .solver
        .dt;
    cfg.time.dt = Some(dt);
    cfg.time.courant = None;
    cfg.boundary.pml_vp_max.get_or_insert(model.vp_range()[1]);
    let mut truth = Config::from_toml(&config([1900.0, 950.0, 1950.0])).unwrap();
    (truth.time, truth.boundary) = (cfg.time.clone(), cfg.boundary.clone());
    let mut observed = Simulation::new(&truth, base).unwrap();
    observed.run(None, false).unwrap();
    let observed = observed.receivers;
    let weights = balancing_weights(&observed.traces, &observed.components);
    let observed = observed.traces;
    let misfit = |traces: &[Vec<f32>]| l2_misfit(traces, &observed, &weights);
    // The L2 misfit's second derivative in the samples: w_c on component c.
    let curvature = |d: &[Vec<f32>]| -> anyhow::Result<Vec<Vec<f32>>> {
        Ok(d.iter()
            .zip(&weights)
            .map(|(d, &w)| d.iter().map(|&x| (w * x as f64) as f32).collect())
            .collect())
    };
    let gradient = |m: Model| {
        let mut sim = Simulation::with_model(&cfg, m.clone()).unwrap();
        let run = gradient_with(&mut sim, |r| misfit(&r.traces), Strategy::default()).unwrap();
        run.gradient.velocities(&m)
    };

    let g = model.grid;
    let [nx, ny, nz] = g.n;
    let mut mask = vec![0.0f64; g.interior_len()];
    for k in 0..nz - pml {
        for j in pml..ny - pml {
            for i in pml..nx - pml {
                mask[g.interior_idx(i, j, k)] = 1.0;
            }
        }
    }
    for (p, seed) in [(0, 23), (1, 24), (2, 25)] {
        let name = ["vp", "vs", "rho"][p];
        let u: Vec<f64> = noise(g.interior_len(), seed)
            .iter()
            .zip(&mask)
            .map(|(v, m)| v * m)
            .collect();
        let values = |m: &Model| [&m.vp, &m.vs, &m.rho][p].clone();
        let start = values(&model);
        let mut speeds = [vec![0.0; u.len()], vec![0.0; u.len()], vec![0.0; u.len()]];
        speeds[p] = u.iter().zip(&start).map(|(u, &m)| u * m as f64).collect();
        let speeds = speeds.each_ref().map(Vec::as_slice);
        let change = Change::from_velocities(&model, speeds);
        let mut sim = Simulation::with_model(&cfg, model.clone()).unwrap();
        let product = hessian_product(&mut sim, &change, |r| misfit(&r.traces), curvature, None)
            .unwrap()
            .velocities(&model, speeds);
        // Without the residuals: the Gauss–Newton product.
        let mut sim = Simulation::with_model(&cfg, model.clone()).unwrap();
        let zeros = |r: &vibe_rations::receiver::Receivers| {
            Ok((0.0, r.traces.iter().map(|t| vec![0.0; t.len()]).collect()))
        };
        let gauss_newton = hessian_product(&mut sim, &change, zeros, curvature, None)
            .unwrap()
            .velocities(&model, speeds);

        let perturbed = |h: f64| {
            let mut m = model.clone();
            let target = [&mut m.vp, &mut m.vs, &mut m.rho]
                .into_iter()
                .nth(p)
                .unwrap();
            for (v, &w) in target.iter_mut().zip(&u) {
                *v = (*v as f64 * (1.0 + h * w)) as f32;
            }
            m
        };
        let difference = |h: f64| -> Vec<Vec<f64>> {
            let (a, b) = (gradient(perturbed(h)), gradient(perturbed(-h)));
            a.iter()
                .zip(&b)
                .map(|(a, b)| a.iter().zip(b).map(|(a, b)| (a - b) / (2.0 * h)).collect())
                .collect()
        };
        let h = 1.6e-2;
        let (coarse, fine) = (difference(h), difference(h / 2.0));
        let extrapolated: Vec<Vec<f64>> = coarse
            .iter()
            .zip(&fine)
            .map(|(c, f)| c.iter().zip(f).map(|(c, f)| (4.0 * f - c) / 3.0).collect())
            .collect();
        // Each parameter's error relative to its own size, outside the layers.
        let error = |a: &[Vec<f64>]| {
            a.iter()
                .zip(&product)
                .map(|(a, b)| {
                    let diff: f64 = a
                        .iter()
                        .zip(b)
                        .zip(&mask)
                        .map(|((a, b), m)| m * (a - b).powi(2))
                        .sum();
                    let size: f64 = b.iter().zip(&mask).map(|(b, m)| m * b * b).sum();
                    (diff / size).sqrt()
                })
                .fold(0.0f64, f64::max)
        };
        println!(
            "Hessian, {label}, {name}: against central differences of the gradient {:.1e} \
             (h {h:.1e}), {:.1e} (h/2), extrapolated {:.1e}; Gauss–Newton {:.1e} apart",
            error(&coarse),
            error(&fine),
            error(&extrapolated),
            error(&gauss_newton)
        );
        assert!(
            error(&extrapolated) < 5e-4,
            "{label}, {name}: {}",
            error(&extrapolated)
        );
        assert!(error(&gauss_newton) > 5e-3, "{label}, {name}");
    }
}

/// A small inversion: two shots of a vertical and a horizontal force, recorded on the free
/// surface, of a model with a faster sphere, inverted for vp and vs from a homogeneous start,
/// with the top 60 m (the sources and the receivers) held fixed and the gradient
/// preconditioned by the illumination. L-BFGS must lower the misfit and bring the sphere's vs
/// closer; vp, whose wavelength (200 m) is twice the sphere, the S-dominated data hardly see.
#[test]
fn an_inversion_finds_a_sphere() {
    use vibe_rations::config::Config;
    use vibe_rations::inversion::{Inversion, Observed};
    use vibe_rations::model::Model;
    use vibe_rations::simulation::Simulation;

    let model = |sphere: &str| {
        format!(
            r#"
            [grid]
            points = [36, 32, 26]
            spacing = 10.0
            [time]
            duration = 0.3
            dt = 4.5e-4
            [boundary]
            pml_width = 8
            pml_vp_max = 3000.0
            [[model.layer]]
            top = 0.0
            vp = 2400.0
            vs = 1300.0
            rho = 2100.0
            {sphere}
            [[source]]
            position = [110.0, 160.0, 30.0]
            force = [0.0, 0.0, 1e9]
            wavelet = {{ kind = "ricker", frequency = 12.0 }}
            [[source]]
            position = [240.0, 150.0, 30.0]
            force = [1e9, 0.0, 0.0]
            wavelet = {{ kind = "ricker", frequency = 12.0 }}
            [receivers]
            components = ["vx", "vy", "vz"]
            [[receivers.line]]
            from = [90.0, 120.0, 0.0]
            to = [260.0, 120.0, 0.0]
            count = 18
            [[receivers.line]]
            from = [90.0, 190.0, 0.0]
            to = [260.0, 190.0, 0.0]
            count = 18
            [[receivers.line]]
            from = [180.0, 90.0, 0.0]
            to = [180.0, 220.0, 0.0]
            count = 14
            [inversion]
            observed = "unused"
            iterations = 3
            fixed_above = 60.0
            "#
        )
    };
    let sphere = "[[model.sphere]]\ncenter = [175.0, 155.0, 115.0]\nradius = 40.0\n\
                  vp = 2700.0\nvs = 1500.0\nrho = 2100.0";
    let base = std::path::Path::new(".");
    let truth = Config::from_toml(&model(sphere)).unwrap();
    let observed: Vec<Observed> = (0..truth.sources.len())
        .map(|s| {
            let mut sim = Simulation::new(&truth.shot(s).unwrap(), base).unwrap();
            sim.run(None, false).unwrap();
            sim.receivers.traces.into()
        })
        .collect();
    let true_model =
        Model::from_config(Simulation::grid(&truth).unwrap(), &truth.model, base).unwrap();

    let cfg = Config::from_toml(&model("")).unwrap();
    let settings = cfg.inversion.clone().unwrap();
    let start = Model::from_config(Simulation::grid(&cfg).unwrap(), &cfg.model, base).unwrap();
    let mut inversion = Inversion::new(&cfg, start.clone(), observed, &settings).unwrap();
    // The model error inside the sphere, relative to the anomaly.
    let inside: Vec<usize> = (0..true_model.vp.len())
        .filter(|&n| true_model.vp[n] != start.vp[n])
        .collect();
    let error = |m: &Model, values: fn(&Model) -> &Vec<f32>| {
        let num: f64 = inside
            .iter()
            .map(|&n| (values(m)[n] - values(&true_model)[n]).powi(2) as f64)
            .sum();
        let den: f64 = inside
            .iter()
            .map(|&n| (values(&start)[n] - values(&true_model)[n]).powi(2) as f64)
            .sum();
        (num / den).sqrt()
    };
    let first = inversion.misfit();
    println!(
        "inversion: {} unknowns, misfit {first:.4e}",
        inversion.unknowns()
    );
    for _ in 0..settings.iterations.unwrap() {
        let it = inversion
            .iterate()
            .unwrap()
            .expect("a step that lowers the misfit");
        println!(
            "inversion, iteration {}: misfit {:.1}% of the start, step {:.2}, {} trial(s), \
             largest change {:.1}%, sphere vp error {:.1}%, vs error {:.1}%",
            it.iteration,
            100.0 * it.misfit / first,
            it.step,
            it.trials,
            100.0 * it.largest_change,
            100.0 * error(&inversion.model, |m| &m.vp),
            100.0 * error(&inversion.model, |m| &m.vs),
        );
    }
    // Measured with the defaults: 17.9% of the misfit left, vs error 91.9%. Before the merge
    // (cubic receivers, Taylor coefficients, the time error kept): 17.9% and 91.8%
    // (5 iterations: 6.6%, 82.7%).
    let last = inversion.misfit() / first;
    let vs = error(&inversion.model, |m| &m.vs);
    assert!(last < 0.2 && vs < 0.93, "misfit {last}, vs error {vs}");
}

/// The inversion's gradient, with respect to its own unknowns (relative changes of vs, ρ and
/// vp, listed in that order, summed over two shots that record velocity and pressure), is the
/// derivative of its misfit, checked down the gradient by extrapolated central differences.
///
/// Measured with the defaults (optimised coefficients, time error removed, band-limited
/// layers): 4.9e-6 at h = 1e-2, between round-off in the single-precision traces at smaller
/// steps (5.7e-7 at 5e-3, 1.4e-5 at 2.5e-3, scattering rather than shrinking) and the
/// extrapolation's h⁴ error at larger ones (7.5e-5 at 2e-2). Node-sampled (until 2026-11-22):
/// 4.0e-6 at 1e-2; 1.4e-5 at 8e-3, 5.1e-5 at 4e-3, 2.7e-5 at 2e-3; 2.5e-5 at 1.6e-2. Before the merge with
/// stencil-long receiver interpolation (cubic receivers, Taylor coefficients, time error
/// kept): 1.1e-5 at 4e-3, with the surface line reaching x = 70 and 220 m; with eight-point
/// interpolation its end receivers read the side layers' damped values and the check gave
/// 3.1e-5, so the line now keeps clear of them. Along a random direction the predicted change
/// cancelled to 4e-4 of this one and round-off held the check at 1e-2.
#[test]
fn the_inversion_gradient_matches_finite_differences() {
    let e = inversion_gradient_check("inversion gradient", "", 1e-2);
    assert!(e < 3e-5, "{e}");
}

/// The other misfits' gradients, through the same inversion and check: the L2 misfit of
/// low-pass filtered traces (a frequency-continuation stage), the cross-correlation time shifts
/// of filtered traces, and the time–frequency phase misfit.
///
/// Measured (extrapolated, at the step given), with band-limited layers: filtered L2 2.4e-6
/// (h = 1e-2; 1.9e-4 at 2e-2, 2.3e-5 at 5e-3, 6.2e-5 at 2.5e-3: h⁴ above, round-off below);
/// phase 6.3e-7 (h = 1e-2); time shifts 1.2e-4 (h = 1e-3). Node-sampled (until 2026-11-22):
/// 4.4e-7, 4.2e-6 and 1.3e-4. The time shifts are checked at a smaller step and against a looser bound: each
/// trace's shift is only piecewise smooth, its derivative jumping where the correlation's peak
/// moves to the next lag, and steps from 2e-3 up cross such jumps (8.5e-4 at 2e-3, 2.0e-3 at
/// 3e-3); below that the check sits on a floor of round-off in the shifts of single-precision
/// traces, 1.3e-4 to 2.2e-4 for h from 2.5e-4 to 1e-3 without a trend. A wrong derivative
/// would be off by tens of per cent. (Before the filter padded the traces' ends: filtered L2
/// 1.5e-5, filtered time shifts 7.5e-6 at 2e-3, by luck clear of the jumps.)
#[test]
fn the_misfits_gradients_match_finite_differences() {
    for (label, settings, h, tolerance) in [
        (
            "filtered l2",
            "[[inversion.stage]]\nmax_frequency = 10.0\niterations = 1",
            1e-2,
            3e-5,
        ),
        (
            "filtered traveltime",
            "misfit = \"traveltime\"\nmax_shift = 0.05\n\
             [[inversion.stage]]\nmax_frequency = 10.0\niterations = 1",
            1e-3,
            5e-4,
        ),
        ("phase", "misfit = \"phase\"", 1e-2, 3e-5),
    ] {
        let e = inversion_gradient_check(label, settings, h);
        assert!(e < tolerance, "{label}: {e}");
    }
}

/// Mutes keep the gradient exact: they weigh observed and synthetic traces alike, before the
/// filter, so the adjoint source goes back through the filter and them in turn. L2 of
/// low-passed traces with the late arrivals down-weighted to a fifth, and the time shifts of
/// one window (which changes what they measure: the predicted change down the gradient is
/// -9.47e-5 against -1.41e-4 without it).
///
/// Measured (extrapolated), with band-limited layers: muted, filtered L2 1.9e-5 (h = 1e-2;
/// 1.3e-4 at 2e-2, 2.5e-5 at 5e-3, 5.3e-5 at 2.5e-3: the extrapolation's h⁴ error above,
/// round-off growing as 1/h below, about equal at 1e-2); windowed time shifts 9.1e-5
/// (h = 1e-3), the shifts' round-off floor. Node-sampled (until 2026-11-22): 8.7e-7, a lucky
/// cancellation within the same band, and 2.6e-5 (1.0e-5 at 2e-3, 1.4e-4 at 5e-4, 3.0e-5 to
/// 9.9e-5 without a window). With the mute after the filter, as first written, the L2 check gave
/// 1.5e-5 at h = 1e-2 and 2.2e-5 to 3.8e-5 at smaller steps; both orders transpose exactly
/// (a unit test), and why the order changed the check's noise was not traced.
#[test]
fn the_mutes_keep_the_gradient_exact() {
    for (label, settings, h, tolerance) in [
        (
            "filtered l2, late arrivals down-weighted",
            "[[inversion.stage]]\nmax_frequency = 10.0\niterations = 1\n\
             [[inversion.mute]]\nkeep = \"before\"\nvelocity = 1300.0\ndelay = 0.15\n\
             taper = 0.03\nweight = 0.2",
            1e-2,
            3e-5,
        ),
        (
            "windowed time shifts",
            "misfit = \"traveltime\"\nmax_shift = 0.05\n\
             [[inversion.mute]]\nkeep = \"after\"\nvelocity = 3000.0\ndelay = 0.06\ntaper = 0.02\n\
             [[inversion.mute]]\nkeep = \"before\"\nvelocity = 1500.0\ndelay = 0.12\ntaper = 0.02",
            1e-3,
            5e-4,
        ),
    ] {
        let e = inversion_gradient_check(label, settings, h);
        assert!(e < tolerance, "{label}: {e}");
    }
}

/// The relative error of the extrapolated central difference of the inversion's misfit down
/// its gradient, at steps `h` and `h`/2, with `settings` added to its [inversion] section.
fn inversion_gradient_check(label: &str, settings: &str, h: f64) -> f64 {
    use vibe_rations::config::Config;
    use vibe_rations::inversion::{Inversion, Observed};
    use vibe_rations::model::Model;
    use vibe_rations::simulation::Simulation;

    let model = |sphere_vs: f64| {
        format!(
            r#"
            [grid]
            points = [30, 28, 24]
            spacing = 10.0
            [time]
            duration = 0.25
            dt = 5e-4
            [boundary]
            pml_width = 6
            pml_vp_max = 3000.0
            [output]
            record_every = 2
            [[model.layer]]
            top = 0.0
            vp = 2200.0
            vs = 1200.0
            rho = 2000.0
            [[model.layer]]
            top = 120.0
            vp = 2600.0
            vs = 1450.0
            rho = 2200.0
            [[model.sphere]]
            center = [150.0, 140.0, 80.0]
            radius = 30.0
            vp = 2300.0
            vs = {sphere_vs}
            rho = 2050.0
            [[source]]
            position = [90.0, 140.0, 40.0]
            force = [0.0, 3e8, 1e9]
            wavelet = {{ kind = "ricker", frequency = 14.0 }}
            [[source]]
            position = [210.0, 120.0, 50.0]
            explosion = 1e12
            wavelet = {{ kind = "ricker", frequency = 14.0 }}
            [receivers]
            components = ["vx", "vz", "p"]
            # The nodes a line of 6 from 100 to 190 m snapped to before positions were
            # exact: exactly, one trace's time shift changes lag at the start (2026-11-22).
            [[receivers.point]]
            position = [100.0, 140.0, 0.0]
            [[receivers.point]]
            position = [120.0, 140.0, 0.0]
            [[receivers.point]]
            position = [140.0, 140.0, 0.0]
            [[receivers.point]]
            position = [150.0, 140.0, 0.0]
            [[receivers.point]]
            position = [170.0, 140.0, 0.0]
            [[receivers.point]]
            position = [190.0, 140.0, 0.0]
            [[receivers.point]]
            position = [150.0, 100.0, 120.0]
            [inversion]
            observed = "unused"
            parameters = ["vs", "rho", "vp"]
            {settings}
            "#
        )
    };
    let base = std::path::Path::new(".");
    let truth = Config::from_toml(&model(1050.0)).unwrap();
    let observed: Vec<Observed> = (0..truth.sources.len())
        .map(|s| {
            let mut sim = Simulation::new(&truth.shot(s).unwrap(), base).unwrap();
            sim.run(None, false).unwrap();
            sim.receivers.traces.into()
        })
        .collect();
    let cfg = Config::from_toml(&model(1200.0)).unwrap();
    for s in 0..cfg.sources.len() {
        let shot = Simulation::new(&cfg.shot(s).unwrap(), base).unwrap();
        assert_eq!(shot.receivers_next_to_layers(), 0, "shot {s}");
    }
    let settings = cfg.inversion.clone().unwrap();
    let start = Model::from_config(Simulation::grid(&cfg).unwrap(), &cfg.model, base).unwrap();
    let inversion = Inversion::new(&cfg, start.clone(), observed, &settings).unwrap();
    let g = inversion.evaluate(&start).unwrap().gradient;
    // Down the gradient, scaled to a largest change of 1: the change it predicts can't cancel
    // between unknowns, as along a random direction, where single-precision round-off in the
    // misfit then swamps it; a gradient with wrong parts still disagrees with the finite
    // difference, which measures g_true · g, not |g|².
    let largest = g.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    let u: Vec<f64> = g.iter().map(|v| -v / largest).collect();
    let predicted: f64 = g.iter().zip(&u).map(|(g, u)| g * u).sum();
    let misfit = |h: f64| {
        let x: Vec<f64> = u.iter().map(|u| h * u).collect();
        inversion.misfit_of(&inversion.model_at(&x)).unwrap()
    };
    let [coarse, fine] = [h, h / 2.0].map(|h| (misfit(h) - misfit(-h)) / (2.0 * h));
    let extrapolated = (4.0 * fine - coarse) / 3.0;
    let e = |fd: f64| (fd - predicted).abs() / predicted.abs();
    println!(
        "{label}: predicted {predicted:.6e}; finite differences off by {:.1e} and {:.1e} \
         (h = {h:.0e}), extrapolated by {:.1e}",
        e(coarse),
        e(fine),
        e(extrapolated)
    );
    e(extrapolated)
}
