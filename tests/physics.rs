//! Whole simulations checked against closed-form solutions and against each other.
//!
//! The homogeneous full-space solutions are from Aki & Richards (2002), §4.2–4.3. Each test
//! keeps receivers clear of the absorbing layers and stops before their (small) reflections
//! could matter, unless reflections are what it measures.

mod common;

use common::*;
use vibe_rations::receiver::Component;
use vibe_rations::source::Wavelet;

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
        println!("explosion, r = {dist:.1} m: misfit {:.2}%", 100.0 * e);
        assert!(e < 0.02, "misfit {e} at {p:?}");
    }
}

#[test]
fn a_point_force_matches_stokes_solution() {
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
            "#
        ),
    ));
    let w = Wavelet::Ricker { f0: 10.0, t0: 0.12 };
    // Forces enter with velocity, at whole steps; the velocity samples are half a step later.
    let times = sim.sample_times();
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
            println!(
                "force, offset {d:?}, {}: misfit {:.2}%",
                c.name(),
                100.0 * e
            );
            assert!(e < 0.01, "misfit {e} in {} at {p:?}", c.name());
        }
    }
}

/// At the epicentre a P wave meets the free surface head-on, and the surface moves twice as
/// much as the same wave would in a full space.
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
    assert!(e < 0.025, "misfit {e}");
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

/// Measured: order 4 at courant 0.85, 1.5e-4 of the peak; order 8 at 0.3 (the default),
/// 3.9e-4 (about -68 dB): longer stencils reach further across the layers' inner edge.
#[test]
fn absorbing_layers_reflect_little() {
    let worst = layer_reflection(false);
    assert!(worst < 4.5e-4, "{worst}");
}

/// Along a free surface the side layers also damp plainly, to stay stable, which costs some
/// reflection. Measured: 1.8e-3 at order 4, 2.1e-3 at order 8 (the default).
#[test]
fn absorbing_layers_reflect_little_along_a_free_surface() {
    let worst = layer_reflection(true);
    assert!(worst < 2.4e-3, "{worst}");
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

/// The adjoint run is the transpose of the forward run. For a force history w at A, recorded
/// as d = F w at receivers in every component, and residuals d̃ at those receivers,
/// ⟨F w, d̃⟩ = ⟨w, Fᵀ d̃⟩, where Fᵀ d̃ is the adjoint velocity read at A after each adjoint step.
/// The scheme is self-adjoint everywhere but in the absorbing layers, so this holds to
/// round-off for receivers outside them, whatever the waves do inside.
#[test]
fn the_adjoint_run_is_the_transpose_of_the_forward_run() {
    for order in [4, 8] {
        for (walls, free_surface, signal, every) in [
            ("rigid", true, "white", 1),
            ("cpml", true, "band-limited", 2),
            ("cpml", false, "white", 1),
        ] {
            let e = dot_product(order, walls, free_surface, signal, every);
            // Measured 2.3e-8 to 8.3e-8: single precision over about 1300 steps. Leaving the
            // adjoint's initial pressure state unmirrored gives 2.4e-7 to 2.8e-6.
            assert!(e < 1.5e-7, "order {order}, {walls}, {signal}: {e}");
        }
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
fn dot_product(order: usize, walls: &str, free_surface: bool, signal: &str, every: usize) -> f64 {
    let pml = if walls == "rigid" { 0 } else { 10 };
    // Without a free surface the top 100 m are an absorbing layer.
    let shallow = if free_surface { 0.0 } else { 110.0 };
    let receivers: String = [
        [180.0, 220.0, shallow],
        [420.0, 260.0, 230.0],
        [300.0, 150.0, shallow + 60.0],
        [250.0, 300.0, 140.0],
        // One node below the free surface, where the stresses' mirror images reach.
        [400.0, 120.0, shallow + 10.0],
    ]
    .iter()
    .map(|p| format!("[[receivers.point]]\nposition = {p:?}\n"))
    .collect();
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
        position = [330.0, 240.0, 200.0]
        force = [0.6, 0.0, 0.8]
        wavelet = {{ kind = "ricker", frequency = 12.0 }}
        [receivers]
        components = ["vx", "vy", "vz", "p"]
        {receivers}
        "#
    ));
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
    let len = sim.receivers.nodes.len() * sim.receivers.samples;
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
    let node = sim.sources[0].node;
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
         every {every}: ⟨Fw, d̃⟩ {forward:.6e}, ⟨w, Fᵀd̃⟩ {backward:.6e}, difference {e:.1e}"
    );
    e
}

/// The gradient is that of the solver's forward run: perturbing vp, vs or ρ at every node
/// outside the absorbing layers by random amounts changes the misfit as the gradient says,
/// to within the finite difference's own error. The model has a free surface, two layers and
/// a sphere; a force and an explosion; receivers on the surface and at depth recording every
/// component, every second step.
///
/// Measured, extrapolated (orders 4 / 8): vp 1.7e-4 / 6.5e-5, vs 1.6e-6 / 5.8e-6, ρ 1.5e-5 /
/// 1.2e-5. The plain central differences are off by 2e-3 to 2e-1 at these steps, falling 4×
/// per halving: their own h² error, with nothing beneath it.
#[test]
fn the_gradient_matches_finite_differences() {
    for order in [4, 8] {
        gradient_check(order);
    }
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
        [[receivers.line]]
        from = [100.0, 200.0, 0.0]
        to = [380.0, 200.0, 0.0]
        count = 8
        [[receivers.point]]
        position = [240.0, 120.0, 180.0]
        [[receivers.point]]
        position = [300.0, 280.0, 90.0]
        "#
    )
}

fn gradient_check(order: usize) {
    use vibe_rations::config::Config;
    use vibe_rations::gradient::{balancing_weights, gradient, l2_misfit};
    use vibe_rations::model::Model;
    use vibe_rations::simulation::Simulation;

    let observed = simulate(&gradient_model(order, [1900.0, 950.0, 1950.0]))
        .receivers
        .traces;
    let weights = balancing_weights(&observed);
    let misfit = |traces: &[Vec<f32>]| l2_misfit(traces, &observed, &weights);
    let cfg = Config::from_toml(&gradient_model(order, [2050.0, 1000.0, 2000.0])).unwrap();
    let base = std::path::Path::new(".");
    let model = Model::from_config(Simulation::grid(&cfg).unwrap(), &cfg.model, base).unwrap();
    let misfit_of = |model: Model| {
        let mut sim = Simulation::with_model(&cfg, model).unwrap();
        sim.run(None, false).unwrap();
        misfit(&sim.receivers.traces).unwrap().0
    };

    let mut sim = Simulation::with_model(&cfg, model.clone()).unwrap();
    let run = gradient(&mut sim, |r| misfit(&r.traces), None).unwrap();
    let [g_vp, g_vs, g_rho] = run.gradient.velocities(&model);

    // Every node outside the layers; vp only below the fastest, which sets dt and the layers.
    let g = model.grid;
    let [nx, ny, nz] = g.n;
    let pml = 8;
    let vp_max = model.vp_range()[1] as f32;
    let inside = |i: usize, j: usize, k: usize| {
        (pml..nx - pml).contains(&i) && (pml..ny - pml).contains(&j) && k < nz - pml
    };
    for (name, grad, seed) in [("vp", &g_vp, 3), ("vs", &g_vs, 4), ("rho", &g_rho, 5)] {
        let u: Vec<f64> = noise(g.interior_len(), seed)
            .iter()
            .map(|v| 2.0 * v - 1.0)
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
                    if inside(i, j, k) && !(name == "vp" && base_values[n] >= vp_max) {
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
            "gradient, order {order}, {name}: predicted {predicted:.6e}; finite differences \
             off by {:.1e} (h {h:.1e}) and {:.1e} (h/2), extrapolated by {:.1e}",
            e(coarse),
            e(fine),
            e(extrapolated)
        );
        assert!(
            e(extrapolated) < 3e-4,
            "{name}, order {order}: {}",
            e(extrapolated)
        );
    }
}
