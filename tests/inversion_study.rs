//! How the misfits, frequency-continuation stages and preconditioners invert small surveys:
//! `just study`, or `just study <name>` for one.
//!
//! Not part of the regular suite (about 30 minutes on 4 cores). Each variant inverts the same
//! observed traces from the same start, and prints after each stage the plain balanced L2
//! misfit (one scale for every variant, relative to the start's) and the model error over the
//! target, relative to the start's error. The research log records the results.
//!
//! - `misfits_on_a_sphere`: `an_inversion_finds_a_sphere`'s survey and start, 6 iterations:
//!   a small anomaly, the start within a fraction of a period of the data.
//! - `density_scaling`: the same, inverted for ρ too, with ρ's steps scaled down.
//! - `receivers_side`: the same with nothing fixed, preconditioned by one or both sides of
//!   the illumination.
//! - `gauss_newton_diagonal`: the sphere preconditioned by the Gauss–Newton Hessian's
//!   diagonal, estimated from random adjoint sources (`hessian_probes`), against the
//!   illumination: for vs, ρ and vp with a denser core; for vp and vs; and with nothing fixed.
//! - `slow_start`: a homogeneous model, started 20% too slow in vp and vs, a line of 43
//!   receivers over 420 m, 8 iterations: at the far offsets the P waves come about half a
//!   period late and the S waves about a full period, where L2 cycle-skips. The target is
//!   every free node below 60 m.

use std::sync::Arc;
use std::time::Instant;

use vibe_rations::config::Config;
use vibe_rations::inversion::{Inversion, Observed};
use vibe_rations::misfit::{Measure, Misfit};
use vibe_rations::model::Model;
use vibe_rations::simulation::Simulation;

/// A variant's name and the settings added to its `[inversion]` section.
type Variant = (&'static str, &'static str);

/// Inverts `observed` from `cfg`'s model with each variant's settings appended to `cfg`'s
/// text, printing the L2 misfit and the target's vs and vp errors after each stage.
fn study(
    survey: &str,
    text: &dyn Fn(&str) -> String,
    observed: &[Observed],
    truth: &Model,
    target: &dyn Fn(&Model) -> Vec<usize>,
    variants: &[Variant],
) {
    let base = std::path::Path::new(".");
    for &(name, settings) in variants {
        let cfg = Config::from_toml(&text(settings)).unwrap();
        let inversion = cfg.inversion.clone().unwrap();
        let start = Model::from_config(Simulation::grid(&cfg).unwrap(), &cfg.model, base).unwrap();
        let nodes = target(&start);
        let error = |m: &Model, values: fn(&Model) -> &Vec<f32>| {
            let sum = |of: &Model| -> f64 {
                nodes
                    .iter()
                    .map(|&n| (values(of)[n] - values(truth)[n]).powi(2) as f64)
                    .sum()
            };
            (sum(m) / sum(&start)).sqrt()
        };
        let t0 = Instant::now();
        let mut inv = Inversion::new(&cfg, start.clone(), observed.to_vec(), &inversion).unwrap();
        // The balanced L2 misfit of a model, the same for every variant.
        let l2 = |m: &Model, inv: &Inversion| -> f64 {
            inv.shots
                .iter()
                .map(|shot| {
                    let mut sim = Simulation::with_model(&shot.config, m.clone()).unwrap();
                    sim.run(None, false).unwrap();
                    let r = &sim.receivers;
                    let interval = sim.info.sample_interval;
                    let measure = Measure::l2(true);
                    let obs = Arc::clone(&shot.observed);
                    Misfit::new(measure, &obs, &r.components, r.samples, interval, &[])
                        .unwrap()
                        .value(&r.traces)
                        .unwrap()
                })
                .sum()
        };
        let first = l2(&start, &inv);
        let stages = inv.stages().to_vec();
        let mut done = 0;
        for (s, stage) in stages.iter().enumerate() {
            if s > 0 {
                inv.begin_stage(s).unwrap();
            }
            for _ in 0..stage.iterations {
                if inv.iterate().unwrap().is_none() {
                    println!("{survey}, {name}: no step lowers the misfit in stage {s}");
                    break;
                }
                done += 1;
            }
            // ρ's error where the truth has a ρ anomaly; otherwise its rms relative change over
            // the target, which is all cross-talk.
            let rho = if nodes.iter().any(|&n| truth.rho[n] != start.rho[n]) {
                format!("rho error {:.1}%", 100.0 * error(&inv.model, |m| &m.rho))
            } else {
                let change: f64 = nodes
                    .iter()
                    .map(|&n| (inv.model.rho[n] as f64 / start.rho[n] as f64 - 1.0).powi(2))
                    .sum();
                format!(
                    "rho change {:.2}%",
                    100.0 * (change / nodes.len() as f64).sqrt()
                )
            };
            println!(
                "{survey}, {name}: after stage {s} ({done} iterations): L2 {:.1}%, vs error \
                 {:.1}%, vp error {:.1}%, {rho}",
                100.0 * l2(&inv.model, &inv) / first,
                100.0 * error(&inv.model, |m| &m.vs),
                100.0 * error(&inv.model, |m| &m.vp),
            );
        }
        println!("{survey}, {name}: {:.0} s", t0.elapsed().as_secs_f64());
    }
}

/// Each shot of `truth` run, its traces as observed data.
fn record(truth: &Config) -> (Vec<Observed>, Model) {
    let base = std::path::Path::new(".");
    let observed = (0..truth.sources.len())
        .map(|s| {
            let mut sim = Simulation::new(&truth.shot(s).unwrap(), base).unwrap();
            sim.run(None, false).unwrap();
            sim.receivers.traces.into()
        })
        .collect();
    let model = Model::from_config(Simulation::grid(truth).unwrap(), &truth.model, base).unwrap();
    (observed, model)
}

fn sphere_survey(sphere: &str, inversion: &str) -> String {
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
        position = [110.0, 155.0, 30.0]
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
        from = [175.0, 90.0, 0.0]
        to = [175.0, 220.0, 0.0]
        count = 14
        [inversion]
        observed = "unused"
        fixed_above = 60.0
        {inversion}
        "#
    )
}

#[test]
#[ignore]
fn misfits_on_a_sphere() {
    let anomaly = "[[model.sphere]]\ncenter = [175.0, 155.0, 115.0]\nradius = 40.0\n\
                   vp = 2700.0\nvs = 1500.0\nrho = 2100.0";
    let (observed, truth) = record(&Config::from_toml(&sphere_survey(anomaly, "")).unwrap());
    study(
        "sphere",
        &|settings| sphere_survey("", settings),
        &observed,
        &truth,
        &|start| {
            (0..truth.vs.len())
                .filter(|&n| truth.vs[n] != start.vs[n])
                .collect()
        },
        &[
            ("l2", "iterations = 6"),
            (
                "l2, 6 Hz then all",
                "[[inversion.stage]]\nmax_frequency = 6.0\niterations = 3\n\
                 [[inversion.stage]]\niterations = 3",
            ),
            (
                "time shifts, then l2",
                "misfit = \"traveltime\"\nmax_shift = 0.05\n[[inversion.stage]]\niterations = 3\n\
                 [[inversion.stage]]\nmisfit = \"l2\"\niterations = 3",
            ),
            ("phase", "misfit = \"phase\"\niterations = 6"),
            (
                "phase, then l2",
                "misfit = \"phase\"\n[[inversion.stage]]\niterations = 3\n\
                 [[inversion.stage]]\nmisfit = \"l2\"\niterations = 3",
            ),
        ],
    );
}

/// The sphere inverted for vs, ρ and vp with density's share of the steps scaled down
/// (`[inversion.scale]`): once with a denser core, and once with none, where any change of ρ
/// is cross-talk.
#[test]
#[ignore]
fn density_scaling() {
    let variants: [Variant; 3] = [
        (
            "rho x 1",
            "parameters = [\"vs\", \"rho\", \"vp\"]\niterations = 5",
        ),
        (
            "rho x 0.3",
            "parameters = [\"vs\", \"rho\", \"vp\"]\niterations = 5\n[inversion.scale]\nrho = 0.3",
        ),
        (
            "rho x 0.1",
            "parameters = [\"vs\", \"rho\", \"vp\"]\niterations = 5\n[inversion.scale]\nrho = 0.1",
        ),
    ];
    for (label, rho) in [("denser core", 2300.0), ("no density anomaly", 2100.0)] {
        let anomaly = format!(
            "[[model.sphere]]\ncenter = [175.0, 155.0, 115.0]\nradius = 40.0\n\
             vp = 2700.0\nvs = 1500.0\nrho = {rho}"
        );
        let (observed, truth) = record(&Config::from_toml(&sphere_survey(&anomaly, "")).unwrap());
        study(
            &format!("sphere, {label}"),
            &|settings| sphere_survey("", settings),
            &observed,
            &truth,
            &|start| {
                (0..truth.vs.len())
                    .filter(|&n| truth.vs[n] != start.vs[n])
                    .collect()
            },
            &variants,
        );
    }
}

/// The sphere with nothing held fixed, where the sources' and receivers' surroundings would
/// take the updates: preconditioned by the sources' side of the illumination alone, then by
/// both sides, estimated from 2 and 8 runs; and with the top fixed, by both.
#[test]
#[ignore]
fn receivers_side() {
    let anomaly = "[[model.sphere]]\ncenter = [175.0, 155.0, 115.0]\nradius = 40.0\n\
                   vp = 2700.0\nvs = 1500.0\nrho = 2100.0";
    let (observed, truth) = record(&Config::from_toml(&sphere_survey(anomaly, "")).unwrap());
    study(
        "receivers' side",
        &|settings| {
            let fixed = settings.contains("fixed");
            let text = sphere_survey("", &settings.replace("fixed", ""));
            if fixed {
                text
            } else {
                text.replace("fixed_above = 60.0", "")
            }
        },
        &observed,
        &truth,
        &|start| {
            (0..truth.vs.len())
                .filter(|&n| truth.vs[n] != start.vs[n])
                .collect()
        },
        &[
            ("sources' side", "iterations = 5"),
            ("both sides, 2 runs", "iterations = 5\nreceiver_runs = 2"),
            ("both sides, 8 runs", "iterations = 5\nreceiver_runs = 8"),
            ("top fixed, sources' side", "iterations = 5\n# fixed"),
            (
                "top fixed, both sides, 8 runs",
                "iterations = 5\nreceiver_runs = 8\n# fixed",
            ),
        ],
    );
}

/// The sphere preconditioned by the Gauss–Newton Hessian's diagonal instead of the
/// illumination: from 8 or 16 probes, smoothed over 0, 10 or 20 m, with the water level at its
/// default or ten times it, and balanced parameter by parameter (`hessian_balance`, so that
/// only its shape differs from the illumination's). Three cases: vs, ρ and vp with a denser core (`density_scaling`'s
/// first), where the diagonal also sets the parameters' shares; vp and vs
/// (`misfits_on_a_sphere`'s); and vp and vs with nothing fixed (`receivers_side`'s), where
/// the diagonal holds the receivers' side.
#[test]
#[ignore]
fn gauss_newton_diagonal() {
    let sphere = |rho: f64| {
        format!(
            "[[model.sphere]]\ncenter = [175.0, 155.0, 115.0]\nradius = 40.0\n\
             vp = 2700.0\nvs = 1500.0\nrho = {rho}"
        )
    };
    let inside = |truth: &Model| {
        let truth = truth.clone();
        move |start: &Model| -> Vec<usize> {
            (0..truth.vs.len())
                .filter(|&n| truth.vs[n] != start.vs[n])
                .collect()
        }
    };
    let (observed, truth) =
        record(&Config::from_toml(&sphere_survey(&sphere(2300.0), "")).unwrap());
    // vs, ρ and vp, 5 iterations, and each variant's settings.
    macro_rules! three {
        ($extra:literal) => {
            concat!(
                "parameters = [\"vs\", \"rho\", \"vp\"]\niterations = 5\n",
                $extra
            )
        };
    }
    let variants: [Variant; 7] = [
        ("illumination", three!("")),
        ("8 probes", three!("hessian_probes = 8")),
        (
            "8 probes, 10 m",
            three!("hessian_probes = 8\nhessian_smoothing = 10.0"),
        ),
        (
            "8 probes, 20 m",
            three!("hessian_probes = 8\nhessian_smoothing = 20.0"),
        ),
        (
            "16 probes, 10 m",
            three!("hessian_probes = 16\nhessian_smoothing = 10.0"),
        ),
        (
            "8 probes, 10 m, water 1e-2",
            three!("hessian_probes = 8\nhessian_smoothing = 10.0\nwater_level = 1e-2"),
        ),
        (
            "8 probes, 10 m, balanced",
            three!("hessian_probes = 8\nhessian_smoothing = 10.0\nhessian_balance = true"),
        ),
    ];
    study(
        "sphere, vs, rho and vp, denser core",
        &|settings| sphere_survey("", settings),
        &observed,
        &truth,
        &inside(&truth),
        &variants,
    );

    let (observed, truth) =
        record(&Config::from_toml(&sphere_survey(&sphere(2100.0), "")).unwrap());
    study(
        "sphere, vp and vs",
        &|settings| sphere_survey("", settings),
        &observed,
        &truth,
        &inside(&truth),
        &[
            ("illumination", "iterations = 5"),
            (
                "8 probes, 10 m",
                "iterations = 5\nhessian_probes = 8\nhessian_smoothing = 10.0",
            ),
            (
                "8 probes, 10 m, balanced",
                "iterations = 5\nhessian_probes = 8\nhessian_smoothing = 10.0\n\
                 hessian_balance = true",
            ),
        ],
    );
    study(
        "sphere, vp and vs, nothing fixed",
        &|settings| sphere_survey("", settings).replace("fixed_above = 60.0", ""),
        &observed,
        &truth,
        &inside(&truth),
        &[
            ("sources' side", "iterations = 5"),
            ("both sides, 8 runs", "iterations = 5\nreceiver_runs = 8"),
            (
                "8 probes, 10 m",
                "iterations = 5\nhessian_probes = 8\nhessian_smoothing = 10.0",
            ),
        ],
    );
}

fn slow_survey(vs: f64, vp: f64, inversion: &str) -> String {
    format!(
        r#"
        [grid]
        points = [60, 36, 28]
        spacing = 10.0
        [time]
        duration = 0.5
        dt = 4.5e-4
        [boundary]
        pml_width = 8
        pml_vp_max = 3000.0
        [[model.layer]]
        top = 0.0
        vp = {vp}
        vs = {vs}
        rho = 2100.0
        [[source]]
        position = [100.0, 180.0, 30.0]
        force = [0.0, 0.0, 1e9]
        wavelet = {{ kind = "ricker", frequency = 12.0 }}
        [[source]]
        position = [500.0, 170.0, 30.0]
        force = [1e9, 0.0, 0.0]
        wavelet = {{ kind = "ricker", frequency = 12.0 }}
        [receivers]
        components = ["vx", "vy", "vz"]
        [[receivers.line]]
        from = [90.0, 160.0, 0.0]
        to = [510.0, 160.0, 0.0]
        count = 43
        [inversion]
        observed = "unused"
        parameters = ["vp", "vs"]
        {inversion}
        "#
    )
}

#[test]
#[ignore]
fn slow_start() {
    let (observed, truth) = record(&Config::from_toml(&slow_survey(1300.0, 2400.0, "")).unwrap());
    study(
        "slow start",
        &|settings| slow_survey(1040.0, 1920.0, settings),
        &observed,
        &truth,
        &|start| {
            let g = start.grid;
            let pml = 8;
            (0..g.n[2] - pml)
                .filter(|&k| k as f64 * g.h[2] >= 60.0)
                .flat_map(|k| {
                    (pml..g.n[1] - pml).flat_map(move |j| {
                        (pml..g.n[0] - pml).map(move |i| g.interior_idx(i, j, k))
                    })
                })
                .collect()
        },
        &[
            ("l2", "iterations = 8"),
            (
                "l2, 4 then 8 Hz then all",
                "[[inversion.stage]]\nmax_frequency = 4.0\niterations = 3\n\
                 [[inversion.stage]]\nmax_frequency = 8.0\niterations = 3\n\
                 [[inversion.stage]]\niterations = 2",
            ),
            (
                "time shifts, then l2",
                "misfit = \"traveltime\"\nmax_shift = 0.1\n[[inversion.stage]]\niterations = 4\n\
                 [[inversion.stage]]\nmisfit = \"l2\"\niterations = 4",
            ),
            (
                "phase, then l2",
                "misfit = \"phase\"\n[[inversion.stage]]\niterations = 4\n\
                 [[inversion.stage]]\nmisfit = \"l2\"\niterations = 4",
            ),
        ],
    );
}
