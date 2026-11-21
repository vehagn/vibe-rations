//! Inverting Marmousi2's elastic data (`examples/marmousi2/central.toml` from its start
//! smoothed over 250 m, against `central-true.toml`'s five shots): where the misfit lies, and
//! how the misfits, stages and preconditioners on main fare. `just marmousi2-study <name>`,
//! after `just marmousi2`. The research log records the results.
//!
//! Not part of the regular suite, except `the_split_adds_up`, which checks the splitting on a
//! small survey.
//!
//! - `misfit_by_arrival`: the start's balanced L2 misfit split by arrival time (before and
//!   after lines t = |Δx|/v + t₀), by component and by offset, and its gradient split at one
//!   of the lines: each part's norm and its cosine with the true update by depth band.
//! - `l2`, `frequency_stages`, `time_shifts_then_l2`, `phase_then_l2`: inversions from the
//!   same start, reporting after each stage the balanced L2 misfit of all the data (one scale
//!   for every variant, relative to the start's), and by depth band the error ratio
//!   ‖m - m_true‖ / ‖m₀ - m_true‖ and the cosine of the update m - m₀ with the true one.

use std::path::{Path, PathBuf};
use std::time::Instant;

use vibe_rations::config::Config;
use vibe_rations::gradient::{Strategy, balancing_weights, gradient_with, l2_misfit};
use vibe_rations::inversion::{Inversion, Observed};
use vibe_rations::model::Model;
use vibe_rations::receiver::Receivers;
use vibe_rations::simulation::Simulation;

/// A survey to study: the start (with `[inversion]`) and the truth, their files under `base`.
struct Case {
    start: Config,
    truth: Config,
    base: PathBuf,
    /// Depth bands (m) to report by, below the fixed part.
    bands: Vec<[f64; 2]>,
}

impl Case {
    /// `examples/marmousi2/central.toml` and `central-true.toml`, the start's `[inversion]`
    /// as `variant` makes it (`central_text`).
    fn central(variant: Variant) -> Case {
        let base = PathBuf::from("examples/marmousi2");
        let read = |name: &str| std::fs::read_to_string(base.join(name)).unwrap();
        let start = central_text(variant);
        let missing = base.join("25m-central/vp.npy");
        assert!(
            missing.is_file(),
            "{} is missing: run `just marmousi2` first",
            missing.display()
        );
        Case {
            start: Config::from_toml(&start).unwrap(),
            truth: Config::from_toml(&read("central-true.toml")).unwrap(),
            base,
            bands: vec![
                [500.0, 750.0],
                [750.0, 1000.0],
                [1000.0, 1500.0],
                [1500.0, 2500.0],
                [2500.0, 3125.0],
            ],
        }
    }

    fn model(&self, cfg: &Config) -> Model {
        Model::from_config(Simulation::grid(cfg).unwrap(), &cfg.model, &self.base).unwrap()
    }

    /// The nodes an inversion of the start may change: outside the absorbing layers and not
    /// above `fixed_above` (the same for vp and vs where the start is solid, as here).
    fn free(&self, start: &Model) -> Vec<usize> {
        let g = start.grid;
        let [nx, ny, nz] = g.n;
        let w = self.start.boundary.pml_width;
        let top = if self.start.boundary.free_surface {
            0
        } else {
            w
        };
        let fixed = self
            .start
            .inversion
            .as_ref()
            .and_then(|i| i.fixed_above)
            .unwrap_or(f64::NEG_INFINITY);
        let mut nodes = Vec::new();
        for k in top..nz - w {
            if (k as f64) * g.h[2] < fixed {
                continue;
            }
            for j in w..ny - w {
                for i in w..nx - w {
                    let n = g.interior_idx(i, j, k);
                    if start.vs[n] > 0.0 {
                        nodes.push(n);
                    }
                }
            }
        }
        nodes
    }

    /// The free nodes in each depth band.
    fn banded(&self, start: &Model) -> Vec<Vec<usize>> {
        let g = start.grid;
        let free = self.free(start);
        self.bands
            .iter()
            .map(|&[top, bottom]| {
                free.iter()
                    .copied()
                    .filter(|&n| {
                        let z = (n / (g.n[0] * g.n[1])) as f64 * g.h[2];
                        z >= top && z < bottom
                    })
                    .collect()
            })
            .collect()
    }

    /// The truth's traces of every shot, as observed data.
    fn record(&self) -> Vec<Observed> {
        (0..self.truth.sources.len())
            .map(|s| {
                let mut sim = Simulation::new(&self.truth.shot(s).unwrap(), &self.base).unwrap();
                sim.run(None, false).unwrap();
                sim.receivers.traces.into()
            })
            .collect()
    }
}

/// A variant of the start's `[inversion]`: keys in place of its `iterations = 5`, and
/// `[[inversion.stage]]` tables after the section, which ends the file.
type Variant = (&'static str, &'static str);

/// `examples/marmousi2/central.toml`'s text with `variant`'s inversion.
fn central_text((keys, tables): Variant) -> String {
    let text = std::fs::read_to_string("examples/marmousi2/central.toml").unwrap();
    assert!(text.contains("\niterations = 5\n"), "central.toml changed");
    let inversion = text
        .rfind("\n[inversion]\n")
        .expect("an [inversion] section");
    assert!(
        !text[inversion + 1..].contains("\n["),
        "[inversion] no longer ends central.toml"
    );
    let text = text.replace("\niterations = 5\n", &format!("\n{keys}\n"));
    format!("{}\n{tables}\n", text.trim_end())
}

const L2: Variant = ("iterations = 5", "");
const FREQUENCY_STAGES: Variant = (
    "",
    "[[inversion.stage]]\nmax_frequency = 0.6\niterations = 2\n\
     [[inversion.stage]]\nmax_frequency = 1.0\niterations = 2\n\
     [[inversion.stage]]\niterations = 1",
);
/// The stages, the best of the first four variants for vs, with the receivers' side of the
/// illumination estimated from eight runs at each stage's start.
const STAGES_WITH_RECEIVERS_SIDE: Variant = ("receiver_runs = 8", FREQUENCY_STAGES.1);
const TIME_SHIFTS_THEN_L2: Variant = (
    "max_shift = 0.8",
    "[[inversion.stage]]\nmisfit = \"traveltime\"\niterations = 3\n\
     [[inversion.stage]]\niterations = 2",
);
const PHASE_THEN_L2: Variant = (
    "",
    "[[inversion.stage]]\nmisfit = \"phase\"\niterations = 3\n\
     [[inversion.stage]]\niterations = 2",
);

/// Every variant's configuration parses and has the stages it means, each with five
/// iterations in all.
#[test]
fn the_variants_parse() {
    for (variant, stages) in [
        (L2, 1),
        (FREQUENCY_STAGES, 3),
        (STAGES_WITH_RECEIVERS_SIDE, 3),
        (TIME_SHIFTS_THEN_L2, 2),
        (PHASE_THEN_L2, 2),
        (("iterations = 0", ""), 1),
    ] {
        let cfg = Config::from_toml(&central_text(variant)).unwrap();
        let inversion = cfg.inversion.unwrap();
        assert_eq!(inversion.fixed_above, Some(500.0));
        let s = inversion.stages().unwrap();
        assert_eq!(s.len(), stages, "{variant:?}");
        let total: usize = s.iter().map(|s| s.iterations).sum();
        assert!(total == 5 || variant.0 == "iterations = 0", "{variant:?}");
    }
}

/// The error ratio ‖m - m_true‖ / ‖m₀ - m_true‖ and the cosine of m - m₀ with m_true - m₀
/// over `nodes`.
fn compare(m: &[f32], m0: &[f32], truth: &[f32], nodes: &[usize]) -> (f64, f64) {
    let (mut e, mut e0, mut dot, mut du, mut dt) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for &n in nodes {
        let (m, m0, t) = (m[n] as f64, m0[n] as f64, truth[n] as f64);
        e += (m - t).powi(2);
        e0 += (m0 - t).powi(2);
        dot += (m - m0) * (t - m0);
        du += (m - m0).powi(2);
        dt += (t - m0).powi(2);
    }
    (
        (e / e0).sqrt(),
        dot / (du * dt).sqrt().max(f64::MIN_POSITIVE),
    )
}

/// A table of `compare` for vp and vs by depth band.
fn print_bands(label: &str, case: &Case, m: &Model, m0: &Model, truth: &Model) {
    println!("{label}: by depth, error ratio and cosine with the true update");
    println!(
        "  {:>14}  {:>8} {:>7}  {:>8} {:>7}",
        "depth (m)", "vp", "cos", "vs", "cos"
    );
    let bands = case.banded(m0);
    let mut all = Vec::new();
    for (band, nodes) in case.bands.iter().zip(&bands) {
        all.extend_from_slice(nodes);
        let (ep, cp) = compare(&m.vp, &m0.vp, &truth.vp, nodes);
        let (es, cs) = compare(&m.vs, &m0.vs, &truth.vs, nodes);
        println!(
            "  {:>14}  {ep:8.3} {cp:+7.2}  {es:8.3} {cs:+7.2}",
            format!("{:.0}–{:.0}", band[0], band[1])
        );
    }
    let (ep, cp) = compare(&m.vp, &m0.vp, &truth.vp, &all);
    let (es, cs) = compare(&m.vs, &m0.vs, &truth.vs, &all);
    println!("  {:>14}  {ep:8.3} {cp:+7.2}  {es:8.3} {cs:+7.2}", "all");
}

/// Where a shot's balanced L2 misfit ½ Σ_c w_c Σ (d - d_obs)² lies, sample by sample.
#[derive(Clone, Debug, Default)]
struct Split {
    total: f64,
    /// Before each line, in the order the lines were given (the rest is after it).
    before: Vec<f64>,
    by_component: Vec<f64>,
    /// By offset bin: [edges[i], edges[i + 1]), the last open.
    by_offset: Vec<f64>,
}

impl Split {
    fn add(&mut self, other: &Split) {
        let sum = |a: &mut Vec<f64>, b: &[f64]| {
            a.resize(b.len().max(a.len()), 0.0);
            a.iter_mut().zip(b).for_each(|(a, b)| *a += b);
        };
        self.total += other.total;
        sum(&mut self.before, &other.before);
        sum(&mut self.by_component, &other.by_component);
        sum(&mut self.by_offset, &other.by_offset);
    }
}

/// Splits the balanced L2 misfit of `synthetic` against `observed` (one shot, its source at
/// `x_source`; receivers at `x`, samples at `times`) by the lines t = |Δx|/v + t₀ in `lines`,
/// by component and by the offset bins starting at `edges`.
fn split(
    synthetic: &Receivers,
    observed: &[Vec<f32>],
    x: &[f64],
    x_source: f64,
    times: &[f64],
    lines: &[(f64, f64)],
    edges: &[f64],
) -> Split {
    let weights = balancing_weights(observed, &synthetic.components);
    let samples = synthetic.samples;
    let mut s = Split {
        before: vec![0.0; lines.len()],
        by_component: vec![0.0; synthetic.components.len()],
        by_offset: vec![0.0; edges.len()],
        ..Split::default()
    };
    for (c, (syn, obs)) in synthetic.traces.iter().zip(observed).enumerate() {
        for (r, &xr) in x.iter().enumerate() {
            let offset = (xr - x_source).abs();
            let bin = edges.iter().rposition(|&e| offset >= e).unwrap_or(0);
            for (m, &t) in times.iter().enumerate() {
                let i = r * samples + m;
                let e = 0.5 * weights[c] * (syn[i] as f64 - obs[i] as f64).powi(2);
                s.total += e;
                s.by_component[c] += e;
                s.by_offset[bin] += e;
                for (b, &(v, t0)) in s.before.iter_mut().zip(lines) {
                    if t < offset / v + t0 {
                        *b += e;
                    }
                }
            }
        }
    }
    s
}

/// The receivers' x in metres.
fn receiver_x(sim: &Simulation) -> Vec<f64> {
    let h = sim.solver.grid.h[0];
    sim.receivers
        .nodes
        .iter()
        .map(|n| n[0] as f64 * h)
        .collect()
}

/// The start's misfit split as `split` does, summed over the shots, then its gradient split
/// at `line` (residuals before it, and the rest): each part's relative gradient (m₀ ∂χ/∂m at
/// the free nodes) by depth band, its norm and the cosine of its descent direction with the
/// true update. Returns the split, to check it.
fn split_by_arrival(
    case: &Case,
    lines: &[(f64, f64)],
    edges: &[f64],
    line: Option<(f64, f64)>,
) -> Split {
    let t = Instant::now();
    let observed = case.record();
    let start = case.model(&case.start);
    let truth = case.model(&case.truth);
    println!(
        "observed {} shot(s) in {:.0} s",
        observed.len(),
        t.elapsed().as_secs_f64()
    );
    let mut sum = Split::default();
    for (s, obs) in observed.iter().enumerate() {
        let shot = case.start.shot(s).unwrap();
        let mut sim = Simulation::with_model(&shot, start.clone()).unwrap();
        sim.run(None, false).unwrap();
        let x_source = shot.sources[0].position[0];
        let one = split(
            &sim.receivers,
            &obs.traces,
            &receiver_x(&sim),
            x_source,
            &sim.sample_times(),
            lines,
            edges,
        );
        let (l2, _) = l2_misfit(
            &sim.receivers.traces,
            &obs.traces,
            &balancing_weights(&obs.traces, &sim.receivers.components),
        )
        .unwrap();
        assert!(
            (one.total - l2).abs() <= 1e-9 * l2.max(1e-300),
            "{} against {l2}",
            one.total
        );
        println!("shot {s} at x = {x_source} m: misfit {:.4e}", one.total);
        sum.add(&one);
    }
    let pc = |v: f64| 100.0 * v / sum.total;
    println!("misfit {:.4e}, all shots", sum.total);
    for (&(v, t0), b) in lines.iter().zip(&sum.before) {
        println!(
            "  before t = |Δx|/{v} + {t0:.2}: {:5.1}%, after: {:5.1}%",
            pc(*b),
            pc(sum.total - b)
        );
    }
    for (c, e) in case
        .start
        .receivers
        .components
        .iter()
        .zip(&sum.by_component)
    {
        println!("  {}: {:5.1}%", c.name(), pc(*e));
    }
    for (i, e) in sum.by_offset.iter().enumerate() {
        let to = edges
            .get(i + 1)
            .map_or("∞".to_string(), |e| format!("{e:.0}"));
        println!("  offset {:.0}–{to} m: {:5.1}%", edges[i], pc(*e));
    }

    let Some((v, t0)) = line else {
        return sum;
    };
    // The gradient is linear in the residuals: the part before the line from residuals masked
    // to it, the part after as the rest.
    let t = Instant::now();
    let n = start.grid.interior_len();
    let mut parts = [
        [vec![0.0f64; n], vec![0.0; n]],
        [vec![0.0; n], vec![0.0; n]],
    ];
    for (s, obs) in observed.iter().enumerate() {
        let shot = case.start.shot(s).unwrap();
        let x_source = shot.sources[0].position[0];
        for (p, before_only) in [(0, true), (1, false)] {
            let mut sim = Simulation::with_model(&shot, start.clone()).unwrap();
            let x = receiver_x(&sim);
            let times = sim.sample_times();
            let samples = sim.receivers.samples;
            let run = gradient_with(
                &mut sim,
                |r| {
                    let w = balancing_weights(&obs.traces, &r.components);
                    let (misfit, mut residuals) = l2_misfit(&r.traces, &obs.traces, &w)?;
                    if before_only {
                        for res in &mut residuals {
                            for (rcv, &xr) in x.iter().enumerate() {
                                let offset = (xr - x_source).abs();
                                for (m, &t) in times.iter().enumerate() {
                                    if t >= offset / v + t0 {
                                        res[rcv * samples + m] = 0.0;
                                    }
                                }
                            }
                        }
                    }
                    Ok((misfit, residuals))
                },
                Strategy::default(),
            )
            .unwrap();
            let [gvp, gvs, _] = run.gradient.velocities(&start);
            for (sum, g) in parts[p].iter_mut().zip([gvp, gvs]) {
                sum.iter_mut().zip(g).for_each(|(a, b)| *a += b);
            }
        }
        println!(
            "gradients of shot {s} done, {:.0} s",
            t.elapsed().as_secs_f64()
        );
    }
    // After the line: all of it less the part before.
    let [before, total] = parts;
    let after: Vec<Vec<f64>> = total
        .iter()
        .zip(&before)
        .map(|(t, b)| t.iter().zip(b).map(|(t, b)| t - b).collect())
        .collect();
    let bands = case.banded(&start);
    println!(
        "gradient split at t = |Δx|/{v} + {t0:.2}: relative (m₀ ∂χ/∂m), its norm by depth band \
         and the cosine of -∇ with the true update"
    );
    for (name, values) in [
        ("before", &before[..]),
        ("after", &after[..]),
        ("all", &total[..]),
    ] {
        for (p, (g, m0, mt)) in [
            (&values[0], &start.vp, &truth.vp),
            (&values[1], &start.vs, &truth.vs),
        ]
        .into_iter()
        .enumerate()
        {
            let cells: Vec<String> = bands
                .iter()
                .map(|nodes| {
                    let (mut norm, mut dot, mut du) = (0.0, 0.0, 0.0);
                    for &n in nodes {
                        let rel = g[n] * m0[n] as f64;
                        let true_rel = mt[n] as f64 / m0[n] as f64 - 1.0;
                        norm += rel * rel;
                        dot += -rel * true_rel;
                        du += true_rel * true_rel;
                    }
                    format!(
                        "{:9.2e} {:+5.2}",
                        norm.sqrt(),
                        dot / (norm * du).sqrt().max(f64::MIN_POSITIVE)
                    )
                })
                .collect();
            println!("  {name:>6} {}: {}", ["vp", "vs"][p], cells.join("  "));
        }
    }
    println!(
        "  depth bands: {:?} m; {:.0} s",
        case.bands,
        t.elapsed().as_secs_f64()
    );
    sum
}

/// Inverts `case` from its start, printing after each stage the balanced L2 misfit of all
/// the data relative to the start's and the model's error and cosine by depth band, and
/// writing the model after each stage to `out/marmousi2/study/<name>`.
fn invert(name: &str, case: &Case) {
    let t = Instant::now();
    let observed = case.record();
    let start = case.model(&case.start);
    let truth = case.model(&case.truth);
    let settings = case.start.inversion.clone().unwrap();
    let mut inv = Inversion::new(&case.start, start.clone(), observed.clone(), &settings).unwrap();
    // The balanced L2 misfit of a model over every shot: one scale for every variant. The
    // shots as the inversion runs them (its time step, layers and operator).
    let shots: Vec<Config> = inv.shots.iter().map(|s| s.config.clone()).collect();
    let l2 = |m: &Model| -> f64 {
        observed
            .iter()
            .enumerate()
            .map(|(s, obs)| {
                let mut sim = Simulation::with_model(&shots[s], m.clone()).unwrap();
                sim.run(None, false).unwrap();
                let w = balancing_weights(&obs.traces, &sim.receivers.components);
                l2_misfit(&sim.receivers.traces, &obs.traces, &w).unwrap().0
            })
            .sum()
    };
    let first = l2(&start);
    println!(
        "{name}: {} unknowns; L2 {first:.4e}; set up in {:.0} s",
        inv.unknowns(),
        t.elapsed().as_secs_f64()
    );
    let out = PathBuf::from("out/marmousi2/study").join(name);
    let stages = inv.stages().to_vec();
    let mut iterations = 0;
    for (s, stage) in stages.iter().enumerate() {
        if s > 0 {
            inv.begin_stage(s).unwrap();
        }
        for _ in 0..stage.iterations {
            match inv.iterate().unwrap() {
                Some(it) => {
                    iterations += 1;
                    println!(
                        "{name}, stage {s}, iteration {}: misfit {:.4e} as the stage measures \
                         it, step {:.2} after {} trial(s), largest change {:.1}%, {:.0} s",
                        it.iteration,
                        it.misfit,
                        it.step,
                        it.trials,
                        100.0 * it.largest_change,
                        t.elapsed().as_secs_f64()
                    );
                }
                None => {
                    println!("{name}: no step lowers stage {s}'s misfit");
                    break;
                }
            }
        }
        println!(
            "{name}: after stage {s} ({}{}{}; {iterations} iterations in all): L2 {:.1}% of \
             the start's",
            stage.misfit.unwrap_or_default().name(),
            stage
                .max_frequency
                .map_or(String::new(), |f| format!(" below {f} Hz")),
            if settings.receiver_runs > 0 {
                format!(", receivers' side from {} runs", settings.receiver_runs)
            } else {
                String::new()
            },
            100.0 * l2(&inv.model) / first
        );
        print_bands(name, case, &inv.model, &start, &truth);
        write_model(&out.join(format!("stage-{s}")), &inv.model);
    }
    println!("{name}: {:.0} s in all", t.elapsed().as_secs_f64());
}

fn write_model(dir: &Path, m: &Model) {
    std::fs::create_dir_all(dir).unwrap();
    let shape = [m.grid.n[2], m.grid.n[1], m.grid.n[0]];
    for (name, values) in [("vp", &m.vp), ("vs", &m.vs), ("rho", &m.rho)] {
        vibe_rations::npy::write(&dir.join(format!("{name}.npy")), &shape, values).unwrap();
    }
}

/// The lines the misfit is split at: S speeds of the seafloor's sediments, from about 300
/// m/s in the true model to 400 in the smoothed one, and one a little faster; t₀ is when a
/// slow wave can first leave the seafloor below the source: the Ricker's peak (1.2/f₀), less
/// half a period, plus the water's 425 m at 1500 m/s.
fn central_lines() -> Vec<(f64, f64)> {
    let f0 = 1.2;
    let t0 = 1.2 / f0 - 0.5 / f0 + 425.0 / 1500.0;
    [300.0, 400.0, 500.0].map(|v| (v, t0)).to_vec()
}

const CENTRAL_OFFSETS: [f64; 4] = [0.0, 1000.0, 2000.0, 3000.0];

#[test]
#[ignore = "a study: `just marmousi2-study misfit_by_arrival` (about 20 minutes)"]
fn misfit_by_arrival() {
    let case = Case::central(("iterations = 0", ""));
    let lines = central_lines();
    split_by_arrival(&case, &lines, &CENTRAL_OFFSETS, Some(lines[1]));
}

#[test]
#[ignore = "a study: `just marmousi2-study l2` (about 40 minutes)"]
fn l2() {
    invert("l2", &Case::central(L2));
}

#[test]
#[ignore = "a study: `just marmousi2-study frequency_stages` (about 50 minutes)"]
fn frequency_stages() {
    invert("frequency-stages", &Case::central(FREQUENCY_STAGES));
}

#[test]
#[ignore = "a study: `just marmousi2-study stages_with_receivers_side` (about an hour)"]
fn stages_with_receivers_side() {
    invert(
        "stages-with-receivers-side",
        &Case::central(STAGES_WITH_RECEIVERS_SIDE),
    );
}

#[test]
#[ignore = "a study: `just marmousi2-study time_shifts_then_l2` (about 45 minutes)"]
fn time_shifts_then_l2() {
    invert("time-shifts-then-l2", &Case::central(TIME_SHIFTS_THEN_L2));
}

#[test]
#[ignore = "a study: `just marmousi2-study phase_then_l2` (about 45 minutes)"]
fn phase_then_l2() {
    invert("phase-then-l2", &Case::central(PHASE_THEN_L2));
}

/// The split on a small two-layer survey adds up: its parts make the misfit `l2_misfit`
/// gives, the gradient's parts the whole (the study prints them), and the cosines are
/// finite.
#[test]
fn the_split_adds_up() {
    let dir = std::env::temp_dir().join(format!("marmousi2-study-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let survey = |vs_top: f64, inversion: &str| {
        format!(
            r#"
            [grid]
            points = [40, 24, 26]
            spacing = 10.0
            [time]
            duration = 0.25
            dt = 5e-4
            [boundary]
            pml_width = 6
            pml_vp_max = 2600.0
            [scheme]
            coefficients = "taylor"
            [[model.layer]]
            top = 0.0
            vp = 2000.0
            vs = {vs_top}
            rho = 2000.0
            [[model.layer]]
            top = 100.0
            vp = 2600.0
            vs = 1400.0
            rho = 2200.0
            [[source]]
            position = [120.0, 120.0, 20.0]
            explosion = 1e12
            wavelet = {{ kind = "ricker", frequency = 20.0 }}
            [[source]]
            position = [270.0, 120.0, 20.0]
            explosion = 1e12
            wavelet = {{ kind = "ricker", frequency = 20.0 }}
            [receivers]
            components = ["vx", "vz"]
            [[receivers.line]]
            from = [100.0, 120.0, 40.0]
            to = [290.0, 120.0, 40.0]
            count = 20
            {inversion}
            "#
        )
    };
    let case = Case {
        start: Config::from_toml(&survey(
            900.0,
            "[inversion]\nobserved = \"unused\"\nfixed_above = 40.0",
        ))
        .unwrap(),
        truth: Config::from_toml(&survey(1000.0, "")).unwrap(),
        base: dir.clone(),
        bands: vec![[40.0, 100.0], [100.0, 200.0]],
    };
    let lines = [(800.0, 0.05), (1500.0, 0.05)];
    let edges = [0.0, 50.0, 100.0];
    let s = split_by_arrival(&case, &lines, &edges, Some(lines[0]));
    let close = |a: f64| (a - s.total).abs() <= 1e-9 * s.total;
    assert!(s.total > 0.0);
    assert!(close(s.by_component.iter().sum()));
    assert!(close(s.by_offset.iter().sum()));
    // A faster line comes earlier: less before it.
    assert!(s.before[1] <= s.before[0] && s.before[0] <= s.total);
    std::fs::remove_dir_all(dir).unwrap();
}
