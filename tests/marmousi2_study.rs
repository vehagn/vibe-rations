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
//! - `misfit_by_mute`: what the mutes of the muted and down-weighted variants keep of that
//!   misfit and of its gradient, by component, by depth band.
//! - `l2_checkpointed`: `l2` with the forward wavefields replayed from checkpoints instead
//!   of rebuilt backwards (the default since 2026-11-24; before, this was `l2_rebuilt`, the
//!   other way round), and how far its model lies from `l2`'s.
//! - `mute_lines`, `west_l2`, `west_muted_l2`, `west_stages_muted_after_the_first`,
//!   `west_misfit_by_mute`: the western 6 km (`Case::part`), held out from the choices made
//!   on the central part, with the slow waves' mute worked out from its own start
//!   (`slow_waves_mute`).
//! - `phase_check`: whether the start predicts the slow waves within half a cycle, from the
//!   data: the phase of the observed traces against the start's after the slow waves' line,
//!   unwrapped along offset shot by shot and side by side, at the stages' frequencies, next
//!   to what the truth predicts. `central_phase_checked_l2`, `west_phase_checked_l2`,
//!   `central_phase_checked_stages`, `west_phase_checked_stages`: inversions muting the slow
//!   waves only beyond the offsets it passes (`phase_mutes`).
//! - `whole_by_region`: the whole-section inversion's result (`examples/marmousi2/whole.toml`)
//!   in 1 km windows along x, next to how far its start's seafloor vs is from the truth's.
//! - `l2`, `frequency_stages`, `time_shifts_then_l2`, `phase_then_l2`,
//!   `stages_with_receivers_side`, and `muted_l2`, `muted_stages`, `down_weighted_l2` and
//!   `down_weighted_stages` with the seafloor's slow waves muted or down-weighted: inversions
//!   from the same start, reporting after each stage the balanced L2 misfit of all the data (one scale
//!   for every variant, relative to the start's), and by depth band the error ratio
//!   ‖m - m_true‖ / ‖m₀ - m_true‖ and the cosine of the update m - m₀ with the true one.

use std::path::{Path, PathBuf};
use std::time::Instant;

use vibe_rations::config::{Config, MuteConfig};
use vibe_rations::gradient::{Strategy, balancing_weights, gradient_with, l2_misfit};
use vibe_rations::inversion::{Inversion, Observed};
use vibe_rations::misfit::{mute_weight, offsets};
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
        Case::part("central", variant)
    }

    /// [`Case::central`] on another 6 km of the section, whose model files `convert.sh` makes
    /// as `25m-<part>` and `25m-<part>-smooth`: the same survey, start and numerics, another
    /// model below the same seafloor (`west`, x = 0–6 km, is held out from the choices made
    /// on the central part).
    fn part(part: &str, variant: Variant) -> Case {
        let base = PathBuf::from("examples/marmousi2");
        let on_part = |text: String| text.replace("25m-central", &format!("25m-{part}"));
        let start = on_part(central_text(variant));
        let truth = on_part(std::fs::read_to_string(base.join("central-true.toml")).unwrap());
        let missing = base.join(format!("25m-{part}/vp.npy"));
        assert!(
            missing.is_file(),
            "{} is missing: run `just marmousi2` first",
            missing.display()
        );
        Case {
            start: Config::from_toml(&start).unwrap(),
            truth: Config::from_toml(&truth).unwrap(),
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
    /// above `fixed_above` (the same for vp and vs where the start is solid, as here), within
    /// `ACROSS` nodes of the receivers' plane.
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
        let line = g
            .nearest_node(self.start.receivers.positions().unwrap()[0])
            .unwrap()[1];
        let mut nodes = Vec::new();
        for k in top..nz - w {
            if (k as f64) * g.h[2] < fixed {
                continue;
            }
            for j in (w..ny - w).filter(|j| j.abs_diff(line) <= ACROSS) {
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

/// How far from the receivers' plane, in nodes, the error ratios and cosines are measured.
/// The parts are 41 nodes across with 16 in each side layer, which leaves these 9 nodes; the
/// whole section (`whole_by_region`) is 201 across, and over all of it the 4 km away from the
/// line, which the survey barely sees, would dilute every number.
const ACROSS: usize = 4;

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
/// What `misfit_by_arrival` found the slow seafloor waves arrive after: a line at 600 m/s from
/// 0.87 s (`central_lines`' t₀), over half a period (0.42 s) of taper, so that the weight is down
/// by the time the start's 400 m/s waves arrive (at 0.5 km offset just then, at 1 km 0.4 s
/// before). The inversion's mutes are the same in every stage.
const SLOW_WAVES_MUTED: &str = "[[inversion.mute]]\nkeep = \"before\"\nvelocity = 600.0\n\
                                delay = 0.87\ntaper = 0.42";
/// The same line, down-weighting to 0.1 (1% in the L2 misfit) instead of removing.
const SLOW_WAVES_DOWN_WEIGHTED: &str = "[[inversion.mute]]\nkeep = \"before\"\n\
                                        velocity = 600.0\ndelay = 0.87\ntaper = 0.42\n\
                                        weight = 0.1";
/// L2 with the slow waves muted, and with them down-weighted; the same with the stages.
fn muted(stages: bool, mute: &'static str) -> Variant {
    muted_with(if stages { "" } else { "iterations = 5" }, stages, mute)
}

/// [`muted`] with other keys in `[inversion]`.
fn muted_with(keys: &'static str, stages: bool, mute: &str) -> Variant {
    let tables = if stages {
        format!("{mute}\n{}", FREQUENCY_STAGES.1)
    } else {
        mute.to_string()
    };
    (keys, leak(tables))
}

/// `Variant` holds `&'static str`s; the tables put together at a test's start live on.
fn leak(s: String) -> &'static str {
    Box::leak(s.into_boxed_str())
}

/// Muted L2 with the receivers' side of the illumination from eight runs: set aside on the
/// stages (`stages_with_receivers_side`) because it moved the deep model further the wrong
/// way, which the mute stopped.
fn muted_l2_with_receivers_side_variant() -> Variant {
    muted_with("iterations = 5\nreceiver_runs = 8", false, SLOW_WAVES_MUTED)
}

/// Muted L2 for twenty iterations, the model reported every five: does it keep approaching
/// the truth, or stall?
fn muted_l2_longer_variant() -> Variant {
    muted_with("iterations = 20", false, SLOW_WAVES_MUTED)
}

/// The stages with the slow waves muted from the second stage on (`[[inversion.stage.mute]]`):
/// at 0.6 Hz, unmuted, they turned vs towards the truth down to 1 km (`frequency_stages`), and
/// muted they didn't (`muted_stages`).
fn stages_muted_after_the_first_variant() -> Variant {
    stages_muted_after_the_first_with(SLOW_WAVES_MUTED)
}

/// [`stages_muted_after_the_first_variant`] with another mute.
fn stages_muted_after_the_first_with(mute: &str) -> Variant {
    let mute = mute.replacen("[[inversion.mute]]", "[[inversion.stage.mute]]", 1);
    (
        "",
        leak(format!(
            "[[inversion.stage]]\nmax_frequency = 0.6\niterations = 2\n\
             [[inversion.stage]]\nmax_frequency = 1.0\niterations = 2\n{mute}\n\
             [[inversion.stage]]\niterations = 1\n{mute}"
        )),
    )
}

/// The leading edge of a window around the early arrivals: kept after a line ahead of the
/// first arrival, faster than any wave in these models (4.7 km/s) from the Ricker's onset (0.2
/// s). On these noise-free traces it removes nothing; with the slow waves' mute it makes the
/// window the time shifts measure in.
const BEFORE_THE_FIRST_ARRIVALS: &str = "[[inversion.mute]]\nkeep = \"after\"\n\
                                         velocity = 4700.0\ndelay = 0.2\ntaper = 0.42";

/// Time shifts in that window for three iterations, then L2 in it for two.
fn windowed_time_shifts_then_l2_variant() -> Variant {
    (
        "max_shift = 0.8",
        leak(format!(
            "{SLOW_WAVES_MUTED}\n{BEFORE_THE_FIRST_ARRIVALS}\n{}",
            TIME_SHIFTS_THEN_L2.1
        )),
    )
}

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

/// L2 with each shot's forward wavefield rebuilt backwards from the band around the layers
/// instead of replayed from checkpoints.
const L2_CHECKPOINTED: Variant = ("iterations = 5\nreconstruct = false", "");

/// Every variant's configuration parses and has the stages, mutes (the inversion's and each
/// stage's) and iterations it means.
#[test]
fn the_variants_parse() {
    for (variant, stages, mutes, stage_mutes, iterations) in [
        (L2, 1, 0, 0, 5),
        (L2_CHECKPOINTED, 1, 0, 0, 5),
        (FREQUENCY_STAGES, 3, 0, 0, 5),
        (STAGES_WITH_RECEIVERS_SIDE, 3, 0, 0, 5),
        (muted(false, SLOW_WAVES_MUTED), 1, 1, 0, 5),
        (muted(true, SLOW_WAVES_MUTED), 3, 1, 0, 5),
        (muted(false, SLOW_WAVES_DOWN_WEIGHTED), 1, 1, 0, 5),
        (muted(true, SLOW_WAVES_DOWN_WEIGHTED), 3, 1, 0, 5),
        (muted_l2_with_receivers_side_variant(), 1, 1, 0, 5),
        (muted_l2_longer_variant(), 1, 1, 0, 20),
        (stages_muted_after_the_first_variant(), 3, 0, 2, 5),
        (windowed_time_shifts_then_l2_variant(), 2, 2, 0, 5),
        (TIME_SHIFTS_THEN_L2, 2, 0, 0, 5),
        (PHASE_THEN_L2, 2, 0, 0, 5),
        (("iterations = 0", ""), 1, 0, 0, 0),
    ] {
        let cfg = Config::from_toml(&central_text(variant)).unwrap();
        let inversion = cfg.inversion.unwrap();
        assert_eq!(inversion.fixed_above, Some(500.0));
        let s = inversion.stages().unwrap();
        assert_eq!(s.len(), stages, "{variant:?}");
        assert_eq!(inversion.mutes.len(), mutes, "{variant:?}");
        let muted_stages = s.iter().filter(|s| !s.mutes.is_empty()).count();
        assert_eq!(muted_stages, stage_mutes, "{variant:?}");
        let total: usize = s.iter().map(|s| s.iterations).sum();
        assert_eq!(total, iterations, "{variant:?}");
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
    sim.receivers.points.iter().map(|p| p[0] * h).collect()
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

/// The slow waves' mute for `case`, worked out from its start and source alone, not from the
/// true model: a line through t₀, the earliest a wave can leave the seafloor below the
/// source (the Ricker's peak 1.2/f₀ less half a period, plus the `water` metres between the
/// sources and the receivers at 1500 m/s), as fast as needs be for a taper of half a period
/// (0.5/f₀) to end, from 500 m offset on, before waves at the start's slowest S speed at the
/// receivers v_s: 1/v_line = 1/v_s - (0.5/f₀)/500 m. The observed slow waves must come
/// later too, which the true models here satisfy (345 m/s at the seafloor) but which these
/// runs didn't check from the data: stacks along lines and f-k spectra of the gathers within
/// 1 km didn't find them under the direct wave's energy. Returns the line's speed and the
/// mute.
fn slow_waves_mute(case: &Case, water: f64) -> (f64, String) {
    let start = case.model(&case.start);
    let grid = Simulation::grid(&case.start).unwrap();
    let slowest = case
        .start
        .receivers
        .positions()
        .unwrap()
        .into_iter()
        .map(|p| {
            let [i, j, k] = grid.nearest_node(p).unwrap();
            start.vs[grid.interior_idx(i, j, k)] as f64
        })
        .fold(f64::INFINITY, f64::min);
    let f0 = 1.2;
    let taper = 0.5 / f0;
    let velocity = 1.0 / (1.0 / slowest - taper / 500.0);
    let delay = 1.2 / f0 - taper + water / 1500.0;
    (
        velocity,
        format!(
            "[[inversion.mute]]\nkeep = \"before\"\nvelocity = {velocity:.1}\n\
             delay = {delay:.2}\ntaper = {taper:.2}"
        ),
    )
}

/// The product of `mutes`' weights on one shot's traces at `offsets`, receiver-major like the
/// traces, at t = m Δ from the run's start, as the inversion's misfits apply them.
fn mute_weights(mutes: &[MuteConfig], offsets: &[f64], samples: usize, interval: f64) -> Vec<f64> {
    let mut weights = Vec::with_capacity(offsets.len() * samples);
    for &x in offsets {
        for m in 0..samples {
            let t = m as f64 * interval;
            weights.push(mutes.iter().map(|mute| mute_weight(mute, t, x)).product());
        }
    }
    weights
}

/// For each depth band (`bands`, free nodes), the norm of the relative gradient m₀ ∂χ/∂m and
/// the cosines of its descent direction with the true update m_true / m₀ - 1, as it is and
/// multiplied by the preconditioner `p` (one value per interior node).
fn gradient_cells(g: &[f64], m0: &[f32], truth: &[f32], p: &[f64], bands: &[Vec<usize>]) -> String {
    let cos = |dot: f64, a: f64, b: f64| dot / (a * b).sqrt().max(f64::MIN_POSITIVE);
    bands
        .iter()
        .map(|nodes| {
            let (mut norm, mut dot, mut norm_p, mut dot_p, mut du) = (0.0, 0.0, 0.0, 0.0, 0.0);
            for &n in nodes {
                let rel = g[n] * m0[n] as f64;
                let true_rel = truth[n] as f64 / m0[n] as f64 - 1.0;
                norm += rel * rel;
                dot -= rel * true_rel;
                norm_p += (p[n] * rel).powi(2);
                dot_p -= p[n] * rel * true_rel;
                du += true_rel * true_rel;
            }
            format!(
                "{:9.2e} {:+5.2} {:+5.2}",
                norm.sqrt(),
                cos(dot, norm, du),
                cos(dot_p, norm_p, du)
            )
        })
        .collect::<Vec<_>>()
        .join("  ")
}

/// What each set of mutes in `sets` (named; a set's weights M multiply) keeps of the start's
/// balanced L2 misfit and of its gradient. The misfit kept is ½ Σ_c w_c Σ (M r)², r = d -
/// d_obs, with the balance w_c of the traces as recorded, so that the shares compare with
/// `misfit_by_arrival`'s; the inversion weighs by w'_c from the muted observed traces instead,
/// and that measure is printed too. Its gradient (residuals w_c M² r) for each component
/// alone, both together, both with the inversion's balance, and the part the mutes remove
/// (the whole less the part kept): by depth band, the norm of m₀ ∂χ/∂m and the cosines of its
/// descent direction with the true update, raw and divided by the source-side illumination
/// with the inversion's water level, as its preconditioner divides. A set without mutes gives
/// the whole misfit and gradient, by component; it must come first.
fn split_by_mute(case: &Case, sets: &[(&str, Vec<MuteConfig>)]) {
    assert!(sets[0].1.is_empty(), "the first set keeps everything");
    let t = Instant::now();
    let observed = case.record();
    let start = case.model(&case.start);
    let truth = case.model(&case.truth);
    println!(
        "observed {} shot(s) in {:.0} s",
        observed.len(),
        t.elapsed().as_secs_f64()
    );
    let n = start.grid.interior_len();
    let names: Vec<&str> = case
        .start
        .receivers
        .components
        .iter()
        .map(|c| c.name())
        .collect();
    let zero = || [vec![0.0f64; n], vec![0.0; n]];
    // Each set's kept gradient (vp, vs) by component and with the inversion's balance, and
    // its kept misfit by component and as the inversion measures it.
    let mut kept: Vec<Vec<[Vec<f64>; 2]>> = sets
        .iter()
        .map(|_| names.iter().map(|_| zero()).collect())
        .collect();
    let mut own: Vec<[Vec<f64>; 2]> = sets.iter().map(|_| zero()).collect();
    let mut kept_misfit = vec![vec![0.0f64; names.len()]; sets.len()];
    let mut own_misfit = vec![0.0f64; sets.len()];
    let mut illumination = vec![0.0f64; n];
    let add = |sum: &mut [Vec<f64>; 2], g: &[Vec<f64>], factor: f64| {
        for (s, g) in sum.iter_mut().zip(g) {
            s.iter_mut().zip(g).for_each(|(a, b)| *a += factor * b);
        }
    };
    for (s, obs) in observed.iter().enumerate() {
        let shot = case.start.shot(s).unwrap();
        let x = offsets(&shot).unwrap();
        let components = &shot.receivers.components;
        let w = balancing_weights(&obs.traces, components);
        for (k, (_, mutes)) in sets.iter().enumerate() {
            for c in 0..names.len() {
                let mut sim = Simulation::with_model(&shot, start.clone()).unwrap();
                let m = mute_weights(mutes, &x, sim.receivers.samples, sim.info.sample_interval);
                let muted: Vec<Vec<f32>> = obs
                    .traces
                    .iter()
                    .map(|o| {
                        o.iter()
                            .zip(&m)
                            .map(|(&v, &m)| (v as f64 * m) as f32)
                            .collect()
                    })
                    .collect();
                let own_w = balancing_weights(&muted, components);
                let mut energy = 0.0;
                let run = gradient_with(
                    &mut sim,
                    |r| {
                        // Component c's kept part alone.
                        let residuals = r
                            .traces
                            .iter()
                            .enumerate()
                            .map(|(d, syn)| {
                                if d != c {
                                    return vec![0.0; syn.len()];
                                }
                                syn.iter()
                                    .zip(&obs.traces[c])
                                    .zip(&m)
                                    .map(|((&a, &b), &m)| {
                                        let r = m * (a as f64 - b as f64);
                                        energy += r * r;
                                        (w[c] * m * r) as f32
                                    })
                                    .collect()
                            })
                            .collect();
                        Ok((0.5 * w[c] * energy, residuals))
                    },
                    Strategy::default(),
                )
                .unwrap();
                kept_misfit[k][c] += run.misfit;
                own_misfit[k] += 0.5 * own_w[c] * energy;
                let g = run.gradient.velocities(&start);
                add(&mut kept[k][c], &g, 1.0);
                if w[c] > 0.0 {
                    add(&mut own[k], &g, own_w[c] / w[c]);
                }
                if k == 0 && c == 0 {
                    illumination
                        .iter_mut()
                        .zip(&run.illumination)
                        .for_each(|(a, b)| *a += b);
                }
            }
        }
        println!(
            "shot {s} at x = {} m done, {:.0} s",
            shot.sources[0].position[0],
            t.elapsed().as_secs_f64()
        );
    }

    // The inversion's preconditioner: 1 / (I / I_max + water level), I_max over the free nodes.
    let bands = case.banded(&start);
    let free = case.free(&start);
    let peak = free.iter().map(|&n| illumination[n]).fold(0.0, f64::max);
    let level = case.start.inversion.as_ref().unwrap().water_level;
    let p: Vec<f64> = illumination
        .iter()
        .map(|&i| 1.0 / (i / peak + level))
        .collect();
    let sum = |parts: &[[Vec<f64>; 2]]| -> [Vec<f64>; 2] {
        let mut total = zero();
        for part in parts {
            add(&mut total, part, 1.0);
        }
        total
    };
    let whole = sum(&kept[0]);
    let all: f64 = kept_misfit[0].iter().sum();
    println!(
        "misfit {all:.4e} in all; by depth band {:?} m, the norm of m₀ ∂χ/∂m, the cosine of -∇ \
         with the true update, and the same divided by the illumination",
        case.bands
    );
    let row = |label: &str, g: &[Vec<f64>; 2]| {
        for (p_name, (g, m0, mt)) in [
            ("vp", (&g[0], &start.vp, &truth.vp)),
            ("vs", (&g[1], &start.vs, &truth.vs)),
        ] {
            println!(
                "  {label:>22} {p_name}: {}",
                gradient_cells(g, m0, mt, &p, &bands)
            );
        }
    };
    for (k, (name, _)) in sets.iter().enumerate() {
        let kept_all: f64 = kept_misfit[k].iter().sum();
        let by_component: Vec<String> = names
            .iter()
            .zip(&kept_misfit[k])
            .map(|(c, e)| format!("{c} {:.1}%", 100.0 * e / all))
            .collect();
        println!(
            "{name}: keeps {:.1}% of the misfit ({}); as the inversion measures it, {:.4e}",
            100.0 * kept_all / all,
            by_component.join(", "),
            own_misfit[k]
        );
        let together = sum(&kept[k]);
        row(&format!("{name}, kept"), &together);
        for (c, g) in names.iter().zip(&kept[k]) {
            row(&format!("{name}, kept {c}"), g);
        }
        if k > 0 {
            row(&format!("{name}, inversion's balance"), &own[k]);
            let removed = [0, 1].map(|i| {
                whole[i]
                    .iter()
                    .zip(&together[i])
                    .map(|(a, b)| a - b)
                    .collect::<Vec<f64>>()
            });
            row(&format!("{name}, removed"), &removed);
        }
    }
    println!("{:.0} s in all", t.elapsed().as_secs_f64());
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
        for i in 1..=stage.iterations {
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
                    // In a long stage, how it goes every five iterations.
                    if i % 5 == 0 && i < stage.iterations {
                        println!(
                            "{name}: after iteration {iterations}: L2 {:.1}% of the start's",
                            100.0 * l2(&inv.model) / first
                        );
                        print_bands(name, case, &inv.model, &start, &truth);
                    }
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

/// How far variant `b`'s model after `stage` lies from variant `a`'s, both as `invert` wrote
/// them: by depth band, ‖m_b - m_a‖ / ‖m_a - m₀‖ (against a's update) and the largest
/// |m_b - m_a| / m_a, for vp and vs.
fn compare_variants(case: &Case, a: &str, b: &str, stage: usize) {
    let path = |name: &str, p: &str| {
        PathBuf::from("out/marmousi2/study")
            .join(name)
            .join(format!("stage-{stage}/{p}.npy"))
    };
    if !path(a, "vp").is_file() {
        println!(
            "no {} to compare with: run `{a}` first",
            path(a, "vp").display()
        );
        return;
    }
    let start = case.model(&case.start);
    let bands = case.banded(&start);
    println!(
        "{b} against {a} after stage {stage}: by depth band {:?} m, ‖m_b - m_a‖ / ‖m_a - m₀‖ and \
         the largest |m_b - m_a| / m_a",
        case.bands
    );
    for (p, m0) in [("vp", &start.vp), ("vs", &start.vs)] {
        let read = |name: &str| vibe_rations::npy::read(&path(name, p)).unwrap().1;
        let (ma, mb) = (read(a), read(b));
        let cells: Vec<String> = bands
            .iter()
            .map(|nodes| {
                let (mut diff, mut update, mut largest) = (0.0f64, 0.0f64, 0.0f64);
                for &n in nodes {
                    let (a, b, m0) = (ma[n] as f64, mb[n] as f64, m0[n] as f64);
                    diff += (b - a).powi(2);
                    update += (a - m0).powi(2);
                    largest = largest.max(((b - a) / a).abs());
                }
                format!("{:9.2e} {largest:9.2e}", (diff / update).sqrt())
            })
            .collect();
        println!("  {p}: {}", cells.join("  "));
    }
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

/// The mutes of a variant's `[inversion]`.
fn mutes_of(variant: Variant) -> Vec<MuteConfig> {
    Config::from_toml(&central_text(variant))
        .unwrap()
        .inversion
        .unwrap()
        .mutes
}

#[test]
#[ignore = "a study: `just marmousi2-study misfit_by_mute` (about 35 minutes)"]
fn misfit_by_mute() {
    split_by_mute(
        &Case::central(("iterations = 0", "")),
        &[
            ("all", Vec::new()),
            ("muted", mutes_of(muted(false, SLOW_WAVES_MUTED))),
            (
                "down-weighted",
                mutes_of(muted(false, SLOW_WAVES_DOWN_WEIGHTED)),
            ),
        ],
    );
}

#[test]
#[ignore = "a study: `just marmousi2-study west_stages_muted_after_the_first` (about an hour)"]
fn west_stages_muted_after_the_first() {
    let (velocity, mute) = slow_waves_mute(&Case::part("west", L2), 425.0);
    println!("the mute line: {velocity:.0} m/s");
    invert(
        "west-stages-muted-after-the-first",
        &Case::part("west", stages_muted_after_the_first_with(&mute)),
    );
}

/// `misfit_by_mute` on the western part, with the mute `slow_waves_mute` works out for it:
/// on the central part the gradient's removed part pointed away from the truth below 750 m,
/// and there the mute helped; here L2 unmuted moved vs towards the truth.
#[test]
#[ignore = "a study: `just marmousi2-study west_misfit_by_mute` (about 20 minutes)"]
fn west_misfit_by_mute() {
    let (velocity, mute) = slow_waves_mute(&Case::part("west", L2), 425.0);
    println!("the mute line: {velocity:.0} m/s");
    split_by_mute(
        &Case::part("west", ("iterations = 0", "")),
        &[
            ("all", Vec::new()),
            (
                "muted",
                mutes_of(muted_with("iterations = 0", false, &mute)),
            ),
        ],
    );
}

/// The start's S speeds under the seafloor against the truth's, for each part: what decides
/// whether the start predicts the slow waves within half a cycle. Over the columns inside the
/// side layers, on the receivers' plane: vs on the seafloor's first solid node (the
/// receivers'), the harmonic mean of vs over the 250 m below it (ten nodes), and the delay a
/// wave at the truth's speeds gathers over 1 km behind one at the start's, 1 km times the
/// difference of those mean slownesses, averaged over the columns and at its largest. Uses
/// the truth: a diagnosis, not a rule an inversion could follow.
#[test]
#[ignore = "a study: `just marmousi2-study seafloor_speeds` (seconds)"]
fn seafloor_speeds() {
    for part in ["central", "west"] {
        let case = Case::part(part, L2);
        let (start, truth) = (case.model(&case.start), case.model(&case.truth));
        let g = start.grid;
        let receiver = case.start.receivers.positions().unwrap()[0];
        let [_, j, top] = g.nearest_node(receiver).unwrap();
        let w = case.start.boundary.pml_width;
        let columns: Vec<usize> = (w..g.n[0] - w).collect();
        let vs = |m: &Model, i: usize, k: usize| m.vs[g.interior_idx(i, j, k)] as f64;
        let slowness =
            |m: &Model, i: usize| (top..top + 10).map(|k| 1.0 / vs(m, i, k)).sum::<f64>() / 10.0;
        let range = |m: &Model| {
            columns
                .iter()
                .map(|&i| vs(m, i, top))
                .fold((f64::INFINITY, 0.0f64), |(lo, hi), v| {
                    (lo.min(v), hi.max(v))
                })
        };
        let mean = |f: &dyn Fn(usize) -> f64| {
            columns.iter().map(|&i| f(i)).sum::<f64>() / columns.len() as f64
        };
        let delay = |i: usize| 1000.0 * (slowness(&truth, i) - slowness(&start, i));
        let largest = columns.iter().map(|&i| delay(i)).fold(0.0f64, f64::max);
        let ((tl, th), (sl, sh)) = (range(&truth), range(&start));
        println!(
            "{part}: vs at {:.0} m, true {tl:.0}–{th:.0} m/s, start {sl:.0}–{sh:.0}; over the \
             250 m below, harmonic means true {:.0}, start {:.0} m/s; delay over 1 km {:.2} s \
             on average, {largest:.2} s at most",
            top as f64 * g.h[2],
            1.0 / mean(&|i| slowness(&truth, i)),
            1.0 / mean(&|i| slowness(&start, i)),
            mean(&delay),
        );
    }
}

/// The whole-section inversion (`examples/marmousi2/whole.toml`, run with `invert`, its
/// model in `out/marmousi2/whole`) region by region: whether its start predicts the slow
/// waves within half a cycle varies along the section, as between the central and western
/// parts (`seafloor_speeds`). In 1 km windows along x, on the receivers' plane: the
/// seafloor's first solid node, vs there (truth and start), the harmonic means of vs over the
/// 250 m below it and the delay over 1 km, as `seafloor_speeds` works them out; then, over
/// every free node of the window (all y), the result's error ratio and cosine with the true
/// update by depth band, and the window's share of the band's added squared error
/// Σ (m - m_true)² - (m₀ - m_true)² (negative where the window gains). First the whole
/// section by depth band, for each iteration kept in `out/marmousi2/whole-it<n>`. Uses the
/// truth: a diagnosis, not a rule an inversion could follow.
#[test]
#[ignore = "a study: `just marmousi2-study whole_by_region` (a minute, after the whole inversion)"]
fn whole_by_region() {
    let base = PathBuf::from("examples/marmousi2");
    let config = |name: &str| Config::from_toml(&std::fs::read_to_string(base.join(name)).unwrap());
    let case = Case {
        start: config("whole.toml").unwrap(),
        truth: config("whole-true.toml").unwrap(),
        base,
        bands: vec![
            [500.0, 750.0],
            [750.0, 1000.0],
            [1000.0, 1500.0],
            [1500.0, 2500.0],
            [2500.0, 3125.0],
        ],
    };
    let (start, truth) = (case.model(&case.start), case.model(&case.truth));
    let result = |dir: &str| -> Option<Model> {
        let path = |p: &str| {
            PathBuf::from("out/marmousi2")
                .join(dir)
                .join(format!("{p}.npy"))
        };
        if !path("vp").is_file() {
            return None;
        }
        let mut m = start.clone();
        m.vp = vibe_rations::npy::read(&path("vp")).unwrap().1;
        m.vs = vibe_rations::npy::read(&path("vs")).unwrap().1;
        assert_eq!(m.vp.len(), start.vp.len(), "{dir} isn't whole.toml's grid");
        Some(m)
    };
    let mut last = None;
    for dir in ["whole-it2", "whole-it3", "whole-it4", "whole-it5", "whole"] {
        if let Some(m) = result(dir) {
            print_bands(dir, &case, &m, &start, &truth);
            last = Some(m);
        }
    }
    let Some(m) = last else {
        println!("no out/marmousi2/whole: run `invert examples/marmousi2/whole.toml` first");
        return;
    };

    let g = start.grid;
    let receiver = case.start.receivers.positions().unwrap()[0];
    let [_, j, _] = g.nearest_node(receiver).unwrap();
    let w = case.start.boundary.pml_width;
    let vs = |m: &Model, i: usize, k: usize| m.vs[g.interior_idx(i, j, k)] as f64;
    let seafloor = |i: usize| {
        (0..g.n[2])
            .find(|&k| truth.vs[g.interior_idx(i, j, k)] > 0.0)
            .unwrap()
    };
    let slowness = |m: &Model, i: usize| {
        let top = seafloor(i);
        (top..top + 10).map(|k| 1.0 / vs(m, i, k)).sum::<f64>() / 10.0
    };
    let x = |i: usize| i as f64 * g.h[0];
    let windows: Vec<[f64; 2]> = (0..16)
        .map(|k| [500.0 + 1000.0 * k as f64, 1500.0 + 1000.0 * k as f64])
        .collect();
    let free = case.banded(&start);
    let excess = |p: &[f32], p0: &[f32], t: &[f32], n: usize| {
        let (m, m0, t) = (p[n] as f64, p0[n] as f64, t[n] as f64);
        (m - t).powi(2) - (m0 - t).powi(2)
    };
    let band_excess: Vec<[f64; 2]> = free
        .iter()
        .map(|nodes| {
            let sum = |p: &[f32], p0: &[f32], t: &[f32]| {
                nodes.iter().map(|&n| excess(p, p0, t, n)).sum::<f64>()
            };
            [
                sum(&m.vp, &start.vp, &truth.vp),
                sum(&m.vs, &start.vs, &truth.vs),
            ]
        })
        .collect();
    println!(
        "by 1 km window along x: seafloor; vs on it, true / start; harmonic means over the 250 m \
         below, true / start; delay over 1 km, mean (largest); then by depth band (m) vp and vs: \
         error ratio, cosine, share of the band's added squared error"
    );
    for [lo, hi] in windows {
        let columns: Vec<usize> = (w..g.n[0] - w)
            .filter(|&i| x(i) >= lo && x(i) < hi)
            .collect();
        let mean = |f: &dyn Fn(usize) -> f64| {
            columns.iter().map(|&i| f(i)).sum::<f64>() / columns.len() as f64
        };
        let delay = |i: usize| 1000.0 * (slowness(&truth, i) - slowness(&start, i));
        let largest = columns
            .iter()
            .map(|&i| delay(i))
            .fold(f64::NEG_INFINITY, f64::max);
        let depth = mean(&|i| seafloor(i) as f64 * g.h[2]);
        println!(
            "{:5.1}–{:4.1} km: seafloor {depth:4.0} m; vs {:3.0} / {:3.0} m/s; 250 m {:3.0} / {:3.0} m/s; \
             delay {:.2} s ({largest:.2})",
            lo / 1000.0,
            hi / 1000.0,
            mean(&|i| vs(&truth, i, seafloor(i))),
            mean(&|i| vs(&start, i, seafloor(i))),
            1.0 / mean(&|i| slowness(&truth, i)),
            1.0 / mean(&|i| slowness(&start, i)),
            mean(&delay),
        );
        for ((band, nodes), [ep_total, es_total]) in
            case.bands.iter().zip(&free).zip(&band_excess).take(3)
        {
            let inside: Vec<usize> = nodes
                .iter()
                .copied()
                .filter(|&n| {
                    let xi = x(n % g.n[0]);
                    xi >= lo && xi < hi
                })
                .collect();
            let (ep, cp) = compare(&m.vp, &start.vp, &truth.vp, &inside);
            let (es, cs) = compare(&m.vs, &start.vs, &truth.vs, &inside);
            let share = |p: &[f32], p0: &[f32], t: &[f32], total: f64| {
                100.0 * inside.iter().map(|&n| excess(p, p0, t, n)).sum::<f64>() / total.abs()
            };
            println!(
                "    {:>9}  vp {ep:6.3} {cp:+6.2} {:+7.1}%   vs {es:6.3} {cs:+6.2} {:+7.1}%",
                format!("{:.0}–{:.0}", band[0], band[1]),
                share(&m.vp, &start.vp, &truth.vp, *ep_total),
                share(&m.vs, &start.vs, &truth.vs, *es_total),
            );
        }
    }
}

/// The slow waves' mute lines `slow_waves_mute` works out for each part, at once.
#[test]
#[ignore = "a study: `just marmousi2-study mute_lines` (seconds)"]
fn mute_lines() {
    for part in ["central", "west"] {
        let (velocity, mute) = slow_waves_mute(&Case::part(part, L2), 425.0);
        println!("{part}: {velocity:.0} m/s\n{mute}");
    }
}

/// Which side of its shot a receiver lies on: towards smaller x, or larger (the shot's own
/// column counts on both).
const SIDES: [&str; 2] = ["left", "right"];

/// How far along one side of a shot the slow-wave phase check lets the slow waves be kept.
#[derive(Clone, Copy, Debug)]
struct Limit {
    /// The offset (m) out to which the phase difference grows continuously from the shot and
    /// stays within ±π (interpolated to where it reaches π).
    offset: f64,
    /// The furthest offset on that side whose window holds anything.
    furthest: f64,
    /// Why it ends: `π` (reached), `jump` (more than π/2 between neighbours 25 m apart, which a
    /// phase drifting with offset never does), `start` (beyond π/2 at the shot already) or
    /// `end` (every offset whose window holds anything).
    why: &'static str,
}

impl std::fmt::Display for Limit {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.why == "end" {
            write!(f, "{:>5.0} all ", self.offset)
        } else {
            write!(f, "{:>5.0} {:<4}", self.offset, self.why)
        }
    }
}

/// `phase` wrapped into (-π, π].
fn wrap(phase: f64) -> f64 {
    let p = (phase + std::f64::consts::PI).rem_euclid(2.0 * std::f64::consts::PI)
        - std::f64::consts::PI;
    if p <= -std::f64::consts::PI {
        p + 2.0 * std::f64::consts::PI
    } else {
        p
    }
}

/// Unwraps `phases` (wrapped, one per receiver, nearest the shot first, at `offsets`) along
/// offset and finds where the unwrapped phase stops being continuous or reaches ±π: the check's
/// limit for one side of a shot, and the unwrapped phases up to it.
fn continuous_limit(offsets: &[f64], phases: &[f64]) -> (Limit, Vec<f64>) {
    let pi = std::f64::consts::PI;
    let furthest = *offsets.last().unwrap();
    let mut unwrapped = vec![phases[0]];
    if phases[0].abs() > 0.5 * pi {
        let limit = Limit {
            offset: offsets[0],
            furthest,
            why: "start",
        };
        return (limit, unwrapped);
    }
    for k in 1..offsets.len() {
        let step = wrap(phases[k] - phases[k - 1]);
        let last = unwrapped[k - 1];
        if step.abs() > 0.5 * pi {
            let limit = Limit {
                offset: offsets[k - 1],
                furthest,
                why: "jump",
            };
            return (limit, unwrapped);
        }
        let next = last + step;
        if next.abs() >= pi {
            let between = (pi - last.abs()) / (next.abs() - last.abs());
            let limit = Limit {
                offset: offsets[k - 1] + between * (offsets[k] - offsets[k - 1]),
                furthest,
                why: "π",
            };
            return (limit, unwrapped);
        }
        unwrapped.push(next);
    }
    let limit = Limit {
        offset: furthest,
        furthest,
        why: "end",
    };
    (limit, unwrapped)
}

/// One trace's spectrum at `f` Hz under the window `window` (one weight per sample):
/// Σ_m w_m d_m e^{-iω t_m}, t_m = m Δ.
fn spectrum(trace: &[f32], window: &[f64], f: f64, interval: f64) -> [f64; 2] {
    let omega = 2.0 * std::f64::consts::PI * f;
    let mut z = [0.0, 0.0];
    for (m, (&d, &w)) in trace.iter().zip(window).enumerate() {
        let (s, c) = (omega * m as f64 * interval).sin_cos();
        let v = w * d as f64;
        z[0] += v * c;
        z[1] -= v * s;
    }
    z
}

/// The slow-wave phase check of one part, from the observed traces and the start's
/// synthetics alone (no truth): in the window after the slow waves' line (the complement of
/// `mute`, the very samples the mute would remove), tapered by half a period of the source's
/// peak frequency before the record's end, each trace's spectrum at each of `frequencies`, and
/// the phase difference Δφ = arg(Σ_c w_c Z_syn,c Z̄_obs,c), w_c the components' balancing
/// weights from the observed traces: the cross-spectrum the balanced L2 misfit's alignment
/// term sums at that frequency, positive where the observed arrivals lag the synthetic ones.
/// Unwrapped along offset from the shot outwards on each side (`continuous_limit`), for both
/// components together and each alone. Returns the limits of both together, by frequency,
/// shot and side, and prints them with each component's and the phase profile every 250 m.
fn phase_limits(
    case: &Case,
    mute: &MuteConfig,
    frequencies: &[f64],
    observed: &[Observed],
) -> Vec<Vec<[Limit; 2]>> {
    let t = Instant::now();
    let f0 = 1.2;
    let mut limits = vec![Vec::new(); frequencies.len()];
    for (s, obs) in observed.iter().enumerate() {
        let shot = case.start.shot(s).unwrap();
        let mut sim = Simulation::new(&shot, &case.base).unwrap();
        sim.run(None, false).unwrap();
        let x = receiver_x(&sim);
        let x_source = shot.sources[0].position[0];
        let samples = sim.receivers.samples;
        let interval = sim.info.sample_interval;
        let duration = samples as f64 * interval;
        let end_taper = 0.5 / f0;
        let components = sim.receivers.components.clone();
        let balance = balancing_weights(&obs.traces, &components);
        let windows: Vec<Vec<f64>> = x
            .iter()
            .map(|&xr| {
                (0..samples)
                    .map(|m| {
                        let t = m as f64 * interval;
                        let end = ((duration - t) / end_taper).min(1.0);
                        let end = 0.5 * (1.0 - (std::f64::consts::PI * end).cos());
                        (1.0 - mute_weight(mute, t, xr - x_source)) * end
                    })
                    .collect()
            })
            .collect();
        // Past the offset where the line meets the record's end the window is empty: there
        // is nothing to judge, and nothing for the mute to remove.
        let filled: Vec<bool> = windows.iter().map(|w| w.iter().any(|&v| v > 0.0)).collect();
        for (fi, &f) in frequencies.iter().enumerate() {
            // Each receiver's cross-spectrum, both components together and each alone.
            let cross: Vec<Vec<[f64; 2]>> = windows
                .iter()
                .enumerate()
                .map(|(r, window)| {
                    let mut sum = [0.0, 0.0];
                    let mut each = Vec::new();
                    let range = r * samples..(r + 1) * samples;
                    let traces = sim.receivers.traces.iter().zip(&obs.traces);
                    for ((syn, obs), w) in traces.zip(&balance) {
                        let zs = spectrum(&syn[range.clone()], window, f, interval);
                        let zo = spectrum(&obs[range.clone()], window, f, interval);
                        let z = [zs[0] * zo[0] + zs[1] * zo[1], zs[1] * zo[0] - zs[0] * zo[1]];
                        sum[0] += w * z[0];
                        sum[1] += w * z[1];
                        each.push(z);
                    }
                    let mut all = vec![sum];
                    all.extend(each);
                    all
                })
                .collect();
            let mut shot_limits = [None, None];
            let mut lines = Vec::new();
            for (side, name) in SIDES.iter().enumerate() {
                // The side's receivers with something in their window, nearest the shot first;
                // not the shot's own, where vx vanishes.
                let mut on_side: Vec<usize> = (0..x.len())
                    .filter(|&r| {
                        filled[r]
                            && if side == 0 {
                                x[r] < x_source
                            } else {
                                x[r] > x_source
                            }
                    })
                    .collect();
                on_side
                    .sort_by(|&a, &b| (x[a] - x_source).abs().total_cmp(&(x[b] - x_source).abs()));
                let offsets: Vec<f64> = on_side.iter().map(|&r| (x[r] - x_source).abs()).collect();
                let mut cells = Vec::new();
                for (k, label) in ["both", "vx", "vz"].iter().enumerate() {
                    let phases: Vec<f64> = on_side
                        .iter()
                        .map(|&r| cross[r][k][1].atan2(cross[r][k][0]))
                        .collect();
                    let (limit, unwrapped) = continuous_limit(&offsets, &phases);
                    cells.push(format!("{label} {limit}"));
                    if k == 0 {
                        shot_limits[side] = Some(limit);
                        let profile: Vec<String> = offsets
                            .iter()
                            .zip(&unwrapped)
                            .filter(|(o, _)| (*o / 250.0 - (*o / 250.0).round()).abs() < 1e-6)
                            .map(|(o, p)| format!("{:.0}:{p:+.2}", o))
                            .collect();
                        lines.push(format!(
                            "    {name:>5} Δφ (rad) by offset (m): {}",
                            profile.join(" ")
                        ));
                    }
                }
                println!(
                    "  {f:.1} Hz, shot {s} at {x_source:.0} m, {name:>5} (window to {:.0} m): {}",
                    offsets.last().unwrap(),
                    cells.join(", ")
                );
            }
            for line in lines {
                println!("{line}");
            }
            limits[fi].push(shot_limits.map(Option::unwrap));
        }
    }
    println!("  the check: {:.0} s", t.elapsed().as_secs_f64());
    limits
}

/// What the truth says the check should find (a check of the check, never part of the rule):
/// the phase 2πf τ a wave along the seafloor gathers between a shot and each receiver, τ the
/// sum over the columns between of the difference of the harmonic mean slownesses over the
/// 250 m below the seafloor, true less the start's (as `seafloor_speeds` works them out), and
/// the offset on each side where it reaches π. By frequency, shot and side.
fn predicted_limits(case: &Case, frequencies: &[f64]) -> Vec<Vec<[f64; 2]>> {
    let (start, truth) = (case.model(&case.start), case.model(&case.truth));
    let g = start.grid;
    let receivers = case.start.receivers.positions().unwrap();
    let [_, j, top] = g.nearest_node(receivers[0]).unwrap();
    let slowness = |m: &Model, i: usize| {
        (top..top + 10)
            .map(|k| 1.0 / m.vs[g.interior_idx(i, j, k)] as f64)
            .sum::<f64>()
            / 10.0
    };
    let excess: Vec<f64> = (0..g.n[0])
        .map(|i| g.h[0] * (slowness(&truth, i) - slowness(&start, i)))
        .collect();
    let furthest = |x_source: f64, side: usize| {
        receivers
            .iter()
            .map(|p| p[0] - x_source)
            .filter(|d| if side == 0 { *d <= 0.0 } else { *d >= 0.0 })
            .fold(0.0f64, |m, d| m.max(d.abs()))
    };
    frequencies
        .iter()
        .map(|&f| {
            case.start
                .sources
                .iter()
                .map(|source| {
                    let i0 = (source.position[0] / g.h[0]).round() as usize;
                    [0, 1].map(|side| {
                        let reach = furthest(source.position[0], side);
                        let mut delay = 0.0;
                        let mut i = i0;
                        loop {
                            let offset = (i as f64 - i0 as f64).abs() * g.h[0];
                            if offset >= reach {
                                return reach;
                            }
                            i = if side == 0 { i - 1 } else { i + 1 };
                            let next = delay + excess[i];
                            if (2.0 * std::f64::consts::PI * f * next).abs() >= std::f64::consts::PI
                            {
                                let crossing = (0.5 / f - delay.abs()) / (next.abs() - delay.abs());
                                return offset + crossing * g.h[0];
                            }
                            delay = next;
                        }
                    })
                })
                .collect()
        })
        .collect()
}

/// The frequencies the check measures at: the stages' corners (`FREQUENCY_STAGES`) and the
/// source's peak frequency, which stands for the unfiltered traces.
const CHECK_FREQUENCIES: [f64; 3] = [0.6, 1.0, 1.2];

/// The slow waves' mutes one frequency's `limits` (by shot and side) call for, `mute` the
/// line, as `[[table]]`s: for each shot with a side cut short, a mute for that shot alone,
/// sparing the offsets out to each side's limit less `spare_taper`, so that what it keeps
/// whole lies within ±π. The offset taper is the distance the line moves over its time taper,
/// velocity × taper (259 m on the central part): as gradual along offset as in time. A side
/// the check passes to its window's end is spared entirely.
fn phase_mutes(table: &str, mute: &MuteConfig, limits: &[[Limit; 2]]) -> String {
    let spare_taper = mute.velocity.unwrap() * mute.taper;
    let end = |l: &Limit| {
        if l.why == "end" {
            1e6
        } else {
            (l.offset - spare_taper).max(0.0)
        }
    };
    limits
        .iter()
        .enumerate()
        .filter(|(_, l)| l.iter().any(|l| l.why != "end"))
        .map(|(s, [left, right])| {
            format!(
                "[[{table}]]\nkeep = \"before\"\nvelocity = {:.1}\ndelay = {:.2}\ntaper = {:.2}\n\
                 shot = {s}\nspare = [{:.1}, {:.1}]\nspare_taper = {spare_taper:.1}",
                mute.velocity.unwrap(),
                mute.delay,
                mute.taper,
                -end(left),
                end(right),
            )
        })
        .collect::<Vec<_>>()
        .join("\n")
}

/// The check on `case` with its own mute line (`slow_waves_mute`): the line's mute, and the
/// limits by frequency (`CHECK_FREQUENCIES`), shot and side. Forward runs of the truth (the
/// observed traces) and of the start only.
fn checked(case: &Case) -> (MuteConfig, Vec<Vec<[Limit; 2]>>) {
    let (velocity, mute) = slow_waves_mute(case, 425.0);
    let mute = mutes_of(muted_with("iterations = 0", false, &mute))[0].clone();
    println!("the window after the line at {velocity:.0} m/s");
    let limits = phase_limits(case, &mute, &CHECK_FREQUENCIES, &case.record());
    (mute, limits)
}

/// L2 for five iterations with the slow waves muted only where the check finds them beyond
/// half a cycle at the source's peak frequency, shot by shot and side by side.
fn phase_checked_l2(part: &str) -> Case {
    let (mute, limits) = checked(&Case::part(part, L2));
    let mutes = phase_mutes("inversion.mute", &mute, &limits[2]);
    println!("the mutes:\n{mutes}");
    Case::part(part, ("iterations = 5", leak(mutes)))
}

/// The frequency stages (0.6 Hz, 1.0 Hz, all) with each stage's mutes from the check at its
/// own frequency.
fn phase_checked_stages(part: &str) -> Case {
    let (mute, limits) = checked(&Case::part(part, L2));
    let stage = |f: usize| phase_mutes("inversion.stage.mute", &mute, &limits[f]);
    let tables = format!(
        "[[inversion.stage]]\nmax_frequency = 0.6\niterations = 2\n{}\n\
         [[inversion.stage]]\nmax_frequency = 1.0\niterations = 2\n{}\n\
         [[inversion.stage]]\niterations = 1\n{}",
        stage(0),
        stage(1),
        stage(2)
    );
    println!("the stages:\n{tables}");
    Case::part(part, ("", leak(tables)))
}

/// `phase_mutes` writes mutes that parse, one per shot cut short, for that shot alone and
/// sparing what the check passed.
#[test]
fn phase_mutes_parse() {
    let mute = mutes_of(muted(false, SLOW_WAVES_MUTED))[0].clone();
    let limit = |offset: f64, why: &'static str| Limit {
        offset,
        furthest: 1500.0,
        why,
    };
    let limits = [
        [limit(1500.0, "end"), limit(1500.0, "end")],
        [limit(800.0, "π"), limit(1500.0, "end")],
        [limit(25.0, "start"), limit(300.0, "jump")],
    ];
    let mutes = phase_mutes("inversion.mute", &mute, &limits);
    let cfg = Config::from_toml(&central_text(("iterations = 5", leak(mutes)))).unwrap();
    let mutes = cfg.inversion.unwrap().mutes;
    assert_eq!(mutes.len(), 2);
    assert_eq!((mutes[0].shot, mutes[1].shot), (Some(1), Some(2)));
    assert_eq!(mutes[0].spare, Some([-548.0, 1e6]));
    assert_eq!(mutes[1].spare, Some([0.0, 48.0]));
    assert_eq!(mutes[0].spare_taper, 252.0);
}

/// The slow-wave phase check on both parts, at 0.6, 1.0 and 1.2 Hz, with the mute line
/// `slow_waves_mute` works out from each part's start, next to what the truth predicts.
#[test]
#[ignore = "a study: `just marmousi2-study phase_check` (about 3 minutes)"]
fn phase_check() {
    for part in ["central", "west"] {
        let case = Case::part(part, L2);
        println!("{part}:");
        let (mute, limits) = checked(&case);
        let predicted = predicted_limits(&case, &CHECK_FREQUENCIES);
        println!(
            "{part}: the limits (m) by shot, left / right, from the data; then from the truth's \
             delays"
        );
        for (fi, f) in CHECK_FREQUENCIES.iter().enumerate() {
            for (s, (l, p)) in limits[fi].iter().zip(&predicted[fi]).enumerate() {
                println!(
                    "  {f:.1} Hz, shot {s}: {} / {} (windows to {:.0} / {:.0}); truth {:.0} / \
                     {:.0}",
                    l[0],
                    l[1],
                    l[0].furthest,
                    l[1].furthest,
                    p[0].min(l[0].furthest),
                    p[1].min(l[1].furthest)
                );
            }
            println!("{}", phase_mutes(&format!("{f} Hz"), &mute, &limits[fi]));
        }
    }
}

#[test]
#[ignore = "a study: `just marmousi2-study central_phase_checked_l2` (about 45 minutes)"]
fn central_phase_checked_l2() {
    invert("phase-checked-l2", &phase_checked_l2("central"));
}

#[test]
#[ignore = "a study: `just marmousi2-study west_phase_checked_l2` (about 40 minutes)"]
fn west_phase_checked_l2() {
    invert("west-phase-checked-l2", &phase_checked_l2("west"));
}

#[test]
#[ignore = "a study: `just marmousi2-study central_phase_checked_stages` (about an hour)"]
fn central_phase_checked_stages() {
    invert("phase-checked-stages", &phase_checked_stages("central"));
}

#[test]
#[ignore = "a study: `just marmousi2-study west_phase_checked_stages` (about an hour)"]
fn west_phase_checked_stages() {
    invert("west-phase-checked-stages", &phase_checked_stages("west"));
}

#[test]
#[ignore = "a study: `just marmousi2-study west_l2` (about 45 minutes)"]
fn west_l2() {
    invert("west-l2", &Case::part("west", L2));
}

#[test]
#[ignore = "a study: `just marmousi2-study west_muted_l2` (about 45 minutes)"]
fn west_muted_l2() {
    let (velocity, mute) = slow_waves_mute(&Case::part("west", L2), 425.0);
    println!("the mute line: {velocity:.0} m/s");
    invert(
        "west-muted-l2",
        &Case::part("west", muted_with("iterations = 5", false, &mute)),
    );
}

#[test]
#[ignore = "a study: `just marmousi2-study l2` (about 40 minutes)"]
fn l2() {
    invert("l2", &Case::central(L2));
}

#[test]
#[ignore = "a study: `just marmousi2-study l2_checkpointed` (about 35 minutes, after `l2`)"]
fn l2_checkpointed() {
    let case = Case::central(L2_CHECKPOINTED);
    invert("l2-checkpointed", &case);
    compare_variants(&case, "l2", "l2-checkpointed", 0);
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
#[ignore = "a study: `just marmousi2-study muted_l2` (about 45 minutes)"]
fn muted_l2() {
    invert("muted-l2", &Case::central(muted(false, SLOW_WAVES_MUTED)));
}

#[test]
#[ignore = "a study: `just marmousi2-study muted_stages` (about an hour)"]
fn muted_stages() {
    invert(
        "muted-stages",
        &Case::central(muted(true, SLOW_WAVES_MUTED)),
    );
}

#[test]
#[ignore = "a study: `just marmousi2-study down_weighted_l2` (about 45 minutes)"]
fn down_weighted_l2() {
    invert(
        "down-weighted-l2",
        &Case::central(muted(false, SLOW_WAVES_DOWN_WEIGHTED)),
    );
}

#[test]
#[ignore = "a study: `just marmousi2-study down_weighted_stages` (about an hour)"]
fn down_weighted_stages() {
    invert(
        "down-weighted-stages",
        &Case::central(muted(true, SLOW_WAVES_DOWN_WEIGHTED)),
    );
}

#[test]
#[ignore = "a study: `just marmousi2-study muted_l2_with_receivers_side` (about 50 minutes)"]
fn muted_l2_with_receivers_side() {
    invert(
        "muted-l2-with-receivers-side",
        &Case::central(muted_l2_with_receivers_side_variant()),
    );
}

#[test]
#[ignore = "a study: `just marmousi2-study muted_l2_longer` (about 2.5 hours)"]
fn muted_l2_longer() {
    invert("muted-l2-longer", &Case::central(muted_l2_longer_variant()));
}

#[test]
#[ignore = "a study: `just marmousi2-study stages_muted_after_the_first` (about an hour)"]
fn stages_muted_after_the_first() {
    invert(
        "stages-muted-after-the-first",
        &Case::central(stages_muted_after_the_first_variant()),
    );
}

#[test]
#[ignore = "a study: `just marmousi2-study windowed_time_shifts_then_l2` (about 55 minutes)"]
fn windowed_time_shifts_then_l2() {
    invert(
        "windowed-time-shifts-then-l2",
        &Case::central(windowed_time_shifts_then_l2_variant()),
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
