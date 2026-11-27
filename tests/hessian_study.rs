//! Hessian probes on the sphere survey, with no change to the solver:
//! `cargo test --release --test hessian_study -- --ignored --nocapture --test-threads 1`.
//!
//! The start is the observed data's own model, so the residuals vanish and the Hessian is
//! the Gauss–Newton one, JᵀWJ (Fichtner & Trampert 2011: with no residual the full and
//! approximate Hessians coincide). The unknowns are the inversion's, relative changes
//! m/m₀ - 1 of vp, vs and ρ at every free node.
//!
//! - `point_spread_functions`: H applied to a Gaussian bump of vp, of vs and of ρ at one
//!   point, from central differences of the exact gradient extrapolated in the step
//!   (Richardson). These are rows of the Hessian, Fichtner & Trampert's Hessian kernels
//!   (Fichtner 2011's resolution and trade-off kernels): how far each spreads, how much of it
//!   lands on the other parameters, and the 3 × 3 block the bumps span, whose symmetry checks
//!   the differences.
//! - `gauss_newton_diagonal`: the diagonal of JᵀWJ for each parameter from adjoint runs with
//!   random residuals: with independent signs z, E[(Jᵀ√W z)²] = diag(JᵀWJ). Checked against the
//!   exact diagonal at the bump's centre, and set beside the source-side illumination the
//!   inversion divides by, depth by depth.

use std::path::Path;
use std::time::Instant;

use vibe_rations::config::{Config, Parameter};
use vibe_rations::gradient::gradient_with;
use vibe_rations::grid::Grid;
use vibe_rations::inversion::{Inversion, Observed};
use vibe_rations::model::Model;
use vibe_rations::simulation::Simulation;

/// Where the bumps sit: the centre of `inversion_study`'s sphere, 85 m below the sources.
const CENTRE: [f64; 3] = [175.0, 155.0, 115.0];

/// `inversion_study`'s sphere survey with no sphere, inverting all three parameters.
fn survey() -> String {
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
    [[source]]
    position = [110.0, 155.0, 30.0]
    force = [0.0, 0.0, 1e9]
    wavelet = { kind = "ricker", frequency = 12.0 }
    [[source]]
    position = [240.0, 150.0, 30.0]
    force = [1e9, 0.0, 0.0]
    wavelet = { kind = "ricker", frequency = 12.0 }
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
    parameters = ["vp", "vs", "rho"]
    iterations = 1
    "#
    .to_string()
}

/// The inversion of the survey's own data from its own model: its shots' configurations
/// (with Δt, the layers' damping and the coefficients' band pinned) record the observed
/// traces, so the residuals are zero to the last bit.
fn setup() -> Inversion {
    let cfg = Config::from_toml(&survey()).unwrap();
    let base = Path::new(".");
    let start = Model::from_config(Simulation::grid(&cfg).unwrap(), &cfg.model, base).unwrap();
    let settings = cfg.inversion.clone().unwrap();
    let record = |shot: &Config| -> Observed {
        let mut sim = Simulation::with_model(shot, start.clone()).unwrap();
        sim.run(None, false).unwrap();
        sim.receivers.traces.into()
    };
    let first: Vec<Observed> = (0..cfg.sources.len())
        .map(|s| record(&cfg.shot(s).unwrap()))
        .collect();
    let inv = Inversion::new(&cfg, start.clone(), first, &settings).unwrap();
    let observed = inv.shots.iter().map(|s| record(&s.config)).collect();
    Inversion::new(&cfg, start, observed, &settings).unwrap()
}

/// Interior node `n`'s position.
fn position(g: &Grid, n: usize) -> [f64; 3] {
    let [nx, ny, _] = g.n;
    g.position([n % nx, (n / nx) % ny, n / (nx * ny)])
}

fn distance2(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| (a[i] - b[i]).powi(2)).sum()
}

/// Where `p`'s unknowns start among all of them.
fn offset(inv: &Inversion, p: Parameter) -> usize {
    let mut o = 0;
    for (q, nodes) in inv.free() {
        if *q == p {
            return o;
        }
        o += nodes.len();
    }
    panic!("{} is not inverted", p.name())
}

/// A relative change of `p` shaped as a Gaussian of standard deviation `sigma` metres about
/// `CENTRE` (cut at 3σ), or, with `sigma` 0, a unit change at the node nearest to it.
fn bump(inv: &Inversion, p: Parameter, sigma: f64) -> Vec<f64> {
    let g = inv.model.grid;
    let mut x = vec![0.0; inv.unknowns()];
    let o = offset(inv, p);
    let nodes = &inv.free().iter().find(|(q, _)| *q == p).unwrap().1;
    if sigma == 0.0 {
        let nearest = (0..nodes.len())
            .min_by(|&a, &b| {
                let d = |i: usize| distance2(position(&g, nodes[i]), CENTRE);
                d(a).total_cmp(&d(b))
            })
            .unwrap();
        x[o + nearest] = 1.0;
        return x;
    }
    for (i, &n) in nodes.iter().enumerate() {
        let r2 = distance2(position(&g, n), CENTRE);
        if r2 < 9.0 * sigma * sigma {
            x[o + i] = (-r2 / (2.0 * sigma * sigma)).exp();
        }
    }
    x
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

/// H v from central differences of the gradient at steps h and h/2, extrapolated to h → 0,
/// and the relative difference between the two before extrapolating (the h² error).
fn hessian_times(inv: &Inversion, v: &[f64], h: f64) -> (Vec<f64>, f64) {
    let gradient = |s: f64| {
        let x: Vec<f64> = v.iter().map(|a| s * a).collect();
        inv.evaluate(&inv.model_at(&x)).unwrap().gradient
    };
    let difference = |h: f64| -> Vec<f64> {
        let (plus, minus) = (gradient(h), gradient(-h));
        plus.iter()
            .zip(&minus)
            .map(|(a, b)| (a - b) / (2.0 * h))
            .collect()
    };
    let (coarse, fine) = (difference(h), difference(h / 2.0));
    let gap: Vec<f64> = coarse.iter().zip(&fine).map(|(a, b)| a - b).collect();
    let hv = coarse
        .iter()
        .zip(&fine)
        .map(|(a, b)| (4.0 * b - a) / 3.0)
        .collect();
    (hv, norm(&gap) / norm(&fine))
}

/// The part of `x` that belongs to `p`.
fn part<'a>(inv: &Inversion, x: &'a [f64], p: Parameter) -> &'a [f64] {
    let o = offset(inv, p);
    let len = inv.free().iter().find(|(q, _)| *q == p).unwrap().1.len();
    &x[o..o + len]
}

/// `x`'s part for `p` as a volume (nz, ny, nx), zero at the fixed nodes.
fn volume(inv: &Inversion, x: &[f64], p: Parameter) -> Vec<f32> {
    let g = inv.model.grid;
    let mut out = vec![0.0f32; g.interior_len()];
    let nodes = &inv.free().iter().find(|(q, _)| *q == p).unwrap().1;
    for (&n, &value) in nodes.iter().zip(part(inv, x, p)) {
        out[n] = value as f32;
    }
    out
}

/// The centroid and rms half-widths along x, y and z of |part|², in metres.
fn spread(inv: &Inversion, x: &[f64], p: Parameter) -> ([f64; 3], [f64; 3]) {
    let g = inv.model.grid;
    let nodes = &inv.free().iter().find(|(q, _)| *q == p).unwrap().1;
    let values = part(inv, x, p);
    let total: f64 = values.iter().map(|v| v * v).sum();
    let mut centre = [0.0; 3];
    for (&n, v) in nodes.iter().zip(values) {
        let r = position(&g, n);
        (0..3).for_each(|a| centre[a] += v * v * r[a] / total);
    }
    let mut width = [0.0; 3];
    for (&n, v) in nodes.iter().zip(values) {
        let r = position(&g, n);
        (0..3).for_each(|a| width[a] += v * v * (r[a] - centre[a]).powi(2) / total);
    }
    (centre, width.map(f64::sqrt))
}

#[test]
#[ignore]
fn point_spread_functions() {
    let t0 = Instant::now();
    let inv = setup();
    println!(
        "misfit at the start {:.3e}, {} unknowns ({:.0} s)",
        inv.misfit(),
        inv.unknowns(),
        t0.elapsed().as_secs_f64()
    );
    let out = Path::new("out/hessian");
    std::fs::create_dir_all(out).unwrap();
    let g = inv.model.grid;
    let shape = g.shape();
    let params = [Parameter::Vp, Parameter::Vs, Parameter::Rho];
    let sigma = 10.0;
    let bumps: Vec<Vec<f64>> = params.iter().map(|&p| bump(&inv, p, sigma)).collect();
    let mut rows = Vec::new();
    for (p, b) in params.iter().zip(&bumps) {
        let t = Instant::now();
        let (hv, gap) = hessian_times(&inv, b, 0.02);
        println!(
            "H on a {} bump (σ {sigma} m at {CENTRE:?}): h vs h/2 differ by {gap:.1e} ({:.0} s)",
            p.name(),
            t.elapsed().as_secs_f64()
        );
        for q in params {
            let (centre, width) = spread(&inv, &hv, q);
            let share = norm(part(&inv, &hv, q)) / norm(part(&inv, &hv, *p));
            println!(
                "    on {:>3}: norm {:.3} of its own; centroid ({:.0}, {:.0}, {:.0}) m; rms \
                 half-widths ({:.0}, {:.0}, {:.0}) m",
                q.name(),
                share,
                centre[0],
                centre[1],
                centre[2],
                width[0],
                width[1],
                width[2]
            );
            let name = format!("{}-on-{}.npy", p.name(), q.name());
            vibe_rations::npy::write(&out.join(name), &shape, &volume(&inv, &hv, q)).unwrap();
        }
        rows.push(hv);
    }
    // The 3 × 3 block bᵢᵀ H bⱼ, its symmetry and its correlations.
    let block: Vec<Vec<f64>> = bumps
        .iter()
        .map(|bi| rows.iter().map(|hj| dot(bi, hj)).collect())
        .collect();
    println!("bᵢᵀ H bⱼ (rows and columns vp, vs, ρ):");
    for (i, row) in block.iter().enumerate() {
        println!(
            "    {:>3}: {}",
            params[i].name(),
            row.iter()
                .map(|v| format!("{v:>11.4e}"))
                .collect::<Vec<_>>()
                .join(" ")
        );
    }
    for i in 0..3 {
        for j in i + 1..3 {
            let scale = (block[i][i] * block[j][j]).sqrt();
            println!(
                "    {}–{}: correlation {:+.3}, asymmetry {:.1e}",
                params[i].name(),
                params[j].name(),
                0.5 * (block[i][j] + block[j][i]) / scale,
                (block[i][j] - block[j][i]).abs() / scale
            );
        }
    }
    println!("{:.0} s", t0.elapsed().as_secs_f64());
}

/// SplitMix64 (Steele, Lea & Flood 2014), for the probes' signs.
fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

/// For each probe, the gradient of every shot with residuals √w ± 1 at every sample (w each
/// trace's weight in the misfit), summed over the shots and taken to the unknowns' relative
/// units: the squares of these average to diag(JᵀWJ) for vp, vs and ρ at every node.
fn probes(inv: &Inversion, count: usize, seed: u64) -> Vec<[Vec<f64>; 3]> {
    let model = &inv.model;
    let n = model.grid.interior_len();
    let mut rng = seed;
    (0..count)
        .map(|_| {
            let mut sum = [vec![0.0; n], vec![0.0; n], vec![0.0; n]];
            for shot in &inv.shots {
                let mut sim = Simulation::with_model(&shot.config, model.clone()).unwrap();
                let weights = shot.misfit.receiver_weights();
                let rng = &mut rng;
                let run = gradient_with(
                    &mut sim,
                    |r| {
                        let residuals = r
                            .traces
                            .iter()
                            .zip(&weights)
                            .map(|(trace, w)| {
                                trace
                                    .chunks(r.samples)
                                    .zip(w)
                                    .flat_map(|(t, &w)| {
                                        t.iter()
                                            .map(|_| {
                                                let sign = if splitmix64(rng) >> 63 == 0 {
                                                    1.0
                                                } else {
                                                    -1.0
                                                };
                                                (sign * w.sqrt()) as f32
                                            })
                                            .collect::<Vec<_>>()
                                    })
                                    .collect()
                            })
                            .collect();
                        Ok((0.0, residuals))
                    },
                    inv.strategy(),
                )
                .unwrap();
                for (s, g) in sum.iter_mut().zip(run.gradient.velocities(model)) {
                    s.iter_mut().zip(g).for_each(|(s, g)| *s += g);
                }
            }
            for (p, s) in Parameter::ALL.iter().zip(sum.iter_mut()) {
                let m0 = model.values(*p);
                s.iter_mut().zip(m0).for_each(|(s, &m)| *s *= m as f64);
            }
            sum
        })
        .collect()
}

#[test]
#[ignore]
fn gauss_newton_diagonal() {
    let t0 = Instant::now();
    let inv = setup();
    let g = inv.model.grid;
    let count = 16;
    let samples = probes(&inv, count, 20_261_127);
    println!(
        "{count} probes of {} shots ({:.0} s)",
        inv.shots.len(),
        t0.elapsed().as_secs_f64()
    );
    let n = g.interior_len();
    // The estimate and its standard error at every node, for each parameter.
    let mut mean = [vec![0.0; n], vec![0.0; n], vec![0.0; n]];
    let mut sq = [vec![0.0; n], vec![0.0; n], vec![0.0; n]];
    for s in &samples {
        for c in 0..3 {
            for i in 0..n {
                let v = s[c][i] * s[c][i];
                mean[c][i] += v / count as f64;
                sq[c][i] += v * v / count as f64;
            }
        }
    }
    // The exact diagonal at the node nearest the centre, from H on a unit change there.
    let mut exact = Vec::new();
    for p in [Parameter::Vp, Parameter::Vs, Parameter::Rho] {
        let e = bump(&inv, p, 0.0);
        let (he, _) = hessian_times(&inv, &e, 0.02);
        let k = e.iter().position(|&v| v == 1.0).unwrap();
        let o = offset(&inv, p);
        let node = inv.free().iter().find(|(q, _)| *q == p).unwrap().1[k - o];
        let c = p as usize;
        let se = ((sq[c][node] - mean[c][node].powi(2)) / count as f64).sqrt();
        println!(
            "{:>3} at {:?}: exact {:.4e}, probes {:.4e} ± {:.1e} ({:+.0}%)",
            p.name(),
            position(&g, node),
            he[k],
            mean[c][node],
            se,
            100.0 * (mean[c][node] / he[k] - 1.0)
        );
        exact.push((p, node, he[k]));
    }
    // Depth by depth over the free nodes: each parameter's mean diagonal, and the source-side
    // illumination the inversion divides by, each relative to its value at the centre's depth.
    let illumination = inv.evaluate(&inv.model).unwrap().illumination;
    let free = &inv.free()[0].1;
    let [nx, ny, nz] = g.n;
    let depth_of = |n: usize| n / (nx * ny);
    let k0 = (CENTRE[2] / g.h[2]).round() as usize;
    let layer_mean = |values: &[f64], k: usize| {
        let at: Vec<f64> = free
            .iter()
            .filter(|&&n| depth_of(n) == k)
            .filter(|&&n| {
                let r = position(&g, n);
                (r[0] - CENTRE[0]).abs() <= 50.0 && (r[1] - CENTRE[1]).abs() <= 50.0
            })
            .map(|&n| values[n])
            .collect();
        at.iter().sum::<f64>() / at.len().max(1) as f64
    };
    println!(
        "depth   vp      vs      rho     illumination   (means within 50 m of the centre's column, / their value at {} m)",
        k0 as f64 * g.h[2]
    );
    let reference: Vec<f64> = (0..3)
        .map(|c| layer_mean(&mean[c], k0))
        .chain([layer_mean(&illumination, k0)])
        .collect();
    for k in 0..nz {
        if !free.iter().any(|&n| depth_of(n) == k) {
            continue;
        }
        let values: Vec<f64> = (0..3)
            .map(|c| layer_mean(&mean[c], k))
            .chain([layer_mean(&illumination, k)])
            .collect();
        println!(
            "{:>5.0}  {}",
            k as f64 * g.h[2],
            values
                .iter()
                .zip(&reference)
                .map(|(v, r)| format!("{:>7.3}", v / r))
                .collect::<Vec<_>>()
                .join(" ")
        );
    }
    let out = Path::new("out/hessian");
    std::fs::create_dir_all(out).unwrap();
    let shape = g.shape();
    for (p, m) in Parameter::ALL.iter().zip(&mean) {
        let v: Vec<f32> = m.iter().map(|&x| x as f32).collect();
        vibe_rations::npy::write(&out.join(format!("diag-{}.npy", p.name())), &shape, &v).unwrap();
    }
    let v: Vec<f32> = illumination.iter().map(|&x| x as f32).collect();
    vibe_rations::npy::write(&out.join("illumination.npy"), &shape, &v).unwrap();
    println!("{:.0} s", t0.elapsed().as_secs_f64());
}
