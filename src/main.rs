use std::path::{Path, PathBuf};
use std::process::ExitCode;

use anyhow::{Context, Result, ensure};
use clap::{Parser, Subcommand};

use vibe_rations::config::{Config, shot_dir};
use vibe_rations::gradient::{Strategy, balancing_weights, gradient_with, l2_misfit};
use vibe_rations::inversion::Inversion;
use vibe_rations::npy;
use vibe_rations::simulation::{self, Simulation};
use vibe_rations::{memory, remote, section, server};

/// 3D elastic seismic wave simulator.
#[derive(Parser)]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Command,
    /// Worker threads; by default one per core.
    #[arg(long, global = true)]
    threads: Option<usize>,
}

#[derive(Subcommand)]
enum Command {
    /// Run a simulation, here or on a machine running `vibe-rations serve`.
    Run {
        /// The TOML configuration.
        config: PathBuf,
        /// Where to write results, instead of the configuration's [output] dir.
        #[arg(long)]
        out: Option<PathBuf>,
        /// Don't report progress.
        #[arg(long, short)]
        quiet: bool,
        /// Run it on another machine: an SSH destination (user@host, or a Host from
        /// ~/.ssh/config) whose `vibe-rations serve` this reaches through an SSH tunnel, or the
        /// server's URL. The run carries on there if this stops; `remote attach` picks it up.
        #[arg(long)]
        remote: Option<String>,
        /// The port `vibe-rations serve` listens on, on the remote machine.
        #[arg(long, default_value_t = 8088)]
        port: u16,
        /// With --remote: start the run and return at once, without following it. `remote
        /// list` shows its progress and `remote fetch` its results.
        #[arg(long, requires = "remote")]
        detach: bool,
        /// Fire the sources one at a time, each shot's results in its own shot-NNN
        /// directory, as an inversion takes them.
        #[arg(long, conflicts_with = "remote")]
        shots: bool,
    },
    /// Keep simulations for clients that come and go: the command line's `--remote` and the
    /// web page. Listens on localhost only by default; reach it through an SSH tunnel.
    Serve {
        /// Address to listen on.
        #[arg(long, default_value = "127.0.0.1:8088")]
        listen: String,
        /// Where runs keep their configurations, model files and results.
        #[arg(long, default_value = "runs")]
        data: PathBuf,
        /// The web page's files (`web/dist`), to serve at /.
        #[arg(long)]
        page: Option<PathBuf>,
        /// Print every request and its answer.
        #[arg(long)]
        log: bool,
    },
    /// Look after the runs on a machine running `vibe-rations serve`.
    Remote {
        /// An SSH destination (user@host, or a Host from ~/.ssh/config), or the server's URL.
        target: String,
        /// The port `vibe-rations serve` listens on, on the remote machine.
        #[arg(long, default_value_t = 8088)]
        port: u16,
        #[command(subcommand)]
        action: RemoteAction,
    },
    /// Set up a simulation without running it: grid, time step, accuracy and memory.
    Check {
        /// The TOML configuration.
        config: PathBuf,
    },
    /// Invert observed traces for the model with L-BFGS, as the configuration's [inversion]
    /// section says: each source is a shot.
    ///
    /// Writes the model as it stands after each iteration (vp.npy, vs.npy, rho.npy, shape nz,
    /// ny, nx, ready for [model.files]) and inversion.json, the settings and every
    /// iteration's misfit.
    Invert {
        /// The TOML configuration, with the starting model and an [inversion] section.
        config: PathBuf,
        /// Where to write results, instead of the configuration's [output] dir.
        #[arg(long)]
        out: Option<PathBuf>,
    },
    /// The gradient of the L2 misfit to observed traces with respect to vp, vs and ρ.
    ///
    /// Runs the configured model, then the adjoint of its receivers, and writes
    /// grad_vp.npy, grad_vs.npy and grad_rho.npy (shape nz, ny, nx, like model files), the
    /// synthetic traces and gradient.json. Exact for the discrete run outside the absorbing
    /// layers, with the time step and layers held fixed.
    Gradient {
        /// The TOML configuration of the model to differentiate around.
        config: PathBuf,
        /// A directory with the observed traces, `<component>.npy` for every recorded
        /// component, shaped as `run` writes them.
        #[arg(long)]
        observed: PathBuf,
        /// Where to write results, instead of the configuration's [output] dir.
        #[arg(long)]
        out: Option<PathBuf>,
        /// Weigh each component by 1/(its observed energy), so pressure and velocity count
        /// alike.
        #[arg(long)]
        balance: bool,
        #[command(flatten)]
        strategy: StrategyArgs,
    },
    /// Make model files from a 2D section: SEG-Y files of vp, vs and ρ, a trace per x position
    /// and a sample per depth, the first at the surface.
    ///
    /// Averages them over cells `spacing` wide around each node (ρ arithmetically, ρvp² and
    /// ρvs² harmonically; a node is fluid if most of its cell is). Writes vp.npy, vs.npy and
    /// rho.npy for [model.files]: sections, shape (nz, nx), which a run repeats along y to
    /// whatever the grid has, or with --ny the whole grid, (nz, ny, nx). Prints the [grid]
    /// they fit.
    Extrude {
        /// P-wave speed, SEG-Y.
        vp: PathBuf,
        /// S-wave speed, SEG-Y; 0 in fluids.
        vs: PathBuf,
        /// Density, SEG-Y.
        rho: PathBuf,
        /// The section's sample spacing in metres, along x and z alike.
        #[arg(long)]
        input_spacing: f64,
        /// The nodes' spacing in metres.
        #[arg(long)]
        spacing: f64,
        /// Repeat the section this many times along y, writing the whole grid; by default
        /// the files are sections, a twentieth or less of the size.
        #[arg(long)]
        ny: Option<usize>,
        /// Keep only this part along x, in metres from the first trace: FROM,TO.
        #[arg(long, value_delimiter = ',', value_name = "FROM,TO")]
        x: Option<Vec<f64>>,
        /// Keep only this part along z, in metres from the first sample, which becomes the
        /// surface: FROM,TO.
        #[arg(long, value_delimiter = ',', value_name = "FROM,TO")]
        z: Option<Vec<f64>>,
        /// Multiply the speeds by this to make them m/s: 1000 for km/s.
        #[arg(long, default_value_t = 1.0)]
        speed_scale: f64,
        /// Multiply the densities by this to make them kg/m³: 1000 for g/cm³.
        #[arg(long, default_value_t = 1.0)]
        rho_scale: f64,
        /// Smooth the solid part with a Gaussian of this standard deviation in metres, leaving
        /// the fluid as it is: a starting model for an inversion.
        #[arg(long, value_name = "METRES")]
        smooth: Option<f64>,
        /// With --smooth, smooth the slownesses 1/vp and 1/vs instead of the speeds (harmonic
        /// means), which keeps the time waves take through thin slow layers.
        #[arg(long, requires = "smooth")]
        slowness: bool,
        /// Raise vs wherever Poisson's ratio exceeds this (vs = vp/√3 at 0.25), sample by
        /// sample before averaging; fluids stay fluid.
        #[arg(long, value_name = "RATIO")]
        max_poisson: Option<f64>,
        /// The directory to write the model files into.
        #[arg(long)]
        out: PathBuf,
    },
    /// Work out the memory a simulation needs, part by part, without setting it up.
    Memory {
        /// The TOML configuration.
        config: PathBuf,
        /// The memory of `gradient` instead, which keeps checkpoints of the forward run and
        /// runs the adjoint beside it.
        #[arg(long, conflicts_with = "invert")]
        gradient: bool,
        /// The memory of `invert` instead: a gradient run per shot, every shot's observed
        /// traces, and L-BFGS's memory.
        #[arg(long)]
        invert: bool,
        #[command(flatten)]
        strategy: StrategyArgs,
    },
}

/// How a gradient run gets the forward wavefield back.
#[derive(clap::Args)]
struct StrategyArgs {
    /// Rebuild the forward wavefield backwards from the band around the absorbing layers
    /// instead of replaying it from checkpoints: less memory, one more pass, and exact only up
    /// to single-precision round-off. With `memory --invert`, as `[inversion] reconstruct =
    /// true` does.
    #[arg(long)]
    reconstruct: bool,
    /// With --reconstruct: steps per segment, each but the last replayed from a checkpoint to
    /// save its band (by default the length that needs least memory; the run's length or more
    /// saves the band for the whole run and replays nothing).
    #[arg(long, requires = "reconstruct")]
    segment: Option<usize>,
}

impl StrategyArgs {
    fn strategy(&self) -> Strategy {
        if self.reconstruct {
            Strategy::Reconstruct {
                segment: self.segment,
            }
        } else {
            Strategy::Checkpoints { every: None }
        }
    }
}

#[derive(Subcommand)]
enum RemoteAction {
    /// The server: version, threads, memory available.
    Info,
    /// Every run on the server, newest first.
    List,
    /// Follow a run until it finishes, then fetch its results.
    Attach {
        /// The run's id, as `list` shows it.
        id: String,
        /// Where to write the results; out/<id> by default.
        #[arg(long)]
        out: Option<PathBuf>,
        /// Don't report progress.
        #[arg(long, short)]
        quiet: bool,
    },
    /// Fetch a finished run's results now, without waiting.
    Fetch {
        /// The run's id, as `list` shows it.
        id: String,
        /// Where to write the results; out/<id> by default.
        #[arg(long)]
        out: Option<PathBuf>,
    },
    /// Stop a run and delete it with its files.
    Cancel {
        /// The run's id, as `list` shows it.
        id: String,
    },
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e:#}");
            ExitCode::FAILURE
        }
    }
}

fn run() -> Result<()> {
    let cli = Cli::parse();
    if let Some(n) = cli.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()?;
    }
    match cli.command {
        Command::Run {
            config,
            out,
            quiet,
            remote: Some(target),
            port,
            detach,
            ..
        } => {
            let cfg = Config::load(&config)?;
            let out = out.unwrap_or_else(|| cfg.output.dir.clone());
            let r = remote::connect(&target, port)?;
            // Nothing here asks for pictures, so the server runs it headless.
            let s = r.submit(&config)?;
            if let Some(info) = &s.info {
                for w in info["warnings"].as_array().into_iter().flatten() {
                    eprintln!("warning: {}", w.as_str().unwrap_or_default());
                }
            }
            if detach {
                println!("{}", s.id);
                eprintln!(
                    "started {id} on {target}, headless. Progress: `vibe-rations remote {target} list`; \
                     results, once done: `vibe-rations remote {target} fetch {id}`",
                    id = s.id
                );
                return Ok(());
            }
            eprintln!(
                "running {id} on {target}; if this stops, the run carries on, and \
                 `vibe-rations remote {target} attach {id}` picks it up",
                id = s.id
            );
            finish(&r, &s.id, &out, quiet)?;
        }
        Command::Run {
            config,
            out,
            quiet,
            shots,
            ..
        } => {
            let (cfg, sim, _) = load(&config, Kind::Run)?;
            let out = out.unwrap_or_else(|| cfg.output.dir.clone());
            report(&sim);
            if shots {
                // Each shot sets up its own; this one would only hold memory.
                drop(sim);
                for s in 0..cfg.sources.len() {
                    let mut shot = Simulation::new(&cfg.shot(s)?, base(&config))?;
                    let dir = out.join(shot_dir(s));
                    if !quiet {
                        eprintln!("shot {} of {}", s + 1, cfg.sources.len());
                    }
                    shot.run(Some(&dir), !quiet)?;
                    shot.write_outputs(&dir)?;
                }
            } else {
                let mut sim = sim;
                sim.run(Some(&out), !quiet)?;
                sim.write_outputs(&out)?;
            }
            eprintln!("wrote {}", out.display());
        }
        Command::Check { config } => {
            let (_, sim, _) = load(&config, Kind::Run)?;
            report(&sim);
        }
        Command::Serve {
            listen,
            data,
            page,
            log,
        } => {
            server::serve(server::Options {
                listen,
                data,
                page,
                log,
            })?;
        }
        Command::Remote {
            target,
            port,
            action,
        } => {
            let r = remote::connect(&target, port)?;
            match action {
                RemoteAction::Info => println!("{}", serde_json::to_string_pretty(&r.info()?)?),
                RemoteAction::List => {
                    for s in r.runs()? {
                        let done = if s.total_steps > 0 {
                            format!("{:.0}%", 100.0 * s.step as f64 / s.total_steps as f64)
                        } else {
                            "-".into()
                        };
                        let left = match (s.phase, remote::seconds_left(&s)) {
                            (server::Phase::Running, Some(t)) => simulation::human_time(t),
                            _ => "-".into(),
                        };
                        println!(
                            "{:<40} {:<11} {:>5} {:>7}  {}{}",
                            s.id,
                            s.phase.name(),
                            done,
                            left,
                            match (&s.error, s.kind) {
                                (Some(e), _) => e.clone(),
                                (None, server::Kind::Inversion) =>
                                    format!("{} (inversion)", s.name),
                                (None, server::Kind::Run) => s.name.clone(),
                            },
                            if s.headless || s.phase.finished() {
                                ""
                            } else {
                                " (watched)"
                            }
                        );
                    }
                }
                RemoteAction::Fetch { id, out } => {
                    let out = out.unwrap_or_else(|| PathBuf::from("out").join(&id));
                    let n = r.fetch(&id, &out)?;
                    eprintln!("wrote {n} files to {}", out.display());
                }
                RemoteAction::Attach { id, out, quiet } => {
                    let out = out.unwrap_or_else(|| PathBuf::from("out").join(&id));
                    finish(&r, &id, &out, quiet)?;
                }
                RemoteAction::Cancel { id } => {
                    r.delete(&id)?;
                    eprintln!("deleted {id}");
                }
            }
        }
        Command::Memory {
            config,
            gradient,
            invert,
            strategy,
        } => {
            let mut cfg = Config::load(&config)?;
            if invert {
                ensure!(
                    strategy.segment.is_none(),
                    "an inversion rebuilds in segments of the length that needs least memory; \
                     --segment is for --gradient"
                );
                let inversion = cfg
                    .inversion
                    .as_mut()
                    .context("the configuration has no [inversion] section")?;
                inversion.reconstruct |= strategy.reconstruct;
            }
            let kind = match (gradient, invert) {
                (true, _) => Kind::Gradient(strategy.strategy()),
                (_, true) => Kind::Invert,
                _ => Kind::Run,
            };
            let m = plan(&cfg, &config, kind)?;
            println!("{}", m.table());
            // The other way to get the forward fields back, for comparison.
            let other = match kind {
                Kind::Run => None,
                Kind::Gradient(Strategy::Checkpoints { .. }) => Some((
                    Kind::Gradient(Strategy::Reconstruct { segment: None }),
                    "with --reconstruct".to_string(),
                )),
                Kind::Gradient(Strategy::Reconstruct { .. }) => Some((
                    Kind::Gradient(Strategy::Checkpoints { every: None }),
                    "without --reconstruct".to_string(),
                )),
                Kind::Invert => {
                    let inversion = cfg.inversion.as_mut().expect("checked above");
                    inversion.reconstruct = !inversion.reconstruct;
                    Some((
                        Kind::Invert,
                        format!("with reconstruct = {}", inversion.reconstruct),
                    ))
                }
            };
            if let Some((kind, flag)) = other {
                let o = plan(&cfg, &config, kind)?;
                println!("  {:>10}  at most {flag}", memory::human(o.peak));
            }
            if let Some(free) = memory::available() {
                println!("  {:>10}  available now", memory::human(free));
            }
        }
        Command::Extrude {
            vp,
            vs,
            rho,
            input_spacing,
            spacing,
            ny,
            x,
            z,
            speed_scale,
            rho_scale,
            smooth,
            slowness,
            max_poisson,
            out,
        } => {
            let pair = |v: Option<Vec<f64>>, name: &str| -> Result<Option<[f64; 2]>> {
                v.map(|v| match v[..] {
                    [from, to] => Ok([from, to]),
                    _ => anyhow::bail!("--{name} takes two numbers, FROM,TO"),
                })
                .transpose()
            };
            let e = section::Extrusion {
                input_spacing,
                spacing,
                ny,
                x: pair(x, "x")?,
                z: pair(z, "z")?,
                speed_scale,
                rho_scale,
                smooth,
                slowness,
                max_poisson,
            };
            let [vp, vs, rho] = [&vp, &vs, &rho].map(|p| section::read_segy(p));
            let m = section::extrude(&vp?, &vs?, &rho?, &e)?;
            std::fs::create_dir_all(&out).with_context(|| format!("creating {}", out.display()))?;
            let [nx, _, nz] = m.n;
            for (name, values) in [("vp", &m.vp), ("vs", &m.vs), ("rho", &m.rho)] {
                npy::write(&out.join(format!("{name}.npy")), &m.shape, values)?;
            }
            let (what, ny) = match ny {
                Some(ny) => (format!("{nx} × {ny} × {nz} nodes"), ny.to_string()),
                None => (
                    format!("sections of {nx} × {nz} nodes"),
                    "<ny: any>".to_string(),
                ),
            };
            let range = |v: &[f32], solid: bool| {
                let it = v.iter().copied().filter(|&x| !solid || x > 0.0);
                let lo = it.clone().fold(f32::INFINITY, f32::min);
                format!("{lo:.0}–{:.0}", it.fold(0.0, f32::max))
            };
            println!(
                "wrote {}/{{vp,vs,rho}}.npy: {what} {spacing} m apart\n\
                 vp {} m/s, vs {} m/s where solid ({} of {} nodes in a y plane fluid), \
                 rho {} kg/m³\n\
                 [grid]\npoints = [{nx}, {ny}, {nz}]\nspacing = {spacing:?}",
                out.display(),
                range(&m.vp, false),
                range(&m.vs, true),
                m.fluid_nodes,
                nx * nz,
                range(&m.rho, false),
            );
        }
        Command::Invert { config, out } => {
            let (cfg, sim, planned) = load(&config, Kind::Invert)?;
            let out = out.unwrap_or_else(|| cfg.output.dir.clone());
            report(&sim);
            eprintln!(
                "inversion {} at most, with a gradient run per shot and L-BFGS's memory; \
                 `vibe-rations memory --invert` itemises it",
                memory::human(planned.peak)
            );
            warn_next_to_layers(&sim);
            drop(sim);
            eprintln!(
                "the starting misfit and gradient, {} shot(s)",
                cfg.sources.len()
            );
            let mut inversion = Inversion::from_config(&cfg, base(&config))?;
            eprintln!(
                "{} unknowns; misfit {:.6e}",
                inversion.unknowns(),
                inversion.misfit()
            );
            inversion.write_outputs(&out, None)?;
            let stages = inversion.stages().to_vec();
            for (s, stage) in stages.iter().enumerate() {
                if s > 0 {
                    inversion.begin_stage(s)?;
                    // So that inversion.json lists every stage begun, even one that takes no
                    // step.
                    inversion.write_outputs(&out, None)?;
                }
                let start = inversion.misfit();
                if stages.len() > 1 {
                    eprintln!(
                        "stage {}: {} misfit{}, {} iteration(s); misfit {start:.6e}",
                        s + 1,
                        stage.misfit.unwrap_or_default().name(),
                        stage
                            .max_frequency
                            .map_or(String::new(), |f| format!(" below {f} Hz")),
                        stage.iterations
                    );
                }
                for _ in 0..stage.iterations {
                    let Some(it) = inversion.iterate()? else {
                        eprintln!("no step lowers the misfit further; ending the stage");
                        break;
                    };
                    eprintln!(
                        "iteration {:>3}: misfit {:.6e} ({:.1}% of the stage's start), step \
                         {:.2} after {} trial(s), largest change {:.1}%{}",
                        it.iteration,
                        it.misfit,
                        100.0 * it.misfit / start,
                        it.step,
                        it.trials,
                        100.0 * it.largest_change,
                        if it.restarted { ", restarted" } else { "" }
                    );
                    inversion.write_outputs(&out, None)?;
                }
            }
            eprintln!("wrote {}", out.display());
        }
        Command::Gradient {
            config,
            observed,
            out,
            balance,
            strategy,
        } => {
            let strategy = strategy.strategy();
            let (cfg, mut sim, planned) = load(&config, Kind::Gradient(strategy))?;
            let out = out.unwrap_or_else(|| cfg.output.dir.clone());
            report(&sim);
            let (how, flag) = match strategy {
                Strategy::Checkpoints { .. } => ("the checkpoints", ""),
                Strategy::Reconstruct { .. } => ("the band around the layers", " --reconstruct"),
            };
            eprintln!(
                "gradient  {} at most, with {how} and the adjoint run; \
                 `vibe-rations memory --gradient{flag}` itemises it",
                memory::human(planned.peak)
            );
            let next_to_layers = warn_next_to_layers(&sim);
            let r = &sim.receivers;
            let want = [r.points.len(), r.samples];
            let data = r
                .components
                .iter()
                .map(|c| {
                    let path = observed.join(format!("{}.npy", c.name()));
                    let (shape, data) =
                        npy::read(&path).with_context(|| format!("reading {}", path.display()))?;
                    ensure!(
                        shape == want,
                        "{} has shape {shape:?}; the receivers record {want:?} \
                         (observed and synthetic traces must share their sampling: \
                         give both configurations the same [time] dt)",
                        path.display()
                    );
                    Ok(data)
                })
                .collect::<Result<Vec<_>>>()?;
            let weights = if balance {
                balancing_weights(&data, &r.components)
            } else {
                vec![1.0; data.len()]
            };
            eprintln!("running forward, then the adjoint back with the forward run replayed");
            let run = gradient_with(
                &mut sim,
                |r| l2_misfit(&r.traces, &data, &weights),
                strategy,
            )?;
            sim.write_outputs(&out)?;
            let shape = sim.solver.grid.shape();
            let [vp, vs, rho] = run.gradient.velocities(&sim.solver.model);
            for (name, values) in [("vp", vp), ("vs", vs), ("rho", rho)] {
                let single: Vec<f32> = values.iter().map(|&v| v as f32).collect();
                npy::write(&out.join(format!("grad_{name}.npy")), &shape, &single)?;
            }
            // What the run held at its peak, with the observed traces: what the plan counts.
            let held = run.memory_bytes + 4 * data.iter().map(Vec::len).sum::<usize>();
            let summary = serde_json::json!({
                "misfit": run.misfit,
                "weights": weights,
                "checkpoint_every": run.checkpoint_every,
                "strategy": run.strategy,
                "memory_bytes": held,
                "receivers_next_to_layers": next_to_layers,
            });
            std::fs::write(
                out.join("gradient.json"),
                serde_json::to_string_pretty(&summary)? + "\n",
            )?;
            let how = match run.strategy {
                Strategy::Checkpoints { .. } => "checkpoints every",
                Strategy::Reconstruct { .. } => "rebuilt in segments of",
            };
            eprintln!(
                "misfit {:.6e}; {how} {} steps; {} at most",
                run.misfit,
                run.checkpoint_every,
                memory::human(held as u64)
            );
            eprintln!("wrote {}", out.display());
        }
    }
    Ok(())
}

/// Follows a remote run to its end and fetches its results into `out`.
fn finish(r: &remote::Remote, id: &str, out: &Path, quiet: bool) -> Result<()> {
    let s = r.follow(id, !quiet)?;
    ensure!(
        s.phase == server::Phase::Done,
        "the run {}: {}",
        s.phase.name(),
        s.error.unwrap_or_default()
    );
    let n = r.download(id, out)?;
    eprintln!("wrote {n} files to {}", out.display());
    Ok(())
}

fn base(path: &Path) -> &Path {
    path.parent().unwrap_or(Path::new("."))
}

/// What a configuration is loaded for, which decides the memory it needs.
#[derive(Clone, Copy)]
enum Kind {
    Run,
    Gradient(Strategy),
    Invert,
}

/// The memory the run `cfg` (loaded from `path`) needs for `kind`.
fn plan(cfg: &Config, path: &Path, kind: Kind) -> Result<memory::Memory> {
    let m = match kind {
        Kind::Run => memory::estimate(cfg, base(path)),
        Kind::Gradient(strategy) => memory::estimate_gradient_with(cfg, base(path), strategy),
        Kind::Invert => memory::estimate_inversion(cfg, base(path)),
    };
    m.with_context(|| format!("planning {}", path.display()))
}

/// Loads and sets up the simulation `path` describes, once the memory `kind` needs, which it
/// returns too, is known to fit.
fn load(path: &Path, kind: Kind) -> Result<(Config, Simulation, memory::Memory)> {
    let cfg = Config::load(path)?;
    let m = plan(&cfg, path, kind)?;
    if let Some(free) = memory::available() {
        ensure!(
            m.peak <= free,
            "{} needs {} but only {} is available:\n{}",
            path.display(),
            memory::human(m.peak),
            memory::human(free),
            m.table()
        );
    }
    let sim = Simulation::new(&cfg, base(path))
        .with_context(|| format!("setting up {}", path.display()))?;
    Ok((cfg, sim, m))
}

fn report(sim: &Simulation) {
    eprintln!("{}", sim.info.summary());
    for w in &sim.info.warnings {
        eprintln!("warning: {w}");
    }
}

/// Warns about receivers whose adjoint sources would reach into the absorbing layers, and
/// returns how many there are.
fn warn_next_to_layers(sim: &Simulation) -> usize {
    let n = sim.receivers_next_to_layers();
    if n > 0 {
        eprintln!(
            "warning: {n} receiver(s) next to the absorbing layers read velocities inside them, \
             which makes the gradient approximate; move them a node further in"
        );
    }
    n
}
