use std::path::{Path, PathBuf};
use std::process::ExitCode;

use anyhow::{Context, Result, ensure};
use clap::{Parser, Subcommand};

use vibe_rations::config::{Config, shot_dir};
use vibe_rations::gradient::{balancing_weights, gradient, l2_misfit};
use vibe_rations::inversion::Inversion;
use vibe_rations::npy;
use vibe_rations::simulation::{self, Simulation};
use vibe_rations::{memory, remote, server};

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
    },
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
                let base = config.parent().unwrap_or(Path::new("."));
                for s in 0..cfg.sources.len() {
                    let mut shot = Simulation::new(&cfg.shot(s)?, base)?;
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
                            format!("{:?}", s.phase).to_lowercase(),
                            done,
                            left,
                            s.error.as_deref().unwrap_or(&s.name),
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
        } => {
            let cfg = Config::load(&config)?;
            let kind = match (gradient, invert) {
                (true, _) => Kind::Gradient,
                (_, true) => Kind::Invert,
                _ => Kind::Run,
            };
            let m = plan(&cfg, &config, kind)?;
            println!("{}", m.table());
            if let Some(free) = memory::available() {
                println!("  {:>10}  available now", memory::human(free));
            }
        }
        Command::Invert { config, out } => {
            let (cfg, sim, planned) = load(&config, Kind::Invert)?;
            let settings = cfg
                .inversion
                .clone()
                .context("the configuration has no [inversion] section")?;
            let out = out.unwrap_or_else(|| cfg.output.dir.clone());
            report(&sim);
            eprintln!(
                "inversion {} at most, with a gradient run per shot and L-BFGS's memory; \
                 `vibe-rations memory --invert` itemises it",
                memory::human(planned.peak)
            );
            let next_to_layers = sim.receivers_next_to_layers();
            if next_to_layers > 0 {
                eprintln!(
                    "warning: {next_to_layers} receiver(s) next to the absorbing layers read \
                     velocities inside them, which makes the gradient approximate"
                );
            }
            drop(sim);
            let base = config.parent().unwrap_or(Path::new("."));
            eprintln!(
                "the starting misfit and gradient, {} shot(s)",
                cfg.sources.len()
            );
            let mut inversion = Inversion::from_config(&cfg, base)?;
            let first = inversion.misfit();
            eprintln!("{} unknowns; misfit {first:.6e}", inversion.unknowns());
            std::fs::create_dir_all(&out).with_context(|| format!("creating {}", out.display()))?;
            save(&out, &inversion, &settings, first)?;
            let stages = inversion.stages().to_vec();
            for (s, stage) in stages.iter().enumerate() {
                if s > 0 {
                    inversion.begin_stage(s)?;
                    // So that inversion.json lists every stage begun, even one that takes no
                    // step.
                    save(&out, &inversion, &settings, first)?;
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
                    save(&out, &inversion, &settings, first)?;
                }
            }
            eprintln!("wrote {}", out.display());
        }
        Command::Gradient {
            config,
            observed,
            out,
            balance,
        } => {
            let (cfg, mut sim, planned) = load(&config, Kind::Gradient)?;
            let out = out.unwrap_or_else(|| cfg.output.dir.clone());
            report(&sim);
            eprintln!(
                "gradient  {} at most, with the checkpoints and the adjoint run; \
                 `vibe-rations memory --gradient` itemises it",
                memory::human(planned.peak)
            );
            let next_to_layers = sim.receivers_next_to_layers();
            if next_to_layers > 0 {
                eprintln!(
                    "warning: {next_to_layers} receiver(s) next to the absorbing layers read \
                     velocities inside them, which makes the gradient approximate; move them a \
                     node further in"
                );
            }
            let r = &sim.receivers;
            let want = [r.nodes.len(), r.samples];
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
            let run = gradient(&mut sim, |r| l2_misfit(&r.traces, &data, &weights), None)?;
            sim.write_outputs(&out)?;
            let g = sim.solver.grid;
            let shape = [g.n[2], g.n[1], g.n[0]];
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
                "memory_bytes": held,
                "receivers_next_to_layers": next_to_layers,
            });
            std::fs::write(
                out.join("gradient.json"),
                serde_json::to_string_pretty(&summary)? + "\n",
            )?;
            eprintln!(
                "misfit {:.6e}; checkpoints every {} steps; {} at most",
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
        format!("{:?}", s.phase).to_lowercase(),
        s.error.unwrap_or_default()
    );
    let n = r.download(id, out)?;
    eprintln!("wrote {n} files to {}", out.display());
    Ok(())
}

fn base(path: &Path) -> &Path {
    path.parent().unwrap_or(Path::new("."))
}

/// Sets a run up, after making sure it fits in the memory available, if that can be known.
/// What a configuration is loaded for, which decides the memory it needs.
#[derive(Clone, Copy)]
enum Kind {
    Run,
    Gradient,
    Invert,
}

/// The memory the run `cfg` (loaded from `path`) needs for `kind`.
fn plan(cfg: &Config, path: &Path, kind: Kind) -> Result<memory::Memory> {
    let m = match kind {
        Kind::Run => memory::estimate(cfg, base(path)),
        Kind::Gradient => memory::estimate_gradient(cfg, base(path), None),
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

/// Writes the model as it stands and the history so far.
fn save(
    out: &Path,
    inversion: &Inversion,
    settings: &vibe_rations::config::InversionConfig,
    first: f64,
) -> Result<()> {
    let m = &inversion.model;
    let shape = [m.grid.n[2], m.grid.n[1], m.grid.n[0]];
    for (name, values) in [("vp", &m.vp), ("vs", &m.vs), ("rho", &m.rho)] {
        npy::write(&out.join(format!("{name}.npy")), &shape, values)?;
    }
    let summary = serde_json::json!({
        "settings": settings,
        "unknowns": inversion.unknowns(),
        "starting_misfit": first,
        "stages": inversion
            .stages()
            .iter()
            .zip(&inversion.stage_misfits)
            .map(|(stage, misfit)| serde_json::json!({
                "max_frequency": stage.max_frequency,
                "misfit": stage.misfit,
                "iterations": stage.iterations,
                "starting_misfit": misfit,
            }))
            .collect::<Vec<_>>(),
        "iterations": inversion.history,
        "memory_bytes": inversion.memory_bytes,
    });
    std::fs::write(
        out.join("inversion.json"),
        serde_json::to_string_pretty(&summary)? + "\n",
    )?;
    Ok(())
}

fn report(sim: &Simulation) {
    eprintln!("{}", sim.info.summary());
    for w in &sim.info.warnings {
        eprintln!("warning: {w}");
    }
}
