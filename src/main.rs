use std::path::{Path, PathBuf};
use std::process::ExitCode;

use anyhow::{Context, Result, ensure};
use clap::{Parser, Subcommand};

use vibe_rations::config::Config;
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
    /// Work out the memory a simulation needs, part by part, without setting it up.
    Memory {
        /// The TOML configuration.
        config: PathBuf,
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
            config, out, quiet, ..
        } => {
            let (cfg, sim) = load(&config)?;
            let mut sim = sim;
            let out = out.unwrap_or_else(|| cfg.output.dir.clone());
            report(&sim);
            sim.run(Some(&out), !quiet)?;
            sim.write_outputs(&out)?;
            eprintln!("wrote {}", out.display());
        }
        Command::Check { config } => {
            let (_, sim) = load(&config)?;
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
        Command::Memory { config } => {
            let cfg = Config::load(&config)?;
            let m = memory::estimate(&cfg, base(&config))
                .with_context(|| format!("planning {}", config.display()))?;
            println!("{}", m.table());
            if let Some(free) = memory::available() {
                println!("  {:>10}  available now", memory::human(free));
            }
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
fn load(path: &Path) -> Result<(Config, Simulation)> {
    let cfg = Config::load(path)?;
    let m = memory::estimate(&cfg, base(path))
        .with_context(|| format!("planning {}", path.display()))?;
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
    Ok((cfg, sim))
}

fn report(sim: &Simulation) {
    eprintln!("{}", sim.info.summary());
    for w in &sim.info.warnings {
        eprintln!("warning: {w}");
    }
}
