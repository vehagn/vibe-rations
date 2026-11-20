use std::path::{Path, PathBuf};
use std::process::ExitCode;

use anyhow::{Context, Result, ensure};
use clap::{Parser, Subcommand};

use vibe_rations::config::Config;
use vibe_rations::memory;
use vibe_rations::simulation::Simulation;

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
    /// Run a simulation.
    Run {
        /// The TOML configuration.
        config: PathBuf,
        /// Where to write results, instead of the configuration's [output] dir.
        #[arg(long)]
        out: Option<PathBuf>,
        /// Don't report progress.
        #[arg(long, short)]
        quiet: bool,
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
        Command::Run { config, out, quiet } => {
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
