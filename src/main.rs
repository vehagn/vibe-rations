use std::path::{Path, PathBuf};
use std::process::ExitCode;

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};

use vibe_rations::config::Config;
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
    }
    Ok(())
}

fn load(path: &Path) -> Result<(Config, Simulation)> {
    let cfg = Config::load(path)?;
    let base = path.parent().unwrap_or(Path::new("."));
    let sim =
        Simulation::new(&cfg, base).with_context(|| format!("setting up {}", path.display()))?;
    Ok((cfg, sim))
}

fn report(sim: &Simulation) {
    eprintln!("{}", sim.info.summary());
    for w in &sim.info.warnings {
        eprintln!("warning: {w}");
    }
}
