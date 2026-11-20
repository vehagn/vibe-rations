use std::path::{Path, PathBuf};
use std::process::ExitCode;

use anyhow::{Context, Result, ensure};
use clap::{Parser, Subcommand};

use vibe_rations::config::Config;
use vibe_rations::gradient::{balancing_weights, gradient, l2_misfit};
use vibe_rations::npy;
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
        Command::Gradient {
            config,
            observed,
            out,
            balance,
        } => {
            let (cfg, mut sim) = load(&config)?;
            let out = out.unwrap_or_else(|| cfg.output.dir.clone());
            report(&sim);
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
                balancing_weights(&data)
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
            let summary = serde_json::json!({
                "misfit": run.misfit,
                "weights": weights,
                "checkpoint_every": run.checkpoint_every,
                "memory_bytes": run.memory_bytes,
                "receivers_next_to_layers": next_to_layers,
            });
            std::fs::write(
                out.join("gradient.json"),
                serde_json::to_string_pretty(&summary)? + "\n",
            )?;
            eprintln!(
                "misfit {:.6e}; checkpoints every {} steps; {:.2} GiB at most",
                run.misfit,
                run.checkpoint_every,
                run.memory_bytes as f64 / (1u64 << 30) as f64
            );
            eprintln!("wrote {}", out.display());
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
