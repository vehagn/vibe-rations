# The container engine for `just dev`, e.g. ENGINE=podman. Any Compose implementation works,
# e.g. COMPOSE=docker-compose.
engine := env("ENGINE", "docker")
compose := env("COMPOSE", engine + " compose") + " --file .devcontainer/compose.yaml"

[private]
default:
    @just --list

# Starts the dev container (building it if need be) and opens a shell in it. `just dev down`
# stops it, keeping the cargo cache and build directory in their volumes.
[doc("Open a shell in the dev container; `just dev down` stops it")]
[group("dev")]
dev command="up":
    #!/usr/bin/env bash
    set -euo pipefail
    case "{{ command }}" in
      up) {{ compose }} up --detach --build --wait; exec {{ compose }} exec dev zsh ;;
      down) exec {{ compose }} down ;;
      *) echo "just dev [up|down]" >&2; exit 1 ;;
    esac

# All tests, optimised: the physics tests run whole simulations.
[group("check")]
test *args:
    cargo test --release {{ args }}

# rustfmt and clippy with -D warnings
[group("check")]
lint:
    cargo fmt --check
    cargo clippy --release --all-targets -- -D warnings

# Run a simulation, e.g. `just run examples/explosion.toml`
[group("simulate")]
run config *args:
    cargo run --release -- run {{ config }} {{ args }}

# Set up a simulation without running it: grid, time step, accuracy and memory
[group("simulate")]
check config:
    cargo run --release -- check {{ config }}
