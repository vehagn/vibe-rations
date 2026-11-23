# Recipes run inside the dev container, except `dev`, `in`, `serve` and the remote ones
# (`remote-up`, `remote-update`, `tunnel`), which run on the host.
#
# The container engine for `just dev` and `just in`, e.g. ENGINE=podman. Any Compose
# implementation works, e.g. COMPOSE=docker-compose.
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

# Runs a recipe in the dev container from the host, starting the container if need be,
# e.g. `just in test` or `just in run examples/explosion.toml`.
[doc("Run a recipe in the dev container from the host, e.g. `just in ci`")]
[group("dev")]
in +recipe:
    #!/usr/bin/env bash
    set -euo pipefail
    {{ compose }} up --detach --wait >/dev/null 2>&1 || {{ compose }} up --detach --build --wait
    # A terminal only when there is one, so agents and scripts can call this too.
    tty=(); [[ -t 0 && -t 1 ]] || tty=(-T)
    exec {{ compose }} exec "${tty[@]}" dev just {{ recipe }}

# Everything that must pass before a commit: formatting, clippy and all tests.
[group("check")]
ci: lint test

# All tests, optimised: the physics tests run whole simulations (about 20 s).
[group("check")]
test *args:
    cargo test --release --workspace {{ args }}

# Unit tests only: every module's tests (the web crate's too) and the CLI and example checks.
[group("check")]
unit *args:
    cargo test --release --workspace --lib --bins --test cli --test examples {{ args }}

# The physics tests, one at a time, printing their misfits and reflections.
[group("check")]
physics *args:
    cargo test --release --test physics {{ args }} -- --nocapture --test-threads 1

# The accuracy study: misfits against exact solutions for every stencil order and time step,
# Lamb's problem on finer grids, and the free surface's excess by frequency and depth, with and
# without `surface_slope` (about 11 minutes; a test's name runs that one alone, e.g.
# `just accuracy lambs_problem_on_finer_grids`).
# The README's "Accuracy" section summarises it.
[group("check")]
accuracy *name:
    cargo test --release --test accuracy {{ name }} -- --ignored --nocapture --test-threads 1

# The gradient accuracy study: the solver's gradients against the exact sensitivity kernels
# of a homogeneous full space, by grid spacing, order, time step and residual (about 40
# minutes; `just gradient-accuracy <name>` for one). The research log records the results.
[group("check")]
gradient-accuracy *args:
    cargo test --release --test gradient_accuracy {{ args }} -- --ignored --nocapture --test-threads 1

# The interface study: plane waves reflected and transmitted at a planar interface, by where it
# falls between the nodes, how layers are put on them and the stencil order (about 3 minutes).
[group("check")]
interfaces:
    cargo test --release --test interfaces -- --ignored --nocapture --test-threads 1

# The inversion study: misfits, frequency-continuation stages and preconditioners on small
# surveys (about 30 minutes; `just study <name>` for one). The research log records the results.
[group("check")]
study *args:
    cargo test --release --test inversion_study {{ args }} -- --ignored --nocapture --test-threads 1

# The reconstruction study: how far the forward wavefield rebuilt backwards drifts, in the
# field and the gradient, by run length, order, time step and boundary (a few minutes).
[group("check")]
drift:
    cargo test --release --test reconstruction -- --ignored --nocapture --test-threads 1

# Builds the command-line program for macOS on Apple silicon into out/mac/vibe-rations, to run
# natively on the host with all its cores and memory: `out/mac/vibe-rations invert ...`. Cross-
# compiled here, linked by zig; rustc's warning that it can't find the macOS SDK is expected.
[group("build")]
mac:
    MACOSX_DEPLOYMENT_TARGET=11.0 cargo zigbuild --release --target aarch64-apple-darwin --bin vibe-rations
    mkdir -p out/mac
    # macOS kills a signed binary rewritten in place (its code-signing cache keeps the old
    # file's), so the new one goes in as a new file.
    rm -f out/mac/vibe-rations
    cp "${CARGO_TARGET_DIR:-target}/aarch64-apple-darwin/release/vibe-rations" out/mac/

# Builds an integration test (tests/<test>.rs) for macOS on Apple silicon into out/mac/<test>,
# to run natively on the host from the repository's root, where its relative paths lead:
# `out/mac/marmousi2_study west_l2 --ignored --exact --nocapture --test-threads 1` (all cores and
# memory). Not cli, examples or remote, which find files by the container's paths.
[group("build")]
mac-test test:
    #!/usr/bin/env bash
    set -euo pipefail
    exe=$(MACOSX_DEPLOYMENT_TARGET=11.0 cargo zigbuild --release --target aarch64-apple-darwin \
        --test {{ test }} --message-format json-render-diagnostics \
        | jq -r 'select(.reason == "compiler-artifact" and .target.name == "{{ test }}") | .executable // empty' | tail -1)
    [ -n "$exe" ] || { echo "no test binary for {{ test }}" >&2; exit 1; }
    mkdir -p out/mac
    rm -f "out/mac/{{ test }}"   # a new file, as in `mac`
    cp "$exe" "out/mac/{{ test }}"
    echo "out/mac/{{ test }}"

# rustfmt and clippy with -D warnings
[group("check")]
lint:
    cargo fmt --check
    cargo clippy --release --workspace --all-targets -- -D warnings
    cargo clippy --release --package vibe-rations-web --target wasm32-unknown-unknown -- -D warnings
    cargo clippy --release --package vibe-rations --target aarch64-apple-darwin --all-targets -- -D warnings

# Formats the code.
[group("check")]
fmt:
    cargo fmt

# Builds the API documentation into the target directory.
[group("check")]
doc:
    cargo doc --no-deps

# Builds the web app into web/dist: the solver as WebAssembly, the page and the examples.
[group("web")]
web:
    web/build.sh

# Serves web/dist on the host with the production server and headers, at
# http://localhost:8088. Build it first with `just in web`.
[group("web")]
serve port="8088":
    {{ engine }} run --rm --name vibe-rations-web-dev -p 127.0.0.1:{{ port }}:8080 \
      -v "{{ justfile_directory() }}/web/dist:/public:ro" \
      -v "{{ justfile_directory() }}/web/sws.toml:/etc/sws.toml:ro" \
      -e SERVER_CONFIG_FILE=/etc/sws.toml \
      ghcr.io/static-web-server/static-web-server:2.44.3

# Run a simulation, e.g. `just run examples/explosion.toml`
[group("simulate")]
run config *args:
    cargo run --release -- run {{ config }} {{ args }}

# The misfit gradient for vp, vs and rho against observed traces, e.g.
# `just gradient guess.toml out/observed --balance`
[group("simulate")]
gradient config observed *args:
    cargo run --release -- gradient {{ config }} --observed {{ observed }} {{ args }}

# Inverts observed traces for the model with L-BFGS, as the configuration's [inversion] says
[group("simulate")]
invert config *args:
    cargo run --release -- invert {{ config }} {{ args }}

# Fetches Marmousi2 (Martin, Wiley & Marfurt 2006; vp, vs and density as SEG-Y, 467 MB) into
# data/marmousi2 unless it is there, checks it, and makes the model files examples/marmousi2
# and the Marmousi2 examples use: sections of the whole model at 25 and 50 m, and of its
# central 6 km as it is and smoothed (examples/marmousi2/convert.sh lists them)
[group("simulate")]
marmousi2:
    examples/marmousi2/fetch.sh
    examples/marmousi2/convert.sh

# A Marmousi2 inversion study (after `just marmousi2`; tests/marmousi2_study.rs): `misfit_by_arrival`,
# `misfit_by_mute`, `l2`, `l2_rebuilt`, `frequency_stages`, `time_shifts_then_l2`,
# `phase_then_l2`, `stages_with_receivers_side`, `muted_l2`, `muted_stages`,
# `down_weighted_l2`, `down_weighted_stages`, `stages_muted_after_the_first`,
# `windowed_time_shifts_then_l2`, `muted_l2_with_receivers_side`, and on the held-out western
# part `west_l2`, `west_muted_l2`, `west_stages_muted_after_the_first`, each 15-65 minutes,
# `muted_l2_longer` (20 iterations, about 2.5 hours), or `mute_lines` (seconds)
[group("simulate")]
marmousi2-study name:
    cargo test --release --test marmousi2_study {{ name }} -- --ignored --exact --nocapture --test-threads 1

# Times forward runs of Marmousi2 and gradient runs of its central part at orders 4 and 8,
# and measures their memory (after `just marmousi2`; about 10 minutes on four cores)
[group("simulate")]
marmousi2-bench:
    examples/marmousi2/bench.sh

# Set up a simulation without running it: grid, time step, accuracy and memory
[group("simulate")]
check config:
    cargo run --release -- check {{ config }}

# Starts `vibe-rations serve` on a remote machine with Docker (or Podman: ENGINE=podman), from
# the published image, if it isn't running already. Its runs and results live in the
# vibe-rations volume there, or where the machine's /etc/vibe-rations.conf says (`data=...`,
# as `deploy/simulation-host/setup.sh` writes it); it listens on that machine's localhost only.
[doc("On the host: start the solver server on a remote machine, e.g. `just remote-up me@box`")]
[group("remote")]
remote-up host image="ghcr.io/vehagn/vibe-rations:main":
    ssh {{ host }} 'data=vibe-rations; [ -r /etc/vibe-rations.conf ] && . /etc/vibe-rations.conf; {{ engine }} start vibe-rations 2>/dev/null || {{ engine }} run -d --name vibe-rations --restart unless-stopped -p 127.0.0.1:8088:8088 -v "$data":/data {{ image }}'

# Replaces the remote server with the image's latest version. Runs in progress there are lost;
# finished runs' results stay where they were.
[doc("On the host: update the remote solver server to the latest image")]
[group("remote")]
remote-update host image="ghcr.io/vehagn/vibe-rations:main":
    ssh {{ host }} 'data=vibe-rations; [ -r /etc/vibe-rations.conf ] && . /etc/vibe-rations.conf; {{ engine }} pull {{ image }} && {{ engine }} rm -f vibe-rations && {{ engine }} run -d --name vibe-rations --restart unless-stopped -p 127.0.0.1:8088:8088 -v "$data":/data {{ image }}'

# Forwards a local port to the remote server for the web page (and for the CLI in the dev
# container, at http://host.docker.internal:<port>). Stop it with Ctrl-C.
[doc("On the host: tunnel to a remote solver server; open http://localhost:8089")]
[group("remote")]
tunnel host port="8089":
    @echo "open http://localhost:{{ port }} (Ctrl-C to close the tunnel)"
    ssh -N -o ExitOnForwardFailure=yes -L {{ port }}:127.0.0.1:8088 {{ host }}

# The memory a simulation needs, part by part, worked out without setting it up
[group("simulate")]
memory config *args:
    cargo run --release -- memory {{ config }} {{ args }}
