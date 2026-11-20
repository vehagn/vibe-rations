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

# The accuracy study: misfits against exact solutions for every stencil order and time step
# (about 1.5 minutes). The README's "Accuracy" section summarises it.
[group("check")]
accuracy:
    cargo test --release --test accuracy -- --ignored --nocapture --test-threads 1

# rustfmt and clippy with -D warnings
[group("check")]
lint:
    cargo fmt --check
    cargo clippy --release --workspace --all-targets -- -D warnings
    cargo clippy --release --package vibe-rations-web --target wasm32-unknown-unknown -- -D warnings

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

# Set up a simulation without running it: grid, time step, accuracy and memory
[group("simulate")]
check config:
    cargo run --release -- check {{ config }}

# Starts `vibe-rations serve` on a remote machine with Docker (or Podman: ENGINE=podman), from
# the published image, if it isn't running already. Its runs and results live in the
# vibe-rations volume there; it listens on that machine's localhost only.
[doc("On the host: start the solver server on a remote machine, e.g. `just remote-up me@box`")]
[group("remote")]
remote-up host image="ghcr.io/vehagn/vibe-rations:main":
    ssh {{ host }} '{{ engine }} start vibe-rations 2>/dev/null || {{ engine }} run -d --name vibe-rations --restart unless-stopped -p 127.0.0.1:8088:8088 -v vibe-rations:/data {{ image }}'

# Replaces the remote server with the image's latest version. Runs in progress there are lost;
# finished runs' results stay in the volume.
[doc("On the host: update the remote solver server to the latest image")]
[group("remote")]
remote-update host image="ghcr.io/vehagn/vibe-rations:main":
    ssh {{ host }} '{{ engine }} pull {{ image }} && {{ engine }} rm -f vibe-rations && {{ engine }} run -d --name vibe-rations --restart unless-stopped -p 127.0.0.1:8088:8088 -v vibe-rations:/data {{ image }}'

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
