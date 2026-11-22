#!/usr/bin/env bash
# Downloads, while a simulation machine is on an unmetered line, what it may need later
# (README.md, "Downloaded while unmetered"). Run as root after setup.sh:
#
#   ssh root@host 'bash -s' < deploy/simulation-host/prefetch.sh
set -uo pipefail
DATA=${DATA:-/srv/data}
export DEBIAN_FRONTEND=noninteractive
step() { echo; echo "=== $(date -u +%H:%M:%S) $*"; }

step "apt: updates; rsync for cheap transfers, tmux and htop for long sessions, time"
apt-get update -qq && apt-get -y -qq upgrade && apt-get install -y -qq rsync tmux htop time
[ -e /var/run/reboot-required ] && echo "reboot required: $(tr '\n' ' ' </var/run/reboot-required.pkgs)"

cd "$DATA/vibe-rations"
compose=(docker compose -f .devcontainer/compose.yaml)
"${compose[@]}" up --detach --wait >/dev/null
step "every crate in Cargo.lock, into the dev container's cache"
"${compose[@]}" exec -T dev cargo fetch
step "Marmousi2: the SEG-Y (467 MB) and its model files"
just in marmousi2
step "the whole section smoothed over 250 m, a starting model (examples/marmousi2/whole.toml)"
"${compose[@]}" exec -T dev cargo run --release --quiet --bin vibe-rations -- extrude \
    data/marmousi2/{vp,vs,density}_marmousi-ii.segy --input-spacing 1.25 \
    --speed-scale 1000 --rho-scale 1000 --spacing 25 --smooth 250 \
    --out examples/marmousi2/25m-smooth
step "the command line built for this CPU, in $DATA/bin"
"${compose[@]}" exec -T -e RUSTFLAGS="-C target-cpu=native" \
    -e CARGO_TARGET_DIR=/home/dev/.cache/target/native dev \
    cargo build --release --quiet --bin vibe-rations
install -d "$DATA/bin"
"${compose[@]}" cp dev:/home/dev/.cache/target/native/release/vibe-rations "$DATA/bin/vibe-rations-native"
step "numpy and matplotlib for figures, in out/venv (the dev container's Python)"
install -d -o 1000 -g 1000 out
"${compose[@]}" exec -T dev sh -c 'python3 -m venv out/venv && out/venv/bin/pip install -q numpy matplotlib'
step "the solver image built here once, for its caches (toolchains, Debian base)"
docker build -q -f Containerfile -t vibe-rations-solver:local .
step done
