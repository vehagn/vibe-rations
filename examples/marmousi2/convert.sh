#!/usr/bin/env bash
# Makes Marmousi2's model files from the SEG-Y fetch.sh fetched, with `vibe-rations extrude`:
# sections, (nz, nx), which a run repeats along y, beside the configurations that use them
# (git-ignored):
#
#   25m/                 the whole section at 25 m        marmousi2.toml
#   25m-central/         its central 6 km, x = 6–12 km     central-true.toml, ../marmousi2-central.toml
#   25m-central-smooth/  the same smoothed over 250 m      central.toml (a starting model)
#   50m/                 the whole section at 50 m        ../marmousi2-50m.toml
set -euo pipefail
cd "$(dirname "$0")/../.."
raw=data/marmousi2
here=examples/marmousi2
extrude=(cargo run --release --quiet --bin vibe-rations -- extrude
  "$raw"/{vp,vs,density}_marmousi-ii.segy --input-spacing 1.25 --speed-scale 1000 --rho-scale 1000)
"${extrude[@]}" --spacing 25 --out "$here/25m"
"${extrude[@]}" --spacing 25 --x 6000,12000 --out "$here/25m-central"
"${extrude[@]}" --spacing 25 --x 6000,12000 --smooth 250 --out "$here/25m-central-smooth"
"${extrude[@]}" --spacing 50 --out "$here/50m"
