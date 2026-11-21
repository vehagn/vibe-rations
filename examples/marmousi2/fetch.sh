#!/usr/bin/env bash
# Fetches Marmousi2 (Martin, Wiley & Marfurt 2006): vp, vs and density as SEG-Y, 467 MB, from
# the Madagascar project's mirror (SEG's own link answers 403), into data/marmousi2 unless it
# is there, and checks it against marmousi2.sha256. `just marmousi2` runs it, then convert.sh;
# so does the web app's image build.
set -euo pipefail
cd "$(dirname "$0")/../.."
dir=data/marmousi2
mkdir -p "$dir"
for f in vp vs density; do
  file="$dir/${f}_marmousi-ii.segy"
  [[ -f $file ]] && continue
  echo "downloading $file (156 MB)" >&2
  curl -fsSL --compressed -o "$file.part" "https://www.ahay.org/data/marm2/${f}_marmousi-ii.segy"
  mv "$file.part" "$file"
done
(cd "$dir" && sha256sum --check --quiet) < examples/marmousi2/marmousi2.sha256
