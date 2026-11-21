#!/usr/bin/env bash
# The Marmousi2 benchmark (`just marmousi2-bench`, in the dev container): forward runs of the
# whole section (marmousi2.toml), and forward and gradient runs of its central 6 km
# (central.toml, against central-true.toml's traces), at orders 4 and 8, each timed, with
# its peak memory measured beside what `memory` planned. Needs the model files `just
# marmousi2` makes. Results: docs/research-log.md, "Marmousi2".
set -euo pipefail
cd "$(dirname "$0")/../.."
here=examples/marmousi2
out=out/marmousi2-bench
mkdir -p "$out"
[[ -f $here/25m/vp.npy && -f $here/25m-central-smooth/vp.npy ]] || {
  echo "no model files: run \`just marmousi2\` first" >&2
  exit 1
}
cargo build --release --quiet
bin="${CARGO_TARGET_DIR:-target}/release/vibe-rations"

# Runs a command with its output in a log, then adds its wall time and peak resident memory.
measure() { # log command...
  local log=$1
  shift
  python3 -c '
import resource, subprocess, sys, time
t = time.monotonic()
code = subprocess.run(sys.argv[1:]).returncode
wall = time.monotonic() - t
peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024
print(f"wall {wall:.1f} s, peak {peak:.0f} MiB", file=sys.stderr)
sys.exit(code)' "$@" 2> "$log" || { cat "$log" >&2; exit 1; }
}

# The configuration at another order and without snapshots, its model files found from
# out/marmousi2-bench, which lies as deep below the repository as examples/marmousi2.
variant() { # config order name
  sed -e '/^\[\[snapshot\]\]/,/^$/d' -e "s/^order = 8/order = $2/" \
    -e "s#^\(vp\|vs\|rho\) = \"#&../../$here/#" "$1" > "$out/$3.toml"
  grep -q "^order = $2" "$out/$3.toml"
  echo "$out/$3.toml"
}

summary=()
# A forward run of `cfg`, timed, into the summary as `name`.
forward() { # cfg name
  local planned log=$out/$2.log
  planned=$("$bin" memory "$1" | awk '/at most/ { print $1, $2 }')
  echo "$2" >&2
  measure "$log" "$bin" run "$1" --out "$out/$2"
  summary+=("$2: $(grep -m1 -oE '[0-9]+ steps' "$log"), \
$(grep -E '^step +([0-9]+)/\1 ' "$log" | grep -oE '[0-9]+ Mcells/s'), $(tail -1 "$log"), \
planned $planned")
}

for order in 8 4; do
  forward "$(variant $here/marmousi2.toml $order forward-$order)" forward-$order
done

echo "observed data" >&2
"$bin" run $here/central-true.toml --quiet --out out/marmousi2/central-true
for order in 8 4; do
  cfg=$(variant $here/central.toml $order gradient-$order)
  forward "$cfg" central-forward-$order
  planned=$("$bin" memory "$cfg" --gradient | awk '/at most/ { print $1, $2 }')
  echo "gradient-$order" >&2
  log="$out/gradient-$order.log"
  measure "$log" "$bin" gradient "$cfg" --observed out/marmousi2/central-true \
    --out "$out/gradient-$order"
  summary+=("gradient-$order: $(grep -m1 -oE '[0-9]+ steps' "$log"), \
$(grep -oE 'checkpoints every [0-9]+ steps' "$log"), $(tail -1 "$log"), planned $planned")
done
printf '%s\n' "${summary[@]}"
