#!/bin/bash
# Runs an inversion, every line of its log timestamped (UTC), keeping each iteration's model:
#
#   docs/analysis/run_timed.sh DIR vibe-rations invert <config> --out DIR/out [--gpu]
#
# writes DIR/log, DIR/iterations/<n>/ (the model and inversion.json after iteration n, 0 the
# start) and DIR/done, which docs/analysis/compare_inversions.py reads.
set -u
d=$1; shift
mkdir -p "$d/iterations"
echo "$(date -u +%FT%T.%3NZ) start: $*" > "$d/log"
"$@" 2>&1 | while IFS= read -r line; do
  printf '%s %s\n' "$(date -u +%FT%T.%3NZ)" "$line"
  case $line in
    *unknowns\;*) n=0 ;;
    iteration*) n=$(echo "$line" | awk '{print $2}' | tr -d :) ;;
    *) n= ;;
  esac
  # invert writes the model and inversion.json right after it prints the line.
  if [ -n "$n" ]; then
    (sleep 5; mkdir -p "$d/iterations/$n"; cp "$d"/out/*.npy "$d/out/inversion.json" "$d/iterations/$n/") &
  fi
done >> "$d/log"
echo "$(date -u +%FT%T.%3NZ) exit ${PIPESTATUS[0]}" >> "$d/log"
touch "$d/done"
