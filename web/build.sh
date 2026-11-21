#!/usr/bin/env bash
# Builds the web app into web/dist: the solver as WebAssembly with its JavaScript bindings, in a
# single-threaded and a threaded variant, the page, and the example configurations. Needs the
# wasm32-unknown-unknown target, the pinned nightly below with rust-src, the wasm-bindgen CLI
# at the version web/Cargo.toml pins, and, optionally, wasm-opt (binaryen).
# `just web` runs it in the dev container; web/Containerfile runs it to build the image.
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
target="${CARGO_TARGET_DIR:-$root/target}"
out="$root/web/dist"
# The threaded build needs a std built with atomics, which takes a nightly toolchain (with
# rust-src). Keep in step with the dev container's and web/Containerfile's.
nightly="${WASM_THREADS_TOOLCHAIN:-nightly-2026-11-18}"

# One thread, on the pinned stable toolchain (SIMD from .cargo/config.toml).
cargo build --release --manifest-path "$root/Cargo.toml" \
  --package vibe-rations-web --target wasm32-unknown-unknown

# Rayon on Web Workers: std rebuilt with atomics, and an imported, shared memory (up to the
# 4 GiB wasm32 allows) with the thread-local storage exports wasm-bindgen's thread support
# needs; rustc no longer adds these by itself. A target directory of its own keeps the two
# builds from rebuilding each other.
threads_flags=(
  -C target-feature=+atomics,+bulk-memory,+mutable-globals,+simd128
  -C link-arg=--shared-memory -C link-arg=--import-memory -C link-arg=--max-memory=4294967296
  -C link-arg=--export=__wasm_init_tls -C link-arg=--export=__tls_size
  -C link-arg=--export=__tls_align -C link-arg=--export=__tls_base
)
RUSTFLAGS="${threads_flags[*]}" \
  cargo "+$nightly" build --release --manifest-path "$root/Cargo.toml" \
  --package vibe-rations-web --features threads --target wasm32-unknown-unknown \
  -Z build-std=panic_abort,std --target-dir "$target/threads"

# Empty the directory rather than replace it: `just serve` bind-mounts it, and a mount of a
# deleted directory stays empty.
mkdir -p "$out"
find "$out" -mindepth 1 -delete
mkdir -p "$out/examples"

features=(--enable-simd --enable-bulk-memory --enable-nontrapping-float-to-int --enable-sign-ext
  --enable-mutable-globals --enable-reference-types --enable-multivalue)
bindgen() { # <directory in dist> <wasm> [wasm-opt flags]
  local dir="$out/$1" wasm="$2"
  shift 2
  wasm-bindgen --target web --no-typescript --out-dir "$dir" "$wasm"
  if command -v wasm-opt >/dev/null; then
    wasm-opt -O3 "${features[@]}" "$@" "$dir/vibe_rations_web_bg.wasm" -o "$dir/vibe_rations_web_bg.wasm"
  fi
}
bindgen pkg "$target/wasm32-unknown-unknown/release/vibe_rations_web.wasm"
bindgen pkg-threads "$target/threads/wasm32-unknown-unknown/release/vibe_rations_web.wasm" --enable-threads

# wasm-bindgen-rayon's worker helper imports its package as a directory, which only bundlers
# resolve; browsers need the module's file.
helper=("$out"/pkg-threads/snippets/wasm-bindgen-rayon-*/src/workerHelpers.js)
sed -i.orig "s#import('../../..')#import('../../../vibe_rations_web.js')#" "${helper[@]}"
rm "${helper[@]/%/.orig}"
if grep -q "import('../../..')" "${helper[@]}"; then
  echo "couldn't point wasm-bindgen-rayon's worker helper at the module" >&2
  exit 1
fi

cp "$root"/web/www/* "$out/"
# The examples, each with the model files it names, at the same paths relative to it, where
# the page fetches them. One whose files aren't there (Marmousi2's, before `just marmousi2`)
# is left out.
listed=()
for toml in "$root"/examples/*.toml; do
  name=$(basename "$toml")
  files=$(awk '/^\[/ { section = $1 }
    section == "[model.files]" && /^ *(vp|vs|rho) *=/ {
      match($0, /"[^"]*"/); print substr($0, RSTART + 1, RLENGTH - 2) }' "$toml")
  missing=""
  for f in $files; do
    [[ $f != *..* && -f "$root/examples/$f" ]] || missing+=" $f"
  done
  if [[ -n $missing ]]; then
    echo "leaving out examples/$name: no$missing" >&2
    continue
  fi
  cp "$toml" "$out/examples/"
  for f in $files; do
    mkdir -p "$(dirname "$out/examples/$f")"
    cp "$root/examples/$f" "$out/examples/$f"
  done
  listed+=("\"$name\"")
done
(IFS=,; printf '[%s]\n' "${listed[*]}") > "$out/examples/index.json"

echo "built $out: $(du -h "$out/pkg/vibe_rations_web_bg.wasm" | cut -f1) single-threaded," \
  "$(du -h "$out/pkg-threads/vibe_rations_web_bg.wasm" | cut -f1) threaded"
