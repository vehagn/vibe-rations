#!/usr/bin/env bash
# Builds the web app into web/dist: the solver as WebAssembly with its JavaScript bindings, the
# page, and the example configurations. Needs the wasm32-unknown-unknown target, the
# wasm-bindgen CLI at the version web/Cargo.toml pins, and, optionally, wasm-opt (binaryen).
# `just web` runs it in the dev container; web/Containerfile runs it to build the image.
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
target="${CARGO_TARGET_DIR:-$root/target}"
out="$root/web/dist"

cargo build --release --manifest-path "$root/Cargo.toml" \
  --package vibe-rations-web --target wasm32-unknown-unknown

# Empty the directory rather than replace it: `just serve` bind-mounts it, and a mount of a
# deleted directory stays empty.
mkdir -p "$out"
find "$out" -mindepth 1 -delete
mkdir -p "$out/pkg" "$out/examples"
wasm-bindgen --target web --no-typescript --out-dir "$out/pkg" \
  "$target/wasm32-unknown-unknown/release/vibe_rations_web.wasm"
if command -v wasm-opt >/dev/null; then
  wasm="$out/pkg/vibe_rations_web_bg.wasm"
  wasm-opt -O3 --enable-simd --enable-bulk-memory --enable-nontrapping-float-to-int \
    --enable-sign-ext --enable-mutable-globals --enable-reference-types --enable-multivalue \
    "$wasm" -o "$wasm"
fi

cp "$root"/web/www/* "$out/"
cp "$root"/examples/*.toml "$out/examples/"
(cd "$out/examples" && printf '[%s]\n' "$(ls -- *.toml | sed 's/.*/"&"/' | paste -sd, -)" > index.json)

echo "built $out ($(du -sh "$out/pkg/vibe_rations_web_bg.wasm" | cut -f1) of WebAssembly)"
