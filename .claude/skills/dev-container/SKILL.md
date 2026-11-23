---
name: dev-container
description: Run cargo, just, tests or the simulator for vibe-rations inside its dev container from the host, and fix container problems. Use before running any build, test or simulation command in this repository; nothing is installed on the host.
---

# Working in the dev container

All building, testing and running happens in the dev container defined in
`.devcontainer/` (Compose project `vibe-rations-dev`, service `dev`, user `dev`, workspace
`/workspaces/vibe-rations`). The host has Docker (Colima) and `just`, but no Rust: never
install a toolchain on the host. Tools the project needs go in `.devcontainer/Containerfile`.

## Running things

From the host, in the repository root:

```bash
just in unit                         # any recipe: unit, physics, test, lint, ci, fmt, run, check
just in run examples/explosion.toml  # recipes with arguments
```

`just in` starts the container if it isn't running and drops the terminal when there is
none, so it works from agents and scripts. For anything that isn't a recipe:

```bash
docker compose -f .devcontainer/compose.yaml exec -T dev cargo tree --depth 1
```

Run commands that take over two minutes in the background, e.g. long simulations.

## What lives where

| Path in the container | What | Persists |
|---|---|---|
| `/workspaces/vibe-rations` | The repository, bind-mounted: edits on the host show up at once | host |
| `/home/dev/.cache/target` | `CARGO_TARGET_DIR`: build output, including the binary | volume `cargo-target` |
| `/usr/local/cargo/registry` | Downloaded crates | volume `cargo-registry` |
| `/home/dev/.local/state/shell` | Shell history | volume `shell-history` |

Simulation output goes where the configuration's `[output] dir` says, relative to the
working directory: under `out/` in the repository (git-ignored), so it is visible on the
host and you can read the PNGs there.

## Toolchain

Rust is pinned twice and both must agree: the `FROM` line of `.devcontainer/Containerfile`
(`rust:1.99.1-trixie`) and `rust-toolchain.toml`, which also lists the targets
(`wasm32-unknown-unknown`, `aarch64-apple-darwin`). zig (`ZIG_VERSION`, with its Linux
wheels' `ZIG_WHEEL_SHA256`) and `CARGO_ZIGBUILD_VERSION` are pinned in
`.devcontainer/Containerfile` for `just mac`; bump the hashes with the version (PyPI's
`ziglang` JSON lists them). The web app adds three more pins, each in
`.devcontainer/Containerfile` and `web/Containerfile` alike: `WASM_BINDGEN_VERSION` (equal to
the `wasm-bindgen` version in `web/Cargo.toml`), `WASM_THREADS_TOOLCHAIN` (also the default in
`web/build.sh`), and the `rust` base image. To upgrade any, bump every copy and rebuild:

```bash
docker compose -f .devcontainer/compose.yaml up --detach --build --wait
```

## Resources and trouble

- The Colima VM has 8 CPUs and 24 GiB of the Mac's 10 and 32 (`colima start --cpu 8
  --memory 24`, since 2026-11-23; it had 4 and 8). `just check <config>` prints a run's
  memory; keep it well under 23 GiB. Order 8 runs `examples/layered.toml` at 176 Mcells/s
  (156 on 4 CPUs: the solver is limited by memory bandwidth, so more CPUs gain little; the
  memory is what grew).
- Working in a git worktree under the repository (e.g. `.claude/worktrees/<name>`): don't run
  `just in` there. Its compose file mounts `..`, the worktree, and `up` under the shared
  project name would recreate the one container with the worktree mounted. Use the main
  checkout's compose file with `exec -w /workspaces/vibe-rations/<worktree path>`, and a
  `CARGO_TARGET_DIR` of its own so the two trees don't keep rebuilding each other.
- "the input device is not a TTY": add `-T` to `docker compose exec` (`just in` does this).
- A stale build after changing the Containerfile: rebuild as above. The volumes survive;
  `docker compose -f .devcontainer/compose.yaml down --volumes` clears them too.
- Git runs on the host (it has the maintainer's identity); the container is for building.
- Compile in the container; running what it builds on the host is fine and often better,
  since the host has more cores. The web app is one case: `just in web`, then `just serve` on
  the host and open http://localhost:8088. The command line is the other: `just in mac`
  cross-compiles it for macOS on Apple silicon into `out/mac/vibe-rations`, which the host
  runs natively, e.g. `out/mac/vibe-rations invert <config>` (all cores and memory, not the
  VM's 8 and 24 GiB; 193 Mcells/s for the same example). A Linux binary from the container won't run on macOS.
