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
(`rust:1.99.1-trixie`) and `rust-toolchain.toml`. To upgrade, bump both and rebuild:

```bash
docker compose -f .devcontainer/compose.yaml up --detach --build --wait
```

## Resources and trouble

- The Colima VM has 4 CPUs and 8 GiB. `just check <config>` prints a run's memory; keep it
  well under 8 GiB. Expect about 200 Mcell updates/s.
- "the input device is not a TTY": add `-T` to `docker compose exec` (`just in` does this).
- A stale build after changing the Containerfile: rebuild as above. The volumes survive;
  `docker compose -f .devcontainer/compose.yaml down --volumes` clears them too.
- Git runs on the host (it has the maintainer's identity); the container is for building.
