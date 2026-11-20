---
name: remote-runs
description: Run vibe-rations simulations on another machine - `vibe-rations serve` there, `run --remote`/`remote` on the command line, the web page served by the server. Use when starting, following, debugging or changing remote runs, the run server's HTTP API (src/server.rs), its client (src/remote.rs), web/www/remote.js, or the solver image (Containerfile, solver.yaml).
---

# Remote runs

`vibe-rations serve` keeps simulations for clients that come and go: the command line
(`run --remote`, `remote`) and the web page (`web/www/remote.js`). It listens on localhost,
and SSH tunnels reach it. The README's "Remote runs" has the user's view; this is for working
on it.

## Try it locally

```bash
cargo run --release -- serve --data /tmp/runs --page web/dist
```

in one dev-container shell (`just dev`), and in another:

```bash
cargo run --release -- run examples/explosion.toml --remote http://127.0.0.1:8088
cargo run --release -- run examples/explosion.toml --remote http://127.0.0.1:8088 --detach
cargo run --release -- remote http://127.0.0.1:8088 list
cargo run --release -- remote http://127.0.0.1:8088 fetch <id>
```

For the page, the browser needs the port: build the image on the host,

```bash
docker build --file Containerfile --tag vibe-rations-solver:test .
```

and start the `solver` entry of `.claude/launch.json` with `preview_start` (it publishes the
container on 127.0.0.1:8090). **Compute** shows "Server", the Threads line says "on the
server", Run steps it and **Runs** lists it; switching to "This browser" runs the threaded
WebAssembly. Afterwards `docker rm -f vibe-rations-solver`: stopping the preview leaves the
container running. `serve --log` prints every request.

**Use a launch.json entry, not `navigate`, for any local origin.** The built-in browser
blocks workers started from inside a worker on origins it didn't start as a preview server.
Threaded WebAssembly needs them, so on such an origin "This browser" hangs at "Switching the
solver" while everything else works. Python's server and static-web-server on another port
fail the same way; a real browser doesn't.

## Rules

- **Never expose the server.** No authentication: localhost plus SSH is the security model.
  Publish the container's port on 127.0.0.1 only. Anything wider needs auth first.
- **Paths from clients stay inside the run** (`server::relative`): model files, outputs and
  page files. Keep the refusal tests in `tests/remote.rs` passing.
- **The page's two solvers answer the same messages** (`load`, `run`, `pause`, `view`,
  `live`, `npy`, `estimate`; `attach` for the server). Change worker.js and remote.js
  together. A `frame` without a section is progress only (live pictures off).
- **Never hold `Run::shared` while stepping, and take `Run::state` through `Run::state()`**
  (not `.state.lock()`), so requests take turns with the scheduler: std's Mutex isn't fair.
  Status must answer in milliseconds during a run; check it with curl's
  `-w "%{time_starttransfer}"` against a headless run.
- **Headless is inferred, not set:** a run is watched while picture requests (`section`,
  `gather`) keep coming, within `WATCH_WINDOW`. Don't add a flag clients must remember to
  clear.
- **Results must match local runs bit for bit** (`a_remote_run_gives_what_a_local_one_does`).
- A run's directory: `config.toml`, model files at their configured relative paths,
  `status.json`, `out/` (what `write_outputs` writes, plus `snapshots/`).

## Tests

`cargo test --release --test remote` starts real server processes on free ports: remote vs
local traces, model-file uploads, pause/resume, reattach from a second connection, deletion,
refusals, and a restart.
