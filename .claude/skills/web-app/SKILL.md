---
name: web-app
description: Change, build, test or deploy the vibe-rations web app - the solver compiled to WebAssembly and run in the visitor's browser, single- or multithreaded, served as static files from a container in Kubernetes. Use for anything under web/, deploy/ or the web image workflow.
---

# The web app

The visitor's browser does all the computing; the server only serves files. The solver runs
in a module Web Worker (`web/www/worker.js`), which steps it in ~50 ms chunks and posts back
ImageBitmaps of the section and the gather. The page (`app.js`) draws them with axes and
handles the controls.

```
web/src/lib.rs   Simulator (wasm-bindgen) → View (plain Rust; unit-tested natively)
web/www/         index.html, app.js, worker.js, style.css (no bundler, no npm)
web/build.sh     → web/dist: pkg/ (1 thread), pkg-threads/ (rayon on workers), page, examples
web/sws.toml     static-web-server: COOP/COEP headers, no-cache, /health
```

## Build, serve, look

```bash
just in web        # in the dev container: both WebAssembly builds, wasm-opt, web/dist
just serve         # on the host: the production server on http://localhost:8088
```

Open the page with the built-in browser (`.claude/launch.json` has a `web` entry for
`preview_start`). Check:

- the status line says "Solver ready" and **Threads** shows the core count (1 means the page
  isn't cross-origin isolated, or `?threads=1`);
- Run produces wavefronts in the section and a gather filling in; there are no console errors;
- a configuration error (e.g. a misspelt key) shows in the red box and disables Run.

Useful URLs: `?example=layered`, `?threads=1` (to compare speeds).

**Speed numbers only mean something with the browser pane visible.** Hidden, it runs as a
background page and is throttled tenfold or more. Expect ~51 Mcell/s on one thread and
~230 on ten (Apple Silicon). Native timings of `View` methods are a reliable alternative:
drawing a section costs ~2 ms.

Driving it from scripts: `#run` toggles run and pause, `#reset` reloads the configuration,
`#editor` and `#apply` hold and apply the TOML, and `#stat-step`, `#stat-speed`,
`#stat-threads` and `#status` show progress.

## Changing it

- New solver features show up in the web app through `vibe-rations` itself. Expose
  anything new the page needs on `View` first, with a native test in `web/src/lib.rs`, then
  wrap it on `Simulator`.
- Keep the worker protocol small: messages `load`, `run`, `pause`, `view`, `npy` in; `ready`,
  `loaded`, `state`, `frame`, `npy`, `error` out.
- `just in unit` runs the web crate's tests; `just in lint` also runs clippy on the
  `wasm32-unknown-unknown` target.
- Examples in `examples/` appear in the page's menu automatically. Keep them small enough
  for a browser: a few million cells.

## Threads

The threaded build fails in subtle ways; `AGENTS.md` (Notes) lists the four requirements and
what each looks like when it's missing. In short: nightly `-Z build-std` with atomics, explicit
shared-memory and TLS linker flags, the rewritten import in wasm-bindgen-rayon's worker
helper (`build.sh` checks it), and the COOP/COEP headers. Pins live in `web/build.sh` and
both Containerfiles; keep them equal.

## Deploying

- `web/Containerfile` (built from the repository root) produces the image; the build stage
  runs on the build platform and the serving stage only copies files, so multi-arch needs no
  emulation. Check a local build with
  `docker build -f web/Containerfile -t vibe-rations-web:dev .` and run it with
  `--read-only --user 65534:65534 --cap-drop ALL -p 127.0.0.1:8099:8080`.
- `.github/workflows/web.yaml` pushes `ghcr.io/<owner>/vibe-rations-web` (branch, version and
  sha tags) outside pull requests.
- `deploy/kubernetes/` is a Kustomization in the maintainer's homelab style (Gateway API
  `HTTPRoute` to the `external` and `internal` gateways). Set the hostname and image tag;
  it needs no secrets or storage.
- Whatever sits in front must keep the `Cross-Origin-Opener-Policy` and
  `Cross-Origin-Embedder-Policy` headers, or visitors get one thread.
