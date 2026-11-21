---
name: web-app
description: Change, build, test or deploy the vibe-rations web app - the solver compiled to WebAssembly and run in the visitor's browser, single- or multithreaded, served as static files from a container in Kubernetes. Use for anything under web/, deploy/ or the web image workflow.
---

# The web app

The visitor's browser does all the computing; the server only serves files. The solver runs
in a module Web Worker (`web/www/worker.js`), which steps it in ~50 ms chunks and posts back
ImageBitmaps of the section and the gather. The page (`app.js`) draws them with axes and
handles the controls. Served by `vibe-rations serve` instead (the `remote-runs` skill), the
page can also use that server's solver through `remote.js`, which answers the same messages;
the pictures come from `src/view.rs` either way.

The configuration's text (`#editor`) is the source of truth. The **Settings** tab's form
reads it and writes into it through `web/src/settings.rs` (`toml_edit`, so comments and
layout survive), loaded on the page's own thread from the single-threaded `pkg/` whichever
solver runs; edits show an Apply/Revert bar until applied. The **FWI** tab runs an inversion
of synthetic data (`view::InversionView`): the configuration's model is the start, a true
model in TOML makes the observed traces, and `settings::synthetic` pins one time step,
`pml_vp_max` and `optimised_for` for both. The **Command** tab shows the command lines that
run the configuration (and the inversion) without the page, and downloads their files.

```
web/src/lib.rs       Simulator → View, Inverter → InversionView (wasm-bindgen; tested natively)
web/src/settings.rs  read/edit the TOML, the commands, the inversion's pair of configurations
web/www/             index.html, app.js (wiring), worker.js, remote.js, config.js (settings on
                     the main thread), form.js, commands.js, fwi.js, plot.js, util.js, style.css
                     (no bundler, no npm)
web/build.sh         → web/dist: pkg/ (1 thread), pkg-threads/ (rayon on workers), page, examples
web/sws.toml         static-web-server: COOP/COEP headers, no-cache, /health
```

## Build, serve, look

```bash
just in web        # in the dev container: both WebAssembly builds, wasm-opt, web/dist
just serve         # on the host: the production server on http://localhost:8088
```

Open the page with the built-in browser (`.claude/launch.json` has a `web` entry for
`preview_start`, and `web-running` to attach when `just serve` is already up: its container
name is taken then). Don't `navigate` to it instead: on an origin it didn't open through
`launch.json` the pane blocks the solver's workers, and the page hangs on "Loading the
solver…". Check:

- served by `just serve`, the Compute switch and the Runs list are hidden (they belong to
  `vibe-rations serve`);
- the status line says "Solver ready" and **Threads** shows the core count (1 means the page
  isn't cross-origin isolated, or `?threads=1`);
- Run produces wavefronts in the section and a gather filling in; there are no console errors;
- a configuration error (e.g. a misspelt key) shows in the red box and disables Run;
- the memory line under the example follows the configuration as it is edited (the worker's
  `estimate` handler, `Simulator.estimateMemory`), and a grid over the 4 GiB WebAssembly can
  address (e.g. `points = [801, 801, 535]`) turns it red and is refused on Apply, before
  anything is allocated (`View::new`, `BROWSER_LIMIT`);
- **Settings:** changing a control rewrites the TOML (comments kept), shows Apply/Revert and
  moves the memory line; a TOML that doesn't parse disables the form with the reason;
- **FWI:** Gradient, then Invert, on the explosion example with the default true model (a
  sphere): the stage line counts shots, the start and iterations, the history table and the
  misfit chart fill in, and the inversion panel draws the update, the true update, the model
  and −gradient for vp, vs, ρ on any plane; Stop waits for the unit in progress;
- **Command:** the commands name the file and machine typed there; Copy copies; the downloads
  give the configuration, and the inversion's start and truth.

Useful URLs: `?example=layered`, `?threads=1` (to compare speeds).

**Speed numbers only mean something with the browser pane visible.** Hidden, it runs as a
background page and is throttled tenfold or more. Expect ~51 Mcell/s on one thread and
~230 on ten (Apple Silicon). Native timings of `View` methods are a reliable alternative:
drawing a section costs ~2 ms.

Driving it from scripts: `#run` toggles run and pause, `#reset` reloads the configuration,
`#editor` and `#apply` hold and apply the TOML (dispatch an `input` event on `#editor` after
setting it, so the form and memory follow), the tabs are `#tab-run`, `#tab-settings`,
`#tab-fwi` and `#tab-command`, `#fwi-gradient`, `#fwi-run` and `#fwi-stop` drive the
inversion, and `#stat-step`, `#stat-speed`, `#stat-threads`, `#status` and `#fwi-stage` show
progress.

## Changing it

- New solver features show up in the web app through `vibe-rations` itself. Expose
  anything new the page needs on `View` first, with a native test in `web/src/lib.rs`, then
  wrap it on `Simulator`.
- Keep the worker protocol small: messages `load`, `run`, `pause`, `view`, `live`, `npy`,
  `estimate` (with `kind: "inversion"` for an inversion's), `invert`, `stopInversion`,
  `resumeInversion`, `inversionView`, `modelNpy` in (`attach` for the server); `ready`,
  `loaded`, `state`, `frame`, `estimate`, `npy`, `inversion`, `inversionFrame`, `modelNpy`,
  `error` out (`inversionAttached` from the server). worker.js and remote.js answer the same.
- Logic worth testing goes in Rust: `settings.rs` (natively tested) for anything about the
  configuration's text or the commands, `view.rs` for pictures. The JavaScript only wires.
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
