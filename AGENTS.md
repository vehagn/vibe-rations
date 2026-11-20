# vibe-rations, for coding agents

A 3D elastic wave simulator (fourth-order staggered-grid velocity–stress finite differences,
CPML, free surface), built as the forward engine for elastic full waveform inversion (FWI).
The [README](README.md) explains the numerics and how they are verified; this file is what
an agent needs to work on it, and the notes below record why things are the way they are.

## Commands

Everything runs in the dev container; nothing is installed on the host. From the host,
`just in <recipe>` runs a recipe in the container (starting it if need be); `just dev`
opens a shell there instead.

| Command | Does | Time |
|---|---|---|
| `just unit` | Module unit tests, the CLI tests and the examples check | seconds |
| `just physics` | Whole simulations against analytic answers, misfits printed | ~20 s |
| `just test` | Everything | ~35 s |
| `just accuracy` | The stencil order and time step study (prints a table) | ~1.5 min |
| `just lint` | rustfmt check and clippy with `-D warnings` | |
| `just ci` | `lint` then `test`: what must pass before every commit | |
| `just fmt` | Format the code | |
| `just check <config>` / `just run <config>` | Set up / run a simulation | |
| `just web` | Build the web app (both WebAssembly variants) into `web/dist` | ~1 min cold |
| `just serve` | On the host: serve `web/dist` at http://localhost:8088 like production | |

## Skills

Task guides in `.claude/skills/*/SKILL.md`; any agent can read them:

| Skill | When |
|---|---|
| `dev-container` | Running cargo, just or the simulator; container trouble |
| `run-simulation` | Writing a configuration, sizing a grid, reading the outputs |
| `verify-solver` | Changing the solver, layers, free surface, sources or receivers; adding a physics test |
| `debug-instability` | A run blows up or doesn't drain |
| `conventional-commit` | Committing work |
| `fwi-literature` | Looking up the theory: adjoints, kernels, misfits, PML, free surface |
| `web-app` | Changing, building, testing or deploying the browser front end |

## Map

```
src/
  main.rs        the CLI: run | check
  config.rs      the TOML configuration (serde, deny_unknown_fields)
  grid.rs        the grid and its two-node halo, indexing
  model.rs       vp, vs, rho on the nodes: layers, .npy files, spheres
  solver.rs      fields, staggered material, the update kernels, free surface, PointOp
  fd.rs          the staggered Stencil (Taylor, order 2–16), a row at a time
  cpml.rs        CPML profiles and memory slabs, the sponge factors
  source.rs      wavelets, moment tensors, how sources enter the grid
  receiver.rs    recording traces through PointOps
  snapshot.rs    sections and volumes during a run
  simulation.rs  a run set up from a configuration; outputs
  image.rs, npy.rs
web/
  src/lib.rs     wasm-bindgen wrapper (Simulator) over a plain-Rust View with native tests
  www/           the page: index.html, app.js (UI), worker.js (runs the solver), style.css
  build.sh       builds web/dist: single-threaded and threaded WebAssembly, page, examples
  sws.toml       static-web-server settings, with the cross-origin isolation headers
  Containerfile  the hosted image
deploy/kubernetes/  Deployment, Service, HTTPRoute (Kustomize)
.github/workflows/  ci.yaml (`just ci`), web.yaml (builds and pushes the image)
tests/
  physics.rs     whole simulations: analytic solutions, layers, reciprocity, stability
  accuracy.rs    the stencil-order study (ignored; `just accuracy`)
  common/        closed-form solutions and setup helpers the two share
  cli.rs         the binary, run as a process
  examples.rs    every example sets up without warnings
examples/        annotated configurations
```

## Conventions

- **Commits:** conventional commits (`feat:`, `fix:`, `test:`, `docs:`, `build:`,
  `refactor:`, `perf:`, `chore:`), unsigned (`git commit --no-gpg-sign`: signing is on in
  the maintainer's global config), one logical change each, and `just ci` green before each.
  Commit as you go, on a branch rather than `main`. The `conventional-commit` skill has the
  details.
- **Tests:** unit tests sit next to the code in `#[cfg(test)] mod tests`; anything that
  runs a whole simulation goes in `tests/physics.rs`. New behaviour gets a test. Don't
  loosen a physics tolerance to make a change pass; the baselines are below.
- **Sources are the transpose of receivers** (`Solver::stagger_weights`, `PointOp`).
  Adjoint FWI depends on it; the reciprocity test checks it to 1e-6.
- SI units; z points down; arrays are x-fastest, padded by `grid::HALO` = 2 nodes.
- Comments say why, not what. Doc comments on public items.

## Notes

Decisions and lessons, newest first. Add to them when you learn something the code
doesn't show.

- **Stencil order and the time step go together** (`just accuracy`, README "Accuracy").
  Order 4 at courant 0.85 benefits from cancelling time and space errors; orders ≥ 8 at
  0.85 are *worse* than order 4, and only beat it at courant ≈ 0.3–0.5. Order 8 at 0.3 needs
  ~4.5 points per S wavelength for ~1% after ten wavelengths, where order 4 needs ~8–9:
  about 3.4× less compute and 5.5× less memory. The halo widens to the stencil's
  half-length; free-surface mirroring stays self-adjoint at every order (reciprocity checked
  at orders 4 and 8). Time–space-optimised coefficients looked good for one speed in a
  von Neumann analysis and failed for other speeds, so they were not adopted.
- **Threads in the browser** need four things, each of which failed once: (1) a std built
  with atomics: nightly `-Z build-std=panic_abort,std`, pinned in `web/build.sh` and both
  Containerfiles; (2) explicit linker flags for an imported, shared memory and the TLS
  exports (`--shared-memory --import-memory --max-memory --export=__wasm_init_tls …`).
  rustc doesn't add them, and without them initThreadPool fails with "Memory could not be
  cloned"; (3) wasm-bindgen-rayon's worker helper imports `'../../..'`, which only bundlers
  resolve, so `build.sh` rewrites it to the module file; (4) a cross-origin-isolated page
  (COOP/COEP in `web/sws.toml`). Without isolation the page falls back to the
  single-threaded build.
- **Browser speeds:** 51 Mcell/s on one thread with SIMD128 (about native per core),
  234 Mcell/s on 10 threads. The built-in browser pane runs as a background page while it
  is hidden and gets throttled (down to 6–7 Mcell/s), so only measure speed with it visible.
- **Web rendering:** sections are interpolated (bilinear, on values) up to ~600 px before
  colouring, because interpolating colours left dotted white zero crossings. The colour
  scale never drops below 5% of the largest it has had this run, so the picture fades as
  energy leaves instead of blowing up the residue. Drawing a section costs ~2 ms.
- **`web/dist` is emptied, not deleted**, on each build: `just serve` bind-mounts it, and a
  mount of a deleted directory serves 404s until the server restarts.
- **Compute:** compile in the dev container; running what it builds on the host is fine (the
  maintainer wants the extra cores), e.g. the web app in a browser on the host.
- **Physics test baselines** (fail if they get worse; update this list if they get better):
  explosion vs analytic 0.7–1.6%; point force vs Stokes 0.3–0.8%; free-surface P doubling
  1.9%; CPML reflection 0.015% of peak without a free surface, 0.18% along one; reciprocity
  ≤ 1e-6; the 9 s soft-sediment run drains to 4e-6 of its peak.
- **Snapshots:** colour ∝ √amplitude, scaled by the 99.5th percentile of |value|. A
  Gaussian moment rate leaves a static stress at the source that swamps pressure (`p`)
  sections; velocity has no static part, so the examples snapshot velocity. In a soft
  near-surface layer the surface amplitudes dominate a vertical section; that is physics.
- **The CPML instability along a free surface.** With a free surface, a soft layer
  (vp/vs ≈ 3.3) and CPML side layers, the field grew exponentially from about 2 s,
  starting in the layer corners at the surface (Festa et al. 2005). The free surface alone
  (no PML) was stable; larger α, κ_max = 4 and a smaller dt didn't help. The fix is a plain
  damping exp(-p·d·x³·dt) in the x and y layers, only when there is a free surface:
  `pml_sponge`, default 0.2. Measured trade-off (reflection at 9 nodes from the layer /
  stable?): d·x⁰ with p 0.1: 0.23% / yes; d·x² 0.1: 0.06–0.15% / yes; d·x⁴ 0.05: — / no;
  d·x³ 0.2 on side layers only: 0.07–0.18% / yes.
- **Sources on the free surface act on half a cell.** vx, vy, σxx, σyy, σzz and σxy on
  the surface plane stand for half a cell (their other half is the mirror image), so their
  source weights are doubled (`Solver::cell_share`). Without it, reciprocity with one end on
  the surface was off by exactly 2.
- **The free surface mirrors velocities evenly.** First version: second-order derivatives
  next to the surface. They were not the adjoint of the stress-imaging operators, energy
  wasn't conserved, and the free-surface misfit was 7.6%. Mirroring vx, vy and vz evenly
  (σ oddly) makes the pair adjoint, keeps fourth order everywhere, and brought it to 1.9%.
  The textbook image method (Fichtner 2011 §3.2.3; Graves 1996) instead derives the
  velocities above the surface from the boundary conditions. It hasn't been tried here; any
  alternative must pass the reciprocity and long-run tests.
- **Cubic interpolation for receivers and sources.** Averaging two staggered values onto a
  node costs ~cos(kh/2), about 2% at 8 points per wavelength; cubic weights
  (-1, 9, 9, -1)/16 cut the misfits roughly threefold.
- **Performance:** about 200–230 Mcell updates/s on the container's 4 cores (Apple
  Silicon, Colima VM with 8 GiB). Memory is 17 padded f32 arrays plus CPML memory: the
  examples need 0.1–0.25 GiB. The work splits by z-plane, each plane a row at a time
  through nine scratch rows.
- **Next:** elastic FWI; the README's roadmap lists the steps.
