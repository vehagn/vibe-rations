# vibe-rations, for coding agents

A 3D elastic wave simulator (staggered-grid velocity–stress finite differences, eighth order
in space by default, CPML, free surface), built as the forward engine for elastic full
waveform inversion (FWI) and as the basis of a paper. The [README](README.md) explains the
numerics and how they are verified; [docs/research-log.md](docs/research-log.md) records the
experiments; this file is what an agent needs to work on it, and the notes below record why
things are the way they are.

## Commands

Everything is built in the dev container; nothing is installed on the host. From the host,
`just in <recipe>` runs a recipe in the container (starting it if need be); `just dev`
opens a shell there instead.

| Command | Does | Time |
|---|---|---|
| `just unit` | Module unit tests, the CLI tests and the examples check | seconds |
| `just physics` | Whole simulations against analytic answers, misfits printed | ~1.5 min |
| `just test` | Everything | ~1.5 min |
| `just accuracy` | The stencil order and time step study (prints a table) | ~2 min |
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
| `fwi-literature` | Looking up the theory: adjoints, kernels, misfits, PML, free surface; citing it |
| `web-app` | Changing, building, testing or deploying the browser front end |

## Map

```
src/
  main.rs        the CLI: run | check
  config.rs      the TOML configuration (serde, deny_unknown_fields); default_courant
  grid.rs        the grid and its halo (as wide as the stencil), indexing
  model.rs       vp, vs, rho on the nodes: layers, .npy files, spheres
  solver.rs      fields, staggered material, the update kernels, free surface, PointOp
  fd.rs          the staggered Stencil (Taylor, order 2–16), a row at a time
  cpml.rs        CPML profiles and memory slabs, the sponge factors
  source.rs      wavelets, moment tensors, how sources enter the grid
  tdt.rs         the time-dispersion transform: pre-warped sources, un-warped traces
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
docs/
  references.bib every work cited anywhere in the repository
  research-log.md experiments, results and decisions, for the paper
  analysis/      von Neumann analysis (pure Python): dispersion, time-space, time dispersion
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
- **Citations:** every work cited anywhere (code, docs, commits, skills) gets an entry in
  `docs/references.bib`, checked against the paper or its Crossref record, with a `note` for
  anything unchecked. Prose cites by author and year; the research log by key, `[@key]`.
- **Research log:** every experiment that informs a decision goes in
  `docs/research-log.md`: question, method, results (with numbers), decision, and how to
  reproduce it. That is the raw material for the paper.
- **Math** in Markdown is LaTeX: `$…$` inline, `$$…$$` displayed. Code comments may use
  Unicode, since rustdoc doesn't render LaTeX.
- **Tests:** unit tests sit next to the code in `#[cfg(test)] mod tests`; anything that
  runs a whole simulation goes in `tests/physics.rs`. New behaviour gets a test. Don't
  loosen a physics tolerance to make a change pass. When a deliberate change trades one
  baseline for another, record both numbers in the test's doc comment, here and in the log.
- **Sources are the transpose of receivers** (`Solver::stagger_weights`, `PointOp`).
  Adjoint FWI depends on it; the reciprocity test checks it to $10^{-6}$.
- SI units; $z$ points down; arrays are $x$-fastest, padded by `grid.halo` nodes.
- Comments say why, not what. Doc comments on public items.

## Notes

Decisions and lessons, newest first. Add to them when you learn something the code
doesn't show.

- **Leapfrog's time error can be removed** (`[scheme] remove_time_dispersion`, `src/tdt.rs`;
  Wang & Xu 2015, Koene et al. 2018). It is a warp of the frequency axis, the same for every
  medium and stencil, so pre-warping the sources and un-warping the traces leaves only the
  spatial error, and $C = 0.85$ costs no accuracy. Order 8 then matches $C = 0.3$ in 2.8 times
  fewer steps; orders 12 and 16, time-limited at 0.3, become 8–12 times more accurate. It is
  opt-in, not the default, until the adjoint run applies the transposed transform to the
  residuals (`Simulation::trace_warp`, `Warp::apply_transpose`). Timing matters: forces enter
  at $n\Delta t$, moments at $(n + \frac12)\Delta t$, velocities are read at
  $(n + \frac12)\Delta t$ and stresses at $(n + 1)\Delta t$; the warp is about $t = 0$.
  Snapshots stay warped, and the layers are only approximately transformed.
- **A research brief** (another agent's literature review and von Neumann analysis,
  2026-11-20) was checked: its numbers reproduce (`docs/analysis/time_dispersion.py`), the
  transform and the stencil-matched interpolation came from it, and the log's entry lists
  what remains, in its open questions. Its references came from memory: check any on
  Crossref before citing it.
- **Error propagation into adjoint fields** is a main research interest for what comes next:
  forward-field error (dispersion, interfaces, free surface, layers), reconstruction error,
  and their effect on gradients. Design the adjoint so each can be measured separately (see
  the log's open questions).
- **Physics test baselines at the default, order 8 with $C = 0.3$** (fail if they get worse;
  update this list if they get better): explosion vs analytic 0.07–0.19%; point force vs
  Stokes 0.03–0.09% (0.001–0.008% at $C = 0.85$ with `remove_time_dispersion`); free-surface
  P doubling 2.4%; CPML reflection $3.9 \times 10^{-4}$ of the peak without a free surface,
  $2.1 \times 10^{-3}$ along one; reciprocity $\le 1.9 \times 10^{-6}$ (orders 4 and 8); the
  9 s soft-sediment run drains to $\le 3 \times 10^{-6}$ of its peak. With cubic receivers
  (before the interpolation matched the stencil) the first two were 0.11–0.26% and
  0.07–0.18%, and reciprocity $\le 1.5 \times 10^{-6}$. At order 4 with $C = 0.85$ they were
  0.7–1.6%, 0.3–0.8%, 1.9%, $1.5 \times 10^{-4}$ and $1.8 \times 10^{-3}$.
- **Order 8 is the default, with $C = 0.3$** (`config::default_courant`). Order 4 at
  $C = 0.85$ benefits from cancelling time and space errors; orders $\ge 8$ at 0.85 are
  *worse* than order 4, and only beat it at $C \approx 0.3$–$0.5$. Order 8 at 0.3 needs about
  4.5 points per S wavelength for about 1% after ten wavelengths, where order 4 needs 8–9:
  about 3.4× less compute and 5.5× less memory. The price: the stress-imaging free surface
  (second order) is now the largest error in its test (2.4%), and layer reflection rose to
  $-68$ dB. Time–space-optimised coefficients looked good for one speed in a von Neumann
  analysis and failed for other speeds, so they were not adopted. Details: the log.
- **Threads in the browser** need four things, each of which failed once: (1) a std built
  with atomics: nightly `-Z build-std=panic_abort,std`, pinned in `web/build.sh` and both
  Containerfiles; (2) explicit linker flags for an imported, shared memory and the TLS
  exports (`--shared-memory --import-memory --max-memory --export=__wasm_init_tls …`).
  rustc doesn't add them, and without them initThreadPool fails with "Memory could not be
  cloned"; (3) wasm-bindgen-rayon's worker helper imports `'../../..'`, which only bundlers
  resolve, so `build.sh` rewrites it to the module file; (4) a cross-origin-isolated page
  (COOP/COEP in `web/sws.toml`). Without isolation the page falls back to the
  single-threaded build.
- **Browser speeds** (order 4): 51 Mcell/s on one thread with SIMD128 (about native per core),
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
- **Snapshots:** colour $\propto \sqrt{|u|}$, scaled by the 99.5th percentile of $|u|$. A
  Gaussian moment rate leaves a static stress at the source that swamps pressure (`p`)
  sections; velocity has no static part, so the examples snapshot velocity. In a soft
  near-surface layer the surface amplitudes dominate a vertical section; that is physics.
- **The CPML instability along a free surface.** With a free surface, a soft layer
  ($v_p/v_s \approx 3.3$) and CPML side layers, the field grew exponentially from about 2 s,
  starting in the layer corners at the surface (Festa et al. 2005). The free surface alone
  (no PML) was stable; larger $\alpha$, $\kappa_{\max} = 4$ and a smaller $\Delta t$ didn't
  help. The fix is a plain damping $\exp(-p\, d\, \xi^3 \Delta t)$ in the $x$ and $y$ layers,
  only when there is a free surface: `pml_sponge`, default $p = 0.2$. The measured trade-off
  is in the log.
- **Sources on the free surface act on half a cell.** $v_x$, $v_y$, $\sigma_{xx}$,
  $\sigma_{yy}$, $\sigma_{zz}$ and $\sigma_{xy}$ on the surface plane stand for half a cell
  (their other half is the mirror image), so their source weights are doubled
  (`Solver::cell_share`). Without it, reciprocity with one end on the surface was off by
  exactly 2.
- **The free surface mirrors velocities evenly.** First version: second-order derivatives
  next to the surface. They were not the adjoint of the stress-imaging operators, energy
  wasn't conserved, and the free-surface misfit was 7.6%. Mirroring $v_x$, $v_y$ and $v_z$
  evenly ($\sigma$ oddly) makes the pair adjoint for any stencil length and brought it to
  1.9%. The textbook image method (Fichtner 2011 §3.2.3; Graves 1996) instead derives the
  velocities above the surface from the boundary conditions. It hasn't been tried here; any
  alternative must pass the reciprocity and long-run tests.
- **Receiver and source interpolation is as long as the stencil** (`Solver::stagger_weights`):
  $2L$-point Lagrange from the staggered values onto the node. Averaging two values costs
  about $\cos(kh/2)$, roughly 2% at 8 points per wavelength; the cubic
  $(-1, 9, 9, -1)/16$ cut the misfits roughly threefold, and at order 8 the eight-point
  weights halved them again (at 4.5 ppw and $2.5 f_0$ the cubic loses 7.6% per receiver,
  eight points 1.0%). Reciprocity rose from $1.5$ to $1.9 \times 10^{-6}$: more terms of
  single-precision round-off, not a loss of the transpose.
- **Performance** on the container's 4 cores (Apple Silicon, Colima VM with 8 GiB): order 4
  about 220 Mcell/s, 8 about 160, 12 about 130, 16 about 105. Memory is 17 padded f32
  arrays plus CPML memory: the examples need 0.1–0.2 GiB. The work splits by $z$-plane, each
  plane a row at a time through nine scratch rows.
- **Next:** elastic FWI; the README's roadmap lists the steps.
