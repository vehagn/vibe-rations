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
| `just physics` | Whole simulations against analytic answers, misfits printed | ~3 min |
| `just test` | Everything | ~3 min |
| `just accuracy` | The stencil order and time step study (prints a table) | ~2 min |
| `just study [name]` | The inversion study: misfits, stages and preconditioning on small surveys | ~30 min |
| `just drift` | How far the rebuilt forward wavefield drifts (`gradient --reconstruct`) | ~3 min |
| `just gradient-accuracy [name]` | The solver's gradients against the exact kernels of a homogeneous full space | ~40 min |
| `just lint` | rustfmt check and clippy with `-D warnings` | |
| `just ci` | `lint` then `test`: what must pass before every commit | |
| `just fmt` | Format the code | |
| `just check <config>` / `just run <config>` | Set up / run a simulation | |
| `just memory <config>` | The memory a run needs, part by part, without building it (`--gradient`, `--invert`: those runs') | instant |
| `just remote-up <host>` / `just tunnel <host>` | On the host: start `vibe-rations serve` on a remote machine (Docker), and forward port 8089 to it | |
| `just gradient <config> <observed>` | Misfit gradient for vp, vs, ρ against observed traces (`--reconstruct`: rebuild the forward field backwards) | |
| `just invert <config>` | L-BFGS inversion as the configuration's `[inversion]` says | |
| `just marmousi2` | Fetch Marmousi2 (467 MB, once) and make its model files, for `examples/marmousi2*` | ~2 min |
| `just marmousi2-bench` | Time forward and gradient runs on it at orders 4 and 8, with their memory | ~10 min |
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
| `remote-runs` | Runs on another machine: `serve`, `run --remote`, the page's server mode, the solver image |

## Map

```
src/
  main.rs        the CLI: run [--remote | --shots] | check | memory | gradient | invert | extrude | serve | remote
  adjoint.rs     the adjoint run: residuals back through the forward solver, reversed
  gradient.rs    misfit gradients for ρ, λ, μ (and vp, vs): checkpointed replay or the field
                 rebuilt backwards (Strategy), correlations
  inversion.rs   preconditioned L-BFGS over shots: relative unknowns, bounds, line search, stages
  misfit.rs      misfits and their exact derivatives: L2, time shifts, phase; zero-phase low-pass,
                 mutes and windows
  config.rs      the TOML configuration (serde, deny_unknown_fields); default_courant
  grid.rs        the grid and its halo (as wide as the stencil), indexing
  memory.rs      the memory a run needs, from its configuration alone
  model.rs       vp, vs, rho on the nodes: layers, .npy files (on disk or in memory), spheres
  solver.rs      fields, staggered material, the update kernels, free surface, PointOp
  fd.rs          the staggered Stencil (Taylor or optimised, order 2–16), a row at a time
  cpml.rs        CPML profiles and memory slabs, the sponge factors
  source.rs      wavelets, moment tensors, how sources enter the grid
  tdt.rs         the time-dispersion transform: pre-warped sources, un-warped traces
  receiver.rs    recording traces through PointOps
  section.rs     model files from 2D SEG-Y sections: cell averages, smoothing, extrusion along y
  server.rs      `vibe-rations serve`: runs kept for clients that come and go, an HTTP API
  remote.rs      its client: submit, follow, fetch; an SSH tunnel through the system's ssh
  view.rs        sections and gathers of a live run, and an inversion of synthetic data stepped
                 a unit at a time with pictures of it (the browser and the server share it)
  snapshot.rs    sections and volumes during a run
  simulation.rs  a run set up from a configuration; outputs
  image.rs, npy.rs
web/
  src/lib.rs     wasm-bindgen wrappers (Simulator, Inverter) over View and InversionView
  src/settings.rs the TOML read and edited for the page's form (toml_edit), the headless
                 commands, an inversion's start and truth; native tests
  www/           the page: index.html, app.js (wiring), worker.js (the solver here),
                 remote.js (the server's solver, same messages), config.js (settings.rs on
                 the page's thread), form.js, commands.js, fwi.js, plot.js, util.js, style.css
  build.sh       builds web/dist: single-threaded and threaded WebAssembly, page, examples
  sws.toml       static-web-server settings, with the cross-origin isolation headers
  Containerfile  the hosted page's image
Containerfile    the solver's image: `vibe-rations serve` with the page, for remote machines
docs/
  references.bib every work cited anywhere in the repository
  research-log.md experiments, results and decisions, for the paper
  analysis/      von Neumann analysis (pure Python): dispersion, time-space, time dispersion
deploy/kubernetes/  Deployment, Service, HTTPRoute (Kustomize)
.github/workflows/  ci.yaml (`just ci`), web.yaml and solver.yaml (build and push the images)
tests/
  physics.rs     whole simulations: analytic solutions, layers, reciprocity, stability
  accuracy.rs    the stencil-order study (ignored; `just accuracy`)
  inversion_study.rs  misfits, stages and preconditioning on small surveys (ignored; `just study`)
  reconstruction.rs  the rebuilt wavefield's drift (ignored; `just drift`)
  gradient_accuracy.rs  gradients against the exact kernels (studies ignored; `just
                 gradient-accuracy`), and a baseline at the coarsest grid
  common/        closed-form solutions and setup helpers the tests share; stokes.rs: a point
                 force's field with its gradient for any time function, and the exact kernels
  cli.rs         the binary, run as a process
  remote.rs      the server and client as processes: same results as local, reattach, restart
  examples.rs    every example sets up without warnings
examples/        annotated configurations; the marmousi2 ones need `just marmousi2`'s model files
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
  Adjoint FWI depends on it; the reciprocity and dot-product tests check it.
- SI units; $z$ points down; arrays are $x$-fastest, padded by `grid.halo` nodes.
- Comments say why, not what. Doc comments on public items.

## Notes

Decisions and lessons, newest first. Add to them when you learn something the code
doesn't show.

- **The gradient against the exact kernels** (`tests/gradient_accuracy.rs`,
  `tests/common/stokes.rs`; the log's 2026-11-21 entry). In a homogeneous full space the
  kernels of an L2 trace misfit are closed-form correlations of Stokes's solution. Compare the
  solver's gradient times the sample interval over the cell volume, not the bare gradient: the
  misfit sums samples, it doesn't integrate them. And compare it with the exact correlations
  taken where the solver takes them (`Problem::grid_kernels`: ρ at the six velocity
  positions, half each; λ and μ's normal part at the node; μ's shear parts at the twelve
  shear positions, a quarter each), or the grid's sampling of the model, $O(h^2)$ and 4% for
  vs at 3.5 points per wavelength (11% for ρ), swamps the wavefields' error (0.26%). Kernels
  are singular at the source and receiver: compare outside half a dominant wavelength. Error
  in the traces and in the gradient differ: the gradient's was a quarter of the trace's at
  the default's coarsest grid, but with order 4 or the time error kept about the recorded
  residual's, 1.3–1.7 times the trace's.
- **Mutes and windows** (`[[inversion.mute]]`, `[[inversion.stage.mute]]`,
  `misfit::mute_weight`, `misfit::offsets`). A mute is a fixed diagonal weight by time and
  horizontal offset, applied before the filter to observed and synthetic traces alike, so the
  misfit is $\chi(FMs, FMd)$ and the adjoint source $MF\,\partial\chi$ (a unit test checks the
  pair transposes to $10^{-12}$). Before, not after, as the codes a research agent read do it
  (DENISE, SeisFlows, LASIF filter first): filtering first smears a muted arrival back across
  the line by about a period of the corner, and in its prototype turned a kept arrival's
  +0.100 s shift into -0.095 s; muting first leaked nothing. Offsets are between the nodes source and receivers snap to,
  worked out from the configuration (`offsets`) only when a stage mutes, and a misfit keeps
  them and its muted copy of the observed traces (both in `Misfit::held_bytes`). Weights and
  quiet traces are judged on the muted traces, so a window that removes a trace's arrivals
  drops the trace. Two mutes make a window; the stage's add to the inversion's. The
  exponential shape is time damping from a line rather than from first-arrival picks, which
  would depend on picking. Times are
  $m\,\Delta$ from the run's start, not the recorded fields' half-step times: a taper hides
  the difference.
- **The misfits meet the web page** (merged 2026-11-21). `InversionConfig.iterations` is now
  an `Option` (stages give their own), so ask `InversionConfig::stages()` how many iterations
  run, never the field. `view::InversionView` runs the stages: a later stage's start
  (`Inversion::begin_stage`: its filtered misfits, L-BFGS afresh, the receivers' side
  re-estimated) is a unit of its own (`Stage::Staging`), `Stage::Iterating` counts iterations
  over every stage, and a stage that finds no step ends that stage, not the inversion. The
  page's Gradient asks for `iterations = 0`; with stages that would be refused, so
  `InversionConfig::gradient_only` turns it into the first stage with no iterations (the
  gradient as that stage measures the misfit), in the view and in
  `memory::estimate_synthetic_inversion` alike (the estimate, made first, refused it once). The FWI tab chooses the misfit; stages are written in the TOML, and the
  Iterations field steps aside when there are any.
- **The web page edits, inverts and explains** (2026-11-21). The TOML stays the source of
  truth: the Settings form reads it (`settings::settings`, defaults filled in from `Config`) and
  writes into it (`settings::edit`, `toml_edit`), so comments and layout survive; new keys are
  typed by name (`INTEGER_KEYS`) and new tables follow their siblings (a table without a
  position follows the one before it). The page loads the single-threaded `pkg/` on its own
  thread for this, whichever solver runs. The FWI tab inverts synthetic data: the
  configuration is the start, a `[model]` block the truth, and `settings::synthetic` writes both
  with one `dt`, `pml_vp_max` and `optimised_for` (`inversion::shared_operator`: the faster
  model's time step and layers, the slower one's band), so the two are sampled alike and the
  downloaded pair, run through `run --shots` and `invert`, gives the page's misfits (checked:
  51.4, 20.7, 14.4, 8.6, 5.5% for `examples/sphere-inversion.toml`, in the browser, on the
  server and on the command line). `view::InversionView` steps it a unit at a time (a shot of
  observed traces, the start, an iteration) for the worker and for `vibe-rations serve`, where
  an inversion is a run of `kind=inversion` in a thread of its own; its pictures
  (`view::Pictures`) are counted in `memory::estimate_synthetic_inversion`, exactly. Two things
  it turned up: `Simulation::run` read the clock on every step, which panics in WebAssembly, so
  it now does only when reporting; and designing optimised coefficients took about 0.7 s
  (order 12, a band the design must shrink), paid by every run an inversion sets up, so
  `Stencil::optimised` remembers its last eight designs (bit for bit the same; 125 ns). The
  built-in browser's previews come from the main checkout's `.claude/launch.json`, so a
  worktree's build needs an entry there pointing at it (and its `out/`, since a separate
  `CARGO_TARGET_DIR` isn't in the shared volume).
- **The receivers' side of the illumination** (`receiver_runs`, `Inversion::receiver_energy`).
  The runs inject through `Solver::adjoint_op`, exactly as the adjoint injects residuals, and
  each trace by `Misfit::receiver_weights`, so velocity and pressure receivers weigh as they
  do in the gradient: exactly for L2 without mutes (with them, by each trace's mean squared
  mute weight, so a trace muted entirely weighs nothing), by each observed (muted) trace's
  energy (or slope's) for the time shifts and phases, which scale every trace by its own
  amplitude. Weighing those by
  component, as at first, let pressure (10⁶ times velocity) take over the estimate. A run
  needs no sources in the configuration; it steps the solver itself. With a Hadamard matrix's signs
  the cross terms cancel exactly (a unit test); with random ones only on average, which is
  enough for a preconditioner. Off by default: with the top fixed it changed little; with
  nothing fixed it helped (the log). The estimate is made at each stage's start and held
  (8 bytes per node, in the memory plan); the source side is still recomputed every
  iteration.
- **Scaling the parameters** (`[inversion.scale]`). A factor $w$ per parameter on the
  preconditioner is exactly a change of variables ($x = Dz$, $D^2 = w$), so it leaves the
  L-BFGS pairs, the Armijo test and the bounds valid; 1 keeps every result bit for bit. Don't
  equalise the classes automatically: relative unknowns already balance $v_s$ and $\rho$, and
  in a research agent's trial equalising amplified $v_p$, which these data barely constrain,
  into cross-talk ($v_p$ error 109–115% after 5 iterations). A
  smaller $\rho$ factor trades density's recovery for less cross-talk.
- **Misfits and frequency continuation** (`src/misfit.rs`, `[inversion] misfit`,
  `[[inversion.stage]]`). Every misfit is a function of the samples with its *exact*
  derivative, so the inversion's finite-difference check covers it
  (`the_misfits_gradients_match_finite_differences`): don't swap in a textbook adjoint source
  that assumes something about the data. Tromp et al.'s time-shift source $\dot s / \int s
  \ddot s$ assumes the observed trace is a shifted copy, and was 45% off for traces of
  different shape. Weights and trace selections come from the observed traces only; one that
  depends on the synthetic makes the misfit jump between models. The time shift is a
  parabola through the correlation's peak: continuous, but its derivative jumps where the
  peak changes lag, so check it at small steps ($h = 10^{-3}$, where round-off in the
  shifts holds it near $10^{-4}$). The phase misfit is
  regularised ($z/\sqrt{|z|^2 + \epsilon^2}$): taken literally, $\Delta\varphi^2$ has a branch
  cut and blows up where the synthetic lacks an arrival. The forward–backward Butterworth
  filter, started at rest, is symmetric, so it is its own adjoint; zeros appended to the end
  keep it so and stop the forward pass's ringing from being cut off (it distorted the last
  ~3/f_c of every trace), but steady-state starts or mirrored padding would not. A stage rebuilds the misfits (balancing weights from the filtered data)
  and clears L-BFGS, whose pairs describe the previous misfit.
- **Rebuilding the forward wavefield backwards** (`gradient_with`, `Strategy::Reconstruct`,
  `gradient --reconstruct`; the log's 2026-11-20 entry). `Solver::step_back_with` undoes a
  step with the same arithmetic, the kernels taking a sign of ±1 (which keeps the forward step
  bit for bit) and skipping the CPML; `Band` supplies what its stencils read from outside the
  nodes outside the layers: velocities on the whole shell L deep, edges included (the
  gradient's strain rates just outside the region read them), stresses on the faces they are
  differentiated across. Records are of the state *before* each step. The band is saved for
  segments, each replayed from a checkpoint at its start, which resets the drift; the whole
  run as one segment needs no replay but more memory than checkpoints except on large grids.
  The drift is single-precision round-off, $10^{-7}$–$10^{-6}$ of the gradient, set by the
  loudest stretch of the field in a segment; the layers and the free surface add nothing
  (rigid walls drift the same), water five times more (not traced). Inside the layers the
  rebuilt field and gradient mean nothing. **Checkpoints stay the default**, also in the
  inversion: switching is a separate decision. When `gradient` gains an array, add it to
  both `Plan::gradient_memory` and `Plan::reconstruction_memory`; when the step gains an
  operation, `step_back_with` must undo it in reverse order (`steps_back_inside_the_band`
  and `the_rebuilt_wavefield_follows_the_run` catch a miss).
- **Marmousi2** (`just marmousi2`, `examples/marmousi2*`, `vibe-rations extrude`; the log's
  "Marmousi2" entries). The model isn't in the repository: `examples/marmousi2/fetch.sh`
  fetches 467 MB of SEG-Y from the Madagascar project's mirror (SEG's own link answers 403)
  into git-ignored `data/` and checks SHA-256 sums; `convert.sh` writes sections, (nz, nx),
  into git-ignored directories beside the configurations, so their paths have no `..` (the run
  server refuses those). `extrude` averages a 2D section onto the nodes (Backus's harmonic
  moduli, a node fluid by majority, so the seafloor stays within half a node). The soft
  sediments under the seafloor ($v_s$ from 303 m/s at 25 m) set the grid and hold most of the
  late energy and the gradient. The whole section's gradient needs 15 GiB, so gradients use
  the central 6 km (4.3 GiB). **Inverting it** (`tests/marmousi2_study.rs`, `just
  marmousi2-study <name>`; the log's Marmousi2 entries of 2026-11-21): from the start smoothed
  over 250 m, the slow seafloor waves carry most of the misfit and turn $v_s$'s gradient the
  wrong way below 750 m; of the tools on main, only muting them (`[[inversion.mute]]` keeping
  what arrives before a 600 m/s line) stops that, and it also removes the low-frequency slow
  waves through which the stages improved $v_s$ (try a mute per stage next). Judge a variant by the
  error ratio and the cosine with the true update by depth band, not by the misfit: L2 fits
  the most data and recovers the least. Each variant is 45–65 minutes and 4.4 GiB, so run
  them one at a time. Without the files, `tests/examples.rs` only parses the
  configurations and `web/build.sh` leaves them out; the web image's `marmousi2` stage fetches
  and converts them itself, so the hosted page has them.
- **Model files in the browser** (`model::Files`). A browser has no files: the worker asks
  `Simulator.modelFiles` what a configuration names, fetches those relative to `examples/`,
  and hands them over as `ModelFiles`; `Model`, `memory` and `Simulation` read them from
  `Files::Memory` as they would from `Files::Dir`. Sections, (nz, nx), keep the downloads
  small: Marmousi2 at 25 m is 1.2 MB as sections against 47 MB as a 41-node-wide grid.
  Inversions in the browser (`Inverter`, `InversionView`) still take layers and spheres only:
  they read models from a directory, which the page hasn't got.
- **The inversion meets this branch's defaults** (merged 2026-11-20). (1) `Inversion::new`
  pins `optimised_for` at the starting model's points per wavelength (the coefficients
  would otherwise change with every model it tries) as well as $\Delta t$ and
  `pml_vp_max`, and takes the vp bound from the stability limit of the stencil it runs,
  not Taylor's. (2) The inversion's finite-difference check failed at
  $3.1 \times 10^{-5}$: its surface line's end receivers, fine with the cubic, read the side
  layers with eight-point interpolation. Clear of them, the check scatters between
  $10^{-7}$ and $5 \times 10^{-5}$ for $h$ from $2 \times 10^{-3}$ to $1.6 \times 10^{-2}$:
  round-off in single-precision traces below $h \approx 10^{-2}$, the extrapolation's $h^4$
  above. It now runs at $h = 10^{-2}$ ($4.0 \times 10^{-6}$) and refuses receivers next to
  the layers. A check that looks off by a few $10^{-5}$ at one $h$ may be noise: try
  several before blaming the gradient. (3) The gradient's illumination array and the
  inversion's own arrays are in the memory plan, exactly (`memory --invert`).
- **Inversion** (`Inversion`, `invert`): L-BFGS on relative changes $m/m_0 - 1$, each source
  a shot, with $\Delta t$, the layers' damping (`pml_vp_max`) and the optimised
  coefficients' band (`optimised_for`) fixed at the starting model's so the misfit is smooth
  and the gradient exact. Unpreconditioned, the updates pile
  up around the sources and the surface receivers, and a sphere below hardly moves (55% of
  the misfit left after 4 iterations, the sphere's $v_s$ error 99%). Dividing by the
  illumination (source side) and holding the top fixed around the sources and receivers
  (`fixed_above`) gave 6.6% and 78% after 6. The line search almost always takes the full
  step at the first try. Forces recorded at the surface barely see $v_p$ in a sphere half
  a P wavelength across. The bounds keep $\lambda \ge 0$ whichever parameters are inverted
  ($v_p$ projected first; $v_p \ge \sqrt2\, v_s$ where $v_s$ is fixed): without that, an
  inversion of $v_p$ alone reached an unphysical model and aborted. A trial the solver can't
  run shortens the step; unknowns pinned at a bound drop out of the direction. Check an
  inversion's gradient *down the gradient*: along a random direction the predicted change
  cancels and round-off swamps it. Details: the log.
- **The adjoint and gradients meet this branch's defaults** (merged 2026-11-20). (1) The
  time-dispersion transform: `gradient` un-warps the traces after the forward run, as `run`
  does, so the misfit is of the same traces as observed data from `run`, and
  `Simulation::recorded_residuals` takes its derivative back through the transform's
  transpose to the samples the solver recorded, which `Adjoint` injects. (2) Optimised
  coefficients follow the model's slowest wave, so a finite-difference check or an inversion
  must pin `[scheme] optimised_for`, or each model gets its own operator (`with_model`
  warns). (3) Receivers interpolate from as many values as the stencil spans, four nodes at
  order 8, so "next to a layer" now means within four nodes: two test receivers that sat
  two nodes in (fine for the cubic) read damped values and made the dot-product test
  2.4e-7 and the gradient check worse, looking like a layer leak until
  `receivers_next_to_layers` named them. The dot-product and gradient checks now assert none.
- **Remote runs** (`vibe-rations serve`, `run --remote`, `remote`; README, "Remote runs").
  The server keeps runs by id in its data directory and steps them in one scheduler thread,
  a chunk of about 0.1 s at a time; request threads draw pictures between chunks under the
  run's lock. It listens on localhost (0.0.0.0 inside its container, published on the host's
  127.0.0.1): SSH is the security, so never publish it more widely without adding
  authentication. The page and the CLI share the API; `web/www/remote.js` answers the
  worker's messages so `app.js` drives either solver the same way. A run the page set up but
  never ran is deleted when replaced. After a server restart, finished runs keep their
  results and set-up ones are interrupted (no checkpoints yet). The dev container has no SSH
  agent, so its CLI reaches remotes through a tunnel opened on the host
  (`http://host.docker.internal:8089`); `--remote user@host` itself (the CLI's own tunnel) is
  untested against a real remote. `tests/remote.rs` checks remote traces match local ones bit
  for bit. In the built-in browser, open local servers through `.claude/launch.json` entries
  (`solver`), not `navigate`: on other origins it blocks nested workers, which threaded
  WebAssembly needs (the `remote-runs` skill). A run is headless (`Status::headless`) when
  nobody has asked for its pictures for 5 s; it then steps in ~2 s chunks instead of 0.1 s
  and draws nothing: 18.7–19.6 s against 21.3–22.1 s watched (4.9 M nodes, 4 cores), about
  12% less. **A run's locks:** std's Mutex isn't fair, and with one lock per run the
  scheduler took it straight back after each chunk, so a status request could wait out a
  whole headless run (and a first measurement, which found no difference, was an artifact
  of picture requests starving the same way). Status lives in `Run::shared`, never held
  while stepping; requests for the wavefield go through `Run::state`, which registers them
  so the scheduler lets them in before its next chunk.
- **Memory is counted before it is allocated** (`memory::estimate`, `vibe-rations memory`).
  It shares `simulation::timing` with `Simulation::new`, takes the model's extremes from
  `Model::extremes` (layers, spheres, or `.npy` files streamed), and must stay exact:
  `counts_what_a_run_allocates` compares it with the arrays a real run holds. When a run
  gains a large array, add it to `memory::Plan` and to `Solver::memory_bytes`. `run` and
  `check` refuse runs larger than `MemAvailable`; the web app refuses runs over 4 GiB. A
  volume snapshot briefly needs 32 bytes per node (a node list, values and shading), which
  could be cut. A gradient run (`Plan::gradient_memory`, `memory --gradient`, checked by
  `gradient` before it allocates) holds the adjoint's solver, residuals, checkpoints and a
  replayed segment besides: about $42\sqrt N$ bytes per padded node at the default interval,
  25 times the forward run for `examples/layered.toml`. `counts_what_a_gradient_run_allocates`
  compares it with `gradient::Run::memory_bytes`, which sums the arrays a real run holds at
  its peak; when `gradient` gains an array, add it to both. `gradient` frees the replay before
  making the gradient, so the output never adds to the peak. An inversion
  (`Plan::inversion_memory`, `memory --invert`, checked by `invert`) adds every shot's
  observed traces, three models, the sums over the shots, and per unknown $8(2m + 5)$ bytes
  for L-BFGS memory $m$; `counts_what_an_inversion_allocates` compares it with
  `Inversion::memory_bytes`, measured the same way. Keep `Inversion::new` from holding a
  solver while it evaluates (it drops its set-up run).
- **Stencil coefficients are optimised by default** (`[scheme] coefficients`, opt out with
  `"taylor"`; `Stencil::optimised`). With the time error removed, the phase error per node,
  $\lvert \mathbf K \rvert - \lvert \mathbf k \rvert h$, is the whole error, so the
  coefficients minimise its largest value over all 3D directions and all wavenumbers up to the
  run's shortest wavelength (Lawson's algorithm; Holberg 1987 did the same for group
  velocity), weighted by $(k/k_{\text{band}})^{-2}$ and capped at $10^{-4}$ per node
  (`PHASE_TOLERANCE`), shrinking the band if need be. Both choices were measured: an
  unweighted minimax, or a cap of $10^{-3}$, spread the error onto long waves, and P waves
  and the explosion test came out 2–4× worse than with Taylor; with them nothing measured
  got worse beyond noise. Gains: S at 3.0 ppw, 19% → 7.3% (order 8), 6.5% → 0.89% (12),
  2.9% → 0.42% (16); order 8 at 4.6 ppw, 1.2% → 0.19%. The band follows the model's slowest
  wave, so an inversion must pin it (`optimised_for`) or the operator changes between
  iterations. The interpolation weights are optimised over the same band.
- **Leapfrog's time error is removed from the traces by default** (`[scheme]
  remove_time_dispersion`, opt out with `false`; `src/tdt.rs`; Wang & Xu 2015, Koene et al.
  2018). It is a warp of the frequency axis, the same for every medium and stencil, so
  pre-warping the sources and un-warping the traces leaves only the spatial error, and
  $C = 0.85$ (the default with it) costs no accuracy. Order 8 then matches $C = 0.3$ in 2.8
  times fewer steps; orders 12 and 16, time-limited at 0.3, become 8–12 times more accurate.
  Our pair is Koene et al.'s eqs. 13–17 (checked against the article: no amplitude factors,
  zero-padding to twice the length); they don't discuss the time origin. Timing matters here:
  forces enter at $n\Delta t$, moments at $(n + \frac12)\Delta t$, velocities are read at
  $(n + \frac12)\Delta t$ and stresses at $(n + 1)\Delta t$; the warp is about $t = 0$.
  Traces are final only after the last step (`Simulation::finish`, which `run` calls);
  snapshots stay warped; the layers are only approximately transformed. **For the adjoint:**
  residuals of finished traces must go through the transposed transform
  (`Simulation::trace_warp(c)` → `Warp::apply_transpose`) before they are injected. Koene et
  al. criticise Wang & Xu for using that adjoint as the forward transform; here it is used
  only for adjoint sources, where it is the exact transpose of the data operator.
- **The free-surface test's 2.5% is mostly its reference.** Doubling the P wave is exact
  only for a plane wave; on finer grids (`free_surface_on_finer_grids`) the misfit tends to
  about 2.2%, so stress imaging's own error is about 0.36% at 8 points per P wavelength
  (0.23% with the time error kept, which partly cancels it). A higher-order free surface is
  therefore less urgent than it looked; measure it on Lamb's problem with an exact solution.
- **A research brief** (another agent's literature review and von Neumann analysis,
  2026-11-20) was checked: its numbers reproduce (`docs/analysis/time_dispersion.py`), the
  transform and the stencil-matched interpolation came from it, and the log's entry lists
  what remains, in its open questions. Its references came from memory: check any on
  Crossref before citing it.
- **Gradients are those of the discrete run** (`gradient::gradient`): the multipliers of the
  forward equations are $W B^{-1} \tilde v$ and $-W \tilde\varepsilon$, with the adjoint
  strain $\tilde\varepsilon$ *summed* from the adjoint strain rates (not $M^{-1}\tilde\sigma$:
  singular in fluids, cancels digits where $\mu \ll \lambda$), then the chain rule through
  the material averaging (once per slot a node fills: the grid's last planes clamp) and the
  free surface's plane-stress moduli. Finite differences agree to $\le 2.8 \times 10^{-4}$
  (CPML at orders 4 and 8; rigid walls with water and uneven spacing) only after
  Richardson's extrapolation: plain central differences in single precision have an $h^2$
  error of percents down to $h \approx 2 \times 10^{-3}$, below which round-off takes over.
  In a finite-difference check, keep $v_p$ below the model's fastest: it sets $\Delta t$ (via
  the Courant fraction) and the layers' damping. The forward field is replayed from
  checkpoints (`Solver::state`, `restore`), bit for bit; a unit test checks every interval
  gives the same gradient. Receiver components may not repeat: a second `p` once went into
  the adjoint run but not the adjoint strain and flipped the gradient's sign.
- **The adjoint is the forward solver run backwards** (`Adjoint`), not separate kernels.
  Transposing the discrete step gives the same step in reversed order with $\tilde v =
  \lambda_v$, $\tilde\sigma = -\lambda_\sigma$. Adjoint step $m$ pairs with forward step
  $n = N - 1 - m$: a velocity residual from forward step $n$ enters during adjoint step $m$
  as a unit force (`Solver::adjoint_op`); a pressure residual, read from $\sigma^{n+1}$, one
  adjoint step earlier as an isotropic stress source $\Delta t\,K/V$, and the last one is the
  initial state, mirrored above a free surface (a receiver a node below the surface is what
  tests that). The dot-product test holds to $8.3 \times 10^{-8}$ (single precision) with
  rigid walls or CPML, orders 4 and 8. Receivers *inside* a layer break it by up to
  $1.4 \times 10^{-3}$, and a velocity receiver whose interpolation reaches into one (the
  first node outside at order 4, the first four at order 8) by up to $2.8 \times 10^{-5}$,
  since it reads damped values (`Simulation::receivers_next_to_layers`; `gradient` warns). For a dot-product test, give
  the residuals the data's signs and weight the components equally: random signs cancelled
  50-fold and made round-off look like $10^{-4}$, and pressure outweighs velocity by $10^6$.
- **Error propagation into adjoint fields** is a main research interest for what comes next:
  forward-field error (dispersion, interfaces, free surface, layers), reconstruction error,
  and their effect on gradients. Design the adjoint so each can be measured separately (see
  the log's open questions).
- **Physics test baselines at the default, order 8 with optimised coefficients at
  $C = 0.85$, the time error removed** (fail if they get worse; update this list if they get
  better): explosion vs analytic 0.005–0.007%; point force vs Stokes 0.001–0.005%;
  free-surface P doubling 2.52% (mostly the reference, above); CPML reflection
  $1.7 \times 10^{-4}$ of the peak without a free surface, $1.9 \times 10^{-3}$ along one;
  reciprocity $\le 8 \times 10^{-7}$ (orders 4 and 8); the adjoint dot-product test
  $\le 8.3 \times 10^{-8}$; gradients against extrapolated finite differences
  $\le 2.8 \times 10^{-4}$ (order 4's $v_p$, whose directional derivative is small; order 8
  $\le 1.6 \times 10^{-5}$); the sphere inversion leaves 17.9% of its misfit and 91.9% of the $v_s$
  anomaly after 3 iterations, and its gradient matches finite differences to
  $4.0 \times 10^{-6}$ (at $h = 10^{-2}$); the rebuilt wavefield retraces the run to
  $\le 1.2 \times 10^{-6}$ of each array's peak and its gradient the replay's to
  $\le 7.1 \times 10^{-7}$, and on the reconstruction path the finite-difference checks give
  the replay's numbers (order 8 $\le 1.6 \times 10^{-5}$); the gradient against the exact
  kernels at 3.5 points per S wavelength (`tests/gradient_accuracy.rs`), the wavefields' part
  0.26% (S) and 0.12% (P) in the Fresnel zone, +0.19% and −0.10% along a smooth change; the 9 s soft-sediment run drains to
  $\le 3 \times 10^{-6}$ of its peak. With Taylor coefficients: 0.002–0.006% (near exact
  for these long waves), 0.001–0.008%, 2.53%, the same layers, $\le 7 \times 10^{-7}$. With
  the time error kept at $C = 0.3$ (the default before): 0.07–0.19%, 0.03–0.09%, 2.40%,
  $3.9 \times 10^{-4}$, $2.1 \times 10^{-3}$ and $1.9 \times 10^{-6}$; with cubic receivers
  before that, the first two were 0.11–0.26% and 0.07–0.18%. At order 4 with $C = 0.85$
  and the time error kept: 0.7–1.6%, 0.3–0.8%, 1.9%, $1.5 \times 10^{-4}$ and
  $1.8 \times 10^{-3}$.
- **Order 8 is the default** (`config::default_courant`). With the time error kept,
  order 4 at $C = 0.85$ benefits from cancelling time and space errors; orders $\ge 8$ at
  0.85 are *worse* than order 4, and only beat it at $C \approx 0.3$–$0.5$, so 0.3 was the
  default until the time error was removed. Order 8 needs about 4.5 points per S wavelength
  for about 1% after ten wavelengths, where order 4 needs 8–9: about 3.4× less compute and
  5.5× less memory. Time–space-optimised coefficients looked good for one speed in a von
  Neumann analysis and failed for other speeds, so they were not adopted. Details: the log.
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
