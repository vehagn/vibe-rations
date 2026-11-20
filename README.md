# vibe-rations

A 3D elastic seismic wave simulator in Rust, built as the forward engine for elastic full
waveform inversion (FWI).

It solves the isotropic elastic wave equation in velocity–stress form on a staggered grid
(Virieux 1986; Levander 1988; Graves 1996), eighth order in space by default (any even order
from 2 to 16) and second order in time, with convolutional PML absorbing layers (Komatitsch &
Martin 2007) and a free surface. The update loops run in parallel across $z$-planes with
rayon.

Citations use the keys of [`docs/references.bib`](docs/references.bib); the
[research log](docs/research-log.md) records the experiments behind the design, for the
paper this simulator is meant to support.

```
step   1783/1783  t 2.5000 s  peak |v| 2.447e-3 m/s  200 Mcells/s  0s left
```

## Getting started

Development happens in the dev container (`.devcontainer/`): Rust 1.99.1, `just`, and a
cargo cache and build directory in volumes. Open the repository in VS Code or RustRover and
reopen it in the container, or open a shell in it from the host:

```bash
just dev
```

Then, in the container:

| Command | Does |
|---|---|
| `just unit` | Unit tests, the CLI tests and the examples check (seconds) |
| `just physics` | Whole simulations against analytic answers, printing misfits (~20 s) |
| `just test` | All tests, optimised |
| `just lint` | rustfmt and clippy with `-D warnings` |
| `just ci` | Lint, then all tests: what must pass before a commit |
| `just check <config>` | Sets a run up without running it: grid, time step, accuracy, memory |
| `just memory <config>` | The memory a run needs, part by part, without setting it up (see "Memory"); `--gradient` or `--invert` for those runs |
| `just remote-up <host>`, `just tunnel <host>` | On the host: start the solver server on a remote machine, and reach it (see "Remote runs") |
| `just run <config>` | Runs a simulation; `--shots` fires the sources one at a time, into `shot-NNN/` |
| `just gradient <config> <observed>` | The misfit gradient for $v_p$, $v_s$ and $\rho$ against traces in `<observed>` |
| `just invert <config>` | Inverts observed traces for the model with L-BFGS, as `[inversion]` says |
| `just web` | Builds the web app into `web/dist` |

```bash
just run examples/explosion.toml
```

From the host, `just in <recipe>` runs any of these in the container without opening a
shell (`just in ci`). `just dev down` stops the container and keeps its volumes.

Working on the code with an agent: [AGENTS.md](AGENTS.md) has the conventions and the
notes behind the numerics, and `.claude/skills/` has task guides.

## Web app

The same solver, compiled to WebAssembly, runs in the browser: pick an example or edit its
TOML, run and pause, and watch a section of any field on any plane while the gather fills in.
The visitor's computer does all the work, in a Web Worker; the server only hands out files,
so it costs next to nothing to host.

```bash
just in web     # build web/dist in the dev container
just serve      # on the host: http://localhost:8088, with the production server and headers
```

When the page is cross-origin isolated (the server sends `Cross-Origin-Opener-Policy:
same-origin` and `Cross-Origin-Embedder-Policy: require-corp`), the solver runs rayon on one
Web Worker per core; otherwise it runs on one thread. `?threads=N` picks the number, and
`?example=layered` an example. Measured on a 10-core Apple Silicon Mac with the explosion
example: one thread 51 Mcell/s, ten threads 234 Mcell/s.

**Hosting.** `web/Containerfile` builds a 4 MB image: static-web-server, non-root, read-only
root filesystem, with the headers above (`web/sws.toml`). GitHub Actions publish it to
`ghcr.io/vehagn/vibe-rations-web` for amd64 and arm64 on every push to `main`
(`.github/workflows/web.yaml`). `deploy/kubernetes/` has a Kustomization with a Deployment,
Service and Gateway API HTTPRoute; set the hostname, then

```bash
kubectl apply -k deploy/kubernetes
```

The threaded build needs a std built with atomics, so it uses a pinned nightly toolchain
(`-Z build-std`) next to the stable one; the single-threaded build stays on stable.

## Configuration

A run is one TOML file; `examples/` has annotated ones. Unknown keys are errors, so typos
don't pass silently. Positions are in metres from the first node, $z$ points down, and with
a free surface $z = 0$ is the surface.

| Section | What it sets |
|---|---|
| `[grid]` | `points = [nx, ny, nz]`, `spacing` (one value or `[dx, dy, dz]`) |
| `[time]` | `duration`; `dt`, or `courant` (fraction of the stability limit; by default 0.85, or 0.3 at orders above 4 without `remove_time_dispersion`) |
| `[scheme]` | `order` of the spatial derivatives: 2, 4, ..., 16 (default 8); `coefficients`, `optimised` (default) or `taylor`, and `optimised_for` (points per wavelength to optimise for; by default the run's own); `remove_time_dispersion` (default true): take leapfrog's time error out of the traces. See "Accuracy" |
| `[boundary]` | `pml_width` (default 20 nodes), `free_surface` (default true), PML tuning |
| `[[model.layer]]` | Horizontal layers: `top` depth, `vp`, `vs`, `rho`. `vs = 0` is a fluid |
| `[model.files]` | Or `.npy` arrays of shape (nz, ny, nx) for `vp`, `vs`, `rho` |
| `[[model.sphere]]` | Spheres painted over the model |
| `[[source]]` | `position`, `wavelet` (Ricker or Gaussian), and one of `explosion`, `moment`, `double_couple` (strike, dip, rake, m0) or `force` |
| `[receivers]` | `components` (`vx`, `vy`, `vz`, `p`), `[[receivers.line]]`s and `[[receivers.point]]`s |
| `[[snapshot]]` | `field`, `plane` (`xy`, `xz`, `yz` or `volume`), `at`, `interval`, `formats` (`png`, `npy`) |
| `[output]` | `dir`, `record_every` |
| `[inversion]` | For `invert`: `observed` traces, `parameters` (`vp`, `vs`, `rho`), `iterations`, `memory`, `initial_step`, `max_step`, `balance`, `fixed_above`, `precondition`, `water_level`, bounds `vp_range`, `vs_range`, `rho_range` |

The wavelet is the moment rate of a moment-tensor source and the force of a force source. A
Ricker wavelet peaks at 1, so its amplitude is a peak rate; a Gaussian has unit area, so its
amplitude is the final moment. A source with a non-zero final moment leaves a static stress
around it, which dominates pressure snapshots near the source. Velocity has no static part.

## Output

In the output directory:

- `vx.npy`, `vy.npy`, `vz.npy`, `p.npy`: traces, shape (receivers, samples), float32
- `time.npy`: the sample times. Velocities lag stresses by half a step, so the first sample
  is at $\Delta t/2$. By default the traces are free of leapfrog's time error (see
  "Numerics"); snapshots are not.
- `receivers.csv`: the receivers' positions after snapping to nodes
- `*_gather.png`: one column per receiver, each trace scaled to its own peak
- `run.json`: the run's settings, time step, accuracy figures and warnings
- `snapshots/`: sections as PNG (blue negative, red positive, colour $\propto \sqrt{|u|}$,
  over the model in grey with the absorbing layers darker) or `.npy`

```python
import numpy as np
vz, t = np.load("out/explosion/vz.npy"), np.load("out/explosion/time.npy")
```

`gradient` compares the configured model's traces with observed ones (`<component>.npy` from
an earlier `run`, sampled the same way: give both configurations the same `[time] dt`) by
the L2 misfit, `--balance` weighting each component by its observed energy, and adds:

- `grad_vp.npy`, `grad_vs.npy`, `grad_rho.npy`: $\partial\chi/\partial v_p$, $\partial\chi /
  \partial v_s$ and $\partial\chi/\partial\rho$ at every node, shape (nz, ny, nx) like model
  files; approximate inside the absorbing layers
- `gradient.json`: the misfit, the component weights, the checkpoint interval, the most
  memory the run held (with the observed traces: what `memory --gradient` works out
  beforehand), and how many receivers sit next to the layers (which makes the gradient
  approximate: move them further in, as many nodes as the stencil reaches, four at order 8)

The misfit is of the traces as `run` writes them, free of leapfrog's time error by default;
the gradient takes it back through the transform's transpose. Give the observed and the
guessed configurations the same `[scheme]`, with `optimised_for` pinned: optimised
coefficients otherwise follow each model's slowest wave, and the two runs would use
different operators.

```bash
just run truth.toml --out out/observed
just gradient guess.toml out/observed --balance --out out/gradient
```

`invert` runs the configuration's `[inversion]`: each source is a shot, whose observed
traces are in `<observed>/shot-NNN/` as `run --shots` writes them (or `<observed>/` itself for
one source). After every iteration it writes the model, `vp.npy`, `vs.npy` and `rho.npy`
(shape (nz, ny, nx), ready for `[model.files]`), and `inversion.json`: the settings and each
iteration's misfit, step and largest change, and the most memory it has held. The time step,
the layers' damping and the optimised coefficients' band (`optimised_for`) are held at the
starting model's, so give the observed run the same `[time] dt` (and `pml_vp_max`, and
`optimised_for` too if the observed traces should come from the same operator); `invert`
checks the shapes and, from `time.npy`, the sample times. `observed` is relative to
the configuration's directory, `--out` (and `[output] dir`) to the working directory.

```bash
just run truth.toml --shots --out out/observed
just invert guess.toml --out out/inverted   # guess.toml: [inversion] observed = "out/observed"
```


## Memory

`just memory <config>` (`vibe-rations memory`) works out what a run needs without building
it, so a model of any size gets an answer at once:

```
      14 GiB  wavefield: 9 arrays on the padded grid
      12 GiB  material: 8 arrays on the padded grid
    4.49 GiB  model: vp, vs and rho
     955 MiB  absorbing layers' memory
      13 KiB  sources' pre-warped time series
    1.79 MiB  traces
      32 GiB  held throughout the run
      32 GiB  at most, plus a few MiB
    7.13 GiB  available now
```

That is `examples/layered.toml` on a $1001 \times 1001 \times 401$ grid. Per node of the grid
padded by the stencil's half-length, the wavefield and material take 68 bytes; the model 12
bytes per interior node; the absorbing layers 24 bytes per layer node; the traces 4 bytes per
sample per receiver and component; with the time error removed, each source 8 bytes per time
step; and a volume snapshot needs 32 bytes per node more while it is written. The time step,
which sets the number of steps and samples, comes from the model's extremes, read from its
layers and spheres, or from its `.npy` files a piece at a time. A test checks that the count
matches what a run allocates, byte for byte.

A gradient run holds much more, and `just memory <config> --gradient` counts it: besides the
forward run, a second solver for the adjoint, the observed traces and the residuals,
$\lceil N/k \rceil$ checkpoints of the wavefield and the layers' memory (36 bytes per padded
node each), $k + 1$ copies of the velocities for the segment being replayed (12 bytes per node
each), the correlation sums and the two strains (112 bytes per node), and the illumination
(8 bytes per interior node). The default interval $k = \lceil\sqrt{3N}\rceil$ about
minimises the checkpoints and the segment together, at about $24\sqrt{3N} \approx 42\sqrt N$
bytes per node, so for
`examples/layered.toml` (1,611 steps) a gradient takes 3.64 GiB against the forward run's
147 MiB, and on the grid above 798 GiB. Rebuilding the forward field backwards instead of
storing it (the roadmap's step 3) is what would cut that. A second test checks this count
against a real gradient run.

`just memory <config> --invert` counts an inversion: one shot's gradient run at a time, every
shot's observed traces, three models (as it stands, the start, and the one stepped to), the
gradient and illumination summed over the shots, and per unknown the L-BFGS pairs ($16m$
bytes for memory $m$) and five more numbers. Unknowns are counted at every free node, so a
model with fluids (no $v_s$ unknown there) needs a little less. A third test checks it against
a real inversion.

`run`, `check`, `gradient` and `invert` refuse a run that needs more than the machine has available
(Linux's `MemAvailable`, or a container's memory limit) before allocating anything. The web
app shows the estimate under the configuration as it is edited, and refuses runs above the
4 GiB a browser gives WebAssembly; the server below has no such limit.

## Remote runs

A machine with more cores and memory can do the computing for both the command line and the
web page. It runs `vibe-rations serve`, which keeps runs going whether or not anyone is
connected, so a long run survives a dropped connection and can be picked up again. The server
listens on that machine's localhost only, and SSH, with your keys and `~/.ssh/config`, does
the authentication and encryption.

On the remote machine, which needs Docker or Podman (`ENGINE=podman`), from your machine:

```bash
just remote-up me@bigbox
```

That starts the image `ghcr.io/vehagn/vibe-rations` there, with the runs and their results
in a volume. `just remote-update me@bigbox` replaces it with the latest image (runs in
progress are lost; results are kept). Without `just`, it is one `docker run`:

```bash
ssh me@bigbox docker run -d --name vibe-rations --restart unless-stopped -p 127.0.0.1:8088:8088 -v vibe-rations:/data ghcr.io/vehagn/vibe-rations:main
```

**Command line.** `--remote` sends the configuration (and any `.npy` model files) over an SSH
tunnel the CLI opens itself, follows the run, and writes its results locally as a local run
would:

```bash
vibe-rations run examples/layered.toml --remote me@bigbox
vibe-rations remote me@bigbox list
vibe-rations remote me@bigbox attach layered-1795184607-3f2a
vibe-rations remote me@bigbox cancel layered-1795184607-3f2a
```

If the CLI stops, the run carries on; `attach` follows it again and fetches its results.
Progress is one line, rewritten in place (a line a minute when not on a terminal): share
done, step, speed and time left. To fire a run off and come back later:

```bash
vibe-rations run big.toml --remote me@bigbox --detach      # prints the run's id
vibe-rations remote me@bigbox list                          # progress and time left
vibe-rations remote me@bigbox fetch big-1795190000-1a2b     # its results, once done
```

The destination can be anything `ssh` accepts, or a URL (`--remote http://localhost:8089`)
for a tunnel opened another way. In the dev container, which has no SSH keys, open the tunnel on the
host (`just tunnel me@bigbox`) and use `--remote http://host.docker.internal:8089`.

**Web page.** The server also serves the page:

```bash
just tunnel me@bigbox     # ssh -N -L 8089:127.0.0.1:8088 me@bigbox
```

and open http://localhost:8089. **Compute** chooses between the server and this browser;
**Runs** lists the server's runs to attach to, including ones started from the command line.
A run started from the page keeps going when the page is closed. The memory line shows what
the server has available. Turning **Live pictures** off leaves progress only: the page asks
for the run's status every 2 s instead of a section and a gather five times a second, which
matters over a slow link; in the browser, the solver stops drawing.

**Headless runs.** The server counts a run as watched while someone asks for its pictures,
and as headless once nobody has for 5 s: runs from the command line, and from the page with
live pictures off or closed. A headless run draws no pictures and steps in chunks of about
2 s instead of 0.1 s (while no watched run shares the machine). On a 4-core machine a
4.9-million-node run took 18.7–19.6 s headless and 21.3–22.1 s watched with pictures five
times a second: about 12% less time, and next to nothing over the network.

**The API** is plain HTTP and JSON (listed in `src/server.rs`): create a run from a
configuration, upload its model files, start, pause, resume or delete it, read its status,
draw sections and gathers as PNG, and fetch its traces and outputs. Each run keeps its
configuration, model files, status and outputs in a directory of its own under the server's
`--data`; after a restart, finished runs' results are still there, and runs that were in
progress are marked interrupted.

## Numerics

**Equations.** With particle velocity $\mathbf{v}$, stress $\boldsymbol{\sigma}$, density
$\rho$ and Lamé parameters $\lambda$, $\mu$, the solver integrates

$$
\rho\,\partial_t v_i = \partial_j \sigma_{ij} + f_i, \qquad
\partial_t \sigma_{ij} = \lambda\,\delta_{ij}\,\partial_k v_k + \mu\,(\partial_i v_j + \partial_j v_i) - \dot{m}_{ij},
$$

with forces $f_i$ and moment-tensor sources $m_{ij}$ (Virieux 1986).

**Grid.** The interior nodes are padded by a halo as wide as the stencil's half-length.
Normal stresses, $\lambda$ and $\mu$ sit on the nodes; $v_x$, $v_y$ and $v_z$ half a node
along $x$, $y$ and $z$; shear stresses half a node along both their axes. Buoyancy is
averaged arithmetically onto the velocity positions and $\mu$ harmonically onto the shear
positions, so shear stress vanishes next to a fluid.

**Derivatives.** A stencil of half-length $L$ (order $2L$) approximates the derivative at a
staggered point $y$, half a node from the nodes around it, as

$$
\partial_x f(y) \approx \frac{1}{h} \sum_{n=1}^{L} c_n \Big[ f\big(y + (n - \tfrac12) h\big) - f\big(y - (n - \tfrac12) h\big) \Big],
$$

with Taylor coefficients in closed form,
$c_n = \frac{1}{2n-1} \prod_{m \ne n} \frac{(2m-1)^2}{(2m-1)^2 - (2n-1)^2}$; order 4 gives
Levander's $c_1 = \tfrac98$, $c_2 = -\tfrac1{24}$. Taylor coefficients are exact for the
longest waves and least accurate for the shortest, which set the grid. By default the
coefficients are instead *optimised* for the run, in the manner of Holberg (1987): for a
plane wave the operator gives $K(\xi) = 2 \sum_n c_n \sin\big((2n - 1)\xi/2\big)$ instead of
$\xi = k h$, and with leapfrog's time error removed (below) the phase error per node,
$\lvert \mathbf K \rvert - \lvert \mathbf k \rvert h$, is the whole error. The coefficients
minimise its largest value over every direction in 3D and every $\lvert \mathbf k \rvert h$ up
to the band $\beta = 2\pi / \text{ppw}$ of the run's shortest wavelength, weighted by
$(\lvert\mathbf k\rvert h/\beta)^{-2}$ so that longer waves keep Taylor-like accuracy, and
capped at $10^{-4}$ per node: if the band can't be had within that, it shrinks. They stay
exact for linear functions. Lawson's algorithm (Rice & Usow 1968) solves the minimax
problem; Zhang & Yao (2013) and Liu (2014) optimise staggered operators in related ways.
Receivers and sources interpolate with weights optimised the same way. The band depends on
the model's slowest wave, so an inversion should fix it with `optimised_for`.

**Stability.** The time step must satisfy

$$
\Delta t \le \frac{1}{v_{p,\max} \, \sum_n |c_n| \, \sqrt{\Delta x^{-2} + \Delta y^{-2} + \Delta z^{-2}}} .
$$

`just check` warns when the grid has fewer points per shortest wavelength than the order
needs (4 at order 8, 5 at order 4, 10 at order 2), taking the slowest wave (S, or P in
fluids) at $2.5 f_0$.

**Free surface.** Stress imaging on the $\sigma_{zz}$ plane (Levander 1988; Graves 1996):
$\sigma_{zz} = 0$ there, $\sigma_{xx}$ and $\sigma_{yy}$ are updated with
$\partial_z v_z = -\frac{\lambda}{\lambda + 2\mu} (\partial_x v_x + \partial_y v_y)$, and
$\sigma_{zz}$, $\sigma_{xz}$, $\sigma_{yz}$ mirror oddly across the surface. Velocities
mirror evenly. With that choice the discrete $z$-derivatives of stress and of velocity stay
negative adjoints of each other next to the surface, for any stencil length, so the scheme
conserves energy there.

**Time dispersion.** Leapfrog's whole time error is a warp of the frequency axis: the
stepped wavefield at frequency $\omega$ is the time-continuous one at
$\omega' = \frac{2}{\Delta t} \sin\frac{\omega \Delta t}{2}$, for every medium and stencil,
because nothing else in a step depends on time. So by default (`remove_time_dispersion`),
each source's samples are pre-warped before the run and each trace is un-warped after it
($\omega = \frac{2}{\Delta t} \arcsin\frac{\omega' \Delta t}{2}$), which leaves only the spatial
error in the traces (Wang & Xu 2015; Koene et al. 2018). On the staggered grid each signal
is warped about its own time stamps: forces enter at $n \Delta t$, moment tensors at
$(n + \frac12)\Delta t$, velocities are read at $(n + \frac12)\Delta t$ and stresses at
$(n + 1)\Delta t$. It is exact outside the absorbing layers; snapshots of the wavefield stay
warped. `remove_time_dispersion = false` opts out.

**Absorbing layers.** CPML (Komatitsch & Martin 2007) with $\kappa = 1$,
$\alpha = \pi f_0$, a quadratic damping profile $d(x)$ and the reflection coefficient of
Collino & Tsogka (2001). Along a free surface, CPMLs can grow without bound when surface
waves in a high-Poisson's-ratio medium run into them (Festa et al. 2005); the side layers then
also damp every field by $\exp(-p\, d\, \xi^3 \Delta t)$ per step, with $\xi \in [0, 1]$ the
depth into the layer and $p = 0.2$ (`pml_sponge`), in the spirit of the multiaxial PML of
Meza-Fajardo & Papageorgiou (2008).

**Sources and receivers.** Receivers interpolate staggered velocities onto nodes from as
many values as the stencil spans, fewer near the edges of the grid, with weights optimised
like the coefficients; with Taylor coefficients they are Lagrange's (the cubic
$(-1, 9, 9, -1)/16$ at order 4).
Sources spread onto the staggered positions with the same weights, divided by the cell's
share of mass (half on the free surface). A source at a point
is therefore exactly the transpose of a receiver there, which is what back-propagating
residuals in FWI needs; reciprocity, $G_{ij}(\mathbf{x}_B, \mathbf{x}_A) = G_{ji}(\mathbf{x}_A, \mathbf{x}_B)$,
holds to $10^{-6}$.

## Accuracy

The scheme is a staggered grid, which is the right family for elastic waves: unlike
conventional (collocated) grids it stays accurate when $v_p/v_s$ is large, as in sediments
(Moczo et al. 2011). Within it, accuracy is mostly a matter of numerical dispersion, and the
choice that matters is the **stencil order together with the time step**.

The leapfrog time stepping is second order. For a plane wave the scheme's dispersion relation is

$$
\sin\!\Big(\frac{\omega \Delta t}{2}\Big) = \frac{v \Delta t}{h} \sqrt{\sum_{i} \Big( \sum_{n} c_n \sin\big((2n - 1) \tfrac{k_i h}{2}\big) \Big)^{2}} ,
$$

and the time error has the opposite sign to the fourth-order spatial error, so at a Courant
fraction of 0.85 the two partly cancel. That is what makes order 4 work well at a large time
step, and why longer stencils only pay off with a smaller one: $2M$-th order in space with
second order in time is second order overall (Liu & Sen 2011). That is, as long as the time
error stays in the traces; by default it doesn't (below).

L2 misfit against the exact solution after about ten dominant wavelengths of travel
(1800 m), measured with `just accuracy`, with the time error kept
(`remove_time_dispersion = false`):

| Waves | Points per shortest wavelength | Order 4, Courant 0.85 | Order 8, 0.4 | Order 8, 0.3 | Order 16, 0.3 |
|---|---|---|---|---|---|
| P (explosion) | 4.0 | 19% | 8.3% | 4.8% | 5.5% |
| S (point force) | 3.5 | 47% | 6.4% | 7.1% | 2.4% |
| S (point force) | 4.6 | 15% | 2.3% | 1.3% | 1.4% |

Higher orders at a Courant fraction of 0.85 are *worse* than order 4 (order 8 there: 40%
for P). For the same 1.3% S misfit, order 4 needs about 8–9 points per wavelength (at 0.5),
which costs about 3.4 times the compute and 5.5 times the memory of order 8 at 4.6. Speeds on
4 cores: order 4 about 220 Mcell/s, 8 about 160, 12 about 130, 16 about 105.

**Removing the time error**, the default (see "Numerics"; Wang & Xu 2015; Koene et al.
2018). The traces carry the spatial error alone, so a large time step costs nothing, and the
longer stencils, which the time error held back, come into their own. Misfit at 1800 m with
Taylor coefficients, Courant 0.85 with the transform unless marked:

| Waves | Points per shortest wavelength | Order 8, 0.3, time error kept | Order 8 (default) | Order 12 | Order 16 |
|---|---|---|---|---|---|
| P (explosion) | 4.0 | 4.8% | 4.0% | 2.9% | 2.9% |
| S (point force) | 3.5 | 7.1% | 8.8% | 2.0% | 0.69% |
| S (point force) | 3.9 | 2.9% | 3.9% | 0.63% | 0.19% |
| S (point force) | 4.6 | 1.3% | 1.2% | 0.17% | 0.12% |

Order 8 is as accurate as with the time error kept at 0.3, in 2.8 times fewer steps (a
little less on the coarsest grids, where its time and space errors had partly cancelled).
Order 16 at 3.5 points per wavelength halves order 8's misfit at 4.6 for about 0.55 of its
compute and under half its memory (and about a fifth of the compute of order 8 at 0.3 with
the time error kept). Snapshots stay warped, and FWI needs the transposed transform in its
adjoint sources.

**Optimised coefficients**, also the default (see "Numerics"). With the time error gone,
the spatial error is all that is left, and coefficients chosen for the grid's shortest
waves rather than its longest cut it severalfold. Misfit at 1800 m, Courant 0.85, time error
removed (Taylor → optimised):

| Waves | Points per shortest wavelength | Order 8 | Order 12 | Order 16 |
|---|---|---|---|---|
| S (point force) | 3.0 | 19% → 7.3% | 6.5% → 0.89% | 2.9% → 0.42% |
| S (point force) | 3.5 | 8.8% → 2.1% | 2.0% → 0.28% | 0.69% → 0.15% |
| S (point force) | 3.9 | 3.9% → 0.61% | 0.63% → 0.13% | 0.19% → 0.11% |
| S (point force) | 4.6 | 1.2% → 0.19% | 0.17% → 0.12% | 0.12% → 0.12% |
| P (explosion) | 3.5 | 7.8% → 2.8% | 2.8% → 2.2% | 2.2% → 2.2% |
| P (explosion) | 4.0 | 4.0% → 3.0% | 2.9% → 3.0% | 2.9% → 2.9% |

Below about 0.1% something other than dispersion sets the floor. The P waves here are
coarser than an elastic grid would make them (they would have $v_p/v_s$ times the S wave's
points per wavelength), and they come out as good as with Taylor's coefficients or better.
Order 12 at 3.0 points per S wavelength is more accurate (0.89%) than order 8 with Taylor
coefficients at 4.6 (1.2%), for about a quarter of the compute and 0.28 of the memory. The
larger $\sum_n \lvert c_n\rvert$ of optimised coefficients costs 4–6% more time steps.

Rules of thumb:

| Purpose | Settings | Points per shortest wavelength |
|---|---|---|
| Waveform fits (FWI); the default | order 8, optimised, Courant 0.85, time error removed | $\ge 3.5$–$4$ |
| Waveform fits, cheapest | order 12–16, optimised, Courant 0.85, time error removed | $\ge 3$ |
| Quick looks | order 4, Courant 0.85 | 6–8 |
| Time error kept (wavefields, not traces) | order 8, Taylor, Courant 0.3 | $\ge 4.5$ |

What higher orders don't fix, from the literature and our measurements:

- **Interfaces.** Sharp interfaces sampled on the grid are only about second-order accurate
  and want a band-limited model (Mittet 2017).
- **The free surface.** Stress imaging is second order, so surface waves need finer grids
  (Kristek et al. 2002); Mittet (2002) gives a free surface suited to long stencils. Our
  free-surface test's 2.5% is mostly its reference, which doubles the P wave as if it were
  plane: on finer grids the misfit tends to about 2.2%, which leaves about 0.4% to the
  scheme at 8 points per P wavelength.
- **Absorbing layers.** Longer stencils reach further across the layers' inner edge:
  with the time error kept, reflection rises from $1.5 \times 10^{-4}$ of the peak at
  order 4 to $3.9 \times 10^{-4}$ (about $-68$ dB) at order 8. At the default it is
  $1.7 \times 10^{-4}$.
- **Source and receiver positions.** They snap to nodes rather than being placed exactly
  (Hicks 2002). A receiver 0.4 of a node off its position cost 26% misfit at 3.9 points per
  wavelength.

Optimised coefficients gain little while the time error remains, which is why they only
became worthwhile with the transform. Coefficients optimised for time and space together
(Liu & Sen 2011) only hold for one wave speed, which breaks down in heterogeneous models.

## Verification

`tests/physics.rs` checks whole simulations against independent answers (Aki & Richards
2002, §4.2–4.3), at the default (order 8, optimised coefficients, Courant 0.85, time error
removed) unless stated:

| Test | Result |
|---|---|
| Explosion in a full space against the analytic solution (near and far field), 100–250 m | 0.005–0.007% L2 misfit (Taylor coefficients 0.002–0.006%) |
| Point force against Stokes' solution: along, across and diagonal to the force | 0.001–0.005% (Taylor 0.001–0.008%) |
| The same with the time error kept, at Courant 0.3 | 0.03–0.09% |
| A P wave doubling at the free surface above a buried explosion | 2.5%, mostly the reference's (see "Accuracy") |
| Absorbing layers against a domain too large to reflect, without / with a free surface | 0.017% / 0.19% of peak |
| Reciprocity in a layered model with a sphere, one end on the free surface (orders 4 and 8) | $\le 8 \times 10^{-7}$ of peak |
| Dot-product test of the adjoint run, $\langle F w, \tilde d \rangle = \langle w, F^{T} \tilde d \rangle$, every component, rigid or absorbing walls (orders 4 and 8) | $\le 8.3 \times 10^{-8}$ |
| Gradients with respect to $v_p$, $v_s$ and $\rho$ against finite differences (Richardson-extrapolated), every node outside the layers perturbed (orders 4 and 8; and rigid walls with a water layer and uneven spacing) | $\le 2.8 \times 10^{-4}$ (order 8: $\le 1.6 \times 10^{-5}$) |
| 9 s with a water layer, soft sediment and a free surface (orders 4 and 8) | drains to $\le 3 \times 10^{-6}$ of peak |

## Towards elastic FWI

Planned next, following Fichtner (2011), ch. 7–11:

1. **Adjoint wavefield** (done: `src/adjoint.rs`). The residuals go back into the forward
   solver, stepped in reversed order, through the transpose of the receivers: a velocity
   residual as a unit force, a pressure residual as an isotropic stress source (Tarantola
   1984; Tromp et al. 2005; Fichtner et al. 2006). Because the scheme is self-adjoint, this
   is the exact transpose of the discrete forward run, not just of the wave equation: the
   dot-product test holds to $10^{-7}$ with the source and receivers clear of the absorbing
   layers.
2. **Gradients** (done: `src/gradient.rs`). Zero-lag correlations of the forward and
   adjoint wavefields for $\lambda$, $\mu$ and $\rho$, reparametrised to $v_p$, $v_s$ and
   $\rho$ (Köhn et al. 2012 on the trade-offs), taken through the solver's own material
   averaging and free surface, so they are the gradient of the discrete run: finite
   differences agree to $10^{-4}$ outside the absorbing layers. The forward field is replayed
   exactly from checkpoints for now.
3. **Forward wavefield during the adjoint run** (done: `gradient --reconstruct`,
   `Strategy::Reconstruct`). The interior scheme is time-reversible, so the field is rebuilt
   backwards (`Solver::step_back_with`) from a band saved around the layers, as deep as the
   stencil reaches (Dussaud et al. 2008; Aaker et al. 2020), instead of replayed from
   checkpoints. Saved for the whole run, the band needs more memory than checkpoints on
   small grids ($B N$ against $\sqrt N$) and less on large ones; saved for segments, each
   replayed from a checkpoint at its start, it halves a gradient's memory on
   `examples/layered.toml` (1.68 against 3.64 GiB) and cuts it by two thirds on a
   $1001 \times 1001 \times 401$ grid, for about 15% more time. It isn't bit-exact, so
   checkpoints stay the default.
4. **Error propagation.** How numerical error in the forward field (dispersion, the
   interface, free-surface and boundary errors above) and in the reconstruction carries into
   the adjoint field, the gradient and the inversion. The reconstruction's part is measured
   (research log, "Rebuilding the forward wavefield backwards"; `just drift`): it is
   single-precision round-off, $10^{-7}$–$10^{-6}$ of the gradient's largest value, growing
   while the rebuilt field is loud rather than with run length, order or time step, and the
   layers and the free surface add nothing to it. The tools: dot-product and reciprocity tests for adjoint
   consistency, gradient checks against finite differences, and the accuracy study repeated
   for gradients. With the time-dispersion transform, the adjoint source is the transposed
   transform of the residuals; since the transformed traces equal the time-continuous ones
   for every model, their gradient should too, which takes the time error out of the
   gradient as well. To be measured.
5. **Misfits and adjoint sources.** L2 waveform first, then cross-correlation time shifts and
   time–frequency phase misfits (Fichtner 2011, ch. 11).
6. **Optimisation** (L-BFGS done: `src/inversion.rs`). Preconditioned L-BFGS (Nocedal 1980;
   Liu & Nocedal 1989) on relative model changes, with bounds and an Armijo line search; each
   source a shot; the illumination (Shin et al. 2001) as preconditioner. Still to come:
   frequency continuation (Fichtner et al. 2013), receiver-side preconditioning, and
   parameter scaling for $\rho$.
7. Later: Born modelling and Hessian-vector products for truncated-Newton FWI (Métivier et
   al. 2013), source estimation, viscoelasticity.

## Licence

BSD 3-Clause; see [LICENSE](LICENSE).
