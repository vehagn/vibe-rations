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
| `[time]` | `duration`; `dt`, or `courant` (fraction of the stability limit; by default 0.3 at orders above 4, 0.85 otherwise) |
| `[scheme]` | `order` of the spatial derivatives: 2, 4, ..., 16 (default 8) |
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
  is at $\Delta t/2$.
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
  memory the run held, and how many receivers sit next to the layers (which makes the
  gradient approximate: move them a node further in)

```bash
just run truth.toml --out out/observed
just gradient guess.toml out/observed --balance --out out/gradient
```

`invert` runs the configuration's `[inversion]`: each source is a shot, whose observed
traces are in `<observed>/shot-NNN/` as `run --shots` writes them (or `<observed>/` itself for
one source). After every iteration it writes the model, `vp.npy`, `vs.npy` and `rho.npy`
(shape (nz, ny, nx), ready for `[model.files]`), and `inversion.json`: the settings and each
iteration's misfit, step and largest change. The time step and the layers' damping are held
at the starting model's, so give the observed run the same `[time] dt` (and `pml_vp_max`).

```bash
just run truth.toml --shots --out out/observed
just invert guess.toml --out out/inverted   # guess.toml: [inversion] observed = "out/observed"
```

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
Levander's $c_1 = \tfrac98$, $c_2 = -\tfrac1{24}$.

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

**Absorbing layers.** CPML (Komatitsch & Martin 2007) with $\kappa = 1$,
$\alpha = \pi f_0$, a quadratic damping profile $d(x)$ and the reflection coefficient of
Collino & Tsogka (2001). Along a free surface, CPMLs can grow without bound when surface
waves in a high-Poisson's-ratio medium run into them (Festa et al. 2005); the side layers then
also damp every field by $\exp(-p\, d\, \xi^3 \Delta t)$ per step, with $\xi \in [0, 1]$ the
depth into the layer and $p = 0.2$ (`pml_sponge`), in the spirit of the multiaxial PML of
Meza-Fajardo & Papageorgiou (2008).

**Sources and receivers.** Receivers interpolate staggered velocities onto nodes with cubic
weights $(-1, 9, 9, -1)/16$. Sources spread onto the staggered positions with the same
weights, divided by the cell's share of mass (half on the free surface). A source at a point
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
second order in time is second order overall (Liu & Sen 2011). Order 8 at 0.3 of the
stability limit is therefore the default.

L2 misfit against the exact solution after about ten dominant wavelengths of travel
(1800 m), measured with `just accuracy`:

| Waves | Points per shortest wavelength | Order 4, Courant 0.85 | Order 8, 0.4 | Order 8, 0.3 | Order 16, 0.3 |
|---|---|---|---|---|---|
| P (explosion) | 4.0 | 19% | 8.5% | 5.2% | 5.8% |
| S (point force) | 3.5 | 47% | 6.4% | 7.1% | 2.4% |
| S (point force) | 4.6 | 15% | 2.3% | 1.3% | 1.4% |

Higher orders at a Courant fraction of 0.85 are *worse* than order 4 (order 8 there: 40%
for P). For the same 1.3% S misfit, order 4 needs about 8–9 points per wavelength (at 0.5),
which costs about 3.4 times the compute and 5.5 times the memory of order 8 at 4.6. Speeds on
4 cores: order 4 about 220 Mcell/s, 8 about 160, 12 about 130, 16 about 105.

Rules of thumb:

| Purpose | Settings | Points per shortest wavelength |
|---|---|---|
| Waveform fits (FWI); the default | order 8, Courant 0.3 | $\ge 4.5$ |
| Quick looks | order 4, Courant 0.85 | 6–8 |
| Coarsest grids | order 12–16, Courant 0.3 | $\ge 3.5$ |

What higher orders don't fix, from the literature and our measurements:

- **Interfaces.** Sharp interfaces sampled on the grid are only about second-order accurate
  and want a band-limited model (Mittet 2017).
- **The free surface.** Stress imaging is second order, so surface waves need finer grids
  (Kristek et al. 2002). It is now the largest error in our free-surface test (2.4% at
  order 8, 1.9% at order 4); Mittet (2002) gives a free surface suited to long stencils.
- **Absorbing layers.** Longer stencils reach further across the layers' inner edge:
  reflection rises from $1.5 \times 10^{-4}$ to $3.9 \times 10^{-4}$ of the peak (about
  $-68$ dB) at order 8.
- **Source and receiver positions.** They snap to nodes rather than being placed exactly
  (Hicks 2002).

Optimised rather than Taylor coefficients (Holberg 1987) gain little at these error levels,
and coefficients optimised for time and space together (Liu & Sen 2011) only hold for one
wave speed, which breaks down in heterogeneous models.

## Verification

`tests/physics.rs` checks whole simulations against independent answers (Aki & Richards
2002, §4.2–4.3), at the default order 8 unless stated:

| Test | Result |
|---|---|
| Explosion in a full space against the analytic solution (near and far field), 100–250 m | 0.1–0.3% L2 misfit |
| Point force against Stokes' solution: along, across and diagonal to the force | 0.07–0.2% |
| A P wave doubling at the free surface above a buried explosion | 2.4% |
| Absorbing layers against a domain too large to reflect, without / with a free surface | 0.039% / 0.21% of peak |
| Reciprocity in a layered model with a sphere, one end on the free surface (orders 4 and 8) | $\le 1.5 \times 10^{-6}$ of peak |
| Dot-product test of the adjoint run, $\langle F w, \tilde d \rangle = \langle w, F^{T} \tilde d \rangle$, every component, rigid or absorbing walls (orders 4 and 8) | $\le 8.3 \times 10^{-8}$ |
| Gradients with respect to $v_p$, $v_s$ and $\rho$ against finite differences (Richardson-extrapolated), every node outside the layers perturbed (orders 4 and 8; and rigid walls with a water layer and uneven spacing) | $\le 2.2 \times 10^{-4}$ |
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
3. **Forward wavefield during the adjoint run.** The interior scheme is time-reversible, so
   it can be rebuilt backward from boundary values saved just inside the layers (Dussaud et
   al. 2008; Aaker et al. 2020) instead of being stored; checkpointing is the fallback.
4. **Error propagation.** How numerical error in the forward field (dispersion, the
   interface, free-surface and boundary errors above) and in the reconstruction (round-off in
   single precision, the non-reversible absorbing layers) carries into the adjoint field, the
   gradient and the inversion. The tools: dot-product and reciprocity tests for adjoint
   consistency, gradient checks against finite differences, and the accuracy study repeated
   for gradients.
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
