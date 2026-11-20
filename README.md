# vibe-rations

A 3D elastic seismic wave simulator in Rust, built as the forward engine for elastic full
waveform inversion (FWI).

It solves the isotropic elastic wave equation in velocity–stress form on a staggered grid,
fourth order in space and second order in time (Virieux 1986; Levander 1988; Graves 1996),
with convolutional PML absorbing layers (Komatitsch & Martin 2007) and a free surface. The
update loops run in parallel across z-planes with rayon.

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
| `just run <config>` | Runs a simulation |
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
don't pass silently. Positions are in metres from the first node, z points down, and with a
free surface z = 0 is the surface.

| Section | What it sets |
|---|---|
| `[grid]` | `points = [nx, ny, nz]`, `spacing` (one value or `[dx, dy, dz]`) |
| `[time]` | `duration`; `dt`, or `courant` (fraction of the stability limit, default 0.85) |
| `[boundary]` | `pml_width` (default 20 nodes), `free_surface` (default true), PML tuning |
| `[[model.layer]]` | Horizontal layers: `top` depth, `vp`, `vs`, `rho`. `vs = 0` is a fluid |
| `[model.files]` | Or `.npy` arrays of shape (nz, ny, nx) for `vp`, `vs`, `rho` |
| `[[model.sphere]]` | Spheres painted over the model |
| `[[source]]` | `position`, `wavelet` (Ricker or Gaussian), and one of `explosion`, `moment`, `double_couple` (strike, dip, rake, m0) or `force` |
| `[receivers]` | `components` (`vx`, `vy`, `vz`, `p`), `[[receivers.line]]`s and `[[receivers.point]]`s |
| `[[snapshot]]` | `field`, `plane` (`xy`, `xz`, `yz` or `volume`), `at`, `interval`, `formats` (`png`, `npy`) |
| `[output]` | `dir`, `record_every` |

The wavelet is the moment rate of a moment-tensor source and the force of a force source. A
Ricker wavelet peaks at 1, so its amplitude is a peak rate; a Gaussian has unit area, so its
amplitude is the final moment. A source with a non-zero final moment leaves a static stress
around it, which dominates pressure snapshots near the source. Velocity has no static part.

## Output

In the output directory:

- `vx.npy`, `vy.npy`, `vz.npy`, `p.npy`: traces, shape (receivers, samples), float32
- `time.npy`: the sample times. Velocities lag stresses by half a step, so the first sample
  is at dt/2.
- `receivers.csv`: the receivers' positions after snapping to nodes
- `*_gather.png`: one column per receiver, each trace scaled to its own peak
- `run.json`: the run's settings, time step, accuracy figures and warnings
- `snapshots/`: sections as PNG (blue negative, red positive, colour ∝ √amplitude, over the
  model in grey with the absorbing layers darker) or `.npy`

```python
import numpy as np
vz, t = np.load("out/explosion/vz.npy"), np.load("out/explosion/time.npy")
```

## Numerics

**Grid.** The interior nodes are padded by a two-node halo for the stencils. Normal stresses,
λ and μ sit on the nodes; vx, vy and vz half a node along x, y and z; shear stresses half a
node along both their axes. Buoyancy is averaged arithmetically onto the velocity positions,
and μ harmonically onto the shear positions, so shear stress vanishes next to a fluid.

**Stability and accuracy.** The time step must satisfy dt ≤ 1 / (vp_max (9/8 + 1/24)
√(1/dx² + 1/dy² + 1/dz²)). `just check` warns if the grid has fewer than 5 points per
shortest wavelength, taking the slowest wave (S, or P in fluids) at 2.5 times the dominant
frequency.

**Free surface.** Stress imaging on the σzz plane: σzz = 0 there, σxx and σyy are updated
with ∂z vz = -λ/(λ+2μ)(∂x vx + ∂y vy), and σzz, σxz, σyz mirror oddly across the surface.
Velocities mirror evenly. With that choice the discrete z-derivatives of stress and of
velocity stay negative adjoints of each other next to the surface, so the scheme conserves
energy there. Every stencil stays fourth order.

**Absorbing layers.** CPML with κ = 1, α = π f_dominant and a quadratic profile; reflection
measured at -76 dB (0.015% of peak). Along a free surface, CPMLs can grow without bound when
surface waves in a high-Poisson's-ratio medium run into them (Festa et al. 2005). The side
layers then also apply a small plain damping, 0.2 · d · x³ (cf. the multiaxial PML of
Meza-Fajardo & Papageorgiou 2008). It keeps them stable at a cost of about -55 dB
reflection, and `pml_sponge` tunes it.

**Sources and receivers.** Receivers interpolate staggered velocities onto nodes with cubic
weights (-1, 9, 9, -1)/16. Sources spread onto the staggered positions with the same
weights, divided by the cell's share of mass (half on the free surface). A source at a point
is therefore exactly the transpose of a receiver there: what FWI needs to back-propagate
residuals.

## Verification

`tests/physics.rs` checks whole simulations against independent answers:

| Test | Result |
|---|---|
| Explosion in a full space against the analytic solution (near and far field), 100–250 m | 0.7–1.6% L2 misfit |
| Point force against Stokes' solution: along, across and diagonal to the force | 0.3–0.8% |
| A P wave doubling at the free surface above a buried explosion | 1.9% |
| Absorbing layers against a domain too large to reflect, without / with a free surface | 0.015% / 0.18% of peak |
| Reciprocity, G_ij(B, A) = G_ji(A, B), layered model with a sphere, one end on the free surface | 1e-6 of peak |
| 9 s with a water layer, soft sediment and a free surface | drains to 4e-6 of peak |

## Towards elastic FWI

Planned next, following Fichtner (2011), *Full Seismic Waveform Modelling and Inversion*,
ch. 7–11:

1. **Adjoint wavefield.** Back-propagate the residuals by injecting them as forces at the
   receivers (the transpose above) in reversed time.
2. **Gradients.** Zero-lag correlations of the forward and adjoint wavefields for λ, μ and
   ρ, reparametrised to vp, vs and ρ (see Köhn et al. 2012 on the trade-offs).
3. **Forward wavefield during the adjoint run.** The interior scheme is time-reversible, so
   it can be rebuilt backward from boundary values saved just inside the layers (Dussaud et
   al. 2008) instead of being stored; checkpointing is the fallback.
4. **Misfits and adjoint sources.** L2 waveform first, then cross-correlation time shifts and
   time–frequency phase misfits (Fichtner ch. 11).
5. **Optimisation.** L-BFGS with a line search, frequency continuation (multiscale), and
   gradient checks against finite differences.
6. Later: Born modelling and Hessian-vector products for truncated-Newton FWI (Métivier et
   al. 2013), source estimation, viscoelasticity.

## Licence

BSD 3-Clause; see [LICENSE](LICENSE).
