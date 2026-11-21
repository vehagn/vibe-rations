---
name: run-simulation
description: Set up and run a vibe-rations simulation - write the TOML configuration, size the grid and time step, run it, and read the traces, gathers and snapshots it writes. Use when asked to model a scenario, make an example, or produce seismograms or wavefield pictures.
---

# Running a simulation

## 1. Start from an example

Copy the closest of `examples/explosion.toml` (homogeneous half-space, explosion),
`examples/layered.toml` (layers, double-couple earthquake) or `examples/sphere-inversion.toml`
(two force shots and an `[inversion]` section, small enough to invert in a browser) into
`out/` for experiments, or
into `examples/` for a new example (`tests/examples.rs` then requires it to set up without
warnings). Every key is documented in `src/config.rs`; unknown keys are errors.

Positions are metres from the first node; z points down; with `free_surface = true`
(default) z = 0 is the surface.

## 2. Size the grid

Check these before running; `just in check <config>` prints all of them:

- **Points per wavelength.** At the default (eighth order, optimised coefficients, time
  error removed), 3.5–4 for about 1% misfit after ten wavelengths (it warns below 3.5):
  $h \le v_{\text{slow}} / (4 \cdot 2.5 f_0)$, with $v_{\text{slow}}$ the slowest S speed (P in
  fluids) and $f_0$ the wavelet frequency. Orders 12–16 get there at 3 (it warns below 3);
  with Taylor coefficients order 8 needs 4.5, and fourth order 8–9 (README, "Accuracy").
  The coefficients are designed for the run's own points per wavelength, which depends on
  the slowest wave; an inversion should pin it with `[scheme] optimised_for`.
- **Absorbing layers** of `pml_width` nodes (default 20, at least 10) line the sides and
  bottom. They are inside the grid, so add 2 × `pml_width` nodes to the region of interest
  and keep sources and receivers out of them (it warns).
- **Time step:** automatic at a fraction $C = 0.85$ of the limit,
  $\Delta t \le 1 / \big(v_{p,\max} \sum_n |c_n| \sqrt{\Delta x^{-2} + \Delta y^{-2} + \Delta z^{-2}}\big)$.
  The fastest $v_p$ sets it. That large a step is only accurate because leapfrog's time error
  is taken out of the traces by default (`[scheme] remove_time_dispersion`; Wang & Xu 2015,
  Koene et al. 2018): traces are final only after the last step, and snapshots stay warped.
  With `remove_time_dispersion = false` the step drops to 0.3 at orders above 4; don't raise
  it there, accuracy collapses. Orders 12–16 with the transform and optimised coefficients
  are the most accurate per unit of compute: order 12 at 3.0 points per wavelength beats
  order 8 with Taylor coefficients at 4.6, for a quarter of the compute (README, "Accuracy").
- **Cost:** `just in memory <config>` gives the memory part by part without building
  anything (add `--gradient` for a gradient run, which needs many times more: checkpoints
  and a replayed segment, about $42\sqrt N$ bytes per padded node for $N$ steps; `--invert`
  for an inversion, which adds L-BFGS's memory), and `run`/`check`/`gradient`/`invert`
  refuse runs larger than the memory available (the container
  has about 7 GiB). Roughly $4 \cdot 17 \cdot (n_x + 2L)(n_y + 2L)(n_z + 2L)$ bytes plus the
  model, the layers and the traces, with $L$ the stencil half-length (4 at order 8). Time
  $\approx$ nodes $\times$ steps $/\, 1.6 \times 10^8$ seconds on the container's 4 cores at
  order 8.

## 3. Sources and wavelets

One of `explosion`, `moment = { xx, yy, zz, xy, xz, yz }`, `double_couple = { strike, dip,
rake, m0 }` (x north, y east, z down) or `force = [fx, fy, fz]`. The wavelet is the moment
rate (or force): a Ricker peaks at 1, so its amplitude is a peak rate; a Gaussian has unit
area, so its amplitude is the final moment. A Gaussian leaves a static stress around the
source; avoid `p` snapshots near it.

## 4. Run

```bash
just in run out/my-run.toml     # progress every 2 s; --quiet to silence it
```

It stops with an error if the wavefield blows up; see the `debug-instability` skill.

## 5. Read the results

In `[output] dir`:

| File | Content |
|---|---|
| `vx.npy` … `p.npy` | float32, shape (receivers, samples) |
| `time.npy` | sample times; the first is dt/2 (velocities lag stresses by half a step) |
| `receivers.csv` | receiver positions after snapping to nodes |
| `run.json` | dt, steps, accuracy figures, warnings |
| `*_gather.png` | a column per receiver, time down, each trace scaled to its own peak |
| `snapshots/<field>_<plane><index>_<step>.png` | blue negative, red positive, colour $\propto \sqrt{|u|}$; grey shows the model, darker in the layers |

Look at the PNGs with your image-reading tool, a mid-run snapshot and the last one at
least. Check that waves leave through the darker layers without coming back, and that the
arrival times in the gather fit distance / speed.

```python
import numpy as np
vz, t = np.load("out/explosion/vz.npy"), np.load("out/explosion/time.npy")
```
