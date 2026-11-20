---
name: run-simulation
description: Set up and run a vibe-rations simulation - write the TOML configuration, size the grid and time step, run it, and read the traces, gathers and snapshots it writes. Use when asked to model a scenario, make an example, or produce seismograms or wavefield pictures.
---

# Running a simulation

## 1. Start from an example

Copy the closest of `examples/explosion.toml` (homogeneous half-space, explosion) or
`examples/layered.toml` (layers, double-couple earthquake) into `out/` for experiments, or
into `examples/` for a new example (`tests/examples.rs` then requires it to set up without
warnings). Every key is documented in `src/config.rs`; unknown keys are errors.

Positions are metres from the first node; z points down; with `free_surface = true`
(default) z = 0 is the surface.

## 2. Size the grid

Check these before running; `just in check <config>` prints all of them:

- **Points per wavelength ≥ 5.** h ≤ v_slow / (5 · 2.5 · f0), with v_slow the slowest S
  speed (P in fluids) and f0 the wavelet frequency. It warns below 5.
- **Absorbing layers** of `pml_width` nodes (default 20, at least 10) line the sides and
  bottom. They are inside the grid, so add 2 × `pml_width` nodes to the region of interest
  and keep sources and receivers out of them (it warns).
- **Time step:** automatic at `courant` (0.85) of the limit
  dt ≤ 1 / (vp_max · 7/6 · √(1/dx² + 1/dy² + 1/dz²)). The fastest vp sets it.
- **Cost:** memory ≈ 4 bytes × 17 × (nx+4)(ny+4)(nz+4) plus the layers' memory; time ≈
  nodes × steps / 2e8 seconds on the container's 4 cores. Keep memory under ~6 GiB.

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
| `snapshots/<field>_<plane><index>_<step>.png` | blue negative, red positive, colour ∝ √amplitude; grey shows the model, darker in the layers |

Look at the PNGs with your image-reading tool, a mid-run snapshot and the last one at
least. Check that waves leave through the darker layers without coming back, and that the
arrival times in the gather fit distance / speed.

```python
import numpy as np
vz, t = np.load("out/explosion/vz.npy"), np.load("out/explosion/time.npy")
```
