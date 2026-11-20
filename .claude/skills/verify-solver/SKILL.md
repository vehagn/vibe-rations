---
name: verify-solver
description: Verify changes to the vibe-rations numerics - the update kernels, absorbing layers, free surface, material averaging, sources or receivers - against analytic solutions, reciprocity and stability, and add new physics tests. Use whenever a change could alter the wavefield, before committing it.
---

# Verifying the solver

Unit tests show that the pieces do what they say. Only whole simulations compared with
independent answers show the physics is right. Any change to `solver.rs`, `cpml.rs`,
`fd.rs`, `source.rs` or the sampling in `receiver.rs` needs this.

## 1. Run the physics tests and compare with the baselines

```bash
just in physics
```

It prints a line per check. Compare with the baselines in the Notes of `AGENTS.md`:

| Test | Baseline | Breaks when |
|---|---|---|
| `an_explosion_matches_the_full_space_solution` | 0.7–1.6% | stencils, stress injection, timing, receiver interpolation |
| `a_point_force_matches_stokes_solution` | 0.3–0.8% | velocity injection, buoyancy, S waves, shear averaging |
| `the_free_surface_doubles_a_vertically_incident_p_wave` | 1.9% | free-surface condition, mirroring, surface sampling |
| `absorbing_layers_reflect_little` (± free surface) | 0.015% / 0.18% | CPML profiles or memory, the sponge |
| `swapping_source_and_receiver_gives_the_same_trace` | ≤ 1e-6 | anything making the scheme non-self-adjoint, or sources no longer the transpose of receivers |
| `a_long_run_stays_stable_and_drains` | 4e-6 | instabilities (see `debug-instability`) |

A number that got worse is a regression even if the test passes. Find out why before
going on. Never loosen a tolerance to get a change through. If a change makes a number
better, tighten the tolerance and update the baseline in `AGENTS.md`.

Changes to the stencils or the time stepping also need `just in accuracy`: compare its table
with the README's "Accuracy" section. Physics tests that exercise the free surface and
stability (`swapping_source_and_receiver_gives_the_same_trace`,
`a_long_run_stays_stable_and_drains`) run at orders 4 and 8; keep both.

## 2. Reading misfits

- A misfit that **grows with distance** is phase error: dispersion (too few points per
  wavelength) or a timing offset. Velocities are at (n + ½)·dt, forces enter at n·dt and
  moments at (n + ½)·dt.
- **The same misfit at every distance** is amplitude or position error: interpolation
  weights, the source scaling (dt/V, buoyancy), or a half-node offset.
- **Exactly a factor 2 or a sign:** half-cell shares on the free surface
  (`Solver::cell_share`), the mirror parities, or the moment-tensor sign (σ -= dt·Ṁ/V).
- **Reciprocity off at 1e-3 or more:** an asymmetry between injection and sampling, or
  between the operators that differentiate stress and velocity.

## 3. Adding a physics test

Put it in `tests/physics.rs`; build the configuration with `homogeneous(...)` and run it
with `simulate(...)`. Keep receivers ~10 nodes clear of the absorbing layers, and stop
before reflections could return unless they are the point. Use at least 8 points per
wavelength so the reference, not dispersion, sets the tolerance. Print the measured
number (`println!`, shown by `just physics`) and set the tolerance just above it.

Closed-form references already in the file (Aki & Richards 2002, §4.2–4.3):

- Explosion with moment rate m0·w(t): radial velocity
  v = m0/(4πρα²)·[w(τ)/r² + w'(τ)/(α r)], τ = t − r/α (`explosion_velocity`).
- Point force f·w(t) along x: Stokes' solution, differentiated in time, with the near-field
  integral by the trapezium rule (`force_velocity`).

Good next references: Lamb's problem (a force on the free surface; Pekeris 1955 for a
Poisson solid), a fluid–solid interface (Scholte waves), and plane-wave reflection
coefficients at a layer interface.
