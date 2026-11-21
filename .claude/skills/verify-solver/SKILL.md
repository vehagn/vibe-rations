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

It prints a line per check. Compare with the baselines in the Notes of `AGENTS.md` (these
are at the default: order 8, optimised coefficients, $C = 0.85$, the time error removed from
the traces):

| Test | Baseline | Breaks when |
|---|---|---|
| `an_explosion_matches_the_full_space_solution` | 0.005–0.007% | stencils, stress injection, timing, receiver interpolation, the time-dispersion transform |
| `a_point_force_matches_stokes_solution` | 0.001–0.005% | velocity injection, buoyancy, S waves, shear averaging, the transform, the optimised coefficients |
| `a_point_force_matches_stokes_solution_with_taylor_coefficients` | 0.001–0.008% | the same with Taylor coefficients |
| `a_point_force_matches_stokes_solution_without_the_transform` | 0.03–0.09% | the same, with the time error kept (opt-out path, $C = 0.3$) |
| `the_free_surface_doubles_a_vertically_incident_p_wave` | 2.52% (about 2.2% is the reference's) | free-surface condition, mirroring, surface sampling |
| `a_planar_interface_reflects_and_transmits_plane_waves` | P direct 0.010%, reflected 0.72%, transmitted 0.078%; S 0.026%, 7.5%, 0.44%; a layer's top on a node reflects $h/v$ early | material averaging (ρ at the velocities, μ harmonically), how layers meet the nodes |
| `absorbing_layers_reflect_little` (± free surface) | $1.7 \times 10^{-4}$ / $1.9 \times 10^{-3}$ of peak | CPML profiles or memory, the sponge |
| `swapping_source_and_receiver_gives_the_same_trace` | $\le 8 \times 10^{-7}$ | anything making the scheme non-self-adjoint, or sources no longer the transpose of receivers |
| `the_adjoint_run_is_the_transpose_of_the_forward_run` | $\le 8.3 \times 10^{-8}$ | the adjoint's step pairing or `Solver::adjoint_op`, or anything that breaks reciprocity; a receiver whose interpolation reaches into a layer (the helper refuses one) |
| `the_inversion_gradient_matches_finite_differences` | $\le 3 \times 10^{-5}$ (measured $4.0 \times 10^{-6}$ at $h = 10^{-2}$) | the inversion's parametrisation, its sum over shots, the pinned operator, or anything the gradient checks catch; a receiver reading the layers (the test refuses one) |
| `the_gradient_matches_finite_differences` (and `…_with_walls_and_water`) | $\le 2.8 \times 10^{-4}$ ($\le 1.6 \times 10^{-4}$) extrapolated; order 8 $\le 1.6 \times 10^{-5}$ | the material averaging, the free surface's moduli, the update's strain rates, sources or receivers that change how the material enters, the transform's transpose (`Simulation::recorded_residuals`) |
| `a_long_run_stays_stable_and_drains` | $\le 3 \times 10^{-6}$ | instabilities (see `debug-instability`) |

A number that got worse is a regression even if the test passes. Find out why before
going on. Never loosen a tolerance to get a change through. If a change makes a number
better, tighten the tolerance and update the baseline in `AGENTS.md`. Record what you
measured in `docs/research-log.md`.

Changes to the stencils or the time stepping also need `just in accuracy`: compare its table
with the README's "Accuracy" section. Physics tests that exercise the free surface and
stability (`swapping_source_and_receiver_gives_the_same_trace`,
`a_long_run_stays_stable_and_drains`) run at orders 4 and 8; keep both.

## 2. Reading misfits

- A misfit that **grows with distance** is phase error: dispersion (too few points per
  wavelength, or too large a time step for the order) or a timing offset. Velocities are at
  $(n + \tfrac12)\Delta t$, forces enter at $n \Delta t$ and moments at
  $(n + \tfrac12)\Delta t$. The time-dispersion transform depends on these times
  (`Source::offset`, `Simulation::trace_warp`).
- **The same misfit at every distance** is amplitude or position error: interpolation
  weights, the source scaling (dt/V, buoyancy), or a half-node offset.
- **Exactly a factor 2 or a sign:** half-cell shares on the free surface
  (`Solver::cell_share`), the mirror parities, or the moment-tensor sign
  ($\sigma_{ij} \mathrel{-}= \Delta t\, \dot{M}_{ij} / V$).
- **Reciprocity off at 1e-3 or more:** an asymmetry between injection and sampling, or
  between the operators that differentiate stress and velocity.

## 3. Adding a physics test

Put it in `tests/physics.rs`; build the configuration with `homogeneous(...)` and run it
with `simulate(...)`. Keep receivers ~10 nodes clear of the absorbing layers, and stop
before reflections could return unless they are the point. Use at least 8 points per
wavelength so the reference, not dispersion, sets the tolerance. Print the measured
number (`println!`, shown by `just physics`) and set the tolerance just above it.

Closed-form references already in the file (Aki & Richards 2002, §4.2–4.3):

- Explosion with moment rate $m_0 w(t)$: radial velocity
  $v_r = \frac{m_0}{4\pi\rho\alpha^2}\big[\frac{w(\tau)}{r^2} + \frac{w'(\tau)}{\alpha r}\big]$,
  $\tau = t - r/\alpha$ (`explosion_velocity` in `tests/common`).
- Point force $f\,w(t)$ along $x$: Stokes' solution, differentiated in time, with the near-field
  integral by the trapezium rule (`force_velocity`).

Plane waves at normal incidence on an interface (`common::plane`): a uniform plane of
sources in a grid whose sides copy the nearest node's values, so the field stays exactly 1D
and the exact answer is $R = (Z_1 - Z_2)/(Z_1 + Z_2)$, $T = 1 + R$. `just interfaces` sweeps
the interface's position between nodes, how layers are sampled and the order.

Good next references: Lamb's problem (a force on the free surface; Pekeris 1955 for a
Poisson solid), a fluid–solid interface (Scholte waves), and reflection at oblique incidence,
where P and S convert.
