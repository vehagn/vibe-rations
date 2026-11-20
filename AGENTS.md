# vibe-rations, for coding agents

A 3D elastic wave simulator (fourth-order staggered-grid velocity–stress finite differences,
CPML, free surface), built as the forward engine for elastic full waveform inversion. The
[README](README.md) explains the numerics and how they are verified.

## Commands

Run these in the dev container (`just dev` on the host opens a shell in it). Nothing is
installed on the host.

| Command | Does |
|---|---|
| `just test` | All tests, optimised (`cargo test --release`) |
| `just lint` | rustfmt check and clippy with `-D warnings` |
| `just check <config>` / `just run <config>` | Set up / run a simulation |

`just test` and `just lint` must pass.

## Map

```
src/
  main.rs        the CLI: run | check
  config.rs      the TOML configuration (serde, deny_unknown_fields)
  grid.rs        the grid and its two-node halo, indexing
  model.rs       vp, vs, rho on the nodes: layers, .npy files, spheres
  solver.rs      fields, staggered material, the update kernels, free surface, PointOp
  fd.rs          fourth-order staggered derivatives along one axis, a row at a time
  cpml.rs        CPML profiles and memory slabs, the sponge factors
  source.rs      wavelets, moment tensors, how sources enter the grid
  receiver.rs    recording traces through PointOps
  snapshot.rs    sections and volumes during a run
  simulation.rs  a run set up from a configuration; outputs
  image.rs, npy.rs
tests/physics.rs whole simulations against analytic solutions, reciprocity, stability
examples/        annotated configurations
```

## Conventions

- Keep sources the exact transpose of receivers (`Solver::stagger_weights`, `PointOp`):
  adjoint-based FWI depends on it, and the reciprocity test checks it.
- Any change to the update loops, the free surface or the CPML must keep every test in
  `tests/physics.rs` passing; don't loosen a tolerance to make a change pass.
- Physical units are SI; z points down; arrays are x-fastest.
