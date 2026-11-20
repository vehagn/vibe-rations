---
name: debug-instability
description: Find out why a vibe-rations run blows up, grows slowly, or doesn't drain - locate where the energy grows and isolate the cause (free surface, absorbing layers, fluids, time step, model). Use when a run stops with "the wavefield blew up" or the long-run test fails.
---

# Debugging an instability

Explicit staggered schemes fail in a few known ways. Don't guess: find where the energy
grows, then switch features off one at a time.

## 1. Where and when does it grow?

Write a throwaway test, `tests/scratch.rs`, marked `#[ignore]`, and **delete it before
committing**. Step the simulation yourself and log the largest velocity with its node
every few hundred steps:

```rust
let mut sim = Simulation::new(&Config::from_toml(TEXT).unwrap(), Path::new(".")).unwrap();
let g = sim.solver.grid;
for n in 0..sim.steps {
    sim.step();
    if n % 500 == 499 {
        let f = &sim.solver.fields;
        let (mut best, mut at) = (0f32, (0, 0, 0, ""));
        for (name, v) in [("vx", &f.vx), ("vy", &f.vy), ("vz", &f.vz)] {
            for k in 0..g.n[2] { for j in 0..g.n[1] { for i in 0..g.n[0] {
                let x = v[g.idx(i, j, k)].abs();
                if !(x <= best) { best = x; at = (i, j, k, name); } // NaN wins too
            }}}
        }
        println!("{n}: {best:.1e} at {at:?}");
    }
}
```

```bash
just in test --test scratch -- --ignored --nocapture
```

Run long enough: 5–10 s of simulated time. Some instabilities only start after a couple
of seconds. A healthy run with absorbing layers decays by many orders of magnitude.

## 2. Read the location

| Grows at | Suspect |
|---|---|
| k = 0–2, inside the x/y layers, near corners | CPML meeting surface waves (high vp/vs): the sponge |
| k = 0–2 everywhere | free-surface update or mirroring |
| inside a layer, away from the surface | CPML profiles or memory indexing |
| at an interface or fluid | material averaging (harmonic μ, arithmetic ρ) |
| everywhere at once, from the start | the time step (CFL), or the update order |

## 3. Isolate by switching features off

Run the same model in variants, one change each:

- `pml_width = 0` (rigid walls). If the energy then stays bounded and doesn't grow, the
  interior and the free surface are fine and the layers are to blame.
- `free_surface = false`.
- A harder model: raise vs (lower vp/vs), or remove fluid layers.
- `courant = 0.5`. If that helps, the step is too long or the limit is wrong.
- Layer settings: `pml_sponge`, `pml_alpha_frequency`, `pml_kappa_max`.

## 4. What it was before

The one found so far: a soft layer (vp/vs = 3.3) with a free surface and CPML grew from
~2 s in the surface corners of the side layers. Rigid walls were stable; α, κ and dt
didn't help. The cure was plain damping in the side layers along the free surface
(`pml_sponge`, AGENTS.md has the measurements). Any fix must keep
`absorbing_layers_reflect_little*` and the reciprocity test within their baselines; run
the `verify-solver` checks afterwards.
