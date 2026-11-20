# Research log

What we tried, what we measured and what we concluded, in order, as raw material for a paper
on this simulator. Each entry says how to reproduce it. Citations are keys of
[`references.bib`](references.bib), written `[@key]`. Math is LaTeX.

Notation: $v_p$, $v_s$ wave speeds; $h$ grid spacing; $\Delta t$ time step; $C$ the Courant
fraction, $\Delta t = C \, \Delta t_{\max}$; $f_0$ the wavelet's dominant frequency; points
per wavelength (ppw) count the shortest wavelength, $v_{\min} / (2.5 f_0)$.

---

## 2026-11-20 — Verification of the forward solver

**Setup.** Velocity–stress staggered grid [@virieux1986], fourth-order Taylor stencils
[@levander1988], leapfrog in time, CPML [@komatitsch2007] with the reflection coefficient of
[@collino2001], stress-imaging free surface [@levander1988; @graves1996].

**Tests** (`tests/physics.rs`, `just physics`), against closed-form solutions for a
homogeneous full space [@aki2002]:

- an explosion with moment rate $m_0 w(t)$, radial velocity
  $v_r = \frac{m_0}{4\pi\rho v_p^2}\big[\frac{w(\tau)}{r^2} + \frac{w'(\tau)}{v_p r}\big]$,
  $\tau = t - r/v_p$, at $r = 100$–$250$ m;
- a point force (Stokes' solution, differentiated in time) along, across and diagonal to it;
- the doubling of a vertically incident P wave at the free surface;
- absorbing-layer reflections, as the difference from a domain too large to reflect;
- reciprocity, $G_{ij}(\mathbf{x}_B, \mathbf{x}_A) = G_{ji}(\mathbf{x}_A, \mathbf{x}_B)$, in a
  layered model with a sphere, one end on the free surface;
- a 9 s run with water, soft sediment and a free surface, which must drain.

**Results, fourth order, $C = 0.85$:** explosion 0.7–1.6% L2 misfit, force 0.3–0.8%, free
surface 1.9%, layer reflection $1.5 \times 10^{-4}$ of the peak ($-76$ dB) without a free
surface, reciprocity $\le 10^{-6}$. The eighth-order results are in the 2026-11-20
accuracy entry below.

## 2026-11-20 — Interpolating receivers and spreading sources

Averaging two staggered values onto a node attenuates a wavenumber $k$ by
$\cos(kh/2)$: about 2% at eight ppw, which dominated the first misfits (2.5–3.7%). Cubic
weights $(-1, 9, 9, -1)/16$ cut them about threefold (force: 0.28–0.76%). Sources use the same
weights, so injection is the exact transpose of sampling.

## 2026-11-20 — A self-adjoint free surface

**Problem.** The first free surface reduced the $z$-derivatives of velocity next to the
surface to second order. That pair of operators is not the negative adjoint of the
stress-imaging operators, so the scheme no longer conserves a discrete energy there; the
free-surface misfit was 7.6%.

**Fix.** Mirror stresses oddly ($\sigma_{zz}$, $\sigma_{xz}$, $\sigma_{yz}$) and velocities
evenly across the surface. With the surface nodes carrying half-cell weights in the energy
inner product, the discrete operators are then negative adjoints for any stencil length, and
every stencil keeps its full order. The free-surface misfit fell to 1.9%. Values on the
surface plane stand for half a cell, so sources there are divided by $\tfrac12$; without that,
reciprocity with one end on the surface was off by exactly 2.

**Reproduce:** `swapping_source_and_receiver_gives_the_same_trace` (orders 4 and 8).

## 2026-11-20 — CPML instability along a free surface

**Observation.** With a soft layer ($v_p/v_s \approx 3.3$), a free surface and CPML side
layers, the field grew exponentially from about 2 s, starting in the layers' corners at the
surface, as [@festa2005] describe. With rigid walls instead of layers the energy stayed
bounded; a larger $\alpha$, $\kappa_{\max} = 4$ and a smaller $\Delta t$ did not help.

**Fix.** Damp every field in the side layers by $\exp(-p\, d\, \xi^3 \Delta t)$ per step, only
with a free surface, cf. the multiaxial PML [@mezafajardo2008]. Reflection against stability,
measured 9 nodes from the layer:

| Profile | $p$ | Reflection (% of peak) | Stable |
|---|---|---|---|
| $d$ | 0.1 | 0.23 | yes |
| $d\,\xi^2$ | 0.1 | 0.06–0.15 | yes |
| $d\,\xi^4$ | 0.05 | — | no |
| $d\,\xi^3$, side layers only | 0.2 | 0.07–0.18 | yes (chosen) |

**Reproduce:** `a_long_run_stays_stable_and_drains`; the `debug-instability` skill.

## 2026-11-20 — Stencil order and time step

**Question.** How much accuracy do longer staggered stencils buy, and at what cost?

**Theory.** For a plane wave in a homogeneous medium the scheme's dispersion relation is

$$
\sin\!\Big(\frac{\omega \Delta t}{2}\Big) = \frac{v \Delta t}{h} \sqrt{\sum_{i} K_i^2},
\qquad K_i = \sum_{n=1}^{L} c_n \sin\big((2n - 1) \tfrac{k_i h}{2}\big),
$$

for P and S alike (`python3 docs/analysis/dispersion.py` and `docs/analysis/time_space.py`
reproduce the numbers below). The leapfrog error makes waves fast and the spatial error makes them slow,
so at fourth order and $C = 0.85$ they partly cancel. Requiring the S group-velocity error to
stay below 0.5% in every direction (von Neumann analysis, $v_p/v_s = \sqrt3$, worst of axis and
diagonal):

| Taylor stencil | ppw at $C = 0.85$ | ppw at $C = 0.5$ | ppw at $C = 0.25$ |
|---|---|---|---|
| 4 points | 7.7 | 8.7 | 9.1 |
| 8 points | 7.0 | 4.6 | 4.9 |
| 16 points | 6.6 | 3.9 | 3.5 |

So longer stencils pay off only with a smaller time step, as [@liu2011] note: $2M$-th order in
space with second order in time is second order overall. Least-squares-optimised spatial
coefficients (in the manner of [@holberg1987]) gained 0.2–0.8 ppw at the 0.5% tolerance and did
worse at stricter ones. Coefficients fitted to
the time–space dispersion relation [@liu2011] reached 4.5 ppw at eighth order for the speed they
were fitted to, but needed 9.9 ppw at half that speed and failed for P waves on the same grid,
so they don't suit heterogeneous models.

**Measurement** (`tests/accuracy.rs`, `just accuracy`). L2 misfit against the closed-form
solutions at 300–1800 m, i.e. up to about ten dominant wavelengths:

| Waves, ppw | Order 4, $C = 0.85$ | Order 8, $C = 0.4$ | Order 8, $C = 0.3$ | Order 16, $C = 0.3$ |
|---|---|---|---|---|
| P, 4.0 | 19% | 8.5% | 5.2% | 5.8% |
| S, 3.5 | 47% | 6.4% | 7.1% | 2.4% |
| S, 4.6 | 15% | 2.3% | 1.3% | 1.4% |

Order 8 at $C = 0.85$ is *worse* than order 4 (40% for P). Order 4 needs about 8–9 ppw (at
$C = 0.5$; 0.83% at 9.2 ppw) for what order 8 at $C = 0.3$ does at 4.6 ppw (1.26%). With
measured speeds (Mcell/s on 4 cores: order 4, 220; 8, 160; 12, 130; 16, 105), the cost per
unit of physical volume and time,
$\propto h^{-4} \, \sum_n |c_n| \, / \, (C \cdot \text{speed})$, is 3.4 times higher for
order 4, and its memory, $\propto h^{-3}$, 5.5 times higher.

**Decision.** Eighth order at $C = 0.3$ is the default; order 4 keeps $C = 0.85$.

**At the default, the verification tests** improve about fivefold: explosion 0.11–0.26%,
force 0.07–0.18%. Two get worse:

- The free surface, 1.9% → 2.4%. Stress imaging is second order [@kristek2002] and is now the
  largest error in that test.
- Layer reflection, $1.5 \times 10^{-4}$ → $3.9 \times 10^{-4}$ ($-68$ dB) without a free
  surface, and $1.8 \times 10^{-3}$ → $2.1 \times 10^{-3}$ with one, as the longer stencils
  reach further across the layers' inner edge.

## 2026-11-20 — The adjoint wavefield by running the solver backwards

**Question.** Can the forward solver, stepped in reversed order, be the exact transpose of the
discrete forward run (the discretise-then-optimise adjoint), rather than only a discretisation
of the continuous adjoint equations [@tarantola1984; @tromp2005; @fichtner2006]? And where does
it stop being one?

**Theory.** One forward step $n$ is

$$
v^{n+\frac12} = v^{n-\frac12} + B D \sigma^{n} + B f^{n}, \qquad
\sigma^{n+1} = \sigma^{n} + M D' v^{n+\frac12} + g^{n},
$$

with $B = \Delta t/\rho$, $M = \Delta t\, C$ and $D' = -W_\sigma^{-1} D^{T} W_v$, where $W$ holds
each value's cell volume (half a cell on the free surface): the staggered derivatives are
negative adjoints, which is also what makes the scheme conserve energy (2026-11-20, the
self-adjoint free surface). Scaling the rows by $W B^{-1}$ and $W M^{-1}$ makes the stepping
symmetric except for the time differences, which transposing turns round. The adjoint
variables obey

$$
\lambda_v^{n+\frac12} = \lambda_v^{n+\frac32} - B D \lambda_\sigma^{n+1} + B W_v^{-1} r_v^{n},
\qquad
\lambda_\sigma^{n} = \lambda_\sigma^{n+1} - M D' \lambda_v^{n+\frac12} + M W_\sigma^{-1} r_\sigma^{n},
$$

with $r = R^{T} \delta$ the residuals $\delta$ spread by the transposed receivers. With
$\tilde v = \lambda_v$ and $\tilde \sigma = -\lambda_\sigma$ this is the forward step taken
backwards: adjoint step $m = N - 1 - n$ reads $\tilde\sigma^{n+1}$ and writes
$\tilde v^{n+\frac12}$, then $\tilde\sigma^{n}$. So

- a velocity residual recorded after forward step $n$ enters during adjoint step $m$ as a
  unit force, $B W_v^{-1} R^{T}$: exactly a force source's operator;
- a pressure residual, which the forward run read from $\sigma^{n+1}$, enters one adjoint step
  earlier as an isotropic stress source $-M W_\sigma^{-1} R^{T} = \Delta t\, K / V$ per unit
  ($K = \lambda + \tfrac23 \mu$; on the free surface the plane-stress modulus loads
  $\sigma_{xx}$ and $\sigma_{yy}$ over half a cell, and $\sigma_{zz}$ not at all). A pressure
  sample from the last forward step is the adjoint run's initial state;
- after adjoint step $m$ the adjoint velocities pair with the forward velocities at
  $(n + \tfrac12)\Delta t$ and the adjoint stresses with the forward stresses at $n\,\Delta t$,
  which is what the gradients will correlate.

Receivers recording every $k$-th step transpose to injecting at those steps only.

**Test** (`tests/physics.rs`). A dot-product test, $\langle F w, \tilde d \rangle = \langle w,
F^{T} \tilde d \rangle$, in the reciprocity test's model (two layers, a sphere, a free surface):
a force history $w$ at 200 m depth, recorded in $v_x$, $v_y$, $v_z$ and $p$ at four receivers
(one on the free surface), 0.5 s (about 1300 steps). $F^{T}\tilde d$ is the adjoint velocity
read at the source, along the force, after each adjoint step.

Two lessons about the test itself. With residuals of random sign, the four receivers' pressure
terms summed to a fiftieth of their size, so round-off in each (agreeing to $10^{-7}$ to
$4 \times 10^{-6}$ per receiver) showed as $10^{-4}$ in the total. And pressure, in Pa,
outweighs velocity by about $10^{6}$, so velocity hardly counted. The residuals are therefore
the data times random positive weights, divided by each component's $\|d_c\|^2$: every term
is positive and every component counts equally.

**Results.** $|\langle F w, \tilde d\rangle - \langle w, F^T \tilde d\rangle| /
\langle F w, \tilde d\rangle$:

| Walls | Signals | Order 4 | Order 8 |
|---|---|---|---|
| Rigid, free surface | white noise, every step | $5.9 \times 10^{-8}$ | $4.2 \times 10^{-8}$ |
| CPML and sponge, free surface | band-limited (12 Hz), every 2nd step | $4.8 \times 10^{-8}$ | $6.0 \times 10^{-8}$ |
| CPML all round | white noise, every step | $3.8 \times 10^{-8}$ | $4.3 \times 10^{-8}$ |

That is single-precision round-off: the time-reversed solver is the transpose of the forward
run, in every component, at the free surface and with decimated recording. The CPML and the
sponge are not self-adjoint in themselves, yet seen from outside the layers they are to this
precision. Inside a layer they are not. With the receivers moved into the top layer (no free
surface, order 4): $1.4 \times 10^{-3}$ at its outer edge, $4.1 \times 10^{-4}$ at 40 m,
$2.3 \times 10^{-4}$ at 80 m, and $9.9 \times 10^{-7}$ at 100 m, the first node outside it
($1.2 \times 10^{-6}$ at order 8; the stencil still reaches in).

**Decision.** The adjoint is the forward solver run backwards (`src/adjoint.rs`), with no
separate adjoint kernels to keep in step with the forward ones. Receivers belong outside the
layers (`check` already warns).

**Reproduce:** `the_adjoint_run_is_the_transpose_of_the_forward_run`, in `just physics`.

## 2026-11-20 — Gradients with respect to ρ, λ and μ

**Question.** The gradient of a misfit with respect to the model, exact for the discrete
forward run rather than for the wave equation, so that what it gets wrong can later be put
down to the physics and the reconstruction of the forward field, not to the adjoint.

**Theory.** Write the forward step as two equations,
$v^{n+\frac12} - v^{n-\frac12} = B\,(D\sigma^{n} + f^{n})$ and
$\sigma^{n+1} - \sigma^{n} = M e^{n+\frac12} + g^{n}$, with $e = D'v$ the strain rate the stress
update applies the stiffness to. Their Lagrange multipliers follow from the weighted adjoint of
the previous entry: $\hat\lambda_v = W B^{-1} \tilde v$ and $\hat\lambda_\sigma = -W \tilde\varepsilon$,
where the adjoint strain

$$
\tilde\varepsilon^{n} = \tilde\varepsilon^{n+1} + D'\tilde v^{n+\frac12} - W^{-1} r_\sigma^{n},
\qquad \tilde\varepsilon^{N} = -W^{-1} r_\sigma^{N},
$$

sums the adjoint strain rates backwards, with the pressure residuals ($-W^{-1} r_\sigma = \delta /
(3W)$ on the normal components); the adjoint stress is $\tilde\sigma = M\tilde\varepsilon$. Since
the forward equations hold, the derivative of the misfit is that of the Lagrangian with respect
to the material alone:

$$
\frac{\partial\chi}{\partial B} = \sum_n \frac{W\, \tilde v^{n+\frac12}\,(v^{n+\frac12} - v^{n-\frac12})}{B^{2}},
\qquad
\frac{\partial\chi}{\partial M} = -\sum_n W\, \tilde\varepsilon^{n+1} \otimes e^{n+\frac12},
$$

the zero-lag correlations of the density and stiffness kernels [@tromp2005; @fichtner2006;
@fichtner2011]. Forces drop out of the first because they enter as $B f$; moment tensors don't
involve the material. The chain rule then goes through the material averaging: $B = 2\Delta t /
(\rho_1 + \rho_2)$, so each node gets $-W \sum_n \tilde v \,\Delta v / (2\Delta t)$; $l2m =
\Delta t(\lambda + 2\mu)$ and $lam = \Delta t\,\lambda$ at the nodes; the shear moduli are
$\Delta t$ times the harmonic mean $H$ of four nodes' $\mu$, $\partial H / \partial \mu_q = H^2 /
(4\mu_q^2)$; on the free surface $\sigma_{xx}$ and $\sigma_{yy}$ use the plane-stress $a = l2m -
lam^2/l2m$ and $b = lam - lam^2/l2m$. The velocity parametrisation follows from
$\lambda = \rho(v_p^2 - 2v_s^2)$, $\mu = \rho v_s^2$.

The adjoint strain is summed rather than recovered from the adjoint stress as
$M^{-1}\tilde\sigma$: $M$ is singular in fluids, and its inverse cancels digits where
$\mu \ll \lambda$.

**Forward field.** The correlations need the forward velocities backwards in time. Here they
are replayed from checkpoints of the solver's state (fields and CPML memory) every $K$ steps,
$K = \lceil\sqrt{3N}\rceil$ by default, which about minimises the checkpoints plus one
segment's velocities. The replay is bit for bit (the updates don't depend on how rayon
splits the work), so it adds no error: a baseline for the boundary-value reconstruction to be
measured against. It costs one more forward run.

**Test** (`tests/physics.rs`). A free surface, two layers and a sphere; a force and an
explosion; ten receivers (eight on the surface) recording $v_x$, $v_y$, $v_z$ and $p$ every
second step for 0.35 s; the misfit $\tfrac12 \sum_c \|d_c - d_c^{obs}\|^2 / \|d_c^{obs}\|^2$
against data from a different sphere. Each of $v_p$, $v_s$ and $\rho$ is perturbed as
$m(1 + h u)$ with $u$ uniform in $[-1, 1]$ at every node outside the layers ($v_p$ only below
the fastest, which sets $\Delta t$ and the layers' damping), and the central difference of the
misfit compared with the gradient's prediction.

**Results.** The plain central differences converge on the gradient exactly as $h^2$: for
$v_p$ at order 4 they are off by 12.8%, 3.2%, 0.79% and 0.18% at $h = 1.6\times10^{-2}$,
$8\times10^{-3}$, $4\times10^{-3}$ and $2\times10^{-3}$, then single-precision round-off takes
over (0.1% at $10^{-3}$ and below). Richardson's extrapolation from $h = 1.6\times10^{-2}$ and
$8 \times 10^{-3}$ removes the $h^2$ term:

| Perturbed | Misfit | Order 4 | Order 8 |
|---|---|---|---|
| every node | all components | $v_p$ $1.7\times10^{-4}$, $v_s$ $1.6\times10^{-6}$, $\rho$ $1.5\times10^{-5}$ | $6.5\times10^{-5}$, $5.8\times10^{-6}$, $1.2\times10^{-5}$ |
| the free surface only | all components | $1.5\times10^{-5}$, $1.2\times10^{-4}$, $2.9\times10^{-6}$ | $2.8\times10^{-4}$, $7.3\times10^{-7}$, $1.0\times10^{-5}$ |
| every node | pressure only | $4.2\times10^{-5}$, $3.2\times10^{-4}$, $1.4\times10^{-5}$ | $1.1\times10^{-4}$, $9.9\times10^{-5}$, $3.4\times10^{-6}$ |

What is left is the extrapolation's own error ($h^4$ terms where the misfit is most nonlinear:
$v_s$ with pressure alone, whose plain differences were 95% off) and round-off. Nothing
suggests an error in the gradient at the $10^{-4}$ level, at the free surface, through the
pressure residuals, or through the force's buoyancy.

**Cost.** On `examples/explosion.toml` (101 × 101 × 67 nodes, 1510 steps, a sphere in the
observed model) a gradient took 5.1 times a forward run: the forward run, its replay, the
adjoint run and two strain-rate passes per step. The checkpoints (every 68 steps) and one
segment's velocities held 1.66 GiB against the forward run's 0.08 GiB: about
$2\sqrt{27 N}$ padded arrays, growing as $\sqrt{N}$. `examples/layered.toml` (4456 steps)
would need about 6 GiB, which is what reconstructing the forward field instead would save.

**Decision.** `gradient::gradient` computes $\partial\chi/\partial\rho$, $\partial\chi/\partial
\lambda$ and $\partial\chi/\partial\mu$ at the nodes (and $v_p$, $v_s$, $\rho$ through
`Gradient::velocities`) with checkpointed replay. Inside the absorbing layers the strain rates
leave out the CPML and the layers aren't self-adjoint, so the gradient there is approximate;
models are usually held fixed there anyway.

**Reproduce:** `the_gradient_matches_finite_differences`, in `just physics`. The surface-only and
pressure-only rows came from the same test with the perturbation restricted to $k = 0$ and the
velocity components' weights set to zero. The cost: `just run` of the example with a sphere
added and `dt = 5.3e-4`, then `just gradient` of the example itself against it.

## Open questions

- **Free surface for long stencils.** Replace stress imaging with a formulation of matching
  order [@mittet2002; @kristek2002] and measure on Lamb's problem.
- **Interfaces.** Band-limit the model to the grid [@mittet2017] and measure reflection and
  transmission at a planar interface against analytic coefficients.
- **Source and receiver positions** off the nodes [@hicks2002].
- **Absorbing layers at high order.** Tune the profile or width so reflection returns to
  $-76$ dB.
- **Error propagation into adjoint fields and gradients.** Forward-field error (dispersion,
  interfaces, free surface, layers), reconstruction error (time-reversed propagation in single
  precision; the non-reversible layers [@dussaud2008; @aaker2020]), and their effect on
  gradients and inversions [@tromp2005; @fichtner2006; @fichtner2011]. Settled: adjoint
  consistency (the dot-product entry) and the gradient of the discrete run against finite
  differences (the gradient entry), with an exactly replayed forward field. Still to
  measure: the reconstruction error of a boundary-value reconstruction against that replay,
  and the accuracy study repeated for gradients (against the gradient of the exact wave
  equation, e.g. a fine-grid reference) as a function of stencil order, $C$ and ppw.
