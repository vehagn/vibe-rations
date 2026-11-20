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

## 2026-11-20 — A research brief, checked

**What.** Another agent reviewed the literature for this simulator and ran its own von Neumann
analysis, ending in nine recommendations. Its references came from memory, so each one cited
here was checked on Crossref first. Its numbers reproduce exactly with
`docs/analysis/time_dispersion.py`, which reuses this repository's coefficients and dispersion
relation: ppw for 0.5% S group error of 4.8 at order 8, $C = 0.3$; 5.0, 4.0 and 3.6 at orders
8, 12 and 16 at $C = 0.85$ with the time error removed; 3.4 for eight-point least-squares
coefficients without time error.

**Verdicts.**

| Recommendation | Verdict |
|---|---|
| Remove leapfrog's time dispersion by transforms [@wang2015; @koene2018] | Confirmed on the 3D solver and adopted, opt-in (next entry). The largest gain. |
| Receiver interpolation as long as the stencil | Adopted (the entry after next). |
| Least-squares coefficients, once the time error is gone [@holberg1987] | Promising: about 0.1 of today's work by the analysis. Adopted as minimax coefficients (2026-11-20, optimised coefficients). |
| Sources and receivers between nodes [@hicks2002; @petersson2016] | Needed for FWI with real acquisition geometry. Open. |
| Volume-averaged material at interfaces | Needs the planar-interface test first. Open. |
| A free surface of the stencil's order, by summation-by-parts closures [@oreilly2017] | Taken to be the largest error in its test (2.4%); most of that turned out to be the test's reference (three entries on). Open; the most work. |
| Absorbing layers: tune at order 8; their transpose | The adjoint run's dot-product test (2026-11-20, the adjoint wavefield) measured the transpose: exact to round-off outside the layers, $1.4 \times 10^{-3}$ inside. Tuning is open. |
| Flush subnormals to zero | Can't be measured here (Apple silicon handles subnormals at full speed). Open for x86. |
| Verification: grid refinement, gradient Taylor tests | The dot-product test is done (the adjoint wavefield); Taylor tests belong with the gradients. |

The brief's caution that the transform pair might need amplitude factors from the paper turned
out unnecessary for this discretisation: the derivation in the next entry has none, and the 1D test agrees
to $10^{-8}$.

## 2026-11-20 — Removing time dispersion

**Question.** The previous entry kept $C = 0.3$ because leapfrog's time error, not the
stencil, limits orders $\ge 8$ at larger steps. Can that error be removed instead of made
small, as [@wang2015; @koene2018; @mittet2019] propose? The research brief (previous entry)
predicted it would cut today's work by 2–2.4× at order 8 and 4.5–6.7× at orders 12–16.

**Theory.** Every time derivative in a step is a centred difference over $\Delta t$. For a
signal $u^n = \hat u\, e^{i \omega n \Delta t}$ (at whatever half-step offset),

$$
\frac{u^{n+\frac12} - u^{n-\frac12}}{\Delta t} = i \omega' u^n,
\qquad \omega' = \frac{2}{\Delta t} \sin\frac{\omega \Delta t}{2}.
$$

Nothing else in a step depends on time: the medium, the stencil, the free surface's mirror
images and the operators by which sources enter are the same every step. So the stepped
wavefield at $\omega$ is the time-continuous (semi-discrete) one at $\omega'$, driven by the
injected samples' spectrum at $\omega$: $\hat u_{\Delta t}(\omega) = G(\omega')\, \hat s_{\Delta t}(\omega)$,
with $G$ the semi-discrete transfer function from source to receiver. Two transforms then
recover $G(\omega') S(\omega')$, the semi-discrete trace:

- before the run, inject samples with spectrum $\hat s_{\Delta t}(\omega) = S(\omega'(\omega))$;
- after it, read the trace's spectrum at $\omega'$ from the recorded one at
  $\omega = \frac{2}{\Delta t} \arcsin\frac{\omega' \Delta t}{2}$.

Two details matter. Phases are about $t = 0$, when the run starts at rest, so each signal
carries its first sample's time: forces enter at $n \Delta t$, moment tensors at
$(n + \frac12)\Delta t$, velocities are read at $(n + \frac12)\Delta t$ and stresses at
$(n + 1)\Delta t$. And the transform is exact wherever the update is a centred difference; the
CPML memory variables and the sponge are recursive filters of their own, so in the layers it
is approximate, which only affects waves being absorbed.

**Implementation** (`src/tdt.rs`, `[scheme] remove_time_dispersion`). Each transform is a real
linear map of a sampled signal: the spectrum at the warped frequencies by direct summation,
then the inverse DFT on a grid twice the signal's length, band-limited to twice the sources'
highest frequency (5 $f_0$, where a Ricker spectrum is $10^{-9}$ of its peak). It costs
$O(N \cdot N_\omega)$ per trace with $N_\omega \approx 10 f_0 T$ a few hundred, which doesn't
show in the run time. The transpose is provided, for adjoint sources.

**Tests.** A 1D staggered velocity–stress run timed like the solver, at $C = 0.9$, against
the same run at $C = 0.1$ (both transformed): 24% of the peak apart without the transform,
$1.1 \times 10^{-8}$ with it. Warping then un-warping a Ricker wavelet gives it back to
$8 \times 10^{-9}$; the transpose passes a dot-product test to $10^{-10}$.

**Measurement** (`just accuracy`, `removing_time_dispersion`, with the stencil-matched
receivers of the entry after this one). L2 misfit at 1800 m, about ten wavelengths, against
the closed-form solutions; "+ T" is $C = 0.85$ with the transform:

| Order | P 4.0 ppw: $C = 0.3$ | + T | S 3.5: $C = 0.3$ | + T | S 3.9: $C = 0.3$ | + T | S 4.6: $C = 0.3$ | + T |
|---|---|---|---|---|---|---|---|---|
| 4 | 33% | 38% | 62% | 64% | 44% | 46% | 26% | 27% |
| 8 | 4.8% | 4.0% | 7.1% | 8.8% | 2.9% | 3.9% | 1.26% | 1.21% |
| 12 | 5.5% | 2.9% | 2.5% | 2.0% | 1.9% | 0.63% | 1.46% | 0.17% |
| 16 | 5.5% | 2.9% | 2.4% | 0.69% | 1.9% | 0.19% | 1.43% | 0.12% |

- Order 8 with the transform matches order 8 at $C = 0.3$ at 4.6 ppw, in 2.8 times fewer
  steps. On the coarsest grids it is up to 1.3 times worse, as the von Neumann analysis
  predicts (5.0 against 4.8 ppw for 0.5%): at $C = 0.3$ the time error, which speeds waves
  up, still cancels part of the spatial error, which slows them down.
- At $C = 0.3$, orders 12 and 16 were limited by the time error. Without it they are 8–12
  times more accurate than the default at 4.6 ppw. Per unit volume and time, with the
  measured speeds, their steps cost 0.45 and 0.57 of the default's, so for the default's
  accuracy at 4.6 ppw (1.26%), order 12 at 3.9 ppw gives 0.63% for 0.23 of its compute and
  0.61 of its memory, and order 16 at 3.5 ppw 0.69% for 0.19 and 0.44. The brief predicted
  0.22 and 0.18.
- $C = 0.95$ gives the same misfits to within 0.06 points, for 11% fewer steps.
- Order 4 gets worse, as predicted: its time and space errors cancelled.
- Limiting the plain traces to the transform's band changes nothing (the same misfits to
  0.07 points), so the gain is the warp's, not the filter's. The transform's own cost doesn't
  show in the measured speeds.
- With the cubic receivers, the transformed P misfits stopped at 3.6% for orders 12 and 16;
  the cubic was what limited them near the source (next entry).

A first try at a 4.0 ppw grid ($h = 17.32$ m) gave 26–33% at 300 and 1500 m for every order
and every time step: the spacing didn't divide the distances, and receivers snapped to nodes
up to 6.9 m (0.4 h) from where they belonged. That is the brief's point about sources and
receivers between nodes, measured by accident; the study now uses $h = 300/17$ m.

**Physics test** (`removing_time_dispersion_allows_a_large_time_step`). The point force
against Stokes' solution at $C = 0.85$ with the transform: 0.001–0.008%, against 0.03–0.09%
at $C = 0.3$ without it. On that grid (6.9 ppw) the time error dominated even at 0.3.

**Decision.** The transform is in, behind `[scheme] remove_time_dispersion = true`, which
also makes 0.85 the default Courant fraction. Not yet the default: the adjoint run must first
apply the transposed transform to the residuals (`Simulation::trace_warp`), and snapshots
stay warped. With it, the best settings for waveform fits become order 12–16 at 3.5–4 ppw.

**Reproduce:** `cargo test --release --lib tdt`; `removes_the_time_error_of_a_staggered_run`;
`just accuracy`; `python3 docs/analysis/time_dispersion.py` for the predictions.

## 2026-11-20 — Receiver interpolation as long as the stencil

**Question.** Reading a staggered value at a node with the cubic $(-1, 9, 9, -1)/16$ keeps
$R(kh) = \sum_n 2 w_n \cos\big((n - \tfrac12) kh\big)$ of a wave's amplitude: it loses 7.6%
at 4.5 ppw and 18.6% at 3.5 (`docs/analysis/time_dispersion.py`), and a source, its
transpose, loses as much again. The brief suggested that this, not dispersion, held back the
longer stencils' P misfits.

**Change.** $2L$-point Lagrange weights, $L$ the stencil's half-length: order 4 keeps the
cubic, order 8 uses $(-5, 49, -245, 1225, 1225, -245, 49, -5)/2048$ (1.0% lost at 4.5 ppw),
and fewer points near the edges of the grid. Folding about the free surface is unchanged, and
sources stay the transpose of receivers.

**Results.** Physics tests at order 8: explosion 0.11–0.26% → 0.07–0.19%, point force
0.07–0.18% → 0.03–0.09%; the free surface (2.40%) and the absorbing layers are unchanged;
reciprocity $1.5 \times 10^{-6}$ → $1.9 \times 10^{-6}$, more terms of single-precision
round-off. In the accuracy study, P at 4.0 ppw, where $v_x$ is interpolated along the path:

| | 300 m, cubic | 300 m, matched | 1800 m, cubic | 1800 m, matched |
|---|---|---|---|---|
| Order 8, $C = 0.3$ | 2.16% | 0.76% | 5.20% | 4.75% |
| Order 16, $C = 0.3$ | 2.22% | 0.81% | 5.84% | 5.47% |
| Order 12, $C = 0.85$ + transform | 1.86% | 0.07% | 3.59% | 2.91% |
| Order 16, $C = 0.85$ + transform | 1.87% | 0.03% | 3.59% | 2.91% |

The S misfits did not change: the S test reads $v_y$ across the path, where the field hardly
varies along $y$. Near the source the cubic was the largest error; far away, dispersion is.

**Decision.** Adopted. The order-4 baselines are untouched (the cubic is its $2L$-point
interpolation).

**Reproduce:** `just physics`; `just accuracy`; `interpolation_matches_the_stencil`.

## 2026-11-20 — Time dispersion removed by default, and the free-surface test's floor

**Decision.** The time-dispersion transform is the default (`remove_time_dispersion = true`,
$C = 0.85$), with `false` to opt out. The adjoint run needs the transposed transform on
residuals of finished traces; the dot-product test of the adjoint run steps the solver
directly and is unaffected.

**The paper.** Checked against [@koene2018]'s text (sections 2–5): our pair is its
eqs. 13–17, the source warped with $\sin$ and the traces with $\arcsin$, no amplitude factors,
zero-padded to twice the length. It doesn't discuss the time origin or staggered time stamps,
absorbing layers, or adjoints; it criticises [@wang2015] for using the inverse transform's
adjoint as the forward one. Here that adjoint serves only as the adjoint source, where it is
the exact transpose of the data operator.

**Physics tests**, order 8, before (time error kept, $C = 0.3$) and after (removed,
$C = 0.85$):

| Test | Before | After |
|---|---|---|
| Explosion against the analytic solution | 0.07–0.19% | 0.002–0.006% |
| Point force against Stokes' solution | 0.03–0.09% | 0.001–0.008% |
| Free-surface P doubling | 2.40% | 2.53% |
| Layer reflection, without / with a free surface | $3.9 \times 10^{-4}$ / $2.1 \times 10^{-3}$ | $1.7 \times 10^{-4}$ / $1.9 \times 10^{-3}$ |
| Reciprocity (orders 4 and 8) | $\le 1.9 \times 10^{-6}$ | $\le 6.9 \times 10^{-7}$ |
| 9 s soft-sediment run, left of its peak | $\le 3.0 \times 10^{-6}$ | $\le 3.0 \times 10^{-6}$ |

Everything improves except the free surface, and the 9 s run with CPML, the sponge and a
free surface stays stable at $C = 0.85$ for order 8.

**The free-surface test's floor.** Why did it get worse? Run on finer grids
(`free_surface_on_finer_grids` in `just accuracy`), its misfit hardly moves:

| $h$ | Order 8, removed | Order 8, kept ($C = 0.3$) | Order 16, removed |
|---|---|---|---|
| 10 m | 2.53% | 2.40% | 2.52% |
| 7.5 m | 2.37% | 2.29% | 2.36% |
| 5 m | 2.26% | 2.23% | 2.26% |

A fit of $R + c\,h^2$ to either column gives $R = 2.17\%$, the same for both, and
predicts the 7.5 m row to 0.01 points. That floor is the reference's: it doubles the incident
P wave as if it were plane, and a spherical wave from 300 m below, 1.5 dominant wavelengths,
is not. The scheme's own
free-surface error is the rest, about 0.36% at $h = 10$ m (8 ppw for P) with the time error
removed and 0.23% with it kept, where the leapfrog error partly cancels it. The test's
tolerance moves from 2.5% to 2.6%, with both numbers in its doc comment. Earlier entries
called stress imaging the largest error in its test; it isn't, and a higher-order free
surface is less urgent than thought. Measuring it needs an exact reference: Lamb's problem.

**Reproduce:** `just physics`; `free_surface_on_finer_grids` in `just accuracy`.

## 2026-11-20 — Optimised coefficients

**Question.** Least-squares coefficients gained only 0.2–0.8 ppw while the time error
remained (2026-11-20, stencil order and time step). With the time error removed, the brief's
analysis predicted eight-point least-squares coefficients would need 3.4 ppw where Taylor's
need 5.0. Do optimised coefficients pay off in 3D runs?

**Criterion.** With the time error removed the trace's error is the spatial one, and for a
plane wave it is a phase error: the operator gives
$K(\xi) = 2 \sum_n c_n \sin\big((2n-1)\xi/2\big)$ for $\xi = kh$ along each axis, and over a
distance $D$ a wave's phase is off by $(\lvert\mathbf K\rvert - \lvert\mathbf k\rvert h)\,D/h$.
So the design minimises the largest phase error per node,
$\max \big\lvert \lvert\mathbf K\rvert - \lvert\mathbf k\rvert h \big\rvert$, over every direction
(91 sampled from axis to diagonal; the error is symmetric under permutations and sign
changes) and every $\lvert\mathbf k\rvert h$ up to the band $\beta = 2\pi/\text{ppw}$ of the run's
shortest wavelength, linearised as $\sum_i d_i\, \delta(\xi d_i)$ with
$\delta(x) = K(x) - x$, which makes it a linear minimax problem. Consistency,
$\sum_n (2n-1)\, c_n = 1$, is kept, so linear functions stay exact. Lawson's algorithm
[@rice1968] solves it; [@holberg1987] bounded the group-velocity error the same way, and
[@zhang2013; @liu2014] optimise staggered operators in related ways. The receivers'
interpolation weights get the same treatment for their amplitude error.

**First try: plain minimax over the run's band.** S misfits fell five- to tenfold, but P
became far worse: on the P test's grids S has only 2.3 ppw, the band was nearly $\pi$, and
the evenly spread error spoiled every wavenumber (order 8, P at 4.0 ppw: 78% against Taylor's
4.0%).

**Two fixes, measured** (`optimised_coefficients` in `just accuracy`, misfit at 1800 m;
Taylor in brackets):

- *Cap the error and shrink the band.* With at most $10^{-3}$ per node, P recovered (order
  8, 4.0 ppw: 4.6%) but stayed worse near the source: at 4.0 ppw, 0.85–1.9% at 300 m against
  Taylor's 0.03–0.46%. A uniform $10^{-3}$ per node is already 1.5% after 15 nodes. $3 \times 10^{-4}$ and
  $10^{-4}$ helped further; $10^{-4}$ was chosen.
- *Weight towards long waves.* An equiripple error reaches the dominant frequency and the P
  waves, where Taylor's error is nearly zero. Weighting the error by $(\xi/\beta)^{-p}$, with
  the cap at the band's edge:

| | Taylor | $p = 0$ | $p = 1$ | $p = 2$ |
|---|---|---|---|---|
| S 3.0 ppw, order 8 / 12 / 16 | 19.5 / 6.5 / 2.9% | 5.4 / 0.72 / 0.38% | 5.8 / 0.75 / 0.39% | 7.3 / 0.89 / 0.42% |
| S 3.5 ppw, order 8 / 12 / 16 | 8.8 / 2.0 / 0.69% | 1.4 / 0.25 / 0.14% | 1.6 / 0.26 / 0.15% | 2.1 / 0.28 / 0.15% |
| P 4.0 ppw, order 8 / 12 / 16 | 4.0 / 2.9 / 2.9% | 3.0 / 3.1 / 3.1% | 3.0 / 3.1 / 3.0% | 3.0 / 3.0 / 2.9% |
| P at 300 m, 4.0 ppw, order 8 / 12 / 16 | 0.46 / 0.07 / 0.03% | 0.11 / 0.14 / 0.20% | 0.09 / 0.08 / 0.09% | 0.08 / 0.04 / 0.04% |
| Explosion physics test (fine grid) | 0.002–0.006% | 0.016–0.021% | 0.009–0.013% | 0.005–0.007% |

$p = 2$ gives up a little at the band's edge (S at 3.0 ppw, order 8) for long waves that are
as good as Taylor's; nothing measured is worse than Taylor's beyond noise. The common floor
of about 3% for P at 1800 m, shared by Taylor's orders 12 and 16, is not dispersion; probably
the narrow domain's side layers at grazing incidence.

**The design** at $p = 2$ (`optimised_coefficients_beat_taylor_over_their_band`): largest
phase error per node over all directions, Taylor → optimised, order 8 at $\beta = 1.8$ (3.5
ppw): $1.3 \times 10^{-2}$ → $4.8 \times 10^{-4}$; order 12 at 2.0: $8.9 \times 10^{-3}$ →
$6.4 \times 10^{-5}$; order 16 at 2.4: $3.0 \times 10^{-2}$ → $1.3 \times 10^{-4}$.
Interpolation errors fall 10–60-fold. $\sum_n \lvert c_n\rvert$ rises 4–6%, and the time step
falls as much.

**Results** at $p = 2$, $10^{-4}$, misfit at 1800 m (Taylor → optimised):

| Waves, ppw | Order 8 | Order 12 | Order 16 |
|---|---|---|---|
| S, 3.0 | 19% → 7.3% | 6.5% → 0.89% | 2.9% → 0.42% |
| S, 3.5 | 8.8% → 2.1% | 2.0% → 0.28% | 0.69% → 0.15% |
| S, 3.9 | 3.9% → 0.61% | 0.63% → 0.13% | 0.19% → 0.11% |
| S, 4.6 | 1.2% → 0.19% | 0.17% → 0.12% | 0.12% → 0.12% |
| P, 2.9 | 19% → 7.8% | 6.9% → 1.7% | 3.5% → 1.4% |
| P, 3.5 | 7.8% → 2.8% | 2.8% → 2.2% | 2.2% → 2.2% |
| P, 4.0 | 4.0% → 3.0% | 2.9% → 3.0% | 2.9% → 2.9% |

Order 12 at 3.0 ppw (0.89%) beats order 8 with Taylor coefficients at 4.6 (1.2%), for about
a quarter of the compute and 0.28 of the memory. Physics tests, Taylor → optimised: explosion
0.002–0.006% → 0.005–0.007% (the one that got worse, on a grid where Taylor is nearly exact
for these waves), point force 0.001–0.008% → 0.001–0.005%, free surface 2.53% → 2.52%, layers
unchanged, reciprocity $\le 6.9$ → $7.9 \times 10^{-7}$; the 9 s run stays stable (any
coefficients keep the scheme self-adjoint).

**Decision.** Optimised coefficients are the default (`coefficients = "optimised"`,
`"taylor"` to opt out), with $p = 2$ and a cap of $10^{-4}$ per node. The band follows the
model's slowest wave, so an inversion must pin it (`optimised_for`) or the operator would
change between iterations.

**Reproduce:** `cargo test --release --lib fd::`; `optimised_coefficients` in
`just accuracy`; `just physics`.

## Open questions

- **Free surface for long stencils.** Replace stress imaging with a formulation of matching
  order [@mittet2002; @kristek2002], e.g. staggered summation-by-parts closures
  [@oreilly2017], which keep the exact discrete adjoint. Its error is about 0.4% at 8 ppw
  for P (the free-surface test's floor), so measure it first on Lamb's problem with an
  exact solution.
- **Interfaces.** Band-limit the model to the grid [@mittet2017] and measure reflection and
  transmission at a planar interface against analytic coefficients.
- **Source and receiver positions** off the nodes [@hicks2002; @petersson2016]: snapping
  a receiver 0.4 h off its position cost 26% misfit at 3.9 ppw (2026-11-20, removing time
  dispersion).
- **Optimised coefficients, further.** Weight the error by the wavelets' spectrum and the
  model's range of speeds instead of $(\xi/\beta)^{-2}$; measure the observed order under grid
  refinement (the band shrinks with $h$, so the operator tends to Taylor's); and check the
  gradients' accuracy with them, since FWI pins the band.
- **The adjoint with the transform.** Residuals of finished traces enter through the
  transposed transform. Since the transformed traces equal the time-continuous ones for every
  model, their gradient should too, taking the time error out of the gradient as well:
  measure it against finite differences. And whether order 12 or 16 should be the default.
- **Absorbing layers at high order.** Tune the profile or width so reflection returns to
  $-76$ dB.
- **Error propagation into adjoint fields and gradients.** Forward-field error (dispersion,
  interfaces, free surface, layers), reconstruction error (time-reversed propagation in single
  precision; the non-reversible layers [@dussaud2008; @aaker2020]), and their effect on
  gradients and inversions [@tromp2005; @fichtner2006; @fichtner2011]. Planned measures:
  dot-product and reciprocity tests of adjoint consistency, gradients against finite
  differences, and the accuracy study repeated for gradients as a function of stencil order,
  $C$ and ppw.
