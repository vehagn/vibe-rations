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
a force history $w$ at 200 m depth, recorded in $v_x$, $v_y$, $v_z$ and $p$ at five receivers
(one on the free surface, one a node below it, where the stresses' mirror images reach),
0.5 s (about 1300 steps). $F^{T}\tilde d$ is the adjoint velocity
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
| Rigid, free surface | white noise, every step | $4.6 \times 10^{-8}$ | $7.0 \times 10^{-8}$ |
| CPML and sponge, free surface | band-limited (12 Hz), every 2nd step | $2.3 \times 10^{-8}$ | $6.6 \times 10^{-8}$ |
| CPML all round | white noise, every step | $4.2 \times 10^{-8}$ | $8.3 \times 10^{-8}$ |

(With four receivers, none a node below the surface, the numbers were similar, but leaving
out the mirroring of the adjoint's initial pressure state went unnoticed; with the fifth it
gives $2.4 \times 10^{-7}$ to $2.8 \times 10^{-6}$, so the test now fails above
$1.5 \times 10^{-7}$.)

That is single-precision round-off: the time-reversed solver is the transpose of the forward
run, in every component, at the free surface and with decimated recording. The CPML and the
sponge are not self-adjoint in themselves, yet seen from outside the layers they are to this
precision. Inside a layer they are not. With the receivers moved into the top layer (no free
surface, order 4): $1.4 \times 10^{-3}$ at its outer edge, $4.1 \times 10^{-4}$ at 40 m,
$2.3 \times 10^{-4}$ at 80 m, and $9.9 \times 10^{-7}$ at 100 m, the first node outside it
($1.2 \times 10^{-6}$ at order 8), where $v_z$ is still interpolated from a value half a
node inside. The same at a side layer, where the sponge damps too: with one surface receiver
on the first node outside it, $2.8 \times 10^{-5}$ (order 4), and $3.9 \times 10^{-8}$ a node
further in.

**Decision.** The adjoint is the forward solver run backwards (`src/adjoint.rs`), with no
separate adjoint kernels to keep in step with the forward ones. Receivers belong outside the
layers (`check` warns), and velocity receivers a node clear of them:
`Simulation::receivers_next_to_layers` counts those that aren't, and `gradient` warns. A
plain forward run isn't warned: its traces hardly change, and the examples' end receivers sit
there.

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
(4\mu_q^2)$, once for every slot the node fills (on the grid's last planes the averaging
clamps the missing neighbours, so a node can fill two or four); on the free surface $\sigma_{xx}$ and $\sigma_{yy}$ use the plane-stress $a = l2m -
lam^2/l2m$ and $b = lam - lam^2/l2m$. The velocity parametrisation follows from
$\lambda = \rho(v_p^2 - 2v_s^2)$, $\mu = \rho v_s^2$.

The adjoint strain is summed rather than recovered from the adjoint stress as
$M^{-1}\tilde\sigma$: $M$ is singular in fluids, and its inverse cancels digits where
$\mu \ll \lambda$.

Where this is not the exact gradient, besides inside the layers:

- At a fluid node next to solid ones, $H$ has a kink in that node's $\mu$ (0 for $\mu \le 0$,
  $\approx 4\mu/k$ above, $k$ the slots it fills), so $\partial\chi/\partial\mu$ there is
  one-sided; `Gradient::mu` reports the part through the normal stresses. The $v_s$
  gradient is exact: $\mu = \rho v_s^2$ smooths the kink out.
- With a free surface the sponge damps the fields after the stress update, so the stored
  velocities the strain rates are taken from are damped where the update's were not, and it
  also reaches the values half a node outside the side layers. Within a stencil's reach of
  the layers' inner edge its factors differ from 1 by $10^{-7}$ to $2 \times 10^{-3}$, and
  enter through the outer stencil coefficients: about $10^{-6}$ in the gradient near the
  layers, below what finite differences resolve.
- Velocity receivers on the first node outside a layer read a damped value (previous entry).

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

A second model (`the_gradient_matches_finite_differences_with_walls_and_water`) covers what
the first leaves out: rigid walls, so that every node is perturbed, the last ones whose
neighbours the averaging clamps included; spacing $11 \times 12 \times 10.5$ m; a water layer
under the free surface; a double couple; every third step recorded. Extrapolated: $v_p$
$4.6 \times 10^{-5}$, $v_s$ $8.9 \times 10^{-7}$, $\rho$ $2.2 \times 10^{-4}$ (the plain
differences for $\rho$ were 44% off). Perturbing only the water layer, $\rho$ agrees to
$2.7 \times 10^{-5}$ and $v_p$ to $1.6 \times 10^{-4}$; only the outermost planes, the plain
difference at $h = 8 \times 10^{-3}$ to $1.1 \times 10^{-4}$–$4.6 \times 10^{-4}$, where
round-off already dominates.

What is left is the extrapolation's own error ($h^4$ terms where the misfit is most nonlinear:
$v_s$ with pressure alone, whose plain differences were 95% off) and round-off, which weighs
most where the predicted change is small. Nothing suggests an error in the gradient at the
$10^{-4}$ level, at the free surface, in fluids, at the grid's edges, through the pressure
residuals, or through the force's buoyancy. A review did find one: listing a component twice
(`p` twice) back-propagated both copies but put only the first into the adjoint strain, which
gave the $v_p$ direction the wrong sign; set-up now refuses duplicates.

**Cost.** On `examples/explosion.toml` (101 × 101 × 67 nodes, 1510 steps, a sphere in the
observed model) a gradient took 4.9 times a forward run (111 s against 23 s on 4 cores): the
forward run, its replay, the adjoint run and two strain-rate passes per step. It held at most
1.91 GiB against the forward run's 0.08 GiB: 1.66 GiB of checkpoints (every 68 steps) and one
segment's velocities, about $2\sqrt{27 N}$ padded arrays and growing as $\sqrt{N}$, and the
rest for the second solver, the strains and the sums. `examples/layered.toml` (4456 steps)
would need over 6 GiB, which is what reconstructing the forward field instead would save. (Both
examples' end receivers sit on the first node outside the layers, so `gradient` warns about
2 receivers there.)

**Decision.** `gradient::gradient` computes $\partial\chi/\partial\rho$, $\partial\chi/\partial
\lambda$ and $\partial\chi/\partial\mu$ at the nodes (and $v_p$, $v_s$, $\rho$ through
`Gradient::velocities`) with checkpointed replay. Inside the absorbing layers the strain rates
leave out the CPML and the layers aren't self-adjoint, so the gradient there is approximate;
models are usually held fixed there anyway.

**Reproduce:** `the_gradient_matches_finite_differences` and
`the_gradient_matches_finite_differences_with_walls_and_water`, in `just physics`. The
restricted rows and checks came from the same tests with the perturbation restricted ($k = 0$,
$k \le 2$, the outermost planes) or the velocity components' weights set to zero. The cost: `just run` of the example with a sphere
added and `dt = 5.3e-4`, then `just gradient` of the example itself against it.
## 2026-11-20 — The adjoint and gradients with the new defaults

**Question.** The adjoint run and the gradients (previous two entries) were built with Taylor
coefficients at $C = 0.3$ and cubic receivers, in parallel with the transform, the optimised
coefficients and the stencil-length receivers. Are they still exact with all of these?

**The transform.** With the time error removed, the run maps the model to $d = W\,u(m)$: the
samples $u$ the solver records, from pre-warped sources, un-warped by the linear map $W$
(`Warp`). A misfit of $d$ has the derivative $W^{T} r$ with respect to $u$, $r$ being its
derivative with respect to $d$. So `gradient` un-warps the traces after its forward run, as
`run` does (the misfit is then of the same traces as observed data written by `run`), and
`Simulation::recorded_residuals` applies $W^{T}$ before the adjoint run injects the
residuals. The gradient is that of the whole discrete map, transforms included.
[@koene2018] don't treat adjoints; [@wang2015]'s use of the inverse transform's adjoint as
the forward transform, which they criticise, is a different thing.

**Optimised coefficients** follow the model's slowest wave, which a finite-difference check
perturbs, so the operator would change with the model. The checks pin `optimised_for`, and
`Simulation::with_model` warns when it isn't.

**Receivers next to the layers.** At first the dot-product test failed at order 8 with CPML
along a free surface: $2.4 \times 10^{-7}$ against $1.5 \times 10^{-7}$; over six seeds
$2 \times 10^{-7}$ to $1.2 \times 10^{-6}$ (Taylor coefficients $4.7 \times 10^{-8}$ to
$5.9 \times 10^{-7}$), against $2$–$5 \times 10^{-8}$ at order 4. Turning the sponge off left
$1.5 \times 10^{-7}$. It looked like the layers, which aren't self-adjoint, leaking through
the longer stencil. It was one receiver, 120 m from a side (node 12, two nodes past the
10-node layer): the cubic read nodes 10–13, eight-point interpolation reads 8–15, two of
them damped. `Simulation::receivers_next_to_layers` named it once the test asserted there
were none; moved to 140 m, that case gives $4.7 \times 10^{-9}$. The gradient check's
receiver line, from 100 m with 8-node layers, had the same problem and now starts at 120 m.

**Results**, extrapolated finite differences against the gradient ($v_p$ / $v_s$ / $\rho$):

| Check | Before: Taylor, $C = 0.3$, cubic | After: defaults |
|---|---|---|
| Dot product, worst of six cases | $8.3 \times 10^{-8}$ | $7.0 \times 10^{-8}$ |
| Gradient, order 4 | $1.7 \times 10^{-4}$ / $1.6 \times 10^{-6}$ / $1.5 \times 10^{-5}$ | $2.8 \times 10^{-4}$ / $2.1 \times 10^{-5}$ / $4.4 \times 10^{-5}$ |
| Gradient, order 8 | $6.5 \times 10^{-5}$ / $5.8 \times 10^{-6}$ / $1.2 \times 10^{-5}$ | $3.9 \times 10^{-6}$ / $1.0 \times 10^{-5}$ / $1.6 \times 10^{-5}$ |
| Walls and water (Taylor coefficients) | $4.6 \times 10^{-5}$ / $8.9 \times 10^{-7}$ / $2.2 \times 10^{-4}$ | $6.4 \times 10^{-5}$ / $1.2 \times 10^{-5}$ / $1.6 \times 10^{-4}$ |

The receiver line moved, so the directional derivatives differ; order 4's $v_p$ one is small
($4.5 \times 10^{-5}$; plain differences 68% off), which is why it sits near the limit. With
optimised coefficients the walls-and-water $\rho$ extrapolates to $6.0 \times 10^{-4}$: its
plain differences are 61% and 15% off, falling fourfold, and smaller steps drown in
round-off ($v_p$ goes from $4.0 \times 10^{-5}$ at $h = 1.6 \times 10^{-2}$ to
$5.2 \times 10^{-4}$ at half that), so that is the finite difference's limit at this model's
small $\rho$ gradient. That check is about the geometry (rigid walls, uneven spacing, water)
and runs with Taylor coefficients; the first one covers optimised ones.

**Decision.** The adjoint and gradients are exact with the defaults; no tolerance changed.
The checks refuse receivers next to the layers, so they can't drift into measuring them.

**Reproduce:** `just physics` (`the_adjoint_run_is_the_transpose_of_the_forward_run`,
`the_gradient_matches_finite_differences*`).


## 2026-11-20 — A gradient run's memory, counted beforehand

**Question.** `run` and `check` refuse runs that won't fit, but `gradient` checked only the
forward run's memory, so it could set out on a run many times larger and fail partway. What
does a gradient run hold, and what dominates it?

**Method.** `memory::Plan::gradient_memory` counts what `gradient::gradient` holds at its
peak, while the forward run is replayed: the forward run, the adjoint's solver (a full copy,
material and model included) and its residuals, a copy of the pressure residuals, the
observed traces, $\lceil N/k \rceil$ checkpoints (the wavefield and the CPML memory), $k + 1$
velocity triples for one segment, the correlation sums (8 doubles per padded node) and the
forward and adjoint strains (12 arrays). `gradient::Run::memory_bytes` now sums the arrays a
real run holds at that moment, and `gradient` frees the replay before it makes the gradient,
so the output never adds to the peak. The test compares the two with and without pressure,
the time-dispersion transform and a free surface (with its sponge), at the default interval
and at one that leaves a short last segment.

**Results.** Exact, byte for byte, in all four. Per padded node the checkpoints and the
segment take $4\,(9 \lceil N/k \rceil + 3(k + 1))$ bytes, about $24\sqrt{3N} \approx
42\sqrt N$ at $k = \lceil\sqrt{3N}\rceil$, against about 108 for the whole forward run.
`examples/layered.toml` with the defaults (1,611 steps, checkpoints every 70): 3.63 GiB against
147 MiB, of which 1.85 GiB checkpoints, 1.33 GiB the segment, 145 MiB the adjoint's solver,
102 MiB the sums and 76 MiB the strains. The same model on $1001 \times 1001 \times 401$
nodes: 795 GiB (the forward run: 32 GiB). With the time error kept at $C = 0.3$ (4,456 steps)
the previous entry's estimate was over 6 GiB, so the new defaults' larger step also cut the
gradient's memory by $\sqrt{N}$.

**Decision.** `gradient` refuses runs that don't fit before allocating anything;
`vibe-rations memory --gradient` itemises them; `gradient.json`'s `memory_bytes` is what the
run held, with the observed traces, and equals the plan. At 25 times the forward run's
memory for 1,611 steps, growing as $\sqrt N$, storing the forward field is what limits gradient
runs on one machine: the case for rebuilding it backwards from boundary values (roadmap step
3), for checkpointing that trades memory for recomputation, or for checkpoints on disk. The
adjoint's solver could also share the forward solver's material and model, 44 bytes per node,
which is small next to them.

**Reproduce:** `just unit counts_what_a_gradient_run_allocates`;
`just memory examples/layered.toml --gradient`.


## 2026-11-20 — L-BFGS inversion

**Question.** Do the gradients drive an inversion that recovers the model, and what does it
need besides the gradient?

**Method** (`src/inversion.rs`). The unknowns are relative changes $x = m / m_0 - 1$ of the
inverted parameters at the free nodes (outside the absorbing layers, below `fixed_above`; $v_s$
stays zero in fluids), so the gradient is $m_0\,\partial\chi/\partial m$. Each source is a shot,
and the balanced misfits and their gradients add up. $\Delta t$ and the layers' damping
(`pml_vp_max`) stay at the starting model's, so the forward map is smooth and the gradient is
its derivative. Each iteration takes the L-BFGS direction [@nocedal1980; @liu1989], the
two-loop recursion with initial matrix $\gamma P$, $\gamma = s^T y / y^T P y$
[@nocedal2006], where $P$ is optionally the source-side pseudo-Hessian, $1 / (I / I_{\max} +
\epsilon)$ with $I = \sum_n |v^{n+\frac12}|^2$ the forward run's energy at the node
[@shin2001]. It scales the direction so that no value changes by more than 2% when there is no
memory and 10% ever, projects onto bounds ($v_p$ below the stability limit of the fixed
$\Delta t$, $v_s \le v_p / \sqrt2$), and backtracks with misfit-only runs until the Armijo
condition holds for the projected step. The pair $(s, y)$ is kept if $s^T y > 0$.

**Experiment.** A homogeneous model ($v_p = 2400$, $v_s = 1300$ m/s, $\rho = 2100$ kg/m³) with
a sphere of radius 45 m, $v_p = 2700$, $v_s = 1500$ m/s, 120 m deep;
40 × 36 × 30 nodes of 10 m, a free surface, 8-node layers; three shots of forces at 30 m depth
(12 Hz Ricker), 31 receivers on the surface recording $v_x$, $v_y$ and $v_z$ for 0.4 s;
inverted for $v_p$ and $v_s$ from the homogeneous model. Errors are the L2 model error over
the sphere's nodes, relative to the anomaly:

| Preconditioner | Fixed | Iterations | Misfit left | $v_s$ error | $v_p$ error |
|---|---|---|---|---|---|
| none | nothing | 4 | 54.8% | 99.1% | 100.0% |
| none | top 60 m | 6 | 13.1% | 84.6% | 100.1% |
| illumination, $\epsilon = 10^{-2}$ | nothing | 6 | 23.8% | 92.5% | 100.1% |
| illumination, $\epsilon = 10^{-3}$ | nothing | 6 | 13.5% | 87.9% | 100.2% |
| illumination, $\epsilon = 10^{-2}$ | top 60 m | 6 | 6.8% | 77.8% | 100.0% |
| illumination, $\epsilon = 10^{-3}$ | top 60 m | 6 | 6.6% | 77.7% | 100.1% |

Without help the gradient is dominated by the sources' and receivers' surroundings, where it is
largest, and L-BFGS spends its steps there. The illumination evens out the source side; the
receivers at the surface still pull, and fixing the top 60 m, which holds both, does most.
With the top fixed, $\epsilon$ hardly matters; without it, $10^{-3}$ is better. Every
iteration here took the full step at the first try (in the smaller test below, one of five
backtracked once), so an iteration costs about one gradient (five forward runs) and one
forward run per shot. $v_p$ didn't move: at 12 Hz its
wavelength (200 m) is twice the sphere's diameter, and the forces radiate mostly S waves
(108 m), which the data are full of.

**Decision.** Precondition by default with $\epsilon = 10^{-3}$; fix the region around the
sources and receivers when they sit near the targets' depth. To come: receiver-side
illumination, frequency continuation [@fichtner2013] to approach the larger anomalies, and
scaling between parameters when $\rho$ is inverted too.

**Review.** An adversarial review found the inversion's safeguards incomplete: with $v_p$
inverted alone, $v_p$ could fall until $\lambda < 0$ and then the model failed validation and
aborted the run; with `parameters = ["vs", "vp"]` $v_s$ was capped by the unprojected $v_p$; a
trial that couldn't be set up ended the line search; at a bound the projected search could
stop while the gradient along the bound's face still pointed downhill; and a nearly silent
component could take over the balanced misfit. Now $v_p$ is projected first and held at or
above $\sqrt2\, v_s$ where $v_s$ isn't inverted; failed trials shorten the step; unknowns
pinned at a bound with the gradient pushing out are left out of the direction (projected
gradient), with steepest descent as a fallback; and a velocity component's weight is floored
at $10^{-4}$ of the loudest one's energy. The inversion's own gradient, with respect to its
unknowns ($v_s$, $\rho$, $v_p$, two shots, velocity and pressure), agrees with extrapolated
central differences down the gradient to $1.1 \times 10^{-5}$
(`the_inversion_gradient_matches_finite_differences`); along a random direction the predicted
change cancelled to $4 \times 10^{-4}$ of that, and round-off held the check at $10^{-2}$. The
other tests' numbers didn't change.

**Reproduce:** `an_inversion_finds_a_sphere` (a smaller version: two shots, 36 × 32 × 26
nodes, 0.3 s, the top 60 m fixed, preconditioned; 3 iterations leave 17.9% of the misfit and a
$v_s$ error of 91.8%, 5 iterations 6.6% and 82.7%). The table came from that test's earlier,
larger setup, with the settings varied.

## 2026-11-20 — The inversion with the new defaults

**Question.** The L-BFGS inversion (previous entry) was built with Taylor coefficients, the
time error kept and cubic receivers. Does it still hold with the defaults this branch
adopted, and what memory does it need?

**Method.** Merged it and reran its tests; varied the scheme, the receivers' positions and
the finite-difference step in its gradient check
(`the_inversion_gradient_matches_finite_differences`: two shots, $v_s$, $\rho$ and $v_p$,
velocity and pressure, extrapolated central differences down the gradient).

**Changes.** Optimised coefficients follow the model's slowest wave, so every trial model
had its own operator; `Inversion::new` now pins `optimised_for` at the start's points per
wavelength, beside $\Delta t$ and `pml_vp_max`. The bound on $v_p$ came from Taylor's
stability limit; it now comes from the stencil the runs use.

**Results.** The check failed at $3.1 \times 10^{-5}$ against its $3 \times 10^{-5}$. The
cause was the geometry, not the defaults: with eight-point interpolation the surface line's
end receivers (x = 70 and 220 m: nodes 7 and 22, whose eight-point reach, nodes 4–11 and
19–26, enters the layers at 0–5 and 24–29) read damped values, which breaks the adjoint's exactness (2026-11-20, "The adjoint and gradients with the
new defaults"); Taylor coefficients with the time error kept gave $7.4 \times 10^{-5}$ there.
With the line at 100–190 m, clear of the layers:

| Scheme | $h$ | Extrapolated error |
|---|---|---|
| optimised, time error kept | $4 \times 10^{-3}$ | $1.6 \times 10^{-7}$ |
| Taylor, time error kept | $4 \times 10^{-3}$ | $1.9 \times 10^{-5}$ |
| Taylor, transform | $4 \times 10^{-3}$ | $2.3 \times 10^{-5}$ |
| defaults | $2 \times 10^{-3}$ | $2.7 \times 10^{-5}$ |
| defaults | $4 \times 10^{-3}$ | $5.1 \times 10^{-5}$ |
| defaults | $6 \times 10^{-3}$ | $1.5 \times 10^{-5}$ |
| defaults | $8 \times 10^{-3}$ | $1.4 \times 10^{-5}$ |
| defaults | $10^{-2}$ | $4.0 \times 10^{-6}$ |
| defaults | $1.2 \times 10^{-2}$ | $1.3 \times 10^{-7}$ |
| defaults | $1.6 \times 10^{-2}$ | $2.5 \times 10^{-5}$ |

The plain differences fall as $h^2$ throughout (2.4% at $10^{-2}$, 0.10% at
$2 \times 10^{-3}$), but the extrapolated ones scatter below $h \approx 10^{-2}$ instead of
shrinking: round-off in the misfit of single-precision traces, which takes over as the
difference shrinks; above it, the extrapolation's $h^4$ error grows. Receivers at 100 m depth
in a larger grid gave $2.0 \times 10^{-5}$ (defaults) and $1.7 \times 10^{-5}$ (Taylor, kept)
at $4 \times 10^{-3}$, the same scatter. The sphere inversion (`an_inversion_finds_a_sphere`)
is unchanged: 17.9% of the misfit and 91.9% of the $v_s$ anomaly left after 3 iterations
(91.8% before), every step taken in one trial.

**Memory.** `memory --invert` counts an inversion exactly: one shot's gradient run, every
shot's observed traces, three models (as it stands, the start, the one stepped to), the
gradient and illumination summed over the shots, and per unknown $8(2m + 5)$ bytes for L-BFGS
memory $m$; `counts_what_an_inversion_allocates` matches a real inversion byte for byte.
`Inversion::new` held its set-up run, a whole solver, through the first gradient; it now
drops it, and `iterate` no longer copies the unknowns and gradient for the line search. The
gradient's illumination adds 8 bytes per interior node: `examples/layered.toml`'s gradient
run needs 3.64 GiB (3.63 before).

**Decision.** The check runs at $h = 10^{-2}$ and refuses receivers next to the layers; its
tolerance stays. The inversion pins the operator. When a check like this is a few
$10^{-5}$ off at one step, scan $h$ before suspecting the gradient.

**Reproduce:** `just physics inversion`; `just unit counts_what_an_inversion_allocates`.

## 2026-11-20 — Misfits beyond L2: filtered traces, time shifts and phases

**Question.** Roadmap step 5: the L2 waveform misfit is the only one, and it cycle-skips once
the starting model is half a period off [@bunks1995]. Which misfits to add, and how to keep
the gradient that of the discrete run, checked by finite differences, as the adjoint and the
gradient are? Each misfit must be a function of the sampled traces with its *exact* derivative
as the adjoint source (`gradient::gradient` takes $(\chi, \partial\chi/\partial s)$).

**Method** (`src/misfit.rs`). Three pieces, each derived for the samples, not for continuous
signals, with a literature and numerical check by research agents (Crossref for every
reference; the papers' equations read where open access).

*Low-pass filter.* A digital Butterworth filter [@butterworth1930] of order $n$: the analog
prototype's poles, prewarped, through the bilinear transform, as second-order sections,
$|H(f)|^2 = 1 / (1 + (\tan \pi f \Delta / \tan \pi f_c \Delta)^{2n})$. Applied forwards and then
backwards, starting at rest, with no padding. A causal recursion started at rest is, on $N$
samples, the lower-triangular Toeplitz matrix $H$ of its impulse response; reversing time
($R$) gives $R H R = H^T$, so the forward–backward filter $F = R H R H = H^T H$ is symmetric
and positive semi-definite. Hence the adjoint source of $\chi(F s, F d)$ is
$F\, \partial\chi/\partial(F s)$: the same filter, no separate transpose. Libraries that pad
the ends or start from steady state (e.g. SciPy's `filtfilt`) are not symmetric. The corner
$f_c$ is where the zero-phase filter halves the amplitude.

*Cross-correlation time shift* [@luo1991]. Per trace, $C(k) = \sum_n s[n + k]\, d[n]$ for
$|k| \le K$ (`max_shift`), $k^*$ its largest, and a parabola through
$a, b, c = C(k^* - 1), C(k^*), C(k^* + 1)$:
$$T = (k^* + \delta)\,\Delta, \qquad \delta = \frac{a - c}{2D}, \quad D = a - 2b + c, \qquad \chi = \tfrac12 \sum T^2 .$$
$T > 0$ when the synthetic is late: the sign of [@tromp2005]'s "predicted minus observed",
opposite to [@luo1991]'s $\Delta\tau$ and [@dahlen2000]'s $\delta T$. With
$\partial C(j)/\partial s[m] = d[m - j]$,
$$\frac{\partial T}{\partial s[m]} = \frac{\Delta}{D}\Big[(\tfrac12 - \delta)\, d[m - k^* + 1] + 2\delta\, d[m - k^*] - (\tfrac12 + \delta)\, d[m - k^* - 1]\Big],$$
the observed trace's slope at the shift found, divided by the correlation's curvature: the
discrete form of [@luo1991]'s eq. 7b. For $|k^*| < K$, $b \ge a, c$ gives $D \le 0$ and
$|\delta| \le \frac12$, and where $k^*$ jumps to the next lag both parabolas give the same half
sample: $T$ is continuous, its derivative jumps. At $|k^*| = K$, $\delta$ is held within
$\pm\frac12$ (no derivative while held), so $T$ saturates at $\pm(K + \frac12)\Delta$ instead of
extrapolating; a synthetic that correlates positively nowhere (a silent one) gets that
saturated value, so a model cannot lower the misfit by silencing a trace. Traces quieter than
$10^{-4}$ of the loudest of their kind (velocity, pressure), judged on the observed data only,
are left out, so the selection does not depend on the model.

*Time–frequency phase* [@fichtner2008]. A Gabor transform with a unit-norm Gaussian window of
width $\sigma = 1/f_0$ ("the dominant period", their §2.2.3), truncated at $\pm 3\sigma$, windows
every $\sigma/2$, frequencies every $1/(4\pi\sigma)$ from the first above zero to $1.2 f_{hi}$
(the Nyquist cell, being real, is left out); $f_0$ is the stage's corner frequency or the
sources' dominant frequency, $f_{hi}$ the corner or $2.5 f_0$. They define
$\chi = \int W^2 \Delta\varphi^2$ with $\Delta\varphi = \arg(\tilde u\, \tilde u_0^*)$ and
$W = |\tilde u_0|$ or $\log(1 + |\tilde u_0|)$ normalised. Taken literally on samples it has two
non-differentiable features: a branch cut at $\Delta\varphi = \pm\pi$, and $1/|\tilde u|$ where
the synthetic lacks an arrival that the data have (any scatterer missing from the starting
model), which in a prototype made the gradient about $10^7$ times larger there and the
finite-difference check fail at every step. We compare instead
$$\chi = \tfrac12 \sum_{\text{cells}} W^2 |u - u_0|^2, \qquad u = \frac{z}{\sqrt{|z|^2 + \epsilon^2}},$$
$z$ the synthetic's coefficient, $u_0$ likewise the observed one's, $\epsilon = 0.05 \max|z_0|$
per trace. Where both are large, $|u - u_0|^2 = 2 - 2\cos\Delta\varphi \approx \Delta\varphi^2$; it
is smooth everywhere and bounded where $|z| \to 0$. $W = \text{taper}(f)\, |z_0| / \max|z_0|$
where that ratio is at least 0.05, else 0 (Kristeková et al. weight by the reference's
$|W_r|$ too [@kristekova2009]), from the observed trace alone, normalised to $\sum W^2 = 1$ per component, like the L2 misfit's
balancing. With $q = |z|^2 + \epsilon^2$ and $e = u - u_0$,
$d\chi = \mathrm{Re}(A\, dz)$, $A = W^2 q^{-1/2} (\bar e - \mathrm{Re}(\bar e z)\, \bar z / q)$, and
$dz = G\, ds$ with $G$ the transform's kernel, so the adjoint source is $\mathrm{Re}(G^T A)$.
LASIF's implementation [@krischer2015] was read and not followed: its weight depends on the
synthetic, its denominator is a water level and its resampling is not transposed, so its
adjoint source is not the derivative of its misfit.

**Results** (`just unit misfit`).

- Filter: $\langle F x, y\rangle = \langle x, F y\rangle$ to $8.9 \times 10^{-15}$ (relative); the measured gain
  matches the squared Butterworth response to $10^{-3}$ at orders 1–6; a pulse inside the band
  keeps its peak's sample.
- Time shift, f64: central differences agree with the derivative to $10^{-6}$ (relative)
  along four directions, between lag jumps. Parabola bias, Ricker at $P$ samples per peak
  period, true shifts across a sample, at most 0.049 samples at $P = 5$ (1.0% of a period), 0.011
at 10 (0.11%), $2.8 \times 10^{-3}$ at 20 and $4.4 \times 10^{-4}$ at 50 ($9 \times 10^{-6}$ of a
period), falling as $P^{-2}$. The bias is odd and zero at whole and half
  samples, so it does not move the minimum for consistent data. Against the familiar
  adjoint source $T\,\dot s / \int s \ddot s$ [@tromp2005; @dahlen2000], which assumes $d$ is a
  shifted copy of $s$: it agrees to $1.5 \times 10^{-7}$ for a whole-sample shift of a copy and
$9.3 \times 10^{-4}$ for 3.3 samples, but for traces of different shape (two Rickers against
one) the two point 27° apart (cosine 0.89) and differ by 45%. Using it would fail the
finite-difference check.
- Phase: $\mathrm{Re}\langle G x, a\rangle = \langle x, G^T a\rangle$ to $10^{-12}$; finite
  differences to $10^{-6}$ with an arrival in the synthetic that the data lack; identical
  traces give 0, a synthetic three times too loud $1.1 \times 10^{-4}$ (from $\epsilon$, in the
  faint cells) against 0.29 for one shifted by a tenth of a period that has an extra arrival.
  Cells at 12 Hz and 1 ms
  sampling: windows every 42 samples, 250 samples either side, 37 frequencies.
- Through `Misfit` (f32 traces, weights, filter), all three misfits, filtered or not, agree
  with central differences to $2 \times 10^{-3}$ (f32 round-off).

**Decision.** Keep the parabola (closed form, local, exact derivative); its bias is below
$10^{-3}$ of a period at the sampling these runs record. Use the regularised phase misfit,
documented as such, rather than $\Delta\varphi^2$ literally. Weights and selections come from the
observed data only, so each misfit is one smooth (time shift: piecewise smooth) function of
the synthetics.

**Reproduce:** `just unit misfit`; the numbers above,
`cargo test --release --lib misfit_study -- --ignored --nocapture`.

## 2026-11-20 — Choosing the misfit, and frequency continuation

**Question.** Roadmap step 6 asks for frequency continuation: invert the low frequencies
first, whose misfit has wider basins, then add the higher ones [@bunks1995]; frequency-domain
codes choose the frequencies by the survey's geometry [@sirgue2004], and Fichtner et al.
broaden the period band as the iterations go [@fichtner2013, §3.2.4; their "multiscale" is
mainly about spatial scales]. How should stages and the new misfits (previous entry) enter
the inversion without losing the exactness of its gradient?

**Method** (`src/inversion.rs`, `[inversion] misfit`, `[[inversion.stage]]`). A stage
low-pass filters the observed and the synthetic traces alike at its `max_frequency`, the
standard way to apply a multiscale strategy in the time domain [@bunks1995], and runs its
`iterations` with its own misfit, or the inversion's. Bunks et al. filtered the source
wavelet and the data instead [as described by @boonyasiriwat2009], which lets the low stages
run on coarser grids; here the grid and $\Delta t$ are pinned for the gradient's exactness, so
that advantage is gone, and a filtered source can only give a causal filter, not the zero-phase
one. Filtering the traces keeps one simulation per shot for every stage, and the filter is its
own transpose, so the adjoint source needs nothing more. At each stage's start
(`Inversion::begin_stage`): each shot's misfit is rebuilt from the observed traces (shared,
not copied, when unfiltered); the L2 misfit's balancing weights are taken from the *filtered*
observed traces; the phase misfit's scale $f_0$ is the stage's corner (the sources' dominant
frequency without one); L-BFGS forgets its pairs, since $(s, y)$ pairs describe the curvature
of another misfit [@nocedal2006] (Modrak & Tromp found restarting at each transition faster in
all nine of their test cases, and without it four failed [@modrak2016, §8.1]); and the misfit, gradient and preconditioner are evaluated
afresh at the model as it stands. The unknowns stay relative to the starting model, so stages
continue one another's model. Every stage's filter is checked against the traces' Nyquist
frequency before the first one runs. Without stages there is one stage on the traces as
recorded, so the defaults reproduce the previous inversion bit for bit: the unfiltered L2
misfit is still `l2_misfit` with the same weights.

**Gradient checks.** The inversion's gradient against extrapolated central differences
down the gradient (`the_inversion_gradient_matches_finite_differences`' model: two shots,
$v_s$, $\rho$ and $v_p$, velocity and pressure), relative error, at steps $h$ and $h/2$:

| Misfit | $h = 2 \times 10^{-3}$ | $5 \times 10^{-3}$ | $10^{-2}$ | $2 \times 10^{-2}$ |
|---|---|---|---|---|
| L2, filtered at 10 Hz | $1.4 \times 10^{-5}$ | $1.0 \times 10^{-5}$ | $1.5 \times 10^{-5}$ | $1.8 \times 10^{-4}$ |
| time shifts | $4.3 \times 10^{-5}$ | $2.0 \times 10^{-4}$ | $4.4 \times 10^{-4}$ | $5.6 \times 10^{-4}$ |
| time shifts, filtered at 10 Hz | $7.5 \times 10^{-6}$ | $3.1 \times 10^{-6}$ | $2.2 \times 10^{-4}$ | $1.7 \times 10^{-3}$ |
| phase | $3.9 \times 10^{-5}$ | | $4.2 \times 10^{-6}$ | |
| phase, filtered at 10 Hz | $3.4 \times 10^{-5}$ | | $1.0 \times 10^{-5}$ | |

The filtered L2 and phase misfits behave like the plain L2 one (previous entries): round-off
in single-precision traces below $h \approx 10^{-2}$, the extrapolation's $h^4$ above. The time
shifts do not converge cleanly with $h$: their plain differences at $5 \times 10^{-3}$ were
$7.5 \times 10^{-5}$ and, at half the step, $1.7 \times 10^{-4}$, worse. Each trace's shift is
continuous but its derivative jumps where the correlation's peak moves to the next lag
(previous entry), and with 42 traces some cross a jump within a few $10^{-3}$ of relative
model change. The jumps are small (in a prototype the slope changed by 4%, 2% and 1% at 5, 10 and
20 samples per period, less at the 70 these traces have, 100 filtered), so the check holds at
small steps; a wrong derivative would be off by tens of per cent. Low-pass filtering makes the traces smoother and
the check cleaner.

**Memory.** A filtering stage keeps a filtered copy of every shot's observed traces; the plan
(`memory --invert`) counts it when any stage filters, and `counts_what_an_inversion_allocates`
still matches a real inversion byte for byte, with and without a filtering stage.

**Decision.** The test checks filtered L2 and phase at $h = 10^{-2}$ and the filtered time
shifts at $2 \times 10^{-3}$, all within the existing $3 \times 10^{-5}$
(`the_misfits_gradients_match_finite_differences`). `iterations` beside stages is refused
rather than guessed at.

**Reproduce:** `just physics misfits_gradients`; `just unit stages`;
`just unit counts_what_an_inversion_allocates`.

## 2026-11-20 — Filtering up to the record's end

**Question.** A research agent's check of the stages (previous entry) against ObsPy's
filtering [@beyreuther2010] and the literature found the filter right (symmetric to
$10^{-14}$, the same operator as ObsPy's zero-phase `lowpass`) but inexact near the end of a
record: the forward pass rings on for several periods of the corner after the last sample,
and the record cut that ringing off, so even a trace at rest by its end came out distorted.
How much, and does padding keep the filter its own transpose?

**Method.** Append $M$ zeros to each trace, filter forwards and backwards, and keep the first
$N$ samples: $G = P^T F P$ with $P$ the padding, symmetric because $F$ is. The slowest pole of
an order-$n$ Butterworth filter decays as $e^{-2\pi f_c \sin(\pi/2n)\, t}$, so
$M = \lceil 2.4 / (f_c \Delta \sin(\pi / 2n)) \rceil$ samples bring it below $10^{-6}$ (2091 samples
at 3 Hz, 1 ms, order 4). The start needs nothing: the forward pass runs first, from rest, so
a trace at rest at its start is filtered exactly there. A 12 Hz Ricker filtered at 3 Hz was
compared with the same pulse in a record ten times as long.

**Results.** Largest error relative to the filtered pulse's peak:

| Pulse ends before the record does | Without zeros | With them |
|---|---|---|
| about $0.5 / f_c$ | 0.56 | $8.1 \times 10^{-14}$ |
| about $3 / f_c$ | $1.4 \times 10^{-3}$ | $3.0 \times 10^{-14}$ |

The padded filter stays symmetric to $10^{-12}$ (unit test). Through the inversion's
finite-difference check (previous entry's model, whose 0.25 s records at a 10 Hz corner are
2.5 periods of the corner long, so their ends mattered), the filtered L2 misfit's gradient now
agrees to $4.4 \times 10^{-7}$ at $h = 10^{-2}$ ($1.5 \times 10^{-5}$ before). The filtered time
shifts changed with the traces, and at $h = 2 \times 10^{-3}$, which had happened to miss every
lag jump before, the coarse difference now crossed one ($8.5 \times 10^{-4}$). Scanning $h$:
$2.2, 1.4, 1.3 \times 10^{-4}$ at $h = 2.5 \times 10^{-4}, 5 \times 10^{-4}, 10^{-3}$, no trend,
which is round-off in the shifts of single-precision traces; $8.5 \times 10^{-4}$ and
$2.0 \times 10^{-3}$ at $2 \times 10^{-3}$ and $3 \times 10^{-3}$, the jumps.

**Decision.** Pad. The time shifts' check runs at $h = 10^{-3}$ against $5 \times 10^{-4}$: their
derivative is exact between lag jumps, a wrong one would be off by tens of per cent, and a
smooth sub-sample estimate (Newton on the correlation's band-limited interpolant, which a
research prototype found unbiased and free of jumps) would cost all $2N - 1$ lags of every
trace without an FFT; it stays an option if the jumps hamper L-BFGS. Records shorter than a
few periods of a stage's corner are still edge-dominated: a cut arrival can't be filtered
correctly by any filter.

**Reproduce:** `just physics misfits_gradients`; the table,
`cargo test --release --lib misfit_study -- --ignored --nocapture`.

## 2026-11-21 — The misfits and stages on two test inversions

**Question.** Do the new misfits and frequency-continuation stages (previous entries) help
the inversion, and where? Frequency continuation and traveltime-like misfits exist to widen
the basin of attraction when the starting model is more than half a period off
[@bunks1995; @luo1991; @fichtner2008]; with a start inside the basin, L2 should do at least as
well.

**Method** (`tests/inversion_study.rs`, `just study`). Two surveys, each inverted from the
same start by every variant; after each stage, the balanced L2 misfit of the model (one scale
for every variant, relative to the start's) and the model error over a target, relative to
the start's error.

- *Sphere*: `an_inversion_finds_a_sphere`'s survey (two force shots at 12 Hz, 50 surface
  receivers, the top 60 m fixed), $v_p$ and $v_s$ from the homogeneous start, 6 iterations;
  the target is the sphere.
- *Slow start*: a homogeneous model ($v_p = 2400$, $v_s = 1300$ m/s), started 20% too slow
  in both, 60 × 36 × 28 nodes of 10 m, two force shots at 30 m depth and 400 m apart, 43
  receivers on a 420 m surface line, 0.5 s, nothing fixed, 8 iterations. At the far offsets
  (about 410 m) the start's P waves come about half a period (43 ms at 12 Hz) late and its S
  waves about a full period (79 ms); the S waves are half a period late at about 220 m. The
  target is every free node below 60 m.

**Results.**

| Sphere, 6 iterations | L2 left | $v_s$ error | $v_p$ error |
|---|---|---|---|
| L2 | 5.4% | 81.1% | 100.9% |
| L2 below 6 Hz (3), then all (3) | 36.1% → 7.1% | 82.8% → 78.9% | 98.9% |
| time shifts (3), then L2 (3) | 61.1% → 9.7% | 95.9% → 86.5% | 100.6% |
| phase (6) | 12.4% | 86.8% | 101.5% |
| phase (3), then L2 (3) | 32.5% → 7.9% | 94.8% → 84.4% | 100.9% |

| Slow start, 8 iterations | L2 left | $v_s$ error below 60 m | $v_p$ error |
|---|---|---|---|
| L2 | 10.5% | 101.7% | 100.3% |
| L2 below 4 Hz (3), 8 Hz (3), then all (2) | 86.3% → 55.4% → 35.1% | 100.0% | 100.2% |
| time shifts (4), then L2 (4) | 67.0% → 16.6% | 99.3% | 99.9% |
| phase (4), then L2 (4) | 88.4% → 29.2% | 100.7% | 100.1% |

Each variant took 70–80 s (sphere) or 170–290 s (slow start) on 4 cores shared with other
work.

On the sphere, whose start is well inside L2's basin, L2 fits the data best; going through
6 Hz first ends with the smallest $v_s$ error (78.9% against 81.1%) at a little more misfit.
The time shifts and phases ignore amplitudes, so they fit less of the L2 misfit, and three
L2 iterations afterwards don't make up for it. On the slow start the test does not
discriminate: no variant changed the deep model in 8 iterations (errors 99–102%), while L2
lowered its misfit to 10.5%: the updates went where the preconditioned gradient is largest,
around the sources and receivers, which nothing held fixed here. The low-passed stages kept
little of the data's energy (a 12 Hz Ricker keeps about 0.5% of its energy below 4 Hz, as a
research agent measured for the filter) and moved slowly.

**Decision.** Keep L2 on the traces as recorded as the default; the misfits and stages are
options, with exact gradients, for starts outside L2's basin. Testing that they widen the
basin needs a survey with more propagation (several wavelengths to the targets) and more
iterations than this machine affords for a test, and a reference for what reaches the deep
model (fixing the near-surface, or the receivers' side of the illumination; next entries).

**Reproduce:** `just study misfits_on_a_sphere`; `just study slow_start` (about 6 and 15
minutes).

## 2026-11-21 — Scaling each parameter's steps

**Question.** Roadmap step 6 asks for parameter scaling for $\rho$. In multiparameter FWI the
classes' gradients differ in size and in how well the data constrain them, and an inversion
"is very sensitive to the parameter scaling" even with an approximate Hessian
[@lavoue2014]; density is the least resolved of the elastic parameters [@kohn2012;
@operto2013]. What should the scaling be here, where the unknowns are already relative
changes $x = m/m_0 - 1$?

**Literature** (read by a research agent on the publishers' pages; equations as images
were not read). Prieux et al. normalise each class by its mean in the starting model, so that
the classes share an order of magnitude, and start L-BFGS from a Gauss–Newton diagonal
[@prieux2013]; that is what relative unknowns already do. Kamei & Pratt scale the second
class's gradient by a factor $\xi$ chosen from trade-off curves (0.025–0.075 for
attenuation), although the gradients were "similar in magnitudes" [@kamei2013]. Köhn et al.
"normalize the material parameters and the gradients" without saying how in the readable
text [@kohn2012].

**Method** (`[inversion.scale]`, `Inversion::preconditioner_for`). A factor $w_c$ per
parameter multiplies the preconditioner: $P' = D^2 P$, $D = \mathrm{diag}(\sqrt{w_c})$. L-BFGS
with initial matrix $\gamma P'$ and pairs $(s, y)$ gives exactly $D$ times the direction of
L-BFGS on $z = D^{-1} x$ with pairs $(D^{-1}s, D y)$, gradient $D g$ and $\gamma P$
[@nocedal2006]: $s^T y$ and $\gamma$ are unchanged, and so are the Armijo test (which uses
$g^T s$) and the bounds, which act on $x$. So a factor is a change of variables for its
class, and the pairs stay valid if the factors change (L-BFGS converges for any bounded,
positive-definite initial matrix [@liu1989]). Without the illumination, the factors alone
make the preconditioner. The first step's largest change in class $c$ becomes
`initial_step` times $w_c \max_c |P g| / \max_{c'} (w_{c'} \max_{c'} |P g|)$. A unit test
checks the change of variables to $10^{-14}$; factors of 1 leave every result bit for bit.

**Experiment.** The sphere survey inverted for $v_s$, $\rho$ and $v_p$, the top 60 m fixed,
5 iterations, with $\rho$'s factor 1, 0.3 and 0.1: once with a denser core too ($\rho$ 2300
against 2100 kg/m³), and once with none, where any change of $\rho$ is cross-talk. Errors
over the sphere, relative to the start's; for the second, the rms relative change of $\rho$
there:

| $\rho$'s factor | L2 left | $v_s$ error | $v_p$ error | $\rho$ error (denser core) |
|---|---|---|---|---|
| 1 | 6.6% | 85.0% | 100.4% | 79.6% |
| 0.3 | 6.1% | 84.1% | 100.0% | 92.3% |
| 0.1 | 6.6% | 83.2% | 99.8% | 97.1% |

| $\rho$'s factor | L2 left | $v_s$ error | $v_p$ error | $\rho$ changed (no anomaly) |
|---|---|---|---|---|
| 1 | 6.2% | 86.0% | 101.1% | 1.78% |
| 0.3 | 7.4% | 82.5% | 101.1% | 0.74% |
| 0.1 | 6.3% | 83.4% | 101.0% | 0.26% |

Where density has no anomaly, a smaller factor cuts its cross-talk sevenfold and the $v_s$
anomaly is recovered a little better; where it has one, the same factor leaves it mostly
unrecovered (97% at 0.1). The factor is a trade-off between the two, as Kamei & Pratt found
for attenuation [@kamei2013], and nothing in the data alone says which case holds.

A research agent also measured, in a copy of this code, that the relative unknowns already
balance the classes: at the start, the largest preconditioned gradient of $\rho$ was 0.70–0.78
of $v_s$'s, and $v_p$'s 1/37 of it (its wavelength is twice the sphere's). Equalising the
classes' largest (or rms) preconditioned gradients automatically made it worse: $v_p$, barely
constrained, was amplified 37-fold into cross-talk ($v_p$ error 109–115%, misfit no better).

**Decision.** Fixed factors, 1 by default, set per inversion (`[inversion.scale] rho = 0.3`,
say) where density takes up what the velocities should explain; no automatic equalisation.
Choosing them is the trade-off-curve exercise of Kamei & Pratt.

**Reproduce:** `just study density_scaling` (about 5 minutes); the unit tests
`preconditioner_factors_change_variables` and `scale_multiplies_each_parameters_preconditioner`.

## 2026-11-21 — The receivers' side of the illumination

**Question.** The preconditioner divides by the sources' side of the illumination, the
forward wavefield's energy $I_s$ at each node [@shin2001], which leaves out the receivers'
Green's functions; without fixing the near-surface, the updates still pile up around the
receivers (the L-BFGS entry). Can the receivers' side be added cheaply, and does it help?

**Method** (`receiver_runs`, `Inversion::receiver_energy`). In the Born approximation the
Gauss–Newton Hessian's diagonal is $\sum_s \sum_a w_a \sum_t (\partial d_{s,a} / \partial m)^2$,
which for pulses of similar spectra is about $I_s(n)\, I_r(n)$ with
$I_r(n) = \sum_a w_a \sum_t |G_a(n, t)|^2$, $G_a$ the waves sent from receiver component $a$
[@plessix2004, whose "type 1" weights, assuming $I_r$ constant, are Shin et al.'s]. Sending
from every component at once with random signs $\varepsilon_a = \pm 1$, as encoded simultaneous
sources [@krebs2009], gives $|\sum_a \varepsilon_a \sqrt{w_a} G_a|^2$, whose expectation is
$I_r$: the cross terms average out. Each run injects the first shot's wavelet through
`Solver::adjoint_op`, as the adjoint run injects residuals, so velocity and pressure
receivers weigh as in the gradient, with $w_a$ each component's weight in the misfit
relative to the largest, averaged over the shots; signs come from SplitMix64 [@steele2014],
seeded by the stage. The estimate, averaged over `receiver_runs` runs, is made at each
stage's start in the model as it stands and held; the preconditioner becomes
$w_c / (I_s I_r / \max + \epsilon)$. Each run costs a forward run, so 8 runs cost about two
fifths of an iteration of this two-shot survey, and less with more shots. A unit test checks
the estimator: with the signs of a 4 × 4 Hadamard matrix's rows, four runs average to the
sum of the four components sent alone to $10^{-5}$ of the largest. A research agent compared
the estimate with the exact sum of 150 single-receiver runs: the log-standard deviation over
the free nodes was 0.36 for 1 run, 0.215 for 4 and 0.099 for 16 in a homogeneous model, and
an analytic $\sum_r 1 / r^2$ missed the depth decay of a layered one by about twofold over
100 m.

**Experiment** (`just study receivers_side`). The sphere survey ($v_p$, $v_s$), 5 iterations:

| | L2 left | $v_s$ error | $v_p$ error |
|---|---|---|---|
| nothing fixed, sources' side | 14.1% | 92.3% | 100.5% |
| nothing fixed, both sides, 2 runs | 8.8% | 87.3% | 100.7% |
| nothing fixed, both sides, 8 runs | 9.5% | 88.6% | 100.7% |
| top 60 m fixed, sources' side | 6.5% | 82.8% | 100.9% |
| top 60 m fixed, both sides, 8 runs | 6.7% | 82.2% | 101.1% |

With nothing fixed, the receivers' side takes over part of what fixing the top does: a third
less misfit and a better $v_s$ anomaly. The estimate's accuracy hardly matters (2 and 8 runs
alike; the research agent found the exact sum and an analytic one within 0.2% of misfit of
each other). With the top fixed, where the receivers' factor varies only 4–8-fold over the
free nodes, it changes nothing beyond noise. The source side is still recomputed every
iteration and the receiver side once per stage.

**Decision.** An option, off by default (`receiver_runs = 0`) so the baselines stand; use it
when the near-surface can't be fixed. 2–8 runs are enough.

**Reproduce:** `just study receivers_side` (about 4 minutes); the unit test
`the_receivers_side_adds_each_receivers_energy`.

## 2026-11-21 — Review of the misfits, stages and preconditioners

**Method.** An adversarial review of this branch's six commits: four reviewers (the misfits'
mathematics; stages and settings; scaling and the receivers' side; documentation and
citations), each finding then handed to a separate agent told to refute it. 21 findings, 19
confirmed (14 distinct, some found twice), 2 refuted.

**Results.** None concerned a gradient: every derivative the reviewers checked, by
derivation or by finite differences, held. Confirmed and fixed:

- *Receivers' side under the time shifts and phases.* Those misfits scale every trace by its
  own amplitude, so their adjoint sources, and the Gauss–Newton diagonal, weigh a trace about
  1/its energy, while the receivers' side weighed each component 1: pressure, $10^6$ times
  velocity, took over the estimate. Each trace now has its own weight
  (`Misfit::receiver_weights`: L2's as before, $1/\sum \dot d^2$ for the time shifts,
  $1 / (E \cdot \text{traces measured})$ for the phases). The previous entry's runs used L2, so
  its numbers stand; its sentence that receivers "weigh as in the gradient" held only for L2.
- *End padding near the Nyquist frequency.* The number of zeros came from the analog
  prototype's decay, which the bilinear transform slows near the Nyquist frequency; it now
  follows the digital filter's largest pole radius $r$, $r^M < 10^{-6}$ with 10% margin (2107
  samples at 3 Hz, 1 ms, order 4, where the entry on filtering up to the record's end has
  2091), and corners above 0.8 of the Nyquist frequency are refused.
- *Memory.* The phase misfit's kernels and the misfits' weights weren't counted; the plan now
  works out each stage's misfit arrays (`Misfit::held_bytes`) and takes the largest, and
  counts the preconditioner only when an inverted parameter's factor differs from 1.
- *Robustness.* One-sample traces with `max_shift` set panicked in any misfit; scale factors
  weren't checked (a negative one would step uphill); a stage that took no step was missing
  from `inversion.json`.
- *Documentation.* Only the scale factors' ratios matter (L-BFGS's $\gamma$ and the first
  step's scaling absorb a common factor), which the docs now say. In the entry on two test
  inversions, the slow start's far-offset S waves come about a full period late (79 ms), not
  half (P waves 43 ms; S half a period at about 220 m), so L2 there faced cycle skipping in
  the S waves, and still lowered its misfit, near the sources and receivers. Corrected in
  place: 42 traces, not 36, in the stages entry; the phase test's shifted synthetic is a
  tenth of a period off and has an extra arrival; the README's "every trace alike" leaves out
  quiet traces; `just study sphere` also ran the density study, so the studies are renamed
  (`misfits_on_a_sphere`, `density_scaling`).

Refuted: that commit b3a9943's "28–39%" contradicts the 45% in the log (the commit quotes a
research prototype's range, the log our own test); that the first misfit entry describes the
filter as unpadded (it was, at the time; the padding entry follows it).

**Reproduce:** the unit tests `receiver_weights_follow_the_adjoint_sources`,
`pads_enough_near_the_nyquist_frequency`, `short_traces_are_refused_only_where_they_must_be`,
`counts_what_it_holds`, `counts_what_an_inversion_allocates` and
`scale_multiplies_each_parameters_preconditioner`; `tests/cli.rs`'s invert test.

## Open questions

- **A test that discriminates the misfits.** The slow start (2026-11-21) left the deep model
  untouched with every misfit in 8 iterations; a survey with more wavelengths to the targets,
  more iterations and the near-surface held (or the receivers' side of the illumination) is
  needed to show whether time shifts, phases or frequency continuation widen the basin here
  [@bunks1995; @luo1991; @fichtner2008]. Windowing single arrivals, which the time shift
  assumes for traces with several, is still to come.
- **Smoother time shifts.** The parabola's shift is continuous but its derivative jumps
  where the correlation's peak changes lag, which limits the finite-difference check to
  $10^{-4}$; Newton's method on the correlation's band-limited interpolant would be smooth and
  unbiased at the cost of every lag of every trace.
- **The phase misfit's scale** is the stage's corner frequency; a low-passed spectrum peaks
  at about 0.85 of it (a research agent's measurement), so $\sigma = 1/f_c$ may be 15% short.
  The receivers' side ignores a stage's filter too.
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
  refinement (the band shrinks with $h$, so the operator tends to Taylor's); and choose the
  band an inversion pins (the starting model's slowest wave, or the lowest speed allowed).
  Gradients with them agree with finite differences (2026-11-20, the adjoint and gradients
  with the new defaults).
- **The time error in gradients.** Residuals enter through the transposed transform, and the
  gradient agrees with finite differences of the transformed run (same entry). Since the
  transformed traces equal the time-continuous ones for every model, the gradient should
  too: compare with a small-$C$ reference. And whether order 12 or 16 should be the
  default.
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
