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

## 2026-11-20 — Rebuilding the forward wavefield backwards

**Question.** A gradient run replays the forward wavefield from checkpoints and keeps a
segment's velocities, about $42\sqrt N$ bytes per padded node: 3.64 GiB for
`examples/layered.toml`, 25 times its forward run. The interior scheme is reversible, so the
field could instead be rebuilt backwards from the final state and values saved around the
absorbing layers [@dussaud2008; @aaker2020; @nguyen2015 compare such schemes in 2D]. How much
does that save, where does it lose, and how far does single precision let it drift? The
alternative it competes with is checkpointing [@griewank2000; @symes2007].

**Step 1: a step back** (`Solver::step_back_with`, `Band`). A forward step is: velocity
update, forces, velocity mirror, stress update (plane-stress moduli and $\sigma_{zz} = 0$ on
the free surface), moments, sponge, stress mirror. The step back undoes these in reverse
order: moments, the stress update, then the band's stresses and the stress mirror, forces,
the velocity update, then the band's velocities and the velocity mirror. Each increment is
computed from the same values as forward and subtracted: the kernels take a sign of $\pm 1$,
which leaves the forward step bit for bit as it was, and a source withdraws exactly what it
added. The CPML and the sponge act only in the layers, which the step back doesn't rebuild,
and $\sigma_{zz}$ on the surface was zero before the step as after it. So in exact
arithmetic the region outside the layers, $I$, comes back exactly, given what the stencils
read from outside it: values within $L$ nodes of $I$ along an axis ($L$ the stencil's
half-length). The band keeps the velocities on that whole shell, edges included, and each
stress only on the faces across which the velocity update differentiates it (6 values per face
node). The velocities on the edges are there for the gradient: the strain rates at positions
just outside $I$ that feed its first nodes (the shear moduli's harmonic means, the density's
arithmetic means) read velocities along the faces. Above a free surface the mirrors rebuild
everything. Records are taken before each forward step, so the step back from $n + 1$ to $n$
reads the record of state $n$. Stepping back 80 steps from the end of a run with a force and an
explosion on a $30^3$ grid, 5-node layers, retraced $I$ to $1.5 \times 10^{-6}$ of each array's
peak, with a free surface and with layers all round alike
(`simulation::tests::steps_back_inside_the_band`).

**Step 2: in the gradient** (`gradient_with`, `Strategy::Reconstruct`). The forward field is
needed in reverse order: the stiffness term pairs $\tilde\varepsilon^{n+1}$ with the strain
rate of $v^{n+\frac12}$, the density term uses $v^{n+\frac12} - v^{n-\frac12}$. So each
backward step copies $v^{n+\frac12}$, steps back, and correlates with the copy and the rebuilt
$v^{n-\frac12}$. Saving the band for the whole run is one way; the other, which the memory
below argues for, saves it for segments of $k$ steps, replaying each segment but the last
from a checkpoint at its start (the last is saved on the way forward) and rebuilding it from
its end, which the replay reaches bit for bit. The error then accumulates over a segment, not
the run, at the cost of one more pass of the solver for all but the last segment. On the
small run of the gradient's unit tests (21 steps, two layers, a force; free surface; 4-node
layers) the rebuilt gradient at every node of $I$ matches the replay's to
$7.4 \times 10^{-8}$ (segments of 1 step), $3.5 \times 10^{-7}$ (7), $5.7 \times 10^{-7}$ (the
default) and $4.0 \times 10^{-6}$ (the whole run) of its largest value
(`gradient::tests::rebuilt_gradients_match_the_replay`). Inside the layers the rebuilt field,
and the gradient, mean nothing; the inversion only updates nodes in $I$.

**Step 3: memory.** A band record holds $3\,|\text{shell}| + 3\,\Sigma\,\text{faces}$ values,
about $6L$ per node of $I$'s surface; a checkpoint 9 arrays and the CPML memory. With
segments of $k$ steps the reconstruction holds $\lceil N/k\rceil - 1$ checkpoints and $k$
records, least at $k \approx \sqrt{N C / B}$ (`reconstruction_segment`); the plan counts it
exactly (`counts_what_a_rebuilding_gradient_run_allocates`). The gradient's peak:

| Grid | Checkpoints | Rebuilt, segments | Rebuilt, whole run |
|---|---|---|---|
| `examples/layered.toml`, $129 \times 129 \times 81$, 1,611 steps | 3.64 GiB | 1.68 GiB ($k = 202$) | 5.67 GiB |
| the same on $1001 \times 1001 \times 401$ | 798 GiB | 245 GiB ($k = 323$) | 466 GiB |

A record of `layered.toml` is 862,968 values (3.3 MiB), against 60 MiB for a checkpoint and
19 MiB for a step of the replayed segment. Saving the band for the whole run loses on the
small grid, where its surface is large next to its volume ($B N$ grows as $N$, the
checkpoints as $\sqrt N$), and wins on the large one; segments win on both, by 54% and 69%,
for one more pass.

**Decision so far.** Offer it as `gradient --reconstruct [--segment N]` beside the default,
checkpoints, which stay the default (also in the inversion) until the drift is measured.

**Step 4: checks** (`tests/physics.rs`). The rebuilt field against the run it retraces, its
largest difference outside the layers over the run relative to each array's peak: order 4
(251 steps) $7.5 \times 10^{-7}$ along a free surface and $6.3 \times 10^{-7}$ with layers all
round; order 8 (283 steps) $6.2 \times 10^{-7}$ and $1.2 \times 10^{-6}$; water under the
surface with rigid walls (213 steps) $9.6 \times 10^{-7}$
(`the_rebuilt_wavefield_follows_the_run`, tolerance $5 \times 10^{-6}$). The rebuilt gradient
against the replay's on the same models: $2.9$–$7.1 \times 10^{-7}$ of each parameter's
largest value, default segments or the whole run (`rebuilt_gradients_match_the_replay`,
$3 \times 10^{-6}$). The finite-difference checks on the reconstruction path give the replay's
numbers within noise (rebuilt, then replayed): order 4 $v_p$ $2.6 \times 10^{-4}$ /
$2.8 \times 10^{-4}$, $v_s$ $2.1 \times 10^{-5}$ / same, $\rho$ $4.4 \times 10^{-5}$ / same;
order 8 $v_p$ $1.4 \times 10^{-6}$ / $3.9 \times 10^{-6}$, $v_s$ $1.0 \times 10^{-5}$, $\rho$
$1.6 \times 10^{-5}$; walls and water $v_p$ $6.3 / 6.4 \times 10^{-5}$, $v_s$
$1.2 \times 10^{-5}$, $\rho$ $1.6 \times 10^{-4}$
(`the_gradient_matches_finite_differences_when_rebuilt`). The other baselines are unchanged:
the forward step is bit for bit what it was.

**Step 5: drift** (`just drift`, `tests/reconstruction.rs`). A $48 \times 40 \times 36$ model
(10 m; two layers, a sphere, a force and an explosion at 14 Hz; 8-node layers;
`optimised_for = 4`), the drift relative to each array's peak over the run after stepping
back the given number of steps from the end. By run length (order 8, $C = 0.85$, free
surface, 3.6 s, 2,905 steps):

| Steps back | 1 | 10 | 100 | 1000 | 2000 | 2905 (to the start) |
|---|---|---|---|---|---|---|
| Drift | $3.3 \times 10^{-12}$ | $4.5 \times 10^{-12}$ | $1.3 \times 10^{-11}$ | $3.5 \times 10^{-11}$ | $9.8 \times 10^{-11}$ | $3.3 \times 10^{-7}$ |

It grows where the field is loud, not with the number of steps. By forward step, with the
field's own largest value on the same scale:

| Step | 1000 | 500 | 300 | 200 | 100 | 50 | 0 |
|---|---|---|---|---|---|---|---|
| Drift | $6.8 \times 10^{-11}$ | $8.4 \times 10^{-10}$ | $5.8 \times 10^{-8}$ | $1.3 \times 10^{-7}$ | $4.8 \times 10^{-7}$ | $3.7 \times 10^{-7}$ | $3.3 \times 10^{-7}$ |
| Field | $3.5 \times 10^{-5}$ | $2.4 \times 10^{-3}$ | $6.7 \times 10^{-2}$ | $0.15$ | $0.81$ | $0.29$ | 0 |

One step back from the exact final state: $6.4 \times 10^{-13}$. The drift rises from
$6 \times 10^{-8}$ to $5 \times 10^{-7}$ while the rebuilt field passes through its loudest
200 steps, about ten units of single-precision round-off ($6 \times 10^{-8}$) at that
amplitude, and stays as an absolute error once the field is quiet again. By stencil order and
Courant fraction, 1.2 s, all the way back: order 4 $7.0 \times 10^{-7}$ (860 steps), 8
$6.2 \times 10^{-7}$ (969), 12 $8.1 \times 10^{-7}$ (1,006), 16 $6.7 \times 10^{-7}$ (1,027);
order 8 at $C = 0.5$ $1.1 \times 10^{-6}$ (1,646 steps) and at $C = 0.3$ $7.6 \times 10^{-7}$
(2,743). Neither matters beyond a factor of two; three times the steps through the same loud
phase add little. By boundary, order 8, 969 steps (100 and 300 steps back, then all the way):

| Boundary | 100 | 300 | All the way |
|---|---|---|---|
| Layers, free surface | $4.2 \times 10^{-11}$ | $1.6 \times 10^{-10}$ | $6.2 \times 10^{-7}$ |
| Layers all round | $9.8 \times 10^{-12}$ | $2.1 \times 10^{-11}$ | $4.7 \times 10^{-7}$ |
| Rigid walls, no layers | $4.6 \times 10^{-8}$ | $7.4 \times 10^{-8}$ | $6.7 \times 10^{-7}$ |
| Rigid walls, free surface | $5.3 \times 10^{-8}$ | $9.4 \times 10^{-8}$ | $8.0 \times 10^{-7}$ |
| Layers, 60 m of water at the surface | $2.7 \times 10^{-9}$ | $9.0 \times 10^{-9}$ | $3.5 \times 10^{-6}$ |

With rigid walls the band is the halo, always zero, and no layers or sponge are involved, so
their drift is the scheme's round-off alone: it is the same at the start as with layers and a
free surface. It is larger at 100–300 steps back only because the walls keep the field loud
throughout. So the layers and the free surface add nothing measurable: in exact arithmetic
the band makes the rebuilt field exact, and the whole drift is single precision. Water is the
exception at five times the rest, $3.5 \times 10^{-6}$, which we haven't traced; in a fluid
the stresses are pressure alone and the velocities large, which changes how round-off in each
array weighs.

**Step 6: in the gradient, and cost.** The rebuilt gradient against the replay's, outside the
layers, on the same model (969 and 2,905 steps), with memory and time on the container's 4
cores:

| | 969 steps | | | 2,905 steps | | |
|---|---|---|---|---|---|---|
| | off by | memory | time | off by | memory | time |
| Checkpoints | | 208 MiB | 4.3 s | | 335 MiB | 9.6 s |
| Segments of 10 | $1.7 \times 10^{-7}$ | 574 MiB | 5.0 s | $8.1 \times 10^{-8}$ | 1.62 GiB | 12.7 s |
| Segments of 100 | $1.8 \times 10^{-7}$ | 125 MiB | 4.9 s | $1.7 \times 10^{-7}$ | 237 MiB | 12.0 s |
| Default (108; 194) | $2.4 \times 10^{-7}$ | 123 MiB | 4.8 s | $2.4 \times 10^{-7}$ | 191 MiB | 12.1 s |
| Whole run | $1.4 \times 10^{-6}$ | 424 MiB | 4.0 s | $1.3 \times 10^{-6}$ | 1.17 GiB | 10.2 s |

(Times from a second run of the study; the first gave 3.5 s for checkpoints, 4.0 s for the
default segments and 3.3 s for the whole run at 969 steps, and 10.2, 12.2 and 10.3 s at 2,905:
on a grid this small they vary by about 20%.) Segments reset the drift at each checkpoint, so the gradient's error is set by one
segment's loudest stretch: $2 \times 10^{-7}$ with the default against $1.4 \times 10^{-6}$ for
the whole run. On `examples/layered.toml` ($129 \times 129 \times 81$, 1,611 steps, a gradient
for the upper crust's $v_s$ lowered from 2600 to 2500 m/s; `gradient --reconstruct`):

| Strategy | Memory | Time | $v_p$ | $v_s$ | $\rho$ |
|---|---|---|---|---|---|
| Checkpoints, every 70 | 3.64 GiB | 70.2 s | | | |
| Rebuilt, segments of 202 | 1.68 GiB | 80.2 s | $2.6 \times 10^{-6}$ | $4.9 \times 10^{-7}$ | $5.7 \times 10^{-7}$ |
| Rebuilt, the whole run | 5.67 GiB | 72.3 s | $3.9 \times 10^{-6}$ | $8.9 \times 10^{-7}$ | $6.9 \times 10^{-7}$ |

The last three columns are the largest difference from the replay's gradient outside the
layers, relative to its largest value (from the written single-precision arrays). The
misfits agree to every printed digit. Segments cost 14% more time here (one more pass over
all but the last segment), the whole run 3% (a step back is a little cheaper than a step:
no CPML), and on the small grid above 12–26% and -7 to +6%.

**Decision.** Keep checkpoints the default everywhere, the inversion included; the
reconstruction is an option (`gradient --reconstruct`), not bit-exact. Its error is
single-precision round-off, $10^{-7}$–$10^{-6}$ of the gradient's largest value, set by how
loud the field gets within a segment, not by run length, order, time step, layers or free
surface; for an inversion that is far below the finite differences' own resolution and
the forward modelling's error. Use segments (the default) to halve the gradient's memory on
grids like `examples/layered.toml`, and to cut it by two thirds on large ones, at about 15%
more time; save the band for the whole run only on large grids, where it needs less memory
than checkpoints and no extra pass. Switching the default is a separate decision, once the
inversion has been run both ways.

**Reproduce:** `just physics rebuilt` (the checks); `just drift` (steps 5 and 6); `just memory
examples/layered.toml --gradient [--reconstruct [--segment 100000]]` (step 3); for the last
table, `just run examples/layered.toml --out out/rebuild/observed`, then `just gradient
guess.toml out/rebuild/observed [--reconstruct [--segment 100000]]` with `guess.toml` the
example with the upper crust's `vs = 2500.0`, and compare the `grad_*.npy`.

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

## 2026-11-20 — Marmousi2

**Question.** How fast, and how large, are runs on a realistic model? Marmousi2
[@martin2006], the elastic and extended successor of Marmousi [@versteeg1994]: forward and
gradient runs at orders 4 and 8, their speed and memory against what `memory` plans.

**Data.** Three SEG-Y files (vp, vs, density), each 155,653,444 bytes: 13,601 traces of
2,801 samples, IBM floats (format 1), km/s and g/cm³, from the Madagascar project's mirror
(`https://www.ahay.org/data/marm2/`); the SEG wiki's S3 link answered 403 and its licence
page could not be read. SHA-256 sums in `examples/marmousi2/marmousi2.sha256`. The headers
don't record the spacing; 1.25 m in both directions gives the published 17 × 3.5 km. Water
fills the first 361 samples (to 450 m) everywhere; $v_p$ spans 1,028–4,700 m/s, $\rho$
1,010–2,627 kg/m³, and $v_s$ is as low as 270 m/s in the soft sediments under the seafloor
(291–500 m/s down to about 1 km, $v_p/v_s \approx 4$–5).

**Method.** `vibe-rations extrude` (`src/section.rs`) averages each property over the
node's cell, $h \times h$ centred on it: $\rho$ arithmetically, the moduli $\rho v_p^2$ and
$\rho v_s^2$ harmonically, the long-wave averages across fine layers [@backus1962]. A node
is fluid when most of its cell is; Backus's average would make fluid every node a water
sample reaches and move the seafloor down by up to half a node. At 25 m the seafloor falls
between the nodes at 450 m (water) and 475 m (the first solid row), the slowest $v_s$ is
303 m/s and $v_p$ keeps its range. The section is repeated along $y$: 41 nodes, 9 of them
outside 16-node layers, enough for receivers at the middle plane. For an inversion,
`--smooth` smooths the solid part by a normalised convolution over the solid nodes only, so
the water stays sharp (a Gaussian, $\sigma = 250$ m, takes $v_s$'s minimum to 407 m/s).

- Whole section (`examples/marmousi2/marmousi2.toml`): $681 \times 41 \times 141$ = 3.94 M
  nodes, an explosion 50 m below the surface at $x = 8.5$ km, Ricker at 1.2 Hz (4.0 points
  per shortest wavelength at 3 Hz), 641 ocean-bottom nodes ($v_x$, $v_z$) on the first solid
  row, 6 s.
- Central 6 km ($x$ = 6–12 km; `central-true.toml`, `central.toml`): $241 \times 41 \times
  141$ = 1.39 M nodes, five shots 1 km apart, 201 receivers, 3.5 s at a pinned
  $\Delta t = 2$ ms (1,750 steps), the operator pinned (`optimised_for = 4.2`) and the layers
  too (`pml_vp_max`), so observed and synthetic traces share their numerics. The gradient is
  of the smoothed model ($v_p$, $v_s$; true $\rho$) against the true model's traces, all five
  shots fired at once.

`just marmousi2-bench` (`examples/marmousi2/bench.sh`) runs each on the dev container's four
cores (Apple M1 Max, Colima VM with 8 GiB, nothing else running) and measures wall time and
peak resident memory; the forward speed is the run's own count of interior nodes times steps
over stepping time. It ran twice, a reboot apart; the ranges below are the two runs, and the
memory was the same both times.

**Results.**

| Run | Nodes | Steps | Stepping | Wall | Peak memory | Planned |
|---|---|---|---|---|---|---|
| Forward, whole section, order 8 | 3.94 M | 3,034 | 178–181 Mcell/s | 67.8–69.1 s | 425 MiB | 471 MiB |
| Forward, whole section, order 4 | 3.94 M | 2,694 | 258–271 Mcell/s | 40.5–42.6 s | 420 MiB | 433 MiB |
| Forward, central, order 8 | 1.39 M | 1,750 | 160–165 Mcell/s | 15.1–15.5 s | 162 MiB | 169 MiB |
| Forward, central, order 4 | 1.39 M | 1,750 | 211–222 Mcell/s | 11.8–12.4 s | 159 MiB | 155 MiB |
| Gradient, central, order 8 | 1.39 M | 1,750 | | 68.0–71.0 s | 4,396 MiB | 4,372 MiB |
| Gradient, central, order 4 | 1.39 M | 1,750 | | 51.7–53.2 s | 3,998 MiB | 3,953 MiB |

- *Speed.* Order 8 steps the whole section at 178–181 Mcell/s and order 4 at 258–271, a
  little above the 160 and 220 the examples gave (README, "Accuracy"). The central part is
  slower, 160–165 and 211–222, because a larger share of it is in the absorbing layers (along $x$, 32 of 241
  nodes against 32 of 681; along $y$, 32 of 41 in both), whose memory variables cost extra.
  Order 4 at $C = 0.85$ takes a longer step (2.23 ms against 1.98), so its run takes 60–62%
  of order 8's time. It is no bargain here: at 4.0 points per shortest wavelength it is outside
  its range (it warns below 5). For order 8's accuracy it needs 8–9 points per wavelength
  where order 8 needs 4.5 (2026-11-20, "Stencil order and time step"): on this 2.5D grid
  $(8.5/4.5)^3 \approx 6.7$ times the node-steps ($n_x$, $n_z$ and the steps grow as $1/h$;
  $n_y$ is set by the layers), which its 1.45–1.5 times faster stepping and 1.13 times longer
  step bring to about 4 times order 8's time.
- *Gradients* take 4.5–4.6 (order 8) and 4.3–4.4 (order 4) forward runs: the forward run, its replay
  from checkpoints, the adjoint run and the strains' correlations, besides setting up two
  solvers and writing the results. They hold 26 times the forward run's memory (checkpoints
  2.23 GiB and the replayed segment 1.50 GiB of the 4.27; 24 checkpoints, every 73 steps). The
  whole section's gradient would need 15 GiB (32 checkpoints, 8.05 GiB; a segment of 97
  steps, 5.45 GiB) against the forward run's 471 MiB, more than the container's 7 GiB:
  `gradient` refuses it before allocating anything.
- *Memory.* Forward runs peak between 10% below the plan and 3% (4 MiB) above it, gradient
  runs 24 and 45 MiB (0.5% and 1.1%) above. The plan counts what is allocated, so pages never
  written (zeros the allocator hands out lazily) aren't resident; the process itself (the
  program, threads' stacks, the allocator's arenas) is, and the plan leaves it out as "a few
  MiB".
- *The wavefield* (`marmousi2.toml`, with snapshots: 72 s) stays bounded: its peak falls
  from $2.9 \times 10^{-2}$ m/s near the source at 0.8 s to $1.3 \times 10^{-5}$ at 6 s. Most
  of the energy left at the end travels as Scholte and S waves in the slow sediments under
  the seafloor.
- *The gradient* sits almost entirely in the top few hundred metres below the seafloor, where
  those waves travel, and for $v_s$ its inner product with the true update,
  $-\nabla_{v_s} \cdot (v_{s,\text{true}} - v_{s,\text{start}})$, is negative.

**Inversion.** Five L-BFGS iterations of `central.toml` ($v_p$ and $v_s$, 395,010 unknowns;
the water and the seafloor's first solid row fixed; illumination preconditioning and
balanced components, the defaults) against `central-true.toml`'s five shots took 40 minutes:
340 s for the starting misfit and gradient, then about 410 s an iteration, five gradient runs
and five forward runs for the line search, which took the full step at the first try every
time. `memory --invert` planned 4.41 GiB. The misfit fell to 97.2%, 83.7%, 73.3%, 64.0% and
57.9% of the start's, but the model came no nearer the truth: over the free nodes
$\lVert m - m_\text{true} \rVert / \lVert m_0 - m_\text{true} \rVert$ is 1.001 for $v_p$ and
1.004 for $v_s$, whose start is off by 11.5% and 21.4% (RMS). By depth, with the cosine
between the update $m - m_0$ and the true one $m_\text{true} - m_0$:

| Depth | $v_p$ error ratio | cosine | $v_s$ error ratio | cosine |
|---|---|---|---|---|
| 500–750 m | 0.900 | +0.44 | 0.973 | +0.25 |
| 750–1000 m | 0.994 | +0.12 | 1.095 | −0.15 |
| 1000–1500 m | 1.020 | −0.15 | 1.070 | −0.17 |
| 1500–2500 m | 1.001 | −0.02 | 0.999 | +0.04 |
| 2500–3100 m | 1.003 | −0.22 | 0.999 | +0.14 |

The updates are a few percent (RMS 2.5% for $v_p$ and 5.6% for $v_s$ in the top band, under
1% below 1.5 km) and blotchy where the true update is layered. In the first 250 m under the
seafloor both speeds move towards the truth; from 750 m to 1.5 km $v_s$ moves away from it;
deeper, nothing moves. That is what one would expect if the slow waves along the seafloor
dominate the misfit and the smoothed start predicts them a cycle or more off: taking 300 and
400 m/s for the sediments' S speed in the true and the smoothed model, their traveltime over
2.5 km differs by about 2 s, at a period of 0.8 s. How much of the misfit they carry was not
measured.

**Decision.** Marmousi2 is now one command away (`just marmousi2`), and
`vibe-rations extrude` makes model files from any 2D SEG-Y section. These numbers are the
benchmark to compare later speed-ups with. Gradients of the whole section need the forward
field rebuilt rather than stored (`feat/reconstruct`, not merged yet) or a machine with
16 GiB; the central part is the gradient and inversion test case meanwhile. The sediments'
303 m/s S waves set the grid, and the seafloor is where the cost and the gradient are. Plain
L2 on all the data from a smoothed start lowers the misfit without recovering the model, so
an inversion of Marmousi2's elastic data needs more: lower frequencies first
[@fichtner2013], the slow seafloor waves muted or down-weighted, a better shallow start, or
a misfit less prone to cycle skipping [@engquist2014; @yang2016] (`feat/misfits`).
`central.toml` is the test case for them, these five iterations the baseline.

**Reproduce:** `just in marmousi2`, then `just in marmousi2-bench` (about 10 minutes);
`just in run examples/marmousi2/marmousi2.toml` for the example with snapshots. The
inversion: `just in run examples/marmousi2/central-true.toml --shots --out
out/marmousi2/observed`, then `just in invert examples/marmousi2/central.toml` (40 minutes,
writing the model after each iteration to `out/marmousi2/central`).

## 2026-11-21 — Marmousi2 in the browser

**Question.** Can the web page run Marmousi2? The browser solver refused model files (a
browser has no files to read), and the whole section at 25 m, 41 nodes wide in $y$, is
47 MB of `.npy` to fetch.

**Method.** Model files may now be sections, shape (nz, nx), which a run repeats along $y$
(`Model::from_files`); `extrude` writes sections unless given `--ny`. `model::Files` says
where a configuration's files come from, a directory or bytes in memory, and `Model`,
`memory` and `Simulation` read either. In the page, the worker asks which files a
configuration names (`Simulator.modelFiles`), fetches them relative to `examples/`, and passes
them in (`ModelFiles`); served by `vibe-rations serve`, the page uploads them with the
configuration instead. `web/build.sh` copies each example's files next to it and leaves out
an example whose files are missing; the web image's `marmousi2` stage fetches Marmousi2 and
converts it first, so the hosted page has them. The model files moved from `data/` to
git-ignored directories beside the configurations, since the run server refuses paths with
`..`. Two examples are for the browser: the whole section decimated to 50 m
(`marmousi2-50m.toml`: cells of $40 \times 40$ samples, $341 \times 41 \times 71$ = 0.99 M
nodes, $v_s$ from 321 m/s, so a Ricker at 0.6 Hz for 4.3 points per wavelength at 1.5 Hz,
8 s in 2,019 steps) and the central 6 km at 25 m (`marmousi2-central.toml`: 1.39 M nodes,
1.2 Hz, 4 s in 2,020 steps). Both ran in the Claude app's built-in browser pane (10 threads,
Apple M1 Max, the page visible), with live pictures, timed from Run to the last step.

**Results.**

| Example | Files | Memory | Speed | Run |
|---|---|---|---|---|
| `marmousi2-50m` | 290 KB | 129 MiB | 162–174 Mcell/s | 14.9 s |
| `marmousi2-central` | 408 KB | 170 MiB | 185 Mcell/s | 20.9 s |
| `marmousi2-50m`, page served by `vibe-rations serve` (4 threads in the container) | uploaded | | 152 Mcell/s | 14.5 s |

As sections the whole section at 25 m is 1.15 MB (47 MB as the grid). Set up from its
sections, `examples/marmousi2/marmousi2.toml` gets the same time step, steps and memory plan
(471 MiB) as from the grid, and a run from sections in memory records the same traces, bit for
bit, as the same model given as a layer (`runs_model_files_it_is_given`). Setting up takes
0.3–0.4 s in the page once the files are fetched. Decimating to 50 m halves the frequency the
grid holds (the sediments' $v_s$ barely changes: 321 m/s at 50 m, 303 at 25 m), so the
whole section at 50 m resolves less than the central part at 25 m, at about the same cost.

**Decision.** Examples may take model files, and the hosted page carries Marmousi2's two:
the decimated whole and the central part at full resolution. The whole section at 25 m
(3.94 M nodes, 3,034 steps) would fit a browser too (471 MiB), at these speeds in about a
minute on ten threads, but is left to the command line. The hosted page redistributes Marmousi2-derived
data; its licence (SEG wiki) could not be read, so that is the maintainer's call, made
2026-11-21.

**Reproduce:** `just in marmousi2`, `just in web`, `just serve`, then open
`http://localhost:8088/?example=marmousi2-50m` (or `marmousi2-central`) and Run; the page
shows speed and steps. `just in unit runs_model_files_it_is_given`.

## 2026-11-21 — Mutes and windows

**Question.** Marmousi2's five L2 iterations lowered the misfit without approaching the model,
and the suspects are the slow waves along the seafloor (2026-11-20, "Marmousi2"); the time
shifts need single arrivals, which whole traces with several don't give (2026-11-21, open
questions). Can the inversion weigh the data by time and offset, mute or window, and keep its
gradient exact?

**Literature** (two research agents; Crossref for every reference, the papers read where
noted in the bibliography). Mutes and time windows before the misfit are old practice: Shipp &
Singh muted multiples behind a line in offset and windowed by offset and time stage by stage
[@shipp2002]. For ocean-bottom data, Warner et al. muted the Scholte waves at short offsets
with "a simple mute because that allows us to be certain that we apply identical processes to
field and synthetic data" [@warner2013]; Prieux et al. saw Scholte waves slower than 350 m/s on
every component of Valhall's data and left them unmodelled [@prieux2013b]. Frequency-domain
codes, which can't window, damp: $e^{-(t - t_0)/\tau}$ from the first arrival
[@brenders2007; @brossier2009; @prieux2013], which "limits the offset range over which surface
waves are seen" [@brossier2009], relaxed over the stages to bring surface waves back in
[@romdhane2011]. Adjoint tomography windows single arrivals automatically, but on the
synthetic and anew for each model [@maggi2009], and leaves the processing's transpose
implicit [@liu2006]. A second agent read the code of DENISE Black-Edition (commit b07d81d;
the code of [@kohn2012]), SeisFlows (main d084262 and legacy 46a9604; [@modrak2018]), LASIF
2.0 (568071f; [@krischer2015]) and Pyadjoint (907cd8a, which has no paper): all filter first
and window afterwards, and none carries both through the adjoint source exactly (DENISE
injects $W r$ where its misfit has $W^2$).

**Method** (`[[inversion.mute]]`, `[[inversion.stage.mute]]`, `misfit::mute_weight`). A mute
keeps one side of a line $t = x / v + t_0$, $x$ the horizontal offset between the nodes source
and receiver snap to, multiplies the other by a floor $w$ (0 removes it), and goes from one to
the other over a taper, a half cosine or an exponential $e^{-\tau/T}$ (the time damping above,
from a line instead of picks); several multiply, and two make a window. The weights $M$ are
fixed by the survey and applied to observed and synthetic traces alike, so the misfit is
$\chi(F M s, F M d)$ with $F$ the stage's filter, and its adjoint source is
$M F\, \partial\chi$. Horizontal offset, not slant distance: in the second agent's runs of this solver the slow
seafloor wave from a source 400 m deep and one 50 m deep lined up in horizontal offset (lags
$-0.11$ to $+0.09$ s at 0.2–1.3 km) where slant distance predicted 0.22–0.85 s. Weights and
quiet traces are judged on the muted traces, from the observed ones only.

*The order.* As first written, the mute came after the filter, as in every code read. The
second agent's prototype found that order leaks, and `mute_study` reproduces it here: a 1.2 Hz
Ricker the synthetic has 0.100 s late, a ten times stronger arrival behind it (the synthetic's
0.3 s early), a mute keeping what comes before 4 s with a 0.42 s taper, a 0.6 Hz corner.
Filtering first, the zero-phase filter smears the strong arrival back across the line, where
the mute can't reach it, and the shift came out $-0.095$ s; muting first, $+0.100$ s. The
prototype also found filtering first lengthens the kept arrival into the taper, which clips it
(4.9% against $5 \times 10^{-10}$). So mutes come first (`the_mute_comes_before_the_filter`).
Even so the line and its taper must end before the muted arrival begins: with a whole period's
taper (0.83 s) the strong arrival's front reached into it, and the shift came out $+0.129$ s
muting first ($-0.153$ filtering first).

**Results.**

- $\langle F M x, y\rangle = \langle x, M F y\rangle$ to $10^{-12}$
  (`muting_then_filtering_is_transposed_exactly`); a change of the synthetic where it is muted, with or
  without a filter, changes neither the misfit nor the adjoint source, and where it is
  down-weighted to $w$ it counts $w^2$ (`muted_samples_do_not_count`); every misfit, filtered
  or not, with a window, agrees with finite differences through `Misfit`.
- Through the inversion (`the_mutes_keep_the_gradient_exact`, the model of the earlier checks):
  filtered L2 with the late arrivals down-weighted to a fifth, $8.7 \times 10^{-7}$ at
  $h = 10^{-2}$ (with the mute after the filter it had been $1.5 \times 10^{-5}$, and up to
  $3.8 \times 10^{-5}$ at smaller steps); time shifts in a window, $2.6 \times 10^{-5}$ at
  $h = 10^{-3}$, their round-off floor. The window changes what they measure: the change
  predicted down the gradient is $-9.47 \times 10^{-5}$ against $-1.41 \times 10^{-4}$ without.
- The memory plan counts the muted copy of the observed traces and the offsets
  (`Misfit::held_bytes`), and still equals what an inversion holds with a mute.

**Decision.** Mutes before the filter, fixed by the configuration, horizontal offsets. The taper
at least half the dominant period: for an arrival on the line, a shift of 0.100 s came out
0.048, 0.068, 0.088, 0.097 and 0.099 s with tapers of 0.01, 0.25, 0.5, 1 and 2 periods as
recorded, and 0.465, 0.221, 0.132, 0.130 and 0.113 s with the 0.6 Hz filter after the mute,
which lengthens the clipped arrival: on filtered traces a mute should not cut through an arrival
whose shift is wanted (`mute_study`). `delay` includes the wavelet's own (1.2/$f_0$ for a
Ricker). For the time shifts a floor must stay far below the kept arrival's amplitude over the
muted one's: behind a kept arrival, a ten times stronger one the same in both traces pulled the
0.100 s shift to 0.099, 0.093, 0.050 and 0.009 s at floors of 0.01, 0.03, 0.1 and 0.3. Windows picked on
the synthetics, as FLEXWIN's, would make the misfit jump between models; picked on the
observed data, or fixed per stage, they would not, which is what these mutes are. Next: the
local session tries them on Marmousi2's seafloor waves.

**Review.** An adversarial review (three reviewers, a skeptic per finding; 12 confirmed, 9 of
them found twice, 3 refuted) found no error in a gradient. Fixed: the L2 misfit's weights on
the receivers' side of the illumination ignored the mutes, so a receiver muted entirely still
sent waves out at full weight; they are now its weight times the trace's mean squared mute
weight. The README's window had a floor of 0.1, which the numbers above warn against (now
0.01); four new entries in the bibliography had a stray comma; docs and the memory report
spoke of a filtered copy where it may be muted, and left out the offsets. It also found, on
main before this branch, that the page's Gradient failed for any configuration with stages:
the memory estimate, made first, refused the iterations the page writes beside them
(`InversionConfig::gradient_only` now serves the estimate and the view alike).

**Reproduce:** `just unit misfit`; `just physics mutes_keep`; the numbers above,
`cargo test --release --lib mute_study -- --ignored --nocapture`.

## 2026-11-21 — Where Marmousi2's misfit lies

**Question.** Five L2 iterations of `examples/marmousi2/central.toml` from its start smoothed
over 250 m lowered the misfit to 57.9% without bringing the model nearer the truth, and moved
$v_s$ away from it between 750 m and 1.5 km (2026-11-20, "Marmousi2"). The suspected cause is
the slow waves along the seafloor (Scholte and S waves in the sediments, about 300 m/s in the
truth, 400 m/s in the start), which the start predicts a cycle or more off. How much of the
misfit, and of its gradient, do they carry?

**Method** (`tests/marmousi2_study.rs`, `misfit_by_arrival`). For the start against the true
model's five shots, the balanced L2 misfit, as the inversion measures it
($\tfrac12 \sum_c w_c \sum (d - d_\text{obs})^2$, $w_c$ the inverse of shot and component $c$'s
observed energy), is summed sample by sample into bins: before and after lines
$t = |\Delta x| / v + t_0$ for $v$ = 300, 400 and 500 m/s, by component, and by horizontal
offset. $t_0 = 0.87$ s is the earliest a slow wave can leave the seafloor below the source: the
Ricker's peak ($1.2 / f_0 = 1.0$ s) less half a period, plus the 425 m of water at 1500 m/s. The
records are 3.5 s long, so such a wave stays in them out to an offset of 0.8 km at 300 m/s and
1.3 km at 500. The gradient is linear in the residuals, so the gradient of the residuals before
the 400 m/s line (the rest zeroed) and the whole gradient give the part after the line as the
difference, two gradient runs per shot. Each part is taken relative ($m_0\, \partial\chi /
\partial m$) at the free nodes, as the inversion's unknowns are, and reported by depth band:
its norm, and the cosine of its descent direction $-\nabla$ with the true update
$m_\text{true} / m_0 - 1$. These are the raw gradients; the inversion divides them by the
source-side illumination [@shin2001], which raises the deep ones.

**Results.** The misfit (0.660 in all; by shot from $x$ = 1 to 5 km 0.057, 0.090, 0.140,
0.255 and 0.119):

| Split | Share |
|---|---|
| before / after $t = \lvert\Delta x\rvert / 300 + 0.87$ | 61.5% / 38.5% |
| before / after $t = \lvert\Delta x\rvert / 400 + 0.87$ | 47.2% / 52.8% |
| before / after $t = \lvert\Delta x\rvert / 500 + 0.87$ | 35.1% / 64.9% |
| $v_x$ / $v_z$ | 76.7% / 23.3% |
| offset 0–1 / 1–2 / 2–3 / over 3 km | 81.4% / 13.5% / 3.8% / 1.3% |

The gradient split at 400 m/s, norm and cosine with the true update by depth band:

| | 500–750 m | 750–1000 m | 1000–1500 m | 1500–2500 m | 2500–3125 m |
|---|---|---|---|---|---|
| $v_p$, before the line | $6.8 \times 10^{-3}$, +0.47 | $1.4 \times 10^{-3}$, +0.28 | $4.8 \times 10^{-4}$, −0.05 | $8.2 \times 10^{-5}$, +0.04 | $4.6 \times 10^{-6}$, −0.14 |
| $v_p$, after it | $3.6 \times 10^{-3}$, −0.04 | $1.2 \times 10^{-3}$, −0.27 | $7.5 \times 10^{-4}$, −0.17 | $7.6 \times 10^{-5}$, +0.01 | $6.6 \times 10^{-6}$, −0.13 |
| $v_p$, all | $8.2 \times 10^{-3}$, +0.37 | $1.8 \times 10^{-3}$, +0.04 | $9.4 \times 10^{-4}$, −0.16 | $1.3 \times 10^{-4}$, +0.03 | $9.9 \times 10^{-6}$, −0.15 |
| $v_s$, before the line | $1.7 \times 10^{-2}$, +0.13 | $2.4 \times 10^{-3}$, +0.01 | $1.1 \times 10^{-3}$, −0.20 | $1.9 \times 10^{-4}$, +0.01 | $3.3 \times 10^{-6}$, +0.15 |
| $v_s$, after it | $1.8 \times 10^{-2}$, +0.03 | $8.6 \times 10^{-3}$, −0.12 | $5.2 \times 10^{-3}$, −0.17 | $4.4 \times 10^{-4}$, −0.05 | $4.1 \times 10^{-6}$, +0.16 |
| $v_s$, all | $3.0 \times 10^{-2}$, +0.10 | $9.5 \times 10^{-3}$, −0.10 | $5.7 \times 10^{-3}$, −0.19 | $5.5 \times 10^{-4}$, −0.03 | $6.5 \times 10^{-6}$, +0.18 |

- *Where the misfit is.* Four fifths of it is within 1 km of the source, three quarters on the
  horizontal component, and more than half after the 400 m/s line: slow, horizontally
  polarised waves at near offsets, which is what Scholte and S waves in soft sediments under a
  fluid look like. The balance weighs components, not offsets, so the loud near offsets
  dominate. More than a third comes even after the 300 m/s line: the truth's slowest
  arrivals, and the coda of reverberations in the sediments.
- *Where its gradient points.* Between 750 m and 1.5 km the residuals after the line carry
  3.5 and 4.8 times the $v_s$ gradient that the earlier ones do, and point away from the truth
  (cosines −0.12 and −0.17): the whole gradient there does too (−0.10, −0.19), which is where
  the baseline's $v_s$ moved away. The earlier residuals' gradient points towards the truth in
  the top kilometre, clearly for $v_p$ (+0.47, +0.28), weakly for $v_s$ (+0.13, +0.01), but
  away for $v_s$ at 1–1.5 km too (−0.20). Below 1.5 km both are 10–1000 times smaller and
  their cosines close to zero: these 1.2 Hz data from 5 km of ocean-bottom nodes barely see
  it.
- The total's cosines match the baseline inversion's by band in sign (2026-11-20 table: $v_s$
  +0.25, −0.15, −0.17 in the top three bands), so its first gradient already decided the
  outcome.

**Decision.** The slow waves are the obstacle, as suspected: they carry most of the misfit and
the wrong-signed part of the gradient below 750 m. Muting or down-weighting them by arrival
time (`feat/mutes`, not yet on main) is the direct remedy; meanwhile the tools on main are run
in the order this suggests (next entries): L2 again in the study's own measure, as the
reference; frequency stages from 0.6 Hz, whose 1.7 s period comes nearer the slow waves'
roughly 2 s of delay [@bunks1995; @fichtner2013]; then time shifts and the phase misfit
[@luo1991; @fichtner2008], whose single-arrival assumption the overlapping slow waves at
near offsets break, so they are expected to do least.

**Reproduce:** `just in marmousi2`, then `just in marmousi2-study misfit_by_arrival` (15
minutes on four cores: 78 s for the observed traces, 706 s for the gradients).

## 2026-11-21 — Marmousi2: frequency stages against L2

**Question.** The slow seafloor waves carry most of the start's misfit, and the part of its
gradient that points away from the truth between 750 m and 1.5 km (previous entry). Does
inverting the low frequencies first, whose longer periods come nearer those waves' delay
[@bunks1995; @fichtner2013], keep the inversion from following them?

**Method** (`tests/marmousi2_study.rs`; `l2`, `frequency_stages`). From the start smoothed over
250 m, five iterations each, $v_p$ and $v_s$, the water and the first solid row fixed,
illumination preconditioning and balanced components (`central.toml`'s settings):

- *L2*: the baseline, rerun in the study, which measures after each stage the balanced L2
  misfit of all the data (the same scale for every variant) and, by depth band, the error
  ratio $\lVert m - m_\text{true} \rVert / \lVert m_0 - m_\text{true} \rVert$ and the cosine of
  $m - m_0$ with $m_\text{true} - m_0$.
- *Stages*: L2 below 0.6 Hz (two iterations), below 1.0 Hz (two), then all (one); the source is
  a 1.2 Hz Ricker, whose spectrum peaks at 1.2 Hz and falls to half by about 0.6 and 2 Hz. The
  0.6 Hz corner's period, 1.7 s, is half the record; each stage starts afresh (L-BFGS forgets
  its pairs, its first step is capped at 2% again) at the cost of a gradient run per shot.

**Results.** Error ratio and cosine by depth band after five iterations:

| | L2 left | | 500–750 m | 750–1000 m | 1000–1500 m | 1500–2500 m | 2500–3125 m |
|---|---|---|---|---|---|---|---|
| L2 | 57.9% | $v_p$ | 0.900, +0.44 | 0.994, +0.12 | 1.020, −0.15 | 1.001, −0.02 | 1.003, −0.24 |
| | | $v_s$ | 0.973, +0.25 | 1.095, −0.15 | 1.070, −0.17 | 0.999, +0.04 | 0.999, +0.15 |
| 0.6 Hz (2) | 97.9% | $v_p$ | 0.992, +0.28 | 0.998, +0.07 | 1.015, −0.21 | 0.999, +0.05 | 1.001, −0.07 |
| | | $v_s$ | 0.970, +0.39 | 0.988, +0.16 | 1.029, −0.12 | 1.009, −0.07 | 1.000, +0.06 |
| then 1.0 Hz (2) | 94.9% | $v_p$ | 0.979, +0.35 | 1.000, +0.01 | 1.019, −0.21 | 0.999, +0.05 | 1.001, −0.09 |
| | | $v_s$ | 0.947, +0.43 | 0.985, +0.17 | 1.043, −0.13 | 1.010, −0.06 | 1.000, +0.06 |
| then all (1) | 93.2% | $v_p$ | 0.976, +0.39 | 1.001, +0.01 | 1.020, −0.22 | 0.999, +0.05 | 1.001, −0.09 |
| | | $v_s$ | 0.946, +0.43 | 0.987, +0.16 | 1.044, −0.14 | 1.010, −0.06 | 1.000, +0.07 |

Over every free node both end at 1.001 for $v_p$; $v_s$ at 1.004 (L2) and 1.007 (stages). The L2
rerun reproduces the baseline's misfits (97.2, 83.7, 73.3, 64.0, 57.9%) and its table to the
last digit (the deepest band now reaches 3125 m, the last free node, where it stopped at
3100). The stages took 3,467 s against L2's 2,634: two more stage starts.

- *Stages keep $v_s$ from turning away at 750–1000 m* (cosine +0.16 against −0.15, error ratio
  0.987 against 1.095) and improve it more in the top 250 m (0.946 against 0.973), but not
  below 1 km: 1000–1500 m still moves away, if less (1.044 against 1.070). Its cosine is the
  same from the first stage on, so the low frequencies already see it the wrong way.
- *They fit much less of the data* (93.2% of the L2 misfit against 57.9%) and do less for $v_p$
  at the top (0.976 against 0.900): of five iterations, four run on 0.6 and 1 Hz, where a 1.2
  Hz Ricker has little energy (after its first iteration the 0.6 Hz stage's own misfit, 1.86,
  is 2.8 times the full band's at the start, because the balance weights each stage's
  filtered data by its own, small, energy),
  and each stage's first step is held to 2%.
- Where the slow waves are more than about a period off (the first entry's 2 s at 0.8 s), even
  0.6 Hz isn't low enough: half a period at 0.6 Hz is 0.8 s, and the data have nothing much
  lower.

**Decision.** Frequency continuation helps where the slow waves are within reach of its lowest
band, and nowhere deeper; it costs fit and time. Five iterations aren't enough to tell where
it would end, but the cosine at 1000–1500 m stays negative through every stage, so more of
them would follow the wrong direction there too. The remedy has to take the slow waves out of
the misfit (`feat/mutes`), or start nearer them; the next entries try the time shifts and the
phase misfit, which this survey's overlapping arrivals at near offsets don't suit.

**Reproduce:** `just in marmousi2-study l2`, `just in marmousi2-study frequency_stages` (44
and 58 minutes on four cores); each stage's model in `out/marmousi2/study/<variant>/`.

## 2026-11-21 — Marmousi2: time shifts, then L2

**Question.** Cross-correlation time shifts measure a phase delay rather than a waveform
difference, so they don't cycle-skip as long as the shift is within the correlation's main
lobe [@luo1991]. Do three iterations of them, then two of L2, get further than L2 alone?

**Method** (`tests/marmousi2_study.rs`, `time_shifts_then_l2`). As in the previous entry, with
`misfit = "traveltime"` and `max_shift = 0.8` s (about a period of the 1.2 Hz Ricker) for three
iterations, then L2 for two.

**Results.**

| | L2 left | | 500–750 m | 750–1000 m | 1000–1500 m | 1500–2500 m | 2500–3125 m |
|---|---|---|---|---|---|---|---|
| time shifts (3) | 99.7% | $v_p$ | 0.997, +0.26 | 0.998, +0.37 | 1.000, −0.17 | 1.000, −0.07 | 1.000, −0.01 |
| | | $v_s$ | 1.001, −0.14 | 1.000, −0.03 | 1.000, +0.01 | 1.000, +0.10 | 1.000, +0.06 |
| then L2 (2) | 83.5% | $v_p$ | 0.958, +0.44 | 0.996, +0.10 | 1.005, −0.18 | 1.000, −0.05 | 1.001, −0.21 |
| | | $v_s$ | 0.990, +0.14 | 1.018, −0.13 | 1.016, −0.17 | 1.000, +0.02 | 1.000, +0.15 |

The time shifts' first iteration took the capped 2% step; in the next two their misfit fell
from 14.611 to 14.596 (0.1%), the line search cutting each step to 0.02 and 0.03 of the
L-BFGS direction after three trials (largest change 0.0% and 0.1%). The two L2 iterations
afterwards behave like L2's first two (83.5% here, 83.7% there, the same signs by band). 3,363 s
in all.

- *The time shifts barely move the model.* Their misfit, $\tfrac12 \sum \tau^2$ over 2,010
  traces, is 14.6, an RMS shift of 0.12 s, and its gradient predicts it only over tiny steps,
  most likely (not measured trace by trace) because each trace's shift is the parabola through its correlation's peak, continuous but with a
  derivative that jumps when the peak moves to another lag (2026-11-20, "Choosing the misfit,
  and frequency continuation"), and near offsets here hold several arrivals of similar size,
  the slow waves among them, whose correlation peaks trade places as the model changes.
- *Their direction is wrong where it matters for $v_s$*: −0.14 in the top band, where L2 and the
  stages go the right way.
- The shift assumes one arrival per trace [@luo1991]; at near offsets, where four fifths of the
  L2 misfit is, these traces have several. Windowing single arrivals (`feat/mutes`'s windows)
  is what the measure needs here.

**Decision.** Time shifts on whole traces don't help on this survey; don't spend iterations on
them before the windows exist. The prediction (2026-11-21, "Where Marmousi2's misfit lies") that
they would do least of the three holds so far; the phase misfit is next.

**Reproduce:** `just in marmousi2-study time_shifts_then_l2` (56 minutes on four cores).

## 2026-11-21 — Marmousi2: the phase misfit, then L2

**Question.** The time–frequency phase misfit compares the phases of the traces' local spectra,
weighted towards where the observed data have energy, and doesn't cycle-skip while each
phase difference is under $\pi$ [@fichtner2008]. Three iterations of it, then two of L2?

**Method** (`tests/marmousi2_study.rs`, `phase_then_l2`). As the previous entries, with
`misfit = "phase"` for three iterations (its scale $f_0$ is the sources' dominant frequency,
1.2 Hz, the stage being unfiltered), then L2 for two.

**Results.**

| | L2 left | | 500–750 m | 750–1000 m | 1000–1500 m | 1500–2500 m | 2500–3125 m |
|---|---|---|---|---|---|---|---|
| phase (3) | 95.4% | $v_p$ | 0.941, +0.50 | 0.949, +0.42 | 0.999, +0.05 | 1.002, −0.06 | 1.000, −0.00 |
| | | $v_s$ | 1.004, −0.06 | 1.005, −0.09 | 1.006, −0.27 | 1.000, +0.06 | 1.000, −0.02 |
| then L2 (2) | 80.0% | $v_p$ | 0.908, +0.53 | 0.951, +0.36 | 1.005, +0.00 | 1.002, −0.06 | 1.001, −0.07 |
| | | $v_s$ | 0.992, +0.12 | 1.026, −0.14 | 1.022, −0.21 | 0.999, +0.05 | 1.000, +0.03 |

Every step was the full one at the first trial; the phase misfit fell from 1.399 after its
first iteration to 1.214 after its third. 3,066 s in all.

- *The phase misfit moves $v_p$, and moves it right*: after three iterations it has the best
  $v_p$ of every variant below 750 m (0.949 at 750–1000 m, where L2 reaches 0.994 in five), and
  its cosines stay positive down to 1.5 km, where L2's and the stages' turn negative.
- *It leaves $v_s$ alone, slightly the wrong way* (−0.06 to −0.27), and two L2 iterations
  afterwards turn $v_s$ the way L2 does (−0.14, −0.21 below 750 m).
- A plausible reading, not tested: the phase weights each time–frequency cell by the observed
  data's energy there, but compares only phases, so the loud, slow, mispredicted arrivals at
  near offsets weigh less than in L2, and the P and early S arrivals, which see $v_p$, count
  for more.

**All four, after five iterations** ($v_p$ / $v_s$ error ratio):

| | L2 left | 500–750 m | 750–1000 m | 1000–1500 m |
|---|---|---|---|---|
| L2 | 57.9% | 0.900 / 0.973 | 0.994 / 1.095 | 1.020 / 1.070 |
| stages, 0.6 → 1 Hz → all | 93.2% | 0.976 / **0.946** | 1.001 / **0.987** | 1.020 / 1.044 |
| time shifts, then L2 | 83.5% | 0.958 / 0.990 | 0.996 / 1.018 | 1.005 / **1.016** |
| phase, then L2 | 80.0% | 0.908 / 0.992 | **0.951** / 1.026 | 1.005 / 1.022 |

None recovers $v_s$ below 1 km; each moves it away there, the time shifts least because they
hardly move anything. Below 1.5 km nothing changes by more than 0.3%.

**Decision.** For $v_s$, which the slow waves spoil, the stages are the best of the four; for
$v_p$, the phase misfit. The receivers' side of the illumination is tried with the stages
(next entry). The phase misfit's way with $v_p$ suggests running it before the stages, or
alone for $v_p$, once the slow waves are muted; that waits for `feat/mutes`.

**Reproduce:** `just in marmousi2-study phase_then_l2` (51 minutes on four cores).

## 2026-11-21 — Marmousi2: the stages with the receivers' side

**Question.** The frequency stages were the best of the four variants for $v_s$ (previous
entry). The receivers' side of the illumination (2026-11-21, "The receivers' side of the
illumination") divides the gradient by the energy that waves sent from the receivers leave at
each node too, which should even out how the 201 ocean-bottom nodes and the five sources weigh
the model. Does it help here?

**Method** (`tests/marmousi2_study.rs`, `stages_with_receivers_side`). The stages of the
"frequency stages" entry with `receiver_runs = 8`: at each stage's start, eight runs with every
receiver firing at once with random signs estimate the receivers' side, which multiplies the
source side in the preconditioner.

**Results.** Error ratio and cosine by depth band:

| | L2 left | | 500–750 m | 750–1000 m | 1000–1500 m | 1500–2500 m | 2500–3125 m |
|---|---|---|---|---|---|---|---|
| 0.6 Hz (2) | 98.6% | $v_p$ | 0.997, +0.28 | 0.999, +0.06 | 1.013, −0.20 | 1.000, −0.00 | 1.000, −0.00 |
| | | $v_s$ | 0.991, +0.19 | 0.992, +0.13 | 1.030, −0.14 | 1.003, −0.07 | 1.000, +0.07 |
| then 1.0 Hz (2) | 96.5% | $v_p$ | 0.986, +0.35 | 1.004, −0.05 | 1.025, −0.20 | 1.000, −0.00 | 1.000, −0.03 |
| | | $v_s$ | 0.981, +0.20 | 0.991, +0.16 | 1.072, −0.15 | 1.006, −0.07 | 1.000, +0.10 |
| then all (1) | 93.9% | $v_p$ | 0.984, +0.38 | 1.005, −0.05 | 1.027, −0.21 | 1.000, +0.00 | 1.000, −0.04 |
| | | $v_s$ | 0.979, +0.22 | 0.996, +0.14 | 1.076, −0.16 | 1.006, −0.07 | 1.000, +0.11 |

Against the stages alone after five iterations (93.2% of the L2 misfit; $v_p$ 0.976, 1.001,
1.020 and $v_s$ 0.946, 0.987, 1.044 in the top three bands), every band that moves is a little
worse. The receivers' side cost 142 s more to set up and 3,930 s in all, against 3,467.

- *It shifts the update downwards*, as it should: the top band, where the receivers sit and
  the source side alone is largest, moves less ($v_s$ 0.979 against 0.946), and 1000–1500 m
  more. But that band's gradient points away from the truth (2026-11-21, "Where Marmousi2's
  misfit lies"), so moving it more makes it worse ($v_s$ 1.076 against 1.044).
- The cosines are nearly the same with and without it (+0.22 against +0.43 at the top is the
  exception); a preconditioner rescales the update but can't turn a gradient that points the
  wrong way.

**Decision.** Not with this misfit: the receivers' side helps where the gradient is right
and the top would otherwise take every update; here the deeper gradient is wrong, because of
the slow waves. Retry it once they are muted (`feat/mutes`). Of the tools on main, the stages
remain the best for $v_s$ and the phase misfit for $v_p$, and neither recovers the model below
1 km in five iterations.

**Reproduce:** `just in marmousi2-study stages_with_receivers_side` (65 minutes on four cores).

## 2026-11-21 — Marmousi2: L2 with the slow seafloor waves muted

**Question.** The slow waves along the seafloor carry most of the start's misfit and turn its
gradient away from the truth below 750 m (2026-11-21, "Where Marmousi2's misfit lies"). Muting
them, as Warner et al. did with the Scholte waves of their ocean-bottom data [@warner2013]
and as Prieux et al. left them out of theirs [@prieux2013b], should leave a misfit whose
gradient points the right way. Does it, with L2 alone?

**Method** (`tests/marmousi2_study.rs`, `muted_l2`; `[[inversion.mute]]`, 2026-11-21, "Mutes
and windows"). As the earlier variants, five L2 iterations, with one mute in every stage: keep
what arrives before $t = |\Delta x| / 600\ \text{m/s} + 0.87$ s, over a half-cosine taper of
0.42 s (half the dominant period, the least the mutes' entry recommends) to zero. $t_0$ is the
earliest a slow wave can leave the seafloor below the source (the first entry); at 600 m/s the
weight is zero by the time the start's 400 m/s waves arrive at 0.5 km (2.12 s) and 0.4 s before
them at 1 km. The mute applies to observed and synthetic traces alike. The study still reports
the L2 misfit of *all* the data, muted or not.

**Results.**

| | L2 left | | 500–750 m | 750–1000 m | 1000–1500 m | 1500–2500 m | 2500–3125 m |
|---|---|---|---|---|---|---|---|
| L2 | 57.9% | $v_p$ | 0.900, +0.44 | 0.994, +0.12 | 1.020, −0.15 | 1.001, −0.02 | 1.003, −0.24 |
| | | $v_s$ | 0.973, +0.25 | 1.095, −0.15 | 1.070, −0.17 | 0.999, +0.04 | 0.999, +0.15 |
| L2, muted | 85.7% | $v_p$ | 0.916, +0.44 | 0.971, +0.35 | 0.999, +0.05 | 1.000, −0.03 | 1.001, −0.29 |
| | | $v_s$ | 0.970, +0.26 | 1.004, −0.05 | 1.001, −0.05 | 1.000, +0.04 | 1.000, +0.13 |

Over every free node: $v_p$ 0.999, $v_s$ 1.000 (L2: 1.001, 1.004). The muted misfit fell from
0.332 after the first iteration to 0.240 after the fifth, every step the full one at the first
trial. 2,769 s in all, as L2 (the mute costs nothing measurable).

- *Muting stops the damage.* $v_s$ below 750 m no longer moves away (1.004 and 1.001 against
  1.095 and 1.070); its cosines there are about zero (−0.05), so it barely moves at all, which
  says the waves left (P, and S waves faster than 600 m/s) see little $v_s$ at these depths.
- *$v_p$ gains.* At 750–1000 m it moves towards the truth as the phase misfit's did (0.971,
  cosine +0.35; the phase misfit's 0.951, L2's 0.994), and at 1000–1500 m it no longer turns away (0.999, +0.05, against 1.020, −0.15):
  the first variant whose $v_p$ is nearer the truth over every free node. The top band is a
  little behind L2's (0.916 against 0.900); the deepest band's cosine is negative in every
  variant (−0.24 to −0.29), a change of 0.1%.
- It fits less of the data (85.7% of the L2 misfit against 57.9%), as it should: the slow
  waves are still mispredicted, and no longer what the inversion fits.

**Decision.** Mute the slow seafloor waves: with them gone, L2 improves $v_p$ down to 1 km and
leaves $v_s$ unharmed. $v_s$ below 750 m needs data that see it, so the next runs try the mute
with the frequency stages (whose low frequencies helped $v_s$ at 750–1000 m before), and
down-weighting instead of removing.

**Reproduce:** `just in marmousi2-study muted_l2` (46 minutes on four cores).

## 2026-11-21 — Marmousi2: the stages with the slow waves muted

**Question.** Muting the slow seafloor waves kept L2 from turning $v_s$ away (previous entry),
and the frequency stages were the variant that improved $v_s$ at 750–1000 m (2026-11-21,
"Marmousi2: frequency stages against L2"). Together?

**Method** (`tests/marmousi2_study.rs`, `muted_stages`). The stages (L2 below 0.6 Hz for two
iterations, below 1 Hz for two, then all for one) with the previous entry's mute in every
stage. The mute comes before each stage's filter (2026-11-21, "Mutes and windows").

**Results.**

| | L2 left | | 500–750 m | 750–1000 m | 1000–1500 m | 1500–2500 m | 2500–3125 m |
|---|---|---|---|---|---|---|---|
| 0.6 Hz (2) | 94.0% | $v_p$ | 0.990, +0.23 | 0.999, +0.08 | 1.001, −0.15 | 1.000, +0.01 | 1.000, −0.12 |
| | | $v_s$ | 0.981, +0.25 | 1.005, −0.08 | 1.002, −0.10 | 1.000, −0.01 | 1.000, +0.10 |
| then 1.0 Hz (2) | 90.4% | $v_p$ | 0.980, +0.27 | 0.998, +0.12 | 1.003, −0.19 | 1.000, −0.01 | 1.000, −0.12 |
| | | $v_s$ | 0.967, +0.28 | 1.008, −0.09 | 1.004, −0.10 | 1.000, −0.01 | 1.000, +0.08 |
| then all (1) | 89.2% | $v_p$ | 0.972, +0.32 | 0.995, +0.19 | 1.003, −0.17 | 1.000, −0.02 | 1.000, −0.14 |
| | | $v_s$ | 0.966, +0.28 | 1.009, −0.09 | 1.004, −0.10 | 1.000, −0.01 | 1.000, +0.09 |

3,487 s in all, as the stages unmuted (3,467).

- *Behind muted L2 everywhere but $v_s$'s top band* ($v_p$ 0.972, 0.995 and 1.003 against
  0.916, 0.971 and 0.999; $v_s$ 0.966 against 0.970), for the same reasons as unmuted: four of
  the five iterations on frequencies where the source has little energy, each stage starting
  again at a 2% step.
- *The stages lose what made them worth having.* Unmuted, the 0.6 Hz stage turned $v_s$ the
  right way at 750–1000 m (0.988, +0.16); muted, it doesn't (1.005, −0.08). So it was the slow
  waves themselves, at frequencies low enough to be within half a cycle of the truth, that
  carried the right $v_s$ there, and the mute takes them out with the rest.

**Decision.** Muting and the stages don't add up as configured. What this points to is a
mute that depends on the stage: keep the slow waves at the lowest frequencies, where they are
within half a cycle and see $v_s$, and mute them above (`[[inversion.stage.mute]]`). That is a
follow-up variant, not run yet. Of the runs so far, muted L2 is the best for $v_p$ and leaves
$v_s$ unharmed; the unmuted stages are the best for $v_s$ above 1 km.

**Reproduce:** `just in marmousi2-study muted_stages` (58 minutes on four cores).

## 2026-11-21 — Marmousi2: down-weighting the slow waves instead

**Question.** Removing the slow seafloor waves helps L2 and takes the stages' low-frequency
$v_s$ information with them (previous two entries). Does keeping a little of them, a floor of
0.1 instead of 0, keep some of that information without the damage?

**Method** (`tests/marmousi2_study.rs`, `down_weighted_l2`, `down_weighted_stages`). The
previous entries' mute with `weight = 0.1`: past the line and its taper the traces are
multiplied by 0.1 rather than zeroed, which the L2 misfit counts as $0.1^2 = 1\%$
(2026-11-21, "Mutes and windows"). L2 for five iterations, and the stages.

**Results.** After five iterations ($v_p$ / $v_s$ error ratio; L2 misfit of all the data):

| | L2 left | 500–750 m | 750–1000 m | 1000–1500 m | 1500–2500 m |
|---|---|---|---|---|---|
| L2, muted | 85.7% | 0.916 / 0.970 | 0.971 / 1.004 | 0.999 / 1.001 | 1.000 / 1.000 |
| L2, down-weighted | 84.9% | 0.915 / 0.970 | 0.971 / 1.004 | 0.999 / 1.002 | 1.000 / 1.000 |
| stages, muted | 89.2% | 0.972 / 0.966 | 0.995 / 1.009 | 1.003 / 1.004 | 1.000 / 1.000 |
| stages, down-weighted | 89.3% | 0.972 / 0.963 | 0.995 / 1.009 | 1.004 / 1.005 | 1.000 / 1.000 |

The cosines differ by 0.03 at most. 2,709 and 3,507 s.

- *A floor of 0.1 is as good as removing them*, for better and worse: the slow waves' 1% share
  of the misfit is too small to steer it, and too small to carry the $v_s$ information the
  unmuted 0.6 Hz stage had. Somewhere between 0.1 and 1 they start to dominate again; that
  isn't measured.

**Decision.** For L2 on this survey, mute (or a floor this low; it makes no difference). The
$v_s$ the slow waves carry at low frequencies needs them in at full weight where they are
within half a cycle and out where they aren't: a mute per stage (`[[inversion.stage.mute]]`,
none in the 0.6 Hz stage), the follow-up the previous entry named.

**Reproduce:** `just in marmousi2-study down_weighted_l2`, `just in marmousi2-study
down_weighted_stages` (45 and 58 minutes on four cores).

## 2026-11-21 — The gradient against the exact kernels

**Question.** The finite-difference checks show that the gradient is that of the discrete
run, to its round-off (2026-11-20, the gradient entries). How far is that from the gradient
of the wave equation itself, and what sets the difference: the wavefields' dispersion, the
time step, the stencil, the residual the adjoint run starts from, or the grid's sampling of
the model? This is the roadmap's step 4 for dispersion, and two open questions: the accuracy
study repeated for gradients, and whether the time-dispersion transform, which takes the time
error out of the traces, takes it out of the gradient too.

**Method** (`tests/gradient_accuracy.rs`, `tests/common/stokes.rs`). In a homogeneous full
space ($v_p$ = 3000 m/s, $v_s$ = 1732 m/s, $\rho$ = 2500 kg/m³) a point force's displacement is
Stokes's solution [@aki2002, eq. 4.23]; `stokes` evaluates it for any time function and
differentiates it in space in closed form. For $\chi = \tfrac12 \int_0^T r(t)^2\, dt$, $r$ the
residual of one velocity component at the receiver, Born's approximation and reciprocity give
the kernels [@tromp2005; @liu2006; @fichtner2011]

$$
K_\rho = -\int \ddot{\mathbf u} \cdot \mathbf a \, d\tau, \qquad
K_\lambda = -\int (\nabla\cdot\mathbf a)(\nabla\cdot\mathbf u)\, d\tau, \qquad
K_\mu = -\int \partial_j a_i\, (\partial_j u_i + \partial_i u_j)\, d\tau,
$$

$\mathbf u$ the forward displacement and $\mathbf a(\tau)$ the velocity at time $T - \tau$ of a
force at the receiver whose time function is the residual reversed. Both are Stokes's
solution, so the kernels are time integrals of closed forms (trapezium rule every 1 ms; the
reversed residual is tabulated every 10 μs with its first three derivatives and its
integral, and interpolated by cubic Hermite polynomials). Checks (`the_reference_is_consistent`):
the velocity agrees with the physics tests' Stokes solution to $10^{-5}$, the spatial
derivatives with finite differences to $10^{-6}$, the table with the residual to $10^{-7}$.

A force of $10^9$ N (Ricker, 10 Hz) is recorded 720 m away: S, a force along $y$ with $v_y$
recorded along $x$ (broadside S, so the $v_s$ kernel); P, a force along $x$ with $v_x$ recorded
(P along its axis, so the $v_p$ kernel). 720 m is 4.2 dominant and 10.4 shortest S
wavelengths. The observed trace is Stokes's solution in a medium 2% faster, so the residual
is what an inversion sees: mostly a shift of 8 ms, a twelfth of a period. The grid reaches
360 m beyond source and receiver and to either side, inside 12 nodes of CPML; record 0.75 s,
after which the residual is below $10^{-6}$ of its peak.

The solver's gradient at a node, times the sample interval over the cell volume (the misfit
sums samples rather than integrating them, and a node stands for a cell), is compared in
three ways:

- *all*: against the exact kernel at the node;
- *fields*: against the exact correlations taken where the solver takes them and averaged as
  it averages the model (`Problem::grid_kernels`): $\rho$'s three components at the six
  velocity positions around the node, half each (arithmetic means); $\lambda$ and the normal
  part of $\mu$ at the node; $\mu$'s shear parts at the twelve shear-stress positions, a
  quarter each (harmonic means; Moczo et al.'s averaging [@moczo2002]). What differs here is
  the solver's wavefields: dispersion, the time step, the residual;
- *grid*: the second against the first, the staggered grid's sampling of the model, with exact
  fields.

On the horizontal plane through source and receiver, outside half a dominant S wavelength
(87 m) around both, where the kernels are singular: in the first Fresnel zone, where the
detour is under half a dominant wavelength [@spetzler2004], and on the whole plane; the
error is $\lVert g - k \rVert / \lVert k \rVert$ for the main parameter ($v_s$ for S, $v_p$ for
P) and the cosine. And along a smooth change, a Gaussian ($\sigma$ = 60 m, cut at $4\sigma$)
midway between them, relative ($v_s(1 + \delta)$): the derivative a line search uses, against
its exact value integrated on a 10 m grid (5 m changes it by $5 \times 10^{-6}$), split in the
same way (`Problem::grid_integral`).

**Results.** *S, the default scheme* (order 8, optimised coefficients, $C$ = 0.85, time error
removed), by grid spacing; points per shortest S wavelength, the trace's misfit and the
recorded residual's, then the gradient's errors for $v_s$:

| Points per wavelength | Trace | Residual | Fresnel: all | cosine | fields | grid | Plane: all | fields | Along the Gaussian: all | fields | grid |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 3.46 | 0.96% | 1.49% | 4.15% | +0.99918 | 0.26% | 4.19% | 5.74% | 0.27% | −0.46% | +0.19% | −0.65% |
| 3.85 | 0.30% | 0.47% | 3.97% | +0.99925 | 0.11% | 4.00% | 4.97% | 0.12% | −0.42% | +0.11% | −0.53% |
| 4.62 | 0.06% | 0.09% | 2.73% | +0.99965 | 0.05% | 2.74% | 3.56% | 0.05% | −0.33% | +0.04% | −0.37% |
| 5.77 | 0.01% | 0.02% | 1.69% | +0.99987 | 0.01% | 1.69% | 2.17% | 0.01% | −0.24% | −0.01% | −0.23% |
| 6.93 | 0.01% | 0.02% | 1.20% | +0.99993 | 0.01% | 1.20% | 1.54% | 0.01% | −0.16% | +0.00% | −0.16% |
| 7.70 | 0.01% | 0.02% | 0.96% | +0.99996 | 0.01% | 0.96% | 1.26% | 0.01% | −0.12% | +0.01% | −0.13% |
| 8.66 | 0.01% | 0.02% | 0.76% | +0.99997 | 0.02% | 0.77% | 1.00% | 0.02% | −0.09% | +0.02% | −0.10% |

The wavefields' error in the gradient (*fields*) follows the trace's down to a floor of
0.01–0.02% (round-off, the layers, the 1 ms quadrature), and is smaller: a quarter of the
trace's misfit at 3.5 points per wavelength and a sixth of the recorded residual's. That is
not general: with the time error kept or with order 4 (below) the gradient is about as wrong
as the residual. Why the default's coarse-grid error weighs so little in the kernel was not
traced. The grid's sampling falls as $h^2$ (5.4 times
from 20 to 8 m, against 6.25) and dominates from the start. For $\rho$, averaged along all
three axes, the whole error is larger (11% at 3.5 points per wavelength, 1.8% at 8.7) and
*fields* about the same (0.18% and 0.02%). Cosines with the exact kernel: 0.9992 at 3.5, 0.99997 at 8.7.

*P*, the same: *fields* 0.12% (3.5 points per S wavelength, 6.0 per P), 0.12%, 0.04% (4.6),
then 0.01%; *grid* 0.00% at every spacing, because $\lambda$, which alone carries the $v_p$
kernel, sits on the node and isn't averaged; along the Gaussian −0.10% at 3.5, −0.04% at 4.6,
under 0.01% from 5.8. For $\rho$, *fields* is 0.35% at 3.5 and 0.10–0.15% from 4.6, and the
whole error falls from 6.1% to 1.2%.

*The time step* (S, 4.6 points per wavelength, order 8, optimised):

| Courant, time error | Trace | Residual | Fresnel: all | cosine | fields | grid | Plane: all | fields | Along the Gaussian: all | fields | grid |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.85, kept | 4.91% | 7.57% | 9.36% | +0.99952 | 8.39% | 2.74% | 9.75% | 8.39% | −7.65% | −7.31% | −0.37% |
| 0.5, kept | 1.69% | 2.61% | 4.40% | +0.99965 | 2.89% | 2.74% | 5.04% | 2.89% | −2.85% | −2.49% | −0.37% |
| 0.3, kept | 0.60% | 0.93% | 3.13% | +0.99966 | 1.02% | 2.74% | 3.92% | 1.02% | −1.24% | −0.88% | −0.37% |
| 0.85, removed | 0.06% | 0.09% | 2.73% | +0.99965 | 0.05% | 2.74% | 3.56% | 0.05% | −0.33% | +0.04% | −0.37% |
| 0.5, removed | 0.06% | 0.09% | 2.73% | +0.99965 | 0.04% | 2.74% | 3.56% | 0.05% | −0.33% | +0.04% | −0.37% |
| 0.3, removed | 0.06% | 0.09% | 2.73% | +0.99965 | 0.04% | 2.74% | 3.57% | 0.04% | −0.33% | +0.03% | −0.37% |

With leapfrog's time error kept, the gradient is wrong by 1.7 times the trace's misfit,
about as much as the recorded residual (8.4% against 7.6% at $C$ = 0.85), and $C$ = 0.3
still leaves 1.0%; removed, the error is 0.04–0.05% at every $C$. The transform
takes the time error out of the gradient as out of the traces, as expected: the transformed
traces are the time-continuous ones for any model, and the adjoint source goes back through
the transform's exact transpose (2026-11-20, "The adjoint and gradients with the new
defaults").

*Orders and coefficients* (S, 4.6 points per wavelength):

| Order, coefficients | Trace | Residual | Fresnel: all | cosine | fields | grid | Plane: all | fields | Along the Gaussian: all | fields | grid |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 4, Taylor | 11.68% | 18.00% | 15.55% | +0.99880 | 16.08% | 2.74% | 15.56% | 16.05% | +12.77% | +13.19% | −0.37% |
| 8, Taylor | 0.48% | 0.75% | 2.67% | +0.99965 | 0.43% | 2.74% | 3.52% | 0.42% | −0.03% | +0.34% | −0.37% |
| 12, Taylor | 0.05% | 0.08% | 2.73% | +0.99965 | 0.03% | 2.74% | 3.58% | 0.03% | −0.35% | +0.02% | −0.37% |
| 16, Taylor | 0.01% | 0.02% | 2.74% | +0.99965 | 0.01% | 2.74% | 3.58% | 0.02% | −0.36% | +0.01% | −0.37% |
| 4, optimised | 8.95% | 13.80% | 10.84% | +0.99912 | 11.24% | 2.74% | 10.96% | 11.23% | +8.66% | +9.06% | −0.37% |
| 8, optimised | 0.06% | 0.09% | 2.73% | +0.99965 | 0.05% | 2.74% | 3.56% | 0.05% | −0.33% | +0.04% | −0.37% |
| 12, optimised | 0.01% | 0.02% | 2.74% | +0.99965 | 0.01% | 2.74% | 3.58% | 0.02% | −0.36% | +0.00% | −0.37% |
| 16, optimised | 0.01% | 0.02% | 2.74% | +0.99965 | 0.02% | 2.74% | 3.58% | 0.02% | −0.36% | +0.00% | −0.37% |

Order 4 gets its gradient 11–16% wrong, 1.3–1.4 times its traces' misfit and a little
less than its residual's (14–18%), and along the Gaussian 9–13%; order 8 with Taylor coefficients 0.43%; from order 8 optimised up, the floor.

*The residual.* With the exact residual as the adjoint source instead of the recorded one,
*fields* falls from 0.26% to 0.12% at 3.5 points per wavelength (along the Gaussian from
+0.19% to −0.02%), and from 0.05% to 0.02% at 4.6: at the coarsest grid about half of the
gradient's error comes in with the residual, whose own error is the recorded trace's
(1.5 times the trace's misfit, the residual being about two thirds of the trace), and half
through the forward and adjoint wavefields.

**Decision.** No change to the defaults. At them, from 3.5 points per shortest S wavelength,
the solver's gradient is the wave equation's to within 0.26% (the cosine 0.9992), less than
its traces' error, and the time step doesn't matter; the gradient is as trustworthy as the
traces, so dispersion doesn't limit an inversion at the grids we use. What separates it from
the continuous kernel is how the staggered grid samples $\rho$ and $\mu$: a node's value stands
for the averages around it, which differs from a point sample by $O(h^2)$. That is the model's
parametrisation, not an error of the gradient: the forward runs use the same averages, and
differentiating the discrete run, as we do, gives the exact gradient for that
parametrisation [@sirkes1997]. Along smooth changes it costs 0.1–0.6%. A baseline test,
`the_gradient_matches_the_exact_kernels`, holds the coarsest grid's numbers (S: 0.26% and
+0.19%; P: 0.12% and −0.10%). Still to measure: interfaces and the free surface, which have
no closed-form kernels, against finer grids; and whether an inversion's model error follows.

**Reproduce:** `just gradient-accuracy` (38 minutes on four cores; `just gradient-accuracy
by_time_step` and the like for one table); the baseline with `cargo test --release --test
gradient_accuracy`.

## 2026-11-21 — A planar interface

**Question.** How accurately does the scheme reflect and transmit at an interface between
two solids, and does it matter where the interface falls between the nodes? The open
question asked to measure it against analytic coefficients and to try a band-limited model
[@mittet2017].

**Method.** Plane waves at normal incidence, where the exact answer is one-dimensional: the
particle velocity is reflected by $R = (Z_1 - Z_2)/(Z_1 + Z_2)$ and transmitted by $T = 1 + R$,
with $Z = \rho v_p$ for P and $\rho v_s$ for S, and nothing converts. A uniform plane of
sources (σzz for P, $x$ forces for S) sends them down; each side of it takes half the jump,
so a moment $m_0$ per node sends $v_z = m_0 w(t)/(2\rho v_p^2 h^2)$ and a force $f$ per node
$v_x = f w(t)/(2\rho v_s h^2)$. The grid is 25 nodes wide; the values outside its sides copy
the nearest node's every half step, so every $z$-plane stays exactly uniform and the side
layers, which act only on $x$ and $y$ derivatives, do nothing. (With ordinary edges, waves
from the sides spoilt the transmitted P wave: 5% low with the receivers 10 nodes from the
side layers, 1% at 30.) Upper medium $v_p = 2000$, $v_s = 1155$ m/s, $\rho = 2000$
kg/m³; lower 3000, 1500, 2500: $R = -0.304$ for P, $-0.238$ for S. $h = 10$ m, a Ricker
at $f_0 = 10$ Hz (P 20 and S 11.6 points per wavelength at $f_0$; S 4.6 at $2.5 f_0$), and at
$f_0 = 5$ Hz with the distances in nodes doubled; the default scheme (order 8, optimised,
$C = 0.85$, time error removed). The interface lies a fraction $\delta$ of a node below a
node. Four ways of putting the two layers on the nodes, each giving node values that the
solver then averages as always (ρ arithmetically onto the velocities, μ harmonically onto
the shear stresses):

- *node*: each node takes the layer it lies in, a node on the interface the one below
  (`[[model.layer]]` now);
- *cell*: each node takes its cell's average over $[z - h/2, z + h/2]$, density
  arithmetically and the moduli $\rho v_p^2$, $\rho v_s^2$ harmonically;
- *band*: the step band-limited to the grid's Nyquist wavenumber (a sinc tapered by a Hann
  window over ±4 nodes, integrated) and sampled, applied to density and the compliances
  $1/\rho v_p^2$, $1/\rho v_s^2$;
- *band moduli*: the same applied to the moduli themselves.

Each arrival is compared with the exact one over its window, the other arrivals' exact
values taken off: the misfit relative to the arrival, the least-squares amplitude ratio, and
the shift from the correlation's peak, in units of $h/v$ of the upper medium (a reflection's
shift is twice the interface's error in position).

**Results.** Reflected waves, $f_0 = 10 \to 5$ Hz (direct waves 0.010% and 0.026%):

| Layers | $\delta$ | P misfit | S misfit | Shift ($h/v$) |
|---|---|---|---|---|
| node | 0 | 35.1% → 17.6% | 57.0% → 29.9% | $-1.00$ |
| node | 0.25 | 17.6% → 8.8% | 30.2% → 15.2% | $+0.50$ |
| node | 0.5 | 0.72% → 0.17% | 7.5% → 1.8% | 0.00 |
| node | 0.75 | 17.6% → 8.8% | 29.4% → 15.1% | $-0.50$ |
| cell | 0, 0.25, 0.75 | 2.9–4.2% → 0.74–1.0% | 18–22% → 4.9–5.9% | $\le 0.08$ |
| band | any | 1.16–1.19% → 0.34–0.35% | 12.4–12.8% → 3.4% | $\le 0.02$ |
| band moduli | 0–0.75 | 3.8–9.0% → 1.9–4.3% | 11–16% → 4.0–6.1% | $-0.10$ to $-0.26$ |

Transmitted waves: node 5.9% (P) and 6.9% (S) at $\delta = 0$, 0.078% and 0.44% at 0.5; band
0.125% and 0.75% (→ 0.04% and 0.21%). Cell averages at $\delta = 0.5$ are the node values.

1. **Node-sampled layers put the interface midway between the last node above it and the
   first below**, wherever it lies: the shifts are $2(\frac12 - \delta)$ to within 0.01. The
   error is first order (it halves with $h/\lambda$). A layer's top on a node, as at both
   interfaces of `examples/layered.toml`, reflects $h/v$ early.
2. **Midway the error is second order** (P ×4.2, S ×4.1 per halving), an amplitude error:
   the S reflection is 6.3% weak at 11.6 points per wavelength at $f_0$, 1.5% at 23.
3. **The stencil's order doesn't matter** (`interfaces_by_order`): orders 8, 12 and 16
   give the same reflection misfits to 0.02%, at $\delta = 0$ and 0.5, node-sampled or
   band-limited, and order 4 differs by at most 0.06% for P and 1.8% for S, where its own
   dispersion shows (its direct S wave is 1.2% off). At interfaces the model's sampling
   sets the error, and longer stencils gain nothing.
4. **Averaging removes the position error but not the second-order one.** Cell averages
   leave shifts of at most 0.08 $h/v$ and larger amplitude errors than midway node sampling.
   Band-limiting makes the result independent of $\delta$ to 0.04% (P) and 0.4% (S), with
   errors between the two. Band-limiting the moduli rather than the compliances leaves a
   $\delta$-dependent position error of 0.05–0.13 $h$, as arithmetic averaging of moduli
   would: the compliances are what to band-limit.
5. **The taper** (`interfaces_by_band_limit_width`): ±2 nodes doubles the error (P 2.3–2.4%,
   S 17–18%); ±8 barely improves on ±4 (P 0.94–0.96%, S 12.5–12.6%).

Relative to the incident wave, band-limited reflections are off by 0.36% (P) and 3% (S) at
these samplings; the S reflection is the weak spot at the point densities FWI runs at.

**Decision.** `a_planar_interface_reflects_and_transmits_plane_waves` in `tests/physics.rs`
checks the interface midway (baselines above) and pins the half-node shift of a layer's top
on a node, so a change to how layers meet the nodes shows up there. Layers stay sampled at
the nodes for now: changing that changes every layered model and every baseline built on
them, and fluids need their own rule (a compliance band-limited into water makes $\mu$ zero
up to four nodes away; the Marmousi2 sections make a node fluid by majority instead).
Band-limited layers, density and compliances over ±4 nodes, are the candidate: they let an
interface lie anywhere at second-order cost. The gradient is with respect to node values,
so this changes only how a configuration's layers become a model, not the gradient. The
sinc and its taper here are ours, not checked against Mittet's paper.

**Reproduce:** `just interfaces` (3 minutes);
`just physics a_planar_interface_reflects_and_transmits_plane_waves`.

## 2026-11-21 — Band-limited layers, opt in

**Question.** The planar-interface entry found that band-limiting density and compliances
puts an interface where it lies. Can layers be built that way in general, including water
and large contrasts, and what would making it the default change?

**Method.** `[model] interfaces = "band-limited"` (`model::layer_profile`): each layer
weighs the difference of the steps at its top and bottom, and density and the compliances
$1/\rho v_p^2$, $1/\rho v_s^2$ are the weighted means. Two problems needed rules. (1) The
step's 7.3% overshoot, across a contrast $r$ in compliance, takes the stiffer side's
compliance $c$ to $c(1 - 0.073(r - 1))$: negative above $r = 14.7$ (a soft sediment on rock,
water on rock), and stiff enough to shorten $\Delta t$ well before. (2) A fluid's S
compliance is infinite, so averaging it would make $\mu$ zero four nodes into the solid.
Candidates for (1), in the study's plane-wave setup (`just interfaces`; 20 and 11.6 points
per wavelength at $f_0$, $\delta$ from 0 to 0.75): the step clamped to $[0, 1]$, and the
logarithms of the moduli band-limited instead. Then a large contrast, a sediment
($v_p = 1800$, $v_s = 400$ m/s, $\rho = 1900$) on rock (4500, 2600, 2600; 48 times the shear
modulus, $R = -0.55$ for P, $-0.80$ for S) at $f_0 = 3.5$ Hz, $\delta$ = 0, 0.25 and 0.5.
Last, `just physics` and the unit tests with the default flipped to band-limited.

**Results.**

- *Clamping* the step costs most of the gain: P reflection 2.5–6.6% (band-limited 1.2%), S
  15–24% (12.6%), no better than cell averages.
- *Logarithms* of the moduli stay positive but leave a position error, shifts of
  $-0.10$ to $-0.12$ $h/v$ for P (3.4–4.5%) and $-0.05$ to $-0.09$ for S: like moduli,
  only compliances put the interface where it lies.
- *At the large contrast* (reflection misfits; transmitted S in brackets):

  | Layers | P, $\delta$ = 0 / 0.25 / 0.5 | S, $\delta$ = 0 / 0.25 / 0.5 |
  |---|---|---|
  | node | 13.5 / 6.8 / 0.07% | 57.7 / 31.2 / 2.4% (25.7 / 14.1 / 3.3%) |
  | sinc, clamped so no modulus stiffens by more than 1.25 | 0.60 / 0.48 / 0.23% | 11.9 / 8.7 / 6.1% (11.5 / 8.5 / 6.8%) |
  | cell averages | 0.37 / 0.26 / 0.07% | 9.1 / 8.0 / 2.4% (9.3 / 7.5 / 3.3%) |

  Cell averages, which never overshoot and also keep the interface in place, beat the
  clamped sinc here; at the moderate contrast the sinc beats them (P 1.2% against
  2.9–4.2%, S 12.6% against 18–22%).
- *The rule*, then: an interface takes the sinc step unless its contrast would let the
  overshoot stiffen a modulus by more than 1.25 (`MAX_STIFFENING`: a compliance ratio above
  3.7, P and, between solids, S), and cell averages otherwise. Water: a node is fluid if the
  layer it lies in is, as with node sampling; a solid node takes the solids' S compliance
  alone, at most half its P modulus where a fluid's P compliance reaches it, so its bulk
  modulus stays positive ($\lambda \ge 0$). Through the configuration the plane-wave tests
  give the study's numbers exactly (`band_limited_interfaces_act_where_they_lie`: P 1.19%
  and 1.17%, S 12.6% and 12.8% with the interface on a node and a quarter of one off).
- *The overshoot speeds up some nodes*: 3022 m/s next to a 3000 m/s layer (2200 above),
  2179 m/s in 2000 m/s rock under water, so $\Delta t$ shrinks by 0.7–8% and
  `pml_vp_max` grows with it.
- *`just physics` with band-limited layers as the default*: the analytic tests, CPML,
  reciprocity, the dot-product test and the sphere inversion are unchanged; the long run
  drains to $2.5 \times 10^{-6}$ and $9.6 \times 10^{-7}$ (node-sampled $2.7 \times 10^{-6}$
  and $1.1 \times 10^{-6}$); the walls-and-water wavefield rebuilt backwards drifts
  $1.8 \times 10^{-6}$ ($9.6 \times 10^{-7}$). Three gradient checks fail, and in the unit
  tests the inversion's set-up test: it expects `pml_vp_max` to be the fastest layer's
  (2000, now 2179). The gradient checks perturb every node except those at the model's
  fastest speed; band-limited, the fastest is an overshoot plane (3022 m/s) and the 3000 m/s
  layer just below it is perturbed past it, which shortens $\Delta t$: the perturbed run has
  1280 samples against the observed 1270, and with water the $v_p$ check is off by 1%. Both
  models of a check have the same $\Delta t$ (checked: 253 steps at order 4, 285 at order 8,
  against 251 and 283 node-sampled), so the gradient isn't at fault; it is with respect to
  node values, which band-limiting only sets.

**Decision.** Band-limited layers are opt-in; node sampling stays the default. Making them
the default needs the gradient checks to spare nodes near the fastest speed (or pin
$\Delta t$ and `pml_vp_max`), the inversion test's expected `pml_vp_max`, and the baselines
above. Spheres are still node-sampled.

**Reproduce:** `just interfaces` (`interfaces_at_a_large_contrast` for the table);
`just physics band_limited_interfaces_act_where_they_lie`; `just unit model::`. For the
default flipped, move `#[default]` in `config::Interfaces` to `BandLimited` and run
`just physics` and `just unit`.

## 2026-11-21 — Marmousi2: what the mute removes

**Question.** Muting the slow seafloor waves stopped L2 from turning $v_s$ away below 750 m,
but left it barely moving there (2026-11-21, "Marmousi2: L2 with the slow seafloor waves
muted" and the two entries after it). How much of the start's misfit and of its gradient
does the mute take away, and does what it leaves between 1 and 1.5 km point towards the
truth? Is the wrong-way part the horizontal component's, which carries three quarters of the
misfit?

**Method** (`tests/marmousi2_study.rs`, `misfit_by_mute`). As "Where Marmousi2's misfit lies",
but split by the mute's weights $M$ instead of a hard line: for the start against the true
model's five shots, the part of the balanced L2 misfit a mute keeps,
$\tfrac12 \sum_c w_c \sum (M r)^2$ with $r = d - d_\text{obs}$, and its gradient (adjoint
source $w_c M^2 r$), each component alone. $w_c$ is the balance of the traces as recorded, so
the shares add up with the earlier entry's; the inversion balances by the muted observed
traces instead, and its gradient (a recombination of the two components') is reported too.
The part removed is the whole less the part kept. By depth band: the norm of the relative
gradient $m_0\, \partial\chi / \partial m$ and the cosine of $-\nabla$ with the true update,
raw and divided by the source-side illumination as the inversion's preconditioner divides
(water level $10^{-3}$). Six gradient runs a shot: the whole, the muted and the down-weighted
mute's kept parts, each component alone.

**Results.** The mute keeps 35.4% of the misfit ($v_x$ 29.3%, $v_z$ 6.1%), so it removes
64.6%; down-weighting to 0.1 keeps 36.5%. Norm, and cosine raw / preconditioned:

| | 500–750 m | 750–1000 m | 1000–1500 m | 1500–2500 m |
|---|---|---|---|---|
| $v_p$, all | $8.2 \times 10^{-3}$, +0.37 / +0.40 | $1.8 \times 10^{-3}$, +0.04 / +0.01 | $9.4 \times 10^{-4}$, −0.16 / −0.17 | $1.3 \times 10^{-4}$, +0.03 / −0.02 |
| $v_p$, kept | $6.6 \times 10^{-3}$, +0.45 / +0.46 | $1.3 \times 10^{-3}$, +0.31 / +0.30 | $3.9 \times 10^{-4}$, +0.02 / +0.01 | $6.5 \times 10^{-5}$, +0.00 / −0.03 |
| $v_p$, removed | $3.8 \times 10^{-3}$, +0.02 / +0.05 | $1.2 \times 10^{-3}$, −0.28 / −0.27 | $8.1 \times 10^{-4}$, −0.20 / −0.20 | $8.9 \times 10^{-5}$, +0.05 / −0.01 |
| $v_s$, all | $3.0 \times 10^{-2}$, +0.10 / +0.07 | $9.5 \times 10^{-3}$, −0.10 / −0.12 | $5.7 \times 10^{-3}$, −0.19 / −0.19 | $5.5 \times 10^{-4}$, −0.03 / −0.01 |
| $v_s$, kept | $1.5 \times 10^{-2}$, +0.20 / +0.17 | $1.5 \times 10^{-3}$, +0.08 / +0.04 | $5.4 \times 10^{-4}$, −0.10 / −0.07 | $6.6 \times 10^{-5}$, −0.01 / +0.01 |
| $v_s$, removed | $2.0 \times 10^{-2}$, −0.01 / −0.02 | $9.1 \times 10^{-3}$, −0.12 / −0.13 | $5.6 \times 10^{-3}$, −0.18 / −0.19 | $5.3 \times 10^{-4}$, −0.03 / −0.01 |

By component, raw cosines at 750–1000 m and 1000–1500 m:

| | $v_p$, $v_x$ | $v_p$, $v_z$ | $v_s$, $v_x$ | $v_s$, $v_z$ |
|---|---|---|---|---|
| all | −0.04, −0.21 | +0.29, +0.10 | −0.11, −0.19 | +0.02, −0.12 |
| kept | +0.24, −0.06 | +0.41, +0.19 | +0.08, −0.09 | +0.01, −0.15 |

With the inversion's own balance the kept gradient's cosines are within 0.04 of these ($v_p$
+0.44, +0.34, +0.05; $v_s$ +0.18, +0.07, −0.11 in the top three bands); its misfit at the start
is 0.338, which `muted_l2`'s first step took to 0.332. The down-weighted mute's numbers are
within 0.02 of the muted ones throughout. 3,968 s, with the parent session's CI sharing the
four cores for part of it (not a benchmark).

- *The mute removes what pointed the wrong way.* Its removed part points away from the truth
  in every band between 750 m and 1.5 km ($v_p$ −0.28 and −0.20, $v_s$ −0.12 and −0.18), as
  the split by arrival time found with a hard 400 m/s line.
- *It removes the slow waves' share of $v_s$, which is most of it below 750 m*: the kept
  $v_s$ gradient is 16% as large as the whole at 750–1000 m and 9% at 1–1.5 km, where the
  kept $v_p$ gradient is 73% and 41% as large. What is left of $v_s$ there is weakly right
  above 1 km (+0.08) and weakly wrong below (−0.10): why muted L2's $v_s$ barely moves.
- *Not only the horizontal component.* Unmuted, $v_x$ carried the wrong sign of $v_p$'s
  gradient (−0.21 at 1–1.5 km, where $v_z$'s is +0.10), and both carried $v_s$'s. Muted, the
  remaining wrong $v_s$ at 1–1.5 km is in both components ($v_x$ −0.09, $v_z$ −0.15), so it isn't
  the horizontal component, nor (muted) the slow waves; it lies in what arrives before the
  line, which this run doesn't split further (S waves and conversions faster than the line, or
  a timing error of the early arrivals from the smoothed start: not measured).
- *The preconditioner changes no cosine by more than 0.05*; it rescales, it doesn't turn.

**Decision.** The mute does what it was for: it takes out the part of the gradient that
pointed away, two thirds of the misfit with it, and leaves a gradient that points right for
$v_p$ down to 1 km. For $v_s$ below 750 m these data, muted, have little to say, and below 1 km
what they say is slightly wrong; the slow waves that carried $v_s$ there only help where they
are within half a cycle (the stages' 0.6 Hz, previous entries), which the next variant keeps.
*Tuned on this dataset:* the line (600 m/s from 0.87 s) was chosen on this part with its true
model's 300 m/s in view; the held-out entry below works it out from the start alone.

**Reproduce:** `just in marmousi2-study misfit_by_mute` (about 50 minutes on four cores
alone).

## 2026-11-21 — Marmousi2: the mute line without the truth, and a held-out part

**Question.** The maintainer's caution: "Be wary of overfitting for one dataset". The slow
waves' mute line, the stages' corners and the ranking of every variant were chosen on the
central 6 km and judged against its own true model. Which of them can be worked out from the
observed data, the start and the source instead, and what can serve as a held-out case?

**Method.** (1) *Measuring the slow waves' speed in the gathers*, observed and the start's
synthetic ones, within 1 km of each source (receivers 25 m apart; the 3.5 s records keep waves
at 300–400 m/s only to 0.8–1 km): the energy (balanced, smoothed over half a period) stacked
along lines $t = \tau + |\Delta x| / v$ for $v$ from 150 to 1000 m/s, slower than sound in
water; the same counting only after the direct wave's peak plus $0.6/f_0$ (a time the water
depth fixes); the f-k energy binned by speed $f/k$; and the loudest wavenumber of speeds
under 1000 m/s at each frequency from 0.3 to 3 Hz under a Hann taper over offset, with and
without the direct wave removed. (2) *A rule from the start and the source alone*
(`slow_waves_mute`, `mute_lines`): the line passes through $t_0$, the earliest a wave can leave
the seafloor below the source (the Ricker's peak $1.2/f_0$ less half a period, plus 425 m of
water at 1500 m/s: 0.87 s), with a taper of half a period ($0.5/f_0$ = 0.42 s), and is as fast
as needs be for the taper to end, from 500 m offset on, before waves at the start's slowest
S speed at the receivers: $1/v_\text{line} = 1/v_s - 0.42\ \text{s} / 500\ \text{m}$. (3) *A
held-out part*: the western 6 km, $x$ = 0–6 km (`convert.sh`'s `25m-west` and its start
smoothed over 250 m, `Case::part("west", …)`), with the same survey, start and numerics. Its
seafloor is the same (first solid node at 475 m, $v_s$ 345 m/s there in the truth); below it
the model differs ($v_p$ from 1028 to 4560 m/s, within the pinned time step's limit).

**Results.**

- *The measurements failed.* The stack picked 1000 m/s, the edge of the search, on every shot,
  observed and synthetic: between 300 m and 1 km the direct wave's hyperbola is close to a
  1000 m/s line. Counting only after the direct wave, the start's gathers still gave 1000 m/s
  and the observed ones 410–1000 m/s by shot; the f-k bins gave 955 m/s throughout and the
  f-k peaks 980–1000 m/s at most frequencies (the direct wave's leakage, falling with
  wavenumber, is largest at the first wavenumbers under 1000 m/s), and with the direct wave removed scattered from 320 to
  1000 m/s by shot and frequency. The direct wave and the early arrivals hold most of the
  recorded energy; the slow waves hold most of the *misfit* (they are mispredicted, the early
  arrivals aren't), not of the energy. The measuring code isn't kept.
- *The rule* gives 616 m/s on the central part (the start's slowest $v_s$ at the receivers,
  407 m/s), next to the 600 its runs used, which were chosen with the true model's 300 m/s and
  the start's 400 in view ("Where Marmousi2's misfit lies"); and 690 m/s on the western part
  (438 m/s).
- *The held-out runs were not done* before the branches were tied up: L2, muted L2 with the
  rule's 690 m/s line, and the per-stage mute, on the western part (`west_l2`,
  `west_muted_l2`, `west_stages_muted_after_the_first`). Open questions.

*What was tuned on which dataset:*

| Choice | From |
|---|---|
| the line's $t_0$ (0.87 s) and taper (0.42 s) | the source and the water's depth, not tuned |
| its speed, 600 m/s | the central part, with its true model in view; the rule gives 616 |
| full mute from 500 m offset, a half-cosine taper | by hand, on the central part |
| the stages' corners, 0.6 and 1.0 Hz | the source's spectrum (half its amplitude by 0.6 Hz) |
| which variant does best | the central part against its own true model, only |

**Decision.** Work the mute out from the start and the source (`slow_waves_mute`), and say
that the observed slow waves must also come after it, which these true models satisfy but
these runs couldn't check from the data. Until the western runs are done, every ranking in
the Marmousi2 entries holds for the central part only.

**Reproduce:** `just in marmousi2` (now also the western part's files), `just in
marmousi2-study mute_lines` (seconds).

## 2026-11-21 — An inversion with its forward wavefields rebuilt

**Question.** Rebuilding a gradient run's forward wavefield backwards from the band around
the layers (2026-11-20, "Rebuilding the forward wavefield backwards") needs less memory than
checkpoints. Can an inversion use it, keeping its memory plan exact?

**Method.** `[inversion] reconstruct = true` makes every gradient run of `Inversion::evaluate`
use `Strategy::Reconstruct` with the default segment, the length that needs least memory.
`Plan::inversion_memory` then counts the reconstruction's parts instead of the checkpoints'
(`counts_what_an_inversion_allocates`, now with the setting too), and `memory --invert`
prints the other way's peak. `rebuilding_the_wavefield_inverts_alike` inverts a small survey
both ways.

**Results.**

- The plan equals what the inversion holds with the setting, as without it.
- `examples/marmousi2/central.toml`: 2.44 GiB against 4.41 with checkpoints. The whole
  section's gradient still needs 8.19 GiB (15 with checkpoints), more than the dev
  container's 7.9.
- On the small survey the starting misfit is the same bit for bit (the forward runs are), the
  gradient differs from the replay's by $2.6 \times 10^{-6}$ of its norm, the misfit after the
  first step by $5.5 \times 10^{-6}$ of it, and $v_s$ by at most $9.4 \times 10^{-8}$.
- *Not run:* the same comparison on Marmousi2's central part over five iterations
  (`l2_rebuilt`, against `l2`'s model), where 425 m of water would drift five times more than
  rock. Open questions. (Run 2026-11-24: the models agree to about $10^{-5}$ of the update;
  "Marmousi2: an inversion with its wavefields rebuilt, against checkpoints".)

**Decision.** An option, off by default: checkpoints stay the default until an inversion of
that size has been compared both ways.

**Reproduce:** `just unit rebuilding_the_wavefield_inverts_alike`, `just unit
counts_what_an_inversion_allocates`, `just memory examples/marmousi2/central.toml --invert
--reconstruct`.

## 2026-11-21 — Marmousi2: the slow waves kept at 0.6 Hz, muted above

**Question.** Unmuted, the stages' 0.6 Hz iterations turned $v_s$ towards the truth at
750–1000 m; muted, they didn't, so it was the slow waves at frequencies low enough to be
within half a cycle that carried it (2026-11-21, "Marmousi2: the stages with the slow waves
muted"). Keep them in the 0.6 Hz stage and mute them in the later ones?

**Method** (`tests/marmousi2_study.rs`, `stages_muted_after_the_first`). The stages (below
0.6 Hz for two iterations, below 1 Hz for two, then all for one) with the slow waves' mute
(600 m/s from 0.87 s, taper 0.42 s) as a stage's own (`[[inversion.stage.mute]]`) in the
second and third stages only. Chosen, like every variant here, on this part against its own
truth (2026-11-21, "Marmousi2: the mute line without the truth, and a held-out part").

**Results.** After each stage:

| | L2 left | | 500–750 m | 750–1000 m | 1000–1500 m | 1500–2500 m | 2500–3125 m |
|---|---|---|---|---|---|---|---|
| 0.6 Hz (2) | 97.9% | $v_p$ | 0.992, +0.28 | 0.998, +0.07 | 1.015, −0.21 | 0.999, +0.05 | 1.001, −0.07 |
| | | $v_s$ | 0.970, +0.39 | 0.988, +0.16 | 1.029, −0.12 | 1.009, −0.07 | 1.000, +0.06 |
| then 1.0 Hz, muted (2) | 92.2% | $v_p$ | 0.979, +0.38 | 0.997, +0.11 | 1.017, −0.23 | 0.999, +0.05 | 1.001, −0.08 |
| | | $v_s$ | 0.949, +0.45 | 0.991, +0.14 | 1.031, −0.12 | 1.009, −0.07 | 1.000, +0.06 |
| then all, muted (1) | 91.3% | $v_p$ | 0.973, +0.41 | 0.995, +0.14 | 1.017, −0.23 | 0.999, +0.05 | 1.001, −0.08 |
| | | $v_s$ | 0.947, +0.45 | 0.991, +0.13 | 1.031, −0.12 | 1.009, −0.07 | 1.000, +0.06 |

The 0.6 Hz stage is the unmuted stages' to the last digit, as it should be. 3,608 s.

After five iterations ($v_p$ / $v_s$ error ratio):

| | L2 left | 500–750 m | 750–1000 m | 1000–1500 m |
|---|---|---|---|---|
| L2, muted | 85.7% | 0.916 / 0.970 | **0.971** / 1.004 | **0.999** / **1.001** |
| stages, unmuted | 93.2% | 0.976 / **0.946** | 1.001 / **0.987** | 1.020 / 1.044 |
| stages, muted | 89.2% | 0.972 / 0.966 | 0.995 / 1.009 | 1.003 / 1.004 |
| stages, muted after 0.6 Hz | 91.3% | 0.973 / 0.947 | 0.995 / 0.991 | 1.017 / 1.031 |

- *It keeps what the low frequencies gave $v_s$ above 1 km* (0.947 and 0.991, the unmuted
  stages' 0.946 and 0.987), which muting every stage lost (0.966 and 1.009).
- *And what they cost below 1 km.* The 0.6 Hz stage already moves $v_s$ at 1–1.5 km away
  (1.029, −0.12), the same slow waves that help above 1 km; the muted stages after it don't
  undo that (1.031), only stop it growing, where the unmuted ones went on to 1.044.
- $v_p$ follows the stages, not muted L2 (0.973, 0.995, 1.017 against 0.916, 0.971, 0.999):
  four of the five iterations are at low frequencies again.

**Decision.** No variant here recovers $v_s$ below 1 km, and each trades it against
something: muted L2 is the safest (nothing gets worse, $v_p$ improves down to 1 km), the
stages with the slow waves kept at 0.6 Hz the best for $v_s$ above 1 km at a cost of 3% below
it. These rankings come from this part and its own truth, and are untested on the held-out
part. The rest of the follow-ups were not run (open questions).

**Reproduce:** `just in marmousi2-study stages_muted_after_the_first` (60 minutes on four
cores).

## 2026-11-22 — Band-limited layers by default

**Question.** Band-limited layers put an interface where it lies (2026-11-21, "Band-limited
layers, opt in"), but flipping the default failed three gradient checks and the inversion's
set-up test. What do the checks need, and what does the flip change elsewhere?

**Method.** `config::Interfaces::BandLimited` made the default. The gradient checks
(`gradient_check` in `tests/physics.rs`) perturbed random signs at every node outside the
layers, but $v_p$ only below the model's fastest speed, which sets $\Delta t$ through the
Courant fraction and the layers' damping. Band-limited, the fastest is an overshoot plane
(3022 m/s) above a 3000 m/s layer, which a perturbation pushes past it. The checks now hold
$\Delta t$ and `pml_vp_max` at the base model's for the observed, base and perturbed runs, as
`Inversion::new` does, and perturb $v_p$ everywhere. The step was swept over
$h = 3.2 \times 10^{-2}$ to $2 \times 10^{-3}$ for $v_p$, with random signs and down the
gradient (random amounts with the gradient's signs, as the inversion's check does), under
both samplings. Then the whole of `just physics` under both samplings with the new checks;
for the misfits' checks (`inversion_gradient_check`), steps from $2 \times 10^{-2}$ to
$2.5 \times 10^{-3}$.

**Results.**

- *Random signs with $\Delta t$ held* fail: $v_p$'s extrapolated error scatters without a
  trend as $h$ shrinks, $2.5 \times 10^{-5}$ to $4.2 \times 10^{-4}$ at order 4 and
  $9.7 \times 10^{-5}$ to $7.3 \times 10^{-4}$ at order 8, and node sampling fails at
  $3.0 \times 10^{-4}$. Perturbing the fast lower layer too, which the old mask had kept
  fixed, leaves the predicted change cancelling between nodes ($4.0 \times 10^{-4}$ per unit
  $h$ against $v_s$'s $6.3 \times 10^{-3}$), and the finite differences' round-off, about
  $10^{-7}$ of a unit of misfit (the same as $v_s$'s), swamps it.
- *Down the gradient* the change can't cancel ($-4.9 \times 10^{-2}$ for $v_p$), and $v_p$'s
  extrapolated error falls as $h^4$, from $1.0 \times 10^{-4}$ at $3.2 \times 10^{-2}$ to a
  floor of about $3 \times 10^{-6}$ from $8 \times 10^{-3}$ down. At $h = 1.6 \times 10^{-2}$,
  band-limited (node-sampled):

  | | order 4 | order 8 | walls and water |
  |---|---|---|---|
  | $v_p$ | $8.5 \times 10^{-6}$ ($1.0 \times 10^{-5}$) | $2.3 \times 10^{-6}$ ($6.3 \times 10^{-6}$) | $5.4 \times 10^{-6}$ |
  | $v_s$ | $6.9 \times 10^{-5}$ ($7.1 \times 10^{-5}$) | $7.4 \times 10^{-5}$ ($7.2 \times 10^{-5}$) | $9.1 \times 10^{-6}$ |
  | $\rho$ | $2.5 \times 10^{-6}$ ($4.8 \times 10^{-6}$) | $1.2 \times 10^{-5}$ ($1.3 \times 10^{-6}$) | $5.1 \times 10^{-6}$ |

  Before, along random signs with $v_p$'s lower layer fixed (node-sampled): $v_p$
  $2.8 \times 10^{-4}$ / $3.9 \times 10^{-6}$, $v_s$ $2.1 \times 10^{-5}$ / $1.0 \times 10^{-5}$,
  $\rho$ $4.4 \times 10^{-5}$ / $1.6 \times 10^{-5}$; walls and water $6.4 \times 10^{-5}$,
  $1.2 \times 10^{-5}$, $1.6 \times 10^{-4}$. $v_s$ is now the largest, the same at both
  orders and samplings: the extrapolation's $h^4$ error along a direction that changes the
  misfit more nonlinearly, not round-off. The rebuilt wavefield's checks give the same
  numbers to $4 \times 10^{-7}$ of the predicted change.
- *Elsewhere in `just physics`*, node-sampled → band-limited: the long run drains to
  $2.7 \times 10^{-6}$ → $2.5 \times 10^{-6}$ (order 4) and $1.1 \times 10^{-6}$ →
  $9.6 \times 10^{-7}$ (order 8); the wavefield rebuilt backwards drifts $7.5$, $6.3$, $6.2$,
  $12 \times 10^{-7}$ → $8.0$, $9.0$, $8.0$, $10 \times 10^{-7}$ (orders 4 and 8, with a
  free surface and with layers all round) and $9.6 \times 10^{-7}$ → $1.8 \times 10^{-6}$ with
  walls and water; its gradient differs from the replay's by at most $7.1 \times 10^{-7}$ →
  $6.2 \times 10^{-7}$. The inversion's checks, down the gradient at one step: inversion
  gradient $4.0 \times 10^{-6}$ → $4.9 \times 10^{-6}$, filtered L2 $4.4 \times 10^{-7}$ →
  $2.4 \times 10^{-6}$, phase $4.2 \times 10^{-6}$ → $6.3 \times 10^{-7}$, filtered time
  shifts $1.3 \times 10^{-4}$ → $1.2 \times 10^{-4}$, muted filtered L2 $8.7 \times 10^{-7}$
  → $1.9 \times 10^{-5}$, windowed time shifts $2.6 \times 10^{-5}$ → $9.1 \times 10^{-5}$.
  Swept, band-limited, the muted L2 check gives $1.3 \times 10^{-4}$, $1.9 \times 10^{-5}$,
  $2.5 \times 10^{-5}$ and $5.3 \times 10^{-5}$ at $h = 2 \times 10^{-2}$, $10^{-2}$,
  $5 \times 10^{-3}$ and $2.5 \times 10^{-3}$: the extrapolation's $h^4$ error above,
  round-off growing as $1/h$ below, about equal at $10^{-2}$. Filtered L2 ($1.9 \times
  10^{-4}$, $2.4 \times 10^{-6}$, $2.3 \times 10^{-5}$, $6.2 \times 10^{-5}$) and the
  inversion gradient ($7.5 \times 10^{-5}$, $4.9 \times 10^{-6}$, $5.7 \times 10^{-7}$,
  $1.4 \times 10^{-5}$) behave the same. A wrong gradient would leave the same error at every
  step, so these moved within the checks' noise; node sampling's $8.7 \times 10^{-7}$ was a
  lucky cancellation in it. Every other test is unchanged (its configuration sets the
  sampling, or has no interface between layers).
- The inversion's set-up test expected `pml_vp_max` to be the rock's 2000 m/s; the start's
  fastest node is now the overshoot beside the water, which it now takes.

**Decision.** Band-limited layers are the default; `interfaces = "nodes"` keeps node
sampling. The gradient checks hold $\Delta t$ and the layers' damping and go down the
gradient, and their tolerance is tightened from $3 \times 10^{-4}$ to $1.5 \times 10^{-4}$,
twice the largest error. The doubled drift of the wavefield rebuilt beside water is recorded,
not traced: band-limiting puts solid nodes with small shear moduli beside the water, a
candidate (the gradient's notes on $\mu \ll \lambda$), untested.

**Reproduce:** `just physics`; for node sampling, move `#[default]` in `config::Interfaces`
back to `Nodes`. The step sweeps need a loop over $h$ in `gradient_check` and
`inversion_gradient_check` (not committed).

## 2026-11-22 — The free surface against Lamb's problem

**Question.** The P-doubling test's 2.5% is mostly its reference (2026-11-20, "The free
surface on finer grids"), which leaves the free surface's own accuracy unmeasured. How
accurate is it for the wave it exists for, the Rayleigh wave, and is the error the free
surface's or the stencil's?

**Method.** Lamb's problem for a Poisson solid: Pekeris's closed form for the vertical
displacement of the surface under a vertical step force on it [@pekeris1955], his eqs.
45–48 (read from the article's scanned pages; `tests/common/lamb.rs`). With
$\tau = \beta t/r$ and $\gamma^2 = (3 + \sqrt3)/4$, the Rayleigh pole, a downward force $F$
moves the surface down by $w = 0$ for $\tau < 1/\sqrt3$;
$$
w = \frac{F}{32\pi\mu r}\left[6 - \frac{\sqrt3}{\sqrt{\tau^2 - \frac14}}
  - \frac{\sqrt{3\sqrt3 + 5}}{\sqrt{\gamma^2 - \tau^2}}
  + \frac{\sqrt{3\sqrt3 - 5}}{\sqrt{\tau^2 - \frac{3 - \sqrt3}{4}}}\right]
$$
for $1/\sqrt3 < \tau < 1$; $w = \frac{F}{16\pi\mu r}\left[6 - \sqrt{3\sqrt3 + 5}/\sqrt{\gamma^2 -
\tau^2}\right]$ for $1 < \tau < \gamma$; and Boussinesq's $3F/(8\pi\mu r)$ after. The
bracket vanishes exactly at the P arrival, both sides at $\tau = 1$ are $(1 - \sqrt3)/8$, the
forms agree with his integrals (eqs. 33–44), and the static value is Boussinesq's for
$\nu = \frac14$: a wrong coefficient breaks one of these. For the solver's Ricker force the
velocity is $w * s''$: Simpson's rule after $t' = t_P + u^2$ and $t' = t_R - u^2$, which
removes the Rayleigh singularity, plus the static part's exact term; 2000 and 4000 intervals
agree to $7 \times 10^{-12}$ of the peak (`lambs_reference_is_consistent`). The survey: $v_p$
3000 m/s, $v_s = v_p/\sqrt3$, $\rho$ 2500 kg/m³; a downward force at 15 Hz on the surface;
$v_z$ on the surface at 150, 240 and 330 m along $x$ and 254.6 m on the diagonal; the
Rayleigh wave ($0.9194\,v_s$) 106 m long at $f_0$, 10.6 points at $h = 10$ m, where the
shortest S wavelength has the default's 4.6. Spacings 10 to 2.5 m at order 8; order 16 and
the time error kept ($C = 0.3$) at 10 and 5 m. Per receiver: the L2 misfit, the amplitude
$\langle v, w\rangle/\langle w, w\rangle$, the largest swing's lateness and its size.

**Results** (order 8; the four receivers' range, or 150 / 240 / 330 m):

| $h$ (m) | misfit | swing's size | lateness (samples) | Rayleigh speed |
|---|---|---|---|---|
| 10 | 15.5 / 21.4 / 28.2% (diagonal 22.5%) | +4.8–5.0% | −1.00 / −1.44 / −1.90 | +0.99% |
| 7.5 | 10.5 / 13.9 / 18.1% | +5.6–5.8% | −0.81 / −1.18 / −1.56 | +0.62% |
| 5 | 6.4 / 7.8 / 9.5% | +5.2–5.4% | −0.57 / −0.83 / −1.11 | +0.30% |
| 3.75 | 4.8 / 5.4 / 6.3% | +4.4–4.5% | −0.42 / −0.62 / −0.84 | +0.17% |
| 2.5 | 3.3 / 3.5 / 3.8% | +3.3% | −0.27 / −0.40 / −0.55 | +0.077% |

The speed comes from the lateness's growth with distance and the sample interval (1.2465 ms
at 10 m, in proportion to $h$). Order 16 is within 1 point of order 8 (15.1 / 20.7 / 27.2% at
10 m, 6.3 / 7.6 / 9.3% at 5 m), and with the time error kept the misfits are within 0.6 point
and the lateness, in time, the same.

- *The free surface, not the stencil, sets the error.* The stencil's order and the time
  error change nothing; at the default's 4.6 points per shortest S wavelength, where the
  interior's S-wave error is about 1% after ten wavelengths, the Rayleigh wave misfits by
  15–28% within three of its own.
- *Two errors.* The Rayleigh wave runs fast, 0.99% at 10 m, falling as $h^{1.8}$ over the
  range: second order, the misfit's part that grows with distance. And its largest swing is
  3–6% too large at every distance, falling only to 3.3% at 2.5 m: first order at the fine
  end (4.4% at 3.75 m, 3.3% at 2.5 m), and at the coarse end below that trend, where the
  phase error also lowers the peak.
- *A candidate for the second, not tested:* the free surface mirrors $v_z$ evenly, so $v_z$
  continued above the surface has a kink, and interpolating across it reads $v_z(0) + c\,h\,
  \partial_z v_z$ with $c = 0.342$ at order 8 (0.375 at 4, 0.329 at 16: hardly the order's,
  as measured). A Rayleigh wave's $\partial_z w / w$ at the surface is $0.227k$ for
  $\nu = \frac14$, so a receiver over-reads by $0.078\,kh$ and the surface force, its
  transpose, over-excites by as much: 2.3% at 5 m and 1.2% at 2.5 m together at $f_0$,
  against 5.2% and 3.3% measured. It accounts for part of it, falling at the same rate.

**Decision.** The free surface is the solver's largest error for surface waves, by far: a
Rayleigh wave needs about 40 points per wavelength at $f_0$ for a few per cent, where body
waves need 4.5 at the shortest. `a_force_on_the_free_surface_matches_lambs_problem` pins
the 10 m numbers (limits just above them). A higher-order free surface [@mittet2002;
@kristek2002] or summation-by-parts closures [@oreilly2017] would address the speed; the
swing's size may be the surface's sampling of $v_z$ (the candidate above), which a free
surface with $v_z$ on the surface plane would also answer. Surface data, Marmousi2's
included, carry this error into the misfit.

**Reproduce:** `just physics lamb`; `just accuracy lambs_problem_on_finer_grids` (about two
and a half minutes).

## 2026-11-22 — Sources and receivers between nodes

**Question.** Sources and receivers snapped to the nearest node, and a receiver 0.4 of a node
off its position cost 26% misfit at 3.9 points per wavelength (2026-11-20, removing time
dispersion). Can they sit where they are, with the scheme's accuracy and the adjoint's
exactness intact?

**Method.** A position becomes grid coordinates, position over spacing, whole numbers on a
node (within $10^{-9}$, so a configuration on nodes reads exactly as before). Each component
lives on a lattice, nodes or half nodes along each axis (`Field::lattice`), and is read at a
point from the $2L$ values of its lattice around it along every axis the point lies off: on
its own lattice point, the value; a staggered value onto a node, the stencil's interpolation as
before; otherwise weights for an offset $\delta \in (0, 1)$ past the value below
(`Stencil::interpolation_at`). With Taylor's coefficients those are Lagrange's; with optimised
ones they minimise the largest complex error $|\sum_j w_j e^{i x (j - \delta)} - 1|$ for
$kh$ up to the stencil's band, weighted by $(k/k_\text{band})^{-2}$ as the node's are, exact
for constants and linear functions (Lagrange's weights plus a null-space combination, Lawson's
algorithm on the modulus). The three axes multiply (`Solver::field_weights`); along $z$, values
above a free surface fold onto their mirror images with the field's parity, odd for $\sigma_{zz}$,
$\sigma_{xz}$ and $\sigma_{yz}$. A source spreads with the same weights, so it stays the
receiver's transpose; Hicks (2002) does the same with Kaiser-windowed sincs [@hicks2002].
Pressure's adjoint
source applies the node's stiffness (plane stress on the surface) to the three normal weights,
and the gradient's pressure residuals go to the adjoint strain with them. Checked against
Stokes's solution, reciprocity, the dot-product test and a gradient check, each with the
source and receivers between nodes, some between the free surface and the first node.

**Results.**

- *Weights.* The design at $\delta = \frac12$ reproduces the node's to $10^{-13}$ (order 4 to
  16). Largest error over the band, Lagrange → optimised: order 4 (band 1.5) $3.3 \times
  10^{-2}$ → $2.0 \times 10^{-2}$ at $\delta = 0.1$, $9.1 \times 10^{-2}$ → $5.7 \times 10^{-2}$ at
  0.37; order 8 (1.8) $1.8 \times 10^{-2}$ → $1.6 \times 10^{-3}$, $5.1 \times 10^{-2}$ →
  $4.7 \times 10^{-3}$; order 16 (2.4) $4.4 \times 10^{-2}$ → $7.2 \times 10^{-4}$,
  $1.3 \times 10^{-1}$ → $2.1 \times 10^{-3}$. Minimising the real and imaginary parts
  separately came out $7 \times 10^{-2}$ away from the node's design at $\delta = \frac12$: part by
  part the optimum isn't unique, and the design picked a needless phase error.
- *Accuracy.* A force at (40.37, 39.62, 40.89) nodes and receivers the same offsets from it
  (6.9 points per S wavelength): 0.001–0.007% against Stokes's solution, as on the nodes
  (0.001–0.005%). With the force on a node and receivers 0.37–0.41 of a node off, snapped:
  10.0%, 17.7%, 7.7% ($v_x$) and 0.14% ($v_y$); exact: 0.003%, 0.005%, 0.007% and 0.001%.
- *Transposes.* Reciprocity, one end on the surface between nodes or 4.1 m below it (between
  the surface and the first node): $\le 6.8 \times 10^{-7}$ at orders 4 and 8, as on the nodes
  ($\le 7.9 \times 10^{-7}$). The dot-product test with every receiver and the source moved by
  (3.7, −2.4, 4.1) m, all four components: $3.8 \times 10^{-10}$ (rigid walls) and
  $5.9 \times 10^{-8}$ (absorbing layers). The first version gave pressure's adjoint the row
  sums of the stiffness, right only where the three normal weights are equal; moved along $z$
  only, pressure broke the test ($2.7 \times 10^{-3}$, $2.1 \times 10^{-2}$), moved sideways it
  didn't ($1.8 \times 10^{-8}$). The gradient check with the gradient model moved by
  (−3.7, −2.4, 4.1) m, pressure recorded 4.1 m below the surface: $v_p$ $2.4 \times 10^{-7}$,
  $v_s$ $7.3 \times 10^{-5}$, $\rho$ $1.1 \times 10^{-5}$, as on the nodes.
- *The existing tests.* Lines whose receivers fell between nodes (and three sources) moved
  them: the sphere inversion left 18.2% of its misfit instead of 17.9%, and the misfits'
  check failed for the time shifts ($6.2 \times 10^{-3}$, $3.9$–$5.4 \times 10^{-3}$ for $h$
  from $2 \times 10^{-3}$ to $1.25 \times 10^{-4}$). One-sided differences show why: forward
  and backward slopes differ by about 1% as $h \to 0$, so one trace's correlation peak changes
  lag at the start itself, where the shift has a kink; filtered L2's errors stay symmetric
  ($\pm 1.8 \times 10^{-2}$ at $1.25 \times 10^{-4}$). The tests now list the positions they
  used to snap to, and every number is as before.

**Decision.** Sources and receivers sit where they are; there is no option to snap. On nodes
nothing changes. A point between nodes reads $(2L)^3$ values (512 at order 8), which costs
nothing next to a step.

**Reproduce:** `just unit fd::tests::interpolates_from_any_offset`; `just physics between`
(the force, reciprocity, dot-product and gradient checks between nodes).

## 2026-11-23 — Marmousi2: the held-out western part

**Question.** Every Marmousi2 ranking so far came from the central 6 km, judged against its
own truth (2026-11-21, "Marmousi2: the mute line without the truth, and a held-out part").
Do the rankings hold on the western 6 km, which none of the choices saw?

**Method** (`tests/marmousi2_study.rs`: `west_l2`, `west_muted_l2`,
`west_stages_muted_after_the_first`, `seafloor_speeds`, `west_misfit_by_mute`). The same
survey, start (smoothed over 250 m) and numerics on $x$ = 0–6 km (`Case::part("west", …)`),
three of the central part's variants:
- L2 for five iterations;
- L2 with the slow waves muted, before the line that `slow_waves_mute` works out from the
  west's own start (690 m/s from 0.87 s, taper 0.42 s);
- the stages with that mute from the second stage on.

The central part's runs used 600 m/s, chosen with its truth in view; the rule gives 616 there.
Every source and receiver lies on a node and the models are files, not layers, so the
defaults changed since those runs (band-limited layers, exact positions) don't apply.
Two diagnostics followed:
- `seafloor_speeds`, which uses the truth, compares the start's S speeds under the seafloor
  with the truth's in both parts.
- `west_misfit_by_mute` splits the west's starting misfit and gradient by its mute, as
  "Marmousi2: what the mute removes" did on the central part.

The inversions ran in the dev container with 8 CPUs, and the split ran natively on the Mac
(`just mac-test`).

**Results.** After five iterations, error ratio and cosine with the true update on the
western part:

| | L2 left | | 500–750 m | 750–1000 m | 1000–1500 m | 1500–2500 m | 2500–3125 m |
|---|---|---|---|---|---|---|---|
| L2 | 17.3% | $v_p$ | 1.227, +0.05 | 0.988, +0.26 | 0.986, +0.19 | 0.999, +0.06 | 1.004, −0.33 |
| | | $v_s$ | 0.714, +0.70 | 1.018, +0.31 | 0.935, +0.36 | 1.000, +0.05 | 0.998, +0.28 |
| L2, muted | 41.3% | $v_p$ | 1.167, +0.13 | 1.055, +0.20 | 1.013, +0.17 | 0.997, +0.10 | 1.001, −0.20 |
| | | $v_s$ | 0.791, +0.61 | 0.978, +0.26 | 0.999, +0.04 | 1.001, −0.17 | 0.999, +0.29 |
| stages, muted after 0.6 Hz | 64.0% | $v_p$ | 1.006, +0.09 | 0.987, +0.24 | 0.992, +0.25 | 1.000, +0.01 | 1.003, −0.39 |
| | | $v_s$ | 0.911, +0.48 | 0.993, +0.13 | 0.983, +0.18 | 1.000, +0.06 | 0.998, +0.38 |

Next to the central part ($v_p$ / $v_s$ error ratio):

| | part | L2 left | 500–750 m | 750–1000 m | 1000–1500 m |
|---|---|---|---|---|---|
| L2 | central | 57.9% | 0.900 / 0.973 | 0.994 / 1.095 | 1.020 / 1.070 |
| | west | 17.3% | 1.227 / 0.714 | 0.988 / 1.018 | 0.986 / 0.935 |
| L2, muted | central | 85.7% | 0.916 / 0.970 | 0.971 / 1.004 | 0.999 / 1.001 |
| | west | 41.3% | 1.167 / 0.791 | 1.055 / 0.978 | 1.013 / 0.999 |
| stages, muted after 0.6 Hz | central | 91.3% | 0.973 / 0.947 | 0.995 / 0.991 | 1.017 / 1.031 |
| | west | 64.0% | 1.006 / 0.911 | 0.987 / 0.993 | 0.992 / 0.983 |

The west's balanced L2 misfit at the start is 0.234, against the central part's 0.660. The
three runs took 2,027, 2,095 and 2,676 s; the central ones took 2,634, 2,769 and 3,608 s on
4 CPUs.

Under the seafloor (`seafloor_speeds`):

| | $v_s$ at 475 m, true / start | 250 m below, harmonic mean, true / start | delay over 1 km, mean (largest) |
|---|---|---|---|
| central | 345 / 407–684 m/s | 414 / 490 m/s | 0.37 s (0.80) |
| west | 345 / 438–448 m/s | 431 / 448 m/s | 0.09 s (0.13) |

Half a period of the 1.2 Hz Ricker is 0.42 s.

- *The central ranking doesn't hold.*
  - In the west, L2 alone recovers $v_s$ best at every depth down to 1.5 km (0.714 with a
    cosine of +0.70 at the top, 0.935 with +0.36 at 1–1.5 km). In the centre it was the one
    variant that turned $v_s$ away below 750 m.
  - Muted L2, the centre's safest, makes $v_p$ worse in the west: 1.167 at 500–750 m and 1.055
    at 750–1000 m.
  - The stages muted after 0.6 Hz are the only variant that makes nothing worse in the west:
    the worst is $v_p$'s 1.006 at the top. In the centre they made $v_s$ worse at 1–1.5 km
    (1.031).
  - No variant is best, or harmless, on both parts.
- *Why: in the west the start predicts the slow waves within half a cycle.*
  - There the start's $v_s$ in the 250 m under the seafloor is within 4% of the truth's.
    Waves at the truth's speeds fall 0.09 s behind the start's over 1 km, and 0.13 s at most.
  - In the centre the start is 18% fast. Smoothing over 250 m has mixed faster rock into the
    soft sediments: the start has 407–684 m/s where the truth has 345 m/s all along the
    seafloor. There the slow waves fall 0.37 s behind per km, and up to 0.80 s, so they pass
    half a cycle within about a kilometre.
  - The slow waves are the same kind of data in both parts. Within half a cycle (the west,
    and the centre's 0.6 Hz stage) their gradient points towards the truth; beyond it (the
    centre at full band) it points away.
- *The split agrees* (`west_misfit_by_mute`). The mute removes about as much of the west's
  misfit as of the centre's (73.7% against 64.6%), but what it removes points the other way.
  Raw cosines with the true update, west / central (the centre's from "Marmousi2: what the
  mute removes"):

  | | 500–750 m | 750–1000 m | 1000–1500 m |
  |---|---|---|---|
  | $v_p$, removed | +0.05 / +0.02 | +0.12 / −0.28 | −0.02 / −0.20 |
  | $v_s$, removed | +0.58 / −0.01 | +0.00 / −0.12 | +0.24 / −0.18 |
  | $v_p$, kept | +0.13 / +0.45 | +0.32 / +0.31 | +0.14 / +0.02 |
  | $v_s$, kept | +0.48 / +0.20 | −0.20 / +0.08 | +0.04 / −0.10 |

  - *The removed part holds almost all of $v_s$'s gradient below 750 m*, as in the centre:
    the kept part is 8% as large at 750–1000 m and 3% at 1–1.5 km. That is why muted $v_s$
    barely moves there (0.978, 0.999).
  - *What the mute keeps points away at 750–1000 m* ($-0.20$).
  - *$v_p$'s whole gradient at 500–750 m is only weakly aligned* (+0.09, against the centre's
    +0.37), so the L2 runs' $v_p$ there drifts away over the iterations rather than starting
    off wrong.
  - $v_x$ carries 94.4% of the west's misfit.
  - The run took 1,060 s natively on the Mac, with the container's CI sharing the cores for
    part of it, so it isn't a benchmark.
- *Fitting more of the data recovers more here.* L2 leaves 17.3% of the west's misfit and
  is its best variant for $v_s$. The central part's "L2 fits the most data and recovers the
  least" held only for that part.
- *$v_p$ at 500–750 m moves away in both L2 runs*, with cosines near zero (+0.05, +0.13): its
  update there is mostly not the true one. The stages, at low frequencies for four of their
  five iterations, leave it alone (1.006).

**Decision.** Whether to mute the slow waves depends on whether the start predicts them
within half a cycle, and the mute rule doesn't ask that. `slow_waves_mute` places the line
ahead of the start's slowest S waves, which is the same in both parts. So no mute becomes a
default, and no ranking from either part is general.

What the decision needs is a test, from the data, of whether the slow waves are
cycle-skipped. One option is the lag between observed and synthetic traces after the line,
offset by offset, measured once at each stage's start and held through the stage, so that
the misfit stays smooth within it. A better shallow start would make the test pass more often.
Neither is written.

The whole-section inversion (`examples/marmousi2/whole.toml`) mutes at 600 m/s throughout.
It covers both parts, so its result should be judged region by region, not only by depth.

**Reproduce:** `just in marmousi2-study west_l2`, `west_muted_l2` and
`west_stages_muted_after_the_first` (35–45 minutes each on 8 CPUs), `seafloor_speeds`
(seconds), `west_misfit_by_mute` (18 minutes natively on the Mac); or natively, `just in mac-test
marmousi2_study`, then `out/mac/marmousi2_study <name> --ignored --exact --nocapture`.

## 2026-11-23 — Marmousi2 in the literature

**Question.** Which published inversions of Marmousi2 [@martin2006] can ours be set
against, how do their setups differ from ours, and which comparison is worth making first?

**Method.** The paper library, web searches and Crossref; the full texts that could be
opened were read (three in the library, three as the publishers' or arXiv's HTML by reading
agents), the rest as abstracts. `docs/marmousi2-literature.md` has the table, the notes on
each study and the proposed comparisons; every work is in the bibliography with how deeply
it was read. Proposed configurations were sized with `vibe-rations memory --invert` (the
macOS build; the 20 m surveys on a stand-in layered model with $\Delta t$ and `pml_vp_max`
pinned, which fixes the plan). The effect of the literature's model change on ours was
measured on `convert.sh`'s files: $v_s$ raised to $v_p/\sqrt3$ wherever it is lower on a
solid node, i.e. Poisson's ratio capped at 0.25.

**Results.**

- *The elastic Marmousi2 inversions read in full raised $v_s$ first*: Poisson's ratio
  capped at 0.25 [@kohn2012], set to 0.25 with the water removed [@jeong2012], set to about
  0.25 [@cao2022]; another runs 20 Hz on a 20 m grid, where the original sediments would get
  under a point per S wavelength, and doesn't say [@sun2022]. None found says it keeps the
  soft sediments. The cap changes $v_s$ at 94.7% of the central part's solid nodes (90.7% of
  the whole section's), by 131%, 98%, 59%, 26% and 9% on average in the bands 475–750,
  750–1000, 1000–1500, 1500–2500 and 2500–3500 m; the slowest $v_s$ goes from 315 to
  886 m/s (start: 407 to 953; western part: 345 to 599; whole section: 303 to 599, set by
  nodes with $v_p$ near 1040 m/s). The capped model has no slow seafloor waves.
- *They invert 8–20 Hz from 6–10 Hz sources on 20 m grids*, lowest bands from 0.2 Hz
  [@jeong2012], 0–5 Hz [@cao2022], up to 2 Hz [@kohn2012], 2–4 Hz [@sun2022], 3–5 Hz
  [@pan2018]; for 40–1300 iterations, from starts as smooth as ours or smoother (Gaussian
  800 m [@kohn2012], a linear gradient [@jeong2012]).
- *One reports our measure*: $E_p = \lVert m_p^\text{est} - m_p^\text{true} \rVert /
  \lVert m_p^0 - m_p^\text{true} \rVert$ was 0.83, 0.72 and 1.04 for $v_p$, $v_s$ and $\rho$
  after 40 L-BFGS iterations at 3–8 Hz on a 3.4 × 1.2 km section with surface acquisition
  [@pan2018]. Cao et al. report relative data residuals of 2.25–2.86% [@cao2022]; the rest,
  pictures. None reports a cosine with the true update or anything by depth.
- *Soft seabeds*: on a Valhall synthetic elastic inversion needed low frequencies, a
  hierarchy of parameters and components [@brossier2009b] and, with noise, $\ell_1$
  [@brossier2010]; ocean-bottom field studies stage $v_p$ from $v_z$ before $v_s$ from PS
  waves [@sears2008; @sears2010], or mute the Scholte waves [@warner2013; @prieux2013b];
  surface waves are tapered out [@yuan2014], fitted by envelope first [@yuan2015], or
  stripped by offset and frequency [@masoni2016].
- *Memory of the proposed comparisons*: our survey on the capped model is `central.toml`'s
  grid and pins, 4.41 GiB (2.44 rebuilt); Köhn et al.'s streamer survey at 20 m ($535
  \times 41 \times 175$ nodes, 3,750 steps, 100 shots) 19 GiB (11 rebuilt), Sun & Demanet's
  ocean-bottom one (50 shots) 18 GiB (10 rebuilt), about 11 minutes per shot and iteration
  on four cores; the whole section 48 GiB (22 rebuilt).

**Decision.** The published elastic Marmousi2 successes are on a model with roughly
doubled shallow $v_s$ and no slow seafloor waves, at frequencies ten times ours, so they
are not yet comparable with our results. First comparison: our central and western parts
with the cap applied to truth and start (start capped before smoothing), everything else
as `central.toml`: five L2 iterations, with and without the mute, error ratio and cosine by
band. If $v_s$ then moves towards the truth below 750 m, the soft seabed is what separates
us from the literature, and the soft-seabed studies above say what to try next; if not,
the band or the iterations do. It needs a cap in `extrude` (applied to the fine samples
before averaging) or a script; no solver change. Then the literature's band on our survey
(the capped model holds about 7.9 Hz at 25 m), Pan et al.'s $E_p$ over all free nodes at
5–40 iterations, and Cao et al.'s residual measure on `whole.toml`'s results. Nothing was
run; this commit changes only `docs/`.

*Revised after the held-out part* (previous entry). The note's "Keep or mute the slow
waves" now puts a check from the data first, and this control second. The check is the
phase of the slow waves against the start's synthetics, offset by offset at 0.6 and
1.2 Hz [@fichtner2011; @shah2012], the latter read as an abstract only. It needs only
forward runs, about 80 s per part. The mute then goes beyond the offset where that phase
stops being continuous.

**Reproduce:** the cap's effect, in numpy (run here with a plain-Python `.npy` reader, the
host having no numpy): `vs_hard = np.where(vs > 0, np.maximum(vs, vp / np.sqrt(3)), 0)` on
`examples/marmousi2/25m-central/{vp,vs}.npy` and the others, then the share of solid nodes
changed and the mean of `vs_hard / vs - 1` by depth band (row × 25 m). The memory:
`vibe-rations memory <config> --invert` on the configurations the note describes.

## 2026-11-23 — The free surface's amplitude excess is the mirrored $v_z$'s kink

**Question.** Lamb's problem found the Rayleigh wave's largest swing 3–6% too large and left
a hypothesis untested (2026-11-22, "The free surface against Lamb's problem"): the free
surface mirrors $v_z$ evenly, so $v_z$ continued above the surface has a kink, and
interpolating $v_z$ onto the surface reads it as if from $c\,h$ below, $c = 0.340$ at order 8
(the mean depth of a surface receiver's weights over $h$, from the solver's own operator). A
surface force, the receiver's transpose, acts from the same depth. A Rayleigh wave's $v_z$
grows with depth at the surface, $\partial_z w / w = 0.227\,k$ for $\nu = \frac14$, so each end
should over-read by about $0.227\,k\,c\,h$: in proportion to frequency at a fixed $h$. And
$v_z$ itself lives at $h/2$, so a source and receivers there, which interpolate nothing along
$z$, should show none of it. (The earlier entry's "2.3% at 5 m together" is that estimate for
one end; both together are 4.6%.)

**Method.** Lamb's survey (`tests/common/lamb.rs`) at a fixed $h = 5$ m and $f_0$ = 7.5, 10, 15,
20 and 30 Hz ($kh$ = 0.15–0.59 for the Rayleigh wave at $f_0$): (a) source and receivers on the
surface, as before; (b) the source on the surface, the receivers at $h/2$; (c) both at $h/2$,
where the exact positions (2026-11-22, "Sources and receivers between nodes") put them on
$v_z$'s own points. The exact answers at $h/2$ are Pekeris's surface traces taken down by the
Rayleigh wave's depth dependence, $(q\,e^{-qkz} - a\,e^{-skz})/(q - a)$ with $q = 0.8475$,
$s = 0.3933$, $a = 2q/(1 + s^2)$, once per end at $h/2$, frequency by frequency: exact for the
Rayleigh wave, which carries the swing, an approximation for the much weaker body waves. Each
end's own excess comes from pairs of runs that differ in that end alone, receiver by receiver:
the receiver's (a) − (b), the source's (b) − (c), so whatever else the free surface gets wrong
cancels. Measured: the largest swing's size and the amplitude $\langle v, e\rangle/\langle e, e
\rangle$ against the exact answer $e$ (`the_free_surface_excess_by_frequency_and_depth` in
`tests/accuracy.rs`).

**Results.** The largest swing's excess (the four receivers' range):

| $f_0$ (Hz) | $kh_R$ | (a) both on the surface | (b) receivers at $h/2$ | (c) both at $h/2$ |
|---|---|---|---|---|
| 7.5 | 0.148 | +3.2 to +4.0% | +1.5 to +2.5% | −0.2 to +0.9% |
| 10 | 0.197 | +4.0 to +4.8% | +1.8 to +2.7% | −0.5 to +0.6% |
| 15 | 0.296 | +5.2 to +5.4% | +1.9 to +2.1% | −1.3 to −1.1% |
| 20 | 0.395 | +5.6 to +6.3% | +1.4 to +2.1% | −2.6 to −1.9% |
| 30 | 0.592 | +4.9 to +6.1% | −0.6 to +0.6% | −5.7 to −4.6% |

Each end's own excess, the same at every receiver to within 0.2 point:

| $f_0$ (Hz) | receiver's end | source's end | per unit $kh_R$ | $0.227\,c\,kh_R$ at $f_0$ |
|---|---|---|---|---|
| 7.5 | +1.55 to +1.71% | +1.54 to +1.68% | 11.1 / 11.0% | 1.15% |
| 10 | +2.12 to +2.26% | +2.08 to +2.21% | 11.3 / 11.0% | 1.53% |
| 15 | +3.27 to +3.29% | +3.17 to +3.19% | 11.1 / 10.8% | 2.29% |
| 20 | +4.17 to +4.20% | +4.01 to +4.04% | 10.6 / 10.2% | 3.04% |
| 30 | +5.44 to +5.49% | +5.15 to +5.20% | 9.3 / 8.8% | 4.52% |

The amplitude measure gives about half a point less per unit $kh_R$ (10.2–10.8% to
20 Hz, 7.8–8.2% at 30 Hz).

- *The hypothesis holds.* Each end on the surface over-reads by the same amount, source and
  receiver alike as reciprocity requires, in proportion to frequency at a fixed $h$: about
  11% per unit $kh_R$, flat to $kh_R = 0.4$ and slightly less at 0.6. With both ends at $h/2$
  the excess is gone at low frequency (+0.4% and +0.1% at 7.5 and 10 Hz, within the body
  waves' approximation). At 15 Hz and 5 m the two ends make +6.5% of Lamb's +5.3%.
- *Its size is about 1.45 times the first-order estimate at $f_0$ up to 20 Hz* (1.2 at
  30 Hz), as if the swing were read at about $1.45 f_0$. That is plausible, since the swing of
  a velocity trace weighs frequencies above $f_0$, but not traced.
- *What is left is another error, negative and growing about as $f^2$:* in (c) $-1.2$%,
  $-2.3$% and $-5.2$% at 15, 20 and 30 Hz. With the kink's excess it partly cancels, which is
  why the total in (a) stopped growing with frequency and why Lamb's excess looked smaller
  than first order at the coarse end. From 15 Hz up, the misfit in (c) is within 1.5 points of
  (a)'s (5.0–8.5% against 6.4–9.5% at 15 Hz): there the Rayleigh wave's speed error, not the
  kink, sets it. At 7.5 and 10 Hz the kink is half the misfit (1.0–1.6% in (c) against
  3.2–3.8% in (a)).

**Decision.** The amplitude excess is not a propagation error: it is how sources and
receivers on the surface read and act on $v_z$ across the even mirror. A better free surface
has two separate things to fix. (1) The surface's sampling of $v_z$: continue $v_z$ above the
surface with the slope the boundary condition gives, $\partial_z v_z =
-\frac{\lambda}{\lambda + 2\mu}(\partial_x v_x + \partial_y v_y)$ from $\sigma_{zz} = 0$,
either in the point operators alone (read $v_z(0) \approx v_z(h/2) - \frac h2\,\partial_z v_z$,
a linear functional whose transpose is the source, so the adjoint stays exact), or in the
mirror itself (the image method, which must then pass the reciprocity and long-run tests).
The point operators are the smaller change and would remove an error of 11% per unit $kh_R$
per end. (2) The Rayleigh wave's speed, about 1% fast at 10 m (second order), and the smaller
negative amplitude error of case (c), which a free surface of higher order must address
[@mittet2002; @kristek2002; @oreilly2017].

**Reproduce:** `just accuracy the_free_surface_excess_by_frequency_and_depth` (about four
minutes on four cores).

## 2026-11-23 — Continuing the mirrored fields with their slope

**Question.** The free surface's amplitude excess is the evenly mirrored $v_z$'s kink, read
by sources and receivers on the surface (the entry before). Can the point operators read
across the surface without it, keeping sources the receivers' transpose and the gradient
exact?

**Method.** Where an interpolation along $z$ folds values from above the surface onto their
mirror images with even parity (velocities, $\sigma_{xx}$, $\sigma_{yy}$, $\sigma_{xy}$),
continue the field with its slope instead of flat: $v(-\zeta) = v(\zeta) - 2\zeta\,\partial_z
v(0)$ to second order, so each folded weight $w$ adds $-2w\zeta$ times the slope, estimated
one-sidedly from the field's first values below (Lagrange's derivative at the surface;
`Solver::continue_with_slope`). Odd fields vanish on the surface and keep their fold. The
operator is linear in the fields and independent of the model, so sources, built from the same
weights, stay the receivers' transpose and the gradient needs nothing new. The slope of the
boundary condition, $\partial_z v_z = -\frac{\lambda}{\lambda+2\mu}(\partial_x v_x +
\partial_y v_y)$, would be exact for every wave but makes the operators depend on the model.
Points for the slope: the reading error per travelling component $e^{i k_z z}$ (true value 1)
for order 8's interpolation onto the surface:

| $k_z h$ | mirror | one-sided, 4 points | 6 points |
|---|---|---|---|
| 0.2 | 6.9% | 0.05% | 0.02% |
| 0.4 | 13.8% | 0.72% | 0.24% |
| 0.6 | 21.0% | 3.5% | 1.7% |
| 0.8 | 28.5% | 10.6% | 7.6% |
| 1.0 | 36.5% | 24.6% | 24.3% |
| 1.2 | 45.4% | 47.9% | 62.2% |
| 1.4 | 55.2% | 82.4% | 134% |

and on an even profile $\cos(k_z z)$, which a vertically incident P wave has (true slope 0,
which the mirror reads exactly): 4 points $-2.3$% at $k_z h = 0.6$, 6 points 0.06%. Six
points, then: better than the mirror up to $k_z h \approx 1.1$, worse above, where at the
default's sampling only nearly vertical S above about $2.2 f_0$ lies. Measured with
`the_free_surface_excess_by_frequency_and_depth` (now with (a′) and (b′), (a) and (b) with
the continuation), Lamb's survey on finer grids and the physics tests, with it on and off.

**Results.**

- *The kink's excess is gone.* Each end's own excess, per unit $kh_R$: 11.1, 11.3, 11.1,
  10.6 and 9.3% at 7.5–30 Hz read across the mirror; 1.9, 2.4, 2.6, 2.5 and 2.5% with the
  slope (at 15 Hz +0.72 to +0.78% per end instead of +3.27 to +3.29%). Four points left 1.0%
  at 7.5 Hz but 5.9% at 30.
- *Lamb's problem* (order 8; misfit at 150, 240, 330 m and the diagonal; the largest swing):

  | $h$ | mirror | swing | slope | swing |
  |---|---|---|---|---|
  | 10 m | 15.5 / 21.4 / 28.2 / 22.6% | +4.8 to +5.0% | 14.6 / 20.8 / 27.2 / 21.9% | −2.9 to −2.7% |
  | 7.5 m | 10.5 / 13.9 / 18.1 / 14.6% | +5.6 to +5.8% | 9.0 / 12.9 / 17.0 / 13.6% | −1.3 to −1.0% |
  | 5 m | 6.4 / 7.8 / 9.5 / 8.1% | +5.2 to +5.4% | 4.3 / 6.2 / 8.2 / 6.5% | +0.2 to +0.4% |
  | 3.75 m | 4.8 / 5.4 / 6.3 / 5.5% | +4.4 to +4.5% | 2.5 / 3.6 / 4.8 / 3.8% | +0.6 to +0.8% |
  | 2.5 m | 3.3 / 3.5 / 3.8 / 3.5% | +3.3% | 1.3 / 1.7 / 2.3 / 1.8% | +0.7 to +0.8% |

  What is left is mostly the Rayleigh wave's speed (the lateness is unchanged). At 10 and
  7.5 m the swing is now too small: the negative second-order error that the kink had
  cancelled shows.
- *P doubling* (Taylor coefficients, 10, 7.5 and 5 m): 2.53, 2.36 and 2.26% mirrored; 2.86,
  2.38 and 2.27% with the slope. The same limit, the reference's 2.2%: this spherical wave
  spreads sideways at the epicentre, so its $v_z$ has a slope there (about 0.75% over the
  3.4 m the mirror read from), and the mirror's error had partly cancelled the plane-wave
  reference's. At the default (optimised, 10 m) 2.52% becomes 2.85%.
- *Sources stay the receivers' transpose*: reciprocity up to $1.1 \times 10^{-6}$ of the peak
  at order 8 ($7.3 \times 10^{-7}$ mirrored; $1.35 \times 10^{-6}$ at order 4), the
  dot-product test $1.6 \times 10^{-8}$ on the nodes and $9.6 \times 10^{-8}$ between them.
  The larger weights add single-precision round-off.
- *The gradient checks with it on*: still exact, but two need other steps. The walls-and-water
  $v_p$ check's direction curves more: extrapolated $4.8 \times 10^{-3}$, $2.9 \times 10^{-4}$,
  $1.2 \times 10^{-5}$ and $3.5 \times 10^{-5}$ at $h$ = 3.2e-2, 1.6e-2, 8e-3 and 4e-3 (the
  $h^4$ of Richardson's extrapolation, then round-off), so 1.6e-2 fails its 1.5e-4. And the
  time-shift check's step of 1e-3 crosses a lag change (1.6e-3), 5e-4 doesn't (2.1e-4).

**Decision.** Opt-in, `[boundary] surface_slope = true`, off by default; every baseline is
unchanged with it off. `the_surface_slope_keeps_sources_the_transpose_of_receivers` and
`the_surface_slope_takes_the_kink_out_of_lambs_problem` test it on. For surface data, a
land survey's especially, it removes the free surface's first-order amplitude error; it isn't
the default because it trades that for an error at the top of the band, reaches six nodes
down (a near-surface interface within them spoils the slope) and would move the baselines
above. Making the boundary condition's slope work instead would fix all three, at the cost of
operators that depend on the model: the gradient would need their derivative (or the ratio
pinned at the starting model, as $\Delta t$ is). The Rayleigh wave's speed is the free
surface's remaining first-order concern for misfits.

**Reproduce:** `just accuracy the_free_surface_excess_by_frequency_and_depth` (about seven
minutes); `just physics surface_slope`; Lamb's finer grids with `surface_slope = true` in
`survey::config`'s `[boundary]`.

## 2026-11-23 — Marmousi2: the whole section, muted L2

**Question.** Every Marmousi2 inversion so far was of the central 6 km, judged against its own
truth (2026-11-21, "Marmousi2: the mute line without the truth, and a held-out part"). With
gauss's 62 GiB the whole section fits. Does muted L2, the safest variant on the central part,
behave the same on the whole section?

**Method** (`examples/marmousi2/whole.toml`, observed traces from `whole-true.toml`). The
whole section, 17 km by 3.5 km, extruded 5 km along $y$: 681 × 201 × 141 nodes at 25 m
(central: 241 × 41 × 141). Nine explosive shots at 50 m depth from $x$ = 0.5 to 16.5 km every
2 km, 1.2 Hz Ricker; 641 receivers recording $v_x$ and $v_z$ on the seafloor (475 m) from 0.5
to 16.5 km; 3.5 s at $\Delta t$ = 2 ms; order 8, pinned as the central part's. Five L2
iterations for $v_p$ and $v_s$, density as it is, the top 500 m fixed, with the slow seafloor
waves muted. **The settings are the central part's, reused, not tuned on the whole
section:** the start (the section smoothed over 250 m, `25m-smooth`), the mute (keep what
arrives before $|\Delta x| / 600\ \text{m/s} + 0.87$ s, a 0.42 s taper), and L2 itself, the
variant chosen on the central part. The rule from the start's slowest seafloor $v_s$
(`slow_waves_mute`) was not worked out for the whole section. Run on gauss with the binary
built for its CPU, 12 threads.

The error ratio and cosine with the true update are measured as for the parts
(`whole_by_region`, which reads the result): over the free nodes below 500 m within four
nodes of the receivers' plane, the strip the parts' side layers leave them (`ACROSS`, now
applied to every case, which changes nothing for the parts). Over the whole 5 km across,
most of which the survey barely sees, every number came out within 1% of 1. The same
study splits the result into 1 km windows along $x$ and sets each beside the start's $v_s$
under the seafloor against the truth's, as `seafloor_speeds` does for the parts.

**Results.** Misfit (the muted L2 misfit, relative to the start's 0.7308; 23,033,010
unknowns):

| Iteration | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| Misfit left | 96.2% | 83.1% | 75.1% | 69.9% | 66.2% |
| Largest change | 2.0% | 10.2% | 11.3% | 12.8% | 14.3% |

Every step was the full one at the first trial. Error ratio and cosine by depth band
(iteration 1's model was overwritten before it could be kept), with the parts' muted L2 after
five iterations:

| | | 500–750 m | 750–1000 m | 1000–1500 m | 1500–2500 m | 2500–3125 m |
|---|---|---|---|---|---|---|
| 2 | $v_p$ | 0.967, +0.41 | 0.992, +0.19 | 0.998, +0.14 | 1.000, −0.01 | 1.000, −0.20 |
| | $v_s$ | 0.993, +0.15 | 1.006, −0.19 | 1.002, −0.19 | 1.000, −0.04 | 1.000, +0.17 |
| 3 | $v_p$ | 0.951, +0.41 | 0.988, +0.20 | 0.996, +0.15 | 1.000, −0.00 | 1.001, −0.20 |
| | $v_s$ | 0.988, +0.17 | 1.009, −0.18 | 1.004, −0.18 | 1.000, −0.04 | 1.000, +0.17 |
| 4 | $v_p$ | 0.942, +0.40 | 0.986, +0.20 | 0.995, +0.16 | 1.000, +0.00 | 1.001, −0.20 |
| | $v_s$ | 0.986, +0.18 | 1.011, −0.18 | 1.004, −0.18 | 1.000, −0.04 | 1.000, +0.17 |
| 5 | $v_p$ | 0.937, +0.40 | 0.984, +0.20 | 0.995, +0.16 | 1.000, +0.01 | 1.001, −0.20 |
| | $v_s$ | 0.984, +0.19 | 1.013, −0.18 | 1.005, −0.18 | 1.000, −0.05 | 1.000, +0.17 |
| central, 5 | $v_p$ | 0.916, +0.44 | 0.971, +0.35 | 0.999, +0.05 | 1.000, −0.03 | 1.001, −0.29 |
| | $v_s$ | 0.970, +0.26 | 1.004, −0.05 | 1.001, −0.05 | 1.000, +0.04 | 1.000, +0.13 |
| west, 5 | $v_p$ | 1.167, +0.13 | 1.055, +0.20 | 1.013, +0.17 | 0.997, +0.10 | 1.001, −0.20 |
| | $v_s$ | 0.791, +0.61 | 0.978, +0.26 | 0.999, +0.04 | 1.001, −0.17 | 0.999, +0.29 |

(The west's line was its own 690 m/s, the others' 600 m/s.) By region, after iteration 5:
the seafloor's $v_s$ (true 345 m/s everywhere), the harmonic mean over the 250 m below it,
true / start, the delay over 1 km, mean (largest), then the error ratio and cosine, and for
750–1000 m the window's share of the band's added squared error
$\sum (m - m_\text{true})^2 - (m_0 - m_\text{true})^2$:

| $x$ (km) | start's $v_s$ | 250 m below | delay (s) | $v_p$ 500–750 m | $v_s$ 500–750 m | $v_s$ 750–1000 m | share |
|---|---|---|---|---|---|---|---|
| 0.5–1.5 | 445 | 427 / 450 | 0.12 (0.13) | 0.988, +0.16 | 0.913, +0.59 | 1.002, −0.05 | +0.1% |
| 1.5–2.5 | 447 | 429 / 453 | 0.12 (0.13) | 0.990, +0.14 | 0.949, +0.53 | 0.997, +0.15 | −0.2% |
| 2.5–3.5 | 442 | 430 / 447 | 0.09 (0.12) | 1.013, −0.02 | 0.945, +0.54 | 0.991, +0.36 | −0.7% |
| 3.5–4.5 | 438 | 429 / 443 | 0.07 (0.08) | 0.998, +0.07 | 0.943, +0.53 | 0.995, +0.17 | −0.3% |
| 4.5–5.5 | 443 | 437 / 447 | 0.05 (0.09) | 1.004, −0.02 | 0.950, +0.51 | 1.004, −0.22 | +0.4% |
| 5.5–6.5 | 428 | 417 / 434 | 0.09 (0.17) | 0.992, +0.13 | 0.962, +0.45 | 0.993, +0.28 | −0.6% |
| 6.5–7.5 | 409 | 388 / 420 | 0.20 (0.28) | 0.931, +0.52 | 0.971, +0.35 | 0.992, +0.30 | −0.7% |
| 7.5–8.5 | 420 | 388 / 440 | 0.30 (0.48) | 0.932, +0.62 | 0.992, +0.18 | 0.998, +0.13 | −1.1% |
| 8.5–9.5 | 483 | 403 / 516 | 0.54 (0.80) | 0.874, +0.70 | 0.984, +0.27 | 1.004, −0.15 | +3.2% |
| 9.5–10.5 | 633 | 567 / 660 | 0.25 (0.77) | 0.971, +0.24 | 0.935, +0.37 | 0.978, +0.33 | −5.2% |
| 10.5–11.5 | 448 | 382 / 483 | 0.55 (0.68) | 0.864, +0.75 | 1.005, −0.13 | 1.012, −0.41 | +13.3% |
| 11.5–12.5 | 472 | 382 / 518 | 0.69 (0.75) | 0.873, +0.71 | 0.995, +0.17 | 1.043, −0.84 | +35.8% |
| 12.5–13.5 | 582 | 498 / 615 | 0.38 (0.78) | 0.992, +0.14 | 0.980, +0.21 | 1.032, −0.38 | +21.6% |
| 13.5–14.5 | 540 | 433 / 573 | 0.56 (0.64) | 0.935, +0.43 | 1.006, −0.10 | 1.031, −0.63 | +45.3% |
| 14.5–15.5 | 607 | 555 / 628 | 0.21 (0.41) | 1.008, +0.04 | 0.949, +0.33 | 0.997, +0.08 | −3.4% |
| 15.5–16.5 | 559 | 444 / 598 | 0.58 (0.79) | 0.908, +0.45 | 1.007, +0.02 | 0.983, +0.19 | −8.3% |

Half a period of the 1.2 Hz Ricker is 0.42 s.

- *Over the section, $v_p$ does nearly as well as on the central part, and $v_s$ worse.* At
  500–750 m $v_p$ ends at 0.937 (central 0.916, west 1.167); its update is a slower band just
  under the seafloor, about 200 m/s at the 99.9th percentile of the centre plane (775 m/s at
  most). It gains most east of 6.5 km, 0.864–0.935 in seven of the ten windows there with
  cosines of +0.43 to +0.75, including the windows where $v_s$ drifts. $v_s$ gains 1.6% at the
  top (central 3.0%) and moves away below 750 m, 1.013 and 1.005 with cosines of −0.18
  (central 1.004 and 1.001, −0.05), a little more each iteration.
- *The $v_s$ drift is where the start is cycle-skipped.* In the four windows from 10.5 to
  14.5 km the start's $v_s$ in the 250 m under the seafloor is 23–36% fast and the slow waves
  fall 0.38–0.69 s behind per km on average, up to 0.78 s, past half a period. There
  $v_s$ at 750–1000 m ends at 1.012–1.043 with cosines of −0.38 to −0.84, and those four
  windows hold 116% of the band's added squared error; 8.5–9.5 km adds 3.2% and the other
  windows together take about 20% off. In the west, 0.5–6.5 km, the start is 2–6% fast and the
  delay 0.05–0.12 s:
  there $v_s$ at 500–750 m improves by 4–9% (0.913–0.962, cosines +0.45 to +0.59) and
  barely moves below, as in the western part's own muted run (0.791 at the top, from five
  iterations of its data alone). The windows with a mean delay of 0.20–0.30 s (6.5–8.5,
  9.5–10.5, 14.5–15.5 km) gain too. 15.5–16.5 km, at 0.58 s, doesn't fit: $v_s$ at
  750–1000 m gains there (0.983), and it gets worse at 1000–1500 m instead (1.022).
- *Below 1.5 km nothing moves*, as on the parts.
- The misfit fell further than the central part's muted misfit (to 66.2% of the start's,
  against 0.332 to 0.240 from the first iteration to the fifth there), but the two are of
  different data and not comparable.

*Cost.* The three observed shots missing after an interrupted run took 13 min 17 s (`run
--shots` on a copy of `whole-true.toml` with only those sources, moved to `shot-006`–`008`;
1.7 GiB). The inversion took 25 h 8 min at 10.6 cores on average: the starting misfit and
gradient 3 h 39 min (24.3 min per shot's gradient), then 4 h 18 min per iteration (a trial
of nine forward runs, about 40 min, and the next gradient). Its peak was 48.3 GiB resident
against the plan's 48 GiB (the inversion's own count, 47.7 GiB); 22 GiB with `reconstruct =
true`. The CPU stayed at 62–71 °C and 3.5 GHz throughout, with no swapping.
`whole_by_region` takes a second.

**Decision.** The whole section repeats, within one survey, what the two parts showed
between them (2026-11-23, "Marmousi2: the held-out western part"): where the start predicts
the slow waves within half a cycle, $v_s$ improves; where it doesn't, $v_s$ below 750 m moves
away. One mute line for the whole section can't serve both. A faster line or the per-stage
mute (this entry's first plan) would take from the west what they take from the east, so
neither is worth 25 hours on gauss yet. The lever is a mute decided from the data by
offset, and here also by source: the slow waves' phase against the start's synthetics,
shot by shot and offset by offset (the literature note's "Keep or mute the slow waves").
Nine shots along 16 km give it the regions to tell apart. Try it first on the two parts,
then on the whole section. Until then this run is the whole section's baseline. Its
settings came from the central part, not this section. The delay uses the truth, so it
explains the result but can't choose a mute.

**Reproduce** (on gauss, or a machine with 50 GiB): `vibe-rations run --shots
examples/marmousi2/whole-true.toml` (about 40 min), then `vibe-rations invert
examples/marmousi2/whole.toml` (25 hours on gauss), then `just in marmousi2-study
whole_by_region`, which reads `out/marmousi2/whole` and any `out/marmousi2/whole-it<n>` kept
beside it (copies of the model after each iteration, which the inversion overwrites).

## 2026-11-24 — The inversion against Born's Gauss–Newton solution

**Question.** The gradient matches the exact kernels (2026-11-21, "The gradient against the
exact kernels"), but does an inversion's *model* follow? FWI has no closed form; its
linearisation about a homogeneous full space has. How does our L-BFGS inversion compare with
the analytic Born / Gauss–Newton solution of the same problem [@tarantola1984; @pratt1998;
@fichtner2011], and which of the nonlinearity, the solver's error and the parametrisation
separates them?

**Method** (`tests/analytic_inversion.rs`). A homogeneous full space ($v_p$ = 3000 m/s,
$v_s$ = 1732 m/s, $\rho$ = 2500 kg/m³) at $h$ = 15 m (4.6 points per shortest S wavelength),
77³ nodes with 12 of CPML, the default scheme, $\Delta t$ and the coefficients' band pinned at
the start's, 0.75 s recorded (402 samples of 1.87 ms). Three forces of $10^9$ N (Ricker,
10 Hz) 300 m from the centre on the $-x$, $-y$ and $-z$ axes, each across its line to the
centre; eleven three-component velocity receivers, 300 m out on the $+x$, $+y$, $+z$ axes and
at the corners of a cube ±195 m. The truth is a relative change of $v_s$ alone, a Gaussian
($\sigma$ = 45 m, a quarter of a dominant S wavelength) sampled at the nodes of a coarse
lattice 45 m apart, 5³ = 125 of them, and trilinearly interpolated between them (support
±135 m): the unknowns are those 125 coefficients, which represent the truth exactly.

The Jacobian $J$ (39 798 samples by 125) is built column by column: a node's change of
$v_s$ ($\delta\mu = 2\rho v_s^2\,\delta$, $\delta\lambda = -2\delta\mu$) scatters as a point
moment source whose moment is $\delta c_{ijkl}\,\partial_l u_k$ of the forward field, so its
velocity at a receiver is $h^3 \delta\mu \sum_{ij} \partial_j G_{ni} * S_{ij}$ with
$S_{ij} = \partial_i u_j + \partial_j u_i - 2\delta_{ij}\nabla\cdot\mathbf u$: Stokes's
solution and its spatial derivatives (`common::stokes`) for $S_{ij}$ tabulated every 0.25 ms
with three derivatives and its integral. Each hat's column sums these over the nodes (λ and
μ's normal part) and the shear positions (μ's shear parts, a quarter of each of four nodes),
as the staggered grid samples the model (`J_grid`), or over the nodes alone with every term
there (`J_point`, the continuous Born integral by the nodes' quadrature). Check
(`the_sensitivities_are_the_exact_kernels`, 2 s): $J^\mathsf T r\,\Delta t$ for a 2% faster
medium's residual equals the adjoint kernels integrated against a hat
(`Problem::grid_integral`, `Problem::kernels`), derived the other way round, to
$6 \times 10^{-10}$. $J^\mathsf T J$'s eigenvalues span $3.8 \times 10^3$ (median
$10^{-2}$ of the largest), computed by Jacobi rotations; the damped solution is
$(J^\mathsf T J + \mu\lambda_{\max} I)^{-1} J^\mathsf T \delta d$.

Data: $\delta d$ from Born ($J_\text{grid}$ times the truth) and from the solver (the
truth's traces less the start's). Runs at $+\varepsilon$ and $-\varepsilon$ split the
solver's difference into its odd part (Born's plus $O(\varepsilon^3)$, so what differs from
Born is the solver's own error) and its even part (the nonlinearity's second order). The
inversions, L2 unbalanced (so all minimise $\tfrac12\lVert d - d_\text{obs}\rVert^2$), from
the homogeneous start against the solver's data at +3%, 8 iterations each:

- *on the hats*: the library's `Lbfgs` and `line_search` on the 125 coefficients, the
  solver's gradient (`Inversion::evaluate`) taken onto the hats, no preconditioner;
- *the inversion*: `Inversion::iterate` as it is, every node outside the layers free,
  divided by the source side's illumination.

Each model is scored over the 17³ nodes of the support: the error ratio
$\lVert m - m_\text{true}\rVert / \lVert m_0 - m_\text{true}\rVert$ and the cosine of
$m - m_0$ with the true change; the misfit left as $\chi/\chi_0$, from solver runs.

**Results.** *Born against the solver's data differences:*

| Peak change | $\lVert\delta d\rVert$ / traces' | Born against the solver's | cosine | odd part against Born | even part / $\lVert\delta d\rVert$ | $J_\text{point}$ against $J_\text{grid}$ |
|---|---|---|---|---|---|---|
| +3% | 1.81% | 4.91% | 0.99919 | 0.29% | 5.06% | 0.83% |
| +10% | 5.67% | 16.3% | 0.99169 | 2.57% | 17.9% | 0.83% |

Born is off by the nonlinearity, about 1.6 times the peak change; the solver's own error in
the scattered waves is 0.29% (at 10% the odd part's 2.6% is mostly $O(\varepsilon^3)$). The
grid's sampling of the model moves Born data by 0.83%.

*Gauss–Newton* (error ratio, cosine; $\chi/\chi_0$ left in the solver's runs):

| Data, $J$ | $\mu$ = 0 | $10^{-4}$ | $10^{-3}$ | $10^{-2}$ | $\chi/\chi_0$ at 0 |
|---|---|---|---|---|---|
| Born, +3% or +10% | 0.00% | 0.25% | 1.95% | 10.4% (0.9955) | — |
| solver's, odd part, +3% | 0.76% (0.99997) | 0.84% | 2.32% | 10.7% | — |
| solver's, +3% | 15.3% (0.98857) | 14.8% | 12.2% | 9.9% (0.99637) | 0.20% |
| solver's, +3%, $J_\text{point}$ | 14.6% (0.98935) | 14.2% | 11.9% | 10.4% | 0.20% |
| solver's, odd part, +10% | 7.6% (0.99725) | 7.5% | 7.6% | 13.4% | — |
| solver's, +10% | 46.4% (0.90665) | 44.9% | 38.2% | 21.8% (0.97669) | 2.05% |

($10^{-8}$ and $10^{-6}$ are as 0.) Born data come back exactly undamped: the problem is well
posed on these unknowns. The solver's error alone costs 0.76% of the model; the nonlinearity,
5% of the data at +3%, costs 15%, amplified by the poorly determined directions (damping
helps only by shrinking them: 9.9% at $10^{-2}$). Whether $J$ samples the model as the grid
does changes little (14.6 against 15.3%).

*The gradient and the first step.* The solver's gradient at the start against
$-J^\mathsf T\delta d$ at the support's nodes: 0.01%, cosine 1.00000 (on the hats the same).
72% of its norm lies in the support; the illumination preconditioner ranges from 189 to 473
over the support (1000 at the water level) and takes that share to 64%. The inversion's first
step against $J^\mathsf T\delta d$ in the support: cosine 0.98796; against the
preconditioned $J^\mathsf T\delta d$, 1.00000. So the first step differs from the analytic
gradient by the preconditioner alone, which varies by a factor 2.5 across the target.
Its largest change, 2.00% (`initial_step`), lies 15 m from a source or receiver, where the
kernels are singular.

*Inversions* (misfit left, error ratio, cosine; on the hats, the error against Gauss–Newton
on the solver's data, $\mu = 10^{-4}$; the inversion, the error over the whole grid):

| Iteration | On the hats: $\chi/\chi_0$ | error | cosine | against GN | Inversion: $\chi/\chi_0$ | error | cosine | whole grid |
|---|---|---|---|---|---|---|---|---|
| 1 | 87.0% | 84.6% | 0.908 | 88.7% | 18.0% | 71.4% | 0.869 | 83.0% |
| 2 | 2.64% | 34.2% | 0.940 | 41.1% | 5.13% | 62.6% | 0.883 | 79.1% |
| 3 | 1.19% | 30.7% | 0.955 | 37.7% | 3.02% | 59.0% | 0.893 | 77.3% |
| 4 | 0.55% | 26.8% | 0.966 | 34.2% | 1.99% | 55.9% | 0.907 | 75.2% |
| 5 | 0.39% | 23.9% | 0.972 | 31.6% | 1.20% | 51.5% | 0.926 | 72.3% |
| 6 | 0.23% | 19.2% | 0.982 | 27.8% | 0.99% | 47.2% | 0.942 | 69.8% |
| 7 | 0.09% | 13.2% | 0.991 | 23.2% | 0.69% | 47.2% | 0.943 | 69.3% |
| 8 | 0.05% | 8.6% | 0.996 | 19.3% | 0.59% | 46.7% | 0.945 | 68.6% |

Every line search took the full step at its first try. On the hats, the first step is held
to 2% (`initial_step`) and barely moves the misfit; then L-BFGS passes the Gauss–Newton
solution (15.3%, $\chi/\chi_0$ 0.20%) by the sixth iteration and reaches 8.6% and 0.05%: it
minimises the nonlinear misfit, which Gauss–Newton only linearises once, so the truth is
nearer the nonlinear minimum than Gauss–Newton's step. The inversion with every node free
fits the data almost as well (0.59%) with the support still 47% wrong and the whole grid
69%: the extra unknowns let it explain the data by changes outside the target (from the
first step on: 77% of its norm lay outside the support, its largest change next to a
source or receiver).

**Decision.** No change to the code. The inversion's gradient is the analytic Born
sensitivity to $10^{-4}$, and its first step is $J^\mathsf T\delta d$ times the preconditioner,
exactly; what separates its model from the truth is not the solver (0.3% of the scattered
data, 0.8% of a Gauss–Newton model) but the nonlinearity, which even at a +3% Gaussian
costs a single Gauss–Newton step 15% and which iterating removes, and above all the
parametrisation: on 125 well-posed unknowns L-BFGS recovers the target to 8.6% in 8
iterations, on every node to 47% with a smaller misfit than Gauss–Newton's. A model's error
is a question of what the data determine, and a misfit that falls says little about it
[@fichtner2011resolution]. This is a reason to give the inversion a coarser or regularised
parametrisation (a basis, or a region), not yet a decision.

**Reproduce:** `cargo test --release --test analytic_inversion -- --ignored --nocapture
--test-threads 1`: `born_and_gauss_newton` (7 minutes on gauss's six cores, $J$ in about
a minute) and `lbfgs_against_gauss_newton` (23 minutes); the check, in `just test`.

## 2026-11-24 — Marmousi2: an inversion with its wavefields rebuilt, against checkpoints

**Question.** An inversion can rebuild each gradient run's forward wavefield backwards
instead of replaying it from checkpoints (`[inversion] reconstruct = true`; 2026-11-21, "An
inversion with its forward wavefields rebuilt"). On a small survey the two agreed to
$10^{-7}$–$10^{-6}$, but checkpoints stayed the default until an inversion of a realistic size
had been compared both ways, since water drifts five times more than rock. Over five
iterations on Marmousi2's central part, does the rebuilt run reach the same model?

**Method** (`tests/marmousi2_study.rs`: `l2`, then `l2_rebuilt`, which compares its model
with `l2`'s in `compare_variants`). Five L2 iterations on the central 6 km from the start
smoothed over 250 m, as the 2026-11-21 `l2`, once with checkpoints and once with the
wavefields rebuilt, everything else equal. Built from main at `ecadcb4` (`just mac-test
marmousi2_study`) and run natively on the maintainer's Mac (M1 Max, 32 GiB), one after the
other, single test thread, with desktop apps running alongside.

**Results.**

| | checkpoints (`l2`) | rebuilt (`l2_rebuilt`) |
|---|---|---|
| Misfit as the stage measures it, iterations 1–5 | 0.64184, 0.55273, 0.48380, 0.42256, 0.38212 | the same, to every printed digit |
| Largest change per step | 2.0, 10.1, 11.2, 11.9, 12.8% | the same |
| L2 left after the stage | 57.9% | 57.9% |
| Error ratio and cosine by band | 0.900 / 0.973 at 500–750 m, 0.994 / 1.095, 1.020 / 1.070 ($v_p$ / $v_s$) | the same, to every printed digit |
| Memory (the plan) | 4.41 GiB | 2.44 GiB |
| Wall time | 1,748 s | 2,032 s (+16%) |

The rebuilt model against the checkpointed one, by depth band, $\lVert m_b - m_a\rVert /
\lVert m_a - m_0\rVert$ and the largest $\lvert m_b - m_a\rvert / m_a$:

| | 500–750 m | 750–1000 m | 1000–1500 m | 1500–2500 m | 2500–3125 m |
|---|---|---|---|---|---|
| $v_p$ | $6.6 \times 10^{-6}$, $3.1 \times 10^{-6}$ | $1.1 \times 10^{-5}$, $1.2 \times 10^{-6}$ | $9.2 \times 10^{-6}$, $8.1 \times 10^{-7}$ | $9.9 \times 10^{-6}$, $2.0 \times 10^{-7}$ | $3.0 \times 10^{-4}$, $1.1 \times 10^{-4}$ |
| $v_s$ | $3.6 \times 10^{-6}$, $2.6 \times 10^{-6}$ | $3.2 \times 10^{-6}$, $1.8 \times 10^{-6}$ | $2.7 \times 10^{-6}$, $1.2 \times 10^{-6}$ | $4.0 \times 10^{-6}$, $3.1 \times 10^{-7}$ | $7.4 \times 10^{-5}$, $1.8 \times 10^{-5}$ |

- *The rebuilt inversion is the checkpointed one* to about $10^{-5}$ of the update down to
  2.5 km, and pointwise to $3.1 \times 10^{-6}$ of the model: single-precision drift, as on
  the small survey, and far below anything the error ratios resolve. The water's larger drift
  doesn't accumulate over five iterations.
- *The deepest band differs most* ($3.0 \times 10^{-4}$ of an update that is itself tiny
  there, error ratios 1.003 and 0.999; pointwise $1.1 \times 10^{-4}$). It lies just above
  the bottom absorbing layer, so drift from the band around the layers is the likely cause;
  not traced.
- *The `l2` run reproduces the 2026-11-21 one to every printed digit*, so neither the defaults
  changed since (band-limited layers and exact positions, which don't apply to models from
  files with sources and receivers on nodes) nor macOS against Linux moves this case.
- It took 16% longer, one run each on a desktop with other work, so indicative only.

**Decision.** For inversions, rebuilding the wavefields is as good as checkpoints at this
size: use `reconstruct = true` wherever memory binds. The whole section, for example, needs
22 GiB rebuilt against 48 with checkpoints, so it fits two at once on gauss, or on a 64 GiB
machine. Making it the default trades roughly 16% more time for 45% less memory. That is the
maintainer's call, not made here. The open question is closed for the central part; the
whole section rebuilt (22 GiB) wasn't run, since the Mac had 13 GiB free.

**Reproduce:** `just in marmousi2-study l2` and then `l2_rebuilt` (about 40 minutes each in
the dev container; 29 and 34 natively on the Mac: `just in mac-test marmousi2_study`, then
`out/mac/marmousi2_study l2 --ignored --exact --nocapture --test-threads 1` and the same for
`l2_rebuilt`).

## 2026-11-24 — Inversions rebuild their wavefields by default

**Question.** Rebuilding the forward wavefields backwards reached the checkpointed model to
about $10^{-5}$ of the update over five iterations on Marmousi2's central part, in 45% less
memory and about 16% more time (previous entry). Should it be the default for inversions,
and does switching move any baseline?

**Method.** `[inversion] reconstruct` now defaults to `true` (`false` replays checkpoints, bit
for bit as before); `gradient` is unchanged and still replays unless given `--reconstruct`.
The unit tests that compared the two ways now name checkpoints explicitly
(`rebuilding_the_wavefield_inverts_alike`, the memory plan's variants in
`counts_what_an_inversion_allocates`), and the Marmousi2 study `l2_rebuilt` became
`l2_checkpointed`, the other way round. `just physics` on main (`3c68b63`) and on the change,
both on gauss, compared line by line with timings set aside; then `just ci`.

**Results.** The 143 lines `just physics` prints are identical, including the sphere
inversion (17.9% of its misfit and 91.9% of the $v_s$ anomaly left after 3 iterations, now
with the wavefields rebuilt) and the inversion's finite-difference check ($4.9 \times
10^{-6}$ at $h = 10^{-2}$). `just ci` passes.

**Decision.** The maintainer's: inversions rebuild by default. The baselines in `AGENTS.md`
stand. A plan made with `memory --invert` now shows the rebuilt peak first and the
checkpointed one second ("with reconstruct = false"). The whole section's inversion needs 22
GiB by default, 48 with checkpoints.

**Reproduce:** `just physics` before and after, compared with
`diff <(grep -v '^test ' before.log | sort) <(grep -v '^test ' after.log | sort)` once
timings are masked; `just memory examples/marmousi2/central.toml --invert`.

## 2026-11-25 — Marmousi2: the slow-wave phase check

**Question.** The slow seafloor waves help $v_s$ where the start predicts them within half a
cycle and hurt it where it doesn't: on the western part L2 alone does best, on the central
part muting them does (2026-11-23, "Marmousi2: the held-out western part"), and the whole
section shows both region by region (2026-11-23, "Marmousi2: the whole section, muted L2").
Can the observed data and the start's synthetics alone, without the truth, decide where to
keep them, and does one such rule give the better variant on both parts? The test is the one
the literature note proposes ("Keep or mute the slow waves", "A test for both parts"): the
phase difference at a low frequency after the slow waves' line, continuous along offset
[@shah2012], within half a cycle [@fichtner2011].

**Method.**

- *The check* (`tests/marmousi2_study.rs`: `phase_check`, `phase_limits`,
  `continuous_limit`). For each part, the observed traces (the truth's runs, as the
  inversions make them) and the start's synthetics, five shots each. The window is the
  complement of the slow waves' mute: what it would remove after the line `slow_waves_mute`
  works out from each part's start and source (616 m/s on the central part, 690 on the
  western, from 0.87 s with a 0.42 s taper), tapered over the record's last 0.42 s. That line,
  not the 600 m/s one the central runs used, because the rule must come from the start alone,
  and the window should hold exactly what the mute removes. For each trace and component
  $c$, $Z_c = \sum_m w_m d_m e^{-i\omega t_m}$; the phase difference is
  $\Delta\varphi = \arg \sum_c w_c Z^\text{syn}_c \bar Z^\text{obs}_c$, with $w_c$ the
  components' balancing weights from the observed traces. That sum is the cross term of the
  balanced L2 misfit at that frequency, so both components count as the misfit weighs them;
  positive means the observed arrivals lag. Along each side of each shot, from the first
  receiver beside the source (25 m; the source's own, where $v_x$ vanishes, is left out)
  outwards, the phase is unwrapped, and the limit is the offset where it reaches $\pm\pi$
  (interpolated), or the last before a step of more than $\pi/2$ between neighbours 25 m
  apart (a drift of even 0.8 s/km moves it 0.15 rad per receiver at 1.2 Hz), or 25 m if
  it already exceeds $\pi/2$ there. Past the offset where the line meets the record's end
  (about 1.6 km on the central part, 1.8 on the western) the window is empty: nothing to judge
  and nothing for the mute to remove. Frequencies: 0.6 and 1.0 Hz (the stages' corners),
  1.2 Hz (the source's peak, standing for the unfiltered traces), and 2.4 Hz (where the
  Ricker's amplitude is a fifth of its peak's). Each component alone is printed too.
- *The check of the check* (`predicted_limits`, uses the truth): the phase $2\pi f\tau$ a
  wave along the seafloor gathers from the shot to each receiver, $\tau$ the sum over the
  columns between of the harmonic mean slowness in the 250 m under the seafloor, true less
  start's (as `seafloor_speeds`), and where it reaches $\pi$.
- *The misfit change* (`[[inversion.mute]]` and `[[inversion.stage.mute]]`: `shot`, `spare`,
  `spare_taper`; `misfit::mute_weight`, `misfit::for_shot`). A mute can act on one shot only
  and spare signed offsets $[a, b]$ (negative towards $-x$; `misfit::offsets` now carries the
  side as its sign, and the line goes by the distance), coming in beyond them over
  `spare_taper` metres with a half cosine: with $s(x)$ rising from 0 to 1, the weight is
  $1 - s(x)\,(1 - w(t, x))$. It is still a fixed weight by time and receiver, decided before
  the stage and held through it, so the misfit stays smooth and its adjoint exact: the
  pair $\langle FMx, y\rangle = \langle x, MFy\rangle$ holds to $10^{-12}$ with spared
  offsets (`spared_offsets_keep_the_transpose`), the weights are checked
  (`mutes_spare_offsets_and_shots`), and without the new keys they are computed as before,
  bit for bit. A mute for another shot weighs 1 rather than being dropped, so every shot
  holds what the memory plan counts. Per shot and per side, because the check's limits
  differ that much (below: 300 to 1,500 m on the central part at 1.2 Hz, the two sides of
  one shot 300 and 1,500 m).
- *The rule* (`phase_mutes`). For each shot with a side cut short, one mute for that shot:
  the slow waves' line, sparing each side out to its limit less the offset taper, which is
  the distance the line moves over its time taper ($616 \times 0.42 = 259$ m; 290 m in the
  west), so that what it keeps whole lies within $\pm\pi$. A side passed to its window's end
  is spared entirely; a shot passed on both sides gets no mute.
- *The inversions*, as `muted_l2` and `west_muted_l2`: five L2 iterations from the start
  smoothed over 250 m, with the rule's mutes from the 1.2 Hz limits (`central_phase_checked_l2`,
  `west_phase_checked_l2`); the frequency stages (0.6 Hz twice, 1.0 Hz twice, all once) with
  each stage's mutes from its own frequency's limits (`central_phase_checked_stages`,
  `west_phase_checked_stages`); and, after the first results, L2 with the 2.4 Hz limits
  (`central_phase_checked_l2_at_2f0`, `west_phase_checked_l2_at_2f0`). The limits are from
  the start's synthetics, once, not re-measured at each stage's start. The wavefields were
  rebuilt backwards (`reconstruct = true`, the default since 2026-11-24, which agrees with
  checkpoints to about $10^{-5}$ of the update). Run on gauss in the dev container, one at a
  time.

**Results: the check.** Limits by shot (source at 1–5 km) and side, left / right, in metres,
for both components together; "all" is every offset whose window holds anything; π: the
phase reached it; jump: a step over $\pi/2$; then what the truth's delays predict, capped at
the window's end.

| part | f (Hz) | shot 0 | shot 1 | shot 2 | shot 3 | shot 4 |
|---|---|---|---|---|---|---|
| central | 0.6 | all / 1100 π | all / 1175 jump | 937 π / 1188 π | 1109 π / all | 1484 π / all |
| | truth | 500 / 1600 | 1500 / 1600 | 1600 / 1600 | 1600 / 1500 | 1600 / 500 |
| | 1.0 | all / 769 π | all / 749 π | 641 π / 627 π | 1074 π / all | 1484 π / all |
| | truth | 500 / 1600 | 1500 / 1178 | 1500 / 1360 | 1104 / 1191 | 1509 / 500 |
| | 1.2 | all / 764 π | all / 300 jump | 609 π / 580 π | 1021 π / 513 π | 391 π / all |
| | truth | 500 / 1591 | 1500 / 1072 | 1201 / 813 | 926 / 1050 | 1197 / 500 |
| | 2.4 | all / 275 jump | all / 575 jump | 350 jump / 550 jump | 550 jump / 487 π | 934 π / all |
| | truth | 500 / 924 | 885 / 650 | 529 / 288 | 640 / 618 | 440 / 348 |
| west | 0.6 | all / all | all / all | all / all | all / all | all / all |
| | 1.0, 1.2 | all / all | all / 1700 jump | all / all | 750 jump / all | all / all |
| | 2.4 | all / all | all / 1700 jump | all / all | 675 jump / all | all / all |
| | truth | all at every frequency, except 1737 m on shot 0's right at 2.4 Hz | | | | |

The windows reach 500 m left of shot 0 and right of shot 4, 1.5 km on the other short sides,
1.6 km (central) and 1.8 km (west) elsewhere. The check took 124–127 s per part (five forward
runs of the start and the spectra), and 8.4 minutes with the truth's runs.

- *It tells the parts apart.* In the west it passes 8 of 10 sides to the window's end at
  every frequency, as the truth predicts; in the centre it cuts 6 of 10 sides short at 1.2
  Hz (300–1,021 m) and 4 even at 0.6 Hz (937–1,484 m), and the truth predicts cuts at 1.0 Hz
  and above.
- *It is conservative against the truth's delays in the centre at 1.2 Hz*: 300–1,021 m where
  the delays give 813–1,591 m. At 1.0 Hz on shots 3 and 4 they agree within 30 m (1,074 against
  1,104, 1,484 against 1,509 m). The window holds whatever arrives after the line, not the
  slow waves alone, and the delays are a 1D estimate of the slow waves only, so neither is
  exact; the literature note's linear estimate (1.1 km at 1.2 Hz, 2.2 km at 0.6) is of the
  same order.
- *Two false cuts in the west*: shot 1's right side at 1.7 km, at the window's end where the
  window is shorter than a period, and shot 3's left at 750 m (the phase goes +0.68 at 500 m,
  +2.08 at 750 m, then jumps), where the truth predicts no cycle skip; the phase there is of
  another arrival, or of interference, not traced.
- *Each component alone* is noisier: $v_x$ alone fails at once on several sides ("start"),
  because near the source it is small and its phase meaningless; summed with the balancing
  weights, $v_z$ carries the near offsets and $v_x$ the far.

**Results: the inversions.** Error ratio and cosine with the true update after five
iterations ($v_p$ / $v_s$ error ratio by band; full rows in the logs):

| | part | L2 left | 500–750 m | 750–1000 m | 1000–1500 m |
|---|---|---|---|---|---|
| L2 | central | 57.9% | 0.900 / 0.973 | 0.994 / 1.095 | 1.020 / 1.070 |
| | west | 17.3% | 1.227 / 0.714 | 0.988 / 1.018 | 0.986 / 0.935 |
| L2, muted | central | 85.7% | 0.916 / 0.970 | 0.971 / 1.004 | 0.999 / 1.001 |
| | west | 41.3% | 1.167 / 0.791 | 1.055 / 0.978 | 1.013 / 0.999 |
| L2, the check's mutes at 1.2 Hz | central | 64.3% | 0.907 / 0.968 | 0.982 / 1.062 | 1.007 / 1.037 |
| | west | 17.4% | 1.226 / 0.714 | 0.988 / 1.013 | 0.985 / 0.938 |
| L2, the check's mutes at 2.4 Hz | central | 63.7% | 0.903 / 0.966 | 0.980 / 1.059 | 1.011 / 1.036 |
| stages (central: unmuted; west: muted after 0.6 Hz) | central | 93.2% | 0.976 / 0.946 | 1.001 / 0.987 | 1.020 / 1.044 |
| | west | 64.0% | 1.006 / 0.911 | 0.987 / 0.993 | 0.992 / 0.983 |
| stages, the check's mutes per stage | central | 92.9% | 0.975 / 0.949 | 1.000 / 0.988 | 1.017 / 1.034 |
| | west | 64.9% | 1.010 / 0.934 | 0.989 / 0.989 | 0.991 / 0.978 |

Cosines, central / west, for the check's L2 at 1.2 Hz: $v_p$ +0.42, +0.19, −0.08 / +0.05,
+0.26, +0.19; $v_s$ +0.26, −0.11, −0.16 / +0.70, +0.31, +0.35. Over every free node $v_p$
1.000 / 1.002, $v_s$ 1.002 / 0.998. Every step was the full one at the first trial. Wall
times: 5,296 s (central L2), 6,498 s (west L2, with other work on gauss for most of it),
6,756 s and 6,943 s (the stages), and 6,986 s (central L2 at 2.4 Hz, gauss shared with other work); the central part's earlier runs took about
2,800 s on the Mac's container, so gauss runs these about twice as long (its 6 cores, and
rebuilding the wavefields).

- *In the west the rule keeps the slow waves, and so does what L2 does there*: it mutes two
  sides of 10, beyond 1.4 km and 460 m, and the result is L2's to within 0.005 (1.226 / 0.714,
  0.988 / 1.013, 0.985 / 0.938). The plain mute's losses ($v_p$ 1.167 at the top, 1.055 at
  750–1000 m; $v_s$ 0.791) don't happen.
- *In the centre it goes part of the way from L2 to the mute.* $v_s$ below 750 m moves away
  less than with L2 (1.062 and 1.037 against 1.095 and 1.070) but still away (cosines −0.11,
  −0.16), where the plain mute holds it (1.004, 1.001); $v_p$ at 750–1000 m gains a little
  (0.982, between 0.994 and 0.971). So the slow waves the check keeps, within $\pm\pi$ at
  1.2 Hz out to 300–1,000 m, still turn $v_s$ the wrong way below 750 m. The full-band misfit
  also fits frequencies above 1.2 Hz, for which those offsets are already beyond half a cycle:
  at 2.4 Hz the check's limits fall to 275–934 m on the cut sides.
- *With the stages the mutes change little*: on the central part the 0.6 and 1.0 Hz limits keep
  most of what the stages see, and the result is the unmuted stages' within 0.01; in the west
  (two sides muted from 1.0 Hz on) it lies near the stages muted after 0.6 Hz, $v_s$ a little
  worse at the top (0.934 against 0.911) and a little better below.
- *Cutting at 2.4 Hz instead changes nothing that matters* on the central part (spares of
  16–675 m on the cut sides instead of 41–762 m): $v_s$ 1.059 and 1.036 below 750 m, against
  1.062 and 1.037; $v_p$ 0.903, 0.980, 1.011. So it isn't that the 1.2 Hz limits were too
  generous: what still turns $v_s$ away is in the slow waves the check passes at every
  frequency, on the sides it keeps whole (shot 1's left to 1.5 km, shots 0's left and 4's right
  to 500 m) or at the shortest offsets, or it isn't in the slow waves. Not split further.
  The western run at 2.4 Hz (`west_phase_checked_l2_at_2f0`) was stopped unfinished: its mutes
  differ from the 1.2 Hz run's on one side by 75 m (spare 385 against 460 m), and with the
  central result unchanged it couldn't decide anything.

**Decision.** The check does what it was built for: from the data and the start alone it
tells the western part, where the start predicts the slow waves within half a cycle, from the
central part, where it doesn't, shot by shot and side by side, and its limits are of the order
the truth's delays give (shorter in the centre at 1.2 Hz). As a rule for the mute it is safe:
on neither part is the result worse than L2's by more than 0.007 in any band, and in the west it keeps what made L2 the best
variant there, which the plain mute lost. But it doesn't give the better variant on both
parts: in the centre it recovers only a third to a half of what the plain mute gained for $v_s$
below 750 m (1.062 / 1.037 against L2's 1.095 / 1.070 and the mute's 1.004 / 1.001), at 1.2 or
2.4 Hz alike. Within half a cycle is necessary for the slow waves to help, not sufficient:
in the centre the slow waves the check passes still push $v_s$ below 750 m the wrong way,
plausibly because the start's shallow $v_s$ is 18% fast and an L2 update of the slow waves
puts part of that error deeper, which no phase test of the data can see. So no mute default
still; the check is worth keeping as a QC (eight minutes) and as the way to place a mute
where one is wanted. The whole section (25 hours) is not yet worth running with this rule:
it would keep the west's gain where the start is good and recover a third to a half of the loss
where it isn't. Next, cheaper: split the central gradient by the check's mutes (as
`misfit_by_mute` did by the line), to find which kept slow waves turn $v_s$ away; and a
better shallow start (slowness averaged, as the literature note suggests), which would make
the check pass more often for the right reason.

**Reproduce:** `just in marmousi2-study phase_check` (8.4 minutes on gauss), then
`central_phase_checked_l2`, `west_phase_checked_l2`, `central_phase_checked_stages`,
`west_phase_checked_stages`, `central_phase_checked_l2_at_2f0` and
`west_phase_checked_l2_at_2f0` (90–115 minutes each on gauss; each runs the check first);
`cargo test --release --lib spare` for the mutes' unit tests.

## 2026-11-25 — Faster updates: AVX2 chosen at run time, bands of rows, one correlation pass

**Question.** The solver is said to be limited by memory bandwidth (gauss ran
`examples/layered.toml` at about 90 Mcells/s from the image, 104 built for its CPU). Where can
the time go without changing a bit of the results: compiling the loops for AVX2 at run time,
so that the default build still runs on any x86-64 CPU; reading each plane fewer times; and
the gradient's own passes?

**Method.** Three changes, each kept only if bit for bit:

1. *AVX2 at run time* (`src/simd.rs`, `solver::sweep`): the update loops and strain rates
   compiled twice, for the target and under `#[target_feature(enable = "avx2")]` without FMA,
   chosen by `is_x86_feature_detected!`; everything called in the loops `#[inline(always)]`,
   so the AVX2 copy compiles it with AVX2. Rust never contracts $a b + c$, so both copies round
   alike. `VIBE_RATIONS_SIMD=off` runs the plain copy; `run.json` and the set-up summary say
   which ran (`kernels`).
2. *Bands of rows* (`solver::blocks`, `block_size`): the updates split the interior into
   blocks of whole rows, four rows a band, each a whole column of planes worked through
   upwards, instead of one plane per task, so that of the $2L$ planes of three fields the
   $z$-derivatives read for one plane all but one were read for the previous plane a moment
   before. The tasks write disjoint rows of shared arrays through a small unsafe wrapper
   (`Shared`, `Outputs::row`), sound because the blocks tile the interior (a unit test checks
   every output and CPML memory value is handed out exactly once).
3. *The gradient's correlations in one pass* (`gradient.rs`): the stiffness terms, the density
   terms and the illumination read and wrote the 64-byte-per-node sums in three passes; the
   adjoint step leaves the adjoint strain alone, so one pass after it does all three, and the
   adjoint strain rate is added into the adjoint strain in place
   (`Solver::add_strain_rate`).

*Bit for bit.* At each commit, on gauss: `just test`; `just physics` with AVX2 and with
`VIBE_RATIONS_SIMD=off`, every line identical to main's (`f23c598`) once timings are masked;
and every output file of seven runs identical byte for byte to main's binary's (`run.json`
apart from its new `kernels` line): `examples/layered.toml` (free surface, snapshots), a
sphere survey `run --shots` and all at once, its `gradient` replayed and rebuilt
(`--reconstruct`), a two-iteration `invert`, and the central Marmousi2 gradient: 110 files,
with AVX2 and with the plain loops, on 12 threads, and the 99 besides the central gradient's
on 1 thread, either way (the same at the earlier commits, but for the plain loops on 1 thread
at the second). Unit tests run both copies and
compare every value: a run with sources, a free surface and layers with their sponge, its
strain rates and its steps back; a gradient, replayed and rebuilt; blocks of several shapes
against whole planes; the strain rate added in place against written first.

*Timing.* On gauss (Xeon E5-1650 v3, 6 cores, 12 threads, 15 MiB L3, 4 × DDR4-2133), idle: a
driver waited until nothing else ran before each run and discarded any run that something
else overlapped. Cases: `examples/layered.toml` as it is (129 × 129 × 81, 1,666 steps,
snapshots); Marmousi2's central part forward (`central-true.toml`, 241 × 41 × 141, 1,750
steps, five sources at once); the whole section's grid for 300 steps (`whole-true.toml` at
0.6 s, 681 × 201 × 141); and the central part's gradient (`central.toml` against
`central-true.toml`'s traces, checkpoints). Mcells/s are the run's own (interior nodes ×
steps / time stepping), the gradient's its wall time. Builds: the default (no `RUSTFLAGS`, as
the image) and the host's (`-C target-cpu=native`), of main and of each commit; three runs
each, medians with the range where it isn't a single value.

**Results.** Forward, 12 threads, Mcells/s:

| Build | `layered.toml` | central part | whole section |
|---|---|---|---|
| main, default | 91 | 102 | 85 (85–86) |
| main, native | 103 (103–104) | 118 | 108 |
| AVX2 at run time | 101 (101–102) | 117 (116–117) | 108 |
| the same, `VIBE_RATIONS_SIMD=off` | 91 | 102 (101–102) | 85 |
| and bands of rows (final), default | **111** | **117** (117–118) | **114** (114–115) |
| the same, `VIBE_RATIONS_SIMD=off` | 99 | 100 | 85 |
| the same, native | 112 | 117 | 114 (114–115) |

On 6 threads, one per core, `layered.toml` runs at 74 (main, default), 95 (main, native) and
106 (final, default or native). The central part's gradient, wall time: main 119.4 s
(default; 119.3–119.4) and 106.9 s (native); AVX2 at run time 107.3 s; and bands of rows
107.1 s (107.1–107.3); and one correlation pass **102.5 s** (102.1–102.5), plain 116.9 s,
native 99.7 s.

- *The default build gains 15–34% forward* (22% on `layered.toml`, 15% on the central part,
  34% on the whole section), 43% on 6 threads, and 14% on a gradient. Built for the CPU it
  is now within 1–3% of that, so a native build no longer pays.
- *AVX2 gives the native build's gain* (11–27%) with the default build: the copy compiled for
  AVX2 is what `target-cpu=native` was buying. So the default build was limited by its
  arithmetic as much as by memory.
- *The bands pay on the grids whose planes leave the cache*: with AVX2 the whole section
  108 → 114 and `layered.toml` 101 → 111; the central part, whose planes are small (about
  50 KB a field) and whose 41 rows make only 11 bands for 12 threads, 117 → 117, and its
  gradient 107.3 → 107.1 s. The whole section is
  now faster than `layered.toml` instead of slower. With the plain loops they gain on
  `layered.toml` (91 → 99) and nothing elsewhere (central 102 → 100, whole 85 → 85).
- *The correlation pass* takes 4.3% off a gradient (107.1 → 102.5 s).

What didn't pay (one run each unless marked; [x, y, z] the block's nodes, whole columns
unless said):

- *Parts of rows* (x-tiles): the whole section 114 (whole rows, bands of 4) against 101 and
  92 for rows cut in 2 and 3, and 89, 82, 77 and 65 cut in 4 to 11 with wider bands; the
  central part 117 against 91 for halves. Shorter rows defeat the hardware prefetchers, and
  every part of a row reads its $x$-halo again.
- *Other bands*: of 2 rows, 108, 111 and 112 on `layered.toml`, the central part and the
  whole section (three runs each) against 111, 117 and 114 for 4; of 6 rows on the whole
  section, 111. Wider bands in whole columns leave too few tasks for 12 threads: 8 rows give
  `layered.toml` 17 and ran at 82, 16 rows in two columns 86, 12 rows in two columns 104; in
  columns of 27 planes, 8 rows ran at 112 against 116 for 4.
- *Shorter columns*, for more tasks: `layered.toml` gains in columns of 27 planes (116
  against 111, three runs each), the central part loses (108 in columns of 47, 105 with bands
  of 3 rows in columns of 71, three runs each, against 117). Whole columns are kept: the
  bigger grids matter more, and a rule by grid would be tuned to three grids.
- Not tried: fusing the velocity and stress sweeps a few planes apart (the sources, the free
  surface's mirrors, the sponge and `step_back_with` make it delicate, and the measured case
  for it isn't there yet); anything that changes the layout, the precision or the rounding
  (FMA), which would change the results.

**Decision.** All three changes are in: results identical bit for bit, the default build as
fast as a native one, 15–34% faster forward and 14% on a gradient on gauss. The solver image
needs no change: built for baseline x86-64, it picks AVX2 itself. Keep every function called
in the loops `#[inline(always)]` (a call left out of line runs the plain code inside the
AVX2 copy, silently slower) and never use `mul_add` or enable FMA there. Apple silicon and
the browser run only the plain loops, with the bands; not measured on the Mac yet, where the
plain loops gained on `layered.toml` and lost 2% on the central part on gauss. gauss's
`/srv/data/bin/vibe-rations-native` was rebuilt from this branch (the old one kept as
`.before-perf`), though the default build now does as well.

**Reproduce.** On gauss, build main and each commit with `cargo build --release` and with
`RUSTFLAGS="-C target-cpu=native"`; then three times each, nothing else running:
`vibe-rations run examples/layered.toml`, `run examples/marmousi2/central-true.toml`, `run` a
copy of `whole-true.toml` with `duration = 0.6`, and `gradient examples/marmousi2/central.toml
--observed <central-true's run>`, with and without `VIBE_RATIONS_SIMD=off`, and with
`RAYON_NUM_THREADS=6`. For identity: `just physics` with and without `VIBE_RATIONS_SIMD=off`
against main's, timings masked, and `cmp` on every output file of the runs above.

## 2026-11-26 — Marmousi2: a start smoothed through its slownesses, the check's split, and a hard seabed

**Question.** In the centre the slow seafloor waves turn $v_s$ below 750 m the wrong way even
where the phase check keeps them (previous Marmousi2 entry). Three questions, each part-sized:
does a start smoothed through its slownesses, which keeps vertical traveltimes, bring the
central seafloor's $v_s$ within 10% of the truth's and stop that (Gate A of the day's plan);
which of the slow waves the check keeps turn $v_s$ away; and does the survey on the model the
published elastic inversions use, Poisson's ratio capped at 0.25 [@kohn2012; @jeong2012;
@cao2022], move $v_s$ towards the truth below 750 m (2026-11-23, "Marmousi2 in the
literature": if so, the soft seabed is what separates us from them).

**Method.**

- *The slowness start* (`vibe-rations extrude --smooth 250 --slowness`;
  `examples/marmousi2/convert.sh`'s `25m-*-slowness`). The solid nodes' $1/v_p$ and $1/v_s$
  are smoothed (the same Gaussian, $\sigma$ = 250 m, over the solid nodes) and inverted
  back: harmonic means instead of arithmetic ones, so a thin slow layer over faster rock
  doesn't come out fast. Density is smoothed as before. Gate A (`seafloor_speeds_by_start`):
  the harmonic mean of $v_s$ over the 250 m below the seafloor within 10% of the truth's on
  the central part.
- *The split* (`central_split_by_check`). The central gradient at the smooth start
  (`iterations = 0`), split as `misfit_by_mute` split it by the line, now by the check's
  mutes at 1.2 Hz (`phase_mutes`, the slow waves' line at 616 m/s): everything; what the
  check keeps; the same with the sides it passes whole muted too ("whole sides muted"); the
  check with every side's spare capped at 300 m ("within 300 m"); the plain line. The
  differences between sets are the gradients of the misfit's parts the bigger set keeps and
  the smaller removes (`split_by_mutes`). For each, by depth band, the norm of
  $m_0\,\partial\chi/\partial m$ and the cosine of $-\nabla\chi$ with the true update, and the
  same divided by the illumination.
- *The hard-seabed control* (`extrude --max-poisson 0.25`; the parts `central-hard` and
  `west-hard`). $v_s$ is raised to $v_p\sqrt{(1 - 2\nu)/(2 - 2\nu)} = v_p/\sqrt3$ on each fine
  sample before the averaging onto the nodes, in the truth and (before smoothing) in the
  start alike. Five L2 iterations as `l2`, on both parts (`hard_l2`, `west_hard_l2`), and on
  the central part with the 600 m/s mute (`hard_muted_l2`, as `muted_l2`), each judged
  against its own (hard) truth.
- Everything ran on gauss, one study at a time (inversions 4,490–4,590 s, of which about
  830 s set-up and recording; the split 5,274 s; the check 451 s).

**Results: the slowness start.**

- *Gate A fails, narrowly.* Under the central seafloor the start's harmonic mean $v_s$ is
  462 m/s against the truth's 414 (+11.4%; the smooth start +18.2%), and the slow waves'
  delay over 1 km falls from 0.37 s on average (0.80 s at most) to 0.25 s (0.60). The west
  was already close: +2.7% against +4.0%, 0.06 s against 0.09 s. The check's line from it
  is 604 m/s in the centre and 676 in the west (616 and 690 from the smooth start), and at
  1.2 Hz it still cuts a side of every central shot (limits of 300–1,475 m) and two western
  sides (muted from 1,416 m on shot 1's right and 466 m on shot 3's left), the same two as
  from the smooth start.
- *L2 from it fits more and turns $v_s$ away more sharply, and shallower*: 53.2% of the
  start's misfit left (57.9% from the smooth start); error ratio and cosine by band:

| start | | 500–750 m | 750–1000 m | 1000–1500 m | 1500–2500 m | 2500–3125 m |
|---|---|---|---|---|---|---|
| slownesses | $v_p$ | 0.940, +0.35 | 0.998, +0.08 | 0.997, +0.08 | 1.000, +0.04 | 1.002, −0.18 |
| | $v_s$ | 0.950, +0.32 | 1.147, −0.41 | 1.003, +0.04 | 0.997, +0.18 | 0.999, +0.20 |
| values (`l2`) | $v_p$ | 0.900, +0.44 | 0.994, +0.12 | 1.020, −0.15 | 1.001, −0.02 | 1.003, −0.24 |
| | $v_s$ | 0.973, +0.25 | 1.095, −0.15 | 1.070, −0.17 | 0.999, +0.04 | 0.999, +0.15 |

  The ratios are each against their own start's error, which is smaller for the slowness
  start, so they don't compare absolute errors; the cosines do compare directions. The
  wrong-way part has moved up into 750–1000 m and strengthened there (−0.41 against −0.15),
  and 1000–1500 m is now neutral. So a start that is less fast under the seafloor moves the
  problem rather than removing it. Since Gate A failed, the slowness start's check-muted runs
  were not run.

**Results: the split.** The central gradient at the smooth start (misfit 0.660; $v_x$ 76.7%,
$v_z$ 23.3%), norm and cosine (divided by the illumination: within 0.05 of these):

| kept | share of the misfit | $v_p$ 750–1000 m | $v_p$ 1000–1500 m | $v_s$ 750–1000 m | $v_s$ 1000–1500 m |
|---|---|---|---|---|---|
| everything | 100% | 1.79e-3, +0.04 | 9.43e-4, −0.16 | 9.53e-3, −0.10 | 5.72e-3, −0.19 |
| the check's mutes | 49.8% | 1.36e-3, +0.26 | 4.52e-4, −0.03 | 4.03e-3, −0.02 | 1.43e-3, −0.11 |
| … and whole sides muted | 36.7% | 1.30e-3, +0.32 | 4.03e-4, +0.01 | 2.02e-3, +0.08 | 5.43e-4, −0.09 |
| … spares capped at 300 m | 36.1% | 1.30e-3, +0.32 | 4.07e-4, +0.01 | 2.15e-3, +0.07 | 5.75e-4, −0.09 |
| the line | 34.3% | 1.30e-3, +0.31 | 3.83e-4, +0.02 | 1.46e-3, +0.09 | 5.09e-4, −0.09 |
| check less whole sides muted | 13.1% | 3.93e-4, −0.16 | 1.89e-4, −0.08 | 2.85e-3, −0.08 | 1.11e-3, −0.10 |
| check less 300 m cap | 13.6% | 3.93e-4, −0.16 | 1.85e-4, −0.08 | 2.75e-3, −0.07 | 1.07e-3, −0.10 |
| 300 m cap less the line | 1.9% | 2.80e-4, +0.05 | 1.45e-4, −0.04 | 1.56e-3, +0.01 | 4.25e-4, −0.01 |

- *What turns the model away is the slow waves beyond 300 m on the sides the check passes
  whole* (shot 1's left to 1.5 km, shots 0's left and 4's right to 500 m): 13% of the
  misfit, pointing away for $v_p$ and $v_s$ alike at 750–1500 m. Muting those sides, or
  capping every spare at 300 m, brings what the check keeps level with the line in every
  band (cosines within 0.03). The two ways differ by 0.6% of the misfit (the cut sides'
  spares beyond 300 m, up to 762 m, against the whole sides' first 300 m) and by at most
  0.01 in cosine.
- *The slow waves within 300 m are harmless*: the 1.9% the cap keeps and the line removes is
  neutral below 750 m (cosines within ±0.05).
- So the check passes, as within half a cycle at every frequency it tries, slow waves whose
  gradient still points the wrong way. A plausible reading, not tested: at 1.2 Hz a slow
  wave's sensitivity reaches a few hundred metres under the seafloor; it sees the start's
  sediments as too fast on average and asks for slower $v_s$ over its whole depth, including
  below 750 m, where the smooth start, averaging the soft sediments with the rock below, may
  already be too slow. That is a depth ambiguity, not cycle skipping, and no phase test
  can see it.

**Results: the hard seabed.** Five L2 iterations against the hard truth (the start's misfit
0.652 in the centre, 0.275 in the west; every step taken whole at the first trial):

| | L2 left | | 500–750 m | 750–1000 m | 1000–1500 m | 1500–2500 m | 2500–3125 m |
|---|---|---|---|---|---|---|---|
| central | 34.3% | $v_p$ | 0.870, +0.50 | 0.923, +0.39 | 1.012, +0.00 | 0.999, +0.04 | 1.004, −0.14 |
| | | $v_s$ | 0.815, +0.63 | 1.039, +0.26 | 1.014, +0.22 | 0.984, +0.18 | 0.995, +0.19 |
| central, muted | 57.2% | $v_p$ | 0.875, +0.49 | 0.922, +0.40 | 0.996, +0.10 | 1.000, +0.01 | 1.001, −0.04 |
| | | $v_s$ | 0.808, +0.59 | 0.945, +0.33 | 1.000, +0.09 | 0.994, +0.12 | 0.999, +0.09 |
| west | 22.6% | $v_p$ | 1.072, +0.03 | 1.057, +0.13 | 1.032, +0.19 | 1.010, −0.11 | 1.012, −0.45 |
| | | $v_s$ | 1.028, +0.31 | 1.523, +0.36 | 1.243, +0.38 | 1.005, +0.04 | 1.002, −0.05 |

- *Without the soft seabed, $v_s$ moves towards the truth at every depth in the centre*:
  cosines +0.26 and +0.22 below 750 m, against −0.15 and −0.17 on the real model, and L2
  fits more (34.3% left against 57.9%). The error ratios below 750 m are still a little
  above 1 (1.039, 1.014): the update points the right way but is too long for its alignment
  there. The mute helps on the hard model too (0.945 and 1.000), and costs nothing above
  750 m.
- *In the west the ratios below 750 m are large, 1.52 and 1.24, with positive cosines*
  (+0.36, +0.38): from $r^2 = 1 - 2s\cos\theta + s^2$, the update there is about 1.6 and 1.2
  times as long as the true update. The hard start's error there is probably small, which
  makes any misdirected update look large; not checked in absolute terms.
- So the answer to the literature entry's question is mostly yes: the soft seabed is what
  turns the central $v_s$ the wrong way below 750 m, and on the published model the same
  survey and inversion move it the right way, though not yet closer than the start below
  750 m without a mute.

**Decision.** Gate A failed and none of the day's runs meets G1 (a data-only start and mute
that keep central $v_s$ at 750–1500 m at or below 1.00 with the west within 0.01 of plain
L2), so Gate B says no whole-section run yet. Capping the spares at 300 m would reproduce the
line's gradient in the centre, but it is a fixed offset chosen here on the central part, and
in the west the slow waves beyond 300 m are what make L2 best; it isn't a recipe. The split
and the hard-seabed control point the same way: the trouble is a soft layer under the
seafloor whose error L2 spreads over the slow waves' whole depth. Next, layer stripping by
depth (Masoni et al. strip by offset and frequency instead [@masoni2016]): `[inversion]
fixed_below`, a new key, holds the model at and below a depth, so three L2 iterations with only the 250 m under
the seafloor free (`fixed_below = 750`) can take up the shallow error first, then five with
everything below 500 m free from that result (`shallow_first_l2`, `west_shallow_first_l2`,
`slowness_shallow_first_l2`; about two hours each on gauss). If the shallow pass absorbs the
error, the second should no longer push $v_s$ below 750 m away.

**Reproduce:** `examples/marmousi2/convert.sh` (writes the slowness and hard directories),
then `just in marmousi2-study <name>` for `seafloor_speeds_by_start` (12 s),
`slowness_phase_check` (7.5 minutes), `slowness_l2`, `hard_l2`, `west_hard_l2`,
`hard_muted_l2` (about 75 minutes each on gauss) and `central_split_by_check` (92 minutes);
`cargo test --release --lib section::`, `--test cli` and `--lib frees_only` for the new
options' tests.

## 2026-11-26 — Marmousi2: layer stripping by depth

**Question.** The previous entry read the central part's wrong-way $v_s$ update below 750 m as
L2 spreading the error of the too fast sediments under the seafloor over the slow waves' whole
depth. If the top 250 m under the seafloor is inverted first, alone, does it take up that
error, and does an inversion of everything from there then leave $v_s$ below 750 m alone
(G1: error ratio at most 1.00 at 750–1500 m in the centre, the west within 0.01 of plain
L2)?

**Method.** `[inversion] fixed_below` (new: the model is held at and below that depth, as
`fixed_above` holds it above): three L2 iterations with only 500–725 m free (`fixed_below =
750`, 10 of the 105 free depths, 37,620 unknowns), then five L2 iterations with everything
below 500 m free (395,010 unknowns) from that result, the misfit otherwise as `l2`
(`shallow_first`: `shallow_first_l2` and `west_shallow_first_l2` from the start smoothed over
250 m, `slowness_shallow_first_l2` from the one smoothed through its slownesses). Judged
after each part and, at the end, against the original start; the misfit reductions multiply
(the second part's is against its own start). On gauss, one at a time: 7,540–7,720 s each.
Masoni et al. strip by offset and frequency instead [@masoni2016]; this is the simplest
version by depth.

**Results.** Error ratio and cosine with the true update, against the original start, beside
five iterations of L2 alone (previous entries):

| | L2 left | | 500–750 m | 750–1000 m | 1000–1500 m | 1500–2500 m | 2500–3125 m |
|---|---|---|---|---|---|---|---|
| central, shallow 3 | 81.2% | $v_p$ | 0.936, +0.42 | 1 | 1 | 1 | 1 |
| | | $v_s$ | 0.990, +0.15 | 1 | 1 | 1 | 1 |
| central, shallow 3 then all 5 | 54.4% | $v_p$ | 0.891, +0.45 | 0.988, +0.16 | 1.018, −0.13 | 1.001, −0.02 | 1.003, −0.24 |
| | | $v_s$ | 0.975, +0.26 | 1.092, −0.17 | 1.066, −0.16 | 0.999, +0.06 | 0.999, +0.15 |
| central, L2 5 | 57.9% | $v_p$ | 0.900, +0.44 | 0.994, +0.12 | 1.020, −0.15 | 1.001, −0.02 | 1.003, −0.24 |
| | | $v_s$ | 0.973, +0.25 | 1.095, −0.15 | 1.070, −0.17 | 0.999, +0.04 | 0.999, +0.15 |
| west, shallow 3 | 45.5% | $v_p$ | 1.038, +0.06 | 1 | 1 | 1 | 1 |
| | | $v_s$ | 0.793, +0.65 | 1 | 1 | 1 | 1 |
| west, shallow 3 then all 5 | 11.0% | $v_p$ | 1.474, +0.04 | 1.133, +0.21 | 1.022, +0.21 | 0.998, +0.07 | 1.008, −0.33 |
| | | $v_s$ | 0.689, +0.73 | 1.070, +0.38 | 0.927, +0.38 | 0.999, +0.06 | 0.997, +0.29 |
| west, L2 5 | 17.3% | $v_p$ | 1.227, +0.05 | 0.988, +0.26 | 0.986, +0.19 | 0.999, +0.06 | 1.004, −0.33 |
| | | $v_s$ | 0.714, +0.70 | 1.018, +0.31 | 0.935, +0.36 | 1.000, +0.05 | 0.998, +0.28 |
| slowness start, shallow 3 | 81.0% | $v_p$ | 0.964, +0.33 | 1 | 1 | 1 | 1 |
| | | $v_s$ | 0.972, +0.24 | 1 | 1 | 1 | 1 |
| slowness start, shallow 3 then all 5 | 40.8% | $v_p$ | 0.931, +0.37 | 0.993, +0.13 | 0.995, +0.10 | 0.999, +0.05 | 1.002, −0.19 |
| | | $v_s$ | 0.936, +0.38 | 1.195, −0.40 | 1.003, +0.07 | 0.996, +0.17 | 0.999, +0.21 |
| slowness start, L2 5 | 53.2% | $v_p$ | 0.940, +0.35 | 0.998, +0.08 | 0.997, +0.08 | 1.000, +0.04 | 1.002, −0.18 |
| | | $v_s$ | 0.950, +0.32 | 1.147, −0.41 | 1.003, +0.04 | 0.997, +0.18 | 0.999, +0.20 |

- *In the centre the shallow part alone barely moves $v_s$* (0.990, cosine +0.15: an update a
  tenth to a fifth as long as the true one), though it fits a fifth of the misfit, mostly
  through $v_p$ (0.936). From the slowness start likewise (0.972). So inverted alone, the
  250 m under the seafloor doesn't take up its own $v_s$ error from these data.
- *The full inversion afterwards ends where L2 alone does*: within 0.01 of `l2` in every band
  in the centre ($v_s$ 1.092 and 1.066 below 750 m, against 1.095 and 1.070), with three more
  iterations and 54.4% of the misfit left against 57.9%. From the slowness start it turns
  $v_s$ at 750–1000 m further away (1.195 against 1.147).
- *In the west, where the shallow part alone does recover $v_s$* (0.793, +0.65), the result
  is mixed: $v_s$ better at the top and at 1000–1500 m (0.689, 0.927 against 0.714, 0.935) but
  worse at 750–1000 m (1.070 against 1.018), $v_p$ much worse above 1 km (1.474, 1.133
  against 1.227, 0.988), with 11.0% of the misfit left against 17.3%.

**Decision.** Layer stripping by depth fails G1 on both counts and is not the recipe: it
changes nothing in the centre and costs the west $v_p$. It also qualifies the previous
entry's reading. The test assumed that the shallow part, held alone, would take up its own
error; in the centre it doesn't (0.990), so the deeper update isn't simply the shallow error
put in the wrong place, which a depth constraint could correct. Where the start predicts the
slow waves within half a cycle (the west), the shallow part alone does recover its $v_s$.
That points back to the start: in the centre the data can't fix the sediments' $v_s$ from this start, with or without
depth constraints, and without the soft sediments (the hard-seabed control) the same
inversion moves the right way. The whole section stays unrun. Next, a better shallow start
from the data (the slow waves' dispersion inverted for $v_s$ under the seafloor, the
literature note's "Scholte-wave dispersion first"), or the soft-seabed remedies of the
literature (staging $v_p$ before $v_s$, $v_z$ before $v_x$); the plan's step 8 meanwhile.

**Reproduce:** `just in marmousi2-study shallow_first_l2` (and `west_shallow_first_l2`,
`slowness_shallow_first_l2`; about 2.1 hours each on gauss), after `convert.sh`;
`cargo test --release --lib frees_only` for `fixed_below`.

## 2026-11-27 — Rebuilding inside a surface

**Question.** Aaker et al. [@aaker2020gji] rebuild the forward wavefield inside a closed
surface one cell wide by injecting the particle velocities and tractions recorded on it: six
values per surface node and step, the velocity–stress form of Mittet's Kirchhoff-integral
reconstruction [@mittet1994], equivalent to Vasmel and Robertsson's multiple point sources
[@vasmel2016]. The band (2026-11-20) saves about $6L$ values per face node instead, so that
each step can be undone exactly. How do the two compare in the rebuilt field's and the
gradient's accuracy, in memory and in time? Aaker et al. name the surface's corners as the
largest error and warn against placing it near a fluid–solid contact.

**Method** (`src/surface.rs`, `Strategy::Surface`, `gradient --surface [--inset N]
[--segment N]`). *The sources, discretely.* The field truncated to the box $\Omega$, $w = Xu$
with $X$ the box's indicator ($\frac12$ on its faces), satisfies the solver's own equations
plus the commutators of $X$ with the derivatives, $B\,[X, D]\,\sigma$ and $C\,[X, D]\,v$,
which live only around the faces. In the continuum they are a force
$\mathbf n \cdot \sigma\,\delta_{\partial\Omega}$ and a stress source
$c : (\mathbf v \otimes \mathbf n)\,\delta_{\partial\Omega}$, which need only $\mathbf v$ and
$\mathbf n \cdot \sigma$ on the surface. On the grid the commutator of an order-$2L$ stencil
spreads over $2L$ nodes across the face, so the one-node surface stands in for it: each
recorded value is read onto the face's node plane with the stencil's own interpolation along
the normal (`Solver::stagger_weights`; each field keeps its own lattice along the face) and
spread back with its transpose, as Aaker et al. shift values with interpolators of the
stencil's length and optimisation. Each value's share is the product of the masks along the
three axes, so edges take $\frac14$ and corners $\frac18$ of a node-lattice value without
special cases. At second order the commutator has exactly this support, so the rebuild is
exact; at higher orders it leaks. *Time reversal.* With $\tilde v = v$ and
$\tilde\sigma = -\sigma$, a forward step of $(\tilde v, \tilde\sigma)$ undoes a step of
$(v, \sigma)$, if the sources are withdrawn (forces) or added (moments) instead, so the
rebuild runs the solver's ordinary step, layers, sponge and mirrors included; what leaks out
of $\Omega$ the layers absorb. Record $R_n$ is the state before forward step $n$: the
reversed step for $n$ takes its forces from $R_n$'s stresses and its stress sources (for step
$n - 1$'s stress update) from $R_n$'s velocities. A segment starts from its end state,
exact, truncated to $\Omega$, with the CPML memory cleared and a lone stress update, so it
needs $k + 1$ records. The physical sources are truncated by the same mask. *Along a free
surface* the box is open at the top: the mirrors carry the truncated field as they carry the
whole one, and on the surface plane the stress sources use the plane-stress moduli.
*Placement.* The surface sits at least $L$ nodes inside the layers (`--inset`, by default
$L$), so that nothing it injects reaches them; the gradient means something only $L + 1$
nodes inside it (`Bounds::region`), where the strain rates and material averages read
nothing on or outside the faces, so the inset and $L + 1$ nodes more (at least $2L + 1$)
inside each absorbing side are lost to it; and
the sources must sit well inside (a force two nodes from a face gave errors of 0.3). The
memory plan counts the records and operators exactly
(`counts_what_a_rebuilding_gradient_run_allocates`).

**Results: second order is exact** (`gradient::tests::surface_gradients_approach_the_replay`;
a $44 \times 40 \times 40$ grid, 10 m, two layers, a sphere, a force and an explosion at
18 Hz, 0.25 s). The largest difference from the replay's gradient, inside the region,
relative to each parameter's largest value: order 2 with Taylor coefficients
$4.5 \times 10^{-8}$–$1.4 \times 10^{-7}$, along a free surface and with layers all round,
default segments or the whole run: round-off, which checks the masks, the signs, the time
pairing, the sources' truncation and the free surface. Order 4 (optimised):
$3.4 \times 10^{-4}$–$1.6 \times 10^{-3}$; order 8: $1.1$–$4.0 \times 10^{-3}$ along a free
surface, $4.5 \times 10^{-4}$–$1.5 \times 10^{-3}$ all round.

**Results: the rebuilt field** (`just drift surface`, `tests/reconstruction.rs`; a
$64 \times 56 \times 48$ grid, 10 m, 8-node layers, two layers, a sphere, a force and an
explosion at 14 Hz, `optimised_for = 4`). The largest difference from the run it retraces,
at nodes one or more inside the faces, relative to each array's peak over the run, after
stepping back the given number of steps from the end. By run length (order 8, $C = 0.85$,
free surface, 2,926 steps):

| Steps back | 1 | 10 | 100 | 300 | 1000 | 2000 | 2926 |
|---|---|---|---|---|---|---|---|
| Band | $5.4 \times 10^{-12}$ | $5.4 \times 10^{-12}$ | $8.7 \times 10^{-12}$ | $1.5 \times 10^{-11}$ | $2.7 \times 10^{-11}$ | $8.3 \times 10^{-11}$ | $4.4 \times 10^{-7}$ |
| Surface | $1.1 \times 10^{-6}$ | $1.1 \times 10^{-6}$ | $6.7 \times 10^{-7}$ | $8.2 \times 10^{-7}$ | $6.2 \times 10^{-7}$ | $1.4 \times 10^{-5}$ | $2.5 \times 10^{-4}$ |

The largest over the whole rebuild: band $4.5 \times 10^{-7}$, surface $6.2 \times 10^{-3}$.
By forward step, with the field's own level:

| Step | 2000 | 1000 | 500 | 300 | 200 | 150 | 100 | 50 | 0 |
|---|---|---|---|---|---|---|---|---|---|
| Surface | $9.0 \times 10^{-7}$ | $1.2 \times 10^{-5}$ | $1.8 \times 10^{-4}$ | $4.3 \times 10^{-3}$ | $4.7 \times 10^{-3}$ | $1.5 \times 10^{-3}$ | $5.5 \times 10^{-4}$ | $3.0 \times 10^{-4}$ | $2.5 \times 10^{-4}$ |
| Field | $1.1 \times 10^{-5}$ | $1.4 \times 10^{-4}$ | $4.3 \times 10^{-3}$ | $6.9 \times 10^{-2}$ | $0.14$ | $0.13$ | $0.84$ | $0.29$ | 0 |

From step 300 on, while the field crosses the faces, the error is 3–9% of the field's level
at that step: what the surface leaks is a share of what crosses it, not round-off. Earlier,
with the field still around the sources, what is left of it is $2.5$–$5.5 \times 10^{-4}$,
since the leakage leaves through the layers. By how far
inside the surface the nodes lie, all the way back: $6.2 \times 10^{-3}$ (one node or more),
$1.4 \times 10^{-3}$ (two), $4.6 \times 10^{-4}$ (four), $3.1 \times 10^{-4}$ (six),
$2.7 \times 10^{-4}$ (nine). By order, 1.2 s, all the way back (band; surface; values a step):
order 4 $5.6 \times 10^{-7}$, $3.3 \times 10^{-3}$, 45,200; 8 $8.7 \times 10^{-7}$,
$6.4 \times 10^{-3}$, 38,064; 12 $9.1 \times 10^{-7}$, $7.9 \times 10^{-3}$, 31,504; 16
$6.4 \times 10^{-7}$, $9.1 \times 10^{-3}$, 25,520 (the box shrinks as the layers' inset
grows). Higher orders leak more, as the commutator the surface stands in for widens. At order
8, the surface six and eight nodes in from the layers rather than four: $6.7 \times 10^{-3}$
and $7.4 \times 10^{-3}$, no better. By the wavelets' frequency, order 8, 1.2 s
(`surface_drift_with_frequency`): 14, 10 and 7 Hz, 3.1, 4.4 and 6.3 points per S wavelength
at $2.5 f_0$, give $6.4$, $6.2$ and $5.2 \times 10^{-3}$ all the way back (the band
$0.9$–$1.7 \times 10^{-6}$). So the leak is not a dispersion error, which would fall fast as
the waves are better sampled: it is the shape of what is injected, one plane read and spread
with interpolation weights, against the commutator's own shape over $2L$ nodes, wrong by a
fixed share of the field at the face whatever its wavelength, and most of it stays within a
few nodes of the faces. By boundary, order 8, 976 steps (100 and 300 steps
back, then all the way; the band all the way):

| Boundary | 100 | 300 | All the way | Band |
|---|---|---|---|---|
| Layers, free surface | $1.6 \times 10^{-5}$ | $7.3 \times 10^{-5}$ | $6.4 \times 10^{-3}$ | $8.7 \times 10^{-7}$ |
| Layers all round | $1.8 \times 10^{-7}$ | $3.4 \times 10^{-7}$ | $4.8 \times 10^{-3}$ | $9.8 \times 10^{-7}$ |
| Rigid walls, free surface | $4.5 \times 10^{-3}$ | $4.9 \times 10^{-3}$ | $6.9 \times 10^{-3}$ | $1.1 \times 10^{-6}$ |
| Layers, 60 m of water at the surface | $5.6 \times 10^{-6}$ | $2.5 \times 10^{-5}$ | $1.2 \times 10^{-2}$ | $2.5 \times 10^{-6}$ |
| All round, seabed 60 m (top face 6 nodes below it, in rock) | $1.1 \times 10^{-6}$ | $8.2 \times 10^{-6}$ | $5.0 \times 10^{-3}$ | $7.3 \times 10^{-7}$ |
| All round, seabed 120 m (on the top face) | $3.3 \times 10^{-6}$ | $6.2 \times 10^{-6}$ | $8.9 \times 10^{-3}$ | $1.2 \times 10^{-6}$ |
| All round, seabed 150 m (top face 3 nodes above it, in water) | $7.2 \times 10^{-4}$ | $8.0 \times 10^{-4}$ | $4.8 \times 10^{-3}$ | $1.3 \times 10^{-6}$ |
| All round, seabed 200 m (top face 8 nodes above it, in water) | $1.2 \times 10^{-3}$ | $2.0 \times 10^{-3}$ | $7.8 \times 10^{-3}$ | $3.4 \times 10^{-6}$ |

(With layers all round the force was moved down to 260 m, 14 nodes below the top face; the
first run, with it two nodes from the face, gave $2.5 \times 10^{-2}$.) Closing the top adds
four edges and four corners and changes little, so with the masks' shares the corners are
not the main error here. With rigid walls the leakage comes back: the late error is
$4.5 \times 10^{-3}$ instead of $1.6 \times 10^{-5}$, so the method needs absorbing layers
around the surface. Water: with the seabed crossing the side faces (water at the top under
the free surface) the worst error doubles; with the top face on the seabed it is the largest
measured; with the top face in the water, the late field (100–300 steps back) is
$10^{-3}$, thousands of times the solid's; with the face in rock six nodes below the seabed, no
different from no water. Aaker et al.'s advice holds: keep the faces in the solid, clear of
the seabed.

**Results: the gradient.** The rebuilt gradients against the replay's at the 20,460 nodes of
the region (order 8, free surface), with the memory and time of each (4 cores):

| | 976 steps | | | | | 2,926 steps | | | | |
|---|---|---|---|---|---|---|---|---|---|---|
| | $v_p$ | $v_s$ | $\rho$ | memory | time | $v_p$ | $v_s$ | $\rho$ | memory | time |
| Checkpoints | | | | 447 MiB | 13.0 s | | | | 725 MiB | 38.4 s |
| Band, default segments (122; 209) | $1.5 \times 10^{-7}$ | $1.4 \times 10^{-7}$ | $1.7 \times 10^{-7}$ | 264 MiB | 15.4 s | $3.7 \times 10^{-8}$ | $1.2 \times 10^{-7}$ | $4.6 \times 10^{-8}$ | 411 MiB | 44.3 s |
| Band, the whole run | $1.3 \times 10^{-7}$ | $8.0 \times 10^{-7}$ | $4.6 \times 10^{-7}$ | 939 MiB | 13.3 s | $1.4 \times 10^{-7}$ | $8.5 \times 10^{-7}$ | $4.6 \times 10^{-7}$ | 2.60 GiB | 38.0 s |
| Surface, default segments (244; 488) | $4.0 \times 10^{-4}$ | $2.3 \times 10^{-3}$ | $1.7 \times 10^{-3}$ | 146 MiB | 15.8 s | $2.8 \times 10^{-4}$ | $4.2 \times 10^{-4}$ | $2.8 \times 10^{-4}$ | 205 MiB | 46.6 s |
| Surface, the whole run | $2.8 \times 10^{-4}$ | $4.2 \times 10^{-4}$ | $2.8 \times 10^{-4}$ | 218 MiB | 13.2 s | $2.8 \times 10^{-4}$ | $4.2 \times 10^{-4}$ | $2.8 \times 10^{-4}$ | 501 MiB | 38.2 s |

(Times from the first of two runs; the second was 5–20% slower throughout, with the same
differences.) With 60 m of water at the top, inside the box, which is open there, and its
sides crossing the seabed (976 steps): the band $9.9 \times 10^{-8}$–$1.1 \times 10^{-7}$ in
default segments and $1.8$–$7.8 \times 10^{-7}$ for the whole run; the surface
$1.6$–$2.6 \times 10^{-3}$ and $2.2$–$9.9 \times 10^{-4}$, up to four times the solid
model's.

The surface's gradient is $3 \times 10^{-4}$–$2 \times 10^{-3}$ off, a thousand times the
band's; shorter segments make it worse, not better as for the band, since each starts from
its end truncated while the field may be loud at the faces. On `examples/layered.toml` (the
2026-11-20 entry's test: the upper crust's $v_s$ lowered from 2600 to 2500 m/s; 1,666 steps;
16-node layers, so the surface at nodes 20 and 108 across and 60 down, and the region nodes
25–103 across and 0–55 down), on the container's 4 cores, with the plan's peak and the
process's largest resident set:

| Strategy | Memory | Resident | Time | $v_p$ | $v_s$ | $\rho$ |
|---|---|---|---|---|---|---|
| Checkpoints, every 71 | 3.66 GiB | 3.70 GiB | 128.8 s | | | |
| Band, segments of 209 | 1.70 GiB | 1.74 GiB | 148.0 s | $4.4 \times 10^{-6}$ | $1.2 \times 10^{-6}$ | $1.0 \times 10^{-6}$ |
| Band, the whole run | 5.85 GiB | 5.83 GiB | 127.1 s | $3.2 \times 10^{-6}$ | $9.3 \times 10^{-7}$ | $7.5 \times 10^{-7}$ |
| Surface, segments of 417 | 1.00 GiB | 1.04 GiB | 164.3 s | $1.0 \times 10^{-4}$ | $2.4 \times 10^{-5}$ | $2.2 \times 10^{-5}$ |
| Surface, the whole run | 1.59 GiB | 1.57 GiB | 159.4 s | $2.9 \times 10^{-5}$ | $8.1 \times 10^{-6}$ | $4.8 \times 10^{-6}$ |

The last three columns are the largest difference from the replay's gradient in the region,
relative to its largest value there (the band's are the same over all its nodes; the
surface's, over all of them, $8.2 \times 10^{-2}$ for $v_p$, since near and outside the faces
it means nothing). The misfits agree to every printed digit. A surface record is 176,266
values against the band's 862,968, a fifth, and the surface's gradient needs 41% less memory
than the band's with segments; but the solvers and the correlation sums, which neither
changes, are most of it. Here the surface is 6–23 times the band's error, far less than on
the small model. Not because its waves are better sampled, which the frequency study rules
out; the small model's sources sit 13–14 nodes from the faces, this one's 20–44, but what sets
the difference is not settled. It
is the slowest: its rebuild is an ordinary step, CPML included, which the band's step back
skips (the layers hold 55% of these nodes), and the injection runs on one thread; 24–28%
over checkpoints against the band's $-1$ to $+15$%. On the same model grown to
$1001 \times 1001 \times 401$ nodes (`memory --gradient`): checkpoints 803 GiB, the band 247
GiB in segments and 478 GiB for the whole run, the surface 174 and 204 GiB.

**Decision.** Keep the band as the inversion's reconstruction; offer the surface as
`gradient --surface`, for comparisons, and not to `invert`. The surface saves memory (a
gradient run needs 41% less than with the band on `examples/layered.toml`, 30% less on a
grid 300 times larger), but its
error is the injection's shape, not round-off: $10^{-5}$–$10^{-3}$ of the gradient's largest
value against the band's $10^{-7}$–$10^{-6}$, growing with the order, largest within a few
nodes of the faces and near a seabed, and the same at any sampling. It is also slower than
the band (its rebuild runs the layers) and constrains the survey: the gradient is lost over
the inset and $L + 1$ nodes more (at least $2L + 1$) inside every absorbing side, the
sources must sit well inside, and the faces must stay in rock clear of the seabed, which a
marine survey can't arrange: along a free surface its sides cross the water and the seabed
(which doubled the error here), and a top face under the seabed would lose the sediments
there, where Marmousi2's gradient is largest. The Marmousi2 configurations, 41 nodes across
with 16-node layers, leave 9 nodes inside, too few for a surface at all. Made exact, the
surface would need everything the stencil reads across each face, $L$ planes deep, which is
what the band records (about $6L$ values per face node): in this scheme the band is the
exact version of the method, and the one-node surface its approximation at a quarter to a
fifth of the records.

**Reproduce:** `just drift surface` (the field and the gradient study); `cargo test --release
--lib surface` (second order); `just memory <config> --gradient --surface [--segment 100000]`;
for `examples/layered.toml`, as in the 2026-11-20 entry, `just run examples/layered.toml
--out out/surface/layered/observed`, then `just gradient guess.toml
out/surface/layered/observed [--reconstruct | --surface] [--segment 100000]` with `guess.toml`
the example with the upper crust's `vs = 2500.0`, and compare the `grad_*.npy` inside the
region.

## 2026-11-27 — The Hessian: its diagonal as a preconditioner, and exact products

**Question.** The maintainer's notes on Hessian kernels (after Fichtner and Trampert
[@fichtner2011hessian] and Fichtner's book [@fichtner2011]) suggest two things for the
inversion: (1) the Gauss–Newton Hessian's diagonal, per parameter and with both the source's
and the receivers' sides in it, as a preconditioner in place of the source-side illumination;
(2) exact Hessian–vector products, the second-order terms included, for Newton-type steps
[@metivier2013] and for resolution analysis. Does (1) help on the sphere survey, and can (2)
be made exact for the discrete run, as the gradient is?

**The Hessian probed** (`tests/hessian_study.rs`, no solver change). On the sphere survey's
own data (zero residual, so $H = J^T W J$), Hessian rows from central differences of the exact
gradient extrapolated in the step, for Gaussian bumps ($\sigma = 10$ m) of $v_p$, $v_s$ and
$\rho$ at the sphere's centre: the 3 × 3 block they span is symmetric to $1.3 \times 10^{-4}$;
$v_s$ and $\rho$ correlate $+0.77$, $v_p$ and $v_s$ $-0.22$, $v_p$ and $\rho$ $+0.06$; a $v_p$
bump's row lands 2.1 times as much on $v_s$ and 2.7 times on $\rho$ as on $v_p$; the rows'
centroids sit 20–30 m above the bump (towards the sources and receivers). The diagonal per
parameter from 16 random adjoint sources (below) is unbiased within its noise at the bump's
centre (±35% per node), and from 60 to 170 m depth it falls about 940-fold for $v_s$, 600 for
$\rho$ and 28 for $v_p$, while the source-side illumination falls 6.5-fold.

**The diagonal from probes** (`Misfit::probe`, `Inversion::gauss_newton_diagonal`). For L2,
$\chi = \frac12 \lVert \sqrt W F M (s - d) \rVert^2$ with $F$ the zero-phase filter, $M$ the
mutes and $W$ the balancing weights, so an adjoint source $z = (FM)^T \sqrt W \xi$ with
independent signs $\xi$ gives $E[(J^T z)^2] = \operatorname{diag}(J^T (FM)^T W FM J)$, the
Gauss–Newton diagonal exactly in expectation [@bekas2007]; the shots' draws are independent,
so one gradient run per shot and probe. A node's estimate scatters by about
$\sqrt{2/\text{probes}}$ of its value, so it can be smoothed within each parameter's free
nodes (`hessian_smoothing`, a Gaussian's $\sigma$ in metres, normalised by the smoothed
indicator so the region's edges keep their level). It replaces the illumination in the
preconditioner, $P = w_c / (e / e_{\max} + \text{water level})$, with $e_{\max}$ the largest
over all parameters, so the diagonal also sets the parameters' shares, or per parameter with
`hessian_balance = true`, so only its shape differs from the illumination's. It holds the
receivers' side already, so `receiver_runs` is refused with it. A unit test checks that
probes sending one sample each sum to the diagonal exactly and that their covariance is the
L2 Hessian in the samples (to $10^{-5}$ with random signs).

*On the sphere survey* (36 × 32 × 26 nodes at 10 m, two surface-near forces at 12 Hz, three
lines of three-component receivers on the surface, the top 60 m fixed, five L-BFGS
iterations; `just study gauss_newton_diagonal`), the error left in the sphere, relative to
its anomaly:

| Case | Preconditioner | L2 left | $v_s$ | $v_p$ | $\rho$ |
|---|---|---|---|---|---|
| $v_s$, $\rho$, $v_p$, denser core | illumination | 6.6% | 85.4% | 100.3% | 79.8% |
| | 8 probes | 7.4% | 101.2% | 106.7% | 86.9% |
| | 8 probes, 10 m | 7.5% | 102.4% | 108.5% | 89.0% |
| | 8 probes, 20 m | 7.2% | 99.9% | 108.3% | 89.2% |
| | 16 probes, 10 m | 7.8% | 103.1% | 107.8% | 89.7% |
| | 8 probes, 10 m, water level $10^{-2}$ | 5.6% | 95.4% | 103.9% | 78.6% |
| | 8 probes, 10 m, balanced | 6.5% | 102.3% | 100.6% | 87.3% |
| $v_p$, $v_s$ | illumination | 6.7% | 83.3% | 100.8% | |
| | 8 probes, 10 m | 9.2% | 99.9% | 114.4% | |
| | 8 probes, 10 m, balanced | 8.8% | 96.4% | 101.9% | |
| $v_p$, $v_s$, nothing fixed | sources' side | 14.1% | 92.6% | 100.5% | |
| | both sides, 8 receiver runs | 9.7% | 89.0% | 100.6% | |
| | 8 probes, 10 m | 5.9% | 91.1% | 103.4% | |

The diagonal is worse than the illumination in every case with the top fixed. Unbalanced, it
gives $v_p$, which the data barely constrain, a large share (its diagonal falls 28-fold with
depth against $v_s$'s 940), which turns into cross-talk ($v_p$ error above 100%), as the
equalised classes did in an earlier trial ("Scaling the parameters"). Balanced, $v_p$ is back
where the illumination leaves it, but $v_s$ still isn't recovered: the diagonal's shape,
falling hundreds of times faster with depth than the illumination, pushes the updates to the
sphere's deep side and below, where the misfit is barely sensitive. The number of probes and
the smoothing change little, so the noise isn't the cause; a higher water level, closer to
the illumination's flatness, helps. With nothing fixed it lowers the misfit most (5.9%
against 9.7% for both sides of the illumination) without recovering the sphere better.

**The Born run** (`src/born.rs`). The tangent-linear model of the discrete run: a second
solver beside the forward one, stepped by the solver's own step with the forward run's
velocity change and strain rate as sources, $(\delta B / B)(v^{n+\frac12} - v^{n-\frac12})$
and $\delta M\, e^{n+\frac12}$, through the forward-differentiated material averaging
(arithmetic for $\rho$, harmonic for $\mu$, the plane-stress moduli on the free surface). Its
sources are the gradient's correlations untransposed, so it is the gradient's exact
transpose. Against central differences of the traces (extrapolated, $\Delta t$ and the layers
held), for random changes of $v_p$, $v_s$ or $\rho$ at every node outside the layers, the
largest component's error is $3.5 \times 10^{-5}$, $1.1 \times 10^{-5}$, $8.1 \times 10^{-5}$
at order 4, $4.3 \times 10^{-5}$, $1.1 \times 10^{-5}$, $8.6 \times 10^{-5}$ at order 8, and
$6.3 \times 10^{-5}$, $3.0 \times 10^{-5}$, $1.6 \times 10^{-4}$ with rigid walls and water;
$\langle J \delta m, r \rangle = \langle \delta m, J^T r \rangle$ to $\le 1.6 \times 10^{-7}$
(`the_born_run_is_the_derivative_and_the_gradient_its_transpose`).

**Gauss–Newton products** (`Inversion::gauss_newton_product`): a Born run, its traces
weighed by the misfit's Hessian in the samples (`Misfit::gauss_newton_source`,
$(FM)^T W FM$, L2 only), and a gradient run driven by them. Unit vectors give the diagonal the
unit probes give, to $10^{-8}$ of its largest entry; against differences of the gradient at
zero residual, entries agree to $1.4 \times 10^{-5}$–$1.7 \times 10^{-4}$ and columns to
$1.5 \times 10^{-5}$–$5.5 \times 10^{-4}$. Products are symmetric to $1.5 \times 10^{-6}$ with
the forward wavefield replayed; rebuilt backwards, its drift near the source puts them
$1.4 \times 10^{-5}$ off ($2.3 \times 10^{-3}$ with the time error kept at $C = 0.3$), so
Hessian work replays checkpoints.

**Full products** (`src/hessian.rs`, `Inversion::hessian_product`). The gradient is
$g = T\,S(u, \tilde u)$, sums bilinear in the forward and adjoint wavefields taken to the nodes
by factors $T$ that depend on the model only through the free surface's plane-stress moduli
and the shear moduli's harmonic means, so
$H \delta m = T\,S(\delta u, \tilde u) + T\,S(u, \delta\tilde u) + \delta T\,S(u, \tilde u)$
(forward over reverse; Fichtner and Trampert's second-order adjoints, discretely). $\delta u$
is the Born run; $\delta\tilde u$ is a second adjoint run driven by the residuals' change (the
Gauss–Newton source of the Born traces) and by the adjoint wavefield scattered by the
material's change, as the Born run scatters the forward one: $(\delta B / B)$ times the
adjoint velocities' change over each step and $\delta M (\tilde\varepsilon^n -
\tilde\varepsilon^{n+1})$, the adjoint strain's, plus $\delta M \tilde\varepsilon^N$ at its
start. The residuals enter through the material too ($B$ for velocity, $M$ for pressure), and
these differences include what they added, so that is taken in. The conversion to $v_p$,
$v_s$, $\rho$ adds its own change applied to the gradient. The forward and Born runs are
replayed together from checkpoints; four runs per shot, and the gradient that comes with it
is `gradient`'s bit for bit. Against central differences of the balanced L2 gradient, at the
gradient checks' guessed sphere against data from another, every node outside the layers
changed (`the_hessian_product_is_the_gradients_derivative`):

| Model | Change of | Against the differences (extrapolated) | Gauss–Newton apart |
|---|---|---|---|
| order 8, free surface, two layers, sphere | $v_p$ | $1.6 \times 10^{-5}$ | $4.4 \times 10^{-2}$ |
| | $v_s$ | $3.6 \times 10^{-6}$ | $1.7 \times 10^{-2}$ |
| | $\rho$ | $1.8 \times 10^{-5}$ | $1.5 \times 10^{-2}$ |
| rigid walls, water, uneven spacing | $v_p$ | $9.5 \times 10^{-6}$ | $1.6 \times 10^{-2}$ |
| | $v_s$ | $8.2 \times 10^{-6}$ | $2.5 \times 10^{-2}$ |
| | $\rho$ | $6.3 \times 10^{-5}$ | $2.5 \times 10^{-2}$ |

(the largest of the three parameters' errors, relative to its size). The plain differences
are $7 \times 10^{-4}$–$5 \times 10^{-3}$ off at $h = 1.6 \times 10^{-2}$ and four times less at
$h/2$: the step's $h^2$. So the second-order terms are right to a few $10^{-3}$ of their size
or better. In the inversion's unknowns, at zero residual the full and Gauss–Newton products
are the same bit for bit; from a start 5% slow in $v_s$ and 3% light, $u^T H v$ and
$v^T H u$ agree to $1.3 \times 10^{-8}$ and $Hu$ is $1.9 \times 10^{-2}$ from $J^T W J u$.

**Against Stokes's full space** (`gauss_newton_products_against_born` in
`tests/analytic_inversion.rs`, 6 minutes). The Jacobian of the 2026-11-24 entry, built from
Stokes's solution for 125 trilinear hats of $v_s$ around the target, gives $J^T J c$; the
inversion's Gauss–Newton product of the hats' nodal values, taken back to the hats
($P^T H P c$), against it, for the truth's Gaussian, the centre's hat and two random sets:
$1.6 \times 10^{-5}$ to $2.8 \times 10^{-5}$ with $J$ sampling the model as the grid does
(the quadratic forms $c^T P^T H P c$ to $10^{-6}$–$2.5 \times 10^{-5}$), and 1.3–3.3% against
the continuous Born integral's ($J$ by the nodes' quadrature), about twice what that sampling
moves the Born data (0.83%). So in a full space the product is the analytic one to the
wavefields' error, and what is left is how the grid samples the model.

**Against Lamb's problem** (`hessian_against_lambs_problem` in `tests/accuracy.rs`, `just
accuracy hessian_against_lambs_problem`, 20 minutes). The survey of "The free surface against
Lamb's problem" (a vertical force on a Poisson solid, $v_z$ on the surface 150–330 m away),
with the half-space's speeds ($v_p$ and $v_s$ together, so it stays a Poisson solid) and its
density changed uniformly, every node and the layers' too. Pekeris's closed form gives the
traces' derivatives (five-point differences in the speeds, in f64; $-d$ for density exactly,
since $d \propto 1/\rho$ at fixed speeds) and the L2 misfit's 2 × 2 Hessian,
$\chi''_{ij} = w \langle \partial_i d, \partial_j d \rangle + w \langle d - d_{\text{obs}},
\partial_{ij} d \rangle$, against Pekeris's data for a solid 2% faster and 1% denser
($w = 1/\lVert d_{\text{obs}} \rVert^2$; the residual is 40.7% of the data). The solver's own
derivatives along the same directions come from five-point differences of its traces
($h = 10^{-2}$; for density the scheme scales exactly as the solid does), and the solver's
Hessian entries from those, with Pekeris's data and with its own for the faster, denser solid.

| | 10 m | 5 m | 5 m, `surface_slope` |
|---|---|---|---|
| Traces against Pekeris's | 15.5–28.2% | 6.4–9.5% | 4.3–8.2% |
| $\partial d / \partial(\text{speeds})$ against Pekeris's | 17.7–34.7% | 6.8–11.5% | 5.4–10.6% |
| Born against the solver's derivative | $0.8$–$1.6 \times 10^{-3}$ | $1.1$–$2.0 \times 10^{-3}$ | $1.1$–$2.1 \times 10^{-3}$ |
| Products against the solver's entries | $2$–$6 \times 10^{-4}$ | $3$–$6 \times 10^{-4}$ | $3$–$6 \times 10^{-4}$ |
| Solver's residual against Pekeris's data | 24.2% | 37.1% | 35.1% |
| $H$ (speeds², speeds–ρ, ρ²) against Pekeris's | +0.27, −33, +4.2% | +5.0, −3.9, +12.7% | −0.64, −11, −3.8% |
| its Gauss–Newton part | +0.17, +4.0, +7.4% | +7.8, +9.4, +9.2% | −2.9, −1.9, −0.52% |
| its second-order part | +2.1, −46, −8.1% | −70, −8.5, +26% | +53, −15, −16% |
| $H$ with the solver's own data | +0.84, +1.5, +6.3% | +8.0, +8.3, +9.1% | −2.7, −2.5, −0.86% |
| its second-order part | +15, +0.63, +2.5% | −2.7, +7.9, +8.6% | +7.8, −2.7, −2.2% |

Four things. (1) The Born run and the products are the solver's own derivatives to
$10^{-3}$ and a few $10^{-4}$; against the $10^{-5}$ of the physics test, the rest is the
layers, which a uniform change fills, and where the Born run's and the second adjoint run's
stress sources lack the CPML's memory terms (density's Born run is $1.1$–$1.5 \times
10^{-2}$ of its derivative, about the same absolute error on a derivative twenty times
smaller). An inversion keeps its unknowns outside the layers, where both are exact. (2) The
free surface's error carries into $J$ about as it is: the derivative is 1.1–1.3 times as far
from Pekeris's as the traces. (3) The Gauss–Newton part is off by about the traces' amplitude
error, squared: +7.8 to +9.4% at 5 m from the mirrored $v_z$'s kink (the Rayleigh wave's
swing 5.2–5.4% too large), within 0.5–2.9% with `surface_slope`; at 10 m it is close only
because the dispersion's loss of peak offsets the kink's larger excess. (4) The second-order
part correlates the residual with $\partial_{ij} d$, and against Pekeris's data the solver's
residual holds its own error, as large as a good part of the true residual: that part is then
8–70% off, and the whole entry up to 33% (10 m, speeds–density). With the solver's own data,
whose residual leaves its error out, the second-order part is within 0.6–15% and the entries
within about the Gauss–Newton part's error. So with real data the second-order terms are only
as good as the forward modelling, which they amplify; the Gauss–Newton part degrades
gracefully.

**Decision.** The illumination stays the default preconditioner; the diagonal is opt-in
(`hessian_probes`, default 0) and not worth its runs on this survey. The products are there
for what needs them: truncated-Newton steps, point-spread functions and resolution. Where the
forward error is a fair part of the residual (near a free surface at the default sampling),
prefer the Gauss–Newton product, or check the second-order part against its size first.

**Reproduce.** `just study gauss_newton_diagonal`; `cargo test --release --test hessian_study
-- --ignored --nocapture --test-threads 1`; `just physics
the_born_run_is_the_derivative_and_the_gradient_its_transpose`, `just physics
the_hessian_product_is_the_gradients_derivative`; `just unit hessian`; `just accuracy hessian_against_lambs_problem`; `cargo test --release --test
analytic_inversion gauss_newton_products_against_born -- --ignored --nocapture`.

## Open questions

- **A test that discriminates the misfits.** The slow start (2026-11-21) left the deep model
  untouched with every misfit in 8 iterations; a survey with more wavelengths to the targets,
  more iterations and the near-surface held (or the receivers' side of the illumination) is
  needed to show whether time shifts, phases or frequency continuation widen the basin here
  [@bunks1995; @luo1991; @fichtner2008]. Windows for single arrivals are now mutes
  (2026-11-21, "Mutes and windows"); choosing them per trace from the observed data is not.
- **Smoother time shifts.** The parabola's shift is continuous but its derivative jumps
  where the correlation's peak changes lag, which limits the finite-difference check to
  $10^{-4}$; Newton's method on the correlation's band-limited interpolant would be smooth and
  unbiased at the cost of every lag of every trace.
- **The phase misfit's scale** is the stage's corner frequency; a low-passed spectrum peaks
  at about 0.85 of it (a research agent's measurement), so $\sigma = 1/f_c$ may be 15% short.
  The receivers' side ignores a stage's filter too.
- **Inverting Marmousi2's elastic data.** From a start smoothed over 250 m, the slow waves
  along the seafloor carry most of the misfit (four fifths within 1 km of the source, more
  than half after a 400 m/s line) and the part of its gradient that turns $v_s$ away from the
  truth below 750 m (2026-11-21, "Where Marmousi2's misfit lies"). Frequency stages, time
  shifts, the phase misfit and the receivers' side don't get round them in five iterations
  (2026-11-21, the four Marmousi2 entries after it): the stages keep $v_s$ right down to 1 km, the phase misfit $v_p$, and
  nothing recovers either below 1 km. Muting them (2026-11-21, "Marmousi2: L2 with the slow
  seafloor waves muted" and the two entries after it) stops the damage: L2 then improves $v_p$
  down to 1 km and leaves $v_s$ alone, but the low-frequency slow waves that had helped $v_s$
  go too; kept at 0.6 Hz and muted above, they help $v_s$ above 1 km and hurt it below
  (2026-11-21, "Marmousi2: the slow waves kept at 0.6 Hz, muted above"). What the mute leaves
  of $v_s$'s gradient at 1–1.5 km points slightly away on both components, from arrivals before
  the line, not split further ("Marmousi2: what the mute removes"). On the held-out western
  part the ranking reverses (2026-11-23, "Marmousi2: the held-out western part"). Its start
  predicts the slow waves within half a cycle, so L2 alone recovers $v_s$ best and the mute
  makes $v_p$ worse near the seafloor. Whether to mute depends on cycle skipping, which
  nothing yet tests from the data. The whole section shows both within one survey: with the
  central part's mute, $v_s$ improves where the start is within half a cycle and drifts
  away at 10.5–14.5 km, where it isn't (2026-11-23, "Marmousi2: the whole section, muted
  L2"). Written but not run (`tests/marmousi2_study.rs`): muted
  L2 with the receivers' side
  (`muted_l2_with_receivers_side`), time shifts in a window around the early arrivals, then L2
  (`windowed_time_shifts_then_l2`) and muted L2 for twenty iterations (`muted_l2_longer`);
  L2 with the wavefields rebuilt (`l2_rebuilt`) has since matched L2 with checkpoints
  (2026-11-24). A test from the data of whether the slow
  waves are within half a cycle (2026-11-25, "Marmousi2: the slow-wave phase check") tells
  the parts apart, and the mutes it decides keep the west's L2 result, but in the centre
  recover only a third to a half of the plain mute's gain for $v_s$ below 750 m: what turns $v_s$ away
  there is not only cycle skipping. Split by the check's mutes (2026-11-26, "Marmousi2: a
  start smoothed through its slownesses, the check's split, and a hard seabed"), it is the
  slow waves beyond 300 m on the sides the check passes whole; a start smoothed through its
  slownesses (+11.4% under the central seafloor against +18.2%) moves the wrong-way update
  up rather than removing it; with Poisson's ratio capped at 0.25, as in the literature,
  $v_s$ moves the right way at every depth. Layer stripping by depth (`fixed_below`; 2026-11-26,
  "Marmousi2: layer stripping by depth") ends where L2 alone does in the centre, whose top
  250 m inverted alone barely moves, and costs the west $v_p$. Not written: the phase misfit with
  the mute; a
  first stage below 0.6 Hz (a 0.45 Hz low-pass keeps 0.8% of a 1.2 Hz Ricker's energy, a 0.6 Hz
  one 2.9%, enough for noise-free data); the slow waves' speed picked from the observed gathers,
  which stacks and f-k spectra within 1 km couldn't do.
- **Free surface for long stencils.** Replace stress imaging with a formulation of matching
  order [@mittet2002; @kristek2002], e.g. staggered summation-by-parts closures
  [@oreilly2017], which keep the exact discrete adjoint. Measured on Lamb's problem
  (2026-11-22, "The free surface against Lamb's problem"): the Rayleigh wave misfits by
  15–28% at the default's sampling, about 1% fast (second order) and 3–6% too large. The
  excess is the evenly mirrored $v_z$'s kink, read by sources and receivers on the surface
  (2026-11-23): about 11% per unit $kh_R$ per end, none with both at $h/2$. Continuing the
  mirrored fields with a one-sided slope (`surface_slope`, opt-in, 2026-11-23) leaves about
  2.5%; the boundary condition's slope, model-dependent, would be exact at every wavenumber.
  Whether either should be the default is open. The horizontal component (Pekeris's eq. 66, elliptic integrals)
  is not measured.
- **Interfaces.** At normal incidence, measured (2026-11-21, "A planar interface"):
  node-sampled layers put an interface midway between nodes (first-order error),
  band-limiting density and compliances [@mittet2017] makes it second order wherever it
  lies, and the stencil's order doesn't matter. Band-limited layers are the default
  (2026-11-21, "Band-limited layers, opt in"; 2026-11-22, "Band-limited layers by default").
  Still open: spheres, and how accurate the water rule is
  against a fluid–solid reference (Scholte waves); oblique incidence
  with converted waves, which needs point sources and an exact reference, or plane waves
  at an angle in a grid whose sides wrap around; and the S reflection's second-order error
  at the point densities FWI runs at.
- **Source and receiver positions.** Between nodes they now interpolate with the stencil's
  own weights (2026-11-22, "Sources and receivers between nodes"). Still open: sources and
  receivers at a fluid–solid interface or near a sharp contrast, where interpolating across
  it smooths what they see [@petersson2016].
- **Optimised coefficients, further.** Weight the error by the wavelets' spectrum and the
  model's range of speeds instead of $(\xi/\beta)^{-2}$; measure the observed order under grid
  refinement (the band shrinks with $h$, so the operator tends to Taylor's); and choose the
  band an inversion pins (the starting model's slowest wave, or the lowest speed allowed).
  Gradients with them agree with finite differences (2026-11-20, the adjoint and gradients
  with the new defaults).
- **Order 12 or 16 by default?** Not for the gradient's sake: from 4.6 points per shortest S
  wavelength every order from 8 up leaves the wavefields' part of the gradient's error at
  0.01–0.05%, and order 8 has 0.26% at 3.5 (2026-11-21, "The gradient against the exact
  kernels"). The case for them is cost on coarse grids (2026-11-20, optimised coefficients).
  (That entry also settled the time error in gradients: the transform takes it out.)
- **Absorbing layers at high order.** Tune the profile or width so reflection returns to
  $-76$ dB.
- **Error propagation into adjoint fields and gradients.** Forward-field error (dispersion,
  interfaces, free surface, layers), reconstruction error (time-reversed propagation in single
  precision; the non-reversible layers [@dussaud2008; @aaker2020]), and their effect on
  gradients and inversions [@tromp2005; @fichtner2006; @fichtner2011]. Settled: adjoint
  consistency (the dot-product entry) and the gradient of the discrete run against finite
  differences (the gradient entry), with an exactly replayed forward field, and the error of
  rebuilding it backwards from the band (2026-11-20: single-precision round-off,
  $10^{-7}$–$10^{-6}$ of the gradient, set by the field's loudest stretch in a segment; the
  layers and the free surface add nothing), against rebuilding it inside a one-node surface
  as Aaker et al. do [@aaker2020gji] (2026-11-27: $10^{-5}$–$10^{-3}$ of the gradient, the
  injection's shape, not round-off). Still open there: why a water layer drifts five
  times more (and with band-limited layers twice that again, 2026-11-22), and whether that
  drift moves an inversion; why the surface's gradient error was 10–100 times smaller on
  `examples/layered.toml` than on the small study model, and whether recording a second
  plane (a dipole term) would close most of its gap to the band for far fewer values. Dispersion in a homogeneous full space
  is settled too, against the exact kernels (2026-11-21, "The gradient against the exact
  kernels"): at the defaults the wavefields put 0.26% into the gradient at 3.5 points per
  shortest S wavelength and $\le 0.05\%$ from 4.6, less than into the traces; with the time
  error kept, or order 4, more. Still to measure: interfaces and the free surface, which have
  no closed-form kernels, against finer grids [@mittet2017]; and whether an inversion's model
  error follows the gradient's (data from a fine grid inverted on coarse ones).
