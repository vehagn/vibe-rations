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
