# Marmousi2 in the FWI literature, and how ours compares

The published full waveform inversions of Marmousi2 (Martin, Wiley & Marfurt 2006, the
elastic successor of Versteeg's 1994 Marmousi) that this repository's elastic inversion of
it can be set against: what each study ran, what it reported, and how a fair comparison
could be made. Gathered 2026-11-23 from the local paper library, web searches and Crossref.
Each study says how deeply it was read: *in full* (the paper itself), *full text, read by an
agent* (an HTML copy read by a reading agent, whose notes are summarised here), or *abstract
only*. Every work cited is in `docs/references.bib`, its fields checked against Crossref.
The research log's 2026-11-23 entry records the decision this note led to.

## Where our setup sits

- **Every elastic Marmousi2 inversion read in full that describes its model raised the
  shear velocities first.** Köhn et al. (2012) recalculated $v_s$ so that Poisson's ratio
  is at most 0.25, replacing the soft seabed with a hard one, since the sediments' low S
  velocities would have needed too fine a grid; Jeong, Lee & Min (2012) replaced the low
  $v_s$ the same way and removed the water; Cao et al. (2022) changed Poisson's ratio to
  about 0.25. Sun & Demanet (2022) don't say, but run 20 Hz on a 20 m grid, where the
  original sediments' $v_s$ would get under one point per wavelength; Pan, Geng & Innanen
  (2018) don't say which part of which Marmousi they used. Marmousi2's Poisson's
  ratio is above 0.25 almost everywhere: capped at 0.25, $v_s$ in our central 6 km rises at
  95% of the solid nodes, by 131% on average at 475–750 m, 98% at 750–1000 m, 59% at
  1–1.5 km and 26% at 1.5–2.5 km, and its slowest value from 315 to 886 m/s (the whole
  section: 91% of nodes, slowest 303 to 599 m/s; measured on `convert.sh`'s files, see the
  log). So those studies invert a model whose $v_s$ is close to $v_p/\sqrt3$, without the
  slow seafloor waves that carry most of our misfit (the log's "Where Marmousi2's misfit
  lies"). None of the studies found says it inverts Marmousi2 with its soft sediments kept.
- **Our frequencies are an order of magnitude lower.** The studies invert up to 8–20 Hz
  from sources of 6–10 Hz, mostly on 20 m grids, starting from 0.2 Hz (Jeong et al.), 0–5
  Hz (Cao et al.), up to 2 Hz (Köhn et al.), 2–4 Hz (Sun & Demanet) or 3–5 Hz (Pan, Geng &
  Innanen 2018). Ours is a 1.2 Hz Ricker (half its peak amplitude at 0.6 and 2 Hz): the
  soft sediments' 303 m/s already need 25 m nodes at 3 Hz. A P wavelength at our peak is
  1.3–2.5 km ($v_p$ of 1.6–3 km/s), several times our depth bands of 250–500 m, while the
  published images resolve structures as small as their 20 m grid (Köhn et al.).
- **Our start is closer to the truth than most.** Ours is the truth smoothed by a Gaussian
  of $\sigma$ = 250 m below the seafloor. Köhn et al. used a Gaussian of 800 m correlation
  length, Jeong et al. a velocity increasing linearly with depth, Sun & Demanet a "crude"
  start, Cao et al. a smoothed truth of unstated length.
- **We stop after 5 iterations; they ran 40–1300.** Pan et al. 40, Sun & Demanet 110, Cao et
  al. 180, Köhn et al. 285, Jeong et al. 800 and 500 more for density.
- **Almost nobody reports a number for the model.** Most judge by pictures and depth
  profiles. Pan, Geng & Innanen (2018) report the relative least-squares error
  $E_p = \lVert m_p^\text{est} - m_p^\text{true} \rVert / \lVert m_p^0 - m_p^\text{true}
  \rVert$, which is exactly our error ratio: 0.83, 0.72 and 1.04 for $v_p$, $v_s$ and
  $\rho$ after 40 iterations. Cao et al. report relative data residuals of 2.25–2.86%.
  Nobody reports the cosine with the true update, or any measure by depth.
- **We hold the true density; most of them invert it.** Köhn, Jeong, Pan and Sun & Demanet
  invert $\rho$ and find it the hardest parameter, carrying the velocities' structure
  (Köhn et al.; Pan et al. trace it to $v_s$, Sun & Demanet to $v_p$ at short offsets);
  Wang et al. (2016) invert it too. Cao et al. fix a smoothed density, the same in the
  data and the inversion. Holding the true $\rho$ makes our problem easier there, which a
  comparison should say.
- **For the slow seafloor waves, our mute is field practice.** On ocean-bottom field data
  Warner et al. (2013) muted the Scholte waves and Prieux et al. (2013b) left them
  unmodelled; Yuan & Simons (2014) tapered surface waves out of their elastic Marmousi
  data. Studies that use such waves fit envelopes first (Yuan, Simons & Bozdağ 2015) or
  strip layers by offset and frequency (Masoni et al. 2016). On a soft-seabed synthetic
  (Valhall), Brossier, Operto & Virieux (2009b, 2010) needed low frequencies, a hierarchy of
  parameters and components, and, with noisy data, an $\ell_1$ misfit to get $v_s$.

## Our setup

| | Central 6 km (`examples/marmousi2/central.toml`) | Whole section (`whole.toml`) |
|---|---|---|
| Model | Marmousi2 averaged onto 25 m nodes (Backus's harmonic moduli), $x$ = 6–12 km, extruded 1 km along $y$; 450 m of water, soft sediments kept ($v_s$ from 315 m/s, 303 in the whole section; $v_p/v_s \approx 4$–5) | all 17 km, extruded 5 km along $y$ |
| Grid | $241 \times 41 \times 141$, order 8 (optimised), CPML 16 nodes, free surface | $681 \times 201 \times 141$ |
| Acquisition | 5 explosions 1 km apart at 50 m depth; 201 ocean-bottom nodes 25 m apart on the first solid row (475 m), $v_x$ and $v_z$; offsets to 4.5 km; 3.5 s | 9 explosions 2 km apart; 641 nodes over 16 km |
| Source | Ricker, 1.2 Hz | the same |
| Start | the truth's $v_p$, $v_s$ smoothed by a Gaussian of $\sigma$ = 250 m over the solid nodes; water and true $\rho$ kept | the same |
| Inverted | $v_p$, $v_s$ as relative changes $m/m_0 - 1$; the top 500 m fixed | the same |
| Misfit | L2, balanced by shot and component; variants with low-pass stages (0.6, 1 Hz), time shifts, phase | L2 |
| Optimiser | L-BFGS, 5 iterations, preconditioned by the source-side illumination | the same |
| Mutes | the slow seafloor waves, before a 600 m/s line from 0.87 s | the same |
| Measures | error ratio and cosine with the true update by depth band (500–750, 750–1000, 1000–1500, 1500–2500, 2500–3125 m), and over all free nodes; the L2 misfit of all the data | (running on gauss) |
| Memory | 4.41 GiB (2.44 GiB rebuilt) | 48 GiB (22 GiB rebuilt) |

The best results so far (the log's 2026-11-21 entries), error ratio for $v_p$ / $v_s$ after
five iterations: muted L2 0.916 / 0.970 at 500–750 m, 0.971 / 1.004 at 750–1000 m, 0.999 /
1.001 at 1–1.5 km; over all free nodes 0.999 / 1.000. No variant recovers $v_s$ below 1 km.
The western 6 km, held out, has not been run yet.

## Elastic FWI of Marmousi2

| Study | Read | Model as used | Acquisition | Source, band (lowest) | Start | Inverted | Misfit, optimiser, iterations |
|---|---|---|---|---|---|---|---|
| Choi, Min & Shin (2008a) | abstract | elastic Marmousi-2 | not in the abstract | frequency domain | not in the abstract | elastic parameters and the wavelet | L2; conjugate gradient scaled by a new pseudo-Hessian |
| Choi, Min & Shin (2008b) | abstract | modified Marmousi-2 (as Köhn et al. cite it), water and solid coupled | streamer: hydrophone, 2C or 1C | frequency domain, also low-cut data | not in the abstract | velocities | L2; conjugate gradient |
| Köhn et al. (2012) | in full | 10 × 3.48 km cut of the complex centre; 500 m of water; Poisson's ratio ≤ 0.25 (hard seabed); 20 m grid | fixed streamer 40 m deep: 400 two-component geophones ($v_x$, $v_z$) 20 m apart, 100 shots at 40 m; 6 s | Ricker 10 Hz; stages up to 2, 5, 10, 20 Hz | the truth smoothed by a 2D Gaussian of 800 m correlation length | $v_p$, $v_s$, $\rho$ (and $\lambda$, $\mu$, $\rho$); $\rho$ bounded to 1000–3000 kg/m³ | L2 (time domain); conjugate gradient with a line search; 285 iterations; gradient zeroed above 380 m, tapered to 480 m, scaled with depth below |
| Jeong, Lee & Min (2012) | full text, read by an agent | 9.2 × 3.04 km, water and the sides removed; low $v_s$ replaced so Poisson's ratio is 0.25; 20 m | receivers at every surface node, 219 shots 40 m apart; 5 s | first derivative of a Gaussian, to 10 Hz; overlapping bands 0.2–2 … 0.2–10 Hz | linear increase with depth | stage 1 $v_p$, $v_s$ with $\rho$ fixed; stage 2 $\lambda$, $\mu$, $\rho$ | L2 (frequency domain); Fletcher–Reeves conjugate gradient scaled by the pseudo-Hessian; 800 + 500 iterations |
| Wang, Dong, Liu & Yang (2016) | abstract | elastic Marmousi2 | not in the abstract | frequency domain | not in the abstract | $v_p$, $v_s$, $\rho$ | nonlinear conjugate gradient preconditioned by a block-diagonal pseudo-Hessian |
| Pan, Liang & Habashy (2016) | abstract | a truncated Marmousi extended into 3D (version not stated) | not in the abstract | frequency domain | several | elastic velocities | Gauss–Newton and preconditioned nonlinear conjugate gradient |
| Pan, Geng & Innanen (2018) | in full | a 3.4 × 1.2 km Marmousi section ($v_p$, $v_s$, $\rho$; which version and part not stated) | 33 sources 100 m apart, 330 receivers 10 m apart, along the top surface | Ricker 6 Hz; bands 3–5 then 3–8 Hz | not described | $v_p$, $v_s$, $\rho$ | L-BFGS, 20 iterations per band; Hessian diagonals as preconditioners |
| Cao et al. (2022) | full text, read by an agent | 3D: Marmousi-II extended along $y$ at 45°, an uneven seabed; Poisson's ratio about 0.25; spectral elements of 60 m | streamer (38 shots, 9 streamers, offsets to 5.8 km) or two ocean-bottom cables of 38 four-component receivers (pressure, 3C) | Ricker 10 Hz; bands 0–5, 0–10, 0–20 Hz | smoothed truth (length not stated) | $v_p$, $v_s$; $\rho$ fixed to a smoothed model | L2 scaled per shot and component; L-BFGS; 60 iterations per band |
| Sun & Demanet (2022) | full text, read by an agent | 500 × 174 nodes at 20 m (about 10 × 3.5 km, the size of Köhn et al.'s grid); 440 m of water; changes to $v_s$ not stated | 50 shots 40 m deep from 0.8 to 8.64 km; 400 receivers ($v_x$, $v_z$) at 460 m, just below the seafloor; free surface; 6 s | Ricker 10 Hz; 2–4 Hz (true or extrapolated by a network), then 4–6, 4–10, 4–20 Hz | "crude", not described | $v_p$, $v_s$, $\rho$ | L2 (time domain, DENISE); L-BFGS; 30 + 30 + 30 + 20 iterations |

**What each found.**

- *Choi, Min & Shin (2008a)* scaled the gradient with a pseudo-Hessian that includes the
  forward wavefields' amplitudes, at no extra cost, and got closer to the true elastic
  Marmousi-2 than with the plain pseudo-Hessian. *Choi, Min & Shin (2008b)* modelled water
  and rock coupled and found the inversion beats an acoustic one, with single- or
  multicomponent and low-cut data; Köhn et al. cite it for the modified Marmousi-II giving
  much the same result from hydrophone, two-component or one-component data.
- *Köhn et al. (2012)*, whose code is DENISE, the closest to our numerics
  (velocity–stress, eighth order, time domain, particle-velocity receivers): with $v_p$,
  $v_s$, $\rho$ the fine layers absent from the start are resolved, the reservoirs found,
  and the residuals very small after 285 iterations. $v_s$ came from streamer data because
  the hard seabed converts P to S; density was recovered but carries $v_p$ and $v_s$
  structure; Lamé parameters gave a worse density with artefacts at layer boundaries.
  Low-frequency stages matter more for density than for the velocities. No numbers.
- *Jeong, Lee & Min (2012)*: inverted together, $v_p$, $v_s$ and $\rho$ give good
  velocities and a badly distorted density; velocities first with $\rho$ fixed, then
  $\lambda$, $\mu$, $\rho$ from them, recovers density with discrepancies at depth. The
  SEG/EAGE salt model, a second test, did less well under the salt. No numbers.
- *Wang et al. (2016)*: the block-diagonal preconditioner reduces cross-talk and speeds
  convergence, more at higher frequencies. A search result mentions a table of starting
  and inverted $v_p$, $v_s$, $\rho$ errors; not seen.
- *Pan, Liang & Habashy (2016)*: a 3D frequency-domain elastic inversion of Marmousi
  extended along $y$, as ours is, which gave good velocity images from several starts.
- *Pan, Geng & Innanen (2018)*: $v_s$ perturbations leak strongly into $\rho$ and, with
  reversed phase, into $v_p$, but receive little leakage themselves. Simultaneous inversion
  reached $E_{v_p}$ = 0.83, $E_{v_s}$ = 0.72, $E_\rho$ = 1.04 after 40 iterations;
  subtracting estimated contamination, 0.76, 0.67 and 0.83 at 2.5 times the cost. They note
  that the weak cross-talk in Köhn et al.'s results came after hundreds of iterations,
  which field data rarely allow.
- *Cao et al. (2022)*: ocean-bottom 3C data recover nearly all of $v_s$'s structure, where
  streamer or ocean-bottom pressure data, which see S only through P-wave amplitudes, merge
  two gas traps into one; $v_p$ from ocean-bottom data loses the streamer's footprint. A
  Scholte wave along the seabed is marked in the pressure data; no mute is described.
  Relative data residuals after 180 iterations: 2.25% (streamer), 2.70% (ocean-bottom
  pressure), 2.86% (3C). Each band took about 4 h on 1520–4560 cores.
- *Sun & Demanet (2022)*: from a crude start, 4–20 Hz data alone gave large errors; 2–4 Hz
  extrapolated from 4–20 Hz by a network gave about the same model as true 2–4 Hz data,
  and the 4–20 Hz stages from either came out close to each other. Density failed at 2–4 Hz
  and lacks long wavelengths at 4–20 Hz. The free surface leaves almost no energy below
  2 Hz. Figures only.

Also found, not read beyond the record: Biondi & Barnier (2020) estimate the elastic
parameters of a gas reservoir in Marmousi2 from pressure data, redatumed to the target;
Wehner et al. (2015) combine elastic waveforms and gravity for the density of Marmousi-II.

## Acoustic FWI of Marmousi with numbers

Acoustic studies of Marmousi or Marmousi2 are many; few report numbers for the model.

- *Métivier et al. (2013)*, in full: the Marmousi II $v_p$ model, 16 × 3.5 km at 25 m, its
  400 m of water held fixed, no free surface; 144 sources and 660 receivers 25 m deep;
  3, 5, 8 and 12 Hz inverted together from the truth smoothed by Seismic Unix's `smooth2`
  ($r_1 = r_2 = 20$); truncated Newton against l-BFGS, compared by convergence. They saw no
  difference between inverting the frequencies together or one after another.
- *Modrak & Tromp (2015)*, in full: Marmousi (version not stated) at 3 and 6 Hz among six
  test problems; L-BFGS needs 0.77 and 0.87 of nonlinear conjugate gradient's wavefield
  simulations with a bracketing line search, and 0.45 and 0.53 with backtracking.
- *Modrak & Tromp (2016)*, full text read by an agent: seven test cases, all acoustic,
  Marmousi and overthrust each on- and offshore (32 shots, 500 hydrophones 10 m deep, 500 m
  of water), starts made by smoothing the truth with Gaussians, L2 without mutes. L-BFGS
  saves 30–50% over nonlinear conjugate gradient; Marmousi converges slowest of the
  exploration cases, its diving waves weak and its reflections strong; and "one or two test
  cases are not enough" to choose settings. Model error is plotted against simulations,
  not tabulated.

## Soft seafloors and slow interface waves

- *Sears, Singh & Barton (2008)*, abstract: elastic inversion of wide-angle multicomponent
  ocean-bottom data; for low-velocity models with a soft seabed, P waves invert for $v_p$
  and both P and converted PS waves for $v_s$. *Sears, Barton & Singh
  (2010)*, abstract, did it on Alba: $v_p$ from $v_z$ (wide angles, then near offsets),
  then $v_s$ from the P waves' amplitude variation, then from PS waves on $v_x$, shallow
  first, below 15 Hz.
- *Brossier, Operto & Virieux (2009b)*, abstract: a Valhall synthetic with a soft seabed,
  where little P converts to S: acoustic inversion of the elastic data gives a good $v_p$;
  elastic inversion needs low frequencies and a hierarchy of parameter classes and data
  components to reach acceptable models. *Brossier, Operto & Virieux (2010)*, abstract:
  with noise, $\ell_1$ (and the Huber and hybrid criteria) recovers $v_s$ in that
  soft-seabed setting where $\ell_2$ does not.
- *Warner et al. (2013)* muted Scholte waves at short offsets in field ocean-bottom data;
  *Prieux et al. (2013b)* saw Scholte waves slower than 350 m/s at Valhall and left them
  out of the modelling (both read earlier; the log's "Mutes and windows").
- *Kugler et al. (2007)*, record only: Scholte waves themselves as data, for the shear
  velocity of shallow-water sediments, a possible shallow start.
- *Yuan & Simons (2014)*, abstract: elastic two-component inversion of the 2D Marmousi
  model (which version is not said), surface waves removed by tapering and filtering, yet
  matched by the final model. *Yuan, Simons & Bozdağ (2014)*, abstract, take an elastic
  Marmousi with its surface waves, and *Yuan, Simons & Bozdağ (2015)*, abstract,
  invert surface waves first with an envelope misfit, then their waveforms, then
  everything, on an industry synthetic model.
- *Masoni et al. (2016)*, abstract: layer stripping for surface waves, low-passing from
  high to low corners with offsets widening, so that each step updates deeper $v_s$.

## Benchmarks for held-out checks

The maintainer's caution against tuning to one dataset is what Modrak & Tromp (2016) argue
for optimisers. Models others have used to check settings:

- *Marmousi2's other parts*: our western 6 km, held out (not run), and the whole section.
- *The elastic overthrust model*, onshore (Aminzadeh et al. 1996 for the SEG/EAGE model;
  Brossier, Operto & Virieux 2009, 2010 elastic): no water, so no seafloor waves, but strong
  surface waves.
- *Valhall synthetic*, offshore with a soft seabed (Brossier et al. 2009b, 2010): the
  closest in physics to our problem; whether the model is available was not checked.
- *GO_3D_OBS* (Górszczyk & Operto 2021), public: a 3D crustal-scale subduction model for
  ocean-bottom seismometer surveys, used elastically by Cao et al. (2022).
- *The SEG/EAGE salt model* (Jeong et al.'s second case) and *BP 2004* (Billette &
  Brandsberg-Dahl 2005): salt, $v_p$-dominated; for $v_p$ settings only.

## The four closest, and how to compare fairly

Memory from `vibe-rations memory --invert` (the macOS build, on configurations in the
scratchpad; a stand-in layered model with $\Delta t$ and the layers' damping pinned, which
fixes the plan, for the 20 m grids). Times scale from `muted_l2`'s 2,769 s for five
iterations of five shots on four cores, starting gradient included: about 110 s per shot
and iteration on the central grid, and in proportion to nodes times steps elsewhere.

1. **Köhn et al. (2012): a hard-seabed control.** Their success, and Jeong's and Cao's, is
   on a model with Poisson's ratio capped at 0.25. The fair first comparison is our own
   survey on that model: `central.toml` and `central-true.toml` unchanged except their
   $v_s$, capped at $v_p/\sqrt3$ on the solid samples before averaging (the start capped
   then smoothed, as Köhn et al. smoothed the modified truth). It needs model files from a
   cap that `extrude` doesn't have (a `--max-poisson` option, or a script over the
   SEG-Y); no solver change. Same grid, steps and pins as `central.toml`, so the same
   4.41 GiB (2.44 GiB rebuilt) and about 45 minutes for five L2 iterations; the western
   part the same again. Report our error ratios and cosines by band, with and without the
   mute. If L2 then moves $v_s$ towards the truth below 750 m, our solver and inversion
   match the literature and the soft seabed is what separates us from it; if not, the
   band or the iterations are. Since the capped $v_s$ is nearly $v_p/\sqrt3$, gains in
   $v_s$ there partly restate $v_p$'s. A streamer variant (receivers 40 m deep in the
   water, Köhn's acquisition on our grid) costs the same. Their own survey at 20 m (100
   shots, 400 receivers, 6 s, 10.7 × 3.48 km with the layers, $535 \times 41 \times 175$
   nodes, 3,750 steps, 5.9 times the central part's node-steps) plans 19 GiB (11 GiB
   rebuilt) and about 11 minutes per shot and iteration on four cores, 18 hours per
   iteration for all shots: a few iterations of every tenth shot at most, on gauss, or on
   the Mac natively with the wavefield rebuilt.
2. **Sun & Demanet (2022): their frequency schedule on our acquisition.** The study whose
   acquisition is most like ours: $v_x$, $v_z$ just below the seafloor, sources 40 m deep,
   a free surface, the same kind of solver and optimiser. Their stages from 2 Hz need the
   hard seabed: with it, 25 m holds about 7.9 Hz at 4.5 points per wavelength (886 m/s),
   where the soft one holds about 3 Hz at the 4.2 that `central.toml` pins. So: the
   control above with a 3 Hz Ricker and stages below 2, 4 and 6 Hz, then all, five
   iterations each: 4.41 GiB, about 3 hours, with
   `optimised_for` re-pinned to the capped truth's points per wavelength at 7.5 Hz (about
   4.7). They report figures only; our measures would put numbers on that schedule.
   Their own survey at 20 m (50 shots, 400 receivers) plans 18 GiB (10 GiB rebuilt).
3. **Pan, Geng & Innanen (2018): the same error measure.** Their $E_p$ is our error ratio.
   Report it over all free nodes, not only by band, after 5, 10, 20 and 40 iterations: the
   written `muted_l2_longer` (20 iterations, about 3 hours, 4.41 GiB) gives the first
   three. Their 0.83 and 0.72 for $v_p$ and $v_s$ after 40 iterations at 3–8 Hz are the
   scale of a survey along the top of a small Marmousi section; their survey is not worth
   rebuilding, since the part of Marmousi it used isn't stated.
4. **Cao et al. (2022): the 3D ocean-bottom analogue.** Theirs is 3D, ocean-bottom,
   $v_p$ and $v_s$ with $\rho$ fixed, l-BFGS in bands: our `whole.toml`, now running on
   gauss, is the nearest of ours. On its results compute their measure, the relative
   residual $\lVert d - d_\text{obs} \rVert / \lVert d_\text{obs} \rVert$ per component
   (theirs: 2.25–2.86% after 180 iterations), beside our error ratios; and compare $v_s$
   logs at two wells through the reservoirs, as they did. Like for like would need the
   capped model too: the same 48 GiB (22 GiB rebuilt), on gauss only.

The soft-seabed studies (Brossier et al. 2009b, 2010 on Valhall; Sears et al. 2008, 2010)
are the closest in physics rather than setup: a soft seabed under ocean-bottom receivers.
Their remedies are the next things to try on our soft seabed, once the control has said how
much of the gap it accounts for: $v_p$ first with $v_s$ held, then $v_s$; the vertical
component before the horizontal; $\ell_1$ for noisy data.

**First:** the hard-seabed control (1), on the central part and the western one. It is the
cheapest (45 minutes and 4.41 GiB each), it makes the one change to the model that every
elastic Marmousi2 inversion read in full that describes its model made, on our own survey
and grid, and its outcome
says which of the literature applies to us: the hard-seabed studies above, or the
soft-seabed ones.

## Not read

Full texts that would sharpen this, behind publishers' pages (the SEG library refused the
reading agent): Choi, Min & Shin (2008a, b), for their Marmousi-2 acquisition and start;
Wang et al. (2016), for their table of errors; Yuan, Simons & Bozdağ (2015), to confirm
which model; Pan, Liang & Habashy (2016), for the 3D extension; and the soft-seabed papers
(Brossier et al. 2009b, 2010; Sears et al. 2008, 2010), for how soft their seabed is and
what they muted. Downloading any of them waits for the maintainer's word.
