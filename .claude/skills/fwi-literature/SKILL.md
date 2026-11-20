---
name: fwi-literature
description: Look up the theory behind vibe-rations and elastic full waveform inversion - adjoint equations, sensitivity kernels, misfit functions, optimisation, PML and free-surface methods - in Fichtner's book and the local paper library. Use before designing or reviewing FWI or numerics work.
---

# Finding the theory

## Fichtner (2011), *Full Seismic Waveform Modelling and Inversion*

The main reference for this project's FWI work. If the maintainer's paper library is on
this machine (`~/claude/papers/`), the book is
`full-seismic-waveform-modeling-and-inversion.pdf`. **PDF page = printed page + 19.**

| Topic | § | Printed page |
|---|---|---|
| Elastic wave equation, formulations | 2.2 | 11 |
| Staggered grid in 3D | 3.2.1 | 39 |
| Free surface: vacuum, image and interior methods | 3.2.3 | 45 |
| Point sources on a grid | 3.4 | 53 |
| PML, convolutional PML | 6.3 | 99 |
| Non-linear minimisation: steepest descent, CG, Newton, multiscale, regularisation | 7 | 113 |
| The time-domain continuous adjoint method | 8 | 141 |
| Adjoint equations for the elastic wave equation | 8.5 | 157 |
| Gradients for isotropic elastic parameters (ρ, λ, μ and others) | 9.1.1 | 165 |
| Hessian kernels | 9.3 | 173 |
| Misfits and adjoint sources: L2, cross-correlation time shifts, time–frequency phase | 11 | 193 |

Extract text with poppler, on the host:

```bash
pdftotext -f 176 -l 182 -layout ~/claude/papers/full-seismic-waveform-modeling-and-inversion.pdf - | less
```

(`-f`/`-l` are PDF pages: §8.5, printed 157–163, is PDF 176–182.) Equations come out
mangled; read the rendered page with your PDF-reading tool when the details matter.

## Papers in the library

Search file names with `ls ~/claude/papers | grep -i <word>`. The most useful:

- Velocity–stress FD: Virieux 1984 and 1986.
- FWI foundations: Tarantola 1984, 1986, 1988; Mora 1987; Gauthier 1986.
- Adjoints and kernels: Tromp, Tape & Liu 2005; Fichtner, Bunge & Igel 2006; Fichtner &
  Trampert 2011 (Hessians, resolution).
- 3D elastic FWI in practice: Raknes, Arntsen & Weibull 2015 (Sleipner); Butzer, Kurzmann &
  Bohlen 2013; Epanomeritakis et al. 2008 (Newton-CG).
- Parametrisation and trade-offs: Köhn et al. 2012; Operto et al. 2013.
- Optimisation: Métivier et al. 2013 (truncated Newton); Pratt, Shin & Hicks 1998.
- Multiscale: Fichtner et al. 2013.
- Misfits: Engquist & Froese 2014, Yang, Engquist & Sun 2016 (Wasserstein).

Cite what you use in doc comments and commit messages, as the code does (e.g. "Levander
1988", "Komatitsch & Martin 2007").
