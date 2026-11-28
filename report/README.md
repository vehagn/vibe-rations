# The solver report

A LaTeX report on vibe-rations: its theory and numerics, how it is verified against analytic
solutions and published work, and how to use it. Its numbers, tables and figures come from
experiments it can rerun, each stamped with the commit and commands that made it, so the
report follows the solver as the solver changes.

```bash
just in report           # the PDF, from the committed data: report/build/vibe-rations-report.pdf
just in report-status    # every experiment: current, stale, dirty or missing
just in report-run       # rerun the quick experiments (names, --tier, --stale, --missing)
```

Building the PDF runs no simulation. It needs the dev container's TeX (pdflatex, BibTeX,
latexmk) and Python's standard library. In another Debian or Ubuntu machine:

```bash
apt-get install --no-install-recommends texlive-latex-base texlive-latex-recommended \
  texlive-latex-extra texlive-pictures texlive-science texlive-fonts-recommended lmodern latexmk
python3 report/report.py pdf
```

## Layout

```
report/
  main.tex               the document: \input of the preamble and of each section
  preamble.tex           packages, styles and the shared macros (colours.tex: the palette)
  sections/<section>.tex abstract, intro, usage, theory, numerics, verification, benchmarks,
                         reproducibility: one file each, one owner each
  figures/<section>/     hand-made figures: TikZ files and static pictures
  experiments/<section>/<name>.toml   an experiment: its title, tier and commands; what it
                         needs sits beside it under the same name: scripts (<name>.py) and
                         a directory of inputs (<name>/, say the solver's configurations)
  data/<section>/<name>/ COMMITTED: the experiment's results, reduced to small text files
                         (CSV, JSON, a little PNG), and provenance.json
  reportlib/             what the scripts share: .npy reader, CSV, values.json and PNG writers
  report.py              the generator: list, run, render, pdf, check
  build/                 git-ignored: raw run outputs (raw/), what LaTeX reads from the
                         generator (gen/), LaTeX's own files (latex/), the PDF
```

LaTeX runs in `report/`, so paths in the sources are relative to it: `\data/...` is
`report/data/...`, `\gen/...` is `report/build/gen/...`, `\repo/...` the repository's root.

## Experiments

Everything in the report that comes from running something is an experiment. Its manifest
names the commands; its data directory holds what they made.

```toml
# report/experiments/verification/stokes.toml
title = "A point force in a homogeneous full space against Stokes's solution"
tier = "quick"     # quick: up to two minutes; study: up to an hour; heavy: longer, or needs
                   # a large machine or downloaded data (say which under `needs`)
depends = ["tests/report_verification.rs", "tests/common"]   # besides the solver itself

[[command]]
run = "cargo test --release --test report_verification stokes -- --ignored --exact --nocapture"

[[command]]
run = "vibe-rations check examples/explosion.toml"
output = "check.txt"   # keep what it prints in {data}/check.txt, to show in the report
```

- Commands run in order from the repository's root, with the freshly built `vibe-rations` first
  on `PATH`. `{data}` and `{raw}` in a command stand for `report/data/<section>/<name>` and
  `report/build/raw/<section>/<name>`; the same paths are in the environment as
  `VIBE_RATIONS_REPORT_DATA` and `VIBE_RATIONS_REPORT_RAW`.
- **Programs write only where those say.** A test or script that finds the variables unset
  writes under `out/report/` instead (`reportlib.data_dir`), so the committed data change
  only through `report.py run`, which clears the data directory first, runs the commands,
  and then writes `provenance.json`: the commit, uncommitted changes in what the experiment
  depends on, the hash of each thing it depends on, the commands with their durations, the
  machine, the toolchain, and every file's size and SHA-256.
- **Committed data are small text.** Thin a curve to a few hundred points
  (`reportlib.thin`), write floats with six significant digits (`reportlib.write_csv`), keep
  rasters to small PNGs. `run` warns about a file over 512 KiB or a directory over 2 MiB.
  Bulky outputs stay in `{raw}`.
- **Status.** An experiment is *missing* without data; *stale* when anything it depends on
  (`src`, `Cargo.toml`, `Cargo.lock`, `rust-toolchain.toml`, its manifest, scripts and inputs,
  and its `depends`) differs from what its data were made from; *dirty* when it ran with uncommitted
  changes in those; otherwise *current*. The PDF marks every figure, table and appendix
  entry whose experiment isn't current, and `report.py check` fails if any isn't. So: commit
  the code, run the experiment, commit its data.
- Python scripts use the standard library only (no NumPy, no matplotlib), so the report
  builds wherever the dev container does. `reportlib.load_npy` reads the solver's arrays.

## Using the results in LaTeX

| In the source | Gives |
|---|---|
| `\addplot table[x=t, y=vz] {\data/verification/stokes/traces.csv};` | a curve, by pgfplots, from the CSV (comma-separated by default) |
| `\pgfplotstabletypeset[...]{\data/verification/stokes/misfits.csv}` | a table from a CSV, with booktabs rules |
| `\addplot graphics[xmin=..,xmax=..,ymin=..,ymax=..] {\data/.../field.png};` | a raster under real axes (`xz section` style: $z$ down, equal scales) |
| `\val{verification/stokes/misfit}` | a value from the experiment's `values.json`, as it stands |
| `\valnum[round-precision=2]{...}`, `\valqty{...}{\percent}` | the same through siunitx's `\num` and `\qty` |
| `\source{verification/stokes}` | "Data: verification/stokes, commit 0dd7470.", linked to the appendix; **end every generated figure's and table's caption with it** |
| `\commands{verification/stokes}` | the experiment's commands as a listing, where showing them teaches something (the appendix lists them all anyway) |

Figures are pgfplots reading the CSVs, so they are redrawn at every build and share the
text's fonts. The solver's curve is `solver` (blue), an exact or published answer
`reference` (black, dashed), their difference `difference`; further curves take the cycle
list's colours in order. Colour maps: `colormap name=vrdiverging` for signed fields,
`vrsequential` for magnitudes; `reportlib.write_png` draws with the same ones.

A missing value or experiment prints as a red **[? name]** rather than stopping the build.

## Conventions

- **Citations:** natbib, `\citet{virieux1986}` in running text and `\citep{...}` in
  parentheses, with the keys of `docs/references.bib`, which the report reads directly: it
  is the repository's one bibliography. Add an entry there, checked against the paper or
  its Crossref record as `AGENTS.md` asks, inside the block of its topic rather than at the
  end of the file (three threads appending at the end is three merge conflicts).
- **Labels:** `sec:<section>:<name>`, `eq:<section>:<name>`, `fig:<section>:<name>`,
  `tab:<section>:<name>`, with `<section>` the file's name; the section itself is
  `sec:<section>`. Refer with `\cref`.
- **Notation:** `\vect{v}`, `\tens{\sigma}`, `\dt`, `\dx` (the spacing $h$), `\vp`, `\vs`,
  `\T` (transpose), `\pd{f}{x}`; units through siunitx v3 (`\qty{10}{m}`, `\num{1e-7}`).
  A macro several sections need goes in `preamble.tex`; one a single section needs, at that
  section's top with `\providecommand`.
- **Code:** `\code{Solver::step_with}` inline (typed as is, underscores and all),
  `\file{src/solver.rs}` for a file (a link to it at the report's commit), `\vr` for the
  program's name; listings in the styles `shell`, `toml`, `rust` and `output`, and
  `\lstinputlisting[style=toml]{\repo/examples/explosion.toml}` to show a file as it is.
- **Boxes:** `keybox` for a point worth remembering, `trybox` for something to try,
  `\pending{...}` for text still to be written (it shows in the PDF).
- **Claims:** a number in the text is a `\val`, or stands next to the experiment it came
  from. "Bit for bit" is claimed only where a test or an experiment checked it.
