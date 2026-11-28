# latexmk's settings for the report: pdflatex and BibTeX, everything written under build/latex.
# `report/report.py pdf` runs it after writing build/gen; `latexmk` here does the same by hand.
use Cwd;

$pdf_mode = 1;
$bibtex_use = 2;
$out_dir = 'build/latex';
$pdflatex = 'pdflatex -interaction=nonstopmode -halt-on-error -file-line-error %O %S';
@default_files = ('main.tex');
# The bibliography is the repository's, docs/references.bib: the generator's copy without the
# `note` fields (which record how each entry was checked) if it is there, else the file itself.
$ENV{'BIBINPUTS'} = join(':', Cwd::abs_path('build/gen') // '', Cwd::abs_path('../docs'),
  $ENV{'BIBINPUTS'} // '');
