#!/bin/bash

# make a pdf file that shows the 
# changes since a specified commit or tag

#old="HEAD"
old="v1.0" # named tag

cd ..
latexdiff-vc --flatten --math-markup 1 --force -r ${old} manuscript/manuscript.tex
#latexdiff-vc -r ${old} manuscript/manuscript.tex
#(don't flatten v0.1, since it had no .bbl file)

cd manuscript/
#pdflatex manuscript-diff${old}.tex
#pdflatex manuscript-diff${old}.tex
#pdflatex manuscript-diff${old}.tex

# if error's use this
pdflatex -interaction scrollmode manuscript-diff${old}.tex
pdflatex -interaction scrollmode manuscript-diff${old}.tex
pdflatex -interaction scrollmode manuscript-diff${old}.tex

for ext in .aux .blg .log Notes.bib; do
  rm -f manuscript-diff${old}${ext}
done
