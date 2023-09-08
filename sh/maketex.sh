#!/bin/bash

name=`echo $1 | sed 's/\.tex//'`

if [ ! -f "refs.bib" ]; then
  echo "no bibliography"
  pdflatex $1
  pdflatex $1
else
#  if [ -f "${name}.aux" ]; then
#    pdflatex $1
#    bibtex ${name}
#    pdflatex $1
#  else
    pdflatex $1 
    bibtex ${name}
    pdflatex $1
    pdflatex $1
#  fi
fi

