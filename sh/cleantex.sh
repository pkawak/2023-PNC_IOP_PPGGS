#!/bin/bash

name=`echo $1 | sed 's/\.tex//'`

rm -f ${name}.aux ${name}.log ${name}Notes.bib 
rm -f ${name}.bbl ${name}.blg ${name}-blx.bib
rm -f ${name}.run.xml
