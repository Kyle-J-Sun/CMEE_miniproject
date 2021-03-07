#!/bin/bash

# to remove the extensions from the input file
name=${1%.*}
texcount -1 -sum $1 > $name.sum
pdflatex $1
pdflatex $1
bibtex $name
pdflatex $1
pdflatex $1
evince $name.pdf &

## Cleanup
rm *~
rm *.aux
rm *.dvi
rm *.log
rm *.nav
rm *.out
rm *.snm
rm *.toc
rm *.bbl
rm *.blg
rm *.fls
rm *.fdb_latexmk
