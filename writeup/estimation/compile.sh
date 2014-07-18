#! /bin/bash
mkdir -p aux

lyx --export pdflatex estimation.lyx

ruby use-template.rb

xelatex -output-directory=aux estimation
bibtex aux/estimation
xelatex -output-directory=aux estimation
xelatex -output-directory=aux estimation

mv aux/estimation.pdf ./