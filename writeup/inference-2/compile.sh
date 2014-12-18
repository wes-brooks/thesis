#! /bin/sh

mkdir -p aux

pdflatex -output-directory=aux inference 
bibtex aux/inference
pdflatex -output-directory=aux inference
pdflatex -output-directory=aux inference

mv aux/inference.pdf ./

