#! /bin/bash
mkdir -p aux
mkdir -p aux-supplement

#Compile the main document:
xelatex -output-directory=aux estimation
bibtex aux/estimation
xelatex -output-directory=aux estimation
xelatex -output-directory=aux estimation
mv aux/estimation.pdf ./

#Compile the supplement
xelatex -output-directory=aux-supplement supplement
bibtex aux-supplement/supplement
xelatex -output-directory=aux-supplement supplement
xelatex -output-directory=aux-supplement supplement
mv aux-supplement/supplement.pdf ./