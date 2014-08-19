#! /bin/bash
mkdir -p aux
mkdir -p aux-supplement

#lyx --export pdflatex estimation.lyx

#Compile the main document:
ruby use-template.rb --source=fromLyx.tex --destination=estimation.tex --template=template.tex
xelatex -output-directory=aux estimation
bibtex aux/estimation
xelatex -output-directory=aux estimation
xelatex -output-directory=aux estimation
mv aux/estimation.pdf ./

#Compile the supplement
ruby use-template.rb --source=rawSupplement.tex --destination=supplement.tex --template=template.tex
xelatex -output-directory=aux-supplement supplement
bibtex aux-supplement/supplement
xelatex -output-directory=aux-supplement supplement
xelatex -output-directory=aux-supplement supplement
mv aux-supplement/supplement.pdf ./