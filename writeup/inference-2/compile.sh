#! /bin/sh

mkdir -p aux

pdflatex -output-directory=aux entry 
biblatex -output-directory=aux entry 
pdflatex -output-directory=aux entry 
pdflatex -output-directory=aux entry 

mv aux/entry.pdf ./


