#! /bin/bash

lyx --export pdflatex glm.lyx
xelatex glm
bibtex glm
xelatex glm