#! /bin/bash

lyx --export pdflatex estimation.lyx
xelatex estimation
bibtex estimation
xelatex estimation
