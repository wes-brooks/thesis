#! /bin/bash

lyx --export pdflatex estimation.lyx
xelatex estimation
biber estimation
xelatex estimation
