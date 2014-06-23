#! /bin/bash

lyx --export pdflatex note.lyx
xelatex note
bibtex note
xelatex note
