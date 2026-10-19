#!/usr/bin/env bash
set -e
cd "$(dirname "$0")"

PDFLATEX="/Library/TeX/texbin/pdflatex"
BIBTEX="/Library/TeX/texbin/bibtex"

echo "=== Compiling paper/manuscript.pdf ==="
$PDFLATEX -interaction=nonstopmode manuscript.tex
$BIBTEX manuscript
$PDFLATEX -interaction=nonstopmode manuscript.tex
$PDFLATEX -interaction=nonstopmode manuscript.tex

echo "=== Compiling paper/supplementary.pdf ==="
$PDFLATEX -interaction=nonstopmode supplementary.tex
$BIBTEX supplementary || true
$PDFLATEX -interaction=nonstopmode supplementary.tex
$PDFLATEX -interaction=nonstopmode supplementary.tex

echo "=== LaTeX Compilation Completed Successfully! ==="
ls -lh manuscript.pdf supplementary.pdf
