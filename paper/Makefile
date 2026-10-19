PDFLATEX = /Library/TeX/texbin/pdflatex
BIBTEX = /Library/TeX/texbin/bibtex

all: manuscript.pdf supplementary.pdf

manuscript.pdf: manuscript.tex references.bib
	$(PDFLATEX) -interaction=nonstopmode manuscript.tex
	$(BIBTEX) manuscript
	$(PDFLATEX) -interaction=nonstopmode manuscript.tex
	$(PDFLATEX) -interaction=nonstopmode manuscript.tex

supplementary.pdf: supplementary.tex references.bib
	$(PDFLATEX) -interaction=nonstopmode supplementary.tex
	$(BIBTEX) supplementary || true
	$(PDFLATEX) -interaction=nonstopmode supplementary.tex
	$(PDFLATEX) -interaction=nonstopmode supplementary.tex

clean:
	rm -f *.aux *.bbl *.blg *.log *.out *.toc *.lot *.lof *.synctex.gz *.pdf

.PHONY: all clean
