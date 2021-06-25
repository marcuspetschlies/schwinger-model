#include Makefile.in

main = main
out  = $(main)

.PHONY: clean

all: pdf

pdf:
	pdflatex ${main}
	pdflatex ${main}
	bibtex ${main}
	bibtex ${main}
	pdflatex ${main}
	pdflatex ${main}

pdfnb:
	pdflatex ${main}
	pdflatex ${main}

clean:
	rm -f *~ *.log *.toc *.bbl *.blg *.lof *.lot *.aux *.out
