# Makefile for the APS (REVTeX) version
# Builds main.pdf and removes all auxiliary compilation files, keeping only the PDF.

MAIN   = main
LATEX  = pdflatex -interaction=nonstopmode -halt-on-error
BIBTEX = bibtex

# Auxiliary / intermediate files produced during compilation (everything but the PDF)
AUX = *.aux *.log *.bbl *.blg *.out *.toc *.lof *.lot *.spl *.end *.brf \
      *.fls *.fdb_latexmk *.synctex.gz *.nav *.snm *.vrb

all: $(MAIN).pdf

$(MAIN).pdf: $(MAIN).tex FRT.bib
	$(LATEX) $(MAIN)
	$(BIBTEX) $(MAIN)
	$(LATEX) $(MAIN)
	$(LATEX) $(MAIN)
	@$(MAKE) --no-print-directory clean

# Remove auxiliary files only (keep the compiled PDF)
clean:
	rm -f $(AUX)

# Remove everything that can be regenerated, including the PDF
distclean: clean
	rm -f $(MAIN).pdf

.PHONY: all clean distclean
