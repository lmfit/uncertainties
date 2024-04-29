# Makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
PAPER         =

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d _build/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .
PDFFILE = uncertainties.pdf

.PHONY:  all html help clean latex pdf

html:
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) _build/html
	@echo
	@echo "Build finished. The HTML pages are in _build/html."

all: html pdf
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) _build/html
	-cp _build/latex/$(PDFFILE) _build/html/.
	@echo
	-cd _build && ln -s html uncertainties_doc && zip -pur uncertainties_doc.zip uncertainties_doc/* && mv uncertainties_doc.zip html && rm -f uncertainties_doc
	@echo "Build finished. The HTML pages are in _build/html."

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  all       to make standalone HTML files with PDF and zipped HTML"
	@echo "  html      to make standalone HTML files"
	@echo "  pdf       to make PDF"
	@echo "  latex     to make LaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  clean     to clean folders"

clean:
	-rm -rf _build/*

latex:
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) _build/latex
	@echo
	@echo "Build finished; the LaTeX files are in _build/latex."
	@echo "Run \`make all-pdf' or \`make all-ps' in that directory to" \
	      "run these through (pdf)latex."

pdf: latex
	(cd _build/latex; $(MAKE) all-pdf)
