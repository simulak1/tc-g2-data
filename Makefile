# Makefile for tcjastrowfitting (flat layout)

# Use your conda Python here (adjust if needed)
PYTHON      := /home/simula/miniconda3/bin/python

MAIN_TEX    := main.tex
MAIN_PDF    := main.pdf

FIG_SCRIPT  := tools/plot_depth_sweep.py
FIG_CSV     := data/heat_tc_energies.csv
FIG_DIR     := figures

.PHONY: all paper figures clean

# Default: just build the paper
all: paper

# --------------------------------------------------------------------
# LaTeX manuscript
# --------------------------------------------------------------------
paper: $(MAIN_PDF)

$(MAIN_PDF): $(MAIN_TEX)
	@echo "==> Building manuscript PDF"
	pdflatex -interaction=nonstopmode $(MAIN_TEX)

# --------------------------------------------------------------------
# Figures from depth sweep script
# --------------------------------------------------------------------
figures:
	@echo "==> Generating figures with plot_depth_sweep.py"
	mkdir -p $(FIG_DIR)
	$(PYTHON) $(FIG_SCRIPT) $(FIG_CSV) --unit mHa --save

# --------------------------------------------------------------------
# Cleanup helper (does NOT touch figures/)
# --------------------------------------------------------------------
clean:
	@echo "==> Cleaning LaTeX aux files"
	rm -f main.aux main.log main.out main.toc main.fls main.fdb_latexmk \
	      figures.aux figures.log
