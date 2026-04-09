.PHONY: all figures tables clean

all: figures tables main.pdf

# --- environment check ---
.venv:
	@echo "No .venv found. Create one with:"
	@echo "  python3 -m venv .venv && .venv/bin/pip install -r requirements.txt"
	@exit 1

# --- data ---
data/atomization_comparison.csv: .venv \
		reference_data/umrigar_ae_reference.csv \
		reference_data/umrigar_molpro_reference.csv \
		data/ccsd_energies.csv
	./scripts.sh/build_ae_data.sh

# --- figures ---
figures: data/atomization_comparison.csv .venv
	./scripts.sh/plot_ae_data.sh
	./scripts.sh/plot_convergence.sh

# --- tables ---
tables: data/atomization_comparison.csv .venv
	./scripts.sh/make_tables.sh
	./scripts.sh/plot_timing.sh

# --- paper ---
main.pdf: main.tex refs.bib figures tables
	latexmk -pdf -bibtex main.tex

clean:
	latexmk -C
