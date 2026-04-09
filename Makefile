.PHONY: all figures clean

all: figures main.pdf

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

# --- paper ---
main.pdf: main.tex figures
	latexmk -pdf main.tex

clean:
	latexmk -C
