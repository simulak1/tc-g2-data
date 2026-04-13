.PHONY: all figures tables clean

all: figures tables main.pdf

# --- environment check ---
.venv:
	@echo "No .venv found. Create one with:"
	@echo "  python3 -m venv .venv && .venv/bin/pip install -r requirements.txt"
	@exit 1

# --- data ---
data/ECP/atomization_comparison.csv: .venv \
		reference_data/umrigar_ae_reference.csv \
		reference_data/umrigar_molpro_reference.csv \
		data/ECP/ccsd_energies.csv
	bash ./scripts.sh/ECP/build_ae_data.sh
	
ALL_ELEC_CSV := $(shell find data/ALL_ELEC -type f -name '*.csv')

data/ALL_ELEC/atomization_energies_H05.csv: .venv $(ALL_ELEC_CSV)
	bash ./scripts.sh/ALL_ELEC/build_ae_data.sh

# --- figures ---
figures: data/ECP/atomization_comparison.csv data/ALL_ELEC/atomization_energies_H05.csv .venv
	bash ./scripts.sh/ECP/plot_ae_data.sh
	bash ./scripts.sh/ECP/plot_convergence.sh
	bash ./scripts.sh/ALL_ELEC/plot.sh

# --- tables ---
tables: data/ECP/atomization_comparison.csv data/ALL_ELEC/atomization_energies_H05.csv .venv
	bash ./scripts.sh/ECP/make_tables.sh
	bash ./scripts.sh/ECP/plot_timing.sh

# --- paper ---
main.pdf: main.tex refs.bib figures tables
	latexmk -pdf -bibtex main.tex

clean:
	latexmk -C
