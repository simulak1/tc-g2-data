.PHONY: all figures clean

all: figures main.pdf

figures: .venv
	# .venv/bin/python tools/script.py

.venv:
	@echo "No .venv found. Create one with:"
	@echo "  python3 -m venv .venv && .venv/bin/pip install -r requirements.txt"
	@exit 1

main.pdf: main.tex
	latexmk -pdf main.tex

clean:
	latexmk -C
