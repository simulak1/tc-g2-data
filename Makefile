.PHONY: all figures clean

all: figures main.pdf

figures:
	# python tools/script.py

main.pdf: main.tex
	latexmk -pdf main.tex

clean:
	latexmk -C
