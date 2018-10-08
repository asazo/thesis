TEX = pdflatex -shell-escape -interaction=nonstopmode -file-line-error
ASY= asy
all: main.pdf

main.pdf : main.tex references.bib
	$(TEX) main.tex
	bibtex main.aux
	$(TEX) main.tex
	$(TEX) main.tex

clean:
	rm -rf *.log *.aux *.blg *.bbl *.synctex.gz *.dvi *.cfg *.idx *.toc *.lof *.lot *.out
