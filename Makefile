all: report.pdf

report.tex: report.md
	pandoc -f markdown+raw_tex -t latex $< --toc --toc-depth=2 --number-sections -s --template=eisvogel.tex --citeproc --natbib --listings -o $@

report.pdf: report.tex
	pdflatex $<
	bibtex report.aux
	pdflatex $<
	pdflatex $<
