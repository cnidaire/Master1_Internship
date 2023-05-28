all: report.pdf

report.tex: report.md
	pandoc -f gfm -t latex+raw_tex $< --toc --toc-depth=2 --number-sections -s --template=eisvogel.tex --listings -o $@

report.pdf: report.tex
	pdflatex $<
