report.pdf: report.md
	pandoc -f gfm -t latex $< --toc --toc-depth=2 --number-sections -o $@
