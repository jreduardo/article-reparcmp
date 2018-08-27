slides-ibc2018: slides-ibc2018.Rnw

	Rscript -e 'knitr::knit("slides-ibc2018.Rnw")'
	pdflatex slides-ibc2018.tex
	-bibtex slides-ibc2018.aux
	pdflatex slides-ibc2018.tex
	pdflatex slides-ibc2018.tex
	make -s clean

clean:
	rm -f *.aux *.bbl *.blg *.brf *.idx *.ilg *.ind *.lof *.log \
	.*lol *.lot *.out *.toc *.synctex.gz *.nav *.snm
