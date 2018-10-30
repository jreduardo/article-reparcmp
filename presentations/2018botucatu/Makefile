slides-mca2018: slides-mca2018.Rnw

	Rscript -e 'knitr::knit("slides-mca2018.Rnw")'
	pdflatex slides-mca2018.tex
	-bibtex slides-mca2018.aux
	pdflatex slides-mca2018.tex
	pdflatex slides-mca2018.tex
	make -s clean

clean:
	rm -f *.aux *.bbl *.blg *.brf *.idx *.ilg *.ind *.lof *.log \
	.*lol *.lot *.out *.toc *.synctex.gz *.nav *.snm
