slides-seminar: slides-seminar.Rnw

	Rscript -e 'knitr::knit("slides-seminar.Rnw")'
	pdflatex slides-seminar.tex
	-bibtex slides-seminar.aux
	pdflatex slides-seminar.tex
	pdflatex slides-seminar.tex
	make -s clean

clean:
	rm -f *.aux *.bbl *.blg *.brf *.idx *.ilg *.ind *.lof *.log \
	.*lol *.lot *.out *.toc *.synctex.gz *.nav *.snm
