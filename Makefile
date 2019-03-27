all: feb2.pdf

feb2.pdf: fardet-ddsnews-feb2.pdf \
	  fardet-t02-feb2.pdf \
	  neardet-ddsnews-feb2.pdf \
	  neardet-ddactivity1-feb2.pdf \
	  neardet-ddenergy-feb2.pdf
	pdftk $^ cat output $@

fardet-ddsnews-feb2.pdf: ligopass2.C ligopass2.sh
	./ligopass2.sh ../ligotest/1549117741-fardet-ddsnews/minbiasfd.ligohists.root $* divide longreadout

fardet-t02-feb2.pdf: ligopass2.C ligopass2.sh
	./ligopass2.sh ../ligotest/1549117741-fardet-t02/minbiasfd.ligohists.root $* divide notlongreadout

neardet-ddsnews-feb2.pdf: ligopass2.C ligopass2.sh
	./ligopass2.sh ../ligotest/1549117741-neardet-ddsnews/minbiasnd.ligohists.root $* divide longreadout

ndactivity1-feb2.pdf: ligopass2.C ligopass2.sh
	./ligopass2.sh ../ligotest/1549117741-neardet-ddsnews/ndactivity.ligohists.root $* dontdivide notlongreadout

fardet-ddenergy-feb2.pdf: ligopass2.C ligopass2.sh
	./ligopass2.sh ../ligotest/1549117741-fardet-ddenergy/ddenergy.ligohists.root $* dontdivide notlongreadout
