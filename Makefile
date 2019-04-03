all: feb2.pdf

feb2.pdf: \
	  neardet-ddsnews-feb2.pdf \
          fardet-ddsnews-feb2.pdf \
	  fardet-t02-feb2.pdf \
	  neardet-ddactivity1-feb2.pdf \
	  fardet-ddenergy-feb2.pdf \
          
	pdftk $^ cat output $@

.SUFFIXES: .pdf

# Don't change the trigger names without looking downstream to see what that effects
fardet-ddsnews-feb2.pdf: ligopass2.C ligopass2.sh
	./ligopass2.sh ../ligotest/1549117741-fardet-ddsnews/minbiasfd.ligohists.root       "FD long readout"   $* divide longreadout

fardet-t02-feb2.pdf: ligopass2.C ligopass2.sh
	./ligopass2.sh ../ligotest/1549117741-fardet-t02/minbiasfd.ligohists.root           "FD 10Hz trigger"   $* divide notlongreadout

neardet-ddsnews-feb2.pdf: ligopass2.C ligopass2.sh
	./ligopass2.sh ../ligotest/1549117741-neardet-ddsnews/minbiasnd.ligohists.root      "ND long readout"   $* divide longreadout

neardet-ddactivity1-feb2.pdf: ligopass2.C ligopass2.sh
	./ligopass2.sh ../ligotest/1549117741-neardet-ddactivity1/ndactivity.ligohists.root "ND energy trigger" $* dontdivide notlongreadout

fardet-ddenergy-feb2.pdf: ligopass2.C ligopass2.sh
	./ligopass2.sh ../ligotest/1549117741-fardet-ddenergy/ddenergy.ligohists.root       "FD energy trigger" $* dontdivide notlongreadout

fardet-ddenergy-feb3.pdf: ligopass2.C ligopass2.sh
	./ligopass2.sh /pnfs/nova/scratch/users/mstrait/ligo/ddenergy.ligohists.root       "FD energy trigger" $* dontdivide notlongreadout
