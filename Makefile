all: ligopass2-fardet-pulser.pdf \
     ligopass2-neardet-bnb.pdf \
     ligopass2-neardet-ddactivity1.pdf

ligopass2-fardet-pulser.pdf: ligopass2.C ligopass2.sh
	./ligopass2.sh ../ligotest/fardet_t02.hists.root ligopass2-fardet-pulser 1

ligopass2-neardet-bnb.pdf: ligopass2.C ligopass2.sh
	./ligopass2.sh ../ligotest/neardet_bnb.hists.root ligopass2-neardet-bnb 1

ligopass2-neardet-ddactivity1.pdf: ligopass2.C ligopass2.sh
	./ligopass2.sh ../ligotest/neardet_ddactivity1.hists.root ligopass2-neardet-ddactivity1 0
