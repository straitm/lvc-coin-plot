all: ligopass2-fardet-pulser.pdf \
     ligopass2-neardet-bnb.pdf \
     ligopass2-neardet-ddactivity1.pdf

ligopass2-fardet-pulser.pdf: ligopass2.C ligopass2.sh
	./ligopass2.sh ../ligotest/fardet_t02.hists.root ligopass2-fardet-pulser 1 0

ligopass2-neardet-bnb.pdf: ligopass2.C ligopass2.sh
	./ligopass2.sh ../ligotest/neardet_bnb.hists.root ligopass2-neardet-bnb 1 0

ligopass2-neardet-ddactivity1.pdf: ligopass2.C ligopass2.sh
	./ligopass2.sh ../ligotest/neardet_ddactivity1.hists.root ligopass2-neardet-ddactivity1 0 0

ndmev-unslicedbighits-10s.pdf: ligopass2.C ligopass2.sh
	./ligopass2.sh ../ligotest/1516535068/ndmev.ligohists.root ndmev 1 0

ddsnews-feb2-unslicedhitpairs-1s.pdf: ligopass2.C ligopass2.sh
	./ligopass2.sh ../ligotest/1549117743/mev.ligohists.root ddsnews-feb2 1 0


feb2: fardet-ddsnews-feb2.pdf neardet-ddsnews-feb2-tracks-1s.pdf \
      ndactivity1-feb2-tracks-1s.pdf ddenergy-feb2-energy_low_cut-10s.pdf \
      t02-feb2.pdf


fardet-ddsnews-feb2.pdf: ligopass2.C ligopass2.sh
	./ligopass2.sh ../ligotest/1549117741-fardet-ddsnews/minbiasfd.ligohists.root fardet-ddsnews-feb2 1 1

neardet-ddsnews-feb2-tracks-1s.pdf: ligopass2.C ligopass2.sh
	./ligopass2.sh ../ligotest/1549117741-neardet-ddsnews/minbiasnd.ligohists.root neardet-ddsnews-feb2 1 1

t02-feb2.pdf: ligopass2.C ligopass2.sh
	./ligopass2.sh ../ligotest/1549117741-fardet-t02/minbiasfd.ligohists.root t02-feb2 1 0

ndactivity1-feb2-tracks-10s.pdf: ligopass2.C ligopass2.sh
	./ligopass2.sh ../ligotest/1549117741-neardet-ddsnews/ndactivity.ligohists.root ndactivity1-feb2 0 0

ddenergy-feb2-energy_low_cut-10s.pdf: ligopass2.C ligopass2.sh
	./ligopass2.sh ../ligotest/1549117741-fardet-ddenergy/ddenergy.ligohists.root ddenergy-feb2 0 0
