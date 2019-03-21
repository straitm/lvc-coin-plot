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


feb2: fardet-ddsnews-feb2-tracks-1s.pdf neardet-ddsnews-feb2-tracks-1s.pdf \
      ndactivity1-feb2-tracks-1s.pdf ddenergy-feb2-energy_low_cut-10s.pdf \
      t02-feb2-set.pdf

t02-feb2-set.pdf: t02-feb2-tracks_point_0-10s.pdf
	@echo making $<
	@pdftk \
	  t02-feb2-rawhits-10s.pdf \
	  t02-feb2-rawhits-1s.pdf \
	  t02-feb2-unslice4ddhits-10s.pdf \
	  t02-feb2-unslice4ddhits-1s.pdf \
	  t02-feb2-unslicedbighits-10s.pdf \
	  t02-feb2-unslicedbighits-1s.pdf \
	  t02-feb2-unslicedhitpairs-10s.pdf \
	  t02-feb2-unslicedhitpairs-1s.pdf \
          t02-feb2-contained_slices-10s.pdf \
	  t02-feb2-contained_slices-1s.pdf \
	  t02-feb2-tracks-10s.pdf \
	  t02-feb2-tracks-1s.pdf \
	  t02-feb2-tracks_point_1-10s.pdf \
	  t02-feb2-tracks_point_1-1s.pdf \
	  t02-feb2-tracks_point_0-10s.pdf \
	  t02-feb2-tracks_point_0-1s.pdf \
	  t02-feb2-halfcontained_tracks-10s.pdf \
	  t02-feb2-halfcontained_tracks-1s.pdf \
	  t02-feb2-halfcontained_tracks_point_0-10s.pdf \
	  t02-feb2-halfcontained_tracks_point_0-1s.pdf \
	  t02-feb2-halfcontained_tracks_point_1-10s.pdf \
	  t02-feb2-halfcontained_tracks_point_1-1s.pdf \
	  t02-feb2-fullycontained_tracks-10s.pdf \
	  t02-feb2-fullycontained_tracks-1s.pdf \
	  t02-feb2-fullycontained_tracks_point_0-10s.pdf \
	  t02-feb2-fullycontained_tracks_point_0-1s.pdf \
	  t02-feb2-fullycontained_tracks_point_1-10s.pdf \
	  t02-feb2-fullycontained_tracks_point_1-1s.pdf \
	  cat output t02-feb2-set.pdf


fardet-ddsnews-feb2-tracks-1s.pdf: ligopass2.C ligopass2.sh
	./ligopass2.sh ../ligotest/1549117741-fardet-ddsnews/minbiasfd.ligohists.root fardet-ddsnews-feb2 1 1

neardet-ddsnews-feb2-tracks-1s.pdf: ligopass2.C ligopass2.sh
	./ligopass2.sh ../ligotest/1549117741-neardet-ddsnews/minbiasnd.ligohists.root neardet-ddsnews-feb2 1 1

t02-feb2-tracks_point_0-10s.pdf: ligopass2.C ligopass2.sh
	./ligopass2.sh ../ligotest/1549117741-fardet-t02/minbiasfd.ligohists.root t02-feb2 1 0

ndactivity1-feb2-tracks-10s.pdf: ligopass2.C ligopass2.sh
	./ligopass2.sh ../ligotest/1549117741-neardet-ddsnews/ndactivity.ligohists.root ndactivity1-feb2 0 0

ddenergy-feb2-energy_low_cut-10s.pdf: ligopass2.C ligopass2.sh
	./ligopass2.sh ../ligotest/1549117741-fardet-ddenergy/ddenergy.ligohists.root ddenergy-feb2 0 0
