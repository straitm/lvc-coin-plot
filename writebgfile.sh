#!/bin/bash

e=$(basename $PWD | cut -d- -f 2)

rm -f *x-neardet-ddactivity1.hadded.root *x-neardet-t02.hadded.root

nnd=$(ls *ddactivity1.hadded.root | wc -l)
nfd=$(ls         *t02.hadded.root | wc -l)

echo $nnd ND files and $nfd FD files

ndhadded=${nnd}x-neardet-ddactivity1.hadded.root
fdhadded=${nfd}x-fardet-t02.hadded.root

if ! hadd $ndhadded *ddactivity1.hadded.root; then exit 1; fi
if ! hadd $fdhadded *t02.hadded.root; then exit 1; fi

(../ligopass2/ligopass2.sh $ndhadded; ../ligopass2/ligopass2.sh $fdhadded) | \
  grep TRACKS | \
  awk '{if(NF == 2) print $0}' | \
  grep -v ^TRACKS_POINT_0_FDMINBIAS_P | \
  grep -v ^HALFCONTAINED_TRACKS_POINT_1_FDMINBIAS_P | \
  sed 's/ *$//' | \
  tee ../ligopass2/$e.bg
