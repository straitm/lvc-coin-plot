#!/bin/bash

(../ligopass2/ligopass2.sh 114x-neardet-ddactivity1.hadded.root; ../ligopass2/ligopass2.sh 109x-fardet-t02.hadded.root) | grep TRACKS | awk '{if(NF == 2) print $0}' | grep -v ^TRACKS_POINT_0_FDMINBIAS_P | grep -v ^HALFCONTAINED_TRACKS_POINT_1_FDMINBIAS_P | sed 's/ *$//'
