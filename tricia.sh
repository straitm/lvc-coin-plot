#!/bin/bash

if [ $# -ne 2 ]; then
  echo Wrong number of arguments $#
  echo Syntax: $(basename $0) bglistfile histfile
  echo Histfile has to have a standard trigger name in it
  exit 1
fi

bgfiles=$1

histfile=$2

outdir=/nova/ana/users/mstrait/ligotriciaplots

mkdir -p $outdir

pdfbase=$outdir/tricia-$(basename $histfile .hadded.root)

if echo $histfile | grep -q fardet-t02; then
  trigname="FD 10Hz trigger"
  livetimediv=0 # because it doesn't work well
elif echo $histfile | grep -q neardet-ddactivity1; then
  trigname="ND energy"
  livetimediv=0
elif echo $histfile | grep -q fardet-ddenergy; then
  trigname="FD energy"
  livetimediv=0
elif echo $histfile | grep -q neardet-ddsnews ||
     echo $histfile | grep -q neardet-ligo; then
  trigname="ND long readout"
  livetimediv=1
elif echo $histfile | grep -q fardet-ddsnews ||
     echo $histfile | grep -q fardet-ligo; then
  trigname="FD long readout"
  livetimediv=1
else
  echo Could not figure out what kind of data $pdfbase is
  exit 1
fi

root -n -l -b -q $(dirname $0)/tricia.C+\
'("'"$bgfiles"'",\
"'"$histfile"'",\
"'"$trigname"'",\
"'"$pdfbase"'",'\
$livetimediv')' \
2> /dev/stdout | tee $pdfbase.log
