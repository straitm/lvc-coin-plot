#!/bin/bash

if [ $# -ne 1 ]; then
  echo Wrong number of arguments $#
  echo Syntax: $(basename $0) histfile
  echo Histfile has to have a standard trigger name in it
  exit 1
fi

histfile=$1

pdfbase=$(basename $histfile .hadded.root)

# If name has the form NNNx-streamname.hadded.root, it's a background
# sample, and NNN is the number of windows.  Otherwise, it's a signal file.
if echo $histfile | grep -q x; then
  nwindows=$(echo $histfile | cut -dx -f 1)
else
  nwindows=1
fi

if echo $histfile | grep fardet-t02; then
  trigname="FD 10Hz trigger"
  livetimediv=1
  longreadout=0
elif echo $histfile | grep neardet-ddactivity1; then
  trigname="ND energy"
  livetimediv=0
  longreadout=0
elif echo $histfile | grep fardet-ddenergy; then
  trigname="FD energy"
  livetimediv=0
  longreadout=0
elif echo $histfile | grep neardet-ddsnews ||
     echo $histfile | grep neardet-ligo; then
  trigname="ND long readout"
  livetimediv=1
  longreadout=1
elif echo $histfile | grep fardet-ddsnews ||
     echo $histfile | grep fardet-ligo; then
  trigname="FD long readout"
  livetimediv=1
  longreadout=1
else
  echo Could not figure out what kind of data $pdfbase is
  exit 1
fi

root -n -l -b -q $(dirname $0)/ligopass2.C+\
'("'"$histfile"'","'"$trigname"'","'"$pdfbase"'",'\
$livetimediv', '$longreadout', '$nwindows')' \
 | tee $pdfbase.log
