#!/bin/bash

if [ $# -ne 4 ]; then
  echo Three arguments: histfile, pdfbase, livetimediv, longreadout
  exit 1
fi

histfile=$1
pdfbase=$2
livetimediv=$3
longreadout=$4

root -n -l -b -q ligopass2.C+'("'$histfile'","'$pdfbase'",'$livetimediv', '$longreadout')'
