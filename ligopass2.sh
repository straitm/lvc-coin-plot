#!/bin/bash

if [ $# -ne 3 ]; then
  echo Three arguments: histfile, pdfbase, livetimediv
  exit 1
fi

histfile=$1
pdfbase=$2
livetimediv=$3

root -n -l -b -q ligopass2.C+'("'$histfile'","'$pdfbase'",'$livetimediv')'
