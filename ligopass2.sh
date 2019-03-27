#!/bin/bash

if [ $# -ne 4 ]; then
  echo Three arguments: histfile, pdfbase, livetimediv, longreadout
  exit 1
fi

histfile=$1
pdfbase=$2

if [ "$3" == divide ]; then
  livetimediv=1
elif [ "$3" == dontdivide ]; then
  livetimediv=0
else
  echo Give \"divide\" or \"dontdivide\" for third argument
  exit 1
fi

if [ "$4" == longreadout ]; then
  longreadout=1
elif [ "$4" == notlongreadout ]; then
  longreadout=0
else
  echo Give \"longreadout\" or \"notlongreadout\" for fourth argument
  exit 1
fi

root -n -l -b -q ligopass2.C+'("'$histfile'","'$pdfbase'",'$livetimediv', '$longreadout')'
