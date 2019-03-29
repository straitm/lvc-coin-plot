#!/bin/bash

if [ $# -ne 5 ]; then
  echo Wrong number of arguments $#
  echo Syntax: $(basename $0) trigname histfile pdfbase '[divide|dontdivide] [longreadout|notlongreadout]'
  exit 1
fi

histfile=$1
trigname=$2
pdfbase=$3

if [ "$4" == divide ]; then
  livetimediv=1
elif [ "$4" == dontdivide ]; then
  livetimediv=0
else
  echo Give \"divide\" or \"dontdivide\" for third argument
  exit 1
fi

if [ "$5" == longreadout ]; then
  longreadout=1
elif [ "$5" == notlongreadout ]; then
  longreadout=0
else
  echo Give \"longreadout\" or \"notlongreadout\" for fourth argument
  exit 1
fi

root -n -l -b -q ligopass2.C+'("'"$histfile"'","'"$trigname"'","'"$pdfbase"'",'$livetimediv', '$longreadout')' \
 | tee $pdfbase.log
