#!/bin/bash

histfile=$1
pdfbase=$2
livetimediv=$3

root -n -l -b -q ligopass2.C+'("'$histfile'","'$pdfbase'",'$livetimediv')'
