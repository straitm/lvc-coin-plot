#!/bin/bash

histfile=$1

root -b -q ligopass2.C'("'$histfile'")'
