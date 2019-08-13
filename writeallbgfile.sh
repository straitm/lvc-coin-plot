#!/bin/bash

cd /nova/ana/users/mstrait

for d in ligobgresults-*; do (cd $d; ../ligopass2/writebgfile.sh); done
