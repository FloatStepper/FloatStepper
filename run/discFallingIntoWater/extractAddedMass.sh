#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

caseDir=${1}
cd $caseDir
stateFile=$(find postProcessing -name addedMassWriter.dat)
grep "^Time = " ${stateFile} | cut -d' ' -f3 > tim
grep -A3 "^Time = " ${stateFile} | grep -B1 "^--" | grep "^(" | cut -d' ' -f2 | cut -d' ' -f1 > A22