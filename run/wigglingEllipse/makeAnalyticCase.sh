#!/bin/sh

#Copying baseCase to generate analytic trajectory for comparison

cp -r baseCase analytic
cd analytic
foamDictionary system/controlDict -entry application -set kelvinKirchhoffEquationSolver
./Allclean
./Allrun
