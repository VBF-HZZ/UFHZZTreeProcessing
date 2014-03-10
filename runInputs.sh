#!/bin/bash

pathTo7TeV=../Histogramming/rootFiles_Legacy
pathTo8TeV=/scratch/osghpc/snowball/UF_HZZ4L_Analysis/FullAnalysis2012/MEKD_Z4L/Histogramming/rootFiles

lumi=19.72

root -l -b -q "makeInput.C(\"${pathTo8TeV}/Data.root\",$lumi,1)"
root -l -b -q "makeInput.C(\"${pathTo8TeV}/Data.root\",$lumi,2)"
root -l -b -q "makeInput.C(\"${pathTo8TeV}/Data.root\",$lumi,3)"

lumi=5.051

#root -l -b -q "makeInput.C(\"${pathTo7TeV}/Data.root\",$lumi,1)"
#root -l -b -q "makeInput.C(\"${pathTo7TeV}/Data.root\",$lumi,2)"
#root -l -b -q "makeInput.C(\"${pathTo7TeV}/Data.root\",$lumi,3)"
