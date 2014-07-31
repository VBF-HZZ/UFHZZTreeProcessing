#!/bin/bash

#pathTo7TeV=../Histogramming/rootFiles_Legacy
#pathTo8TeV=/scratch/osghpc/snowball/UF_HZZ4L_Analysis/FullAnalysis2012/MEKD_Z4L/Histogramming/rootFiles
pathTo7TeV=/scratch/osg/lihengne/cms/csa14/analyzer8tev/treeProcess/UFHZZTreeProcessing/rootFiles_7TeV
pathTo8TeV=/scratch/osg/lihengne/cms/csa14/analyzer8tev/treeProcess/UFHZZTreeProcessing/rootFiles_8TeV

lumi=19.7

root -l -b -q "makeInput.C(\"${pathTo8TeV}/Data_8TeV.root\",$lumi,1)"
root -l -b -q "makeInput.C(\"${pathTo8TeV}/Data_8TeV.root\",$lumi,2)"
root -l -b -q "makeInput.C(\"${pathTo8TeV}/Data_8TeV.root\",$lumi,3)"

lumi=5.1

root -l -b -q "makeInput.C(\"${pathTo7TeV}/Data_7TeV.root\",$lumi,1)"
root -l -b -q "makeInput.C(\"${pathTo7TeV}/Data_7TeV.root\",$lumi,2)"
root -l -b -q "makeInput.C(\"${pathTo7TeV}/Data_7TeV.root\",$lumi,3)"
