#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <boost/regex.hpp>

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TBranch.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TString.h>



void addDataMCScalingFactors(TString filename,TString treeDir, TString treename, int myYear);
double getAllWeight(double LepPt, double LepEta, int year, int LepID);
double getWeight( int year,double lepEtas[],double lepPhis[], double lepPts[],int lepIDs[]);

const TString pathToScaleFactors = "scaleFactors/";

const TString filePath7TeV = "/scratch/osg/lihengne/cms/csa14/analyzer8tev/treeProcess/UFHZZTreeProcessing/rootFiles_7TeV";
const TString filePath8TeV = "/scratch/osg/lihengne/cms/csa14/analyzer8tev/treeProcess/UFHZZTreeProcessing/rootFiles_8TeV";

const TString treeDirectory7 = "AnaAfterHlt";
const TString treeDirectory8 = "AnaAfterHlt";

const TString treeName7 = "passedEvents";
const TString treeName8 = "passedEvents";

//const int nFiles7 = 10;
const int nFiles7 = 0;
//const int nFiles8 = 13;
const int nFiles8 = 1;


const TString files7[nFiles7];
//const TString files7[nFiles7] = {
//"ZZ_2e2mu.root",
//"ZZ_4e.root",
//"ZZ_4mu.root",
//"ZZto2e2tau.root",
//"ZZto2mu2tau.root",
//"ZZto4tau.root",
//"ggZZ_2e2mu.root",
//"ggZZ_4l.root",
//"mH_350_VBF.root",
//"mH_350_ggH.root"
// }; 


const TString files8[nFiles8] = {
//"ZZTo2e2mu_8TeV_ext.root"
//"ZZTo2e2tau_8TeV_ext.root",
//"ZZTo2mu2tau_8TeV_ext.root",
//"ZZTo4e_8TeV_ext.root",
//"ZZTo4mu_8TeV_ext.root",
//"ZZTo4tau_8TeV_ext.root",
//"mH_125_TTH.root",
"mH_126_WH.root"//,
//"mH_125_ZH.root",
//"mH_126_SMH.root",
//"mH_126_TTH.root",
//"mH_126_WH.root",
//"mH_126_ZH.root"
//"ZZJetsTo4L.root"
};
//const TString files8[nFiles8] =  {
//"ZZ_2e2mu.root",
//"ZZ_4e.root",
//"ZZ_4mu.root",
//"ZZto2e2tau.root",
//"ZZto2mu2tau.root",
//"ZZto4tau.root",
//"ggZZ_2e2mu.root",
//"ggZZ_4l.root",
//"mH_126_VBF.root",
//"mH_126_ggH.root",
//"mH_126_ggH_powheg15.root",
//"mH_350_VBF.root",
//"mH_350_ggH.root"
//};

const bool applySystMu = false;
const bool applySystEle = false;
