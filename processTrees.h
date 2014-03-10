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

const TString filePath7TeV = "/scratch/osghpc/snowball/UF_HZZ4L_Analysis/FullAnalysis2012/Histogramming/rootFiles_Legacy_dataMC/";
//const TString filePath8TeV = "/scratch/osghpc/snowball/UF_HZZ4L_Analysis/FullAnalysis2012/Histogramming_8TeV/rootFiles_Legacy_dataMC/";
const TString filePath8TeV = "/scratch/osghpc/snowball/UF_HZZ4L_Analysis/FullAnalysis2012/MEKD_Z4L/Histogramming/rootFiles/";

const TString treeDirectory7 = "AnaAfterHlt";
const TString treeDirectory8 = "AnaAfterHlt";

const TString treeName7 = "passedEvents";
const TString treeName8 = "passedEvents";

const int nFiles7 = 0;//74
const int nFiles8 = 10;//116

const TString files7[nFiles7] = {};

/* const TString files7[nFiles7] = { */


/*   "ggZZ_2e2mu.root", */
/*   "ggZZ_4l.root", */
/*   "Graviton2MHToZZTo4L_M-126_7TeV-JHUGenV2-pythia6_7TeV.root", */
/*   "Graviton2PBToZZTo4L_M-126_7TeV-JHUGenV2-pythia6_7TeV.root", */
/*   "Graviton2PHToZZTo4L_M-126_7TeV-JHUGenV2-pythia6_7TeV.root", */
/*   "Graviton2PMqqbarToZZTo4L_M-126_7TeV-JHUgenV2-PYTHIA6_Tauola_7TeV.root", */
/*   "Graviton2PMToZZTo4L_M-126_7TeV_ext-JHUgenV2-PYTHIA6_Tauola_7TeV.root", */
/*   "Higgs0MToZZTo4L_M-126_7TeV_ext-JHUgenV2-pythia6_7TeV.root", */
/*   "Higgs0PHToZZTo4L_M-126_7TeV_ext-JHUgenV2-pythia6_7TeV.root", */
/*   "mH_115_ggH.root", */
/*   "mH_115_VBF.root", */
/*   "mH_120_ggH.root", */
/*   "mH_120_VBF.root", */
/*   "mH_130_ggH.root", */
/*   "mH_130_VBF.root", */
/*   "mH_140_ggH.root", */
/*   "mH_140_VBF.root", */
/*   "mH_150_ggH.root", */
/*   "mH_150_VBF.root", */
/*   "mH_160_ggH.root", */
/*   "mH_160_VBF.root", */
/*   "mH_170_ggH.root", */
/*   "mH_170_VBF.root", */
/*   "mH_180_ggH.root", */
/*   "mH_180_VBF.root", */
/*   "mH_190_ggH.root", */
/*   "mH_190_VBF.root", */
/*   "mH_200_ggH.root", */
/*   "mH_200_VBF.root", */
/*   "mH_210_ggH.root", */
/*   "mH_210_VBF.root", */
/*   "mH_220_ggH.root", */
/*   "mH_220_VBF.root", */
/*   "mH_230_ggH.root", */
/*   "mH_230_VBF.root", */
/*   "mH_250_ggH.root", */
/*   "mH_250_VBF.root", */
/*   "mH_275_ggH.root", */
/*   "mH_275_VBF.root", */
/*   "mH_300_ggH.root", */
/*   "mH_300_VBF.root", */
/*   "mH_325_ggH.root", */
/*   "mH_325_VBF.root", */
/*   "mH_350_ggH.root", */
/*   "mH_350_VBF.root", */
/*   "mH_375_ggH.root", */
/*   "mH_375_VBF.root", */
/*   "mH_400_ggH.root", */
/*   "mH_400_VBF.root", */
/*   "mH_425_ggH.root", */
/*   "mH_425_VBF.root", */
/*   "mH_450_ggH.root", */
/*   "mH_450_VBF.root", */
/*   "mH_475_ggH.root", */
/*   "mH_475_VBF.root", */
/*   "mH_500_ggH.root", */
/*   "mH_500_VBF.root", */
/*   "mH_525_ggH.root", */
/*   "mH_525_VBF.root", */
/*   "mH_550_ggH.root", */
/*   "mH_550_VBF.root", */
/*   "mH_575_ggH.root", */
/*   "mH_575_VBF.root", */
/*   "mH_600_ggH.root", */
/*   "mH_600_VBF.root", */
/*   "SMHiggsToZZTo4L_M-126_7TeV_ext-JHUgenV2-PYTHIA6_Tauola_7TeV.root", */
/*   "Vector1MToZZTo4L_M-126_7TeV-JHUgenV2-PYTHIA6_Tauola_7TeV.root", */
/*   "Vector1PToZZTo4L_M-126_7TeV-JHUgenV2-PYTHIA6_Tauola_7TeV.root", */
/*   "ZZ_2e2mu.root", */
/*   "ZZ_4e.root", */
/*   "ZZ_4mu.root", */
/*   "ZZto2e2tau.root", */
/*   "ZZto2mu2tau.root", */
/*   "ZZto4tau.root" */
/* }; */


const TString files8[nFiles8] =  {
  //   "Higgs0MToZZTo4L_M-126_8TeV_ext-JHUgenV2-pythia6_8TeV.root",
   "ZZTo2e2mu.root",
   "ZZTo2e2tau.root",
   "ZZTo2mu2tau.root",
   "ZZTo4e.root",
   "ZZTo4mu.root",
   "ZZTo4tau.root",
   "ggH_91.2.root",
   //   "mH_126_ggH.root",
   "ZZ2e2mu_tchan.root",
   "ZZ4e_tchan.root",
   "ZZ4mu_tchan.root"
};
  /*
const TString files8[nFiles8] =  { 
  "ggZZ_2e2mu.root",
  "ggZZ_4l.root",
  "Graviton2MHToZZTo4L_M-126_8TeV-JHUGenV2-pythia6_8TeV.root",
  "Graviton2PBToZZTo4L_M-126_8TeV-JHUGenV2-pythia6_8TeV.root",
  "Graviton2PHToZZTo4L_M-126_8TeV-JHUGenV2-pythia6_8TeV.root",
  "Graviton2PMqqbarToZZTo4L_M-126_8TeV-JHUgenV2-PYTHIA6_Tauola_8TeV.root",
  "Graviton2PMToZZTo4L_M-126_8TeV_ext-JHUgenV2-PYTHIA6_Tauola_8TeV.root",
  "Higgs0MToZZTo4L_M-126_8TeV_ext-JHUgenV2-pythia6_8TeV.root",
  "Higgs0PHToZZTo4L_M-126_8TeV-JHUgenV2-pythia6_8TeV.root",
  "mH_1000_ggH.root",
  "mH_1000_VBF.root",
  "mH_115_ggH.root",
  "mH_115_VBF.root",
  "mH_116_ggH.root",
  "mH_116_VBF.root",
  "mH_117_ggH.root",
  "mH_117_VBF.root",
  "mH_118_ggH.root",
  "mH_118_VBF.root",
  "mH_119_ggH.root",
  "mH_119_VBF.root",
  "mH_120_ggH.root",
  "mH_120_VBF.root",
  "mH_121_ggH.root",
  "mH_121_VBF.root",
  "mH_122_ggH.root",
  "mH_122_VBF.root",
  "mH_123_ggH.root",
  "mH_123_VBF.root",
  "mH_124_ggH.root",
  "mH_124_VBF.root",
  "mH_125_ggH.root",
  "mH_125_VBF.root",
  "mH_126_ggH.root",
  "mH_126_VBF.root",
  "mH_127_ggH.root",
  "mH_127_VBF.root",
  "mH_128_ggH.root",
  "mH_128_VBF.root",
  "mH_129_ggH.root",
  "mH_129_VBF.root",
  "mH_130_ggH.root",
  "mH_130_VBF.root",
  "mH_135_ggH.root",
  "mH_135_VBF.root",
  "mH_140_ggH.root",
  "mH_140_VBF.root",
  "mH_145_ggH.root",
  "mH_145_VBF.root",
  "mH_150_ggH.root",
  "mH_150_VBF.root",
  "mH_160_ggH.root",
  "mH_160_VBF.root",
  "mH_170_ggH.root",
  "mH_170_VBF.root",
  "mH_180_ggH.root",
  "mH_180_VBF.root",
  "mH_190_ggH.root",
  "mH_190_VBF.root",
  "mH_200_ggH.root",
  "mH_200_VBF.root",
  "mH_220_ggH.root",
  "mH_220_VBF.root",
  "mH_250_ggH.root",
  "mH_250_VBF.root",
  "mH_275_ggH.root",
  "mH_275_VBF.root",
  "mH_300_ggH.root",
  "mH_300_VBF.root",
  "mH_325_ggH.root",
  "mH_325_VBF.root",
  "mH_350_ggH.root",
  "mH_350_VBF.root",
  "mH_375_ggH.root",
  "mH_375_VBF.root",
  "mH_400_ggH.root",
  "mH_400_VBF.root",
  "mH_425_ggH.root",
  "mH_425_VBF.root",
  "mH_450_ggH.root",
  "mH_450_VBF.root",
  "mH_475_ggH.root",
  "mH_475_VBF.root",
  "mH_500_ggH.root",
  "mH_500_VBF.root",
  "mH_525_ggH.root",
  "mH_525_VBF.root",
  "mH_550_ggH.root",
  "mH_550_VBF.root",
  "mH_575_ggH.root",
  "mH_575_VBF.root",
  "mH_600_ggH.root",
  "mH_600_VBF.root",
  "mH_650_ggH.root",
  "mH_650_VBF.root",
  "mH_700_ggH.root",
  "mH_700_VBF.root",
  "mH_750_ggH.root",
  "mH_750_VBF.root",
  "mH_800_ggH.root",
  "mH_800_VBF.root",
  "mH_850_ggH.root",
  "mH_850_VBF.root",
  "mH_900_ggH.root",
  "mH_900_VBF.root",
  "mH_950_ggH.root",
  "mH_950_VBF.root",
  "SMHiggsToZZTo4L_M-126_8TeV_ext-JHUgenV2-pythia6_8TeV.root",
  "Vector1MToZZTo4L_M-126_8TeV-JHUgenV2-PYTHIA6_Tauola_8TeV.root",
  "Vector1PToZZTo4L_M-126_8TeV-JHUgenV2-PYTHIA6_Tauola_8TeV.root",
  "ZZ_2e2mu.root",
  "ZZ_4e.root",
  "ZZ_4mu.root",
  "ZZto2e2tau.root",
  "ZZto2mu2tau.root",
  "ZZto4tau.root"
  
  };


  */



const bool applySystMu = false;
const bool applySystEle = false;
