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
#include "processTrees.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraphPolar.h"
#include "TSystem.h"
#include "TPaletteAxis.h"

using namespace std;



//Scale factors for data/MC efficiency                                                                                                                     
TFile *fMuWeight = new TFile(pathToScaleFactors+"/MuonScaleFactors_2011_2012.root","READ");
TFile *fMuWeight12 = new TFile(pathToScaleFactors+"/SF2013_V5_old2011.root","READ");
//TFile *fElWeight = new TFile(pathToScaleFactors+"/scale_factors_ele2011.root","READ");
TFile *fElWeight = new TFile(pathToScaleFactors+"/CombinedMethod_ScaleFactors_RecoIdIsoSip_2011.root","READ");
TFile *fElWeight12 = new TFile(pathToScaleFactors+"/scale_factors_ele2012.root","READ");

TH2D *hTH2D_Mu_All_2011A = (TH2D*)fMuWeight->Get("TH2D_ALL_2011A");
TH2D *hTH2D_Mu_All_2011B = (TH2D*)fMuWeight->Get("TH2D_ALL_2011B");
TH2D *hTH2D_Mu_All_2012  = (TH2D*)fMuWeight12->Get("TH2D_ALL_2012");

//TH2D *hTH2D_El_All_2011A = (TH2D*)fElWeight->Get("h_electron_scale_factor_RECO_ID_ISO_SIP");
//TH2D *hTH2D_El_All_2011B = (TH2D*)fElWeight->Get("h_electron_scale_factor_RECO_ID_ISO_SIP");
TH2D *hTH2D_El_All_2011A = (TH2D*)fElWeight->Get("h_electronScaleFactor_RecoIdIsoSip");
TH2D *hTH2D_El_All_2011B = (TH2D*)fElWeight->Get("h_electronScaleFactor_RecoIdIsoSip");
//TH2D *hTH2D_El_All_2012  = (TH2D*)fElWeight12->Get("h_electron_scale_factor_RECO_ID_ISO_SIP");

//TH2D *hTH2D_El_All_2011A = (TH2D*)fElWeight->Get("heff_electron_selection");
//TH2D *hTH2D_El_All_2011B = (TH2D*)fElWeight->Get("heff_electron_selection");
TH2D *hTH2D_El_All_2012  = (TH2D*)fElWeight12->Get("h_electronScaleFactor_RecoIdIsoSip");



int main(){

  gSystem->Load("libHistPainter");  

  if(fMuWeight == NULL){ cout << "Cannot open " << fMuWeight->GetName() << endl; return -1;}
  if(fMuWeight12 == NULL){ cout << "Cannot open " << fMuWeight12->GetName() <<endl; return -1;}
  if(fElWeight == NULL){ cout << "Cannot open " << fElWeight->GetName() <<endl; return -1;}
  if(fElWeight12 == NULL){ cout << "Cannot open " << fElWeight12->GetName() <<endl; return -1;}


  for(int i = 0; i < nFiles7; i++)
    {
      TString newFileName = filePath7TeV+"/"+files7[i];
      addDataMCScalingFactors(newFileName,treeDirectory7,treeName7,2011);
    }

  
  for(int i = 0; i < nFiles8; i++)
    {
      TString newFileName = filePath8TeV+"/"+files8[i];
      addDataMCScalingFactors(newFileName,treeDirectory8,treeName8,2012);
    }

  fMuWeight->Close();
  fElWeight->Close();
  fMuWeight12->Close();
  fElWeight12->Close();

  return 0;
}




void addDataMCScalingFactors(TString filename,TString treeDir, TString treename, int myYear)
{

  int idL1,idL2,idL3,idL4;
  double pTL1,pTL2,pTL3,pTL4;
  double etaL1,etaL2,etaL3,etaL4;
  double phiL1,phiL2,phiL3,phiL4;

  cout << filename << endl;

  TFile *file = new TFile(filename,"UPDATE");
  
  TTree *oldtree = (TTree*)file->Get(treeDir+"/"+treename);
  if(oldtree==NULL)
    {
      cout << "Could not find tree " << treeDir << "/" << treename << endl
	   << "in file " << file->GetName() << endl;
      return;
    }
  
  oldtree->SetBranchAddress("idL1",&idL1);
  oldtree->SetBranchAddress("idL2",&idL2);
  oldtree->SetBranchAddress("idL3",&idL3);
  oldtree->SetBranchAddress("idL4",&idL4);
  oldtree->SetBranchAddress("pTL1",&pTL1);
  oldtree->SetBranchAddress("pTL2",&pTL2);
  oldtree->SetBranchAddress("pTL3",&pTL3);
  oldtree->SetBranchAddress("pTL4",&pTL4);
  oldtree->SetBranchAddress("etaL1",&etaL1);
  oldtree->SetBranchAddress("etaL2",&etaL2);
  oldtree->SetBranchAddress("etaL3",&etaL3);
  oldtree->SetBranchAddress("etaL4",&etaL4);
  oldtree->SetBranchAddress("phiL1",&phiL1);
  oldtree->SetBranchAddress("phiL2",&phiL2);
  oldtree->SetBranchAddress("phiL3",&phiL3);
  oldtree->SetBranchAddress("phiL4",&phiL4);

  TH1F *hNevents = (TH1F*)file->Get(treeDir+"/nEvents");
  double dataMCweight = 1, nEvents = 0;
  TBranch *branch = oldtree->Branch("dataMC_weight",&dataMCweight,"dataMC_weight/D");
  TBranch *nevbranch = oldtree->Branch("nEvents",&nEvents,"nEvents/D");
  double nEventsTmp = hNevents->GetBinContent(1);  

  double pts[4],phis[4],etas[4];
  int lepIDs[4];

  for(int i = 0; i < oldtree->GetEntries(); i++)
    {
      oldtree->GetEntry(i);
      
      /////////////////
      lepIDs[0] = idL1;
      lepIDs[1] = idL2;
      lepIDs[2] = idL3;
      lepIDs[3] = idL4;

      pts[0] = pTL1;
      pts[1] = pTL2;
      pts[2] = pTL3;
      pts[3] = pTL4;

      phis[0] = phiL1;
      phis[1] = phiL2;
      phis[2] = phiL3;
      phis[3] = phiL4;

      etas[0] = etaL1;
      etas[1] = etaL2;
      etas[2] = etaL3;
      etas[3] = etaL4;
      /////////////////

      dataMCweight = getWeight(myYear,etas,phis,pts,lepIDs);
      nEvents = nEventsTmp;
      branch->Fill();
      nevbranch->Fill();
    }

  file->cd();
  oldtree->CloneTree()->Write(treename+"_dataMC", TObject::kOverwrite);
  file->Close();
  
}




double getAllWeight(double LepPt,double LepEta,int year, int LepID)
{

  int scaleFactYear = -1;
  double Run2011AFraction = 0.465;

  if( year == 2011 )
    {
      TRandom3 rand;
      double whatPeriod = rand.Rndm(1234); 
      if(whatPeriod < Run2011AFraction) scaleFactYear = 0;
      else if(whatPeriod > Run2011AFraction) scaleFactYear = 1;
      else {cout << "Random number problem!!" << endl; return -1;}
    }
  else if( year == 2012 ){scaleFactYear = 2;}
  else {cout << "Unknown year " << year << endl; return -1;}

  double weight  = 1.; 
  double errCorr = 0.;

  double myLepPt = LepPt;
  double myLepEta = LepEta;
  double myLepID = abs(LepID);
  
  //avoid to go out of the TH boundary
  if(abs(LepID) == 13 && myLepPt > 99.) myLepPt = 99.;
  if(abs(LepID) == 11 && myLepPt > 199.) myLepPt = 199.;
  if(abs(LepID) == 11) myLepEta = fabs(myLepEta);

  if(scaleFactYear == 0)
    {
      if( myLepID == 13)
	{                                               
	  weight  = hTH2D_Mu_All_2011A->GetBinContent(hTH2D_Mu_All_2011A->GetXaxis()->FindBin(myLepPt),hTH2D_Mu_All_2011A->GetYaxis()->FindBin(LepEta));
	  errCorr = hTH2D_Mu_All_2011A->GetBinError(hTH2D_Mu_All_2011A->GetXaxis()->FindBin(myLepPt),hTH2D_Mu_All_2011A->GetYaxis()->FindBin(LepEta));
	  //if (weight == 0) cout << "muon " << myLepPt << "   " << LepEta << endl;

	}
      else if( myLepID == 11)
	{   
	  weight  = hTH2D_El_All_2011A->GetBinContent(hTH2D_El_All_2011A->GetXaxis()->FindBin(myLepPt),hTH2D_El_All_2011A->GetYaxis()->FindBin(myLepEta));
	  errCorr = hTH2D_El_All_2011A->GetBinError(hTH2D_El_All_2011A->GetXaxis()->FindBin(myLepPt),hTH2D_El_All_2011A->GetYaxis()->FindBin(myLepEta));  
	  //if (weight ==0) cout<< "elec " << myLepPt << "   " << LepEta << endl; 
	}
      else {
	return -1;
      }    
    }
  else if(scaleFactYear == 1)
    {
      if( myLepID == 13)
	{                                               
	  weight  = hTH2D_Mu_All_2011B->GetBinContent(hTH2D_Mu_All_2011B->GetXaxis()->FindBin(myLepPt),hTH2D_Mu_All_2011B->GetYaxis()->FindBin(LepEta));
	  errCorr = hTH2D_Mu_All_2011B->GetBinError(hTH2D_Mu_All_2011B->GetXaxis()->FindBin(myLepPt),hTH2D_Mu_All_2011B->GetYaxis()->FindBin(LepEta));
	  //if (weight == 0) cout << "muon " << myLepPt << "   " << LepEta << endl;
	}
      else if( myLepID == 11)
	{   
	  weight  = hTH2D_El_All_2011B->GetBinContent(hTH2D_El_All_2011B->GetXaxis()->FindBin(myLepPt),hTH2D_El_All_2011B->GetYaxis()->FindBin(myLepEta));
	  errCorr = hTH2D_El_All_2011B->GetBinError(hTH2D_El_All_2011B->GetXaxis()->FindBin(myLepPt),hTH2D_El_All_2011B->GetYaxis()->FindBin(myLepEta));
	  //if (weight ==0) cout<< "elec " << myLepPt << "   " << LepEta << endl;
	}
      else {
	return -1;
      }    
    }
  else if(scaleFactYear == 2)
    {
      if( myLepID == 13)
	{
	  weight  = hTH2D_Mu_All_2012->GetBinContent(hTH2D_Mu_All_2012->GetXaxis()->FindBin(myLepPt),hTH2D_Mu_All_2012->GetYaxis()->FindBin(LepEta));
	  errCorr = hTH2D_Mu_All_2012->GetBinError(hTH2D_Mu_All_2012->GetXaxis()->FindBin(myLepPt),hTH2D_Mu_All_2012->GetYaxis()->FindBin(LepEta));
	  //if (weight == 0) cout << "muon " << myLepPt << "   " << LepEta << endl;
	}
      else if( myLepID == 11)
	{
	  weight  = hTH2D_El_All_2012->GetBinContent(hTH2D_El_All_2012->GetXaxis()->FindBin(myLepPt),hTH2D_El_All_2012->GetYaxis()->FindBin(myLepEta));
	  errCorr = hTH2D_El_All_2012->GetBinError(hTH2D_El_All_2012->GetXaxis()->FindBin(myLepPt),hTH2D_El_All_2012->GetYaxis()->FindBin(myLepEta));
          //if (weight ==0) cout<< "elec " << myLepPt << "   " << LepEta << endl;
	}
      else {
	return -1;
      }
    }
  else return -1;
  
  //add the systematics on T&P corrections
  if(myLepID == 13)
    {
      if(myLepPt >= 15.) errCorr = sqrt(errCorr*errCorr + 0.005*0.005);
      else{
	if(myLepEta < 1.2) errCorr = sqrt(errCorr*errCorr + 0.01*0.01);
	else errCorr = sqrt(errCorr*errCorr + 0.014*0.014);
      }
    }

  //Safety
  if(myLepPt < 5. && myLepID == 13) weight = 1;
  if(myLepPt < 7. && myLepID == 11) weight = 1;
  
  //if(weight < 0.001 || weight > 10.)
  // {
  //  cout << "pt = " << myLepPt << " eta = " << myLepEta << " weight = " << weight << endl;
  //  return -1;
  //}
  
  int Seeder = 0;//always 0 for system studies
  TRandom3 randomToss;
  if( (applySystMu && myLepID == 13) || (applySystEle && myLepID == 11) ){
    weight = randomToss.Gaus(weight,errCorr);
    randomToss.SetSeed(Seeder);
  }

  return weight;


}


double getWeight(int year,double lepEtas[],double lepPhis[], double lepPts[],int lepIDs[])
{
  double eff_weight = 1.;
  
  for(int nLep=0;nLep<4;nLep++)
    {
      if(fabs(lepIDs[nLep]) == 13) eff_weight *= getAllWeight(lepPts[nLep], lepEtas[nLep], year, lepIDs[nLep]);
      if(fabs(lepIDs[nLep]) == 11) eff_weight *= getAllWeight(lepPts[nLep], lepEtas[nLep], year, lepIDs[nLep]);
    }

  return eff_weight;
  
}





