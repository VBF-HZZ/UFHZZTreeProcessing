#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TMath.h>
#include <TROOT.h>
#include "TRandom.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <stdio.h>
#include <math.h>

using namespace std;

int main() {
  makeInput();
  return 0;
}

void makeInput(TString file, double lumi,int finalState)
{

  TFile *f = new TFile(file,"READ");
  TTree *oldtree = (TTree*)f->Get("AnaAfterHlt/passedEvents");
  if(oldtree==NULL)oldtree = (TTree*)f->Get("Ana/passedEvents");
  if(oldtree==NULL)oldtree = (TTree*)f->Get("AnaCorrector/passedEvents");

  if(!f || !oldtree){cout << "Cannot find file or tree" << endl; return;}

  double mzz,m4lErr,mekdLD,mekdLDNew,mz1,mz2;
  double melaLD, mekdSig, mekdBg;
  int idL1,idL2,idL3,idL4, _fs;
  double pTL1,pTL2,pTL3,pTL4;
  double pXL1,pXL2,pXL3,pXL4;
  double pYL1,pYL2,pYL3,pYL4;
  double pZL1,pZL2,pZL3,pZL4;
  double EL1,EL2,EL3,EL4;
  ULong64_t Run, Event, LumiSect;
  bool passedQCDSelection;
  double pseudoKD, graviKD, p0hplusKD;
  double p1plusKD, p1minusKD, qqgraviKD;
  double PIp1plusKD, PIp1minusKD;
  double PIp2plusKD, p2bplusKD, p2hplusKD, p2hminusKD;
  double me_qqZZtchan, me_H, me_qqZZschan, me_qqZZschan2;
  double mass4lNoFSR;
  bool FSRZ1, FSRZ2;
  double Z4l_minMass2l;

  std::vector<unsigned long> runVec, lumiVec, eventVec;

  oldtree->SetBranchAddress("mass4l",&mzz);
  oldtree->SetBranchAddress("mass4lNoFSR",&mass4lNoFSR);
  oldtree->SetBranchAddress("massZ1",&mz1);
  oldtree->SetBranchAddress("massZ2",&mz2);
  oldtree->SetBranchAddress("FSR_Z1",&FSRZ1);
  oldtree->SetBranchAddress("FSR_Z2",&FSRZ2);
  oldtree->SetBranchAddress("Z4l_minMass2l",&Z4l_minMass2l);


  oldtree->SetBranchAddress("JHUKD_H_qqZZ_noPDF",&melaLD);
  oldtree->SetBranchAddress("MEKD_noPDF_noFSR",&mekdLD);
  oldtree->SetBranchAddress("MEKD_ME_H_noPDF_noFSR",&mekdSig);
  oldtree->SetBranchAddress("MEKD_ME_ZZ_noPDF_noFSR",&mekdBg);

  
  oldtree->SetBranchAddress("JHUKD_H_h0M_noPDF",&pseudoKD);
  oldtree->SetBranchAddress("JHUKD_H_h0P_noPDF",&p0hplusKD);
  oldtree->SetBranchAddress("JHUKD_H_h1M_noPDF",&p1minusKD);
  oldtree->SetBranchAddress("JHUKD_H_h1P_noPDF",&p1plusKD);
  oldtree->SetBranchAddress("JHUKD_H_ggh2P_noPDF",&graviKD);
  oldtree->SetBranchAddress("JHUKD_H_qqh2P_noPDF",&qqgraviKD);
  oldtree->SetBranchAddress("JHUKD_H_h2hP_noPDF",&p2hplusKD);
  oldtree->SetBranchAddress("JHUKD_H_h2hM_noPDF",&p2hminusKD);
  oldtree->SetBranchAddress("JHUKD_H_h2bP_noPDF",&p2bplusKD);
  oldtree->SetBranchAddress("JHUKD_H_h2P_prodInd_noPDF",&PIp2plusKD);
  oldtree->SetBranchAddress("JHUKD_H_h1P_prodInd_noPDF",&PIp1plusKD);
  oldtree->SetBranchAddress("JHUKD_H_h1M_prodInd_noPDF",&PIp1minusKD);

  oldtree->SetBranchAddress("MEKD_ME_H_noPDF",&me_H);
  oldtree->SetBranchAddress("MEKD_ME_qqZ4l_bkg_noPDF",&me_qqZZtchan);
  oldtree->SetBranchAddress("MEKD_ME_qqZ4l_sig_noPDF",&me_qqZZschan);
  oldtree->SetBranchAddress("MEKD_ME_qqZ4l_sig_noPDF_noFSR",&me_qqZZschan2);
  
  /*
  oldtree->SetBranchAddress("MEKD_h0M_ZZ_noPDF_noFSR",&pseudoKD);
  oldtree->SetBranchAddress("MEKD_h0P_ZZ_noPDF_noFSR",&p0hplusKD);
  oldtree->SetBranchAddress("MEKD_h1M_ZZ_noPDF_noFSR",&p1minusKD);
  oldtree->SetBranchAddress("MEKD_h1P_ZZ_noPDF_noFSR",&p1plusKD);
  oldtree->SetBranchAddress("MEKD_qqh2P_ZZ_noPDF_noFSR",&qqgraviKD);
  oldtree->SetBranchAddress("MEKD_ggh2P_ZZ_noPDF_noFSR",&graviKD);
  oldtree->SetBranchAddress("MEKD_h2hP_ZZ_noPDF_noFSR",&p2hplusKD);
  oldtree->SetBranchAddress("MEKD_h2hM_ZZ_noPDF_noFSR",&p2hminusKD);
  oldtree->SetBranchAddress("MEKD_h2bP_ZZ_noPDF_noFSR",&p2bplusKD);
  oldtree->SetBranchAddress("MEKD_h2P_ZZ_prodInd_noPDF_noFSR",&PIp2plusKD);
  oldtree->SetBranchAddress("MEKD_h1P_ZZ_prodInd_noPDF_noFSR",&PIp1plusKD);
  oldtree->SetBranchAddress("MEKD_h1M_ZZ_prodInd_noPDF_noFSR",&PIp1minusKD);
  */

  oldtree->SetBranchAddress("finalState",&_fs);
  oldtree->SetBranchAddress("massErrorUFCorr",&m4lErr);
  oldtree->SetBranchAddress("idL1",&idL1);
  oldtree->SetBranchAddress("idL2",&idL2);
  oldtree->SetBranchAddress("idL3",&idL3);
  oldtree->SetBranchAddress("idL4",&idL4);
  oldtree->SetBranchAddress("pTL1",&pTL1);
  oldtree->SetBranchAddress("pTL2",&pTL2);
  oldtree->SetBranchAddress("pTL3",&pTL3);
  oldtree->SetBranchAddress("pTL4",&pTL4);

  oldtree->SetBranchAddress("pXL1",&pXL1);
  oldtree->SetBranchAddress("pXL2",&pXL2);
  oldtree->SetBranchAddress("pXL3",&pXL3);
  oldtree->SetBranchAddress("pXL4",&pXL4);

  oldtree->SetBranchAddress("pYL1",&pYL1);
  oldtree->SetBranchAddress("pYL2",&pYL2);
  oldtree->SetBranchAddress("pYL3",&pYL3);
  oldtree->SetBranchAddress("pYL4",&pYL4);

  oldtree->SetBranchAddress("pZL1",&pZL1);
  oldtree->SetBranchAddress("pZL2",&pZL2);
  oldtree->SetBranchAddress("pZL3",&pZL3);
  oldtree->SetBranchAddress("pZL4",&pZL4);

  oldtree->SetBranchAddress("EL1",&EL1);
  oldtree->SetBranchAddress("EL2",&EL2);
  oldtree->SetBranchAddress("EL3",&EL3);
  oldtree->SetBranchAddress("EL4",&EL4);

  oldtree->SetBranchAddress("Run",&Run);
  oldtree->SetBranchAddress("Event",&Event);
  oldtree->SetBranchAddress("LumiSect",&LumiSect);
  oldtree->SetBranchAddress("passedQCDcut",&passedQCDSelection);


  double globalPtMuCut = 5, globalPtElCut = 7;
  double globalMZ1Low = 40, globalMZ2Low = 4, globalM4lCut = 40;

  std::string fsn;
  if (finalState == 1) fsn = "4mu";
  if (finalState == 2) fsn = "4e";
  if (finalState == 3) fsn = "2e2mu";
  

  char Name[192];
  if( lumi > 10) sprintf(Name,"StatisticalInputTrees/hzz%s_%.2f.root",fsn.c_str(),lumi);
  else sprintf(Name,"StatisticalInputTrees/hzz%s_%.3f.root",fsn.c_str(),lumi);


  TFile *hzzF = new TFile(Name,"RECREATE");
  TTree *tree = new TTree("data_obs","data_obs");

  double mass, massErr, mela, mekd, mekdLog;
  double CMS_zz4l_pseudoKD, CMS_zz4l_graviKD, CMS_zz4l_p0hplusKD;
  double CMS_zz4l_p1plusKD, CMS_zz4l_p1minusKD, CMS_zz4l_qqgraviKD;
  double CMS_zz4l_PIp1plusKD, CMS_zz4l_PIp1minusKD;
  double CMS_zz4l_PIp2plusKD, CMS_zz4l_p2bplusKD;
  double CMS_zz4l_p2hplusKD, CMS_zz4l_p2hminusKD;


  tree->Branch("CMS_zz4l_mass",&mass,"CMS_zz4l_mass/D");
  //tree->Branch("CMS_zz4l_massRelErr",&massErr,"CMS_zz4l_massRelErr/D");
  //tree->Branch("melaLD",&mela,"melaLD/D");
  //tree->Branch("mekdLD",&mekd,"mekdLD/D");
  //tree->Branch("mekdLLD",&mekdLog,"mekdLLD/D");

  tree->Branch("CMS_zz4l_smd",&mela,"CMS_zz4l_smd/D");
  tree->Branch("CMS_zz4l_pseudoKD",&CMS_zz4l_pseudoKD,"CMS_zz4l_pseudoKD/D");
  //tree->Branch("CMS_zz4l_p0hplusKD",&CMS_zz4l_p0hplusKD,"CMS_zz4l_p0hplusKD/D");
  //tree->Branch("CMS_zz4l_p1plusKD",&CMS_zz4l_p1plusKD,"CMS_zz4l_p1plusKD/D");
  //tree->Branch("CMS_zz4l_p1minusKD",&CMS_zz4l_p1minusKD,"CMS_zz4l_p1minusKD/D");
  //tree->Branch("CMS_zz4l_qqgraviKD",&CMS_zz4l_qqgraviKD,"CMS_zz4l_qqgraviKD/D");
  //tree->Branch("CMS_zz4l_graviKD",&CMS_zz4l_graviKD,"CMS_zz4l_graviKD/D");
  //tree->Branch("CMS_zz4l_PIp1plusKD",&CMS_zz4l_PIp1plusKD,"CMS_zz4l_PIp1plusKD/D");
  //tree->Branch("CMS_zz4l_PIp1minusKD",&CMS_zz4l_PIp1minusKD,"CMS_zz4l_PIp1minusKD/D");
  //tree->Branch("CMS_zz4l_PIp2plusKD",&CMS_zz4l_PIp2plusKD,"CMS_zz4l_PIp2plusKD/D");
  //tree->Branch("CMS_zz4l_p2bplusKD",&CMS_zz4l_p2bplusKD,"CMS_zz4l_p2bplusKD/D");
  //tree->Branch("CMS_zz4l_p2hplusKD",&CMS_zz4l_p2hplusKD,"CMS_zz4l_p2hplusKD/D");
  //tree->Branch("CMS_zz4l_p2hminusKD",&CMS_zz4l_p2hminusKD,"CMS_zz4l_p2hminusKD/D");




  for (int i = 0; i < oldtree->GetEntries(); i++)
    {

      oldtree->GetEntry(i);

      bool notDuplicateEvent = true;

      ULong64_t runId = Run;
      ULong64_t eventId = Event;
      ULong64_t lumiId = LumiSect;


      for (unsigned int n = 0; n < runVec.size(); n++)
        {
          if(runId == runVec[n] && lumiId == lumiVec[n] && eventId == eventVec[n]){notDuplicateEvent = false;}
        }

      if(!notDuplicateEvent) continue;

      //if(!passedQCDSelection) continue;

      if( mz1 < globalMZ1Low ) continue;
      if( mz2 < globalMZ2Low ) continue;
      if( mzz < 86 ) continue;
      if( mzz > 96 ) continue;
      if( finalState == 3)
	{
	  if (_fs != 3 && _fs != 4) continue;
	}
      else{
	if( _fs != finalState) continue;
      }
      if(Z4l_minMass2l < 4) continue;

      if( abs(idL1) == 13 && pTL1 < globalPtMuCut ) continue;
      if( abs(idL2) == 13 && pTL2 < globalPtMuCut ) continue;
      if( abs(idL3) == 13 && pTL3 < globalPtMuCut ) continue;
      if( abs(idL4) == 13 && pTL4 < globalPtMuCut ) continue;

      if( abs(idL1) == 11 && pTL1 < globalPtElCut ) continue;
      if( abs(idL2) == 11 && pTL2 < globalPtElCut ) continue;
      if( abs(idL3) == 11 && pTL3 < globalPtElCut ) continue;
      if( abs(idL4) == 11 && pTL4 < globalPtElCut ) continue;

      mass = mzz;
      //if(mekdLD == NULL) mekd == mekdLDNew;
      //else mekd = mekdLD;
      mekdLog = TMath::Log(mekdSig) - TMath::Log(mekdBg);
      mekd = mekdLD;
      //mekd = mekdLDNew;
      //mela = melaLD;
      massErr = (double)m4lErr/mass;
      /*
      CMS_zz4l_pseudoKD = pseudoKD;
      CMS_zz4l_graviKD = graviKD;
      CMS_zz4l_p0hplusKD = p0hplusKD;
      CMS_zz4l_p1plusKD = p1plusKD;
      CMS_zz4l_p1minusKD = p1minusKD;
      CMS_zz4l_qqgraviKD = qqgraviKD;
      CMS_zz4l_PIp1plusKD = PIp1plusKD;
      CMS_zz4l_PIp1minusKD = PIp1minusKD;
      CMS_zz4l_PIp2plusKD = PIp2plusKD;
      CMS_zz4l_p2bplusKD = p2bplusKD;
      CMS_zz4l_p2hplusKD = p2hplusKD;
      CMS_zz4l_p2hminusKD = p2hminusKD;
      */

      mela = probRatioMEKD(3e-06, me_H, me_qqZZtchan);
      CMS_zz4l_pseudoKD = probRatioMEKD(3e-08, me_H, me_qqZZschan);
      if (CMS_zz4l_pseudoKD != CMS_zz4l_pseudoKD){ CMS_zz4l_pseudoKD = probRatioMEKD(3e-08, me_H, me_qqZZschan2);}
      CMS_zz4l_graviKD = -1;
      CMS_zz4l_p0hplusKD = -1;
      CMS_zz4l_p1plusKD = -1;
      CMS_zz4l_p1minusKD = -1;
      CMS_zz4l_qqgraviKD = -1;
      CMS_zz4l_PIp1plusKD = -1;
      CMS_zz4l_PIp1minusKD = -1;
      CMS_zz4l_PIp2plusKD = -1;
      CMS_zz4l_p2bplusKD = -1;
      CMS_zz4l_p2hplusKD = -1;
      CMS_zz4l_p2hminusKD = -1;

      //cout << " m" << fsn << "="<< mzz << " mZ1=" << mz1 << " mZ2=" << mz2 << " m4lErr=" << m4lErr << " mela=" << melaLD << " mekd=" << mekd <<"," 
      //<< mekdLog << endl
      //<< "     0-=" << pseudoKD << " gg2+m=" << graviKD << " 0h+=" << p0hplusKD << " qq1+=" << p1plusKD << " qq1-=" << p1minusKD << endl
      //<< "     qq2+=" << qqgraviKD << " PI1-=" << PIp1minusKD << " PI1+=" << PIp1plusKD << " PI2+=" << PIp2plusKD << " 2b+=" << p2bplusKD 
      //<< " 2h+=" << p2hplusKD << " 2h-=" << p2hminusKD << endl;

      cout << Event <<  " m" << fsn << "="<< mzz << " mZ1=" << mz1 << " mZ2=" << mz2 << " Xaxis=" << mela << " Yaxis=" << CMS_zz4l_pseudoKD << endl;

      cout << mass4lNoFSR << "  " << FSRZ1 << "  " << FSRZ2 << "  " << me_H << "   " << me_qqZZtchan << "  " << me_qqZZschan << "  " << me_qqZZschan2 << endl;
      cout << idL1 << "  " <<  EL1 << "  " <<  pXL1 << "  " << pYL1 << "  " << pZL1 << endl;
      cout << idL2 << "  " <<  EL2 << "  " <<  pXL2 << "  " << pYL2 << "  " << pZL2 << endl;
      cout << idL3 << "  " <<  EL3 << "  " <<  pXL3 << "  " << pYL3 << "  " << pZL3 << endl;
      cout << idL4 << "  " <<  EL4 << "  " <<  pXL4 << "  " << pYL4 << "  " << pZL4 << endl;


      tree->Fill();

      runVec.push_back(runId);
      lumiVec.push_back(lumiId);
      eventVec.push_back(eventId);

    }

  hzzF->cd();
  tree->Write();
  hzzF->Close();



}

double probRatioMEKD(double c, double me2processA, double me2processB){
  if (me2processA + c * me2processB == 0) return -999.;
  return me2processA/( me2processA + c * me2processB );
}
