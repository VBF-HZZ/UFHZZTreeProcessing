#include <Riostream.h>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "Math/DistFunc.h"

//#include "printData.C"

int scanTrees(TString treeFile="", TString treeName = "")
{

  //if(treeFile == "") treeFile = "sync_SMHiggs_JHU";
  //if(treeFile == "") treeFile = "test";
  if(treeName == "") treeName = "AnaAfterHlt/passedEvents";
  if(treeName == "") treeName = "passedEvents_dataMC";

  TFile *f = new TFile(treeFile+".root","READ");
  TTree *t = (TTree*)f->Get(treeName);

  //TCut scanCuts = "passedFullSelection && passedQCDcut && mass4l > 70";
  TCut scanCuts = "passedZ4lSelection && mass4l > 86 && mass4l < 96";
  //TCut scanCuts = "passedFullSelection && passedQCDcut && (VBFJet1 && VBFJet2) && mass4l > 100";
  
  TString treeFileTxt = treeFile+".txt";
  ((TTreePlayer*)(t->GetPlayer()))->SetScanRedirect(true);
  ((TTreePlayer*)(t->GetPlayer()))->SetScanFileName(treeFileTxt);
  t->Scan("Event:mass4l:massZ1:massZ2:massErrorUFCorr:FSRPhot1_Pt:FSRPhot2_Pt:JHUKD_H_qqZZ_noPDF:dataMC_weight",scanCuts);
  //t->Scan("Event:pTL1:pTL2:pTL3:pTL4:FSRPhot1_eta:FSRPhot2_eta:JHUKD_H_qqZZ_noPDF:pTVBFJet1:etaVBFJet1:pTVBFJet2:etaVBFJet2:VBFDiJetMass:VBFDeltaEta:FisherDiscrim",scanCuts);
  //t->Scan("Event:mass4l:massZ1:massZ2:idL1:pTL1:etaL1:idL2:pTL2:etaL2:idL3:pTL3:etaL3:idL4:pTL4:etaL4:massErrorUFCorr:melaLD:FSRPhot1_Pt:FSRPhot1_eta:FSRPhot2_Pt:FSRPhot2_eta:pTVBFJet1:etaVBFJet1:pTVBFJet2:etaVBFJet2:VBFDiJetMass:FisherDiscrim",scanCuts);

  treeFileTxt = treeFile+"_MEMs.txt";
  ((TTreePlayer*)(t->GetPlayer()))->SetScanRedirect(true);
  ((TTreePlayer*)(t->GetPlayer()))->SetScanFileName(treeFileTxt);
  t->Scan("Event:mass4l:MEKD_noPDF:MEKD_ME_H_noPDF:MEKD_ME_ZZ_noPDF:MEKD_h0M_ZZ_noPDF_noFSR:MEKD_ME_h0M_noPDF_noFSR:",scanCuts);


  TTree *t1 = (TTree*)f->Get("AnaAfterHlt/jets/jetDumpTree");
  treeFileTxt = treeFile+"_jets.txt";
  ((TTreePlayer*)(t1->GetPlayer()))->SetScanRedirect(true);
  ((TTreePlayer*)(t1->GetPlayer()))->SetScanFileName(treeFileTxt);
  t1->Scan("Event:pT:eta:phi:PFid:PUIdPass:PUIdFlag:PUMVA");
 
  TTree *t2 = (TTree*)f->Get("AnaAfterHlt/muons/muonDumpTree");
  treeFileTxt = treeFile+"_looseMuons.txt";
  ((TTreePlayer*)(t2->GetPlayer()))->SetScanRedirect(true);
  ((TTreePlayer*)(t2->GetPlayer()))->SetScanFileName(treeFileTxt);
  t2->Scan("Event:pT:eta:phi:isPFMuon:SIP:relIso:isoNH:isoCH:isoPhot:dB");

  TTree *t3 = (TTree*)f->Get("AnaAfterHlt/electrons/electronDumpTree");
  treeFileTxt = treeFile+"_looseElectrons.txt";
  ((TTreePlayer*)(t3->GetPlayer()))->SetScanRedirect(true);
  ((TTreePlayer*)(t3->GetPlayer()))->SetScanFileName(treeFileTxt);
  t3->Scan("Event:pT:eta:phi:mvaID:SIP:relIso:isoNH:isoCH:isoPhot:rho");

  TTree *t4 = (TTree*)f->Get("AnaAfterHlt/photons/photonDumpTree");
  treeFileTxt = treeFile+"_fsrPhotons.txt";
  ((TTreePlayer*)(t4->GetPlayer()))->SetScanRedirect(true);
  ((TTreePlayer*)(t4->GetPlayer()))->SetScanFileName(treeFileTxt);
  t4->Scan("Event:pT:eta:phi:deltaR:relIso:isoNH:isoCH:isoCHPU:isoPhot");

  TTree *t5 = (TTree*)f->Get("AnaAfterHlt/finalLeptons/finalLepDumpTree");
  treeFileTxt = treeFile+"_finalLeptons.txt";
  ((TTreePlayer*)(t5->GetPlayer()))->SetScanRedirect(true);
  ((TTreePlayer*)(t5->GetPlayer()))->SetScanFileName(treeFileTxt);
  t5->Scan("Event:pdgid:pT:eta:phi:mvaID:dz:dxy:relIso:isoNH:isoCH:isoPhot:rho:dB");

  TString treeFileTxt = treeFile+"_MassErrs.txt";
  ((TTreePlayer*)(t->GetPlayer()))->SetScanRedirect(true);
  ((TTreePlayer*)(t->GetPlayer()))->SetScanFileName(treeFileTxt);
  t->Scan("Run:Event:mass4l:massZ1:massZ2:idL1:idL3:massErrorUFCorr:massErrorUF:massErrorUCSD",scanCuts);


  f->Close();

  treeFileTxt = treeFile+"_syncShort.txt";
  std::string tmp = treeFileTxt;
  printDataSyncShort(treeFile,treeName,tmp,40,12,5,7,70,1000);
  treeFileTxt = treeFile+"_syncLong.txt";
  tmp = treeFileTxt;
  printDataSyncLong(treeFile,treeName,tmp,40,12,5,7,70,1000);

  return 0;
 
}





void printData(TString file, TString treeName, std::string outFile,double massZ1Cut, double massZ2Cut, double ptMuCut, double ptElCut, double m4lLowCut, double m4lHighCut)
{

  using namespace std;

  ULong64_t Run, Event, LumiSect;
  double mass4l, mass4lNoFSR, pT4l, massZ1, massZ2;
  double eta4l, phi4l;
  int idL1, idL2, idL3, idL4;
  double IPL1, IPL2, IPL3, IPL4;
  double dIPL1, dIPL2, dIPL3, dIPL4;
  double EL1, EL2, EL3, EL4;
  double EZ1, EZ2;
  double pTVBFJet1, pTVBFJet2, VBFDiJetMass, VBFDeltaEta, FisherDiscrim; 
  double pTL1, pTL2, pTL3, pTL4;
  double pXL1, pXL2, pXL3, pXL4;
  double pYL1, pYL2, pYL3, pYL4;
  double pZL1, pZL2, pZL3, pZL4;
  double pTZ1, pTZ2;
  double pXZ1, pXZ2;
  double pYZ1, pYZ2;
  double pZZ1, pZZ2;
  double etaL1, etaL2, etaL3, etaL4;
  double phiL1, phiL2, phiL3, phiL4;
  double isoCHL1, isoCHL2, isoCHL3, isoCHL4;
  double isoNHL1, isoNHL2, isoNHL3, isoNHL4;
  double isoPhotL1, isoPhotL2, isoPhotL3, isoPhotL4;
  double RelIsoL1, RelIsoL2, RelIsoL3, RelIsoL4;
  double RelIsoUCL1, RelIsoUCL2, RelIsoUCL3, RelIsoUCL4;
  double muRhoCor, elRhoCor;
  double worstIso, worstSIP;
  double cosTheta1, cosTheta2, Phi, cosThetaStar, phiStar1, phiStar2, Phi1, Phi2;
  int nJets, nVtx, nPhotons;
  double metVal, rawRho;
  int extraLep_id[10];
  double extraLep_pT[10], extraLep_iso[10], extraLep_e[10];
  double extraLep_pX[10], extraLep_pY[10], extraLep_pZ[10];
  double extraLep_eta[10], extraLep_phi[10], extraLep_sip[10];
  double extraLep_chIso[10], extraLep_nhIso[10], extraLep_phIso[10];
  double massError;
  double minM3l, minDeltR, maxP;
  double melaLD;
  double FSRPhot1_Pt,FSRPhot2_Pt,FSRPhot1_eta,FSRPhot2_eta,FSRPhot1_phi,FSRPhot2_phi;
  bool FSR_Z1,FSR_Z2;

  TBranch *b_extraLep_id;
  TBranch *b_extraLep_pT, *b_extraLep_iso, *b_extraLep_e;
  TBranch *b_extraLep_pX, *b_extraLep_pY, *b_extraLep_pZ;
  TBranch *b_extraLep_eta, *b_extraLep_phi, *b_extraLep_sip;
  TBranch *b_extraLep_nhIso, *b_extraLep_chIso, *b_extraLep_phIso;

  double thetaPhoton_deg, theta12_deg, theta13_deg, theta14_deg, minMass2Lep, maxMass2Lep;
  double thetaPhotonZ_deg;

  TFile *f = new TFile(file+".root","READ");
  TTree *tree = (TTree*)f->Get(treeName);
  tree->SetBranchAddress("Run",&Run);
  tree->SetBranchAddress("Event",&Event);
  tree->SetBranchAddress("LumiSect",&LumiSect);
  tree->SetBranchAddress("mass4l",&mass4l);
  tree->SetBranchAddress("mass4lNoFSR",&mass4lNoFSR);
  tree->SetBranchAddress("minMass3l",&minM3l);
  tree->SetBranchAddress("minDeltaR",&minDeltR);
  tree->SetBranchAddress("Z4l_maxP",&maxP);
  tree->SetBranchAddress("pT4l",&pT4l);
  tree->SetBranchAddress("Z4l_thetaPhoton_deg",&thetaPhoton_deg);
  tree->SetBranchAddress("Z4l_thetaPhotonZ_deg",&thetaPhotonZ_deg);
  tree->SetBranchAddress("Z4l_theta12_deg",&theta12_deg);
  tree->SetBranchAddress("Z4l_theta13_deg",&theta13_deg);
  tree->SetBranchAddress("Z4l_theta14_deg",&theta14_deg);
  tree->SetBranchAddress("Z4l_minMass2l",&minMass2Lep);
  tree->SetBranchAddress("Z4l_maxMass2l",&maxMass2Lep);
  tree->SetBranchAddress("massZ1",&massZ1);
  tree->SetBranchAddress("massZ2",&massZ2);
  tree->SetBranchAddress("eta4l",&eta4l);
  tree->SetBranchAddress("phi4l",&phi4l);
  tree->SetBranchAddress("EL1",&EL1);
  tree->SetBranchAddress("EL2",&EL2);
  tree->SetBranchAddress("EL3",&EL3);
  tree->SetBranchAddress("EL4",&EL4);
  tree->SetBranchAddress("EZ1",&EZ1);
  tree->SetBranchAddress("EZ2",&EZ2);
  tree->SetBranchAddress("idL1",&idL1);
  tree->SetBranchAddress("idL2",&idL2);
  tree->SetBranchAddress("idL3",&idL3);
  tree->SetBranchAddress("idL4",&idL4);
  tree->SetBranchAddress("pTL1",&pTL1);
  tree->SetBranchAddress("pTL2",&pTL2);
  tree->SetBranchAddress("pTL3",&pTL3);
  tree->SetBranchAddress("pTL4",&pTL4);
  tree->SetBranchAddress("pXL1",&pXL1);
  tree->SetBranchAddress("pXL2",&pXL2);
  tree->SetBranchAddress("pXL3",&pXL3);
  tree->SetBranchAddress("pXL4",&pXL4);
  tree->SetBranchAddress("pYL1",&pYL1);
  tree->SetBranchAddress("pYL2",&pYL2);
  tree->SetBranchAddress("pYL3",&pYL3);
  tree->SetBranchAddress("pYL4",&pYL4);
  tree->SetBranchAddress("pZL1",&pZL1);
  tree->SetBranchAddress("pZL2",&pZL2);
  tree->SetBranchAddress("pZL3",&pZL3);
  tree->SetBranchAddress("pZL4",&pZL4);
  tree->SetBranchAddress("pTZ1",&pTZ1);
  tree->SetBranchAddress("pTZ2",&pTZ2);
  tree->SetBranchAddress("pXZ1",&pXZ1);
  tree->SetBranchAddress("pXZ2",&pXZ2);
  tree->SetBranchAddress("pYZ1",&pYZ1);
  tree->SetBranchAddress("pYZ2",&pYZ2);
  tree->SetBranchAddress("pZZ1",&pZZ1);
  tree->SetBranchAddress("pZZ2",&pZZ2);
  tree->SetBranchAddress("etaL1",&etaL1);
  tree->SetBranchAddress("etaL2",&etaL2);
  tree->SetBranchAddress("etaL3",&etaL3);
  tree->SetBranchAddress("etaL4",&etaL4);
  tree->SetBranchAddress("phiL1",&phiL1);
  tree->SetBranchAddress("phiL2",&phiL2);
  tree->SetBranchAddress("phiL3",&phiL3);
  tree->SetBranchAddress("phiL4",&phiL4);
  tree->SetBranchAddress("IPL1",&IPL1);
  tree->SetBranchAddress("IPL2",&IPL2);
  tree->SetBranchAddress("IPL3",&IPL3);
  tree->SetBranchAddress("IPL4",&IPL4);
  tree->SetBranchAddress("dIPL1",&dIPL1);
  tree->SetBranchAddress("dIPL2",&dIPL2);
  tree->SetBranchAddress("dIPL3",&dIPL3);
  tree->SetBranchAddress("dIPL4",&dIPL4);
  tree->SetBranchAddress("isoNHL1",&isoNHL1);
  tree->SetBranchAddress("isoNHL2",&isoNHL2);
  tree->SetBranchAddress("isoNHL3",&isoNHL3);
  tree->SetBranchAddress("isoNHL4",&isoNHL4);
  tree->SetBranchAddress("isoCHL1",&isoCHL1);
  tree->SetBranchAddress("isoCHL2",&isoCHL2);
  tree->SetBranchAddress("isoCHL3",&isoCHL3);
  tree->SetBranchAddress("isoCHL4",&isoCHL4);
  tree->SetBranchAddress("isoPhotL1",&isoPhotL1);
  tree->SetBranchAddress("isoPhotL2",&isoPhotL2);
  tree->SetBranchAddress("isoPhotL3",&isoPhotL3);
  tree->SetBranchAddress("isoPhotL4",&isoPhotL4);
  tree->SetBranchAddress("RelIsoL1",&RelIsoL1);
  tree->SetBranchAddress("RelIsoL2",&RelIsoL2);
  tree->SetBranchAddress("RelIsoL3",&RelIsoL3);
  tree->SetBranchAddress("RelIsoL4",&RelIsoL4);
  tree->SetBranchAddress("RelIsoUCL1",&RelIsoUCL1);
  tree->SetBranchAddress("RelIsoUCL2",&RelIsoUCL2);
  tree->SetBranchAddress("RelIsoUCL3",&RelIsoUCL3);
  tree->SetBranchAddress("RelIsoUCL4",&RelIsoUCL4);
  tree->SetBranchAddress("muRhoCor",&muRhoCor);
  tree->SetBranchAddress("elRhoCor",&elRhoCor);
  tree->SetBranchAddress("worstIso",&worstIso);
  tree->SetBranchAddress("worstSIP",&worstSIP);
  tree->SetBranchAddress("cosTheta1",&cosTheta1);
  tree->SetBranchAddress("cosTheta2",&cosTheta2);
  tree->SetBranchAddress("Phi",&Phi);
  tree->SetBranchAddress("cosThetaStar",&cosThetaStar);
  tree->SetBranchAddress("phiStar1",&phiStar1);
  tree->SetBranchAddress("phiStar2",&phiStar2);
  tree->SetBranchAddress("Phi1",&Phi1);
  tree->SetBranchAddress("Phi2",&Phi2);
  tree->SetBranchAddress("nJets",&nJets);
  tree->SetBranchAddress("nVtx",&nVtx);
  tree->SetBranchAddress("nPhotons",&nPhotons);
  tree->SetBranchAddress("metVal",&metVal);
  tree->SetBranchAddress("rawRho",&rawRho);
  tree->SetBranchAddress("extraLep_id", extraLep_id, &b_extraLep_id);
  tree->SetBranchAddress("extraLep_pT", extraLep_pT, &b_extraLep_pT);
  tree->SetBranchAddress("extraLep_pX", extraLep_pX, &b_extraLep_pX);
  tree->SetBranchAddress("extraLep_pY", extraLep_pY, &b_extraLep_pY);
  tree->SetBranchAddress("extraLep_pZ", extraLep_pZ, &b_extraLep_pZ);
  tree->SetBranchAddress("extraLep_e",extraLep_e, &b_extraLep_e);
  tree->SetBranchAddress("extraLep_iso", extraLep_iso, &b_extraLep_iso);
  tree->SetBranchAddress("extraLep_chIso",extraLep_chIso, &b_extraLep_chIso);
  tree->SetBranchAddress("extraLep_nhIso", extraLep_nhIso, &b_extraLep_nhIso);
  tree->SetBranchAddress("extraLep_phIso",extraLep_phIso, &b_extraLep_phIso);
  tree->SetBranchAddress("extraLep_eta",extraLep_eta, &b_extraLep_eta);
  tree->SetBranchAddress("extraLep_phi",extraLep_phi, &b_extraLep_phi);
  tree->SetBranchAddress("extraLep_sip",extraLep_sip, &b_extraLep_sip);
  tree->SetBranchAddress("massErrorUFCorr",&massError);
  tree->SetBranchAddress("melaLD",&melaLD);
  tree->SetBranchAddress("FSRPhot1_Pt",&FSRPhot1_Pt);
  tree->SetBranchAddress("FSRPhot2_Pt",&FSRPhot2_Pt);
  tree->SetBranchAddress("FSRPhot1_eta",&FSRPhot1_eta);
  tree->SetBranchAddress("FSRPhot2_eta",&FSRPhot2_eta);
  tree->SetBranchAddress("FSRPhot1_phi",&FSRPhot1_phi);
  tree->SetBranchAddress("FSRPhot2_phi",&FSRPhot2_phi);
  tree->SetBranchAddress("FSR_Z1",&FSR_Z1);
  tree->SetBranchAddress("FSR_Z2",&FSR_Z2);

  tree->SetBranchAddress("pTVBFJet1",&pTVBFJet1);
  tree->SetBranchAddress("pTVBFJet2",&pTVBFJet2);
  tree->SetBranchAddress("VBFDiJetMass",&VBFDiJetMass);
  tree->SetBranchAddress("VBFDeltaEta",&VBFDeltaEta);
  tree->SetBranchAddress("FisherDiscrim",&FisherDiscrim);


  std::vector<ULong64_t> runVec, lumiVec, eventVec;


  ofstream out;
  out.open(outFile.c_str());

  int counter = 0, counter_4mu = 0, counter_4e = 0, counter_2e2mu = 0, VBFCounter = 0;

  tree->BuildIndex("Run");
  TTreeIndex *index = (TTreeIndex*)tree->GetTreeIndex();

  for( int i = 0; i < index->GetN(); i++ )
    {

      ULong64_t local = tree->LoadTree( index->GetIndex()[i] );
      tree->GetEntry(local);

      bool notDuplicateEvent = true;

      //ULong64_t runId = Run;
      //ULong64_t eventId = Event;
      //ULong64_t lumiId = LumiSect;

      /*
      for (unsigned int n = 0; n < runVec.size(); n++)
	{
	  if(runId == runVec[n] && lumiId == lumiVec[n] && eventId == eventVec[n]){notDuplicateEvent = false;}
	}
      */

      if (notDuplicateEvent)
	{

	  if(massZ1 < massZ1Cut) continue;
	  if(massZ2 < massZ2Cut) continue;
	  if(mass4l < m4lLowCut) continue;
	  if(mass4l > m4lHighCut) continue;

	  if( abs(idL1) == 13 && pTL1 < ptMuCut ) continue;
	  if( abs(idL2) == 13 && pTL2 < ptMuCut ) continue;
	  if( abs(idL3) == 13 && pTL3 < ptMuCut ) continue;
	  if( abs(idL4) == 13 && pTL4 < ptMuCut ) continue;

	  if( abs(idL1) == 11 && pTL1 < ptElCut ) continue;
	  if( abs(idL2) == 11 && pTL2 < ptElCut ) continue;
	  if( abs(idL3) == 11 && pTL3 < ptElCut ) continue;
	  if( abs(idL4) == 11 && pTL4 < ptElCut ) continue;
	  

	  //runVec.push_back(runId);
	  //lumiVec.push_back(lumiId);
	  //eventVec.push_back(eventId);
	  
	  counter++;
	  if( abs(idL1) == 13 && abs(idL2) == 13 && abs(idL3) == 13 && abs(idL4) == 13 ){counter_4mu++;}
	  if( abs(idL1) == 11 && abs(idL2) == 11 && abs(idL3) == 11 && abs(idL4) == 11 ){counter_4e++;}
	  if( abs(idL1) == 13 && abs(idL2) == 13 && abs(idL3) == 11 && abs(idL4) == 11 ){counter_2e2mu++;}
	  if( abs(idL1) == 11 && abs(idL2) == 11 && abs(idL3) == 13 && abs(idL4) == 13 ){counter_2e2mu++;}


		  

	  out << "########" << endl
	      << counter << endl
	      << "########" << endl;

	      out.precision(7);
	      out << left
		  << "Run:  " << Run          << "    Event: " << Event    << "   LumiSec: " << LumiSect << endl
		  << "m4l:  " << mass4l   << "   m4lError: " << massError   << "      mZ1: " << massZ1   << "       mZ2: " << massZ2  << endl
		  << "FSR_Z1  " << FSR_Z1 << " pt: " << FSRPhot1_Pt << "    FSR_Z2: " << FSR_Z2 << " pt: " << FSRPhot2_Pt << endl  
		  << "pTVBFJet1: " << pTVBFJet1 << "  pTVBFJet2: " << pTVBFJet2 << "  jjmass/deta: " << VBFDiJetMass <<"/"<<VBFDeltaEta
		  << "  Fisher: " << FisherDiscrim << endl 
		  << "cos(theta1): " << cosTheta1 << "  cos(theta2): " << cosTheta2 << "  cos(theta*): " << cosThetaStar << endl
		  << "phi: " << Phi << "  Phi1: " << Phi1 << "  Phi2: " << Phi2 << "  Phi*1: " << phiStar1 << "  Phi*2: " << phiStar2 << endl
		  << " MuonRhoCorr: " << muRhoCor << "      ElecRhoCor: " << elRhoCor << " Rho: " << rawRho << endl
		  << "MELALD: " << melaLD << "   m4lNoFSR: " << mass4lNoFSR << endl
		  << "minM3l:  " << minM3l << "   maxP:  " << maxP << "  minDeltaR:  " << minDeltR << endl
		  << "thetaPhotonZ:  " << thetaPhotonZ_deg << "  thetaPhoton:  " << thetaPhoton_deg << endl 
		  << "maxMass2L:  " << maxMass2Lep << "  minMass2L:  " << minMass2Lep << endl
		  << "pT4l: " << pT4l         << "    eta4l: " << eta4l    << "     phi4l: " << phi4l   << endl
		  << "worstIso: " << worstIso << "    worstSIP: " << worstSIP << endl
		  << "nJets: " << nJets << "    nPhotons: " << nPhotons << "   nVtx: " << nVtx << endl
		  << "MET: " << metVal << endl
		  << "L1: " << idL1         << "       pT: " << pTL1 << "  eta: " << etaL1 << "  phi: " << phiL1
		  << "  iso: CH " << isoCHL1 << " NH " << isoNHL1 << " Ph " << isoPhotL1 << "  RelIso: " << RelIsoL1 << "," << RelIsoUCL1 << endl
		  << "L2: " << idL2         << "       pT: " << pTL2 << "  eta: " << etaL2 << "  phi: " << phiL2
		  << "  iso: CH " << isoCHL2 << " NH " << isoNHL2 << " Ph " << isoPhotL2 << "  RelIso: " << RelIsoL2 << "," << RelIsoUCL2 << endl
		  << "L3: " << idL3         << "       pT: " << pTL3 << "  eta: " << etaL3 << "  phi: " << phiL3
		  << "  iso: CH " << isoCHL3 << " NH " << isoNHL3 << " Ph " << isoPhotL3 << "  RelIso: " << RelIsoL3 << "," << RelIsoUCL3 << endl
		  << "L4: " << idL4         << "       pT: " << pTL4 << "  eta: " << etaL4 << "  phi: " << phiL4
		  << "  iso: CH " << isoCHL4 << " NH " << isoNHL4 << " Ph " << isoPhotL4 << "  RelIso: " << RelIsoL4 << "," << RelIsoUCL4 << endl
		  << "  ID             px             py              pz             e        IP       deltaIP  "  << endl
		  << "-------------------------------------------------------------------------------------------" << endl
		  << right
		  << "  Z1    " << setprecision(4) << setw(12) << pXZ1 << "   " << setprecision(4) << setw(12) << pYZ1 << "   " 
		  << setprecision(4) << setw(12) << pZZ1 << "  " << setprecision(4) << setw(12) <<  EZ1 << endl
		  << "  Z2    " << setprecision(4) << setw(12) << pXZ2 << "   " << setprecision(4) << setw(12) <<  pYZ2 << "   " 
		  << setprecision(4) << setw(12) << pZZ2 << "  " << setprecision(4) << setw(12) << EZ2 << endl
		  << "  L1    " << setprecision(4) << setw(12) << pXL1 << "   " << setprecision(4) << setw(12) << pYL1 << "   " 
		  << setprecision(4) << setw(12) << pZL1 << "  " << setprecision(4) << setw(12) << EL1 
		  << "   " << setprecision(4) << setw(12) << IPL1 << "  " << setprecision(4) << setw(12) << dIPL1 << endl
		  << "  L2    " << setprecision(4) << setw(12) << pXL2 << "   " << setprecision(4) << setw(12) << pYL2 << "   " 
		  << setprecision(4) << setw(12) << pZL2 << "  " << setprecision(4) << setw(12) << EL2 
                  << "   " << setprecision(4) << setw(12) << IPL2 << "  " << setprecision(4) << setw(12) << dIPL2 << endl
		  << "  L3    " << setprecision(4) << setw(12) << pXL3 << "   " << setprecision(4) << setw(12) << pYL3 << "   " 
		  << setprecision(4) << setw(12) << pZL3 << "  " << setprecision(4) << setw(12) << EL3 
                  << "   " << setprecision(4) << setw(12) << IPL3 << "  " << setprecision(4) << setw(12) << dIPL3 << endl
		  << "  L4    " << setprecision(4) << setw(12) << pXL4 << "   " << setprecision(4) << setw(12) << pYL4 << "   " 
		  << setprecision(4) << setw(12) << pZL4 << "  " << setprecision(4) << setw(12) << EL4 
                  << "   " << setprecision(4) << setw(12) << IPL4 << "  " << setprecision(4) << setw(12) << dIPL4 << endl;
	      

	      out << "  Other Leptons" << endl
		  << "---------------" << endl;

	      for( int k = 0; k < 10; k++ )
		{
		  if( abs(extraLep_id[k]) == 11 || abs(extraLep_id[k]) == 13 )
		    {
		      if( extraLep_pT[k] != pTL1 && extraLep_pT[k] != pTL2 && extraLep_pT[k] != pTL3 && extraLep_pT[k] != pTL4 )
			{
		      
		      out << " " << extraLep_id[k] << setprecision(4) << setw(12) << extraLep_pX[k] << "   " << setprecision(4) << setw(12) << extraLep_pY[k] << "   "
			  << setprecision(4) << setw(12) << extraLep_pZ[k] << "  " << setprecision(4) << setw(12) <<  extraLep_e[k] << "  "  << endl
			  << "      pT: " << extraLep_pT[k] << "   eta: " << extraLep_eta[k] << "   phi: " << extraLep_phi[k]  
			  << "  SIP: " << extraLep_sip[k] << "  iso: " << extraLep_iso[k]   
			  << "  CH " << extraLep_chIso[k] << " NH " << extraLep_nhIso[k] << " h " << extraLep_phIso[k] << endl;
			}
		    }
		  
		}

	      
	      out << endl << endl;
	      
	      if(pTVBFJet1 > 0 && pTVBFJet2 > 0 && VBFDiJetMass > 300 && VBFDeltaEta > 3) VBFCounter++;
	      
	}
    
    }

  out << counter << " event candidates found" << endl;
  out << counter_4mu << " 4mu events " << endl;
  out << counter_4e << " 4e events " << endl;
  out << counter_2e2mu << " 2e2mu events" << endl;
  out << VBFCounter << " VBF tagged events" << endl;

  //runVec.clear();
  //lumiVec.clear();
  //eventVec.clear();  
  
  out.close();
}


void printDataShort(TString file, TString treeName, std::string outFile,double massZ1Cut, double massZ2Cut, double ptMuCut, double ptElCut, double m4lLowCut, double m4lHighCut)
{

  using namespace std;

  std::vector<ULong64_t> runVec, lumiVec, eventVec;

  ULong64_t Run, Event, LumiSect;
  double mass4l,  massZ1, massZ2;
  int idL1, idL2, idL3, idL4;
  double pTL1, pTL2, pTL3, pTL4;
  double melaLD;
  string type;

  TFile *f = new TFile(file+".root","READ");
  TTree *tree = (TTree*)f->Get(treeName);
  tree->SetBranchAddress("Run",&Run);
  tree->SetBranchAddress("Event",&Event);
  tree->SetBranchAddress("LumiSect",&LumiSect);
  tree->SetBranchAddress("mass4l",&mass4l);
  tree->SetBranchAddress("massZ1",&massZ1);
  tree->SetBranchAddress("massZ2",&massZ2);
  tree->SetBranchAddress("idL1",&idL1);
  tree->SetBranchAddress("idL2",&idL2);
  tree->SetBranchAddress("idL3",&idL3);
  tree->SetBranchAddress("idL4",&idL4);
  tree->SetBranchAddress("pTL1",&pTL1);
  tree->SetBranchAddress("pTL2",&pTL2);
  tree->SetBranchAddress("pTL3",&pTL3);
  tree->SetBranchAddress("pTL4",&pTL4);
  tree->SetBranchAddress("melaLD",&melaLD);


  ofstream out;
  out.open(outFile.c_str());

  int counter = 0, counter_4mu = 0, counter_4e = 0, counter_2e2mu = 0;

  tree->BuildIndex("Run");
  TTreeIndex *index = (TTreeIndex*)tree->GetTreeIndex();

  for( int i = 0; i < index->GetN(); i++ )
    {

      ULong64_t local = tree->LoadTree( index->GetIndex()[i] );
      tree->GetEntry(local);

      bool notDuplicateEvent = true;

      //ULong64_t runId = Run;
      //ULong64_t eventId = Event;
      //ULong64_t lumiId = LumiSect;
      /*
      for (unsigned int n = 0; n < runVec.size(); n++)
	{
	  if(runId == runVec[n] && lumiId == lumiVec[n] && eventId == eventVec[n]){notDuplicateEvent = false;}
	}
      */

      if (notDuplicateEvent)
	{

	  if(massZ1 < massZ1Cut) continue;
	  if(massZ2 < massZ2Cut) continue;
	  if(mass4l < m4lLowCut) continue;
	  if(mass4l > m4lHighCut) continue;

	  if( abs(idL1) == 13 && pTL1 < ptMuCut ) continue;
	  if( abs(idL2) == 13 && pTL2 < ptMuCut ) continue;
	  if( abs(idL3) == 13 && pTL3 < ptMuCut ) continue;
	  if( abs(idL4) == 13 && pTL4 < ptMuCut ) continue;

	  if( abs(idL1) == 11 && pTL1 < ptElCut ) continue;
	  if( abs(idL2) == 11 && pTL2 < ptElCut ) continue;
	  if( abs(idL3) == 11 && pTL3 < ptElCut ) continue;
	  if( abs(idL4) == 11 && pTL4 < ptElCut ) continue;

	  //runVec.push_back(runId);
	  //lumiVec.push_back(lumiId);
	  //eventVec.push_back(eventId);
	  
	  counter++;
	  if( abs(idL1) == 13 && abs(idL2) == 13 && abs(idL3) == 13 && abs(idL4) == 13 ){counter_4mu++;type="4mu";}
	  if( abs(idL1) == 11 && abs(idL2) == 11 && abs(idL3) == 11 && abs(idL4) == 11 ){counter_4e++;type="4e";}
	  if( abs(idL1) == 13 && abs(idL2) == 13 && abs(idL3) == 11 && abs(idL4) == 11 ){counter_2e2mu++;type="2mu2e";}
	  if( abs(idL1) == 11 && abs(idL2) == 11 && abs(idL3) == 13 && abs(idL4) == 13 ){counter_2e2mu++;type="2e2mu";}
		  
	      out.precision(7);
	      out << left
		  << setw(10) << Run  <<  "   " << setw(10) <<  Event  << "  " << setw(10) << LumiSect << "  " << setw(5) << type  
		  << "  " << setw(8) <<  mass4l   << "  " << setw(8) <<  massZ1  << "  " << setw(8) <<  massZ2  << "  " << melaLD << endl;
		
	}
    
    }

  out << endl << endl;
  out << counter << " event candidates found" << endl;
  out << counter_4mu << " 4mu events " << endl;
  out << counter_4e << " 4e events " << endl;
  out << counter_2e2mu << " 2e2mu events" << endl;

  //runVec.clear();
  //lumiVec.clear();
  //eventVec.clear();  
  
  out.close();
}


void printDataSyncShort(TString file, TString treeName, std::string outFile,double massZ1Cut, double massZ2Cut, double ptMuCut, double ptElCut, double m4lLowCut, double m4lHighCut)
{

  using namespace std;

  std::vector<ULong64_t> runVec, lumiVec, eventVec;

  ULong64_t Run, Event, LumiSect;
  double mass4l,  massZ1, massZ2,massErrRaw,massErrCorr;
  int idL1, idL2, idL3, idL4;
  double pTL1, pTL2, pTL3, pTL4;
  double melaLD, pT4l, FisherDiscrim;
  string type;
  bool VBFJet1, VBFJet2;
  double pTVBFJet1, pTVBFJet2;
  double JHUKD_H_qqZZ_noPDF,JHUKD_H_h0M_noPDF,JHUKD_H_h0P_noPDF;
  double JHUKD_H_h1P_noPDF,JHUKD_H_h1M_noPDF,JHUKD_H_ggh2P_noPDF;
  double JHUKD_H_qqh2P_noPDF;

  double VBFDiJetMass, VBFDeltaEta;

  int nJets = 0;

  TFile *f = new TFile(file+".root","READ");
  TTree *tree = (TTree*)f->Get(treeName);
  tree->SetBranchAddress("Run",&Run);
  tree->SetBranchAddress("Event",&Event);
  tree->SetBranchAddress("LumiSect",&LumiSect);
  tree->SetBranchAddress("mass4l",&mass4l);
  tree->SetBranchAddress("massZ1",&massZ1);
  tree->SetBranchAddress("massZ2",&massZ2);
  tree->SetBranchAddress("idL1",&idL1);
  tree->SetBranchAddress("idL2",&idL2);
  tree->SetBranchAddress("idL3",&idL3);
  tree->SetBranchAddress("idL4",&idL4);
  tree->SetBranchAddress("pTL1",&pTL1);
  tree->SetBranchAddress("pTL2",&pTL2);
  tree->SetBranchAddress("pTL3",&pTL3);
  tree->SetBranchAddress("pTL4",&pTL4);
  tree->SetBranchAddress("pT4l",&pT4l);
  tree->SetBranchAddress("melaLD",&melaLD);
  tree->SetBranchAddress("massErrorUF",&massErrRaw);
  tree->SetBranchAddress("massErrorUFCorr",&massErrCorr);

  tree->SetBranchAddress("VBFJet1",&VBFJet1);
  tree->SetBranchAddress("pTVBFJet1",&pTVBFJet1);
  tree->SetBranchAddress("VBFJet2",&VBFJet2);
  tree->SetBranchAddress("pTVBFJet2",&pTVBFJet2);
  tree->SetBranchAddress("VBFDiJetMass",&VBFDiJetMass);
  tree->SetBranchAddress("VBFDeltaEta",&VBFDeltaEta);
  tree->SetBranchAddress("FisherDiscrim",&FisherDiscrim);

  tree->SetBranchAddress("JHUKD_H_qqZZ_noPDF",&JHUKD_H_qqZZ_noPDF);
  tree->SetBranchAddress("JHUKD_H_h0M_noPDF",&JHUKD_H_h0M_noPDF);
  tree->SetBranchAddress("JHUKD_H_h0P_noPDF",&JHUKD_H_h0P_noPDF);
  tree->SetBranchAddress("JHUKD_H_h1P_noPDF",&JHUKD_H_h1P_noPDF);
  tree->SetBranchAddress("JHUKD_H_h1M_noPDF",&JHUKD_H_h1M_noPDF);
  tree->SetBranchAddress("JHUKD_H_ggh2P_noPDF",&JHUKD_H_ggh2P_noPDF);
  tree->SetBranchAddress("JHUKD_H_qqh2P_noPDF",&JHUKD_H_qqh2P_noPDF);
	  

  ofstream out;
  out.open(outFile.c_str());

  tree->BuildIndex("Run");
  TTreeIndex *index = (TTreeIndex*)tree->GetTreeIndex();

  for( int i = 0; i < index->GetN(); i++ )
    {

      ULong64_t local = tree->LoadTree( index->GetIndex()[i] );
      tree->GetEntry(local);

      bool notDuplicateEvent = true;

      
      //ULong64_t runId = Run;
      //ULong64_t eventId = Event;
      //ULong64_t lumiId = LumiSect;
      /*
      for (unsigned int n = 0; n < runVec.size(); n++)
	{
	  if(runId == runVec[n] && lumiId == lumiVec[n] && eventId == eventVec[n]){notDuplicateEvent = false;}
	}
    
      */
      if (notDuplicateEvent)
	{

	  if(massZ1 < massZ1Cut) continue;
	  if(massZ2 < massZ2Cut) continue;
	  if(mass4l < m4lLowCut) continue;
	  if(mass4l > m4lHighCut) continue;

	  //if(mass4l > 110 && mass4l < 140) continue;
	  //if(mass4l > 300) continue;
	  
	  if( abs(idL1) == 13 && pTL1 < ptMuCut ) continue;
	  if( abs(idL2) == 13 && pTL2 < ptMuCut ) continue;
	  if( abs(idL3) == 13 && pTL3 < ptMuCut ) continue;
	  if( abs(idL4) == 13 && pTL4 < ptMuCut ) continue;

	  if( abs(idL1) == 11 && pTL1 < ptElCut ) continue;
	  if( abs(idL2) == 11 && pTL2 < ptElCut ) continue;
	  if( abs(idL3) == 11 && pTL3 < ptElCut ) continue;
	  if( abs(idL4) == 11 && pTL4 < ptElCut ) continue;

	  //runVec.push_back(runId);
	  //lumiVec.push_back(lumiId);
	  //eventVec.push_back(eventId);

	  //out <<Run<<":"<<LumiSect<<":"<<Event<<":";
	  //out.setf(ios::fixed,ios::floatfield);
	  //out.precision(2);
	  //out << mass4l <<":"<< massZ1 <<":"<< massZ2 <<":"
	  //   << massErrRaw <<":"<< massErrCorr <<":";
	  //out.setf(ios::fixed,ios::floatfield);
	  //out.precision(3);
	  //out << melaLD <<endl;
	  nJets = 0;

	  if((VBFJet1 && !VBFJet2) || (VBFJet2 && !VBFJet1)) nJets = 1;
	  if(!VBFJet2 && !VBFJet1) nJets = 0;
	  if(VBFJet1 && VBFJet2) nJets = 2;

	  if(!VBFJet1){pTVBFJet1 = 0; VBFDiJetMass = 0; VBFDeltaEta = 0;}
	  if(!VBFJet2){pTVBFJet2 = 0; VBFDiJetMass = 0; VBFDeltaEta = 0;}

	  out <<Run<<":"<<LumiSect<<":"<<Event<<":";
	  out.setf(ios::fixed,ios::floatfield);
	  out.precision(2);
	  out << mass4l <<":"<< massZ1 <<":"<< massZ2 <<":"
	      << massErrRaw <<":"<< massErrCorr <<":";
	  out.setf(ios::fixed,ios::floatfield);
	  out.precision(3);
	  out << JHUKD_H_qqZZ_noPDF<<":";
	  out.precision(2);
	  out << pT4l << ":";
	  out.precision(0);
	  out << nJets << ":";
	  out.precision(2);
	  out << pTVBFJet1 << ":" << pTVBFJet2 << ":" << VBFDiJetMass << ":";
	  out.precision(3);
	  out << VBFDeltaEta << ":" << FisherDiscrim << endl;


	}
    
    }


  //runVec.clear();
  //lumiVec.clear();
  //eventVec.clear();  
  
  out.close();
}


void printDataSyncLong(TString file, TString treeName, std::string outFile,double massZ1Cut, double massZ2Cut, double ptMuCut, double ptElCut, double m4lLowCut, double m4lHighCut)
{

  using namespace std;

  std::vector<ULong64_t> runVec, lumiVec, eventVec;

  ULong64_t Run, Event, LumiSect;
  double mass4l,  massZ1, massZ2,massErrRaw,massErrCorr;
  int idL1, idL2, idL3, idL4;
  double pTL1, pTL2, pTL3, pTL4;
  double melaLD,pT4l;
  string type;
  bool VBFJet1, VBFJet2;
  double pTVBFJet1, pTVBFJet2;
  double JHUKD_H_qqZZ_noPDF,JHUKD_H_h0M_noPDF,JHUKD_H_h0P_noPDF;
  double JHUKD_H_h1P_noPDF,JHUKD_H_h1M_noPDF,JHUKD_H_ggh2P_noPDF;
  double JHUKD_H_qqh2P_noPDF;

  double VBFDiJetMass, VBFDeltaEta;

  int nJets = 0;

  TFile *f = new TFile(file+".root","READ");
  TTree *tree = (TTree*)f->Get(treeName);
  tree->SetBranchAddress("Run",&Run);
  tree->SetBranchAddress("Event",&Event);
  tree->SetBranchAddress("LumiSect",&LumiSect);
  tree->SetBranchAddress("mass4l",&mass4l);
  tree->SetBranchAddress("massZ1",&massZ1);
  tree->SetBranchAddress("massZ2",&massZ2);
  tree->SetBranchAddress("idL1",&idL1);
  tree->SetBranchAddress("idL2",&idL2);
  tree->SetBranchAddress("idL3",&idL3);
  tree->SetBranchAddress("idL4",&idL4);
  tree->SetBranchAddress("pTL1",&pTL1);
  tree->SetBranchAddress("pTL2",&pTL2);
  tree->SetBranchAddress("pTL3",&pTL3);
  tree->SetBranchAddress("pTL4",&pTL4);
  tree->SetBranchAddress("pT4l",&pT4l);
  tree->SetBranchAddress("melaLD",&melaLD);
  tree->SetBranchAddress("massErrorUF",&massErrRaw);
  tree->SetBranchAddress("massErrorUFCorr",&massErrCorr);

  tree->SetBranchAddress("VBFJet1",&VBFJet1);
  tree->SetBranchAddress("pTVBFJet1",&pTVBFJet1);
  tree->SetBranchAddress("VBFJet2",&VBFJet2);
  tree->SetBranchAddress("pTVBFJet2",&pTVBFJet2);
  tree->SetBranchAddress("VBFDiJetMass",&VBFDiJetMass);
  tree->SetBranchAddress("VBFDeltaEta",&VBFDeltaEta);

  tree->SetBranchAddress("JHUKD_H_qqZZ_noPDF",&JHUKD_H_qqZZ_noPDF);
  tree->SetBranchAddress("JHUKD_H_h0M_noPDF",&JHUKD_H_h0M_noPDF);
  tree->SetBranchAddress("JHUKD_H_h0P_noPDF",&JHUKD_H_h0P_noPDF);
  tree->SetBranchAddress("JHUKD_H_h1P_noPDF",&JHUKD_H_h1P_noPDF);
  tree->SetBranchAddress("JHUKD_H_h1M_noPDF",&JHUKD_H_h1M_noPDF);
  tree->SetBranchAddress("JHUKD_H_ggh2P_noPDF",&JHUKD_H_ggh2P_noPDF);
  tree->SetBranchAddress("JHUKD_H_qqh2P_noPDF",&JHUKD_H_qqh2P_noPDF);
	  

  ofstream out;
  out.open(outFile.c_str());

  tree->BuildIndex("Run");
  TTreeIndex *index = (TTreeIndex*)tree->GetTreeIndex();

  for( int i = 0; i < index->GetN(); i++ )
    {

      ULong64_t local = tree->LoadTree( index->GetIndex()[i] );
      tree->GetEntry(local);

      bool notDuplicateEvent = true;
      
      //ULong64_t runId = Run;
      //ULong64_t eventId = Event;
      //ULong64_t lumiId = LumiSect;
      /*
      for (unsigned int n = 0; n < runVec.size(); n++)
	{
	  if(runId == runVec[n] && lumiId == lumiVec[n] && eventId == eventVec[n]){notDuplicateEvent = false;}
	}
      */

      if (notDuplicateEvent)
	{

	  if(massZ1 < massZ1Cut) continue;
	  if(massZ2 < massZ2Cut) continue;
	  if(mass4l < m4lLowCut) continue;
	  if(mass4l > m4lHighCut) continue;

	  //if(mass4l > 110 && mass4l < 140) continue;
	  //if(mass4l > 300) continue;
	  
	  if( abs(idL1) == 13 && pTL1 < ptMuCut ) continue;
	  if( abs(idL2) == 13 && pTL2 < ptMuCut ) continue;
	  if( abs(idL3) == 13 && pTL3 < ptMuCut ) continue;
	  if( abs(idL4) == 13 && pTL4 < ptMuCut ) continue;

	  if( abs(idL1) == 11 && pTL1 < ptElCut ) continue;
	  if( abs(idL2) == 11 && pTL2 < ptElCut ) continue;
	  if( abs(idL3) == 11 && pTL3 < ptElCut ) continue;
	  if( abs(idL4) == 11 && pTL4 < ptElCut ) continue;

	  //runVec.push_back(runId);
	  //lumiVec.push_back(lumiId);
	  //eventVec.push_back(eventId);

	  //out <<Run<<":"<<LumiSect<<":"<<Event<<":";
	  //out.setf(ios::fixed,ios::floatfield);
	  //out.precision(2);
	  //out << mass4l <<":"<< massZ1 <<":"<< massZ2 <<":"
	  //   << massErrRaw <<":"<< massErrCorr <<":";
	  //out.setf(ios::fixed,ios::floatfield);
	  //out.precision(3);
	  //out << melaLD <<endl;
	  nJets = 0;

	  if((VBFJet1 && !VBFJet2) || (VBFJet2 && !VBFJet1)) nJets = 1;
	  if(!VBFJet2 && !VBFJet1) nJets = 0;
	  if(VBFJet1 && VBFJet2) nJets = 2;

	  if(!VBFJet1){pTVBFJet1 = 0; VBFDiJetMass = 0; VBFDeltaEta = 0;}
	  if(!VBFJet2){pTVBFJet2 = 0; VBFDiJetMass = 0; VBFDeltaEta = 0;}

	  out <<Run<<":"<<LumiSect<<":"<<Event<<":";
	  out.setf(ios::fixed,ios::floatfield);
	  out.precision(2);
	  out << mass4l <<":"<< massZ1 <<":"<< massZ2 <<":"
	      << massErrRaw <<":"<< massErrCorr <<":";
	  out.setf(ios::fixed,ios::floatfield);
	  out.precision(3);
	  out << JHUKD_H_qqZZ_noPDF<<":";
	  out.precision(2);
	  out << pT4l << ":";
	  out.precision(0);
	  out << nJets << ":";
	  out.precision(2);
	  out << pTVBFJet1 << ":" << pTVBFJet2 << ":" << VBFDiJetMass << ":";
	  out.precision(3);
	  out << VBFDeltaEta << ":" << JHUKD_H_h0M_noPDF << ":" << JHUKD_H_h0P_noPDF << ":" <<JHUKD_H_h1P_noPDF 
	      << ":" << JHUKD_H_h1M_noPDF << ":" << JHUKD_H_ggh2P_noPDF << ":" << JHUKD_H_qqh2P_noPDF <<endl;

	}
    
    }


  //runVec.clear();
  //lumiVec.clear();
  //eventVec.clear();  
  
  out.close();
}

