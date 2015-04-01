#include <TF1.h>
#include <TProfile.h>
#include <TChain.h>
#include <TTree.h>
#include <TEventList.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom3.h>

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <cstring>
#include <sstream>
#include <vector>
#include <algorithm>


using namespace std;

#ifdef _MAKECINT_
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#endif



static const int nbins_pt = 39;
static const double boundaries_pt[nbins_pt+1] = {
  3, 4, 5, 7, 9, 12,
  15, 18, 21, 24, 28,
  32, 37, 43, 49, 56,
  64, 74, 84, 97, 114,
  133, 153, 174, 196,
  220, 245, 272, 300,
  330, 362, 395, 430,
  468, 507, 548, 592,
  638, 686, 1000
};


int trigger_turnon(){

  TFile *fin = new TFile("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PbPb/2011/data/PbPb_MinBiasUPC_ntuple_SuperNovaRejected_akPuCalo_20150311.root","r");


  TTree *trig = (TTree*)fin->Get("trigger_info");

  int   run_value;	 
  int   evt_value;	 
  int   lumi_value;	 
  int   Jet80;		 
  int   Jet80_prescl;	 
  int   Jet65;		 
  int   Jet65_prescl;	 
  int   Jet55;		 
  int   Jet55_prescl;	 
  int   HIMinBias;	 
  int   HIMinBias_prescl;

  int   nref;		 
  float jetpt[100];		 
  float jeteta[100];		 
    

  trig->SetBranchAddress("run_value",&run_value);	 
  trig->SetBranchAddress("evt_value",&evt_value);	 
  trig->SetBranchAddress("lumi_value",&lumi_value);	 
  trig->SetBranchAddress("Jet80",&Jet80);		 
  trig->SetBranchAddress("Jet80_prescl",&Jet80_prescl);	 
  trig->SetBranchAddress("Jet65",&Jet65);		 
  trig->SetBranchAddress("Jet65_prescl",&Jet65_prescl);	 
  trig->SetBranchAddress("Jet55",&Jet55);		 
  trig->SetBranchAddress("Jet55_prescl",&Jet55_prescl);	 
  trig->SetBranchAddress("HIMinBias",&HIMinBias);	 
  trig->SetBranchAddress("HIMinBias_prescl",&HIMinBias_prescl);
  trig->SetBranchAddress("nref",&nref);		 
  trig->SetBranchAddress("jetpt",jetpt);		 
  trig->SetBranchAddress("jeteta",jeteta);		 

  TFile *fout = new TFile("trigger_histos_pbpb.root","recreate");

  TH1D *hevt = new TH1D("hevt","# of events",30,-0.5,30-0.5);
  hevt->Sumw2();
  TH2D *hjet = new TH2D("hjet","jet pt for different trigger comb",30,-0.5,30-0.5,nbins_pt,boundaries_pt);
  hjet->Sumw2();


  Long64_t nentries = trig->GetEntries();
  Long64_t nbytes=0;
  for(int iev=0; iev<nentries; iev++){
    nbytes += trig->GetEntry(iev);

    if(HIMinBias)hevt->Fill(0);
    if(Jet55)hevt->Fill(1);
    if(Jet65)hevt->Fill(2);
    if(Jet80)hevt->Fill(3);

    if(HIMinBias_prescl)hevt->Fill(4);
    if(Jet55_prescl)hevt->Fill(5);
    if(Jet65_prescl)hevt->Fill(6);
    if(Jet80_prescl)hevt->Fill(7);
    
    if(HIMinBias && Jet55)hevt->Fill(8);
    if(HIMinBias && Jet65)hevt->Fill(9);
    if(HIMinBias && Jet80)hevt->Fill(10);


    if(HIMinBias && Jet55_prescl)hevt->Fill(11);
    if(HIMinBias && Jet65_prescl)hevt->Fill(12);
    if(HIMinBias && Jet80_prescl)hevt->Fill(13);

    if(HIMinBias_prescl && Jet55_prescl)hevt->Fill(14);
    if(HIMinBias_prescl && Jet65_prescl)hevt->Fill(15);
    if(HIMinBias_prescl && Jet80_prescl)hevt->Fill(16);

    if(HIMinBias && Jet55 && Jet55_prescl)hevt->Fill(17);
    if(HIMinBias && Jet65 && Jet65_prescl)hevt->Fill(18);
    if(HIMinBias && Jet80 && Jet80_prescl)hevt->Fill(19);

    if(HIMinBias && (Jet55 || Jet55_prescl))hevt->Fill(20);
    if(HIMinBias && (Jet65 || Jet65_prescl))hevt->Fill(21);
    if(HIMinBias && (Jet80 || Jet80_prescl))hevt->Fill(22);

    if(HIMinBias && Jet55 && Jet65)hevt->Fill(23);
    if(HIMinBias && Jet65 && Jet80)hevt->Fill(24);
    if(HIMinBias && Jet55 && !Jet65)hevt->Fill(25);
    if(HIMinBias && Jet65 && !Jet80)hevt->Fill(26);

    if(HIMinBias && HIMinBias_prescl)hevt->Fill(27);    


    for(int i=0; i<nref; i++){ 

      if(HIMinBias)hjet->Fill(0.,jetpt[i]);
      if(HIMinBias_prescl)hjet->Fill(1.,jetpt[i]);

      if(HIMinBias && Jet80)hjet->Fill(2.,jetpt[i]);
      if(HIMinBias_prescl && Jet80)hjet->Fill(3.,jetpt[i]);
      if(HIMinBias && Jet80_prescl)hjet->Fill(4.,jetpt[i]);
      if(HIMinBias_prescl && Jet80_prescl)hjet->Fill(5.,jetpt[i]);
      if(HIMinBias && Jet80 && Jet80_prescl)hjet->Fill(6.,jetpt[i]);
      if(HIMinBias && (Jet80 || Jet80_prescl))hjet->Fill(7.,jetpt[i]);
      if(HIMinBias && Jet65)hjet->Fill(8.,jetpt[i]);
      

      if(HIMinBias_prescl && Jet65)hjet->Fill(9.,jetpt[i]);
      if(HIMinBias && Jet65_prescl)hjet->Fill(10.,jetpt[i]);
      if(HIMinBias_prescl && Jet65_prescl)hjet->Fill(11.,jetpt[i]);
      if(HIMinBias && Jet65 && Jet65_prescl)hjet->Fill(12.,jetpt[i]);
      if(HIMinBias && (Jet65 || Jet65_prescl))hjet->Fill(13.,jetpt[i]);


      if(HIMinBias && Jet55)hjet->Fill(14.,jetpt[i]);
      if(HIMinBias_prescl && Jet55)hjet->Fill(15.,jetpt[i]);
      if(HIMinBias && Jet55_prescl)hjet->Fill(16.,jetpt[i]);
      if(HIMinBias_prescl && Jet55_prescl)hjet->Fill(17.,jetpt[i]);
      if(HIMinBias && Jet55 && Jet55_prescl)hjet->Fill(18.,jetpt[i]);
      if(HIMinBias && (Jet55 || Jet55_prescl))hjet->Fill(19.,jetpt[i]);

      if(HIMinBias && HIMinBias_prescl)hjet->Fill(20.,jetpt[i]);


      if(HIMinBias && Jet65 && !Jet80)hjet->Fill(21.,jetpt[i]);
      if(HIMinBias && Jet55 && !Jet65 && !Jet80)hjet->Fill(22.,jetpt[i]);

    }
  }


  fout->cd();
  fout->Write();
  fout->Close();


  return 0;

}
