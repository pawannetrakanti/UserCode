//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Dec 18 07:26:09 2015 by ROOT version 6.05/03
// from TTree HltTree/
// found on file: root://eoscms//eos/cms/store/group/phys_heavyions/azsigmon/HiForestAODpp5TeV/HiForestAOD_withTupel_pp_MC_Z30mumuJet_v1.root
//////////////////////////////////////////////////////////
#include "commonSetup.h"
#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

class Skims {
public :
   Skims(){};
   ~Skims(){};

   // Declaration of leaf types
   Int_t           ana_step;
   Int_t           pHBHENoiseFilterResultProducer;
   Int_t           pPAprimaryVertexFilter;
   Int_t           pBeamScrapingFilter;
   Int_t           pVertexFilterCutG;
   Int_t           pVertexFilterCutGloose;
   Int_t           pVertexFilterCutGtight;
   Int_t           pVertexFilterCutGplus;
   Int_t           pVertexFilterCutE;
   Int_t           pVertexFilterCutEandG;
   Int_t           HBHENoiseFilterResultRun1;
   Int_t           HBHENoiseFilterResultRun2Loose;
   Int_t           HBHENoiseFilterResultRun2Tight;
   Int_t           HBHENoiseFilterResult;
   Int_t           HBHEIsoNoiseFilterResult;

   // List of branches
   TBranch        *b_ana_step;   //!
   TBranch        *b_pHBHENoiseFilterResultProducer;   //!
   TBranch        *b_pPAprimaryVertexFilter;   //!
   TBranch        *b_pBeamScrapingFilter;   //!
   TBranch        *b_pVertexFilterCutG;   //!
   TBranch        *b_pVertexFilterCutGloose;   //!
   TBranch        *b_pVertexFilterCutGtight;   //!
   TBranch        *b_pVertexFilterCutGplus;   //!
   TBranch        *b_pVertexFilterCutE;   //!
   TBranch        *b_pVertexFilterCutEandG;   //!
   TBranch        *b_HBHENoiseFilterResultRun1;   //!
   TBranch        *b_HBHENoiseFilterResultRun2Loose;   //!
   TBranch        *b_HBHENoiseFilterResultRun2Tight;   //!
   TBranch        *b_HBHENoiseFilterResult;   //!
   TBranch        *b_HBHEIsoNoiseFilterResult;   //!

};


void setupSkimTree(TTree *t,Skims &tSkims,bool doCheck = 0)
{
   // Set branch addresses and branch pointers
   t->SetBranchAddress("ana_step", &tSkims.ana_step, &tSkims.b_ana_step);
   t->SetBranchAddress("pHBHENoiseFilterResultProducer", &tSkims.pHBHENoiseFilterResultProducer, &tSkims.b_pHBHENoiseFilterResultProducer);
   t->SetBranchAddress("pPAprimaryVertexFilter", &tSkims.pPAprimaryVertexFilter, &tSkims.b_pPAprimaryVertexFilter);
   t->SetBranchAddress("pBeamScrapingFilter", &tSkims.pBeamScrapingFilter, &tSkims.b_pBeamScrapingFilter);
   t->SetBranchAddress("pVertexFilterCutG", &tSkims.pVertexFilterCutG, &tSkims.b_pVertexFilterCutG);
   t->SetBranchAddress("pVertexFilterCutGloose", &tSkims.pVertexFilterCutGloose, &tSkims.b_pVertexFilterCutGloose);
   t->SetBranchAddress("pVertexFilterCutGtight", &tSkims.pVertexFilterCutGtight, &tSkims.b_pVertexFilterCutGtight);
   t->SetBranchAddress("pVertexFilterCutGplus", &tSkims.pVertexFilterCutGplus, &tSkims.b_pVertexFilterCutGplus);
   t->SetBranchAddress("pVertexFilterCutE", &tSkims.pVertexFilterCutE, &tSkims.b_pVertexFilterCutE);
   t->SetBranchAddress("pVertexFilterCutEandG", &tSkims.pVertexFilterCutEandG, &tSkims.b_pVertexFilterCutEandG);
   t->SetBranchAddress("HBHENoiseFilterResultRun1", &tSkims.HBHENoiseFilterResultRun1, &tSkims.b_HBHENoiseFilterResultRun1);
   t->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &tSkims.HBHENoiseFilterResultRun2Loose, &tSkims.b_HBHENoiseFilterResultRun2Loose);
   t->SetBranchAddress("HBHENoiseFilterResultRun2Tight", &tSkims.HBHENoiseFilterResultRun2Tight, &tSkims.b_HBHENoiseFilterResultRun2Tight);
   t->SetBranchAddress("HBHENoiseFilterResult", &tSkims.HBHENoiseFilterResult, &tSkims.b_HBHENoiseFilterResult);
   t->SetBranchAddress("HBHEIsoNoiseFilterResult", &tSkims.HBHEIsoNoiseFilterResult, &tSkims.b_HBHEIsoNoiseFilterResult);
   if (doCheck) {
   }
}

