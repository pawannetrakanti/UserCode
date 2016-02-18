//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Dec 18 07:26:09 2015 by ROOT version 6.05/03
// from TTree pfTree/dijet tree
// found on file: root://eoscms//eos/cms/store/group/phys_heavyions/azsigmon/HiForestAODpp5TeV/HiForestAOD_withTupel_pp_MC_Z30mumuJet_v1.root
//////////////////////////////////////////////////////////
#include "commonSetup.h"
#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

using namespace std;

class PFCands {
public :
   PFCands(){};
   ~PFCands(){};

   // Declaration of leaf types
   Int_t           nPFpart;
   Int_t           pfId[MaxEntry];   //[nPFpart]
   Float_t         pfPt[MaxEntry];   //[nPFpart]
   Float_t         pfVsPtInitial[MaxEntry];   //[nPFpart]
   Float_t         pfEta[MaxEntry];   //[nPFpart]
   Float_t         pfPhi[MaxEntry];   //[nPFpart]
   Float_t         pfEnergy[MaxEntry];   //[nPFpart]
   Float_t         vn  [5][15];
   Float_t         psin[5][15];
   Float_t         sumpt[15];

   // List of branches
   TBranch        *b_nPFpart;   //!
   TBranch        *b_pfId;   //!
   TBranch        *b_pfPt;   //!
   TBranch        *b_pfVsPtInitial;   //!
   TBranch        *b_pfEta;   //!
   TBranch        *b_pfPhi;   //!
   TBranch        *b_pfEnergy;   //!
   TBranch        *b_vn;   //!
   TBranch        *b_vpsi;   //!
   TBranch        *b_sumpt;   //!

};


void setupPFCandTree(TTree *t,PFCands &tPFCands,bool doCheck = 0)
{
   // Set branch addresses and branch pointers
   t->SetBranchAddress("nPFpart", &tPFCands.nPFpart, &tPFCands.b_nPFpart);
   t->SetBranchAddress("pfId", tPFCands.pfId, &tPFCands.b_pfId);
   t->SetBranchAddress("pfPt", tPFCands.pfPt, &tPFCands.b_pfPt);
   t->SetBranchAddress("pfVsPtInitial", tPFCands.pfVsPtInitial, &tPFCands.b_pfVsPtInitial);
   t->SetBranchAddress("pfEta", tPFCands.pfEta, &tPFCands.b_pfEta);
   t->SetBranchAddress("pfPhi", tPFCands.pfPhi, &tPFCands.b_pfPhi);
   t->SetBranchAddress("pfEnergy", tPFCands.pfEnergy, &tPFCands.b_pfEnergy);
   t->SetBranchAddress("vn", tPFCands.vn, &tPFCands.b_vn);
   t->SetBranchAddress("psin", tPFCands.psin, &tPFCands.b_vpsi);
   t->SetBranchAddress("sumpt", tPFCands.sumpt, &tPFCands.b_sumpt);
   if (doCheck) {
     if (t->GetMaximum("nPFpart")>MaxEntry)std::cout <<"FATAL ERROR: Arrary size of nPFpart too small!!!  "<<t->GetMaximum("nPFpart")<<std::endl;
   }
}

