//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Dec 14 11:54:48 2015 by ROOT version 6.04/00
// from TTree trackTree/v1
// found on file: ../Data/HiForest-HighPtJet80-RECO753p1.root
//////////////////////////////////////////////////////////
#include "commonSetup.h"
#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

using namespace std;

class Tracks {
public :
   Tracks(){};
   ~Tracks(){};

   // Declaration of leaf types
   Int_t           nEv;
   Int_t           nLumi;
   Int_t           nBX;
   Int_t           nRun;
   Int_t           N;
   Int_t           nVtx;
   Int_t           nTrk;
   Int_t           maxPtVtx;
   Int_t           maxMultVtx;
   Int_t           nTrkVtx[50];   //[nVtx]
   Float_t         normChi2Vtx[50];   //[nVtx]
   Float_t         sumPtVtx[50];   //[nVtx]
   Float_t         xVtx[50];   //[nVtx]
   Float_t         yVtx[50];   //[nVtx]
   Float_t         zVtx[50];   //[nVtx]
   Float_t         xVtxErr[50];   //[nVtx]
   Float_t         yVtxErr[50];   //[nVtx]
   Float_t         zVtxErr[50];   //[nVtx]
   Float_t         vtxDist2D[50];   //[nVtx]
   Float_t         vtxDist2DErr[50];   //[nVtx]
   Float_t         vtxDist2DSig[50];   //[nVtx]
   //Float_t         vtxDist3D[50];   //[nVtx]
   Float_t         vtxDist3DErr[50];   //[nVtx]
   Float_t         vtxDist3DSig[50];   //[nVtx]
   Int_t           nVtxSim;
   Float_t         xVtxSim[2];   //[nVtxSim]
   Float_t         yVtxSim[2];   //[nVtxSim]
   Float_t         zVtxSim[2];   //[nVtxSim]
   Float_t         trkPt[MaxEntry];   //[nTrk]
   Float_t         trkPtError[MaxEntry];   //[nTrk]
   UChar_t         trkNHit[MaxEntry];   //[nTrk]
   UChar_t         trkNlayer[MaxEntry];   //[nTrk]
   Float_t         trkEta[MaxEntry];   //[nTrk]
   Float_t         trkPhi[MaxEntry];   //[nTrk]
   Int_t           trkCharge[MaxEntry];   //[nTrk]
   UChar_t         trkNVtx[MaxEntry];   //[nTrk]
   Int_t           nTrkTimesnVtx;
   Bool_t          trkAssocVtx[MaxNumber];   //[nTrkTimesnVtx]
   Bool_t          highPurity[MaxEntry];   //[nTrk]
   Bool_t          tight[MaxEntry];   //[nTrk]
   Bool_t          loose[MaxEntry];   //[nTrk]
   Float_t         trkChi2[MaxEntry];   //[nTrk]
   UChar_t         trkNdof[MaxEntry];   //[nTrk]
   Float_t         trkDxy1[MaxEntry];   //[nTrk]
   Float_t         trkDxyError1[MaxEntry];   //[nTrk]
   Float_t         trkDz1[MaxEntry];   //[nTrk]
   Float_t         trkDzError1[MaxEntry];   //[nTrk]
   Bool_t          trkFake[MaxEntry];   //[nTrk]
   UChar_t         trkAlgo[MaxEntry];   //[nTrk]
   UChar_t         trkOriginalAlgo[MaxEntry];   //[nTrk]
   Float_t         trkMVA[MaxEntry];   //[nTrk]
   Int_t           pfType[MaxEntry];   //[nTrk]
   Float_t         pfCandPt[MaxEntry];   //[nTrk]
   Float_t         pfEcal[MaxEntry];   //[nTrk]
   Float_t         pfHcal[MaxEntry];   //[nTrk]

   // List of branches
   TBranch        *b_nEv;   //!
   TBranch        *b_nLumi;   //!
   TBranch        *b_nBX;   //!
   TBranch        *b_nRun;   //!
   TBranch        *b_N;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_nTrk;   //!
   TBranch        *b_maxPtVtx;   //!
   TBranch        *b_maxMultVtx;   //!
   TBranch        *b_nTrkVtx;   //!
   TBranch        *b_normChi2Vtx;   //!
   TBranch        *b_sumPtVtx;   //!
   TBranch        *b_xVtx;   //!
   TBranch        *b_yVtx;   //!
   TBranch        *b_zVtx;   //!
   TBranch        *b_xVtxErr;   //!
   TBranch        *b_yVtxErr;   //!
   TBranch        *b_zVtxErr;   //!
   TBranch        *b_vtxDist2D;   //!
   TBranch        *b_vtxDist2DErr;   //!
   TBranch        *b_vtxDist2DSig;   //!
   //TBranch        *b_vtxDist3D;   //!
   TBranch        *b_vtxDist3DErr;   //!
   TBranch        *b_vtxDist3DSig;   //!
   TBranch        *b_nVtxSim;   //!
   TBranch        *b_xVtxSim;   //!
   TBranch        *b_yVtxSim;   //!
   TBranch        *b_zVtxSim;   //!
   TBranch        *b_trkPt;   //!
   TBranch        *b_trkPtError;   //!
   TBranch        *b_trkNHit;   //!
   TBranch        *b_trkNlayer;   //!
   TBranch        *b_trkEta;   //!
   TBranch        *b_trkPhi;   //!
   TBranch        *b_trkCharge;   //!
   TBranch        *b_trkNVtx;   //!
   TBranch        *b_nTrkTimesnVtx;   //!
   TBranch        *b_trkAssocVtx;   //!
   TBranch        *b_highPurity;   //!
   TBranch        *b_tight;   //!
   TBranch        *b_loose;   //!
   TBranch        *b_trkChi2;   //!
   TBranch        *b_trkNdof;   //!
   TBranch        *b_trkDxy1;   //!
   TBranch        *b_trkDxyError1;   //!
   TBranch        *b_trkDz1;   //!
   TBranch        *b_trkDzError1;   //!
   TBranch        *b_trkFake;   //!
   TBranch        *b_trkAlgo;   //!
   TBranch        *b_trkOriginalAlgo;   //!
   TBranch        *b_trkMVA;   //!
   TBranch        *b_pfType;   //!
   TBranch        *b_pfCandPt;   //!
   TBranch        *b_pfEcal;   //!
   TBranch        *b_pfHcal;   //!

};


void setupTrackTree(TTree *t,Tracks &tTracks,bool doCheck = 0)
{
   // Set branch addresses and branch pointers
   t->SetBranchAddress("nEv", &tTracks.nEv, &tTracks.b_nEv);
   t->SetBranchAddress("nLumi", &tTracks.nLumi, &tTracks.b_nLumi);
   t->SetBranchAddress("nBX", &tTracks.nBX, &tTracks.b_nBX);
   t->SetBranchAddress("nRun", &tTracks.nRun, &tTracks.b_nRun);
   t->SetBranchAddress("N", &tTracks.N, &tTracks.b_N);
   t->SetBranchAddress("nVtx", &tTracks.nVtx, &tTracks.b_nVtx);
   t->SetBranchAddress("nTrk", &tTracks.nTrk, &tTracks.b_nTrk);
   t->SetBranchAddress("maxPtVtx", &tTracks.maxPtVtx, &tTracks.b_maxPtVtx);
   t->SetBranchAddress("maxMultVtx", &tTracks.maxMultVtx, &tTracks.b_maxMultVtx);
   t->SetBranchAddress("nTrkVtx", tTracks.nTrkVtx, &tTracks.b_nTrkVtx);
   t->SetBranchAddress("normChi2Vtx", tTracks.normChi2Vtx, &tTracks.b_normChi2Vtx);
   t->SetBranchAddress("sumPtVtx", tTracks.sumPtVtx, &tTracks.b_sumPtVtx);
   t->SetBranchAddress("xVtx", tTracks.xVtx, &tTracks.b_xVtx);
   t->SetBranchAddress("yVtx", tTracks.yVtx, &tTracks.b_yVtx);
   t->SetBranchAddress("zVtx", tTracks.zVtx, &tTracks.b_zVtx);
   t->SetBranchAddress("xVtxErr", tTracks.xVtxErr, &tTracks.b_xVtxErr);
   t->SetBranchAddress("yVtxErr", tTracks.yVtxErr, &tTracks.b_yVtxErr);
   t->SetBranchAddress("zVtxErr", tTracks.zVtxErr, &tTracks.b_zVtxErr);
   t->SetBranchAddress("vtxDist2D", tTracks.vtxDist2D, &tTracks.b_vtxDist2D);
   t->SetBranchAddress("vtxDist2DErr", tTracks.vtxDist2DErr, &tTracks.b_vtxDist2DErr);
   t->SetBranchAddress("vtxDist2DSig", tTracks.vtxDist2DSig, &tTracks.b_vtxDist2DSig);
   //t->SetBranchAddress("vtxDist3D", tTracks.vtxDist3D, &tTracks.b_vtxDist3D);
   t->SetBranchAddress("vtxDist3DErr", tTracks.vtxDist3DErr, &tTracks.b_vtxDist3DErr);
   t->SetBranchAddress("vtxDist3DSig", tTracks.vtxDist3DSig, &tTracks.b_vtxDist3DSig);
   t->SetBranchAddress("nVtxSim", &tTracks.nVtxSim, &tTracks.b_nVtxSim);
   t->SetBranchAddress("xVtxSim", &tTracks.xVtxSim, &tTracks.b_xVtxSim);
   t->SetBranchAddress("yVtxSim", &tTracks.yVtxSim, &tTracks.b_yVtxSim);
   t->SetBranchAddress("zVtxSim", &tTracks.zVtxSim, &tTracks.b_zVtxSim);
   t->SetBranchAddress("trkPt", tTracks.trkPt, &tTracks.b_trkPt);
   t->SetBranchAddress("trkPtError", tTracks.trkPtError, &tTracks.b_trkPtError);
   t->SetBranchAddress("trkNHit", tTracks.trkNHit, &tTracks.b_trkNHit);
   t->SetBranchAddress("trkNlayer", tTracks.trkNlayer, &tTracks.b_trkNlayer);
   t->SetBranchAddress("trkEta", tTracks.trkEta, &tTracks.b_trkEta);
   t->SetBranchAddress("trkPhi", tTracks.trkPhi, &tTracks.b_trkPhi);
   t->SetBranchAddress("trkCharge", tTracks.trkCharge, &tTracks.b_trkCharge);
   t->SetBranchAddress("trkNVtx", tTracks.trkNVtx, &tTracks.b_trkNVtx);
   t->SetBranchAddress("nTrkTimesnVtx", &tTracks.nTrkTimesnVtx, &tTracks.b_nTrkTimesnVtx);
   t->SetBranchAddress("trkAssocVtx", tTracks.trkAssocVtx, &tTracks.b_trkAssocVtx);
   t->SetBranchAddress("highPurity", tTracks.highPurity, &tTracks.b_highPurity);
   t->SetBranchAddress("tight", tTracks.tight, &tTracks.b_tight);
   t->SetBranchAddress("loose", tTracks.loose, &tTracks.b_loose);
   t->SetBranchAddress("trkChi2", tTracks.trkChi2, &tTracks.b_trkChi2);
   t->SetBranchAddress("trkNdof", tTracks.trkNdof, &tTracks.b_trkNdof);
   t->SetBranchAddress("trkDxy1", tTracks.trkDxy1, &tTracks.b_trkDxy1);
   t->SetBranchAddress("trkDxyError1", tTracks.trkDxyError1, &tTracks.b_trkDxyError1);
   t->SetBranchAddress("trkDz1", tTracks.trkDz1, &tTracks.b_trkDz1);
   t->SetBranchAddress("trkDzError1", tTracks.trkDzError1, &tTracks.b_trkDzError1);
   t->SetBranchAddress("trkFake", tTracks.trkFake, &tTracks.b_trkFake);
   t->SetBranchAddress("trkAlgo", tTracks.trkAlgo, &tTracks.b_trkAlgo);
   t->SetBranchAddress("trkOriginalAlgo", tTracks.trkOriginalAlgo, &tTracks.b_trkOriginalAlgo);
   t->SetBranchAddress("trkMVA", tTracks.trkMVA, &tTracks.b_trkMVA);
   t->SetBranchAddress("pfType", tTracks.pfType, &tTracks.b_pfType);
   t->SetBranchAddress("pfCandPt", tTracks.pfCandPt, &tTracks.b_pfCandPt);
   t->SetBranchAddress("pfEcal", tTracks.pfEcal, &tTracks.b_pfEcal);
   t->SetBranchAddress("pfHcal", tTracks.pfHcal, &tTracks.b_pfHcal);
   if (doCheck) {
     if (t->GetMaximum("nVtx")>50) std::cout <<"FATAL ERROR: Arrary size of nVtx too small!!!  "<<t->GetMaximum("nVtx")<<std::endl;
      if (t->GetMaximum("nVtxSim")>2) std::cout <<"FATAL ERROR: Arrary size of nVtxSim too small!!!  "<<t->GetMaximum("nVtxSim")<<std::endl;
      if (t->GetMaximum("nTrk")>MaxEntry) std::cout <<"FATAL ERROR: Arrary size of nTrk too small!!!  "<<t->GetMaximum("nTrk")<<std::endl;
      if (t->GetMaximum("nTrkTimesnVtx")>MaxNumber) std::cout <<"FATAL ERROR: Arrary size of nTrkTimesnVtx too small!!!  "<<t->GetMaximum("nTrkTimesnVtx")<<std::endl;
   }
}

