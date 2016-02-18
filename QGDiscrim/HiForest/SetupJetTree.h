//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Dec 18 07:MaxJets:09 2015 by ROOT version 6.05/03
// from TTree t/ak4PFpatJetsWithBtagging Jet Analysis Tree
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

class Jets {
public :
   Jets(){};
   ~Jets(){};

   // Declaration of leaf types
   Int_t           evt;
   Float_t         b;
   Int_t           nref;
   Float_t         rawpt[MaxJets];   //[nref]
   Float_t         jtpt[MaxJets];   //[nref]
   Float_t         jteta[MaxJets];   //[nref]
   Float_t         jty[MaxJets];   //[nref]
   Float_t         jtphi[MaxJets];   //[nref]
   Float_t         jtpu[MaxJets];   //[nref]
   Float_t         jtenergy[MaxJets];   //[nref]
   Float_t         jtm[MaxJets];   //[nref]
   Float_t         discr_fr01[MaxJets];   //[nref]
   Float_t         trackMax[MaxJets];   //[nref]
   Float_t         trackSum[MaxJets];   //[nref]
   Int_t           trackN[MaxJets];   //[nref]
   Float_t         trackHardSum[MaxJets];   //[nref]
   Int_t           trackHardN[MaxJets];   //[nref]
   Float_t         chargedMax[MaxJets];   //[nref]
   Float_t         chargedSum[MaxJets];   //[nref]
   Int_t           chargedN[MaxJets];   //[nref]
   Float_t         chargedHardSum[MaxJets];   //[nref]
   Int_t           chargedHardN[MaxJets];   //[nref]
   Float_t         photonMax[MaxJets];   //[nref]
   Float_t         photonSum[MaxJets];   //[nref]
   Int_t           photonN[MaxJets];   //[nref]
   Float_t         photonHardSum[MaxJets];   //[nref]
   Int_t           photonHardN[MaxJets];   //[nref]
   Float_t         neutralMax[MaxJets];   //[nref]
   Float_t         neutralSum[MaxJets];   //[nref]
   Int_t           neutralN[MaxJets];   //[nref]
   Float_t         eMax[MaxJets];   //[nref]
   Float_t         eSum[MaxJets];   //[nref]
   Int_t           eN[MaxJets];   //[nref]
   Float_t         muMax[MaxJets];   //[nref]
   Float_t         muSum[MaxJets];   //[nref]
   Int_t           muN[MaxJets];   //[nref]
   Float_t         discr_ssvHighEff[MaxJets];   //[nref]
   Float_t         discr_ssvHighPur[MaxJets];   //[nref]
   Float_t         discr_csvMva[MaxJets];   //[nref]
   Float_t         discr_csvSimple[MaxJets];   //[nref]
   Float_t         discr_muByIp3[MaxJets];   //[nref]
   Float_t         discr_muByPt[MaxJets];   //[nref]
   Float_t         discr_prob[MaxJets];   //[nref]
   Float_t         discr_probb[MaxJets];   //[nref]
   Float_t         discr_tcHighEff[MaxJets];   //[nref]
   Float_t         discr_tcHighPur[MaxJets];   //[nref]
   Float_t         ndiscr_ssvHighEff[MaxJets];   //[nref]
   Float_t         ndiscr_ssvHighPur[MaxJets];   //[nref]
   Float_t         ndiscr_csvSimple[MaxJets];   //[nref]
   Float_t         ndiscr_muByPt[MaxJets];   //[nref]
   Float_t         ndiscr_prob[MaxJets];   //[nref]
   Float_t         ndiscr_probb[MaxJets];   //[nref]
   Float_t         ndiscr_tcHighEff[MaxJets];   //[nref]
   Float_t         ndiscr_tcHighPur[MaxJets];   //[nref]
   Float_t         pdiscr_csvSimple[MaxJets];   //[nref]
   Float_t         pdiscr_prob[MaxJets];   //[nref]
   Float_t         pdiscr_probb[MaxJets];   //[nref]
   Int_t           nsvtx[MaxJets];   //[nref]
   Int_t           svtxntrk[MaxJets];   //[nref]
   Float_t         svtxdl[MaxJets];   //[nref]
   Float_t         svtxdls[MaxJets];   //[nref]
   Float_t         svtxdl2d[MaxJets];   //[nref]
   Float_t         svtxdls2d[MaxJets];   //[nref]
   Float_t         svtxm[MaxJets];   //[nref]
   Float_t         svtxpt[MaxJets];   //[nref]
   Int_t           nIPtrk[MaxJets];   //[nref]
   Int_t           nselIPtrk[MaxJets];   //[nref]
   Float_t         mue[MaxJets];   //[nref]
   Float_t         mupt[MaxJets];   //[nref]
   Float_t         mueta[MaxJets];   //[nref]
   Float_t         muphi[MaxJets];   //[nref]
   Float_t         mudr[MaxJets];   //[nref]
   Float_t         muptrel[MaxJets];   //[nref]
   Int_t           muchg[MaxJets];   //[nref]
   Int_t           beamId1;
   Int_t           beamId2;
   Float_t         pthat;
   Float_t         refpt[MaxJets];   //[nref]
   Float_t         refeta[MaxJets];   //[nref]
   Float_t         refy[MaxJets];   //[nref]
   Float_t         refphi[MaxJets];   //[nref]
   Float_t         refdphijt[MaxJets];   //[nref]
   Float_t         refdrjt[MaxJets];   //[nref]
   Float_t         refparton_pt[MaxJets];   //[nref]
   Int_t           refparton_flavor[MaxJets];   //[nref]
   Int_t           refparton_flavorForB[MaxJets];   //[nref]
   Float_t         genChargedSum[MaxJets];   //[nref]
   Float_t         genHardSum[MaxJets];   //[nref]
   Float_t         signalChargedSum[MaxJets];   //[nref]
   Float_t         signalHardSum[MaxJets];   //[nref]
   Int_t           ngen;
   Int_t           genmatchindex[MaxJets];   //[ngen]
   Float_t         genpt[MaxJets];   //[ngen]
   Float_t         geneta[MaxJets];   //[ngen]
   Float_t         geny[MaxJets];   //[ngen]
   Float_t         genphi[MaxJets];   //[ngen]
   Float_t         gendphijt[MaxJets];   //[ngen]
   Float_t         gendrjt[MaxJets];   //[ngen]

   // List of branches
   TBranch        *b_evt;   //!
   TBranch        *b_b;   //!
   TBranch        *b_nref;   //!
   TBranch        *b_rawpt;   //!
   TBranch        *b_jtpt;   //!
   TBranch        *b_jteta;   //!
   TBranch        *b_jty;   //!
   TBranch        *b_jtphi;   //!
   TBranch        *b_jtenergy;   //!
   TBranch        *b_jtpu;   //!
   TBranch        *b_jtm;   //!
   TBranch        *b_discr_fr01;   //!
   TBranch        *b_trackMax;   //!
   TBranch        *b_trackSum;   //!
   TBranch        *b_trackN;   //!
   TBranch        *b_trackHardSum;   //!
   TBranch        *b_trackHardN;   //!
   TBranch        *b_chargedMax;   //!
   TBranch        *b_chargedSum;   //!
   TBranch        *b_chargedN;   //!
   TBranch        *b_chargedHardSum;   //!
   TBranch        *b_chargedHardN;   //!
   TBranch        *b_photonMax;   //!
   TBranch        *b_photonSum;   //!
   TBranch        *b_photonN;   //!
   TBranch        *b_photonHardSum;   //!
   TBranch        *b_photonHardN;   //!
   TBranch        *b_neutralMax;   //!
   TBranch        *b_neutralSum;   //!
   TBranch        *b_neutralN;   //!
   TBranch        *b_eMax;   //!
   TBranch        *b_eSum;   //!
   TBranch        *b_eN;   //!
   TBranch        *b_muMax;   //!
   TBranch        *b_muSum;   //!
   TBranch        *b_muN;   //!
   TBranch        *b_discr_ssvHighEff;   //!
   TBranch        *b_discr_ssvHighPur;   //!
   TBranch        *b_discr_csvMva;   //!
   TBranch        *b_discr_csvSimple;   //!
   TBranch        *b_discr_muByIp3;   //!
   TBranch        *b_discr_muByPt;   //!
   TBranch        *b_discr_prob;   //!
   TBranch        *b_discr_probb;   //!
   TBranch        *b_discr_tcHighEff;   //!
   TBranch        *b_discr_tcHighPur;   //!
   TBranch        *b_ndiscr_ssvHighEff;   //!
   TBranch        *b_ndiscr_ssvHighPur;   //!
   TBranch        *b_ndiscr_csvSimple;   //!
   TBranch        *b_ndiscr_muByPt;   //!
   TBranch        *b_ndiscr_prob;   //!
   TBranch        *b_ndiscr_probb;   //!
   TBranch        *b_ndiscr_tcHighEff;   //!
   TBranch        *b_ndiscr_tcHighPur;   //!
   TBranch        *b_pdiscr_csvSimple;   //!
   TBranch        *b_pdiscr_prob;   //!
   TBranch        *b_pdiscr_probb;   //!
   TBranch        *b_nsvtx;   //!
   TBranch        *b_svtxntrk;   //!
   TBranch        *b_svtxdl;   //!
   TBranch        *b_svtxdls;   //!
   TBranch        *b_svtxdl2d;   //!
   TBranch        *b_svtxdls2d;   //!
   TBranch        *b_svtxm;   //!
   TBranch        *b_svtxpt;   //!
   TBranch        *b_nIPtrk;   //!
   TBranch        *b_nselIPtrk;   //!
   TBranch        *b_mue;   //!
   TBranch        *b_mupt;   //!
   TBranch        *b_mueta;   //!
   TBranch        *b_muphi;   //!
   TBranch        *b_mudr;   //!
   TBranch        *b_muptrel;   //!
   TBranch        *b_muchg;   //!
   TBranch        *b_beamId1;   //!
   TBranch        *b_beamId2;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_refpt;   //!
   TBranch        *b_refeta;   //!
   TBranch        *b_refy;   //!
   TBranch        *b_refphi;   //!
   TBranch        *b_refdphijt;   //!
   TBranch        *b_refdrjt;   //!
   TBranch        *b_refparton_pt;   //!
   TBranch        *b_refparton_flavor;   //!
   TBranch        *b_refparton_flavorForB;   //!
   TBranch        *b_genChargedSum;   //!
   TBranch        *b_genHardSum;   //!
   TBranch        *b_signalChargedSum;   //!
   TBranch        *b_signalHardSum;   //!
   TBranch        *b_ngen;   //!
   TBranch        *b_genmatchindex;   //!
   TBranch        *b_genpt;   //!
   TBranch        *b_geneta;   //!
   TBranch        *b_geny;   //!
   TBranch        *b_genphi;   //!
   TBranch        *b_gendphijt;   //!
   TBranch        *b_gendrjt;   //!

};


void setupJetTree(TTree *t,Jets &tJets,bool doCheck = 0)
{
   // Set branch addresses and branch pointers
   t->SetBranchAddress("evt", &tJets.evt, &tJets.b_evt);
   t->SetBranchAddress("b", &tJets.b, &tJets.b_b);
   t->SetBranchAddress("nref", &tJets.nref, &tJets.b_nref);
   t->SetBranchAddress("rawpt", tJets.rawpt, &tJets.b_rawpt);
   t->SetBranchAddress("jtpt", tJets.jtpt, &tJets.b_jtpt);
   t->SetBranchAddress("jteta", tJets.jteta, &tJets.b_jteta);
   t->SetBranchAddress("jty", tJets.jty, &tJets.b_jty);
   t->SetBranchAddress("jtphi", tJets.jtphi, &tJets.b_jtphi);
   t->SetBranchAddress("jtenergy", tJets.jtenergy, &tJets.b_jtenergy);
   t->SetBranchAddress("jtpu", tJets.jtpu, &tJets.b_jtpu);
   t->SetBranchAddress("jtm", tJets.jtm, &tJets.b_jtm);
   t->SetBranchAddress("discr_fr01", tJets.discr_fr01, &tJets.b_discr_fr01);
   t->SetBranchAddress("trackMax", tJets.trackMax, &tJets.b_trackMax);
   t->SetBranchAddress("trackSum", tJets.trackSum, &tJets.b_trackSum);
   t->SetBranchAddress("trackN", tJets.trackN, &tJets.b_trackN);
   t->SetBranchAddress("trackHardSum", tJets.trackHardSum, &tJets.b_trackHardSum);
   t->SetBranchAddress("trackHardN", tJets.trackHardN, &tJets.b_trackHardN);
   t->SetBranchAddress("chargedMax", tJets.chargedMax, &tJets.b_chargedMax);
   t->SetBranchAddress("chargedSum", tJets.chargedSum, &tJets.b_chargedSum);
   t->SetBranchAddress("chargedN", tJets.chargedN, &tJets.b_chargedN);
   t->SetBranchAddress("chargedHardSum", tJets.chargedHardSum, &tJets.b_chargedHardSum);
   t->SetBranchAddress("chargedHardN", tJets.chargedHardN, &tJets.b_chargedHardN);
   t->SetBranchAddress("photonMax", tJets.photonMax, &tJets.b_photonMax);
   t->SetBranchAddress("photonSum", tJets.photonSum, &tJets.b_photonSum);
   t->SetBranchAddress("photonN", tJets.photonN, &tJets.b_photonN);
   t->SetBranchAddress("photonHardSum", tJets.photonHardSum, &tJets.b_photonHardSum);
   t->SetBranchAddress("photonHardN", tJets.photonHardN, &tJets.b_photonHardN);
   t->SetBranchAddress("neutralMax", tJets.neutralMax, &tJets.b_neutralMax);
   t->SetBranchAddress("neutralSum", tJets.neutralSum, &tJets.b_neutralSum);
   t->SetBranchAddress("neutralN", tJets.neutralN, &tJets.b_neutralN);
   t->SetBranchAddress("eMax", tJets.eMax, &tJets.b_eMax);
   t->SetBranchAddress("eSum", tJets.eSum, &tJets.b_eSum);
   t->SetBranchAddress("eN", tJets.eN, &tJets.b_eN);
   t->SetBranchAddress("muMax", tJets.muMax, &tJets.b_muMax);
   t->SetBranchAddress("muSum", tJets.muSum, &tJets.b_muSum);
   t->SetBranchAddress("muN", tJets.muN, &tJets.b_muN);
   t->SetBranchAddress("discr_ssvHighEff", tJets.discr_ssvHighEff, &tJets.b_discr_ssvHighEff);
   t->SetBranchAddress("discr_ssvHighPur", tJets.discr_ssvHighPur, &tJets.b_discr_ssvHighPur);
   t->SetBranchAddress("discr_csvMva", tJets.discr_csvMva, &tJets.b_discr_csvMva);
   t->SetBranchAddress("discr_csvSimple", tJets.discr_csvSimple, &tJets.b_discr_csvSimple);
   t->SetBranchAddress("discr_muByIp3", tJets.discr_muByIp3, &tJets.b_discr_muByIp3);
   t->SetBranchAddress("discr_muByPt", tJets.discr_muByPt, &tJets.b_discr_muByPt);
   t->SetBranchAddress("discr_prob", tJets.discr_prob, &tJets.b_discr_prob);
   t->SetBranchAddress("discr_probb", tJets.discr_probb, &tJets.b_discr_probb);
   t->SetBranchAddress("discr_tcHighEff", tJets.discr_tcHighEff, &tJets.b_discr_tcHighEff);
   t->SetBranchAddress("discr_tcHighPur", tJets.discr_tcHighPur, &tJets.b_discr_tcHighPur);
   t->SetBranchAddress("ndiscr_ssvHighEff", tJets.ndiscr_ssvHighEff, &tJets.b_ndiscr_ssvHighEff);
   t->SetBranchAddress("ndiscr_ssvHighPur", tJets.ndiscr_ssvHighPur, &tJets.b_ndiscr_ssvHighPur);
   t->SetBranchAddress("ndiscr_csvSimple", tJets.ndiscr_csvSimple, &tJets.b_ndiscr_csvSimple);
   t->SetBranchAddress("ndiscr_muByPt", tJets.ndiscr_muByPt, &tJets.b_ndiscr_muByPt);
   t->SetBranchAddress("ndiscr_prob", tJets.ndiscr_prob, &tJets.b_ndiscr_prob);
   t->SetBranchAddress("ndiscr_probb", tJets.ndiscr_probb, &tJets.b_ndiscr_probb);
   t->SetBranchAddress("ndiscr_tcHighEff", tJets.ndiscr_tcHighEff, &tJets.b_ndiscr_tcHighEff);
   t->SetBranchAddress("ndiscr_tcHighPur", tJets.ndiscr_tcHighPur, &tJets.b_ndiscr_tcHighPur);
   t->SetBranchAddress("pdiscr_csvSimple", tJets.pdiscr_csvSimple, &tJets.b_pdiscr_csvSimple);
   t->SetBranchAddress("pdiscr_prob", tJets.pdiscr_prob, &tJets.b_pdiscr_prob);
   t->SetBranchAddress("pdiscr_probb", tJets.pdiscr_probb, &tJets.b_pdiscr_probb);
   t->SetBranchAddress("nsvtx", tJets.nsvtx, &tJets.b_nsvtx);
   t->SetBranchAddress("svtxntrk", tJets.svtxntrk, &tJets.b_svtxntrk);
   t->SetBranchAddress("svtxdl", tJets.svtxdl, &tJets.b_svtxdl);
   t->SetBranchAddress("svtxdls", tJets.svtxdls, &tJets.b_svtxdls);
   t->SetBranchAddress("svtxdl2d", tJets.svtxdl2d, &tJets.b_svtxdl2d);
   t->SetBranchAddress("svtxdls2d", tJets.svtxdls2d, &tJets.b_svtxdls2d);
   t->SetBranchAddress("svtxm", tJets.svtxm, &tJets.b_svtxm);
   t->SetBranchAddress("svtxpt", tJets.svtxpt, &tJets.b_svtxpt);
   t->SetBranchAddress("nIPtrk", tJets.nIPtrk, &tJets.b_nIPtrk);
   t->SetBranchAddress("nselIPtrk", tJets.nselIPtrk, &tJets.b_nselIPtrk);
   t->SetBranchAddress("mue", tJets.mue, &tJets.b_mue);
   t->SetBranchAddress("mupt", tJets.mupt, &tJets.b_mupt);
   t->SetBranchAddress("mueta", tJets.mueta, &tJets.b_mueta);
   t->SetBranchAddress("muphi", tJets.muphi, &tJets.b_muphi);
   t->SetBranchAddress("mudr", tJets.mudr, &tJets.b_mudr);
   t->SetBranchAddress("muptrel", tJets.muptrel, &tJets.b_muptrel);
   t->SetBranchAddress("muchg", tJets.muchg, &tJets.b_muchg);
   t->SetBranchAddress("beamId1", &tJets.beamId1, &tJets.b_beamId1);
   t->SetBranchAddress("beamId2", &tJets.beamId2, &tJets.b_beamId2);
   t->SetBranchAddress("pthat", &tJets.pthat, &tJets.b_pthat);
   t->SetBranchAddress("refpt", tJets.refpt, &tJets.b_refpt);
   t->SetBranchAddress("refeta", tJets.refeta, &tJets.b_refeta);
   t->SetBranchAddress("refy", tJets.refy, &tJets.b_refy);
   t->SetBranchAddress("refphi", tJets.refphi, &tJets.b_refphi);
   t->SetBranchAddress("refdphijt", tJets.refdphijt, &tJets.b_refdphijt);
   t->SetBranchAddress("refdrjt", tJets.refdrjt, &tJets.b_refdrjt);
   t->SetBranchAddress("refparton_pt", tJets.refparton_pt, &tJets.b_refparton_pt);
   t->SetBranchAddress("refparton_flavor", tJets.refparton_flavor, &tJets.b_refparton_flavor);
   t->SetBranchAddress("refparton_flavorForB", tJets.refparton_flavorForB, &tJets.b_refparton_flavorForB);
   t->SetBranchAddress("genChargedSum", tJets.genChargedSum, &tJets.b_genChargedSum);
   t->SetBranchAddress("genHardSum", tJets.genHardSum, &tJets.b_genHardSum);
   t->SetBranchAddress("signalChargedSum", tJets.signalChargedSum, &tJets.b_signalChargedSum);
   t->SetBranchAddress("signalHardSum", tJets.signalHardSum, &tJets.b_signalHardSum);
   t->SetBranchAddress("ngen", &tJets.ngen, &tJets.b_ngen);
   t->SetBranchAddress("genmatchindex", tJets.genmatchindex, &tJets.b_genmatchindex);
   t->SetBranchAddress("genpt", tJets.genpt, &tJets.b_genpt);
   t->SetBranchAddress("geneta", tJets.geneta, &tJets.b_geneta);
   t->SetBranchAddress("geny", tJets.geny, &tJets.b_geny);
   t->SetBranchAddress("genphi", tJets.genphi, &tJets.b_genphi);
   t->SetBranchAddress("gendphijt", tJets.gendphijt, &tJets.b_gendphijt);
   t->SetBranchAddress("gendrjt", tJets.gendrjt, &tJets.b_gendrjt);
   if (doCheck) {
      if (t->GetMaximum("nref")>MaxJets) std::cout <<"FATAL ERROR: Arrary size of nref too small!!!  "<<t->GetMaximum("nref")<<std::endl;
      if (t->GetMaximum("ngen")>MaxJets) std::cout <<"FATAL ERROR: Arrary size of ngen too small!!!  "<<t->GetMaximum("ngen")<<std::endl;
   }
}

