//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Dec 15 07:55:04 2015 by ROOT version 6.05/03
// from TTree EventTree/ EventTree
// found on file: root://eoscms//eos/cms//store/group/phys_heavyions/azsigmon/HiForestAODPbPb5TeV/HiForestAOD_withTupel_PbPb_MC_Z30mumuJet_v1.root
//////////////////////////////////////////////////////////

#ifndef EventTree_h
#define EventTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <TObject.h>

// Header file for the classes stored in the TTree if any.
#include <iostream>
#include <vector>

using namespace std;

class EventTree : public TObject {

public :

  EventTree(const char *name);
  EventTree(const char *name,const char *filename);
  EventTree(const char *name, TTree *inTree);
  virtual ~EventTree();


  virtual Long64_t GetEntries();
  virtual Int_t    GetEntry(Long64_t entry);
  virtual void     Init(TTree *tree);


  TFile          *fInfile;
  TTree          *fTree;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Fixed size dimensions of array or collections stored in the TTree if any.

  // Declaration of leaf types
  Int_t           EvtIsRealData;
  UInt_t          EvtNum;
  UInt_t          EvtRunNum;
  Int_t           EvtLumiNum;
  Int_t           EvtBxNum;
  Int_t           EvtVtxCnt;
  Int_t           EvtPuCnt;
  Int_t           EvtPuCntTruth;
  std::vector<double>  *EvtWeights;
  UInt_t          TrigHlt;
  ULong64_t       TrigHltMu;
  ULong64_t       TrigHltDiMu;
  ULong64_t       TrigHltEl;
  ULong64_t       TrigHltDiEl;
  ULong64_t       TrigHltElMu;
  std::vector<float>   *METPt;
  std::vector<float>   *METPx;
  std::vector<float>   *METPy;
  std::vector<float>   *METPz;
  std::vector<float>   *METE;
  std::vector<float>   *METsigx2;
  std::vector<float>   *METsigxy;
  std::vector<float>   *METsigy2;
  std::vector<float>   *METsig;
  std::vector<float>   *GLepDr01Pt;
  std::vector<float>   *GLepDr01Eta;
  std::vector<float>   *GLepDr01Phi;
  std::vector<float>   *GLepDr01E;
  std::vector<int>     *GLepDr01Id;
  std::vector<int>     *GLepDr01St;
  std::vector<int>     *GLepDr01MomId;
  std::vector<float>   *GLepBarePt;
  std::vector<float>   *GLepBareEta;
  std::vector<float>   *GLepBarePhi;
  std::vector<float>   *GLepBareE;
  std::vector<int>     *GLepBareId;
  std::vector<int>     *GLepBareSt;
  std::vector<int>     *GLepBareMomId;
  std::vector<float>   *GLepSt3Pt;
  std::vector<float>   *GLepSt3Eta;
  std::vector<float>   *GLepSt3Phi;
  std::vector<float>   *GLepSt3E;
  std::vector<int>     *GLepSt3Id;
  std::vector<int>     *GLepSt3St;
  std::vector<int>     *GLepSt3Mother0Id;
  std::vector<int>     *GLepSt3MotherCnt;
  std::vector<float>   *GPhotPt;
  std::vector<float>   *GPhotEta;
  std::vector<float>   *GPhotPhi;
  std::vector<float>   *GPhotE;
  std::vector<int>     *GPhotMotherId;
  std::vector<int>     *GPhotSt;
  std::vector<float>   *GPhotIsoEDR03;
  std::vector<float>   *GPhotIsoEDR04;
  std::vector<float>   *GPhotIsoEDR05;
  std::vector<float>   *GPhotIsoSumPtDR03;
  std::vector<float>   *GPhotIsoSumPtDR04;
  std::vector<float>   *GPhotIsoSumPtDR05;
  std::vector<float>   *GLepClosePhotPt;
  std::vector<float>   *GLepClosePhotEta;
  std::vector<float>   *GLepClosePhotPhi;
  std::vector<float>   *GLepClosePhotE;
  std::vector<int>     *GLepClosePhotId;
  std::vector<int>     *GLepClosePhotMother0Id;
  std::vector<int>     *GLepClosePhotMotherCnt;
  std::vector<int>     *GLepClosePhotSt;
  std::vector<float>   *GJetAk04Pt;
  std::vector<float>   *GJetAk04Eta;
  std::vector<float>   *GJetAk04Phi;
  std::vector<float>   *GJetAk04E;
  std::vector<float>   *GJetAk04ChFrac;
  std::vector<int>     *GJetAk04ConstCnt;
  std::vector<int>     *GJetAk04ConstId;
  std::vector<float>   *GJetAk04ConstPt;
  std::vector<float>   *GJetAk04ConstEta;
  std::vector<float>   *GJetAk04ConstPhi;
  std::vector<float>   *GJetAk04ConstE;
  std::vector<int>     *GPdfId1;
  std::vector<int>     *GPdfId2;
  std::vector<float>   *GPdfx1;
  std::vector<float>   *GPdfx2;
  std::vector<float>   *GPdfScale;
  Float_t         GBinningValue;
  Int_t           GNup;
  std::vector<float>   *MuPt;
  std::vector<float>   *MuEta;
  std::vector<float>   *MuPhi;
  std::vector<float>   *MuE;
  std::vector<unsigned int> *MuId;
  std::vector<unsigned int> *MuIdTight;
  std::vector<float>   *MuCh;
  std::vector<float>   *MuVtxZ;
  std::vector<float>   *MuDxy;
  std::vector<float>   *MuPfIso;
  std::vector<float>   *MuType;
  std::vector<float>   *MuIsoTkIsoAbs;
  std::vector<float>   *MuIsoTkIsoRel;
  std::vector<float>   *MuIsoCalAbs;
  std::vector<float>   *MuIsoCombRel;
  std::vector<float>   *MuTkNormChi2;
  std::vector<int>     *MuTkHitCnt;
  std::vector<int>     *MuMatchedStationCnt;
  std::vector<float>   *MuDz;
  std::vector<int>     *MuPixelHitCnt;
  std::vector<int>     *MuTkLayerCnt;
  std::vector<float>   *MuPfIsoChHad;
  std::vector<float>   *MuPfIsoNeutralHad;
  std::vector<float>   *MuPfIsoRawRel;
  std::vector<unsigned int> *MuHltMatch;
  std::vector<float>   *ElPt;
  std::vector<float>   *ElEta;
  std::vector<float>   *ElEtaSc;
  std::vector<float>   *ElPhi;
  std::vector<float>   *ElE;
  std::vector<unsigned int> *ElId;
  std::vector<float>   *ElCh;
  std::vector<float>   *ElMvaTrig;
  std::vector<float>   *ElMvaNonTrig;
  std::vector<float>   *ElMvaPresel;
  std::vector<float>   *ElDEtaTkScAtVtx;
  std::vector<float>   *ElDPhiTkScAtVtx;
  std::vector<float>   *ElHoE;
  std::vector<float>   *ElSigmaIetaIeta;
  std::vector<float>   *ElSigmaIetaIetaFull5x5;
  std::vector<float>   *ElEinvMinusPinv;
  std::vector<float>   *ElD0;
  std::vector<float>   *ElDz;
  std::vector<int>     *ElExpectedMissingInnerHitCnt;
  std::vector<int>     *ElPassConvVeto;
  std::vector<unsigned int> *ElHltMatch;
  std::vector<float>   *ElPfIsoChHad;
  std::vector<float>   *ElPfIsoNeutralHad;
  std::vector<float>   *ElPfIsoIso;
  std::vector<float>   *ElPfIsoPuChHad;
  std::vector<float>   *ElPfIsoRaw;
  std::vector<float>   *ElPfIsoDbeta;
  std::vector<float>   *ElAEff;
  std::vector<float>   *charged;
  std::vector<float>   *photon;
  std::vector<float>   *neutral;
  std::vector<float>   *charged_Tom;
  std::vector<float>   *photon_Tom;
  std::vector<float>   *neutral_Tom;
  std::vector<float>   *PhotPt;
  std::vector<float>   *PhotEta;
  std::vector<float>   *PhotPhi;
  std::vector<float>   *PhotScRawE;
  std::vector<float>   *PhotScEta;
  std::vector<float>   *PhotScPhi;
  std::vector<float>   *PhotIsoEcal;
  std::vector<float>   *PhotIsoHcal;
  std::vector<float>   *PhotIsoTk;
  std::vector<float>   *PhotPfIsoChHad;
  std::vector<float>   *PhotPfIsoNeutralHad;
  std::vector<float>   *PhotPfIsoPhot;
  std::vector<float>   *PhotPfIsoPuChHad;
  std::vector<float>   *PhotPfIsoEcalClus;
  std::vector<float>   *PhotPfIsoHcalClus;
  std::vector<float>   *PhotE3x3;
  std::vector<float>   *PhotE1x5;
  std::vector<float>   *PhotE2x5;
  std::vector<float>   *PhotE5x5;
  std::vector<float>   *PhotSigmaIetaIeta;
  std::vector<float>   *PhotEtaWidth;
  std::vector<float>   *PhotPhiWidth;
  std::vector<float>   *PhotHoE;
  std::vector<unsigned int> *PhotId;
  std::vector<bool>    *PhotHasPixelSeed;
  std::vector<float>   *JetAk04Pt;
  std::vector<float>   *JetAk04Eta;
  std::vector<float>   *JetAk04Phi;
  std::vector<float>   *JetAk04E;
  std::vector<float>   *JetAk04Id;
  std::vector<bool>    *JetAk04PuId;
  std::vector<float>   *JetAk04PuMva;
  std::vector<float>   *JetAk04RawPt;
  std::vector<float>   *JetAk04RawE;
  std::vector<float>   *JetAk04HfHadE;
  std::vector<float>   *JetAk04HfEmE;
  std::vector<float>   *JetAk04ChHadFrac;
  std::vector<float>   *JetAk04NeutralHadAndHfFrac;
  std::vector<float>   *JetAk04ChEmFrac;
  std::vector<float>   *JetAk04NeutralEmFrac;
  std::vector<float>   *JetAk04ChMult;
  std::vector<float>   *JetAk04ConstCnt;
  std::vector<float>   *JetAk04JetBeta;
  std::vector<float>   *JetAk04JetBetaClassic;
  std::vector<float>   *JetAk04JetBetaStar;
  std::vector<float>   *JetAk04JetBetaStarClassic;
  std::vector<float>   *JetAk04BTagCsv;
  std::vector<float>   *JetAk04BTagCsvV1;
  std::vector<float>   *JetAk04BTagCsvSLV1;
  std::vector<float>   *JetAk04BDiscCisvV2;
  std::vector<float>   *JetAk04BDiscJp;
  std::vector<float>   *JetAk04BDiscBjp;
  std::vector<float>   *JetAk04BDiscTche;
  std::vector<float>   *JetAk04BDiscTchp;
  std::vector<float>   *JetAk04BDiscSsvhe;
  std::vector<float>   *JetAk04BDiscSsvhp;
  std::vector<float>   *JetAk04PartFlav;
  std::vector<float>   *JetAk04JecUncUp;
  std::vector<float>   *JetAk04JecUncDwn;
  std::vector<int>     *JetAk04ConstId;
  std::vector<float>   *JetAk04ConstPt;
  std::vector<float>   *JetAk04ConstEta;
  std::vector<float>   *JetAk04ConstPhi;
  std::vector<float>   *JetAk04ConstE;
  std::vector<int>     *JetAk04GenJet;

  // List of branches
  TBranch        *b_EvtIsRealData;   //!
  TBranch        *b_EvtNum;   //!
  TBranch        *b_EvtRunNum;   //!
  TBranch        *b_EvtLumiNum;   //!
  TBranch        *b_EvtBxNum;   //!
  TBranch        *b_EvtVtxCnt;   //!
  TBranch        *b_EvtPuCnt;   //!
  TBranch        *b_EvtPuCntTruth;   //!
  TBranch        *b_EvtWeights;   //!
  TBranch        *b_TrigHlt;   //!
  TBranch        *b_TrigHltMu;   //!
  TBranch        *b_TrigHltDiMu;   //!
  TBranch        *b_TrigHltEl;   //!
  TBranch        *b_TrigHltDiEl;   //!
  TBranch        *b_TrigHltElMu;   //!
  TBranch        *b_METPt;   //!
  TBranch        *b_METPx;   //!
  TBranch        *b_METPy;   //!
  TBranch        *b_METPz;   //!
  TBranch        *b_METE;   //!
  TBranch        *b_METsigx2;   //!
  TBranch        *b_METsigxy;   //!
  TBranch        *b_METsigy2;   //!
  TBranch        *b_METsig;   //!
  TBranch        *b_GLepDr01Pt;   //!
  TBranch        *b_GLepDr01Eta;   //!
  TBranch        *b_GLepDr01Phi;   //!
  TBranch        *b_GLepDr01E;   //!
  TBranch        *b_GLepDr01Id;   //!
  TBranch        *b_GLepDr01St;   //!
  TBranch        *b_GLepDr01MomId;   //!
  TBranch        *b_GLepBarePt;   //!
  TBranch        *b_GLepBareEta;   //!
  TBranch        *b_GLepBarePhi;   //!
  TBranch        *b_GLepBareE;   //!
  TBranch        *b_GLepBareId;   //!
  TBranch        *b_GLepBareSt;   //!
  TBranch        *b_GLepBareMomId;   //!
  TBranch        *b_GLepSt3Pt;   //!
  TBranch        *b_GLepSt3Eta;   //!
  TBranch        *b_GLepSt3Phi;   //!
  TBranch        *b_GLepSt3E;   //!
  TBranch        *b_GLepSt3Id;   //!
  TBranch        *b_GLepSt3St;   //!
  TBranch        *b_GLepSt3Mother0Id;   //!
  TBranch        *b_GLepSt3MotherCnt;   //!
  TBranch        *b_GPhotPt;   //!
  TBranch        *b_GPhotEta;   //!
  TBranch        *b_GPhotPhi;   //!
  TBranch        *b_GPhotE;   //!
  TBranch        *b_GPhotMotherId;   //!
  TBranch        *b_GPhotSt;   //!
  TBranch        *b_GPhotIsoEDR03;   //!
  TBranch        *b_GPhotIsoEDR04;   //!
  TBranch        *b_GPhotIsoEDR05;   //!
  TBranch        *b_GPhotIsoSumPtDR03;   //!
  TBranch        *b_GPhotIsoSumPtDR04;   //!
  TBranch        *b_GPhotIsoSumPtDR05;   //!
  TBranch        *b_GLepClosePhotPt;   //!
  TBranch        *b_GLepClosePhotEta;   //!
  TBranch        *b_GLepClosePhotPhi;   //!
  TBranch        *b_GLepClosePhotE;   //!
  TBranch        *b_GLepClosePhotId;   //!
  TBranch        *b_GLepClosePhotMother0Id;   //!
  TBranch        *b_GLepClosePhotMotherCnt;   //!
  TBranch        *b_GLepClosePhotSt;   //!
  TBranch        *b_GJetAk04Pt;   //!
  TBranch        *b_GJetAk04Eta;   //!
  TBranch        *b_GJetAk04Phi;   //!
  TBranch        *b_GJetAk04E;   //!
  TBranch        *b_GJetAk04ChFrac;   //!
  TBranch        *b_GJetAk04ConstCnt;   //!
  TBranch        *b_GJetAk04ConstId;   //!
  TBranch        *b_GJetAk04ConstPt;   //!
  TBranch        *b_GJetAk04ConstEta;   //!
  TBranch        *b_GJetAk04ConstPhi;   //!
  TBranch        *b_GJetAk04ConstE;   //!
  TBranch        *b_GPdfId1;   //!
  TBranch        *b_GPdfId2;   //!
  TBranch        *b_GPdfx1;   //!
  TBranch        *b_GPdfx2;   //!
  TBranch        *b_GPdfScale;   //!
  TBranch        *b_GBinningValue;   //!
  TBranch        *b_GNup;   //!
  TBranch        *b_MuPt;   //!
  TBranch        *b_MuEta;   //!
  TBranch        *b_MuPhi;   //!
  TBranch        *b_MuE;   //!
  TBranch        *b_MuId;   //!
  TBranch        *b_MuIdTight;   //!
  TBranch        *b_MuCh;   //!
  TBranch        *b_MuVtxZ;   //!
  TBranch        *b_MuDxy;   //!
  TBranch        *b_MuPfIso;   //!
  TBranch        *b_MuType;   //!
  TBranch        *b_MuIsoTkIsoAbs;   //!
  TBranch        *b_MuIsoTkIsoRel;   //!
  TBranch        *b_MuIsoCalAbs;   //!
  TBranch        *b_MuIsoCombRel;   //!
  TBranch        *b_MuTkNormChi2;   //!
  TBranch        *b_MuTkHitCnt;   //!
  TBranch        *b_MuMatchedStationCnt;   //!
  TBranch        *b_MuDz;   //!
  TBranch        *b_MuPixelHitCnt;   //!
  TBranch        *b_MuTkLayerCnt;   //!
  TBranch        *b_MuPfIsoChHad;   //!
  TBranch        *b_MuPfIsoNeutralHad;   //!
  TBranch        *b_MuPfIsoRawRel;   //!
  TBranch        *b_MuHltMatch;   //!
  TBranch        *b_ElPt;   //!
  TBranch        *b_ElEta;   //!
  TBranch        *b_ElEtaSc;   //!
  TBranch        *b_ElPhi;   //!
  TBranch        *b_ElE;   //!
  TBranch        *b_ElId;   //!
  TBranch        *b_ElCh;   //!
  TBranch        *b_ElMvaTrig;   //!
  TBranch        *b_ElMvaNonTrig;   //!
  TBranch        *b_ElMvaPresel;   //!
  TBranch        *b_ElDEtaTkScAtVtx;   //!
  TBranch        *b_ElDPhiTkScAtVtx;   //!
  TBranch        *b_ElHoE;   //!
  TBranch        *b_ElSigmaIetaIeta;   //!
  TBranch        *b_ElSigmaIetaIetaFull5x5;   //!
  TBranch        *b_ElEinvMinusPinv;   //!
  TBranch        *b_ElD0;   //!
  TBranch        *b_ElDz;   //!
  TBranch        *b_ElExpectedMissingInnerHitCnt;   //!
  TBranch        *b_ElPassConvVeto;   //!
  TBranch        *b_ElHltMatch;   //!
  TBranch        *b_ElPfIsoChHad;   //!
  TBranch        *b_ElPfIsoNeutralHad;   //!
  TBranch        *b_ElPfIsoIso;   //!
  TBranch        *b_ElPfIsoPuChHad;   //!
  TBranch        *b_ElPfIsoRaw;   //!
  TBranch        *b_ElPfIsoDbeta;   //!
  TBranch        *b_ElAEff;   //!
  TBranch        *b_charged;   //!
  TBranch        *b_photon;   //!
  TBranch        *b_neutral;   //!
  TBranch        *b_charged_Tom;   //!
  TBranch        *b_photon_Tom;   //!
  TBranch        *b_neutral_Tom;   //!
  TBranch        *b_PhotPt;   //!
  TBranch        *b_PhotEta;   //!
  TBranch        *b_PhotPhi;   //!
  TBranch        *b_PhotScRawE;   //!
  TBranch        *b_PhotScEta;   //!
  TBranch        *b_PhotScPhi;   //!
  TBranch        *b_PhotIsoEcal;   //!
  TBranch        *b_PhotIsoHcal;   //!
  TBranch        *b_PhotIsoTk;   //!
  TBranch        *b_PhotPfIsoChHad;   //!
  TBranch        *b_PhotPfIsoNeutralHad;   //!
  TBranch        *b_PhotPfIsoPhot;   //!
  TBranch        *b_PhotPfIsoPuChHad;   //!
  TBranch        *b_PhotPfIsoEcalClus;   //!
  TBranch        *b_PhotPfIsoHcalClus;   //!
  TBranch        *b_PhotE3x3;   //!
  TBranch        *b_PhotE1x5;   //!
  TBranch        *b_PhotE2x5;   //!
  TBranch        *b_PhotE5x5;   //!
  TBranch        *b_PhotSigmaIetaIeta;   //!
  TBranch        *b_PhotEtaWidth;   //!
  TBranch        *b_PhotPhiWidth;   //!
  TBranch        *b_PhotHoE;   //!
  TBranch        *b_PhotId;   //!
  TBranch        *b_PhotHasPixelSeed;   //!
  TBranch        *b_JetAk04Pt;   //!
  TBranch        *b_JetAk04Eta;   //!
  TBranch        *b_JetAk04Phi;   //!
  TBranch        *b_JetAk04E;   //!
  TBranch        *b_JetAk04Id;   //!
  TBranch        *b_JetAk04PuId;   //!
  TBranch        *b_JetAk04PuMva;   //!
  TBranch        *b_JetAk04RawPt;   //!
  TBranch        *b_JetAk04RawE;   //!
  TBranch        *b_JetAk04HfHadE;   //!
  TBranch        *b_JetAk04HfEmE;   //!
  TBranch        *b_JetAk04ChHadFrac;   //!
  TBranch        *b_JetAk04NeutralHadAndHfFrac;   //!
  TBranch        *b_JetAk04ChEmFrac;   //!
  TBranch        *b_JetAk04NeutralEmFrac;   //!
  TBranch        *b_JetAk04ChMult;   //!
  TBranch        *b_JetAk04ConstCnt;   //!
  TBranch        *b_JetAk04JetBeta;   //!
  TBranch        *b_JetAk04JetBetaClassic;   //!
  TBranch        *b_JetAk04JetBetaStar;   //!
  TBranch        *b_JetAk04JetBetaStarClassic;   //!
  TBranch        *b_JetAk04BTagCsv;   //!
  TBranch        *b_JetAk04BTagCsvV1;   //!
  TBranch        *b_JetAk04BTagCsvSLV1;   //!
  TBranch        *b_JetAk04BDiscCisvV2;   //!
  TBranch        *b_JetAk04BDiscJp;   //!
  TBranch        *b_JetAk04BDiscBjp;   //!
  TBranch        *b_JetAk04BDiscTche;   //!
  TBranch        *b_JetAk04BDiscTchp;   //!
  TBranch        *b_JetAk04BDiscSsvhe;   //!
  TBranch        *b_JetAk04BDiscSsvhp;   //!
  TBranch        *b_JetAk04PartFlav;   //!
  TBranch        *b_JetAk04JecUncUp;   //!
  TBranch        *b_JetAk04JecUncDwn;   //!
  TBranch        *b_JetAk04ConstId;   //!
  TBranch        *b_JetAk04ConstPt;   //!
  TBranch        *b_JetAk04ConstEta;   //!
  TBranch        *b_JetAk04ConstPhi;   //!
  TBranch        *b_JetAk04ConstE;   //!
  TBranch        *b_JetAk04GenJet;   //!

 private:
  const char *fname;

  ClassDef(EventTree,0)
};
#endif
