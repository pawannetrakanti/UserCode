//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Dec 18 07:26:09 2015 by ROOT version 6.05/03
// from TTree EventTree/ EventTree
// found on file: root://eoscms//eos/cms/store/group/phys_heavyions/azsigmon/HiForestAODpp5TeV/HiForestAOD_withTupel_pp_MC_Z30mumuJet_v1.root
//////////////////////////////////////////////////////////
#include "commonSetup.h"
#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <vector>


#ifdef __MAKECINT__
#pragma link C++ class vector<double>+;
#pragma link C++ class vector<float>+;
#pragma link C++ class vector<int>+;
#pragma link C++ class vector<unsigned int>+;
#endif


class Tuples {
public :
   Tuples(){};
   ~Tuples(){};

   // Declaration of leaf types
   Int_t           EvtIsRealData;
   UInt_t          EvtNum;
   UInt_t          EvtRunNum;
   Int_t           EvtLumiNum;
   Int_t           EvtBxNum;
   Int_t           EvtVtxCnt;
   Int_t           EvtPuCnt;
   Int_t           EvtPuCntTruth;

   UInt_t          TrigHlt;
   ULong64_t       TrigHltMu;
   ULong64_t       TrigHltDiMu;
   ULong64_t       TrigHltEl;
   ULong64_t       TrigHltDiEl;
   ULong64_t       TrigHltElMu;
   Float_t         GBinningValue;
   Int_t           GNup;

   std::vector<double>  *EvtWeights=0;
   std::vector<float>   *METPt=0;
   std::vector<float>   *METPx=0;
   std::vector<float>   *METPy=0;
   std::vector<float>   *METPz=0;
   std::vector<float>   *METE=0;
   std::vector<float>   *METsigx2=0;
   std::vector<float>   *METsigxy=0;
   std::vector<float>   *METsigy2=0;
   std::vector<float>   *METsig=0;
   std::vector<float>   *GLepDr01Pt=0;
   std::vector<float>   *GLepDr01Eta=0;
   std::vector<float>   *GLepDr01Phi=0;
   std::vector<float>   *GLepDr01E=0;
   std::vector<int>     *GLepDr01Id=0;
   std::vector<int>     *GLepDr01St=0;
   std::vector<int>     *GLepDr01MomId=0;
   std::vector<float>   *GLepBarePt=0;
   std::vector<float>   *GLepBareEta=0;
   std::vector<float>   *GLepBarePhi=0;
   std::vector<float>   *GLepBareE=0;
   std::vector<int>     *GLepBareId=0;
   std::vector<int>     *GLepBareSt=0;
   std::vector<int>     *GLepBareMomId=0;
   std::vector<float>   *GLepSt3Pt=0;
   std::vector<float>   *GLepSt3Eta=0;
   std::vector<float>   *GLepSt3Phi=0;
   std::vector<float>   *GLepSt3E=0;
   std::vector<int>     *GLepSt3Id=0;
   std::vector<int>     *GLepSt3St=0;
   std::vector<int>     *GLepSt3Mother0Id=0;
   std::vector<int>     *GLepSt3MotherCnt=0;
   std::vector<float>   *GPhotPt=0;
   std::vector<float>   *GPhotEta=0;
   std::vector<float>   *GPhotPhi=0;
   std::vector<float>   *GPhotE=0;
   std::vector<int>     *GPhotMotherId=0;
   std::vector<int>     *GPhotSt=0;
   std::vector<float>   *GPhotIsoEDR03=0;
   std::vector<float>   *GPhotIsoEDR04=0;
   std::vector<float>   *GPhotIsoEDR05=0;
   std::vector<float>   *GPhotIsoSumPtDR03=0;
   std::vector<float>   *GPhotIsoSumPtDR04=0;
   std::vector<float>   *GPhotIsoSumPtDR05=0;
   std::vector<float>   *GLepClosePhotPt=0;
   std::vector<float>   *GLepClosePhotEta=0;
   std::vector<float>   *GLepClosePhotPhi=0;
   std::vector<float>   *GLepClosePhotE=0;
   std::vector<int>     *GLepClosePhotId=0;
   std::vector<int>     *GLepClosePhotMother0Id=0;
   std::vector<int>     *GLepClosePhotMotherCnt=0;
   std::vector<int>     *GLepClosePhotSt=0;
   std::vector<float>   *GJetAk04Pt=0;
   std::vector<float>   *GJetAk04Eta=0;
   std::vector<float>   *GJetAk04Phi=0;
   std::vector<float>   *GJetAk04E=0;
   std::vector<float>   *GJetAk04ChFrac=0;
   std::vector<int>     *GJetAk04ConstCnt=0;
   std::vector<int>     *GJetAk04ConstId=0;
   std::vector<float>   *GJetAk04ConstPt=0;
   std::vector<float>   *GJetAk04ConstEta=0;
   std::vector<float>   *GJetAk04ConstPhi=0;
   std::vector<float>   *GJetAk04ConstE=0;
   std::vector<int>     *GPdfId1=0;
   std::vector<int>     *GPdfId2=0;
   std::vector<float>   *GPdfx1=0;
   std::vector<float>   *GPdfx2=0;
   std::vector<float>   *GPdfScale=0;

   std::vector<float>   *MuPt=0;
   std::vector<float>   *MuEta=0;
   std::vector<float>   *MuPhi=0;
   std::vector<float>   *MuE=0;
   std::vector<unsigned int> *MuId=0;
   std::vector<unsigned int> *MuIdTight=0;
   std::vector<float>   *MuCh=0;
   std::vector<float>   *MuVtxZ=0;
   std::vector<float>   *MuDxy=0;
   std::vector<float>   *MuPfIso=0;
   std::vector<float>   *MuType=0;
   std::vector<float>   *MuIsoTkIsoAbs=0;
   std::vector<float>   *MuIsoTkIsoRel=0;
   std::vector<float>   *MuIsoCalAbs=0;
   std::vector<float>   *MuIsoCombRel=0;
   std::vector<float>   *MuTkNormChi2=0;
   std::vector<int>     *MuTkHitCnt=0;
   std::vector<int>     *MuMatchedStationCnt=0;
   std::vector<float>   *MuDz=0;
   std::vector<int>     *MuPixelHitCnt=0;
   std::vector<int>     *MuTkLayerCnt=0;
   std::vector<float>   *MuPfIsoChHad=0;
   std::vector<float>   *MuPfIsoNeutralHad=0;
   std::vector<float>   *MuPfIsoRawRel=0;
   std::vector<unsigned int> *MuHltMatch=0;
   std::vector<float>   *ElPt=0;
   std::vector<float>   *ElEta=0;
   std::vector<float>   *ElEtaSc=0;
   std::vector<float>   *ElPhi=0;
   std::vector<float>   *ElE=0;
   std::vector<unsigned int> *ElId=0;
   std::vector<float>   *ElCh=0;
   std::vector<float>   *ElMvaTrig=0;
   std::vector<float>   *ElMvaNonTrig=0;
   std::vector<float>   *ElMvaPresel=0;
   std::vector<float>   *ElDEtaTkScAtVtx=0;
   std::vector<float>   *ElDPhiTkScAtVtx=0;
   std::vector<float>   *ElHoE=0;
   std::vector<float>   *ElSigmaIetaIeta=0;
   std::vector<float>   *ElSigmaIetaIetaFull5x5=0;
   std::vector<float>   *ElEinvMinusPinv=0;
   std::vector<float>   *ElD0=0;
   std::vector<float>   *ElDz=0;
   std::vector<int>     *ElExpectedMissingInnerHitCnt=0;
   std::vector<int>     *ElPassConvVeto=0;
   std::vector<unsigned int> *ElHltMatch=0;
   std::vector<float>   *ElPfIsoChHad=0;
   std::vector<float>   *ElPfIsoNeutralHad=0;
   std::vector<float>   *ElPfIsoIso=0;
   std::vector<float>   *ElPfIsoPuChHad=0;
   std::vector<float>   *ElPfIsoRaw=0;
   std::vector<float>   *ElPfIsoDbeta=0;
   std::vector<float>   *ElAEff=0;
   std::vector<float>   *charged=0;
   std::vector<float>   *photon=0;
   std::vector<float>   *neutral=0;
   std::vector<float>   *charged_Tom=0;
   std::vector<float>   *photon_Tom=0;
   std::vector<float>   *neutral_Tom=0;
   std::vector<float>   *PhotPt=0;
   std::vector<float>   *PhotEta=0;
   std::vector<float>   *PhotPhi=0;
   std::vector<float>   *PhotScRawE=0;
   std::vector<float>   *PhotScEta=0;
   std::vector<float>   *PhotScPhi=0;
   std::vector<float>   *PhotIsoEcal=0;
   std::vector<float>   *PhotIsoHcal=0;
   std::vector<float>   *PhotIsoTk=0;
   std::vector<float>   *PhotPfIsoChHad=0;
   std::vector<float>   *PhotPfIsoNeutralHad=0;
   std::vector<float>   *PhotPfIsoPhot=0;
   std::vector<float>   *PhotPfIsoPuChHad=0;
   std::vector<float>   *PhotPfIsoEcalClus=0;
   std::vector<float>   *PhotPfIsoHcalClus=0;
   std::vector<float>   *PhotE3x3=0;
   std::vector<float>   *PhotE1x5=0;
   std::vector<float>   *PhotE2x5=0;
   std::vector<float>   *PhotE5x5=0;
   std::vector<float>   *PhotSigmaIetaIeta=0;
   std::vector<float>   *PhotEtaWidth=0;
   std::vector<float>   *PhotPhiWidth=0;
   std::vector<float>   *PhotHoE=0;
   std::vector<unsigned int> *PhotId=0;
   std::vector<bool>    *PhotHasPixelSeed=0;
   std::vector<float>   *JetAk04Pt=0;
   std::vector<float>   *JetAk04Eta=0;
   std::vector<float>   *JetAk04Phi=0;
   std::vector<float>   *JetAk04E=0;
   std::vector<float>   *JetAk04Id=0;
   std::vector<bool>    *JetAk04PuId=0;
   std::vector<float>   *JetAk04PuMva=0;
   std::vector<float>   *JetAk04RawPt=0;
   std::vector<float>   *JetAk04RawE=0;
   std::vector<float>   *JetAk04HfHadE=0;
   std::vector<float>   *JetAk04HfEmE=0;
   std::vector<float>   *JetAk04ChHadFrac=0;
   std::vector<float>   *JetAk04NeutralHadAndHfFrac=0;
   std::vector<float>   *JetAk04ChEmFrac=0;
   std::vector<float>   *JetAk04NeutralEmFrac=0;
   std::vector<float>   *JetAk04ChMult=0;
   std::vector<float>   *JetAk04ConstCnt=0;
   std::vector<float>   *JetAk04JetBeta=0;
   std::vector<float>   *JetAk04JetBetaClassic=0;
   std::vector<float>   *JetAk04JetBetaStar=0;
   std::vector<float>   *JetAk04JetBetaStarClassic=0;
   std::vector<float>   *JetAk04BTagCsv=0;
   std::vector<float>   *JetAk04BTagCsvV1=0;
   std::vector<float>   *JetAk04BTagCsvSLV1=0;
   std::vector<float>   *JetAk04BDiscCisvV2=0;
   std::vector<float>   *JetAk04BDiscJp=0;
   std::vector<float>   *JetAk04BDiscBjp=0;
   std::vector<float>   *JetAk04BDiscTche=0;
   std::vector<float>   *JetAk04BDiscTchp=0;
   std::vector<float>   *JetAk04BDiscSsvhe=0;
   std::vector<float>   *JetAk04BDiscSsvhp=0;
   std::vector<float>   *JetAk04PartFlav=0;
   std::vector<float>   *JetAk04JecUncUp=0;
   std::vector<float>   *JetAk04JecUncDwn=0;
   std::vector<int>     *JetAk04ConstId=0;
   std::vector<float>   *JetAk04ConstPt=0;
   std::vector<float>   *JetAk04ConstEta=0;
   std::vector<float>   *JetAk04ConstPhi=0;
   std::vector<float>   *JetAk04ConstE=0;
   std::vector<int>     *JetAk04GenJet=0;

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
};
void setupTupleTree(TTree *t,Tuples *tTuples,bool doCheck = 0)
{

  // Set branch addresses and branch pointers
  /* t->SetBranchAddress("EvtIsRealData", &tTuples.EvtIsRealData, &tTuples.b_EvtIsRealData); */
  /* t->SetBranchAddress("EvtNum", &tTuples.EvtNum, &tTuples.b_EvtNum); */
  /* t->SetBranchAddress("EvtRunNum", &tTuples.EvtRunNum, &tTuples.b_EvtRunNum); */
  /* t->SetBranchAddress("EvtLumiNum", &tTuples.EvtLumiNum, &tTuples.b_EvtLumiNum); */
  /* t->SetBranchAddress("EvtBxNum", &tTuples.EvtBxNum, &tTuples.b_EvtBxNum); */
  /* t->SetBranchAddress("EvtVtxCnt", &tTuples.EvtVtxCnt, &tTuples.b_EvtVtxCnt); */
  /* t->SetBranchAddress("EvtPuCnt", &tTuples.EvtPuCnt, &tTuples.b_EvtPuCnt); */
  /* t->SetBranchAddress("EvtPuCntTruth", &tTuples.EvtPuCntTruth, &tTuples.b_EvtPuCntTruth); */
  /* t->SetBranchAddress("GBinningValue", &tTuples.GBinningValue, &tTuples.b_GBinningValue); */
  /* t->SetBranchAddress("GNup", &tTuples.GNup, &tTuples.b_GNup); */
  /* t->SetBranchAddress("TrigHlt", &tTuples.TrigHlt, &tTuples.b_TrigHlt); */
  /* t->SetBranchAddress("TrigHltMu", &tTuples.TrigHltMu, &tTuples.b_TrigHltMu); */
  /* t->SetBranchAddress("TrigHltDiMu", &tTuples.TrigHltDiMu, &tTuples.b_TrigHltDiMu); */
  /* t->SetBranchAddress("TrigHltEl", &tTuples.TrigHltEl, &tTuples.b_TrigHltEl); */
  /* t->SetBranchAddress("TrigHltDiEl", &tTuples.TrigHltDiEl, &tTuples.b_TrigHltDiEl); */
  /* t->SetBranchAddress("TrigHltElMu", &tTuples.TrigHltElMu, &tTuples.b_TrigHltElMu); */

  t->SetBranchAddress("EvtWeights", tTuples->EvtWeights, &tTuples->b_EvtWeights);
  /* t->SetBranchAddress("METPt", tTuples.METPt, &tTuples.b_METPt); */
  /* t->SetBranchAddress("METPx", tTuples.METPx, &tTuples.b_METPx); */
  /* t->SetBranchAddress("METPy", tTuples.METPy, &tTuples.b_METPy); */
  /* t->SetBranchAddress("METPz", tTuples.METPz, &tTuples.b_METPz); */
  /* t->SetBranchAddress("METE", tTuples.METE, &tTuples.b_METE); */
  /* t->SetBranchAddress("METsigx2", tTuples.METsigx2, &tTuples.b_METsigx2); */
  /* t->SetBranchAddress("METsigxy", tTuples.METsigxy, &tTuples.b_METsigxy); */
  /* t->SetBranchAddress("METsigy2", tTuples.METsigy2, &tTuples.b_METsigy2); */
  /* t->SetBranchAddress("METsig", tTuples.METsig, &tTuples.b_METsig); */
  /* t->SetBranchAddress("GLepDr01Pt", tTuples.GLepDr01Pt, &tTuples.b_GLepDr01Pt); */
  /* t->SetBranchAddress("GLepDr01Eta", tTuples.GLepDr01Eta, &tTuples.b_GLepDr01Eta); */
  /* t->SetBranchAddress("GLepDr01Phi", tTuples.GLepDr01Phi, &tTuples.b_GLepDr01Phi); */
  /* t->SetBranchAddress("GLepDr01E", tTuples.GLepDr01E, &tTuples.b_GLepDr01E); */
  /* t->SetBranchAddress("GLepDr01Id", tTuples.GLepDr01Id, &tTuples.b_GLepDr01Id); */
  /* t->SetBranchAddress("GLepDr01St", tTuples.GLepDr01St, &tTuples.b_GLepDr01St); */
  /* t->SetBranchAddress("GLepDr01MomId", tTuples.GLepDr01MomId, &tTuples.b_GLepDr01MomId); */
  /* t->SetBranchAddress("GLepBarePt", tTuples.GLepBarePt, &tTuples.b_GLepBarePt); */
  /* t->SetBranchAddress("GLepBareEta", tTuples.GLepBareEta, &tTuples.b_GLepBareEta); */
  /* t->SetBranchAddress("GLepBarePhi", tTuples.GLepBarePhi, &tTuples.b_GLepBarePhi); */
  /* t->SetBranchAddress("GLepBareE", tTuples.GLepBareE, &tTuples.b_GLepBareE); */
  /* t->SetBranchAddress("GLepBareId", tTuples.GLepBareId, &tTuples.b_GLepBareId); */
  /* t->SetBranchAddress("GLepBareSt", tTuples.GLepBareSt, &tTuples.b_GLepBareSt); */
  /* t->SetBranchAddress("GLepBareMomId", tTuples.GLepBareMomId, &tTuples.b_GLepBareMomId); */
  /* t->SetBranchAddress("GLepSt3Pt", tTuples.GLepSt3Pt, &tTuples.b_GLepSt3Pt); */
  /* t->SetBranchAddress("GLepSt3Eta", tTuples.GLepSt3Eta, &tTuples.b_GLepSt3Eta); */
  /* t->SetBranchAddress("GLepSt3Phi", tTuples.GLepSt3Phi, &tTuples.b_GLepSt3Phi); */
  /* t->SetBranchAddress("GLepSt3E", tTuples.GLepSt3E, &tTuples.b_GLepSt3E); */
  /* t->SetBranchAddress("GLepSt3Id", tTuples.GLepSt3Id, &tTuples.b_GLepSt3Id); */
  /* t->SetBranchAddress("GLepSt3St", tTuples.GLepSt3St, &tTuples.b_GLepSt3St); */
  /* t->SetBranchAddress("GLepSt3Mother0Id", tTuples.GLepSt3Mother0Id, &tTuples.b_GLepSt3Mother0Id); */
  /* t->SetBranchAddress("GLepSt3MotherCnt", tTuples.GLepSt3MotherCnt, &tTuples.b_GLepSt3MotherCnt); */
  /* t->SetBranchAddress("GPhotPt", tTuples.GPhotPt, &tTuples.b_GPhotPt); */
  /* t->SetBranchAddress("GPhotEta", tTuples.GPhotEta, &tTuples.b_GPhotEta); */
  /* t->SetBranchAddress("GPhotPhi", tTuples.GPhotPhi, &tTuples.b_GPhotPhi); */
  /* t->SetBranchAddress("GPhotE", tTuples.GPhotE, &tTuples.b_GPhotE); */
  /* t->SetBranchAddress("GPhotMotherId", tTuples.GPhotMotherId, &tTuples.b_GPhotMotherId); */
  /* t->SetBranchAddress("GPhotSt", tTuples.GPhotSt, &tTuples.b_GPhotSt); */
  /* t->SetBranchAddress("GPhotIsoEDR03", tTuples.GPhotIsoEDR03, &tTuples.b_GPhotIsoEDR03); */
  /* t->SetBranchAddress("GPhotIsoEDR04", tTuples.GPhotIsoEDR04, &tTuples.b_GPhotIsoEDR04); */
  /* t->SetBranchAddress("GPhotIsoEDR05", tTuples.GPhotIsoEDR05, &tTuples.b_GPhotIsoEDR05); */
  /* t->SetBranchAddress("GPhotIsoSumPtDR03", tTuples.GPhotIsoSumPtDR03, &tTuples.b_GPhotIsoSumPtDR03); */
  /* t->SetBranchAddress("GPhotIsoSumPtDR04", tTuples.GPhotIsoSumPtDR04, &tTuples.b_GPhotIsoSumPtDR04); */
  /* t->SetBranchAddress("GPhotIsoSumPtDR05", tTuples.GPhotIsoSumPtDR05, &tTuples.b_GPhotIsoSumPtDR05); */
  /* t->SetBranchAddress("GLepClosePhotPt", tTuples.GLepClosePhotPt, &tTuples.b_GLepClosePhotPt); */
  /* t->SetBranchAddress("GLepClosePhotEta", tTuples.GLepClosePhotEta, &tTuples.b_GLepClosePhotEta); */
  /* t->SetBranchAddress("GLepClosePhotPhi", tTuples.GLepClosePhotPhi, &tTuples.b_GLepClosePhotPhi); */
  /* t->SetBranchAddress("GLepClosePhotE", tTuples.GLepClosePhotE, &tTuples.b_GLepClosePhotE); */
  /* t->SetBranchAddress("GLepClosePhotId", tTuples.GLepClosePhotId, &tTuples.b_GLepClosePhotId); */
  /* t->SetBranchAddress("GLepClosePhotMother0Id", tTuples.GLepClosePhotMother0Id, &tTuples.b_GLepClosePhotMother0Id); */
  /* t->SetBranchAddress("GLepClosePhotMotherCnt", tTuples.GLepClosePhotMotherCnt, &tTuples.b_GLepClosePhotMotherCnt); */
  /* t->SetBranchAddress("GLepClosePhotSt", tTuples.GLepClosePhotSt, &tTuples.b_GLepClosePhotSt); */
  /* t->SetBranchAddress("GJetAk04Pt", tTuples.GJetAk04Pt, &tTuples.b_GJetAk04Pt); */
  /* t->SetBranchAddress("GJetAk04Eta", tTuples.GJetAk04Eta, &tTuples.b_GJetAk04Eta); */
  /* t->SetBranchAddress("GJetAk04Phi", tTuples.GJetAk04Phi, &tTuples.b_GJetAk04Phi); */
  /* t->SetBranchAddress("GJetAk04E", tTuples.GJetAk04E, &tTuples.b_GJetAk04E); */
  /* t->SetBranchAddress("GJetAk04ChFrac", tTuples.GJetAk04ChFrac, &tTuples.b_GJetAk04ChFrac); */
  /* t->SetBranchAddress("GJetAk04ConstCnt", tTuples.GJetAk04ConstCnt, &tTuples.b_GJetAk04ConstCnt); */
  /* t->SetBranchAddress("GJetAk04ConstId", tTuples.GJetAk04ConstId, &tTuples.b_GJetAk04ConstId); */
  /* t->SetBranchAddress("GJetAk04ConstPt", tTuples.GJetAk04ConstPt, &tTuples.b_GJetAk04ConstPt); */
  /* t->SetBranchAddress("GJetAk04ConstEta", tTuples.GJetAk04ConstEta, &tTuples.b_GJetAk04ConstEta); */
  /* t->SetBranchAddress("GJetAk04ConstPhi", tTuples.GJetAk04ConstPhi, &tTuples.b_GJetAk04ConstPhi); */
  /* t->SetBranchAddress("GJetAk04ConstE", tTuples.GJetAk04ConstE, &tTuples.b_GJetAk04ConstE); */
  /* t->SetBranchAddress("GPdfId1", tTuples.GPdfId1, &tTuples.b_GPdfId1); */
  /* t->SetBranchAddress("GPdfId2", tTuples.GPdfId2, &tTuples.b_GPdfId2); */
  /* t->SetBranchAddress("GPdfx1", tTuples.GPdfx1, &tTuples.b_GPdfx1); */
  /* t->SetBranchAddress("GPdfx2", tTuples.GPdfx2, &tTuples.b_GPdfx2); */
  /* t->SetBranchAddress("GPdfScale", tTuples.GPdfScale, &tTuples.b_GPdfScale); */
  /* t->SetBranchAddress("MuPt", tTuples.MuPt, &tTuples.b_MuPt); */
  /* t->SetBranchAddress("MuEta", tTuples.MuEta, &tTuples.b_MuEta); */
  /* t->SetBranchAddress("MuPhi", tTuples.MuPhi, &tTuples.b_MuPhi); */
  /* t->SetBranchAddress("MuE", tTuples.MuE, &tTuples.b_MuE); */
  /* t->SetBranchAddress("MuId", tTuples.MuId, &tTuples.b_MuId); */
  /* t->SetBranchAddress("MuIdTight", tTuples.MuIdTight, &tTuples.b_MuIdTight); */
  /* t->SetBranchAddress("MuCh", tTuples.MuCh, &tTuples.b_MuCh); */
  /* t->SetBranchAddress("MuVtxZ", tTuples.MuVtxZ, &tTuples.b_MuVtxZ); */
  /* t->SetBranchAddress("MuDxy", tTuples.MuDxy, &tTuples.b_MuDxy); */
  /* t->SetBranchAddress("MuPfIso", tTuples.MuPfIso, &tTuples.b_MuPfIso); */
  /* t->SetBranchAddress("MuType", tTuples.MuType, &tTuples.b_MuType); */
  /* t->SetBranchAddress("MuIsoTkIsoAbs", tTuples.MuIsoTkIsoAbs, &tTuples.b_MuIsoTkIsoAbs); */
  /* t->SetBranchAddress("MuIsoTkIsoRel", tTuples.MuIsoTkIsoRel, &tTuples.b_MuIsoTkIsoRel); */
  /* t->SetBranchAddress("MuIsoCalAbs", tTuples.MuIsoCalAbs, &tTuples.b_MuIsoCalAbs); */
  /* t->SetBranchAddress("MuIsoCombRel", tTuples.MuIsoCombRel, &tTuples.b_MuIsoCombRel); */
  /* t->SetBranchAddress("MuTkNormChi2", tTuples.MuTkNormChi2, &tTuples.b_MuTkNormChi2); */
  /* t->SetBranchAddress("MuTkHitCnt", tTuples.MuTkHitCnt, &tTuples.b_MuTkHitCnt); */
  /* t->SetBranchAddress("MuMatchedStationCnt", tTuples.MuMatchedStationCnt, &tTuples.b_MuMatchedStationCnt); */
  /* t->SetBranchAddress("MuDz", tTuples.MuDz, &tTuples.b_MuDz); */
  /* t->SetBranchAddress("MuPixelHitCnt", tTuples.MuPixelHitCnt, &tTuples.b_MuPixelHitCnt); */
  /* t->SetBranchAddress("MuTkLayerCnt", tTuples.MuTkLayerCnt, &tTuples.b_MuTkLayerCnt); */
  /* t->SetBranchAddress("MuPfIsoChHad", tTuples.MuPfIsoChHad, &tTuples.b_MuPfIsoChHad); */
  /* t->SetBranchAddress("MuPfIsoNeutralHad", tTuples.MuPfIsoNeutralHad, &tTuples.b_MuPfIsoNeutralHad); */
  /* t->SetBranchAddress("MuPfIsoRawRel", tTuples.MuPfIsoRawRel, &tTuples.b_MuPfIsoRawRel); */
  /* t->SetBranchAddress("MuHltMatch", tTuples.MuHltMatch, &tTuples.b_MuHltMatch); */
  /* t->SetBranchAddress("ElPt", tTuples.ElPt, &tTuples.b_ElPt); */
  /* t->SetBranchAddress("ElEta", tTuples.ElEta, &tTuples.b_ElEta); */
  /* t->SetBranchAddress("ElEtaSc", tTuples.ElEtaSc, &tTuples.b_ElEtaSc); */
  /* t->SetBranchAddress("ElPhi", tTuples.ElPhi, &tTuples.b_ElPhi); */
  /* t->SetBranchAddress("ElE", tTuples.ElE, &tTuples.b_ElE); */
  /* t->SetBranchAddress("ElId", tTuples.ElId, &tTuples.b_ElId); */
  /* t->SetBranchAddress("ElCh", tTuples.ElCh, &tTuples.b_ElCh); */
  /* t->SetBranchAddress("ElMvaTrig", tTuples.ElMvaTrig, &tTuples.b_ElMvaTrig); */
  /* t->SetBranchAddress("ElMvaNonTrig", tTuples.ElMvaNonTrig, &tTuples.b_ElMvaNonTrig); */
  /* t->SetBranchAddress("ElMvaPresel", tTuples.ElMvaPresel, &tTuples.b_ElMvaPresel); */
  /* t->SetBranchAddress("ElDEtaTkScAtVtx", tTuples.ElDEtaTkScAtVtx, &tTuples.b_ElDEtaTkScAtVtx); */
  /* t->SetBranchAddress("ElDPhiTkScAtVtx", tTuples.ElDPhiTkScAtVtx, &tTuples.b_ElDPhiTkScAtVtx); */
  /* t->SetBranchAddress("ElHoE", tTuples.ElHoE, &tTuples.b_ElHoE); */
  /* t->SetBranchAddress("ElSigmaIetaIeta", tTuples.ElSigmaIetaIeta, &tTuples.b_ElSigmaIetaIeta); */
  /* t->SetBranchAddress("ElSigmaIetaIetaFull5x5", tTuples.ElSigmaIetaIetaFull5x5, &tTuples.b_ElSigmaIetaIetaFull5x5); */
  /* t->SetBranchAddress("ElEinvMinusPinv", tTuples.ElEinvMinusPinv, &tTuples.b_ElEinvMinusPinv); */
  /* t->SetBranchAddress("ElD0", tTuples.ElD0, &tTuples.b_ElD0); */
  /* t->SetBranchAddress("ElDz", tTuples.ElDz, &tTuples.b_ElDz); */
  /* t->SetBranchAddress("ElExpectedMissingInnerHitCnt", tTuples.ElExpectedMissingInnerHitCnt, &tTuples.b_ElExpectedMissingInnerHitCnt); */
  /* t->SetBranchAddress("ElPassConvVeto", tTuples.ElPassConvVeto, &tTuples.b_ElPassConvVeto); */
  /* t->SetBranchAddress("ElHltMatch", tTuples.ElHltMatch, &tTuples.b_ElHltMatch); */
  /* t->SetBranchAddress("ElPfIsoChHad", tTuples.ElPfIsoChHad, &tTuples.b_ElPfIsoChHad); */
  /* t->SetBranchAddress("ElPfIsoNeutralHad", tTuples.ElPfIsoNeutralHad, &tTuples.b_ElPfIsoNeutralHad); */
  /* t->SetBranchAddress("ElPfIsoIso", tTuples.ElPfIsoIso, &tTuples.b_ElPfIsoIso); */
  /* t->SetBranchAddress("ElPfIsoPuChHad", tTuples.ElPfIsoPuChHad, &tTuples.b_ElPfIsoPuChHad); */
  /* t->SetBranchAddress("ElPfIsoRaw", tTuples.ElPfIsoRaw, &tTuples.b_ElPfIsoRaw); */
  /* t->SetBranchAddress("ElPfIsoDbeta", tTuples.ElPfIsoDbeta, &tTuples.b_ElPfIsoDbeta); */
  /* t->SetBranchAddress("ElAEff", tTuples.ElAEff, &tTuples.b_ElAEff); */
  /* t->SetBranchAddress("charged", tTuples.charged, &tTuples.b_charged); */
  /* t->SetBranchAddress("photon", tTuples.photon, &tTuples.b_photon); */
  /* t->SetBranchAddress("neutral", tTuples.neutral, &tTuples.b_neutral); */
  /* t->SetBranchAddress("charged_Tom", tTuples.charged_Tom, &tTuples.b_charged_Tom); */
  /* t->SetBranchAddress("photon_Tom", tTuples.photon_Tom, &tTuples.b_photon_Tom); */
  /* t->SetBranchAddress("neutral_Tom", tTuples.neutral_Tom, &tTuples.b_neutral_Tom); */
  /* t->SetBranchAddress("PhotPt", tTuples.PhotPt, &tTuples.b_PhotPt); */
  /* t->SetBranchAddress("PhotEta", tTuples.PhotEta, &tTuples.b_PhotEta); */
  /* t->SetBranchAddress("PhotPhi", tTuples.PhotPhi, &tTuples.b_PhotPhi); */
  /* t->SetBranchAddress("PhotScRawE", tTuples.PhotScRawE, &tTuples.b_PhotScRawE); */
  /* t->SetBranchAddress("PhotScEta", tTuples.PhotScEta, &tTuples.b_PhotScEta); */
  /* t->SetBranchAddress("PhotScPhi", tTuples.PhotScPhi, &tTuples.b_PhotScPhi); */
  /* t->SetBranchAddress("PhotIsoEcal", tTuples.PhotIsoEcal, &tTuples.b_PhotIsoEcal); */
  /* t->SetBranchAddress("PhotIsoHcal", tTuples.PhotIsoHcal, &tTuples.b_PhotIsoHcal); */
  /* t->SetBranchAddress("PhotIsoTk", tTuples.PhotIsoTk, &tTuples.b_PhotIsoTk); */
  /* t->SetBranchAddress("PhotPfIsoChHad", tTuples.PhotPfIsoChHad, &tTuples.b_PhotPfIsoChHad); */
  /* t->SetBranchAddress("PhotPfIsoNeutralHad", tTuples.PhotPfIsoNeutralHad, &tTuples.b_PhotPfIsoNeutralHad); */
  /* t->SetBranchAddress("PhotPfIsoPhot", tTuples.PhotPfIsoPhot, &tTuples.b_PhotPfIsoPhot); */
  /* t->SetBranchAddress("PhotPfIsoPuChHad", tTuples.PhotPfIsoPuChHad, &tTuples.b_PhotPfIsoPuChHad); */
  /* t->SetBranchAddress("PhotPfIsoEcalClus", tTuples.PhotPfIsoEcalClus, &tTuples.b_PhotPfIsoEcalClus); */
  /* t->SetBranchAddress("PhotPfIsoHcalClus", tTuples.PhotPfIsoHcalClus, &tTuples.b_PhotPfIsoHcalClus); */
  /* t->SetBranchAddress("PhotE3x3", tTuples.PhotE3x3, &tTuples.b_PhotE3x3); */
  /* t->SetBranchAddress("PhotE1x5", tTuples.PhotE1x5, &tTuples.b_PhotE1x5); */
  /* t->SetBranchAddress("PhotE2x5", tTuples.PhotE2x5, &tTuples.b_PhotE2x5); */
  /* t->SetBranchAddress("PhotE5x5", tTuples.PhotE5x5, &tTuples.b_PhotE5x5); */
  /* t->SetBranchAddress("PhotSigmaIetaIeta", tTuples.PhotSigmaIetaIeta, &tTuples.b_PhotSigmaIetaIeta); */
  /* t->SetBranchAddress("PhotEtaWidth", tTuples.PhotEtaWidth, &tTuples.b_PhotEtaWidth); */
  /* t->SetBranchAddress("PhotPhiWidth", tTuples.PhotPhiWidth, &tTuples.b_PhotPhiWidth); */
  /* t->SetBranchAddress("PhotHoE", tTuples.PhotHoE, &tTuples.b_PhotHoE); */
  /* t->SetBranchAddress("PhotId", tTuples.PhotId, &tTuples.b_PhotId); */
  /* t->SetBranchAddress("PhotHasPixelSeed", tTuples.PhotHasPixelSeed, &tTuples.b_PhotHasPixelSeed); */
  /* t->SetBranchAddress("JetAk04Pt", tTuples.JetAk04Pt, &tTuples.b_JetAk04Pt); */
  /* t->SetBranchAddress("JetAk04Eta", tTuples.JetAk04Eta, &tTuples.b_JetAk04Eta); */
  /* t->SetBranchAddress("JetAk04Phi", tTuples.JetAk04Phi, &tTuples.b_JetAk04Phi); */
  /* t->SetBranchAddress("JetAk04E", tTuples.JetAk04E, &tTuples.b_JetAk04E); */
  /* t->SetBranchAddress("JetAk04Id", tTuples.JetAk04Id, &tTuples.b_JetAk04Id); */
  /* t->SetBranchAddress("JetAk04PuId", tTuples.JetAk04PuId, &tTuples.b_JetAk04PuId); */
  /* t->SetBranchAddress("JetAk04PuMva", tTuples.JetAk04PuMva, &tTuples.b_JetAk04PuMva); */
  /* t->SetBranchAddress("JetAk04RawPt", tTuples.JetAk04RawPt, &tTuples.b_JetAk04RawPt); */
  /* t->SetBranchAddress("JetAk04RawE", tTuples.JetAk04RawE, &tTuples.b_JetAk04RawE); */
  /* t->SetBranchAddress("JetAk04HfHadE", tTuples.JetAk04HfHadE, &tTuples.b_JetAk04HfHadE); */
  /* t->SetBranchAddress("JetAk04HfEmE", tTuples.JetAk04HfEmE, &tTuples.b_JetAk04HfEmE); */
  /* t->SetBranchAddress("JetAk04ChHadFrac", tTuples.JetAk04ChHadFrac, &tTuples.b_JetAk04ChHadFrac); */
  /* t->SetBranchAddress("JetAk04NeutralHadAndHfFrac", tTuples.JetAk04NeutralHadAndHfFrac, &tTuples.b_JetAk04NeutralHadAndHfFrac); */
  /* t->SetBranchAddress("JetAk04ChEmFrac", tTuples.JetAk04ChEmFrac, &tTuples.b_JetAk04ChEmFrac); */
  /* t->SetBranchAddress("JetAk04NeutralEmFrac", tTuples.JetAk04NeutralEmFrac, &tTuples.b_JetAk04NeutralEmFrac); */
  /* t->SetBranchAddress("JetAk04ChMult", tTuples.JetAk04ChMult, &tTuples.b_JetAk04ChMult); */
  /* t->SetBranchAddress("JetAk04ConstCnt", tTuples.JetAk04ConstCnt, &tTuples.b_JetAk04ConstCnt); */
  /* t->SetBranchAddress("JetAk04JetBeta", tTuples.JetAk04JetBeta, &tTuples.b_JetAk04JetBeta); */
  /* t->SetBranchAddress("JetAk04JetBetaClassic", tTuples.JetAk04JetBetaClassic, &tTuples.b_JetAk04JetBetaClassic); */
  /* t->SetBranchAddress("JetAk04JetBetaStar", tTuples.JetAk04JetBetaStar, &tTuples.b_JetAk04JetBetaStar); */
  /* t->SetBranchAddress("JetAk04JetBetaStarClassic", tTuples.JetAk04JetBetaStarClassic, &tTuples.b_JetAk04JetBetaStarClassic); */
  /* t->SetBranchAddress("JetAk04BTagCsv", tTuples.JetAk04BTagCsv, &tTuples.b_JetAk04BTagCsv); */
  /* t->SetBranchAddress("JetAk04BTagCsvV1", tTuples.JetAk04BTagCsvV1, &tTuples.b_JetAk04BTagCsvV1); */
  /* t->SetBranchAddress("JetAk04BTagCsvSLV1", tTuples.JetAk04BTagCsvSLV1, &tTuples.b_JetAk04BTagCsvSLV1); */
  /* t->SetBranchAddress("JetAk04BDiscCisvV2", tTuples.JetAk04BDiscCisvV2, &tTuples.b_JetAk04BDiscCisvV2); */
  /* t->SetBranchAddress("JetAk04BDiscJp", tTuples.JetAk04BDiscJp, &tTuples.b_JetAk04BDiscJp); */
  /* t->SetBranchAddress("JetAk04BDiscBjp", tTuples.JetAk04BDiscBjp, &tTuples.b_JetAk04BDiscBjp); */
  /* t->SetBranchAddress("JetAk04BDiscTche", tTuples.JetAk04BDiscTche, &tTuples.b_JetAk04BDiscTche); */
  /* t->SetBranchAddress("JetAk04BDiscTchp", tTuples.JetAk04BDiscTchp, &tTuples.b_JetAk04BDiscTchp); */
  /* t->SetBranchAddress("JetAk04BDiscSsvhe", tTuples.JetAk04BDiscSsvhe, &tTuples.b_JetAk04BDiscSsvhe); */
  /* t->SetBranchAddress("JetAk04BDiscSsvhp", tTuples.JetAk04BDiscSsvhp, &tTuples.b_JetAk04BDiscSsvhp); */
  /* t->SetBranchAddress("JetAk04PartFlav", tTuples.JetAk04PartFlav, &tTuples.b_JetAk04PartFlav); */
  /* t->SetBranchAddress("JetAk04JecUncUp", tTuples.JetAk04JecUncUp, &tTuples.b_JetAk04JecUncUp); */
  /* t->SetBranchAddress("JetAk04JecUncDwn", tTuples.JetAk04JecUncDwn, &tTuples.b_JetAk04JecUncDwn); */
  /* t->SetBranchAddress("JetAk04ConstId", tTuples.JetAk04ConstId, &tTuples.b_JetAk04ConstId); */
  /* t->SetBranchAddress("JetAk04ConstPt", tTuples.JetAk04ConstPt, &tTuples.b_JetAk04ConstPt); */
  /* t->SetBranchAddress("JetAk04ConstEta", tTuples.JetAk04ConstEta, &tTuples.b_JetAk04ConstEta); */
  /* t->SetBranchAddress("JetAk04ConstPhi", tTuples.JetAk04ConstPhi, &tTuples.b_JetAk04ConstPhi); */
  /* t->SetBranchAddress("JetAk04ConstE", tTuples.JetAk04ConstE, &tTuples.b_JetAk04ConstE); */
  /* t->SetBranchAddress("JetAk04GenJet", tTuples.JetAk04GenJet, &tTuples.b_JetAk04GenJet); */
  if (doCheck) {
  }
}
