#define EventTree_cxx

#include "EventTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

ClassImp(EventTree)


EventTree::EventTree(const char *name)
{
  std::cout <<"  Tuple Name : " << name << std::endl;
  std::cout <<"  Needs to define other paramters : " << name << std::endl;
}
EventTree::EventTree(const char *name, const char *infName) :  fInfile(0), fTree(0), fname(name)
{
  fInfile = TFile::Open(infName);
  if( !fInfile ) { std::cout <<" No input file found " << std::endl; exit(1); }

  fTree = (TTree*)fInfile->Get("tupel/EventTree");

  Init(fTree);
}
EventTree::EventTree(const char *name, TTree *inTree) : fname(name)
{
  fTree = inTree;
  Init(fTree);
}
EventTree::~EventTree()
{
  if (!fTree) return;
  delete fTree->GetCurrentFile();
}
Int_t EventTree::GetEntry(Long64_t entry)
{
  // Read contents of entry.                                                                                                            
  if (!fTree) return 0;
  return fTree->GetEntry(entry);
}
Long64_t EventTree::GetEntries()
{
  //!   # of entries
  if (!fTree) return 0;
  return fTree->GetEntries();
}
void EventTree::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize 
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.       
  // It is normally not necessary to make changes to the generated   
  // code, but the routine can be extended by the user if needed.  
  // Init() will be called many times when running on PROOF        
  // (once per file to be processed).                 

  // Set object pointer              
  EvtWeights = 0;
  METPt = 0;
  METPx = 0;
  METPy = 0;
  METPz = 0;
  METE = 0;
  METsigx2 = 0;
  METsigxy = 0;
  METsigy2 = 0;
  METsig = 0;
  GLepDr01Pt = 0;
  GLepDr01Eta = 0;
  GLepDr01Phi = 0;
  GLepDr01E = 0;
  GLepDr01Id = 0;
  GLepDr01St = 0;
  GLepDr01MomId = 0;
  GLepBarePt = 0;
  GLepBareEta = 0;
  GLepBarePhi = 0;
  GLepBareE = 0;
  GLepBareId = 0;
  GLepBareSt = 0;
  GLepBareMomId = 0;
  GLepSt3Pt = 0;
  GLepSt3Eta = 0;
  GLepSt3Phi = 0;
  GLepSt3E = 0;
  GLepSt3Id = 0;
  GLepSt3St = 0;
  GLepSt3Mother0Id = 0;
  GLepSt3MotherCnt = 0;
  GPhotPt = 0;
  GPhotEta = 0;
  GPhotPhi = 0;
  GPhotE = 0;
  GPhotMotherId = 0;
  GPhotSt = 0;
  GPhotIsoEDR03 = 0;
  GPhotIsoEDR04 = 0;
  GPhotIsoEDR05 = 0;
  GPhotIsoSumPtDR03 = 0;
  GPhotIsoSumPtDR04 = 0;
  GPhotIsoSumPtDR05 = 0;
  GLepClosePhotPt = 0;
  GLepClosePhotEta = 0;
  GLepClosePhotPhi = 0;
  GLepClosePhotE = 0;
  GLepClosePhotId = 0;
  GLepClosePhotMother0Id = 0;
  GLepClosePhotMotherCnt = 0;
  GLepClosePhotSt = 0;
  GJetAk04Pt = 0;
  GJetAk04Eta = 0;
  GJetAk04Phi = 0;
  GJetAk04E = 0;
  GJetAk04ChFrac = 0;
  GJetAk04ConstCnt = 0;
  GJetAk04ConstId = 0;
  GJetAk04ConstPt = 0;
  GJetAk04ConstEta = 0;
  GJetAk04ConstPhi = 0;
  GJetAk04ConstE = 0;
  GPdfId1 = 0;
  GPdfId2 = 0;
  GPdfx1 = 0;
  GPdfx2 = 0;
  GPdfScale = 0;
  MuPt = 0;
  MuEta = 0;
  MuPhi = 0;
  MuE = 0;
  MuId = 0;
  MuIdTight = 0;
  MuCh = 0;
  MuVtxZ = 0;
  MuDxy = 0;
  MuPfIso = 0;
  MuType = 0;
  MuIsoTkIsoAbs = 0;
  MuIsoTkIsoRel = 0;
  MuIsoCalAbs = 0;
  MuIsoCombRel = 0;
  MuTkNormChi2 = 0;
  MuTkHitCnt = 0;
  MuMatchedStationCnt = 0;
  MuDz = 0;
  MuPixelHitCnt = 0;
  MuTkLayerCnt = 0;
  MuPfIsoChHad = 0;
  MuPfIsoNeutralHad = 0;
  MuPfIsoRawRel = 0;
  MuHltMatch = 0;
  ElPt = 0;
  ElEta = 0;
  ElEtaSc = 0;
  ElPhi = 0;
  ElE = 0;
  ElId = 0;
  ElCh = 0;
  ElMvaTrig = 0;
  ElMvaNonTrig = 0;
  ElMvaPresel = 0;
  ElDEtaTkScAtVtx = 0;
  ElDPhiTkScAtVtx = 0;
  ElHoE = 0;
  ElSigmaIetaIeta = 0;
  ElSigmaIetaIetaFull5x5 = 0;
  ElEinvMinusPinv = 0;
  ElD0 = 0;
  ElDz = 0;
  ElExpectedMissingInnerHitCnt = 0;
  ElPassConvVeto = 0;
  ElHltMatch = 0;
  ElPfIsoChHad = 0;
  ElPfIsoNeutralHad = 0;
  ElPfIsoIso = 0;
  ElPfIsoPuChHad = 0;
  ElPfIsoRaw = 0;
  ElPfIsoDbeta = 0;
  ElAEff = 0;
  charged = 0;
  photon = 0;
  neutral = 0;
  charged_Tom = 0;
  photon_Tom = 0;
  neutral_Tom = 0;
  PhotPt = 0;
  PhotEta = 0;
  PhotPhi = 0;
  PhotScRawE = 0;
  PhotScEta = 0;
  PhotScPhi = 0;
  PhotIsoEcal = 0;
  PhotIsoHcal = 0;
  PhotIsoTk = 0;
  PhotPfIsoChHad = 0;
  PhotPfIsoNeutralHad = 0;
  PhotPfIsoPhot = 0;
  PhotPfIsoPuChHad = 0;
  PhotPfIsoEcalClus = 0;
  PhotPfIsoHcalClus = 0;
  PhotE3x3 = 0;
  PhotE1x5 = 0;
  PhotE2x5 = 0;
  PhotE5x5 = 0;
  PhotSigmaIetaIeta = 0;
  PhotEtaWidth = 0;
  PhotPhiWidth = 0;
  PhotHoE = 0;
  PhotId = 0;
  PhotHasPixelSeed = 0;
  JetAk04Pt = 0;
  JetAk04Eta = 0;
  JetAk04Phi = 0;
  JetAk04E = 0;
  JetAk04Id = 0;
  JetAk04PuId = 0;
  JetAk04PuMva = 0;
  JetAk04RawPt = 0;
  JetAk04RawE = 0;
  JetAk04HfHadE = 0;
  JetAk04HfEmE = 0;
  JetAk04ChHadFrac = 0;
  JetAk04NeutralHadAndHfFrac = 0;
  JetAk04ChEmFrac = 0;
  JetAk04NeutralEmFrac = 0;
  JetAk04ChMult = 0;
  JetAk04ConstCnt = 0;
  JetAk04JetBeta = 0;
  JetAk04JetBetaClassic = 0;
  JetAk04JetBetaStar = 0;
  JetAk04JetBetaStarClassic = 0;
  JetAk04BTagCsv = 0;
  JetAk04BTagCsvV1 = 0;
  JetAk04BTagCsvSLV1 = 0;
  JetAk04BDiscCisvV2 = 0;
  JetAk04BDiscJp = 0;
  JetAk04BDiscBjp = 0;
  JetAk04BDiscTche = 0;
  JetAk04BDiscTchp = 0;
  JetAk04BDiscSsvhe = 0;
  JetAk04BDiscSsvhp = 0;
  JetAk04PartFlav = 0;
  JetAk04JecUncUp = 0;
  JetAk04JecUncDwn = 0;
  JetAk04ConstId = 0;
  JetAk04ConstPt = 0;
  JetAk04ConstEta = 0;
  JetAk04ConstPhi = 0;
  JetAk04ConstE = 0;
  JetAk04GenJet = 0;



  // Set branch addresses and branch pointers                                                                                                                            
  if (!tree) return;
  fCurrent = -1;
  tree->SetMakeClass(1);

  tree->SetBranchAddress("EvtIsRealData", &EvtIsRealData, &b_EvtIsRealData);
  tree->SetBranchAddress("EvtNum", &EvtNum, &b_EvtNum);
  tree->SetBranchAddress("EvtRunNum", &EvtRunNum, &b_EvtRunNum);
  tree->SetBranchAddress("EvtLumiNum", &EvtLumiNum, &b_EvtLumiNum);
  tree->SetBranchAddress("EvtBxNum", &EvtBxNum, &b_EvtBxNum);
  tree->SetBranchAddress("EvtVtxCnt", &EvtVtxCnt, &b_EvtVtxCnt);
  tree->SetBranchAddress("EvtPuCnt", &EvtPuCnt, &b_EvtPuCnt);
  tree->SetBranchAddress("EvtPuCntTruth", &EvtPuCntTruth, &b_EvtPuCntTruth);
  tree->SetBranchAddress("EvtWeights", &EvtWeights, &b_EvtWeights);
  tree->SetBranchAddress("TrigHlt", &TrigHlt, &b_TrigHlt);
  tree->SetBranchAddress("TrigHltMu", &TrigHltMu, &b_TrigHltMu);
  tree->SetBranchAddress("TrigHltDiMu", &TrigHltDiMu, &b_TrigHltDiMu);
  tree->SetBranchAddress("TrigHltEl", &TrigHltEl, &b_TrigHltEl);
  tree->SetBranchAddress("TrigHltDiEl", &TrigHltDiEl, &b_TrigHltDiEl);
  tree->SetBranchAddress("TrigHltElMu", &TrigHltElMu, &b_TrigHltElMu);
  tree->SetBranchAddress("METPt", &METPt, &b_METPt);
  tree->SetBranchAddress("METPx", &METPx, &b_METPx);
  tree->SetBranchAddress("METPy", &METPy, &b_METPy);
  tree->SetBranchAddress("METPz", &METPz, &b_METPz);
  tree->SetBranchAddress("METE", &METE, &b_METE);
  tree->SetBranchAddress("METsigx2", &METsigx2, &b_METsigx2);
  tree->SetBranchAddress("METsigxy", &METsigxy, &b_METsigxy);
  tree->SetBranchAddress("METsigy2", &METsigy2, &b_METsigy2);
  tree->SetBranchAddress("METsig", &METsig, &b_METsig);
  tree->SetBranchAddress("GLepDr01Pt", &GLepDr01Pt, &b_GLepDr01Pt);
  tree->SetBranchAddress("GLepDr01Eta", &GLepDr01Eta, &b_GLepDr01Eta);
  tree->SetBranchAddress("GLepDr01Phi", &GLepDr01Phi, &b_GLepDr01Phi);
  tree->SetBranchAddress("GLepDr01E", &GLepDr01E, &b_GLepDr01E);
  tree->SetBranchAddress("GLepDr01Id", &GLepDr01Id, &b_GLepDr01Id);
  tree->SetBranchAddress("GLepDr01St", &GLepDr01St, &b_GLepDr01St);
  tree->SetBranchAddress("GLepDr01MomId", &GLepDr01MomId, &b_GLepDr01MomId);
  tree->SetBranchAddress("GLepBarePt", &GLepBarePt, &b_GLepBarePt);
  tree->SetBranchAddress("GLepBareEta", &GLepBareEta, &b_GLepBareEta);
  tree->SetBranchAddress("GLepBarePhi", &GLepBarePhi, &b_GLepBarePhi);
  tree->SetBranchAddress("GLepBareE", &GLepBareE, &b_GLepBareE);
  tree->SetBranchAddress("GLepBareId", &GLepBareId, &b_GLepBareId);
  tree->SetBranchAddress("GLepBareSt", &GLepBareSt, &b_GLepBareSt);
  tree->SetBranchAddress("GLepBareMomId", &GLepBareMomId, &b_GLepBareMomId);
  tree->SetBranchAddress("GLepSt3Pt", &GLepSt3Pt, &b_GLepSt3Pt);
  tree->SetBranchAddress("GLepSt3Eta", &GLepSt3Eta, &b_GLepSt3Eta);
  tree->SetBranchAddress("GLepSt3Phi", &GLepSt3Phi, &b_GLepSt3Phi);
  tree->SetBranchAddress("GLepSt3E", &GLepSt3E, &b_GLepSt3E);
  tree->SetBranchAddress("GLepSt3Id", &GLepSt3Id, &b_GLepSt3Id);
  tree->SetBranchAddress("GLepSt3St", &GLepSt3St, &b_GLepSt3St);
  tree->SetBranchAddress("GLepSt3Mother0Id", &GLepSt3Mother0Id, &b_GLepSt3Mother0Id);
  tree->SetBranchAddress("GPhotSt", &GPhotSt, &b_GPhotSt);
  tree->SetBranchAddress("GPhotIsoEDR03", &GPhotIsoEDR03, &b_GPhotIsoEDR03);
  tree->SetBranchAddress("GPhotIsoEDR04", &GPhotIsoEDR04, &b_GPhotIsoEDR04);
  tree->SetBranchAddress("GPhotIsoEDR05", &GPhotIsoEDR05, &b_GPhotIsoEDR05);
  tree->SetBranchAddress("GPhotIsoSumPtDR03", &GPhotIsoSumPtDR03, &b_GPhotIsoSumPtDR03);
  tree->SetBranchAddress("GPhotIsoSumPtDR04", &GPhotIsoSumPtDR04, &b_GPhotIsoSumPtDR04);
  tree->SetBranchAddress("GPhotIsoSumPtDR05", &GPhotIsoSumPtDR05, &b_GPhotIsoSumPtDR05);
  tree->SetBranchAddress("GLepClosePhotPt", &GLepClosePhotPt, &b_GLepClosePhotPt);
  tree->SetBranchAddress("GLepClosePhotEta", &GLepClosePhotEta, &b_GLepClosePhotEta);
  tree->SetBranchAddress("GLepClosePhotPhi", &GLepClosePhotPhi, &b_GLepClosePhotPhi);
  tree->SetBranchAddress("GLepClosePhotE", &GLepClosePhotE, &b_GLepClosePhotE);
  tree->SetBranchAddress("GLepClosePhotId", &GLepClosePhotId, &b_GLepClosePhotId);
  tree->SetBranchAddress("GLepClosePhotMother0Id", &GLepClosePhotMother0Id, &b_GLepClosePhotMother0Id);
  tree->SetBranchAddress("GLepClosePhotMotherCnt", &GLepClosePhotMotherCnt, &b_GLepClosePhotMotherCnt);
  tree->SetBranchAddress("GLepClosePhotSt", &GLepClosePhotSt, &b_GLepClosePhotSt);
  tree->SetBranchAddress("GJetAk04Pt", &GJetAk04Pt, &b_GJetAk04Pt);
  tree->SetBranchAddress("GJetAk04Eta", &GJetAk04Eta, &b_GJetAk04Eta);
  tree->SetBranchAddress("GJetAk04Phi", &GJetAk04Phi, &b_GJetAk04Phi);
  tree->SetBranchAddress("GJetAk04E", &GJetAk04E, &b_GJetAk04E);
  tree->SetBranchAddress("GJetAk04ChFrac", &GJetAk04ChFrac, &b_GJetAk04ChFrac);
  tree->SetBranchAddress("GJetAk04ConstCnt", &GJetAk04ConstCnt, &b_GJetAk04ConstCnt);
  tree->SetBranchAddress("GJetAk04ConstId", &GJetAk04ConstId, &b_GJetAk04ConstId);
  tree->SetBranchAddress("GJetAk04ConstPt", &GJetAk04ConstPt, &b_GJetAk04ConstPt);
  tree->SetBranchAddress("GJetAk04ConstEta", &GJetAk04ConstEta, &b_GJetAk04ConstEta);
  tree->SetBranchAddress("GJetAk04ConstPhi", &GJetAk04ConstPhi, &b_GJetAk04ConstPhi);
  tree->SetBranchAddress("GJetAk04ConstE", &GJetAk04ConstE, &b_GJetAk04ConstE);
  tree->SetBranchAddress("GPdfId1", &GPdfId1, &b_GPdfId1);
  tree->SetBranchAddress("GPdfId2", &GPdfId2, &b_GPdfId2);
  tree->SetBranchAddress("GPdfx1", &GPdfx1, &b_GPdfx1);
  tree->SetBranchAddress("GPdfx2", &GPdfx2, &b_GPdfx2);
  tree->SetBranchAddress("GPdfScale", &GPdfScale, &b_GPdfScale);
  tree->SetBranchAddress("GBinningValue", &GBinningValue, &b_GBinningValue);
  tree->SetBranchAddress("GNup", &GNup, &b_GNup);
  tree->SetBranchAddress("MuPt", &MuPt, &b_MuPt);
  tree->SetBranchAddress("MuEta", &MuEta, &b_MuEta);
  tree->SetBranchAddress("MuPhi", &MuPhi, &b_MuPhi);
  tree->SetBranchAddress("MuE", &MuE, &b_MuE);
  tree->SetBranchAddress("MuId", &MuId, &b_MuId);
  tree->SetBranchAddress("MuIdTight", &MuIdTight, &b_MuIdTight);
  tree->SetBranchAddress("MuCh", &MuCh, &b_MuCh);
  tree->SetBranchAddress("MuVtxZ", &MuVtxZ, &b_MuVtxZ);
  tree->SetBranchAddress("MuDxy", &MuDxy, &b_MuDxy);
  tree->SetBranchAddress("MuPfIso", &MuPfIso, &b_MuPfIso);
  tree->SetBranchAddress("MuType", &MuType, &b_MuType);
  tree->SetBranchAddress("MuIsoTkIsoAbs", &MuIsoTkIsoAbs, &b_MuIsoTkIsoAbs);
  tree->SetBranchAddress("MuIsoTkIsoRel", &MuIsoTkIsoRel, &b_MuIsoTkIsoRel);
  tree->SetBranchAddress("MuIsoCalAbs", &MuIsoCalAbs, &b_MuIsoCalAbs);
  tree->SetBranchAddress("MuIsoCombRel", &MuIsoCombRel, &b_MuIsoCombRel);
  tree->SetBranchAddress("MuTkNormChi2", &MuTkNormChi2, &b_MuTkNormChi2);
  tree->SetBranchAddress("MuTkHitCnt", &MuTkHitCnt, &b_MuTkHitCnt);
  tree->SetBranchAddress("MuMatchedStationCnt", &MuMatchedStationCnt, &b_MuMatchedStationCnt);
  tree->SetBranchAddress("MuDz", &MuDz, &b_MuDz);
  tree->SetBranchAddress("MuPixelHitCnt", &MuPixelHitCnt, &b_MuPixelHitCnt);
  tree->SetBranchAddress("MuTkLayerCnt", &MuTkLayerCnt, &b_MuTkLayerCnt);
  tree->SetBranchAddress("MuPfIsoChHad", &MuPfIsoChHad, &b_MuPfIsoChHad);
  tree->SetBranchAddress("MuPfIsoNeutralHad", &MuPfIsoNeutralHad, &b_MuPfIsoNeutralHad);
  tree->SetBranchAddress("MuPfIsoRawRel", &MuPfIsoRawRel, &b_MuPfIsoRawRel);
  tree->SetBranchAddress("MuHltMatch", &MuHltMatch, &b_MuHltMatch);
  tree->SetBranchAddress("ElPt", &ElPt, &b_ElPt);
  tree->SetBranchAddress("ElEta", &ElEta, &b_ElEta);
  tree->SetBranchAddress("ElEtaSc", &ElEtaSc, &b_ElEtaSc);
  tree->SetBranchAddress("ElPhi", &ElPhi, &b_ElPhi);
  tree->SetBranchAddress("ElE", &ElE, &b_ElE);
  tree->SetBranchAddress("ElId", &ElId, &b_ElId);
  tree->SetBranchAddress("ElCh", &ElCh, &b_ElCh);
  tree->SetBranchAddress("ElMvaTrig", &ElMvaTrig, &b_ElMvaTrig);
  tree->SetBranchAddress("ElMvaNonTrig", &ElMvaNonTrig, &b_ElMvaNonTrig);
  tree->SetBranchAddress("ElMvaPresel", &ElMvaPresel, &b_ElMvaPresel);
  tree->SetBranchAddress("ElDEtaTkScAtVtx", &ElDEtaTkScAtVtx, &b_ElDEtaTkScAtVtx);
  tree->SetBranchAddress("ElDPhiTkScAtVtx", &ElDPhiTkScAtVtx, &b_ElDPhiTkScAtVtx);
  tree->SetBranchAddress("ElHoE", &ElHoE, &b_ElHoE);
  tree->SetBranchAddress("ElSigmaIetaIeta", &ElSigmaIetaIeta, &b_ElSigmaIetaIeta);
  tree->SetBranchAddress("ElSigmaIetaIetaFull5x5", &ElSigmaIetaIetaFull5x5, &b_ElSigmaIetaIetaFull5x5);
  tree->SetBranchAddress("ElEinvMinusPinv", &ElEinvMinusPinv, &b_ElEinvMinusPinv);
  tree->SetBranchAddress("ElD0", &ElD0, &b_ElD0);
  tree->SetBranchAddress("ElDz", &ElDz, &b_ElDz);
  tree->SetBranchAddress("ElExpectedMissingInnerHitCnt", &ElExpectedMissingInnerHitCnt, &b_ElExpectedMissingInnerHitCnt);
  tree->SetBranchAddress("ElPassConvVeto", &ElPassConvVeto, &b_ElPassConvVeto);
  tree->SetBranchAddress("ElHltMatch", &ElHltMatch, &b_ElHltMatch);
  tree->SetBranchAddress("ElPfIsoChHad", &ElPfIsoChHad, &b_ElPfIsoChHad);
  tree->SetBranchAddress("ElPfIsoNeutralHad", &ElPfIsoNeutralHad, &b_ElPfIsoNeutralHad);
  tree->SetBranchAddress("ElPfIsoIso", &ElPfIsoIso, &b_ElPfIsoIso);
  tree->SetBranchAddress("ElPfIsoPuChHad", &ElPfIsoPuChHad, &b_ElPfIsoPuChHad);
  tree->SetBranchAddress("ElPfIsoRaw", &ElPfIsoRaw, &b_ElPfIsoRaw);
  tree->SetBranchAddress("ElPfIsoDbeta", &ElPfIsoDbeta, &b_ElPfIsoDbeta);
  tree->SetBranchAddress("ElAEff", &ElAEff, &b_ElAEff);
  tree->SetBranchAddress("charged", &charged, &b_charged);
  tree->SetBranchAddress("photon", &photon, &b_photon);
  tree->SetBranchAddress("neutral", &neutral, &b_neutral);
  tree->SetBranchAddress("charged_Tom", &charged_Tom, &b_charged_Tom);
  tree->SetBranchAddress("photon_Tom", &photon_Tom, &b_photon_Tom);
  tree->SetBranchAddress("neutral_Tom", &neutral_Tom, &b_neutral_Tom);
  tree->SetBranchAddress("PhotPt", &PhotPt, &b_PhotPt);
  tree->SetBranchAddress("PhotEta", &PhotEta, &b_PhotEta);
  tree->SetBranchAddress("PhotPhi", &PhotPhi, &b_PhotPhi);
  tree->SetBranchAddress("PhotScRawE", &PhotScRawE, &b_PhotScRawE);
  tree->SetBranchAddress("PhotScEta", &PhotScEta, &b_PhotScEta);
  tree->SetBranchAddress("PhotScPhi", &PhotScPhi, &b_PhotScPhi);
  tree->SetBranchAddress("PhotIsoEcal", &PhotIsoEcal, &b_PhotIsoEcal);
  tree->SetBranchAddress("PhotIsoHcal", &PhotIsoHcal, &b_PhotIsoHcal);
  tree->SetBranchAddress("PhotIsoTk", &PhotIsoTk, &b_PhotIsoTk);
  tree->SetBranchAddress("PhotPfIsoChHad", &PhotPfIsoChHad, &b_PhotPfIsoChHad);
  tree->SetBranchAddress("PhotPfIsoNeutralHad", &PhotPfIsoNeutralHad, &b_PhotPfIsoNeutralHad);
  tree->SetBranchAddress("PhotPfIsoPhot", &PhotPfIsoPhot, &b_PhotPfIsoPhot);
  tree->SetBranchAddress("PhotPfIsoPuChHad", &PhotPfIsoPuChHad, &b_PhotPfIsoPuChHad);
  tree->SetBranchAddress("PhotPfIsoEcalClus", &PhotPfIsoEcalClus, &b_PhotPfIsoEcalClus);
  tree->SetBranchAddress("PhotPfIsoHcalClus", &PhotPfIsoHcalClus, &b_PhotPfIsoHcalClus);
  tree->SetBranchAddress("PhotE3x3", &PhotE3x3, &b_PhotE3x3);
  tree->SetBranchAddress("PhotE1x5", &PhotE1x5, &b_PhotE1x5);
  tree->SetBranchAddress("PhotE2x5", &PhotE2x5, &b_PhotE2x5);
  tree->SetBranchAddress("PhotE5x5", &PhotE5x5, &b_PhotE5x5);
  tree->SetBranchAddress("PhotSigmaIetaIeta", &PhotSigmaIetaIeta, &b_PhotSigmaIetaIeta);
  tree->SetBranchAddress("PhotEtaWidth", &PhotEtaWidth, &b_PhotEtaWidth);
  tree->SetBranchAddress("PhotPhiWidth", &PhotPhiWidth, &b_PhotPhiWidth);
  tree->SetBranchAddress("PhotHoE", &PhotHoE, &b_PhotHoE);
  tree->SetBranchAddress("PhotId", &PhotId, &b_PhotId);
  tree->SetBranchAddress("PhotHasPixelSeed", &PhotHasPixelSeed, &b_PhotHasPixelSeed);
  tree->SetBranchAddress("JetAk04Pt", &JetAk04Pt, &b_JetAk04Pt);
  tree->SetBranchAddress("JetAk04Eta", &JetAk04Eta, &b_JetAk04Eta);
  tree->SetBranchAddress("JetAk04Phi", &JetAk04Phi, &b_JetAk04Phi);
  tree->SetBranchAddress("JetAk04E", &JetAk04E, &b_JetAk04E);
  tree->SetBranchAddress("JetAk04Id", &JetAk04Id, &b_JetAk04Id);
  tree->SetBranchAddress("JetAk04PuId", &JetAk04PuId, &b_JetAk04PuId);
  tree->SetBranchAddress("JetAk04PuMva", &JetAk04PuMva, &b_JetAk04PuMva);
  tree->SetBranchAddress("JetAk04RawPt", &JetAk04RawPt, &b_JetAk04RawPt);
  tree->SetBranchAddress("JetAk04RawE", &JetAk04RawE, &b_JetAk04RawE);
  tree->SetBranchAddress("JetAk04HfHadE", &JetAk04HfHadE, &b_JetAk04HfHadE);
  tree->SetBranchAddress("JetAk04HfEmE", &JetAk04HfEmE, &b_JetAk04HfEmE);
  tree->SetBranchAddress("JetAk04ChHadFrac", &JetAk04ChHadFrac, &b_JetAk04ChHadFrac);
  tree->SetBranchAddress("JetAk04NeutralHadAndHfFrac", &JetAk04NeutralHadAndHfFrac, &b_JetAk04NeutralHadAndHfFrac);
  tree->SetBranchAddress("JetAk04ChEmFrac", &JetAk04ChEmFrac, &b_JetAk04ChEmFrac);
  tree->SetBranchAddress("JetAk04NeutralEmFrac", &JetAk04NeutralEmFrac, &b_JetAk04NeutralEmFrac);
  tree->SetBranchAddress("JetAk04ChMult", &JetAk04ChMult, &b_JetAk04ChMult);
  tree->SetBranchAddress("JetAk04ConstCnt", &JetAk04ConstCnt, &b_JetAk04ConstCnt);
  tree->SetBranchAddress("JetAk04JetBeta", &JetAk04JetBeta, &b_JetAk04JetBeta);
  tree->SetBranchAddress("JetAk04JetBetaClassic", &JetAk04JetBetaClassic, &b_JetAk04JetBetaClassic);
  tree->SetBranchAddress("JetAk04JetBetaStar", &JetAk04JetBetaStar, &b_JetAk04JetBetaStar);
  tree->SetBranchAddress("JetAk04JetBetaStarClassic", &JetAk04JetBetaStarClassic, &b_JetAk04JetBetaStarClassic);
  tree->SetBranchAddress("JetAk04BTagCsv", &JetAk04BTagCsv, &b_JetAk04BTagCsv);
  tree->SetBranchAddress("JetAk04BTagCsvV1", &JetAk04BTagCsvV1, &b_JetAk04BTagCsvV1);
  tree->SetBranchAddress("JetAk04BTagCsvSLV1", &JetAk04BTagCsvSLV1, &b_JetAk04BTagCsvSLV1);
  tree->SetBranchAddress("JetAk04BDiscCisvV2", &JetAk04BDiscCisvV2, &b_JetAk04BDiscCisvV2);
  tree->SetBranchAddress("JetAk04BDiscJp", &JetAk04BDiscJp, &b_JetAk04BDiscJp);
  tree->SetBranchAddress("JetAk04BDiscBjp", &JetAk04BDiscBjp, &b_JetAk04BDiscBjp);
  tree->SetBranchAddress("JetAk04BDiscTche", &JetAk04BDiscTche, &b_JetAk04BDiscTche);
  tree->SetBranchAddress("JetAk04BDiscTchp", &JetAk04BDiscTchp, &b_JetAk04BDiscTchp);
  tree->SetBranchAddress("JetAk04BDiscSsvhe", &JetAk04BDiscSsvhe, &b_JetAk04BDiscSsvhe);
  tree->SetBranchAddress("JetAk04BDiscSsvhp", &JetAk04BDiscSsvhp, &b_JetAk04BDiscSsvhp);
  tree->SetBranchAddress("JetAk04PartFlav", &JetAk04PartFlav, &b_JetAk04PartFlav);
  tree->SetBranchAddress("JetAk04JecUncUp", &JetAk04JecUncUp, &b_JetAk04JecUncUp);
  tree->SetBranchAddress("JetAk04JecUncDwn", &JetAk04JecUncDwn, &b_JetAk04JecUncDwn);
  tree->SetBranchAddress("JetAk04ConstId", &JetAk04ConstId, &b_JetAk04ConstId);
  tree->SetBranchAddress("JetAk04ConstPt", &JetAk04ConstPt, &b_JetAk04ConstPt);
  tree->SetBranchAddress("JetAk04ConstEta", &JetAk04ConstEta, &b_JetAk04ConstEta);
  tree->SetBranchAddress("JetAk04ConstPhi", &JetAk04ConstPhi, &b_JetAk04ConstPhi);
  tree->SetBranchAddress("JetAk04ConstE", &JetAk04ConstE, &b_JetAk04ConstE);
  tree->SetBranchAddress("JetAk04GenJet", &JetAk04GenJet, &b_JetAk04GenJet);
}
