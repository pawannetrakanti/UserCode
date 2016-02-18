#include <iostream>

#include "hiForest.h"


ClassImp(HiForest)
HiForest::HiForest(const char *name, const char *infName, bool ispp, bool ishireco, bool ismc):
  tree(0),
  pp(ispp),
  hireco(ishireco),
  mc(ismc),
  nEntries(0),
  currentEvent(0),
  fname(name)
{
  tree = new TTree("tree",fname);

  //! Initialize the jet collection
  mJets = new Jets[ndir];
  
  // Input file
  inf = TFile::Open(infName);
  //inf = new TFile(infName,"r");

  hltTree       = (TTree*) inf->Get("hltanalysis/HltTree");
  hasHltTree    = (hltTree   != 0);
  if (hasHltTree) {
    hltTree->SetName("hlt");
    if (tree == 0) tree = hltTree; else tree->AddFriend(hltTree);
    setupHltTree(hltTree,hlt);
    //std::cout<<"SetupHlt Tree : " <<std::endl;
  }

  skimTree      = (TTree*) inf->Get("skimanalysis/HltTree");
  hasSkimTree   = (skimTree   != 0);
  if (hasSkimTree) {
    skimTree->SetName("skim");
    if (tree == 0) tree = skimTree; else tree->AddFriend(skimTree);
    setupSkimTree(skimTree,skim);
    //std::cout<<"SetupSkim Tree : " <<std::endl;
  }

  evtTree       = (TTree*) inf->Get("hiEvtAnalyzer/HiTree");
  hasEvtTree    = (evtTree   != 0);
  if (hasEvtTree) {
    evtTree->SetName("event");
    if (tree == 0) tree = evtTree; else tree->AddFriend(evtTree);
    setupEvtTree(evtTree,evt);
    //std::cout<<"SetupEvt Tree : " <<std::endl;
  }

  pfTree        = (TTree*) inf->Get("pfcandAnalyzer/pfTree");
  hasPFCandTree = (pfTree   != 0);
  if (hasPFCandTree) {
    pfTree->SetName("pfcand");
    if (tree == 0) tree = pfTree; else tree->AddFriend(pfTree);
    setupPFCandTree(pfTree,pfcand);
    //std::cout<<"SetupPFCand Tree : " <<std::endl;
  }

  if( pp && !hireco)trackTree  = (TTree*) inf->Get("ppTrack/trackTree");
  else trackTree     = (TTree*) inf->Get("anaTrack/trackTree");

  hasTrackTree  = (trackTree   != 0);
  if (hasTrackTree) {
    trackTree->SetName("track");
    if (tree == 0) tree = trackTree; else tree->AddFriend(trackTree);
    setupTrackTree(trackTree,track);
    //std::cout<<"SetupTrack Tree : " <<std::endl;
  }

  //! Jet trees
  for(int idir=0;idir<ndir;idir++){
    JetTree[idir] = (TTree*)inf->Get(Form("%sJetAnalyzer/t",calgo[idir]));
    hasJetTree[idir] = (JetTree[idir]!=0);
    if(hasJetTree[idir]){
      JetTree[idir]->SetName(Form("%s",calgo[idir]));
      if (tree == 0) tree = JetTree[idir]; else tree->AddFriend(JetTree[idir]);
      setupJetTree(JetTree[idir],mJets[idir]);      
      //std::cout<<Form("SetupJet Tree : %s",calgo[idir])<<"\t # of entries : "<<JetTree[idir]->GetEntries()<<std::endl;
    }
  }


  tupleTree     = (TTree*) inf->Get("tupel/EventTree");
  hasTupleTree  = (tupleTree   != 0);
  if (hasTupleTree) {
    tupleTree->SetName("tuple");
    if (tree == 0) tree = tupleTree; else tree->AddFriend(tupleTree);
    tuple = new EventTree("tupleObj",tupleTree);
    //setupTupleTree(tupleTree,tuple);
    //std::cout<<"SetupTuple Tree : " <<std::endl;
  }

  tree->SetMarkerStyle(20);
  std::cout<<"All Trees set up for " << name << std::endl;

  // Print the status of thre forest
  PrintStatus();

  //! Fill map for association of jet algo and the trees
  FillMap();
}

HiForest::~HiForest()
{
  delete [] mJets;
  EmptyMap();
}

Long64_t HiForest::GetEntry(int i)
{
  Long64_t nb=0;
  currentEvent = i;
  //std::cout << " Entry " << i <<std::endl;
  // get the entry of the available trees
  if (hasHltTree){      
    nb += hltTree      ->GetEntry(i);
    //std::cout << " \t Entry from  hlt " << std::endl;
  }
  if (hasSkimTree){
    nb += skimTree     ->GetEntry(i);
    //std::cout << " \t Entry from  skim " << std::endl;
  }
  if (hasPFCandTree){
    nb += pfTree       ->GetEntry(i);
    //std::cout << " \t Entry from PFCand " << std::endl;
  }
  if (hasTrackTree){
    nb += trackTree    ->GetEntry(i);
    //std::cout << " \t Entry from trackTree " << std::endl;
  }
  if (hasEvtTree){
    nb += evtTree      ->GetEntry(i);
    //std::cout << " \t Entry from  evtTree " << std::endl;
  }
  if (hasTupleTree){
    //std::cout << " \t ********************* Entry from  tuple " << std::endl;
    nb += tupleTree    ->GetEntry(i);
    //std::cout << " \t Entry from  tuple " << std::endl;
  }
  for(int idir=0;idir<ndir;idir++){
    if (hasJetTree[idir]){
      nb += JetTree[idir]->GetEntry(i);
      //std::cout << " \t Entry from  JetTree  " << calgo[idir] << std::endl;
    }
  }
  return nb;
}

Long64_t HiForest::GetEntries()
{
  // get the entries of the available trees
  return nEntries;
}

void HiForest::CheckTree(TTree *t,const char *title)
{
  Long64_t entries = t->GetEntries();
  if (nEntries==0) nEntries = entries;
  std::cout <<title<<"\t : "<<entries<<" entries loaded.";
  if (entries != nEntries) {
    std::cout <<" Inconsistent number of entries!!"<<std::endl;
  } else {
    std::cout <<std::endl;
  }
}
void HiForest::PrintStatus()
{
  if (hasHltTree)      CheckTree(hltTree,      "HltTree");
  if (hasSkimTree)     CheckTree(skimTree,     "SkimTree");
  if (hasPFCandTree)   CheckTree(pfTree,       "PFCandTree");
  if (hasTrackTree)    CheckTree(trackTree,    "TrackTree");
  if (hasEvtTree)      CheckTree(evtTree,      "EvtTree");
  if (hasTupleTree)    CheckTree(tupleTree,    "TupleTree");
  for(int i=0;i<ndir;i++){
    if(hasJetTree[i])CheckTree(JetTree[i],Form("%sJetTree",calgo[i]));
    //std::cout<<"Jet Tree : "<<calgo[i]<<std::endl;
  }
}
// ====================== Jet related functions ==================
int HiForest::GetNAlgo(){
  //return ndir;
  return JetContainer.size();
}
const char *HiForest::GetAlgoName(int i){
  return calgo[i];
}
Jets *HiForest::GetJet(int i){
  return &(mJets[i]);
}

void HiForest::FillMap(){
  std::cout <<std::endl;
  std::cout <<std::endl;
  for(int idir=0;idir<ndir;idir++){
    if(!hasJetTree[idir])continue;
    JetContainer.insert(JetPair(calgo[idir],idir));
  }
  std::cout<<" FillMap() : Size of the Jet Container : " <<JetContainer.size()<<std::endl;
}
Jets *HiForest::GetJetByAlgo(const char *algo){
  int iFound=-1;
  for (it = JetContainer.begin(); it != JetContainer.end(); ++it) {
    if(strcmp((*it).first,algo)==0){
      iFound=(*it).second;
      break;
    }
  }
  return &(mJets[iFound]);
}
void HiForest::SelectJetAlgo(const char *algo[],const int len){

  for(int i=0;i<ndir;i++){
    hasJetTree[i]=0;
    for(int j=0;j<len;j++){
      if(strcmp(algo[j],calgo[i])==0){
	hasJetTree[i]=1;
	CheckTree(JetTree[i],Form("%sJetTree",calgo[i]));
	std::cout<<"Jet Tree : "<<calgo[i]<<"\t "<<algo[j]<<std::endl;
      }
    }
  }
}
void HiForest::SelectBranches(const char *tname,const char *blst[],const int len){
  
  if(strcmp(tname,"JetTree")==0){
    for(int i=0;i<ndir;i++){
      //std::cout<<"Jet Tree status : "<<calgo[i]<<"\t"<<hasJetTree[i]<<std::endl;
      if(!hasJetTree[i])continue;
      JetTree[i]->SetBranchStatus("*",0);
      //std::cout<<calgo[i]<<"\t # of entries : "<<JetTree[i]->GetEntries()<<std::endl;
      for(int j=0;j<len;j++){  
	JetTree[i]->SetBranchStatus(blst[j],1);
      //std::cout<<"\t branches : "<<blst[j]<<"\t status : "<<JetTree[i]->GetBranchStatus(blst[j])<<std::endl;
      }  
    }
  }
  if(strcmp(tname,"evtTree")==0){
    if (hasEvtTree){
      evtTree->SetBranchStatus("*",0);
      for(int j=0;j<len;j++){        
	evtTree->SetBranchStatus(blst[j],1);
      }
    }else std::cout<<"eveTree not available : "<<std::endl;
  }
}
void HiForest::SelectTupleBranches(std::vector<std::string> sblist)
{
  std::cout <<"  SelectTupleBranches " << sblist.size() << std::endl;
  if (hasTupleTree){
    tupleTree->SetBranchStatus("*",0,0);  
    tupleTree->SetBranchStatus("EvtIsRealData",1);
    tupleTree->SetBranchStatus("EvtNum",1);
    tupleTree->SetBranchStatus("EvtRunNum",1);
    tupleTree->SetBranchStatus("EvtLumiNum",1);
    tupleTree->SetBranchStatus("EvtBxNum",1);
    tupleTree->SetBranchStatus("EvtVtxCnt",1);
    tupleTree->SetBranchStatus("EvtPuCnt",1);
    tupleTree->SetBranchStatus("EvtPuCntTruth",1);
    tupleTree->SetBranchStatus("EvtWeights",1);
    tupleTree->SetBranchStatus("TrigHlt",1);
    tupleTree->SetBranchStatus("TrigHltMu",1);
    tupleTree->SetBranchStatus("TrigHltDiMu",1);
    tupleTree->SetBranchStatus("TrigHltEl",1);
    tupleTree->SetBranchStatus("TrigHltDiEl",1);
    tupleTree->SetBranchStatus("TrigHltElMu",1);

    std::cout <<" I  ok switched  Evet and Trig one "  << std::endl;

    for(std::vector<std::string>::const_iterator istr = sblist.begin(); istr != sblist.end(); ++istr){
      std::string cont_type = (*istr); 
      if( cont_type == "Zmumu" ){
	tupleTree->SetBranchStatus("MuPt",1);
	tupleTree->SetBranchStatus("MuEta",1);
	tupleTree->SetBranchStatus("MuPhi",1);
	tupleTree->SetBranchStatus("MuE",1);      
	tupleTree->SetBranchStatus("MuId",1);
	tupleTree->SetBranchStatus("MuIdTight",1);
	tupleTree->SetBranchStatus("MuCh",1);
	tupleTree->SetBranchStatus("MuVtxZ",1);
	tupleTree->SetBranchStatus("MuDxy",1);
	tupleTree->SetBranchStatus("MuPfIso",1);
	tupleTree->SetBranchStatus("MuType",1);
	tupleTree->SetBranchStatus("MuIsoTkIsoAbs",1);
	tupleTree->SetBranchStatus("MuIsoTkIsoRel",1);
	tupleTree->SetBranchStatus("MuIsoCalAbs",1);
	tupleTree->SetBranchStatus("MuIsoCombRel",1);
	tupleTree->SetBranchStatus("MuTkNormChi2",1);
	tupleTree->SetBranchStatus("MuTkHitCnt",1);
	tupleTree->SetBranchStatus("MuMatchedStationCnt",1);
	tupleTree->SetBranchStatus("MuDz",1);
	tupleTree->SetBranchStatus("MuPixelHitCnt",1);
	tupleTree->SetBranchStatus("MuTkLayerCnt",1);
	tupleTree->SetBranchStatus("MuPfIsoChHad",1);
	tupleTree->SetBranchStatus("MuPfIsoNeutralHad",1);
	tupleTree->SetBranchStatus("MuPfIsoRawRel",1);
	tupleTree->SetBranchStatus("MuHltMatch",1);
	std::cout <<" II  ok Zmumu  "  << std::endl;
	continue;
      }
      else if( cont_type == "Zee" ){
	tupleTree->SetBranchStatus("ElPt",1);
	tupleTree->SetBranchStatus("ElEta",1);
	tupleTree->SetBranchStatus("ElEtaSc",1);
	tupleTree->SetBranchStatus("ElPhi",1);
	tupleTree->SetBranchStatus("ElE",1);
	tupleTree->SetBranchStatus("ElId",1);
	tupleTree->SetBranchStatus("ElCh",1);
	tupleTree->SetBranchStatus("ElMvaTrig",1);
	tupleTree->SetBranchStatus("ElMvaNonTrig",1);
	tupleTree->SetBranchStatus("ElMvaPresel",1);
	tupleTree->SetBranchStatus("ElDEtaTkScAtVtx",1);
	tupleTree->SetBranchStatus("ElDPhiTkScAtVtx",1);
	tupleTree->SetBranchStatus("ElHoE",1);
	tupleTree->SetBranchStatus("ElSigmaIetaIeta",1);
	tupleTree->SetBranchStatus("ElSigmaIetaIetaFull5x5",1);
	tupleTree->SetBranchStatus("ElEinvMinusPinv",1);
	tupleTree->SetBranchStatus("ElD0",1);
	tupleTree->SetBranchStatus("ElDz",1);
	tupleTree->SetBranchStatus("ElExpectedMissingInnerHitCnt",1);
	tupleTree->SetBranchStatus("ElPassConvVeto",1);
	tupleTree->SetBranchStatus("ElHltMatch",1);
	tupleTree->SetBranchStatus("ElPfIsoChHad",1);
	tupleTree->SetBranchStatus("ElPfIsoNeutralHad",1);
	tupleTree->SetBranchStatus("ElPfIsoIso",1);
	tupleTree->SetBranchStatus("ElPfIsoPuChHad",1);
	tupleTree->SetBranchStatus("ElPfIsoRaw",1);
	tupleTree->SetBranchStatus("ElPfIsoDbeta",1);
	tupleTree->SetBranchStatus("ElAEff",1);
	tupleTree->SetBranchStatus("charged",1);
	tupleTree->SetBranchStatus("photon",1);
	tupleTree->SetBranchStatus("neutral",1);
	tupleTree->SetBranchStatus("charged_Tom",1);
	tupleTree->SetBranchStatus("photon_Tom",1);
	tupleTree->SetBranchStatus("neutral_Tom",1);
	std::cout <<" III  ok Zee  "  << std::endl;
	continue;
      }
      else if( cont_type == "Jet" ){
	tupleTree->SetBranchStatus("JetAk04Pt",1);
	tupleTree->SetBranchStatus("JetAk04Eta",1);
	tupleTree->SetBranchStatus("JetAk04Phi",1);
	tupleTree->SetBranchStatus("JetAk04E",1);
	tupleTree->SetBranchStatus("JetAk04Id",1);
	tupleTree->SetBranchStatus("JetAk04PuId",1);
	tupleTree->SetBranchStatus("JetAk04PuMva",1);
	tupleTree->SetBranchStatus("JetAk04RawPt",1);
	tupleTree->SetBranchStatus("JetAk04RawE",1);
	tupleTree->SetBranchStatus("JetAk04HfHadE",1);
	tupleTree->SetBranchStatus("JetAk04HfEmE",1);
	tupleTree->SetBranchStatus("JetAk04ChHadFrac",1);
	tupleTree->SetBranchStatus("JetAk04NeutralHadAndHfFrac",1);
	tupleTree->SetBranchStatus("JetAk04ChEmFrac",1);
	tupleTree->SetBranchStatus("JetAk04NeutralEmFrac",1);
	tupleTree->SetBranchStatus("JetAk04ChMult",1);
	tupleTree->SetBranchStatus("JetAk04ConstCnt",1);
	tupleTree->SetBranchStatus("JetAk04JetBeta",1);
	tupleTree->SetBranchStatus("JetAk04JetBetaClassic",1);
	tupleTree->SetBranchStatus("JetAk04JetBetaStar",1);
	tupleTree->SetBranchStatus("JetAk04JetBetaStarClassic",1);
	tupleTree->SetBranchStatus("JetAk04BTagCsv",1);
	tupleTree->SetBranchStatus("JetAk04BTagCsvV1",1);
	tupleTree->SetBranchStatus("JetAk04BTagCsvSLV1",1);
	tupleTree->SetBranchStatus("JetAk04BDiscCisvV2",1);
	tupleTree->SetBranchStatus("JetAk04BDiscJp",1);
	tupleTree->SetBranchStatus("JetAk04BDiscBjp",1);
	tupleTree->SetBranchStatus("JetAk04BDiscTche",1);
	tupleTree->SetBranchStatus("JetAk04BDiscTchp",1);
	tupleTree->SetBranchStatus("JetAk04BDiscSsvhe",1);
	tupleTree->SetBranchStatus("JetAk04BDiscSsvhp",1);
	tupleTree->SetBranchStatus("JetAk04PartFlav",1);
	tupleTree->SetBranchStatus("JetAk04JecUncUp",1);
	tupleTree->SetBranchStatus("JetAk04JecUncDwn",1);
	tupleTree->SetBranchStatus("JetAk04ConstId",1);
	tupleTree->SetBranchStatus("JetAk04ConstPt",1);
	tupleTree->SetBranchStatus("JetAk04ConstEta",1);
	tupleTree->SetBranchStatus("JetAk04ConstPhi",1);
	tupleTree->SetBranchStatus("JetAk04ConstE",1);
	if( mc )tupleTree->SetBranchStatus("JetAk04GenJet",1);
	continue;
      }
      else if( cont_type == "Photon" ){
	tupleTree->SetBranchStatus("PhotPt",1);
	tupleTree->SetBranchStatus("PhotEta",1);
	tupleTree->SetBranchStatus("PhotPhi",1);
	tupleTree->SetBranchStatus("PhotScRawE",1);
	tupleTree->SetBranchStatus("PhotScEta",1);
	tupleTree->SetBranchStatus("PhotScPhi",1);
	tupleTree->SetBranchStatus("PhotIsoEcal",1);
	tupleTree->SetBranchStatus("PhotIsoHcal",1);
	tupleTree->SetBranchStatus("PhotIsoTk",1);
	tupleTree->SetBranchStatus("PhotPfIsoChHad",1);
	tupleTree->SetBranchStatus("PhotPfIsoNeutralHad",1);
	tupleTree->SetBranchStatus("PhotPfIsoPhot",1);
	tupleTree->SetBranchStatus("PhotPfIsoPuChHad",1);
	tupleTree->SetBranchStatus("PhotPfIsoEcalClus",1);
	tupleTree->SetBranchStatus("PhotPfIsoHcalClus",1);
	tupleTree->SetBranchStatus("PhotE3x3",1);
	tupleTree->SetBranchStatus("PhotE1x5",1);
	tupleTree->SetBranchStatus("PhotE2x5",1);
	tupleTree->SetBranchStatus("PhotE5x5",1);
	tupleTree->SetBranchStatus("PhotSigmaIetaIeta",1);
	tupleTree->SetBranchStatus("PhotEtaWidth",1);
	tupleTree->SetBranchStatus("PhotPhiWidth",1);
	tupleTree->SetBranchStatus("PhotHoE",1);
	tupleTree->SetBranchStatus("PhotId",1);
	tupleTree->SetBranchStatus("PhotHasPixelSeed",1);
	std::cout <<" V  ok Photon  "  << std::endl;
	continue;
      }
      else if( cont_type == "PhotonJet" ){
	tupleTree->SetBranchStatus("GLepDr01Pt",1);
	tupleTree->SetBranchStatus("GLepDr01Eta",1);
	tupleTree->SetBranchStatus("GLepDr01Phi",1);
	tupleTree->SetBranchStatus("GLepDr01E",1);
	tupleTree->SetBranchStatus("GLepDr01Id",1);
	tupleTree->SetBranchStatus("GLepDr01St",1);
	tupleTree->SetBranchStatus("GLepDr01MomId",1);
	tupleTree->SetBranchStatus("GLepBarePt",1);
	tupleTree->SetBranchStatus("GLepBareEta", 1);
	tupleTree->SetBranchStatus("GLepBarePhi", 1);
	tupleTree->SetBranchStatus("GLepBareE"  , 1);
	tupleTree->SetBranchStatus("GLepBareId" , 1);
	tupleTree->SetBranchStatus("GLepBareSt" , 1);
	tupleTree->SetBranchStatus("GLepBareMomId",1);
	tupleTree->SetBranchStatus("GLepSt3Pt",1);
	tupleTree->SetBranchStatus("GLepSt3Eta",1);
	tupleTree->SetBranchStatus("GLepSt3Phi",1);
	tupleTree->SetBranchStatus("GLepSt3E",1);
	tupleTree->SetBranchStatus("GLepSt3Id",1);
	tupleTree->SetBranchStatus("GLepSt3St",1);
	tupleTree->SetBranchStatus("GLepSt3Mother0Id",1);

	tupleTree->SetBranchStatus("GPhotSt",1);
	tupleTree->SetBranchStatus("GPhotIsoEDR03",1);
	tupleTree->SetBranchStatus("GPhotIsoEDR04",1);
	tupleTree->SetBranchStatus("GPhotIsoEDR05",1);
	tupleTree->SetBranchStatus("GPhotIsoSumPtDR03",1);
	tupleTree->SetBranchStatus("GPhotIsoSumPtDR04",1);
	tupleTree->SetBranchStatus("GPhotIsoSumPtDR05",1);
	tupleTree->SetBranchStatus("GLepClosePhotPt" ,1);
	tupleTree->SetBranchStatus("GLepClosePhotEta",1);
	tupleTree->SetBranchStatus("GLepClosePhotPhi",1);
	tupleTree->SetBranchStatus("GLepClosePhotE" ,1);
	tupleTree->SetBranchStatus("GLepClosePhotId",1);
	tupleTree->SetBranchStatus("GLepClosePhotMother0Id",1);
	tupleTree->SetBranchStatus("GLepClosePhotMotherCnt",1);
	tupleTree->SetBranchStatus("GLepClosePhotSt",1);
	tupleTree->SetBranchStatus("GJetAk04Pt" ,1);
	tupleTree->SetBranchStatus("GJetAk04Eta",1);
	tupleTree->SetBranchStatus("GJetAk04Phi",1);
	tupleTree->SetBranchStatus("GJetAk04E"  ,1);
	tupleTree->SetBranchStatus("GJetAk04ChFrac",1);
	tupleTree->SetBranchStatus("GJetAk04ConstCnt",1);
	tupleTree->SetBranchStatus("GJetAk04ConstId",1);
	tupleTree->SetBranchStatus("GJetAk04ConstPt",1);
	tupleTree->SetBranchStatus("GJetAk04ConstEta",1);
	tupleTree->SetBranchStatus("GJetAk04ConstPhi",1);
	tupleTree->SetBranchStatus("GJetAk04ConstE",1);
	tupleTree->SetBranchStatus("GPdfId1",1);
	tupleTree->SetBranchStatus("GPdfId2",1);
	tupleTree->SetBranchStatus("GPdfx1" ,1);
	tupleTree->SetBranchStatus("GPdfx2" ,1);
	tupleTree->SetBranchStatus("GPdfScale",1);
	tupleTree->SetBranchStatus("GBinningValue",1);
	tupleTree->SetBranchStatus("GNup",1);
	std::cout <<" VI  ok PhotonJet  "  << std::endl;
	continue;
      }
      else if( cont_type == "MET" ){
	tupleTree->SetBranchStatus("METPt",1);
	tupleTree->SetBranchStatus("METPx",1);
	tupleTree->SetBranchStatus("METPy",1);
	tupleTree->SetBranchStatus("METPz",1);
	tupleTree->SetBranchStatus("METE",1);
	tupleTree->SetBranchStatus("METsigx2",1);
	tupleTree->SetBranchStatus("METsigxy",1);
	tupleTree->SetBranchStatus("METsigy2",1);
	tupleTree->SetBranchStatus("METsig",1);
	std::cout <<" VII  ok MET  "  << std::endl;
	continue;
      }else std::cout <<" No option found "    <<std::endl;
    } 
  }else  std::cout<<"tupleTree not available : "<<std::endl;
}

void HiForest::EmptyMap(){
  if (JetContainer.empty()) return;
  JetContainer.clear();
}
// ====================== Event Utilities ========================
