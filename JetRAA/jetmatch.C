
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
#include <map>
#include <set>
#include <utility>

using namespace std;

#ifdef _MAKECINT_
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#endif

void AddInputFiles(TChain */*ch*/, string /*iname*/, string /*inputTree*/);

int GetPhiBin(float /*phi*/);
int GetEtaBin(float /*eta*/);
int GetPtBin(float /*pt*/);
int GetCentBin(int /*hiBin*/);
float delphi(float /*phi1*/, float /*phi2*/);
float deltaR(float /*eta1*/, float /*phi1*/, 
	      float /*eta2*/, float /*phi2*/);

struct Jet{
  int id;
  float pt;
  float eta;
  float phi;
};

bool compare_pt(Jet jet1, Jet jet2);
bool compare_pt(Jet jet1, Jet jet2){
  return jet1.pt > jet2.pt;
}

typedef std::pair< Jet, Jet > CaloPFJetPair;
struct CompareMatchedJets {
  bool operator()(const CaloPFJetPair &A1, const CaloPFJetPair &A2){
    Jet pf1 = A1.first; //! PFJet   1st pair
    Jet cj1 = A1.second;//! CaloJet 1st pair

    Jet pf2 = A2.first; //! PFJet   2nd pair
    Jet cj2 = A2.second;//! CaloJet 2nd pair

    float delr1 = deltaR(pf1.eta, pf1.phi, cj1.eta, cj1.phi);
    float delr2 = deltaR(pf2.eta, pf2.phi, cj2.eta, cj2.phi);

    //float delpt1 = fabs(pf1.pt - cj1.pt);
    //float delpt2 = fabs(pf2.pt - cj2.pt);

    //return ((delpt1 < delpt2) && (delr1 < delr2) && (pf1.pt > pf2.pt) && cj1.pt > cj2.pt);
    
    return ((delr1 < delr2) && (pf1.pt > pf2.pt));

    //return (delr1 < delr2);
  }
};

typedef std::multiset< CaloPFJetPair, CompareMatchedJets > CaloPFMatchedJets;
typedef std::multiset< CaloPFJetPair >::iterator CPFItr;
//typedef std::multiset< Jet >::value_type MatchedJet;



const int ncen=7;
const char *cdir [ncen] = {"05","510","1030","3050","5070","7090","90100"};
const char *ccent[ncen] = {"0-5%","5-10%","10-30%","30-50%","50-70%","70-90%","90-100%"};

const int ncand=5;
const char *ccand[ncand] = {"h^{#pm}","e^{#pm}","#mu^{#pm}","#gamma","h0"};

//! constants
int iYear=2015;
const double pi=acos(-1.);

const double ketacut=2.0;
const double kptrawcut =0.0;
const double kptrecocut=30.0;
const double kdelrcut=0.2;

//! pt binning
double ptbins[] ={24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 
		  114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 
		  395, 430, 468, 507, 548, 592, 638, 686, 1000};
const int nbins = sizeof(ptbins)/sizeof(double) - 1;

const double ptbins_wide[]={30,50,80,120,200,340,548};
const int nbins_wide = sizeof(ptbins_wide)/sizeof(double) - 1;

const double etabins[] = {-2.000, -1.4000, -0.4500, 0.000, 0.4500, 1.400, 2.000};
const int neta = sizeof(etabins)/sizeof(double) - 1;
const double phibins[] = {-3.141,-2.100,-1.500,-0.800,-0.300, 
 			  0.300,0.800,1.500,2.100,3.141
};
const int nphi = sizeof(phibins)/sizeof(double) - 1;

//! PuPtMin values 
// akPu1PFJets.puPtMin = 10
// akPu2PFJets.puPtMin = 10
// akPu3PFJets.puPtMin = 15
// akPu4PFJets.puPtMin = 20
// akPu5PFJets.puPtMin = 25
// akPu6PFJets.puPtMin = 30
// akPu6PFJets.puPtMin = 35


TStopwatch timer;
const char *ksp="pbpb";
int jetmatch(std::string kAlgName="akPu3",
		    //std::string finame="HiForest_jet55or65or80_JetRAA_v1_lumi1_Part18.root",
		    std::string fileList="/mnt/hadoop/cms/store/user/belt/HiForest_jet55or65or80_JetRAA_v1_final/HiForest_jet55or65or80_JetRAA_v1_lumi9_Part8.root,/mnt/hadoop/cms/store/user/belt/HiForest_jet55or65or80_JetRAA_v1_final/HiForest_jet55or65or80_JetRAA_v1_lumi9_Part80.root,/mnt/hadoop/cms/store/user/belt/HiForest_jet55or65or80_JetRAA_v1_final/HiForest_jet55or65or80_JetRAA_v1_lumi9_Part81.root",
		    std::string foname="outputhisto_data.root" 
		    )
{
  
  timer.Start();


  //tr_jet = (TTree*)fin->Get("akPu3PFJetAnalyzer/t");
  TChain *tch_pfjet = new TChain(Form("%sPFJetAnalyzer/t",kAlgName.c_str()));
  AddInputFiles(tch_pfjet,fileList,Form("%sPFJetAnalyzer/t",kAlgName.c_str()));
  cout <<" # of events in PFJet Tree : " <<  tch_pfjet->GetEntries() <<endl;

  //tr_jet = (TTree*)fin->Get("akPu3CaloJetAnalyzer/t");
  TChain *tch_calojet = new TChain(Form("%sCaloJetAnalyzer/t",kAlgName.c_str()));
  AddInputFiles(tch_calojet,fileList,Form("%sCaloJetAnalyzer/t",kAlgName.c_str()));
  cout <<" # of events in CaloJet Tree : " <<  tch_calojet->GetEntries() <<endl;

  //tr_ev = (TTree*)fin->Get("hiEvtAnalyzer/HiTree");
  TChain *tch_ev = new TChain("hiEvtAnalyzer/HiTree");
  AddInputFiles(tch_ev,fileList,"hiEvtAnalyzer/HiTree");
  cout <<" # of events in Event Tree : " <<  tch_ev->GetEntries() <<endl;

  //tr_hlt = (TTree*)fin->Get("hltanalysis/HltTree");
  TChain *tch_hlt = new TChain("hltanalysis/HltTree");
  AddInputFiles(tch_hlt,fileList,"hltanalysis/HltTree");
  cout <<" # of events in HLT Tree : " <<  tch_hlt->GetEntries() <<endl;

  //tr_skim = (TTree*)fin->Get("skimanalysis/HltTree");
  TChain *tch_skim = new TChain("skimanalysis/HltTree");
  AddInputFiles(tch_skim,fileList,"skimanalysis/HltTree");
  cout <<" # of events in Skim Tree : " <<  tch_skim->GetEntries() <<endl;

  //tr_trobj = (TTree*)fin->Get("hltobject/jetObjTree");  
  // TChain *tch_trgobj = new TChain("hltobject/jetObjTree");
  // AddInputFiles(tch_trgobj,fileList,"hltobject/jetObjTree");
  // cout <<" # of events in TrigObj Tree : " <<  tch_trgobj->GetEntries() <<endl;
  // cout <<endl;


  //! Event Tree
  int run_value;
  int evt_value;
  int lumi_value;
  int hiNpix;
  int hiBin;
  float vz;

  tch_ev->SetBranchAddress("run",&run_value);  
  tch_ev->SetBranchAddress("evt",&evt_value);  
  tch_ev->SetBranchAddress("lumi",&lumi_value);  
  tch_ev->SetBranchAddress("hiBin",&hiBin);
  tch_ev->SetBranchAddress("hiNpix",&hiNpix);  
  tch_ev->SetBranchAddress("vz",&vz);


  //! HLT tree
  int jet55;
  int jet55_prescl;
  int jet65;
  int jet65_prescl;
  int jet80;
  int jet80_prescl;
  tch_hlt->SetBranchAddress("HLT_HIJet55_v1",&jet55);
  tch_hlt->SetBranchAddress("HLT_HIJet55_v1_Prescl",&jet55_prescl);
  tch_hlt->SetBranchAddress("HLT_HIJet65_v1",&jet65);
  tch_hlt->SetBranchAddress("HLT_HIJet65_v1_Prescl",&jet65_prescl);
  tch_hlt->SetBranchAddress("HLT_HIJet80_v1",&jet80);
  tch_hlt->SetBranchAddress("HLT_HIJet80_v1_Prescl",&jet80_prescl);

  // int L1_sj36_1;
  // int L1_sj36_p_1;
  // int L1_sj52_1;
  // int L1_sj52_p_1;
  // tch_hlt->SetBranchAddress("L1_SingleJet36_BptxAND",&L1_sj36_1);
  // tch_hlt->SetBranchAddress("L1_SingleJet36_BptxAND_Prescl",&L1_sj36_p_1);
  // tch_hlt->SetBranchAddress("L1_SingleJet52_BptxAND",&L1_sj52_1);
  // tch_hlt->SetBranchAddress("L1_SingleJet52_BptxAND_Prescl",&L1_sj52_p_1);


  //! Skim Tree
  int pcollisionEventSelection;
  int pHBHENoiseFilter;
  tch_skim->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
  tch_skim->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);

  //! Trigger object tree
  // float trgObj_id;
  // float trgObj_pt;
  // float trgObj_eta;
  // float trgObj_phi;
  // tch_trgobj->SetBranchAddress("id",&trgObj_id);
  // tch_trgobj->SetBranchAddress("pt",&trgObj_pt);
  // tch_trgobj->SetBranchAddress("eta",&trgObj_eta);
  // tch_trgobj->SetBranchAddress("phi",&trgObj_phi);


  //! CaloJet tree
  //Declaration of leaves types
  int   nref_calo;
  float jtpt_calo[1000];
  float rawpt_calo[1000];
  float jtpu_calo[1000];
  float jteta_calo[1000];
  float jtphi_calo[1000];
  float hcalSum_calo[1000];
  float ecalSum_calo[1000];

  tch_calojet->SetBranchAddress("nref",&nref_calo);
  tch_calojet->SetBranchAddress("rawpt",rawpt_calo);
  tch_calojet->SetBranchAddress("jtpt" ,jtpt_calo);
  tch_calojet->SetBranchAddress("jtpu" ,jtpu_calo);
  tch_calojet->SetBranchAddress("jteta",jteta_calo);
  tch_calojet->SetBranchAddress("jtphi",jtphi_calo);
  tch_calojet->SetBranchAddress("hcalSum",hcalSum_calo);
  tch_calojet->SetBranchAddress("ecalSum",ecalSum_calo);

  //! PFJet tree
  //Declaration of leaves types
  int   nref;
  float jtpt[1000];
  float rawpt[1000];
  float jteta[1000];
  float jtpu [1000];
  float jtphi[1000];
  float neutralSum[1000];
  float chargedSum[1000];
  float photonSum[1000];
  float eleSum[1000];
  float muonSum[1000];
  float neutralMax[1000];
  float chargedMax[1000];
  float photonMax[1000];
  float eleMax[1000];
  float muonMax[1000];
  float hcalSum_pf[1000];
  float ecalSum_pf[1000];


  tch_pfjet->SetBranchAddress("nref",&nref);
  tch_pfjet->SetBranchAddress("rawpt",rawpt);
  tch_pfjet->SetBranchAddress("jtpt" ,jtpt);
  tch_pfjet->SetBranchAddress("jtpu" ,jtpu);
  tch_pfjet->SetBranchAddress("jteta",jteta);
  tch_pfjet->SetBranchAddress("jtphi",jtphi);
  tch_pfjet->SetBranchAddress("neutralSum",neutralSum);
  tch_pfjet->SetBranchAddress("chargedSum",chargedSum);
  tch_pfjet->SetBranchAddress("photonSum",photonSum);
  tch_pfjet->SetBranchAddress("eSum",eleSum);
  tch_pfjet->SetBranchAddress("muSum",muonSum);
  tch_pfjet->SetBranchAddress("neutralMax",neutralMax);
  tch_pfjet->SetBranchAddress("chargedMax",chargedMax);
  tch_pfjet->SetBranchAddress("photonMax",photonMax);
  tch_pfjet->SetBranchAddress("eMax",eleMax);
  tch_pfjet->SetBranchAddress("muMax",muonMax);
  tch_pfjet->SetBranchAddress("hcalSum",hcalSum_pf);
  tch_pfjet->SetBranchAddress("ecalSum",ecalSum_pf);
  
  tch_pfjet->AddFriend(tch_ev);
  tch_pfjet->AddFriend(tch_hlt);
  tch_pfjet->AddFriend(tch_skim);
  //tch_pfjet->AddFriend(tch_trgobj);
  tch_pfjet->AddFriend(tch_calojet);

  //! Disable branches 
  //! Jet Tree
  tch_pfjet->SetBranchStatus("*",0,0);
  tch_pfjet->SetBranchStatus("nref" ,1,0);
  tch_pfjet->SetBranchStatus("rawpt",1,0);
  tch_pfjet->SetBranchStatus("jtpt" ,1,0);
  tch_pfjet->SetBranchStatus("jtpu" ,1,0);
  tch_pfjet->SetBranchStatus("jteta",1,0);
  tch_pfjet->SetBranchStatus("jtphi",1,0);
  tch_pfjet->SetBranchStatus("neutralSum",1,0);
  tch_pfjet->SetBranchStatus("chargedSum",1,0);
  tch_pfjet->SetBranchStatus("photonSum",1,0);
  tch_pfjet->SetBranchStatus("eSum",1,0);
  tch_pfjet->SetBranchStatus("muSum",1,0);
  tch_pfjet->SetBranchStatus("neutralMax",1,0);
  tch_pfjet->SetBranchStatus("chargedMax",1,0);
  tch_pfjet->SetBranchStatus("photonMax",1,0);
  tch_pfjet->SetBranchStatus("eMax",1,0);
  tch_pfjet->SetBranchStatus("muMax",1,0);
  tch_pfjet->SetBranchStatus("hcalSum",1,0);
  tch_pfjet->SetBranchStatus("ecalSum",1,0);

  tch_pfjet->SetBranchStatus("run",1,0);
  tch_pfjet->SetBranchStatus("evt",1,0);
  tch_pfjet->SetBranchStatus("lumi",1,0);
  tch_pfjet->SetBranchStatus("hiNpix",1,0);
  tch_pfjet->SetBranchStatus("hiBin",1,0);
  tch_pfjet->SetBranchStatus("vz",1,0);

  // tch_pfjet->SetBranchStatus("id",1,0);
  // tch_pfjet->SetBranchStatus("pt",1,0);
  // tch_pfjet->SetBranchStatus("eta",1,0);
  // tch_pfjet->SetBranchStatus("phi",1,0);


  tch_pfjet->SetBranchStatus("HLT_HIJet55_v1",1,0);
  tch_pfjet->SetBranchStatus("HLT_HIJet55_v1_Prescl",1,0);
  tch_pfjet->SetBranchStatus("HLT_HIJet65_v1",1,0);
  tch_pfjet->SetBranchStatus("HLT_HIJet65_v1_Prescl",1,0);
  tch_pfjet->SetBranchStatus("HLT_HIJet80_v1",1,0);
  tch_pfjet->SetBranchStatus("HLT_HIJet80_v1_Prescl",1,0);
  // tch_pfjet->SetBranchStatus("L1_SingleJet36_BptxAND",1,0);
  // tch_pfjet->SetBranchStatus("L1_SingleJet36_BptxAND_Prescl",1,0);
  // tch_pfjet->SetBranchStatus("L1_SingleJet52_BptxAND",1,0);
  // tch_pfjet->SetBranchStatus("L1_SingleJet52_BptxAND_Prescl",1,0);

  tch_pfjet->SetBranchStatus("pcollisionEventSelection",1,0);
  tch_pfjet->SetBranchStatus("pHBHENoiseFilter",1,0);


  std::string outdir="";
  std::string outname=outdir+foname;
  TFile *fout = new TFile(outname.c_str(),"RECREATE");


  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"**************************************************** "<<std::endl;
  std::cout<<Form("PbPbData : %s %s ",ksp, kAlgName.c_str())<<std::endl;
  std::cout<<Form("Outfile : %s",outname.c_str())<<std::endl;
  std::cout<<Form("eta cut for %0.3f ",ketacut)<<std::endl;
  std::cout<<"**************************************************** "<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;


  //! 
  //! Define histograms here
  fout->mkdir(Form("%sJetAnalyzer",kAlgName.c_str()));
  fout->cd(Form("%sJetAnalyzer"   ,kAlgName.c_str()));


  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  Float_t deltar;
  Float_t calopt, pfpt;
  Float_t calorawpt, pfrawpt;
  Float_t caloeta, pfeta;
  Float_t calophi, pfphi;
  Float_t calopu, pfpu;
  Float_t chMax, phMax, neMax,  chSum, phSum, neSum,  eSum, muSum, muMax, eMax, hcalSum, ecalSum;
  /*Int_t hiBin jet80, jet80_prescl, jet65, jet65_prescl, jet55, jet55_prescl;*/
  /*Int_t evt_value, run_value, lumi_value;*/

  TTree* matchJets = new TTree("matchJets","Ntuple containing important information about matched jets");
  matchJets->Branch("hiBin",&hiBin,"hiBin/I");
  matchJets->Branch("run_value",&run_value,"run_value/I");
  matchJets->Branch("evt_value",&evt_value,"evt_value/I");
  matchJets->Branch("lumi_value",&lumi_value,"lumi_value/I");
  matchJets->Branch("jet80",&jet80,"jet80/I"); matchJets->Branch("jet80_prescl",&jet80_prescl,"jet80_prescl/I");
  matchJets->Branch("jet65",&jet65,"jet65/I"); matchJets->Branch("jet65_prescl",&jet65_prescl,"jet65_prescl/I");
  matchJets->Branch("jet55",&jet55,"jet55/I"); matchJets->Branch("jet55_prescl",&jet55_prescl,"jet55_prescl/I");
  matchJets->Branch("vz",&vz,"vz/F");
  matchJets->Branch("deltar",&deltar,"deltar/F"); 
  matchJets->Branch("calopt",&calopt,"calopt/F");    matchJets->Branch("pfpt",&pfpt,"pfpt/F");       
  matchJets->Branch("calorawpt",&calorawpt,"calorawpt/F");    matchJets->Branch("pfrawpt",&pfrawpt,"pfrawpt/F");       
  matchJets->Branch("caloeta",&caloeta,"caloeta/F"); matchJets->Branch("pfeta",&pfeta,"pfeta/F");       
  matchJets->Branch("calophi",&calophi,"calophi/F"); matchJets->Branch("pfphi",&pfphi,"pfphi/F");       
  matchJets->Branch("calopu",&calopu,"calopu/F");    matchJets->Branch("pfpu",&pfpu,"pfpu/F");       
  matchJets->Branch("chMax",&chMax,"chMax/F");       matchJets->Branch("chSum",&chSum,"chSum/F");    
  matchJets->Branch("phMax",&phMax,"phMax/F");       matchJets->Branch("phSum",&phSum,"phSum/F");
  matchJets->Branch("neMax",&neMax,"neMax/F");       matchJets->Branch("neSum",&neSum,"neSum/F");
  matchJets->Branch("muMax",&muMax,"muMax/F");       matchJets->Branch("muSum",&muSum,"muSum/F");
  matchJets->Branch("eMax",&eMax,"eMax/F");          matchJets->Branch("eSum",&eSum,"eSum/F");
  matchJets->Branch("hcalSum",&hcalSum,"hcalSum/F"); matchJets->Branch("ecalSum",&ecalSum,"ecalSum/F");


  TTree* unmatchPFJets = new TTree("unmatchPFJets","Ntuple containing important information about unmatched PF jets");
  unmatchPFJets->Branch("hiBin",&hiBin,"hiBin/I");
  unmatchPFJets->Branch("run_value",&run_value,"run_value/I");
  unmatchPFJets->Branch("evt_value",&evt_value,"evt_value/I");
  unmatchPFJets->Branch("lumi_value",&lumi_value,"lumi_value/I");
  unmatchPFJets->Branch("jet80",&jet80,"jet80/I"); unmatchPFJets->Branch("jet80_prescl",&jet80_prescl,"jet80_prescl/I");
  unmatchPFJets->Branch("jet65",&jet65,"jet65/I"); unmatchPFJets->Branch("jet65_prescl",&jet65_prescl,"jet65_prescl/I");
  unmatchPFJets->Branch("jet55",&jet55,"jet55/I"); unmatchPFJets->Branch("jet55_prescl",&jet55_prescl,"jet55_prescl/I");
  unmatchPFJets->Branch("vz",&vz,"vz/F");
  unmatchPFJets->Branch("pfpt",&pfpt,"pfpt/F");       
  unmatchPFJets->Branch("pfrawpt",&pfrawpt,"pfrawpt/F");       
  unmatchPFJets->Branch("pfeta",&pfeta,"pfeta/F");       
  unmatchPFJets->Branch("pfphi",&pfphi,"pfphi/F");       
  unmatchPFJets->Branch("pfpu",&pfpu,"pfpu/F");       
  unmatchPFJets->Branch("chMax",&chMax,"chMax/F");       unmatchPFJets->Branch("chSum",&chSum,"chSum/F");    
  unmatchPFJets->Branch("phMax",&phMax,"phMax/F");       unmatchPFJets->Branch("phSum",&phSum,"phSum/F");
  unmatchPFJets->Branch("neMax",&neMax,"neMax/F");       unmatchPFJets->Branch("neSum",&neSum,"neSum/F");
  unmatchPFJets->Branch("muMax",&muMax,"muMax/F");       unmatchPFJets->Branch("muSum",&muSum,"muSum/F");
  unmatchPFJets->Branch("eMax",&eMax,"eMax/F");          unmatchPFJets->Branch("eSum",&eSum,"eSum/F");

  TTree* unmatchCaloJets = new TTree("unmatchCaloJets","Ntuple containing important information about unmatched Calo jets");
  unmatchCaloJets->Branch("hiBin",&hiBin,"hiBin/I");
  unmatchCaloJets->Branch("run_value",&run_value,"run_value/I");
  unmatchCaloJets->Branch("evt_value",&evt_value,"evt_value/I");
  unmatchCaloJets->Branch("lumi_value",&lumi_value,"lumi_value/I");
  unmatchCaloJets->Branch("jet80",&jet80,"jet80/I"); unmatchCaloJets->Branch("jet80_prescl",&jet80_prescl,"jet80_prescl/I");
  unmatchCaloJets->Branch("jet65",&jet65,"jet65/I"); unmatchCaloJets->Branch("jet65_prescl",&jet65_prescl,"jet65_prescl/I");
  unmatchCaloJets->Branch("jet55",&jet55,"jet55/I"); unmatchCaloJets->Branch("jet55_prescl",&jet55_prescl,"jet55_prescl/I");
  unmatchCaloJets->Branch("vz",&vz,"vz/F");
  unmatchCaloJets->Branch("calopt",&calopt,"calopt/F");
  unmatchCaloJets->Branch("calorawpt",&calorawpt,"calorawpt/F");
  unmatchCaloJets->Branch("caloeta",&caloeta,"caloeta/F"); 
  unmatchCaloJets->Branch("calophi",&calophi,"calophi/F"); 
  unmatchCaloJets->Branch("calopu",&calopu,"calopu/F");    
  unmatchCaloJets->Branch("hcalSum",&hcalSum,"hcalSum/F"); unmatchCaloJets->Branch("ecalSum",&ecalSum,"ecalSum/F");

  TH1F *hEvents_Total = new TH1F("hEvents_Total","Total # of events ",10,0.,1.);
  hEvents_Total->Sumw2();
  TH1F *hEvents_pCollEvent = new TH1F("hEvents_pCollEvent","# of events ",10,0.,1.);
  hEvents_pCollEvent->Sumw2();
  TH1F *hEvents_pHBHENoise = new TH1F("hEvents_pHBHENoise","# of events ",10,0.,1.);
  hEvents_pHBHENoise->Sumw2();
  TH1F *hEvents_Vzcut = new TH1F("hEvents_Vzcut","# of events ",10,0.,1.);
  hEvents_Vzcut->Sumw2();
  TH1F *hEvents_Cent = new TH1F("hEvents_Cent","Cent # of events ",10,0.,1.);
  hEvents_Cent->Sumw2();
  fout->cd("../");
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout<<"Initialized the histograms " <<std::endl;

  Long64_t nbytes=0;
  Long64_t nentries = tch_pfjet->GetEntries();
  std::cout<<Form("# of entries in TTree for %s %s : ",kAlgName.c_str(),ksp)<<nentries<<std::endl;

  for (Long64_t i=0; i<nentries;i++) {
    nbytes += tch_pfjet->GetEntry(i);

    //if(i==500)break;
    
    float rndm=gRandom->Rndm();

    hEvents_Total->Fill(rndm);

    if(pcollisionEventSelection)hEvents_pCollEvent->Fill(rndm);
    if(pcollisionEventSelection && pHBHENoiseFilter)hEvents_pHBHENoise->Fill(rndm);
    if(pcollisionEventSelection && pHBHENoiseFilter && fabs(vz)<15. )hEvents_Vzcut->Fill(rndm);
    if(pcollisionEventSelection==0 || pHBHENoiseFilter==0 || fabs(vz) > 15.)continue;

    //if(!jet55_1 || !jet65_1 || !jet80_1)continue;
    
    int iCent = GetCentBin(hiBin);
    if(iCent<0 || iCent>=ncen)continue;

    //! SuperNovae events
    int lowjetCounter=0;
    int jetCounter=0;
    for(int g = 0;g<nref;g++){
      if( fabs(jteta[g]) < 2. && jtpt[g]>=kptrecocut ){ //to select inside
	lowjetCounter++;
      }
      if( fabs(jteta[g]) < 2. && jtpt[g]>=50. ){ //to select inside
	jetCounter++;
      }//eta selection cut
    }// jet loop

    // apply the correct supernova selection cut rejection here:
    if( hiNpix > 38000 - 500*jetCounter )continue;
    if( nref==0 && nref_calo==0 )continue;
    if( lowjetCounter == 0 )continue;
    
    hEvents_Cent->Fill(rndm);

    //std::cout<<" ********** Event # " <<i<<"\t vz  : "<<vz<<"\t hiBin : "<<hiBin<<"\t  # of PF jets "<<nref<<" # of calojets  "<<nref_calo<<std::endl;

    if(i%50000==0)std::cout<<" ********** Event # " <<i<<"\t vz  : "<<vz<<"\t hiBin : "<<hiBin<<"\t  # of PF jets "<<nref<< "  # of Calo jets " <<nref_calo<<std::endl;

    int pfjets=0;
    int calojets=0;

    std::vector <Jet> vPFJets, vCaloJets;
    std::vector <int> pfid(nref), caloid(nref_calo);

    for(int pj=0; pj<nref; pj++){ //! PFjet
      
      if( rawpt[pj] < kptrawcut || jtpt[pj] < kptrecocut ) continue;
      if( fabs(jteta[pj]) > ketacut ) continue;
      //if( (eMax[pj]/jtpt[pj])>=0.6 || (chMax[pj]/jtpt[pj])<=0.02 )continue;
      
      Jet pfj;
      pfj.id  = pj;
      pfj.eta = jteta[pj];
      pfj.phi = jtphi[pj];
      pfj.pt  = jtpt [pj];

      vPFJets.push_back(pfj);
      pfjets++;
    }    

    for(int cj=0; cj<nref_calo; cj++){ //! CaloJet
      
      if( rawpt_calo[cj] < kptrawcut || jtpt_calo[cj] < kptrecocut) continue;
      if( fabs(jteta_calo[cj]) > ketacut ) continue;
      
      Jet clj;
      clj.id  = cj;
      clj.eta = jteta_calo[cj];
      clj.phi = jtphi_calo[cj];
      clj.pt  = jtpt_calo[cj];

      vCaloJets.push_back(clj);
      calojets++;
    }//! calo jet loop
    

    bool onlyCalo   = (pfjets==0 && calojets >0) ? true : false;
    bool onlyPF     = (pfjets>0  && calojets==0) ? true : false;
    bool bothPFCalo = (pfjets>0  && calojets >0) ? true : false;

    int matchedJets=0;
    int unmatchedPFJets=0;
    int unmatchedCaloJets=0;

    std::vector < Jet >::const_iterator iJet, jJet;

    if( onlyPF ){
      
      for(iJet = vPFJets.begin(); iJet != vPFJets.end(); ++iJet){ //! PFjet

	int pj = (*iJet).id; 
	 	
    	deltar = -999;
    	calopt = -999;
    	pfpt   = jtpt [pj];
    	pfrawpt= rawpt[pj];
    	pfeta  = jteta[pj];
    	pfphi  = jtphi[pj];
    	pfpu   = jtpu [pj];
	
    	chMax  = chargedMax[pj];
    	phMax  = photonMax [pj];
    	neMax  = neutralMax[pj];
    	muMax  = muonMax   [pj];
    	eMax   = eleMax    [pj];
	
    	chSum  = chargedSum[pj];
    	phSum  = photonSum [pj];
    	neSum  = neutralSum[pj];
    	muSum  = muonSum   [pj];
    	eSum   = eleSum    [pj];
	
    	hcalSum = hcalSum_pf[pj];
    	ecalSum = ecalSum_pf[pj];

    	unmatchedPFJets++;
    	unmatchPFJets->Fill();
      }
    }
    else if( onlyCalo ){
      for(iJet = vCaloJets.begin(); iJet != vCaloJets.end(); ++iJet){ //! Calojet

	int cj = (*iJet).id; 	

    	deltar    = -999;
    	calopt    = jtpt_calo [cj];
    	pfpt      = -999;
    	calorawpt = rawpt_calo[cj];
    	caloeta   = jteta_calo[cj];
    	calophi   = jtphi_calo[cj];
    	calopu    = jtpu_calo [cj];
	
    	hcalSum   = hcalSum_calo[cj];
    	ecalSum   = ecalSum_calo[cj];
	
    	unmatchedCaloJets++;
    	unmatchCaloJets->Fill();
      }
    }
    else if( bothPFCalo ){

      CaloPFMatchedJets mCaloPFMatchedJets;
      for(iJet = vPFJets.begin(); iJet != vPFJets.end(); ++iJet){ //! PFjet

	for(jJet = vCaloJets.begin(); jJet != vCaloJets.end(); ++jJet){ //! Calojet      

	  mCaloPFMatchedJets.insert(std::make_pair(*iJet,*jJet));

	}//! calo jet loop
      }//! PF jet loop
    
      CPFItr itr;
      //! Matched jets (PF jet matched to Calo jet)
      for(itr = mCaloPFMatchedJets.begin(); itr != mCaloPFMatchedJets.end(); ++itr){
	CaloPFJetPair jetpair = (*itr);
	Jet pfj = jetpair.first;
	Jet clj = jetpair.second;
	float delr  = deltaR(pfj.eta, pfj.phi, clj.eta, clj.phi);
	//float delpt = fabs(pfj.pt - clj.pt);
	if( delr < kdelrcut && pfid[pfj.id]==0 ){
	
	  deltar = delr;
	
	  calopt    = jtpt_calo [clj.id];
	  pfpt      = jtpt      [pfj.id];
	  calorawpt = rawpt_calo[clj.id];
	  pfrawpt   = rawpt     [pfj.id];
	  caloeta   = jteta_calo[clj.id];
	  pfeta     = jteta     [pfj.id];
	  calophi   = jtphi_calo[clj.id];
	  pfphi     = jtphi     [pfj.id];
	  calopu    = jtpu_calo [clj.id];
	  pfpu      = jtpu      [pfj.id];
	
	  chMax  = chargedMax[pfj.id];
	  phMax  = photonMax [pfj.id];
	  neMax  = neutralMax[pfj.id];
	  muMax  = muonMax   [pfj.id];
	  eMax   = eleMax    [pfj.id];
	
	  chSum  = chargedSum[pfj.id];
	  phSum  = photonSum [pfj.id];
	  neSum  = neutralSum[pfj.id];
	  muSum  = muonSum   [pfj.id];
	  eSum   = eleSum    [pfj.id];
	
	  hcalSum = hcalSum_pf[clj.id];
	  ecalSum = ecalSum_pf[clj.id];
	  
	  matchedJets++;	  
	  matchJets->Fill();
	
	  // cout <<"\t *** Matched jet " << " id : " << pfj.id << " " << clj.id 
	  //      << " pfpt :  " << pfj.pt << " calopt : " << clj.pt 
	  //      << " delr :  " << delr   << " delpt  : " << delpt <<endl;
	  
	  pfid  [pfj.id] = 1;
	  caloid[clj.id] = 1;
	
	}
      }

      // //! Unmatched jets
      for(itr = mCaloPFMatchedJets.begin(); itr != mCaloPFMatchedJets.end(); ++itr){
	CaloPFJetPair jetpair = (*itr);
	Jet pfj = jetpair.first;
	Jet clj = jetpair.second;
	//float delr  = deltaR(pfj.eta, pfj.phi, clj.eta, clj.phi);
	//float delpt = fabs(pfj.pt - clj.pt);
	//if( pfid[pfj.id]==1 || caloid[clj.id]==1 )continue;
      
	if( pfid[pfj.id] == 0 ){//! Unmatched PF jet

	  pfpt      = jtpt      [pfj.id];
	  pfrawpt   = rawpt     [pfj.id];
	  pfeta     = jteta     [pfj.id];
	  pfphi     = jtphi     [pfj.id];
	  pfpu      = jtpu      [pfj.id];
	
	  chMax  = chargedMax[pfj.id];
	  phMax  = photonMax [pfj.id];
	  neMax  = neutralMax[pfj.id];
	  muMax  = muonMax   [pfj.id];
	  eMax   = eleMax    [pfj.id];
	
	  chSum  = chargedSum[pfj.id];
	  phSum  = photonSum [pfj.id];
	  neSum  = neutralSum[pfj.id];
	  muSum  = muonSum   [pfj.id];
	  eSum   = eleSum    [pfj.id];
	
	  hcalSum = hcalSum_pf[pfj.id];
	  ecalSum = ecalSum_pf[pfj.id];
	
	  unmatchedPFJets++;	  	  
	  unmatchPFJets->Fill();
	  pfid  [pfj.id] = 1;	  
	  // cout <<"\t %%%  UnMatched PF jet " << " id     : " << pfj.id << " " << clj.id 
	  //      << " pfpt :  " << pfj.pt   << " calopt : " << clj.pt 
	  //      << " delr : "  << delr     << " delpt  : " << delpt <<endl;
	}
      
	if( caloid[clj.id] == 0 ){//! Unmatched Calo jet

	  calopt    = jtpt_calo [clj.id];
	  calorawpt = rawpt_calo[clj.id];
	  caloeta   = jteta_calo[clj.id];
	  calophi   = jtphi_calo[clj.id];
	  calopu    = jtpu_calo [clj.id];
	
	  hcalSum = hcalSum_pf[clj.id];
	  ecalSum = ecalSum_pf[clj.id];
	
	  unmatchedCaloJets++;	  	  
	  unmatchCaloJets->Fill();
	  caloid[clj.id] = 1;	  
	  // cout <<"\t XXX  UnMatched Calo jet " << " id     : " << pfj.id << " " << clj.id 
	  //      << " pfpt :  " << pfj.pt   << " calopt : " << clj.pt 
	  //      << " delr : "  << delr     << " delpt  : " << delpt <<endl;
	}
      }
    }//! bothPFCalo
    
    // if( bothPFCalo    )std::cout<<" ****** Both PFCalo Event # " <<i;/*<<" vz  : "<<vz<<" hiBin : "<<hiBin;*/
    // else if( onlyCalo )std::cout<<" ****** Only Calo   Event # " <<i;/*<<" vz  : "<<vz<<" hiBin : "<<hiBin;*/
    // else if( onlyPF   )std::cout<<" ****** Only PF     Event # " <<i;/*<<" vz  : "<<vz<<" hiBin : "<<hiBin;*/
    // cout << " " 
    //   //<<" PF "<<nref<<" CALO  "<<nref_calo
    // 	 <<" PF : " << pfjets << " Calo : "  << calojets 
    // 	 <<" mjets : "<<  matchedJets 
    // 	 <<" umPF  : "<<  unmatchedPFJets 
    // 	 <<" umCalo: "<<  unmatchedCaloJets 
    // 	 <<std::endl;
    // cout << endl;

    //std::cout<<"Completed event #  "<<ievt<<std::endl; 
  }//! event loop ends

  //! Write to output file
  fout->cd();
  fout->Write();
  fout->Close();

  
  //! Check
  timer.Stop();
  double rtime  = timer.RealTime();
  double ctime  = timer.CpuTime();
  
  std::cout<<std::endl;
  std::cout<<Form("RealTime=%f seconds, CpuTime=%f seconds",rtime,ctime)<<std::endl;
  std::cout<<std::endl;
  std::cout<<"Good bye : " <<"\t"<<std::endl;
  return 1;
}

int GetEtaBin(float eta)
{
  for(int ix=0;ix<neta;ix++){
    if(eta>=etabins[ix] && eta<etabins[ix+1]){
      return ix;
    }
  }
  return -1;
}
int GetPhiBin(float phi)
{
  for(int ix=0;ix<nphi;ix++){
    if(phi>=phibins[ix] && phi<phibins[ix+1]){
      return ix;
    }
  }
  return -1;
}

int GetPtBin(float pt)
{
  for(int ix=0;ix<nbins;ix++){
    if(pt>=ptbins[ix] && pt<ptbins[ix+1]){
      return ix;
    }
  }
  return -1;
}
float delphi(float phi1, float phi2)
{
  float dphi = phi2 - phi1;
  dphi = fabs(atan2(sin(dphi),cos(dphi)));
  return dphi;
}
int GetCentBin(int bin)
{
  int ibin=-1;
  //! centrality is defined as 0.5% bins of cross section
  //! in 0-200 bins               
  if(bin<=10)ibin=0; //! 0-5%
  else if(bin>10  && bin<=20 )ibin=1; //! 5-10%
  else if(bin>20  && bin<=60 )ibin=2;  //! 10-30%
  else if(bin>60  && bin<=100)ibin=3;  //! 30-50%
  else if(bin>100 && bin<=140)ibin=4;  //! 50-70%
  else if(bin>140 && bin<=180)ibin=5;  //! 70-90%
  else if(bin>180 && bin<=200)ibin=6;  //! 90-100%
  return ibin;
}
float deltaR(float eta1, float phi1, float eta2, float phi2)
{
  float deta = eta1 - eta2;
  float dphi = fabs(phi1 - phi2);
  if(dphi > pi)dphi -= 2*pi;
  float dr = sqrt(pow(deta,2) + pow(dphi,2));
  return dr;
}
//void AddInputFiles(TChain *tch, string inputname, vector<string> inputTrees)
void AddInputFiles(TChain *tch, string inputname, string inputTree)
{
  //cout << " Add input files " << tch->GetName() << "  " << inputname.c_str() <<endl;
  //cout<<endl;
  std::stringstream tmpstr(inputname.c_str());
  std::string segment;
  std::vector<string> infiles;

  while(std::getline(tmpstr, segment, ',')){
    infiles.push_back(segment);
  }

  int im=0;
  for(std::vector<string>::iterator it=infiles.begin(); it!=infiles.end(); ++it, im++){
    std::string sFile = (*it);
    //   cout <<im<<"  "<< sFile.c_str() << endl;
    string stree = sFile+"/"+inputTree;
    tch->Add(stree.c_str());
  }
  cout<<endl;
}
