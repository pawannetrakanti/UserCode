
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

double GetXsec(double /*maxpthat*/);
void GetCentWeight(TH1F */*hCentWeight*/);


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
  //! Calo-PF match
  bool operator()(const CaloPFJetPair &A1, const CaloPFJetPair &A2){
    
    Jet cj1 = A1.first;  //! CaloJet 1st pair
    Jet pf1 = A1.second; //! PFJet   1st pair
    
    Jet cj2 = A2.first;  //! CaloJet 2nd pair
    Jet pf2 = A2.second; //! PFJet   2nd pair

    float delr1 = deltaR(pf1.eta, pf1.phi, cj1.eta, cj1.phi);
    float delr2 = deltaR(pf2.eta, pf2.phi, cj2.eta, cj2.phi);

    //float delpt1 = fabs(pf1.pt - cj1.pt);
    //float delpt2 = fabs(pf2.pt - cj2.pt);
    //return ((delpt1 < delpt2) && (delr1 < delr2) && (pf1.pt > pf2.pt) && cj1.pt > cj2.pt);
    
    return ((delr1 < delr2) && (cj1.pt > cj2.pt));
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
#define pi 3.14159265

const double ketacut=2.0;
const double kptrawcut =0.0;
const double kptrecocut=30.0;
const double kdelrmatch=0.2;
const double kdelrcut=0.3;

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

const double xsec[12][3] ={{2.034e-01,15  ,30}, //! 15                                                                    
                           {1.075e-02,30  ,50}, //! 30                                                                    
                           {1.025e-03,50  ,80}, //! 50                                                                    
                           {9.865e-05,80  ,120}, //! 80                                                                   
                           {1.129e-05,120 ,170}, //! 120                                                                  
                           {1.465e-06,170 ,220}, //! 170                                                                  
                           {2.837e-07,220 ,280}, //! 220                                                                  
                           {5.323e-08,280 ,370}, //! 280                                                                  
                           {5.934e-09,370 ,9999}, //! 370                                                                 
                           {0.0000000,9999,9999}
};


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
		std::string kDataset="data",
		//std::string fileList="/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat80_Track9_Jet30_matchEqR_merged_forest_0.root",
		std::string fileList="/mnt/hadoop/cms/store/user/belt/HiForest_jet55or65or80_JetRAA_v1_final/HiForest_jet55or65or80_JetRAA_v1_lumi9_Part8.root,/mnt/hadoop/cms/store/user/belt/HiForest_jet55or65or80_JetRAA_v1_final/HiForest_jet55or65or80_JetRAA_v1_lumi9_Part80.root,/mnt/hadoop/cms/store/user/belt/HiForest_jet55or65or80_JetRAA_v1_final/HiForest_jet55or65or80_JetRAA_v1_lumi9_Part81.root",
		//std::string foname="outputhisto_mc.root", 
		std::string foname="outputhisto_data.root", 
		double maxpthat=120
		)
{
  
  timer.Start();

  bool printDebug=false;


  //tr_jet = (TTree*)fin->Get("akPu3PFJetAnalyzer/t");
  TChain *tch_pfjet = new TChain(Form("%sPFJetAnalyzer/t",kAlgName.c_str()));
  AddInputFiles(tch_pfjet,fileList,Form("%sPFJetAnalyzer/t",kAlgName.c_str()));
  cout <<" # of events in PFJet Tree : " <<  tch_pfjet->GetEntries() <<endl;

  //tr_jet = (TTree*)fin->Get("akPu3CaloJetAnalyzer/t");
  //TChain *tch_calojet = new TChain(Form("%sCaloJetAnalyzer/t",kAlgName.c_str()));
  //AddInputFiles(tch_calojet,fileList,Form("%sCaloJetAnalyzer/t",kAlgName.c_str()));
  TChain *tch_calojet = new TChain("akPu3CaloJetAnalyzer/t");
  AddInputFiles(tch_calojet,fileList,"akPu3CaloJetAnalyzer/t");
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
  
  int sid_calo[1000];
  float refpt_calo[1000];
  float refeta_calo[1000];
  float refphi_calo[1000];
  float refdrjt_calo[1000];

  tch_calojet->SetBranchAddress("nref",&nref_calo);
  tch_calojet->SetBranchAddress("rawpt",rawpt_calo);
  tch_calojet->SetBranchAddress("jtpt" ,jtpt_calo);
  tch_calojet->SetBranchAddress("jtpu" ,jtpu_calo);
  tch_calojet->SetBranchAddress("jteta",jteta_calo);
  tch_calojet->SetBranchAddress("jtphi",jtphi_calo);
  tch_calojet->SetBranchAddress("hcalSum",hcalSum_calo);
  tch_calojet->SetBranchAddress("ecalSum",ecalSum_calo);
  if( kDataset == "mc" ){
    tch_calojet->SetBranchAddress("subid" ,sid_calo);
    tch_calojet->SetBranchAddress("refpt" ,refpt_calo);
    tch_calojet->SetBranchAddress("refeta",refeta_calo);
    tch_calojet->SetBranchAddress("refphi",refphi_calo);
    tch_calojet->SetBranchAddress("refdrjt",refdrjt_calo);
  }

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

  float pthat;
  int   sid[1000];
  float pfrefpt[1000];
  float pfrefeta[1000];
  float pfrefphi[1000];
  float pfrefdrjt[1000];

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
  if( kDataset == "mc" ){
    tch_pfjet->SetBranchAddress("pthat",&pthat);    
    tch_pfjet->SetBranchAddress("subid" ,sid);
    tch_pfjet->SetBranchAddress("refpt" ,pfrefpt);
    tch_pfjet->SetBranchAddress("refeta",pfrefeta);
    tch_pfjet->SetBranchAddress("refphi",pfrefphi);
    tch_pfjet->SetBranchAddress("refdrjt",pfrefdrjt);
  }
  
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

  if( kDataset == "mc"){
    tch_pfjet->SetBranchStatus("pthat",1,0);    
    tch_pfjet->SetBranchStatus("subid" ,1,0);
    tch_pfjet->SetBranchStatus("refpt" ,1,0);
    tch_pfjet->SetBranchStatus("refeta",1,0);
    tch_pfjet->SetBranchStatus("refphi",1,0);
    tch_pfjet->SetBranchStatus("refdrjt",1,0);
  }


  std::string outdir="";
  std::string outname=outdir+foname;
  TFile *fout = new TFile(outname.c_str(),"RECREATE");


  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"**************************************************** "<<std::endl;
  std::cout<<Form("%s %s %s ",kDataset.c_str(), ksp, kAlgName.c_str())<<std::endl;
  std::cout<<Form("Outfile : %s",outname.c_str())<<std::endl;
  std::cout<<Form("pT cut : %0.3f ; eta cut :  %0.3f ",kptrecocut, ketacut)<<std::endl;
  std::cout<<"**************************************************** "<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;

  //! Vertex re-weighting 
  TF1 *fVz=0;
  fVz = new TF1("fVz","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
  fVz->SetParameters(9.86748e-01, -8.91367e-03, 5.35416e-04, 2.67665e-06, -2.01867e-06);

  TH1F *hCentWeight = new TH1F("hCentWeight","Centrality weight",200,0,200);
  GetCentWeight(hCentWeight);


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

  Int_t subid;
  Float_t weight;
  Float_t refpt, refeta, refphi, refdrjt; 

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
  if( kDataset == "mc"){
    matchJets->Branch("subid",&subid,"subid/I");    
    matchJets->Branch("weight",&weight,"weight/F");    matchJets->Branch("pthat",&pthat,"pthat/F");       
    matchJets->Branch("refpt",&refpt,"refpt/F");       
    matchJets->Branch("refeta",&refeta,"refeta/F");       
    matchJets->Branch("refphi",&refphi,"refphi/F");       
    matchJets->Branch("refdrjt",&refdrjt,"refdrjt/F");       
  }

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
  if( kDataset == "mc"){
    unmatchPFJets->Branch("subid",&subid,"subid/I");    
    unmatchPFJets->Branch("weight",&weight,"weight/F");    unmatchPFJets->Branch("pthat",&pthat,"pthat/F");       
    unmatchPFJets->Branch("refpt",&refpt,"refpt/F");       
    unmatchPFJets->Branch("refeta",&refeta,"refeta/F");       
    unmatchPFJets->Branch("refphi",&refphi,"refphi/F");       
    unmatchPFJets->Branch("refdrjt",&refdrjt,"refdrjt/F");       
  }

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
  if( kDataset == "mc" ){
    unmatchCaloJets->Branch("subid",&subid,"subid/I");    
    unmatchCaloJets->Branch("weight",&weight,"weight/F");    unmatchCaloJets->Branch("pthat",&pthat,"pthat/F");
    unmatchCaloJets->Branch("refpt",&refpt,"refpt/F");    
    unmatchCaloJets->Branch("refeta",&refeta,"refeta/F"); 
    unmatchCaloJets->Branch("refphi",&refphi,"refphi/F"); 
    unmatchCaloJets->Branch("refdrjt",&refdrjt,"refdrjt/F"); 
  }

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

  weight=1.;
  double wxs=1.;
  if( kDataset == "mc" ){
    TEventList* el = new TEventList("el","el");
    stringstream selection; selection<<"pthat<="<<maxpthat;
    tch_pfjet->Draw(">>el",selection.str().c_str());
    double fentries = (double)el->GetN();
    std::cout<<"tree entries  :  "<<kAlgName.c_str()<<" algorithm : " << nentries<<" elist: "<< fentries <<std::endl;
    delete el;
    double tmpw = GetXsec(maxpthat);
    wxs = tmpw/(fentries/100000.);
  }

  //! Start event loop
  for (Long64_t i=0; i<nentries;i++) {
    nbytes += tch_pfjet->GetEntry(i);

    if(printDebug && i==20)break;

    float rndm=gRandom->Rndm();

    hEvents_Total->Fill(rndm);

    if(pcollisionEventSelection)hEvents_pCollEvent->Fill(rndm);
    if(pcollisionEventSelection && pHBHENoiseFilter)hEvents_pHBHENoise->Fill(rndm);
    if(pcollisionEventSelection && pHBHENoiseFilter && fabs(vz)<15. )hEvents_Vzcut->Fill(rndm);
    if(pcollisionEventSelection==0 || pHBHENoiseFilter==0 || fabs(vz) > 15.)continue;

    //if(!jet55_1 || !jet65_1 || !jet80_1)continue;
    
    int iCent = GetCentBin(hiBin);
    if(iCent<0 || iCent>=ncen)continue;

    //! SuperNovae events use calo jets
    //int lowjetCounter=0;
    int jetCounter=0;
    for(int g = 0;g<nref_calo;g++){
      // if( fabs(jteta_calo[g]) < 2. && jtpt_calo[g]>=kptrecocut ){ //to select inside
      // 	lowjetCounter++;
      // }
      if( fabs(jteta_calo[g]) < 2. && jtpt_calo[g]>=50. ){ //to select inside
	jetCounter++;
      }//eta selection cut
    }// jet loop

    // apply the correct supernova selection cut rejection here:
    if( hiNpix > 38000 - 500*jetCounter )continue;
    if( nref==0 && nref_calo==0 )continue;
    //if( lowjetCounter == 0 )continue;

    if( kDataset == "mc" && pthat > maxpthat )continue;

    hEvents_Cent->Fill(rndm);

    if(printDebug)std::cout << "------------------------------------Start Event # " << i <<"------------------------------------------------------------------ " << std::endl;
    //std::cout<<" ***** Event # " <<i<<"\t vz  : "<<vz<<"\t hiBin : "<<hiBin<<"\t  # of PF jets "<<nref<<" # of calojets  "<<nref_calo<<std::endl;
    if(i%50000==0)std::cout<<" ********** Event # " <<i<<"\t vz  : "<<vz<<"\t hiBin : "<<hiBin<<"\t  # of PF jets "<<nref<< "  # of Calo jets " <<nref_calo<<std::endl;

    int pfjets=0;
    int calojets=0;

    std::vector <Jet> vPFJets, vCaloJets;
    std::vector <int> pfid(nref), caloid(nref_calo);

    for(int pj=0; pj<nref; pj++){ //! PFjet

      if( rawpt[pj] < kptrawcut || jtpt[pj] < kptrecocut ) continue;
      if( fabs(jteta[pj]) > ketacut ) continue;

      if ( kDataset == "mc" ){
	if ( pfrefpt[pj] > 3.*pthat )continue;
      }
      //if( (eleMax[pj]/jtpt[pj])>=0.6 || (chargedMax[pj]/jtpt[pj])<=0.02 )continue;
      
      Jet pfj;
      pfj.id  = pj;
      pfj.eta = jteta[pj];
      pfj.phi = jtphi[pj];
      pfj.pt  = jtpt [pj];

      if(printDebug){
	std::cout << " PF jets : " << std::endl;
	if( kDataset == "mc" )std::cout <<"\t" << pj << "  pt : " << jtpt[pj] << " eta : "  << jteta[pj] << " subid : " << sid[pj] << std::endl;
	else std::cout <<"\t"<< pj << "  pt : " << jtpt[pj] << " eta : "  << jteta[pj] << std::endl;
      }
      vPFJets.push_back(pfj);
      pfjets++;
    }    
    if(printDebug)std::cout << std::endl;
    for(int cj=0; cj<nref_calo; cj++){ //! CaloJet
      
      if( rawpt_calo[cj] < kptrawcut || jtpt_calo[cj] < kptrecocut) continue;
      if( fabs(jteta_calo[cj]) > ketacut ) continue;

      if( kDataset == "mc" ){
	if ( refpt_calo[cj] > 3.*pthat )continue;
      }
      
      Jet clj;
      clj.id  = cj;
      clj.eta = jteta_calo[cj];
      clj.phi = jtphi_calo[cj];
      clj.pt  = jtpt_calo[cj];

      if(printDebug){
	std::cout << " Calo jets : " << std::endl;
	if( kDataset == "mc" )std::cout <<"\t" << cj << "  pt : " << jtpt_calo[cj] << " eta : "  << jteta_calo[cj] << " subid : " << sid_calo[cj] << std::endl;
	else  std::cout <<"\t" << cj << "  pt : " << jtpt_calo[cj] << " eta : "  << jteta_calo[cj] << std::endl;
      }
      vCaloJets.push_back(clj);
      calojets++;
    }//! calo jet loop

    if(printDebug)std::cout << std::endl;
    if(pfjets==0 && calojets==0){
      if(printDebug){
	std::cout <<" XXXXXXXXXXX  No Calo and PF jets passed the cuts " << std::endl;
	std::cout << "------------------------------------End Event # " << i <<"------------------------------------------------------------------ " << "\n"<<std::endl;
      }
      continue;
    }

    bool onlyCalo   = (pfjets==0 && calojets >0) ? true : false;
    bool onlyPF     = (pfjets>0  && calojets==0) ? true : false;
    bool bothPFCalo = (pfjets>0  && calojets >0) ? true : false;

    int matchedJets=0;
    int unmatchedPFJets=0;
    int unmatchedCaloJets=0;

    double tmpwt = 1.;
    if( kDataset == "mc" ){
      double wvz  = fVz->Eval(vz);
      double wcen = hCentWeight->GetBinContent(hCentWeight->FindBin(hiBin));
      tmpwt = (wxs*wvz*wcen);
    }else tmpwt = 1.;

    if(printDebug){
      std::cout <<" Total # of PF jets : " << pfjets << " Total # of calojets : "  << calojets <<"\n"<<std::endl;
      std::cout << std::endl;    
    }
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

	if( kDataset == "mc" ){
	  weight = tmpwt;
	  subid  = sid[pj];
	  refpt  = pfrefpt[pj];
	  refeta = pfrefeta[pj];
	  refphi = pfrefphi[pj];
	  refdrjt=pfrefdrjt[pj];
	}
	if(printDebug)std::cout <<" unmatched pf jets w ncalo=0 : " << unmatchedPFJets << "  pt : " << jtpt[pj] << " eta : "  << jteta[pj] << std::endl;
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

	if( kDataset == "mc" ){
	  weight = tmpwt;
	  subid  = sid_calo[cj];
	  refpt  = refpt_calo[cj];
	  refeta = refeta_calo[cj];
	  refphi = refphi_calo[cj];
	  refdrjt=refdrjt_calo[cj];
	}
	if(printDebug)std::cout <<" unmatched calo jets w npf=0 : " << unmatchedCaloJets << "  pt : " << jtpt_calo[cj] << " eta : "  << jteta_calo[cj] << std::endl;	
    	unmatchedCaloJets++;
    	unmatchCaloJets->Fill();
      }
    }
    else if( bothPFCalo ){

      CaloPFMatchedJets mCaloPFMatchedJets;
      for(iJet = vCaloJets.begin(); iJet != vCaloJets.end(); ++iJet){ //! Calojet      
	
	for(jJet = vPFJets.begin(); jJet != vPFJets.end(); ++jJet){ //! PFjet
	  
	  mCaloPFMatchedJets.insert(std::make_pair(*iJet,*jJet));
	  
	}//! calo jet loop
      }//! PF jet loop
      
      CPFItr itr;
      //! Matched jets (Calo jet matched to PF jet)
      for(itr = mCaloPFMatchedJets.begin(); itr != mCaloPFMatchedJets.end(); ++itr){

	CaloPFJetPair jetpair = (*itr);
	Jet clj = jetpair.first;
	Jet pfj = jetpair.second;

	float delr  = deltaR(pfj.eta, pfj.phi, clj.eta, clj.phi);
	//float delpt = fabs(pfj.pt - clj.pt);
	if( delr < kdelrmatch && caloid[clj.id]==0 ){
	
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
	  

	  if( kDataset == "mc" ){
	    weight = tmpwt;
	    subid  = sid [pfj.id];
	    refpt  = pfrefpt [pfj.id];
	    refeta = pfrefeta[pfj.id];
	    refphi = pfrefphi[pfj.id];
	    refdrjt=pfrefdrjt[pfj.id];
	  }

	  if(printDebug)std::cout <<"\t *** Matched jet " << " id : " << pfj.id << " " << clj.id << " pfpt :  " << pfj.pt << " calopt : " << clj.pt <<std::endl;
	  
	  matchedJets++;	  
	  matchJets->Fill();

	  pfid  [pfj.id] = 1;
	  caloid[clj.id] = 1;
	}
      }

      // //! Unmatched jets
      for(itr = mCaloPFMatchedJets.begin(); itr != mCaloPFMatchedJets.end(); ++itr){
	CaloPFJetPair jetpair = (*itr);

	Jet clj = jetpair.first;
	Jet pfj = jetpair.second;

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
	
	  if( kDataset == "mc" ){
	    weight = tmpwt;
	    subid  = sid[pfj.id];
	    refpt  = pfrefpt[pfj.id];
	    refeta = pfrefeta[pfj.id];
	    refphi = pfrefphi[pfj.id];
	    refdrjt= pfrefdrjt[pfj.id]; 
	  }

	  if(printDebug)std::cout <<"\t %%%  UnMatched PF   jet " << " id     : " << pfj.id << " pfpt :  " << pfj.pt << std::endl;
	  
	  unmatchedPFJets++;	  	  
	  unmatchPFJets->Fill();
	  pfid  [pfj.id] = 1;	  
	}
      
	if( caloid[clj.id] == 0 ){//! Unmatched Calo jet

	  calopt    = jtpt_calo [clj.id];
	  calorawpt = rawpt_calo[clj.id];
	  caloeta   = jteta_calo[clj.id];
	  calophi   = jtphi_calo[clj.id];
	  calopu    = jtpu_calo [clj.id];
	
	  hcalSum = hcalSum_pf[clj.id];
	  ecalSum = ecalSum_pf[clj.id];

	  if( kDataset == "mc" ){
	    weight = tmpwt;
	    subid  = sid_calo [clj.id];
	    refpt  = refpt_calo [clj.id];
	    refeta = refeta_calo[clj.id];
	    refphi = refphi_calo[clj.id];
	    refdrjt= refdrjt_calo[clj.id];
	  }

	  if(printDebug)std::cout <<"\t XXX  UnMatched Calo jet " << " id     : " << clj.id  << " calopt : " << clj.pt << std::endl;
	
	  unmatchedCaloJets++;	  	  
	  unmatchCaloJets->Fill();
	  caloid[clj.id] = 1;	  
	}
      }
    }//! bothPFCalo
    if(printDebug){
      std::cout << std::endl;
      if( bothPFCalo    )std::cout<<" ****** Both PFCalo Event # " <<i;/*<<" vz  : "<<vz<<" hiBin : "<<hiBin;*/
      else if( onlyCalo )std::cout<<" ****** Only Calo   Event # " <<i;/*<<" vz  : "<<vz<<" hiBin : "<<hiBin;*/
      else if( onlyPF   )std::cout<<" ****** Only PF     Event # " <<i;/*<<" vz  : "<<vz<<" hiBin : "<<hiBin;*/
      std::cout << " " 
	   <<" All PF  " << nref <<" CALO  "<< nref_calo << ";"  
	   <<" After cut PF : " << pfjets << " Calo : "  << calojets << ";"
	   <<" mjets : "<<  matchedJets 
	   <<" umPF  : "<<  unmatchedPFJets 
	   <<" umCalo: "<<  unmatchedCaloJets 
	   <<std::endl;
      std::cout << "------------------------------------End Event # " << i <<"------------------------------------------------------------------ " << "\n"<<std::endl;
    }
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
double GetXsec(double maxpt)
{
  //std::cout << " GetXsec() ::: max pt hat : " << maxpt << std::endl;                                                                          
  double effxsec=0;
  for(int i=0; i<11; i++){
    if(fabs(maxpt - xsec[i][2]) < 1e-08){
      effxsec = xsec[i][0] - xsec[i+1][0];
      //effxsec = xsec[i][0];                                                                                                                   
      //std::cout <<"\t \t  effective xsec : " << effxsec << "\t"<<  xsec[i][0] << "\t pthat : "<< xsec[i][1] << std::endl;                     
      return effxsec;
    }
  }
  return  1;
}
void GetCentWeight(TH1F *hCentWeight)
{
  hCentWeight->SetBinContent(1,6.7634);
  hCentWeight->SetBinContent(2,8.9638);
  hCentWeight->SetBinContent(3,8.32666);
  hCentWeight->SetBinContent(4,8.85033);
  hCentWeight->SetBinContent(5,7.48557);
  hCentWeight->SetBinContent(6,7.07842);
  hCentWeight->SetBinContent(7,7.17439);
  hCentWeight->SetBinContent(8,7.39451);
  hCentWeight->SetBinContent(9,7.03825);
  hCentWeight->SetBinContent(10,8.04546);
  hCentWeight->SetBinContent(11,6.32941);
  hCentWeight->SetBinContent(12,4.84289);
  hCentWeight->SetBinContent(13,7.05322);
  hCentWeight->SetBinContent(14,5.6361);
  hCentWeight->SetBinContent(15,5.59001);
  hCentWeight->SetBinContent(16,4.90395);
  hCentWeight->SetBinContent(17,5.13768);
  hCentWeight->SetBinContent(18,4.98226);
  hCentWeight->SetBinContent(19,3.76756);
  hCentWeight->SetBinContent(20,4.44141);
  hCentWeight->SetBinContent(21,4.01054);
  hCentWeight->SetBinContent(22,3.29702);
  hCentWeight->SetBinContent(23,3.21606);
  hCentWeight->SetBinContent(24,3.60559);
  hCentWeight->SetBinContent(25,3.36325);
  hCentWeight->SetBinContent(26,2.6244);
  hCentWeight->SetBinContent(27,3.17479);
  hCentWeight->SetBinContent(28,2.6614);
  hCentWeight->SetBinContent(29,2.1703);
  hCentWeight->SetBinContent(30,2.5898);
  hCentWeight->SetBinContent(31,2.56079);
  hCentWeight->SetBinContent(32,2.4732);
  hCentWeight->SetBinContent(33,2.02533);
  hCentWeight->SetBinContent(34,2.03333);
  hCentWeight->SetBinContent(35,1.75553);
  hCentWeight->SetBinContent(36,1.61111);
  hCentWeight->SetBinContent(37,1.55114);
  hCentWeight->SetBinContent(38,1.63142);
  hCentWeight->SetBinContent(39,1.57826);
  hCentWeight->SetBinContent(40,1.45919);
  hCentWeight->SetBinContent(41,1.5094);
  hCentWeight->SetBinContent(42,1.32344);
  hCentWeight->SetBinContent(43,1.36376);
  hCentWeight->SetBinContent(44,1.00481);
  hCentWeight->SetBinContent(45,0.924943);
  hCentWeight->SetBinContent(46,1.04942);
  hCentWeight->SetBinContent(47,0.976863);
  hCentWeight->SetBinContent(48,0.864318);
  hCentWeight->SetBinContent(49,0.772694);
  hCentWeight->SetBinContent(50,0.864502);
  hCentWeight->SetBinContent(51,0.829401);
  hCentWeight->SetBinContent(52,0.75895);
  hCentWeight->SetBinContent(53,0.702901);
  hCentWeight->SetBinContent(54,0.545314);
  hCentWeight->SetBinContent(55,0.571537);
  hCentWeight->SetBinContent(56,0.447445);
  hCentWeight->SetBinContent(57,0.565426);
  hCentWeight->SetBinContent(58,0.411747);
  hCentWeight->SetBinContent(59,0.381474);
  hCentWeight->SetBinContent(60,0.389027);
  hCentWeight->SetBinContent(61,0.345071);
  hCentWeight->SetBinContent(62,0.39263);
  hCentWeight->SetBinContent(63,0.340061);
  hCentWeight->SetBinContent(64,0.363676);
  hCentWeight->SetBinContent(65,0.342351);
  hCentWeight->SetBinContent(66,0.311693);
  hCentWeight->SetBinContent(67,0.2503);
  hCentWeight->SetBinContent(68,0.258714);
  hCentWeight->SetBinContent(69,0.269137);
  hCentWeight->SetBinContent(70,0.278774);
  hCentWeight->SetBinContent(71,0.254098);
  hCentWeight->SetBinContent(72,0.198022);
  hCentWeight->SetBinContent(73,0.217276);
  hCentWeight->SetBinContent(74,0.233573);
  hCentWeight->SetBinContent(75,0.233931);
  hCentWeight->SetBinContent(76,0.205296);
  hCentWeight->SetBinContent(77,0.187256);
  hCentWeight->SetBinContent(78,0.206262);
  hCentWeight->SetBinContent(79,0.192841);
  hCentWeight->SetBinContent(80,0.174259);
  hCentWeight->SetBinContent(81,0.157487);
  hCentWeight->SetBinContent(82,0.1807);
  hCentWeight->SetBinContent(83,0.135957);
  hCentWeight->SetBinContent(84,0.143054);
  hCentWeight->SetBinContent(85,0.158412);
  hCentWeight->SetBinContent(86,0.158663);
  hCentWeight->SetBinContent(87,0.130637);
  hCentWeight->SetBinContent(88,0.105144);
  hCentWeight->SetBinContent(89,0.109533);
  hCentWeight->SetBinContent(90,0.115536);
  hCentWeight->SetBinContent(91,0.103691);
  hCentWeight->SetBinContent(92,0.0988995);
  hCentWeight->SetBinContent(93,0.0899957);
  hCentWeight->SetBinContent(94,0.091202);
  hCentWeight->SetBinContent(95,0.0947045);
  hCentWeight->SetBinContent(96,0.0990303);
  hCentWeight->SetBinContent(97,0.074485);
  hCentWeight->SetBinContent(98,0.0904833);
  hCentWeight->SetBinContent(99,0.0745771);
  hCentWeight->SetBinContent(100,0.0746246);
  hCentWeight->SetBinContent(101,0.0666776);
  hCentWeight->SetBinContent(102,0.0631808);
  hCentWeight->SetBinContent(103,0.0645528);
  hCentWeight->SetBinContent(104,0.0721828);
  hCentWeight->SetBinContent(105,0.0640522);
  hCentWeight->SetBinContent(106,0.0544978);
  hCentWeight->SetBinContent(107,0.0602298);
  hCentWeight->SetBinContent(108,0.052432);
  hCentWeight->SetBinContent(109,0.0499806);
  hCentWeight->SetBinContent(110,0.05452);
  hCentWeight->SetBinContent(111,0.0456856);
  hCentWeight->SetBinContent(112,0.0464227);
  hCentWeight->SetBinContent(113,0.0389109);
  hCentWeight->SetBinContent(114,0.0429926);
  hCentWeight->SetBinContent(115,0.0423068);
  hCentWeight->SetBinContent(116,0.0436439);
  hCentWeight->SetBinContent(117,0.032317);
  hCentWeight->SetBinContent(118,0.0351724);
  hCentWeight->SetBinContent(119,0.0378572);
  hCentWeight->SetBinContent(120,0.0356574);
  hCentWeight->SetBinContent(121,0.0300515);
  hCentWeight->SetBinContent(122,0.0294732);
  hCentWeight->SetBinContent(123,0.0279459);
  hCentWeight->SetBinContent(124,0.0275134);
  hCentWeight->SetBinContent(125,0.0274872);
  hCentWeight->SetBinContent(126,0.0262874);
  hCentWeight->SetBinContent(127,0.0228082);
  hCentWeight->SetBinContent(128,0.0268362);
  hCentWeight->SetBinContent(129,0.0235638);
  hCentWeight->SetBinContent(130,0.019708);
  hCentWeight->SetBinContent(131,0.0203582);
  hCentWeight->SetBinContent(132,0.0191097);
  hCentWeight->SetBinContent(133,0.0169256);
  hCentWeight->SetBinContent(134,0.018112);
  hCentWeight->SetBinContent(135,0.0175009);
  hCentWeight->SetBinContent(136,0.0144258);
  hCentWeight->SetBinContent(137,0.0155731);
  hCentWeight->SetBinContent(138,0.0135958);
  hCentWeight->SetBinContent(139,0.0129593);
  hCentWeight->SetBinContent(140,0.0134124);
  hCentWeight->SetBinContent(141,0.0102854);
  hCentWeight->SetBinContent(142,0.00902376);
  hCentWeight->SetBinContent(143,0.00938477);
  hCentWeight->SetBinContent(144,0.00979958);
  hCentWeight->SetBinContent(145,0.00981297);
  hCentWeight->SetBinContent(146,0.00830205);
  hCentWeight->SetBinContent(147,0.00828065);
  hCentWeight->SetBinContent(148,0.0075616);
  hCentWeight->SetBinContent(149,0.00721783);
  hCentWeight->SetBinContent(150,0.00742391);
  hCentWeight->SetBinContent(151,0.00668121);
  hCentWeight->SetBinContent(152,0.00490303);
  hCentWeight->SetBinContent(153,0.00689083);
  hCentWeight->SetBinContent(154,0.00620564);
  hCentWeight->SetBinContent(155,0.00501006);
  hCentWeight->SetBinContent(156,0.00467418);
  hCentWeight->SetBinContent(157,0.00358751);
  hCentWeight->SetBinContent(158,0.0043082);
  hCentWeight->SetBinContent(159,0.00353042);
  hCentWeight->SetBinContent(160,0.00356054);
  hCentWeight->SetBinContent(161,0.00277187);
  hCentWeight->SetBinContent(162,0.00259774);
  hCentWeight->SetBinContent(163,0.0026294);
  hCentWeight->SetBinContent(164,0.00266786);
  hCentWeight->SetBinContent(165,0.00251157);
  hCentWeight->SetBinContent(166,0.00218918);
  hCentWeight->SetBinContent(167,0.00229047);
  hCentWeight->SetBinContent(168,0.00178743);
  hCentWeight->SetBinContent(169,0.00182462);
  hCentWeight->SetBinContent(170,0.00204086);
  hCentWeight->SetBinContent(171,0.00189708);
  hCentWeight->SetBinContent(172,0.00203718);
  hCentWeight->SetBinContent(173,0.0020711);
  hCentWeight->SetBinContent(174,0.00180765);
  hCentWeight->SetBinContent(175,0.00159439);
  hCentWeight->SetBinContent(176,0.00216191);
  hCentWeight->SetBinContent(177,0.00136735);
  hCentWeight->SetBinContent(178,0.00182475);
  hCentWeight->SetBinContent(179,0.00160661);
  hCentWeight->SetBinContent(180,0.00138471);
  hCentWeight->SetBinContent(181,0.00156103);
  hCentWeight->SetBinContent(182,0.00200855);
  hCentWeight->SetBinContent(183,0.0023071);
  hCentWeight->SetBinContent(184,0.00211314);
  hCentWeight->SetBinContent(185,0.00155022);
  hCentWeight->SetBinContent(186,0.00204334);
  hCentWeight->SetBinContent(187,0.00180985);
  hCentWeight->SetBinContent(188,0.00165799);
  hCentWeight->SetBinContent(189,0.00253497);
  hCentWeight->SetBinContent(190,0.00271872);
  hCentWeight->SetBinContent(191,0.00223219);
  hCentWeight->SetBinContent(192,0.00272361);
  hCentWeight->SetBinContent(193,0.00296343);
  hCentWeight->SetBinContent(194,0.00455219);
  hCentWeight->SetBinContent(195,0.00947736);
  hCentWeight->SetBinContent(196,0.0159602);
  hCentWeight->SetBinContent(197,0.0463495);
  hCentWeight->SetBinContent(198,0.156464);
  hCentWeight->SetBinContent(199,0);
  hCentWeight->SetBinContent(200,0);
}


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
