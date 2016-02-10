
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
#include <TCut.h>

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


//! constants
#define iYear 2015

#define pi 3.14159265

#define  ketacut    2.0
#define  kptrawcut  0.0
#define  kptrecocut 0.0
#define  kdelrmatch 0.2
#define  kptmatch   20.0
#define  kdelrcut   0.3
#define  kvzcut     15.0

const int ntrig=3;
double pi2 = 2*pi -1.;

int rbins=50;
double rbinl=0.0;
double rbinh=2.0;


void AddInputFiles(TChain */*ch*/, string /*iname*/, string /*inputTree*/);

int GetPhiBin(float /*phi*/);
int GetEtaBin(float /*eta*/);
int GetEtaBinWide(float /*eta*/);
int GetPtBin(float /*pt*/);
int GetPtBinWide(float /*pt*/);
int GetCentBin(int /*hiBin*/);
float delPhi(float /*phi1*/, float /*phi2*/);
float deltaR(float /*eta1*/, float /*phi1*/, 
	      float /*eta2*/, float /*phi2*/);

TF1 *fSmear=new TF1("fSmear","pol4",0.0,800.0);
double GetSmearedPt(int /*icent*/, double /*recopt*/, double /*refpt*/);


double GetXsecWt(double /*pthat*/, string /*kspecies*/);
double GetXsec(double /*maxpthat*/);
void GetCentWeight(TH1D */*hCentWeight*/);


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

    float delr1 = deltaR(cj1.eta, cj1.phi, pf1.eta, pf1.phi);
    float delr2 = deltaR(cj2.eta, cj2.phi, pf2.eta, pf2.phi);
    
    //float delpt1 = fabs(cj1.pt - pf1.pt);
    //float delpt2 = fabs(cj2.pt - pf2.pt);
    
    return ( (delr1 < delr2) && (cj1.pt > cj2.pt) );
    //return ( (delr1 < delr2) && (cj1.pt > cj2.pt) && (delpt1 < delpt2) );
  }
};

typedef std::multiset< CaloPFJetPair, CompareMatchedJets > CaloPFMatchedJets;
typedef std::multiset< CaloPFJetPair >::iterator CPFItr;
//typedef std::multiset< CaloPFJetPair >::value_type CPFJet;

const int ncen=8;
const char *cdir [ncen] = {"05","510","1030","3050","5070","7090","90100","pp"};
const char *ccent[ncen] = {"0-5%","5-10%","10-30%","30-50%","50-70%","70-90%","90-100%","pp"};

const int ncand=5;
const char *ccand[ncand] = {"h^{#pm}","e^{#pm}","#mu^{#pm}","#gamma","h0"};


//! pt binning
// double ptbins[] ={24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 
// 		  114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 
// 		  395, 430, 468, 507, 548, 592, 638, 686, 1000};
double ptbins[]={40, 50 ,60 ,70 ,80 ,90 ,100, 110, 120, 130, 140, 160, 200, 250, 300, 400, 548};
const int nbins = sizeof(ptbins)/sizeof(double) - 1;

const double ptbins_wide[]={30,50,80,120,200,340,548};
const int nbins_wide = sizeof(ptbins_wide)/sizeof(double) - 1;

//const double etabins[] = {-2.000, -1.4000, -0.4500, 0.000, 0.4500, 1.400, 2.000};
const double etabins[] = {-2.000, -1.4000, -0.4500, 0.000, 0.4500, 1.400, 2.000};
const int neta = sizeof(etabins)/sizeof(double) - 1;

const double etabins_wide[] = {0.00, 1.00, 1.80, 2.00};
const int neta_wide = sizeof(etabins_wide)/sizeof(double) - 1;

const double phibins[] = {-3.141,-2.100,-1.500,-0.800,-0.300, 
 			  0.300,0.800,1.500,2.100,3.141
};
const int nphi = sizeof(phibins)/sizeof(double) - 1;

double xsec[12][3] ={{2.034e-01,15.  ,30.},   //! 15     0   
		     {1.075e-02,30.  ,50.},   //! 30     1   
		     {1.025e-03,50.  ,80.},   //! 50     2   
		     {9.865e-05,80.  ,120.},  //! 80     3  
		     {1.129e-05,120. ,170.},  //! 120    4  
		     {1.465e-06,170. ,220.},  //! 170    5  
		     {2.837e-07,220. ,280.},  //! 220    6  
		     {5.323e-08,280. ,370.},  //! 280    7  
		     {5.934e-09,370. ,460.},  //! 370    8 
		     {8.125e-10,460. ,540.},  //! 460    9
		     {1.467e-10,540. ,9999.}, //! 540    10
		     {0.0000000,9999.,9999.}  //         11 
};

//! PbPb Weight factors
//!                         pthat effEv    effxec       wt
double pthatwt_pbpb[9][4]={{15,   315269,  0.19265   ,  6.11065e-07},
			   {30,   243155,  0.009725  ,  3.99951e-08},
			   {50,   379848,  0.00092635,  2.43874e-09},
			   {80,   362476,  8.736e-05 ,  2.41009e-10},
			   {120,  359581,  9.825e-06 ,  2.73235e-11},
			   {170,  359272,  1.1813e-06,  3.28804e-12},
			   {220,  217873,  2.3047e-07,  1.05782e-12},
			   {280,   88991,  4.7296e-08,  5.31469e-13},
			   {370,   25252,  5.934e-09 ,  2.34991e-13}
};
			   
//! pp Weight factors
double pthatwt_pp[11][4] ={{15 ,  4570974,  0.19265   ,   4.21464e-08},
			   {30 ,  4665419,  0.009725  ,   2.08449e-09},
			   {50 ,  4902822,  0.00092635,   1.88942e-10},
			   {80 ,  4871373,  8.736e-05 ,   1.79333e-11},
			   {120,  4840555,  9.825e-06 ,   2.02973e-12},
			   {170,  4603464,  1.1813e-06,   2.56611e-13},
			   {220,  4913932,  2.3047e-07,   4.69013e-14},
			   {280,  5384459,  4.7296e-08,   8.78380e-15},
			   {370,  4870894,  5.1215e-09,   1.05145e-15},
			   {460,  4655293,  6.658e-10 ,   1.43020e-16},
			   {540,  7171716,  1.467e-10 ,   2.04554e-17}
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

// int jetmatch(std::string kSpecies="pbpb",
// 	     std::string kAlgName="akPu3",
// 	     std::string kDataset="mc",
// 	     //! PbPb
// 	     std::string kFileList="/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/HiForest_PYTHIA_HYDJET_pthat80_Track9_Jet30_matchEqR_merged_forest_0.root",
// 	     std::string kFoname="output_tree_mc.root", 
// 	     double kMaxpthat=120.
// 	     )


int jetmatch(std::string kSpecies="pp",
 	     std::string kAlgName="ak3",
 	     std::string kDataset="mc",
 	     std::string kFileList="/mnt/hadoop/cms/store/user/velicanu/HiForest_pp_Offical_MC_pthat_80_53X_STARTHI53_V29_5_3_20_correctJEC_pawan_30Nov2014_hadd/0.root",
 	     std::string kFoname="output_tree_mc.root", 
 	     double kMaxpthat=120.
 	     )
{
  
  timer.Start();

  std::cout << std::endl;


  Float_t Jet55_prescl = 1;//2.0475;
  Float_t Jet40_prescl = 1;//9.275;
  int iSet=1;
  if( kDataset =="data"){
    iSet=0;
    if( kSpecies == "pbpb")Jet55_prescl = 2.0475;
    else if( kSpecies == "pp" )Jet40_prescl = 9.275;
  }
  
  bool printDebug=false;

  //! For PbPb MC
  if( (kSpecies == "pbpb" || kSpecies == "pbpb_mb") && kDataset == "mc" ){
    xsec[8][2]=9999; //! 370 in PbPb
    xsec[9][0]=0.00000;  xsec[9][1]=9999;  xsec[9][2]=9999;
  }

  //tr_jet = (TTree*)fin->Get("akPu3PFJetAnalyzer/t");
  TChain *tch_pfjet = new TChain(Form("%sPFJetAnalyzer/t",kAlgName.c_str()));
  AddInputFiles(tch_pfjet,kFileList,Form("%sPFJetAnalyzer/t",kAlgName.c_str()));
  cout <<" # of events in PFJet   Tree : " <<  tch_pfjet->GetEntries() <<endl;

  //tr_jet = (TTree*)fin->Get("akPu3CaloJetAnalyzer/t");
  //TChain *tch_calojet = new TChain(Form("%sCaloJetAnalyzer/t",kAlgName.c_str()));
  //AddInputFiles(tch_calojet,kFileList,Form("%sCaloJetAnalyzer/t",kAlgName.c_str()));
  TChain *tch_calojet=0;
  if( kSpecies == "pbpb" || kSpecies == "pbpb_mb" ){
    tch_calojet = new TChain("akPu3CaloJetAnalyzer/t");
    AddInputFiles(tch_calojet,kFileList,"akPu3CaloJetAnalyzer/t");
    cout <<" # of events in CaloJet Tree : " <<  tch_calojet->GetEntries() <<endl;
  }else{
    tch_calojet = new TChain("ak3CaloJetAnalyzer/t");
    AddInputFiles(tch_calojet,kFileList,"ak3CaloJetAnalyzer/t");
    cout <<" # of events in CaloJet Tree : " <<  tch_calojet->GetEntries() <<endl;
  }
  //tr_ev = (TTree*)fin->Get("hiEvtAnalyzer/HiTree");
  TChain *tch_ev = new TChain("hiEvtAnalyzer/HiTree");
  AddInputFiles(tch_ev,kFileList,"hiEvtAnalyzer/HiTree");
  cout <<" # of events in Event   Tree : " <<  tch_ev->GetEntries() <<endl;

  //tr_hlt = (TTree*)fin->Get("hltanalysis/HltTree");
  TChain *tch_hlt = new TChain("hltanalysis/HltTree");
  AddInputFiles(tch_hlt,kFileList,"hltanalysis/HltTree");
  cout <<" # of events in HLT     Tree : " <<  tch_hlt->GetEntries() <<endl;

  //tr_skim = (TTree*)fin->Get("skimanalysis/HltTree");
  TChain *tch_skim = new TChain("skimanalysis/HltTree");
  AddInputFiles(tch_skim,kFileList,"skimanalysis/HltTree");
  cout <<" # of events in Skim    Tree : " <<  tch_skim->GetEntries() <<endl;

  //tr_trobj = (TTree*)fin->Get("hltobject/jetObjTree");  
  TChain *tch_trgobj = new TChain("hltobject/jetObjTree");
  if( kDataset == "data" ){
    AddInputFiles(tch_trgobj,kFileList,"hltobject/jetObjTree");
    cout <<" # of events in TrigObj Tree : " <<  tch_trgobj->GetEntries() <<endl;
    cout <<endl;
  }

  //! Event Tree
  int run_value;
  int evt_value;
  int lumi_value;
  int hiNpix;
  int hiBin;
  float vz;

  float hiHF;
  float hiZDC;
  float hiZDCminus;
  tch_ev->SetBranchAddress("run",&run_value);  
  tch_ev->SetBranchAddress("evt",&evt_value);  
  tch_ev->SetBranchAddress("lumi",&lumi_value);  
  tch_ev->SetBranchAddress("hiBin",&hiBin);
  tch_ev->SetBranchAddress("hiNpix",&hiNpix);  
  tch_ev->SetBranchAddress("hiHF",&hiHF);  
  tch_ev->SetBranchAddress("hiZDC",&hiZDC);  
  tch_ev->SetBranchAddress("hiZDCminus",&hiZDCminus);  
  tch_ev->SetBranchAddress("vz",&vz);


  //! HLT tree
  //! HI
  int jet55;
  int jet55_prescl;
  int jet65;
  int jet65_prescl;
  int jet80;
  int jet80_prescl;
  //! PP
  int jet40;
  int jet40_prescl;
  int jet60;
  int jet60_prescl;

  //! MinBias
  int L1_MB;
  int L1_MB_p;
  int L1_sj36;
  int L1_sj52;
  int jetMB;
  int jetMB_prescl;
  int L1_sj36_p;
  int L1_sj52_p;

  if ( kSpecies == "pbpb" || kSpecies == "pbpb_mb" ){
    if( kDataset == "data" ){
      tch_hlt->SetBranchAddress("HLT_HIJet55_v1",&jet55);
      tch_hlt->SetBranchAddress("HLT_HIJet55_v1_Prescl",&jet55_prescl);
      tch_hlt->SetBranchAddress("HLT_HIJet65_v1",&jet65);
      tch_hlt->SetBranchAddress("HLT_HIJet65_v1_Prescl",&jet65_prescl);
      tch_hlt->SetBranchAddress("HLT_HIJet80_v1",&jet80);
      tch_hlt->SetBranchAddress("HLT_HIJet80_v1_Prescl",&jet80_prescl);
      tch_hlt->SetBranchAddress("HLT_HIMinBiasHfOrBSC_v1",&jetMB);
      tch_hlt->SetBranchAddress("HLT_HIMinBiasHfOrBSC_v1_Prescl",&jetMB_prescl);
    }else{
      tch_hlt->SetBranchAddress("HLT_HIJet55_v7",&jet55);
      tch_hlt->SetBranchAddress("HLT_HIJet55_v7_Prescl",&jet55_prescl);
      tch_hlt->SetBranchAddress("HLT_HIJet65_v7",&jet65);
      tch_hlt->SetBranchAddress("HLT_HIJet65_v7_Prescl",&jet65_prescl);
      tch_hlt->SetBranchAddress("HLT_HIJet80_v7",&jet80);
      tch_hlt->SetBranchAddress("HLT_HIJet80_v7_Prescl",&jet80_prescl);
      tch_hlt->SetBranchAddress("HLT_HIMinBiasHfOrBSC_v4",&jetMB);
      tch_hlt->SetBranchAddress("HLT_HIMinBiasHfOrBSC_v4_Prescl",&jetMB_prescl);
    }
    tch_hlt->SetBranchAddress("L1_ZeroBias",&L1_MB);
    tch_hlt->SetBranchAddress("L1_ZeroBias_Prescl",&L1_MB_p);
    tch_hlt->SetBranchAddress("L1_SingleJet36_BptxAND",&L1_sj36);
    tch_hlt->SetBranchAddress("L1_SingleJet36_BptxAND_Prescl",&L1_sj36_p);
    tch_hlt->SetBranchAddress("L1_SingleJet52_BptxAND",&L1_sj52);
    tch_hlt->SetBranchAddress("L1_SingleJet52_BptxAND_Prescl",&L1_sj52_p);
  }else {
    tch_hlt->SetBranchAddress("HLT_PAJet40_NoJetID_v1",&jet40);
    tch_hlt->SetBranchAddress("HLT_PAJet40_NoJetID_v1_Prescl",&jet40_prescl);
    tch_hlt->SetBranchAddress("HLT_PAJet60_NoJetID_v1",&jet60);
    tch_hlt->SetBranchAddress("HLT_PAJet60_NoJetID_v1_Prescl",&jet60_prescl);
    tch_hlt->SetBranchAddress("HLT_PAJet80_NoJetID_v1",&jet80);
    tch_hlt->SetBranchAddress("HLT_PAJet80_NoJetID_v1_Prescl",&jet80_prescl);
  }
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
  if( kSpecies == "pbpb" || kSpecies == "pbpb_mb" )tch_skim->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
  else tch_skim->SetBranchAddress("pPAcollisionEventSelectionPA",&pcollisionEventSelection);
  tch_skim->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);

  //! Trigger object tree
  float trgObj_id=0;
  float trgObj_pt=0;
  float trgObj_eta=0;
  float trgObj_phi=0;
  if( kDataset == "data" ){
    tch_trgobj->SetBranchAddress("id",&trgObj_id);
    tch_trgobj->SetBranchAddress("pt",&trgObj_pt);
    tch_trgobj->SetBranchAddress("eta",&trgObj_eta);
    tch_trgobj->SetBranchAddress("phi",&trgObj_phi);
  }

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
  int   nref_pf;
  float jtpt[1000];
  float rawpt[1000];
  float jteta[1000];
  float jtpu [1000];
  float jtphi[1000];
  float trackSum[1000];  
  float chargedSum[1000];
  float neutralSum[1000];
  float photonSum[1000];
  float eleSum[1000];
  float muonSum[1000];
  float trackMax[1000];  
  float chargedMax[1000];
  float neutralMax[1000];
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
  int   pfrefparton_flavor[1000];

  tch_pfjet->SetBranchAddress("nref",&nref_pf);
  tch_pfjet->SetBranchAddress("rawpt",rawpt);
  tch_pfjet->SetBranchAddress("jtpt" ,jtpt);
  tch_pfjet->SetBranchAddress("jtpu" ,jtpu);
  tch_pfjet->SetBranchAddress("jteta",jteta);
  tch_pfjet->SetBranchAddress("jtphi",jtphi);
  tch_pfjet->SetBranchAddress("trackSum",trackSum);
  tch_pfjet->SetBranchAddress("chargedSum",chargedSum);
  tch_pfjet->SetBranchAddress("neutralSum",neutralSum);
  tch_pfjet->SetBranchAddress("photonSum",photonSum);
  tch_pfjet->SetBranchAddress("eSum",eleSum);
  tch_pfjet->SetBranchAddress("muSum",muonSum);
  tch_pfjet->SetBranchAddress("trackMax",trackMax);
  tch_pfjet->SetBranchAddress("chargedMax",chargedMax);
  tch_pfjet->SetBranchAddress("neutralMax",neutralMax);
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
    tch_pfjet->SetBranchAddress("refparton_flavor",pfrefparton_flavor);
  }
  
  tch_pfjet->AddFriend(tch_ev);
  tch_pfjet->AddFriend(tch_hlt);
  tch_pfjet->AddFriend(tch_skim);
  if( kDataset == "data")tch_pfjet->AddFriend(tch_trgobj);
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
  tch_pfjet->SetBranchStatus("trackSum",1,0);
  tch_pfjet->SetBranchStatus("chargedSum",1,0);
  tch_pfjet->SetBranchStatus("neutralSum",1,0);
  tch_pfjet->SetBranchStatus("photonSum",1,0);
  tch_pfjet->SetBranchStatus("eSum",1,0);
  tch_pfjet->SetBranchStatus("muSum",1,0);
  tch_pfjet->SetBranchStatus("trackMax",1,0);
  tch_pfjet->SetBranchStatus("chargedMax",1,0);
  tch_pfjet->SetBranchStatus("neutralMax",1,0);
  tch_pfjet->SetBranchStatus("photonMax",1,0);
  tch_pfjet->SetBranchStatus("eMax",1,0);
  tch_pfjet->SetBranchStatus("muMax",1,0);
  tch_pfjet->SetBranchStatus("hcalSum",1,0);
  tch_pfjet->SetBranchStatus("ecalSum",1,0);

  tch_pfjet->SetBranchStatus("run",1,0);
  tch_pfjet->SetBranchStatus("evt",1,0);
  tch_pfjet->SetBranchStatus("lumi",1,0);
  tch_pfjet->SetBranchStatus("hiNpix",1,0);
  tch_pfjet->SetBranchStatus("hiHF",1,0);
  tch_pfjet->SetBranchStatus("hiZDC",1,0);
  tch_pfjet->SetBranchStatus("hiZDCminus",1,0);
  tch_pfjet->SetBranchStatus("hiBin",1,0);
  tch_pfjet->SetBranchStatus("vz",1,0);



  if( kSpecies == "pbpb" || kSpecies == "pbpb_mb" ){
    if ( kDataset == "data" ){
      tch_pfjet->SetBranchStatus("HLT_HIJet55_v1",1,0);
      tch_pfjet->SetBranchStatus("HLT_HIJet55_v1_Prescl",1,0);
      tch_pfjet->SetBranchStatus("HLT_HIJet65_v1",1,0);
      tch_pfjet->SetBranchStatus("HLT_HIJet65_v1_Prescl",1,0);
      tch_pfjet->SetBranchStatus("HLT_HIJet80_v1",1,0);
      tch_pfjet->SetBranchStatus("HLT_HIJet80_v1_Prescl",1,0);
      tch_pfjet->SetBranchStatus("HLT_HIMinBiasHfOrBSC_v1",1,0);
      tch_pfjet->SetBranchStatus("HLT_HIMinBiasHfOrBSC_v1_Prescl",1,0);

      tch_pfjet->SetBranchStatus("id",1,0);
      tch_pfjet->SetBranchStatus("pt",1,0);
      tch_pfjet->SetBranchStatus("eta",1,0);
      tch_pfjet->SetBranchStatus("phi",1,0);

    }else{
      tch_pfjet->SetBranchStatus("HLT_HIJet55_v7",1,0);
      tch_pfjet->SetBranchStatus("HLT_HIJet55_v7_Prescl",1,0);
      tch_pfjet->SetBranchStatus("HLT_HIJet65_v7",1,0);
      tch_pfjet->SetBranchStatus("HLT_HIJet65_v7_Prescl",1,0);
      tch_pfjet->SetBranchStatus("HLT_HIJet80_v7",1,0);
      tch_pfjet->SetBranchStatus("HLT_HIJet80_v7_Prescl",1,0);
      tch_pfjet->SetBranchStatus("HLT_HIMinBiasHfOrBSC_v4",1,0);
      tch_pfjet->SetBranchStatus("HLT_HIMinBiasHfOrBSC_v4_Prescl",1,0);
    }
    tch_pfjet->SetBranchStatus("pcollisionEventSelection",1,0);
    tch_pfjet->SetBranchStatus("L1_ZeroBias",1,0);
    tch_pfjet->SetBranchStatus("L1_ZeroBias_Prescl",1,0);
    tch_pfjet->SetBranchStatus("L1_SingleJet36_BptxAND",1,0);
    tch_pfjet->SetBranchStatus("L1_SingleJet36_BptxAND_Prescl",1,0);
    tch_pfjet->SetBranchStatus("L1_SingleJet52_BptxAND",1,0);
    tch_pfjet->SetBranchStatus("L1_SingleJet52_BptxAND_Prescl",1,0);
  }else{
    tch_pfjet->SetBranchStatus("HLT_PAJet40_NoJetID_v1",1,0);
    tch_pfjet->SetBranchStatus("HLT_PAJet40_NoJetID_v1_Prescl",1,0);
    tch_pfjet->SetBranchStatus("HLT_PAJet60_NoJetID_v1",1,0);
    tch_pfjet->SetBranchStatus("HLT_PAJet60_NoJetID_v1_Prescl",1,0);
    tch_pfjet->SetBranchStatus("HLT_PAJet80_NoJetID_v1",1,0);
    tch_pfjet->SetBranchStatus("HLT_PAJet80_NoJetID_v1_Prescl",1,0);
    tch_pfjet->SetBranchStatus("pPAcollisionEventSelectionPA",1,0);
  }
  tch_pfjet->SetBranchStatus("pHBHENoiseFilter",1,0);
  
  if( kDataset == "mc"){
    tch_pfjet->SetBranchStatus("pthat",1,0);    
    tch_pfjet->SetBranchStatus("subid" ,1,0);
    tch_pfjet->SetBranchStatus("refpt" ,1,0);
    tch_pfjet->SetBranchStatus("refeta",1,0);
    tch_pfjet->SetBranchStatus("refphi",1,0);
    tch_pfjet->SetBranchStatus("refdrjt",1,0);
    tch_pfjet->SetBranchStatus("refparton_flavor",1,0);
  }

  //! Vertex re-weighting 
  TF1 *fVz=0;
  fVz = new TF1("fVz","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
  if( kSpecies == "pbpb" || kSpecies=="pbpb_mb" )fVz->SetParameters(9.86748e-01, -8.91367e-03, 5.35416e-04, 2.67665e-06, -2.01867e-06);
  else fVz->SetParameters(8.41684e-01,-2.58609e-02,4.86550e-03,-3.10581e-04,2.07918e-05);

  //! Centrality re-weighting 
  TH1D *hCentWeight = new TH1D("hCentWeight","Centrality weight",200,0,200);
  GetCentWeight(hCentWeight);

  //! Create output file
  std::string outdir="";
  std::string outfile=outdir+kFoname;
  TFile *fout = new TFile(outfile.c_str(),"RECREATE");

  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"**************************************************** "<<std::endl;
  std::cout<<Form("Dataset : %s, Species : %s, Jet Algorithm :  %s ",kDataset.c_str(), kSpecies.c_str(), kAlgName.c_str())<<std::endl;
  std::cout<<Form("Outfile : %s",outfile.c_str())<<std::endl;
  std::cout<<Form("vertex z (c.m.) cut : %0.3f ",kvzcut)<<std::endl;
  std::cout<<Form("Reco pT cut : %0.3f ; Reco eta cut : %0.3f ",kptrecocut, ketacut)<<std::endl;
  std::cout<<Form("CALO-PF jet matching delta R cut   : %0.3f ",kdelrmatch)<<std::endl;
  std::cout<<"**************************************************** "<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;

  //! 
  //! Define histograms here
  fout->mkdir(Form("%sJetAnalyzer",kAlgName.c_str()));
  fout->cd(Form("%sJetAnalyzer"   ,kAlgName.c_str()));


  TH1D *hEvents_Total = new TH1D("hEvents_Total","Total # of events Total ",10,0.,1.);
  hEvents_Total->Sumw2();
  TH1D *hEvents_pCollEvent = new TH1D("hEvents_pCollEvent","# of events pCollisionEvent",10,0.,1.);
  hEvents_pCollEvent->Sumw2();
  TH1D *hEvents_pHBHENoise = new TH1D("hEvents_pHBHENoise","# of events HBHENoise",10,0.,1.);
  hEvents_pHBHENoise->Sumw2();
  TH1D *hEvents_Vzcut = new TH1D("hEvents_Vzcut","# of events vz cut ",10,0.,1.);
  hEvents_Vzcut->Sumw2();
  TH1D *hEvents_Good = new TH1D("hEvents_Good","# of events of good events",10,0.,1.);
  hEvents_Good->Sumw2();
  TH1D *hEvents_bad = new TH1D("hEvents_bad","# of events bad ",10,0.,1.);
  hEvents_bad->Sumw2();
  TH1D *hEvents_supernova = new TH1D("hEvents_supernova","supernova # of events ",10,0.,1.);
  hEvents_supernova->Sumw2();
  TH1D *hEvents_nopfcalo = new TH1D("hEvents_nopfcalo","nopfcalo # of events ",10,0.,1.);
  hEvents_nopfcalo->Sumw2();
  TH1D *hEvents_maxpthat = new TH1D("hEvents_maxpthat","maxpthat # of events ",10,0.,1.);
  hEvents_maxpthat->Sumw2();
  TH1D *hEvents_pileup = new TH1D("hEvents_pileup","pileup # of events ",10,0.,1.);
  hEvents_pileup->Sumw2();

  TH1D *hEvents_jet40 = new TH1D("hEvents_jet40","# of events jet40",10,0.,1.);
  hEvents_jet40->Sumw2();
  TH1D *hEvents_jet60 = new TH1D("hEvents_jet60","# of events jet60",10,0.,1.);
  hEvents_jet60->Sumw2();
  TH1D *hEvents_jet55 = new TH1D("hEvents_jet55","# of events jet55",10,0.,1.);
  hEvents_jet55->Sumw2();
  TH1D *hEvents_jet65 = new TH1D("hEvents_jet65","# of events jet65",10,0.,1.);
  hEvents_jet65->Sumw2();
  TH1D *hEvents_jet80 = new TH1D("hEvents_jet80","# of events jet80",10,0.,1.);
  hEvents_jet80->Sumw2();

  TH1D *hEvents_jet40_nojet60_nojet80 = new TH1D("hEvents_jet40_nojet60_nojet80","# of events jet40 && !jet60 && !jet80",10,0.,1.);
  hEvents_jet40_nojet60_nojet80->Sumw2();
  TH1D *hEvents_jet60_nojet80 = new TH1D("hEvents_jet60_nojet80","# of events jet60 && !jet80",10,0.,1.);
  hEvents_jet60_nojet80->Sumw2();

  TH2F *hEvents_jet40_prescl = new TH2F("hEvents_jet40_prescl","prescaled # of events jet40",50,-0.5,50-0.5,10,0.,1.);
  hEvents_jet40_prescl->Sumw2();
  TH2F *hEvents_jet60_prescl = new TH2F("hEvents_jet60_prescl","prescaled # of events jet60",50,-0.5,50-0.5,10,0.,1.);
  hEvents_jet60_prescl->Sumw2();				
  TH2F *hEvents_jet55_prescl = new TH2F("hEvents_jet55_prescl","prescaled # of events jet55",50,-0.5,50-0.5,10,0.,1.);
  hEvents_jet55_prescl->Sumw2();				
  TH2F *hEvents_jet65_prescl = new TH2F("hEvents_jet65_prescl","prescaled # of events jet65",50,-0.5,50-0.5,10,0.,1.);
  hEvents_jet65_prescl->Sumw2();				
  TH2F *hEvents_jet80_prescl = new TH2F("hEvents_jet80_prescl","prescaled # of events jet80",50,-0.5,50-0.5,10,0.,1.);
  hEvents_jet80_prescl->Sumw2();
  
  TH1D *hEvents_jet55_nojet65_nojet80 = new TH1D("hEvents_jet55_nojet65_nojet80","# of events jet55 && !jet65 && !jet80",10,0.,1.);
  hEvents_jet55_nojet65_nojet80->Sumw2();
  TH1D *hEvents_jet65_nojet80 = new TH1D("hEvents_jet65_nojet80","# of events jet60 && !jet80",10,0.,1.);
  hEvents_jet65_nojet80->Sumw2();

  TH1D *hVz = new TH1D("hVz","Vertex z (in c.m.)",60,-15.,15.);
  hVz->Sumw2();

  TH1D *hBin = new TH1D("hBin","hiBin",200,-0.5,199.5);
  hBin->Sumw2();

  TH1D *hpthat = new TH1D("hpthat","pthat",1200,0.0,12000.0);
  hpthat->Sumw2();

  //! 0 : w/o jet id and 1:w jet id
  //! Resposnse 
  TH1D *hrecopt_genm[2][ncen][nbins];
  TH1D *hrescrpt_genm [2][ncen][nbins], *hresrrpt_genm [2][ncen][nbins];
  TH1D *hrescrpt_sm_genm[2][ncen][nbins];
  TH1D *hrescrpt_genm_eta [2][ncen][nbins][neta], *hresrrpt_genm_eta [2][ncen][nbins][neta];
  TH1D *hrescrpt_genm_phi [2][ncen][nbins][nphi], *hresrrpt_genm_phi [2][ncen][nbins][nphi];
  TH1D *hrescrpt_wide_genm_eta [2][ncen][nbins_wide][neta], *hresrrpt_wide_genm_eta [2][ncen][nbins_wide][neta];
  TH1D *hrescrpt_wide_genm_phi [2][ncen][nbins_wide][nphi], *hresrrpt_wide_genm_phi [2][ncen][nbins_wide][nphi];

  //!  quark jet energy scale
  TH1D *hrescrpt_q_genm [2][ncen][nbins], *hresrrpt_q_genm [2][ncen][nbins];
  TH1D *hrescrpt_q_genm_eta [2][ncen][nbins][neta], *hresrrpt_q_genm_eta [2][ncen][nbins][neta];

  //!  gluon jet energy scale
  TH1D *hrescrpt_g_genm [2][ncen][nbins], *hresrrpt_g_genm [2][ncen][nbins];
  TH1D *hrescrpt_g_genm_eta [2][ncen][nbins][neta], *hresrrpt_g_genm_eta [2][ncen][nbins][neta];

  //! Flavour
  TH1D *hrescrpt_u_genm [2][ncen][nbins], *hresrrpt_u_genm [2][ncen][nbins];
  TH1D *hrescrpt_u_genm_eta [2][ncen][nbins][neta], *hresrrpt_u_genm_eta [2][ncen][nbins][neta];
  TH1D *hrescrpt_d_genm [2][ncen][nbins], *hresrrpt_d_genm [2][ncen][nbins];
  TH1D *hrescrpt_d_genm_eta [2][ncen][nbins][neta], *hresrrpt_d_genm_eta [2][ncen][nbins][neta];
  TH1D *hrescrpt_s_genm [2][ncen][nbins], *hresrrpt_s_genm [2][ncen][nbins];
  TH1D *hrescrpt_s_genm_eta [2][ncen][nbins][neta], *hresrrpt_s_genm_eta [2][ncen][nbins][neta];

  //! Background
  TH1D *havbkgpt_comb[2][ncen][nbins];
  TH1D *havbkgpt_comb_eta[2][ncen][nbins][neta_wide];

  TH1D *havbkgpt_ind[2][ncen][ntrig][nbins];
  TH1D *havbkgpt_ind_eta[2][ncen][ntrig][nbins][neta_wide];
  for(int k=0; k<2; k++){
    for(int ic=0; ic<ncen; ic++){
      for(int it=0;it<ntrig;it++){
	for(int ip=0; ip<nbins; ip++){
	  havbkgpt_ind[k][ic][it][ip]= new TH1D(Form("havbkgpt_ind_%sPF_%d_%d_%d_%d_%d",kAlgName.c_str(),iSet,k,ic,it,ip),
						//Form("(<bkg>) jet p_{T} %sPF %s %0.0f < p_{T}^{REF} < %0.0f",kAlgName.c_str(),ccent[ic],ptbins[ip],ptbins[ip+1]),50,0.,250.);
						Form("(<bkg>) jet p_{T} %sPF %s %0.0f < p_{T}^{REF} < %0.0f",kAlgName.c_str(),ccent[ic],ptbins[ip],ptbins[ip+1]),150,0.,550.);
	  havbkgpt_ind[k][ic][it][ip]->Sumw2();
	  for(int ie=0; ie<neta_wide; ie++){
	    havbkgpt_ind_eta[k][ic][it][ip][ie]= new TH1D(Form("havbkgpt_ind_eta_%sPF_%d_%d_%d_%d_%d_%d",kAlgName.c_str(),iSet,k,ic,it,ip,ie),
							  Form("(<bkg>) jet p_{T} %sPF %s %0.0f < p_{T}^{REF} < %0.0f",
							       //kAlgName.c_str(),ccent[ic],ptbins[ip],ptbins[ip+1]),50,0.,250.);
							       kAlgName.c_str(),ccent[ic],ptbins[ip],ptbins[ip+1]),150,0.,550.);
	    havbkgpt_ind_eta[k][ic][it][ip][ie]->Sumw2();
	  }
	}
      }
    }
  }


  TH1D *htrig_deta[ncen][ntrig][nbins];
  TH1D *htrig_dphi[ncen][ntrig][nbins];
  TH1D *htrig_dpt [ncen][ntrig][nbins];
  for(int ic=0;ic<ncen; ic++){
    for(int it=0;it<ntrig;it++){
      for(int ip=0; ip<nbins; ip++){
	htrig_deta[ic][it][ip] = new TH1D(Form("htrig_deta_%sPF_%d_%d_%d_%d",kAlgName.c_str(),iSet,ic,it,ip),
					  Form("htrig_deta_%sPF_%d_%d_%d_%d",kAlgName.c_str(),iSet,ic,it,ip),42,-4.2,4.2);
	htrig_deta[ic][it][ip]->Sumw2();
	htrig_dphi[ic][it][ip] = new TH1D(Form("htrig_dphi_%sPF_%d_%d_%d_%d",kAlgName.c_str(),iSet,ic,it,ip),
					  Form("htrig_dphi_%sPF_%d_%d_%d_%d",kAlgName.c_str(),iSet,ic,it,ip),72,-1.0,pi2);
	htrig_dphi[ic][it][ip]->Sumw2();
	htrig_dpt [ic][it][ip] = new TH1D(Form("htrig_dpt_%sPF_%d_%d_%d_%d",kAlgName.c_str(),iSet,ic,it,ip),
					  Form("htrig_dphi_%sPF_%d_%d_%d_%d",kAlgName.c_str(),iSet,ic,it,ip),200,-200.0,200.0);
	htrig_dpt[ic][it][ip]->Sumw2();
      }
    }
  }

  for(int k=0; k<2; k++){
    for(int ic=0; ic<ncen; ic++){
      for(int ip=0; ip<nbins; ip++){
	havbkgpt_comb[k][ic][ip]= new TH1D(Form("havbkgpt_comb_%sPF_%d_%d_%d_%d",kAlgName.c_str(),iSet,k,ic,ip),
					   //Form("(<bkg>) jet p_{T} %sPF %s %0.0f < p_{T}^{REF} < %0.0f",kAlgName.c_str(),ccent[ic],ptbins[ip],ptbins[ip+1]),50,0.,250.);
					   Form("(<bkg>) jet p_{T} %sPF %s %0.0f < p_{T}^{REF} < %0.0f",kAlgName.c_str(),ccent[ic],ptbins[ip],ptbins[ip+1]),150,0.,550.);
	havbkgpt_comb[k][ic][ip]->Sumw2();
	for(int ie=0; ie<neta_wide; ie++){
	  havbkgpt_comb_eta[k][ic][ip][ie]= new TH1D(Form("havbkgpt_comb_eta_%sPF_%d_%d_%d_%d_%d",kAlgName.c_str(),iSet,k,ic,ip,ie),
						     Form("(<bkg>) jet p_{T} %sPF %s %0.0f < p_{T}^{REF} < %0.0f",
							  //kAlgName.c_str(),ccent[ic],ptbins[ip],ptbins[ip+1]),50,0.,250.);
							  kAlgName.c_str(),ccent[ic],ptbins[ip],ptbins[ip+1]),150,0.,550.);
	  havbkgpt_comb_eta[k][ic][ip][ie]->Sumw2();
	}
	
	hrecopt_genm [k][ic][ip]= new TH1D(Form("hrecopt_genm%d_%d_%d",k,ic,ip),
					   Form("(Recopt ) jet p_{T} %s %s %0.0f < p_{T}^{REF} < %0.0f",kAlgName.c_str(),ccent[ic],ptbins[ip],ptbins[ip+1]),
					   200,0,1000);
	hrecopt_genm [k][ic][ip]->Sumw2();
	hrescrpt_genm [k][ic][ip]= new TH1D(Form("hrescrpt_genm%d_%d_%d",k,ic,ip),
					    Form("(Reco/Gen) jet p_{T} %s %s %0.0f < p_{T}^{REF} < %0.0f",kAlgName.c_str(),ccent[ic],ptbins[ip],ptbins[ip+1]),
					    rbins,rbinl,rbinh);
	hrescrpt_genm [k][ic][ip]->Sumw2();
	hrescrpt_sm_genm [k][ic][ip]= new TH1D(Form("hrescrpt_sm_genm%d_%d_%d",k,ic,ip),
					       Form("(smear Reco/Gen) jet p_{T} %s %s %0.0f < p_{T}^{REF} < %0.0f",kAlgName.c_str(),ccent[ic],ptbins[ip],ptbins[ip+1]),
					       rbins,rbinl,rbinh);
	hrescrpt_sm_genm [k][ic][ip]->Sumw2();
	hresrrpt_genm [k][ic][ip]= new TH1D(Form("hresrrpt_genm%d_%d_%d",k,ic,ip),
					    Form("(Raw/Gen) jet p_{T} %s %s %0.0f < p_{T}^{REF} < %0.0f",kAlgName.c_str(),ccent[ic],ptbins[ip],ptbins[ip+1]),
					    rbins,rbinl,rbinh);
	hresrrpt_genm [k][ic][ip]->Sumw2();

	
	hrescrpt_q_genm [k][ic][ip]= new TH1D(Form("hrescrpt_q_genm%d_%d_%d",k,ic,ip),
					    Form("(Reco/Gen) quark jet p_{T} %s %s %0.0f < p_{T}^{REF} < %0.0f",kAlgName.c_str(),ccent[ic],ptbins[ip],ptbins[ip+1]),
					    rbins,rbinl,rbinh);
	hrescrpt_q_genm [k][ic][ip]->Sumw2();
	hresrrpt_q_genm [k][ic][ip]= new TH1D(Form("hresrrpt_q_genm%d_%d_%d",k,ic,ip),
					    Form("(Raw/Gen) quark jet p_{T} %s %s %0.0f < p_{T}^{REF} < %0.0f",kAlgName.c_str(),ccent[ic],ptbins[ip],ptbins[ip+1]),
					    rbins,rbinl,rbinh);
	hresrrpt_q_genm [k][ic][ip]->Sumw2();

	hrescrpt_g_genm [k][ic][ip]= new TH1D(Form("hrescrpt_g_genm%d_%d_%d",k,ic,ip),
					    Form("(Reco/Gen) gluon jet p_{T} %s %s %0.0f < p_{T}^{REF} < %0.0f",kAlgName.c_str(),ccent[ic],ptbins[ip],ptbins[ip+1]),
					    rbins,rbinl,rbinh);
	hrescrpt_g_genm [k][ic][ip]->Sumw2();
	hresrrpt_g_genm [k][ic][ip]= new TH1D(Form("hresrrpt_g_genm%d_%d_%d",k,ic,ip),
					      Form("(Raw/Gen) gluon jet p_{T} %s %s %0.0f < p_{T}^{REF} < %0.0f",kAlgName.c_str(),ccent[ic],ptbins[ip],ptbins[ip+1]),
					      rbins,rbinl,rbinh);
	hresrrpt_g_genm [k][ic][ip]->Sumw2();

	hrescrpt_u_genm [k][ic][ip]= new TH1D(Form("hrescrpt_u_genm%d_%d_%d",k,ic,ip),
					    Form("(Reco/Gen) u-quark jet p_{T} %s %s %0.0f < p_{T}^{REF} < %0.0f",kAlgName.c_str(),ccent[ic],ptbins[ip],ptbins[ip+1]),
					    rbins,rbinl,rbinh);
	hrescrpt_u_genm [k][ic][ip]->Sumw2();
	hresrrpt_u_genm [k][ic][ip]= new TH1D(Form("hresrrpt_u_genm%d_%d_%d",k,ic,ip),
					    Form("(Raw/Gen) u-quark jet p_{T} %s %s %0.0f < p_{T}^{REF} < %0.0f",kAlgName.c_str(),ccent[ic],ptbins[ip],ptbins[ip+1]),
					    rbins,rbinl,rbinh);
	hresrrpt_u_genm [k][ic][ip]->Sumw2();

	hrescrpt_d_genm [k][ic][ip]= new TH1D(Form("hrescrpt_d_genm%d_%d_%d",k,ic,ip),
					      Form("(Reco/Gen) d-quark jet p_{T} %s %s %0.0f < p_{T}^{REF} < %0.0f",kAlgName.c_str(),ccent[ic],ptbins[ip],ptbins[ip+1]),
					      rbins,rbinl,rbinh);
	hrescrpt_d_genm [k][ic][ip]->Sumw2();
	hresrrpt_d_genm [k][ic][ip]= new TH1D(Form("hresrrpt_d_genm%d_%d_%d",k,ic,ip),
					    Form("(Raw/Gen) d-quark jet p_{T} %s %s %0.0f < p_{T}^{REF} < %0.0f",kAlgName.c_str(),ccent[ic],ptbins[ip],ptbins[ip+1]),
					    rbins,rbinl,rbinh);
	hresrrpt_d_genm [k][ic][ip]->Sumw2();

	hrescrpt_s_genm [k][ic][ip]= new TH1D(Form("hrescrpt_s_genm%d_%d_%d",k,ic,ip),
					    Form("(Reco/Gen) s-quark jet p_{T} %s %s %0.0f < p_{T}^{REF} < %0.0f",kAlgName.c_str(),ccent[ic],ptbins[ip],ptbins[ip+1]),
					    rbins,rbinl,rbinh);
	hrescrpt_s_genm [k][ic][ip]->Sumw2();
	hresrrpt_s_genm [k][ic][ip]= new TH1D(Form("hresrrpt_s_genm%d_%d_%d",k,ic,ip),
					    Form("(Raw/Gen) s-quark jet p_{T} %s %s %0.0f < p_{T}^{REF} < %0.0f",kAlgName.c_str(),ccent[ic],ptbins[ip],ptbins[ip+1]),
					    rbins,rbinl,rbinh);
	hresrrpt_s_genm [k][ic][ip]->Sumw2();


	for(int ie=0;ie<neta;ie++){
          hrescrpt_genm_eta [k][ic][ip][ie]= new TH1D(Form("hrescrpt_genm_eta%d_%d_%d_%d",k,ic,ip,ie),
						      Form("(Reco/Gen) jet p_{T} %s %s %d %d",kAlgName.c_str(),ccent[ic],ip,ie),
						      rbins,rbinl,rbinh);
          hrescrpt_genm_eta [k][ic][ip][ie]->Sumw2();
          hresrrpt_genm_eta [k][ic][ip][ie]= new TH1D(Form("hresrrpt_genm_eta%d_%d_%d_%d",k,ic,ip,ie),
						      Form("(Raw/Gen) jet p_{T} %s %s %d %d",kAlgName.c_str(),ccent[ic],ip,ie),
						      rbins,rbinl,rbinh);
          hresrrpt_genm_eta [k][ic][ip][ie]->Sumw2();


          hrescrpt_q_genm_eta [k][ic][ip][ie]= new TH1D(Form("hrescrpt_q_genm_eta%d_%d_%d_%d",k,ic,ip,ie),
						      Form("(Reco/Gen) quark jet p_{T} %s %s %d %d",kAlgName.c_str(),ccent[ic],ip,ie),
						      rbins,rbinl,rbinh);
          hrescrpt_q_genm_eta [k][ic][ip][ie]->Sumw2();
          hresrrpt_q_genm_eta [k][ic][ip][ie]= new TH1D(Form("hresrrpt_q_genm_eta%d_%d_%d_%d",k,ic,ip,ie),
						      Form("(Raw/Gen) quark jet p_{T} %s %s %d %d",kAlgName.c_str(),ccent[ic],ip,ie),
						      rbins,rbinl,rbinh);
          hresrrpt_q_genm_eta [k][ic][ip][ie]->Sumw2();


          hrescrpt_g_genm_eta [k][ic][ip][ie]= new TH1D(Form("hrescrpt_g_genm_eta%d_%d_%d_%d",k,ic,ip,ie),
							Form("(Reco/Gen) gluon jet p_{T} %s %s %d %d",kAlgName.c_str(),ccent[ic],ip,ie),
							rbins,rbinl,rbinh);
          hrescrpt_g_genm_eta [k][ic][ip][ie]->Sumw2();
          hresrrpt_g_genm_eta [k][ic][ip][ie]= new TH1D(Form("hresrrpt_g_genm_eta%d_%d_%d_%d",k,ic,ip,ie),
							Form("(Raw/Gen) gluon jet p_{T} %s %s %d %d",kAlgName.c_str(),ccent[ic],ip,ie),
							rbins,rbinl,rbinh);
          hresrrpt_g_genm_eta [k][ic][ip][ie]->Sumw2();
	  
	  
          hrescrpt_u_genm_eta [k][ic][ip][ie]= new TH1D(Form("hrescrpt_u_genm_eta%d_%d_%d_%d",k,ic,ip,ie),
							Form("(Reco/Gen) uquark  jet p_{T} %s %s %d %d",kAlgName.c_str(),ccent[ic],ip,ie),
							rbins,rbinl,rbinh);
          hrescrpt_u_genm_eta [k][ic][ip][ie]->Sumw2();
          hresrrpt_u_genm_eta [k][ic][ip][ie]= new TH1D(Form("hresrrpt_u_genm_eta%d_%d_%d_%d",k,ic,ip,ie),
							Form("(Raw/Gen) u-quark jet p_{T} %s %s %d %d",kAlgName.c_str(),ccent[ic],ip,ie),
							rbins,rbinl,rbinh);
          hresrrpt_u_genm_eta [k][ic][ip][ie]->Sumw2();

          hrescrpt_d_genm_eta [k][ic][ip][ie]= new TH1D(Form("hrescrpt_d_genm_eta%d_%d_%d_%d",k,ic,ip,ie),
							Form("(Reco/Gen) d-quark  jet p_{T} %s %s %d %d",kAlgName.c_str(),ccent[ic],ip,ie),
							rbins,rbinl,rbinh);
          hrescrpt_d_genm_eta [k][ic][ip][ie]->Sumw2();
          hresrrpt_d_genm_eta [k][ic][ip][ie]= new TH1D(Form("hresrrpt_d_genm_eta%d_%d_%d_%d",k,ic,ip,ie),
							Form("(Raw/Gen) d-quark jet p_{T} %s %s %d %d",kAlgName.c_str(),ccent[ic],ip,ie),
							rbins,rbinl,rbinh);
          hresrrpt_d_genm_eta [k][ic][ip][ie]->Sumw2();


          hrescrpt_s_genm_eta [k][ic][ip][ie]= new TH1D(Form("hrescrpt_s_genm_eta%d_%d_%d_%d",k,ic,ip,ie),
							Form("(Reco/Gen) s-quark  jet p_{T} %s %s %d %d",kAlgName.c_str(),ccent[ic],ip,ie),
							rbins,rbinl,rbinh);
          hrescrpt_s_genm_eta [k][ic][ip][ie]->Sumw2();
          hresrrpt_s_genm_eta [k][ic][ip][ie]= new TH1D(Form("hresrrpt_s_genm_eta%d_%d_%d_%d",k,ic,ip,ie),
							Form("(Raw/Gen) s-quark jet p_{T} %s %s %d %d",kAlgName.c_str(),ccent[ic],ip,ie),
							rbins,rbinl,rbinh);
          hresrrpt_s_genm_eta [k][ic][ip][ie]->Sumw2();
	}
	
        for(int ij=0;ij<nphi;ij++){
          hrescrpt_genm_phi [k][ic][ip][ij]= new TH1D(Form("hrescrpt_genm_phi%d_%d_%d_%d",k,ic,ip,ij),
						      Form("(Reco/Gen) jet p_{T} %s %s %d %d",kAlgName.c_str(),ccent[ic],ip,ij),
						      rbins,rbinl,rbinh);
          hrescrpt_genm_phi [k][ic][ip][ij]->Sumw2();
          hresrrpt_genm_phi [k][ic][ip][ij]= new TH1D(Form("hresrrpt_genm_phi%d_%d_%d_%d",k,ic,ip,ij),
						      Form("(Raw/Gen) jet p_{T} %s %s %d %d",kAlgName.c_str(),ccent[ic],ip,ij),
						      rbins,rbinl,rbinh);
          hresrrpt_genm_phi [k][ic][ip][ij]->Sumw2();
        }
      }

      //! coarse pt bin
      for(int ip=0;ip<nbins_wide;ip++){
        for(int ie=0;ie<neta;ie++){
          hrescrpt_wide_genm_eta [k][ic][ip][ie]= new TH1D(Form("hrescrpt_wide_genm_eta%d_%d_%d_%d",k,ic,ip,ie),Form("(Reco/Gen) jet p_{T} %s %d %d",ccent[ic],ip,ie),rbins,rbinl,rbinh);
          hrescrpt_wide_genm_eta [k][ic][ip][ie]->Sumw2();
          hresrrpt_wide_genm_eta [k][ic][ip][ie]= new TH1D(Form("hresrrpt_wide_genm_eta%d_%d_%d_%d",k,ic,ip,ie),Form("(Raw/Gen) jet p_{T} %s %d %d",ccent[ic],ip,ie),rbins,rbinl,rbinh);
          hresrrpt_wide_genm_eta [k][ic][ip][ie]->Sumw2();
        }
        for(int ij=0;ij<nphi;ij++){
          hrescrpt_wide_genm_phi [k][ic][ip][ij]= new TH1D(Form("hrescrpt_wide_genm_phi%d_%d_%d_%d",k,ic,ip,ij),Form("(Reco/Gen) jet p_{T} %s %d %d",ccent[ic],ip,ij),rbins,rbinl,rbinh);
          hrescrpt_wide_genm_phi [k][ic][ip][ij]->Sumw2();
          hresrrpt_wide_genm_phi [k][ic][ip][ij]= new TH1D(Form("hresrrpt_wide_genm_phi%d_%d_%d_%d",k,ic,ip,ij),Form("(Raw/Gen) jet p_{T} %s %d %d",ccent[ic],ip,ij),rbins,rbinl,rbinh);
          hresrrpt_wide_genm_phi [k][ic][ip][ij]->Sumw2();
	}
      }
    }
  }

  fout->cd("../");
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout<<"Initialized the histograms " <<std::endl;

  Long64_t nbytes=0;
  Long64_t nentries = tch_pfjet->GetEntries();
  std::cout<<Form("# of entries in TTree for %s %s : ",kAlgName.c_str(),kSpecies.c_str())<<nentries<<std::endl;

  double wxs=1.;
  if( kDataset == "mc" ){
    if( kSpecies == "pbpb" || kSpecies == "pp" ){
    }else if( kSpecies == "pbpb_mb" )wxs=1.0;
  }else wxs=1.0;
  //return 0;

  //! Start event loop
  for (Long64_t i=0; i<nentries;i++) {
    nbytes += tch_pfjet->GetEntry(i);
    if(printDebug && i==10)break;
    //if( i == 10 )break; 
    float rndm=gRandom->Rndm();
    

    int TrigFill=0;
    float prescl=1.0;
    int iTrig=-1;
    if( kSpecies == "pp" ){
      if( jet40==1  && jet60==0 && jet80==0 ){TrigFill=1;iTrig=0;prescl=Jet40_prescl;}
      else if( jet60==1 && jet80==0 ){TrigFill=1;iTrig=1;prescl=1.0;}
      else if( jet80==1 ){TrigFill=1;iTrig=2;prescl=1.0;}
    }else{
      if( jet55==1  && jet65==0 && jet80==0 ){TrigFill=1; iTrig=0; prescl=Jet55_prescl;}
      else if( jet65==1 && jet80==0 ){TrigFill=1; iTrig=1; prescl=1.0;}
      else if( jet80==1 ){TrigFill=1; iTrig=2; prescl=1.0;}
    }
    if( iTrig < 0 )continue;

    hEvents_Total->Fill(rndm);
    if( kSpecies == "pp" ){
      if(jet40)hEvents_jet40->Fill(rndm);
      if(jet40_prescl)hEvents_jet40_prescl->Fill(jet40_prescl,rndm);
      if(jet60)hEvents_jet60->Fill(rndm);
      if(jet60_prescl)hEvents_jet60_prescl->Fill(jet60_prescl,rndm);
      if(jet80)hEvents_jet80->Fill(rndm);
      if(jet80_prescl)hEvents_jet80_prescl->Fill(jet80_prescl,rndm);

      if(jet40==1 && jet60==0 && jet80==0)hEvents_jet40_nojet60_nojet80->Fill(rndm);
      if(jet60==1 && jet80==0)hEvents_jet60_nojet80->Fill(rndm);
    }else{
      if(jet55)hEvents_jet55->Fill(rndm);
      if(jet55_prescl)hEvents_jet55_prescl->Fill(jet55_prescl,rndm);
      if(jet65)hEvents_jet65->Fill(rndm);
      if(jet65_prescl)hEvents_jet65_prescl->Fill(jet65_prescl,rndm);
      if(jet80)hEvents_jet80->Fill(rndm);
      if(jet80_prescl)hEvents_jet80_prescl->Fill(jet80_prescl,rndm);

      if(jet55==1 && jet65==0 && jet80==0)hEvents_jet55_nojet65_nojet80->Fill(rndm);
      if(jet65==1 && jet80==0)hEvents_jet65_nojet80->Fill(rndm);
    }

    if( kDataset == "data" ){
      if(pcollisionEventSelection)hEvents_pCollEvent->Fill(rndm);
      if(pcollisionEventSelection && pHBHENoiseFilter)hEvents_pHBHENoise->Fill(rndm);
      if(pcollisionEventSelection && pHBHENoiseFilter && fabs(vz)<kvzcut )hEvents_Vzcut->Fill(rndm);
      if(pcollisionEventSelection==0 || pHBHENoiseFilter==0 || fabs(vz) > kvzcut){
	hEvents_bad->Fill(rndm);
	continue;
      }
    }else if( kDataset == "mc" ){//! HBHENoiseFilter
      if(pcollisionEventSelection)hEvents_pCollEvent->Fill(rndm);
      if(pcollisionEventSelection)hEvents_pHBHENoise->Fill(rndm);
      if(pcollisionEventSelection && fabs(vz)<kvzcut )hEvents_Vzcut->Fill(rndm);
      if(pcollisionEventSelection==0 || fabs(vz) > kvzcut){
	hEvents_bad->Fill(rndm);
	continue;
      }
    }
    
    int iCent = -1;
    if( kSpecies == "pbpb" || kSpecies == "pbpb_mb")iCent = GetCentBin(hiBin);
    else iCent = ncen-1;

    if(iCent<0 || iCent>=ncen)continue;

    if( kDataset == "data" ){
      if ( kSpecies == "pbpb" || kSpecies == "pbpb_mb" ){

	//! pileupcut
	// if((run_value > 182323) && ((float)hiZDC/1350 >= (400 - (float)0.091 * hiHF))){
	//   hEvents_pileup->Fill(rndm);
	//   continue;
	// }
	// if((run_value <= 182323) && ((float)hiZDCminus/1350 >= (65 - (float)0.0155 * hiHF))){
	//   hEvents_pileup->Fill(rndm);
	//   continue;        
	// }

	if((float)hiZDCminus/1350 >= (90 - 0.0204 * hiHF)) continue;


	//! SuperNovae events use calo jets
	//int lowjetCounter=0;
	int jetCounter=0;
	for(int g = 0;g<nref_calo;g++){
	  // if( fabs(jteta_calo[g]) < 2. && jtpt_calo[g]>=kptrecocut ){ //to select inside
	  // 	lowjetCounter++;
	  // }
	  if( fabs(jteta_calo[g]) < ketacut && jtpt_calo[g]>=50. ){ //to select inside
	    jetCounter++;
	  }//eta selection cut
	}// jet loop
	// apply the correct supernova selection cut rejection here:
	if( hiNpix > 38000 - 500*jetCounter ){
	  hEvents_supernova->Fill(rndm);
	  continue;
	}
      }
    }

    if( nref_pf==0 && nref_calo==0 ){
      hEvents_nopfcalo->Fill(rndm);
      continue;
    }
    //if( lowjetCounter == 0 )continue;

    //!  pt hat weight factor taking all events
    wxs=1.0;
    if( kDataset == "mc" ){
      wxs = GetXsecWt( pthat, kSpecies );
      hpthat->Fill( pthat, wxs );
    }
    hEvents_maxpthat->Fill(rndm);

    //! ----------------------------------------------------------------------------
    // int iStat=0;
    // for(int pj=0; pj<nref; pj++){ //! PFjet
    //   if( pfrefdrjt[pj] > 0.3 && sid[pj]==0 && pfrefpt[pj] > 50. && fabs(jteta[pj]) < ketacut){
    // 	std::cout << "\t\t Found an interesting cand : "<< i << " dr : " << pfrefdrjt[pj] << " refpt : " << pfrefpt[pj] << " recopt : " << jtpt[pj] << std::endl;
    // 	iStat=1;
    //   }
    // }
    // if( iStat )printDebug=true;
    // else printDebug=false;
    //! ----------------------------------------------------------------------------    

    if(printDebug){
      std::cout << "------------------------------------Start Event # " << i <<"------------------------------------------------------------------ " << std::endl;
      std::cout<<" ***** Event # " <<i<<" "<<run_value<<" " <<evt_value<<" " <<lumi_value<<"\t vz  : "<<vz<<"\t hiBin : "<<hiBin<<"\t  # of PF jets "<<nref_pf<<" # of calojets  "<<nref_calo<<std::endl;
    }
    //if(i%50000==0)std::cout<<" ********** Event # " <<i<<"\t vz  : "<<vz<<"\t hiBin : "<<hiBin<<"\t  # of PF jets "<<nref<< "  # of Calo jets " <<nref_calo<<std::endl;
    
    int pfjets=0;
    int calojets=0;
    
    std::vector < Jet > vPFJets, vCaloJets;
    std::vector < int > pfid(nref_pf), caloid(nref_calo);
    
    if(printDebug)std::cout << " PF jets : " << std::endl;
 
    for(int pj=0; pj<nref_pf; pj++){ //! PFjet
      //if( rawpt[pj] < kptrawcut || jtpt[pj] < kptrecocut ) continue;
      //if( fabs(jteta[pj]) > ketacut ) continue;

      // if( kSpecies == "pbpb" || kSpecies == "pbpb_mb" ){
      // 	if( trackMax[pj]/jtpt[pj] > 0.02)
      // 	  if( muonMax[pj]/(neutralMax[pj]+photonMax[pj]+chargedMax[pj]+muonMax[pj]+eleMax[pj])<0.975 )
      // }else if( kSpecies == "pp" ){
      // }

      //if( pj == 53 )cout << " jet pT "<< jtpt[pj] << endl;

      Jet pfj;
      pfj.id  = pj;
      pfj.eta = jteta[pj];
      pfj.phi = jtphi[pj];
      pfj.pt  = jtpt [pj];

      //if( abs(pfrefparton_flavor[pj]) <= 21 )cout << pj << "  ref parton : " << pfrefparton_flavor[pj] << endl;


      if(printDebug){
	if( kDataset == "mc" )std::cout <<"\t" << pj << "  pt : " << jtpt[pj] << " eta : "  << jteta[pj] << " phi : " << jtphi[pj] << " subid : " << sid[pj] << " dr : " << pfrefdrjt[pj] << std::endl;
	else std::cout <<"\t"<< pj << "  pt : " << jtpt[pj] << " eta : "  << jteta[pj] << " phi : " << jtphi[pj] << std::endl;
      }
      vPFJets.push_back(pfj);
      pfjets++;
    }    
    
    if(printDebug){
      std::cout << std::endl;
      std::cout << " Calo jets : " << std::endl;
    }
    for(int cj=0; cj<nref_calo; cj++){ //! CaloJet
      //if( rawpt_calo[cj] < kptrawcut || jtpt_calo[cj] < kptrecocut) continue;
      //if( fabs(jteta_calo[cj]) > ketacut ) continue;
      
      Jet clj;
      clj.id  = cj;
      clj.eta = jteta_calo[cj];
      clj.phi = jtphi_calo[cj];
      clj.pt  = jtpt_calo[cj];

      if(printDebug){
	if( kDataset == "mc" )std::cout <<"\t" << cj << "  pt : " << jtpt_calo[cj] << " eta : "  << jteta_calo[cj] << " phi : " << jtphi_calo[cj] << " subid : " << sid_calo[cj] << " dr : " << refdrjt_calo[cj] << std::endl;
	else  std::cout <<"\t" << cj << "  pt : " << jtpt_calo[cj] << " eta : "  << jteta_calo[cj] << " phi : " << jtphi_calo[cj] << std::endl;
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
      hEvents_nopfcalo->Fill(rndm);
      continue;
    }

    hEvents_Good->Fill(rndm);
    hVz->Fill(vz);
    hBin->Fill(hiBin);

    bool onlyCalo   = (pfjets==0 && calojets >0) ? true : false;
    bool onlyPF     = (pfjets>0  && calojets==0) ? true : false;
    bool bothPFCalo = (pfjets>0  && calojets >0) ? true : false;

    int matchedJets=0;
    int unmatchedPFJets=0;
    int unmatchedCaloJets=0;
    int multimatchedPFJets=0;

    //! for jetTree
    double weight = 1.;
    if( kDataset == "mc" ){
      double wvz  = fVz->Eval(vz);
      double wcen = hCentWeight->GetBinContent(hCentWeight->FindBin(hiBin));
      if ( kSpecies == "pbpb" )weight = (wxs*wvz*wcen);
      else if( kSpecies == "pbpb_mb" )weight = (wvz*wcen);
      else weight = (wxs*wvz);
    }else weight = 1.;

    if(printDebug)std::cout <<" Total ==>  # of PF jets : " << pfjets << ", # of calojets : "  << calojets <<"\n"<<std::endl;

    std::vector < Jet >::const_iterator iJet, jJet;

    if( onlyPF ){
      
      for(iJet = vPFJets.begin(); iJet != vPFJets.end(); ++iJet){ //! PFjet

	int pj = (*iJet).id; 

	Int_t PFElecCut=0;
	Int_t MuCut=0;
	Int_t TrCut=0;
	if( kSpecies == "pbpb" || kSpecies == "pbpb_mb" ){
	  Float_t Sumcand  = chargedSum[pj] +  photonSum[pj] + neutralSum[pj] + muonSum[pj];
	  if( eleMax[pj]/Sumcand < 0.05 )PFElecCut=1;
	  //if( muonMax[pj]/(neutralMax[pj]+photonMax[pj]+chargedMax[pj]+muonMax[pj]+eleMax[pj])<0.975 )MuCut = 1;
	  //if( trackMax[pj]/jtpt[pj] > 0.02)TrCut=1;
	  MuCut=1;
	  TrCut=1;
	}else if( kSpecies == "pp" ){
	  PFElecCut=1;
	  MuCut=1;
	  TrCut=1;
	}

	bool isel = true;	
	if( kDataset=="mc"){
	  if(sid[pj] != 0)isel=false;
	  if(pfrefdrjt[pj] > kdelrcut)isel=false;
	  //if( fabs(pfrefpt[pj]) > 2.*pthat )isel=false;
	  if( fabs(pfrefeta[pj]) > ketacut )isel=false;

	  if( isel ){
	    double resp_corr =  jtpt[pj]  / pfrefpt[pj];
	    double resp_raw  =  rawpt[pj] / pfrefpt[pj];

	    int iphi     = GetPhiBin(pfrefphi[pj]);
	    int ieta     = GetEtaBin(pfrefeta[pj]);
	    int ipt_mc   = GetPtBin (pfrefpt[pj]);
	    int ipt_wide = GetPtBinWide (pfrefpt[pj]);

	    if( ipt_mc>=0 && ipt_mc<nbins ){
	      hrecopt_genm  [0][iCent][ipt_mc]->Fill(jtpt[pj],weight*prescl);
	      hrescrpt_genm [0][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
	      hresrrpt_genm [0][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);

	      if( iCent == ncen-1 ){ //! smearing only for pp
		for(int ic=0; ic< ncen; ic++){
		  double smpt = GetSmearedPt(ic,jtpt[pj], pfrefpt[pj]);
		  double resp_smear=  smpt/pfrefpt[pj];
		  hrescrpt_sm_genm [0][ic][ipt_mc]->Fill(resp_smear,weight*prescl);
		}
	      }

	      //! Flavor dependence
	      if( abs(pfrefparton_flavor[pj]) < 5){
		hrescrpt_q_genm [0][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		hresrrpt_q_genm [0][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);

		if( abs(pfrefparton_flavor[pj]) == 1 ){
		  hrescrpt_d_genm [0][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		  hresrrpt_d_genm [0][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);
		}else if( abs(pfrefparton_flavor[pj]) == 2 ){
		  hrescrpt_u_genm [0][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		  hresrrpt_u_genm [0][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);
		}else if( abs(pfrefparton_flavor[pj]) == 3 ){
		  hrescrpt_s_genm [0][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		  hresrrpt_s_genm [0][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);
		}
	      }else if( abs(pfrefparton_flavor[pj]) == 21 ){
		hrescrpt_g_genm [0][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		hresrrpt_g_genm [0][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);
		//cout << " ONLY PF found a gluon jet : " << resp_corr << endl;
	      }

	      if( ieta>=0 && ieta<neta ){
		hrescrpt_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
		hresrrpt_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);

		//! Flavor dependence
		if( abs(pfrefparton_flavor[pj]) < 5){
		  hrescrpt_q_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
		  hresrrpt_q_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		  
		  if( abs(pfrefparton_flavor[pj]) == 1 ){
		    hrescrpt_d_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
		    hresrrpt_d_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		  }else if( abs(pfrefparton_flavor[pj]) == 2 ){
		    hrescrpt_u_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
		    hresrrpt_u_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		  }else if( abs(pfrefparton_flavor[pj]) == 3 ){
		    hrescrpt_s_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
		    hresrrpt_s_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		  }
		}else if( abs(pfrefparton_flavor[pj]) == 21 ){
		  hrescrpt_g_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
		  hresrrpt_g_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		}


		if( ipt_wide>=0 && ipt_wide <nbins_wide){
		  hrescrpt_wide_genm_eta [0][iCent][ipt_wide][ieta]->Fill(resp_corr,weight*prescl);
		  hresrrpt_wide_genm_eta [0][iCent][ipt_wide][ieta]->Fill(resp_raw ,weight*prescl);
		}
	      }
	      if( iphi>=0 && iphi<nphi ){
		hrescrpt_genm_phi [0][iCent][ipt_mc][iphi]->Fill(resp_corr,weight*prescl);
		hresrrpt_genm_phi [0][iCent][ipt_mc][iphi]->Fill(resp_raw ,weight*prescl);
		if( ipt_wide>=0 && ipt_wide<nbins_wide){
		  hrescrpt_wide_genm_phi [0][iCent][ipt_wide][iphi]->Fill(resp_corr,weight*prescl);
		  hresrrpt_wide_genm_phi [0][iCent][ipt_wide][iphi]->Fill(resp_raw ,weight*prescl);
		}
	      }
	    }
	    if( (PFElecCut*MuCut*TrCut) == 1 ){
	      if( ipt_mc>=0 && ipt_mc<nbins ){
		hrecopt_genm  [1][iCent][ipt_mc]->Fill(jtpt[pj],weight*prescl);
		hrescrpt_genm [1][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		hresrrpt_genm [1][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);
		
		if( iCent == ncen-1 ){ //! smearing only for pp
		  for(int ic=0; ic< ncen; ic++){
		    double smpt = GetSmearedPt(ic,jtpt[pj],pfrefpt[pj]);
		    double resp_smear=  smpt/pfrefpt[pj];
		    hrescrpt_sm_genm [1][ic][ipt_mc]->Fill(resp_smear,weight*prescl);
		  }
		}

		//! Flavor dependence
		if( abs(pfrefparton_flavor[pj]) < 5){
		  hrescrpt_q_genm [1][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		  hresrrpt_q_genm [1][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);
		  
		  if( abs(pfrefparton_flavor[pj]) == 1 ){
		    hrescrpt_d_genm [1][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		    hresrrpt_d_genm [1][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);
		  }else if( abs(pfrefparton_flavor[pj]) == 2 ){
		    hrescrpt_u_genm [1][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		    hresrrpt_u_genm [1][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);
		  }else if( abs(pfrefparton_flavor[pj]) == 3 ){
		    hrescrpt_s_genm [1][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		    hresrrpt_s_genm [1][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);
		  }
		}else if( abs(pfrefparton_flavor[pj]) == 21 ){
		  hrescrpt_g_genm [1][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		  hresrrpt_g_genm [1][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);
		  //cout << " ONLY PF PFElec cut found a gluon jet : " << resp_corr << endl;
		}


		if( ieta>=0 && ieta<neta ){
		  hrescrpt_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
		  hresrrpt_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);

		  //! Flavor dependence
		  if( abs(pfrefparton_flavor[pj]) < 5){
		    hrescrpt_q_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
		    hresrrpt_q_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		  
		    if( abs(pfrefparton_flavor[pj]) == 1 ){
		      hrescrpt_d_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
		      hresrrpt_d_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		    }else if( abs(pfrefparton_flavor[pj]) == 2 ){
		      hrescrpt_u_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
		      hresrrpt_u_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		    }else if( abs(pfrefparton_flavor[pj]) == 3 ){
		      hrescrpt_s_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
		      hresrrpt_s_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		    }
		  }else if( abs(pfrefparton_flavor[pj]) == 21 ){
		    hrescrpt_g_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
		    hresrrpt_g_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		  }

		  if( ipt_wide>=0 && ipt_wide <nbins_wide){
		    hrescrpt_wide_genm_eta [1][iCent][ipt_wide][ieta]->Fill(resp_corr,weight*prescl);
		    hresrrpt_wide_genm_eta [1][iCent][ipt_wide][ieta]->Fill(resp_raw ,weight*prescl);
		  }
		}
		if( iphi>=0 && iphi<nphi ){
		  hrescrpt_genm_phi [1][iCent][ipt_mc][iphi]->Fill(resp_corr,weight*prescl);
		  hresrrpt_genm_phi [1][iCent][ipt_mc][iphi]->Fill(resp_raw ,weight*prescl);
		  if( ipt_wide>=0 && ipt_wide<nbins_wide){
		    hrescrpt_wide_genm_phi [1][iCent][ipt_wide][iphi]->Fill(resp_corr,weight*prescl);
		    hresrrpt_wide_genm_phi [1][iCent][ipt_wide][iphi]->Fill(resp_raw ,weight*prescl);
		  }
		}
	      }
	    }
	  }
	}

	float sumPF = (chargedSum[pj] + neutralSum[pj] + photonSum[pj] + muonSum[pj] + eleSum[pj]);
	float bkgd  = sumPF - rawpt[pj];

	int ipt = GetPtBin(jtpt[pj]);
	int ieta_wide = GetEtaBinWide(fabs(jteta[pj]));
	
	double trig_deta = 0;
	double trig_dphi = 0;
	double trig_dpt  = 0;
	if( kDataset == "data" ){
	  trig_deta = trgObj_eta - jteta[pj];
	  trig_dphi = delPhi(trgObj_phi,jtphi[pj]);
	  trig_dpt  = trgObj_pt  - jtpt[pj];
	}

	if( TrigFill && isel ){
	  if(ipt >=0 && ipt < nbins){
	    if(ieta_wide >= 0 && ieta_wide < neta_wide){

	      htrig_deta[iCent][iTrig][ipt]->Fill(trig_deta);
	      htrig_dphi[iCent][iTrig][ipt]->Fill(trig_dphi);
	      htrig_dpt [iCent][iTrig][ipt]->Fill(trig_dpt);

	      havbkgpt_comb    [0][iCent][ipt]->Fill(bkgd, weight*prescl);
	      havbkgpt_comb_eta[0][iCent][ipt][ieta_wide]->Fill(bkgd, weight*prescl);
	      havbkgpt_ind     [0][iCent][iTrig][ipt]->Fill(bkgd, weight*prescl);
	      havbkgpt_ind_eta [0][iCent][iTrig][ipt][ieta_wide]->Fill(bkgd, weight*prescl);
	      if( (PFElecCut*MuCut*TrCut) == 1 ){
		havbkgpt_comb    [1][iCent][ipt]->Fill(bkgd, weight*prescl);
		havbkgpt_comb_eta[1][iCent][ipt][ieta_wide]->Fill(bkgd, weight*prescl);
		havbkgpt_ind     [1][iCent][iTrig][ipt]->Fill(bkgd, weight*prescl);
		havbkgpt_ind_eta [1][iCent][iTrig][ipt][ieta_wide]->Fill(bkgd, weight*prescl);
	      }
	    }
	  }
	}
  
	if(printDebug){
	  if( kDataset == "mc" )std::cout <<" unmatched PF jets w ncalo=0 : " << unmatchedPFJets << "  pt : " << jtpt[pj] << " eta : "  << jteta[pj] << " phi : " << jtphi[pj] << " dr : " << pfrefdrjt[pj] <<std::endl;
	  else std::cout <<" unmatched PF jets w ncalo=0 : " << unmatchedPFJets << "  pt : " << jtpt[pj] << " eta : "  << jteta[pj] << " phi : " << jtphi[pj] << std::endl;
	}
    	unmatchedPFJets++;
      }
    }
    else if( onlyCalo ){
      for(iJet = vCaloJets.begin(); iJet != vCaloJets.end(); ++iJet){ //! Calojet
	int cj = (*iJet).id; 
	if(printDebug){
	  if( kDataset == "mc" )std::cout <<" unmatched CALO jets w npf=0 : " << unmatchedCaloJets << "  pt : " << jtpt_calo[cj] << " eta : "  << jteta_calo[cj] << " phi : " << jtphi_calo[cj] << " dr : " << refdrjt_calo[cj] << std::endl;
	  else std::cout <<" unmatched CALO jets w npf=0 : " << unmatchedCaloJets << "  pt : " << jtpt_calo[cj] << " eta : "  << jteta_calo[cj] << " phi : " << jtphi_calo[cj] << std::endl;
	}
    	unmatchedCaloJets++;
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
      //! Matched jets (PF jet matched to Calo jet)
      for(itr = mCaloPFMatchedJets.cbegin(); itr != mCaloPFMatchedJets.cend(); ++itr){

	CaloPFJetPair jetpair = (*itr);
	Jet clj = jetpair.first;
	Jet pfj = jetpair.second;
	
	float delr  = deltaR(clj.eta, clj.phi, pfj.eta, pfj.phi);
	//float delpt = fabs(clj.pt - pfj.pt);
	//if( delr < kdelrmatch && caloid[clj.id]==0 && pfid[pfj.id]==0 ){ //! Removes the multimatch
	if( delr < kdelrmatch && caloid[clj.id]==0 ){//! keep the multimatch with tag 

	  if( pfid[pfj.id] > 0 ){
	    multimatchedPFJets++;
	    continue;
	  }

	  Int_t PFElecCut=0;
	  Int_t MuCut=0;
	  Int_t TrCut=0;
	  if( kSpecies == "pbpb" || kSpecies == "pbpb_mb" ){
	    Float_t Sumcand  = chargedSum[pfj.id] +  photonSum[pfj.id] + neutralSum[pfj.id] + muonSum[pfj.id];
	    Float_t caloBypf = jtpt_calo[clj.id]/jtpt[pfj.id];
	    Float_t ePFSel   = (18./7.*caloBypf) - 9./7.;

	    if( caloBypf > 0.5 && caloBypf <= 0.85 && eleMax[pfj.id]/Sumcand < ePFSel ) PFElecCut = 1;
	    if( caloBypf > 0.85 )PFElecCut = 1;
	    if( caloBypf <= 0.5 && eleMax[pfj.id]/Sumcand < 0.05) PFElecCut = 1;

	    //if( muonMax[pfj.id]/(neutralMax[pfj.id]+photonMax[pfj.id]+chargedMax[pfj.id]+muonMax[pfj.id]+eleMax[pfj.id])<0.975 )MuCut = 1;
	    //if( trackMax[pfj.id]/jtpt[pfj.id] > 0.02)TrCut=1;
	    MuCut=1;
	    TrCut=1;
	    //if( pfj.id == 53 )cout <<pfj.id<<" matched  jet pT "<< jtpt[pfj.id] <<" calo pt : "<< jtpt_calo[clj.id] << " Sumcand : "<< Sumcand <<" caloBypf  : " << caloBypf <<"  ePFSel : " << ePFSel << "  PFElecCut : " << PFElecCut << endl;
	  }else if ( kSpecies == "pp" ){
	    PFElecCut=1;
	    MuCut=1;
	    TrCut=1;
	  }


	  bool isel = true;	
	  if( kDataset=="mc"){
	    if(sid[pfj.id] != 0)isel=false;
	    if(pfrefdrjt[pfj.id] > kdelrcut)isel=false;
	    //if( fabs(pfrefpt[pfj.id]) > 2.*pthat )isel=false;
	    if( fabs(pfrefeta[pfj.id]) > ketacut )isel=false;

	    if( isel ){
	      double resp_corr =  jtpt[pfj.id]  / pfrefpt[pfj.id];
	      double resp_raw  =  rawpt[pfj.id] / pfrefpt[pfj.id];

	      int iphi     = GetPhiBin(pfrefphi[pfj.id]);
	      int ieta     = GetEtaBin(pfrefeta[pfj.id]);
	      int ipt_mc   = GetPtBin (pfrefpt[pfj.id]);
	      int ipt_wide = GetPtBinWide (pfrefpt[pfj.id]);
	      if( ipt_mc>=0 && ipt_mc<nbins ){
		hrecopt_genm  [0][iCent][ipt_mc]->Fill(jtpt[pfj.id],weight*prescl);
		hrescrpt_genm [0][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		hresrrpt_genm [0][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);

		if( iCent == ncen-1 ){ //! smearing only for pp
		  for(int ic=0; ic< ncen; ic++){
		    double smpt = GetSmearedPt(ic,jtpt[pfj.id],pfrefpt[pfj.id]);
		    double resp_smear=  smpt/pfrefpt[pfj.id];
		    hrescrpt_sm_genm [0][ic][ipt_mc]->Fill(resp_smear,weight*prescl);
		    //cout << " \t smear  " << ic << "  " << jtpt[pfj.id] << "  " << resp_corr << "  "<<  smpt <<  "   " << resp_smear << endl;
		  }
		}
		
		//! Flavor dependence
		if( abs(pfrefparton_flavor[pfj.id]) < 5){
		  hrescrpt_q_genm [0][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		  hresrrpt_q_genm [0][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);
		  
		  if( abs(pfrefparton_flavor[pfj.id]) == 1 ){
		    hrescrpt_d_genm [0][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		    hresrrpt_d_genm [0][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);
		  }else if( abs(pfrefparton_flavor[pfj.id]) == 2 ){
		    hrescrpt_u_genm [0][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		    hresrrpt_u_genm [0][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);
		  }else if( abs(pfrefparton_flavor[pfj.id]) == 3 ){
		    hrescrpt_s_genm [0][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		    hresrrpt_s_genm [0][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);
		  }
		}else if( abs(pfrefparton_flavor[pfj.id]) == 21 ){
		  hrescrpt_g_genm [0][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		  hresrrpt_g_genm [0][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);
		  //cout << " PF Matched  found a gluon jet : " << resp_corr << endl;
		}

		if( ieta>=0 && ieta<neta ){
		  hrescrpt_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
		  hresrrpt_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);

		  //! Flavor dependence
		  if( abs(pfrefparton_flavor[pfj.id]) < 5){
		    hrescrpt_q_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
		    hresrrpt_q_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		  
		    if( abs(pfrefparton_flavor[pfj.id]) == 1 ){
		      hrescrpt_d_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
		      hresrrpt_d_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		    }else if( abs(pfrefparton_flavor[pfj.id]) == 2 ){
		      hrescrpt_u_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
		      hresrrpt_u_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		    }else if( abs(pfrefparton_flavor[pfj.id]) == 3 ){
		      hrescrpt_s_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
		      hresrrpt_s_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		    }
		  }else if( abs(pfrefparton_flavor[pfj.id]) == 21 ){
		    hrescrpt_g_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
		    hresrrpt_g_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		  }


		  if( ipt_wide>=0 && ipt_wide <nbins_wide){
		    hrescrpt_wide_genm_eta [0][iCent][ipt_wide][ieta]->Fill(resp_corr,weight*prescl);
		    hresrrpt_wide_genm_eta [0][iCent][ipt_wide][ieta]->Fill(resp_raw ,weight*prescl);
		  }
		}
		if( iphi>=0 && iphi<nphi ){
		  hrescrpt_genm_phi [0][iCent][ipt_mc][iphi]->Fill(resp_corr,weight*prescl);
		  hresrrpt_genm_phi [0][iCent][ipt_mc][iphi]->Fill(resp_raw ,weight*prescl);
		  if( ipt_wide>=0 && ipt_wide<nbins_wide){
		    hrescrpt_wide_genm_phi [0][iCent][ipt_wide][iphi]->Fill(resp_corr,weight*prescl);
		    hresrrpt_wide_genm_phi [0][iCent][ipt_wide][iphi]->Fill(resp_raw ,weight*prescl);
		  }
		}
	      }
	      if( (PFElecCut*MuCut*TrCut) == 1 ){
		if( ipt_mc>=0 && ipt_mc<nbins ){
		  hrecopt_genm  [1][iCent][ipt_mc]->Fill(jtpt[pfj.id],weight*prescl);
		  hrescrpt_genm [1][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		  hresrrpt_genm [1][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);

		  if( iCent == ncen-1 ){ //! smearing only for pp
		    for(int ic=0; ic< ncen; ic++){
		      double smpt = GetSmearedPt(ic,jtpt[pfj.id],pfrefpt[pfj.id]);
		      double resp_smear=  smpt/pfrefpt[pfj.id];
		      hrescrpt_sm_genm [1][ic][ipt_mc]->Fill(resp_smear,weight*prescl);
		    }
		  }

		  //! Flavor dependence
		  if( abs(pfrefparton_flavor[pfj.id]) < 5){
		    hrescrpt_q_genm [1][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		    hresrrpt_q_genm [1][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);
		  
		    if( abs(pfrefparton_flavor[pfj.id]) == 1 ){
		      hrescrpt_d_genm [1][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		      hresrrpt_d_genm [1][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);
		    }else if( abs(pfrefparton_flavor[pfj.id]) == 2 ){
		      hrescrpt_u_genm [1][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		      hresrrpt_u_genm [1][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);
		    }else if( abs(pfrefparton_flavor[pfj.id]) == 3 ){
		      hrescrpt_s_genm [1][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		      hresrrpt_s_genm [1][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);
		    }
		  }else if( abs(pfrefparton_flavor[pfj.id]) == 21 ){
		    hrescrpt_g_genm [1][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		    hresrrpt_g_genm [1][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);
		    //cout << " PF Matched PFElec  found a gluon jet : " << resp_corr << endl;
		  }


		  if( ieta>=0 && ieta<neta ){
		    hrescrpt_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
		    hresrrpt_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);

		    //! Flavor dependence
		    if( abs(pfrefparton_flavor[pfj.id]) < 5){
		      hrescrpt_q_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
		      hresrrpt_q_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		  
		      if( abs(pfrefparton_flavor[pfj.id]) == 1 ){
			hrescrpt_d_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
			hresrrpt_d_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		      }else if( abs(pfrefparton_flavor[pfj.id]) == 2 ){
			hrescrpt_u_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
			hresrrpt_u_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		      }else if( abs(pfrefparton_flavor[pfj.id]) == 3 ){
			hrescrpt_s_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
			hresrrpt_s_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		      }
		    }else if( abs(pfrefparton_flavor[pfj.id]) == 21 ){
		      hrescrpt_g_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
		      hresrrpt_g_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		    }

		    if( ipt_wide>=0 && ipt_wide <nbins_wide){
		      hrescrpt_wide_genm_eta [1][iCent][ipt_wide][ieta]->Fill(resp_corr,weight*prescl);
		      hresrrpt_wide_genm_eta [1][iCent][ipt_wide][ieta]->Fill(resp_raw ,weight*prescl);
		    }
		  }
		  if( iphi>=0 && iphi<nphi ){
		    hrescrpt_genm_phi [1][iCent][ipt_mc][iphi]->Fill(resp_corr,weight*prescl);
		    hresrrpt_genm_phi [1][iCent][ipt_mc][iphi]->Fill(resp_raw ,weight*prescl);
		    if( ipt_wide>=0 && ipt_wide<nbins_wide){
		      hrescrpt_wide_genm_phi [1][iCent][ipt_wide][iphi]->Fill(resp_corr,weight*prescl);
		      hresrrpt_wide_genm_phi [1][iCent][ipt_wide][iphi]->Fill(resp_raw ,weight*prescl);
		    }
		  }
		}
	      }
	    }
	  }
	  
	  float sumPF = (chargedSum[pfj.id] + neutralSum[pfj.id] + photonSum[pfj.id] + muonSum[pfj.id] + eleSum[pfj.id]);
	  float bkgd  = sumPF - rawpt[pfj.id];
	  
	  int ipt = GetPtBin(jtpt[pfj.id]);
	  int ieta_wide = GetEtaBinWide(fabs(jteta[pfj.id]));

	  float trig_deta = 0;
	  float trig_dphi = 0;
	  float trig_dpt  = 0;
	  if( kDataset == "data" ){
	    trig_deta = trgObj_eta - jteta[pfj.id];
	    trig_dphi = delPhi(trgObj_phi,jtphi[pfj.id]);
	    trig_dpt  = trgObj_pt  - jtpt[pfj.id];
	  }
	  
	  if( TrigFill && isel ){
	    if(ipt >=0 && ipt < nbins){
	      if(ieta_wide >= 0 && ieta_wide < neta_wide){

		htrig_deta[iCent][iTrig][ipt]->Fill(trig_deta);
		htrig_dphi[iCent][iTrig][ipt]->Fill(trig_dphi);
		htrig_dpt [iCent][iTrig][ipt]->Fill(trig_dpt);

		havbkgpt_comb    [0][iCent][ipt]->Fill(bkgd, weight*prescl);
		havbkgpt_comb_eta[0][iCent][ipt][ieta_wide]->Fill(bkgd, weight*prescl);
		havbkgpt_ind     [0][iCent][iTrig][ipt]->Fill(bkgd, weight*prescl);
		havbkgpt_ind_eta [0][iCent][iTrig][ipt][ieta_wide]->Fill(bkgd, weight*prescl);
		if( (PFElecCut*MuCut*TrCut) == 1 ){
		  havbkgpt_comb    [1][iCent][ipt]->Fill(bkgd, weight*prescl);
		  havbkgpt_comb_eta[1][iCent][ipt][ieta_wide]->Fill(bkgd, weight*prescl);
		  havbkgpt_ind     [1][iCent][iTrig][ipt]->Fill(bkgd, weight*prescl);
		  havbkgpt_ind_eta [1][iCent][iTrig][ipt][ieta_wide]->Fill(bkgd, weight*prescl);
		}
	      }
	    }
	  }
	  
	  if(printDebug){
	    if( kDataset == "data" )std::cout <<"\t *** Matched" << " id : " << pfj.id << " " << clj.id << " pfpt :  " << pfj.pt << " calopt : " << clj.pt << std::endl;
	    else {
	      std::cout <<"\t *** Matched" << " id : " << pfj.id << " " << clj.id << " pfpt :  " << pfj.pt << " calopt : " << clj.pt << " spf : " << sid[pfj.id] << " scj : " << sid [clj.id] << " dr : " << pfrefdrjt[pfj.id] << std::endl;
	    }
	  }
	  //if( sid[pfj.id]!=0 && pfpt>20 )std::cout <<"\t *** Matched" << " pthat : " << pthat << " refpt : " << refpt << " pfpt :  " << pfj.pt << " calopt : " << clj.pt << " spf : " << sid[pfj.id] << " scj : " << sid [clj.id] << std::endl;	  
	  matchedJets++;	  
	  pfid  [pfj.id] = 1;
	  caloid[clj.id] = 1;
	}
      }//! itr loop
      

      //! Unmatched jets
      for(itr = mCaloPFMatchedJets.cbegin(); itr != mCaloPFMatchedJets.cend(); ++itr){
	CaloPFJetPair jetpair = (*itr);
	
	Jet clj = jetpair.first;
	Jet pfj = jetpair.second;
	
	//float delr  = deltaR(pfj.eta, pfj.phi, clj.eta, clj.phi);
	//float delpt = fabs(pfj.pt - clj.pt);
	//if( pfid[pfj.id]==1 || caloid[clj.id]==1 )continue;
	
	if( pfid[pfj.id] == 0 ){//! Unmatched PF jet
	  
	  Int_t PFElecCut = 0;
	  Int_t MuCut=0;
	  Int_t TrCut=0;
	  Float_t Sumcand = chargedSum[pfj.id] +  photonSum[pfj.id] + neutralSum[pfj.id] + muonSum[pfj.id];
	  if( kSpecies == "pbpb" || kSpecies == "pbpb_mb" ){
	    if( eleMax[pfj.id]/Sumcand < 0.05 )PFElecCut=1;      
	    //if( muonMax[pfj.id]/(neutralMax[pfj.id]+photonMax[pfj.id]+chargedMax[pfj.id]+muonMax[pfj.id]+eleMax[pfj.id])<0.975 )MuCut = 1;
	    //if( trackMax[pfj.id]/jtpt[pfj.id] > 0.02)TrCut=1;
	    MuCut=1;
	    TrCut=1;
	  }else if( kSpecies == "pp" ){
	    PFElecCut=1;
	    MuCut=1;
	    TrCut=1;
	  }

	  bool isel = true;	
	  if( kDataset=="mc"){
	    if(sid[pfj.id] != 0)isel=false;
	    if(pfrefdrjt[pfj.id] > kdelrcut)isel=false;
	    //if( fabs(pfrefpt[pfj.id]) > 2.*pthat )isel=false;
	    if( fabs(pfrefpt[pfj.id]) > ketacut )isel=false;

	    if( isel ){
	      double resp_corr =  jtpt[pfj.id]  / pfrefpt[pfj.id];
	      double resp_raw  =  rawpt[pfj.id] / pfrefpt[pfj.id];
	      int iphi     = GetPhiBin(pfrefphi[pfj.id]);
	      int ieta     = GetEtaBin(pfrefeta[pfj.id]);
	      int ipt_mc   = GetPtBin (pfrefpt[pfj.id]);
	      int ipt_wide = GetPtBinWide (pfrefpt[pfj.id]);
	      if( ipt_mc>=0 && ipt_mc<nbins ){
		hrecopt_genm  [0][iCent][ipt_mc]->Fill(jtpt[pfj.id],weight*prescl);
		hrescrpt_genm [0][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		hresrrpt_genm [0][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);

		if( iCent == ncen-1 ){ //! smearing only for pp
		  for(int ic=0; ic< ncen; ic++){
		    double smpt = GetSmearedPt(ic,jtpt[pfj.id], pfrefpt[pfj.id]);
		    double resp_smear=  smpt/pfrefpt[pfj.id];
		    hrescrpt_sm_genm [0][ic][ipt_mc]->Fill(resp_smear,weight*prescl);
		  }
		}


		//! Flavor dependence
		if( abs(pfrefparton_flavor[pfj.id]) < 5){
		  hrescrpt_q_genm [0][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		  hresrrpt_q_genm [0][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);
		  
		  if( abs(pfrefparton_flavor[pfj.id]) == 1 ){
		    hrescrpt_d_genm [0][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		    hresrrpt_d_genm [0][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);
		  }else if( abs(pfrefparton_flavor[pfj.id]) == 2 ){
		    hrescrpt_u_genm [0][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		    hresrrpt_u_genm [0][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);
		  }else if( abs(pfrefparton_flavor[pfj.id]) == 3 ){
		    hrescrpt_s_genm [0][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		    hresrrpt_s_genm [0][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);
		  }
		}else if( abs(pfrefparton_flavor[pfj.id]) == 21 ){
		  hrescrpt_g_genm [0][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		  hresrrpt_g_genm [0][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);
		}


		if( ieta>=0 && ieta<neta ){
		  hrescrpt_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
		  hresrrpt_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);

		  //! Flavor dependence
		  if( abs(pfrefparton_flavor[pfj.id]) < 5){
		    hrescrpt_q_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
		    hresrrpt_q_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		    
		    if( abs(pfrefparton_flavor[pfj.id]) == 1 ){
		      hrescrpt_d_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
		      hresrrpt_d_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		    }else if( abs(pfrefparton_flavor[pfj.id]) == 2 ){
		      hrescrpt_u_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
		      hresrrpt_u_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		    }else if( abs(pfrefparton_flavor[pfj.id]) == 3 ){
		      hrescrpt_s_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
		      hresrrpt_s_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		    }
		  }else if( abs(pfrefparton_flavor[pfj.id]) == 21 ){
		    hrescrpt_g_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
		    hresrrpt_g_genm_eta [0][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		  }

		  if( ipt_wide>=0 && ipt_wide <nbins_wide){
		    hrescrpt_wide_genm_eta [0][iCent][ipt_wide][ieta]->Fill(resp_corr,weight*prescl);
		    hresrrpt_wide_genm_eta [0][iCent][ipt_wide][ieta]->Fill(resp_raw ,weight*prescl);
		  }
		}
		if( iphi>=0 && iphi<nphi ){
		  hrescrpt_genm_phi [0][iCent][ipt_mc][iphi]->Fill(resp_corr,weight*prescl);
		  hresrrpt_genm_phi [0][iCent][ipt_mc][iphi]->Fill(resp_raw ,weight*prescl);
		  if( ipt_wide>=0 && ipt_wide<nbins_wide){
		    hrescrpt_wide_genm_phi [0][iCent][ipt_wide][iphi]->Fill(resp_corr,weight*prescl);
		    hresrrpt_wide_genm_phi [0][iCent][ipt_wide][iphi]->Fill(resp_raw ,weight*prescl);
		  }
		}
	      }
	      if( (PFElecCut*MuCut*TrCut) == 1 ){
		if( ipt_mc>=0 && ipt_mc<nbins ){
		  hrecopt_genm  [1][iCent][ipt_mc]->Fill(jtpt[pfj.id],weight*prescl);
		  hrescrpt_genm [1][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		  hresrrpt_genm [1][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);

		  if( iCent == ncen-1 ){ //! smearing only for pp
		    for(int ic=0; ic< ncen; ic++){
		      double smpt = GetSmearedPt(ic,jtpt[pfj.id], pfrefpt[pfj.id]);
		      double resp_smear=  smpt/pfrefpt[pfj.id];
		      hrescrpt_sm_genm [1][ic][ipt_mc]->Fill(resp_smear,weight*prescl);
		    }
		  }

		  //! Flavor dependence
		  if( abs(pfrefparton_flavor[pfj.id]) < 5){
		    hrescrpt_q_genm [1][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		    hresrrpt_q_genm [1][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);
		  
		    if( abs(pfrefparton_flavor[pfj.id]) == 1 ){
		      hrescrpt_d_genm [1][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		      hresrrpt_d_genm [1][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);
		    }else if( abs(pfrefparton_flavor[pfj.id]) == 2 ){
		      hrescrpt_u_genm [1][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		      hresrrpt_u_genm [1][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);
		    }else if( abs(pfrefparton_flavor[pfj.id]) == 3 ){
		      hrescrpt_s_genm [1][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		      hresrrpt_s_genm [1][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);
		    }
		  }else if( abs(pfrefparton_flavor[pfj.id]) == 21 ){
		    hrescrpt_g_genm [1][iCent][ipt_mc]->Fill(resp_corr,weight*prescl);
		    hresrrpt_g_genm [1][iCent][ipt_mc]->Fill(resp_raw ,weight*prescl);
		  }

		  if( ieta>=0 && ieta<neta ){
		    hrescrpt_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
		    hresrrpt_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		    //! Flavor dependence
		    if( abs(pfrefparton_flavor[pfj.id]) < 5){
		      hrescrpt_q_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
		      hresrrpt_q_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		      
		      if( abs(pfrefparton_flavor[pfj.id]) == 1 ){
			hrescrpt_d_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
			hresrrpt_d_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		      }else if( abs(pfrefparton_flavor[pfj.id]) == 2 ){
			hrescrpt_u_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
			hresrrpt_u_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		      }else if( abs(pfrefparton_flavor[pfj.id]) == 3 ){
			hrescrpt_s_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
			hresrrpt_s_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		      }
		    }else if( abs(pfrefparton_flavor[pfj.id]) == 21 ){
		      hrescrpt_g_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_corr,weight*prescl);
		      hresrrpt_g_genm_eta [1][iCent][ipt_mc][ieta]->Fill(resp_raw ,weight*prescl);
		    }

		    if( ipt_wide>=0 && ipt_wide <nbins_wide){
		      hrescrpt_wide_genm_eta [1][iCent][ipt_wide][ieta]->Fill(resp_corr,weight*prescl);
		      hresrrpt_wide_genm_eta [1][iCent][ipt_wide][ieta]->Fill(resp_raw ,weight*prescl);
		    }
		  }
		  if( iphi>=0 && iphi<nphi ){
		    hrescrpt_genm_phi [1][iCent][ipt_mc][iphi]->Fill(resp_corr,weight*prescl);
		    hresrrpt_genm_phi [1][iCent][ipt_mc][iphi]->Fill(resp_raw ,weight*prescl);
		    if( ipt_wide>=0 && ipt_wide<nbins_wide){
		      hrescrpt_wide_genm_phi [1][iCent][ipt_wide][iphi]->Fill(resp_corr,weight*prescl);
		      hresrrpt_wide_genm_phi [1][iCent][ipt_wide][iphi]->Fill(resp_raw ,weight*prescl);
		    }
		  }
		}
	      }
	    }
	  }
	  
	  float sumPF = (chargedSum[pfj.id] + neutralSum[pfj.id] + photonSum[pfj.id] + muonSum[pfj.id] + eleSum[pfj.id]);
	  float bkgd  = sumPF - rawpt[pfj.id];
	  
	  int ipt = GetPtBin(jtpt[pfj.id]);
	  int ieta_wide = GetEtaBinWide(fabs(jteta[pfj.id]));
	  
	  float trig_deta = 0;
	  float trig_dphi = 0;
	  float trig_dpt  = 0;
	  if( kDataset == "data" ){
	    trig_deta = trgObj_eta - jteta[pfj.id];
	    trig_dphi = delPhi(trgObj_phi,jtphi[pfj.id]);
	    trig_dpt  = trgObj_pt  - jtpt[pfj.id];
	  }
	  
	  if( TrigFill && isel ){
	    if(ipt >=0 && ipt < nbins){
	      if(ieta_wide >= 0 && ieta_wide < neta_wide){

		htrig_deta[iCent][iTrig][ipt]->Fill(trig_deta);
		htrig_dphi[iCent][iTrig][ipt]->Fill(trig_dphi);
		htrig_dpt [iCent][iTrig][ipt]->Fill(trig_dpt);

		havbkgpt_comb    [0][iCent][ipt]->Fill(bkgd, weight*prescl);
		havbkgpt_comb_eta[0][iCent][ipt][ieta_wide]->Fill(bkgd, weight*prescl);
		havbkgpt_ind     [0][iCent][iTrig][ipt]->Fill(bkgd, weight*prescl);
		havbkgpt_ind_eta [0][iCent][iTrig][ipt][ieta_wide]->Fill(bkgd, weight*prescl);
		if( (PFElecCut*MuCut*TrCut) == 1 ){
		  havbkgpt_comb    [1][iCent][ipt]->Fill(bkgd, weight*prescl);
		  havbkgpt_comb_eta[1][iCent][ipt][ieta_wide]->Fill(bkgd, weight*prescl);
		  havbkgpt_ind     [1][iCent][iTrig][ipt]->Fill(bkgd, weight*prescl);
		  havbkgpt_ind_eta [1][iCent][iTrig][ipt][ieta_wide]->Fill(bkgd, weight*prescl);
		}
	      }
	    }
	  }
	  
	  if(printDebug){
	    if( kDataset == "data" )std::cout <<"\t %%% UnMatched PF" << " id  : " << pfj.id << " pfpt :  " << pfj.pt << std::endl;
	    else {
	      std::cout <<"\t %%% UnMatched PF" << " id  : " << pfj.id << " pfpt :  " << pfj.pt << " subid : " << sid[pfj.id] << " dr : " << pfrefdrjt[pfj.id] << std::endl;
	    }
	  }
	  unmatchedPFJets++;	  	  
	  pfid  [pfj.id] = 1;	  
	}
	
	if( caloid[clj.id] == 0 ){//! Unmatched Calo jet
	  if(printDebug){
	    if( kDataset == "mc" )std::cout <<"\t XXX UnMatched Calo" << " id : " << clj.id  << " calopt : " << clj.pt << " subid : " << sid_calo [clj.id] << " dr : " << refdrjt_calo[clj.id] << std::endl;
	    else std::cout <<"\t XXX UnMatched Calo" << " id : " << clj.id  << " calopt : " << clj.pt << std::endl;
	  }
	  unmatchedCaloJets++;	  	  
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
		<<" All ==> PF : " << nref_pf <<", CALO : "<< nref_calo << ";"  
		<<" After cut ==> PF : " << pfjets << ", Calo : "  << calojets << ";"
		<<" mCALOPF : "<< matchedJets <<", multiMatched : "<< multimatchedPFJets << ", unmPF : "<<  unmatchedPFJets <<", unmCalo : "<<  unmatchedCaloJets 
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
int GetEtaBinWide(float eta)
{
  for(int ix=0;ix<neta_wide;ix++){
    if(eta>=etabins_wide[ix] && eta<etabins_wide[ix+1]){
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
int GetPtBinWide(float pt)
{
  for(int ix=0;ix<nbins_wide;ix++){
    if(pt>=ptbins_wide[ix] && pt<ptbins_wide[ix+1]){
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
  // if(bin<=10)ibin=0; //! 0-5%
  // else if(bin>10  && bin<=20 )ibin=1; //! 5-10%
  // else if(bin>20  && bin<=60 )ibin=2;  //! 10-30%
  // else if(bin>60  && bin<=100)ibin=3;  //! 30-50%
  // else if(bin>100 && bin<=140)ibin=4;  //! 50-70%
  // else if(bin>140 && bin<=180)ibin=5;  //! 70-90%
  // else if(bin>180 && bin<=200)ibin=6;  //! 90-100%


  if(bin>=0 && bin<10)ibin=0; //! 0-5%
  else if(bin>=10  && bin<20 )ibin=1; //! 5-10%
  else if(bin>=20  && bin<60 )ibin=2;  //! 10-30%
  else if(bin>=60  && bin<100)ibin=3;  //! 30-50%
  else if(bin>=100 && bin<140)ibin=4;  //! 50-70%
  else if(bin>=140 && bin<180)ibin=5;  //! 70-90%
  else if(bin>=180 && bin<200)ibin=6;  //! 90-100%

  return ibin;
}
float delPhi(float phi1, float phi2)
{
  float dphi = phi2 - phi1;
  //dphi = fabs(atan2(sin(dphi),cos(dphi)));
  if(dphi < -1.0) dphi+=2*pi;
  if(dphi > pi2 ) dphi-=2*pi;
  return dphi;
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

double GetXsecWt(double pthat, string ksp)
{
  double wt=1.0;
  if( ksp == "pp" ){
    for( int i=0; i<11; i++){
      if( pthat >=  pthatwt_pp[i][0] )wt = pthatwt_pp[i][3];
    }
  }else if( ksp == "pbpb" ){
    for( int i=0; i<9; i++){
      if( pthat >=  pthatwt_pbpb[i][0] )wt = pthatwt_pbpb[i][3];
    }
  }else if( ksp == "pbpb_mb" ){
    wt = 1.0;
  }
  return wt;
}


void GetCentWeight(TH1D *hCentWeight)
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
  //  cout<<endl;
}
double GetSmearedPt(int icent, double recopt, double refpt)
{
  // if( recopt < 60. || recopt > 300. )return recopt;

  switch(icent){
  case 0: //! 0-5%
    fSmear->SetParameters(0.378558,-0.0045491,2.85944e-05,-8.58614e-08);
    break;
  case 1: //! 5-10%
    fSmear->SetParameters(0.285197,-0.00255444,1.11539e-05,-1.71385e-08);
    break;
  case 2: //! 10-30%
    fSmear->SetParameters(0.344794,-0.00443056,2.92684e-05,-9.1955e-08);
    break;
  case 3: //! 30-50%
    fSmear->SetParameters(0.200427,-0.00139927,3.16825e-06,6.15753e-09);
    break;
  case 4: //! 50-70%
    fSmear->SetParameters(0.134483,3.4815e-05,-1.02942e-05,6.06084e-08);
    break;
  case 5: //!70-90%
    fSmear->SetParameters(0.174897,-0.00175978,1.05952e-05,-3.18472e-08);
    break;
  case 6: //! 90-100%
    fSmear->SetParameters(0.190651,-0.00234196,1.5059e-05,-4.26833e-08);
    break;
  case 7: //! pp
    fSmear->SetParameters(0,0,0,0);    
  default : //! No smear
    fSmear->SetParameters(0,0,0,0);
    break;
  }

  double smpt = gRandom->Gaus(recopt,fSmear->Eval(refpt));
  return smpt;
}
