
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
double delphi(double /*phi1*/, double /*phi2*/);
double deltaR(float /*eta1*/, float /*phi1*/, 
	      float /*eta2*/, float /*phi2*/);
struct Jet{
  float pt;
  float eta;
  float phi;
  int id;
};
bool compare_pt(Jet jet1, Jet jet2);
bool compare_pt(Jet jet1, Jet jet2){
  return jet1.pt > jet2.pt;
}

//! Triggers  55 || 65 || 80
//!  no difference between 55_p_1, 65_p_1 and 80_p_1 all are exactly same
const int ntrig=18;
const char *ctrig[ntrig] = {
  " jet80_1",
  " jet80_p_1",
  " jet65_1",
  " jet65_p_1",
  " jet55_1",
  " jet55_p_1",

  " jet80_1 && !jet80_p_1",
  " jet65_1 && !jet65_p_1",
  " jet55_1 && !jet55_p_1",
  " jet55_1 &&  jet65_1",
  " jet65_1 &&  jet80_1",
  " jet55_1 &&  jet80_1",
  " jet65_1 && !jet80_1",
  "!jet55_1 &&  jet65_1",
  " jet55_1 && !jet65_1"	  
};


const int ncen=6;
const char *cdir[ncen]  = {"05","510","1030","3050","5070","70100"};
const char *ccent[ncen] = {"0-5%","5-10%","10-30%","30-50%","50-70%","70-90%"};

const int ncand=5;
const char *ccand[ncand] = {"h^{#pm}","e^{#pm}","#mu^{#pm}","#gamma","h0"};

//! constants
int iYear=2015;
const double pi=acos(-1.);

const double ketacut=2.0;
const double kptrawcut =0.0;
const double kptrecocut=20.0;

//! pt binning
//double ptbins[] ={24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 
// 		  114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 
// 		  395, 430, 468, 507, 548, 592, 638, 686, 1000};
double ptbins[] ={30,40,50,60,70,80,90,100,110,120,130,140,160,200,250,300,400,548};
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
int CondorPbPb_Data(std::string kAlgName="akPu3",
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
  TChain *tch_trgobj = new TChain("hltobject/jetObjTree");
  AddInputFiles(tch_trgobj,fileList,"hltobject/jetObjTree");
  cout <<" # of events in TrigObj Tree : " <<  tch_trgobj->GetEntries() <<endl;
  cout <<endl;


  //! Event Tree
  float vz;
  int hiBin;
  int hiNpix;
  tch_ev->SetBranchAddress("vz",&vz);
  tch_ev->SetBranchAddress("hiBin",&hiBin);
  tch_ev->SetBranchAddress("hiNpix",&hiNpix);
  
  //! HLT tree
  int jet55_1;
  int jet55_p_1;
  int jet65_1;
  int jet65_p_1;
  int jet80_1;
  int jet80_p_1;
  int L1_sj36_1;
  int L1_sj36_p_1;
  int L1_sj52_1;
  int L1_sj52_p_1;
  tch_hlt->SetBranchAddress("HLT_HIJet55_v1",&jet55_1);
  tch_hlt->SetBranchAddress("HLT_HIJet55_v1_Prescl",&jet55_p_1);
  tch_hlt->SetBranchAddress("HLT_HIJet65_v1",&jet65_1);
  tch_hlt->SetBranchAddress("HLT_HIJet65_v1_Prescl",&jet65_p_1);
  tch_hlt->SetBranchAddress("HLT_HIJet80_v1",&jet80_1);
  tch_hlt->SetBranchAddress("HLT_HIJet80_v1_Prescl",&jet80_p_1);
  tch_hlt->SetBranchAddress("L1_SingleJet36_BptxAND",&L1_sj36_1);
  tch_hlt->SetBranchAddress("L1_SingleJet36_BptxAND_Prescl",&L1_sj36_p_1);
  tch_hlt->SetBranchAddress("L1_SingleJet52_BptxAND",&L1_sj52_1);
  tch_hlt->SetBranchAddress("L1_SingleJet52_BptxAND_Prescl",&L1_sj52_p_1);

  //! Skim Tree
  int pcollisionEventSelection;
  int pHBHENoiseFilter;
  tch_skim->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
  tch_skim->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);

  //! Trigger object tree
  float trgObj_id;
  float trgObj_pt;
  float trgObj_eta;
  float trgObj_phi;
  tch_trgobj->SetBranchAddress("id",&trgObj_id);
  tch_trgobj->SetBranchAddress("pt",&trgObj_pt);
  tch_trgobj->SetBranchAddress("eta",&trgObj_eta);
  tch_trgobj->SetBranchAddress("phi",&trgObj_phi);


  //! CaloJet tree
  //Declaration of leaves types
  int   nref_calo;
  float jtpt_calo[1000];
  float rawpt_calo[1000];
  float jteta_calo[1000];
  float jtphi_calo[1000];

  tch_calojet->SetBranchAddress("nref",&nref_calo);
  tch_calojet->SetBranchAddress("rawpt",rawpt_calo);
  tch_calojet->SetBranchAddress("jtpt" ,jtpt_calo);
  tch_calojet->SetBranchAddress("jteta",jteta_calo);
  tch_calojet->SetBranchAddress("jtphi",jtphi_calo);


  //! PFJet tree
  //Declaration of leaves types
  int   nref;
  float jtpt[1000];
  float rawpt[1000];
  float jteta[1000];
  float jtphi[1000];
  float neutralSum[1000];
  float chargedSum[1000];
  float photonSum[1000];
  float eSum[1000];
  float muSum[1000];
  float neutralMax[1000];
  float chargedMax[1000];
  float photonMax[1000];
  float eMax[1000];
  float muMax[1000];

  tch_pfjet->SetBranchAddress("nref",&nref);
  tch_pfjet->SetBranchAddress("rawpt",rawpt);
  tch_pfjet->SetBranchAddress("jtpt" ,jtpt);
  tch_pfjet->SetBranchAddress("jteta",jteta);
  tch_pfjet->SetBranchAddress("jtphi",jtphi);
  tch_pfjet->SetBranchAddress("neutralSum",neutralSum);
  tch_pfjet->SetBranchAddress("chargedSum",chargedSum);
  tch_pfjet->SetBranchAddress("photonSum",photonSum);
  tch_pfjet->SetBranchAddress("eSum",eSum);
  tch_pfjet->SetBranchAddress("muSum",muSum);
  tch_pfjet->SetBranchAddress("neutralMax",neutralMax);
  tch_pfjet->SetBranchAddress("chargedMax",chargedMax);
  tch_pfjet->SetBranchAddress("photonMax",photonMax);
  tch_pfjet->SetBranchAddress("muMax",muMax);
  tch_pfjet->SetBranchAddress("eMax",eMax);
  
  tch_pfjet->AddFriend(tch_ev);
  tch_pfjet->AddFriend(tch_hlt);
  tch_pfjet->AddFriend(tch_skim);
  tch_pfjet->AddFriend(tch_trgobj);
  tch_pfjet->AddFriend(tch_calojet);

  //! Disable branches 
  //! Jet Tree
  tch_pfjet->SetBranchStatus("*",0,0);
  tch_pfjet->SetBranchStatus("nref" ,1,0);
  tch_pfjet->SetBranchStatus("rawpt",1,0);
  tch_pfjet->SetBranchStatus("jtpt" ,1,0);
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
  tch_pfjet->SetBranchStatus("muMax",1,0);
  tch_pfjet->SetBranchStatus("eMax",1,0);

  tch_pfjet->SetBranchStatus("vz",1,0);
  tch_pfjet->SetBranchStatus("hiBin",1,0);
  tch_pfjet->SetBranchStatus("hiNpix",1,0);

  tch_pfjet->SetBranchStatus("id",1,0);
  tch_pfjet->SetBranchStatus("pt",1,0);
  tch_pfjet->SetBranchStatus("eta",1,0);
  tch_pfjet->SetBranchStatus("phi",1,0);

  tch_pfjet->SetBranchStatus("HLT_HIJet55_v1",1,0);
  tch_pfjet->SetBranchStatus("HLT_HIJet55_v1_Prescl",1,0);
  tch_pfjet->SetBranchStatus("HLT_HIJet65_v1",1,0);
  tch_pfjet->SetBranchStatus("HLT_HIJet65_v1_Prescl",1,0);
  tch_pfjet->SetBranchStatus("HLT_HIJet80_v1",1,0);
  tch_pfjet->SetBranchStatus("HLT_HIJet80_v1_Prescl",1,0);
  tch_pfjet->SetBranchStatus("L1_SingleJet36_BptxAND",1,0);
  tch_pfjet->SetBranchStatus("L1_SingleJet36_BptxAND_Prescl",1,0);
  tch_pfjet->SetBranchStatus("L1_SingleJet52_BptxAND",1,0);
  tch_pfjet->SetBranchStatus("L1_SingleJet52_BptxAND_Prescl",1,0);

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

  TH1F *hBin_full  = new TH1F("hBin_full","Centrality bin",200,-0.5,200-0.5);
  hBin_full->Sumw2();

  TH1F *hVz_full  = new TH1F("hVz_full","Vz (in cm)",80,-20.0,20.0);
  hVz_full->Sumw2();

  TH2F *hBin_NJets = new TH2F("hBin_NJets","hiBin vs. # of pfjets / #of calojets",200,-0.5,200-0.5,100,0.,5.0);
  hBin_NJets->Sumw2();

  TH1F *hEvts[ntrig];
  TH1F *hVz[ntrig];
  TH1F *hBin[ntrig];

  TProfile *hchargedfrc[ntrig][ncen], *hphotonfrc[ntrig][ncen], *hneutralfrc[ntrig][ncen], *hmuonfrc[ntrig][ncen], *helecfrc[ntrig][ncen];
  TProfile *hchargedsum[ntrig][ncen][2], *hphotonsum[ntrig][ncen][2], *hneutralsum[ntrig][ncen][2], *hmuonsum[ntrig][ncen][2], *helecsum[ntrig][ncen][2];
  TProfile *hchargedmax[ntrig][ncen][2], *hphotonmax[ntrig][ncen][2], *hneutralmax[ntrig][ncen][2], *hmuonmax[ntrig][ncen][2], *helecmax[ntrig][ncen][2];

  //! Gen matched jets
  TH2F *hrawrecpt[ntrig][ncen][2];
  TH1F *hrecopt [ntrig][ncen][2], *hrawpt[ntrig][ncen][2], *hjeteta[ntrig][ncen][2], *hjetphi[ntrig][ncen][2];
  TH2F *havbkgpt[ntrig][ncen][2], *havbkgpt_eta[ntrig][ncen][2][neta], *havbkgpt_phi[ntrig][ncen][2][nphi];


  //! Data-driven JEC's
  TH2F *hrescrpt_datadr [2][ncen];
  for(int i=0; i<2; i++){	
    for(int ic=0; ic<ncen; ic++){
      hrescrpt_datadr [i][ic] = new TH2F(Form("hrescrpt_datadr_%d_%d",i,ic),Form("jet p_{T} %d %s %s ",i,kAlgName.c_str(),ccent[ic]), 30, 50, 650, 50,-1.50,1.50);
      hrescrpt_datadr [i][ic]->Sumw2();
    }
  }
  
  for(int it=0; it<ntrig; it++){
    hEvts[it]  = new TH1F(Form("hEvts%d",it),Form("# of events, %s",ctrig[it]),ntrig,-0.5,ntrig-0.5);
    hEvts[it]->Sumw2();

    hVz[it]  = new TH1F(Form("hVz%d",it),Form("Vz (in cm), %s",ctrig[it]),80,-20.,20.);
    hVz[it]->Sumw2();

    hBin[it]  = new TH1F(Form("hBin%d",it),Form("Centrality bin, %s",ctrig[it]),200,-0.5,200-0.5);
    hBin[it]->Sumw2();


    for(int ic=0; ic<ncen; ic++){
      
      hchargedfrc[it][ic] = new TProfile(Form("hchargedfrc%d_%d",it,ic),"chargedMax/recopt",nbins,ptbins);
      hchargedfrc[it][ic]->Sumw2();
      hphotonfrc [it][ic] = new TProfile(Form("hphotonfrc%d_%d",it,ic),"photonMax/recopt",nbins,ptbins);
      hphotonfrc [it][ic]->Sumw2();
      hneutralfrc[it][ic] = new TProfile(Form("hneutralfrc%d_%d",it,ic),"neutralMax/recopt",nbins,ptbins);
      hneutralfrc[it][ic]->Sumw2();
      hmuonfrc[it][ic] = new TProfile(Form("hmuonfrc%d_%d",it,ic),"muonMax/recopt",nbins,ptbins);
      hmuonfrc[it][ic]->Sumw2();
      helecfrc[it][ic] = new TProfile(Form("helecfrc%d_%d",it,ic),"elecMax/recopt",nbins,ptbins);
      helecfrc[it][ic]->Sumw2();


      
      for(int i=0; i<2; i++){	
	
	hrawrecpt[it][ic][i] = new TH2F(Form("hrecrawpt%d_%d_%d",it,ic,i),Form("Rec-Raw matched gen p_{T} distribution jet centb %s %s %s",ctrig[it],kAlgName.c_str(),ccent[ic]),200,0,1000,200,0,1000);
	hrawrecpt[it][ic][i]->Sumw2();

	hrecopt[it][ic][i] = new TH1F(Form("hrecopt%d_%d_%d",it,ic,i),Form("Reco p_{T} distribution jet %s %s %s",ctrig[it],kAlgName.c_str(),ccent[ic]),200,0,1000);
	hrecopt[it][ic][i]->Sumw2();
	hrawpt [it][ic][i] = new TH1F(Form("hrawpt%d_%d_%d",it,ic,i),Form("Raw p_{T} distribution jet  %s %s %s",ctrig[it],kAlgName.c_str(),ccent[ic]),200,0,1000);
	hrawpt [it][ic][i]->Sumw2();
      
	hjeteta[it][ic][i] = new TH1F(Form("hjeteta%d_%d_%d",it,ic,i),Form("jet eta distribution jet %s %s %s",ctrig[it],kAlgName.c_str(),ccent[ic]),36,-ketacut,ketacut);
	hjeteta[it][ic][i]->Sumw2();
	hjetphi[it][ic][i] = new TH1F(Form("hjetphi%d_%d_%d",it,ic,i),Form("jet phi distribution jet %s %s %s",ctrig[it],kAlgName.c_str(),ccent[ic]),36,-pi,pi);
	hjetphi[it][ic][i]->Sumw2();


	hchargedmax[it][ic][i] = new TProfile(Form("hchargedmax%d_%d_%d",it,ic,i),"chargedMax",nbins,ptbins);
	hchargedmax[it][ic][i]->Sumw2();
	hphotonmax [it][ic][i] = new TProfile(Form("hphotonmax%d_%d_%d",it,ic,i),"photonMax",nbins,ptbins);
	hphotonmax [it][ic][i]->Sumw2();
	hneutralmax[it][ic][i] = new TProfile(Form("hneutralmax%d_%d_%d",it,ic,i),"neutralMax",nbins,ptbins);
	hneutralmax[it][ic][i]->Sumw2();
	hmuonmax[it][ic][i] = new TProfile(Form("hmuonmax%d_%d_%d",it,ic,i),"muonMax",nbins,ptbins);
	hmuonmax[it][ic][i]->Sumw2();
	helecmax[it][ic][i] = new TProfile(Form("helecmax%d_%d_%d",it,ic,i),"elecMax",nbins,ptbins);
	helecmax[it][ic][i]->Sumw2();

	hchargedsum[it][ic][i] = new TProfile(Form("hchargedsum%d_%d_%d",it,ic,i),"chargedSum",nbins,ptbins);
	hchargedsum[it][ic][i]->Sumw2();
	hphotonsum [it][ic][i] = new TProfile(Form("hphotonsum%d_%d_%d",it,ic,i),"photonSum",nbins,ptbins);
	hphotonsum [it][ic][i]->Sumw2();
	hneutralsum[it][ic][i] = new TProfile(Form("hneutralsum%d_%d_%d",it,ic,i),"neutralSum",nbins,ptbins);
	hneutralsum[it][ic][i]->Sumw2();
	hmuonsum[it][ic][i] = new TProfile(Form("hmuonsum%d_%d_%d",it,ic,i),"muonSum",nbins,ptbins);
	hmuonsum[it][ic][i]->Sumw2();
	helecsum[it][ic][i] = new TProfile(Form("helecsum%d_%d_%d",it,ic,i),"elecSum",nbins,ptbins);
	helecsum[it][ic][i]->Sumw2();

	havbkgpt[it][ic][i] = new TH2F(Form("havbkgpt%d_%d_%d",it,ic,i),Form("(<bkg>) jet p_{T} %s %s %s",ctrig[it],kAlgName.c_str(),ccent[ic]),nbins,ptbins,200,0.,1000.);
	havbkgpt[it][ic][i]->Sumw2();
    
	for(int ie=0;ie<neta;ie++){      
	  havbkgpt_eta[it][ic][i][ie] = new TH2F(Form("havbkgpt_eta%d_%d_%d_%d",it,ic,i,ie),Form("(<bkg>) jet p_{T} %s %s %s",ctrig[it],kAlgName.c_str(),ccent[ic]),nbins,ptbins,200,0.,1000.);
	  havbkgpt_eta[it][ic][i][ie]->Sumw2();
	}
	
	for(int ij=0;ij<nphi;ij++){      
	  havbkgpt_phi[it][ic][i][ij] = new TH2F(Form("havbkgpt_phi%d_%d_%d_%d",it,ic,i,ij),Form("(<bkg>) jet p_{T} %s %s %s",ctrig[it],kAlgName.c_str(),ccent[ic]),nbins,ptbins,200,0.,1000.);
	  havbkgpt_phi[it][ic][i][ij]->Sumw2();
	}
      }
    }//! ic
  }//! itrig
  fout->cd("../");
  std::cout<<"Initialized the histograms " <<std::endl;


  // //std::vector<double> vJets;
  Long64_t nbytes=0;
  Long64_t nentries = tch_pfjet->GetEntries();
  std::cout<<Form("# of entries in TTree for %s %s : ",kAlgName.c_str(),ksp)<<nentries<<std::endl;

  double nevt [ntrig][ncen]={0};
  double njets[ntrig][ncen]={0};

  for (Long64_t i=0; i<nentries;i++) {
    nbytes += tch_pfjet->GetEntry(i);

    //if(i==10)break;
    
    double wxs=1;
    double wcen=1;
    double wvz=1;

    float rndm=gRandom->Rndm();

    hEvents_Total->Fill(rndm);
    //std::cout<<" ********** Event # " <<i<<"\t vz  : "<<vz<<"\t hiBin : "<<hiBin<<"\t  # of PF jets "<<nref<<" # of calojets  "<<nref_calo<<std::endl;

    if( pcollisionEventSelection )hEvents_pCollEvent->Fill(rndm);
    if( pcollisionEventSelection && pHBHENoiseFilter )hEvents_pHBHENoise->Fill(rndm);
    if( pcollisionEventSelection && pHBHENoiseFilter && fabs(vz)<15. )hEvents_Vzcut->Fill(rndm);

    if( pcollisionEventSelection==0 || pHBHENoiseFilter==0 || fabs(vz) > 15. )continue;

    int jetCounter = 0;//counts jets which are going to be used in the supernova cut rejection.
    for(int g = 0;g<nref; g++){
      if(fabs(jteta[g]) < 2){ //to select inside
	if(jtpt[g]>=50) jetCounter++;
      }//eta selection cut
    }// jet loop
    // apply the correct supernova selection cut rejection here:
    if(hiNpix > 38000 - 500*jetCounter)continue;

    if(jet80_1 || jet55_1 || jet65_1){
      hBin_full->Fill(hiBin,wxs*wcen*wvz);
      hVz_full->Fill(vz,wxs*wcen*wvz);
    }

    int iCent = GetCentBin(hiBin);
    if( iCent<0 || iCent>=ncen )continue;

    hEvents_Cent->Fill(rndm);

    if(i%50000==0)std::cout<<" ********** Event # " <<i<<"\t vz  : "<<vz<<"\t hiBin : "<<hiBin<<"\t  # of PF jets "<<nref<< "  # of Calo jets " <<nref_calo<<std::endl;
    
    int isFired[ntrig]={0};

    if( jet80_1   )isFired[0]=1;
    if( jet80_p_1 )isFired[1]=1;

    if( jet65_1   )isFired[2]=1;
    if( jet65_p_1 )isFired[3]=1;

    if( jet55_1   )isFired[4]=1;
    if( jet55_p_1 )isFired[5]=1;

    if( jet80_1 && !jet80_p_1 )isFired[6]=1;
    if( jet65_1 && !jet65_p_1 )isFired[7]=1;
    if( jet55_1 && !jet55_p_1 )isFired[8]=1;

    if( jet55_1 &&  jet65_1)isFired[9]=1;
    if( jet65_1 &&  jet80_1)isFired[10]=1;
    if( jet55_1 &&  jet80_1)isFired[11]=1;

    if( jet65_1 && !jet80_1)isFired[12]=1;
    if(!jet55_1 &&  jet65_1)isFired[13]=1;
    if( jet55_1 && !jet65_1)isFired[14]=1;	  

    // if( jet65_1 && !jet80_1)isFired[4]=1;
    // if( jet65_p_1 && !jet80_p_1)isFired[5]=1;
    // if( jet55_1 && !jet65_1 && !jet80_1)isFired[8]=1;
    // if( jet55_p_1 && !jet65_p_1 && !jet80_p_1)isFired[9]=1;
    // if(!jet55_1 &&  jet65_1 && !jet80_1 )isFired[10]=1;
    // if( jet55_1 &&  jet65_1 && !jet80_1 )isFired[11]=1;
    // if(!jet55_1 &&  jet65_1 &&  jet80_1 )isFired[12]=1;
    // if( jet55_1 &&  jet65_1 &&  jet80_1 )isFired[13]=1;
    // if( jet80_1 &&  L1_sj52_1)isFired[14]=1;
    // if( jet65_1 &&  L1_sj36_1 && !jet80_1)isFired[15]=1;
    // if( jet55_1 &&  L1_sj36_1 && !jet65_1 && !jet80_1)isFired[16]=1;
    // if( jet55_1 ||  jet65_1 || jet80_1)isFired[17]=1;



    std::vector < Jet > pfjet_coll;
    if(isFired[17]==1){
      if(nref_calo>0)hBin_NJets->Fill(hiBin,(float)nref/nref_calo);
      else if(nref_calo==0)hBin_NJets->Fill(hiBin, -0.5);
      else if(nref==0)hBin_NJets->Fill(hiBin, 100.5);

      for(int gj=0; gj<nref; gj++){
	
	if( rawpt[gj] < kptrawcut || jtpt[gj] > 600. ) continue;
	if( fabs(jteta[gj]) > ketacut ) continue;
	
	if( (eMax[gj]/jtpt[gj])>=0.6 || (chargedMax[gj]/jtpt[gj])<=0.02 )continue;
	
	float recopt  = jtpt[gj];
	float recoeta = jteta[gj];
	float recophi = jtphi[gj];

	Jet pfjet;
	pfjet.pt  = recopt;
        pfjet.eta = recoeta;
        pfjet.phi = recophi;
        pfjet.id  = gj;
	pfjet_coll.push_back(pfjet);
      }
    }

    std::vector <Jet>::const_iterator ijet;
    //! Find data driven JEC
    if(pfjet_coll.size()>1){
      std::sort(pfjet_coll.begin(), pfjet_coll.end(), compare_pt);
      ijet =pfjet_coll.begin();
      Jet ljet  =  *(ijet);
      Jet sljet =  *(ijet + 1);
      double dphi=-9999;
      if ( ljet.pt > 60. && sljet.pt > 60. ){
        dphi = delphi(ljet.phi, sljet.phi);
        if ( dphi > 2.*pi/3. ){ //! matched dijet
          double ptdij = (ljet.pt + sljet.pt)/2.;
          int mstat=1;
          if (pfjet_coll.size()>2){
            Jet thjet =  *(ijet + 2);
            if(thjet.pt/ptdij > 0.2)mstat=0;
          }
          if(mstat){
            double B=-9999;
            double rn1 = gRandom->Rndm();
            double rn2 = gRandom->Rndm();
            //cout  << " leading jet : "  << ljet.pt << "  subleading jet : "  <<sljet.pt << endl; 
            if(rn1 > rn2){
              B = (ljet.pt - sljet.pt)/(ljet.pt + sljet.pt);
            }else{
              B = (sljet.pt - ljet.pt)/(sljet.pt + ljet.pt);
            }
            hrescrpt_datadr[0][iCent]->Fill(ptdij, B, wxs*wcen*wvz);
          }
        }
      }
    }




    for(int it=0; it<ntrig; it++){
      if(isFired[it]==0)continue;
      
      hBin [it]->Fill(hiBin,wxs*wcen*wvz);
      hVz[it]->Fill(vz,wxs*wcen*wvz);
      hEvts[it]->Fill(it);
      nevt [it][iCent]++;

      for(int gj=0; gj<nref; gj++){
	
      if( rawpt[gj] < kptrawcut || jtpt[gj] > 600. ) continue;
      if( fabs(jteta[gj]) > ketacut ) continue;

      int jetId=1;
      if( (eMax[gj]/jtpt[gj])>=0.6 || (chargedMax[gj]/jtpt[gj])<=0.02 )jetId=0;

      float recopt  = jtpt[gj];
      float recoeta = jteta[gj];
      float recophi = jtphi[gj];

      double rchmax = chargedMax[gj] / recopt;
      double rnemax = neutralMax[gj] / recopt;
      double rphmax = photonMax[gj] / recopt;
      double rmumax = muMax[gj] / recopt;
      double relmax = eMax[gj] / recopt;

      hchargedfrc[it][iCent]->Fill(recopt, rchmax,wxs*wcen*wvz);
      hneutralfrc[it][iCent]->Fill(recopt, rnemax,wxs*wcen*wvz);
      hphotonfrc [it][iCent]->Fill(recopt, rphmax,wxs*wcen*wvz);
      hmuonfrc   [it][iCent]->Fill(recopt, rmumax,wxs*wcen*wvz);
      helecfrc   [it][iCent]->Fill(recopt, relmax,wxs*wcen*wvz);

      hchargedmax[it][iCent][jetId]->Fill(recopt, chargedMax[gj],wxs*wcen*wvz);
      hneutralmax[it][iCent][jetId]->Fill(recopt, neutralMax[gj],wxs*wcen*wvz);
      hphotonmax [it][iCent][jetId]->Fill(recopt, photonMax[gj],wxs*wcen*wvz);
      hmuonmax   [it][iCent][jetId]->Fill(recopt, muMax[gj],wxs*wcen*wvz);
      helecmax   [it][iCent][jetId]->Fill(recopt, eMax[gj],wxs*wcen*wvz);

      hchargedsum[it][iCent][jetId]->Fill(recopt, chargedSum[gj],wxs*wcen*wvz);
      hneutralsum[it][iCent][jetId]->Fill(recopt, neutralSum[gj],wxs*wcen*wvz);
      hphotonsum [it][iCent][jetId]->Fill(recopt, photonSum[gj],wxs*wcen*wvz);
      hmuonsum   [it][iCent][jetId]->Fill(recopt, muSum[gj],wxs*wcen*wvz);
      helecsum   [it][iCent][jetId]->Fill(recopt, eSum[gj],wxs*wcen*wvz);
      
      njets[it][iCent]++;
      
      
      int iphi = GetPhiBin(recophi);
      int ieta = GetEtaBin(recoeta);
      
      double sumPF = ( chargedSum[gj] + neutralSum[gj] + photonSum[gj] + eSum[gj] + muSum[gj] );
      double bkgd  = sumPF - rawpt[gj];

      havbkgpt[it][iCent][jetId]->Fill( recopt , bkgd, wxs*wcen*wvz );
      if( ieta>=0 && ieta<neta)havbkgpt_eta[it][iCent][jetId][ieta]->Fill( recopt , bkgd, wxs*wcen*wvz);
      if( iphi>=0 && iphi<nphi)havbkgpt_phi[it][iCent][jetId][iphi]->Fill( recopt , bkgd, wxs*wcen*wvz);

      hrecopt   [it][iCent][jetId]->Fill(recopt,wxs*wcen*wvz);
      hrawpt    [it][iCent][jetId]->Fill(rawpt[gj],wxs*wcen*wvz);
      hrawrecpt [it][iCent][jetId]->Fill(rawpt[gj],recopt,wxs*wcen*wvz);
      
      hjeteta   [it][iCent][jetId]->Fill(recoeta,wxs*wcen*wvz);
      hjetphi   [it][iCent][jetId]->Fill(recophi,wxs*wcen*wvz);
      
      }//! irec loop
    }//! ntrig loop
    //std::cout<<"Completed event #  "<<ievt<<std::endl; 
  }//! event loop ends
  

  std::cout<<std::endl;
  std::cout<<kAlgName.c_str() << std::endl;
  for(int it=0; it<ntrig; it++){
    cout <<ctrig[it]<<endl;
    for(int ic=0;ic<ncen;ic++){
      std::cout<<"\t cent : " << ccent[ic] << "  # of events PbPb : " << nevt[it][ic] << " # of jets  " << njets[it][ic] << std::endl;      
    }
  }
  std::cout<<std::endl;

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
double delphi(double phi1, double phi2)
{
  double dphi = phi2 - phi1;
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
  return ibin;
}
double deltaR(float eta1, float phi1, float eta2, float phi2)
{
  float deta = eta1 - eta2;
  float dphi = fabs(phi1 - phi2);
  if(dphi > pi)dphi -= 2*pi;
  double dr = sqrt(pow(deta,2) + pow(dphi,2));
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
