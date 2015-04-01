
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


using namespace std;

#ifdef _MAKECINT_
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#endif

int GetPhiBin(float /*phi*/);
int GetEtaBin(float /*eta*/);
int GetPtBin(float /*pt*/);
int GetPtBinWide(float /*pt*/);
int GetCentBin(int /*hiBin*/);
double delphi(double /*phi1*/, double /*phi2*/);
double GetXsec(double /*maxpthat*/);
double deltaR(float /*eta1*/, float /*phi1*/, 
	      float /*eta2*/, float /*phi2*/);
void GetCentWeight(TH1F */*hCentWeight*/);

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

//! pp 2.76 TeV PYTHIA di-jet cross sections
// const double xsec[12][3] ={{2.034e-01,15,30}, //! 15
//                            {1.075e-02,30,50}, //! 30
//                            {1.025e-03,50,80}, //! 50
//                            {9.865e-05,80,120}, //! 80
//                            {1.129e-05,120,170}, //! 120
//                            {1.465e-06,170,220}, //! 170
//                            {2.837e-07,220,280}, //! 220
//                            {5.323e-08,280,370}, //! 280
//                            {5.934e-09,370,460}, //! 370
//                            {8.125e-10,460,540}, //! 460
//                            {1.467e-10,540,9999}, //! 540 
//                            {0,9999,9999}
// };

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

const int ncen=6;
const char *cdir[ncen]  = {"05","510","1030","3050","5070","7090"};
const char *ccent[ncen] = {"0-5%","5-10%","10-30%","30-50%","50-70%","70-90%"};


//! constants
int iYear=2014;
const double pi=acos(-1.);
//const double pi2=2*pi -1;

bool wJetId = true;
const double ketacut=2.0;
const double kptrawcut =0.0;
const double kptrecocut=0.0;
const double kptgencut =20.0;
const double kdRCut=0.30;

int rbins=50;
double rbinl=0.0;
double rbinh=2.0;

//! pt binning
double ptbins[] ={30,40,50,60,70,80,90,100,110,120,130,140,160,200,250,300,400,548};
const int nbins = sizeof(ptbins)/sizeof(double) - 1;

//! data pt binning
//double ptbins_data[] ={100, 110, 120, 130, 140, 150, 160, 170, 180, 200, 240, 300};
//const int dbins = sizeof(ptbins_data)/sizeof(double) - 1;

const double ptbins_wide[]={30,50,80,120,200,340,548};
const int nbins_wide = sizeof(ptbins_wide)/sizeof(double) - 1;

//const double etabins[] = {-3.000, -2.400, -2.000, -1.3000, 0.000, 1.300, 2.000, 2.400, 3.000};
const double etabins[] = {-2.000, -1.4000, -0.4500, 0.000, 0.4500, 1.400, 2.000};
const int neta = sizeof(etabins)/sizeof(double) - 1;

//const double etabins[]={-3.0,-2.4,-1.8,-1.4,-1.0,-0.8,-0.6,-0.4,-0.2,
//			0.0,0.2,0.4,0.6,0.8,1.0,1.4,1.8,2.4,3.0};
//const int neta = sizeof(etabins)/sizeof(double) - 1;


//const double phibins[] = {-3.141,-2.700,-2.100,-1.500,-0.900,-0.300, 
// 			  0.300,0.900,1.500,2.100,2.700,3.141
//};
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
int CondorPbPb_MC(std::string kAlgName="akPu3",
		  std::string finame="HiForest_PYTHIA_HYDJET_pthat220_Track9_Jet30_matchEqR_merged_forest_0.root",
		  std::string foname="akPu3_Histo_PYTHIA_HYDJET_pthat220_Track9_Jet30_matchEqR_merged_forest_0.root", 
		  double maxpthat=280
		  )
{
  
  timer.Start();



  //! Doga's
  std::string indir="/mnt/hadoop/cms/store/user/dgulhan/PYTHIA_HYDJET_Track9_Jet30_Pyquen_DiJet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_merged/";
  
  //! Marguerite
  //std::string indir="/mnt/hadoop/cms/store/user/belt/Validation53X/Pyquen_Dijet_TuneZ2_Unquenched_Hydjet1p8_2760GeV_Track9_Jet30_v15_full/";

  std::string inname=indir+finame;


  TFile *fin = new TFile(inname.c_str(),"r");
  //std::string outdir="/net/hidsk0001/d00/scratch/pawan/condorfiles/pbpb/Doga/";
  //std::string outdir="/net/hidsk0001/d00/scratch/pawan/condorfiles/pbpb/Marguerite/";
  std::string outdir="";
  std::string outname=outdir+foname;
  TFile *fout = new TFile(outname.c_str(),"RECREATE");

  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"**************************************************** "<<std::endl;
  std::cout<<Form("PbPbCalJec : %s %s ",ksp, kAlgName.c_str())<<std::endl;
  std::cout<<Form("Outfile : %s",outname.c_str())<<std::endl;
  std::cout<<Form("eta cut for %0.3f ",ketacut)<<std::endl;
  std::cout<<"**************************************************** "<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;

  //! 
  //! Define histograms here
  //   TH1::SetDefaultSumw2();
  //   TH2::SetDefaultSumw2();
  //   TH3::SetDefaultSumw2();
  //  TProfile::SetDefaultSumw2();

  fout->mkdir(Form("%sJetAnalyzer",kAlgName.c_str()));
  fout->cd(Form("%sJetAnalyzer"   ,kAlgName.c_str()));

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  TH1F *hVz   = new TH1F("hVz","vertex z",80,-20.0,20.0);
  hVz->Sumw2();
  TH1F *hBin    = new TH1F("hBin","Centrality bin",200,-0.5,200-0.5);
  hBin->Sumw2();
  TH1F *hpthat = new TH1F("hpthat","pt-hat distribution",200,0,1000);
  hpthat->Sumw2();

  TH1F *hchargedmax[ncen];
  TH1F *hphotonmax[ncen];
  TH1F *hneutralmax[ncen];
  TH1F *hmuonmax[ncen];
  TH1F *helecmax[ncen];
  

  //! Gen matched jets
  TH1F *hgenpt_genm [2][ncen], *hrecopt_genm[2][ncen], *hrawpt_genm[2][ncen];
  TH1F *hgeneta     [2][ncen], *hgenphi     [2][ncen];
  TH1F *hjeteta     [2][ncen], *hjetphi     [2][ncen];

  TH2F *hgenrecpt_genm[2][ncen], *hgenrawpt_genm[2][ncen];

  TH1F *havbkgpt_genm[2][ncen][nbins];
  TH1F *havbkgpt_eta_genm[2][ncen][nbins][2];

  //! Resposnse
  TH1F *hrescrpt_genm [2][ncen][nbins], *hresrrpt_genm [2][ncen][nbins];
  TH1F *hrescrpt_genm_eta [2][ncen][nbins][neta], *hresrrpt_genm_eta [2][ncen][nbins][neta];
  TH1F *hrescrpt_genm_phi [2][ncen][nbins][nphi], *hresrrpt_genm_phi [2][ncen][nbins][nphi];
  TH1F *hrescrpt_wide_genm_eta [2][ncen][nbins_wide][neta], *hresrrpt_wide_genm_eta [2][ncen][nbins_wide][neta];
  TH1F *hrescrpt_wide_genm_phi [2][ncen][nbins_wide][nphi], *hresrrpt_wide_genm_phi [2][ncen][nbins_wide][nphi];
 
  TH1F *havrefpt_genm[2][ncen][nbins], *havjetpt_genm[2][ncen][nbins], *havrawpt_genm[2][ncen][nbins];
  TProfile *hrawpt_genpt [2][ncen], *hrecopt_genpt[2][ncen];

  //!
  TH1F *hPtAll [2][ncen], *hPtEff [2][ncen];
  TH1F *hEtaAll[2][ncen], *hEtaEff[2][ncen];
  TH1F *hPhiAll[2][ncen], *hPhiEff[2][ncen];

  TH1F *hPtFakeAll[2][ncen], *hPtFake[2][ncen];
  TH1F *hEtaFakeAll_20[2][ncen],*hEtaFakeAll_30[2][ncen], *hEtaFakeAll_40[2][ncen],*hEtaFakeAll_45[2][ncen], *hEtaFakeAll_50[2][ncen];
  TH1F *hPhiFakeAll_20[2][ncen],*hPhiFakeAll_30[2][ncen], *hPhiFakeAll_40[2][ncen],*hPhiFakeAll_45[2][ncen], *hPhiFakeAll_50[2][ncen];
  TH1F *hEtaFake_20[2][ncen],*hEtaFake_30[2][ncen], *hEtaFake_40[2][ncen],*hEtaFake_45[2][ncen], *hEtaFake_50[2][ncen];
  TH1F *hPhiFake_20[2][ncen],*hPhiFake_30[2][ncen], *hPhiFake_40[2][ncen],*hPhiFake_45[2][ncen], *hPhiFake_50[2][ncen];

  //! 
  TH1F *hsubid[ncen];
  TH1F *hPt[ncen][5], *hRawPt[ncen][5];

  //! Data-driven JEC's
  TH2F *hrescrpt_datadr[2][ncen];


  for(int ic=0; ic<ncen; ic++){
    
    hsubid[ic] = new TH1F(Form("hsubid_%d",ic),"Sumid :",100,-2.0-0.5,103-0.5);
    hsubid[ic]->Sumw2();
    for(int i=0; i<5;i++){
      hPt[ic][i] = new TH1F(Form("hPt%d_%d",ic,i),Form(" pT for algorithm %s %s %d",kAlgName.c_str(),ccent[ic],i),50,20,620);
      hPt[ic][i]->Sumw2();
      hRawPt[ic][i] = new TH1F(Form("hRawPt%d_%d",ic,i),Form("Raw  pT for algorithm %s %s %d",kAlgName.c_str(),ccent[ic],i),50,20,620);
      hRawPt[ic][i]->Sumw2();
    }
    hchargedmax[ic] = new TH1F(Form("hchargedmax_%d",ic),"chargedMax/recopt",100,0.,5.);
    hchargedmax[ic]->Sumw2();
    hphotonmax[ic] = new TH1F(Form("hphotonmax_%d",ic),"photonMax/recopt",100,0.,5.);
    hphotonmax[ic]->Sumw2();
    hneutralmax[ic] = new TH1F(Form("hneutralmax_%d",ic),"neutralMax/recopt",100,0.,5.);
    hneutralmax[ic]->Sumw2();
    hmuonmax[ic] = new TH1F(Form("hmuonmax_%d",ic),"muonMax/recopt",100,0.,5.);
    hmuonmax[ic]->Sumw2();
    helecmax[ic] = new TH1F(Form("helecmax_%d",ic),"elecMax/recopt",100,0.,5.);
    helecmax[ic]->Sumw2();
  }

  for(int id=0; id<2; id++){//! 0:PF and 1:Calo
    for(int ic=0; ic<ncen; ic++){
      hgenrecpt_genm [id][ic] = new TH2F(Form("hgenrecpt_genm%d_%d",id,ic),Form("Gen-Rec matched gen p_{T} distribution jet centb %s %s",kAlgName.c_str(),ccent[ic]),200,0,1000,200,0,1000);
      hgenrecpt_genm [id][ic]->Sumw2();
      hgenrawpt_genm [id][ic] = new TH2F(Form("hgenrawpt_genm%d_%d",id,ic),Form("Gen-Raw matched gen p_{T} distribution jet centb %s %s",kAlgName.c_str(),ccent[ic]),200,0,1000,200,0,1000);
      hgenrawpt_genm [id][ic]->Sumw2();
      
      //hgenrecpt_genm [id][ic] = new TH2F(Form("hgenrecpt_genm%d_%d",id,ic),Form("Gen-Rec matched gen p_{T} distribution jet centb %s %s",kAlgName.c_str(),ccent[ic]),nbins,ptbins,nbins,ptbins);
      //hgenrecpt_genm [id][ic]->Sumw2();
  
      hgenpt_genm [id][ic] = new TH1F(Form("hgenpt_genm%d_%d",id,ic),Form("Gen matched gen p_{T} distribution jet centb %s %s",kAlgName.c_str(),ccent[ic]),200,0,1000);
      hgenpt_genm [id][ic]->Sumw2();
      hrecopt_genm[id][ic] = new TH1F(Form("hrecopt_genm%d_%d",id,ic),Form("Gen matched reco p_{T} distribution jet %s %s",kAlgName.c_str(),ccent[ic]),200,0,1000);
      hrecopt_genm[id][ic]->Sumw2();
      hrawpt_genm [id][ic] = new TH1F(Form("hrawpt_genm%d_%d",id,ic),Form("Gen matched raw p_{T} distribution jet  %s %s",kAlgName.c_str(),ccent[ic]),200,0,1000);
      hrawpt_genm [id][ic]->Sumw2();
  
      hjeteta[id][ic] = new TH1F(Form("hjeteta%d_%d",id,ic),Form("jet eta distribution jet %s %s",kAlgName.c_str(),ccent[ic]),36,-ketacut,ketacut);
      hjeteta[id][ic]->Sumw2();
      hjetphi[id][ic] = new TH1F(Form("hjetphi%d_%d",id,ic),Form("jet phi distribution jet %s %s",kAlgName.c_str(),ccent[ic]),36,-pi,pi);
      hjetphi[id][ic]->Sumw2();
  
      hgeneta[id][ic] = new TH1F(Form("hgeneta%d_%d",id,ic),Form("gen jet eta distribution jet %s %s",kAlgName.c_str(),ccent[ic]),36,-ketacut,ketacut);
      hgeneta[id][ic]->Sumw2();
      hgenphi[id][ic] = new TH1F(Form("hgenphi%d_%d",id,ic),Form("gen jet phi distribution jet %s %s",kAlgName.c_str(),ccent[ic]),36,-pi,pi);
      hgenphi[id][ic]->Sumw2();
  
      hrawpt_genpt [id][ic] = new TProfile(Form("hrawpt_genpt%d_%d",id,ic),Form("reco p_{T} gen p_{T} %s %s",kAlgName.c_str(),ccent[ic]),nbins,ptbins);
      hrawpt_genpt [id][ic]->Sumw2();
      hrecopt_genpt[id][ic] = new TProfile(Form("hrecopt_genpt%d_%d",id,ic),Form("reco p_{T} gen p_{T} %s %s",kAlgName.c_str(),ccent[ic]),nbins,ptbins);
      hrecopt_genpt[id][ic]->Sumw2();

      //! efficiency histograms
      //hPtAll [ic] = new TH1F(Form("hPtAll_%d_%d",id,ic),Form("Denominator pT for algorithm %s %s",kAlgName.c_str(),ccent[ic]),nbins,ptbins);
      hPtAll [id][ic] = new TH1F(Form("hPtAll_%d_%d",id,ic),Form("Denominator pT for algorithm %s %s",kAlgName.c_str(),ccent[ic]),50,20,620);
      hPtAll [id][ic]->Sumw2();
      hEtaAll[id][ic] = new TH1F(Form("hEtaAll_%d_%d",id,ic),Form("Denominator eta  for algorithm %s %s",kAlgName.c_str(),ccent[ic]), 20, -2.0,2.0);
      hEtaAll[id][ic]->Sumw2();
      hPhiAll[id][ic] = new TH1F(Form("hPhiAll_%d_%d",id,ic),Form("Denominator  phi  for algorithm %s %s",kAlgName.c_str(),ccent[ic]),18,-pi,pi);
      hPhiAll[id][ic]->Sumw2();

      hPtFakeAll [id][ic] = new TH1F(Form("hPtFakeAll_%d_%d",id,ic),Form("Fake All Denominator pT for algorithm %s %s",kAlgName.c_str(),ccent[ic]),50,20,620);
      hPtFakeAll [id][ic]->Sumw2();

      hEtaFakeAll_20[id][ic] = new TH1F(Form("hEtaFakeAll_20_%d_%d",id,ic),Form("Fake All 20 GeV/c Denominator eta  for algorithm %s %s",kAlgName.c_str(),ccent[ic]), 20, -2.0,2.0);
      hEtaFakeAll_20[id][ic]->Sumw2();
      hPhiFakeAll_20[id][ic] = new TH1F(Form("hPhiFakeAll_20_%d_%d",id,ic),Form("Fake All 20 GeV/c Denominator  phi  for algorithm %s %s",kAlgName.c_str(),ccent[ic]),18,-pi,pi);
      hPhiFakeAll_20[id][ic]->Sumw2();


      hEtaFakeAll_30[id][ic] = new TH1F(Form("hEtaFakeAll_30_%d_%d",id,ic),Form("Fake All 30 GeV/c Denominator eta  for algorithm %s %s",kAlgName.c_str(),ccent[ic]), 20, -2.0,2.0);
      hEtaFakeAll_30[id][ic]->Sumw2();
      hPhiFakeAll_30[id][ic] = new TH1F(Form("hPhiFakeAll_30_%d_%d",id,ic),Form("Fake All 30 GeV/c Denominator  phi  for algorithm %s %s",kAlgName.c_str(),ccent[ic]),18,-pi,pi);
      hPhiFakeAll_30[id][ic]->Sumw2();


      hEtaFakeAll_40[id][ic] = new TH1F(Form("hEtaFakeAll_40_%d_%d",id,ic),Form("Fake All 40 GeV/c Denominator eta  for algorithm %s %s",kAlgName.c_str(),ccent[ic]), 20, -2.0,2.0);
      hEtaFakeAll_40[id][ic]->Sumw2();
      hPhiFakeAll_40[id][ic] = new TH1F(Form("hPhiFakeAll_40_%d_%d",id,ic),Form("Fake All 40 GeV/c Denominator  phi  for algorithm %s %s",kAlgName.c_str(),ccent[ic]),18,-pi,pi);
      hPhiFakeAll_40[id][ic]->Sumw2();


      hEtaFakeAll_45[id][ic] = new TH1F(Form("hEtaFakeAll_45_%d_%d",id,ic),Form("Fake All 45 GeV/c Denominator eta  for algorithm %s %s",kAlgName.c_str(),ccent[ic]), 20, -2.0,2.0);
      hEtaFakeAll_45[id][ic]->Sumw2();
      hPhiFakeAll_45[id][ic] = new TH1F(Form("hPhiFakeAll_45_%d_%d",id,ic),Form("Fake All 45 GeV/c Denominator  phi  for algorithm %s %s",kAlgName.c_str(),ccent[ic]),18,-pi,pi);
      hPhiFakeAll_45[id][ic]->Sumw2();

      hEtaFakeAll_50[id][ic] = new TH1F(Form("hEtaFakeAll_50_%d_%d",id,ic),Form("Fake All 50 GeV/c Denominator eta  for algorithm %s %s",kAlgName.c_str(),ccent[ic]), 20, -2.0,2.0);
      hEtaFakeAll_50[id][ic]->Sumw2();
      hPhiFakeAll_50[id][ic] = new TH1F(Form("hPhiFakeAll_50_%d_%d",id,ic),Form("Fake All 50 GeV/c Denominator  phi  for algorithm %s %s",kAlgName.c_str(),ccent[ic]),18,-pi,pi);
      hPhiFakeAll_50[id][ic]->Sumw2();

      hPtEff [id][ic] = new TH1F(Form("hPtEff_%d_%d",id,ic),Form("Neunominator eff pT for algorithm %s %s",kAlgName.c_str(),ccent[ic]),50,20,620);
      hPtEff [id][ic]->Sumw2();
      hEtaEff[id][ic] = new TH1F(Form("hEtaEff_%d_%d",id,ic),Form("Neunominator eff eta  for algorithm %s %s",kAlgName.c_str(),ccent[ic]), 20, -2.0,2.0);
      hEtaEff[id][ic]->Sumw2();
      hPhiEff[id][ic] = new TH1F(Form("hPhiEff_%d_%d",id,ic),Form("Neunominator  eff phi  for algorithm %s %s",kAlgName.c_str(),ccent[ic]),18,-pi,pi);
      hPhiEff[id][ic]->Sumw2();

      hPtFake [id][ic] = new TH1F(Form("hPtFake_%d_%d",id,ic),Form("Neunominator fake pT for algorithm %s %s",kAlgName.c_str(),ccent[ic]),50,20,620);
      hPtFake [id][ic]->Sumw2();
    
      hEtaFake_20[id][ic] = new TH1F(Form("hEtaFake_20_%d_%d",id,ic),Form("Neunominator 20 fake eta  for algorithm %s %s",kAlgName.c_str(),ccent[ic]), 20, -2.0,2.0);
      hEtaFake_20[id][ic]->Sumw2();
      hPhiFake_20[id][ic] = new TH1F(Form("hPhiFake_20_%d_%d",id,ic),Form("Neunominator 20 fake phi  for algorithm %s %s",kAlgName.c_str(),ccent[ic]),18,-pi,pi);
      hPhiFake_20[id][ic]->Sumw2();

      hEtaFake_30[id][ic] = new TH1F(Form("hEtaFake_30_%d_%d",id,ic),Form("Neunominator 30 fake eta  for algorithm %s %s",kAlgName.c_str(),ccent[ic]), 20, -2.0,2.0);
      hEtaFake_30[id][ic]->Sumw2();
      hPhiFake_30[id][ic] = new TH1F(Form("hPhiFake_30_%d_%d",id,ic),Form("Neunominator 30 fake phi  for algorithm %s %s",kAlgName.c_str(),ccent[ic]),18,-pi,pi);
      hPhiFake_30[id][ic]->Sumw2();

      hEtaFake_40[id][ic] = new TH1F(Form("hEtaFake_40_%d_%d",id,ic),Form("Neunominator 40 fake eta  for algorithm %s %s",kAlgName.c_str(),ccent[ic]), 20, -2.0,2.0);
      hEtaFake_40[id][ic]->Sumw2();
      hPhiFake_40[id][ic] = new TH1F(Form("hPhiFake_40_%d_%d",id,ic),Form("Neunominator 40 fake phi  for algorithm %s %s",kAlgName.c_str(),ccent[ic]),18,-pi,pi);
      hPhiFake_40[id][ic]->Sumw2();

      hEtaFake_45[id][ic] = new TH1F(Form("hEtaFake_45_%d_%d",id,ic),Form("Neunominator 45 fake eta  for algorithm %s %s",kAlgName.c_str(),ccent[ic]), 20, -2.0,2.0);
      hEtaFake_45[id][ic]->Sumw2();
      hPhiFake_45[id][ic] = new TH1F(Form("hPhiFake_45_%d_%d",id,ic),Form("Neunominator 45 fake phi  for algorithm %s %s",kAlgName.c_str(),ccent[ic]),18,-pi,pi);
      hPhiFake_45[id][ic]->Sumw2();

      hEtaFake_50[id][ic] = new TH1F(Form("hEtaFake_50_%d_%d",id,ic),Form("Neunominator 50 fake eta  for algorithm %s %s",kAlgName.c_str(),ccent[ic]), 20, -2.0,2.0);
      hEtaFake_50[id][ic]->Sumw2();
      hPhiFake_50[id][ic] = new TH1F(Form("hPhiFake_50_%d_%d",id,ic),Form("Neunominator 50 fake phi  for algorithm %s %s",kAlgName.c_str(),ccent[ic]),18,-pi,pi);
      hPhiFake_50[id][ic]->Sumw2();

      hrescrpt_datadr [id][ic] = new TH2F(Form("hrescrpt_datadr_%d_%d",id,ic),Form("(Reco/Gen) jet p_{T} %s %s ",kAlgName.c_str(),ccent[ic]), 30, 50, 650, 50,-1.50,1.50);
      hrescrpt_datadr [id][ic]->Sumw2();


      for(int ip=0;ip<nbins;ip++){
	havrefpt_genm [id][ic][ip]= new TH1F(Form("havrefpt_genm%d_%d_%d",id,ic,ip),Form("(<Ref>) jet p_{T} %s %s %0.0f < p_{T}^{REF} < %0.0f",kAlgName.c_str(),ccent[ic],ptbins[ip],ptbins[ip+1]),500,0,1000);
	havjetpt_genm [id][ic][ip]= new TH1F(Form("havjetpt_genm%d_%d_%d",id,ic,ip),Form("(<Reco>) jet p_{T} %s %s %0.0f < p_{T}^{REF} < %0.0f",kAlgName.c_str(),ccent[ic],ptbins[ip],ptbins[ip+1]),500,0,1000);
	havjetpt_genm [id][ic][ip]->Sumw2();
	havrawpt_genm [id][ic][ip]= new TH1F(Form("havrawpt_genm%d_%d_%d",id,ic,ip),Form("(<Raw>) jet p_{T} %s %s %0.0f < p_{T}^{REF} < %0.0f",kAlgName.c_str(),ccent[ic],ptbins[ip],ptbins[ip+1]),500,0,1000);
	havrawpt_genm [id][ic][ip]->Sumw2();
    
	havbkgpt_genm [id][ic][ip]= new TH1F(Form("havbkgpt_genm%d_%d_%d",id,ic,ip),Form("(<bkg>) jet p_{T} %s %s %0.0f < p_{T}^{REF} < %0.0f",kAlgName.c_str(),ccent[ic],ptbins[ip],ptbins[ip+1]),100,0,300);
	havbkgpt_genm [id][ic][ip]->Sumw2();
	for(int ie=0; ie<2;ie++){
	  havbkgpt_eta_genm [id][ic][ip][ie]= new TH1F(Form("havbkgpt_eta_genm%d_%d_%d_%d",id,ic,ip,ie),Form("(<bkg>) jet p_{T} %s %s %0.0f < p_{T}^{RAW} < %0.0f",kAlgName.c_str(),ccent[ic],ptbins[ip],ptbins[ip+1]),100,0,300);
	  havbkgpt_eta_genm [id][ic][ip][ie]->Sumw2();
	}
    
	//! Gen matched Response and resolution
	hrescrpt_genm [id][ic][ip]= new TH1F(Form("hrescrpt_genm%d_%d_%d",id,ic,ip),Form("(Reco/Gen) jet p_{T} %s %s %0.0f < p_{T}^{REF} < %0.0f",kAlgName.c_str(),ccent[ic],ptbins[ip],ptbins[ip+1]),rbins,rbinl,rbinh);
	hrescrpt_genm [id][ic][ip]->Sumw2();
	hresrrpt_genm [id][ic][ip]= new TH1F(Form("hresrrpt_genm%d_%d_%d",id,ic,ip),Form("(Raw/Gen) jet p_{T} %s %s %0.0f < p_{T}^{REF} < %0.0f",kAlgName.c_str(),ccent[ic],ptbins[ip],ptbins[ip+1]),rbins,rbinl,rbinh);
	hresrrpt_genm [id][ic][ip]->Sumw2();

    
	for(int ie=0;ie<neta;ie++){      
	  hrescrpt_genm_eta [id][ic][ip][ie]= new TH1F(Form("hrescrpt_genm_eta%d_%d_%d_%d",id,ic,ip,ie),Form("(Reco/Gen) jet p_{T} %s %s %d %d",kAlgName.c_str(),ccent[ic],ip,ie),rbins,rbinl,rbinh);
	  hrescrpt_genm_eta [id][ic][ip][ie]->Sumw2();
	  hresrrpt_genm_eta [id][ic][ip][ie]= new TH1F(Form("hresrrpt_genm_eta%d_%d_%d_%d",id,ic,ip,ie),Form("(Raw/Gen) jet p_{T} %s %s %d %d",kAlgName.c_str(),ccent[ic],ip,ie),rbins,rbinl,rbinh);
	  hresrrpt_genm_eta [id][ic][ip][ie]->Sumw2();
	}
    
	for(int ij=0;ij<nphi;ij++){      
	  hrescrpt_genm_phi [id][ic][ip][ij]= new TH1F(Form("hrescrpt_genm_phi%d_%d_%d_%d",id,ic,ip,ij),Form("(Reco/Gen) jet p_{T} %s %s %d %d",kAlgName.c_str(),ccent[ic],ip,ij),rbins,rbinl,rbinh);
	  hrescrpt_genm_phi [id][ic][ip][ij]->Sumw2();
	  hresrrpt_genm_phi [id][ic][ip][ij]= new TH1F(Form("hresrrpt_genm_phi%d_%d_%d_%d",id,ic,ip,ij),Form("(Raw/Gen) jet p_{T} %s %s %d %d",kAlgName.c_str(),ccent[ic],ip,ij),rbins,rbinl,rbinh);
	  hresrrpt_genm_phi [id][ic][ip][ij]->Sumw2();
	}
      }

      //! coarse pt bin 
      for(int ip=0;ip<nbins_wide;ip++){
	for(int ie=0;ie<neta;ie++){
	  hrescrpt_wide_genm_eta [id][ic][ip][ie]= new TH1F(Form("hrescrpt_wide_genm_eta%d_%d_%d_%d",id,ic,ip,ie),Form("(Reco/Gen) jet p_{T} %s %d %d",ccent[ic],ip,ie),rbins,rbinl,rbinh);
	  hrescrpt_wide_genm_eta [id][ic][ip][ie]->Sumw2();
	  hresrrpt_wide_genm_eta [id][ic][ip][ie]= new TH1F(Form("hresrrpt_wide_genm_eta%d_%d_%d_%d",id,ic,ip,ie),Form("(Raw/Gen) jet p_{T} %s %d %d",ccent[ic],ip,ie),rbins,rbinl,rbinh);
	  hresrrpt_wide_genm_eta [id][ic][ip][ie]->Sumw2();
	}
	for(int ij=0;ij<nphi;ij++){
	  hrescrpt_wide_genm_phi [id][ic][ip][ij]= new TH1F(Form("hrescrpt_wide_genm_phi%d_%d_%d_%d",id,ic,ip,ij),Form("(Reco/Gen) jet p_{T} %s %d %d",ccent[ic],ip,ij),rbins,rbinl,rbinh);
	  hrescrpt_wide_genm_phi [id][ic][ip][ij]->Sumw2();
	  hresrrpt_wide_genm_phi [id][ic][ip][ij]= new TH1F(Form("hresrrpt_wide_genm_phi%d_%d_%d_%d",id,ic,ip,ij),Form("(Raw/Gen) jet p_{T} %s %d %d",ccent[ic],ip,ij),rbins,rbinl,rbinh);
	  hresrrpt_wide_genm_phi [id][ic][ip][ij]->Sumw2();
	}
      }
    }//! ic
  }//! id
  fout->cd("../");

  std::cout<<"Initialized the histograms " <<std::endl;
  ///////////////////////////////////////////////////////////////////////////////////////// 

  Long64_t nbytes=0;

  TTree *tr_ev = 0;
  float vz;
  int hiBin;

  tr_ev = (TTree*)fin->Get("hiEvtAnalyzer/HiTree");
  tr_ev->SetBranchAddress("vz",&vz);
  tr_ev->SetBranchAddress("hiBin",&hiBin);
  tr_ev->SetBranchStatus("*",0,0);
  tr_ev->SetBranchStatus("vz",1,0);
  tr_ev->SetBranchStatus("hiBin",1,0);

  TTree *tr_skim=0;
  int pcollisionEventSelection;
  tr_skim = (TTree*)fin->Get("skimanalysis/HltTree");
  tr_skim->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
  tr_skim->SetBranchStatus("*",0,0);
  tr_skim->SetBranchStatus("pcollisionEventSelection",1,0);

  TTree *tr_in=0;
  tr_in = (TTree*)fin->Get(Form("%sPFJetAnalyzer/t",kAlgName.c_str()));
  Long64_t nentries = tr_in->GetEntries();
  //std::cout<<Form("# of entries in TTree for %s %s : ",kAlgName.c_str(),ksp)<<nentries<<std::endl;
    
  //Declaration of leaves types PF Jets
  //int hiBin=0;
  int   nref;
  float pthat;
  float weight=1;
  float jtpt[1000];
  float rawpt[1000];
  float jteta[1000];
  float jtphi[1000];
  float neutralSum[1000];
  float chargedSum[1000];
  float photonSum [1000];
  float neutralMax[1000];
  float chargedMax[1000];
  float photonMax[1000];
  float eMax[1000];
  float muMax[1000];
  float refpt[1000];
  float refeta[1000];
  float refphi[1000];
  float refdrjt[1000];
  int subid[1000];
  
  int ngen;
  int genmatchindex[500];
  int gensubid[500];

  tr_in->SetBranchAddress("nref",&nref);
  tr_in->SetBranchAddress("pthat",&pthat);
  tr_in->SetBranchAddress("rawpt",rawpt);
  tr_in->SetBranchAddress("jtpt" ,jtpt);
  tr_in->SetBranchAddress("jteta",jteta);
  tr_in->SetBranchAddress("jtphi",jtphi);
  tr_in->SetBranchAddress("neutralSum",neutralSum);
  tr_in->SetBranchAddress("chargedSum",chargedSum);
  tr_in->SetBranchAddress("photonSum",photonSum);
  tr_in->SetBranchAddress("neutralMax",neutralMax);
  tr_in->SetBranchAddress("chargedMax",chargedMax);
  tr_in->SetBranchAddress("photonMax",photonMax);
  tr_in->SetBranchAddress("muMax",muMax);
  tr_in->SetBranchAddress("eMax",eMax);
  tr_in->SetBranchAddress("refpt",refpt);
  tr_in->SetBranchAddress("refphi",refphi);
  tr_in->SetBranchAddress("refeta",refeta);
  tr_in->SetBranchAddress("refdrjt",refdrjt);
  tr_in->SetBranchAddress("subid",subid);

  tr_in->SetBranchAddress("ngen",&ngen);
  tr_in->SetBranchAddress("genmatchindex",genmatchindex);
  tr_in->SetBranchAddress("gensubid",gensubid);
  
  tr_in->SetBranchStatus("*",0,0);
  tr_in->SetBranchStatus("nref" ,1,0);
  tr_in->SetBranchStatus("pthat",1,0);
  tr_in->SetBranchStatus("rawpt",1,0);
  tr_in->SetBranchStatus("jtpt" ,1,0);
  tr_in->SetBranchStatus("jteta",1,0);
  tr_in->SetBranchStatus("jtphi",1,0);
  tr_in->SetBranchStatus("neutralSum",1,0);
  tr_in->SetBranchStatus("chargedSum",1,0);
  tr_in->SetBranchStatus("photonSum",1,0);
  tr_in->SetBranchStatus("neutralMax",1,0);
  tr_in->SetBranchStatus("chargedMax",1,0);
  tr_in->SetBranchStatus("photonMax",1,0);
  tr_in->SetBranchStatus("muMax",1,0);
  tr_in->SetBranchStatus("eMax",1,0);
  tr_in->SetBranchStatus("refpt",1,0);
  tr_in->SetBranchStatus("refphi",1,0);
  tr_in->SetBranchStatus("refeta",1,0);
  tr_in->SetBranchStatus("refdrjt",1,0);
  tr_in->SetBranchStatus("subid",1,0);
  
  tr_in->SetBranchStatus("ngen",1,0);
  tr_in->SetBranchStatus("genmatchindex",1,0);
  tr_in->SetBranchStatus("gensubid",1,0);


  //! Calo jets
  //Declaration of leaves types                                                                                                               
  int   nref_calo;
  float jtpt_calo[1000];
  float rawpt_calo[1000];
  float jteta_calo[1000];
  float jtphi_calo[1000];
  float hcalSum_calo[1000];
  float ecalSum_calo[1000];
  float refpt_calo[1000];
  float refeta_calo[1000];
  float refphi_calo[1000];
  float refdrjt_calo[1000];
  int subid_calo[1000];
  
  int ngen_calo;
  int genmatchindex_calo[500];
  int gensubid_calo[500];

  TTree *tr_calo=0;
  tr_calo = (TTree*)fin->Get(Form("%sCaloJetAnalyzer/t",kAlgName.c_str()));
  tr_calo->SetBranchAddress("nref",&nref_calo);
  tr_calo->SetBranchAddress("rawpt",rawpt_calo);
  tr_calo->SetBranchAddress("jtpt" ,jtpt_calo);
  tr_calo->SetBranchAddress("jteta",jteta_calo);
  tr_calo->SetBranchAddress("jtphi",jtphi_calo);
  tr_calo->SetBranchAddress("hcalSum",hcalSum_calo);
  tr_calo->SetBranchAddress("ecalSum",ecalSum_calo);
  tr_calo->SetBranchAddress("refpt",refpt_calo);
  tr_calo->SetBranchAddress("refphi",refphi_calo);
  tr_calo->SetBranchAddress("refeta",refeta_calo);
  tr_calo->SetBranchAddress("refdrjt",refdrjt_calo);
  tr_calo->SetBranchAddress("subid",subid_calo);

  tr_calo->SetBranchAddress("ngen",&ngen_calo);
  tr_calo->SetBranchAddress("genmatchindex",genmatchindex_calo);
  tr_calo->SetBranchAddress("gensubid",gensubid_calo);
  
  tr_calo->SetBranchStatus("*",0,0);
  tr_calo->SetBranchStatus("nref" ,1,0);
  tr_calo->SetBranchStatus("rawpt",1,0);
  tr_calo->SetBranchStatus("jtpt" ,1,0);
  tr_calo->SetBranchStatus("jteta",1,0);
  tr_calo->SetBranchStatus("jtphi",1,0);
  tr_calo->SetBranchStatus("hcalSum",1,0);
  tr_calo->SetBranchStatus("ecalSum",1,0);
  tr_calo->SetBranchStatus("refpt",1,0);
  tr_calo->SetBranchStatus("refphi",1,0);
  tr_calo->SetBranchStatus("refeta",1,0);
  tr_calo->SetBranchStatus("refdrjt",1,0);
  tr_calo->SetBranchStatus("subid",1,0);
  
  tr_calo->SetBranchStatus("ngen",1,0);
  tr_calo->SetBranchStatus("genmatchindex",1,0);
  tr_calo->SetBranchStatus("gensubid",1,0);



  TEventList* el = new TEventList("el","el");
  stringstream selection; selection<<"pthat<"<<maxpthat;
  tr_in->Draw(">>el",selection.str().c_str());
  double fentries = el->GetN();
  cout<<"tree entries  :  "<<kAlgName.c_str()<<" algorithm : " << nentries<<" elist: "<< fentries <<endl;
  delete el;

  //! weight  for the merging of the samples for different pT hat bins
  weight  = GetXsec(maxpthat);
  double wxs = weight/(fentries/100000.);
    
  //! Add Friends to the TTree
  tr_in->AddFriend(tr_ev);
  tr_in->AddFriend(tr_skim);
  tr_in->AddFriend(tr_calo);

  //! Vertex re-weighting
  TF1 *fVz=0;
  fVz = new TF1("fVz","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
  fVz->SetParameters(9.86748e-01, -8.91367e-03, 5.35416e-04, 2.67665e-06, -2.01867e-06);


  TH1F *hCentWeight = new TH1F("hCentWeight","Centrality weight",200,0,200);
  GetCentWeight(hCentWeight);  


  double nevt [ncen]={0};
  double njets[ncen]={0};
  double njets_calo[ncen]={0};

  Int_t iEvent=0;     
  for (Long64_t i=0; i<nentries;i++) {
    //for (Long64_t i=0; i<500;i++) {
    nbytes += tr_in->GetEntry(i);
    //return 0;
    
    double wcen= 1;//hCentWeight->GetBinContent(hCentWeight->FindBin(hiBin));
    double wvz = 1;//fVz->Eval(vz);

    if(pthat > maxpthat || !pcollisionEventSelection || fabs(vz) > 15.)continue;
    
    //std::cout<<" ********** Event # " <<i<<"\t weight  : "<<wxs<<"\t pthat : "<<pthat<<" nref_pf : "<< nref << " nref_calo : " << nref_calo <<  std::endl;

    hVz->Fill(vz);
    hBin->Fill(hiBin);

    //! xsec-weight
    hpthat->Fill(pthat,wxs*wcen*wvz);

    int iCent = GetCentBin(hiBin);
    if(iCent<0 || iCent>=ncen)continue;

    
    //if(i%50000==0)std::cout<<" ********** Event # " <<i<<"\t weight  : "<<wxs<<"\t pthat : "<<pthat<<std::endl;
    nevt[iCent]++;


    for(int irec=0; irec<nref; irec++){
      if( fabs(jteta[irec]) > ketacut || jtpt[irec] < 20 || refpt[irec] > 3.*pthat )continue;

      hsubid[iCent]->Fill(subid[irec],wxs*wcen*wvz);

      hPt[iCent][0]->Fill(jtpt[irec],wxs*wcen*wvz);
      hRawPt[iCent][0]->Fill(rawpt[irec],wxs*wcen*wvz);
      if(subid[irec]>=0){
	hPt[iCent][1]->Fill(jtpt[irec],wxs*wcen*wvz);
	hRawPt[iCent][1]->Fill(rawpt[irec],wxs*wcen*wvz);
      }
      if(subid[irec]==0){
	hPt[iCent][2]->Fill(jtpt[irec],wxs*wcen*wvz);
	hRawPt[iCent][2]->Fill(rawpt[irec],wxs*wcen*wvz);
      }
      if(subid[irec]> 0){
	hPt[iCent][3]->Fill(jtpt[irec],wxs*wcen*wvz);
	hRawPt[iCent][3]->Fill(rawpt[irec],wxs*wcen*wvz);
      }
      if(subid[irec]< 0){
	hPt[iCent][4]->Fill(jtpt[irec],wxs*wcen*wvz);
	hRawPt[iCent][4]->Fill(rawpt[irec],wxs*wcen*wvz);
      }
    }


    //! PF
    for(int irec=0; irec<nref; irec++){
      if( subid[irec] !=0 )continue;
      if(fabs(refeta[irec]) > ketacut || jtpt[irec] < kptrecocut)continue;
      if( wJetId ){
	if( (eMax[irec]/jtpt[irec])>=0.6 || (chargedMax[irec]/jtpt[irec])<=0.02 )continue;
      }

      hPtAll  [0][iCent]->Fill(refpt [irec],wxs*wcen*wvz);
      hEtaAll [0][iCent]->Fill(refeta[irec],wxs*wcen*wvz);
      hPhiAll [0][iCent]->Fill(refphi[irec],wxs*wcen*wvz);

      //! eff
      if( refdrjt[irec] < kdRCut ){
	hPtEff [0][iCent]->Fill(refpt [irec],wxs*wcen*wvz);
	hEtaEff[0][iCent]->Fill(refeta[irec],wxs*wcen*wvz);
	hPhiEff[0][iCent]->Fill(refphi[irec],wxs*wcen*wvz);
      }
    }

    //! Calo
    for(int irec=0; irec<nref_calo; irec++){
      if( subid_calo[irec] !=0 )continue;
      if(fabs(refeta_calo[irec]) > ketacut || jtpt_calo[irec] < kptrecocut)continue;

      hPtAll  [1][iCent]->Fill(refpt_calo [irec],wxs*wcen*wvz);
      hEtaAll [1][iCent]->Fill(refeta_calo[irec],wxs*wcen*wvz);
      hPhiAll [1][iCent]->Fill(refphi_calo[irec],wxs*wcen*wvz);

      //! eff
      if( refdrjt_calo[irec] < kdRCut ){
	hPtEff [1][iCent]->Fill(refpt_calo[irec],wxs*wcen*wvz);
	hEtaEff[1][iCent]->Fill(refeta_calo[irec],wxs*wcen*wvz);
	hPhiEff[1][iCent]->Fill(refphi_calo[irec],wxs*wcen*wvz);
      }
    }



    //! Fake
    for(int irec=0; irec<nref; irec++){
      if(fabs(jteta[irec]) > ketacut)continue;
      if(jtpt[irec] < 20 /*|| jtpt[irec] > 3.*pthat*/ ) continue;
      if( wJetId ){
	if( (eMax[irec]/jtpt[irec])>=0.6 || (chargedMax[irec]/jtpt[irec])<=0.02 )continue;
      }
      hPtFakeAll  [0][iCent]->Fill(jtpt [irec]/*,wxs*wcen*wvz*/);
      if( jtpt[irec] > 20) {
	hEtaFakeAll_20 [0][iCent]->Fill(jteta[irec]/*,wxs*wcen*wvz*/);
	hPhiFakeAll_20 [0][iCent]->Fill(jtphi[irec]/*,wxs*wcen*wvz*/);
      }
      if( jtpt[irec] > 30) {
	hEtaFakeAll_30 [0][iCent]->Fill(jteta[irec]/*,wxs*wcen*wvz*/);
	hPhiFakeAll_30 [0][iCent]->Fill(jtphi[irec]/*,wxs*wcen*wvz*/);
      }
      if( jtpt[irec] > 40) {
	hEtaFakeAll_40 [0][iCent]->Fill(jteta[irec]/*,wxs*wcen*wvz*/);
	hPhiFakeAll_40 [0][iCent]->Fill(jtphi[irec]/*,wxs*wcen*wvz*/);
      }
      if( jtpt[irec] > 45) {
	hEtaFakeAll_45 [0][iCent]->Fill(jteta[irec]/*,wxs*wcen*wvz*/);
	hPhiFakeAll_45 [0][iCent]->Fill(jtphi[irec]/*,wxs*wcen*wvz*/);
      }
      if( jtpt[irec] > 50) {
	hEtaFakeAll_50 [0][iCent]->Fill(jteta[irec]/*,wxs*wcen*wvz*/);
	hPhiFakeAll_50 [0][iCent]->Fill(jtphi[irec]/*,wxs*wcen*wvz*/);
      }


      //! Fake
      if(refpt[irec] < 0){
	//cout <<"  Fake : " << jtpt[irec] << endl;
	hPtFake [0][iCent]->Fill(jtpt [irec]/*,wxs*wcen*wvz*/);
	if(jtpt[irec] > 20){
	  hEtaFake_20[0][iCent]->Fill(jteta[irec]/*,wxs*wcen*wvz*/);
	  hPhiFake_20[0][iCent]->Fill(jtphi[irec]/*,wxs*wcen*wvz*/);
	}
	if(jtpt[irec] > 30){
	  hEtaFake_30[0][iCent]->Fill(jteta[irec]/*,wxs*wcen*wvz*/);
	  hPhiFake_30[0][iCent]->Fill(jtphi[irec]/*,wxs*wcen*wvz*/);
	}
	if(jtpt[irec] > 40){
	  hEtaFake_40[0][iCent]->Fill(jteta[irec]/*,wxs*wcen*wvz*/);
	  hPhiFake_40[0][iCent]->Fill(jtphi[irec]/*,wxs*wcen*wvz*/);
	}
	if(jtpt[irec] > 45){
	  hEtaFake_45[0][iCent]->Fill(jteta[irec]/*,wxs*wcen*wvz*/);
	  hPhiFake_45[0][iCent]->Fill(jtphi[irec]/*,wxs*wcen*wvz*/);
	}
	if(jtpt[irec] > 50){
	  hEtaFake_50[0][iCent]->Fill(jteta[irec]/*,wxs*wcen*wvz*/);
	  hPhiFake_50[0][iCent]->Fill(jtphi[irec]/*,wxs*wcen*wvz*/);
	}
      }
    }

    //!
    //! Fake Calo
    for(int irec=0; irec<nref_calo; irec++){
      if(fabs(jteta_calo[irec]) > ketacut)continue;
      if(jtpt_calo[irec] < 20 /*|| jtpt[irec] > 3.*pthat*/ ) continue;

      hPtFakeAll  [1][iCent]->Fill(jtpt_calo[irec]/*,wxs*wcen*wvz*/);
      if( jtpt_calo[irec] > 20) {
	hEtaFakeAll_20 [1][iCent]->Fill(jteta_calo[irec]/*,wxs*wcen*wvz*/);
	hPhiFakeAll_20 [1][iCent]->Fill(jtphi_calo[irec]/*,wxs*wcen*wvz*/);
      }
      if( jtpt_calo[irec] > 30) {
	hEtaFakeAll_30 [1][iCent]->Fill(jteta_calo[irec]/*,wxs*wcen*wvz*/);
	hPhiFakeAll_30 [1][iCent]->Fill(jtphi_calo[irec]/*,wxs*wcen*wvz*/);
      }
      if( jtpt_calo[irec] > 40) {
	hEtaFakeAll_40 [1][iCent]->Fill(jteta_calo[irec]/*,wxs*wcen*wvz*/);
	hPhiFakeAll_40 [1][iCent]->Fill(jtphi_calo[irec]/*,wxs*wcen*wvz*/);
      }
      if( jtpt_calo[irec] > 45) {
	hEtaFakeAll_45 [1][iCent]->Fill(jteta_calo[irec]/*,wxs*wcen*wvz*/);
	hPhiFakeAll_45 [1][iCent]->Fill(jtphi_calo[irec]/*,wxs*wcen*wvz*/);
      }
      if( jtpt_calo[irec] > 50) {
	hEtaFakeAll_50 [1][iCent]->Fill(jteta_calo[irec]/*,wxs*wcen*wvz*/);
	hPhiFakeAll_50 [1][iCent]->Fill(jtphi_calo[irec]/*,wxs*wcen*wvz*/);
      }
      //! Fake
      if(refpt_calo[irec] < 0){
	//cout <<"  Fake : " << jtpt[irec] << endl;
	hPtFake [1][iCent]->Fill(jtpt_calo[irec]/*,wxs*wcen*wvz*/);
	if(jtpt_calo[irec] > 20){
	  hEtaFake_20[1][iCent]->Fill(jteta_calo[irec]/*,wxs*wcen*wvz*/);
	  hPhiFake_20[1][iCent]->Fill(jtphi_calo[irec]/*,wxs*wcen*wvz*/);
	}
	if(jtpt_calo[irec] > 30){
	  hEtaFake_30[1][iCent]->Fill(jteta_calo[irec]/*,wxs*wcen*wvz*/);
	  hPhiFake_30[1][iCent]->Fill(jtphi_calo[irec]/*,wxs*wcen*wvz*/);
	}
	if(jtpt_calo[irec] > 40){
	  hEtaFake_40[1][iCent]->Fill(jteta_calo[irec]/*,wxs*wcen*wvz*/);
	  hPhiFake_40[1][iCent]->Fill(jtphi_calo[irec]/*,wxs*wcen*wvz*/);
	}
	if(jtpt_calo[irec] > 45){
	  hEtaFake_45[1][iCent]->Fill(jteta_calo[irec]/*,wxs*wcen*wvz*/);
	  hPhiFake_45[1][iCent]->Fill(jtphi_calo[irec]/*,wxs*wcen*wvz*/);
	}
	if(jtpt_calo[irec] > 50){
	  hEtaFake_50[1][iCent]->Fill(jteta_calo[irec]/*,wxs*wcen*wvz*/);
	  hPhiFake_50[1][iCent]->Fill(jtphi_calo[irec]/*,wxs*wcen*wvz*/);
	}
      }
    }


    std::vector < Jet > pfjet_coll;

    //! Gen matched jets loop PF
    for(int igen=0; igen<ngen; igen++){
      if( gensubid[igen] !=0 ) continue;
      int gj = genmatchindex[igen];

      //for(int igen=0; igen<nref; igen++){
      //int gj=igen;
      
      if( rawpt[gj] < kptrawcut  ) continue;
      if( jtpt[gj]  < kptrecocut ) continue;
      if( refpt[gj] < kptgencut  ) continue;
      if( fabs(refeta[gj]) > ketacut ) continue;
      if( refdrjt[gj] > kdRCut )continue;
      if( fabs(refpt[gj]) > 2.*pthat )continue;

      if( wJetId ){
	if( (eMax[gj]/jtpt[gj])>=0.6 || (chargedMax[gj]/jtpt[gj])<=0.02 )continue;
      }

      float recopt  = jtpt[gj];
      float recoeta = jteta[gj];
      float recophi = jtphi[gj];
      //float delr    = refdrjt[gj];
      //double ratio  = recopt/refpt[gj];

      double rchmax = chargedMax[gj] / recopt;
      double rnemax = neutralMax[gj] / recopt;
      double rphmax = photonMax[gj]  / recopt;
      double rmumax = muMax[gj]      / recopt;
      double relmax = eMax[gj]       / recopt;

      double resp_corr =  recopt    / refpt[gj];
      double resp_raw  =  rawpt[gj] / refpt[gj];
      
      int iphi = GetPhiBin(refphi[gj]);
      int ieta = GetEtaBin(refeta[gj]);
      int ipt  = GetPtBin (refpt[gj]);
      int ipt_wide = GetPtBinWide (refpt[gj]);

      double sumPF     = (chargedSum[gj] + neutralSum[gj] + photonSum[gj] + eMax[gj] + muMax[gj]);
      double bkgd      = sumPF - rawpt[gj];

      hchargedmax[iCent]->Fill(rchmax,wxs*wcen*wvz);
      hneutralmax[iCent]->Fill(rnemax,wxs*wcen*wvz);
      hphotonmax [iCent]->Fill(rphmax,wxs*wcen*wvz);
      hmuonmax   [iCent]->Fill(rmumax,wxs*wcen*wvz);
      helecmax   [iCent]->Fill(relmax,wxs*wcen*wvz);
      
      njets[iCent]++;

      hrawpt_genpt [0][iCent]->Fill(refpt[gj],rawpt[gj],wxs*wcen*wvz);
      hrecopt_genpt[0][iCent]->Fill(refpt[gj],recopt,wxs*wcen*wvz);

      hgenpt_genm [0][iCent]->Fill(refpt[gj],wxs*wcen*wvz);
      hrecopt_genm[0][iCent]->Fill(recopt,wxs*wcen*wvz);
      hrawpt_genm [0][iCent]->Fill(rawpt[gj],wxs*wcen*wvz);
      
      hjeteta[0][iCent]->Fill(recoeta,wxs*wcen*wvz);
      hjetphi[0][iCent]->Fill(recophi,wxs*wcen*wvz);
      
      hgeneta[0][iCent]->Fill(refeta[gj],wxs*wcen*wvz);
      hgenphi[0][iCent]->Fill(refphi[gj],wxs*wcen*wvz);

      hgenrecpt_genm[0][iCent]->Fill(refpt[gj], recopt, wxs*wcen*wvz);
      hgenrawpt_genm[0][iCent]->Fill(refpt[gj], rawpt[gj], wxs*wcen*wvz);


      //! Response in pt and eta
      if( ipt>=0 && ipt<nbins ){

	Jet pfjet;
	pfjet.pt  = recopt;
	pfjet.eta = recoeta;
	pfjet.phi = recophi;
	pfjet.id  = gj;

	pfjet_coll.push_back(pfjet);

	havrefpt_genm [0][iCent][ipt]->Fill(refpt[gj],wxs*wcen*wvz);
	havjetpt_genm [0][iCent][ipt]->Fill(recopt,wxs*wcen*wvz);
	havrawpt_genm [0][iCent][ipt]->Fill(rawpt[gj],wxs*wcen*wvz);

	havbkgpt_genm [0][iCent][ipt]->Fill(bkgd,wxs*wcen*wvz);
	if(fabs(recoeta) < 1.3)havbkgpt_eta_genm [0][iCent][ipt][0]->Fill(bkgd,wxs*wcen*wvz);
	else havbkgpt_eta_genm [0][iCent][ipt][1]->Fill(bkgd,wxs*wcen*wvz);

	hrescrpt_genm [0][iCent][ipt]->Fill(resp_corr,wxs*wcen*wvz);
	hresrrpt_genm [0][iCent][ipt]->Fill(resp_raw ,wxs*wcen*wvz);
	
	if( ieta>=0 && ieta<neta ){
	  hrescrpt_genm_eta [0][iCent][ipt][ieta]->Fill(resp_corr,wxs*wcen*wvz);
	  hresrrpt_genm_eta [0][iCent][ipt][ieta]->Fill(resp_raw ,wxs*wcen*wvz);
	  if( ipt_wide>=0 ){
	    hrescrpt_wide_genm_eta [0][iCent][ipt_wide][ieta]->Fill(resp_corr,wxs*wcen*wvz);
	    hresrrpt_wide_genm_eta [0][iCent][ipt_wide][ieta]->Fill(resp_raw ,wxs*wcen*wvz);
	  }
	}
	
	if( iphi>=0 && iphi<nphi ){
	  hrescrpt_genm_phi [0][iCent][ipt][iphi]->Fill(resp_corr,wxs*wcen*wvz);
	  hresrrpt_genm_phi [0][iCent][ipt][iphi]->Fill(resp_raw ,wxs*wcen*wvz);
	  if( ipt_wide>=0 ){
	    hrescrpt_wide_genm_phi [0][iCent][ipt_wide][iphi]->Fill(resp_corr,wxs*wcen*wvz);
	    hresrrpt_wide_genm_phi [0][iCent][ipt_wide][iphi]->Fill(resp_raw ,wxs*wcen*wvz);
	  }
	}
      }
    }//! igen loop

    std::vector < Jet > calojet_coll;
    //! Calo jet loop
    //! Gen matched jets loop
    for(int igen=0; igen<ngen_calo; igen++){
      if( gensubid_calo[igen] !=0 ) continue;
      int gj = genmatchindex_calo[igen];

      //for(int igen=0; igen<nref; igen++){
      //int gj=igen;
      
      if( rawpt_calo[gj] < kptrawcut  ) continue;
      if( jtpt_calo[gj]  < kptrecocut ) continue;
      if( refpt_calo[gj] < kptgencut  ) continue;
      if( fabs(refeta_calo[gj]) > ketacut ) continue;
      if( refdrjt_calo[gj] > kdRCut )continue;
      if( fabs(refpt_calo[gj]) > 2.*pthat )continue;

      float recopt  = jtpt_calo[gj];
      float recoeta = jteta_calo[gj];
      float recophi = jtphi_calo[gj];
      //float delr    = refdrjt_calo[gj];
      //double ratio  = recopt/refpt_calo[gj];

      double resp_corr =  recopt    / refpt_calo[gj];
      double resp_raw  =  rawpt_calo[gj] / refpt_calo[gj];
      
      int iphi = GetPhiBin(refphi_calo[gj]);
      int ieta = GetEtaBin(refeta_calo[gj]);
      int ipt  = GetPtBin (refpt_calo[gj]);
      int ipt_wide = GetPtBinWide (refpt_calo[gj]);

      double sumCalo   = (hcalSum_calo[gj] + ecalSum_calo[gj]);
      double bkgd      = sumCalo - rawpt_calo[gj];
      
      njets_calo[iCent]++;

      hrawpt_genpt [1][iCent]->Fill(refpt_calo[gj],rawpt_calo[gj],wxs*wcen*wvz);
      hrecopt_genpt[1][iCent]->Fill(refpt_calo[gj],recopt,wxs*wcen*wvz);
      
      hgenpt_genm [1][iCent]->Fill(refpt_calo[gj],wxs*wcen*wvz);
      hrecopt_genm[1][iCent]->Fill(recopt,wxs*wcen*wvz);
      hrawpt_genm [1][iCent]->Fill(rawpt_calo[gj],wxs*wcen*wvz);
      
      hjeteta[1][iCent]->Fill(recoeta,wxs*wcen*wvz);
      hjetphi[1][iCent]->Fill(recophi,wxs*wcen*wvz);
      
      hgeneta[1][iCent]->Fill(refeta_calo[gj],wxs*wcen*wvz);
      hgenphi[1][iCent]->Fill(refphi_calo[gj],wxs*wcen*wvz);

      hgenrecpt_genm[1][iCent]->Fill(refpt_calo[gj], recopt, wxs*wcen*wvz);
      hgenrawpt_genm[1][iCent]->Fill(refpt_calo[gj], rawpt_calo[gj], wxs*wcen*wvz);


      //! Response in pt and eta
      if( ipt>=0 && ipt<nbins ){

	Jet calojet;
	calojet.pt  = recopt;
	calojet.eta = recoeta;
	calojet.phi = recophi;
	calojet.id  = gj;

	calojet_coll.push_back(calojet);

	havrefpt_genm [1][iCent][ipt]->Fill(refpt_calo[gj],wxs*wcen*wvz);
	havjetpt_genm [1][iCent][ipt]->Fill(recopt,wxs*wcen*wvz);
	havrawpt_genm [1][iCent][ipt]->Fill(rawpt_calo[gj],wxs*wcen*wvz);

	havbkgpt_genm [1][iCent][ipt]->Fill(bkgd,wxs*wcen*wvz);
	if(fabs(recoeta) < 1.3)havbkgpt_eta_genm [1][iCent][ipt][0]->Fill(bkgd,wxs*wcen*wvz);
	else havbkgpt_eta_genm [1][iCent][ipt][1]->Fill(bkgd,wxs*wcen*wvz);

	hrescrpt_genm [1][iCent][ipt]->Fill(resp_corr,wxs*wcen*wvz);
	hresrrpt_genm [1][iCent][ipt]->Fill(resp_raw ,wxs*wcen*wvz);
	
	if( ieta>=0 && ieta<neta ){
	  hrescrpt_genm_eta [1][iCent][ipt][ieta]->Fill(resp_corr,wxs*wcen*wvz);
	  hresrrpt_genm_eta [1][iCent][ipt][ieta]->Fill(resp_raw ,wxs*wcen*wvz);
	  if( ipt_wide>=0 ){
	    hrescrpt_wide_genm_eta [1][iCent][ipt_wide][ieta]->Fill(resp_corr,wxs*wcen*wvz);
	    hresrrpt_wide_genm_eta [1][iCent][ipt_wide][ieta]->Fill(resp_raw ,wxs*wcen*wvz);
	  }
	}
	
	if( iphi>=0 && iphi<nphi ){
	  hrescrpt_genm_phi [1][iCent][ipt][iphi]->Fill(resp_corr,wxs*wcen*wvz);
	  hresrrpt_genm_phi [1][iCent][ipt][iphi]->Fill(resp_raw ,wxs*wcen*wvz);
	  if( ipt_wide>=0 ){
	    hrescrpt_wide_genm_phi [1][iCent][ipt_wide][iphi]->Fill(resp_corr,wxs*wcen*wvz);
	    hresrrpt_wide_genm_phi [1][iCent][ipt_wide][iphi]->Fill(resp_raw ,wxs*wcen*wvz);
	  }
	}
      }
    }//! igen loop

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

    //! Find data driven JEC Calo jet
    if(calojet_coll.size()>1){
      std::sort(calojet_coll.begin(), calojet_coll.end(), compare_pt);
      ijet = calojet_coll.begin();
      Jet ljet  =  *(ijet);
      Jet sljet =  *(ijet + 1);
      double dphi=-9999;
      if ( ljet.pt > 60. && sljet.pt > 60. ){
	dphi = delphi(ljet.phi, sljet.phi);
	if ( dphi > 2.*pi/3. ){ //! matched dijet
	  double ptdij = (ljet.pt + sljet.pt)/2.;  
	  int mstat=1;
	  if (calojet_coll.size()>2){
	    Jet thjet =  *(ijet + 2);
	    if(thjet.pt/ptdij > 0.2)mstat=0;	
	  }
	  if(mstat){
	    double B=-9999;
	    double rn1 = gRandom->Rndm();
	    double rn2 = gRandom->Rndm();
	    if(rn1 > rn2){
	      B = (ljet.pt - sljet.pt)/(ljet.pt + sljet.pt);
	    }else{
	      B = (sljet.pt - ljet.pt)/(sljet.pt + ljet.pt);
	    }
	    hrescrpt_datadr[1][iCent]->Fill(ptdij, B, wxs*wcen*wvz);
	  }
	}
      }
    }

    iEvent++;
    //std::cout<<"Completed event #  "<<ievt<<std::endl; 
  }//! event loop ends


  std::cout<<std::endl;
  std::cout<<kAlgName.c_str() << std::endl;
  for(int ic=0;ic<ncen;ic++){
    std::cout<<"\t cent : " << ccent[ic] << "  # of events : "<< " PbPb : " << nevt[ic] << " PFJets  " << njets[ic] << " CaloJets : " << njets_calo[ic] << std::endl;
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
int GetPtBinWide(float pt)
{
  if(pt<ptbins_wide[0] || pt>ptbins_wide[nbins_wide])return -1;
  for(int ix=0;ix<nbins_wide;ix++){
    if(pt>=ptbins_wide[ix] && pt<ptbins_wide[ix+1]){
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
double deltaR(float eta1, float phi1, float eta2, float phi2)
{
  float deta = eta1 - eta2;
  float dphi = fabs(phi1 - phi2);
  if(dphi > pi)dphi -= 2*pi;
  double dr = sqrt(pow(deta,2) + pow(dphi,2));
  return dr;
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
