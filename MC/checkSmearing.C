#include "hiForest.h"
#include "SmearingFactors.h"
#include <TF1.h>
#include <TProfile.h>
#include <TChain.h>
#include <TTree.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TFile.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <cstdlib>
#include <cmath>
#include <iostream>
using namespace std;


#ifdef _MAKECINT_
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ class HiForest;
#pragma link C++ class Evts;
#pragma link C++ class Jets;
#endif


int GetCentBin(int /*bin*/);
int GetEtaBin(float /*eta*/);
int GetDetEtaBin(float /*eta*/);
int GetPhiBin(float /*phi*/);
int GetPtBin(float /*pt*/);
bool selecJet(Jets */*jetc*/, int /*ij*/);
void FindLeadSubLeadJets(Jets */*jetc*/, int */*ljet*/);

//! pt binning
//double ptbins[] ={30,40,50,60,70,80,90, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400};
double ptbins[] ={30,50,70,90,120,160,200,240,280,340,400};
const int bins = sizeof(ptbins)/sizeof(double) - 1;

//! constants
const int iYear=2012;
const double pi=acos(-1.);
const double pi2=2*pi -1;
const int ncen=7;
const int ketar=4;
const int maxe=5;
const int maxph=5;
const float ketacut=2.;
const double kptrecocut=15.;
const double kptgencut =0.;
const double kdRcut=0.3;
const double kVzcut=15.;

const int rbins=50;
double rbinl=0.0;
double rbinh=5.00;
/* Full list 
const int knj=25;
const char *calgo[knj]= {
  "icPu5Calo",
  "ak1PF","ak2PF","ak3PF","ak4PF","ak5PF","ak6PF",
  "ak1Calo","ak2Calo","ak3Calo","ak4Calo","ak5Calo","ak6Calo",
  "akPu1PF","akPu2PF","akPu3PF","akPu4PF","akPu5PF","akPu6PF",
  "akPu1Calo","akPu2Calo","akPu3Calo","akPu4Calo","akPu5Calo","akPu6Calo"
};
*/


const int knj = 7;
//const char *cjets[knj] = {"icPu5","ak2PF","ak3PF","ak4PF","akPu2PF","akPu3PF","akPu4PF"};

  

TStopwatch timer;
int checkSmearing(double kPt=80,const char *ksp="pp")
{

  timer.Start();

  //! Load Lib
  gSystem->Load("hiForest_h.so");
  
  //! Load the hiforest input root file
  HiForest *c = new HiForest(Form("/net/hisrv0001/home/icali/hadoop/Pythia/Z2/ppDijet_merged/pp276Dijet%0.0f_merged.root",kPt),"pp2012",1,1);

  //! Load Array for smearing
  LoadParameters();

  double xsection=0;
  double xup=0;
  double xsub=0;
  double maxpthat=9999;
  if(kPt==15){
    maxpthat=30;
    xup=1.566e-01;
    xsub=1.079e-02;
  }
  else if(kPt==30){
    maxpthat=50;
    xup=1.079e-02;
    xsub=1.021e-03;
  }
  else if(kPt==50){
    maxpthat=80;
    xup=1.021e-03;
    xsub=9.913e-05;
  }
  else if(kPt==80){
    maxpthat=120;
    xup=9.913e-05;
    xsub=1.128e-05;
  }
  else if(kPt==120){
    maxpthat=170;
    xup=1.128e-05;
    xsub=1.470e-06;
  }
  else if(kPt==170){
    maxpthat=200;
    xup=1.470e-06;
    xsub=5.310e-07;
  }
  else if(kPt==200){
    maxpthat=250;
    xup=5.310e-07;
    xsub=1.192e-07;
  }
  else if(kPt==250){
    maxpthat=300;
    xup=1.192e-07;
    xsub=3.176e-08;
  }
  else if(kPt==300){
    maxpthat=9999;
    xup=3.176e-08;
    xsub=0;
  }
  xsection = xup-xsub;

  std::cout<<std::endl;
  std::cout<<std::endl;


  //! Don't want to loop over trees which is not used in the analysis
  //! event trees
  //c->hasEvtTree=0;

  //! jet trees
  //! Switch on only the jet trees which you  require
  const char *cjets[knj] = {"icPu5","ak2PF","ak3PF","ak4PF","akPu2PF","akPu3PF","akPu4PF"};
  c->SelectJetAlgo(cjets,knj);
  

  //! photon tree
  //c->hasPhotonTree=0;


  //! Select only the jet branches which you are going to use
  //! This increases the speed for running over trees with lot of branches.
  //! This is currently for only jet algos and evtTree  uniformly applied to all the 

  //! Event Tree 
  const char *evlist[]={"hiBin","vz","hiHF"};
  const int kevbr = sizeof(evlist)/sizeof(const char *);
  c->SelectBranches("evtTree",evlist,kevbr);

  //! jet Tree algorithms
  const char *jtlist[]={"nref","pthat","rawpt","jtpt","jteta","jtphi","jtpu","refpt","refeta","refphi","refdrjt","refparton_flavor",
			"ngen","gensubid","genmatchindex"};

  const int kjtbr = sizeof(jtlist)/sizeof(const char *);
  c->SelectBranches("JetTree",jtlist,kjtbr);

  std::cout<<"Selected the branches of need from evtTree and JetTrees : "<<std::endl;

  //! To get the jet object from hiforest
  Jets *iJet=0;
  //const int knj = 7;//c->GetNAlgo(); //! # of jet algorithms in this hiforest 
  std::cout<<"Loaded all tree variables and # of jet algorithms : "<<knj<<std::endl;
  std::cout<<"\t"<<std::endl;


  //! Open a output file for histos
  TFile *fout = new TFile(Form("check/test_Smear_%0.0fGeV_%s_%d.root",kPt,ksp,iYear),"RECREATE");

  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"**************************************************** "<<std::endl;
  std::cout<<Form("Running for %s ",ksp)<<std::endl;
  std::cout<<Form("pT  cut for %0.3f ",kptrecocut)<<std::endl;
  std::cout<<Form("eta cut for %0.3f ",ketacut)<<std::endl;
  std::cout<<"My hiForest TTree : " <<c->GetName()<<std::endl;
  std::cout<<"Output file  "<<fout->GetName()<<std::endl;
  std::cout<<"**************************************************** "<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;


  //! 
  //! Define histograms here
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TProfile::SetDefaultSumw2();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  TH1F *hEvt    = new TH1F("hEvt","# of events ",4,0,4);
  TH1F *hVz     = new TH1F("hVz","# of events ",80,-20,20);
  TH1F *hBin    = new TH1F("hBin","Centrality bin",40,0,40);
  TH1F *hHF     = new TH1F("hHF","Centrality variable from HF",600,0,6000);
  TH1F *hTotEve = new TH1F("hTotEve","# of events in the skimmed files",4,0,4);

  TH1F *hgenpt [knj][ncen], *hrecopt[knj][ncen], *hrawpt[knj][ncen], *hsmpt[knj][ncen];
  TH1F *hgenptC [knj][ncen], *hrecoptC[knj][ncen], *hrawptC[knj][ncen], *hsmptC[knj][ncen];
  TH1F *hgeneta[knj][ncen], *hrecoeta[knj][ncen];
  TH1F *hgenphi[knj][ncen], *hrecophi[knj][ncen];

  //! Ratios of the pt distributions
  TProfile *hrecogen[knj][ncen], *hrecoraw[knj][ncen], *hrawgen[knj][ncen], *hrecoraw_ref[knj][ncen];

  //! Resposnse
  TH2F *hcorrptrefpt[knj][ncen], *hrawptrefpt[knj][ncen], *hcorrptrawpt[knj][ncen];
  TH2F *hrescrpt[knj][ncen], *hresrrpt[knj][ncen], *hresrcrpt[knj][ncen], *hresorpt[knj][ncen];
  TH2F *hratiorawrefpt[knj][ncen], *hratiocorrrefpt[knj][ncen], *hratiocorrrawpt[knj][ncen];
  TH2F *hratiorawrefpt_eta[knj][ncen][2], *hratiocorrrefpt_eta[knj][ncen][2];

  TH2F *hratiocorrrefpt_genm[knj][ncen];
  TH2F *hratiocorrrefpt_lead[knj][ncen], *hratiocorrrefpt_slead[knj][ncen], *hratiocorrrefpt_remain[knj][ncen];
  TH2F *hratiocorrrefpt_pthat[knj][ncen];
  TH2F *hratiocorrorgpt[knj][ncen];

  TH2F *hpteta[knj][ncen][maxe], *hptphi[knj][ncen][maxph];

  TH2F *hgenjrecoj[knj][ncen];
  TH2F *hgenjrecoj_pflavor[knj][ncen];
  
  TH2F *hFracjets[knj][ncen];
  TH2F *hAvjets[knj][ncen];
  TH2F *hMeanPtjets[knj][ncen];
  TH1F *hNjets[knj][ncen];
  TH1F *hNevt [knj][ncen];

  //! Pileup effect study
  TH2F *hjetptpu[knj][ncen];             //! centrality                                                                       
  TH2F *hjetptpu_etab[knj][ncen][ketar]; //! eta dependence             

  //! Efficency histos
  TH1F *hPtAll [knj][ncen], *hPtSel[knj][ncen];
  TH1F *hEtaAll[knj][ncen], *hEtaSel[knj][ncen];
  TH1F *hPhiAll[knj][ncen], *hPhiSel[knj][ncen];

  //! Response vs deltar
  TH1F *hRspVsDeltaR[knj][ncen][25];

  //! DeltaR efficiency                                                                                                                                                              
  TH1F *hDeltaR[knj][ncen];
  TH1F *hDeltaRAll[knj][ncen];
  TH1F *hDeltaRSel[knj][ncen];

  for(int nj=0;nj<knj;nj++){
    //const char *algoname = c->GetAlgoName(nj);
    //char *algoname = cjets[nj];
    //strcpy(algoname,cjets[nj]);

    for(int icen=0;icen<ncen;icen++){
      hFracjets[nj][icen]   = new TH2F(Form("hFracjets%d_%d",nj,icen),Form("Fraction of jets in given pt hat cent %d %s",icen,cjets[nj]),500,0,1000,500,0.,1000.);
      hAvjets[nj][icen]     = new TH2F(Form("hAvjets%d_%d",nj,icen),Form("<#> of jets cent %d %s",icen,cjets[nj]),500,0,1000,30,0,30);
      hMeanPtjets[nj][icen] = new TH2F(Form("hMeanPtjets%d_%d",nj,icen),Form("<pT> of jets cent %d %s",icen,cjets[nj]),500,0,1000,500,0,1000);

      hNjets[nj][icen] = new TH1F(Form("hNjets%d_%d",nj,icen),Form("# of jets cent %d jets %s",icen,cjets[nj]),3,0.0,3.0);
      hNevt [nj][icen] = new TH1F(Form("hNevt%d_%d",nj,icen),Form("# of events cent %d %s",icen,cjets[nj]),40,-40,40);
      hgeneta[nj][icen] = new TH1F(Form("hgeneta%d_%d",nj,icen),Form("gen eta distribution jet centb %d %s",icen,cjets[nj]),60,-3.0,3.0);
      hrecoeta[nj][icen] = new TH1F(Form("hrecoeta%d_%d",nj,icen),Form("reco eta distribution jet centb %d %s",icen,cjets[nj]),60,-3.0,3.0);

      hgenphi[nj][icen] = new TH1F(Form("hgenphi%d_%d",nj,icen),Form("gen phi distribution jet centb %d %s",icen,cjets[nj]),36,-pi,pi);
      hrecophi[nj][icen] = new TH1F(Form("hrecophi%d_%d",nj,icen),Form("reco phi distribution jet centb %d %s",icen,cjets[nj]),36,-pi,pi);

      hgenpt[nj][icen]  = new TH1F(Form("hgenpt%d_%d",nj,icen),Form("gen p_{T} distribution jet centb %d %s",icen,cjets[nj]),500,0,1000);
      hrecopt[nj][icen] = new TH1F(Form("hrecopt%d_%d",nj,icen),Form("reco p_{T} distribution jet centb %d %s",icen,cjets[nj]),500,0,1000);
      hrawpt[nj][icen]  = new TH1F(Form("hrawpt%d_%d",nj,icen),Form("raw p_{T} distribution jet centb %d %s",icen,cjets[nj]),500,0,1000);
      hsmpt[nj][icen] = new TH1F(Form("hsmpt%d_%d",nj,icen),Form("smeared p_{T} distribution jet centb %d %s",icen,cjets[nj]),500,0,1000);


      //! with pt bins
      hgenptC[nj][icen]  = new TH1F(Form("hgenptC%d_%d",nj,icen),Form("gen p_{T} distribution jet centb %d %s",icen,cjets[nj]),bins,ptbins);
      hrecoptC[nj][icen] = new TH1F(Form("hrecoptC%d_%d",nj,icen),Form("reco p_{T} distribution jet centb %d %s",icen,cjets[nj]),bins,ptbins);
      hrawptC[nj][icen]  = new TH1F(Form("hrawptC%d_%d",nj,icen),Form("raw p_{T} distribution jet centb %d %s",icen,cjets[nj]),bins,ptbins);
      hsmptC[nj][icen] = new TH1F(Form("hsmptC%d_%d",nj,icen),Form("smeared p_{T} distribution jet centb %d %s",icen,cjets[nj]),bins,ptbins);

      //! Ratios
      hrecogen[nj][icen] = new TProfile(Form("hrecogen%d_%d",nj,icen),Form("reco/gen p_{T} distribution jet centb %d %s",icen,cjets[nj]),500,0,1000);
      hrecoraw[nj][icen] = new TProfile(Form("hrecoraw%d_%d",nj,icen),Form("reco/raw p_{T} distribution jet centb %d %s",icen,cjets[nj]),500,0,1000);
      hrecoraw_ref[nj][icen]  = new TProfile(Form("hrecoraw_ref%d_%d",nj,icen),Form("reco/raw : ref p_{T} distribution jet centb %d %s",icen,cjets[nj]),500,0,1000);
      hrawgen[nj][icen]  = new TProfile(Form("hrawgen%d_%d",nj,icen),Form("raw/gen p_{T} distribution jet centb %d %s",icen,cjets[nj]),500,0,1000);


      hcorrptrefpt[nj][icen]= new TH2F(Form("hcorrptrefpt%d_%d",nj,icen),Form("Gen jet:Reco jet p_{T} distribution jet centb %d %s",icen,cjets[nj]),500,0,1000,500,0,1000);
      hrawptrefpt[nj][icen]= new TH2F(Form("hrawptrefpt%d_%d",nj,icen),Form("Gen jet:Raw jet p_{T}  distribution jet centb %d %s",icen,cjets[nj]),500,0,1000,500,0,1000);
      hcorrptrawpt[nj][icen]= new TH2F(Form("hcorrptrawpt%d_%d",nj,icen),Form("Reco jet Corr:Raw jet p_{T}  distribution jet centb %d %s",icen,cjets[nj]),500,0,1000,500,0,1000);

      hrescrpt[nj][icen]= new TH2F(Form("hrescrpt%d_%d",nj,icen),Form("Gen jet:(Reco/Gen) jet p_{T} distribution jet centb %d %s",icen,cjets[nj]),500,0,1000,150,rbinl,rbinh);
      hresrrpt[nj][icen]= new TH2F(Form("hresrrpt%d_%d",nj,icen),Form("Gen jet:(Raw/Gen) jet p_{T}  distribution jet centb %d %s",icen,cjets[nj]),500,0,1000,150,rbinl,rbinh);
      hresrcrpt[nj][icen]= new TH2F(Form("hresrcrpt%d_%d",nj,icen),Form("Reco jet:(Reco/Raw) jet p_{T} distribution jet centb %d %s",icen,cjets[nj]),500,0,1000,150,rbinl,rbinh);
      hresorpt[nj][icen]= new TH2F(Form("hresorpt%d_%d",nj,icen),Form("Gen jet:(Smeared/Reco) jet p_{T} distribution jet centb %d %s",icen,cjets[nj]),500,0,1000,150,rbinl,rbinh);

      hratiorawrefpt[nj][icen]= new TH2F(Form("hratiorawrefpt%d_%d",nj,icen),Form("Raw jet / Gen jet p_{T} (raw) distribution jet centb %d %s",icen,cjets[nj]),
					 bins,ptbins,rbins,rbinl,rbinh);
      hratiocorrrefpt[nj][icen]= new TH2F(Form("hratiocorrrefpt%d_%d",nj,icen),Form("Reco jet / Gen jet p_{T} (corr.) distribution jet centb %d %s",icen,cjets[nj]),
					  bins,ptbins,rbins,rbinl,rbinh);
      hratiocorrrawpt[nj][icen]= new TH2F(Form("hratiocorrrawpt%d_%d",nj,icen),Form("Correc. jet / Raw jet jet p_{T} distribution jet centb %d %s",icen,cjets[nj]),
                                          bins,ptbins,rbins,rbinl,rbinh);


      hratiocorrrefpt_lead[nj][icen]= new TH2F(Form("hratiocorrrefpt_lead%d_%d",nj,icen),Form("Leading jet Reco jet / Gen jet p_{T} (corr.) distribution jet centb %d %s",icen,cjets[nj]),
                                               bins,ptbins,rbins,rbinl,rbinh);
      hratiocorrrefpt_slead[nj][icen]= new TH2F(Form("hratiocorrrefpt_slead%d_%d",nj,icen),Form("sub-Leading jet Reco jet / Gen jet p_{T} (corr.) distribution jet centb %d %s",icen,cjets[nj]),
                                                bins,ptbins,rbins,rbinl,rbinh);
      hratiocorrrefpt_remain[nj][icen]= new TH2F(Form("hratiocorrrefpt_remain%d_%d",nj,icen),Form("Remaing jet Reco jet / Gen jet p_{T} (corr.) distribution jet centb %d %s",icen,cjets[nj]),
                                                 bins,ptbins,rbins,rbinl,rbinh);
      hratiocorrrefpt_pthat[nj][icen]= new TH2F(Form("hratiocorrrefpt_pthat%d_%d",nj,icen),Form("pT hat jet Reco jet / Gen jet p_{T} (corr.) distribution jet centb %d %s",icen,cjets[nj]),
						bins,ptbins,rbins,rbinl,rbinh);

      hratiocorrrefpt_genm[nj][icen]= new TH2F(Form("hratiocorrrefpt_genm%d_%d",nj,icen),Form("Gen matched jet Reco jet / Gen jet p_{T} (corr.) distribution jet centb %d %s",icen,cjets[nj]),
                                               bins,ptbins,rbins,rbinl,rbinh);
      hratiocorrorgpt[nj][icen]= new TH2F(Form("hratiocorrorgpt%d_%d",nj,icen),Form("Smeared Reco jet / Reco jet p_{T} (corr.) distribution jet centb %d %s",icen,cjets[nj]),
					  bins,ptbins,rbins,rbinl,rbinh);
      


      for(int ie=0;ie<2;ie++){
        hratiorawrefpt_eta[nj][icen][ie]= new TH2F(Form("hratiorawrefpt_eta%d_%d_%d",nj,icen,ie),
						   Form("Raw jet / Gen jet p_{T} (raw) distribution jet centb %d %s etabin%d",icen,cjets[nj],ie),
						   bins,ptbins,rbins,rbinl,rbinh);
        hratiocorrrefpt_eta[nj][icen][ie]= new TH2F(Form("hratiocorrrefpt_eta%d_%d_%d",nj,icen,ie),
						    Form("Reco jet / Gen jet p_{T} (raw) distribution jet centb %d %s etabin%d",icen,cjets[nj],ie),
						    bins,ptbins,rbins,rbinl,rbinh);
      }

      hjetptpu[nj][icen] = new TH2F(Form("hjetptpu%d_%d",nj,icen),Form("jet(pt:pu) distribution jet centb %d %s",icen,cjets[nj]),bins,ptbins,100,0,300);
      for(int ie=0;ie<ketar;ie++){
        hjetptpu_etab[nj][icen][ie] = new TH2F(Form("hjetptpu_etab%d_%d_%d",nj,icen,ie),Form("jet(pt:pu) distribution jet %s etabin %d cen %d",cjets[nj],icen,ie),
					       bins,ptbins,100,0,300);
      }
      
      for(int m=0;m<maxe;m++){
        hpteta[nj][icen][m] = new TH2F(Form("hpteta%d_%d_%d",nj,icen,m),Form("resolution  pt(eta) distribution cent %d jet %s etabin%d",icen,cjets[nj],m),
				       bins,ptbins,rbins,rbinl,rbinh);
      }
      
      for(int m=0;m<maxph;m++){
	hptphi[nj][icen][m] = new TH2F(Form("hptphi%d_%d_%d",nj,icen,m),Form("resolution pt(phi) distribution cent %d jet %s phibin%d",icen,cjets[nj],m),
				       bins,ptbins,rbins,rbinl,rbinh);
      }

      hgenjrecoj[nj][icen] =new TH2F(Form("hgenjrecoj%d_%d",nj,icen),Form("gen jet2 : reco jet2 %s cent %d",cjets[nj],icen),500,0.,1000.,500,0.,1000.);
      hgenjrecoj_pflavor[nj][icen] =new TH2F(Form("hgenjrecoj_pflavor%d_%d",nj,icen),Form("gen jet : reco jet with parton flavor %s cent %d",cjets[nj],icen),500,0.,1000.,500,0.,1000.);

      //! efficency histograms
      hPtAll [nj][icen] = new TH1F(Form("hPtAll%d_%d",nj,icen),Form("Denominator pT for algorithm %s cent %d",cjets[nj],icen),bins,ptbins);
      hEtaAll[nj][icen] = new TH1F(Form("hEtaAll%d_%d",nj,icen),Form("Denominator eta  for algorithm %s cent %d",cjets[nj],icen),20,-2.0,2.0);
      hPhiAll[nj][icen] = new TH1F(Form("hPhiAll%d_%d",nj,icen),Form("Denominator  phi  for algorithm %s cent %d",cjets[nj],icen),20,-pi,pi);
      
      hPtSel [nj][icen] = new TH1F(Form("hPtSel%d_%d",nj,icen),Form("Numerator pT for algorithm %s cent %d",cjets[nj],icen),bins,ptbins);
      hEtaSel[nj][icen] = new TH1F(Form("hEtaSel%d_%d",nj,icen),Form("Numerator eta  for algorithm %s cent %d",cjets[nj],icen),20,-2.0,2.0);
      hPhiSel[nj][icen] = new TH1F(Form("hPhiSel%d_%d",nj,icen),Form("Numerator  phi  for algorithm %s cent %d",cjets[nj],icen),20,-pi,pi);

      hDeltaR[nj][icen]    = new TH1F(Form("hDeltaR%d_%d",nj,icen),Form("#DeltaR for algorithm %s cent %d",cjets[nj],icen),100,0,1);
      hDeltaRAll[nj][icen] = new TH1F(Form("hDeltaRAll%d_%d",nj,icen),Form("#DeltaR (all) for algorithm %s cent %d",cjets[nj],icen),100,0,1);
      hDeltaRSel[nj][icen] = new TH1F(Form("hDeltaRSel%d_%d",nj,icen),Form("#DeltaR (sel) for algorithm %s cent %d",cjets[nj],icen),100,0,1);

      for(int ir=0;ir<25;ir++){
	//! Response vs DeltaR
        hRspVsDeltaR[nj][icen][ir] = new TH1F(Form("hRspVsDeltaR%d_%d_%d",nj,icen,ir),Form(" <recopt/refpt> vs. #DeltaR (%d) algorithm %s cent %d",ir,cjets[nj],icen),rbins,rbinl,rbinh);
      }
    }//! icen
  }//! nj
  std::cout<<"Initialized the histograms " <<std::endl;
  ///////////////////////////////////////////////////////////////////////////////////////// 


  //! Centrality reweighting function
  //TF1* fcen = new TF1("fcen","exp(-1.0*pow(x+1.11957e+01,2)/pow(1.34120e+01,2)/2)",0,40);
  TF1 *fcen = new TF1("fcen","[0]*exp([1]+[2]*x+[3]*x*x+[4]*x*x*x)",0,40);
  fcen->SetParameters(2.10653e-02,5.61607,-1.41493e-01,1.00586e-03,-1.32625e-04);

  //! vertex z reqeighting
  TF1 *fVz = new TF1("fVz","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
  if(strcmp(ksp,"pp")==0) fVz->SetParameters(8.41684e-01,-2.58609e-02,4.86550e-03,-3.10581e-04,2.07918e-05);
  else  fVz->SetParameters(7.62788e-01,-1.13228e-02,5.85199e-03,-3.04550e-04,4.43440e-05);

       
  Long64_t nb = 0;
  Long64_t nentries = c->GetEntries();
  std::cout<<Form("# of entries in TTree for %s : ",ksp)<<nentries<<std::endl;
  std::cout<<std::endl;
  hEvt->Fill(2,nentries);

  //! weight  for the merging of the samples for different pT hat bins
  Float_t wxs = xsection/(nentries/100000.);
  
  Int_t iEvent=0; 
  for (Long64_t ievt=0; ievt<nentries;ievt++) {//! event loop
    //for (Long64_t ievt=0; ievt<200;ievt++) {//! event loop
    //! load the hiForest event
    nb += c->GetEntry(ievt);

    int hiBin       = c->evt.hiBin;
    float vz        = c->evt.vz;
    float hiHF      = c->evt.hiHF;

    //! testing
    //if(hiBin>4 && strcmp(ksp,"pbpb")==0)continue;

    //! apply vertex cut
    if(fabs(vz)>kVzcut)continue;


    //! Centrality bin
    if(hiBin<0 || hiBin>39)continue;
    int centb=-1;
    double wcen=1;
    double wvz=fVz->Eval(vz);
    if(strcmp(ksp,"pbpb")==0){
      centb=GetCentBin(hiBin);
      wcen = fcen->Eval(hiBin);
    }else{
      centb=ncen-1; //! pp
      wcen=1;
    }

    if(ievt%10000==0)std::cout<<" ********** Event # " <<ievt<<std::endl;
    //std::cout<<" ********** Event # " <<ievt<<"\t vz : "<<vz<<"\t hiBin : "<<hiBin<<"\t wxs : "<<wxs<<"\t centb : "<<centb<<std::endl;
    

    //! Centrality from 0-90% 
    if(centb==-1 || centb==ncen)continue;


    //! Methods to get the jet objects of a given algorithm
    /* Get the Jet object of a given type of algorithm according to the number from Full list
       iJet = c->GetJet(15);
       std::cout<<"\t GetJet()       : "<<c->GetCjets[Nj](15)<<"\t # of Jets  : "<<iJet->nref<<std::endl;
    */
    /* Get the Jet object of a given type of algorithm according to the NAME from Full list
    iJet = c->GetJetByAlgo("akPu3PF");
    std::cout<<"\t GetJetByAlgo() : "<<c->GetCjets[Nj](15)<<"\t # of Jets  : "<<iJet->nref<<std::endl;
    */

    int istat=0;
    for(int nj=0;nj<knj;nj++){ //! loop over different jet algorithms
      //for(int nj=15;nj<16;nj++){ //! loop over different jet algorithms
      //! Get the jet object
      //iJet = c->GetJet(nj);
      iJet = c->GetJetByAlgo(cjets[nj]);

      //! xsec-weight
      double pthat = iJet->pthat;
      if(pthat > maxpthat)continue;
      istat=1;
      
      //std::cout<<"\t Jet Algorithm : "<<c->GetCjets[Nj](nj)<<"\t # of Jets  : "<<iJet->nref<<"\t pthat : "<<pthat<<std::endl;
      if(nj==0)hTotEve->Fill(1); //! akPu3PF      

      int *ljet = new int[2];
      FindLeadSubLeadJets(iJet,ljet);
      if(ljet[0]>=0){
	hratiocorrrefpt_lead[nj][centb]->Fill(iJet->refpt[ljet[0]],iJet->jtpt[ljet[0]]/iJet->refpt[ljet[0]],wxs*wcen*wvz);
      }
      if(ljet[1]>=0){
	hratiocorrrefpt_slead[nj][centb]->Fill(iJet->refpt[ljet[1]],iJet->jtpt[ljet[1]]/iJet->refpt[ljet[1]],wxs*wcen*wvz);
      }

      if(nj>3){
	//! for gen matched jets
	for(int igen=0; igen<iJet->ngen; igen++){
	  if( iJet->gensubid[igen] != 0) continue;
	  int gj = iJet->genmatchindex[igen];
	  
	  float refpt   = iJet->refpt[gj];
	  float recopt  = iJet->jtpt[gj];
	  float recoeta = iJet->jteta[gj];
	  float delr    = iJet->refdrjt[gj];
	  
	  if(recopt<kptrecocut || refpt<kptgencut || refpt==0 || fabs(recoeta)>ketacut || fabs(delr)>kdRcut)continue;
	  hratiocorrrefpt_genm[nj][centb]->Fill(refpt,recopt/refpt,wxs*wcen*wvz);
	}
      }

      float njets=0.;
      float meanpt=0.;

      
      //int istat=0;
      for(int ij=0; ij<iJet->nref; ij++){
	//if(!selecJet(iJet,ij))continue;

	//! ic:6 is the corrected pT without any smearing
	//! ic:0-5 pp smeared pT with heavy-ion resolution

	float refpt   = iJet->refpt[ij];
	float recopt  =iJet->jtpt[ij];
	float rawpt   = iJet->rawpt[ij];
	float refeta  = iJet->refeta[ij];
	float recoeta = iJet->jteta[ij];
	float refphi  = iJet->refphi[ij];
	float recophi = iJet->jtphi[ij];
	int refparton_flavor = iJet->refparton_flavor[ij];


	if(refpt<40)continue;

	
	//! Matching criteria 
	float delr = iJet->refdrjt[ij];
	
	if(recopt<kptrecocut || refpt<kptgencut || refpt==0)continue;
	
	if(fabs(refeta)<ketacut && refpt>80){
	  //! Denominator for matching efficiency
	  hPtAll [nj][centb]->Fill(refpt,wxs*wcen*wvz);
	  hEtaAll[nj][centb]->Fill(refeta,wxs*wcen*wvz);
	  hPhiAll[nj][centb]->Fill(refphi,wxs*wcen*wvz);
	  
	  //! Response
	  for (int idr=0;idr<25;idr++) {
	    double drcut = 0.0+idr*(0.25-0.00)/(25-1);
	    if (fabs(delr)>drcut) continue;
	    hRspVsDeltaR[nj][centb][idr]->Fill(recopt/refpt,wxs*wcen*wvz);
	  }
	  
	  //! DeltaR efficiency
	  hDeltaR[nj][centb]->Fill(delr,wxs*wcen*wvz);
	  for (int idrmax=0;idrmax<100;idrmax++) {
	    float drmax = idrmax*0.01+0.005;
	    hDeltaRAll[nj][centb]->Fill(drmax,wxs*wcen*wvz);
	    if (delr<drmax) hDeltaRSel[nj][centb]->Fill(drmax,wxs*wcen*wvz);
	  }
	}
	
	//! Matching cut for gen-jet and reco-jet
	if(fabs(delr) > kdRcut)continue;
	
	
	if(fabs(refeta)<ketacut && refpt>80){
	  //! Numerator for matching efficiency
	  hPtSel [nj][centb]->Fill(refpt,wxs*wcen*wvz);
	  hEtaSel[nj][centb]->Fill(refeta,wxs*wcen*wvz);
	  hPhiSel[nj][centb]->Fill(refphi,wxs*wcen*wvz);
	}
	
	//! pile up eta dependence
	//int ebin = GetDetEtaBin(fabs(recoeta));
	
	//! jet selction cut
	if(fabs(recoeta)>ketacut)continue;

	//istat=1;
	//! 1D distributions
	njets++;
	meanpt += refpt;
	hNjets [nj][centb]->Fill(1.);
	hgenpt [nj][centb]->Fill(refpt,wxs*wcen*wvz);
	hgeneta[nj][centb]->Fill(refeta,wxs*wcen*wvz);
	hgenphi[nj][centb]->Fill(refphi,wxs*wcen*wvz);
	hgenptC[nj][centb]->Fill(refpt,wxs*wcen*wvz);	
	hrecoeta[nj][centb]->Fill(recoeta,wxs*wcen*wvz);
	hrecophi[nj][centb]->Fill(recophi,wxs*wcen*wvz);
	hrawpt  [nj][centb]->Fill(rawpt,wxs*wcen*wvz);
	hrawptC [nj][centb]->Fill(rawpt,wxs*wcen*wvz);	

	hFracjets[nj][centb]->Fill(pthat,refpt,wxs*wcen*wvz);

	int ieta=-1;
	if(fabs(recoeta)<1.3)ieta=0; //! barrel region
	else ieta=1; //! HCAL region

	  //! Response in different eta and phi bins
	int etabin = GetEtaBin(fabs(refeta));
	int phibin = GetPhiBin(refphi);
	if(etabin < 0 || etabin>=maxe || phibin < 0 || phibin>=maxph)continue;
	
	//! pileup study
	hjetptpu[nj][centb]->Fill(recopt,iJet->jtpu[ij],wxs*wcen*wvz);
	
	for(int ic=0;ic<ncen;ic++){
	  
	  float smpt = GetSmearedPtMC(nj,ic,recopt,refpt);
	  //if(nj==2)std::cout<<"\t \t "<<ij<<"\t refpt : "<<refpt<<"\t recopt : "<<recopt<<"\t raw pT : "<<rawpt<<"\t smpt "<<smpt<<std::endl;
	  
	  /*
	  //! For different centrality binning
	  int nic =-1;
	  if(ic==0 || ic==1)nic=0;
	  else if(ic==2)nic=1;
	  else if(ic==3)nic=2;
	  else if(ic==4 || ic==5)nic=3;

	  //!  and replace ic with nic in your code
	  //! make sure that your bins are coarser or equal to these centrality bins
	  //!   0     1      2       3       4       5 
	  //!  0-5%   5-10%  10-30%  30-50%  50-70%  70-90%
	  */

	  hrecopt [nj][ic]->Fill(recopt,wxs*wcen*wvz);
	  hsmpt   [nj][ic]->Fill(smpt,wcen*wvz*wxs);	  
	  
	  hrecoptC[nj][ic]->Fill(recopt,wxs*wcen*wvz);
	  hsmptC  [nj][ic]->Fill(smpt,wcen*wvz*wxs);
	  
	  hrecogen[nj][ic]->Fill(refpt,smpt/refpt,wxs*wcen*wvz);
	  hrawgen [nj][ic]->Fill(refpt,rawpt/refpt,wxs*wcen*wvz);
	  hrecoraw[nj][ic]->Fill(rawpt,smpt/rawpt,wxs*wcen*wvz);
	  hrecoraw_ref [nj][ic]->Fill(refpt,smpt/rawpt,wxs*wcen*wvz);
	  
	  //! Response  (ratio of recopt/refpt)
	  hratiocorrrefpt[nj][ic]->Fill(refpt,smpt/refpt,wxs*wcen*wvz);
	  hratiorawrefpt [nj][ic]->Fill(refpt,rawpt/refpt,wxs*wcen*wvz);
	  hratiocorrrawpt[nj][ic]->Fill(recopt,smpt/rawpt,wxs*wcen*wvz);
	  
	  //! remaing jets
	  bool iRemain = ij!=ljet[0] || ij!=ljet[1];
	  if(iRemain)hratiocorrrefpt_remain[nj][ic]->Fill(refpt,smpt/refpt,wxs*wcen*wvz);
	  
	  //! pT hat
	  hratiocorrrefpt_pthat[nj][ic]->Fill(pthat,smpt/refpt,wxs*wcen*wvz);	
	  
	  //! 2D correlation between refpt and recopt
	  hcorrptrefpt[nj][ic]->Fill(refpt,smpt,wxs*wcen*wvz);
	  hrawptrefpt [nj][ic]->Fill(refpt,rawpt,wxs*wcen*wvz);
	  hcorrptrawpt[nj][ic]->Fill(rawpt,smpt,wxs*wcen*wvz);
	  
	  //! Very fine bin in ref pt
	  hrescrpt[nj][ic]->Fill(refpt,smpt/refpt,wxs*wcen*wvz);
	  hresrrpt[nj][ic]->Fill(refpt,rawpt/refpt,wxs*wcen*wvz);
	  hresrcrpt[nj][ic]->Fill(recopt,smpt/rawpt,wxs*wcen*wvz);
	  hresorpt[nj][ic]->Fill(refpt,smpt/recopt,wcen*wvz*wxs); //! finner bins

	  //! Ratio of original pT and the smeared pT
	  hratiocorrorgpt[nj][ic]->Fill(recopt,smpt/recopt,wcen*wvz*wxs);

	  if(ieta>=0){
	    hratiocorrrefpt_eta[nj][ic][ieta]->Fill(refpt,smpt/refpt,wxs*wcen*wvz);
	    hratiorawrefpt_eta [nj][ic][ieta]->Fill(refpt,rawpt/refpt,wxs*wcen*wvz);
	  }

	  //! Response in eta and phi bins
	  hpteta[nj][ic][etabin]->Fill(refpt,smpt/refpt,wxs*wcen*wvz);
	  hptphi[nj][ic][phibin]->Fill(refpt,smpt/refpt,wxs*wcen*wvz);
	  
	  //! Gen:Reco correlation plots
	  hgenjrecoj [nj][ic]->Fill(refpt,smpt,wxs*wcen*wvz);
	  if(fabs(refparton_flavor)<=21)hgenjrecoj_pflavor [nj][ic]->Fill(refpt,smpt,wxs*wcen*wvz);

	}//! ic loop
      }//! ij loop

      hNevt[nj][centb]->Fill(vz);
      if(njets>0){
	hAvjets[nj][centb]->Fill(pthat,njets,wxs*wcen*wvz);
	hMeanPtjets[nj][centb]->Fill(pthat,meanpt/njets,wxs*wcen*wvz);
      }
      delete [] ljet;
    }//! nj jet loop
    
    if(istat){
      hBin->Fill(hiBin);
      hVz->Fill(vz);
      hHF->Fill(hiHF,wxs*wcen*wvz);
      iEvent++;
    }
    //std::cout<<"Completed event #  "<<ievt<<std::endl; 
  }//! event loop ends
  
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;

  std::cout<<"Events which passed the pT hat cut : "<<hTotEve->Integral()<<" out of  : "<<hEvt->Integral()
	   <<" efficiency of all the cuts : " <<hTotEve->Integral()/hEvt->Integral()<<std::endl;
  std::cout<<std::endl;


  if(strcmp(ksp,"pp")==0){
    for(int nj=0;nj<knj;nj++){
      //std::cout<<"# of Events for : "<<c->GetCjets[Nj](nj)<<"\t"<<hNevt[nj][0]->Integral()<<"\t # of Jets : "<<hNjets[nj][0]->Integral()<<std::endl;
      std::cout<<"# of Events for : "<<cjets[nj]<<"\t"<<hNevt[nj][ncen-1]->Integral()<<"\t # of Jets : "<<hNjets[nj][ncen-1]->Integral()<<std::endl;
    }
  }else{
    for(int nj=0;nj<knj;nj++){
      for(int icen=0;icen<ncen;icen++){
	std::cout<<"# of Events for : "<<cjets[nj]<<"\t icen : "<<icen<<"\t"<<hNevt[nj][icen]->Integral()<<"\t # of Jets : "<<hNjets[nj][icen]->Integral()<<std::endl;
      }
      std::cout<<std::endl;
    }
  }

  //! Write to output file
  fout->cd();
  fout->Write();
  fout->Close();


  //! Check
  timer.Stop();
  float  mbytes = 0.000001*nb;
  double rtime  = timer.RealTime();
  double ctime  = timer.CpuTime();

  std::cout<<std::endl;
  std::cout<<Form("RealTime=%f seconds, CpuTime=%f seconds",rtime,ctime)<<std::endl;
  std::cout<<Form("You read %f Mbytes/RealTime seconds",mbytes/rtime)<<std::endl;
  std::cout<<Form("You read %f Mbytes/CpuTime  seconds",mbytes/ctime)<<std::endl;
  std::cout<<std::endl;
  std::cout<<"Good bye : " <<"\t"<<std::endl;
  return 1;
}
int GetCentBin(int bin)
{
  int ibin=-1;
  //! centrality is defined as 2.5% bins of cross section
  //! in 0-39 bins

  if(bin<2)ibin=0; //! 0-5%
  else if(bin>=2 && bin<4)ibin=1;   //! 5-10%
  else if(bin>=4 && bin<12)ibin=2;  //! 10-30%
  else if(bin>=12&& bin<20)ibin=3;  //! 30-50%
  else if(bin>=20&& bin<28)ibin=4;  //! 50-70%
  else if(bin>=28&& bin<36)ibin=5;  //! 70-90%

  return ibin;
}
int GetEtaBin(float eta)
{
  int ibin=-1;
  float min =0.0;
  ibin = (eta - min)*maxe/2.;
  return ibin;
}
int GetPhiBin(float phi)
{
  int ibin=-1;
  float min = -pi;
  ibin = (phi - min)*maxph/2*pi;
  return ibin;
}
bool selecJet(Jets *jetc,int it)
{
  bool goodJet = jetc->refpt[it]!=-999 && jetc->refphi[it]!=-999 && jetc->refeta[it]!=-999 && jetc->jtpt[it]!=-999 && jetc->jtphi[it]!=-999 && jetc->jteta[it]!=-999;
  return goodJet;
}
int GetDetEtaBin(float eta)
{
  int ibin=-1;
  if(eta>=0 && eta<1.3)ibin=0;       //! barrel
  else if(eta>=1.3 && eta<2.0)ibin=1;//! endcap+tracks
  else if(eta>=2.0 && eta<3.0)ibin=2;//! endcap-notracks
  else if(eta>=3.0 && eta<5.0)ibin=3;//! forward

  return ibin;
}
int GetPtBin(float pt)
{
  for(int ix=0;ix<bins;ix++){
    if(pt>=ptbins[ix] && pt<ptbins[ix+1]){
      return ix;
    }
  }
  return -1;
}
void FindLeadSubLeadJets(Jets *jetc, int *ljet)
{
  ljet[0]=-1; ljet[1]=-2;

  float tempt1=-9;
  float tempt2=-9;

  //! Get the leading jet                                                                                                                                                                                      
  for(int ij=0; ij<jetc->nref; ij++){
    //if(!selecJet(jetc,ij))continue;
    if(fabs(jetc->refdrjt[ij])>kdRcut || fabs(jetc->jteta[ij])>ketacut || jetc->jtpt[ij]<kptrecocut)continue;

    float jetpt = jetc->jtpt[ij];
    if(jetpt > tempt1){
      tempt1 = jetpt;
      ljet[0] = ij;
    }
  }

  for(int ij=0; ij<jetc->nref; ij++){
    //if(!selecJet(jetc,ij) || ij==ljet[0])continue;
    if(ij==ljet[0])continue;
    if(fabs(jetc->refdrjt[ij])>kdRcut || fabs(jetc->jteta[ij])>ketacut || jetc->jtpt[ij]<kptrecocut)continue;

    float jetpt = jetc->jtpt[ij];
    //float dphi  = jetc->jtphi[ij] - jetc->jtphi[ljet[0]];
    //if (dphi > pi ) dphi = dphi - 2 * pi;
    //if (dphi < -pi) dphi = dphi + 2 * pi;
    //if(dphi < 2*pi/3.)continue;
    if (jetpt > tempt2){
      tempt2 = jetpt;
      ljet[1] = ij;
    }
  }
}

