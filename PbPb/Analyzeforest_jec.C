#include <TF1.h>
#include <TProfile.h>
#include <TChain.h>
#include <TTree.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TFile.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom3.h>

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <cstring>
#include <vector>
#include <algorithm>


#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

using namespace std;


#ifdef _MAKECINT_
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#endif

int GetEtaBin(float /*eta*/);
int GetDetEtaBin(float /*eta*/);
int GetPhiBin  (float /*phi*/);
int GetCentBin(int /*hiBin*/);
int GetPtBin(float /*pt*/);
double delphi(double /*phi1*/, double /*phi2*/);


//! pt binning
double ptbins[] ={10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,280,300,320,340,360,400,450,500,550,600,835};
const int bins = sizeof(ptbins)/sizeof(double) - 1;

//! data pt binning
//double ptbins_data[] ={100, 110, 120, 130, 140, 150, 160, 170, 180, 200, 240, 300};
//const int dbins = sizeof(ptbins_data)/sizeof(double) - 1;

//! constants
int iYear=2014;
const double pi=acos(-1.);
const double pi2=2*pi -1;

const int ncent=6; //! for pp fill in 0, for ppb 0-4 and pbpb 0-5
const int ketar=4;
const int maxe =5;
const int maxph=5;
const float ketacut=2.0;
const double kptrecocut=15.;
const double kptgencut =0.;
const double kdRcut=0.3;


const int rbins=50;
double rbinl=0.0;
double rbinh=5.0;


 const int kAlgos = 18;
 const char *calgo[kAlgos] = {"akVs3PF","akVs4PF","akVs5PF",
  			     "akVs3Calo","akVs4Calo","akVs5Calo",
  			     "akPu3PF","akPu4PF","akPu5PF",
  			     "akPu3Calo","akPu4Calo","akPu5Calo",
  			     "ak3PF","ak4PF","ak5PF",
  			     "ak3Calo","ak4Calo","ak5Calo"
 };
 string corrFileName[kAlgos]= {"AKVs3PF","AKVs4PF","AKVs5PF",
  			      "AKVs3Calo","AKVs4Calo","AKVs5Calo",
   			      "AKPu3PF","AKPu4PF","AKPu5PF",
  			      "AKPu3Calo","AKPu4Calo","AKPu5Calo",
  			      "AK3PF","AK4PF","AK5PF",
  			      "AK3Calo","AK4Calo","AK5Calo",
  // 			      "AK3PF","AK4PF","AK5PF",
  // 			      "AK3Calo","AK4Calo","AK5Calo"
 };

 const float kDelR[kAlgos]  = {0.3,0.4,0.4,
    			      0.3,0.4,0.5,
    			      0.3,0.4,0.5,
    			      0.3,0.4,0.5,
    			      0.3,0.4,0.5,
    			      0.3,0.4,0.5
 };



// const int kAlgos = 4;
// const char *calgo[kAlgos] = {"akVs3PF","akVs3Calo",
//  			     "akPu3PF","akPu3Calo",
// };
// string corrFileName[kAlgos]= {"AKVs3PF","AKVs3Calo",
//   			      "AKPu3PF","AKPu3Calo",
// };

// const float kDelR[kAlgos]  = {0.3,
//    			      0.3,
//    			      0.3,
//    			      0.3,
// };

const char *ccent[ncent] = {"0-5%","5-10%","10-30%","30-50%","50-70%","70-90%"};
TStopwatch timer;
int Analyzeforest_jec(const char *ksp="pp")
{

  timer.Start();

  std::string inname="";
  int knj = kAlgos;
  if(strcmp(ksp,"pp")==0)inname="dijet_pp_mergedpthatbins_Track8_Jet22MC.root";
  else {
    inname="dijet_pp_mergedpthatbins_embedpp.root";
    knj = 12;
    //knj = 4;
  }
  TFile *fin = new TFile(inname.c_str(),"r");

  std::string outname="";  
  if(strcmp(ksp,"pp")==0)outname="JetResponse_histos_JECv14_pp.root";  
  else outname="JetResponse_histos_JECv14_embeded_pp_cent.root";
  TFile *fout = new TFile(outname.c_str(),"RECREATE");
  
  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"**************************************************** "<<std::endl;
  std::cout<<Form("Analyzeforest_jec : %s",ksp)<<std::endl;
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

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  TH1F *hBin    = new TH1F("hBin","Centrality bin",200,-0.5,200-0.5);
  hBin->Sumw2();
  TH1F *hpthat  = new TH1F("hpthat","pt-hat distribution",500,0,1000);
  hpthat->Sumw2();

  TH2F *hrajptrpt [knj][ncent];

  //! Gen matched jets
  TH1F *hgenpt_genm [knj][ncent], *hrecopt_genm[knj][ncent], *hrawpt_genm[knj][ncent];
  TH1F *hjeteta     [knj][ncent], *hjetphi     [knj][ncent];
  TH2F *hjetpteta   [knj][ncent], *hjetptphi   [knj][ncent], *hjetetaphi [knj][ncent];

  //! Ratios of the pt distributions
  TProfile *hrecogen[knj][ncent], *hrecoraw[knj][ncent], *hrawgen[knj][ncent];

  //! Resposnse
  TH2F *hrescrpt_genm [knj][ncent], *hresrrpt_genm [knj][ncent], *hresrcrpt_genm [knj][ncent];
  TH3F *hrescreta_genm[knj][ncent], *hresrreta_genm[knj][ncent], *hresrcreta_genm[knj][ncent];
  TH3F *hrescrphi_genm[knj][ncent], *hresrrphi_genm[knj][ncent], *hresrcrphi_genm[knj][ncent];
  TH2F *hratiorawrefpt_eta[knj][ncent][2], *hratiocorrrefpt_eta[knj][ncent][2];

  TH2F *hratiocorrrefpt_genm[knj][ncent];

  TH2F *hpteta[knj][ncent][maxe] ;
  TH2F *hptphi[knj][ncent][maxph] ;

  TH2F *hgenjrecoj[knj][ncent];
  TH2F *hgenjrawj [knj][ncent];

  //! Background jet pt
  TH2F *hjbkgComb[knj][ncent], *hjbkgComb_match[knj][ncent];


  //! For comparison with data
  TH2F *hJetEnergyScale[knj][ncent];


  //! Efficency histos 
  TH1F *hdRAll [knj][ncent], *hdRSel[knj][ncent];
  TH1F *hPtAll [knj][ncent], *hPtSel[knj][ncent];
  TH1F *hEtaAll[knj][ncent], *hEtaSel[knj][ncent];
  TH1F *hPhiAll[knj][ncent], *hPhiSel[knj][ncent];


  TH1F *hFakePtAll [knj][ncent], *hFakePtSel[knj][ncent];
  TH1F *hFakeEtaAll[knj][ncent], *hFakeEtaSel[knj][ncent];
  TH1F *hFakePhiAll[knj][ncent], *hFakePhiSel[knj][ncent];


  for(int nj=0;nj<knj;nj++){
    for(int ic=0;ic<ncent;ic++){

      hgenpt_genm [nj][ic] = new TH1F(Form("hgenpt_genm%d_%d",nj,ic),Form("Gen matched gen p_{T} distribution jet centb %d %s",ic,calgo[nj]),500,0,1000);
      hgenpt_genm [nj][ic]->Sumw2();
      hrecopt_genm[nj][ic] = new TH1F(Form("hrecopt_genm%d_%d",nj,ic),Form("Gen matched reco p_{T} distribution jet centb %d %s",ic,calgo[nj]),500,0,1000);
      hrecopt_genm[nj][ic]->Sumw2();
      hrawpt_genm [nj][ic] = new TH1F(Form("hrawpt_genm%d_%d",nj,ic),Form("Gen matched raw p_{T} distribution jet centb %d %s",ic,calgo[nj]),500,0,1000);
      hrawpt_genm [nj][ic]->Sumw2();

      hjbkgComb[nj][ic] = new TH2F(Form("hjbkgComb%d_%d",nj,ic),Form("jet background distribution jet centb %d %s",ic,calgo[nj]),500,0,1000,500,0.,250.);
      hjbkgComb[nj][ic]->Sumw2();

      hjbkgComb_match[nj][ic] = new TH2F(Form("hjbkgComb_match%d_%d",nj,ic),Form("match jet background distribution jet centb %d %s",ic,calgo[nj]),500,0,1000,500,0.,250.);
      hjbkgComb_match[nj][ic]->Sumw2();

      //! Ratios
      hrecogen[nj][ic] = new TProfile(Form("hrecogen%d_%d",nj,ic),Form("reco/gen p_{T} distribution jet centb %d %s",ic,calgo[nj]),500,0,1000);
      hrecogen[nj][ic]->Sumw2();
      hrecoraw[nj][ic] = new TProfile(Form("hrecoraw%d_%d",nj,ic),Form("reco/raw p_{T} distribution jet centb %d %s",ic,calgo[nj]),500,0,1000);
      hrecoraw[nj][ic]->Sumw2();
      hrawgen[nj][ic]  = new TProfile(Form("hrawgen%d_%d",nj,ic),Form("raw/gen p_{T} distribution jet centb %d %s",ic,calgo[nj]),500,0,1000);
      hrawgen[nj][ic]->Sumw2();

      //! Gen matched Response and resolution
      hrescrpt_genm[nj][ic]= new TH2F(Form("hrescrpt_genm%d_%d",nj,ic),Form("Gen jet:(Reco/Gen) jet p_{T} distribution jet centb %d %s",ic,calgo[nj]),500,0,1000,150,rbinl,rbinh);
      hrescrpt_genm[nj][ic]->Sumw2();
      hresrrpt_genm[nj][ic]= new TH2F(Form("hresrrpt_genm%d_%d",nj,ic),Form("Gen jet:(Raw/Gen) jet p_{T}  distribution jet centb %d %s",ic,calgo[nj]),500,0,1000,150,rbinl,rbinh);
      hresrrpt_genm[nj][ic]->Sumw2();
      hresrcrpt_genm[nj][ic]= new TH2F(Form("hresrcrpt_genm%d_%d",nj,ic),Form("Reco jet:(Reco/Raw) jet p_{T} distribution jet centb %d %s",ic,calgo[nj]),500,0,1000,150,rbinl,rbinh);
      hresrcrpt_genm[nj][ic]->Sumw2();

      //! Gen matched response and resolutions in different eta bins
      hrescreta_genm[nj][ic]  = new TH3F(Form("hrescreta_genm%d_%d",nj,ic),Form("Gen jet:(Reco/Gen) jet p_{T} distribution jet centb %d %s",ic,calgo[nj])  ,50,-2.5,2.5,500,0,1000,rbins,rbinl,rbinh);
      hrescreta_genm[nj][ic]->Sumw2();
      hresrreta_genm[nj][ic]  = new TH3F(Form("hresrreta_genm%d_%d",nj,ic),Form("Gen jet:(Raw/Gen) jet p_{T}  distribution jet centb %d %s",ic,calgo[nj])  ,50,-2.5,2.5,500,0,1000,rbins,rbinl,rbinh);
      hresrreta_genm[nj][ic]->Sumw2();
      hresrcreta_genm[nj][ic] = new TH3F(Form("hresrcreta_genm%d_%d",nj,ic),Form("Reco jet:(Reco/Raw) jet p_{T} distribution jet centb %d %s",ic,calgo[nj]),50,-2.5,2.5,500,0,1000,rbins,rbinl,rbinh);
      hresrcreta_genm[nj][ic]->Sumw2();

      //! Gen matched response and resolutions in different phi bins
      hrescrphi_genm[nj][ic]  = new TH3F(Form("hrescrphi_genm%d_%d",nj,ic),Form("Gen jet:(Reco/Gen) jet p_{T} distribution jet centb %d %s",ic,calgo[nj])  ,72,-pi,pi,500,0,1000,rbins,rbinl,rbinh);
      hrescrphi_genm[nj][ic]->Sumw2();
      hresrrphi_genm[nj][ic]  = new TH3F(Form("hresrrphi_genm%d_%d",nj,ic),Form("Gen jet:(Raw/Gen) jet p_{T}  distribution jet centb %d %s",ic,calgo[nj])  ,72,-pi,pi,500,0,1000,rbins,rbinl,rbinh);
      hresrrphi_genm[nj][ic]->Sumw2();
      hresrcrphi_genm[nj][ic] = new TH3F(Form("hresrcrphi_genm%d_%d",nj,ic),Form("Reco jet:(Reco/Raw) jet p_{T} distribution jet centb %d %s",ic,calgo[nj]),72,-pi,pi,500,0,1000,rbins,rbinl,rbinh);
      hresrcrphi_genm[nj][ic]->Sumw2();

      hjeteta[nj][ic] = new TH1F(Form("hjeteta%d_%d",nj,ic),Form("jet eta distribution jet centb %d %s",ic,calgo[nj]),72,-ketacut,ketacut);
      hjeteta[nj][ic]->Sumw2();
      hjetphi[nj][ic] = new TH1F(Form("hjetphi%d_%d",nj,ic),Form("jet phi distribution jet centb %d %s",ic,calgo[nj]),72,-pi,pi);
      hjetphi[nj][ic]->Sumw2();

      hjetetaphi[nj][ic] = new TH2F(Form("hjetetaphi%d_%d",nj,ic),Form("jet eta-phi distribution jet centb %d %s",ic,calgo[nj]),72,-ketacut,ketacut,72,-pi,pi);
      hjetetaphi[nj][ic] ->Sumw2();
      hjetpteta[nj][ic] = new TH2F(Form("hjetpteta%d_%d",nj,ic),Form("jet pt-eta distribution jet centb %d %s",ic,calgo[nj]),500,0,1000,72,-ketacut,ketacut);
      hjetpteta[nj][ic]->Sumw2();
      hjetptphi[nj][ic] = new TH2F(Form("hjetptphi%d_%d",nj,ic),Form("jet pt-phi distribution jet centb %d %s",ic,calgo[nj]),500,0,1000,72,-pi,pi);
      hjetptphi[nj][ic]->Sumw2();

      hratiocorrrefpt_genm[nj][ic]= new TH2F(Form("hratiocorrrefpt_genm%d_%d",nj,ic),Form("Gen matched jet Reco jet / Gen jet p_{T} (corr.) distribution jet centb %d %s",ic,calgo[nj]),
                                               500,0,1000,rbins,rbinl,rbinh);
      hratiocorrrefpt_genm[nj][ic]->Sumw2();

      for(int ie=0;ie<2;ie++){
        hratiorawrefpt_eta[nj][ic][ie]= new TH2F(Form("hratiorawrefpt_eta%d_%d_%d",nj,ic,ie),
						   Form("Raw jet / Gen jet p_{T} (raw) distribution jet centb %d %s etabin%d",ic,calgo[nj],ie),
						   500,0,1000,rbins,rbinl,rbinh);
	hratiorawrefpt_eta[nj][ic][ie]->Sumw2();
        hratiocorrrefpt_eta[nj][ic][ie]= new TH2F(Form("hratiocorrrefpt_eta%d_%d_%d",nj,ic,ie),
						    Form("Reco jet / Gen jet p_{T} (raw) distribution jet centb %d %s etabin%d",ic,calgo[nj],ie),
						    500,0,1000,rbins,rbinl,rbinh);
	hratiocorrrefpt_eta[nj][ic][ie]->Sumw2();
      }
      
      for(int m=0;m<maxe;m++){
        hpteta[nj][ic][m] = new TH2F(Form("hpteta%d_%d_%d",nj,ic,m),Form("resolution  pt(eta) distribution cent %d jet %s etabin%d",ic,calgo[nj],m),
				       500,0,1000,rbins,rbinl,rbinh);
	hpteta[nj][ic][m]->Sumw2();
      }
      
      for(int m=0;m<maxph;m++){
	hptphi[nj][ic][m] = new TH2F(Form("hptphi%d_%d_%d",nj,ic,m),Form("resolution pt(phi) distribution cent %d jet %s phibin%d",ic,calgo[nj],m),
				       500,0,1000,rbins,rbinl,rbinh);
	hptphi[nj][ic][m]->Sumw2();
      }

      hgenjrecoj[nj][ic] =new TH2F(Form("hgenjrecoj%d_%d",nj,ic),Form("gen jet2 : reco jet2 %s cent %d",calgo[nj],ic),500,0.,1000.,500,0.,1000.);
      hgenjrecoj[nj][ic]->Sumw2();


      hgenjrawj[nj][ic] =new TH2F(Form("hgenjrawj%d_%d",nj,ic),Form("gen jet2 : raw jet2 %s cent %d",calgo[nj],ic),500,0.,1000.,500,0.,1000.);
      hgenjrawj[nj][ic]->Sumw2();

      hrajptrpt[nj][ic] = new TH2F(Form("hrajptrpt%d_%d",nj,ic),Form("corr pT / jet(raw pt) p_{T} distribution %d jet %s",ic,calgo[nj]),500,0,1000,50,0,10);
      hrajptrpt[nj][ic]->Sumw2();

      hJetEnergyScale[nj][ic] = new TH2F(Form("hJetEnergyScale%d_%d",nj,ic),Form("hJetEnergyScale%d_%d",nj,ic),500,0,1000,50,-1.00,1.00);
      hJetEnergyScale[nj][ic]->Sumw2();

      //! efficcy histograms
      hdRAll [nj][ic] = new TH1F(Form("hdRAll%d_%d",nj,ic),Form("dR Denominator pT for algorithm %s cent %d",calgo[nj],ic),500,0,1000);
      hdRAll [nj][ic]->Sumw2();
      hPtAll [nj][ic] = new TH1F(Form("hPtAll%d_%d",nj,ic),Form("Denominator pT for algorithm %s cent %d",calgo[nj],ic),500,0,1000);
      hPtAll [nj][ic]->Sumw2();
      hEtaAll[nj][ic] = new TH1F(Form("hEtaAll%d_%d",nj,ic),Form("Denominator eta  for algorithm %s cent %d",calgo[nj],ic),20,-ketacut,ketacut);
      hEtaAll[nj][ic]->Sumw2();
      hPhiAll[nj][ic] = new TH1F(Form("hPhiAll%d_%d",nj,ic),Form("Denominator  phi  for algorithm %s cent %d",calgo[nj],ic),20,-pi,pi);
      hPhiAll[nj][ic]->Sumw2();

      hdRSel [nj][ic] = new TH1F(Form("hdRSel%d_%d",nj,ic),Form("dR Numerator pT for algorithm %s cent %d",calgo[nj],ic),500,0,1000);
      hdRSel [nj][ic]->Sumw2();
      hPtSel [nj][ic] = new TH1F(Form("hPtSel%d_%d",nj,ic),Form("Numerator pT for algorithm %s cent %d",calgo[nj],ic),500,0,1000);
      hPtSel [nj][ic]->Sumw2();
      hEtaSel[nj][ic] = new TH1F(Form("hEtaSel%d_%d",nj,ic),Form("Numerator eta  for algorithm %s cent %d",calgo[nj],ic),20,-ketacut,ketacut);
      hEtaSel[nj][ic]->Sumw2();
      hPhiSel[nj][ic] = new TH1F(Form("hPhiSel%d_%d",nj,ic),Form("Numerator  phi  for algorithm %s cent %d",calgo[nj],ic),20,-pi,pi);
      hPhiSel[nj][ic]->Sumw2();

      hFakePtAll [nj][ic] = new TH1F(Form("hFakePtAll%d_%d",nj,ic),Form("Fake Denominator pT for algorithm %s cent %d",calgo[nj],ic),500,0,1000);
      hFakePtAll [nj][ic]->Sumw2();
      hFakeEtaAll[nj][ic] = new TH1F(Form("hFakeEtaAll%d_%d",nj,ic),Form("Fake Denominator eta  for algorithm %s cent %d",calgo[nj],ic),20,-ketacut,ketacut);
      hFakeEtaAll[nj][ic]->Sumw2();
      hFakePhiAll[nj][ic] = new TH1F(Form("hFakePhiAll%d_%d",nj,ic),Form("Fake Denominator  phi  for algorithm %s cent %d",calgo[nj],ic),20,-pi,pi);
      hFakePhiAll[nj][ic]->Sumw2();

      hFakePtSel [nj][ic] = new TH1F(Form("hFakePtSel%d_%d",nj,ic),Form("Fake Numerator pT for algorithm %s cent %d",calgo[nj],ic),500,0,1000);
      hFakePtSel [nj][ic]->Sumw2();
      hFakeEtaSel[nj][ic] = new TH1F(Form("hFakeEtaSel%d_%d",nj,ic),Form("Fake Numerator eta  for algorithm %s cent %d",calgo[nj],ic),20,-ketacut,ketacut);
      hFakeEtaSel[nj][ic]->Sumw2();
      hFakePhiSel[nj][ic] = new TH1F(Form("hFakePhiSel%d_%d",nj,ic),Form("Fake Numerator  phi  for algorithm %s cent %d",calgo[nj],ic),20,-pi,pi);
      hFakePhiSel[nj][ic]->Sumw2();

    }//! ic
  }//! nj
  std::cout<<"Initialized the histograms " <<std::endl;
  ///////////////////////////////////////////////////////////////////////////////////////// 


  TTree *tr_in=0;
  std::vector<double> vJets;
  Long64_t nbytes=0;

  double nevt [kAlgos][ncent]={{0}};
  double njets[kAlgos][ncent]={{0}};

  for(int nj=0; nj<knj; nj++){

    tr_in = (TTree*)fin->Get(Form("%sJetAnalyzer/t",calgo[nj]));
    Long64_t nentries = tr_in->GetEntries();
    std::cout<<Form("# of entries in TTree for %s %s : ",calgo[nj],ksp)<<nentries<<std::endl;
    std::cout<<std::endl;
    
    //Declaration of leaves types
    int hiBin;
    int   nref;
    float pthat;
    float weight;
    float corrpt[1000];
    float jtpt[1000];
    float rawpt[1000];
    float jteta[1000];
    float jtphi[1000];
    float refpt[1000];
    float refeta[1000];
    float refphi[1000];
    float refdrjt[1000];
    float chargedMax[1000];
    float chargedSum[1000];
    float photonSum [1000];
    float neutralSum[1000];
    float refparton_pt[1000];
    int subid[1000];

    tr_in->SetBranchAddress("hiBin",&hiBin);
    tr_in->SetBranchAddress("nref",&nref);
    tr_in->SetBranchAddress("pthat",&pthat);
    tr_in->SetBranchAddress("weight",&weight);
    tr_in->SetBranchAddress("rawpt",rawpt);
    tr_in->SetBranchAddress("corrpt",corrpt);
    tr_in->SetBranchAddress("jtpt",jtpt);  //! this is also raw pt
    tr_in->SetBranchAddress("jteta",jteta);
    tr_in->SetBranchAddress("jtphi",jtphi);
    tr_in->SetBranchAddress("refpt",refpt);
    tr_in->SetBranchAddress("refphi",refphi);
    tr_in->SetBranchAddress("refeta",refeta);
    tr_in->SetBranchAddress("refdrjt",refdrjt);
    tr_in->SetBranchAddress("refparton_pt",refparton_pt);
    tr_in->SetBranchAddress("subid",subid);
    if(strcmp(ksp,"pp")!=0){
      tr_in->SetBranchAddress("chargedMax",chargedMax);
      tr_in->SetBranchAddress("chargedSum",chargedSum);
      tr_in->SetBranchAddress("photonSum",photonSum);
      tr_in->SetBranchAddress("neutralSum",neutralSum);
    }

    //! Load the jet energy correction factors on fly
    string L2Name, L3Name;
    JetCorrectorParameters *parHI_l2=0, *parHI_l3=0;
    vector<JetCorrectorParameters> vpar_HI;
    FactorizedJetCorrector *_JEC_HI = new FactorizedJetCorrector(vpar_HI);//JR                                                                         

    
     L2Name = "/net/hisrv0001/home/pawan/Validation/Track8_Jet19/combinePtHatBins/ppJEC2014/JECv14/HI_PythiaZ2_2760GeV_5316_v14_L2Relative_"+corrFileName[nj]+"_offline.txt";
     cout<<"**** ++++++  L2Name "<<L2Name<<endl;
     L3Name = "/net/hisrv0001/home/pawan/Validation/Track8_Jet19/combinePtHatBins/ppJEC2014/JECv14/HI_PythiaZ2_2760GeV_5316_v14_L3Absolute_"+corrFileName[nj]+"_offline.txt";
     cout<<"**** ++++++  L3Name "<<L3Name<<endl;
    
//     L2Name = "/net/hisrv0001/home/pawan/Validation/Track8_Jet19/CMSSW_5_3_16/src/JEC_base/pp2014/txtfiles/JEC_STARTHI53_LV1_Track8_Jet22_dijet_L2Relative_"+corrFileName[nj]+".txt";
//     cout<<"**** ++++++  L2Name "<<L2Name<<endl;
//     L3Name = "/net/hisrv0001/home/pawan/Validation/Track8_Jet19/CMSSW_5_3_16/src/JEC_base/pp2014/txtfiles/JEC_STARTHI53_LV1_Track8_Jet22_dijet_L3Absolute_"+corrFileName[nj]+".txt";
//     cout<<"**** ++++++  L3Name "<<L3Name<<endl;

//!  Improved at Low pT
//      L2Name = "/net/hisrv0001/home/pawan/Validation/Track8_Jet19/CMSSW_5_3_16/src/JEC_base/pp2014/JEC_STARTHI53_LV1_Track8_Jet22_dijet_L2Relative_"+corrFileName[nj]+".txt";
//      cout<<"**** ++++++  L2Name "<<L2Name<<endl;
//      L3Name = "/net/hisrv0001/home/pawan/Validation/Track8_Jet19/CMSSW_5_3_16/src/JEC_base/pp2014/JEC_STARTHI53_LV1_Track8_Jet22_dijet_L3Absolute_"+corrFileName[nj]+".txt";
//      cout<<"**** ++++++  L3Name "<<L3Name<<endl;


    parHI_l2 = new JetCorrectorParameters(L2Name.c_str());
    parHI_l3 = new JetCorrectorParameters(L3Name.c_str());
    
    vpar_HI.push_back(*parHI_l2);
    vpar_HI.push_back(*parHI_l3);
    _JEC_HI = new FactorizedJetCorrector(vpar_HI);
    
    
    Int_t iEvent=0;     
    for (Long64_t i=0; i<nentries;i++) {
      nbytes += tr_in->GetEntry(i);
      //return 0;

      int icent=-1;
      double wcen=1;
      double wvz=1;
      //! weight  for the merging of the samples for different pT hat bins
      double wxs = weight;
      
      if(strcmp(ksp,"pp")==0){
	icent=0;
      } else {
	icent = GetCentBin(hiBin);
      }
      //! MinBias
      //icent=0;
      
      if(i%10000==0)std::cout<<" ********** Event # " <<i<<"\t weight  : "<<weight<<"\t pthat : "<<pthat<<std::endl;
      hBin->Fill(hiBin,wxs*wcen*wvz);

      //! Centrality
      if(icent<0 || icent>=ncent ){
	//std::cout<<" Something wrong !!!!!  : "  << icent << "\t hiBin : "<< hiBin << " \t event " << i << std::endl;
	continue;
      }
      //! xsec-weight
      hpthat->Fill(pthat,wxs*wcen*wvz);
      nevt[nj][icent]++;

//       //! Jet energy scale comparison with data      
//       vJets.clear();
//       for(int igen=0; igen<nref; igen++){
// 	int gj = igen;
	
// 	if(subid[gj] != 0 || refparton_pt[gj]==-999)continue;
	
// 	_JEC_HI->setJetEta(jteta[gj]);
// 	_JEC_HI->setJetPt (rawpt[gj]);
	
// 	//float recopt  = corrpt[gj];
// 	float recopt  = rawpt[gj]*_JEC_HI->getCorrection();  //! correction with JECv14
// 	if(recopt<80 && fabs(jteta[gj])>1.4 && fabs(refdrjt[gj])>kDelR[nj])continue;	
// 	vJets.push_back(recopt);
//       }
//       if(vJets.size()>=2){
// 	std::sort(vJets.begin(),vJets.end());
// 	double B=-9999;
// 	double rn1 = gRandom->Rndm();
// 	double rn2 = gRandom->Rndm();
// 	double ptdij = (vJets[0] + vJets[1])/2.;
// 	if(rn1 > rn2){
// 	  B = (vJets[0] - vJets[1])/(vJets[0] + vJets[1]);
// 	}else{
// 	  B = (vJets[1] - vJets[0])/(vJets[1] + vJets[0]);
// 	}
// 	if(B!=-9999)hJetEnergyScale[nj][icent]->Fill(ptdij,B,wxs*wcen*wvz);
//       }


      //! Gen matched jets loop
      for(int igen=0; igen<nref; igen++){

	int gj = igen;

	if(fabs(refeta[gj]) > ketacut)continue;

	_JEC_HI->setJetEta(jteta[gj]);
	_JEC_HI->setJetPt (rawpt[gj]);


	float recopt  = rawpt[gj]*_JEC_HI->getCorrection();  //! correction with JECv14
	float recoeta = jteta[gj];
	float recophi = jtphi[gj];
	float delr    = refdrjt[gj];
	double jetbkgd = (chargedSum[gj] + photonSum[gj] + neutralSum[gj]) - rawpt[gj];

	
	if(subid[gj] == 0 && fabs(refeta[gj]) < ketacut && refpt[gj]>kptgencut){
	  hdRAll[nj][icent]->Fill(refpt[gj],wxs*wcen*wvz);
	  if(fabs(delr)<kdRcut)hdRSel[nj][icent]->Fill(refpt[gj],wxs*wcen*wvz);
	}


	//! Fake 
	if(fabs(delr)<kdRcut && fabs(refeta[gj]) < ketacut){
	  
	  hFakePtAll [nj][icent]->Fill(recopt,wxs*wcen*wvz);
	  hFakeEtaAll[nj][icent]->Fill(recoeta,wxs*wcen*wvz);
	  hFakePhiAll[nj][icent]->Fill(recophi,wxs*wcen*wvz);

	  if( refpt[gj] < 0 ){
	    hFakePtSel [nj][icent]->Fill(recopt,wxs*wcen*wvz);
	    hFakeEtaSel[nj][icent]->Fill(recoeta,wxs*wcen*wvz);
	    hFakePhiSel[nj][icent]->Fill(recophi,wxs*wcen*wvz);
	  }
	  if(strcmp(ksp,"pp")!=0){
	    hjbkgComb[nj][icent]->Fill(recopt,jetbkgd,wxs*wcen*wvz);
	  }
	}

	//! Reconstruction 
	if(fabs(refeta[gj]) < ketacut && refpt[gj]>kptgencut){

	  //! Denominator for reconstruction efficiency
	  hPtAll [nj][icent]->Fill(refpt[gj],wxs*wcen*wvz);
	  hEtaAll[nj][icent]->Fill(refeta[gj],wxs*wcen*wvz);
	  hPhiAll[nj][icent]->Fill(refphi[gj],wxs*wcen*wvz);

	  if(fabs(delr)<kdRcut){
	    //! Numerator for reconstrunction efficiency
	    hPtSel [nj][icent]->Fill(refpt[gj],wxs*wcen*wvz);
	    hEtaSel[nj][icent]->Fill(refeta[gj],wxs*wcen*wvz);
	    hPhiSel[nj][icent]->Fill(refphi[gj],wxs*wcen*wvz);
	  }
	}

	if(subid[gj] != 0)continue;

	//float recopt  = corrpt[gj];
	if(recopt<kptrecocut || refpt[gj]<kptgencut || refpt[gj]==0 || fabs(recoeta)>ketacut || fabs(delr)>kdRcut)continue;	
	
	if(strcmp(ksp,"pp")!=0){
	  hjbkgComb_match[nj][icent]->Fill(recopt,jetbkgd,wxs*wcen*wvz);
	}

	
	hgenjrecoj [nj][icent]->Fill(refpt[gj],recopt,wxs*wcen*wvz);
	hgenjrawj  [nj][icent]->Fill(refpt[gj],rawpt[gj],wxs*wcen*wvz);
	
	hratiocorrrefpt_genm[nj][icent]->Fill(refpt[gj],recopt/refpt[gj],wxs*wcen*wvz);
	hrecogen[nj][icent]->Fill(refpt[gj],recopt/refpt[gj],wxs*wcen*wvz);
	hrecoraw[nj][icent]->Fill(recopt,recopt/rawpt[gj],wxs*wcen*wvz);
	hrawgen [nj][icent]->Fill(refpt[gj],rawpt[gj]/refpt[gj],wxs*wcen*wvz);


	njets[nj][icent]++;
	
	
	int ieta=-1;
	if(fabs(recoeta)<1.3)ieta=0; //! barrel region
	else ieta=1; //! HCAL region
	
	
	//! Response in eta
	hratiocorrrefpt_eta[nj][icent][ieta]->Fill(refpt[gj],recopt/refpt[gj],wxs*wcen*wvz);
	hratiorawrefpt_eta [nj][icent][ieta]->Fill(refpt[gj],rawpt[gj]/refpt[gj],wxs*wcen*wvz);
	
	//! Ratio of recopT / rawpT
	hrajptrpt[nj][icent] ->Fill(recopt,recopt/rawpt[gj],wxs*wcen*wvz);
	
	//! Jet eta, phi, pt, eta-pt, eta-phi and pt-phi
	hjeteta   [nj][icent]->Fill(jteta[gj],wxs*wcen*wvz);
	hjetphi   [nj][icent]->Fill(jtphi[gj],wxs*wcen*wvz);
	hjetpteta [nj][icent]->Fill(recopt,jteta[gj],wxs*wcen*wvz);
	hjetptphi [nj][icent]->Fill(jtpt[gj],jtphi[gj],wxs*wcen*wvz);
	hjetetaphi[nj][icent]->Fill(jteta[gj],jtphi[gj],wxs*wcen*wvz);
	
	hgenpt_genm [nj][icent]->Fill(refpt[gj],wxs*wcen*wvz);
	hrecopt_genm[nj][icent]->Fill(recopt,wxs*wcen*wvz);	  
	hrawpt_genm [nj][icent]->Fill(rawpt[gj],wxs*wcen*wvz);
	
	//! Very fine bin in ref pt used for response
	hrescrpt_genm [nj][icent]->Fill(refpt[gj],recopt/refpt[gj],wxs*wcen*wvz);
	hresrrpt_genm [nj][icent]->Fill(refpt[gj],rawpt[gj]/refpt[gj],wxs*wcen*wvz);
	hresrcrpt_genm[nj][icent]->Fill(recopt,recopt/rawpt[gj],wxs*wcen*wvz);
	
	//! Very fine bin in ref eta used for response
	hrescreta_genm [nj][icent]->Fill(refeta[gj],refpt[gj],recopt/refpt[gj],wxs*wcen*wvz);
	hresrreta_genm [nj][icent]->Fill(refeta[gj],refpt[gj],rawpt[gj]/refpt[gj],wxs*wcen*wvz);
	hresrcreta_genm[nj][icent]->Fill(recoeta,refpt[gj],recopt/rawpt[gj],wxs*wcen*wvz);

	//! Very fine bin in ref phi used for response
	hrescrphi_genm [nj][icent]->Fill(refphi[gj],refpt[gj],recopt/refpt[gj],wxs*wcen*wvz);
	hresrrphi_genm [nj][icent]->Fill(refphi[gj],refpt[gj],rawpt[gj]/refpt[gj],wxs*wcen*wvz);
	hresrcrphi_genm[nj][icent]->Fill(recophi,refpt[gj],recopt/rawpt[gj],wxs*wcen*wvz);
	
	//! Response in different eta and phi bins
	int etabin = GetEtaBin(fabs(refeta[gj]));
	int phibin = GetPhiBin(refphi[gj]);
	
	//! Response in eta and phi bins
	if(etabin >= 0 && etabin<maxe){
	  hpteta[nj][icent][etabin]->Fill(refpt[gj],recopt/refpt[gj],wxs*wcen*wvz);
	}
	if(phibin >= 0 && phibin<maxph){
	  hptphi[nj][icent][phibin]->Fill(refpt[gj],recopt/refpt[gj],wxs*wcen*wvz);
	}

      }//! igen loop


      iEvent++;
      //std::cout<<"Completed event #  "<<ievt<<std::endl; 

    }//! event loop ends

    delete parHI_l2;
    delete parHI_l3;
    delete _JEC_HI;

  }//! jet loop

  std::cout<<std::endl;
  for(int nj=0;nj<knj;nj++){
    std::cout<<calgo[nj] << std::endl;
    if(strcmp(ksp,"pp")==0){
      std::cout<<"\t # of events : "<< " pp " << "  " << nevt[nj][0] << "  " << njets[nj][0] << std::endl;      
    }
    else{
      for(int ic=0;ic<ncent;ic++){
	std::cout<<"\t # of events : "<< ccent[ic] << "  " << nevt[nj][ic] << "  " << njets[nj][ic] << std::endl;      
      }
    }
    std::cout<<std::endl;
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
  int ibin=-1;
  float min =0.0;
  ibin = (eta - min)*maxe/ketacut;
  return ibin;
}
int GetPhiBin(float phi)
{
  int ibin=-1;
  float min = -pi;
  ibin = (phi - min)*maxph/2*pi;
  return ibin;
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
double delphi(double phi1, double phi2)
{
  double dphi = phi2 - phi1;
  dphi = fabs(atan2(sin(dphi),cos(dphi)));
  return dphi;
}
int GetCentBin(int bin)
{
  //! PbPb
  int ibin=-1;
  //! centrality is defined as 2.5% bins of cross section
  //! in 0-39 bins

//   if(bin<2)ibin=0; //! 0-5% 
//   else if(bin>=2 && bin<4)ibin=1;   //! 5-10%
//   else if(bin>=4 && bin<12)ibin=2;  //! 10-30%
//   else if(bin>=12&& bin<20)ibin=3;  //! 30-50% 
//   else if(bin>=20&& bin<28)ibin=4;  //! 50-70% 
//   else if(bin>=28&& bin<36)ibin=5;  //! 70-90%

  //! centrality bins as 0.5% 0-200
  if(bin<10)ibin=0; //! 0-5%
  else if(bin>=10  && bin<20)ibin=1;   //! 5-10%
  else if(bin>=20  && bin<60)ibin=2;   //! 10-30%
  else if(bin>=60  && bin<100)ibin=3;  //! 30-50%
  else if(bin>=100 && bin<140)ibin=4;  //! 50-70%
  else if(bin>=140 && bin<180)ibin=5;  //! 70-90%   
  return ibin;
}
