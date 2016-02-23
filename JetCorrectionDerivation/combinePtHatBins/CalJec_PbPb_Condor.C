
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
#include <fstream>
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


// struct Jet{
//   double frefpt;
//   double frawpt;
//   double fcorpt;
//   double fparpt;
// };
// bool compare_pt (float pt1, float pt2) { return (pt1 < pt2); }



int GetCentBin(int /*hiBin*/);
int GetPhiBin(float /*phi*/);
int GetEtaBin(float /*eta*/);
int GetPtHatBin(float /*pt*/);
int GetPtBin(float /*pt*/);
int GetPtBinWide(float /*pt*/);
double delphi(double /*phi1*/, double /*phi2*/);
double GetXsecWt(float /*pthat*/);

//const int npthat=11;
const int npthat=11;
double pthatwt[npthat][4] ={
//!  pthat  effEv    effxec     wt  
     // {15   ,241077  ,0.49235   ,2.0422935e-06},
     // {30   ,176205  ,0.030482  ,1.7299168e-07},
     // {50   ,222893  ,0.003572  ,1.6026075e-08},
     // {80   ,244756  ,0.000425  ,1.7361781e-09},
     // {120  ,250177  ,5.873e-05 ,2.3475379e-10},
     // {170  ,196962  ,9.199e-06 ,4.6704440e-11},
     // {220  ,190772  ,2.2564e-06,1.1827732e-11},
     // {280  ,263915  ,6.336e-07 ,2.4007798e-12},
     // {370  ,232644  ,1.0884e-07,4.6783927e-13},
     // {460  ,192713  ,2.215e-08 ,1.1493776e-13},
     // {540  ,310286  ,1.001e-08 ,3.2260560e-14}


  //! Embedded sample v9
     // {15 ,280123,0.49235   ,1.75762076e-06},
     // {30 ,194098,0.030482  ,1.57044380e-07},
     // {50 ,242561,0.003572  ,1.47266048e-08},
     // {80 ,284872,0.000425  ,1.49168749e-09},
     // {120,288696,5.873e-05 ,2.03431984e-10},
     // {170,232095,9.199e-06 ,3.96346324e-11},
     // {220,244596,2.2564e-06,9.225007768e-12},
     // {280,313811,6.336e-07 ,2.019049682e-12},
     // {370,244899,1.0884e-07,4.444281112e-13},
     // {460,224013,2.215e-08 ,9.887819011e-14},
     // {540,350236,1.001e-08 ,2.858072842e-14}

  {15,161479,0.49235    ,3.049003e-06},
  {30,405739,0.030482   ,7.512711e-08},
  {50,259546,0.003572   ,1.376288e-08},
  {80,256506,0.00042494 ,1.656647e-09},
  {120,255111,5.873e-05 ,2.302135e-10},
  {170,220926,9.199e-06 ,4.163838e-11},
  {220,184475,2.2564e-06,1.223147e-11},
  {280,50928,6.336e-07  ,1.244109e-11},
  {370,195420,1.0884e-07,5.569542e-13},
  {460,206363,2.215e-08 ,1.073351e-13},
  {540,345207,1.001e-08 ,2.899709e-14}
};

const int ncen=10;
const char *centbin[ncen] = {"0-5%","5-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-100%"};

//! pt binning
double ptbins[] ={17, 22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 400, 550, 790, 1000};
const int nbins = sizeof(ptbins)/sizeof(double) - 1;

double pthatbins[] ={15, 30, 50, 80, 120, 170, 220, 280, 370, 460, 540, 1000};
const int kbins = sizeof(pthatbins)/sizeof(double) - 1;

//! constants
int iYear=2015;
const double pi=acos(-1.);
//const double pi2=2*pi -1;

const double kvzcut    =15.0;
const float  ketacut   =2.0;
const double krawptcut =10.0;
const double kptrecocut=0.0;
const double kptgencut =0.0;
const double kdRCut=0.30;

int rbins=150;
double rbinl=0.0;
double rbinh=2.0;

const double ptbins_wide[]={20,50,80,120,200,350,800};
const int nbins_wide = sizeof(ptbins_wide)/sizeof(double) - 1;

double etabins[] = {-3.000, -2.400, -2.000, -1.3000, 0.000, 1.300, 2.000, 2.400, 3.000};
const int neta = sizeof(etabins)/sizeof(double) - 1;

const double phibins[] = {-3.141,-2.700,-2.100,-1.500,-0.900,-0.300, 
			  0.300,0.900,1.500,2.100,2.700,3.141
};
const int nphi = sizeof(phibins)/sizeof(double) - 1;

bool is_file(const char *fileName);
bool is_file(const char *fileName)
{
  std::ifstream infile(fileName,ios::in);
  return infile.good();
}

TStopwatch timer;

int CalJec_PbPb_Condor(const char *runalgo="akVs1PF", 
		       string inname="Pythia8_Dijet_pp_TuneCUETP8M1_Hydjet_MinBias_5020GeV_HiForest_PrivMC_HiForestPbPbJECv4LocalVoronoi_merged_forest_50_0.root",
		       string outfile="test.root", 
		       const char *ksp="PbPb")
{

  timer.Start();

  int nj=0;
  string corrFileName = "";
  if( strcmp(runalgo,"akPu1PF") == 0){corrFileName="AK1PF";nj=0;}
  else if( strcmp(runalgo,"akPu2PF") == 0){corrFileName="AK2PF";nj=1;}
  else if( strcmp(runalgo,"akPu3PF") == 0){corrFileName="AK3PF";nj=2;}
  else if( strcmp(runalgo,"akPu4PF") == 0){corrFileName="AK4PF";nj=3;}
  else if( strcmp(runalgo,"akPu5PF") == 0){corrFileName="AK5PF";nj=4;}
  else if( strcmp(runalgo,"akPu6PF") == 0){corrFileName="AK6PF";nj=5;}
  else if( strcmp(runalgo,"akPu1Calo") == 0){corrFileName="AK1Calo";nj=6;}
  else if( strcmp(runalgo,"akPu2Calo") == 0){corrFileName="AK2Calo";nj=7;}
  else if( strcmp(runalgo,"akPu3Calo") == 0){corrFileName="AK3Calo";nj=8;}
  else if( strcmp(runalgo,"akPu4Calo") == 0){corrFileName="AK4Calo";nj=9;}
  else if( strcmp(runalgo,"akPu5Calo") == 0){corrFileName="AK5Calo";nj=10;}
  else if( strcmp(runalgo,"akPu6Calo") == 0){corrFileName="AK6Calo";nj=11;}
  else if( strcmp(runalgo,"akVs1PF") == 0){corrFileName="AK1PF";nj=12;}
  else if( strcmp(runalgo,"akVs2PF") == 0){corrFileName="AK2PF";nj=13;}
  else if( strcmp(runalgo,"akVs3PF") == 0){corrFileName="AK3PF";nj=14;}
  else if( strcmp(runalgo,"akVs4PF") == 0){corrFileName="AK4PF";nj=15;}
  else if( strcmp(runalgo,"akVs5PF") == 0){corrFileName="AK5PF";nj=16;}
  else if( strcmp(runalgo,"akVs6PF") == 0){corrFileName="AK6PF";nj=17;}
  else if( strcmp(runalgo,"akVs1Calo") == 0){corrFileName="AK1Calo";nj=18;}
  else if( strcmp(runalgo,"akVs2Calo") == 0){corrFileName="AK2Calo";nj=19;}
  else if( strcmp(runalgo,"akVs3Calo") == 0){corrFileName="AK3Calo";nj=20;}
  else if( strcmp(runalgo,"akVs4Calo") == 0){corrFileName="AK4Calo";nj=21;}
  else if( strcmp(runalgo,"akVs5Calo") == 0){corrFileName="AK5Calo";nj=22;}
  else if( strcmp(runalgo,"akVs6Calo") == 0){corrFileName="AK6Calo";nj=23;}

  std::string basedir ="";
  //std::string dirname ="JEC_HIReco/";
  std::string dirname ="JEC_HIReco_HcalRespCorrs_v4_00_mc/";
  std::string jecera  ="JEC_pp_PYTHIA_TuneCUETP8M1_5020GeV_patch3";



  //std::string inDir=Form("/mnt/hadoop/cms/store/user/pawan/Pythia8_Dijet%d_pp_TuneCUETP8M1_Hydjet_MinBias_5020GeV_HiForest_PrivMC_TAGHiForestPbPbJECv9/",kpthat);
  std::string inDir="";
  //std::string inDir="/mnt/hadoop/cms/store/user/dgulhan/mergedForest/Pythia8_Dijet_pp_TuneCUETP8M1_Hydjet_MinBias_5020GeV_HiForest_PrivMC_HiForestPbPbJECv4LocalVoronoi/";
  TFile *fin = new TFile((inDir+inname).c_str(),"r");

  TFile *fout = new TFile(outfile.c_str(),"RECREATE");
  
  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"**************************************************** "<<std::endl;
  std::cout<<Form("Analyzeforest_jec : %s  %s",ksp,dirname.c_str())<<std::endl;
  std::cout<<Form("eta cut for %0.3f ",ketacut)<<std::endl;
  std::cout<<"**************************************************** "<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;



  //! Before running check all the txt files are available
  std::string L2Name = basedir+dirname+jecera+"_L2Relative_"+corrFileName+".txt";
  if(!is_file(L2Name.c_str())){
    cout<<"**** ++++++  L2Name does not exists  "<<L2Name<<endl;
    return  2;
  }
  std::string L3Name = basedir+dirname+jecera+"_L3Absolute_"+corrFileName+".txt";
  if(!is_file(L3Name.c_str())){
    cout<<"**** ++++++  L3Name does not exists  "<<L3Name<<endl;
    return 3;
  }
  cout <<"  Files are there " << endl;
  //! 
  //! Define histograms here
  //   TH1::SetDefaultSumw2();
  //   TH2::SetDefaultSumw2();
  //   TH3::SetDefaultSumw2();
  //  TProfile::SetDefaultSumw2();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  //! Gen matched jets
  // TH1D *hgenpt_genm[ncen] , *hrecopt_genm[ncen], *hrawpt_genm[ncen];
  // TH1D *hgeneta[ncen]     , *hgenphi[ncen]     ;
  // TH1D *hjeteta[ncen]     , *hjetphi[ncen]     ;

  //! Resposnse
  TH1D *hrescrpt_genm [ncen][nbins], *hresrrpt_genm [ncen][nbins], *hratio [ncen][nbins];
  //TH1D *hrescrpt_genm_eta [ncen][nbins][neta], *hresrrpt_genm_eta [ncen][nbins][neta], *hratio_eta [ncen][nbins][neta];
  //TH1D *hrescrpt_genm_phi [ncen][nbins][nphi], *hresrrpt_genm_phi [ncen][nbins][nphi], *hratio_phi [ncen][nbins][nphi];
  TH1D *hrescrpt_wide_genm_eta [ncen][nbins_wide][neta], *hresrrpt_wide_genm_eta [ncen][nbins_wide][neta], *hratio_wide_eta [ncen][nbins_wide][neta];
  TH1D *hrescrpt_wide_genm_phi [ncen][nbins_wide][nphi], *hresrrpt_wide_genm_phi [ncen][nbins_wide][nphi], *hratio_wide_phi [ncen][nbins_wide][nphi];

  // TH1D *havjetpt_genm[ncen][nbins], *havrawpt_genm[ncen][nbins];
  // TProfile *hcorr_rawpt[ncen], *hcorr_genpt[ncen], *hcorr_recopt[ncen];
  // TProfile *hl3corr_rawpt[ncen], *hl3corr_genpt[ncen], *hl3corr_recopt[ncen];
  // TProfile *hl2corr_rawpt[ncen], *hl2corr_genpt[ncen], *hl2corr_recopt[ncen];
  // TProfile *hrawpt_genpt[ncen], *hrecopt_genpt[ncen];

  TH1D *hBin;
  TH1D *hpthat[ncen];


  fout->mkdir(Form("%sJetAnalyzer",runalgo));
  fout->cd(Form("%sJetAnalyzer",runalgo));

  hBin = new TH1D(Form("hBin_%d",nj),"Centrality bin",200,-0.5,200-0.5);
  hBin->Sumw2();

  for(int ic=0; ic<ncen; ic++){
    hpthat[ic]  = new TH1D(Form("hpthat%d_%d",nj,ic),Form("%s pt-hat distribution",centbin[ic]),500,0.,1000.);
    hpthat[ic]->Sumw2();
    
    // hgenpt_genm[ic]  = new TH1D(Form("hgenpt_genm%d_%d",nj,ic),Form("Gen matched gen p_{T} distribution jet %s %s",runalgo,centbin[ic]),500,0.,1000.);
    // hgenpt_genm[ic]->Sumw2();
    // hrecopt_genm[ic] = new TH1D(Form("hrecopt_genm%d_%d",nj,ic),Form("Gen matched reco p_{T} distribution jet %s %s",runalgo,centbin[ic]),500,0.,1000.);
    // hrecopt_genm[ic]->Sumw2();
    // hrawpt_genm[ic]  = new TH1D(Form("hrawpt_genm%d_%d",nj,ic),Form("Gen matched raw p_{T} distribution jet %s %s",runalgo,centbin[ic]),500,0.,1000.);
    // hrawpt_genm[ic]->Sumw2();

    // hjeteta[ic] = new TH1D(Form("hjeteta%d_%d",nj,ic),Form("jet eta distribution jet %s %s",runalgo,centbin[ic]),72,-ketacut,ketacut);
    // hjeteta[ic]->Sumw2();
    // hjetphi[ic] = new TH1D(Form("hjetphi%d_%d",nj,ic),Form("jet phi distribution jet %s %s",runalgo,centbin[ic]),72,-pi,pi);
    // hjetphi[ic]->Sumw2();
    
    // hgeneta[ic] = new TH1D(Form("hgeneta%d_%d",nj,ic),Form("gen jet eta distribution jet %s %s",runalgo,centbin[ic]),72,-ketacut,ketacut);
    // hgeneta[ic]->Sumw2();
    // hgenphi[ic] = new TH1D(Form("hgenphi%d_%d",nj,ic),Form("gen jet phi distribution jet %s %s",runalgo,centbin[ic]),72,-pi,pi);
    // hgenphi[ic]->Sumw2();
  

    // hcorr_genpt[ic]  = new TProfile(Form("hcorr_genpt%d_%d",nj,ic),Form("Corr Fac gen p_{T} %s %s",runalgo,centbin[ic]),nbins,ptbins);
    // hcorr_genpt[ic]->Sumw2();
    // hcorr_rawpt[ic]  = new TProfile(Form("hcorr_rawpt%d_%d",nj,ic),Form("Corr Fac raw p_{T} %s %s",runalgo,centbin[ic]),nbins,ptbins);
    // hcorr_rawpt[ic]->Sumw2();
    // hcorr_recopt[ic] = new TProfile(Form("hcorr_recopt%d_%d",nj,ic),Form("Corr Fac reco p_{T} %s %s",runalgo,centbin[ic]),nbins,ptbins);
    // hcorr_recopt[ic]->Sumw2();
  
    // hl2corr_genpt[ic]  = new TProfile(Form("hl2corr_genpt%d_%d",nj,ic),Form("l2Corr Fac gen p_{T} %s %s",runalgo,centbin[ic]),nbins,ptbins);
    // hl2corr_genpt[ic]->Sumw2();
    // hl2corr_rawpt[ic]  = new TProfile(Form("hl2corr_rawpt%d_%d",nj,ic),Form("l2Corr Fac raw p_{T} %s %s",runalgo,centbin[ic]),nbins,ptbins);
    // hl2corr_rawpt[ic]->Sumw2();
    // hl2corr_recopt[ic]  = new TProfile(Form("hl2corr_recopt%d_%d",nj,ic),Form("l2Corr Fac reco p_{T} %s %s",runalgo,centbin[ic]),nbins,ptbins);
    // hl2corr_recopt[ic]->Sumw2();
  
    // hl3corr_genpt[ic]  = new TProfile(Form("hl3corr_genpt%d_%d",nj,ic),Form("l3Corr Fac gen p_{T} %s %s",runalgo,centbin[ic]),nbins,ptbins);
    // hl3corr_genpt[ic] ->Sumw2();
    // hl3corr_rawpt[ic]  = new TProfile(Form("hl3corr_rawpt%d_%d",nj,ic),Form("l3Corr Fac raw p_{T} %s %s",runalgo,centbin[ic]),nbins,ptbins);
    // hl3corr_rawpt[ic] ->Sumw2();
    // hl3corr_recopt[ic]  = new TProfile(Form("hl3corr_recopt%d_%d",nj,ic),Form("l3Corr Fac reco p_{T} %s %s",runalgo,centbin[ic]),nbins,ptbins);
    // hl3corr_recopt[ic] ->Sumw2();
  
    // hrawpt_genpt[ic]  = new TProfile(Form("hrawpt_genpt%d_%d",nj,ic),Form("reco p_{T} gen p_{T} %s %s",runalgo,centbin[ic]),nbins,ptbins);
    // hrawpt_genpt[ic]->Sumw2();
    // hrecopt_genpt[ic]  = new TProfile(Form("hrecopt_genpt%d_%d",nj,ic),Form("reco p_{T} gen p_{T} %s %s",runalgo,centbin[ic]),nbins,ptbins);
    // hrecopt_genpt[ic]->Sumw2();

    for(int ip=0;ip<nbins;ip++){
    
      // havjetpt_genm [ic][ip]= new TH1D(Form("havjetpt_genm%d_%d_%d",nj,ic,ip),Form("(<Reco>) jet p_{T} %s %0.0f < p_{T}^{REF} < %0.0f",runalgo,ptbins[ip],ptbins[ip+1]),500,0.,1000.);
      // havjetpt_genm [ic][ip]->Sumw2();
      // havrawpt_genm [ic][ip]= new TH1D(Form("havrawpt_genm%d_%d_%d",nj,ic,ip),Form("(<Raw>) jet p_{T} %s %0.0f < p_{T}^{REF} < %0.0f",runalgo,ptbins[ip],ptbins[ip+1]),500,0.,1000.);
      // havrawpt_genm [ic][ip]->Sumw2();
    
      //! Gen matched Response and resolution
      hrescrpt_genm [ic][ip]= new TH1D(Form("hrescrpt_genm%d_%d_%d",nj,ic,ip),Form("(Reco/Gen) jet p_{T} %s %0.0f < p_{T}^{REF} < %0.0f",runalgo,ptbins[ip],ptbins[ip+1]),rbins,rbinl,rbinh);
      hrescrpt_genm [ic][ip]->Sumw2();
      hresrrpt_genm [ic][ip]= new TH1D(Form("hresrrpt_genm%d_%d_%d",nj,ic,ip),Form("(Raw/Gen) jet p_{T} %s %0.0f < p_{T}^{REF} < %0.0f",runalgo,ptbins[ip],ptbins[ip+1]),rbins,rbinl,rbinh);
      hresrrpt_genm [ic][ip]->Sumw2();
      hratio[ic][ip]= new TH1D(Form("hratio%d_%d_%d",nj,ic,ip),Form("(Reco/Raw) jet p_{T} %s %0.0f < p_{T}^{REF} < %0.0f",runalgo,ptbins[ip],ptbins[ip+1]),rbins,rbinl,rbinh);
      hratio[ic][ip]->Sumw2();
    
      // for(int ie=0;ie<neta;ie++){      
      // 	hrescrpt_genm_eta [ic][ip][ie]= new TH1D(Form("hrescrpt_genm_eta%d_%d_%d_%d",nj,ic,ip,ie),Form("(Reco/Gen) jet p_{T} %s %d %d",runalgo,ip,ie),rbins,rbinl,rbinh);
      // 	hrescrpt_genm_eta [ic][ip][ie]->Sumw2();
      // 	hresrrpt_genm_eta [ic][ip][ie]= new TH1D(Form("hresrrpt_genm_eta%d_%d_%d_%d",nj,ic,ip,ie),Form("(Raw/Gen) jet p_{T} %s %d %d",runalgo,ip,ie),rbins,rbinl,rbinh);
      // 	hresrrpt_genm_eta [ic][ip][ie]->Sumw2();
      // 	hratio_eta[ic][ip][ie]= new TH1D(Form("hratio_eta%d_%d_%d_%d",nj,ic,ip,ie),Form("(Reco/Raw) jet p_{T} %s %d %d",runalgo,ip,ie),rbins,rbinl,rbinh);
      // 	hratio_eta[ic][ip][ie]->Sumw2();
      // }
      
    //   for(int ij=0;ij<nphi;ij++){      
    // 	hrescrpt_genm_phi [ic][ip][ij]= new TH1D(Form("hrescrpt_genm_phi%d_%d_%d_%d",nj,ic,ip,ij),Form("(Reco/Gen) jet p_{T} %s %d %d",runalgo,ip,ij),rbins,rbinl,rbinh);
    // 	hrescrpt_genm_phi [ic][ip][ij]->Sumw2();
    // 	hresrrpt_genm_phi [ic][ip][ij]= new TH1D(Form("hresrrpt_genm_phi%d_%d_%d_%d",nj,ic,ip,ij),Form("(Raw/Gen) jet p_{T} %s %d %d",runalgo,ip,ij),rbins,rbinl,rbinh);
    // 	hresrrpt_genm_phi [ic][ip][ij]->Sumw2();
    // 	hratio_phi[ic][ip][ij]= new TH1D(Form("hratio_phi%d_%d_%d_%d",nj,ic,ip,ij),Form("(Reco/Raw) jet p_{T} %s %d %d",runalgo,ip,ij),rbins,rbinl,rbinh);
    // 	hratio_phi[ic][ip][ij]->Sumw2();
    //   }

    }

    
    //! coarse pt bin
    for(int ip=0;ip<nbins_wide;ip++){
      for(int ie=0;ie<neta;ie++){      
	hrescrpt_wide_genm_eta [ic][ip][ie]= new TH1D(Form("hrescrpt_wide_genm_eta%d_%d_%d_%d",nj,ic,ip,ie),Form("(Reco/Gen) jet p_{T} %s %d %d",runalgo,ip,ie),rbins,rbinl,rbinh);
	hrescrpt_wide_genm_eta [ic][ip][ie]->Sumw2();
	hresrrpt_wide_genm_eta [ic][ip][ie]= new TH1D(Form("hresrrpt_wide_genm_eta%d_%d_%d_%d",nj,ic,ip,ie),Form("(Raw/Gen) jet p_{T} %s %d %d",runalgo,ip,ie),rbins,rbinl,rbinh);
	hresrrpt_wide_genm_eta [ic][ip][ie]->Sumw2();
	hratio_wide_eta[ic][ip][ie]= new TH1D(Form("hratio_wide_eta%d_%d_%d_%d",nj,ic,ip,ie),Form("(Reco/Raw) jet p_{T} %s %d %d",runalgo,ip,ie),rbins,rbinl,rbinh);
	hratio_wide_eta[ic][ip][ie]->Sumw2();
      }
      for(int ij=0;ij<nphi;ij++){      
	hrescrpt_wide_genm_phi [ic][ip][ij]= new TH1D(Form("hrescrpt_wide_genm_phi%d_%d_%d_%d",nj,ic,ip,ij),Form("(Reco/Gen) jet p_{T} %s %d %d",runalgo,ip,ij),rbins,rbinl,rbinh);
	hrescrpt_wide_genm_phi [ic][ip][ij]->Sumw2();
	hresrrpt_wide_genm_phi [ic][ip][ij]= new TH1D(Form("hresrrpt_wide_genm_phi%d_%d_%d_%d",nj,ic,ip,ij),Form("(Raw/Gen) jet p_{T} %s %d %d",runalgo,ip,ij),rbins,rbinl,rbinh);
	hresrrpt_wide_genm_phi [ic][ip][ij]->Sumw2();
	hratio_wide_phi[ic][ip][ij]= new TH1D(Form("hratio_wide_phi%d_%d_%d_%d",nj,ic,ip,ij),Form("(Reco/Raw) jet p_{T} %s %d %d",runalgo,ip,ij),rbins,rbinl,rbinh);
	hratio_wide_phi[ic][ip][ij]->Sumw2();
      }
    }
  }

  fout->cd("../");
  std::cout<<"Initialized the histograms " <<std::endl;
  ///////////////////////////////////////////////////////////////////////////////////////// 


  TTree *tr_in=0;
  //std::vector<double> vJets;
  Long64_t nbytes=0;

  double nevt[ncen]={0};
  double njets[ncen]={0};
  
  tr_in = (TTree*)fin->Get(Form("%sJetAnalyzer/t",runalgo));
  Long64_t nentries = tr_in->GetEntries();
  std::cout<<Form("# of entries in TTree for %s %s : ",runalgo,ksp)<<nentries<<std::endl;
  
  
  //Declaration of leaves types
  TTree *tr_ev = 0;
  float vz;
  int hiBin;
  tr_ev = (TTree*)fin->Get("hiEvtAnalyzer/HiTree");
  tr_ev->SetBranchAddress("vz",&vz);
  tr_ev->SetBranchAddress("hiBin",&hiBin);
  tr_ev->SetBranchStatus("*",0,0);
  tr_ev->SetBranchStatus("hiBin",1);
  tr_ev->SetBranchStatus("vz",1);


  int ngen;
  int genmatchindex[1000];

  int   nref;
  float pthat;
  float jtpt[1000];
  float rawpt[1000];
  float jteta[1000];
  float jtphi[1000];
  float refpt[1000];
  float refeta[1000];
  float refphi[1000];
  float refdrjt[1000];
  //float refparton_pt[1000];
  int refparton_flavor[1000];
  int subid[1000];

  tr_in->SetBranchAddress("nref",&nref);
  tr_in->SetBranchAddress("pthat",&pthat);
  tr_in->SetBranchAddress("rawpt",rawpt);
  tr_in->SetBranchAddress("jtpt",jtpt);
  tr_in->SetBranchAddress("jteta",jteta);
  tr_in->SetBranchAddress("jtphi",jtphi);
  tr_in->SetBranchAddress("refpt",refpt);
  tr_in->SetBranchAddress("refphi",refphi);
  tr_in->SetBranchAddress("refeta",refeta);
  tr_in->SetBranchAddress("refdrjt",refdrjt);
  //tr_in->SetBranchAddress("refparton_pt",refparton_pt);
  tr_in->SetBranchAddress("refparton_flavor",refparton_flavor);
  tr_in->SetBranchAddress("subid",subid);
  tr_in->SetBranchAddress("ngen",&ngen);
  tr_in->SetBranchAddress("genmatchindex",genmatchindex);

  tr_in->SetBranchStatus("*",0,0);
  tr_in->SetBranchStatus("nref",1);
  tr_in->SetBranchStatus("pthat",1);
  tr_in->SetBranchStatus("rawpt",1);
  tr_in->SetBranchStatus("jtpt" ,1);
  tr_in->SetBranchStatus("jteta",1);
  tr_in->SetBranchStatus("jtphi",1);
  tr_in->SetBranchStatus("refpt",1);
  tr_in->SetBranchStatus("refphi",1);
  tr_in->SetBranchStatus("refeta",1);
  //tr_in->SetBranchStatus("refdrjt",1);
  tr_in->SetBranchStatus("subid",1);
  tr_in->SetBranchStatus("refparton_flavor",1);
  tr_in->SetBranchStatus("ngen",1);
  tr_in->SetBranchStatus("genmatchindex",1);  

  //! Add Friends to the TTree
  tr_in->AddFriend(tr_ev);

  //! Load the jet energy correction factors on fly
  JetCorrectorParameters *parHI_l2=0, *parHI_l3=0;
  vector<JetCorrectorParameters> vpar_HI;
  FactorizedJetCorrector *_JEC_HI = new FactorizedJetCorrector(vpar_HI);//JR
    
    
  parHI_l2 = new JetCorrectorParameters(L2Name.c_str());
  parHI_l3 = new JetCorrectorParameters(L3Name.c_str());

  //cout <<   L3Name.c_str() <<  endl;
    
  vpar_HI.push_back(*parHI_l2);
  vpar_HI.push_back(*parHI_l3);
  _JEC_HI = new FactorizedJetCorrector(vpar_HI);
    
    
  Int_t iEvent=0;     
  for (Long64_t i=0; i<nentries;i++) {
    //for (Long64_t i=0; i<50;i++) {
    nbytes += tr_in->GetEntry(i);

    if( fabs(vz) > kvzcut )continue;

    double wcen=1.0;
    double wvz =1.0;

    //! weight  for the merging of the samples for different pT hat bins
    double wxs = GetXsecWt(pthat);
    if( wxs < 0 )continue;

    if(i%10000==0)std::cout<<" ********** Event # " <<i<<"\t weight  : "<<wxs<<"\t pthat : "<<pthat<<"\t hiBin : " <<hiBin<<std::endl;
    //std::cout<<" ********** Event # " <<i<<"\t weight  : "<<wxs<<"\t pthat : "<<pthat<<"\t hiBin : " << hiBin << std::endl;
    hBin->Fill(hiBin,wxs*wcen*wvz);

    int iCent = GetCentBin(hiBin);
    if( iCent < 0 || iCent >= ncen )continue;

    //std::cout<<" ********** Event # " <<i<<"\t weight  : "<<weight<<"\t pthat : "<<pthat<<std::endl;

    //! xsec-weight
    hpthat[iCent]->Fill(pthat,wxs*wcen*wvz);
    nevt[iCent]++;

    //! Gen matched jets loop
    for(int igen=0; igen<ngen; igen++){

      //! Only gen matched jets
      if( genmatchindex[igen] < 0 )continue;

      int gj = igen;

      if( fabs(refeta[gj]) > ketacut || refpt[gj] < 0 || subid[gj]!=0 )continue;
	
      if( rawpt[gj] < krawptcut ) continue;

      //if( abs(refparton_flavor[gj]) > 21 )continue;

      float delr    = refdrjt[gj];

      if(delr > kdRCut)continue;

      //if( subid[gj]!=0 || fabs(refeta[gj]) > ketacut || rawpt[gj] < krawptcut  
      //  ||  refpt[gj] < 0 ||  refpt[gj] < kptgencut || rawpt[gj]/refpt[gj] > 1)continue;
	
      _JEC_HI->setJetEta(jteta[gj]);
      _JEC_HI->setJetPt (rawpt[gj]);
	
      float corrfac = _JEC_HI->getCorrection();

      std::vector <float> subcorr = _JEC_HI->getSubCorrections();
      float l2corr = subcorr[0];
      float l3corr = subcorr[1];
      float recopt  = rawpt[gj]*corrfac;  
      //float recopt_l2corr = rawpt[gj]*l2corr;
      //float recopt_l3corr = rawpt[gj]*l3corr;

      //std::cout << " \t" << runalgo << " " << pthat << " refpt : " << refpt[gj] <<  " rawpt :  " << rawpt[gj] << "\t recopt "<< recopt << " l3 corr "<< l3corr << " l2 corr :  "<< recopt_l2corr << " corr : " << corrfac << std::endl; 
	
      if( (abs(corrfac - recopt/rawpt[gj]) > 0.1) || corrfac > 5 ){
	std::cout <<" WARNING " << runalgo << " refpt : " << refpt[gj] <<  " rawpt :  " << rawpt[gj] << "\t recopt "<< recopt << " corr : " << corrfac << std::endl; 
	//continue;
      }
      njets[iCent]++;

      // hcorr_genpt[iCent]->Fill(refpt[gj],corrfac,wxs*wcen*wvz);
      // hl2corr_genpt[iCent]->Fill(refpt[gj],l2corr,wxs*wcen*wvz);
      // hl3corr_genpt[iCent]->Fill(refpt[gj],l3corr,wxs*wcen*wvz);

      // hcorr_rawpt[iCent]->Fill(rawpt[gj],corrfac,wxs*wcen*wvz);
      // hl2corr_rawpt[iCent]->Fill(rawpt[gj],l2corr,wxs*wcen*wvz);
      // hl3corr_rawpt[iCent]->Fill(rawpt[gj],l3corr,wxs*wcen*wvz);

      // hcorr_recopt[iCent]->Fill(recopt,corrfac,wxs*wcen*wvz);
      // hl2corr_recopt[iCent]->Fill(recopt,l2corr,wxs*wcen*wvz);
      // hl3corr_recopt[iCent]->Fill(recopt,l3corr,wxs*wcen*wvz);

      // hrawpt_genpt[iCent]->Fill(refpt[gj],rawpt[gj],wxs*wcen*wvz);
      // hrecopt_genpt[iCent]->Fill(refpt[gj],recopt,wxs*wcen*wvz);


      float recoeta = jteta[gj];
      float recophi = jtphi[gj];
	
      double resp_corr =  recopt    / refpt[gj];
      double resp_raw  =  rawpt[gj] / refpt[gj];
      double ratio     =  recopt    / rawpt[gj];

      int iphi = GetPhiBin(refphi[gj]);
      int ieta = GetEtaBin(refeta[gj]);
      int ipt  = GetPtBin (refpt[gj]);
      int ipt_wide = GetPtBinWide (refpt[gj]);

      // hgenpt_genm[iCent]->Fill(refpt[gj],wxs*wcen*wvz);
      // hrecopt_genm[iCent]->Fill(recopt,wxs*wcen*wvz);
      // hrawpt_genm[iCent]->Fill(rawpt[gj],wxs*wcen*wvz);

      // hjeteta[iCent]->Fill(recoeta,wxs*wcen*wvz);
      // hjetphi[iCent]->Fill(recophi,wxs*wcen*wvz);

      // hgeneta[iCent]->Fill(refeta[gj],wxs*wcen*wvz);
      // hgenphi[iCent]->Fill(refphi[gj],wxs*wcen*wvz);


      //! Response in pt and eta
      if( ipt>=0 ){
	// havjetpt_genm[iCent] [ipt]->Fill(recopt,wxs*wcen*wvz);
	// havrawpt_genm[iCent] [ipt]->Fill(rawpt[gj],wxs*wcen*wvz);

	hrescrpt_genm [iCent][ipt]->Fill(resp_corr,wxs*wcen*wvz);
	hresrrpt_genm [iCent][ipt]->Fill(resp_raw ,wxs*wcen*wvz);
	hratio        [iCent][ipt]->Fill(ratio    ,wxs*wcen*wvz);

	if( ieta>=0 ){
	  // hrescrpt_genm_eta [iCent][ipt][ieta]->Fill(resp_corr,wxs*wcen*wvz);
	  // hresrrpt_genm_eta [iCent][ipt][ieta]->Fill(resp_raw ,wxs*wcen*wvz);
	  // hratio_eta        [iCent][ipt][ieta]->Fill(ratio    ,wxs*wcen*wvz);
	  if( ipt_wide>=0 ){
	    hrescrpt_wide_genm_eta [iCent][ipt_wide][ieta]->Fill(resp_corr,wxs*wcen*wvz);
	    hresrrpt_wide_genm_eta [iCent][ipt_wide][ieta]->Fill(resp_raw ,wxs*wcen*wvz);
	    hratio_wide_eta        [iCent][ipt_wide][ieta]->Fill(ratio    ,wxs*wcen*wvz);
	  }
	}

	if( iphi>=0 ){
	  // hrescrpt_genm_phi [iCent][ipt][iphi]->Fill(resp_corr,wxs*wcen*wvz);
	  // hresrrpt_genm_phi [iCent][ipt][iphi]->Fill(resp_raw ,wxs*wcen*wvz);
	  // hratio_phi        [iCent][ipt][iphi]->Fill(ratio    ,wxs*wcen*wvz);
	  if( ipt_wide>=0 ){
	    hrescrpt_wide_genm_phi [iCent][ipt_wide][iphi]->Fill(resp_corr,wxs*wcen*wvz);
	    hresrrpt_wide_genm_phi [iCent][ipt_wide][iphi]->Fill(resp_raw ,wxs*wcen*wvz);
	    hratio_wide_phi        [iCent][ipt_wide][iphi]->Fill(ratio    ,wxs*wcen*wvz);
	  }
	}
      }
    }//! igen loop
    iEvent++;
    //std::cout<<"Completed event #  "<<i<<std::endl; 
  }//! event loop ends
  
  delete parHI_l2;
  delete parHI_l3;
  delete _JEC_HI;

  std::cout<<std::endl;
  std::cout<<runalgo << std::endl;
  for(int ic=0; ic<ncen; ic++){
    std::cout<<"\t # of events : "<< centbin[ic] << "  " << nevt[ic] << " # of jets : " << njets[ic] << std::endl;      
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
  if(eta<etabins[0] || eta>etabins[neta])return -1;
  for(int ix=0;ix<neta;ix++){
    if(eta>=etabins[ix] && eta<etabins[ix+1]){
      return ix;
    }
  }
  return -1;
}
int GetPhiBin(float phi)
{
  if(phi<phibins[0] || phi>phibins[nphi])return -1;
  for(int ix=0;ix<nphi;ix++){
    if(phi>=phibins[ix] && phi<phibins[ix+1]){
      return ix;
    }
  }
  return -1;
}

int GetPtBin(float pt)
{
  if(pt<ptbins[0] || pt>ptbins[nbins])return -1;
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

int GetPtHatBin(float pt)
{
  if(pt<pthatbins[0] || pt>pthatbins[kbins])return -1;
  for(int ix=0;ix<kbins;ix++){    
    if(pt>=pthatbins[ix] && pt<pthatbins[ix+1]){
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
double GetXsecWt(float pthat)
{
  double wt=-1.0;
  for( int i=0; i<npthat; i++){
    if( pthat >  pthatwt[i][0] )wt = pthatwt[i][3];
  }
  return wt;
}
int GetCentBin(int bin)
{
  int ibin=-1;
  //! centrality is defined as 0.5% bins of cross section
  //const char *centbin[ncen] = {"0-5%","5-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-100%"};
  if(bin>=0 && bin<10)ibin=0; //! 0-5% 
  else if(bin>=10  && bin<20 )ibin=1; //! 5-10%
  else if(bin>=20  && bin<40 )ibin=2; //! 10-20%
  else if(bin>=40  && bin<60 )ibin=3; //! 20-30%
  else if(bin>=60  && bin<80 )ibin=4; //! 30-40%
  else if(bin>=80  && bin<100)ibin=5; //! 40-50%
  else if(bin>=100 && bin<120)ibin=6; //! 50-60%
  else if(bin>=120 && bin<140)ibin=7; //! 60-70%
  else if(bin>=140 && bin<160)ibin=8; //! 70-80%
  else if(bin>=160 && bin<200)ibin=9; //! 80-100%

  return ibin;
}
