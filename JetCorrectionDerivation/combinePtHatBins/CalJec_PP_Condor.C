
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


int GetPhiBin(float /*phi*/);
int GetEtaBin(float /*eta*/);
int GetPtHatBin(float /*pt*/);
int GetPtBin(float /*pt*/);
int GetPtBinWide(float /*pt*/);
double delphi(double /*phi1*/, double /*phi2*/);
//

//! pt binning
//double ptbins[] ={17, 22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 400, 550, 790, 1200};
double ptbins[] ={15, 17, 22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 400, 550, 790, 1000};
const int nbins = sizeof(ptbins)/sizeof(double) - 1;

double pthatbins[] ={15, 30, 50, 80, 120, 170, 220, 280, 370, 460, 540, 1000};
const int kbins = sizeof(pthatbins)/sizeof(double) - 1;

//! data pt binning
//double ptbins_data[] ={100, 110, 120, 130, 140, 150, 160, 170, 180, 200, 240, 300};
//const int dbins = sizeof(ptbins_data)/sizeof(double) - 1;

//! constants
int iYear=2014;
const double pi=acos(-1.);
//const double pi2=2*pi -1;

const float  ketacut   =3.2;
const double krawptcut =1.0;
const double kptrecocut=0.0;
const double kptgencut =0.0;
const double kdRCut=0.30;
int rbins=150;
double rbinl=0.0;
double rbinh=2.0;

//const double etabins[]={-3.000,
//			      -2.500, -2.043, -1.930, -1.740, -1.653, -1.566, -1.392,  
//			      -1.218, -1.131, -0.957, -0.879, -0.783, -0.609, -0.522, 
//			      -0.435, -0.348, -0.261, -0.087, 	    	    	    
//			      +0.000, 	    	    	    	    	    	    
//			      +0.087, +0.261, +0.348, +0.435, +0.522, +0.609, +0.783, 
//			      +0.879, +0.957, +1.131, +1.218, +1.392, +1.566, 
//			      +1.653, +1.740, +1.930, +2.043, +2.500,  
//			      +3.000
//};

const int nbinsdr = 25;
const double drmin = 0.05;
const double drmax = 5.00;

const double ptbins_wide[]={20,40,60,80,110,200,350,800};
const int nbins_wide = sizeof(ptbins_wide)/sizeof(double) - 1;

//const double etabins[]={-3.0,-2.4,-1.8,-1.4,-1.0,-0.8,-0.6,-0.4,-0.2,
//			0.0,0.2,0.4,0.6,0.8,1.0,1.4,1.8,2.4,3.0};

// double etabins[] = {-3.000,
// 		    -2.500, -1.740, -1.392,  
// 		    -1.131, -0.879, -0.435, -0.087, 
// 		    +0.000, 
// 		    +0.087, +0.435, +0.879, +1.131, +1.392,  
// 		    +1.740, +2.500,
// 		    +3.000
// };
//const double etabins[]={-3.0,-2.0,-1.4,-0.8,0.0,0.8,1.4,2.0,3.0};
//double etabins[] = {-3.000, -2.000, -1.392, -0.435, -0.087, 0.000, 0.087, 0.435, 1.392, 2.000, 3.000};


 
// double etabins[]={-3.000,
// 		  -2.500, -2.043, -1.740, -1.392,  
// 		  -1.131, -0.879, -0.609, -0.435, -0.261, -0.087, 
// 		  +0.000, 
// 		  +0.087, +0.261, +0.435, +0.609,
// 		  +0.879, +1.131, +1.392,  
// 		  +1.740, +2.043, +2.500,  
// 		  +3.000
// };
//double etabins[] = {-3.000, -2.400, -2.000, -1.3000, 0.000, 1.300, 2.000, 2.400, 3.000};


// double etabins[] ={-3.000,
// 		   -2.500, -2.043, -1.740, -1.653, -1.566, -1.392,  
// 		   -1.218, -1.131, -0.957, -0.879, -0.783, -0.609, 
// 		   -0.522, -0.435, -0.348, -0.261, -0.087, 
// 		   +0.000, 
// 		   +0.087, +0.261, +0.348, +0.435, +0.522, +0.609, 
// 		   +0.783, +0.879, +0.957, +1.131, +1.218, +1.392, 
// 		   +1.566, +1.653, +1.740, +1.930, +2.043, +2.500,
// 		   +3.000
// };


// double etabins[] ={-5.191, -4.363, -3.664, 
// 		   -3.139, -2.853, -2.500, -2.322, -2.172, 
// 		   -2.043, -1.930, -1.830, -1.740, -1.653, -1.566, -1.479, -1.392, -1.305, 
// 		   -1.218, -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, 
// 		   -0.435, -0.348, -0.261, -0.174, -0.087, 
// 		   +0.000,                            
// 		   +0.087, +0.174, +0.261, +0.348, +0.435, +0.522, +0.609, +0.696, +0.783, 
// 		   +0.879, +0.957, +1.044, +1.131, +1.218, +1.305, +1.392, +1.479, +1.566, 
// 		   +1.653, +1.740, +1.830, +1.930, +2.043, +2.172, +2.322, +2.500,  
// 		   +2.853, +3.139, +3.664,   
// 		   +4.363, +5.191
// };

// double etabins[] ={-4.363, -3.664, 
// 		   -3.139, -2.853, -2.500, -2.322, -2.172, 
// 		   -2.043, -1.930, -1.830, -1.740, -1.653, -1.566, -1.479, -1.392, -1.305, 
// 		   -1.218, -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, 
// 		   -0.435, -0.348, -0.261, -0.174, -0.087, 
// 		   +0.000,                            
// 		   +0.087, +0.174, +0.261, +0.348, +0.435, +0.522, +0.609, +0.696, +0.783, 
// 		   +0.879, +0.957, +1.044, +1.131, +1.218, +1.305, +1.392, +1.479, +1.566, 
// 		   +1.653, +1.740, +1.830, +1.930, +2.043, +2.172, +2.322, +2.500,  
// 		   +2.853, +3.139, +3.664, +4.363
// };
// double etabins[]={-3.664, -3.139, -2.853,
// 		  -2.500, -2.043, -1.740, -1.392,  
// 		  -1.131, -0.879, -0.609, -0.435, -0.261, -0.087, 
// 		  +0.000, 
// 		  +0.087, +0.261, +0.435, +0.609,  
// 		  +0.879, +1.131, +1.392,  
// 		  +1.740, +2.043, +2.500, +2.853, 
// 		  +3.139, +3.664
// };


double etabins[] = {-4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, 
		    -3.489, -3.314, -3.139, -2.964, -2.853, -2.650, -2.500, -2.322, -2.172, 
		    -2.043, -1.930, -1.830, -1.740, -1.653, -1.566, -1.479, -1.392, -1.305, 
		    -1.218, -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, 
		    -0.435, -0.348, -0.261, -0.174, -0.087, 
		    +0.000, 
		    +0.087, +0.174, +0.261, +0.348, +0.435, +0.522, +0.609, +0.696, +0.783, 
		    +0.879, +0.957, +1.044, +1.131, +1.218, +1.305, +1.392, +1.479, +1.566, 
		    +1.653, +1.740, +1.830, +1.930, +2.043, +2.172, +2.322, +2.500, +2.650, 
		    +2.853, +2.964, +3.139, +3.314, +3.489, +3.664, +3.839, +4.013, +4.191, 
		    +4.363, +4.538, +4.716
};
const int neta = sizeof(etabins)/sizeof(double) - 1;

const double phibins[] = {-3.141,-2.700,-2.100,-1.500,-0.900,-0.300, 
			  0.300,0.900,1.500,2.100,2.700,3.141
};
const int nphi = sizeof(phibins)/sizeof(double) - 1;

//const int kAlgos = 1;
//string salgo       [kAlgos]= {"ak1PF","ak2PF","ak3PF","ak4PF","ak5PF","ak6PF","ak1Calo","ak2Calo","ak3Calo","ak4Calo","ak5Calo","ak6Calo"};
//string corrFileName[kAlgos]= {"AK1PF","AK2PF","AK3PF","AK4PF","AK5PF","AK6PF","AK1Calo","AK2Calo","AK3Calo","AK4Calo","AK5Calo","AK6Calo"};
//string corrFileName[kAlgos]={"AK3PF"};

bool is_file(const char *fileName);
bool is_file(const char *fileName)
{
  std::ifstream infile(fileName,ios::in);
  return infile.good();
}



TStopwatch timer;
int CalJec_PP_Condor(const char *runalgo="ak1Calo", const char *outfile="test.root", const char *ksp="pp")
{

  timer.Start();

  int nj=0;
  string corrFileName = "";
  if( strcmp(runalgo,"ak1PF") == 0){corrFileName="AK1PF";nj=0;}
  else if( strcmp(runalgo,"ak2PF") == 0){corrFileName="AK2PF";nj=1;}
  else if( strcmp(runalgo,"ak3PF") == 0){corrFileName="AK3PF";nj=2;}
  else if( strcmp(runalgo,"ak4PF") == 0){corrFileName="AK4PF";nj=3;}
  else if( strcmp(runalgo,"ak5PF") == 0){corrFileName="AK5PF";nj=4;}
  else if( strcmp(runalgo,"ak6PF") == 0){corrFileName="AK6PF";nj=5;}
  else if( strcmp(runalgo,"ak1Calo") == 0){corrFileName="AK1Calo";nj=6;}
  else if( strcmp(runalgo,"ak2Calo") == 0){corrFileName="AK2Calo";nj=7;}
  else if( strcmp(runalgo,"ak3Calo") == 0){corrFileName="AK3Calo";nj=8;}
  else if( strcmp(runalgo,"ak4Calo") == 0){corrFileName="AK4Calo";nj=9;}
  else if( strcmp(runalgo,"ak5Calo") == 0){corrFileName="AK5Calo";nj=10;}
  else if( strcmp(runalgo,"ak6Calo") == 0){corrFileName="AK6Calo";nj=11;}


  //std::string basedir ="/net/hisrv0001/home/pawan/Validation/CMSSW_7_1_1/src/JEC_base/pp5020/";
  //std::string basedir ="/net/hisrv0001/home/pawan/Validation/CMSSW_7_1_1/src/combinePtHatBins/ppJEC2015/";

  // std::string basedir ="";
  // std::string dirname ="JEC/";
  // std::string jecera  ="JEC_pp_PYTHIA_TuneCUETP8M1_5020GeV_patch3";

  // std::string basedir ="";
  // std::string dirname ="JEC_HIReco/";
  // std::string jecera  ="JEC_pp_PYTHIA_TuneCUETP8M1_5020GeV_patch3";

  // std::string basedir ="";
  // //std::string dirname ="JEC_ppReco_757p1/";
  // std::string dirname ="test_jec/";
  // std::string jecera  ="JEC_pp_PYTHIA_TuneCUETP8M1_5020GeV_patch3";

  std::string basedir ="";
  std::string dirname ="JEC_HIReco_HcalRespCorrs_v4_00_mc/";
  std::string jecera  ="JEC_pp_PYTHIA_TuneCUETP8M1_5020GeV_patch3";

  // std::string basedir ="/net/hisrv0001/home/pawan/Validation/CMSSW_7_1_1/src/JEC_base/pp2014/";
  // //std::string basedir ="";
  // std::string dirname ="official/";
  // std::string jecera  ="JEC_pp_PYTHIAZ2_2760GeV_5320_v28";


  //std::string dirname = "official";
  //std::string dirname = "AKCalo";

  std::string inname="";
  inname="/mnt/hadoop/cms/store/user/pawan/pp5020/HiForest_PYTHIA_QCD_merged_TuneCUETP8M1_cfi_5020GeV_20151229_ppReco_757p1_HcalRespCorrs_v4_00_mc.root";

  //inname="/mnt/hadoop/cms/store/user/pawan/pp5020/HiForest_PYTHIA_QCD_merged_TuneCUETP8M1_cfi_5020GeV_757p1_NominalHICollisions2015_20151207_ppReco.root";

  //inname="/mnt/hadoop/cms/store/user/pawan/pp5020/HiForest_PYTHIA_QCD_merged_TuneCUETP8M1_cfi_5020GeV_757p1_NominalHICollisions2015_20151207_ppReco.root";

  //inname="/mnt/hadoop/cms/store/user/pawan/pp5020/HiForest_PYTHIA_QCD_merged_TuneCUETP8M1_cfi_5020GeV_patch3_NominalHICollisions2015_20151130_ppReco.root";

  //inname="/mnt/hadoop/cms/store/user/pawan/pp5020/HIReco_ntuple/HiForest_PYTHIA_QCD_merged_TuneCUETP8M1_cfi_5020GeV_patch3_NominalHICollisions2015_20151126_HIReco.root"; //! is with pp reco

  //inname="/mnt/hadoop/cms/store/user/pawan/pp5020/HiForest_PYTHIA_QCD_merged_TuneCUETP8M1_cfi_5020GeV_patch3_NominalHICollisions2015_20151112_ppReco_nocut.root";
  //inname="/mnt/hadoop/cms/store/user/pawan/pp5020/HiForest_PYTHIA_QCD_merged_TuneCUETP8M1_cfi_5020GeV_patch3_NominalHICollisions2015_20151010_ppReco.root";
  //inname="/mnt/hadoop/cms/store/user/pawan/pp5020/HiForest_PYTHIA_QCD_merged_TuneCUETP8M1_cfi_5020GeV_patch3_NominalHICollisions2015_20151010_HIReco.root"; //! is with pp reco
  //outname="JetResponse_histos_ppSignal_ppReco_PYTHIA_TuneCUETP8M1_5020GeV_patch3.root";  

  TFile *fin = new TFile(inname.c_str(),"r");

  TFile *fout = new TFile(outfile,"RECREATE");
  
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
  TH1F *hgenpt_genm , *hrecopt_genm, *hrawpt_genm;
  TH1F *hgeneta     , *hgenphi     ;
  TH1F *hjeteta     , *hjetphi     ;


  //! Resposnse
  TH1F *hrescrpt_genm [nbins], *hresrrpt_genm [nbins];
  TH1F *hrescrpt_genm_eta [nbins][neta], *hresrrpt_genm_eta [nbins][neta];
  TH1F *hrescrpt_genm_phi [nbins][nphi], *hresrrpt_genm_phi [nbins][nphi];
  TH1F *hrescrpt_wide_genm_eta [nbins_wide][neta], *hresrrpt_wide_genm_eta [nbins_wide][neta];
  TH1F *hrescrpt_wide_genm_phi [nbins_wide][nphi], *hresrrpt_wide_genm_phi [nbins_wide][nphi];
    
  fout->mkdir(Form("%sJetAnalyzer",runalgo));
  fout->cd(Form("%sJetAnalyzer",runalgo));

  TH1F *hBin    = new TH1F(Form("hBin_%d",nj),"Centrality bin",200,-0.5,200-0.5);
  hBin->Sumw2();
  TH1F *hpthat  = new TH1F(Form("hpthat_%d",nj),"pt-hat distribution",500,0,1000);
  hpthat->Sumw2();

  
  hgenpt_genm  = new TH1F(Form("hgenpt_genm%d",nj),Form("Gen matched gen p_{T} distribution jet centb %s",runalgo),500,0,1000);
  hgenpt_genm ->Sumw2();
  hrecopt_genm = new TH1F(Form("hrecopt_genm%d",nj),Form("Gen matched reco p_{T} distribution jet %s",runalgo),500,0,1000);
  hrecopt_genm->Sumw2();
  hrawpt_genm  = new TH1F(Form("hrawpt_genm%d",nj),Form("Gen matched raw p_{T} distribution jet  %s",runalgo),500,0,1000);
  hrawpt_genm ->Sumw2();
  
  hjeteta = new TH1F(Form("hjeteta%d",nj),Form("jet eta distribution jet %s",runalgo),72,-ketacut,ketacut);
  hjeteta->Sumw2();
  hjetphi = new TH1F(Form("hjetphi%d",nj),Form("jet phi distribution jet %s",runalgo),72,-pi,pi);
  hjetphi->Sumw2();
  
  hgeneta = new TH1F(Form("hgeneta%d",nj),Form("gen jet eta distribution jet %s",runalgo),72,-ketacut,ketacut);
  hgeneta->Sumw2();
  hgenphi = new TH1F(Form("hgenphi%d",nj),Form("gen jet phi distribution jet %s",runalgo),72,-pi,pi);
  hgenphi->Sumw2();
  
  for(int ip=0;ip<nbins;ip++){
    //! Gen matched Response and resolution
    hrescrpt_genm [ip]= new TH1F(Form("hrescrpt_genm%d_%d",nj,ip),Form("(Reco/Gen) jet p_{T} %s %0.0f < p_{T}^{REF} < %0.0f",runalgo,ptbins[ip],ptbins[ip+1]),rbins,rbinl,rbinh);
    hrescrpt_genm [ip]->Sumw2();
    hresrrpt_genm [ip]= new TH1F(Form("hresrrpt_genm%d_%d",nj,ip),Form("(Raw/Gen) jet p_{T} %s %0.0f < p_{T}^{REF} < %0.0f",runalgo,ptbins[ip],ptbins[ip+1]),rbins,rbinl,rbinh);
    hresrrpt_genm [ip]->Sumw2();
    
    for(int ie=0;ie<neta;ie++){      
      hrescrpt_genm_eta [ip][ie]= new TH1F(Form("hrescrpt_genm_eta%d_%d_%d",nj,ip,ie),Form("(Reco/Gen) jet p_{T} %s %d %d",runalgo,ip,ie),rbins,rbinl,rbinh);
      hrescrpt_genm_eta [ip][ie]->Sumw2();
      hresrrpt_genm_eta [ip][ie]= new TH1F(Form("hresrrpt_genm_eta%d_%d_%d",nj,ip,ie),Form("(Raw/Gen) jet p_{T} %s %d %d",runalgo,ip,ie),rbins,rbinl,rbinh);
      hresrrpt_genm_eta [ip][ie]->Sumw2();
    }
    
    for(int ij=0;ij<nphi;ij++){      
      hrescrpt_genm_phi [ip][ij]= new TH1F(Form("hrescrpt_genm_phi%d_%d_%d",nj,ip,ij),Form("(Reco/Gen) jet p_{T} %s %d %d",runalgo,ip,ij),rbins,rbinl,rbinh);
      hrescrpt_genm_phi [ip][ij]->Sumw2();
      hresrrpt_genm_phi [ip][ij]= new TH1F(Form("hresrrpt_genm_phi%d_%d_%d",nj,ip,ij),Form("(Raw/Gen) jet p_{T} %s %d %d",runalgo,ip,ij),rbins,rbinl,rbinh);
      hresrrpt_genm_phi [ip][ij]->Sumw2();
    }
  }
  
  
  //! coarse pt bin
  for(int ip=0;ip<nbins_wide;ip++){
    for(int ie=0;ie<neta;ie++){      
      hrescrpt_wide_genm_eta [ip][ie]= new TH1F(Form("hrescrpt_wide_genm_eta%d_%d_%d",nj,ip,ie),Form("(Reco/Gen) jet p_{T} %s %d %d",runalgo,ip,ie),rbins,rbinl,rbinh);
      hrescrpt_wide_genm_eta [ip][ie]->Sumw2();
      hresrrpt_wide_genm_eta [ip][ie]= new TH1F(Form("hresrrpt_wide_genm_eta%d_%d_%d",nj,ip,ie),Form("(Raw/Gen) jet p_{T} %s %d %d",runalgo,ip,ie),rbins,rbinl,rbinh);
      hresrrpt_wide_genm_eta [ip][ie]->Sumw2();
    }
    for(int ij=0;ij<nphi;ij++){      
      hrescrpt_wide_genm_phi [ip][ij]= new TH1F(Form("hrescrpt_wide_genm_phi%d_%d_%d",nj,ip,ij),Form("(Reco/Gen) jet p_{T} %s %d %d",runalgo,ip,ij),rbins,rbinl,rbinh);
      hrescrpt_wide_genm_phi [ip][ij]->Sumw2();
      hresrrpt_wide_genm_phi [ip][ij]= new TH1F(Form("hresrrpt_wide_genm_phi%d_%d_%d",nj,ip,ij),Form("(Raw/Gen) jet p_{T} %s %d %d",runalgo,ip,ij),rbins,rbinl,rbinh);
      hresrrpt_wide_genm_phi [ip][ij]->Sumw2();
    }
  }
  fout->cd("../");
  std::cout<<"Initialized the histograms " <<std::endl;
  ///////////////////////////////////////////////////////////////////////////////////////// 


  TTree *tr_in=0;
  //std::vector<double> vJets;
  Long64_t nbytes=0;

  double nevt =0;
  double njets=0;
  
  tr_in = (TTree*)fin->Get(Form("%sJetAnalyzer/t",runalgo));
  Long64_t nentries = tr_in->GetEntries();
  std::cout<<Form("# of entries in TTree for %s %s : ",runalgo,ksp)<<nentries<<std::endl;
  
  
  //Declaration of leaves types
  // float vz;
  // int hiBin;
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
  float refparton_pt[1000];
  int refparton_flavor[1000];
  int subid[1000];

  // tr_in->SetBranchAddress("vz",&vz);
  // tr_in->SetBranchAddress("hiBin",&hiBin);
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
  tr_in->SetBranchAddress("refparton_flavor",refparton_flavor);
  tr_in->SetBranchAddress("subid",subid);

  //! Load the jet energy correction factors on fly
  JetCorrectorParameters *parHI_l2=0, *parHI_l3=0;
  vector<JetCorrectorParameters> vpar_HI;
    
  parHI_l2 = new JetCorrectorParameters(L2Name.c_str());
  parHI_l3 = new JetCorrectorParameters(L3Name.c_str());

  //cout <<   L3Name.c_str() <<  endl;
    
  vpar_HI.push_back(*parHI_l2);
  vpar_HI.push_back(*parHI_l3);
  //_JEC_HI = new FactorizedJetCorrector(vpar_HI);
  FactorizedJetCorrector *_JEC_HI = new FactorizedJetCorrector(vpar_HI);//JR                                                                             
    
  Int_t iEvent=0;     
  for (Long64_t i=0; i<nentries;i++) {
    nbytes += tr_in->GetEntry(i);
    
    float wcen=1;
    float wvz=1;
    //! weight  for the merging of the samples for different pT hat bins
    float wxs = weight;

    //if(i%50000==0)std::cout<<" ********** Event # " <<i<<"\t weight  : "<<weight<<"\t pthat : "<<pthat<<std::endl;
    //if(nj==3)hBin->Fill(hiBin,wxs*wcen*wvz);

    //std::cout<<" ********** Event # " <<i<<"\t weight  : "<<weight<<"\t pthat : "<<pthat<<std::endl;

    if( pthat < 17. )continue;


    //! xsec-weight
    hpthat->Fill(pthat,wxs*wcen*wvz);
    nevt++;

    //! Gen matched jets loop
    for(int igen=0; igen<nref; igen++){
      int gj = igen;

      if( fabs(refeta[gj]) > ketacut || refpt[gj] < 0 || subid[gj]!=0 )continue;
	
      if( rawpt[gj] < krawptcut ) continue;

      if( fabs(refdrjt[gj]) > kdRCut)continue;

      //if( subid[gj]!=0 || fabs(refeta[gj]) > ketacut || rawpt[gj] < krawptcut  
      //  ||  refpt[gj] < 0 ||  refpt[gj] < kptgencut || rawpt[gj]/refpt[gj] > 1)continue;
	
      _JEC_HI->setJetEta(jteta[gj]);
      _JEC_HI->setJetPt (rawpt[gj]);
	
      float corrfac = _JEC_HI->getCorrection();

      std::vector <float> subcorr = _JEC_HI->getSubCorrections();
      float l2corr = subcorr[0];
      float l3corr = subcorr[1];
      float recopt  = rawpt[gj]*corrfac;  
      // float recopt_l2corr = rawpt[gj]*l2corr;
      // float recopt_l3corr = rawpt[gj]*l3corr;

      //std::cout << " \t" << runalgo << " " << pthat << " refpt : " << refpt[gj] <<  " rawpt :  " << rawpt[gj] << "\t recopt "<< recopt << " l3 corr "<< l3corr << " l2 corr :  "<< recopt_l2corr << " corr : " << corrfac << std::endl; 
	
      if( (abs(corrfac - recopt/rawpt[gj]) > 0.1) || corrfac > 5 ){
	std::cout <<" WARNING " << runalgo << " refpt : " << refpt[gj] <<  " rawpt :  " << rawpt[gj] << "\t recopt "<< recopt 
		  << " corr : " << corrfac << "\t l3corr : "<< l3corr << "\t l2corr : " << l2corr << std::endl; 
	//continue;
      }
      njets++;


      float recoeta = jteta[gj];
      float recophi = jtphi[gj];
      //float delr    = refdrjt[gj];
	
      double resp_corr =  recopt    / refpt[gj];
      double resp_raw  =  rawpt[gj] / refpt[gj];

      int iphi = GetPhiBin(refphi[gj]);
      int ieta = GetEtaBin(refeta[gj]);
      int ipt  = GetPtBin (refpt[gj]);
      int ipt_wide = GetPtBinWide (refpt[gj]);

      hgenpt_genm ->Fill(refpt[gj],wxs*wcen*wvz);
      hrecopt_genm->Fill(recopt,wxs*wcen*wvz);
      hrawpt_genm ->Fill(rawpt[gj],wxs*wcen*wvz);

      hjeteta->Fill(recoeta,wxs*wcen*wvz);
      hjetphi->Fill(recophi,wxs*wcen*wvz);

      hgeneta->Fill(refeta[gj],wxs*wcen*wvz);
      hgenphi->Fill(refphi[gj],wxs*wcen*wvz);

      //! Response in pt and eta
      if( ipt>=0 ){
	hrescrpt_genm [ipt]->Fill(resp_corr,wxs*wcen*wvz);
	hresrrpt_genm [ipt]->Fill(resp_raw ,wxs*wcen*wvz);
	if( ieta>=0 ){
	  hrescrpt_genm_eta [ipt][ieta]->Fill(resp_corr,wxs*wcen*wvz);
	  hresrrpt_genm_eta [ipt][ieta]->Fill(resp_raw ,wxs*wcen*wvz);
	  if( ipt_wide>=0 ){
	    hrescrpt_wide_genm_eta [ipt_wide][ieta]->Fill(resp_corr,wxs*wcen*wvz);
	    hresrrpt_wide_genm_eta [ipt_wide][ieta]->Fill(resp_raw ,wxs*wcen*wvz);
	  }
	}

	if( iphi>=0 ){
	  hrescrpt_genm_phi [ipt][iphi]->Fill(resp_corr,wxs*wcen*wvz);
	  hresrrpt_genm_phi [ipt][iphi]->Fill(resp_raw ,wxs*wcen*wvz);
	  if( ipt_wide>=0 ){
	    hrescrpt_wide_genm_phi [ipt_wide][iphi]->Fill(resp_corr,wxs*wcen*wvz);
	    hresrrpt_wide_genm_phi [ipt_wide][iphi]->Fill(resp_raw ,wxs*wcen*wvz);
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
  std::cout<<"\t # of events : "<< " pp " << "  " << nevt << "  " << njets << std::endl;      
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
