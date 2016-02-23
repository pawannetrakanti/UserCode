
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


//! pt binning
double ptbins[] ={17, 22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 400, 550, 790, 1000};
//double ptbins[]={20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 250, 300, 540, 1000};
const int nbins = sizeof(ptbins)/sizeof(double) - 1;

double pthatbins[] ={15, 30, 50, 80, 120, 170, 220, 280, 370, 460, 540, 1000};
const int kbins = sizeof(pthatbins)/sizeof(double) - 1;

//! data pt binning
//double ptbins_data[] ={100, 110, 120, 130, 140, 150, 160, 170, 180, 200, 240, 300};
//const int dbins = sizeof(ptbins_data)/sizeof(double) - 1;

//! constants
int iYear=2015;
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

const int nbinsdr = 25;
const double drmin = 0.05;
const double drmax = 5.00;

const double ptbins_wide[]={20,50,80,120,200,350,800};
const int nbins_wide = sizeof(ptbins_wide)/sizeof(double) - 1;

//double etabins[] = {-3.000, -2.400, -2.000, -1.3000, 0.000, 1.300, 2.000, 2.400, 3.000};
double etabins[]={-3.139, -2.853,
		  -2.500, -2.043, -1.740, -1.392,  
		  -1.131, -0.879, -0.609, -0.435, -0.261, -0.087, 
		  +0.000, 
		  +0.087, +0.261, +0.435, +0.609,  
		  +0.879, +1.131, +1.392,  
		  +1.740, +2.043, +2.500, +2.853, 
		  +3.139
};
const int neta = sizeof(etabins)/sizeof(double) - 1;

const double phibins[] = {-3.141,-2.700,-2.100,-1.500,-0.900,-0.300, 
			  0.300,0.900,1.500,2.100,2.700,3.141
};
const int nphi = sizeof(phibins)/sizeof(double) - 1;

const int kAlgos = 12;
const char *calgo  [kAlgos]= {"ak1PF","ak2PF","ak3PF","ak4PF","ak5PF","ak6PF","ak1Calo","ak2Calo","ak3Calo","ak4Calo","ak5Calo","ak6Calo"};
string corrFileName[kAlgos]= {"AK1PF","AK2PF","AK3PF","AK4PF","AK5PF","AK6PF","AK1Calo","AK2Calo","AK3Calo","AK4Calo","AK5Calo","AK6Calo"};

//float kDelRCut=0.3;

bool is_file(const char *fileName);
bool is_file(const char *fileName)
{
  std::ifstream infile(fileName,ios::in);
  return infile.good();
}



TStopwatch timer;
int CalJec(string inname="/net/hidsk0001/d00/scratch/pawan/combinedPtHat/pp5020/HiForest_PYTHIA_QCD_merged_TuneCUETP8M1_cfi_5020GeV_757p1_NominalHICollisions2015_20151207_ppReco.root", 
	   string outname="Histo_jec_ppReco_757p1.root",const char *ksp="pp")
{

  timer.Start();


  int knj = kAlgos;

  //std::string basedir ="/net/hisrv0001/home/pawan/Validation/CMSSW_7_1_1/src/JEC_base/pp5020/";
  // std::string basedir ="";
  // std::string dirname ="JEC/";
  // std::string jecera  ="JEC_pp_PYTHIA_TuneCUETP8M1_5020GeV_patch3";


  // std::string basedir ="";
  // std::string dirname ="JEC/";
  // std::string jecera  ="JEC_pp_PYTHIAZ2_2760GeV_5320_v28";


  // std::string basedir ="/net/hisrv0001/home/pawan/Validation/CMSSW_7_1_1/src/JEC_base/pp2014/";
  // std::string dirname ="official/";
  // std::string jecera  ="JEC_pp_PYTHIAZ2_2760GeV_5320_v28";

  //! 757p1
  std::string basedir ="";
  std::string dirname ="JEC_ppReco_757p1/";
  std::string jecera  ="JEC_pp_PYTHIA_TuneCUETP8M1_5020GeV_patch3";


  TFile *fin = new TFile(inname.c_str(),"r");

  TFile *fout = new TFile(outname.c_str(),"RECREATE");
  
  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"**************************************************** "<<std::endl;
  std::cout<<Form("Analyzeforest_jec : %s  %s",ksp,dirname.c_str())<<std::endl;
  std::cout<<Form("eta cut for %0.3f ",ketacut)<<std::endl;
  std::cout<<"**************************************************** "<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;




  //! Before running check all the txt files are available
  for(int nj=0; nj<knj; nj++){
    std::string L2Name = basedir+dirname+jecera+"_L2Relative_"+corrFileName[nj]+".txt";
    if(!is_file(L2Name.c_str())){
      cout<<"**** ++++++  L2Name does not exists  "<<L2Name<<endl;
      return  2;
    }
    std::string L3Name = basedir+dirname+jecera+"_L3Absolute_"+corrFileName[nj]+".txt";
    if(!is_file(L3Name.c_str())){
      cout<<"**** ++++++  L3Name does not exists  "<<L3Name<<endl;
      return 3;
    }
  }

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

  //! Gen matched jets
  TH1F *hgenpt_genm [knj], *hrecopt_genm[knj], *hrawpt_genm[knj];
  TH1F *hgeneta     [knj], *hgenphi     [knj];
  TH1F *hjeteta     [knj], *hjetphi     [knj];


  //! Resposnse
  TH1F *hrescrpt_genm [knj][nbins], *hresrrpt_genm [knj][nbins];
  TH1F *hrescrpt_genm_eta [knj][nbins][neta], *hresrrpt_genm_eta [knj][nbins][neta];
  TH1F *hrescrpt_genm_phi [knj][nbins][nphi], *hresrrpt_genm_phi [knj][nbins][nphi];
  TH1F *hrescrpt_wide_genm_eta [knj][nbins_wide][neta], *hresrrpt_wide_genm_eta [knj][nbins_wide][neta];
  TH1F *hrescrpt_wide_genm_phi [knj][nbins_wide][nphi], *hresrrpt_wide_genm_phi [knj][nbins_wide][nphi];

  for(int nj=0;nj<knj;nj++){
    
    fout->mkdir(Form("%sJetAnalyzer",calgo[nj]));
    fout->cd(Form("%sJetAnalyzer",calgo[nj]));

    hgenpt_genm [nj] = new TH1F(Form("hgenpt_genm%d",nj),Form("Gen matched gen p_{T} distribution jet centb %s",calgo[nj]),500,0,1000);
    hgenpt_genm [nj]->Sumw2();
    hrecopt_genm[nj] = new TH1F(Form("hrecopt_genm%d",nj),Form("Gen matched reco p_{T} distribution jet %s",calgo[nj]),500,0,1000);
    hrecopt_genm[nj]->Sumw2();
    hrawpt_genm [nj] = new TH1F(Form("hrawpt_genm%d",nj),Form("Gen matched raw p_{T} distribution jet  %s",calgo[nj]),500,0,1000);
    hrawpt_genm [nj]->Sumw2();

    hjeteta[nj] = new TH1F(Form("hjeteta%d",nj),Form("jet eta distribution jet %s",calgo[nj]),72,-ketacut,ketacut);
    hjeteta[nj]->Sumw2();
    hjetphi[nj] = new TH1F(Form("hjetphi%d",nj),Form("jet phi distribution jet %s",calgo[nj]),72,-pi,pi);
    hjetphi[nj]->Sumw2();

    hgeneta[nj] = new TH1F(Form("hgeneta%d",nj),Form("gen jet eta distribution jet %s",calgo[nj]),72,-ketacut,ketacut);
    hgeneta[nj]->Sumw2();
    hgenphi[nj] = new TH1F(Form("hgenphi%d",nj),Form("gen jet phi distribution jet %s",calgo[nj]),72,-pi,pi);
    hgenphi[nj]->Sumw2();

    for(int ip=0;ip<nbins;ip++){
      //! Gen matched Response and resolution
      hrescrpt_genm [nj][ip]= new TH1F(Form("hrescrpt_genm%d_%d",nj,ip),Form("(Reco/Gen) jet p_{T} %s %0.0f < p_{T}^{REF} < %0.0f",calgo[nj],ptbins[ip],ptbins[ip+1]),rbins,rbinl,rbinh);
      hrescrpt_genm [nj][ip]->Sumw2();
      hresrrpt_genm [nj][ip]= new TH1F(Form("hresrrpt_genm%d_%d",nj,ip),Form("(Raw/Gen) jet p_{T} %s %0.0f < p_{T}^{REF} < %0.0f",calgo[nj],ptbins[ip],ptbins[ip+1]),rbins,rbinl,rbinh);
      hresrrpt_genm [nj][ip]->Sumw2();

      for(int ie=0;ie<neta;ie++){      
	hrescrpt_genm_eta [nj][ip][ie]= new TH1F(Form("hrescrpt_genm_eta%d_%d_%d",nj,ip,ie),Form("(Reco/Gen) jet p_{T} %s %d %d",calgo[nj],ip,ie),rbins,rbinl,rbinh);
	hrescrpt_genm_eta [nj][ip][ie]->Sumw2();
	hresrrpt_genm_eta [nj][ip][ie]= new TH1F(Form("hresrrpt_genm_eta%d_%d_%d",nj,ip,ie),Form("(Raw/Gen) jet p_{T} %s %d %d",calgo[nj],ip,ie),rbins,rbinl,rbinh);
	hresrrpt_genm_eta [nj][ip][ie]->Sumw2();
      }

      for(int ij=0;ij<nphi;ij++){      
	hrescrpt_genm_phi [nj][ip][ij]= new TH1F(Form("hrescrpt_genm_phi%d_%d_%d",nj,ip,ij),Form("(Reco/Gen) jet p_{T} %s %d %d",calgo[nj],ip,ij),rbins,rbinl,rbinh);
	hrescrpt_genm_phi [nj][ip][ij]->Sumw2();
	hresrrpt_genm_phi [nj][ip][ij]= new TH1F(Form("hresrrpt_genm_phi%d_%d_%d",nj,ip,ij),Form("(Raw/Gen) jet p_{T} %s %d %d",calgo[nj],ip,ij),rbins,rbinl,rbinh);
	hresrrpt_genm_phi [nj][ip][ij]->Sumw2();
      }
    }


    //! coarse pt bin
    for(int ip=0;ip<nbins_wide;ip++){
      for(int ie=0;ie<neta;ie++){      
	hrescrpt_wide_genm_eta [nj][ip][ie]= new TH1F(Form("hrescrpt_wide_genm_eta%d_%d_%d",nj,ip,ie),Form("(Reco/Gen) jet p_{T} %s %d %d",calgo[nj],ip,ie),rbins,rbinl,rbinh);
	hrescrpt_wide_genm_eta [nj][ip][ie]->Sumw2();
	hresrrpt_wide_genm_eta [nj][ip][ie]= new TH1F(Form("hresrrpt_wide_genm_eta%d_%d_%d",nj,ip,ie),Form("(Raw/Gen) jet p_{T} %s %d %d",calgo[nj],ip,ie),rbins,rbinl,rbinh);
	hresrrpt_wide_genm_eta [nj][ip][ie]->Sumw2();
      }
      for(int ij=0;ij<nphi;ij++){      
	hrescrpt_wide_genm_phi [nj][ip][ij]= new TH1F(Form("hrescrpt_wide_genm_phi%d_%d_%d",nj,ip,ij),Form("(Reco/Gen) jet p_{T} %s %d %d",calgo[nj],ip,ij),rbins,rbinl,rbinh);
	hrescrpt_wide_genm_phi [nj][ip][ij]->Sumw2();
	hresrrpt_wide_genm_phi [nj][ip][ij]= new TH1F(Form("hresrrpt_wide_genm_phi%d_%d_%d",nj,ip,ij),Form("(Raw/Gen) jet p_{T} %s %d %d",calgo[nj],ip,ij),rbins,rbinl,rbinh);
	hresrrpt_wide_genm_phi [nj][ip][ij]->Sumw2();
      }
    }
    fout->cd("../");
  }//! nj
  std::cout<<"Initialized the histograms " <<std::endl;
  ///////////////////////////////////////////////////////////////////////////////////////// 


  TTree *tr_in=0;
  //std::vector<double> vJets;
  Long64_t nbytes=0;

  double nevt [kAlgos]={0};
  double njets[kAlgos]={0};

  for(int nj=0; nj<knj; nj++){

    tr_in = (TTree*)fin->Get(Form("%sJetAnalyzer/t",calgo[nj]));
    Long64_t nentries = tr_in->GetEntries();
    std::cout<<Form("# of entries in TTree for %s %s : ",calgo[nj],ksp)<<nentries<<std::endl;

    
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
    // float neutralMax[1000];
    // float chargedMax[1000];
    // float photonMax[1000];
    float refpt[1000];
    float refeta[1000];
    float refphi[1000];
    float refdrjt[1000];
    //float refparton_pt[1000];
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
    // tr_in->SetBranchAddress("neutralMax",neutralMax);
    // tr_in->SetBranchAddress("chargedMax",chargedMax);
    // tr_in->SetBranchAddress("photonMax",photonMax);
    tr_in->SetBranchAddress("refpt",refpt);
    tr_in->SetBranchAddress("refphi",refphi);
    tr_in->SetBranchAddress("refeta",refeta);
    tr_in->SetBranchAddress("refdrjt",refdrjt);
    //tr_in->SetBranchAddress("refparton_pt",refparton_pt);
    tr_in->SetBranchAddress("refparton_flavor",refparton_flavor);
    tr_in->SetBranchAddress("subid",subid);

    //! Load the jet energy correction factors on fly
    string L2Name, L3Name;
    JetCorrectorParameters *parHI_l2=0, *parHI_l3=0;
    vector<JetCorrectorParameters> vpar_HI;
    FactorizedJetCorrector *_JEC_HI = new FactorizedJetCorrector(vpar_HI);//JR                                                                         
    
    //! Sample from which JEC was derived
    //! Official pp sample JEC
    // L2Name = "/net/hisrv0001/home/pawan/Validation/CMSSW_7_1_1/src/JEC_base/pp2014/"+dirname+"/JEC_2013RECO_STARTHI53_LV1_5_3_16_Track8_Jet28_L2Relative_"+corrFileName[nj]+".txt";
    // cout<<"**** ++++++  L2Name "<<L2Name<<endl;
    // L3Name = "/net/hisrv0001/home/pawan/Validation/CMSSW_7_1_1/src/JEC_base/pp2014/"+dirname+"/JEC_2013RECO_STARTHI53_LV1_5_3_16_Track8_Jet28_L3Absolute_"+corrFileName[nj]+".txt";

    //if(nj<knj/2)dirname="AKPF";
    //else  dirname="AKCalo";

    L2Name = basedir+dirname+jecera+"_L2Relative_"+corrFileName[nj]+".txt";
    cout<<"**** ++++++  L2Name "<<L2Name<<endl;
    L3Name = basedir+dirname+jecera+"_L3Absolute_"+corrFileName[nj]+".txt";
    cout<<"**** ++++++  L3Name "<<L3Name<<endl;
    std::cout<<std::endl;
    
    
    parHI_l2 = new JetCorrectorParameters(L2Name.c_str());
    parHI_l3 = new JetCorrectorParameters(L3Name.c_str());
    
    vpar_HI.push_back(*parHI_l2);
    vpar_HI.push_back(*parHI_l3);
    _JEC_HI = new FactorizedJetCorrector(vpar_HI);
    
    
    Int_t iEvent=0;     
    for (Long64_t i=0; i<nentries;i++) {
      //for (Long64_t i=0; i<100/*nentries*/;i++) {
      nbytes += tr_in->GetEntry(i);

      //if(fabs(vz)>15.)continue;

      double wcen=1;
      double wvz=1;

      if( pthat < 17. )continue;

      //! weight  for the merging of the samples for different pT hat bins
      double  wxs = weight;

      
      //if(i%50000==0)std::cout<<" ********** Event # " <<i<<"\t weight  : "<<weight<<"\t pthat : "<<pthat<<std::endl;
      //if(nj==3)hBin->Fill(hiBin,wxs*wcen*wvz);

      //std::cout<<" ********** Event # " <<i<<"\t weight  : "<<weight<<"\t pthat : "<<pthat<<std::endl;

      //! xsec-weight
      if(nj==3)hpthat->Fill(pthat,wxs*wcen*wvz);
      nevt[nj]++;

      //std::vector<Jet> vJet;
      //! Gen matched jets loop
      for(int igen=0; igen<nref; igen++){
	int gj = igen;

	if(fabs(refeta[gj]) > ketacut || refpt[gj] < 0 || subid[gj]!=0 )continue;
	
	if( rawpt[gj] < krawptcut ) continue;

	if( fabs(refdrjt[gj]) > kdRCut)continue;

	//if( abs(refparton_flavor[gj]) > 21 )continue;
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

	if( (abs(corrfac - recopt/rawpt[gj]) > 0.1) || corrfac > 5 ){
	  std::cout <<" WARNING " << calgo[nj] << " refpt : " << refpt[gj] <<  " rawpt :  " << rawpt[gj] 
		    << "\t recopt "<< recopt 
		    << " corr : " << corrfac << " l3corr : " << l3corr  << "\t l2corr : " << l2corr << std::endl; 
	  //continue;
	}
	
	njets[nj]++;

	// Jet selJet;
	// selJet.frefpt = refpt[gj];
	// selJet.frawpt = rawpt[gj];
	// selJet.fcorpt = recopt;
	// selJet.fparpt = refparton_pt[gj];
	// vJet.push_back(selJet);	


	float recoeta = jteta[gj];
	float recophi = jtphi[gj];
	
	double resp_corr =  recopt    / refpt[gj];
	double resp_raw  =  rawpt[gj] / refpt[gj];

	int iphi = GetPhiBin(refphi[gj]);
	int ieta = GetEtaBin(refeta[gj]);
	int ipt  = GetPtBin (refpt[gj]);
	int ipt_wide = GetPtBinWide (refpt[gj]);

	hgenpt_genm[nj]->Fill(refpt[gj],wxs*wcen*wvz);
	hrecopt_genm[nj]->Fill(recopt,wxs*wcen*wvz);
	hrawpt_genm[nj]->Fill(rawpt[gj],wxs*wcen*wvz);

	hjeteta[nj]->Fill(recoeta,wxs*wcen*wvz);
	hjetphi[nj]->Fill(recophi,wxs*wcen*wvz);

	hgeneta[nj]->Fill(refeta[gj],wxs*wcen*wvz);
	hgenphi[nj]->Fill(refphi[gj],wxs*wcen*wvz);


	//! Response in pt and eta
	if( ipt>=0 ){
	  hrescrpt_genm [nj][ipt]->Fill(resp_corr,wxs*wcen*wvz);
	  hresrrpt_genm [nj][ipt]->Fill(resp_raw ,wxs*wcen*wvz);
	  if( ieta>=0 ){
	    hrescrpt_genm_eta [nj][ipt][ieta]->Fill(resp_corr,wxs*wcen*wvz);
	    hresrrpt_genm_eta [nj][ipt][ieta]->Fill(resp_raw ,wxs*wcen*wvz);
	    if( ipt_wide>=0 ){
	      hrescrpt_wide_genm_eta [nj][ipt_wide][ieta]->Fill(resp_corr,wxs*wcen*wvz);
	      hresrrpt_wide_genm_eta [nj][ipt_wide][ieta]->Fill(resp_raw ,wxs*wcen*wvz);
	    }
	  }

	  if( iphi>=0 ){
	    hrescrpt_genm_phi [nj][ipt][iphi]->Fill(resp_corr,wxs*wcen*wvz);
	    hresrrpt_genm_phi [nj][ipt][iphi]->Fill(resp_raw ,wxs*wcen*wvz);
	    if( ipt_wide>=0 ){
	      hrescrpt_wide_genm_phi [nj][ipt_wide][iphi]->Fill(resp_corr,wxs*wcen*wvz);
	      hresrrpt_wide_genm_phi [nj][ipt_wide][iphi]->Fill(resp_raw ,wxs*wcen*wvz);
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

    //if(nj==0)break;

  }//! jet loop

  std::cout<<std::endl;
  for(int nj=0;nj<knj;nj++){
    std::cout<<calgo[nj] << std::endl;
    std::cout<<"\t # of events : "<< " pp " << "  " << nevt[nj] << "  " << njets[nj] << std::endl;      
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
