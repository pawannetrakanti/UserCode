#include <iostream>
#include <fstream>
#include <limits>
#include <iomanip>

#include "TCanvas.h"
#include "TError.h"
#include "TPad.h"
#include "TString.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TMath.h"
#include "TF1.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TChain.h"

#include "TCut.h"
#include "TNtuple.h"

#include "THStack.h"
#include <TGraph.h>

using namespace std;

int xSection(string recotype="HI"){

  TH1::SetDefaultSumw2();

  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiForest2015#5_02_TeV_pp_MC

  double pthats[] = {15       , 30       , 50       , 80       , 120      , 170      , 220      , 280      , 370      , 460      , 540      , 9999};
  double xs502[]  = {5.269e-01, 3.455e-02, 4.068e-03, 4.959e-04, 7.096e-05, 1.223e-05, 3.031e-06, 7.746e-07, 1.410e-07, 3.216e-08, 1.001e-08, 0.0};

  //double xs276[]  = {2.034e-01, 1.075e-02, 1.025e-03, 9.865e-05, 1.129e-05, 1.465e-06, 2.837e-07, 5.323e-08, 5.934e-09, 8.125e-10, 1.467e-10, 0.0};


  //double pthats[] = {15       , 30       , 50       , 80       , 120      , 170      , 220      , 280      , 9999};
  //double xs502[]  = {5.269e-01, 3.455e-02, 4.068e-03, 4.959e-04, 7.096e-05, 1.223e-05, 3.031e-06, 7.746e-07, 0.0};

  //double pthats[] = {15       , 30       , 50       , 80       , 120      , 170      , 220      , 280      , 460      , 540      , 9999};
  //double xs502 [] = {5.269e-01, 3.455e-02, 4.068e-03, 4.959e-04, 7.096e-05, 1.223e-05, 3.031e-06, 7.746e-07, 3.216e-08, 1.001e-08, 0.0};


  const int Npt = sizeof(pthats)/sizeof(double) - 1;
  double n[Npt]={0};

  // double n  [Npt] = {279441,909299,939442,790651,482031,120627,422180,764369,1143728};
  // double wtf[Npt] = {1.07364e-07,3.66964e-09,4.04208e-10,6.48706e-11,2.06051e-11,
  double *xs = xs502;
  for(int i=0; i<Npt; i++){
    cout << pthats[i] << "\t" << xs[i] << " " << xs502[i] << endl;
  }
  cout << "Npt :  " << Npt << endl;
  //return 0;

  //! Read the files and get the effective # of events

  TChain* nt;
  //TChain* evt;
  // TFile* outf[Nfiles];
  
  nt  = new TChain("ak3PFJetAnalyzer/t");
  //nt  = new TChain("akPu3PFJetAnalyzer/t");

  std::string infile_Forest;
  infile_Forest="filelist_PP_Signal_758_HcalRespCorr_v4.txt";
  //infile_Forest="filelist_PbPb_hcalv4.txt";
  //infile_Forest="filelist_ppJEC2015_"+recotype+"Reco_757p1_HcalRespCorrs_v4_00_mc.txt";
  std::ifstream infile(infile_Forest.c_str(),std::ifstream::in);
  std::string filename_Forest;

  while(1){
    infile>>filename_Forest;
    if( !infile.good() )break;
    //cout << filename_Forest.c_str() << endl;
    nt->AddFile(filename_Forest.c_str());
  }

  for(int i = 0; i < Npt; i++){
    TCut pthatCut(Form("pthat >= %f && pthat < %f",pthats[i],pthats[i+1]));
    n[i] = 1.0*nt->GetEntries(pthatCut);
    //cout<<"no of events in pthat = "<<pthats[i]<<" = "<<n[i]<<endl;
    //xs[i] = xs502[i] - xs502[i+1];
    xs[i]-= xs[i+1];
    cout << "{" << pthats[i] << "," << std::setprecision (std::numeric_limits<double>::digits10 + 1) 
	 << n[i]  << "," << std::setprecision (std::numeric_limits<double>::digits10 + 1) 
	 << xs[i] << "," << std::setprecision (std::numeric_limits<double>::digits10 + 1)  
	 << xs[i]/n[i] << "}," 
	 << endl;
  }

  //TFile *fout = new TFile(Form("xsec_weight_PbPb_%sreco_757p1_HcalRespCorrs_v4_00_mc.root",recotype.c_str()),"RECREATE");
  //TFile *fout = new TFile("xsec_weight_PbPb.root","RECREATE");
  TFile *fout = new TFile(Form("xsec_weight_pp_%sreco_758_HcalRespCorrs_v4_mc.root",recotype.c_str()),"RECREATE");
  TH1D *hpthat_comb = new TH1D("hpthat","pt hat combined",500,0,1500);
  hpthat_comb->Sumw2();

  float pthat;
  nt->SetBranchAddress("pthat",&pthat);
  nt->SetBranchStatus("*",0,0);
  nt->SetBranchStatus("pthat",1);  

  for(int i=0; i<nt->GetEntries(); i++){
    nt->GetEntry(i);
    double pthatwt=-1.0;
    for(int j = 0; j < Npt; j++){
      if(n[j] > 0 && pthat >= pthats[j]){
	pthatwt = xs[j]/n[j];
      }
    }
    if( pthatwt > 0 )hpthat_comb->Fill(pthat,pthatwt);
  }

  fout->cd();
  fout->Write();
  fout->Close();

  return 0;
}

