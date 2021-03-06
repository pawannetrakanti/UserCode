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
#include <TCanvas.h>
#include <TMath.h>

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <cstring>
#include <sstream>
#include <vector>
#include <algorithm>

#include "RooUnfold-1.1.1/src/RooUnfold.h"
#include "RooUnfold-1.1.1/src/RooUnfoldResponse.h"
#include "RooUnfold-1.1.1/src/RooUnfoldBayes.h"
#include "RooUnfold-1.1.1/src/RooUnfoldSvd.h"
#include "RooUnfold-1.1.1/src/RooUnfoldBinByBin.h"

#include "utilities.h"

using namespace std;

const int ncen=6;
const char *cdir [ncen] = {"05","510","1030","3050","5070","7090"};
const char *ccent[ncen] = {"0-5%","5-10%","10-30%","30-50%","50-70%","70-90%"};

const double ptbins[] = {
  3, 4, 5, 7, 9, 12,
  15, 18, 21, 24, 28,
  32, 37, 43, 49, 56,
  64, 74, 84, 97, 114,
  133, 153, 174, 196,
  220, 245, 300,
  330, 362, 395,
  430, 468, 507,
  548
};
const int nbins = sizeof(ptbins)/sizeof(double) - 1;


int plotUnfold_PbPb(int kRadii=3, int kJetID=0)
{
  std::string kSpecies = "pbpb";
  bool doRebin=true;
  int kRebin=20;

  gSystem->Load("RooUnfold-1.1.1/libRooUnfold.so");

  TFile *fin = new TFile(Form("Outputhist_unfold_akPu%dPF.root",kRadii),"r");

  double lumi_scale=1.;
  //double lumi_scale=1./(145.156 * 1e6);
  // if( kSpecies == "pbpb" ){
  //   lumi_scale = 1./(145.156 * 1e6);
  // } else {
  //   lumi_scale = 1./(5.3 * 1e9);
  // }


  TH1F *hgen[ncen], *hrec[ncen], *hgen_check[ncen];
  TH1D *hunf[ncen] , *hunf_check[ncen];
  TH1D *hratio[ncen], *hratio_check[ncen];

  const int kIter=4;

  for(int ic=0; ic<ncen; ic++){
    hgen[ic]  = (TH1F*)fin->Get(Form("akPu%dJetAnalyzer/hgen_f_pbpb_akPu%d_%d_%d",kRadii,kRadii,kJetID,ic));
    hgen[ic]->SetMarkerColor(kBlack);
    hgen[ic]->SetLineColor(kBlack);
    hgen[ic]->SetMarkerStyle(20);
    hgen[ic]->SetMarkerSize(0.9);
    if(doRebin)hgen[ic]->Rebin(kRebin);

    hrec[ic]  = (TH1F*)fin->Get(Form("akPu%dJetAnalyzer/hrec_f_pbpb_akPu%d_%d_%d",kRadii,kRadii,kJetID,ic));
    hrec[ic]->SetMarkerColor(kRed);
    hrec[ic]->SetLineColor(kRed);
    hrec[ic]->SetMarkerStyle(20);
    hrec[ic]->SetMarkerSize(0.9);
    if(doRebin)hrec[ic]->Rebin(kRebin);

    hunf[ic]  = (TH1D*)fin->Get(Form("akPu%dJetAnalyzer/hunf_f_pbpb_akPu%d_%d_%d"  ,kRadii,kRadii,kJetID,ic));
    hunf[ic]->SetMarkerColor(kGreen+3);
    hunf[ic]->SetLineColor(kGreen+3);
    hunf[ic]->SetMarkerStyle(20);
    hunf[ic]->SetMarkerSize(0.9);
    if(doRebin)hunf[ic]->Rebin(kRebin);

    hratio[ic] = (TH1D*)hunf[ic]->Clone(Form("hratio_f_pbpb_akPu%d_%d_%d", kRadii,kJetID,ic));
    hratio[ic]->Divide(hgen[ic]);
    hratio[ic]->SetMarkerStyle(20);
    hratio[ic]->SetMarkerColor(kBlack);
    hratio[ic]->SetMarkerSize(1.2);
    hratio[ic]->SetLineColor(kBlack);

    // //! Checking closure
    hgen_check[ic]  = (TH1F*)fin->Get(Form("akPu%dJetAnalyzer/hgen_c_f_pbpb_akPu%d_%d_%d",kRadii,kRadii,kJetID,ic));
    hgen_check[ic]->SetMarkerColor(kBlack);
    hgen_check[ic]->SetLineColor(kBlack);
    hgen_check[ic]->SetMarkerStyle(20);
    hgen_check[ic]->SetMarkerSize(0.9);
    if(doRebin)hgen_check[ic]->Rebin(kRebin);

    hunf_check[ic]  = (TH1D*)fin->Get(Form("akPu%dJetAnalyzer/hunf_c_f_pbpb_akPu%d_%d_%d"  ,kRadii,kRadii,kJetID,ic));
    hunf_check[ic]->SetMarkerColor(kMagenta+3);
    hunf_check[ic]->SetLineColor(kMagenta+3);
    hunf_check[ic]->SetMarkerStyle(20);
    hunf_check[ic]->SetMarkerSize(0.9);
    if(doRebin)hunf_check[ic]->Rebin(kRebin);
    
    
    hratio_check[ic] = (TH1D*)hunf_check[ic]->Clone(Form("hratio_c_f_pbpb_akPu%d_%d_%d", kRadii,kJetID,ic));
    hratio_check[ic]->Divide(hgen_check[ic]);
    hratio_check[ic]->SetMarkerStyle(20);
    hratio_check[ic]->SetMarkerColor(kRed);
    hratio_check[ic]->SetMarkerSize(1.2);
    hratio_check[ic]->SetLineColor(kRed);
  }

  TLegend *leg0 = getLegend(0.05, 0.75, 0.40, 0.95);
  leg0->SetHeader("Bayesian Unfolding");
  leg0->SetTextSize(0.07);
  leg0->AddEntry(hratio[0],"Train","p");
  leg0->AddEntry(hratio_check[0],"Check","p");

  TLine *line = new TLine(20.1,1,340.1,1);
  line->SetLineStyle(2); 
  line->SetLineWidth(2);

  TH1D *hDum = GetDummyHist(20.0, 340.0, 0.5, 1.5,"hDum","Gen Jet p_{T} (GeV/c)","PbPb Unfolding Closure", false);
  hDum->GetYaxis()->SetNdivisions(610); 
  hDum->GetYaxis()->SetLabelFont(43);
  hDum->GetYaxis()->SetTitleFont(43);
  hDum->GetYaxis()->SetLabelSize(20);
  hDum->GetYaxis()->SetTitleSize(22);
  hDum->GetYaxis()->SetTitleOffset(2.6);
  hDum->GetXaxis()->SetLabelFont(43);
  hDum->GetXaxis()->SetTitleFont(43);
  hDum->GetXaxis()->SetLabelSize(20);
  hDum->GetXaxis()->SetTitleSize(22);
  hDum->GetXaxis()->SetTitleOffset(2.3);
  hDum->GetXaxis()->SetNoExponent();
  hDum->GetXaxis()->SetMoreLogLabels();

  int ipad=0;
  TCanvas *c1 = new TCanvas("c1","Spectra",1100,770);
  makeMultiPanelCanvasWithGap(c1,3,2,0.01,0.01,0.16,0.2,0.04,0.04); 

  for(int ic=ncen-1;ic>=0;ic--){
    c1->cd(++ipad);
    hDum->Draw("hist");
    hratio[ic]->Draw("psame");
    hratio_check[ic]->Draw("psame");
    drawText(ccent[ic],0.648173,0.8459761,22);
    line->Draw();
    
    if(ipad==1){
      if(kJetID)drawText("wJetID",0.22,0.78,21);
      else drawText("w/o JetID",0.22,0.78,21);
      drawText("PYTHIA+HYDJET",0.22,0.88,21);
    }else if(ipad==2)leg0->Draw();
  }
  c1->SaveAs("Unfolding_Bayesian_PbPb_MC.pdf");
  return 0;
}

