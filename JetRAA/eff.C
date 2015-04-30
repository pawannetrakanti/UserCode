#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cstring>

#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TMarker.h>
#include <TString.h>
#include <TVirtualFitter.h>
#include "../Utils/utilities.h"

//void CorrectBinWidth(TH1F *&/*h1*/);
//void setBinomialErrors(TGraphAsymmErrors* grEff,const TH1F* hEnum, const TH1F* hDenom);

const double etabins[] = {0, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0};
const int neta = sizeof(etabins)/sizeof(double) - 1;
const char *seta[neta] = {"|#eta|<1.0","1.0<|#eta|<1.1","1.1<|#eta|<1.2","1.2<|#eta|<1.3","1.3<|#eta|<1.4",
			  "1.4<|#eta|<1.5","1.5<|#eta|<1.6","1.6<|#eta|<1.7","1.7<|#eta|<1.8","1.8<|#eta|<1.9"
			  ,"1.9<|#eta|<2.0"};

int eff(int wJetID=1)
{
  LoadStyle();
  gStyle->SetErrorX(0);

  
  cout << " neta : " << neta << endl;
  //return 0;

  int iSave=1;

  int plotCent=0;
  int plotRadii=3;

  std::string  plotAlgo="PF";
  int id=wJetID;  
  
  const int ncen=6;
  const char *ccent[ncen] = {"0-5%","5-10%","10-30%","30-50%","50-70%","70-90%"};
  
  const int knj=3;
  std::string calgo_pbpb[knj]={"akPu2","akPu3","akPu4"};

  //! Input file PbPb
  TFile *fin_pbpb = NULL;

  fin_pbpb = new TFile("OutputHist_ntuples_short_PbPb.root","r");

  TGraphAsymmErrors *gr_eff_pt [knj][ncen], *gr_eff_eta [knj][ncen], *gr_eff_phi [knj][ncen];
  TGraphAsymmErrors *gr_eff_pt_eta[knj][ncen][neta];

  Color_t icol[ncen] ={kViolet+2, kGreen-3, kRed-3, kBlue, kMagenta, kOrange-2};
  int isty[ncen]     ={24       , 25      , 26    , 30   , 32      , 27       };
  TH1F *hden_eff_pt =0, *hnum_eff_pt =0;
  TH1F *hden_eff_eta=0, *hnum_eff_eta=0;
  TH1F *hden_eff_phi=0, *hnum_eff_phi=0;

  TFile *fin=0;
  std::string algname="";
  for(int nj=0;nj<knj;nj++){

    for(int ic=0;ic<ncen; ic++){

      fin     = fin_pbpb;     
      algname = calgo_pbpb[nj];

      hden_eff_pt = (TH1F*)fin->Get(Form("hPtAll_%s_%d_%d" ,(algname+plotAlgo).c_str(),id,ic));
      hnum_eff_pt = (TH1F*)fin->Get(Form("hPtEff_%s_%d_%d" ,(algname+plotAlgo).c_str(),id,ic));
      gr_eff_pt[nj][ic] = new TGraphAsymmErrors(hnum_eff_pt, hden_eff_pt,"cl=0.683 b(1,1) mode");
      gr_eff_pt[nj][ic]->SetName(Form("gr_eff_pt_%d_%d",nj,ic));
      gr_eff_pt[nj][ic]->SetMarkerStyle(isty[ic]);
      gr_eff_pt[nj][ic]->SetMarkerColor(icol[ic]);
      gr_eff_pt[nj][ic]->SetLineColor(icol[ic]);
      gr_eff_pt[nj][ic]->SetMarkerSize(1.6);

      hden_eff_eta = (TH1F*)fin->Get(Form("hEtaAll_%s_%d_%d" ,(algname+plotAlgo).c_str(),id,ic));
      hnum_eff_eta = (TH1F*)fin->Get(Form("hEtaEff_%s_%d_%d" ,(algname+plotAlgo).c_str(),id,ic));
      gr_eff_eta[nj][ic] = new TGraphAsymmErrors(hnum_eff_eta, hden_eff_eta,"cl=0.683 b(1,1) mode");
      gr_eff_eta[nj][ic]->SetName(Form("gr_eff_eta_%d_%d",nj,ic));
      gr_eff_eta[nj][ic]->SetMarkerStyle(isty[ic]);
      gr_eff_eta[nj][ic]->SetMarkerColor(icol[ic]);
      gr_eff_eta[nj][ic]->SetLineColor(icol[ic]);
      gr_eff_eta[nj][ic]->SetMarkerSize(1.6);

      hden_eff_phi = (TH1F*)fin->Get(Form("hPhiAll_%s_%d_%d" ,(algname+plotAlgo).c_str(),id,ic));
      hnum_eff_phi = (TH1F*)fin->Get(Form("hPhiEff_%s_%d_%d" ,(algname+plotAlgo).c_str(),id,ic));
      gr_eff_phi[nj][ic] = new TGraphAsymmErrors(hnum_eff_phi, hden_eff_phi,"cl=0.683 b(1,1) mode");
      gr_eff_phi[nj][ic]->SetName(Form("gr_eff_eta_%d_%d",nj,ic));
      gr_eff_phi[nj][ic]->SetMarkerStyle(isty[ic]);
      gr_eff_phi[nj][ic]->SetMarkerColor(icol[ic]);
      gr_eff_phi[nj][ic]->SetLineColor(icol[ic]);
      gr_eff_phi[nj][ic]->SetMarkerSize(1.6);



      for(int ie=0; ie<neta; ie++){
	hden_eff_pt = (TH1F*)fin->Get(Form("hPtAll_etabin_%s_%d_%d_%d" ,(algname+plotAlgo).c_str(),id,ic,ie));
	hnum_eff_pt = (TH1F*)fin->Get(Form("hPtEff_etabin_%s_%d_%d_%d" ,(algname+plotAlgo).c_str(),id,ic,ie));
	gr_eff_pt_eta[nj][ic][ie] = new TGraphAsymmErrors(hnum_eff_pt, hden_eff_pt,"cl=0.683 b(1,1) mode");
	gr_eff_pt_eta[nj][ic][ie]->SetName(Form("gr_eff_pt_eta_%s_%d_%d",(algname+plotAlgo).c_str(),ic,ie));
	gr_eff_pt_eta[nj][ic][ie]->SetMarkerStyle(isty[ic]);
	gr_eff_pt_eta[nj][ic][ie]->SetMarkerColor(icol[ic]);
	gr_eff_pt_eta[nj][ic][ie]->SetLineColor(icol[ic]);
	gr_eff_pt_eta[nj][ic][ie]->SetMarkerSize(1.6);
      }

      // //setBinomialErrors(gr_eff_pt[nj][ic], hnum_eff, hden_eff);
      // //! Binwidth correction
      // // CorrectBinWidth(hden_eff);
      // // CorrectBinWidth(hnum_eff);
      // // CorrectBinWidth(hden_fak);
      // // CorrectBinWidth(hnum_fak);

    }
  }//! nj

  //gr_eff_pt_eta[1][0][0]->Draw("ap");
  //return 0;


  TLegend *leg = new TLegend(0.4312431,0.2531486,0.6325633,0.6209068,"","BRNDC");
  if(wJetID)leg->SetHeader("w JetID");
  else leg->SetHeader("w/o JetID");
  leg->SetBorderSize(0);
  leg->SetLineColor(10);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04);
  for(int ic=0; ic<ncen; ic++){
    leg->AddEntry(gr_eff_pt[knj-1][ic],ccent[ic],"p");
  }

  TF1 *fitfunc = new TF1("fitfunc","([0]/2.)*(1 + TMath::Erf((x - [1])/[2]))",20,620);
  fitfunc->SetParameters(0.92,0.4,1.0);

  //976,863
  int ipad=0;
  int maxc=4;
  int maxr=3;
  TCanvas *c1[knj], *c11[knj], *c111[knj], *c2[knj], *c22[knj], *c222[knj];
  TCanvas *ceta[knj];
  for(int nj=0; nj<knj; nj++){
    c1[nj] = new TCanvas(Form("c1_%d",nj),Form("%s Eff ",(calgo_pbpb[nj]+plotAlgo).c_str()),975,865);
    c1[nj]->cd();
    c1[nj]->SetLeftMargin(0.13);
    c1[nj]->SetBottomMargin(0.18);
    //gPad->SetLogx();
    for(int ic=0; ic<ncen; ic++){
      gr_eff_pt[nj][ic]->SetMaximum(1.015);
      //gr_eff_pt[nj][ic]->SetMinimum(0.);
      gr_eff_pt[nj][ic]->SetMinimum(0.74);
      gr_eff_pt[nj][ic]->SetTitle("");
      gr_eff_pt[nj][ic]->GetXaxis()->SetRangeUser(21,330);
      gr_eff_pt[nj][ic]->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
      gr_eff_pt[nj][ic]->GetXaxis()->SetLabelSize(0.05);
      gr_eff_pt[nj][ic]->GetXaxis()->SetTitleSize(0.05);
      gr_eff_pt[nj][ic]->GetXaxis()->CenterTitle();
      gr_eff_pt[nj][ic]->GetXaxis()->SetDecimals(true);    
      gr_eff_pt[nj][ic]->GetXaxis()->SetNoExponent();
      gr_eff_pt[nj][ic]->GetXaxis()->SetMoreLogLabels();
      
      gr_eff_pt[nj][ic]->GetYaxis()->SetTitle("Jet reconstruction efficiency");
      gr_eff_pt[nj][ic]->GetYaxis()->CenterTitle();
      gr_eff_pt[nj][ic]->GetYaxis()->SetTitleOffset(1.26);
      gr_eff_pt[nj][ic]->GetYaxis()->SetLabelSize(0.05);
      gr_eff_pt[nj][ic]->GetYaxis()->SetTitleSize(0.05);

      string drawOpt = (ic==0) ? "apz" : "psamez";
      gr_eff_pt[nj][ic]->Draw(drawOpt.c_str());
      if(ic==0){
	TLine *line = new TLine(21,1.0,330,1.0);
	line->SetLineStyle(2);
	line->Draw();
      }
    }
    leg->Draw();
    drawText2(Form("%s, PbPb (PYTHIA+HYDJET), |#eta|<2",(calgo_pbpb[nj]+plotAlgo).c_str()),0.22,0.93,29);
    if(iSave){
      c1[nj]->SaveAs(Form("ANPlots/Eff/EffPt_%d_PbPb_%s.gif",wJetID,(calgo_pbpb[nj]+plotAlgo).c_str()));
      c1[nj]->SaveAs(Form("ANPlots/Eff/EffPt_%d_PbPb_%s.pdf",wJetID,(calgo_pbpb[nj]+plotAlgo).c_str()));
    }


    //! in different eta bins
    ceta[nj] = new TCanvas(Form("ceta_%d",nj),Form("%s",(calgo_pbpb[nj]+plotAlgo).c_str()),1237,839);
    ceta[nj]->Divide(4,3);
    //makeMultiPanelCanvas(ceta[nj],maxc,maxr,0.0,0.0,0.22,0.22,0.02,0);
    ipad=0;
    for(int ic=0; ic<ncen; ic++){
      ipad=0;
      for(int ie=0;ie<neta;ie++){
	ceta[nj]->cd(++ipad);
	gr_eff_pt_eta[nj][ic][ie]->SetMaximum(1.015);
	gr_eff_pt_eta[nj][ic][ie]->SetMinimum(0.74);
	gr_eff_pt_eta[nj][ic][ie]->SetTitle("");
	gr_eff_pt_eta[nj][ic][ie]->GetXaxis()->SetRangeUser(21,330);
	gr_eff_pt_eta[nj][ic][ie]->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
	gr_eff_pt_eta[nj][ic][ie]->GetXaxis()->SetLabelSize(0.05);
	gr_eff_pt_eta[nj][ic][ie]->GetXaxis()->SetTitleSize(0.05);
	gr_eff_pt_eta[nj][ic][ie]->GetXaxis()->CenterTitle();
	gr_eff_pt_eta[nj][ic][ie]->GetXaxis()->SetDecimals(true);    
	gr_eff_pt_eta[nj][ic][ie]->GetXaxis()->SetNoExponent();
	gr_eff_pt_eta[nj][ic][ie]->GetXaxis()->SetMoreLogLabels();
	
	gr_eff_pt_eta[nj][ic][ie]->GetYaxis()->SetTitle("Jet reconstruction efficiency");
	gr_eff_pt_eta[nj][ic][ie]->GetYaxis()->CenterTitle();
	gr_eff_pt_eta[nj][ic][ie]->GetYaxis()->SetTitleOffset(1.26);
	gr_eff_pt_eta[nj][ic][ie]->GetYaxis()->SetLabelSize(0.05);
	gr_eff_pt_eta[nj][ic][ie]->GetYaxis()->SetTitleSize(0.05);

	string drawOpt = (ic==0) ? "apz" : "psamez";
	gr_eff_pt_eta[nj][ic][ie]->Draw(drawOpt.c_str());
	if(ic==0 & ipad==1){
	  drawText2(Form("%s",(calgo_pbpb[nj]+plotAlgo).c_str()),0.35,0.73,21);
	  drawText2("PYTHIA+HYDJET",0.35,0.63,21);
	}
	if(ic==0)drawText2(seta[ie],0.35,0.53,24);
      }
      if(iSave){
	ceta[nj]->SaveAs(Form("ANPlots/Eff/EffPt_etabins_%d_PbPb_%s.gif",wJetID,(calgo_pbpb[nj]+plotAlgo).c_str()));
	ceta[nj]->SaveAs(Form("ANPlots/Eff/EffPt_etabins_%d_PbPb_%s.pdf",wJetID,(calgo_pbpb[nj]+plotAlgo).c_str()));
      }
    }

    // /// eta
    // c11[nj] = new TCanvas(Form("c11_%d",nj),Form("%s Eta Eff ",(calgo_pbpb[nj]+plotAlgo).c_str()),975,865);
    // c11[nj]->cd();
    // c11[nj]->SetLeftMargin(0.13);
    // c11[nj]->SetBottomMargin(0.18);
    // for(int ic=0; ic<ncen; ic++){
    //   gr_eff_eta[nj][ic]->SetMaximum(1.015);
    //   gr_eff_eta[nj][ic]->SetMinimum(0.0);
    //   //gr_eff_eta[nj][ic]->SetMinimum(0.60);
    //   gr_eff_eta[nj][ic]->SetTitle("");
    //   gr_eff_eta[nj][ic]->GetXaxis()->SetRangeUser(-2.0,2.0);
    //   gr_eff_eta[nj][ic]->GetXaxis()->SetTitle("GenJet #eta");
    //   gr_eff_eta[nj][ic]->GetXaxis()->SetLabelSize(0.05);
    //   gr_eff_eta[nj][ic]->GetXaxis()->SetTitleSize(0.05);
    //   gr_eff_eta[nj][ic]->GetXaxis()->CenterTitle();
    //   gr_eff_eta[nj][ic]->GetXaxis()->SetDecimals(true);    
    //   gr_eff_eta[nj][ic]->GetXaxis()->SetNoExponent();
    //   gr_eff_eta[nj][ic]->GetXaxis()->SetMoreLogLabels();
    //   gr_eff_eta[nj][ic]->GetXaxis()->SetNdivisions(505);
    //   gr_eff_eta[nj][ic]->GetYaxis()->SetTitle("Jet reconstruction efficiency");
    //   gr_eff_eta[nj][ic]->GetYaxis()->CenterTitle();
    //   gr_eff_eta[nj][ic]->GetYaxis()->SetTitleOffset(1.26);
    //   gr_eff_eta[nj][ic]->GetYaxis()->SetLabelSize(0.05);
    //   gr_eff_eta[nj][ic]->GetYaxis()->SetTitleSize(0.05);
    //   gr_eff_eta[nj][ic]->GetYaxis()->SetNdivisions(505);

    //   string drawOpt = (ic==0) ? "apz" : "psamez";
    //   gr_eff_eta[nj][ic]->Draw(drawOpt.c_str());
    //   if(ic==0){
    // 	TLine *line = new TLine(-2.0,1.0,2.0,1.0);
    // 	line->SetLineStyle(2);
    // 	line->Draw();
    //   }
    // }
    // leg->Draw();
    // drawText2(Form("%s, PbPb (PYTHIA+HYDJET)",(calgo_pbpb[nj]+plotAlgo).c_str()),0.22,0.93,29);
    // //c11[nj]->Close();
    // if(iSave){
    //   c11[nj]->SaveAs(Form("../ANPlots//EffEta_%d_PbPb_%s.gif",wJetID,(calgo_pbpb[nj]+plotAlgo).c_str()));
    //   c11[nj]->SaveAs(Form("../ANPlots//EffEta_%d_PbPb_%s.pdf",wJetID,(calgo_pbpb[nj]+plotAlgo).c_str()));
    // }


    // /// phi
    // c111[nj] = new TCanvas(Form("c111_%d",nj),Form("%s phi Eff ",(calgo_pbpb[nj]+plotAlgo).c_str()),975,865);
    // c111[nj]->cd();
    // c111[nj]->SetLeftMargin(0.13);
    // c111[nj]->SetBottomMargin(0.18);
    // for(int ic=0; ic<ncen; ic++){
    //   gr_eff_phi[nj][ic]->SetMaximum(1.015);
    //   gr_eff_phi[nj][ic]->SetMinimum(0.0);
    //   //gr_eff_phi[nj][ic]->SetMinimum(0.60);
    //   gr_eff_phi[nj][ic]->SetTitle("");
    //   gr_eff_phi[nj][ic]->GetXaxis()->SetRangeUser(-3.14,3.14);
    //   gr_eff_phi[nj][ic]->GetXaxis()->SetTitle("GenJet #phi");
    //   gr_eff_phi[nj][ic]->GetXaxis()->SetLabelSize(0.05);
    //   gr_eff_phi[nj][ic]->GetXaxis()->SetTitleSize(0.05);
    //   gr_eff_phi[nj][ic]->GetXaxis()->CenterTitle();
    //   gr_eff_phi[nj][ic]->GetXaxis()->SetDecimals(true);    
    //   gr_eff_phi[nj][ic]->GetXaxis()->SetNoExponent();
    //   gr_eff_phi[nj][ic]->GetXaxis()->SetMoreLogLabels();

    //   gr_eff_phi[nj][ic]->GetYaxis()->SetTitle("Jet reconstruction efficiency");
    //   gr_eff_phi[nj][ic]->GetYaxis()->CenterTitle();
    //   gr_eff_phi[nj][ic]->GetYaxis()->SetTitleOffset(1.26);
    //   gr_eff_phi[nj][ic]->GetYaxis()->SetLabelSize(0.05);
    //   gr_eff_phi[nj][ic]->GetYaxis()->SetTitleSize(0.05);
    //   gr_eff_phi[nj][ic]->GetYaxis()->SetNdivisions(507);

    //   string drawOpt = (ic==0) ? "apz" : "psamez";
    //   gr_eff_phi[nj][ic]->Draw(drawOpt.c_str());
    //   if(ic==0){
    // 	TLine *line = new TLine(-3.14,1.0,3.14,1.0);
    // 	line->SetLineStyle(2);
    // 	line->Draw();
    //   }
    // }
    // leg->Draw();
    // drawText2(Form("%s, PbPb (PYTHIA+HYDJET), |#eta|<2.0",(calgo_pbpb[nj]+plotAlgo).c_str()),0.22,0.93,29);
    // //c111[nj]->Close();
    // if(iSave){
    //   c111[nj]->SaveAs(Form("../ANPlots//EffPhi_%d_PbPb_%s.gif",wJetID,(calgo_pbpb[nj]+plotAlgo).c_str()));
    //   c111[nj]->SaveAs(Form("../ANPlots//EffPhi_%d_PbPb_%s.pdf",wJetID,(calgo_pbpb[nj]+plotAlgo).c_str()));
    // }
  }
  return 0; 
}
// void CorrectBinWidth(TH1F *&h1)
// {
//   for(int ix=1; ix<=h1->GetNbinsX(); ix++){
//     h1->SetBinContent(ix, h1->GetBinContent(ix)/h1->GetBinWidth(ix));
//   }
// }
// void setBinomialErrors(TGraphAsymmErrors* grEff,const TH1F* hnum, const TH1F* hden)
// {
//   for (int i=0;i<grEff->GetN();i++) {
//     float nenum =hnum ->GetBinContent(i+1);
//     float ndenom=hden->GetBinContent(i+1);
//     float eeff=(ndenom>0.0) ? std::sqrt(nenum/(ndenom*ndenom)*(1-nenum/ndenom)):0.0;
//     grEff->SetPointError(i,0,0,eeff,eeff);
//   }
// }
