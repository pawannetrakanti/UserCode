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

#include "../Utils/utilities.h"

const int ncen=8;
const char *cdir [ncen] = {"05","510","1030","3050","5070","7090","90100","pp"};
const char *ccent[ncen] = {"0-5%","5-10%","10-30%","30-50%","50-70%","70-90%","90-100%","pp"};
const int dsty   [ncen] = {20, 24, 21, 25, 33, 27, 29, 30};
const int dcol   [ncen] = {kRed, kBlue, kBlack, kAzure+8, kGreen+3, kViolet+1, kOrange+1, kGray+1};
double ptbins[] ={40, 50 ,60 ,70 ,80 ,90 ,100, 110, 120, 130, 140, 160, 200, 250, 300, 400, 548};
const int nbins = sizeof(ptbins)/sizeof(double) - 1;

const int knj=3;
std::string srad[knj] = {"2", "3", "4"};

const char *run[2] = {"DATA","MC"};

void MakeIndHist(TH1F *&/*h1*/,float /*ymax*/,float /*ymin*/);

int DrawBkg(int drawBin=0)
{

  int iSave=1;
  TFile *fin = new TFile("OutputHist_ntuples_28042015.root","r");

  TH1F *havbkgpt_wjid[2][knj][ncen][nbins];

  //! 0: data , 1: MC
  TH1F *hbmean[2][knj][ncen], *hbrms[2][knj][ncen];

  std::string algname="";
  TH1F *hist=0;
  for(int in=0; in<2; in++){
    for(int nj=0; nj<knj; nj++){
      for(int ic=0; ic<ncen; ic++){    
	if(ic==ncen-1)algname="ak"+srad[nj]+"PF";
	else algname="akPu"+srad[nj]+"PF";
      
	//! bkg mean from data
	hbmean[in][nj][ic] = new TH1F(Form("hbmean_%s_%s_%d",run[in],algname.c_str(),ic),Form("hbmean_%s_%s_%s",run[in],algname.c_str(),ccent[ic]),nbins,ptbins);
	hbmean[in][nj][ic]->Sumw2();
	MakeHist(hbmean[in][nj][ic],"Reco Jet p_{T} (GeV/c)","< Background p_{T} > (GeV/c)");
	hbmean[in][nj][ic]->SetMarkerStyle(dsty[ic]);
	hbmean[in][nj][ic]->SetMarkerColor(dcol[ic]);
	hbmean[in][nj][ic]->SetLineColor(dcol[ic]);
	hbmean[in][nj][ic]->SetMarkerSize(1.2);

	if(nj==0){
	  hbmean[in][nj][ic]->SetMinimum(-0.8);
	  hbmean[in][nj][ic]->SetMaximum(84);
	}else if(nj==1){
	  hbmean[in][nj][ic]->SetMinimum(-0.8);
	  hbmean[in][nj][ic]->SetMaximum(114);
	}else if(nj==2){
	  hbmean[in][nj][ic]->SetMinimum(-0.8);
	  hbmean[in][nj][ic]->SetMaximum(184);
	}
	hbmean[in][nj][ic]->GetXaxis()->SetTitleSize(0.06);
	hbmean[in][nj][ic]->GetYaxis()->SetTitleSize(0.06);
	hbmean[in][nj][ic]->GetYaxis()->SetTitleOffset(1.13);

	hbrms[in][nj][ic] = new TH1F(Form("hbrms_data_%s_%s_%d",run[in],algname.c_str(),ic),Form("hbrms_%s_%s_%s",run[in],algname.c_str(),ccent[ic]),nbins,ptbins);
	hbrms[in][nj][ic]->Sumw2();
	MakeHist(hbrms[in][nj][ic],"Reco Jet p_{T} (GeV/c)","RMS ( Background p_{T} ) (GeV/c)");
	hbrms[in][nj][ic]->SetMarkerStyle(dsty[ic]);
	hbrms[in][nj][ic]->SetMarkerColor(dcol[ic]);
	hbrms[in][nj][ic]->SetLineColor(dcol[ic]);
	hbrms[in][nj][ic]->SetMarkerSize(1.2);
	if(nj==0){
	  hbrms[in][nj][ic]->SetMinimum(-0.2);
	  hbrms[in][nj][ic]->SetMaximum(12);
	}else if(nj==1){
	  hbrms[in][nj][ic]->SetMinimum(-0.2);
	  hbrms[in][nj][ic]->SetMaximum(22);
	}else{
	  hbrms[in][nj][ic]->SetMinimum(-0.2);
	  hbrms[in][nj][ic]->SetMaximum(34);
	}
	hbrms[in][nj][ic]->GetXaxis()->SetTitleSize(0.06);
	hbrms[in][nj][ic]->GetYaxis()->SetTitleSize(0.06);
	hbrms[in][nj][ic]->GetYaxis()->SetTitleOffset(0.86);
	
	for(int ip=0; ip<nbins; ip++){    
	  hist = (TH1F*)fin->Get(Form("havbkgpt_wjid_comb_%s_%d_%d_%d",algname.c_str(),in,ic,ip));
	  hbmean[in][nj][ic]->SetBinContent(ip+1, hist->GetMean());
	  hbmean[in][nj][ic]->SetBinError  (ip+1, hist->GetMeanError());
	  hbrms[in][nj][ic]->SetBinContent (ip+1, hist->GetRMS());
	  hbrms[in][nj][ic]->SetBinError   (ip+1, hist->GetRMSError());


	  havbkgpt_wjid[in][nj][ic][ip] = (TH1F*)fin->Get(Form("havbkgpt_wjid_comb_%s_%d_%d_%d",algname.c_str(),in,ic,ip));
	  if(havbkgpt_wjid[in][nj][ic][ip]->Integral()>0){
	    havbkgpt_wjid[in][nj][ic][ip]->Scale(1./havbkgpt_wjid[in][nj][ic][ip]->Integral());
	  }else{
	    havbkgpt_wjid[in][nj][ic][ip]->Scale(0.);
	  }
	  if(in==0){//! Data
	    havbkgpt_wjid[in][nj][ic][ip]->SetLineColor(1);
	    havbkgpt_wjid[in][nj][ic][ip]->SetMarkerColor(1);
	    havbkgpt_wjid[in][nj][ic][ip]->SetMarkerStyle(20);
	  }else{
	    havbkgpt_wjid[in][nj][ic][ip]->SetLineColor(2);
	    havbkgpt_wjid[in][nj][ic][ip]->SetMarkerColor(2);
	    havbkgpt_wjid[in][nj][ic][ip]->SetMarkerStyle(20);
	    havbkgpt_wjid[in][nj][ic][ip]->SetMarkerSize(1.2);
	    havbkgpt_wjid[in][nj][ic][ip]->SetFillStyle(3005);
	    havbkgpt_wjid[in][nj][ic][ip]->SetFillColor(2);
	    havbkgpt_wjid[in][nj][ic][ip]->SetLineWidth(2);
	  }
	}
      }
    }
  }

  gStyle->SetErrorX(0);
  TLegend   *leg0 = new TLegend(0.1733478,0.5879662,0.4349987,0.9085542,NULL,"brNDC");
  leg0->SetBorderSize(1);
  leg0->SetTextSize(0.04);
  leg0->SetLineColor(10);
  leg0->SetLineStyle(1);
  leg0->SetLineWidth(1);
  leg0->SetFillColor(10);
  leg0->SetFillStyle(1001);
  leg0->AddEntry(hbmean[0][1][0],"PbPb","");
  for(int ic=0; ic<ncen-2; ic++){
    leg0->AddEntry(hbmean[0][1][ic],ccent[ic],"p");
  }
  leg0->AddEntry(hbmean[1][1][0],"PYTHIA + HYDJET","l");


  TCanvas *cbkg[knj];
  for(int nj=0; nj<knj; nj++){
    cbkg[nj] = new TCanvas(Form("cbkg_%d",nj),Form("%s <bkgd> pT",srad[nj].c_str()),1270,550);
    cbkg[nj]->Divide(2,1);

    cbkg[nj]->cd(1);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.06);
    gPad->SetTopMargin(0.03);
    gPad->SetBottomMargin(0.20);



    for(int ic=ncen-3; ic>=0; ic--){
      hbmean[0][nj][ic]->GetXaxis()->SetRangeUser(40,300);
      hbmean[1][nj][ic]->GetXaxis()->SetRangeUser(40,300);
      
      if(ic==ncen-3){
	hbmean[1][nj][ic]->Draw("hist");
	leg0->Draw();

	drawText2("(a)  CMS",0.19,0.92,22);
	drawText2(Form("Anti-k_{T}, PF Jets, R = 0.%s",srad[nj].c_str()),0.40,0.88,21);
	drawText2("| #eta | < 2.0",0.40,0.82,21);
      }else {
	hbmean[1][nj][ic]->Draw("histsame");
      }
      hbmean[0][nj][ic]->Draw("psame0");
    }
    

    cbkg[nj]->cd(2);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.06);
    gPad->SetTopMargin(0.03);
    gPad->SetBottomMargin(0.20);

    for(int ic=ncen-3; ic>=0; ic--){
      hbrms[0][nj][ic]->GetXaxis()->SetRangeUser(40,300);
      hbrms[1][nj][ic]->GetXaxis()->SetRangeUser(40,300);
      
      if(ic==ncen-3){
	hbrms[1][nj][ic]->Draw("hist");
	drawText2("(b)  CMS",0.19,0.92,22);
	drawText2(Form("Anti-k_{T}, PF Jets, R = 0.%s",srad[nj].c_str()),0.40,0.88,21);
	drawText2("| #eta | < 2.0",0.40,0.82,21);
      }else {
	hbrms[1][nj][ic]->Draw("histsame");
      }
      hbrms[0][nj][ic]->Draw("psame0");
    }
    if(iSave){
      cbkg[nj]->SaveAs(Form("ANPlots/New/Bkgd_Mean_RMS_R%s.pdf",srad[nj].c_str()));
      cbkg[nj]->SaveAs(Form("ANPlots/New/Bkgd_Mean_RMS_R%s.png",srad[nj].c_str()));
    }
  }

  int maxc=3;
  int maxr=2;

  //int drawBin=6;
  int ipad=0;

  TLegend *leg1 = getLegend(0.6245365,0.475705,0.8578539,0.7242597);
  leg1->SetHeader("");
  leg1->SetTextSize(0.09);
  leg1->AddEntry(havbkgpt_wjid[0][1][0][drawBin],"Data","p");
  leg1->AddEntry(havbkgpt_wjid[1][1][0][drawBin],"MC"  ,"lf");

  TCanvas *cind[knj];  
  for(int nj=0; nj<knj; nj++){
    cind[nj] = new TCanvas(Form("cind_%d",nj),Form("R=%s w jetid Match Bakcground pT",srad[nj].c_str()),1097,609);
    makeMultiPanelCanvas(cind[nj],maxc,maxr,0.0,0.0,0.22,0.22,0.02,0);
    ipad=0;
    for(int ic=ncen-3; ic>=0; ic--){
      cind[nj]->cd(++ipad);
      gPad->SetLogy();

      MakeIndHist(havbkgpt_wjid[0][nj][ic][drawBin], 5.25, 2.25e-08);
      MakeIndHist(havbkgpt_wjid[1][nj][ic][drawBin], 5.25, 2.25e-08);


      if(ipad==1){
	havbkgpt_wjid[1][nj][ic][drawBin]->GetYaxis()->SetTitleSize(0.09);
	havbkgpt_wjid[1][nj][ic][drawBin]->GetYaxis()->SetTitleOffset(1.14);
	havbkgpt_wjid[1][nj][ic][drawBin]->GetYaxis()->SetLabelSize(0.09);
	havbkgpt_wjid[1][nj][ic][drawBin]->GetYaxis()->SetLabelOffset(0.005);
      }

      havbkgpt_wjid[1][nj][ic][drawBin]->Draw("hist");
      havbkgpt_wjid[0][nj][ic][drawBin]->Draw("psame");
      if(ipad==1){
	drawText2(Form("Anti-k_{T}, PF Jets, R = 0.%s",srad[nj].c_str()),0.30,0.78,21);
	drawText2(" | #eta | < 2.0",0.31,0.66,21);
	leg1->Draw();
      }
      drawText2(ccent[ic],0.30,0.90,21);
      if(ipad==2){
	std::ostringstream strs; 
	strs << ptbins[drawBin] << " < p_{T} (GeV/c) < " << ptbins[drawBin+1];
	std::string spt = strs.str();
	drawText2(spt.c_str(),0.30,0.80,21);
      }
    }
    gPad->Update();
    if(iSave){
      cind[nj]->SaveAs(Form("ANPlots/New/AvbkgdpT_wjid_R%s_%0.0f_pT_%0.0f.pdf",srad[nj].c_str(),ptbins[drawBin],ptbins[drawBin+1]));
      cind[nj]->SaveAs(Form("ANPlots/New/AvbkgdpT_wjid_R%s_%0.0f_pT_%0.0f.png",srad[nj].c_str(),ptbins[drawBin],ptbins[drawBin+1]));
    }
  }


  return 0;
}
void MakeIndHist(TH1F *&h1,float ymax,float ymin)
{
  h1->SetMaximum(ymax);
  h1->SetMinimum(ymin);
  h1->SetTitle("");
  //h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->GetXaxis()->SetTitle("Jet Background p_{T} (GeV/c)");
  h1->GetXaxis()->CenterTitle(true);
  h1->GetXaxis()->SetMoreLogLabels();
  h1->GetXaxis()->SetNoExponent();
  h1->GetXaxis()->SetTitleFont(42);
  h1->GetXaxis()->SetLabelFont(42);
  h1->GetXaxis()->SetTitleSize(0.07);
  h1->GetXaxis()->SetTitleOffset(1.15);
  h1->GetXaxis()->SetLabelSize(0.07);
  h1->GetXaxis()->SetLabelOffset(0.005);
  h1->GetXaxis()->SetNdivisions(507);
  h1->GetYaxis()->SetTitle("Event Fraction");
  h1->GetYaxis()->CenterTitle(true);
  h1->GetYaxis()->SetTitleSize(0.07);
  h1->GetYaxis()->SetTitleOffset(1.50);
  h1->GetYaxis()->SetLabelSize(0.07);
  h1->GetYaxis()->SetNdivisions(507);
  h1->GetYaxis()->SetDecimals(true);
  h1->GetYaxis()->SetTitleFont(42);
  h1->GetYaxis()->SetLabelFont(42);
}
