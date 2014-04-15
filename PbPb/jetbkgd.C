#include <string.h>
#include <iostream>
#include "utilities.h"
#include <vector>

#include <TVirtualPad.h>
#include <TPaveStats.h>
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>

using namespace std;
int GetPtBin(float /*pt*/);
void MakePad(TVirtualPad */*cpad*/,int /*logx*/,int /*logy*/);
void GetPileup(TH2 */*h1*/,TH1 */*hmean*/,TH1 */*hrms*/);

const int kfiles=2;
const char *version[kfiles] = {"data","mc"};
const char *cleg[kfiles]    = {"Data","MC"};
Color_t dcol[kfiles]={kBlack,kRed};
int     dsty[kfiles]={20,24};

const int knj=1;
const char *calgo[knj]  = {"akPu3PF"}; 
const char *algn [knj]  = {"Anti-k_{T}, PF, R = 0.3"};
const char *cnpbpb[knj] = {"akPu3PF"}; 
const int iAlg[knj]     = {6};

// const char *calgo[knj]={"akPu3PF","akVs3PF"}; 
// const char *algn[knj] = {"Anti-k_{T}, PU PF, R = 0.3","Anti-k_{T}, VS PF, R = 0.3"};
// const char *cnpbpb[knj]  =  {"akPu3PF","akPu3VS"}; 

double ptbins[] = {30,40,50,60,70,75,80,85,90,100};
//double ptbins[] = {30,50,80,100,120,140,160,180,200,250};
//double ptbins[] = {10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,240,280,300,346};
const int b1  = sizeof(ptbins)/sizeof(Double_t) - 1;
const int bins = b1;

const int ncen=6;
const char *ccent[ncen] ={"0-5%","5-10%","10-30%","30-50%","50-70%","70-90%"};
const char *mcent[ncen] ={"05","510","1030","3050","5070","7090"};
int pileup(int isel=0)
{
  LoadStyle();

  cout<<" # of pT bins : " << bins << endl;

  TFile *fin[kfiles];
  fin[0] = new TFile("input/jetbkg_data.root","r"); 
  fin[1] = new TFile("input/jetbkg_mc_genmatched.root","r"); 
  //fin[1] = new TFile("input/jetbkg_mc_incl.root","r"); 

  //! Background for jets
  TH2F *hjetptbkgd[kfiles][knj][ncen];
  TH1F *hbmean[kfiles][knj][ncen];
  TH1F *hbrms [kfiles][knj][ncen];
  TH1F *hjet1D[kfiles][knj][ncen][bins];

  const char *ctrig[] = {"Combined",
			 "Jet80",
			 "Jet65 && !Jet80",
			 "#splitline{Jet55_ps && Jet55}{&& !Jet65 && !Jet80}"
  };

  for(int in=0; in<kfiles; in++){
    for(int nj=0;nj<knj;nj++){
      for(int ic=0;ic<ncen;ic++){
	if(in==0){ //! data
	  if(isel==0){//! combined
	    hjetptbkgd[in][nj][ic] = (TH2F*)fin[in]->Get(Form("hjbkgComb_%d",ic)); 
	  }else if(isel==1){//! HLT_HIJet80_v1
	    hjetptbkgd[in][nj][ic] = (TH2F*)fin[in]->Get(Form("hjbkg%d_0",ic));
	  }else if(isel==2){//! HLT_HIJet65_v1 && !HLT_HIJet80_v1
	    hjetptbkgd[in][nj][ic] = (TH2F*)fin[in]->Get(Form("hjbkg%d_1",ic));
	  }else{ //! HLT_HIJet55_v1_Prescl && HLT_HIJet55_v1 && !HLT_HIJet65_v1 && !HLT_HIJet80_v1
	    hjetptbkgd[in][nj][ic] = (TH2F*)fin[in]->Get(Form("hjbkg%d_2",ic));
	  }
	} else { //! mc
	  hjetptbkgd[in][nj][ic] = (TH2F*)fin[in]->Get(Form("hjbkgComb%d_%d",iAlg[nj],ic));//! Inclusive no cut
	  //hjetptbkgd[in][nj][ic] = (TH2F*)fin[in]->Get(Form("hjbkgComb_match%d_%d",iAlg[nj],ic));//! genmatched
	}
	hjetptbkgd[in][nj][ic]->SetName(Form("hjetptbkgd_%s_%s_%d",version[in],calgo[nj],ic));

	//! pileup mean from data
	hbmean[in][nj][ic] = new TH1F(Form("hbmean_%s_%s_%d",version[in],calgo[nj],ic),Form("hbmean_%s_%s_%s",version[in],calgo[nj],ccent[ic]),bins,ptbins);
	hbmean[in][nj][ic]->Sumw2();
	MakeHist(hbmean[in][nj][ic],"Reco Jet p_{T} (GeV/c)","< Jet bkgd >");
	hbmean[in][nj][ic]->SetMarkerStyle(dsty[in]);
	hbmean[in][nj][ic]->SetMarkerColor(dcol[in]);
	hbmean[in][nj][ic]->SetLineColor(dcol[in]);

	hbrms[in][nj][ic] = new TH1F(Form("hbrms_data_%s_%s_%d",version[in],calgo[nj],ic),Form("hbrms_%s_%s_%s",version[in],calgo[nj],ccent[ic]),bins,ptbins);
	hbrms[in][nj][ic]->Sumw2();
	MakeHist(hbrms[in][nj][ic],"Reco Jet p_{T} (GeV/c)","#sigma( Jet bkgd )");
	hbrms[in][nj][ic]->SetMarkerStyle(dsty[in]);
	hbrms[in][nj][ic]->SetMarkerColor(dcol[in]);
	hbrms[in][nj][ic]->SetLineColor(dcol[in]);

	GetPileup(hjetptbkgd[in][nj][ic],hbmean[in][nj][ic],hbrms[in][nj][ic]);

	for(int ix=0;ix<bins;ix++){
	  
	  //int lbin = (int)(ptbins[ix]/binw + 1); 
	  //int hbin = (int)(ptbins[ix+1]/binw + 1); 

	  int lbin = hjetptbkgd[in][nj][ic]->GetXaxis()->FindBin(ptbins[ix]); 
	  int hbin = hjetptbkgd[in][nj][ic]->GetXaxis()->FindBin(ptbins[ix+1]); 

	  hjet1D[in][nj][ic][ix] = (TH1F*)hjetptbkgd[in][nj][ic]->ProjectionY(Form("hjet1D_%s_%d_%d_%d",version[in],nj,ic,ix),lbin,hbin);
	  hjet1D[in][nj][ic][ix]->Rebin(5);
	  hjet1D[in][nj][ic][ix]->Scale(1./hjet1D[in][nj][ic][ix]->Integral());
	  hjet1D[in][nj][ic][ix]->SetMarkerColor(dcol[in]);
	  hjet1D[in][nj][ic][ix]->SetLineColor(dcol[in]);
	  hjet1D[in][nj][ic][ix]->SetMarkerStyle(dsty[in]);
	  if(in==1){
	    hjet1D[in][nj][ic][ix]->SetFillColor(dcol[in]);
	    hjet1D[in][nj][ic][ix]->SetFillStyle(3004);
	  }
	  MakeHist(hjet1D[in][nj][ic][ix],"Bkgd Jet p_{T} (GeV/c)","Event Fraction");
	  hjet1D[in][nj][ic][ix]->GetXaxis()->SetRangeUser(0,100);
	}
      }
    }
  }

  // hjet1D[0][0][0][0]->Draw("p");
  // hjet1D[1][0][0][0]->Draw("psame");
  // return 0;

  TLegend *leg = new TLegend(0.7345273,0.6274441,0.9071231,0.7950419,"","BRNDC");
  leg->SetHeader("");
  leg->SetBorderSize(0);
  leg->SetLineColor(10);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(10);
  leg->SetTextSize(0.07);
  leg->AddEntry(hbmean[0][0][0],cleg[0],"p");
  leg->AddEntry(hbmean[1][0][0],cleg[1],"p");


  TLegend *leg1 = new TLegend(0.6009514,0.5476372,0.9383781,0.7985264,"","BRNDC");
  leg1->SetHeader("");
  leg1->SetBorderSize(0);
  leg1->SetLineColor(10);
  leg1->SetLineStyle(1);
  leg1->SetLineWidth(1);
  leg1->SetFillColor(10);
  leg1->SetTextSize(0.07);
  leg1->AddEntry(hjet1D[0][0][0][0],cleg[0],"p");
  leg1->AddEntry(hjet1D[1][0][0][0],cleg[1],"f");


  int ipad=0;
  // TCanvas *c0= new TCanvas("c0","Jetbkgd distributions",1030,923);
  // c0->Divide(2,2);
  // for(int in=0; in<kfiles; in++){
  //   ipad=0;
  //   for(int nj=0;nj<knj;nj++){

  //     c0->cd(++ipad);
  //     MakePad(gPad,0,0); 

  //     if(in==0){
  // 	hbmean[in][nj][0]->Draw("p");
  // 	if(nj==0)drawText2("HYDJET 1.8 (MB)",0.28,0.82,21);
  //     }
  //     else {
  // 	hbmean[in][nj][0]->Draw("psame");
  //     }
  //     if(in==0){
  // 	drawText2("0-100%",0.75,0.82,21);
  // 	drawText2(Form("%s",calgo[nj]),0.28,0.75,21); 
  //     }

  //     c0->cd(ipad+2);
  //     MakePad(gPad,0,0); 
  //     if(in==0){
  // 	hbrms[in][nj][0]->Draw("p");
  // 	if(nj==0)drawText2("HYDJET 1.8 (MB)",0.28,0.82,21);
  //     }
  //     else {
  // 	hbrms[in][nj][0]->Draw("psame");
  //     }
  //     if(in==0){
  // 	drawText2("0-100%",0.75,0.82,21);
  // 	drawText2(Form("%s",calgo[nj]),0.28,0.75,21); 
  //     }
  //   }
  // }
  // c0->cd(1);
  // leg->Draw();
  // return 0;

  TCanvas *c1[ncen];
  for(int nj=0;nj<knj;nj++){
    for(int ic=0;ic<ncen;ic++){      
      c1[ic] = new TCanvas(Form("c1%d",ic),Form("%s 1D Jetbkgd",ccent[ic]),1576,684);
      c1[ic]->Divide(4,2);
      ipad=0;
      for(int ix=0; ix<8; ix++){
	c1[ic]->cd(++ipad);
	MakePad(gPad,0,0); 
	//hjet1D[0][nj][ic][ix]->SetMaximum(1.2);
	hjet1D[0][nj][ic][ix]->SetMaximum(0.243);
	hjet1D[0][nj][ic][ix]->SetMinimum(1e-04);
	hjet1D[0][nj][ic][ix]->Draw("p");
	hjet1D[1][nj][ic][ix]->Draw("ehistsame");
	hjet1D[0][nj][ic][ix]->Draw("psame");

	drawText2(Form("%0.0f < p_{T}^{Jet} (GeV/c) < %0.0f",ptbins[ix],ptbins[ix+1]),0.28,0.80,21); 
	drawText2(calgo[nj],0.28,0.91,21);
	drawText2("|#eta_{jet}| < 2",0.78,0.91,21);
	if(ix==0)drawText2(ccent[ic],0.35,0.65,21);
	if(ipad==1)leg1->Draw();
	if(ipad==2){
	  drawText2(ctrig[isel],0.35,0.65,21);
	}
      }
      //c1[ic]->SaveAs(Form("anplots/jetbkgd_%s_%s.pdf",calgo[nj],mcent[ic]));
      //c1[ic]->SaveAs(Form("anplots/jetbkgd_%s_%s.gif",calgo[nj],mcent[ic]));
    }
  }

  return 0;
}
void GetPileup(TH2 *h2,TH1 *hmean, TH1*hrms)
{
  for(int ix=1;ix<=h2->GetNbinsX();ix++){
    //! pile up pt distribution
    TH1F *hpu = (TH1F*)h2->ProjectionY("hpu",ix,ix);
    hpu->SetName("hpu");
    //hpu->Scale(binwx);
    
    hmean->SetBinContent(ix,hpu->GetMean());
    hmean->SetBinError  (ix,hpu->GetMeanError());
      
    hrms->SetBinContent(ix,hpu->GetRMS());
    hrms->SetBinError  (ix,hpu->GetRMSError());
    delete hpu;
  }//! ix
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

void MakePad(TVirtualPad *cpad, int xlog, int ylog)
{
  cpad->SetLeftMargin(0.22);
  cpad->SetRightMargin(0.05);
  cpad->SetBottomMargin(0.17);
  cpad->SetTopMargin(0.12);

  if(xlog && ylog){
    cpad->SetLogx();
    cpad->SetLogy();
  }else if (xlog){
    cpad->SetLogx();
  }else if (ylog){
    cpad->SetLogy();
  }
}
