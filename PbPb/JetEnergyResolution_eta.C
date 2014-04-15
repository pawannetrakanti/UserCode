
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>

#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
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
#include "TCanvas.h"
#include <TMarker.h>
#include <TString.h>


using namespace std;

const double pi=acos(-1.);
const double pi2=2*pi -1;

const int knj = 12;
const char *calgo[knj] = {"akVs3PF","akVs4PF","akVs5PF",
			  "akVs3Calo","akVs4Calo","akVs5Calo",
			  "ak3PF","ak4PF","ak5PF",
			  "ak3Calo","ak4Calo","ak5Calo"
};
const char *algn [knj] = {"Anti-k_{T}, VS PF, R = 0.3","Anti-k_{T}, VS PF, R = 0.4","Anti-k_{T}, VS PF, R = 0.5",
			  "Anti-k_{T}, VS Calo, R = 0.3","Anti-k_{T}, VS Calo, R = 0.4","Anti-k_{T}, VS Calo, R = 0.5",
			  "Anti-k_{T}, PF, R = 0.3","Anti-k_{T}, PF, R = 0.4""Anti-k_{T}, PF, R = 0.5",
			  "Anti-k_{T}, Calo, R = 0.3","Anti-k_{T}, Calo, R = 0.4","Anti-k_{T}, Calo, R = 0.5",
};

//! Note :  for pp we are using the jet algorithm with out pu subtraction ie. ak3PF, ak4PF and ak5PF

const int ncen=1;
//!                          0 
const char *ccent[ncen] = {"0-100%"};
const char *ksp  [ncen] = {"pbpb"};

//const char *fopt="MLRQ+"; 
int iFit=0; 
const char *fopt="RQ+";
const int knpx=2000;
float fitmin=0.00;
float fitmax=5.00;

double xmin=20.;
double xmax=300.;
double xfitmin=20;
double xfitmax=300;
int maxEntry=5;
int statop=1;


int GetPtBin(float /*pt*/);
void MakeHist(TH1 */*hist*/,int /*istat*/, const char */*xname*/, const char */*yname*/);
void MakeHistRMS(TH1 */*hRMS*/,float /*max*/,float /*min*/);
void MakeHistMean(TH1 */*Mean*/,float /*max*/,float /*min*/);
void CleanHist(TH1 */*h1D*/,float /*lxrange*/,float /*hxrange*/);

void DefineHisto(const char */*cver*/, int /*nj*/, int /*icen*/, int /*ip*/, TH1D *&/*h1*/, const char */*name*/, const char */*caption*/);
void FillMeanSigma(int /*ip*/,TH1 */*h1*/,TH1 */*ArM*/,TH1 */*RMS*/,TH1 */*Mean*/,TH1 */*Sigma*/);

void LoadStyle();
void drawText(const char */*text*/, float /*xp*/, float /*yp*/, int /*size*/);
void drawText2(const char */*text*/, float /*xp*/, float /*yp*/, int /*size*/);
void makeMultiPanelCanvas(TCanvas*& /*canv*/,
                          const Int_t /*columns*/,
                          const Int_t /*rows*/,
                          const Float_t /*leftOffset*/,
                          const Float_t /*bottomOffset*/,
                          const Float_t /*leftMargin*/,
                          const Float_t /*bottomMargin*/,
                          const Float_t /*edge*/,const Float_t /*asyoffset*/); 

void MakeZero(TH1 */*hist*/);


//! Histogram binned
int bins =500;
float minval = 0;
float maxval = 1000;
float binw   = 2;
//! used for fit
double ptbins[]   = {15,30,70,110,150,190,250,300,400,1000};
const int b1  = sizeof(ptbins)/sizeof(Double_t) - 1;
const int nbins = b1;

float minval_eta = -2.5;
float maxval_eta = 2.5;
float binw_eta   = 0.1;
double etabins[] ={-2.0,-1.8,-1.5,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,
		   0,0.2,0.4,0.6,0.8,1.0,1.2,1.5,1.8,2.0};  
const int b2  = sizeof(etabins)/sizeof(Double_t) - 1;
const int neta = b2;

const int kfiles=1;
const char *version[kfiles] = {"53x-80-Track8-Jet22"};
const char *cleg[kfiles]    = {"Track8-Jet22"};
const int icol[kfiles] ={1};
const int isty[kfiles] ={20};
const int jsty[kfiles] ={20};

int JetEnergyResolution_eta(const char *reta = "eta2")
{
  
  LoadStyle();
  
  int rfit=9;

  int kRebin=1;

  float ketacut=2.0;
  if(strcmp(reta,"eta3.0")==0)ketacut=3;
  bool iSigma=false;

  if(rfit==0){fitmin=0.01;fitmax=2.00;}   
  else if(rfit==1){fitmin=0.00;fitmax=1.25;}
  else if(rfit==2){fitmin=0.00;fitmax=1.50;}
  else if(rfit==3){fitmin=0.00;fitmax=1.75;}
  else if(rfit==4){fitmin=0.00;fitmax=2.50;}
  else if(rfit==5){fitmin=0.00;fitmax=3.50;} 
  else if(rfit==6){fitmin=0.50;fitmax=1.25;}
  else if(rfit==7){fitmin=0.50;fitmax=1.50;}
  else if(rfit==8){fitmin=0.50;fitmax=1.75;}
  else if(rfit==9){fitmin=0.50;fitmax=2.50;}
  else if(rfit==10){fitmin=0.50;fitmax=3.50;}
  else if(rfit==11){fitmin=0.75;fitmax=1.25;}
  else if(rfit==12){fitmin=0.75;fitmax=1.50;}
  else if(rfit==13){fitmin=0.75;fitmax=1.75;}
  else if(rfit==14){fitmin=0.75;fitmax=2.50;}
  else if(rfit==15){fitmin=0.75;fitmax=3.50;}
  else if(rfit==16){fitmin=1.00;fitmax=1.25;}
  else if(rfit==17){fitmin=1.00;fitmax=1.50;}
  else if(rfit==18){fitmin=1.00;fitmax=1.75;}
  else if(rfit==19){fitmin=1.00;fitmax=2.50;}
  else if(rfit==20){fitmin=1.00;fitmax=3.50;}

  if(kRebin){
    if(bins%kRebin!=0){
      cout<<"Cannot be divided in these bins chose another combination : "<<endl;
      return 0;
    }
    bins /= kRebin;
    binw = (maxval  - minval)/(1.0*bins);
    cout<<"kRebin : "<<kRebin<<"\t bins : "<<bins<<"\t binw  : "<<binw<<endl;
  }
 
  //! Input files 
  TFile *fin_pbpb[kfiles];

  cout<<"\t"<<endl
      <<"rfit : "<<rfit<<"\t fitmin : "<<fitmin<<"\t fitmax : "<<fitmax<<endl
      <<"\t"<<endl;

  cout<<"# of pt bins : "<<nbins<<endl;
  //return 0;


  xmin=ptbins[0];
  xmax=ptbins[nbins];
  xfitmin=xmin;
  xfitmax=xmax;

  int maxr=2;
  int ipad=0;

  //! 0-100 PbPb Resposnse 
  //! Response corrected pT  (jtpt) / gen pT vs gen pT
  TH1D *hMean[kfiles][knj][ncen][nbins], *hSigma[kfiles][knj][ncen][nbins], *hRMS[kfiles][knj][ncen][nbins], *hArM[kfiles][knj][ncen][nbins];

  //! Response raw pT / gen pT  vs gen pT
  TH1D *hMean_r[kfiles][knj][ncen][nbins], *hSigma_r[kfiles][knj][ncen][nbins], *hRMS_r[kfiles][knj][ncen][nbins], *hArM_r[kfiles][knj][ncen][nbins];

  TFile *fout = new TFile("output.root","RECREATE");

  for(int in=0; in<kfiles; in++){
    fin_pbpb[in] = new TFile("input/JetResponse_histos_pp.root","r");
    for(int nj=0;nj<knj;nj++){

      cout<<"nj : "<<nj<<Form("\t %s",calgo[nj])<<endl;

      for(int icen=0;icen<ncen;icen++){
	
	//cout<<"icen : "<<ccent[icen]<<endl;
	
	TH3F *hreco = (TH3F*)fin_pbpb[in]->Get(Form("hrescreta_genm%d_%d",nj,icen));
	TH3F *hraw  = (TH3F*)fin_pbpb[in]->Get(Form("hresrreta_genm%d_%d",nj,icen));

	for(int ip=0; ip<nbins; ip++){	  	

	  DefineHisto(cleg[in], nj, icen, ip, hMean   [in][nj][icen][ip],"hMean"    ,"<reco p_{T}/gen p_{T}>");
	  DefineHisto(cleg[in], nj, icen, ip, hMean_r [in][nj][icen][ip],"hMean_r"  ,"<reco p_{T}/gen p_{T}>");
	  DefineHisto(cleg[in], nj, icen, ip, hArM    [in][nj][icen][ip],"hArm"     ,"<reco p_{T}/gen p_{T}>");
	  DefineHisto(cleg[in], nj, icen, ip, hArM_r  [in][nj][icen][ip],"hArm_r"   ,"<reco p_{T}/gen p_{T}>");
	  DefineHisto(cleg[in], nj, icen, ip, hSigma  [in][nj][icen][ip],"hSigma"   ,"#sigma(reco p_{T}/gen p_{T})");
	  DefineHisto(cleg[in], nj, icen, ip, hSigma_r[in][nj][icen][ip],"hSigma_r" ,"#sigma(reco p_{T}/gen p_{T})");
	  DefineHisto(cleg[in], nj, icen, ip, hRMS    [in][nj][icen][ip],"hRMS"     ,"#sigma(reco p_{T}/gen p_{T})");
	  DefineHisto(cleg[in], nj, icen, ip, hRMS_r  [in][nj][icen][ip],"hRMS_r"   ,"#sigma(reco p_{T}/gen p_{T})");

	  int lbin = hreco->GetYaxis()->FindBin(ptbins[ip]);  //(int)(ptbins[ip]   - minval)/binw + 1;
	  int hbin = hreco->GetYaxis()->FindBin(ptbins[ip+1]);//(int)(ptbins[ip+1] - minval)/binw + 1;

	  //if(nj==0)cout<<" pt : " << ptbins[ip] << "\t lbin : "<< lbin << " " << hreco->GetYaxis()->FindBin(ptbins[ip]) << " \t hbin : " << hbin << " " << hreco->GetYaxis()->FindBin(ptbins[ip+1]) << endl;

	  for(int ie=0; ie<neta; ie++){

	    int leta = hreco->GetXaxis()->FindBin(etabins[ie]);  //(int)(etabins[ie]   - minval_eta)/binw_eta + 1;
	    int heta = hreco->GetXaxis()->FindBin(etabins[ie+1]);//(int)(etabins[ie+1] - minval_eta)/binw_eta + 1;
	    
	    TH1D *h1reco = (TH1D*)hreco->ProjectionZ(Form("hpz_reco_%d_%d_%d",nj,ip,ie),leta,heta,lbin,hbin,"e");
	    // if(nj==1){
	    //   //TCanvas *c0 = new TCanvas("c0","",1);
	    //   //c0->cd();
	    //   //h1reco->SetMarkerStyle(20);
	    //   //h1reco->Draw("p");
	    //   //c0->SaveAs("hist1d_ak3VsCalo_0_0.gif");
	    //   fout->cd();
	    //   h1reco->Write();
	    //   //return 0;
	    // }
	    FillMeanSigma(ie, h1reco, hArM[in][nj][icen][ip], hRMS[in][nj][icen][ip], hMean[in][nj][icen][ip], hSigma[in][nj][icen][ip]);

	    TH1D *h1raw = (TH1D*)hraw->ProjectionZ("hpz_raw",leta,heta,lbin,hbin,"e");
	    FillMeanSigma(ie, h1raw, hArM_r[in][nj][icen][ip], hRMS_r[in][nj][icen][ip], hMean_r[in][nj][icen][ip], hSigma_r[in][nj][icen][ip]);

	    delete h1reco;
	    delete h1raw;
	  }//! ie
	}//! ip
	delete hreco;
	delete hraw;
      }//! icen loop ends
    }//! nj loop ends
  }//! kfiles version 
  //return 0;
  //fout->Close();

  for(int in=0; in<kfiles; in++){
    for(int nj=0;nj<knj;nj++){
      for(int ic=0;ic<ncen;ic++){
	for(int ip=0; ip<nbins; ip++){	  	
	  MakeHistMean(hArM  [in][nj][ic][ip],1.15,0.75);
	  MakeHistMean(hMean [in][nj][ic][ip],1.15,0.75);
	  
	  MakeHistRMS (hRMS  [in][nj][ic][ip],0.43,0.001);
	  MakeHistRMS (hSigma[in][nj][ic][ip],0.43,0.001);
	}      
      }
    }
  }

  TCanvas *c3[4][nbins];
  TLegend *l4 = new TLegend(0.384227,0.7972266,0.8904035,0.967033,NULL,"BRNDC");
  l4->SetHeader("CMSSW 5_3_16");
  l4->SetBorderSize(0);
  l4->SetTextFont(42);
  l4->SetTextSize(0.07);
  l4->SetLineColor(1);
  l4->SetLineStyle(1);
  l4->SetLineWidth(1);
  l4->SetFillColor(10);
  l4->SetFillStyle(1001);
  l4->AddEntry(hRMS[0][0][0][0],cleg[0],"p"); 


  TLine *line = new TLine(etabins[0],1.00,etabins[b2],1.00);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->SetLineColor(1);

  int maxc=3;
  maxr=2;
  const char *alg[4] = {"VsPF","VsCalo","PF","Calo"};
  int plot[knj]={0,1,2,3,4,5,6,7,8,9,10,11};

  ipad=0;
  for(int in=0; in<kfiles; in++){
    ipad=0;

    for(int ic=0;ic<ncen;ic++){
      for(int ip=0; ip<nbins;ip++){      
	ipad=0;
	int ik=-1;
	for(int nj=0; nj<knj;nj++){
	  if(nj%3==0){
	    ++ik;
	    c3[ik][ip] = new TCanvas(Form("c3_%d_%d",ik,ip),Form("%s JES JER",alg[ik]),1211,652);
	    makeMultiPanelCanvas(c3[ik][ip],maxc,maxr,0.0,0.0,0.22,0.22,0.02,0);
	    ipad=0;
	  }
	  hRMS[in][plot[nj]][ic][ip]->SetMarkerStyle(isty[in]);
	  hRMS[in][plot[nj]][ic][ip]->SetMarkerColor(icol[in]);
	  hRMS[in][plot[nj]][ic][ip]->SetLineColor(icol[in]);
	  hRMS[in][plot[nj]][ic][ip]->SetMarkerSize(1.3);

	  hArM[in][plot[nj]][ic][ip]->SetMarkerStyle(isty[in]);
	  hArM[in][plot[nj]][ic][ip]->SetMarkerColor(icol[in]);
	  hArM[in][plot[nj]][ic][ip]->SetLineColor(icol[in]);
	  hArM[in][plot[nj]][ic][ip]->SetMarkerSize(1.3);

	  hArM_r[in][plot[nj]][ic][ip]->SetMarkerStyle(jsty[in]);
	  hArM_r[in][plot[nj]][ic][ip]->SetMarkerColor(icol[in]);
	  hArM_r[in][plot[nj]][ic][ip]->SetLineColor(icol[in]);
	  hArM_r[in][plot[nj]][ic][ip]->SetMarkerSize(1.3);

	  hSigma[in][plot[nj]][ic][ip]->SetMarkerStyle(isty[in]);
	  hSigma[in][plot[nj]][ic][ip]->SetMarkerColor(icol[in]);
	  hSigma[in][plot[nj]][ic][ip]->SetLineColor(icol[in]);
	  hSigma[in][plot[nj]][ic][ip]->SetMarkerSize(1.3);

	  hMean[in][plot[nj]][ic][ip]->SetMarkerStyle(isty[in]);
	  hMean[in][plot[nj]][ic][ip]->SetMarkerColor(icol[in]);
	  hMean[in][plot[nj]][ic][ip]->SetLineColor(icol[in]);
	  hMean[in][plot[nj]][ic][ip]->SetMarkerSize(1.3);

	  hMean_r[in][plot[nj]][ic][ip]->SetMarkerStyle(jsty[in]);
	  hMean_r[in][plot[nj]][ic][ip]->SetMarkerColor(icol[in]);
	  hMean_r[in][plot[nj]][ic][ip]->SetLineColor(icol[in]);
	  hMean_r[in][plot[nj]][ic][ip]->SetMarkerSize(1.3);

	  c3[ik][ip]->cd(++ipad);

	  //cout<<" nj : " << nj << "\t ip : "<< ip << "\t ipad : "<< ipad << "\t ik : " <<ik <<endl;
	
	  hRMS[in][plot[nj]][ic][ip]->SetMaximum(0.532);
	  hSigma[in][plot[nj]][ic][ip]->SetMaximum(0.532);
	  if(in==0)hRMS[in][plot[nj]][ic][ip]->Draw("p");
	  else hRMS[in][plot[nj]][ic][ip]->Draw("psame");
	  //if(in==0)hSigma[in][plot[nj]][ic][ip]->Draw("p");
	  //else hSigma[in][plot[nj]][ic][ip]->Draw("psame");

	  if(ipad!=1){
	    MakeZero(hRMS[in][plot[nj]][ic][ip]);
	    MakeZero(hSigma[in][plot[nj]][ic][ip]);
	    MakeZero(hMean[in][plot[nj]][ic][ip]);
	    MakeZero(hArM[in][plot[nj]][ic][ip]);
	  }
	  if(ipad==1){
	    drawText2("CMS Simulation", 0.28, 0.86, 22);
	    //drawText2("HYDJET 1.8",0.28,0.75,21);
	    drawText2("PYTHIA Z2",0.28,0.75,21);
	  }
	  else if(ipad==2){
	    drawText2("#sqrt{s_{NN}} = 2.76 TeV",0.25,0.85,21);
	    drawText2(Form("%0.0f < p_{T}^{gen} (GeV/c) < %0.0f",ptbins[ip],ptbins[ip+1]),0.15,0.72,21);
	  }
	  drawText2(Form("%s",calgo[plot[nj]]),0.28,0.60,21);  

	  c3[ik][ip]->cd(ipad+3);

	  hArM[in][plot[nj]][ic][ip]->SetMaximum(1.08);      
	  hArM[in][plot[nj]][ic][ip]->SetMinimum(0.92);      
	  
	  hMean[in][plot[nj]][ic][ip]->SetMaximum(1.08);      
	  hMean[in][plot[nj]][ic][ip]->SetMinimum(0.92);      
	  
	  if(in==0){
	    hArM[in][plot[nj]][ic][ip]->Draw("p");
	    //hArM_r[in][plot[nj]][ic][ip]->Draw("psame");
	  }else{
	    hArM[in][plot[nj]][ic][ip]->Draw("psame");
	    //hArM_r[in][plot[nj]][ic][ip]->Draw("psame");
	  }
	  
	  // if(in==0){
	  //   hMean[in][plot[nj]][ic][ip]->Draw("p");
	  //   //hMean_r[in][plot[nj]][ic][ip]->Draw("psame");
	  // }else{
	  //   hMean[in][plot[nj]][ic][ip]->Draw("psame");
	  //   //hMean_r[in][plot[nj]][ic][ip]->Draw("psame");
	  // }
	  if(ipad==2)l4->Draw();
	  line->Draw();
	  //c3[ik][ip]->SaveAs(Form("JESJER_etadep_%0.0fpT%0.0f_%s.gif",ptbins[ip],ptbins[ip+1],alg[ik]));
	}//! ip
      }//! icen
    }//! knj
  }//! in


  return 0;
  
}
void DefineHisto(const char *cver, int nj, int ic, int ip, TH1D *&h1, const char *name, const char *caption)
{
  h1 = new TH1D(Form("%s_%s_%d_%d_%d",name,cver,nj,ic,ip),Form("CMSSW %s %s %s %s %d %d",cver,caption,ksp[ic],calgo[nj],ic,ip),neta,etabins);
  MakeHist(h1,statop,"#eta","<p_{T}^{RecoJet}/p_{T}^{GenJet}>");
}
void MakeHist(TH1 *histo,int istat,const char *xname, const char *yname)
{
  histo->SetStats(istat);
  histo->SetMarkerStyle(24);
  histo->SetMarkerColor(1);
  histo->SetLineColor(1);
  histo->SetLineStyle(1);
  histo->GetXaxis()->SetTitle(xname);
  histo->GetXaxis()->CenterTitle(true);
  histo->GetYaxis()->SetTitle(yname);
  histo->GetYaxis()->CenterTitle(true);
}
void FillMeanSigma(int ip,TH1 *h1F,TH1 *harm,TH1 *hrms,TH1 *hm,TH1 *hs)
{

  TF1 *f1 = new TF1("f1","(([0]/(2*pi*[1]*[1]))*exp(-1.*((x-[2])*(x-[2])/(2*[1]*[1]))))",0,5);
  f1->SetParameters(1,0.1,1);
  f1->SetParNames("A","#sigma","mean");
  f1->SetParLimits(0,0,20);     //! A
  f1->SetParLimits(1,0.,10.0);  //! Sigma
  f1->SetParLimits(2,fitmin,fitmax); //! Mean
  f1->SetLineWidth(1);
  f1->SetNpx(knpx);
  
  float mm=0,ss=0,p0=0;  
  if(h1F->GetEntries()<maxEntry){
    h1F->Scale(0.);
    harm  ->SetBinContent(ip+1,-9999);
    harm  ->SetBinError  (ip+1,0);
    hrms  ->SetBinContent(ip+1,-9999);
    hrms  ->SetBinError  (ip+1,0);
  }
  if(h1F->Integral()>0){

    h1F->Scale(1./h1F->Integral());
    if(iFit==0){
      //h1F->Fit("gaus",fopt,"",fitmin,fitmax);
      h1F->Fit("gaus","MRQ0+","",fitmin,fitmax);
      TF1* f2 = (TF1*)h1F->GetFunction("gaus");
      f2->SetLineWidth(1);
      f2->SetLineStyle(2);
      f2->SetNpx(knpx);
      hm->SetBinContent(ip+1,f2->GetParameter(1));
      hs->SetBinContent(ip+1,f2->GetParameter(2));
      
      if(strcmp(fopt,"RQ0+")==0){
	hm->SetBinError  (ip+1,h1F->GetMeanError());
	hs->SetBinError  (ip+1,h1F->GetRMSError());
      }else{
	hm->SetBinError(ip+1,f2->GetParError(1));
	hs->SetBinError(ip+1,f2->GetParError(2));
      }
    }else{
      mm = h1F->GetMean();
      ss = h1F->GetRMS();
      p0 = h1F->GetMaximum();
      f1->SetParameters(p0,ss,mm);
      f1->SetParLimits(0,0,2*p0);
      f1->SetParLimits(1,0,2*ss);
      f1->SetParLimits(2,fitmin,fitmax);
      //f1->SetParLimits(2,mm-2.5*ss,mm+2.5*ss);
      
      //h1F->Fit("f1",fopt,"",fitmin,fitmax);
      h1F->Fit("f1","RQ0+","",fitmin,fitmax);
      hm->SetBinContent(ip+1,f1->GetParameter(2));
      hs->SetBinContent(ip+1,f1->GetParameter(1));
      
      if(strcmp(fopt,"RQ0+")==0){
	hm->SetBinError  (ip+1,h1F->GetMeanError());
	hs->SetBinError  (ip+1,h1F->GetRMSError());
      }else{
	hm->SetBinError  (ip+1,f1->GetParError(2));
	hs->SetBinError  (ip+1,f1->GetParError(1));
      }
    }
    harm->SetBinContent(ip+1,h1F->GetMean());
    harm->SetBinError  (ip+1,h1F->GetMeanError());
    hrms->SetBinContent(ip+1,h1F->GetRMS());
    hrms->SetBinError  (ip+1,h1F->GetRMSError());
  }
  delete f1;
}
void MakeHistRMS(TH1 *h1,float ymax,float ymin)
{

  h1->SetTitle("");
  h1->SetMaximum(ymax);
  h1->SetMinimum(ymin);
  h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->GetXaxis()->SetTitle("GenJet #eta");
  h1->GetXaxis()->CenterTitle(true);
  h1->GetXaxis()->SetMoreLogLabels();
  h1->GetXaxis()->SetNoExponent();
  h1->GetXaxis()->SetNdivisions(507);
  h1->GetXaxis()->SetLabelFont(42);
  h1->GetXaxis()->SetLabelOffset(0.01);
  h1->GetXaxis()->SetLabelSize(0.07);
  h1->GetXaxis()->SetTitleSize(0.07);
  h1->GetXaxis()->SetTitleOffset(1.15);
  h1->GetXaxis()->SetTitleFont(42);
  h1->GetYaxis()->SetTitle("#sigma (RecoJet p_{T} / GenJet p_{T})");
  h1->GetYaxis()->CenterTitle(true);
  h1->GetYaxis()->SetNdivisions(507);
  h1->GetYaxis()->SetLabelFont(42);
  h1->GetYaxis()->SetLabelOffset(0.01);
  h1->GetYaxis()->SetLabelSize(0.09);
  h1->GetYaxis()->SetTitleSize(0.09);
  h1->GetYaxis()->SetTitleOffset(1.12);
  h1->GetYaxis()->SetTitleFont(42);
  h1->GetYaxis()->SetDecimals(true);

}
void MakeZero(TH1 *h1)
{
  h1->GetYaxis()->SetLabelSize(0);
  h1->GetYaxis()->SetTitleSize(0);
}
void MakeHistMean(TH1 *h1,float ymax,float ymin)
{
  h1->SetMaximum(ymax);
  h1->SetMinimum(ymin);
  h1->SetTitle("");
  h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->GetXaxis()->SetTitle("GenJet #eta");
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
  h1->GetYaxis()->SetTitle("<RecoJet p_{T} / GenJet p_{T}>");
  h1->GetYaxis()->SetTitleSize(0.07);
  h1->GetYaxis()->SetTitleOffset(1.50);
  h1->GetYaxis()->SetLabelSize(0.07);
  h1->GetYaxis()->SetNdivisions(507);
  h1->GetYaxis()->SetDecimals(true);
  h1->GetYaxis()->SetTitleFont(42);
  h1->GetYaxis()->SetLabelFont(42);
}
int GetPtBin(float pt)
{
  int ibin=-1;
  ibin = (int)(pt + minval)/binw;
  return ibin;
}
void CleanHist(TH1 *h1F,float lxval, float hxval)
{
  for(int ix=1;ix<=h1F->GetNbinsX();ix++){
    double val = h1F->GetBinCenter(ix);
    if(val<lxval || val>hxval){
      h1F->SetBinContent(ix,0);
      h1F->SetBinError(ix,0);
    }
  }
}
void makeMultiPanelCanvas(TCanvas*& canv,
                          const Int_t columns,
                          const Int_t rows,
                          const Float_t leftOffset,
                          const Float_t bottomOffset,
                          const Float_t leftMargin,
                          const Float_t bottomMargin,
                          const Float_t edge,const Float_t asyoffset) {
  if (canv==0) {
    cout<<"makeMultiPanelCanvas :  Got null canvas."<<endl;
    return;
  }
  canv->Clear();
  
  TPad* pad[columns][rows];
  
  Float_t Xlow[columns];
  Float_t Xup[columns];
  Float_t Ylow[rows];
  Float_t Yup[rows];
  
  Float_t PadWidth =
    (1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
		      (1.0/(1.0-edge))+(Float_t)columns-2.0);
  Float_t PadHeight =
    (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
			(1.0/(1.0-edge))+(Float_t)rows-2.0);
  
  Xlow[0] = leftOffset - asyoffset;
  Xup[0]  = leftOffset + PadWidth/(1.0-leftMargin);
  Xup[columns-1] = 1;
  Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);
  
  Yup[0] = 1;
  Ylow[0] = 1.0-PadHeight/(1.0-edge);
  Ylow[rows-1] = bottomOffset;
  Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);
  
  for(Int_t i=1;i<columns-1;i++) {
    Xlow[i] = Xup[0] + (i-1)*PadWidth;
    Xup[i] = Xup[0] + (i)*PadWidth;
  }
  Int_t ct = 0;
  for(Int_t i=rows-2;i>0;i--) {
    Ylow[i] = Yup[rows-1] + ct*PadHeight;
    Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
    ct++;
  }
  
  TString padName;
  for(Int_t i=0;i<columns;i++) {
    for(Int_t j=0;j<rows;j++) {
      canv->cd();
      padName = Form("p_%d_%d",i,j);
      pad[i][j] = new TPad(padName.Data(),padName.Data(),
			   Xlow[i],Ylow[j],Xup[i],Yup[j]);
      
      // this is hacked version to create aysmmetric pads around low 
      if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
      else pad[i][j]->SetLeftMargin(0);
      
      if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
      else pad[i][j]->SetRightMargin(0);
      
      if(j==0) pad[i][j]->SetTopMargin(edge);
      else pad[i][j]->SetTopMargin(0);
      
      if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
      else pad[i][j]->SetBottomMargin(0);
      
      pad[i][j]->Draw();
      pad[i][j]->cd();
      pad[i][j]->SetNumber(columns*j+i+1);
    }
  }
}
void LoadStyle()
{
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetErrorX(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetFillColor(10);
  gStyle->SetFrameFillColor(10);
  gStyle->SetPadColor(10);
  gStyle->SetStatColor(10);
  gStyle->SetTitleFillColor(10);
  gStyle->SetPadBorderSize(0);
  gStyle->SetLegendBorderSize(1);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetPalette(1);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
}
void drawText(const char *text, float xp, float yp, int size){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(63);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  //tex->SetTextFont(42);
  tex->SetNDC();
  tex->Draw();
}

void drawText2(const char *text, float xp, float yp, int size){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(43);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw();
}
