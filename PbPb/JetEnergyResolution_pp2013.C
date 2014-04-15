
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

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
#include "TCanvas.h"
#include <TMarker.h>
#include <TString.h>
#include <TVirtualFitter.h>

using namespace std;

const double pi=acos(-1.);
const double pi2=2*pi -1;

int bins =500;
float minval = 0;
float maxval = 1000;
float binw   = 2;

const int knj = 6;
const char *calgo[knj] = {"ak3PF","ak4PF","ak5PF",
			  "ak3Calo","ak4Calo","ak5Calo"
};
const char *algn [knj] = {"Anti-k_{T}, PF, R = 0.3","Anti-k_{T}, PF, R = 0.4""Anti-k_{T}, PF, R = 0.5",
			  "Anti-k_{T}, Calo, R = 0.3","Anti-k_{T}, Calo, R = 0.4","Anti-k_{T}, Calo, R = 0.5"
};
//! Note :  for pp we are using the jet algorithm with out pu subtraction ie. ak3PF, ak4PF and ak5PF

const int ktype=2;
const char *alg[ktype] = {"PF","Calo"};
int plot[knj]={0,1,2,3,4,5};

const int ncen=1;
//!                          0 
const char *ccent[ncen] = {"0-100%"};
const char *ksp  [ncen] = {"pp"};

//const char *fopt="MLLRQ0+"; 
const char *fopt="RQ+";
int iFit=1; 

const int knpx=2000;
float fitmin=0.00;
float fitmax=5.00;

double xmin=15.;
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

void DefineHisto(const char */*cver*/, int /*nj*/, int /*icen*/, TH1F *&/*h1*/, const char */*name*/, const char */*caption*/);
void FillMeanSigma(int /*ip*/,TH1F *&/*h1*/,TH1F */*ArM*/,TH1F */*RMS*/,TH1F */*Mean*/,TH1F */*Sigma*/);
void CalJesJer(const char */*cver*/, int /*nj*/, int /*icen*/, TH2F */*h2d*/, TH1F **/*h1d*/, TH1F */*ham*/, TH1F */*hr*/, TH1F */*hm*/, TH1F */*hs*/);

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

//! Fit functions
void fit_gaussian(TH1F *&/*hrsp*/, const double /*nsigma*/, const int /*niter*/);
int  fit_dscb(TH1F *&/*hrsp*/, const double /*nsigma*/, const int /*niter*/);

/// double sided crystal ball function definition
double fnc_dscb(double */*xx*/, double */*pp*/);
void adjust_fitrange(TH1F */*h*/, double &/*min*/, double &/*max*/);


void MakeZero(TH1 */*hist*/);


//! used for fit
//double ptbins[] = {80,90,100,110,120,130,140,150,160,170,180,190,200,286};
//double ptbins[] = {18,28,38,48,58,68,78,88,98,108,118,128,138,148,158,168,178,188,198,208,248,288,308,358,408,458,508,648};
//double ptbins[] = {14,18,22,26,30,34,38,42,44,48,52,56,60};
double ptbins[]   = {15,20,25,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,220,240,260,280,300,320,340,360,400,450,500,550,600,848};
//double ptbins[]   = {20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,200,220,240,260,300,386};
const int b1  = sizeof(ptbins)/sizeof(Double_t) - 1;
const int nbins = b1;

const char *version  = "prod22 v81";
const char *cleg  = "prod22 v81";
const int icol =1;
const int isty =20;
const int jsty =20;

//TFile *fout=0;

int JetEnergyResolution_pp2013(const char *reta = "eta2")
{
  LoadStyle();

  std::string inname="input/JetResponse_histos_prod22v81_pp.root";  
  int rfit=21;
  int kRebin=1;

  float ketacut=2.0;
  if(strcmp(reta,"eta3.0")==0)ketacut=3;
  bool iSigma=false;

  if(iSigma){
    gStyle->SetOptStat(1);
    gStyle->SetOptFit(11111);
  }

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
  else if(rfit==21){fitmin=0.20;fitmax=1.70;}


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
  TFile *fin_pp=0;

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

  //! 0-100 Pp Resposnse 
  //! Response corrected pT  (jtpt) / gen pT vs gen pT
  TH1F *hgenpt[knj][ncen], *hrecpt[knj][ncen], *hrawpt[knj][ncen];
  
  
  TH2F *hratiocorrrefpt[knj][ncen];
  TH1F *hratiocorrrefpt1D[knj][ncen][nbins];  
  TH1F *hMean[knj][ncen], *hSigma[knj][ncen], *hRMS[knj][ncen], *hArM[knj][ncen];

  //! Response raw pT / gen pT  vs gen pT
  TH2F *hratiorawrefpt[knj][ncen];
  TH1F *hratiorawrefpt1D[knj][ncen][nbins];  
  TH1F *hMean_r[knj][ncen], *hSigma_r[knj][ncen], *hRMS_r[knj][ncen], *hArM_r[knj][ncen];

  //! Ratio of Fit and geometric mean/RMS
  //TH1F *hRatio_Mean[knj][ncen], *hRatio_RMS[knj][ncen];

  TH2F *hgenjrecoj[knj][ncen];

  fin_pp = new TFile(inname.c_str(),"r"); 

  for(int nj=0;nj<knj;nj++){
    //cout<<"nj : "<<nj<<Form("\t %s",calgo[nj])<<endl;
    for(int icen=0;icen<ncen;icen++){
      
      DefineHisto(cleg, nj, icen, hMean   [nj][icen],"hMean"    ,"<reco p_{T}/gen p_{T}>");
      DefineHisto(cleg, nj, icen, hMean_r [nj][icen],"hMean_r"  ,"<reco p_{T}/gen p_{T}>");
      DefineHisto(cleg, nj, icen, hArM    [nj][icen],"hArm"     ,"<reco p_{T}/gen p_{T}>");
      DefineHisto(cleg, nj, icen, hArM_r  [nj][icen],"hArm_r"   ,"<reco p_{T}/gen p_{T}>");
      DefineHisto(cleg, nj, icen, hSigma  [nj][icen],"hSigma"   ,"#sigma(reco p_{T}/gen p_{T})");
      DefineHisto(cleg, nj, icen, hSigma_r[nj][icen],"hSigma_r" ,"#sigma(reco p_{T}/gen p_{T})");
      DefineHisto(cleg, nj, icen, hRMS    [nj][icen],"hRMS"     ,"#sigma(reco p_{T}/gen p_{T})");
      DefineHisto(cleg, nj, icen, hRMS_r  [nj][icen],"hRMS_r"   ,"#sigma(reco p_{T}/gen p_{T})");
      
      //cout<<"icen : "<<ccent[icen]<<endl;
      //! Pp  starts here
      hgenjrecoj[nj][icen] = (TH2F*)fin_pp->Get(Form("hgenjrecoj%d_%d",nj,icen));
      
      hratiocorrrefpt[nj][icen] = (TH2F*)fin_pp->Get(Form("hrescrpt_genm%d_%d",nj,icen));
      hratiocorrrefpt[nj][icen]->SetName(Form("hratiocorrrefpt_pp_%s_%d_%d",version,nj,icen));
      hratiorawrefpt [nj][icen] = (TH2F*)fin_pp->Get(Form("hresrrpt_genm%d_%d",nj,icen));
      hratiorawrefpt [nj][icen]->SetName(Form("hratiorawrefpt_pp_%s_%d_%d",version,nj,icen));
      
      //if(nj==1)cout<<nj<<"\t icen : "<<icen<<"\t hMean : "<<hMean[nj][icen]->GetName()<<endl;
      //! Calculate jes and jer
      //! Raw / Gen
      CalJesJer(version, nj, icen, 
		hratiorawrefpt [nj][icen], hratiorawrefpt1D [nj][icen], 
		hArM_r[nj][icen], hRMS_r[nj][icen], 
		hMean_r[nj][icen],hSigma_r[nj][icen]);
      
      //! Reco / Gen
      CalJesJer(version, nj, icen, 
		hratiocorrrefpt[nj][icen], hratiocorrrefpt1D[nj][icen], 
		hArM[nj][icen], hRMS[nj][icen], 
		hMean[nj][icen],hSigma[nj][icen]);
      
      
      hgenpt[nj][icen] = (TH1F*)fin_pp->Get(Form("hgenpt_genm%d_%d",nj,icen));
      hgenpt[nj][icen]->SetMarkerStyle(24);
      hgenpt[nj][icen]->SetMarkerSize(1.2);
      hgenpt[nj][icen]->SetMarkerColor(1);
      hgenpt[nj][icen]->Rebin(5);
      hgenpt[nj][icen]->SetTitle(Form("%s",calgo[nj]));
      hgenpt[nj][icen]->GetXaxis()->SetRangeUser(ptbins[0],ptbins[b1]);
      hgenpt[nj][icen]->GetXaxis()->SetTitle("jet p_{T} (GeV/c)");
      hgenpt[nj][icen]->GetXaxis()->CenterTitle();
      hgenpt[nj][icen]->GetXaxis()->SetTitleSize(0.07);
      hgenpt[nj][icen]->GetXaxis()->SetLabelSize(0.07);
      hgenpt[nj][icen]->GetYaxis()->SetRangeUser(0.1,6e+02);
      hgenpt[nj][icen]->GetYaxis()->SetTitle("Counts");
      hgenpt[nj][icen]->GetYaxis()->CenterTitle();
      hgenpt[nj][icen]->GetYaxis()->SetTitleSize(0.07);
      hgenpt[nj][icen]->GetYaxis()->SetLabelSize(0.07);
      
      hrecpt[nj][icen] = (TH1F*)fin_pp->Get(Form("hrecopt_genm%d_%d",nj,icen));
      hrecpt[nj][icen]->SetMarkerStyle(20);
      hrecpt[nj][icen]->SetMarkerSize(1.2);
      hrecpt[nj][icen]->SetMarkerColor(2);
      hrecpt[nj][icen]->Rebin(5);
      //hrecpt[nj][icen]->GetXaxis()->SetRangeUser(ptbins[0],ptbins[b1]);
      
      hrawpt[nj][icen] = (TH1F*)fin_pp->Get(Form("hrawpt_genm%d_%d",nj,icen));
      hrawpt[nj][icen]->SetMarkerStyle(25);
      hrawpt[nj][icen]->SetMarkerSize(1.2);
      hrawpt[nj][icen]->SetMarkerColor(4);
      hrawpt[nj][icen]->Rebin(5);
      //hrawpt[nj][icen]->GetXaxis()->SetRangeUser(ptbins[0],ptbins[b1]);
      
      //! Pp ends
    }//! icen loop ends
  }//! nj loop ends
  //return 0;

  for(int nj=0;nj<knj;nj++){
    for(int ic=0;ic<ncen;ic++){
      //MakeHistMean(hArM  [nj][ic],1.19,0.84);
      //MakeHistMean(hMean [nj][ic],1.19,0.84);
      MakeHistMean(hArM  [nj][ic],1.05,0.95);
      MakeHistMean(hMean [nj][ic],1.05,0.95);
      
      MakeHistRMS (hRMS  [nj][ic],0.43,0.001);
      MakeHistRMS (hSigma[nj][ic],0.43,0.001);
    }      
  }

  
  
   ipad=0;
   TCanvas *c99[knj];
   for(int nj=0;nj<1;nj++){
     c99[nj] = new TCanvas(Form("c99_%d",nj),Form("%s Fitting plots",calgo[nj]),100,102,1399,942);
     c99[nj]->Divide(4,3,0,0);
     ipad=0;
  
     for(int ip=0;ip<12/*nbins*/;ip++){      
       c99[nj]->cd(++ipad);
       if(ipad%4==0)gPad->SetRightMargin(0.02);
       if(ipad==1 || ipad==5 || ipad==9)gPad->SetLeftMargin(0.15);
       gPad->SetBottomMargin(0.15);
       //gPad->SetLogy();
    
       //hratiocorrrefpt1D[nj][0][ip]->SetMaximum(25.634);
       hratiocorrrefpt1D[nj][0][ip]->SetMaximum(0.164);
       hratiocorrrefpt1D[nj][0][ip]->SetMinimum(1e-09);
       hratiocorrrefpt1D[nj][0][ip]->SetTitle(0);
       hratiocorrrefpt1D[nj][0][ip]->GetXaxis()->SetTitle("<reco jet p_{T} / gen jet p_{T}>");
       hratiocorrrefpt1D[nj][0][ip]->GetXaxis()->SetTitleFont(42);
       hratiocorrrefpt1D[nj][0][ip]->GetXaxis()->SetLabelFont(42);
       hratiocorrrefpt1D[nj][0][ip]->GetXaxis()->SetLabelSize(0.08);
       hratiocorrrefpt1D[nj][0][ip]->GetXaxis()->SetTitleSize(0.07);
       hratiocorrrefpt1D[nj][0][ip]->GetYaxis()->SetTitle("");
       hratiocorrrefpt1D[nj][0][ip]->GetYaxis()->SetTitleFont(42);
       hratiocorrrefpt1D[nj][0][ip]->GetYaxis()->SetLabelFont(42);
       hratiocorrrefpt1D[nj][0][ip]->GetYaxis()->SetLabelSize(0.08);
       hratiocorrrefpt1D[nj][0][ip]->SetMarkerStyle(24);
       hratiocorrrefpt1D[nj][0][ip]->SetMarkerColor(1);
       hratiocorrrefpt1D[nj][0][ip]->SetLineColor(1);
       hratiocorrrefpt1D[nj][0][ip]->SetMarkerSize(1.1);
       hratiocorrrefpt1D[nj][0][ip]->Draw("p");  

       if(ipad==1){drawText2(Form("%s",calgo[nj]),0.28,0.90,21);      
	 drawText2(Form("%0.0f < p_{T} (GeV/c) < %0.0f",ptbins[ip], ptbins[ip+1]),0.22,0.80,21);		
       }else 	 drawText2(Form("%0.0f < p_{T} (GeV/c) < %0.0f",ptbins[ip], ptbins[ip+1]),0.17,0.80,21);		
     }
   }
   //return 0;



  TLegend *l3 = new TLegend(0.05976299,0.2574956,0.3741464,0.5714286,NULL,"BRNDC");
  l3->SetHeader("");
  l3->SetBorderSize(0);
  l3->SetTextFont(42);
  l3->SetTextSize(0.09);
  l3->SetLineColor(1);
  l3->SetLineStyle(1);
  l3->SetLineWidth(1);
  l3->SetFillColor(10);
  l3->SetFillStyle(1001);
  l3->AddEntry(hgenpt[0][0],"Gen","p"); 
  l3->AddEntry(hrecpt[0][0],"Reco","p"); 
  l3->AddEntry(hrawpt[0][0],"Raw","p"); 


  //TLegend *l4 = new TLegend(0.2527884,0.6276059,0.5060804,0.963133,NULL,"BRNDC");
  TLegend *l4 = new TLegend(0.2419691,0.8166283,0.8611024,0.9792,NULL,"BRNDC");
  l4->SetHeader("");
  l4->SetBorderSize(0);
  l4->SetTextFont(42);
  l4->SetTextSize(0.07);
  l4->SetLineColor(1);
  l4->SetLineStyle(1);
  l4->SetLineWidth(1);
  l4->SetFillColor(10);
  l4->SetFillStyle(1001);
  l4->AddEntry(hRMS[0][0],"PYTHIA Z2","p"); 


  // TLegend *l4 = new TLegend(0.3628797,0.2534562,0.9777591,0.5745008,NULL,"BRNDC");
  // l4->SetHeader("UnCorr.  Corr.");
  // l4->SetNColumns(2);
  // l4->SetColumnSeparation(0.1);
  // l4->SetBorderSize(0);
  // l4->SetTextFont(42);
  // l4->SetTextSize(0.06);
  // l4->SetLineColor(1);
  // l4->SetLineStyle(1);
  // l4->SetLineWidth(1);
  // l4->SetFillColor(10);
  // l4->SetFillStyle(1001);
  // hArM_r[0][0][0]->SetTitle("   ");
  // l4->AddEntry(hArM_r[0][0][0]," ","p"); 
  // l4->AddEntry(hArM  [0][0][0],cleg[0],"p"); 
  // hArM_r[1][0][0]->SetTitle("   ");
  // l4->AddEntry(hArM_r[1][0][0]," ","p"); 
  // l4->AddEntry(hArM[1][0][0],cleg[1],"p"); 
  // hArM_r[2][0][0]->SetTitle("  ");
  // l4->AddEntry(hArM_r[2][0][0]," ","p"); 
  // l4->AddEntry(hArM[2][0][0],cleg[2],"p"); 


  TLine *line = new TLine(ptbins[0],1.00,ptbins[b1],1.00);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->SetLineColor(1);

  int maxc=3;
  maxr=2;

  TCanvas *c3[ktype], *c4=0, *c5;
  //std::cout<<std::endl;
  ipad=0;
  int ik=-1;
  for(int nj=0; nj<knj;nj++){
    if(nj%3==0){
      ++ik;
      c3[ik] = new TCanvas(Form("c3_%d",ik),Form("%s JES JER",alg[ik]),1211,652);
      makeMultiPanelCanvas(c3[ik],maxc,maxr,0.0,0.0,0.22,0.22,0.02,0);
      ipad=0;	
      //cout<<" nj : " <<nj <<endl;
    }
    
    for(int ic=0;ic<ncen;ic++){
      
      hRMS[plot[nj]][ic]->SetMarkerStyle(isty);
      hRMS[plot[nj]][ic]->SetMarkerColor(icol);
      hRMS[plot[nj]][ic]->SetLineColor(icol);
      hRMS[plot[nj]][ic]->SetMarkerSize(1.3);
      
      hArM[plot[nj]][ic]->SetMarkerStyle(isty);
      hArM[plot[nj]][ic]->SetMarkerColor(icol);
      hArM[plot[nj]][ic]->SetLineColor(icol);
      hArM[plot[nj]][ic]->SetMarkerSize(1.3);
      
      hArM_r[plot[nj]][ic]->SetMarkerStyle(jsty+4);
      //hArM_r[plot[nj]][ic]->SetMarkerStyle(24);
      hArM_r[plot[nj]][ic]->SetMarkerColor(icol);
      hArM_r[plot[nj]][ic]->SetLineColor(icol);
      hArM_r[plot[nj]][ic]->SetMarkerSize(1.3);
      
      hSigma[plot[nj]][ic]->SetMarkerStyle(isty);
      hSigma[plot[nj]][ic]->SetMarkerColor(icol);
      hSigma[plot[nj]][ic]->SetLineColor(icol);
      hSigma[plot[nj]][ic]->SetMarkerSize(1.3);
      
      hMean[plot[nj]][ic]->SetMarkerStyle(isty);
      hMean[plot[nj]][ic]->SetMarkerColor(icol);
      hMean[plot[nj]][ic]->SetLineColor(icol);
      hMean[plot[nj]][ic]->SetMarkerSize(1.3);
      
      hMean_r[plot[nj]][ic]->SetMarkerStyle(jsty+4);
      hMean_r[plot[nj]][ic]->SetMarkerColor(icol);
      hMean_r[plot[nj]][ic]->SetLineColor(icol);
      hMean_r[plot[nj]][ic]->SetMarkerSize(1.3);
      
      
      c3[ik]->cd(++ipad);
      gPad->SetLogx();
      
      hRMS[plot[nj]][ic]->SetMaximum(0.532);
      hSigma[plot[nj]][ic]->SetMaximum(0.532);
      if(iSigma){
	hSigma[plot[nj]][ic]->Draw("p");
      }else{
	hRMS[plot[nj]][ic]->Draw("p");
      }
      
      if(ipad!=1){
	MakeZero(hRMS[plot[nj]][ic]);
	MakeZero(hSigma[plot[nj]][ic]);
	MakeZero(hMean[plot[nj]][ic]);
	MakeZero(hArM[plot[nj]][ic]);
      }
      
      if(ipad==1){
	drawText2("CMS Simulation", 0.28, 0.86, 22);
	//drawText2("HYDJET 1.8",0.28,0.75,21);
	drawText2("PYTHIA Z2 (prod22 v81)",0.28,0.75,21);
      }
      else if(ipad==2){
	drawText2("#sqrt{s_{NN}} = 2.76 TeV",0.45,0.85,21);
	drawText2(Form("|#eta| < %0.0f",ketacut),0.15,0.85,21);
	//l4->Draw();
      }
      
      //drawText2("0-100%",0.75,0.82,21);
      drawText2(Form("%s",calgo[plot[nj]]),0.28,0.65,21);      
      
      c3[ik]->cd(ipad+3);
      
      //cout<<"pad : "<<(ipad+(ncen-1))<<"\t "<<calgo[nj]<<"\t centrality : "<<ccent[ic]<<"\t name : "<<hArM[nj][ic]->GetName()<<endl;
      gPad->SetLogx();
      
      hArM[plot[nj]][ic]->SetMaximum(1.05);      
      hArM[plot[nj]][ic]->SetMinimum(0.95);      
      hMean[plot[nj]][ic]->SetMaximum(1.05);      
      hMean[plot[nj]][ic]->SetMinimum(0.95);      
      
      //! Low pT
      // hArM[plot[nj]][ic]->SetMaximum(1.184);      
      // hArM[plot[nj]][ic]->SetMinimum(0.95);      
      // hMean[plot[nj]][ic]->SetMaximum(1.18);      
      // hMean[plot[nj]][ic]->SetMinimum(0.95);      
      
      //! With Raw/gen 
      // hArM[plot[nj]][ic]->SetMaximum(1.05);      
      // hArM[plot[nj]][ic]->SetMinimum(0.15);      
      // hMean[plot[nj]][ic]->SetMaximum(1.05);      
      // hMean[plot[nj]][ic]->SetMinimum(0.15);      
      
      if(iSigma){
	hMean[plot[nj]][ic]->Draw("p");
	//hMean_r[plot[nj]][ic]->Draw("psame");
      }else{
	hArM[plot[nj]][ic]->Draw("p");
	//hArM_r[plot[nj]][ic]->Draw("psame");
      }
      if(ipad==2)l4->Draw();
      line->Draw();
    }//! icen
    c3[ik]->SaveAs(Form("JESJER_pT_full_%s_prod22v81_pp2013.pdf",alg[ik]));
    c3[ik]->SaveAs(Form("JESJER_pT_full_%s_prod22v81_pp2013.gif",alg[ik]));
  }//! knj

  // TFile *fout1 = new TFile("ppJESJER.root","RECREATE");
  // fout1->cd();
  // for(int nj=0;nj<12;nj++){
  //   hArM[nj][0]->SetName(Form("hArm_%s",calgo[nj]));
  //   hArM[nj][0]->Write();
  //   hRMS[nj][0]->SetName(Form("hRMS_%s",calgo[nj]));
  //   hRMS[nj][0]->Write();
  // }
  // fout1->Close();

  return 0;

  c5 = new TCanvas("c5","Genj-Recoj",1302,842);
  c5->Divide(3,2);
  ipad=0;
  for(int nj=0; nj<knj;nj++){
    for(int ic=0;ic<ncen;ic++){
      
      if(nj==0)ipad=1;
      else if(nj==1)ipad=4;
      
      if(nj==0)ipad=2;
      else if(nj==1)ipad=5;
      
      if(nj==0)ipad=3;
      else if(nj==1)ipad=6;
      
      c5->cd(ipad);
      
      gPad->SetBottomMargin(0.2);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.05);
      
      hgenjrecoj[nj][ic]->SetTitle("");
      hgenjrecoj[nj][ic]->GetXaxis()->SetRangeUser(0,180);
      hgenjrecoj[nj][ic]->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
      hgenjrecoj[nj][ic]->GetXaxis()->CenterTitle(true);
      hgenjrecoj[nj][ic]->GetXaxis()->SetMoreLogLabels();
      hgenjrecoj[nj][ic]->GetXaxis()->SetNoExponent();
      hgenjrecoj[nj][ic]->GetXaxis()->SetNdivisions(507);
      hgenjrecoj[nj][ic]->GetXaxis()->SetLabelFont(42);
      hgenjrecoj[nj][ic]->GetXaxis()->SetLabelOffset(0.01);
      hgenjrecoj[nj][ic]->GetXaxis()->SetLabelSize(0.06);
      hgenjrecoj[nj][ic]->GetXaxis()->SetTitleSize(0.06);
      hgenjrecoj[nj][ic]->GetXaxis()->SetTitleOffset(1.15);
      hgenjrecoj[nj][ic]->GetXaxis()->SetTitleFont(42);
      
      hgenjrecoj[nj][ic]->GetYaxis()->SetTitle("RecoJet p_{T} (GeV/c)");
      hgenjrecoj[nj][ic]->GetYaxis()->SetRangeUser(0,180);
      hgenjrecoj[nj][ic]->GetYaxis()->CenterTitle(true);
      hgenjrecoj[nj][ic]->GetYaxis()->SetNdivisions(507);
      hgenjrecoj[nj][ic]->GetYaxis()->SetLabelFont(42);
      hgenjrecoj[nj][ic]->GetYaxis()->SetLabelOffset(0.01);
      hgenjrecoj[nj][ic]->GetYaxis()->SetLabelSize(0.06);
      hgenjrecoj[nj][ic]->GetYaxis()->SetTitleSize(0.06);
      hgenjrecoj[nj][ic]->GetYaxis()->SetTitleOffset(1.12);
      hgenjrecoj[nj][ic]->GetYaxis()->SetTitleFont(42);
      hgenjrecoj[nj][ic]->GetYaxis()->SetDecimals(true);
      hgenjrecoj[nj][ic]->Draw("colz");
      
      if(ipad==1){
	drawText2("HYDJET 1.8 MB",0.28,0.83,21);
	//drawText2("PYTHIA + HYDJET 1.8",0.28,0.85,21);
      }
      else if(ipad==2){
	drawText2(Form("|#eta| < %0.0f",ketacut),0.75,0.83,21);
      }
      drawText2(Form("%s",cleg),0.25,0.92,21);
      drawText2("0-100%",0.75,0.92,21);
      drawText2(Form("%s",calgo[nj]),0.28,0.78,21);      
    }
  }
  return 0;

  c4 = new TCanvas("c4","pt",100,142,1528,658);
  c4->Divide(2,2,0,0);
  //std::cout<<std::endl;
  ipad=0;
  for(int nj=0; nj<knj;nj++){
    for(int ic=0;ic<ncen;ic++){
      c4->cd(++ipad);
      gPad->SetLogy();
      gPad->SetBottomMargin(0.2);
      if(ipad==1)gPad->SetLeftMargin(0.15);
      else if(ipad==4)gPad->SetRightMargin(0.05);
      
      hgenpt[nj][ic]->SetMaximum(1e+04);
      hgenpt[nj][ic]->Draw("p");
      hrecpt[nj][ic]->Draw("psame");
      hrawpt[nj][ic]->Draw("psame");
      
      if(ipad==1){
	drawText2("CMS Preliminary", 0.28, 0.90, 22);
	drawText2("HYDJET 1.8 MB",0.28,0.85,21);
	//drawText2("PYTHIA + HYDJET 1.8",0.28,0.85,21);
      }
      else if(ipad==2){
	drawText2(Form("|#eta| < %0.0f",ketacut),0.75,0.85,21);
      }else if(ipad==3)l3->Draw();
      drawText2("0-100%",0.75,0.90,21);
      drawText2(Form("%s",calgo[nj]),0.58,0.75,21);      
    }//! icen
  }//! knj
  return 0;

}
void DefineHisto(const char *cver, int nj, int ic, TH1F *&h1, const char *name, const char *caption)
{
  h1 = new TH1F(Form("%s_%s_%d_%d",name,cver,nj,ic),Form("CMSSW %s %s %s %s %d",cver,caption,ksp[ic],calgo[nj],ic),nbins,ptbins);
  MakeHist(h1,statop,"p_{T}^{GenJet} (GeV/c)","<p_{T}^{RecoJet}/p_{T}^{GenJet}>");
}

void CalJesJer(const char *cver, int nj, int ic, TH2F *h2d, TH1F **h1d, TH1F *ham, TH1F *hr, TH1F *hm, TH1F *hs)
{

  // if(nj==0){
  //   fout = new TFile("output.root","RECREATE");
  // }
  for(int ip=0;ip<nbins;ip++){
    int lbin = h2d->GetXaxis()->FindBin(ptbins[ip]);
    int hbin = h2d->GetXaxis()->FindBin(ptbins[ip+1]);
    h1d[ip]  = (TH1F*)h2d->ProjectionY(Form("hratiocorrrefpt1D_%d_%d_%d",nj,ic,ip),lbin,hbin,"e");
    if(h1d[ip]->GetEntries()<maxEntry)continue;
    //h1d[ip]->Rebin(5);
    // if(nj==0){
    //   fout->cd();
    //   h1d[ip]->Write();
    // }
    FillMeanSigma(ip,h1d[ip],ham,hr,hm,hs);
  }
}
// void FillMeanSigma(int ip,TH1F *&h1F, TH1F *harm,TH1F *hrms,TH1F *hm,TH1F *hs)
// {
//   if(h1F->GetEntries() < maxEntry || 0==h1F){
//     h1F->Scale(0.);
//     harm  ->SetBinContent(ip+1,-9999);
//     harm  ->SetBinError  (ip+1,0);
//     hrms  ->SetBinContent(ip+1,-9999);
//     hrms  ->SetBinError  (ip+1,0);
//     std::cout<<" WARNING :::::  " << " Not enough entries, please change the binning :: "   << std::endl;
//     return;
//   }else{
//     h1F->Scale(1./h1F->Integral());
//   }
  
//   //! Fill arthematic mean and RMS
//   harm->SetBinContent(ip+1,h1F->GetMean());
//   harm->SetBinError  (ip+1,h1F->GetMeanError());
//   hrms->SetBinContent(ip+1,h1F->GetRMS());
//   hrms->SetBinError  (ip+1,h1F->GetRMSError());
  
//   double nsigma=2.0;
//   double niter=10;
//   if(iFit==1){//! Gaussian Fit
//     fit_gaussian(h1F, nsigma, niter);
//   }
//   else if(iFit==2){//! Crystal Ball
//     int fitstatus = fit_dscb(h1F, nsigma, niter);
//   }
//   TF1 *f1 = (TF1*)h1F->GetListOfFunctions()->Last();
//   if(0==f1){
//     delete f1; return;
//   }else{
//     hm->SetBinContent(ip+1,f1->GetParameter(1));
//     hs->SetBinContent(ip+1,f1->GetParameter(2));      
//     hm->SetBinError  (ip+1,f1->GetParError(1));
//     hs->SetBinError  (ip+1,f1->GetParError(2));
//     delete f1;
//   }
// }

void fit_gaussian(TH1F *&hrsp, const double nsigma, const int niter)
{
  if (0==hrsp) {
    cout<<"ERROR: Empty pointer to fit_gaussian()"<<endl;return;
  }
  
  double mean     = hrsp->GetMean();
  double rms      = hrsp->GetRMS();
  double rspMax   = 1.0;

  double norm  = hrsp->GetMaximumStored();
  double peak  = mean;
  double sigma = rms;

  double xmin_tmp  = hrsp->GetXaxis()->GetXmin();
  double xmax_tmp  = hrsp->GetXaxis()->GetXmax();
  TF1* fitfnc(0); int fitstatus(-1);
  for (int iiter=0;iiter<niter;iiter++) {
    vector<double> vv;
    vv.push_back(rspMax);
    vv.push_back(xmin_tmp);
    vv.push_back(peak-nsigma*sigma);
    double fitrange_min = *std::max_element(vv.begin(),vv.end());
    double fitrange_max = std::min(xmax_tmp,peak+nsigma*sigma);
    fitrange_min = 0.2;                                 
    fitrange_max = 2.0;
    adjust_fitrange(hrsp,fitrange_min,fitrange_max);
    fitfnc = new TF1("fgaus","gaus",fitrange_min,fitrange_max);
    fitfnc->SetParNames("N","#mu","#sigma");
    fitfnc->SetParameter(0,norm);
    fitfnc->SetParameter(1,peak);
    fitfnc->SetParameter(2,sigma);
    fitstatus = hrsp->Fit(fitfnc,"RQ+");
    delete fitfnc;
    fitfnc = hrsp->GetFunction("fgaus");
    //fitfnc->ResetBit(TF1::kNotDraw);
    norm  = fitfnc->GetParameter(0);
    peak  = fitfnc->GetParameter(1);
    sigma = fitfnc->GetParameter(2);
  }
  if(hrsp->GetFunction("fgaus")==0){
    cout << "No function recorded in histogram " << hrsp->GetName() << endl;
  }
  if (0!=fitstatus){
    cout<< " fit_gaussian() to " << hrsp->GetName()  
        <<" failed. Fitstatus: "<<fitstatus
        <<" - FNC deleted."<<endl;
    hrsp->GetListOfFunctions()->Delete();
  }
}
int fit_dscb(TH1F *&hrsp, const double nsigma, const int niter)
{
  if (0==hrsp) {
    cout<<"ERROR: Empty pointer to fit_dscb()"<<endl;return -1;
  }

  // first use a gaussian to constrain crystal ball gaussian core  
  fit_gaussian(hrsp,nsigma,niter);

  TF1* fgaus = hrsp->GetFunction("fgaus");
  if (0==fgaus) {
    hrsp->GetListOfFunctions()->Delete();
    return -1;
  }

  // implementation of the low pt bias threshold
  double fitrange_min = 0.2;
  double fitrange_max = 2.0;
  
  adjust_fitrange(hrsp,fitrange_min,fitrange_max);
  TF1* fdscb = new TF1("fdscb",fnc_dscb,fitrange_min,fitrange_max,7);
  
  // set the std values
  double norm = fgaus->GetParameter(0);
  double mean = fgaus->GetParameter(1);
  double sigma= fgaus->GetParameter(2);
  double aone=2.0, atwo=2.0, pone=5.0, ptwo=5.0;
  //double aone=5.0, atwo=5.0, pone=15.0, ptwo=15.0;
  TVirtualFitter::SetDefaultFitter("Minuit2");

  int fitstatus(0);
  for (unsigned i=0;i<4;i++) {
    
    fdscb->SetParameter(0,norm); // N 
    fdscb->SetParameter(1,mean); // mean
    fdscb->SetParameter(2,sigma);// sigma
    fdscb->SetParameter(3,aone); // a1
    fdscb->SetParameter(4,pone); // p1
    fdscb->SetParameter(5,atwo); // a2
    fdscb->SetParameter(6,ptwo); // p2
    
    fdscb->FixParameter(1,mean);
    fdscb->FixParameter(2,sigma);

    if (i>0) fdscb->FixParameter(3,aone);
    else fdscb->SetParLimits(3,1.,5.);

    if (i>1) fdscb->FixParameter(5,atwo);
    else fdscb->SetParLimits(5,1.,5.);

    fdscb->SetParLimits(4,0.,25.);
    fdscb->SetParLimits(6,0.,25.);

    fitstatus = hrsp->Fit(fdscb,"RQB+");
    if (0==fitstatus) i=999;
    
    delete fdscb;
    fdscb = hrsp->GetFunction("fdscb");
    
    if (0==fdscb) return -1;

    norm  = fdscb->GetParameter(0);
    aone  = fdscb->GetParameter(3);
    pone  = fdscb->GetParameter(4);
    atwo  = fdscb->GetParameter(5);
    ptwo  = fdscb->GetParameter(6);

  }

  // reset sigma and mean to gauss values...
  fdscb->SetParameter(1,fgaus->GetParameter(1));
  fdscb->SetParError (1,fgaus->GetParError(1));
  fdscb->SetParameter(2,fgaus->GetParameter(2));
  fdscb->SetParError (2,fgaus->GetParError(2));
  if (0!=fitstatus){
    cout<<"fit_fdscb() to "<<hrsp->GetName()
        <<" failed. Fitstatus: "<<fitstatus<<endl;
    hrsp->GetFunction("fdscb")->Delete();
  }
  
  return fitstatus;
}


double fnc_dscb(double *xx,double *pp)
{
  double x   = xx[0];
  // gaussian core 
  double N   = pp[0];//norm
  double mu  = pp[1];//mean
  double sig = pp[2];//variance
  // transition parameters
  double a1  = pp[3];
  double p1  = pp[4];
  double a2  = pp[5];
  double p2  = pp[6];

  double u   = (x-mu)/sig;
  double A1  = pow(p1/fabs(a1),p1)*exp(-a1*a1/2);
  double A2  = pow(p2/fabs(a2),p2)*exp(-a2*a2/2);
  double B1  = p1/fabs(a1) - fabs(a1);
  double B2  = p2/fabs(a2) - fabs(a2);

  double result(N);
  if      (u<-a1) result *= A1*pow(B1-u,-p1);
  else if (u<a2)  result *= exp(-u*u/2);
  else            result *= A2*pow(B2+u,-p2);
  return result;
}
void adjust_fitrange(TH1F *h,double& min,double& max)
{
  int imin=1; while (h->GetBinLowEdge(imin)<min) imin++;
  int imax=1; while (h->GetBinLowEdge(imax)<max) imax++;
  while ((imax-imin)<8) {
    if (imin>1) {imin--; min = h->GetBinCenter(imin); }
    if (imax<h->GetNbinsX()-1) { imax++; max=h->GetBinCenter(imax); }
  }
}

 void FillMeanSigma(int ip,TH1F *&h1F,TH1F *harm,TH1F *hrms,TH1F *hm,TH1F *hs)
 {
   TF1 *f1 = new TF1("f1","(([0]/(2*pi*[1]*[1]))*exp(-1.*((x-[2])*(x-[2])/(2*[1]*[1]))))",fitmin,fitmax);
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
       h1F->Fit("gaus",fopt,"",fitmin,fitmax);
       TF1* f2 = (TF1*)h1F->GetFunction("gaus");
       f2->SetLineWidth(1);
       f2->SetLineStyle(2);
       f2->SetNpx(knpx);
       hm->SetBinContent(ip+1,f2->GetParameter(1));
       hs->SetBinContent(ip+1,f2->GetParameter(2));
    
       if(strcmp(fopt,"MLLRQ0+")==0){
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
    
       h1F->Fit("f1",fopt,"",fitmin,fitmax);
       //h1F->Fit("f1","MLLRQ0+","",fitmin,fitmax);
       hm->SetBinContent(ip+1,f1->GetParameter(2));
       hs->SetBinContent(ip+1,f1->GetParameter(1));
    
       if(strcmp(fopt,"MLLRQ+")==0){
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

void MakeHistRMS(TH1 *h1,float ymax,float ymin)
{

  h1->SetTitle("");
  h1->SetMaximum(ymax);
  h1->SetMinimum(ymin);
  h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
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
  h1->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
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
