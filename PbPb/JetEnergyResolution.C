
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>

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


using namespace std;

const double pi=acos(-1.);
const double pi2=2*pi -1;

int bins =500;
float minval = 0;
float maxval = 1000;
float binw   = 2;

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
void FillMeanSigma(int /*ip*/,TH1 */*h1*/,TH1 */*ArM*/,TH1 */*RMS*/,TH1 */*Mean*/,TH1 */*Sigma*/);
void CalJesJer(const char */*cver*/, int /*nj*/, int /*icen*/, TH2 */*h2d*/, TH1F **/*h1d*/, TH1 */*ham*/, TH1 */*hr*/, TH1 */*hm*/, TH1 */*hs*/);

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


//! used for fit
//double ptbins[]   = {80,90,100,110,120,130,140,150,160,170,180,190,200,286};
//double ptbins[]     = {15,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,240,280,300,350,400,450,500,645};
double ptbins[]     = {14,16,18,20,22,24,26,28,30,32,34,36,38,40};
const int b1  = sizeof(ptbins)/sizeof(Double_t) - 1;
const int nbins = b1;

// const int kfiles=3;
// const char *version[kfiles] = {"44x-80","53x-Track7-80_v2","53x-LV1-Track7-80"};
// const char *cleg[kfiles]    = {"44x","53x_V28-Track7","53x_LV1-Track7"};
// const int icol[kfiles] ={1,2,4};
// const int isty[kfiles] ={24,25,30};
// const int jsty[kfiles] ={20,21,29};


const int kfiles=1;
const char *version[kfiles] = {"53x-80-Track8-Jet22"};
const char *cleg[kfiles]    = {"Track8-Jet22"};
const int icol[kfiles] ={1};
const int isty[kfiles] ={20};
const int jsty[kfiles] ={20};

TFile *fout=0;

int JetEnergyResolution(const char *reta = "eta2")
{
  
  LoadStyle();
  
  int rfit=4;

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
  TH1F *hgenpt[kfiles][knj][ncen], *hrecpt[kfiles][knj][ncen], *hrawpt[kfiles][knj][ncen];
  
  
  TH2F *hratiocorrrefpt[kfiles][knj][ncen];
  TH1F *hratiocorrrefpt1D[kfiles][knj][ncen][nbins];  
  TH1F *hMean[kfiles][knj][ncen], *hSigma[kfiles][knj][ncen], *hRMS[kfiles][knj][ncen], *hArM[kfiles][knj][ncen];

  //! Response raw pT / gen pT  vs gen pT
  TH2F *hratiorawrefpt[kfiles][knj][ncen];
  TH1F *hratiorawrefpt1D[kfiles][knj][ncen][nbins];  
  TH1F *hMean_r[kfiles][knj][ncen], *hSigma_r[kfiles][knj][ncen], *hRMS_r[kfiles][knj][ncen], *hArM_r[kfiles][knj][ncen];

  //! Ratio of Fit and geometric mean/RMS
  //TH1F *hRatio_Mean[knj][ncen], *hRatio_RMS[knj][ncen];

  TH2F *hgenjrecoj[kfiles][knj][ncen];

  for(int in=0; in<kfiles; in++){
    //fin_pbpb[in] = new TFile(Form("input/testoutput_%s.root",version[in]),"r");
    fin_pbpb[in] = new TFile("input/JetResponse_histos_pp.root","r");
    //fin_pbpb[in] = new TFile("input/JetResponse_histos_vsalgo_pp.root","r");
    //fin_pbpb[in] = new TFile("JetResponse_histos_pp_embedded.root","r");
    for(int nj=0;nj<knj;nj++){
      //cout<<"nj : "<<nj<<Form("\t %s",calgo[nj])<<endl;
      for(int icen=0;icen<ncen;icen++){
	
	DefineHisto(cleg[in], nj, icen, hMean   [in][nj][icen],"hMean"    ,"<reco p_{T}/gen p_{T}>");
	DefineHisto(cleg[in], nj, icen, hMean_r [in][nj][icen],"hMean_r"  ,"<reco p_{T}/gen p_{T}>");
	DefineHisto(cleg[in], nj, icen, hArM    [in][nj][icen],"hArm"     ,"<reco p_{T}/gen p_{T}>");
	DefineHisto(cleg[in], nj, icen, hArM_r  [in][nj][icen],"hArm_r"   ,"<reco p_{T}/gen p_{T}>");
	DefineHisto(cleg[in], nj, icen, hSigma  [in][nj][icen],"hSigma"   ,"#sigma(reco p_{T}/gen p_{T})");
	DefineHisto(cleg[in], nj, icen, hSigma_r[in][nj][icen],"hSigma_r" ,"#sigma(reco p_{T}/gen p_{T})");
	DefineHisto(cleg[in], nj, icen, hRMS    [in][nj][icen],"hRMS"     ,"#sigma(reco p_{T}/gen p_{T})");
	DefineHisto(cleg[in], nj, icen, hRMS_r  [in][nj][icen],"hRMS_r"   ,"#sigma(reco p_{T}/gen p_{T})");
	
	//cout<<"icen : "<<ccent[icen]<<endl;
	//! PbPb  starts here
	hgenjrecoj[in][nj][icen] = (TH2F*)fin_pbpb[in]->Get(Form("hgenjrecoj%d_%d",nj,icen));

	hratiocorrrefpt[in][nj][icen] = (TH2F*)fin_pbpb[in]->Get(Form("hrescrpt_genm%d_%d",nj,icen));
	hratiocorrrefpt[in][nj][icen]->SetName(Form("hratiocorrrefpt_pbpb_%s_%d_%d",version[in],nj,icen));
	hratiorawrefpt [in][nj][icen] = (TH2F*)fin_pbpb[in]->Get(Form("hresrrpt_genm%d_%d",nj,icen));
	hratiorawrefpt [in][nj][icen]->SetName(Form("hratiorawrefpt_pbpb_%s_%d_%d",version[in],nj,icen));
      
	//if(nj==1)cout<<nj<<"\t icen : "<<icen<<"\t hMean : "<<hMean[nj][icen]->GetName()<<endl;
	//! Calculate jes and jer
	//! Raw / Gen
	CalJesJer(version[in], nj, icen, 
		  hratiorawrefpt [in][nj][icen], hratiorawrefpt1D [in][nj][icen], 
		  hArM_r[in][nj][icen], hRMS_r[in][nj][icen], 
		  hMean_r[in][nj][icen],hSigma_r[in][nj][icen]);
      
	//! Reco / Gen
	CalJesJer(version[in], nj, icen, 
		  hratiocorrrefpt[in][nj][icen], hratiocorrrefpt1D[in][nj][icen], 
		  hArM[in][nj][icen], hRMS[in][nj][icen], 
		  hMean[in][nj][icen],hSigma[in][nj][icen]);


	hgenpt[in][nj][icen] = (TH1F*)fin_pbpb[in]->Get(Form("hgenpt_genm%d_%d",nj,icen));
	hgenpt[in][nj][icen]->SetMarkerStyle(24);
	hgenpt[in][nj][icen]->SetMarkerSize(1.2);
	hgenpt[in][nj][icen]->SetMarkerColor(1);
	hgenpt[in][nj][icen]->Rebin(5);
	hgenpt[in][nj][icen]->SetTitle(Form("%s",calgo[nj]));
	hgenpt[in][nj][icen]->GetXaxis()->SetRangeUser(ptbins[0],ptbins[b1]);
	hgenpt[in][nj][icen]->GetXaxis()->SetTitle("jet p_{T} (GeV/c)");
	hgenpt[in][nj][icen]->GetXaxis()->CenterTitle();
	hgenpt[in][nj][icen]->GetXaxis()->SetTitleSize(0.07);
	hgenpt[in][nj][icen]->GetXaxis()->SetLabelSize(0.07);
	hgenpt[in][nj][icen]->GetYaxis()->SetRangeUser(0.1,6e+02);
	hgenpt[in][nj][icen]->GetYaxis()->SetTitle("Counts");
	hgenpt[in][nj][icen]->GetYaxis()->CenterTitle();
	hgenpt[in][nj][icen]->GetYaxis()->SetTitleSize(0.07);
	hgenpt[in][nj][icen]->GetYaxis()->SetLabelSize(0.07);

	hrecpt[in][nj][icen] = (TH1F*)fin_pbpb[in]->Get(Form("hrecopt_genm%d_%d",nj,icen));
	hrecpt[in][nj][icen]->SetMarkerStyle(20);
	hrecpt[in][nj][icen]->SetMarkerSize(1.2);
	hrecpt[in][nj][icen]->SetMarkerColor(2);
	hrecpt[in][nj][icen]->Rebin(5);
	//hrecpt[in][nj][icen]->GetXaxis()->SetRangeUser(ptbins[0],ptbins[b1]);

	hrawpt[in][nj][icen] = (TH1F*)fin_pbpb[in]->Get(Form("hrawpt_genm%d_%d",nj,icen));
	hrawpt[in][nj][icen]->SetMarkerStyle(25);
	hrawpt[in][nj][icen]->SetMarkerSize(1.2);
	hrawpt[in][nj][icen]->SetMarkerColor(4);
	hrawpt[in][nj][icen]->Rebin(5);
	//hrawpt[in][nj][icen]->GetXaxis()->SetRangeUser(ptbins[0],ptbins[b1]);
      
	//! PbPb ends
      }//! icen loop ends
    }//! nj loop ends
  }//! version 
  //return 0;

  for(int in=0; in<kfiles; in++){
    for(int nj=0;nj<knj;nj++){
      for(int ic=0;ic<ncen;ic++){
	//MakeHistMean(hArM  [in][nj][ic],1.19,0.84);
	//MakeHistMean(hMean [in][nj][ic],1.19,0.84);
	MakeHistMean(hArM  [in][nj][ic],1.05,0.95);
	MakeHistMean(hMean [in][nj][ic],1.05,0.95);

	MakeHistRMS (hRMS  [in][nj][ic],0.43,0.001);
	MakeHistRMS (hSigma[in][nj][ic],0.43,0.001);
      }      
    }
  }

  TCanvas *c3[4], *c4=0, *c5;
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
  l3->AddEntry(hgenpt[0][0][0],"Gen","p"); 
  l3->AddEntry(hrecpt[0][0][0],"Reco","p"); 
  l3->AddEntry(hrawpt[0][0][0],"Raw","p"); 


  //TLegend *l4 = new TLegend(0.2527884,0.6276059,0.5060804,0.963133,NULL,"BRNDC");
  TLegend *l4 = new TLegend(0.1044898,0.7233131,0.7247845,0.9396888,NULL,"BRNDC");
  l4->SetHeader("CMSSW 5_3_16");
  l4->SetBorderSize(0);
  l4->SetTextFont(42);
  l4->SetTextSize(0.07);
  l4->SetLineColor(1);
  l4->SetLineStyle(1);
  l4->SetLineWidth(1);
  l4->SetFillColor(10);
  l4->SetFillStyle(1001);
  l4->AddEntry(hRMS[0][0][0],cleg[0],"p"); 
  //l4->AddEntry(hRMS[1][0][0],cleg[1],"p"); 
  //l4->AddEntry(hRMS[2][0][0],cleg[2],"p"); 


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

  const char *alg[4] = {"VsPF","VsCalo","PF","Calo"};
  int plot[knj]={0,1,2,3,4,5,6,7,8,9,10,11};
  //std::cout<<std::endl;
  ipad=0;
  for(int in=0; in<kfiles; in++){
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
	
	hRMS[in][plot[nj]][ic]->SetMarkerStyle(isty[in]);
	hRMS[in][plot[nj]][ic]->SetMarkerColor(icol[in]);
	hRMS[in][plot[nj]][ic]->SetLineColor(icol[in]);
	hRMS[in][plot[nj]][ic]->SetMarkerSize(1.3);

	hArM[in][plot[nj]][ic]->SetMarkerStyle(isty[in]);
	hArM[in][plot[nj]][ic]->SetMarkerColor(icol[in]);
	hArM[in][plot[nj]][ic]->SetLineColor(icol[in]);
	hArM[in][plot[nj]][ic]->SetMarkerSize(1.3);

	//hArM_r[in][plot[nj]][ic]->SetMarkerStyle(jsty[in]);
	hArM_r[in][plot[nj]][ic]->SetMarkerStyle(24);
	hArM_r[in][plot[nj]][ic]->SetMarkerColor(icol[in]);
	hArM_r[in][plot[nj]][ic]->SetLineColor(icol[in]);
	hArM_r[in][plot[nj]][ic]->SetMarkerSize(1.3);

	hSigma[in][plot[nj]][ic]->SetMarkerStyle(isty[in]);
	hSigma[in][plot[nj]][ic]->SetMarkerColor(icol[in]);
	hSigma[in][plot[nj]][ic]->SetLineColor(icol[in]);
	hSigma[in][plot[nj]][ic]->SetMarkerSize(1.3);

	hMean[in][plot[nj]][ic]->SetMarkerStyle(isty[in]);
	hMean[in][plot[nj]][ic]->SetMarkerColor(icol[in]);
	hMean[in][plot[nj]][ic]->SetLineColor(icol[in]);
	hMean[in][plot[nj]][ic]->SetMarkerSize(1.3);

	hMean_r[in][plot[nj]][ic]->SetMarkerStyle(jsty[in]);
	hMean_r[in][plot[nj]][ic]->SetMarkerColor(icol[in]);
	hMean_r[in][plot[nj]][ic]->SetLineColor(icol[in]);
	hMean_r[in][plot[nj]][ic]->SetMarkerSize(1.3);


	c3[ik]->cd(++ipad);
	//gPad->SetLogx();
	
	hRMS[in][plot[nj]][ic]->SetMaximum(0.532);
	hSigma[in][plot[nj]][ic]->SetMaximum(0.532);
	if(iSigma){
	  if(in==0)hSigma[in][plot[nj]][ic]->Draw("p");
	  else hSigma[in][plot[nj]][ic]->Draw("psame");
	}else{
	  if(in==0)hRMS[in][plot[nj]][ic]->Draw("p");
	  else hRMS[in][plot[nj]][ic]->Draw("psame");
	}
	
	if(ipad!=1){
	  MakeZero(hRMS[in][plot[nj]][ic]);
	  MakeZero(hSigma[in][plot[nj]][ic]);
	  MakeZero(hMean[in][plot[nj]][ic]);
	  MakeZero(hArM[in][plot[nj]][ic]);
	}

	if(ipad==1){
	  drawText2("CMS Simulation", 0.28, 0.86, 22);
          //drawText2("HYDJET 1.8",0.28,0.75,21);
	  drawText2("PYTHIA Z2",0.28,0.75,21);
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
	//gPad->SetLogx();
 
	// hArM[in][plot[nj]][ic]->SetMaximum(1.05);      
	// hArM[in][plot[nj]][ic]->SetMinimum(0.95);      
	// hMean[in][plot[nj]][ic]->SetMaximum(1.05);      
	// hMean[in][plot[nj]][ic]->SetMinimum(0.95);      

	//! Low pT
	hArM[in][plot[nj]][ic]->SetMaximum(1.384);      
	hArM[in][plot[nj]][ic]->SetMinimum(0.95);      
	hMean[in][plot[nj]][ic]->SetMaximum(1.15);      
	hMean[in][plot[nj]][ic]->SetMinimum(0.95);      


	// hArM[in][plot[nj]][ic]->SetMaximum(1.05);      
	// hArM[in][plot[nj]][ic]->SetMinimum(0.15);      
	// hMean[in][plot[nj]][ic]->SetMaximum(1.05);      
	// hMean[in][plot[nj]][ic]->SetMinimum(0.15);      
	
	if(iSigma){
	  if(in==0){
	    hMean[in][plot[nj]][ic]->Draw("p");
	    //hMean_r[in][plot[nj]][ic]->Draw("psame");
	  }else{
	    hMean[in][plot[nj]][ic]->Draw("psame");
	    //hMean_r[in][plot[nj]][ic]->Draw("psame");
	  }
	}else{
	  if(in==0){
	    hArM[in][plot[nj]][ic]->Draw("p");
	    //hArM_r[in][plot[nj]][ic]->Draw("psame");
	  }else{
	    hArM[in][plot[nj]][ic]->Draw("psame");
	    //hArM_r[in][plot[nj]][ic]->Draw("psame");
	  }
	}
	if(ipad==2)l4->Draw();
	line->Draw();
      }//! icen
      //c3[ik]->SaveAs(Form("JESJER_pT_full_%s.gif",alg[ik]));
    }//! knj
  }//! in
  return 0;

  c5 = new TCanvas("c5","Genj-Recoj",1302,842);
  c5->Divide(3,2);
  ipad=0;
  for(int in=0; in<kfiles; in++){
    for(int nj=0; nj<knj;nj++){
      for(int ic=0;ic<ncen;ic++){
	
	if(in==0 && nj==0)ipad=1;
	else if(in==0 && nj==1)ipad=4;

	if(in==1 && nj==0)ipad=2;
	else if(in==1 && nj==1)ipad=5;

	if(in==2 && nj==0)ipad=3;
	else if(in==2 && nj==1)ipad=6;

	c5->cd(ipad);
	
	gPad->SetBottomMargin(0.2);
	gPad->SetLeftMargin(0.15);
	gPad->SetRightMargin(0.05);
	
	hgenjrecoj[in][nj][ic]->SetTitle("");
	hgenjrecoj[in][nj][ic]->GetXaxis()->SetRangeUser(0,180);
	hgenjrecoj[in][nj][ic]->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
	hgenjrecoj[in][nj][ic]->GetXaxis()->CenterTitle(true);
	hgenjrecoj[in][nj][ic]->GetXaxis()->SetMoreLogLabels();
	hgenjrecoj[in][nj][ic]->GetXaxis()->SetNoExponent();
	hgenjrecoj[in][nj][ic]->GetXaxis()->SetNdivisions(507);
	hgenjrecoj[in][nj][ic]->GetXaxis()->SetLabelFont(42);
	hgenjrecoj[in][nj][ic]->GetXaxis()->SetLabelOffset(0.01);
	hgenjrecoj[in][nj][ic]->GetXaxis()->SetLabelSize(0.06);
	hgenjrecoj[in][nj][ic]->GetXaxis()->SetTitleSize(0.06);
	hgenjrecoj[in][nj][ic]->GetXaxis()->SetTitleOffset(1.15);
	hgenjrecoj[in][nj][ic]->GetXaxis()->SetTitleFont(42);

	hgenjrecoj[in][nj][ic]->GetYaxis()->SetTitle("RecoJet p_{T} (GeV/c)");
	hgenjrecoj[in][nj][ic]->GetYaxis()->SetRangeUser(0,180);
	hgenjrecoj[in][nj][ic]->GetYaxis()->CenterTitle(true);
	hgenjrecoj[in][nj][ic]->GetYaxis()->SetNdivisions(507);
	hgenjrecoj[in][nj][ic]->GetYaxis()->SetLabelFont(42);
	hgenjrecoj[in][nj][ic]->GetYaxis()->SetLabelOffset(0.01);
	hgenjrecoj[in][nj][ic]->GetYaxis()->SetLabelSize(0.06);
	hgenjrecoj[in][nj][ic]->GetYaxis()->SetTitleSize(0.06);
	hgenjrecoj[in][nj][ic]->GetYaxis()->SetTitleOffset(1.12);
	hgenjrecoj[in][nj][ic]->GetYaxis()->SetTitleFont(42);
	hgenjrecoj[in][nj][ic]->GetYaxis()->SetDecimals(true);
	hgenjrecoj[in][nj][ic]->Draw("colz");

	if(ipad==1){
	  drawText2("HYDJET 1.8 MB",0.28,0.83,21);
	  //drawText2("PYTHIA + HYDJET 1.8",0.28,0.85,21);
	}
	else if(ipad==2){
	  drawText2(Form("|#eta| < %0.0f",ketacut),0.75,0.83,21);
	}
	drawText2(Form("%s",cleg[in]),0.25,0.92,21);
	drawText2("0-100%",0.75,0.92,21);
	drawText2(Form("%s",calgo[nj]),0.28,0.78,21);      
      }
    }
  }
  return 0;

  c4 = new TCanvas("c4","pt",100,142,1528,658);
  c4->Divide(2,2,0,0);
  //std::cout<<std::endl;
  for(int in=0; in<kfiles; in++){
    ipad=0;
    for(int nj=0; nj<knj;nj++){
      for(int ic=0;ic<ncen;ic++){
	c4->cd(++ipad);
	gPad->SetLogy();
	gPad->SetBottomMargin(0.2);
	if(ipad==1)gPad->SetLeftMargin(0.15);
	else if(ipad==4)gPad->SetRightMargin(0.05);
	
	hgenpt[in][nj][ic]->SetMaximum(1e+04);
	hgenpt[in][nj][ic]->Draw("p");
	hrecpt[in][nj][ic]->Draw("psame");
	hrawpt[in][nj][ic]->Draw("psame");
	
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
  }//! in

  return 0;
    


  
  // ipad=0;
  // TCanvas *c99[knj][ncen];
  // maxr=3;
  // for(int nj=1;nj<2;nj++){
  //   for(int ic=0;ic<ncen-1;ic++){
  //     c99[nj][ic] = new TCanvas(Form("c99_%d_%d",nj,ic),Form("%s Fitting plots %s",calgo[nj],ccent[ic]),100,102,1399,942);
  //     //c99[nj][ic]->Divide(6,maxr,0,0);
  //     c99[nj][ic]->Divide(7,4,0,0);
  //     ipad=0;
      
  //     for(int ip=0;ip<nbins;ip++){      
  // 	c99[nj][ic]->cd(++ipad);
  // 	if(ipad%7==0)gPad->SetRightMargin(0.02);
  // 	if(ipad==1 || ipad==8 || ipad==15 || ipad==22)gPad->SetLeftMargin(0.15);
  // 	gPad->SetBottomMargin(0.15);
  // 	gPad->SetLogy();

  // 	hratiocorrrefpt1D[nj][ic][ip]->SetMaximum(25.634);
  // 	hratiocorrrefpt1D[nj][ic][ip]->SetMinimum(1e-09);
  // 	hratiocorrrefpt1D[nj][ic][ip]->SetTitle(0);
  // 	hratiocorrrefpt1D[nj][ic][ip]->GetXaxis()->SetTitle("<reco jet p_{T} / gen jet p_{T}>");
  // 	hratiocorrrefpt1D[nj][ic][ip]->GetXaxis()->SetTitleFont(42);
  // 	hratiocorrrefpt1D[nj][ic][ip]->GetXaxis()->SetLabelFont(42);
  // 	hratiocorrrefpt1D[nj][ic][ip]->GetXaxis()->SetLabelSize(0.08);
  // 	hratiocorrrefpt1D[nj][ic][ip]->GetXaxis()->SetTitleSize(0.07);
	
  // 	hratiocorrrefpt1D[nj][ic][ip]->GetYaxis()->SetTitle("");
  // 	hratiocorrrefpt1D[nj][ic][ip]->GetYaxis()->SetTitleFont(42);
  // 	hratiocorrrefpt1D[nj][ic][ip]->GetYaxis()->SetLabelFont(42);
  // 	hratiocorrrefpt1D[nj][ic][ip]->GetYaxis()->SetLabelSize(0.08);

  // 	hratiocorrrefpt1D[nj][ic][ip]->SetMarkerStyle(24);
  // 	hratiocorrrefpt1D[nj][ic][ip]->SetMarkerColor(1);
  // 	hratiocorrrefpt1D[nj][ic][ip]->SetLineColor(1);
  // 	hratiocorrrefpt1D[nj][ic][ip]->SetMarkerSize(1.1);
  // 	hratiocorrrefpt1D[nj][ic][ip]->Draw("p");  


  // 	// c99[nj][ic]->Update();
  // 	// TPaveStats *ps = (TPaveStats*)  hratiocorrrefpt1D[nj][ic][ip]->GetListOfFunctions()->FindObject("stats");
  // 	// ps->SetX1NDC(0.50);
  // 	// ps->SetY1NDC(0.41);       
  // 	// ps->SetX2NDC(0.95);
  // 	// ps->SetY2NDC(0.79);
  // 	// ps->SetTextFont(42);
  // 	// ps->Draw();

  // 	//pp
  // 	if(ic==0){
  // 	  hratiocorrrefpt1D[nj][ncen-1][ip]->SetMaximum(25.634);
  // 	  hratiocorrrefpt1D[nj][ncen-1][ip]->SetMinimum(1e-09);
  // 	  hratiocorrrefpt1D[nj][ncen-1][ip]->SetTitle(0);
  // 	  hratiocorrrefpt1D[nj][ncen-1][ip]->GetXaxis()->SetTitle("<reco jet p_{T} / gen jet p_{T}>");
  // 	  hratiocorrrefpt1D[nj][ncen-1][ip]->GetXaxis()->SetTitleFont(42);
  // 	  hratiocorrrefpt1D[nj][ncen-1][ip]->GetXaxis()->SetLabelFont(42);
  // 	  hratiocorrrefpt1D[nj][ncen-1][ip]->GetXaxis()->SetLabelSize(0.08);
  // 	  hratiocorrrefpt1D[nj][ncen-1][ip]->GetXaxis()->SetTitleSize(0.07);
	  
  // 	  hratiocorrrefpt1D[nj][ncen-1][ip]->GetYaxis()->SetTitle("");
  // 	  hratiocorrrefpt1D[nj][ncen-1][ip]->GetYaxis()->SetTitleFont(42);
  // 	  hratiocorrrefpt1D[nj][ncen-1][ip]->GetYaxis()->SetLabelFont(42);
  // 	  hratiocorrrefpt1D[nj][ncen-1][ip]->GetYaxis()->SetLabelSize(0.08);
	  
  // 	  hratiocorrrefpt1D[nj][ncen-1][ip]->SetMarkerStyle(24);
  // 	  hratiocorrrefpt1D[nj][ncen-1][ip]->SetMarkerColor(2);
  // 	  hratiocorrrefpt1D[nj][ncen-1][ip]->SetLineColor(2);
  // 	  hratiocorrrefpt1D[nj][ncen-1][ip]->SetMarkerSize(0.8);

  // 	}	  

  // 	// TPaveStats *ps = (TPaveStats*)  hratiocorrrefpt1D[nj][ic][ip]->GetListOfFunctions()->FindObject("stats");
  // 	// ps->SetX1NDC(0.50);
  // 	// ps->SetY1NDC(0.41);       
  // 	// ps->SetX2NDC(0.95);
  // 	// ps->SetY2NDC(0.79);
  // 	// ps->SetTextFont(42);
  // 	// ps->Draw();

  // 	hratiocorrrefpt1D[nj][ncen-1][ip]->Draw("psame");  
  // 	c99[nj][ic]->Update();

  // 	TPaveText *pt   = new TPaveText(0.4524683,0.8914759,0.7023389,0.9597512,"brNDC");
  // 	pt->SetBorderSize(0);
  // 	pt->SetFillColor(10);
  // 	pt->SetTextFont(42);
  // 	TText *text = pt->AddText(Form("%0.0f < p_{T} (GeV/c) < %0.0f",ptbins[ip], ptbins[ip+1]));
  // 	text->SetTextSize(0.07);
  // 	pt->Draw();


	
	// if(ic<ncen-1){
	//   hratiocorrrefpt1D[nj][ic][ip]->SetMaximum(25.634);
	//   hratiocorrrefpt1D[nj][ic][ip]->SetMinimum(1e-09);
	//   hratiocorrrefpt1D[nj][ic][ip]->SetTitle(0);
	//   hratiocorrrefpt1D[nj][ic][ip]->GetXaxis()->SetTitle("<reco jet p_{T} / gen jet p_{T}>");
	//   hratiocorrrefpt1D[nj][ic][ip]->GetXaxis()->SetTitleFont(42);
	//   hratiocorrrefpt1D[nj][ic][ip]->GetXaxis()->SetLabelFont(42);
	//   hratiocorrrefpt1D[nj][ic][ip]->GetXaxis()->SetLabelSize(0.08);
	//   hratiocorrrefpt1D[nj][ic][ip]->GetXaxis()->SetTitleSize(0.07);
	
	//   hratiocorrrefpt1D[nj][ic][ip]->GetYaxis()->SetTitle("");
	//   hratiocorrrefpt1D[nj][ic][ip]->GetYaxis()->SetTitleFont(42);
	//   hratiocorrrefpt1D[nj][ic][ip]->GetYaxis()->SetLabelFont(42);
	//   hratiocorrrefpt1D[nj][ic][ip]->GetYaxis()->SetLabelSize(0.08);

	//   hratiocorrrefpt1D[nj][ic][ip]->SetMarkerStyle(24);
	//   hratiocorrrefpt1D[nj][ic][ip]->SetMarkerColor(1);
	//   hratiocorrrefpt1D[nj][ic][ip]->SetLineColor(1);
	//   hratiocorrrefpt1D[nj][ic][ip]->SetMarkerSize(0.9);
	//   hratiocorrrefpt1D[nj][ic][ip]->Draw("p");  

	//   hratiorawrefpt1D[nj][ic][ip]->SetMarkerStyle(24);
	//   hratiorawrefpt1D[nj][ic][ip]->SetMarkerColor(4);
	//   hratiorawrefpt1D[nj][ic][ip]->SetLineColor(4);
	//   hratiorawrefpt1D[nj][ic][ip]->SetMarkerSize(0.9);
	//   //hratiorawrefpt1D[nj][ic][ip]->Draw("psame");  
	
	//   c99[nj][ic]->Update();
	//   TPaveStats *ps = (TPaveStats*)  hratiocorrrefpt1D[nj][ic][ip]->GetListOfFunctions()->FindObject("stats");
	//   ps->SetX1NDC(0.50);
	//   ps->SetY1NDC(0.41);       
	//   ps->SetX2NDC(0.95);
	//   ps->SetY2NDC(0.79);
	//   ps->SetTextFont(42);
	//   ps->Draw();

	// }else{
	  
	

	
	// TPaveText *pt   = new TPaveText(0.4524683,0.8914759,0.7023389,0.9597512,"brNDC");
	// pt->SetBorderSize(0);
	// pt->SetFillColor(10);
	// pt->SetTextFont(42);
	// TText *text = pt->AddText(Form("%0.0f < p_{T} (GeV/c) < %0.0f",ptbins[ip], ptbins[ip+1]));
	// text->SetTextSize(0.07);
	// pt->Draw();
	
	// if(ipad==1){
	//   TPaveText *pt1 = new TPaveText(0.6044166,0.2194909,0.8644171,0.3668644,"brNDC");
	//   pt1->SetBorderSize(0);
	//   pt1->SetFillColor(10);
	//   pt1->SetTextFont(42);
	//   TText *text1 = 0;
	//   if(ic==ncen-1)text1 = pt1->AddText("ak3PF");
	//   else text1 = pt1->AddText(calgo[nj]);
	//   text1->SetTextSize(0.09);	
	//   pt1->Draw();
	  
	//   TPaveText *pt2 = new TPaveText(0.60104,0.8160025,0.8475339,0.9142515,"brNDC");
	//   pt2->SetBorderSize(0);
	//   pt2->SetFillColor(10);
	//   pt->SetTextFont(42);
	//   TText *text2 = pt2->AddText(Form("%s",ccent[ic]));
	//   text2->SetTextSize(0.08);
	//   pt2->Draw();
	// }
  //     }
  //   }
  // }
  // return 0;


  /*
  //! Difference between Gaussian fits mean/sigma and ArM/RMS
  cout<<"Sigma and RMS ....." <<endl;
  ipad=0;
  ialg=1;
  TCanvas *c7 = new TCanvas("c7","Sigma and RMSFit resolution",15,131,1885,546);
  makeMultiPanelCanvas(c7,ncen,2,0.0,0.0,0.22,0.22,0.02,0);
  //c7->SetGridx();
  TLegend *l7 = new TLegend(0.2453033,0.6247655,0.4838159,0.9399625,NULL,"BRNDC");
  l7->SetHeader("");
  l7->SetBorderSize(0);
  l7->SetTextFont(42);
  l7->SetTextSize(0.09);
  l7->SetLineColor(1);
  l7->SetLineStyle(1);
  l7->SetLineWidth(1);
  l7->SetFillColor(10);
  l7->SetFillStyle(1001);
  l7->SetHeader("");			  
  l7->AddEntry(hRMS[ialg][0],"RMS","p");
  l7->AddEntry(hSigma[ialg][0],"Gaussian Sigma","p");

  for(int icen=ncen-1;icen>=0;icen--){
    c7->cd(++ipad); 
    //gPad->SetGridx();
    //gPad->SetLogx();

    if(ipad==1){//! pp

      hSigma[ialg][icen]->SetStats(0);
      hSigma[ialg][icen]->SetMarkerStyle(20);
      hSigma[ialg][icen]->SetMarkerColor(1);
      hSigma[ialg][icen]->SetLineColor(1);
      hSigma[ialg][icen]->SetMarkerSize(1.0);

      hRMS[ialg][icen]->SetStats(0);      
      hRMS[ialg][icen]->SetMarkerStyle(24);
      hRMS[ialg][icen]->SetMarkerColor(1);
      hRMS[ialg][icen]->SetLineColor(1);
      hRMS[ialg][icen]->SetMarkerSize(1.2);
      
      hRMS[ialg][icen]->Draw("p");
      hSigma[ialg][icen]->Draw("psame");
      hRatio_RMS[ialg][icen] = (TH1F*)hRMS[ialg][icen]->Clone(Form("hRatio_RMS%d_%d",1,icen));
      hRatio_RMS[ialg][icen]->Divide(hSigma[ialg][icen]);

    }else{//! pbpb

      hSigma[ialg][icen]->SetStats(0);      
      hSigma[ialg][icen]->SetMarkerStyle(20);
      hSigma[ialg][icen]->SetMarkerColor(1);
      hSigma[ialg][icen]->SetLineColor(1);
      hSigma[ialg][icen]->SetMarkerSize(1.0);

      hRMS[ialg][icen]->SetStats(0);            
      hRMS[ialg][icen]->SetMarkerStyle(24);
      hRMS[ialg][icen]->SetMarkerColor(1);
      hRMS[ialg][icen]->SetLineColor(1);
      hRMS[ialg][icen]->SetMarkerSize(1.2);
      
      hRMS[ialg][icen]->Draw("p");
      hSigma[ialg][icen]->Draw("psame");
      


      hRatio_RMS[ialg][icen] = (TH1F*)hRMS[ialg][icen]->Clone(Form("hRatio_RMS%d_%d",1,icen));
      hRatio_RMS[ialg][icen]->Divide(hSigma[ialg][icen]);
    }
  
    gPad->Update();
    
    TPaveText *pt1 = new TPaveText(0.2846445,0.06345905,0.4451356,0.2028512,"brNDC");
    pt1->SetBorderSize(0);
    pt1->SetFillColor(10);
    pt1->SetTextFont(42);
    TText *text1 = pt1->AddText(Form("%s",ccent[icen]));
    text1->SetTextSize(0.08);
    pt1->Draw();

    if(ipad==1)l7->Draw();
    
    gPad->Update();

    c7->cd(ipad+ncen);
    //gPad->SetLogx();
    //gPad->SetGridx();
    
    TLine *l2 = new TLine(xmin,1,xmax+250,1);
    l2->SetLineWidth(1);
    l2->SetLineStyle(2);

    hRatio_RMS[ialg][icen]->SetStats(0);
    hRatio_RMS[ialg][icen]->SetMaximum(1.546);
    hRatio_RMS[ialg][icen]->SetMinimum(0.846);
    hRatio_RMS[ialg][icen]->GetYaxis()->SetTitle("RMS/Sigma");
    hRatio_RMS[ialg][icen]->Draw("p");
    l2->Draw();
  }

  // c7->SaveAs("AN/JER/Fits_RMSSigma_JER.png");
  // c7->SaveAs("AN/JER/Fits_RMSSigma_JER.C");
  // c7->SaveAs("AN/JER/Fits_RMSSigma_JER.eps");
  // c7->SaveAs("AN/JER/Fits_RMsSigma_JER.pdf");


  cout<<endl;
  cout<<"ArM and Mean ....." <<endl;
  //! Mean
  ipad=0;
  TCanvas *c8 = new TCanvas("c8","ArM and Mean",15,131,1885,546);
  makeMultiPanelCanvas(c8,ncen,2,0.0,0.0,0.22,0.22,0.02,0);

  TLegend *l8 = new TLegend(0.2453033,0.6247655,0.4838159,0.9399625,NULL,"BRNDC");
  l8->SetHeader("");
  l8->SetBorderSize(0);
  l8->SetTextFont(42);
  l8->SetTextSize(0.09);
  l8->SetLineColor(1);
  l8->SetLineStyle(1);
  l8->SetLineWidth(1);
  l8->SetFillColor(10);
  l8->SetFillStyle(1001);
  l8->SetHeader("");			  
  l8->AddEntry(hArM[ialg][0],"Ar Mean","p");
  l8->AddEntry(hMean[ialg][0],"Gaussian Mean","p");

  for(int icen=ncen-1;icen>=0;icen--){
    
    c8->cd(++ipad); 

    if(ipad==1){//! pp
      hMean[ialg][icen]->SetStats(0);
      hMean[ialg][icen]->SetMarkerStyle(20);
      hMean[ialg][icen]->SetMarkerColor(1);
      hMean[ialg][icen]->SetLineColor(1);
      hMean[ialg][icen]->SetMarkerSize(1.0);
      
      hArM[ialg][icen]->SetStats(0);
      hArM[ialg][icen]->SetMarkerStyle(24);
      hArM[ialg][icen]->SetMarkerColor(1);
      hArM[ialg][icen]->SetLineColor(1);
      hArM[ialg][icen]->SetMarkerSize(1.2);

      hArM[ialg][icen]->Draw("p");
      hMean[ialg][icen]->Draw("psame");
      
      hRatio_Mean[ialg][icen] = (TH1F*)hArM[ialg][icen]->Clone(Form("hRatio_ArM%d_%d",ialg,icen));
      hRatio_Mean[ialg][icen]->Divide(hMean[ialg][icen]);
      
    }else{

      hMean[ialg][icen]->SetStats(0);
      hMean[ialg][icen]->SetMarkerStyle(20);
      hMean[ialg][icen]->SetMarkerColor(1);
      hMean[ialg][icen]->SetLineColor(1);
      hMean[ialg][icen]->SetMarkerSize(1.0);

      hArM[ialg][icen]->SetStats(0);
      hArM[ialg][icen]->SetMarkerStyle(24);
      hArM[ialg][icen]->SetMarkerColor(1);
      hArM[ialg][icen]->SetLineColor(1);
      hArM[ialg][icen]->SetMarkerSize(1.2);

      hArM[ialg][icen]->Draw("p");      
      hMean[ialg][icen]->Draw("psame");

      hRatio_Mean[ialg][icen] = (TH1F*)hArM[ialg][icen]->Clone(Form("hRatio_Mean%d_%d",ialg,icen));
      hRatio_Mean[ialg][icen]->Divide(hMean[ialg][icen]);
    }      

    TPaveText *pt1 = new TPaveText(0.2846445,0.06345905,0.4451356,0.2028512,"brNDC");
    pt1->SetBorderSize(0);
    pt1->SetFillColor(10);
    pt1->SetTextFont(42);
    TText *text1 = pt1->AddText(Form("%s",ccent[icen]));
    text1->SetTextSize(0.08);
    pt1->Draw();

    if(ipad==1)l8->Draw();

    c8->cd(ipad+ncen);
    TLine *l2 = new TLine(xmin,1,xmax+250,1);
    l2->SetLineWidth(1);
    l2->SetLineStyle(2);

    hRatio_Mean[ialg][icen]->SetStats(0);
    hRatio_Mean[ialg][icen]->SetMaximum(1.046);
    hRatio_Mean[ialg][icen]->SetMinimum(0.946);
    hRatio_Mean[ialg][icen]->GetYaxis()->SetTitle("ArM/Mean");
    hRatio_Mean[ialg][icen]->Draw("p");
    l2->Draw();
  }
  // c8->SaveAs("AN/JER/Fits_RMSSigma_JES.png");
  // c8->SaveAs("AN/JER/Fits_RMSSigma_JES.C");
  // c8->SaveAs("AN/JER/Fits_RMSSigma_JES.eps");
  // c8->SaveAs("AN/JER/Fits_RMsSigma_JES.pdf");
  */
  return 0;

}
void DefineHisto(const char *cver, int nj, int ic, TH1F *&h1, const char *name, const char *caption)
{
  h1 = new TH1F(Form("%s_%s_%d_%d",name,cver,nj,ic),Form("CMSSW %s %s %s %s %d",cver,caption,ksp[ic],calgo[nj],ic),nbins,ptbins);
  MakeHist(h1,statop,"p_{T}^{GenJet} (GeV/c)","<p_{T}^{RecoJet}/p_{T}^{GenJet}>");
}

void CalJesJer(const char *cver, int nj, int ic, TH2 *h2d, TH1F **h1d, TH1 *ham, TH1 *hr, TH1 *hm, TH1 *hs)
{

  if(nj==0){
    fout = new TFile("output.root","RECREATE");
  }
  for(int ip=0;ip<nbins;ip++){
    int lbin = h2d->GetXaxis()->FindBin(ptbins[ip]);
    int hbin = h2d->GetXaxis()->FindBin(ptbins[ip+1]);
    h1d[ip]  = (TH1F*)h2d->ProjectionY(Form("hratiocorrrefpt1D_%d_%d_%d",nj,ic,ip),lbin,hbin,"e");
    if(h1d[ip]->GetEntries()<maxEntry)continue;
    //h1d[ip]->Rebin(5);
    if(nj==0){
      fout->cd();
      h1d[ip]->Write();
    }
    FillMeanSigma(ip,h1d[ip],ham,hr,hm,hs);
  }
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
      h1F->Fit("gaus","MLLRQ0+","",fitmin,fitmax);
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
      
      //h1F->Fit("f1",fopt,"",fitmin,fitmax);
      h1F->Fit("f1","MLLRQ0+","",fitmin,fitmax);
      hm->SetBinContent(ip+1,f1->GetParameter(2));
      hs->SetBinContent(ip+1,f1->GetParameter(1));
      
      if(strcmp(fopt,"MLLRQ0+")==0){
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
