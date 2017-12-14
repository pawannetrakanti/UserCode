#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TKey.h>
#include <TH1F.h>
#include <TF1.h>
#include <TVirtualFitter.h>
#include <TPaveStats.h>
#include <TMath.h>
#include <TSpectrum.h>
#include <TCanvas.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include <cmath>
#include <algorithm>
#include <Math/MinimizerOptions.h>

#include "utilities.h"

using namespace std;

void fit_double_gaussian(TH1F*& hrsp);

/// default fit with gaussian in niter iteration of mean  
void fit_gaussian(TH1F*& hrsp,
                  const double nsigma,
                  const double jtptmin,
                  const int niter,
		  const string alg,
		  const int stat
);

/// optional double sided crystal ball fit to response distributions
int fit_dscb(TH1F*& hrsp,
             const double nsigma,
             const double jtptmin,
             const int niter,
             const string alg,
	     const int stat
);

/// double sided crystal ball function definition  
double fnc_dscb(double*xx,double*pp);

double fit_rsp(double *x, double *p);
double fit_rsp(double *x, double *p)
{
  double xval = x[0];
  // if(xval > p[7] && xval < p[8]){
  //   TF1::RejectPoint();
  //   return 0;
  // }
  double result = p[0]-p[1]/(pow(log10(xval)-p[2],2)+p[2])-p[3]*exp((-p[4]*(log10(xval)-p[5])*(log10(xval)-p[5]))+(p[6]*(log10(xval)-p[5])));
  return result;

}

void adjust_fitrange(TH1* h,double& min,double& max);
void adjust_fitrange_new(TH1* h,double& min,double& max);

template <class T>
bool from_string(T& t,
                 const std::string& s,
                 std::ios_base& (*f)(std::ios_base&))
{
  std::istringstream iss(s);
  return !(iss >> f >> t).fail();
}

double GetHistMax(TH1 *h);

double ptbins[] ={
  1,5,10,15,20,25,30,35,40,45,57,72,90,120,150,200,300,400,550,750,1000,1500
};


const int nbins = sizeof(ptbins)/sizeof(double) - 1;

const int knj = 1;
const char *calgo[knj]={"ak4"};

int PFL3(const int scene=5, std::string ctype="pf")
{

  bool iSave=false;
  LoadStyle();
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);

  double minx=1.0;
  string outdir="txtfiles/";
  cout << " # of bins in pt : "  << nbins << endl;
  //return 0;

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000); 
  TVirtualFitter::SetDefaultFitter("Minuit");
  TH1F *hrsp[knj][nbins];
  TGraphErrors *grsp[knj], *gcor[knj];
  TH1F *hMean[knj], *hRMS[knj];
  
  TFile *fin = 0;
  // if( scene == 3 )fin = new TFile("input/jra_ecal-aged-hcal-aged-sipm-nonaged_scenario3_f.root","r");
  // else if( scene == 5 )fin = new TFile("input/jra_ecal-aged-hcal-noaged-sipm-aged_scenario5.root","r");

  // if( scene == 3 )fin = new TFile("input/jra_ecal-aged-hcal-aged-sipm-nonaged_scenario3.root","r");
  // if( scene == 3 )fin = new TFile("input/jra_ecal-aged-hcal-aged-sipm-nonaged_scenario3_jec.root","r");
  // else if( scene == 5)fin = new TFile("input/jra_ecal-aged-hcal-noaged-sipm-aged_scenario5_f.root","r");

  //fin = new TFile("input/jra_ecal-aged-hcal-aged-sipm-nonaged_scenario3_jec.root","r");

  fin = new TFile("input/jra_ecal-aged-hcal-noaged-sipm-aged_scenario5.root","r");

  double refptmin=1;
  double refptmax=1500;

  TLegend *l0 = getLegend(0.6349664,0.3265565,0.8491835,0.7128335);
  l0->SetHeader(Form("Anti-k_{T}, %s",ctype.c_str()));
  TLegend *l1 = getLegend(0.1500646,0.5052258,0.3653091,0.7490942);
  l1->SetHeader("");

  int njmin=0;
  int njmax=knj;

  TF1 *fgaus=0;
  TH1F *hrefpt=0, *hjetpt=0;
  for(int i=njmin; i<njmax; i++){

    grsp[i] = new TGraphErrors(nbins);
    gcor[i] = new TGraphErrors(nbins);

    grsp[i]->SetName(Form("L3RspVsRefPt_%s",calgo[i]));
    grsp[i]->SetTitle(Form("L3RspVsRefPt %s",calgo[i]));
    grsp[i]->SetLineColor(1);
    grsp[i]->SetMarkerColor(1);
    grsp[i]->SetMarkerStyle(20);
    grsp[i]->SetMarkerSize(1.2);

    gcor[i]->SetName(Form("L3CorVsJetPt_%s",calgo[i]));
    gcor[i]->SetTitle(Form("L3CorVsJetPt %s",calgo[i]));
    gcor[i]->SetLineColor(1);
    gcor[i]->SetMarkerColor(1);
    gcor[i]->SetMarkerStyle(20);
    gcor[i]->SetMarkerSize(1.2);

    hMean [i] = new TH1F(Form("hMean_%d",i),Form("%s%s",calgo[i],ctype.c_str()),nbins,ptbins);
    hMean [i]->SetMarkerColor(1);
    hMean [i]->SetMarkerStyle(20);
    hMean [i]->SetLineColor(1);
    hMean [i]->GetXaxis()->SetMoreLogLabels();
    hMean [i]->GetXaxis()->SetNoExponent();

    hRMS [i] = new TH1F(Form("hRMS_%d",i),Form("%s%s",calgo[i],ctype.c_str()),nbins,ptbins);
    hRMS [i]->SetMarkerColor(1);
    hRMS [i]->SetMarkerStyle(20);
    hRMS [i]->SetLineColor(1);


    std::ostringstream strs;
    //strs << calgo[i] << "pf/";
    strs << calgo[i] << ctype.c_str() << "/";

    cout << " Running for " << strs.str().c_str() << endl;
    
    
    for(int j=0; j<nbins; j++){
      //if(j<2)continue;
      std::ostringstream spt; 
      spt << ptbins[j] << "to" << ptbins[j+1];
      
      std::string sone  (strs.str()+"RefPt_Barrel_RefPt"+spt.str());
      std::string stwo  (strs.str()+"JetPt_Barrel_RefPt"+spt.str());
      std::string sthree(strs.str()+"RelRsp_Barrel_RefPt"+spt.str());

      hrefpt     = (TH1F*)fin->Get(Form("%s",sone.c_str()));
      hjetpt     = (TH1F*)fin->Get(Form("%s",stwo.c_str()));
      hrsp[i][j] = (TH1F*)fin->Get(Form("%s",sthree.c_str()));

      if(hrsp[i][j]->Integral()==0) continue;

      assert(hrefpt->GetEntries()>0 && hjetpt->GetEntries()>0);
      
      double refpt   =hrefpt->GetMean();
      double erefpt  =hrefpt->GetMeanError();
      double jetpt   =hjetpt->GetMean();
      double ejetpt  =hjetpt->GetMeanError();

      hrsp[i][j]->SetLineColor(1);
      hrsp[i][j]->SetMarkerColor(1);
      hrsp[i][j]->SetMarkerStyle(20);
      if(j==0)l0->AddEntry(hrsp[i][j],Form("R = %0.1f",(i+1)*0.1),"p");
      
      double rms   = hrsp[i][j]->GetRMS(); 
      double peak  = hrsp[i][j]->GetMean();
      double epeak = hrsp[i][j]->GetMeanError();

      //if(i==0)cout <<"  refpt : " << refpt << " rawpt : " << jetpt << endl;
      fit_gaussian(hrsp[i][j], 1.0, 1, 4, Form("%s%s",calgo[i],ctype.c_str()), 0);
      
      //if ( ptbins[j] < 23 ) fit_dscb(hrsp[i][j], 1.1, 15, 4, Form("%s%s",calgo[i],ctype.c_str()), 0);
      //else fit_dscb(hrsp[i][j], 1.0, 15, 4, Form("%s%s",calgo[i],ctype.c_str()), 0);
      //if( j != 0 )fit_dscb(hrsp[i][j], 1.00, 1, 4, Form("%s%s",calgo[i],ctype.c_str()), 0);
      

      TF1*  frelrsp = (TF1*)hrsp[i][j]->GetListOfFunctions()->Last();
      peak    = (frelrsp==0) ? hrsp[i][j]->GetMean()      :  frelrsp->GetParameter(1);
      epeak   = (frelrsp==0) ? hrsp[i][j]->GetMeanError() :  frelrsp->GetParError(1);

      //cout << ptbins[j] << " peak : " << peak << endl;
      // if( j == 1 ){
      // 	peak = 0.675;
      // 	epeak = 3.97703e-03;
      // }

      // if( j==0 ){
      // 	peak  = 9.12848e-01; //hrsp[i][j]->GetMean();
      // 	epeak = 1.26090e-03; //hrsp[i][j]->GetMeanError();	
      // }
      // peak    = (frelrsp==0) ? GetHistMax(hrsp[i][j])  :  frelrsp->GetParameter(1);
      // epeak   = (frelrsp==0) ? 0.02*peak               :  frelrsp->GetParError(1);

      double rsp     =peak;
      double ersp    =epeak;
      double cor     =1.0/rsp;
      double ecor    =cor*cor*ersp;

      // double sigma  = (frelrsp==0) ? hrsp[i][j]->GetRMS() :  frelrsp->GetParameter(2);
      // double esigma = (frelrsp==0) ? hrsp[i][j]->GetRMSError() :  frelrsp->GetParError(2);
      // hMean[i]->SetBinContent(j+1,peak);
      // hMean[i]->SetBinError(j+1,epeak);
      // hRMS[i]->SetBinContent(j+1,sigma/peak);
      // hRMS[i]->SetBinError(j+1,esigma);

      //assert(grsp[i]->GetN()==gcor[i]->GetN());
      //int ij = grsp[i]->GetN();

      grsp[i]->SetPoint     (j, refpt,  peak);
      grsp[i]->SetPointError(j, erefpt, epeak);
      gcor[i]->SetPoint     (j, jetpt,  cor);
      gcor[i]->SetPointError(j, ejetpt, ecor);
    }
  }
  //hrsp[2][15]->Draw("p");
  //return 0;
  if(gPad)gPad->Close();

  if( 0 ){
    //! 0 < |eta| < 0.5
    //! Anti k_{T} R=0.4 (PF);
    TCanvas *cClosureM = new TCanvas("cClosureM","Mean Closure",950,800);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.08);
    gPad->SetBottomMargin(0.15);
    gPad->SetLogx();
    hMean[0]->SetStats(0);
    hMean[0]->SetTitle("Scenario3 HCAL-AGED ak4pfl2l3");
    hMean[0]->GetXaxis()->SetRangeUser(20,1500);
    hMean[0]->GetXaxis()->SetTitle("p_{T}^{REF} [GeV]");
    hMean[0]->GetXaxis()->SetTitleFont(42);
    hMean[0]->GetXaxis()->SetLabelFont(42);
    hMean[0]->GetXaxis()->SetLabelSize(0.07);
    hMean[0]->GetXaxis()->SetTitleSize(0.07);
    hMean[0]->GetXaxis()->SetNdivisions(507);
    hMean[0]->GetXaxis()->SetNoExponent();
    hMean[0]->GetXaxis()->SetMoreLogLabels();
    hMean[0]->GetYaxis()->SetTitle("p_{T}/p_{T}^{REF}");
    hMean[0]->GetYaxis()->SetRangeUser(0.8,1.2);
    hMean[0]->GetYaxis()->SetTitleFont(42);
    hMean[0]->GetYaxis()->SetLabelFont(42);
    hMean[0]->GetYaxis()->SetLabelSize(0.07);
    hMean[0]->GetYaxis()->SetNdivisions(507);
    
    hMean[0]->Draw("p");
    TLine *line1 = new TLine(20,0.99,1500,0.99);
    line1->SetLineStyle(2);
    line1->Draw();
    TLine *line2 = new TLine(20,1.01,1500,1.01);
    line2->SetLineStyle(2);
    line2->Draw();
    
    const string s_sigma="sqrt([0]*abs([0])/(x*x) + [1]*[1]*pow(x,[3]) + [2]*[2])";
    TF1 *fres = new TF1("fres",s_sigma.c_str(),20,1250);
    fres->SetLineWidth(2);
    fres->SetLineColor(2);  
    fres->SetParameters(3.5,0.5,0.03,-1);
    
    TCanvas *cClosureS = new TCanvas("cClosureS","Sigma Closure",950,800);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.08);
    gPad->SetBottomMargin(0.15);
    gStyle->SetOptStat(1);
    gStyle->SetOptFit(1);
    hRMS[0]->SetStats(1);
    hRMS[0]->SetTitle("Scenario3 HCAL-AGED ak4pfl2l3");
    hRMS[0]->GetXaxis()->SetRangeUser(20,1500);
    hRMS[0]->GetXaxis()->SetTitle("p_{T}^{REF} [GeV]");
    hRMS[0]->GetXaxis()->SetTitleFont(42);
    hRMS[0]->GetXaxis()->SetLabelFont(42);
    hRMS[0]->GetXaxis()->SetLabelSize(0.07);
    hRMS[0]->GetXaxis()->SetTitleSize(0.07);
    hRMS[0]->GetXaxis()->SetNdivisions(507);
    hRMS[0]->GetXaxis()->SetNoExponent();
    hRMS[0]->GetXaxis()->SetMoreLogLabels();
    hRMS[0]->GetYaxis()->SetTitle("p_{T}/p_{T}^{REF}");
    hRMS[0]->GetYaxis()->SetRangeUser(0.0,0.4);
    hRMS[0]->GetYaxis()->SetTitleFont(42);
    hRMS[0]->GetYaxis()->SetLabelFont(42);
    hRMS[0]->GetYaxis()->SetLabelSize(0.07);
    hRMS[0]->GetYaxis()->SetNdivisions(507);

    hRMS[0]->Fit(fres,"R+");
    hRMS[0]->Draw("p");
    
    return 0;
  }
  



  int ipad=0;
  //! 0 - 20 GeV
  TCanvas *c98[knj];
  for(int nj=njmin;nj<njmax;nj++){
    c98[nj] = new TCanvas(Form("c99_%d",nj),Form("%s%s Fitting plots",calgo[nj],ctype.c_str()),100,102,1399,942);
    c98[nj]->Divide(5,4,0,0);
    ipad=0;
    for(int ip=0;ip<nbins;ip++){      
      c98[nj]->cd(++ipad);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.01);
      gPad->SetBottomMargin(0.15);

      hrsp[nj][ip]->SetTitle(0);
      hrsp[nj][ip]->GetXaxis()->SetTitle("< REC p_{T} / REF p_{T}>");
      hrsp[nj][ip]->GetXaxis()->SetTitleFont(42);
      hrsp[nj][ip]->GetXaxis()->SetLabelFont(42);
      hrsp[nj][ip]->GetXaxis()->SetLabelSize(0.08);
      hrsp[nj][ip]->GetXaxis()->SetTitleSize(0.07);
      hrsp[nj][ip]->GetXaxis()->SetNdivisions(507);
      hrsp[nj][ip]->GetYaxis()->SetTitle("");
      hrsp[nj][ip]->GetYaxis()->SetTitleFont(42);
      hrsp[nj][ip]->GetYaxis()->SetLabelFont(42);
      hrsp[nj][ip]->GetYaxis()->SetLabelSize(0.08);
      hrsp[nj][ip]->GetYaxis()->SetNdivisions(507);
      
      hrsp[nj][ip]->SetMarkerStyle(20);
      hrsp[nj][ip]->SetMarkerColor(1);
      hrsp[nj][ip]->SetLineColor(1);
      hrsp[nj][ip]->SetMarkerSize(1.1);
      hrsp[nj][ip]->Draw("hist");  
      
      TF1 *fdscb = (TF1*)hrsp[nj][ip]->GetFunction("fdscb");
      if(fdscb){
   	fdscb->SetLineWidth(2);
   	fdscb->SetLineColor(2);
   	fdscb->Draw("lsame");
      }
      TF1 *fgaus = (TF1*)hrsp[nj][ip]->GetFunction("fgaus");
      if(fgaus){
   	fgaus->SetLineWidth(2);
   	fgaus->SetLineColor(4);
   	fgaus->Draw("lsame");
      }
      
      std::ostringstream strs; 
      strs << ptbins[ip] << "< p_{T} (GeV/c) <" << ptbins[ip+1];
      std::string spt = strs.str();
      
      if(ipad==1){drawText2(Form("%s%s",calgo[nj],ctype.c_str()),0.28,0.90,18);      
	
   	drawText2(spt.c_str(),0.22,0.80,15);		
      } else drawText2(spt.c_str(),0.17,0.80,15);		
    }
    //if( strcmp(calgo[nj],"ak1") !=0 )c98[nj]->Close();
  }
  //return 0;

  //cout << endl;
  //cout << endl;

  cout << "Fitting the L3 Response "  << endl;

  //cout << endl;
  //cout << endl;


  TF1 *fitrsp=0, *fitcor=0;
  //std::string fcn_rsp="[0]-[1]/(pow(log10(x),2)+[2])-[3]*exp((-[4]*(log10(x)-[5])*(log10(x)-[5]))+([6]*(log10(x)-[5])*(log10(x)-[5])))";
  //std::string fcn_cor="[0]+[1]/(pow(log10(x),2)+[2])+[3]*exp((-[4]*(log10(x)-[5])*(log10(x)-[5]))+([6]*(log10(x)-[5])*(log10(x)-[5])))";

  //! Orgx
   std::string fcn_rsp="[0]-[1]/(pow(log10(x),2)+[2])-[3]*exp((-[4]*(log10(x)-[5])*(log10(x)-[5]))+([6]*(log10(x)-[5])))";
   std::string fcn_cor="[0]+[1]/(pow(log10(x),2)+[2])+[3]*exp((-[4]*(log10(x)-[5])*(log10(x)-[5]))+([6]*(log10(x)-[5])))";

  TCanvas *c11 = new TCanvas("c11","L3 Response",950,800);

  for(int i=njmin; i<njmax; i++){  
    c11->cd(i+1);
   //gPad->SetLogx();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->Update();
    grsp[i]->GetXaxis()->SetTitle("p_{T}^{REF}");
    grsp[i]->GetXaxis()->SetTitleFont(42);
    grsp[i]->GetXaxis()->SetLabelFont(42);
    grsp[i]->GetXaxis()->SetLabelSize(0.06);
    grsp[i]->GetXaxis()->SetTitleSize(0.06);
    grsp[i]->GetXaxis()->SetNdivisions(507);
    grsp[i]->GetXaxis()->SetNoExponent();
    grsp[i]->GetXaxis()->SetMoreLogLabels();
    grsp[i]->GetYaxis()->SetTitle("p_{T}/p_{T}^{REF}");
    grsp[i]->GetYaxis()->SetTitleFont(42);
    grsp[i]->GetYaxis()->SetLabelFont(42);
    grsp[i]->GetYaxis()->SetLabelSize(0.06);
    grsp[i]->GetYaxis()->SetTitleSize(0.06);
    grsp[i]->GetYaxis()->SetNdivisions(507);
    grsp[i]->GetYaxis()->SetRangeUser(0,1.160);

    //fitrsp = new TF1("fitrsp",fcn_rsp.c_str(),1.0,grsp[i]->GetX()[grsp[i]->GetN()-1]);
    fitrsp = new TF1("fitrsp",fcn_rsp.c_str(),grsp[i]->GetX()[2], grsp[i]->GetX()[grsp[i]->GetN()-1]);

    //fitrsp = new TF1("fitrsp",fcn_rsp.c_str(),20, grsp[i]->GetX()[grsp[i]->GetN()-1]);
    
    // for(int ij=0; ij<grsp[i]->GetN(); ij++){
    //   cout << " fit_min " << grsp[i]->GetX()[ij] << endl;
    // }

    //! Org //! used for Scenario3
    // fitrsp->SetParameter(0, 0.96);
    // fitrsp->SetParameter(1, 0.033);
    // fitrsp->SetParameter(2,-0.70);
    // fitrsp->SetParameter(3, 0.02);
    // fitrsp->SetParameter(4, 1.02);
    // fitrsp->SetParameter(5, 2.70);
    // fitrsp->SetParameter(6, 0.016);

    //! Scenario5
    fitrsp->SetParameter(0, 0.80);
    fitrsp->SetParameter(1, 0.033);
    fitrsp->SetParameter(2,-0.70);
    fitrsp->SetParameter(3, 0.02);
    fitrsp->SetParameter(4, 1.02);
    fitrsp->SetParameter(5, 2.70);
    fitrsp->SetParameter(6, 0.016);
    
    // fitrsp->SetParLimits(2,-0.90,1.10);
    // fitrsp->SetParLimits(3,0.003,0.08);
    fitrsp->SetLineWidth(2);

    //grsp[i]->Fit(fitrsp,"QR","",15.0, grsp[i]->GetX()[grsp[i]->GetN()-1]);
    grsp[i]->Fit(fitrsp,"QR");

    grsp[i]->SetMaximum(1.2);
    grsp[i]->SetMinimum(0.2);


    grsp[i]->Draw("ap");
    drawText2(Form("%s%s",calgo[i],ctype.c_str()),0.20,0.78,21);


    gPad->Update();
    TPaveStats *strsp = (TPaveStats*)grsp[i]->FindObject("stats");
    strsp->SetX1NDC(0.46);
    strsp->SetY1NDC(0.18);
    strsp->SetX2NDC(0.87);
    strsp->SetY2NDC(0.59);
    strsp->Draw();
    //delete fitrsp;
    //continue;
  }    

  
  TCanvas *c12 = new TCanvas("c12","L3 Corrections",950,800);
  for(int i=njmin; i<njmax; i++){  
    c12->cd(i+1);
    // cout << "Fitting the L3 Corrections "  << endl;
    // cout << endl;
  
    fitcor = new TF1("fitcor",fcn_cor.c_str(),
    		     gcor[i]->GetX()[2],gcor[i]->GetX()[gcor[i]->GetN()-1]);
    //fitcor = new TF1("fitcor",fcn_cor.c_str(),
    //		     1.0,gcor[i]->GetX()[gcor[i]->GetN()-1]);
    
    //! Org Scenario3
    // fitcor->SetParameter(0,1.04);
    // fitcor->SetParameter(1,.033);
    // fitcor->SetParameter(2,-0.70);
    // fitcor->SetParameter(3,0.02);
    // fitcor->SetParameter(4,1.02);
    // fitcor->SetParameter(5,2.7);
    // fitcor->SetParameter(6,0.016);

    // fitcor->SetParameter(0,1.04);
    // fitcor->SetParameter(1,0.033);
    // fitcor->SetParameter(2,-0.70);
    // fitcor->SetParameter(3,0.020);
    // fitcor->SetParameter(4,1.020);
    // fitcor->SetParameter(5,2.10);
    // fitcor->SetParameter(6,0.013);

    //! Scenario5
    fitcor->SetParameter(0, 1.20);
    fitcor->SetParameter(1, 0.033);
    fitcor->SetParameter(2,-0.70);
    fitcor->SetParameter(3, 0.02);
    fitcor->SetParameter(4, 1.02);
    fitcor->SetParameter(5, 2.70);
    fitcor->SetParameter(6, 0.016);

    // fitcor->SetParLimits(2,-0.90,1.10);
    // fitcor->SetParLimits(3,0.003,0.08);

    //gcor[i]->Fit(fitcor,"QR+","",20,gcor[i]->GetX()[gcor[i]->GetN()-1]);
    gcor[i]->Fit(fitcor,"QR");
    
    string algn="4";
    string era=Form("PhaseIISummer16_25nsV3_MC_Scenario%d",scene);
    string txtfilename = outdir+era+"_L3Absolute_AK"+algn+"PF.txt";
    ofstream outf(txtfilename.c_str());
    outf.setf(ios::left);
    outf<<"{1 JetEta 1 JetPt "<<fcn_cor<<" Correction L3Absolute}"<<endl;
    outf<<setw(12)<<-5.191                // eta_min
     	<<setw(12)<<+5.191                // eta_max
     	<<setw(12)<<fitcor->GetNpar()+2   // number of parameters + 2
     	<<setw(12)<<4.0                   // minimum pT
     	<<setw(12)<<5000.0;               // maximum pT
    for(int p=0; p<fitcor->GetNpar(); p++){
      outf<<setw(12)<<fitcor->GetParameter(p); // p0-p6
    }
    outf.close();

    gcor[i]->SetMaximum(1.35);
    gcor[i]->SetMinimum(0.98);
    gcor[i]->Draw("ap");
    //gPad->SetLogx();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->Update();

    gcor[i]->GetXaxis()->SetTitle("p_{T}");
    gcor[i]->GetXaxis()->SetTitleFont(42);
    gcor[i]->GetXaxis()->SetLabelFont(42);
    gcor[i]->GetXaxis()->SetLabelSize(0.06);
    gcor[i]->GetXaxis()->SetTitleSize(0.06);
    gcor[i]->GetXaxis()->SetNdivisions(507);
    gcor[i]->GetXaxis()->SetNoExponent();
    gcor[i]->GetXaxis()->SetMoreLogLabels();
    gcor[i]->GetYaxis()->SetTitle("L3 Correction");
    gcor[i]->GetYaxis()->SetTitleFont(42);
    gcor[i]->GetYaxis()->SetLabelFont(42);
    gcor[i]->GetYaxis()->SetLabelSize(0.06);
    gcor[i]->GetYaxis()->SetTitleSize(0.06);
    gcor[i]->GetYaxis()->SetNdivisions(507);

    gcor[i]->GetXaxis()->SetMoreLogLabels();
    gcor[i]->GetXaxis()->SetNoExponent();
    drawText2(Form("%s%s",calgo[i],ctype.c_str()),0.20,0.78,21);


    TPaveStats *stcor = (TPaveStats*)gcor[i]->FindObject("stats");
    stcor->SetX1NDC(0.48);
    stcor->SetY1NDC(0.53);
    stcor->SetX2NDC(0.89);
    stcor->SetY2NDC(0.89);
    stcor->Draw();
  
    // cout << endl;
    //delete fitcor;
  }

  // if(iSave){
  //   c12->SaveAs("CorrectionPlots/L3Absolute_Corrections_UpgradeJEC_ak4pf.pdf");
  // }

  TFile *fout = new TFile(Form("l3pf_scene%d_input.root",scene),"RECREATE");
  for(int nj=njmin;nj<njmax;nj++){
    fout->mkdir(Form("%s%s",calgo[nj],ctype.c_str()),Form("%s%s",calgo[nj],ctype.c_str()));
    fout->cd(Form("%s%s",calgo[nj],ctype.c_str()));
    grsp[nj]->Write();
    gcor[nj]->Write();
    fout->cd("../");
  }
  fout->Close();

  return 0;
}
int fit_dscb(TH1F*& hrsp,
             const double nsigma,
             const double jtptmin,
             const int niter,
             const string alg,
	     const int stat
)
{
  if (0==hrsp) {
    cout<<"ERROR: Empty pointer to fit_dscb()"<<endl;return -1;
  }


  // first use a gaussian to constrain crystal ball gaussian core
  fit_gaussian(hrsp,nsigma,jtptmin,niter,alg,stat);
  TF1* fgaus = hrsp->GetFunction("fgaus");
  if (0==fgaus) {
    hrsp->GetListOfFunctions()->Delete();
    return -1;
  }

  // implementation of the low pt bias threshold 
  string histname = hrsp->GetName();
  double ptRefMax(1.0),rspMax(0.0);
  int pos1     = histname.find("RefPt");
  int pos2     = histname.find("to",pos1);
  string ss    = histname.substr(pos1+5,pos2);
  if (from_string(ptRefMax,ss,std::dec)) {
    if (histname.find("RelRsp")==0)
      rspMax = jtptmin/ptRefMax;
    if (histname.find("AbsRsp")==0)
      rspMax = jtptmin-ptRefMax;
  }

  double fitrange_min(0.10);
  double fitrange_max(6.00);

  adjust_fitrange_new(hrsp,fitrange_min,fitrange_max);
  
  TF1* fdscb = new TF1("fdscb",fnc_dscb,fitrange_min,fitrange_max,7);
  fdscb->SetLineWidth(2);
  fdscb->SetLineStyle(2);

  double norm = fgaus->GetParameter(0);
  double mean = fgaus->GetParameter(1);
  double sigma= fgaus->GetParameter(2);

  // double norm = 2*hrsp->GetMaximumStored();
  // double mean = hrsp->GetMean();
  // double sigma= hrsp->GetRMS();

  //double aone(2.0),atwo(2.0),pone(10.0),ptwo(10.0);
  //double aone(4.0),atwo(4.0),pone(15.0),ptwo(15.0);

  double aone(2.0),atwo(2.0),pone(5.0),ptwo(5.0);
  TVirtualFitter::SetDefaultFitter("Minuit");
  
  int fitstatus(0);
  for (int i=0;i<niter;i++) {
    fdscb->SetParameter(0,norm); // N
    fdscb->SetParameter(1,mean); // mean
    fdscb->SetParameter(2,sigma);// sigma
    fdscb->SetParameter(3,aone); // a1
    fdscb->SetParameter(4,pone); // p1
    fdscb->SetParameter(5,atwo); // a2
    fdscb->SetParameter(6,ptwo); // p2                

    fdscb->FixParameter(1,mean);
    fdscb->FixParameter(2,sigma);

    //! HI reco and pp reco
    if (i>0) fdscb->FixParameter(3,aone);
    else fdscb->SetParLimits(3,0.,20.);
    
    if (i>1) fdscb->FixParameter(5,atwo);
    else fdscb->SetParLimits(5,0.,20.);
    
    fdscb->SetParLimits(4,0.,60.);
    fdscb->SetParLimits(6,0.,60.);


    fitstatus = hrsp->Fit(fdscb,"RQ+");
    if (0==fitstatus) i=999;
    delete fdscb;
    fdscb = hrsp->GetFunction("fdscb");
    if (0==fdscb) return -1;

    norm  = fdscb->GetParameter(0);
    aone  = fdscb->GetParameter(3);
    pone  = fdscb->GetParameter(4);
    atwo  = fdscb->GetParameter(5);
    ptwo  = fdscb->GetParameter(6);

    fdscb->ResetBit(TF1::kNotDraw);
    // reset sigma and mean to gauss values...
    fdscb->SetParameter(1,fgaus->GetParameter(1));
    fdscb->SetParError (1,fgaus->GetParError(1));
    fdscb->SetParameter(2,fgaus->GetParameter(2));
    fdscb->SetParError (2,fgaus->GetParError(2));
  }
  if (0!=fitstatus){
    cout<<alg.c_str()<<" fit_fdscb() to "<<hrsp->GetName()
        <<" failed. Fitstatus: "<<fitstatus<<endl;
    hrsp->GetFunction("fdscb")->Delete();
  }
  return fitstatus;
}
void fit_gaussian(TH1F*& hrsp,
                  const double nsigma,
                  const double jtptmin,
                  const int niter,
		  const string alg,
		  const int stat
)
{
  if (0==hrsp) {
    cout<<alg.c_str()<<" ERROR: Empty pointer to fit_gaussian()"<<endl;return;
  }
  
  string histname = hrsp->GetName();
  double mean     = GetHistMax(hrsp);//hrsp->GetMean();
  double rms      = hrsp->GetRMS();
  double ptRefMax(1.0),rspMax(0.0);

  //cout << histname.c_str()<< " Mean : "  << mean << " \t  RMS  : " << rms << endl;

  double norm  = 2*hrsp->GetMaximumStored();
  double peak  = mean;
  double sigma = rms;
  int pos1     = histname.find("RefPt");
  int pos2     = histname.find("to",pos1);
  string ss    = histname.substr(pos1+5,pos2);
  if (from_string(ptRefMax,ss,std::dec)) {
    if (histname.find("RelRsp")==0)
      rspMax = jtptmin/ptRefMax;
    if (histname.find("AbsRsp")==0)
      rspMax = jtptmin-ptRefMax;
  }
  double xmin  = hrsp->GetXaxis()->GetXmin();
  double xmax  = hrsp->GetXaxis()->GetXmax();
  TF1* fitfnc(0); int fitstatus(-1);
  for (int iiter=0;iiter<niter;iiter++) {
    vector<double> vv;
    vv.push_back(rspMax);
    vv.push_back(xmin);
    vv.push_back(peak-nsigma*sigma);

    //double fitrange_min = peak - 0.01*peak;
    //double fitrange_max = peak + 0.01*peak;

    double fitrange_min = *std::max_element(vv.begin(),vv.end());
    double fitrange_max = std::min(xmax,peak+nsigma*sigma);

    // if( stat ){
    //   fitrange_min = 0.90;
    //   fitrange_max = 1.10;
    // }
    // else adjust_fitrange(hrsp,fitrange_min,fitrange_max);

    //adjust_fitrange(hrsp,fitrange_min,fitrange_max);

    fitfnc = new TF1("fgaus","gaus",fitrange_min,fitrange_max);
    fitfnc->SetParNames("N","#mu","#sigma");
    fitfnc->SetParameter(0,norm);
    fitfnc->SetParameter(1,peak);
    fitfnc->SetParameter(2,sigma);
    
    fitstatus = hrsp->Fit(fitfnc,"RQ");
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
  // if (0==fitstatus){
  //   cout<<"fit_gaussian() to "<<hrsp->GetName()    <<"  sucessful : " <<endl;
  // }          

  if (0!=fitstatus){
    cout<<alg.c_str()<<"  fit_gaussian() to "<<hrsp->GetName()
        <<" failed. Fitstatus: "<<fitstatus
        <<" - FNC deleted."<<endl;
    hrsp->GetListOfFunctions()->Delete();
  }
}


double fnc_dscb(double*xx,double*pp)
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
  double A1  = TMath::Power(p1/TMath::Abs(a1),p1)*TMath::Exp(-a1*a1/2);
  double A2  = TMath::Power(p2/TMath::Abs(a2),p2)*TMath::Exp(-a2*a2/2);
  double B1  = p1/TMath::Abs(a1) - TMath::Abs(a1);
  double B2  = p2/TMath::Abs(a2) - TMath::Abs(a2);

  double result(N);
  if      (u<-a1) result *= A1*TMath::Power(B1-u,-p1);
  else if (u<a2)  result *= TMath::Exp(-u*u/2);
  else            result *= A2*TMath::Power(B2+u,-p2);
  return result;
}


void fit_double_gaussian(TH1F*& hrsp)
{
  if (0==hrsp) {
    cout<<"ERROR: Empty pointer to fit_double_gaussian()"<<endl;return;
  }
  
  string histname = hrsp->GetName();
  double mean     = hrsp->GetMean();
  double rms      = hrsp->GetRMS();

  TSpectrum *spec = new TSpectrum(10);
  int nfound = spec->Search(hrsp,1,"new");
  if(nfound!=2)return;

  double *xpeaks = spec->GetPositionX();

  double peak1  = xpeaks[0];
  double bin1   = hrsp->FindBin(peak1);
  double norm1  = hrsp->GetBinContent(bin1);
  double sigma1 = 0.2;

  double peak2  = xpeaks[1];
  double bin2   = hrsp->FindBin(peak2);
  double norm2  = hrsp->GetBinContent(bin2);
  double sigma2 = 0.3;


  cout << " Mean  : "  << mean  << " \t  RMS  : " << rms    << endl;
  cout << " norm1 : "  << norm1 << " \t  norm2 : " << norm2 << endl;
  cout << " peak1 : "  << peak1 << " \t  sig1 : " << sigma1 << endl;
  cout << " peak2 : "  << peak2 << " \t  sig2 : " << sigma2 << endl;

  double fitrange_min = 0.2;
  double fitrange_max = 1.7;

  TF1* fitfnc(0); int fitstatus(-1);
  TF1 *fitg1(0), *fitg2(0);
  fitfnc = new TF1("fdgaus","gaus(0)+gaus(3)",fitrange_min,fitrange_max);
  fitfnc->SetLineColor(1);
  fitfnc->SetLineStyle(2);
  fitfnc->SetLineWidth(2);

  fitfnc->SetParNames("N_{1}", "#mu_{1}", "#sigma_{1}",
		      "N_{2}", "#mu_{2}", "#sigma_{2}");
  fitfnc->SetParameters(norm1, peak1, sigma1, 
   			norm2, peak2, sigma2); 


  //fitfnc->SetParLimits(0,0.0,2.0*norm1);
  //fitfnc->SetParLimits(1,peak1-3*sigma1,peak1+3*sigma1);
  //fitfnc->SetParLimits(2,0.1,1.0);

  //fitfnc->SetParLimits(3,0.0,1.5*norm2);
  //fitfnc->SetParLimits(4,peak2-3*sigma2,peak2+3*sigma2);
  //fitfnc->SetParLimits(5,0.01,1.0);

  fitstatus = hrsp->Fit(fitfnc,"RQ");
  // if (0!=fitstatus){
  //   fitfnc->SetParLimits(4,0.2,1.7);
  //   fitfnc->SetParLimits(5,2.0,10.0);
  //   //cout <<" Not able to Fit this pt bin " << hrsp->GetName() << endl;
  // }

  fitstatus = hrsp->Fit(fitfnc,"RQ");
  hrsp->SetMaximum(norm1+0.2*norm1);
  fitg1 = new TF1("fg1","gaus(0)",fitrange_min,fitrange_max);
  fitg1->SetParameters(fitfnc->GetParameter(0),
		       fitfnc->GetParameter(1),
		       fitfnc->GetParameter(2));
  fitg1->SetLineColor(2);
  fitg1->SetLineStyle(2);
  hrsp->GetListOfFunctions()->Add(fitg1);

  fitg2 = new TF1("fg2","gaus(0)",fitrange_min,fitrange_max);
  fitg2->SetParameters(fitfnc->GetParameter(3),
		       fitfnc->GetParameter(4),
		       fitfnc->GetParameter(5));
  fitg2->SetLineColor(4);
  fitg2->SetLineStyle(4);
  hrsp->GetListOfFunctions()->Add(fitg2);

  if(hrsp->GetFunction("fdgaus")==0){
    cout << "No function recorded in histogram " << hrsp->GetName() << endl;
  }
  if (0!=fitstatus){
    cout<<"fit_double_gaussian() to "<<hrsp->GetName()
        <<" failed. Fitstatus: "<<fitstatus
        <<" - FNC deleted."<<endl;
    hrsp->GetListOfFunctions()->Delete();
  }
}

void adjust_fitrange(TH1* h,double& min,double& max)
{
  int imin=1; while (h->GetBinLowEdge(imin)<min) imin++;
  int imax=1; while (h->GetBinLowEdge(imax)<max) imax++;
  while ((imax-imin)<8) {
    if (imin>1) {imin--; min = h->GetBinCenter(imin); }
    if (imax<h->GetNbinsX()-1) { imax++; max=h->GetBinCenter(imax); }
  }
}
void adjust_fitrange_new(TH1* h,double& min,double& max)
{
  int imin=-1;
  for(int ix=1; ix<=h->GetNbinsX(); ix++){
    if( h->GetBinContent(ix) > 0 ){
      imin=ix;
      break;
    }
  }
  int imax=-1;
  double intg = h->Integral();
  double sumval=0;
  for(int ix=1; ix<=h->GetNbinsX(); ix++){
    sumval += h->GetBinContent(ix);
    if( (sumval / intg) >= 0.98){
      imax = ix;
      break;
    }
  }
  min = h->GetBinCenter(imin);
  max = h->GetBinCenter(imax);
}
double GetHistMax(TH1 *h)
{
  double max=-999;
  int ibin=-1;
  for(int ix=1; ix<=h->GetNbinsX(); ix++){
    if( h->GetBinContent(ix) > max ){
      max=h->GetBinContent(ix);
      ibin=ix;
    }
  }
  return (ibin > 0) ? h->GetBinCenter(ibin) : -1.0;
}
