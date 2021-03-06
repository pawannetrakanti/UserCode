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
TLegend *getLegend(double /*x1*/, double /*y1*/, double /*x2*/, double /*y2*/);

template <class T>
bool from_string(T& t,
                 const std::string& s,
                 std::ios_base& (*f)(std::ios_base&))
{
  std::istringstream iss(s);
  return !(iss >> f >> t).fail();
}
//double ptbins[] ={12, 17, 22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 400, 550, 790, 1000};
//double ptbins[] ={22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 350, 400, 550, 790, 1000};
//double ptbins[] ={17, 22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 350, 400, 550, 790, 1000};

double ptbins[] ={12,  12.5,  13,  13.5,  14,  15,  17,  20,  23,  27,  30,  35,  40,  45,  57,  
		  72,  90,  120,  150,  200,  300,  400,  600,  1000
};


const int nbins = sizeof(ptbins)/sizeof(double) - 1;

const int knj = 6;
const char *calgo[knj]={"ak1","ak2","ak3","ak4","ak5","ak6"};

int PFL3(const char *ctype="PF")
{

  bool iSave=false;
  LoadStyle();
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);

  double minx=1.0;
  string outdir="txtfiles/";
  //string outdir="junk/";
  cout << " # of bins in pt : "  << nbins << endl;
  //return 0;

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000); 
  TVirtualFitter::SetDefaultFitter("Minuit");
  TH1F *hrsp[knj][nbins];
  TGraphErrors *grsp[knj], *gcor[knj];
  TH1F *hMean[knj], *hRMS[knj];

  //TFile *fin = new TFile("jra_hiF_akpf_dijet_HIreco.root","r");
  //TFile *fin = new TFile("jra_hiF_akpf_dijet_ppreco_nocut.root","r");

  //TFile *fin = new TFile("jra_hiF_akpf_dijet_hireco.root","r");
  //TFile *fin = new TFile("jra_hiF_akpf_ppReco_new_dijet.root","r");

  //TFile *fin = new TFile("jra_hiF_akpf_ppReco_757p1_dijet.root","r");

  //TFile *fin = new TFile("jra_hiF_akpf_ppReco_757p1_dijet_fulleta_fullpt.root","r");

  TFile *fin = new TFile("jra_hiF_akpf_ppReco_757p1_HcalRespCorrs_v4_00_mc_dijet.root","r");

  double refptmin=0;
  double refptmax=1000;

  TLegend *l0 = getLegend(0.6349664,0.3265565,0.8491835,0.7128335);
  l0->SetHeader(Form("Anti-k_{T}, %s",ctype));
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

    hMean [i] = new TH1F(Form("hMean_%d",i),Form("%s%s",calgo[i],ctype),nbins,ptbins);
    hMean [i]->SetMarkerColor(1);
    hMean [i]->SetMarkerStyle(20);
    hMean [i]->SetLineColor(1);
    hMean [i]->GetXaxis()->SetMoreLogLabels();
    hMean [i]->GetXaxis()->SetNoExponent();

    hRMS [i] = new TH1F(Form("hRMS_%d",i),Form("%s%s",calgo[i],ctype),nbins,ptbins);
    hRMS [i]->SetMarkerColor(1);
    hRMS [i]->SetMarkerStyle(20);
    hRMS [i]->SetLineColor(1);


    std::ostringstream strs;
     strs << calgo[i] << "PFJetAnalyzer/";
     cout << " Running for " << strs.str().c_str() << endl;

    // if(strcmp(ctype,"PF")==0)strs << "ak" << i+1 << "PFJetAnalyzer/";
    // else if(strcmp(ctype,"Calo"  )==0)strs << "ak"   << i+1 << "CaloJetAnalyzer/";
    // else if(strcmp(ctype,"PuPF"  )==0)strs << "akPu" << i+1 << "PFJetAnalyzer/";
    // else if(strcmp(ctype,"PuCalo")==0)strs << "akPu" << i+1 << "CaloJetAnalyzer/";
    // else if(strcmp(ctype,"VsPF"  )==0)strs << "akVs" << i+1 << "PFJetAnalyzer/";
    // else if(strcmp(ctype,"VsCalo")==0)strs << "akVs" << i+1 << "CaloJetAnalyzer/";

    // if( njmin !=0 && njmax != knj ){
    //   cout << " Running for " << strs.str().c_str() << endl;
    //   cout << endl;
    // }

    
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

      // if(i==0){
      // 	//fit_dscb(hrsp[i][j], 1.75, 1.0, 5, Form("ak%d%s",i+1,ctype), 1);      
      // 	if( j < 8)fit_dscb(hrsp[i][j], 0.80, 1.0, 5, Form("ak%d%s",i+1,ctype), 1);      
      // 	else fit_dscb(hrsp[i][j], 1.00, 1.0, 5, Form("ak%d%s",i+1,ctype), 2);      

      // }else fit_dscb(hrsp[i][j], 1.00, 1.0, 8, Form("ak%d%s",i+1,ctype), 0);      

      
      //fit_dscb(hrsp[i][j], 1.0, minx, 4, Form("%s%s",calgo[i],ctype), 0);      

      
      // if( i<2 )fit_dscb(hrsp[i][j], 0.95, minx, 5, Form("ak%d%s",i+1,ctype),0);      
      // else if(i==2)fit_dscb(hrsp[i][j], 1.45, minx, 5, Form("ak%d%s",i+1,ctype),0);      
      // else if(i==3 || i==4)fit_dscb(hrsp[i][j], 1.50, minx, 4, Form("ak%d%s",i+1,ctype),0);
      // else fit_dscb(hrsp[i][j], 1.30, minx, 5, Form("ak%d%s",i+1,ctype),0);

      //! pp reco
      // if( i == 0 )fit_dscb(hrsp[i][j], 1.25, minx, 4, Form("ak%d%s",i+1,ctype),0); //! ak1 
      // else if( i == 1 )fit_dscb(hrsp[i][j], 1.00, minx, 5, Form("%s%s",calgo[i],ctype),0); //! ak2 
      // else if( i == 2 )fit_dscb(hrsp[i][j], 1.04, minx, 5, Form("%s%s",calgo[i],ctype),0); //! ak3 
      // else if( i == 3 )fit_dscb(hrsp[i][j], 1.00, minx, 5, Form("%s%s",calgo[i],ctype),0); //! ak4
      // else if( i == 4 )fit_dscb(hrsp[i][j], 1.00, minx, 5, Form("%s%s",calgo[i],ctype),0); //! ak5
      // else if( i == 5 )fit_dscb(hrsp[i][j], 0.97, minx, 5, Form("%s%s",calgo[i],ctype),0); //! ak6


      //! HI reco new
      // if( i == 5 )fit_dscb(hrsp[i][j], 0.97, minx, 4, Form("ak%d%s",i+1,ctype),0); //! ak6
      // else if( i == 4 )fit_dscb(hrsp[i][j], 1.00, minx, 5, Form("ak%d%s",i+1,ctype),0); //! ak5
      // else if( i == 3 )fit_dscb(hrsp[i][j], 1.00, minx, 5, Form("ak%d%s",i+1,ctype),0); //! ak4
      // else if( i == 2 )fit_dscb(hrsp[i][j], 1.04, minx, 5, Form("ak%d%s",i+1,ctype),0); //! ak3 
      // else if( i == 1 )fit_dscb(hrsp[i][j], 0.97, minx, 5, Form("ak%d%s",i+1,ctype),0); //! ak2 
      // else if( i == 0 )fit_dscb(hrsp[i][j], 1.02, minx, 5, Form("ak%d%s",i+1,ctype),0); //! ak2 


      // if(i!=0)fit_dscb(hrsp[i][j], 1.15, minx, 4, Form("ak%d%s",i+1,ctype), 0);      
      // if(i==2)fit_dscb(hrsp[i][j], 0.95, minx, 5, Form("ak%d%s",i+1,ctype), 0);      
      // else fit_dscb(hrsp[i][j], 1.00, minx, 5, Form("ak%d%s",i+1,ctype), 0);      

      //fit_dscb(hrsp[i][j], 1.00, minx, 5, Form("ak%d%s",i+1,ctype), 0);      

      TF1*  frelrsp = (TF1*)hrsp[i][j]->GetListOfFunctions()->Last();
      peak    = (frelrsp==0) ? hrsp[i][j]->GetMean()      :  frelrsp->GetParameter(1);
      epeak   = (frelrsp==0) ? hrsp[i][j]->GetMeanError() :  frelrsp->GetParError(1);
      peak    = (frelrsp==0) ? hrsp[i][j]->GetMean()      :  frelrsp->GetParameter(1);
      epeak   = (frelrsp==0) ? hrsp[i][j]->GetMeanError() :  frelrsp->GetParError(1);


      double cor  = 1.0/peak;
      double ecor = cor*cor*epeak;
      
      grsp[i]->SetPoint     (j, refpt,  peak);
      grsp[i]->SetPointError(j, erefpt, epeak);
      gcor[i]->SetPoint     (j, jetpt,  cor);
      gcor[i]->SetPointError(j, ejetpt, ecor);
    }
  }
  //hrsp[2][15]->Draw("p");
  //return 0;
  if(gPad)gPad->Close();

  int ipad=0;
  //! 0 - 20 GeV
  TCanvas *c98[knj];
  for(int nj=njmin;nj<njmax;nj++){
    c98[nj] = new TCanvas(Form("c99_%d",nj),Form("Fine %s%s Fitting plots",calgo[nj],ctype),100,102,1399,942);
    //c98[nj]->Divide(6,4,0,0);
    c98[nj]->Divide(7,4,0,0);
    ipad=0;
    for(int ip=0;ip<nbins;ip++){      
      c98[nj]->cd(++ipad);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.01);
      gPad->SetBottomMargin(0.15);
      //if(ipad%10==0)gPad->SetRightMargin(0.02);
      //if(ipad%(ip+1)==0)gPad->SetLeftMargin(0.15);

      gPad->SetLogy();
      
      // hrsp[nj][ip]->SetMaximum((hrsp[nj][ip]->GetMaximum() + 0.25*hrsp[nj][ip]->GetMaximum()));
      // hrsp[nj][ip]->SetMinimum(1e-12);
      hrsp[nj][ip]->SetMaximum(hrsp[nj][ip]->Integral()*2e-00);
      hrsp[nj][ip]->SetMinimum(hrsp[nj][ip]->Integral()*1e-07);
      hrsp[nj][ip]->SetTitle(0);
      hrsp[nj][ip]->GetXaxis()->SetTitle("<reco jet p_{T} / gen jet p_{T}>");
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
      // TF1 *fdg = (TF1*)hrsp[nj][ip]->GetFunction("fdgaus");
      // if(fdg){
      // 	TF1 *fg1 = (TF1*)hrsp[nj][ip]->GetFunction("fg1");      
      // 	TF1 *fg2 = (TF1*)hrsp[nj][ip]->GetFunction("fg2");      
      // 	fdg->Draw("lsame");
      // 	fg1->Draw("lsame");
      // 	fg2->Draw("lsame");
      // }
      
      std::ostringstream strs; 
      strs << ptbins[ip] << "< p_{T} (GeV/c) <" << ptbins[ip+1];
      std::string spt = strs.str();
      
      if(ipad==1){drawText2(Form("%s%s",calgo[nj],ctype),0.28,0.90,18);      
	
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
  std::string fcn_rsp="";
  std::string fcn_cor="";

  //fcn_rsp="[0]-[1]/(pow(log10(x)-[2],2)+[2])-[3]*exp((-[4]*(log10(x)-[5])*(log10(x)-[5]))+([6]*(log10(x)-[5])))";
  //fcn_cor="[0]+[1]/(pow(log10(x)-[2],2)+[2])+[3]*exp((-[4]*(log10(x)-[5])*(log10(x)-[5]))+([6]*(log10(x)-[5])))";

  //! Used this one
  fcn_rsp="[0]-[1]/(pow(log10(x)-[2],2)+[2])-[3]*exp((-[4]*(log10(x)-[5])*(log10(x)-[5]))+([6]*(log10(x)-[5])))";
  fcn_cor="[0]+[1]/(pow(log10(x)-[2],2)+[2])+[3]*exp((-[4]*(log10(x)-[5])*(log10(x)-[5]))+([6]*(log10(x)-[5])))";

  //double xmin[knj]={10,5e+03,18,18,35,35,45};
  //double xmax[knj]={27,5e+04,40,40,72,72,90};
  TCanvas *c11 = new TCanvas("c11","L3 Response",1234,730/*1531,683*/);
  c11->Divide(3,2);
  TCanvas *c12 = new TCanvas("c12","L3 Corrections",1234,730/*1531,683*/);
  c12->Divide(3,2);
  for(int i=njmin; i<njmax; i++){  
    c11->cd(i+1);

    fitrsp = new TF1("fitrsp",fcn_rsp.c_str(),1.0,grsp[i]->GetX()[grsp[i]->GetN()-1]);

    fitrsp->SetParameter(0,0.96);
    fitrsp->SetParameter(1,0.033);
    fitrsp->SetParameter(2,-0.70);
    fitrsp->SetParameter(3,0.020);
    fitrsp->SetParameter(4,1.020);
    fitrsp->SetParameter(5,2.10);
    fitrsp->SetParameter(6,0.013);


    // fitrsp->SetParameter(0,0.96);
    // fitrsp->SetParameter(1,0.033);
    // fitrsp->SetParameter(2,-0.70);
    // fitrsp->SetParameter(3,0.020);
    // fitrsp->SetParameter(4,1.020);
    // fitrsp->SetParameter(5,1.90);
    // fitrsp->SetParameter(6,0.010);

    grsp[i]->Fit(fitrsp,"QR");
    grsp[i]->SetMaximum(1.2);
    grsp[i]->SetMinimum(0.2);
    grsp[i]->Draw("ap");
    gPad->SetLogx();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->Update();

    grsp[i]->GetXaxis()->SetTitle("< Ref p_{T} > (GeV/c)");
    grsp[i]->GetXaxis()->SetTitleFont(42);
    grsp[i]->GetXaxis()->SetLabelFont(42);
    grsp[i]->GetXaxis()->SetLabelSize(0.06);
    grsp[i]->GetXaxis()->SetTitleSize(0.06);
    grsp[i]->GetXaxis()->SetNdivisions(507);
    grsp[i]->GetXaxis()->SetNoExponent();
    grsp[i]->GetXaxis()->SetMoreLogLabels();
    grsp[i]->GetYaxis()->SetTitle("L3 Response");
    grsp[i]->GetYaxis()->SetTitleFont(42);
    grsp[i]->GetYaxis()->SetLabelFont(42);
    grsp[i]->GetYaxis()->SetLabelSize(0.06);
    grsp[i]->GetYaxis()->SetTitleSize(0.06);
    grsp[i]->GetYaxis()->SetNdivisions(507);

    grsp[i]->GetYaxis()->SetRangeUser(0,1.160);
    //grsp[i]->GetXaxis()->SetRangeUser(15,600);
    drawText2(Form("%s%s",calgo[i],ctype),0.20,0.78,21);

    TPaveStats *strsp = (TPaveStats*)grsp[i]->FindObject("stats");
    strsp->SetX1NDC(0.46);
    strsp->SetY1NDC(0.18);
    strsp->SetX2NDC(0.87);
    strsp->SetY2NDC(0.59);
    strsp->Draw();
    //delete fitrsp;
    //continue;
    
    c12->cd(i+1);
    //cout << endl;
    //cout << endl;
    
    // cout << "Fitting the L3 Corrections "  << endl;
    // cout << endl;
  
    fitcor = new TF1("fitcor",fcn_cor.c_str(),
		     1.0,gcor[i]->GetX()[gcor[i]->GetN()-1]);

    fitcor->SetParameter(0,1.04);
    fitcor->SetParameter(1,0.033);
    fitcor->SetParameter(2,-0.70);
    fitcor->SetParameter(3,0.020);
    fitcor->SetParameter(4,1.020);
    fitcor->SetParameter(5,2.10);
    fitcor->SetParameter(6,0.013);

    // fitcor->SetParameter(0,1.04);
    // fitcor->SetParameter(1,0.033);
    // fitcor->SetParameter(2,-0.70);
    // fitcor->SetParameter(3,0.020);
    // fitcor->SetParameter(4,1.020);
    // fitcor->SetParameter(5,1.90);
    // fitcor->SetParameter(6,0.010);

    gcor[i]->Fit(fitcor,"QR");

    string algn=Form("%d",i+1);
    string era="JEC_pp_PYTHIA_TuneCUETP8M1_5020GeV_patch3";
    string txtfilename = outdir + era+"_L3Absolute_AK"+algn+"PF.txt";
    ofstream outf(txtfilename.c_str());
    outf.setf(ios::left);

    //cout<<"{1 JetEta 1 JetPt "<<fcn_cor<<" Correction L3Absolute}"<<endl;
    outf<<"{1 JetEta 1 JetPt "<<fcn_cor<<" Correction L3Absolute}"<<endl;

    // cout<<-3.000                // eta_min
    // 	<<setw(12)<<+3.000                // eta_max
    // 	<<setw(12)<<fitcor->GetNpar()+2   // number of parameters + 2
    // 	<<setw(12)<<4.0                   // minimum pT
    // 	<<setw(12)<<1000.0;               // maximum pT

    outf<<setw(12)<<-3.000                // eta_min
     	<<setw(12)<<+3.000                // eta_max
     	<<setw(12)<<fitcor->GetNpar()+2   // number of parameters + 2
     	<<setw(12)<<4.0                   // minimum pT
     	<<setw(12)<<1000.0;               // maximum pT


    for(int p=0; p<fitcor->GetNpar(); p++){
      //cout<<setw(12)<<fitcor->GetParameter(p); // p0-p4
      outf<<setw(12)<<fitcor->GetParameter(p); // p0-p4
    }
    outf.close();
    //cout<<endl;

    // if(i==0){
    //   double *rawpt = gcor[i]->GetX();
    //   for(int ix=0; ix<gcor[i]->GetN(); ix++){
    // 	double cor_fac = fitcor->Eval(rawpt[ix]);
    // 	double corrpt  =  cor_fac*rawpt[ix];
    // 	std::cout << " rawpt : " << rawpt[ix] << " corrpt : "  << corrpt << " cor_fac : "  << cor_fac << std::endl;
    //   }
    // }

    // gcor[i]->GetXaxis()->SetRangeUser(15,600);
    gcor[i]->SetMaximum(1.35);
    gcor[i]->SetMinimum(0.98);
    gcor[i]->Draw("ap");
    gPad->SetLogx();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->Update();

    gcor[i]->GetXaxis()->SetTitle("< raw jet p_{T} > (GeV/c)");
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
    drawText2(Form("%s%s",calgo[i],ctype),0.20,0.78,21);


    TPaveStats *stcor = (TPaveStats*)gcor[i]->FindObject("stats");
    stcor->SetX1NDC(0.48);
    stcor->SetY1NDC(0.53);
    stcor->SetX2NDC(0.89);
    stcor->SetY2NDC(0.89);
    stcor->Draw();
  
    // cout << endl;
    //delete fitcor;
  }

  if(iSave){
    c12->SaveAs("CorrectionPlots/L3Absolute_Corrections_ppReco_PFJets_757p1.pdf");
  }

  TFile *fout = new TFile("l3pf_input.root","RECREATE");
  for(int nj=njmin;nj<njmax;nj++){
    fout->mkdir(Form("%sPFJetAnalyzer",calgo[nj]),Form("%sPFJetAnalyzer",calgo[nj]));
    fout->cd(Form("%sPFJetAnalyzer",calgo[nj]));
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
  double fitrange_max(2.00);

  //adjust_fitrange(hrsp,fitrange_min,fitrange_max);
  
  TF1* fdscb = new TF1("fdscb",fnc_dscb,fitrange_min,fitrange_max,7);
  fdscb->SetLineWidth(2);
  fdscb->SetLineStyle(2);

  double norm = 2.0*fgaus->GetParameter(0);
  double mean = fgaus->GetParameter(1);
  double sigma= fgaus->GetParameter(2);

  // double norm = 2*hrsp->GetMaximumStored();
  // double mean = hrsp->GetMean();
  // double sigma= hrsp->GetRMS();

  double aone(2.0),atwo(2.0),pone(10.0),ptwo(10.0);
  //double aone(4.0),atwo(4.0),pone(15.0),ptwo(15.0);

  //double aone(2.0),atwo(2.0),pone(5.0),ptwo(5.0);
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


    // if( alg == "ak1PF" ){
    //   //! pp reco with no cut
    //   if (i>0) fdscb->FixParameter(3,aone);
    //   else fdscb->SetParLimits(3,0.,10.);
      
    //   if (i>1) fdscb->FixParameter(5,atwo);
    //   else fdscb->SetParLimits(5,0.,10.);
      
    //   fdscb->SetParLimits(4,0.,30.);
    //   fdscb->SetParLimits(6,0.,30.);
    // }else if( alg == "ak3PF" ){
    //   if (i>0) fdscb->FixParameter(3,aone);
    //   else fdscb->SetParLimits(3,0.,20.);
      
    //   if (i>1) fdscb->FixParameter(5,atwo);
    //   else fdscb->SetParLimits(5,0.,20.);
      
    //   fdscb->SetParLimits(4,0.,50.);
    //   fdscb->SetParLimits(6,0.,50.);
    // }
    // else{
    //   //! pp reco with no cut
    //   if (i>0) fdscb->FixParameter(3,aone);
    //   else fdscb->SetParLimits(3,0.,40.);
      
    //   if (i>1) fdscb->FixParameter(5,atwo);
    //   else fdscb->SetParLimits(5,0.,40.);
      
    //   fdscb->SetParLimits(4,0.,80.);
    //   fdscb->SetParLimits(6,0.,80.);
    // }


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
  double mean     = hrsp->GetMean();
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
    double fitrange_min = *std::max_element(vv.begin(),vv.end());
    double fitrange_max = std::min(xmax,peak+nsigma*sigma);

    // if( stat ){
    //   fitrange_min = 0.90;
    //   fitrange_max = 1.10;
    // }
    // else adjust_fitrange(hrsp,fitrange_min,fitrange_max);

    adjust_fitrange(hrsp,fitrange_min,fitrange_max);

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
TLegend *getLegend(double x1, double y1, double x2, double y2)
{
  TLegend *leg = new TLegend(x1,y1,x2,y2,NULL,"BRNDC");
  leg->SetHeader("");
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(10);
  leg->SetFillStyle(1001);
  return leg;
}
