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

#include "utilities.h"

using namespace std;

void fit_double_gaussian(TH1F*& hrsp);

/// default fit with gaussian in niter iteration of mean  
void fit_gaussian(TH1F*& hrsp,
                  const double nsigma,
                  const double jtptmin,
                  const int niter);

/// optional double sided crystal ball fit to response distributions
int fit_dscb(TH1F*& hrsp,
             const double nsigma,
             const double jtptmin,
             const int niter,
             const string alg);

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
double ptbins[] = {10,15,20,27,35,45,57,72,90,120,150,200,300,400,550,750,1000}; 
//double ptbins[]={10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15, 16, 17, 18, 19, 20, 22, 25, 27, 30, 32, 35, 37, 40, 42, 45, 47,
//  		 50, 52, 55, 57, 60, 62, 65, 67, 70, 72, 75, 77, 80, 85, 90, 95, 100,
//  		 120, 150, 200, 300, 400, 550, 750, 1000
//};
const int nbins = sizeof(ptbins)/sizeof(double) - 1;

const int knj = 7;
const char *calgo[knj]={"ak1","ak2","ak3","ak4","ak5","ak6","ak7"};

int cal_l3corr(const char *ctype="PF")
{
  LoadStyle();
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);

  cout << " # of bins in pt : "  << nbins << endl;
  //return 0;

  //TVirtualFitter::SetDefaultFitter("Minuit");
  TH1F *hrsp[knj][nbins];
  TGraphErrors *grsp[knj], *gcor[knj];
  TH1F *hMean[knj], *hRMS[knj];
  TH1F *hMean1[knj], *hSigma1[knj];
  TH1F *hMean2[knj], *hSigma2[knj];

  //TFile *fin = new TFile("jra_hiF_ak_dijet_full_finebin_lowpt.root","r");
  //TFile *fin = new TFile("jra_hiF_ak_dijet_final_lowpt.root","r");
  //TFile *fin = new TFile("jra_hiF_ak_dijet_test_lowpt.root","r");

  TFile *fin = new TFile("jra_hiF_ak_dijet_official_lowpt.root","r");

  double refptmin=0;
  double refptmax=1000;

  TLegend *l0 = getLegend(0.6349664,0.3265565,0.8491835,0.7128335);
  l0->SetHeader(Form("Anti-k_{T}, %s",ctype));
  TLegend *l1 = getLegend(0.1500646,0.5052258,0.3653091,0.7490942);
  l1->SetHeader("");

  TH1F *hrefpt=0, *hjetpt=0;
  //TSpectrum *spec = new TSpectrum(10);
  for(int i=0; i<knj; i++){

    grsp[i] = new TGraphErrors(nbins);
    gcor[i] = new TGraphErrors(nbins);

    grsp[i]->SetName(Form("L3RspVsRefPt_%d",i));
    grsp[i]->SetTitle(Form("L3RspVsRefPt %s",calgo[i]));
    grsp[i]->SetLineColor(1);
    grsp[i]->SetMarkerColor(1);
    grsp[i]->SetMarkerStyle(20);
    grsp[i]->SetMarkerSize(1.2);

    gcor[i]->SetName(Form("L3CorVsJetPt_%d",i));
    gcor[i]->SetTitle(Form("L3CorVsJetPt %s",calgo[i]));
    gcor[i]->SetLineColor(1);
    gcor[i]->SetMarkerColor(1);
    gcor[i]->SetMarkerStyle(20);
    gcor[i]->SetMarkerSize(1.2);


    hMean [i] = new TH1F(Form("hMean_%d",i),Form("ak%d%s",i+1,ctype),nbins,ptbins);
    hMean [i]->SetMarkerColor(1);
    hMean [i]->SetMarkerStyle(20);
    hMean [i]->SetLineColor(1);
    hMean [i]->GetXaxis()->SetMoreLogLabels();
    hMean [i]->GetXaxis()->SetNoExponent();

    hRMS [i] = new TH1F(Form("hRMS_%d",i),Form("ak%d%s",i+1,ctype),nbins,ptbins);
    hRMS [i]->SetMarkerColor(1);
    hRMS [i]->SetMarkerStyle(20);
    hRMS [i]->SetLineColor(1);

    hMean1 [i] = new TH1F(Form("hMean1_%d",i),Form("ak%d%s",i+1,ctype),nbins,ptbins);
    hMean1 [i]->SetMarkerColor(2);
    hMean1 [i]->SetMarkerStyle(20);
    hMean1 [i]->SetLineColor(2);

    hSigma1[i] = new TH1F(Form("hSigma1_%d",i),Form("ak%d%s",i+1,ctype),nbins,ptbins);
    hSigma1[i]->SetMarkerColor(2);
    hSigma1[i]->SetMarkerStyle(20);
    hSigma1[i]->SetLineColor(2);

    hMean2 [i] = new TH1F(Form("hMean2_%d",i),Form("ak%d%s",i+1,ctype),nbins,ptbins);
    hMean2 [i]->SetMarkerColor(4);
    hMean2 [i]->SetMarkerStyle(20);
    hMean2 [i]->SetLineColor(4);

    hSigma2[i] = new TH1F(Form("hSigma2_%d",i),Form("ak%d%s",i+1,ctype),nbins,ptbins);
    hSigma2[i]->SetMarkerColor(4);
    hSigma2[i]->SetMarkerStyle(20);
    hSigma2[i]->SetLineColor(4);

    std::ostringstream strs;
    if(strcmp(ctype,"PF")==0)strs << "ak" << i+1 << "PFJetAnalyzer/";
    else if(strcmp(ctype,"Calo"  )==0)strs << "ak"   << i+1 << "CaloJetAnalyzer/";
    else if(strcmp(ctype,"PuPF"  )==0)strs << "akPu" << i+1 << "PFJetAnalyzer/";
    else if(strcmp(ctype,"PuCalo")==0)strs << "akPu" << i+1 << "CaloJetAnalyzer/";
    else if(strcmp(ctype,"VsPF"  )==0)strs << "akVs" << i+1 << "PFJetAnalyzer/";
    else if(strcmp(ctype,"VsCalo")==0)strs << "akVs" << i+1 << "CaloJetAnalyzer/";
    
    for(int j=0; j<nbins; j++){
      
      //if(j<2)continue;
      std::ostringstream spt; 
      spt << ptbins[j] << "to" << ptbins[j+1];

      std::string sone  (strs.str()+"RefPt_Barrel_RefPt"+spt.str());
      std::string stwo  (strs.str()+"JetPt_Barrel_RefPt"+spt.str());
      std::string sthree(strs.str()+"RelRsp_Barrel_RefPt"+spt.str());

      // if(j==nbins-1){
      // 	cout << " hrefpt : "  << sone.c_str() << endl;
      // 	cout << " hjetpt : "  << stwo.c_str() << endl;
      // 	cout << " hrelsp : "  << sthree.c_str() << endl;
      // }
      hrefpt     = (TH1F*)fin->Get(Form("%s",sone.c_str()));
      hjetpt     = (TH1F*)fin->Get(Form("%s",stwo.c_str()));
      hrsp[i][j] = (TH1F*)fin->Get(Form("%s",sthree.c_str()));
      //hrsp[i][j]->Rebin(5);


      assert(hrefpt->GetEntries()>0 && hjetpt->GetEntries()>0);
      
      double refpt   =hrefpt->GetMean();
      double erefpt  =hrefpt->GetMeanError();
      double jetpt   =hjetpt->GetMean();
      double ejetpt  =hjetpt->GetMeanError();

      // if(i==2){
      // 	cout << spt.str().c_str() << setw(15) <<  "pt : "<< setw(10) << ptval << setw(10) << "refpt : "<< setw(10) << refpt << setw(15) << "jetpt : " << setw(10) << jetpt << endl;
      // }


      if(j==0){
	refptmin = refpt;
      }else if(j==nbins-2){
	refptmax = refpt;
      }

      hrsp[i][j]->SetLineColor(1);
      hrsp[i][j]->SetMarkerColor(1);
      hrsp[i][j]->SetMarkerStyle(20);
      //hrsp[i][j]->Scale(1./hrsp[i][j]->Integral());
      //hrsp[i][j]->SetMinimum(0);
      if(j==0)l0->AddEntry(hrsp[i][j],Form("R = %0.1f",(i+1)*0.1),"p");
      
      double peak  = hrsp[i][j]->GetMean();
      double epeak = hrsp[i][j]->GetMeanError();

      // if(strcmp(ctype,"PuPF")==0 && ptbins[j]>20){
      // 	int nfound = spec->Search(hrsp[i][j],1.2,"new");
      // 	//if(nfound==2 && i==2){
      // 	if(nfound>1){
      // 	  float *xpeaks = spec->GetPositionX();
      // 	  cout << calgo[i] << ctype << " " << j << "  " << spt.str().c_str() << "  peak1 : " << xpeaks[0] << " peak2 : " << xpeaks[1] << endl;
      // 	  for(int mp=0; mp<nfound; mp++){
      // 	    peak += xpeaks[mp];
      // 	  }
      // 	  peak /= nfound;
      // 	}
      // }


      //if(i==0){
      //fit_dscb(hrsp[i][j], 2.5, 1.0, 4, Form("ak%dPF",i+1));
      //fit_gaussian(hrsp[i][j], 1.5, 1.0, 4);
      //fit_double_gaussian(hrsp[i][j]);
      //fit_npeaks(hrsp[i][j]);
      
      //TF1*  frelrsp = (TF1*)hrsp[i][j]->GetListOfFunctions()->Last();
      //peak    = (frelrsp==0) ? hrsp[i][j]->GetMean()      :  frelrsp->GetParameter(1);
      //epeak   = (frelrsp==0) ? hrsp[i][j]->GetMeanError() :  frelrsp->GetParError(1);
      //peak    = (frelrsp==0) ? hrsp[i][j]->GetMean()      : ((frelrsp->GetParameter(1)+hrsp[i][j]->GetMean())*0.5);
      //epeak   = (frelrsp==0) ? hrsp[i][j]->GetMeanError() : 0.5*sqrt(pow(hrsp[i][j]->GetMeanError(),2)+pow(frelrsp->GetParError(1),2));
      //}

      double cor  = 1.0/peak;
      double ecor = cor*cor*epeak;
      
      grsp[i]->SetPoint     (j, refpt, peak);
      grsp[i]->SetPointError(j, erefpt, epeak);
      gcor[i]->SetPoint     (j, jetpt, cor);
      gcor[i]->SetPointError(j, ejetpt, ecor);

      // grsp[i]->SetPoint     (j, ptval, peak);
      // grsp[i]->SetPointError(j, 0    , epeak);
      // gcor[i]->SetPoint     (j, ptval, cor);
      // gcor[i]->SetPointError(j, 0, ecor);
      
    }
  }
  //hrsp[2][15]->Draw("p");
  //int nfound = spec->Search(hrsp[2][15],1.0,"new",0.01);
  //cout << " $$$$$$$$$$ " << nfound << endl;
  //return 0;

  int ipad=0;
  TCanvas *c99[knj];
  for(int nj=0;nj<knj;nj++){
    //c99[nj] = new TCanvas(Form("c99_%d",nj),Form("%s Fitting plots",calgo[nj]),100,102,1399,942);
    c99[nj] = new TCanvas(Form("c99_%d",nj),Form("%s Fitting plots",calgo[nj]),1854,866);
    //c99[nj]->Divide(10,5,0,0);
    c99[nj]->Divide(4,4,0,0);
    ipad=0;
    for(int ip=0;ip<nbins;ip++){      
      c99[nj]->cd(++ipad);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.01);
      gPad->SetBottomMargin(0.15);
      //if(ipad%10==0)gPad->SetRightMargin(0.02);
      //if(ipad%(ip+1)==0)gPad->SetLeftMargin(0.15);

      
      hrsp[nj][ip]->SetMaximum((hrsp[nj][ip]->GetMaximum() + 0.25*hrsp[nj][ip]->GetMaximum()));
      hrsp[nj][ip]->SetMinimum(0);
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

      TF1 *fdg = (TF1*)hrsp[nj][ip]->GetFunction("fdgaus");
      if(fdg){
	TF1 *fg1 = (TF1*)hrsp[nj][ip]->GetFunction("fg1");      
	TF1 *fg2 = (TF1*)hrsp[nj][ip]->GetFunction("fg2");      
	
	fdg->Draw("lsame");
	fg1->Draw("lsame");
	fg2->Draw("lsame");
      }

      std::ostringstream strs; 
      strs << ptbins[ip] << "< p_{T} (GeV/c) <" << ptbins[ip+1];
      std::string spt = strs.str();
      
      if(ipad==1){drawText2(Form("%s%s",calgo[nj],ctype),0.28,0.90,18);      
	
	drawText2(spt.c_str(),0.22,0.80,15);		
      } else drawText2(spt.c_str(),0.17,0.80,15);		
    }
    //if(nj!=0)c99[nj]->Close();
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

  if(strcmp(ctype,"PF")==0 || strcmp(ctype,"PuPF")==0 || strcmp(ctype,"VsPF")==0){
    fcn_rsp="[0]-[1]/(pow(log10(x)-[2],2)+[2])-[3]*exp((-[4]*(log10(x)-[5])*(log10(x)-[5]))+([6]*(log10(x)-[5])))";
    fcn_cor="[0]+[1]/(pow(log10(x)-[2],2)+[2])+[3]*exp((-[4]*(log10(x)-[5])*(log10(x)-[5]))+([6]*(log10(x)-[5])))";
  }else if(strcmp(ctype,"Calo")==0 || strcmp(ctype,"PuCalo")==0 || strcmp(ctype,"VsCalo")==0){
    fcn_rsp="[0]-[1]/(pow(log10(x),[2])+[3])+[4]/x";
    fcn_cor="[0]+[1]/(pow(log10(x),[2])+[3])-[4]/x";
  }

  //double xmin[knj]={10,5e+03,18,18,35,35,45};
  //double xmax[knj]={27,5e+04,40,40,72,72,90};
  TCanvas *c11 = new TCanvas("c11","L3 Response",1531,683);
  c11->Divide(4,2);
  TCanvas *c12 = new TCanvas("c12","L3 Corrections",1531,683);
  c12->Divide(4,2);
  for(int i=0; i<knj; i++){  
    c11->cd(i+1);

    fitrsp = new TF1("fitrsp",fcn_rsp.c_str(),1.0,grsp[i]->GetX()[grsp[i]->GetN()-1]);

    if(strcmp(ctype,"PF")==0){

      //! PbPb
      // fitrsp->SetParameter(0,0.98);
      // fitrsp->SetParameter(1,0.033);
      // fitrsp->SetParameter(2,-0.70);
      // fitrsp->SetParameter(3,0.020);
      // fitrsp->SetParameter(4,1.020);
      // fitrsp->SetParameter(5,2.10);
      // fitrsp->SetParameter(6,0.013);

      //! pp
      fitrsp->SetParameter(0,0.96);
      fitrsp->SetParameter(1,0.033);
      fitrsp->SetParameter(2,-0.70);
      fitrsp->SetParameter(3,0.020);
      fitrsp->SetParameter(4,1.020);
      fitrsp->SetParameter(5,2.10);
      fitrsp->SetParameter(6,0.013);

    }
    else if(strcmp(ctype,"PuPF")==0){
      fitrsp->SetParameter(0,0.98);
      fitrsp->SetParameter(1,0.044);
      fitrsp->SetParameter(2,-0.70);
      fitrsp->SetParameter(3,0.020);
      fitrsp->SetParameter(4,1.020);
      fitrsp->SetParameter(5,2.10);
      fitrsp->SetParameter(6,0.025);

     }
    else if(strcmp(ctype,"VsPF")==0){
      fitrsp->SetParameter(0,0.96);
      fitrsp->SetParameter(1,0.033);
      fitrsp->SetParameter(2,-0.70);
      fitrsp->SetParameter(3,0.020);
      fitrsp->SetParameter(4,1.020);
      fitrsp->SetParameter(5,2.10);
      fitrsp->SetParameter(6,0.013);
    }
    else if(strcmp(ctype,"Calo")==0 || strcmp(ctype,"PuCalo")==0 || strcmp(ctype,"VsCalo")==0){

      fitrsp->SetParameter(0,1.0);
      fitrsp->SetParameter(1,1.0);
      fitrsp->SetParameter(2,1.0);
      fitrsp->SetParameter(3,1.0);
      fitrsp->SetParameter(4,1.0);

    }
    
    //if(strcmp(ctype,"PuPF")==0)grsp[i]->Fit(fitrsp,"R","",20,grsp[i]->GetX()[gcor[i]->GetN()-1]);
    //else grsp[i]->Fit(fitrsp,"R");

    grsp[i]->Fit(fitrsp,"R");
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

    //grsp[i]->GetXaxis()->SetRangeUser(15,600);
    drawText2(Form("%s%s",calgo[i],ctype),0.20,0.78,21);

    TPaveStats *strsp = (TPaveStats*)grsp[i]->FindObject("stats");
    strsp->SetX1NDC(0.46);
    strsp->SetY1NDC(0.18);
    strsp->SetX2NDC(0.87);
    strsp->SetY2NDC(0.59);
    strsp->Draw();
    delete fitrsp;
    //continue;
    
    c12->cd(i+1);
    //cout << endl;
    //cout << endl;
    
    // cout << "Fitting the L3 Corrections "  << endl;
    // cout << endl;
  
    fitcor = new TF1("fitcor",fcn_cor.c_str(),
		     1.0,gcor[i]->GetX()[gcor[i]->GetN()-1]);

    if(strcmp(ctype,"PF")==0){

      //! PbPb
      // fitcor->SetParameter(0,1.02);
      // fitcor->SetParameter(1,0.033);
      // fitcor->SetParameter(2,-0.70);
      // fitcor->SetParameter(3,0.020);
      // fitcor->SetParameter(4,1.020);
      // fitcor->SetParameter(5,2.10);
      // fitcor->SetParameter(6,0.013);
      
      //! pp
      fitcor->SetParameter(0,1.04);
      fitcor->SetParameter(1,0.033);
      fitcor->SetParameter(2,-0.70);
      fitcor->SetParameter(3,0.020);
      fitcor->SetParameter(4,1.020);
      fitcor->SetParameter(5,2.10);
      fitcor->SetParameter(6,0.013);
  
    }else if(strcmp(ctype,"PuPF")==0){
      fitcor->SetParameter(0,1.04);
      fitcor->SetParameter(1,0.044);
      fitcor->SetParameter(2,-0.7);
      fitcor->SetParameter(3,0.02);
      fitcor->SetParameter(4,1.02);
      fitcor->SetParameter(5,2.10);
      fitcor->SetParameter(6,0.025);

    }else if(strcmp(ctype,"VsPF")==0){
      fitcor->SetParameter(0,1.04);
      fitcor->SetParameter(1,0.033);
      fitcor->SetParameter(2,-0.70);
      fitcor->SetParameter(3,0.020);
      fitcor->SetParameter(4,1.020);
      fitcor->SetParameter(5,2.10);
      fitcor->SetParameter(6,0.013);

    }else if(strcmp(ctype,"Calo")==0 || strcmp(ctype,"PuCalo")==0 || strcmp(ctype,"VsCalo")==0){
      fitcor->SetParameter(0,1.0);
      fitcor->SetParameter(1,1.0);
      fitcor->SetParameter(2,1.0);
      fitcor->SetParameter(3,1.0);
      fitcor->SetParameter(4,1.0);

    }

    //if(strcmp(ctype,"PuPF")==0)gcor[i]->Fit(fitcor,"R","",20,gcor[i]->GetX()[gcor[i]->GetN()-1]);
    //else gcor[i]->Fit(fitcor,"R");

    gcor[i]->Fit(fitcor,"R");

    // gcor[i]->GetXaxis()->SetRangeUser(15,600);
    gcor[i]->SetMaximum(2.45);
    gcor[i]->SetMinimum(0.89);
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
    delete fitcor;
  }
  //c12->Close();

  // TFile *fo = new TFile("test_ak3PF.root","RECREATE");
  // fo->cd();

  // for(int i=13; i<22; i++){
  //   for(int j=13; j<22; i++){
  //   hrsp[2][i]->Write();
  // }
  // fo->Close();
  // return 0;


  //int isel=4;
  // new TCanvas();
  // grsp[isel]->SetMaximum(1.2);
  // grsp[isel]->SetMinimum(0.2);
  // grsp[isel]->Draw("ap");
  // gPad->SetLogx();
  return 0;
}
int fit_dscb(TH1F*& hrsp,
             const double nsigma,
             const double jtptmin,
             const int niter,
             const string alg)
{
  if (0==hrsp) {
    cout<<"ERROR: Empty pointer to fit_dscb()"<<endl;return -1;
  }


  // first use a gaussian to constrain crystal ball gaussian core
  fit_gaussian(hrsp,nsigma,jtptmin,niter);
  TF1* fgaus = hrsp->GetFunction("fgaus");
  if (0==fgaus) {
    hrsp->GetListOfFunctions()->Delete();
    return -1;
  }

  // implementation of the low pt bias threshold 
  string histname = hrsp->GetName();
  //double ptRefMax(1.0),rspMax(0.0);
  double fitrange_min(0.2);
  double fitrange_max(1.7);

  TF1* fdscb = new TF1("fdscb",fnc_dscb,fitrange_min,fitrange_max,7);
  fdscb->SetLineWidth(2);
  fdscb->SetLineStyle(2);

  double norm = fgaus->GetParameter(0);
  double mean = fgaus->GetParameter(1);
  double sigma= fgaus->GetParameter(2);

  double aone(2.0),atwo(2.0),pone(10.0),ptwo(10.0);
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
    
    if (i>0) fdscb->FixParameter(3,aone);
    else fdscb->SetParLimits(3,1.,50.);
    
    if (i>1) fdscb->FixParameter(5,atwo);
    else fdscb->SetParLimits(5,1.,50.);

    fdscb->SetParLimits(4,0.,100.);
    fdscb->SetParLimits(6,0.,100.);

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

    // reset sigma and mean to gauss values...
    fdscb->SetParameter(1,fgaus->GetParameter(1));
    fdscb->SetParError (1,fgaus->GetParError(1));
    fdscb->SetParameter(2,fgaus->GetParameter(2));
    fdscb->SetParError (2,fgaus->GetParError(2));
  }
  if (0!=fitstatus){
    cout<<"fit_fdscb() to "<<hrsp->GetName()
        <<" failed. Fitstatus: "<<fitstatus<<endl;
    hrsp->GetFunction("fdscb")->Delete();
  }
  return fitstatus;
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

void fit_gaussian(TH1F*& hrsp,
                  const double nsigma,
                  const double jtptmin,
                  const int niter)
{
  if (0==hrsp) {
    cout<<"ERROR: Empty pointer to fit_gaussian()"<<endl;return;
  }
  
  string histname = hrsp->GetName();
  double mean     = hrsp->GetMean();
  double rms      = hrsp->GetRMS();
  double ptRefMax(1.0),rspMax(0.0);

  //cout << " Mean : "  << mean << " \t  RMS  : " << rms << endl;

  double norm  = hrsp->GetMaximumStored();
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
    fitrange_min = 0.4;
    fitrange_max = 1.2;
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
  //   if (0==fitstatus){
  //     cout<<"fit_gaussian() to "<<hrsp->GetName()    <<"  sucessful : " <<endl;
  //   }          

  if (0!=fitstatus){
    cout<<"fit_gaussian() to "<<hrsp->GetName()
        <<" failed. Fitstatus: "<<fitstatus
        <<" - FNC deleted."<<endl;
    hrsp->GetListOfFunctions()->Delete();
  }
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

  float *xpeaks = spec->GetPositionX();

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
