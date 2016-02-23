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
#include "TCanvas.h"
#include <TMarker.h>
#include <TString.h>
#include <TVirtualFitter.h>

#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

using namespace std;

double ptbins[] ={17, 22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 400, 550, 790, 1000};
const int nbins = sizeof(ptbins)/sizeof(double) - 1;
double xmin=ptbins[2];
double xmax=ptbins[nbins];

const int ncen=10;
const char *centbin[ncen] = {"0-5%","5-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-100%"};

void  LoadStyle();
void drawText2(const char */*text*/, float /*xp*/, float /*yp*/, int /*size*/);
void makeMultiPanelCanvas(TCanvas*& /*canv*/,
                          const Int_t /*columns*/,
                          const Int_t /*rows*/,
                          const Float_t /*leftOffset*/,
                          const Float_t /*bottomOffset*/,
                          const Float_t /*leftMargin*/,
                          const Float_t /*bottomMargin*/,
                          const Float_t /*edge*/,const Float_t /*asyoffset*/); 
void MakeHist(TH1D *&/*hist*/,int /*istat*/, const char */*xname*/, const char */*yname*/);
void MakeHistRMS(TH1D *&/*hRMS*/,float /*max*/,float /*min*/);
void MakeHistMean(TH1D *&/*Mean*/,float /*max*/,float /*min*/);
void MakeZero(TH1D *&/*hist*/);



void fit_double_gaussian(TH1D *&hrsp);
int fit_dscb(TH1D *&hrsp,
             const double nsigma,
             const double jtptmin,
             const int niter,
             const string alg);

/// double sided crystal ball function definition  
double fnc_dscb(double*xx,double*pp);
template <class T>
bool from_string(T& t,
                 const std::string& s,
                 std::ios_base& (*f)(std::ios_base&))
{
  std::istringstream iss(s);
  return !(iss >> f >> t).fail();
}


void fit_gaussian(TH1D*& hrsp,
                  const double nsigma,
                  const double jtptmin,
                  const int niter);
void adjust_fitrange(TH1 *h,double& min,double& max);

void set_range_truncatedRMS(TH1D *&hist,float frac);


int PbPb_JRsp(string radii="4", string algtype="PuPF")
{

  LoadStyle();

  int iSave=0;
  double fracRMS = 1.000;
  cout <<" npt : " << nbins << endl;

  string calgo;
  int ist=-1;
  if(algtype == "PuPF"){ ist=-1; calgo = "akPu"+radii+"PF";  }
  else if(algtype == "PuCalo"){ ist= 5; calgo = "akPu"+radii+"Calo";}
  else if(algtype == "VsPF"  ){ ist=11; calgo = "akVs"+radii+"PF";  }
  else if(algtype == "VsCalo"){ ist=17; calgo = "akVs"+radii+"Calo";}

  if( radii == "1" )ist += 1;
  else if( radii == "2" )ist += 2;
  else if( radii == "3" )ist += 3;
  else if( radii == "4" )ist += 4;
  else if( radii == "5" )ist += 5;

  std::string inname="Histo_pbpb_hireco.root";//"input/Histo_PbPb_merged_trim_JEC.root";
  TFile *fin = new TFile(inname.c_str(),"r");
  
  TH1D *hrsp [ncen][nbins];
  TH1D *hMean[ncen], *hSigma[ncen];
  TF1 *fgaus=0;
  for(int ic=0; ic<ncen; ic++){

    hMean[ic] = new TH1D(Form("hMean_%d_%d",ist,ic),Form("Mean %s",calgo.c_str()),nbins,ptbins);
    hMean[ic]->SetMarkerColor(1);
    hMean[ic]->SetMarkerStyle(20);
    hMean[ic]->SetLineColor(1);
    hMean[ic]->SetMarkerSize(1.3);
    //MakeHistMean(hMean[ic],1.072,0.980); 
    MakeHistMean(hMean[ic],1.054,0.948); 

    //MakeHistMean(hMean[ic],1.012,0.780); 

    hSigma[ic] = new TH1D(Form("hSigma_%d_%d",ist,ic),Form("Sigma %s",calgo.c_str()),nbins,ptbins);
    hSigma[ic]->SetMarkerColor(1);
    hSigma[ic]->SetMarkerStyle(20);
    hSigma[ic]->SetLineColor(1);
    hSigma[ic]->SetMarkerSize(1.3);
    MakeHistRMS(hSigma[ic],0.563,0.001); 


    for(int ip=0; ip<nbins; ip++){
      
      hrsp[ic][ip] = (TH1D*)fin->Get(Form("%sJetAnalyzer/hrescrpt_genm%d_%d_%d",calgo.c_str(),ist,ic,ip));

      //set_range_truncatedRMS(hrsp[ic][ip], fracRMS);

      double norm  = hrsp[ic][ip]->GetMaximumStored();
      double mean  = hrsp[ic][ip]->GetMean();
      double emean = hrsp[ic][ip]->GetMeanError();
      // double sig   = hrsp[ic][ip]->GetRMS();
      // double esig  = hrsp[nj][ip]->GetRMSError();

      double sig   = hrsp[ic][ip]->GetRMS()/mean;
      double err   = hrsp[ic][ip]->GetRMSError();
      double esig = (pow(1/mean,2)*pow(err,2))+(pow(-sig/pow(mean,2),2)*pow(emean,2));
      esig = sqrt(esig);
      
      fgaus = new TF1("fgaus","gaus",0.0,2.00);
      fgaus->SetParameters(norm,1.00,0.2);
      //fgaus->SetParLimits(1,0.65,1.20);
      //fgaus->SetParLimits(2,0.1,1.0);

      int fitstatus = hrsp[ic][ip]->Fit(fgaus,"RQ","",0.55,1.25);
      //fitstatus=-1;
      mean  = (fitstatus!=0) ? hrsp[ic][ip]->GetMean()     :  fgaus->GetParameter(1);
      emean = (fitstatus!=0) ? hrsp[ic][ip]->GetMeanError():  fgaus->GetParError(1);
      sig   = (fitstatus!=0) ? hrsp[ic][ip]->GetRMS()      :  fgaus->GetParameter(2);
      esig  = (fitstatus!=0) ? hrsp[ic][ip]->GetRMSError() :  fgaus->GetParError(2);
      
      
      // fit_gaussian(hrsp[ic][ip], 1.0, 1.0, 5);      
      // TF1*  fgaus = (TF1*) hrsp[ic][ip]->GetListOfFunctions()->Last();
      // mean  = (fgaus==0) ? hrsp[ic][ip]->GetMean()     :  fgaus->GetParameter(1);      
      // emean = (fgaus==0) ? hrsp[ic][ip]->GetMeanError():  fgaus->GetParError(1);
      // sig   = (fgaus==0) ? hrsp[ic][ip]->GetRMS()/mean :  fgaus->GetParameter(2)/mean;
      // esig  = (fgaus==0) ? hrsp[ic][ip]->GetRMSError() :  fgaus->GetParError(2);

      //mean  = (fgaus==0) ? hrsp[ic][ip]->GetMean()     :  ((fgaus->GetParameter(1)+hrsp[ic][ip]->GetMean())*0.5);
      //emean = (fgaus==0) ? hrsp[ic][ip]->GetMeanError():   0.5*sqrt(pow(hrsp[ic][ip]->GetMeanError(),2)+pow(fgaus->GetParError(1),2));//fgaus->GetParError(1);
      //esig  = (fgaus==0) ? sqrt((pow(1/mean,2)*pow(hrsp[ic][ip]->GetRMSError(),2))+(pow(-hrsp[ic][ip]->GetRMS()/pow(mean,2),2)*pow(emean,2))) : sqrt((pow(1./fgaus->GetParameter(1),2)*pow(fgaus->GetParError(2),2))+pow(-fgaus->GetParameter(2)/pow(fgaus->GetParameter(1),2),2)*pow(fgaus->GetParError(1),2));


      //! Double sided crystal ball
       // fit_dscb(hrsp[ic][ip], 1.0, 1.0, 15, Form("ak%d%s",nj+1,ccalo.c_str()));
       // TF1*  frelrsp = (TF1*) hrsp[ic][ip]->GetListOfFunctions()->Last();
       // //mean  = (frelrsp==0) ? hrsp[ic][ip]->GetMean()     :  frelrsp->GetParameter(1);
       // //emean = (frelrsp==0) ? hrsp[ic][ip]->GetMeanError():  frelrsp->GetParError(1);
       // mean  = (frelrsp==0) ? hrsp[ic][ip]->GetMean()     :   ((frelrsp->GetParameter(1)+hrsp[ic][ip]->GetMean())*0.5);
       // emean = (frelrsp==0) ? hrsp[ic][ip]->GetMeanError():   0.5*sqrt(pow(hrsp[ic][ip]->GetMeanError(),2)+pow(frelrsp->GetParError(1),2));
       // sig   = (frelrsp==0) ? hrsp[nj][ip]->GetRMS()/mean :  frelrsp->GetParameter(2)/frelrsp->GetParameter(1);
       // esig  = (frelrsp==0) ? sqrt((pow(1/mean,2)*pow(hrsp[ic][ip]->GetRMSError(),2))+(pow(-hrsp[ic][ip]->GetRMS()/pow(mean,2),2)*pow(emean,2))) 
       // 	: sqrt((pow(1./frelrsp->GetParameter(1),2)*pow(frelrsp->GetParError(2),2))+pow(-frelrsp->GetParameter(2)/pow(frelrsp->GetParameter(1),2),2)*pow(frelrsp->GetParError(1),2));

      //sig   = (frelrsp==0) ? hrsp[ic][ip]->GetRMS()/mean :  frelrsp->GetParameter(2);
      //esig  = (frelrsp==0) ? hrsp[ic][ip]->GetRMSError() :  frelrsp->GetParError(2);

      //peak    =(frelrsp==0) ? hrsp[ic][j]->GetMean()     : ((frelrsp->GetParameter(1)+hrsp[ic][j]->GetMean())*0.5);
      //epeak   =(frelrsp==0) ? hrsp[ic][j]->GetMeanError(): 0.5*sqrt(pow(hrsp[ic][j]->GetMeanError(),2)+pow(frelrsp->GetParError(1),2));


      //fit_double_gaussian(hrsp[ic][ip]);
      // TF1*  frelrsp = (TF1*) hrsp[ic][ip]->GetFunction("fgaus");
      // mean         = (frelrsp==0) ? hrsp[ic][ip]->GetMean()     :  frelrsp->GetParameter(1);
      // double emean = (frelrsp==0) ? hrsp[ic][ip]->GetMeanError():  frelrsp->GetParError(1);
      // double sig   = (frelrsp==0) ? hrsp[ic][ip]->GetRMS()      :  frelrsp->GetParameter(2);
      // double esig  = (frelrsp==0) ? hrsp[ic][ip]->GetRMSError() :  frelrsp->GetParError(2);

      hMean [ic]->SetBinContent (ip+1, mean);
      hMean [ic]->SetBinError   (ip+1, emean);
      hSigma[ic]->SetBinContent (ip+1, sig);
      hSigma[ic]->SetBinError   (ip+1, esig);
    }
  }
  //return 0;

  int maxc=10;
  int maxr=2;
  int ipad=0;

  TLine *l0 = new TLine(xmin,1.00,xmax,1.0);
  l0->SetLineStyle(2);
  l0->SetLineWidth(2);
  TLine *l1 = new TLine(xmin,0.99,xmax,0.99);
  l1->SetLineStyle(2);
  TLine *l2 = new TLine(xmin,1.01,xmax,1.01);
  l2->SetLineStyle(2);

  TCanvas *c11 = new TCanvas("c11",Form("%s JES JER",calgo.c_str()),1828,408);
  makeMultiPanelCanvas(c11,maxc,maxr,0.0,0.0,0.22,0.22,0.02,0);
  for(int ic=ncen-1; ic>=0; ic--){
    c11->cd(++ipad);
    gPad->SetLogx();
    
    hSigma[ic]->GetXaxis()->SetRangeUser(xmin,xmax);
    hSigma[ic]->Draw("p");

    if(ipad==1){
      drawText2("PYTHIA + HYDJET",0.35,0.85,14);      
    }
    if(ipad==2){
      drawText2(calgo.c_str(),0.28,0.78,14);
    }
    if(ipad==3){
      drawText2("#sqrt{s_{NN}} = 5.02 TeV",0.25,0.85,14);
      drawText2("|#eta| < 2",0.25,0.76,14);
    }

    drawText2(centbin[ic],0.33,0.65,14);
    c11->cd(ipad+maxc);
    gPad->SetLogx();

    hMean[ic]->GetXaxis()->SetRangeUser(xmin,xmax);
    hMean[ic]->Draw("p");

    l1->Draw();
    l0->Draw();
    l2->Draw();
  }
  if(iSave){
    c11->SaveAs(Form("FinalPlots/JEC_%s.gif",calgo.c_str()));
    c11->SaveAs(Form("FinalPlots/JEC_%s.pdf",calgo.c_str()));
  }


  ipad=0;
  TCanvas *c99[ncen];
  for(int ic=0;ic<ncen; ic++){
    c99[ic] = new TCanvas(Form("c99_%d",ic),Form("%s %s Fitting FinalPlots",calgo.c_str(),centbin[ic]),100,102,1718,641);
    //c99[ic]->Divide(7,3,0,0);
    c99[ic]->Divide(6,4,0,0);
    ipad=0;
    for(int ip=0;ip<nbins;ip++){      
      c99[ic]->cd(++ipad);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.15);
      gPad->SetRightMargin(0.1);

      hrsp[ic][ip]->SetMaximum((hrsp[ic][ip]->GetMaximum() + 1.25*hrsp[ic][ip]->GetMaximum()));
      hrsp[ic][ip]->SetMinimum(0);
      hrsp[ic][ip]->SetTitle(0);
      hrsp[ic][ip]->GetXaxis()->SetTitle("<reco jet p_{T} / gen jet p_{T}>");
      hrsp[ic][ip]->GetXaxis()->SetTitleFont(42);
      hrsp[ic][ip]->GetXaxis()->SetLabelFont(42);
      hrsp[ic][ip]->GetXaxis()->SetLabelSize(0.06);
      hrsp[ic][ip]->GetXaxis()->SetTitleSize(0.06);
      hrsp[ic][ip]->GetXaxis()->SetNdivisions(505);
      hrsp[ic][ip]->GetYaxis()->SetTitle("");
      hrsp[ic][ip]->GetYaxis()->SetTitleFont(42);
      hrsp[ic][ip]->GetYaxis()->SetLabelFont(42);
      hrsp[ic][ip]->GetYaxis()->SetLabelSize(0.06);
      hrsp[ic][ip]->GetYaxis()->SetNdivisions(505);
      
      hrsp[ic][ip]->SetMarkerStyle(20);
      hrsp[ic][ip]->SetMarkerColor(1);
      hrsp[ic][ip]->SetLineColor(1);
      hrsp[ic][ip]->SetMarkerSize(1.1);
      hrsp[ic][ip]->Draw("p");  

      std::ostringstream strs; 
      strs << ptbins[ip] << "< p_{T} (GeV/c) <" << ptbins[ip+1];
      std::string spt = strs.str();
      
      if(ipad==1){drawText2(calgo.c_str(),0.28,0.90,19);      
	drawText2(spt.c_str(),0.22,0.80,15);		
      } else drawText2(spt.c_str(),0.17,0.80,15);		
    }
  }
  return 0;
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
void drawText2(const char *text, float xp, float yp, int size){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(43);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw();
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
  Float_t Xup [columns];
  Float_t Ylow[rows];
  Float_t Yup [rows];
  
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
void MakeHist(TH1D *&histo,int istat,const char *xname, const char *yname)
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
void MakeHistMean(TH1D *&h1,float ymax,float ymin)
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
  h1->GetYaxis()->SetTitle("#mu");
  h1->GetYaxis()->CenterTitle(true);
  h1->GetYaxis()->SetTitleSize(0.07);
  h1->GetYaxis()->SetTitleOffset(1.50);
  h1->GetYaxis()->SetLabelSize(0.07);
  h1->GetYaxis()->SetNdivisions(507);
  h1->GetYaxis()->SetDecimals(true);
  h1->GetYaxis()->SetTitleFont(42);
  h1->GetYaxis()->SetLabelFont(42);
}
void MakeHistRMS(TH1D *&h1,float ymax,float ymin)
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
  h1->GetYaxis()->SetTitle("#sigma / #mu");
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

void fit_gaussian(TH1D *&hrsp,
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

  rspMax =jtptmin/ptRefMax;

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

    adjust_fitrange(hrsp,fitrange_min,fitrange_max);
    fitfnc = new TF1("fgaus","gaus",fitrange_min,fitrange_max);
    fitfnc->SetParNames("N","#mu","#sigma");
    fitfnc->SetParameter(0,norm);
    fitfnc->SetParameter(1,peak);
    fitfnc->SetParameter(2,sigma);
    fitstatus = hrsp->Fit(fitfnc,"RQ0");
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
void adjust_fitrange(TH1 *h,double& min,double& max)
{
  int imin=1; while (h->GetBinLowEdge(imin)<min) imin++;
  int imax=1; while (h->GetBinLowEdge(imax)<max) imax++;
  while ((imax-imin)<8) {
    if (imin>1) {imin--; min = h->GetBinCenter(imin); }
    if (imax<h->GetNbinsX()-1) { imax++; max=h->GetBinCenter(imax); }
  }
}
void MakeZero(TH1D *&h1)
{
  h1->GetYaxis()->SetLabelSize(0);
  h1->GetYaxis()->SetTitleSize(0);
}
int fit_dscb(TH1D *&hrsp,
             const double nsigma,
             const double jtptmin,
             const int niter,
             const string alg)
{
  if (0==hrsp) {
    cout<<"ERROR: Empty pointer to fit_dscb()"<<endl;return -1;
  }


  // first use a gaussian to constrain crystal ball gaussian core
  fit_gaussian(hrsp, nsigma, jtptmin, niter);
  TF1* fgaus = hrsp->GetFunction("fgaus");
  if (0==fgaus) {
    hrsp->GetListOfFunctions()->Delete();
    return -1;
  }

  // implementation of the low pt bias threshold 
  string histname = hrsp->GetName();
  //double ptRefMax(1.0),rspMax(0.0);

  double fitrange_min(0.4);
  double fitrange_max(1.6);

  //cout <<" \t \t  xmin : "  << fitrange_min << "\t" << fitrange_max << endl;


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
  double A1  = pow(p1/fabs(a1),p1)*exp(-a1*a1/2);
  double A2  = pow(p2/abs(a2),p2)*exp(-a2*a2/2);
  double B1  = p1/fabs(a1) - fabs(a1);
  double B2  = p2/fabs(a2) - fabs(a2);

  double result(N);
  if      (u<-a1) result *= A1*pow(B1-u,-p1);
  else if (u<a2)  result *= exp(-u*u/2);
  else            result *= A2*pow(B2+u,-p2);
  return result;
}
void fit_double_gaussian(TH1D *&hrsp)
{
  if (0==hrsp) {
    cout<<"ERROR: Empty pointer to fit_double_gaussian()"<<endl;return;
  }
  
  string histname = hrsp->GetName();
  hrsp->Scale(1./hrsp->Integral());
  double mean     = hrsp->GetMean();
  double rms      = hrsp->GetRMS();

  int maxbin    = hrsp->GetMaximumBin();
  double norm1  = hrsp->GetBinContent(maxbin);
  double peak1  = hrsp->GetBinCenter(maxbin);
  double sigma1 = 0.04;

  double norm2  = norm1/2.0;
  double peak2  = mean;
  double sigma2 = 2*rms;

  //cout << " Mean  : "  << mean  << " \t  RMS  : " << rms    << endl;
  //cout << " norm1 : "  << norm1 << " \t  norm2 : " << norm2 << endl;
  //cout << " peak1 : "  << peak1 << " \t  sig1 : " << sigma1 << endl;
  //cout << " peak2 : "  << peak2 << " \t  sig2 : " << sigma2 << endl;

  double fitrange_min = 0.4;
  double fitrange_max = 1.4;

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
  fitstatus = hrsp->Fit(fitfnc,"RQ");

  fitfnc->SetParLimits(0,0.01,50*norm1);
  // fitfnc->SetParLimits(1,0.7,1.2);
  // fitfnc->SetParLimits(2,0.01,5.0);

  //fitfnc->SetParLimits(3,0.0,2*norm2);
  // fitfnc->SetParLimits(4,0.2,1.7);
  // fitfnc->SetParLimits(5,1.0,10.0);

  //fitfnc->SetParLimits(4,peak2-3.0*sigma2,peak2+3.0*sigma2);
  //fitfnc->SetParLimits(5,0.10,2.0*sigma2);


  // if (0!=fitstatus){
  //   fitfnc->SetParLimits(4,0.2,1.7);
  //   fitfnc->SetParLimits(5,2.5,20.0);
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

  // if(hrsp->GetFunction("fdgaus")==0){
  //   cout << "No function recorded in histogram " << hrsp->GetName() << endl;
  // }
  // if (0!=fitstatus){
  //   cout<<"fit_double_gaussian() to "<<hrsp->GetName()
  //       <<" failed. Fitstatus: "<<fitstatus
  //       <<" - FNC deleted."<<endl;
  //   hrsp->GetListOfFunctions()->Delete();
  // }
}
void set_range_truncatedRMS(TH1D *&hist,float frac)
{
  if (0==hist) return;

  const float nevts = hist->Integral(); if (0==nevts) return;
  const int   nbins = hist->GetNbinsX();

  if (frac<=0.0 || frac==1.) return;

  for (int ibin=1;ibin<nbins;++ibin) {
    int binx1   = ibin;
    int binx2   = nbins+1-ibin;
    float ievts = hist->Integral(binx1,binx2);

    if ( (ievts/nevts)>frac ) continue;
    else { hist->GetXaxis()->SetRange(binx1,binx2); break; }
  }
  return;
}
