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

//#include "utilities.h"

using namespace std;

const int knj = 6;
const char *calgo[knj] = {"ak1","ak2","ak3","ak4","ak5","ak6"};
double ptbins[] ={15, 17, 22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 400, 550, 790, 1000};
//double ptbins[] ={ 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 250, 300, 540, 1000};
//const double ptbins[]={20,50,80,120,200,350,800};
const int nbins = sizeof(ptbins)/sizeof(double) - 1;
double xmin=ptbins[0];
double xmax=ptbins[nbins-1];

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
void MakeHist(TH1F *&/*hist*/,int /*istat*/, const char */*xname*/, const char */*yname*/);
void MakeHistRMS(TH1F *&/*hRMS*/,float /*max*/,float /*min*/);
void MakeHistMean(TH1F *&/*Mean*/,float /*max*/,float /*min*/);
void MakeZero(TH1F *&/*hist*/);
TGraphErrors* TgraphIt(TH1F */*hist*/);


void fit_double_gaussian(TH1F *&hrsp);
int fit_dscb(TH1F *&hrsp,
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


void fit_gaussian(TH1F*& hrsp,
                  const double nsigma,
                  const double jtptmin,
                  const int niter);
void adjust_fitrange(TH1 *h,double& min,double& max);

void set_range_truncatedRMS(TH1F *&hist,float frac);


int JRsp(const char *ctype="Calo",std::string ceta="eta3")
{
  int iSave=0;
  double fracRMS = 1.0;
  cout <<" npt : " << nbins << endl;

  if(strcmp(ctype,"Calo")==0)xmin=ptbins[2];


  TVirtualFitter::SetDefaultFitter("Minuit");
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-04); 

  int ist=0;
  if(     strcmp(ctype,"Calo")  ==0)ist=knj;
  else if(strcmp(ctype,"PuPF")  ==0)ist=2*knj;
  else if(strcmp(ctype,"PuCalo")==0)ist=3*knj;
  else if(strcmp(ctype,"VsPF")  ==0)ist=4*knj;
  else if(strcmp(ctype,"VsCalo")==0)ist=5*knj;

  LoadStyle();
  std::string inname="";
  std::string etaname="";
  if( ceta == "eta2" ){
    //inname="JetResponse_histos_ppSignal_ppReco_PYTHIA_TuneCUETP8M1_5020GeV_patch3_eta2.root";
    //inname="Histo_jec_ppreco.root";
    //inname="Histo_calo_ppreco.root";

    //inname="Histo_ppreco.root";
    //inname="Histo_hireco.root";
    //inname="JetResponse_histos_ppSignal_ppReco_PYTHIA_TuneCUETP8M1_5020GeV_patch3_eta2_ppreco_nocut.root";

    //inname="Histo_hireco_26112015.root";
    //inname="Histo_ppreconew_30112015.root";

    //inname="Histo_ppreco_JEC753_757p1.root";

    inname="Histo_test_akCalo.root";

  }else{
    //
    //inname="Histo_ppReco_757p1.root";

    inname="Histo_ak26_hireco.root";
    etaname="|#eta| < 3.2";
  }
  TFile *fin = new TFile(inname.c_str(),"r");
  fin->ls();
  //return 0;

  TH1F *hrsp [knj][nbins];
  TH1F *hMean[knj], *hSigma[knj];
  TF1 *fgaus=0;
  for(int nj=0; nj<knj; nj++){

    hMean[nj] = new TH1F(Form("hMean_%d",nj),Form("Mean %s %s",calgo[nj],ctype),nbins,ptbins);
    hMean[nj]->SetMarkerColor(1);
    hMean[nj]->SetMarkerStyle(20);
    hMean[nj]->SetLineColor(1);
    hMean[nj]->SetMarkerSize(1.3);
    //MakeHistMean(hMean[nj],1.072,0.980); 
    //MakeHistMean(hMean[nj],1.022,0.978); 
    //MakeHistMean(hMean[nj],1.022,0.908); 

    MakeHistMean(hMean[nj],1.082,0.908); 

    hSigma[nj] = new TH1F(Form("hSigma_%d",nj),Form("Sigma %s %s",calgo[nj],ctype),nbins,ptbins);
    hSigma[nj]->SetMarkerColor(1);
    hSigma[nj]->SetMarkerStyle(20);
    hSigma[nj]->SetLineColor(1);
    hSigma[nj]->SetMarkerSize(1.3);
    MakeHistRMS(hSigma[nj],0.763,0.001); 

    if( nj==0 ) continue;


    std::ostringstream algn;
    if(strcmp(ctype,"PF"         )==0)algn << "ak"   << nj+1 << "PFJetAnalyzer";
    else if(strcmp(ctype,"Calo"  )==0)algn << "ak"   << nj+1 << "CaloJetAnalyzer";
    else if(strcmp(ctype,"PuPF"  )==0)algn << "akPu" << nj+1 << "PFJetAnalyzer";
    else if(strcmp(ctype,"PuCalo")==0)algn << "akPu" << nj+1 << "CaloJetAnalyzer";
    else if(strcmp(ctype,"VsPF"  )==0)algn << "akVs" << nj+1 << "PFJetAnalyzer";
    else if(strcmp(ctype,"VsCalo")==0)algn << "akVs" << nj+1 << "CaloJetAnalyzer";


    for(int ip=0; ip<nbins; ip++){
      
      hrsp[nj][ip] = (TH1F*)fin->Get(Form("%s/hrescrpt_genm%d_%d",algn.str().c_str(),nj+ist,ip));
      //hrsp[nj][ip]->Rebin(3);
      
      //! JEC closure from produced samples
      //hrsp[nj][ip] = (TH1F*)fin->Get(Form("%s/hrescrpt_genm%d_0_%d",algn.str().c_str(),nj+ist,ip));
      //hrsp[nj][ip]->Rebin(3);

      //hrsp[nj][ip] = (TH1F*)fin->Get(Form("%s%sJetAnalyzer/hresrrpt_genm%d_%d",calgo[nj],ctype,nj+ist,ip));
      //cout << Form("%s/hrescrpt_genm%d_%d",algn.str().c_str(),nj+ist,ip) << endl;

      set_range_truncatedRMS(hrsp[nj][ip], fracRMS);

      double norm  = hrsp[nj][ip]->GetMaximumStored();
      double mean  = hrsp[nj][ip]->GetMean();
      double emean = hrsp[nj][ip]->GetMeanError();
      double rms   = hrsp[nj][ip]->GetRMS();
      // double sig   = hrsp[nj][ip]->GetRMS();
      // double esig  = hrsp[nj][ip]->GetRMSError();

      double sig   = hrsp[nj][ip]->GetRMS()/mean;
      double err   = hrsp[nj][ip]->GetRMSError();
      double esig = (pow(1/mean,2)*pow(err,2))+(pow(-sig/pow(mean,2),2)*pow(emean,2));
      esig = sqrt(esig);

      int fitstatus = 0;
      //fgaus = new TF1("fgaus","gaus", 0.50, 1.50);
      fgaus = new TF1("fgaus","gaus", 0.50, 1.50);
      fgaus->SetParameters(norm, 0.9999, 0.1);
      //fgaus->SetParLimits(1,0.95,1.05);
      //fgaus->SetParLimits(1,0.85,1.25);
      if( strcmp(ctype,"PF") == 0 ){
	fgaus->SetParLimits(1,0.75,1.30);
	fgaus->SetParLimits(2,0.01,0.20);
	//fgaus->SetParLimits(1,0.90,1.10);
	//fgaus->SetParLimits(2,0.01,0.20);
	fitstatus = hrsp[nj][ip]->Fit(fgaus,"RQ");
	//fitstatus=-1;
      }else if( strcmp(ctype,"Calo") == 0 ){
	fgaus->SetParLimits(1,0.85,1.25);
	fgaus->SetParLimits(2,0.01,0.60);
	fitstatus = hrsp[nj][ip]->Fit(fgaus,"RQ","",0.85,1.15);
	//if( strcmp(calgo[nj],"ak1")==0)fitstatus=-1;

	//if( ip == 0 )fitstatus=-1;
	//fitstatus=-1;
      }
      
      //fitstatus = hrsp[nj][ip]->Fit(fgaus,"RQ");
      //cout <<" \t Fit Mean : " << fgaus->GetParameter(1) << " Hist Mean : " << mean << endl;
      // if(nj==2){
      // 	cout << "*******START*************** "<< ptbins[ip]  << " to " << ptbins[ip+1] << endl;
      // 	fitstatus = hrsp[nj][ip]->Fit(fgaus,"RQ");
      // 	cout <<" \t Fit Mean : " << fgaus->GetParameter(1) << " Hist Mean : " << mean << endl;
      // 	cout << "*******END***************** "<< ptbins[ip]  << " to " << ptbins[ip+1] << endl;
      // 	cout << endl;
      // }
      // else fitstatus = hrsp[nj][ip]->Fit(fgaus,"RQ");

      
      // fgaus = new TF1("fgaus","gaus", mean - 1.5*rms, mean + 1.5*rms);
      // fgaus->SetParameters(norm, 1.0, 0.25);
      // int fitstatus = hrsp[nj][ip]->Fit(fgaus,"RQ");
      //fitstatus=-1;

      mean  = (fitstatus!=0) ? hrsp[nj][ip]->GetMean()     :  fgaus->GetParameter(1);
      emean = (fitstatus!=0) ? hrsp[nj][ip]->GetMeanError():  fgaus->GetParError(1);
      sig   = (fitstatus!=0) ? hrsp[nj][ip]->GetRMS()/mean :  fgaus->GetParameter(2)/fgaus->GetParameter(1);
      esig  = (fitstatus!=0) ? hrsp[nj][ip]->GetRMSError() :  fgaus->GetParError(2);


      // mean  = (fitstatus!=0) ? hrsp[nj][ip]->GetMean()     :  ((fgaus->GetParameter(1)+hrsp[nj][ip]->GetMean())*0.5);
      // emean = (fitstatus!=0) ? hrsp[nj][ip]->GetMeanError():   0.5*sqrt(pow(hrsp[nj][ip]->GetMeanError(),2)+pow(fgaus->GetParError(1),2));//fgaus->GetParError(1);
      // sig   = (fitstatus!=0) ? hrsp[nj][ip]->GetRMS()/mean :  fgaus->GetParameter(2)/fgaus->GetParameter(1);
      // esig  = (fitstatus!=0) ? sqrt((pow(1/mean,2)*pow(hrsp[nj][ip]->GetRMSError(),2))+(pow(-hrsp[nj][ip]->GetRMS()/pow(mean,2),2)*pow(emean,2))) : sqrt((pow(1./fgaus->GetParameter(1),2)*pow(fgaus->GetParError(2),2))+pow(-fgaus->GetParameter(2)/pow(fgaus->GetParameter(1),2),2)*pow(fgaus->GetParError(1),2));


      // fit_gaussian(hrsp[nj][ip], 1.0, 1.0, 5);      
      // TF1*  fgaus = (TF1*) hrsp[nj][ip]->GetListOfFunctions()->Last();
      // mean  = (fgaus==0) ? hrsp[nj][ip]->GetMean()     :  fgaus->GetParameter(1);      
      // emean = (fgaus==0) ? hrsp[nj][ip]->GetMeanError():  fgaus->GetParError(1);
      // sig   = (fgaus==0) ? hrsp[nj][ip]->GetRMS()/mean :  fgaus->GetParameter(2)/mean;
      // esig  = (fgaus==0) ? hrsp[nj][ip]->GetRMSError() :  fgaus->GetParError(2);

      //mean  = (fgaus==0) ? hrsp[nj][ip]->GetMean()     :  ((fgaus->GetParameter(1)+hrsp[nj][ip]->GetMean())*0.5);
      //emean = (fgaus==0) ? hrsp[nj][ip]->GetMeanError():   0.5*sqrt(pow(hrsp[nj][ip]->GetMeanError(),2)+pow(fgaus->GetParError(1),2));//fgaus->GetParError(1);
      //esig  = (fgaus==0) ? sqrt((pow(1/mean,2)*pow(hrsp[nj][ip]->GetRMSError(),2))+(pow(-hrsp[nj][ip]->GetRMS()/pow(mean,2),2)*pow(emean,2))) : sqrt((pow(1./fgaus->GetParameter(1),2)*pow(fgaus->GetParError(2),2))+pow(-fgaus->GetParameter(2)/pow(fgaus->GetParameter(1),2),2)*pow(fgaus->GetParError(1),2));


      //! Double sided crystal ball
      // fit_dscb(hrsp[nj][ip], 1.0, 1.0, 15, Form("ak%d%s",nj+1,ctype));
      // TF1*  frelrsp = (TF1*) hrsp[nj][ip]->GetListOfFunctions()->Last();
      // mean  = (frelrsp==0) ? hrsp[nj][ip]->GetMean()     :  frelrsp->GetParameter(1);
      // emean = (frelrsp==0) ? hrsp[nj][ip]->GetMeanError():  frelrsp->GetParError(1);
      // // mean  = (frelrsp==0) ? hrsp[nj][ip]->GetMean()     :   ((frelrsp->GetParameter(1)+hrsp[nj][ip]->GetMean())*0.5);
      // // emean = (frelrsp==0) ? hrsp[nj][ip]->GetMeanError():   0.5*sqrt(pow(hrsp[nj][ip]->GetMeanError(),2)+pow(frelrsp->GetParError(1),2));
      // sig   = (frelrsp==0) ? hrsp[nj][ip]->GetRMS()/mean :  frelrsp->GetParameter(2)/frelrsp->GetParameter(1);
      // esig  = (frelrsp==0) ? sqrt((pow(1/mean,2)*pow(hrsp[nj][ip]->GetRMSError(),2))+(pow(-hrsp[nj][ip]->GetRMS()/pow(mean,2),2)*pow(emean,2))) 
      //  	: sqrt((pow(1./frelrsp->GetParameter(1),2)*pow(frelrsp->GetParError(2),2))+pow(-frelrsp->GetParameter(2)/pow(frelrsp->GetParameter(1),2),2)*pow(frelrsp->GetParError(1),2));

      //sig   = (frelrsp==0) ? hrsp[nj][ip]->GetRMS()/mean :  frelrsp->GetParameter(2);
      //esig  = (frelrsp==0) ? hrsp[nj][ip]->GetRMSError() :  frelrsp->GetParError(2);
      //peak    =(frelrsp==0) ? hrsp[i][j]->GetMean()     : ((frelrsp->GetParameter(1)+hrsp[i][j]->GetMean())*0.5);
      //epeak   =(frelrsp==0) ? hrsp[i][j]->GetMeanError(): 0.5*sqrt(pow(hrsp[i][j]->GetMeanError(),2)+pow(frelrsp->GetParError(1),2));


      // fit_double_gaussian(hrsp[nj][ip]);
      // TF1*  frelrsp = (TF1*) hrsp[nj][ip]->GetFunction("fgaus");
      // mean  = (frelrsp==0) ? hrsp[nj][ip]->GetMean()     :  frelrsp->GetParameter(1);
      // emean = (frelrsp==0) ? hrsp[nj][ip]->GetMeanError():  frelrsp->GetParError(1);
      // sig   = (frelrsp==0) ? hrsp[nj][ip]->GetRMS()      :  frelrsp->GetParameter(2);
      // esig  = (frelrsp==0) ? hrsp[nj][ip]->GetRMSError() :  frelrsp->GetParError(2);

      hMean[nj]->SetBinContent (ip+1, mean);
      hMean[nj]->SetBinError   (ip+1, emean);
      hSigma[nj]->SetBinContent(ip+1, sig);
      hSigma[nj]->SetBinError  (ip+1, esig);
    }
  }
  if( gPad )gPad->Close();
  //return 0;

  int maxc=6;
  int maxr=2;
  int ipad=0;

  TLine *l0 = new TLine(xmin,1.00,xmax,1.0);
  l0->SetLineStyle(2);
  l0->SetLineWidth(2);
  TLine *l1 = new TLine(xmin,0.99,xmax,0.99);
  l1->SetLineStyle(2);
  TLine *l2 = new TLine(xmin,1.01,xmax,1.01);
  l2->SetLineStyle(2);

  TCanvas *c11 = new TCanvas("c11","1-3 JES JER",1488,500/*1211,652*/);
  makeMultiPanelCanvas(c11,maxc,maxr,0.0,0.0,0.22,0.22,0.02,0);
  for(int nj=0; nj<knj; nj++){
    c11->cd(++ipad);


    gPad->SetLogx();
    //if( nj==0 ) continue;
    
    if(ipad!=1)MakeZero(hSigma[nj]);
    hSigma[nj]->GetXaxis()->SetRangeUser(xmin,xmax);
    hSigma[nj]->Draw("p");
    
    if(ipad==1){
      drawText2("PYTHIA TuneCUETP8M1",0.25,0.85,18);      
      drawText2("pp Reco (757p1)",0.55,0.75,18);
      drawText2("#sqrt{s} = 5.02 TeV",0.55,0.62,18);
      drawText2(etaname.c_str(),0.55,0.47,19);
    }
    if(ipad==2){
      //drawText2("#sqrt{s_{NN}} = 5.02 TeV",0.25,0.85,18);
    }
    if(ipad==3){
      drawText2("HI Reco",0.35,0.85,18);
      //drawText2("pp Reco (757p1)",0.35,0.85,18);
    }
    if(ipad==4){
      //drawText2(etaname.c_str(),0.25,0.85,21);
    }
    drawText2(Form("%s%s",calgo[nj],ctype),0.28,0.69,18);
    c11->cd(ipad+maxc);

    gPad->SetLogx();

    if(ipad!=1)MakeZero(hMean[nj]);
    hMean[nj]->GetXaxis()->SetRangeUser(xmin,xmax);
    hMean[nj]->Draw("p");

    l1->Draw();
    l0->Draw();
    l2->Draw();
  }
  if(iSave){
    c11->SaveAs(Form("FinalPlots/ppOfficial/JEC_ak123_%s.gif",ctype));
    c11->SaveAs(Form("FinalPlots/ppOfficial/JEC_ak123_%s.pdf",ctype));
  }

  // TGraphErrors *grMean[knj], *grSigma[knj];
  // TFile *fout = new TFile("JEClosure_757p1_with753JEC.root","RECREATE");
  // for(int nj=0;nj<knj;nj++){
  //   fout->mkdir(Form("%s%sJetAnalyzer",calgo[nj],ctype),Form("%s%sJetAnalyzer",calgo[nj],ctype));
  //   fout->cd(Form("%s%sJetAnalyzer",calgo[nj],ctype));

  //   grMean[nj]  = TgraphIt(hMean[nj]);
  //   grSigma[nj] = TgraphIt(hSigma[nj]);

  //   grMean[nj]->SetName(Form("grMean_%s%s",calgo[nj],ctype));
  //   grMean[nj]->SetTitle(Form("grMean_%s%s",calgo[nj],ctype));
  //   grMean[nj]->Write();

  //   grSigma[nj]->SetName(Form("grSigma_%s%s",calgo[nj],ctype));
  //   grSigma[nj]->SetTitle(Form("grSigma_%s%s",calgo[nj],ctype));
  //   grSigma[nj]->Write();


  //   // hMean [nj]->Write();
  //   // hSigma[nj]->Write();
  //   fout->cd("../");
  // }
  // fout->Close();


  return 0;

  // ipad=0;
  // TCanvas *c12 = new TCanvas("c12","4-6 JES JER",1211,652);
  // makeMultiPanelCanvas(c12,maxc,maxr,0.0,0.0,0.22,0.22,0.02,0);
  // for(int nj=3; nj<knj; nj++){
  //   c12->cd(++ipad);
  //   gPad->SetLogx();
  //   if(ipad!=1)MakeZero(hSigma[nj]);
  //   hSigma[nj]->GetXaxis()->SetRangeUser(xmin,xmax);
  //   hSigma[nj]->Draw("p");

  //   if(ipad==1){
  //     drawText2("PYTHIA TuneCUETP8M1",0.35,0.85,21);      
  //   }
  //   if(ipad==2){
  //     drawText2("#sqrt{s_{NN}} = 5.02 TeV",0.55,0.85,21);
  //     drawText2(etaname.c_str(),0.55,0.70,21);
  //   }
  //   drawText2(Form("%s%s",calgo[nj],ctype),0.28,0.78,21);
  //   c12->cd(ipad+maxc);
  //   gPad->SetLogx();
  //   if(ipad!=1)MakeZero(hMean[nj]);
  //   hMean[nj]->GetXaxis()->SetRangeUser(xmin,xmax);
  //   hMean[nj]->Draw("p");
  //   l1->Draw();
  //   l0->Draw();
  //   l2->Draw();
  // }
  // if(iSave){
  //   c12->SaveAs(Form("FinalPlots/ppOfficial/JEC_ak456_%s.gif",ctype));
  //   c12->SaveAs(Form("FinalPlots/ppOfficial/JEC_ak456_%s.pdf",ctype));
  // }

  // ipad=0;
  // TCanvas *c13 = new TCanvas("c13","7 JES JER",620,627);
  // makeMultiPanelCanvas(c13,1,2,0.0,0.0,0.22,0.22,0.02,0);
  // c13->cd(++ipad);
  // gPad->SetLogx();
  // hSigma[knj-1]->GetXaxis()->SetRangeUser(xmin,xmax);
  // hSigma[knj-1]->Draw("p");
  // drawText2("PYTHIA",0.55,0.85,21);      
  // drawText2("#sqrt{s_{NN}} = 2.76 TeV",0.55,0.75,21);
  // drawText2(etaname.c_str(),0.55,0.65,21);
  // drawText2(Form("%s%s",calgo[knj-1],ctype),0.28,0.78,21);
  // c13->cd(++ipad);
  // gPad->SetLogx();
  // hMean[knj-1]->GetXaxis()->SetRangeUser(xmin,xmax);
  // hMean[knj-1]->Draw("p");
  // l1->Draw();
  // l0->Draw();
  // l2->Draw();
  // if(iSave){
  //   c13->SaveAs(Form("FinalPlots/ppOfficial/JEC_ak7_%s.gif",ctype));
  //   c13->SaveAs(Form("FinalPlots/ppOfficial/JEC_ak7_%s.pdf",ctype));
  // }
  // //return 0;

  ipad=0;
  TCanvas *c99[knj];
  for(int nj=0;nj<knj;nj++){
    c99[nj] = new TCanvas(Form("c99_%d",nj),Form("%s Fitting FinalPlots",calgo[nj]),100,102,1656,566/*1399,942*/);
    c99[nj]->Divide(7,3,0,0);
    ipad=0;
    for(int ip=0;ip<nbins;ip++){      
      c99[nj]->cd(++ipad);
      // if(ipad%4==0)gPad->SetRightMargin(0.02);
      // if(ipad==1 || ipad==5 || ipad==9 || ipad==13)gPad->SetLeftMargin(0.15);

      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.15);
      gPad->SetRightMargin(0.1);

      hrsp[nj][ip]->SetMaximum((hrsp[nj][ip]->GetMaximum() + 0.45*hrsp[nj][ip]->GetMaximum()));
      hrsp[nj][ip]->SetMinimum(0);
      hrsp[nj][ip]->SetTitle(0);
      hrsp[nj][ip]->GetXaxis()->SetTitle("<reco jet p_{T} / gen jet p_{T}>");
      hrsp[nj][ip]->GetXaxis()->SetTitleFont(42);
      hrsp[nj][ip]->GetXaxis()->SetLabelFont(42);
      hrsp[nj][ip]->GetXaxis()->SetLabelSize(0.08);
      hrsp[nj][ip]->GetXaxis()->SetTitleSize(0.07);
      hrsp[nj][ip]->GetYaxis()->SetTitle("");
      hrsp[nj][ip]->GetYaxis()->SetTitleFont(42);
      hrsp[nj][ip]->GetYaxis()->SetLabelFont(42);
      hrsp[nj][ip]->GetYaxis()->SetLabelSize(0.08);

      hrsp[nj][ip]->SetMarkerStyle(20);
      hrsp[nj][ip]->SetMarkerColor(1);
      hrsp[nj][ip]->SetLineColor(1);
      hrsp[nj][ip]->SetMarkerSize(1.1);
      hrsp[nj][ip]->Draw("p");  

      std::ostringstream strs; 
      strs << ptbins[ip] << "< p_{T} (GeV/c) <" << ptbins[ip+1];
      std::string spt = strs.str();
      
      if(ipad==1){drawText2(Form("%s%s",calgo[nj],ctype),0.28,0.90,21);      
	drawText2(spt.c_str(),0.22,0.80,18);		
      } else drawText2(spt.c_str(),0.17,0.80,18);		
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
void MakeHist(TH1F *&histo,int istat,const char *xname, const char *yname)
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
void MakeHistMean(TH1F *&h1,float ymax,float ymin)
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
void MakeHistRMS(TH1F *&h1,float ymax,float ymin)
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

void fit_gaussian(TH1F *&hrsp,
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
void MakeZero(TH1F *&h1)
{
  h1->GetYaxis()->SetLabelSize(0);
  h1->GetYaxis()->SetTitleSize(0);
}
int fit_dscb(TH1F *&hrsp,
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
void fit_double_gaussian(TH1F *&hrsp)
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
void set_range_truncatedRMS(TH1F *&hist,float frac)
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
TGraphErrors* TgraphIt(TH1F* hist){

  TGraphErrors *tg;
  int nbins = hist->GetNbinsX();

  const int nlines = nbins;

  float pt[nlines], xsec[nlines];
  float pterr[nlines], xsecerr[nlines];

  for(int i = 0; i<nbins; i++ ){
    pt[i] = hist->GetBinCenter(i+1);
    xsec[i] = hist->GetBinContent(i+1);
    xsecerr[i] = hist->GetBinError(i+1);
    pterr[i] = 0;
  }

  tg = new TGraphErrors(nlines,pt,xsec,pterr,xsecerr);
  return tg;
}
