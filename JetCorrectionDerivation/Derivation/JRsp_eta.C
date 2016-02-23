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

using namespace std;

const int knj = 6;
const char *calgo[knj] = {"ak1","ak2","ak3","ak4","ak5","ak6"};
//double ptbins[] ={ 22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 400, 550, 790, 1000};
//const double ptbins[]={(15.+30.)/2.,(30.+50.)/2.,(50.+80.)/2.,(80.+120.)/2.,(120.+170.)/2.,(170.+220.)/2.,(220.+280.)/2.,(280.+370.)/2.,(370.+460.)/2.,(460.+540.)/2.,(540.+1000.)/2.};
const double ptbins[]={20,40,60,80,110,200,350,800};
const int nbins = sizeof(ptbins)/sizeof(double) - 1;
double xmin=ptbins[0];
double xmax=ptbins[nbins];

Int_t color[13]={kViolet+2,kBlue,kAzure+6,kGreen-3,kOrange-5,kOrange-3,kOrange+4,kRed-3,kRed,kRed+2,kMagenta+1,kRed+1,kGreen+3};

// const double etabins[]={-3.000,
//  			-2.500, -2.043, -1.930, -1.740, -1.653, -1.566, -1.392,
//  			-1.218, -1.131, -0.957, -0.879, -0.783, -0.609, -0.522,
//  			-0.435, -0.348, -0.261, -0.087,
//  			+0.000,
//  			+0.087, +0.261, +0.348, +0.435, +0.522, +0.609, +0.783,
//  			+0.879, +0.957, +1.131, +1.218, +1.392, +1.566,
//  			+1.653, +1.740, +1.930, +2.043, +2.500,
//  			+3.000
// };



// double etabins[]={-3.000,
//                   -2.500, -2.043, -1.740, -1.392,
//                   -1.131, -0.879, -0.609, -0.435, -0.261, -0.087,
//                   +0.000,
//                   +0.087, +0.261, +0.435, +0.609,
//                   +0.879, +1.131, +1.392,
//                   +1.740, +2.043, +2.500,
//                   +3.000
// };



//const double etabins[]={-3.0,-2.4,-1.8,-1.4,-1.0,-0.8,-0.6,-0.4,-0.2,
//			0.0,0.2,0.4,0.6,0.8,1.0,1.4,1.8,2.4,3.0};

//double etabins[] = {-3.000, -2.000, -1.3000, 0.000, 1.300, 2.000, 3.000};
//double etabins[] = {-3.000, -2.400, -2.000, -1.3000, 0.000, 1.300, 2.000, 2.400, 3.000};



// double etabins[] ={-3.000,
// 		   -2.500, -2.043, -1.740, -1.653, -1.566, -1.392,  
// 		   -1.218, -1.131, -0.957, -0.879, -0.783, -0.609, 
// 		   -0.522, -0.435, -0.348, -0.261, -0.087, 
// 		   +0.000, 
// 		   +0.087, +0.261, +0.348, +0.435, +0.522, +0.609, 
// 		   +0.783, +0.879, +0.957, +1.131, +1.218, +1.392, 
// 		   +1.566, +1.653, +1.740, +1.930, +2.043, +2.500,
// 		   +3.000
// };

// double etabins[] ={-4.363, -3.664,
//                    -3.139, -2.853, -2.500, -2.322, -2.172,
//                    -2.043, -1.930, -1.830, -1.740, -1.653, -1.566, -1.479, -1.392, -1.305,
//                    -1.218, -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522,
//                    -0.435, -0.348, -0.261, -0.174, -0.087,
//                    +0.000,
//                    +0.087, +0.174, +0.261, +0.348, +0.435, +0.522, +0.609, +0.696, +0.783,
//                    +0.879, +0.957, +1.044, +1.131, +1.218, +1.305, +1.392, +1.479, +1.566,
//                    +1.653, +1.740, +1.830, +1.930, +2.043, +2.172, +2.322, +2.500,
//                    +2.853, +3.139, +3.664, +4.363,
// };


// double etabins[]={-3.139, -2.853,
//                   -2.500, -2.043, -1.740, -1.392,
//                   -1.131, -0.879, -0.609, -0.435, -0.261, -0.087,
//                   +0.000,
//                   +0.087, +0.261, +0.435, +0.609,
//                   +0.879, +1.131, +1.392,
//                   +1.740, +2.043, +2.500, +2.853,
//                   +3.139
// };

double etabins[] = {-4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664,
                    -3.489, -3.314, -3.139, -2.964, -2.853, -2.650, -2.500, -2.322, -2.172,
                    -2.043, -1.930, -1.830, -1.740, -1.653, -1.566, -1.479, -1.392, -1.305,
                    -1.218, -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522,
                    -0.435, -0.348, -0.261, -0.174, -0.087,
                    +0.000,
                    +0.087, +0.174, +0.261, +0.348, +0.435, +0.522, +0.609, +0.696, +0.783,
                    +0.879, +0.957, +1.044, +1.131, +1.218, +1.305, +1.392, +1.479, +1.566,
                    +1.653, +1.740, +1.830, +1.930, +2.043, +2.172, +2.322, +2.500, +2.650,
                    +2.853, +2.964, +3.139, +3.314, +3.489, +3.664, +3.839, +4.013, +4.191,
                    +4.363, +4.538, +4.716
};
const int neta = sizeof(etabins)/sizeof(double) - 1;
// double xmineta=etabins[0];
// double xmaxeta=etabins[neta-1];
double xmineta=-2.964;
double xmaxeta=+2.964;


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

int JRsp_eta(const char *ctype="PF",const char *csel="lowpt")
{

  int iSave=0;

  int ist=0;
  if(     strcmp(ctype,"Calo")  ==0)ist=knj;
  else if(strcmp(ctype,"PuPF")  ==0)ist=2*knj;
  else if(strcmp(ctype,"PuCalo")==0)ist=3*knj;
  else if(strcmp(ctype,"VsPF")  ==0)ist=4*knj;
  else if(strcmp(ctype,"VsCalo")==0)ist=5*knj;

  LoadStyle();

  //! pp Reco JEC
  //std::string inname="JetResponse_histos_ppSignal_ppReco_official_53X_STARTHI53_V28_5_3_20_MC_eta3.root";
  //std::string inname="Histo_hireco_26112015.root";
  
  //std::string inname="Histo_ppreconew_30112015.root";
  
  //std::string inname="Histo_ppReco_757p1.root";

  //std::string inname="Histo_test_akCalo.root";

  std::string inname="Histo_ak26_hireco.root";  
  TFile *fin = new TFile(inname.c_str(),"r");

  cout << " npt : " << nbins << " neta "<< neta << endl;
  //return 0;

  TH1F *hrsp [knj][nbins][neta];
  TH1F *hMean[knj][nbins], *hSigma[knj][nbins];
  TF1 *fgaus=0;
  for(int nj=0; nj<knj; nj++){

    for(int ip=0; ip<nbins; ip++){
      
      hMean[nj][ip] = new TH1F(Form("hMean_%d_%d",nj,ip),Form("Mean %s %d",calgo[nj],ip),neta,etabins);
      hMean[nj][ip]->SetMarkerStyle(20);
      hMean[nj][ip]->SetMarkerColor(1);
      hMean[nj][ip]->SetLineColor(1);
      hMean[nj][ip]->SetMarkerSize(1.3);
      //MakeHistMean(hMean[nj][ip],1.024,0.974); 
      //MakeHistMean(hMean[nj][ip],1.054,0.954); 
      MakeHistMean(hMean[nj][ip],1.084,0.914); 
      
      hSigma[nj][ip] = new TH1F(Form("hSigma_%d_%d",nj,ip),Form("Sigma %s %d",calgo[nj],ip),neta,etabins);
      hSigma[nj][ip]->SetMarkerColor(1);
      hSigma[nj][ip]->SetMarkerStyle(20);
      hSigma[nj][ip]->SetLineColor(1);
      hSigma[nj][ip]->SetMarkerSize(1.3);
      //MakeHistRMS(hSigma[nj][ip],0.563,0.001); 
      MakeHistRMS(hSigma[nj][ip],0.863,0.001); 

      if( nj == 0 )continue;

      std::ostringstream algn;
      if(strcmp(ctype,"PF"         )==0)algn << "ak"   << nj+1 << "PFJetAnalyzer";
      else if(strcmp(ctype,"Calo"  )==0)algn << "ak"   << nj+1 << "CaloJetAnalyzer";
      else if(strcmp(ctype,"PuPF"  )==0)algn << "akPu" << nj+1 << "PFJetAnalyzer";
      else if(strcmp(ctype,"PuCalo")==0)algn << "akPu" << nj+1 << "CaloJetAnalyzer";
      else if(strcmp(ctype,"VsPF"  )==0)algn << "akVs" << nj+1 << "PFJetAnalyzer";
      else if(strcmp(ctype,"VsCalo")==0)algn << "akVs" << nj+1 << "CaloJetAnalyzer";

      for(int ie=0; ie<neta; ie++){
	//hrsp[nj][ip][ie] = (TH1F*)fin->Get(Form("%s/hrescrpt_genm_eta%d_%d_%d",algn.str().c_str(),nj+ist,ip,ie));

	//hrsp[nj][ip][ie] = (TH1F*)fin->Get(Form("%s/hrescrpt_wide_genm_eta%d_0_%d_%d",algn.str().c_str(),nj+ist,ip,ie));


	hrsp[nj][ip][ie] = (TH1F*)fin->Get(Form("%s/hrescrpt_wide_genm_eta%d_%d_%d",algn.str().c_str(),nj+ist,ip,ie));
	hrsp[nj][ip][ie]->SetMarkerStyle(20);
	//hrsp[nj][ip][ie]->Rebin(3);
	//hrsp[nj][ip][ie] = (TH1F*)fin->Get(Form("%s%sJetAnalyzer/hrescrpt_genm_eta%d_%d_%d",calgo[nj],ctype,nj+ist,ip,ie));
      
	if( hrsp[nj][ip][ie]->GetEntries() < 5 )continue;

	double norm  = hrsp[nj][ip][ie]->GetMaximumStored();
	double mean  = hrsp[nj][ip][ie]->GetMean();
	double emean = hrsp[nj][ip][ie]->GetMeanError();
	double sig   = hrsp[nj][ip][ie]->GetRMS();
	double esig  = hrsp[nj][ip][ie]->GetRMSError();

	// fgaus = new TF1("fgaus","gaus", mean - 2.*sig, mean + 2*sig);
	// fgaus->SetParameters(norm,1.0,1.5);
	// int fitstatus = hrsp[nj][ip][ie]->Fit(fgaus,"RQ");
	// mean  = (fitstatus!=0) ? hrsp[nj][ip][ie]->GetMean()     :  fgaus->GetParameter(1);
	// emean = (fitstatus!=0) ? hrsp[nj][ip][ie]->GetMeanError():  fgaus->GetParError(1);
	// sig   = (fitstatus!=0) ? hrsp[nj][ip][ie]->GetRMS()      :  fgaus->GetParameter(2);
	// esig  = (fitstatus!=0) ? hrsp[nj][ip][ie]->GetRMSError() :  fgaus->GetParError(2);
	// if(0!=fitstatus){
	//   cout << ctype << " \t  pt bin "  << ptbins[ip] << "  eta : " << etabins[ie] << endl;
	// }

	int fitstatus = 0;
	fgaus = new TF1("fgaus","gaus", 0.50, 1.50);
	fgaus->SetParameters(norm, 0.9999, 0.1);
	//fgaus->SetParLimits(1,0.95,1.05);
	//fgaus->SetParLimits(1,0.85,1.25);
	if( strcmp(ctype,"PF") == 0 ){
	  // fgaus->SetParLimits(1,0.75,1.10);
	  // fgaus->SetParLimits(2,0.01,0.20);
	  
	  //fgaus->SetParLimits(1,0.90,1.10);
	  //fgaus->SetParLimits(2,0.01,0.20);
	  
	  fitstatus = hrsp[nj][ip][ie]->Fit(fgaus,"RQ");
	  //fitstatus=-1;
	}else if( strcmp(ctype,"Calo") == 0 ){
	  // fgaus->SetParLimits(1,0.90,1.10);
	  // fgaus->SetParLimits(2,0.01,0.60);
	  fitstatus = hrsp[nj][ip][ie]->Fit(fgaus,"RQ");
	  //fitstatus=-1;
	}


	mean  = (fitstatus!=0) ? hrsp[nj][ip][ie]->GetMean()     :  fgaus->GetParameter(1);
	emean = (fitstatus!=0) ? hrsp[nj][ip][ie]->GetMeanError():  fgaus->GetParError(1);
	sig   = (fitstatus!=0) ? hrsp[nj][ip][ie]->GetRMS()/mean :  fgaus->GetParameter(2)/fgaus->GetParameter(1);
	esig  = (fitstatus!=0) ? hrsp[nj][ip][ie]->GetRMSError() :  fgaus->GetParError(2);

	hMean [nj][ip]->SetBinContent (ie+1, mean);
	hMean [nj][ip]->SetBinError   (ie+1, emean);
	hSigma[nj][ip]->SetBinContent (ie+1, sig);
	hSigma[nj][ip]->SetBinError   (ie+1, esig);
      }
    }
  }
  if ( gPad )gPad->Close();

  //hrsp[4][3][4]->Draw("p");
  //return 0;
  
  int maxc=6;
  int maxr=2;
  int ipad=0;

  int npt = 1; //! lowpt
  if(strcmp(csel,"lowpt")!=0)npt = 3; //! highpt

  std::ostringstream spt; 
  spt << ptbins[npt] << "< p_{T} (GeV/c) <" << ptbins[npt+1];
  ipad=0;
  TCanvas *c99[knj];
  for(int nj=0;nj<knj;nj++){
    if( nj == 0 )continue;
    c99[nj] = new TCanvas(Form("c99_%d",nj),Form("%s Fitting FinalPlots",calgo[nj]),100,102,1399,942);
    //c99[nj]->Divide(7,5,0,0);
    //c99[nj]->Divide(6,4,0,0);
    c99[nj]->Divide(6,10,0,0);
    ipad=0;
    for(int ie=0;ie<neta;ie++){      
      c99[nj]->cd(++ipad);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.15);
      gPad->SetRightMargin(0.1);

      hrsp[nj][npt][ie]->SetMaximum((hrsp[nj][npt][ie]->GetMaximum() + 0.25*hrsp[nj][npt][ie]->GetMaximum()));
      hrsp[nj][npt][ie]->SetMinimum(0);
      hrsp[nj][npt][ie]->SetTitle(0);
      hrsp[nj][npt][ie]->GetXaxis()->SetTitle("");
      hrsp[nj][npt][ie]->GetXaxis()->SetTitleFont(42);
      hrsp[nj][npt][ie]->GetXaxis()->SetLabelFont(42);
      hrsp[nj][npt][ie]->GetXaxis()->SetLabelSize(0.08);
      hrsp[nj][npt][ie]->GetXaxis()->SetTitleSize(0.07);
      hrsp[nj][npt][ie]->GetYaxis()->SetTitle("");
      hrsp[nj][npt][ie]->GetYaxis()->SetTitleFont(42);
      hrsp[nj][npt][ie]->GetYaxis()->SetLabelFont(42);
      hrsp[nj][npt][ie]->GetYaxis()->SetLabelSize(0.08);

      hrsp[nj][npt][ie]->SetMarkerStyle(20);
      hrsp[nj][npt][ie]->SetMarkerColor(1);
      hrsp[nj][npt][ie]->SetLineColor(1);
      hrsp[nj][npt][ie]->SetMarkerSize(1.1);
      hrsp[nj][npt][ie]->Draw("p");  

      std::ostringstream strs; 
      strs << etabins[ie] << "< #eta <" << etabins[ie+1];
      std::string spt = strs.str();
      
      if(ipad==1){drawText2(Form("%s%s",calgo[nj],ctype),0.28,0.90,17);      
	
	drawText2(spt.c_str(),0.22,0.80,21);		
      } else drawText2(spt.c_str(),0.17,0.80,21);		
    }
  }



  TLine *l0 = new TLine(xmineta,1.00,xmaxeta,1.0);
  l0->SetLineStyle(2);
  l0->SetLineWidth(2);
  TLine *l1 = new TLine(xmineta,0.99,xmaxeta,0.99);
  l1->SetLineStyle(2);
  TLine *l2 = new TLine(xmineta,1.01,xmaxeta,1.01);
  l2->SetLineStyle(2);

  ipad=0;
  TCanvas *c11 = new TCanvas("c11","JES JER",1488,500/*1211,652*/);
  makeMultiPanelCanvas(c11,maxc,maxr,0.0,0.0,0.22,0.22,0.02,0);
  for(int nj=0; nj<knj; nj++){
    c11->cd(++ipad);


    if(ipad!=1)MakeZero(hSigma[nj][npt]);
    hSigma[nj][npt]->Draw("p");
    if(ipad==1){
      drawText2("PYTHIA TuneCUETP8M1",0.25,0.85,18);      
      drawText2("pp Reco (757p1)",0.55,0.75,18);
      drawText2("#sqrt{s} = 5.02 TeV",0.55,0.62,18);
    }
    //if( nj == 0 )continue;

    if(ipad==2){
      //drawText2("#sqrt{s_{NN}} = 5.02 TeV",0.25,0.85,15);
      drawText2(spt.str().c_str(),0.25,0.66,18);
    }
    if(ipad==3){
      drawText2("HI Reco",0.35,0.85,18);
      //drawText2("pp Reco",0.35,0.85,18);
    }
    drawText2(Form("%s%s",calgo[nj],ctype),0.28,0.78,18);

    c11->cd(ipad+maxc);
    if(ipad!=1)MakeZero(hMean[nj][npt]);
    hMean[nj][npt]->Draw("p");
    // hMean[nj][0]->Draw("p");
    // hMean[nj][1]->Draw("psame");
    // hMean[nj][2]->Draw("psame");
    // hMean[nj][3]->Draw("psame");
    // hMean[nj][4]->Draw("psame");
    // hMean[nj][5]->Draw("psame");
    l1->Draw();
    l0->Draw();
    l2->Draw();
  }
  if(iSave){
    c11->SaveAs(Form("FinalPlots/ppOfficial/JEC_ak123_%s_eta_%s.gif",ctype,csel));
    c11->SaveAs(Form("FinalPlots/ppOfficial/JEC_ak123_%s_eta_%s.pdf",ctype,csel));
  }
  //return 0;





  // ipad=0;
  // TCanvas *c11 = new TCanvas("c11","1-3 JES JER",1211,652);
  // makeMultiPanelCanvas(c11,maxc,maxr,0.0,0.0,0.22,0.22,0.02,0);
  // for(int nj=0; nj<3; nj++){
  //   c11->cd(++ipad);
  //   if(ipad!=1)MakeZero(hSigma[nj][npt]);
  //   hSigma[nj][npt]->Draw("p");

  //   if(ipad==1){
  //     drawText2("PYTHIA",0.55,0.85,21);      
  //   }
  //   if(ipad==2){
  //     drawText2("#sqrt{s_{NN}} = 2.76 TeV",0.55,0.85,21);
  //     drawText2(spt.str().c_str(),0.35,0.66,21);
  //   }
  //   drawText2(Form("%s%s",calgo[nj],ctype),0.28,0.78,21);
  //   c11->cd(ipad+maxc);
  //   if(ipad!=1)MakeZero(hMean[nj][npt]);
  //   hMean[nj][npt]->Draw("p");
  // }
  // if(iSave){
  //   c11->SaveAs(Form("FinalPlots/ppOfficial/JEC_ak123_%s_eta_%s.gif",ctype,csel));
  //   c11->SaveAs(Form("FinalPlots/ppOfficial/JEC_ak123_%s_eta_%s.pdf",ctype,csel));
  // }
  // //return 0;

  // ipad=0;
  // TCanvas *c12 = new TCanvas("c12","4-6 JES JER",1211,652);
  // makeMultiPanelCanvas(c12,maxc,maxr,0.0,0.0,0.22,0.22,0.02,0);
  // for(int nj=3; nj<6; nj++){
  //   c12->cd(++ipad);
  //   if(ipad!=1)MakeZero(hSigma[nj][npt]);
  //   hSigma[nj][npt]->Draw("p");

  //   if(ipad==1){
  //     drawText2("PYTHIA",0.55,0.85,21);      
  //   }
  //   if(ipad==2){
  //     drawText2("#sqrt{s_{NN}} = 2.76 TeV",0.55,0.85,21);
  //     drawText2(spt.str().c_str(),0.35,0.66,21);
  //   }
  //   drawText2(Form("%s%s",calgo[nj],ctype),0.28,0.78,21);
  //   c12->cd(ipad+maxc);
  //   if(ipad!=1)MakeZero(hMean[nj][npt]);
  //   hMean[nj][npt]->Draw("p");
  // }
  // if(iSave){
  //   c12->SaveAs(Form("FinalPlots/ppOfficial/JEC_ak456_%s_eta_%s.gif",ctype,csel));
  //   c12->SaveAs(Form("FinalPlots/ppOfficial/JEC_ak456_%s_eta_%s.pdf",ctype,csel));
  // }
  // ipad=0;
  // TCanvas *c13 = new TCanvas("c13","0.7 JES JER",620,627);
  // makeMultiPanelCanvas(c13,1,2,0.0,0.0,0.22,0.22,0.02,0);
  // c13->cd(++ipad);
  // hSigma[knj-1][npt]->Draw("p");
  // drawText2("PYTHIA",0.55,0.85,21);      
  // drawText2("#sqrt{s_{NN}} = 2.76 TeV",0.55,0.75,21);
  // drawText2(spt.str().c_str(),0.35,0.66,21);
  // drawText2(Form("%s%s",calgo[knj-1],ctype),0.28,0.78,21);
  // c13->cd(++ipad);
  // hMean[knj-1][npt]->Draw("p");
  // if(iSave){
  //   c13->SaveAs(Form("FinalPlots/ppOfficial/JEC_ak7_%s_eta_%s.gif",ctype,csel));
  //   c13->SaveAs(Form("FinalPlots/ppOfficial/JEC_ak7_%s_eta_%s.pdf",ctype,csel));
  // }
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
  h1->GetXaxis()->SetRangeUser(xmineta,xmaxeta);
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
  h1->GetYaxis()->SetTitle("#mu");
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
  h1->GetXaxis()->SetRangeUser(xmineta,xmaxeta);
  h1->GetXaxis()->SetTitle("GenJet #eta ");
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
void MakeZero(TH1F *&h1)
{
  h1->GetYaxis()->SetLabelSize(0);
  h1->GetYaxis()->SetTitleSize(0);
}

