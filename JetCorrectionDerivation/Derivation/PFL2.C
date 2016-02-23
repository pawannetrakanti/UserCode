#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TString.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

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

void perform_smart_fit(int nj, TGraphErrors * gabscor, TF1 * fabscor);
//double ptbins[] ={22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 400, 550, 790, 1000};
//double ptbins[] ={17, 22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 350, 400, 550, 790, 1000};
double ptbins[] ={12,  12.5,  13,  13.5,  14,  15,  17,  20,  23,  27,  30,  35,  40,  45,  57,  
		  72,  90,  120,  150,  200,  300,  400,  600,  1000
};
const int npt = sizeof(ptbins)/sizeof(double) - 1;

// double etabins[] = {-3.000,
// 		    -2.500, -2.043, -1.930, -1.740, -1.653, -1.566, -1.392,  
// 		    -1.218, -1.131, -0.957, -0.879, -0.783, -0.609, -0.522, 
// 		    -0.435, -0.348, -0.261, -0.087, 
// 		    +0.000, 		  	  
// 		    +0.087, +0.261, +0.348, +0.435, +0.522, +0.609, +0.783, 
// 		    +0.879, +0.957, +1.131, +1.218, +1.392, +1.566, 
// 		    +1.653, +1.740, +1.930, +2.043, +2.500,  
// 		    +3.000
// };

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

// double etabins[] = {-3.000,
// 		    -2.500, -2.043, -1.740, -1.392,  
// 		    -1.131, -0.879, -0.609, -0.435, -0.261, -0.087, 
// 		    +0.000, 	  	  
// 		    +0.087, +0.261, +0.435, +0.609,  
// 		    +0.879, +1.131, +1.392,  
// 		    +1.740, +2.043, +2.500,  
// 		    +3.000,
// };


//! uisng 
// double etabins[] ={-5.191, -4.363, -3.664, 
// 		 -3.139, -2.853, -2.500, -2.322, -2.172, 
// 		 -2.043, -1.930, -1.830, -1.740, -1.653, -1.566, -1.479, -1.392, -1.305, 
// 		 -1.218, -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, 
// 		 -0.435, -0.348, -0.261, -0.174, -0.087, 
// 		 +0.000,       	       	       	       
// 		 +0.087, +0.174, +0.261, +0.348, +0.435, +0.522, +0.609, +0.696, +0.783, 
// 		 +0.879, +0.957, +1.044, +1.131, +1.218, +1.305, +1.392, +1.479, +1.566, 
// 		 +1.653, +1.740, +1.830, +1.930, +2.043, +2.172, +2.322, +2.500,  
// 		 +2.853, +3.139, +3.664,   
// 		 +4.363, +5.191
// };

const double ketacut=4.73;
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
// double etabins[]={-3.139, -2.853,
// 		  -2.500, -2.043, -1.740, -1.392,  
// 		  -1.131, -0.879, -0.609, -0.435, -0.261, -0.087, 
// 		  +0.000, 		
// 		  +0.087, +0.261, +0.435, +0.609,  
// 		  +0.879, +1.131, +1.392,  
// 		  +1.740, +2.043, +2.500, +2.853, 
// 		  +3.139
// };
const int neta = sizeof(etabins)/sizeof(double) - 1;


const int knj = 6;
const char *calgo[knj] = {"ak1","ak2","ak3","ak4","ak5","ak6"};

int PFL2(const char *ctype="PF")
{
  bool iSave=false;
  double minx=1.0;

  LoadStyle();
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(11111);


  string outdir="txtfiles/";
  //string outdir="junk/";
  cout << "npt : " << npt << "  " << "neta : " << neta <<endl;
  //return 0;


  //TFile *fin = new TFile("jra_hiF_akpf_dijet_HIreco.root","r");
  //TFile *fin = new TFile("jra_hiF_akpf_dijet_ppreco_nocut.root","r");

  //TFile *fin = new TFile("jra_hiF_akpf_dijet_hireco.root","r");

  //TFile *fin = new TFile("jra_hiF_akpf_ppReco_new_dijet.root","r");

  //TFile *fin = new TFile("jra_hiF_akpf_ppReco_757p1_dijet.root","r");

   TFile *fin = new TFile("jra_hiF_akpf_ppReco_757p1_HcalRespCorrs_v4_00_mc_dijet.root","r");
  TFile *fl3 = new TFile("l3pf_input.root","r");

  TGraphErrors *grsp[knj]; //! l3 input

  TH1F *hl_rsp[knj][npt][neta], *hl_refpt[knj][npt][neta], *hl_jetpt[knj][npt][neta];

  TGraphErrors *grelcor[knj][neta], *gabscor[knj][neta], *gabsrsp[knj][neta];
  TF1 *fabscor=0, *frelcor=0, *frsp=0, *fl3rsp=0;
  
  int njmin = 0;
  int njmax = knj;
  for(int i=njmin; i<njmax; i++){

    if( strcmp(calgo[i],"ak1")==0 || strcmp(calgo[i],"ak2")==0 )minx=4.0;
    else minx=4.0;
    
    std::ostringstream algn;
    algn << calgo[i] << "PFJetAnalyzer/";    
    cout << " Running for "  << algn.str().c_str() << endl;

    grsp[i] = (TGraphErrors*)fl3->Get(Form("%sL3RspVsRefPt_%s",algn.str().c_str(),calgo[i]));
    fl3rsp = (TF1*)grsp[i]->GetListOfFunctions()->First();

    for(int j=0; j<npt; j++){ //! pt
      std::ostringstream spt; 
      spt << ptbins[j] << "to" << ptbins[j+1];
      for(int k=0; k<neta; k++){ //! eta
	std::ostringstream seta; 
	seta << etabins[k] << "to" << etabins[k+1];
	
	std::string s1(algn.str()+"RefPt_JetEta"+seta.str()+"_RefPt"+spt.str());
	std::string s2(algn.str()+"JetPt_JetEta"+seta.str()+"_RefPt"+spt.str());
	std::string s3(algn.str()+"RelRsp_JetEta"+seta.str()+"_RefPt"+spt.str());

	hl_refpt[i][j][k] = (TH1F*)fin->Get(Form("%s",s1.c_str()));
	hl_jetpt[i][j][k] = (TH1F*)fin->Get(Form("%s",s2.c_str()));
	hl_rsp  [i][j][k] = (TH1F*)fin->Get(Form("%s",s3.c_str()));

	if(j==0){ //! for each eta bin get the abs corr/rsp as a function of pT
	  gabscor[i][k] = new TGraphErrors();
	  gabscor[i][k]->SetName(("AbsCorVsJetPt_JetEta"+seta.str()).c_str());
	  gabsrsp[i][k] = new TGraphErrors();
	  gabsrsp[i][k]->SetName(("AbsRspVsRefPt_JetEta"+seta.str()).c_str());	  
	  grelcor[i][k] = new TGraphErrors();
	  grelcor[i][k]->SetName(("RelCorVsJetPt_JetEta"+seta.str()).c_str());
	}

	if (hl_rsp[i][j][k]->GetEntries() > 10) {
	  frsp    = (TF1*)hl_rsp[i][j][k]->GetListOfFunctions()->Last();
	  assert(hl_refpt[i][j][k]->GetEntries()>0 && hl_jetpt[i][j][k]->GetEntries()>0);

	  double refpt  =hl_refpt[i][j][k]->GetMean();
	  double erefpt =hl_refpt[i][j][k]->GetMeanError();
	  double jetpt  =hl_jetpt[i][j][k]->GetMean();
	  double ejetpt =hl_jetpt[i][j][k]->GetMeanError();

	  double peak;
	  double epeak;

	  peak  = (frsp==0) ? hl_rsp[i][j][k]->GetMean()     : frsp->GetParameter(1);
	  epeak = (frsp==0) ? hl_rsp[i][j][k]->GetMeanError(): frsp->GetParError(1);

	  double absrsp = peak;
	  double eabsrsp = epeak;
	  double abscor = 0.0;
	  double eabscor = 0.0;

	  if (absrsp > 0){
	    abscor  =1.0/absrsp;
	    eabscor = abscor*abscor*epeak;
	  }
	  if ((abscor>0) && (absrsp>0) && (eabscor>1e-5) && (eabscor/abscor<0.5) && (eabsrsp>1e-4) && (eabsrsp/absrsp<0.5)){
	    gabsrsp[i][k]->SetPoint(j,refpt,absrsp);
	    gabsrsp[i][k]->SetPointError(j,erefpt,eabsrsp);	    
	    gabscor[i][k]->SetPoint(j,jetpt,abscor);
	    gabscor[i][k]->SetPointError(j,ejetpt,eabscor);	    
	  } 
	}
      }
    }

    //! All the graph filling is done... now try to fit
    for(int k=0; k<neta; k++){ //! eta    
      if( fabs(etabins[k]) > ketacut ) continue;
      int npoints = gabscor[i][k]->GetN();
      double xmin(1.0),xmax(100.0);
      if (npoints > 0){
	xmin = gabscor[i][k]->GetX()[0];
	xmax = gabscor[i][k]->GetX()[gabscor[i][k]->GetN()-1];
      }
      //if (xmin < 7) xmin=7;

      if (npoints < 4) {
	gabscor[i][k]->SetPoint     (0, 10.0,1.0);
	gabscor[i][k]->SetPointError(0,  0.0,0.0);
	gabscor[i][k]->SetPoint     (1,100.0,1.0);
	gabscor[i][k]->SetPointError(1,  0.0,0.0);
	fabscor = new TF1("fit","[0]",10.0,100.0);
	fabscor->FixParameter(0,1.0);
      }else if(npoints > 4 && npoints < 10){
	fabscor=new TF1("fit","[0]+[1]*log10(x)+[2]*pow(log10(x),2)",xmin,xmax);
	fabscor->SetParameter(0,1.0);
	fabscor->SetParameter(1,0.0);
	fabscor->SetParameter(2,0.0);
      }else{
	fabscor = new TF1("fit","[0]+[1]/(pow(log10(x),2)+[2])+[3]*exp(-[4]*(log10(x)-[5])*(log10(x)-[5]))",xmin,xmax);
	fabscor->SetParameter(0,0.5);
	fabscor->SetParameter(1,9.0);
	fabscor->SetParameter(2,8.0);
	fabscor->SetParameter(3,-0.3);
	fabscor->SetParameter(4,0.6);
	fabscor->SetParameter(5,1.0);

	fabscor->SetParLimits(2,0.1,100);
	fabscor->SetParLimits(3,-100,0);
	fabscor->SetParLimits(4,0,100);



	// fabscor = new TF1("fit","[0]+[1]/(pow(log10(x),2)+[2])+[3]*exp(-[4]*(log10(x)-[5])*(log10(x)-[5])) + [6]/pow(log10(x),2)",xmin,xmax);
	// fabscor->SetParameter(0,0.5);
	// fabscor->SetParameter(1,9.0);
	// fabscor->SetParameter(2,8.0);
	// fabscor->SetParameter(3,-0.3);
	// fabscor->SetParameter(4,0.8);
	// fabscor->SetParameter(5,2.0);
	// fabscor->SetParameter(6,1.0);

 	// fabscor->SetParLimits(2,0.1,100);
	// fabscor->SetParLimits(3,-100,0);
	// fabscor->SetParLimits(4,0,100);

      }
      perform_smart_fit(i, gabscor[i][k],fabscor);
      gabscor[i][k]->GetListOfFunctions()->First()->ResetBit(TF1::kNotDraw);
      gabsrsp[i][k]->SetMarkerStyle(20);
      gabscor[i][k]->SetMarkerStyle(20);
    }

    for(int j=0; j<npt; j++){ //! pt
      for(int k=0; k<neta; k++){ //! eta    
	if(hl_jetpt[i][j][k]->Integral()!=0) {
	  fabscor         = gabscor[i][k]->GetFunction("fit");
	  double jetpt    = hl_jetpt[i][j][k]->GetMean();
	  double refpt    = jetpt*fabscor->Eval(jetpt);
	  double l3cor    = fl3rsp->Eval(refpt);
	  double controlpt=refpt*l3cor;
	  double relcor   =controlpt/jetpt;
	  if (relcor > 5 && jetpt>15)
	    cout<<"WARNING !!! suspicious point: "<<hl_jetpt[i][j][k]->GetName()
		<<", jet pt = "<<jetpt<<", ref pt = "<<refpt<<", control pt =  "<<controlpt<<", relcor = "<<relcor<<endl;
	  else {
	    grelcor[i][k]->SetPoint(j, jetpt, relcor);
	  }
	}
      }
    }

    string fnc_as_str;
    fnc_as_str  = "([0]+[1]/(pow(log10(x),2)+[2])+[3]*exp(-[4]*(log10(x)-[5])*(log10(x)-[5])))";
    string alg=Form("%d",i+1);
    string era="JEC_pp_PYTHIA_TuneCUETP8M1_5020GeV_patch3";
    TString txtfilename = outdir + era + "_L2Relative_AK" + alg + "PF.txt";
    ofstream outf(txtfilename);
    outf.setf(ios::right);

    for(int k=0; k<neta; k++){ //! eta  
      double  xmin    = grelcor[i][k]->GetX()[0];
      double  xmax    = grelcor[i][k]->GetX()[grelcor[i][k]->GetN()-1];
      frelcor = new TF1("fitrelcor",fnc_as_str.c_str(),xmin,xmax);
      frelcor->SetParameter(0,0.0);
      frelcor->SetParameter(1,0.0);
      frelcor->SetParameter(2,0.0);
      frelcor->SetParameter(3,0.0);
      frelcor->SetParameter(4,0.0);
      frelcor->SetParameter(5,0.0);
      frelcor->SetParameter(6,0.0);
      if (grelcor[i][k]->GetN()<2) {
	grelcor[i][k]->SetPoint(0,10,1.0);
	grelcor[i][k]->SetPoint(1,100,1.0);
	frelcor->FixParameter(1,0.0);
	frelcor->FixParameter(2,0.0);
	frelcor->FixParameter(3,0.0);
	frelcor->FixParameter(4,0.0);
	frelcor->FixParameter(5,0.0);
	frelcor->FixParameter(6,0.0);
	cout<<" Unfortuante enough to be only with N<2" <<endl;
      }
      else if (grelcor[i][k]->GetN()==2) {
	frelcor->FixParameter(2,0.0);
	frelcor->FixParameter(3,0.0);
	frelcor->FixParameter(4,0.0);
	frelcor->FixParameter(5,0.0);
	frelcor->FixParameter(6,0.0);
	cout<<" Unfortuante enough to be only with N==2" <<endl;
      }
      grelcor[i][k]->Fit(frelcor,"QRB0");
      grelcor[i][k]->GetListOfFunctions()->First()->ResetBit(TF1::kNotDraw);
      grelcor[i][k]->SetMarkerStyle(20);
     
      if(k==0){
	//cout<<"{1 JetEta 1 JetPt "<<fnc_as_str<<" Correction L2Relative}"<<endl;
	outf<<"{1 JetEta 1 JetPt "<<fnc_as_str<<" Correction L2Relative}"<<endl;
      }
      double  etamin  = etabins[k];
      double  etamax  = etabins[k+1];
      double  ptmin   = grelcor[i][k]->GetX()[0];
      double  ptmax   = grelcor[i][k]->GetX()[grelcor[i][k]->GetN()-1];
      // cout<<setw(11)<<etamin
      // 	  <<setw(11)<<etamax
      // 	  <<setw(11)<<(int)(frelcor->GetNpar()+2) //Number of parameters + 2
      // 	  <<setw(12)<<ptmin
      // 	  <<setw(12)<<ptmax;

       outf<<setw(11)<<etamin
	  <<setw(11)<<etamax
	  <<setw(11)<<(int)(frelcor->GetNpar()+2) //Number of parameters + 2
	  <<setw(12)<<ptmin
	  <<setw(12)<<ptmax;

      for(int p=0; p<frelcor->GetNpar(); p++){
	//cout<<setw(13)<<frelcor->GetParameter(p);
	outf<<setw(13)<<frelcor->GetParameter(p);
      }
      //cout<<endl;
      outf<<endl;
    }
    outf.close();
  }//! i jet algo loop
  //return 0;

   int ipad=0;
   //! Absolute L2 corrections
   TCanvas *c99[knj], *c98[knj];
   for(int nj=njmin; nj<njmax; nj++){
     c99[nj] = new TCanvas(Form("c99_%d",nj),Form("%s Abs L2 Fitting plots",calgo[nj]),1854,866);
     //c99[nj]->Divide(9,7,0,0);
     c99[nj]->Divide(6,4,0,0);
     //c99[nj]->Divide(8,5,0,0);
     ipad=0;
     for(int ie=0;ie<neta;ie++){      
       c99[nj]->cd(++ipad);
       gPad->SetLogx();
       gPad->SetLeftMargin(0.15);
       gPad->SetRightMargin(0.01);
       gPad->SetBottomMargin(0.15);
       
       gabscor[nj][ie]->SetMaximum(3.84);
       gabscor[nj][ie]->SetMinimum(0.34);
       gabscor[nj][ie]->SetTitle(0);
       //gabscor[nj][ie]->GetXaxis()->SetRangeUser(30,1000);      
       gabscor[nj][ie]->GetXaxis()->SetTitle("< raw p_{T} > (GeV/c)");
       gabscor[nj][ie]->GetXaxis()->SetTitleFont(42);
       gabscor[nj][ie]->GetXaxis()->SetLabelFont(42);
       gabscor[nj][ie]->GetXaxis()->SetLabelSize(0.08);
       gabscor[nj][ie]->GetXaxis()->SetTitleSize(0.07);
       gabscor[nj][ie]->GetXaxis()->SetNdivisions(507);
       gabscor[nj][ie]->GetXaxis()->SetNoExponent();
       gabscor[nj][ie]->GetXaxis()->SetMoreLogLabels();
       gabscor[nj][ie]->GetYaxis()->SetTitle("Abs L2 Corr.");
       gabscor[nj][ie]->GetYaxis()->SetTitleFont(42);
       gabscor[nj][ie]->GetYaxis()->SetLabelFont(42);
       gabscor[nj][ie]->GetYaxis()->SetLabelSize(0.08);
       gabscor[nj][ie]->GetYaxis()->SetTitleSize(0.07);
       gabscor[nj][ie]->GetYaxis()->SetNdivisions(507);
       
       gabscor[nj][ie]->Draw("ap");  
       
       std::ostringstream seta; 
       seta << etabins[ie] << " < #eta < " << etabins[ie+1];
       drawText2(seta.str().c_str(),0.27,0.80,15);
       if(ipad==1){
	 drawText2(Form("%s%s",calgo[nj],ctype),0.27,0.90,18);      
       }
     }
     //c99[nj]->Close();
   }
   //return 0;

   //! Relative L2 Correction
   ipad=0;
   for(int nj=njmin; nj<njmax; nj++){
     c98[nj] = new TCanvas(Form("c98_%d",nj),Form("%s Rel L2 Fitting plots",calgo[nj]),1854,866);
     //c98[nj]->Divide(9,7,0,0);
     c98[nj]->Divide(6,4,0,0);
     ipad=0;
     for(int ie=0;ie<neta;ie++){      

       c98[nj]->cd(++ipad);
       gPad->SetLogx();
       gPad->SetLeftMargin(0.15);
       gPad->SetRightMargin(0.01);
       gPad->SetBottomMargin(0.15);
    
       grelcor[nj][ie]->SetMaximum(3.84);
       grelcor[nj][ie]->SetMinimum(0.28);
       grelcor[nj][ie]->SetTitle(0);
       //grelcor[nj][ie]->GetXaxis()->SetRangeUser(30,1000);      
       grelcor[nj][ie]->GetXaxis()->SetTitle("< raw p_{T} > (GeV/c)");
       grelcor[nj][ie]->GetXaxis()->SetTitleFont(42);
       grelcor[nj][ie]->GetXaxis()->SetLabelFont(42);
       grelcor[nj][ie]->GetXaxis()->SetLabelSize(0.08);
       grelcor[nj][ie]->GetXaxis()->SetTitleSize(0.07);
       grelcor[nj][ie]->GetXaxis()->SetNdivisions(507);
       grelcor[nj][ie]->GetXaxis()->SetNoExponent();
       grelcor[nj][ie]->GetXaxis()->SetMoreLogLabels();
       grelcor[nj][ie]->GetYaxis()->SetTitle("Rel L2 Corr.");
       grelcor[nj][ie]->GetYaxis()->SetTitleFont(42);
       grelcor[nj][ie]->GetYaxis()->SetLabelFont(42);
       grelcor[nj][ie]->GetYaxis()->SetLabelSize(0.08);
       grelcor[nj][ie]->GetYaxis()->SetTitleSize(0.07);
       grelcor[nj][ie]->GetYaxis()->SetNdivisions(507);

       grelcor[nj][ie]->SetMarkerStyle(20);
       grelcor[nj][ie]->Draw("ap");  

       std::ostringstream seta; 
       seta << etabins[ie] << " < #eta < " << etabins[ie+1];
       drawText2(seta.str().c_str(),0.27,0.80,15);
       if(ipad==1){
	 drawText2(Form("%s%s",calgo[nj],ctype),0.27,0.90,18);      
       }
     }
     if(iSave){
       c98[nj]->SaveAs(Form("CorrectionPlots/L2Relative_Corrections_%s%s_ppReco_PFJets_757p1.pdf",calgo[nj],ctype));
     } 
   }
   
  return 0;
}
void perform_smart_fit(int nj, TGraphErrors * gabscor, TF1 * fabscor) {

  int maxFitIter = 150;
  int fitIter = 0;
  vector<double> bestPars;
  double bestRChi2 = 0;
  do {
    //
    // do the fit, get the results and the parameters of the fitted function
    //

    TFitResultPtr fitResPtr = gabscor->Fit(fabscor,"RQ0S");
    vector<double> auxPars = fitResPtr.Get()->Parameters();

    //
    // compute the reduced chi2 of this fit and if it is the best fit so far
    // then save the parameters
    //
    double rchi2 = fitResPtr.Get()->Chi2()/ fitResPtr.Get()->Ndf();
    if (fitResPtr.Get()->Ndf() == 0) rchi2 = 0;
    if (rchi2 > 0 && (rchi2<bestRChi2 || bestRChi2==0)){
      bestRChi2 = rchi2;
      bestPars  = auxPars;
    }
    
    //
    // increment the counter
    //
    fitIter++;
  }while(( bestRChi2 > 2 || bestRChi2 == 0 ) && fitIter < maxFitIter);
  
  //
  // set the best parameters and chi2 to the fit function
  //

  TF1 * ffh = gabscor->GetFunction("fit");
  for (unsigned int np=0;np < bestPars.size() ; np++){
    ffh->SetParameter(np,bestPars[np]);
    fabscor->SetParameter(np,bestPars[np]);
  }
  fabscor->SetChisquare(bestRChi2 * fabscor->GetNDF());
  ffh->SetChisquare(bestRChi2 * fabscor->GetNDF());

  //
  // warn if the fit diverges at low pt
  //
  // if (fabscor->Integral(0,10) > 25)
  //   cout << "\t***ERROR***, fit for histo " << calgo[nj] << " " << gabscor->GetName() << " diverges at low pt (<10 GeV/c)" << endl;
  
  //
  // check for failed fits
  // a chi2 of zero is symptomatic of a failed fit.
  //
  
if (bestRChi2 < 0.001){
  cout<<"\t***ERROR***, FIT HAS FAILED for histo "<< calgo[nj] << "  " << gabscor->GetName()
      <<" which has a reduced chi2="<<bestRChi2
      <<" after "<<fitIter<<" iterations. "<<endl;
  }

  //
  // check for large reduced chi2's
  // above 10 is a plain error; between 5 and 10 is a warning
  //

 if (bestRChi2 > 5){
    if (bestRChi2 > 10)
      cout<<"\t***ERROR***,";
    else
      cout<<"\tWARNING,";

    cout<<" fit for histo "<< calgo[nj] << " " << gabscor->GetName()
        <<" has a reduced chi2="<<bestRChi2
        <<" after "<<fitIter<<" iterations"<<endl;
  }
}

