#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TLegend.h>
#include <TString.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TVirtualFitter.h>

#include "Math/GSLIntegrator.h"

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
/// default fit with gaussian in niter iteration of mean  
void fit_gaussian(TH1F*& hrsp,
                  const double nsigma,
                  const double jtptmin,
                  const int niter,
		  const string alg
);

int fit_dscb(TH1F*& hrsp,
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
void adjust_fitrange(TH1* h,double& min,double& max);
double GetPeak(TH1F *);

void perform_smart_fit(TGraphErrors * gabscor, TF1 * fabscor);
//double ptbins[] ={17, 22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 350, 400, 550, 790, 1000};

double ptbins[] ={15, 17, 20, 23, 27, 30, 35, 40, 45, 57, 72, 90, 120, 150, 200, 300, 400, 600, 1000};
const int npt = sizeof(ptbins)/sizeof(double) - 1;


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


//! Using for 757p1
// double etabins[]={-3.139, -2.853,
// 		  -2.500, -2.043, -1.740, -1.392,  
// 		  -1.131, -0.879, -0.609, -0.435, -0.261, -0.087, 
// 		  +0.000, 		
// 		  +0.087, +0.261, +0.435, +0.609,  
// 		  +0.879, +1.131, +1.392,  
// 		  +1.740, +2.043, +2.500, +2.853, 
// 		  +3.139
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
// 		    -2.500, -2.043, -1.930, -1.740, -1.653, -1.566, -1.392,  
// 		    -1.218, -1.131, -0.957, -0.879, -0.783, -0.609, -0.522, 
// 		    -0.435, -0.348, -0.261, -0.087, 
// 		    +0.000, 		  	  
// 		    +0.087, +0.261, +0.348, +0.435, +0.522, +0.609, +0.783, 
// 		    +0.879, +0.957, +1.131, +1.218, +1.392, +1.566, 
// 		    +1.653, +1.740, +1.930, +2.043, +2.500,  
// 		    +3.000
// };



// double etabins[] = {-3.000,
//    		    -2.500, -2.043, -1.740, -1.392,  
//    		    -1.131, -0.879, -0.609, -0.435, -0.261, -0.087, 
//    		    +0.000, 	  	  
//    		    +0.087, +0.261, +0.435, +0.609,  
//    		    +0.879, +1.131, +1.392,  
//    		    +1.740, +2.043, +2.500,  
//    		    +3.000,
// };



// double etabins[] = {-3.000,
//    		    -2.400, -1.392,  
//    		    -0.879, -0.609, -0.087, 
//    		    +0.000, 	  	  
//    		    +0.087, +0.609,  
//    		    +0.879, +1.392,  
//    		    +2.400,  
//    		    +3.000,
// };

const int neta = sizeof(etabins)/sizeof(double) - 1;


const int knj = 6;
string calgo[knj] = {"ak1","ak2","ak3","ak4","ak5","ak6"};

//const int knj = 1;
//string calgo[knj] = {"ak6"};

int CaloL2(const char *ctype="Calo")
{

  bool iSave=false;
  double minx=1.0;

  LoadStyle();
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(11111);

  ROOT::Math::GSLIntegrator GLSIntegrator(double absTol = 1.0000000000000001E-9, 
					  double relTol = 9.9999999999999995E-7, 
					  size_t size = 10000);

  string outdir="txtfiles/";
  //string outdir="junk/";

  cout << "npt : " << npt << "  " << "neta : " << neta <<endl;
  //return 0;
  
  //TFile *fin = new TFile("jra_hiF_akcalo_ppReco_757p1_dijet.root","r");

  //TFile *fin = new TFile("jra_hiF_akcalo_ppReco_757p1_dijet_fulleta_fullpt.root","r");

  TFile *fin = new TFile("jra_hiF_akcalo_ppReco_757p1_HcalRespCorrs_v4_00_mc_dijet.root","r");
  TFile *fl3 = new TFile("l3calo_input.root","r");

  TGraphErrors *grsp[knj]; //! l3 input

  TH1F *hl_rsp[knj][npt][neta], *hl_refpt[knj][npt][neta], *hl_jetpt[knj][npt][neta];

  TGraphErrors *grelcor[knj][neta], *gabscor[knj][neta], *gabsrsp[knj][neta];
  TF1 *fabscor=0, *frelcor=0, *frsp=0, *fl3rsp=0;
  TF1 *fgaus=0;

  int njmin = 0;
  int njmax = knj;
  for(int i=njmin; i<njmax; i++){
    // if( calgo[i] == "ak1" || calgo[i] == "ak2" )minx=10.0;
    // else minx=7.0;

    // if( calgo[i] == "ak1" || calgo[i] == "ak2" )minx=5.0;
    // else minx=5.0;
 
    std::ostringstream algn;
    algn << calgo[i] << "CaloJetAnalyzer/";    
    cout << " Running for "  << algn.str().c_str() << endl;

    // stringstream strValue;
    // strValue << calgo[i];
    // int iAlgo;
    // strValue >> iAlgo;
    //cout << " Running for algo "  << algn.str().c_str() << " iAlgo : "<< iAlgo <<endl;

    grsp[i] = (TGraphErrors*)fl3->Get(Form("%sL3RspVsRefPt_%s",algn.str().c_str(),calgo[i].c_str()));
    fl3rsp = (TF1*)grsp[i]->GetListOfFunctions()->First();

    for(int j=0; j<npt; j++){ //! pt
      std::ostringstream spt; 
      spt << ptbins[j] << "to" << ptbins[j+1];
      //cout << "\t  " << spt.str().c_str() <<endl;
      for(int k=0; k<neta; k++){ //! eta
	std::ostringstream seta; 
	seta << etabins[k] << "to" << etabins[k+1];
	//cout << "\t \t " << seta.str().c_str() <<endl;	
	std::string s1(algn.str()+"RefPt_JetEta"+seta.str()+"_RefPt"+spt.str());
	std::string s2(algn.str()+"JetPt_JetEta"+seta.str()+"_RefPt"+spt.str());
	std::string s3(algn.str()+"RelRsp_JetEta"+seta.str()+"_RefPt"+spt.str());

	hl_refpt[i][j][k] = (TH1F*)fin->Get(Form("%s",s1.c_str()));
	hl_jetpt[i][j][k] = (TH1F*)fin->Get(Form("%s",s2.c_str()));
	hl_rsp  [i][j][k] = (TH1F*)fin->Get(Form("%s",s3.c_str()));

	double tmppeak = GetPeak(hl_rsp[i][j][k]);

	if(j==0){ //! for each eta bin get the abs corr/rsp as a function of pT
	  gabscor[i][k] = new TGraphErrors();
	  gabscor[i][k]->SetName(("AbsCorVsJetPt_JetEta"+seta.str()).c_str());
	  gabsrsp[i][k] = new TGraphErrors();
	  gabsrsp[i][k]->SetName(("AbsRspVsRefPt_JetEta"+seta.str()).c_str());
	  grelcor[i][k] = new TGraphErrors();
	  grelcor[i][k]->SetName(("RelCorVsJetPt_JetEta"+seta.str()).c_str());
	}


	if (hl_rsp[i][j][k]->GetEntries() > 5 ) {
	  double norm = hl_rsp[i][j][k]->GetMaximumStored();
	  double rms  = hl_rsp[i][j][k]->GetRMS();

	  //! fit the with a fuction
	  //int fitstatus=-1;
	  //fgaus = new TF1("fgaus","gaus", tmppeak - 1.5*rms, tmppeak + 1.5*rms);
	  //fgaus->SetParameters(norm, tmppeak, 0.15);
	  //fitstatus = hl_rsp[i][j][k]->Fit(fgaus,"QR");      
	  //fit_dscb(hl_rsp[i][j][k], 1.03, minx, 5, Form("ak%d%s",i+1,ctype));      
	  
	  frsp   = (TF1*)hl_rsp[i][j][k]->GetListOfFunctions()->Last();
	  assert(hl_refpt[i][j][k]->GetEntries()>0 && hl_jetpt[i][j][k]->GetEntries()>0);
	
	  double refpt  =hl_refpt[i][j][k]->GetMean();
	  double erefpt =hl_refpt[i][j][k]->GetMeanError();
	  double jetpt  =hl_jetpt[i][j][k]->GetMean();
	  double ejetpt =hl_jetpt[i][j][k]->GetMeanError();

	  double peak;
	  double epeak;

	  //if( frsp!=0 )// cout << " yes the fit is there" << endl;
	  //cout << " frsp : " << frsp->GetParameter(1) <<endl;

	  peak  = (frsp==0) ? hl_rsp[i][j][k]->GetMean()     : frsp->GetParameter(1);
	  epeak = (frsp==0) ? hl_rsp[i][j][k]->GetMeanError(): frsp->GetParError(1);

	  double absrsp  = peak;
	  double eabsrsp = epeak;
	  double abscor  = 0.0;
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
	}else{
	  //cout << "  Hist has very few entries ::::: " << hl_rsp[i][j][k]->GetName() << " \t " << hl_rsp[i][j][k]->GetEntries() <<endl;
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
      if (xmin < minx) xmin=minx;

      // if( fabs(xmin - xmax) < 2.0 ){
      //   cout << " xmin is much smaller than xmax" <<endl;
      // 	xmin=4.0;
      // }

      if (npoints < 3) {
	gabscor[i][k]->SetPoint     (0, 10.0,1.0);
	gabscor[i][k]->SetPointError(0,  0.0,0.0);
	gabscor[i][k]->SetPoint     (1,100.0,1.0);
	gabscor[i][k]->SetPointError(1,  0.0,0.0);
	fabscor = new TF1("fit","[0]",10.0,100.0);
	fabscor->FixParameter(0,1.0);
      }else if(npoints > 2 && npoints < 10){
	//fabscor=new TF1("fit","[0]+[1]*log10(x)+[2]*(pow(log10(x),2))",xmin,xmax);
	fabscor=new TF1("fit","[0]+[1]*log10(x)+[2]/x",xmin,xmax);
	fabscor->SetParameter(0,1.0);
	fabscor->SetParameter(1,0.0);
	fabscor->SetParameter(2,0.0);
      }else{
	if (xmin<minx) xmin=minx;
	//! org
	//fabscor=new TF1("fit","[0] + [1]/(pow(log10(x),[2])+[3])",xmin,xmax);
	
	//fabscor=new TF1("fit","[0]+[1]/(pow(log10(x),[2])+[3])-[4]/x",xmin,xmax);
	fabscor=new TF1("fit","[0]+[1]/(pow(log10(x),[2])+[3])-[4]/x-[5]/x/x",xmin,xmax);
	fabscor->SetParameter(0,1.0);
	fabscor->SetParameter(1,5.0);
	fabscor->SetParameter(2,3.0);
	fabscor->SetParameter(3,3.0);
	fabscor->SetParLimits(3,0,100);
	fabscor->SetParameter(4,1.0);
	fabscor->SetParameter(5,1.0);
      }
      perform_smart_fit(gabscor[i][k],fabscor);
      gabscor[i][k]->GetListOfFunctions()->First()->ResetBit(TF1::kNotDraw);
      gabsrsp[i][k]->SetMarkerStyle(20);
      gabscor[i][k]->SetMarkerStyle(20);
    }


    for(int j=0; j<npt; j++){ //! pt
      for(int k=0; k<neta; k++){ //! eta    
	if( fabs(etabins[k]) > ketacut ) continue;
	if(hl_jetpt[i][j][k]->Integral()!=0) {
	  fabscor         = gabscor[i][k]->GetFunction("fit");
	  double jetpt    = hl_jetpt[i][j][k]->GetMean();
	  double refpt    = jetpt*fabscor->Eval(jetpt);
	  double l3cor    = fl3rsp->Eval(refpt);
	  double controlpt= refpt*l3cor;
	  double relcor   = controlpt/jetpt;
	  
	  // if(controlpt < 20) {
	  //   cout <<"  control pt : " << controlpt << " jetpt : " << jetpt << "\t  refpt : "<< refpt << " refptbin : " << hl_refpt[i][j][k]->GetMean() << endl;
	  // }
	  if (relcor > 5 && jetpt > minx)
	    cout<<"WARNING !!! suspicious point: "<<hl_jetpt[i][j][k]->GetName()
		<<", jet pt = "<<jetpt<<", ref pt = "<<refpt<<", control pt =  "<<controlpt<<", relcor = "<<relcor<<endl;
	  else {
	    grelcor[i][k]->SetPoint(j, jetpt, relcor);
	  }
	}
      }
    }

    string fnc_as_str;
    //fnc_as_str =  "([0] - [1]/(pow(log10(x),[2])+[3]))"; //! org
    //! using this one
    //fnc_as_str   = "[0]-[1]/(pow(log10(x),[2])+[3])+[4]/x";

    fnc_as_str   = "[0]-[1]/(pow(log10(x),[2])+[3])+[4]/x+[5]/x/x";

    string alg=Form("%d",i+1);
    string era="JEC_pp_PYTHIA_TuneCUETP8M1_5020GeV_patch3";
    TString txtfilename = outdir+era+"_L2Relative_AK"+alg+"Calo.txt";
    ofstream outf(txtfilename);
    //outf.setf(ios::right);
    for(int k=0; k<neta; k++){ //! eta  
      if( fabs(etabins[k]) > ketacut ) continue;
      double  xmin    = grelcor[i][k]->GetX()[0];
      double  xmax    = grelcor[i][k]->GetX()[grelcor[i][k]->GetN()-1];
      frelcor = new TF1("fitrelcor",fnc_as_str.c_str(),xmin,xmax);
      frelcor->SetParameter(0,0.0);
      frelcor->SetParameter(1,0.0);
      frelcor->SetParameter(2,0.0);
      frelcor->SetParameter(3,0.0);
      frelcor->SetParameter(4,0.0);
      frelcor->SetParameter(1,0.0);
      if (grelcor[i][k]->GetN()<2) {
	grelcor[i][k]->SetPoint(0,10,1.0);
	grelcor[i][k]->SetPoint(1,100,1.0);
	frelcor->FixParameter(1,0.0);
	frelcor->FixParameter(2,0.0);
	frelcor->FixParameter(3,0.0);
	frelcor->FixParameter(4,0.0);
	frelcor->FixParameter(5,0.0);
	cout<<" Unfortuante enough to be only with N<2" <<endl;
      }
      else if (grelcor[i][k]->GetN()==2) {
	frelcor->FixParameter(2,0.0);
	frelcor->FixParameter(3,0.0);
	frelcor->FixParameter(4,0.0);
	frelcor->FixParameter(5,0.0);
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
     c99[nj] = new TCanvas(Form("c99_%d",nj),Form("%sCalo Abs L2 Fitting plots",calgo[nj].c_str()),1854,866);
     c99[nj]->Divide(10,8,0,0);
     //c99[nj]->Divide(9,7,0,0);
     //c99[nj]->Divide(6,4,0,0);
     //c99[nj]->Divide(4,3,0,0);
     ipad=0;
     for(int ie=0;ie<neta;ie++){      
       c99[nj]->cd(++ipad);
       gPad->SetLogx();
       gPad->SetLeftMargin(0.15);
       gPad->SetRightMargin(0.01);
       gPad->SetBottomMargin(0.15);
       
       gabscor[nj][ie]->SetMaximum(3.14);
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
	 drawText2(Form("%s%s",calgo[nj].c_str(),ctype),0.27,0.90,18);      
       }
     }
     //c99[nj]->Close();
   }
   //return 0;

   //! Relative L2 Correction
   ipad=0;
   for(int nj=njmin; nj<njmax; nj++){
     c98[nj] = new TCanvas(Form("c98_%d",nj),Form("%sCalo Rel L2 Fitting plots",calgo[nj].c_str()),1854,866);
     c98[nj]->Divide(10,8,0,0);
     //c98[nj]->Divide(8,5,0,0);
     //c98[nj]->Divide(6,4,0,0);
     //c98[nj]->Divide(4,3,0,0);
     ipad=0;
     for(int ie=0;ie<neta;ie++){      

       c98[nj]->cd(++ipad);
       gPad->SetLogx();
       gPad->SetLeftMargin(0.15);
       gPad->SetRightMargin(0.01);
       gPad->SetBottomMargin(0.15);
    
       grelcor[nj][ie]->SetMaximum(2.38);
       grelcor[nj][ie]->SetMinimum(0.38);
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
	 drawText2(Form("%s%s",calgo[nj].c_str(),ctype),0.27,0.90,18);      
       }
     }
     if(iSave){
       c98[nj]->SaveAs(Form("CorrectionPlots/L2Relative_Corrections_%s%s_ppReco_CaloJets_757p1.pdf",calgo[nj].c_str(),ctype));
     }
   }
   
  return 0;
}
void perform_smart_fit(TGraphErrors * gabscor, TF1 * fabscor) {

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
  //   cout << "\t***ERROR***, fit for histo " << gabscor->GetName() << " diverges at low pt (<10 GeV/c)" << endl;
  
  //
  // check for failed fits
  // a chi2 of zero is symptomatic of a failed fit.
  //
  
  if (bestRChi2 < 0.001){
    cout<<"\t***ERROR***, FIT HAS FAILED for histo "<<gabscor->GetName()
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

    cout<<" fit for histo "<<gabscor->GetName()
        <<" has a reduced chi2="<<bestRChi2
        <<" after "<<fitIter<<" iterations"<<endl;
  }
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
  if( hrsp->GetEntries() < 10 )return -1;


  // first use a gaussian to constrain crystal ball gaussian core
  fit_gaussian(hrsp, nsigma, jtptmin, niter, alg);
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

  double fitrange_min(0.0);
  double fitrange_max(2.0);

  fitrange_min = std::max(rspMax,0.01);
  adjust_fitrange(hrsp,fitrange_min,fitrange_max);

  TF1* fdscb = new TF1("fdscb",fnc_dscb,fitrange_min,fitrange_max,7);
  fdscb->SetLineWidth(2);
  fdscb->SetLineStyle(2);

  double norm = 2.*fgaus->GetParameter(0);
  double mean = fgaus->GetParameter(1);
  double sigma= fgaus->GetParameter(2);

  //cout << hrsp->GetName() << "  fgaus "<< mean << " \t "  << hrsp->GetMean() << endl;

  // double norm = hrsp->GetMaximum();
  // double mean = hrsp->GetMean();
  // double sigma= hrsp->GetRMS();

  double aone(2.0),atwo(2.0),pone(10.0),ptwo(10.0);
  //double aone(5.0),atwo(5.0),pone(15.0),ptwo(15.0);
  TVirtualFitter::SetDefaultFitter("Minuit");
  
  int fitstatus(0);
  for (int i=0;i<niter;i++) {
    fdscb->SetParameter(0,2*norm); // N
    fdscb->SetParameter(1,mean); // mean
    fdscb->SetParameter(2,sigma);// sigma
    fdscb->SetParameter(3,aone); // a1
    fdscb->SetParameter(4,pone); // p1
    fdscb->SetParameter(5,atwo); // a2
    fdscb->SetParameter(6,ptwo); // p2                

    // fdscb->FixParameter(1,mean);
    // fdscb->FixParameter(2,sigma);
    

    if (i>0) fdscb->FixParameter(3,aone);
    else fdscb->SetParLimits(3,1.,20.);
    
    if (i>1) fdscb->FixParameter(5,atwo);
    else fdscb->SetParLimits(5,1.,20.);

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

    //reset sigma and mean to gauss values...
    fdscb->SetParameter(1,fgaus->GetParameter(1));
    fdscb->SetParError (1,fgaus->GetParError(1));
    fdscb->SetParameter(2,fgaus->GetParameter(2));
    fdscb->SetParError (2,fgaus->GetParError(2));
  }

  if (0!=fitstatus){
    cout<<"fit_fdscb() to "<<alg.c_str()<<"  " <<hrsp->GetName()
        <<" failed. Fitstatus: "<<fitstatus<<endl;
    hrsp->GetFunction("fdscb")->Delete();
  }
  else fdscb->ResetBit(TF1::kNotDraw);
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
                  const int niter,
		  const string alg)
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
    adjust_fitrange(hrsp,fitrange_min,fitrange_max);
    fitfnc = new TF1("fgaus","gaus",fitrange_min,fitrange_max);
    fitfnc->SetParNames("N","#mu","#sigma");
    fitfnc->SetParameter(0,norm);
    fitfnc->SetParameter(1,peak);
    fitfnc->SetParameter(2,sigma);
    fitstatus = hrsp->Fit(fitfnc,"RQ");
    delete fitfnc;
    fitfnc = hrsp->GetFunction("fgaus");
    fitfnc->ResetBit(TF1::kNotDraw);
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
    cout<<"fit_gaussian() to "<< alg.c_str() << "  " << hrsp->GetName()
        <<" failed. Fitstatus: "<<fitstatus
        <<" - FNC deleted."<<endl;
    hrsp->GetListOfFunctions()->Delete();
  }
}
double GetPeak(TH1F *hrsp)
{
  double peak=0;
  double max=-999;
  for(int ix=1; ix<=hrsp->GetNbinsX(); ix++){
    if(hrsp->GetBinContent(ix)>max){
      max  = hrsp->GetBinContent(ix);
      peak = hrsp->GetBinCenter(ix); 
    }
  }
  return peak;
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
