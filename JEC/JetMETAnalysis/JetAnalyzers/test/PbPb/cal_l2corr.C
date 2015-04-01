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

void perform_smart_fit(TGraphErrors * gabscor, TF1 * fabscor);


double ptbins[] = {10,15,20,27,35,45,57,72,90,120,150,200,300,400,550,750,1000};
// double ptbins[]={10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15, 16, 17, 18, 19, 20, 22, 25, 27, 30, 32, 35, 37, 40, 42, 45, 47,
//  		 50, 52, 55, 57, 60, 62, 65, 67, 70, 72, 75, 77, 80, 85, 90, 95, 100,
//  		 120, 150, 200, 300, 400, 550, 750, 1000
//};
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

double etabins[] = {-3.000,
		    -2.500, -2.043, -1.740, -1.392,  
		    -1.131, -0.879, -0.609, -0.435, -0.261, -0.087, 
		    +0.000, 	  	  
		    +0.087, +0.261, +0.435, +0.609,  
		    +0.879, +1.131, +1.392,  
		    +1.740, +2.043, +2.500,  
		    +3.000,
};
const int neta = sizeof(etabins)/sizeof(double) - 1;


const int knj = 7;
const char *calgo[knj] = {"ak1","ak2","ak3","ak4","ak5","ak6","ak7"};

int cal_l2corr(const char *ctype="PF")
{


  cout << "npt : " << npt << "  " << "neta : " << neta <<endl;

  // for(int i=0; i<npt; i++){
  //   cout <<" pt : "<< ptbins[i] << setw(15) << "  log10(pt)  : " << setw(15) << log10(ptbins[i]) << endl;
  // }
  // return 0;
  
  LoadStyle();
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(11111);

  //TFile* fin = new TFile("jra_hiF_ak_l2_dijet_full_finebin_lowpt.root","r");

  //TFile* fin = new TFile("jra_hiF_ak_l2_dijet_final_lowpt.root","r");
  //TFile* fin = new TFile("jra_hiF_ak_l2_dijet_official_lowpt.root","r");
  TFile* fin = new TFile("jra_hiF_ak_l2_dijet_test_lowpt.root","r");

  TGraphErrors *grelcor[knj][neta], *gabscor[knj][neta];
  TF1 *fabscor=0;

  int njmin = 2;
  int njmax = knj;
  //for(int i=njmin; i<njmin+1/*knj*/; i++){
  for(int i=njmin; i<njmax; i++){

    std::ostringstream algn;
    if(strcmp(ctype,"PF"         )==0)algn << "ak"   << i+1 << "PFJetAnalyzer/";
    else if(strcmp(ctype,"Calo"  )==0)algn << "ak"   << i+1 << "CaloJetAnalyzer/";
    else if(strcmp(ctype,"PuPF"  )==0)algn << "akPu" << i+1 << "PFJetAnalyzer/";
    else if(strcmp(ctype,"PuCalo")==0)algn << "akPu" << i+1 << "CaloJetAnalyzer/";
    else if(strcmp(ctype,"VsPF"  )==0)algn << "akVs" << i+1 << "PFJetAnalyzer/";
    else if(strcmp(ctype,"VsCalo")==0)algn << "akVs" << i+1 << "CaloJetAnalyzer/";

    cout << " Doing fits for ...." << algn.str().c_str() << endl;

    for(int j=0; j<neta; j++){
      std::ostringstream seta; 
      seta << etabins[j] << "to" << etabins[j+1];

      std::string sone(algn.str()+"AbsCorVsJetPt_JetEta"+seta.str());
      std::string stwo(algn.str()+"RelCorVsJetPt_JetEta"+seta.str());
      
      gabscor[i][j] = (TGraphErrors*)fin->Get(Form("%s",sone.c_str()));
      grelcor[i][j] = (TGraphErrors*)fin->Get(Form("%s",stwo.c_str()));

      assert(gabscor[i][j]->GetN()>0 && grelcor[i][j]->GetN()>0);
      int npoints = gabscor[i][j]->GetN();
      double xmin(1.0),xmax(100.0);
      if (npoints > 0){
	xmin = gabscor[i][j]->GetX()[0];
	xmax = gabscor[i][j]->GetX()[gabscor[i][j]->GetN()-1];
      }


      //if (xmin < 15) xmin=15;
      //cout << i << " " << j << "  xmin :  "<<  xmin << " xmax : " << xmax << endl;

      if (npoints < 3) {
	gabscor[i][j]->SetPoint     (0, 10.0,1.0);
	gabscor[i][j]->SetPointError(0,  0.0,0.0);
	gabscor[i][j]->SetPoint     (1,100.0,1.0);
	gabscor[i][j]->SetPointError(1,  0.0,0.0);
	fabscor = new TF1("fit","[0]",10.0,100.0);
	fabscor->FixParameter(0,1.0);
      }
      else{
	//fabscor = new TF1("fit","[0]+[1]/(pow(log10(x),2)+[2])+[3]*exp(-[4]*(log10(x)-[5])*(log10(x)-[5]))",xmin,xmax);

	if(strcmp(ctype,"PF")==0 || strcmp(ctype,"PuPF")==0 || strcmp(ctype,"VsPF")==0){
	  fabscor = new TF1("fit","[0]+[1]/(pow(log10(x),2)+[2])+[3]*exp(-[4]*(log10(x)-[5])*(log10(x)-[5]))",xmin,xmax);
	  fabscor->SetParameter(0,0.5);
	  fabscor->SetParameter(1,9.0);
	  fabscor->SetParameter(2,8.0);
	  fabscor->SetParameter(3,-0.3);
	  fabscor->SetParameter(4,0.6);
	  fabscor->SetParameter(5,1.00);
	  //fabscor->SetParameter(5,1.25);

	  fabscor->SetParLimits(2,0.1,100);
	  fabscor->SetParLimits(3,-100,0);
	  fabscor->SetParLimits(4,0,100);

	  if (xmax < 15) {
	    fabscor->FixParameter(0,1.0);
	    fabscor->FixParameter(1,0.0);
	    fabscor->FixParameter(2,0.0);
	    fabscor->FixParameter(3,0.0);
	    fabscor->FixParameter(4,0.0);
	    fabscor->FixParameter(5,0.0);
	    fabscor->FixParameter(6,0.0);
	    cout<<" !!!!!!!!!!  xmax less than 20 " <<xmax<<endl;
	  }
	}else{
	  fabscor=new TF1("fit","[0]+[1]/(pow(log10(x),[2])+[3])",xmin,xmax);
	  //fabscor=new TF1("fit","[0]+[1]/(pow(log10(x),[2]))+[3]*exp(-[4]*pow(log10(x)-[5],2))",xmin,xmax);
	  fabscor->SetParameter(0,1.0);
	  fabscor->SetParameter(1,5.0);
	  fabscor->SetParameter(2,3.0);
	  fabscor->SetParameter(3,3.0);
	  fabscor->SetParLimits(3,0,100);
	  //fabscor->SetParameter(4,1.0);
	  //fabscor->SetParameter(5,1.0);
	}
      }
      perform_smart_fit(gabscor[i][j],fabscor);
      gabscor[i][j]->GetListOfFunctions()->First()->ResetBit(TF1::kNotDraw);
      gabscor[i][j]->SetMarkerStyle(20);
      //continue;

      string fnc_as_str; 
      if (strcmp(ctype,"PF")==0 || 
	  strcmp(ctype,"PuPF")==0 || 
	  strcmp(ctype,"VsPF")==0 
	  ) fnc_as_str  = "(x>=15)*([0]+[1]/(pow(log10(x)-[2],2)+[2])+[3]*exp(-[4]*(log10(x)-[5])*(log10(x)-[5])))+(x<15)*1";
      else  if (strcmp(ctype,"Calo")==0 ||
		strcmp(ctype,"PuCalo")==0 ||
		strcmp(ctype,"VsCalo")==0
		) {
	fnc_as_str =  "(x>=15)*([0]+[1]/(pow(log10(x),[2])+[3]))+(x<15)*1";
	//fnc_as_str =  "(x>=15)*([0]+[1]/(pow(log10(x),[2]))+[3]*exp(-[4]*pow(log10(x)-[5],2)))+(x<15)*1";
      }
      xmin    = grelcor[i][j]->GetX()[0];
      xmax    = grelcor[i][j]->GetX()[grelcor[i][j]->GetN()-1];
      TF1*    frelcor = new TF1("fit",fnc_as_str.c_str(),xmin,xmax);
      frelcor->SetParameter(0,0.0);
      frelcor->SetParameter(1,0.0);
      frelcor->SetParameter(2,0.0);
      frelcor->SetParameter(3,0.0);
      frelcor->SetParameter(4,0.0);
      frelcor->SetParameter(5,0.0);

      if (grelcor[i][j]->GetN()<2) {
	grelcor[i][j]->SetPoint(0,10,1.0);
	grelcor[i][j]->SetPoint(1,100,1.0);
	frelcor->FixParameter(1,0.0);
	frelcor->FixParameter(2,0.0);
	frelcor->FixParameter(3,0.0);
	frelcor->FixParameter(4,0.0);
	frelcor->FixParameter(5,0.0);
	//cout<<" Unfortuante enough to be only with N<2" <<endl;  
      }else if (grelcor[i][j]->GetN()==2) {
	frelcor->FixParameter(2,0.0);
	frelcor->FixParameter(3,0.0);
	frelcor->FixParameter(4,0.0);
	frelcor->FixParameter(5,0.0);
	//cout<<" Unfortuante enough to be only with N==2" <<endl;  
      }
      //grelcor[i][j]->Fit(frelcor,"QRB0");
      grelcor[i][j]->Fit(frelcor,"RB");
      grelcor[i][j]->GetListOfFunctions()->First()->ResetBit(TF1::kNotDraw);
      grelcor[i][j]->SetMarkerStyle(20);

      //cout <<" \t \t " << sone.c_str() << " " << gacor[i][j]->GetN() << "\t \t  " << stwo.c_str() << endl;
    }//! j eta
    cout<<endl;
  }//! i knj
  //grcor[1][35]->Draw("ap");
  //return 0;






  int ipad=0;
  //! neta =38;
  //! Absolute L2 corrections
  TCanvas *c99[knj], *c98[knj];
  for(int nj=njmin; nj<njmax; nj++){
    c99[nj] = new TCanvas(Form("c99_%d",nj),Form("%s Abs L2 Fitting plots",calgo[nj]),1854,866);
    //c99[nj]->Divide(4,5,0,0);
    //c99[nj]->Divide(8,5,0,0);
    c99[nj]->Divide(6,4,0,0);
    ipad=0;
    for(int ie=0;ie<neta;ie++){      
      c99[nj]->cd(++ipad);
      gPad->SetLogx();
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.01);
      gPad->SetBottomMargin(0.15);
      
      gabscor[nj][ie]->SetMaximum(2.14);
      gabscor[nj][ie]->SetMinimum(0.64);
      gabscor[nj][ie]->SetTitle(0);
      gabscor[nj][ie]->GetXaxis()->SetRangeUser(30,1000);      
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
      if(ipad==1)drawText2(ctype,0.27,0.90,15);      
    }
    //c99[nj]->Close();
  }
  //return 0;

  //! Relative L2 Correction
  ipad=0;
  for(int nj=njmin; nj<njmax; nj++){
    c98[nj] = new TCanvas(Form("c98_%d",nj),Form("%s Rel L2 Fitting plots",calgo[nj]),1854,866);
    //c98[nj]->Divide(4,5,0,0);
    //c98[nj]->Divide(8,5,0,0);
    c98[nj]->Divide(6,4,0,0);
    ipad=0;
    for(int ie=0;ie<neta;ie++){      

      c98[nj]->cd(++ipad);
      gPad->SetLogx();
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.01);
      gPad->SetBottomMargin(0.15);
      
      grelcor[nj][ie]->SetMaximum(1.38);
      grelcor[nj][ie]->SetMinimum(0.68);
      grelcor[nj][ie]->SetTitle(0);
      grelcor[nj][ie]->GetXaxis()->SetRangeUser(30,1000);      
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
      if(ipad==1)drawText2(ctype,0.27,0.90,15);      
    }
  }

  return 0;
}
void perform_smart_fit(TGraphErrors * gabscor, TF1 * fabscor) {

  int maxFitIter = 50;
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
  //if (fabscor->Integral(0,10) > 25)
  //cout << "\t***ERROR***, fit for histo " << gabscor->GetName() << " diverges at low pt (<10 GeV/c)" << endl;
  
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

