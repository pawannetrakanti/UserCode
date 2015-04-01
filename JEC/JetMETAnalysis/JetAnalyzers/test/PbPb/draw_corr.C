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


double ptbins[] = {10,15,20,27,35,45,57,72,90,120,150,200,300,400,550,750,1000}; 
const int npt = sizeof(ptbins)/sizeof(double) - 1;
const int knj = 7;
const char *calgo[knj]={"ak1","ak2","ak3","ak4","ak5","ak6","ak7"};
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

TLegend *getLegend(double /*x1*/, double /*y1*/, double /*x2*/, double /*y2*/);

int draw_corr(const char *ctype="PF")
{

  int iSave=0;

  LoadStyle();
  //gStyle->SetOptStat(1);
  //gStyle->SetOptFit(1);

  cout << "npt : " << npt << "  " << "neta : " << neta <<endl;
  //return 0;

  double refmin=15;
  double refmax=1000;

  TFile *fin_l3 = new TFile("jra_hiF_ak_l3_dijet_final_lowpt.root","r");
  TFile *fin_l2 = new TFile("jra_hiF_ak_l2_dijet_final_lowpt.root","r");

  //TFile *fin_l3    = new TFile("jra_hiF_ak_l3_dijet_test_lowpt.root","r");
  //TFile *fin_l2 = new TFile("jra_hiF_ak_l2_dijet_test_lowpt.root","r");

  //TFile *fin_l3    = new TFile("jra_hiF_ak_l3_dijet_official_lowpt.root","r");
  //TFile *fin_l2 = new TFile("jra_hiF_ak_l2_dijet_official_lowpt.root","r");

  TGraphErrors *grsp_l3[knj], *gcor_l3[knj];
  TGraphErrors *grelcor_l2[knj][neta], *gabscor_l2[knj][neta];

  //TLegend *l0 = getLegend(0.6349664,0.3265565,0.8491835,0.7128335);
  //l0->SetHeader(Form("Anti-k_{T}, %s",ctype));
  //TLegend *l1 = getLegend(0.1500646,0.5052258,0.3653091,0.7490942);
  //l1->SetHeader("");


  for(int i=0; i<knj; i++){
     std::ostringstream cnj;
     if(strcmp(ctype,"PF")==0)cnj << "ak" << i+1 << "PFJetAnalyzer/";
     else if(strcmp(ctype,"Calo"  )==0)cnj << "ak"   << i+1 << "CaloJetAnalyzer/";
     else if(strcmp(ctype,"PuPF"  )==0)cnj << "akPu" << i+1 << "PFJetAnalyzer/";
     else if(strcmp(ctype,"PuCalo")==0)cnj << "akPu" << i+1 << "CaloJetAnalyzer/";
     else if(strcmp(ctype,"VsPF"  )==0)cnj << "akVs" << i+1 << "PFJetAnalyzer/";
     else if(strcmp(ctype,"VsCalo")==0)cnj << "akVs" << i+1 << "CaloJetAnalyzer/";
     
     grsp_l3[i] = (TGraphErrors*)fin_l3->Get(Form("%sL3RspVsRefPt",cnj.str().c_str()));
     grsp_l3[i]->GetXaxis()->SetTitle("< Ref p_{T} > (GeV/c)");

     grsp_l3[i]->GetXaxis()->SetTitleFont(42);
     grsp_l3[i]->GetXaxis()->SetLabelFont(42);
     grsp_l3[i]->GetXaxis()->SetLabelSize(0.06);
     grsp_l3[i]->GetXaxis()->SetTitleSize(0.06);
     grsp_l3[i]->GetXaxis()->SetNdivisions(507);
     grsp_l3[i]->GetYaxis()->SetTitle("L3 Response");
     grsp_l3[i]->GetYaxis()->SetTitleFont(42);
     grsp_l3[i]->GetYaxis()->SetLabelFont(42);
     grsp_l3[i]->GetYaxis()->SetLabelSize(0.06);
     grsp_l3[i]->GetYaxis()->SetTitleSize(0.06);
     grsp_l3[i]->GetYaxis()->SetNdivisions(507);

     gcor_l3[i] = (TGraphErrors*)fin_l3->Get(Form("%sL3CorVsJetPt",cnj.str().c_str()));
     gcor_l3[i]->GetXaxis()->SetTitle("< raw jet p_{T} > (GeV/c)");
     //gcor_l3[i]->GetXaxis()->SetRangeUser(refmin,refmax);
     gcor_l3[i]->GetXaxis()->SetTitleFont(42);
     gcor_l3[i]->GetXaxis()->SetLabelFont(42);
     gcor_l3[i]->GetXaxis()->SetLabelSize(0.06);
     gcor_l3[i]->GetXaxis()->SetTitleSize(0.06);
     gcor_l3[i]->GetXaxis()->SetNdivisions(507);
     gcor_l3[i]->GetXaxis()->SetNoExponent();
     gcor_l3[i]->GetXaxis()->SetMoreLogLabels();
     gcor_l3[i]->GetYaxis()->SetTitle("L3 Correction");
     gcor_l3[i]->GetYaxis()->SetTitleFont(42);
     gcor_l3[i]->GetYaxis()->SetLabelFont(42);
     gcor_l3[i]->GetYaxis()->SetLabelSize(0.06);
     gcor_l3[i]->GetYaxis()->SetTitleSize(0.06);
     gcor_l3[i]->GetYaxis()->SetNdivisions(507);

     
     for(int j=0; j<neta; j++){
       std::ostringstream seta; 
       seta << etabins[j] << "to" << etabins[j+1];
       
       std::string sone(cnj.str()+"AbsCorVsJetPt_JetEta"+seta.str());
       std::string stwo(cnj.str()+"RelCorVsJetPt_JetEta"+seta.str());
       
       gabscor_l2[i][j] = (TGraphErrors*)fin_l2->Get(Form("%s",sone.c_str()));
       gabscor_l2[i][j]->GetXaxis()->SetRangeUser(refmin,refmax);      
       gabscor_l2[i][j]->GetXaxis()->SetTitle("< raw p_{T} > (GeV/c)");
       gabscor_l2[i][j]->GetXaxis()->SetTitleFont(42);
       gabscor_l2[i][j]->GetXaxis()->SetLabelFont(42);
       gabscor_l2[i][j]->GetXaxis()->SetLabelSize(0.08);
       gabscor_l2[i][j]->GetXaxis()->SetTitleSize(0.07);
       gabscor_l2[i][j]->GetXaxis()->SetNdivisions(507);
       gabscor_l2[i][j]->GetXaxis()->SetNoExponent();
       gabscor_l2[i][j]->GetXaxis()->SetMoreLogLabels();
       gabscor_l2[i][j]->GetYaxis()->SetTitle("Abs. L2 Corr.");
       gabscor_l2[i][j]->GetYaxis()->SetTitleFont(42);
       gabscor_l2[i][j]->GetYaxis()->SetLabelFont(42);
       gabscor_l2[i][j]->GetYaxis()->SetLabelSize(0.08);
       gabscor_l2[i][j]->GetYaxis()->SetTitleSize(0.07);
       gabscor_l2[i][j]->GetYaxis()->SetNdivisions(507);



       grelcor_l2[i][j] = (TGraphErrors*)fin_l2->Get(Form("%s",stwo.c_str()));
       grelcor_l2[i][j]->SetTitle(0);
       //grelcor_l2[i][j]->GetXaxis()->SetRangeUser(refmin,refmax);      
       grelcor_l2[i][j]->GetXaxis()->SetTitle("< raw p_{T} > (GeV/c)");
       grelcor_l2[i][j]->GetXaxis()->SetTitleFont(42);
       grelcor_l2[i][j]->GetXaxis()->SetLabelFont(42);
       grelcor_l2[i][j]->GetXaxis()->SetLabelSize(0.08);
       grelcor_l2[i][j]->GetXaxis()->SetTitleSize(0.07);
       grelcor_l2[i][j]->GetXaxis()->SetNdivisions(507);
       grelcor_l2[i][j]->GetXaxis()->SetNoExponent();
       grelcor_l2[i][j]->GetXaxis()->SetMoreLogLabels();
       grelcor_l2[i][j]->GetYaxis()->SetTitle("Rel. L2 Corr.");
       grelcor_l2[i][j]->GetYaxis()->SetTitleFont(42);
       grelcor_l2[i][j]->GetYaxis()->SetLabelFont(42);
       grelcor_l2[i][j]->GetYaxis()->SetLabelSize(0.08);
       grelcor_l2[i][j]->GetYaxis()->SetTitleSize(0.07);
       grelcor_l2[i][j]->GetYaxis()->SetNdivisions(507);

     }
     

  }

  TCanvas *c11 = new TCanvas("c11","L3 Response",1531,683);
  c11->Divide(4,2);
  TCanvas *c12 = new TCanvas("c12","L3 Corrections",1531,683);
  c12->Divide(4,2);
  for(int i=0; i<knj; i++){  
    c11->cd(i+1);
    gPad->SetLogx();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);

    grsp_l3[i]->GetXaxis()->SetNoExponent();
    grsp_l3[i]->GetXaxis()->SetMoreLogLabels();
    grsp_l3[i]->GetXaxis()->SetRangeUser(refmin,refmax);
    grsp_l3[i]->Draw("ap");
    gPad->Update();

    drawText2(Form("%s%s",calgo[i],ctype),0.20,0.78,21);

    // TPaveStats *strsp = (TPaveStats*)grsp_l3[i]->FindObject("stats");
    // strsp->SetX1NDC(0.46);
    // strsp->SetY1NDC(0.18);
    // strsp->SetX2NDC(0.87);
    // strsp->SetY2NDC(0.59);
    // strsp->Draw();
    
    c12->cd(i+1);
    gPad->SetLogx();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->Update();
    //gcor_l3[i]->GetXaxis()->SetRangeUser(15,1000);
    gcor_l3[i]->GetXaxis()->SetNoExponent();
    gcor_l3[i]->GetXaxis()->SetMoreLogLabels();
    gcor_l3[i]->Draw("ap");

    drawText2(Form("%s%s",calgo[i],ctype),0.60,0.78,21);
    
    // TPaveStats *stcor = (TPaveStats*)gcor[i]->FindObject("stats");
    // stcor->SetX1NDC(0.48);
    // stcor->SetY1NDC(0.53);
    // stcor->SetX2NDC(0.89);
    // stcor->SetY2NDC(0.89);
    // stcor->Draw();
  }

  c11->cd(8);
  drawText2("PYTHIA pp Signal",0.17,0.80,21);
  drawText2("#sqrt{s_{NN}} = 2.76 TeV",0.17,0.70,21);    
  drawText2("15 < #hat{p_{T}} (GeV/c) < 540",0.17,0.60,21);    
  drawText2("| #eta^{Ref}_{Jet} | < 3.0",0.17,0.50,21);    
  drawText2("CMSSW_5_3_16" ,0.17,0.40,21);    
  drawText2("STARTHI53_LV1",0.17,0.30,21);    
  drawText2("Track8_Jet29" ,0.17,0.20,21);    
  if(iSave){
    c11->SaveAs(Form("FinalPlots/L3Response_%s.gif",ctype));
    c11->SaveAs(Form("FinalPlots/L3Response_%s.pdf",ctype));
  }

  c12->cd(8);
  drawText2("PYTHIA",0.17,0.80,21);
  drawText2("#sqrt{s_{NN}} = 2.76 TeV",0.17,0.70,21);    
  drawText2("15 < #hat{p_{T}} (GeV/c) < 540",0.17,0.60,21);    
  drawText2("| #eta^{Ref}_{Jet} | < 3.0",0.17,0.50,21);    
  drawText2("CMSSW_5_3_16" ,0.17,0.40,21);    
  drawText2("STARTHI53_LV1",0.17,0.30,21);    
  drawText2("Track8_Jet29" ,0.17,0.20,21);    
  if(iSave){
    c12->SaveAs(Form("FinalPlots/L3Correction_%s.gif",ctype));
    c12->SaveAs(Form("FinalPlots/L3Correction_%s.pdf",ctype));
  }


  int ipad=0;
  //! Absolute L2 corrections
  TCanvas *c99[knj], *c98[knj];
  for(int nj=0; nj<knj; nj++){
    c99[nj] = new TCanvas(Form("c99_%d",nj),Form("%s Abs L2 Fitting plots",calgo[nj]),1854,866);
    c99[nj]->Divide(6,4,0,0);
    ipad=0;
    for(int ie=0;ie<neta;ie++){      
      c99[nj]->cd(++ipad);
      gPad->SetLogx();
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.01);
      gPad->SetBottomMargin(0.15);
      
      gabscor_l2[nj][ie]->SetMaximum(2.14);
      gabscor_l2[nj][ie]->SetMinimum(0.64);
      gabscor_l2[nj][ie]->SetTitle(0);
      gabscor_l2[nj][ie]->GetXaxis()->SetRangeUser(refmin,refmax);      
      gabscor_l2[nj][ie]->GetXaxis()->SetTitle("< raw p_{T} > (GeV/c)");
      gabscor_l2[nj][ie]->GetXaxis()->SetTitleFont(42);
      gabscor_l2[nj][ie]->GetXaxis()->SetLabelFont(42);
      gabscor_l2[nj][ie]->GetXaxis()->SetLabelSize(0.08);
      gabscor_l2[nj][ie]->GetXaxis()->SetTitleSize(0.07);
      gabscor_l2[nj][ie]->GetXaxis()->SetNdivisions(507);
      gabscor_l2[nj][ie]->GetXaxis()->SetNoExponent();
      gabscor_l2[nj][ie]->GetXaxis()->SetMoreLogLabels();
      gabscor_l2[nj][ie]->GetYaxis()->SetTitle("Abs L2 Corr.");
      gabscor_l2[nj][ie]->GetYaxis()->SetTitleFont(42);
      gabscor_l2[nj][ie]->GetYaxis()->SetLabelFont(42);
      gabscor_l2[nj][ie]->GetYaxis()->SetLabelSize(0.08);
      gabscor_l2[nj][ie]->GetYaxis()->SetTitleSize(0.07);
      gabscor_l2[nj][ie]->GetYaxis()->SetNdivisions(507);

      gabscor_l2[nj][ie]->Draw("ap");  

      std::ostringstream seta; 
      seta << etabins[ie] << " < #eta < " << etabins[ie+1];
      drawText2(seta.str().c_str(),0.27,0.80,15);
      if(ipad==1)drawText2(Form("%s%s",calgo[nj],ctype),0.20,0.87,21);
    }

    c99[nj]->cd(23);
    drawText2("PYTHIA pp Signal",0.17,0.80,21);
    drawText2("#sqrt{s_{NN}} = 2.76 TeV",0.17,0.60,21);    
    drawText2("15 < #hat{p_{T}} (GeV/c) < 540",0.17,0.40,21);    
    drawText2("| #eta^{Ref}_{Jet} | < 3.0",0.17,0.20,21);    


    c99[nj]->cd(24);
    drawText2("CMSSW_5_3_16" ,0.17,0.70,21);    
    drawText2("STARTHI53_LV1",0.17,0.50,21);    
    drawText2("Track8_Jet29" ,0.17,0.30,21);    

    if(iSave){
      c99[nj]->SaveAs(Form("FinalPlots/L2AbsCorr_%s%s.gif",calgo[nj],ctype));
      c99[nj]->SaveAs(Form("FinalPlots/L2AbsCorr_%s%s.pdf",calgo[nj],ctype));
    }

    //c99[nj]->Close();
  }
  //return 0;
  //! Relative L2 Correction
  ipad=0;
  for(int nj=0; nj<knj; nj++){
    c98[nj] = new TCanvas(Form("c98_%d",nj),Form("%s Rel L2 Fitting plots",calgo[nj]),1854,866);
    c98[nj]->Divide(6,4,0,0);
    ipad=0;
    for(int ie=0;ie<neta;ie++){      
      
      c98[nj]->cd(++ipad);
      gPad->SetLogx();
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.01);
      gPad->SetBottomMargin(0.15);
      
      grelcor_l2[nj][ie]->SetMaximum(1.38);
      grelcor_l2[nj][ie]->SetMinimum(0.68);
      grelcor_l2[nj][ie]->SetTitle(0);
      grelcor_l2[nj][ie]->GetXaxis()->SetRangeUser(20,1000);      
      grelcor_l2[nj][ie]->GetXaxis()->SetTitle("< raw p_{T} > (GeV/c)");
      grelcor_l2[nj][ie]->GetXaxis()->SetTitleFont(42);
      grelcor_l2[nj][ie]->GetXaxis()->SetLabelFont(42);
      grelcor_l2[nj][ie]->GetXaxis()->SetLabelSize(0.08);
      grelcor_l2[nj][ie]->GetXaxis()->SetTitleSize(0.07);
      grelcor_l2[nj][ie]->GetXaxis()->SetNdivisions(507);
      grelcor_l2[nj][ie]->GetXaxis()->SetNoExponent();
      grelcor_l2[nj][ie]->GetXaxis()->SetMoreLogLabels();
      grelcor_l2[nj][ie]->GetYaxis()->SetTitle("Rel L2 Corr.");
      grelcor_l2[nj][ie]->GetYaxis()->SetTitleFont(42);
      grelcor_l2[nj][ie]->GetYaxis()->SetLabelFont(42);
      grelcor_l2[nj][ie]->GetYaxis()->SetLabelSize(0.08);
      grelcor_l2[nj][ie]->GetYaxis()->SetTitleSize(0.07);
      grelcor_l2[nj][ie]->GetYaxis()->SetNdivisions(507);

      grelcor_l2[nj][ie]->SetMarkerStyle(20);
      grelcor_l2[nj][ie]->Draw("ap");  
      std::ostringstream seta; 
      seta << etabins[ie] << " < #eta < " << etabins[ie+1];
      drawText2(seta.str().c_str(),0.27,0.80,15);
      if(ipad==1)drawText2(Form("%s%s",calgo[nj],ctype),0.20,0.87,21);
    }
    c98[nj]->cd(23);
    drawText2("PYTHIA pp Signal",0.17,0.80,21);
    drawText2("#sqrt{s_{NN}} = 2.76 TeV",0.17,0.60,21);    
    drawText2("15 < #hat{p_{T}} (GeV/c) < 540",0.17,0.40,21);    
    drawText2("| #eta^{Ref}_{Jet} | < 3.0",0.17,0.20,21);    
    c98[nj]->cd(24);
    drawText2("CMSSW_5_3_16" ,0.17,0.70,21);    
    drawText2("STARTHI53_LV1",0.17,0.50,21);    
    drawText2("Track8_Jet29" ,0.17,0.30,21);    

    if(iSave){
      c98[nj]->SaveAs(Form("FinalPlots/L2RelativeCorr_%s%s.gif",calgo[nj],ctype));
      c98[nj]->SaveAs(Form("FinalPlots/L2RelativeCorr_%s%s.pdf",calgo[nj],ctype));
    }
  }

 return 0;
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
