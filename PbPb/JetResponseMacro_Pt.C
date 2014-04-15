
#include <iostream>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <string.h>

#include <TROOT.h>
#include <TVector3.h> 
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include "utils.h"


using namespace std;

#define PI 3.14159;
const double leadingjetcut = 20. ;
bool SavePlot = true ;


const char* AlgoNames[] = {"ak3PF","ak4PF","ak5PF",
			   "ak3Calo","ak4Calo","ak5Calo",
			   "akPu3PF","akPu4PF","akPu5PF",
			   "akPu3Calo","akPu4Calo","akPu5Calo",
			   "akVs3PF","akVs4PF","akVs5PF",
			   "akVs3Calo","akVs4Calo","akVs5Calo",			   
};
const Int_t nBin = 18 ;

const Double_t jetPtBin[] = {10, 20, 30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,400};
const int nJetPtBin = sizeof(jetPtBin)/sizeof(Double_t)-1 ;

void makeMultiPanelCanvas(TCanvas*& canv, const Int_t columns, 
			  const Int_t rows, const Float_t leftOffset=0.,
                          const Float_t bottomOffset=0., 
			  const Float_t leftMargin=0.2, 
			  const Float_t bottomMargin=0.2,
                          const Float_t edge=0.05);

void JetResponseMacro_Pt(){

  TFile *FileA = TFile::Open(Form("dijet_pp_mergedpthatbins_Track8_Jet21MC_corrected.root"));
  TString outname = "Plots/JetResponse_Plots.root"; 
  TFile* outf = new TFile(outname,"recreate");


  TString canv_name = "c1";
  const Double_t kw = 1300;
  const Double_t kh = 480;
  TCanvas *c1 = new TCanvas(canv_name," ",10,10,kw,kh);
  makeMultiPanelCanvas(c1,nBin/2,2,0.0,0.0,0.1,0.2,0.05);  
  
  gStyle->SetOptStat(0);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadTopMargin   (0.025);
  gStyle->SetPadLeftMargin  (0.15);
  gStyle->SetPadRightMargin (0.025);
  gStyle->SetPadTickX       (1);
  gStyle->SetPadTickY       (1);
  
  
  TLegend *t1=new TLegend(0.2,0.75,0.8,0.92);
  t1->SetFillColor(0);
  t1->SetBorderSize(0);
  t1->SetFillStyle(0);
  t1->SetTextFont(63);
  t1->SetTextSize(15);
  TLegend *t2=new TLegend(0.20,0.45,0.35,0.6);
  t2->SetFillColor(0);
  t2->SetBorderSize(0);
  t2->SetFillStyle(0);
  t2->SetTextFont(63);
  t2->SetTextSize(17);
  
  TH1F * dummy = new TH1F("dummy", "dummy", 500, 0., 500.);
  dummy->SetAxisRange(0., 200., "X") ;
  dummy->GetXaxis()->SetTitle("p_{T}^{ref} (GeV/c)");
  
  const int cCanDiv = nBin;
  int nCanvDiv = cCanDiv;
  const int CnXbins = nJetPtBin;
  int nXbins = CnXbins;
  
  TH2D* iHistoCorrPt[cCanDiv];
  TH2D* iHistoRefPt[cCanDiv];
  TH2D* iHistoJtPt[cCanDiv];
  TH2D* iHistoRawPt[cCanDiv];
  TH1D* h1[cCanDiv][CnXbins];
  
  TF1* Gauss6[cCanDiv][CnXbins];
  double mean[cCanDiv][CnXbins];
  double sigma[cCanDiv][CnXbins];
  double meanErr[cCanDiv][CnXbins];
  double sigmaErr[cCanDiv][CnXbins];
  double xPoint[cCanDiv][CnXbins];
  double xPointErr[cCanDiv][CnXbins];
  TLegend*  leg[cCanDiv];
  
  for(int i =1; i<nCanvDiv; i++){
    c1->cd(i+1);
    
    TTree* t      = (TTree*)FileA->Get(Form("%sJetAnalyzer/t", AlgoNames[i]));
    
    leg[i]= new TLegend(0.22,0.63,0.45,0.90);//top right
    leg[i]->SetFillColor(0);
    leg[i]->SetTextSize(0.05);
    leg[i]->SetBorderSize(0);
    
    iHistoCorrPt[i] = new TH2D(Form("%s_CorrPt",AlgoNames[i]),Form("%s_CorrPt",AlgoNames[i]),nXbins,0,400,500,0,5);
    t->Draw(Form("corrpt/refpt:refpt>>%s_CorrPt",AlgoNames[i]),"weight*((corrpt>0 && TMath::Abs(jteta)<2.0))","goff");
    //t->Draw(Form("rawpt/refpt:refpt>>%s_CorrPt",AlgoNames[i]),Form("rawpt>0 && TMath::Abs(jteta)<2.0"),"goff");
    iHistoCorrPt[i]->Draw("colz");
    leg[i]->AddEntry(iHistoCorrPt[i],Form("%s |#eta|<2.0",AlgoNames[i]),"");
    leg[i]->AddEntry(iHistoCorrPt[i],Form("Jet p_{T}^{corrected}/p_{T}^{ref}"),"");
    format1Dhisto(*iHistoCorrPt[i],250,-1,2,20,2,1,"Jet p_{T}^{ref} (GeV/c)","Corrected p_{T}/ ref p_{T}");//Red      
    leg[i]->Draw();
    for(int iX=iHistoCorrPt[i]->GetXaxis()->GetFirst(); iX<iHistoCorrPt[i]->GetXaxis()->GetLast(); iX++ ){
      //cout<<"bin ( "<<iX<<" ) --limits ["<<iHistoCorrPt[i]->GetXaxis()->GetBinLowEdge(iX)<<" , "<<iHistoCorrPt[i]->GetXaxis()->GetBinUpEdge(iX)<<" ] "<<endl;
      h1[i][iX] = new TH1D(Form("h1_%d_%d_py",i,iX),Form("h1_%d_%d_py",i,iX),500,0,5);
      h1[i][iX] =(TH1D*)iHistoCorrPt[i]->ProjectionY(Form("h1_%d_%d_py",i,iX),iX,iX,"e");
      //h1[i][iX] =(TH1D*)iHistoCorrPt[i]->ProjectionY(Form("h1_%s_Pt%f_%f_py",AlgoNames[i],jetPtBin[iX],jetPtBin[iX+1]),iHistoCorrPt[i]->GetXaxis()->FindBin(jetPtBin[iX]),iHistoCorrPt[i]->GetXaxis()->FindBin(jetPtBin[iX+1])-1,"e");
      
      Gauss6[i][iX] = new TF1(Form("F6_c%d_d%d",i,iX),"gaus",0,3);
      Gauss6[i][iX]->SetParLimits(1,0.7,1.5);
      h1[i][iX]->Fit(Form("F6_c%d_d%d",i,iX),"0Q");
      //   h1[i][iX]->Fit(Gauss6[i][iX],"0Q");
      mean[i][iX] = Gauss6[i][iX]->GetParameter(1);
      sigma[i][iX] = Gauss6[i][iX]->GetParameter(2);
      meanErr[i][iX] = Gauss6[i][iX]->GetParError(1);
      sigmaErr[i][iX] = Gauss6[i][iX]->GetParError(2);
      xPoint[i][iX] =iHistoCorrPt[i]->GetXaxis()->GetBinCenter(iX);
      xPointErr[i][iX] = 0;
      
      //       if(i==2 && iX==3){//Just to test 
      // 	cout<<"bin ( "<<iX<<" ) --limits ["<<iHistoCorrPt[i]->GetXaxis()->GetBinLowEdge(iX)<<" , "<<iHistoCorrPt[i]->GetXaxis()->GetBinUpEdge(iX)<<" ] "<<endl;
      // 	TCanvas *c11 = new TCanvas("c11","c11",500,400);
      // 	c11->cd();	    
      // 	h1[i][iX]->Draw();
      // 	Gauss6[i][iX]->Draw("same");
      //      }	  
    }
  }    
  
  if(SavePlot){
    c1->Print("Plots/CorrectedOverRefJet.gif");
  }
  
  c1->Update();
  
  TLegend*  leg2B[cCanDiv];
  TLegend*  leg2A[cCanDiv];
  TGraphErrors *CorrPt_Mean[cCanDiv];
  TCanvas *c2 = new TCanvas("c2","c2",10,10,kw,kh);
  makeMultiPanelCanvas(c2,nBin/2,2,0.0,0.0,0.1,0.2,0.05);  
  
  TGraphErrors *CorrPt_Sigma[cCanDiv];
  TCanvas *c3 = new TCanvas("c3","c3",10,10,kw,kh);
  makeMultiPanelCanvas(c3,nBin/2,2,0.0,0.0,0.1,0.2,0.05);  
  
  for(int i2 =1; i2<nCanvDiv; i2++){
    leg2A[i2]= new TLegend(0.22,0.73,0.45,0.95);//top right
    leg2A[i2]->SetFillColor(0);
    leg2A[i2]->SetTextSize(0.05);
    leg2A[i2]->SetBorderSize(0);
    
    leg2B[i2]= new TLegend(0.22,0.73,0.45,0.95);//top right
    leg2B[i2]->SetFillColor(0);
    leg2B[i2]->SetTextSize(0.05);
    leg2B[i2]->SetBorderSize(0);
    
    double xLoc[CnXbins];
    double y2Loc[CnXbins];
    double y3Loc[CnXbins];
    double xLocErr[CnXbins];
    double y2LocErr[CnXbins];
    double y3LocErr[CnXbins];
    for(int iBin=0; iBin<nXbins; iBin++){	  
      xLoc[iBin] = xPoint[i2][iBin];
      xLocErr[iBin] = xPointErr[i2][iBin];
      y2Loc[iBin] = mean[i2][iBin];
      y2LocErr[iBin] = meanErr[i2][iBin];
      y3Loc[iBin] = sigma[i2][iBin];
      y3LocErr[iBin] = sigmaErr[i2][iBin];
    }
    
    c2->cd(i2+1);
    CorrPt_Mean[i2] = new TGraphErrors(nXbins,xLoc,y2Loc,xLocErr,y2LocErr);
    formatTGraph(*CorrPt_Mean[i2],250,-1,1,20,1,1,"p_{T}^{ref} (GeV/c)","#mu ");
    leg2A[i2]->AddEntry(CorrPt_Mean[i2],Form("%s |#eta|<2.0",AlgoNames[i2]),"");
    leg2A[i2]->AddEntry(CorrPt_Mean[i2],Form("p_{T}^{corrected}/p_{T}^{ref}"),"lp");
    CorrPt_Mean[i2]->SetMaximum(1.2);
    CorrPt_Mean[i2]->SetMinimum(0.8);
    CorrPt_Mean[i2]->Draw("AP");
    leg2A[i2]->Draw();
    
    c3->cd(i2+1);
    CorrPt_Sigma[i2] = new TGraphErrors(nXbins,xLoc,y3Loc,xLocErr,y3LocErr);
    formatTGraph(*CorrPt_Sigma[i2],250,-1,1,20,1,1,"p_{T}^{ref} (GeV/c)","#sigma ");
    leg2B[i2]->AddEntry(CorrPt_Sigma[i2],Form("%s |#eta|<2.0",AlgoNames[i2]),"");
    leg2B[i2]->AddEntry(CorrPt_Sigma[i2],Form("p_{T}^{corrected}/p_{T}^{ref}"),"lp");
    CorrPt_Sigma[i2]->SetMaximum(0.4);
    CorrPt_Sigma[i2]->SetMinimum(0.);
    CorrPt_Sigma[i2]->Draw("AP");
    leg2B[i2]->Draw();
    
  }
  if(SavePlot){
    c2->Print("Plots/CorrectedOverRefJetMean.gif");  
    c3->Print("Plots/CorrectedOverRefJetSigma.gif") ; 
  }
  c2->Update();
  c3->Update();

  outf->cd();
  outf->Write();
  outf->Close();
}




///////////////////////////////////////////////////////////////////////
//         TOOL BOX
//////////////////////////////////////////////////////////////////////
void makeMultiPanelCanvas(TCanvas*& canv,
                          const Int_t columns,
                          const Int_t rows,
                          const Float_t leftOffset,
                          const Float_t bottomOffset,
                          const Float_t leftMargin,
                          const Float_t bottomMargin,
                          const Float_t edge) {
   if (canv==0) {
      Error("makeMultiPanelCanvas","Got null canvas.");
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
   Xlow[0] = leftOffset;
   Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
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
void drawText(const char *text, float xp, float yp){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(63);
  tex->SetTextSize(17);
  //tex->SetTextSize(0.05);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw();
}
