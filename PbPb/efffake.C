#include "utilities.h"

int efffake()
{
  LoadStyle();
  
  const int ncen=6;
  const char *ccent[ncen] = {"0-5%","5-10%","10-30%","30-50%","50-70%","70-90%"};
  
  const int knj=2;
  const char *calgo[knj]={"akVs3PF","akPu3PF"};
  const int iAlg[knj] = {0,2};
  //! Input file PbPb
  TFile *fin_pbpb = NULL;
  fin_pbpb = new TFile("../input/JetResponse_histos_JECv14_embeded_pp_cent.root","r");

  TFile *fin_pp=NULL;
  fin_pp = new TFile("../input/JetResponse_histos_JECv14_pp.root","r");

  //! pp
  TGraphAsymmErrors *gr_fake_pp_pt[knj];
  TGraphAsymmErrors *gr_eff_pp_pt[knj];

  //! Heavy-ion
  TGraphAsymmErrors  *gr_fake_pbpb_pt[knj][ncen];
  TGraphAsymmErrors *gr_eff_pbpb_pt [knj][ncen];

  Color_t icol[ncen] = {kRed, kBlue, kGreen+3, kOrange+3, kMagenta, kBlack};
  int isty[ncen]     = {24  , 25   , 26      , 30       , 32      , 27    };
  TH1F *hnum=0, *hden=0;
  for(int nj=0;nj<knj;nj++){
    
    TH1F *hden = (TH1F*)fin_pp->Get(Form("hPtAll%d_0",iAlg[nj]));
    hden->Rebin(5);
    TH1F *hnum = (TH1F*)fin_pp->Get(Form("hPtSel%d_0",iAlg[nj]));
    hnum->Rebin(5);

    gr_eff_pp_pt[nj] = new TGraphAsymmErrors(hnum, hden,"cl=0.683 b(1,1) mode");
    gr_eff_pp_pt[nj]->SetMarkerStyle(20);
    gr_eff_pp_pt[nj]->SetMarkerColor(1);
    gr_eff_pp_pt[nj]->SetMarkerSize(1.2);

    delete hden;
    delete hnum;

    TH1F *hden = (TH1F*)fin_pp->Get(Form("hFakePtAll%d_0",iAlg[nj]));
    TH1F *hnum = (TH1F*)fin_pp->Get(Form("hFakePtSel%d_0",iAlg[nj]));    

    gr_fake_pp_pt[nj] = new TGraphAsymmErrors(hnum, hden,"cl=0.683 b(1,1) mode");
    gr_fake_pp_pt[nj]->SetMarkerStyle(24);
    gr_fake_pp_pt[nj]->SetMarkerColor(1);
    gr_fake_pp_pt[nj]->SetMarkerSize(1.3);

    delete hden;
    delete hnum;

    for(int ic=0;ic<ncen; ic++){
      TH1F *hden = (TH1F*)fin_pbpb->Get(Form("hPtAll%d_%d",iAlg[nj],ic));
      hden->Rebin(5);
      TH1F *hnum = (TH1F*)fin_pbpb->Get(Form("hPtSel%d_%d",iAlg[nj],ic));
      hnum->Rebin(5);

      gr_eff_pbpb_pt[nj][ic] = new TGraphAsymmErrors(hnum, hden,"cl=0.683 b(1,1) mode");
      gr_eff_pbpb_pt[nj][ic]->SetName(Form("gr_eff_pbpb_pt_%d_%d",nj,ic));
      gr_eff_pbpb_pt[nj][ic]->SetMarkerStyle(isty[ic]);
      gr_eff_pbpb_pt[nj][ic]->SetMarkerColor(icol[ic]);
      gr_eff_pbpb_pt[nj][ic]->SetLineColor(icol[ic]);
      gr_eff_pbpb_pt[nj][ic]->SetMarkerSize(1.6);

      delete hden;
      delete hnum;

      TH1F *hden = (TH1F*)fin_pp->Get(Form("hFakePtAll%d_%d",iAlg[nj],ic));
      TH1F *hnum = (TH1F*)fin_pp->Get(Form("hFakePtSel%d_%d",iAlg[nj],ic));    
      
      gr_fake_pbpb_pt[nj][ic] = new TGraphAsymmErrors(hnum, hden,"cl=0.683 b(1,1) mode");
      gr_fake_pbpb_pt[nj][ic]->SetName(Form("gr_fake_pbpb_pt_%d_%d",nj,ic));
      gr_fake_pbpb_pt[nj][ic]->SetMarkerStyle(isty[ic]);
      gr_fake_pbpb_pt[nj][ic]->SetMarkerColor(icol[ic]);
      gr_fake_pbpb_pt[nj][ic]->SetLineColor(icol[ic]);
      gr_fake_pbpb_pt[nj][ic]->SetMarkerSize(1.6);
      
      delete hden;
      delete hnum;
      
    }
  }//! nj


  TLegend *leg = new TLegend(0.3771486,0.2040302,0.7583418,0.5906801,"","BRNDC");
  leg->SetHeader("");
  leg->SetBorderSize(0);
  leg->SetLineColor(10);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04);
  leg->AddEntry(gr_eff_pp_pt[0],"pp","p");
  leg->AddEntry(gr_eff_pbpb_pt[0][0],ccent[0],"p");
  leg->AddEntry(gr_eff_pbpb_pt[0][1],ccent[1],"p");
  leg->AddEntry(gr_eff_pbpb_pt[0][2],ccent[2],"p");
  leg->AddEntry(gr_eff_pbpb_pt[0][3],ccent[3],"p");
  leg->AddEntry(gr_eff_pbpb_pt[0][4],ccent[4],"p");
  leg->AddEntry(gr_eff_pbpb_pt[0][5],ccent[5],"p");

  TCanvas *c1[knj], *c2[knj];
  for(int nj=0; nj<knj; nj++){
    c1[nj] = new TCanvas(Form("c1_%d",nj),Form("%s Eff ",calgo[nj]),913,822);
    c1[nj]->cd();
    c1[nj]->SetLeftMargin(0.13);
    c1[nj]->SetBottomMargin(0.18);
    gr_eff_pp_pt[nj]->SetMaximum(1.04);
    gr_eff_pp_pt[nj]->SetMinimum(0.56);
    gr_eff_pp_pt[nj]->SetTitle("");
    gr_eff_pp_pt[nj]->GetXaxis()->SetRangeUser(15,150);
    gr_eff_pp_pt[nj]->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
    gr_eff_pp_pt[nj]->GetXaxis()->SetLabelSize(0.05);
    gr_eff_pp_pt[nj]->GetXaxis()->SetTitleSize(0.05);
    gr_eff_pp_pt[nj]->GetXaxis()->CenterTitle();
    gr_eff_pp_pt[nj]->GetYaxis()->SetTitle("Jet reconstruction efficiency");
    gr_eff_pp_pt[nj]->GetYaxis()->CenterTitle();
    gr_eff_pp_pt[nj]->GetYaxis()->SetTitleOffset(1.26);
    gr_eff_pp_pt[nj]->GetYaxis()->SetLabelSize(0.05);
    gr_eff_pp_pt[nj]->GetYaxis()->SetTitleSize(0.05);
    gr_eff_pp_pt[nj]->Draw("ap");
    for(int ic=0; ic<ncen; ic++){
      gr_eff_pbpb_pt[nj][ic]->Draw("psame");
    }
    leg->Draw();
    drawText2(Form("%s |#eta|<2",calgo[nj]),0.45,0.60,26);

    //
    c2[nj] = new TCanvas(Form("c2_%d",nj),Form("%s Fake ",calgo[nj]),913,822);
    c2[nj]->cd();
    c2[nj]->SetLeftMargin(0.13);
    c2[nj]->SetBottomMargin(0.18);
    gr_fake_pp_pt[nj]->SetMaximum(1.04);
    gr_fake_pp_pt[nj]->SetMinimum(0.0);
    gr_fake_pp_pt[nj]->SetTitle("");
    gr_fake_pp_pt[nj]->GetXaxis()->SetRangeUser(15,150);
    gr_fake_pp_pt[nj]->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
    gr_fake_pp_pt[nj]->GetXaxis()->SetLabelSize(0.05);
    gr_fake_pp_pt[nj]->GetXaxis()->SetTitleSize(0.05);
    gr_fake_pp_pt[nj]->GetXaxis()->CenterTitle();
    gr_fake_pp_pt[nj]->GetYaxis()->SetTitle("Fake Rate");
    gr_fake_pp_pt[nj]->GetYaxis()->CenterTitle();
    gr_fake_pp_pt[nj]->GetYaxis()->SetTitleOffset(1.26);
    gr_fake_pp_pt[nj]->GetYaxis()->SetLabelSize(0.05);
    gr_fake_pp_pt[nj]->GetYaxis()->SetTitleSize(0.05);
    gr_fake_pp_pt[nj]->Draw("ap");
    for(int ic=0; ic<ncen; ic++){
      gr_fake_pbpb_pt[nj][ic]->Draw("psame");
    }
    leg->Draw();
    drawText2(Form("%s |#eta|<2",calgo[nj]),0.45,0.60,26);
  }
  return 0; 
}
