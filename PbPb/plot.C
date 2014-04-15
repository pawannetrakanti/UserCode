const int nRad=3;
const int nalg=6;
const char *calg [nRad][nalg] = {{"akVs3Calo","ak3Calo","akPu3Calo","akVs3PF","ak3PF","akPu3PF"},
				 {"akVs4Calo","ak4Calo","akPu4Calo","akVs4PF","ak4PF","akPu4PF"},
				 {"akVs5Calo","ak5Calo","akPu5Calo","akVs5PF","ak5PF","akPu5PF"}
};
int plot(int iRad=0,const char *file="corr")
{

  gStyle->SetOptStat(0);

  TFile *fin = 0;
  if(strcmp(file,"corr")==0)fin = TFile::Open(Form("dijet_pp_mergedpthatbins_Track8_Jet21MC_corrected.root",file));
  else fin  = TFile::Open(Form("dijet_pp_mergedpthatbins_Track8_Jet21MC.root",file));

  TH1D *hpthat = new TH1D("hpthat","#hat{p_{T}} distribution",500,0,1000);
  hpthat->GetXaxis()->SetTitle("#hat{p_{T}} (GeV/c)");
  hpthat->GetXaxis()->CenterTitle();
  hpthat->GetYaxis()->SetTitle("Normalized Yield");
  hpthat->GetYaxis()->CenterTitle();
  hpthat->SetMarkerStyle(20);
  hpthat->Sumw2();

  TH2D *hjes[nalg];
  TProfile *hprof_jes[nalg];
  TTree *tr=0;

  for(int i=0; i<nalg; i++){
    hjes[i] = new TH2D(Form("hjes_%d",i),Form("Jet energy scale %s",calg[iRad][i]),200,0,1000,100,0,2);
    hjes[i]->GetXaxis()->SetTitle("Gen p_{T} (GeV/c)");
    hjes[i]->GetXaxis()->CenterTitle();
    hjes[i]->GetXaxis()->SetRangeUser(20,600);
    hjes[i]->GetYaxis()->SetTitle("Reco p_{T} / Gen p_{T} ");
    hjes[i]->GetYaxis()->CenterTitle();

    hjes[i]->Sumw2();

    tr = (TTree*)fin->Get(Form("%sJetAnalyzer/t",calg[iRad][i]));

    if(i==0)tr->Draw("pthat>>hpthat","weight*(TMath::Abs(refeta)<2. && refpt>10)");

    if(strcmp(file,"corr")==0){
      tr->Draw(Form("corrpt/refpt:refpt>>hjes_%d",i),"weight*(TMath::Abs(refeta)<2. && refpt>10 && subid==0)");
    }
    else tr->Draw(Form("corrpt/refpt:refpt>>hjes_%d",i),"weight*(TMath::Abs(refeta)<2. && refpt>10 && subid==0)");

    hprof_jes[i] = (TProfile*)hjes[i]->ProfileX(Form("hprof_jes_%d",i),1,hjes[i]->GetNbinsX());
    hprof_jes[i]->SetMarkerStyle(20);
  }



  TCanvas *cpThat = new TCanvas("cpThat","pT Hat",639,567);
  cpThat->cd();
  gPad->SetLogy();
  hpthat->Draw("p");
  cpThat->SaveAs("pTHat.gif");

  TCanvas *cJes = new TCanvas("cJes","Jet Energy Scale",1685,894);
  cJes->Divide(3,2);
  for(int i=0; i<nalg; i++){
    cJes->cd(i+1);
    hjes[i]->Draw("colz");
    hprof_jes[i]->Draw("psame");
  }
  cJes->SaveAs(Form("jes_%d_%s.gif",iRad,file));
  return 1;
}
