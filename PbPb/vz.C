void vz()
{
  TFile *f1 = TFile::Open("/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet50_FOREST_Track8_Jet24_FixedPtHat_v0_mergedpkurt/0.root");
  TTree *tr_mc = (TTree*)f1->Get("hiEvtAnalyzer/HiTree");

  TFile *f2 = TFile::Open("/mnt/hadoop/cms/store/user/velicanu/hiForest_Jet80or95_GR_R_53_LV6_12Mar2014_0000CET_Track8_Jet21/0.root");

  TTree *tr_da = (TTree*)f2->Get("hiEvtAnalyzer/HiTree");
  TTree *skmpbpb = (TTree*)f2->Get("skimanalysis/HltTree");
  tr_da->AddFriend(skmpbpb);

  TH1D *hvz_mc = new TH1D("hvz_mc","vz",200,-20,20);
  hvz_mc->Sumw2();


  TH1D *hvz_data = new TH1D("hvz_data","vz data",200,-20,20);
  hvz_data->Sumw2();

  tr_mc->Project("hvz_mc","vz","TMath::Abs(vz)<15");
  tr_da->Project("hvz_data","vz","TMath::Abs(vz)<15&&pHBHENoiseFilter&&pcollisionEventSelection");

  TFile *fout = new TFile("vz_datamc.root","RECREATE");
  fout->cd();
  hvz_mc->Write();
  hvz_data->Write();
  fout->Close();

}
