#include <iostream>
#include <stdio.h>

#include <TRandom.h>
#include <TH1F.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TProfile.h>
#include <TStopwatch.h>
#include <TCut.h>
#include <cstdlib>
#include <cmath>
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TLine.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"


const Double_t boundaries_jetPtBin[]={0,4,8,14,24,34,44,54,64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 429, 507, 592, 1000};
const int nbins_jetPtBin = 25;

const int ncen=6;
int GetCentBin(int /*hiBin*/);

void pfCand_mc(){
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  // number convension:
  // 1 - 55 or 65
  // 2 - 80 or 95
  // 80 is the unprescaled trigger - yes
  //
  
  //data files - PbPb
  TFile *fpbpb=0;
  const int nfiles=5;
  const char *infile[nfiles] ={
    "/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet30_FOREST_Track8_Jet24_FixedPtHat_v0_mergedpkurt/0.root",
    "/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet50_FOREST_Track8_Jet24_FixedPtHat_v0_mergedpkurt/0.root",
    "/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet80_FOREST_Track8_Jet24_FixedPtHat_v0_mergedpkurt/0.root",
    "/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet100_FOREST_Track8_Jet24_FixedPtHat_v0/0.root",
    "/mnt/hadoop/cms/store/user/velicanu/HydjetDrum_Pyquen_Dijet120_FOREST_Track8_Jet24_FixedPtHat_v0/0.root"
  };

  double xsec[nfiles]     = {1.075e-02,1.025e-03,9.865e-05,3.069e-05,1.129E-05};
  double maxpthat[nfiles] = {50       ,80       ,100      ,120      , 9999    };

  for(int i=0;i<nfiles-1;i++){
    xsec[i] -= xsec[i+1];
  }

  
  //scale the PbPb histograms before adding them
  //we have to scale them according to the lumi of the Jet80 file.
  // HLT file | Lumi inverse micro barns
  // HLT_80 | 149.382
  // HLT_65 | 3.195
  // HLT_55 | 2.734

  //! Add the Jet Trees
  TTree *jetpbpb=0, *evtpbpb=0, *hltpbpb=0, *skmpbpb=0, *pfcpbpb=0;
  static const int nbins = 29;
  static const double boundaries[nbins+1] = {22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 790, 967};

  TFile *fout = new TFile("test_new.root","RECREATE");
  TH1F *hpbpbComb[ncen];
  TH2F *hjbkgComb[ncen];
  
  TH2F *hpfptComb[ncen];
  TH2F *hpfvsptComb[ncen];
  TH2F *hpfvsptiniComb[ncen];
  TH2F *hdiffpfvsptComb[ncen];

  for(int ic=0; ic<ncen; ic++){
    hpbpbComb[ic] = new TH1F(Form("hpbpbComb_%d",ic),"pbpbComb",500,0,1000);
    hjbkgComb[ic] = new TH2F(Form("hjbkgComb_%d",ic),"jbkgComb",500,0,1000,100,-5.0,5.0);
    hpfptComb[ic] = new TH2F(Form("hpfptComb_%d",ic),Form("pfptComb %d",ic),500,0,1000,500,0,100);
    hpfvsptComb[ic] = new TH2F(Form("hpfvsptComb_%d",ic),Form("pfvsptComb %d",ic),500,0,1000,500,0,100);
    hpfvsptiniComb[ic] = new TH2F(Form("hpfvsptiniComb_%d",ic),Form("pfvsptiniComb %d",ic),500,0,1000,100,-10.0,10.0);
    hdiffpfvsptComb[ic] = new TH2F(Form("hdiffpfvsptComb_%d",ic),Form("hdiffpfvspt Comb %d",ic),500,0,1000,500,-100,100);
  }    


  //! Load the jet energy correction factors on fly
  string L2Name, L3Name;
  JetCorrectorParameters *parHI_l2=0, *parHI_l3=0;
  vector<JetCorrectorParameters> vpar_HI;
  FactorizedJetCorrector *_JEC_HI = new FactorizedJetCorrector(vpar_HI);//JR
                                                                                                                         

  //L2Name = "/net/hisrv0001/home/pawan/Validation/Track8_Jet19/combinePtHatBins/ppJEC2014/JECv14/HI_PythiaZ2_2760GeV_5316_v14_L2Relative_"+corrFileName[nj]+"_offline.txt";
  L2Name = "/net/hisrv0001/home/pawan/Validation/Track8_Jet19/combinePtHatBins/ppJEC2014/JECv14/HI_PythiaZ2_2760GeV_5316_v14_L2Relative_AKVs3PF_offline.txt";
  cout<<"**** ++++++  L2Name "<<L2Name<<endl;
  //L3Name = "/net/hisrv0001/home/pawan/Validation/Track8_Jet19/combinePtHatBins/ppJEC2014/JECv14/HI_PythiaZ2_2760GeV_5316_v14_L3Absolute_"+corrFileName[nj]+"_offline.txt";
  L3Name = "/net/hisrv0001/home/pawan/Validation/Track8_Jet19/combinePtHatBins/ppJEC2014/JECv14/HI_PythiaZ2_2760GeV_5316_v14_L3Absolute_AKVs3PF_offline.txt";
  cout<<"**** ++++++  L3Name "<<L3Name<<endl;

  parHI_l2 = new JetCorrectorParameters(L2Name.c_str());
  parHI_l3 = new JetCorrectorParameters(L3Name.c_str());
  
  vpar_HI.push_back(*parHI_l2);
  vpar_HI.push_back(*parHI_l3);
  _JEC_HI = new FactorizedJetCorrector(vpar_HI);


  Long64_t nbytes=0;
  for(int it=0; it<nfiles; it++){

    fpbpb = TFile::Open(infile[it]);
    std::cout<<" Reading file # "<<infile[it]<<std::endl;
    
    jetpbpb = (TTree*)fpbpb->Get("akVs3PFJetAnalyzer/t");
    evtpbpb = (TTree*)fpbpb->Get("hiEvtAnalyzer/HiTree");
    pfcpbpb = (TTree*)fpbpb->Get("pfcandAnalyzer/pfTree");

    jetpbpb->AddFriend(evtpbpb);
    jetpbpb->AddFriend(pfcpbpb);

    //! Event
    Int_t hiBin;
    Float_t vz;
    jetpbpb->SetBranchAddress("hiBin",&hiBin);
    jetpbpb->SetBranchAddress("vz",&vz);

    //! Jet
    Int_t nref;
    Float_t pthat;
    Float_t rawpt[1000];
    Float_t jtpt[1000];
    Float_t jtphi[1000];
    Float_t jteta[1000];
    Float_t chargedMax[1000];
    Float_t chargedSum[1000];
    Float_t photonSum [1000];
    Float_t neutralSum[1000];
    Float_t refdrjt[1000];
    Int_t subid[1000];

    jetpbpb->SetBranchAddress("nref",&nref);
    jetpbpb->SetBranchAddress("pthat",&pthat);
    jetpbpb->SetBranchAddress("rawpt",rawpt);
    jetpbpb->SetBranchAddress("jtpt",jtpt);
    jetpbpb->SetBranchAddress("jteta",jteta);
    jetpbpb->SetBranchAddress("chargedMax",chargedMax);
    jetpbpb->SetBranchAddress("chargedSum",chargedSum);
    jetpbpb->SetBranchAddress("photonSum",photonSum);
    jetpbpb->SetBranchAddress("neutralSum",neutralSum);
    jetpbpb->SetBranchAddress("subid",subid);
    jetpbpb->SetBranchAddress("refdrjt",refdrjt);


    
    //! PF cand
    Int_t nPFpart;
    Int_t pfId[10000];
    Float_t pfPt[10000];
    Float_t pfVsPt[10000];
    Float_t pfVsPtInitial[10000];
    Float_t pfEta[10000];
    Float_t pfPhi[10000];

    jetpbpb->SetBranchAddress("nPFpart",&nPFpart);
    jetpbpb->SetBranchAddress("pfId",pfId);
    jetpbpb->SetBranchAddress("pfPt",pfPt);
    jetpbpb->SetBranchAddress("pfVsPt",pfVsPt);
    jetpbpb->SetBranchAddress("pfVsPtInitial",pfVsPtInitial);
    jetpbpb->SetBranchAddress("pfEta",pfEta);    
    jetpbpb->SetBranchAddress("pfPhi",pfPhi);


    std::cout<<" Added the necessary trees " <<endl;

    double wxs = xsec[it]/jetpbpb->GetEntries()/1000.;
 
    for(int i=0; i<jetpbpb->GetEntries(); i++){
      nbytes += jetpbpb->GetEntry(i);

      if(i%1000==0)std::cout<<" \t  events processed # "<< i <<std::endl;

      bool selEvent = fabs(vz)<15 && pthat < maxpthat[it];
      if(!selEvent)continue;

      int icent = GetCentBin(hiBin);
      if(icent <0 || icent>=ncen)continue;



      double leadEta=-999;
      double leadPhi=-999;
      double leadpT=-999;
      for(int ij=0; ij<nref;ij++){

	if(subid[ij]!=0)continue;

	_JEC_HI->setJetEta(jteta[ij]);
        _JEC_HI->setJetPt (rawpt[ij]);
	
	//double recopt = jtpt[ij];
	double recopt = rawpt[ij]*_JEC_HI->getCorrection();

	bool selJet = fabs(jteta[ij])<2 &&  chargedMax[ij]/recopt > 0.01 && recopt>20. && rawpt[ij]>10; 
	if(!selJet)continue;
	
	double jetbkg = (chargedSum[ij] + neutralSum[ij] + photonSum[ij]) - rawpt[ij];
	hpbpbComb[icent]->Fill(recopt,wxs);

	if(recopt > leadpT){
	  leadpT      = recopt;
	  leadPhi     = jtphi[ij];
	  leadEta     = jteta[ij];
	}
      }
      double sumpt=0;
      //! pf cand loop
      for(int ip=0; ip<nPFpart;ip++){
	double dr = sqrt(pow(leadEta - pfEta[ip],2) + pow(leadPhi - pfPhi[ip],2));
	if(dr>0.3)continue;
	hpfptComb[icent]->Fill(leadpT,pfPt[ip],wxs);
	hpfvsptComb[icent]->Fill(leadpT,pfVsPt[ip],wxs);
	hpfvsptiniComb[icent]->Fill(leadpT,pfVsPtInitial[ip],wxs);
	hdiffpfvsptComb[icent]->Fill(leadpT,pfVsPtInitial[ip] - pfVsPt[ip],wxs);
	sumpt += pfVsPtInitial[ip];
      }
      hjbkgComb[icent]->Fill(leadpT,sumpt,wxs);
    }
  }
  delete parHI_l2;
  delete parHI_l3;
  delete _JEC_HI;

  fout->cd();
  fout->Write();
  fout->Close();
}
int GetCentBin(int bin)
{
  //! PbPb
  int ibin=-1;
  //! centrality bins as 0.5% 0-200
  if(bin<10)ibin=0; //! 0-5%
  else if(bin>=10  && bin<20)ibin=1;   //! 5-10% 
  else if(bin>=20  && bin<60)ibin=2;   //! 10-30%
  else if(bin>=60  && bin<100)ibin=3;  //! 30-50%
  else if(bin>=100 && bin<140)ibin=4;  //! 50-70%
  else if(bin>=140 && bin<180)ibin=5;  //! 70-90%
  return ibin;
}
