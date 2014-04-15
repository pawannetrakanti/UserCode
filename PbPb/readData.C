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

void readData(){
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();

  // number convension:
  // 1 - 55 or 65
  // 2 - 80 or 95
  // 80 is the unprescaled trigger - yes
  //
  
  //data files - PbPb
  TFile *fpbpb=0;
  const char *infile[2] ={"/mnt/hadoop/cms/store/user/velicanu/hiForest_Jet55or65_GR_R_53_LV6_12Mar2014_0000CET_Track8_Jet21/0.root",
			  "/mnt/hadoop/cms/store/user/velicanu/hiForest_Jet80or95_GR_R_53_LV6_12Mar2014_0000CET_Track8_Jet21/0.root"};
  

  
  //scale the PbPb histograms before adding them
  //we have to scale them according to the lumi of the Jet80 file.
  // HLT file | Lumi inverse micro barns
  // HLT_80 | 149.382
  // HLT_65 | 3.195
  // HLT_55 | 2.734

  //! Add the Jet Trees
  TTree *jetpbpb=0, *evtpbpb=0, *hltpbpb=0, *skmpbpb=0;
  static const int nbins = 29;
  static const double boundaries[nbins+1] = {22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 790, 967};

  TFile *fout = new TFile("output_PbPb_data.root","RECREATE");
  TH1F *hBin    = new TH1F("hBin","Centrality bin",200,-0.5,200-0.5);
  TH1F *hpbpb[ncen][3], *hpbpbComb[ncen];
  TH2F *hjbkg[ncen][3], *hjbkgComb[ncen];

  TH2F *hchargedMax[ncen][3], *hchargedMaxComb[ncen];
  TH2F *hphotonMax[ncen][3], *hphotonMaxComb[ncen];
  TH2F *hneutralMax[ncen][3], *hneutralMaxComb[ncen];

  TH2F *hchargedMaxSp[ncen][3], *hchargedMaxSpComb[ncen];
  TH2F *hphotonMaxSp[ncen][3], *hphotonMaxSpComb[ncen];
  TH2F *hneutralMaxSp[ncen][3], *hneutralMaxSpComb[ncen];

  TH3F *hjbkgvz [ncen][3], *hjbkgvzComb[ncen];
  TH3F *hjbkgeta[ncen][3], *hjbkgetaComb[ncen];
  TH3F *hjbkgraw[ncen][3], *hjbkgrawComb[ncen];

  for(int ic=0; ic<ncen; ic++){
    for(int in=0;in<3;in++){
      hpbpb[ic][in] = new TH1F(Form("hpbpb%d_%d",ic,in),Form("pbpb %d %d",ic,in),100,0,1000);
      hjbkg[ic][in] = new TH2F(Form("hjbkg%d_%d",ic,in),Form("jbkg %d %d",ic,in),100,0,1000,500,0.,250.);
      hchargedMax[ic][in] = new TH2F(Form("hchargedMax%d_%d",ic,in),Form("chargedMax/jtpt %d %d",ic,in),100,0,1000,250,0.,5.);      
      hchargedMaxSp[ic][in] = new TH2F(Form("hchargedMaxSp%d_%d",ic,in),Form("chargedMaxSp/jtpt %d %d",ic,in),100,0,1000,250,0.,5.);      
      hphotonMax[ic][in] = new TH2F(Form("hphotonMax%d_%d",ic,in),Form("photonMax/jtpt %d %d",ic,in),100,0,1000,250,0.,5.);      
      hphotonMaxSp[ic][in] = new TH2F(Form("hphotonMaxSp%d_%d",ic,in),Form("photonMaxSp/jtpt %d %d",ic,in),100,0,1000,250,0.,5.);      
      hneutralMax[ic][in] = new TH2F(Form("hneutralMax%d_%d",ic,in),Form("neutralMax/jtpt %d %d",ic,in),100,0,1000,250,0.,5.);      
      hneutralMaxSp[ic][in] = new TH2F(Form("hneutralMaxSp%d_%d",ic,in),Form("neutralMaxSp/jtpt %d %d",ic,in),100,0,1000,250,0.,5.);      

      hjbkgvz[ic][in]  = new TH3F(Form("hjbkgvz%d_%d",ic,in),Form("vz jbkg %d %d",ic,in),6,-15.,15,100,0,1000,500,0.,250.);
      hjbkgeta[ic][in] = new TH3F(Form("hjbkgeta%d_%d",ic,in),Form("eta jbkg %d %d",ic,in),8,-2.0,2.0,100,0,1000,500,0.,250.);
      hjbkgraw[ic][in] = new TH3F(Form("hjbkgraw%d_%d",ic,in),Form("raw pt jbkg %d %d",ic,in),100,0,1000,100,0,1000,500,0.,250.);
    }
    hpbpbComb[ic] = new TH1F(Form("hpbpbComb_%d",ic),"pbpbComb",100,0,1000);
    hjbkgComb[ic] = new TH2F(Form("hjbkgComb_%d",ic),"jbkgComb",100,0,1000,500,0.,250);
    hchargedMaxComb[ic] = new TH2F(Form("hchargedMaxComb_%d",ic),Form("chargedMax/jtpt %d",ic),100,0,1000,250,0.,5.);      
    hchargedMaxSpComb[ic] = new TH2F(Form("hchargedMaxSpComb_%d",ic),Form("chargedMaxSp/jtpt %d",ic),100,0,1000,250,0.,5.);      
    hphotonMaxComb[ic] = new TH2F(Form("hphotonMaxComb_%d",ic),Form("photonMax/jtpt %d ",ic),100,0,1000,250,0.,5.);      
    hphotonMaxSpComb[ic] = new TH2F(Form("hphotonMaxSpComb_%d",ic),Form("photonMaxSp/jtpt %d",ic),100,0,1000,250,0.,5.);      
    hneutralMaxComb[ic] = new TH2F(Form("hneutralMaxComb_%d",ic),Form("neutralMax/jtpt %d ",ic),100,0,1000,250,0.,5.);      
    hneutralMaxSpComb[ic] = new TH2F(Form("hneutralMaxSpcomb_%d",ic),Form("neutralMaxSp/jtpt %d",ic),100,0,1000,250,0.,5.);      

    hjbkgvzComb[ic]  = new TH3F(Form("hjbkgvzComb_%d",ic),Form("vz jbkg %d",ic),6,-15.,15, 100,0,1000, 500,0.,250.);
    hjbkgetaComb[ic] = new TH3F(Form("hjbkgetaComb_%d",ic),Form("eta jbkg %d",ic),8,-2.0,2.0, 100,0,1000, 500,0.,250.);
    hjbkgrawComb[ic] = new TH3F(Form("hjbkgrawComb_%d",ic),Form("raw pt jbkg %d ",ic),100,0,1000,100,0,1000,500,0.,250.);
  }    


  //! Load the jet energy correction factors on fly
  string L2Name, L3Name;
  JetCorrectorParameters *parHI_l2=0, *parHI_l3=0;
  vector<JetCorrectorParameters> vpar_HI;
  FactorizedJetCorrector *_JEC_HI = new FactorizedJetCorrector(vpar_HI);//JR
                                                                                                                         

  //L2Name = "/net/hisrv0001/home/pawan/Validation/Track8_Jet19/combinePtHatBins/ppJEC2014/JECv14/HI_PythiaZ2_2760GeV_5316_v14_L2Relative_"+corrFileName[nj]+"_offline.txt";
  L2Name = "/net/hisrv0001/home/pawan/Validation/Track8_Jet19/combinePtHatBins/ppJEC2014/JECv14/HI_PythiaZ2_2760GeV_5316_v14_L2Relative_AKPu3PF_offline.txt";
  cout<<"**** ++++++  L2Name "<<L2Name<<endl;
  //L3Name = "/net/hisrv0001/home/pawan/Validation/Track8_Jet19/combinePtHatBins/ppJEC2014/JECv14/HI_PythiaZ2_2760GeV_5316_v14_L3Absolute_"+corrFileName[nj]+"_offline.txt";
  L3Name = "/net/hisrv0001/home/pawan/Validation/Track8_Jet19/combinePtHatBins/ppJEC2014/JECv14/HI_PythiaZ2_2760GeV_5316_v14_L3Absolute_AKPu3PF_offline.txt";
  cout<<"**** ++++++  L3Name "<<L3Name<<endl;

  parHI_l2 = new JetCorrectorParameters(L2Name.c_str());
  parHI_l3 = new JetCorrectorParameters(L3Name.c_str());
  
  vpar_HI.push_back(*parHI_l2);
  vpar_HI.push_back(*parHI_l3);
  _JEC_HI = new FactorizedJetCorrector(vpar_HI);


  Long64_t nbytes=0;
  for(int it=0; it<2; it++){

    fpbpb = TFile::Open(infile[it]);
    std::cout<<" Reading file # "<<infile[it]<<std::endl;
    
    jetpbpb = (TTree*)fpbpb->Get("akPu3PFJetAnalyzer/t");
    evtpbpb = (TTree*)fpbpb->Get("hiEvtAnalyzer/HiTree");
    hltpbpb = (TTree*)fpbpb->Get("hltanalysis/HltTree");
    skmpbpb = (TTree*)fpbpb->Get("skimanalysis/HltTree");


    jetpbpb->AddFriend(evtpbpb);
    jetpbpb->AddFriend(hltpbpb);
    jetpbpb->AddFriend(skmpbpb);

    //jetpbpb->SetAlias("jetbkg","(chargedSum+neutralSum+photonSum) - rawpt");

    //! HLT
    Int_t HLT_HIJet80_v1;
    Int_t HLT_HIJet65_v1;
    Int_t HLT_HIJet55_v1;
    Int_t HLT_HIJet55_v1_Prescl;

    //! Skim
    Int_t pHBHENoiseFilter;
    Int_t pcollisionEventSelection;
    
    //! Event
    Int_t run;
    Int_t evt;
    Int_t hiBin;
    Float_t vz;
    
    Int_t nref;
    Float_t rawpt[1000];
    Float_t jtpt[1000];
    Float_t jteta[1000];
    Float_t chargedMax[1000];
    Float_t chargedSum[1000];
    Float_t photonSum [1000];
    Float_t photonMax [1000];
    Float_t neutralSum[1000];
    Float_t neutralMax[1000];
    
    jetpbpb->SetBranchAddress("HLT_HIJet80_v1",&HLT_HIJet80_v1);
    jetpbpb->SetBranchAddress("HLT_HIJet65_v1",&HLT_HIJet65_v1);
    jetpbpb->SetBranchAddress("HLT_HIJet55_v1",&HLT_HIJet55_v1);
    jetpbpb->SetBranchAddress("HLT_HIJet55_v1_Prescl",&HLT_HIJet55_v1_Prescl);

    jetpbpb->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
    jetpbpb->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);

    jetpbpb->SetBranchAddress("run",&run);
    jetpbpb->SetBranchAddress("evt",&evt);
    jetpbpb->SetBranchAddress("hiBin",&hiBin);
    jetpbpb->SetBranchAddress("vz",&vz);

    jetpbpb->SetBranchAddress("nref",&nref);
    jetpbpb->SetBranchAddress("rawpt",rawpt);
    jetpbpb->SetBranchAddress("jtpt",jtpt);
    jetpbpb->SetBranchAddress("jteta",jteta);
    jetpbpb->SetBranchAddress("chargedMax",chargedMax);
    jetpbpb->SetBranchAddress("chargedSum",chargedSum);
    jetpbpb->SetBranchAddress("photonMax",photonMax);
    jetpbpb->SetBranchAddress("photonSum",photonSum);
    jetpbpb->SetBranchAddress("neutralMax",neutralMax);
    jetpbpb->SetBranchAddress("neutralSum",neutralSum);

    //std::cout<<" Added the necessary trees " <<endl;


 
    for(int i=0; i<jetpbpb->GetEntries(); i++){
      nbytes += jetpbpb->GetEntry(i);

      //if(i%1000==0)std::cout<<" \t  events processed # "<< i <<std::endl;

      bool selEvent = fabs(vz)<15 && pcollisionEventSelection && pHBHENoiseFilter;
      if(!selEvent)continue;

      int iTrig=-1;

      if(it==0){
	if(HLT_HIJet65_v1 && !HLT_HIJet80_v1)iTrig=1;
	if(HLT_HIJet55_v1_Prescl && HLT_HIJet55_v1 && !HLT_HIJet65_v1 && !HLT_HIJet80_v1)iTrig=2;
      }else if(it==1){
	if(HLT_HIJet80_v1)iTrig=0;
      }

      if(iTrig<0)continue;

      //TCut pbpb1 = "abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_HIJet80_v1&&chargedMax/jtpt>0.01&&hiBin>=0 && hiBin<10";
      //TCut pbpb2 = "abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_HIJet65_v1&&!HLT_HIJet80_v1&&chargedMax/jtpt>0.01&&hiBin>0&&hiBin<10";
      //TCut pbpb3 = "abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_HIJet55_v1&&!HLT_HIJet65_v1&&!HLT_HIJet80_v1&&chargedMax/jtpt>0.01&&hiBin>0&&hiBin<10";
      hBin->Fill(hiBin);

      int icent = GetCentBin(hiBin);
      if(icent <0 || icent>=ncen)continue;

//       //! Jet energy scale comparison with data 
//       vJets.clear();
//       for(int igen=0; igen<nref; igen++){ 
// 	int gj = igen; 
// 	_JEC_HI->setJetEta(jteta[gj]);
// 	_JEC_HI->setJetPt (rawpt[gj]);
// 	float recopt  = rawpt[gj]*_JEC_HI->getCorrection();  //! correction with JECv14
// 	if(recopt<80 && fabs(jteta[gj])>1.4 && chargedMax[ij]/recopt > 0.01)continue;	
// 	vJets.push_back(recopt);
//       }                         
      
//       if(vJets.size()>=2){
// 	std::sort(vJets.begin(),vJets.end());
// 	double B=-9999;
// 	double rn1 = gRandom->Rndm();
// 	double rn2 = gRandom->Rndm();
// 	double ptdij = (vJets[0] + vJets[1])/2.;
// 	if(rn1 > rn2){
// 	  B = (vJets[0] - vJets[1])/(vJets[0] + vJets[1]); 
// 	}else{ 
// 	  B = (vJets[1] - vJets[0])/(vJets[1] + vJets[0]);
// 	}
// 	if(B!=-9999)hJetEnergyScale[nj][icent]->Fill(ptdij,B,wxs*wcen*wvz);
//       }       
      
      
      for(int ij=0; ij<nref;ij++){

	_JEC_HI->setJetEta(jteta[ij]);
        _JEC_HI->setJetPt (rawpt[ij]);
	
	//double recopt = jtpt[ij];
	double recopt = rawpt[ij]*_JEC_HI->getCorrection();

	bool selJet = fabs(jteta[ij])<2 &&  chargedMax[ij]/recopt > 0.01 && recopt>20. && rawpt[ij]>10; 
	if(!selJet)continue;
	
	double jetbkg = (chargedSum[ij] + neutralSum[ij] + photonSum[ij]) - rawpt[ij];
	hpbpb[icent][iTrig]->Fill(recopt);
	hjbkg[icent][iTrig]->Fill(recopt,jetbkg);

	hchargedMax[icent][iTrig]->Fill(recopt,chargedMax[ij]/recopt);
	hphotonMax [icent][iTrig]->Fill(recopt,photonMax [ij]/recopt);
	hneutralMax[icent][iTrig]->Fill(recopt,neutralMax[ij]/recopt);


	hjbkgvz [icent][iTrig]->Fill(vz,recopt,jetbkg);
	hjbkgeta[icent][iTrig]->Fill(jteta[ij],recopt,jetbkg);
	hjbkgraw[icent][iTrig]->Fill(rawpt[ij],recopt,jetbkg);

	if( jetbkg < 40 ){
	  hchargedMaxSp[icent][iTrig]->Fill(recopt,chargedMax[ij]/recopt);
	  hphotonMaxSp [icent][iTrig]->Fill(recopt,photonMax [ij]/recopt);
	  hneutralMaxSp[icent][iTrig]->Fill(recopt,neutralMax[ij]/recopt);

	  //	  if(iTrig==2 && icent==0 && recopt>50 && recopt<100){
	  //std::cout<<"  \t  Run # : " << run << " evt : " <<  evt << "i : " <<ij<<" recopt : "<< recopt << " rawpt : " << rawpt[ij] <<" chsum : "<<chargedSum[ij]<<" phsum : "<<photonSum[ij]<< " nusum : " << neutralSum[ij] << " chMax/recpt : " << chargedMax[ij]/recopt << " phMax/recpt : " << photonMax[ij]/recopt << " nuMax/recpt : " << neutralMax[ij]/recopt << " bkgd : " << jetbkg << std::endl;
	  //	  }
	}
      }
    }
  }

  //respective lumi seen by the trigger all in inverse micro barns
  double lumi[3]={149.382e6,3.195e6,2.734e6};
  for(int ic=0;ic<ncen;ic++){
    for(int in=0;in<3;in++){
      hpbpb[ic][in]->Scale(1./lumi[in]/4.);
      hpbpbComb[ic]->Add(hpbpb[ic][in]);

      hjbkg[ic][in]->Scale(1./lumi[in]);
      hjbkgComb[ic]->Add(hjbkg[ic][in]);

      hjbkgvz[ic][in]->Scale(1./lumi[in]);
      hjbkgvzComb[ic]->Add(hjbkgvz[ic][in]);

      hjbkgeta[ic][in]->Scale(1./lumi[in]);
      hjbkgetaComb[ic]->Add(hjbkgeta[ic][in]);

      hjbkgraw[ic][in]->Scale(1./lumi[in]);
      hjbkgrawComb[ic]->Add(hjbkgraw[ic][in]);


      
      hchargedMax[ic][in]->Scale(1./lumi[in]);
      hchargedMaxComb[ic]->Add(hchargedMax[ic][in]);
      hphotonMax [ic][in]->Scale(1./lumi[in]);
      hphotonMaxComb[ic]->Add(hphotonMax[ic][in]);
      hneutralMax[ic][in]->Scale(1./lumi[in]);
      hneutralMaxComb[ic]->Add(hneutralMax[ic][in]);      

      hchargedMaxSp[ic][in]->Scale(1./lumi[in]);
      hchargedMaxSpComb[ic]->Add(hchargedMaxSp[ic][in]);
      hphotonMaxSp [ic][in]->Scale(1./lumi[in]);
      hphotonMaxSpComb[ic]->Add(hphotonMaxSp[ic][in]);
      hneutralMaxSp[ic][in]->Scale(1./lumi[in]);
      hneutralMaxSpComb[ic]->Add(hneutralMaxSp[ic][in]);      

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
