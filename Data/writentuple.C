// From the prompt reco
#include "/afs/cern.ch/user/p/pawan/scratch0/CMSSW_6_2_0/src/work/pPb/HiForest/V3/hiForest.h"
#include <TChain.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TF1.h>
#include <TFile.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom.h>
#include <cstdlib>
#include <cmath>
#include <iostream>
using namespace std;

const float ketacut=2.0;
const float kptrecocut=30.;

void LoadLib();
void ShutoffBranches(HiForest */*hi*/);
void FindLeadSubLeadJets(Jets */*mJets*/, int */*ljet*/);

TStopwatch timer;
int writentuple(char *ksp="ppJet40")
{

  timer.Start();

  LoadLib();

  TString inname="";
  if(strcmp(ksp,"ppJet40")==0)inname = "root://eoscms//eos/cms/store/group/phys_heavyions/yjlee/pp2013/promptReco/PP2013_HiForest_PromptReco_JSon_Jet40Jet60_ppTrack_forestv84.root";
  else if(strcmp(ksp,"ppJet80")==0)inname = "root://eoscms//eos/cms/store/caf/user/yjlee/pp2013/promptReco/PP2013_HiForest_PromptReco_JsonPP_Jet80_PPReco_forestv82.root";

  //! Load Lib
  //gSystem->Load("/afs/cern.ch/user/p/pawan/scratch0/CMSSW_6_2_0/src/work/pPb/HiForest/V3/hiForest_h.so");

  //! Define the input file and HiForest
  //! CMSSW_5_3_3
  HiForest *c = new HiForest(inname,Form("Forest%s",ksp),cPP);
  cout<<"Loaded the hiforest tree : "<<c->GetName()<<endl;
  ShutoffBranches(c);

  //! Output file
  //! HIHighPt
  TFile *fout = new TFile(Form("ntuple_2013_%s.root",ksp),"RECREATE");  

  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"**************************************************** "<<std::endl;
  std::cout<<Form("Running for %s ",ksp)<<std::endl;
  std::cout<<Form("pT  cut for %0.3f ",kptrecocut)<<std::endl;
  std::cout<<Form("eta cut for %0.3f ",ketacut)<<std::endl;
  std::cout<<"My hiForest Tree : " <<c->GetName()<<"\t Entries "<<c->GetEntries()<<std::endl;
  std::cout<<"Output file  "<<fout->GetName()<<std::endl;
  std::cout<<"**************************************************** "<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;

  //! shut off jet trees
  c->hasAk2CaloJetTree=0;
  c->hasAk4CaloJetTree=0;
  c->hasAk3CaloJetTree=0;
  c->hasAk5CaloJetTree=0;

  c->hasAkPu2CaloJetTree=0;
  c->hasAkPu4CaloJetTree=0;
  c->hasAkPu3CaloJetTree=0;
  c->hasAkPu5CaloJetTree=0;

  c->hasAk2PFJetTree=0;
  c->hasAk4PFJetTree=0;
  c->hasAk5PFJetTree=0;

  c->hasAkPu2PFJetTree=0;
  c->hasAkPu4PFJetTree=0;
  c->hasAkPu5PFJetTree=0;

  c->hasTrackTree=0;

  //! For jets
  Jets *mJets=0;
  Long64_t nentries = c->GetEntries();
  std::cout<<Form("# of entries in TTree for %s : ",ksp)<<nentries<<std::endl;

  string jetVars = "";
  jetVars += "evt:run:vz:trig:jet40:jet60:jet80:jet100:ntrk:pt1:raw1:eta1:phi1:chMax1:chSum1:phSum1:neSum1:pt2:raw2:eta2:phi2:chMax2:chSum2:phSum2:neSum2:pt3:raw3:eta3:phi3:chMax3:chSum3:phSum3:neSum3";
  TNtuple *ntjet=0;
  ntjet = new TNtuple("ntjet","",jetVars.data());

  for (Long64_t ievt=0; ievt<nentries;ievt++) {//! event loop
  //for (Long64_t ievt=0; ievt<100;ievt++) {//! event loop
    //! load the hiForest event
    c->GetEntry(ievt);
    
    //! events with Single vertex
    bool evSel = false;
    float trig=-9;
    if(strcmp(ksp,"ppJet40")==0){
      evSel = fabs(c->evt.vz)<15. && c->skim.pHBHENoiseFilter  && c->skim.pPAcollisionEventSelectionPA && (c->hlt.HLT_PAJet40_NoJetID_v1 || c->hlt.HLT_PAJet60_NoJetID_v1); 
      trig=1;
    }
    else if(strcmp(ksp,"ppJet80")==0){
      evSel = fabs(c->evt.vz)<15. && c->skim.pHBHENoiseFilter  && c->skim.pPAcollisionEventSelectionPA && (c->hlt.HLT_PAJet80_NoJetID_v1 || c->hlt.HLT_PAJet100_NoJetID_v1); 
      trig=2;
    }
    if(!evSel)continue;

    float pt1 = -9, pt2    = -9, pt3    = -9,
      raw1    = -9, raw2   = -9, raw3   = -9,
      eta1    = -9, eta2   = -9, eta3   = -9,
      phi1    = -9, phi2   = -9, phi3   = -9,
      chMax1  = -9, chMax2 = -9, chMax3 = -9,
      chSum1  = -9, chSum2 = -9, chSum3 = -9,
      phSum1  = -9, phSum2 = -9, phSum3 = -9,
      neSum1  = -9, neSum2 = -9, neSum3 = -9;

    
    float run   = c->evt.run;
    float evt   = c->evt.evt;
    float vz    = c->evt.vz;
    
    float jet40  = c->hlt.HLT_PAJet40_NoJetID_v1;
    float jet60  = c->hlt.HLT_PAJet60_NoJetID_v1;
    float jet80  = c->hlt.HLT_PAJet80_NoJetID_v1;
    float jet100 = c->hlt.HLT_PAJet100_NoJetID_v1;
    
    float ntrk    = c->evt.hiNtracks;

    if(ievt%10000==0)std::cout<<" ********** Event # " <<ievt<<"\t Run : "<<run<<std::endl;
    
    mJets = &(c->ak3PF);
    
    int *ljet = new int[3];
    FindLeadSubLeadJets(mJets,ljet);
    
    int jtLead = -1, jtSubLead = -1, jtThird = -1;
    
    if(ljet[0] >=0 ) jtLead    = ljet[0];
    if(ljet[1] >=0 ) jtSubLead = ljet[1];    
    if(ljet[2] >=0 ) jtThird   = ljet[2];    
    
    if(jtLead<0)continue;
    
    if(jtLead > -1){
      pt1     = mJets->jtpt[jtLead];
      eta1    = mJets->jteta[jtLead];
      phi1    = mJets->jtphi[jtLead];
      raw1    = mJets->rawpt[jtLead];
      chMax1  = mJets->chargedMax[jtLead];
      chSum1  = mJets->chargedSum[jtLead];
      phSum1  = mJets->photonSum[jtLead];
      neSum1  = mJets->neutralSum[jtLead];
    }
    
    if(jtSubLead > -1){
      pt2     = mJets->jtpt[jtSubLead];
      eta2    = mJets->jteta[jtSubLead];
      phi2    = mJets->jtphi[jtSubLead];
      raw2    = mJets->rawpt[jtSubLead];
      chMax2  = mJets->chargedMax[jtSubLead];
      chSum2  = mJets->chargedSum[jtSubLead];
      phSum2  = mJets->photonSum [jtSubLead];
      neSum2  = mJets->neutralSum[jtSubLead];
    }
    
    if(jtThird > -1){
      pt3     = mJets->jtpt[jtThird];
      eta3    = mJets->jteta[jtThird];
      phi3    = mJets->jtphi[jtThird];
      raw3    = mJets->rawpt[jtThird];
      chMax3  = mJets->chargedMax[jtThird];
      chSum3  = mJets->chargedSum[jtThird];
      phSum3  = mJets->photonSum [jtThird];
      neSum3  = mJets->neutralSum[jtThird];
    }
    
    float jentry[] = {evt,run,vz,trig,jet40,jet60,jet80,jet100,ntrk,
		       pt1,raw1,eta1,phi1,chMax1,chSum1,phSum1,neSum1,
		       pt2,raw2,eta2,phi2,chMax2,chSum2,phSum2,neSum2,
		       pt3,raw3,eta3,phi3,chMax3,chSum3,phSum3,neSum3
    };
    
    ntjet->Fill(jentry);
    delete [] ljet;
  }//! event loop ends


  //! Write to output file
  fout->cd();
  fout->Write();
  fout->Close();

  //! Check
  timer.Stop();
  float rtime  = timer.RealTime();
  float ctime  = timer.CpuTime();

  std::cout<<"\t"<<std::endl;
  std::cout<<Form("RealTime=%f seconds, CpuTime=%f seconds",rtime,ctime)<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"Good bye : " <<"\t"<<std::endl;

  return 0;
}
void FindLeadSubLeadJets(Jets *jetc, int *ljet)
{
  ljet[0]=-1; ljet[1]=-2; ljet[2]=-3;
  
  float tempt=-9;
  //! Get the leading jet
  for(int ij=0; ij<jetc->nref; ij++){
    if(fabs(jetc->jteta[ij])>ketacut || jetc->jtpt[ij]<kptrecocut || jetc->rawpt[ij]<15)continue;
    float jetpt = jetc->jtpt[ij];
    if(jetpt > tempt){
      tempt = jetpt;
      ljet[0] = ij;
    }
  }
  if(ljet[0]>=0){
    // Subleading
    tempt=-9;
    for(int ij=0; ij<jetc->nref; ij++){
      if(ij==ljet[0])continue;
      if(fabs(jetc->jteta[ij])>ketacut || jetc->jtpt[ij]>jetc->jtpt[ljet[0]] || jetc->rawpt[ij]<15)continue;
      float jetpt = jetc->jtpt[ij];
      if (jetpt > tempt){
	tempt = jetpt;
	ljet[1] = ij;
      }
    }

    if(ljet[1]>=0){
      // third jet
      tempt=-9;
      for(int ij=0; ij<jetc->nref; ij++){
	if(ij==ljet[0] || ij==ljet[1])continue;
	if(fabs(jetc->jteta[ij])>ketacut || jetc->jtpt[ij]>jetc->jtpt[ljet[0]] || jetc->rawpt[ij]<15)continue;
	float jetpt = jetc->jtpt[ij];
	if (jetpt > tempt){
	  tempt = jetpt;
	  ljet[2] = ij;
	}
      }
    }//! third jet
  }
}
void LoadLib()
{
  gSystem->Load("/afs/cern.ch/user/p/pawan/scratch0/CMSSW_6_2_0/src/work/pPb/HiForest/V3/hiForest_h.so");
}
void ShutoffBranches(HiForest *hi)
{

  //! added by pawan
  //! Select only the branches you want to use for the analysis
  //! This increases the speed for running

  //! For Hlt
  hi->hltTree->SetBranchStatus("*",0,0);
  hi->hltTree->SetBranchStatus("HLT_PAJet40_NoJetID_v1",1,0);
  hi->hltTree->SetBranchStatus("HLT_PAJet60_NoJetID_v1",1,0);
  hi->hltTree->SetBranchStatus("HLT_PAJet80_NoJetID_v1",1,0);
  hi->hltTree->SetBranchStatus("HLT_PAJet100_NoJetID_v1",1,0);

  //! for Skim Tree
  hi->skimTree->SetBranchStatus("*",0,0);
  hi->skimTree->SetBranchStatus("pcollisionEventSelection",1,0);
  hi->skimTree->SetBranchStatus("pPAcollisionEventSelectionPA",1,0);
  hi->skimTree->SetBranchStatus("pHBHENoiseFilter",1,0);

  //! Evt tree
  hi->evtTree->SetBranchStatus("*",0,0);
  hi->evtTree->SetBranchStatus("run",1,0);
  hi->evtTree->SetBranchStatus("evt",1,0);
  hi->evtTree->SetBranchStatus("vz",1,0);
  hi->evtTree->SetBranchStatus("hiNtracks",1,0);

  hi->ak3PFJetTree->SetBranchStatus("*",0,0);
  hi->ak3PFJetTree->SetBranchStatus("nref",1,0);
  hi->ak3PFJetTree->SetBranchStatus("rawpt",1,0);
  hi->ak3PFJetTree->SetBranchStatus("jtpt",1,0);
  hi->ak3PFJetTree->SetBranchStatus("jteta",1,0);
  hi->ak3PFJetTree->SetBranchStatus("jtphi",1,0);
  hi->ak3PFJetTree->SetBranchStatus("chargedMax",1,0);
  hi->ak3PFJetTree->SetBranchStatus("chargedSum",1,0);
  hi->ak3PFJetTree->SetBranchStatus("photonSum",1,0);
  hi->ak3PFJetTree->SetBranchStatus("neutralSum",1,0);



  //! Track tree
  //hi->trackTree->SetBranchStatus("*",0,0);
  //hi->trackTree->SetBranchStatus("nTrk",1,0);
  // hi->trackTree->SetBranchStatus("trkPt",1,0);
  // hi->trackTree->SetBranchStatus("trkEta",1,0);
  // hi->trackTree->SetBranchStatus("trkPhi",1,0);
  // hi->trackTree->SetBranchStatus("highPurity",1,0);
  // hi->trackTree->SetBranchStatus("trkDz1",1,0);
  // hi->trackTree->SetBranchStatus("trkDzError1",1,0);
  // hi->trackTree->SetBranchStatus("trkDxy1",1,0);
  // hi->trackTree->SetBranchStatus("trkDxyError1",1,0);

  /*
  //! PF jet tree
  hi->ak2PFJetTree->SetBranchStatus("*",0,0);
  hi->ak2PFJetTree->SetBranchStatus("nref",1,0);
  hi->ak2PFJetTree->SetBranchStatus("rawpt",1,0);
  hi->ak2PFJetTree->SetBranchStatus("jtpt",1,0);
  hi->ak2PFJetTree->SetBranchStatus("jteta",1,0);
  hi->ak2PFJetTree->SetBranchStatus("jtphi",1,0);
  hi->ak2PFJetTree->SetBranchStatus("trackMax",1,0);
  */
  /*
  hi->ak3PFJetTree->SetBranchStatus("*",0,0);
  hi->ak3PFJetTree->SetBranchStatus("nref",1,0);
  hi->ak3PFJetTree->SetBranchStatus("rawpt",1,0);
  hi->ak3PFJetTree->SetBranchStatus("jtpt",1,0);
  hi->ak3PFJetTree->SetBranchStatus("jteta",1,0);
  hi->ak3PFJetTree->SetBranchStatus("jtphi",1,0);
  hi->ak3PFJetTree->SetBranchStatus("trackMax",1,0);
  */
  /*
  hi->ak4PFJetTree->SetBranchStatus("*",0,0);
  hi->ak4PFJetTree->SetBranchStatus("nref",1,0);
  hi->ak4PFJetTree->SetBranchStatus("rawpt",1,0);
  hi->ak4PFJetTree->SetBranchStatus("jtpt",1,0);
  hi->ak4PFJetTree->SetBranchStatus("jteta",1,0);
  hi->ak4PFJetTree->SetBranchStatus("jtphi",1,0);
  hi->ak4PFJetTree->SetBranchStatus("trackMax",1,0);
  */
  /*
  hi->ak5PFJetTree->SetBranchStatus("*",0,0);
  hi->ak5PFJetTree->SetBranchStatus("nref",1,0);
  hi->ak5PFJetTree->SetBranchStatus("rawpt",1,0);
  hi->ak5PFJetTree->SetBranchStatus("jtpt",1,0);
  hi->ak5PFJetTree->SetBranchStatus("jteta",1,0);
  hi->ak5PFJetTree->SetBranchStatus("jtphi",1,0);
  hi->ak5PFJetTree->SetBranchStatus("trackMax",1,0);
  */

  //
  /*
  hi->akPu2PFJetTree->SetBranchStatus("*",0,0);
  hi->akPu2PFJetTree->SetBranchStatus("nref",1,0);
  hi->akPu2PFJetTree->SetBranchStatus("rawpt",1,0);
  hi->akPu2PFJetTree->SetBranchStatus("jtpt",1,0);
  hi->akPu2PFJetTree->SetBranchStatus("jteta",1,0);
  hi->akPu2PFJetTree->SetBranchStatus("jtphi",1,0);
  hi->akPu2PFJetTree->SetBranchStatus("trackMax",1,0);
  */
  /*
  hi->akPu3PFJetTree->SetBranchStatus("*",0,0);
  hi->akPu3PFJetTree->SetBranchStatus("nref",1,0);
  hi->akPu3PFJetTree->SetBranchStatus("rawpt",1,0);
  hi->akPu3PFJetTree->SetBranchStatus("jtpt",1,0);
  hi->akPu3PFJetTree->SetBranchStatus("jteta",1,0);
  hi->akPu3PFJetTree->SetBranchStatus("jtphi",1,0);
  hi->akPu3PFJetTree->SetBranchStatus("trackMax",1,0);
  */
  /*
  hi->akPu4PFJetTree->SetBranchStatus("*",0,0);
  hi->akPu4PFJetTree->SetBranchStatus("nref",1,0);
  hi->akPu4PFJetTree->SetBranchStatus("rawpt",1,0);
  hi->akPu4PFJetTree->SetBranchStatus("jtpt",1,0);
  hi->akPu4PFJetTree->SetBranchStatus("jteta",1,0);
  hi->akPu4PFJetTree->SetBranchStatus("jtphi",1,0);
  hi->akPu4PFJetTree->SetBranchStatus("trackMax",1,0);
  */
  /*
  hi->akPu5PFJetTree->SetBranchStatus("*",0,0);
  hi->akPu5PFJetTree->SetBranchStatus("nref",1,0);
  hi->akPu5PFJetTree->SetBranchStatus("rawpt",1,0);
  hi->akPu5PFJetTree->SetBranchStatus("jtpt",1,0);
  hi->akPu5PFJetTree->SetBranchStatus("jteta",1,0);
  hi->akPu5PFJetTree->SetBranchStatus("jtphi",1,0);
  hi->akPu5PFJetTree->SetBranchStatus("trackMax",1,0);
  */
  //


  //! Calo jet trees
  /*
  hi->ak2CaloJetTree->SetBranchStatus("*",0,0);
  hi->ak2CaloJetTree->SetBranchStatus("nref",1,0);
  hi->ak2CaloJetTree->SetBranchStatus("rawpt",1,0);
  hi->ak2CaloJetTree->SetBranchStatus("jtpt",1,0);
  hi->ak2CaloJetTree->SetBranchStatus("jteta",1,0);
  hi->ak2CaloJetTree->SetBranchStatus("jtphi",1,0);
  hi->ak2CaloJetTree->SetBranchStatus("trackMax",1,0);

  hi->ak3CaloJetTree->SetBranchStatus("*",0,0);
  hi->ak3CaloJetTree->SetBranchStatus("nref",1,0);
  hi->ak3CaloJetTree->SetBranchStatus("rawpt",1,0);
  hi->ak3CaloJetTree->SetBranchStatus("jtpt",1,0);
  hi->ak3CaloJetTree->SetBranchStatus("jteta",1,0);
  hi->ak3CaloJetTree->SetBranchStatus("jtphi",1,0);
  hi->ak3CaloJetTree->SetBranchStatus("trackMax",1,0);

  hi->ak4CaloJetTree->SetBranchStatus("*",0,0);
  hi->ak4CaloJetTree->SetBranchStatus("nref",1,0);
  hi->ak4CaloJetTree->SetBranchStatus("rawpt",1,0);
  hi->ak4CaloJetTree->SetBranchStatus("jtpt",1,0);
  hi->ak4CaloJetTree->SetBranchStatus("jteta",1,0);
  hi->ak4CaloJetTree->SetBranchStatus("jtphi",1,0);
  hi->ak4CaloJetTree->SetBranchStatus("trackMax",1,0);

  hi->ak5CaloJetTree->SetBranchStatus("*",0,0);
  hi->ak5CaloJetTree->SetBranchStatus("nref",1,0);
  hi->ak5CaloJetTree->SetBranchStatus("rawpt",1,0);
  hi->ak5CaloJetTree->SetBranchStatus("jtpt",1,0);
  hi->ak5CaloJetTree->SetBranchStatus("jteta",1,0);
  hi->ak5CaloJetTree->SetBranchStatus("jtphi",1,0);
  hi->ak5CaloJetTree->SetBranchStatus("trackMax",1,0);
  */

  ////
  /*
  hi->akPu2CaloJetTree->SetBranchStatus("*",0,0);
  hi->akPu2CaloJetTree->SetBranchStatus("nref",1,0);
  hi->akPu2CaloJetTree->SetBranchStatus("rawpt",1,0);
  hi->akPu2CaloJetTree->SetBranchStatus("jtpt",1,0);
  hi->akPu2CaloJetTree->SetBranchStatus("jteta",1,0);
  hi->akPu2CaloJetTree->SetBranchStatus("jtphi",1,0);
  hi->akPu2CaloJetTree->SetBranchStatus("trackMax",1,0);

  hi->akPu3CaloJetTree->SetBranchStatus("*",0,0);
  hi->akPu3CaloJetTree->SetBranchStatus("nref",1,0);
  hi->akPu3CaloJetTree->SetBranchStatus("rawpt",1,0);
  hi->akPu3CaloJetTree->SetBranchStatus("jtpt",1,0);
  hi->akPu3CaloJetTree->SetBranchStatus("jteta",1,0);
  hi->akPu3CaloJetTree->SetBranchStatus("jtphi",1,0);
  hi->akPu3CaloJetTree->SetBranchStatus("trackMax",1,0);

  hi->akPu4CaloJetTree->SetBranchStatus("*",0,0);
  hi->akPu4CaloJetTree->SetBranchStatus("nref",1,0);
  hi->akPu4CaloJetTree->SetBranchStatus("rawpt",1,0);
  hi->akPu4CaloJetTree->SetBranchStatus("jtpt",1,0);
  hi->akPu4CaloJetTree->SetBranchStatus("jteta",1,0);
  hi->akPu4CaloJetTree->SetBranchStatus("jtphi",1,0);
  hi->akPu4CaloJetTree->SetBranchStatus("trackMax",1,0);

  hi->akPu5CaloJetTree->SetBranchStatus("*",0,0);
  hi->akPu5CaloJetTree->SetBranchStatus("nref",1,0);
  hi->akPu5CaloJetTree->SetBranchStatus("rawpt",1,0);
  hi->akPu5CaloJetTree->SetBranchStatus("jtpt",1,0);
  hi->akPu5CaloJetTree->SetBranchStatus("jteta",1,0);
  hi->akPu5CaloJetTree->SetBranchStatus("jtphi",1,0);
  hi->akPu5CaloJetTree->SetBranchStatus("trackMax",1,0);
  */
  ///
  std::cout<<"Loaded all tree variables "<<std::endl;

}
