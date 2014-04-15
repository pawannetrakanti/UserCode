#include "TROOT.h"
#include "TString.h"
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h" 
using namespace std;

int weight_ppTracking(const char *infile="/mnt/hadoop/cms/store/user/belt/Validations53X/Track8_Jet21_cut1/hiForest_QCDpT80_STARTHI53_LV1_Track8_Jet21_9_1GeVcut.root",
		      const char *outfile="dijet_pp_merged_corrected_Track8_Jet21MC.root", 
		      float maxpthat=120., 
		      float  xSection=9.865e-05,
		      bool updateCorrections=false
		      )
  
{

  //Reset ROOT and connect tree file
  gROOT->Reset();

  TFile *fin = TFile::Open(infile, "readonly");
  cout<<infile<<endl;

  cout<<endl;
  cout<<endl;
  fin->ls();
  cout<<endl;
  cout<<endl;

  const int ndir=12;
  
  string inputDirName[ndir]= {
    "akPu3CaloJetAnalyzer",
    "akPu3PFJetAnalyzer",
    "ak3CaloJetAnalyzer",
    "ak3PFJetAnalyzer",
    "akPu4CaloJetAnalyzer",
    "akPu4PFJetAnalyzer",
    "ak4CaloJetAnalyzer",
    "ak4PFJetAnalyzer",
    "akPu5CaloJetAnalyzer",
    "akPu5PFJetAnalyzer",
    "ak5CaloJetAnalyzer",
    "ak5PFJetAnalyzer"
  };

  string outputDirName[ndir]= {
    "akPu3CaloJetAnalyzer",
    "akPu3PFJetAnalyzer",
    "ak3CaloJetAnalyzer",
    "ak3PFJetAnalyzer",
    "akPu4CaloJetAnalyzer",
    "akPu4PFJetAnalyzer",
    "ak4CaloJetAnalyzer",
    "ak4PFJetAnalyzer",
    "akPu5CaloJetAnalyzer",
    "akPu5PFJetAnalyzer",
    "ak5CaloJetAnalyzer",
    "ak5PFJetAnalyzer"
  };
 
  string corrFileName[ndir]= {
    "AKPu3Calo",
    "AKPu3PF",
    "AK3Calo",
    "AK3PF",
    "AKPu4Calo",
    "AKPu4PF",
    "AK4Calo",
    "AK4PF",
    "AKPu5Calo",
    "AKPu5PF",
    "AK5Calo",
    "AK5PF"
  };



  TFile *fout=new TFile(outfile,"RECREATE");

  TTree *tr_in=0, *tr_out=0;
  Int_t           hiBin;  

  TTree *tr_ev = 0;
  if(!updateCorrections){
    tr_ev = (TTree*)fin->Get("hiEvtAnalyzer/HiTree");
    tr_ev->SetBranchAddress("hiBin",&hiBin);
    tr_ev->SetBranchStatus("*",0,0);
    tr_ev->SetBranchStatus("hiBin",1,0);
  }

  for(Int_t idir=0;idir<ndir;idir++){

    cout <<"idir =" << idir <<" JetName ="<< inputDirName[idir].c_str() <<endl ;
    tr_in = (TTree*)fin->Get(Form("%s/t",inputDirName[idir].c_str()));
    
    float fentries = (float)tr_in->GetEntries();
    float weight = xSection/(fentries/1000.);
    cout<<" weight "<<weight<<" \t " << tr_in->GetName() << endl;
    
    //Declaration of leaves types
    int   nref;
    float pthat;
    float corrpt[100];
    float jtpt[100];
    float rawpt[100];
    float jteta[100];
    float jtphi[100];
    float jty[100];
    float jtpu[100];
    float refpt[100];
    float refeta[100];
    float refphi[100];
    float refdphijt[100];       
    float refdrjt[100];       
    float refparton_pt[100];
    int refparton_flavor[100];
    int subid[100];

    tr_in->SetBranchAddress("nref",&nref);
    tr_in->SetBranchAddress("pthat",&pthat);
    tr_in->SetBranchAddress("rawpt",rawpt);
    tr_in->SetBranchAddress("jtpt",jtpt);
    tr_in->SetBranchAddress("jteta",jteta);
    tr_in->SetBranchAddress("jty",jty);
    tr_in->SetBranchAddress("jtphi",jtphi);
    tr_in->SetBranchAddress("jtpu",jtpu);
    tr_in->SetBranchAddress("refpt",refpt);
    tr_in->SetBranchAddress("refphi",refphi);
    tr_in->SetBranchAddress("refeta",refeta);
    tr_in->SetBranchAddress("refdphijt",refdphijt);
    tr_in->SetBranchAddress("refdrjt",refdrjt);
    tr_in->SetBranchAddress("refparton_pt",refparton_pt);
    tr_in->SetBranchAddress("refparton_flavor",refparton_flavor);
    tr_in->SetBranchAddress("subid",subid);

    if(updateCorrections){
      tr_in->SetBranchAddress("hiBin",&hiBin);
      tr_in->SetBranchAddress("weight",&weight);

      maxpthat=9999;
    }
    cout<<"get jet trees!!! "<<endl;

    //! Add Friends to the TTree
    if(!updateCorrections)tr_in->AddFriend(tr_ev);

    fout->mkdir(outputDirName[idir].c_str());
    fout->cd(outputDirName[idir].c_str());

    tr_out = new TTree("t","Jet  Response Analyzer");
    tr_out->SetMaxTreeSize(4200000000);


    // Set output branch addresses.
    tr_out->Branch("hiBin",&hiBin,"hiBin/I");
    tr_out->Branch("nref",&nref,"nref/I");
    tr_out->Branch("pthat",&pthat,"pthat/F");
    tr_out->Branch("weight",&weight,"weight/F");
    tr_out->Branch("jtpt",jtpt,"jtpt[nref]/F");
    tr_out->Branch("rawpt",rawpt,"rawpt[nref]/F");
    tr_out->Branch("corrpt",corrpt,"corrpt[nref]/F");
    tr_out->Branch("jtpu",jtpu,"jtpu[nref]/F");
    tr_out->Branch("jteta",jteta,"jteta[nref]/F");
    tr_out->Branch("jty",jty,"jty[nref]/F");
    tr_out->Branch("jtphi",jtphi,"jtphi[nref]/F");
    tr_out->Branch("refpt",refpt,"refpt[nref]/F");
    tr_out->Branch("refeta",refeta,"refeta[nref]/F");
    tr_out->Branch("refphi",refphi,"refphi[nref]/F");
    tr_out->Branch("refdphijt",refdphijt,"refdphijt[nref]/F");
    tr_out->Branch("refdrjt",refdrjt,"refdrjt[nref]/F");
    tr_out->Branch("refparton_pt",refparton_pt,"refparton_pt[nref]/F");
    tr_out->Branch("refparton_flavor",refparton_flavor,"refparton_flavor[nref]/I");
    tr_out->Branch("subid",subid,"subid[nref]/I");


    Long64_t nbytes = 0;
    //! grab the JEC's  
    string L2Name, L3Name;
    JetCorrectorParameters *parHI_l2=0, *parHI_l3=0;
    vector<JetCorrectorParameters> vpar_HI;
    FactorizedJetCorrector *_JEC_HI = new FactorizedJetCorrector(vpar_HI);//JR 

    if(updateCorrections){
      L2Name = "/net/hisrv0001/home/pawan/Validation/Track8_Jet19/CMSSW_5_3_16/src/JEC_base/pp2014/JEC_pp2014ReReco53x_dijet_L2Relative_"+corrFileName[idir]+".txt";
      cout<<"**** ++++++  L2Name "<<L2Name<<endl;
      L3Name = "/net/hisrv0001/home/pawan/Validation/Track8_Jet19/CMSSW_5_3_16/src/JEC_base/pp2014/JEC_pp2014ReReco53x_dijet_L3Absolute_"+corrFileName[idir]+".txt";
      cout<<"**** ++++++  L3Name "<<L3Name<<endl;
      
      parHI_l2 = new JetCorrectorParameters(L2Name.c_str());
      parHI_l3 = new JetCorrectorParameters(L3Name.c_str());
      
      vpar_HI.push_back(*parHI_l2);
      vpar_HI.push_back(*parHI_l3);
      _JEC_HI = new FactorizedJetCorrector(vpar_HI);
    }

    
    for (Long64_t i=0; i<fentries;i++) {
      nbytes += tr_in->GetEntry(i);


      if(pthat>maxpthat)continue;	 

      //! jet loop
      for(int iref=0;iref<nref;iref++){

 	corrpt[iref] = jtpt[iref];
 	jtpt[iref]   = rawpt[iref];

	if(updateCorrections){
          _JEC_HI->setJetEta(jteta[iref]);
          _JEC_HI->setJetPt(rawpt[iref]);
          corrpt[iref] = rawpt[iref]*_JEC_HI->getCorrection();
	  //if(TMath::Abs(corrpt[iref]-(rawpt[iref]*_JEC_HI->getCorrection()))>1. && rawpt[iref]>10.) cout<<"*****!!!!!"<<endl;
	  //if(TMath::Abs(_JEC_HI->getCorrection()-(corrpt[iref]/rawpt[iref]))>0.2 && rawpt[iref]>10.) cout<<"iref: "<<iref<<" correction is: "<<(_JEC_HI->getCorrection())<<" ratio of corr/raw: "<<(corrpt    [iref]/rawpt[iref])<<" for rawpt: "<<rawpt[iref]<<" corrpt: "<<corrpt[iref]<<" and rawpT*correction: "<<(rawpt[iref]*_JEC_HI->getCorrection())<<" and correction*rawpt: "<<(_JEC_HI->getCorrection()*rawpt[iref])<<" and eta: "<<jteta[iref]<<endl;
	}
      }//ref for loop
      tr_out->Fill();
    }//nentries loop

    tr_out->Write();

    delete parHI_l2;
    delete parHI_l3;
    delete _JEC_HI;
    
  }    

  cout <<"finish the code!!"<<endl ;
  fout->Close();
  cout<<"XXXX"<<endl;

  return 0;
}
