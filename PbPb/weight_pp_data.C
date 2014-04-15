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

int weight_pp_data(const char *infile="/mnt/hadoop/cms/store/user/velicanu/hiForest_Jet80or95_GR_R_53_LV6_02Mar2014_1300CET_Track8_Jet15/0.root",
		   const char *outfile="/net/hidsk0001/d00/scratch/pawan/ppData/dijet_pp_data_Jet80or95_corrected.root")
  
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
    "akVs3CaloJetAnalyzer",
    "akVs3PFJetAnalyzer",
    "akVs4CaloJetAnalyzer",
    "akVs4PFJetAnalyzer",
    "akVs5CaloJetAnalyzer",
    "akVs5PFJetAnalyzer",

//      "ak3CaloJetAnalyzer",
//      "ak3PFJetAnalyzer",
//      "ak4CaloJetAnalyzer",
//      "ak4PFJetAnalyzer",
//      "ak5CaloJetAnalyzer",
//      "ak5PFJetAnalyzer",

    "akPu3CaloJetAnalyzer",
    "akPu3PFJetAnalyzer",
    "akPu4CaloJetAnalyzer",
    "akPu4PFJetAnalyzer",
    "akPu5CaloJetAnalyzer",
    "akPu5PFJetAnalyzer"
  };

  string outputDirName[ndir]= {
    "akVs3CaloJetAnalyzer",
    "akVs3PFJetAnalyzer",
    "akVs4CaloJetAnalyzer",
    "akVs4PFJetAnalyzer",
    "akVs5CaloJetAnalyzer",
    "akVs5PFJetAnalyzer",

//     "ak3CaloJetAnalyzer",
//     "ak3PFJetAnalyzer",
//     "ak4CaloJetAnalyzer",
//     "ak4PFJetAnalyzer",
//     "ak5CaloJetAnalyzer",
//     "ak5PFJetAnalyzer",
    
    "akPu3CaloJetAnalyzer",
    "akPu3PFJetAnalyzer",
    "akPu4CaloJetAnalyzer",
    "akPu4PFJetAnalyzer",
    "akPu5CaloJetAnalyzer",
    "akPu5PFJetAnalyzer"
  };
 
  string corrFileName[ndir]= {
    "AKVs3Calo",
    "AKVs3PF",
    "AKVs4Calo",
    "AKVs4PF",
    "AKVs5Calo",
    "AKVs5PF",

//      "AK3Calo",
//      "AK3PF",
//      "AK4Calo",
//      "AK4PF",
//      "AK5Calo",
//      "AK5PF",

     "AK3Calo",
     "AK3PF",
     "AK4Calo",
     "AK4PF",
     "AK5Calo",
     "AK5PF"

//      "AKPu3Calo",
//      "AKPu3PF",
//      "AKPu4Calo",
//      "AKPu4PF",
//      "AKPu5Calo",
//      "AKPu5PF"
  };



  TFile *fout=new TFile(outfile,"RECREATE");

  TTree *tr_in=0, *tr_out=0;
  Int_t           hiBin;  
  double  vz=0;
  TTree *tr_ev = 0;
  tr_ev = (TTree*)fin->Get("hiEvtAnalyzer/HiTree");
  tr_ev->SetBranchAddress("hiBin",&hiBin);
  tr_ev->SetBranchStatus("*",0,0);
  tr_ev->SetBranchStatus("hiBin",1,0);
  tr_ev->SetBranchStatus("vz",1,0);

  for(Int_t idir=0;idir<ndir;idir++){

    cout <<"idir =" << idir <<" JetName ="<< inputDirName[idir].c_str() <<endl ;
    tr_in = (TTree*)fin->Get(Form("%s/t",inputDirName[idir].c_str()));
    
    float fentries = (float)tr_in->GetEntries();
    //cout<<" weight "<<weight<<" \t " << tr_in->GetName() << " entries : " << fentries << endl;
    
    //Declaration of leaves types
    int   nref;
    float corrpt[1000];
    float jtpt[1000];
    float rawpt[1000];
    float jteta[1000];
    float jtphi[1000];
    float jty[1000];
    float jtpu[1000];
    float chargedMax[1000];
    float chargedSum[1000];
    float photonMax[1000];
    float photonSum[1000];
    float neutralMax[1000];
    float neutralSum[1000];

    tr_in->SetBranchAddress("nref",&nref);
    tr_in->SetBranchAddress("rawpt",rawpt);
    tr_in->SetBranchAddress("jtpt",jtpt);
    tr_in->SetBranchAddress("jteta",jteta);
    tr_in->SetBranchAddress("jty",jty);
    tr_in->SetBranchAddress("jtphi",jtphi);
    tr_in->SetBranchAddress("jtpu",jtpu);
    tr_in->SetBranchAddress("chargedMax",chargedMax);
    tr_in->SetBranchAddress("chargedSum",chargedSum);
    tr_in->SetBranchAddress("photonMax",photonMax);
    tr_in->SetBranchAddress("photonSum",photonSum);
    tr_in->SetBranchAddress("neutralMax",neutralMax);
    tr_in->SetBranchAddress("neutralSum",neutralSum);

    cout<<"get jet trees!!! "<<endl;

    //! Add Friends to the TTree
    tr_in->AddFriend(tr_ev);


    fout->mkdir(outputDirName[idir].c_str());
    fout->cd(outputDirName[idir].c_str());

    tr_out = new TTree("t","Jet  Response Analyzer");
    tr_out->SetMaxTreeSize(4200000000);


    // Set output branch addresses.
    tr_out->Branch("hiBin",&hiBin,"hiBin/I");
    tr_out->Branch("nref",&nref,"nref/I");
    tr_out->Branch("jtpt",jtpt,"jtpt[nref]/F");
    tr_out->Branch("rawpt",rawpt,"rawpt[nref]/F");
    tr_out->Branch("corrpt",corrpt,"corrpt[nref]/F");
    tr_out->Branch("jtpu",jtpu,"jtpu[nref]/F");
    tr_out->Branch("jteta",jteta,"jteta[nref]/F");
    tr_out->Branch("jty",jty,"jty[nref]/F");
    tr_out->Branch("jtphi",jtphi,"jtphi[nref]/F");
    tr_out->Branch("chargedMax",chargedMax,"chargedMax[nref]/F");
    tr_out->Branch("chargedSum",chargedSum,"chargedSum[nref]/F");
    tr_out->Branch("photonMax" ,photonMax,"photonMax[nref]/F");
    tr_out->Branch("photonSum" ,photonSum,"photonSum[nref]/F");
    tr_out->Branch("neutralMax",neutralMax,"neutralMax[nref]/F");
    tr_out->Branch("neutralSum",neutralSum,"neutralSum[nref]/F");

    Long64_t nbytes = 0;
    //! grab the JEC's  
    string L2Name, L3Name;
    JetCorrectorParameters *parHI_l2=0, *parHI_l3=0;
    vector<JetCorrectorParameters> vpar_HI;
    FactorizedJetCorrector *_JEC_HI = new FactorizedJetCorrector(vpar_HI);//JR 

    L2Name = "/net/hisrv0001/home/pawan/Validation/Track8_Jet19/CMSSW_5_3_16/src/JEC_base/pp2014/JEC_STARTHI53_LV1_Track8_Jet22_dijet_L2Relative_"+corrFileName[idir]+".txt";
    cout<<"**** ++++++  L2Name "<<L2Name<<endl;
    L3Name = "/net/hisrv0001/home/pawan/Validation/Track8_Jet19/CMSSW_5_3_16/src/JEC_base/pp2014/JEC_STARTHI53_LV1_Track8_Jet22_dijet_L3Absolute_"+corrFileName[idir]+".txt";
    cout<<"**** ++++++  L3Name "<<L3Name<<endl;
    
    parHI_l2 = new JetCorrectorParameters(L2Name.c_str());
    parHI_l3 = new JetCorrectorParameters(L3Name.c_str());
    
    vpar_HI.push_back(*parHI_l2);
    vpar_HI.push_back(*parHI_l3);
    _JEC_HI = new FactorizedJetCorrector(vpar_HI);
    //continue;
    
    for (Long64_t i=0; i<fentries;i++) {
      nbytes += tr_in->GetEntry(i);
      if(fabs(vz)>15.)continue;

      //! jet loop
      for(int iref=0;iref<nref;iref++){
	
	_JEC_HI->setJetEta(jteta[iref]);
	_JEC_HI->setJetPt(rawpt[iref]);
	corrpt[iref] = rawpt[iref]*_JEC_HI->getCorrection();
	
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
