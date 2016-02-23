#include <iostream>
#include <sstream>
#include <string>

#include "TROOT.h"
#include "TString.h"
#include <TEventList.h>
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1.h"
using namespace std;

// //! pp 2.76 TeV PYTHIA cross sections
// const double xsec[12][3] ={{2.034e-01,15,30}, //! 15
// 			   {1.075e-02,30,50}, //! 30
// 			   {1.025e-03,50,80}, //! 50
// 			   {9.865e-05,80,120}, //! 80
// 			   {1.129e-05,120,170}, //! 120
// 			   {1.465e-06,170,220}, //! 170
// 			   {2.837e-07,220,280}, //! 220
// 			   {5.323e-08,280,370}, //! 280
// 			   {5.934e-09,370,460}, //! 370
// 			   {8.125e-10,460,540}, //! 460
// 			   {1.467e-10,540,9999}, //! 540
// 			   {0,9999,9999}
// };


const int npthat=11;
//const int npthat=10; //! HI reco new one with 460 missing 
//const int npthat=10; //! pp reco new one with 370 757p1 missing 
//! pp 5.02 TeV Weight factors
//!                         pthat effEv  effxec       wt 

float pthatwt_pp[npthat][4] ={
  //! pp reco
// {15,928810,0.49235,5.30087e-07},
// {30,690508,0.030482,4.41443e-08},
// {50,956852,0.0035721,3.73318e-09},
// {80,945232,0.00042494,4.49562e-10},
// {120,791451,5.873e-05,7.42055e-11},
// {170,631722,9.199e-06,1.45618e-11},
// {220,307586,2.2564e-06,7.33583e-12},
// {280,1.2073e+06,6.336e-07,5.24806e-13},
// {370,511586,1.0884e-07,2.1275e-13},
// {460,783032,2.215e-08,2.82875e-14},
// {540,1.29072e+06,1.001e-08,7.75537e-15}


//! pp Reco w/o pT cut
// {15,930691  ,0.49235    ,5.290155e-07},
// {30,867066  ,0.030482   ,3.515534e-08},
// {50,978072  ,0.0035721  ,3.652185e-09},
// {80,950280  ,0.00042494 ,4.471735e-10},
// {120,876776 ,5.873e-05  ,6.698404e-11},
// {170,803916 ,9.199e-06  ,1.144274e-11},
// {220,947250 ,2.2564e-06 ,2.382053e-12},
// {280,1074188,6.336e-07  ,5.898409e-13},
// {370,895608 ,1.0884e-07 ,1.215264e-13},
// {460,863923 ,2.215e-08  ,2.563886e-14},
// {540,1327230,1.001e-08  ,7.5420236e-15}

//! HI reco (actually was with pp tracking only)
// {15,929040,0.49235,5.29956e-07},
// {30,356393,0.030482,8.55292e-08},
// {50,679355,0.0035721,5.25808e-09},
// {80,913874,0.00042494,4.64988e-10},
// {120,128147,5.873e-05,4.58302e-10},
// {170,539051,9.199e-06,1.70652e-11},
// {220,880261,2.2564e-06,2.56333e-12},
// {280,244502,6.336e-07,2.59139e-12},
// {370,706424,1.0884e-07,1.54072e-13},
// {460,560562,2.215e-08,3.95139e-14},
// {540,651991,1.001e-08,1.5353e-14}


//! HI reco with HI tracking
  // {15 ,1868518  ,0.492350  ,2.634976e-07},
  // {30 ,999064   ,0.030482  ,3.051056e-08},
  // {50 ,994627   ,0.003572  ,3.591396e-09},
  // {80 ,973474   ,0.00042494,4.365191e-10},
  // {120,965725   ,5.873e-05 ,6.081441e-11},
  // {170,2170421  ,9.199e-06 ,4.238348e-12},
  // {220,1364257  ,2.2564e-06,1.653941e-12},
  // {280,1197433  ,6.336e-07 ,5.291319e-13},
  // {370,1173556  ,1.3099e-07,1.116180e-13},
  // {540,1090925  ,1.001e-08 ,9.175699e-15}



  //! pp with new RECO 757p1 for hcal and ecal timing slews
  // {15,1046617 ,0.492350  ,4.7042041e-07},
  // {30,930294  ,0.030482  ,3.2765986e-08},
  // {50,1091974 ,0.003572  ,3.2712317e-09},
  // {80,1087400 ,0.00042494,3.90785360e-10},
  // {120,1081058,5.873e-05 ,5.43264099e-11},
  // {170,795845 ,9.199e-06 ,1.15587834e-11},
  // {220,722028 ,2.2564e-06,3.12508656e-12},
  // {280,1118212,6.336e-07 ,5.66618852e-13},
  // {370,1048958,1.0884e-07,1.03760112e-13},
  // {460,991650 ,2.215e-08 ,2.23365099e-14},
  // {540,1425964,1.001e-08 ,7.01981256e-15}


  //! pp with new HI RECO 757p1 for hcal and ecal timing slews v4_mc_00
  {15,1046429,0.49235    ,4.70504927e-07},
  {30,930110,0.030482    ,3.27724680e-08},
  {50,1091940,0.003572   ,3.27133359e-09},
  {80,1087095,0.00042494 ,3.90894999e-10},
  {120,1080993,5.873e-05 ,5.43296765e-11},
  {170,795680,9.199e-06  ,1.156118040e-11},
  {220,931164,2.2564e-06 ,2.423203646e-12},
  {280,1176166,6.336e-07 ,5.386994693e-13},
  {370,1059159,1.0884e-07,1.027607753e-13},
  {460,992914,2.215e-08  ,2.230807502e-14},
  {540,1426350,1.001e-08 ,7.017912854e-15}

};


//float GetXsecwt(float/*maxpthat*/);
float GetXsecWt(float pthat);

int weight(std::string infile="/mnt/hadoop/cms/store/user/pawan/PYTHIA_QCD30_TuneCUETP8M1_cfi_RECODEBUGHI_757p1_HcalRespCorrs_v4_00_mc/6.root",
	   std::string outfile="HiForest_PYTHIA_QCD_30_TuneCUETP8M1_cfi_RECODEBUGHI_757p1_HcalRespCorrs_v4_00_mc_6.root",
	   std::string radii="4"

)
{

  TFile *fin = TFile::Open(infile.c_str(), "r");
  //cout<<infile<<endl;

  // cout<<endl;
  // cout<<endl;
  // fin->ls();
  // cout<<endl;
  // cout<<endl;

  const int ndir=4;
  string DirName[ndir]= {
    "ak"+radii+"CaloJetAnalyzer",
    "ak"+radii+"PFJetAnalyzer",
    "akPu"+radii+"CaloJetAnalyzer",
    "akPu"+radii+"PFJetAnalyzer"
    // "akVs"+radii+"CaloJetAnalyzer",
    // "akVs"+radii+"PFJetAnalyzer"
  };

  TFile *fout=new TFile(outfile.c_str(),"RECREATE");

  TTree *tr_in=0, *tr_out=0;
  TH1D *hpthat[ndir], *hgenpt[ndir], *hrawpt[ndir];

  for(Int_t idir=0;idir<ndir;idir++){

    cout <<"idir =" << idir <<" JetName ="<< DirName[idir].c_str() <<endl ;
    tr_in = (TTree*)fin->Get(Form("%s/t",DirName[idir].c_str()));
    
    //Declaration of leaves types
    int   nref;
    float pthat;
    float weight;
    float corrpt[1000];
    float jtpt[1000];
    float rawpt[1000];
    float jteta[1000];
    float jtphi[1000];
    float jty[1000];
    float jtpu[1000];
    // float neutralMax[1000];
    // float chargedMax[1000];
    // float photonMax[1000];
    // float neutralSum[1000];
    // float chargedSum[1000];
    // float photonSum[1000];
    float refpt[1000];
    float refeta[1000];
    float refphi[1000];
    float refdphijt[1000];       
    float refdrjt[1000];       
    float refparton_pt[1000];
    int refparton_flavor[1000];
    int subid[1000];
    

    tr_in->SetBranchAddress("nref",&nref);
    tr_in->SetBranchAddress("pthat",&pthat);
    tr_in->SetBranchAddress("rawpt",rawpt);
    tr_in->SetBranchAddress("jtpt",jtpt);
    tr_in->SetBranchAddress("jteta",jteta);
    tr_in->SetBranchAddress("jty",jty);
    tr_in->SetBranchAddress("jtphi",jtphi);
    tr_in->SetBranchAddress("jtpu",jtpu);
    // tr_in->SetBranchAddress("neutralMax",neutralMax);
    // tr_in->SetBranchAddress("chargedMax",chargedMax);
    // tr_in->SetBranchAddress("photonMax",photonMax);
    // tr_in->SetBranchAddress("neutralSum",neutralSum);
    // tr_in->SetBranchAddress("chargedSum",chargedSum);
    // tr_in->SetBranchAddress("photonSum",photonSum);
    tr_in->SetBranchAddress("refpt",refpt);
    tr_in->SetBranchAddress("refphi",refphi);
    tr_in->SetBranchAddress("refeta",refeta);
    tr_in->SetBranchAddress("refdphijt",refdphijt);
    tr_in->SetBranchAddress("refdrjt",refdrjt);
    tr_in->SetBranchAddress("refparton_pt",refparton_pt);
    tr_in->SetBranchAddress("refparton_flavor",refparton_flavor);
    tr_in->SetBranchAddress("subid",subid);
    cout<<"get jet trees!!! "<<endl;

    //! Add Friends to the TTree
    //tr_in->AddFriend(tr_hlt);
    //tr_in->AddFriend(tr_ev);

    fout->mkdir(DirName[idir].c_str());
    fout->cd(DirName[idir].c_str());
    
    tr_out = new TTree("t","Jet  Response Analyzer");
    tr_out->SetMaxTreeSize(4200000000);

    // Set output branch addresses.
    //tr_out->Branch("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA,"pPAcollisionEventSelectionPA/I");
    //tr_out->Branch("vz",&vz,"vz/F");
    //tr_out->Branch("hiBin",&hiBin,"hiBin/I");
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
    // tr_out->Branch("neutralMax",neutralMax,"neutralMax[nref]/F");
    // tr_out->Branch("chargedMax",chargedMax,"chargedMax[nref]/F");
    // tr_out->Branch("photonMax",photonMax,"photonMax[nref]/F");
    // tr_out->Branch("neutralSum",neutralSum,"neutralSum[nref]/F");
    // tr_out->Branch("chargedSum",chargedSum,"chargedSum[nref]/F");
    // tr_out->Branch("photonSum",photonSum,"photonSum[nref]/F");
    tr_out->Branch("refpt",refpt,"refpt[nref]/F");
    tr_out->Branch("refeta",refeta,"refeta[nref]/F");
    tr_out->Branch("refphi",refphi,"refphi[nref]/F");
    tr_out->Branch("refdphijt",refdphijt,"refdphijt[nref]/F");
    tr_out->Branch("refdrjt",refdrjt,"refdrjt[nref]/F");
    tr_out->Branch("refparton_pt",refparton_pt,"refparton_pt[nref]/F");
    tr_out->Branch("refparton_flavor",refparton_flavor,"refparton_flavor[nref]/I");
    tr_out->Branch("subid",subid,"subid[nref]/I");

    hpthat[idir] = new TH1D(Form("pthhat_%d",idir),Form("pthat %s",DirName[idir].c_str()),400,0.,2000.);
    hpthat[idir]->Sumw2();
    hgenpt[idir] = new TH1D(Form("hgenpt_%d",idir),Form("genpt %s",DirName[idir].c_str()),400,0.,2000.);
    hgenpt[idir]->Sumw2();
    hrawpt[idir] = new TH1D(Form("hrawopt_%d",idir),Form("rawopt %s",DirName[idir].c_str()),400,0.,2000.);
    hrawpt[idir]->Sumw2();


    Long64_t nbytes = 0;
    cout<< tr_in->GetName() << " total entries : " << tr_in->GetEntries() << endl;
    
    for (Long64_t i=0; i<tr_in->GetEntries();i++) {
      nbytes += tr_in->GetEntry(i);

      weight = GetXsecWt(pthat);
      hpthat[idir]->Fill(pthat,weight);

      //! jet loop
      for(int iref=0;iref<nref;iref++){
	corrpt[iref] = jtpt[iref];
	jtpt  [iref] = rawpt[iref];

	// if( abs(refparton_flavor[iref]) <= 21 ){
	//   subid[iref]=0;
	if( subid[iref] == 0){
	  hgenpt[idir]->Fill(refpt[iref],weight);
	  hrawpt[idir]->Fill(rawpt[iref],weight);
	}
	// }else subid[iref]=-1;

      }//ref for loop

      tr_out->Fill();
    }//nentries loop

    tr_out->Write();
    hpthat[idir]->Write();
    hgenpt[idir]->Write();
    hrawpt[idir]->Write();
    fout->cd("../");
  }    
  cout <<"finish the code!!"<<endl ;
  fout->Close();
  cout<<"XXXX"<<endl;

  return 0;
}
// float GetXsecWt(float maxpt)
// {
//   float effxsec=0;
//   for(int i=0; i<11; i++){
//     //for(int i=0; i<10; i++){
//     if(fabs(maxpt - xsec[i][2]) < 1e-08){
//       effxsec = xsec[i][0] - xsec[i+1][0];
//       //effxsec = xsec[i][0];
//       return effxsec;
//     }
//   }
//   return  1;
// }
float GetXsecWt(float pthat)
{
  float wt=1.0;
  for( int i=0; i<npthat; i++){
    if( pthat >  pthatwt_pp[i][0] )wt = pthatwt_pp[i][3];
  }
  return wt;
}
