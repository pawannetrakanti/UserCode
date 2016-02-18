#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <utility>

#include <TSystem.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStopwatch.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include "HiForest/hiForest.h"


using namespace std;

#ifdef _MAKECINT_
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ class Jets+;
#pragma link C++ class Hlts+;
#pragma link C++ class Skims+;
#pragma link C++ class EventTree+;
#endif

#define zmass 91.1876

const float kvzcut=24.;

// const float kmuetacut=2.4;
// const float kmuptcut=20.;
// const float kzptcut=40.;
// const float kzmassmin=70.;
// const float kzmassmax=110.;

const float kjetetacut=2.4;
const float kjetptcut=30.;
const float kdrcut=0.3;
float kdelrcut=0.4;

// struct Muon{
//   TLorentzVector muVec;
//   int much;
// };
// struct DiMuon{
//   float OppSign;
//   float M;
//   float Pt;
//   float Eta;
//   float Rap;
//   float Phi;
//   float E;
// };
// struct Jet{
//   float Pt;
//   float Eta;
//   float Phi;
//   float E;
// };

float deltaR(float /*eta1*/, float /*phi1*/, float /*eta2*/, float /*phi2*/);
float deltaEta(float /*eta1*/, float /*eta2*/);
float deltaPhi(float /*phi1*/, float /*phi2*/);

// typedef std::pair< DiMuon, Jet > ZJetPair;
// struct CompareZJetPairs {
//   bool operator()(const ZJetPair &A1, const ZJetPair &A2){
//     DiMuon dimuOne = A1.first;
//     Jet    jetOne  = A1.second;

//     DiMuon dimuTwo = A2.first;
//     Jet    jetTwo  = A2.second;
    
//     // float delrOne = deltaR( dimuOne.dimuEta, dimuOne.dimuPhi, jetOne.jetEta, jetOne.jetPhi );
//     // float delrTwo = deltaR( dimuTwo.dimuEta, dimuTwo.dimuPhi, jetTwo.jetEta, jetTwo.jetPhi );

//     float delPhiOne = deltaPhi( dimuOne.Phi, jetOne.Phi );
//     float delPhiTwo = deltaPhi( dimuTwo.Phi, jetTwo.Phi );

//     return ( fabs( delPhiOne - pi ) < fabs( delPhiTwo - pi ));
//   }
// };
// typedef std::multiset< ZJetPair, CompareZJetPairs > ZJetCollection;
// typedef std::multiset< ZJetPair >::iterator ZJItr;


const int npthat=11;
float pthatwt_pp[npthat][4] ={
  {15,932778,0.49235    ,5.278319171e-07},
  {30,903567,0.030482   ,3.373518510e-08},
  {50,983531,0.0035721  ,3.631913992e-09},
  {80,1820782,0.00042494,2.333832387e-10},
  {120,1080554,5.873e-05,5.435174920e-11},
  {170,836152,9.199e-06 ,1.100158823e-11},
  {220,954396,2.2564e-06,2.364217790e-12},
  {280,1083994,6.336e-07,5.84505080e-13},
  {370,948240,1.0884e-07,1.147810683e-13},
  {460,1558268,2.215e-08,1.421449971e-14},
  {540,2597338,1.001e-08,3.85394585e-15}

};
float  GetXsecWt(float pthat);

bool is_file(const char *fileName);
bool is_file(const char *fileName)
{
  std::ifstream infile(fileName,ios::in);
  return infile.good();
}

TStopwatch timer;
int bookJetTree(std::string runAlgo="ak4PF",
		std::string inputFile="/mnt/hadoop/cms/store/user/abaty/mergedForests/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV_ppSignal/Pythia8_Dijet80_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/0.root",
		std::string outputFile="test.root",
		std::string kDataType="mc"
)
{
  timer.Start();
  //gSystem->Load("HiForest/EventTree_cxx.so");
  //gSystem->Load("HiForest/hiForest_cxx.so");

  bool debug=false;

  int  ispp=1;
  bool ismc=1;
  bool ishireco=0;
  if( kDataType == "mc" )ismc=1;
  HiForest *eve = new HiForest("QGAna",inputFile.c_str(), ispp, ishireco, ismc);

  std::cout << "# of entries : " << eve->GetEntries() <<  std::endl;
  Long64_t nEntries = ( debug == true ) ? 1000 : eve->GetEntries();


  //eve->hasHltTree=0;
  //eve->hasPFCandTree=0;

  const int njetTrees=6;

  //! Only reading ak4PF jets
  //! 0 : ak3PF, 1 : ak4PF, 2 : ak5PF, 3 : akPu3PF, 4 : akPu4PF, 5 : akPu5PF
  for(int i=0; i<njetTrees; i++){
    eve->hasJetTree[i]=0; //! ak3PF
  }

  if( runAlgo == "ak3PF" ){
    eve->hasJetTree[0]=1; //! ak3PF
  }else if( runAlgo ==  "ak4PF" ){
    eve->hasJetTree[1]=1; //! ak4PF
  }else if( runAlgo == "ak5PF" ){
    eve->hasJetTree[2]=1; //! ak5PF
  }else if( runAlgo == "akPu3PF" ){
    eve->hasJetTree[3]=1; //! akPu3PF
  }else if ( runAlgo == "akPu4PF" ){
    eve->hasJetTree[4]=1; //! akPu4PF
  }else if( runAlgo == "akPu5PF" ){
    eve->hasJetTree[5]=1; //! akPu5PF
  }

  if( debug )std::cout << " All the ncecessary algos initialized " << std::endl;


  eve->hasTupleTree=0;

  eve->skimTree->SetBranchStatus("*",0,0);
  eve->skimTree->SetBranchStatus("pHBHENoiseFilterResultProducer",1);

  eve->trackTree->SetBranchStatus("*",0,0);
  eve->trackTree->SetBranchStatus("nVtx",1);  
  eve->trackTree->SetBranchStatus("zVtx",1);  

  if( debug )std::cout << " TrackTree branches set " << std::endl;

  //std::vector<std::string> tupleList;
  ///! Possible options for now are :
  ///!  Zmumu, Zee, Jet (Zmumu and Zee + Jet), BJet, Photon,  PhotonJet, MET
  //tupleList.push_back("Zmumu");
  //tupleList.push_back("Jet");
  //eve->SelectTupleBranches(tupleList);

  //! Jet Tree for discriminatnts
  float vz;
  float pthat;
  float weight;
  int njets;
  float jetpt [200];
  float jeteta[200];
  float jetphi[200];
  float jetmass[200];
  int jetflav [200];

  int nPFCand [200];
  int nCharged[200];
  int nNeutral[200];
  int nPFMult [200];
  float jetMByPt[200];
  float RMSCand  [200];
  float Axis1[200];
  float Axis2[200];
  float Sigma[200];
  float R   [200];
  float pTD [200];
  float pull[200];

  
  TFile *fout = new TFile(outputFile.c_str(),"RECREATE");
  TTree *jetTree = new TTree("jetTree","QG Discriminant Tree");
  jetTree->Branch("vz",&vz,"vz/F");
  jetTree->Branch("pthat",&pthat,"pthat/F");
  jetTree->Branch("weight",&weight,"weight/F");
  jetTree->Branch("njets",&njets,"njets/I");
  jetTree->Branch("jetpt",jetpt,"jetpt[njets]/F");
  jetTree->Branch("jeteta",jeteta,"jeteta[njets]/F");
  jetTree->Branch("jetphi",jetphi,"jetphi[njets]/F");
  jetTree->Branch("jetmass",jetmass,"jetmass[njets]/F");
  jetTree->Branch("jetflav",jetflav,"jetflav[njets]/I");
  jetTree->Branch("nPFCand",nPFCand,"nPFCand[njets]/I");  
  jetTree->Branch("nCharged",nCharged,"nCharged[njets]/I");  
  jetTree->Branch("nNeutral",nNeutral,"nNeutral[njets]/I");  
  jetTree->Branch("nPFMult",nPFMult,"nPFMult[njets]/I");  
  jetTree->Branch("jetMByPt",jetMByPt,"jetMByPt[njets]/F");  
  jetTree->Branch("RMSCand",RMSCand,"RMSCand[njets]/F");  
  jetTree->Branch("Axis1",Axis1,"Axis1[njets]/F");  
  jetTree->Branch("Axis2",Axis2,"Axis2[njets]/F");  
  jetTree->Branch("Sigma",Sigma,"Sigma[njets]/F");  
  jetTree->Branch("R",R,"R[njets]/F");  
  jetTree->Branch("pTD",pTD,"pTD[njets]/F");  
  jetTree->Branch("pull",pull,"pull[njets]/F");  

  Jets *iJet=0;
  Long64_t nb=0;
  for(int iev=0; iev < nEntries; iev++){
    nb += eve->GetEntry(iev);
    //if( debug )std::cout<<"Processing " << iev << std::endl;
    if( iev % 10000 == 0 )std::cout<<"Processing " << iev << std::endl;
    //std::cout << iev << " muon Size : " << eve->MuPt->size() << " electron size : " << eve->ElPt->size() << " photon size " << eve->PhotPt->size() << "  jet size : " << eve->JetA//k04Pt->size() << std::endl;
    //if( eve->MuPt->size() != 2 ) continue;

    if( !eve->skim.pHBHENoiseFilterResultProducer && kDataType=="data" ) continue;
    if( eve->track.nVtx == 0 )continue;
    if( fabs( eve->track.zVtx[0] ) > kvzcut )continue;

    if( eve->pfcand.nPFpart < 2 )continue;

    vz    = eve->track.zVtx[0];
    //std::cout<<"Processing " << iev << std::endl;

    iJet  = eve->GetJetByAlgo(runAlgo.c_str());
    if( kDataType == "mc" ){
      pthat  = iJet->pthat;
      weight = GetXsecWt(pthat);
    }else weight=1;

    njets=0;
    for( int ig=0; ig < iJet->ngen; ig++ ){
      if( iJet->genmatchindex[ig] < 0 )continue;
      int ij = iJet->genmatchindex[ig];

      float recopt = iJet->jtpt[ij];
      float jtmass = iJet->jtm[ij];
      if( kDataType == "mc" ){
	if( iJet->refpt[ij] < 0 || abs(iJet->refdrjt[ij]) > kdrcut )continue; 
	if( abs(iJet->refparton_flavor[ij]) > 21 )continue;
      }else{
	if( fabs(iJet->refeta[ij]) > kjetetacut || 
	    recopt < kjetptcut )continue;      
	//iJet->jtpt[ij] < kjetptcut )continue;      
      }
      
      jetpt [njets]  = recopt; //iJet->jtpt[ij];
      jeteta[njets]  = iJet->jteta[ij];
      jetphi[njets]  = iJet->jtphi[ij];
      jetmass[njets] = jtmass;
      if( kDataType == "mc" )jetflav[njets] = iJet->refparton_flavor[ij];
      else jetflav[njets] = 0;

      TLorentzVector jetvec;
      jetvec.SetPtEtaPhiE(iJet->jtpt[ij], iJet->jteta[ij], iJet->jtphi[ij], iJet->jtenergy[ij]);
      float jety = jetvec.Rapidity();

      nPFCand[njets]=0;
      nCharged[njets]=0;
      nNeutral[njets]=0;

      nPFMult[njets]=0;
      int nCh=0;
      int nEl=0;
      int nMu=0;
      int nPh=0;
      int nNh=0;

      float tot_wt=0;

      RMSCand[njets]=0;
      float rmscand_n=0;
      float rmscand_d=0;

      //! Axes
      float M11=0,M22=0,M12=0,M21=0;
      Axis1[njets]=0;
      Axis2[njets]=0;
      Sigma[njets]=0;

      float maxCandPt=-999;
      float sumCandPt=0;
      float sumCandPtSq=0;
      R[njets]=0;

      pTD[njets]=0;
      
      TVector2 t_Vect(0,0);
      TVector2 r(0,0);
      pull[njets]=0;

      for(int ipf=0; ipf < eve->pfcand.nPFpart; ipf++){
	if( eve->pfcand.pfId[ipf] == 0 )continue;
	if( eve->pfcand.pfPt[ipf] < 1.0 )continue;
	//if( fabs(eve->pfcand.pfEta[ipf] > kjetetacut )continue;

	float dr = deltaR(eve->pfcand.pfEta[ipf], eve->pfcand.pfPhi[ipf], iJet->jteta[ij], iJet->jtphi[ij]);
	if( fabs( dr ) > kdelrcut ) continue;

	float deta = deltaEta(eve->pfcand.pfEta[ipf], iJet->jteta[ij]);
	float dphi = deltaPhi(eve->pfcand.pfPhi[ipf], iJet->jtphi[ij]);
	
	if( eve->pfcand.pfPt[ipf] > maxCandPt )maxCandPt=eve->pfcand.pfPt[ipf];
	sumCandPt   += eve->pfcand.pfPt[ipf];
	sumCandPtSq += pow(eve->pfcand.pfPt[ipf],2);

	nPFCand[njets]++;
	switch( eve->pfcand.pfId[ipf] ){
	case 1 : 
	  nCh++; break;
	case 2 : 
	  nEl++; break;
	case 3 :
	  nMu++; break;
	case 4 :
	  nPh++; break;
	case 5 :
	  nNh++; break;
	default :
	  break;
	}
	
	float wt = eve->pfcand.pfPt[ipf] * eve->pfcand.pfPt[ipf];
	tot_wt += wt;

	//! RMSCand
	rmscand_n += wt*dr*dr;
	rmscand_d += wt; 
	
	M11 += wt*deta*deta;
	M22 += wt*dphi*dphi;
	M12 += wt*deta*dphi;
	M21 += wt*deta*dphi;

	//! Pull
	TLorentzVector pfvec;
	pfvec.SetPtEtaPhiE(eve->pfcand.pfPt[ipf], eve->pfcand.pfEta[ipf], eve->pfcand.pfPhi[ipf], eve->pfcand.pfEnergy[ipf]);
	float pfy = pfvec.Rapidity();
	float dy = pfy - jety;
	r.Set( dy, dphi );
	float r_mag = r.Mod();
	t_Vect += ( eve->pfcand.pfPt[ipf] /  recopt ) * r_mag *r;
	//t_Vect += ( eve->pfcand.pfPt[ipf] /  iJet->jtpt[ij] ) * r_mag *r;

      }//! PF cand loop ends

      bool noPFCand = (tot_wt==0);
      if( noPFCand )continue;

      nPFMult[njets] = min(9,nMu)
	+ 10*min(99,nCh)
	+ 1000*min(99,nNh)
	+ 100000*min(99,nEl)
	+ 10000000*min(99,nPh);
      
      nCharged[njets] = (nCh + nMu + nEl);
      nNeutral[njets] = (nNh + nPh);
      
      if( min(9,nMu)  == 9  )std::cout << iev << " jetpt : " << recopt << " More Muons    in this jet " << nMu <<std::endl;
      if( min(99,nCh) == 99 )std::cout << iev << " jetpt : " << recopt << " More Charge   in this jet " << nCh <<std::endl;
      if( min(99,nNh) == 99 )std::cout << iev << " jetpt : " << recopt << " More Neutral  in this jet " << nNh <<std::endl;
      if( min(99,nEl) == 99 )std::cout << iev << " jetpt : " << recopt << " More Electron in this jet " << nEl <<std::endl;
      if( min(99,nPh) == 99 )std::cout << iev << " jetpt : " << recopt << " More Photons  in this jet " << nPh <<std::endl;
      
      RMSCand[njets] = sqrt ( rmscand_n / rmscand_d );
      
      M12 = -1.*M12;
      M21 = -1.*M21;

      //! eign values
      float trace = M11 + M22;
      float detrm = (M11*M22) - (M12*M21);
      
      float lam1 = trace/2. + sqrt( pow(trace,2)/4. - detrm );
      float lam2 = trace/2. - sqrt( pow(trace,2)/4. - detrm );
      
      Axis1[njets] = sqrt( lam1 / tot_wt );
      Axis2[njets] = sqrt( lam2 / tot_wt );

      Sigma[njets] = sqrt( pow(Axis1[njets],2) + pow(Axis2[njets],2) );
      
      R[njets] = maxCandPt / sumCandPt;
      
      pTD[njets] = sqrt( sumCandPtSq ) / sumCandPt;
      
      pull[njets] = t_Vect.Mod();
      
      jetMByPt[njets] = jtmass / recopt;
      
      njets++;
    }//! jet loop ends
    if( njets > 0 )jetTree->Fill();
  }

  fout->cd();
  fout->Write();
  fout->Close();

  timer.Stop();
  float mbytes = 1.0e-06*nb;
  float rtime = timer.RealTime();
  float ctime = timer.CpuTime();
  std::cout<<std::endl;
  std::cout<<Form("RealTime=%f seconds, CpuTime=%f seconds",rtime,ctime)<<std::endl;
  std::cout<<Form("You read %f Mbytes/RealTime seconds",mbytes/rtime)<<std::endl;
  std::cout<<Form("You read %f Mbytes/CpuTime seconds",mbytes/ctime)<<std::endl;
  std::cout<<std::endl;
  std::cout<<"Good bye : " <<"\t"<<std::endl;

  return 0;

}
float deltaEta(float eta1, float eta2)
{
  float deta = eta1 - eta2;
  return deta;
}
float deltaPhi(float phi1, float phi2)
{
  float dphi = fabs(phi1 - phi2);
  if(dphi > pi)dphi -= 2*pi;
  return dphi;
}
float deltaR(float eta1, float phi1, float eta2, float phi2)
{
  float deta = eta1 - eta2;
  float dphi = deltaPhi(phi1, phi2);
  float dr = sqrt(pow(deta,2) + pow(dphi,2));
  return dr;
}
float GetXsecWt(float pthat)
{
  float wt=1.0;
  for( int i=0; i<npthat; i++){
    if( pthat >  pthatwt_pp[i][0] )wt = pthatwt_pp[i][3];
  }
  return wt;
}
