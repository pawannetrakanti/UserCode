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

const float kvzcut=24.;
const float kjetetacut=2.4;
const float kjetptcut=15.;
const float kdrcut=0.3;
float kdelrcut=0.4;


float deltaR(float /*eta1*/, float /*phi1*/, float /*eta2*/, float /*phi2*/);
float deltaEta(float /*eta1*/, float /*eta2*/);
float deltaPhi(float /*phi1*/, float /*phi2*/);
int GetPtBin(float /*pt*/);

//! pt binning                                                                                             
double ptbins[] ={17, 22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 
		  196, 220, 245, 272, 300, 400, 550, 790, 1000
};
const int nptbins = sizeof(ptbins)/sizeof(double) - 1;


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

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();



  bool debug=false;
  bool fillTree=false;

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
    eve->hasJetTree[i]=0; //! 
  }
  if( runAlgo == "ak3PF" ){
    eve->hasJetTree[0]=1; //! ak3PF
    kdelrcut=0.3;
  }else if( runAlgo ==  "ak4PF" ){
    eve->hasJetTree[1]=1; //! ak4PF
    kdelrcut=0.4;
  }else if( runAlgo == "ak5PF" ){
    eve->hasJetTree[2]=1; //! ak5PF
    kdelrcut=0.5;
  }else if( runAlgo == "akPu3PF" ){
    eve->hasJetTree[3]=1; //! akPu3PF
  }else if ( runAlgo == "akPu4PF" ){
    eve->hasJetTree[4]=1; //! akPu4PF
  }else if( runAlgo == "akPu5PF" ){
    eve->hasJetTree[5]=1; //! akPu5PF
  }


  if( debug )std::cout << " All the ncecessary algos initialized " << std::endl;

  //eve->hasTupleTree=0;

  eve->skimTree->SetBranchStatus("*",0,0);
  eve->skimTree->SetBranchStatus("pHBHENoiseFilterResultProducer",1);

  eve->trackTree->SetBranchStatus("*",0,0);
  eve->trackTree->SetBranchStatus("nVtx",1);  
  eve->trackTree->SetBranchStatus("zVtx",1);  

  if( debug )std::cout << " TrackTree branches set " << std::endl;

  //! Jet Tree for discriminatnts
  float vz;
  float pthat;
  float weight;
  int njets;
  float jetpt [500];
  float jeteta[500];
  float jetphi[500];
  float jetmass[500];
  //float jetarea[500];
  int jetflav  [500];

  int nPFCand [500];
  int nCharged[500];
  int nNeutral[500];
  int nPFMult [500];
  float jetMByPt[500];
  float RMSCand [500];
  float Axis1[500];
  float Axis2[500];
  float Sigma[500];
  float R   [500];
  float pTD [500];
  float pull[500];

  
  TFile *fout = new TFile(outputFile.c_str(),"RECREATE");
  fout->mkdir(Form("%sJetAnalyzer",runAlgo.c_str()));
  fout->cd(Form("%sJetAnalyzer",runAlgo.c_str()));
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


  //! Histograms
  TH1D *hPull[2][nptbins], *hPtD[2][nptbins], *hRMSCand[2][nptbins], *hSigma[2][nptbins],
    *hR[2][nptbins], *hAxis1[2][nptbins], *hAxis2[2][nptbins],  *hJetMByPt[2][nptbins];

  for(int i=0; i<2; i++){
    for(int k=0; k<3; k++){
      for(int j=0; j<nptbins; j++){
	hRMSCand[i][k][j] = new TH1D(Form("hRMSCand_%s_%d_%d_%d",runAlgo.c_str(),i,k,j),"RMS Cand ",200,0,5);    
	hAxis1[i][k][j]   = new TH1D(Form("hAxis1_%s_%d_%d_%d"  ,runAlgo.c_str(),i,k,j),"Axis1",200,0,5); 
	hAxis2[i][k][j]   = new TH1D(Form("hAxis2_%s_%d_%d_%d"  ,runAlgo.c_str(),i,k,j),"Axis2",200,0,5); 
	hSigma[i][k][j]   = new TH1D(Form("hSigma_%s_%d_%d_%d"  ,runAlgo.c_str(),i,k,j),"Sigma",500,0,1); 
	hPull[i][k][j]    = new TH1D(Form("hPull_%s_%d_%d_%d"   ,runAlgo.c_str(),i,k,j),"Pull ",500,0,5);    
	hR[i][k][j]       = new TH1D(Form("hR_%s_%d_%d_%d"      ,runAlgo.c_str(),i,k,j),"R ",500,0,5.0);  
	hPtD[i][k][j]     = new TH1D(Form("hPtD_%s_%d_%d_%d"    ,runAlgo.c_str(),i,k,j),"p_{T}D",500,0,5.0);
	hJetMByPt[i][k][j] = new TH1D(Form("hJetMByPt_%s_%d_%d_%d",runAlgo.c_str(),i,k,j),"JetMass/Jet pT",200,0,1.0);
      }
    }
  }
  fout->cd("../");

  if( debug )cout << " Histograms defined " << endl;

  Jets *iJet=0;
  iJet  = eve->GetJetByAlgo(runAlgo.c_str());
  if( debug )cout<<" Got Jet algo " <<endl;
  
  EventTree *itupel=0;
  itupel = eve->GetTupelTree();
  if( debug )cout<<" Got Tupel " <<endl;

  Long64_t nb=0;
  for(int iev=0; iev < nEntries; iev++){
    nb += eve->GetEntry(iev);
    if( debug )std::cout<<"Processing " << iev << std::endl;
    if( iev % 10000 == 0 )std::cout<<"Processing " << iev << std::endl;

    if( debug )cout << " # of jets in tupel " <<itupel->JetAk04Pt->size() 
		    << " # of jets in jettree : " << iJet->nref <<endl;

    if( !eve->skim.pHBHENoiseFilterResultProducer && kDataType=="data" ) continue;
    if( eve->track.nVtx == 0 )continue;
    if( fabs( eve->track.zVtx[0] ) > kvzcut )continue;

    vz    = eve->track.zVtx[0];
    //std::cout<<"Processing " << iev << std::endl;

    if( kDataType == "mc" ){
      pthat  = iJet->pthat;
      weight = GetXsecWt(pthat);
    }else weight=1;


    njets=0;
    for( int ij=0; ij < (int)itupel->JetAk04Pt->size(); ij++ ){

      int imatch = itupel->JetAk04GenJet->at(ij);
      //if( debug ) cout << "  imatch : " << imatch << endl;
      if( imatch < 0 )continue;

      float recoeta = itupel->JetAk04Eta->at(ij);
      float recopt  = itupel->JetAk04Pt->at(ij);
      float recophi = itupel->JetAk04Phi->at(ij);      
      float recoen  = itupel->JetAk04E->at(ij);      
      float partflav= itupel->JetAk04PartFlav->at(ij);
      
      if ( fabs (recoeta) > kjetetacut || ( recopt < kjetptcut ) )continue;

      //! Jet ID 
      //if(chf>0 && nhf<0.99 && cmult>0.0 && cemf<0.99 && nemf<0.99 && nconst>1) tempJetID=1;
      // float jetid =  itupel->JetAk04Id->at(ij);
      // if( debug ) cout << "  jetid : " << jetid << endl;
      // if( jetid != 1 )continue;
      
      TLorentzVector jetvec;
      jetvec.SetPtEtaPhiE(recopt,recoeta,recophi,recoen);
      
      float jtmass = jetvec.M();
      if( kDataType == "mc" ){
    	if( fabs(partflav) > 21 )continue;
      }

      int iFlav=-1;
      if( fabs(partflav) < 6 )iFlav=1;
      if( fabs(partflav) == 21 )iFlav=2;

      jetpt [njets]  = recopt;
      jeteta[njets]  = recoeta;
      jetphi[njets]  = recophi;
      jetmass[njets] = jtmass;
      if( kDataType == "mc" )jetflav[njets] = (int)partflav;
      else jetflav[njets] = 0;
      
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

      int nConstituents = itupel->JetAk04ConstCnt->at(ij);

      for(int ipf=0; ipf < nConstituents; ipf++){
	//cout << itupel->JetAk04ConstPt->at(ipf) <<endl;
    	if(  itupel->JetAk04ConstId->at(ipf) == 0 )continue;
    	if(  itupel->JetAk04ConstPt->at(ipf) < 1.0 )continue;

    	int   pfid  = fabs(itupel->JetAk04ConstId->at(ipf));
	



	

    	float pfeta = itupel->JetAk04ConstEta->at(ipf);
    	float pfphi = itupel->JetAk04ConstPhi->at(ipf);
    	float pfpt  = itupel->JetAk04ConstPt->at(ipf);
    	float pfen  = itupel->JetAk04ConstE->at(ipf);
    	//if( fabs(pfeta) > kjetetacut )continue;

    	float deta = deltaEta(pfeta, recoeta);
    	float dphi = deltaPhi(pfphi, recophi);
	
    	if( pfpt > maxCandPt )maxCandPt=pfpt;
    	sumCandPt   += pfpt;
    	sumCandPtSq += pow(pfpt,2);
	
    	nPFCand[njets]++;
    	float wt = pow(pfpt,2);
    	tot_wt += wt;
	
    	//! RMSCand
    	float dr = deltaR(pfeta, pfphi, recoeta, recophi);
    	rmscand_n += wt*dr*dr;
    	rmscand_d += wt; 
	
    	M11 += wt*deta*deta;
    	M22 += wt*dphi*dphi;
    	M12 += wt*deta*dphi;
    	M21 += wt*deta*dphi;
	
    	//! Pull
    	TLorentzVector pfvec;
    	pfvec.SetPtEtaPhiE(pfpt, pfeta, pfphi, pfen);
    	float pfy = pfvec.Rapidity();
    	float dy = pfy - jety;
    	r.Set( dy, dphi );
    	float r_mag = r.Mod();
    	t_Vect += ( pfpt /  recopt ) * r_mag *r;
      }//! PF cand loop ends

      bool noPFCand = (tot_wt==0);
      if( noPFCand )continue;
      
      int iBin = GetPtBin(recopt);
      
      //if( debug )cout << " \t tupel " << ij << "\t" << recopt << " nConst :  "<< nConstituents << " \t iBin " << iBin <<  "  njets " << njets  <<endl;
      
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
      
      if( iBin >=0 && iBin < nptbins){
	hRMSCand [0][0][iBin]->Fill(RMSCand[njets],weight);
       	hAxis1   [0][0][iBin]->Fill(Axis1[njets],weight);
       	hAxis2   [0][0][iBin]->Fill(Axis2[njets],weight);
       	hSigma   [0][0][iBin]->Fill(Sigma[njets],weight);
      	hR       [0][0][iBin]->Fill(R[njets],weight);
      	hPtD     [0][0][iBin]->Fill(pTD[njets],weight);
      	hPull    [0][0][iBin]->Fill(pull[njets],weight);
      	hJetMByPt[0][0][iBin]->Fill(jetMByPt[njets],weight);

	if( iFlav>0){
	  hRMSCand [0][iFlav][iBin]->Fill(RMSCand[njets],weight);
	  hAxis1   [0][iFlav][iBin]->Fill(Axis1[njets],weight);
	  hAxis2   [0][iFlav][iBin]->Fill(Axis2[njets],weight);
	  hSigma   [0][iFlav][iBin]->Fill(Sigma[njets],weight);
	  hR       [0][iFlav][iBin]->Fill(R[njets],weight);
	  hPtD     [0][iFlav][iBin]->Fill(pTD[njets],weight);
	  hPull    [0][iFlav][iBin]->Fill(pull[njets],weight);
	  hJetMByPt[0][iFlav][iBin]->Fill(jetMByPt[njets],weight);
	}
      }
      njets++;
    }//! jet loop
    if( njets>0 && fillTree )jetTree->Fill();

    // //! ref jets
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

      int iFlav=-1;
      if( abs(iJet->refparton_flavor[ij]) < 6 )iFlav=1;
      if( abs(iJet->refparton_flavor[ij]) == 21 )iFlav=2;


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
        switch( abs(eve->pfcand.pfId[ipf]) ){
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
        pfvec.SetPtEtaPhiE(eve->pfcand.pfPt[ipf], eve->pfcand.pfEta[ipf], eve->pfcand.pfPhi[ipf], 
    			   eve->pfcand.pfEnergy[ipf]);
        float pfy = pfvec.Rapidity();
        float dy = pfy - jety;
    	r.Set( dy, dphi );
        float r_mag = r.Mod();
        t_Vect += ( eve->pfcand.pfPt[ipf] /  recopt ) * r_mag *r;
      }//! PF cand loop ends
      bool noPFCand = (tot_wt==0);
      if( noPFCand )continue;

      int iBin = GetPtBin(recopt);

      nPFMult[njets] = min(9,nMu)
        + 10*min(99,nCh)
    	+ 1000*min(99,nNh)
        + 100000*min(99,nEl)
        + 10000000*min(99,nPh);

      nCharged[njets] = (nCh + nMu + nEl);
      nNeutral[njets] = (nNh + nPh);
      if( min(9,nMu)  == 9  )std::cout << iev << " jetpt : " << recopt << " More Muons    in this jet " 
    				       << nMu <<std::endl;
      if( min(99,nCh) == 99 )std::cout << iev << " jetpt : " << recopt << " More Charge   in this jet " 
    				       << nCh <<std::endl;
      if( min(99,nNh) == 99 )std::cout << iev << " jetpt : " << recopt << " More Neutral  in this jet " 
    				       << nNh <<std::endl;
      if( min(99,nEl) == 99 )std::cout << iev << " jetpt : " << recopt << " More Electron in this jet "
    				       << nEl <<std::endl;
      if( min(99,nPh) == 99 )std::cout << iev << " jetpt : " << recopt << " More Photons  in this jet " 
    				       << nPh <<std::endl;

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

      if( iBin >=0 && iBin <nptbins){
    	hRMSCand[1][0][iBin]->Fill(RMSCand[njets],weight);
    	hAxis1[1][0][iBin]->Fill(Axis1[njets],weight);
    	hAxis2[1][0][iBin]->Fill(Axis2[njets],weight);
    	hSigma[1][0][iBin]->Fill(Sigma[njets],weight);
    	hR    [1][0][iBin]->Fill(R[njets],weight);
    	hPtD  [1][0][iBin]->Fill(pTD[njets],weight);
    	hPull [1][0][iBin]->Fill(pull[njets],weight);
    	hJetMByPt [1][0][iBin]->Fill(jetMByPt[njets],weight);

	if( iFlav>0 ){
	  hRMSCand [1][iFlav][iBin]->Fill(RMSCand[njets],weight);
	  hAxis1   [1][iFlav][iBin]->Fill(Axis1[njets],weight);
	  hAxis2   [1][iFlav][iBin]->Fill(Axis2[njets],weight);
	  hSigma   [1][iFlav][iBin]->Fill(Sigma[njets],weight);
	  hR       [1][iFlav][iBin]->Fill(R[njets],weight);
	  hPtD     [1][iFlav][iBin]->Fill(pTD[njets],weight);
	  hPull    [1][iFlav][iBin]->Fill(pull[njets],weight);
	  hJetMByPt[1][iFlav][iBin]->Fill(jetMByPt[njets],weight);
	}


      }
      njets++;
      if( debug )cout << iev <<"\t Ended jet "<< njets << endl;
    }//! jet loop ends
  }

  fout->cd(Form("%sJetAnalyzer",runAlgo.c_str()));
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
int GetPtBin(float pt)
{
  if(pt<ptbins[0] || pt>ptbins[nptbins])return -1;
  for(int ix=0;ix<nptbins;ix++){
    if(pt>=ptbins[ix] && pt<ptbins[ix+1]){
      return ix;
    }
  }
  return -1;
}
