// -*- C++ -*-
//
// Package:    HOAnalyzer
// Class:      HOAnalyzer
// 
/**\class HOAnalyzer HOAnalyzer.cc UseOfHOTower/HOAnalyzer/src/HOAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Pawan Kumar Netrakanti,40 4-B28,+41227671589,
//         Created:  Wed Apr 23 12:03:42 CEST 2014
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/JetMatching/interface/JetMatchedPartons.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"

#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/JetID.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "FWCore/Framework/interface/GenericHandle.h"



#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "CLHEP/Vector/LorentzVector.h"

#include <string>
#include <iostream>
#include <fstream>

using namespace std;
using namespace edm;
using namespace pat;
using namespace reco;
using namespace CLHEP;

//
// class declaration
//


//
// constants, enums and typedefs
//
const int nchnmx = 10;
//
// static data member definitions
//
static const int nphimx = 72;
static const int njetmx = 50;
static const int  hcmx  = 1000;
static const int  vtmax = 20;

const int netahotower=15;
float etahotowers[netahotower+1]={0.000,  0.086,  0.173,  0.259,  0.341,  0.431,  0.518,  0.605,  0.692,  0.779,  0.887,  0.958,  1.044,  1.131,  1.218, 1.305/*1.262*/
};



class HOAnalyzer : public edm::EDAnalyzer {
   public:
  explicit HOAnalyzer(const edm::ParameterSet&);
  ~HOAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  int getieta(float tmpeta);
  int getiphi(float tmpphi);
  int getring(int ieta);
  float getetaval(int ieta);
  float getphival(int iphi);

  private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  

  
  // ----------member data ---------------------------
  bool isHistFill;
  bool isTrigger;
  bool isRECO;
  bool isMC;
  //bool m_digiInput;
  //std::string hoLabel_;
  //std::string hoLabel2_;
  edm::InputTag towerLabel_;


  std::vector<edm::InputTag>  jettags_;


  //! Histograms
  TH1D  *hist_nprim;
  TH1D  *hist_prim_alltrk;
  TH1D  *hist_prim_seltrk;
  TH1D  *hist_prim_goodtrk;
  TH1D  *hist_prim_dx;
  TH1D  *hist_prim_dy;
  TH2D  *hist_prim_dxy;
  TH1D  *hist_prim_dz;
  TH1D  *hist_prim_prob;

  TH1D *hist_nJet;          
  TH1D *hist_jetEt;         
  TH1D *hist_jetPt;        
  TH2D *hist_jetEtaVsPhi;   
  TH1D *hist_jetNTracks;    
  TH2D *hist_jetNTracksVsPt;
  TH1D *hist_jetEMF;        
  TH1D *hist_jetHAF;        
  TH1D *hist_jetGenEmE;     
  TH1D *hist_jetGenHadE;    
  TH1D *hist_jetGenEMF;     
  TH1D *hist_jetGenHAF;     
  TH1D *hist_jetEoverGenE;  
  TH2D *hist_jetPtoverGenPt;  
  TH1D *hist_jetCorr;       
  TH1D *hist_n90Hits;       
  TH1D *hist_fHPD;          
  TH1D *hist_nConstituents; 
  TH1D *hist_jetCHF;        
  TH1D *hist_deltar;

  TH1D *hist_pf_nJet;            
  TH1D *hist_pf_jetEt;          
  TH1D *hist_pf_jetPt;          
  TH2D *hist_pf_jetEtaVsPhi;    
  TH1D *hist_pf_jetNTracks;     
  TH2D *hist_pf_jetNTracksVsPt; 
  TH1D *hist_pf_jetCHF;         
  TH1D *hist_pf_jetNHF;         
  TH1D *hist_pf_jetCEF;         
  TH1D *hist_pf_jetNEF;         
  TH1D *hist_pf_jetGenEmE;      
  TH1D *hist_pf_jetGenHadE;     
  TH1D *hist_pf_jetGenEMF;      
  TH1D *hist_pf_jetGenHAF;     
  TH1D *hist_pf_jetEoverGenE;   
  TH2D *hist_pf_jetPtoverGenPt;   
  TH1D *hist_pf_jetCorr;        
  TH1D *hist_pf_nConstituents;  
  TH1D *hist_pf_deltar;

  //! Tree
  TTree* T1;
  unsigned int Nevents; // # of analyzed events 
  int irun, ilumi, ievt;
  float /*wtfact,*/ qscale;
  //float bestVtxX, bestVtxY, bestVtxZ;

  int ncalojets;
  //int calojet_cons[njetmx];
  float calojet_uce [njetmx], calojet_e   [njetmx], calojet_et [njetmx], calojet_pt[njetmx], calojet_ucpt[njetmx];
  float calojet_eta [njetmx], calojet_phi [njetmx];
  int   calojet_ieta[njetmx], calojet_iphi[njetmx];
  float calojet_emf [njetmx], calojet_hadf[njetmx];
  int   calojet_n90 [njetmx];
  float calojet_fHPD[njetmx], calojet_fRBX[njetmx];
  int   calojet_n90hits[njetmx];
  float calojet_towa   [njetmx];
  float calojet_hadenhb[njetmx]; 
  //float calojet_hadenhe[njetmx];
  //float calojet_hadenhf[njetmx];
  float calojet_hadenho[njetmx];
  
  float calotowE    [njetmx][500];
  float calotowHadE [njetmx][500];
  float calotowEmE  [njetmx][500];
  float calotowHoE  [njetmx][500];
  float calotowEt   [njetmx][500];
  float calotowHadEt[njetmx][500];
  float calotowEmEt [njetmx][500];
  float calotowHoEt [njetmx][500];
  float calotoweta  [njetmx][500];
  float calotowphi  [njetmx][500];
  int   calotowieta [njetmx][500];
  int   calotowiphi [njetmx][500];

  int ncalotow      [njetmx];
  float sumCaloE    [njetmx];
  float sumCaloEmE  [njetmx];
  float sumCaloHadE [njetmx];
  float sumCaloHoE  [njetmx];
  float sumCaloEt   [njetmx];
  float sumCaloEmEt [njetmx];
  float sumCaloHadEt[njetmx];
  float sumCaloHoEt [njetmx];

  int calojet_ntrks           [njetmx];
  int calojet_trkqual         [njetmx][500];
  float  calojet_trkpt        [njetmx][500];
  float  calojet_trketa       [njetmx][500];
  float  calojet_trkphi       [njetmx][500];
  float  calojet_trkPtError   [njetmx][500];
  int calojet_trkCharge       [njetmx][500];
  int calojet_trkNHit         [njetmx][500];
  float calojet_trkDxy        [njetmx][500];
  float calojet_trkDxyError   [njetmx][500];
  float calojet_trkDz         [njetmx][500];
  float calojet_trkDzError    [njetmx][500];
  float calojet_trkChi2       [njetmx][500];
  int   calojet_trkNdof       [njetmx][500];
  float calojet_trkVx         [njetmx][500];
  float calojet_trkVy         [njetmx][500];
  float calojet_trkVz         [njetmx][500];
  float calojet_trkDxyBS      [njetmx][500];
  float calojet_trkDxyErrorBS [njetmx][500];

  int npfjets;
  int pfjet_mult [njetmx];
  float pfjet_uce [njetmx];
  float pfjet_e   [njetmx];
  float pfjet_et  [njetmx];
  float pfjet_pt  [njetmx];
  float pfjet_ucpt[njetmx];
  int   pfjet_cons[njetmx];
  float pfjet_eta [njetmx];
  float pfjet_phi [njetmx];
  int   pfjet_ieta[njetmx];
  int   pfjet_iphi[njetmx];
	 
  float pfjet_chhaden[njetmx];
  float pfjet_nehaden[njetmx];
  float pfjet_phen   [njetmx];
  float pfjet_elen   [njetmx];
  float pfjet_muen   [njetmx];
  //float pfjet_hfhaden[njetmx];
  //float pfjet_hfemen [njetmx];

  //float pfjet_chemf  [njetmx];
  //float pfjet_chmuenf[njetmx];
  float pfjet_chenf[njetmx];
  float pfjet_nuenf[njetmx];
  float pfjet_phenf[njetmx];
  float pfjet_elenf[njetmx];
  float pfjet_muenf[njetmx];
  //float pfjet_hfhadenf[njetmx];
  //float pfjet_hfemenf [njetmx];

  float pfjet_tracksumecal[njetmx];
  float pfjet_tracksumhcal[njetmx];
  float pfjet_tracksumho  [njetmx];

  int pfjet_candId   [njetmx][500];
  float pfjet_cande  [njetmx][500];
  float pfjet_candet [njetmx][500];
  float pfjet_candpt [njetmx][500];
  float pfjet_candeta[njetmx][500];
  float pfjet_candphi[njetmx][500];

  int nhb;
  float hb_en [hcmx];
  float hb_et [hcmx];
  float hb_ti [hcmx];
  int hb_ieta [hcmx];
  int hb_iphi [hcmx];
  int hb_depth[hcmx];
  int hb_ring [hcmx];
  float hb_eta[hcmx];
  float hb_phi[hcmx];


  int nho;
  float ho_en [hcmx];
  float ho_et [hcmx];
  float ho_ti [hcmx];
  int ho_ieta [hcmx];
  int ho_iphi [hcmx];
  int ho_depth[hcmx];
  int ho_ring [hcmx];
  float ho_eta[hcmx];
  float ho_phi[hcmx];

  int   nref;
  float refpt  [njetmx];
  float refen  [njetmx];
  float refet  [njetmx];
  float refeta [njetmx];
  float refphi [njetmx];
  float refemen [njetmx];
  float refemet [njetmx];
  float refhaden[njetmx];
  float refhadet[njetmx];

  int nvtx;
  float vertx[vtmax];
  float verty[vtmax];
  float vertz[vtmax];
  float vertchisqp[vtmax];
  int   ntkpm[vtmax];

};




//
// constructors and destructor
//
HOAnalyzer::HOAnalyzer(const edm::ParameterSet& iConfig)
{

  isHistFill = iConfig.getUntrackedParameter<bool>("HistFill", true);
  isTrigger  = iConfig.getUntrackedParameter<bool>("Trigger", false);
  isRECO     = iConfig.getUntrackedParameter<bool>("RECO", false);
  isMC       = iConfig.getUntrackedParameter<bool>("MonteCarlo", true);

  //hoLabel_    = iConfig.getUntrackedParameter<string>("hoInput");
  //hoLabel2_   = iConfig.getUntrackedParameter<string>("hoInput2");
  towerLabel_ = iConfig.getParameter<edm::InputTag>("towerInput");


  //theRootFileName = iConfig.getUntrackedParameter<string>("RootFileName");
  
  jettags_    = iConfig.getParameter<std::vector<edm::InputTag> >( "jetTags" );

  //now do what ever initialization is needed
  edm::Service<TFileService> fs;

  //theFile = new TFile(theRootFileName.c_str(), "RECREATE");
  //theFile->cd();

  //T1 = new TTree("T1", "UseOfHO");
  T1 = fs->make<TTree>("T1","UseOfHO");

  T1->Branch("nvtx", &nvtx, "nvtx/I");
  T1->Branch("vertx", vertx, "vertx[nvtx]/F");
  T1->Branch("verty", verty, "verty[nvtx]/F");
  T1->Branch("vertz", vertz, "vertz[nvtx]/F");
  T1->Branch("vertchisqp", vertchisqp, "vertchisqp[nvtx]/F");
  T1->Branch("ntkpm", ntkpm, "ntkpm[nvtx]/I");

  T1->Branch("ncalojets", &ncalojets, "ncalojets/I");
  //T1->Branch("calojet_cons",calojet_cons,"calojet_cons[ncalojets]/I");
  T1->Branch("calojet_n90",calojet_n90,"calojet_n90[ncalojets]/I");
  T1->Branch("calojet_n90hits",calojet_n90hits,"calojet_n90hits[ncalojets]/I");
  T1->Branch("calojet_fHPD",calojet_fHPD,"calojet_fHPD[ncalojets]/F");
  T1->Branch("calojet_fRBX",calojet_fRBX,"calojet_fRBX[ncalojets]/F");

  T1->Branch("calojet_phi",calojet_phi,"calojet_phi[ncalojets]/F");
  T1->Branch("calojet_eta",calojet_eta,"calojet_eta[ncalojets]/F");
  T1->Branch("calojet_iphi",calojet_iphi,"calojet_iphi[ncalojets]/I");
  T1->Branch("calojet_ieta",calojet_ieta,"calojet_ieta[ncalojets]/I");
  T1->Branch("calojet_uce",calojet_uce,"calojet_uce[ncalojets]/F");
  T1->Branch("calojet_e",calojet_e,"calojet_e[ncalojets]/F");
  T1->Branch("calojet_et",calojet_et,"calojet_et[ncalojets]/F");
  T1->Branch("calojet_ucpt",calojet_ucpt,"calojet_ucpt[ncalojets]/F");
  T1->Branch("calojet_pt",calojet_pt,"calojet_pt[ncalojets]/F");
  T1->Branch("calojet_hadenhb",calojet_hadenhb,"calojet_hadenhb[ncalojets]/F");	 
  //T1->Branch("calojet_hadenhe",calojet_hadenhe,"calojet_hadenhe[ncalojets]/F");	 
  //T1->Branch("calojet_hadenhf",calojet_hadenhf,"calojet_hadenhf[ncalojets]/F");	   
  T1->Branch("calojet_hadenho",calojet_hadenho,"calojet_hadenho[ncalojets]/F");	 

  T1->Branch("calojet_towa",calojet_towa,"calojet_towa[ncalojets]/F");	 
  T1->Branch("calojet_hadf",calojet_hadf,"calojet_hadf[ncalojets]/F");	
  T1->Branch("calojet_emf",calojet_emf,"calojet_emf[ncalojets]/F");

  
  T1->Branch("ncalotow",ncalotow,"ncalotow[ncalojets]/I");	     
  T1->Branch("calotowE",calotowE,"calotowE[ncalojets][500]/F");	   
  T1->Branch("calotowHadE",calotowHadE,"calotowHadE[ncalojets][500]/F");	     
  T1->Branch("calotowEmE",calotowEmE,"calotowEmE[ncalojets][500]/F");	     
  T1->Branch("calotowHoE",calotowHoE,"calotowHoE[ncalojets][500]/F");	     
  T1->Branch("calotowEt",calotowEt,"calotowE[ncalojets][500]/F");	   
  T1->Branch("calotowHadEt",calotowHadEt,"calotowHadEt[ncalojets][500]/F");	     
  T1->Branch("calotowEmEt",calotowEmEt,"calotowEmEt[ncalojets][500]/F");	     
  T1->Branch("calotowHoEt",calotowHoEt,"calotowHoEt[ncalojets][500]/F");	     
  T1->Branch("calotowieta",calotowieta,"calotowieta[ncalojets][500]/I");	     
  T1->Branch("calotowiphi",calotowiphi,"calotowiphi[ncalojets][500]/I");	     
  T1->Branch("calotoweta",calotoweta,"calotoweta[ncalojets][500]/F");	     
  T1->Branch("calotowphi",calotowphi,"calotowphi[ncalojets][500]/F");	     


  T1->Branch("sumCaloE",sumCaloE,"sumCaloE[ncalojets]/F");	     
  T1->Branch("sumCaloEmE",sumCaloEmE,"sumCaloEmE[ncalojets]/F");	     
  T1->Branch("sumCaloHadE",sumCaloHadE,"sumCaloHadE[ncalojets]/F");	     
  T1->Branch("sumCaloHoE",sumCaloHoE,"sumCaloHoE[ncalojets]/F");	     


  T1->Branch("sumCaloEt",sumCaloEt,"sumCaloEt[ncalojets]/F");	     
  T1->Branch("sumCaloEmEt",sumCaloEmEt,"sumCaloEmEt[ncalojets]/F");	     
  T1->Branch("sumCaloHadEt",sumCaloHadEt,"sumCaloHadEt[ncalojets]/F");	     
  T1->Branch("sumCaloHoEt",sumCaloHoEt,"sumCaloHoEt[ncalojets]/F");	     


  T1->Branch("calojet_ntrks"          ,calojet_ntrks        ,"calojet_ntrks         [ncalojets]/I");  
  T1->Branch("calojet_trkqual"        ,calojet_trkqual      ,"calojet_trkqual       [ncalojets][500]/I");
  T1->Branch("calojet_trkpt"          ,calojet_trkpt        ,"calojet_trkpt         [ncalojets][500]/F");
  T1->Branch("calojet_trketa"         ,calojet_trketa       ,"calojet_trketa        [ncalojets][500]/F");
  T1->Branch("calojet_trkphi"         ,calojet_trkphi       ,"calojet_trkphi        [ncalojets][500]/F");
  T1->Branch("calojet_trkPtError"     ,calojet_trkPtError   ,"calojet_trkPtError    [ncalojets][500]/F");
  T1->Branch("calojet_trkCharge"      ,calojet_trkCharge    ,"calojet_trkCharge     [ncalojets][500]/I");
  T1->Branch("calojet_trkNHit"        ,calojet_trkNHit      ,"calojet_trkNHit       [ncalojets][500]/I");
  T1->Branch("calojet_trkDxy"         ,calojet_trkDxy       ,"calojet_trkDxy        [ncalojets][500]/F");
  T1->Branch("calojet_trkDxyError"    ,calojet_trkDxyError  ,"calojet_trkDxyError   [ncalojets][500]/F");
  T1->Branch("calojet_trkDz"          ,calojet_trkDz        ,"calojet_trkDz         [ncalojets][500]/F");
  T1->Branch("calojet_trkDzError"     ,calojet_trkDzError   ,"calojet_trkDzError    [ncalojets][500]/F");
  T1->Branch("calojet_trkChi2"        ,calojet_trkChi2      ,"calojet_trkChi2       [ncalojets][500]/F");
  T1->Branch("calojet_trkNdof"        ,calojet_trkNdof      ,"calojet_trkNdof       [ncalojets][500]/I");
  T1->Branch("calojet_trkVx"          ,calojet_trkVx        ,"calojet_trkVx         [ncalojets][500]/F");
  T1->Branch("calojet_trkVy"          ,calojet_trkVy        ,"calojet_trkVy         [ncalojets][500]/F");
  T1->Branch("calojet_trkVz"          ,calojet_trkVz        ,"calojet_trkVz         [ncalojets][500]/F");
  T1->Branch("calojet_trkDxyBS"       ,calojet_trkDxyBS     ,"calojet_trkDxyBS      [ncalojets][500]/F");
  T1->Branch("calojet_trkDxyErrorBS"  ,calojet_trkDxyErrorBS,"calojet_trkDxyErrorBS [ncalojets][500]/F");


  T1->Branch("npfjets", &npfjets, "npfjets/I");
  T1->Branch("pfjet_cons",pfjet_cons,"pfjet_cons[npfjets]/I");
  T1->Branch("pfjet_mult",pfjet_mult,"pfjet_mult[npfjets]/I");
  T1->Branch("pfjet_uce",pfjet_uce,"pfjet_uce[npfjets]/F");
  T1->Branch("pfjet_e",pfjet_e,"pfjet_e[npfjets]/F");
  T1->Branch("pfjet_et",pfjet_et,"pfjet_et[npfjets]/F");
  T1->Branch("pfjet_ucpt",pfjet_ucpt,"pfjet_ucpt[npfjets]/F");
  T1->Branch("pfjet_pt",pfjet_pt,"pfjet_pt[npfjets]/F");
  T1->Branch("pfjet_phi",pfjet_phi,"pfjet_phi[npfjets]/F");
  T1->Branch("pfjet_eta",pfjet_eta,"pfjet_eta[npfjets]/F");
  T1->Branch("pfjet_iphi",pfjet_iphi,"pfjet_iphi[npfjets]/I");
  T1->Branch("pfjet_ieta",pfjet_ieta,"pfjet_ieta[npfjets]/I");

  T1->Branch("pfjet_chhaden",pfjet_chhaden,"pfjet_chhaden[npfjets]/F");	 
  T1->Branch("pfjet_nehaden",pfjet_nehaden,"pfjet_nehaden[npfjets]/F");	 
  T1->Branch("pfjet_phen",pfjet_phen,"pfjet_phen[npfjets]/F");	 
  T1->Branch("pfjet_elen",pfjet_elen,"pfjet_elen[npfjets]/F");	 
  T1->Branch("pfjet_muen",pfjet_muen,"pfjet_muen[npfjets]/F");	 
  //T1->Branch("pfjet_hfhaden",pfjet_hfhaden,"pfjet_hfhaden[npfjets]/F");	 
  //T1->Branch("pfjet_hfemen",pfjet_hfemen,"pfjet_hfemen[npfjets]/F");	 
  //T1->Branch("pfjet_chemf",pfjet_chemf,"pfjet_chemf[npfjets]/F");	 
  //T1->Branch("pfjet_chmuenf",pfjet_chmuenf,"pfjet_chmuef[npfjets]/F");	 

  T1->Branch("pfjet_chenf",pfjet_chenf,"pfjet_chenf[npfjets]/F");	   
  T1->Branch("pfjet_nuenf",pfjet_nuenf,"pfjet_nuenf[npfjets]/F");	   
  T1->Branch("pfjet_phenf",pfjet_phenf,"pfjet_phenf[npfjets]/F");	   
  T1->Branch("pfjet_elenf",pfjet_elenf,"pfjet_elenf[npfjets]/F");	   
  T1->Branch("pfjet_muenf",pfjet_muenf,"pfjet_muenf[npfjets]/F");	   
  //T1->Branch("pfjet_hfhadenf",pfjet_hfhadenf,"pfjet_hfhadenf[npfjets]/F");	   
  //T1->Branch("pfjet_hfemenf",pfjet_hfemenf,"pfjet_hfemenf[npfjets]/F");	   

  T1->Branch("pfjet_tracksumecal",pfjet_tracksumecal,"pfjet_tracksumecal[npfjets]/F");
  T1->Branch("pfjet_tracksumhcal",pfjet_tracksumhcal,"pfjet_tracksumhcal[npfjets]/F");
  T1->Branch("pfjet_tracksumho",pfjet_tracksumho,"pfjet_tracksumho[npfjets]/F");

  T1->Branch("pfjet_candId" ,    pfjet_candId ,"pfjet_candId [npfjets][500]/I");
  T1->Branch("pfjet_cande"  ,    pfjet_cande  ,"pfjet_cande  [npfjets][500]/F");
  T1->Branch("pfjet_candet" ,    pfjet_candet ,"pfjet_candet [npfjets][500]/F");
  T1->Branch("pfjet_candpt" ,    pfjet_candpt ,"pfjet_candpt [npfjets][500]/F");
  T1->Branch("pfjet_candeta",    pfjet_candeta,"pfjet_candeta[npfjets][500]/F");
  T1->Branch("pfjet_candphi",    pfjet_candphi,"pfjet_candphi[npfjets][500]/F");

  T1->Branch("nhb",&nhb,"nhb/I");	   
  T1->Branch("hb_en"   ,hb_en    ,"hb_en[nhb]/F");	   
  T1->Branch("hb_et"   ,hb_et    ,"hb_et[nhb]/F");	   
  T1->Branch("hb_ti"   ,hb_ti    ,"hb_ti[nhb]/F");	   
  T1->Branch("hb_ieta" ,hb_ieta  ,"hb_ieta[nhb]/I");	   
  T1->Branch("hb_iphi" ,hb_iphi  ,"hb_iphi[nhb]/I");	   
  T1->Branch("hb_depth",hb_depth ,"hb_depth[nhb]/I");	   
  T1->Branch("hb_ring" ,hb_ring  ,"hb_ring[nhb]/I");	   
  T1->Branch("hb_eta"  ,hb_eta   ,"hb_eta[nhb]/F");	   
  T1->Branch("hb_phi"  ,hb_phi   ,"hb_phi[nhb]/F");	   


  T1->Branch("nho",&nho,"nho/I");	   
  T1->Branch("ho_en",ho_en,"ho_en[nho]/F");	   
  T1->Branch("ho_et",ho_et,"ho_et[nho]/F");	   
  T1->Branch("ho_ti",ho_ti,"ho_ti[nho]/F");	   
  T1->Branch("ho_ieta",ho_ieta,"ho_ieta[nho]/I");	   
  T1->Branch("ho_iphi",ho_iphi,"ho_iphi[nho]/I");	   
  T1->Branch("ho_depth",ho_depth,"ho_depth[nho]/I");	   
  T1->Branch("ho_ring",ho_ring,"ho_ring[nho]/I");	   
  T1->Branch("ho_eta",ho_eta,"ho_eta[nho]/F");	   
  T1->Branch("ho_phi",ho_phi,"ho_phi[nho]/F");	   
  
  if(isMC){
    T1->Branch("qscale", &qscale, "qscale/F");
    //T1->Branch("wtfact", &wtfact, "wtfact/F");
    T1->Branch("nref",  &nref, "nref/I");
    T1->Branch("refpt",  refpt, "refpt[nref]/F");
    T1->Branch("refen",  refen, "refen[nref]/F");
    T1->Branch("refet",  refet, "refet[nref]/F");
    T1->Branch("refeta", refeta, "refeta[nref]/F");
    T1->Branch("refphi", refphi, "refphi[nref]/F");
    T1->Branch("refemen", refemen, "refemen[nref]/F");
    T1->Branch("refhaden", refhaden, "refhaden[nref]/F");
    T1->Branch("refemet", refemet, "refemet[nref]/F");
    T1->Branch("refhadet", refhadet, "refhadet[nref]/F");
  }




  if (isHistFill) {
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    TFileDirectory histDir = fs->mkdir("Histos");

    hist_nprim            = histDir.make<TH1D>("hist_nprim"                    , "# of primary vtx", 50, -0.5, 49.5);
    hist_prim_alltrk      = histDir.make<TH1D>("hist_primall"                  , "All tracks in primary vtx", 240, -0.5, 239.5);
    hist_prim_goodtrk     = histDir.make<TH1D>("hist_primgood"                 , "Good tracks in primary vtx", 240, -0.5, 239.5);
    hist_prim_seltrk      = histDir.make<TH1D>("hist_primsel"                  , "Selected tracks in primary vtx", 240, -0.5, 239.5);
    hist_prim_dx          = histDir.make<TH1D>("hist_primdx"                   , "#Delta x of prim wrt beam spot", 120, -2.4, 2.4);
    hist_prim_dy          = histDir.make<TH1D>("hist_primdy"                   , "#Delta y of prim wrt beam spot", 120, -2.4, 2.4);
    hist_prim_dxy         = histDir.make<TH2D>("hist_primdxy"                  , "#Delta y vs #Delta x of prim", 60, -0.18, 0.18, 60, -0.18, 0.18);
    hist_prim_dz          = histDir.make<TH1D>("hist_primdz"                   , "#Delta z of prim wrt beam spot", 120, -30.0, 30.0);
    hist_prim_prob        = histDir.make<TH1D>("hist_primprob"                 , "vertex fit prob", 400, 0.0, 1.0);

    
    hist_nJet             = histDir.make<TH1D>( "hist_nJet"                    , "Number of Calo Jets", 10, 0, 10 ) ;
    hist_jetEt            = histDir.make<TH1D>( "hist_jetEt"                   , "log10(1+Et_{calojets})",240,1.0, 7.0);
    hist_jetPt            = histDir.make<TH1D>( "hist_jetPt"                   , "log10(1+pt {calojets})", 600, 1.0, 7.0) ;
    hist_jetEtaVsPhi      = histDir.make<TH2D>( "hist_jetEtaVsPhi"             , "Jet #phi versus #eta;#eta;#phi", 50, -5.0, 5.0, 50, -M_PI, M_PI ) ;
    hist_jetNTracks       = histDir.make<TH1D>( "hist_jetNTracks"              , "Jet N_{TRACKS}", 20, 0, 20 ) ;
    hist_jetNTracksVsPt   = histDir.make<TH2D>( "hist_jetNTracksVsPt"          , "Number of Tracks versus Jet p_{T};Jet p_{T}(GeV/c) ;N_{Tracks}",600, 1.0, 7.0, 20, 0, 20 ) ;
    hist_jetEMF           = histDir.make<TH1D>( "hist_jetEMF"                  , "Jet EMF", 200, 0, 1) ;
    hist_jetHAF           = histDir.make<TH1D>( "hist_jetHAF"                  , "Jet HADF", 200, 0, 1) ;
    hist_jetGenEmE        = histDir.make<TH1D>( "hist_jetGenEmE"               , "Gen Jet EM Energy", 200, 0, 200 ) ;
    hist_jetGenHadE       = histDir.make<TH1D>( "hist_jetGenHadE"              , "Gen Jet HAD Energy", 200, 0, 200 ) ;
    hist_jetGenEMF        = histDir.make<TH1D>( "hist_jetGenEMF"               , "Gen Jet EMF", 200, 0, 1) ;
    hist_jetGenHAF        = histDir.make<TH1D>( "hist_jetGenHAF"               , "Gen Jet HADF", 200, 0, 1) ;
    hist_jetEoverGenE     = histDir.make<TH1D>( "hist_jetEoverGenE"            , "Energy of reco Jet / Energy of gen Jet", 200, 0, 2.0) ;
    hist_jetPtoverGenPt   = histDir.make<TH2D>( "hist_jetPtoverGenPt"          , "Pt of reco Jet / Pt of gen Jet", 600, 1.0, 7.0, 200, 0, 2.0) ;
    hist_jetCorr          = histDir.make<TH1D>( "hist_jetCorr"                 , "Jet Correction Factor", 200, 0, 1.0 ) ;
    hist_n90Hits          = histDir.make<TH1D>( "hist_n90Hits"                 , "Jet n90Hits", 20, 0, 20) ;
    hist_fHPD             = histDir.make<TH1D>( "hist_fHPD"                    , "Jet fHPD", 200, 0, 1) ;
    hist_nConstituents    = histDir.make<TH1D>( "hist_nConstituents"           , "Jet nConstituents", 20, 0, 20 ) ;
    hist_jetCHF           = histDir.make<TH1D>( "hist_jetCHF"                  , "Jet Charged Hadron Fraction", 200, 0, 1.0) ;
    hist_deltar           = histDir.make<TH1D>( "hist_deltar"                  , "Gen-Reco DeltaR", 200, 0, 1.0) ;

    hist_pf_nJet             = histDir.make<TH1D>( "hist_pf_nJet"              , "Number of PF Jets", 10, 0, 10 ) ;
    hist_pf_jetEt            = histDir.make<TH1D>( "hist_jetEt"                , "log10(1+Et_{pfjets})",240,1.0, 7.0);
    hist_pf_jetPt            = histDir.make<TH1D>( "hist_pf_jetPt"             , "log10(1+p_{T} {pfjets})", 600, 1.0, 7.0) ;
    hist_pf_jetEtaVsPhi      = histDir.make<TH2D>( "hist_pf_jetEtaVsPhi"       , "PFJet #phi versus #eta;#eta;#phi", 50, -5.0, 5.0, 50, -M_PI, M_PI ) ;
    hist_pf_jetNTracks       = histDir.make<TH1D>( "hist_pf_jetNTracks"        , "PFJet N_{TRACKS}", 20, 0, 20 ) ;
    hist_pf_jetNTracksVsPt   = histDir.make<TH2D>( "hist_pf_jetNTracksVsPt"    , "Number of Tracks versus Jet p_{T};Jet p_{T}(GeV/c) ;N_{Tracks}",600, 1.0, 7.0, 20, 0, 20 ) ;
    hist_pf_jetCHF           = histDir.make<TH1D>( "hist_pf_jetCHF"            , "PFJet CHF", 200, 0, 1) ;
    hist_pf_jetNHF           = histDir.make<TH1D>( "hist_pf_jetNHF"            , "PFJet NHF", 200, 0, 1) ;
    hist_pf_jetCEF           = histDir.make<TH1D>( "hist_pf_jetCEF"            , "PFJet CEF", 200, 0, 1) ;
    hist_pf_jetNEF           = histDir.make<TH1D>( "hist_pf_jetNEF"            , "PFJet NEF", 200, 0, 1) ;
    hist_pf_jetGenEmE        = histDir.make<TH1D>( "hist_pf_jetGenEmE"         , "Gen Jet EM Energy", 200, 0, 200 ) ;
    hist_pf_jetGenHadE       = histDir.make<TH1D>( "hist_pf_jetGenHadE"        , "Gen Jet HAD Energy", 200, 0, 200 ) ;
    hist_pf_jetGenEMF        = histDir.make<TH1D>( "hist_pf_jetGenEMF"         , "Gen Jet EMF", 200, 0, 1) ;
    hist_pf_jetGenHAF        = histDir.make<TH1D>( "hist_pf_jetGenHAF"         , "Gen Jet HADF",200, 0, 1) ;
    hist_pf_jetEoverGenE     = histDir.make<TH1D>( "hist_pf_jetEoverGenE"      , "Energy of reco Jet / Energy of gen Jet", 200, 0, 2.0) ;
    hist_pf_jetPtoverGenPt   = histDir.make<TH2D>( "hist_pf_jetPtoverGenPt"    , "Pt of reco Jet / Pt of gen Jet", 600, 1.0, 7.0, 200, 0, 2.0) ;
    hist_pf_jetCorr          = histDir.make<TH1D>( "hist_pf_jetCorr"           , "PFJet Correction Factor", 200, 0, 1.0 ) ;
    hist_pf_nConstituents    = histDir.make<TH1D>( "hist_pf_nConstituents"     , "PFJet nConstituents", 20, 0, 20 ) ;
    hist_pf_deltar           = histDir.make<TH1D>( "hist_pf_deltar"            , "PFJet Gen-Reco DeltaR", 200, 0, 1.0) ;
  }

  Nevents = 0;
}


HOAnalyzer::~HOAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
int HOAnalyzer::getieta(float tmpeta) {
  int ietabin = -999;
  for(int ix=0; ix<netahotower; ix++){
    if(fabs(tmpeta) >= etahotowers[ix] && fabs(tmpeta) < etahotowers[ix+1]){ietabin = ix+1; break;}
  }
  return (tmpeta >= 0) ? ietabin : -ietabin;
}  
int HOAnalyzer::getiphi(float tmpphi){
  if (tmpphi<0) tmpphi +=2*M_PI;
  int iphibn = int(tmpphi*nphimx/(2.*M_PI)) + 1;
  if (iphibn>nphimx) iphibn -=nphimx;
  if (iphibn<=0) iphibn +=nphimx;
  return iphibn;
}
int HOAnalyzer::getring(int ieta) {
  int iring=-999;
  if( abs(ieta) <= 5) iring=0;
  else if( abs(ieta) <= 10) iring=1;
  else if( abs(ieta) <= 15) iring=2;
  return (ieta > 0) ? iring : -iring;
}
float HOAnalyzer::getetaval(int ieta) {
  if (ieta==0 || abs(ieta)>netahotower) return -999;
  if (ieta>0) {
    return 0.5*(etahotowers[ieta-1]+etahotowers[ieta]);
  } else {
    return -0.5*(etahotowers[abs(ieta)-1]+etahotowers[abs(ieta)]);
  }
}
float HOAnalyzer::getphival(int iphi) {
  float phi = (iphi - 0.5)*(2*M_PI)/nphimx;
  if (phi > M_PI) phi -=2*M_PI;
  if (phi <-M_PI) phi +=2*M_PI;
  return phi;
}



// ------------ method called for each event  ------------
void
HOAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   Nevents++;

   irun  = iEvent.id().run();
   ilumi = iEvent.luminosityBlock();
   ievt  = iEvent.id().event();

   float weight = 1.0;
   if (!iEvent.isRealData()) {
     edm::Handle< GenEventInfoProduct > geneventinfo;
     iEvent.getByLabel("generator", geneventinfo);
     qscale = -1.;
     //wtfact = -1.;
     if (geneventinfo.isValid()) {
       //wtfact = geneventinfo->weight();
       qscale = geneventinfo->qScale();
     }
   }

   reco::TrackBase::Point beamPoint(0,0, 0);
   edm::Handle<reco::BeamSpot> beamSpotH;
   iEvent.getByLabel("offlineBeamSpot",beamSpotH);
   if (beamSpotH.isValid()){
     beamPoint = beamSpotH->position();
   }
   
   
   Handle<VertexCollection> recVtxs;
   iEvent.getByLabel("offlinePrimaryVertices", recVtxs);
   //! Not the best choice but for time been
   //const reco::Vertex* v = &(*recVtxs->begin()); 
   //const reco::Vertex *bestVertex=0; 
   //float chi2ndf = 99999;
   nvtx=0;
   if (recVtxs.isValid()) {
     for (reco::VertexCollection::const_iterator vert=recVtxs->begin(); vert<recVtxs->end(); vert++) {
       if (vert->isValid() && !vert->isFake() && vert->ndof()>=4) {
	 int ngoodtrk=0;
	 int nseltrk = 0;
	 float chisqp = ChiSquaredProbability(vert->chi2(),vert->ndof());
	 // if(chisqp < chi2ndf){
	 //   chi2ndf = chisqp;
	 //   bestVertex = &(*vert);
	 // }
	 //if(chisqp > 0.01){
	 // bestVertex = &(*vert);	   
	 //}
	 for (Vertex::trackRef_iterator reftrk =vert->tracks_begin(); reftrk<vert->tracks_end(); reftrk++) {
	   if ((*reftrk)->quality(TrackBase::highPurity) && vert->trackWeight(*reftrk)>0) {
	     ngoodtrk++;
	     if ((*reftrk)->normalizedChi2()<100000 &&
		 abs((*reftrk)->dxy()) < 10000 &&
		 (*reftrk)->pt() >0.50) {nseltrk++; }
	   }
	 }
	 vertx[nvtx] = vert->position().x();
	 verty[nvtx] = vert->position().y();
	 vertz[nvtx] = vert->position().z();
	 vertchisqp[nvtx] = chisqp ;
	 ntkpm[nvtx] = 1000*(1000*min(int(vert->tracksSize()),999) + min(ngoodtrk,999)) + min(999, nseltrk);

	 if(isHistFill){
	   hist_prim_alltrk->Fill(vert->tracksSize());
	   hist_prim_goodtrk->Fill(ngoodtrk);
	   hist_prim_seltrk->Fill(nseltrk);
	   hist_prim_dx->Fill(vert->position().x() - beamPoint.x());
	   hist_prim_dy->Fill(vert->position().y() - beamPoint.y());
	   hist_prim_dxy->Fill(vert->position().x() - beamPoint.x(), vert->position().y() - beamPoint.y());
	   hist_prim_dz->Fill(vert->position().z() - beamPoint.z());
	   hist_prim_prob->Fill(chisqp);
	   //hist_prim_prob->Fill(max(float(-30.0), log10(chisqp)));
	 }
	 //std::cout << Nevents << " vz : " << vert->position().z() << "\t  chisqp : "<< chisqp << "\t # of tracks : " << vert->tracksSize() << " " << ngoodtrk << " " << nseltrk << std::endl;
	 nvtx++;
       }
       if(isHistFill)hist_nprim->Fill(nvtx);       
     }
   }
   //bestVtxX = bestVertex->position().x();
   //bestVtxY = bestVertex->position().y();
   //bestVtxZ = bestVertex->position().z();
   //std::cout << " Total  # of vertices  : "  << recVtxs->size() << "\t best vz : " <<  bestVtxZ << std::endl;

   //edm::Handle<vector<reco::GenJet> >genjets;
   //iEvent.getByLabel(genjet_, genjets);


   //! Track collection
   //edm::Handle<reco::TrackCollection> tracks;
   //iEvent.getByLabel("generalTracks", tracks);


   if(jettags_.empty()){
     std::cout << " No jet algorithms are defined  please check your cfg's ..... " << std::endl;
   }

   ncalojets=0;
   npfjets=0;
   nref=0;
   int il=0;
   for (std::vector<edm::InputTag>::const_iterator jetLabel = jettags_.begin(); jetLabel != jettags_.end(); ++jetLabel, il++) {
     
     std::string jLabel = jetLabel->label();
     edm::Handle<pat::JetCollection> jets;
     iEvent.getByLabel(jLabel, jets);
     
     // if(il==0){
     //   std::cout<<" iEvent : " << Nevents 
     //  		<< "\t patjets  : " << jLabel.c_str() << "\t" << jets->size() 
     //  		<< std::endl;
     // }


     for(unsigned int ij=0; ij< jets->size(); ++ij){
       const pat::Jet &jet = (*jets)[ij];
       
       if(fabs(jet.eta()) > 2.6 || jet.pt() <= 20.0)continue;
       
       //float rawpt  = jet.correctedJet("Uncorrected").pt();
       //float rawe   = jet.correctedJet("Uncorrected").energy();
       //float corrpt = jet.pt();

       const reco::GenJet *genjet = jet.genJet();
       const reco::TrackRefVector & jetTracks = jet.associatedTracks();

       double dr = 9999;
       if(il==0){
	 if(genjet){
	   refpt   [nref] = genjet->pt();
	   refeta  [nref] = genjet->eta();
	   refphi  [nref] = genjet->phi();
	   refen   [nref] = genjet->energy();
	   refemen [nref] = genjet->emEnergy(); 
	   refhaden[nref] = genjet->hadEnergy(); 	   
	   
	   refet  [nref]  = refen[nref]/cosh(genjet->eta());
	   refemet[nref]  = refemen[nref]/cosh(genjet->eta());
	   refhadet[nref] = refhaden[nref]/cosh(genjet->eta());

	   dr = reco::deltaR(jet.eta(), jet.phi(), genjet->eta(), genjet->phi());

	   //std::cout<<"\t refet  : " << refet[nref] << "\t ref emEt : "<< refemet[nref] << "\t ref hadEt : "<< refhadet[nref] << std::endl;
	   nref++;
	 }
       }
       //if(dr < 0.3)continue;

       //int ipass=0;
       float jetemf=-999;
       int istat=0;
       if( jet.isCaloJet() ){
	 istat=1;
	 if(jet.getCaloConstituents().size() < 2)continue;
	 //if (ncalojets > 1 && (fabs(jet.eta()) > 4.5 || jet.pt() < 30.0)) continue;
	 
	 jetemf = jet.emEnergyFraction();
	 if(jetemf < 0.01 || jetemf >=1.0)continue;
	 
	 reco::JetID jetid = jet.jetID();
	 //if (jetid.n90Hits <= 1) ipass = 0;
	 //if (jetid.fHPD  > 0.98) ipass = 0;
	 //if (!ipass)continue;


	 if(isHistFill){
	   hist_jetEt ->Fill(log10(1+jet.et()),weight);
	   hist_jetPt->Fill(log10(1+jet.pt()),weight);
	   hist_jetEtaVsPhi->Fill(jet.eta(),jet.phi(),weight);
	   hist_jetNTracks->Fill(jetTracks.size(),weight);
	   hist_jetNTracksVsPt->Fill(log10(1+jet.pt()),jetTracks.size(),weight);
	   hist_jetEMF->Fill(jet.emEnergyFraction(),weight ); 
	   hist_jetHAF->Fill(jet.energyFractionHadronic(),weight ); 
	   hist_jetCorr->Fill(jet.jecFactor("Uncorrected"));
	   hist_n90Hits->Fill(jetid.n90Hits,weight);
	   hist_fHPD->Fill(jetid.fHPD,weight);
	   hist_nConstituents->Fill(jet.getCaloConstituents().size(),weight);

	   if(genjet){
	     hist_jetGenEmE->Fill(genjet->emEnergy(),weight);
	     hist_jetGenHadE->Fill(genjet->hadEnergy(),weight);
	     hist_jetEoverGenE->Fill(jet.energy()/genjet->energy(),weight);
	     hist_jetPtoverGenPt->Fill(log10(1+jet.pt()), jet.pt()/genjet->pt(),weight);
	     hist_jetGenEMF->Fill(genjet->emEnergy()/genjet->energy(),weight);
	     hist_jetGenHAF->Fill(genjet->hadEnergy()/genjet->energy(),weight);
	     hist_deltar->Fill(dr,weight);
	   }


	   calojet_ntrks[ncalojets]=0;
	   TLorentzVector p4_tracks(0,0,0,0);	   
	   for ( reco::TrackRefVector::const_iterator itrk = jetTracks.begin(),
		   itrkEnd = jetTracks.end(); itrk != itrkEnd; ++itrk ) {
	     TLorentzVector p4_trk;
	     double M_PION = 0.140;
	     p4_trk.SetPtEtaPhiM( (*itrk)->pt(), (*itrk)->eta(), (*itrk)->phi(), M_PION );
	     p4_tracks += p4_trk;
	     
	     calojet_trkqual    [ncalojets][calojet_ntrks[ncalojets]]  = (*itrk)->quality(TrackBase::highPurity);
	     calojet_trkpt      [ncalojets][calojet_ntrks[ncalojets]]  = (*itrk)->pt();
	     calojet_trketa     [ncalojets][calojet_ntrks[ncalojets]]  = (*itrk)->eta();
	     calojet_trkphi     [ncalojets][calojet_ntrks[ncalojets]]  = (*itrk)->phi();
	     calojet_trkPtError [ncalojets][calojet_ntrks[ncalojets]]  = (*itrk)->ptError();
	     calojet_trkCharge  [ncalojets][calojet_ntrks[ncalojets]]  = (*itrk)->charge();
	     calojet_trkNHit    [ncalojets][calojet_ntrks[ncalojets]]  = (*itrk)->numberOfValidHits();
	     calojet_trkDxy     [ncalojets][calojet_ntrks[ncalojets]]  = (*itrk)->dxy();
	     calojet_trkDxyError[ncalojets][calojet_ntrks[ncalojets]]  = (*itrk)->dxyError();
	     calojet_trkDz      [ncalojets][calojet_ntrks[ncalojets]]  = (*itrk)->dz();
	     calojet_trkDzError [ncalojets][calojet_ntrks[ncalojets]]  = (*itrk)->dzError();
	     calojet_trkChi2    [ncalojets][calojet_ntrks[ncalojets]]  = (*itrk)->chi2();
	     calojet_trkNdof    [ncalojets][calojet_ntrks[ncalojets]]  = (*itrk)->ndof();
	     calojet_trkVx      [ncalojets][calojet_ntrks[ncalojets]]  = (*itrk)->vx();
	     calojet_trkVy      [ncalojets][calojet_ntrks[ncalojets]]  = (*itrk)->vy();
	     calojet_trkVz      [ncalojets][calojet_ntrks[ncalojets]]  = (*itrk)->vz();
	     calojet_trkDxyBS   [ncalojets][calojet_ntrks[ncalojets]]  = (*itrk)->dxy(beamPoint);
	     calojet_trkDxyErrorBS [ncalojets][calojet_ntrks[ncalojets]]  = sqrt(pow(calojet_trkDxyError[ncalojets][calojet_ntrks[ncalojets]],2) + beamSpotH->BeamWidthX()*beamSpotH->BeamWidthY());

	     calojet_ntrks[ncalojets]++;

	   }
	   hist_jetCHF->Fill(p4_tracks.Energy()/jet.energy(),weight);
	 }

	 calojet_uce [ncalojets] = jet.correctedJet("Uncorrected").energy(); 
	 calojet_e   [ncalojets] = jet.energy(); 
	 calojet_et  [ncalojets] = jet.et();
	 calojet_pt  [ncalojets] = jet.pt();
	 calojet_ucpt[ncalojets] = jet.correctedJet("Uncorrected").pt();
	 calojet_n90 [ncalojets] = jet.n90();
	 calojet_n90hits [ncalojets] = jetid.n90Hits;
	 calojet_fHPD[ncalojets] = jetid.fHPD;
	 calojet_fRBX[ncalojets] = jetid.fRBX;
	 //calojet_cons[ncalojets] = jet.getCaloConstituents().size(); 
	 calojet_eta [ncalojets] = jet.eta();
	 calojet_phi [ncalojets] = jet.phi();
	 
	 calojet_hadenhb[ncalojets] = jet.hadEnergyInHB(); 
	 //calojet_hadenhe[ncalojets] = jet.hadEnergyInHE();
	 //calojet_hadenhf[ncalojets] = jet.hadEnergyInHF();
	 calojet_hadenho[ncalojets] = jet.hadEnergyInHO();
	 
	 calojet_towa[ncalojets] = jet.towersArea();
	 
	 calojet_hadf[ncalojets] = jet.energyFractionHadronic();
	 calojet_emf [ncalojets] = jetemf; 
	 
	 calojet_ieta [ncalojets] = getieta(jet.eta());
	 calojet_iphi [ncalojets] = getiphi(jet.phi());
	 
	 //std::cout <<"\t \t Calo Jets  "<< "\t corrpt : "<< jet.pt() <<"\t raw pt : "<< calojet_ucpt[ncalojets] < "\t ieta : " << calojet_ieta [ncalojets] << std::endl;
	 //const reco::TrackRefVector & jetTracks = jet.associatedTracks();
	 
	 
	 //! Try to read Calo Constituents
	 sumCaloEt[ncalojets]=sumCaloEmEt[ncalojets]=sumCaloHadEt[ncalojets]=sumCaloHoEt[ncalojets]=0;
	 sumCaloE[ncalojets]=sumCaloEmE[ncalojets]=sumCaloHadE[ncalojets]=sumCaloHoE[ncalojets]=0;
	 ncalotow[ncalojets]=0;
	 for (unsigned ic = 0; ic < jet.getCaloConstituents().size(); ic++) {
	   CaloTowerPtr tower = (jet.getCaloConstituents())[ic];
	   CaloTowerDetId id = tower->id();

	   calotowE    [ncalojets][ncalotow[ncalojets]] = tower->energy(); 
	   calotowHadE [ncalojets][ncalotow[ncalojets]] = tower->hadEnergy();
	   calotowEmE  [ncalojets][ncalotow[ncalojets]] = tower->emEnergy();
	   calotowHoE  [ncalojets][ncalotow[ncalojets]] = tower->outerEnergy();
	   calotowEt   [ncalojets][ncalotow[ncalojets]] = tower->et(); 
	   calotowHadEt[ncalojets][ncalotow[ncalojets]] = tower->hadEt();
	   calotowEmEt [ncalojets][ncalotow[ncalojets]] = tower->emEt();
	   calotowHoEt [ncalojets][ncalotow[ncalojets]] = tower->outerEt();
	   calotoweta  [ncalojets][ncalotow[ncalojets]] = tower->eta();
	   calotowphi  [ncalojets][ncalotow[ncalojets]] = tower->phi();
	   calotowieta [ncalojets][ncalotow[ncalojets]] = id.ieta();
	   calotowiphi [ncalojets][ncalotow[ncalojets]] = id.iphi();

	   ncalotow[ncalojets]++;

	   //float aeta = tower->eta();
	   //float aphi = tower->phi();
	   //float theta   = tower->p4(bestVtxZ).theta();
	   //float caleta  = tower->p4(bestVtxZ).eta();
	   //float calphi  = tower->p4(bestVtxZ).phi();
	   //std::cout << "\t \t \t \t ###################  : "<< theta << "\t caleta : "<< caleta << " aeta : " << aeta << "\t calphi : "<< calphi << " aphi : "<< aphi << std::endl;

	   //int ieta_calo = id.ieta();
	   //int iphi_calo = id.iphi();
	   
	   sumCaloEt   [ncalojets] += tower->et();
	   sumCaloEmEt [ncalojets] += tower->emEt();
	   sumCaloHadEt[ncalojets] += tower->hadEt();
	   sumCaloHoEt [ncalojets] += tower->outerEt();

	   sumCaloE   [ncalojets]   += tower->energy();
	   sumCaloEmE [ncalojets]   += tower->emEnergy();
	   sumCaloHadE[ncalojets]   += tower->hadEnergy();
	   sumCaloHoE [ncalojets]   += tower->outerEnergy();


	   // std:: cout<< " \t " << ic  <<"  " <<  "\t caltowEt : " 
	   // 	       << caltowEt << "\t hadet : "<< calohadEt << "\t caloEmEt " << caloEmEt  <<" \t caloHoEt "<< caloHoEt 
	   // 	       <<"\t eta : "<< caleta << "\t calphi : "<< calphi << std::endl;
	   
	   //std::cout << "\t \t \t \t ###################  : "<< "\t caleta : "<< caleta << " ieta : " << ieta_calo << "\t calphi : "<< calphi << " iphi : "<< iphi_calo << std::endl;

	   if(ncalotow[ncalojets]>500){
	     std::cout<<"##########################################################################" <<std::endl;
	     std::cout<<"##########################################################################" <<std::endl;
	     std::cout<<"#################### ncalotow > 500 please check #########################" <<std::endl;
	     std::cout<<"##########################################################################" <<std::endl;
	     std::cout<<"##########################################################################" <<std::endl;
	   }

	 }   

	 //std::cout <<"\t \t Calo Jets  "<< "\t  jet Et : "<< calojetet  [ncalojets]  << "   " << ncal << "\t SumEt "<< sumEt <<"\t EmEt : "<< sumEmEt << "\t sumHadEt : "<< sumHadEt << "\t sumHoEt : " << sumHoEt << std::endl;

	 ncalojets++;
       }
       else if( jet.isPFJet() ){
	 
	 istat=2;
	 
	 // float jetenergy_uncorrected = jet.chargedHadronEnergy()
	 //   + jet.neutralHadronEnergy()  
	 //   + jet.photonEnergy()
	 //   + jet.electronEnergy()
	 //   + jet.muonEnergy()
	 //   + jet.HFHadronEnergy()
	 //   + jet.HFEMEnergy();
	 
	 pfjet_mult[npfjets] = pow(10,0)*min(9,jet.muonMultiplicity())
	   + pow(10,1)*min(99,jet.chargedHadronMultiplicity())
	   + pow(10,3)*min(99,jet.neutralHadronMultiplicity())
	   + pow(10,5)*min(99,jet.electronMultiplicity())
	   + pow(10,7)*min(99,jet.photonMultiplicity());
	 
	 // if(Nevents==1){
	 //   std::cout << " pfmult : " << pfjetmult[npfjets] 
	 // 	     << "\t mu   : " << jet.muonMultiplicity()
	 // 	     << "\t ch   : " << jet.chargedHadronMultiplicity()
	 // 	     << "\t ne   : " << jet.neutralHadronMultiplicity()
	 // 	     << "\t el   : " << jet.electronMultiplicity()
	 // 	     << "\t ph   : " << jet.photonMultiplicity()
	 // 	     << std::endl;
	 // }

	 if(isHistFill){
	   hist_pf_jetEt->Fill(log10(1+jet.et()),weight);
	   hist_pf_jetPt->Fill(log10(1+jet.pt()),weight );
	   hist_pf_jetEtaVsPhi->Fill(jet.eta(),jet.phi(),weight);
	   hist_pf_jetNTracks->Fill(jetTracks.size(),weight);
	   hist_pf_jetNTracksVsPt->Fill(log10(1+jet.pt()),jetTracks.size(),weight);
	   hist_pf_nConstituents->Fill(jet.getPFConstituents().size(),weight);
	   hist_pf_jetCEF->Fill(jet.chargedEmEnergyFraction(),weight);
	   hist_pf_jetNEF->Fill(jet.neutralEmEnergyFraction(),weight);
	   hist_pf_jetCHF->Fill(jet.chargedHadronEnergyFraction(),weight);
	   hist_pf_jetNHF->Fill(jet.neutralHadronEnergyFraction(),weight);
	   hist_pf_jetCorr->Fill(jet.jecFactor("Uncorrected"));

	   if (genjet) {
	     hist_pf_jetGenEmE->Fill(genjet->emEnergy(),weight);
	     hist_pf_jetGenHadE->Fill(genjet->hadEnergy(),weight);
	     hist_pf_jetEoverGenE->Fill(jet.energy() / genjet->energy(),weight);
	     hist_pf_jetPtoverGenPt->Fill(log10(1+jet.pt()), jet.pt()/genjet->pt(),weight);
	     hist_pf_jetGenEMF->Fill(genjet->emEnergy()/genjet->energy(),weight);
	     hist_pf_jetGenHAF->Fill(genjet->hadEnergy()/genjet->energy(),weight);
	     hist_pf_deltar->Fill((reco::deltaR(jet.eta(), jet.phi(), genjet->eta(), genjet->phi())),weight);
	   }
	 }
	 
	 pfjet_uce [npfjets] = jet.correctedJet("Uncorrected").energy(); 
	 pfjet_et  [npfjets] = jet.et(); 
	 pfjet_e   [npfjets] = jet.energy(); 
	 pfjet_pt  [npfjets] = jet.pt();
	 pfjet_ucpt[npfjets] = jet.correctedJet("Uncorrected").pt();
	 pfjet_cons[npfjets] = jet.getPFConstituents().size(); 
	 pfjet_eta [npfjets] = jet.eta();
	 pfjet_phi [npfjets] = jet.phi();
	 
	 pfjet_chhaden[npfjets] = jet.chargedHadronEnergy();
	 pfjet_nehaden[npfjets] = jet.neutralHadronEnergy();
	 pfjet_phen   [npfjets] = jet.photonEnergy();
	 pfjet_elen   [npfjets] = jet.electronEnergy();
	 pfjet_muen   [npfjets] = jet.muonEnergy();
	 //pfjet_hfhaden[npfjets] = jet.HFHadronEnergy();
	 //pfjet_hfemen [npfjets] = jet.HFEMEnergy(); 
	 //pfjet_chemf  [npfjets] = jet.chargedEmEnergyFraction();
	 //pfjet_chmuenf[npfjets] = jet.chargedMuEnergyFraction();

	 pfjet_chenf[npfjets] = jet.chargedHadronEnergyFraction();
	 pfjet_nuenf[npfjets] = jet.neutralHadronEnergyFraction();
	 pfjet_phenf[npfjets] = jet.photonEnergyFraction();
	 pfjet_elenf[npfjets] = jet.electronEnergyFraction();
	 pfjet_muenf[npfjets] = jet.muonEnergyFraction();
	 //pfjet_hfhadenf[npfjets] = jet.HFHadronEnergyFraction();
	 //pfjet_hfemenf [npfjets] = jet.HFEMEnergyFraction();
	 
	 pfjet_ieta [npfjets] = getieta(jet.eta());
	 pfjet_iphi [npfjets] = getiphi(jet.phi());

	 //! pf candidites
	 for (int ix=0; ix<pfjet_cons[npfjets]; ix++) {
	   //const reco::PFCandidatePtr candPtr = jet.getPFConstituent(ix);
	   //reco::PFCandidate cand(candPtr);

	   reco::PFCandidate cand  = (jet.getPFConstituents())[ix];

	   pfjet_candId [npfjets][ix] = (int)cand.particleId();
	   pfjet_cande  [npfjets][ix] = cand.energy();
	   pfjet_candet [npfjets][ix] = cand.et();
	   pfjet_candpt [npfjets][ix] = cand.pt();
	   pfjet_candeta[npfjets][ix] = cand.eta();
	   pfjet_candphi[npfjets][ix] = cand.phi();
	   pfjet_tracksumecal[npfjets] += 0;
	   pfjet_tracksumhcal[npfjets] += 0;
	   pfjet_tracksumho[npfjets] += 0;

	   //only charged hadrons and leptons can be asscociated with a track
	    float cand_type = cand.particleId();
	    if(cand_type == PFCandidate::h || cand_type == PFCandidate::e || cand_type == PFCandidate::mu){

	      for(unsigned iblock=0; iblock<cand.elementsInBlocks().size(); iblock++) {
	        PFBlockRef blockRef   = cand.elementsInBlocks()[iblock].first;
	        unsigned indexInBlock = cand.elementsInBlocks()[iblock].second;

	        const edm::OwnVector<  reco::PFBlockElement>&  elements = (*blockRef).elements();

	        switch (elements[indexInBlock].type()) {
	    	 //This tells you what type of element it is:
	    	 //cout<<" block type"<<elements[indexInBlock].type()<<endl; 
	        case PFBlockElement::ECAL: {
	    	 reco::PFClusterRef clusterRef = elements[indexInBlock].clusterRef();
	    	 double eet = clusterRef->energy()/cosh(clusterRef->eta());
	    	 //cout<<" ecal energy "<<clusterRef->energy()<<endl;
	    	 pfjet_tracksumecal[npfjets] += eet;
	    	 break;
	        }
	        case PFBlockElement::HCAL: {
	    	 reco::PFClusterRef clusterRef = elements[indexInBlock].clusterRef();
	    	 double eet = clusterRef->energy()/cosh(clusterRef->eta());
	    	 //cout<<" hcal energy "<<clusterRef->energy()<<endl;
	    	 pfjet_tracksumhcal[npfjets] += eet;
	    	 break;
	        }
	        case PFBlockElement::HO: {
	    	 reco::PFClusterRef clusterRef = elements[indexInBlock].clusterRef();
	    	 double eet = clusterRef->energy()/cosh(clusterRef->eta());
	    	 //cout<<" ho energy "<<clusterRef->energy()<<endl;
	    	 pfjet_tracksumho[npfjets] += eet;
	    	 break;
	        }
	        case PFBlockElement::TRACK: {
	    	 //This is just the reference to the track itself, since tracks can never be linked to other tracks
	    	 break;
	        }
	        default:
	    	 break;
	        }
   	     }
	   }
	 }
	 //std::cout <<"\t \t PF Jets  "<< "\t corrpt : "<< jet.pt() <<"\t raw pt : "<< pfjet_ucpt[npfjets] << "\t ieta : " << pfjet_ieta [npfjets] << std::endl;
	 npfjets++;
       }
       if(isHistFill){
	 if(istat==1)hist_nJet->Fill(ncalojets);
	 else if(istat==2)hist_pf_nJet->Fill(npfjets);
       }
     }
   }//! jet algos completed


   // Calo Geometry
   edm::ESHandle<CaloGeometry> pG;
   const CaloGeometry* caloGeom;

   iSetup.get<CaloGeometryRecord>().get(pG);
   caloGeom = pG.product();

   //! HBHE RecHit Collection
   edm::Handle<HBHERecHitCollection> hbheht;
   iEvent.getByLabel("hbhereco","",hbheht);
   nhb=0;
   if (hbheht.isValid()) {
     if ((*hbheht).size()>0) {
       for (HBHERecHitCollection::const_iterator hbheit=(*hbheht).begin(); hbheit!=(*hbheht).end(); ++hbheit){

	 HcalDetId id  = (*hbheit).id();
	 if( id.subdet() != HcalBarrel ) continue;

	 float henr    = (*hbheit).energy();
	 float htime   = (*hbheit).time();
	 int hieta     = id.ieta();
	 int hiphi     = id.iphi();
	 int hdepth    = id.depth();

	 GlobalPoint pos = caloGeom->getPosition(hbheit->detid());
	 math::XYZPoint posV(pos.x(), pos.y(), pos.z());
	 
	 hb_en  [nhb]  = henr;
	 hb_et  [nhb]  = henr*sin(posV.theta());
	 hb_ti  [nhb]  = htime;
	 hb_ieta[nhb]  = hieta; 
	 hb_iphi[nhb]  = hiphi;
	 hb_depth[nhb] = hdepth;
	 hb_ring[nhb]  = getring(hieta);
	 hb_eta [nhb]  = pos.eta();
	 hb_phi [nhb]  = pos.phi();
	 nhb++;
       }
     }
   }

   //! HO Rechit Collection
   edm::Handle<HORecHitCollection> hoht;
   iEvent.getByLabel("horeco","",hoht);
   //iEvent.getByLabel(hoLabel_,hoLabel2_,hoht);

   nho=0;
   if (hoht.isValid()) {
     if ((*hoht).size()>0) {
       for (HORecHitCollection::const_iterator hoit=(*hoht).begin(); hoit!=(*hoht).end(); ++hoit){

	 HcalDetId id  = (*hoit).id();
	 float hoenr   = (*hoit).energy();
	 float hotime  = (*hoit).time();
	 int ietaho    = id.ieta();
	 int iphiho    = id.iphi();
	 int idepth    = id.depth();

	 GlobalPoint pos = caloGeom->getPosition(hoit->detid());

	 //if (ietaho== 5 && (iphiho==18 || iphiho==19)) continue;
	 //if (ietaho==-5 && (iphiho>=11 && iphiho<=14)) continue;
	 //if ((ietaho<-10) || (ietaho>10 && iphiho<59)) continue;

	 //math::XYZPoint posV(pos.x() - bestVtxX, pos.y() - bestVtxY, pos.z() - bestVtxZ);
	 math::XYZPoint posV(pos.x(), pos.y(), pos.z());

	 ho_en  [nho] = hoenr;
	 ho_et  [nho] = hoenr*sin(posV.theta());
	 ho_ti  [nho] = hotime;
	 ho_ieta[nho] = ietaho; 
	 ho_iphi[nho] = iphiho;
	 ho_depth[nho] = idepth;
	 ho_ring[nho] = getring(ietaho);
	 ho_eta [nho] = pos.eta();
	 ho_phi [nho] = pos.phi();
	 nho++;
       }
     }
   }
   if(npfjets>0 || ncalojets>0)T1->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
HOAnalyzer::beginJob( )
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HOAnalyzer::endJob() 
{
  // theFile->cd();
  // theFile->Write();
  // theFile->Close();
  cout <<"End of HOAnalyzer "<<Nevents<<endl;

}

// ------------ method called when starting to processes a run  ------------
void 
HOAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
HOAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
HOAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
HOAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HOAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HOAnalyzer);
