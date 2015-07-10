
// -*- C++ -*-
//
// Package:    HOinPFAlgo
// Class:      HOinPFAlgo
// 
/**\class HOinPFAlgo HOinPFAlgo.cc Test/HOinPFAlgo/src/HOinPFAlgo.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Gobinda Majumder,22 1-031,+41227679681,
//         Created:  Thu Dec  1 13:09:13 CET 2011
// $Id$
//
//


// system include files
#include <memory>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <functional>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"


#include "DataFormats/ParticleFlowReco/src/PFBlock.cc"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrackFwd.h"


#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETFwd.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/JetReco/interface/JetID.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "CLHEP/Vector/LorentzVector.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
//#include "CalibFormats/HcalObjects/interface/HcalCalibrationWidths.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
//#include "CondFormats/HcalObjects/interface/HcalQIECoder.h"
//#include "CondFormats/HcalObjects/interface/HcalPedestal.h"
//#include "CondFormats/HcalObjects/interface/HcalPedestalWidth.h"
//#include "CondFormats/HcalObjects/interface/HcalElectronicsMap.h"
#include "FWCore/Framework/interface/ESHandle.h"

//#include "DataFormats/HcalCalibObjects/interface/HBHERecTime.h"
//#include "DataFormats/HcalCalibObjects/interface/HBHERecTimeCollection.h"

using namespace std;
using namespace edm;
using namespace reco;  
using namespace CLHEP;

static const int njetmx = 500;

struct JRA{
  int irun, ievt, ilumi;
  int nprim;

  int ncalojet;
  int calojetmul  [njetmx];
  int calojetmat  [njetmx];
  float calojetmom[njetmx];
  float calojetthe[njetmx];
  float calojeteta[njetmx];
  float calojetphi[njetmx];
  float calojeten [njetmx];
  float calojetumm[njetmx]; 
  float calojetemf[njetmx];
  float calojeteb [njetmx];
  float calojethb [njetmx];
  float calojetho [njetmx];

  int npfjet;
  int pfjetmul[njetmx];
  int pfjetmat[njetmx];
  float pfjetmom[njetmx];
  float pfjetthe[njetmx];
  float pfjeteta[njetmx];
  float pfjetphi[njetmx];
  float pfjeten [njetmx];
  float pfjetumm[njetmx];
  float pfjetunc[njetmx];
  float pfchghad[njetmx];
  float pfneuhad[njetmx];
  float pfneuemf[njetmx]; 

  float pfjetsecal[njetmx];
  float pfjetshcal[njetmx];
  float pfjetsho  [njetmx];

  float pfjetsecal_ch[njetmx];
  float pfjetshcal_ch[njetmx];
  float pfjetsho_ch  [njetmx];

  float pfjetsecal_ph[njetmx];
  float pfjetshcal_ph[njetmx];
  float pfjetsho_ph  [njetmx];

  float pfjetsecal_el[njetmx];
  float pfjetshcal_el[njetmx];
  float pfjetsho_el  [njetmx];

  float pfjetsecal_nh[njetmx];
  float pfjetshcal_nh[njetmx];
  float pfjetsho_nh  [njetmx];

  float pfjetsecal_mu[njetmx];
  float pfjetshcal_mu[njetmx];
  float pfjetsho_mu  [njetmx];

  //! Raw energies
  float pfjetsecal_r[njetmx];
  float pfjetshcal_r[njetmx];
  float pfjetsho_r  [njetmx];

  float pfjetsecal_ch_r[njetmx];
  float pfjetshcal_ch_r[njetmx];
  float pfjetsho_ch_r  [njetmx];

  float pfjetsecal_ph_r[njetmx];
  float pfjetshcal_ph_r[njetmx];
  float pfjetsho_ph_r  [njetmx];

  float pfjetsecal_el_r[njetmx];
  float pfjetshcal_el_r[njetmx];
  float pfjetsho_el_r  [njetmx];

  float pfjetsecal_nh_r[njetmx];
  float pfjetshcal_nh_r[njetmx];
  float pfjetsho_nh_r  [njetmx];

  float pfjetsecal_mu_r[njetmx];
  float pfjetshcal_mu_r[njetmx];
  float pfjetsho_mu_r  [njetmx];



  int   ngenjet;
  int   ngenpar [njetmx]; //number of particle with E>1 GeV;
  float genjetmom[njetmx];
  float genjetthe[njetmx];
  float genjeteta[njetmx];
  float genjetphi[njetmx];
  float genjeten [njetmx];

  int ngenlep;
  float genmet; 
  float genmetphi;
  float gencalomet; 
  float gencalometphi;
  float pfmet;
  float pfmetphi;
  float pfmetsig;
  float pfchmet;
  float pfchmetphi;
  float pfchmetsig;
  float calomet;
  float calometphi;
  float tcmet;
  float tcmetphi;

};



//
// class declaration
//

class HOinPFAlgo : public edm::EDAnalyzer {
   public:
      explicit HOinPFAlgo(const edm::ParameterSet&);
      ~HOinPFAlgo();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  // ----------member data ---------------------------
  edm::Service<TFileService> fs;


  int Nevt;
  bool doGenJets_;
  bool doGenMets_;
  edm::InputTag genjetLabel_;
  edm::InputTag genmetLabel_;
  bool doCaloJets_;
  bool doCaloMets_;
  bool doPFJets_;
  bool doPFMets_;
  bool doTCMets_;
  edm::InputTag jetidLabel_;
  edm::InputTag calojetLabel_;
  edm::InputTag calometLabel_;
  edm::InputTag pfjetLabel_;
  edm::InputTag pfmetLabel_;
  edm::InputTag pfchmetLabel_;
  edm::InputTag tcmetLabel_;
  int tagHO_;
  double jtEtaCut_;
  double jtPtCut_;
  double genPtCut_;
  double genEtaCut_;
  double dRmatch_;
  bool   doHist_;
  bool   doGenMatch_;

  TH1F *h_calojtpt ;
  TH1F *h_calojteta;
  TH1F *h_calojtphi;

  TH1F *h_pfjtpt ;
  TH1F *h_pfjteta;
  TH1F *h_pfjtphi;

  TH1F *h_genjtpt ;
  TH1F *h_genjteta;
  TH1F *h_genjtphi;

  TH1F *h_genmet ;
  TH1F *h_genmetphi;
  TH1F *h_gencalomet ;
  TH1F *h_gencalometphi;

  TH1F *h_calomet ;
  TH1F *h_calometphi;
  TH1F *h_tcmet ;
  TH1F *h_tcmetphi;

  TH1F *h_pfmet ;
  TH1F *h_pfmetphi;
  TH1F *h_pfchmet ;
  TH1F *h_pfchmetphi;

  TTree* T1;

  JRA jets_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HOinPFAlgo::HOinPFAlgo(const edm::ParameterSet& pset)

{
  //now do what ever initialization is needed
  doGenJets_    = pset.getParameter<bool>("doGenJets");
  doGenMets_    = pset.getParameter<bool>("doGenMets");
  genjetLabel_  = pset.getUntrackedParameter<InputTag>("GenJetTag");
  genmetLabel_  = pset.getUntrackedParameter<InputTag>("GenMetTag");
  doCaloJets_   = pset.getParameter<bool>("doCaloJets");
  doCaloMets_   = pset.getParameter<bool>("doCaloMets");
  doPFJets_     = pset.getParameter<bool>("doPFJets");
  doPFMets_     = pset.getParameter<bool>("doPFMets");
  doTCMets_     = pset.getParameter<bool>("doTCMets");
  jetidLabel_   = pset.getParameter<InputTag>("JetIdTag");
  calojetLabel_ = pset.getUntrackedParameter<InputTag>("CaloJetTag");
  calometLabel_ = pset.getUntrackedParameter<InputTag>("CaloMetTag");
  pfjetLabel_   = pset.getUntrackedParameter<InputTag>("PFJetTag");
  pfmetLabel_   = pset.getUntrackedParameter<InputTag>("PFMetTag");
  pfchmetLabel_ = pset.getUntrackedParameter<InputTag>("PFChMetTag");
  tcmetLabel_   = pset.getUntrackedParameter<InputTag>("TCMetTag");
  jtPtCut_      = pset.getParameter<double>("JetPtCut");
  jtEtaCut_     = pset.getParameter<double>("JetEtaCut");
  doHist_       = pset.getParameter<bool>("doHist");
  doGenMatch_   = pset.getParameter<bool>("doGenMatch");
  genPtCut_     = pset.getParameter<double>("GenPtCut");
  genEtaCut_    = pset.getParameter<double>("GenEtaCut");

  dRmatch_ = 0.2;

  Nevt=0;

}


HOinPFAlgo::~HOinPFAlgo()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
HOinPFAlgo::analyze(const edm::Event& iEvent, const edm::EventSetup& pset)
{

  Nevt++;
  using namespace edm;
  
  bool isTreeFill = true;

  //std::cout <<"\t HOinPFAlgo::analyze "<< Nevt<<" "<<iEvent.id().run()<<" "<<iEvent.id().event()<<" "<<tagHO_<<" "<<jtPtCut_<< std::endl;
  //if (Nevt%100==1) cout <<"HOinPFAlgo::analyze "<< Nevt<<" "<<iEvent.id().run()<<" "<<iEvent.id().event()<<" "<<tagHO_<<" "<<jtPtCut_<<endl;
  
  jets_.irun  = iEvent.id().run();
  jets_.ievt  = iEvent.id().event();
  jets_.ilumi = iEvent.id().luminosityBlock();


  reco::TrackBase::Point beamPoint(0,0, 0);
  edm::Handle<reco::BeamSpot> beamSpotH;
  iEvent.getByLabel("offlineBeamSpot",beamSpotH);
  if (beamSpotH.isValid()){
    beamPoint = beamSpotH->position();
  }


  Handle<VertexCollection> primaryVertices;  
   try {
     iEvent.getByLabel("offlinePrimaryVertices",primaryVertices);
   } catch ( cms::Exception &iEvent ) {;}
   if (primaryVertices.isValid()) {
     int ndofct=0;
     int nchict=0;
     int nvert = 0;
     for (reco::VertexCollection::const_iterator vert=primaryVertices->begin(); vert<primaryVertices->end(); vert++) {
       nvert++;
       if (vert->isValid() && !vert->isFake()) {
	 if (vert->ndof()>7) { 
	   ndofct++;
	   if (vert->normalizedChi2()<5) nchict++;
	 }
       }
     }
     jets_.nprim = min(99,nvert) + 100*min(99,ndofct) + 10000*min(99,nchict);
   } else { jets_.nprim = 0;}

   //std::cout <<" \t # of primary vertices  : " << jets_.nprim << std::endl;

   if( doGenJets_ || doGenMatch_ ){
     Handle<reco::GenJetCollection> GenJetColl;
     iEvent.getByLabel(genjetLabel_, GenJetColl);
     jets_.ngenjet = 0;
     if (GenJetColl.isValid()) {
       reco::GenJetCollection::const_iterator genjet;
       for(genjet = GenJetColl->begin(); genjet != GenJetColl->end(); ++genjet) {
	 if ( fabs( genjet->eta() ) > genEtaCut_ ) continue;

	 Hep3Vector tmp3v( genjet->px(), genjet->py(), genjet->pz());
	 double jetpt=tmp3v.perp();
	 if ( jetpt < genPtCut_ ) continue;

	 h_genjtpt ->Fill(genjet->pt());
         h_genjteta->Fill(genjet->eta());
         h_genjtphi->Fill(genjet->phi());
	 
	 jets_.ngenpar  [jets_.ngenjet] = genjet->getGenConstituents().size();
	 jets_.genjetmom[jets_.ngenjet] = tmp3v.mag();
	 jets_.genjetthe[jets_.ngenjet] = tmp3v.theta();
	 jets_.genjeteta[jets_.ngenjet] = genjet->eta();
	 jets_.genjetphi[jets_.ngenjet] = tmp3v.phi();
	 jets_.genjeten [jets_.ngenjet] = genjet->energy();
	 jets_.ngenjet++;
	 
	 if (jets_.ngenjet >= njetmx){
	   isTreeFill=false;
	   break;
	 }
       }
     }
     //std::cout <<" \t # of gen jets  : " << jets_.ngenjet << std::endl;
   }

   if( doGenMets_ ){
     jets_.ngenlep = 0;

     Handle <GenParticleCollection> GenParticleColl;
     iEvent.getByLabel("genParticles",GenParticleColl);
     if (GenParticleColl.isValid()) {
       for (GenParticleCollection::const_iterator cand = GenParticleColl->begin(); cand < GenParticleColl->end(); cand++) {
	 if (cand->status()!=1) continue;
	 if (cand->pt()<5) continue;
	 int pdgid = abs(cand->pdgId())-11;
	 if (pdgid<0 || pdgid>5) continue;
	 jets_.ngenlep +=pow(10., pdgid);
       }
     }

     jets_.gencalomet = jets_.gencalometphi = -100;
     edm::Handle<GenMETCollection> GenCaloMetColl;
     iEvent.getByLabel(genmetLabel_ , GenCaloMetColl);
     if (GenCaloMetColl.isValid()) {
       jets_.gencalomet    = GenCaloMetColl->begin()->et();
       jets_.gencalometphi = GenCaloMetColl->begin()->phi();
       h_gencalomet->Fill(jets_.gencalomet);
       h_gencalometphi->Fill(jets_.gencalometphi);
     }

     jets_.genmet = jets_.genmetphi = -100;
     edm::Handle<GenMETCollection> GenMetColl;
     iEvent.getByLabel("genMetTrue" , GenMetColl);
     if (GenMetColl.isValid()) {
       jets_.genmet    = GenMetColl->begin()->et();
       jets_.genmetphi = GenMetColl->begin()->phi();
       h_genmet->Fill(jets_.genmet);
       h_genmetphi->Fill(jets_.genmetphi);
     }
   }


   if( doCaloMets_ ){
     edm::Handle<CaloMETCollection> CaloMetColl;
     //iEvent.getByLabel("corMetGlobalMuons18", metsis);
     iEvent.getByLabel(calometLabel_, CaloMetColl);
     jets_.calomet = jets_.calometphi = -100;
     if (CaloMetColl.isValid()) {
       jets_.calomet  = CaloMetColl->begin()->et();
       jets_.calometphi = CaloMetColl->begin()->phi();
       h_calomet ->Fill(jets_.calomet);
       h_calometphi->Fill(jets_.calometphi);
     }
   }


   if( doTCMets_ ){
     edm::Handle<METCollection> TCMetColl;
     iEvent.getByLabel(tcmetLabel_, TCMetColl);
     jets_.tcmet = jets_.tcmetphi = -100;
     if (TCMetColl.isValid()) {
       jets_.tcmet  = TCMetColl->begin()->et();
       jets_.tcmetphi = TCMetColl->begin()->phi();
       h_tcmet ->Fill(jets_.tcmet);
       h_tcmetphi->Fill(jets_.tcmetphi);
     }
   }

   if( doPFMets_ ){
     edm::Handle<PFMETCollection> PFMetColl;
     iEvent.getByLabel(pfmetLabel_, PFMetColl);
     jets_.pfmet = jets_.pfmetphi = -100;
     if (PFMetColl.isValid()) {
       jets_.pfmet     = PFMetColl->begin()->et();
       jets_.pfmetphi  = PFMetColl->begin()->phi();
       jets_.pfmetsig  = PFMetColl->begin()->significance();
       h_pfmet ->Fill(jets_.pfmet);
       h_pfmetphi->Fill(jets_.pfmetphi);
     }
     edm::Handle<PFMETCollection> PFChMetColl;
     iEvent.getByLabel(pfchmetLabel_, PFChMetColl);
     jets_.pfchmet = jets_.pfchmetphi = -100;
     if (PFChMetColl.isValid()) {
       jets_.pfchmet    = PFChMetColl->begin()->et();
       jets_.pfchmetphi = PFChMetColl->begin()->phi();
       jets_.pfchmetsig = PFChMetColl->begin()->significance();
       h_pfchmet ->Fill(jets_.pfchmet);
       h_pfchmetphi->Fill(jets_.pfchmetphi);
     }
   }


   //bool hasHO=false;
   if ( doCaloJets_ ) {

     //! ak4CaloJet
     edm::Handle<edm::ValueMap<reco::JetID> > caloIds;
     iEvent.getByLabel(jetidLabel_, caloIds);
     
     edm::Handle<reco::CaloJetCollection>CaloJetColl;
     iEvent.getByLabel(calojetLabel_, CaloJetColl);
     
     jets_.ncalojet=0;
     if (CaloJetColl.isValid()) {
       reco::CaloJetCollection::const_iterator calojet;
       for ( calojet = CaloJetColl->begin(); calojet != CaloJetColl->end(); ++calojet ) {
	 
	 if (fabs(calojet->eta()) > jtEtaCut_ || calojet->pt() < jtPtCut_) continue;

	 int ipass = 1;
	 if (calojet->nConstituents() < 2) continue;
	 if (fabs(calojet->eta())<genEtaCut_ && calojet->emEnergyFraction() <0.01) ipass = 0; //continue;
	 int ipasstight = ipass;
	 if (fabs(calojet->eta())<genEtaCut_ && calojet->emEnergyFraction() >=1.0) ipasstight = 0; //continue;


	 h_calojtpt ->Fill(calojet->pt());
         h_calojteta->Fill(calojet->eta());
         h_calojtphi->Fill(calojet->phi());


      	 Hep3Vector tmp3v1(calojet->px(), calojet->py(), calojet->pz());
      	 jets_.calojetmom  [jets_.ncalojet] = tmp3v1.mag();
      	 jets_.calojetthe  [jets_.ncalojet] = tmp3v1.theta();
      	 jets_.calojeteta  [jets_.ncalojet] = calojet->eta();
      	 jets_.calojetphi  [jets_.ncalojet] = tmp3v1.phi();
      	 //jets_.calojeten   [jets_.ncalojet] = (2*ipass - 1)*calojet->energy();
	 jets_.calojeten   [jets_.ncalojet] = calojet->energy();
      	 jets_.calojetmul  [jets_.ncalojet] = calojet->nConstituents();
      	 if (ipasstight==0) jets_.calojetmul[jets_.ncalojet] = -jets_.calojetmul[jets_.ncalojet];
     	 
      	 jets_.calojetemf[jets_.ncalojet] = calojet->emEnergyFraction();
      	 jets_.calojeteb [jets_.ncalojet] = calojet->emEnergyInEB();
      	 jets_.calojethb [jets_.ncalojet] = calojet->hadEnergyInHB();
      	 jets_.calojetho [jets_.ncalojet] = calojet->hadEnergyInHO();

	 double energy = calojet->energy();
	 //double energy = 0;
      	 // for (int ij=0; ij< calojet->nConstituents(); ij++) {
	 //   energy += calojet->getCaloConstituent(ij)->p();
      	 // }
      	 jets_.calojetumm[jets_.ncalojet] = energy;
	 jets_.calojetmat[jets_.ncalojet] = -1;
	 // if( calojet->hadEnergyInHO() > 0 ){
	 //   //std::cout << "\t \t " <<  jets_.ncalojet <<  " Calo Jet uncorrected  energy : " << energy << "  energy :  "  << calojet->energy() << std::endl;
	 //   std::cout << "\t \t Calojet : " <<  jets_.ncalojet <<  " umm  energy : " << energy << "  HO energy :  "  << calojet->hadEnergyInHO() << std::endl;
	 // }
      	 jets_.ncalojet++;
      	 if (jets_.ncalojet>=njetmx) {
	   isTreeFill = false;
	   break;
	 }
       }
     }
     //std::cout <<" \t # of calo jets  : " << jets_.ncalojet << std::endl;
   }

   
   if( doPFJets_ ){     

     //! ak4PFJets
     edm::Handle<reco::PFJetCollection> PFJetColl;
     jets_.npfjet = 0;
     iEvent.getByLabel(pfjetLabel_, PFJetColl);

     if (PFJetColl.isValid()) {
       //cout <<"  Jet container : " << pfjetLabel_.label() << endl;
       reco::PFJetCollection::const_iterator pfjet;
       for ( pfjet = PFJetColl->begin(); pfjet != PFJetColl->end(); ++pfjet ) {	 

         if ( fabs(pfjet->eta()) > jtEtaCut_ || pfjet->pt() < jtPtCut_ ) continue;

	 h_pfjtpt ->Fill(pfjet->pt());
         h_pfjteta->Fill(pfjet->eta());
         h_pfjtphi->Fill(pfjet->phi());
	 
	 //if(hasHO)std::cout <<"   \t " <<  pfjetLabel_.label() << "  pt : "  << pfjet->pt() << " en : " << pfjet->energy() << std::endl;
	 //std::cout <<" \t \t PFJet " <<  pfjetLabel_.label() << "  pt : "  << pfjet->pt() << " en : " << pfjet->energy() << std::endl;
	 
         HepLorentzVector jet4v( pfjet->px(), pfjet->py(), pfjet->pz(), pfjet->p() );

	 jets_.pfjetunc[jets_.npfjet] = 0.0;
         jets_.pfjetumm[jets_.npfjet] = 0.0;
         jets_.pfjetmom[jets_.npfjet] = jet4v.rho();
         jets_.pfjetthe[jets_.npfjet] = jet4v.theta();
         jets_.pfjeteta[jets_.npfjet] = pfjet->eta();
         jets_.pfjetphi[jets_.npfjet] = jet4v.phi();
	 
         jets_.pfneuemf[jets_.npfjet] = pfjet->neutralEmEnergyFraction();
         jets_.pfchghad[jets_.npfjet] = pfjet->chargedHadronEnergyFraction();
         jets_.pfneuhad[jets_.npfjet] = pfjet->neutralHadronEnergyFraction();

	 jets_.pfjeten[jets_.npfjet] = pfjet->energy();
	 
	 jets_.pfjetmul[jets_.npfjet] = min(9,pfjet->muonMultiplicity())
           + 10*min(99,pfjet->chargedHadronMultiplicity())
           + 1000*min(99,pfjet->neutralHadronMultiplicity())
           + 100000*min(99,pfjet->electronMultiplicity())
           + 10000000*min(99,pfjet->photonMultiplicity ());
	 
         HepLorentzVector const4v(0,0,0,0);
	 float sumHO=0.0, sumHO_raw=0.0;
	 float sumHCAL=0.0, sumHCAL_raw=0.0;
	 float sumECAL=0.0, sumECAL_raw=0.0;

	 float sumHO_ch=0.0, sumHO_ch_raw=0.0;
	 float sumHCAL_ch=0.0, sumHCAL_ch_raw=0.0;
	 float sumECAL_ch=0.0, sumECAL_ch_raw=0.0;

	 float sumHO_ph=0.0, sumHO_ph_raw=0.0;
	 float sumHCAL_ph=0.0, sumHCAL_ph_raw=0.0;
	 float sumECAL_ph=0.0, sumECAL_ph_raw=0.0;

	 float sumHO_nh=0.0, sumHO_nh_raw=0.0;
	 float sumHCAL_nh=0.0, sumHCAL_nh_raw=0.0;
	 float sumECAL_nh=0.0, sumECAL_nh_raw=0.0;

	 float sumHO_el=0.0, sumHO_el_raw=0.0;
	 float sumHCAL_el=0.0, sumHCAL_el_raw=0.0;
	 float sumECAL_el=0.0, sumECAL_el_raw=0.0;

	 float sumHO_mu=0.0, sumHO_mu_raw=0.0;
	 float sumHCAL_mu=0.0, sumHCAL_mu_raw=0.0;
	 float sumECAL_mu=0.0, sumECAL_mu_raw=0.0;

	 for (int ix=0; ix< pfjet->nConstituents(); ix++) {
           const reco::PFCandidatePtr pfcand = pfjet->getPFConstituent (ix);
	   if(pfcand->particleId() == reco::PFCandidate::X || pfcand->energy() == 0)continue;
           const4v += HepLorentzVector(pfcand->px(), pfcand->py(), pfcand->pz(), pfcand->energy());

	   sumECAL += pfcand->ecalEnergy();
	   sumHCAL += pfcand->hcalEnergy();
	   sumHO   += pfcand->hoEnergy();

	   sumECAL_raw += pfcand->rawEcalEnergy();
	   sumHCAL_raw += pfcand->rawHcalEnergy();
	   sumHO_raw   += pfcand->rawHoEnergy();

	   if( pfcand->particleId() == reco::PFCandidate::h ){
	     sumECAL_ch += pfcand->ecalEnergy();
	     sumHCAL_ch += pfcand->hcalEnergy();
	     sumHO_ch   += pfcand->hoEnergy();
	     
	     sumECAL_ch_raw += pfcand->rawEcalEnergy();
	     sumHCAL_ch_raw += pfcand->rawHcalEnergy();
	     sumHO_ch_raw   += pfcand->rawHoEnergy();
	   }else if(pfcand->particleId() == reco::PFCandidate::h0){
	     sumECAL_nh += pfcand->ecalEnergy();
	     sumHCAL_nh += pfcand->hcalEnergy();
	     sumHO_nh   += pfcand->hoEnergy();
	     
	     sumECAL_nh_raw += pfcand->rawEcalEnergy();
	     sumHCAL_nh_raw += pfcand->rawHcalEnergy();
	     sumHO_nh_raw   += pfcand->rawHoEnergy();
	   }else if(pfcand->particleId() == reco::PFCandidate::e){
	     sumECAL_el += pfcand->ecalEnergy();
	     sumHCAL_el += pfcand->hcalEnergy();
	     sumHO_el   += pfcand->hoEnergy();
	     
	     sumECAL_el_raw += pfcand->rawEcalEnergy();
	     sumHCAL_el_raw += pfcand->rawHcalEnergy();
	     sumHO_el_raw   += pfcand->rawHoEnergy();
	   }else if(pfcand->particleId() == reco::PFCandidate::gamma){
	     sumECAL_ph += pfcand->ecalEnergy();
	     sumHCAL_ph += pfcand->hcalEnergy();
	     sumHO_ph   += pfcand->hoEnergy();
	     
	     sumECAL_ph_raw += pfcand->rawEcalEnergy();
	     sumHCAL_ph_raw += pfcand->rawHcalEnergy();
	     sumHO_ph_raw   += pfcand->rawHoEnergy();
	   }else if(pfcand->particleId() == reco::PFCandidate::mu){
	     sumECAL_mu += pfcand->ecalEnergy();
	     sumHCAL_mu += pfcand->hcalEnergy();
	     sumHO_mu   += pfcand->hoEnergy();
	     
	     sumECAL_mu_raw += pfcand->rawEcalEnergy();
	     sumHCAL_mu_raw += pfcand->rawHcalEnergy();
	     sumHO_mu_raw   += pfcand->rawHoEnergy();
	   }
	 }

         jets_.pfjetunc[jets_.npfjet] =const4v.e();
         jets_.pfjetumm[jets_.npfjet] =const4v.rho();
	 jets_.pfjetmat[jets_.npfjet] = -1;

	 jets_.pfjetsecal[jets_.npfjet] = sumECAL;
	 jets_.pfjetshcal[jets_.npfjet] = sumHCAL;
	 jets_.pfjetsho  [jets_.npfjet] = sumHO;

	 jets_.pfjetsecal_r[jets_.npfjet] = sumECAL_raw;
	 jets_.pfjetshcal_r[jets_.npfjet] = sumHCAL_raw;
	 jets_.pfjetsho_r  [jets_.npfjet] = sumHO_raw;

	 jets_.pfjetsecal_ch[jets_.npfjet] = sumECAL_ch;
	 jets_.pfjetshcal_ch[jets_.npfjet] = sumHCAL_ch;
	 jets_.pfjetsho_ch  [jets_.npfjet] = sumHO_ch;

	 jets_.pfjetsecal_ch_r[jets_.npfjet] = sumECAL_ch_raw;
	 jets_.pfjetshcal_ch_r[jets_.npfjet] = sumHCAL_ch_raw;
	 jets_.pfjetsho_ch_r  [jets_.npfjet] = sumHO_ch_raw;

	 jets_.pfjetsecal_el[jets_.npfjet] = sumECAL_el;
	 jets_.pfjetshcal_el[jets_.npfjet] = sumHCAL_el;
	 jets_.pfjetsho_el  [jets_.npfjet] = sumHO_el;

	 jets_.pfjetsecal_el_r[jets_.npfjet] = sumECAL_el_raw;
	 jets_.pfjetshcal_el_r[jets_.npfjet] = sumHCAL_el_raw;
	 jets_.pfjetsho_el_r  [jets_.npfjet] = sumHO_el_raw;

	 jets_.pfjetsecal_ph[jets_.npfjet] = sumECAL_ph;
	 jets_.pfjetshcal_ph[jets_.npfjet] = sumHCAL_ph;
	 jets_.pfjetsho_ph  [jets_.npfjet] = sumHO_ph;

	 jets_.pfjetsecal_ph_r[jets_.npfjet] = sumECAL_ph_raw;
	 jets_.pfjetshcal_ph_r[jets_.npfjet] = sumHCAL_ph_raw;
	 jets_.pfjetsho_ph_r  [jets_.npfjet] = sumHO_ph_raw;

	 jets_.pfjetsecal_nh[jets_.npfjet] = sumECAL_nh;
	 jets_.pfjetshcal_nh[jets_.npfjet] = sumHCAL_nh;
	 jets_.pfjetsho_nh  [jets_.npfjet] = sumHO_nh;

	 jets_.pfjetsecal_nh_r[jets_.npfjet] = sumECAL_nh_raw;
	 jets_.pfjetshcal_nh_r[jets_.npfjet] = sumHCAL_nh_raw;
	 jets_.pfjetsho_nh_r  [jets_.npfjet] = sumHO_nh_raw;

	 jets_.pfjetsecal_mu[jets_.npfjet] = sumECAL_mu;
	 jets_.pfjetshcal_mu[jets_.npfjet] = sumHCAL_mu;
	 jets_.pfjetsho_mu  [jets_.npfjet] = sumHO_mu;

	 jets_.pfjetsecal_mu_r[jets_.npfjet] = sumECAL_mu_raw;
	 jets_.pfjetshcal_mu_r[jets_.npfjet] = sumHCAL_mu_raw;
	 jets_.pfjetsho_mu_r  [jets_.npfjet] = sumHO_mu_raw;

	 jets_.npfjet++;

         if (jets_.npfjet >= njetmx) {
	   isTreeFill = false;
	   break;
	 }
       }
     }
     //std::cout <<" \t # of PF jets  : " << jets_.npfjet << std::endl;
   }


   if( doGenMatch_ ){

     //cout << "  # of calojets : " << jets_.ncalojet << " # of pfjets : " << jets_.npfjet << " # of gen jets : "<< jets_.ngenjet << endl;
     for(int igen  = 0; igen  < jets_.ngenjet; igen++){
       double dRmin=1000.0;
       int imatch = -1;
        for(int icalo = 0; icalo < jets_.ncalojet; icalo++){
	  float delr = deltaR(jets_.calojeteta[icalo], jets_.calojetphi[icalo], jets_.genjeteta[igen], jets_.genjetphi[igen]);
	  if( delr < dRmin && delr < dRmatch_ ){
	    dRmin = delr;
	    imatch = icalo;
	  }
        }
	if( imatch > -1 )jets_.calojetmat[imatch] = igen;
        //if(imatch>=0 && dRmin < dRmatch_ ) cout << " $$$$  Calo Gen-Matched jets : " <<  igen << "  gen : " << ( jets_.genjeten[igen] * sin(jets_.genjetthe[igen]) )
	// 					    //					       << " rec  : " << ( jets_.calojetmom[imatch] * sin(jets_.pfjetthe[imatch]) ) 
       // 					       << " rec  : " << jets_.calojeten[imatch]  
       // 					       << " eb   : " << jets_.calojeteb[imatch] 
       // 					       << " hb   : " << jets_.calojethb[imatch] 
       // 					       << " ho   : " << jets_.calojetho[imatch] << endl;
       dRmin=1000.0;
       imatch = -1;
       for(int ipf = 0; ipf < jets_.npfjet; ipf++){
	 float delr = deltaR(jets_.pfjeteta[ipf], jets_.pfjetphi[ipf], jets_.genjeteta[igen], jets_.genjetphi[igen]);
	 if( delr < dRmin && delr < dRmatch_ ){
	   dRmin = delr;
	   imatch = ipf;
	 }
       }
       // if(imatch>=0 && dRmin < dRmatch_ ) cout << " ****  PF   Gen-Matched jets : "   <<  igen 
       // 					       << "  gen : " <<  jets_.genjeten[igen] << "  rec : " <<  jets_.pfjeten [imatch] 
       // 					       << endl;
       if( imatch > -1 )jets_.pfjetmat[imatch] = igen;
     }
   }


   //if( jets_.ncalojet == 0 && jets_.npfjet == 0 )isTreeFill=false;
   if( jets_.ngenjet == 0 )isTreeFill=false;
   if( isTreeFill )T1->Fill();

   //cout << "\t " << endl;

}


// ------------ method called once each job just before starting event loop  ------------
void 
HOinPFAlgo::beginJob()
{
  
  string jetTitle = calojetLabel_.label();
  jetTitle.erase(3,4);
  jetTitle += " HOinPF";

  T1 = fs->make<TTree>("T1",jetTitle.c_str());
  T1->Branch("irun" , &jets_.irun , "irun/I");
  T1->Branch("ievt" , &jets_.ievt , "ievt/I");
  T1->Branch("ilumi", &jets_.ilumi, "ilumi/I");
  T1->Branch("nprim", &jets_.nprim, "nprim/I");
  
  if( doCaloJets_ ){
    //! Calo
    T1->Branch("ncalojet" , &jets_.ncalojet,"ncalojet/I");
    T1->Branch("calojetmul",jets_.calojetmul ,"calojetmul[ncalojet]/I");
    T1->Branch("calojetmat",jets_.calojetmat ,"calojetmat[ncalojet]/I");
    T1->Branch("calojetmom",jets_.calojetmom ,"calojetmom[ncalojet]/F");
    T1->Branch("calojetthe",jets_.calojetthe ,"calojeteta[ncalojet]/F");
    T1->Branch("calojeteta",jets_.calojeteta ,"calojetthe[ncalojet]/F");
    T1->Branch("calojetphi",jets_.calojetphi ,"calojetphi[ncalojet]/F");
    T1->Branch("calojetemf",jets_.calojetemf ,"calojetemf[ncalojet]/F");
    T1->Branch("calojeten" ,jets_.calojeten  ,"calojeten[ncalojet]/F");
    T1->Branch("calojetumm",jets_.calojetumm ,"calojetumm[ncalojet]/F");
    T1->Branch("calojeteb" ,jets_.calojeteb  ,"calojeteb[ncalojet]/F");
    T1->Branch("calojethb" ,jets_.calojethb  ,"calojethb[ncalojet]/F");
    T1->Branch("calojetho" ,jets_.calojetho  ,"calojetho[ncalojet]/F");
  }
  
  if( doPFJets_ ){
    //! PF jets
    T1->Branch("npfjet", &jets_.npfjet , "npfjet/I");
    T1->Branch("pfjetmul",jets_.pfjetmul ,"pfjetmul[npfjet]/I");
    T1->Branch("pfjetmat",jets_.pfjetmat ,"pfjetmat[npfjet]/I");
    T1->Branch("pfjeten" ,jets_.pfjeten  ,"pfjeten [npfjet]/F");
    T1->Branch("pfjetunc",jets_.pfjetunc ,"pfjetunc[npfjet]/F");
    T1->Branch("pfjetumm",jets_.pfjetumm ,"pfjetumm[npfjet]/F");
    T1->Branch("pfjetmom",jets_.pfjetmom ,"pfjetmom[npfjet]/F");
    T1->Branch("pfjetthe",jets_.pfjetthe ,"pfjetthe[npfjet]/F");
    T1->Branch("pfjeteta",jets_.pfjeteta ,"pfjeteta[npfjet]/F");
    T1->Branch("pfjetphi",jets_.pfjetphi ,"pfjetphi[npfjet]/F");
    T1->Branch("pfchghad",jets_.pfchghad ,"pfchghad[npfjet]/F");
    T1->Branch("pfneuemf",jets_.pfneuemf ,"pfneuemf[npfjet]/F");
    T1->Branch("pfneuhad",jets_.pfneuhad ,"pfneuhad[npfjet]/F");

    T1->Branch("pfjetsecal",jets_.pfjetsecal,"pfjetsecal[npfjet]/F");
    T1->Branch("pfjetshcal",jets_.pfjetshcal,"pfjetshcal[npfjet]/F");
    T1->Branch("pfjetsho",jets_.pfjetsho ,"pfjetsho[npfjet]/F");
    T1->Branch("pfjetsecal_r",jets_.pfjetsecal_r,"pfjetsecal_r[npfjet]/F");
    T1->Branch("pfjetshcal_r",jets_.pfjetshcal_r,"pfjetshcal_r[npfjet]/F");
    T1->Branch("pfjetsho_r",jets_.pfjetsho_r ,"pfjetsho_r[npfjet]/F");

    T1->Branch("pfjetsecal_ch",jets_.pfjetsecal_ch,"pfjetsecal_ch[npfjet]/F");
    T1->Branch("pfjetshcal_ch",jets_.pfjetshcal_ch,"pfjetshcal_ch[npfjet]/F");
    T1->Branch("pfjetsho_ch",jets_.pfjetsho_ch ,"pfjetsho_ch[npfjet]/F");
    T1->Branch("pfjetsecal_ch_r",jets_.pfjetsecal_ch_r,"pfjetsecal_ch_r[npfjet]/F");
    T1->Branch("pfjetshcal_ch_r",jets_.pfjetshcal_ch_r,"pfjetshcal_ch_r[npfjet]/F");
    T1->Branch("pfjetsho_ch_r",jets_.pfjetsho_ch_r ,"pfjetsho_ch_r[npfjet]/F");

    T1->Branch("pfjetsecal_nh",jets_.pfjetsecal_nh,"pfjetsecal_nh[npfjet]/F");
    T1->Branch("pfjetshcal_nh",jets_.pfjetshcal_nh,"pfjetshcal_nh[npfjet]/F");
    T1->Branch("pfjetsho_nh",jets_.pfjetsho_nh ,"pfjetsho_nh[npfjet]/F");
    T1->Branch("pfjetsecal_nh_r",jets_.pfjetsecal_nh_r,"pfjetsecal_nh_r[npfjet]/F");
    T1->Branch("pfjetshcal_nh_r",jets_.pfjetshcal_nh_r,"pfjetshcal_nh_r[npfjet]/F");
    T1->Branch("pfjetsho_nh_r",jets_.pfjetsho_nh_r ,"pfjetsho_nh_r[npfjet]/F");

    T1->Branch("pfjetsecal_ph",jets_.pfjetsecal_ph,"pfjetsecal_ph[npfjet]/F");
    T1->Branch("pfjetshcal_ph",jets_.pfjetshcal_ph,"pfjetshcal_ph[npfjet]/F");
    T1->Branch("pfjetsho_ph",jets_.pfjetsho_ph ,"pfjetsho_ph[npfjet]/F");
    T1->Branch("pfjetsecal_ph_r",jets_.pfjetsecal_ph_r,"pfjetsecal_ph_r[npfjet]/F");
    T1->Branch("pfjetshcal_ph_r",jets_.pfjetshcal_ph_r,"pfjetshcal_ph_r[npfjet]/F");
    T1->Branch("pfjetsho_ph_r",jets_.pfjetsho_ph_r ,"pfjetsho_ph_r[npfjet]/F");

    T1->Branch("pfjetsecal_el",jets_.pfjetsecal_el,"pfjetsecal_el[npfjet]/F");
    T1->Branch("pfjetshcal_el",jets_.pfjetshcal_el,"pfjetshcal_el[npfjet]/F");
    T1->Branch("pfjetsho_el",jets_.pfjetsho_el ,"pfjetsho_el[npfjet]/F");
    T1->Branch("pfjetsecal_el_r",jets_.pfjetsecal_el_r,"pfjetsecal_el_r[npfjet]/F");
    T1->Branch("pfjetshcal_el_r",jets_.pfjetshcal_el_r,"pfjetshcal_el_r[npfjet]/F");
    T1->Branch("pfjetsho_el_r",jets_.pfjetsho_el_r ,"pfjetsho_el_r[npfjet]/F");

    T1->Branch("pfjetsecal_mu",jets_.pfjetsecal_mu,"pfjetsecal_mu[npfjet]/F");
    T1->Branch("pfjetshcal_mu",jets_.pfjetshcal_mu,"pfjetshcal_mu[npfjet]/F");
    T1->Branch("pfjetsho_mu",jets_.pfjetsho_mu ,"pfjetsho_mu[npfjet]/F");
    T1->Branch("pfjetsecal_mu_r",jets_.pfjetsecal_mu_r,"pfjetsecal_mu_r[npfjet]/F");
    T1->Branch("pfjetshcal_mu_r",jets_.pfjetshcal_mu_r,"pfjetshcal_mu_r[npfjet]/F");
    T1->Branch("pfjetsho_mu_r",jets_.pfjetsho_mu_r ,"pfjetsho_mu_r[npfjet]/F");
  }
  
  if( doGenJets_ ){
    T1->Branch("ngenjet"  , &jets_.ngenjet,"ngenjet/I");
    T1->Branch("ngenpar" ,jets_.ngenpar ,"ngenpar[ngenjet]/I");
    T1->Branch("genjetmom",jets_.genjetmom,"genjetmom[ngenjet]/F");
    T1->Branch("genjetthe",jets_.genjetthe,"genjetthe[ngenjet]/F");
    T1->Branch("genjeteta",jets_.genjeteta,"genjeteta[ngenjet]/F");
    T1->Branch("genjetphi",jets_.genjetphi,"genjetphi[ngenjet]/F");
    T1->Branch("genjeten" ,jets_.genjeten ,"genjeten[ngenjet]/F");
  }


  if( doGenMets_ ){
    T1->Branch("ngenlep"  , &jets_.ngenlep  ,"ngenlep/I");
    T1->Branch("genmet"  , &jets_.genmet  ,"genmet/F");
    T1->Branch("genmetphi", &jets_.genmetphi,"genmetphi/F");
    T1->Branch("gencalomet"  , &jets_.gencalomet  ,"gencalomet/F");
    T1->Branch("gencalometphi", &jets_.gencalometphi,"gencalometphi/F");
  }

  if( doCaloMets_ ){
    T1->Branch("calomet", &jets_.calomet    ,"calomet/F");
    T1->Branch("calometphi" , &jets_.calometphi ,"calometphi/F");
  }
  if( doTCMets_ ){
    T1->Branch("tcmet", &jets_.tcmet    ,"tcmet/F");
    T1->Branch("tcmetphi" , &jets_.tcmetphi ,"tcmetphi/F");
  }

  if( doPFMets_ ){
    T1->Branch("pfmet"  , &jets_.pfmet  ,"pfmet/F");
    T1->Branch("pfmetphi" , &jets_.pfmetphi ,"pfmetphi/F");
    T1->Branch("pfmetsig" , &jets_.pfmetsig ,"pfmetsig/F");

    T1->Branch("pfchmet"  , &jets_.pfchmet  ,"pfchmet/F");
    T1->Branch("pfchmetphi" , &jets_.pfchmetphi ,"pfchmetphi/F");
    T1->Branch("pfchmetsig" , &jets_.pfchmetsig ,"pfchmetsig/F");
  }

  if( doHist_ ){

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
    TFileDirectory histDir = fs->mkdir("Histos");
    h_calojtpt  = histDir.make<TH1F>("hcalojtpt" , "Calojet pT ", 800, jtPtCut_ - 20., 4500);
    h_calojteta = histDir.make<TH1F>("hcalojteta", "Calojet eta", 72, -jtEtaCut_, jtEtaCut_);
    h_calojtphi = histDir.make<TH1F>("hcalojtphi", "Calojet phi", 72, -M_PI, M_PI);

    h_pfjtpt  = histDir.make<TH1F>("hpfjtpt" , "PFJet pT ", 800, jtPtCut_ - 20., 4500);
    h_pfjteta = histDir.make<TH1F>("hpfjteta", "PFJet eta", 72, -jtEtaCut_, jtEtaCut_);
    h_pfjtphi = histDir.make<TH1F>("hpfjtphi", "PFJet phi", 72, -M_PI, M_PI);

    h_genjtpt  = histDir.make<TH1F>("hgenjtpt" , "GenJet pT  ", 800, 100., 4500);
    h_genjteta = histDir.make<TH1F>("hgenjteta", "GenJet eta ", 72, -jtEtaCut_, jtEtaCut_);
    h_genjtphi = histDir.make<TH1F>("hgenjtphi", "GenJet phi ", 72, -M_PI, M_PI);


    h_genmet    = histDir.make<TH1F>("hgenmet"   , "True GenMet   ", 800, 0., 4500);
    h_genmetphi = histDir.make<TH1F>("hgenmetphi", "True GenMet phi ", 72, -M_PI, M_PI);
    h_gencalomet    = histDir.make<TH1F>("hgencalomet"   , "Calo GenMet   ", 800, 0., 4500);
    h_gencalometphi = histDir.make<TH1F>("hgencalometphi", "Calo GenMet phi ", 72, -M_PI, M_PI);

    h_calomet    = histDir.make<TH1F>("hcalomet" , "Calomet  ", 800, 0., 4500);
    h_calometphi = histDir.make<TH1F>("hcalometphi", "Calomet phi", 72, -M_PI, M_PI);

    h_tcmet    = histDir.make<TH1F>("htcmet" , "TCMET  ", 800, 0., 4500);
    h_tcmetphi = histDir.make<TH1F>("htcmetphi", "TCMET phi", 72, -M_PI, M_PI);

    h_pfmet  = histDir.make<TH1F>("hpfmet" , "PFmet ", 800, 0., 4500);
    h_pfmetphi = histDir.make<TH1F>("hpfmetphi", "PFmet phi", 72, -M_PI, M_PI);

    h_pfchmet  = histDir.make<TH1F>("hpfchmet" , "PFChmet  ", 800, 0., 4500);
    h_pfchmetphi = histDir.make<TH1F>("hpfchmetphi", "PFChmet phi", 72, -M_PI, M_PI);

  }



}

// ------------ method called once each job just after ending the event loop  ------------
void 
HOinPFAlgo::endJob() 
{
  cout<<"End of HOinPFAlgo with event "<<tagHO_<< "  " << Nevt<<endl;
}

// ------------ method called when starting to processes a run  ------------
void 
HOinPFAlgo::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
HOinPFAlgo::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
HOinPFAlgo::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
HOinPFAlgo::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HOinPFAlgo::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HOinPFAlgo);
