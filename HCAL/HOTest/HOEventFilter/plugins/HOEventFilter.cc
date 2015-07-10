// -*- C++ -*-
//
// Package:    HOTest/HOEventFilter
// Class:      HOEventFilter
// 
/**\class HOEventFilter HOEventFilter.cc HOTest/HOEventFilter/plugins/HOEventFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Fri, 13 Mar 2015 10:49:52 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

#include "FWCore/Framework/interface/ESHandle.h"


using namespace std;
using namespace edm;



//
// class declaration
//

class HOEventFilter : public edm::EDFilter {
   public:
      explicit HOEventFilter(const edm::ParameterSet&);
      ~HOEventFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
  virtual void beginJob() override;
  virtual bool filter(edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  // ----------member data ---------------------------
  edm::InputTag caloJetLabel_;
  edm::InputTag hoLabel_;
  double hoenthresh_;
  int nEvent;
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
HOEventFilter::HOEventFilter(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed
  caloJetLabel_   = iConfig.getUntrackedParameter<InputTag>("CaloJetTag");
  hoLabel_        = iConfig.getUntrackedParameter<InputTag>("hoInput");
  hoenthresh_     = iConfig.getParameter<double>("hoEnThresh");
  nEvent=0;
}


HOEventFilter::~HOEventFilter()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
HOEventFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  nEvent++;

  edm::Handle<reco::CaloJetCollection>CaloJetColl;
  iEvent.getByLabel(caloJetLabel_, CaloJetColl);
  if (CaloJetColl.isValid()) {
    reco::CaloJetCollection::const_iterator calojet;
    for ( calojet = CaloJetColl->begin(); calojet != CaloJetColl->end(); ++calojet ) {
      if ( calojet->hadEnergyInHO() >0 ) {
	std::cout <<" \t \t  HOEventFilter ::  processing  event "  << nEvent << " " << iEvent.id().event() 
		  << " HO en : " <<  calojet->hadEnergyInHO() << std::endl;
	return true;
      }
    }
  }


  // edm::Handle<HORecHitCollection> hoht;
  // iEvent.getByLabel(hoLabel_.label(),"",hoht);
  // if ( !hoht.isValid() )return false;
  // if ((*hoht).size()==0)return false;
  
  // float sumen=0.0;
  // for (HORecHitCollection::const_iterator hoit=(*hoht).begin(); hoit!=(*hoht).end(); ++hoit){
  //   sumen        += (*hoit).energy();
  // }
  // if( sumen <= hoenthresh_ )return false;
  
  return false;
}

// ------------ method called once each job just before starting event loop  ------------
void 
HOEventFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HOEventFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
HOEventFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
HOEventFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
HOEventFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
HOEventFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HOEventFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(HOEventFilter);
