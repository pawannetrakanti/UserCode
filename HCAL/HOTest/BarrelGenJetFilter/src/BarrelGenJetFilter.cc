/*
cd /cvmfs/cms.cern.ch/slc6_amd64_gcc472/cms/cmssw/CMSSW_5_3_20/src/
diff RecoParticleFlow/PFClusterProducer/python/particleFlowRecHitHO_cfi.py  /cvmfs/cms.cern.ch/slc6_amd64_gcc481/cms/cmssw/CMSSW_7_2_1/src/RecoParticleFlow/PFClusterProducer/python/particleFlowRecHitHO_cfi.py  
diff RecoParticleFlow/PFClusterProducer/python/particleFlowCluster_cff.py   /cvmfs/cms.cern.ch/slc6_amd64_gcc481/cms/cmssw/CMSSW_7_2_1/src/RecoParticleFlow/PFClusterProducer/python/particleFlowCluster_cff.py   
diff RecoParticleFlow/PFProducer/python/particleFlowBlock_cff.py	    /cvmfs/cms.cern.ch/slc6_amd64_gcc481/cms/cmssw/CMSSW_7_2_1/src/RecoParticleFlow/PFProducer/python/particleFlowBlock_cff.py            
diff RecoParticleFlow/PFProducer/python/particleFlow_cff.py		    /cvmfs/cms.cern.ch/slc6_amd64_gcc481/cms/cmssw/CMSSW_7_2_1/src/RecoParticleFlow/PFProducer/python/particleFlow_cff.py                 
diff RecoJets/JetProducers/python/ak5PFJets_cfi.py			    /cvmfs/cms.cern.ch/slc6_amd64_gcc481/cms/cmssw/CMSSW_7_2_1/src/RecoJets/JetProducers/python/ak5PFJets_cfi.py                          
diff RecoMET/Configuration/python/RecoPFMET_cff.py			    /cvmfs/cms.cern.ch/slc6_amd64_gcc481/cms/cmssw/CMSSW_7_2_1/src/RecoMET/Configuration/python/RecoPFMET_cff.py                          
diff RecoParticleFlow/PFProducer/python/pfLinker_cff.py			    /cvmfs/cms.cern.ch/slc6_amd64_gcc481/cms/cmssw/CMSSW_7_2_1/src/RecoParticleFlow/PFProducer/python/pfLinker_cff.py                     
diff RecoParticleFlow/Configuration/python/RecoParticleFlow_cff.py	    /cvmfs/cms.cern.ch/slc6_amd64_gcc481/cms/cmssw/CMSSW_7_2_1/src/RecoParticleFlow/Configuration/python/RecoParticleFlow_cff.py 


diff RecoParticleFlow/PFClusterProducer/python/particleFlowRecHitHO_cfi.py  /home/gobinda/analho/CMSSW_5_3_20/src/RecoParticleFlow/PFClusterProducer/python/particleFlowRecHitHO_cfi.py  
diff RecoParticleFlow/PFClusterProducer/python/particleFlowCluster_cff.py   /home/gobinda/analho/CMSSW_5_3_20/src/RecoParticleFlow/PFClusterProducer/python/particleFlowCluster_cff.py   
diff RecoParticleFlow/PFProducer/python/particleFlowBlock_cff.py	    /home/gobinda/analho/CMSSW_5_3_20/src/RecoParticleFlow/PFProducer/python/particleFlowBlock_cff.py            
diff RecoParticleFlow/PFProducer/python/particleFlow_cff.py		    /home/gobinda/analho/CMSSW_5_3_20/src/RecoParticleFlow/PFProducer/python/particleFlow_cff.py                 
diff RecoJets/JetProducers/python/ak5PFJets_cfi.py			    /home/gobinda/analho/CMSSW_5_3_20/src/RecoJets/JetProducers/python/ak5PFJets_cfi.py                          
diff RecoMET/Configuration/python/RecoPFMET_cff.py			    /home/gobinda/analho/CMSSW_5_3_20/src/RecoMET/Configuration/python/RecoPFMET_cff.py                          
diff RecoParticleFlow/PFProducer/python/pfLinker_cff.py			    /home/gobinda/analho/CMSSW_5_3_20/src/RecoParticleFlow/PFProducer/python/pfLinker_cff.py                     
diff RecoParticleFlow/Configuration/python/RecoParticleFlow_cff.py	    /home/gobinda/analho/CMSSW_5_3_20/src/RecoParticleFlow/Configuration/python/RecoParticleFlow_cff.py 






*/

// -*- C++ -*-
//
// Package:    BarrelGenJetFilter
// Class:      BarrelGenJetFilter
// 
/**\class BarrelGenJetFilter BarrelGenJetFilter.cc Test/BarrelGenJetFilter/src/BarrelGenJetFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Wed Nov  5 19:04:47 IST 2014
// $Id$
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
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

//
// class declaration
//

class BarrelGenJetFilter : public edm::EDFilter {
   public:
      explicit BarrelGenJetFilter(const edm::ParameterSet&);
      ~BarrelGenJetFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

  // ----------member data ---------------------------
  edm::InputTag genColl_;
  double genPtCut_;
  double genEtaCut_;

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
BarrelGenJetFilter::BarrelGenJetFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  genColl_       = iConfig.getUntrackedParameter<edm::InputTag>("GenJetColl");
  genPtCut_      = iConfig.getParameter<double>("GenPtCut");
  genEtaCut_     = iConfig.getParameter<double>("GenEtaCut");
}


BarrelGenJetFilter::~BarrelGenJetFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
BarrelGenJetFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace reco;
  
  Handle<reco::GenJetCollection> genjets;
  iEvent.getByLabel(genColl_, genjets);
  if (genjets.isValid()) {
    for(GenJetCollection::const_iterator jet = genjets->begin(); jet != genjets->end(); ++jet) {
      if (jet->getGenConstituents().size() < 2) continue;
      if (jet->pt() < genPtCut_ ) continue;
      if (fabs(jet->eta()) > genEtaCut_ ) continue;
      return true;
    }
  }
  return false;
}

// ------------ method called once each job just before starting event loop  ------------
void 
BarrelGenJetFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
BarrelGenJetFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
BarrelGenJetFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
BarrelGenJetFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
BarrelGenJetFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
BarrelGenJetFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BarrelGenJetFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(BarrelGenJetFilter);
