import FWCore.ParameterSet.VarParsing as VarParsing

ivars = VarParsing.VarParsing('python')
#ivars.inputFiles = 'file:00A8FBA5-41E1-E111-934D-0025901D4C3E.root'
ivars.inputFiles = 'file:./RECO/RECO_3000_3500_8TeV_DIGIRAWRECO_1.root'
ivars.outputFile = './Output_hoanalyzer.root'
ivars.parseArguments()


import FWCore.ParameterSet.Config as cms  

process = cms.Process ("hoana") 
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff") 

# Reconstruction geometry services 
#  Tracking Geometry 
#process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi") 

#Tracker 
#process.load("RecoTracker.GeometryESProducer.TrackerRecoGeometryESProducer_cfi") 
#process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi") 
#process.load("TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff") 


process.MessageLogger = cms.Service("MessageLogger", 
 cout = cms.untracked.PSet( 
  default = cms.untracked.PSet( ## kill all messages in the log 
   limit = cms.untracked.int32(0) 
  ), 
  FwkJob = cms.untracked.PSet( ## but FwkJob category - those unlimitted 
   limit = cms.untracked.int32(-1) 
  ) 
 ), 
 categories = cms.untracked.vstring('FwkJob'), 
 destinations = cms.untracked.vstring('cout') 
) 

process.source = cms.Source("PoolSource", 
     fileNames = cms.untracked.vstring(ivars.inputFiles) 
) 

#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck", 
#  ignoreTotal = cms.untracked.int32(1) ## default is one 
#) 


process.load('Configuration.StandardSequences.Reconstruction_cff') 
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

process.load('RecoJets.Configuration.GenJetParticles_cff')
from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
process.ak5GenJetsNoNu = ak5GenJets.clone(src = cms.InputTag("genParticlesForJetsNoNu"))

process.ak5GenJetsSequence = cms.Sequence(process.genParticlesForJetsNoNu*process.ak5GenJetsNoNu)

from RecoJets.Configuration.RecoGenJets_cff import *
from RecoJets.Configuration.RecoJets_cff import *
from RecoJets.Configuration.RecoPFJets_cff import *
process.kt6PFJets.doRhoFastjet = True 
process.ak5PFJets.doAreaFastjet = True 

process.load("Configuration.StandardSequences.Services_cff") 
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff") 
process.GlobalTag.globaltag = 'START53_V29B::All' 
#process.GlobalTag.globaltag = 'START53_V7A::All' 


from RecoJets.JetProducers.JetIDParams_cfi import * 
process.ak5CaloJetID = cms.EDProducer('JetIDProducer', JetIDParams, 
                                      src = cms.InputTag('ak5CaloJets') 
) 
from RecoJets.JetProducers.JetIDParams_cfi import * 
process.ak5PFJetID = cms.EDProducer('JetIDProducer', JetIDParams, 
                                    src = cms.InputTag('ak5CaloJets') 
) 


from PhysicsTools.PatAlgos.mcMatchLayer0.jetMatch_cfi import *
process.ak5CaloPartonMatch = patJetPartonMatch.clone(
    src         = cms.InputTag("ak5CaloJets"),   # RECO objects to match
    matched     = cms.InputTag("genParticles"),   # mc-truth particle collection
    maxDeltaR   = cms.double(0.4),                # Minimum deltaR for the match
    maxDPtRel   = cms.double(3.0),                # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),       # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(False),# False = just match input in order; True = pick lowest deltaR pair first
)
process.ak5CaloGenMatchJets = patJetGenJetMatch.clone(
    src = cms.InputTag("ak5CaloJets"),
    matched = cms.InputTag("ak5GenJetsNoNu"),
    maxDeltaR = cms.double(0.4),
    maxDPtRel = cms.double(3.0),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(False)
)
process.ak5PFPartonMatch = patJetPartonMatch.clone(
    src         = cms.InputTag("ak5PFJets"),   # RECO objects to match
    matched     = cms.InputTag("genParticles"),   # mc-truth particle collection
    maxDeltaR   = cms.double(0.4),                # Minimum deltaR for the match
    maxDPtRel   = cms.double(3.0),                # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),       # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(False),# False = just match input in order; True = pick lowest deltaR pair first
)
process.ak5PFGenMatchJets = patJetGenJetMatch.clone(
    src = cms.InputTag("ak5PFJets"),
    matched = cms.InputTag("ak5GenJetsNoNu"),
    maxDeltaR = cms.double(0.4),
    maxDPtRel = cms.double(3.0),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(False)
)

from PhysicsTools.PatAlgos.recoLayer0.jetCorrFactors_cfi import patJetCorrFactors
process.ak5CaloJetCorrFactors = patJetCorrFactors.clone(
    emf = cms.bool(False),
    src = cms.InputTag("ak5CaloJets"),
    payload = cms.string('AK5Calo'),
    levels = cms.vstring('L1FastJet', 'L2Relative', 'L3Absolute',#'L5Flavor', 'L7Parton'
                         ),
)
process.ak5PFJetCorrFactors = patJetCorrFactors.clone(
    emf = cms.bool(False),
    src = cms.InputTag("ak5PFJets"),
    payload = cms.string('AK5PF'),
    levels = cms.vstring('L1FastJet', 'L2Relative', 'L3Absolute',#'L5Flavor', 'L7Parton'
                         ),
)


from  PhysicsTools.PatAlgos.mcMatchLayer0.jetFlavourId_cff import *
process.ak5CaloJetPartons = patJetPartons.clone(
    src = cms.InputTag("genParticles")
)
process.ak5CaloJetPartonAssociation = patJetPartonAssociation.clone(
    jets    = cms.InputTag("ak5CaloJets"),
    partons = cms.InputTag("ak5CaloJetPartons"),
    coneSizeToAssociate = cms.double(0.3),
)
process.ak5CaloJetFlavourAssociation = patJetFlavourAssociation.clone(
    srcByReference = cms.InputTag("ak5CaloJetPartonAssociation"),
)
process.ak5CaloJetFlavourId = cms.Sequence(process.ak5CaloJetPartons 
                                           * process.ak5CaloJetPartonAssociation
                                           * process.ak5CaloJetFlavourAssociation
)
process.ak5PFJetPartons = patJetPartons.clone(
    src = cms.InputTag("genParticles")
)
process.ak5PFJetPartonAssociation = patJetPartonAssociation.clone(
    jets    = cms.InputTag("ak5PFJets"),
    partons = cms.InputTag("ak5PFJetPartons"),
    coneSizeToAssociate = cms.double(0.3),
)
process.ak5PFJetFlavourAssociation = patJetFlavourAssociation.clone(
    srcByReference = cms.InputTag("ak5PFJetPartonAssociation"),
)
process.ak5PFJetFlavourId = cms.Sequence(process.ak5PFJetPartons 
                                         * process.ak5PFJetPartonAssociation
                                         * process.ak5PFJetFlavourAssociation
)

from RecoJets.JetAssociationProducers.j2tParametersCALO_cfi import *
from RecoJets.JetAssociationProducers.j2tParametersVX_cfi import *
process.ak5JetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
                                                        j2tParametersVX,
                                                        jets = cms.InputTag("ak5CaloJets")
)
 
process.ak5JetTracksAssociatorAtVertexPF = cms.EDProducer("JetTracksAssociatorAtVertex",
                                                          j2tParametersVX,
                                                          jets = cms.InputTag("ak5PFJets")
)

from  PhysicsTools.PatAlgos.recoLayer0.jetTracksCharge_cff import patJetCharge
process.ak5CaloJetCharge = patJetCharge.clone(
    src = cms.InputTag("ak5JetTracksAssociatorAtVertex"), ## a reco::JetTracksAssociation::Container
)
from  PhysicsTools.PatAlgos.recoLayer0.jetTracksCharge_cff import patJetCharge
process.ak5PFJetCharge = patJetCharge.clone(
    src = cms.InputTag("ak5JetTracksAssociatorAtVertexPF"), ## a reco::JetTracksAssociation::Container
)


from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import patJets
process.ak5CalopatJets = patJets.clone(jetSource = cms.InputTag("ak5CaloJets"),
                                       addGenPartonMatch   = cms.bool(True),
                                       embedGenPartonMatch = cms.bool(True),          
                                       embedGenJetMatch    = cms.bool(True),
                                       addGenJetMatch      = cms.bool(True),
                                       genPartonMatch = cms.InputTag("ak5CaloPartonMatch"),
                                       
                                       genJetMatch = cms.InputTag("ak5CaloGenMatchJets"),
                                       
                                       jetIDMap = cms.InputTag("ak5CaloJetID"),

                                       addJetCorrFactors    = cms.bool(True),
                                       jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak5CaloJetCorrFactors") ),
                                       
                                       getJetMCFlavour    = cms.bool(True),
                                       JetPartonMapSource = cms.InputTag("ak5CaloJetFlavourAssociation"),
                                       
                                       addAssociatedTracks    = cms.bool(True),
                                       trackAssociationSource = cms.InputTag("ak5JetTracksAssociatorAtVertex"),
                                       
                                       addJetCharge    = cms.bool(True),
                                       jetChargeSource = cms.InputTag("ak5CaloJetCharge"),
                                       
                                       # btag information
                                       addBTagInfo          = cms.bool(False),   ## master switch
                                       addDiscriminators    = cms.bool(False),   ## addition btag discriminator

                                       # embedding of RECO items (do not use on AOD input!)
                                       embedCaloTowers = cms.bool(True)
)


process.ak5PFpatJets = patJets.clone(jetSource = cms.InputTag("ak5PFJets"),
                                     addGenPartonMatch   = cms.bool(True),
                                     embedGenPartonMatch = cms.bool(True),          
                                     embedGenJetMatch    = cms.bool(True),
                                     addGenJetMatch      = cms.bool(True),
                                     genPartonMatch = cms.InputTag("ak5PFPartonMatch"),
                                     
                                     genJetMatch = cms.InputTag("ak5PFGenMatchJets"),
                                     
                                     jetIDMap = cms.InputTag("ak5JetID"),
                                     
                                     addJetCorrFactors    = cms.bool(True),
                                     jetCorrFactorsSource = cms.VInputTag(cms.InputTag("ak5PFJetCorrFactors") ),
                                     
                                     getJetMCFlavour    = cms.bool(True),
                                     JetPartonMapSource = cms.InputTag("ak5PFJetFlavourAssociation"),
                                     
                                     addAssociatedTracks    = cms.bool(True),
                                     trackAssociationSource = cms.InputTag("ak5JetTracksAssociatorAtVertexPF"),
                                       
                                     addJetCharge    = cms.bool(True),
                                     jetChargeSource = cms.InputTag("ak5PFJetCharge"),
                                     
                                     # btag information
                                     addBTagInfo          = cms.bool(False),   ## master switch
                                     addDiscriminators    = cms.bool(False),   ## addition btag discriminator

                                     embedPFCandidates = cms.bool(True)
)


process.ak5CalopatJetSequence = cms.Sequence(process.ak5CaloJetCorrFactors
                                             * process.ak5JetTracksAssociatorAtVertex
                                             * process.ak5CaloJetCharge
                                             * process.ak5CaloPartonMatch
                                             * process.ak5CaloGenMatchJets
                                             * process.ak5CaloJetFlavourId
                                             * process.ak5CaloJetID
                                             * process.ak5CalopatJets
)                                             

process.ak5PFpatJetSequence = cms.Sequence(process.ak5PFJetCorrFactors
                                           * process.ak5JetTracksAssociatorAtVertexPF
                                           * process.ak5PFJetCharge
                                           * process.ak5PFPartonMatch
                                           * process.ak5PFGenMatchJets
                                           * process.ak5PFJetFlavourId
                                           * process.ak5PFJetID
                                           * process.ak5PFpatJets
)                                             


process.hotower = cms.EDAnalyzer("HOAnalyzer",
                                 HistFill = cms.untracked.bool(True),
                                 Trigger = cms.untracked.bool(False),
                                 RECO = cms.untracked.bool(True),
                                 MonteCarlo =  cms.untracked.bool(True), 
                                 towerInput = cms.InputTag('towerMaker'),
                                 jetTags = cms.VInputTag(cms.InputTag('ak5CalopatJets'),
                                                         cms.InputTag('ak5PFpatJets')
                                                         )
)
 
process.TFileService = cms.Service("TFileService", 
#                                   fileName = cms.string('Output_hoanalyzer_whistos.root'), 
                                   fileName = cms.string(ivars.outputFile), 
) 

process.p = cms.Path(process.kt6PFJets
                     * process.ak5GenJetsSequence
                     * process.ak5CalopatJetSequence
                     * process.ak5PFpatJetSequence
                     * process.hotower
)
#process.e = cms.EndPath(process.o) 
#process.p = cms.Path(process.patDefaultSequence*process.kt6PFJets*process.ak5PFJetsL1FastL2L3*process.hotower) 
#process.out.SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p'))

