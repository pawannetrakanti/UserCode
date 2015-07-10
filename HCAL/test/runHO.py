import FWCore.ParameterSet.VarParsing as VarParsing

ivars = VarParsing.VarParsing('python')
#ivars.inputFiles='file:./step3_0E886EA8-798E-E411-8CE9-0025905B85F6.root'
#ivars.inputFiles='/store/user/pawan/Phys14dr/step3/set0/step3_0E56F5ED-848E-E411-979F-0025905A611E.root'
ivars.inputFiles='file:./Output_step3_Phys14dr_crab_75.root'
ivars.outputFile = './Output_hoanalyzer_Phys14dr.root'

ivars.parseArguments()

import FWCore.ParameterSet.Config as cms  

process = cms.Process ("hoana") 
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')

process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.Simulation_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('CommonTools.ParticleFlow.EITopPAG_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.load("RecoLocalTracker.Configuration.RecoLocalTracker_cff")
process.load("CondCore.DBCommon.CondDBCommon_cfi")

process.load('FWCore.MessageService.MessageLogger_cfi')
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


#process.load('FWCore.MessageService.MessageLogger_cfi')
##process.MessageLogger.cerr.INFO.limit =  cms.untracked.int32(0)
#process.MessageLogger.cerr.FwkReport.reportEvery =  cms.untracked.int32(100)
#process.MessageLogger.cerr.threshold = cms.untracked.string('ERROR')
##process.MessageLogger.suppressInfo = cms.untracked.vstring('MeasurementTrackerEventProducer:MeasurementTrackerEvent',
##                                                              'CSCHaloDataProducer:CSCHaloData',
##                                                              'SeedGeneratorFromRegionHitsEDProducer:pixelPairStepSeeds'
##                                                           )
##process.MessageLogger.suppressWarning = cms.untracked.vstring('MeasurementTrackerEventProducer:MeasurementTrackerEvent',
##                                                              'CSCHaloDataProducer:CSCHaloData',
##                                                              'SeedGeneratorFromRegionHitsEDProducer:pixelPairStepSeeds'
##                                                              )
##
#
process.source = cms.Source("PoolSource", 
     fileNames = cms.untracked.vstring(ivars.inputFiles) 
) 
process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string(ivars.outputFile) 
) 
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1),
#    output = cms.untracked.int32(20)
)

#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck", 
#  ignoreTotal = cms.untracked.int32(5) ## default is one 
#) 

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc_GRun', '')

process.load('RecoJets.JetProducers.ak4PFJets_cfi')
process.load('RecoJets.JetProducers.ak4CaloJets_cfi')
process.load('RecoMET.Configuration.RecoTCMET_cff')
process.load('RecoMET.Configuration.RecoMET_cff')

# These are the original threshold used for HO
process.particleFlowRecHitHO.threshold_ring0 = cms.double(0.05)
process.particleFlowRecHitHO.threshold_rig12 = cms.double(0.08)

#process.dump = cms.EDAnalyzer("EventContentAnalyzer")

# Original
process.caloTowerForTrk.hbheInput = cms.InputTag("hbheprereco")

#process.localReReco = cms.Sequence(process.trackerlocalreco
#                                   +process.muonlocalreco
#                                   +process.calolocalreco
#                                   )
process.localReReco = cms.Sequence(process.trackerlocalreco)

process.load('RecoJets.JetProducers.ak4JetID_cfi')

process.globalReReco = cms.Sequence(process.offlineBeamSpot
                                    * process.MeasurementTrackerEventPreSplitting  # unclear where to put this
                                    * process.siPixelClusterShapeCachePreSplitting # unclear where to put this
                                    * process.standalonemuontracking
                                    * process.trackingGlobalReco
                                    * process.vertexreco                         
                                    * process.hcalGlobalRecoSequence
                                    * process.particleFlowCluster
                                    * process.ecalClusters
                                    * process.caloTowersRec
                                    * process.egammaGlobalReco
                                    * process.jetGlobalReco
                                    * process.muonGlobalReco
                                    * process.pfTrackingGlobalReco
                                    * process.muoncosmicreco
                                    * process.metreco
)


# Particle Flow re-processing                            
process.pfReReco = cms.Sequence(process.egammaHighLevelRecoPrePF
                                * process.particleFlowReco
                                * process.egammaHighLevelRecoPostPF
                                * process.regionalCosmicTracksSeq
                                * process.muoncosmichighlevelreco
                                * process.muonshighlevelreco 
                                * process.particleFlowLinks
                                * process.jetHighLevelReco
#                                * process.recoPFJets
                                * process.recoPFMET
)
#


## Gen Info re-processing
process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
process.load("RecoJets.Configuration.GenJetParticles_cff")
process.load("RecoJets.Configuration.RecoGenJets_cff")
process.load("RecoMET.Configuration.GenMETParticles_cff")
process.load("RecoMET.Configuration.RecoGenMET_cff")
process.load("RecoParticleFlow.PFProducer.particleFlowSimParticle_cff")
process.load("RecoParticleFlow.Configuration.HepMCCopy_cfi")

process.genMetCaloAndNonPrompt.src = cms.InputTag("genParticlesForJetsNoNu")

process.genReReco = cms.Sequence(process.genParticles
                                 * process.genParticlesForJetsNoNu
                                 * process.ak4GenJetsNoNu
                                 * process.genMETParticles
                                 * process.recoGenMET
)

process.load('HOTest.HOinPFAlgo.ho_00_cfi')
process.load('HOTest.HOinPFAlgo.ho_01_cfi')
process.load('HOTest.HOinPFAlgo.ho_02_cfi')
#
process.load('HOTest.HOinPFAlgo.hoinpfalgo_cfi')


process.primaryVertexFilter = cms.EDFilter(
    "VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
    filter = cms.bool(True)
)


#Filter to keep event in PoolOutputModule (first MET fileter>300 GeV,
#The HO tower energy >40 GeV
#process.metstorefilter = cms.EDFilter("HighMetFilter")
process.load('RecoMET.METFilters.metFilters_cff')
process.trackingFailureFilter.JetSource = cms.InputTag("ak4PFJets")
process.metFilters_sel = cms.Sequence( process.HBHENoiseFilter
#                                       * process.hcalLaserEventFilter
#                                       * process.EcalDeadCellTriggerPrimitiveFilter 
                                       * process.goodVertices 
                                       * process.trackingFailureFilter
#                                       * process.eeBadScFilter
#                                       * process.ecalLaserCorrFilter 
                                       * process.trkPOGFilters
)

process.barrelGenJets = cms.EDFilter("BarrelGenJetFilter",
                                     GenJetColl= cms.untracked.InputTag("ak4GenJetsNoNu"),
                                     GenPtCut  = cms.double(500.0),
                                     GenEtaCut = cms.double(2.0),
)
process.hoHitFilter = cms.EDFilter("HOEventFilter",
                                   hoInput    = cms.untracked.InputTag("horeco"),
                                   hoEnThresh = cms.double(0.0),
                                   CaloJetTag = cms.untracked.InputTag("ak4CaloJet")
)


process.allfilters = cms.Sequence(process.primaryVertexFilter
                                  * process.barrelGenJets
#                                  * process.hoHitFilter
#                                  * process.metFilters_sel
)

# PFClusterAnalyzerHO
process.pfClusterAnalyzerHO = cms.EDAnalyzer("PFClusterAnalyzer",
                                             PFClusters = cms.InputTag("particleFlowClusterHO"),
                                             verbose = cms.untracked.bool(True),
                                             printBlocks = cms.untracked.bool(False)
)
# PFClusterAnalyzerHO00
process.pfClusterAnalyzerHO00 = cms.EDAnalyzer("PFClusterAnalyzer",
                                             PFClusters = cms.InputTag("particleFlowClusterHO00"),
                                             verbose = cms.untracked.bool(True),
                                             printBlocks = cms.untracked.bool(False)
)


process.pfClusterCompare00 = cms.EDAnalyzer("PFClusterComparator",
                                             PFClusters = cms.InputTag("particleFlowClusterHO"),
                                             PFClustersCompare = cms.InputTag("particleFlowClusterHO00"),
                                             verbose = cms.untracked.bool(True)
)

process.pfcCheck00 = cms.EDAnalyzer("PFCandidateChecker",
                                    pfCandidatesReco = cms.InputTag("particleFlow","","ReRECO"),
                                    pfCandidatesReReco = cms.InputTag("particleFlow00","",""),
                                    pfJetsReco = cms.InputTag("ak4PFJets","","ReRECO"),
                                    pfJetsReReco = cms.InputTag("ak4PFJets00","",""),
                                    deltaEMax = cms.double(1E-5),
                                    deltaEtaMax = cms.double(1E-5),
                                    deltaPhiMax = cms.double(1E-5),
                                    printBlocks = cms.untracked.bool(False),
                                    rankByPt = cms.untracked.bool(True)
)

process.raw2digi_step = cms.Path(process.RawToDigi)
process.p = cms.Path(process.genReReco
#                     * process.hoinpf
#                     * process.allfilters
                     * (process.localReReco    * process.globalReReco    * process.pfReReco    * process.hoinpf   )    
                     * (process.localReReco00  * process.globalReReco00  * process.pfReReco00  * process.hoinpf00 )  
                     * (process.localReReco01  * process.globalReReco01  * process.pfReReco01  * process.hoinpf01 )  
                     * (process.localReReco02  * process.globalReReco02  * process.pfReReco02  * process.hoinpf02 )  
#                    * (process.localReReco00  + process.globalReReco00  + process.pfReReco00  + process.hoinpf00 + process.pfClusterCompare00 + process.pfcCheck00 )  
)

process.schedule = cms.Schedule(process.raw2digi_step,process.p)
##process.e = cms.EndPath(process.out) 
#process.p = cms.Path(process.patDefaultSequence*process.kt6PFJets*process.ak5PFJetsL1FastL2L3*process.hotower) 
#process.out.SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p'))


