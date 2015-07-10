
import FWCore.ParameterSet.Config as cms  

process = cms.Process ("Combined") 
process.load('Configuration/StandardSequences/Reconstruction_cff') 
#process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")

process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
#Tracker 
process.load('RecoTracker.GeometryESProducer.TrackerRecoGeometryESProducer_cfi')
process.load('TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi') 
process.load('TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff')
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff') 
#process.GlobalTag.globaltag = 'START42_V12::All'
process.GlobalTag.globaltag = 'GR_R_53_V21::All'
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('RecoJets.JetProducers.ak5PFJets_cfi')
process.load('RecoJets.JetProducers.ak5CaloJets_cfi')
process.load('RecoMET.Configuration.RecoTCMET_cff')

#GMA CHECK THIS
process.particleFlowRecHitHO.thresh_Barrel=0.4
process.particleFlowRecHitHO.thresh_Endcap=1.0

process.particleFlowClusterHO.thresh_Seed_Barrel=0.5
process.particleFlowClusterHO.thresh_Barrel=0.5
process.particleFlowClusterHO.thresh_Seed_Endcap=1.0
process.particleFlowClusterHO.thresh_Endcap=1.0

process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(-1)
) 
process.MessageLogger = cms.Service('MessageLogger', 
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

process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(
#'/store/relval/CMSSW_5_0_0_pre3/RelValQCD_Pt_3000_3500/GEN-SIM-RECO/START50_V3-v1/0010/82AC7FB1-C2FA-E011-9A33-0026189438DB.root'
#'file:/localdata/gobinda/anal/mc/singleneutron500_reco.root'
#'/store/data/Run2012A/Jet/RECO/PromptReco-v1/000/191/701/AE2DE611-818B-E111-A42D-0025901D6272.root'
#'/store/data/Run2012A/MET/RECO/PromptReco-v1/000/193/336/0E2E715B-9597-E111-828C-BCAEC53296F7.root'
#'/store/data/Run2012A/ElectronHad/RECO/PromptReco-v1/000/193/336/F6E9EA47-9F97-E111-8B15-001D09F2B2CF.root'
#'/store/data/Run2012A/MuHad/RECO/PromptReco-v1/000/193/336/F84DF932-A497-E111-82C5-001D09F2447F.root'
#'file:/localdata/gobinda/data/jet_2012a_promptv1_AE2DE611-818B-E111-A42D-0025901D6272.root'
#'/store/data/Run2012A/SingleMu/RECO/PromptReco-v1/000/192/257/4AFB2FD4-BC8E-E111-BAF5-003048F117B4.root'
'/store/data/Run2012A/Jet/RAW-RECO/20Nov2012-v2/00000/EEE23F90-AC37-E211-8923-001BFCDBD15E.root'
    )
) 

process.SimpleMemoryCheck = cms.Service('SimpleMemoryCheck',
  ignoreTotal = cms.untracked.int32(1) ## default is one 
) 

process.TFileService = cms.Service("TFileService",
# fileName = cms.string('hist_met2012a_mod.root'),
# fileName = cms.string('hist_jet2012a_mod.root'),
# fileName = cms.string('hist_ht2012a_mod.root'),
# fileName = cms.string('hist_jetmon2012b_mod.root'),
 fileName = cms.string('hist_jetht2012b_mod.root'),
)

process.hoinpf = cms.EDAnalyzer("HOinPFAlgo",
  MonteCarlo =  cms.untracked.bool(False),
  HOTagged = cms.untracked.int32(1),
  OnlyDIGI =  cms.untracked.bool(False),
  ReRECO =  cms.untracked.bool(True),                             
  Tagged = cms.untracked.int32(17),

#  RootFileName = cms.untracked.string('hoinpf_met2012a_mod.root'),
#  RootFileName = cms.untracked.string('hoinpf_jet2012a_mod.root'),
#  RootFileName = cms.untracked.string('hoinpf_ht2012a_mod.root'),
#  RootFileName = cms.untracked.string('hoinpf_jetmon2012b_mod.root'),
  RootFileName = cms.untracked.string('hoinpf_jetht2012b_mod.root'),
  HistFill = cms.untracked.bool(False)
)

process.dump = cms.EDAnalyzer("EventContentAnalyzer")
process.oout = cms.OutputModule("PoolOutputModule",
   outputCommands = cms.untracked.vstring('keep *'),
                                
#    fileName = cms.untracked.string('event_met2012a_mod.root'),
#    fileName = cms.untracked.string('event_jet2012a_mod.root'),
#    fileName = cms.untracked.string('event_ht2012a_mod.root'),
#    fileName = cms.untracked.string('event_jetmon2012b_mod.root'),
    fileName = cms.untracked.string('event_jetht2012b_mod.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p1')
  )
)

# Local re-reco: Produce tracker rechits, pf rechits and pf clusters

process.localReReco = cms.Sequence(process.siPixelRecHits+
                                   process.siStripMatchedRecHits+
                                   #process.hbhereflag+
##                                   process.particleFlowCluster+
                                   process.particleFlowClusterWithoutHO+
                                   process.pfClusteringHO+
                                   process.ecalClusters)

process.localReReco00 = cms.Sequence(process.pfClusterHO00)

process.localReReco01 = cms.Sequence(process.pfClusterHO01)

process.localReReco02 = cms.Sequence(process.pfClusterHO02)

process.localReReco03 = cms.Sequence(process.pfClusterHO03)

process.localReReco04 = cms.Sequence(process.pfClusterHO04)

process.localReReco05 = cms.Sequence(process.pfClusterHO05)

process.localReReco06 = cms.Sequence(process.pfClusterHO06)

process.localReReco07 = cms.Sequence(process.pfClusterHO07)

process.localReReco08 = cms.Sequence(process.pfClusterHO08)

process.localReReco09 = cms.Sequence(process.pfClusterHO09)

process.localReReco10 = cms.Sequence(process.pfClusterHO10)
process.localReReco11 = cms.Sequence(process.pfClusterHO11)
process.localReReco12 = cms.Sequence(process.pfClusterHO12)
process.localReReco13 = cms.Sequence(process.pfClusterHO13)
process.localReReco14 = cms.Sequence(process.pfClusterHO14)
process.localReReco15 = cms.Sequence(process.pfClusterHO15)
process.localReReco16 = cms.Sequence(process.pfClusterHO16)


# Track re-reco
process.globalReReco =  cms.Sequence(
#    process.offlineBeamSpot+
#    process.recopixelvertexing+
    process.ckftracks+
    process.caloTowersRec+
#    process.vertexreco+
#    process.recoJets+
    process.ak5CaloJets+
    process.muonrecoComplete+
    process.muoncosmicreco+
    process.egammaGlobalReco+
    process.pfTrackingGlobalReco+
    process.egammaHighLevelRecoPrePF+
    process.muoncosmichighlevelreco+
    process.metreco
#    process.met+process.tcMet
    )

process.globalReReco00 =  cms.Sequence(
    process.towerMakerWithHO00+
    process.ak5CaloJets00+
    process.metHO00+
#    process.muonrecoComplete+
#    process.muoncosmicreco+
#    process.egammaGlobalReco+
#    process.pfTrackingGlobalReco+
#    process.egammaHighLevelRecoPrePF+
#    process.muoncosmichighlevelreco+
    process.muonMETValueMapProducer00+
    process.corMetGlobalMuons00+
    process.muonTCMETValueMapProducer00+
    process.tcMet00
    )

process.globalReReco01 =  cms.Sequence(
    process.towerMakerWithHO01+
    process.ak5CaloJets01+
    process.metHO01+
    process.muonMETValueMapProducer01+
    process.corMetGlobalMuons01+
    process.muonTCMETValueMapProducer01+   
    process.tcMet01
    )

#process.globalReReco02 =  cms.Sequence(
#    process.towerMakerWithHO02+
#    process.ak5CaloJets02+
#    process.metHO02+
#    process.muonMETValueMapProducer02+
#    process.corMetGlobalMuons02+
#    process.muonTCMETValueMapProducer02+
#    process.tcMet02
#    )

process.globalReReco02 =  cms.Sequence(
    process.towerMakerWithHO02+
    process.ak5CaloJets02+
    process.metHO02+
    process.muonMETValueMapProducer02+
    process.corMetGlobalMuons02+
    process.muonTCMETValueMapProducer02+
    process.tcMet02
    )

process.globalReReco03 =  cms.Sequence(
    process.towerMakerWithHO03+
    process.ak5CaloJets03+
    process.metHO03+
    process.muonMETValueMapProducer03+
    process.corMetGlobalMuons03+
    process.muonTCMETValueMapProducer03+
    process.tcMet03
    )

process.globalReReco04 =  cms.Sequence(
    process.towerMakerWithHO04+
    process.ak5CaloJets04+
    process.metHO04+
    process.muonMETValueMapProducer04+
    process.corMetGlobalMuons04+
    process.muonTCMETValueMapProducer04+
    process.tcMet04
    )

process.globalReReco05 =  cms.Sequence(
    process.towerMakerWithHO05+
    process.ak5CaloJets05+
    process.metHO05+
    process.muonMETValueMapProducer05+
    process.corMetGlobalMuons05+
    process.muonTCMETValueMapProducer05+
    process.tcMet05
    )

process.globalReReco06 =  cms.Sequence(
    process.towerMakerWithHO06+
    process.ak5CaloJets06+
    process.metHO06+
    process.muonMETValueMapProducer06+
    process.corMetGlobalMuons06+
    process.muonTCMETValueMapProducer06+
    process.tcMet06
    )

process.globalReReco07 =  cms.Sequence(
    process.towerMakerWithHO07+
    process.ak5CaloJets07+
    process.metHO07+
    process.muonMETValueMapProducer07+
    process.corMetGlobalMuons07+
    process.muonTCMETValueMapProducer07+
    process.tcMet07
    )

process.globalReReco08 =  cms.Sequence(
    process.towerMakerWithHO08+
    process.ak5CaloJets08+
    process.metHO08+
    process.muonMETValueMapProducer08+
    process.corMetGlobalMuons08+
    process.muonTCMETValueMapProducer08+
    process.tcMet08
    )

process.globalReReco09 =  cms.Sequence(
    process.towerMakerWithHO09+
    process.ak5CaloJets09+
    process.metHO09+
    process.muonMETValueMapProducer09+
    process.corMetGlobalMuons09+
    process.muonTCMETValueMapProducer09+
    process.tcMet09
    )

process.globalReReco10 =  cms.Sequence(
    process.towerMakerWithHO10+
    process.ak5CaloJets10+
    process.metHO10+
    process.muonMETValueMapProducer10+
    process.corMetGlobalMuons10+
    process.muonTCMETValueMapProducer10+
    process.tcMet10
    )

process.globalReReco11 =  cms.Sequence(
    process.towerMakerWithHO11+
    process.ak5CaloJets11+
    process.metHO11+
    process.muonMETValueMapProducer11+
    process.corMetGlobalMuons11+
    process.muonTCMETValueMapProducer11+
    process.tcMet11
    )

process.globalReReco12 =  cms.Sequence(
    process.towerMakerWithHO12+
    process.ak5CaloJets12+
    process.metHO12+
    process.muonMETValueMapProducer12+
    process.corMetGlobalMuons12+
    process.muonTCMETValueMapProducer12+
    process.tcMet12
    )

process.globalReReco13 =  cms.Sequence(
    process.towerMakerWithHO13+
    process.ak5CaloJets13+
    process.metHO13+
    process.muonMETValueMapProducer13+
    process.corMetGlobalMuons13+
    process.muonTCMETValueMapProducer13+
    process.tcMet13
    )

process.globalReReco14 =  cms.Sequence(
    process.towerMakerWithHO14+
    process.ak5CaloJets14+
    process.metHO14+
    process.muonMETValueMapProducer14+
    process.corMetGlobalMuons14+
    process.muonTCMETValueMapProducer14+
    process.tcMet14
    )

process.globalReReco15 =  cms.Sequence(
    process.towerMakerWithHO15+
    process.ak5CaloJets15+
    process.metHO15+
    process.muonMETValueMapProducer15+
    process.corMetGlobalMuons15+
    process.muonTCMETValueMapProducer15+
    process.tcMet15
    )

process.globalReReco16 =  cms.Sequence(
    process.towerMakerWithHO16+
    process.ak5CaloJets16+
    process.metHO16+
    process.muonMETValueMapProducer16+
    process.corMetGlobalMuons16+
    process.muonTCMETValueMapProducer16+
    process.tcMet16
    )


# Particle Flow re-processing
process.pfReReco = cms.Sequence(process.particleFlowReco+
                                process.egammaHighLevelRecoPostPF+
                                process.muonshighlevelreco+
                                process.particleFlowLinks+
                                process.recoPFJets+
                                process.recoPFMET)
#                                process.PFTau)
 
# Gen Info re-processing
process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
process.load("RecoJets.Configuration.GenJetParticles_cff")
process.load("RecoJets.Configuration.RecoGenJets_cff")
process.load("RecoMET.Configuration.GenMETParticles_cff")
process.load("RecoMET.Configuration.RecoGenMET_cff")
process.load("RecoParticleFlow.PFProducer.particleFlowSimParticle_cff")
process.load("RecoParticleFlow.Configuration.HepMCCopy_cfi")

#process.genReReco = cms.Sequence(process.generator+
#                                 process.genParticles+
#                                 process.genJetParticles+
#                                 process.recoGenJets+process.genParticlesForJetsNoNu*process.ak5GenJetsNoNu+
#                                 process.genMETParticles+
#                                 process.recoGenMET+
#                                 process.particleFlowSimParticle)

process.genNoNuJet =cms.Sequence(process.genParticles*
                                 process.genJetParticles*
                                 process.genParticlesForJetsNoNu*
                                 process.ak5GenJetsNoNu)

process.pfReReco00 = cms.Sequence(process.particleFlowReco00+
                                process.particleFlowLinks00+
                                process.ak5PFJets00+
                                process.recoPFMET00)

process.pfReReco01 = cms.Sequence(process.particleFlowReco01+
                                 process.particleFlowLinks01+
                                 process.ak5PFJets01+
                                 process.recoPFMET01)

process.pfReReco02 = cms.Sequence(process.particleFlowReco02+
                                 process.particleFlowLinks02+
                                 process.ak5PFJets02+
                                 process.recoPFMET02)

process.pfReReco03 = cms.Sequence(process.particleFlowReco03+
                                 process.particleFlowLinks03+
                                 process.ak5PFJets03+
                                 process.recoPFMET03)


process.pfReReco04 = cms.Sequence(process.particleFlowReco04+
                                 process.particleFlowLinks04+
                                 process.ak5PFJets04+
                                 process.recoPFMET04)

process.pfReReco05 = cms.Sequence(process.particleFlowReco05+
                                 process.particleFlowLinks05+
                                 process.ak5PFJets05+
                                 process.recoPFMET05)

process.pfReReco06 = cms.Sequence(process.particleFlowReco06+
                                 process.particleFlowLinks06+
                                 process.ak5PFJets06+
                                 process.recoPFMET06)

process.pfReReco07 = cms.Sequence(process.particleFlowReco07+
                                 process.particleFlowLinks07+
                                 process.ak5PFJets07+
                                 process.recoPFMET07)

process.pfReReco08 = cms.Sequence(process.particleFlowReco08+
                                 process.particleFlowLinks08+
                                 process.ak5PFJets08+
                                 process.recoPFMET08)

process.pfReReco09 = cms.Sequence(process.particleFlowReco09+
                                 process.particleFlowLinks09+
                                 process.ak5PFJets09+
                                 process.recoPFMET09)

process.pfReReco10 = cms.Sequence(process.particleFlowReco10+
                                 process.particleFlowLinks10+
                                 process.ak5PFJets10+
                                 process.recoPFMET10)

process.pfReReco11 = cms.Sequence(process.particleFlowReco11+
                                 process.particleFlowLinks11+
                                 process.ak5PFJets11+
                                 process.recoPFMET11)

process.pfReReco12 = cms.Sequence(process.particleFlowReco12+
                                 process.particleFlowLinks12+
                                 process.ak5PFJets12+
                                 process.recoPFMET12)

process.pfReReco13 = cms.Sequence(process.particleFlowReco13+
                                 process.particleFlowLinks13+
                                 process.ak5PFJets13+
                                 process.recoPFMET13)

process.pfReReco14 = cms.Sequence(process.particleFlowReco14+
                                 process.particleFlowLinks14+
                                 process.ak5PFJets14+
                                 process.recoPFMET14)

process.pfReReco15 = cms.Sequence(process.particleFlowReco15+
                                 process.particleFlowLinks15+
                                 process.ak5PFJets15+
                                 process.recoPFMET15)

process.pfReReco16 = cms.Sequence(process.particleFlowReco16+
                                 process.particleFlowLinks16+
                                 process.ak5PFJets16+
                                 process.recoPFMET16)


process.hojetfilter = cms.EDFilter("BarrelJetFilter",
  Ptcut = cms.untracked.double(99.0),
  Etacut = cms.untracked.double(1.1),
  HOcut = cms.untracked.double(5.0),
  OthHistFill = cms.untracked.bool(False)                                 
)

# Filter to keep event in PoolOutputModule (first MET fileter>300 GeB,
#The HO tower energy >40 GeV
process.metstorefilter = cms.EDFilter("HighMetFilter")

process.primaryVertexFilter = cms.EDFilter(
    "VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
    filter = cms.bool(True)
    )

process.allfilters = cms.Sequence(
   process.hojetfilter* 
   process.primaryVertexFilter 
)

process.p1 = cms.Path(
    process.allfilters*
    (process.localReReco+process.globalReReco+process.pfReReco)+
    (process.localReReco00+process.globalReReco00+process.pfReReco00)+
    (process.localReReco01+process.globalReReco01+process.pfReReco01)*
    (process.localReReco02+process.globalReReco02+process.pfReReco02)*
    (process.localReReco03+process.globalReReco03+process.pfReReco03)*
    (process.localReReco04+process.globalReReco04+process.pfReReco04)*
    (process.localReReco05+process.globalReReco05+process.pfReReco05)*
    (process.localReReco06+process.globalReReco06+process.pfReReco06)*
    (process.localReReco07+process.globalReReco07+process.pfReReco07)*
    (process.localReReco08+process.globalReReco08+process.pfReReco08)*
    (process.localReReco09+process.globalReReco09+process.pfReReco09)*
    (process.localReReco10+process.globalReReco10+process.pfReReco10)*
    (process.localReReco11+process.globalReReco11+process.pfReReco11)*
    (process.localReReco12+process.globalReReco12+process.pfReReco12)*
    (process.localReReco13+process.globalReReco13+process.pfReReco13)*
    (process.localReReco14+process.globalReReco14+process.pfReReco14)*
    (process.localReReco15+process.globalReReco15+process.pfReReco15)*
    (process.localReReco16+process.globalReReco16+process.pfReReco16)*
        
    (process.hoinpf*process.metstorefilter ))

#process.RECOSIMoutput_step = cms.EndPath(process.oout)
#process.schedule = cms.Schedule(process.p1,process.RECOSIMoutput_step)

#process.e = cms.EndPath(process.oout) 

#End of BarrelJetFilter with event 9223 Jetpassed 6389 passed 7248
#End of HOinPFAlgo with event 2796
#End of HighMetFilter 2796 43 1
