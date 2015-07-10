import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Reconstruction_cff import *
from RecoLocalCalo.CaloTowersCreator.calotowermaker_cfi import *
from RecoJets.Configuration.CaloTowersRec_cff import *

from RecoParticleFlow.PFProducer.pfLinker_cfi import *

from RecoJets.JetProducers.ak4PFJets_cfi import *
from RecoJets.JetProducers.ak4CaloJets_cfi import *
from RecoMET.Configuration.RecoTCMET_cff import *
from RecoMET.METProducers.CaloMET_cfi import *
from RecoMET.METProducers.MetMuonCorrections_cff import *
from RecoMET.METProducers.MuonTCMETValueMapProducer_cff import *

from RecoMET.METProducers.PFMET_cfi import *
from RecoMET.METProducers.pfChMet_cfi import *

_useho        = cms.bool(True)
_tagho        = cms.int32(0) # if _useho is False _tagho is redundant

_threshold0  = cms.double(0.08)
_threshold12 = cms.double(0.08)

_seedingThreshold0     = cms.double(0.08)
_gatheringThreshold0   = cms.double(0.05)
_rechitEnergyNorm0     = cms.double(0.05)

_seedingThreshold12   = cms.double(0.08)
_gatheringThreshold12 = cms.double(0.05)
_rechitEnergyNorm12   = cms.double(0.05)

_pcallogWeightDenominator = cms.double(0.05)
_allcellpcalclogWeightDenominator = cms.double(0.05)


#Default
# threshold_ring0 = cms.double(0.4),
# threshold_rig12 = cms.double(1.0)
particleFlowRecHitHO00 = particleFlowRecHitHO.clone(
    threshold_ring0 = _threshold0,
    threshold_rig12 = _threshold12
)

#Default
# RING0   seedingThreshold = cms.double(1.0),
# RING1   seedingThreshold = cms.double(3.1),
from RecoParticleFlow.PFClusterProducer.particleFlowClusterHO_cfi import *
particleFlowClusterHO00 = particleFlowClusterHO.clone(
    seedFinder = cms.PSet(
        algoName = cms.string("LocalMaximumSeedFinder"),    
        thresholdsByDetector = cms.VPSet(                   
            cms.PSet( detector = cms.string("HCAL_BARREL2_RING0"),
                      seedingThreshold = _seedingThreshold0,
                      seedingThresholdPt = cms.double(0.0)
                      ),                             
            cms.PSet( detector = cms.string("HCAL_BARREL2_RING1"), 
                      seedingThreshold = _seedingThreshold12,
                      seedingThresholdPt = cms.double(0.0)
                      ) 
            ),
        nNeighbours = cms.int32(4) 
        ),

    initialClusteringStep = cms.PSet(
        algoName = cms.string("Basic2DGenericTopoClusterizer"),
        thresholdsByDetector = cms.VPSet(
            cms.PSet( detector = cms.string("HCAL_BARREL2_RING0"),
                      gatheringThreshold = _gatheringThreshold0,
                      gatheringThresholdPt = cms.double(0.0)
                      ),
            cms.PSet( detector = cms.string("HCAL_BARREL2_RING1"),
                      gatheringThreshold   = _gatheringThreshold12,
                      gatheringThresholdPt = cms.double(0.0)
                      )
            ),
        useCornerCells = cms.bool(True)
        ),

    pfClusterBuilder = cms.PSet(
        algoName = cms.string("Basic2DGenericPFlowClusterizer"),
        #pf clustering parameters
        minFractionToKeep = cms.double(1e-7),
        positionCalc = cms.PSet(  
            algoName = cms.string("Basic2DGenericPFlowPositionCalc"),    ## 
            minFractionInCalc = cms.double(1e-9),
            posCalcNCrystals = cms.int32(5),
            logWeightDenominator    = _pcallogWeightDenominator, # same as gathering threshold
            minAllowedNormalization = cms.double(1e-9)
            ),
        allCellsPositionCalc = cms.PSet(  
            algoName = cms.string("Basic2DGenericPFlowPositionCalc"),    ## 
            minFractionInCalc = cms.double(1e-9),
            posCalcNCrystals = cms.int32(-1),
            logWeightDenominator    = _allcellpcalclogWeightDenominator, # same as gathering threshold
            minAllowedNormalization = cms.double(1e-9)
            ),
        showerSigma = cms.double(00.0),
        stoppingTolerance = cms.double(1e-8),
        maxIterations = cms.uint32(50),
        excludeOtherSeeds = cms.bool(True),
        minFracTot = cms.double(1e-20), ## numerical stabilization
        recHitEnergyNorms = cms.VPSet(
            cms.PSet( detector = cms.string("HCAL_BARREL2_RING0"),
                      recHitEnergyNorm = _rechitEnergyNorm0 
                      ),
            cms.PSet( detector = cms.string("HCAL_BARREL2_RING1"),
                      recHitEnergyNorm = _rechitEnergyNorm12
                      )
            )
        )
)

pfClusteringHO00 = cms.Sequence( particleFlowRecHitHO00 * particleFlowClusterHO00 )
localReReco00    = cms.Sequence( pfClusteringHO00 )

# Default
#HOThreshold0 = cms.double(1.1),
#HOThresholdPlus1 = cms.double(3.5),
#HOThresholdMinus1 = cms.double(3.5),
#HOThresholdPlus2 = cms.double(3.5),
#HOThresholdMinus2 = cms.double(3.5),
towerMakerWithHO00 = towerMakerWithHO.clone(
    UseHO = _useho,
    HOThreshold0      = _threshold0,
    HOThresholdPlus1  = _threshold12,
    HOThresholdMinus1 = _threshold12,
    HOThresholdPlus2  = _threshold12,
    HOThresholdMinus2 = _threshold12,
    
    hbheInput = cms.InputTag("hbhereco"),
    hoInput   = cms.InputTag("horeco")
)
ak4CaloJets00 = ak4CaloJets.clone(
    src = cms.InputTag("towerMakerWithHO00")
)

# RecoMET/METProducers/python/CaloMET_cfi.py
#caloMetHO00 = caloMetBEFO.clone(
#    src = cms.InputTag("towerMakerWithHO00"),
#    alias = cms.string("caloMetBEFO00")
#)
caloMetHO00 = caloMet.clone(
    src = cms.InputTag("towerMakerWithHO00"),
    alias = cms.string("caloMetHO00")
)

# RecoMET/METProducers/python/MuonMETValueMapProducer_cff.py
muonMETValueMapProducer00 = muonMETValueMapProducer.clone(
   useHO = _useho
)
# RecoMET/METProducers/python/MetMuonCorrections_cff.py 
corMetGlobalMuons00 = corMetGlobalMuons.clone(
    uncorMETInputTag = cms.InputTag("caloMetHO00"),
    muonMETDepositValueMapInputTag = cms.InputTag("muonMETValueMapProducer00","muCorrData","")
)

# RecoMET/METProducers/python/MuonTCMETValueMapProducer_cff.py 
muonTCMETValueMapProducer00 = muonTCMETValueMapProducer.clone(
)

tcMet00 = tcMet.clone(
    alias  = cms.string('tcMet00'),
    metInputTag       = cms.InputTag("caloMetHO00"),
    muonDepValueMap   = cms.InputTag("muonMETValueMapProducer00"  , "muCorrData"),
    tcmetDepValueMap  = cms.InputTag("muonTCMETValueMapProducer00", "muCorrData"),
)

globalReReco00 =  cms.Sequence(
    towerMakerWithHO00
    + ak4CaloJets00
    + caloMetHO00
    + muonMETValueMapProducer00
    + corMetGlobalMuons00
    + muonTCMETValueMapProducer00
    + tcMet00
)

particleFlowBlock00 = particleFlowBlock.clone(
#   PFClustersHO = cms.InputTag("particleFlowClusterHO00"),
    useHO   = _useho,
#   verbose = cms.untracked.bool(True),
#   debug   = cms.untracked.bool(True)
    elementImporters = cms.VPSet(
        cms.PSet( importerName = cms.string("GSFTrackImporter"),
                  source = cms.InputTag("pfTrackElec"),
                  gsfsAreSecondary = cms.bool(False),
                  superClustersArePF = cms.bool(True) ),
        cms.PSet( importerName = cms.string("ConvBremTrackImporter"),
                  source = cms.InputTag("pfTrackElec") ),
        cms.PSet( importerName = cms.string("SuperClusterImporter"),
                  source_eb = cms.InputTag("particleFlowSuperClusterECAL:particleFlowSuperClusterECALBarrel"),
                  source_ee = cms.InputTag("particleFlowSuperClusterECAL:particleFlowSuperClusterECALEndcapWithPreshower"),
                  source_towers = cms.InputTag("towerMaker"),
                  maximumHoverE = cms.double(0.5),
                  minSuperClusterPt = cms.double(00.0),
                  minPTforBypass = cms.double(100.0),
                  superClustersArePF = cms.bool(True) ),
        cms.PSet( importerName = cms.string("ConversionTrackImporter"),
                  source = cms.InputTag("pfConversions") ),
        # V0's not actually used in particle flow block building so far
        #cms.PSet( importerName = cms.string("V0TrackImporter"),
        #          source = cms.InputTag("pfV0") ),
        #NuclearInteraction's also come in Loose and VeryLoose varieties
        cms.PSet( importerName = cms.string("NuclearInteractionTrackImporter"),
                  source = cms.InputTag("pfDisplacedTrackerVertex") ),
        #for best timing GeneralTracksImporter should come after 
        # all secondary track importers                                                           
        cms.PSet( importerName = cms.string("GeneralTracksImporter"),
                  source = cms.InputTag("pfTrack"),
                  muonSrc = cms.InputTag("muons1stStep"),
                  cleanBadConvertedBrems = cms.bool(True),
                  useIterativeTracking = cms.bool(True),
                  DPtOverPtCuts_byTrackAlgo = cms.vdouble(-1.0,-1.0,-1.0,
                                                           1.0,1.0),
                  NHitCuts_byTrackAlgo = cms.vuint32(3,3,3,3,3)
                  ),
        # secondary GSF tracks are also turned off
        #cms.PSet( importerName = cms.string("GSFTrackImporter"), 
        #          source = cms.InputTag("pfTrackElec:Secondary"),
        #          gsfsAreSecondary = cms.bool(True),
        #          superClustersArePF = cms.bool(True) ),
        # to properly set SC based links you need to run ECAL importer
        # after you've imported all SCs to the block
        cms.PSet( importerName = cms.string("ECALClusterImporter"),
                  source = cms.InputTag("particleFlowClusterECAL"),
                  BCtoPFCMap = cms.InputTag('particleFlowSuperClusterECAL:PFClusterAssociationEBEE') ),
        cms.PSet( importerName = cms.string("GenericClusterImporter"),
                  source = cms.InputTag("particleFlowClusterHCAL") ),
        cms.PSet( importerName = cms.string("GenericClusterImporter"),
                  source = cms.InputTag("particleFlowClusterHO00"),
#                  useHO = _useho,
                  ),
         cms.PSet( importerName = cms.string("GenericClusterImporter"),
                   source = cms.InputTag("particleFlowClusterHF") ),
        cms.PSet( importerName = cms.string("GenericClusterImporter"),
                  source = cms.InputTag("particleFlowClusterPS") ),
        ),
)

particleFlowTmp00 = particleFlowTmp.clone(
    blocks = cms.InputTag("particleFlowBlock00"),
    useHO  = _useho,
    tagHO  = _tagho,  # using get_HO_weight() from GMA and Swagata
)

particleFlow00 = pfLinker.clone(
    PFCandidate = cms.VInputTag(cms.InputTag("particleFlowTmp00"))
)

from RecoParticleFlow.PFProducer.particleFlowTmpPtrs_cfi import *
particleFlowPtrs00 = particleFlowTmpPtrs.clone(
    src = cms.InputTag("particleFlow00")
)

from RecoEgamma.EgammaIsolationAlgos.pfBlockBasedIsolation_cfi import *
particleBasedIsolation00 = particleBasedIsolation.clone(
    pfCandidates = cms.InputTag("particleFlow00"),
)

ak4PFJets00 = ak4PFJets.clone(
    rParam       = cms.double(0.5),
    src =cms.InputTag("particleFlow00")
)

pfMet00 = pfMet.clone(
    src = cms.InputTag("particleFlow00"),
    alias = cms.string('pfMet00'),
    jets = cms.InputTag("ak4PFJets00")
)
particleFlowForChargedMET00 = particleFlowForChargedMET.clone(
    PFCollectionLabel = cms.InputTag("particleFlow00"),
    PVCollectionLabel = cms.InputTag("offlinePrimaryVertices"),
    dzCut = cms.double(0.2),
    neutralEtThreshold = cms.double(-1.0)
)
 
pfChMet00 = pfChMet.clone(
    src = cms.InputTag("particleFlowForChargedMET00"),
    alias = cms.string('pfChMet00'),
    globalThreshold = cms.double(0.0),
    calculateSignificance = cms.bool(False),
)


particleFlowReco00 = cms.Sequence(particleFlowBlock00 + particleFlowTmp00)
particleFlowLinks00 = cms.Sequence( particleFlow00 * particleFlowPtrs00 * particleBasedIsolation00 )

recoPFMET00 = cms.Sequence(pfMet00 + particleFlowForChargedMET00 + pfChMet00)


pfReReco00 = cms.Sequence(particleFlowReco00
                          + particleFlowLinks00
                          + ak4PFJets00
                          + recoPFMET00
)
