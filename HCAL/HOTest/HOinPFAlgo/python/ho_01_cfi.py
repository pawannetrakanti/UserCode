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
_tagho        = cms.int32(1) # if _useho is False _tagho is redundant

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
particleFlowRecHitHO01 = particleFlowRecHitHO.clone(
    threshold_ring0 = _threshold0,
    threshold_rig12 = _threshold12
)

#Default
# RING0   seedingThreshold = cms.double(1.0),
# RING1   seedingThreshold = cms.double(3.1),
from RecoParticleFlow.PFClusterProducer.particleFlowClusterHO_cfi import *
particleFlowClusterHO01 = particleFlowClusterHO.clone(
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
        showerSigma = cms.double(01.0),
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

pfClusteringHO01 = cms.Sequence( particleFlowRecHitHO01 * particleFlowClusterHO01 )
localReReco01    = cms.Sequence( pfClusteringHO01 )

# Default
#HOThreshold0 = cms.double(1.1),
#HOThresholdPlus1 = cms.double(3.5),
#HOThresholdMinus1 = cms.double(3.5),
#HOThresholdPlus2 = cms.double(3.5),
#HOThresholdMinus2 = cms.double(3.5),
towerMakerWithHO01 = towerMakerWithHO.clone(
    UseHO = _useho,
    HOThreshold0      = _threshold0,
    HOThresholdPlus1  = _threshold12,
    HOThresholdMinus1 = _threshold12,
    HOThresholdPlus2  = _threshold12,
    HOThresholdMinus2 = _threshold12,
    
    hbheInput = cms.InputTag("hbhereco"),
    hoInput   = cms.InputTag("horeco")
)
ak4CaloJets01 = ak4CaloJets.clone(
    src = cms.InputTag("towerMakerWithHO01")
)

# RecoMET/METProducers/python/CaloMET_cfi.py
#caloMetHO01 = caloMetBEFO.clone(
#    src = cms.InputTag("towerMakerWithHO01"),
#    alias = cms.string("caloMetBEFO01")
#)
caloMetHO01 = caloMet.clone(
    src = cms.InputTag("towerMakerWithHO01"),
    alias = cms.string("caloMetHO01")
)

# RecoMET/METProducers/python/MuonMETValueMapProducer_cff.py
muonMETValueMapProducer01 = muonMETValueMapProducer.clone(
   useHO = _useho
)
# RecoMET/METProducers/python/MetMuonCorrections_cff.py 
corMetGlobalMuons01 = corMetGlobalMuons.clone(
    uncorMETInputTag = cms.InputTag("caloMetHO01"),
    muonMETDepositValueMapInputTag = cms.InputTag("muonMETValueMapProducer01","muCorrData","")
)

# RecoMET/METProducers/python/MuonTCMETValueMapProducer_cff.py 
muonTCMETValueMapProducer01 = muonTCMETValueMapProducer.clone(
)

tcMet01 = tcMet.clone(
    alias  = cms.string('tcMet01'),
    metInputTag       = cms.InputTag("caloMetHO01"),
    muonDepValueMap   = cms.InputTag("muonMETValueMapProducer01"  , "muCorrData"),
    tcmetDepValueMap  = cms.InputTag("muonTCMETValueMapProducer01", "muCorrData"),
)

globalReReco01 =  cms.Sequence(
    towerMakerWithHO01
    + ak4CaloJets01
    + caloMetHO01
    + muonMETValueMapProducer01
    + corMetGlobalMuons01
    + muonTCMETValueMapProducer01
    + tcMet01
)

particleFlowBlock01 = particleFlowBlock.clone(
#   PFClustersHO = cms.InputTag("particleFlowClusterHO01"),
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
                  minSuperClusterPt = cms.double(01.0),
                  minPTforBypass = cms.double(101.0),
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
                  source = cms.InputTag("particleFlowClusterHO01"),
#                  useHO = _useho,
                  ),
         cms.PSet( importerName = cms.string("GenericClusterImporter"),
                   source = cms.InputTag("particleFlowClusterHF") ),
        cms.PSet( importerName = cms.string("GenericClusterImporter"),
                  source = cms.InputTag("particleFlowClusterPS") ),
        ),
)

particleFlowTmp01 = particleFlowTmp.clone(
    blocks = cms.InputTag("particleFlowBlock01"),
    useHO  = _useho,
    tagHO  = _tagho,  # using get_HO_weight() from GMA and Swagata
)

particleFlow01 = pfLinker.clone(
    PFCandidate = cms.VInputTag(cms.InputTag("particleFlowTmp01"))
)

from RecoParticleFlow.PFProducer.particleFlowTmpPtrs_cfi import *
particleFlowPtrs01 = particleFlowTmpPtrs.clone(
    src = cms.InputTag("particleFlow01")
)

from RecoEgamma.EgammaIsolationAlgos.pfBlockBasedIsolation_cfi import *
particleBasedIsolation01 = particleBasedIsolation.clone(
    pfCandidates = cms.InputTag("particleFlow01"),
)

ak4PFJets01 = ak4PFJets.clone(
    rParam       = cms.double(0.5),
    src =cms.InputTag("particleFlow01")
)

pfMet01 = pfMet.clone(
    src = cms.InputTag("particleFlow01"),
    alias = cms.string('pfMet01'),
    jets = cms.InputTag("ak4PFJets01")
)
particleFlowForChargedMET01 = particleFlowForChargedMET.clone(
    PFCollectionLabel = cms.InputTag("particleFlow01"),
    PVCollectionLabel = cms.InputTag("offlinePrimaryVertices"),
    dzCut = cms.double(0.2),
    neutralEtThreshold = cms.double(-1.0)
)
 
pfChMet01 = pfChMet.clone(
    src = cms.InputTag("particleFlowForChargedMET01"),
    alias = cms.string('pfChMet01'),
    globalThreshold = cms.double(0.0),
    calculateSignificance = cms.bool(False),
)


particleFlowReco01 = cms.Sequence(particleFlowBlock01 + particleFlowTmp01)
particleFlowLinks01 = cms.Sequence( particleFlow01 * particleFlowPtrs01 * particleBasedIsolation01 )

recoPFMET01 = cms.Sequence(pfMet01 + particleFlowForChargedMET01 + pfChMet01)


pfReReco01 = cms.Sequence(particleFlowReco01
                          + particleFlowLinks01
                          + ak4PFJets01
                          + recoPFMET01
)
