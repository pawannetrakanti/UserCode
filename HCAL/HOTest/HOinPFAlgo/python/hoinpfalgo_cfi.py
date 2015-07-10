import FWCore.ParameterSet.Config as cms


#demo = cms.EDAnalyzer('HOinPFAlgo'
#)

hoinpf = cms.EDAnalyzer("HOinPFAlgo",
                        doHist     = cms.bool(True),  
                        doGenJets  = cms.bool(True),   
                        doGenMets  = cms.bool(True),   
                        doGenMatch = cms.bool(True),   
                        GenJetTag  = cms.untracked.InputTag("ak4GenJetsNoNu"),
                        GenMetTag  = cms.untracked.InputTag("genMetCalo"),
                        doCaloJets = cms.bool(True),   
                        doCaloMets = cms.bool(True),   
                        doPFJets   = cms.bool(True),                        
                        doPFMets   = cms.bool(True),                        
                        doTCMets   = cms.bool(True),                        
                        JetIdTag   = cms.InputTag("ak4JetID"),
                        CaloJetTag = cms.untracked.InputTag("ak4CaloJets"),
                        CaloMetTag = cms.untracked.InputTag("caloMetBEFO"),
                        PFJetTag   = cms.untracked.InputTag("ak4PFJets"),
                        PFMetTag   = cms.untracked.InputTag("pfMet"),
                        PFChMetTag = cms.untracked.InputTag("pfChMet"),
                        TCMetTag   = cms.untracked.InputTag("tcMet"),
                        JetPtCut   = cms.double(200.),
                        JetEtaCut  = cms.double(1.4),
                        GenPtCut   = cms.double(500.),
                        GenEtaCut  = cms.double(1.4),
)
hoinpf00 = hoinpf.clone(
#    doGenJets  = cms.bool(False),   
    CaloJetTag = cms.untracked.InputTag("ak4CaloJets00"),
    PFJetTag   = cms.untracked.InputTag("ak4PFJets00"),
    CaloMetTag = cms.untracked.InputTag("caloMetHO00"),
    TCMetTag   = cms.untracked.InputTag("tcMet00"),
    PFMetTag   = cms.untracked.InputTag("pfMet00"),
    PFChMetTag = cms.untracked.InputTag("pfChMet00"),
)
hoinpf01 = hoinpf.clone(
#    doGenJets  = cms.bool(False),   
    CaloJetTag = cms.untracked.InputTag("ak4CaloJets01"),
    PFJetTag   = cms.untracked.InputTag("ak4PFJets01"),
    CaloMetTag = cms.untracked.InputTag("caloMetHO01"),
    TCMetTag   = cms.untracked.InputTag("tcMet01"),
    PFMetTag   = cms.untracked.InputTag("pfMet01"),
    PFChMetTag = cms.untracked.InputTag("pfChMet01"),
)
hoinpf02 = hoinpf.clone(
#    doGenJets  = cms.bool(False),   
    CaloJetTag = cms.untracked.InputTag("ak4CaloJets02"),
    PFJetTag   = cms.untracked.InputTag("ak4PFJets02"),
    CaloMetTag = cms.untracked.InputTag("caloMetHO02"),
    TCMetTag   = cms.untracked.InputTag("tcMet02"),
    PFMetTag   = cms.untracked.InputTag("pfMet02"),
    PFChMetTag = cms.untracked.InputTag("pfChMet02"),
)

