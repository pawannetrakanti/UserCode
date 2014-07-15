import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.load("CondCore.DBCommon.CondDBCommon_cfi")

process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1)
        )

process.source = cms.Source("EmptySource")

process.PoolDBESSource = cms.ESSource("PoolDBESSource",
      DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0)
        ),
      timetype = cms.string('runnumber'),
      toGet = cms.VPSet(
      cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AK3PF'),
            label  = cms.untracked.string('AK3PFLocal')
            ),
      cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AK4PF'),
            label  = cms.untracked.string('AK4PFLocal')
            ),
      cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AKVs3PF'),
            label  = cms.untracked.string('AKVs3PFLocal')
            ),
      cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AKVs4PF'),
            label  = cms.untracked.string('AKVs4PFLocal')
            ),                                                                                
      cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AK4Calo'),
            label  = cms.untracked.string('AK4CaloLocal')
            ),                                                                                
      cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AK3Calo'),
            label  = cms.untracked.string('AK3CaloLocal')
            ),                                                                                
      cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AKVs3Calo'),
            label  = cms.untracked.string('AKVs3CaloLocal')
            ),

      cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AKVs4Calo'),
            label  = cms.untracked.string('AKVs4CaloLocal')
            ),                                                                                
      ),

      connect = cms.string('sqlite:JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29.db')
)


process.demo1 = cms.EDAnalyzer('JetCorrectorDBReader', 
        payloadName    = cms.untracked.string('AK3CaloLocal'),
        printScreen    = cms.untracked.bool(False),
        createTextFile = cms.untracked.bool(True),
        globalTag      = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29')
)


process.demo2 = cms.EDAnalyzer('JetCorrectorDBReader', 
        payloadName    = cms.untracked.string('AK4PFLocal'),
        printScreen    = cms.untracked.bool(False),
        createTextFile = cms.untracked.bool(True),
        globalTag      = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29')
)

process.demo3 = cms.EDAnalyzer('JetCorrectorDBReader', 
        payloadName    = cms.untracked.string('AKVs3PFLocal'),
        printScreen    = cms.untracked.bool(False),
        createTextFile = cms.untracked.bool(True),
        globalTag      = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29')
)

process.demo4 = cms.EDAnalyzer('JetCorrectorDBReader', 
        payloadName    = cms.untracked.string('AKVs3CaloLocal'),
        printScreen    = cms.untracked.bool(False),
        createTextFile = cms.untracked.bool(True),
        globalTag      = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29')                               
)

process.p = cms.Path(process.demo1 * process.demo2 * process.demo3 * process.demo4 )
