import FWCore.ParameterSet.Config as cms 
process = cms.Process('jecdb') 
process.load('CondCore.DBCommon.CondDBCommon_cfi') 
process.CondDBCommon.connect = 'sqlite_file:JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29.db' 
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1)) 
process.source = cms.Source('EmptySource') 
process.PoolDBOutputService = cms.Service('PoolDBOutputService', 
   process.CondDBCommon, 
   toPut = cms.VPSet( 
    cms.PSet(
    record = cms.string('AK1PF'), 
    tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AK1PF'), 
    label  = cms.string('AK1PF') 
    ),
    cms.PSet(
    record = cms.string('AK2PF'), 
    tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AK2PF'), 
    label  = cms.string('AK2PF') 
    ),
    cms.PSet(
    record = cms.string('AK3PF'), 
    tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AK3PF'), 
    label  = cms.string('AK3PF') 
    ),
    cms.PSet(
    record = cms.string('AK4PF'), 
    tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AK4PF'), 
    label  = cms.string('AK4PF') 
    ),
    cms.PSet(
    record = cms.string('AK5PF'), 
    tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AK5PF'), 
    label  = cms.string('AK5PF') 
    ),
    cms.PSet(
    record = cms.string('AK6PF'), 
    tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AK6PF'), 
    label  = cms.string('AK6PF') 
    ),
    cms.PSet(
    record = cms.string('AK7PF'), 
    tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AK7PF'), 
    label  = cms.string('AK7PF') 
    ),
    cms.PSet(
    record = cms.string('AK1Calo'), 
    tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AK1Calo'), 
    label  = cms.string('AK1Calo') 
    ),
    cms.PSet(
    record = cms.string('AK2Calo'), 
    tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AK2Calo'), 
    label  = cms.string('AK2Calo') 
    ),
    cms.PSet(
    record = cms.string('AK3Calo'), 
    tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AK3Calo'), 
    label  = cms.string('AK3Calo') 
    ),
    cms.PSet(
    record = cms.string('AK4Calo'), 
    tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AK4Calo'), 
    label  = cms.string('AK4Calo') 
    ),
    cms.PSet(
    record = cms.string('AK5Calo'), 
    tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AK5Calo'), 
    label  = cms.string('AK5Calo') 
    ),
    cms.PSet(
    record = cms.string('AK6Calo'), 
    tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AK6Calo'), 
    label  = cms.string('AK6Calo') 
    ),
    cms.PSet(
    record = cms.string('AK7Calo'), 
    tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AK7Calo'), 
    label  = cms.string('AK7Calo') 
    ),
    cms.PSet(
    record = cms.string('AKVs1PF'), 
    tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AKVs1PF'), 
    label  = cms.string('AKVs1PF') 
    ),
    cms.PSet(
    record = cms.string('AKVs2PF'), 
    tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AKVs2PF'), 
    label  = cms.string('AKVs2PF') 
    ),
    cms.PSet(
    record = cms.string('AKVs3PF'), 
    tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AKVs3PF'), 
    label  = cms.string('AKVs3PF') 
    ),
    cms.PSet(
    record = cms.string('AKVs4PF'), 
    tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AKVs4PF'), 
    label  = cms.string('AKVs4PF') 
    ),
    cms.PSet(
    record = cms.string('AKVs5PF'), 
    tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AKVs5PF'), 
    label  = cms.string('AKVs5PF') 
    ),
    cms.PSet(
    record = cms.string('AKVs6PF'), 
    tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AKVs6PF'), 
    label  = cms.string('AKVs6PF') 
    ),
    cms.PSet(
    record = cms.string('AKVs7PF'), 
    tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AKVs7PF'), 
    label  = cms.string('AKVs7PF') 
    ),
    cms.PSet(
    record = cms.string('AKVs1Calo'), 
    tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AKVs1Calo'), 
    label  = cms.string('AKVs1Calo') 
    ),
    cms.PSet(
    record = cms.string('AKVs2Calo'), 
    tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AKVs2Calo'), 
    label  = cms.string('AKVs2Calo') 
    ),
    cms.PSet(
    record = cms.string('AKVs3Calo'), 
    tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AKVs3Calo'), 
    label  = cms.string('AKVs3Calo') 
    ),
    cms.PSet(
    record = cms.string('AKVs4Calo'), 
    tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AKVs4Calo'), 
    label  = cms.string('AKVs4Calo') 
    ),
    cms.PSet(
    record = cms.string('AKVs5Calo'), 
    tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AKVs5Calo'), 
    label  = cms.string('AKVs5Calo') 
    ),
    cms.PSet(
    record = cms.string('AKVs6Calo'), 
    tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AKVs6Calo'), 
    label  = cms.string('AKVs6Calo') 
    ),
    cms.PSet(
    record = cms.string('AKVs7Calo'), 
    tag    = cms.string('JetCorrectorParametersCollection_JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29_AKVs7Calo'), 
    label  = cms.string('AKVs7Calo') 
    )
  )
) 
 

process.dbWriterAK1PF = cms.EDAnalyzer('JetCorrectorDBWriter', 
   era    = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29'), 
   algo   = cms.untracked.string('AK1PF') 
) 
process.dbWriterAK2PF = cms.EDAnalyzer('JetCorrectorDBWriter', 
   era    = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29'), 
   algo   = cms.untracked.string('AK2PF') 
) 
process.dbWriterAK3PF = cms.EDAnalyzer('JetCorrectorDBWriter', 
   era    = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29'), 
   algo   = cms.untracked.string('AK3PF') 
) 
process.dbWriterAK4PF = cms.EDAnalyzer('JetCorrectorDBWriter', 
   era    = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29'), 
   algo   = cms.untracked.string('AK4PF') 
) 
process.dbWriterAK5PF = cms.EDAnalyzer('JetCorrectorDBWriter', 
   era    = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29'), 
   algo   = cms.untracked.string('AK5PF') 
) 
process.dbWriterAK6PF = cms.EDAnalyzer('JetCorrectorDBWriter', 
   era    = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29'), 
   algo   = cms.untracked.string('AK6PF') 
) 
process.dbWriterAK7PF = cms.EDAnalyzer('JetCorrectorDBWriter', 
   era    = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29'), 
   algo   = cms.untracked.string('AK7PF') 
) 


process.dbWriterAK1Calo = cms.EDAnalyzer('JetCorrectorDBWriter', 
   era    = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29'), 
   algo   = cms.untracked.string('AK1Calo') 
) 
process.dbWriterAK2Calo = cms.EDAnalyzer('JetCorrectorDBWriter', 
   era    = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29'), 
   algo   = cms.untracked.string('AK2Calo') 
) 
process.dbWriterAK3Calo = cms.EDAnalyzer('JetCorrectorDBWriter', 
   era    = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29'), 
   algo   = cms.untracked.string('AK3Calo') 
) 
process.dbWriterAK4Calo = cms.EDAnalyzer('JetCorrectorDBWriter', 
   era    = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29'), 
   algo   = cms.untracked.string('AK4Calo') 
) 
process.dbWriterAK5Calo = cms.EDAnalyzer('JetCorrectorDBWriter', 
   era    = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29'), 
   algo   = cms.untracked.string('AK5Calo') 
) 
process.dbWriterAK6Calo = cms.EDAnalyzer('JetCorrectorDBWriter', 
   era    = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29'), 
   algo   = cms.untracked.string('AK6Calo') 
) 
process.dbWriterAK7Calo = cms.EDAnalyzer('JetCorrectorDBWriter', 
   era    = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29'), 
   algo   = cms.untracked.string('AK7Calo') 
) 


process.dbWriterAKVs1PF = cms.EDAnalyzer('JetCorrectorDBWriter', 
   era    = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29'), 
   algo   = cms.untracked.string('AKVs1PF') 
) 
process.dbWriterAKVs2PF = cms.EDAnalyzer('JetCorrectorDBWriter', 
   era    = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29'), 
   algo   = cms.untracked.string('AKVs2PF') 
) 
process.dbWriterAKVs3PF = cms.EDAnalyzer('JetCorrectorDBWriter', 
   era    = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29'), 
   algo   = cms.untracked.string('AKVs3PF') 
) 
process.dbWriterAKVs4PF = cms.EDAnalyzer('JetCorrectorDBWriter', 
   era    = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29'), 
   algo   = cms.untracked.string('AKVs4PF') 
) 
process.dbWriterAKVs5PF = cms.EDAnalyzer('JetCorrectorDBWriter', 
   era    = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29'), 
   algo   = cms.untracked.string('AKVs5PF') 
) 
process.dbWriterAKVs6PF = cms.EDAnalyzer('JetCorrectorDBWriter', 
   era    = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29'), 
   algo   = cms.untracked.string('AKVs6PF') 
) 
process.dbWriterAKVs7PF = cms.EDAnalyzer('JetCorrectorDBWriter', 
   era    = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29'), 
   algo   = cms.untracked.string('AKVs7PF') 
) 

process.dbWriterAKVs1Calo = cms.EDAnalyzer('JetCorrectorDBWriter', 
   era    = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29'), 
   algo   = cms.untracked.string('AKVs1Calo') 
) 
process.dbWriterAKVs2Calo = cms.EDAnalyzer('JetCorrectorDBWriter', 
   era    = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29'), 
   algo   = cms.untracked.string('AKVs2Calo') 
) 
process.dbWriterAKVs3Calo = cms.EDAnalyzer('JetCorrectorDBWriter', 
   era    = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29'), 
   algo   = cms.untracked.string('AKVs3Calo') 
) 
process.dbWriterAKVs4Calo = cms.EDAnalyzer('JetCorrectorDBWriter', 
   era    = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29'), 
   algo   = cms.untracked.string('AKVs4Calo') 
) 
process.dbWriterAKVs5Calo = cms.EDAnalyzer('JetCorrectorDBWriter', 
   era    = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29'), 
   algo   = cms.untracked.string('AKVs5Calo') 
) 
process.dbWriterAKVs6Calo = cms.EDAnalyzer('JetCorrectorDBWriter', 
   era    = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29'), 
   algo   = cms.untracked.string('AKVs6Calo') 
) 
process.dbWriterAKVs7Calo = cms.EDAnalyzer('JetCorrectorDBWriter', 
   era    = cms.untracked.string('JEC_2011RECO_STARTHI53_LV1_5_3_16_Track8_Jet29'), 
   algo   = cms.untracked.string('AKVs7Calo') 
) 



process.p = cms.Path( 
    process.dbWriterAK1PF *
    process.dbWriterAK2PF *
    process.dbWriterAK3PF *
    process.dbWriterAK4PF *
    process.dbWriterAK5PF *
    process.dbWriterAK6PF *
    process.dbWriterAK7PF *
    process.dbWriterAK1Calo *
    process.dbWriterAK2Calo *
    process.dbWriterAK3Calo *
    process.dbWriterAK4Calo *
    process.dbWriterAK5Calo *
    process.dbWriterAK6Calo *
    process.dbWriterAK7Calo *
    process.dbWriterAKVs1PF *
    process.dbWriterAKVs2PF *
    process.dbWriterAKVs3PF *
    process.dbWriterAKVs4PF *
    process.dbWriterAKVs5PF *
    process.dbWriterAKVs6PF *
    process.dbWriterAKVs7PF *
    process.dbWriterAKVs1Calo *
    process.dbWriterAKVs2Calo *
    process.dbWriterAKVs3Calo *
    process.dbWriterAKVs4Calo *
    process.dbWriterAKVs5Calo *
    process.dbWriterAKVs6Calo *
    process.dbWriterAKVs7Calo
) 
