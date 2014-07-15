import FWCore.ParameterSet.Config as cms
process = cms.Process("jectxt")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('HeavyIonsAnalysis.Configuration.CommonFunctions_cff')
# define your favorite global tag
process.GlobalTag.globaltag = 'STARTHI53_LV1::All'
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))
process.source = cms.Source("EmptySource")
process.readAK3PF    = cms.EDAnalyzer('JetCorrectorDBReader',  
      # below is the communication to the database 
      payloadName    = cms.untracked.string('AK3PF'),
      # this is used ONLY for the name of the printed txt files. You can use any name that you like, 
      # but it is recommended to use the GT name that you retrieved the files from.
      globalTag      = cms.untracked.string('STARTHI53_LV1'),  
      printScreen    = cms.untracked.bool(False),
      createTextFile = cms.untracked.bool(True)
)

process.readAK3Calo   = process.readAK3PF.clone(payloadName = cms.untracked.string('AKVs3Calo'))
process.readAKVs3PF   = process.readAK3PF.clone(payloadName = cms.untracked.string('AKVs3PF'))
process.readAKVs3Calo = process.readAK3PF.clone(payloadName = cms.untracked.string('AKVs3Calo'))

from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import *
overrideJEC_PbPb2760_Test(process)

process.p = cms.Path( process.readAK3PF * process.readAK3Calo * process.readAKVs3PF * process.readAKVs3Calo)
