#from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

#config.section_("General")
config.General.requestName = 'QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8_Scenario3'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
#config.General.transferLogs = True

#config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_JRA_cfg.py'

#config.section_("Data")
config.Data.inputDataset = '/QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8/snabili-Data_analysis-BIGSamples-step3-scenario3-ECAL_aged-HCAL_aged-SiPM_nonaged_FEVTDEBUGHLToutput-c44fb673c801dfcf4c2a42605690049d/USER'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 50
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8_Scenario3'
config.Data.ignoreLocality = True

#config.section_("Site")
config.Site.storageSite = 'T2_IN_TIFR'
#config.Site.whitelist = ["T2_IN_TIFR"]
#config.Site.whitelist = ["T2_FR_IPHC","T2_IN_TIFR","T2_PT_NCG_Lisbon"]
#config.Site.blacklist = ["T0_*","T1_US_FNAL_Disk"]

#config.section_("User")
#config.User.voGroup = 'cms'
