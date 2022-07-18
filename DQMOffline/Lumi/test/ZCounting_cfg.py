import FWCore.ParameterSet.Config as cms

process = cms.Process('RECODQM')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/GeometryDB_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')
process.load('TrackingTools/TransientTrack/TransientTrackBuilder_cfi')
process.load('Configuration.StandardSequences.EDMtoMEAtRunEnd_cff')

# load DQM
process.load("DQMServices.Core.DQM_cfg")
process.load("DQMServices.Components.DQMEnvironment_cfi")
process.load('DQMOffline.Configuration.DQMOffline_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.GlobalTag.globaltag = '123X_dataRun3_Prompt_v12' 



# trigger filter
process.load('HLTrigger/HLTfilters/hltHighLevel_cfi')
process.hltHighLevel.throw = cms.bool(False)
process.hltHighLevel.HLTPaths = cms.vstring()


from CondCore.CondDB.CondDB_cfi import *

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # 'file:/eos/cms/store/data/Run2018C/SingleMuon/AOD/12Nov2019_UL2018-v3/240001/AFEC3E96-9E21-9A43-89E4-A8F03C2F4EC0.root'
        'file:/eos/cms/tier0/store/data/Run2022B/SingleMuon/AOD/PromptReco-v1/000/355/680/00000/469e6457-60af-4f1b-aee7-361fdd7d232f.root'
    )
)

process.source.inputCommands = cms.untracked.vstring("keep *",
                                                         "drop *_MEtoEDMConverter_*_*")

process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True),
  Rethrow     = cms.untracked.vstring('ProductNotFound'),
  fileMode    = cms.untracked.string('FULLMERGE')
  )


from DQMOffline.Lumi.ZCounting_cff import zcounting

process.load("DQMOffline.Lumi.ZCounting_cff") 

process.zcounting.MuonTriggerNames = cms.vstring("HLT_IsoMu24_v*")
process.zcounting.MuonTriggerObjectNames = cms.vstring("hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p08")


#process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
#                                     fileName = cms.untracked.string("OUT_step1.root"))

#process.out1 = cms.OutputModule("PoolOutputModule",
#    outputCommands = cms.untracked.vstring('drop *'),
#    fileName = cms.untracked.string('dqm.root')
#)

# Path and EndPath definitions
process.dqmoffline_step = cms.Path(process.zcounting)
process.dqmsave_step = cms.Path(process.DQMSaver)
#process.DQMoutput_step = cms.EndPath(process.DQMoutput)



# Schedule definition
process.schedule = cms.Schedule(
    process.dqmoffline_step,
#    process.DQMoutput_step,
    process.dqmsave_step
    )

# process.dqmSaver.workflow = '/SingleMuon/Run2017B-PromptReco-v1/RECO'
