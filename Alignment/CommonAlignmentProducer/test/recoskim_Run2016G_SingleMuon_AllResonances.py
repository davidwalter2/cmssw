# Auto generated configuration file
# using:
# Revision: 1.19
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v
# with command line options: RECO -s RAW2DIGI,L1Reco,RECO,SKIM:LogError+LogErrorMonitor,ALCA:TkAlKsToPiPi+TkAlLambdaToProtonPi+TkAlDstToD0Pi+TkAlJpsiMuMu,EI,DQM:@rerecoCommon --runUnscheduled --nThreads 8 --data --era Run2_2016 --scenario pp --conditions 106X_dataRun2_v27 --eventcontent AOD,DQM --datatier AOD,DQMIO --customise Configuration/DataProcessing/RecoTLR.customisePostEra_Run2_2016,Configuration/DataProcessing/Utils.addMonitoring,Alignment/CommonAlignmentProducer/alcarecoSplitLevel.setAlcaRecoSplitLevel --filein /store/data/Run2016G/SingleMuon/RAW/v1/000/279/654/00000/8245F1D7-D46B-E611-A448-02163E014137.root -n 1000 --python_filename Alignment/CommonAlignmentProducer/test/recoskim_Run2016G_SingleMuon_AllResonances.py --no_exec
#
# Extension of recoskim_Run2016G_SingleMuon_KsLambdaDst.py adding the
# TkAlJpsiMuMu ALCAREco stream so the new TwoBodyDecayCandidateProducer +
# VertexCompositeCandidateRemapper plumbing can be exercised on real RAW
# data alongside the existing V0 and D* streams.
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run2_2016_cff import Run2_2016

process = cms.Process('RECO',Run2_2016)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.Skims_cff')
process.load('Configuration.StandardSequences.AlCaRecoStreams_cff')
process.load('CommonTools.ParticleFlow.EITopPAG_cff')
process.load('DQMServices.Core.DQMStoreNonLegacy_cff')
process.load('DQMOffline.Configuration.DQMOffline_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/data/Run2016G/SingleMuon/RAW/v1/000/279/654/00000/8245F1D7-D46B-E611-A448-02163E014137.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('RECO nevts:1000'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.AODoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('AOD'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(31457280),
    fileName = cms.untracked.string('RECO_RAW2DIGI_L1Reco_RECO_SKIM_ALCA_EI_DQM.root'),
    outputCommands = process.AODEventContent.outputCommands
)

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('DQMIO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('RECO_RAW2DIGI_L1Reco_RECO_SKIM_ALCA_EI_DQM_inDQM.root'),
    outputCommands = process.DQMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition
process.ALCARECOStreamTkAlDstToD0Pi = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pathALCARECOTkAlDstToD0Pi')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('ALCARECO'),
        filterName = cms.untracked.string('TkAlDstToD0Pi')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('TkAlDstToD0Pi.root'),
    outputCommands = process.OutALCARECOTkAlDstToD0Pi.outputCommands,
)
process.ALCARECOStreamTkAlKsToPiPi = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pathALCARECOTkAlKsToPiPi')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('ALCARECO'),
        filterName = cms.untracked.string('TkAlKsToPiPi')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('TkAlKsToPiPi.root'),
    outputCommands = process.OutALCARECOTkAlKsToPiPi.outputCommands,
)
process.ALCARECOStreamTkAlLambdaToProtonPi = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pathALCARECOTkAlLambdaToProtonPi')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('ALCARECO'),
        filterName = cms.untracked.string('TkAlLambdaToProtonPi')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('TkAlLambdaToProtonPi.root'),
    outputCommands = process.OutALCARECOTkAlLambdaToProtonPi.outputCommands,
)
process.ALCARECOStreamTkAlJpsiMuMu = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pathALCARECOTkAlJpsiMuMu')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('ALCARECO'),
        filterName = cms.untracked.string('TkAlJpsiMuMu')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('TkAlJpsiMuMu.root'),
    outputCommands = process.OutALCARECOTkAlJpsiMuMu.outputCommands,
)
process.SKIMStreamLogError = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pathlogerror')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('RAW-RECO'),
        filterName = cms.untracked.string('LogError')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('LogError.root'),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep edmTriggerResults_*_*_*',
        'keep *_logErrorHarvester_*_*',
    )
)
process.SKIMStreamLogErrorMonitor = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pathlogerrormonitor')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('USER'),
        filterName = cms.untracked.string('LogErrorMonitor')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('LogErrorMonitor.root'),
    outputCommands = cms.untracked.vstring(
        'drop *_*_*_*',
        'keep edmErrorSummaryEntrys_*_*_*'
    )
)

# Other statements
process.ALCARECOEventContent.outputCommands.extend(process.OutALCARECOTkAlLambdaToProtonPi_noDrop.outputCommands)
process.ALCARECOEventContent.outputCommands.extend(process.OutALCARECOTkAlDstToD0Pi_noDrop.outputCommands)
process.ALCARECOEventContent.outputCommands.extend(process.OutALCARECOTkAlKsToPiPi_noDrop.outputCommands)
process.ALCARECOEventContent.outputCommands.extend(process.OutALCARECOTkAlJpsiMuMu_noDrop.outputCommands)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v27', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.eventinterpretaion_step = cms.Path(process.EIsequence)
process.dqmoffline_step = cms.EndPath(process.DQMOfflineCommon)
process.dqmoffline_1_step = cms.EndPath(process.DQMOfflineMuon)
process.dqmoffline_2_step = cms.EndPath(process.DQMOfflineHcal)
process.dqmoffline_3_step = cms.EndPath(process.DQMOfflineJetMET)
process.dqmoffline_4_step = cms.EndPath(process.DQMOfflineEcal)
process.dqmoffline_5_step = cms.EndPath(process.DQMOfflineEGamma)
process.dqmoffline_6_step = cms.EndPath(process.DQMOfflineL1TMuon)
process.dqmoffline_7_step = cms.EndPath(process.DQMOfflineL1TEgamma)
process.dqmoffline_8_step = cms.EndPath(process.DQMOfflineCTPPS)
process.dqmoffline_9_step = cms.EndPath(process.DQMOfflineL1TMonitoring)
process.dqmofflineOnPAT_step = cms.EndPath(process.PostDQMOffline)
process.AODoutput_step = cms.EndPath(process.AODoutput)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)
process.ALCARECOStreamTkAlDstToD0PiOutPath = cms.EndPath(process.ALCARECOStreamTkAlDstToD0Pi)
process.ALCARECOStreamTkAlKsToPiPiOutPath = cms.EndPath(process.ALCARECOStreamTkAlKsToPiPi)
process.ALCARECOStreamTkAlLambdaToProtonPiOutPath = cms.EndPath(process.ALCARECOStreamTkAlLambdaToProtonPi)
process.ALCARECOStreamTkAlJpsiMuMuOutPath = cms.EndPath(process.ALCARECOStreamTkAlJpsiMuMu)
process.SKIMStreamLogErrorOutPath = cms.EndPath(process.SKIMStreamLogError)
process.SKIMStreamLogErrorMonitorOutPath = cms.EndPath(process.SKIMStreamLogErrorMonitor)

# Schedule definition
process.schedule = cms.Schedule(
    process.raw2digi_step,
    process.L1Reco_step,
    process.reconstruction_step,
    process.pathlogerror,
    process.pathlogerrormonitor,
    process.pathALCARECOTkAlLambdaToProtonPi,
    process.pathALCARECOTkAlDstToD0Pi,
    process.pathALCARECOTkAlKsToPiPi,
    process.pathALCARECOTkAlJpsiMuMu,
    process.eventinterpretaion_step,
    process.dqmoffline_step,
    process.dqmoffline_1_step,
    process.dqmoffline_2_step,
    process.dqmoffline_3_step,
    process.dqmoffline_4_step,
    process.dqmoffline_5_step,
    process.dqmoffline_6_step,
    process.dqmoffline_7_step,
    process.dqmoffline_8_step,
    process.dqmoffline_9_step,
    process.dqmofflineOnPAT_step,
    process.AODoutput_step,
    process.DQMoutput_step,
    process.ALCARECOStreamTkAlDstToD0PiOutPath,
    process.ALCARECOStreamTkAlKsToPiPiOutPath,
    process.ALCARECOStreamTkAlLambdaToProtonPiOutPath,
    process.ALCARECOStreamTkAlJpsiMuMuOutPath,
    process.SKIMStreamLogErrorOutPath,
    process.SKIMStreamLogErrorMonitorOutPath,
)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(16)
process.options.numberOfStreams=cms.untracked.uint32(0)
process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1)

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.RecoTLR
from Configuration.DataProcessing.RecoTLR import customisePostEra_Run2_2016

#call to customisation function customisePostEra_Run2_2016 imported from Configuration.DataProcessing.RecoTLR
process = customisePostEra_Run2_2016(process)

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# Automatic addition of the customisation function from Alignment.CommonAlignmentProducer.alcarecoSplitLevel
from Alignment.CommonAlignmentProducer.alcarecoSplitLevel import setAlcaRecoSplitLevel

#call to customisation function setAlcaRecoSplitLevel imported from Alignment.CommonAlignmentProducer.alcarecoSplitLevel
process = setAlcaRecoSplitLevel(process)

# End of customisation functions
#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)


# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
