## Stage-2 CVH refit driver for J/psi ALCARECO (data, Run 2 2016).
## Reads ALCARECOTkAlJpsiMuMu directly; the producer falls back to the
## legacy in-module track-pair loop when srcCandidates is empty, so no
## prior candidate-producer step is required.
##
## Knobs are exposed through VarParsing('analysis'); see the opts.register
## calls below for the full list. Driven by calibration_studies/slurm/.
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
from Configuration.AlCa.GlobalTag import GlobalTag

opts = VarParsing.VarParsing('analysis')
opts.register('input', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'absolute path or root:// URL of ALCARECO file')
opts.register('nEvents', 500, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, 'number of events to process (-1 = all)')
opts.register('fillJac', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool, 'store per-track Jacobians')
opts.register('fillGrads', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool, 'store per-event gradient + packed Hessian')
opts.register('doMassConstraint', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool, 'apply J/psi mass constraint in the two-track fit')
opts.register('useIdealGeometry', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool, 'use ideal (uncorrected) tracker geometry')
opts.register('goldenJson', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'optional Golden JSON file to filter run/lumi pre-processing; empty = no filter')
opts.parseArguments()

JPSI_TRIGGERS = [
    "HLT_Dimuon0_Jpsi_Muon",
    "HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing",
    "HLT_Dimuon0er16_Jpsi_NoVertexing",
    "HLT_Dimuon10_Jpsi_Barrel",
    "HLT_Dimuon13_PsiPrime",
    "HLT_Dimuon16_Jpsi",
    "HLT_Dimuon20_Jpsi",
    "HLT_Dimuon8_PsiPrime_Barrel",
    "HLT_DoubleMu4_3_Bs",
    "HLT_DoubleMu4_3_Jpsi_Displaced",
    "HLT_DoubleMu4_JpsiTrk_Displaced",
    "HLT_DoubleMu4_PsiPrimeTrk_Displaced",
    "HLT_Mu7p5_Track2_Jpsi",
    "HLT_Mu7p5_Track3p5_Jpsi",
    "HLT_Mu7p5_Track7_Jpsi",
]

process = cms.Process("BENCH", Run2_2016)

process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.EndOfProcess_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.GeometrySimDB_cff")

process.GlobalTag = GlobalTag(process.GlobalTag, "auto:run2_data", "")
process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(
        record=cms.string("GeometryFileRcd"),
        tag=cms.string("XMLFILE_Geometry_2016_81YV1_Extended2016_mc"),
        label=cms.untracked.string("Extended"),
    ),
)
process.XMLFromDBSource.label = cms.string("Extended")

process.load("TrackPropagation.Geant4e.geantRefit_cff")

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(opts.nEvents))

assert opts.input, "must set input=<path> on the cmsRun command line"
# Accept local paths (prepend "file:") or xrootd URLs as-is.
_url = opts.input if opts.input.startswith(("root://", "file:")) else "file:" + opts.input
process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(_url),
    secondaryFileNames=cms.untracked.vstring(),
)

# Golden-JSON pre-filter: drop run/lumi pairs that aren't certified.
# Applied at the source so the framework never delivers those events to the
# CVH module, saving CPU on bad lumis.
if opts.goldenJson:
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(
        filename=opts.goldenJson).getVLuminosityBlockRange()

process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(1),
    numberOfStreams=cms.untracked.uint32(1),
    numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1),
)

# Reduce log spam (every 100 events instead of every event).
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")

# HLT pre-filter: drop events that don't pass any of the J/psi paths we
# also store decisions for. Saves the Geant4e/CVH cost on triggers we'd
# never analyse. throw=False so the filter tolerates menu changes across
# eras (any path missing in a given menu is silently skipped).
process.hltFilter = cms.EDFilter(
    "HLTHighLevel",
    HLTPaths=cms.vstring(*[t + "_v*" for t in JPSI_TRIGGERS]),
    eventSetupPathsKey=cms.string(""),
    andOr=cms.bool(True),     # OR over the path list
    throw=cms.bool(False),
    TriggerResultsTag=cms.InputTag("TriggerResults", "", "HLT"),
)

process.globalCor = cms.EDProducer(
    "ResidualGlobalCorrectionMakerTwoTrackG4e",
    src=cms.InputTag("ALCARECOTkAlJpsiMuMu"),
    fitFromGenParms=cms.bool(False),
    fitFromSimParms=cms.bool(False),
    fillTrackTree=cms.bool(True),
    fillGrads=cms.bool(bool(opts.fillGrads)),
    fillJac=cms.bool(bool(opts.fillJac)),
    fillRunTree=cms.bool(True),
    doGen=cms.bool(False),
    doSim=cms.bool(False),
    requireGen=cms.bool(False),
    doMuons=cms.bool(False),
    doMuonAssoc=cms.bool(False),
    doTrigger=cms.bool(True),
    doRes=cms.bool(False),
    useIdealGeometry=cms.bool(bool(opts.useIdealGeometry)),
    bsConstraint=cms.bool(False),
    applyHitQuality=cms.bool(True),
    doVtxConstraint=cms.bool(False),
    doMassConstraint=cms.bool(bool(opts.doMassConstraint)),
    massConstraint=cms.double(3.0969),
    massConstraintWidth=cms.double(1e-5),
    corFiles=cms.vstring(),
    triggers=cms.vstring(*JPSI_TRIGGERS),
    MagneticFieldLabel=cms.string(""),
    outprefix=cms.untracked.string("globalcor"),
)

process.reconstruction_step = cms.Path(
    process.hltFilter * process.geopro * process.offlineBeamSpot * process.globalCor
)
process.schedule = cms.Schedule(process.reconstruction_step)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
