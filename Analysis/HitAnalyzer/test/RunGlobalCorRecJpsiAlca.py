import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
from Configuration.AlCa.GlobalTag import GlobalTag


ONIA_MODE = "jpsi"


def _onia_settings(mode):
    key = mode.lower()
    settings = {
        "jpsi": {
            "src": "ALCARECOTkAlJpsiMuMu",
            "input_files": [
                "/store/mc/RunIISummer20UL16RECO/JPsiToMuMu_Pt8toInf-pythia8/ALCARECO/TkAlJpsiMuMu-106X_mcRun2_asymptotic_v13-v2/2820000/02536347-8E42-C74C-9ED0-7EAF4290973A.root",
            ],
            "triggers": [
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
            ],
        },
        "upsilon": {
            "src": "ALCARECOTkAlUpsilonMuMu",
            "input_files": [
                "/store/mc/RunIISummer20UL16RECO/Upsilon_MultiResonanceFilter_Pt8to13-pythia8/ALCARECO/TkAlUpsilonMuMu-106X_mcRun2_asymptotic_v13-v2/30000/0387C9FE-F2A1-D14A-9064-72B99FF3B170.root",
            ],
            "triggers": [
                "HLT_Dimuon0_Phi_Barrel",
                "HLT_Dimuon0_Upsilon_Muon",
                "HLT_Dimuon13_Upsilon",
                "HLT_Dimuon8_Upsilon_Barrel",
                "HLT_Mu16_TkMu0_dEta18_Onia",
                "HLT_Mu16_TkMu0_dEta18_Phi",
                "HLT_Mu25_TkMu0_dEta18_Onia",
                "HLT_Mu7p5_L2Mu2_Upsilon",
                "HLT_Mu7p5_Track2_Upsilon",
                "HLT_Mu7p5_Track3p5_Upsilon",
                "HLT_Mu7p5_Track7_Upsilon",
                "HLT_QuadMuon0_Dimuon0_Upsilon",
            ],
        },
    }
    if key not in settings:
        raise RuntimeError("Unsupported ONIA_MODE '{}'; use 'jpsi' or 'upsilon'.".format(mode))
    return settings[key]


onia = _onia_settings(ONIA_MODE)


process = cms.Process("RECO2", Run2_2016)

process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.EndOfProcess_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.load("Configuration.StandardSequences.GeometrySimDB_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, "auto:run2_mc", "")
process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(
        record=cms.string("GeometryFileRcd"),
        tag=cms.string("XMLFILE_Geometry_2016_81YV1_Extended2016_mc"),
        label=cms.untracked.string("Extended"),
    ),
)
process.XMLFromDBSource.label = cms.string("Extended")

process.load("TrackPropagation.Geant4e.geantRefit_cff")

process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(*onia["input_files"]),
    secondaryFileNames=cms.untracked.vstring(),
)

process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(1),
    numberOfStreams=cms.untracked.uint32(1),
    numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1),
)

process.configurationMetadata = cms.untracked.PSet(
    annotation=cms.untracked.string("globalCor onia ALCARECO cfg"),
    name=cms.untracked.string("Applications"),
    version=cms.untracked.string("$Revision: 1.19 $"),
)

process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")

process.globalCor = cms.EDProducer(
    "ResidualGlobalCorrectionMakerTwoTrackG4e",
    src=cms.InputTag(onia["src"]),
    fitFromGenParms=cms.bool(False),
    fitFromSimParms=cms.bool(False),
    fillTrackTree=cms.bool(True),
    fillGrads=cms.bool(False),
    fillJac=cms.bool(True),
    fillRunTree=cms.bool(True),
    doGen=cms.bool(False),
    doSim=cms.bool(False),
    requireGen=cms.bool(False),
    doMuons=cms.bool(False),
    doMuonAssoc=cms.bool(False),
    doTrigger=cms.bool(True),
    doRes=cms.bool(False),
    useIdealGeometry=cms.bool(True),
    bsConstraint=cms.bool(False),
    applyHitQuality=cms.bool(True),
    doVtxConstraint=cms.bool(False),
    doMassConstraint=cms.bool(False),
    massConstraint=cms.double(3.0969 if ONIA_MODE.lower() == "jpsi" else 9.4603),
    massConstraintWidth=cms.double(1e-5),
    corFiles=cms.vstring(),
    triggers=cms.vstring(*onia["triggers"]),
    MagneticFieldLabel=cms.string(""),
    scalarPotentialLmax=cms.uint32(5),
    scalarPotentialExtra=cms.vstring(),
    outprefix=cms.untracked.string("globalcor"),
)

process.reconstruction_step = cms.Path(process.geopro * process.offlineBeamSpot * process.globalCor)
process.schedule = cms.Schedule(process.reconstruction_step)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask

associatePatAlgosToolsTask(process)

from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands

process = customiseLogErrorHarvesterUsingOutputCommands(process)

from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete

process = customiseEarlyDelete(process)
