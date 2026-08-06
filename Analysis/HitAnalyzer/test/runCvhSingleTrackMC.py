## Single-track CVH refit driver for JPsiToMuMu MC ALCARECO (UL2016,
## 106X production, repacked to split-1 on ceph — see the repack notes).
## Mirror of runCvhSingleTrack.py with MC conditions and gen matching:
## doGen=True attaches the dR<0.1 same-charge status-1 gen muon to each
## refit track (genPt/genEta/... branches). Primary use: gen-matched
## kink-finder null on MC muons (muons do not decay -> kinkDchisq must be
## chi2-like; deviations = Q-model mismodeling), and the MC reference for
## the pion/kaon decay-tagging studies.
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os

from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
from Configuration.AlCa.GlobalTag import GlobalTag

opts = VarParsing.VarParsing('analysis')
opts.register('input', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'comma-separated absolute paths or root:// URLs of MC ALCARECO files')
opts.register('nEvents', 500, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, 'number of events to process (-1 = all)')
opts.register('useIdealGeometry', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'ideal (uncorrected) tracker geometry; default False = MC-truth '
              'alignment from the MC GlobalTag. NOTE: for the kink null keep '
              'False — the MC GT alignment is the geometry the MC was '
              'simulated with')
opts.register('scalarPot3DInitFile', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'coefficient dump file produced by mfs/dump_coeffs_for_cmssw.py (required)')
opts.register('numberOfThreads', 1, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'framework numberOfThreads (numberOfStreams follows)')
opts.register('propagationDirection', 'anyDirection', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'Geant4ePropagator PropagationDirection')
_defaultGroupsFile = os.path.join(os.environ.get('CMSSW_BASE', ''),
                                  'src/Analysis/HitAnalyzer/data/materialGroups50.txt')
opts.register('materialGroupsFile', _defaultGroupsFile, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'global material model grouping-tier rules file')
opts.register('doKinkFinder', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'per-material-step decay-in-flight score test')
opts.register('kinkInjectLayer', -1, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'closure test: inject a synthetic kink at this step (-1 = off)')
opts.register('kinkInjectDqop', 0.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float, 'injected q/p step [1/GeV]')
opts.register('kinkInjectDxdz', 0.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float, 'injected dx/dz kink')
opts.register('kinkInjectDydz', 0.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float, 'injected dy/dz kink')
opts.parseArguments()
if not opts.scalarPot3DInitFile:
    raise SystemExit("scalarPot3DInitFile=<path> is required (coefficient dump file)")

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

# Conditions the MC was produced with (106X UL2016 MC production chain).
process.GlobalTag = GlobalTag(process.GlobalTag, "106X_mcRun2_asymptotic_v17", "")
process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(
        record=cms.string("GeometryFileRcd"),
        tag=cms.string("XMLFILE_Geometry_2016_81YV1_Extended2016_mc"),
        label=cms.untracked.string("Extended"),
    ),
)
process.XMLFromDBSource.label = cms.string("Extended")

process.load("TrackPropagation.Geant4e.geantRefit_cff")
from TrackPropagation.Geant4e.cvhMaster_cfi import CvhMasterPSet

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(opts.nEvents))

_paths = [p.strip() for p in opts.input.split(',') if p.strip()]
assert _paths, "must set input=<paths> on the cmsRun command line"
_urls = [p if p.startswith(("root://", "file:")) else "file:" + p for p in _paths]
process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(*_urls),
    secondaryFileNames=cms.untracked.vstring(),
)

process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(int(opts.numberOfThreads)),
    numberOfStreams=cms.untracked.uint32(int(opts.numberOfThreads)),
    numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.RandomNumberGeneratorService.globalCor = cms.PSet(
    initialSeed=cms.untracked.uint32(123456789),
    engineName=cms.untracked.string('HepJamesRandom'),
)

process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")

process.globalCor = cms.EDProducer(
    "ResidualGlobalCorrectionMakerG4e",
    src=cms.InputTag("ALCARECOTkAlJpsiMuMu"),
    fitFromGenParms=cms.bool(False),
    fitFromSimParms=cms.bool(False),
    fillTrackTree=cms.bool(True),
    fillGrads=cms.bool(False),
    fillJac=cms.bool(False),
    fillRunTree=cms.bool(True),
    doGen=cms.bool(True),
    genParticles=cms.InputTag("genParticles"),
    pileupInfo=cms.InputTag("addPileupInfo"),
    doSim=cms.bool(False),
    requireGen=cms.bool(False),
    doMuons=cms.bool(False),
    doMuonAssoc=cms.bool(False),
    doTrigger=cms.bool(False),
    doRes=cms.bool(False),
    useIdealGeometry=cms.bool(bool(opts.useIdealGeometry)),
    bsConstraint=cms.bool(False),
    applyHitQuality=cms.bool(True),
    corFiles=cms.vstring(),
    triggers=cms.vstring(),
    trackParticleName=cms.string('mu'),
    MagneticFieldLabel=cms.string(""),
    scalarPotentialInitFile=cms.string(opts.scalarPot3DInitFile),
    materialGroupsFile=cms.string(opts.materialGroupsFile),
    globalMaterialModel=cms.bool(True),
    perStepFieldModes=cms.bool(True),
    skipHitlessSurfaces=cms.bool(True),
    doKinkFinder=cms.bool(bool(opts.doKinkFinder)),
    kinkInjectLayer=cms.int32(int(opts.kinkInjectLayer)),
    kinkInjectDqop=cms.double(float(opts.kinkInjectDqop)),
    kinkInjectDxdz=cms.double(float(opts.kinkInjectDxdz)),
    kinkInjectDydz=cms.double(float(opts.kinkInjectDydz)),
    outprefix=cms.untracked.string("globalcor_singlemc"),
    CvhMaster=CvhMasterPSet.clone(Particles=cms.vstring("mu+", "mu-")),
)

from MagneticField.ParametrizedEngine.parametrizedMagneticField_ScalarPot3D_cfi \
    import ParametrizedMagneticFieldProducer as ScalarPot3DMagneticFieldProducer
process.ScalarPot3DMagneticFieldProducer = ScalarPot3DMagneticFieldProducer.clone()
process.ScalarPot3DMagneticFieldProducer.parameters.InitFile = opts.scalarPot3DInitFile
fieldlabel = "ScalarPot3DMf"
process.ScalarPot3DMagneticFieldProducer.label = fieldlabel

process.geopro.MagneticFieldLabel = fieldlabel
process.Geant4ePropagator.MagneticFieldLabel = fieldlabel
process.Geant4ePropagator.ForCVH = cms.bool(True)
process.Geant4ePropagator.PropagationDirection = cms.string(opts.propagationDirection)
process.globalCor.MagneticFieldLabel = cms.string(fieldlabel)
from TrackPropagation.Geant4e.cvhMasterESProducer_cfi import cvhMasterESProducer
process.cvhMasterESProducer = cvhMasterESProducer.clone()
process.cvhMasterESProducer.MagneticFieldLabel = cms.string(fieldlabel)

process.reconstruction_step = cms.Path(
    process.offlineBeamSpot * process.globalCor
)
process.schedule = cms.Schedule(process.reconstruction_step)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
