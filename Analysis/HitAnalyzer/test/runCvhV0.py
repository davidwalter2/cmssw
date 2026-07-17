## CVH two-track refit driver for V0 ALCARECO (K_S -> pi pi and
## Lambda -> p pi), data, Run 2 2016. Candidate-driven: reads the
## persisted *Resonances collections from the stage-1 V0 ALCARECO
## production. First consumer of the displaced-track channels the global
## material model was built for -- the origin-independent corrections
## apply to these tracks with no code change.
##
## V0 lesson inherited from the 10_6 drivers: the propagator's momentum
## floor must come down (default 0.5 GeV) -- soft V0 daughters dominate
## the failure accounting otherwise. Exposed as ptotLimit (default
## 0.05 GeV).
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os

from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
from Configuration.AlCa.GlobalTag import GlobalTag

opts = VarParsing.VarParsing('analysis')
opts.register('mode', 'ks', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'ks or lambda')
opts.register('input', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'absolute path or root:// URL of V0 ALCARECO file')
opts.register('nEvents', -1, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, 'number of events to process (-1 = all)')
opts.register('numberOfThreads', 1, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'framework numberOfThreads (numberOfStreams follows)')
opts.register('scalarPot3DInitFile', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'coefficient dump file (always required; registers the parmtype-14 modes)')
opts.register('ptotLimit', 0.05, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'Geant4ePropagator momentum floor in GeV (V0 daughters are soft; '
              'the 0.5 GeV default rejects most of them)')
opts.register('useIdealGeometry', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'ideal (uncorrected) tracker geometry; default False = aligned from GT')
opts.register('propagationDirection', 'anyDirection', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'Geant4ePropagator PropagationDirection')
_defaultGroupsFile = os.path.join(os.environ.get('CMSSW_BASE', ''),
                                  'src/Analysis/HitAnalyzer/data/materialGroups50.txt')
opts.register('materialGroupsFile', _defaultGroupsFile, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'global material model grouping-tier rules file')
opts.register('globalMaterialModel', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'parmtype-15 global material groups (default True -- the whole '
              'point for displaced tracks); False = legacy per-module block')
opts.register('skipHitlessSurfaces', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'hit-to-hit propagation (effective only with globalMaterialModel)')
opts.register('perStepFieldModes', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'per-step field-mode application/attribution')
opts.register('useStartingState', 'perigee', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              "iteration-0 reference state. Default 'perigee': measured BETTER than "
              "'midPropagated' on V0s (KS failures 3.65% vs 7.20% on the 2016 "
              "repacked sample) -- the midPropagated mode, although built for "
              "displaced starting states, degrades soft V0 daughters; to be "
              "understood before changing this default")
opts.register('eventsToProcess', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'comma-separated run:event list; empty = all')
opts.parseArguments()
if not opts.scalarPot3DInitFile:
    raise SystemExit("scalarPot3DInitFile=<path> is required (coefficient dump file)")
if opts.mode not in ('ks', 'lambda'):
    raise SystemExit("mode must be 'ks' or 'lambda'")

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
from TrackPropagation.Geant4e.cvhMaster_cfi import CvhMasterPSet

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(opts.nEvents))

assert opts.input, "must set input=<path> on the cmsRun command line"
_url = opts.input if opts.input.startswith(("root://", "file:")) else "file:" + opts.input
process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(_url),
    secondaryFileNames=cms.untracked.vstring(),
)
if opts.eventsToProcess:
    process.source.eventsToProcess = cms.untracked.VEventRange(
        *[s.strip() for s in opts.eventsToProcess.split(',') if s.strip()])

process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(int(opts.numberOfThreads)),
    numberOfStreams=cms.untracked.uint32(int(opts.numberOfThreads)),
    numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.RandomNumberGeneratorService.globalCor = cms.PSet(
    initialSeed=cms.untracked.uint32(123456789),
    engineName=cms.untracked.string('HepJamesRandom'),
)

process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")

if opts.mode == 'ks':
    from Analysis.HitAnalyzer.ResidualGlobalCorrectionMakerTwoTrackPiPiG4e_cfi import globalCorKs as _v0maker
    _particles = ["pi+", "pi-"]
else:
    from Analysis.HitAnalyzer.ResidualGlobalCorrectionMakerTwoTrackProtonPiG4e_cfi import globalCorLambda as _v0maker
    _particles = ["proton", "anti_proton", "pi+", "pi-"]

process.globalCor = _v0maker.clone(
    useIdealGeometry=cms.bool(bool(opts.useIdealGeometry)),
    fillRunTree=cms.bool(True),
    fillJac=cms.bool(True),
    scalarPotentialInitFile=cms.string(opts.scalarPot3DInitFile),
    materialGroupsFile=cms.string(opts.materialGroupsFile),
    globalMaterialModel=cms.bool(bool(opts.globalMaterialModel)),
    skipHitlessSurfaces=cms.bool(bool(opts.skipHitlessSurfaces) and bool(opts.globalMaterialModel)),
    perStepFieldModes=cms.bool(bool(opts.perStepFieldModes)),
    useStartingState=cms.string(opts.useStartingState),
    outprefix=cms.untracked.string("globalcor_" + opts.mode),
    CvhMaster=CvhMasterPSet.clone(Particles=cms.vstring(*_particles)),
)

# Field label routing (ScalarPot3D as in the J/psi drivers)
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
process.Geant4ePropagator.PropagationPtotLimit = cms.double(float(opts.ptotLimit))
process.globalCor.MagneticFieldLabel = cms.string(fieldlabel)
process.globalCor.CvhMaster.MagneticFieldLabel = cms.string(fieldlabel)

process.reconstruction_step = cms.Path(process.offlineBeamSpot * process.globalCor)
process.schedule = cms.Schedule(process.reconstruction_step)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
