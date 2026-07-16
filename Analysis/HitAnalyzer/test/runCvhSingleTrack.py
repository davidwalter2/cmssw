## Single-track CVH refit driver for J/psi ALCARECO (data, Run 2 2016).
## Uses ResidualGlobalCorrectionMakerG4e (NOT TwoTrack) so there is no
## kinematic-vertex fit / mass constraint -- each muon track is refit
## independently. Useful for validating the single-track CVH refit and
## as a lighter-weight cross-check of the propagator / residual chain.
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
opts.register('useIdealGeometry', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool, 'use ideal (uncorrected) tracker geometry')
opts.register('useScalarPot3D', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool, 'use the ScalarPot3D field model')
opts.register('useOpera3D', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'use the full 3D TOSCA grid as baseline field (takes precedence over useScalarPot3D)')
opts.register('scalarPot3DInitFile', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'coefficient dump file produced by mfs/dump_coeffs_for_cmssw.py (always required)')
opts.register('numberOfThreads', 1, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'framework numberOfThreads (numberOfStreams follows the same value)')
opts.register('propagationDirection', 'anyDirection', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'Geant4ePropagator PropagationDirection (anyDirection = per-leg '
              'forward/backward choice, default; alongMomentum = legacy '
              'forward-only)')
opts.register('nIters', 10, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'Gauss-Newton iteration cap (default 10 = baseline)')
opts.register('edmConvergence', 1e-5, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'EDM convergence threshold on the reference-state block (default 1e-5)')
opts.register('debugPerIterDump', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'store per-iteration chi2/EDM trajectory vectors in the tree')
opts.register('gnDampAfter', 0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'damp Gauss-Newton steps from this iteration on (0 = off); '
              'collapses limit cycles between chi2-degenerate states')
opts.register('gnDampFactor', 0.5, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'damping factor applied to the GN step when gnDampAfter is active')
opts.parseArguments()
if not opts.scalarPot3DInitFile:
    raise SystemExit(
        "scalarPot3DInitFile=<path> is required (coefficient dump file)")

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
from TrackPropagation.Geant4e.cvhMaster_cfi import CvhMasterPSet

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(opts.nEvents))

assert opts.input, "must set input=<path> on the cmsRun command line"
_url = opts.input if opts.input.startswith(("root://", "file:")) else "file:" + opts.input
process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(_url),
    secondaryFileNames=cms.untracked.vstring(),
)

process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(int(opts.numberOfThreads)),
    numberOfStreams=cms.untracked.uint32(int(opts.numberOfThreads)),
    numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Per-stream CLHEP engine for the residual-maker (required by the MT-safe
# CvhMaster path; seeds are derived deterministically per stream).
process.RandomNumberGeneratorService.globalCor = cms.PSet(
    initialSeed=cms.untracked.uint32(123456789),
    engineName=cms.untracked.string('HepJamesRandom'),
)

process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")

process.hltFilter = cms.EDFilter(
    "HLTHighLevel",
    HLTPaths=cms.vstring(*[t + "_v*" for t in JPSI_TRIGGERS]),
    eventSetupPathsKey=cms.string(""),
    andOr=cms.bool(True),
    throw=cms.bool(False),
    TriggerResultsTag=cms.InputTag("TriggerResults", "", "HLT"),
)

# Single-track CVH refit. No kinematic-vertex fit, no mass constraint --
# each ALCARECOTkAlJpsiMuMu track is refit by itself.
process.globalCor = cms.EDProducer(
    "ResidualGlobalCorrectionMakerG4e",
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
    corFiles=cms.vstring(),
    triggers=cms.vstring(*JPSI_TRIGGERS),
    MagneticFieldLabel=cms.string(""),
    scalarPotentialInitFile=cms.string(opts.scalarPot3DInitFile),
    nIters=cms.uint32(int(opts.nIters)),
    edmConvergence=cms.double(float(opts.edmConvergence)),
    debugPerIterDump=cms.bool(bool(opts.debugPerIterDump)),
    gnDampAfter=cms.uint32(int(opts.gnDampAfter)),
    gnDampFactor=cms.double(float(opts.gnDampFactor)),
    outprefix=cms.untracked.string("globalcor_single"),
    # MT G4Error master (GlobalCache) -- owns the G4 world / master field.
    # Muon-only particle set: this driver refits muon tracks only.
    CvhMaster=CvhMasterPSet.clone(Particles=cms.vstring("mu+", "mu-")),
)

if opts.useOpera3D:
    from MagneticField.Engine.volumeBasedMagneticField_160812_cfi import \
        VolumeBasedMagneticFieldESProducer as Opera3DMagneticFieldProducer
    from MagneticField.Engine.volumeBasedMagneticField_160812_cfi import magfield as MagneticFieldGeometry
    process.magfield = MagneticFieldGeometry
    process.es_prefer_magfield_cvhrefit = cms.ESPrefer("XMLIdealGeometryESSource", "magfield")
    process.Opera3DMagneticFieldProducer = Opera3DMagneticFieldProducer
    fieldlabel = "grid_160812_3_8t"
    process.Opera3DMagneticFieldProducer.label = fieldlabel
    process.Opera3DMagneticFieldProducer.useParametrizedTrackerField = cms.bool(False)
elif not opts.useScalarPot3D:
    raise RuntimeError("useScalarPot3D=False not supported; use ScalarPot3D or Opera3D")
else:
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
process.globalCor.CvhMaster.MagneticFieldLabel = cms.string(fieldlabel)

process.reconstruction_step = cms.Path(
    # geopro removed: CvhMasterThread (GlobalCache) owns the G4 world setup.
    process.hltFilter * process.offlineBeamSpot * process.globalCor
)
process.schedule = cms.Schedule(process.reconstruction_step)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
