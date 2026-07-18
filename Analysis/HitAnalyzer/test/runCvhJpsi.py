## Stage-2 CVH refit driver for J/psi ALCARECO (data, Run 2 2016).
## Reads ALCARECOTkAlJpsiMuMu directly; the producer falls back to the
## legacy in-module track-pair loop when srcCandidates is empty, so no
## prior candidate-producer step is required.
##
## Knobs are exposed through VarParsing('analysis'); see the opts.register
## calls below for the full list. Driven by calibration_studies/slurm/.
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os

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
opts.register('fillGradsFactored', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'store per-event gradient + low-rank factored Hessian (H = B^T B)')
opts.register('doMassConstraint', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool, 'apply J/psi mass constraint in the two-track fit')
opts.register('useIdealGeometry', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'Default False (btojpsik option (B), aligned geometry from GT). Set True '
              'only for the AN Stage-1 broken baseline (Stage-2 corrections not applied '
              'here). See openspec/finalize-cvh-producer-15-0-19.')
opts.register('goldenJson', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'optional Golden JSON file to filter run/lumi pre-processing; empty = no filter')
opts.register('useScalarPot3D', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'use the spherical-harmonic scalar-potential field  '
              'in the CVH refit (default; only model supported in this port)')
opts.register('useOpera3D', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'use the full 3D TOSCA volumetric grid (160812) as the baseline '
              'field for the propagator + geopro + globalCor; takes precedence '
              'over useScalarPot3D when True.')
opts.register('scalarPot3DInitFile', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'coefficient dump file produced by mfs/dump_coeffs_for_cmssw.py. '
              'Always required: the residual-correction maker uses it to register '
              'parmtype-14 modes and seed their initial coefficients. Also reused '
              'as the field producer init file when useScalarPot3D=True.')
opts.register('runFDClosure', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'Numerical-FD closure of the per-mode chain rule '
              '(debug; runs once on the first chain-rule site)')
opts.register('epsilonFDClosure', 1e-4, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'eps for the FD closure (used as eps * dB_perMode for each test mode)')
opts.register('numberOfThreads', 1, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'framework numberOfThreads (numberOfStreams follows the same value)')
opts.register('debugPerIterDump', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'per-iteration debug trace (tree vectors + stdout dbgSeed/dbgIter '
              'lines); use together with eventsToProcess on a few events')
opts.register('eventsToProcess', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'comma-separated run:event list to select specific events '
              '(e.g. 278769:15462343,278769:16101980); empty = all')
opts.register('nIters', 10, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'Gauss-Newton iteration cap per constraint phase (default 10 = baseline)')
opts.register('edmConvergence', 1e-5, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'EDM convergence threshold on the reference-state block (default 1e-5; '
              '0 disables early stopping, e.g. for per-iteration trajectory studies)')
opts.register('keepPixelEdgeHits', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'keep pixel hits whose cluster touches the sensor boundary '
              '(isOnEdge) in the fit instead of demoting them to inactive; '
              'the pixelMinSizeX CPE-quality cut applies independently '
              '(default False = baseline)')
opts.register('pixelMinSizeX', 2, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'minimum pixel cluster size in x for a hit to stay in the fit '
              '(default 2 = baseline sizeX>1 cut; 1 admits all clusters)')
_defaultGroupsFile = os.path.join(os.environ.get('CMSSW_BASE', ''),
                                  'src/Analysis/HitAnalyzer/data/materialGroups50.txt')
opts.register('materialGroupsFile', _defaultGroupsFile, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'global material model grouping-tier rules file '
              '(Analysis/HitAnalyzer/data/materialGroups{50,100}.txt); empty = off')
opts.register('perStepFieldModes', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'apply the scalar-potential correction and attribute the per-mode '
              'derivatives per Geant4 step instead of piecewise-constant per leg '
              '(leg-structure-free field attribution; default True)')
opts.register('skipHitlessSurfaces', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'drop hitless module surfaces (dead-module placeholders, '
              'quality-demoted hits) from the fit; propagation goes hit to hit. '
              'Default True; effective only with globalMaterialModel=True '
              '(auto-disabled otherwise)')
opts.register('globalMaterialModel', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'replace the per-module material parameters (parmtype 7) with the '
              'parmtype-15 global material groups of materialGroupsFile '
              '(exclusive switch). Default True (tier-50 groups file from the '
              'release); set False for the legacy per-module parameterisation')
opts.register('propagationDirection', 'anyDirection', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'Geant4ePropagator PropagationDirection. "anyDirection" (default) '
              'picks forward/backward per leg from the target-plane geometry, '
              'recovering legs whose target plane is marginally behind the '
              'state (runaway-leg failure mode); "alongMomentum" is the '
              'legacy forward-only behaviour (bit-identical for all fits '
              'that do not fail with it).')
opts.parseArguments()
if not opts.scalarPot3DInitFile:
    raise SystemExit(
        "scalarPot3DInitFile=<path> is required (coefficient dump file): "
        "the basis evaluator in globalCor needs it for chain-rule columns "
        "even when useOpera3D=True swaps the baseline field model.")

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

if opts.eventsToProcess:
    process.source.eventsToProcess = cms.untracked.VEventRange(
        *[s.strip() for s in opts.eventsToProcess.split(',') if s.strip()])

process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(int(opts.numberOfThreads)),
    numberOfStreams=cms.untracked.uint32(int(opts.numberOfThreads)),
    numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1),
)

# Per-stream CLHEP engine for the residual-maker. The Tier-3 MT path calls
# setG4RandomEngineForStream() at the top of every produce() to wire this
# engine into Geant4's thread-local RNG. Reproducible across thread counts
# because the framework derives per-stream seeds deterministically from the
# initialSeed below + stream index.
process.RandomNumberGeneratorService.globalCor = cms.PSet(
    initialSeed=cms.untracked.uint32(123456789),
    engineName=cms.untracked.string('HepJamesRandom'),
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
    # Low-rank factored Hessian storage (H = B^T B, nRank x nParms):
    # ~9x smaller than hesspackedv at 360 field modes; see
    # ResidualGlobalCorrectionMakerTwoTrackG4e.cc for the rank argument.
    fillGradsFactored=cms.untracked.bool(bool(opts.fillGradsFactored)),
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
    keepPixelEdgeHits=cms.bool(bool(opts.keepPixelEdgeHits)),
    pixelMinSizeX=cms.int32(int(opts.pixelMinSizeX)),
    doVtxConstraint=cms.bool(False),
    doMassConstraint=cms.bool(bool(opts.doMassConstraint)),
    massConstraint=cms.double(3.0969),
    massConstraintWidth=cms.double(1e-5),
    corFiles=cms.vstring(),
    triggers=cms.vstring(*JPSI_TRIGGERS),
    MagneticFieldLabel=cms.string(""),
    # Scalar-potential B-field correction (parmtype-14, absolute-field
    # model). Initial coefficients + basis structure are loaded from a
    # coefficient dump file (mfs/dump_coeffs_for_cmssw.py output). The
    # dump's mode count determines nFieldModes -- use a 50-mode
    # ("custom50": lphi5-base + l=6,m=1; see mfs/CLAUDE.md) dump to keep
    # the per-event Hessian workspace small in MT runs.
    scalarPotentialInitFile=cms.string(opts.scalarPot3DInitFile),
    # Numerical-FD closure (debug only).
    runFDClosure=cms.bool(bool(opts.runFDClosure)),
    epsilonFDClosure=cms.double(float(opts.epsilonFDClosure)),
    debugPerIterDump=cms.bool(bool(opts.debugPerIterDump)),
    nIters=cms.uint32(int(opts.nIters)),
    edmConvergence=cms.double(float(opts.edmConvergence)),
    materialGroupsFile=cms.string(opts.materialGroupsFile),
    globalMaterialModel=cms.bool(bool(opts.globalMaterialModel)),
    perStepFieldModes=cms.bool(bool(opts.perStepFieldModes)),
    skipHitlessSurfaces=cms.bool(bool(opts.skipHitlessSurfaces) and bool(opts.globalMaterialModel)),
    outprefix=cms.untracked.string("globalcor"),
    # MT G4Error master: GlobalCache config for CvhMasterThread. The master
    # spawns a dedicated thread in initializeGlobalCache that builds DDDWorld
    # + master magnetic field BEFORE any TBB worker starts. Each per-stream
    # CvhWorker then attaches per-thread G4 state to it on first produce().
    # This replaces the geopro side-effect dependency that blocked
    # numberOfThreads >= 2 previously.
    #
    # Narrowed to muons -- this runner only propagates J/psi -> mu mu
    # daughters, so the rest of the canonical CVH particle set
    # (gamma, e+-, pi+-, K+-, p, anti_p) is skipped at physics-list
    # construction. Saves the per-thread ProcessManager + process
    # allocations for ~9 unused particles.
    CvhMaster=CvhMasterPSet.clone(Particles=cms.vstring("mu+", "mu-")),
)

# Bring up the labelled 3D field producer and rewire the consumers
# present in this driver (geopro, Geant4ePropagator, and our
# globalCor analyzer). Uses the scalar-potential ScalarPot3D model
# from scalar-potential field model. Independent from
# nano_cff.setup3DFieldForRefit (which assumes the full set of seven
# CVH-side consumers from the NanoAOD configuration).
if opts.useOpera3D:
    # Use the full 3D TOSCA volumetric grid (160812) as the baseline field
    # for the propagator / geopro / globalCor instead of the scalar-potential
    # ScalarPot3D model. Provided as an alternative field-model option for
    # cross-checks and B-field studies.
    from MagneticField.Engine.volumeBasedMagneticField_160812_cfi import \
        VolumeBasedMagneticFieldESProducer as Opera3DMagneticFieldProducer
    from MagneticField.Engine.volumeBasedMagneticField_160812_cfi import magfield as MagneticFieldGeometry
    process.magfield = MagneticFieldGeometry
    process.es_prefer_magfield_cvhrefit = cms.ESPrefer("XMLIdealGeometryESSource", "magfield")
    process.Opera3DMagneticFieldProducer = Opera3DMagneticFieldProducer
    fieldlabel = "grid_160812_3_8t"
    process.Opera3DMagneticFieldProducer.label = fieldlabel
    process.Opera3DMagneticFieldProducer.useParametrizedTrackerField = cms.bool(False)
    # Route the labelled field into the CPEs as in the production data refit
    # (nano_cff.nanoAOD_customizeData) and the 10_6 cross-release driver:
    # the Lorentz-drift in the hit re-evaluation then uses the same field.
    for _cpe in ("stripCPEESProducer", "StripCPEfromTrackAngleESProducer",
                 "siPixelTemplateDBObjectESProducer", "templates"):
        if hasattr(process, _cpe):
            getattr(process, _cpe).MagneticFieldLabel = fieldlabel
elif not opts.useScalarPot3D:
    raise RuntimeError(
        "useScalarPot3D=False is no longer supported; the legacy non-thread-safe "
        "wrapper class is not part of this port. Use the ScalarPot3D model.")
else:
    if not opts.scalarPot3DInitFile:
        raise RuntimeError(
            "useScalarPot3D=True requires scalarPot3DInitFile to point "
            "at a coefficient dump file produced by mfs/dump_coeffs_for_cmssw.py")
    from MagneticField.ParametrizedEngine.parametrizedMagneticField_ScalarPot3D_cfi \
        import ParametrizedMagneticFieldProducer as ScalarPot3DMagneticFieldProducer
    process.ScalarPot3DMagneticFieldProducer = ScalarPot3DMagneticFieldProducer.clone()
    process.ScalarPot3DMagneticFieldProducer.parameters.InitFile = opts.scalarPot3DInitFile
    fieldlabel = "ScalarPot3DMf"
    process.ScalarPot3DMagneticFieldProducer.label = fieldlabel
process.geopro.MagneticFieldLabel = fieldlabel
process.Geant4ePropagator.MagneticFieldLabel = fieldlabel
# The Geant4 master is now the shared EventSetup product from
# cvhMasterESProducer (CvhMasterRecord), consumed by globalCor via esConsumes.
# It builds its master G4 field via SimG4Core's FieldBuilder on top of the same
# labelled magnetic field the propagator consumes.
from TrackPropagation.Geant4e.cvhMasterESProducer_cfi import cvhMasterESProducer
process.cvhMasterESProducer = cvhMasterESProducer.clone()
process.cvhMasterESProducer.MagneticFieldLabel = cms.string(fieldlabel)
# Activate the CVH-specific propagator path: instantiates the custom fluct
# (G4UniversalFluctuationForExtrapolator) and routes its table pointer via
# SetParticleAndCharge. Without this, computeErrorIoni dereferences a null
# fluct->table on the first event.
process.Geant4ePropagator.ForCVH = cms.bool(True)
process.Geant4ePropagator.PropagationDirection = cms.string(opts.propagationDirection)
process.globalCor.MagneticFieldLabel = cms.string(fieldlabel)

process.reconstruction_step = cms.Path(
    # geopro is removed: CvhMasterThread (residual-maker GlobalCache) now
    # owns the G4 world / master magnetic field setup in an MT-safe way.
    # See TrackPropagation/Geant4e/{interface,src}/CvhMaster*.
    process.hltFilter * process.offlineBeamSpot * process.globalCor
)
process.schedule = cms.Schedule(process.reconstruction_step)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
