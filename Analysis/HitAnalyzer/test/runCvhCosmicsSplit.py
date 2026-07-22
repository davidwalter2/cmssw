## Single-track CVH refit driver for cosmic-muon ALCARECO (data, Run 2 2016).
## Smoke test of the CVH machinery on cosmic topology: tracks that do not
## originate from the beamline and cross the upper tracker hemisphere
## outside-in. Enabled by the origin-independent global material model
## (parmtype-15) and the per-step scalar-potential field application;
## anyDirection propagation picks forward/backward per leg.
##
## Input: TkAlCosmicsInCollisions ALCARECO (NoBPTX PD, cosmics reconstructed
## in collision runs -> guaranteed 3.8 T) or TkAlCosmics0T ALCARECO from the
## Cosmics PD (mixed 0 T / 3.8 T runs -- select field-on runs upstream!).
## Both keep full siPixelClusters/siStripClusters; 10_6-written files must be
## repacked at split level 0 before reading here (ROOT #19773).
##
## Differences w.r.t. runCvhSingleTrack.py:
##  - src = cosmic track collection (default ALCARECOTkAlCosmicsInCollisions)
##  - no HLT filter, doTrigger=False (cosmics have no J/psi trigger menu)
##  - offlineBeamSpot still produced from DB (not in the cosmics ALCARECO);
##    it only serves as the PCA reference axis, bsConstraint stays False.
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os

from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
from Configuration.AlCa.GlobalTag import GlobalTag

opts = VarParsing.VarParsing('analysis')
opts.register('input', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'absolute path or root:// URL of ALCARECO file')
opts.register('trackSrc', 'ALCARECOTkAlCosmicsInCollisions', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'cosmic track collection (ALCARECOTkAlCosmicsInCollisions, or '
              'ALCARECOTkAlCosmicsCTF0T etc. for the Cosmics-PD streams)')
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
opts.register('keepPixelEdgeHits', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'keep pixel hits whose cluster touches the sensor boundary '
              '(isOnEdge) in the fit instead of demoting them to inactive')
opts.register('pixelMinSizeX', 2, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'minimum pixel cluster size in x for a hit to stay in the fit')
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
              'damp Gauss-Newton steps from this iteration on (0 = off)')
opts.register('gnDampFactor', 0.5, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'damping factor applied to the GN step when gnDampAfter is active')
opts.register('allowChargeFlipAboveP', 1e9, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'momentum [GeV] above which the Gauss-Newton step may cross '
              'q/p=0 (a genuine charge flip). Crossing q/p=0 is p->inf, benign '
              'at high p where the seed charge is ambiguous; the sign-flip '
              'clamp stays a divergence guard for stiffer tracks. The hard '
              'p>=2 GeV momentum floor always applies. Default 1e9 = never '
              '(legacy). For high-p cosmics try e.g. 100.')
opts.register('seedChargeSign', 1, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'manual two-hypothesis charge test: +1 = nominal seed charge, '
              '-1 = seed the OPPOSITE charge (fit converges to the opposite '
              'minimum). For diagnostics; the automatic one-pass resolution '
              '(twoHypothesisCharge) supersedes this.')
opts.register('twoHypothesisCharge', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'one-pass two-hypothesis charge resolution (default True): when '
              'the nominal fit clamp catches a q/p sign crossing (high-p '
              'charge mis-ID signature), re-fit the opposite charge and keep '
              'the lower-chi2 result. Inert for well-measured tracks.')
opts.register('cosmicSeedFromEntry', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'seed cosmic legs from the muon-ENTRY (higher-y) end and iterate '
              'hits downward so the whole leg propagates FORWARD (default False '
              '= inner seed, upper half backward). Moves the report point to '
              'the entry end; used to validate the backward energy-loss model.')
opts.register('forceCosmic', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'force the outside-in cosmic handling regardless of track algo '
              '(needed for split-track legs whose algo is not recognised as '
              'ctf/cosmics). Default False.')
_defaultGroupsFile = os.path.join(os.environ.get('CMSSW_BASE', ''),
                                  'src/Analysis/HitAnalyzer/data/materialGroups50.txt')
opts.register('materialGroupsFile', _defaultGroupsFile, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'global material model grouping-tier rules file; empty = off')
opts.register('runFDClosure', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'one-shot numerical-FD closure of the field-mode chain rule')
opts.register('epsilonFDClosure', 1e-4, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'eps for the FD closure')
opts.register('perStepFieldModes', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'apply the scalar-potential correction and attribute the per-mode '
              'derivatives per Geant4 step (leg-structure-free; default True)')
opts.register('skipHitlessSurfaces', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'drop hitless module surfaces from the fit; propagation goes hit '
              'to hit. Effective only with globalMaterialModel=True')
opts.register('globalMaterialModel', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'parmtype-15 global material groups (exclusive switch); '
              'REQUIRED in spirit for cosmics -- the legacy per-module '
              'parameterisation assumes prompt tracks from the IP')
opts.register('materialFDGroup', -1, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'one-shot V1 FD closure for this material group (-1 = off)')
opts.register('materialFDEps', 1e-3, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'epsilon for the material FD closure')
opts.register('goodRunsFile', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'text file of run numbers (one per line, # comments ok) to KEEP; '
              'source is restricted to these whole runs via lumisToProcess. '
              'REQUIRED for cosmics from the Cosmics/NoBPTX PDs because the '
              'solenoid was ramped down for parts of era G -- the 3.8T '
              'ScalarPot3D coefficients must not be fit to reduced-field runs '
              '(see repack/goodruns_cosmics_2016GH.txt). Empty = no filter.')
opts.parseArguments()
if not opts.scalarPot3DInitFile:
    raise SystemExit(
        "scalarPot3DInitFile=<path> is required (coefficient dump file)")

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

# ---- Split-track chain -------------------------------------------------
# Refit the cosmic to regenerate the trajectory + association, split into
# top/bottom half-tracks at closest approach (CosmicTrackSplitter), refit each
# leg -> two INDEPENDENT reco::Tracks per cosmic. Both go to the single-track
# CVH fit; comparing their momenta tests resolution + any coherent hit bias.
process.load("RecoTracker.Configuration.RecoTrackerP5_cff")
process.load("RecoTracker.TrackProducer.TrackRefitterP5_cfi")
if not hasattr(process, "MeasurementTrackerEvent"):
    process.load("RecoTracker.MeasurementDet.MeasurementTrackerEventProducer_cfi")
process.MeasurementTrackerEvent.pixelClusterProducer = 'siPixelClusters'
process.MeasurementTrackerEvent.stripClusterProducer = 'siStripClusters'
process.MeasurementTrackerEvent.inactivePixelDetectorLabels = cms.VInputTag()
process.MeasurementTrackerEvent.inactiveStripDetectorLabels = cms.VInputTag()
process.TrackRefitterP5.src = "ALCARECOTkAlCosmicsInCollisions"
process.TrackRefitterP5.TrajectoryInEvent = cms.bool(True)
process.cosmicTrackSplitting.tracks = "TrackRefitterP5"
process.cosmicTrackSplitting.tjTkAssociationMapTag = "TrackRefitterP5"
process.splittedTracksP5.src = "cosmicTrackSplitting"
process.splittedTracksP5.TrajectoryInEvent = cms.bool(True)

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(opts.nEvents))

assert opts.input, "must set input=<path> on the cmsRun command line"
_url = opts.input if opts.input.startswith(("root://", "file:")) else "file:" + opts.input
process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(_url),
    secondaryFileNames=cms.untracked.vstring(),
)

# Field-status filter: keep only whole 3.8T runs. A cosmics ALCARECO file can
# mix good and ramped/off runs (era G), so filter at the source (lumi 0 = "to
# end of run") rather than skip whole files.
if opts.goodRunsFile:
    _goodruns = []
    for _line in open(opts.goodRunsFile):
        _line = _line.split('#', 1)[0].strip()
        if _line:
            _goodruns.append(int(_line))
    if not _goodruns:
        raise SystemExit(f"goodRunsFile {opts.goodRunsFile} contains no run numbers")
    process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(
        *[f"{r}:1-{r}:0" for r in _goodruns])
    print(f"[runCvhCosmics] field filter: keeping {len(_goodruns)} runs "
          f"from {opts.goodRunsFile}")

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

# The cosmics ALCARECO does not keep offlineBeamSpot -- produce it from the
# DB. It only provides the PCA reference axis for the track parameterisation.
process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")

# Single-track CVH refit of cosmic tracks. No HLT filter: the ALCARECO
# selection (pathALCARECOTkAlCosmicsInCollisions) already ran at production.
process.globalCor = cms.EDProducer(
    "ResidualGlobalCorrectionMakerG4e",
    src=cms.InputTag("splittedTracksP5"),
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
    doTrigger=cms.bool(False),
    doRes=cms.bool(False),
    useIdealGeometry=cms.bool(bool(opts.useIdealGeometry)),
    bsConstraint=cms.bool(False),
    applyHitQuality=cms.bool(True),
    keepPixelEdgeHits=cms.bool(bool(opts.keepPixelEdgeHits)),
    pixelMinSizeX=cms.int32(int(opts.pixelMinSizeX)),
    corFiles=cms.vstring(),
    triggers=cms.vstring(),
    MagneticFieldLabel=cms.string(""),
    scalarPotentialInitFile=cms.string(opts.scalarPot3DInitFile),
    nIters=cms.uint32(int(opts.nIters)),
    edmConvergence=cms.double(float(opts.edmConvergence)),
    debugPerIterDump=cms.bool(bool(opts.debugPerIterDump)),
    gnDampAfter=cms.uint32(int(opts.gnDampAfter)),
    gnDampFactor=cms.double(float(opts.gnDampFactor)),
    allowChargeFlipAboveP=cms.double(float(opts.allowChargeFlipAboveP)),
    seedChargeSign=cms.int32(int(opts.seedChargeSign)),
    twoHypothesisCharge=cms.bool(bool(opts.twoHypothesisCharge)),
    cosmicSeedFromEntry=cms.bool(bool(opts.cosmicSeedFromEntry)),
    forceCosmic=cms.bool(True),  # split legs need cosmic handling
    materialGroupsFile=cms.string(opts.materialGroupsFile),
    runFDClosure=cms.bool(bool(opts.runFDClosure)),
    epsilonFDClosure=cms.double(float(opts.epsilonFDClosure)),
    globalMaterialModel=cms.bool(bool(opts.globalMaterialModel)),
    perStepFieldModes=cms.bool(bool(opts.perStepFieldModes)),
    skipHitlessSurfaces=cms.bool(bool(opts.skipHitlessSurfaces) and bool(opts.globalMaterialModel)),
    materialFDGroup=cms.int32(int(opts.materialFDGroup)),
    materialFDEps=cms.double(float(opts.materialFDEps)),
    outprefix=cms.untracked.string("globalcor_cosmicsplit"),
    # MT G4Error master (GlobalCache) -- owns the G4 world / master field.
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
    for _cpe in ("stripCPEESProducer", "StripCPEfromTrackAngleESProducer",
                 "siPixelTemplateDBObjectESProducer", "templates"):
        if hasattr(process, _cpe):
            getattr(process, _cpe).MagneticFieldLabel = fieldlabel
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
    process.offlineBeamSpot * process.MeasurementTrackerEvent
    * process.TrackRefitterP5 * process.cosmicTrackSplitting * process.splittedTracksP5
    * process.globalCor
)
process.schedule = cms.Schedule(process.reconstruction_step)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
