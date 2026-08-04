## Single-track CVH gen-closure driver with resolution parameters (doRes).
##
## Re-runs the 2022 resolution-correction closure test with the current
## (2026) hit selection and fit infrastructure: refits gen-matched muon
## tracks of the inclusive B->J/psi+X MC ALCARECO (2016 postVFP, split=1,
## gen kept) with fitFromGenParms=True (reference state frozen to gen ->
## no weak modes, no mass constraint needed) and doRes=True, which
## registers the log-variance resolution parameter families
##   parmtype  8: local-x/phi hit resolution scale (per module)
##   parmtype  9: local-y hit resolution scale (per pixel module)
##   parmtype 10: multiple-scattering process-noise scale (per glued detid)
##   parmtype 11: ionization fluctuation scale (per glued detid)
## and emits their gradient/Hessian rows (incl. the log-det terms) into
## the fillGrads output. The closure criterion: solving the global system
## with any resolution family floated (calibration_studies/
## global_corrections/fit_global_grads.py --parmtypes 10) must return
## parameters consistent with zero. The 2022 attempt failed this because
## non-Gaussian straggling/MS tails pulled the second-moment estimator;
## see Documents/Resolution/NOTES.md.
##
## Assembled from runCvhSingleTrack.py (single-track maker block) and
## runCvhJpsiGenMC.py (MC conditions, source handling, gen anchoring).
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os

from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
from Configuration.AlCa.GlobalTag import GlobalTag

opts = VarParsing.VarParsing('analysis')
opts.register('input', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'comma-separated absolute paths or root:// URLs of ALCARECO files')
opts.register('inputFileList', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'text file with one input path per line (the MC has ~10 events '
              'per file, so runs typically need many files); combined with '
              'input= if both are given')
opts.register('nEvents', 500, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, 'number of events to process (-1 = all)')
opts.register('particle', 'mu', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'track hypothesis + gen-match species: mu (default), pi, kaon, '
              'proton. Sets trackParticleName, the gen-match |pdgId| and the '
              'CvhMaster G4 particle list. The TkAlJpsiX MC track collection '
              'contains the non-muon B daughters, so pi/kaon closures run on '
              'the same input files.')
opts.register('doRes', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'register the resolution parameter families (parmtypes 8-11) '
              'and emit their gradient/Hessian rows (default True; set False '
              'for a doRes-off reference run with identical selection)')
opts.register('doSimHits', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'read tracker PSimHits and match them to hits (input must keep them)')
opts.register('fitSimHitPositions', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'rung-E closure: fit simulated hit positions instead of cluster '
              'positions (covariances unchanged); requires doSimHits=True')
opts.register('trackSrc', 'ALCARECOTkAlJpsiX', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'input track collection (ALCARECOTkAlJpsiX for the custom '
              'B->JpsiX ALCARECO, ALCARECOTkAlJpsiMuMu for standard TkAl)')
opts.register('fitFromGenParms', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'freeze the reference parameters to the gen-muon values '
              '(gen-closure mode, default True)')
opts.register('fillJac', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'store per-track Jacobians (not needed for the grads solve)')
opts.register('fillGrads', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'store per-track gradient + packed Hessian (default True; '
              'input to fit_global_grads.py)')
opts.register('useIdealGeometry', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'use ideal (uncorrected) tracker geometry instead of the '
              'MC-production alignment from the GT (default False)')
opts.register('useScalarPot3D', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'use the spherical-harmonic scalar-potential field model (default)')
opts.register('useOpera3D', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'use the full 3D TOSCA volumetric grid (160812) as baseline field '
              '(takes precedence over useScalarPot3D)')
opts.register('scalarPot3DInitFile', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'coefficient dump file produced by mfs/dump_coeffs_for_cmssw.py '
              '(always required: parmtype-14 registration needs it)')
opts.register('numberOfThreads', 1, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'framework numberOfThreads (numberOfStreams follows the same value)')
opts.register('ioniTruncationAlpha', 0.999, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'delta-electron truncation of the ionization variance in the '
              'G4e error propagation (default 0.999 = historical baseline). '
              'Scan 0.995-0.999: fitted resolution parameters drifting with '
              'this knob is the 2022 failure signature')
opts.register('propagationDirection', 'anyDirection', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'Geant4ePropagator PropagationDirection (anyDirection = per-leg '
              'forward/backward choice, default; alongMomentum = legacy '
              'forward-only)')
opts.register('keepPixelEdgeHits', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'keep pixel hits whose cluster touches the sensor boundary '
              '(default False = baseline exclusion; the 2022 attempt predates '
              'this cut, so the baseline is the interesting configuration)')
opts.register('pixelMinSizeX', 2, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'minimum pixel cluster size in x for a hit to stay in the fit '
              '(default 2 = baseline sizeX>1 cut; 1 admits all clusters)')
opts.register('corFile', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'optional correction file (parmtree/x) applied via corparms_, '
              'e.g. fitted pixel-hit class corrections')
opts.register('nIters', 10, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'Gauss-Newton iteration cap (default 10; gen-anchored fits run 1)')
opts.register('edmConvergence', 1e-5, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'EDM convergence threshold on the reference-state block')
opts.register('gnDampAfter', 0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'damp Gauss-Newton steps from this iteration on (0 = off)')
opts.register('gnDampFactor', 0.5, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'damping factor applied to the GN step when gnDampAfter is active')
opts.register('debugPerIterDump', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'per-iteration debug trace; use with eventsToProcess on a few events')
opts.register('eventsToProcess', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'comma-separated run:event list to select specific events (empty = all)')
_defaultGroupsFile = os.path.join(os.environ.get('CMSSW_BASE', ''),
                                  'src/Analysis/HitAnalyzer/data/materialGroups50.txt')
opts.register('materialGroupsFile', _defaultGroupsFile, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'global material model grouping-tier rules file; empty = off')
opts.register('globalMaterialModel', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'parmtype-15 global material groups instead of per-module '
              'parmtype 7 (default True = current production model)')
opts.register('perStepFieldModes', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'per-Geant4-step field-mode attribution (default True)')
opts.register('skipHitlessSurfaces', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'drop hitless module surfaces from the fit (default True; '
              'effective only with globalMaterialModel=True)')
opts.parseArguments()
if not opts.scalarPot3DInitFile:
    raise SystemExit(
        "scalarPot3DInitFile=<path> is required (coefficient dump file): "
        "the basis evaluator in globalCor needs it for parmtype-14 registration.")

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

# Conditions the MC was produced with (CMSSW_10_6_20_patch1 production
# chain) -- alignment/CPE/beamspot consistent with the simulated detector,
# which is what a gen-closure fit must use.
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
if opts.inputFileList:
    with open(opts.inputFileList) as _f:
        _paths += [l.strip() for l in _f if l.strip() and not l.startswith('#')]
assert _paths, "must set input=<paths> and/or inputFileList=<file> on the cmsRun command line"
_urls = [p if p.startswith(("root://", "file:")) else "file:" + p for p in _paths]
process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(*_urls),
    secondaryFileNames=cms.untracked.vstring(),
    # The condor MC production has a small tail of corrupt files; skip them
    # instead of aborting (pre-scan the filelist for the
    # FormatIncompatibility cases, see calibration_studies/pixelhits).
    skipBadFiles=cms.untracked.bool(True),
    # Every condor job numbers its events from the same (run=1, lumi=1)
    # range, so distinct physics events collide in (run, lumi, event) and
    # the default duplicate check silently drops most of the sample.
    duplicateCheckMode=cms.untracked.string('noDuplicateCheck'),
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

# Per-stream CLHEP engine for the residual-maker (MT-safe CvhMaster path).
process.RandomNumberGeneratorService.globalCor = cms.PSet(
    initialSeed=cms.untracked.uint32(123456789),
    engineName=cms.untracked.string('HepJamesRandom'),
)

process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")

_PARTICLES = {
    'mu':     dict(pdg=13,   g4=["mu+", "mu-"]),
    'pi':     dict(pdg=211,  g4=["pi+", "pi-"]),
    'kaon':   dict(pdg=321,  g4=["kaon+", "kaon-"]),
    'proton': dict(pdg=2212, g4=["proton", "anti_proton"]),
}
if opts.particle not in _PARTICLES:
    raise SystemExit(f"unknown particle={opts.particle} (use mu/pi/kaon/proton)")
_pcfg = _PARTICLES[opts.particle]

# Single-track CVH refit of the gen-matched muon tracks. requireGen with
# the dR<0.1 status-1 |pdgId|==13 matching drops the other B daughters
# (kaon/pion tracks) of the TkAlJpsiX track collection.
process.globalCor = cms.EDProducer(
    "ResidualGlobalCorrectionMakerG4e",
    src=cms.InputTag(opts.trackSrc),
    fitFromGenParms=cms.bool(bool(opts.fitFromGenParms)),
    fitFromSimParms=cms.bool(False),
    fillTrackTree=cms.bool(True),
    fillGrads=cms.bool(bool(opts.fillGrads)),
    fillJac=cms.bool(bool(opts.fillJac)),
    fillRunTree=cms.bool(True),
    doGen=cms.bool(True),
    genParticles=cms.InputTag("genParticles"),
    pileupInfo=cms.InputTag("addPileupInfo"),
    doSim=cms.bool(bool(opts.doSimHits)),
    fitSimHitPositions=cms.untracked.bool(bool(opts.fitSimHitPositions)),
    requireGen=cms.bool(True),
    doMuons=cms.bool(False),
    doMuonAssoc=cms.bool(False),
    doTrigger=cms.bool(False),
    doRes=cms.bool(bool(opts.doRes)),
    useIdealGeometry=cms.bool(bool(opts.useIdealGeometry)),
    bsConstraint=cms.bool(False),
    applyHitQuality=cms.bool(True),
    keepPixelEdgeHits=cms.bool(bool(opts.keepPixelEdgeHits)),
    pixelMinSizeX=cms.int32(int(opts.pixelMinSizeX)),
    corFiles=cms.vstring(*( [opts.corFile] if opts.corFile else [] )),
    triggers=cms.vstring(),
    trackParticleName=cms.string(opts.particle),
    genMatchPdgId=cms.int32(_pcfg['pdg']),
    MagneticFieldLabel=cms.string(""),
    scalarPotentialInitFile=cms.string(opts.scalarPot3DInitFile),
    nIters=cms.uint32(int(opts.nIters)),
    edmConvergence=cms.double(float(opts.edmConvergence)),
    debugPerIterDump=cms.bool(bool(opts.debugPerIterDump)),
    gnDampAfter=cms.uint32(int(opts.gnDampAfter)),
    gnDampFactor=cms.double(float(opts.gnDampFactor)),
    runFDClosure=cms.bool(False),
    epsilonFDClosure=cms.double(1e-4),
    materialGroupsFile=cms.string(opts.materialGroupsFile),
    globalMaterialModel=cms.bool(bool(opts.globalMaterialModel)),
    perStepFieldModes=cms.bool(bool(opts.perStepFieldModes)),
    skipHitlessSurfaces=cms.bool(bool(opts.skipHitlessSurfaces) and bool(opts.globalMaterialModel)),
    outprefix=cms.untracked.string("globalcor_resclosure"),
    # MT G4Error master (GlobalCache); particle set follows the hypothesis.
    CvhMaster=CvhMasterPSet.clone(Particles=cms.vstring(*_pcfg['g4'])),
)

# B-field model wiring (identical to the other 15_0 drivers).
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
    raise RuntimeError(
        "useScalarPot3D=False is no longer supported; use the ScalarPot3D model.")
else:
    from MagneticField.ParametrizedEngine.parametrizedMagneticField_ScalarPot3D_cfi \
        import ParametrizedMagneticFieldProducer as ScalarPot3DMagneticFieldProducer
    process.ScalarPot3DMagneticFieldProducer = ScalarPot3DMagneticFieldProducer.clone()
    process.ScalarPot3DMagneticFieldProducer.parameters.InitFile = opts.scalarPot3DInitFile
    fieldlabel = "ScalarPot3DMf"
    process.ScalarPot3DMagneticFieldProducer.label = fieldlabel
process.geopro.MagneticFieldLabel = fieldlabel
process.Geant4ePropagator.MagneticFieldLabel = fieldlabel
from TrackPropagation.Geant4e.cvhMasterESProducer_cfi import cvhMasterESProducer
process.cvhMasterESProducer = cvhMasterESProducer.clone()
process.cvhMasterESProducer.MagneticFieldLabel = cms.string(fieldlabel)
# CVH-specific propagator path: instantiates the custom fluct
# (G4UniversalFluctuationForExtrapolator). Without this, computeErrorIoni
# dereferences a null fluct->table on the first event.
process.Geant4ePropagator.ForCVH = cms.bool(True)
process.Geant4ePropagator.PropagationDirection = cms.string(opts.propagationDirection)
process.Geant4ePropagator.IoniTruncationAlpha = cms.double(float(opts.ioniTruncationAlpha))
process.globalCor.MagneticFieldLabel = cms.string(fieldlabel)

process.reconstruction_step = cms.Path(process.offlineBeamSpot * process.globalCor)
process.schedule = cms.Schedule(process.reconstruction_step)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
