## Stage-2 unified CVH refit driver for B+ -> J/psi K+.
##
## Composes: JpsiKCandidateSplitter -> two-track CVH refit of the J/psi
## (mass-constraint ON, dimuon-side) and single-track CVH refit of the
## bachelor kaon (kaon hypothesis), then writes each maker's sidecar
## TFile separately (offline-join recipe; see §5).
##
## Driver design follows runCvhJpsi.py for the framework/B-field plumbing
## and runCvhJpsiXSmoke.py for the cross-release single-track maker
## conventions (no geopro on the path -- CvhMasterThread owns the G4
## master MT-safely; no HLT filter -- the ALCAReco is already
## trigger-skimmed at Stage 1).
##
## VarParsing knobs:
##   input               : single ALCARECO sl1 file
##   nEvents             : events to process (-1 = all)
##   fillJac             : store per-track Jacobians (default True)
##   scalarPot3DInitFile : coefficient dump path (default = David W. mfs file)
##   mode                : "dimuon" (default) or "kaon" -- ONE maker per
##                         cmsRun job. The two makers cannot co-exist in the
##                         same process: each `CvhMasterThread` constructs
##                         its own `G4MTRunManagerKernel`, which trips the
##                         G4 single-master singleton (Geant4's
##                         `G4Region` ctor segfaults on the second instance,
##                         confirmed empirically). "both" is still accepted
##                         for future-proofing if CvhMasterThread ever gets
##                         a shared-singleton refactor, but at present it
##                         WILL crash; run two jobs instead.
import os

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
from Configuration.AlCa.GlobalTag import GlobalTag

# AN2021_131_v8 §3.3-3.4 canonical correction file. Located relative to
# $WREM_BASE (set by setup.sh); a clear EnvironmentError is raised at
# driver-parse time if the env var is unset or the file is missing.
def _resolve_default_corfile():
    wrem_base = os.environ.get('WREM_BASE')
    if not wrem_base:
        raise EnvironmentError(
            'WREM_BASE environment variable is not set. Source setup.sh at the '
            'repo root before invoking cmsRun so the default corFiles path can '
            'be resolved. (Explicit override: pass corFiles=<path> on the CLI.)'
        )
    rel = 'wremnants-data/data/calibration/correctionResults_v721_recjpsidata.root'
    p = os.path.join(wrem_base, rel)
    if not os.path.isfile(p):
        raise EnvironmentError(
            'Default correction file not found at: {}. Either restore the '
            'wremnants-data submodule (git submodule update --init) or override '
            'via corFiles=<path> on the CLI.'.format(p)
        )
    return p

opts = VarParsing.VarParsing('analysis')
opts.register('input', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'absolute path or root:// URL of the Stage-1 ALCARECO sl1 file')
opts.register('nEvents', 200, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, 'events to process')
opts.register('fillJac', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool, 'store per-track Jacobians')
opts.register('scalarPot3DInitFile',
              '/work/submit/david_w/ZMass/mfs/data/fitresults/polyfit3d_full_coeffs_lmax18_cmsswnorm.txt',
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'scalar-potential coefficient dump')
opts.register('mode', 'dimuon', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'which maker to schedule: "dimuon" or "kaon" (one per cmsRun; '
              '"both" reserved, currently crashes on G4 singleton)')
opts.register('runFDClosure', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'enable the maker\'s finite-difference Jacobian closure test '
              '(once per job; prints "===== Numerical-FD closure =====" to stdout)')
opts.register('epsilonFDClosure', 1e-4, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'epsilon for the FD-closure perturbation (default 1e-4)')
opts.register('useScalarPot3D', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'when True (default), use the scalar-potential 3D field via the '
              '`ScalarPot3DMagneticFieldProducer`; when False, fall back to the '
              'standard CMSSW field from `MagneticField_cff` (used for the '
              '§9.4.b A/B test of the kaon q/p anomaly)')
opts.register('kaonAsMuon', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'when True, override `trackParticleName` on the kaon-side maker '
              'to "mu" (mass-hypothesis A/B test; §9.4.c). Inputs remain the '
              'bachelor kaon tracks, only the propagation hypothesis changes. '
              'NOTE (openspec add-jpsi-x-muons-and-preprod-refinements): the '
              'JpsiKCandidateSplitter now emits `bachelorPdgId` / `muon0PdgId` '
              '/ `muon1PdgId` branches derived from daughter->pdgId(). This '
              'B+ config still uses the hard-coded kaon mass hypothesis; the '
              'multi-channel Stage-2 (Bc/K*0/phi/Ks/Lambda/psi2S) will consume '
              'the pdgId branches to pick the mass hypothesis per event.')
opts.register('plimit', 1.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'Geant4ePropagator PropagationPtotLimit [GeV/c]. Default 1.0 '
              'matches the cfi; §9.8 ramps this down to 0.05 to recover the '
              'soft-bachelor tail. Lowering it on the dimuon side is a no-op '
              '(muons clear 1.0 trivially); the knob affects the kaon mode.')
# openspec/improve-cvh-refit-convergence: CVH joint-refit convergence knobs
# (default values reproduce the published baseline bit-identically).
opts.register('nIters', 10, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'CVH Gauss-Newton iteration cap per icons phase. Default 10 '
              'matches the previous hard-coded constant; matrix sweeps {10, 20, 50, 100}.')
opts.register('edmConvergence', 1e-5, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'CVH convergence threshold on edmval = -delta(chi^2). Default 1e-5 '
              'matches the previous hard-coded constant; matrix sweeps {1e-5, 1e-3, 1e-2}.')
opts.register('useStartingState', 'perigee', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'CVH iter-0 reference state. Currently only "perigee" is implemented; '
              '"midPropagated" is wired but throws cms::Exception until the propagation '
              'helper lands (openspec/improve-cvh-refit-convergence task 1.4).')
opts.register('debug', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'when True, write per-iter vector branches (chisqval_iter, edmval_iter, '
              'deltachisqval_iter, mu_qoverp_iter, Jpsi_mass_iter) for the dimuon-side '
              'maker. Use only for the matrix per-iter deep dive; bloats output ~80 B/event.')
opts.register('useIdealGeometry', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'Default False (btojpsik option (B)): propagate against the aligned geometry '
              'loaded from the CMSSW GlobalTag. Set True for the AN2021_131_v8 canonical '
              'config, in which case populate corFiles (below) with a correctionResults*.root '
              'whose parmset matches this producer build; the WMass-era v721 file does NOT '
              'match this build (idxmaptree remapping not implemented in the current maker), '
              'so True + corFiles=[v721] currently crashes at the parmset-size assert. '
              'See openspec/enable-an-canonical-corrections (follow-up) for the (A) path.')
opts.register('useIdealGeometryMuon', -1, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'per-leg override: -1 (default) inherits from useIdealGeometry, 0 forces False, '
              '1 forces True (int-bool because VarParsing bools cannot express "unset").')
opts.register('useIdealGeometryKaon', -1, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'per-leg override for the bachelor-kaon maker; semantics as useIdealGeometryMuon.')
opts.register('corFiles', [], VarParsing.VarParsing.multiplicity.list,
              VarParsing.VarParsing.varType.string,
              'AN2021_131_v8 Stage-2 correction file list. Default [] (btojpsik option (B), '
              'physics-quality via aligned geometry with no Stage-2 corrections). To attempt '
              'AN-canonical (option (A)) pass a correctionResults*.root whose parmtree row '
              'count matches this producer build (WMass v721 does NOT match; see follow-up).')
opts.register('corFilesMuon', [], VarParsing.VarParsing.multiplicity.list,
              VarParsing.VarParsing.varType.string,
              'per-leg override: empty (default) inherits from corFiles.')
opts.register('corFilesKaon', [], VarParsing.VarParsing.multiplicity.list,
              VarParsing.VarParsing.varType.string,
              'per-leg override for the bachelor-kaon maker; semantics as corFilesMuon.')
opts.register('disableCorFiles', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'explicit knob to disable Stage-2 corrections (both legs). Useful for the '
              'broken-baseline point P3 and for the (B) side of the (A)-vs-(B) overlay.')
opts.parseArguments()

# Resolve corFiles: empty on the CLI --> driver default (v721); non-empty --> take as-is.
# `disableCorFiles=True` --> force empty on both legs regardless of other opts.
def _resolve_corfiles(opts_value, base_default):
    if opts.disableCorFiles:
        return []
    if opts_value:
        return list(opts_value)
    return list(base_default)

_base_corfiles = [] if opts.disableCorFiles else list(opts.corFiles)
_muon_corfiles = _resolve_corfiles(opts.corFilesMuon, _base_corfiles)
_kaon_corfiles = _resolve_corfiles(opts.corFilesKaon, _base_corfiles)

# Resolve per-leg useIdealGeometry: -1 inherits, 0/1 override.
def _resolve_ideal_geom(per_leg, base):
    if per_leg < 0:
        return bool(base)
    return bool(per_leg)

_muon_ideal_geom = _resolve_ideal_geom(opts.useIdealGeometryMuon, opts.useIdealGeometry)
_kaon_ideal_geom = _resolve_ideal_geom(opts.useIdealGeometryKaon, opts.useIdealGeometry)

# Echo the resolved config to stdout so the CMSSW log carries an obvious
# summary line even before the plugin's own LogInfo runs.
print('[runCvhBplusJpsiK.py] resolved CVH config:', flush=True)
print('  muon: useIdealGeometry={}, corFiles={}'.format(_muon_ideal_geom, _muon_corfiles), flush=True)
print('  kaon: useIdealGeometry={}, corFiles={}'.format(_kaon_ideal_geom, _kaon_corfiles), flush=True)
assert opts.input, 'must set input=<path>'
assert opts.mode in ('both', 'dimuon', 'kaon'), \
    f'mode must be both|dimuon|kaon, got {opts.mode!r}'

process = cms.Process('CVHBPLUS', Run2_2016)

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')

process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(
        record=cms.string('GeometryFileRcd'),
        tag=cms.string('XMLFILE_Geometry_2016_81YV1_Extended2016_mc'),
        label=cms.untracked.string('Extended'),
    ),
)
process.XMLFromDBSource.label = cms.string('Extended')

process.load('TrackPropagation.Geant4e.geantRefit_cff')
from TrackPropagation.Geant4e.cvhMaster_cfi import CvhMasterPSet

# ---- source / general options ---------------------------------------------
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(opts.nEvents))
_url = opts.input if opts.input.startswith(('root://', 'file:')) else 'file:' + opts.input
process.source = cms.Source(
    'PoolSource',
    fileNames=cms.untracked.vstring(_url),
    secondaryFileNames=cms.untracked.vstring(),
)
process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(1),
    numberOfStreams=cms.untracked.uint32(1),
    numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 50

# ---- splitter --------------------------------------------------------------
from Analysis.HitAnalyzer.JpsiKCandidateSplitter_cfi import jpsiKCandidateSplitter
process.jpsiKCandidateSplitter = jpsiKCandidateSplitter.clone()

# ---- beamspot --------------------------------------------------------------
process.offlineBeamSpot = cms.EDProducer('BeamSpotProducer')

# ---- two-track maker (J/psi side, mass-constraint ON) ---------------------
from Analysis.HitAnalyzer.ResidualGlobalCorrectionMakerTwoTrackJpsiKMuMuG4e_cfi \
    import globalCorJpsiK
process.globalCorJpsiK = globalCorJpsiK.clone(
    scalarPotentialInitFile=cms.string(opts.scalarPot3DInitFile),
    fillJac=cms.bool(bool(opts.fillJac)),
    runFDClosure=cms.bool(bool(opts.runFDClosure)),
    epsilonFDClosure=cms.double(float(opts.epsilonFDClosure)),
    useIdealGeometry=cms.bool(_muon_ideal_geom),
    corFiles=cms.vstring(*_muon_corfiles),
    # Defaults (10, 1e-5, "perigee", False) reproduce the baseline bit-identically.
    nIters=cms.uint32(int(opts.nIters)),
    edmConvergence=cms.double(float(opts.edmConvergence)),
    useStartingState=cms.string(str(opts.useStartingState)),
    debugPerIterDump=cms.bool(bool(opts.debug)),
    # The CvhMasterThread Particles list is the union of what either maker
    # may propagate: muons (for the dimuon) + kaons (for the bachelor).
    # Same superset on both cfis when scheduled in the same process; the
    # individual per-track Geant4 hypothesis is controlled by the maker
    # (massConstraint=true here, trackParticleName=kaon below).
    CvhMaster=CvhMasterPSet.clone(
        Particles=cms.vstring('mu+', 'mu-', 'kaon+', 'kaon-')),
)
process.RandomNumberGeneratorService.globalCorJpsiK = cms.PSet(
    initialSeed=cms.untracked.uint32(123456789),
    engineName=cms.untracked.string('HepJamesRandom'),
)

# ---- single-track maker (bachelor kaon) -----------------------------------
from Analysis.HitAnalyzer.ResidualGlobalCorrectionMakerJpsiKSingleTrackKaonG4e_cfi \
    import globalCorJpsiKKaon
process.globalCorJpsiKKaon = globalCorJpsiKKaon.clone(
    scalarPotentialInitFile=cms.string(opts.scalarPot3DInitFile),
    fillJac=cms.bool(bool(opts.fillJac)),
    runFDClosure=cms.bool(bool(opts.runFDClosure)),
    epsilonFDClosure=cms.double(float(opts.epsilonFDClosure)),
    useIdealGeometry=cms.bool(_kaon_ideal_geom),
    corFiles=cms.vstring(*_kaon_corfiles),
    # §9.4.c knob: override the propagation hypothesis from "kaon" to "mu"
    # while leaving the input track collection (bachelor kaons) untouched.
    trackParticleName=cms.string('mu' if opts.kaonAsMuon else 'kaon'),
    CvhMaster=CvhMasterPSet.clone(
        Particles=cms.vstring('mu+', 'mu-', 'kaon+', 'kaon-')),
)
process.RandomNumberGeneratorService.globalCorJpsiKKaon = cms.PSet(
    # HepJamesRandom requires the seed in [0, 900_000_000]; pick a distinct
    # value from the dimuon-maker's seed above to keep the two streams
    # independent (per-stream seeds are derived deterministically).
    initialSeed=cms.untracked.uint32(234567890),
    engineName=cms.untracked.string('HepJamesRandom'),
)

# ---- magnetic field --------------------------------------------------------
# Default: load the scalar-potential 3D field producer with a unique label and
# point Geant4ePropagator + makers at it.
# `useScalarPot3D=False` (§9.4.b A/B test) skips the producer entirely;
# `MagneticField_cff` (loaded above) provides the default CMSSW field, and the
# propagator + makers leave `MagneticFieldLabel` empty (== default ESProducer).
if opts.useScalarPot3D:
    from MagneticField.ParametrizedEngine.parametrizedMagneticField_ScalarPot3D_cfi \
        import ParametrizedMagneticFieldProducer as ScalarPot3DMagneticFieldProducer
    process.ScalarPot3DMagneticFieldProducer = ScalarPot3DMagneticFieldProducer.clone()
    process.ScalarPot3DMagneticFieldProducer.parameters.InitFile = opts.scalarPot3DInitFile
    fieldlabel = 'ScalarPot3DMf'
    process.ScalarPot3DMagneticFieldProducer.label = fieldlabel
    process.Geant4ePropagator.MagneticFieldLabel = fieldlabel
    for m in (process.globalCorJpsiK, process.globalCorJpsiKKaon):
        m.MagneticFieldLabel = cms.string(fieldlabel)
        m.CvhMaster.MagneticFieldLabel = cms.string(fieldlabel)
else:
    # Standard CMSSW field. Leave MagneticFieldLabel at its cfi default (empty
    # string -> default ESProducer). ForCVH on the propagator stays on.
    pass
process.Geant4ePropagator.ForCVH = cms.bool(True)
process.Geant4ePropagator.PropagationPtotLimit = cms.double(opts.plimit)

# ---- path / schedule -------------------------------------------------------
# geopro is intentionally NOT on the path: CvhMasterThread (residual-maker
# GlobalCache) owns the G4 world / master magnetic field MT-safely.
if opts.mode == 'both':
    process.reconstruction_step = cms.Path(
        process.offlineBeamSpot
        * process.jpsiKCandidateSplitter
        * process.globalCorJpsiK
        * process.globalCorJpsiKKaon
    )
elif opts.mode == 'dimuon':
    process.reconstruction_step = cms.Path(
        process.offlineBeamSpot
        * process.jpsiKCandidateSplitter
        * process.globalCorJpsiK
    )
else:  # kaon
    process.reconstruction_step = cms.Path(
        process.offlineBeamSpot
        * process.jpsiKCandidateSplitter
        * process.globalCorJpsiKKaon
    )
process.schedule = cms.Schedule(process.reconstruction_step)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
