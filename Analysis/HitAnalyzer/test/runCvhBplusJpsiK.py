## Stage-2 unified CVH refit driver for B+ -> J/psi K+.
##
## Composes the two-track CVH refit of the J/psi (mass-constraint ON) and the
## single-track CVH refit of the bachelor kaon (kaon hypothesis). Since the
## shared Geant4 master landed (e33c8d, `cvhMasterESProducer` -> one master per
## job consumed by every maker), both makers run in ONE cmsRun job -- see
## mode="both", now the default. The legacy per-maker modes are kept for A/B
## work; they are what required the offline join.
##
## The two-track maker reads the nested stage-1 candidate directly
## (`ALCARECOTkAlJpsiXBPlusResonances`, whose daughter(0) is the J/psi VCC):
## the maker descends composite daughters itself, so no candidate splitter is
## needed on the dimuon side.
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
##   mode                : "both" (default), "dimuon", or "kaon".
##                         "both" schedules the two makers in one process.
##                         This used to be rejected: each maker built its own
##                         `CvhMasterThread` / `G4MTRunManagerKernel` and
##                         tripped the G4 single-master singleton. That is
##                         fixed -- `cvhMasterESProducer` now supplies ONE
##                         master as an EventSetup product and both makers
##                         consume it via esConsumes.
import os

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
from Configuration.AlCa.GlobalTag import GlobalTag

# AN2021_131_v8 §3.3-3.4 canonical correction file. Located relative to
# $WREM_BASE (set by setup.sh); a clear EnvironmentError is raised at
# driver-parse time if the env var is unset or the file is missing.
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
opts.register('mode', 'both', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'which maker(s) to schedule: "both" (default, one job -- the '
              'makers share the one CvhMaster ES product), "dimuon", or "kaon"')
opts.register('emitRefitTracks', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'single-track maker emits refit reco::Tracks (+refitOk map) for '
              'the downstream constrained B-vertex fit')
opts.register('refitMaxRelPtErr', -1.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'reject a refit track whose ptError/pt exceeds this (refitOk=0). '
              'Negative (default) disables the cut, leaving the covariance tail '
              'available for the follow-up investigation.')
opts.register('jpsiConstraint', 'inFit', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'B-candidate kinematic fit mode: inFit (Bmm5, constraint applied '
              'inside the mother fit) | upstream (muons already CVH '
              'J/psi-constrained, no second constraint) | cascade (dimuon fit '
              'first, then combined with the bachelor, vertex floating)')
opts.register('nanoOut', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'when set, write a NanoAOD to this path: candidate table with '
              'raw VCC columns plus the refit ValueMaps as externalVariables, '
              'and Track / PV / Muon tables. Turns on produceValueMaps and '
              'turns off the sidecar TTree.')
opts.register('srcTracks', 'ALCARECOTkAlJpsiX',
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'persisted AlCaReco track collection (also the dE/dx ValueMap '
              'key and the NanoAOD Track table source)')
opts.register('srcMuons', 'ALCARECOTkAlJpsiXLooseMuons',
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'persisted reco::Muon collection for the NanoAOD Muon table')
opts.register('srcCandidates', 'ALCARECOTkAlJpsiXBPlusResonances',
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'stage-1 candidate collection for the two-track maker. Default is '
              'the nested B+ VCC read directly (the maker descends composite '
              'daughters); set to "jpsiKCandidateSplitter:dimuon" for the '
              'legacy pre-split input')
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
opts.register('plimit', 0.05, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'Geant4ePropagator PropagationPtotLimit [GeV/c]. Default 0.05: '
              'B+ bachelor kaons are soft (p ~ 0.3-0.9 GeV) and the old cfi '
              'default of 1.0 aborted ~2/3 of them with fail[plimit]. '
              'The cfi value was 1.0; §9.8 ramps this down to 0.05 to recover the '
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
              'CVH iter-0 reference state: "perigee" (default) or "midPropagated" '
              '(analytical extrapolation of both daughter perigee states to their '
              'midpoint, with per-event perigee fallback).')
opts.register('debug', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'when True, write per-iter vector branches (chisqval_iter, edmval_iter, '
              'deltachisqval_iter, mu_qoverp_iter, Jpsi_mass_iter) for the dimuon-side '
              'maker. Use only for the matrix per-iter deep dive; bloats output ~80 B/event.')
# openspec/add-btojpsik-cvh-global-calibration: global-correction calibration
# output. Defaults are all OFF/empty so the diagnostic offline-join behaviour
# (the only prior use of this driver) stays bit-identical.
opts.register('fillGrads', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'write the per-candidate gradient/Hessian sidecar branches '
              '(nParms, globalidxv, gradv, hesspackedv/hessfactorv) needed by '
              'the global-correction aggregation. Default False = diagnostic '
              'offline-join output only.')
opts.register('fillGradsFactored', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'store the per-candidate Hessian as its matrix square root '
              '(hessfactorv, H = B^T B) instead of the packed upper triangle. '
              'Only meaningful when fillGrads=True.')
opts.register('fillRunTree', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'write the runtree global-parameter catalog (one row per global '
              'parameter). Needed once per production; the aggregation and '
              'solve stages read it to interpret globalidxv.')
opts.register('globalMaterialModel', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'fit material as ~50 grouped scales (parmtype 15) instead of the '
              'per-module energy-loss form (parmtype 7). Requires '
              'materialGroupsFile. Implies skipHitlessSurfaces=True.')
opts.register('materialGroupsFile', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'path to the material grouping rules file (e.g. materialGroups50.txt). '
              'Required when globalMaterialModel=True. MUST be identical across '
              'every maker whose grads are summed, or the parameter catalogs '
              'will not line up.')
opts.register('globalTag', 'auto:run2_data', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'CMSSW GlobalTag string. Default auto:run2_data for the R2016H ALCARECO. '
              'For MC ALCARECO pass the MC GT that produced the sample, e.g. '
              '106X_mcRun2_asymptotic_v17 for 2016 postVFP.')
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

# Resolve corFiles: empty on the CLI --> no corrections; non-empty --> take as-is.
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

process.GlobalTag = GlobalTag(process.GlobalTag, opts.globalTag, '')
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

# ---- calibration (global-correction) output -------------------------------
# Shared by BOTH makers: the dimuon and kaon grads are summed downstream in a
# common global-parameter index space, so the material model and grouping file
# must be byte-identical between them (and with any other production whose
# grads enter the same fit).
if opts.globalMaterialModel and not opts.materialGroupsFile:
    raise ValueError(
        'globalMaterialModel=True requires materialGroupsFile=<path> '
        '(e.g. materialGroups50.txt)')
if opts.materialGroupsFile and not os.path.isfile(opts.materialGroupsFile):
    raise EnvironmentError(
        f'materialGroupsFile not found: {opts.materialGroupsFile}')

_calib_pset = dict(
    fillGrads=cms.bool(bool(opts.fillGrads)),
    # NOTE: the maker reads this one with getUntrackedParameter, so it must be
    # cms.untracked.bool -- a tracked bool lands in a different namespace and
    # would silently fall back to False (packed instead of factored Hessian).
    fillGradsFactored=cms.untracked.bool(bool(opts.fillGradsFactored)),
    fillRunTree=cms.bool(bool(opts.fillRunTree)),
    globalMaterialModel=cms.bool(bool(opts.globalMaterialModel)),
    materialGroupsFile=cms.string(str(opts.materialGroupsFile)),
)
if opts.globalMaterialModel:
    # the maker defaults skipHitlessSurfaces to globalMaterialModel, but set it
    # explicitly so the configuration is self-documenting in the job log
    _calib_pset['skipHitlessSurfaces'] = cms.bool(True)

print('[runCvhBplusJpsiK] calibration output: '
      f'fillGrads={bool(opts.fillGrads)} '
      f'fillGradsFactored={bool(opts.fillGradsFactored)} '
      f'fillRunTree={bool(opts.fillRunTree)} '
      f'globalMaterialModel={bool(opts.globalMaterialModel)} '
      f'materialGroupsFile="{opts.materialGroupsFile}"')

# ---- bachelor track source -------------------------------------------------
# Splitter removal: the two-track maker reads the nested candidate directly and
# CandidateLeafTrackProducer extracts the bachelor (direct leaf-daughter) tracks
# generically -- so JpsiKCandidateSplitter is only needed for the legacy
# pre-split A/B path (srcCandidates=jpsiKCandidateSplitter:*).
_uses_splitter = 'jpsiKCandidateSplitter' in opts.srcCandidates
if _uses_splitter:
    from Analysis.HitAnalyzer.JpsiKCandidateSplitter_cfi import jpsiKCandidateSplitter
    process.jpsiKCandidateSplitter = jpsiKCandidateSplitter.clone()
    _bach_src = cms.InputTag('jpsiKCandidateSplitter', 'bachelor')
    _bach_idx = cms.InputTag('jpsiKCandidateSplitter', 'bachelorBCandIdx')
else:
    process.bplusBachelorTracks = cms.EDProducer(
        'CandidateLeafTrackProducer', src=cms.InputTag(opts.srcCandidates))
    _bach_src = cms.InputTag('bplusBachelorTracks')
    _bach_idx = cms.InputTag('bplusBachelorTracks', 'candIdx')

# ---- beamspot --------------------------------------------------------------
process.offlineBeamSpot = cms.EDProducer('BeamSpotProducer')

# ---- two-track maker (J/psi side, mass-constraint ON) ---------------------
from Analysis.HitAnalyzer.ResidualGlobalCorrectionMakerTwoTrackJpsiKMuMuG4e_cfi \
    import globalCorJpsiK
_src_cands = cms.InputTag(*opts.srcCandidates.split(':')) \
    if ':' in opts.srcCandidates else cms.InputTag(opts.srcCandidates)
# Reading the nested stage-1 B+ VCC directly: the maker descends daughter(0)
# (the J/psi composite) itself, so bCandIdx bookkeeping is unnecessary -- the
# ValueMaps key straight to this collection.
process.globalCorJpsiK = globalCorJpsiK.clone(
    srcCandidates=_src_cands,
    bCandIdxSrc=(globalCorJpsiK.bCandIdxSrc if _uses_splitter else cms.InputTag('')),
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
    **_calib_pset,
)
process.RandomNumberGeneratorService.globalCorJpsiK = cms.PSet(
    initialSeed=cms.untracked.uint32(123456789),
    engineName=cms.untracked.string('HepJamesRandom'),
)

# ---- single-track maker (bachelor kaon) -----------------------------------
from Analysis.HitAnalyzer.ResidualGlobalCorrectionMakerJpsiKSingleTrackKaonG4e_cfi \
    import globalCorJpsiKKaon
process.globalCorJpsiKKaon = globalCorJpsiKKaon.clone(
    src=_bach_src,
    bCandIdxSrc=_bach_idx,
    scalarPotentialInitFile=cms.string(opts.scalarPot3DInitFile),
    fillJac=cms.bool(bool(opts.fillJac)),
    runFDClosure=cms.bool(bool(opts.runFDClosure)),
    epsilonFDClosure=cms.double(float(opts.epsilonFDClosure)),
    useIdealGeometry=cms.bool(_kaon_ideal_geom),
    corFiles=cms.vstring(*_kaon_corfiles),
    # §9.4.c knob: override the propagation hypothesis from "kaon" to "mu"
    # while leaving the input track collection (bachelor kaons) untouched.
    trackParticleName=cms.string('mu' if opts.kaonAsMuon else 'kaon'),
    emitRefitTracks=cms.bool(bool(opts.emitRefitTracks)),
    refitMaxRelPtErr=cms.double(float(opts.refitMaxRelPtErr)),
    CvhMaster=CvhMasterPSet.clone(
        Particles=cms.vstring('mu+', 'mu-', 'kaon+', 'kaon-')),
    **_calib_pset,
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
else:
    # Standard CMSSW field. Leave MagneticFieldLabel at its cfi default (empty
    # string -> default ESProducer). ForCVH on the propagator stays on.
    fieldlabel = ''

# Shared CVH G4 master (EventSetup product on CvhMasterRecord), consumed by
# every maker via esConsumes -- so it must exist whichever field is in use, and
# it is what lets both makers run in one job on one G4 master.
from TrackPropagation.Geant4e.cvhMasterESProducer_cfi import cvhMasterESProducer
process.cvhMasterESProducer = cvhMasterESProducer.clone()
process.cvhMasterESProducer.MagneticFieldLabel = cms.string(fieldlabel)
process.Geant4ePropagator.ForCVH = cms.bool(True)
process.Geant4ePropagator.PropagationPtotLimit = cms.double(opts.plimit)

# ---- path / schedule -------------------------------------------------------
# geopro is intentionally NOT on the path: CvhMasterThread (residual-maker
# GlobalCache) owns the G4 world / master magnetic field MT-safely.
#
# mode="both" runs the two makers in ONE process: `cvhMasterESProducer`
# supplies a single Geant4 master as an EventSetup product and both makers
# consume it, so there is no second G4 kernel and no offline join.
#
# Splitter removed from the default path: the two-track maker descends the
# nested candidate, and the kaon maker's bachelor tracks come from
# CandidateLeafTrackProducer. The splitter only reappears on the legacy
# pre-split A/B path (_uses_splitter).
_seq = process.offlineBeamSpot
if _uses_splitter:
    _seq = _seq * process.jpsiKCandidateSplitter
elif opts.mode in ('both', 'kaon'):
    _seq = _seq * process.bplusBachelorTracks
if opts.mode in ('both', 'dimuon'):
    _seq = _seq * process.globalCorJpsiK
if opts.mode in ('both', 'kaon'):
    _seq = _seq * process.globalCorJpsiKKaon
process.reconstruction_step = cms.Path(_seq)

# ---- NanoAOD output --------------------------------------------------------
# The refit's corrected quantities are already EDM ValueMaps keyed to
# srcCandidates, so they attach to the candidate table as externalVariables --
# the same mechanism the dE/dx maps use on the Track table. Raw VCC columns and
# corrected columns therefore sit side by side, with no offline join.
if opts.nanoOut:
    from PhysicsTools.NanoAOD.common_cff import Var, ExtVar, P3Vars

    process.globalCorJpsiK.produceValueMaps = cms.bool(True)
    process.globalCorJpsiK.fillTrackTree = cms.bool(False)
    process.globalCorJpsiKKaon.fillTrackTree = cms.bool(False)

    _cor = 'globalCorJpsiK'
    process.bplusTable = cms.EDProducer(
        'SimpleVertexCompositeCandidateFlatTableProducer',
        src=_src_cands,
        cut=cms.string(''),
        name=cms.string('BuJpsiK'),
        doc=cms.string('B+ -> J/psi K+ stage-1 candidates, raw + CVH-refit'),
        singleton=cms.bool(False), extension=cms.bool(False),
        variables=cms.PSet(
            P3Vars,
            mass=Var('mass', float, doc='raw candidate mass'),
            charge=Var('charge', 'int16', doc='charge'),
            vertexChi2=Var('vertexChi2', float, doc='raw vertex chi2'),
            nDau=Var('numberOfDaughters', 'int16', doc='n daughters'),
            jpsiPt=Var('daughter(0).pt', float, doc='J/psi (daughter0) pt'),
            jpsiPdgId=Var('daughter(0).pdgId', int, doc='J/psi pdgId'),
            kaonPt=Var('daughter(1).pt', float, doc='bachelor pt'),
            kaonEta=Var('daughter(1).eta', float, doc='bachelor eta'),
            kaonPhi=Var('daughter(1).phi', float, doc='bachelor phi'),
            kaonPdgId=Var('daughter(1).pdgId', int, doc='bachelor signed pdgId'),
        ),
        externalVariables=cms.PSet(
            corMass=ExtVar(cms.InputTag(_cor, 'corMass'), float, doc='CVH-refit dimuon mass'),
            corMassErr=ExtVar(cms.InputTag(_cor, 'corMassErr'), float, doc='CVH-refit mass error'),
            corPt=ExtVar(cms.InputTag(_cor, 'corPt'), float, doc='CVH-refit dimuon pt'),
            corEta=ExtVar(cms.InputTag(_cor, 'corEta'), float, doc='CVH-refit dimuon eta'),
            corPhi=ExtVar(cms.InputTag(_cor, 'corPhi'), float, doc='CVH-refit dimuon phi'),
            corMuPlusPt=ExtVar(cms.InputTag(_cor, 'muPlusPt'), float, doc='CVH-refit mu+ pt'),
            corMuMinusPt=ExtVar(cms.InputTag(_cor, 'muMinusPt'), float, doc='CVH-refit mu- pt'),
            # edmval < 0 marks a candidate whose dimuon leg did not converge:
            # the orphan flag that used to live in the offline join.
            corEdmval=ExtVar(cms.InputTag(_cor, 'edmval'), float, doc='CVH fit EDM (<0 = not refit)'),
            # Fitted mother candidate. Distinct cvh* names so these can never
            # be confused with BParking's own bkmm_jpsimc_* / bkmm_nomc_*.
            cvhFitMass=ExtVar(cms.InputTag('bplusFit', 'fitMass'), float, doc='fitted m(mumuK)'),
            cvhFitMassErr=ExtVar(cms.InputTag('bplusFit', 'fitMassErr'), float, doc='fitted mass error'),
            cvhFitPt=ExtVar(cms.InputTag('bplusFit', 'fitPt'), float, doc='fitted pt'),
            cvhFitVtxChi2=ExtVar(cms.InputTag('bplusFit', 'fitVtxChi2'), float, doc='fit vertex chi2'),
            cvhFitVtxProb=ExtVar(cms.InputTag('bplusFit', 'fitVtxProb'), float, doc='fit vertex prob'),
            cvhFitOk=ExtVar(cms.InputTag('bplusFit', 'fitOk'), int, doc='1 = kinematic fit succeeded'),
            # Cross-links into the Track table (-1 = no match). Enables e.g.
            # Track_dedxHarmonic2[BuJpsiK_kaonTrackIdx[i]] downstream.
            mu0TrackIdx=ExtVar(cms.InputTag('bplusLeafIdx', 'mu0TrackIdx'), int, doc='J/psi mu0 row in Track'),
            mu1TrackIdx=ExtVar(cms.InputTag('bplusLeafIdx', 'mu1TrackIdx'), int, doc='J/psi mu1 row in Track'),
            kaonTrackIdx=ExtVar(cms.InputTag('bplusLeafIdx', 'bach0TrackIdx'), int, doc='bachelor kaon row in Track'),
        ),
    )

    # Track -> Muon / Track -> PV cross-links, inverting the persisted
    # associations into row indices (-1 = none) on the Track table.
    process.trackMuonIdx = cms.EDProducer(
        'TrackToMuonIndexProducer',
        trackSrc=cms.InputTag(opts.srcTracks),
        association=cms.InputTag(opts.srcTracks + 'TrackToMuon'))
    process.trackPvIdx = cms.EDProducer(
        'TrackToVertexIndexProducer',
        trackSrc=cms.InputTag(opts.srcTracks),
        association=cms.InputTag('offlinePrimaryVertices'))

    process.trackTable = cms.EDProducer(
        'SimpleTrackFlatTableProducer',
        src=cms.InputTag(opts.srcTracks),
        cut=cms.string(''), name=cms.string('Track'),
        doc=cms.string('AlCaReco cloned alignment tracks'),
        singleton=cms.bool(False), extension=cms.bool(False),
        variables=cms.PSet(
            P3Vars,
            charge=Var('charge', 'int16', doc='charge'),
            dxy=Var('dxy', float, doc='dxy'), dz=Var('dz', float, doc='dz'),
            normChi2=Var('normalizedChi2', float, doc='chi2/ndof'),
            nValidHits=Var('numberOfValidHits', 'int16', doc='n valid hits'),
        ),
        externalVariables=cms.PSet(
            dedxHarmonic2=ExtVar(cms.InputTag(opts.srcTracks + 'DeDxHarmonic2'), float, doc='dE/dx harmonic2'),
            dedxPixelHarmonic2=ExtVar(cms.InputTag(opts.srcTracks + 'DeDxPixelHarmonic2'), float, doc='dE/dx pixel harmonic2'),
            originalIndex=ExtVar(cms.InputTag(opts.srcTracks, 'originalIndex'), 'uint', doc='index into generalTracks'),
            muonIdx=ExtVar(cms.InputTag('trackMuonIdx'), int, doc='row in Muon (-1 = none)'),
            pvIdx=ExtVar(cms.InputTag('trackPvIdx'), int, doc='row in PV (-1 = none)'),
        ),
    )

    process.muonTable = cms.EDProducer(
        'SimpleMuonFlatTableProducer',
        src=cms.InputTag(opts.srcMuons),
        cut=cms.string(''), name=cms.string('Muon'),
        doc=cms.string('AlCaReco persisted reco::Muon'),
        singleton=cms.bool(False), extension=cms.bool(False),
        variables=cms.PSet(
            P3Vars,
            charge=Var('charge', 'int16', doc='charge'),
            isGlobal=Var('isGlobalMuon', bool, doc='is global muon'),
            isTracker=Var('isTrackerMuon', bool, doc='is tracker muon'),
            nMatches=Var('numberOfMatches', 'int16', doc='n matched stations'),
        ),
    )

    process.pvTable = cms.EDProducer(
        'SimpleVertexFlatTableProducer',
        src=cms.InputTag('offlinePrimaryVertices'),
        cut=cms.string(''), name=cms.string('PV'),
        doc=cms.string('offline primary vertices'),
        singleton=cms.bool(False), extension=cms.bool(False),
        variables=cms.PSet(
            x=Var('x', float, doc='x'), y=Var('y', float, doc='y'), z=Var('z', float, doc='z'),
            chi2=Var('chi2', float, doc='chi2'), ndof=Var('ndof', float, doc='ndof'),
        ),
    )

    # Detector conditions: magnet current is what the CVH calibration wants.
    # DcsStatus is tiny (two floats + a 25-bit partition mask).
    process.dcsTable = cms.EDProducer(
        'SimpleDcsStatusFlatTableProducer',
        src=cms.InputTag('scalersRawToDigi'),
        cut=cms.string(''), name=cms.string('Dcs'),
        doc=cms.string('DcsStatus from scalersRawToDigi'),
        singleton=cms.bool(False), extension=cms.bool(False),
        variables=cms.PSet(
            magnetCurrent=Var('magnetCurrent', float, doc='magnet current [A]'),
            magnetTemperature=Var('magnetTemperature', float, doc='magnet temperature'),
            ready=Var('ready', 'uint', doc='per-partition ready bitmask'),
        ),
    )

    # Legacy L1: finalOR + packed 128-bit algo / 64-bit tech decision words.
    process.l1Table = cms.EDProducer(
        'L1LegacyDecisionTableProducer',
        src=cms.InputTag('gtDigis'),
        name=cms.string('L1'),
    )

    # Refit bachelor tracks, when the maker emits them. These are the inputs
    # a constrained B-vertex fit needs; refitOk flags legs that were not refit.
    # Fitted mother candidate (the Bmm5 chain). The stage-1 mass is a raw
    # four-vector sum; this is the actual vertex/kinematic fit.
    process.bplusFit = cms.EDProducer(
        'JpsiXKinematicFitProducer',
        src=_src_cands,
        srcTracks=cms.InputTag('globalCorJpsiKKaon', 'refit') if opts.emitRefitTracks
                  else cms.InputTag(''),
        srcRefitOk=cms.InputTag('globalCorJpsiKKaon', 'refitOk'),
        jpsiConstraint=cms.string(str(opts.jpsiConstraint)),
        jpsiMass=cms.double(3.0969),
        maxChi2=cms.double(-1.),
    )

    # Candidate-daughter -> Track row cross-links (flat-tree join keys).
    process.bplusLeafIdx = cms.EDProducer(
        'CandidateLeafTrackIndexProducer',
        src=_src_cands,
        trackSrc=cms.InputTag(opts.srcTracks),
    )

    _extra_tables = []
    if opts.emitRefitTracks:
        process.refitTrackTable = cms.EDProducer(
            'SimpleTrackFlatTableProducer',
            src=cms.InputTag('globalCorJpsiKKaon', 'refit'),
            cut=cms.string(''), name=cms.string('RefitTrack'),
            doc=cms.string('CVH-refit bachelor tracks (aligned with input tracks)'),
            singleton=cms.bool(False), extension=cms.bool(False),
            variables=cms.PSet(
                P3Vars,
                charge=Var('charge', 'int16', doc='charge'),
                dxy=Var('dxy', float, doc='dxy'), dz=Var('dz', float, doc='dz'),
                normChi2=Var('normalizedChi2', float, doc='chi2/ndof'),
                ptErr=Var('ptError', float, doc='pt uncertainty from the refit 5x5'),
            ),
            externalVariables=cms.PSet(
                refitOk=ExtVar(cms.InputTag('globalCorJpsiKKaon', 'refitOk'),
                               int, doc='1 = refit succeeded, 0 = input copy'),
            ),
        )
        _extra_tables.append(process.refitTrackTable)

    process.nanoTables = cms.Task(
        process.bplusTable, process.trackTable, process.muonTable,
        process.pvTable, process.dcsTable, process.l1Table, process.bplusFit,
        process.bplusLeafIdx, process.trackMuonIdx, process.trackPvIdx,
        *_extra_tables)
    process.nano_step = cms.Path(process.nanoTables)

    process.nanoOutput = cms.OutputModule(
        'NanoAODOutputModule',
        fileName=cms.untracked.string(opts.nanoOut if opts.nanoOut.startswith('file:')
                                      else 'file:' + opts.nanoOut),
        # NanoAODOutputModule turns edm::TriggerResults into HLT_* decision
        # branches on its own, so keeping the product is all the HLT bits need.
        outputCommands=cms.untracked.vstring(
            'drop *',
            'keep nanoaodFlatTable_*_*_*',
            'keep edmTriggerResults_*_*_*',
        ),
        compressionLevel=cms.untracked.int32(9),
        compressionAlgorithm=cms.untracked.string('LZMA'),
    )
    process.nano_out_step = cms.EndPath(process.nanoOutput)
process.schedule = cms.Schedule(process.reconstruction_step)
if opts.nanoOut:
    process.schedule.extend([process.nano_step, process.nano_out_step])

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
