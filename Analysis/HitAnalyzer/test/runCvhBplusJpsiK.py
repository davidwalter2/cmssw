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
              'makers share the one CvhMaster ES product), "dimuon", "kaon", or '
              '"joint" (the CALIBRATION path: preselector + joint N-body maker '
              'writing a grads sidecar, no nano, no single-track makers)')
opts.register('selTight', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'mode=joint only: use the Bmm5 ANALYSIS-level kinematic cuts in '
              'the preselector instead of the loose nominal. This is the '
              'labelled systematic variation -- it keeps 18%% of bachelors '
              '(median bachelor pT ~0.53 GeV) and so removes exactly the soft '
              'phase space the channel exists to probe.')
opts.register('motherConstraintWidth', -1.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'sigma of the MOTHER (B) mass constraint, GeV. Negative keeps the '
              'maker default of 1e-3. NOT the B natural width (~4e-13): the '
              'constraint compares a MODEL PREDICTION to the true mass, so its '
              'width must cover the model error, not just the true spread.')
opts.register('gradsPass', 'allcons', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'mode=joint only: which constraint pass the stored grad/Hess come '
              'from. "allcons" (default, J/psi AND B mass constrained) is the '
              'headline arm and needs no maker change -- the extraction block '
              'sits after the icons loop, so it already stores the last pass. '
              '"subcons" (J/psi only) is the null control and truncates the '
              'pass list. "free" drops all mass rows.')
opts.register('emitRefitTracks', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'single-track maker emits refit reco::Tracks (+refitOk map) for '
              'the downstream constrained B-vertex fit')
opts.register('jointCvh', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'schedule the joint N-body CVH arm (jointCvh* columns). '
              'Measured: adds ~0.87 s/event to a ~3.4 s/event job, i.e. '
              'about +50%%. False drops the arm and its columns entirely.')
opts.register('twoTrackArm', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'schedule the TWO-TRACK dimuon CVH maker (globalCorJpsiK) and '
              'its cor* columns -- the J/psi-mass-constrained dimuon refit. '
              'This is the single biggest cost in the job: it issues ~254 '
              'Geant4e propagations per fit against ~31 for a single track, '
              'about 89%% of the propagation work. Production sets False: the '
              'three arms that are wanted are kvfRaw*, kvfCvh* and jointCvh*, '
              'and the joint maker applies its own J/psi constraint anyway. '
              'True (default) keeps the historical four-arm output for A/B '
              'work.')
opts.register('emitRefJacobian', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'store the reference-point jacobian d(track state)/d(global '
              'params) and the global-parameter indices, flattened into '
              'companion tables. Lets updated global corrections be applied '
              'downstream by a linearised update WITHOUT reproducing the '
              'NanoAOD (AN2021_131 sec. 420-433). This is NOT the calibration '
              'payload -- no Hessian, no mass jacobians; those stay behind '
              'fillGradsFactored. Costs ~3-4x the nano size.')
opts.register('timing', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'enable the Timing service with per-module breakdown for '
              'resource cost. Off in production -- it prints per module.')
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
              'A/B test of the kaon q/p anomaly)')
opts.register('kaonAsMuon', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'when True, override `trackParticleName` on the kaon-side maker '
              'to "mu" (mass-hypothesis A/B test). Inputs remain the '
              'bachelor kaon tracks, only the propagation hypothesis changes. '
              'NOTE: the '
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
              'The cfi value was 1.0; this default ramps it down to 0.05 to recover the '
              'soft-bachelor tail. Lowering it on the dimuon side is a no-op '
              '(muons clear 1.0 trivially); the knob affects the kaon mode.')
# CVH joint-refit convergence knobs
# (default values reproduce the published baseline bit-identically).
opts.register('clampMomentumFloor', -1.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'Gauss-Newton momentum floor [GeV] of the refit step clamp. '
              '<0 (default) = derive it from plimit as 1.25*plimit. THE FLOOR '
              'MUST STAY ABOVE THE PROPAGATION LIMIT (the clamp exists only to '
              'keep the Gauss-Newton state out of the propagator refusal '
              'region) and NOT ABOVE THE PHYSICAL MOMENTUM SPECTRUM: a floor '
              'above it pins soft tracks at the floor and, where p_ref is '
              'already below it, scales the step to zero and freezes the fit '
              'at its seed. The built-in 2.0 GeV default of the maker would '
              'do exactly that to the 0.3-0.9 GeV bachelor kaon, which is why '
              'this driver always sets the floor from plimit.')
opts.register('maxMomentumStepFactor', 2.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'RELATIVE Gauss-Newton step damping: the max factor by which a '
              'track momentum may change in one iteration (default 2 = p may '
              'at most halve or double). Effective floor '
              'max(clampMomentumFloor, p_ref/f) plus the symmetric cap '
              'p_ref*f, so the bound is ALWAYS strictly inside p_ref: neither '
              'the soft bachelor kaon nor any other leg can be pinned at a '
              'fixed momentum or frozen by a zero step. <=1 switches it off '
              'and leaves clampMomentumFloor as the only bound.')
opts.register('stepBacktracking', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'chi2-based (Armijo) retroactive step backtracking: if the chi2 '
              'realized by the previous step fails the sufficient-decrease '
              'test, restore that linearization, halve the step and redo the '
              'iteration. No extra propagation on the accept path.')
opts.register('stepBacktrackFromIter', 2, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'first Gauss-Newton iteration at which the chi2 backtracking test '
              'may fire (default 2). NOT 1: at iteration 0 the GBL '
              'propagation/kink residuals are identically zero by construction, '
              'so the iteration-0 chi2 is a different objective from every later '
              'one and comparing across that boundary would backtrack every '
              'candidate.')
opts.register('maxChi2Backtrack', 4, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'max chi2-backtracking halvings per accepted step (default 4). '
              'Distinct from maxBacktracks, the failed-leg retry budget.')
opts.register('armijoC', 1.e-4, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'Armijo sufficient-decrease coefficient c1 (default 1e-4).')
opts.register('armijoSlack', 1.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'relative chi2 slack in the Armijo test (default 1.0 = the chi2 may not more than DOUBLE in one iteration). This is a DIVERGENCE TRAP, not a line-search tolerance: the CVH/GBL iteration does not monotonically decrease r^T Vinv r, so a textbook 1e-4..1e-3 makes it halve the step forever at 6x the propagation cost for no change in the result.')
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
# Global-correction calibration output. Defaults are all OFF/empty so the diagnostic offline-join behaviour
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
opts.register('isMC', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'MC input: emit a Gen table + generator weight and gen-match the B '
              'candidate (genB* columns). Off for data (schema unchanged). When '
              'True, also pass the MC globalTag, e.g. 106X_mcRun2_asymptotic_v17.')
opts.register('useIdealGeometry', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'Default False (btojpsik option (B)): propagate against the aligned geometry '
              'loaded from the CMSSW GlobalTag. Set True for the AN2021_131_v8 canonical '
              'config, in which case populate corFiles (below) with a correctionResults*.root '
              'whose parmset matches this producer build; the WMass-era v721 file does NOT '
              'match this build (idxmaptree remapping not implemented in the current maker), '
              'so True + corFiles=[v721] currently crashes at the parmset-size assert. '
              'The (A) path (idxmaptree remapping) is a follow-up.')
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
assert opts.mode in ('both', 'dimuon', 'kaon', 'joint'), \
    f'mode must be both|dimuon|kaon|joint, got {opts.mode!r}'
assert opts.gradsPass in ('allcons', 'subcons', 'free'), \
    f'gradsPass must be allcons|subcons|free, got {opts.gradsPass!r}'
if opts.mode == 'joint':
    # The calibration path writes a grads sidecar, not a nano: the joint maker's
    # tree IS the output. Forcing nanoOut off here rather than asking the caller
    # to remember keeps a grads production from silently also paying for (and
    # writing) a NanoAOD.
    if opts.nanoOut:
        print('[runCvhBplusJpsiK] mode=joint: forcing nanoOut=False '
              '(the grads sidecar is the output)')
        opts.nanoOut = False
    if opts.emitRefitTracks:
        raise ValueError('mode=joint does not use the single-track makers; '
                         'emitRefitTracks must be False')
    # NB: _uses_splitter is defined further down, so test opts directly here.
    if 'jpsiKCandidateSplitter' in opts.srcCandidates:
        raise ValueError('mode=joint reads the nested candidate directly; '
                         'srcCandidates must not be a jpsiKCandidateSplitter '
                         'instance')

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
if opts.timing:
    # Per-module CPU/real time. summaryOnly=False gives the per-module table,
    # which is what separates the joint CVH arm's cost from the rest.
    process.Timing = cms.Service(
        'Timing',
        summaryOnly=cms.untracked.bool(False),
        excessiveTimeThreshold=cms.untracked.double(0.),
        useJobReport=cms.untracked.bool(True))

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
        'CandidateLeafTrackProducer', src=cms.InputTag(opts.srcCandidates),
        mode=cms.string('bachelor'))
    _bach_src = cms.InputTag('bplusBachelorTracks')
    _bach_idx = cms.InputTag('bplusBachelorTracks', 'candIdx')

# The refit ARM of the B fit also needs the muon legs refit (add-cvh-refit-tracks-
# into-bfit). The muons are the leaves of daughter(0) (the J/psi), so a second
# CandidateLeafTrackProducer in "jpsi" mode extracts them; a cloned single-track
# maker (trackParticleName='mu') refits them. Leaf-keyed via candIdx/leafIdx.
# Requires the non-splitter path (the splitter emits no leafIdx).
if opts.emitRefitTracks and _uses_splitter:
    raise ValueError('emitRefitTracks=True requires the non-splitter path '
                     '(srcCandidates must not be jpsiKCandidateSplitter:*): the '
                     'refit-arm mapping needs CandidateLeafTrackProducer leafIdx')
if opts.emitRefitTracks:
    process.bplusJpsiMuonTracks = cms.EDProducer(
        'CandidateLeafTrackProducer', src=cms.InputTag(opts.srcCandidates),
        mode=cms.string('jpsi'))

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
    # Mass-hypothesis knob: override the propagation hypothesis from "kaon" to "mu"
    # while leaving the input track collection (bachelor kaons) untouched.
    trackParticleName=cms.string('mu' if opts.kaonAsMuon else 'kaon'),
    emitRefitTracks=cms.bool(bool(opts.emitRefitTracks)),
    refitMaxRelPtErr=cms.double(float(opts.refitMaxRelPtErr)),
    # Track-keyed reference-point jacobian (5 x nParms) + global indices. The
    # muon-keyed jacRef map this maker has always produced is written inside an
    # `if (muonref.isNonnull())` gate and is therefore EMPTY here -- the
    # bachelor and J/psi-leg makers run on bare track collections with no muon
    # association at all.
    emitRefJacobian=cms.bool(bool(opts.emitRefJacobian)),
    # Muon-keyed jacRef rows: the full 5-row reference state whenever the
    # joint N-body arm or the jacobian export runs (displacement cuts need
    # dxy/dsz); 3 momentum rows otherwise. Only matters with muon association.
    jacRefRows=cms.int32(5 if (opts.jointCvh or opts.emitRefJacobian) else 3),
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

# ---- single-track maker (J/psi muons, refit arm) --------------------------
# Clone of the single-track kaon maker but fed the two muon legs (from the
# "jpsi"-mode CandidateLeafTrackProducer) and propagated with the muon
# hypothesis (trackParticleName='mu'). Emits refit reco::Tracks + refitOk for
# the constrained B-vertex fit's REFIT arm. Only scheduled when emitRefitTracks.
if opts.emitRefitTracks:
    process.globalCorJpsiKMuon = globalCorJpsiKKaon.clone(
        src=cms.InputTag('bplusJpsiMuonTracks'),
        bCandIdxSrc=cms.InputTag('bplusJpsiMuonTracks', 'candIdx'),
        scalarPotentialInitFile=cms.string(opts.scalarPot3DInitFile),
        fillJac=cms.bool(bool(opts.fillJac)),
        runFDClosure=cms.bool(bool(opts.runFDClosure)),
        epsilonFDClosure=cms.double(float(opts.epsilonFDClosure)),
        useIdealGeometry=cms.bool(_muon_ideal_geom),
        corFiles=cms.vstring(*_muon_corfiles),
        trackParticleName=cms.string('mu'),
        emitRefitTracks=cms.bool(True),
        refitMaxRelPtErr=cms.double(float(opts.refitMaxRelPtErr)),
        emitRefJacobian=cms.bool(bool(opts.emitRefJacobian)),
        # Muon-keyed jacRef rows: the full 5-row reference state whenever the
        # joint N-body arm or the jacobian export runs (displacement cuts need
        # dxy/dsz); 3 momentum rows otherwise. Only matters with muon association.
        jacRefRows=cms.int32(5 if (opts.jointCvh or opts.emitRefJacobian) else 3),
        CvhMaster=CvhMasterPSet.clone(
            Particles=cms.vstring('mu+', 'mu-', 'kaon+', 'kaon-')),
        **_calib_pset,
    )
    process.RandomNumberGeneratorService.globalCorJpsiKMuon = cms.PSet(
        initialSeed=cms.untracked.uint32(345678901),
        engineName=cms.untracked.string('HepJamesRandom'),
    )

# ---- magnetic field --------------------------------------------------------
# Default: load the scalar-potential 3D field producer with a unique label and
# point Geant4ePropagator + makers at it.
# `useScalarPot3D=False` (field A/B test) skips the producer entirely;
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
    _field_makers = [process.globalCorJpsiKKaon]
    if opts.twoTrackArm:
        _field_makers.append(process.globalCorJpsiK)
    if opts.emitRefitTracks:
        _field_makers.append(process.globalCorJpsiKMuon)
    for m in _field_makers:
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

# Gauss-Newton momentum floor of every residual maker scheduled in this
# process, derived from the propagation limit unless given explicitly. Applied
# generically (the maker set depends on `mode`), after all modules exist.
_clampFloor = (float(opts.clampMomentumFloor) if float(opts.clampMomentumFloor) > 0.
               else 1.25 * float(opts.plimit))
if _clampFloor <= float(opts.plimit):
    raise RuntimeError(
        "clampMomentumFloor (%g GeV) must be ABOVE PropagationPtotLimit (%g GeV)"
        % (_clampFloor, float(opts.plimit)))

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
# The refit arm needs BOTH single-track makers (kaon + muon), so emitRefitTracks
# implies mode must schedule the kaon maker (mode in both|kaon). The default
# mode="both" satisfies this; guard the pathological combination explicitly.
if opts.emitRefitTracks and opts.mode not in ('both', 'kaon'):
    raise ValueError('emitRefitTracks=True needs the single-track makers; use '
                     'mode="both" (default) or "kaon", not mode="%s"' % opts.mode)
_seq = process.offlineBeamSpot
if _uses_splitter:
    _seq = _seq * process.jpsiKCandidateSplitter
elif opts.mode in ('both', 'kaon'):
    _seq = _seq * process.bplusBachelorTracks
if opts.emitRefitTracks:
    _seq = _seq * process.bplusJpsiMuonTracks
# twoTrackArm=False drops the dimuon TWO-TRACK maker. It is ~89% of the
# Geant4e propagation work, and the three arms production wants (kvfRaw*,
# kvfCvh*, jointCvh*) do not consume it -- the joint maker applies its own
# J/psi constraint, and the KVF arm applies jpsiConstraint=inFit.
_two_track = opts.twoTrackArm and opts.mode in ('both', 'dimuon')
if opts.mode == 'dimuon' and not opts.twoTrackArm:
    raise ValueError('mode="dimuon" with twoTrackArm=False schedules no maker '
                     'at all -- pick one.')
if _two_track:
    _seq = _seq * process.globalCorJpsiK
if opts.mode in ('both', 'kaon'):
    _seq = _seq * process.globalCorJpsiKKaon
if opts.emitRefitTracks:
    _seq = _seq * process.globalCorJpsiKMuon

# mode="joint" -- the CALIBRATION path.
#
# Deliberately minimal: beamspot -> raw KVF arm -> preselector -> joint N-body
# maker. Three things it does NOT schedule, each for a reason:
#   * no single-track makers -- refitLegs is empty, so the KVF arm runs on the
#     stage-1 tracks only. The preselector must cut on RAW quantities (a cut on
#     CVH-refit quantities would make the selection a function of the
#     corrections being fitted), so the refit arm has no consumer here.
#   * no CandidateVertexGeometryProducer -- the nominal selection uses the
#     DIMUON flight significance (bplusFit:rawDimuonSl3d), not the mother's, so
#     the vertex-geometry block has no consumer either.
#   * no nano -- the maker's grads tree is the output (forced above).
if opts.mode == 'joint':
    process.bplusFit = cms.EDProducer(
        'JpsiXKinematicFitProducer',
        src=_src_cands,
        refitLegs=cms.VPSet(),          # raw arm only, by design (see above)
        jpsiConstraint=cms.string(str(opts.jpsiConstraint)),
        jpsiMass=cms.double(3.0969),
        maxChi2=cms.double(-1.),
        beamSpot=cms.InputTag('offlineBeamSpot'),
        primaryVertices=cms.InputTag('offlinePrimaryVertices'),
        genParticles=cms.InputTag('genParticles') if opts.isMC else cms.InputTag(''),
    )

    from Analysis.HitAnalyzer.jpsiXCandidatePreselectorForCorrections_cfi import (
        jpsiXCandidatePreselectorForCorrections, jpsiXCandidatePreselectorForCorrectionsTight)
    _presel = jpsiXCandidatePreselectorForCorrectionsTight if opts.selTight \
        else jpsiXCandidatePreselectorForCorrections
    process.bplusPreselect = _presel.clone(src=_src_cands)

    # The joint maker runs over the SURVIVORS, not the parent collection.
    _sel_cands = cms.InputTag('bplusPreselect')
    process.jointCvhGrads = globalCorJpsiK.clone(
        src=cms.InputTag(opts.srcTracks),
        dedxSourceTracks=cms.InputTag(opts.srcTracks),
        srcCandidates=_sel_cands,
        bCandIdxSrc=cms.InputTag(''),
        useIdealGeometry=cms.bool(False),
        corFiles=cms.vstring(),
        # The grads tree is the whole point, so unlike the nano-path clone this
        # one has fillTrackTree ON and produceValueMaps OFF.
        fillTrackTree=cms.bool(True),
        produceValueMaps=cms.bool(False),
        fillJac=cms.bool(bool(opts.fillJac)),
        gradsPass=cms.string(str(opts.gradsPass)),
        outprefix=cms.untracked.string('globalcor_jpsik_joint'),
        scalarPotentialInitFile=cms.string(opts.scalarPot3DInitFile),
        nIters=cms.uint32(int(opts.nIters)),
        edmConvergence=cms.double(float(opts.edmConvergence)),
        useStartingState=cms.string(str(opts.useStartingState)),
        debugPerIterDump=cms.bool(bool(opts.debug)),
        CvhMaster=CvhMasterPSet.clone(
            Particles=cms.vstring('mu+', 'mu-', 'kaon+', 'kaon-')),
        **_calib_pset,
    )
    if opts.motherConstraintWidth > 0.:
        process.jointCvhGrads.motherConstraintWidth = \
            cms.double(float(opts.motherConstraintWidth))
    process.jointCvhGrads.__dict__['_TypedParameterizable__type'] = \
        'ResidualGlobalCorrectionMakerNTrackG4e'
    process.RandomNumberGeneratorService.jointCvhGrads = cms.PSet(
        initialSeed=cms.untracked.uint32(123456789),
        engineName=cms.untracked.string('HepJamesRandom'))

    _seq = _seq * process.bplusFit * process.bplusPreselect * process.jointCvhGrads
    print('[runCvhBplusJpsiK] mode=joint: gradsPass={} selection={}'.format(
        opts.gradsPass, 'TIGHT (analysis, systematic)' if opts.selTight
        else 'loose nominal'), flush=True)

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
            # ---- Fitted mother, three arms (naming per proposal 2026-07-28) ----
            # METHOD first, track source second:
            #   kvfRaw*    KVF on the raw stage-1 tracks
            #   kvfCvh*    KVF on the CVH single-track refit  ("route 1")
            #   jointCvh*  joint N-body CVH -- fit and vertex in one
            # The old cvhFit*/dimuon* UNSUFFIXED aliases are deliberately gone:
            # a column meaning "whichever arm is currently default" does not
            # survive the decision that set it. `cvhFit*` was doubly bad -- it
            # read as "the CVH fit" but meant "a KVF on CVH-refit tracks".
            kvfRawMass=ExtVar(cms.InputTag('bplusFit', 'rawFitMass'), float, doc='m(mumuK), KVF on raw tracks'),
            kvfRawMassErr=ExtVar(cms.InputTag('bplusFit', 'rawFitMassErr'), float, doc='mass error, KVF on raw tracks'),
            kvfRawPt=ExtVar(cms.InputTag('bplusFit', 'rawFitPt'), float, doc='pt, KVF on raw tracks'),
            kvfRawVtxChi2=ExtVar(cms.InputTag('bplusFit', 'rawFitVtxChi2'), float, doc='vertex chi2, KVF on raw tracks'),
            kvfRawVtxProb=ExtVar(cms.InputTag('bplusFit', 'rawFitVtxProb'), float, doc='vertex prob, KVF on raw tracks'),
            kvfRawOk=ExtVar(cms.InputTag('bplusFit', 'rawFitOk'), int, doc='1 = KVF on raw tracks succeeded'),

            kvfCvhMass=ExtVar(cms.InputTag('bplusFit', 'refFitMass'), float, doc='m(mumuK), KVF on CVH-refit tracks'),
            kvfCvhMassErr=ExtVar(cms.InputTag('bplusFit', 'refFitMassErr'), float, doc='mass error, KVF on CVH-refit tracks'),
            kvfCvhPt=ExtVar(cms.InputTag('bplusFit', 'refFitPt'), float, doc='pt, KVF on CVH-refit tracks'),
            kvfCvhVtxChi2=ExtVar(cms.InputTag('bplusFit', 'refFitVtxChi2'), float, doc='vertex chi2, KVF on CVH-refit tracks'),
            kvfCvhVtxProb=ExtVar(cms.InputTag('bplusFit', 'refFitVtxProb'), float, doc='vertex prob, KVF on CVH-refit tracks'),
            kvfCvhOk=ExtVar(cms.InputTag('bplusFit', 'refFitOk'), int, doc='1 = KVF on CVH-refit tracks succeeded'),
            nLegsRefit=ExtVar(cms.InputTag('bplusFit', 'nLegsRefit'), int, doc='number of legs (0-3) using a CVH refit track'),

            # ---- Dimuon (J/psi) VERTEX quantities -----------------------
            # These come from route 1's SEPARATE dimuon fit. The joint arm has
            # no counterpart: all N tracks share ONE vertex by construction, so
            # there is no distinct dimuon vertex to characterise.
            kvfRawDimuonVtxProb=ExtVar(cms.InputTag('bplusFit', 'rawDimuonVtxProb'), float, doc='dimuon vertex prob, raw tracks'),
            kvfRawDimuonAlphaBS=ExtVar(cms.InputTag('bplusFit', 'rawDimuonAlphaBS'), float, doc='dimuon XY pointing angle wrt BS, raw tracks'),
            kvfRawDimuonSxy=ExtVar(cms.InputTag('bplusFit', 'rawDimuonSxy'), float, doc='dimuon 2D Lxy significance wrt BS, raw tracks'),
            kvfRawDimuonSl3d=ExtVar(cms.InputTag('bplusFit', 'rawDimuonSl3d'), float, doc='dimuon 3D flight significance wrt closest-z PV, raw tracks'),
            kvfCvhDimuonVtxProb=ExtVar(cms.InputTag('bplusFit', 'refDimuonVtxProb'), float, doc='dimuon vertex prob, CVH-refit tracks'),
            kvfCvhDimuonAlphaBS=ExtVar(cms.InputTag('bplusFit', 'refDimuonAlphaBS'), float, doc='dimuon XY pointing angle wrt BS, CVH-refit tracks'),
            kvfCvhDimuonSxy=ExtVar(cms.InputTag('bplusFit', 'refDimuonSxy'), float, doc='dimuon 2D Lxy significance wrt BS, CVH-refit tracks'),
            kvfCvhDimuonSl3d=ExtVar(cms.InputTag('bplusFit', 'refDimuonSl3d'), float, doc='dimuon 3D flight significance wrt closest-z PV, CVH-refit tracks'),
            # UNCONSTRAINED dimuon mass (gap G11). The dimuon fit above applies
            # no mass constraint, so this is a free J/psi mass -- the natural
            # cross-check on the muon momentum scale, and the quantity tying
            # this channel to the J/psi calibration. Not to be confused with
            # jointCvhJpsiMass, which is the joint fit's CONSTRAINED parameter
            # and is pinned to the PDG value for every candidate.
            kvfRawDimuonMass=ExtVar(cms.InputTag('bplusFit', 'rawDimuonMass'), float, doc='UNCONSTRAINED dimuon mass, raw tracks'),
            kvfRawDimuonMassErr=ExtVar(cms.InputTag('bplusFit', 'rawDimuonMassErr'), float, doc='uncertainty on the unconstrained dimuon mass, raw tracks'),
            kvfCvhDimuonMass=ExtVar(cms.InputTag('bplusFit', 'refDimuonMass'), float, doc='UNCONSTRAINED dimuon mass, CVH-refit tracks'),
            kvfCvhDimuonMassErr=ExtVar(cms.InputTag('bplusFit', 'refDimuonMassErr'), float, doc='uncertainty on the unconstrained dimuon mass, CVH-refit tracks'),

            # ---- B-vertex geometry, one geometry module per arm ----------
            # Same implementation for every arm: pure
            # arithmetic on each arm's fitted vertex block. No fit happens here.
            kvfRawAlphaBS=ExtVar(cms.InputTag('vtxGeomKvfRaw', 'alphaBS'), float, doc='B XY pointing angle wrt BS, KVF raw'),
            kvfRawLxy=ExtVar(cms.InputTag('vtxGeomKvfRaw', 'lxy'), float, doc='B transverse flight from BS, KVF raw'),
            kvfRawSxy=ExtVar(cms.InputTag('vtxGeomKvfRaw', 'sxy'), float, doc='B Lxy significance, KVF raw'),
            kvfRawL3d=ExtVar(cms.InputTag('vtxGeomKvfRaw', 'l3d'), float, doc='B 3D flight from closest-z PV, KVF raw'),
            kvfRawSl3d=ExtVar(cms.InputTag('vtxGeomKvfRaw', 'sl3d'), float, doc='B 3D flight significance, KVF raw'),
            kvfRawAlpha3dPV=ExtVar(cms.InputTag('vtxGeomKvfRaw', 'alpha3dPV'), float, doc='B 3D pointing angle wrt PV, KVF raw'),
            kvfCvhAlphaBS=ExtVar(cms.InputTag('vtxGeomKvfCvh', 'alphaBS'), float, doc='B XY pointing angle wrt BS, KVF CVH'),
            kvfCvhLxy=ExtVar(cms.InputTag('vtxGeomKvfCvh', 'lxy'), float, doc='B transverse flight from BS, KVF CVH'),
            kvfCvhSxy=ExtVar(cms.InputTag('vtxGeomKvfCvh', 'sxy'), float, doc='B Lxy significance, KVF CVH'),
            kvfCvhL3d=ExtVar(cms.InputTag('vtxGeomKvfCvh', 'l3d'), float, doc='B 3D flight from closest-z PV, KVF CVH'),
            kvfCvhSl3d=ExtVar(cms.InputTag('vtxGeomKvfCvh', 'sl3d'), float, doc='B 3D flight significance, KVF CVH'),
            kvfCvhAlpha3dPV=ExtVar(cms.InputTag('vtxGeomKvfCvh', 'alpha3dPV'), float, doc='B 3D pointing angle wrt PV, KVF CVH'),
            # Cross-links into the Track table (-1 = no match). Enables e.g.
            # Track_dedxHarmonic2[BuJpsiK_kaonTrackIdx[i]] downstream.
            mu0TrackIdx=ExtVar(cms.InputTag('bplusLeafIdx', 'mu0TrackIdx'), int, doc='J/psi mu0 row in Track'),
            mu1TrackIdx=ExtVar(cms.InputTag('bplusLeafIdx', 'mu1TrackIdx'), int, doc='J/psi mu1 row in Track'),
            kaonTrackIdx=ExtVar(cms.InputTag('bplusLeafIdx', 'bach0TrackIdx'), int, doc='bachelor kaon row in Track'),
        ),
    )

    # Two-track dimuon CVH arm (cor*). Added here rather than inline so the arm
    # can be switched off without leaving dangling InputTags -- production runs
    # with twoTrackArm=False, and an ExtVar pointing at an unscheduled module
    # is a configuration error, not a missing column.
    if opts.twoTrackArm:
        _ev = process.bplusTable.externalVariables
        _ev.corMass = ExtVar(cms.InputTag(_cor, 'corMass'), float, doc='CVH-refit dimuon mass')
        _ev.corMassErr = ExtVar(cms.InputTag(_cor, 'corMassErr'), float, doc='CVH-refit mass error')
        _ev.corPt = ExtVar(cms.InputTag(_cor, 'corPt'), float, doc='CVH-refit dimuon pt')
        _ev.corEta = ExtVar(cms.InputTag(_cor, 'corEta'), float, doc='CVH-refit dimuon eta')
        _ev.corPhi = ExtVar(cms.InputTag(_cor, 'corPhi'), float, doc='CVH-refit dimuon phi')
        _ev.corMuPlusPt = ExtVar(cms.InputTag(_cor, 'muPlusPt'), float, doc='CVH-refit mu+ pt')
        _ev.corMuMinusPt = ExtVar(cms.InputTag(_cor, 'muMinusPt'), float, doc='CVH-refit mu- pt')
        # edmval < 0 marks a candidate whose dimuon leg did not converge:
        # the orphan flag that used to live in the offline join.
        _ev.corEdmval = ExtVar(cms.InputTag(_cor, 'edmval'), float, doc='CVH fit EDM (<0 = not refit)')

    # Joint-CVH arm columns. Added here rather than inline so the arm can be
    # switched off without leaving dangling InputTags.
    if opts.jointCvh:
        _j, _jg = 'jointCvhBu', 'vtxGeomJointCvh'
        _ev = process.bplusTable.externalVariables
        _ev.jointCvhMass = ExtVar(cms.InputTag(_j, 'motherMass'), float, doc='m(mumuK), joint N-body CVH (subcons pass)')
        _ev.jointCvhMassErr = ExtVar(cms.InputTag(_j, 'motherMassErr'), float, doc='mass error, joint N-body CVH')
        _ev.jointCvhPt = ExtVar(cms.InputTag(_j, 'motherPt'), float, doc='pt, joint N-body CVH')
        _ev.jointCvhEta = ExtVar(cms.InputTag(_j, 'motherEta'), float, doc='eta, joint N-body CVH')
        _ev.jointCvhPhi = ExtVar(cms.InputTag(_j, 'motherPhi'), float, doc='phi, joint N-body CVH')
        # allcons pass: the B mass is FORCED to PDG by construction. Never use
        # this as a mass measurement -- it is a constraint, not a result.
        _ev.jointCvhConsMass = ExtVar(cms.InputTag(_j, 'motherConsMass'), float, doc='m(mumuK) with the mother mass constraint APPLIED (not a measurement)')
        _ev.jointCvhChisq = ExtVar(cms.InputTag(_j, 'chisq'), float, doc='joint CVH fit chi2')
        _ev.jointCvhNdof = ExtVar(cms.InputTag(_j, 'ndof'), float, doc='joint CVH fit ndof')
        # REFERENCE-BLOCK edm. The full-state edmval is re-zeroed each
        # re-linearisation and says nothing about convergence.
        _ev.jointCvhEdmRef = ExtVar(cms.InputTag(_j, 'edmvalRef'), float, doc='joint CVH reference-block EDM (converged: < 1e-5)')
        _ev.jointCvhNiter = ExtVar(cms.InputTag(_j, 'niter'), int, doc='joint CVH iterations')
        # 1 = the fit COMPLETED. NOT a physics-quality flag:
        # a runaway fit can complete. Cut on chisq/ndof and edmRef.
        _ev.jointCvhOk = ExtVar(cms.InputTag(_j, 'fitOk'), int, doc='1 = joint CVH fit completed (NOT a quality flag)')
        # CONSTRAINED-SUBSYSTEM (J/psi) quantities from the SAME joint fit that
        # produced the mother. The maker computes these already and puts them
        # on its corMass/corPt/... maps -- they were simply never consumed,
        # because the cor* columns were wired to the two-track maker instead.
        # Summed over massconstrainttracks, so this is the J/psi for a B+ and
        # whatever subsystem carries the constraint for another channel.
        # Unlike jointCvhConsMass, this is a MEASUREMENT: it comes from the
        # unconstrained (icons==0) pass, so the J/psi mass is not forced.
        _ev.jointCvhJpsiMass = ExtVar(cms.InputTag(_j, 'corMass'), float, doc='m(mumu) from the joint N-body CVH fit (unconstrained pass)')
        _ev.jointCvhJpsiMassErr = ExtVar(cms.InputTag(_j, 'corMassErr'), float, doc='mass error on m(mumu), joint N-body CVH')
        _ev.jointCvhJpsiPt = ExtVar(cms.InputTag(_j, 'corPt'), float, doc='pt of the mumu subsystem, joint N-body CVH')
        _ev.jointCvhJpsiEta = ExtVar(cms.InputTag(_j, 'corEta'), float, doc='eta of the mumu subsystem, joint N-body CVH')
        _ev.jointCvhJpsiPhi = ExtVar(cms.InputTag(_j, 'corPhi'), float, doc='phi of the mumu subsystem, joint N-body CVH')
        _ev.jointCvhVtxX = ExtVar(cms.InputTag(_j, 'vtxX'), float, doc='joint CVH common vertex x')
        _ev.jointCvhVtxY = ExtVar(cms.InputTag(_j, 'vtxY'), float, doc='joint CVH common vertex y')
        _ev.jointCvhVtxZ = ExtVar(cms.InputTag(_j, 'vtxZ'), float, doc='joint CVH common vertex z')
        _ev.jointCvhAlphaBS = ExtVar(cms.InputTag(_jg, 'alphaBS'), float, doc='B XY pointing angle wrt BS, joint CVH')
        _ev.jointCvhLxy = ExtVar(cms.InputTag(_jg, 'lxy'), float, doc='B transverse flight from BS, joint CVH')
        _ev.jointCvhSxy = ExtVar(cms.InputTag(_jg, 'sxy'), float, doc='B Lxy significance, joint CVH')
        _ev.jointCvhL3d = ExtVar(cms.InputTag(_jg, 'l3d'), float, doc='B 3D flight from closest-z PV, joint CVH')
        _ev.jointCvhSl3d = ExtVar(cms.InputTag(_jg, 'sl3d'), float, doc='B 3D flight significance, joint CVH')
        _ev.jointCvhAlpha3dPV = ExtVar(cms.InputTag(_jg, 'alpha3dPV'), float, doc='B 3D pointing angle wrt PV, joint CVH')
        # Per-leg FITTED momenta, decomposition leaf order (mu-, mu+, K for a
        # B+). Fixed-index scalar columns because a nano flat table cannot hold
        # a variable-length ValueMap<vector<float>>; absent legs are -99.
        # Leaf order matches the mu0/mu1/kaonTrackIdx cross-links above.
        for _i in range(3):
            setattr(_ev, f'jointCvhLeg{_i}Pt',
                    ExtVar(cms.InputTag(_j, f'leg{_i}Pt'), float, doc=f'leg {_i} fitted pt, joint CVH'))
            setattr(_ev, f'jointCvhLeg{_i}Eta',
                    ExtVar(cms.InputTag(_j, f'leg{_i}Eta'), float, doc=f'leg {_i} fitted eta, joint CVH'))
            setattr(_ev, f'jointCvhLeg{_i}Phi',
                    ExtVar(cms.InputTag(_j, f'leg{_i}Phi'), float, doc=f'leg {_i} fitted phi, joint CVH'))

    # Gen-truth columns (MC only): the closest last-copy b-hadron to the
    # candidate. genBPdgId ~ +-521 with genBDR < ~0.1 tags a true B+ -> J/psi K.
    if opts.isMC:
        _gf = 'bplusFit'
        process.bplusTable.externalVariables.genBMass = ExtVar(cms.InputTag(_gf, 'genBMass'), float, doc='matched gen b-hadron mass')
        process.bplusTable.externalVariables.genBPt = ExtVar(cms.InputTag(_gf, 'genBPt'), float, doc='matched gen b-hadron pt')
        process.bplusTable.externalVariables.genBEta = ExtVar(cms.InputTag(_gf, 'genBEta'), float, doc='matched gen b-hadron eta')
        process.bplusTable.externalVariables.genBPhi = ExtVar(cms.InputTag(_gf, 'genBPhi'), float, doc='matched gen b-hadron phi')
        process.bplusTable.externalVariables.genBDR = ExtVar(cms.InputTag(_gf, 'genBDR'), float, doc='dR(candidate, matched gen b-hadron)')
        process.bplusTable.externalVariables.genBPdgId = ExtVar(cms.InputTag(_gf, 'genBPdgId'), int, doc='matched gen b-hadron pdgId (0 = none)')
        process.bplusTable.externalVariables.genPartIdx = ExtVar(cms.InputTag(_gf, 'genBIdx'), int, doc='row in Gen of the matched b-hadron (-1 = none)')

        # Per-leg generator match (leg 0/1 = muons, leg 2 = bachelor), plus the
        # common ancestor of the three matched legs. The ancestor's species is
        # the truth category: a B+ is signal, another b-hadron is a
        # mis-reconstructed b, none found is genuine combinatorial.
        #
        # Matching follows Bmm5: dR < 0.02 AND |dpt|/pt_gen < 0.1, both
        # required. An unmatched leg keeps -1 and is never assigned a nearest
        # neighbour, because a wrong match mislabels the very categories this
        # exists to separate.
        for _l in (0, 1, 2):
            setattr(process.bplusTable.externalVariables, f'leg{_l}GenIdx',
                    ExtVar(cms.InputTag(_gf, f'leg{_l}GenIdx'), int,
                           doc=f'row in Gen matched to leg {_l} (-1 = none)'))
            setattr(process.bplusTable.externalVariables, f'leg{_l}GenPdgId',
                    ExtVar(cms.InputTag(_gf, f'leg{_l}GenPdgId'), int,
                           doc=f'pdgId of the Gen particle matched to leg {_l} (0 = none)'))
            setattr(process.bplusTable.externalVariables, f'leg{_l}GenMotherPdgId',
                    ExtVar(cms.InputTag(_gf, f'leg{_l}GenMotherPdgId'), int,
                           doc=f'pdgId of that particle\'s mother (0 = none)'))
            setattr(process.bplusTable.externalVariables, f'leg{_l}GenPt',
                    ExtVar(cms.InputTag(_gf, f'leg{_l}GenPt'), float,
                           doc=f'generated pt of the particle matched to leg {_l}'))
            setattr(process.bplusTable.externalVariables, f'leg{_l}GenDR',
                    ExtVar(cms.InputTag(_gf, f'leg{_l}GenDR'), float,
                           doc=f'dR(leg {_l}, its Gen match)'))
        process.bplusTable.externalVariables.genAncestorPdgId = ExtVar(
            cms.InputTag(_gf, 'genAncestorPdgId'), int,
            doc='pdgId of the common ancestor of all matched legs (0 = none found)')
        process.bplusTable.externalVariables.genAncestorIdx = ExtVar(
            cms.InputTag(_gf, 'genAncestorIdx'), int,
            doc='row in Gen of that common ancestor (-1 = none)')
        process.bplusTable.externalVariables.genAncestorPt = ExtVar(
            cms.InputTag(_gf, 'genAncestorPt'), float, doc='common ancestor pt')
        process.bplusTable.externalVariables.genAncestorMass = ExtVar(
            cms.InputTag(_gf, 'genAncestorMass'), float, doc='common ancestor mass')
        process.bplusTable.externalVariables.nLegsGenMatched = ExtVar(
            cms.InputTag(_gf, 'nLegsGenMatched'), int,
            doc='number of legs with a Gen match (ancestor requires all three)')

    # Track -> Muon / Track -> PV cross-links, inverting the persisted
    # associations into row indices (-1 = none) on the Track table.
    process.trackMuonIdx = cms.EDProducer(
        'TrackToMuonIndexProducer',
        trackSrc=cms.InputTag(opts.srcTracks),
        association=cms.InputTag(opts.srcTracks + 'TrackToMuon'),
        # Track->Muon association is keyed to the alignment tracks directly, so
        # no originalIndex bridge is needed (leave the bridge inputs empty).
        originalIndex=cms.InputTag(''),
        pvSrc=cms.InputTag(''))
    process.trackPvIdx = cms.EDProducer(
        'TrackToVertexIndexProducer',
        trackSrc=cms.InputTag(opts.srcTracks),
        association=cms.InputTag('offlinePrimaryVertices'),
        # Track->PV association is keyed to generalTracks, so resolve via the
        # persisted originalIndex map (alignment track -> generalTracks row)
        # and the primary vertices' own track refs.
        originalIndex=cms.InputTag(opts.srcTracks, 'originalIndex'),
        pvSrc=cms.InputTag('offlinePrimaryVertices'))

    # Impact parameters w.r.t. the associated primary vertex. reco::Track's own
    # dxy/dz are measured from the origin, so the stored columns are unusable as
    # impact parameters; this producer supplies the corrected ones.
    process.trackImpactParameter = cms.EDProducer(
        'TrackImpactParameterProducer',
        trackSrc=cms.InputTag(opts.srcTracks),
        pvSrc=cms.InputTag('offlinePrimaryVertices'),
        pvIdx=cms.InputTag('trackPvIdx'),
        sentinel=cms.double(-99.))

    process.trackTable = cms.EDProducer(
        'SimpleTrackFlatTableProducer',
        src=cms.InputTag(opts.srcTracks),
        cut=cms.string(''), name=cms.string('Track'),
        doc=cms.string('AlCaReco cloned alignment tracks'),
        singleton=cms.bool(False), extension=cms.bool(False),
        variables=cms.PSet(
            P3Vars,
            charge=Var('charge', 'int16', doc='charge'),
            # dxy/dz here are reco::Track's own accessors, measured from the
            # ORIGIN (0,0,0), not from the primary vertex. They are kept so the
            # defect stays visible; d0/dzPV below are the ones to use.
            dxy=Var('dxy', float, doc='dxy w.r.t. the ORIGIN, not the PV'),
            dz=Var('dz', float, doc='dz w.r.t. the ORIGIN, not the PV'),
            normChi2=Var('normalizedChi2', float, doc='chi2/ndof'),
            nValidHits=Var('numberOfValidHits', 'int16', doc='n valid hits'),
            # Quality block. Verified present on the cloned alignment tracks:
            # the covariance survived cloning (ptError/dxyError/dzError finite
            # and positive for 100% of 752 tracks over 40 events) and so did the
            # hit pattern (pixel hits > 0 for 99.3%, tracker layers for 100%).
            ptErr=Var('ptError', float, doc='pt uncertainty'),
            dxyErr=Var('dxyError', float, doc='dxy uncertainty'),
            dzErr=Var('dzError', float, doc='dz uncertainty'),
            nValidPixelHits=Var('hitPattern().numberOfValidPixelHits()', 'int16',
                                doc='n valid pixel hits'),
            trackerLayers=Var('hitPattern().trackerLayersWithMeasurement()', 'int16',
                              doc='tracker layers with measurement'),
            pixelLayers=Var('hitPattern().pixelLayersWithMeasurement()', 'int16',
                            doc='pixel layers with measurement'),
            highPurity=Var('quality("highPurity")', bool, doc='highPurity quality flag'),
        ),
        externalVariables=cms.PSet(
            d0=ExtVar(cms.InputTag('trackImpactParameter', 'd0'), float,
                      doc='transverse impact parameter w.r.t. the associated PV'),
            dzPV=ExtVar(cms.InputTag('trackImpactParameter', 'dzPV'), float,
                        doc='longitudinal impact parameter w.r.t. the associated PV'),
            d0Err=ExtVar(cms.InputTag('trackImpactParameter', 'd0Err'), float,
                         doc='uncertainty on d0, track and PV contributions'),
            dzPVErr=ExtVar(cms.InputTag('trackImpactParameter', 'dzPVErr'), float,
                           doc='uncertainty on dzPV, track and PV contributions'),
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
    # Refit-arm leg sources (leaf-keyed). Empty -> the refit arm equals the raw
    # arm. Each leg: the single-track maker's refit collection + refitOk, plus
    # the input-track producer's candIdx/leafIdx that key refit track j back to
    # (candidate, leaf). Muon legs from the muon maker, bachelor from the kaon.
    _refit_legs = cms.VPSet()
    if opts.emitRefitTracks:
        _refit_legs = cms.VPSet(
            cms.PSet(
                tracks=cms.InputTag('globalCorJpsiKMuon', 'refit'),
                refitOk=cms.InputTag('globalCorJpsiKMuon', 'refitOk'),
                candIdx=cms.InputTag('bplusJpsiMuonTracks', 'candIdx'),
                leafIdx=cms.InputTag('bplusJpsiMuonTracks', 'leafIdx'),
            ),
            cms.PSet(
                tracks=cms.InputTag('globalCorJpsiKKaon', 'refit'),
                refitOk=cms.InputTag('globalCorJpsiKKaon', 'refitOk'),
                candIdx=cms.InputTag('bplusBachelorTracks', 'candIdx'),
                leafIdx=cms.InputTag('bplusBachelorTracks', 'leafIdx'),
            ),
        )
    process.bplusFit = cms.EDProducer(
        'JpsiXKinematicFitProducer',
        src=_src_cands,
        refitLegs=_refit_legs,
        # Per-leg gen matching thresholds. Bmm5 (GenBmmProducer.cc:74) requires
        # dR < 0.02 AND |dpt|/pt < 0.1 together, and that is what this produced
        # until 2026-08-18.
        #
        # The pT window is now DISABLED, deliberately, because it biases the one
        # quantity this analysis measures. A badly measured kaon fails the window,
        # so its candidate loses its common ancestor and leaves the signal
        # category: measured, 3.6% of B+-like candidates lose their bachelor leg
        # that way, and those are the broad ones -- p16-p84 of 5.126-5.411 against
        # 5.250-5.298 for the matched. The signal template's mass tails were being
        # cut by the truth definition rather than by the detector, and the kaon
        # momentum response was hard-truncated at |resp - 1| = 0.0999.
        #
        # dR is tightened 0.02 -> 0.01 to pay for it. Track direction is measured
        # orders of magnitude better than momentum, so a dR-only match biases a
        # variable the calibration never uses. The Bmm5 warning still stands in
        # principle -- dR alone can match the wrong track in a dense b jet -- so
        # the wrong-match rate is measured rather than assumed: the bachelor's
        # matched gen particle must be a daughter of the same ancestor as the
        # muons, and the fraction that is not is the contamination this buys.
        # dR back to 0.02, having measured that 0.01 was too aggressive: at 0.01
        # the three-leg match rate fell 78.5% -> 69.4%, which cost more than
        # dropping the pT window gained. And the wrong-match rate -- the
        # bachelor's matched particle not traceable to the ancestor -- was
        # **0.4% either way**, so tightening dR bought nothing it was meant to.
        genLegMaxDR=cms.untracked.double(0.02),
        genLegMaxRelDPt=cms.untracked.double(1e9),
        # Depth budget for the common-ancestor search. 10 and 30 were compared
        # on 1068 candidates and give byte-identical results, so the residual
        # ~13% of three-leg-matched signal that finds no ancestor is NOT the
        # budget running out. Kept at 10; raising it only costs permutations.
        genAncestorMaxDepth=cms.untracked.uint32(10),
        jpsiConstraint=cms.string(str(opts.jpsiConstraint)),
        jpsiMass=cms.double(3.0969),
        maxChi2=cms.double(-1.),
        beamSpot=cms.InputTag('offlineBeamSpot'),
        # PV collection for the true 3D dimuon flight-length significance
        # (dimuonSl3d). offlinePrimaryVertices is persisted in the AlCaReco.
        primaryVertices=cms.InputTag('offlinePrimaryVertices'),
        # Gen-matching (MC only): match the candidate to the closest last-copy
        # b-hadron. Empty on data -> genB* columns stay sentinel.
        genParticles=cms.InputTag('genParticles') if opts.isMC else cms.InputTag(''),
    )

    # ---- Joint N-body CVH arm --------------------------------------------
    # One fit for all N legs through a common vertex, emitting m, sigma_m,
    # jacMass and the vertex block directly -- no post-hoc vertex fit.
    # Cost (measured): ~0.87 s/event plus ~32 s startup, against a
    # ~3.4 s/event job with it on -- so it is roughly +50%, NOT the 10x that
    # was assumed before measuring. `jointCvh=False` drops the arm entirely.
    if opts.jointCvh:
        process.jointCvhBu = globalCorJpsiK.clone(
            src=cms.InputTag(opts.srcTracks),
            dedxSourceTracks=cms.InputTag(opts.srcTracks),
            srcCandidates=_src_cands,
            bCandIdxSrc=cms.InputTag(''),
            useIdealGeometry=cms.bool(False),
            corFiles=cms.vstring(),
            fillTrackTree=cms.bool(False),
            fillJac=cms.bool(False),
            produceValueMaps=cms.bool(True),
            # Reference-point jacobian WITHOUT the calibration payload: this
            # emits globalIdxs + jacRefLeg{i} + jacRefVtx and no Hessian.
            emitRefJacobian=cms.bool(bool(opts.emitRefJacobian)),
            # The clone source carries fillGradsFactored from _calib_pset,
            # whose driver default is True -- so this arm was running a
            # SelfAdjointEigenSolver over the full ~nParms x nParms Hessian for
            # EVERY candidate and emitting motherJacMass/jpsiJacMass/hessFactor
            # that nothing in this config reads. Off here. The calibration path
            # is mode="joint", which forces nanoOut off and builds its own
            # maker; it is unaffected.
            fillGradsFactored=cms.untracked.bool(False),
            fillGrads=cms.bool(False),
            outprefix=cms.untracked.string('jointcvh'),
            scalarPotentialInitFile=cms.string(opts.scalarPot3DInitFile),
            CvhMaster=CvhMasterPSet.clone(
                Particles=cms.vstring('mu+', 'mu-', 'kaon+', 'kaon-')),
        )
        process.jointCvhBu.__dict__['_TypedParameterizable__type'] = \
            'ResidualGlobalCorrectionMakerNTrackG4e'
        process.RandomNumberGeneratorService.jointCvhBu = cms.PSet(
            initialSeed=cms.untracked.uint32(123456789),
            engineName=cms.untracked.string('HepJamesRandom'))

    # ---- Vertex geometry, one instance per vertex block ---------------------
    # Every fitter emits the same block; this module is instantiated per
    # emitter rather than each fitter deriving its own geometry. It does NO
    # fitting; see CandidateVertexGeometryProducer.
    def _vtxGeom(src, prefix, okTag):
        return cms.EDProducer(
            'CandidateVertexGeometryProducer',
            srcCandidates=_src_cands,
            vertexSrc=cms.InputTag(src),
            blockPrefix=cms.string(prefix),
            fitOk=okTag,
            beamSpot=cms.InputTag('offlineBeamSpot'),
            primaryVertices=cms.InputTag('offlinePrimaryVertices'),
        )
    process.vtxGeomKvfRaw = _vtxGeom('bplusFit', 'raw', cms.InputTag('bplusFit', 'rawFitOk'))
    process.vtxGeomKvfCvh = _vtxGeom('bplusFit', 'ref', cms.InputTag('bplusFit', 'refFitOk'))
    if opts.jointCvh:
        process.vtxGeomJointCvh = _vtxGeom('jointCvhBu', '', cms.InputTag('jointCvhBu', 'fitOk'))

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
        # Refit muon legs (the J/psi leaves), aligned with bplusJpsiMuonTracks.
        process.refitMuTrackTable = cms.EDProducer(
            'SimpleTrackFlatTableProducer',
            src=cms.InputTag('globalCorJpsiKMuon', 'refit'),
            cut=cms.string(''), name=cms.string('RefitMuTrack'),
            doc=cms.string('CVH-refit J/psi muon tracks (aligned with input tracks)'),
            singleton=cms.bool(False), extension=cms.bool(False),
            variables=cms.PSet(
                P3Vars,
                charge=Var('charge', 'int16', doc='charge'),
                dxy=Var('dxy', float, doc='dxy'), dz=Var('dz', float, doc='dz'),
                normChi2=Var('normalizedChi2', float, doc='chi2/ndof'),
                ptErr=Var('ptError', float, doc='pt uncertainty from the refit 5x5'),
            ),
            externalVariables=cms.PSet(
                refitOk=ExtVar(cms.InputTag('globalCorJpsiKMuon', 'refitOk'),
                               int, doc='1 = refit succeeded, 0 = input copy'),
            ),
        )
        _extra_tables.append(process.refitMuTrackTable)

    if opts.isMC:
        # Full genParticles as a browsable Gen table; BuJpsiK_genPartIdx indexes
        # into it. Plus the per-event generator weight (GenEventInfoProduct).
        # DO NOT ADD A `cut` HERE.
        #
        # `cut=''` is what makes table row i equal genParticles[i]. Both
        # `genPartIdxMother` below and `BuJpsiK_genPartIdx` (the candidate's
        # matched b-hadron row) are raw indices into this collection. A cut
        # would renumber the rows while leaving both columns pointing at the
        # old numbering -- wrong answers, no error, nothing to notice.
        #
        # If the table ever must be slimmed, use a GenParticlePruner (as
        # standard NanoAOD does via `finalGenParticles`), which fixes up the
        # mother references to point within its own output. A `cut` on this
        # producer does not.
        process.genTable = cms.EDProducer(
            'SimpleGenParticleFlatTableProducer',
            src=cms.InputTag('genParticles'),
            cut=cms.string(''), name=cms.string('Gen'),
            doc=cms.string('generator particles (full genParticles, unpruned)'),
            singleton=cms.bool(False), extension=cms.bool(False),
            variables=cms.PSet(
                P3Vars,
                mass=Var('mass', float, doc='mass'),
                pdgId=Var('pdgId', int, doc='PDG id'),
                status=Var('status', 'int16', doc='status'),
                charge=Var('charge', 'int16', doc='charge'),
                # Row index of the first mother, -1 at the record's roots. Valid
                # only because cut='' keeps rows 1:1 with genParticles.
                genPartIdxMother=Var('?numberOfMothers>0?motherRef(0).key():-1',
                                     'int16', doc='row in Gen of the first mother '
                                                  '(-1 = none)'),
            ),
        )
        process.genWeightTable = cms.EDProducer(
            'SimpleGenEventFlatTableProducer',
            src=cms.InputTag('generator'),
            name=cms.string('GenEvt'),
            doc=cms.string('generator event info'),
            extension=cms.bool(False),
            variables=cms.PSet(
                weight=Var('weight', float, doc='generator event weight'),
            ),
        )
        # Pileup. `PileupSummaryInfos_addPileupInfo` is persisted in the MC
        # AlCaReco (via the keepPileupSummaryInfo customise), and CMSSW's own
        # NPUTablesProducer consumes vector<PileupSummaryInfo> directly, so no
        # new plugin is needed. Data carries no such product and gets no table.
        process.puTable = cms.EDProducer(
            'NPUTablesProducer',
            src=cms.InputTag('addPileupInfo'),
            pvsrc=cms.InputTag('offlinePrimaryVertices'),
            zbins=cms.vdouble([0.0, 1.7, 2.6, 3.0, 3.5, 4.2, 5.2, 6.0, 7.5, 9.0, 12.0]),
            savePtHatMax=cms.bool(False),
        )
        _extra_tables += [process.genTable, process.genWeightTable, process.puTable]

    # ---- Reference-point jacobian tables ---------------------------------
    # A nano flat-table column must be scalar, and these jacobians are
    # 5 (or 3) rows x nParms with nParms varying per track. The stock NanoAOD
    # FlattenedValueMapVectorTableProducer is the mechanism for exactly this --
    # muons_cff.py already uses it to store cvhJacRef/cvhmergedGlobalIdxs on
    # the Muon table. It emits per-object counts plus one concatenated values
    # table, which is the standard jagged-nano idiom and far more compact than
    # one row per matrix element.
    #
    # Constraint: within a table, all vector maps of the SAME TYPE must have
    # equal length per object, or the producer throws. int and float lengths
    # are tracked separately, so globalIdxs (nParms, int) rides along with the
    # float jacobians. momCov (9 floats) gets its own table for that reason.
    #
    # Two keying schemes, because the two fits differ:
    #   BuJpsiKJac  joint N-body arm, keyed to the CANDIDATE. Every float map
    #               is 3*nParms: 3 momentum rows per leg, plus the shared
    #               vertex jacobian. There is no 5-parameter per-leg state in
    #               a fit whose vertex is imposed by the parameterisation.
    #   *TrackJac   single-track arms, keyed to the input TRACK collection.
    #               5*nParms -- the full (qop, lambda, phi, dxy, dsz)
    #               reference state, i.e. the AN's dxref.
    if opts.emitRefJacobian:
        _JPREC = 12   # matches muons_cff.py's cvhJacRef
        _IPREC = 16   # global indices need 16 bits today
        if opts.jointCvh:
            _jvars = cms.PSet(
                cvhGlobalIdxs=ExtVar(cms.InputTag('jointCvhBu', 'globalIdxs'),
                                     'std::vector<int>',
                                     doc='global correction-parameter indices', precision=_IPREC),
                cvhJacRefVtx=ExtVar(cms.InputTag('jointCvhBu', 'jacRefVtx'),
                                    'std::vector<float>',
                                    doc='d(common vertex xyz)/d(globalparms), 3 x nParms row-major',
                                    precision=_JPREC),
            )
            for _i in range(3):
                setattr(_jvars, 'cvhJacRefLeg%d' % _i,
                        ExtVar(cms.InputTag('jointCvhBu', 'jacRefLeg%d' % _i),
                               'std::vector<float>',
                               doc='d(leg %d momentum)/d(globalparms), 3 x nParms row-major' % _i,
                               precision=_JPREC))
            process.bplusJacTable = cms.EDProducer(
                'FlattenedCandValueMapVectorTableProducer',
                name=cms.string('BuJpsiKJac'),
                src=_src_cands,
                cut=cms.string(''),
                doc=cms.string('joint N-body CVH reference-point jacobian'),
                variables=_jvars,
            )
            _extra_tables.append(process.bplusJacTable)

        def _trackJacTable(src, maker, name, doc):
            return cms.EDProducer(
                'FlattenedTrackValueMapVectorTableProducer',
                name=cms.string(name), src=src, cut=cms.string(''),
                doc=cms.string(doc),
                variables=cms.PSet(
                    cvhGlobalIdxs=ExtVar(cms.InputTag(maker, 'trackGlobalIdxs'),
                                         'std::vector<int>',
                                         doc='global correction-parameter indices',
                                         precision=_IPREC),
                    cvhJacRef=ExtVar(cms.InputTag(maker, 'trackJacRef'),
                                     'std::vector<float>',
                                     doc='d(qop,lambda,phi,dxy,dsz)/d(globalparms), '
                                         '5 x nParms row-major',
                                     precision=_JPREC),
                ),
            )
        process.kaonJacTable = _trackJacTable(
            _bach_src, 'globalCorJpsiKKaon', 'KaonTrackJac',
            'bachelor single-track CVH reference-point jacobian')
        _extra_tables.append(process.kaonJacTable)
        if opts.emitRefitTracks:
            process.muJacTable = _trackJacTable(
                cms.InputTag('bplusJpsiMuonTracks'), 'globalCorJpsiKMuon',
                'MuTrackJac', 'J/psi muon single-track CVH reference-point jacobian')
            _extra_tables.append(process.muJacTable)

    process.nanoTables = cms.Task(
        process.bplusTable, process.trackTable, process.muonTable,
        process.pvTable, process.dcsTable, process.l1Table, process.bplusFit,
        process.bplusLeafIdx, process.trackMuonIdx, process.trackPvIdx,
        process.trackImpactParameter,
        process.vtxGeomKvfRaw, process.vtxGeomKvfCvh,
        *( [process.jointCvhBu, process.vtxGeomJointCvh] if opts.jointCvh else [] ),
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

# Momentum floor -> every residual-maker instance actually scheduled (see the
# `_clampFloor` block next to PropagationPtotLimit above).
_clamped = []
for _n, _m in process.producers.items():
    if str(_m.type_()).startswith('ResidualGlobalCorrectionMaker'):
        _m.clampMomentumFloor = cms.double(_clampFloor)
        # Relative step damping + chi2 backtracking on top of the absolute
        # momentum floor; see the option help.
        _m.maxMomentumStepFactor = cms.double(float(opts.maxMomentumStepFactor))
        _m.stepBacktracking = cms.bool(bool(opts.stepBacktracking))
        _m.maxChi2Backtrack = cms.uint32(int(opts.maxChi2Backtrack))
        _m.stepBacktrackFromIter = cms.uint32(int(opts.stepBacktrackFromIter))
        _m.armijoC = cms.double(float(opts.armijoC))
        _m.armijoSlack = cms.double(float(opts.armijoSlack))
        _clamped.append(_n)
print("[cvh] effective: PropagationPtotLimit=%g GeV, clampMomentumFloor=%g GeV on %s"
      % (float(opts.plimit), _clampFloor, ",".join(_clamped) or "NO maker"))
print("[cvh] effective: maxMomentumStepFactor=%g, stepBacktracking=%s "
      "(fromIter=%d, maxChi2Backtrack=%d, armijoC=%g, armijoSlack=%g)"
      % (float(opts.maxMomentumStepFactor), bool(opts.stepBacktracking),
         int(opts.stepBacktrackFromIter), int(opts.maxChi2Backtrack),
         float(opts.armijoC), float(opts.armijoSlack)))

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
