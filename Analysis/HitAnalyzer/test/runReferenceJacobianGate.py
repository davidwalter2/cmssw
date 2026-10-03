## Reference-jacobian gate: validate hybrid2curvJacobianD against
## twoTrackPca2curvJacobianD.
##
## The N-body common-vertex reference jacobian is claimed to be a SCATTER of
## the existing base-class helper hybrid2curvJacobianD (5x6, hybrid params
## (qop, lam, phi, x, y, z)) rather than new SymPy. That helper is dead code
## in the tree -- declared, defined, called nowhere -- so it has to be shown
## correct before the state-layout rewrite leans on it.
##
## This config runs ResidualGlobalCorrectionMakerNTrackG4e (currently an exact
## fork of the two-track maker) with validateRefJacobian=True. On the first
## candidate the maker builds both matrices on synthetic states with d0 forced
## to 0 and prints the max deviation. PASS means the scatter reproduces the
## SymPy jacobian to double precision on the 9 shared columns.
##
## Usage:
##   cmsRun runD13Gate.py input=<alcareco.root> [nEvents=20]

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
from Configuration.AlCa.GlobalTag import GlobalTag

opts = VarParsing.VarParsing('analysis')
opts.register('input', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'input ALCARECO file')
opts.register('nEvents', 20, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, 'events to process')
opts.register('globalTag', 'auto:run2_data',
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'CMSSW GlobalTag')
opts.register('srcTracks', 'ALCARECOTkAlJpsiX',
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'track collection')
opts.register('srcCandidates', 'ALCARECOTkAlJpsiXBPlusResonances',
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'candidate collection')
opts.register('validateRefJacobian', True,
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'Run the one-shot reference-jacobian and N-track-mass self-test gates. Off for production '
              'runs: they print per candidate and are diagnostics, not output.')
opts.register('motherConstraintWidth', -1.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'N-body maker only: sigma of the MOTHER mass constraint. Negative '
              'leaves the maker default (1e-3 GeV). The two-track maker uses its '
              'cfi massConstraintWidth = 9.29e-5 (the PDG J/psi natural width), '
              'so the two disagree by 10.8x out of the box -- which weights the '
              'mass row 116x differently and is the prime suspect for the N=2 '
              'payload non-closure.')
opts.register('doMassConstraint', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'mass constraint on BOTH makers. False removes the mass row and '
              'with it the second-order convolution bias correction -- the '
              'isolating test for the payload-closure difference, since it is '
              'the only term the two makers compute by different code paths.')
opts.register('gradsPass', 'allcons', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'N-body maker only (the two-track maker ignores it); must be '
              '"free" when doMassConstraint=False')
opts.register('nIters', 10, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'Gauss-Newton iteration cap per constraint pass, applied to BOTH '
              'makers. Relevant to the payload closure: grad/Hess are evaluated '
              'wherever the fit halted, so a different stopping point shows up '
              'as a coherent payload offset.')
opts.register('edmConvergence', 1e-5, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'EDM convergence threshold, applied to BOTH makers (see nIters)')
opts.register('fillGrads', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'store the packed per-candidate global Hessian on BOTH makers -- '
              'the N=2 CALIBRATION-PAYLOAD closure. Earlier work validated '
              'the FITTED quantities to 6.9 keV; nothing has '
              'ever compared grad/Hess, and a silent grad bug would poison all '
              '92 fitted globals with no visible symptom.')
opts.register('fillGradsFactored', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'store the factored Hessian too (H = B^T B)')
opts.register('fillRunTree', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'write the runtree parameter catalog')
opts.register('globalMaterialModel', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'grouped material scales (parmtype 15) instead of per-module '
              'energy loss (parmtype 7) -- must match the reference channels')
opts.register('materialGroupsFile', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'materialGroups tier file, required with globalMaterialModel')
opts.register('refVtxConstraint', True,
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'doVtxConstraint on the TWO-TRACK reference. True is the like-for-like '
              'setting: with d0 free the reference carries a parameter the '
              'common-vertex layout structurally lacks.')
opts.register('scalarPot3DInitFile',
              '/work/submit/david_w/ZMass/mfs/data/fitresults/'
              'polyfit3d_full_coeffs_lmax18_cmsswnorm.txt',
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'scalar-potential coefficient dump')
opts.parseArguments()

process = cms.Process('D13GATE', Run2_2016)

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
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.offlineBeamSpot = cms.EDProducer('BeamSpotProducer')

# Calibration payload, applied IDENTICALLY to both makers -- the whole point of
# the closure is that only the module type differs. fillGradsFactored is read
# with getUntrackedParameter, so a tracked cms.bool would silently fall back to
# False and we would compare a factored Hessian against nothing.
_grads_pset = dict(
    fillGrads=cms.bool(bool(opts.fillGrads)),
    fillGradsFactored=cms.untracked.bool(bool(opts.fillGradsFactored)),
    fillRunTree=cms.bool(bool(opts.fillRunTree)),
    globalMaterialModel=cms.bool(bool(opts.globalMaterialModel)),
    materialGroupsFile=cms.string(str(opts.materialGroupsFile)),
)
if opts.globalMaterialModel:
    if not opts.materialGroupsFile:
        raise ValueError('globalMaterialModel=True requires materialGroupsFile')
    _grads_pset['skipHitlessSurfaces'] = cms.bool(True)

# The NTrack maker takes the two-track cfi verbatim -- it is currently an exact
# fork, which is exactly what makes the N=2 comparison meaningful.
from Analysis.HitAnalyzer.ResidualGlobalCorrectionMakerTwoTrackJpsiKMuMuG4e_cfi \
    import globalCorJpsiK

process.d13gate = globalCorJpsiK.clone(
    src=cms.InputTag(opts.srcTracks),
    dedxSourceTracks=cms.InputTag(opts.srcTracks),
    srcCandidates=cms.InputTag(opts.srcCandidates),
    bCandIdxSrc=cms.InputTag(''),
    useIdealGeometry=cms.bool(False),   # aligned geometry (see the memo)
    corFiles=cms.vstring(),
    fillTrackTree=cms.bool(True),
    fillJac=cms.bool(False),
    validateRefJacobian=cms.bool(bool(opts.validateRefJacobian)),
    outprefix=cms.untracked.string('d13_ntrack'),
    # Real 360-mode scalar-potential dump, matching the production driver.
    # A nonzero dB makes the gate STRONGER: the field correction enters both
    # jacobians, so a mishandled field term cannot cancel out of the diff.
    scalarPotentialInitFile=cms.string(opts.scalarPot3DInitFile),
    CvhMaster=CvhMasterPSet.clone(
        Particles=cms.vstring('mu+', 'mu-', 'kaon+', 'kaon-')),
    nIters=cms.uint32(int(opts.nIters)),
    edmConvergence=cms.double(float(opts.edmConvergence)),
    doMassConstraint=cms.bool(bool(opts.doMassConstraint)),
    gradsPass=cms.string(str(opts.gradsPass)),
    **_grads_pset,
)
if opts.motherConstraintWidth > 0.:
    process.d13gate.motherConstraintWidth = cms.double(float(opts.motherConstraintWidth))
process.d13gate.__dict__['_TypedParameterizable__type'] = \
    'ResidualGlobalCorrectionMakerNTrackG4e'

process.RandomNumberGeneratorService.d13gate = cms.PSet(
    initialSeed=cms.untracked.uint32(123456789),
    engineName=cms.untracked.string('HepJamesRandom'),
)

# Shared CVH G4 master (EventSetup product on CvhMasterRecord), consumed by the
# maker via esConsumes. Default field label (empty) -- the gate does not need
# the ScalarPot3D field *map*, only the scalar-potential dB correction, which
# the maker applies itself from scalarPotentialInitFile.
from TrackPropagation.Geant4e.cvhMasterESProducer_cfi import cvhMasterESProducer
process.cvhMasterESProducer = cvhMasterESProducer.clone()
process.cvhMasterESProducer.MagneticFieldLabel = cms.string('')
process.Geant4ePropagator.ForCVH = cms.bool(True)
process.Geant4ePropagator.PropagationPtotLimit = cms.double(0.05)

# ---- N=2 closure reference: the unmodified two-track maker, same events ----
# Same cfi, same config, only the module type differs -- so any difference in
# the fitted output is the state-layout rewrite and nothing else.
process.d13ref = process.d13gate.clone(
    validateRefJacobian=cms.bool(False),
    doVtxConstraint=cms.bool(bool(opts.refVtxConstraint)),
    outprefix=cms.untracked.string('d13_tworef'),
)
process.d13ref.__dict__['_TypedParameterizable__type'] = \
    'ResidualGlobalCorrectionMakerTwoTrackG4e'
process.RandomNumberGeneratorService.d13ref = cms.PSet(
    initialSeed=cms.untracked.uint32(123456789),
    engineName=cms.untracked.string('HepJamesRandom'),
)

process.p = cms.Path(process.offlineBeamSpot * process.d13gate * process.d13ref)
process.schedule = cms.Schedule(process.p)
