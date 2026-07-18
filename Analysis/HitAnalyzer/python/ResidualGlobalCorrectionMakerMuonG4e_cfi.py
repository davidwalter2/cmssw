import FWCore.ParameterSet.Config as cms

import os

# Single-muon-track CVH refit for the WMass custom NanoAOD.
#
# This is the `trackrefit` producer wired into muons_cff.py: it refits every
# muon inner track (delivered as a reco::Track + pat::Muon association by the
# `tracksfrommuons` TrackProducerFromPatMuons) and emits the per-muon
# ValueMaps the muon tables consume (corPt/corEta/corPhi/corCharge/edmval/
# nValidHits/nValidPixelHits and the vector maps globalIdxs/jacRef/momCov).
#
# `doMuonAssoc = True` keys the output ValueMaps to the pat::Muon collection
# via the association carried on `src`. `doMuons/doGen/doSim/doTrigger` are
# all False for production NanoAOD (no debug TTree, no gen/sim/trigger paths),
# so the guarded `muons`/`genParticles`/`pileupInfo`/`triggers` inputs are not
# consumed.
#
# The B-field model + coefficient file and the CvhMaster field label are set
# at customise time by nano_cff.setup3DFieldForRefit (which routes the
# "ScalarPot3DMf" labelled field into MagneticFieldLabel here and into
# CvhMaster.MagneticFieldLabel). scalarPotentialInitFile is likewise filled
# by the helper.
#
# Only the nominal, real-geometry refit is provided here (bsConstraint =
# False, useIdealGeometry = False). The ideal-geometry (trackrefitideal) and
# beamspot-constrained (trackrefitbs) MC-only variants are deferred: they
# require a second producer instance in the same job, which is blocked by the
# one-G4-master-per-job (G4MTRunManagerKernel singleton) constraint of the
# per-producer CvhMasterThread GlobalCache. See the Stage-A2 notes in the
# migration plan.
_materialGroups50 = os.path.join(
    os.environ.get("CMSSW_BASE", ""),
    "src/Analysis/HitAnalyzer/data/materialGroups50.txt")

ResidualGlobalCorrectionMakerMuonG4e = cms.EDProducer(
    'ResidualGlobalCorrectionMakerG4e',
    src = cms.InputTag('tracksfrommuons'),
    fitFromGenParms = cms.bool(False),
    fitFromSimParms = cms.bool(False),
    fillTrackTree = cms.bool(False),
    fillGrads = cms.bool(False),
    fillJac = cms.bool(True),
    fillRunTree = cms.bool(False),
    doGen = cms.bool(False),
    doSim = cms.bool(False),
    requireGen = cms.bool(False),
    doMuons = cms.bool(False),
    doMuonAssoc = cms.bool(True),
    doTrigger = cms.bool(False),
    doRes = cms.bool(False),
    useIdealGeometry = cms.bool(False),
    bsConstraint = cms.bool(False),
    applyHitQuality = cms.bool(True),
    keepPixelEdgeHits = cms.bool(False),
    pixelMinSizeX = cms.int32(2),
    corFiles = cms.vstring(),
    triggers = cms.vstring(),
    trackParticleName = cms.string('mu'),
    # Set by nano_cff.setup3DFieldForRefit at customise time.
    MagneticFieldLabel = cms.string(''),
    scalarPotentialInitFile = cms.string(''),
    # Convergence / debug knobs (production defaults).
    nIters = cms.uint32(10),
    edmConvergence = cms.double(1e-5),
    debugPerIterDump = cms.bool(False),
    runFDClosure = cms.bool(False),
    epsilonFDClosure = cms.double(1e-4),
    # Global material model (parmtype-15 groups) + per-step field attribution.
    globalMaterialModel = cms.bool(True),
    perStepFieldModes = cms.bool(True),
    skipHitlessSurfaces = cms.bool(True),
    materialGroupsFile = cms.string(_materialGroups50),
    outprefix = cms.untracked.string('globalcor_muon'),
    # The Geant4 master is no longer per-producer: it is the shared
    # EventSetup product from cvhMasterESProducer (CvhMasterRecord), consumed
    # via esConsumes by the maker. See nano_cff.setup3DFieldForRefit.
)
