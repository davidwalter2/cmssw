import FWCore.ParameterSet.Config as cms

# Stage-2 CVH refit of the bachelor kaon track of B+ -> J/psi K+. Runs
# independently of the dimuon refit (no link to the J/psi vertex inside
# this fit -- option-a; the B+ mass/vertex constraint, if ever needed,
# enters as a Lagrange term at the downstream GBL stage).
#
# `src` is the splitter's `bachelor` track collection (one Track per
# surviving B+; positionally aligned with `bachelorBCandIdx`). The
# `Track::TrackExtraRef` is preserved through the splitter's value-copy
# of `daughter(1)->track()`, so hits resolve via the persisted
# `vector<reco::TrackExtra>` in the input event.
#
# `trackParticleName = "kaon"` engages the kaon mass + G4 particle name
# in propagation. Already verified: `getParticleProperties("kaon")` =
# 0.4937 GeV (`plugins/ParticleProperties.cc:14`); `g4ParticleName(
# "kaon", +-1) = "kaon+" / "kaon-"`, both pre-registered in the default
# `CvhMasterPSet.Particles` list (`TrackPropagation/Geant4e/python/
# cvhMaster_cfi.py:46-53`).
globalCorJpsiKKaon = cms.EDProducer(
    'ResidualGlobalCorrectionMakerG4e',
    src = cms.InputTag('jpsiKCandidateSplitter', 'bachelor'),
    bCandIdxSrc = cms.InputTag('jpsiKCandidateSplitter', 'bachelorBCandIdx'),
    fitFromGenParms = cms.bool(False),
    fitFromSimParms = cms.bool(False),
    fillTrackTree = cms.bool(True),
    fillGrads = cms.bool(False),
    fillJac = cms.bool(True),
    fillRunTree = cms.bool(False),
    doGen = cms.bool(False),
    doSim = cms.bool(False),
    requireGen = cms.bool(False),
    doMuons = cms.bool(False),
    doMuonAssoc = cms.bool(False),
    doTrigger = cms.bool(False),
    doRes = cms.bool(False),
    # False by default: True together with the default empty corFiles is the
    # broken Stage-1-only configuration (ideal geometry, no corrections)
    # the maker itself warns about. The driver overrides per leg.
    useIdealGeometry = cms.bool(False),
    bsConstraint = cms.bool(False),
    applyHitQuality = cms.bool(True),
    corFiles = cms.vstring(),
    triggers = cms.vstring(),
    trackParticleName = cms.string('kaon'),
    MagneticFieldLabel = cms.string(''),
    # Decay-in-flight kink finder (per-material-step score test); off by
    # default, enable for the K -> mu nu / K -> pi pi0 contamination veto.
    doKinkFinder = cms.bool(False),
    outprefix = cms.untracked.string('globalcor_jpsik_kaon'),
)
