import FWCore.ParameterSet.Config as cms

# Configuration of the existing ResidualGlobalCorrectionMakerTwoTrackG4e
# C++ class for KS->pi+pi-: both daughter masses set to the charged pion.
# The class itself is already mass-parameterised, so this is a pure
# configuration clone (no new C++ needed).
globalCorKs = cms.EDProducer(
    'ResidualGlobalCorrectionMakerTwoTrackG4e',
    src = cms.InputTag('ALCARECOTkAlKsToPiPi'),
    # Default fast path: read the persisted *Resonances candidates produced
    # by VertexCompositeCandidateRemapper at stage-1. One tree row per
    # candidate, no re-pairing. The ALCAReco track collection is still
    # consumed (via `src`) for muon-matching and dE/dx lookups; daughter
    # TrackRefs in the candidate already point at it.
    srcCandidates = cms.InputTag('ALCARECOTkAlKsToPiPiResonances'),
    dedxSourceTracks   = cms.InputTag('ALCARECOTkAlKsToPiPi'),
    dedxHarmonic2      = cms.InputTag('ALCARECOTkAlKsToPiPiDeDxHarmonic2'),
    dedxPixelHarmonic2 = cms.InputTag('ALCARECOTkAlKsToPiPiDeDxPixelHarmonic2'),
    dedxAllHarmonic2   = cms.InputTag('ALCARECOTkAlKsToPiPiDeDxAllHarmonic2'),
    fitFromGenParms = cms.bool(False),
    fitFromSimParms = cms.bool(False),
    fillTrackTree = cms.bool(True),
    fillGrads = cms.bool(False),
    fillJac = cms.bool(False),
    fillRunTree = cms.bool(False),
    doGen = cms.bool(False),
    genParticles = cms.InputTag('genParticles'),
    pileupInfo = cms.InputTag('addPileupInfo'),
    doSim = cms.bool(False),
    requireGen = cms.bool(False),
    doMuons = cms.bool(False),
    muons = cms.InputTag('muons'),
    doMuonAssoc = cms.bool(False),
    doTrigger = cms.bool(False),
    triggers = cms.vstring(),
    doL1Trigger = cms.bool(False),
    l1Results = cms.InputTag('gtDigis', '', 'RECO'),
    l1Triggers = cms.vstring(),
    doRes = cms.bool(False),
    useIdealGeometry = cms.bool(True),
    bsConstraint = cms.bool(False),
    applyHitQuality = cms.bool(True),
    doVtxConstraint = cms.bool(False),
    doMassConstraint = cms.bool(False),
    massConstraint = cms.double(0.497611),       # KS mass
    massConstraintWidth = cms.double(7.351e-15),  # natural width Gamma = hbar/tau (KS), GeV
    # Per-daughter Geant4 particle base. Mass + mass uncertainty are
    # looked up from the PDG table in
    # Analysis/HitAnalyzer/interface/ParticleProperties.h .
    daughterParticleName1 = cms.string("pi"),
    daughterParticleName2 = cms.string("pi"),
    corFiles = cms.vstring(),
    MagneticFieldLabel = cms.string(''),
    # Scalar-potential B-field correction (absolute-field model). Path to a
    # Phase-A.5 dump file produced by mfs/dump_coeffs_for_cmssw.py; the basis
    # structure (l_max, mode list, Schmidt convention) and the initial
    # coefficients come from the file.
    scalarPotentialInitFile = cms.string(''),
    outprefix = cms.untracked.string('globalcor_ks'),
)
