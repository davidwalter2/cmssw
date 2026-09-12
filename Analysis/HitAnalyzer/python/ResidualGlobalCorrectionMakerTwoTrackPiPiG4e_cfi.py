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
    doVtxConstraint = cms.bool(True),
    # Minimum size of a pair, required BEFORE the fit (see
    # ResidualGlobalCorrectionMakerTwoTrackG4e.cc). ndof = nvalid +
    # nvalidpixel - 10 (+3 beamspot, +1 pointing, +1 vertex constraint): one
    # measurement coordinate per strip hit, two per pixel hit, against the ten
    # state parameters the common vertex costs. minNdof = 1 therefore requires
    # more than NINE measurement coordinates with the vertex constraint on and
    # more than TEN with it off -- at ndof == 0 the fit is exactly determined
    # (chi2 identically zero, chi2/ndof undefined) and the factored-Hessian
    # export indexes past the end of its eigenvalue vector and aborts the
    # process. minPairHits is the same requirement read on VALID HITS rather
    # than on measurement coordinates; -1 = auto = 10 (constraint on) / 11
    # (off). 0 disables either.
    minNdof = cms.int32(1),
    minPairHits = cms.int32(-1),
    # minimum valid hits on the WEAKER leg; 0 = off (see the .cc). Default 8:
    # a thin leg is background (82-93 % of what it removes is `dup`/`unmatched`
    # by gen truth, at 0.9983 signal efficiency) and it is what every
    # non-finite mass-resolution export has in common.
    minLegHits = cms.int32(8),
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
    # coefficient dump file produced by mfs/dump_coeffs_for_cmssw.py; the basis
    # structure (l_max, mode list, Schmidt convention) and the initial
    # coefficients come from the file.
    scalarPotentialInitFile = cms.string(''),
    outprefix = cms.untracked.string('globalcor_ks'),
)
