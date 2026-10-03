import FWCore.ParameterSet.Config as cms

# CVH refit for Z -> mu+ mu-, candidate-driven (stage-1 selection).
globalCorZ = cms.EDProducer(
    'ResidualGlobalCorrectionMakerTwoTrackG4e',
    src = cms.InputTag('ALCARECOTkAlZMuMu'),
    srcCandidates = cms.InputTag('ALCARECOTkAlZMuMuResonances'),
    dedxSourceTracks   = cms.InputTag('ALCARECOTkAlZMuMu'),
    dedxHarmonic2      = cms.InputTag(''),
    dedxPixelHarmonic2 = cms.InputTag(''),
    dedxAllHarmonic2   = cms.InputTag(''),
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
    # THE LUMINOUS REGION AS A GAUSSIAN NOISE BLOCK: three rows on the common
    # vertex with the beam-width covariance, ONE per pair (see the beam-line
    # block in ResidualGlobalCorrectionMakerBase.h). ON here because the Z is PROMPT: it decays at the primary vertex, so the
    # common vertex IS in the luminous region.
    bsConstraint = cms.bool(True),
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
    massConstraint = cms.double(91.1876),       # Z mass, GeV
    massConstraintWidth = cms.double(2.4952),   # Z natural width, GeV
    daughterParticleName1 = cms.string("mu"),
    daughterParticleName2 = cms.string("mu"),
    corFiles = cms.vstring(),
    MagneticFieldLabel = cms.string(''),
    outprefix = cms.untracked.string('globalcor_z'),
)
