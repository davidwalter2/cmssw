import FWCore.ParameterSet.Config as cms

# Configuration of the existing ResidualGlobalCorrectionMakerTwoTrackG4e
# C++ class for Lambda0->p pi-: asymmetric, slot 0 holds the proton.
# The class itself is already mass-parameterised, so this is a pure
# configuration clone (no new C++ needed). Track ordering matters --
# track[0] gets daughterMass1, track[1] gets daughterMass2 -- and is
# guaranteed by V0Producer's convention (preserved through stage-1's
# VertexCompositeCandidateRemapper): daughter(0) = baryon (p / pbar),
# daughter(1) = pion.
globalCorLambda = cms.EDProducer(
    'ResidualGlobalCorrectionMakerTwoTrackG4e',
    src = cms.InputTag('ALCARECOTkAlLambdaToProtonPi'),
    # Default fast path: persisted Lambda candidates (positive-charge-first
    # convention from V0Producer means daughter[0] = p / pbar, daughter[1] =
    # pi). See PiPi cfi for details.
    srcCandidates = cms.InputTag('ALCARECOTkAlLambdaToProtonPiResonances'),
    dedxSourceTracks   = cms.InputTag('ALCARECOTkAlLambdaToProtonPi'),
    dedxHarmonic2      = cms.InputTag('ALCARECOTkAlLambdaToProtonPiDeDxHarmonic2'),
    dedxPixelHarmonic2 = cms.InputTag('ALCARECOTkAlLambdaToProtonPiDeDxPixelHarmonic2'),
    dedxAllHarmonic2   = cms.InputTag('ALCARECOTkAlLambdaToProtonPiDeDxAllHarmonic2'),
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
    massConstraint = cms.double(1.115683),        # Lambda mass
    massConstraintWidth = cms.double(2.502e-15),  # natural width Gamma = hbar/tau (Lambda), GeV
    # Per-daughter Geant4 particle base (track[0] = baryon = proton).
    # Mass + mass uncertainty are looked up from the PDG table in
    # Analysis/HitAnalyzer/interface/ParticleProperties.h .
    daughterParticleName1 = cms.string("proton"),
    daughterParticleName2 = cms.string("pi"),
    corFiles = cms.vstring(),
    MagneticFieldLabel = cms.string(''),
    # Scalar-potential B-field correction (absolute-field model); path to a
    # coefficient dump file produced by mfs/dump_coeffs_for_cmssw.py.
    scalarPotentialInitFile = cms.string(''),
    outprefix = cms.untracked.string('globalcor_lambda'),
)
