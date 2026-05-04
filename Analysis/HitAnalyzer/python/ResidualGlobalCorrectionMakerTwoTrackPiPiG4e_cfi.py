import FWCore.ParameterSet.Config as cms

# Configuration of the existing ResidualGlobalCorrectionMakerTwoTrackG4e
# C++ class for KS->pi+pi-: both daughter masses set to the charged pion.
# The class itself is already mass-parameterised, so this is a pure
# configuration clone (no new C++ needed).
globalCorKs = cms.EDProducer(
    'ResidualGlobalCorrectionMakerTwoTrackG4e',
    src = cms.InputTag('ALCARECOTkAlKsToPiPi'),
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
    massConstraintWidth = cms.double(1.e-5),
    daughterMass1 = cms.double(0.139570),              # both daughters: pion
    daughterMass2 = cms.double(0.139570),
    daughterMass1Err = cms.double(1.e-6),
    daughterMass2Err = cms.double(1.e-6),
    minPairMass = cms.double(0.40),               # KS V0Producer-like window
    maxPairMass = cms.double(0.60),
    respectTrackOrder = cms.bool(False),          # symmetric: order doesn't matter
    corFiles = cms.vstring(),
    MagneticFieldLabel = cms.string(''),
    outprefix = cms.untracked.string('globalcor_ks'),
)
