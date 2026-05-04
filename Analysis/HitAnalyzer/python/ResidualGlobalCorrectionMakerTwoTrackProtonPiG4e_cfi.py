import FWCore.ParameterSet.Config as cms

# Configuration of the existing ResidualGlobalCorrectionMakerTwoTrackG4e
# C++ class for Lambda0->p pi-: asymmetric, kaon position holds the proton.
# The class itself is already mass-parameterised, so this is a pure
# configuration clone (no new C++ needed). Track ordering matters --
# track[0] gets daughterMass1, track[1] gets daughterMass2 -- so the upstream
# V0CandidateProducer (with tryBothAssignments=True) is responsible for
# emitting (proton-candidate, pion-candidate) in that order.
globalCorLambda = cms.EDProducer(
    'ResidualGlobalCorrectionMakerTwoTrackG4e',
    src = cms.InputTag('ALCARECOTkAlLambdaToProtonPi'),
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
    massConstraintWidth = cms.double(1.e-5),
    daughterMass1 = cms.double(0.938272),              # proton (in track[0] slot)
    daughterMass2 = cms.double(0.139570),
    daughterMass1Err = cms.double(1.e-6),
    daughterMass2Err = cms.double(1.e-6),
    minPairMass = cms.double(1.05),               # Lambda window
    maxPairMass = cms.double(1.18),
    respectTrackOrder = cms.bool(True),           # asymmetric: track[0]=proton, track[1]=pion
    corFiles = cms.vstring(),
    MagneticFieldLabel = cms.string(''),
    outprefix = cms.untracked.string('globalcor_lambda'),
)
