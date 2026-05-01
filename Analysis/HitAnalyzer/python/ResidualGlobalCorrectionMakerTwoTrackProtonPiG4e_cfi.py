import FWCore.ParameterSet.Config as cms

# Configuration of the existing ResidualGlobalCorrectionMakerTwoTrackKPiG4e
# C++ class for Lambda0->p pi-: asymmetric, kaon position holds the proton.
# The class itself is already mass-parameterised, so this is a pure
# configuration clone (no new C++ needed). Track ordering matters --
# track[0] gets kaonMass, track[1] gets pionMass -- so the upstream
# V0CandidateProducer (with tryBothAssignments=True) is responsible for
# emitting (proton-candidate, pion-candidate) in that order.
globalCorLambda = cms.EDProducer(
    'ResidualGlobalCorrectionMakerTwoTrackKPiG4e',
    src = cms.InputTag('ALCARECOTkAlLambdaToProtonPi'),
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
    kaonMass = cms.double(0.938272),              # proton (in track[0] slot)
    pionMass = cms.double(0.139570),
    kaonMassErr = cms.double(1.e-6),
    pionMassErr = cms.double(1.e-6),
    minPairMass = cms.double(1.05),               # Lambda window
    maxPairMass = cms.double(1.18),
    respectTrackOrder = cms.bool(True),           # asymmetric: track[0]=proton, track[1]=pion
    corFiles = cms.vstring(),
    MagneticFieldLabel = cms.string(''),
    outprefix = cms.untracked.string('globalcor_lambda'),
)
