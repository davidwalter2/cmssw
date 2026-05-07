import FWCore.ParameterSet.Config as cms

# CVH refit for J/psi -> mu+ mu-, driven directly by the persisted
# VertexCompositeCandidateCollection from stage-1 ALCAReco
# (ALCARECOTkAlJpsiMuMuResonances). One tree row per candidate; no re-pairing
# at stage-2 (selection is fully delegated to stage-1).
globalCorJpsi = cms.EDProducer(
    'ResidualGlobalCorrectionMakerTwoTrackG4e',
    src = cms.InputTag('ALCARECOTkAlJpsiMuMu'),
    srcCandidates = cms.InputTag('ALCARECOTkAlJpsiMuMuResonances'),
    dedxSourceTracks   = cms.InputTag('ALCARECOTkAlJpsiMuMu'),
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
    bsConstraint = cms.bool(False),
    applyHitQuality = cms.bool(True),
    doVtxConstraint = cms.bool(False),
    doMassConstraint = cms.bool(False),
    massConstraint = cms.double(3.0969),          # J/psi mass, GeV
    massConstraintWidth = cms.double(9.29e-5),    # J/psi natural width, GeV
    # Both daughters are muons. The TwoBodyDecayCandidateProducer at stage-1
    # uses the positive-charge-first convention so daughter[0] is mu+; that
    # ordering does not affect the symmetric mass refit but is preserved for
    # downstream consumers.
    daughterParticleName1 = cms.string("mu"),
    daughterParticleName2 = cms.string("mu"),
    corFiles = cms.vstring(),
    MagneticFieldLabel = cms.string(''),
    outprefix = cms.untracked.string('globalcor_jpsi'),
)
