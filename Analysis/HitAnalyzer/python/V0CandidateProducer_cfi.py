import FWCore.ParameterSet.Config as cms

# Defaults: KS -> pi+pi- selection. For Lambda, override
# daughterMass1 = 0.93827, expectedV0Mass = 1.115683, tryBothAssignments = True,
# minV0Mass / maxV0Mass to the Lambda window.
V0CandidateProducer = cms.EDProducer(
    'V0CandidateProducer',
    # inputs
    tracks   = cms.InputTag('ALCARECOTkAlKsToPiPi'),
    beamSpot = cms.InputTag('offlineBeamSpot'),
    # mass hypotheses (track[0] gets daughterMass1, track[1] gets daughterMass2 downstream)
    daughterMass1 = cms.double(0.139570),     # pion
    daughterMass2 = cms.double(0.139570),     # pion
    daughterMassErr = cms.double(1.e-6),
    tryBothAssignments = cms.bool(False),     # symmetric for KS; True for Lambda
    expectedV0Mass = cms.double(0.497611),    # used to disambiguate when tryBothAssignments=True
    # post-fit V0 mass window
    minV0Mass = cms.double(0.40),
    maxV0Mass = cms.double(0.60),
    # vertex / pointing / flight cuts (V0Producer-like; loosen as needed)
    pvalMin           = cms.double(0.0),      # vertex chi2 p-value > pvalMin (0 = no cut)
    cosThetaXYMin     = cms.double(0.998),    # 2D pointing angle
    LxyOverSigmaMin   = cms.double(15.0),     # transverse flight significance
    # basic per-track cut (redundant with ALCARECO; cheap safety net)
    minTrackPt = cms.double(0.1),
    # opposite-sign requirement
    applyChargeFilter = cms.bool(True),
    charge            = cms.int32(0),
)
