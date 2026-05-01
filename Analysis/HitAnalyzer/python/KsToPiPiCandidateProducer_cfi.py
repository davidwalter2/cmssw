import FWCore.ParameterSet.Config as cms

# KS -> pi+pi- candidate finder for step-2 of the alignment workflow.
# Configuration of the generic V0CandidateProducer C++ class with both
# daughter masses set to the charged pion mass; symmetric, no
# proton/pion-assignment ambiguity.
KsToPiPiCandidateProducer = cms.EDProducer(
    'V0CandidateProducer',
    # inputs
    tracks   = cms.InputTag('ALCARECOTkAlKsToPiPi'),
    beamSpot = cms.InputTag('offlineBeamSpot'),
    # mass hypotheses
    daughterMass1 = cms.double(0.139570),     # pi+
    daughterMass2 = cms.double(0.139570),     # pi-
    daughterMassErr = cms.double(1.e-6),
    tryBothAssignments = cms.bool(False),     # symmetric: nothing to swap
    expectedV0Mass = cms.double(0.497611),    # KS PDG mass
    # post-fit V0 mass window
    minV0Mass = cms.double(0.40),
    maxV0Mass = cms.double(0.60),
    # vertex / pointing / flight cuts (V0Producer-like)
    pvalMin           = cms.double(0.0),      # vertex chi2 p-value > pvalMin
    cosThetaXYMin     = cms.double(0.998),    # 2D pointing angle
    LxyOverSigmaMin   = cms.double(15.0),     # transverse flight significance
    # per-track minimum (redundant safety net)
    minTrackPt = cms.double(0.1),
    # opposite-sign requirement
    applyChargeFilter = cms.bool(True),
    charge            = cms.int32(0),
)
