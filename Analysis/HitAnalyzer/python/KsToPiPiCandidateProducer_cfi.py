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
    # Daughter Geant4 particle base (both pions). Mass + mass uncertainty
    # are looked up from the PDG table in
    # Analysis/HitAnalyzer/interface/ParticleProperties.h .
    daughterParticleName1 = cms.string("pi"),
    daughterParticleName2 = cms.string("pi"),
    tryBothAssignments = cms.bool(False),     # symmetric: nothing to swap
    expectedV0Mass = cms.double(0.497611),    # KS PDG mass
    # Post-fit V0 mass window: KS PDG +/-60 MeV. Tightened from the previous
    # +/-100 MeV window to suppress combinatorial background.
    minV0Mass = cms.double(0.44),
    maxV0Mass = cms.double(0.56),
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
