import FWCore.ParameterSet.Config as cms

# Lambda0 -> p pi- candidate finder for step-2 of the alignment workflow.
# Configuration of the generic V0CandidateProducer C++ class for the
# asymmetric proton + pion case. tryBothAssignments=True means both
# (m_p, m_pi) and (m_pi, m_p) hypotheses are tried per track pair and the
# assignment whose post-fit mass is closer to the nominal Lambda mass is
# kept. Output ordering is (proton-candidate, pion-candidate) so the
# downstream ntuplizer with respectTrackOrder=True assigns the correct
# masses (track[0]=proton-mass, track[1]=pion-mass).
LambdaToProtonPiCandidateProducer = cms.EDProducer(
    'V0CandidateProducer',
    # inputs
    tracks   = cms.InputTag('ALCARECOTkAlLambdaToProtonPi'),
    beamSpot = cms.InputTag('offlineBeamSpot'),
    # mass hypotheses (track[0] = baryon, track[1] = pion)
    daughterMass1 = cms.double(0.938272),     # proton (or anti-proton)
    daughterMass2 = cms.double(0.139570),     # pion
    daughterMassErr = cms.double(1.e-6),
    tryBothAssignments = cms.bool(True),      # asymmetric: disambiguate per pair
    expectedV0Mass = cms.double(1.115683),    # Lambda PDG mass
    # Post-fit V0 mass window: asymmetric around the Lambda PDG mass
    # (-46 MeV / +44 MeV). Tightened from the previous +/-65 MeV window
    # to suppress combinatorial background.
    minV0Mass = cms.double(1.07),
    maxV0Mass = cms.double(1.16),
    # vertex / pointing / flight cuts (V0Producer-like)
    pvalMin           = cms.double(0.0),
    cosThetaXYMin     = cms.double(0.998),
    LxyOverSigmaMin   = cms.double(15.0),
    # per-track minimum (redundant safety net)
    minTrackPt = cms.double(0.1),
    # opposite-sign requirement
    applyChargeFilter = cms.bool(True),
    charge            = cms.int32(0),
)
