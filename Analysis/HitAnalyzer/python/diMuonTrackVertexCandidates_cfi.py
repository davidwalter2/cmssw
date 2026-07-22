import FWCore.ParameterSet.Config as cms

# Opposite-sign muon-track pairs (from TrackProducerFromPatMuons) in a mass
# window, as reco::VertexCompositeCandidate with RecoChargedCandidate daughters
# -- the input form the CVH two-track refit consumes. Resonance-agnostic: set
# the window for Z / J/psi / Upsilon; clone for multiple resonances.
diMuonTrackVertexCandidates = cms.EDProducer(
    "DiMuonTrackVertexCandidateProducer",
    src = cms.InputTag("tracksfrommuons"),
    massMin = cms.double(50.0),
    massMax = cms.double(150.0),
    oppositeSign = cms.bool(True),
)
