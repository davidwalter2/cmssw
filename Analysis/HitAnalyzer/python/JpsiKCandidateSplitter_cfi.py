import FWCore.ParameterSet.Config as cms

# Reads the B+ -> J/psi K+ candidates produced upstream by
# `JpsiXCandidateProducer`. Each B+ is a VertexCompositeCandidate whose
# daughter(0) is the J/psi VCC and daughter(1) is the bachelor kaon
# RecoChargedCandidate. The splitter emits:
#   - `dimuon`             : VertexCompositeCandidateCollection (J/psi sub-cands)
#   - `bachelor`           : TrackCollection (kaon tracks)
#   - `dimuonBCandIdx`     : std::vector<int>, parallel to `dimuon`,
#                            element = source B+ candidate index
#   - `bachelorBCandIdx`   : std::vector<int>, parallel to `bachelor`
#
# Channel-species disambiguation comes from the per-leg pdgId branches
# emitted by the splitter (`bachelorPdgId`, `muon0PdgId`, `muon1PdgId`);
# the downstream single-track maker cfi sets `trackParticleName="kaon"`.
jpsiKCandidateSplitter = cms.EDProducer(
    'JpsiKCandidateSplitter',
    src=cms.InputTag('ALCARECOTkAlJpsiXBPlusResonances'),
)
