import FWCore.ParameterSet.Config as cms

TkAlGoodIdMuonSelector = cms.EDFilter("MuonSelector",
    src = cms.InputTag('muons'),
    cut = cms.string('isGlobalMuon &'
                     'isTrackerMuon &'
                     'numberOfMatches > 1 &'
                     'globalTrack.hitPattern.numberOfValidMuonHits > 0 &'
                     'abs(eta) < 2.5 &'
                     'globalTrack.normalizedChi2 < 20.'),
    filter = cms.bool(True)
)

# Loose selector for the J/psi + X AlCaReco stream. Drops the AND on
# tracker/global muon so tracker-only muons are kept, and drops the two
# globalTrack.* sub-cuts (they dereference a null globalTrack on
# tracker-only muons, which would break expression evaluation). By
# construction this is a superset of TkAlGoodIdMuonSelector.
TkAlLooseIdMuonSelector = cms.EDFilter("MuonSelector",
    src = cms.InputTag('muons'),
    cut = cms.string('(isGlobalMuon | isTrackerMuon) &'
                     'abs(eta) < 2.5 &'
                     'numberOfMatches > 1'),
    filter = cms.bool(True)
)

TkAlRelCombIsoMuonSelector = cms.EDFilter("MuonSelector",
    src = cms.InputTag(''),
    cut = cms.string('(isolationR03().sumPt + isolationR03().emEt + isolationR03().hadEt)/pt  < 0.15'),
    filter = cms.bool(True)
)
