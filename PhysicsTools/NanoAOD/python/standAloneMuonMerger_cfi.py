import FWCore.ParameterSet.Config as cms

mergedStandAloneMuons = cms.EDProducer("StandAloneMuonMerger",
                                       standAlone = cms.InputTag("standAloneMuons"),
                                       standAloneUpdatedAtVtx = cms.InputTag("standAloneMuons", "UpdatedAtVtx"),
                                       muons = cms.InputTag("linkedMuons"),
                                       )
                                       
