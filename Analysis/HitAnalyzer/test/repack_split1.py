import FWCore.ParameterSet.Config as cms
process = cms.Process("REPACK")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(30))
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(
    "file:/ceph/submit/data/user/d/david_w/ZMass/alcareco/260506_AllResonances_charmonium_full/TkAlJpsiMuMu.root"))
process.out = cms.OutputModule("PoolOutputModule",
    fileName=cms.untracked.string("file:/work/submit/david_w/ZMass/BtoJpsiX_MCprod/TkAlJpsiMuMu_split1.root"),
    splitLevel=cms.untracked.int32(1),
    outputCommands=cms.untracked.vstring("keep *"),
)
process.e = cms.EndPath(process.out)
