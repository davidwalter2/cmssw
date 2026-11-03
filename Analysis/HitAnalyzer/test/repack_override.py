import FWCore.ParameterSet.Config as cms
process = cms.Process("REPACK")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(3))
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(
    "file:/ceph/submit/data/user/d/david_w/mz/mc/inclusive_btojpsix_2016postvfp/mc_inclusive_btojpsix_c3228824_p0.root"))
process.out = cms.OutputModule("PoolOutputModule",
    fileName=cms.untracked.string("file:/work/submit/david_w/ZMass/BtoJpsiX_MCprod/mc_split1_override.root"),
    splitLevel=cms.untracked.int32(1),
    overrideInputFileSplitLevels=cms.untracked.bool(True),   # <-- the missing piece
    outputCommands=cms.untracked.vstring("keep *"),
)
process.e = cms.EndPath(process.out)
