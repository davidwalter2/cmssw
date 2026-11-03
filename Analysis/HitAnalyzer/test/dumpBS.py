import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
from Configuration.AlCa.GlobalTag import GlobalTag
process = cms.Process("BS", eras.Run2_2016)
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, "auto:run2_data", "")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(
    "file:/ceph/submit/data/user/d/david_w/ZMass/alcareco/260506_AllResonances_charmonium_full/TkAlJpsiMuMu.root"))
process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")
process.dump = cms.EDAnalyzer("BeamSpotAnalyzer") if False else cms.EDAnalyzer("EventContentAnalyzer")
# simple: print BS via a python-free analyzer isn't trivial; use a PSetless approach
process.p = cms.Path(process.offlineBeamSpot)
process.add_(cms.Service("Tracer", dumpPathsAndConsumes=cms.untracked.bool(False)))
