## Same as runGlobalCorRecKsFromAlca.py but pointing at the Charmonium-based
## ALCAREco we generated in /ceph/.../charmonium_full/ for apples-to-apples
## comparison with the J/psi variants (also from Charmonium).
import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
from Configuration.AlCa.GlobalTag import GlobalTag

process = cms.Process('TEST', eras.Run2_2016)

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.load('TrackPropagation.Geant4e.geantRefit_cff')

process.Geant4ePropagator.PropagationPtotLimit = cms.double(0.05)
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.source = cms.Source(
    'PoolSource',
    fileNames=cms.untracked.vstring(
        'file:/ceph/submit/data/user/d/david_w/ZMass/alcareco/260506_AllResonances_charmonium_100files/merged/TkAlKsToPiPi.root'
    )
)

process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(1),
    numberOfStreams=cms.untracked.uint32(0)
)

process.load('Analysis.HitAnalyzer.ResidualGlobalCorrectionMakerTwoTrackPiPiG4e_cfi')
process.globalCorKs = process.globalCorKs.clone(
    useIdealGeometry = False,
    outprefix = 'globalcor_ks',
)

process.offlineBeamSpot = cms.EDProducer('BeamSpotProducer')

process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(
        record=cms.string('GeometryFileRcd'),
        tag=cms.string('XMLFILE_Geometry_2016_81YV1_Extended2016_mc'),
        label=cms.untracked.string('Extended'),
    ),
)
process.XMLFromDBSource.label = cms.string('Extended')

process.p = cms.Path(
    process.geopro *
    process.offlineBeamSpot *
    process.globalCorKs
)
