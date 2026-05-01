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

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

# Use the existing 8h-test ALCARECO output as input. Edit this to point at
# whichever TkAlKsToPiPi.root you want to ntuplise.
process.source = cms.Source(
    'PoolSource',
    fileNames=cms.untracked.vstring(
        'file:/work/submit/david_w/ZMass/test_output_8h/TkAlKsToPiPi_8h.root'
    )
)

process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(1),
    numberOfStreams=cms.untracked.uint32(0)
)

# Step-2 V0 candidate finder: re-pairs the deduplicated daughter tracks,
# applies V0Producer-like cuts, picks the best candidate per event.
process.load('Analysis.HitAnalyzer.V0CandidateProducer_cfi')
process.selectedKsTracks = process.V0CandidateProducer.clone(
    tracks = cms.InputTag('ALCARECOTkAlKsToPiPi'),
    daughterMass1 = 0.139570,        # both pions
    daughterMass2 = 0.139570,
    tryBothAssignments = False,      # symmetric
    expectedV0Mass = 0.497611,       # KS PDG mass
    minV0Mass = 0.40,
    maxV0Mass = 0.60,
    cosThetaXYMin = 0.998,
    LxyOverSigmaMin = 15.0,
)

# CVH 2-track refit + flat tree, KS configuration (both pions).
process.load('Analysis.HitAnalyzer.ResidualGlobalCorrectionMakerTwoTrackPiPiG4e_cfi')
process.globalCorKs = process.globalCorKs.clone(
    src = 'selectedKsTracks',
    useIdealGeometry = False,
    respectTrackOrder = True,         # exactly two tracks per event in pair order
    outprefix = 'globalcor_ks',
)

# BeamSpot is consumed by V0CandidateProducer. The ALCARECO output does not
# keep offlineBeamSpot, so produce a fresh one from the standard service.
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
    process.selectedKsTracks *
    process.globalCorKs
)
