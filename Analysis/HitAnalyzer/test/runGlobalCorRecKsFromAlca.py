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

# Override the Geant4e propagator's momentum threshold (default 0.5 GeV) so
# low-pT V0 daughters (KS pions can have p < 0.5 GeV) are not rejected at the
# `plimit` exit. ~70% of CVH propagation failures on the V0 channels were
# from this cut; lowering it to 0.05 GeV recovers them while staying above
# the regime where Geant4 step-finding becomes unreliable.
process.Geant4ePropagator.PropagationPtotLimit = cms.double(0.05)

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

# Edit this to point at whichever TkAlKsToPiPi.root you want to ntuplise.
process.source = cms.Source(
    'PoolSource',
    fileNames=cms.untracked.vstring(
        'file:/work/submit/david_w/ZMass/test_output_multifile/TkAlKsToPiPi.root'
    )
)

process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(1),
    numberOfStreams=cms.untracked.uint32(0)
)

# Step-2 V0 candidate finder: re-pairs the deduplicated daughter tracks,
# applies V0Producer-like cuts, picks the best candidate per event.
process.load('Analysis.HitAnalyzer.KsToPiPiCandidateProducer_cfi')
process.selectedKsTracks = process.KsToPiPiCandidateProducer.clone(
    tracks = cms.InputTag('ALCARECOTkAlKsToPiPi'),
    # other parameters (mass hypotheses, tryBothAssignments, V0 mass window,
    # cosThetaXYMin, LxyOverSigmaMin) are taken from the cfi defaults
    # (V0Producer-style: KS PDG +/-70 MeV, cosThetaXY > 0.998, Lxy/sigma > 15).
)

# CVH 2-track refit + flat tree, KS configuration (both pions).
process.load('Analysis.HitAnalyzer.ResidualGlobalCorrectionMakerTwoTrackPiPiG4e_cfi')
process.globalCorKs = process.globalCorKs.clone(
    src = 'selectedKsTracks',
    useIdealGeometry = False,
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
