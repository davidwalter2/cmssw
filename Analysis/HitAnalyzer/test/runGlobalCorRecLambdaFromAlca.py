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

process.source = cms.Source(
    'PoolSource',
    fileNames=cms.untracked.vstring(
        'file:/work/submit/david_w/ZMass/test_output_multifile/TkAlLambdaToProtonPi.root'
    )
)

process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(1),
    numberOfStreams=cms.untracked.uint32(0)
)

# Step-2 V0 candidate finder: re-pairs the deduplicated daughter tracks,
# applies V0Producer-like cuts, tries both (p,pi) and (pi,p) assignments per
# pair, picks the assignment closer to m_Lambda, then keeps the best candidate
# per event in (proton-track, pion-track) order.
process.load('Analysis.HitAnalyzer.LambdaToProtonPiCandidateProducer_cfi')
process.selectedLambdaTracks = process.LambdaToProtonPiCandidateProducer.clone(
    tracks = cms.InputTag('ALCARECOTkAlLambdaToProtonPi'),
    # other parameters (mass hypotheses, tryBothAssignments, V0 mass window,
    # cosThetaXYMin, LxyOverSigmaMin) are taken from the cfi defaults
    # (V0Producer-style: Lambda PDG +/-50 MeV, cosThetaXY > 0.998,
    # Lxy/sigma > 15).
)

# CVH 2-track refit + flat tree, Lambda configuration (proton + pion).
process.load('Analysis.HitAnalyzer.ResidualGlobalCorrectionMakerTwoTrackProtonPiG4e_cfi')
process.globalCorLambda = process.globalCorLambda.clone(
    src = 'selectedLambdaTracks',
    useIdealGeometry = False,
    respectTrackOrder = True,         # track[0]=proton, track[1]=pion
    outprefix = 'globalcor_lambda',
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
    process.selectedLambdaTracks *
    process.globalCorLambda
)
