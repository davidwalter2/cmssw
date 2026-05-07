## CVH refit, candidate-driven fast path, J/psi -> mu+ mu-.
## Reads the persisted ALCARECOTkAlJpsiMuMuResonances VertexCompositeCandidate
## collection produced by the new stage-1 ALCAREco -- one tree row per
## candidate, no track-pair re-combinatorics or stage-2 vertex pre-fit.
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

# Match the V0 driver: lower the Geant4e momentum cutoff to recover
# low-pT muons (J/psi muons can be soft).
process.Geant4ePropagator.PropagationPtotLimit = cms.double(0.05)

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.source = cms.Source(
    'PoolSource',
    fileNames=cms.untracked.vstring(
        'file:/ceph/submit/data/user/d/david_w/ZMass/alcareco/260506_AllResonances_charmonium_100files/merged/TkAlJpsiMuMu.root'
    ),
)

process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(1),
    numberOfStreams=cms.untracked.uint32(0),
)

# CVH 2-track refit, candidate-driven. The cfi already sets
# srcCandidates = ALCARECOTkAlJpsiMuMuResonances.
process.load('Analysis.HitAnalyzer.ResidualGlobalCorrectionMakerTwoTrackJpsiMuMuG4e_cfi')
process.globalCorJpsi = process.globalCorJpsi.clone(
    useIdealGeometry = False,
    outprefix = 'globalcor_jpsi',
)

# CVH base consumes offlineBeamSpot; ALCAREco doesn't keep it, so produce one.
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
    process.globalCorJpsi
)
