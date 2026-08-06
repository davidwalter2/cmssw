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
from TrackPropagation.Geant4e.cvhMaster_cfi import CvhMasterPSet
process.globalCorKs = process.globalCorKs.clone(
    useIdealGeometry = False,
    outprefix = 'globalcor_ks',
    # MT G4Error master — same wiring as the J/psi driver. Particles set
    # to pi+/pi- since the Ks daughters propagate as pions.
    CvhMaster = CvhMasterPSet.clone(Particles=cms.vstring('pi+', 'pi-')),
    # Scalar-pot init file: needed even with the nominal CMSSW field
    # because the basis evaluator inside globalCor registers parmtype-14
    # modes from it. Default = David's polyfit3d coeff dump.
    scalarPotentialInitFile = cms.string('/work/submit/david_w/ZMass/mfs/data/fitresults/polyfit3d_full_coeffs_lmax18_cmsswnorm.txt'),
)

# Per-stream CLHEP engine; the maker calls setG4RandomEngineForStream() at
# the top of every produce() to wire this into Geant4's thread-local RNG.
process.RandomNumberGeneratorService = cms.Service('RandomNumberGeneratorService',
    globalCorKs = cms.PSet(
        initialSeed = cms.untracked.uint32(123456789),
        engineName = cms.untracked.string('HepJamesRandom'),
    ),
)
process.Geant4ePropagator.ForCVH = cms.bool(True)

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
    # geopro removed: CvhMasterThread (the maker's GlobalCache) owns the G4
    # world / master magnetic field init in an MT-safe way; running geopro
    # AS WELL double-inits G4 and segfaults in G4Region::G4Region.
    process.offlineBeamSpot *
    process.globalCorKs
)
