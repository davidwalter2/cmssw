## Ditrack (two-track, candidate-driven) CVH refit on the produced B->J/psi+X MC.
## Reads J/psi->mu mu candidates from ALCARECOTkAlJpsiXJpsiOnlyResonances and
## refits the two daughter tracks together. Field/ES-master setup mirrors the
## validated single-track driver (ScalarPot3D + shared cvhMasterESProducer).
import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
from Configuration.AlCa.GlobalTag import GlobalTag

INPUT = "file:/ceph/submit/data/user/d/david_w/mz/mc/inclusive_btojpsix_2016postvfp/mc_inclusive_btojpsix_c3228824_p52.root"
COEFF = "/work/submit/david_w/ZMass/mfs/data/fitresults/polyfit3d_full_coeffs_lmax18_cmsswnorm.txt"

process = cms.Process("TEST", eras.Run2_2016)
process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.GeometrySimDB_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("TrackPropagation.Geant4e.geantRefit_cff")

process.Geant4ePropagator.PropagationPtotLimit = cms.double(0.05)  # keep soft J/psi muons
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(INPUT))
process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(1),
    numberOfStreams=cms.untracked.uint32(0),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.RandomNumberGeneratorService.globalCorJpsi = cms.PSet(
    initialSeed=cms.untracked.uint32(223456789),
    engineName=cms.untracked.string("HepJamesRandom"),
)

# Candidate-driven two-track J/psi maker, repointed at the TkAlJpsiX collections.
process.load("Analysis.HitAnalyzer.ResidualGlobalCorrectionMakerTwoTrackJpsiMuMuG4e_cfi")
process.globalCorJpsi = process.globalCorJpsi.clone(
    useIdealGeometry=False,
    outprefix="globalcor_ditrack_mc",
    src=cms.InputTag("ALCARECOTkAlJpsiX"),
    srcCandidates=cms.InputTag("ALCARECOTkAlJpsiXJpsiOnlyResonances"),
    dedxSourceTracks=cms.InputTag("ALCARECOTkAlJpsiX"),
    scalarPotentialInitFile=cms.string(COEFF),
)

process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")

process.GlobalTag = GlobalTag(process.GlobalTag, "auto:run2_mc", "")
process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(
        record=cms.string("GeometryFileRcd"),
        tag=cms.string("XMLFILE_Geometry_2016_81YV1_Extended2016_mc"),
        label=cms.untracked.string("Extended"),
    ),
)
process.XMLFromDBSource.label = cms.string("Extended")

# ScalarPot3D field + shared CVH G4 master (mirror the single-track driver).
from MagneticField.ParametrizedEngine.parametrizedMagneticField_ScalarPot3D_cfi \
    import ParametrizedMagneticFieldProducer as ScalarPot3DMagneticFieldProducer
process.ScalarPot3DMagneticFieldProducer = ScalarPot3DMagneticFieldProducer.clone()
process.ScalarPot3DMagneticFieldProducer.parameters.InitFile = COEFF
fieldlabel = "ScalarPot3DMf"
process.ScalarPot3DMagneticFieldProducer.label = fieldlabel
process.geopro.MagneticFieldLabel = fieldlabel
process.Geant4ePropagator.MagneticFieldLabel = fieldlabel
process.Geant4ePropagator.ForCVH = cms.bool(True)
process.Geant4ePropagator.PropagationDirection = cms.string("anyDirection")
process.globalCorJpsi.MagneticFieldLabel = cms.string(fieldlabel)
from TrackPropagation.Geant4e.cvhMasterESProducer_cfi import cvhMasterESProducer
process.cvhMasterESProducer = cvhMasterESProducer.clone()
process.cvhMasterESProducer.MagneticFieldLabel = cms.string(fieldlabel)

process.p = cms.Path(process.offlineBeamSpot * process.globalCorJpsi)
