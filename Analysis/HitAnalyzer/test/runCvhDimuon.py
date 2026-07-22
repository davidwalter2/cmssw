## Runtime test of the dimuon (two-track) CVH refit EDM ValueMap output, on
## J/psi ALCARECO (which carries tracker hits). Feeds ALCARECOTkAlJpsiMuMu
## tracks to the new diMuonTrackVertexCandidates producer (J/psi window) ->
## trackrefitdimuon (produceValueMaps=True) -> writes the ValueMaps to EDM.
## Inspect corMass: it should peak at the J/psi mass with filled per-muon momenta.
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
from Configuration.AlCa.GlobalTag import GlobalTag

opts = VarParsing.VarParsing('analysis')
opts.register('input', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'J/psi ALCARECO file')
opts.register('nEvents', 50, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, 'events')
opts.register('numberOfThreads', 2, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, 'threads/streams')
opts.register('scalarPot3DInitFile', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'scalar-potential coeff dump')
opts.register('massMin', 2.6, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float, 'min mu-mu mass')
opts.register('massMax', 3.5, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float, 'max mu-mu mass')
opts.register('storeJacobians', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool, 'store the global-fit payload (fillGradsFactored)')
opts.parseArguments()
assert opts.input, "set input=<jpsi ALCARECO>"
assert opts.scalarPot3DInitFile, "set scalarPot3DInitFile=<coeff dump>"

process = cms.Process("CVHDIMU", Run2_2016)
process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.GeometrySimDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag = GlobalTag(process.GlobalTag, "auto:run2_data", "")
process.GlobalTag.toGet = cms.VPSet(cms.PSet(
    record=cms.string("GeometryFileRcd"),
    tag=cms.string("XMLFILE_Geometry_2016_81YV1_Extended2016_mc"),
    label=cms.untracked.string("Extended"),
))
process.XMLFromDBSource.label = cms.string("Extended")
process.load("TrackPropagation.Geant4e.geantRefit_cff")

process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(opts.nEvents))
_url = opts.input if opts.input.startswith(("root://", "file:")) else "file:" + opts.input
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(_url))
process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(int(opts.numberOfThreads)),
    numberOfStreams=cms.untracked.uint32(int(opts.numberOfThreads)),
    numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1),
)
process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")

# Labelled ScalarPot3D field + shared G4 master.
from MagneticField.ParametrizedEngine.parametrizedMagneticField_ScalarPot3D_cfi \
    import ParametrizedMagneticFieldProducer as ScalarPot3DMagneticFieldProducer
process.ScalarPot3DMagneticFieldProducer = ScalarPot3DMagneticFieldProducer.clone()
process.ScalarPot3DMagneticFieldProducer.parameters.InitFile = opts.scalarPot3DInitFile
fieldlabel = "ScalarPot3DMf"
process.ScalarPot3DMagneticFieldProducer.label = fieldlabel
for _c in ("geopro", "Geant4ePropagator", "stripCPEESProducer",
           "StripCPEfromTrackAngleESProducer", "siPixelTemplateDBObjectESProducer",
           "templates"):
    if hasattr(process, _c):
        getattr(process, _c).MagneticFieldLabel = cms.string(fieldlabel)
process.Geant4ePropagator.ForCVH = cms.bool(True)
process.Geant4ePropagator.PropagationDirection = cms.string("anyDirection")
from TrackPropagation.Geant4e.cvhMasterESProducer_cfi import cvhMasterESProducer
process.cvhMasterESProducer = cvhMasterESProducer.clone()
process.cvhMasterESProducer.MagneticFieldLabel = cms.string(fieldlabel)

# Dimuon candidate producer (J/psi window) + generic dimuon refit.
from Analysis.HitAnalyzer.diMuonTrackVertexCandidates_cfi import diMuonTrackVertexCandidates
from Analysis.HitAnalyzer.ResidualGlobalCorrectionMakerDiMuonG4e_cfi import ResidualGlobalCorrectionMakerDiMuonG4e
process.diMuonTrackVertexCandidates = diMuonTrackVertexCandidates.clone(
    src="ALCARECOTkAlJpsiMuMu", massMin=opts.massMin, massMax=opts.massMax)
process.trackrefitdimuon = ResidualGlobalCorrectionMakerDiMuonG4e.clone(
    src=cms.InputTag("ALCARECOTkAlJpsiMuMu"),
    srcCandidates=cms.InputTag("diMuonTrackVertexCandidates"),
    MagneticFieldLabel=cms.string(fieldlabel),
    scalarPotentialInitFile=cms.string(opts.scalarPot3DInitFile),
    fillGradsFactored=cms.untracked.bool(bool(opts.storeJacobians)),
)
process.RandomNumberGeneratorService.trackrefitdimuon = cms.PSet(
    initialSeed=cms.untracked.uint32(423456789), engineName=cms.untracked.string('HepJamesRandom'))

process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName=cms.untracked.string("file:cvhdimuon.root"),
    outputCommands=cms.untracked.vstring(
        "drop *",
        "keep *_diMuonTrackVertexCandidates_*_*",
        "keep *_trackrefitdimuon_*_*",
    ),
)
process.p = cms.Path(process.offlineBeamSpot * process.diMuonTrackVertexCandidates * process.trackrefitdimuon)
process.e = cms.EndPath(process.out)
process.schedule = cms.Schedule(process.p, process.e)
