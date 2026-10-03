## A2 coexistence test: run a single-track maker (ResidualGlobalCorrectionMakerG4e)
## AND a two-track maker (ResidualGlobalCorrectionMakerTwoTrackG4e) in the SAME
## cmsRun process, both consuming the ONE shared CVH Geant4 master
## (cvhMasterESProducer -> CvhMasterRecord). Before A2 this aborted: each maker
## owned its own edm::GlobalCache<CvhMasterThread> -> a second G4MTRunManagerKernel
## -> Geant4 singleton abort. It must now build G4 once and run both makers.
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
from Configuration.AlCa.GlobalTag import GlobalTag

opts = VarParsing.VarParsing('analysis')
opts.register('input', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'J/psi ALCARECO file')
opts.register('nEvents', 20, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, 'events')
opts.register('numberOfThreads', 2, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, 'threads/streams')
opts.register('scalarPot3DInitFile', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'scalar-potential coeff dump')
opts.parseArguments()
assert opts.input, "set input=<jpsi ALCARECO>"
assert opts.scalarPot3DInitFile, "set scalarPot3DInitFile=<coeff dump>"

process = cms.Process("COEXIST", Run2_2016)
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

process.MessageLogger.cerr.FwkReport.reportEvery = 5
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(opts.nEvents))
_url = opts.input if opts.input.startswith(("root://", "file:")) else "file:" + opts.input
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(_url))
process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(int(opts.numberOfThreads)),
    numberOfStreams=cms.untracked.uint32(int(opts.numberOfThreads)),
    numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1),
)
process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")

# Labelled 3D scalar-potential field.
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

# The ONE shared CVH Geant4 master.
from TrackPropagation.Geant4e.cvhMasterESProducer_cfi import cvhMasterESProducer
process.cvhMasterESProducer = cvhMasterESProducer.clone()
process.cvhMasterESProducer.MagneticFieldLabel = cms.string(fieldlabel)

# Two-track maker (J/psi -> mu mu).
from Analysis.HitAnalyzer.ResidualGlobalCorrectionMakerTwoTrackJpsiMuMuG4e_cfi import globalCorJpsi
process.globalCorTwo = globalCorJpsi.clone(
    srcCandidates=cms.InputTag(''),          # fall back to the j>i pair loop
    scalarPotentialInitFile=cms.string(opts.scalarPot3DInitFile),
    MagneticFieldLabel=cms.string(fieldlabel),
    outprefix=cms.untracked.string('coexist_two'),
)

# Single-track maker on the same ALCARECO tracks (different plugin class ->
# the case GlobalCache could never share).
from Analysis.HitAnalyzer.ResidualGlobalCorrectionMakerMuonG4e_cfi import ResidualGlobalCorrectionMakerMuonG4e
_single = ResidualGlobalCorrectionMakerMuonG4e.clone(
    src=cms.InputTag('ALCARECOTkAlJpsiMuMu'),
    doMuonAssoc=cms.bool(False),
    scalarPotentialInitFile=cms.string(opts.scalarPot3DInitFile),
    MagneticFieldLabel=cms.string(fieldlabel),
)
# The MC 3-instance scenario: nominal + ideal + beamspot single-track makers,
# all sharing the one G4 master, PLUS the two-track maker (4 instances, 2
# plugin classes). None of this could coexist before A2.
process.globalCorOne   = _single.clone(useIdealGeometry=cms.bool(True),  bsConstraint=cms.bool(False), outprefix=cms.untracked.string('coexist_one'))
process.globalCorIdeal = _single.clone(useIdealGeometry=cms.bool(True),  bsConstraint=cms.bool(False), outprefix=cms.untracked.string('coexist_ideal'))
process.globalCorBs    = _single.clone(useIdealGeometry=cms.bool(True),  bsConstraint=cms.bool(True),  outprefix=cms.untracked.string('coexist_bs'))

for _m, _seed in (("globalCorTwo", 123456789), ("globalCorOne", 234567890),
                  ("globalCorIdeal", 334567890), ("globalCorBs", 434567890)):
    setattr(process.RandomNumberGeneratorService, _m, cms.PSet(
        initialSeed=cms.untracked.uint32(_seed), engineName=cms.untracked.string('HepJamesRandom')))

process.p = cms.Path(process.offlineBeamSpot * process.globalCorTwo
                     * process.globalCorOne * process.globalCorIdeal * process.globalCorBs)
process.schedule = cms.Schedule(process.p)
