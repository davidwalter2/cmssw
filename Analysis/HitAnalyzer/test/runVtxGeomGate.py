## Gate: the maker emits the joint fit's vertex + covariance
## as ValueMaps (no PV token), and CandidateVertexGeometryProducer turns those
## into PV/BS-relative geometry without doing any fitting.
##
## Usage:  cmsRun runVtxGeomGate.py input=<alcareco.root> [nEvents=50]

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
from Configuration.AlCa.GlobalTag import GlobalTag

opts = VarParsing.VarParsing('analysis')
opts.register('input', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'input ALCARECO file')
opts.register('nEvents', 50, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, 'events to process')
opts.register('globalTag', 'auto:run2_data',
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'CMSSW GlobalTag')
opts.register('srcTracks', 'ALCARECOTkAlJpsiX',
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'track collection')
opts.register('srcCandidates', 'ALCARECOTkAlJpsiXBPlusResonances',
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'candidate collection')
opts.register('scalarPot3DInitFile',
              '/work/submit/david_w/ZMass/mfs/data/fitresults/'
              'polyfit3d_full_coeffs_lmax18_cmsswnorm.txt',
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'scalar-potential dump')
opts.parseArguments()

process = cms.Process('VTXGEOM', Run2_2016)
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')

process.GlobalTag = GlobalTag(process.GlobalTag, opts.globalTag, '')
process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(record=cms.string('GeometryFileRcd'),
             tag=cms.string('XMLFILE_Geometry_2016_81YV1_Extended2016_mc'),
             label=cms.untracked.string('Extended')))
process.XMLFromDBSource.label = cms.string('Extended')

process.load('TrackPropagation.Geant4e.geantRefit_cff')
from TrackPropagation.Geant4e.cvhMaster_cfi import CvhMasterPSet

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(opts.nEvents))
_url = opts.input if opts.input.startswith(('root://', 'file:')) else 'file:' + opts.input
process.source = cms.Source('PoolSource', fileNames=cms.untracked.vstring(_url))
process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(1),
    numberOfStreams=cms.untracked.uint32(1))
process.MessageLogger.cerr.FwkReport.reportEvery = 25
process.offlineBeamSpot = cms.EDProducer('BeamSpotProducer')

from Analysis.HitAnalyzer.ResidualGlobalCorrectionMakerTwoTrackJpsiKMuMuG4e_cfi \
    import globalCorJpsiK

process.cvhJoint = globalCorJpsiK.clone(
    src=cms.InputTag(opts.srcTracks),
    dedxSourceTracks=cms.InputTag(opts.srcTracks),
    srcCandidates=cms.InputTag(opts.srcCandidates),
    bCandIdxSrc=cms.InputTag(''),
    useIdealGeometry=cms.bool(False),
    corFiles=cms.vstring(),
    fillTrackTree=cms.bool(True),
    fillJac=cms.bool(False),
    # fillGradsFactored is UNTRACKED and separate from fillJac; it is what
    # gates the global-fit payload, including motherJacMass/jpsiJacMass.
    fillGradsFactored=cms.untracked.bool(True),
    produceValueMaps=cms.bool(True),
    outprefix=cms.untracked.string('vtxgeom'),
    scalarPotentialInitFile=cms.string(opts.scalarPot3DInitFile),
    CvhMaster=CvhMasterPSet.clone(
        Particles=cms.vstring('mu+', 'mu-', 'kaon+', 'kaon-')),
)
process.cvhJoint.__dict__['_TypedParameterizable__type'] = \
    'ResidualGlobalCorrectionMakerNTrackG4e'
process.RandomNumberGeneratorService.cvhJoint = cms.PSet(
    initialSeed=cms.untracked.uint32(123456789),
    engineName=cms.untracked.string('HepJamesRandom'))

from TrackPropagation.Geant4e.cvhMasterESProducer_cfi import cvhMasterESProducer
process.cvhMasterESProducer = cvhMasterESProducer.clone()
process.cvhMasterESProducer.MagneticFieldLabel = cms.string('')
process.Geant4ePropagator.ForCVH = cms.bool(True)
process.Geant4ePropagator.PropagationPtotLimit = cms.double(0.05)

# The middle module. No fitting: it reads the maker's vertex + covariance
# and the PV/BS, and emits pointing / flight-significance variables.
process.cvhJointGeom = cms.EDProducer(
    'CandidateVertexGeometryProducer',
    srcCandidates=cms.InputTag(opts.srcCandidates),
    vertexSrc=cms.InputTag('cvhJoint'),
    fitOk=cms.InputTag('cvhJoint', 'fitOk'),
    # The N-body maker emits ONE block, so its instance names are unprefixed.
    # JpsiXKinematicFitProducer emits one per arm -> blockPrefix 'raw' / 'ref'.
    blockPrefix=cms.string(''),
    beamSpot=cms.InputTag('offlineBeamSpot'),
    primaryVertices=cms.InputTag('offlinePrimaryVertices'),
)

process.out = cms.OutputModule(
    'PoolOutputModule',
    fileName=cms.untracked.string('vtxgeom_products.root'),
    outputCommands=cms.untracked.vstring(
        'drop *',
        'keep *_cvhJoint_*_*',
        'keep *_cvhJointGeom_*_*',
        'keep *_%s_*_*' % opts.srcCandidates,
    ),
)
process.p = cms.Path(process.offlineBeamSpot * process.cvhJoint * process.cvhJointGeom)
process.e = cms.EndPath(process.out)
process.schedule = cms.Schedule(process.p, process.e)
