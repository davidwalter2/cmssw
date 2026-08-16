## Single-track CVH refit driver for V0 (K_S -> pi pi / Lambda -> p pi)
## ALCARECO daughter tracks (data, Run 2 2016). Unlike runCvhV0.py this uses
## ResidualGlobalCorrectionMakerG4e (NOT TwoTrack): every selected V0-daughter
## track is refit independently, with no vertex / mass constraint. Purpose:
## per-track studies on a decay-in-flight-rich hadron sample — primarily the
## kink-finder score test (doKinkFinder), where pions (cTau 7.8 m) provide
## real in-tracker decays at the ~0.5-1% level vs the muon null.
##
## The single-track maker carries a hard p >= 2 GeV momentum floor (GN clamp),
## so V0 daughters are preselected with a p > pMin string cut (default 3 GeV)
## via a generic TrackSelector; the soft part of the spectrum is deliberately
## dropped rather than fit against the clamp.
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os

from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
from Configuration.AlCa.GlobalTag import GlobalTag

opts = VarParsing.VarParsing('analysis')
opts.register('input', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'comma-separated absolute paths or root:// URLs of V0 ALCARECO files')
opts.register('nEvents', -1, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, 'number of events to process (-1 = all)')
opts.register('collection', 'ALCARECOTkAlKsToPiPi', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'input track collection (ALCARECOTkAlKsToPiPi or '
              'ALCARECOTkAlLambdaToProtonPi)')
opts.register('particle', 'pi', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'track particle hypothesis: pi (KS daughters), kaon, proton')
opts.register('pMin', 3.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'track momentum preselection in GeV (single-track maker momentum '
              'floor is 2 GeV; stay above it)')
opts.register('useIdealGeometry', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'ideal (uncorrected) tracker geometry; default False = aligned from GT. '
              'NOTE: ideal geometry is NOT a kink-finder null (misalignment = real kinks)')
opts.register('scalarPot3DInitFile', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'coefficient dump file produced by mfs/dump_coeffs_for_cmssw.py (required)')
opts.register('numberOfThreads', 1, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'framework numberOfThreads (numberOfStreams follows)')
opts.register('propagationDirection', 'anyDirection', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'Geant4ePropagator PropagationDirection')
_defaultGroupsFile = os.path.join(os.environ.get('CMSSW_BASE', ''),
                                  'src/Analysis/HitAnalyzer/data/materialGroups50.txt')
opts.register('materialGroupsFile', _defaultGroupsFile, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'global material model grouping-tier rules file')
opts.register('doKinkFinder', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'per-material-step decay-in-flight score test (default True — '
              'this driver exists for the kink study)')
opts.register('kinkInjectLayer', -1, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'closure test: inject a synthetic kink at this step (-1 = off)')
opts.register('kinkInjectDqop', 0.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float, 'injected q/p step [1/GeV]')
opts.register('kinkInjectDxdz', 0.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float, 'injected dx/dz kink')
opts.register('kinkInjectDydz', 0.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float, 'injected dy/dz kink')
opts.register('fillJac', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'export the per-track reference Jacobian jacrefv (needed to see a '
              'reference change downstream of the fit)')
opts.register('fillGrads', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'export the Millepede objects gradv/hesspackedv (implies fillJac)')
opts.parseArguments()
if not opts.scalarPot3DInitFile:
    raise SystemExit("scalarPot3DInitFile=<path> is required (coefficient dump file)")

_g4particles = {'pi': ('pi+', 'pi-'), 'kaon': ('kaon+', 'kaon-'),
                'proton': ('proton', 'anti_proton'), 'mu': ('mu+', 'mu-')}
assert opts.particle in _g4particles, "particle must be one of %s" % list(_g4particles)

process = cms.Process("BENCH", Run2_2016)

process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.EndOfProcess_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.GeometrySimDB_cff")

process.GlobalTag = GlobalTag(process.GlobalTag, "auto:run2_data", "")
process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(
        record=cms.string("GeometryFileRcd"),
        tag=cms.string("XMLFILE_Geometry_2016_81YV1_Extended2016_mc"),
        label=cms.untracked.string("Extended"),
    ),
)
process.XMLFromDBSource.label = cms.string("Extended")

process.load("TrackPropagation.Geant4e.geantRefit_cff")
from TrackPropagation.Geant4e.cvhMaster_cfi import CvhMasterPSet

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(opts.nEvents))

_paths = [p.strip() for p in opts.input.split(',') if p.strip()]
assert _paths, "must set input=<paths> on the cmsRun command line"
_urls = [p if p.startswith(("root://", "file:")) else "file:" + p for p in _paths]
process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(*_urls),
    secondaryFileNames=cms.untracked.vstring(),
)

process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(int(opts.numberOfThreads)),
    numberOfStreams=cms.untracked.uint32(int(opts.numberOfThreads)),
    numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.RandomNumberGeneratorService.globalCor = cms.PSet(
    initialSeed=cms.untracked.uint32(123456789),
    engineName=cms.untracked.string('HepJamesRandom'),
)

process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")

# Momentum preselection: the single-track maker clamps Gauss-Newton steps at
# p = 2 GeV, so tracks below pMin are dropped upstream rather than fit
# against the clamp.
process.selectedV0Tracks = cms.EDFilter(
    "TrackSelector",
    src=cms.InputTag(opts.collection),
    cut=cms.string("p > %f" % float(opts.pMin)),
)

process.globalCor = cms.EDProducer(
    "ResidualGlobalCorrectionMakerG4e",
    src=cms.InputTag("selectedV0Tracks"),
    fitFromGenParms=cms.bool(False),
    fitFromSimParms=cms.bool(False),
    fillTrackTree=cms.bool(True),
    fillGrads=cms.bool(bool(opts.fillGrads)),
    fillJac=cms.bool(bool(opts.fillJac) or bool(opts.fillGrads)),
    fillRunTree=cms.bool(True),
    doGen=cms.bool(False),
    doSim=cms.bool(False),
    requireGen=cms.bool(False),
    doMuons=cms.bool(False),
    doMuonAssoc=cms.bool(False),
    doTrigger=cms.bool(False),
    doRes=cms.bool(False),
    useIdealGeometry=cms.bool(bool(opts.useIdealGeometry)),
    bsConstraint=cms.bool(False),
    applyHitQuality=cms.bool(True),
    corFiles=cms.vstring(),
    triggers=cms.vstring(),
    trackParticleName=cms.string(opts.particle),
    MagneticFieldLabel=cms.string(""),
    scalarPotentialInitFile=cms.string(opts.scalarPot3DInitFile),
    materialGroupsFile=cms.string(opts.materialGroupsFile),
    globalMaterialModel=cms.bool(True),
    perStepFieldModes=cms.bool(True),
    skipHitlessSurfaces=cms.bool(True),
    doKinkFinder=cms.bool(bool(opts.doKinkFinder)),
    kinkInjectLayer=cms.int32(int(opts.kinkInjectLayer)),
    kinkInjectDqop=cms.double(float(opts.kinkInjectDqop)),
    kinkInjectDxdz=cms.double(float(opts.kinkInjectDxdz)),
    kinkInjectDydz=cms.double(float(opts.kinkInjectDydz)),
    outprefix=cms.untracked.string("globalcor_v0single"),
    CvhMaster=CvhMasterPSet.clone(Particles=cms.vstring(*_g4particles[opts.particle])),
)

from MagneticField.ParametrizedEngine.parametrizedMagneticField_ScalarPot3D_cfi \
    import ParametrizedMagneticFieldProducer as ScalarPot3DMagneticFieldProducer
process.ScalarPot3DMagneticFieldProducer = ScalarPot3DMagneticFieldProducer.clone()
process.ScalarPot3DMagneticFieldProducer.parameters.InitFile = opts.scalarPot3DInitFile
fieldlabel = "ScalarPot3DMf"
process.ScalarPot3DMagneticFieldProducer.label = fieldlabel

process.geopro.MagneticFieldLabel = fieldlabel
process.Geant4ePropagator.MagneticFieldLabel = fieldlabel
process.Geant4ePropagator.ForCVH = cms.bool(True)
process.Geant4ePropagator.PropagationDirection = cms.string(opts.propagationDirection)
process.globalCor.MagneticFieldLabel = cms.string(fieldlabel)
from TrackPropagation.Geant4e.cvhMasterESProducer_cfi import cvhMasterESProducer
process.cvhMasterESProducer = cvhMasterESProducer.clone()
process.cvhMasterESProducer.MagneticFieldLabel = cms.string(fieldlabel)

process.reconstruction_step = cms.Path(
    process.offlineBeamSpot * process.selectedV0Tracks * process.globalCor
)
process.schedule = cms.Schedule(process.reconstruction_step)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
