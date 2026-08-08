## Single-track CVH refit driver for the inclusive B -> J/psi + X MC
## ALCARECO (ALCARECOTkAlJpsiX, ongoing production on the ceph CMS group
## store). Fits gen-matched tracks of a chosen species through the
## single-track maker — no vertex / mass constraint.
##
## Primary use: kink-finder decay-in-flight rates by TRUE species in the
## SAME sample: gen-matched kaons must show ~7x the pion decay-tag excess
## at equal momentum (lambda_K/lambda_pi = (ctau_K/m_K)/(ctau_pi/m_pi) =
## 0.135), muons are the null. Species selection = gen matching
## (status 1, |pdgId| = genMatchPdgId, same charge, dR < 0.1, pT within
## 50%) with requireGen=True. NOTE the 50% pT window biases against
## hard early decays (reconstructed pT far from gen) — tag rates from
## this sample are a lower bound for those.
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os

from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
from Configuration.AlCa.GlobalTag import GlobalTag

opts = VarParsing.VarParsing('analysis')
opts.register('input', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'comma-separated paths of B->J/psi+X MC ALCARECO files')
opts.register('inputFileList', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'text file with one input path per line (the production files '
              'hold ~15 events each, so lists of hundreds of files are normal); '
              'concatenated with input= if both are given')
opts.register('nEvents', -1, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, 'number of events to process (-1 = all)')
opts.register('particle', 'kaon', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'species pass: kaon (genMatchPdgId 321), pi (211), mu (13)')
opts.register('pMin', 3.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'track momentum preselection in GeV (maker momentum floor is 2 GeV)')
opts.register('useIdealGeometry', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'ideal tracker geometry; default False = MC-truth alignment from GT')
opts.register('scalarPot3DInitFile', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'coefficient dump file (required)')
opts.register('numberOfThreads', 1, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, 'framework threads/streams')
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
              'per-material-step decay-in-flight score test')
opts.register('genMatchDR', 0.1, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'dR search window of the gen match')
opts.register('doSimDecayTruth', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'Geant4 decay/interaction truth for the gen-matched particle '
              '(simTrk*/simVtx*/simDau* branches). Requires the v3+ production, '
              'which keeps SimTracks+SimVertices; earlier campaigns do not have '
              'them and the job will fail on the missing product')
opts.register('genMatchPtWindow', 10.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'relative pT window of the gen match. Default 10 = effectively '
              'off (dR + charge only), REQUIRED for decay studies: with the '
              'legacy 0.5 the reco pT of a decayed kaon (daughter mu can '
              'carry 5% of p) fails the match and decays are silently '
              'removed from the requireGen sample')
opts.parseArguments()
if not opts.scalarPot3DInitFile:
    raise SystemExit("scalarPot3DInitFile=<path> is required (coefficient dump file)")

_species = {'kaon': (321, ('kaon+', 'kaon-')),
            'pi':   (211, ('pi+', 'pi-')),
            'mu':   (13,  ('mu+', 'mu-'))}
assert opts.particle in _species, "particle must be one of %s" % list(_species)
_pdgid, _g4parts = _species[opts.particle]

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

# Conditions of the 106X UL2016 MC production chain (as runCvhJpsiGenMC.py).
process.GlobalTag = GlobalTag(process.GlobalTag, "106X_mcRun2_asymptotic_v17", "")
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
if opts.inputFileList:
    with open(opts.inputFileList) as _f:
        _paths += [l.strip() for l in _f if l.strip() and not l.startswith('#')]
assert _paths, "must set input=<paths> and/or inputFileList=<file>"
_urls = [p if p.startswith(("root://", "file:")) else "file:" + p for p in _paths]
process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(*_urls),
    secondaryFileNames=cms.untracked.vstring(),
    # The production writes each condor job with Run 1 / Lumi 1 and event
    # numbers restarting from 1, so (run, lumi, event) COLLIDE across files
    # while the events are genuinely distinct (different seeds). The default
    # duplicate check silently drops all but the first file's worth of
    # events (~15 events kept out of thousands). NOTE the flip side: tree
    # rows cannot be uniquely keyed by run/lumi/event for this sample.
    duplicateCheckMode=cms.untracked.string("noDuplicateCheck"),
    # The production is written by thousands of condor jobs and a handful of
    # outputs are zero-length (job died during the copy out). One such file
    # aborts the whole cmsRun with a FileOpenError, so skip rather than die.
    skipBadFiles=cms.untracked.bool(True),
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

process.selectedJpsiXTracks = cms.EDFilter(
    "TrackSelector",
    src=cms.InputTag("ALCARECOTkAlJpsiX"),
    cut=cms.string("p > %f" % float(opts.pMin)),
)

process.globalCor = cms.EDProducer(
    "ResidualGlobalCorrectionMakerG4e",
    src=cms.InputTag("selectedJpsiXTracks"),
    fitFromGenParms=cms.bool(False),
    fitFromSimParms=cms.bool(False),
    fillTrackTree=cms.bool(True),
    fillGrads=cms.bool(False),
    fillJac=cms.bool(False),
    fillRunTree=cms.bool(True),
    doGen=cms.bool(True),
    genParticles=cms.InputTag("genParticles"),
    pileupInfo=cms.InputTag("addPileupInfo"),
    genMatchPdgId=cms.int32(_pdgid),
    genMatchPtWindow=cms.double(float(opts.genMatchPtWindow)),
    # Every stable charged species competes for the dR match and the winner
    # must be the species of this pass. Without this, opening the pT window
    # (mandatory for decay studies) lets a soft gen hadron steal the match to
    # an unrelated J/psi muon track: on the v3 MC that mislabelled ~2/3 of
    # the "decayed kaon" sample and diluted the ROC from AUC 0.78 to 0.59.
    genMatchPdgIds=cms.vint32(11, 13, 211, 321, 2212),
    genMatchDR=cms.double(float(opts.genMatchDR)),
    doSim=cms.bool(False),
    # doSim (PSimHit-based) stays off: the ALCARECO keeps SimTracks and
    # SimVertices but NOT the TrackerHits*LowTof collections.
    doSimDecayTruth=cms.bool(bool(opts.doSimDecayTruth)),
    requireGen=cms.bool(True),
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
    outprefix=cms.untracked.string("globalcor_jpsix_%s" % opts.particle),
    CvhMaster=CvhMasterPSet.clone(Particles=cms.vstring(*_g4parts)),
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
    process.offlineBeamSpot * process.selectedJpsiXTracks * process.globalCor
)
process.schedule = cms.Schedule(process.reconstruction_step)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
