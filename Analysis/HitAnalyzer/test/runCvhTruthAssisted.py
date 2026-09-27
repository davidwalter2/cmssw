## Truth-assisted module-level corrections on realistic-alignment MC, and the
## single-track closure that tests them.
##
##   mode=derive   single-track CVH fit of gen-matched muons with the five
##                 reference parameters frozen to the generator values
##                 (fitFromGenParms): no weak modes, the hit residuals measure
##                 the reconstruction geometry against the true trajectory.
##                 Stores gradv / hesspackedv / globalidxv for the global solve
##                 (the single-track maker has no factored Hessian).
##   mode=closure  the FREE single-track fit of the same kind of tracks, storing
##                 the reference-point jacobian jacrefv, so the fitted
##                 corrections can be applied linearly:
##                     (q/p)_cor = (q/p)_ref + sum_j jacrefv[0, j] x[globalidxv[j]]
##
##   inputType=alcareco  J/psi MC TkAlJpsiMuMu ALCARECO tracks (trackSrc)
##   inputType=miniaod   DY MiniAODv2: slimmedMuons -> TrackProducerFromPatMuons
##                       (innerTrackOnly=False), the chain the W-mass nano uses
##
## Fixed physics configuration, the one the nano v15 production refits with:
##   * module-level model: one Bz (parmtype 6) + one material (7) parameter per
##     module, hitless (dead / inactive) surfaces kept in the fit;
##   * pixel edge / single-pixel hits re-admitted with their class-correction
##     parameters 16-21 (pixelClassHits, default True);
##   * the SIMULATION's field: the default MC field of the process = the 160812
##     volume map with the OAE parametrisation inside the tracker, which is what
##     the UL16 SIM propagated through (the 160812 3D grid differs from it by
##     -1.8 mT at the origin up to +5.7 mT at the tracker ends);
##   * realistic reconstruction alignment from the GT (useIdealGeometry=False;
##     TrackerAlignment_2016_ultralegacymc_v1), the SIM being ideal.
##
## example:
##   cmsRun runCvhTruthAssisted.py mode=derive inputType=alcareco \
##       inputFileList=files.txt nEvents=-1 numberOfThreads=4
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
from Configuration.AlCa.GlobalTag import GlobalTag

opts = VarParsing.VarParsing('analysis')
opts.register('input', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'comma-separated input files')
opts.register('inputFileList', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'text file with one input file per line')
opts.register('nEvents', 500, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, 'events to process (-1 = all)')
opts.register('skipEvents', 0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, 'events to skip')
opts.register('numberOfThreads', 1, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, 'threads (= streams; one output file per stream)')
opts.register('mode', 'derive', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'derive | closure')
opts.register('inputType', 'alcareco', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'alcareco | miniaod')
opts.register('trackSrc', 'ALCARECOTkAlJpsiMuMu', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'track collection (inputType=alcareco)')
opts.register('muonSrc', 'slimmedMuons', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'pat::Muon collection (inputType=miniaod)')
opts.register('pixelClassHits', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              're-admit pixel edge / single-pixel hits + class-correction parameters 16-21')
opts.register('useIdealGeometry', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'reconstruct with the ideal geometry (the SIM geometry) instead of the '
              'realistic MC alignment')
opts.register('globalTag', '106X_mcRun2_asymptotic_v17', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'UL16 MC conditions')
opts.register('propagationPtotLimit', 0.2, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float, 'Geant4e momentum floor (GeV)')
opts.register('outprefix', 'globalcor_truth', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'output file prefix (<prefix>_<stream>.root)')
opts.parseArguments()

assert opts.mode in ('derive', 'closure'), 'mode must be derive|closure'
assert opts.inputType in ('alcareco', 'miniaod'), 'inputType must be alcareco|miniaod'
_derive = opts.mode == 'derive'
_mini = opts.inputType == 'miniaod'

process = cms.Process("TRUTH", Run2_2016)
process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.GeometrySimDB_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, opts.globalTag, "")
process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(
        record=cms.string("GeometryFileRcd"),
        tag=cms.string("XMLFILE_Geometry_2016_81YV1_Extended2016_mc"),
        label=cms.untracked.string("Extended"),
    ),
)
process.XMLFromDBSource.label = cms.string("Extended")
process.load("TrackPropagation.Geant4e.geantRefit_cff")

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(opts.nEvents))
_paths = [p.strip() for p in opts.input.split(',') if p.strip()]
if opts.inputFileList:
    with open(opts.inputFileList) as _f:
        _paths += [l.strip() for l in _f if l.strip() and not l.startswith('#')]
assert _paths, "set input=<paths> and/or inputFileList=<file>"
process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(*[p if p.startswith(("root://", "file:")) else "file:" + p
                                      for p in _paths]),
    skipEvents=cms.untracked.uint32(int(opts.skipEvents)),
    skipBadFiles=cms.untracked.bool(True),
    duplicateCheckMode=cms.untracked.string('noDuplicateCheck'),
)
process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(int(opts.numberOfThreads)),
    numberOfStreams=cms.untracked.uint32(int(opts.numberOfThreads)),
    numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.RandomNumberGeneratorService.globalCor = cms.PSet(
    initialSeed=cms.untracked.uint32(123456789),
    engineName=cms.untracked.string('HepJamesRandom'),
)

_seq = []
if _mini:
    # MiniAOD carries offlineBeamSpot (the tag the maker consumes): no producer.
    process.tracksfrommuons = cms.EDProducer(
        "TrackProducerFromPatMuons",
        src=cms.InputTag(opts.muonSrc),
        innerTrackOnly=cms.bool(False),
        ptMin=cms.double(-1.),
    )
    _seq.append(process.tracksfrommuons)
    _tracks, _gen, _pu = "tracksfrommuons", "prunedGenParticles", "slimmedAddPileupInfo"
else:
    process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")
    _seq.append(process.offlineBeamSpot)
    _tracks, _gen, _pu = opts.trackSrc, "genParticles", "addPileupInfo"

process.globalCor = cms.EDProducer(
    "ResidualGlobalCorrectionMakerG4e",
    src=cms.InputTag(_tracks),
    fitFromGenParms=cms.bool(_derive),
    fitFromSimParms=cms.bool(False),
    fillTrackTree=cms.bool(True),
    fillGrads=cms.bool(_derive),
    fillJac=cms.bool(not _derive),
    fillRunTree=cms.bool(True),
    doGen=cms.bool(True),
    genParticles=cms.InputTag(_gen),
    pileupInfo=cms.InputTag(_pu),
    requireGen=cms.bool(True),
    genMatchPdgId=cms.int32(13),
    doSim=cms.bool(False),
    doMuons=cms.bool(False),
    doMuonAssoc=cms.bool(False),
    doTrigger=cms.bool(False),
    doRes=cms.bool(False),
    useIdealGeometry=cms.bool(bool(opts.useIdealGeometry)),
    bsConstraint=cms.bool(False),
    applyHitQuality=cms.bool(True),
    keepPixelEdgeHits=cms.bool(bool(opts.pixelClassHits)),
    pixelMinSizeX=cms.int32(1 if opts.pixelClassHits else 2),
    pixelHitClassCorrections=cms.bool(bool(opts.pixelClassHits)),
    perModuleBfield=cms.bool(True),
    globalMaterialModel=cms.bool(False),
    perStepFieldModes=cms.bool(False),
    skipHitlessSurfaces=cms.bool(False),
    materialGroupsFile=cms.string(""),
    scalarPotentialInitFile=cms.string(""),
    corFiles=cms.vstring(),
    triggers=cms.vstring(),
    trackParticleName=cms.string('mu'),
    MagneticFieldLabel=cms.string(""),
    nIters=cms.uint32(10),
    edmConvergence=cms.double(1e-5),
    outprefix=cms.untracked.string(opts.outprefix),
)
_seq.append(process.globalCor)

# The simulation's field: the default (unlabelled) MC field of the process.
process.geopro.MagneticFieldLabel = ""
process.Geant4ePropagator.MagneticFieldLabel = ""
process.Geant4ePropagator.ForCVH = cms.bool(True)
process.Geant4ePropagator.PropagationDirection = cms.string("anyDirection")
process.Geant4ePropagator.PropagationPtotLimit = cms.double(float(opts.propagationPtotLimit))
from TrackPropagation.Geant4e.cvhMasterESProducer_cfi import cvhMasterESProducer
process.cvhMasterESProducer = cvhMasterESProducer.clone(MagneticFieldLabel=cms.string(""))
process.globalCor.clampMomentumFloor = cms.double(1.25 * float(opts.propagationPtotLimit))
process.globalCor.maxMomentumStepFactor = cms.double(2.0)
process.globalCor.stepBacktracking = cms.bool(True)

_path = _seq[0]
for _m in _seq[1:]:
    _path = _path * _m
process.reconstruction_step = cms.Path(_path)
process.schedule = cms.Schedule(process.reconstruction_step)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)
