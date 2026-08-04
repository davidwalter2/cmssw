## CVH gen-closure refit driver for the inclusive B->J/psi+X MC ALCARECO
## (2016 postVFP, produced with CMSSW_10_6_20_patch1, split=1, gen kept).
## Adapted from runCvhJpsi.py (data driver) for the pixel edge / single-
## column hit study: fits J/psi->mumu candidates with fitFromGenParms=True
## (reference state frozen to gen kinematics -> no mass constraint, no
## weak modes; validated unbiased in the past), so hit-pathology biases
## can be measured directly against gen truth.
##
## Default input collections are the ALCARECOTkAlJpsiX labels of the MC
## (tracks + JpsiOnlyResonances candidates); the Golden-JSON and HLT
## filters of the data driver are dropped / made opt-in.
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os

from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
from Configuration.AlCa.GlobalTag import GlobalTag

opts = VarParsing.VarParsing('analysis')
opts.register('input', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'comma-separated absolute paths or root:// URLs of ALCARECO files')
opts.register('inputFileList', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'text file with one input path per line (the MC has ~10 events '
              'per file, so runs typically need many files); combined with '
              'input= if both are given')
opts.register('fitFromGenParms', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'freeze the 10 vertex/kinematic reference parameters to the '
              'gen-muon values (gen-closure mode, default True)')
opts.register('applyHltFilter', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'require one of the J/psi HLT paths (default False for MC '
              'closure; the ALCARECO selection already ran)')
opts.register('useLegacyPairLoop', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'ignore the JpsiOnlyResonances candidates and pair all tracks '
              'in the module (legacy fallback; default False)')
opts.register('fillHitDiagnostics', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'store per-pixel-hit diagnostic branches (hitdiag_*): local '
              'residuals + side-resolved pathology class')
opts.register('deweightPathoHits', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'deweight pathological pixel hits (x1e-6) instead of using '
              'them: keeps surface+state so hitdiag residuals are unbiased '
              'w.r.t. the rest of the fit; combine with keepPixelEdgeHits='
              'True pixelMinSizeX=1 and fitFromGenParms=True')
opts.register('pixelHitClassCorrections', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'register the per-pixel-module pathology-class correction '
              'parameters (parmtypes 16-21: edge-x-mean/diff, edge-y-mean/'
              'diff, sizeX1, sizeY1) and emit their Jacobian columns; use '
              'with keepPixelEdgeHits=True pixelMinSizeX=1')
opts.register('corFile', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'optional correction file (parmtree/x, one entry per catalog '
              'parameter) applied via corparms_, e.g. the fitted class '
              'corrections from write_classcorr_corfile.py')
opts.register('pixelLorentzParam', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'physics mode: replace edge-x-mean/sizeX1 by a per-module '
              'delta-tan(thetaL) parameter (parmtype 22) with a Jacobian '
              'column on EVERY valid pixel hit (weights: size-1 = 1, '
              'x-edge = lorentzWedge, regular = lorentzWclean)')
opts.register('lorentzWclean', 0.67, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'dtanLA response weight of regular (clean) pixel hits '
              '(measured: digitizer twin-sample study)')
opts.register('lorentzWsize1', 0.02, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'dtanLA response weight of size-1 pixel hits (measured ~0: '
              'pixel-center quantisation)')
opts.register('lorentzWedge', 0.45, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'dtanLA response weight of x-edge pixel hits (low-stats '
              'measurement, +-0.4)')
opts.register('injectLorentzTan', 0., VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'closure test: inject a true Lorentz-angle mismatch of this '
              'size (shifts every valid pixel hit local-x by t/2*w*value)')
opts.register('injectLorentzWclean', -999., VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'clean-hit response weight used for the INJECTION (response-'
              'model error study); < -900 = same as lorentzWclean')
opts.register('nEvents', 500, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, 'number of events to process (-1 = all)')
opts.register('fillJac', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool, 'store per-track Jacobians')
opts.register('fillGrads', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool, 'store per-event gradient + packed Hessian')
opts.register('fillGradsFactored', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'store per-event gradient + low-rank factored Hessian (H = B^T B)')
opts.register('doTrigger', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'read TriggerResults::HLT (off for private samples without HLT)')
opts.register('doVtxConstraint', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'apply the common-vertex constraint in the two-track fit')
opts.register('doSimHits', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'read tracker PSimHits (input must retain them)')
opts.register('fitSimHitPositions', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'rung-E closure: fit simulated hit positions (measured '
              'coordinates only, covariances unchanged); needs doSimHits=True')
opts.register('propagationPtotLimit', 0.2, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'G4e propagation momentum floor [GeV]; cfi default was 1.0')
opts.register('doRes', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'register resolution families and export the per-candidate '
              'mass-CF ingredients (dV blocks, step records, mass-projected '
              'influence weights)')
opts.register('trackSrc', 'ALCARECOTkAlJpsiX', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'input track collection (ALCARECOTkAlJpsiMuMu for the standard '
              'TkAl ALCARECO of the JPsiToMuMu MC; pair with useLegacyPairLoop=True)')
opts.register('doMassConstraint', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool, 'apply J/psi mass constraint in the two-track fit')
opts.register('useIdealGeometry', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'Default False (btojpsik option (B), aligned geometry from GT). Set True '
              'only for the AN Stage-1 broken baseline (Stage-2 corrections not applied '
              'here). See openspec/finalize-cvh-producer-15-0-19.')
opts.register('useScalarPot3D', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'use the spherical-harmonic scalar-potential field  '
              'in the CVH refit (default; only model supported in this port)')
opts.register('useOpera3D', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'use the full 3D TOSCA volumetric grid (160812) as the baseline '
              'field for the propagator + geopro + globalCor; takes precedence '
              'over useScalarPot3D when True.')
opts.register('scalarPot3DInitFile', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'coefficient dump file produced by mfs/dump_coeffs_for_cmssw.py. '
              'Always required: the residual-correction maker uses it to register '
              'parmtype-14 modes and seed their initial coefficients. Also reused '
              'as the field producer init file when useScalarPot3D=True.')
opts.register('runFDClosure', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'Numerical-FD closure of the per-mode chain rule '
              '(debug; runs once on the first chain-rule site)')
opts.register('epsilonFDClosure', 1e-4, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'eps for the FD closure (used as eps * dB_perMode for each test mode)')
opts.register('numberOfThreads', 1, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'framework numberOfThreads (numberOfStreams follows the same value)')
opts.register('debugPerIterDump', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'per-iteration debug trace (tree vectors + stdout dbgSeed/dbgIter '
              'lines); use together with eventsToProcess on a few events')
opts.register('eventsToProcess', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'comma-separated run:event list to select specific events '
              '(e.g. 278769:15462343,278769:16101980); empty = all')
opts.register('nIters', 10, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'Gauss-Newton iteration cap per constraint phase (default 10 = baseline)')
opts.register('edmConvergence', 1e-5, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'EDM convergence threshold on the reference-state block (default 1e-5; '
              '0 disables early stopping, e.g. for per-iteration trajectory studies)')
opts.register('keepPixelEdgeHits', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'keep pixel hits whose cluster touches the sensor boundary '
              '(isOnEdge) in the fit instead of demoting them to inactive; '
              'the pixelMinSizeX CPE-quality cut applies independently '
              '(default False = baseline)')
opts.register('pixelMinSizeX', 2, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'minimum pixel cluster size in x for a hit to stay in the fit '
              '(default 2 = baseline sizeX>1 cut; 1 admits all clusters)')
_defaultGroupsFile = os.path.join(os.environ.get('CMSSW_BASE', ''),
                                  'src/Analysis/HitAnalyzer/data/materialGroups50.txt')
opts.register('materialGroupsFile', _defaultGroupsFile, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'global material model grouping-tier rules file '
              '(Analysis/HitAnalyzer/data/materialGroups{50,100}.txt); empty = off')
opts.register('perStepFieldModes', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'apply the scalar-potential correction and attribute the per-mode '
              'derivatives per Geant4 step instead of piecewise-constant per leg '
              '(leg-structure-free field attribution; default True)')
opts.register('skipHitlessSurfaces', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'drop hitless module surfaces (dead-module placeholders, '
              'quality-demoted hits) from the fit; propagation goes hit to hit. '
              'Default True; effective only with globalMaterialModel=True '
              '(auto-disabled otherwise)')
opts.register('globalMaterialModel', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'replace the per-module material parameters (parmtype 7) with the '
              'parmtype-15 global material groups of materialGroupsFile '
              '(exclusive switch). Default True (tier-50 groups file from the '
              'release); set False for the legacy per-module parameterisation')
opts.register('propagationDirection', 'anyDirection', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'Geant4ePropagator PropagationDirection. "anyDirection" (default) '
              'picks forward/backward per leg from the target-plane geometry, '
              'recovering legs whose target plane is marginally behind the '
              'state (runaway-leg failure mode); "alongMomentum" is the '
              'legacy forward-only behaviour (bit-identical for all fits '
              'that do not fail with it).')
opts.parseArguments()
if not opts.scalarPot3DInitFile:
    raise SystemExit(
        "scalarPot3DInitFile=<path> is required (coefficient dump file): "
        "the basis evaluator in globalCor needs it for chain-rule columns "
        "even when useOpera3D=True swaps the baseline field model.")

JPSI_TRIGGERS = [
    "HLT_Dimuon0_Jpsi_Muon",
    "HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing",
    "HLT_Dimuon0er16_Jpsi_NoVertexing",
    "HLT_Dimuon10_Jpsi_Barrel",
    "HLT_Dimuon13_PsiPrime",
    "HLT_Dimuon16_Jpsi",
    "HLT_Dimuon20_Jpsi",
    "HLT_Dimuon8_PsiPrime_Barrel",
    "HLT_DoubleMu4_3_Bs",
    "HLT_DoubleMu4_3_Jpsi_Displaced",
    "HLT_DoubleMu4_JpsiTrk_Displaced",
    "HLT_DoubleMu4_PsiPrimeTrk_Displaced",
    "HLT_Mu7p5_Track2_Jpsi",
    "HLT_Mu7p5_Track3p5_Jpsi",
    "HLT_Mu7p5_Track7_Jpsi",
]

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

# Conditions the MC was produced with (CMSSW_10_6_20_patch1 production
# chain) -- alignment/CPE/beamspot consistent with the simulated detector,
# which is what a gen-closure fit must use.
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
assert _paths, "must set input=<paths> and/or inputFileList=<file> on the cmsRun command line"
# Accept local paths (prepend "file:") or xrootd URLs as-is.
_urls = [p if p.startswith(("root://", "file:")) else "file:" + p for p in _paths]
process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(*_urls),
    secondaryFileNames=cms.untracked.vstring(),
    # The condor MC production has a small tail of corrupt files (garbled
    # embedded provenance -> FormatIncompatibility at readFile_); skip
    # them instead of aborting the whole many-file job. Skipped files are
    # reported in the log. NOTE: this does NOT cover the
    # FormatIncompatibility case -- pre-scan the filelist (see
    # calibration_studies/pixelhits) and exclude those files.
    skipBadFiles=cms.untracked.bool(True),
    # Every condor job numbers its events from the same (run=1, lumi=1)
    # range, so distinct physics events collide in (run, lumi, event) and
    # the default duplicate check silently drops most of the sample.
    duplicateCheckMode=cms.untracked.string('noDuplicateCheck'),
)

if opts.eventsToProcess:
    process.source.eventsToProcess = cms.untracked.VEventRange(
        *[s.strip() for s in opts.eventsToProcess.split(',') if s.strip()])

process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(int(opts.numberOfThreads)),
    numberOfStreams=cms.untracked.uint32(int(opts.numberOfThreads)),
    numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1),
)

# Per-stream CLHEP engine for the residual-maker. The Tier-3 MT path calls
# setG4RandomEngineForStream() at the top of every produce() to wire this
# engine into Geant4's thread-local RNG. Reproducible across thread counts
# because the framework derives per-stream seeds deterministically from the
# initialSeed below + stream index.
process.RandomNumberGeneratorService.globalCor = cms.PSet(
    initialSeed=cms.untracked.uint32(123456789),
    engineName=cms.untracked.string('HepJamesRandom'),
)

# Reduce log spam (every 100 events instead of every event).
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")

# HLT pre-filter: drop events that don't pass any of the J/psi paths we
# also store decisions for. Saves the Geant4e/CVH cost on triggers we'd
# never analyse. throw=False so the filter tolerates menu changes across
# eras (any path missing in a given menu is silently skipped).
process.hltFilter = cms.EDFilter(
    "HLTHighLevel",
    HLTPaths=cms.vstring(*[t + "_v*" for t in JPSI_TRIGGERS]),
    eventSetupPathsKey=cms.string(""),
    andOr=cms.bool(True),     # OR over the path list
    throw=cms.bool(False),
    TriggerResultsTag=cms.InputTag("TriggerResults", "", "HLT"),
)

process.globalCor = cms.EDProducer(
    "ResidualGlobalCorrectionMakerTwoTrackG4e",
    src=cms.InputTag(opts.trackSrc),
    fitFromGenParms=cms.bool(bool(opts.fitFromGenParms)),
    fitFromSimParms=cms.bool(False),
    fillTrackTree=cms.bool(True),
    fillGrads=cms.bool(bool(opts.fillGrads)),
    # Low-rank factored Hessian storage (H = B^T B, nRank x nParms):
    # ~9x smaller than hesspackedv at 360 field modes; see
    # ResidualGlobalCorrectionMakerTwoTrackG4e.cc for the rank argument.
    fillGradsFactored=cms.untracked.bool(bool(opts.fillGradsFactored)),
    fillJac=cms.bool(bool(opts.fillJac)),
    fillRunTree=cms.bool(True),
    doGen=cms.bool(True),
    genParticles=cms.InputTag("genParticles"),
    pileupInfo=cms.InputTag("addPileupInfo"),
    doSim=cms.bool(bool(opts.doSimHits)),
    fitSimHitPositions=cms.untracked.bool(bool(opts.fitSimHitPositions)),
    # Gen matching (dR < 0.1, same charge, status-1 muons) is required both
    # to anchor fitFromGenParms and to reject combinatorial pairs.
    requireGen=cms.bool(True),
    doMuons=cms.bool(False),
    doMuonAssoc=cms.bool(False),
    doTrigger=cms.bool(bool(opts.doTrigger)),
    doRes=cms.bool(bool(opts.doRes)),
    useIdealGeometry=cms.bool(bool(opts.useIdealGeometry)),
    bsConstraint=cms.bool(False),
    applyHitQuality=cms.bool(True),
    keepPixelEdgeHits=cms.bool(bool(opts.keepPixelEdgeHits)),
    pixelMinSizeX=cms.int32(int(opts.pixelMinSizeX)),
    fillHitDiagnostics=cms.bool(bool(opts.fillHitDiagnostics)),
    deweightPathoHits=cms.bool(bool(opts.deweightPathoHits)),
    pixelHitClassCorrections=cms.bool(bool(opts.pixelHitClassCorrections)),
    pixelLorentzParam=cms.bool(bool(opts.pixelLorentzParam)),
    lorentzWclean=cms.double(float(opts.lorentzWclean)),
    lorentzWsize1=cms.double(float(opts.lorentzWsize1)),
    lorentzWedge=cms.double(float(opts.lorentzWedge)),
    injectLorentzTan=cms.double(float(opts.injectLorentzTan)),
    injectLorentzWclean=cms.double(float(opts.injectLorentzWclean)),
    doVtxConstraint=cms.bool(bool(opts.doVtxConstraint)),
    doMassConstraint=cms.bool(bool(opts.doMassConstraint)),
    massConstraint=cms.double(3.0969),
    massConstraintWidth=cms.double(1e-5),
    corFiles=cms.vstring(*( [opts.corFile] if opts.corFile else [] )),
    triggers=cms.vstring(*JPSI_TRIGGERS),
    MagneticFieldLabel=cms.string(""),
    # Scalar-potential B-field correction (parmtype-14, absolute-field
    # model). Initial coefficients + basis structure are loaded from a
    # coefficient dump file (mfs/dump_coeffs_for_cmssw.py output). The
    # dump's mode count determines nFieldModes -- use a 50-mode
    # ("custom50": lphi5-base + l=6,m=1; see mfs/CLAUDE.md) dump to keep
    # the per-event Hessian workspace small in MT runs.
    scalarPotentialInitFile=cms.string(opts.scalarPot3DInitFile),
    # Numerical-FD closure (debug only).
    runFDClosure=cms.bool(bool(opts.runFDClosure)),
    epsilonFDClosure=cms.double(float(opts.epsilonFDClosure)),
    debugPerIterDump=cms.bool(bool(opts.debugPerIterDump)),
    nIters=cms.uint32(int(opts.nIters)),
    edmConvergence=cms.double(float(opts.edmConvergence)),
    materialGroupsFile=cms.string(opts.materialGroupsFile),
    globalMaterialModel=cms.bool(bool(opts.globalMaterialModel)),
    perStepFieldModes=cms.bool(bool(opts.perStepFieldModes)),
    skipHitlessSurfaces=cms.bool(bool(opts.skipHitlessSurfaces) and bool(opts.globalMaterialModel)),
    outprefix=cms.untracked.string("globalcor"),
    # MT G4Error master: GlobalCache config for CvhMasterThread. The master
    # spawns a dedicated thread in initializeGlobalCache that builds DDDWorld
    # + master magnetic field BEFORE any TBB worker starts. Each per-stream
    # CvhWorker then attaches per-thread G4 state to it on first produce().
    # This replaces the geopro side-effect dependency that blocked
    # numberOfThreads >= 2 previously.
    #
    # Narrowed to muons -- this runner only propagates J/psi -> mu mu
    # daughters, so the rest of the canonical CVH particle set
    # (gamma, e+-, pi+-, K+-, p, anti_p) is skipped at physics-list
    # construction. Saves the per-thread ProcessManager + process
    # allocations for ~9 unused particles.
    CvhMaster=CvhMasterPSet.clone(Particles=cms.vstring("mu+", "mu-")),
)

# Candidate-driven pair building: use the persisted J/psi->mumu candidates
# of the TkAlJpsiX ALCARECO instead of the all-pairs legacy loop (the MC
# track collection also contains the other B daughters, e.g. the kaon).
if not opts.useLegacyPairLoop:
    process.globalCor.srcCandidates = cms.InputTag("ALCARECOTkAlJpsiXJpsiOnlyResonances")

# Bring up the labelled 3D field producer and rewire the consumers
# present in this driver (geopro, Geant4ePropagator, and our
# globalCor analyzer). Uses the scalar-potential ScalarPot3D model
# from scalar-potential field model. Independent from
# nano_cff.setup3DFieldForRefit (which assumes the full set of seven
# CVH-side consumers from the NanoAOD configuration).
if opts.useOpera3D:
    # Use the full 3D TOSCA volumetric grid (160812) as the baseline field
    # for the propagator / geopro / globalCor instead of the scalar-potential
    # ScalarPot3D model. Provided as an alternative field-model option for
    # cross-checks and B-field studies.
    from MagneticField.Engine.volumeBasedMagneticField_160812_cfi import \
        VolumeBasedMagneticFieldESProducer as Opera3DMagneticFieldProducer
    from MagneticField.Engine.volumeBasedMagneticField_160812_cfi import magfield as MagneticFieldGeometry
    process.magfield = MagneticFieldGeometry
    process.es_prefer_magfield_cvhrefit = cms.ESPrefer("XMLIdealGeometryESSource", "magfield")
    process.Opera3DMagneticFieldProducer = Opera3DMagneticFieldProducer
    fieldlabel = "grid_160812_3_8t"
    process.Opera3DMagneticFieldProducer.label = fieldlabel
    process.Opera3DMagneticFieldProducer.useParametrizedTrackerField = cms.bool(False)
    # Route the labelled field into the CPEs as in the production data refit
    # (nano_cff.nanoAOD_customizeData) and the 10_6 cross-release driver:
    # the Lorentz-drift in the hit re-evaluation then uses the same field.
    for _cpe in ("stripCPEESProducer", "StripCPEfromTrackAngleESProducer",
                 "siPixelTemplateDBObjectESProducer", "templates"):
        if hasattr(process, _cpe):
            getattr(process, _cpe).MagneticFieldLabel = fieldlabel
elif not opts.useScalarPot3D:
    raise RuntimeError(
        "useScalarPot3D=False is no longer supported; the legacy non-thread-safe "
        "wrapper class is not part of this port. Use the ScalarPot3D model.")
else:
    if not opts.scalarPot3DInitFile:
        raise RuntimeError(
            "useScalarPot3D=True requires scalarPot3DInitFile to point "
            "at a coefficient dump file produced by mfs/dump_coeffs_for_cmssw.py")
    from MagneticField.ParametrizedEngine.parametrizedMagneticField_ScalarPot3D_cfi \
        import ParametrizedMagneticFieldProducer as ScalarPot3DMagneticFieldProducer
    process.ScalarPot3DMagneticFieldProducer = ScalarPot3DMagneticFieldProducer.clone()
    process.ScalarPot3DMagneticFieldProducer.parameters.InitFile = opts.scalarPot3DInitFile
    fieldlabel = "ScalarPot3DMf"
    process.ScalarPot3DMagneticFieldProducer.label = fieldlabel
process.geopro.MagneticFieldLabel = fieldlabel
process.Geant4ePropagator.MagneticFieldLabel = fieldlabel
# The Geant4 master is now the shared EventSetup product from
# cvhMasterESProducer (CvhMasterRecord), consumed by globalCor via esConsumes.
# It builds its master G4 field via SimG4Core's FieldBuilder on top of the same
# labelled magnetic field the propagator consumes.
from TrackPropagation.Geant4e.cvhMasterESProducer_cfi import cvhMasterESProducer
process.cvhMasterESProducer = cvhMasterESProducer.clone()
process.cvhMasterESProducer.MagneticFieldLabel = cms.string(fieldlabel)
# Activate the CVH-specific propagator path: instantiates the custom fluct
# (G4UniversalFluctuationForExtrapolator) and routes its table pointer via
# SetParticleAndCharge. Without this, computeErrorIoni dereferences a null
# fluct->table on the first event.
process.Geant4ePropagator.ForCVH = cms.bool(True)
process.Geant4ePropagator.PropagationDirection = cms.string(opts.propagationDirection)
process.Geant4ePropagator.PropagationPtotLimit = cms.double(float(opts.propagationPtotLimit))
process.globalCor.MagneticFieldLabel = cms.string(fieldlabel)

# geopro is removed: CvhMasterThread (residual-maker GlobalCache) now
# owns the G4 world / master magnetic field setup in an MT-safe way.
# See TrackPropagation/Geant4e/{interface,src}/CvhMaster*.
if opts.applyHltFilter:
    process.reconstruction_step = cms.Path(
        process.hltFilter * process.offlineBeamSpot * process.globalCor
    )
else:
    process.reconstruction_step = cms.Path(
        process.offlineBeamSpot * process.globalCor
    )
process.schedule = cms.Schedule(process.reconstruction_step)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
