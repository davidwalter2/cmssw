import FWCore.ParameterSet.Config as cms
import os

# Generic dimuon (two-track) CVH refit for the WMass custom NanoAOD, run on
# MINIAOD-derived muon-track pairs (diMuonTrackVertexCandidates). Resonance-
# agnostic: a plain common-vertex fit, no mass constraint, so the SAME config
# serves Z / J/psi / Upsilon -- the resonance is selected only by the candidate
# producer's mass window. produceValueMaps=True emits per-candidate EDM
# ValueMaps (keyed to srcCandidates) that the Dimuon NanoAOD table consumes;
# kinematics always, the global-fit payload (globalIdxs/jacRef/jacMass/Hessian)
# only when fillGradsFactored is turned on.
#
# The Geant4 master is the shared EventSetup product (cvhMasterESProducer);
# MagneticFieldLabel / scalarPotentialInitFile are set at customise time by
# nano_cff.setup3DFieldForRefit (ScalarPot3D for data, default field for MC).
_materialGroups50 = os.path.join(os.environ.get('CMSSW_BASE', ''),
                                 'src/Analysis/HitAnalyzer/data/materialGroups50.txt')

ResidualGlobalCorrectionMakerDiMuonG4e = cms.EDProducer(
    'ResidualGlobalCorrectionMakerTwoTrackG4e',
    src = cms.InputTag('tracksfrommuons'),
    srcCandidates = cms.InputTag('diMuonTrackVertexCandidates'),
    dedxSourceTracks   = cms.InputTag(''),
    dedxHarmonic2      = cms.InputTag(''),
    dedxPixelHarmonic2 = cms.InputTag(''),
    dedxAllHarmonic2   = cms.InputTag(''),
    fitFromGenParms = cms.bool(False),
    fitFromSimParms = cms.bool(False),
    fillTrackTree = cms.bool(False),
    fillGrads = cms.bool(False),
    # Flag: gates the (large) global-fit ValueMap payload
    # (globalIdxs / jacRef / jacMass / factored Hessian). Off by default ->
    # only the dimuon kinematics are stored.
    fillGradsFactored = cms.untracked.bool(False),
    fillJac = cms.bool(False),
    fillRunTree = cms.bool(False),
    doGen = cms.bool(False),
    genParticles = cms.InputTag('genParticles'),
    pileupInfo = cms.InputTag('addPileupInfo'),
    doSim = cms.bool(False),
    requireGen = cms.bool(False),
    doMuons = cms.bool(False),
    muons = cms.InputTag('muons'),
    doMuonAssoc = cms.bool(False),
    doTrigger = cms.bool(False),
    triggers = cms.vstring(),
    doL1Trigger = cms.bool(False),
    l1Results = cms.InputTag('gtDigis', '', 'RECO'),
    l1Triggers = cms.vstring(),
    doRes = cms.bool(False),
    useIdealGeometry = cms.bool(False),
    bsConstraint = cms.bool(False),
    applyHitQuality = cms.bool(True),
    keepPixelEdgeHits = cms.bool(False),
    pixelMinSizeX = cms.int32(2),
    doVtxConstraint = cms.bool(True),
    # Minimum size of a pair, required BEFORE the fit (see
    # ResidualGlobalCorrectionMakerTwoTrackG4e.cc). ndof = nvalid +
    # nvalidpixel - 10 (+3 beamspot, +1 pointing, +1 vertex constraint): one
    # measurement coordinate per strip hit, two per pixel hit, against the ten
    # state parameters the common vertex costs. minNdof = 1 therefore requires
    # more than NINE measurement coordinates with the vertex constraint on and
    # more than TEN with it off -- at ndof == 0 the fit is exactly determined
    # (chi2 identically zero, chi2/ndof undefined) and the factored-Hessian
    # export indexes past the end of its eigenvalue vector and aborts the
    # process. minPairHits is the same requirement read on VALID HITS rather
    # than on measurement coordinates; -1 = auto = 10 (constraint on) / 11
    # (off). 0 disables either.
    minNdof = cms.int32(1),
    minPairHits = cms.int32(-1),
    # minimum valid hits on the WEAKER leg; 0 = off (see the .cc).
    minLegHits = cms.int32(0),
    doMassConstraint = cms.bool(False),
    massConstraint = cms.double(91.1876),        # unused (no mass constraint)
    massConstraintWidth = cms.double(2.4952),
    daughterParticleName1 = cms.string('mu'),
    daughterParticleName2 = cms.string('mu'),
    corFiles = cms.vstring(),
    MagneticFieldLabel = cms.string(''),
    scalarPotentialInitFile = cms.string(''),
    nIters = cms.uint32(10),
    edmConvergence = cms.double(1e-5),
    globalMaterialModel = cms.bool(True),
    perStepFieldModes = cms.bool(True),
    skipHitlessSurfaces = cms.bool(True),
    materialGroupsFile = cms.string(_materialGroups50),
    produceValueMaps = cms.bool(True),
    outprefix = cms.untracked.string('globalcor_dimuon'),
)
