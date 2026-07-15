import FWCore.ParameterSet.Config as cms

# Stage-2 CVH refit of the J/psi -> mu+ mu- subsystem of B+ -> J/psi K+,
# composed with the bachelor single-track refit (see
# ResidualGlobalCorrectionMakerJpsiKSingleTrackKaonG4e_cfi.py).
#
# Differences vs. the standalone J/psi cfi
# (`ResidualGlobalCorrectionMakerTwoTrackJpsiMuMuG4e_cfi.py`):
#
#   - `srcCandidates` points at the JpsiKCandidateSplitter's `dimuon`
#     output (one VCC per surviving B+, daughters = the two original J/psi
#     muon-RCCs);
#   - `bCandIdxSrc` points at the splitter's `dimuonBCandIdx` companion
#     vector<int>, enabling the additive `bCandIdx` branch in the maker
#     so the offline join can pair the dimuon refit row to the kaon
#     refit row on (run, lumi, event, bCandIdx);
#   - `doMassConstraint = True` (the dimuon mass is pulled to the PDG
#     J/psi mass inside the joint two-track CVH vertex fit);
#   - `fillJac = True` (per-track Jacobians for the downstream GBL
#     handoff, even though 2a stops at the m(mu mu K) closure).
#
# `src` keeps pointing at the persisted track collection
# (`ALCARECOTkAlJpsiX` -- the cross-release name); the maker iterates the
# fast path over `srcCandidates`, reaching tracks via each VCC daughter's
# `RecoChargedCandidate::track()` ref into `src`.
globalCorJpsiK = cms.EDProducer(
    'ResidualGlobalCorrectionMakerTwoTrackG4e',
    src = cms.InputTag('ALCARECOTkAlJpsiX'),
    srcCandidates = cms.InputTag('jpsiKCandidateSplitter', 'dimuon'),
    bCandIdxSrc = cms.InputTag('jpsiKCandidateSplitter', 'dimuonBCandIdx'),
    dedxSourceTracks   = cms.InputTag('ALCARECOTkAlJpsiX'),
    dedxHarmonic2      = cms.InputTag(''),
    dedxPixelHarmonic2 = cms.InputTag(''),
    dedxAllHarmonic2   = cms.InputTag(''),
    fitFromGenParms = cms.bool(False),
    fitFromSimParms = cms.bool(False),
    fillTrackTree = cms.bool(True),
    fillGrads = cms.bool(False),
    fillJac = cms.bool(True),
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
    # False by default: True together with the default empty corFiles is the
    # broken Stage-1-only configuration (ideal geometry, no corrections)
    # the maker itself warns about. The driver overrides per leg.
    useIdealGeometry = cms.bool(False),
    bsConstraint = cms.bool(False),
    applyHitQuality = cms.bool(True),
    doVtxConstraint = cms.bool(False),
    doMassConstraint = cms.bool(True),
    massConstraint = cms.double(3.0969),          # PDG J/psi mass, GeV
    massConstraintWidth = cms.double(9.29e-5),    # PDG J/psi natural width, GeV
    daughterParticleName1 = cms.string("mu"),
    daughterParticleName2 = cms.string("mu"),
    corFiles = cms.vstring(),
    MagneticFieldLabel = cms.string(''),
    outprefix = cms.untracked.string('globalcor_jpsik_dimuon'),
)
