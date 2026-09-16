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
    # The luminous-region rows (see the beam-line block in
    # ResidualGlobalCorrectionMakerBase.h) are OFF here because the B vertex is DISPLACED by construction; the
    # pointing constraint is this channel's analogue of the beam line.
    bsConstraint = cms.bool(False),
    applyHitQuality = cms.bool(True),
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
    # minimum valid hits on the WEAKER leg; 0 = off (see the .cc). Default 8:
    # a thin leg is background (82-93 % of what it removes is `dup`/`unmatched`
    # by gen truth, at 0.9983 signal efficiency) and it is what every
    # non-finite mass-resolution export has in common.
    minLegHits = cms.int32(8),
    doMassConstraint = cms.bool(True),
    massConstraint = cms.double(3.0969),          # PDG J/psi mass, GeV
    massConstraintWidth = cms.double(9.29e-5),    # PDG J/psi natural width, GeV
    daughterParticleName1 = cms.string("mu"),
    daughterParticleName2 = cms.string("mu"),
    corFiles = cms.vstring(),
    MagneticFieldLabel = cms.string(''),
    outprefix = cms.untracked.string('globalcor_jpsik_dimuon'),
)
