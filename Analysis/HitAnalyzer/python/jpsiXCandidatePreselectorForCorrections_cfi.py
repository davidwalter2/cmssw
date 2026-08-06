import FWCore.ParameterSet.Config as cms

# Candidate preselection for the CVH grads path
# Driver that consumes it: test/runCvhBplusJpsiK.py.
#
# Defaults are the LOOSE NOMINAL: loosen the kinematics, keep the displacement.
# The displacement/vertex-quality cuts suppress combinatorial background without
# touching bachelor kinematics; the analysis-level kinematic cuts sculpt exactly
# the soft/forward bachelor phase space the channel exists to probe. The AlCaReco
# itself (ALCARECOTkAlJpsiX_cff.py) already imposes minBachelorPt = 0.1,
# maxBachelorEta = 2.5, minJpsiPt = 3.0, minMotherPt = 5.0, so these floors sit
# at the dataset's own limits rather than above them.
#
# The tight (Bmm5 analysis) variant is a labelled systematic, supplied by the
# driver's selTight=True switch -- not the nominal.
#
# All ValueMap tags MUST be wired to the RAW arm of JpsiXKinematicFitProducer.
# `ref*` is computed from CVH-refit tracks, i.e. from the corrections being
# fitted; selecting on it would make the selection a function of the answer.
jpsiXCandidatePreselectorForCorrections = cms.EDProducer(
    'JpsiXCandidatePreselectorForCorrections',
    src=cms.InputTag('ALCARECOTkAlJpsiXBPlusResonances'),

    # --- raw-arm ValueMaps -------------------------------------------------
    fitOk=cms.InputTag('bplusFit', 'rawFitOk'),
    fitMass=cms.InputTag('bplusFit', 'rawFitMass'),
    fitVtxProb=cms.InputTag('bplusFit', 'rawFitVtxProb'),
    dimuonVtxProb=cms.InputTag('bplusFit', 'rawDimuonVtxProb'),
    dimuonAlphaBS=cms.InputTag('bplusFit', 'rawDimuonAlphaBS'),
    dimuonSl3d=cms.InputTag('bplusFit', 'rawDimuonSl3d'),

    # --- mother mass window ------------------------------------------------
    # +/-100 MeV about the PDG B+ mass. Its job is to BOUND THE PULL, not to
    # provide sidebands (no sideband subtraction is applied): at the
    # measured lever arm S_K = d ln m_B / d ln p_K = 0.248, a 100 MeV window
    # caps the momentum displacement the B-mass row can impose on a background
    # candidate at (0.1/5.279)/0.248 ~ 8%.
    massCentre=cms.double(5.27934),
    massHalfWindow=cms.double(0.100),

    # --- displacement / vertex quality: KEPT at analysis strength ----------
    fitVtxProbMin=cms.double(0.1),
    dimuonVtxProbMin=cms.double(0.1),
    dimuonAlphaBSMax=cms.double(0.4),
    dimuonSl3dMin=cms.double(4.0),

    # --- kinematics: LOOSENED to the AlCaReco floors -----------------------
    muonPtMin=cms.double(3.0),
    muonEtaMax=cms.double(2.4),
    bachelorPtMin=cms.double(0.1),
    bachelorPtMax=cms.double(1e4),   # analysis uses 8.0; effectively off here
    bachelorEtaMax=cms.double(2.4),
    dimuonPtMin=cms.double(0.0),   # analysis uses 7.0; dropped nominally

    # --- correctness, not tuning -------------------------------------------
    # At ~1.7 candidates/event the same dimuon's hits would otherwise enter the
    # summed global Hessian repeatedly and the sum would no longer be the
    # information of an independent sample.
    bestPerEvent=cms.bool(True),
)

# Bmm5 get_bkmm_selections values, for the systematic variation. Retained
# explicitly so the comparison is reproducible: on our sample these keep 18% of
# bachelors (median bachelor pT ~0.53 GeV, 82% below 1 GeV).
jpsiXCandidatePreselectorForCorrectionsTight = jpsiXCandidatePreselectorForCorrections.clone(
    muonPtMin=4.0,
    muonEtaMax=1.4,
    bachelorPtMin=1.0,
    bachelorPtMax=8.0,
    bachelorEtaMax=1.4,
    dimuonPtMin=7.0,
)
