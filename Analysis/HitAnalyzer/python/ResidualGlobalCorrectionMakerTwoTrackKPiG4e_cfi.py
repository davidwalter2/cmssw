import FWCore.ParameterSet.Config as cms

globalCorD0 = cms.EDProducer(
    'ResidualGlobalCorrectionMakerTwoTrackG4e',
    src=cms.InputTag('ALCARECOTkAlDstToD0Pi'),
    # Default fast path: persisted D* candidates (positive-charge-first
    # convention; daughter[0] = K, daughter[1] = pi from D0). See PiPi
    # cfi for details. Requires stage-1 to provide the *Resonances
    # collection (ThreeBodyDecayCandidateProducer + remapper); leave the
    # legacy in-selector pair-finding chain disabled at stage-1 for this
    # to apply.
    srcCandidates=cms.InputTag('ALCARECOTkAlDstToD0PiResonances'),
    dedxSourceTracks  =cms.InputTag('ALCARECOTkAlDstToD0Pi'),
    dedxHarmonic2     =cms.InputTag('ALCARECOTkAlDstToD0PiDeDxHarmonic2'),
    dedxPixelHarmonic2=cms.InputTag('ALCARECOTkAlDstToD0PiDeDxPixelHarmonic2'),
    dedxAllHarmonic2  =cms.InputTag('ALCARECOTkAlDstToD0PiDeDxAllHarmonic2'),
    fitFromGenParms=cms.bool(False),
    fitFromSimParms=cms.bool(False),
    fillTrackTree=cms.bool(True),
    fillGrads=cms.bool(False),
    fillJac=cms.bool(False),
    fillRunTree=cms.bool(False),
    doGen=cms.bool(False),
    genParticles=cms.InputTag('genParticles'),
    pileupInfo=cms.InputTag('addPileupInfo'),
    doSim=cms.bool(False),
    requireGen=cms.bool(False),
    doMuons=cms.bool(False),
    muons=cms.InputTag('muons'),
    doMuonAssoc=cms.bool(False),
    doTrigger=cms.bool(False),
    triggers=cms.vstring(),
    doL1Trigger=cms.bool(False),
    l1Results=cms.InputTag('gtDigis', '', 'RECO'),
    l1Triggers=cms.vstring(),
    doRes=cms.bool(False),
    useIdealGeometry=cms.bool(True),
    bsConstraint=cms.bool(False),
    applyHitQuality=cms.bool(True),
    doVtxConstraint=cms.bool(True),
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
    doMassConstraint=cms.bool(False),
    massConstraint=cms.double(1.86483),
    massConstraintWidth=cms.double(1.605e-12),  # natural width Gamma = hbar/tau (D0), GeV
    # Per-daughter Geant4 particle base (track[0] = K, track[1] = pi).
    # Mass + mass uncertainty are looked up from the PDG table in
    # Analysis/HitAnalyzer/interface/ParticleProperties.h .
    daughterParticleName1=cms.string("kaon"),
    daughterParticleName2=cms.string("pi"),
    corFiles=cms.vstring(),
    MagneticFieldLabel=cms.string(''),
    # Scalar-potential B-field correction (absolute-field model); path to a
    # coefficient dump file produced by mfs/dump_coeffs_for_cmssw.py.
    scalarPotentialInitFile=cms.string(''),
    outprefix=cms.untracked.string('globalcor_d0')
)
