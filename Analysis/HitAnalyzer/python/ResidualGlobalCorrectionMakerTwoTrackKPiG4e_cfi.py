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
    doVtxConstraint=cms.bool(False),
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
