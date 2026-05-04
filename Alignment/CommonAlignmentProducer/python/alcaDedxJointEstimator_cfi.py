import FWCore.ParameterSet.Config as cms

# Joint Harmonic-2 truncated dE/dx estimator (Ih) over both strip and
# pixel hits, equivalent to the standard CMSSW
# `dedxPixelAndStripHarmonic2T085` (defined in
# RecoTracker/DeDx/python/dedxEstimators_cff.py). Defined here under a
# unique label so the V0 ALCARECO skims (KS/Lambda/D*) can share a
# single producer instance without clashing with the upstream RECO
# sequence's identically-named module.
alcaDedxJointEstimator = cms.EDProducer('DeDxEstimatorProducer',
    tracks          = cms.InputTag('generalTracks'),
    estimator       = cms.string('genericTruncated'),
    fraction        = cms.double(-0.15),     # drop the lowest 15% of hits
    exponent        = cms.double(-2.0),      # Harmonic-2
    UseStrip        = cms.bool(True),
    UsePixel        = cms.bool(True),
    ShapeTest       = cms.bool(True),
    MeVperADCStrip  = cms.double(3.61e-06*265),
    MeVperADCPixel  = cms.double(3.61e-06),
    Reccord         = cms.string('SiStripDeDxMip_3D_Rcd'),
    ProbabilityMode = cms.string('Accumulation'),
    UseCalibration  = cms.bool(False),
    calibrationPath = cms.string(''),
)
