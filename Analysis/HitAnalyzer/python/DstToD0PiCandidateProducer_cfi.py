import FWCore.ParameterSet.Config as cms

DstToD0PiCandidateProducer = cms.EDProducer(
    'DstToD0PiCandidateProducer',
    src=cms.InputTag('ALCARECOTkAlDstToD0Pi'),
    kaonMass=cms.double(0.493677),
    pionMass=cms.double(0.139570),
    softPionMass=cms.double(0.139570),
    minD0Mass=cms.double(1.70),
    maxD0Mass=cms.double(2.00),
    minDstMass=cms.double(1.75),
    maxDstMass=cms.double(2.30),
    minDeltaMass=cms.double(0.135),
    maxDeltaMass=cms.double(0.160),
    pvalMin=cms.double(0.0),
    applyChargeFilter=cms.bool(True),
    useUnsignedCharge=cms.bool(True),
    charge=cms.int32(1),
)
