import FWCore.ParameterSet.Config as cms

DstToD0PiCandidateProducer = cms.EDProducer(
    'DstToD0PiCandidateProducer',
    src=cms.InputTag('ALCARECOTkAlDstToD0Pi'),
    kaonMass=cms.double(0.493677),
    pionMass=cms.double(0.139570),
    softPionMass=cms.double(0.139570),
    # Mass windows harmonized with the Stage-1 AlignmentThreeBodyDecayTrackSelector
    # cuts in ALCARECOTkAlDstToD0Pi_cff.py (CMSSW_10_6_17_patch1) so step-2 does
    # not re-admit configurations the ALCARECO already rejected.
    minD0Mass=cms.double(1.81483),     # was 1.70 -- tightened to D0 +/-50 MeV
    maxD0Mass=cms.double(1.91483),     # was 2.00
    minDstMass=cms.double(1.860),      # was 1.75
    maxDstMass=cms.double(2.160),      # was 2.30
    minDeltaMass=cms.double(0.14243),  # was 0.135 -- tightened to delta-m +/-3 MeV
    maxDeltaMass=cms.double(0.14843),  # was 0.160
    pvalMin=cms.double(0.0),
    applyChargeFilter=cms.bool(True),
    useUnsignedCharge=cms.bool(True),
    charge=cms.int32(1),
)
