import FWCore.ParameterSet.Config as cms

import DPGAnalysis.Skims.skim_detstatus_cfi
ALCARECOTkAlDstToD0PiDCSFilter = DPGAnalysis.Skims.skim_detstatus_cfi.dcsstatus.clone(
    DetectorType = cms.vstring('TIBTID','TOB','TECp','TECm','BPIX','FPIX',
                               'DT0','DTp','DTm','CSCp','CSCm'),
    ApplyFilter = cms.bool(True),
    AndOr = cms.bool(True),
    DebugOn = cms.untracked.bool(False)
)

import Alignment.CommonAlignmentProducer.AlignmentTrackSelectorWithIndexMap_cfi
ALCARECOTkAlDstToD0Pi = Alignment.CommonAlignmentProducer.AlignmentTrackSelectorWithIndexMap_cfi.AlignmentTrackSelectorWithIndexMap.clone()
ALCARECOTkAlDstToD0Pi.filter = True
ALCARECOTkAlDstToD0Pi.src = 'generalTracks'

ALCARECOTkAlDstToD0Pi.applyBasicCuts = True
ALCARECOTkAlDstToD0Pi.ptMin = 0.5
ALCARECOTkAlDstToD0Pi.etaMin = -3.5
ALCARECOTkAlDstToD0Pi.etaMax = 3.5
ALCARECOTkAlDstToD0Pi.nHitMin = 0

ALCARECOTkAlDstToD0Pi.GlobalSelector.applyGlobalMuonFilter = False
ALCARECOTkAlDstToD0Pi.GlobalSelector.applyIsolationtest = False

ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.applyMassrangeFilter = True
ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.minXMass = 1.860
ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.maxXMass = 2.160
ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.applyIntermediateMassrangeFilter = True
ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.minIntermediateMass = 1.81483
ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.maxIntermediateMass = 1.91483
ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.applyMassDifferenceFilter = True
ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.minMassDifference = 0.14243
ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.maxMassDifference = 0.14843
ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.firstDaughterMass = 0.493677
ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.secondDaughterMass = 0.139570
ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.thirdDaughterMass = 0.139570
ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.applyChargeFilter = True
ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.charge = 1
ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.useUnsignedCharge = True
ALCARECOTkAlDstToD0Pi.ThreeBodyDecaySelector.numberOfCandidates = 1

# Persist per-track dE/dx (Harmonic2 strip + pixel-only + joint strip+pixel)
# for the selected D* daughter tracks, re-keyed onto the cloned
# ALCARECOTkAlDstToD0Pi track collection. The projection uses each cloned
# Track's preserved TrackExtraRef.key() to look up the original
# generalTracks-keyed value.
from Alignment.CommonAlignmentProducer.alcaDedxJointEstimator_cfi import alcaDedxJointEstimator
ALCARECOTkAlDstToD0PiDeDxHarmonic2 = cms.EDProducer('DeDxValueMapProjector',
    selectedTracks     = cms.InputTag('ALCARECOTkAlDstToD0Pi'),
    intermediateTracks = cms.InputTag('generalTracks'),  # selector took generalTracks directly
    sourceTracks       = cms.InputTag('generalTracks'),
    sourceValueMap     = cms.InputTag('dedxHarmonic2'),
    originalIndexMap   = cms.InputTag('ALCARECOTkAlDstToD0Pi', 'originalIndex'),
)
ALCARECOTkAlDstToD0PiDeDxPixelHarmonic2 = ALCARECOTkAlDstToD0PiDeDxHarmonic2.clone(
    sourceValueMap = cms.InputTag('dedxPixelHarmonic2'),
)
ALCARECOTkAlDstToD0PiDeDxAllHarmonic2 = ALCARECOTkAlDstToD0PiDeDxHarmonic2.clone(
    sourceValueMap = cms.InputTag('alcaDedxJointEstimator'),
)

seqALCARECOTkAlDstToD0Pi = cms.Sequence(ALCARECOTkAlDstToD0PiDCSFilter+
                                        ALCARECOTkAlDstToD0Pi+
                                        alcaDedxJointEstimator+
                                        ALCARECOTkAlDstToD0PiDeDxHarmonic2+
                                        ALCARECOTkAlDstToD0PiDeDxPixelHarmonic2+
                                        ALCARECOTkAlDstToD0PiDeDxAllHarmonic2)
