# AlCaReco for track based alignment using Lambda0->proton+pi- events
import FWCore.ParameterSet.Config as cms

# DCS status filter -- require tracker in good operational state
import DPGAnalysis.Skims.skim_detstatus_cfi
ALCARECOTkAlLambdaToProtonPiDCSFilter = DPGAnalysis.Skims.skim_detstatus_cfi.dcsstatus.clone(
    DetectorType = cms.vstring('TIBTID','TOB','TECp','TECm','BPIX','FPIX',
                               'DT0','DTp','DTm','CSCp','CSCm'),
    ApplyFilter  = cms.bool(True),
    AndOr        = cms.bool(True),
    DebugOn      = cms.untracked.bool(False)
)

# Local V0Producer clone with a lower track-pT cut (see ALCARECOTkAlV0Candidates_cff.py).
# Shared between the KS and Lambda ALCARECOs (deduplicated by the framework so it
# runs at most once per event even if both ALCARECO paths reference it).
from Alignment.CommonAlignmentProducer.ALCARECOTkAlV0Candidates_cff import ALCARECOTkAlV0Candidates

# Event pre-filter: require >=1 reconstructed Lambda0 candidate in our private
# V0 collection. All V0Producer cuts (incl. post-fit mass +/-50 MeV) are applied.
ALCARECOTkAlLambdaToProtonPiV0Filter = cms.EDFilter('CandViewCountFilter',
    src = cms.InputTag('ALCARECOTkAlV0Candidates', 'Lambda'),
    minNumber = cms.uint32(1)
)

# Extract the daughter tracks of the V0 candidates as a small TrackCollection.
# By V0Producer convention, daughter(0) = baryon (p or pbar) and daughter(1) =
# pion (pi- or pi+); both Lambda0 and anti-Lambda0 are accepted automatically.
# The TrackExtraRefs in each copied Track still point back to generalTracks, so
# the downstream AlignmentTrackSelectorModule + TrackCollectionStoreManager
# clones tracks + extras + hits + clusters into the ALCARECO output.
ALCARECOTkAlLambdaToProtonPiV0Tracks = cms.EDProducer('V0DaughterTrackProducer',
    src = cms.InputTag('ALCARECOTkAlV0Candidates', 'Lambda'),
)

# Standard alignment track selector. Track-quality cuts mirror our V0Producer
# clone's so every accepted daughter track passes; they're nominally redundant
# and serve mainly as documentation + a failsafe. TwoBodyDecaySelector is not
# used (left at defaults: all switches off) -- the V0 candidates have already
# passed V0Producer's tighter post-fit mass cut so re-pairing the tracks here
# would be redundant and could spuriously reject candidates due to small
# post-fit vs raw mass shifts.
import Alignment.CommonAlignmentProducer.AlignmentTrackSelectorWithIndexMap_cfi
ALCARECOTkAlLambdaToProtonPi = Alignment.CommonAlignmentProducer.AlignmentTrackSelectorWithIndexMap_cfi.AlignmentTrackSelectorWithIndexMap.clone(
    src = cms.InputTag('ALCARECOTkAlLambdaToProtonPiV0Tracks'),
    filter = True,
    applyBasicCuts = True,
    ptMin   = 0.1,    ## matches our local V0Producer clone tkPtCut
    etaMin  = -3.5,
    etaMax  = 3.5,
    nHitMin = 3,      ## matches V0Producer tkNHitsCut
)
ALCARECOTkAlLambdaToProtonPi.GlobalSelector.applyGlobalMuonFilter = False
ALCARECOTkAlLambdaToProtonPi.GlobalSelector.applyIsolationtest    = False

# Persist per-track dE/dx (Harmonic2 strip + pixel-only + joint strip+pixel)
# for the selected V0 daughters, re-keyed onto the cloned
# ALCARECOTkAlLambdaToProtonPi track collection. The projection uses each
# cloned Track's preserved TrackExtraRef.key() to look up the original
# generalTracks-keyed value.
from Alignment.CommonAlignmentProducer.alcaDedxJointEstimator_cfi import alcaDedxJointEstimator
ALCARECOTkAlLambdaToProtonPiDeDxHarmonic2 = cms.EDProducer('DeDxValueMapProjector',
    selectedTracks     = cms.InputTag('ALCARECOTkAlLambdaToProtonPi'),
    intermediateTracks = cms.InputTag('ALCARECOTkAlLambdaToProtonPiV0Tracks'),
    sourceTracks       = cms.InputTag('generalTracks'),
    sourceValueMap     = cms.InputTag('dedxHarmonic2'),
    originalIndexMap   = cms.InputTag('ALCARECOTkAlLambdaToProtonPi', 'originalIndex'),
)
ALCARECOTkAlLambdaToProtonPiDeDxPixelHarmonic2 = ALCARECOTkAlLambdaToProtonPiDeDxHarmonic2.clone(
    sourceValueMap = cms.InputTag('dedxPixelHarmonic2'),
)
ALCARECOTkAlLambdaToProtonPiDeDxAllHarmonic2 = ALCARECOTkAlLambdaToProtonPiDeDxHarmonic2.clone(
    sourceValueMap = cms.InputTag('alcaDedxJointEstimator'),
)

# Re-key the V0 candidate collection's daughter TrackRefs onto the cloned
# AlignmentTrackSelector output so downstream consumers can navigate
# candidate -> daughter -> track without dereferencing generalTracks.
# Candidates whose daughters were dropped by AlignmentTrackSelector are
# silently removed. Drives off the AlignmentTrackSelectorWithIndexMapModule
# side-channel ValueMap of source-track indices.
ALCARECOTkAlLambdaToProtonPiResonances = cms.EDProducer('VertexCompositeCandidateRemapper',
    srcCandidates      = cms.InputTag('ALCARECOTkAlV0Candidates', 'Lambda'),
    selectedTracks     = cms.InputTag('ALCARECOTkAlLambdaToProtonPi'),
    intermediateTracks = cms.InputTag('ALCARECOTkAlLambdaToProtonPiV0Tracks'),
    originalIndexMap   = cms.InputTag('ALCARECOTkAlLambdaToProtonPi', 'originalIndex'),
)

seqALCARECOTkAlLambdaToProtonPi = cms.Sequence(
    ALCARECOTkAlLambdaToProtonPiDCSFilter +
    ALCARECOTkAlV0Candidates +
    ALCARECOTkAlLambdaToProtonPiV0Filter +
    ALCARECOTkAlLambdaToProtonPiV0Tracks +
    ALCARECOTkAlLambdaToProtonPi +
    ALCARECOTkAlLambdaToProtonPiResonances +
    alcaDedxJointEstimator +
    ALCARECOTkAlLambdaToProtonPiDeDxHarmonic2 +
    ALCARECOTkAlLambdaToProtonPiDeDxPixelHarmonic2 +
    ALCARECOTkAlLambdaToProtonPiDeDxAllHarmonic2
)
