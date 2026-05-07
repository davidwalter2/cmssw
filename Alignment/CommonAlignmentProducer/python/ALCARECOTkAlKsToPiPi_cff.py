# AlCaReco for track based alignment using KS->pi+pi- events
import FWCore.ParameterSet.Config as cms

# DCS status filter -- require tracker in good operational state
import DPGAnalysis.Skims.skim_detstatus_cfi
ALCARECOTkAlKsToPiPiDCSFilter = DPGAnalysis.Skims.skim_detstatus_cfi.dcsstatus.clone(
    DetectorType = cms.vstring('TIBTID','TOB','TECp','TECm','BPIX','FPIX',
                               'DT0','DTp','DTm','CSCp','CSCm'),
    ApplyFilter  = cms.bool(True),
    AndOr        = cms.bool(True),
    DebugOn      = cms.untracked.bool(False)
)

# Local V0Producer clone with a lower track-pT cut (see ALCARECOTkAlV0Candidates_cff.py).
# Standalone clone so the standard generalV0Candidates collection (consumed by DQM,
# MINIAOD, PF, ...) is left untouched.
from Alignment.CommonAlignmentProducer.ALCARECOTkAlV0Candidates_cff import ALCARECOTkAlV0Candidates

# Event pre-filter: require >=1 reconstructed KS candidate in our private V0 collection.
# All V0Producer cuts (flight significance, pointing angle, IP significance,
# post-fit mass +/-70 MeV) are applied. No HLT filter so the ALCARECO is usable
# from any trigger path (single muon, dimuon, etc.) and is suitable for multiple
# physics channels including W/Z+KS, Psi(2S)->J/psiKS, and inclusive KS production.
ALCARECOTkAlKsToPiPiV0Filter = cms.EDFilter('CandViewCountFilter',
    src = cms.InputTag('ALCARECOTkAlV0Candidates', 'Kshort'),
    minNumber = cms.uint32(1)
)

# Extract the daughter tracks of the V0 candidates as a small TrackCollection.
# The TrackExtraRefs in each copied Track still point back to generalTracks, so
# the downstream AlignmentTrackSelectorModule + TrackCollectionStoreManager
# clones tracks + extras + hits + clusters into the ALCARECO output exactly as
# it does when reading directly from generalTracks -- but only for V0 daughters.
ALCARECOTkAlKsToPiPiV0Tracks = cms.EDProducer('V0DaughterTrackProducer',
    src = cms.InputTag('ALCARECOTkAlV0Candidates', 'Kshort'),
)

# Standard alignment track selector. Track-quality cuts mirror V0Producer's so
# every V0Producer-accepted daughter track passes; they're nominally redundant
# with the V0Producer cuts and serve mainly as documentation + a failsafe.
# TwoBodyDecaySelector is not used (left at defaults: all switches off) -- the
# V0 candidates have already passed V0Producer's tighter post-fit mass cut so
# re-pairing the tracks here would be redundant and could spuriously reject
# candidates due to small post-fit vs raw mass shifts.
import Alignment.CommonAlignmentProducer.AlignmentTrackSelectorWithIndexMap_cfi
ALCARECOTkAlKsToPiPi = Alignment.CommonAlignmentProducer.AlignmentTrackSelectorWithIndexMap_cfi.AlignmentTrackSelectorWithIndexMap.clone(
    src = cms.InputTag('ALCARECOTkAlKsToPiPiV0Tracks'),
    filter = True,
    applyBasicCuts = True,
    ptMin   = 0.1,    ## matches our local V0Producer clone tkPtCut
    etaMin  = -3.5,
    etaMax  = 3.5,
    nHitMin = 3,      ## matches V0Producer tkNHitsCut
)
ALCARECOTkAlKsToPiPi.GlobalSelector.applyGlobalMuonFilter = False
ALCARECOTkAlKsToPiPi.GlobalSelector.applyIsolationtest    = False

# Persist per-track dE/dx (Harmonic2 strip + pixel-only + joint strip+pixel)
# for the selected V0 daughters, re-keyed onto the cloned
# ALCARECOTkAlKsToPiPi track collection. The projection uses each cloned
# Track's preserved TrackExtraRef.key() to look up the original
# generalTracks-keyed value.
from Alignment.CommonAlignmentProducer.alcaDedxJointEstimator_cfi import alcaDedxJointEstimator
ALCARECOTkAlKsToPiPiDeDxHarmonic2 = cms.EDProducer('DeDxValueMapProjector',
    selectedTracks     = cms.InputTag('ALCARECOTkAlKsToPiPi'),
    intermediateTracks = cms.InputTag('ALCARECOTkAlKsToPiPiV0Tracks'),
    sourceTracks       = cms.InputTag('generalTracks'),
    sourceValueMap     = cms.InputTag('dedxHarmonic2'),
    originalIndexMap   = cms.InputTag('ALCARECOTkAlKsToPiPi', 'originalIndex'),
)
ALCARECOTkAlKsToPiPiDeDxPixelHarmonic2 = ALCARECOTkAlKsToPiPiDeDxHarmonic2.clone(
    sourceValueMap = cms.InputTag('dedxPixelHarmonic2'),
)
ALCARECOTkAlKsToPiPiDeDxAllHarmonic2 = ALCARECOTkAlKsToPiPiDeDxHarmonic2.clone(
    sourceValueMap = cms.InputTag('alcaDedxJointEstimator'),
)

# Re-key the V0 candidate collection's daughter TrackRefs onto the cloned
# AlignmentTrackSelector output so downstream consumers can navigate
# candidate -> daughter -> track without dereferencing generalTracks.
# Candidates whose daughters were dropped by AlignmentTrackSelector are
# silently removed. Drives off the AlignmentTrackSelectorWithIndexMapModule
# side-channel ValueMap of source-track indices.
ALCARECOTkAlKsToPiPiResonances = cms.EDProducer('VertexCompositeCandidateRemapper',
    srcCandidates      = cms.InputTag('ALCARECOTkAlV0Candidates', 'Kshort'),
    selectedTracks     = cms.InputTag('ALCARECOTkAlKsToPiPi'),
    intermediateTracks = cms.InputTag('ALCARECOTkAlKsToPiPiV0Tracks'),
    originalIndexMap   = cms.InputTag('ALCARECOTkAlKsToPiPi', 'originalIndex'),
)

seqALCARECOTkAlKsToPiPi = cms.Sequence(
    ALCARECOTkAlKsToPiPiDCSFilter +
    ALCARECOTkAlV0Candidates +
    ALCARECOTkAlKsToPiPiV0Filter +
    ALCARECOTkAlKsToPiPiV0Tracks +
    ALCARECOTkAlKsToPiPi +
    ALCARECOTkAlKsToPiPiResonances +
    alcaDedxJointEstimator +
    ALCARECOTkAlKsToPiPiDeDxHarmonic2 +
    ALCARECOTkAlKsToPiPiDeDxPixelHarmonic2 +
    ALCARECOTkAlKsToPiPiDeDxAllHarmonic2
)
