# AlCaReco for track based alignment using J/Psi->MuMu events
import FWCore.ParameterSet.Config as cms

import HLTrigger.HLTfilters.hltHighLevel_cfi
ALCARECOTkAlJpsiMuMuHLT = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    andOr = True, ## choose logical OR between Triggerbits
    eventSetupPathsKey = 'TkAlJpsiMuMu',
    throw = False # tolerate triggers stated above, but not available
    )

# DCS partitions
# "EBp","EBm","EEp","EEm","HBHEa","HBHEb","HBHEc","HF","HO","RPC"
# "DT0","DTp","DTm","CSCp","CSCm","CASTOR","TIBTID","TOB","TECp","TECm"
# "BPIX","FPIX","ESp","ESm"
import DPGAnalysis.Skims.skim_detstatus_cfi
ALCARECOTkAlJpsiMuMuDCSFilter = DPGAnalysis.Skims.skim_detstatus_cfi.dcsstatus.clone(
    DetectorType = cms.vstring('TIBTID','TOB','TECp','TECm','BPIX','FPIX',
                               'DT0','DTp','DTm','CSCp','CSCm'),
    ApplyFilter  = cms.bool(True),
    AndOr        = cms.bool(True),
    DebugOn      = cms.untracked.bool(False)
)

import Alignment.CommonAlignmentProducer.TkAlMuonSelectors_cfi
ALCARECOTkAlJpsiMuMuGoodMuons = Alignment.CommonAlignmentProducer.TkAlMuonSelectors_cfi.TkAlGoodIdMuonSelector.clone()

# Build the J/psi candidate collection upstream of AlignmentTrackSelector so
# the candidate object survives downstream (V0-pattern, mirrors KsToPiPi).
# All passing pairs are emitted (not the one-best of the legacy
# AlignmentTwoBodyDecayTrackSelector); per-pair charge & mass cuts here mirror
# the legacy TwoBodyDecaySelector configuration. Muon-id filter (replacing
# applyGlobalMuonFilter) is enforced via muonSrc.
ALCARECOTkAlJpsiMuMuCandidates = cms.EDProducer('TwoBodyDecayCandidateProducer',
    src     = cms.InputTag('generalTracks'),
    muonSrc = cms.InputTag('ALCARECOTkAlJpsiMuMuGoodMuons'),
    minMass        = cms.double(2.7),  ## GeV
    maxMass        = cms.double(3.4),  ## GeV
    daughterMass   = cms.double(0.105),
    daughterPdgId  = cms.int32(13),    ## mu-
    motherPdgId    = cms.int32(443),   ## J/psi
    applyChargeFilter      = cms.bool(False),
    charge                 = cms.int32(0),
    useUnsignedCharge      = cms.bool(True),
    applyAcoplanarityFilter = cms.bool(False),
    acoplanarDistance      = cms.double(1.0),
)

# Extract the daughter tracks of the candidates as a small TrackCollection.
# TrackExtraRefs in each copied Track still point back to generalTracks, so
# the downstream AlignmentTrackSelectorModule + TrackCollectionStoreManager
# clones tracks + extras + hits + clusters into the ALCAREco output for the
# J/psi daughter tracks only.
ALCARECOTkAlJpsiMuMuTracks = cms.EDProducer('V0DaughterTrackProducer',
    src = cms.InputTag('ALCARECOTkAlJpsiMuMuCandidates'),
)

import Alignment.CommonAlignmentProducer.AlignmentTrackSelectorWithIndexMap_cfi
# Drop-in replacement for AlignmentTrackSelectorModule: same cloned outputs
# (tracks/extras/hits/clusters) plus a ValueMap<unsigned int>
# (instance label "originalIndex") of source-track indices, computed by
# pointer arithmetic on the selector chain's Track* output. Downstream
# remapping (candidates, dE/dx) consumes this index map directly so no
# kinematic fingerprinting is required.
ALCARECOTkAlJpsiMuMu = Alignment.CommonAlignmentProducer.AlignmentTrackSelectorWithIndexMap_cfi.AlignmentTrackSelectorWithIndexMap.clone(
    src = cms.InputTag('ALCARECOTkAlJpsiMuMuTracks'),
    filter = True, ##do not store empty events
    applyBasicCuts = True,
    ptMin  = 0.8, ##GeV
    etaMin = -3.5,
    etaMax = 3.5,
    nHitMin = 0,
)
# Muon-id and pair-finding moved upstream into the candidate producer; here
# we only apply per-track quality cuts.
ALCARECOTkAlJpsiMuMu.GlobalSelector.applyGlobalMuonFilter = False
ALCARECOTkAlJpsiMuMu.GlobalSelector.applyIsolationtest    = False

# Re-key the candidate collection's daughter TrackRefs onto the cloned
# AlignmentTrackSelector output so downstream consumers can navigate
# candidate -> daughter -> track without dereferencing generalTracks.
# Candidates whose daughters were dropped by AlignmentTrackSelector are
# silently removed.
ALCARECOTkAlJpsiMuMuResonances = cms.EDProducer('VertexCompositeCandidateRemapper',
    srcCandidates      = cms.InputTag('ALCARECOTkAlJpsiMuMuCandidates'),
    selectedTracks     = cms.InputTag('ALCARECOTkAlJpsiMuMu'),
    intermediateTracks = cms.InputTag('ALCARECOTkAlJpsiMuMuTracks'),
    originalIndexMap   = cms.InputTag('ALCARECOTkAlJpsiMuMu', 'originalIndex'),
)

seqALCARECOTkAlJpsiMuMu = cms.Sequence(
    ALCARECOTkAlJpsiMuMuHLT +
    ALCARECOTkAlJpsiMuMuDCSFilter +
    ALCARECOTkAlJpsiMuMuGoodMuons +
    ALCARECOTkAlJpsiMuMuCandidates +
    ALCARECOTkAlJpsiMuMuTracks +
    ALCARECOTkAlJpsiMuMu +
    ALCARECOTkAlJpsiMuMuResonances
)

## customizations for the pp_on_AA eras
from Configuration.Eras.Modifier_pp_on_XeXe_2017_cff import pp_on_XeXe_2017
from Configuration.Eras.Modifier_pp_on_AA_2018_cff import pp_on_AA_2018
(pp_on_XeXe_2017 | pp_on_AA_2018).toModify(ALCARECOTkAlJpsiMuMuHLT,
                                           eventSetupPathsKey='TkAlJpsiMuMuHI'
                                           )
