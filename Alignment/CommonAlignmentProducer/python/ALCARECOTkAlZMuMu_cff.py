import FWCore.ParameterSet.Config as cms

import HLTrigger.HLTfilters.hltHighLevel_cfi
# AlCaReco for track based alignment using ZMuMu events
ALCARECOTkAlZMuMuHLT = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    andOr = True, ## choose logical OR between Triggerbits
    eventSetupPathsKey = 'TkAlZMuMu',
    throw = False # tolerate triggers stated above, but not available
    )

# DCS partitions
# "EBp","EBm","EEp","EEm","HBHEa","HBHEb","HBHEc","HF","HO","RPC"
# "DT0","DTp","DTm","CSCp","CSCm","CASTOR","TIBTID","TOB","TECp","TECm"
# "BPIX","FPIX","ESp","ESm"
import DPGAnalysis.Skims.skim_detstatus_cfi
ALCARECOTkAlZMuMuDCSFilter = DPGAnalysis.Skims.skim_detstatus_cfi.dcsstatus.clone(
    DetectorType = cms.vstring('TIBTID','TOB','TECp','TECm','BPIX','FPIX',
                               'DT0','DTp','DTm','CSCp','CSCm'),
    ApplyFilter  = cms.bool(True),
    AndOr        = cms.bool(True),
    DebugOn      = cms.untracked.bool(False)
)

import Alignment.CommonAlignmentProducer.TkAlMuonSelectors_cfi
ALCARECOTkAlZMuMuGoodMuons = Alignment.CommonAlignmentProducer.TkAlMuonSelectors_cfi.TkAlGoodIdMuonSelector.clone()
ALCARECOTkAlZMuMuRelCombIsoMuons = Alignment.CommonAlignmentProducer.TkAlMuonSelectors_cfi.TkAlRelCombIsoMuonSelector.clone(
    src = 'ALCARECOTkAlZMuMuGoodMuons'
)

# Z candidates from all opposite-charge muon-track pairs in the mass window.
# V0-pattern pre-stage that replaces the legacy in-selector pairing.
ALCARECOTkAlZMuMuCandidates = cms.EDProducer('TwoBodyDecayCandidateProducer',
    src     = cms.InputTag('generalTracks'),
    muonSrc = cms.InputTag('ALCARECOTkAlZMuMuRelCombIsoMuons'),
    minMass        = cms.double(65.0),
    maxMass        = cms.double(115.0),
    daughterMass   = cms.double(0.105),
    daughterPdgId  = cms.int32(13),
    motherPdgId    = cms.int32(23),    ## Z
    applyChargeFilter      = cms.bool(True),
    charge                 = cms.int32(0),
    useUnsignedCharge      = cms.bool(True),
    applyAcoplanarityFilter = cms.bool(False),
    acoplanarDistance      = cms.double(1.0),
)

ALCARECOTkAlZMuMuTracks = cms.EDProducer('V0DaughterTrackProducer',
    src = cms.InputTag('ALCARECOTkAlZMuMuCandidates'),
)

import Alignment.CommonAlignmentProducer.AlignmentTrackSelectorWithIndexMap_cfi
ALCARECOTkAlZMuMu = Alignment.CommonAlignmentProducer.AlignmentTrackSelectorWithIndexMap_cfi.AlignmentTrackSelectorWithIndexMap.clone(
    src = cms.InputTag('ALCARECOTkAlZMuMuTracks'),
    filter = True, ##do not store empty events
    applyBasicCuts = True,
    ptMin  = 15.0, ##GeV
    etaMin = -3.5,
    etaMax = 3.5,
    nHitMin = 0,
)
ALCARECOTkAlZMuMu.GlobalSelector.applyGlobalMuonFilter = False
ALCARECOTkAlZMuMu.GlobalSelector.applyIsolationtest    = False

ALCARECOTkAlZMuMuResonances = cms.EDProducer('VertexCompositeCandidateRemapper',
    srcCandidates      = cms.InputTag('ALCARECOTkAlZMuMuCandidates'),
    selectedTracks     = cms.InputTag('ALCARECOTkAlZMuMu'),
    intermediateTracks = cms.InputTag('ALCARECOTkAlZMuMuTracks'),
    originalIndexMap   = cms.InputTag('ALCARECOTkAlZMuMu', 'originalIndex'),
)

seqALCARECOTkAlZMuMu = cms.Sequence(
    ALCARECOTkAlZMuMuHLT +
    ALCARECOTkAlZMuMuDCSFilter +
    ALCARECOTkAlZMuMuGoodMuons +
    ALCARECOTkAlZMuMuRelCombIsoMuons +
    ALCARECOTkAlZMuMuCandidates +
    ALCARECOTkAlZMuMuTracks +
    ALCARECOTkAlZMuMu +
    ALCARECOTkAlZMuMuResonances
)

## customizations for the pp_on_AA eras
from Configuration.Eras.Modifier_pp_on_XeXe_2017_cff import pp_on_XeXe_2017
from Configuration.Eras.Modifier_pp_on_AA_2018_cff import pp_on_AA_2018
(pp_on_XeXe_2017 | pp_on_AA_2018).toModify(ALCARECOTkAlZMuMuHLT,
                                           eventSetupPathsKey='TkAlZMuMuHI'
                                           )
