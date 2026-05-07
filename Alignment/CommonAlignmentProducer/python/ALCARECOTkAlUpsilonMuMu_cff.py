# AlCaReco for track based alignment using Upsilon->MuMu events
import FWCore.ParameterSet.Config as cms

import HLTrigger.HLTfilters.hltHighLevel_cfi
ALCARECOTkAlUpsilonMuMuHLT = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    andOr = True, ## choose logical OR between Triggerbits
    eventSetupPathsKey = 'TkAlUpsilonMuMu',
    throw = False # tolerate triggers stated above, but not available
    )

# DCS partitions
# "EBp","EBm","EEp","EEm","HBHEa","HBHEb","HBHEc","HF","HO","RPC"
# "DT0","DTp","DTm","CSCp","CSCm","CASTOR","TIBTID","TOB","TECp","TECm"
# "BPIX","FPIX","ESp","ESm"
import DPGAnalysis.Skims.skim_detstatus_cfi
ALCARECOTkAlUpsilonMuMuDCSFilter = DPGAnalysis.Skims.skim_detstatus_cfi.dcsstatus.clone(
    DetectorType = cms.vstring('TIBTID','TOB','TECp','TECm','BPIX','FPIX',
                               'DT0','DTp','DTm','CSCp','CSCm'),
    ApplyFilter  = cms.bool(True),
    AndOr        = cms.bool(True),
    DebugOn      = cms.untracked.bool(False)
)

import Alignment.CommonAlignmentProducer.TkAlMuonSelectors_cfi
ALCARECOTkAlUpsilonMuMuGoodMuons = Alignment.CommonAlignmentProducer.TkAlMuonSelectors_cfi.TkAlGoodIdMuonSelector.clone()
ALCARECOTkAlUpsilonMuMuRelCombIsoMuons = Alignment.CommonAlignmentProducer.TkAlMuonSelectors_cfi.TkAlRelCombIsoMuonSelector.clone(
    src = 'ALCARECOTkAlUpsilonMuMuGoodMuons',
    cut = '(isolationR03().sumPt + isolationR03().emEt + isolationR03().hadEt)/pt  < 0.3'

)

# Upsilon candidates from all opposite-charge muon-track pairs in the mass
# window. V0-pattern pre-stage that replaces the legacy in-selector pairing.
ALCARECOTkAlUpsilonMuMuCandidates = cms.EDProducer('TwoBodyDecayCandidateProducer',
    src     = cms.InputTag('generalTracks'),
    muonSrc = cms.InputTag('ALCARECOTkAlUpsilonMuMuRelCombIsoMuons'),
    minMass        = cms.double(8.9),
    maxMass        = cms.double(9.9),
    daughterMass   = cms.double(0.105),
    daughterPdgId  = cms.int32(13),
    motherPdgId    = cms.int32(553),   ## Upsilon(1S)
    applyChargeFilter      = cms.bool(True),
    charge                 = cms.int32(0),
    useUnsignedCharge      = cms.bool(True),
    applyAcoplanarityFilter = cms.bool(False),
    acoplanarDistance      = cms.double(1.0),
)

ALCARECOTkAlUpsilonMuMuTracks = cms.EDProducer('V0DaughterTrackProducer',
    src = cms.InputTag('ALCARECOTkAlUpsilonMuMuCandidates'),
)

import Alignment.CommonAlignmentProducer.AlignmentTrackSelectorWithIndexMap_cfi
ALCARECOTkAlUpsilonMuMu = Alignment.CommonAlignmentProducer.AlignmentTrackSelectorWithIndexMap_cfi.AlignmentTrackSelectorWithIndexMap.clone(
    src = cms.InputTag('ALCARECOTkAlUpsilonMuMuTracks'),
    filter = True, ##do not store empty events
    applyBasicCuts = True,
    ptMin  = 3.,  ##GeV
    etaMin = -3.5,
    etaMax = 3.5,
    nHitMin = 0,
)
ALCARECOTkAlUpsilonMuMu.GlobalSelector.applyGlobalMuonFilter = False
ALCARECOTkAlUpsilonMuMu.GlobalSelector.applyIsolationtest    = False

ALCARECOTkAlUpsilonMuMuResonances = cms.EDProducer('VertexCompositeCandidateRemapper',
    srcCandidates      = cms.InputTag('ALCARECOTkAlUpsilonMuMuCandidates'),
    selectedTracks     = cms.InputTag('ALCARECOTkAlUpsilonMuMu'),
    intermediateTracks = cms.InputTag('ALCARECOTkAlUpsilonMuMuTracks'),
    originalIndexMap   = cms.InputTag('ALCARECOTkAlUpsilonMuMu', 'originalIndex'),
)

seqALCARECOTkAlUpsilonMuMu = cms.Sequence(
    ALCARECOTkAlUpsilonMuMuHLT +
    ALCARECOTkAlUpsilonMuMuDCSFilter +
    ALCARECOTkAlUpsilonMuMuGoodMuons +
    ALCARECOTkAlUpsilonMuMuRelCombIsoMuons +
    ALCARECOTkAlUpsilonMuMuCandidates +
    ALCARECOTkAlUpsilonMuMuTracks +
    ALCARECOTkAlUpsilonMuMu +
    ALCARECOTkAlUpsilonMuMuResonances
)

## customizations for the pp_on_AA eras
from Configuration.Eras.Modifier_pp_on_XeXe_2017_cff import pp_on_XeXe_2017
from Configuration.Eras.Modifier_pp_on_AA_2018_cff import pp_on_AA_2018
(pp_on_XeXe_2017 | pp_on_AA_2018).toModify(ALCARECOTkAlUpsilonMuMuHLT,
                                           eventSetupPathsKey='TkAlUpsilonMuMuHI'
                                           )
