# AlCaReco for track-based alignment using J/psi + X events
# (B+->J/psiK, B0->J/psiK*0, B0->J/psiKs, Bs->J/psiphi, Lb->J/psiLambda,
#  psi(2S)->J/psiKs, Bc->J/psipi)
#
# Five-stage single-sequence architecture:
#   Stage 0: DCS filter + good-muon selector. No HLT filter in the default
#            sequence; ALCARECOTkAlJpsiXHLT is defined below but excluded
#            from seqALCARECOTkAlJpsiX. Prepend it for online/express use.
#   Stage 1: Shared candidate construction (J/psi, V0s, K*0, phi).
#   Stage 2: Per-channel B-meson / quarkonium candidates (7 instances of
#            JpsiXCandidateProducer).
#   Stage 3: ALCARECOTkAlJpsiXAllTracks -- CompositeDaughterTrackProducer
#            that merges and deduplicates all leaf tracks from all 7 VCC
#            collections into one TrackCollection.
#   Stage 4: ALCARECOTkAlJpsiX -- AlignmentTrackSelectorWithIndexMapModule
#            on the merged collection (clones tracks + extras + hits +
#            clusters, emits originalIndex ValueMap) + 3 DeDxValueMapProjectors.
#   Stage 5: 7 x VertexCompositeCandidateRemapper re-keying daughter
#            TrackRefs onto the cloned ALCARECOTkAlJpsiX collection.
#
# TrackRef chain (correctness argument):
#   JpsiXCandidateProducer stores TrackRef(generalTracks, K) for every leaf.
#   CompositeDaughterTrackProducer shallow-copies each track; extra().key()
#   is still K.  AlignmentTrackSelectorWithIndexMapModule deep-clones tracks
#   and emits originalIndex[j] = origIdx (cloned->intermediate).
#   VertexCompositeCandidateRemapper builds origToSel[K] = j and rewrites
#   daughter.track().key() K -> j in the cloned collection.

import os
import FWCore.ParameterSet.Config as cms

# ---------------------------------------------------------------------------
# Selection preset switching (TkAlJpsiX -- non-V0 channels only)
# ---------------------------------------------------------------------------
# Set the env var TKALJPSIX_SELECTION_PRESET before cmsDriver/cmsRun:
#   A = mass windows only (channel-defining cuts, no kinematic refinement).
#   B = mass windows + kinematic + geometric, low-pT bachelor / daughter floor
#       (Phase-1-locked values; default; bachelor pT > 0.1 GeV matches the
#       central V0 clone's tkPtCut).
#   C = preset B + B-level Kalman vertex fit on the 4 non-V0 channels
#       (minBVtxProb, maxMotherAlphaBS, minBLxyOverSigma activated). 3-body
#       Kalman fit in track mode (B+, Bc); 4-body Kalman fit on all leaf tracks
#       in VCC mode (B0->K*0, Bs->phi).
# V0-mode channels (B0->Ks, Lambda_b, psi(2S)) are preset-invariant (they ride
# on the central V0Producer which already gives clean candidates at ~1/event).
_TKALJPSIX_SELECTION_PRESET = os.environ.get('TKALJPSIX_SELECTION_PRESET', 'B')
if _TKALJPSIX_SELECTION_PRESET not in ('A', 'B', 'C'):
    raise ValueError(
        "TKALJPSIX_SELECTION_PRESET must be 'A', 'B', or 'C'; got "
        + repr(_TKALJPSIX_SELECTION_PRESET))

# Sentinel that disables a cut (used when preset A wants no cut on a field
# that preset B uses).
_DISABLED = 1e9

# Preset-C Kalman cut values shared across BPlus / Bc / B0Kstar / BsPhi.
# alpha_BS is interpreted as a max-angle (radians) cut; the C++ implementation
# computes acos(cos_alpha) in 3D from the Kalman-fitted vertex relative to the
# beamspot. Lxy/sigma cut uses the first offlinePrimaryVertices entry.
#
# **Iteration 1 (2026-06-11)**: Phase-2 v1 with the proposal's initial values
# (vtxProb>0.01, cos(alpha_xy)>0.99 -> 0.142 rad, Lxy/sigma>3) was 10x more
# aggressive than Phase-1 closed-form vertex predicted -- preset C killed all
# non-V0 channels (Bs->phi went to zero cands). Relaxing two of three cuts in
# one iteration:
#   minBLxyOverSigma  3.0 -> 1.0   (probable dominant kill; Bc-lifetime aware)
#   maxMotherAlphaBS  acos(0.99) ~ 0.142 rad -> acos(0.95) ~ 0.317 rad
#                                   (compensates 3D vs xy convention)
#   minBVtxProb       0.01 unchanged (standard B-physics cut)
import math as _math
_PRESETC_KALMAN = dict(
    minBVtxProb       = 0.01,
    maxMotherAlphaBS  = _math.acos(0.95),    # ~ 0.3176 rad (was acos(0.99) ~ 0.1415)
    minBLxyOverSigma  = 1.0,                 # was 3.0
)

# Per-channel non-V0 kinematic + geometric cuts. Mass windows are NOT keyed by
# preset -- they are channel-defining. Phase-1 lock: A vs B on B+ MC settled
# preset B at (Muon_pt>4, |kaon_eta|<2.5, kaon-mu-DOCA<0.03 cm, raw J/psi pT > 3,
# raw B pT > 5). Bachelor / daughter pT lowered from 1.5 / 1.0 GeV to 0.1 GeV
# to align with the V0 clone's tkPtCut and serve the alignment use case
# (openspec change add-jpsi-x-vertex-fit-and-low-pt).
# Bachelor / daughter eta cap raised to the tracker acceptance edge (2.5)
# under every preset (openspec change add-jpsi-x-muons-and-preprod-refinements).
_NON_V0_PRESETS = {
    'A': {
        'BPlus':   dict(minBachelorPt=0.5,  maxBachelorEta=2.5, minJpsiPt=0., minMotherPt=0.,
                        maxBachelorMuTrackDOCA=_DISABLED),
        'Bc':      dict(minBachelorPt=0.3,  maxBachelorEta=2.5, minJpsiPt=0., minMotherPt=0.,
                        maxBachelorMuTrackDOCA=_DISABLED),
        'B0Kstar': dict(minJpsiPt=0., minMotherPt=0., maxBachelorMuTrackDOCA=_DISABLED),
        'BsPhi':   dict(minJpsiPt=0., minMotherPt=0., maxBachelorMuTrackDOCA=_DISABLED),
        # K*0 / phi sub-resonance daughter cuts (TwoBodyDecayCandidateProducer)
        'Kstar':   dict(minDaughterPt=0.0, maxDaughterEta=2.5),
        'Phi':     dict(minDaughterPt=0.0, maxDaughterEta=2.5),
    },
    'B': {
        'BPlus':   dict(minBachelorPt=0.1,  maxBachelorEta=2.5, minJpsiPt=3., minMotherPt=5.,
                        maxBachelorMuTrackDOCA=0.03),
        'Bc':      dict(minBachelorPt=0.1,  maxBachelorEta=2.5, minJpsiPt=3., minMotherPt=5.,
                        maxBachelorMuTrackDOCA=0.03),
        'B0Kstar': dict(minJpsiPt=3., minMotherPt=5., maxBachelorMuTrackDOCA=0.03),
        'BsPhi':   dict(minJpsiPt=3., minMotherPt=5., maxBachelorMuTrackDOCA=0.03),
        # K*0 / phi sub-resonance daughter cuts
        'Kstar':   dict(minDaughterPt=0.1, maxDaughterEta=2.5),
        'Phi':     dict(minDaughterPt=0.1, maxDaughterEta=2.5),
    },
    'C': {
        # Inherits updated preset B kinematics; adds Kalman vertex-fit cuts.
        'BPlus':   dict(minBachelorPt=0.1,  maxBachelorEta=2.5, minJpsiPt=3., minMotherPt=5.,
                        maxBachelorMuTrackDOCA=0.03, **_PRESETC_KALMAN),
        'Bc':      dict(minBachelorPt=0.1,  maxBachelorEta=2.5, minJpsiPt=3., minMotherPt=5.,
                        maxBachelorMuTrackDOCA=0.03, **_PRESETC_KALMAN),
        'B0Kstar': dict(minJpsiPt=3., minMotherPt=5., maxBachelorMuTrackDOCA=0.03,
                        **_PRESETC_KALMAN),
        'BsPhi':   dict(minJpsiPt=3., minMotherPt=5., maxBachelorMuTrackDOCA=0.03,
                        **_PRESETC_KALMAN),
        # K*0 / phi sub-resonance daughter cuts (same as B; no sub-resonance
        # Kalman fit added here -- the B-level multi-body fit in
        # JpsiXCandidateProducer already constrains all leaf tracks to a common
        # B vertex).
        'Kstar':   dict(minDaughterPt=0.1, maxDaughterEta=2.5),
        'Phi':     dict(minDaughterPt=0.1, maxDaughterEta=2.5),
    },
}[_TKALJPSIX_SELECTION_PRESET]


# ---------------------------------------------------------------------------
# Stage 0 -- DCS filter + good-muon selector
# ---------------------------------------------------------------------------

# HLT filter: defined but NOT included in seqALCARECOTkAlJpsiX.
# Targets reprocessing / calibration use cases where all J/psi candidates
# are wanted regardless of trigger. Prepend for online/express production.
import HLTrigger.HLTfilters.hltHighLevel_cfi
ALCARECOTkAlJpsiXHLT = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    andOr               = True,
    eventSetupPathsKey  = 'TkAlJpsiMuMu',   # reuse existing key; no HLT clone needed
    throw               = False,
)

import DPGAnalysis.Skims.skim_detstatus_cfi
ALCARECOTkAlJpsiXDCSFilter = DPGAnalysis.Skims.skim_detstatus_cfi.dcsstatus.clone(
    DetectorType = cms.vstring('TIBTID','TOB','TECp','TECm','BPIX','FPIX',
                               'DT0','DTp','DTm','CSCp','CSCm'),
    ApplyFilter  = cms.bool(True),
    AndOr        = cms.bool(True),
    DebugOn      = cms.untracked.bool(False),
)

import Alignment.CommonAlignmentProducer.TkAlMuonSelectors_cfi
ALCARECOTkAlJpsiXGoodMuons = Alignment.CommonAlignmentProducer.TkAlMuonSelectors_cfi.TkAlGoodIdMuonSelector.clone()

# Looser muon selector for the J/psi-only production channel below.
# Superset of ALCARECOTkAlJpsiXGoodMuons by construction (relaxes the
# tracker&global AND to OR, drops the two globalTrack.* sub-cuts). Also
# doubles as the persisted-muon source (ALCARECOTkAlJpsiXMuons) since a
# superset of the tight selector's output.
ALCARECOTkAlJpsiXLooseMuons = Alignment.CommonAlignmentProducer.TkAlMuonSelectors_cfi.TkAlLooseIdMuonSelector.clone(
    filter = cms.bool(False)   ## don't gate the event on loose-muon presence;
                               ## the tight selector gates it, which is a
                               ## subset so loose is already non-empty when
                               ## tight is
)

# ---------------------------------------------------------------------------
# Stage 1 -- Shared candidate construction
# ---------------------------------------------------------------------------

# J/psi -> mu+mu- candidates (muon filter applied via muonSrc).
ALCARECOTkAlJpsiXJpsiCandidates = cms.EDProducer('TwoBodyDecayCandidateProducer',
    src           = cms.InputTag('generalTracks'),
    muonSrc       = cms.InputTag('ALCARECOTkAlJpsiXGoodMuons'),
    minMass       = cms.double(2.95),   ## GeV  +/-5sigma at 30 MeV resolution
    maxMass       = cms.double(3.25),   ## GeV
    daughterMass  = cms.double(0.105),  ## muon
    daughterPdgId = cms.int32(13),      ## mu-
    motherPdgId   = cms.int32(443),     ## J/psi
    applyChargeFilter        = cms.bool(True),
    charge                   = cms.int32(0),
    useUnsignedCharge        = cms.bool(True),
    applyAcoplanarityFilter  = cms.bool(False),
    acoplanarDistance        = cms.double(1.0),
)

# Ks and Lambda V0 candidates: shared with TkAlKsToPiPi and TkAlLambdaToProtonPi
# via CMS framework deduplication (module runs once regardless of how many
# sequences reference it in the same job).
from Alignment.CommonAlignmentProducer.ALCARECOTkAlV0Candidates_cff import ALCARECOTkAlV0Candidates

# K*0(892) -> K+pi- candidates: asymmetric mass mode with both charge
# assignments tried per pair so both K+pi- and K-pi+ are emitted when
# both fall in the mass window.
ALCARECOTkAlJpsiXKstarCandidates = cms.EDProducer('TwoBodyDecayCandidateProducer',
    src                     = cms.InputTag('generalTracks'),
    muonSrc                 = cms.InputTag(''),       ## no muon filter
    daughterMass            = cms.double(0.493677),   ## required base parameter; overridden by first/secondDaughterMass
    daughterPdgId           = cms.int32(321),         ## required base parameter; overridden by first/secondDaughterPdgId
    firstDaughterMass       = cms.double(0.493677),   ## K
    secondDaughterMass      = cms.double(0.139570),   ## pi
    firstDaughterPdgId      = cms.int32(321),         ## K+
    secondDaughterPdgId     = cms.int32(211),         ## pi+
    motherPdgId             = cms.int32(313),         ## K*0
    minMass                 = cms.double(0.80),   ## +/-2 Gamma around K*0(892)
    maxMass                 = cms.double(0.99),
    minDaughterPt           = cms.double(_NON_V0_PRESETS['Kstar']['minDaughterPt']),
    maxDaughterEta          = cms.double(_NON_V0_PRESETS['Kstar']['maxDaughterEta']),
    tryBothChargeAssignments = cms.bool(True),
    applyChargeFilter       = cms.bool(True),
    charge                  = cms.int32(0),
    useUnsignedCharge       = cms.bool(True),
    applyAcoplanarityFilter = cms.bool(False),
    acoplanarDistance       = cms.double(1.0),
    applyVertexFit          = cms.bool(False),   ## no Kalman fit at AlCaReco (cost)
    minVtxProb              = cms.double(0.0),
)

# Dipion pair for psi(2S) -> J/psi pi+pi- (BR ~34.7 %). This replaces
# the previous psi(2S) wiring which sourced pi+pi- from the V0 Ks
# output -- that was a physics bug (psi(2S) does not have a Ks-quality
# displaced pi+pi- vertex; the dipion is at the psi(2S) decay vertex).
# psi(2S) itself can be prompt (~80 % of the total psi(2S) rate at CMS)
# or non-prompt from B decays (~20 %); this producer accepts both since
# preset B applies no displacement cut. See openspec change
# add-jpsi-x-muons-and-preprod-refinements.
#
# Mass window (0.28, 0.65) covers the physical dipion mass range: 2*m_pi
# threshold at 0.279 GeV up to (m_psi(2S) - m_J/psi) ~ 0.589 GeV. The
# spectrum peaks near 0.5 GeV due to chiral dynamics; window has small
# sideband margin above the phase-space edge.
#
# Combinatorial control via maxTrackTrackDOCA (openspec item 6):
# static straight-line 3D DCA between the two pion tracks. Real pi+pi-
# from a common psi(2S) vertex have DCAs at the beamspot + tracking
# resolution scale (few tens of um), well below the 0.03 cm cut.
# Random-track combinatorics have DCAs distributed up to cm-scale and
# are rejected. Matches the physical scale of maxBachelorMuTrackDOCA
# used on the four non-V0 preset-B JpsiXCandidateProducer instances.
# Not a Kalman fit -- purely geometric.
ALCARECOTkAlJpsiXPiPiCandidates = cms.EDProducer('TwoBodyDecayCandidateProducer',
    src           = cms.InputTag('generalTracks'),
    muonSrc       = cms.InputTag(''),           ## no muon filter
    daughterMass  = cms.double(0.139570),       ## pion
    daughterPdgId = cms.int32(211),             ## pi+
    motherPdgId   = cms.int32(100443),          ## psi(2S) (candidate tag only)
    minMass       = cms.double(0.28),
    maxMass       = cms.double(0.65),
    minDaughterPt            = cms.double(0.1), ## matches K*0 / phi under all presets
    maxDaughterEta           = cms.double(2.5), ## matches K*0 / phi under all presets
    applyChargeFilter        = cms.bool(True),
    charge                   = cms.int32(0),
    useUnsignedCharge        = cms.bool(True),
    applyAcoplanarityFilter  = cms.bool(False),
    acoplanarDistance        = cms.double(1.0),
    applyVertexFit           = cms.bool(False), ## no Kalman fit at AlCaReco
    minVtxProb               = cms.double(0.0),
    maxTrackTrackDOCA        = cms.double(0.03),## 3D DCA between the two pion tracks (cm)
)

# phi(1020) -> K+K- candidates: symmetric mode (existing interface).
ALCARECOTkAlJpsiXPhiCandidates = cms.EDProducer('TwoBodyDecayCandidateProducer',
    src           = cms.InputTag('generalTracks'),
    muonSrc       = cms.InputTag(''),
    daughterMass  = cms.double(0.493677),  ## K
    daughterPdgId = cms.int32(321),        ## K+
    motherPdgId   = cms.int32(333),        ## phi
    minMass       = cms.double(0.990),  ## at K+K- threshold (2mK = 0.987 GeV)
    maxMass       = cms.double(1.040),  ## tightened upper edge
    minDaughterPt            = cms.double(_NON_V0_PRESETS['Phi']['minDaughterPt']),
    maxDaughterEta           = cms.double(_NON_V0_PRESETS['Phi']['maxDaughterEta']),
    applyChargeFilter        = cms.bool(True),
    charge                   = cms.int32(0),
    useUnsignedCharge        = cms.bool(True),
    applyAcoplanarityFilter  = cms.bool(False),
    acoplanarDistance        = cms.double(1.0),
    applyVertexFit           = cms.bool(False),   ## no Kalman fit at AlCaReco (cost)
    minVtxProb               = cms.double(0.0),
)

# ---------------------------------------------------------------------------
# Stage 2 -- Per-channel B-meson / quarkonium candidate producers
# ---------------------------------------------------------------------------

# B+ -> J/psi K+/-  (track mode: bachelor = kaon from generalTracks)
ALCARECOTkAlJpsiXBPlusCandidates = cms.EDProducer('JpsiXCandidateProducer',
    xMode          = cms.string('track'),
    jpsiSrc        = cms.InputTag('ALCARECOTkAlJpsiXJpsiCandidates'),
    trackSrc       = cms.InputTag('generalTracks'),
    minBachelorPt  = cms.double(_NON_V0_PRESETS['BPlus']['minBachelorPt']),
    bachelorMass   = cms.double(0.493677),  ## kaon
    bachelorPdgId  = cms.int32(321),
    motherPdgId    = cms.int32(521),        ## B+
    minMotherMass  = cms.double(5.0),
    maxMotherMass  = cms.double(5.5),
    minJpsiPt      = cms.double(_NON_V0_PRESETS['BPlus']['minJpsiPt']),
    minMotherPt    = cms.double(_NON_V0_PRESETS['BPlus']['minMotherPt']),
    maxBachelorEta = cms.double(_NON_V0_PRESETS['BPlus']['maxBachelorEta']),
    maxBachelorMuTrackDOCA   = cms.double(_NON_V0_PRESETS['BPlus']['maxBachelorMuTrackDOCA']),
    applyJpsiMassConstraint  = cms.bool(False),  ## flipped to True under preset C only
    ## maxBachelorIPToJpsiVertex / maxMotherAlphaBS / minBVtxProb / minBLxyOverSigma
    ## removed: no Kalman B-vertex fit at AlCaReco (cost). Phase-1-tuned cuts
    ## live on maxBachelorMuTrackDOCA above.
)

# B0 -> J/psi K*0  (VCC mode: intermediate resonance = K*0)
ALCARECOTkAlJpsiXB0KstarCandidates = cms.EDProducer('JpsiXCandidateProducer',
    xMode         = cms.string('vcc'),
    jpsiSrc       = cms.InputTag('ALCARECOTkAlJpsiXJpsiCandidates'),
    xSrc          = cms.InputTag('ALCARECOTkAlJpsiXKstarCandidates'),
    motherPdgId   = cms.int32(511),    ## B0
    minMotherMass = cms.double(5.0),
    maxMotherMass = cms.double(5.5),
    minJpsiPt     = cms.double(_NON_V0_PRESETS['B0Kstar']['minJpsiPt']),
    minMotherPt   = cms.double(_NON_V0_PRESETS['B0Kstar']['minMotherPt']),
    maxBachelorMuTrackDOCA = cms.double(_NON_V0_PRESETS['B0Kstar']['maxBachelorMuTrackDOCA']),
    applyJpsiMassConstraint = cms.bool(False),  ## flipped to True under preset C only
)

# B0 -> J/psi Ks  (VCC mode: Ks from shared V0Producer clone)
ALCARECOTkAlJpsiXB0KsCandidates = cms.EDProducer('JpsiXCandidateProducer',
    xMode         = cms.string('vcc'),
    jpsiSrc       = cms.InputTag('ALCARECOTkAlJpsiXJpsiCandidates'),
    xSrc          = cms.InputTag('ALCARECOTkAlV0Candidates', 'Kshort'),
    motherPdgId   = cms.int32(511),    ## B0
    minMotherMass = cms.double(5.0),
    maxMotherMass = cms.double(5.5),
    minJpsiPt     = cms.double(3.0),
    applyJpsiMassConstraint = cms.bool(False),  ## V0-mode invariant: False under every preset
    ## maxMotherAlphaBS removed: requires Kalman B-vertex fit (forbidden at
    ## AlCaReco). V0-mode channels are preset-invariant; rely on V0Producer
    ## upstream-RECO quality + J/psi mass+pT cuts.
)

# Bs -> J/psi phi  (VCC mode: phi)
ALCARECOTkAlJpsiXBsPhiCandidates = cms.EDProducer('JpsiXCandidateProducer',
    xMode         = cms.string('vcc'),
    jpsiSrc       = cms.InputTag('ALCARECOTkAlJpsiXJpsiCandidates'),
    xSrc          = cms.InputTag('ALCARECOTkAlJpsiXPhiCandidates'),
    motherPdgId   = cms.int32(531),    ## Bs
    minMotherMass = cms.double(5.2),
    maxMotherMass = cms.double(5.6),
    minJpsiPt     = cms.double(_NON_V0_PRESETS['BsPhi']['minJpsiPt']),
    minMotherPt   = cms.double(_NON_V0_PRESETS['BsPhi']['minMotherPt']),
    maxBachelorMuTrackDOCA = cms.double(_NON_V0_PRESETS['BsPhi']['maxBachelorMuTrackDOCA']),
    applyJpsiMassConstraint = cms.bool(False),  ## flipped to True under preset C only
)

# Lambda_b -> J/psi Lambda  (VCC mode: Lambda from shared V0Producer clone)
ALCARECOTkAlJpsiXLambdabCandidates = cms.EDProducer('JpsiXCandidateProducer',
    xMode         = cms.string('vcc'),
    jpsiSrc       = cms.InputTag('ALCARECOTkAlJpsiXJpsiCandidates'),
    xSrc          = cms.InputTag('ALCARECOTkAlV0Candidates', 'Lambda'),
    motherPdgId   = cms.int32(5122),   ## Lambda_b
    minMotherMass = cms.double(5.3),
    maxMotherMass = cms.double(6.0),
    minJpsiPt     = cms.double(3.0),
    applyJpsiMassConstraint = cms.bool(False),  ## V0-mode invariant: False under every preset
    ## maxMotherAlphaBS removed: see B0->Ks comment.
)

# psi(2S) -> J/psi pi+pi-  (VCC mode: prompt pi+pi- from ALCARECOTkAlJpsiXPiPiCandidates).
# Historically wired to V0 Ks output (WRONG -- psi(2S) does not have a
# displaced Ks-quality pi+pi- vertex). Fixed in openspec
# add-jpsi-x-muons-and-preprod-refinements: the dipion producer above
# supplies a prompt-relative-to-psi(2S) pi+pi- pair via
# TwoBodyDecayCandidateProducer. psi(2S) itself may be prompt (from PV,
# dominant fraction) or non-prompt (from B decays, ~20 %); no
# displacement cut is applied under preset B, so both are accepted.
ALCARECOTkAlJpsiXPsi2SCandidates = cms.EDProducer('JpsiXCandidateProducer',
    xMode         = cms.string('vcc'),
    jpsiSrc       = cms.InputTag('ALCARECOTkAlJpsiXJpsiCandidates'),
    xSrc          = cms.InputTag('ALCARECOTkAlJpsiXPiPiCandidates'),
    motherPdgId   = cms.int32(100443), ## psi(2S)
    minMotherMass = cms.double(3.5),
    maxMotherMass = cms.double(3.9),
    minJpsiPt     = cms.double(3.0),
    minMotherPt   = cms.double(3.0),
    applyJpsiMassConstraint = cms.bool(False),  ## preset-invariant: no dimuon constraint
    ## maxMotherAlphaBS removed: no Kalman B-vertex fit at AlCaReco (preset B).
    ## Under preset C, alpha_BS + Lxy sig can be added here to filter to
    ## non-prompt (from B decays) psi(2S) only.
)

# Bc+ -> J/psi pi+/-  (track mode: bachelor = pion from generalTracks)
ALCARECOTkAlJpsiXBcCandidates = cms.EDProducer('JpsiXCandidateProducer',
    xMode          = cms.string('track'),
    jpsiSrc        = cms.InputTag('ALCARECOTkAlJpsiXJpsiCandidates'),
    trackSrc       = cms.InputTag('generalTracks'),
    minBachelorPt  = cms.double(_NON_V0_PRESETS['Bc']['minBachelorPt']),
    bachelorMass   = cms.double(0.139570),  ## pion
    bachelorPdgId  = cms.int32(211),
    motherPdgId    = cms.int32(541),        ## Bc+
    minMotherMass  = cms.double(5.9),
    maxMotherMass  = cms.double(6.6),
    minJpsiPt      = cms.double(_NON_V0_PRESETS['Bc']['minJpsiPt']),
    minMotherPt    = cms.double(_NON_V0_PRESETS['Bc']['minMotherPt']),
    maxBachelorEta = cms.double(_NON_V0_PRESETS['Bc']['maxBachelorEta']),
    maxBachelorMuTrackDOCA   = cms.double(_NON_V0_PRESETS['Bc']['maxBachelorMuTrackDOCA']),
    applyJpsiMassConstraint  = cms.bool(False),  ## flipped to True under preset C only
)

# ---------------------------------------------------------------------------
# Eighth production channel: J/psi-only, dimuon-only, run over the LOOSE muon
# selector (ALCARECOTkAlJpsiXLooseMuons) so it captures muon populations the
# tight selector would reject (tracker-only muons, muons with no globalTrack
# info). Preset-invariant: the muon-quality axis is orthogonal to the
# selection preset. See openspec change add-jpsi-x-muons-and-preprod-refinements.
# ---------------------------------------------------------------------------
ALCARECOTkAlJpsiXJpsiOnlyCandidates = cms.EDProducer('TwoBodyDecayCandidateProducer',
    src           = cms.InputTag('generalTracks'),
    muonSrc       = cms.InputTag('ALCARECOTkAlJpsiXLooseMuons'),
    minMass       = cms.double(2.95),   ## same window as the tight-selector J/psi
    maxMass       = cms.double(3.25),
    daughterMass  = cms.double(0.105),  ## muon
    daughterPdgId = cms.int32(13),      ## mu-
    motherPdgId   = cms.int32(443),     ## J/psi
    applyChargeFilter        = cms.bool(True),
    charge                   = cms.int32(0),
    useUnsignedCharge        = cms.bool(True),
    applyAcoplanarityFilter  = cms.bool(False),
    acoplanarDistance        = cms.double(1.0),
    applyVertexFit           = cms.bool(False),
    minVtxProb               = cms.double(0.0),
)

# ---------------------------------------------------------------------------
# Preset C: inject Kalman vertex-fit cut parameters into the 4 non-V0 producers.
# The C++ JpsiXCandidateProducer activates the Kalman branch when ANY of
# minBVtxProb, maxMotherAlphaBS, or minBLxyOverSigma is present in the python
# config (existence-check, not value-check), so under presets A/B no Kalman
# parameter is set and the producer's fit branch stays inactive. Under preset
# C the parameters are set on all 4 non-V0 channels (BPlus and Bc -> 3-body
# fit on (mu, mu, bachelor); B0Kstar and BsPhi -> 4-body fit on all leaf
# tracks at a common B vertex). V0-mode channels are untouched.
# ---------------------------------------------------------------------------
if _TKALJPSIX_SELECTION_PRESET == 'C':
    for _prod, _channel in (
        (ALCARECOTkAlJpsiXBPlusCandidates,   'BPlus'),
        (ALCARECOTkAlJpsiXBcCandidates,      'Bc'),
        (ALCARECOTkAlJpsiXB0KstarCandidates, 'B0Kstar'),
        (ALCARECOTkAlJpsiXBsPhiCandidates,   'BsPhi'),
    ):
        _prod.minBVtxProb             = cms.double(_NON_V0_PRESETS[_channel]['minBVtxProb'])
        _prod.maxMotherAlphaBS        = cms.double(_NON_V0_PRESETS[_channel]['maxMotherAlphaBS'])
        _prod.minBLxyOverSigma        = cms.double(_NON_V0_PRESETS[_channel]['minBLxyOverSigma'])
        _prod.applyJpsiMassConstraint = cms.bool(True)
    ## V0-mode producers (ALCARECOTkAlJpsiX{B0Ks,Lambdab,Psi2S}Candidates) keep
    ## applyJpsiMassConstraint=False under preset C: they are preset-invariant
    ## by the V0 spec requirement.

# ---------------------------------------------------------------------------
# Stage 3 -- Merged deduplicated leaf-track collection
# ---------------------------------------------------------------------------

ALCARECOTkAlJpsiXAllTracks = cms.EDProducer('CompositeDaughterTrackProducer',
    srcs = cms.VInputTag(
        cms.InputTag('ALCARECOTkAlJpsiXBPlusCandidates'),
        cms.InputTag('ALCARECOTkAlJpsiXB0KstarCandidates'),
        cms.InputTag('ALCARECOTkAlJpsiXB0KsCandidates'),
        cms.InputTag('ALCARECOTkAlJpsiXBsPhiCandidates'),
        cms.InputTag('ALCARECOTkAlJpsiXLambdabCandidates'),
        cms.InputTag('ALCARECOTkAlJpsiXPsi2SCandidates'),
        cms.InputTag('ALCARECOTkAlJpsiXBcCandidates'),
        cms.InputTag('ALCARECOTkAlJpsiXJpsiOnlyCandidates'),
    ),
)

# ---------------------------------------------------------------------------
# Stage 4 -- Track selector + dE/dx projectors
# ---------------------------------------------------------------------------

import Alignment.CommonAlignmentProducer.AlignmentTrackSelectorWithIndexMap_cfi
ALCARECOTkAlJpsiX = Alignment.CommonAlignmentProducer.AlignmentTrackSelectorWithIndexMap_cfi.AlignmentTrackSelectorWithIndexMap.clone(
    src            = cms.InputTag('ALCARECOTkAlJpsiXAllTracks'),
    filter         = True,   ## reject events with no candidates
    applyBasicCuts = True,
    ptMin          = 0.1,    ## GeV; matches every upstream pT cut + V0 clone tkPtCut
    etaMin         = -3.5,
    etaMax         = 3.5,
    nHitMin        = 0,
)
ALCARECOTkAlJpsiX.GlobalSelector.applyGlobalMuonFilter = False
ALCARECOTkAlJpsiX.GlobalSelector.applyIsolationtest    = False

# dE/dx projected onto all tracks in the merged collection.
# Physically meaningful PID information: kaon channels (B+, B0->K*0, Bs)
# and proton channel (Lambda_b). Pion and muon tracks also get values
# but these carry no useful PID information.
from Alignment.CommonAlignmentProducer.alcaDedxJointEstimator_cfi import alcaDedxJointEstimator
ALCARECOTkAlJpsiXDeDxHarmonic2 = cms.EDProducer('DeDxValueMapProjector',
    selectedTracks     = cms.InputTag('ALCARECOTkAlJpsiX'),
    intermediateTracks = cms.InputTag('ALCARECOTkAlJpsiXAllTracks'),
    sourceTracks       = cms.InputTag('generalTracks'),
    sourceValueMap     = cms.InputTag('dedxHarmonic2'),
    originalIndexMap   = cms.InputTag('ALCARECOTkAlJpsiX', 'originalIndex'),
)
ALCARECOTkAlJpsiXDeDxPixelHarmonic2 = ALCARECOTkAlJpsiXDeDxHarmonic2.clone(
    sourceValueMap = cms.InputTag('dedxPixelHarmonic2'),
)
ALCARECOTkAlJpsiXDeDxAllHarmonic2 = ALCARECOTkAlJpsiXDeDxHarmonic2.clone(
    sourceValueMap = cms.InputTag('alcaDedxJointEstimator'),
)

# ---------------------------------------------------------------------------
# Stage 5 -- Per-channel candidate remappers
# ---------------------------------------------------------------------------

ALCARECOTkAlJpsiXBPlusResonances = cms.EDProducer('VertexCompositeCandidateRemapper',
    srcCandidates      = cms.InputTag('ALCARECOTkAlJpsiXBPlusCandidates'),
    selectedTracks     = cms.InputTag('ALCARECOTkAlJpsiX'),
    intermediateTracks = cms.InputTag('ALCARECOTkAlJpsiXAllTracks'),
    originalIndexMap   = cms.InputTag('ALCARECOTkAlJpsiX', 'originalIndex'),
)
ALCARECOTkAlJpsiXB0KstarResonances = cms.EDProducer('VertexCompositeCandidateRemapper',
    srcCandidates      = cms.InputTag('ALCARECOTkAlJpsiXB0KstarCandidates'),
    selectedTracks     = cms.InputTag('ALCARECOTkAlJpsiX'),
    intermediateTracks = cms.InputTag('ALCARECOTkAlJpsiXAllTracks'),
    originalIndexMap   = cms.InputTag('ALCARECOTkAlJpsiX', 'originalIndex'),
)
ALCARECOTkAlJpsiXB0KsResonances = cms.EDProducer('VertexCompositeCandidateRemapper',
    srcCandidates      = cms.InputTag('ALCARECOTkAlJpsiXB0KsCandidates'),
    selectedTracks     = cms.InputTag('ALCARECOTkAlJpsiX'),
    intermediateTracks = cms.InputTag('ALCARECOTkAlJpsiXAllTracks'),
    originalIndexMap   = cms.InputTag('ALCARECOTkAlJpsiX', 'originalIndex'),
)
ALCARECOTkAlJpsiXBsPhiResonances = cms.EDProducer('VertexCompositeCandidateRemapper',
    srcCandidates      = cms.InputTag('ALCARECOTkAlJpsiXBsPhiCandidates'),
    selectedTracks     = cms.InputTag('ALCARECOTkAlJpsiX'),
    intermediateTracks = cms.InputTag('ALCARECOTkAlJpsiXAllTracks'),
    originalIndexMap   = cms.InputTag('ALCARECOTkAlJpsiX', 'originalIndex'),
)
ALCARECOTkAlJpsiXLambdabResonances = cms.EDProducer('VertexCompositeCandidateRemapper',
    srcCandidates      = cms.InputTag('ALCARECOTkAlJpsiXLambdabCandidates'),
    selectedTracks     = cms.InputTag('ALCARECOTkAlJpsiX'),
    intermediateTracks = cms.InputTag('ALCARECOTkAlJpsiXAllTracks'),
    originalIndexMap   = cms.InputTag('ALCARECOTkAlJpsiX', 'originalIndex'),
)
ALCARECOTkAlJpsiXPsi2SResonances = cms.EDProducer('VertexCompositeCandidateRemapper',
    srcCandidates      = cms.InputTag('ALCARECOTkAlJpsiXPsi2SCandidates'),
    selectedTracks     = cms.InputTag('ALCARECOTkAlJpsiX'),
    intermediateTracks = cms.InputTag('ALCARECOTkAlJpsiXAllTracks'),
    originalIndexMap   = cms.InputTag('ALCARECOTkAlJpsiX', 'originalIndex'),
)
ALCARECOTkAlJpsiXBcResonances = cms.EDProducer('VertexCompositeCandidateRemapper',
    srcCandidates      = cms.InputTag('ALCARECOTkAlJpsiXBcCandidates'),
    selectedTracks     = cms.InputTag('ALCARECOTkAlJpsiX'),
    intermediateTracks = cms.InputTag('ALCARECOTkAlJpsiXAllTracks'),
    originalIndexMap   = cms.InputTag('ALCARECOTkAlJpsiX', 'originalIndex'),
)
ALCARECOTkAlJpsiXJpsiOnlyResonances = cms.EDProducer('VertexCompositeCandidateRemapper',
    srcCandidates      = cms.InputTag('ALCARECOTkAlJpsiXJpsiOnlyCandidates'),
    selectedTracks     = cms.InputTag('ALCARECOTkAlJpsiX'),
    intermediateTracks = cms.InputTag('ALCARECOTkAlJpsiXAllTracks'),
    originalIndexMap   = cms.InputTag('ALCARECOTkAlJpsiX', 'originalIndex'),
)

# Track -> reco::Muon association keyed on the persisted deduplicated track
# collection ALCARECOTkAlJpsiX and valued into ALCARECOTkAlJpsiXLooseMuons (a
# superset of the tight selector's output, so a single association covers
# both). Non-muon tracks (bachelors, V0 daughters, K*0 / phi daughters) map
# to null MuonRef. See plugin AlignmentTrackToMuonAssociator.
ALCARECOTkAlJpsiXTrackToMuon = cms.EDProducer('AlignmentTrackToMuonAssociator',
    selectedTracks     = cms.InputTag('ALCARECOTkAlJpsiX'),
    intermediateTracks = cms.InputTag('ALCARECOTkAlJpsiXAllTracks'),
    originalIndexMap   = cms.InputTag('ALCARECOTkAlJpsiX', 'originalIndex'),
    muons              = cms.InputTag('ALCARECOTkAlJpsiXLooseMuons'),
)

# ---------------------------------------------------------------------------
# Combined sequence (Stage 0 -> 1 -> 2 -> 3 -> 4 -> 5)
# Note: ALCARECOTkAlJpsiXHLT is defined above but intentionally excluded.
# ---------------------------------------------------------------------------
seqALCARECOTkAlJpsiX = cms.Sequence(
    ALCARECOTkAlJpsiXDCSFilter +
    ALCARECOTkAlJpsiXGoodMuons +
    ALCARECOTkAlJpsiXLooseMuons +
    ALCARECOTkAlJpsiXJpsiCandidates +
    ALCARECOTkAlJpsiXJpsiOnlyCandidates +
    ALCARECOTkAlV0Candidates +
    ALCARECOTkAlJpsiXKstarCandidates +
    ALCARECOTkAlJpsiXPhiCandidates +
    ALCARECOTkAlJpsiXPiPiCandidates +
    ALCARECOTkAlJpsiXBPlusCandidates +
    ALCARECOTkAlJpsiXB0KstarCandidates +
    ALCARECOTkAlJpsiXB0KsCandidates +
    ALCARECOTkAlJpsiXBsPhiCandidates +
    ALCARECOTkAlJpsiXLambdabCandidates +
    ALCARECOTkAlJpsiXPsi2SCandidates +
    ALCARECOTkAlJpsiXBcCandidates +
    ALCARECOTkAlJpsiXAllTracks +
    ALCARECOTkAlJpsiX +
    ALCARECOTkAlJpsiXTrackToMuon +
    alcaDedxJointEstimator +
    ALCARECOTkAlJpsiXDeDxHarmonic2 +
    ALCARECOTkAlJpsiXDeDxPixelHarmonic2 +
    ALCARECOTkAlJpsiXDeDxAllHarmonic2 +
    ALCARECOTkAlJpsiXBPlusResonances +
    ALCARECOTkAlJpsiXB0KstarResonances +
    ALCARECOTkAlJpsiXB0KsResonances +
    ALCARECOTkAlJpsiXBsPhiResonances +
    ALCARECOTkAlJpsiXLambdabResonances +
    ALCARECOTkAlJpsiXPsi2SResonances +
    ALCARECOTkAlJpsiXBcResonances +
    ALCARECOTkAlJpsiXJpsiOnlyResonances
)
