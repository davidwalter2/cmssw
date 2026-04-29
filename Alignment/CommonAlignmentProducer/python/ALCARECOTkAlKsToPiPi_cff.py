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

# Event pre-filter: require >=1 reconstructed KS candidate.
# generalV0Candidates:Kshort is produced by the standard V0Producer in RECO and
# already has displaced-vertex cuts applied (flight significance, pointing angle,
# IP significance). This suppresses prompt-track combinatorial background without
# requiring a new C++ module or any specific HLT trigger.
# No HLT filter is applied so the ALCARECO is usable from any trigger path
# (single muon, dimuon, etc.) and is suitable for multiple physics channels
# including W/Z+KS, Psi(2S)->J/psiKS, and inclusive KS production.
ALCARECOTkAlKsToPiPiV0Filter = cms.EDFilter('CandViewCountFilter',
    src = cms.InputTag('generalV0Candidates', 'Kshort'),
    minNumber = cms.uint32(1)
)

# Store the pion tracks of every reconstructed KS candidate, with cloned hits
# and clusters. AlignmentTracksFromV0Selector extracts the unique daughter
# TrackRefs from generalV0Candidates:Kshort and uses the standard
# helper::TrackCollectionStoreManager machinery to clone tracks + extras + hits
# + clusters into self-contained collections. This avoids the combinatorial
# track-pair search done by AlignmentTwoBodyDecayTrackSelector and ensures all
# V0Producer KS daughters (with their displaced-vertex cuts already applied)
# are preserved. filter() always returns true; gating is done upstream by
# the V0 CandViewCountFilter.
ALCARECOTkAlKsToPiPi = cms.EDFilter('AlignmentTracksFromV0Selector',
    src    = cms.InputTag('generalTracks'),                         ## main input (ObjectSelector convention)
    v0src  = cms.InputTag('generalV0Candidates', 'Kshort'),         ## V0 candidates whose daughter tracks to keep
    filter = cms.bool(True),                                         ## drop event if zero daughters survive
)

seqALCARECOTkAlKsToPiPi = cms.Sequence(
    ALCARECOTkAlKsToPiPiDCSFilter +
    ALCARECOTkAlKsToPiPiV0Filter +
    ALCARECOTkAlKsToPiPi
)
