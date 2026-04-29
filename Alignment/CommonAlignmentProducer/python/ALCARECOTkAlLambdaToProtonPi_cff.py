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

# Event pre-filter: require >=1 reconstructed Lambda0 candidate.
# generalV0Candidates:Lambda is produced by the standard V0Producer in RECO with
# full displaced-vertex cuts already applied.
ALCARECOTkAlLambdaToProtonPiV0Filter = cms.EDFilter('CandViewCountFilter',
    src = cms.InputTag('generalV0Candidates', 'Lambda'),
    minNumber = cms.uint32(1)
)

# Store the proton + pion tracks of every reconstructed Lambda0 candidate, with
# cloned hits and clusters. AlignmentTracksFromV0Selector extracts the unique
# daughter TrackRefs from generalV0Candidates:Lambda and uses the standard
# helper::TrackCollectionStoreManager machinery to clone tracks + extras + hits
# + clusters into self-contained collections. By V0Producer convention,
# daughter(0) = baryon (p or pbar) and daughter(1) = pion (pi- or pi+); both
# Lambda0 and anti-Lambda0 are accepted automatically. filter() always returns
# true; gating is done upstream by the V0 CandViewCountFilter.
ALCARECOTkAlLambdaToProtonPi = cms.EDFilter('AlignmentTracksFromV0Selector',
    src    = cms.InputTag('generalTracks'),                         ## main input (ObjectSelector convention)
    v0src  = cms.InputTag('generalV0Candidates', 'Lambda'),         ## V0 candidates whose daughter tracks to keep
    filter = cms.bool(True),                                         ## drop event if zero daughters survive
)

seqALCARECOTkAlLambdaToProtonPi = cms.Sequence(
    ALCARECOTkAlLambdaToProtonPiDCSFilter +
    ALCARECOTkAlLambdaToProtonPiV0Filter +
    ALCARECOTkAlLambdaToProtonPi
)
