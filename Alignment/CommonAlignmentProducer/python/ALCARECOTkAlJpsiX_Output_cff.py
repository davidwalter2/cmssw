import copy
import FWCore.ParameterSet.Config as cms

# AlCaReco output for TkAlJpsiX: B-meson and quarkonium alignment channels.
#
# Saved per event:
#   - Cloned leaf tracks (+ extras + hits + clusters) from all 7 channels,
#     deduplicated into a single ALCARECOTkAlJpsiX collection.
#   - 7 per-channel candidate collections with daughter TrackRefs re-keyed
#     onto ALCARECOTkAlJpsiX (no dangling generalTracks refs).
#   - 3 dE/dx ValueMaps projected onto the cloned track collection
#     (physically useful for kaon channels B+/B0->K*0/Bs and proton
#     channel Lambda_b; saved for all tracks for uniformity).
#   - Standard auxiliary: trigger results, DCS status, primary vertices.

OutALCARECOTkAlJpsiX_noDrop = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pathALCARECOTkAlJpsiX')
    ),
    outputCommands = cms.untracked.vstring(
        'keep *_ALCARECOTkAlJpsiX_*_*',                       ## cloned tracks + extras + hits + clusters
        'keep *_ALCARECOTkAlJpsiXDeDxHarmonic2_*_*',          ## dE/dx strip Harmonic2
        'keep *_ALCARECOTkAlJpsiXDeDxPixelHarmonic2_*_*',     ## dE/dx pixel Harmonic2
        'keep *_ALCARECOTkAlJpsiXDeDxAllHarmonic2_*_*',       ## dE/dx joint strip+pixel
        'keep *_ALCARECOTkAlJpsiXBPlusResonances_*_*',        ## B+ -> J/psi K
        'keep *_ALCARECOTkAlJpsiXB0KstarResonances_*_*',      ## B0 -> J/psi K*0
        'keep *_ALCARECOTkAlJpsiXB0KsResonances_*_*',         ## B0 -> J/psi Ks
        'keep *_ALCARECOTkAlJpsiXBsPhiResonances_*_*',        ## Bs -> J/psi phi
        'keep *_ALCARECOTkAlJpsiXLambdabResonances_*_*',      ## Lambda_b -> J/psi Lambda
        'keep *_ALCARECOTkAlJpsiXPsi2SResonances_*_*',        ## psi(2S) -> J/psi Ks
        'keep *_ALCARECOTkAlJpsiXBcResonances_*_*',           ## Bc -> J/psi pi
        'keep L1AcceptBunchCrossings_*_*_*',
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*',
        'keep *_TriggerResults_*_*',
        'keep DcsStatuss_scalersRawToDigi_*_*',
        'keep *_offlinePrimaryVertices_*_*',
    )
)

OutALCARECOTkAlJpsiX = copy.deepcopy(OutALCARECOTkAlJpsiX_noDrop)
OutALCARECOTkAlJpsiX.outputCommands.insert(0, 'drop *')
