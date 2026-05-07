# Output definition for ALCARECOTkAlLambdaToProtonPi
#
# Minimal output: the cloned V0 daughter tracks (with their extras, hits and
# clusters), the per-track dE/dx scalar estimators (Harmonic2 strip + pixel)
# re-keyed onto the cloned tracks, and standard auxiliary collections. The
# V0 candidates and full generalTracks collection are NOT kept -- step 2
# re-pairs the stored tracks (with vertex constraint) and can recompute any
# additional per-track quantities directly from the cloned hits + clusters.
import copy
import FWCore.ParameterSet.Config as cms

OutALCARECOTkAlLambdaToProtonPi_noDrop = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pathALCARECOTkAlLambdaToProtonPi')
    ),
    outputCommands = cms.untracked.vstring(
        'keep *_ALCARECOTkAlLambdaToProtonPi_*_*',                      ## V0 daughter tracks + extras + hits + clusters
        'keep *_ALCARECOTkAlLambdaToProtonPiResonances_*_*',            ## Lambda candidates with daughter refs into the cloned tracks
        'keep *_ALCARECOTkAlLambdaToProtonPiDeDxHarmonic2_*_*',         ## per-track dE/dx Harmonic2 (strip)
        'keep *_ALCARECOTkAlLambdaToProtonPiDeDxPixelHarmonic2_*_*',    ## per-track dE/dx Harmonic2 (pixel-only)
        'keep *_ALCARECOTkAlLambdaToProtonPiDeDxAllHarmonic2_*_*',      ## per-track dE/dx Harmonic2-truncated (strip+pixel joint)
        'keep L1AcceptBunchCrossings_*_*_*',
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*',
        'keep *_TriggerResults_*_*',                                     ## downstream HLT selection
        'keep DcsStatuss_scalersRawToDigi_*_*',
        'keep *_offlinePrimaryVertices_*_*',                             ## downstream PV constraints
    )
)
OutALCARECOTkAlLambdaToProtonPi = copy.deepcopy(OutALCARECOTkAlLambdaToProtonPi_noDrop)
OutALCARECOTkAlLambdaToProtonPi.outputCommands.insert(0, 'drop *')
