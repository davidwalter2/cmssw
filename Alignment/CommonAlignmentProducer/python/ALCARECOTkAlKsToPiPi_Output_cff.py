# Output definition for ALCARECOTkAlKsToPiPi
import copy
import FWCore.ParameterSet.Config as cms

OutALCARECOTkAlKsToPiPi_noDrop = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pathALCARECOTkAlKsToPiPi')
    ),
    outputCommands = cms.untracked.vstring(
        'keep *_ALCARECOTkAlKsToPiPi_*_*',
        'keep *_generalV0Candidates_Kshort_*',
        'keep L1AcceptBunchCrossings_*_*_*',
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*',
        'keep *_TriggerResults_*_*',
        'keep DcsStatuss_scalersRawToDigi_*_*',
        'keep recoBeamSpot_offlineBeamSpot_*_*',
        'keep *_offlinePrimaryVertices_*_*',
        # generalTracks + dE/dx: needed to access dE/dx info for V0 daughter tracks
        # via the daughter().track() TrackRef -> dE/dx ValueMap key chain.
        'keep *_generalTracks_*_*',
        'keep *_dedxHarmonic2_*_*',
        'keep *_dedxPixelHarmonic2_*_*',
        'keep *_dedxHitInfo_*_*',
    )
)
OutALCARECOTkAlKsToPiPi = copy.deepcopy(OutALCARECOTkAlKsToPiPi_noDrop)
OutALCARECOTkAlKsToPiPi.outputCommands.insert(0, 'drop *')
