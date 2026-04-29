# Output definition for ALCARECOTkAlLambdaToProtonPi
import copy
import FWCore.ParameterSet.Config as cms

OutALCARECOTkAlLambdaToProtonPi_noDrop = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pathALCARECOTkAlLambdaToProtonPi')
    ),
    outputCommands = cms.untracked.vstring(
        'keep *_ALCARECOTkAlLambdaToProtonPi_*_*',
        'keep *_generalV0Candidates_Lambda_*',
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
OutALCARECOTkAlLambdaToProtonPi = copy.deepcopy(OutALCARECOTkAlLambdaToProtonPi_noDrop)
OutALCARECOTkAlLambdaToProtonPi.outputCommands.insert(0, 'drop *')
