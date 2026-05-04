import FWCore.ParameterSet.Config as cms

# AlCaReco for track based alignment using D*->D0 pi events
OutALCARECOTkAlDstToD0Pi_noDrop = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pathALCARECOTkAlDstToD0Pi')
    ),
    outputCommands = cms.untracked.vstring(
        'keep *_ALCARECOTkAlDstToD0Pi_*_*',
        'keep *_ALCARECOTkAlDstToD0PiDeDxHarmonic2_*_*',         ## per-track dE/dx Harmonic2 (strip)
        'keep *_ALCARECOTkAlDstToD0PiDeDxPixelHarmonic2_*_*',    ## per-track dE/dx Harmonic2 (pixel-only)
        'keep *_ALCARECOTkAlDstToD0PiDeDxAllHarmonic2_*_*',      ## per-track dE/dx Harmonic2-truncated (strip+pixel joint)
        'keep L1AcceptBunchCrossings_*_*_*',
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*',
        'keep *_TriggerResults_*_*',
        'keep DcsStatuss_scalersRawToDigi_*_*',
        'keep *_offlinePrimaryVertices_*_*')
)

import copy
OutALCARECOTkAlDstToD0Pi = copy.deepcopy(OutALCARECOTkAlDstToD0Pi_noDrop)
OutALCARECOTkAlDstToD0Pi.outputCommands.insert(0, "drop *")
