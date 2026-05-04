# Local V0Producer clone for the TkAlKsToPiPi and TkAlLambdaToProtonPi
# ALCARECOs.
#
# generalV0Candidates is consumed by many downstream modules (DQM, MINIAOD
# slimming, PF candidate-PV association, Onia analyses, the standard Geant4e
# refit, ...). Modifying its cuts globally would change behaviour everywhere.
# Instead we keep the standard generalV0Candidates untouched and run a second
# instance under a private label, with looser cuts tuned for the ALCARECO
# (in particular a lower track-pT threshold to capture the soft V0 tail).
#
# Both KS and Lambda are produced by the same V0Producer instance
# (doKShorts/doLambdas both True), so this single clone serves both ALCARECOs.
# The framework deduplicates module execution per event, so even though both
# seqALCARECOTkAlKsToPiPi and seqALCARECOTkAlLambdaToProtonPi reference it,
# it runs at most once per event.

import FWCore.ParameterSet.Config as cms
from RecoVertex.V0Producer.generalV0Candidates_cff import generalV0Candidates

ALCARECOTkAlV0Candidates = generalV0Candidates.clone(
    tkPtCut = 0.1,   ## lower than the default 0.35 GeV
)
