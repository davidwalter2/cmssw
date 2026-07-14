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
# (doKShorts/doLambdas both True), so this single clone serves both ALCARECOs
# and the three J/psi + X V0-mode channels (B0->Ks, Lambda_b, psi(2S)).
# The framework deduplicates module execution per event, so even though both
# seqALCARECOTkAlKsToPiPi and seqALCARECOTkAlLambdaToProtonPi reference it,
# it runs at most once per event.
#
# openspec change add-jpsi-x-muons-and-preprod-refinements: loosen the vertex
# and IP-significance cuts on the local clone (only) to recover V0 yield for
# the alignment use case. All three V0-mode J/psi + X channels get the boost
# together since they share this handle. A future proposal may scan for a
# better working point; this proposal just lands a defensible loosened
# default.

import FWCore.ParameterSet.Config as cms
from RecoVertex.V0Producer.generalV0Candidates_cff import generalV0Candidates

ALCARECOTkAlV0Candidates = generalV0Candidates.clone(
    tkPtCut          = 0.1,     ## lower than the default 0.35 GeV
    vtxChi2Cut       = 15.,     ## loosened from default 6.63
    vtxDecaySigXYCut = 10.,     ## loosened from default 15
    cosThetaXYCut    = 0.995,   ## loosened from default 0.998
    tkIPSigXYCut     = 1.,      ## loosened from default 2
)
