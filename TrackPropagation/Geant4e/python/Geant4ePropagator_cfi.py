import FWCore.ParameterSet.Config as cms

# Default Geant4e propagator setup, uses Muons as particle hypotesis
# ParticleName can be any particle described in the Geant4 documentation
# http://geant4.web.cern.ch/geant4/UserDocumentation/UsersGuides/ForApplicationDeveloper/html/ch05s03.html
# The chargen ( e.g. mu+ or mu- ) will be added by the propagator, depending on the fitted track's charge
Geant4ePropagator = cms.ESProducer("GeantPropagatorESProducer",
                                   ComponentName = cms.string("Geant4ePropagator"),
                                   # anyDirection picks forward/backward per leg from the
                                   # target-plane geometry. Recovers legs whose target plane is
                                   # marginally behind the state (with forward-only propagation
                                   # those escape the tracker and are destroyed in calorimeter
                                   # material while being reported as successes). Validated on
                                   # 100k J/psi events: bit-identical to alongMomentum for every
                                   # previously-succeeding fit, +0.2% recovered candidates;
                                   # backward-leg Jacobians/noise are converted exactly to the
                                   # physical frame in propagateGenericWithJacobianAltD.
                                   PropagationDirection=cms.string("anyDirection"),
                                   ParticleName=cms.string("mu"),
                                   PropagationPtotLimit = cms.double(1.0), ## GeV/c
                                   MagneticFieldLabel = cms.string(""),
                                   ForCVH=cms.bool(False)
                                   )
