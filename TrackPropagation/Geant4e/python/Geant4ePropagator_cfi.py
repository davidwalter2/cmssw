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
                                   ForCVH=cms.bool(False),
                                   # CDF fraction of the delta-electron spectrum kept in the
                                   # truncated ionization variance (G4UniversalFluctuationFor-
                                   # Extrapolator). 0.999 = historical CVH baseline. Scanned by
                                   # the resolution-closure diagnostic (0.995-0.999): fitted
                                   # resolution parameters must not depend on it.
                                   IoniTruncationAlpha=cms.double(0.999),

                                   # The CVH energy-loss corrections.  These
                                   # were CVH_* ENVIRONMENT VARIABLES until
                                   # 2026-08-18; they are configuration now, so
                                   # that they land in the output provenance
                                   # (edmProvDump recovers exactly which
                                   # corrections produced a file), so a typo is
                                   # a configuration error instead of a silent
                                   # default, and so the defaults live in one
                                   # place.  cvhcgf::configure reads them in
                                   # GeantPropagatorESProducer's constructor.
                                   #
                                   # DEFAULT ON: each exists because the model
                                   # was missing something Geant4 actually
                                   # runs, so the default configuration should
                                   # model the simulation.  Turn one off to
                                   # ATTRIBUTE it, not to get the baseline.
                                   IoniExactDelta=cms.bool(True),
                                   IoniKokoulin=cms.bool(True),
                                   ReferenceChargeAware=cms.bool(True),
                                   ReferenceSpeciesDedx=cms.bool(True),
                                   ReferenceHadronRadiative=cms.bool(True),

                                   # Diagnostics.  DEFAULT OFF -- these do NOT
                                   # move the model toward the simulation.
                                   # IoniUrban2021 in particular was built for
                                   # a prediction NOTES_DELTASPEC s10.4
                                   # falsified.
                                   ReferenceIonizationOnly=cms.bool(False),
                                   IoniUrban2021=cms.bool(False),
                                   # Dump the G4EmParameters block (the MODEL
                                   # half of the SIM/MODEL pair; the SIM half is
                                   # ProcessActivationWatcher's own parameter),
                                   # and set the sim's measured EM values in the
                                   # model job.  Both live in Geant4 classes
                                   # with no PSet of their own, so they ride
                                   # here.
                                   DumpEmParameters=cms.bool(False),
                                   EmHarmonise=cms.bool(False),
                                   # dE/dx table scale; a probe for the J/psi
                                   # mass bias, not a tune.  1.0 = unscaled.
                                   DedxScale=cms.double(1.0),

                                   # Simpson interval counts; must be even and
                                   # >= 2, enforced in cvhcgf::configure.
                                   ReferenceSpeciesDedxNbin=cms.int32(16),
                                   IoniKokoulinNbin=cms.int32(96),
                                   # MeV; 0 means "use e0".  A scan knob for the
                                   # exact-delta channel's lower bound, not a
                                   # tune -- nothing is fitted to it.
                                   IoniExactDeltaT0=cms.double(0.0)
                                   )
