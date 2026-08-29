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
                                   # The CDF fraction of the delta-electron
                                   # spectrum kept when the ionization variance
                                   # is formed.  This is a convention, and it
                                   # exists only because a Gaussian weight needs
                                   # a second moment while the Urban delta
                                   # channel's does not converge.  It is read
                                   # ONLY under `CgfQoPMode = 0` (the legacy
                                   # truncated-Q refit); under the default
                                   # Fisher weight the cut has no meaning and
                                   # this value is never used.  For the record:
                                   # scanning it under the legacy weight moves
                                   # the fitted q/p by rms 1.2e-5, i.e. the
                                   # Z-mass target itself.
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
                                   # Log-T buckets of the Kokoulin term in the
                                   # BLOCK CGF (cvhcgf), which is evaluated at
                                   # every point of the inversion grid and is
                                   # therefore the one place where this count
                                   # costs real time -- unlike IoniKokoulinNbin
                                   # above, which is a once-per-step Simpson
                                   # rule.  Measured on the pT = 3 reference
                                   # ray: against nbin = 96, the block 1/I
                                   # moves by 1.6e-4 (8), 3.8e-5 (16), 1.6e-5
                                   # (24), 3.3e-6 (48) relative, while the
                                   # correction ITSELF is +1.2e-3 .. +3.3e-3 of
                                   # 1/I.  16 therefore reproduces the
                                   # correction to ~1 % of itself at a sixth of
                                   # the cost.  Not Simpson: any value >= 1.
                                   #
                                   # DEFAULT 0 = OMIT THE TERM, a cost decision
                                   # with the physics measured both ways: 25
                                   # events with the CGF weight take 96.6 s with
                                   # the exact-delta channel and no Kokoulin --
                                   # indistinguishable from the 95.2 s of the
                                   # tree-level channel -- and more than 20x
                                   # that with 16 buckets, because the bucketing
                                   # is evaluated at every point of every
                                   # block's inversion grid.  Against that, the
                                   # term moves 1/I by +1.2e-3 .. +3.3e-3 at
                                   # pT = 3 and +2.1e-3 .. +8.6e-3 at pT = 40
                                   # (f_K grows logarithmically with E).
                                   #
                                   # MEASURED IN THE FIT, 48 tracks, nbin 8 vs
                                   # off: the fitted q/p moves by mean
                                   # +2.6e-07 with an rms of 2.4e-06 -- and the
                                   # rms is NOT the channel. A control that
                                   # changes 1/I by 1e-16 (the same sum
                                   # regrouped) already gives an rms of
                                   # 2.4e-07, and nbin 8 vs 4 -- a far smaller
                                   # physics change than 8 vs off -- gives a
                                   # LARGER 2.7e-06. The rms at this level is
                                   # the fit's conditioning; only the mean is
                                   # the correction. 2.6e-07 is 40x below the
                                   # Z-mass target.
                                   #
                                   # Cost, same machine, nbin 8 with the
                                   # shared-edge form: 162 s -> 560 s, i.e.
                                   # +245 %. That is the trade: a factor 3.5 in
                                   # wall clock for 2.6e-07 on a momentum.
                                   IoniKokoulinCgfNbin=cms.int32(0),

                                   # The RADIATIVE (bremsstrahlung + pair)
                                   # channel of the block CGF.  Implemented and
                                   # validated against the offline reference
                                   # (S(t) to 1.6e-9, 1/I to 8.5e-7), and
                                   # DEFAULT OFF on a measurement:
                                   #
                                   #   d(1/I) from the channel, per leg
                                   #     pT = 3     +1.3e-5 .. +6.4e-5
                                   #     pT = 40    +6.3e-5 .. +8.1e-4
                                   #     pT = 100   +6.1e-5 .. +3.0e-4
                                   #
                                   # and, measured in the fit, a q/p shift of
                                   # mean -1.4e-07 with an rms of 9.4e-07 (the
                                   # rms is largely the conditioning floor --
                                   # see IoniKokoulinCgfNbin above -- so the
                                   # channel is worth ~1e-07 on a momentum).  Note this is NOT the 18.6 % of
                                   # the q/p kappa2 that radiative carries at
                                   # pT = 40: Fisher information is
                                   # core-dominated, so a tail channel with a
                                   # fifth of the VARIANCE carries a thousandth
                                   # of the INFORMATION.  What it does move is
                                   # the block's MODE (+0.3 % at pT = 3), which
                                   # is what the mode-3 re-centring uses -- so
                                   # it exists, it is validated, and it is one
                                   # switch away.
                                   CgfRadiativeChannel=cms.bool(False),

                                   # THE PROCESS-NOISE WEIGHT FOR q/p.
                                   #
                                   #   0 = THE LEGACY CVH REFIT.  The Gaussian
                                   #       weight on the delta-truncated
                                   #       variance, with `IoniTruncationAlpha`
                                   #       live -- i.e. the estimator this fit
                                   #       was built and validated on, restored
                                   #       bit for bit (verified: 95/95 branches
                                   #       identical to a build at 663639e, the
                                   #       commit before the CGF default).  It
                                   #       is a supported configuration and not
                                   #       merely a diagnostic limit: it carries
                                   #       an unphysical convention the Fisher
                                   #       weight does not, but it costs
                                   #       114 ms/track against the CGF's
                                   #       1562 ms/track -- a factor 13.7 --
                                   #       and the CGF has NO measured accuracy
                                   #       advantage over it against gen truth
                                   #       (NOTES_CGFFIT s88, s89).
                                   #   1 = the FISHER INFORMATION of the block,
                                   #       by exact inversion of its
                                   #       characteristic function.  THE
                                   #       DEFAULT.
                                   #   2 = 1 plus the per-leg diagnostic print.
                                   #   3 = 1 plus the IRLS re-centring, a
                                   #       prototype whose fixed point is still
                                   #       schedule-dependent at 1.6e-4
                                   #       (NOTES_CGFFIT s58).
                                   #
                                   # Why the default moved to 1: the Urban
                                   # delta channel has a 1/E^2 single-collision
                                   # spectrum, so its second moment does not
                                   # exist without a cut, and that cut was a
                                   # convention worth ~1e-5 on the fitted q/p
                                   # -- the Z-mass target itself.  The Fisher
                                   # information of the same distribution needs
                                   # no cut.  Measured: schedule-independent to
                                   # 1.4e-6, no convergence cost (8.60 against
                                   # 8.56 mean iterations), +11 % wall clock.
                                   CgfQoPMode=cms.int32(1),

                                   # Recompute the block every N Gauss-Newton
                                   # sweeps; 0 = freeze after the first.
                                   # Freezing is justified, not assumed: the
                                   # fixed point is the same to 1.4e-6 (five
                                   # orders on the EDM target only shrink the
                                   # difference, i.e. it is a stopping residue,
                                   # not a second fixed point) and it is 14x
                                   # cheaper.
                                   CgfQoPRefresh=cms.int32(0),

                                   # Shrinkage on the mode-3 re-centring.
                                   # 1.0 = full size (what stage 7 built and
                                   # what section 85 measured), 0 = none.
                                   # Against gen truth the correction is
                                   # anti-correlated with the residual --
                                   # 2Cov/Var = -5.4e-3, the only thing in this
                                   # programme that is -- while adding +6.8e-3
                                   # of scatter, so at full size it is a wash.
                                   # The variance-minimising value is
                                   # -Cov/Var(delta) = 0.40, predicted to buy
                                   # -1.1e-3 on the residual variance.
                                   CgfRecentreDamping=cms.double(1.0),
                                   # MeV; 0 means "use e0".  A scan knob for the
                                   # exact-delta channel's lower bound, not a
                                   # tune -- nothing is fitted to it.
                                   IoniExactDeltaT0=cms.double(0.0)
                                   )
