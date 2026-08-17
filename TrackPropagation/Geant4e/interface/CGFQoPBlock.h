#ifndef TrackPropagation_Geant4e_CGFQoPBlock_h
#define TrackPropagation_Geant4e_CGFQoPBlock_h

// Cumulant-generating-function (CGF) description of the q/p process-noise
// block, and its FISHER INFORMATION.
//
// WHAT THIS IS FOR
// ----------------
// The CVH fit weights each process-noise block by the inverse of a VARIANCE
// (Geant4ePropagator::computeErrorIoni, accumulated into errMSIout(0,0)).
// For the Urban ionization model that variance does not exist without a
// truncation: the delta-ray channel has a 1/E^2 single-collision spectrum, so
// the second moment is dominated by the hard edge and is only finite because
// G4UniversalFluctuationForExtrapolator truncates the spectrum at a fixed
// QUANTILE alpha = 0.999. That truncation is a convention, it is not additive
// under step subdivision, and it sets the fit's block weights.
//
// The statistically correct weight is the FISHER INFORMATION of the block,
//
//     I = INT (dp/dr)^2 / p dr ,      sigma2_eff = 1/I
//
// which is finite with no truncation anywhere (heavy-tailed ML: the estimator
// variance is 1/I even where the data variance is divergent). I is a property
// of the block's DISTRIBUTION, not of the realised residual, so it is a
// data-independent weight -- the Fisher-scoring, not the Newton, choice.
//
// HOW I IS OBTAINED
// -----------------
// By EXACT FFT INVERSION of the block characteristic function, never from a
// saddlepoint. That is not a preference: the saddlepoint Fisher route was
// measured to be 5-38 % wrong and, worse, unpredictably so (adding a channel
// carrying 5e-11 of the variance moved it by 33 %), because I is dominated by
// the delta-ray HARD EDGE which an exponential-tilting expansion places in the
// wrong location. See Documents/Resolution/NOTES_XXII_msrad.md section 4.
//
// SCOPE OF THIS IMPLEMENTATION
// ----------------------------
//   * IONIZATION channel only. On the reference pT = 3 track the multiple
//     scattering and radiative channels move 1/I of the q/p block by +0.0 %,
//     +0.05 % and +0.15 % at planes 0 / 9 / 18 (up to +1.8 % at pT = 100).
//     That is a known, bounded and one-signed omission, quantified offline,
//     NOT a claim that they are negligible in general.
//   * SCALAR: the q/p (curvilinear component 0) marginal of the block. The
//     joint 5x5 CGF is not attempted.
//
// UNITS AND THE STANDARDIZATION CONVENTION
// ----------------------------------------
// Per step the Urban record carries energies in MeV; `gs` is the full linear
// map from that step's energy loss to the block's residual coordinate,
//
//     gs = w_s * cs_s * 1e-3 / sigma        [ (q/p units) / MeV ]
//
// with w_s the transport weight of the step's noise to the end of the block
// and cs = E/p^3 [GeV^-2] the propagator's own dE -> d(q/p) map. `sigma` is a
// purely INTERNAL numerical scale that keeps the inversion grid conditioned;
// the returned inverse information is in units of sigma^2, and the caller
// multiplies it back. It is not a modelling choice and the answer must not
// depend on it (checked by scanning it).
//
// Note the transport weight is folded into `gs`, i.e. into the RESIDUAL
// coordinate, and never into the CGF parameters themselves. That is the
// convention NOTES_EXPORTS.md section 4.3 requires: it is what makes
// d eta_b/d a = 0 exactly for the alignment and B-field global parameters,
// and therefore what keeps the existing 50-mode scalar-potential production
// valid.

#include <cstddef>
#include <complex>
#include <vector>

namespace cvhcgf {

  // TRI-STATE ENVIRONMENT FLAG -- the mechanism that lets a switch have a
  // non-trivial DEFAULT and still be operable in both directions.
  //
  //     unset                                   -> `dflt`
  //     "0" / "false" / "off" / "no" / ""       -> false
  //     anything else (including "1")           -> true
  //
  // The original readers were `getenv(name) != nullptr`, i.e. presence-only.
  // That is fine for a default-OFF knob but cannot express "off" once the
  // default is ON, and a correction that cannot be turned off cannot be
  // ATTRIBUTED: the factorial comparison (all four off vs all four on) that
  // the global fit needs would be impossible to set up. Every switch below
  // therefore goes through this, whatever its default, so that one convention
  // covers both directions and `CVH_X=0` always means off.
  //
  // Case-insensitive on the false words. Note "" (set but empty) is FALSE,
  // not true: `export CVH_X=` reads as "I want it off" far more often than as
  // "I want it on".
  bool envFlag(const char *name, bool dflt);

  // THE SINGLE READER of CVH_IONONLY.
  //
  // The reference dE/dx table drops the radiative mean when this is set, so
  // the block model must supply that mean instead. The two were measured to
  // differ by EXACTLY a translation, which means each alone is a bias and the
  // pair is a no-op: at pT = 3 the reference moves by 1.03e-5 relative on q/p
  // at the outermost plane, i.e. at the Z-mass target precision.
  //
  // Having ONE function produce the value, called by both the energy-loss
  // table and the block model, is what makes the coupling structural: there is
  // no second place that could read a different answer. Do not add another
  // getenv("CVH_IONONLY") anywhere.
  bool referenceIsIonOnly();

  // THE SINGLE READER of CVH_IONI_KOKOULIN.
  //
  // Geant4's G4MuBetheBlochModel multiplies the knock-on cross section by
  // R. Kokoulin's radiative correction
  //
  //     f_K(T) = 1 + (alpha/2pi) a1 (a3 - a1),
  //     a1 = ln(1 + 2T/m_e),  a3 = ln(4 E (E - T)/M^2),
  //
  // for T above 100 keV and muons above 1 GeV. It reaches +6 % (pT = 3) to
  // +9.4 % (pT = 40) at the hard end of the spectrum.
  //
  // The MEAN already carries it: Geant4's own ComputeDEDXPerVolume /
  // CrossSectionPerVolume include it, so the dE/dx table -- and therefore the
  // reference trajectory -- is already correct (matched to 5 significant
  // digits, NOTES_SAMPLERGAP section 2b). What does NOT carry it is the
  // FLUCTUATION: the exact-delta channel of
  // G4UniversalFluctuationForExtrapolator is the tree-level PDG spectrum, so
  // both the variance it returns (the fit's Q) and the offline analytic CF
  // built from its record are missing the correction.
  //
  // Having ONE function name -- and ONE environment variable -- govern both
  // is the same structural coupling as referenceIsIonOnly above: the C++
  // variance path calls this, and the offline consumer
  // (cf_track_resolution.IONI_KOKOULIN) derives its default from the SAME
  // variable CVH_IONI_KOKOULIN. "On" therefore means on in the fit's Q and in
  // the model the closure is measured with, and there is no convention for a
  // second site to get wrong. Do not add another getenv("CVH_IONI_KOKOULIN")
  // anywhere.
  //
  // DEFAULT OFF (see the block comment on the four defaults in
  // CGFQoPBlock.cc; briefly default-ON 2026-08-16, reverted). Unset or
  // CVH_IONI_KOKOULIN=0 means not one line of the correction runs and the
  // export is bit-identical to the historical one (proven at 27/27 branches,
  // NOTES_QVALID s3). CVH_IONI_KOKOULIN=1 enables it.
  bool ioniKokoulinEnabled();

  // Number of Simpson intervals of the log-T quadrature that the fluctuation
  // model uses for the Kokoulin variance integral (CVH_IONI_KOKOULIN_NBIN,
  // default 96). A numerical-accuracy knob only -- exposed so convergence can
  // be MEASURED rather than asserted.
  int ioniKokoulinNbin();

  // THE SINGLE READER of CVH_REF_CHARGEAWARE.
  //
  // The reference trajectory's mean energy loss is charge-ODD in nature and
  // charge-EVEN in G4TablesForExtrapolatorForCVH, because that class builds
  // ONE muon table -- from G4MuonPlus -- and ONE hadron table -- from
  // G4Proton -- and then hands the muon table to both signs and scales the
  // hadron table by q*q. `muonMinus` is a member of the class, is assigned in
  // its constructor, and is never used. So the reference is not charge-blind
  // at the charge AVERAGE: it is pinned to the POSITIVE particle, and every
  // negative track carries the whole of the error.
  //
  // The physics it drops is the charge-odd part of Geant4's high-order
  // stopping-power block, G4EmCorrections::HighOrderCorrections =
  // 2*(Barkas + Bloch) + Mott. Bloch is a function of q^2 and cannot
  // contribute; Barkas (Ashley-Ritchie, ~1.29 z, FALLS with beta) and Mott
  // (Ahlen's z^3, pi*alpha*beta*z, GROWS with beta) are odd in z. At beta*gamma
  // 3.3-29.7 the odd part is 0.309-0.314 % of the restricted dE/dx, Mott
  // supplying 99.3 % of it and Barkas 0.65-0.75 %; below beta*gamma ~ 1 the two
  // exchange roles. See Documents/Resolution/NOTES_BARKAS.md.
  //
  // With this set, the tables are built for the NEGATIVE particle as well --
  // ComputeMuonDEDX(muonMinus, ...) and ComputeProtonDEDX(antiProton, ...),
  // the SAME functions with the other argument -- and ComputeDEDX,
  // ComputeRange and ComputeEnergy select between the two on the sign of the
  // track's own G4ParticleDefinition charge. No term is written out by hand
  // anywhere: the correction is whatever Geant4's own models say the negative
  // particle's stopping power is, so it is right at every beta, for every
  // singly-charged species, and it cannot drift away from Geant4 as the models
  // are revised.
  //
  // ALL THREE of ComputeDEDX / ComputeRange / ComputeEnergy must switch
  // together: EnergyAfterStep picks between `step * ComputeDEDX` (short step)
  // and `ComputeEnergy(range - step)` (long step) on linLossLimit, so a
  // partial fix would make the two branches disagree by exactly this term.
  //
  // THIS MOVES THE REFERENCE TRAJECTORY, i.e. the MEAN -- unlike
  // referenceIsIonOnly's companion knobs CVH_IONI_EXACTDELTA and
  // CVH_IONI_KOKOULIN, which move the FLUCTUATION and preserve the mean by
  // construction. Everything downstream of the reference -- the transport
  // Jacobians, the exported jacref/gradv/hesspackedv, the fitted momenta -- is
  // therefore affected. It is identically zero for a POSITIVE track, which is
  // the cheapest available regression test: any movement of a positive track
  // is a bug.
  //
  // The fluctuation model's own `meanLoss = length * dedx` reads the SAME
  // tables through G4UniversalFluctuationForExtrapolator::SetParticleAndCharge
  // and is charge-aware with it, so the reference's mean and the noise model's
  // mean cannot disagree by the charge-odd term (they did until 2026-08-16;
  // NOTES_SPECIESDEDX s8 diagnosed it, NOTES_DEFAULTON s1 closed it -- that
  // fix is a CORRECTNESS fix and is kept independently of the default).
  //
  // DEFAULT OFF (see the block comment on the four defaults in
  // CGFQoPBlock.cc; briefly default-ON 2026-08-16, reverted -- this is the
  // only charge-ODD member of the four and therefore the one degenerate with
  // the calibration's `M`). Unset or CVH_REF_CHARGEAWARE=0 means not one extra
  // table is built, no dispatch changes, and the export is bit-identical to
  // the historical one (proven at 27/27 branches, NOTES_CHARGEODD s3.2).
  // CVH_REF_CHARGEAWARE=1 enables it. Do not add another
  // getenv("CVH_REF_CHARGEAWARE") anywhere.
  bool referenceIsChargeAware();

  // THE SINGLE READER of CVH_REF_SPECIESDEDX.
  //
  // `G4TablesForExtrapolatorForCVH` builds NO hadron table except the
  // proton's. `G4EnergyLossForExtrapolatorForCVH::ComputeDEDX` serves every
  // other hadron from it at the SCALED kinetic energy
  //
  //     e = ekin * proton_mass_c2 / m ,
  //
  // and `ComputeRange` / `ComputeEnergy` do the same with the matching
  // `massratio` factors. That scaling is exact for a stopping power that is a
  // function of beta*gamma alone: it preserves beta, gamma, the mean
  // excitation energy I, the density-effect delta, the shell correction and
  // the whole Barkas/Bloch/Mott block. It does NOT preserve `Tmax`, whose
  // recoil denominator carries the PROJECTILE mass,
  //
  //     Tmax(m) = 2 m_e (b g)^2 / (1 + 2 gamma m_e/m + (m_e/m)^2) .
  //
  // The reference therefore integrates a stopping power built on the PROTON's
  // Tmax at the species' beta*gamma. Unrestricted Bethe-Bloch is
  // dE = xi [ln(2 m_e bg^2 Tmax/I^2) - 2 beta^2 - delta], so
  // d(dE)/d ln Tmax = xi EXACTLY and the whole error is
  //
  //     d(dE/dx) = xi_perlength * ln( Tmax_species / Tmax_table ) ,
  //     xi_perlength = twopi_mc2_rcl2 * n_el * z^2 / beta^2
  //
  // -- the same xi this class already forms in EnergyDispersion. Measured
  // against Geant4's own unrestricted dE/dx for the species: +5.24e-3 (pi),
  // +2.86e-4 (K), and identically zero for the proton (it IS the table) and
  // for the muon (its own table). See Documents/Resolution/NOTES_PION.md.
  //
  // WHY A LOOKUP CORRECTION AND NOT A PER-SPECIES TABLE. Unlike
  // referenceIsChargeAware -- where "call the same function with the other
  // particle" was strictly better than transcribing Barkas and Mott -- a
  // per-species table is NOT the cheap option here:
  //   * the tables live on ONE energy grid, [1 MeV, 100 TeV]. That grid is a
  //     grid in beta*gamma only after the mass scaling. A pion table on the
  //     same grid starts at bg = 0.12, below G4BetheBlochModel's validity, and
  //     the RANGE table is a cumulative integral from the bottom of the grid,
  //     so the invalid region would corrupt the range at every energy. The
  //     mass scaling is the mechanism that covers every hadron with one grid;
  //     it is not an approximation to be removed.
  //   * what is transcribed here is not a MODEL. `d(dE/dx)/d ln Tmax = xi` is
  //     the exact analytic derivative of the leading term, and Tmax is exact
  //     two-body recoil kinematics. Geant4 cannot revise either without
  //     ceasing to be Bethe-Bloch -- which is precisely NOT true of the
  //     Barkas parameterization or the Mott series.
  //   * the species set is open (the propagator builds its particle from a
  //     configured name), so a per-species table needs either a PDG list or a
  //     lazily mutated process-wide static.
  // Sub-leading mass dependence is left alone and is BELOW the table's own
  // interpolation floor: the only other mass-dependent term of
  // G4BetheBlochModel at fixed bg is the spin-1/2 recoil term
  // (0.5 Tmax/E_tot)^2, worth 5.8e-6 relative for the pion and 5.3e-7 for the
  // kaon against a spline floor of ~1e-5 (measured, NOTES_SPECIESDEDX s2).
  //
  // ALL THREE of ComputeDEDX / ComputeRange / ComputeEnergy switch together,
  // for the same reason as referenceIsChargeAware: EnergyAfterStep picks
  // between `step * ComputeDEDX` and `ComputeEnergy(range - step)` on
  // linLossLimit. The range and inverse range carry the correction as the
  // exact perturbation integral of the corrected stopping power (see
  // G4EnergyLossForExtrapolatorForCVH::speciesRangeDefect), so that
  // ComputeEnergy remains the numerical inverse of ComputeRange and the two
  // branches lose the same energy over the same step.
  // G4UniversalFluctuationForExtrapolator's own `meanLoss = length * dedx`
  // reads the same proton table and moves with it.
  //
  // THIS MOVES THE REFERENCE TRAJECTORY, i.e. the MEAN, and it is
  // charge-EVEN: unlike the Barkas/Mott term it does not cancel between charge
  // conjugates. It is identically zero for a MUON and for a PROTON, which is
  // the cheapest available regression test: any movement of either is a bug.
  //
  // DEFAULT OFF (see the block comment on the four defaults in
  // CGFQoPBlock.cc; briefly default-ON 2026-08-16, reverted). Unset or
  // CVH_REF_SPECIESDEDX=0 means not one line of the correction runs and the
  // export is bit-identical to the historical one (proven at 27/27 branches,
  // NOTES_SPECIESDEDX s3.2). CVH_REF_SPECIESDEDX=1 enables it. Do not add
  // another getenv("CVH_REF_SPECIESDEDX") anywhere.
  bool referenceIsSpeciesDedx();

  // CVH_REF_HADRAD. The hadron reference is the PROTON dE/dx table looked up at
  // e = ekin * m_p/m, which preserves beta*gamma -- correct for ionization,
  // which is a function of beta*gamma, and WRONG for radiative loss, which
  // carries explicit mass dependence. ComputeProtonDEDX builds only
  // G4BetheBlochModel, so a hadron's reference subtracts no radiative mean at
  // all while the simulation runs hBrems/hPairProd.
  //
  // Enabling this adds the radiative mean, computed per species from
  // G4hBremsstrahlungModel/G4hPairProductionModel at the particle's OWN
  // kinetic energy (no proton scaling), and -- in the propagator -- the
  // matching radiative fluctuation block.
  //
  // BOTH HALVES OR NEITHER. Measured: the missing mean and the missing
  // fluctuation cancel to ~90%, so the fluctuation alone is 0.00049 against
  // 0.00005 for the complete correction -- a 10x DEGRADATION. This single
  // switch drives both so they cannot be enabled separately.
  bool referenceHasHadronRadiative();

  // Number of Simpson intervals of the range-defect quadrature
  // (CVH_REF_SPECIESDEDX_NBIN, default 16, forced even and >= 2). A
  // numerical-accuracy knob only -- exposed, exactly as ioniKokoulinNbin is,
  // so convergence can be MEASURED rather than asserted.
  int speciesDedxNbin();

  // The dE/dx defect of a proton-table lookup, to be ADDED to the table value.
  //
  //     xi_perlength * ln( Tmax(mass) / Tmax(refMass) )   at the SAME bg
  //
  // CLHEP units throughout: `ekin`, `mass`, `refMass` in MeV, `electronDensity`
  // in 1/mm^3, the result in MeV/mm. `refMass` is the PDG mass of the particle
  // the table was actually built from (G4Proton / G4AntiProton), NOT
  // CLHEP::proton_mass_c2 -- the two differ by 7.5e-5 MeV (8.0e-8 relative),
  // which would make the correction 1.4e-11 rather than 0 for the proton and
  // destroy the exact null. Written so that mass == refMass returns +0.0 BIT
  // FOR BIT: only the recoil denominators are formed, the common
  // 2 m_e bg^2 numerator cancels in the ratio, and log(1.0) == 0.0.
  //
  // Deliberately G4-free and in this namespace so that there is exactly ONE
  // implementation, shared by the energy-loss extrapolator and by the
  // fluctuation's meanLoss. A second copy is how the two would drift apart.
  double speciesTmaxDedx(double ekin, double mass, double refMass, double charge2, double electronDensity);

  // Per-step Urban ionization record, already scaled and weighted.
  // Mirrors G4UniversalFluctuationForExtrapolator::UrbanFluctRecord with the
  // `scaling` factor applied to the energies and the transport weight folded
  // into `gs`.
  struct IoniStep {
    // 0 = Gaussian regime, 1 = Urban compound Poisson.
    //
    // regime 2/3 (the exact knock-on cross section, CVH_IONI_EXACTDELTA) is
    // REFUSED by blockExponent/blockKappa2 -- see below. It is not that the
    // struct cannot hold one: it is that in that regime `a3` means something
    // else, and there is no field for `beta2`/`etot` either.
    int regime = -1;
    double gsig2 = 0.; // regime-0 variance [MeV^2] (the ALPHA-TRUNCATED one)
    double a1 = 0., e1 = 0.;    // excitation channel 1: count, energy [MeV]
    double a2 = 0., e2 = 0.;    // excitation channel 2
    // Delta-ray collision COUNT -- in regime 1 only. The Urban record reuses
    // this slot for `xi` (an energy) in regime 2/3, which is why regime 2/3
    // must be refused rather than read.
    double a3 = 0.;
    double e0 = 0., tmax = 0.;  // 1/E^2 spectrum support [MeV]
    double gs = 0.;             // (residual units) per MeV, incl. transport
  };

  struct Config {
    // |phi| = e^{lncut} sets the end of the t grid. The contribution of the
    // truncated tail to p(z) is bounded by (1/pi) INT |phi| dt, so -60
    // (|phi| = 9e-27) is far below any relative floor used downstream.
    double lncut = -60.;
    // uniform t samples of the CF; fixes z_max = pi * nt / t_max (aliasing)
    int nt = 1 << 13;
    // zero padding; fixes dz = 2 pi / (npad * t_max) (peak resolution)
    int npad = 32;
    // RELATIVE density floor selecting the contiguous support around the mode.
    // An ABSOLUTE floor on an oversized grid was a measured 57x error in this
    // study (NOTES XV) -- do not change the sense of this.
    double floor = 1e-8;
    // Number of points of the resampled score table (see Result::psi), and
    // the half-width of the WINDOW it spans, below the mode, in standardized
    // units.
    //
    // The window is not a detail. The contiguous support of a block runs to
    // z ~ -2300 (the 1/E^2 loss tail), while all of psi's structure lives in
    // z ~ [-10, +7] and it diverges at the hard edge. Resampling npsi points
    // over the WHOLE support gives dz = 0.29 there and gets psi wrong by 38 %
    // one z unit from the edge -- measured, and it is this study's recurring
    // "match the grid to the thing being scanned" failure in a new place.
    // Spanning [zmode - psiWindow, zhi] instead gives dz ~ the native FFT
    // spacing. Below the window psi is continued by its exact 1/|z| power-law
    // asymptote rather than clamped.
    int npsi = 8192;
    double psiWindow = 60.;
    // Rigid TRANSLATION of the block density, in the same standardized units
    // as the residual. This is how the uncentred convention is carried: with
    // the radiative mean removed from the reference dE/dx table the block must
    // supply it, and the centred and uncentred radiative CFs were measured to
    // differ by EXACTLY `i t dz_rad` (6.6e-19 relative), so the whole change
    // is a translation and nothing about the shape moves. Applied by offsetting
    // the z axis, which is exact and free.
    double meanShift = 0.;
  };

  struct Result {
    bool ok = false;
    double invFisher = 0.;  // 1/I, in units of sigma^2 (i.e. of the residual)
    double kappa2 = 0.;     // block variance from the CF, same units, for ref
    double tmax = 0.;       // end of the t grid actually used
    double mass = 0.;       // probability mass on the contiguous support
    double massFrac = 0.;   // that mass / total mass on the grid
    double zmode = 0.;      // mode of the density
    double zlo = 0., zhi = 0.;  // support edges
    int nsteps = 0;

    // The SCORE psi(z) = -d ln p/dz = -p'(z)/p(z), resampled uniformly on the
    // contiguous support: psi[i] at z = psiZ0 + i*psiDz.
    //
    // It comes from the SAME two transforms that produce I -- p from phi and
    // p' from -i t phi -- so it costs nothing extra and no density (or log
    // density) is ever differenced numerically.
    //
    // IT IS DELIBERATELY NOT THE SADDLEPOINT FORM. `psi_SPA = thetahat +
    // K'''/(2 K''^2)` was measured against exact inversion at 0.52-1.44x over
    // r in [-3, +4] -- i.e. wrong by up to a factor 2 exactly where the fit
    // sits (NOTES_EXPORTS section 5.2). It also needs a saddlepoint solve per
    // evaluation, which this does not.
    double psiZ0 = 0., psiDz = 0.;
    std::vector<double> psi;
    // Diagnostic: sign changes of psi on the NATIVE grid, restricted to where
    // the density carries mass (p > 1e-4 * pmax). The restriction is
    // essential: over the full support psi = -p'/p is pure inversion noise in
    // the deep tail, where it flips sign thousands of times and any
    // unrestricted counter reports garbage (measured: 4322 crossings at
    // plane 0, all of them noise).
    //
    // The MODE is a zero of psi and there should be exactly one in the
    // mass-carrying region. The "second zero crossing" reported for the
    // SADDLEPOINT score is that density's hard-edge artefact -- K'' -> 0 makes
    // -0.5 ln(2 pi K'') diverge, so ln p_SPA turns back up -- and this counter
    // is how one finds out whether the exact score inherits it.
    int nZeroCross = 0;
    double psiWinLo = 0.;  // low edge of the score window
    int nBelowWindow = 0;  // unused; reserved for the maker's clamp counter
  };

  // psi at an arbitrary z by linear interpolation of Result::psi. psi is
  // interpolated as a FUNCTION of z; it is never inverted (it is not
  // monotonic). Outside the support the value is clamped to the nearest edge
  // and `clamped` is set -- a residual outside the support means the block
  // has essentially zero likelihood there, and the fit must not see an
  // infinity.
  double scoreAt(const Result &r, double z, bool *clamped = nullptr);

  // sine and cosine integrals of a real non-negative argument.
  // ci is returned as Cin(x) = gamma + ln x - Ci(x), which is what the
  // delta-ray term needs and which stays finite and small as x -> 0.
  void siCin(double x, double &si, double &cin);

  // <(e^{i a E} - 1 - i a E)/E^2> over the 1/E^2 spectrum on [1, w] in units
  // of e0, i.e. exactly cf_track_resolution._delta_term_2d(a, w).
  std::complex<double> deltaTerm(double a, double w);

  // Block CF exponent S(t): phi(t) = exp(S(t)). Centred (mean subtracted),
  // which is the convention the propagator's mean-loss table implies.
  //
  // THROWS std::runtime_error if any step has `regime >= 2`. That regime's
  // `a3` slot holds `xi` (an energy) rather than a collision count, so reading
  // it here is wrong by ~1e-5 in the delta channel -- and, being a WEIGHT, it
  // would not fail, it would silently change the answer. The refusal replaced
  // a one-shot warning, and is KEPT now that CVH_IONI_EXACTDELTA is
  // default-off again: a warning is not a safe guard for a silent-wrong-answer
  // hazard at any default, and the two switches can still be set together by
  // hand.
  std::complex<double> blockExponent(const std::vector<IoniStep> &steps, double t);

  // Gaussian-limit variance of the block in residual units, i.e. -S''(0).
  // Throws on regime >= 2 for the same reason as blockExponent.
  double blockKappa2(const std::vector<IoniStep> &steps);

  // t at which Re S(t) = lncut, bisected on a coarse log scan.
  double tauReach(const std::vector<IoniStep> &steps, double lncut);

  // The whole thing: 1/I by exact inversion.
  Result inverseFisher(const std::vector<IoniStep> &steps, const Config &cfg = Config());

  // in-place radix-2 FFT, forward sign convention e^{-2 pi i j k / N}
  void fftInPlace(std::vector<std::complex<double>> &a);

}  // namespace cvhcgf

#endif
