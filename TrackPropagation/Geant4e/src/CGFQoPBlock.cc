#include "TrackPropagation/Geant4e/interface/CGFQoPBlock.h"

#include "FWCore/Utilities/interface/Exception.h"

#include <CLHEP/Units/PhysicalConstants.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>

namespace cvhcgf {

  namespace {
    Switches g_switches;
    bool g_configured = false;
  }

  void configure(const edm::ParameterSet &pset) {
    Switches s;
    s.ioniExactDelta = pset.getParameter<bool>("IoniExactDelta");
    s.ioniKokoulin = pset.getParameter<bool>("IoniKokoulin");
    s.referenceChargeAware = pset.getParameter<bool>("ReferenceChargeAware");
    s.referenceSpeciesDedx = pset.getParameter<bool>("ReferenceSpeciesDedx");
    s.referenceHadRad = pset.getParameter<bool>("ReferenceHadronRadiative");
    s.referenceIonOnly = pset.getParameter<bool>("ReferenceIonizationOnly");
    s.ioniUrban2021 = pset.getParameter<bool>("IoniUrban2021");
    s.dumpEmParameters = pset.getParameter<bool>("DumpEmParameters");
    s.emHarmonise = pset.getParameter<bool>("EmHarmonise");
    s.dedxScale = pset.getParameter<double>("DedxScale");
    s.speciesDedxNbin = pset.getParameter<int>("ReferenceSpeciesDedxNbin");
    s.ioniKokoulinNbin = pset.getParameter<int>("IoniKokoulinNbin");
    s.ioniKokoulinCgfNbin = pset.getParameter<int>("IoniKokoulinCgfNbin");
    s.cgfRadiative = pset.getParameter<bool>("CgfRadiativeChannel");
    s.cgfQoPMode = pset.getParameter<int>("CgfQoPMode");
    s.cgfQoPRefresh = pset.getParameter<int>("CgfQoPRefresh");
    if (s.cgfQoPMode < 0 || s.cgfQoPMode > 3) {
      throw cms::Exception("Configuration")
          << "CgfQoPMode must be 0 (legacy Gaussian weight, diagnostic only), 1 (the Fisher weight), "
             "2 (1 + per-leg print) or 3 (1 + IRLS re-centring); got "
          << s.cgfQoPMode;
    }
    if (s.cgfQoPRefresh < 0) {
      throw cms::Exception("Configuration")
          << "CgfQoPRefresh must be >= 0 (0 = freeze after the first sweep); got " << s.cgfQoPRefresh;
    }
    s.ioniExactDeltaT0 = pset.getParameter<double>("IoniExactDeltaT0");

    for (auto const &nb : {std::make_pair("ReferenceSpeciesDedxNbin", s.speciesDedxNbin),
                           std::make_pair("IoniKokoulinNbin", s.ioniKokoulinNbin)}) {
      if (nb.second < 2 || (nb.second % 2) != 0) {
        throw cms::Exception("Configuration")
            << nb.first << " must be even and >= 2 (Simpson needs an even interval count); got "
            << nb.second;
      }
    }

    if (s.ioniKokoulinCgfNbin < 0) {
      throw cms::Exception("Configuration")
          << "IoniKokoulinCgfNbin must be >= 0 (0 = omit the term in the block CGF, N = N log-T "
             "buckets); got "
          << s.ioniKokoulinCgfNbin;
    }

    if (g_configured) {
      // The readers below configure PROCESS-GLOBAL Geant4 model classes, so two
      // propagators asking for different physics cannot both be served.
      // Honouring whichever was constructed first would be a silent, ordering-
      // dependent answer -- exactly the failure the move off getenv removes.
      const bool same = s.ioniExactDelta == g_switches.ioniExactDelta &&
                        s.ioniKokoulin == g_switches.ioniKokoulin &&
                        s.referenceChargeAware == g_switches.referenceChargeAware &&
                        s.referenceSpeciesDedx == g_switches.referenceSpeciesDedx &&
                        s.referenceHadRad == g_switches.referenceHadRad &&
                        s.referenceIonOnly == g_switches.referenceIonOnly &&
                        s.ioniUrban2021 == g_switches.ioniUrban2021 &&
                        s.dumpEmParameters == g_switches.dumpEmParameters &&
                        s.emHarmonise == g_switches.emHarmonise &&
                        s.dedxScale == g_switches.dedxScale &&
                        s.speciesDedxNbin == g_switches.speciesDedxNbin &&
                        s.ioniKokoulinNbin == g_switches.ioniKokoulinNbin &&
                        s.ioniExactDeltaT0 == g_switches.ioniExactDeltaT0 &&
                        s.cgfRadiative == g_switches.cgfRadiative &&
                        s.cgfQoPMode == g_switches.cgfQoPMode &&
                        s.cgfQoPRefresh == g_switches.cgfQoPRefresh &&
                        s.ioniKokoulinCgfNbin == g_switches.ioniKokoulinCgfNbin;
      if (!same) {
        throw cms::Exception("Configuration")
            << "cvhcgf::configure called twice with different values. The CVH energy-loss "
            << "switches configure process-global Geant4 model classes, so every "
            << "GeantPropagatorESProducer in a job must declare the same ones.";
      }
      return;
    }
    g_switches = s;
    g_configured = true;
  }

  const Switches &switches() {
    if (!g_configured) {
      throw cms::Exception("Configuration")
          << "cvhcgf::switches() used before cvhcgf::configure(). The CVH energy-loss switches "
          << "are set from the Geant4ePropagator ESProducer's ParameterSet; a job that reaches "
          << "this code without one is misconfigured. (These were environment variables before "
          << "and silently defaulted, which is the behaviour this replaces.)";
    }
    return g_switches;
  }

  bool referenceIsIonOnly() {
    return switches().referenceIonOnly;
  }

  // The energy-loss corrections are DEFAULT-ON, so that the MODEL represents
  // what the SIMULATION actually does.  Each of them exists because the model
  // was missing something Geant4 runs: the exact knock-on cross section, the
  // Kokoulin term, a charge-aware mean loss, a species-specific dE/dx, and
  // hadron radiative loss (the sim runs hBrems/hPairProd while a hadron's
  // reference subtracted no radiative mean at all).  Leaving them off makes the
  // default configuration a model of a simulation nobody runs, and makes every
  // unqualified number a measurement of a known omission.
  //
  // History: flipped ON 2026-08-16 (NOTES_DEFAULTON), reverted to OFF the same
  // week by 4b24be8, and restored here.
  //
  // THE ATTRIBUTION GATE THAT 4b24be8 NAMED IS NOT CLOSED BY THIS COMMIT, and
  // whoever consumes these Jacobians needs to know it: the global fit has still
  // not been run both ways, and `CVH_REF_CHARGEAWARE` is the one charge-ODD
  // member -- degenerate in a single fit with the calibration's `M`, which is
  // read as physical misalignment.  Default-ON changes which way that risk
  // points (a real charge-odd loss now enters the reference instead of being
  // absorbed into M) but does not remove it.  Run the fit both ways before
  // quoting M.
  //
  // `CVH_IONI_URBAN2021` deliberately stays OFF: NOTES_DELTASPEC s10.4
  // falsified the prediction it was built for, so it is a diagnostic, not a
  // correction, and turning it on would NOT move the model toward the sim.
  bool ioniKokoulinEnabled() {
    return switches().ioniKokoulin;
  }

  bool referenceIsChargeAware() {
    return switches().referenceChargeAware;
  }

  bool referenceIsSpeciesDedx() {
    return switches().referenceSpeciesDedx;
  }

  bool referenceHasHadronRadiative() {
    return switches().referenceHadRad;
  }

  int speciesDedxNbin() {
    // Validated in `configure` (Simpson needs an even, positive interval
    // count), where a bad value is a configuration ERROR rather than being
    // silently replaced by the default the way the getenv reader did.
    return switches().speciesDedxNbin;
  }

  double speciesTmaxDedx(double ekin, double mass, double refMass, double charge2, double electronDensity) {
    // Guard rather than compute nonsense: every one of these is a programming
    // error upstream, and a silent 0 here is a silently uncorrected reference.
    if (!(ekin > 0.) || !(mass > 0.) || !(refMass > 0.) || !(electronDensity > 0.) || !(charge2 > 0.)) {
      return 0.;
    }
    const double tau = ekin / mass;
    const double gam = tau + 1.0;
    const double bg2 = tau * (tau + 2.0);
    const double beta2 = bg2 / (gam * gam);
    if (!(beta2 > 0.)) {
      return 0.;
    }
    // Tmax(m) = 2 m_e bg^2 / den(m). The numerator is common to both masses at
    // fixed bg (and gamma = sqrt(1+bg^2) is a function of bg alone), so the
    // ratio is den(refMass)/den(mass) and the numerator never has to be formed.
    // This is what makes mass == refMass give exactly log(1.0) = +0.0.
    const double r = CLHEP::electron_mass_c2 / mass;
    const double rref = CLHEP::electron_mass_c2 / refMass;
    const double den = 1.0 + 2.0 * gam * r + r * r;
    const double denref = 1.0 + 2.0 * gam * rref + rref * rref;
    const double xi = CLHEP::twopi_mc2_rcl2 * electronDensity * charge2 / beta2;
    return xi * std::log(denref / den);
  }

  int ioniKokoulinNbin() {
    return switches().ioniKokoulinNbin;
  }

  int ioniKokoulinCgfNbin() {
    return switches().ioniKokoulinCgfNbin;
  }

  bool cgfRadiativeEnabled() {
    return switches().cgfRadiative;
  }

  int cgfQoPMode() {
    return switches().cgfQoPMode;
  }

  int cgfQoPRefresh() {
    return switches().cgfQoPRefresh;
  }

  namespace {
    constexpr double kEuler = 0.5772156649015328606;
    constexpr double kPi = 3.14159265358979323846;
    // Series/closed-form crossover of the delta-ray term, and the series
    // length. Identical to cf_track_resolution._DT_SER / _DT_NSER so the two
    // implementations select the same branch at the same argument.
    constexpr double kDtSer = 2.0;
    constexpr int kDtNser = 80;
  }  // namespace

  //--------------------------------------------------------------------------
  // Si(x) and Cin(x) = gamma + ln x - Ci(x), real x >= 0.
  //
  // Numerical Recipes `cisi`: power series below x = 2, Lentz continued
  // fraction for E1(ix) above. Cin rather than Ci is returned because the
  // delta-ray term needs Cin -- which is O(x^2/4) at small x, where Ci is
  // logarithmically divergent and would have to be cancelled by hand.
  //--------------------------------------------------------------------------
  void siCin(double x, double &si, double &cin) {
    constexpr double kEps = 1e-17;
    constexpr int kMaxIt = 200;
    constexpr double kFpMin = 1e-300;
    constexpr double kTMin = 2.0;

    const double t = std::fabs(x);
    if (t == 0.) {
      si = 0.;
      cin = 0.;
      return;
    }
    if (t < 1e-150) {
      // leading terms; below this the k = 2 term underflows and NR's own
      // err = term/|sum| goes 0/0
      si = t;
      cin = 0.25 * t * t;
      if (x < 0.)
        si = -si;
      return;
    }
    if (t <= kTMin) {
      // power series: Si = sum (-1)^k x^{2k+1}/((2k+1)(2k+1)!)
      //               Cin = sum (-1)^{k+1} x^{2k}/(2k (2k)!)   (k >= 1)
      double sum = 0., sums = 0., sumc = 0.;
      double sign = 1., fact = 1.;
      bool odd = true;
      int k = 1;
      for (; k <= kMaxIt; ++k) {
        fact *= t / k;
        const double term = fact / k;
        sum += sign * term;
        const double err = term / std::fabs(sum);
        if (odd) {
          sign = -sign;
          sums = sum;
          sum = sumc;
        } else {
          sumc = sum;
          sum = sums;
        }
        if (err < kEps)
          break;
        odd = !odd;
      }
      si = sums;
      cin = -sumc;  // the alternating sum above accumulates -Cin
    } else {
      // continued fraction for E1(i t): exp(i t) * (sum) = ...
      std::complex<double> b(1.0, t), c(1.0 / kFpMin, 0.0), d, h;
      d = h = 1.0 / b;
      int i = 2;
      for (; i <= kMaxIt; ++i) {
        const double a = -static_cast<double>(i - 1) * static_cast<double>(i - 1);
        b += 2.0;
        d = 1.0 / (a * d + b);
        c = b + a / c;
        const std::complex<double> del = c * d;
        h *= del;
        if (std::fabs(del.real() - 1.0) + std::fabs(del.imag()) < kEps)
          break;
      }
      h *= std::complex<double>(std::cos(t), -std::sin(t));
      const double ci = -h.real();
      si = kPi / 2. + h.imag();
      cin = kEuler + std::log(t) - ci;
    }
    if (x < 0.)
      si = -si;
  }

  //--------------------------------------------------------------------------
  // The delta-ray term.
  //
  //   <(e^{iaE} - 1 - iaE)/E^2> * N = expm1(y) - expm1(x)/w + y (F(y) - F(x))
  //   y = i a,  x = i a w,  N = 1 - 1/w,  F(s) = -sum_{k>=1} s^k/(k k!)
  //
  // F on the imaginary axis is elementary: F(iu) = Ein(-iu) = Cin(u) - i Si(u)
  // for u > 0, Cin even and Si odd. That removes any need for a complex
  // exponential integral, and both pieces are then evaluated in the branch
  // where they do not cancel.
  //
  // THE EXPANSION PARAMETER IS a*w, NOT a. The support reaches E = w, and for
  // muons w = tmax/e0 ~ 1e8-1e10, so guarding the series on |a| alone selects
  // it in a regime where it is wrong by orders of magnitude (a documented
  // failure: it turned the model CF into a pure oscillation at small t).
  //--------------------------------------------------------------------------
  std::complex<double> deltaTerm(double a, double w) {
    const double N = 1.0 - 1.0 / w;
    const double aw = a * w;

    if (std::fabs(aw) <= kDtSer) {
      // sum_{k>=2} (x^k/w - y^k) / (k! (k-1)),  x = i a w, y = i a
      const std::complex<double> x(0., aw), y(0., a);
      std::complex<double> tx(1., 0.), ty(1., 0.), acc(0., 0.);
      for (int k = 1; k < kDtNser; ++k) {
        tx *= x / static_cast<double>(k);
        ty *= y / static_cast<double>(k);
        if (k >= 2) {
          acc += (tx / w - ty) / static_cast<double>(k - 1);
          if (std::abs(tx) < 1e-19 * std::max(std::abs(acc), 1e-300))
            break;
        }
      }
      return acc / N;
    }

    // closed form. expm1(y) - expm1(x)/w cancels the 1/w analytically; the
    // grouped form e^{ia} - (1 - 1/w) loses relative precision linearly in w
    // (7e-5 at w = 1e9), which matters because the result is only O(a^2 w).
    const double sa = std::sin(a), ca = std::cos(a);
    const double saw = std::sin(aw), caw = std::cos(aw);
    const std::complex<double> em1y(ca - 1.0, sa);
    const std::complex<double> em1x(caw - 1.0, saw);

    double siA, cinA, siW, cinW;
    siCin(std::fabs(a), siA, cinA);
    siCin(std::fabs(aw), siW, cinW);
    // y (F(y) - F(x)) with F(iu) = Cin(|u|) - i sgn(u) Si(|u|), y = ia, x = iaw.
    //
    // The sgn is already IN the expression rather than carried as a variable,
    // which is why the real part has a fabs and the imaginary part does not.
    // With w > 0 -- always, w is a ratio of energies and N = 1 - 1/w -- both
    // arguments share s = sgn(a), so
    //
    //   y (F(y) - F(x)) = ia [ (CinA - CinW) - i s (SiA - SiW) ]
    //                   = a s (SiA - SiW)  +  i a (CinA - CinW)
    //
    // and a s = a sgn(a) = |a|. (A `sgn` local was computed here and never
    // used; it was this factor, already folded in.)
    const std::complex<double> ydF(std::fabs(a) * (siA - siW), a * (cinA - cinW));

    return (em1y - em1x / w + ydF) / N;
  }

  //--------------------------------------------------------------------------
  // THE EXACT KNOCK-ON CHANNEL (regime 2/3).
  //
  // With CVH_IONI_EXACTDELTA the Urban 1/E^2 delta channel is replaced by the
  // cross section Geant4's own G4BetheBlochModel samples,
  //
  //     dN/dT = (xi / T^2) [ 1 - beta^2 T/Tmax + (spin 1/2) T^2/(2 E^2) ] ,
  //
  // on [e0, tmax], with `xi` -- the record's `a3` slot in this regime -- the
  // step's Landau energy scale in MeV. Its centred log-CF is
  //
  //     S = (xi/e0) [ J0 - (beta^2 e0/tmax) J1 ]
  //         + (spin 1/2) (xi/e0) (e0^2/2E^2) J2 ,
  //
  //     J0 = INT_1^w (e^{yu} - 1 - yu)/u^2 du       y = i a,  a = gs e0 t
  //     J1 = INT_1^w (e^{yu} - 1 - yu)/u   du       w = tmax/e0
  //     J2 = INT_1^w (e^{yu} - 1 - yu)     du
  //
  // -- the same three integrals the tree-level 1/E^2 term already needs (it is
  // J0/N), so the only new mathematics is J1 and J2.
  //
  // UNTIL 2026-08-20 THIS THREW. The refusal was correct while nothing
  // implemented the channel -- reading `a3` as a collision count when it holds
  // xi is a ~1e-5 error in a WEIGHT, i.e. one that changes the answer without
  // failing -- but it made the in-fit CGF and the exact-delta correction
  // mutually exclusive, and the exact delta is default-ON (e232c20) precisely
  // because it is what the simulation does. A weight that cannot be evaluated
  // on the physics the fit runs is not a weight.
  //
  // This is a port of `cf_track_resolution.exact_delta_exponent` and
  // `_delta_terms_exact`, the offline routines every published closure number
  // in NOTES_DELTASPEC rests on, and it is validated against them step for
  // step (testCGFQoPBlock + cgf_cxx_validate.py).
  //
  // `deltaTerm` above is deliberately NOT refactored to share this code. It
  // could be -- J0/N is exactly its return value -- but the series branches
  // differ in their termination test, so sharing would perturb the regime-1
  // path in its last bits, and bit-identity of what already ships is worth
  // more than the fifteen lines.
  //--------------------------------------------------------------------------
  void deltaTermsExact(double a, double w, std::complex<double> &J0, std::complex<double> &J1,
                       std::complex<double> &J2) {
    const std::complex<double> y(0., a), x(0., a * w);

    // THE EXPANSION PARAMETER IS a*w, NOT a -- same trap as `deltaTerm`: the
    // support reaches u = w ~ 1e8, so a series guarded on |a| is selected in a
    // regime where it is wrong by orders of magnitude.
    if (std::fabs(a * w) <= kDtSer) {
      std::complex<double> tx(1., 0.), ty(1., 0.), s0(0., 0.), s1(0., 0.), s2(0., 0.);
      for (int k = 1; k < kDtNser; ++k) {
        tx *= x / static_cast<double>(k);
        ty *= y / static_cast<double>(k);
        if (k >= 2) {
          s0 += (tx / w - ty) / static_cast<double>(k - 1);
          s1 += (tx - ty) / static_cast<double>(k);
          s2 += (tx * w - ty) / static_cast<double>(k + 1);
          // J2 is the largest of the three by a factor w, so it sets the
          // termination; the |tx| * w on the left is the same quantity.
          if (std::abs(tx) * w < 1e-19 * std::max(std::abs(s2), 1e-300))
            break;
        }
      }
      J0 = s0;
      J1 = s1;
      J2 = s2;
      return;
    }

    // Closed form. Ein(-i u) = Cin(|u|) - i sgn(u) Si(|u|), Cin even, Si odd,
    // so the exponential integral never needs a complex argument.
    double siA, cinA, siW, cinW;
    siCin(std::fabs(a), siA, cinA);
    siCin(std::fabs(a * w), siW, cinW);
    const double sgn = (a >= 0.) ? 1. : -1.;
    const std::complex<double> ein(cinA - cinW, -sgn * (siA - siW));

    const double sa = std::sin(a), ca = std::cos(a);
    const double saw = std::sin(a * w), caw = std::cos(a * w);
    const std::complex<double> em1y(ca - 1.0, sa);
    const std::complex<double> em1x(caw - 1.0, saw);
    const std::complex<double> ey(ca, sa);
    const std::complex<double> ex(caw, saw);

    J0 = em1y - em1x / w + y * ein;
    J1 = ein - y * (w - 1.0);
    J2 = (ex - ey) / y - (w - 1.0) - y * (w * w - 1.0) / 2.0;
  }

  // One step's exact-delta exponent. `at` = gs * t, i.e. the conjugate
  // variable in 1/MeV BEFORE the e0 scaling that `deltaTermsExact` works in --
  // which is what lets the Kokoulin bucketing below call this on a sub-range
  // [lo, hi] with the same `at`.
  std::complex<double> exactDeltaExponent(
      double xi, double e0, double tmax, double beta2, double etot, double at, bool spinHalf) {
    // `xi` is only required to be NONZERO, not positive: the Kokoulin
    // bucketing calls this with xi * (f_K - 1), and while Geant4's f_K is >= 1
    // everywhere it is used, a guard that silently drops negative weights would
    // turn a sign error into a small answer instead of a visible one.
    if (!(e0 > 0.) || !(tmax > e0) || xi == 0.)
      return std::complex<double>(0., 0.);
    const double w = tmax / e0;
    std::complex<double> J0, J1, J2;
    deltaTermsExact(at * e0, w, J0, J1, J2);
    const double pref = xi / e0;
    const double c1 = beta2 * e0 / tmax;
    std::complex<double> S = pref * (J0 - c1 * J1);
    if (spinHalf && etot > 0.)
      S += pref * (e0 * e0 / (2.0 * etot * etot)) * J2;
    return S;
  }

  //--------------------------------------------------------------------------
  // KOKOULIN, in the CGF.
  //
  // G4MuBetheBlochModel multiplies the knock-on cross section by
  // f_K(T) = 1 + (alpha/2pi) a1 (a3 - a1) above T = 100 keV for muons above
  // 1 GeV, and `G4UniversalFluctuationForExtrapolator` already puts it into
  // the VARIANCE the legacy weight uses. The block CGF has to carry it too, or
  // switching the fit from the variance to the Fisher information would
  // silently DROP a correction that is +6 % at the hard end -- exactly where
  // the Fisher information of this block lives.
  //
  // f_K - 1 is smooth in ln T and e^{iaT} is not, so f_K - 1 is made piecewise
  // constant on a log grid and each bin's oscillatory integral is done in
  // CLOSED FORM by `exactDeltaExponent` on the sub-range. The substitution
  // beta^2 -> beta^2 * hi/tmax is what turns that routine's "upper limit == the
  // beta^2 denominator" convention into a genuine sub-range integral whose
  // beta^2 term still carries the KINEMATIC Tmax; summing contiguous bins
  // reproduces the full-range closed form (7e-15 offline).
  //
  // The muon guard is Geant4's, mirrored rather than approximated: the factor
  // exists in G4MuBetheBlochModel and nowhere else, so a kaon's regime-3 record
  // must not receive it. The mass comes from the record by the identity
  // m = E sqrt(1 - beta^2) -- no PDG column needed, and no assumption.
  //--------------------------------------------------------------------------
  constexpr double kKokAlphaPrime = 1.0 / (2.0 * 3.14159265358979323846 * 137.035999084);
  constexpr double kKokTMin = 0.1;      // G4MuBetheBlochModel::limitKinEnergy, MeV
  constexpr double kKokMuMin = 1000.0;  // G4MuBetheBlochModel::lowestKinEnergy, MeV
  constexpr double kMuMass = 105.6583745;
  constexpr double kKokMassTol = 1.0;   // separates mu from its neighbour pi (34 MeV away)

  double kokoulinFactor(double T, double etot) {
    if (!(T > kKokTMin) || !(T < etot - kMuMass))
      return 1.0;
    const double a1 = std::log(1.0 + 2.0 * T / CLHEP::electron_mass_c2);
    const double a3 = std::log(4.0 * etot * (etot - T) / (kMuMass * kMuMass));
    return 1.0 + kKokAlphaPrime * a1 * (a3 - a1);
  }

  // Antiderivative of the exact-delta integrand, for the bucketed sum below:
  //
  //   F(T) = B0(T) - (beta2/tmax) B1(T) + inv2E2 B2(T),
  //   B0 = INT (e^{icT}-1-icT)/T^2 dT = -G/T - i c Ein(-icT)
  //   B1 = INT (e^{icT}-1-icT)/T   dT = -Ein(-icT) - i c T
  //   B2 = INT (e^{icT}-1-icT)     dT = expm1(icT)/(ic) - T - i c T^2/2
  //
  // so a sub-range [lo, hi] costs F(hi) - F(lo) and CONSECUTIVE BUCKETS SHARE
  // AN EDGE: n+1 evaluations instead of 2n. That is the whole point -- the
  // Kokoulin term is evaluated at every point of the inversion grid, so its
  // cost is what decides whether the correction is affordable at all.
  //
  // Below |cT| = 2 the closed forms cancel (B2's expm1/(ic) is 1/(ic) plus the
  // terms that are subtracted off), so the series are used instead; they are
  // the same expansions with the cancelling orders removed analytically:
  //   B1 = sum_{k>=2} u^k/(k k!),        u = i c T
  //   B2 = (1/(ic)) sum_{k>=3} u^k/k!
  //   B0 = (1/T) sum_{j>=2} u^j [1/((j-1)(j-1)!) - 1/j!]
  std::complex<double> exactDeltaAnti(double c, double T, double beta2, double tmax, double inv2E2) {
    if (c == 0. || T <= 0.)
      return std::complex<double>(0., 0.);
    const std::complex<double> ic(0., c);
    const double u = c * T;
    std::complex<double> B0, B1, B2;

    if (std::fabs(u) <= kDtSer) {
      const std::complex<double> uu(0., u);
      std::complex<double> tk(1., 0.);  // u^k / k!  built up
      std::complex<double> s0(0., 0.), s1(0., 0.), s2(0., 0.);
      double fact = 1.;
      for (int k = 1; k < kDtNser; ++k) {
        tk *= uu / static_cast<double>(k);  // = u^k / k!
        fact *= k;                          // = k!
        if (k >= 2) {
          // u^k/(k k!)
          s1 += tk / static_cast<double>(k);
          // u^j [1/((j-1)(j-1)!) - 1/j!] = u^j/j! * [j!/((j-1)(j-1)!) - 1]
          const double jf = fact / (static_cast<double>(k - 1) * (fact / static_cast<double>(k)));
          s0 += tk * (jf - 1.);
        }
        if (k >= 3)
          s2 += tk;
        if (std::abs(tk) < 1e-19 * std::max(std::abs(s2), 1e-300))
          break;
      }
      B0 = s0 / T;
      B1 = s1;
      B2 = s2 / ic;
    } else {
      double si, cin;
      siCin(std::fabs(u), si, cin);
      const double sgn = (u >= 0.) ? 1. : -1.;
      const std::complex<double> ein(cin, -sgn * si);  // Ein(-i u)
      const double su = std::sin(u), cu = std::cos(u);
      const std::complex<double> em1(cu - 1.0, su);    // e^{iu} - 1
      const std::complex<double> G = em1 - std::complex<double>(0., u);
      B0 = -G / T - ic * ein;
      B1 = -ein - std::complex<double>(0., u);
      B2 = em1 / ic - T - ic * (T * T) * 0.5;
    }
    return B0 - (beta2 / tmax) * B1 + inv2E2 * B2;
  }

  std::complex<double> kokoulinExponent(const IoniStep &s, double at, int nbin) {
    const double m = s.etot * std::sqrt(std::max(1.0 - s.beta2, 0.0));
    if (std::fabs(m - kMuMass) >= kKokMassTol)  // not a muon: Geant4 does not correct it
      return std::complex<double>(0., 0.);
    if (!(s.etot - kMuMass > kKokMuMin))        // below G4MuIonisation's lowestKinEnergy
      return std::complex<double>(0., 0.);
    const double lo0 = std::max(s.e0, kKokTMin);
    if (!(lo0 < s.tmax))
      return std::complex<double>(0., 0.);
    const bool spinHalf = (s.regime == 2);
    // Edges by `pow(ratio, i/nbin)` rather than by repeated multiplication:
    // the offline reference builds them that way (`lo0 * (tmax/lo0)**fr`), and
    // the two differ in the last bits otherwise, which is exactly the sort of
    // difference a like-for-like validation should not have to absorb.
    const double ratio = s.tmax / lo0;
    const double inv2E2 = (spinHalf && s.etot > 0.) ? 1.0 / (2.0 * s.etot * s.etot) : 0.;
    // One antiderivative per EDGE, reused by the two buckets that share it.
    std::complex<double> Flo = exactDeltaAnti(at, lo0, s.beta2, s.tmax, inv2E2);
    std::complex<double> S(0., 0.);
    for (int i = 0; i < nbin; ++i) {
      const double lo = lo0 * std::pow(ratio, static_cast<double>(i) / nbin);
      const double hi = lo0 * std::pow(ratio, static_cast<double>(i + 1) / nbin);
      const std::complex<double> Fhi = exactDeltaAnti(at, hi, s.beta2, s.tmax, inv2E2);
      const double kap = kokoulinFactor(std::sqrt(lo * hi), s.etot) - 1.0;
      if (kap != 0.)
        S += (s.a3 * kap) * (Fhi - Flo);
      Flo = Fhi;
    }
    return S;
  }

  // The exact channel's second cumulant, INT T^2 dN/dT dT, in closed form:
  //   xi [ (tmax - e0) - beta^2 (tmax^2 - e0^2)/(2 tmax)
  //        + (spin 1/2) (tmax^3 - e0^3)/(6 E^2) ] .
  // This is the same integrand `kokoulinVarIntegral` reweights in the
  // fluctuation model, which is what keeps kappa2 here and the variance there
  // the same object.
  double exactDeltaKappa2(const IoniStep &s) {
    if (!(s.e0 > 0.) || !(s.tmax > s.e0) || !(s.a3 > 0.))
      return 0.;
    const double t0 = s.e0, t1 = s.tmax;
    double m2 = (t1 - t0) - s.beta2 * (t1 * t1 - t0 * t0) / (2.0 * t1);
    if (s.regime == 2 && s.etot > 0.)
      m2 += (t1 * t1 * t1 - t0 * t0 * t0) / (6.0 * s.etot * s.etot);
    return s.a3 * m2;
  }

  // The Kokoulin excess of the same second moment. Composite Simpson in ln T,
  // mirroring `G4UniversalFluctuationForExtrapolator::kokoulinVarIntegral`
  // term for term -- the two must not drift apart, because one is the weight
  // the fit used to use and the other is the weight it uses now.
  double kokoulinKappa2(const IoniStep &s, int nbin) {
    const double m = s.etot * std::sqrt(std::max(1.0 - s.beta2, 0.0));
    if (std::fabs(m - kMuMass) >= kKokMassTol || !(s.etot - kMuMass > kKokMuMin))
      return 0.;
    const double lo = std::max(s.e0, kKokTMin);
    if (!(s.tmax > lo))
      return 0.;
    const double du = std::log(s.tmax / lo) / nbin;
    const double inv2E2 = (s.regime == 2 && s.etot > 0.) ? 1.0 / (2.0 * s.etot * s.etot) : 0.;
    double acc = 0.;
    for (int i = 0; i <= nbin; ++i) {
      const double T = lo * std::exp(i * du);
      const double wgt = 1.0 - s.beta2 * T / s.tmax + inv2E2 * T * T;
      const double f = (kokoulinFactor(T, s.etot) - 1.0) * wgt * T;  // dT = T du
      const double c = (i == 0 || i == nbin) ? 1. : ((i % 2) ? 4. : 2.);
      acc += c * f;
    }
    return s.a3 * acc * du / 3.0;
  }

  //--------------------------------------------------------------------------
  // THE RADIATIVE CHANNEL.
  //
  // Bremsstrahlung and pair production are a compound Poisson exactly as
  // ionization is, so the exponent simply adds
  //
  //     S_rad(t) = sum_steps INT dv (dN/dv) (e^{i a v E} - 1 - i a v E),
  //     a = gs * t ,
  //
  // with the same centring `-1 - i a eps` that the ionization channels use.
  // The centring is not a convention here, it is a REQUIREMENT: the reference
  // trajectory already subtracts the radiative mean (the dE/dx table is built
  // with ionOnly = false), so the block must describe the fluctuation about
  // that mean and nothing else. `makeRadSpectrum` normalizes the spectrum to
  // the same mean the table subtracted, which is what makes S'(0) = 0 exact
  // rather than approximate.
  //
  // The integral is a trapezoid on the shared v grid -- the same rule the
  // offline reference uses, so the two agree to their common discretization
  // rather than to a difference of two quadratures.
  //--------------------------------------------------------------------------
  void makeRadSpectrum(const double *v,
                       const double *shapeBrem,
                       const double *shapePair,
                       double dEBrem,
                       double dEPair,
                       double etot,
                       int nv,
                       double *out) {
    for (int i = 0; i < nv; ++i)
      out[i] = 0.;
    if (nv < 2 || v == nullptr || out == nullptr)
      return;
    for (int proc = 0; proc < 2; ++proc) {
      const double *shape = proc == 0 ? shapeBrem : shapePair;
      const double dE = proc == 0 ? dEBrem : dEPair;
      if (shape == nullptr || !(dE > 0.))
        continue;
      // norm = INT v E shape dv, the mean this shape would produce at unit
      // normalization; trapezoid on the same grid as the exponent's integral.
      double norm = 0.;
      for (int i = 0; i + 1 < nv; ++i) {
        const double dv = v[i + 1] - v[i];
        norm += 0.5 * dv * (v[i] * shape[i] + v[i + 1] * shape[i + 1]);
      }
      norm *= etot;
      if (!(norm > 0.))
        continue;
      const double f = dE / norm;
      for (int i = 0; i < nv; ++i)
        out[i] += shape[i] * f;
    }
  }

  namespace {
    // Trapezoid weights of the shared v grid, formed once per call rather than
    // per step: the integral is sum_i w_i (dN/dv)_i f(v_i), and folding the
    // rule into `w` turns each step into one dot product.
    void trapWeights(const double *v, int nv, std::vector<double> &w) {
      w.assign(nv, 0.);
      for (int i = 0; i + 1 < nv; ++i) {
        const double dv = 0.5 * (v[i + 1] - v[i]);
        w[i] += dv;
        w[i + 1] += dv;
      }
    }
  }  // namespace

  std::complex<double> blockExponent(const Block &blk, double t) {
    double sre = 0., sim = 0.;
    for (const IoniStep &s : blk.ioni) {
      if (s.regime == 0) {
        sre += -0.5 * t * t * s.gsig2 * s.gs * s.gs;
        continue;
      }
      // excitations: Poisson at FIXED energy, centred
      if (s.a1 > 0. && s.e1 > 0.) {
        const double th = s.gs * s.e1 * t;
        sre += s.a1 * (std::cos(th) - 1.);
        sim += s.a1 * (std::sin(th) - th);
      }
      if (s.a2 > 0. && s.e2 > 0.) {
        const double th = s.gs * s.e2 * t;
        sre += s.a2 * (std::cos(th) - 1.);
        sim += s.a2 * (std::sin(th) - th);
      }
      // delta rays. Regime 1: `a3` collisions from 1/E^2 on [e0, tmax].
      // Regime 2/3: the exact knock-on cross section normalized by xi = `a3`,
      // plus Geant4's Kokoulin radiative correction when that is on.
      if (s.a3 > 0. && s.tmax > s.e0 && s.e0 > 0.) {
        if (s.regime >= 2) {
          std::complex<double> d = exactDeltaExponent(
              s.a3, s.e0, s.tmax, s.beta2, s.etot, s.gs * t, s.regime == 2);
          if (s.kokNbin > 0)
            d += kokoulinExponent(s, s.gs * t, s.kokNbin);
          sre += d.real();
          sim += d.imag();
        } else {
          const std::complex<double> d = deltaTerm(s.gs * s.e0 * t, s.tmax / s.e0);
          sre += s.a3 * d.real();
          sim += s.a3 * d.imag();
        }
      }
    }

    // The radiative channel, on the shared v grid.
    if (!blk.rad.empty() && blk.v != nullptr && blk.nv > 1) {
      std::vector<double> w;
      trapWeights(blk.v, blk.nv, w);
      for (const RadStep &r : blk.rad) {
        if (r.dNdv == nullptr || !(r.etot > 0.) || r.gs == 0.)
          continue;
        const double a = r.gs * t;
        for (int i = 0; i < blk.nv; ++i) {
          const double q = w[i] * r.dNdv[i];
          if (q == 0.)
            continue;
          const double x = a * blk.v[i] * r.etot;
          // e^{ix} - 1 - ix for REAL x: no complex transcendental is needed,
          // and `cos x - 1` is written as -2 sin^2(x/2), the cancellation-free
          // form. Below |x| = 1e-4 the leading terms are used instead, for the
          // same reason -- and the guard is on x itself, never on `a`, which
          // is the trap the delta-ray term documents.
          if (std::fabs(x) < 1e-4) {
            sre += q * (-0.5 * x * x);
            sim += q * (x * x * x / 6.);
          } else {
            const double sh = std::sin(0.5 * x);
            sre += q * (-2.0 * sh * sh);
            sim += q * (std::sin(x) - x);
          }
        }
      }
    }
    return std::complex<double>(sre, sim);
  }

  //--------------------------------------------------------------------------
  double blockKappa2(const Block &blk) {
    // kappa2 = sum over channels of (count) * <E^2> * gs^2. For the 1/E^2
    // spectrum on [e0, tmax], <E^2>/N = (tmax - e0) * e0 * tmax / (tmax - e0)
    // ... written out: the normalized density is (e0 tmax/(tmax-e0)) / E^2,
    // so <E^2> = e0 tmax (tmax - e0)/(tmax - e0) = e0 * tmax.
    // IONIZATION ONLY, DELIBERATELY. The radiative second moment is dominated
    // by v -> 1 -- it is the catastrophic radiator, not the block -- so
    // including it here would hand the standardization scale and the
    // Gaussian-shape carrier `dQ2u` to the rarest events on the track. The
    // radiative channel belongs in the CF (blockExponent) and in the Fisher
    // information built from it, and that is where it is.
    double k2 = 0.;
    for (const IoniStep &s : blk.ioni) {
      if (s.regime == 0) {
        k2 += s.gsig2 * s.gs * s.gs;
        continue;
      }
      const double g2 = s.gs * s.gs;
      if (s.a1 > 0. && s.e1 > 0.)
        k2 += s.a1 * s.e1 * s.e1 * g2;
      if (s.a2 > 0. && s.e2 > 0.)
        k2 += s.a2 * s.e2 * s.e2 * g2;
      if (s.a3 > 0. && s.tmax > s.e0 && s.e0 > 0.) {
        if (s.regime >= 2) {
          double m2 = exactDeltaKappa2(s);
          if (s.kokNbin > 0)
            m2 += kokoulinKappa2(s, s.kokNbin);
          k2 += m2 * g2;
        } else {
          k2 += s.a3 * s.e0 * s.tmax * g2;
        }
      }
    }
    return k2;
  }

  //--------------------------------------------------------------------------
  double tauReach(const Block &blk, double lncut) {
    constexpr double lo = 1e-3, hi = 1e6;
    constexpr int n = 400;
    const double r = std::pow(hi / lo, 1.0 / (n - 1));
    double tprev = lo;
    double t = lo;
    for (int i = 0; i < n; ++i, t = lo * std::pow(r, i)) {
      if (blockExponent(blk, t).real() < lncut) {
        if (i == 0)
          return lo;
        double a = tprev, b = t;
        for (int k = 0; k < 24; ++k) {
          const double m = std::sqrt(a * b);
          if (blockExponent(blk, m).real() < lncut)
            b = m;
          else
            a = m;
        }
        return std::sqrt(a * b);
      }
      tprev = t;
    }
    return hi;
  }

  //--------------------------------------------------------------------------
  void fftInPlace(std::vector<std::complex<double>> &a) {
    const size_t n = a.size();
    if (n < 2)
      return;
    // bit reversal
    for (size_t i = 1, j = 0; i < n; ++i) {
      size_t bit = n >> 1;
      for (; j & bit; bit >>= 1)
        j ^= bit;
      j ^= bit;
      if (i < j)
        std::swap(a[i], a[j]);
    }
    for (size_t len = 2; len <= n; len <<= 1) {
      const double ang = -2. * kPi / static_cast<double>(len);
      const std::complex<double> wl(std::cos(ang), std::sin(ang));
      for (size_t i = 0; i < n; i += len) {
        std::complex<double> w(1., 0.);
        for (size_t k = 0; k < len / 2; ++k) {
          const std::complex<double> u = a[i + k];
          const std::complex<double> v = a[i + k + len / 2] * w;
          a[i + k] = u + v;
          a[i + k + len / 2] = u - v;
          w *= wl;
        }
      }
    }
  }

  //--------------------------------------------------------------------------
  Result inverseFisher(const Block &blk, const Config &cfg) {
    Result out;
    out.nsteps = static_cast<int>(blk.ioni.size());
    if (blk.ioni.empty() && blk.rad.empty())
      return out;

    out.kappa2 = blockKappa2(blk);
    const double tmax = tauReach(blk, cfg.lncut) * 1.3;
    out.tmax = tmax;
    if (!(tmax > 0.) || !std::isfinite(tmax))
      return out;

    const int nt = cfg.nt;
    const size_t N = static_cast<size_t>(nt) * static_cast<size_t>(cfg.npad);
    const double dt = tmax / static_cast<double>(nt);

    // phi(t_j) * trapezoid weight, zero padded. Two transforms with the same
    // phi: multiplier 1 gives p, multiplier -i t gives dp/dz. Differencing a
    // density (or its log) numerically never enters.
    std::vector<std::complex<double>> cp(N, std::complex<double>(0., 0.));
    std::vector<std::complex<double>> cd(N, std::complex<double>(0., 0.));
    for (int j = 0; j < nt; ++j) {
      const double t = dt * j;
      const std::complex<double> S = blockExponent(blk, t);
      if (S.real() < -745.)  // exp underflows; the tail is already negligible
        continue;
      const std::complex<double> phi = std::exp(S);
      const double wgt = (j == 0) ? 0.5 * dt : dt;
      cp[j] = phi * wgt;
      cd[j] = phi * wgt * std::complex<double>(0., -t);
    }
    fftInPlace(cp);
    fftInPlace(cd);

    // z_k = 2 pi k /(N dt) for k < N/2, 2 pi (k-N)/(N dt) above: the FFT
    // output is already in "sorted" order once the second half is read as
    // negative z, so walk it as [N/2 .. N-1] then [0 .. N/2-1].
    const double dz = 2. * kPi / (static_cast<double>(N) * dt);
    const size_t half = N / 2;
    // cfg.meanShift translates the density rigidly; applying it to the z axis
    // is exact and costs nothing (see Config::meanShift).
    auto zOf = [&](size_t idx) {
      const size_t k = (idx < half) ? (idx + half) : (idx - half);
      const double kk = (k < half) ? static_cast<double>(k) : (static_cast<double>(k) - static_cast<double>(N));
      return kk * dz + cfg.meanShift;
    };
    auto pOf = [&](size_t idx) {
      const size_t k = (idx < half) ? (idx + half) : (idx - half);
      return cp[k].real() / kPi;
    };
    auto dpOf = [&](size_t idx) {
      const size_t k = (idx < half) ? (idx + half) : (idx - half);
      return cd[k].real() / kPi;
    };

    // contiguous support around the mode, RELATIVE floor
    size_t imode = 0;
    double pmax = -std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < N; ++i) {
      const double p = pOf(i);
      if (p > pmax) {
        pmax = p;
        imode = i;
      }
    }
    if (!(pmax > 0.))
      return out;
    const double thr = cfg.floor * pmax;
    size_t ilo = imode, ihi = imode;
    while (ilo > 0 && pOf(ilo - 1) > thr)
      --ilo;
    while (ihi + 1 < N && pOf(ihi + 1) > thr)
      ++ihi;
    if (ihi - ilo < 10)
      return out;

    // I = INT (p')^2/p, mass = INT p, both by trapezoid on the support
    double I = 0., mass = 0., tot = 0.;
    for (size_t i = ilo; i < ihi; ++i) {
      const double p0 = pOf(i), p1 = pOf(i + 1);
      const double d0 = dpOf(i), d1 = dpOf(i + 1);
      const double f0 = (p0 > 0.) ? d0 * d0 / p0 : 0.;
      const double f1 = (p1 > 0.) ? d1 * d1 / p1 : 0.;
      I += 0.5 * dz * (f0 + f1);
      mass += 0.5 * dz * (p0 + p1);
    }
    for (size_t i = 0; i + 1 < N; ++i)
      tot += 0.5 * dz * (std::max(pOf(i), 0.) + std::max(pOf(i + 1), 0.));

    if (!(I > 0.) || !std::isfinite(I) || !(mass > 0.))
      return out;

    out.invFisher = mass / I;  // I is normalized by the mass on the support
    out.mass = mass;
    out.massFrac = (tot > 0.) ? mass / tot : 0.;
    out.zmode = zOf(imode);
    out.zlo = zOf(ilo);
    out.zhi = zOf(ihi);

    // Score table. psi = -p'/p on the native grid, resampled uniformly onto
    // npsi points spanning the support. Both p and p' come from transforms of
    // the same phi, so nothing is differenced.
    const int npsi = std::max(16, cfg.npsi);
    const double pfl = cfg.floor * pmax;
    auto psiOf = [&](size_t idx) {
      const double p = pOf(idx);
      return (p > pfl) ? -dpOf(idx) / p : 0.;
    };
    out.psi.assign(npsi, 0.);
    out.psiZ0 = std::max(out.zlo, out.zmode - cfg.psiWindow);
    out.psiWinLo = out.psiZ0;
    out.psiDz = (out.zhi - out.psiZ0) / static_cast<double>(npsi - 1);
    {
      size_t j = ilo;
      for (int i = 0; i < npsi; ++i) {
        const double zt = out.psiZ0 + i * out.psiDz;
        while (j + 1 < ihi && zOf(j + 1) < zt)
          ++j;
        const double za = zOf(j), zb = zOf(j + 1);
        const double fa = psiOf(j), fb = psiOf(j + 1);
        const double u = (zb > za) ? (zt - za) / (zb - za) : 0.;
        out.psi[i] = fa + u * (fb - fa);
      }
    }
    // sign changes where the density carries mass -- see Result::nZeroCross
    {
      const double sig = 1e-4 * pmax;
      int prev = 0;
      for (size_t i = ilo; i <= ihi; ++i) {
        if (pOf(i) <= sig)
          continue;
        const double v = psiOf(i);
        const int s = (v < 0.) ? -1 : ((v > 0.) ? 1 : 0);
        if (s != 0 && prev != 0 && s != prev)
          ++out.nZeroCross;
        if (s != 0)
          prev = s;
      }
    }

    out.ok = true;
    return out;
  }

  //--------------------------------------------------------------------------
  double scoreAt(const Result &r, double z, bool *clamped) {
    if (clamped)
      *clamped = false;
    if (r.psi.empty())
      return 0.;
    const double x = (z - r.psiZ0) / r.psiDz;
    if (x <= 0.) {
      // Below the window the density is the 1/E^2 loss tail, where
      // p ~ C/z^2 exactly, hence psi = -p'/p -> -2/z. Continue by that law
      // anchored on the window edge rather than clamping: clamping would
      // OVERSTATE |psi| far out, and psi falling off like 1/|z| is exactly
      // the property that forces the forward-expansion safeguard.
      if (clamped)
        *clamped = true;
      const double z0 = r.psiZ0;
      if (z < 0. && z0 < 0.)
        return r.psi.front() * (z0 / z);
      return r.psi.front();
    }
    const size_t n = r.psi.size();
    if (x >= static_cast<double>(n - 1)) {
      if (clamped)
        *clamped = true;
      return r.psi.back();
    }
    const size_t i = static_cast<size_t>(x);
    const double u = x - static_cast<double>(i);
    return r.psi[i] + u * (r.psi[i + 1] - r.psi[i]);
  }

}  // namespace cvhcgf
