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
    s.ioniExactDeltaT0 = pset.getParameter<double>("IoniExactDeltaT0");

    for (auto const &nb : {std::make_pair("ReferenceSpeciesDedxNbin", s.speciesDedxNbin),
                           std::make_pair("IoniKokoulinNbin", s.ioniKokoulinNbin)}) {
      if (nb.second < 2 || (nb.second % 2) != 0) {
        throw cms::Exception("Configuration")
            << nb.first << " must be even and >= 2 (Simpson needs an even interval count); got "
            << nb.second;
      }
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
                        s.ioniExactDeltaT0 == g_switches.ioniExactDeltaT0;
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
    const double sgn = (a >= 0. ? 1. : -1.);
    // y (F(y) - F(x)) with F(iu) = Cin(|u|) - i sgn(u) Si(|u|)
    const std::complex<double> ydF(std::fabs(a) * (siA - siW), a * (cinA - cinW));

    return (em1y - em1x / w + ydF) / N;
  }

  //--------------------------------------------------------------------------
  // REGIME 2/3 IS REFUSED, NOT MISREAD.
  //
  // blockExponent/blockKappa2 below branch on `regime == 0` and treat
  // everything else as regime 1, i.e. they read `a3` as a COLLISION COUNT. In
  // regime 2/3 (CVH_IONI_EXACTDELTA) that slot holds `xi`, an ENERGY -- ~0.072
  // MeV where the count is ~7.6 -- so the delta channel comes out wrong by
  // ~1e-5. Being a WEIGHT, it does not fail: it silently changes the answer.
  //
  // That was tolerable while CVH_IONI_EXACTDELTA was default-off. It is not
  // tolerable now that it is default-ON: anyone enabling the CGF prototype
  // would inherit the wrong answer without doing anything. The block CGF has
  // no exact-delta channel (and no Kokoulin term either), so the only correct
  // behaviours are "implement it" or "refuse". Until the CGF workstream ports
  // cf_track_resolution.exact_delta_exponent -- which also needs `beta2` and
  // `etot`, which IoniStep does not carry -- this refuses.
  //
  // Throwing std::runtime_error rather than cms::Exception keeps this
  // translation unit free of framework headers, as it has always been; the
  // in-fit call site (Geant4ePropagator) converts to a cms::Exception before
  // any record reaches here.
  namespace {
    void refuseExactDelta(const std::vector<IoniStep> &steps) {
      for (const IoniStep &s : steps) {
        if (s.regime >= 2) {
          throw std::runtime_error(
              "cvhcgf: the block CGF has no regime-2/3 (exact-delta) channel. The record's `a3` slot holds xi "
              "(an energy), not a collision count, so evaluating it here would be wrong by ~1e-5 in the delta "
              "channel WITHOUT failing. Run with CVH_IONI_EXACTDELTA=0, or give blockExponent/blockKappa2 the "
              "exact knock-on cross section (cf_track_resolution.exact_delta_exponent).");
        }
      }
    }
  }  // namespace

  std::complex<double> blockExponent(const std::vector<IoniStep> &steps, double t) {
    refuseExactDelta(steps);
    double sre = 0., sim = 0.;
    for (const IoniStep &s : steps) {
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
      // delta rays: 1/E^2 compound Poisson on [e0, tmax]
      if (s.a3 > 0. && s.tmax > s.e0 && s.e0 > 0.) {
        const std::complex<double> d = deltaTerm(s.gs * s.e0 * t, s.tmax / s.e0);
        sre += s.a3 * d.real();
        sim += s.a3 * d.imag();
      }
    }
    return std::complex<double>(sre, sim);
  }

  //--------------------------------------------------------------------------
  double blockKappa2(const std::vector<IoniStep> &steps) {
    // kappa2 = sum over channels of (count) * <E^2> * gs^2. For the 1/E^2
    // spectrum on [e0, tmax], <E^2>/N = (tmax - e0) * e0 * tmax / (tmax - e0)
    // ... written out: the normalized density is (e0 tmax/(tmax-e0)) / E^2,
    // so <E^2> = e0 tmax (tmax - e0)/(tmax - e0) = e0 * tmax.
    refuseExactDelta(steps);  // see the note above blockExponent
    double k2 = 0.;
    for (const IoniStep &s : steps) {
      if (s.regime == 0) {
        k2 += s.gsig2 * s.gs * s.gs;
        continue;
      }
      const double g2 = s.gs * s.gs;
      if (s.a1 > 0. && s.e1 > 0.)
        k2 += s.a1 * s.e1 * s.e1 * g2;
      if (s.a2 > 0. && s.e2 > 0.)
        k2 += s.a2 * s.e2 * s.e2 * g2;
      if (s.a3 > 0. && s.tmax > s.e0 && s.e0 > 0.)
        k2 += s.a3 * s.e0 * s.tmax * g2;
    }
    return k2;
  }

  //--------------------------------------------------------------------------
  double tauReach(const std::vector<IoniStep> &steps, double lncut) {
    constexpr double lo = 1e-3, hi = 1e6;
    constexpr int n = 400;
    const double r = std::pow(hi / lo, 1.0 / (n - 1));
    double tprev = lo;
    double t = lo;
    for (int i = 0; i < n; ++i, t = lo * std::pow(r, i)) {
      if (blockExponent(steps, t).real() < lncut) {
        if (i == 0)
          return lo;
        double a = tprev, b = t;
        for (int k = 0; k < 24; ++k) {
          const double m = std::sqrt(a * b);
          if (blockExponent(steps, m).real() < lncut)
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
  Result inverseFisher(const std::vector<IoniStep> &steps, const Config &cfg) {
    Result out;
    out.nsteps = static_cast<int>(steps.size());
    if (steps.empty())
      return out;

    out.kappa2 = blockKappa2(steps);
    const double tmax = tauReach(steps, cfg.lncut) * 1.3;
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
      const std::complex<double> S = blockExponent(steps, t);
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
