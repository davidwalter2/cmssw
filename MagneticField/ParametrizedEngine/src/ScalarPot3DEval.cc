#include "MagneticField/ParametrizedEngine/interface/ScalarPot3DEval.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

// ---------------------------------------------------------------------------
// Associated Legendre polynomial P_l^m(x), scipy.special.lpmv convention
// (Condon-Shortley phase (-1)^m INCLUDED — empirically verified against
// scipy.special.lpmv despite docstring saying otherwise).  Computes the
// full table up to (l_max, l_max) and stores in plm[l*(l_max+1) + m].
//
// Recurrence (l-recurrence at fixed m):
//   P_m^m(x) = (-1)^m (2m-1)!! (1-x^2)^(m/2)    [m >= 0; 0!!:=1, (-1)!!:=1]
//   P_{m+1}^m(x) = x (2m+1) P_m^m(x)
//   (l - m + 1) P_{l+1}^m(x) = (2l+1) x P_l^m(x) - (l + m) P_{l-1}^m(x)
//
// Returns zero entries for m > l (well-defined; lpmv returns 0 in that case).
// ---------------------------------------------------------------------------
inline unsigned int plm_index(unsigned int l, unsigned int m, unsigned int lmax) {
  return l * (lmax + 1) + m;
}

void compute_plm_table(double x, unsigned int lmax, std::vector<double>& plm) {
  plm.assign((lmax + 1) * (lmax + 1), 0.0);
  if (lmax == 0) {
    plm[0] = 1.0;
    return;
  }
  const double sin_t = std::sqrt(std::max(0.0, 1.0 - x * x));

  // Diagonal: P_m^m = (-1)^m (2m-1)!! (1-x^2)^(m/2)  (CS phase included)
  plm[plm_index(0, 0, lmax)] = 1.0;
  double pmm = 1.0;
  double dfact = 1.0;  // (2m-1)!!, starts at (-1)!!=1 for m=0.
  double cs_sign = 1.0;
  for (unsigned int m = 1; m <= lmax; ++m) {
    dfact *= static_cast<double>(2 * m - 1);
    cs_sign = -cs_sign;  // (-1)^m
    pmm = cs_sign * dfact * std::pow(sin_t, static_cast<double>(m));
    plm[plm_index(m, m, lmax)] = pmm;
  }

  // First off-diagonal: P_{m+1}^m
  for (unsigned int m = 0; m + 1 <= lmax; ++m) {
    plm[plm_index(m + 1, m, lmax)] =
        x * static_cast<double>(2 * m + 1) * plm[plm_index(m, m, lmax)];
  }

  // Upward recurrence in l at fixed m
  for (unsigned int m = 0; m + 2 <= lmax; ++m) {
    for (unsigned int l = m + 2; l <= lmax; ++l) {
      const double a = (2.0 * l - 1.0) * x * plm[plm_index(l - 1, m, lmax)];
      const double b = (static_cast<double>(l) + m - 1.0)
                       * plm[plm_index(l - 2, m, lmax)];
      plm[plm_index(l, m, lmax)] = (a - b) / static_cast<double>(l - m);
    }
  }
}

// log(n!) via Stirling-tied evaluation; we need this for Schmidt factors only
// once at construction, so just use the std::lgamma path directly.
inline double schmidt_factor(unsigned int l, unsigned int m) {
  if (m == 0) return 1.0;
  // sqrt((l-m)! / (l+m)!) = exp(0.5 * (lgamma(l-m+1) - lgamma(l+m+1)))
  const double lg = 0.5 * (std::lgamma(static_cast<double>(l - m + 1))
                          - std::lgamma(static_cast<double>(l + m + 1)));
  return std::exp(lg);
}

// Trim leading whitespace.
inline std::string ltrim(const std::string& s) {
  size_t i = 0;
  while (i < s.size() && (s[i] == ' ' || s[i] == '\t')) ++i;
  return s.substr(i);
}

}  // namespace

namespace magfieldparam {

ScalarPot3DEval::ScalarPot3DEval(const std::string& dump_path) {
  std::ifstream f(dump_path);
  if (!f) {
    throw std::runtime_error("ScalarPot3DEval: cannot open " + dump_path);
  }

  bool basis_ok = false;
  unsigned int n_modes_declared = 0;
  std::string line;
  std::vector<std::string> body;
  body.reserve(512);

  while (std::getline(f, line)) {
    const std::string s = ltrim(line);
    if (s.empty()) continue;
    if (s[0] == '#') {
      std::istringstream iss(s.substr(1));
      std::string key;
      iss >> key;
      if (key == "basis_type") {
        std::string v;
        iss >> v;
        if (v != "harmonic") {
          throw std::runtime_error("ScalarPot3DEval: unsupported basis_type='"
                                   + v + "' (only 'harmonic' supported)");
        }
        basis_ok = true;
      } else if (key == "l_max") {
        int v;
        iss >> v;
        l_max_ = static_cast<unsigned int>(v);
      } else if (key == "r_scale") {
        iss >> r_scale_;
      } else if (key == "z0") {
        iss >> z0_;
      } else if (key == "cmssw_norm") {
        int v;
        iss >> v;
        cmssw_norm_ = (v != 0);
      } else if (key == "nmodes") {
        iss >> n_modes_declared;
      }
      // Other header keys (sha1_npz, column legend, ...) are ignored.
      continue;
    }
    body.push_back(s);
  }

  if (!basis_ok) {
    throw std::runtime_error("ScalarPot3DEval: missing '# basis_type harmonic' header in " + dump_path);
  }
  if (z0_ != 0.0) {
    throw std::runtime_error("ScalarPot3DEval: z0!=0 not supported (CMSSW loader assumes z0=0)");
  }
  if (n_modes_declared != 0 && body.size() != n_modes_declared) {
    throw std::runtime_error("ScalarPot3DEval: nmodes header says "
                             + std::to_string(n_modes_declared)
                             + " but found " + std::to_string(body.size()) + " body lines");
  }

  params_.reserve(body.size());
  init_coeffs_.reserve(body.size());
  schmidt_.reserve(body.size());
  for (const auto& bs : body) {
    std::istringstream iss(bs);
    int idx, l, m;
    char cs;
    double c;
    if (!(iss >> idx >> l >> m >> cs >> c)) {
      throw std::runtime_error("ScalarPot3DEval: malformed body line: '" + bs + "'");
    }
    if (cs != 'c' && cs != 's') {
      throw std::runtime_error("ScalarPot3DEval: bad cs flag '"
                               + std::string(1, cs) + "' (must be 'c' or 's')");
    }
    if (l < 0 || m < 0 || static_cast<unsigned int>(l) > l_max_
        || static_cast<unsigned int>(m) > static_cast<unsigned int>(l)) {
      throw std::runtime_error("ScalarPot3DEval: bad (l,m) pair in body: l="
                               + std::to_string(l) + " m=" + std::to_string(m));
    }
    Param p{static_cast<unsigned int>(l), static_cast<unsigned int>(m), cs};
    params_.push_back(p);
    init_coeffs_.push_back(c);
    schmidt_.push_back(cmssw_norm_ ? schmidt_factor(p.l, p.m) : 1.0);
  }
}

// ---------------------------------------------------------------------------
// Build the per-call geometry cache: trig table cos/sin(m*phi), Plm table at
// cos theta, and the Rn_pow[k] = (R/r_scale)^k power table (k=0..l_max).
// All three are O(l_max) work, so the per-mode loop pays no std::pow cost.
// ---------------------------------------------------------------------------
void ScalarPot3DEval::fillGeomCache(const GlobalPoint& gp,
                                    GeomCache& g) const {
  const double x = gp.x();
  const double y = gp.y();
  const double zshift = gp.z() - z0_;
  g.r = std::sqrt(x * x + y * y);
  g.R = std::sqrt(g.r * g.r + zshift * zshift);

  const double phi = std::atan2(y, x);
  g.cphi = std::cos(phi);
  g.sphi = std::sin(phi);

  g.cosm.resize(l_max_ + 1);
  g.sinm.resize(l_max_ + 1);
  for (unsigned int m = 0; m <= l_max_; ++m) {
    g.cosm[m] = std::cos(static_cast<double>(m) * phi);
    g.sinm[m] = std::sin(static_cast<double>(m) * phi);
  }

  if (g.R > 0.0) {
    g.cos_t = zshift / g.R;
    g.sin_t = g.r / g.R;
  } else {
    g.cos_t = (zshift >= 0.0 ? 1.0 : -1.0);
    g.sin_t = 0.0;
  }
  compute_plm_table(g.cos_t, l_max_, g.plm);

  // (R/r_scale)^k built by repeated multiply: 1 division + l_max
  // multiplies, replacing the std::pow call inside each mode's loop body.
  g.Rn_pow.resize(l_max_ + 1);
  g.Rn_pow[0] = 1.0;
  if (l_max_ >= 1) {
    const double Rn = g.R / r_scale_;
    for (unsigned int k = 1; k <= l_max_; ++k) {
      g.Rn_pow[k] = g.Rn_pow[k - 1] * Rn;
    }
  }
}

void ScalarPot3DEval::evaluateBasisFromCache(const GeomCache& g,
                                             std::vector<double>& bz,
                                             std::vector<double>& br,
                                             std::vector<double>& bphi) const {
  const unsigned int N = nModes();
  bz.assign(N, 0.0);
  br.assign(N, 0.0);
  bphi.assign(N, 0.0);

  const double inv_rscale = 1.0 / r_scale_;
  const double inv_sin_t = (g.sin_t > 1e-12) ? 1.0 / g.sin_t : 0.0;
  const double inv_r = (g.r > 1e-12) ? 1.0 / g.r : 0.0;

  for (unsigned int i = 0; i < N; ++i) {
    const Param& p = params_[i];
    const unsigned int l = p.l;
    const unsigned int m = p.m;
    const double S = schmidt_[i];

    if (l == 0) continue;  // gauge mode, all components zero

    // (R/r_scale)^(l-1) / r_scale  -- table lookup, no std::pow
    const double R_pow = g.Rn_pow[l - 1] * inv_rscale;
    const double R_l   = g.Rn_pow[l];

    const double phi_factor = (p.cs == 'c' ? g.cosm[m] : g.sinm[m]);

    // ---------- Bz: (l+m) R^{l-1}/r_scale * P_{l-1}^m * phi_factor ---------
    double plm_lm1 = 0.0;
    if (m <= l - 1) {
      plm_lm1 = g.plm[plm_index(l - 1, m, l_max_)];
    }
    bz[i] = static_cast<double>(l + m) * R_pow * plm_lm1 * phi_factor * S;

    // ---------- Br: R^{l-1}/r_scale * [l Plm - cos*(l+m)*Pl-1m] / sin -----
    double plm_l = 0.0;
    if (m <= l) {
      plm_l = g.plm[plm_index(l, m, l_max_)];
    }
    if (g.sin_t > 1e-12) {
      const double numerator = static_cast<double>(l) * plm_l
                             - g.cos_t * static_cast<double>(l + m) * plm_lm1;
      br[i] = R_pow * numerator * inv_sin_t * phi_factor * S;
    }  // else br[i] stays 0 (L'Hopital on axis)

    // ---------- Bphi: (m/r) (R/r_scale)^l * P_l^m * (∓sin|+cos)(m phi) ----
    if (m > 0 && g.r > 1e-12) {
      const double phi_dphi = (p.cs == 'c' ? -g.sinm[m] : g.cosm[m]);
      bphi[i] = static_cast<double>(m) * inv_r * R_l * plm_l * phi_dphi * S;
    }  // else bphi[i] stays 0
  }
}

// ---------------------------------------------------------------------------
// Per-mode basis evaluation (public): fills bz, br, bphi vectors at gp.
// ---------------------------------------------------------------------------
void ScalarPot3DEval::evaluateBasisAt(const GlobalPoint& gp,
                                      std::vector<double>& bz,
                                      std::vector<double>& br,
                                      std::vector<double>& bphi) const {
  GeomCache g;
  fillGeomCache(gp, g);
  evaluateBasisFromCache(g, bz, br, bphi);
}

// ---------------------------------------------------------------------------
// Single-component basis getters (delegate to evaluateBasisAt + projection).
// ---------------------------------------------------------------------------
void ScalarPot3DEval::getBzBasisAt(const GlobalPoint& gp,
                                   std::vector<double>& out) const {
  std::vector<double> br, bphi;
  evaluateBasisAt(gp, out, br, bphi);
}

void ScalarPot3DEval::getBrBasisAt(const GlobalPoint& gp,
                                   std::vector<double>& out) const {
  std::vector<double> bz, bphi;
  evaluateBasisAt(gp, bz, out, bphi);
}

void ScalarPot3DEval::getBphiBasisAt(const GlobalPoint& gp,
                                     std::vector<double>& out) const {
  std::vector<double> bz, br;
  evaluateBasisAt(gp, bz, br, out);
}

void ScalarPot3DEval::getBxBasisAt(const GlobalPoint& gp,
                                   std::vector<double>& out) const {
  std::vector<double> bz, br, bphi;
  evaluateBasisAt(gp, bz, br, bphi);
  const double phi = std::atan2(gp.y(), gp.x());
  const double cphi = std::cos(phi);
  const double sphi = std::sin(phi);
  const unsigned int N = nModes();
  out.resize(N);
  for (unsigned int i = 0; i < N; ++i) {
    out[i] = br[i] * cphi - bphi[i] * sphi;
  }
}

void ScalarPot3DEval::getByBasisAt(const GlobalPoint& gp,
                                   std::vector<double>& out) const {
  std::vector<double> bz, br, bphi;
  evaluateBasisAt(gp, bz, br, bphi);
  const double phi = std::atan2(gp.y(), gp.x());
  const double cphi = std::cos(phi);
  const double sphi = std::sin(phi);
  const unsigned int N = nModes();
  out.resize(N);
  for (unsigned int i = 0; i < N; ++i) {
    out[i] = br[i] * sphi + bphi[i] * cphi;
  }
}

// ---------------------------------------------------------------------------
// Aggregate field at gp = Σ c_i × basis_i, returned as GlobalVector (Tesla).
// Hot path during Geant4 stepping — accumulates in three scalars rather
// than allocating per-mode bz/br/bphi vectors and summing them in a
// second pass. Geometry tables (Plm, Rn_pow, trig) are built once via
// fillGeomCache; the per-mode body is one fused multiply-add per
// component.
// ---------------------------------------------------------------------------
GlobalVector ScalarPot3DEval::evaluateAbsoluteAt(
    const GlobalPoint& gp, const std::vector<double>& corparms) const {
  const unsigned int N = nModes();
  if (corparms.size() != N) {
    throw std::runtime_error("ScalarPot3DEval::evaluateAbsoluteAt: corparms size "
                             + std::to_string(corparms.size())
                             + " != nModes() " + std::to_string(N));
  }

  GeomCache g;
  fillGeomCache(gp, g);

  const double inv_rscale = 1.0 / r_scale_;
  const double inv_sin_t  = (g.sin_t > 1e-12) ? 1.0 / g.sin_t : 0.0;
  const double inv_r      = (g.r     > 1e-12) ? 1.0 / g.r     : 0.0;

  double sumBz = 0.0, sumBr = 0.0, sumBphi = 0.0;
  for (unsigned int i = 0; i < N; ++i) {
    const Param& p = params_[i];
    const unsigned int l = p.l;
    if (l == 0) continue;
    const unsigned int m = p.m;

    // Coefficient absorbs both the user value and the Schmidt factor in
    // one multiply (saves one per mode).
    const double c = corparms[i] * schmidt_[i];
    if (c == 0.0) continue;

    const double R_pow = g.Rn_pow[l - 1] * inv_rscale;  // R^{l-1}/r_scale
    const double R_l   = g.Rn_pow[l];                   // R^l (for Bphi)

    const double phi_factor = (p.cs == 'c' ? g.cosm[m] : g.sinm[m]);

    const double plm_lm1 = (m <= l - 1) ? g.plm[plm_index(l - 1, m, l_max_)] : 0.0;
    const double plm_l   = (m <= l)     ? g.plm[plm_index(l,     m, l_max_)] : 0.0;

    // Bz
    sumBz += c * static_cast<double>(l + m) * R_pow * plm_lm1 * phi_factor;

    // Br (zero on axis by L'Hopital; numerator vanishes there too)
    if (g.sin_t > 1e-12) {
      const double numerator = static_cast<double>(l) * plm_l
                             - g.cos_t * static_cast<double>(l + m) * plm_lm1;
      sumBr += c * R_pow * numerator * inv_sin_t * phi_factor;
    }

    // Bphi (zero on axis: Plm ~ sin^m, m>=1 -> 0; m=0 -> dPhi/dphi=0)
    if (m > 0 && g.r > 1e-12) {
      const double phi_dphi = (p.cs == 'c' ? -g.sinm[m] : g.cosm[m]);
      sumBphi += c * static_cast<double>(m) * inv_r * R_l * plm_l * phi_dphi;
    }
  }

  const double Bx = sumBr * g.cphi - sumBphi * g.sphi;
  const double By = sumBr * g.sphi + sumBphi * g.cphi;
  return GlobalVector(Bx, By, sumBz);
}

}  // namespace magfieldparam
