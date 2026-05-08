#ifndef MagneticField_ParametrizedEngine_ScalarPot3DEval_h
#define MagneticField_ParametrizedEngine_ScalarPot3DEval_h

/** \class ScalarPot3DEval
 *
 *  Spherical-harmonic scalar-potential basis evaluator for the CMS
 *  tracker-volume B-field.  C++ port of mfs/harmonic_basis.py.
 *
 *  Phi(r, phi, z) = Σ c_{l,m,cs} * (R/r_scale)^l * P_l^m(cos θ)
 *                                * {cos|sin}(m*phi)              (z' = z - z0)
 *  Bz   = ∂Phi/∂z      = Σ c (l+m) (R/r_scale)^{l-1}/r_scale
 *                          P_{l-1}^m(cos θ) * phi_factor
 *  Br   = ∂Phi/∂r      = Σ c (R/r_scale)^{l-1}/r_scale
 *                          [l Plm − cos θ (l+m) P_{l-1}^m] / sin θ * phi_factor
 *  Bphi = (1/r) ∂Phi/∂φ = Σ c (m/r) (R/r_scale)^l Plm
 *                            * (∓sin|cos)(m*phi)
 *
 *  Bx = Br cos φ − Bphi sin φ
 *  By = Br sin φ + Bphi cos φ
 *
 *  Plm uses the scipy.special.lpmv convention (Condon-Shortley phase
 *  (-1)^m INCLUDED — verified empirically; the scipy docstring is
 *  misleading on this point).  When cmssw_norm is true, the per-mode
 *  Schmidt factor
 *  S(l,m) = sqrt((l-m)!/(l+m)!) — the cumulative LadderUp factor of
 *  CMSSW HarmBasis3DCyl — multiplies the *whole* mode (Bz, Br, Bphi)
 *  uniformly.  See replicated-bouncing-cloud.md A.0.
 *
 *  The evaluator is constructed from a flat-text dump file produced
 *  by mfs/dump_coeffs_for_cmssw.py (Phase A.5 format).
 *
 *  \author David Walter
 */

#include <string>
#include <vector>

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

namespace magfieldparam {

class ScalarPot3DEval {
public:
  /// Construct from a Phase-A.5 flat-text dump file.
  /// Throws std::runtime_error on malformed input.
  explicit ScalarPot3DEval(const std::string& dump_path);

  /// Evaluate the absolute field at gp [cm] in Tesla, summing
  /// c_i × basis_i over the basis with the given coefficient vector.
  /// corparms.size() must equal nModes().
  GlobalVector evaluateAbsoluteAt(const GlobalPoint& gp,
                                  const std::vector<double>& corparms) const;

  /// Per-mode basis values at gp; out is resized to nModes().
  void getBzBasisAt  (const GlobalPoint& gp, std::vector<double>& out) const;
  void getBxBasisAt  (const GlobalPoint& gp, std::vector<double>& out) const;
  void getByBasisAt  (const GlobalPoint& gp, std::vector<double>& out) const;
  void getBrBasisAt  (const GlobalPoint& gp, std::vector<double>& out) const;
  void getBphiBasisAt(const GlobalPoint& gp, std::vector<double>& out) const;

  /// Compute Bz, Br, Bphi per mode in one pass (cheaper than three
  /// separate calls).  Each output vector resized to nModes().
  void evaluateBasisAt(const GlobalPoint& gp,
                       std::vector<double>& bz,
                       std::vector<double>& br,
                       std::vector<double>& bphi) const;

  /// Header info (read from the dump file).
  unsigned int nModes() const { return params_.size(); }
  unsigned int lMax() const { return l_max_; }
  double rScale() const { return r_scale_; }
  double z0() const { return z0_; }
  bool cmsswNorm() const { return cmssw_norm_; }

  /// Initial coefficients from the dump (size = nModes()).
  /// In the convention indicated by cmsswNorm().
  const std::vector<double>& initCoeffs() const { return init_coeffs_; }

  /// Mode descriptor (l, m, cs).
  struct Param {
    unsigned int l;
    unsigned int m;
    char cs;        // 'c' or 's'
  };
  const std::vector<Param>& params() const { return params_; }

private:
  // Header.
  unsigned int l_max_   = 0;
  double r_scale_       = 1.0;
  double z0_            = 0.0;
  bool   cmssw_norm_    = false;

  // Basis.
  std::vector<Param>  params_;
  std::vector<double> init_coeffs_;
  // Per-mode Schmidt factor sqrt((l-m)!/(l+m)!), or 1 if cmssw_norm_=false.
  // Pre-baked once at construction so per-evaluation cost is one multiply.
  std::vector<double> schmidt_;
};

}  // namespace magfieldparam

#endif
