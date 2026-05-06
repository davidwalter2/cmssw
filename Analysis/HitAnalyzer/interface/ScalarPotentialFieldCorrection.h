#ifndef Analysis_HitAnalyzer_ScalarPotentialFieldCorrection_h
#define Analysis_HitAnalyzer_ScalarPotentialFieldCorrection_h

// Wraps magfieldparam::HarmBasis3DCyl as a per-hit field-correction provider
// for the CVH global fit. Replaces the per-module dBz parameter with a global
// block of ~50 spherical-harmonic coefficients of the magnetic scalar potential.
//
// Modes: all (L,M) with L in [1, lmaxFull], plus optional extras at given
// (L, |M|). Each (L, M=0) contributes 1 mode; each (L, |M|>=1) contributes 2
// modes (cos and sin azimuthal phases).
//
// Sentinel parmset key: (parmtype = ParmTypeBfieldGlobal, DetId(modeIdx)) so
// the modes piggyback on the existing detidparms infrastructure.

#include <Eigen/Dense>
#include <map>
#include <memory>
#include <set>
#include <utility>
#include <vector>

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "MagneticField/ParametrizedEngine/interface/HarmBasis3DCyl.h"

namespace ana_hitanalyzer {

class ScalarPotentialFieldCorrection {
public:
  // Sentinel parmtype for global B-field block in (parmtype, DetId) keys.
  // Distinct from per-module types 0-7 and resolution types 8-11.
  static constexpr int ParmTypeBfieldGlobal = 14;

  // lmaxFull: include all (L,M) modes with L in [1, lmaxFull].
  // extras: list of (L, |M|) pairs to include in addition. For |M|=0 one mode
  //   is added; for |M|>=1 two modes (cos and sin) are added.
  ScalarPotentialFieldCorrection(unsigned int lmaxFull,
                                 const std::vector<std::pair<int, int>> &extras);

  unsigned int nModes() const { return nModes_; }

  // Register sentinel parmset entries — call once during beginRun.
  void appendParmsetEntries(std::set<std::pair<int, DetId>> &parmset) const;

  // Cache the global index for each mode after detidparms is built.
  void resolveGlobalIndices(
      const std::map<std::pair<int, DetId>, unsigned int> &detidparms);

  // Global index of mode i in corparms / dxdparms. Valid after resolve.
  unsigned int basisGlobalIdx(unsigned int i) const { return basisGlobalIdx_[i]; }

  // Field correction (Bx, By, Bz) [Tesla] at a hit, contracted against the
  // current corparms vector. Reads the per-mode coefficients from
  // corparms[basisGlobalIdx(i)].
  Eigen::Vector3d getCorrectionAt(const GlobalPoint &pos,
                                  const std::vector<double> &corparms) const;

  // Per-mode Bz basis value at a hit, for use in the transport-Jacobian
  // chain rule: J_col(mode i) = J_col(dBz) * bzPerMode[i]. Output size nModes_.
  void getBzBasisAt(const GlobalPoint &pos,
                    std::vector<double> &bzPerMode) const;

  // Full 3 x N Jacobian d(Bx, By, Bz) / d(c_i) at a hit, for use in the
  // alignment-Jacobian chain rule. Each column is one basis mode evaluated at
  // the position, rotated to Cartesian. Output is sized (3, nModes_).
  void getJacobianAt(const GlobalPoint &pos,
                     Eigen::Matrix<double, 3, Eigen::Dynamic> &dBdc) const;

private:
  // Build the active flat-index list (into the underlying HarmBasis3DCyl) from
  // (lmaxFull, extras).
  void buildActiveIndices(unsigned int lmaxFull,
                          const std::vector<std::pair<int, int>> &extras);

  unsigned int nModes_;
  std::vector<unsigned int> activeIdx_;          // flat indices into basis_
  std::vector<unsigned int> basisGlobalIdx_;     // global index in corparms_
  mutable magfieldparam::HarmBasis3DCyl basis_;  // mutable for SetPoint/EvalXxx
};

}  // namespace ana_hitanalyzer

#endif
