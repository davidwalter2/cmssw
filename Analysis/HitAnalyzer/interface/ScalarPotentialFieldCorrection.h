#ifndef Analysis_HitAnalyzer_ScalarPotentialFieldCorrection_h
#define Analysis_HitAnalyzer_ScalarPotentialFieldCorrection_h

// Per-hit field-correction provider for the CVH global fit, backed by
// magfieldparam::ScalarPot3DEval. Replaces the per-module dBz parameter
// with a global block of spherical-harmonic coefficients of the magnetic
// scalar potential. The basis structure (l_max, mode list, Schmidt
// convention) is loaded from a Phase-A.5 dump file produced by
// mfs/dump_coeffs_for_cmssw.py.
//
// Sentinel parmset key: (parmtype = ParmTypeBfieldGlobal, DetId(modeIdx)).

#include <Eigen/Dense>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

namespace magfieldparam { class ScalarPot3DEval; }

namespace ana_hitanalyzer {

class ScalarPotentialFieldCorrection {
public:
  // Sentinel parmtype for global B-field block in (parmtype, DetId) keys.
  static constexpr int ParmTypeBfieldGlobal = 14;

  // Construct from a Phase-A.5 dump file (path produced by
  // mfs/dump_coeffs_for_cmssw.py). The file determines the basis
  // structure and supplies the initial coefficients (initCoeffs()).
  explicit ScalarPotentialFieldCorrection(const std::string& dumpPath);

  ~ScalarPotentialFieldCorrection();

  unsigned int nModes() const;

  // Initial coefficient values from the dump file, in the convention
  // declared by the dump's cmssw_norm flag. Size = nModes().
  const std::vector<double>& initCoeffs() const;

  // Register sentinel parmset entries -- call once during beginRun.
  void appendParmsetEntries(std::set<std::pair<int, DetId>>& parmset) const;

  // Cache the global index for each mode after detidparms is built.
  void resolveGlobalIndices(
      const std::map<std::pair<int, DetId>, unsigned int>& detidparms);

  // Global index of mode i in corparms / dxdparms. Valid after resolve.
  unsigned int basisGlobalIdx(unsigned int i) const { return basisGlobalIdx_[i]; }

  // Aggregate field at a hit, contracted against the current corparms
  // vector. Reads per-mode coefficients from corparms[basisGlobalIdx(i)].
  Eigen::Vector3d getCorrectionAt(const GlobalPoint& pos,
                                  const std::vector<double>& corparms) const;

  // Per-mode basis values at a hit, for use in the transport-Jacobian
  // chain rule. Each output sized to nModes().
  void getBzBasisAt(const GlobalPoint& pos, std::vector<double>& out) const;
  void getBxBasisAt(const GlobalPoint& pos, std::vector<double>& out) const;
  void getByBasisAt(const GlobalPoint& pos, std::vector<double>& out) const;

  // Full 3 x N Jacobian d(Bx, By, Bz) / d(c_i) at a hit, for the
  // alignment-Jacobian chain rule. Output sized (3, nModes()).
  void getJacobianAt(const GlobalPoint& pos,
                     Eigen::Matrix<double, 3, Eigen::Dynamic>& dBdc) const;

private:
  std::unique_ptr<magfieldparam::ScalarPot3DEval> eval_;
  std::vector<unsigned int> basisGlobalIdx_;  // global index in corparms_
};

}  // namespace ana_hitanalyzer

#endif
