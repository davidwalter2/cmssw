#include "Analysis/HitAnalyzer/interface/ScalarPotentialFieldCorrection.h"

#include <cmath>
#include <stdexcept>

namespace ana_hitanalyzer {

namespace {

// HarmBasis3DCyl constructor builds modes for L=1..N in flat order:
// for each L, k advances by 2L+1: M=0, then M=+1, M=-1, M=+2, M=-2, ...
// Returns the flat index of (L, M=0) in the basis.
unsigned int firstIdxForL(unsigned int L) {
  // sum_{l=1}^{L-1} (2l+1) = (L-1)(L+1) = L^2 - 1
  return L * L - 1;
}

// Required basis dimension to cover lmaxFull plus any extras.
unsigned int requiredBasisDim(unsigned int lmaxFull,
                              const std::vector<std::pair<int, int>> &extras) {
  unsigned int N = lmaxFull;
  for (const auto &lm : extras) {
    if (static_cast<unsigned int>(lm.first) > N)
      N = static_cast<unsigned int>(lm.first);
  }
  if (N < 1) N = 1;
  return N;
}

}  // namespace

ScalarPotentialFieldCorrection::ScalarPotentialFieldCorrection(
    unsigned int lmaxFull, const std::vector<std::pair<int, int>> &extras)
    : basis_(requiredBasisDim(lmaxFull, extras)) {
  buildActiveIndices(lmaxFull, extras);
}

void ScalarPotentialFieldCorrection::buildActiveIndices(
    unsigned int lmaxFull, const std::vector<std::pair<int, int>> &extras) {
  activeIdx_.clear();

  // Block 1: all (L, M) with L in [1, lmaxFull].
  // The HarmBasis3DCyl flat ordering is contiguous for L=1..lmaxFull,
  // covering the first lmaxFull*(lmaxFull+2) flat indices.
  const unsigned int nFull = lmaxFull * (lmaxFull + 2);
  activeIdx_.reserve(nFull + 2 * extras.size());
  for (unsigned int k = 0; k < nFull; ++k) {
    activeIdx_.push_back(k);
  }

  // Block 2: extras. Each (L, |M|) with L > lmaxFull adds 1 mode if |M|=0,
  // or 2 modes (cos and sin phases) if |M|>=1.
  for (const auto &lm : extras) {
    const int L = lm.first;
    const int absM = std::abs(lm.second);
    if (L < 1 || static_cast<unsigned int>(L) > basis_.GetDim()) {
      throw std::out_of_range("ScalarPotentialFieldCorrection: extra L out of basis range");
    }
    if (static_cast<unsigned int>(L) <= lmaxFull) {
      // already in block 1
      continue;
    }
    if (absM > L) {
      throw std::out_of_range("ScalarPotentialFieldCorrection: extra |M| > L");
    }
    const unsigned int kBase = firstIdxForL(static_cast<unsigned int>(L));
    if (absM == 0) {
      activeIdx_.push_back(kBase);
    } else {
      // M=+|absM| at kBase + 2|absM| - 1; M=-|absM| at kBase + 2|absM|
      activeIdx_.push_back(kBase + 2 * absM - 1);
      activeIdx_.push_back(kBase + 2 * absM);
    }
  }

  nModes_ = activeIdx_.size();
  basisGlobalIdx_.assign(nModes_, 0u);  // resolved later
}

void ScalarPotentialFieldCorrection::appendParmsetEntries(
    std::set<std::pair<int, DetId>> &parmset) const {
  for (unsigned int i = 0; i < nModes_; ++i) {
    parmset.emplace(ParmTypeBfieldGlobal, DetId(i));
  }
}

void ScalarPotentialFieldCorrection::resolveGlobalIndices(
    const std::map<std::pair<int, DetId>, unsigned int> &detidparms) {
  for (unsigned int i = 0; i < nModes_; ++i) {
    auto it = detidparms.find(std::make_pair(ParmTypeBfieldGlobal, DetId(i)));
    if (it == detidparms.end()) {
      throw std::runtime_error(
          "ScalarPotentialFieldCorrection: mode global index not in detidparms");
    }
    basisGlobalIdx_[i] = it->second;
  }
}

namespace {

// Convert (Br, Bphi, Bz) at azimuth phi to Cartesian (Bx, By, Bz).
inline Eigen::Vector3d cylToCart(double Br, double Bphi, double Bz, double phi) {
  const double cphi = std::cos(phi);
  const double sphi = std::sin(phi);
  return Eigen::Vector3d(Br * cphi - Bphi * sphi,
                         Br * sphi + Bphi * cphi,
                         Bz);
}

}  // namespace

Eigen::Vector3d ScalarPotentialFieldCorrection::getCorrectionAt(
    const GlobalPoint &pos, const std::vector<double> &corparms) const {
  const double r = std::hypot(pos.x(), pos.y());
  const double z = pos.z();
  const double phi = std::atan2(pos.y(), pos.x());

  basis_.SetPoint(r, z, phi);
  basis_.EvalBr();
  basis_.EvalBz();
  basis_.EvalBphi();

  double Br = 0., Bphi = 0., Bz = 0.;
  for (unsigned int i = 0; i < nModes_; ++i) {
    const unsigned int k = activeIdx_[i];
    const double c = corparms[basisGlobalIdx_[i]];
    Br   += c * basis_.GetBr_k(k);
    Bz   += c * basis_.GetBz_k(k);
    Bphi += c * basis_.GetBphi_k(k);
  }
  return cylToCart(Br, Bphi, Bz, phi);
}

void ScalarPotentialFieldCorrection::getBzBasisAt(
    const GlobalPoint &pos, std::vector<double> &bzPerMode) const {
  const double r = std::hypot(pos.x(), pos.y());
  const double z = pos.z();
  const double phi = std::atan2(pos.y(), pos.x());

  basis_.SetPoint(r, z, phi);
  basis_.EvalBz();

  bzPerMode.resize(nModes_);
  for (unsigned int i = 0; i < nModes_; ++i) {
    bzPerMode[i] = basis_.GetBz_k(activeIdx_[i]);
  }
}

void ScalarPotentialFieldCorrection::getJacobianAt(
    const GlobalPoint &pos,
    Eigen::Matrix<double, 3, Eigen::Dynamic> &dBdc) const {
  const double r = std::hypot(pos.x(), pos.y());
  const double z = pos.z();
  const double phi = std::atan2(pos.y(), pos.x());

  basis_.SetPoint(r, z, phi);
  basis_.EvalBr();
  basis_.EvalBz();
  basis_.EvalBphi();

  dBdc.resize(3, nModes_);
  const double cphi = std::cos(phi);
  const double sphi = std::sin(phi);
  for (unsigned int i = 0; i < nModes_; ++i) {
    const unsigned int k = activeIdx_[i];
    const double br = basis_.GetBr_k(k);
    const double bphi = basis_.GetBphi_k(k);
    const double bz = basis_.GetBz_k(k);
    dBdc(0, i) = br * cphi - bphi * sphi;
    dBdc(1, i) = br * sphi + bphi * cphi;
    dBdc(2, i) = bz;
  }
}

}  // namespace ana_hitanalyzer
