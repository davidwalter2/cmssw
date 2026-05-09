#include "Analysis/HitAnalyzer/interface/ScalarPotentialFieldCorrection.h"

#include <stdexcept>

#include "MagneticField/ParametrizedEngine/interface/ScalarPot3DEval.h"

namespace ana_hitanalyzer {

ScalarPotentialFieldCorrection::ScalarPotentialFieldCorrection(
    const std::string& dumpPath)
    : eval_(std::make_unique<magfieldparam::ScalarPot3DEval>(dumpPath)),
      basisGlobalIdx_(eval_->nModes(), 0u) {}

ScalarPotentialFieldCorrection::~ScalarPotentialFieldCorrection() = default;

unsigned int ScalarPotentialFieldCorrection::nModes() const {
  return eval_->nModes();
}

const std::vector<double>& ScalarPotentialFieldCorrection::initCoeffs() const {
  return eval_->initCoeffs();
}

void ScalarPotentialFieldCorrection::appendParmsetEntries(
    std::set<std::pair<int, DetId>>& parmset) const {
  const unsigned int N = nModes();
  for (unsigned int i = 0; i < N; ++i) {
    parmset.emplace(ParmTypeBfieldGlobal, DetId(i));
  }
}

void ScalarPotentialFieldCorrection::resolveGlobalIndices(
    const std::map<std::pair<int, DetId>, unsigned int>& detidparms) {
  const unsigned int N = nModes();
  basisGlobalIdx_.assign(N, 0u);
  for (unsigned int i = 0; i < N; ++i) {
    auto it = detidparms.find(std::make_pair(ParmTypeBfieldGlobal, DetId(i)));
    if (it == detidparms.end()) {
      throw std::runtime_error(
          "ScalarPotentialFieldCorrection: mode global index not in detidparms");
    }
    basisGlobalIdx_[i] = it->second;
  }
}

Eigen::Vector3d ScalarPotentialFieldCorrection::getCorrectionAt(
    const GlobalPoint& pos, const std::vector<double>& corparms) const {
  // Gather per-mode coefficient slice from the global corparms vector
  // at the resolved basis indices. ScalarPot3DEval::evaluateAbsoluteAt
  // takes a vector sized to nModes() (one entry per mode).
  const unsigned int N = nModes();
  std::vector<double> c(N);
  for (unsigned int i = 0; i < N; ++i) {
    c[i] = corparms[basisGlobalIdx_[i]];
  }
  const GlobalVector B = eval_->evaluateAbsoluteAt(pos, c);
  return Eigen::Vector3d(B.x(), B.y(), B.z());
}

void ScalarPotentialFieldCorrection::getBzBasisAt(
    const GlobalPoint& pos, std::vector<double>& out) const {
  eval_->getBzBasisAt(pos, out);
}

void ScalarPotentialFieldCorrection::getBxBasisAt(
    const GlobalPoint& pos, std::vector<double>& out) const {
  eval_->getBxBasisAt(pos, out);
}

void ScalarPotentialFieldCorrection::getByBasisAt(
    const GlobalPoint& pos, std::vector<double>& out) const {
  eval_->getByBasisAt(pos, out);
}

void ScalarPotentialFieldCorrection::getJacobianAt(
    const GlobalPoint& pos,
    Eigen::Matrix<double, 3, Eigen::Dynamic>& dBdc) const {
  const unsigned int N = nModes();
  std::vector<double> bx, by, bz;
  eval_->getBxBasisAt(pos, bx);
  eval_->getByBasisAt(pos, by);
  eval_->getBzBasisAt(pos, bz);
  dBdc.resize(3, N);
  for (unsigned int i = 0; i < N; ++i) {
    dBdc(0, i) = bx[i];
    dBdc(1, i) = by[i];
    dBdc(2, i) = bz[i];
  }
}

}  // namespace ana_hitanalyzer
