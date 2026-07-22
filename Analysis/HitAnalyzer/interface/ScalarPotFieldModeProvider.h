#ifndef Analysis_HitAnalyzer_ScalarPotFieldModeProvider_h
#define Analysis_HitAnalyzer_ScalarPotFieldModeProvider_h

// Adapter implementing sim::FieldModeProvider on top of the
// scalar-potential basis evaluator, for the per-step (leg-structure-free)
// field-mode application/attribution in the Geant4e propagator. Owned per
// residual-maker stream instance; the mutable basis buffers make it
// single-stream only. The corparms pointer refers to the maker's live
// corrections vector, so calibration iterations are picked up
// automatically. setInjection adds eps * basis_i to the applied
// correction field for the per-step FD closure.

#include "Analysis/HitAnalyzer/interface/ScalarPotentialFieldCorrection.h"
#include "SimG4Core/MagneticField/interface/Field.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include <vector>

namespace ana_hitanalyzer {

class ScalarPotFieldModeProvider : public sim::FieldModeProvider {
public:
  ScalarPotFieldModeProvider(const ScalarPotentialFieldCorrection *fc,
                             const std::vector<double> *corparms)
      : fc_(fc), corparms_(corparms) {}

  unsigned int nModes() const override { return fc_->nModes(); }

  void setInjection(int mode, double eps) {
    injMode_ = mode;
    injEps_ = eps;
  }

  void sampleAt(double x_cm, double y_cm, double z_cm, double b[3],
                const double *&bx, const double *&by, const double *&bz) const override {
    const GlobalPoint p(x_cm, y_cm, z_cm);
    fc_->getBxBasisAt(p, bx_);
    fc_->getByBasisAt(p, by_);
    fc_->getBzBasisAt(p, bz_);
    bx = bx_.data();
    by = by_.data();
    bz = bz_.data();
    // applied correction = coefficients (by global index) contracted with
    // the basis, plus the FD injection
    double vx = 0., vy = 0., vz = 0.;
    const unsigned int n = fc_->nModes();
    for (unsigned int i = 0; i < n; ++i) {
      const double c = (*corparms_)[fc_->basisGlobalIdx(i)];
      vx += c * bx_[i];
      vy += c * by_[i];
      vz += c * bz_[i];
    }
    if (injMode_ >= 0) {
      vx += injEps_ * bx_[injMode_];
      vy += injEps_ * by_[injMode_];
      vz += injEps_ * bz_[injMode_];
    }
    b[0] = vx;
    b[1] = vy;
    b[2] = vz;
  }

private:
  const ScalarPotentialFieldCorrection *fc_;
  const std::vector<double> *corparms_;
  int injMode_ = -1;
  double injEps_ = 0.;
  mutable std::vector<double> bx_, by_, bz_;
};

}  // namespace ana_hitanalyzer

#endif
