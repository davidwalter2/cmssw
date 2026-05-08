/** \file
 *
 *  Phase A.1 of replicated-bouncing-cloud.md.
 *
 *  \author David Walter
 */

#include "ScalarPot3DMagneticField.h"

#include "MagneticField/ParametrizedEngine/interface/ScalarPot3DEval.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <cmath>
#include <stdexcept>

ScalarPot3DMagneticField::ScalarPot3DMagneticField(const edm::ParameterSet& parameters)
    : eval_(nullptr),
      validity_radius_cm_(parameters.exists("ValidityRadius")
                          ? parameters.getParameter<double>("ValidityRadius")
                          : 320.0) {
  const std::string init_file = parameters.getParameter<std::string>("InitFile");
  if (init_file.empty()) {
    throw cms::Exception("Configuration")
        << "ScalarPot3DMagneticField requires a non-empty InitFile parameter";
  }
  eval_ = std::make_unique<magfieldparam::ScalarPot3DEval>(init_file);
  edm::LogInfo("MagneticField|ScalarPot3D")
      << "Loaded ScalarPot3DEval from " << init_file
      << ":  nModes=" << eval_->nModes()
      << "  l_max=" << eval_->lMax()
      << "  r_scale=" << eval_->rScale() << " cm"
      << "  cmssw_norm=" << eval_->cmsswNorm()
      << "  validity_radius=" << validity_radius_cm_ << " cm";
}

ScalarPot3DMagneticField::~ScalarPot3DMagneticField() = default;

GlobalVector ScalarPot3DMagneticField::inTesla(const GlobalPoint& gp) const {
  if (isDefined(gp)) {
    return inTeslaUnchecked(gp);
  }
  return GlobalVector();
}

GlobalVector ScalarPot3DMagneticField::inTeslaUnchecked(const GlobalPoint& gp) const {
  // For Phase A.1 we use the dump file's stored coefficients directly
  // as the basis amplitudes.  Per-event coefficient updates (A.4)
  // will replace this with a thread-safe shared accessor.
  return eval_->evaluateAbsoluteAt(gp, eval_->initCoeffs());
}

bool ScalarPot3DMagneticField::isDefined(const GlobalPoint& gp) const {
  // Sphere |R| <= validity_radius_cm_, with R = sqrt(x^2 + y^2 + z^2).
  // The basis is fitted on a similar sphere; outside, R^L extrapolation
  // becomes unreliable.  A.2 will replace this hard cutoff with a
  // C^1 cosine blend to the underlying CMSSW field.
  const double R = std::sqrt(gp.x() * gp.x() + gp.y() * gp.y()
                             + gp.z() * gp.z());
  return R <= validity_radius_cm_;
}
