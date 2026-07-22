#ifndef MagneticField_ParametrizedEngine_ScalarPot3DMagneticField_h
#define MagneticField_ParametrizedEngine_ScalarPot3DMagneticField_h

/** \class ScalarPot3DMagneticField
 *
 *  CMSSW MagneticField subclass that evaluates the spherical-harmonic
 *  scalar-potential basis (mfs/harmonic_basis.py) for the CMS tracker
 *  volume.  Published by the `ParametrizedMagneticFieldProducer`
 *  ESProducer at a labelled `IdealMagneticFieldRecord` (default label
 *  `"ScalarPot3DMf"`), consumed by the CVH-side `MagneticFieldLabel`
 *  overrides in `nano_cff.py` / `runCvhJpsi.py`.  Outside the validity
 *  sphere `inTesla` returns the zero vector.
 *
 *  scalar-potential field model.  Per-event coefficient
 *  updates (per-event corFiles updates) are deferred — at this stage the field
 *  uses the dump file's stored coefficients directly.
 *
 *  \author David Walter
 */

#include <memory>
#include <string>

#include "MagneticField/Engine/interface/MagneticField.h"

namespace edm { class ParameterSet; }

namespace magfieldparam {
class ScalarPot3DEval;
}

class ScalarPot3DMagneticField : public MagneticField {
public:
  /// Construct from cfi PSet. Required keys:
  /// InitFile (string) — path to a coefficient dump
  /// ValidityRadius (double, cm) — sphere radius outside which inTesla
  /// returns zero (default 320 cm)
  explicit ScalarPot3DMagneticField(const edm::ParameterSet& parameters);

  ~ScalarPot3DMagneticField() override;

  GlobalVector inTesla(const GlobalPoint& gp) const override;
  GlobalVector inTeslaUnchecked(const GlobalPoint& gp) const override;

  bool isDefined(const GlobalPoint& gp) const override;

private:
  std::unique_ptr<magfieldparam::ScalarPot3DEval> eval_;
  double validity_radius_cm_;  // sphere radius for isDefined (cm)
};

#endif
