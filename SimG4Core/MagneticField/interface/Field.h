#ifndef SimG4Core_Field_H
#define SimG4Core_Field_H

#include "G4MagneticField.hh"

#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include <CLHEP/Units/SystemOfUnits.h>

class G4LogicalVolume;

namespace sim {
  // Volume-resolved material-scaling provider (global material model,
  // Analysis/HitAnalyzer/doc/global-material-model-plan.md). When a
  // provider is attached to the field wrapper, the CVH energy-loss
  // process adds materialOffset(volume, r, z) of the step's volume to
  // the leg-constant dxi. Implemented by
  // TrackPropagation/Geant4e MaterialGroupModel; the abstract base
  // lives here so SimG4Core does not depend on TrackPropagation.
  class MaterialOffsetProvider {
  public:
    virtual ~MaterialOffsetProvider() = default;
    virtual double materialOffset(const G4LogicalVolume *lv, double r_cm, double z_cm) const = 0;
  };

  // Per-step field-mode provider (leg-structure-free B-field correction):
  // supplies the current correction field and the per-mode basis values at
  // arbitrary points, so the propagator can apply the correction and
  // attribute the per-mode derivatives per Geant4 step instead of
  // piecewise-constant per leg. Implemented maker-side on top of the
  // scalar-potential basis evaluator; the abstract base lives here so
  // SimG4Core/TrackPropagation do not depend on Analysis/HitAnalyzer.
  class FieldModeProvider {
  public:
    virtual ~FieldModeProvider() = default;
    virtual unsigned int nModes() const = 0;
    // One evaluation per step: the applied correction field b [T]
    // (coefficients contracted with the basis, plus any FD injection) AND
    // the per-mode basis values at the same global point (cm). Using one
    // sample for both the application and the derivative columns makes
    // them consistent by construction. The returned pointers stay valid
    // until the next call on the same provider instance.
    virtual void sampleAt(double x_cm, double y_cm, double z_cm, double b[3],
                          const double *&bx, const double *&by, const double *&bz) const = 0;
  };

  class Field final : public G4MagneticField {
  public:
    Field(const MagneticField *f, double d);
    ~Field() override;
    inline void GetFieldValue(const G4double p[], G4double b[3]) const override;
    void SetOffset(double x, double y, double z);
    void SetMaterialOffset(double offset) { dxi = offset; }
    double GetMaterialOffset() const { return dxi; }
    void SetMaterialOffsetProvider(const MaterialOffsetProvider *p) { offsetProvider = p; }
    const MaterialOffsetProvider *GetMaterialOffsetProvider() const { return offsetProvider; }

  private:
    const MagneticField *theCMSMagneticField;
    double theDelta;

    mutable double oldx[3];
    mutable double oldb[3];

    double offset[3];
    double dxi;
    const MaterialOffsetProvider *offsetProvider = nullptr;
  };
};  // namespace sim

void sim::Field::GetFieldValue(const G4double xyz[], G4double bfield[3]) const {
  if (std::abs(oldx[0] - xyz[0]) > theDelta || std::abs(oldx[1] - xyz[1]) > theDelta ||
      std::abs(oldx[2] - xyz[2]) > theDelta) {
    constexpr float lunit = (1.0 / CLHEP::cm);
    GlobalPoint ggg((float)(xyz[0]) * lunit, (float)(xyz[1]) * lunit, (float)(xyz[2]) * lunit);
    GlobalVector v = theCMSMagneticField->inTesla(ggg);

    constexpr float btesla = CLHEP::tesla;
    oldb[0] = (v.x() * btesla);
    oldb[1] = (v.y() * btesla);
    oldb[2] = (v.z() * btesla);
    oldx[0] = xyz[0];
    oldx[1] = xyz[1];
    oldx[2] = xyz[2];
  }

  bfield[0] = oldb[0] + offset[0];
  bfield[1] = oldb[1] + offset[1];
  bfield[2] = oldb[2] + offset[2];
}

#endif
