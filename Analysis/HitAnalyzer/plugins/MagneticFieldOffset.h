#ifndef HitAnalyzer_MagneticFieldOffset_h
#define HitAnalyzer_MagneticFieldOffset_h

/** \class MagneticFieldOffset
 *
 *  Field engine providing interpolation within the full CMS region.
 *
 *  \author N. Amapane - CERN
 */

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/VolumeBasedEngine/interface/MagGeometry.h"

// Class for testing MagneticFieldOffset
class testMagneticField;
class testMagGeometryAnalyzer;

class MagneticFieldOffset : public MagneticField {

 public:
   MagneticFieldOffset(const MagneticField *field) : field_(field), offset_(0.) {}
//   ~MagneticFieldOffset() override;

  /// Copy constructor implement a shallow copy (ie no ownership of actual engines)
//   MagneticFieldOffset(const MagneticFieldOffset& vbf);

  /// Returns a shallow copy.
  MagneticField* clone() const override;

  GlobalVector inTesla ( const GlobalPoint& g) const override;

  GlobalVector inTeslaUnchecked ( const GlobalPoint& g) const override;

  bool isDefined(const GlobalPoint& gp) const override { return field_->isDefined(gp); }
  
  void setOffset(float offset) { offset_ = offset; }

//   bool isZSymmetric() const;


 private:
  const MagneticField* field_;
  float offset_;
};

#endif
