#ifndef HitAnalyzer_OffsetMagneticField_h
#define HitAnalyzer_OffsetMagneticField_h

#include "MagneticField/Engine/interface/MagneticField.h"

class OffsetMagneticField : public MagneticField {
  
public:
  OffsetMagneticField() {};
  OffsetMagneticField(const MagneticField* nominal, GlobalVector const& offset) :
    field_(nominal), offset_(offset) {}
    
  virtual ~OffsetMagneticField() {}
  
  GlobalVector const& offset() const { return offset_; }
  void setOffset(GlobalVector const& offset) { offset_ = offset; }
  
  virtual GlobalVector inTesla(const GlobalPoint& gp) const { return field_->inTesla(gp) + offset_; }
  virtual GlobalVector inTeslaUnchecked(const GlobalPoint& gp) const { return field_->inTeslaUnchecked(gp) + offset_; }
    
private:
  const MagneticField* field_;
  GlobalVector offset_;
  
};
#endif
