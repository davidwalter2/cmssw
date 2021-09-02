#include "MagneticFieldOffset.h"

MagneticField* MagneticFieldOffset::clone() const {
  return new MagneticFieldOffset(*this);
}

GlobalVector MagneticFieldOffset::inTesla (const GlobalPoint& gp) const {
  GlobalVector fieldval = field_->inTesla(gp);
  const float mag = fieldval.mag();
  return mag > 0. ? (mag + offset_)/mag*fieldval : fieldval;
}

GlobalVector MagneticFieldOffset::inTeslaUnchecked(const GlobalPoint& gp) const{
  //same as above, but do not check range
  GlobalVector fieldval = field_->inTeslaUnchecked(gp);
  const float mag = fieldval.mag();
  return mag > 0. ? (mag + offset_)/mag*fieldval : fieldval;
}
