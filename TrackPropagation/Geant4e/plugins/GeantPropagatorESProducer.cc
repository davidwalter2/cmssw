#include "GeantPropagatorESProducer.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackPropagation/Geant4e/interface/Geant4ePropagator.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ModuleFactory.h"

#include <memory>
#include <string>

using namespace edm;

GeantPropagatorESProducer::GeantPropagatorESProducer(const edm::ParameterSet &p)
    : fieldlabel_(p.getParameter<std::string>("MagneticFieldLabel")),
      magFieldToken_(setWhatProduced(this, p.getParameter<std::string>("ComponentName"))
                         .consumesFrom<MagneticField, IdealMagneticFieldRecord>(edm::ESInputTag("", fieldlabel_))) {
  pset_ = p;
  plimit_ = pset_.getParameter<double>("PropagationPtotLimit");
  forCVH_ = pset_.getParameter<bool>("ForCVH");
}

GeantPropagatorESProducer::~GeantPropagatorESProducer() {}

std::unique_ptr<Propagator> GeantPropagatorESProducer::produce(const TrackingComponentsRecord &iRecord) {
  std::string pdir = pset_.getParameter<std::string>("PropagationDirection");
  std::string particleName = pset_.getParameter<std::string>("ParticleName");

  PropagationDirection dir = alongMomentum;

  if (pdir == "oppositeToMomentum") {
    dir = oppositeToMomentum;
  } else if (pdir == "alongMomentum") {
    dir = alongMomentum;
  } else if (pdir == "anyDirection") {
    dir = anyDirection;
  }

  // Delta-electron truncation of the ionization variance; optional so
  // hand-written PSets predating the parameter keep working.
  const double ioniTruncAlpha = pset_.existsAs<double>("IoniTruncationAlpha")
                                    ? pset_.getParameter<double>("IoniTruncationAlpha")
                                    : 0.999;

  return std::make_unique<Geant4ePropagator>(
      &(iRecord.get(magFieldToken_)), particleName, dir, plimit_, forCVH_, ioniTruncAlpha);
}
