#include "GeantPropagatorESProducer.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackPropagation/Geant4e/interface/Geant4ePropagator.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TrackPropagation/Geant4e/interface/CGFQoPBlock.h"

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
  // The CVH energy-loss switches used to be read from CVH_* environment
  // variables inside the Geant4 model classes, which left them out of the
  // output file's provenance entirely.  They are ParameterSet parameters now,
  // pushed once into the process-global readers here -- in the CONSTRUCTOR, so
  // they are set before any physics list or fluctuation model is built.
  cvhcgf::configure(p);

  // A knob that belongs to the OTHER weight is silently inert -- the one cost
  // of having both estimators in one build: a job configured with
  // `IoniTruncationAlpha = 0.995` under the Fisher weight, or with
  // `CgfRecentreDamping = 0.4` under the legacy one, runs happily and
  // produces the DEFAULT physics while its provenance says otherwise.  Say so
  // once, at construction, rather than let a scan of a dead parameter be
  // reported as a null result.
  {
    const int mode = cvhcgf::cgfQoPMode();
    const double alpha = p.existsAs<double>("IoniTruncationAlpha")
                             ? p.getParameter<double>("IoniTruncationAlpha")
                             : 0.999;
    std::string dead;
    if (mode == 0) {
      if (p.existsAs<bool>("CgfRadiativeChannel") && p.getParameter<bool>("CgfRadiativeChannel"))
        dead += " CgfRadiativeChannel";
      if (p.existsAs<int>("CgfQoPRefresh") && p.getParameter<int>("CgfQoPRefresh") != 0)
        dead += " CgfQoPRefresh";
      if (p.existsAs<int>("IoniKokoulinCgfNbin") && p.getParameter<int>("IoniKokoulinCgfNbin") != 0)
        dead += " IoniKokoulinCgfNbin";
      if (p.existsAs<double>("CgfRecentreDamping") && p.getParameter<double>("CgfRecentreDamping") != 1.)
        dead += " CgfRecentreDamping";
      if (!dead.empty())
        edm::LogWarning("Geant4e") << "CgfQoPMode = 0 (legacy truncated-Gaussian Q): the CGF "
                                      "parameter(s)" << dead << " are set but have NO effect.";
    } else if (alpha != 0.999) {
      edm::LogWarning("Geant4e") << "CgfQoPMode = " << mode << " (Fisher weight): IoniTruncationAlpha = "
                                 << alpha << " is set but has NO effect -- the block CGF needs no "
                                    "delta-ray cut. Set CgfQoPMode = 0 to make it live.";
    }
  }
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
  // hand-written PSets predating the parameter keep working.  Read
  // unconditionally, but USED only when `CgfQoPMode = 0` selects the legacy
  // truncated-Gaussian Q; the default Fisher weight ignores it.
  const double ioniTruncationAlpha = pset_.existsAs<double>("IoniTruncationAlpha")
                                         ? pset_.getParameter<double>("IoniTruncationAlpha")
                                         : 0.999;

  // Geant4e maximum step [mm]. Default 10.0 reproduces the previous hard-coded
  // behaviour exactly, so existing configs are unaffected.
  const double stepLengthLimit = pset_.existsAs<double>("StepLengthLimit")
                                     ? pset_.getParameter<double>("StepLengthLimit")
                                     : 10.0;

  return std::make_unique<Geant4ePropagator>(
      &(iRecord.get(magFieldToken_)), particleName, dir, plimit_, forCVH_,
      ioniTruncationAlpha, stepLengthLimit);
}
