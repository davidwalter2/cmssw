// HadElasticKernelWatcher -- record the deflection the SIMULATION actually
// applies to the primary at each `hadElastic` collision.
//
// WHY THIS EXISTS
// ---------------
// The nuclear elastic CF channel (calibration_studies/resolution/
// cf_nucel_exact.py) builds its angular kernel from a standalone driver that
// calls the G4 model's ApplyYourself() directly.  Every check on that kernel so
// far has been INFERENCE FROM THE DRIVER ALONE -- cross sections, model class,
// dispatch, energy ranges -- and all of it says the driver reproduces the
// physics list.  Yet the antiproton over-corrects by ~10x in locx, and four
// hypotheses have been tested and refuted (the Coulomb branch inside
// G4AntiNuclElastic at P=5e-4; acceptance truncation of its hard tail, 300x too
// small; the XSFactor*Elastic scale factors, all no-ops; the dropped recoil
// energy loss, which fails its own qop prediction).
//
// This is the measurement nobody has made: what the SIM draws, per collision,
// compared against what the driver draws.  It settles rate and shape at once,
// and it validates the pi/K/p kernels too rather than leaving them resting on a
// closure that happens to work.
//
// WHAT IS RECORDED, AND THE CONTAMINATION CONTROL
// -----------------------------------------------
// The deflection is taken as the angle between the pre- and post-step momentum
// directions of a step whose PostStep process is `hadElastic`.  That angle is
// NOT purely elastic: transportation bends the track in the field along the
// step and msc adds its own kick at the end of the same step.  Rather than
// disable either -- which would change the stepping and therefore the physics
// the closure is judged on -- every step is recorded with a TAG, and steps that
// did NOT end in hadElastic are written (prescaled) as a control.  The control
// distribution IS the field+msc contamination, measured under exactly the same
// conditions, so the signal can be read against it instead of being assumed
// clean.  At pT = 3 the expected contamination is ~1.5 mrad of field bend per
// ~4 mm step against a 25-35 mrad elastic kick.
//
// Record: 6 x float32 = [tag, theta_rad, ekin_pre_MeV, dE_MeV, steplen_mm, r_cm]
//   tag = 1 -> the step ended in `process` (default hadElastic)
//   tag = 0 -> control step, prescaled by `controlPrescale`
// `r` is kept so the offline side can cut to the toy tracker (r < 123 cm) and
// discard everything the primary does after it leaves for the real CMS
// calorimeters -- which is where most of its elastic collisions happen and none
// of them are in the modelled path.
//
// Nothing here is enabled by default; the watcher must be added explicitly.

#include "SimG4Core/Watcher/interface/SimWatcher.h"
#include "SimG4Core/Watcher/interface/SimWatcherFactory.h"
#include "SimG4Core/Notification/interface/Observer.h"
#include "SimG4Core/Notification/interface/EndOfRun.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

class HadElasticKernelWatcher : public SimWatcher,
                                public Observer<const G4Step *>,
                                public Observer<const EndOfRun *> {
public:
  explicit HadElasticKernelWatcher(const edm::ParameterSet &p)
      : out_(p.getUntrackedParameter<std::string>("output", "hadelastic.bin")),
        proc_(p.getUntrackedParameter<std::string>("process", "hadElastic")),
        prescale_(p.getUntrackedParameter<int>("controlPrescale", 1000)),
        rmax_(p.getUntrackedParameter<double>("rmax", 123.0)) {
    fh_.open(out_.c_str(), std::ios::binary | std::ios::trunc);
    std::cout << "[hadelkern] writing " << out_ << "  process=" << proc_
              << "  controlPrescale=" << prescale_ << "  rmax=" << rmax_ << " cm" << std::endl;
  }

  ~HadElasticKernelWatcher() override { close(); }

  void update(const G4Step *s) override {
    if (s == nullptr) return;
    const G4Track *t = s->GetTrack();
    if (t == nullptr || t->GetParentID() != 0) return;   // PRIMARY only

    const G4StepPoint *pre = s->GetPreStepPoint();
    const G4StepPoint *post = s->GetPostStepPoint();
    if (pre == nullptr || post == nullptr) return;

    const G4VProcess *dp = post->GetProcessDefinedStep();
    const bool sig = (dp != nullptr && dp->GetProcessName() == proc_);

    ++nseen_;
    if (!sig) {
      // control: prescaled, and only where it can be compared like for like
      if (prescale_ <= 0 || (nseen_ % prescale_) != 0) return;
    }

    const G4ThreeVector &d0 = pre->GetMomentumDirection();
    const G4ThreeVector &d1 = post->GetMomentumDirection();
    double c = d0.dot(d1);
    c = std::max(-1.0, std::min(1.0, c));
    const double theta = std::acos(c);

    const G4ThreeVector &x = post->GetPosition();
    const double r = std::sqrt(x.x() * x.x() + x.y() * x.y()) / CLHEP::cm;
    if (r > rmax_) return;          // outside the toy tracker: not modelled

    rec_[0] = sig ? 1.0f : 0.0f;
    rec_[1] = (float)theta;
    rec_[2] = (float)(pre->GetKineticEnergy() / CLHEP::MeV);
    rec_[3] = (float)((pre->GetKineticEnergy() - post->GetKineticEnergy()) / CLHEP::MeV);
    rec_[4] = (float)(s->GetStepLength() / CLHEP::mm);
    rec_[5] = (float)r;
    fh_.write(reinterpret_cast<const char *>(rec_.data()), rec_.size() * sizeof(float));
    if (sig) ++nsig_;
    else ++nctl_;
  }

  void update(const EndOfRun *) override { close(); }

private:
  void close() {
    if (closed_) return;
    closed_ = true;
    if (fh_.is_open()) fh_.close();
    std::cout << "[hadelkern] wrote " << nsig_ << " " << proc_ << " records and " << nctl_
              << " control records (of " << nseen_ << " primary steps) to " << out_ << std::endl;
  }

  const std::string out_;
  const std::string proc_;
  const int prescale_;
  const double rmax_;
  std::ofstream fh_;
  bool closed_ = false;
  long nseen_ = 0, nsig_ = 0, nctl_ = 0;
  std::vector<float> rec_ = std::vector<float>(6, 0.f);
};

DEFINE_SIMWATCHER(HadElasticKernelWatcher);
