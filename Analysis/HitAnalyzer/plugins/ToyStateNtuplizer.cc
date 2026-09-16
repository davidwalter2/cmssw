// Records the TRUE 5D track state where the primary crosses each toy shell
// boundary, replacing PSimHits for the homogeneous clean-propagation test.
//
// WHY THIS EXISTS. SimHitStateNtuplizer flattens PSimHits, which only exist
// where the tracker sensitive detector runs -- that is what ties the clean
// propagation test to the real tracker geometry and its DetIds. In the toy the
// tracking volume is one uniform medium (Analysis/HitAnalyzer/data/tracker.xml)
// so there are no sensors and no PSimHits, and the state has to be taken from
// Geant4 directly.
//
// WHY BOUNDARIES RATHER THAN INTERPOLATION. The toy volume is sliced into
// concentric shells of IDENTICAL material. Geant4 always terminates a step on a
// geometry boundary, so the post-step point of a boundary-limited step lies
// EXACTLY on the surface: no interpolation, no sagitta error. Because the
// material is the same on both sides, the slicing is pure bookkeeping and
// changes no physics -- except that it truncates steps at the boundaries, which
// is also what module boundaries do in the real detector.
//
// WHAT IS STORED, AND WHY IT IS GLOBAL. The state is written in GLOBAL
// coordinates (position and momentum), not in a per-surface local frame.
//
// An earlier version built the tangent basis from each ray's OWN azimuth and
// set locx = pos.eU. That is identically zero -- a position has no component
// along its own tangent -- and the output confirmed it (|locx| < 1e-14 for
// every entry). The real lesson is that the local frame is a property of the
// SURFACE, fixed by the reference trajectory, and the watcher does not know the
// reference. Choosing a frame here can only be wrong or redundant.
//
// So the frame is applied offline, where the reference IS known, and where the
// same plane definitions are handed to the model side. Global quantities are
// sufficient: any local frame is a rotation away.

#include "SimG4Core/Watcher/interface/SimProducer.h"
#include "SimG4Core/Watcher/interface/SimWatcherFactory.h"
#include "SimG4Core/Notification/interface/Observer.h"
#include "SimG4Core/Notification/interface/BeginOfEvent.h"
#include "SimG4Core/Notification/interface/EndOfEvent.h"
#include "SimG4Core/Notification/interface/EndOfRun.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <TFile.h>
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4VProcess.hh"

#include <TTree.h>
#include <cmath>
#include <vector>

class ToyStateNtuplizer : public SimProducer,
                          public Observer<const BeginOfEvent *>,
                          public Observer<const G4Step *>,
                          public Observer<const EndOfEvent *>,
                          public Observer<const EndOfRun *> {
public:
  explicit ToyStateNtuplizer(const edm::ParameterSet &);

  ~ToyStateNtuplizer() override { writeOut(); }

  void produce(edm::Event &, const edm::EventSetup &) override {}

private:
  void update(const BeginOfEvent *) override;
  void update(const G4Step *) override;
  void update(const EndOfEvent *) override;
  void update(const EndOfRun *) override;
  // Written at EndOfRun, NOT in the destructor: the destructor runs after ROOT
  // has torn down and the file came out with no keys at all.
  void writeOut();

  // target radii [cm], ascending
  std::vector<double> radii_;
  double tol_;   // how close to a radius counts as "on" it [cm]

  // A plain TFile, NOT TFileService: the watcher is constructed on a Geant4
  // worker thread where the edm service registry is not reachable, and calling
  // fs->make<TTree> there segfaults. Owning the file keeps the watcher
  // self-contained; the job is pinned to one thread/stream anyway.
  TFile *file_ = nullptr;
  TTree *tree_ = nullptr;
  std::vector<int> ishell_;
  std::vector<float> qop_;
  std::vector<float> globx_, globy_, globz_;      // position [cm]
  std::vector<float> globpx_, globpy_, globpz_;   // momentum [GeV]
  std::vector<float> pabs_, eloss_, globr_;
  int nhit_ = 0;
  float p0_ = 0.f;   // momentum at the first recorded crossing, for eloss
};

ToyStateNtuplizer::ToyStateNtuplizer(const edm::ParameterSet &p)
    : radii_(p.getParameter<std::vector<double>>("radii")),
      tol_(p.getParameter<double>("tolerance")) {
  file_ = TFile::Open(p.getParameter<std::string>("output").c_str(), "RECREATE");
  file_->cd();
  tree_ = new TTree("simstates", "true 5D local states at toy shell boundaries");
  tree_->Branch("detid", &ishell_);   // shell index, standing in for a DetId
  tree_->Branch("qop", &qop_);
  tree_->Branch("globx", &globx_);
  tree_->Branch("globy", &globy_);
  tree_->Branch("globz", &globz_);
  tree_->Branch("globpx", &globpx_);
  tree_->Branch("globpy", &globpy_);
  tree_->Branch("globpz", &globpz_);
  tree_->Branch("pabs", &pabs_);
  tree_->Branch("eloss", &eloss_);
  tree_->Branch("globr", &globr_);
  tree_->Branch("nhit", &nhit_);
}

void ToyStateNtuplizer::update(const BeginOfEvent *) {
  ishell_.clear();
  qop_.clear();
  globx_.clear(); globy_.clear(); globz_.clear();
  globpx_.clear(); globpy_.clear(); globpz_.clear();
  pabs_.clear(); eloss_.clear(); globr_.clear();
  nhit_ = 0;
  p0_ = 0.f;
}

void ToyStateNtuplizer::update(const G4Step *step) {
  const G4Track *trk = step->GetTrack();
  // primary only: secondaries (delta rays) are not the propagated particle
  if (trk->GetParentID() != 0)
    return;

  const G4StepPoint *post = step->GetPostStepPoint();
  // a boundary-limited step ends exactly ON a surface; anything else is a
  // physics-limited step in the middle of a shell and is not a crossing
  if (post->GetStepStatus() != fGeomBoundary)
    return;

  const G4ThreeVector pos = post->GetPosition() / CLHEP::cm;
  const double r = std::hypot(pos.x(), pos.y());

  int ishell = -1;
  for (size_t k = 0; k < radii_.size(); ++k) {
    if (std::abs(r - radii_[k]) < tol_) {
      ishell = static_cast<int>(k);
      break;
    }
  }
  if (ishell < 0)
    return;   // a boundary, but not one of the surfaces we score on

  const G4ThreeVector mom = post->GetMomentum() / CLHEP::GeV;
  const double p = mom.mag();
  if (p <= 0.)
    return;

  if (p0_ == 0.f)
    p0_ = static_cast<float>(p);

  ishell_.push_back(ishell);
  qop_.push_back(static_cast<float>(trk->GetDefinition()->GetPDGCharge() > 0 ? 1. / p : -1. / p));
  globx_.push_back(static_cast<float>(pos.x()));
  globy_.push_back(static_cast<float>(pos.y()));
  globz_.push_back(static_cast<float>(pos.z()));
  globpx_.push_back(static_cast<float>(mom.x()));
  globpy_.push_back(static_cast<float>(mom.y()));
  globpz_.push_back(static_cast<float>(mom.z()));
  pabs_.push_back(static_cast<float>(p));
  eloss_.push_back(static_cast<float>(p0_ - p));
  globr_.push_back(static_cast<float>(r));
  ++nhit_;
}

void ToyStateNtuplizer::update(const EndOfEvent *) { tree_->Fill(); }

void ToyStateNtuplizer::update(const EndOfRun *) { writeOut(); }

void ToyStateNtuplizer::writeOut() {
  if (!file_)
    return;
  file_->cd();
  tree_->Write();
  file_->Close();
  file_ = nullptr;
  tree_ = nullptr;
}

DEFINE_SIMWATCHER(ToyStateNtuplizer);
