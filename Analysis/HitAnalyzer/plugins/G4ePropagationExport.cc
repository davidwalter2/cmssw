// G4ePropagationExport
//
// Model side of the "clean propagation test": propagate ONE fixed initial
// state deterministically with the CVH Geant4e propagator through the same
// sequence of sensor surfaces the full simulation crosses, and export
// everything needed to predict -- analytically, in characteristic-function
// space -- the PDF of the propagated 5D state that Geant4 would produce.
//
// Per leg (surface k-1 -> surface k) it writes:
//   * the deterministic reference state at surface k, in the DetUnit local
//     parameterization (q/p, dx/dz, dy/dz, x, y) -- directly comparable with
//     the PSimHit-derived true states from SimHitStateNtuplizer;
//   * the leg transport Jacobian F and the noise matrices Q, dQMS, dQI
//     (curvilinear), which give the exact Gaussian limit;
//   * the per-step Urban (ionization) and Moliere (multiple-scattering)
//     physics records -- identical layout to the maker's ioniurbanv/msmoliv
//     exports, so the existing offline CF code reads them unchanged;
//   * the per-step CUMULATIVE transport Jacobian, so that the noise of every
//     individual Geant4 step can be transported to any target surface
//     exactly, rather than through one RMS-matched scalar weight per pooled
//     block as the track-fit exports must do. In the fit that pooling is
//     unavoidable; here it would be a confound, because the whole point of
//     the test is the tails.
//
// One cmsRun event = one full propagation. The output is tiny and is consumed
// by calibration_studies/resolution/cf_propagation_test.py.

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"

#include "TrackPropagation/Geant4e/interface/Geant4ePropagator.h"
#include "TrackPropagation/Geant4e/interface/CvhMasterThread.h"
#include "TrackPropagation/Geant4e/interface/CvhMasterRecord.h"
#include "TrackPropagation/Geant4e/interface/CvhWorker.h"

#include "CLHEP/Random/RandomEngine.h"
#include "Randomize.hh"

#include <Eigen/Dense>

#include "TTree.h"

#include <cmath>
#include <memory>
#include <vector>

class G4ePropagationExport : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit G4ePropagationExport(const edm::ParameterSet &);

private:
  void analyze(const edm::Event &, const edm::EventSetup &) override;

  // local copy of ResidualGlobalCorrectionMakerBase::surfaceToDouble (private there)
  static GloballyPositioned<double> surfaceToDouble(const Surface &surface);

  edm::ESGetToken<Propagator, TrackingComponentsRecord> g4ePropToken_;
  edm::ESGetToken<CvhMasterThread, CvhMasterRecord> cvhMasterToken_;
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;
  // PSimHit local coordinates are always expressed in the IDEAL geometry (the
  // simulation knows nothing about the alignment applied at reco), so the
  // reference surfaces of the model side must be the ideal ones too or the
  // two sides are compared in different frames.
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomIdealToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> fieldToken_;
  bool useIdealGeometry_;

  std::unique_ptr<CvhWorker> worker_;
  std::unique_ptr<Geant4ePropagator> prop_;

  std::vector<double> initPos_, initMom_;
  double charge_;
  std::string particleName_;
  std::vector<unsigned int> targetDetIds_;
  std::vector<double> targetLocalZ_;

  bool done_ = false;

  TTree *tree_ = nullptr;

  unsigned int ileg_ = 0;
  unsigned int detid_ = 0;
  bool ok_ = false;
  double zoff_ = 0.;
  // deterministic reference state at the target plane
  double refqop_ = 0., refdxdz_ = 0., refdydz_ = 0., reflocx_ = 0., reflocy_ = 0., reflocz_ = 0.;
  double refglobr_ = 0., refglobz_ = 0., refp_ = 0., refpt_ = 0.;
  double dEdxlast_ = 0.;
  std::vector<double> F_, Q_, dQMS_, dQI_;
  // per-step physics records (same layout as the maker's exports)
  std::vector<double> msmoliv_;
  std::vector<double> ioniurbanv_;
  // per-step cumulative transport (25 doubles) + log-size markers
  std::vector<double> stepjacc_;
  std::vector<int> stepnms_, stepnioni_;
};

G4ePropagationExport::G4ePropagationExport(const edm::ParameterSet &iConfig)
    : g4ePropToken_(esConsumes(edm::ESInputTag("", "Geant4ePropagator"))),
      cvhMasterToken_(esConsumes()),
      geomToken_(esConsumes()),
      geomIdealToken_(esConsumes(edm::ESInputTag("", "idealForDigi"))),
      fieldToken_(esConsumes()),
      useIdealGeometry_(iConfig.getParameter<bool>("useIdealGeometry")),
      initPos_(iConfig.getParameter<std::vector<double>>("initialPosition")),
      initMom_(iConfig.getParameter<std::vector<double>>("initialMomentum")),
      charge_(iConfig.getParameter<double>("charge")),
      particleName_(iConfig.getParameter<std::string>("particleName")),
      targetDetIds_(iConfig.getParameter<std::vector<unsigned int>>("targetDetIds")),
      targetLocalZ_(iConfig.getParameter<std::vector<double>>("targetLocalZ")) {
  usesResource("TFileService");
  if (initPos_.size() != 3 || initMom_.size() != 3) {
    throw cms::Exception("Configuration") << "initialPosition and initialMomentum must have 3 entries";
  }
  if (targetDetIds_.size() != targetLocalZ_.size()) {
    throw cms::Exception("Configuration") << "targetDetIds and targetLocalZ must have the same length";
  }
  worker_ = std::make_unique<CvhWorker>();

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("legs", "deterministic reference propagation, per leg");
  tree_->Branch("ileg", &ileg_);
  tree_->Branch("detid", &detid_);
  tree_->Branch("ok", &ok_);
  tree_->Branch("zoff", &zoff_);
  tree_->Branch("refqop", &refqop_);
  tree_->Branch("refdxdz", &refdxdz_);
  tree_->Branch("refdydz", &refdydz_);
  tree_->Branch("reflocx", &reflocx_);
  tree_->Branch("reflocy", &reflocy_);
  tree_->Branch("reflocz", &reflocz_);
  tree_->Branch("refglobr", &refglobr_);
  tree_->Branch("refglobz", &refglobz_);
  tree_->Branch("refp", &refp_);
  tree_->Branch("refpt", &refpt_);
  tree_->Branch("dEdxlast", &dEdxlast_);
  tree_->Branch("F", &F_);
  tree_->Branch("Q", &Q_);
  tree_->Branch("dQMS", &dQMS_);
  tree_->Branch("dQI", &dQI_);
  tree_->Branch("msmoliv", &msmoliv_);
  tree_->Branch("ioniurbanv", &ioniurbanv_);
  tree_->Branch("stepjacc", &stepjacc_);
  tree_->Branch("stepnms", &stepnms_);
  tree_->Branch("stepnioni", &stepnioni_);
}

GloballyPositioned<double> G4ePropagationExport::surfaceToDouble(const Surface &surface) {
  const Point3DBase<double, GlobalTag> pos = surface.position();
  auto const &gy = surface.rotation().y();
  auto const &gz = surface.rotation().z();
  const Eigen::Matrix<double, 3, 1> vy(gy.x(), gy.y(), gy.z());
  const Eigen::Matrix<double, 3, 1> vz(gz.x(), gz.y(), gz.z());
  const Eigen::Matrix<double, 3, 1> uz = vz.normalized();
  const Eigen::Matrix<double, 3, 1> ux = vy.cross(vz).normalized();
  const Eigen::Matrix<double, 3, 1> uy = uz.cross(ux);
  const TkRotation<double> tkrot(ux[0], ux[1], ux[2], uy[0], uy[1], uy[2], uz[0], uz[1], uz[2]);
  return GloballyPositioned<double>(pos, tkrot);
}

void G4ePropagationExport::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
  if (done_) {
    return;
  }
  done_ = true;

  // Per-thread G4 setup, then bind the stream's CLHEP engine into Geant4's
  // thread-local engine, then take a private clone of the ES propagator so we
  // can switch the step logs on (the setters are non-const). Same ordering as
  // ResidualGlobalCorrectionMakerG4e::produce -- it matters.
  worker_->ensureInitialized(iSetup.getData(cvhMasterToken_).cvhMaster());

  edm::Service<edm::RandomNumberGenerator> rng;
  if (!rng.isAvailable()) {
    throw cms::Exception("Configuration") << "G4ePropagationExport requires the RandomNumberGeneratorService";
  }
  CLHEP::HepRandomEngine &engine = rng->getEngine(iEvent.streamID());
  G4Random::setTheEngine(&engine);

  if (!prop_) {
    const Geant4ePropagator *templateProp =
        dynamic_cast<const Geant4ePropagator *>(iSetup.getHandle(g4ePropToken_).product());
    if (templateProp == nullptr) {
      throw cms::Exception("Configuration") << "ESProducer 'Geant4ePropagator' is not a Geant4ePropagator";
    }
    if (!templateProp->GetForCVH()) {
      throw cms::Exception("Configuration") << "Geant4ePropagator must be configured with ForCVH = True";
    }
    prop_.reset(templateProp->clone());
    prop_->setIoniStepLogging(true);
    prop_->setStepTransportLogging(true);
  }

  const TrackerGeometry *geom =
      useIdealGeometry_ ? &iSetup.getData(geomIdealToken_) : &iSetup.getData(geomToken_);

  Eigen::Matrix<double, 7, 1> state;
  state << initPos_[0], initPos_[1], initPos_[2], initMom_[0], initMom_[1], initMom_[2], charge_;

  edm::LogPrint("G4ePropagationExport")
      << "start state: x=(" << state[0] << "," << state[1] << "," << state[2] << ") cm  p=(" << state[3] << ","
      << state[4] << "," << state[5] << ") GeV  q=" << state[6] << "  particle=" << particleName_ << "  targets="
      << targetDetIds_.size();

  for (size_t k = 0; k < targetDetIds_.size(); ++k) {
    const DetId did(targetDetIds_[k]);
    const GeomDet *det = geom->idToDet(did);
    if (det == nullptr) {
      throw cms::Exception("Configuration") << "detid " << targetDetIds_[k] << " not in the tracker geometry";
    }
    // reference frame = the DetUnit frame (the frame PSimHit local coordinates
    // live in); the propagation TARGET is that plane displaced along local z
    // to the sensor entry face, so the leg stops just before the silicon.
    const GloballyPositioned<double> base = surfaceToDouble(det->surface());
    const double zoff = targetLocalZ_[k];
    const GloballyPositioned<double> target(base.toGlobal(GloballyPositioned<double>::LocalPoint(0., 0., zoff)),
                                            base.rotation());

    auto const &res = prop_->propagateGenericWithJacobianAltD(state, target, Eigen::Vector3d::Zero(), 0., 0., 0., -1.,
                                                              particleName_);

    ileg_ = k;
    detid_ = targetDetIds_[k];
    ok_ = std::get<0>(res);
    zoff_ = zoff;

    F_.assign(25, 0.);
    Q_.assign(25, 0.);
    dQMS_.assign(25, 0.);
    dQI_.assign(25, 0.);
    msmoliv_.clear();
    ioniurbanv_.clear();
    stepjacc_.clear();
    stepnms_.clear();
    stepnioni_.clear();

    if (!ok_) {
      edm::LogWarning("G4ePropagationExport") << "propagation FAILED on leg " << k << " to detid " << detid_;
      tree_->Fill();
      break;
    }

    const Eigen::Matrix<double, 7, 1> endState = std::get<1>(res);
    const Eigen::Matrix<double, 5, 5> Q = std::get<2>(res);
    const Eigen::Matrix<double, 5, 9> FdF = std::get<3>(res);
    dEdxlast_ = std::get<4>(res);
    const Eigen::Matrix<double, 5, 5> dQMS = std::get<5>(res);
    const Eigen::Matrix<double, 5, 5> dQI = std::get<6>(res);

    for (int i = 0; i < 5; ++i) {
      for (int j = 0; j < 5; ++j) {
        F_[5 * i + j] = FdF(i, j);
        Q_[5 * i + j] = Q(i, j);
        dQMS_[5 * i + j] = dQMS(i, j);
        dQI_[5 * i + j] = dQI(i, j);
      }
    }

    // reference state expressed in the DetUnit local frame
    const GloballyPositioned<double>::GlobalPoint gp(endState[0], endState[1], endState[2]);
    const GloballyPositioned<double>::GlobalVector gv(endState[3], endState[4], endState[5]);
    const GloballyPositioned<double>::LocalPoint lp = base.toLocal(gp);
    const GloballyPositioned<double>::LocalVector lv = base.toLocal(gv);
    refp_ = gv.mag();
    // cos(lambda) = pt/p: needed offline to put the two curvilinear angle
    // noises on a common footing before using azimuthal isotropy
    refpt_ = std::sqrt(endState[3] * endState[3] + endState[4] * endState[4]);
    refqop_ = refp_ > 0. ? endState[6] / refp_ : 0.;
    refdxdz_ = lv.x() / lv.z();
    refdydz_ = lv.y() / lv.z();
    reflocx_ = lp.x();
    reflocy_ = lp.y();
    reflocz_ = lp.z();
    refglobr_ = std::sqrt(gp.x() * gp.x() + gp.y() * gp.y());
    refglobz_ = gp.z();

    // drain the per-step logs of THIS leg before any further propagate call
    for (auto const &ms : prop_->msStepLog()) {
      msmoliv_.push_back(ms.effZ);
      msmoliv_.push_back(ms.effA);
      msmoliv_.push_back(ms.xg);
      msmoliv_.push_back(ms.pGeV);
      msmoliv_.push_back(ms.beta);
      msmoliv_.push_back(ms.thp2);
      msmoliv_.push_back(ms.dOverX0);
      msmoliv_.push_back(ms.stepGroup);
    }
    for (auto const &us : prop_->ioniStepLog()) {
      ioniurbanv_.push_back(us.rec.regime);
      ioniurbanv_.push_back(us.rec.gsig2);
      ioniurbanv_.push_back(us.rec.a1);
      ioniurbanv_.push_back(us.rec.e1);
      ioniurbanv_.push_back(us.rec.a2);
      ioniurbanv_.push_back(us.rec.e2);
      ioniurbanv_.push_back(us.rec.a3);
      ioniurbanv_.push_back(us.rec.e0r);
      ioniurbanv_.push_back(us.rec.tmaxr);
      ioniurbanv_.push_back(us.rec.scaling);
      ioniurbanv_.push_back(us.cs);
    }
    for (auto const &st : prop_->stepTransportLog()) {
      for (int i = 0; i < 25; ++i) {
        stepjacc_.push_back(st.jacc[i]);
      }
      stepnms_.push_back(st.nMs);
      stepnioni_.push_back(st.nIoni);
    }

    edm::LogPrint("G4ePropagationExport")
        << "leg " << k << " detid " << detid_ << "  r=" << refglobr_ << " z=" << refglobz_ << " cm  p=" << refp_
        << " GeV  steps=" << stepnms_.size() << " ms=" << msmoliv_.size() / 8 << " ioni=" << ioniurbanv_.size() / 11;

    tree_->Fill();
    state = endState;
  }
}

DEFINE_FWK_MODULE(G4ePropagationExport);
