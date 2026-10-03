// SimHitStateNtuplizer
//
// Ground-truth side of the "clean propagation test" (Josh Bendavid's
// suggestion, 2026-08-04): propagate ONE particle with a FIXED initial state
// many times through the full Geant4 simulation and record the true final 5D
// helix state on sensor surfaces, so that the analytic (characteristic
// function) prediction of the propagated-state PDF can be tested against
// Geant4 with no hits, no track fit, no FSR, no selection and no background.
//
// A PSimHit already carries exactly the CMSSW 5D local trajectory
// parameterization at the point where the track enters the sensitive silicon:
//   (q/p, dx/dz, dy/dz, x, y)   in the DetUnit local frame
// built from entryPoint(), thetaAtEntry()/phiAtEntry() and pabs(). So no
// custom G4 stepping action is needed -- we simply flatten the PSimHits of the
// primary track into a tree, one entry per (event, crossed module).
//
// The matching deterministic reference state at the SAME surfaces, together
// with the per-step Urban/Moliere material records, is produced by
// G4ePropagationExport (model side); the two are compared offline in
// characteristic-function space by calibration_studies/resolution/cf_propagation_test.py.
//
// Note the surface convention: the reference plane is the sensor ENTRY face,
// i.e. the DetUnit plane displaced by -thickness/2 along local z. entryLocalZ
// is written out so that this can be verified to be constant event-by-event
// rather than assumed.

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <vector>

class SimHitStateNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit SimHitStateNtuplizer(const edm::ParameterSet &);
  static void fillDescriptions(edm::ConfigurationDescriptions &);

private:
  void analyze(const edm::Event &, const edm::EventSetup &) override;

  std::vector<edm::EDGetTokenT<edm::PSimHitContainer>> simHitTokens_;
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;

  int pdgId_;
  unsigned int trackId_;

  TTree *tree_ = nullptr;

  unsigned int event_ = 0;
  unsigned int nhit_ = 0;
  // one entry per crossed module, ordered along the trajectory (by tof)
  std::vector<unsigned int> detid_;
  std::vector<float> tof_;
  std::vector<float> qop_;    // q/|p| [1/GeV]  -- local param 0
  std::vector<float> dxdz_;   // local param 1
  std::vector<float> dydz_;   // local param 2
  std::vector<float> locx_;   // local param 3 [cm]
  std::vector<float> locy_;   // local param 4 [cm]
  std::vector<float> locz_;   // entry local z [cm] (surface-convention check)
  std::vector<float> pabs_;   // [GeV]
  std::vector<float> eloss_;  // [GeV] deposited in this sensor
  std::vector<float> pathlen_;
  std::vector<float> globr_;
  std::vector<float> globz_;
  std::vector<unsigned short> proctype_;

public:
  ~SimHitStateNtuplizer() override = default;
};

SimHitStateNtuplizer::SimHitStateNtuplizer(const edm::ParameterSet &iConfig)
    : geomToken_(esConsumes()),
      pdgId_(iConfig.getParameter<int>("pdgId")),
      trackId_(iConfig.getParameter<unsigned int>("trackId")) {
  usesResource("TFileService");
  for (const auto &tag : iConfig.getParameter<std::vector<edm::InputTag>>("simHitTags")) {
    simHitTokens_.push_back(consumes<edm::PSimHitContainer>(tag));
  }

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("simstates", "true 5D local states at sensor entry faces");
  tree_->Branch("event", &event_);
  tree_->Branch("nhit", &nhit_);
  tree_->Branch("detid", &detid_);
  tree_->Branch("tof", &tof_);
  tree_->Branch("qop", &qop_);
  tree_->Branch("dxdz", &dxdz_);
  tree_->Branch("dydz", &dydz_);
  tree_->Branch("locx", &locx_);
  tree_->Branch("locy", &locy_);
  tree_->Branch("locz", &locz_);
  tree_->Branch("pabs", &pabs_);
  tree_->Branch("eloss", &eloss_);
  tree_->Branch("pathlen", &pathlen_);
  tree_->Branch("globr", &globr_);
  tree_->Branch("globz", &globz_);
  tree_->Branch("proctype", &proctype_);
}

void SimHitStateNtuplizer::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
  const TrackerGeometry *geom = &iSetup.getData(geomToken_);

  detid_.clear();
  tof_.clear();
  qop_.clear();
  dxdz_.clear();
  dydz_.clear();
  locx_.clear();
  locy_.clear();
  locz_.clear();
  pabs_.clear();
  eloss_.clear();
  pathlen_.clear();
  globr_.clear();
  globz_.clear();
  proctype_.clear();

  // collect the primary's hits from every tracker container, then order them
  // along the trajectory. Time of flight is monotonic along the track for a
  // primary and is unaffected by which container (Low/HighTof) the hit landed
  // in, so it is the robust ordering key.
  std::vector<const PSimHit *> hits;
  for (const auto &token : simHitTokens_) {
    edm::Handle<edm::PSimHitContainer> handle;
    iEvent.getByToken(token, handle);
    if (!handle.isValid()) {
      continue;
    }
    for (const auto &hit : *handle) {
      if (hit.trackId() != trackId_) {
        continue;
      }
      if (pdgId_ != 0 && std::abs(hit.particleType()) != std::abs(pdgId_)) {
        continue;
      }
      hits.push_back(&hit);
    }
  }
  std::sort(hits.begin(), hits.end(), [](const PSimHit *a, const PSimHit *b) { return a->tof() < b->tof(); });

  for (const PSimHit *hit : hits) {
    const float theta = hit->thetaAtEntry();
    const float phi = hit->phiAtEntry();
    const float tanTheta = std::tan(theta);
    // charge sign: PSimHit stores the PDG id; for leptons the charge is
    // -sign(pdg) (pdg 13 = mu-), for the rest we take sign(pdg).
    const int pdg = hit->particleType();
    const double charge = (std::abs(pdg) == 11 || std::abs(pdg) == 13 || std::abs(pdg) == 15)
                              ? (pdg > 0 ? -1. : 1.)
                              : (pdg > 0 ? 1. : -1.);

    const Local3DPoint entry = hit->entryPoint();
    const GeomDet *det = geom->idToDet(DetId(hit->detUnitId()));
    GlobalPoint gp(0., 0., 0.);
    if (det != nullptr) {
      gp = det->surface().toGlobal(entry);
    }

    detid_.push_back(hit->detUnitId());
    tof_.push_back(hit->tof());
    qop_.push_back(hit->pabs() > 0.f ? charge / hit->pabs() : 0.f);
    dxdz_.push_back(tanTheta * std::cos(phi));
    dydz_.push_back(tanTheta * std::sin(phi));
    locx_.push_back(entry.x());
    locy_.push_back(entry.y());
    locz_.push_back(entry.z());
    pabs_.push_back(hit->pabs());
    eloss_.push_back(hit->energyLoss());
    pathlen_.push_back((hit->exitPoint() - hit->entryPoint()).mag());
    globr_.push_back(gp.perp());
    globz_.push_back(gp.z());
    proctype_.push_back(hit->processType());
  }

  event_ = iEvent.id().event();
  nhit_ = detid_.size();
  tree_->Fill();
}

void SimHitStateNtuplizer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::vector<edm::InputTag>>("simHitTags", {});
  desc.add<int>("pdgId", 13);
  desc.add<unsigned int>("trackId", 1);
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(SimHitStateNtuplizer);
