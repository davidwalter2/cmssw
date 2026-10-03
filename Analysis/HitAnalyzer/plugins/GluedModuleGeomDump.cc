// One-shot dump of surface centers + normals for a configurable list of
// tracker detids, plus their glued-partner relationships. Diagnostic for
// the CVH runaway-leg localization study: the failing propagation targets
// cluster on specific glued-module faces whose center separation looks
// anomalous in the trace prints; this dumps the same quantities straight
// from the TrackerGeometry (aligned or ideal depending on the config).

#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GluedGeomDet.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "CalibFormats/SiStripObjects/interface/SiStripQuality.h"
#include "CalibTracker/Records/interface/SiStripDependentRecords.h"

#include <map>

class GluedModuleGeomDump : public edm::one::EDAnalyzer<> {
public:
  explicit GluedModuleGeomDump(const edm::ParameterSet& cfg)
      : geomToken_(esConsumes()),
        topoToken_(esConsumes()),
        detids_(cfg.getParameter<std::vector<unsigned int>>("detIds")),
        scanAllGlued_(cfg.getParameter<bool>("scanAllGlued")),
        useQuality_(cfg.getParameter<bool>("useQuality")),
        compareIdeal_(cfg.getParameter<bool>("compareIdeal")) {
    if (useQuality_) {
      qualityToken_ = esConsumes<SiStripQuality, SiStripQualityRcd>(edm::ESInputTag(""));
    }
    if (compareIdeal_) {
      idealGeomToken_ = esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>(edm::ESInputTag("", "idealForDigi"));
    }
  }

private:
  void analyze(const edm::Event&, const edm::EventSetup& es) override {
    const TrackerGeometry& geom = es.getData(geomToken_);
    const TrackerTopology& topo = es.getData(topoToken_);
    const SiStripQuality* quality = useQuality_ ? &es.getData(qualityToken_) : nullptr;

    for (unsigned int rawid : detids_) {
      const DetId id(rawid);
      const GeomDet* det = geom.idToDet(id);
      if (det == nullptr) {
        std::cout << "detid " << rawid << " NOT FOUND" << std::endl;
        continue;
      }
      const auto& pos = det->surface().position();
      const auto normal = det->surface().toGlobal(LocalVector(0., 0., 1.));
      std::cout << "detid " << rawid
                << "  glued=" << topo.glued(id)
                << "  isStereo=" << (topo.isStereo(id) ? 1 : 0)
                << "  pos(r,phi,z)=(" << pos.perp() << ", " << pos.phi() << ", " << pos.z() << ")"
                << "  normal=(" << normal.x() << ", " << normal.y() << ", " << normal.z() << ")";
      if (quality != nullptr && id.subdetId() >= 3) {  // strip subdetectors only
        std::cout << "  IsModuleBad=" << quality->IsModuleBad(rawid);
      }
      std::cout << std::endl;
    }

    if (useQuality_) {
      // Global bad-module census from the conditions (merged quality:
      // bad channels + fibers + modules + DCS/RunInfo for this run).
      unsigned long nBadModules = 0;
      for (auto const* du : geom.detUnits()) {
        const DetId id = du->geographicalId();
        if (id.subdetId() >= 3 && quality->IsModuleBad(id.rawId())) {
          ++nBadModules;
        }
      }
      std::cout << "SiStripQuality: total bad strip modules this run = " << nBadModules << std::endl;
    }

    if (scanAllGlued_) {
      // Sweep all glued pairs: flag pairs whose aligned faces are not
      // parallel (|sin(angle)| between normals beyond threshold). A healthy
      // glued module has exactly (anti)parallel faces; a large relative
      // rotation is an alignment-payload artifact (unconstrained module).
      std::map<uint32_t, std::pair<const GeomDet*, const GeomDet*>> pairs;
      for (auto const* du : geom.detUnits()) {
        const DetId id = du->geographicalId();
        const uint32_t glued = topo.glued(id);
        if (glued == 0) continue;
        auto& pr = pairs[glued];
        if (topo.isStereo(id)) pr.first = du;
        else pr.second = du;
      }
      unsigned int nflagged = 0;
      for (auto const& [glued, pr] : pairs) {
        if (pr.first == nullptr || pr.second == nullptr) continue;
        const auto n1 = pr.first->surface().toGlobal(LocalVector(0., 0., 1.));
        const auto n2 = pr.second->surface().toGlobal(LocalVector(0., 0., 1.));
        const double cosang = n1.x()*n2.x() + n1.y()*n2.y() + n1.z()*n2.z();
        const double sinang = std::sqrt(std::max(0., 1. - cosang*cosang));
        if (sinang > 0.01) {  // > 10 mrad relative face rotation
          ++nflagged;
          const auto& p1 = pr.first->surface().position();
          const auto& p2 = pr.second->surface().position();
          std::cout << "NONPARALLEL glued=" << glued
                    << "  stereo=" << pr.first->geographicalId().rawId()
                    << "  mono=" << pr.second->geographicalId().rawId()
                    << "  sin(angle)=" << sinang
                    << "  stereo(r,phi,z)=(" << p1.perp() << "," << p1.phi() << "," << p1.z() << ")"
                    << "  mono(r,phi,z)=(" << p2.perp() << "," << p2.phi() << "," << p2.z() << ")";
          if (quality != nullptr) {
            std::cout << "  badStereo=" << quality->IsModuleBad(pr.first->geographicalId().rawId())
                      << "  badMono=" << quality->IsModuleBad(pr.second->geographicalId().rawId());
          }
          std::cout << std::endl;
        }
      }
      std::cout << "glued pairs scanned = " << pairs.size()
                << "  non-parallel (>10 mrad) = " << nflagged << std::endl;
    }

    if (compareIdeal_) {
      // Absolute aligned-vs-ideal sweep over EVERY tracker sensor (pixels
      // and strips, single- and double-sided): flag modules whose aligned
      // orientation or position departs from the design far beyond any
      // physical alignment correction (normals > 20 mrad or centers > 5 mm).
      // This catches unconstrained-alignment ("garbage") entries also for
      // single-sided modules, which the glued-pair parallelism check above
      // cannot see.
      const TrackerGeometry& ideal = es.getData(idealGeomToken_);
      unsigned int nchecked = 0, nflagged = 0;
      for (auto const* du : geom.detUnits()) {
        const DetId id = du->geographicalId();
        const GeomDet* dui = ideal.idToDet(id);
        if (dui == nullptr) continue;
        ++nchecked;
        const auto na = du->surface().toGlobal(LocalVector(0., 0., 1.));
        const auto ni = dui->surface().toGlobal(LocalVector(0., 0., 1.));
        const double cosang = std::abs(na.x()*ni.x() + na.y()*ni.y() + na.z()*ni.z());
        const double sinang = std::sqrt(std::max(0., 1. - cosang*cosang));
        const auto& pa = du->surface().position();
        const auto& pi = dui->surface().position();
        const double dx = pa.x()-pi.x(), dy = pa.y()-pi.y(), dz = pa.z()-pi.z();
        const double dpos = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (sinang > 0.02 || dpos > 0.5) {
          ++nflagged;
          std::cout << "OUTLIER detid=" << id.rawId()
                    << "  subdet=" << id.subdetId()
                    << "  glued=" << topo.glued(id)
                    << "  sin(tilt)=" << sinang
                    << "  dpos=" << dpos << " cm"
                    << "  aligned(r,phi,z)=(" << pa.perp() << "," << pa.phi() << "," << pa.z() << ")"
                    << "  ideal(r,phi,z)=(" << pi.perp() << "," << pi.phi() << "," << pi.z() << ")";
          if (quality != nullptr && id.subdetId() >= 3) {
            std::cout << "  IsModuleBad=" << quality->IsModuleBad(id.rawId());
          }
          std::cout << std::endl;
        }
      }
      std::cout << "aligned-vs-ideal sweep: sensors checked = " << nchecked
                << "  outliers (tilt > 20 mrad or |dpos| > 5 mm) = " << nflagged << std::endl;
    }
  }

  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> idealGeomToken_;
  edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
  edm::ESGetToken<SiStripQuality, SiStripQualityRcd> qualityToken_;
  std::vector<unsigned int> detids_;
  bool scanAllGlued_;
  bool useQuality_;
  bool compareIdeal_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GluedModuleGeomDump);
