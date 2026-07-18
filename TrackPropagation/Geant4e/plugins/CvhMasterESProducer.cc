// CvhMasterESProducer: builds the single, shared CVH Geant4 master
// (CvhMasterThread) as an EventSetup product on CvhMasterRecord.
//
// The master owns a dedicated G4 thread + G4MTRunManagerKernel (a Geant4
// process-global singleton), so exactly one may exist per job. Producing it
// here -- and consuming it via esConsumes in every CVH residual maker -- lets
// any number/mix of maker instances share the one master, replacing the old
// per-producer edm::GlobalCache<CvhMasterThread> (which is per-module-label
// and collides on the second instance).
//
// The geometry (IdealGeometryRecord) and labelled master field
// (IdealMagneticFieldRecord) are consumed from the dependent record and handed
// to CvhMasterThread::ensureG4Started, which starts G4 on the master thread on
// the first produce() and is a guarded no-op thereafter. Both dependencies are
// single-IOV over a normal job, so the master is built exactly once.

#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Utilities/interface/ESInputTag.h"

#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "DetectorDescription/DDCMS/interface/DDCompactView.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "TrackPropagation/Geant4e/interface/CvhMasterThread.h"
#include "TrackPropagation/Geant4e/interface/CvhMasterRecord.h"

#include <memory>
#include <string>

class CvhMasterESProducer : public edm::ESProducer {
public:
  explicit CvhMasterESProducer(const edm::ParameterSet& iConfig)
      : masterPSet_(iConfig),
        geoFromDD4hep_(iConfig.getParameter<bool>("g4GeometryDD4hepSource")),
        useMagneticField_(iConfig.getParameter<bool>("UseMagneticField")) {
    auto cc = setWhatProduced(this);
    if (geoFromDD4hep_) {
      dd4hepToken_ = cc.consumes();
    } else {
      dddToken_ = cc.consumes();
    }
    if (useMagneticField_) {
      mfToken_ = cc.consumes(edm::ESInputTag("", iConfig.getParameter<std::string>("MagneticFieldLabel")));
    }
  }

  std::shared_ptr<CvhMasterThread> produce(const CvhMasterRecord& iRecord) {
    // Build once; the framework caches the returned shared_ptr and (for the
    // normal single-IOV job) calls produce() exactly once.
    if (!master_) {
      master_ = std::make_shared<CvhMasterThread>(masterPSet_);
    }
    const DDCompactView* ddd = nullptr;
    const cms::DDCompactView* dd4hep = nullptr;
    const MagneticField* mf = nullptr;
    if (geoFromDD4hep_) {
      dd4hep = iRecord.getTransientHandle(dd4hepToken_).product();
    } else {
      ddd = iRecord.getTransientHandle(dddToken_).product();
    }
    if (useMagneticField_) {
      mf = &iRecord.get(mfToken_);
    }
    master_->ensureG4Started(ddd, dd4hep, mf);
    return master_;
  }

private:
  const edm::ParameterSet masterPSet_;
  const bool geoFromDD4hep_;
  const bool useMagneticField_;
  edm::ESGetToken<DDCompactView, IdealGeometryRecord> dddToken_;
  edm::ESGetToken<cms::DDCompactView, IdealGeometryRecord> dd4hepToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> mfToken_;
  std::shared_ptr<CvhMasterThread> master_;
};

DEFINE_FWK_EVENTSETUP_MODULE(CvhMasterESProducer);
