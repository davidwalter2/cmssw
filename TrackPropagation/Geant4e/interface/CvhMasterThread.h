#ifndef TrackPropagation_Geant4e_CvhMasterThread_h
#define TrackPropagation_Geant4e_CvhMasterThread_h

// CvhMasterThread: dedicated G4 master thread for the CVH G4Error flow.
//
// Slim analogue of SimG4Core/Application/OscarMTMasterThread. Owns a
// dedicated std::thread that constructs and drives a CvhMaster (the G4 master
// kernel + DDDWorld + master magnetic field). Held in an edm::GlobalCache by
// the CVH residual-maker; initializeGlobalCache spawns the master thread once
// per job, globalBeginRun forwards BeginRun signals, globalEndJob joins it.
//
// EventSetup reads happen on the framework's main thread (callConsumes
// registers the tokens during the residual-maker's ctor; beginRun() reads
// them in); the G4 work itself happens on the master thread, signalled via a
// pair of condition variables and a small state enum.

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "DetectorDescription/DDCMS/interface/DDCompactView.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include <memory>
#include <string>
#include <thread>
#include <mutex>
#include <condition_variable>

namespace edm {
  class EventSetup;
  class ConsumesCollector;
}  // namespace edm

class CvhMaster;

class CvhMasterThread {
public:
  explicit CvhMasterThread(const edm::ParameterSet&);
  ~CvhMasterThread();

  // Called from the residual-maker's ctor (which receives the GlobalCache
  // pointer) to register the ES tokens this master needs. Idempotent.
  void callConsumes(edm::ConsumesCollector&& iC) const;

  void beginRun(const edm::EventSetup&) const;
  void endRun() const;
  void stopThread();

  inline CvhMaster& cvhMaster() const { return *m_cvhMaster; }
  inline CvhMaster* cvhMasterPtr() const { return m_cvhMaster.get(); }

private:
  enum class ThreadState { NotExist = 0, BeginRun = 1, EndRun = 2, Destruct = 3 };

  const bool m_pGeoFromDD4hep;
  const bool m_pUseMagneticField;
  const std::string m_magFieldLabel;

  std::shared_ptr<CvhMaster> m_cvhMaster;
  std::thread m_masterThread;

  mutable const DDCompactView* m_pDDD{nullptr};
  mutable const cms::DDCompactView* m_pDD4hep{nullptr};
  mutable const MagneticField* m_pMF{nullptr};
  mutable edm::ESGetToken<DDCompactView, IdealGeometryRecord> m_DDD;
  mutable edm::ESGetToken<cms::DDCompactView, IdealGeometryRecord> m_DD4hep;
  mutable edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> m_MagField;

  mutable std::mutex m_protectMutex;
  mutable std::mutex m_threadMutex;
  mutable std::condition_variable m_notifyMasterCv;
  mutable std::condition_variable m_notifyMainCv;

  mutable ThreadState m_masterThreadState;

  mutable bool m_hasToken{false};
  mutable bool m_masterCanProceed{false};
  mutable bool m_mainCanProceed{false};
  mutable bool m_firstRun{true};
  mutable bool m_stopped{false};
};

#endif
