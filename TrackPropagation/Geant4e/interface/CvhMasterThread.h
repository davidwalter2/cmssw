#ifndef TrackPropagation_Geant4e_CvhMasterThread_h
#define TrackPropagation_Geant4e_CvhMasterThread_h

// CvhMasterThread: dedicated G4 master thread for the CVH G4Error flow.
//
// Slim analogue of SimG4Core/Application/OscarMTMasterThread. Owns a
// dedicated std::thread that constructs and drives a CvhMaster (the G4 master
// kernel + DDDWorld + master magnetic field).
//
// Shared across all CVH residual makers in a job as an EventSetup product on
// CvhMasterRecord (produced by CvhMasterESProducer). The ESProducer builds one
// instance per job, fetches the geometry + labelled field from the EventSetup,
// and calls ensureG4Started(products) (first call starts G4, later calls are
// no-ops guarded against a genuine product change). The master thread is joined
// by ~CvhMasterThread (via stopThread) when the ESProducer releases it at end
// of job. The G4 work itself happens on the master thread, signalled via a
// pair of condition variables and a small state enum.

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "DetectorDescription/DDCMS/interface/DDCompactView.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include <memory>
#include <string>
#include <thread>
#include <mutex>
#include <condition_variable>

class CvhMaster;

class CvhMasterThread {
public:
  explicit CvhMasterThread(const edm::ParameterSet&);
  ~CvhMasterThread();

  // The G4 world is a JOB-scoped resource with a lazy start: it needs the
  // geometry + field, so it cannot be built at construction. The
  // CvhMasterESProducer fetches those from the EventSetup and passes them
  // here; the first call initializes G4 exactly once, later calls are no-ops
  // (and refuse a genuinely changed geometry/field product, since G4 holds the
  // raw pointers for the whole job and cannot be rebuilt). Pass the DDD or the
  // DD4hep view (whichever matches g4GeometryDD4hepSource); the other is null.
  // stopThread() tears G4 down at end of job.
  void ensureG4Started(const DDCompactView* ddd,
                       const cms::DDCompactView* dd4hep,
                       const MagneticField* mf) const;
  void stopThread();

  inline CvhMaster& cvhMaster() const { return *m_cvhMaster; }
  inline CvhMaster* cvhMasterPtr() const { return m_cvhMaster.get(); }

private:
  enum class ThreadState { NotExist = 0, BeginRun = 1, EndRun = 2, Destruct = 3 };

  // Signals the state loop to tear down G4 (EndRun handshake). Only used
  // from stopThread(); no-op if G4 was never started.
  void stopG4() const;

  const bool m_pGeoFromDD4hep;
  const bool m_pUseMagneticField;
  const std::string m_magFieldLabel;

  std::shared_ptr<CvhMaster> m_cvhMaster;
  std::thread m_masterThread;

  // The geometry/field products the G4 world captured at first start; kept to
  // detect (and refuse) a mid-job change -- G4 cannot be rebuilt, so stale
  // conditions must fail loudly.
  mutable const DDCompactView* m_pDDD{nullptr};
  mutable const cms::DDCompactView* m_pDD4hep{nullptr};
  mutable const MagneticField* m_pMF{nullptr};

  mutable std::mutex m_protectMutex;
  mutable std::mutex m_threadMutex;
  mutable std::condition_variable m_notifyMasterCv;
  mutable std::condition_variable m_notifyMainCv;

  mutable ThreadState m_masterThreadState;

  mutable bool m_masterCanProceed{false};
  mutable bool m_mainCanProceed{false};
  mutable bool m_g4Started{false};
  mutable bool m_stopped{false};
};

#endif
