#include "TrackPropagation/Geant4e/interface/CvhMasterThread.h"
#include "TrackPropagation/Geant4e/interface/CvhMaster.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "G4PhysicalVolumeStore.hh"

// Slim port of SimG4Core/Application/OscarMTMasterThread, adapted for the
// CVH G4Error flow. The state-loop / signalling pattern is preserved
// verbatim (it's load-bearing for the master/worker handshake); what
// differs is what initG4() actually does -- see CvhMaster.

CvhMasterThread::CvhMasterThread(const edm::ParameterSet& iConfig)
    : m_pGeoFromDD4hep(iConfig.getParameter<bool>("g4GeometryDD4hepSource")),
      m_pUseMagneticField(iConfig.getParameter<bool>("UseMagneticField")),
      m_magFieldLabel(iConfig.getParameter<std::string>("MagneticFieldLabel")),
      m_masterThreadState(ThreadState::NotExist) {
  std::unique_lock<std::mutex> lk(m_threadMutex);

  edm::LogVerbatim("Geant4e") << "CvhMasterThread: creating master thread (DD4hep="
                              << m_pGeoFromDD4hep << ")";

  m_masterThread = std::thread([&]() {
    std::unique_lock<std::mutex> lk2(m_threadMutex);
    edm::LogVerbatim("Geant4e") << "CvhMasterThread: initializing CvhMaster on dedicated thread";

    // CvhMaster construction touches G4 singletons; must happen on this
    // dedicated thread so the G4MTRunManagerKernel binds correctly.
    m_cvhMaster = std::make_shared<CvhMaster>(iConfig);

    bool isG4Alive = false;
    while (true) {
      m_mainCanProceed = true;
      m_notifyMainCv.notify_one();

      m_masterCanProceed = false;
      m_notifyMasterCv.wait(lk2, [&] { return m_masterCanProceed; });

      if (m_masterThreadState == ThreadState::BeginRun) {
        edm::LogVerbatim("Geant4e") << "CvhMasterThread: BeginRun -> initG4";
        m_cvhMaster->initG4(m_pDDD, m_pDD4hep, m_pMF);
        isG4Alive = true;
      } else if (m_masterThreadState == ThreadState::EndRun) {
        edm::LogVerbatim("Geant4e") << "CvhMasterThread: EndRun -> stopG4";
        m_cvhMaster->stopG4();
        G4PhysicalVolumeStore::Clean();
        isG4Alive = false;
      } else if (m_masterThreadState == ThreadState::Destruct) {
        if (isG4Alive) {
          throw cms::Exception("LogicError")
              << "CvhMasterThread: Geant4 is still alive at Destruct; state must transition through EndRun first";
        }
        break;
      } else {
        throw cms::Exception("LogicError")
            << "CvhMasterThread: illegal master thread state " << static_cast<int>(m_masterThreadState);
      }
    }

    // Destruct on this thread to keep G4 singleton teardown on the
    // construction thread (mirrors OscarMTMasterThread's pattern).
    m_cvhMaster.reset();
    lk2.unlock();
    edm::LogVerbatim("Geant4e") << "CvhMasterThread: master thread finished";
  });

  // Wait for the master thread to finish constructing CvhMaster.
  m_mainCanProceed = false;
  m_notifyMainCv.wait(lk, [&]() { return m_mainCanProceed; });
  lk.unlock();
  edm::LogVerbatim("Geant4e") << "CvhMasterThread: constructed";
}

CvhMasterThread::~CvhMasterThread() {
  if (!m_stopped) {
    stopThread();
  }
}

void CvhMasterThread::callConsumes(edm::ConsumesCollector&& iC) const {
  if (m_hasToken) {
    return;
  }
  if (m_pGeoFromDD4hep) {
    m_DD4hep = iC.esConsumes<cms::DDCompactView, IdealGeometryRecord, edm::Transition::BeginRun>();
  } else {
    m_DDD = iC.esConsumes<DDCompactView, IdealGeometryRecord, edm::Transition::BeginRun>();
  }
  if (m_pUseMagneticField) {
    m_MagField = iC.esConsumes<MagneticField, IdealMagneticFieldRecord, edm::Transition::BeginRun>(
        edm::ESInputTag("", m_magFieldLabel));
  }
  m_hasToken = true;
}

void CvhMasterThread::ensureG4Started(const edm::EventSetup& iSetup) const {
  std::lock_guard<std::mutex> lk(m_protectMutex);

  // Job-scoped lazy start: the first call builds the G4 world; every later
  // call (subsequent runs) is a no-op. The ES products captured here stay in
  // use for the whole job -- valid only while their IOVs are unchanged, so
  // a mid-job conditions change is refused loudly below (G4 cannot be
  // rebuilt; silently keeping the stale field/geometry pointers would be
  // physics corruption or a dangling pointer).
  if (m_g4Started) {
    if (iSetup.get<IdealGeometryRecord>().cacheIdentifier() != m_geometryCacheId) {
      throw cms::Exception("Conditions")
          << "CvhMasterThread: the geometry payload changed at a run boundary, but the G4 world "
          << "was built from the previous one and cannot be rebuilt. Process one geometry IOV per job.";
    }
    if (m_pUseMagneticField &&
        iSetup.get<IdealMagneticFieldRecord>().cacheIdentifier() != m_magFieldCacheId) {
      // The record-level cache identifier rolls at every run boundary where
      // ANY product in IdealMagneticFieldRecord has a run-limited IOV (e.g.
      // AutoMagneticFieldESProducer following RunInfo). That alone does not
      // mean OUR field product changed: a file-based labelled field (e.g.
      // ScalarPot3D) is produced once and the EventSetup keeps returning the
      // same cached object. Re-fetch the product and accept the boundary iff
      // it is literally the same object the G4 world (master + workers) was
      // built with; only a genuinely rebuilt product is fatal, since G4
      // holds the raw pointer for the whole job and cannot be re-pointed.
      const MagneticField* mfNew = &iSetup.getData(m_MagField);
      if (mfNew != m_pMF) {
        throw cms::Exception("Conditions")
            << "CvhMasterThread: the magnetic-field product was rebuilt at a run boundary (e.g. a "
            << "magnet-current change picked up by AutoMagneticFieldESProducer), but the G4 world "
            << "was built with the previous field object and cannot be rebuilt. Process one field "
            << "IOV per job.";
      }
      m_magFieldCacheId = iSetup.get<IdealMagneticFieldRecord>().cacheIdentifier();
      edm::LogInfo("Geant4e") << "CvhMasterThread: magnetic-field record IOV rolled at a run "
                              << "boundary but the field product is unchanged; continuing.";
    }
    return;
  }

  std::unique_lock<std::mutex> lk2(m_threadMutex);

  if (m_pGeoFromDD4hep) {
    m_pDD4hep = &(*iSetup.getTransientHandle(m_DD4hep));
  } else {
    m_pDDD = &(*iSetup.getTransientHandle(m_DDD));
  }
  if (m_pUseMagneticField) {
    m_pMF = &iSetup.getData(m_MagField);
    m_magFieldCacheId = iSetup.get<IdealMagneticFieldRecord>().cacheIdentifier();
  }
  m_geometryCacheId = iSetup.get<IdealGeometryRecord>().cacheIdentifier();

  m_masterThreadState = ThreadState::BeginRun;
  m_masterCanProceed = true;
  m_mainCanProceed = false;
  m_notifyMasterCv.notify_one();
  m_notifyMainCv.wait(lk2, [&]() { return m_mainCanProceed; });
  m_g4Started = true;
}

void CvhMasterThread::stopG4() const {
  std::lock_guard<std::mutex> lk(m_protectMutex);
  if (!m_g4Started) {
    return;
  }
  std::unique_lock<std::mutex> lk2(m_threadMutex);
  m_masterThreadState = ThreadState::EndRun;
  m_mainCanProceed = false;
  m_masterCanProceed = true;
  m_notifyMasterCv.notify_one();
  m_notifyMainCv.wait(lk2, [&]() { return m_mainCanProceed; });
  m_g4Started = false;
}

void CvhMasterThread::stopThread() {
  if (m_stopped) {
    return;
  }
  edm::LogVerbatim("Geant4e") << "CvhMasterThread::stopThread";
  // The state loop requires the EndRun handshake before Destruct.
  stopG4();
  std::unique_lock<std::mutex> lk2(m_threadMutex);
  m_masterThreadState = ThreadState::Destruct;
  m_masterCanProceed = true;
  m_notifyMasterCv.notify_one();
  lk2.unlock();
  m_masterThread.join();
  m_stopped = true;
}
