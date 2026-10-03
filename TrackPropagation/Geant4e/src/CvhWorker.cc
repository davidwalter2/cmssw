#include "TrackPropagation/Geant4e/interface/CvhWorker.h"
#include "TrackPropagation/Geant4e/interface/CvhMaster.h"

#include "SimG4Core/Geometry/interface/DDDWorld.h"
#include "SimG4Core/MagneticField/interface/FieldBuilder.h"
#include "SimG4Core/MagneticField/interface/CMSFieldManager.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "G4Threading.hh"
#include "G4WorkerThread.hh"
#include "G4WorkerRunManagerKernel.hh"
#include "G4TransportationManager.hh"
#include "G4StateManager.hh"
#include "G4UImanager.hh"
#include "G4VUserPhysicsList.hh"
#include "G4ErrorPropagatorData.hh"
#include "G4ErrorPropagatorManager.hh"

#include <atomic>
#include <mutex>

namespace {
  // Single counter so each TBB worker thread gets a unique G4 thread id
  // (matches RunManagerMTWorker's `thread_counter` pattern).
  std::atomic<int> g_threadCounter{0};
  int newThreadIndex() { return g_threadCounter.fetch_add(1); }

  // Serialise the per-thread G4 setup. G4WorkerThread::BuildGeometryAndPhysics
  // Vector and G4WorkerRunManagerKernel construction touch shared G4 stores
  // even though the bulk of their state ends up in G4ThreadLocal; concurrent
  // first-time invocations on different threads have historically caused
  // store corruption (the same reason the existing geant4eInitMutex pattern
  // exists in Geant4ePropagator.cc).
  std::mutex& cvhWorkerInitMutex() {
    static std::mutex m;
    return m;
  }
}  // namespace

void CvhWorker::ensureInitialized(CvhMaster& master) {
  // Per-thread "did we run" flag. Each TBB worker thread runs this block
  // exactly once for its lifetime; subsequent produce()s skip directly past.
  static thread_local bool initializedOnThisThread = false;
  if (initializedOnThisThread) {
    return;
  }
  std::lock_guard<std::mutex> lk(cvhWorkerInitMutex());
  // Re-check under lock (another thread won the race -- not possible since
  // flag is per-thread, but cheap to defend).
  if (initializedOnThisThread) {
    return;
  }

  const int thisID = newThreadIndex();
  G4Threading::G4SetThreadId(thisID);
  edm::LogVerbatim("Geant4e") << "CvhWorker::ensureInitialized on thread " << thisID;

  // G4 11 stores per-thread workspace (G4ThreadLocal arrays for physical
  // volume rotation/translation indices, navigator safety helper, etc.)
  // that must be allocated before any geometry method is called on this
  // thread. Without this, G4VPhysicalVolume::GetTranslation segfaults
  // when SetWorldForTracking dereferences the world.
  G4WorkerThread::BuildGeometryAndPhysicsVector();

  auto* existing = G4WorkerRunManagerKernel::GetRunManagerKernel();
  G4WorkerRunManagerKernel* kernel =
      (existing != nullptr) ? dynamic_cast<G4WorkerRunManagerKernel*>(existing)
                            : new G4WorkerRunManagerKernel();
  if (kernel == nullptr) {
    throw cms::Exception("Geant4e")
        << "CvhWorker: existing G4 worker kernel is not a G4WorkerRunManagerKernel; "
        << "another G4 client on this thread has registered an incompatible kernel.";
  }

  G4VPhysicalVolume* worldPV = master.world().GetWorldVolume();
  kernel->WorkerDefineWorldVolume(worldPV);
  G4TransportationManager* tM = G4TransportationManager::GetTransportationManager();
  tM->SetWorldForTracking(worldPV);

  // Per-thread magnetic field. The MagneticField object itself is shared
  // (read-only); each thread gets its own CMSFieldManager + stepper.
  if (master.useMagneticField() && master.magneticField() != nullptr) {
    sim::FieldBuilder fieldBuilder(master.magneticField(), master.fieldPSet());
    CMSFieldManager* fieldManager = new CMSFieldManager();
    tM->SetFieldManager(fieldManager);
    fieldBuilder.build(fieldManager, tM->GetPropagatorInField());
  }

  // Share the master's G4ErrorPhysicsListForCVH with this worker.
  // InitializeWorker rebuilds the per-thread physics tables; the
  // physics-list ConstructProcess is one-shot per worker thread (guarded
  // by a thread_local flag) so the worker side does NOT double-register
  // processes.
  G4StateManager::GetStateManager()->SetNewState(G4State_Init);
  G4VUserPhysicsList* physicsList = master.physicsListForWorker();
  physicsList->InitializeWorker();
  kernel->SetPhysics(physicsList);
  kernel->InitializePhysics();
  if (!kernel->RunInitialization()) {
    throw cms::Exception("Geant4e")
        << "CvhWorker: G4 worker kernel RunInitialization failed on thread " << thisID;
  }

  // Re-apply the step-length limit on this thread's G4UImanager (the
  // master only set it on its own G4UI; UImanager is G4ThreadLocal).
  G4UImanager::GetUIpointer()->ApplyCommand("/geant4e/limits/stepLength 10.0 mm");

  // Leave G4 in PreInit so the propagator's first-call init (which calls
  // ensureGeant4eIsInitilizedForCVH(true)) is allowed to drive
  // G4ErrorRunManagerHelper through its own setup -- which is what
  // allocates G4ErrorPropagator (theG4ErrorPropagator) and its internal
  // navigator. Without that allocation InitTrackPropagation segfaults.
  // Process re-registration cannot happen: G4ErrorPhysicsListForCVH
  // ConstructProcess is one-shot per worker thread (thread_local flag),
  // and all physics-list instances in the job share the same particle
  // set (job-wide registry in G4ErrorPhysicsListForCVH), so the skipped
  // invocations concern exactly the particles already configured.
  G4StateManager::GetStateManager()->SetNewState(G4State_PreInit);

  initializedOnThisThread = true;
  edm::LogVerbatim("Geant4e") << "CvhWorker::ensureInitialized done on thread " << thisID;
}
