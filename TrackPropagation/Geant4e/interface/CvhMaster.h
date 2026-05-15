#ifndef TrackPropagation_Geant4e_CvhMaster_h
#define TrackPropagation_Geant4e_CvhMaster_h

// CvhMaster: G4MT master-side runtime for the CVH G4Error propagator.
//
// Slim analogue of SimG4Core/Application/RunManagerMT, restricted to what the
// CVH residual-makers need: a master G4MTRunManagerKernel, a DDDWorld
// (populates G4PhysicalVolumeStore once for the whole job), and the master
// magnetic field. No sensitive detectors, no run actions, no event-generator
// hooks. Constructed on a dedicated thread by CvhMasterThread; workers
// (CvhWorker) attach to it via DDDWorld::GetWorldVolume() + the
// master-supplied magnetic field.
//
// The CVH-specific physics (G4ErrorPhysicsListForCVH) is *not* registered
// here. The existing Geant4ePropagator::ensureGeant4eIsInitilizedForCVH
// handles that per-thread on the first propagate() call, layered on top of
// the master/worker G4 setup populated below. The mismatch is why pre-MT
// runs aborted on the first event in worker threads.

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <memory>

namespace cms {
  class DDCompactView;
}
class DDCompactView;
class MagneticField;
class G4StateManager;
class G4MTRunManagerKernel;
class DDDWorld;
class CustomUIsession;
class SensitiveDetectorCatalog;

class G4VUserPhysicsList;

class CvhMaster {
public:
  explicit CvhMaster(const edm::ParameterSet&);
  ~CvhMaster();

  // Called from CvhMasterThread's state loop on the dedicated G4 master
  // thread. Builds DDDWorld + master G4MTRunManagerKernel + master magnetic
  // field + master G4ErrorPhysicsListForCVH; brings the kernel through
  // PreInit -> Init -> Idle so the master's per-particle ProcessManagers
  // have transportation / eLoss / step-length-limit processes attached.
  // Workers copy this state via G4WorkerThread::BuildGeometryAndPhysicsVector.
  void initG4(const DDCompactView* pDD, const cms::DDCompactView* pDD4hep, const MagneticField* pMF);

  // Symmetric: tear down so the master thread can join cleanly.
  void stopG4();

  inline const DDDWorld& world() const { return *m_world; }
  inline SensitiveDetectorCatalog& catalog() { return *m_catalog; }
  inline const MagneticField* magneticField() const { return m_pMagField; }
  inline const edm::ParameterSet& fieldPSet() const { return m_pField; }
  inline bool useMagneticField() const { return m_pUseMagneticField; }
  inline bool geoFromDD4hep() const { return m_geoFromDD4hep; }
  // Non-const accessor: workers call InitializeWorker on this pointer.
  // The physics list is owned by m_kernel (after SetPhysics).
  inline G4VUserPhysicsList* physicsListForWorker() const { return m_physicsList; }

private:
  void terminateRun();

  G4MTRunManagerKernel* m_kernel{nullptr};
  G4StateManager* m_stateManager{nullptr};
  CustomUIsession* m_UIsession{nullptr};
  std::unique_ptr<DDDWorld> m_world;
  std::unique_ptr<SensitiveDetectorCatalog> m_catalog;
  G4VUserPhysicsList* m_physicsList{nullptr};

  const bool m_geoFromDD4hep;
  const bool m_pUseMagneticField;
  edm::ParameterSet m_pField;

  const MagneticField* m_pMagField{nullptr};
  bool m_managerInitialized{false};
  bool m_runTerminated{false};
};

#endif
