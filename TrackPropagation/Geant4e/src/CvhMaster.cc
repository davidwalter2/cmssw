#include "TrackPropagation/Geant4e/interface/CvhMaster.h"

#include "SimG4Core/Geometry/interface/DDDWorld.h"
#include "SimG4Core/Geometry/interface/CustomUIsession.h"
#include "SimG4Core/Geometry/interface/SensitiveDetectorCatalog.h"
#include "SimG4Core/MagneticField/interface/FieldBuilder.h"
#include "SimG4Core/MagneticField/interface/CMSFieldManager.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "DetectorDescription/DDCMS/interface/DDCompactView.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "G4StateManager.hh"
#include "G4ApplicationState.hh"
#include "G4UImanager.hh"
#include "G4TransportationManager.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4RegionStore.hh"
#include "G4MTRunManagerKernel.hh"
#include "G4PhysListUtil.hh"
#include "G4LossTableManager.hh"

#include "TrackPropagation/Geant4e/interface/G4ErrorPhysicsListForCVH.h"

CvhMaster::CvhMaster(const edm::ParameterSet& p)
    : m_geoFromDD4hep(p.getParameter<bool>("g4GeometryDD4hepSource")),
      m_pUseMagneticField(p.getParameter<bool>("UseMagneticField")),
      m_pField(p.getParameter<edm::ParameterSet>("MagneticField")) {
  m_kernel = new G4MTRunManagerKernel();
  m_stateManager = G4StateManager::GetStateManager();
  m_UIsession = new CustomUIsession();
  G4UImanager::GetUIpointer()->SetCoutDestination(m_UIsession);
  G4UImanager::GetUIpointer()->SetMasterUIManager(true);
  G4PhysListUtil::InitialiseParameters();
  G4LossTableManager::Instance();
  m_catalog = std::make_unique<SensitiveDetectorCatalog>();
}

CvhMaster::~CvhMaster() {
  if (!m_runTerminated) {
    terminateRun();
  }
  delete m_UIsession;
  m_UIsession = nullptr;
  // m_kernel: deleted by G4StateManager / G4MTRunManager destruction order;
  // we do not delete it here to match RunManagerMT's pattern.
}

void CvhMaster::initG4(const DDCompactView* pDD,
                       const cms::DDCompactView* pDD4hep,
                       const MagneticField* pMF) {
  if (m_managerInitialized) {
    edm::LogWarning("Geant4e") << "CvhMaster::initG4 already done; ignoring";
    return;
  }
  edm::LogVerbatim("Geant4e") << "CvhMaster::initG4 begin (DD4hep=" << m_geoFromDD4hep
                              << ", UseMagneticField=" << m_pUseMagneticField << ")";

  // Build the geometry. DDDWorld populates the global G4 stores
  // (G4PhysicalVolumeStore, G4LogicalVolumeStore, G4RegionStore) -- this is
  // the side-effect previously delivered by `geopro`.
  // verb=0, cuts=false, pcut=false: CVH does not need region cuts.
  m_world = std::make_unique<DDDWorld>(pDD, pDD4hep, *m_catalog, /*verb=*/0,
                                       /*cuts=*/false, /*pcut=*/false);
  G4VPhysicalVolume* world = m_world->GetWorldVolume();
  m_kernel->DefineWorldVolume(world, /*topologyIsChanged=*/true);

  const G4PhysicalVolumeStore* pvs = G4PhysicalVolumeStore::GetInstance();
  const G4LogicalVolumeStore* lvs = G4LogicalVolumeStore::GetInstance();
  const G4RegionStore* rgs = G4RegionStore::GetInstance();
  edm::LogVerbatim("Geant4e") << "CvhMaster: " << pvs->size() << " physical volumes; "
                              << lvs->size() << " logical volumes; " << rgs->size() << " regions";

  // Master magnetic field. Workers set up their own per-thread field
  // manager attached to the same MagneticField object in CvhWorker.
  if (m_pUseMagneticField) {
    m_pMagField = pMF;
    sim::FieldBuilder fieldBuilder(m_pMagField, m_pField);
    CMSFieldManager* fieldManager = new CMSFieldManager();
    G4TransportationManager* tM = G4TransportationManager::GetTransportationManager();
    tM->SetFieldManager(fieldManager);
    fieldBuilder.build(fieldManager, tM->GetPropagatorInField());
  }

  // Register the master physics list (G4ErrorPhysicsListForCVH) and bring
  // the kernel through PreInit -> Init -> Idle so that the master's
  // per-particle ProcessManagers get transportation / eLoss / step-length /
  // B-field limit processes attached. Workers later copy these process
  // managers via G4WorkerThread::BuildGeometryAndPhysicsVector (called from
  // CvhWorker) -- without this master-side setup G4SteppingManager would
  // see NULL ProcessManager for mu- on workers and abort (Tracking0011).
  //
  // G4ErrorPhysicsListForCVH::ConstructProcess is MT-safe: it skips
  // particles whose ProcessManager already has a Transportation process,
  // so InitializeWorker re-invocations on worker threads do not
  // double-register.
  m_stateManager->SetNewState(G4State_PreInit);
  m_physicsList = new G4ErrorPhysicsListForCVH();
  m_kernel->SetPhysics(m_physicsList);

  m_stateManager->SetNewState(G4State_Init);
  m_kernel->InitializePhysics();
  m_kernel->SetUpDecayChannels();
  if (!m_kernel->RunInitialization()) {
    throw cms::Exception("Geant4e") << "CvhMaster::initG4: master G4 kernel RunInitialization failed";
  }
  m_physicsList->SetCutsWithDefault();

  // 10 mm step-length limit -- the existing single-thread propagator
  // applies this in ensureGeant4eIsInitilizedForCVH. Apply on the master
  // here; workers re-apply on their own G4UImanager in CvhWorker.
  G4UImanager::GetUIpointer()->ApplyCommand("/geant4e/limits/stepLength 10.0 mm");

  m_stateManager->SetNewState(G4State_GeomClosed);
  m_managerInitialized = true;
  edm::LogVerbatim("Geant4e") << "CvhMaster::initG4 done";
}

void CvhMaster::stopG4() {
  edm::LogVerbatim("Geant4e") << "CvhMaster::stopG4";
  G4GeometryManager::GetInstance()->OpenGeometry();
  m_stateManager->SetNewState(G4State_Quit);
  if (!m_runTerminated) {
    terminateRun();
  }
}

void CvhMaster::terminateRun() {
  if (!m_runTerminated && m_kernel != nullptr && m_managerInitialized) {
    m_kernel->RunTermination();
  }
  m_runTerminated = true;
}
