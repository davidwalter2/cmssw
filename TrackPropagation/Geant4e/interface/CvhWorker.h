#ifndef TrackPropagation_Geant4e_CvhWorker_h
#define TrackPropagation_Geant4e_CvhWorker_h

// CvhWorker: per-stream / per-thread G4 setup for the CVH G4Error flow.
//
// Owned by ResidualGlobalCorrectionMakerBase (one per stream). ensureInitialized
// is idempotent within a thread (gated by a static thread_local flag): the
// first produce() call on each TBB worker thread does the work
// (G4WorkerThread::BuildGeometryAndPhysicsVector +
//  G4WorkerRunManagerKernel + WorkerDefineWorldVolume +
//  G4TransportationManager::SetWorldForTracking + per-thread magnetic field
//  via sim::FieldBuilder + CMSFieldManager); subsequent calls are a TLS-load
// no-op.
//
// After this runs, Geant4ePropagator::ensureGeant4eIsInitilizedForCVH layers
// G4Error physics on top -- the navigator's world is now set so the chained
// G4SafetyHelper::InitialiseHelper / G4WentzelVIModel::Initialise no longer
// abort.

class CvhMaster;

class CvhWorker {
public:
  CvhWorker() = default;
  ~CvhWorker() = default;

  CvhWorker(const CvhWorker&) = delete;
  CvhWorker& operator=(const CvhWorker&) = delete;

  // Idempotent per-thread G4 setup. Safe to call from every produce(); only
  // the first call on each TBB worker thread does work.
  void ensureInitialized(CvhMaster& master);
};

#endif
