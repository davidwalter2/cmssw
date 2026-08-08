#include <algorithm>
#include <sstream>

// Geant4e
#include "TrackPropagation/Geant4e/interface/ConvertFromToCLHEP.h"
#include "TrackPropagation/Geant4e/interface/Geant4ePropagator.h"
#include "TrackPropagation/Geant4e/interface/G4ErrorPhysicsListForCVH.h"
#include "TrackPropagation/Geant4e/interface/MaterialGroupModel.h"

// CMSSW
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/TrajectoryState/interface/SurfaceSideDefinition.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/AnalyticalJacobians/interface/AnalyticalCurvilinearJacobian.h"

// Geant4
#include "G4Box.hh"
#include "G4ErrorCylSurfaceTarget.hh"
#include "G4ErrorFreeTrajState.hh"
#include "G4ErrorPlaneSurfaceTarget.hh"
#include "G4ErrorPropagatorData.hh"
#include "G4ErrorRunManagerHelper.hh"
#include "G4EventManager.hh"
#include "G4Field.hh"
#include "G4FieldManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4SteppingControl.hh"
#include "G4TransportationManager.hh"
#include "G4Tubs.hh"
#include "G4UImanager.hh"
#include "G4PathFinder.hh"
#include "G4ErrorPropagationNavigator.hh"
#include "G4RunManagerKernel.hh"
#include "G4StateManager.hh"

// CLHEP
#include <CLHEP/Units/SystemOfUnits.h>

#include "TrackingTools/AnalyticalJacobians/interface/JacobianCurvilinearToLocal.h"
#include "TrackingTools/TrajectoryParametrization/interface/CurvilinearTrajectoryParameters.h"

#include <Eigen/Geometry>

#include "SimG4Core/MagneticField/interface/Field.h"

#include "Math/ProbFuncMathCore.h"
#include "Math/QuantFuncMathCore.h"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4MuBremsstrahlungModel.hh"
#include "G4MuPairProductionModel.hh"
#include "G4DataVector.hh"
#include "G4Element.hh"

#include <mutex>

namespace {
  // Process-wide gate for first-time Geant4e initialisation per thread.
  // G4 11.1 makes G4ErrorPropagatorManager / G4ErrorPropagatorData
  // G4ThreadLocal, but constructing the TLS instance still mutates the
  // global G4 geometry store (G4RegionStore, world volume, etc.). Two TBB
  // workers landing in InitGeant4e at the same time corrupt that store —
  // so we serialise the first-call init while leaving subsequent (cheap,
  // TLS-pointer-only) accesses unlocked. This is single-acquire-per-thread.
  std::mutex& geant4eInitMutex() {
    static std::mutex m;
    return m;
  }
  bool& geant4eInitDoneForThread() {
    static thread_local bool done = false;
    return done;
  }
}

/** Constructor.
 */
Geant4ePropagator::Geant4ePropagator(const MagneticField *field,
                                     std::string particleName,
                                     PropagationDirection dir,
                                     double plimit,
                                     bool forCVH,
                                     double ioniTruncAlpha)
    : Propagator(dir),
      theField(field),
      theParticleName(particleName),
      plimit_(plimit),
      forCVH_(forCVH),
      ioniTruncAlpha_(ioniTruncAlpha) {
  LogDebug("Geant4e") << "Geant4e Propagator initialized";

  // G4 init is deferred: in MT mode the ESProducer's produce() runs eagerly
  // on a TBB worker, before any EDProducer (in particular `geopro`, which
  // populates G4's world volume) has run. Calling InitGeant4e here aborts
  // with "No world defined in your geometry". The first propagate() call
  // (which can only happen after geopro has run for the current event)
  // does the init then under geant4eInitMutex(). fluct allocation is
  // deferred too -- the G4UniversalFluctuationForExtrapolator ctor chains
  // into G4TablesForExtrapolatorForCVH -> G4WentzelVIModel::Initialise ->
  // G4SafetyHelper::InitialiseHelper, which aborts on a worker thread that
  // hasn't had its world set up yet. fluct is allocated lazily in the
  // first-call init block of propagateGenericWithJacobianAltD, AFTER the
  // CvhWorker has installed the world on this thread and
  // ensureGeant4eIsInitilizedForCVH has registered the G4Error physics.
}

/** Deep-copy constructor.
 *
 * The default copy ctor would shallow-copy `fluct` and let two propagator
 * instances share the same G4UniversalFluctuationForExtrapolator — and then
 * the dtor of one would delete it from under the other. clone() (used by
 * the residual-makers in Tier-3 MT) goes through this ctor, so each per-
 * stream copy gets its own fluct, configured for the same particle.
 */
Geant4ePropagator::Geant4ePropagator(const Geant4ePropagator &other)
    : Propagator(other.propagationDirection()),
      theField(other.theField),
      theParticleName(other.theParticleName),
      plimit_(other.plimit_),
      forCVH_(other.forCVH_),
      ioniTruncAlpha_(other.ioniTruncAlpha_),
      ioniStepLogging_(other.ioniStepLogging_),
      stepTransportLogging_(other.stepTransportLogging_) {
  // fluct allocation is deferred to the first propagate() call on this
  // thread (under geant4eInitMutex), AFTER the per-thread G4 world has
  // been set up by CvhWorker. Allocating it eagerly in the deep-copy ctor
  // aborts on a worker thread that hasn't yet had its navigator initialised.
}

/** Destructor.
 */
Geant4ePropagator::~Geant4ePropagator() {
  LogDebug("Geant4e") << "Geant4ePropagator::~Geant4ePropagator()" << std::endl;
  // One-line summary of propagateGenericWithJacobianAltD failures over the
  // lifetime of this propagator. Counters are zero unless the function was
  // called (CVH ntuplizer use case).
  if (propTotalCalls_ > 0ULL) {
    const unsigned long long fail_total =
        propFailCounts_[0] + propFailCounts_[1] + propFailCounts_[2] + propFailCounts_[3] + propFailCounts_[4] +
        propFailCounts_[5];
    std::cout << "Geant4ePropagator::propagateGenericWithJacobianAltD summary"
              << "  calls="        << propTotalCalls_
              << "  failures="     << fail_total
              << " (" << (100. * fail_total / propTotalCalls_) << "% )"
              << "   exit1[plimit]=" << propFailCounts_[0]
              << "   exit2[ierr]="   << propFailCounts_[1]
              << "   exit3[maxlen]=" << propFailCounts_[2]
              << "   exit4[pdrain]=" << propFailCounts_[3]
              << "   exit5[offsurface]=" << propFailCounts_[4]
              << "   exit6[fieldbound]=" << propFailCounts_[5]
              << "   backwardLegs=" << propBackwardLegs_
              << std::endl;
  }

  // don't close the g4 Geometry here, because the propagator might have been
  // cloned
  // but there is only one, globally-shared Geometry
  delete fluct;
}

//
////////////////////////////////////////////////////////////////////////////
//

/** Propagate from a free state (e.g. position and momentum in
 *  in global cartesian coordinates) to a plane.
 */

void Geant4ePropagator::ensureGeant4eIsInitilized(bool) const {
  LogDebug("Geant4ePropagator") << "G4 propagator starts isInitialized, theField: " << theField;

  auto *theG4eManager = G4ErrorPropagatorManager::GetErrorPropagatorManager();
  auto man = G4RunManagerKernel::GetRunManagerKernel();
  if (G4StateManager::GetStateManager()->GetCurrentState() == G4State_PreInit) {
    man->SetVerboseLevel(0);
    theG4eManager->InitGeant4e();

    // define 10 mm step limit for propagator
    G4UImanager::GetUIpointer()->ApplyCommand("/geant4e/limits/stepLength 10.0 mm");
  }
  const G4Field *field = G4TransportationManager::GetTransportationManager()->GetFieldManager()->GetDetectorField();
  if (field == nullptr) {
    edm::LogError("Geant4e") << "No G4 magnetic field defined";
  }
  LogDebug("Geant4ePropagator") << "G4 propagator initialized; field: " << field;
}

void Geant4ePropagator::ensureGeant4eIsInitilizedForCVH(bool forceInit) const {
  LogDebug("Geant4e") << "ensureGeant4eIsInitilized called" << std::endl;
  if ((G4ErrorPropagatorData::GetErrorPropagatorData()->GetState() == G4ErrorState_PreInit) || forceInit) {
    LogDebug("Geant4e") << "Initializing G4 propagator" << std::endl;

    auto *theG4eManager = G4ErrorPropagatorManager::GetErrorPropagatorManager();
    theG4eManager->SetUserInitialization(new G4ErrorPhysicsListForCVH());
    theG4eManager->InitGeant4e();

    const G4Field *field = G4TransportationManager::GetTransportationManager()->GetFieldManager()->GetDetectorField();
    if (field == nullptr) {
      edm::LogError("Geant4e") << "No G4 magnetic field defined";
    }
    LogDebug("Geant4e") << "G4 propagator initialized" << std::endl;
  } else {
    LogDebug("Geant4e") << "G4 not in preinit state: " << G4ErrorPropagatorData::GetErrorPropagatorData()->GetState()
                        << std::endl;
  }
  // define 10 mm step limit for propagator
  G4UImanager::GetUIpointer()->ApplyCommand("/geant4e/limits/stepLength 10.0 mm");
}

template <>
Geant4ePropagator::ErrorTargetPair Geant4ePropagator::transformToG4SurfaceTarget(const Plane &pDest,
                                                                                 bool moveTargetToEndOfSurface) const {
  //* Get position and normal (orientation) of the destination plane
  GlobalPoint posPlane = pDest.toGlobal<double>(LocalPoint(0, 0, 0));
  GlobalVector normalPlane = pDest.toGlobal<double>(LocalVector(0, 0, 1.));
  normalPlane = normalPlane.unit();

  //* Transform this into HepGeom::Point3D<double>  and
  // HepGeom::Normal3D<double>  that define a plane for
  //  Geant4e.
  //  CMS uses cm and GeV while Geant4 uses mm and MeV
  HepGeom::Point3D<double> surfPos = TrackPropagation::globalPointToHepPoint3D(posPlane);
  HepGeom::Normal3D<double> surfNorm = TrackPropagation::globalVectorToHepNormal3D(normalPlane);

  //* Set the target surface
  return ErrorTargetPair(false, std::make_shared<G4ErrorPlaneSurfaceTarget>(surfNorm, surfPos));
}

template <>
Geant4ePropagator::ErrorTargetPair Geant4ePropagator::transformToG4SurfaceTargetD(
    const GloballyPositioned<double> &pDest, bool moveTargetToEndOfSurface) const {
  //* Get position and normal (orientation) of the destination plane
  Vector3DBase<double, GlobalTag> normalPlane = pDest.toGlobal(Vector3DBase<double, LocalTag>(0, 0, 1.));

  const Point3DBase<double, GlobalTag> &posPlane = pDest.position();

  //* Transform this into HepGeom::Point3D<double>  and
  // HepGeom::Normal3D<double>  that define a plane for
  //  Geant4e.
  //  CMS uses cm and GeV while Geant4 uses mm and MeV
  HepGeom::Point3D<double> surfPos(posPlane.x() * CLHEP::cm, posPlane.y() * CLHEP::cm, posPlane.z() * CLHEP::cm);
  HepGeom::Normal3D<double> surfNorm(normalPlane.x(), normalPlane.y(), normalPlane.z());

  //* Set the target surface
  return ErrorTargetPair(false, std::make_shared<G4ErrorPlaneSurfaceTarget>(surfNorm, surfPos));
}

template <>
Geant4ePropagator::ErrorTargetPair Geant4ePropagator::transformToG4SurfaceTarget(const Cylinder &pDest,
                                                                                 bool moveTargetToEndOfSurface) const {
  // Get Cylinder parameters.
  // CMS uses cm and GeV while Geant4 uses mm and MeV.
  // - Radius
  G4float radCyl = pDest.radius() * CLHEP::cm;
  // - Position: PositionType & GlobalPoint are Basic3DPoint<float,GlobalTag>
  G4ThreeVector posCyl = TrackPropagation::globalPointToHep3Vector(pDest.position());
  // - Rotation: Type in CMSSW is RotationType == TkRotation<T>, T=float
  G4RotationMatrix rotCyl = TrackPropagation::tkRotationFToHepRotation(pDest.rotation());

  // DEBUG
  TkRotation<float> rotation = pDest.rotation();
  LogDebug("Geant4e") << "G4e -  TkRotation" << rotation;
  LogDebug("Geant4e") << "G4e -  G4Rotation" << rotCyl << "mm";

  return ErrorTargetPair(!moveTargetToEndOfSurface, std::make_shared<G4ErrorCylSurfaceTarget>(radCyl, posCyl, rotCyl));
}

template <>
std::string Geant4ePropagator::getSurfaceType(Cylinder const &c) const {
  return "Cylinder";
}

template <>
std::string Geant4ePropagator::getSurfaceType(Plane const &c) const {
  return "Plane";
}

std::string Geant4ePropagator::generateParticleName(int charge) const {
  std::string particleName = theParticleName;

  if (charge > 0) {
    particleName += "+";
  }
  if (charge < 0) {
    particleName += "-";
  }

  LogDebug("Geant4e") << "G4e -  Particle name: " << particleName;

  return particleName;
}

template <>
bool Geant4ePropagator::configureAnyPropagation(G4ErrorMode &mode,
                                                GloballyPositioned<double> const &pDest,
                                                GlobalPoint const &cmsInitPos,
                                                GlobalVector const &cmsInitMom) const {
  if (cmsInitMom.mag() < plimit_)
    return false;
  if (pDest.toLocal(cmsInitPos).z() * pDest.toLocal(cmsInitMom).z() < 0) {
    mode = G4ErrorMode_PropForwards;
    LogDebug("Geant4e") << "G4e -  Propagator mode is \'forwards\' indirect "
                           "via the Any direction"
                        << std::endl;
  } else {
    mode = G4ErrorMode_PropBackwards;
    LogDebug("Geant4e") << "G4e -  Propagator mode is \'backwards\' indirect "
                           "via the Any direction"
                        << std::endl;
  }

  return true;
}

template <>
bool Geant4ePropagator::configureAnyPropagation(G4ErrorMode &mode,
                                                Plane const &pDest,
                                                GlobalPoint const &cmsInitPos,
                                                GlobalVector const &cmsInitMom) const {
  if (cmsInitMom.mag() < plimit_)
    return false;
  if (pDest.localZ(cmsInitPos) * pDest.localZ(cmsInitMom) < 0) {
    mode = G4ErrorMode_PropForwards;
    LogDebug("Geant4e") << "G4e -  Propagator mode is \'forwards\' indirect "
                           "via the Any direction"
                        << std::endl;
  } else {
    mode = G4ErrorMode_PropBackwards;
    LogDebug("Geant4e") << "G4e -  Propagator mode is \'backwards\' indirect "
                           "via the Any direction"
                        << std::endl;
  }

  return true;
}

template <>
bool Geant4ePropagator::configureAnyPropagation(G4ErrorMode &mode,
                                                Cylinder const &pDest,
                                                GlobalPoint const &cmsInitPos,
                                                GlobalVector const &cmsInitMom) const {
  if (cmsInitMom.mag() < plimit_)
    return false;
  //------------------------------------
  // For cylinder assume outside is backwards, inside is along
  // General use for particles from collisions
  LocalPoint lpos = pDest.toLocal(cmsInitPos);
  Surface::Side theSide = pDest.side(lpos, 0);
  if (theSide == SurfaceOrientation::positiveSide) {  // outside cylinder
    mode = G4ErrorMode_PropBackwards;
    LogDebug("Geant4e") << "G4e -  Propagator mode is \'backwards\' indirect "
                           "via the Any direction";
  } else {  // inside cylinder
    mode = G4ErrorMode_PropForwards;
    LogDebug("Geant4e") << "G4e -  Propagator mode is \'forwards\' indirect "
                           "via the Any direction";
  }

  return true;
}

template <class SurfaceType>
bool Geant4ePropagator::configurePropagation(G4ErrorMode &mode,
                                             SurfaceType const &pDest,
                                             GlobalPoint const &cmsInitPos,
                                             GlobalVector const &cmsInitMom) const {
  if (cmsInitMom.mag() < plimit_)
    return false;
  if (propagationDirection() == oppositeToMomentum) {
    mode = G4ErrorMode_PropBackwards;
    LogDebug("Geant4e") << "G4e -  Propagator mode is \'backwards\' " << std::endl;
  } else if (propagationDirection() == alongMomentum) {
    mode = G4ErrorMode_PropForwards;
    LogDebug("Geant4e") << "G4e -  Propagator mode is \'forwards\'" << std::endl;
  } else if (propagationDirection() == anyDirection) {
    if (configureAnyPropagation(mode, pDest, cmsInitPos, cmsInitMom) == false)
      return false;
  } else {
    edm::LogError("Geant4e") << "G4e - Unsupported propagation mode";
    return false;
  }
  return true;
}

template <class SurfaceType>
std::pair<TrajectoryStateOnSurface, double> Geant4ePropagator::propagateGeneric(const FreeTrajectoryState &ftsStart,
                                                                                const SurfaceType &pDest) const {
  // Deferred per-thread Geant4e init: do the first GetErrorPropagatorManager()
  // call on this thread under geant4eInitMutex() so concurrent first-time
  // inits don't corrupt G4's global geometry store. Subsequent calls on the
  // same thread skip the mutex (the bool is thread_local).
  if (!geant4eInitDoneForThread()) {
    std::lock_guard<std::mutex> lk(geant4eInitMutex());
    if (!geant4eInitDoneForThread()) {
      if (forCVH_) {
        // forceInit=true to drive G4ErrorRunManagerHelper through its
        // init -- this is what allocates G4ErrorPropagator's internal
        // navigator + transportation, needed before InitTrackPropagation
        // can run. The physics-list ConstructProcess is one-shot per
        // worker thread (thread_local flag), and this instance shares the
        // job-wide particle set registered by CvhMaster's list, so
        // re-invoking it after CvhWorker's master/worker InitializeWorker
        // neither double-registers processes nor defines extra
        // process-less particles.
        ensureGeant4eIsInitilizedForCVH(true);
      } else {
        ensureGeant4eIsInitilized(true);
      }
      geant4eInitDoneForThread() = true;
    }
  }
  // fluct is per-INSTANCE while the init flag above is per-THREAD: another
  // stream's propagator may already have initialized this thread, so the
  // allocation must not live inside the guarded block. It is deferred from
  // the ctors to here so it runs AFTER the world is installed on this thread
  // (CvhWorker) and the G4Error physics list is registered (guaranteed once
  // geant4eInitDoneForThread() is true) -- only then is
  // G4WentzelVIModel::Initialise's chained G4SafetyHelper::InitialiseHelper
  // safe to call, and only then does G4ParticleTable know about mu+ etc.
  if (forCVH_ && !fluct) {
    fluct = new G4UniversalFluctuationForExtrapolator();
    const G4ParticleDefinition *partdef =
        G4ParticleTable::GetParticleTable()->FindParticle(generateParticleName(1));
    fluct->SetParticleAndCharge(partdef, 1.);
    fluct->SetIoniTruncationAlpha(ioniTruncAlpha_);
  }
  auto *theG4eManager = G4ErrorPropagatorManager::GetErrorPropagatorManager();
  auto *theG4eData = G4ErrorPropagatorData::GetErrorPropagatorData();

  ///////////////////////////////
  // Construct the target surface
  //
  //* Set the target surface

  ErrorTargetPair g4eTarget_center = transformToG4SurfaceTarget(pDest, false);

  // * Get the starting point and direction and convert them to
  // CLHEP::Hep3Vector
  //   for G4. CMS uses cm and GeV while Geant4 uses mm and MeV
  GlobalPoint cmsInitPos = ftsStart.position();
  GlobalVector cmsInitMom = ftsStart.momentum();
  bool flipped = false;
  if (propagationDirection() == oppositeToMomentum) {
    // flip the momentum vector as Geant4 will not do this
    // on it's own in a backward propagation
    cmsInitMom = -cmsInitMom;
    flipped = true;
  }

  // Set the mode of propagation according to the propagation direction
  G4ErrorMode mode = G4ErrorMode_PropForwards;
  if (!configurePropagation(mode, pDest, cmsInitPos, cmsInitMom))
    return TsosPP(TrajectoryStateOnSurface(), 0.0f);

  // re-check propagation direction chosen in case of AnyDirection
  if (mode == G4ErrorMode_PropBackwards && !flipped)
    cmsInitMom = -cmsInitMom;

  CLHEP::Hep3Vector g4InitPos = TrackPropagation::globalPointToHep3Vector(cmsInitPos);
  CLHEP::Hep3Vector g4InitMom = TrackPropagation::globalVectorToHep3Vector(cmsInitMom * CLHEP::GeV);

  debugReportTrackState("intitial", cmsInitPos, g4InitPos, cmsInitMom, g4InitMom, pDest);

  // Set the mode of propagation according to the propagation direction
  // G4ErrorMode mode = G4ErrorMode_PropForwards;

  // if (!configurePropagation(mode, pDest, cmsInitPos, cmsInitMom))
  //	return TsosPP(TrajectoryStateOnSurface(), 0.0f);

  ///////////////////////////////
  // Set the error and trajectories, and finally propagate
  //
  G4ErrorTrajErr g4error(5, 1);
  if (ftsStart.hasError()) {
    CurvilinearTrajectoryError initErr;
    initErr = ftsStart.curvilinearError();
    g4error = TrackPropagation::algebraicSymMatrix55ToG4ErrorTrajErr(initErr, ftsStart.charge());
    LogDebug("Geant4e") << "CMS -  Error matrix: " << std::endl << initErr.matrix();
  } else {
    LogDebug("Geant4e") << "No error matrix available" << std::endl;
    return TsosPP(TrajectoryStateOnSurface(), 0.0f);
  }

  LogDebug("Geant4e") << "G4e -  Error matrix: " << std::endl << g4error;

  // in CMSSW, the state errors are deflated when performing the backward
  // propagation
  if (mode == G4ErrorMode_PropForwards) {
    G4ErrorPropagatorData::GetErrorPropagatorData()->SetStage(G4ErrorStage_Inflation);
  } else if (mode == G4ErrorMode_PropBackwards) {
    G4ErrorPropagatorData::GetErrorPropagatorData()->SetStage(G4ErrorStage_Deflation);
  }

  G4ErrorFreeTrajState g4eTrajState(generateParticleName(ftsStart.charge()), g4InitPos, g4InitMom, g4error);
  LogDebug("Geant4e") << "G4e -  Traj. State: " << (g4eTrajState);

  //////////////////////////////
  // Propagate
  int iterations = 0;
  double finalPathLength = 0;

  HepGeom::Point3D<double> finalRecoPos;

  G4ErrorPropagatorData::GetErrorPropagatorData()->SetMode(mode);

  theG4eData->SetTarget(g4eTarget_center.second.get());
  LogDebug("Geant4e") << "Running Propagation to the RECO surface" << std::endl;

  theG4eManager->InitTrackPropagation();

  if (!forCVH_) {
    // re-initialize navigator to avoid mismatches and/or segfaults
    theG4eManager->GetErrorPropagationNavigator()->LocateGlobalPointAndSetup(
        g4InitPos, &g4InitMom, /*pRelativeSearch = */ false, /*ignoreDirection = */ false);
  }

  bool continuePropagation = true;
  while (continuePropagation) {
    iterations++;
    LogDebug("Geant4e") << std::endl << "step count " << iterations << " step length " << finalPathLength;

    if (!forCVH_) {
      // re-initialize navigator to avoid mismatches and/or segfaults
      theG4eManager->GetErrorPropagationNavigator()->LocateGlobalPointWithinVolume(g4eTrajState.GetPosition());
    }

    const int ierr = theG4eManager->PropagateOneStep(&g4eTrajState, mode);

    if (ierr != 0) {
      // propagation failed, return invalid track state
      return TsosPP(TrajectoryStateOnSurface(), 0.0f);
    }

    const float thisPathLength = TrackPropagation::g4doubleToCmsDouble(g4eTrajState.GetG4Track()->GetStepLength());

    LogDebug("Geant4e") << "step Length was " << thisPathLength << " cm, current global position: "
                        << TrackPropagation::hepPoint3DToGlobalPoint(g4eTrajState.GetPosition()) << std::endl;

    finalPathLength += thisPathLength;

    // Cap deliberately loosened from the stock/10_6 value of 200 cm: low-pT
    // curling tracks (e.g. B->J/psiK bachelor kaons, helix radius ~50 cm) can
    // legitimately accumulate >2 m of path before intersecting the target
    // surface. 100 m still bounds pathological loopers.
    if (std::fabs(finalPathLength) > 10000.0f) {
      LogDebug("Geant4e") << "ERROR: Quitting propagation: path length mega large" << std::endl;
      theG4eManager->GetPropagator()->InvokePostUserTrackingAction(g4eTrajState.GetG4Track());
      continuePropagation = false;
      LogDebug("Geant4e") << "WARNING: Quitting propagation: max path length "
                             "exceeded, returning invalid state"
                          << std::endl;

      // reached maximum path length, bail out
      return TsosPP(TrajectoryStateOnSurface(), 0.0f);
    }

    if (theG4eManager->GetPropagator()->CheckIfLastStep(g4eTrajState.GetG4Track())) {
      theG4eManager->GetPropagator()->InvokePostUserTrackingAction(g4eTrajState.GetG4Track());
      continuePropagation = false;
    }
  }

  // CMSSW Tracking convention, backward propagations have negative path length
  if (propagationDirection() == oppositeToMomentum)
    finalPathLength = -finalPathLength;

  // store the correct location for the hit on the RECO surface
  LogDebug("Geant4e") << "Position on the RECO surface" << g4eTrajState.GetPosition() << std::endl;
  finalRecoPos = g4eTrajState.GetPosition();

  theG4eManager->EventTermination();

  LogDebug("Geant4e") << "Final position of the Track :" << g4eTrajState.GetPosition() << std::endl;

  //////////////////////////////
  // Retrieve the state in the end from Geant4e, convert them to CMS vectors
  // and points, and build global trajectory parameters.
  // CMS uses cm and GeV while Geant4 uses mm and MeV
  //
  const HepGeom::Vector3D<double> momEnd = g4eTrajState.GetMomentum();

  // use the hit on the the RECO plane as the final position to be d'accor with
  // the RecHit measurements
  const GlobalPoint posEndGV = TrackPropagation::hepPoint3DToGlobalPoint(finalRecoPos);
  GlobalVector momEndGV = TrackPropagation::hep3VectorToGlobalVector(momEnd) / CLHEP::GeV;

  debugReportTrackState("final", posEndGV, finalRecoPos, momEndGV, momEnd, pDest);

  // Get the error covariance matrix from Geant4e. It comes in curvilinear
  // coordinates so use the appropiate CMS class
  G4ErrorTrajErr g4errorEnd = g4eTrajState.GetError();

  CurvilinearTrajectoryError curvError(
      TrackPropagation::g4ErrorTrajErrToAlgebraicSymMatrix55(g4errorEnd, ftsStart.charge()));

  if (mode == G4ErrorMode_PropBackwards) {
    GlobalTrajectoryParameters endParm(
        posEndGV, momEndGV, ftsStart.parameters().charge(), &ftsStart.parameters().magneticField());

    // flip the momentum direction because it has been flipped before running
    // G4's backwards prop
    momEndGV = -momEndGV;
  }

  LogDebug("Geant4e") << "G4e -  Error matrix after propagation: " << std::endl << g4errorEnd;

  LogDebug("Geant4e") << "CMS -  Error matrix after propagation: " << std::endl << curvError.matrix();

  GlobalTrajectoryParameters tParsDest(posEndGV, momEndGV, ftsStart.charge(), theField);

  SurfaceSideDefinition::SurfaceSide side;

  side = propagationDirection() == alongMomentum ? SurfaceSideDefinition::afterSurface
                                                 : SurfaceSideDefinition::beforeSurface;

  return TsosPP(TrajectoryStateOnSurface(tParsDest, curvError, pDest, side), finalPathLength);
}

std::tuple<bool,
           Eigen::Matrix<double, 7, 1>,
           Eigen::Matrix<double, 5, 5>,
           Eigen::Matrix<double, 5, 9>,
           double,
           Eigen::Matrix<double, 5, 5>,
           Eigen::Matrix<double, 5, 5>,
           double,
           double>
Geant4ePropagator::propagateGenericWithJacobianAltD(const Eigen::Matrix<double, 7, 1> &ftsStart,
                                                    const GloballyPositioned<double> &pDest,
                                                    const Eigen::Vector3d &dB,
                                                    double dxi,
                                                    double dms,
                                                    double dioni,
                                                    double pforced,
                                                    const std::string &particleNameOverride,
                                                    const MaterialGroupModel *matGroups,
                                                    std::vector<std::pair<int, Eigen::Matrix<double, 5, 1>>>
                                                        *groupJacOut,
                                                    const sim::FieldModeProvider *fieldModes,
                                                    std::vector<Eigen::Matrix<double, 5, 1>> *modeJacOut) const {
  using namespace Eigen;

  // Urban + Moliere step logs: one propagate call = one leg; clearing here
  // keeps the logs in sync with the call whose noise matrices the caller
  // consumes.
  if (ioniStepLogging_) {
    ioniStepLog_.clear();
    msStepLog_.clear();
    radStepLog_.clear();
  }
  if (stepTransportLogging_) {
    stepTransportLog_.clear();
  }

  // Deferred per-thread Geant4e init under mutex (see propagateGeneric).
  if (!geant4eInitDoneForThread()) {
    std::lock_guard<std::mutex> lk(geant4eInitMutex());
    if (!geant4eInitDoneForThread()) {
      if (forCVH_) {
        // forceInit=true to drive G4ErrorRunManagerHelper through its
        // init -- this is what allocates G4ErrorPropagator's internal
        // navigator + transportation, needed before InitTrackPropagation
        // can run. The physics-list ConstructProcess is one-shot per
        // worker thread (thread_local flag), and this instance shares the
        // job-wide particle set registered by CvhMaster's list, so
        // re-invoking it after CvhWorker's master/worker InitializeWorker
        // neither double-registers processes nor defines extra
        // process-less particles.
        ensureGeant4eIsInitilizedForCVH(true);
      } else {
        ensureGeant4eIsInitilized(true);
      }
      geant4eInitDoneForThread() = true;
    }
  }
  // Per-instance fluct allocation, outside the per-thread guard.
  // See propagateGeneric for the rationale.
  if (forCVH_ && !fluct) {
    fluct = new G4UniversalFluctuationForExtrapolator();
    const G4ParticleDefinition *partdef =
        G4ParticleTable::GetParticleTable()->FindParticle(generateParticleName(1));
    fluct->SetParticleAndCharge(partdef, 1.);
    fluct->SetIoniTruncationAlpha(ioniTruncAlpha_);
  }
  auto *theG4eManager = G4ErrorPropagatorManager::GetErrorPropagatorManager();
  auto *theG4eData = G4ErrorPropagatorData::GetErrorPropagatorData();

  const G4Field *field = G4TransportationManager::GetTransportationManager()->GetFieldManager()->GetDetectorField();
  //FIXME check thread safety of this
  sim::Field *cmsField = const_cast<sim::Field *>(static_cast<const sim::Field *>(field));

  cmsField->SetOffset(dB.x(), dB.y(), dB.z());
  cmsField->SetMaterialOffset(dxi);
  // Global material model: volume-resolved k_g applied per step by
  // G4ErrorEnergyLossForCVH on top of the leg-constant dxi.
  cmsField->SetMaterialOffsetProvider(matGroups);
  if (groupJacOut != nullptr) {
    groupJacOut->clear();
  }
  if (fieldModes != nullptr && modeJacOut != nullptr) {
    modeJacOut->assign(fieldModes->nModes(), Matrix<double, 5, 1>::Zero());
  }

  auto retDefault = [cmsField]() {
    cmsField->SetOffset(0., 0., 0.);
    cmsField->SetMaterialOffset(0.);
    cmsField->SetMaterialOffsetProvider(nullptr);
    return std::tuple<bool,
                      Matrix<double, 7, 1>,
                      Matrix<double, 5, 5>,
                      Matrix<double, 5, 9>,
                      double,
                      Matrix<double, 5, 5>,
                      Matrix<double, 5, 5>,
                      double,
                      double>(false,
                              Matrix<double, 7, 1>::Zero(),
                              Matrix<double, 5, 5>::Zero(),
                              Matrix<double, 5, 9>::Zero(),
                              0.,
                              Matrix<double, 5, 5>::Zero(),
                              Matrix<double, 5, 5>::Zero(),
                              0.,
                              0.);
  };

  ///////////////////////////////
  // Construct the target surface
  //
  //* Set the target surface

  //TODO fix precision of this
  ErrorTargetPair g4eTarget_center = transformToG4SurfaceTargetD(pDest, false);

  CLHEP::Hep3Vector g4InitPos(ftsStart[0] * CLHEP::cm, ftsStart[1] * CLHEP::cm, ftsStart[2] * CLHEP::cm);
  CLHEP::Hep3Vector g4InitMom(ftsStart[3] * CLHEP::GeV, ftsStart[4] * CLHEP::GeV, ftsStart[5] * CLHEP::GeV);

  // * Get the starting point and direction and convert them to
  // CLHEP::Hep3Vector
  //   for G4. CMS uses cm and GeV while Geant4 uses mm and MeV
  GlobalPoint cmsInitPos(ftsStart[0], ftsStart[1], ftsStart[2]);
  GlobalVector cmsInitMom(ftsStart[3], ftsStart[4], ftsStart[5]);
  bool flipped = false;
  if (propagationDirection() == oppositeToMomentum) {
    // flip the momentum vector as Geant4 will not do this
    // on it's own in a backward propagation
    cmsInitMom = -cmsInitMom;
    g4InitMom = -g4InitMom;
    flipped = true;
  }

  const double charge = ftsStart[6];

  // bookkeeping: count every entry to this function
  ++propTotalCalls_;

  // Set the mode of propagation according to the propagation direction
  G4ErrorMode mode = G4ErrorMode_PropForwards;
  if (!configurePropagation(mode, pDest, cmsInitPos, cmsInitMom)) {
    ++propFailCounts_[0];
    std::cout << "Geant4e fail[plimit]"
              << "  p="        << cmsInitMom.mag()
              << "  plimit="   << plimit_
              << "  charge="   << charge
              << "  particle=" << (particleNameOverride.empty()
                                       ? generateParticleName(static_cast<int>(charge))
                                       : particleNameOverride)
              << std::endl;
    return retDefault();
  }

  // re-check propagation direction chosen in case of AnyDirection
  if (mode == G4ErrorMode_PropBackwards && !flipped) {
    cmsInitMom = -cmsInitMom;
    g4InitMom = -g4InitMom;
  }

  ///////////////////////////////
  // Set the error and trajectories, and finally propagate
  //
  G4ErrorTrajErr g4error(5, 0);

  LogDebug("Geant4e") << "G4e -  Error matrix: " << std::endl << g4error;

  // in CMSSW, the state errors are deflated when performing the backward
  // propagation
  if (mode == G4ErrorMode_PropForwards) {
    G4ErrorPropagatorData::GetErrorPropagatorData()->SetStage(G4ErrorStage_Inflation);
  } else if (mode == G4ErrorMode_PropBackwards) {
    G4ErrorPropagatorData::GetErrorPropagatorData()->SetStage(G4ErrorStage_Deflation);
  }

  // Per-call particle-name override (used by the V0 CVH ntuplizer to
  // propagate pions / protons / kaons rather than muons). Empty string =
  // use the propagator's stored particle name (default "mu+"/"mu-"). The
  // override string must already include any +/- sign (or be a baryon name
  // like "proton" / "anti_proton") since it is passed to G4 verbatim.
  const std::string g4ParticleName = particleNameOverride.empty()
                                         ? generateParticleName(charge)
                                         : particleNameOverride;
  G4ErrorFreeTrajState g4eTrajState(g4ParticleName, g4InitPos, g4InitMom, g4error);
  LogDebug("Geant4e") << "G4e -  Traj. State: " << (g4eTrajState);

  //////////////////////////////
  // Propagate
  int iterations = 0;
  double finalPathLength = 0;

  HepGeom::Point3D<double> finalRecoPos;

  G4ErrorPropagatorData::GetErrorPropagatorData()->SetMode(mode);

  theG4eData->SetTarget(g4eTarget_center.second.get());
  LogDebug("Geant4e") << "Running Propagation to the RECO surface" << std::endl;

  theG4eManager->InitTrackPropagation();

  // re-initialize navigator to avoid mismatches and/or segfaults
  theG4eManager->GetErrorPropagationNavigator()->LocateGlobalPointAndSetup(g4InitPos, &g4InitMom, false, false);

  // initial jacobian is the identity matrix for the state components,
  // and 0 for b-field and material variations
  Matrix<double, 5, 9> jac;
  jac.leftCols<5>() = Matrix<double, 5, 5>::Identity();
  jac.rightCols<4>() = Matrix<double, 5, 4>::Zero();

  const G4ErrorTrajErr errnull(5, 0);

  Matrix<double, 5, 5> g4errorEnd = Matrix<double, 5, 5>::Zero();
  Matrix<double, 5, 5> dQ = Matrix<double, 5, 5>::Zero();
  Matrix<double, 5, 5> dQ2 = Matrix<double, 5, 5>::Zero();

  double dEdxlast = 0.;

  Matrix<double, 5, 5> dErrorDxLast = Matrix<double, 5, 5>::Zero();

  double RItotal = 0.;
  double deltaTotal = 0.;
  double wTotal = 0.;

  bool continuePropagation = true;
  while (continuePropagation) {
    // re-initialize navigator to avoid mismatches and/or segfaults
    theG4eManager->GetErrorPropagationNavigator()->LocateGlobalPointWithinVolume(g4eTrajState.GetPosition());

    iterations++;
    LogDebug("Geant4e") << std::endl << "step count " << iterations << " step length " << finalPathLength;

    Matrix<double, 7, 1> statepre;
    statepre[0] = g4eTrajState.GetPosition().x() / CLHEP::cm;
    statepre[1] = g4eTrajState.GetPosition().y() / CLHEP::cm;
    statepre[2] = g4eTrajState.GetPosition().z() / CLHEP::cm;
    statepre[3] = g4eTrajState.GetMomentum().x() / CLHEP::GeV;
    statepre[4] = g4eTrajState.GetMomentum().y() / CLHEP::GeV;
    statepre[5] = g4eTrajState.GetMomentum().z() / CLHEP::GeV;
    statepre[6] = charge;

    // Per-step field-mode application: one basis sample per step at the
    // predicted step midpoint (half the 10 mm step-length cap ahead along
    // the momentum -- steps in the tracker essentially always run at the
    // cap). The SAME sample provides the applied offset for this step and,
    // below, the per-mode derivative columns, so application and
    // derivatives are consistent by construction. The constant dB argument
    // stays available on top for FD perturbations.
    Eigen::Vector3d dBstep = dB;
    const double *modeBx = nullptr;
    const double *modeBy = nullptr;
    const double *modeBz = nullptr;
    if (fieldModes != nullptr) {
      const double pmag = statepre.segment<3>(3).norm();
      const double look = 0.5;  // cm, half the 10 mm step cap
      const Eigen::Vector3d peval =
          statepre.head<3>() + (look / pmag) * statepre.segment<3>(3);
      double cstep[3];
      fieldModes->sampleAt(peval[0], peval[1], peval[2], cstep, modeBx, modeBy, modeBz);
      dBstep += Eigen::Vector3d(cstep[0], cstep[1], cstep[2]);
      cmsField->SetOffset(dBstep.x(), dBstep.y(), dBstep.z());
    }

    //set the error matrix to null to disentangle MS and ionization contributions
    g4eTrajState.SetError(errnull);
    const int ierr = theG4eManager->PropagateOneStep(&g4eTrajState, mode);

    if (ierr != 0) {
      // propagation failed, return invalid track state
      ++propFailCounts_[1];
      std::cout << "Geant4e fail[ierr=" << ierr << "]"
                << "  pT="     << cmsInitMom.perp()
                << "  eta="    << cmsInitMom.eta()
                << "  phi="    << cmsInitMom.phi()
                << "  charge=" << charge
                << "  r0="     << std::hypot(ftsStart[0], ftsStart[1])
                << "  z0="     << ftsStart[2]
                << "  surf_r=" << std::hypot(pDest.position().x(), pDest.position().y())
                << "  surf_z=" << pDest.position().z()
                << "  iter="   << iterations
                << "  pathLen="<< finalPathLength
                << "  particle=" << g4ParticleName
                << std::endl;
      return retDefault();
    }

    // Field-model validity bound: every tracker module (and therefore every
    // legitimate module-to-module leg) lies strictly inside the field
    // model's defined region (ScalarPot3D: sphere R = 320 cm vs. outermost
    // module corners at R ~ 295 cm). A state outside it is unambiguous
    // proof of a runaway/wrong-way leg, caught here within one G4 step
    // (<= 10 mm) of leaving -- metres earlier and thousands of steps
    // cheaper than the momentum-drain or off-surface checks below. For
    // field models defined everywhere (volume-based map) isDefined() is
    // always true and this guard never fires.
    {
      const GlobalPoint curPos(g4eTrajState.GetPosition().x() / CLHEP::cm,
                               g4eTrajState.GetPosition().y() / CLHEP::cm,
                               g4eTrajState.GetPosition().z() / CLHEP::cm);
      if (!theField->isDefined(curPos)) {
        ++propFailCounts_[5];
        std::cout << "Geant4e fail[fieldbound]"
                  << "  r="        << curPos.perp()
                  << "  z="        << curPos.z()
                  << "  p="        << g4eTrajState.GetMomentum().mag() / CLHEP::GeV
                  << "  pT0="      << cmsInitMom.perp()
                  << "  eta0="     << cmsInitMom.eta()
                  << "  charge="   << charge
                  << "  surf_r="   << std::hypot(pDest.position().x(), pDest.position().y())
                  << "  surf_z="   << pDest.position().z()
                  << "  iter="     << iterations
                  << "  pathLen="  << finalPathLength
                  << "  particle=" << g4ParticleName
                  << std::endl;
        theG4eManager->GetPropagator()->InvokePostUserTrackingAction(g4eTrajState.GetG4Track());
        return retDefault();
      }
    }

    // In-flight momentum floor: a leg whose momentum drains below plimit_
    // can only end at the extrapolator table floor (Ekin ~ 1 MeV) after
    // grinding through meters of dense material -- a runaway/wrong-way leg
    // that will never reach the intended (bounded) module. Without this
    // check such a leg is returned as a *success* with p ~ 15 MeV, and the
    // 1/p^n transport-Jacobian terms poison the downstream fit (NaN state
    // updates). Fail fast instead, mirroring the entry-point plimit check.
    if (g4eTrajState.GetMomentum().mag() / CLHEP::GeV < plimit_) {
      ++propFailCounts_[3];
      std::cout << "Geant4e fail[pdrain]"
                << "  p="        << g4eTrajState.GetMomentum().mag() / CLHEP::GeV
                << "  plimit="   << plimit_
                << "  pT0="      << cmsInitMom.perp()
                << "  eta0="     << cmsInitMom.eta()
                << "  charge="   << charge
                << "  r="        << g4eTrajState.GetPosition().perp() / CLHEP::cm
                << "  z="        << g4eTrajState.GetPosition().z() / CLHEP::cm
                << "  surf_r="   << std::hypot(pDest.position().x(), pDest.position().y())
                << "  surf_z="   << pDest.position().z()
                << "  iter="     << iterations
                << "  pathLen="  << finalPathLength
                << "  particle=" << g4ParticleName
                << std::endl;
      theG4eManager->GetPropagator()->InvokePostUserTrackingAction(g4eTrajState.GetG4Track());
      return retDefault();
    }

    const double thisPathLength = TrackPropagation::g4doubleToCmsDouble(g4eTrajState.GetG4Track()->GetStepLength());

    const double ePre = g4eTrajState.GetG4Track()->GetStep()->GetPreStepPoint()->GetTotalEnergy() / CLHEP::GeV;
    const double ePost = g4eTrajState.GetG4Track()->GetStep()->GetPostStepPoint()->GetTotalEnergy() / CLHEP::GeV;
    const double dEdx = (ePost - ePre) / thisPathLength;
    const double mass = g4eTrajState.GetG4Track()->GetDynamicParticle()->GetMass() / CLHEP::GeV;

    if (std::abs(thisPathLength) > 0.) {
      dEdxlast = dEdx;
    }

    const Matrix<double, 5, 9> transportJac = transportJacobianBxByBzD(statepre, thisPathLength, dEdx, mass, dBstep);

    // transport contribution to error
    g4errorEnd = (transportJac.leftCols<5>() * g4errorEnd * transportJac.leftCols<5>().transpose()).eval();

    dQ = (transportJac.leftCols<5>() * dQ * transportJac.leftCols<5>().transpose()).eval();
    dQ2 = (transportJac.leftCols<5>() * dQ2 * transportJac.leftCols<5>().transpose()).eval();

    // Global material model, step group + coherent (M1) scaling: the same
    // per-step k_g that scales the mean energy loss (applied by the eloss
    // process via the provider) coherently scales the step's MS and
    // ionization-fluctuation variances -- e^k more material means e^k more
    // scattering and fluctuation power. Only the CURRENT k values enter
    // (two-step scheme: the fit never differentiates through the weights).
    int stepGroup = -1;
    double matStepFact = 1.;
    if (matGroups != nullptr) {
      const G4Step *stpm = g4eTrajState.GetG4Track()->GetStep();
      const G4ThreeVector midm =
          0.5 * (stpm->GetPreStepPoint()->GetPosition() + stpm->GetPostStepPoint()->GetPosition());
      const G4LogicalVolume *lvm =
          stpm->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
      stepGroup = matGroups->classify(lvm, midm.perp() / CLHEP::cm, midm.z() / CLHEP::cm);
      matStepFact = std::exp(matGroups->offsetOf(stepGroup));
    }

    Matrix<double, 5, 5> errMSIout = PropagateErrorMSC(g4eTrajState.GetG4Track(), pforced);

    // scaling only affects MS
    const double msfact = std::exp(dms) * matStepFact;

    errMSIout *= msfact;

    const G4Material *mate = g4eTrajState.GetG4Track()->GetVolume()->GetLogicalVolume()->GetMaterial();
    const double X0 = mate->GetRadlen() / CLHEP::cm;
    RItotal += msfact * thisPathLength / X0;

    if (ioniStepLogging_ && thisPathLength > 0.) {
      const G4Track *trk = g4eTrajState.GetG4Track();
      const double pGeV = (pforced > 0. ? pforced : trk->GetMomentum().mag() / CLHEP::GeV);
      const double mass = trk->GetDynamicParticle()->GetMass() / CLHEP::GeV;
      const double beta = pGeV / std::sqrt(pGeV * pGeV + mass * mass);
      MoliereMsStep ms;
      CalculateEffectiveZandA(mate, ms.effZ, ms.effA);
      // per-element Moliere sums (effZ/effA are mass averages and both
      // Moliere parameters are non-linear in Z -- see MoliereMsStep)
      CalculateMoliereSums(mate, beta, ms.zzp1OverA, ms.lnScreenW);
      // areal density rho*d in g/cm^2 (GetDensity in G4 internal units)
      ms.xg = (mate->GetDensity() / (CLHEP::g / CLHEP::cm3)) * thisPathLength;
      ms.pGeV = pGeV;
      ms.beta = beta;
      // projected-angle variance as it enters Q (post-msfact):
      // curvilinear (qop, lambda, phi, yT, zT) -> (lambda, lambda)
      ms.thp2 = errMSIout(1, 1);
      ms.dOverX0 = thisPathLength / X0;
      ms.stepGroup = stepGroup;
      msStepLog_.push_back(ms);

      // radiative (brems + pair) record for the offline CF, one per step and
      // aligned with the Moliere log so an offline reader can zip them
      RadiativeStep rs;
      rs.effZ = ms.effZ;
      rs.effA = ms.effA;
      rs.xg = ms.xg;
      rs.etotGeV = std::sqrt(pGeV * pGeV + mass * mass);
      rs.pGeV = pGeV;
      rs.dOverX0 = ms.dOverX0;
      rs.stepCm = thisPathLength;
      computeRadiativeDEDX(trk, rs.dedxBrem, rs.dedxPair);
      rs.dedxRad = rs.dedxBrem + rs.dedxPair;
      rs.cs = rs.etotGeV / (pGeV * pGeV * pGeV);
      fillRadiativeSpectrum(trk, rs);
      radStepLog_.push_back(rs);
    }

    const double ionifact = std::exp(dioni) * matStepFact;

    errMSIout(0, 0) = ionifact * computeErrorIoni(g4eTrajState.GetG4Track(), pforced);

    // separate scaling for ionization and MS

    g4errorEnd += errMSIout;

    if (std::abs(thisPathLength) > 0.) {
      dErrorDxLast = errMSIout / thisPathLength;
    }

    // transport  + nominal energy loss contribution to jacobian
    jac = (transportJac.leftCols<5>() * jac).eval();

    // Per-step cumulative transport, logged AFTER the update so that it is
    // the transport from the leg start to the point where this step's noise
    // was just injected into g4errorEnd/dQ/dQ2 above (see StepTransport).
    if (stepTransportLogging_) {
      StepTransport st;
      Map<Matrix<double, 5, 5, RowMajor>>(st.jacc) = jac.leftCols<5>();
      st.nMs = static_cast<int>(msStepLog_.size());
      st.nIoni = static_cast<int>(ioniStepLog_.size());
      stepTransportLog_.push_back(st);
    }

    //TODO assess relevance of approximations (does the order matter? position/momentum before or after step?)
    //b-field (dBx, dBy, dBz) and material (dxi) contributions to jacobian
    jac.rightCols<4>() += transportJac.rightCols<4>();

    // Per-step field modes: transport the accumulated per-mode columns and
    // add this step's contribution, scaled by the SAME basis sample the
    // applied offset used above -- exact consistency between application
    // and derivative (the FD closure probes precisely this).
    if (fieldModes != nullptr && modeJacOut != nullptr) {
      for (auto &mc : *modeJacOut) {
        mc = (transportJac.leftCols<5>() * mc).eval();
      }
      const unsigned int nmodes = fieldModes->nModes();
      for (unsigned int im = 0; im < nmodes; ++im) {
        (*modeJacOut)[im] += transportJac.col(5) * modeBx[im] + transportJac.col(6) * modeBy[im] +
                             transportJac.col(7) * modeBz[im];
      }
    }

    // Global material model: per-group split of the dxi column. Same
    // recursion as the integrated column above (transport previous
    // accumulations by this step's state Jacobian, add this step's dxi
    // column to the group of the step's volume), so the sum over groups
    // equals jac.col(8) exactly.
    if (matGroups != nullptr && groupJacOut != nullptr) {
      for (auto &gc : *groupJacOut) {
        gc.second = (transportJac.leftCols<5>() * gc.second).eval();
      }
      // stepGroup classified above (M1 block)
      auto it = std::find_if(groupJacOut->begin(), groupJacOut->end(),
                             [stepGroup](auto const &e) { return e.first == stepGroup; });
      if (it == groupJacOut->end()) {
        groupJacOut->emplace_back(stepGroup, transportJac.col(8));
      } else {
        it->second += transportJac.col(8);
      }
    }

    Matrix<double, 5, 5> errMS = errMSIout;
    errMS(0, 0) = 0.;

    Matrix<double, 5, 5> errI = errMSIout;
    errI.bottomRightCorner<4, 4>() *= 0.;

    dQ += errMS;
    dQ2 += errI;

    LogDebug("Geant4e") << "step Length was " << thisPathLength << " cm, current global position: "
                        << TrackPropagation::hepPoint3DToGlobalPoint(g4eTrajState.GetPosition()) << std::endl;

    finalPathLength += thisPathLength;

    // Cap deliberately loosened from the stock/10_6 value of 200 cm: low-pT
    // curling tracks (e.g. B->J/psiK bachelor kaons, helix radius ~50 cm) can
    // legitimately accumulate >2 m of path before intersecting the target
    // surface. 100 m still bounds pathological loopers.
    if (std::fabs(finalPathLength) > 10000.0f) {
      LogDebug("Geant4e") << "ERROR: Quitting propagation: path length mega large" << std::endl;
      theG4eManager->GetPropagator()->InvokePostUserTrackingAction(g4eTrajState.GetG4Track());
      continuePropagation = false;
      LogDebug("Geant4e") << "WARNING: Quitting propagation: max path length "
                             "exceeded, returning invalid state"
                          << std::endl;

      // reached maximum path length, bail out
      ++propFailCounts_[2];
      std::cout << "Geant4e fail[maxlen]"
                << "  iter="     << iterations
                << "  pathLen="  << finalPathLength
                << "  pT="       << cmsInitMom.perp()
                << "  eta="      << cmsInitMom.eta()
                << "  phi="      << cmsInitMom.phi()
                << "  charge="   << charge
                << "  r0="       << std::hypot(ftsStart[0], ftsStart[1])
                << "  z0="       << ftsStart[2]
                << "  surf_r="   << std::hypot(pDest.position().x(), pDest.position().y())
                << "  surf_z="   << pDest.position().z()
                << "  particle=" << g4ParticleName
                << std::endl;

      return retDefault();
    }

    if (theG4eManager->GetPropagator()->CheckIfLastStep(g4eTrajState.GetG4Track())) {
      theG4eManager->GetPropagator()->InvokePostUserTrackingAction(g4eTrajState.GetG4Track());
      continuePropagation = false;
    }
  }

  // CMSSW Tracking convention, backward propagations have negative path length
  if (propagationDirection() == oppositeToMomentum)
    finalPathLength = -finalPathLength;

  // store the correct location for the hit on the RECO surface
  LogDebug("Geant4e") << "Position on the RECO surface" << g4eTrajState.GetPosition() << std::endl;
  finalRecoPos = g4eTrajState.GetPosition();

  theG4eManager->EventTermination();

  LogDebug("Geant4e") << "Final position of the Track :" << g4eTrajState.GetPosition() << std::endl;

  // On-surface closure check: G4's target-reached logic (CheckIfLastStep)
  // can report success without the state being anywhere near the target
  // plane -- observed when a leg whose target plane is behind the state
  // escapes to the G4 world boundary (r ~ 17.5 m) and every subsequent leg
  // returns a zero-step "success" pinned at the world edge. Such states
  // poison the fit with astronomic residuals. Require the final position
  // to lie on the destination plane within 1 mm.
  {
    const Point3DBase<double, GlobalTag> finalPosCms(finalRecoPos.x() / CLHEP::cm,
                                                     finalRecoPos.y() / CLHEP::cm,
                                                     finalRecoPos.z() / CLHEP::cm);
    const double distToPlane = pDest.toLocal(finalPosCms).z();
    if (std::abs(distToPlane) > 0.1) {
      ++propFailCounts_[4];
      std::cout << "Geant4e fail[offsurface]"
                << "  dPlane="   << distToPlane
                << "  pT0="      << cmsInitMom.perp()
                << "  eta0="     << cmsInitMom.eta()
                << "  charge="   << charge
                << "  r="        << finalPosCms.perp()
                << "  z="        << finalPosCms.z()
                << "  surf_r="   << std::hypot(pDest.position().x(), pDest.position().y())
                << "  surf_z="   << pDest.position().z()
                << "  iter="     << iterations
                << "  pathLen="  << finalPathLength
                << "  particle=" << g4ParticleName
                << std::endl;
      return retDefault();
    }
  }

  //////////////////////////////
  // Retrieve the state in the end from Geant4e, convert them to CMS vectors
  // and points, and build global trajectory parameters.
  // CMS uses cm and GeV while Geant4 uses mm and MeV
  //

  Matrix<double, 7, 1> ftsEnd;
  ftsEnd[0] = g4eTrajState.GetPosition().x() / CLHEP::cm;
  ftsEnd[1] = g4eTrajState.GetPosition().y() / CLHEP::cm;
  ftsEnd[2] = g4eTrajState.GetPosition().z() / CLHEP::cm;
  ftsEnd[3] = g4eTrajState.GetMomentum().x() / CLHEP::GeV;
  ftsEnd[4] = g4eTrajState.GetMomentum().y() / CLHEP::GeV;
  ftsEnd[5] = g4eTrajState.GetMomentum().z() / CLHEP::GeV;
  ftsEnd[6] = charge;

  // Backward legs (anyDirection mode picking PropBackwards): the momentum
  // was flipped before running G4's backward propagation; flip it back so
  // the returned state has physical (along-track) momentum, mirroring the
  // TSOS-based propagateGeneric.
  //
  // Exact frame conversion of the derivative/noise outputs: everything was
  // accumulated on the momentum-FLIPPED trajectory, i.e. in the curvilinear
  // frame of the reversed direction. Under momentum reversal the CMSSW
  // curvilinear basis transforms as q/p even, lambda odd, dphi even (phi
  // shifts by the constant pi), xT odd (U = Z x T flips with T), yT even
  // (V = T x U is invariant) -- the flip map has Jacobian
  // P = diag(1, -1, 1, -1, 1), which is its own inverse. Hence:
  //   J_phys(state block) = P * J_flip * P   (rows: end frame, cols: start frame)
  //   J_phys(dB/dxi cols) = P * cols         (the perturbations are global-frame)
  //   Q_phys              = P * Q_flip * P   (both indices in the end frame)
  //   dEdx_phys           = -dEdx_flip       (retraced path gains energy where
  //                                           the physical track loses it)
  // With this conversion the per-leg gradients and noise matrices are exact
  // for backward legs of any length, not only the mm-scale marginal legs.
  if (mode == G4ErrorMode_PropBackwards) {
    ftsEnd[3] = -ftsEnd[3];
    ftsEnd[4] = -ftsEnd[4];
    ftsEnd[5] = -ftsEnd[5];

    Matrix<double, 5, 5> Pflip = Matrix<double, 5, 5>::Identity();
    Pflip(1, 1) = -1.;
    Pflip(3, 3) = -1.;

    jac.leftCols<5>() = (Pflip * jac.leftCols<5>() * Pflip).eval();
    jac.rightCols<4>() = (Pflip * jac.rightCols<4>()).eval();
    // the cumulative per-step transports carry one index in the leg-start
    // frame and one in the step frame, so they flip on both sides like jac
    for (auto &st : stepTransportLog_) {
      Map<Matrix<double, 5, 5, RowMajor>> m(st.jacc);
      m = (Pflip * m * Pflip).eval();
    }
    // per-group dxi columns transform like the integrated dxi column
    if (groupJacOut != nullptr) {
      for (auto &gc : *groupJacOut) {
        gc.second = (Pflip * gc.second).eval();
      }
    }
    // per-mode field columns transform like the dB columns
    if (modeJacOut != nullptr) {
      for (auto &mc : *modeJacOut) {
        mc = (Pflip * mc).eval();
      }
    }
    g4errorEnd = (Pflip * g4errorEnd * Pflip).eval();
    dQ = (Pflip * dQ * Pflip).eval();
    dQ2 = (Pflip * dQ2 * Pflip).eval();
    dEdxlast = -dEdxlast;

    ++propBackwardLegs_;
  }

  cmsField->SetOffset(0., 0., 0.);
  cmsField->SetMaterialOffset(0.);
  cmsField->SetMaterialOffsetProvider(nullptr);

  return std::tuple<bool,
                    Matrix<double, 7, 1>,
                    Matrix<double, 5, 5>,
                    Matrix<double, 5, 9>,
                    double,
                    Matrix<double, 5, 5>,
                    Matrix<double, 5, 5>,
                    double,
                    double>(true, ftsEnd, g4errorEnd, jac, dEdxlast, dQ, dQ2, deltaTotal, wTotal);
}

//
////////////////////////////////////////////////////////////////////////////
//

/** The methods propagateWithPath() are identical to the corresponding
 *  methods propagate() in what concerns the resulting
 *  TrajectoryStateOnSurface, but they provide in addition the
 *  exact path length along the trajectory.
 */

std::pair<TrajectoryStateOnSurface, double> Geant4ePropagator::propagateWithPath(const FreeTrajectoryState &ftsStart,
                                                                                 const Plane &pDest) const {
  // Finally build the pair<...> that needs to be returned where the second
  // parameter is the exact path length. Currently calculated with a stepping
  // action that adds up the length of every step
  return propagateGeneric(ftsStart, pDest);
}

std::pair<TrajectoryStateOnSurface, double> Geant4ePropagator::propagateWithPath(const FreeTrajectoryState &ftsStart,
                                                                                 const Cylinder &cDest) const {
  // Finally build the pair<...> that needs to be returned where the second
  // parameter is the exact path length.
  return propagateGeneric(ftsStart, cDest);
}

std::pair<TrajectoryStateOnSurface, double> Geant4ePropagator::propagateWithPath(
    const TrajectoryStateOnSurface &tsosStart, const Plane &pDest) const {
  // Finally build the pair<...> that needs to be returned where the second
  // parameter is the exact path length.
  const FreeTrajectoryState ftsStart = *tsosStart.freeState();
  return propagateGeneric(ftsStart, pDest);
}

std::pair<TrajectoryStateOnSurface, double> Geant4ePropagator::propagateWithPath(
    const TrajectoryStateOnSurface &tsosStart, const Cylinder &cDest) const {
  const FreeTrajectoryState ftsStart = *tsosStart.freeState();
  // Finally build the pair<...> that needs to be returned where the second
  // parameter is the exact path length.
  return propagateGeneric(ftsStart, cDest);
}

void Geant4ePropagator::debugReportPlaneSetup(GlobalPoint const &posPlane,
                                              HepGeom::Point3D<double> const &surfPos,
                                              GlobalVector const &normalPlane,
                                              HepGeom::Normal3D<double> const &surfNorm,
                                              const Plane &pDest) const {
  LogDebug("Geant4e") << "G4e -  Destination CMS plane position:" << posPlane << "cm\n"
                      << "G4e -                  (Ro, eta, phi): (" << posPlane.perp() << " cm, " << posPlane.eta()
                      << ", " << posPlane.phi().degrees() << " deg)\n"
                      << "G4e -  Destination G4  plane position: " << surfPos << " mm, Ro = " << surfPos.perp()
                      << " mm";
  LogDebug("Geant4e") << "G4e -  Destination CMS plane normal  : " << normalPlane << "\n"
                      << "G4e -  Destination G4  plane normal  : " << normalPlane;
  LogDebug("Geant4e") << "G4e -  Distance from plane position to plane: " << pDest.localZ(posPlane) << " cm";
}

template <class SurfaceType>
void Geant4ePropagator::debugReportTrackState(std::string const &currentContext,
                                              GlobalPoint const &cmsInitPos,
                                              CLHEP::Hep3Vector const &g4InitPos,
                                              GlobalVector const &cmsInitMom,
                                              CLHEP::Hep3Vector const &g4InitMom,
                                              const SurfaceType &pDest) const {
  LogDebug("Geant4e") << "G4e - Current Context: " << currentContext;
  LogDebug("Geant4e") << "G4e -  CMS point position:" << cmsInitPos << "cm\n"
                      << "G4e -              (Ro, eta, phi): (" << cmsInitPos.perp() << " cm, " << cmsInitPos.eta()
                      << ", " << cmsInitPos.phi().degrees() << " deg)\n"
                      << "G4e -   G4  point position: " << g4InitPos << " mm, Ro = " << g4InitPos.perp() << " mm";
  LogDebug("Geant4e") << "G4e -   CMS momentum      :" << cmsInitMom << "GeV\n"
                      << " pt: " << cmsInitMom.perp() << "G4e -  G4  momentum      : " << g4InitMom << " MeV";
}

//------------------------------------------------------------------------
Eigen::Matrix<double, 5, 5> Geant4ePropagator::PropagateErrorMSC(const G4Track *aTrack, double pforced) const {
  G4ThreeVector vpPre = aTrack->GetMomentum() / CLHEP::GeV;
  G4double mass = aTrack->GetDynamicParticle()->GetMass() / CLHEP::GeV;
  G4double pPre = pforced > 0. ? pforced : aTrack->GetMomentum().mag() / CLHEP::GeV;
  G4double Etot = sqrt(pPre * pPre + mass * mass);
  G4double beta = pPre / Etot;
  G4double pBeta = pPre * beta;
  G4double stepLengthCm = aTrack->GetStep()->GetStepLength() / CLHEP::cm;

  G4Material *mate = aTrack->GetVolume()->GetLogicalVolume()->GetMaterial();
  G4double effZ, effA;
  CalculateEffectiveZandA(mate, effZ, effA);

#ifdef G4EVERBOSE
  if (iverbose >= 4)
    G4cout << "material " << mate->GetName() << " effZ:" << effZ << " effA:" << effA
           << " dens(g/mole):" << mate->GetDensity() / CLHEP::g * CLHEP::mole
           << " Radlen/cm:" << mate->GetRadlen() / CLHEP::cm << " nuclLen/cm"
           << mate->GetNuclearInterLength() / CLHEP::cm << G4endl;
#endif

  G4double RI = stepLengthCm / (mate->GetRadlen() / CLHEP::cm);
#ifdef G4EVERBOSE
  if (iverbose >= 4)
    G4cout << std::setprecision(6) << std::setw(6) << "G4EP:MSC: RI=X/X0 " << RI << " stepLengthCm " << stepLengthCm
           << " radlen/cm " << (mate->GetRadlen() / CLHEP::cm) << " RI*1.e10:" << RI * 1.e10 << G4endl;
#endif
  G4double charge = aTrack->GetDynamicParticle()->GetCharge();
  G4double X0 = mate->GetRadlen() / CLHEP::cm;
  G4double Xs = X0 * (effZ + 1.) * std::log(287. / std::sqrt(effZ)) / std::log(159. * std::pow(effZ, -1. / 3.)) / effZ;

  G4double DD = 2.25e-4 * stepLengthCm * (charge / pBeta * charge / pBeta) / Xs;

#ifdef G4EVERBOSE
  if (iverbose >= 3)
    G4cout << "G4EP:MSC: D*1E6= " << DD * 1.E6 << " pBeta " << pBeta << G4endl;
#endif
  G4double S1 = DD * stepLengthCm * stepLengthCm / 3.;
  G4double S2 = DD;
  G4double S3 = DD * stepLengthCm / 2.;

  G4double CLA = std::sqrt(vpPre.x() * vpPre.x() + vpPre.y() * vpPre.y()) / pPre;
#ifdef G4EVERBOSE
  if (iverbose >= 2)
    G4cout << std::setw(6) << "G4EP:MSC: RI " << RI << " S1 " << S1 << " S2 " << S2 << " S3 " << S3 << " CLA " << CLA
           << G4endl;
#endif
  Eigen::Matrix<double, 5, 5> res = Eigen::Matrix<double, 5, 5>::Zero();
  res(1, 1) = S2;
  res(1, 4) = -S3;
  res(2, 2) = S2 / CLA / CLA;
  res(2, 3) = S3 / CLA;
  res(3, 3) = S1;
  res(4, 4) = S1;

  res(4, 1) = res(1, 4);
  res(3, 2) = res(2, 3);

#ifdef G4EVERBOSE
  if (iverbose >= 2)
    G4cout << "G4EP:MSC: error matrix propagated msc " << fError << G4endl;
#endif

  return res;
}

//------------------------------------------------------------------------
std::pair<double, double> Geant4ePropagator::computeLandau(const G4Track *aTrack) const {
  G4double stepLengthCm = aTrack->GetStep()->GetStepLength() / CLHEP::cm;

  if (stepLengthCm <= 0.) {
    return std::make_pair(0., 0.);
  }

  //  *     Calculate xi factor (KeV).
  G4Material *mate = aTrack->GetVolume()->GetLogicalVolume()->GetMaterial();
  G4double effZ, effA;
  CalculateEffectiveZandA(mate, effZ, effA);

  G4double mass = aTrack->GetDynamicParticle()->GetMass() / CLHEP::GeV;
  G4double pPre = aTrack->GetMomentum().mag() / CLHEP::GeV;
  G4double Etot = sqrt(pPre * pPre + mass * mass);
  G4double beta = pPre / Etot;
  G4double gamma = Etot / mass;

  // *     Maximum energy transfer to atomic electron (KeV).
  G4double eta = beta * gamma;
  G4double etasq = eta * eta;
  G4double eMass = 0.51099906 / CLHEP::GeV;
  G4double massRatio = eMass / mass;
  G4double F1 = 2 * eMass * etasq;
  G4double F2 = 1. + 2. * massRatio * gamma + massRatio * massRatio;
  G4double Emax = 1.E+6 * F1 / F2;  // now in keV

  G4double Emaxmev = Emax * 1e-3;
  G4double stepLengthmm = stepLengthCm * 10;
  const double ekinmev = aTrack->GetStep()->GetPreStepPoint()->GetKineticEnergy();

  const double ePremev = aTrack->GetStep()->GetPreStepPoint()->GetTotalEnergy();
  const double ePostmev = aTrack->GetStep()->GetPostStepPoint()->GetTotalEnergy();
  const double elossmev = ePremev - ePostmev;

  const unsigned int nsamples = 1;

  std::vector<double> elossv;
  elossv.reserve(nsamples);
  for (unsigned int isample = 0; isample < nsamples; ++isample) {
    double eloss =
        1e-3 * fluct->SampleFluctuations2(mate, aTrack->GetDynamicParticle(), Emaxmev, stepLengthmm, ekinmev, elossmev);

    elossv.push_back(eloss);
  }

  std::sort(elossv.begin(), elossv.end());

  unsigned int imode = 0.24672310 * nsamples;
  unsigned int imed = nsamples / 2;

  const double mode = elossv[imode];
  const double med = elossv[imed];

  constexpr double k = -0.22278;
  constexpr double km = 1.35578;

  const double deltap = mode;
  const double c = std::max(0., 0.5 * M_PI * (med - mode) / (km - k));

  const double ascale = 2. * c / M_PI;
  const double mshift = deltap - ascale * k;

  const double mu = c > 0. ? mshift - ascale * std::log(ascale) : deltap;

  return std::make_pair(mu, c);
}

//------------------------------------------------------------------------
namespace {
  // ComputeDMicroscopicCrossSection is protected in both models; a trivial
  // derived class exposes it without touching Geant4. Using G4's own
  // differential cross sections (rather than reimplementing Petrukhin-
  // Shestakov / Kelner-Kokoulin-Petrukhin offline) removes any re-derivation
  // risk AND guarantees the spectrum is consistent with the dE/dx table and
  // with the full simulation, which run the same models.
  struct BremProbe : public G4MuBremsstrahlungModel {
    explicit BremProbe(const G4ParticleDefinition *p) : G4MuBremsstrahlungModel(p) {}
    using G4MuBremsstrahlungModel::ComputeDMicroscopicCrossSection;
  };
  struct PairProbe : public G4MuPairProductionModel {
    explicit PairProbe(const G4ParticleDefinition *p) : G4MuPairProductionModel(p) {}
    using G4MuPairProductionModel::ComputeDMicroscopicCrossSection;
  };
}  // namespace

void Geant4ePropagator::radVGrid(double *v) {
  const double lo = std::log(kRadVMin), hi = std::log(kRadVMax);
  for (int i = 0; i < kNRadV; ++i) {
    v[i] = std::exp(lo + (hi - lo) * double(i) / double(kNRadV - 1));
  }
}

void Geant4ePropagator::fillRadiativeSpectrum(const G4Track *aTrack, RadiativeStep &rs) const {
  const G4ParticleDefinition *part = aTrack->GetDynamicParticle()->GetParticleDefinition();
  if (std::abs(part->GetPDGEncoding()) != 13) {
    return;
  }
  static thread_local BremProbe *brem = nullptr;
  static thread_local PairProbe *pair = nullptr;
  static thread_local const G4ParticleDefinition *muPlus = nullptr;
  if (brem == nullptr) {
    muPlus = G4MuonPlus::MuonPlus();
    G4DataVector cuts(std::max<size_t>(G4Material::GetNumberOfMaterials(), 1), DBL_MAX);
    brem = new BremProbe(muPlus);
    pair = new PairProbe(muPlus);
    brem->Initialise(muPlus, cuts);
    pair->Initialise(muPlus, cuts);
    brem->SetUseBaseMaterials(false);
    pair->SetUseBaseMaterials(false);
  }

  const G4Material *mate = aTrack->GetVolume()->GetLogicalVolume()->GetMaterial();
  const double ePre = aTrack->GetStep()->GetPreStepPoint()->GetKineticEnergy();
  const double ePost = aTrack->GetStep()->GetPostStepPoint()->GetKineticEnergy();
  const double tkin = 0.5 * (ePre + ePost);
  const double etot = tkin + part->GetPDGMass();
  const double stepLen = aTrack->GetStep()->GetStepLength();

  double v[kNRadV];
  radVGrid(v);

  // sum over the material's ACTUAL elements with their atom densities -- the
  // cross sections are per atom and go as Z^2, so an effZ shortcut would bias
  // mixtures. n[i] is in 1/mm^3, stepLen in mm, dsigma/deps in mm^2/MeV, and
  // deps = E dv, so n * L * dsigma/deps * E is dimensionless: dN/dv.
  const G4ElementVector *elems = mate->GetElementVector();
  const double *natoms = mate->GetVecNbOfAtomsPerVolume();
  const size_t nel = mate->GetNumberOfElements();

  for (int i = 0; i < kNRadV; ++i) {
    const double eps = v[i] * etot;
    double sb = 0., sp = 0.;
    // above the kinematic limit the models are not meaningful; leave zero
    if (eps > 0. && eps < tkin) {
      for (size_t ie = 0; ie < nel; ++ie) {
        const double Z = (*elems)[ie]->GetZ();
        const double w = natoms[ie] * stepLen * etot;
        sb += w * brem->ComputeDMicroscopicCrossSection(tkin, Z, eps);
        sp += w * pair->ComputeDMicroscopicCrossSection(tkin, Z, eps);
      }
    }
    rs.dNdvBrem[i] = (sb > 0. && std::isfinite(sb)) ? sb : 0.;
    rs.dNdvPair[i] = (sp > 0. && std::isfinite(sp)) ? sp : 0.;
  }
}

void Geant4ePropagator::computeRadiativeDEDX(const G4Track *aTrack, double &dedxBrem, double &dedxPair) const {
  dedxBrem = 0.;
  dedxPair = 0.;
  const G4ParticleDefinition *part = aTrack->GetDynamicParticle()->GetParticleDefinition();
  // radiative loss ~ 1/m^2: a muon-only effect at tracker momenta
  if (std::abs(part->GetPDGEncoding()) != 13) {
    return;
  }

  // Build against muonPlus with unrestricted cuts, EXACTLY as
  // G4TablesForExtrapolatorForCVH::ComputeMuonDEDX does (that table is built
  // for muonPlus and used for both charges), so what is exported here is the
  // radiative part of the mean the propagator actually subtracts. Cached per
  // thread: model construction is expensive and this runs per Geant4 step.
  static thread_local G4MuPairProductionModel *pairModel = nullptr;
  static thread_local G4MuBremsstrahlungModel *bremModel = nullptr;
  static thread_local const G4ParticleDefinition *muPlus = nullptr;
  if (pairModel == nullptr) {
    muPlus = G4MuonPlus::MuonPlus();
    G4DataVector cuts(std::max<size_t>(G4Material::GetNumberOfMaterials(), 1), DBL_MAX);
    pairModel = new G4MuPairProductionModel(muPlus);
    bremModel = new G4MuBremsstrahlungModel(muPlus);
    pairModel->Initialise(muPlus, cuts);
    bremModel->Initialise(muPlus, cuts);
    pairModel->SetUseBaseMaterials(false);
    bremModel->SetUseBaseMaterials(false);
  }

  const G4Material *mate = aTrack->GetVolume()->GetLogicalVolume()->GetMaterial();
  const double ePre = aTrack->GetStep()->GetPreStepPoint()->GetKineticEnergy();
  const double ePost = aTrack->GetStep()->GetPostStepPoint()->GetKineticEnergy();
  const double ekin = 0.5 * (ePre + ePost);

  // ComputeDEDXPerVolume(material, particle, kineticEnergy, cut); the table
  // passes e for both energy and cut, i.e. unrestricted -- matched here.
  // Kept SEPARATE so each tabulated shape can be normalized to its own
  // process mean; the sum is exactly what the mean-loss table adds.
  const double u = CLHEP::GeV / CLHEP::cm;
  dedxBrem = bremModel->ComputeDEDXPerVolume(mate, muPlus, ekin, ekin) / u;
  dedxPair = pairModel->ComputeDEDXPerVolume(mate, muPlus, ekin, ekin) / u;
}

double Geant4ePropagator::computeErrorIoni(const G4Track *aTrack, double pforced) const {
  G4double stepLengthCm = aTrack->GetStep()->GetStepLength() / CLHEP::cm;
#ifdef G4EVERBOSE
  G4double DEDX2;
  if (stepLengthCm < 1.E-7) {
    DEDX2 = 0.;
  }
#endif
  //  *     Calculate xi factor (KeV).
  G4Material *mate = aTrack->GetVolume()->GetLogicalVolume()->GetMaterial();
  G4double effZ, effA;
  CalculateEffectiveZandA(mate, effZ, effA);

  G4double mass = aTrack->GetDynamicParticle()->GetMass() / CLHEP::GeV;
  G4double pPre = pforced > 0. ? pforced : aTrack->GetMomentum().mag() / CLHEP::GeV;
  G4double Etot = sqrt(pPre * pPre + mass * mass);
  G4double beta = pPre / Etot;
  G4double gamma = Etot / mass;

  // *     Calculate xi factor (keV).
  G4double XI = 153.5 * effZ * stepLengthCm * (mate->GetDensity() / CLHEP::mg * CLHEP::mole) / (effA * beta * beta);

#ifdef G4EVERBOSE
  if (iverbose >= 2) {
    G4cout << "G4EP:IONI: XI/keV " << XI << " beta " << beta << " gamma " << gamma << G4endl;
    G4cout << " density " << (mate->GetDensity() / CLHEP::mg * CLHEP::mole) << " effA " << effA << " step "
           << stepLengthCm << G4endl;
  }
#endif
  // *     Maximum energy transfer to atomic electron (KeV).
  G4double eta = beta * gamma;
  G4double etasq = eta * eta;
  G4double eMass = 0.51099906 / CLHEP::GeV;
  G4double massRatio = eMass / mass;
  G4double F1 = 2 * eMass * etasq;
  G4double F2 = 1. + 2. * massRatio * gamma + massRatio * massRatio;
  G4double Emax = 1.E+6 * F1 / F2;  // now in keV

  /*The above  formula for var(1/p) good for dens scatterers. However, for MIPS
    passing through a gas it leads to overestimation. Further more for incident
    electrons the Emax is almost equal to incident energy. This leads  to
    k=Xi/Emax  as small as e-6  and gradually the cov matrix  explodes.
    http://www2.pv.infn.it/~rotondi/kalman_1.pdf
    Since I do not have enough info at the moment to implement Landau &
    sub-Landau models for k=Xi/Emax <0.01 I'll saturate k at this value for now
  */

  G4double Emaxmev = Emax * 1e-3;
  G4double stepLengthmm = stepLengthCm * 10;
  const double ePre = aTrack->GetStep()->GetPreStepPoint()->GetKineticEnergy() / CLHEP::GeV;
  const double ePost = aTrack->GetStep()->GetPostStepPoint()->GetKineticEnergy() / CLHEP::GeV;
  G4double ekinmev = 0.5 * (ePre + ePost) * 1e3;

  G4double dedxsqurban =
      1e-6 * fluct->SampleFluctuations(mate, aTrack->GetDynamicParticle(), Emaxmev, stepLengthmm, ekinmev);

  G4double dedxSq = dedxsqurban;

  if (ioniStepLogging_ && fluct->lastRecordValid()) {
    const double etotGeV = aTrack->GetTotalEnergy() / CLHEP::GeV;
    const double pGeV = aTrack->GetStep()->GetPreStepPoint()->GetMomentum().mag() / CLHEP::GeV;
    UrbanIoniStep stepRec;
    stepRec.rec = fluct->lastRecord();
    stepRec.cs = etotGeV / (pGeV * pGeV * pGeV);
    ioniStepLog_.push_back(stepRec);
  }

#ifdef G4EVERBOSE
  if (iverbose >= 2)
    G4cout << "G4EP:IONI: DEDX^2(GeV^2) " << dedxSq << " emass/GeV: " << eMass << " Emax/keV: " << Emax
           << "  k=Xi/Emax=" << XI / Emax << G4endl;

#endif

  Etot = aTrack->GetTotalEnergy() / CLHEP::GeV;

  G4double pPre6 = (aTrack->GetStep()->GetPreStepPoint()->GetMomentum() / CLHEP::GeV).mag();
  pPre6 = std::pow(pPre6, 6);
  // Apply it to error
  const double res = Etot * Etot * dedxSq / pPre6;
#ifdef G4EVERBOSE
  if (iverbose >= 2)
    G4cout << "G4:IONI Etot/GeV: " << Etot << " err_dedx^2/GeV^2: " << dedxSq << " p^6: " << pPre6 << G4endl;
  if (iverbose >= 2)
    G4cout << "G4EP:IONI: error2_from_ionisation " << (Etot * Etot * dedxSq) / pPre6 << G4endl;
#endif

  return res;
}

//------------------------------------------------------------------------
void Geant4ePropagator::CalculateMoliereSums(const G4Material *mate, double beta,
                                             double &zzp1OverA, double &lnScreenW) {
  // Both Moliere parameters are non-linear in Z, and Geant4 evaluates them
  // PER ELEMENT (G4WentzelOKandVIxSection::SetupTarget is called with each
  // element's Z, and the cross sections are summed). Averaging Z first --
  // which is what effZ/effA do -- is therefore wrong for compounds.
  //
  //   chi_c,i^2      ~ w_i Z_i(Z_i+1)/A_i                       (scattering power)
  //   chi_a,i^2      ~ Z_i^(2/3) (1.13 + 3.76 (alpha Z_i/beta)^2)
  //                      * (1 + exp(-Z_i^2/1000))               (G4 screening radius,
  //                                G4WentzelOKandVIxSection.cc:154)
  //
  // The exponent depends on chi_a only logarithmically, so the correct
  // effective value is the scattering-power-weighted GEOMETRIC mean.
  constexpr double kAlpha = 1.0 / 137.035999084;
  zzp1OverA = 0.;
  lnScreenW = 0.;
  double wsum = 0.;
  const G4int nelem = mate->GetNumberOfElements();
  const G4double *fracVec = mate->GetFractionVector();   // MASS fractions
  const double b = (beta > 1e-6) ? beta : 1e-6;
  for (G4int ii = 0; ii < nelem; ++ii) {
    const double Z = mate->GetElement(ii)->GetZ();
    const double A = mate->GetElement(ii)->GetA() / CLHEP::g * CLHEP::mole;
    if (Z <= 0. || A <= 0.) continue;
    const double w = fracVec[ii] * Z * (Z + 1.) / A;     // chi_c,i^2 weight
    if (w <= 0.) continue;
    const double az = kAlpha * Z / b;
    const double scr = std::pow(Z, 2. / 3.) * (1.13 + 3.76 * az * az)
                       * (1. + std::exp(-Z * Z * 1.0e-3));
    zzp1OverA += w;
    lnScreenW += w * std::log(scr);
    wsum += w;
  }
  if (wsum > 0.) lnScreenW /= wsum;
}

void Geant4ePropagator::CalculateEffectiveZandA(const G4Material *mate, G4double &effZ, G4double &effA) {
  effZ = 0.;
  effA = 0.;
  G4int ii, nelem = mate->GetNumberOfElements();
  const G4double *fracVec = mate->GetFractionVector();
  for (ii = 0; ii < nelem; ii++) {
    effZ += mate->GetElement(ii)->GetZ() * fracVec[ii];
    effA += mate->GetElement(ii)->GetA() * fracVec[ii] / CLHEP::g * CLHEP::mole;
  }
}

Eigen::Matrix<double, 5, 9> Geant4ePropagator::transportJacobianBxByBzD(
    const Eigen::Matrix<double, 7, 1> &start, double s, double dEdx, double mass, const Eigen::Vector3d &dB) const {
  // Body (CSE + derivative + Eigen result blocks below) generated by
  // TrackPropagation/Geant4e/python/gen_transport_jacobian.py.
  // 5x9 columns: (qop0, lam0, phi0, xt0, yt0, dBx, dBy, dBz, dxi).
  // The dBx/dBy/dBz columns are scaled by kTeslaToInvGeV at the end so the
  // derivatives are per-Tesla; dxi is unscaled.

  if (s == 0.) {
    Eigen::Matrix<double, 5, 9> res;
    res.leftCols<5>() = Eigen::Matrix<double, 5, 5>::Identity();
    res.rightCols<4>() = Eigen::Matrix<double, 5, 4>::Zero();
    return res;
  }

  const GlobalPoint pos(start[0], start[1], start[2]);
  const GlobalVector &bfield = theField->inInverseGeV(pos);

  // Apply the 3D field offset (Tesla) so the integrated trajectory uses the
  // corrected field. kTeslaToInvGeV converts dB (Tesla) to the inInverseGeV
  // units in which `bfield` is already expressed.
  const double Bx = bfield.x() + MagneticField::kTeslaToInvGeV * dB.x();
  const double By = bfield.y() + MagneticField::kTeslaToInvGeV * dB.y();
  const double Bz = bfield.z() + MagneticField::kTeslaToInvGeV * dB.z();

  const double M0x = start[0];
  const double M0y = start[1];
  const double M0z = start[2];

  const Eigen::Matrix<double, 3, 1> W0 = start.segment<3>(3).normalized();
  const double W0x = W0[0];
  const double W0y = W0[1];
  const double W0z = W0[2];

  const double q = start[6];
  const double qop0 = q / start.segment<3>(3).norm();

  const double x0 = std::pow(mass, 2);
  const double x1 = std::pow(qop0, -2);
  const double x2 = std::sqrt(std::pow(q, 2) * x1 + x0);
  const double x3 = dEdx * s;
  const double x4 = x2 + x3;
  const double x5 = std::pow(-x0 + std::pow(x4, 2), -3.0 / 2.0);
  const double x6 = std::pow(W0z, 2);
  const double x7 = std::pow(W0x, 2);
  const double x8 = std::pow(W0y, 2);
  const double x9 = x7 + x8;
  const double x10 = 1.0/x9;
  const double x11 = std::pow(x10*x6 + 1, -1.0/2.0);
  const double x12 = std::pow(Bx, 2);
  const double x13 = std::pow(By, 2);
  const double x14 = std::pow(Bz, 2);
  const double x15 = x12 + x13 + x14;
  const double x16 = std::pow(x15, 5.0/2.0);
  const double x17 = std::sqrt(x15);
  const double x18 = s*x17;
  const double x19 = qop0*x18;
  const double x20 = std::sin(x19);
  const double x21 = x16*x20;
  const double x22 = x11*x21;
  const double x23 = std::sqrt(x9);
  const double x24 = 1.0/x23;
  const double x25 = W0z*x24;
  const double x26 = x22*x25;
  const double x27 = Bx*W0y;
  const double x28 = By*W0x;
  const double x29 = x24*x27 - x24*x28;
  const double x30 = std::pow(x15, 2);
  const double x31 = std::cos(x19);
  const double x32 = x31 - 1;
  const double x33 = x30*x32;
  const double x34 = x11*x33;
  const double x35 = x29*x34;
  const double x36 = std::pow(x15, 3);
  const double x37 = M0x*W0x;
  const double x38 = M0y*W0y;
  const double x39 = M0z*W0z + x37 + x38;
  const double x40 = M0z*(x24*x7 + x24*x8) - x25*x37 - x25*x38;
  const double x41 = W0z*x39 + x23*x40;
  const double x42 = x36*x41;
  const double x43 = x19 - x20;
  const double x44 = std::pow(x15, 3.0/2.0);
  const double x45 = Bx*x11;
  const double x46 = W0x*x24;
  const double x47 = x45*x46;
  const double x48 = By*x11;
  const double x49 = W0y*x24;
  const double x50 = x48*x49;
  const double x51 = Bz*x11;
  const double x52 = x25*x51;
  const double x53 = x47 + x50 + x52;
  const double x54 = x44*x53;
  const double x55 = x43*x54;
  const double x56 = Bz*x55 + qop0*x42 + x26 + x35;
  const double x57 = 1.0/x36;
  const double x58 = x1*x57;
  const double x59 = x56*x58;
  const double x60 = x31*x36;
  const double x61 = x11*x60;
  const double x62 = x25*x61;
  const double x63 = x22*x29;
  const double x64 = s*x17 - x18*x31;
  const double x65 = Bz*x54;
  const double x66 = 1.0/qop0;
  const double x67 = x57*x66;
  const double x68 = x67*(s*x62 - s*x63 + x42 + x64*x65);
  const double x69 = x33*x53;
  const double x70 = -Bz*x69 + W0z*x11*x24*x31*x36 - x63;
  const double x71 = x57*x70;
  const double x72 = x23*x39;
  const double x73 = x36*x72;
  const double x74 = W0y*x73;
  const double x75 = W0y*x22;
  const double x76 = x25*x45;
  const double x77 = x46*x51;
  const double x78 = x76 - x77;
  const double x79 = x23*x33;
  const double x80 = -M0x*x49 + M0y*W0x*x24;
  const double x81 = W0z*x40;
  const double x82 = W0x*x80 - W0y*x81;
  const double x83 = x36*x82;
  const double x84 = x23*x55;
  const double x85 = By*x84 + qop0*x74 + qop0*x83 + x75 - x78*x79;
  const double x86 = x24*x58;
  const double x87 = s*x61;
  const double x88 = x21*x78;
  const double x89 = s*x23;
  const double x90 = By*x23;
  const double x91 = x54*x64;
  const double x92 = W0y*x87 + x74 + x83 + x88*x89 + x90*x91;
  const double x93 = x24*x67;
  const double x94 = x49*x61;
  const double x95 = -By*x69 + x88 + x94;
  const double x96 = x57*x95;
  const double x97 = W0x*x73;
  const double x98 = W0x*x22;
  const double x99 = x25*x48;
  const double x100 = x49*x51;
  const double x101 = -x100 + x99;
  const double x102 = W0x*x81 + W0y*x80;
  const double x103 = x102*x36;
  const double x104 = Bx*x84 - qop0*x103 + qop0*x97 + x101*x79 + x98;
  const double x105 = x101*x21;
  const double x106 = Bx*x23;
  const double x107 = W0x*x87 - x103 - x105*x89 + x106*x91 + x97;
  const double x108 = -Bx*x69 + W0x*x11*x24*x31*x36 - x105;
  const double x109 = x108*x57;
  const double x110 = -x109*(-x104*x86 + x107*x93) - x71*(-x59 + x68) - x96*(-x85*x86 + x92*x93);
  const double x111 = q*x4*x5;
  const double x112 = dEdx*x111;
  const double x113 = W0z*x10;
  const double x114 = W0x*x45;
  const double x115 = W0y*x113;
  const double x116 = Bz*x11 - x113*x114 - x115*x48;
  const double x117 = Bz*x116;
  const double x118 = x43*x44;
  const double x119 = x117*x118 + x22 - x25*x35;
  const double x120 = std::pow(x15, -6);
  const double x121 = x120*x66;
  const double x122 = x121*x70;
  const double x123 = W0z*x46;
  const double x124 = x115*x51 + x48;
  const double x125 = x106*x116*x118 - x123*x22 + x124*x79;
  const double x126 = x121*x24;
  const double x127 = x108*x126;
  const double x128 = W0z*x49;
  const double x129 = W0x*x51;
  const double x130 = x113*x129 + x45;
  const double x131 = By*x116*x23*x43*x44 - x128*x22 - x130*x79;
  const double x132 = x126*x95;
  const double x133 = -x119*x122 - x125*x127 - x131*x132;
  const double x134 = Bx*W0x;
  const double x135 = By*W0y;
  const double x136 = x134*x24 + x135*x24;
  const double x137 = x45*x49;
  const double x138 = By*W0x*x11*x24 - x137;
  const double x139 = Bz*x138;
  const double x140 = x118*x139 + x136*x34;
  const double x141 = W0y*x51;
  const double x142 = x118*x138*x90 - x141*x33 + x98;
  const double x143 = Bx*x138*x23*x43*x44 - x129*x33 - x75;
  const double x144 = -x122*x140 - x127*x143 - x132*x142;
  const double x145 = W0y*x108*x24*x57 - x46*x96;
  const double x146 = x109*x123 + x128*x96 - x23*x71;
  const double x147 = 6/std::pow(x15, 4);
  const double x148 = Bx*x147;
  const double x149 = x56*x66;
  const double x150 = x148*x149;
  const double x151 = x20*x44;
  const double x152 = 5*x151;
  const double x153 = x30*x31;
  const double x154 = s*x153;
  const double x155 = qop0*x154;
  const double x156 = x15*x32;
  const double x157 = 4*x156;
  const double x158 = x157*x29;
  const double x159 = qop0*s;
  const double x160 = x151*x159;
  const double x161 = x160*x29;
  const double x162 = Bx*qop0;
  const double x163 = 6*x30;
  const double x164 = x162*x163;
  const double x165 = Bx*x53;
  const double x166 = 3*x17*x43;
  const double x167 = Bz*x166;
  const double x168 = 1.0/x17;
  const double x169 = s*x168;
  const double x170 = Bx*qop0*s*x168 - x162*x169*x31;
  const double x171 = x118*x77 + x152*x76 + x155*x76 + x158*x45 - x161*x45 + x164*x41 + x165*x167 + x170*x65 + x34*x49;
  const double x172 = x24*x66;
  const double x173 = x148*x172;
  const double x174 = qop0*x163;
  const double x175 = x174*x72;
  const double x176 = W0y*x152;
  const double x177 = W0z*x34;
  const double x178 = W0y*x155;
  const double x179 = x118*x48;
  const double x180 = x106*x157;
  const double x181 = x162*x78;
  const double x182 = x151*x89;
  const double x183 = x165*x166*x90;
  const double x184 = x170*x54;
  const double x185 = W0x*x179 + x164*x82 + x175*x27 + x176*x45 - x177 + x178*x45 - x180*x78 + x181*x182 + x183 + x184*x90;
  const double x186 = x12*x53;
  const double x187 = x166*x23;
  const double x188 = -x101*x162*x182 + x101*x180 - x102*x164 + x106*x184 + x114*x118 + x114*x152 + x114*x155 + x134*x175 + x186*x187 + x84;
  const double x189 = -x109*(-x104*x173 + x188*x24*x57*x66) - x71*(-x150 + x171*x57*x66) - x96*(-x173*x85 + x185*x24*x57*x66);
  const double x190 = By*x147;
  const double x191 = x149*x190;
  const double x192 = By*x174;
  const double x193 = By*x53;
  const double x194 = qop0*x169;
  const double x195 = By*qop0*s*x168 - By*x194*x31;
  const double x196 = x100*x118 + x152*x99 + x155*x99 + x158*x48 - x161*x48 + x167*x193 + x192*x41 + x195*x65 - x34*x46;
  const double x197 = x172*x190;
  const double x198 = x152*x48;
  const double x199 = x155*x48;
  const double x200 = x118*x45;
  const double x201 = x157*x90;
  const double x202 = qop0*x101;
  const double x203 = s*x151;
  const double x204 = x203*x90;
  const double x205 = x195*x54;
  const double x206 = W0x*x198 + W0x*x199 + W0y*x200 + x101*x201 - x102*x192 + x106*x205 + x175*x28 + x177 + x183 - x202*x204;
  const double x207 = qop0*x78;
  const double x208 = x13*x53;
  const double x209 = W0y*x179 + x135*x175 + x176*x48 + x178*x48 + x187*x208 + x192*x82 - x201*x78 + x204*x207 + x205*x90 + x84;
  const double x210 = -x109*(-x104*x197 + x206*x24*x57*x66) - x71*(-x191 + x196*x57*x66) - x96*(-x197*x85 + x209*x24*x57*x66);
  const double x211 = Bz*x147;
  const double x212 = x149*x211;
  const double x213 = Bz*x174;
  const double x214 = x14*x53;
  const double x215 = Bz*qop0*s*x168 - Bz*x194*x31;
  const double x216 = x118*x52 + x152*x52 + x155*x52 + x158*x51 - x161*x51 + x166*x214 + x213*x41 + x215*x65 + x55;
  const double x217 = x172*x211;
  const double x218 = Bz*W0y;
  const double x219 = Bz*x157;
  const double x220 = x219*x23;
  const double x221 = Bz*x182;
  const double x222 = x167*x53;
  const double x223 = x215*x54;
  const double x224 = W0x*x34 + W0z*x179 + x141*x152 + x141*x155 + x175*x218 + x207*x221 + x213*x82 - x220*x78 + x222*x90 + x223*x90;
  const double x225 = Bz*W0x;
  const double x226 = -W0y*x34 + W0z*x200 + x101*x220 - x102*x213 + x106*x222 + x106*x223 + x129*x152 + x129*x155 + x175*x225 - x202*x221;
  const double x227 = -x109*(-x104*x217 + x226*x24*x57*x66) - x71*(-x212 + x216*x57*x66) - x96*(-x217*x85 + x224*x24*x57*x66);
  const double x228 = std::pow(x95, 2);
  const double x229 = std::pow(x108, 2);
  const double x230 = x120*x228 + x120*x229;
  const double x231 = 1.0/x230;
  const double x232 = x120*x231;
  const double x233 = 1.0/(x232*std::pow(x70, 2) + 1);
  const double x234 = s*x20;
  const double x235 = x11*std::pow(x15, 7.0/2.0);
  const double x236 = x235*x25;
  const double x237 = std::pow(x230, -1.0/2.0);
  const double x238 = x237*x57;
  const double x239 = x234*x235;
  const double x240 = x239*x49;
  const double x241 = s*x60;
  const double x242 = x241*x78;
  const double x243 = x21*x53;
  const double x244 = By*s*x243;
  const double x245 = (1.0/2.0)*x120;
  const double x246 = x245*x95;
  const double x247 = x239*x46;
  const double x248 = x101*x241;
  const double x249 = x108*x245;
  const double x250 = x71/std::pow(x230, 3.0/2.0);
  const double x251 = qop0*x20;
  const double x252 = qop0*x61;
  const double x253 = x235*x251;
  const double x254 = x253*x49;
  const double x255 = qop0*x60;
  const double x256 = x255*x78;
  const double x257 = By*qop0*x243;
  const double x258 = x253*x46;
  const double x259 = x101*x255;
  const double x260 = x233*(x238*(Bz*qop0*x16*x20*x53 - x236*x251 - x252*x29) + x250*(-x246*(-2*x254 + 2*x256 + 2*x257) - x249*(2*Bx*qop0*x16*x20*x53 - 2*x258 - 2*x259)));
  const double x261 = x115*x61;
  const double x262 = x116*x33;
  const double x263 = By*x262;
  const double x264 = W0x*x113;
  const double x265 = x264*x61;
  const double x266 = x124*x21;
  const double x267 = Bx*x262;
  const double x268 = x46*x61;
  const double x269 = x100*x21;
  const double x270 = 2*x269;
  const double x271 = x138*x33;
  const double x272 = By*x271;
  const double x273 = Bx*x271;
  const double x274 = x237*x70;
  const double x275 = x22*x49;
  const double x276 = x152*x29;
  const double x277 = x155*x29;
  const double x278 = W0z*x21;
  const double x279 = x159*x24*x278;
  const double x280 = -Bx*Bz*qop0*s*x20*x44*x53 + x165*x219;
  const double x281 = std::pow(x15, -7);
  const double x282 = 6*x281;
  const double x283 = Bx*x282;
  const double x284 = 2*x26;
  const double x285 = 12*x153;
  const double x286 = x46*x48;
  const double x287 = x286*x33;
  const double x288 = x159*x21;
  const double x289 = x137*x288;
  const double x290 = 10*x151;
  const double x291 = Bx*x290;
  const double x292 = x154*x181;
  const double x293 = 8*x156;
  const double x294 = x165*x293;
  const double x295 = By*x294;
  const double x296 = By*x162*x234*x54;
  const double x297 = x33*x47;
  const double x298 = x288*x47;
  const double x299 = x101*x154;
  const double x300 = x162*x299;
  const double x301 = 2*x69;
  const double x302 = x22*x46;
  const double x303 = -x302;
  const double x304 = x193*x219;
  const double x305 = By*x282;
  const double x306 = x137*x33;
  const double x307 = x286*x288;
  const double x308 = By*x290;
  const double x309 = By*qop0;
  const double x310 = x299*x309;
  const double x311 = x33*x50;
  const double x312 = x288*x50;
  const double x313 = x154*x309;
  const double x314 = x313*x78;
  const double x315 = qop0*x234;
  const double x316 = x13*x315*x54;
  const double x317 = Bz*x282;
  const double x318 = x33*x99;
  const double x319 = Bz*x290;
  const double x320 = Bz*s*x153;
  const double x321 = x207*x320;
  const double x322 = By*x315*x65;
  const double x323 = x33*x76;
  const double x324 = x21*x77;
  const double x325 = x202*x320;
  const double x326 = x108*x232;
  const double x327 = x232*x95;
  const double x328 = x326*(-x254 + x256 + x257) - x327*(Bx*qop0*x16*x20*x53 - x258 - x259);
  const double x329 = x163*x31;
  const double x330 = Bx*x152;
  const double x331 = By*x157;
  const double x332 = x165*x331;
  const double x333 = x109*x231;
  const double x334 = x231*x96;
  const double x335 = By*x152;
  const double x336 = Bz*x152;
  const double x337 = Bx*W0z;
  const double x338 = -x225 + x337;
  const double x339 = Bz*W0z + x134 + x135;
  const double x340 = x33*x339;
  const double x341 = -By*x340 + W0y*x60 + x21*x338;
  const double x342 = std::pow(x341, 2);
  const double x343 = x6 + x9;
  const double x344 = 1.0/x343;
  const double x345 = x120*x344;
  const double x346 = x342*x345;
  const double x347 = By*W0z;
  const double x348 = -x218 + x347;
  const double x349 = -Bx*x340 + W0x*x31*x36 - x21*x348;
  const double x350 = std::pow(x349, 2);
  const double x351 = x345*x350;
  const double x352 = x346 + x351;
  const double x353 = std::pow(x352, -1.0/2.0);
  const double x354 = x341*x353;
  const double x355 = x1*x104;
  const double x356 = std::pow(x10*x343, -1.0/2.0);
  const double x357 = x10*x356;
  const double x358 = x120*x357;
  const double x359 = x349*x353;
  const double x360 = x1*x85;
  const double x361 = x121*x357;
  const double x362 = x359*x361;
  const double x363 = x354*x361;
  const double x364 = qop0*x23;
  const double x365 = -qop0*x17*x31 + qop0*x17;
  const double x366 = x365*x54;
  const double x367 = W0y*x252 + x364*x88 + x366*x90;
  const double x368 = W0x*x252 - x105*x364 + x106*x366;
  const double x369 = x362*x367 - x363*x368;
  const double x370 = x357*x57;
  const double x371 = x356*x57;
  const double x372 = 12*x281;
  const double x373 = Bx*x372;
  const double x374 = x357*x66;
  const double x375 = x373*x374;
  const double x376 = x104*x354;
  const double x377 = x359*x85;
  const double x378 = x28*x33;
  const double x379 = x27*x288;
  const double x380 = x154*x162;
  const double x381 = x338*x380;
  const double x382 = Bx*x339;
  const double x383 = x331*x382;
  const double x384 = By*x339;
  const double x385 = x162*x203*x384;
  const double x386 = x27*x329 + x278 + x330*x338 - x378 - x379 + x381 - x383 + x385;
  const double x387 = x353*x361;
  const double x388 = x104*x387;
  const double x389 = x134*x33;
  const double x390 = x134*x288;
  const double x391 = x348*x380;
  const double x392 = x12*x339;
  const double x393 = 6*Bx*W0x*x30*x31 + qop0*s*x12*x20*x339*x44 - x157*x392 - x330*x348 - x340 - x389 - x390 - x391;
  const double x394 = x387*x85;
  const double x395 = std::pow(x352, -3.0/2.0);
  const double x396 = x283*x344;
  const double x397 = 2*x278;
  const double x398 = By*x293*x382;
  const double x399 = x27*x285 + x291*x338 - 2*x378 - 2*x379 + 2*x381 + 2*x385 + x397 - x398;
  const double x400 = (1.0/2.0)*x345;
  const double x401 = x341*x400;
  const double x402 = 2*x340;
  const double x403 = 12*Bx*W0x*x30*x31 + 2*qop0*s*x12*x20*x339*x44 - x291*x348 - x293*x392 - 2*x389 - 2*x390 - 2*x391 - x402;
  const double x404 = x349*x400;
  const double x405 = x395*(x342*x396 + x350*x396 - x399*x401 - x403*x404);
  const double x406 = x361*x405;
  const double x407 = x104*x341;
  const double x408 = x349*x85;
  const double x409 = By*x372;
  const double x410 = x374*x409;
  const double x411 = x27*x33;
  const double x412 = x28*x288;
  const double x413 = x313*x348;
  const double x414 = Bx*By*qop0*s*x20*x339*x44 + 6*By*W0x*x30*x31 - x278 - x335*x348 - x383 - x411 - x412 - x413;
  const double x415 = x135*x33;
  const double x416 = x135*x288;
  const double x417 = x313*x338;
  const double x418 = x13*x339;
  const double x419 = x160*x418;
  const double x420 = x135*x329 - x157*x418 + x335*x338 - x340 - x415 - x416 + x417 + x419;
  const double x421 = x305*x344;
  const double x422 = 2*Bx*By*qop0*s*x20*x339*x44 + 12*By*W0x*x30*x31 - x308*x348 - x397 - x398 - 2*x411 - 2*x412 - 2*x413;
  const double x423 = x135*x285 - x293*x418 + x308*x338 - x402 - 2*x415 - 2*x416 + 2*x417 + 2*x419;
  const double x424 = x395*(x342*x421 + x350*x421 - x401*x423 - x404*x422);
  const double x425 = x361*x424;
  const double x426 = x317*x374;
  const double x427 = x27 - x28;
  const double x428 = -Bz*x340 + W0z*x31*x36 - x21*x427;
  const double x429 = std::pow(x15, -9);
  const double x430 = x344*x428*x429;
  const double x431 = x24*x430;
  const double x432 = x172*x430;
  const double x433 = x354*x432;
  const double x434 = x359*x432;
  const double x435 = x346*x353 + x351*x353;
  const double x436 = -x367*x433 - x368*x434 + x435*x57*x66*(qop0*x62 - qop0*x63 + x365*x65);
  const double x437 = x435*x67;
  const double x438 = x345*x353;
  const double x439 = x341*x438;
  const double x440 = x428*x439;
  const double x441 = x349*x438;
  const double x442 = x428*x441;
  const double x443 = std::pow(x15, -10);
  const double x444 = x172*x429;
  const double x445 = x444*x85;
  const double x446 = x344*x353;
  const double x447 = x428*x446;
  const double x448 = x445*x447;
  const double x449 = x344*(Bx*Bz*qop0*s*x20*x339*x44 + 6*Bx*W0z*x30*x31 - W0y*x21 - x219*x382 - x225*x33 - x288*x337 - x330*x427 - x380*x427);
  const double x450 = x104*x444;
  const double x451 = x359*x450;
  const double x452 = x354*x445;
  const double x453 = x447*x450;
  const double x454 = x405*x432;
  const double x455 = x104*x349;
  const double x456 = x341*x85;
  const double x457 = x373*x446;
  const double x458 = x344*(By*Bz*qop0*s*x20*x339*x44 + 6*By*W0z*x30*x31 + W0x*x16*x20 - x218*x33 - x219*x384 - x288*x347 - x313*x427 - x335*x427);
  const double x459 = x424*x432;
  const double x460 = x409*x446;
  const double x461 = 6*Bz*x172*x344*x428*x443;
  const double dqopdqop0 = std::pow(q, 3)*x4*x5/(std::pow(qop0, 3)*x2) - x110*x112;
  const double dqopdlam0 = -x112*x133;
  const double dqopdphi0 = -x112*x144;
  const double dqopdxt0 = -x112*x145;
  const double dqopdyt0 = -x112*x146;
  const double dqopdBx = -x112*x189;
  const double dqopdBy = -x112*x210;
  const double dqopdBz = -x112*x227;
  const double dqopdxi = -x111*x3;
  const double dlamdqop0 = x110*x260 + x233*(x238*(Bz*s*x16*x20*x53 - x234*x236 - x29*x87) + x250*(-x246*(-2*x240 + 2*x242 + 2*x244) - x249*(2*Bx*s*x16*x20*x53 - 2*x247 - 2*x248)));
  const double dlamdlam0 = x133*x260 + x233*(x238*(-x117*x33 + x25*x63 + x61) + x250*(-x246*(2*x130*x16*x20 - 2*x261 - 2*x263) - x249*(-2*x265 - 2*x266 - 2*x267)));
  const double dlamdphi0 = x144*x260 + x233*(x238*(-x136*x22 - x139*x33) + x250*(-x246*(2*x268 + x270 - 2*x272) - x249*(2*Bz*W0x*x11*x16*x20*x24 - 2*x273 - 2*x94)));
  const double dlamdxt0 = x145*x260;
  const double dlamdyt0 = x146*x260;
  const double dlamdBx = x189*x260 + x233*(-x148*x274 + x238*(6*Bx*W0z*x11*x24*x30*x31 - x275 - x276*x45 - x277*x45 - x279*x45 - x280 - x33*x77) + x250*(x228*x283 + x229*x283 - x246*(x137*x285 + x284 - 2*x287 - 2*x289 + x291*x78 + 2*x292 - x295 + 2*x296) - x249*(12*Bx*W0x*x11*x24*x30*x31 + 2*qop0*s*x12*x20*x44*x53 - x101*x291 - x186*x293 - 2*x297 - 2*x298 - 2*x300 - x301)));
  const double dlamdBy = x210*x260 + x233*(-x190*x274 + x238*(By*Bz*qop0*s*x20*x44*x53 + 6*By*W0z*x11*x24*x30*x31 - x100*x33 - x198*x29 - x199*x29 - x279*x48 - x303 - x304) + x250*(x228*x305 + x229*x305 - x246*(-x208*x293 + x285*x50 - x301 + x308*x78 - 2*x311 - 2*x312 + 2*x314 + 2*x316) - x249*(2*Bx*By*qop0*s*x20*x44*x53 + 12*By*W0x*x11*x24*x30*x31 - x101*x308 - x284 - x295 - 2*x306 - 2*x307 - 2*x310)));
  const double dlamdBz = x227*x260 + x233*(-x211*x274 + x238*(6*Bz*W0z*x11*x24*x30*x31 + qop0*s*x14*x20*x44*x53 - x157*x214 - x276*x51 - x277*x51 - x279*x51 - x33*x52 - x69) + x250*(x228*x317 + x229*x317 - x246*(-Bz*x193*x293 + x100*x285 - x159*x270 - 2*x302 - 2*x318 + x319*x78 + 2*x321 + 2*x322) - x249*(2*Bx*Bz*qop0*s*x20*x44*x53 + 12*Bz*W0x*x11*x24*x30*x31 - Bz*x294 + 2*W0y*x11*x16*x20*x24 - x101*x319 - 2*x159*x324 - 2*x323 - 2*x325)));
  const double dlamdxi = 0;
  const double dphidqop0 = x110*x328 + x326*(-x240 + x242 + x244) - x327*(Bx*s*x16*x20*x53 - x247 - x248);
  const double dphidlam0 = x133*x328 + x326*(x130*x16*x20 - x261 - x263) - x327*(-x265 - x266 - x267);
  const double dphidphi0 = x144*x328 + x326*(x268 + x269 - x272) - x327*(Bz*W0x*x11*x16*x20*x24 - x273 - x94);
  const double dphidxt0 = x145*x328;
  const double dphidyt0 = x146*x328;
  const double dphidBx = x189*x328 + x333*(-x148*x95 + x57*(x137*x329 + x26 - x287 - x289 + x292 + x296 + x330*x78 - x332)) - x334*(-x108*x148 + x57*(6*Bx*W0x*x11*x24*x30*x31 + qop0*s*x12*x20*x44*x53 - x101*x330 - x157*x186 - x297 - x298 - x300 - x69));
  const double dphidBy = x210*x328 + x333*(-x190*x95 + x57*(-x157*x208 - x311 - x312 + x314 + x316 + x329*x50 + x335*x78 - x69)) - x334*(-x108*x190 + x57*(Bx*By*qop0*s*x20*x44*x53 + 6*By*W0x*x11*x24*x30*x31 - x101*x335 - x26 - x306 - x307 - x310 - x332));
  const double dphidBz = x227*x328 + x333*(-x211*x95 + x57*(x100*x329 - x159*x269 + x303 - x304 - x318 + x321 + x322 + x336*x78)) - x334*(-x108*x211 + x57*(6*Bz*W0x*x11*x24*x30*x31 - x101*x336 - x159*x324 + x275 - x280 - x323 - x325));
  const double dphidxi = 0;
  const double dxtdqop0 = -x107*x363 + x110*x369 + x354*x355*x358 - x358*x359*x360 + x362*x92;
  const double dxtdlam0 = -x125*x363 + x131*x362 + x133*x369;
  const double dxtdphi0 = x142*x362 - x143*x363 + x144*x369;
  const double dxtdxt0 = W0x*x359*x370 + W0y*x354*x370 + x145*x369;
  const double dxtdyt0 = -x115*x359*x371 + x146*x369 + x264*x354*x371;
  const double dxtdBx = x185*x362 - x188*x363 + x189*x369 + x375*x376 - x375*x377 - x386*x388 + x393*x394 - x406*x407 + x406*x408;
  const double dxtdBy = -x206*x363 + x209*x362 + x210*x369 + x376*x410 - x377*x410 - x388*x420 + x394*x414 - x407*x425 + x408*x425;
  const double dxtdBz = x224*x362 - x226*x363 + x227*x369 + x376*x426 - x377*x426;
  const double dxtdxi = 0;
  const double dytdqop0 = -x107*x434 + x110*x436 + x354*x360*x431 + x355*x359*x431 - x433*x92 - x435*x59 + x435*x68;
  const double dytdlam0 = x119*x437 - x125*x434 - x131*x433 + x133*x436;
  const double dytdphi0 = x140*x437 - x142*x433 - x143*x434 + x144*x436;
  const double dytdxt0 = x145*x436 - x440*x46 + x442*x49;
  const double dytdyt0 = x123*x442 + x128*x440 + x146*x436 + x23*x435;
  const double dytdBx = 18*Bx*x104*x24*x344*x349*x353*x428*x443*x66 + 18*Bx*x24*x341*x344*x353*x428*x443*x66*x85 - x150*x435 + x171*x435*x57*x66 - x185*x433 - x188*x434 + x189*x436 - x386*x448 - x393*x453 - x449*x451 - x449*x452 - x454*x455 - x454*x456 + x56*x57*x66*(-x342*x457 + x346*x405 - x350*x457 + x351*x405 + x399*x439 + x403*x441);
  const double dytdBy = 18*By*x104*x24*x344*x349*x353*x428*x443*x66 + 18*By*x24*x341*x344*x353*x428*x443*x66*x85 - x191*x435 + x196*x435*x57*x66 - x206*x434 - x209*x433 + x210*x436 - x414*x453 - x420*x448 - x451*x458 - x452*x458 - x455*x459 - x456*x459 + x56*x57*x66*(-x342*x460 + x346*x424 - x350*x460 + x351*x424 + x422*x441 + x423*x439);
  const double dytdBz = x104*x359*x461 - x212*x435 + x216*x435*x67 - x224*x433 - x226*x434 + x227*x436 + x354*x461*x85;
  const double dytdxi = 0;
  Eigen::Matrix<double, 5, 9> res;
  res(0,0) = dqopdqop0;
  res(0,1) = dqopdlam0;
  res(0,2) = dqopdphi0;
  res(0,3) = dqopdxt0;
  res(0,4) = dqopdyt0;
  res(0,5) = dqopdBx;
  res(0,6) = dqopdBy;
  res(0,7) = dqopdBz;
  res(0,8) = dqopdxi;
  res(1,0) = dlamdqop0;
  res(1,1) = dlamdlam0;
  res(1,2) = dlamdphi0;
  res(1,3) = dlamdxt0;
  res(1,4) = dlamdyt0;
  res(1,5) = dlamdBx;
  res(1,6) = dlamdBy;
  res(1,7) = dlamdBz;
  res(1,8) = dlamdxi;
  res(2,0) = dphidqop0;
  res(2,1) = dphidlam0;
  res(2,2) = dphidphi0;
  res(2,3) = dphidxt0;
  res(2,4) = dphidyt0;
  res(2,5) = dphidBx;
  res(2,6) = dphidBy;
  res(2,7) = dphidBz;
  res(2,8) = dphidxi;
  res(3,0) = dxtdqop0;
  res(3,1) = dxtdlam0;
  res(3,2) = dxtdphi0;
  res(3,3) = dxtdxt0;
  res(3,4) = dxtdyt0;
  res(3,5) = dxtdBx;
  res(3,6) = dxtdBy;
  res(3,7) = dxtdBz;
  res(3,8) = dxtdxi;
  res(4,0) = dytdqop0;
  res(4,1) = dytdlam0;
  res(4,2) = dytdphi0;
  res(4,3) = dytdxt0;
  res(4,4) = dytdyt0;
  res(4,5) = dytdBx;
  res(4,6) = dytdBy;
  res(4,7) = dytdBz;
  res(4,8) = dytdxi;

  res.middleCols(5, 3) *= MagneticField::kTeslaToInvGeV;

  return res;

}
