#include <sstream>

// Geant4e
#include "TrackPropagation/Geant4e/interface/ConvertFromToCLHEP.h"
#include "TrackPropagation/Geant4e/interface/Geant4ePropagator.h"
#include "TrackPropagation/Geant4e/interface/G4ErrorPhysicsListForCVH.h"

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

/** Constructor.
 */
Geant4ePropagator::Geant4ePropagator(
    const MagneticField *field, std::string particleName, PropagationDirection dir, double plimit, bool forCVH)
    : Propagator(dir),
      theField(field),
      theParticleName(particleName),
      theG4eManager(G4ErrorPropagatorManager::GetErrorPropagatorManager()),
      theG4eData(G4ErrorPropagatorData::GetErrorPropagatorData()),
      plimit_(plimit),
      forCVH_(forCVH) {
  LogDebug("Geant4e") << "Geant4e Propagator initialized";

  // has to be called here, doing it later will not load the G4 physics list
  // properly when using the G4 ES Producer. Reason: unclear
  if (forCVH_) {
    ensureGeant4eIsInitilizedForCVH(true);
    const G4ParticleDefinition *partdef = G4ParticleTable::GetParticleTable()->FindParticle(generateParticleName(1));

    fluct = new G4UniversalFluctuationForExtrapolator();
    fluct->SetParticleAndCharge(partdef, 1.);
  } else
    ensureGeant4eIsInitilized(true);
}

/** Destructor.
 */
Geant4ePropagator::~Geant4ePropagator() {
  LogDebug("Geant4e") << "Geant4ePropagator::~Geant4ePropagator()" << std::endl;

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

    // if (std::fabs(finalPathLength) > 10000.0f)
    if (std::fabs(finalPathLength) > 200.0f) {
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
           Eigen::Matrix<double, 5, 7>,
           double,
           Eigen::Matrix<double, 5, 5>,
           Eigen::Matrix<double, 5, 5>,
           double,
           double>
Geant4ePropagator::propagateGenericWithJacobianAltD(const Eigen::Matrix<double, 7, 1> &ftsStart,
                                                    const GloballyPositioned<double> &pDest,
                                                    double dBz,
                                                    double dxi,
                                                    double dms,
                                                    double dioni,
                                                    double pforced) const {
  using namespace Eigen;

  const G4Field *field = G4TransportationManager::GetTransportationManager()->GetFieldManager()->GetDetectorField();
  //FIXME check thread safety of this
  sim::Field *cmsField = const_cast<sim::Field *>(static_cast<const sim::Field *>(field));

  cmsField->SetOffset(0., 0., dBz);
  cmsField->SetMaterialOffset(dxi);

  auto retDefault = [cmsField]() {
    cmsField->SetOffset(0., 0., 0.);
    cmsField->SetMaterialOffset(0.);
    return std::tuple<bool,
                      Matrix<double, 7, 1>,
                      Matrix<double, 5, 5>,
                      Matrix<double, 5, 7>,
                      double,
                      Matrix<double, 5, 5>,
                      Matrix<double, 5, 5>,
                      double,
                      double>(false,
                              Matrix<double, 7, 1>::Zero(),
                              Matrix<double, 5, 5>::Zero(),
                              Matrix<double, 5, 7>::Zero(),
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

  // Set the mode of propagation according to the propagation direction
  G4ErrorMode mode = G4ErrorMode_PropForwards;
  if (!configurePropagation(mode, pDest, cmsInitPos, cmsInitMom))
    return retDefault();

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

  G4ErrorFreeTrajState g4eTrajState(generateParticleName(charge), g4InitPos, g4InitMom, g4error);
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

  Matrix<double, 5, 7> jac;
  jac.leftCols<5>() = Matrix<double, 5, 5>::Identity();
  jac.rightCols<2>() = Matrix<double, 5, 2>::Zero();

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

    //set the error matrix to null to disentangle MS and ionization contributions
    g4eTrajState.SetError(errnull);
    const int ierr = theG4eManager->PropagateOneStep(&g4eTrajState, mode);

    if (ierr != 0) {
      // propagation failed, return invalid track state
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

    const Matrix<double, 5, 7> transportJac = transportJacobianBzD(statepre, thisPathLength, dEdx, mass, dBz);

    // transport contribution to error
    g4errorEnd = (transportJac.leftCols<5>() * g4errorEnd * transportJac.leftCols<5>().transpose()).eval();

    dQ = (transportJac.leftCols<5>() * dQ * transportJac.leftCols<5>().transpose()).eval();
    dQ2 = (transportJac.leftCols<5>() * dQ2 * transportJac.leftCols<5>().transpose()).eval();

    Matrix<double, 5, 5> errMSIout = PropagateErrorMSC(g4eTrajState.GetG4Track(), pforced);

    // scaling only affects MS
    const double msfact = std::exp(dms);

    errMSIout *= msfact;

    const G4Material *mate = g4eTrajState.GetG4Track()->GetVolume()->GetLogicalVolume()->GetMaterial();
    const double X0 = mate->GetRadlen() / CLHEP::cm;
    RItotal += msfact * thisPathLength / X0;

    const double ionifact = std::exp(dioni);

    errMSIout(0, 0) = ionifact * computeErrorIoni(g4eTrajState.GetG4Track(), pforced);

    // separate scaling for ionization and MS

    g4errorEnd += errMSIout;

    if (std::abs(thisPathLength) > 0.) {
      dErrorDxLast = errMSIout / thisPathLength;
    }

    // transport  + nominal energy loss contribution to jacobian
    jac = (transportJac.leftCols<5>() * jac).eval();

    //TODO assess relevance of approximations (does the order matter? position/momentum before or after step?)
    //b-field and material contributions to jacobian
    jac.rightCols<2>() += transportJac.rightCols<2>();

    Matrix<double, 5, 5> errMS = errMSIout;
    errMS(0, 0) = 0.;

    Matrix<double, 5, 5> errI = errMSIout;
    errI.bottomRightCorner<4, 4>() *= 0.;

    dQ += errMS;
    dQ2 += errI;

    LogDebug("Geant4e") << "step Length was " << thisPathLength << " cm, current global position: "
                        << TrackPropagation::hepPoint3DToGlobalPoint(g4eTrajState.GetPosition()) << std::endl;

    finalPathLength += thisPathLength;

    if (std::fabs(finalPathLength) > 200.0f) {
      LogDebug("Geant4e") << "ERROR: Quitting propagation: path length mega large" << std::endl;
      theG4eManager->GetPropagator()->InvokePostUserTrackingAction(g4eTrajState.GetG4Track());
      continuePropagation = false;
      LogDebug("Geant4e") << "WARNING: Quitting propagation: max path length "
                             "exceeded, returning invalid state"
                          << std::endl;

      // reached maximum path length, bail out
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

  cmsField->SetOffset(0., 0., 0.);
  cmsField->SetMaterialOffset(0.);

  return std::tuple<bool,
                    Matrix<double, 7, 1>,
                    Matrix<double, 5, 5>,
                    Matrix<double, 5, 7>,
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

Eigen::Matrix<double, 5, 7> Geant4ePropagator::transportJacobianBzD(
    const Eigen::Matrix<double, 7, 1> &start, double s, double dEdx, double mass, double dBz) const {
  if (s == 0.) {
    Eigen::Matrix<double, 5, 7> res;
    res.leftCols<5>() = Eigen::Matrix<double, 5, 5>::Identity();
    res.rightCols<2>() = Eigen::Matrix<double, 5, 2>::Zero();
    return res;
  }

  const GlobalPoint pos(start[0], start[1], start[2]);
  const GlobalVector &bfield = theField->inInverseGeV(pos);

  const double Bx = bfield.x();
  const double By = bfield.y();
  const double Bz = bfield.z() + 2.99792458e-3 * dBz;

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
  const double x10 = 1.0 / x9;
  const double x11 = std::pow(x10 * x6 + 1, -1.0 / 2.0);
  const double x12 = std::pow(Bz, 2);
  const double x13 = std::pow(Bx, 2) + std::pow(By, 2) + x12;
  const double x14 = std::pow(x13, 5.0 / 2.0);
  const double x15 = std::sqrt(x13);
  const double x16 = s * x15;
  const double x17 = qop0 * x16;
  const double x18 = std::sin(x17);
  const double x19 = x14 * x18;
  const double x20 = x11 * x19;
  const double x21 = std::sqrt(x9);
  const double x22 = 1.0 / x21;
  const double x23 = W0z * x22;
  const double x24 = Bx * W0y;
  const double x25 = By * W0x;
  const double x26 = x22 * x24 - x22 * x25;
  const double x27 = std::pow(x13, 2);
  const double x28 = std::cos(x17);
  const double x29 = x28 - 1;
  const double x30 = x27 * x29;
  const double x31 = x11 * x30;
  const double x32 = x26 * x31;
  const double x33 = std::pow(x13, 3);
  const double x34 = M0x * W0x;
  const double x35 = M0y * W0y;
  const double x36 = M0z * W0z + x34 + x35;
  const double x37 = M0z * (x22 * x7 + x22 * x8) - x23 * x34 - x23 * x35;
  const double x38 = W0z * x36 + x21 * x37;
  const double x39 = x33 * x38;
  const double x40 = x17 - x18;
  const double x41 = std::pow(x13, 3.0 / 2.0);
  const double x42 = Bx * x11;
  const double x43 = W0x * x22;
  const double x44 = By * x11;
  const double x45 = W0y * x22;
  const double x46 = Bz * x11;
  const double x47 = x23 * x46;
  const double x48 = x42 * x43 + x44 * x45 + x47;
  const double x49 = x41 * x48;
  const double x50 = x40 * x49;
  const double x51 = Bz * x50 + qop0 * x39 + x20 * x23 + x32;
  const double x52 = 1.0 / x33;
  const double x53 = x1 * x52;
  const double x54 = x51 * x53;
  const double x55 = x28 * x33;
  const double x56 = x11 * x55;
  const double x57 = x23 * x56;
  const double x58 = x20 * x26;
  const double x59 = s * x15 - x16 * x28;
  const double x60 = Bz * x49;
  const double x61 = 1.0 / qop0;
  const double x62 = x52 * x61;
  const double x63 = x62 * (s * x57 - s * x58 + x39 + x59 * x60);
  const double x64 = x30 * x48;
  const double x65 = -Bz * x64 + W0z * x11 * x22 * x28 * x33 - x58;
  const double x66 = x52 * x65;
  const double x67 = x21 * x36;
  const double x68 = x33 * x67;
  const double x69 = W0y * x68;
  const double x70 = W0y * x20;
  const double x71 = x23 * x42;
  const double x72 = x43 * x46;
  const double x73 = x71 - x72;
  const double x74 = x21 * x30;
  const double x75 = -M0x * x45 + M0y * W0x * x22;
  const double x76 = W0z * x37;
  const double x77 = W0x * x75 - W0y * x76;
  const double x78 = x33 * x77;
  const double x79 = By * x21;
  const double x80 = qop0 * x69 + qop0 * x78 + x50 * x79 + x70 - x73 * x74;
  const double x81 = x22 * x53;
  const double x82 = s * x56;
  const double x83 = x19 * x73;
  const double x84 = s * x21;
  const double x85 = x49 * x59;
  const double x86 = W0y * x82 + x69 + x78 + x79 * x85 + x83 * x84;
  const double x87 = x22 * x62;
  const double x88 = x45 * x56;
  const double x89 = -By * x64 + x83 + x88;
  const double x90 = x52 * x89;
  const double x91 = W0x * x68;
  const double x92 = W0x * x20;
  const double x93 = x23 * x44;
  const double x94 = x45 * x46;
  const double x95 = x93 - x94;
  const double x96 = W0x * x76 + W0y * x75;
  const double x97 = x33 * x96;
  const double x98 = Bx * x21;
  const double x99 = qop0 * x91 - qop0 * x97 + x50 * x98 + x74 * x95 + x92;
  const double x100 = x19 * x95;
  const double x101 = W0x * x82 - x100 * x84 + x85 * x98 + x91 - x97;
  const double x102 = -Bx * x64 + W0x * x11 * x22 * x28 * x33 - x100;
  const double x103 = x102 * x52;
  const double x104 = -x103 * (x101 * x87 - x81 * x99) - x66 * (-x54 + x63) - x90 * (-x80 * x81 + x86 * x87);
  const double x105 = q * x4 * x5;
  const double x106 = dEdx * x105;
  const double x107 = W0z * x10;
  const double x108 = W0x * x107;
  const double x109 = W0y * x107;
  const double x110 = Bz * x11 - x108 * x42 - x109 * x44;
  const double x111 = Bz * x110;
  const double x112 = x40 * x41;
  const double x113 = x111 * x112 + x20 - x23 * x32;
  const double x114 = std::pow(x13, -6);
  const double x115 = x114 * x61;
  const double x116 = x115 * x65;
  const double x117 = W0z * x43;
  const double x118 = x109 * x46 + x44;
  const double x119 = x110 * x112 * x98 - x117 * x20 + x118 * x74;
  const double x120 = x115 * x22;
  const double x121 = x102 * x120;
  const double x122 = W0z * x45;
  const double x123 = x108 * x46 + x42;
  const double x124 = By * x110 * x21 * x40 * x41 - x122 * x20 - x123 * x74;
  const double x125 = x120 * x89;
  const double x126 = -x113 * x116 - x119 * x121 - x124 * x125;
  const double x127 = Bx * W0x;
  const double x128 = By * W0y;
  const double x129 = x127 * x22 + x128 * x22;
  const double x130 = By * W0x * x11 * x22 - x42 * x45;
  const double x131 = Bz * x130;
  const double x132 = x112 * x131 + x129 * x31;
  const double x133 = x30 * x46;
  const double x134 = -W0y * x133 + x112 * x130 * x79 + x92;
  const double x135 = Bx * x130 * x21 * x40 * x41 - W0x * x133 - x70;
  const double x136 = -x116 * x132 - x121 * x135 - x125 * x134;
  const double x137 = W0y * x102 * x22 * x52 - x43 * x90;
  const double x138 = x103 * x117 + x122 * x90 - x21 * x66;
  const double x139 = 6 * Bz;
  const double x140 = x139 / std::pow(x13, 4);
  const double x141 = x140 * x61;
  const double x142 = x141 * x51;
  const double x143 = x18 * x41;
  const double x144 = 5 * x143;
  const double x145 = x27 * x28;
  const double x146 = qop0 * s;
  const double x147 = x145 * x146;
  const double x148 = x13 * x29;
  const double x149 = 4 * x148;
  const double x150 = x26 * x46;
  const double x151 = x146 * x150;
  const double x152 = Bz * qop0;
  const double x153 = 6 * x27;
  const double x154 = x152 * x153;
  const double x155 = x12 * x48;
  const double x156 = 3 * x15 * x40;
  const double x157 = 1.0 / x15;
  const double x158 = s * x152;
  const double x159 = Bz * qop0 * s * x157 - x157 * x158 * x28;
  const double x160 =
      x112 * x47 - x143 * x151 + x144 * x47 + x147 * x47 + x149 * x150 + x154 * x38 + x155 * x156 + x159 * x60 + x50;
  const double x161 = x141 * x22;
  const double x162 = Bz * W0y;
  const double x163 = qop0 * x153 * x67;
  const double x164 = x144 * x46;
  const double x165 = x147 * x46;
  const double x166 = W0z * x112;
  const double x167 = Bz * x73;
  const double x168 = x149 * x21;
  const double x169 = x158 * x73;
  const double x170 = x143 * x21;
  const double x171 = Bz * x48;
  const double x172 = x156 * x171;
  const double x173 = x159 * x49;
  const double x174 = W0x * x31 + W0y * x164 + W0y * x165 + x154 * x77 + x162 * x163 + x166 * x44 - x167 * x168 +
                      x169 * x170 + x172 * x79 + x173 * x79;
  const double x175 = Bz * W0x;
  const double x176 = Bz * x95;
  const double x177 = x158 * x95;
  const double x178 = W0x * x164 + W0x * x165 - W0y * x31 - x154 * x96 + x163 * x175 + x166 * x42 + x168 * x176 -
                      x170 * x177 + x172 * x98 + x173 * x98;
  const double x179 = -x103 * (-x161 * x99 + x178 * x22 * x52 * x61) - x66 * (-x142 + x160 * x52 * x61) -
                      x90 * (-x161 * x80 + x174 * x22 * x52 * x61);
  const double x180 = std::pow(x89, 2);
  const double x181 = std::pow(x102, 2);
  const double x182 = x114 * x180 + x114 * x181;
  const double x183 = 1.0 / x182;
  const double x184 = x114 * x183;
  const double x185 = 1.0 / (x184 * std::pow(x65, 2) + 1);
  const double x186 = s * x18;
  const double x187 = x11 * std::pow(x13, 7.0 / 2.0);
  const double x188 = x187 * x23;
  const double x189 = std::pow(x182, -1.0 / 2.0);
  const double x190 = x189 * x52;
  const double x191 = x186 * x187;
  const double x192 = x191 * x45;
  const double x193 = s * x55;
  const double x194 = x193 * x73;
  const double x195 = x19 * x48;
  const double x196 = By * s * x195;
  const double x197 = (1.0 / 2.0) * x114;
  const double x198 = x197 * x89;
  const double x199 = x191 * x43;
  const double x200 = x193 * x95;
  const double x201 = x102 * x197;
  const double x202 = x66 / std::pow(x182, 3.0 / 2.0);
  const double x203 = qop0 * x18;
  const double x204 = qop0 * x56;
  const double x205 = x187 * x203;
  const double x206 = x205 * x45;
  const double x207 = qop0 * x55;
  const double x208 = x207 * x73;
  const double x209 = By * qop0 * x195;
  const double x210 = x205 * x43;
  const double x211 = x207 * x95;
  const double x212 = x185 * (x190 * (Bz * qop0 * x14 * x18 * x48 - x188 * x203 - x204 * x26) +
                              x202 * (-x198 * (-2 * x206 + 2 * x208 + 2 * x209) -
                                      x201 * (2 * Bx * qop0 * x14 * x18 * x48 - 2 * x210 - 2 * x211)));
  const double x213 = x109 * x56;
  const double x214 = x110 * x30;
  const double x215 = By * x214;
  const double x216 = x108 * x56;
  const double x217 = x118 * x19;
  const double x218 = Bx * x214;
  const double x219 = x43 * x56;
  const double x220 = x19 * x94;
  const double x221 = 2 * x220;
  const double x222 = x130 * x30;
  const double x223 = By * x222;
  const double x224 = Bx * x222;
  const double x225 = x139 / std::pow(x13, 7);
  const double x226 = x20 * x43;
  const double x227 = x30 * x93;
  const double x228 = 10 * x143;
  const double x229 = x145 * x169;
  const double x230 = By * x171;
  const double x231 = 8 * x148;
  const double x232 = By * qop0 * x186 * x60;
  const double x233 = x30 * x71;
  const double x234 = x19 * x72;
  const double x235 = x145 * x177;
  const double x236 = Bx * x171;
  const double x237 = x102 * x184;
  const double x238 = x184 * x89;
  const double x239 = x237 * (-x206 + x208 + x209) - x238 * (Bx * qop0 * x14 * x18 * x48 - x210 - x211);
  const double x240 = x30 * (Bz * W0z + x127 + x128);
  const double x241 = -By * x240 + W0y * x55 + x19 * (Bx * W0z - x175);
  const double x242 = x6 + x9;
  const double x243 = 1.0 / x242;
  const double x244 = x114 * x243;
  const double x245 = std::pow(x241, 2) * x244;
  const double x246 = -Bx * x240 + W0x * x28 * x33 - x19 * (By * W0z - x162);
  const double x247 = x244 * std::pow(x246, 2);
  const double x248 = std::pow(x245 + x247, -1.0 / 2.0);
  const double x249 = x241 * x248;
  const double x250 = x1 * x99;
  const double x251 = std::pow(x10 * x242, -1.0 / 2.0);
  const double x252 = x10 * x251;
  const double x253 = x114 * x252;
  const double x254 = x246 * x248;
  const double x255 = x1 * x80;
  const double x256 = x115 * x252;
  const double x257 = x254 * x256;
  const double x258 = x249 * x256;
  const double x259 = qop0 * x21;
  const double x260 = -qop0 * x15 * x28 + qop0 * x15;
  const double x261 = x260 * x49;
  const double x262 = W0y * x204 + x259 * x83 + x261 * x79;
  const double x263 = W0x * x204 - x100 * x259 + x261 * x98;
  const double x264 = x257 * x262 - x258 * x263;
  const double x265 = x252 * x52;
  const double x266 = x251 * x52;
  const double x267 = x225 * x252;
  const double x268 = x249 * x61;
  const double x269 = x254 * x61;
  const double x270 = -Bz * x240 + W0z * x28 * x33 - x19 * (x24 - x25);
  const double x271 = x22 * x243 * x270 / std::pow(x13, 9);
  const double x272 = x249 * x271;
  const double x273 = x272 * x61;
  const double x274 = x269 * x271;
  const double x275 = x245 * x248 + x247 * x248;
  const double x276 = -x262 * x273 - x263 * x274 + x275 * x52 * x61 * (qop0 * x57 - qop0 * x58 + x260 * x60);
  const double x277 = x275 * x62;
  const double x278 = x244 * x270;
  const double x279 = x249 * x278;
  const double x280 = x254 * x278;
  const double x281 = x139 * x22 * x243 * x270 / std::pow(x13, 10);
  const double dqopdqop0 = std::pow(q, 3) * x4 * x5 / (std::pow(qop0, 3) * x2) - x104 * x106;
  const double dqopdlam0 = -x106 * x126;
  const double dqopdphi0 = -x106 * x136;
  const double dqopdxt0 = -x106 * x137;
  const double dqopdyt0 = -x106 * x138;
  const double dqopdBz = -x106 * x179;
  const double dqopdxi = -x105 * x3;
  const double dlamdqop0 = x104 * x212 + x185 * (x190 * (Bz * s * x14 * x18 * x48 - x186 * x188 - x26 * x82) +
                                                 x202 * (-x198 * (-2 * x192 + 2 * x194 + 2 * x196) -
                                                         x201 * (2 * Bx * s * x14 * x18 * x48 - 2 * x199 - 2 * x200)));
  const double dlamdlam0 =
      x126 * x212 +
      x185 * (x190 * (-x111 * x30 + x23 * x58 + x56) +
              x202 * (-x198 * (2 * x123 * x14 * x18 - 2 * x213 - 2 * x215) - x201 * (-2 * x216 - 2 * x217 - 2 * x218)));
  const double dlamdphi0 =
      x136 * x212 + x185 * (x190 * (-x129 * x20 - x131 * x30) +
                            x202 * (-x198 * (2 * x219 + x221 - 2 * x223) -
                                    x201 * (2 * Bz * W0x * x11 * x14 * x18 * x22 - 2 * x224 - 2 * x88)));
  const double dlamdxt0 = x137 * x212;
  const double dlamdyt0 = x138 * x212;
  const double dlamdBz =
      x179 * x212 + x185 * (-x140 * x189 * x65 +
                            x190 * (6 * Bz * W0z * x11 * x22 * x27 * x28 + qop0 * s * x12 * x18 * x41 * x48 -
                                    x144 * x150 - x145 * x151 - x146 * x19 * x47 - x149 * x155 - x30 * x47 - x64) +
                            x202 * (x180 * x225 + x181 * x225 -
                                    x198 * (12 * x145 * x94 - x146 * x221 + x167 * x228 - 2 * x226 - 2 * x227 +
                                            2 * x229 - x230 * x231 + 2 * x232) -
                                    x201 * (2 * Bx * Bz * qop0 * s * x18 * x41 * x48 +
                                            12 * Bz * W0x * x11 * x22 * x27 * x28 + 2 * W0y * x11 * x14 * x18 * x22 -
                                            2 * x146 * x234 - x176 * x228 - x231 * x236 - 2 * x233 - 2 * x235)));
  const double dlamdxi = 0;
  const double dphidqop0 = x104 * x239 + x237 * (-x192 + x194 + x196) - x238 * (Bx * s * x14 * x18 * x48 - x199 - x200);
  const double dphidlam0 = x126 * x239 + x237 * (x123 * x14 * x18 - x213 - x215) - x238 * (-x216 - x217 - x218);
  const double dphidphi0 =
      x136 * x239 + x237 * (x219 + x220 - x223) - x238 * (Bz * W0x * x11 * x14 * x18 * x22 - x224 - x88);
  const double dphidxt0 = x137 * x239;
  const double dphidyt0 = x138 * x239;
  const double dphidBz =
      x103 * x183 *
          (-x140 * x89 +
           x52 * (x144 * x167 - x146 * x220 - x149 * x230 + x153 * x28 * x94 - x226 - x227 + x229 + x232)) +
      x179 * x239 -
      x183 * x90 *
          (-x102 * x140 + x52 * (Bx * Bz * qop0 * s * x18 * x41 * x48 + 6 * Bz * W0x * x11 * x22 * x27 * x28 +
                                 W0y * x11 * x14 * x18 * x22 - x144 * x176 - x146 * x234 - x149 * x236 - x233 - x235));
  const double dphidxi = 0;
  const double dxtdqop0 = -x101 * x258 + x104 * x264 + x249 * x250 * x253 - x253 * x254 * x255 + x257 * x86;
  const double dxtdlam0 = -x119 * x258 + x124 * x257 + x126 * x264;
  const double dxtdphi0 = x134 * x257 - x135 * x258 + x136 * x264;
  const double dxtdxt0 = W0x * x254 * x265 + W0y * x249 * x265 + x137 * x264;
  const double dxtdyt0 = x108 * x249 * x266 - x109 * x254 * x266 + x138 * x264;
  const double dxtdBz = x174 * x257 - x178 * x258 + x179 * x264 + x267 * x268 * x99 - x267 * x269 * x80;
  const double dxtdxi = 0;
  const double dytdqop0 =
      -x101 * x274 + x104 * x276 + x250 * x254 * x271 + x255 * x272 - x273 * x86 - x275 * x54 + x275 * x63;
  const double dytdlam0 = x113 * x277 - x119 * x274 - x124 * x273 + x126 * x276;
  const double dytdphi0 = x132 * x277 - x134 * x273 - x135 * x274 + x136 * x276;
  const double dytdxt0 = x137 * x276 - x279 * x43 + x280 * x45;
  const double dytdyt0 = x117 * x280 + x122 * x279 + x138 * x276 + x21 * x275;
  const double dytdBz = -x142 * x275 + x160 * x275 * x62 - x174 * x273 - x178 * x274 + x179 * x276 + x268 * x281 * x80 +
                        x269 * x281 * x99;
  const double dytdxi = 0;
  Eigen::Matrix<double, 5, 7> res;
  res(0, 0) = dqopdqop0;
  res(0, 1) = dqopdlam0;
  res(0, 2) = dqopdphi0;
  res(0, 3) = dqopdxt0;
  res(0, 4) = dqopdyt0;
  res(0, 5) = dqopdBz;
  res(0, 6) = dqopdxi;
  res(1, 0) = dlamdqop0;
  res(1, 1) = dlamdlam0;
  res(1, 2) = dlamdphi0;
  res(1, 3) = dlamdxt0;
  res(1, 4) = dlamdyt0;
  res(1, 5) = dlamdBz;
  res(1, 6) = dlamdxi;
  res(2, 0) = dphidqop0;
  res(2, 1) = dphidlam0;
  res(2, 2) = dphidphi0;
  res(2, 3) = dphidxt0;
  res(2, 4) = dphidyt0;
  res(2, 5) = dphidBz;
  res(2, 6) = dphidxi;
  res(3, 0) = dxtdqop0;
  res(3, 1) = dxtdlam0;
  res(3, 2) = dxtdphi0;
  res(3, 3) = dxtdxt0;
  res(3, 4) = dxtdyt0;
  res(3, 5) = dxtdBz;
  res(3, 6) = dxtdxi;
  res(4, 0) = dytdqop0;
  res(4, 1) = dytdlam0;
  res(4, 2) = dytdphi0;
  res(4, 3) = dytdxt0;
  res(4, 4) = dytdyt0;
  res(4, 5) = dytdBz;
  res(4, 6) = dytdxi;

  res.col(5) *= 2.99792458e-3;

  return res;
}
