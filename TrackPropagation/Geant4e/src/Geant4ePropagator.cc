#include <sstream>

// Geant4e
#include "TrackPropagation/Geant4e/interface/ConvertFromToCLHEP.h"
#include "TrackPropagation/Geant4e/interface/Geant4ePropagator.h"
#include "TrackPropagation/Geant4e/interface/G4ErrorPhysicsListCustom.h"

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
#include "G4ErrorPropagationNavigator.hh"
#include "G4RunManagerKernel.hh"
#include "G4StateManager.hh"
#include "G4PathFinder.hh"

// CLHEP
#include <CLHEP/Units/SystemOfUnits.h>

#include "TrackingTools/AnalyticalJacobians/interface/JacobianCurvilinearToLocal.h"
#include "TrackingTools/TrajectoryParametrization/interface/CurvilinearTrajectoryParameters.h"

#include <Eigen/Geometry>

#include "SimG4Core/MagneticField/interface/Field.h"


/** Constructor.
 */
Geant4ePropagator::Geant4ePropagator(const MagneticField *field,
                                     std::string particleName,
                                     PropagationDirection dir,
                                     double plimit)
    : Propagator(dir),
      theField(field),
      theParticleName(particleName),
      theG4eManager(G4ErrorPropagatorManager::GetErrorPropagatorManager()),
      theG4eData(G4ErrorPropagatorData::GetErrorPropagatorData()),
      plimit_(plimit) {
  LogDebug("Geant4e") << "Geant4e Propagator initialized";

  // has to be called here, doing it later will not load the G4 physics list
  // properly when using the G4 ES Producer. Reason: unclear
  ensureGeant4eIsInitilized(true);
  fluct = new G4UniversalFluctuationForExtrapolator();
  fluct->SetParticleAndCharge(G4ParticleTable::GetParticleTable()->FindParticle(generateParticleName(1)), 1);
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
    //G4UImanager::GetUIpointer()->ApplyCommand("/exerror/setField -10. kilogauss");

//     G4ProductionCutsTable::GetProductionCutsTable()->UpdateCoupleTable(G4RunManagerKernel::GetRunManagerKernel()->GetCurrentWorld());

    theG4eManager->SetUserInitialization(new G4ErrorPhysicsListCustom());
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
Geant4ePropagator::ErrorTargetPair Geant4ePropagator::transformToG4SurfaceTargetD(const GloballyPositioned<double> &pDest,
                                                                                 bool moveTargetToEndOfSurface) const {
  //* Get position and normal (orientation) of the destination plane
//   Point3DBase<double, GlobalTag> posPlane = pDest.toGlobal<double>(Point3DBase<double, LocalTag>(0, 0, 0));
  Vector3DBase<double, GlobalTag> normalPlane = pDest.toGlobal(Vector3DBase<double, LocalTag>(0, 0, 1.));
                                                                                   
                                                                                   
  Point3DBase<double, GlobalTag> posPlane = pDest.position();
//   TkRotation<double> rot = pDest.rotation();
  
//   Vector3DBase<double, LocalTag> lz(0, 0, 1.);
//   Vector3DBase<double, GlobalTag> normalPlane = Vector3DBase<double, GlobalTag>(rot.multiplyInverse(lz.basicVector()));

  //* Transform this into HepGeom::Point3D<double>  and
  // HepGeom::Normal3D<double>  that define a plane for
  //  Geant4e.
  //  CMS uses cm and GeV while Geant4 uses mm and MeV
  HepGeom::Point3D<double> surfPos(posPlane.x()*cm, posPlane.y()*cm, posPlane.z()*cm);
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

  const G4Field *field = G4TransportationManager::GetTransportationManager()->GetFieldManager()->GetDetectorField();
//   //FIXME check thread safety of this
//   sim::Field *cmsField = const_cast<sim::Field*>(static_cast<const sim::Field*>(field));
//   
// //   cmsField->SetOffset(0., 0., dBz);
// //   cmsField->SetMaterialOffset(dxi);
//   cmsField->SetMaterialOffset(std::log(1e-10));
                                                                                  
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

  // re-initialize navigator to avoid mismatches and/or segfaults
  theG4eManager->GetErrorPropagationNavigator()->LocateGlobalPointAndSetup(
      g4InitPos, &g4InitMom, /*pRelativeSearch = */ false, /*ignoreDirection = */ false);

  bool continuePropagation = true;
  while (continuePropagation) {
    iterations++;
    LogDebug("Geant4e") << std::endl << "step count " << iterations << " step length " << finalPathLength;

    // re-initialize navigator to avoid mismatches and/or segfaults
    theG4eManager->GetErrorPropagationNavigator()->LocateGlobalPointWithinVolume(g4eTrajState.GetPosition());

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

 std::tuple<bool, Eigen::Matrix<double, 7, 1>, Eigen::Matrix<double, 5, 5>, Eigen::Matrix<double, 5, 7>, double, Eigen::Matrix<double, 5, 5>> Geant4ePropagator::propagateGenericWithJacobianAltD(const Eigen::Matrix<double, 7, 1> &ftsStart,
                                                                                const GloballyPositioned<double> &pDest, double dBz, double dxi, double pforced) const {
                          
  using namespace Eigen;

            
//   std::cout << "start propgation, ftsStart:\n" << ftsStart << std::endl;
  
  const G4Field *field = G4TransportationManager::GetTransportationManager()->GetFieldManager()->GetDetectorField();
  //FIXME check thread safety of this
  sim::Field *cmsField = const_cast<sim::Field*>(static_cast<const sim::Field*>(field));
  
  cmsField->SetOffset(0., 0., dBz);
  cmsField->SetMaterialOffset(dxi);
//   cmsField->SetMaterialOffset(std::log(1e-6));
//   cmsField->SetMaterialOffset(std::log(1e-10));
  
  auto retDefault = [cmsField]() {
    cmsField->SetOffset(0., 0., 0.);
    cmsField->SetMaterialOffset(0.);
//     return std::tuple<TrajectoryStateOnSurface, AlgebraicMatrix57, AlgebraicMatrix55, double>();
    return std::tuple<bool, Matrix<double, 7, 1>, Matrix<double, 5, 5>, Matrix<double, 5, 7>, double, Matrix<double, 5, 5>>(false, Matrix<double, 7, 1>::Zero(), Matrix<double, 5, 5>::Zero(), Matrix<double, 5, 7>::Zero(), 0., Matrix<double, 5, 5>::Zero());
  };
  


  
//   std::cout << "start propagation" << std::endl;
  
  ///////////////////////////////
  // Construct the target surface
  //
  //* Set the target surface

  //TODO fix precision of this
  ErrorTargetPair g4eTarget_center = transformToG4SurfaceTargetD(pDest, false);

  CLHEP::Hep3Vector g4InitPos(ftsStart[0]*cm, ftsStart[1]*cm, ftsStart[2]*cm);
  CLHEP::Hep3Vector g4InitMom(ftsStart[3]*GeV, ftsStart[4]*GeV, ftsStart[5]*GeV);
  
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

  
  
//   double posarr[4] = { ftsStart[0]*cm, ftsStart[1]*cm, ftsStart[2]*cm, 0. };
//   double fieldval[10];
//   field->GetFieldValue(posarr, fieldval);
//   
//   std::cout << "cmsInitPos:" << std::endl;
//   std::cout << cmsInitPos << std::endl;
//   std::cout << "fieldval " << fieldval[0] << " " << fieldval[1] << " " << fieldval[2] << std::endl;
  
//   CLHEP::Hep3Vector g4InitPos = TrackPropagation::globalPointToHep3Vector(cmsInitPos);
//   CLHEP::Hep3Vector g4InitMom = TrackPropagation::globalVectorToHep3Vector(cmsInitMom * GeV);
  


//   debugReportTrackState("intitial", cmsInitPos, g4InitPos, cmsInitMom, g4InitMom, pDest);

  // Set the mode of propagation according to the propagation direction
  // G4ErrorMode mode = G4ErrorMode_PropForwards;

  // if (!configurePropagation(mode, pDest, cmsInitPos, cmsInitMom))
  //	return TsosPP(TrajectoryStateOnSurface(), 0.0f);

  ///////////////////////////////
  // Set the error and trajectories, and finally propagate
  //
  G4ErrorTrajErr g4error(5, 0);
//   if (ftsStart.hasError()) {
//     CurvilinearTrajectoryError initErr;
//     initErr = ftsStart.curvilinearError();
//     g4error = TrackPropagation::algebraicSymMatrix55ToG4ErrorTrajErr(initErr, ftsStart.charge());
//     LogDebug("Geant4e") << "CMS -  Error matrix: " << std::endl << initErr.matrix();
//   } else {
//     LogDebug("Geant4e") << "No error matrix available" << std::endl;
//     return retDefault();
//   }

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
//   G4TransportationManager::GetTransportationManager()->GetPropagatorInField()->GetNavigatorForPropagating()->LocateGlobalPointAndSetup(g4InitPos, &g4InitMom);
//   theG4eManager->GetErrorPropagationNavigator()->LocateGlobalPointAndSetup(g4InitPos, &g4InitMom);
  theG4eManager->GetErrorPropagationNavigator()->LocateGlobalPointAndSetup(g4InitPos, &g4InitMom, false, false);
//   std::cout << "pathfinder = " << G4PathFinder::GetInstance() << std::endl;
//   G4PathFinder::GetInstance()->Locate(g4InitPos, g4InitMom, false);
  
//   std::cout << "pathfinder nav0 = " << G4PathFinder::GetInstance()->GetNavigator(0) << " prop navigator = " << theG4eManager->GetErrorPropagationNavigator() << std::endl;
  
  // initial jacobian is the identity matrix for the state components,
  // and 0 for b-field and material variations
//   G4ErrorMatrix jac(5, 7);
//   for (unsigned int i = 0; i < 5; ++i) {
//     for (unsigned int j = 0; j < 7; ++j) {
//       jac(i+1, j+1) = i == j ? 1. : 0.;
//     }
//   }
  
  Matrix<double, 5, 7> jac;
  jac.leftCols<5>() = Matrix<double, 5, 5>::Identity();
  jac.rightCols<2>() = Matrix<double, 5, 2>::Zero();
  
  const G4ErrorTrajErr errnull(5, 0);
  
  
//   const G4ErrorTrajErr g4errorEnd = g4eTrajState.GetError();
//   const AlgebraicMatrix55 errstart = ftsStart.curvilinearError().matrix();
  Matrix<double, 5, 5> g4errorEnd = Matrix<double, 5, 5>::Zero();
  
//   Matrix<double, 7, 1> ftsEnd = ftsStart;

  
//   G4ErrorTrajErr dQ(5, 0);
  Matrix<double, 5, 5> dQ = Matrix<double, 5, 5>::Zero();

  double dEdxlast = 0.;
//   double masslast = 0.;
  
  Matrix<double, 5, 5> dErrorDxLast = Matrix<double, 5, 5>::Zero();
  
  bool firstStep = true;

  bool continuePropagation = true;
  while (continuePropagation) {
//     if (iterations > 0) {
//       G4PathFinder::GetInstance()->ReLocate(g4eTrajState.GetPosition());
//     }
    
    // re-initialize navigator to avoid mismatches and/or segfaults
    theG4eManager->GetErrorPropagationNavigator()->LocateGlobalPointWithinVolume(g4eTrajState.GetPosition());
//     G4PathFinder::GetInstance()->ReLocate(g4eTrajState.GetPosition());
//     const G4ThreeVector momstep = g4eTrajState.GetMomentum();
//     theG4eManager->GetErrorPropagationNavigator()->LocateGlobalPointAndSetup(g4eTrajState.GetPosition(), &momstep, false, false);
//     theG4eManager->GetErrorPropagationNavigator()
//     G4PathFinder::GetInstance()->Locate(g4InitPos, g4InitMom, false);

    // brute force re-initialization at every step to be safe
//     const G4ThreeVector momstep = g4eTrajState.GetMomentum();
//     theG4eManager->GetErrorPropagationNavigator()->LocateGlobalPointAndSetup(g4eTrajState.GetPosition(), &momstep, false, false);
    
    
    iterations++;
    LogDebug("Geant4e") << std::endl << "step count " << iterations << " step length " << finalPathLength;

//     const GlobalPoint pos = TrackPropagation::hepPoint3DToGlobalPoint(g4eTrajState.GetPosition());
//     const GlobalVector mom = TrackPropagation::hep3VectorToGlobalVector(g4eTrajState.GetMomentum() / GeV);
//     constexpr double mass = 0.1056583745;
    
//     const FreeTrajectoryState statepre(pos, mom, ftsStart.charge(), &ftsStart.parameters().magneticField());
    
    Matrix<double, 7, 1> statepre;
    statepre[0] = g4eTrajState.GetPosition().x()/cm;
    statepre[1] = g4eTrajState.GetPosition().y()/cm;
    statepre[2] = g4eTrajState.GetPosition().z()/cm;
    statepre[3] = g4eTrajState.GetMomentum().x()/GeV;
    statepre[4] = g4eTrajState.GetMomentum().y()/GeV;
    statepre[5] = g4eTrajState.GetMomentum().z()/GeV;
    statepre[6] = charge;
    
//     std::cout << "prestep: " << g4eTrajState.GetG4Track()->GetStep()->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName() << std::endl;
//     std::cout << "prestep: " << g4eTrajState.GetG4Track() << std::endl;
    
    //set the error matrix to null to disentangle MS and ionization contributions
    g4eTrajState.SetError(errnull);
    const int ierr = theG4eManager->PropagateOneStep(&g4eTrajState, mode);

//     std::cout << "poststep: " << g4eTrajState.GetG4Track()->GetStep()->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << std::endl;
    
    if (ierr != 0) {
      // propagation failed, return invalid track state
      return retDefault();
    }

    const double thisPathLength = TrackPropagation::g4doubleToCmsDouble(g4eTrajState.GetG4Track()->GetStepLength());

//     std::cout << "thisPathLength = " << thisPathLength << std::endl;
    
    const double pPre = g4eTrajState.GetMomentum().mag()/GeV;
    const double ePre = g4eTrajState.GetG4Track()->GetStep()->GetPreStepPoint()->GetTotalEnergy() / GeV;
    const double ePost = g4eTrajState.GetG4Track()->GetStep()->GetPostStepPoint()->GetTotalEnergy() / GeV;
    const double dEdx = (ePost - ePre)/thisPathLength;
    const double mass  = g4eTrajState.GetG4Track()->GetDynamicParticle()->GetMass() / GeV;

    if (std::abs(thisPathLength) > 0.) {
      dEdxlast = dEdx;
//       masslast = mass;
    }
    

    Matrix<double, 7, 1> statepost;
    statepost[0] = g4eTrajState.GetPosition().x()/cm;
    statepost[1] = g4eTrajState.GetPosition().y()/cm;
    statepost[2] = g4eTrajState.GetPosition().z()/cm;
    statepost[3] = g4eTrajState.GetMomentum().x()/GeV;
    statepost[4] = g4eTrajState.GetMomentum().y()/GeV;
    statepost[5] = g4eTrajState.GetMomentum().z()/GeV;
    statepost[6] = charge;

//     const Matrix<double, 7, 1> dstate = statepost - statepre;
//     std::cout << "dstate:\n" << dstate << std::endl;

//     std::cout << "thisPathLength = " << thisPathLength << std::endl;
//     const Matrix<double, 5, 7> transportJac = transportJacobian(statepre, thisPathLength, dEdx, mass);

//     const Matrix<double, 5, 7> transportJac = transportJacobianBzD(statepre, thisPathLength, dEdx, mass, dBz);
    const Matrix<double, 5, 7> transportJac = transportJacobianBzMidpointD(statepre, statepost, thisPathLength, dEdx, mass, dBz);

//     Matrix<double, 5, 7> transportJac = transportJacobianBzD(statepre, thisPathLength, dEdx, mass, dBz);

//     transportJac(0, 1) = 0.;
//     transportJac(0, 2) = 0.;
//     transportJac(0, 3) = 0.;
//     transportJac(0, 4) = 0.;
//     transportJac(0, 5) = 0.;
    
//     const Matrix<double, 5, 7> transportJacAlt = transportJacobianD(statepre, thisPathLength, dEdx, mass);
    
//     const double dh = 1e-4;
    
//     const Matrix<double, 6, 1> dest = transportResultD(statepre, thisPathLength, dEdx, mass, 0.);
//     const Matrix<double, 6, 1> destalt = transportResultD(statepre, thisPathLength, dEdx, mass, dh);
    
//     const Matrix<double, 6, 1> ddest = (destalt-dest)/dh;
//     
//     Matrix<double, 7, 1> statepost;
//     statepost[0] = g4eTrajState.GetPosition().x()/cm;
//     statepost[1] = g4eTrajState.GetPosition().y()/cm;
//     statepost[2] = g4eTrajState.GetPosition().z()/cm;
//     statepost[3] = g4eTrajState.GetMomentum().x()/GeV;
//     statepost[4] = g4eTrajState.GetMomentum().y()/GeV;
//     statepost[5] = g4eTrajState.GetMomentum().z()/GeV;
//     statepost[6] = charge;
//     
//     const Matrix<double, 3, 1> Wpost = statepost.segment<3>(3).normalized();
//     const Matrix<double, 3, 1> khat(0., 0., 1.);
//     
//     const Matrix<double, 3, 1> Upost = khat.cross(Wpost).normalized();
// //     const Matrix<double, 3, 1> Upost = Upostpre.normalized();
//     
//     const double ddxt = ddest.head<3>().dot(Upost);
    
//     std::cout <<
    
    
    
//     const Matrix<double, 5, 7> transportJac = transportJacobianBzAdvanced(statepre, thisPathLength, dEdx, mass);
    
//     const Matrix<double, 5, 7> transportJacAlt = transportJacobianBz(statepre, thisPathLength, dEdx, mass);

//     std::cout << "transportJac" << std::endl;
//     std::cout << transportJac << std::endl;
// //     std::cout << "transportJacAlt" << std::endl;
// //     std::cout << transportJacAlt << std::endl;
//     std::cout << "ddest numerical" << std::endl;
//     std::cout << ddest.transpose() << std::endl;
//     std::cout << "ddxt jac = " << transportJac(3, 5) << " ddxt numerical = " << ddxt << std::endl;
    
//     bool isnan = false;
//     for (unsigned int i=0; i<transportJac.rows(); ++i) {
//       for (unsigned int j=0; j<transportJac.cols(); ++j) {
//         if (std::isnan(transportJac(i, j))) {
//           isnan = true;
//         }
//       }
//     }
//     
//     if (isnan) {
//       std::cout << "thisPathLength = " << thisPathLength << " ePre = " << ePre << " ePost = " << ePost << " dEdx = " << dEdx << " mass = " << mass << std::endl;
//       std::cout << "transportJac" << std::endl;
//       std::cout << transportJac << std::endl;
//     }
    
    // transport contribution to error
    g4errorEnd = (transportJac.leftCols<5>()*g4errorEnd*transportJac.leftCols<5>().transpose()).eval();
    
    dQ = (transportJac.leftCols<5>()*dQ*transportJac.leftCols<5>().transpose()).eval();
    
//     const G4ErrorMatrix &transferMatrix = g4eTrajState.GetTransfMat();
    
    // transport contribution to error
//     g4errorEnd = g4errorEnd.similarity(transferMatrix).T();
//     dQ = dQ.similarity(transferMatrix).T();
    
    // MS and ionization contributions to error
    
//     G4ErrorTrajErr errMSI = g4eTrajState.GetError();
//     
//     // recompute ionization contribution (without truncation for the moment)
//     errMSI(0+1, 0+1) = computeErrorIoni(g4eTrajState.GetG4Track());
//     
//     Matrix<double, 5, 5> errMSIout = Matrix<double, 5, 5>::Zero();
//     
// //     std::cout << "errMSI" << std::endl;
// //     std::cout << errMSI << std::endl;
//     
//     for (unsigned int i = 0; i < 5; i++) {
//       for (unsigned int j = 0; j < 5; j++) {
//         double w = 1.;
//         if (i==0) {
//           w *= charge;
//         }
//         if (j==0) {
//           w *= charge;
//         }
//         errMSIout(i, j) += w*errMSI(i+1, j+1);
//       }
//     }
    
    Matrix<double, 5, 5> errMSIout = PropagateErrorMSC(g4eTrajState.GetG4Track(), pforced);
    errMSIout(0, 0) = computeErrorIoni(g4eTrajState.GetG4Track(), pforced);
      
      
    
    g4errorEnd += errMSIout;
    
    

    
    if (std::abs(thisPathLength) > 0.) {
      dErrorDxLast = errMSIout/thisPathLength;
    }
    
    // transport  + nominal energy loss contribution to jacobian
//     jac = transferMatrix*jac;
    jac = (transportJac.leftCols<5>()*jac).eval();

    //TODO assess relevance of approximations (does the order matter? position/momentum before or after step?)
    //b-field and material contributions to jacobian
    jac.rightCols<2>() += transportJac.rightCols<2>();
        
    
//     const double betaPre = pPre/ePre;    
//     G4ErrorTrajErr errMS = errMSI;
//     errMS(0+1, 0+1) = 0.;
// //     std::cout << "ppre = " << pPre << " betaPre = " << betaPre << " mass = " << mass << " ePre  = " << ePre << std::endl;
// //     G4ErrorTrajErr dQpre = 2.*pPre*(betaPre*betaPre + 2.*mass*mass/ePre/ePre)*errMS;
//     G4ErrorTrajErr dQpre = (2.*pPre + 4.*mass*mass/pPre)*betaPre*betaPre*errMS;
//     
//     for (unsigned int i = 0; i < 5; i++) {
//       for (unsigned int j = 0; j < 5; j++) {
//         double w = charge;
//         if (i==0) {
//           w *= charge;
//         }
//         if (j==0) {
//           w *= charge;
//         }
//         w *= jac(0,0);
//         dQ(i, j) += w*dQpre(i+1, j+1);
//       }
//     }
    
//     Matrix<double, 5, 5> errMS = errMSIout;
//     errMS(0,0) = 0.;
    
//     Matrix<double, 5, 5> errIoni = Matrix<double, 5, 5>::Zero();
//     errIoni(0,0) = errMSIout(0,0);
//     
//     const double q = charge;
//     const double qop = charge/pPre;
//     const double eMass = 510.99906e-6;
//                                             
//     const double x0 = std::pow(q, 2);
//     const double x1 = std::pow(mass, 2)*std::pow(qop, 2);
//     const double x2 = 1./(qop*(x0 + x1));
//     const Matrix<double, 5, 5> dQms = 2.*errMS*x2*(x0 + 2.*x1);
//     const Matrix<double, 5, 5> dQion = errIoni*x2*(4.0*x0 + 10.0*x1);



//     std::cout << "dQionNorm" << std::endl;
//     std::cout << dQion/errIoni(0,0) << std::endl;

    
//      const double dQmsNorm = 2*x0*(x1 + 2*x4)/(x1 + x4);
//      std::cout << "dQmsNorm = " << dQmsNorm << std::endl;

//     std::cout << "errMSIout" << std::endl;
//     std::cout << errMSIout << std::endl;
//     
//     std::cout << "dQms" << std::endl;
//     std::cout << dQms << std::endl;
//     
//     std::cout << "dQion" << std::endl;
//     std::cout << dQion << std::endl;
    
    



    
//     dQ += 0.*jac(0,0)*dQion;
//     dQ += jac(0,0)*dQms;
    
    
    Matrix<double, 5, 5> errMS = errMSIout;
    errMS(0,0) = 0.;
    
    dQ += errMS;
    
    
    
    
    
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
      return retDefault();
    }

    if (theG4eManager->GetPropagator()->CheckIfLastStep(g4eTrajState.GetG4Track())) {
      theG4eManager->GetPropagator()->InvokePostUserTrackingAction(g4eTrajState.GetG4Track());
      continuePropagation = false;
    }
  }

//   std::cout << "propgationDirection = " << propagationDirection() << "iterations = " << iterations << " finalPathLength = " << finalPathLength << std::endl;
//   if (finalPathLength < 0.) {
//     throw std::runtime_error("negative path length!");
//   }
  
  // CMSSW Tracking convention, backward propagations have negative path length
  if (propagationDirection() == oppositeToMomentum)
    finalPathLength = -finalPathLength;

//   std::cout << "finalPathLength = " << finalPathLength << std::endl;
  
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
  ftsEnd[0] = g4eTrajState.GetPosition().x()/cm;
  ftsEnd[1] = g4eTrajState.GetPosition().y()/cm;
  ftsEnd[2] = g4eTrajState.GetPosition().z()/cm;
  ftsEnd[3] = g4eTrajState.GetMomentum().x()/GeV;
  ftsEnd[4] = g4eTrajState.GetMomentum().y()/GeV;
  ftsEnd[5] = g4eTrajState.GetMomentum().z()/GeV;
  ftsEnd[6] = charge;
  
  cmsField->SetOffset(0., 0., 0.);
  cmsField->SetMaterialOffset(0.);
  
  return std::tuple<bool, Matrix<double, 7, 1>, Matrix<double, 5, 5>, Matrix<double, 5, 7>, double, Matrix<double, 5, 5>>(true, ftsEnd, g4errorEnd, jac, dEdxlast, dQ);
  
//   return std::tuple<TrajectoryStateOnSurface, AlgebraicMatrix57, AlgebraicMatrix55, double>(TrajectoryStateOnSurface(localparms, localerr, pDest, theField, side), jacfinal, dQfinal, dEdxlast);
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
Eigen::Matrix<double, 5, 5> Geant4ePropagator::PropagateErrorMSC( const G4Track* aTrack, double pforced) const
{ 
  G4ThreeVector vpPre = aTrack->GetMomentum()/GeV;
//   G4double pPre = vpPre.mag();
//   G4double pBeta = pPre*pPre / (aTrack->GetTotalEnergy()/GeV);
  G4double mass  = aTrack->GetDynamicParticle()->GetMass() / GeV;
  G4double pPre = pforced > 0. ? pforced : aTrack->GetMomentum().mag() / GeV;
  G4double Etot = sqrt(pPre*pPre + mass*mass);
  G4double beta = pPre/Etot;
  G4double pBeta = pPre*beta;
  G4double  stepLengthCm = aTrack->GetStep()->GetStepLength()/cm;

  G4Material* mate = aTrack->GetVolume()->GetLogicalVolume()->GetMaterial();
  G4double effZ, effA;
  CalculateEffectiveZandA( mate, effZ, effA );

#ifdef G4EVERBOSE
  if( iverbose >= 4 ) G4cout << "material " << mate->GetName() 
                     //<< " " << mate->GetZ() << " "  << mate->GetA() 
                        << " effZ:" << effZ << " effA:" << effA
                        << " dens(g/mole):"  << mate->GetDensity()/g*mole << " Radlen/cm:" << mate->GetRadlen()/cm << " nuclLen/cm" << mate->GetNuclearInterLength()/cm << G4endl;
#endif

  G4double RI = stepLengthCm / (mate->GetRadlen()/cm);
#ifdef G4EVERBOSE
  if( iverbose >= 4 ) G4cout << std::setprecision(6) << std::setw(6) << "G4EP:MSC: RI=X/X0 " << RI << " stepLengthCm " << stepLengthCm << " radlen/cm " << (mate->GetRadlen()/cm) << " RI*1.e10:" << RI*1.e10 << G4endl;
#endif
  G4double charge = aTrack->GetDynamicParticle()->GetCharge();
  G4double DDold = 1.8496E-4*RI*(charge/pBeta * charge/pBeta );
  G4double X0 = mate->GetRadlen()/cm;
  G4double Xs = X0*(effZ + 1.)*std::log(287./std::sqrt(effZ))/std::log(159.*std::pow(effZ, -1./3.))/effZ;
  G4double DD = 2.25e-4*stepLengthCm*(charge/pBeta * charge/pBeta )/Xs;
//   std::cout << "DDold = " << DDold << " DD = " << DD << std::endl;
#ifdef G4EVERBOSE
  if( iverbose >= 3 ) G4cout << "G4EP:MSC: D*1E6= " << DD*1.E6 <<" pBeta " << pBeta << G4endl;
#endif
  G4double S1 = DD*stepLengthCm*stepLengthCm/3.;
  G4double S2 = DD;
  G4double S3 = DD*stepLengthCm/2.;

  G4double CLA = std::sqrt( vpPre.x() * vpPre.x() + vpPre.y() * vpPre.y() )/pPre;
#ifdef G4EVERBOSE
  if( iverbose >= 2 ) G4cout << std::setw(6) << "G4EP:MSC: RI " << RI << " S1 " << S1 << " S2 "  << S2 << " S3 "  << S3 << " CLA " << CLA << G4endl;
#endif
  Eigen::Matrix<double, 5, 5> res = Eigen::Matrix<double, 5, 5>::Zero();
  res(1, 1) = S2;
  res(1, 4) = -S3;
  res(2, 2) = S2/CLA/CLA;
  res(2, 3) = S3/CLA;
  res(3, 3) = S1;
  res(4, 4) = S1;
  
  res(4, 1) = res(1, 4);
  res(3, 2) = res(2, 3);

#ifdef G4EVERBOSE
  if( iverbose >= 2 ) G4cout << "G4EP:MSC: error matrix propagated msc " << fError << G4endl;
#endif

  return res;
}

//------------------------------------------------------------------------
double Geant4ePropagator::computeErrorIoni(const G4Track* aTrack, double pforced) const
{
  G4double stepLengthCm = aTrack->GetStep()->GetStepLength() / cm;
#ifdef G4EVERBOSE
  G4double DEDX2;
  if(stepLengthCm < 1.E-7)
  {
    DEDX2 = 0.;
  }
#endif
  //  *     Calculate xi factor (KeV).
  G4Material* mate = aTrack->GetVolume()->GetLogicalVolume()->GetMaterial();
  G4double effZ, effA;
  CalculateEffectiveZandA(mate, effZ, effA);

//   G4double Etot  = aTrack->GetTotalEnergy() / GeV;
//   G4double beta  = aTrack->GetMomentum().mag() / GeV / Etot;
  
  G4double mass  = aTrack->GetDynamicParticle()->GetMass() / GeV;
  G4double pPre = pforced > 0. ? pforced : aTrack->GetMomentum().mag() / GeV;
  G4double Etot = sqrt(pPre*pPre + mass*mass);
  G4double beta = pPre/Etot;
  G4double gamma = Etot / mass;

  // *     Calculate xi factor (keV).
  G4double XI = 153.5 * effZ * stepLengthCm * (mate->GetDensity() / mg * mole) /
                (effA * beta * beta);
                

#ifdef G4EVERBOSE
  if(iverbose >= 2)
  {
    G4cout << "G4EP:IONI: XI/keV " << XI << " beta " << beta << " gamma "
           << gamma << G4endl;
    G4cout << " density " << (mate->GetDensity() / mg * mole) << " effA "
           << effA << " step " << stepLengthCm << G4endl;
  }
#endif
  // *     Maximum energy transfer to atomic electron (KeV).
  G4double eta       = beta * gamma;
  G4double etasq     = eta * eta;
  G4double eMass     = 0.51099906 / GeV;
  G4double massRatio = eMass / mass;
  G4double F1        = 2 * eMass * etasq;
  G4double F2        = 1. + 2. * massRatio * gamma + massRatio * massRatio;
  G4double Emax      = 1.E+6 * F1 / F2;  // now in keV

//   std::cout << "eMass = " << eMass << " mass = " << mass << std::endl;
  //  * *** and now sigma**2  in GeV
//   G4double dedxSq =
//     XI * Emax * (1. - (beta * beta / 2.)) * 1.E-12;  // now in GeV^2
  /*The above  formula for var(1/p) good for dens scatterers. However, for MIPS
    passing through a gas it leads to overestimation. Further more for incident
    electrons the Emax is almost equal to incident energy. This leads  to
    k=Xi/Emax  as small as e-6  and gradually the cov matrix  explodes.
    http://www2.pv.infn.it/~rotondi/kalman_1.pdf
    Since I do not have enough info at the moment to implement Landau &
    sub-Landau models for k=Xi/Emax <0.01 I'll saturate k at this value for now
  */
  

//   const double ePre = aTrack->GetStep()->GetPreStepPoint()->GetTotalEnergy() / GeV;
//   const double ePost = aTrack->GetStep()->GetPostStepPoint()->GetTotalEnergy() / GeV;
//   const double dEdxactual = -(ePost - ePre)/stepLengthCm;
//   
// //   const double poti = 16.e-9 * 10.75;                 // = 16 eV * Z**0.9, for Si Z=14
// //   const double eplasma = 28.816e-9 * sqrt(2.33 * 0.498);  // 28.816 eV * sqrt(rho*(Z/A)) for Si
//   const double poti = 16.e-9 * pow(effZ, 0.9);                // = 16 eV * Z**0.9, for Si Z=14
//   const double density = mate->GetDensity() / mg * mole;
//   const double eplasma = 28.816e-9 * sqrt(density*effZ/effA);  // 28.816 eV * sqrt(rho*(Z/A)) for Si
// //   const double delta0 = 2 * log(eplasma / poti) - 1.;
//   const double delta0 = 2 * log(eplasma / poti) - 1. + 2.*log(beta*gamma);
//   const double beta2 = beta*beta;
// //   const double dEdx = XI*1e-6 * (log(2. * eMass * Emax*1e-6 / (poti * poti)) - 2. * (beta2)-delta0)/stepLengthCm;
//   const double dEdx = XI*1e-6 * (log(2. * eMass * Emax*1e-6*beta2*gamma*gamma / (poti * poti)) - 2. * (beta2)-delta0)/stepLengthCm;
// 
//   std::cout << "material = " << mate->GetName() <<" ePre = " << ePre << " dEdx simple = " << dEdx << " dEdx actual = " << dEdxactual << " actual/simple = " << dEdxactual/dEdx << std::endl;
  
//   dedxSq *= 0.;
  
  // Implementation based on PANDA Report PV/01-07 Section 7
  // TODO understand implications for thick scatterer split into small steps

#if 0
  G4double dedxSq;
  
  const double kappa = XI/Emax;
//   std::cout << "steplength = " << stepLengthCm << " xi = " << XI*1e-6 << " xi/dx = " << XI*1e-6/stepLengthCm << " emax = " << Emax*1e-6 << " kappa = " << kappa << " effZ = " << effZ << std::endl;
  
  if (kappa > 0.005) {
    //vavilov distribution (or gaussian limit which is equivalent for the variance)
  
    dedxSq = XI * Emax * (1. - (beta * beta / 2.)) * 1.E-12;  // now in GeV^2
    
    std::cout << "vavilov: dedxSq = " << dedxSq << std::endl;
  }
  else {
  
  
    const double I = 16.*pow(effZ,0.9);
    const double f2 = effZ <= 2. ? 0. : 2./effZ;
    const double f1 = 1. - f2;
    const double e2 = 10.*effZ*effZ;
    const double e1 = pow(I/pow(e2,f2),1./f1);
    const double r = 0.4;
//     const double emaxev = Emax*1e6;
    const double emaxev = Emax*1e3; // keV -> eV
    const double massev = mass*1e9;
    
    const double ePre = aTrack->GetStep()->GetPreStepPoint()->GetTotalEnergy() / GeV;
    const double ePost = aTrack->GetStep()->GetPostStepPoint()->GetTotalEnergy() / GeV;
    const double C = (ePre-ePost)/stepLengthCm*1e9;
    
    const double sigma1 = C*f1*(log(2.*massev*beta*beta*gamma*gamma/e1) - beta*beta)/e1/(log(2.*massev*beta*beta*gamma*gamma/I) - beta*beta)*(1.-r);
    
    const double sigma2 = C*f1*(log(2.*massev*beta*beta*gamma*gamma/e2) - beta*beta)/e2/(log(2.*massev*beta*beta*gamma*gamma/I) - beta*beta)*(1.-r);
    
    const double sigma3 = C*emaxev/I/(emaxev+I)/log((emaxev+I)/I)*r;
    
    const double Nc = (sigma1 + sigma2 + sigma3)*stepLengthCm;
    
//     std::cout << "Nc = " << Nc << std::endl;
    
    if (Nc > 50.) {
//     if (false) {
      //landau
      constexpr double sigalpha = 15.76; //corresponds to 0.996 quantile
//       constexpr double sigalpha = 22.33; //corresponds to 0.998 quantile
//       constexpr double sigalpha = 31.59; //corresponds to 0.999 quantile
      
      dedxSq = sigalpha*sigalpha*XI*XI*1e-12;
      
      const double dedxsqvavilov = XI * Emax * (1. - (beta * beta / 2.)) * 1.E-12;
      
      std::cout << "dedxsq: vavilov = " << dedxsqvavilov << " landau = " << dedxSq << std::endl;
//       std::cout << "landau\n";
    }
    else {
      //sub-landau
      const double alpha = 0.996;
//       const double alpha = 0.998;
//       const double alpha = 0.999;
//       const double alpha = 1.;
      const double ealpha = I/(1. - alpha*emaxev/(emaxev + I));
      const double e3 = I*(emaxev +I)*log(ealpha/I)/emaxev;
      const double e3sq = I*(emaxev + I)*(ealpha - I)/emaxev;
      const double sigmae3sq = e3sq  - e3*e3;
      
      dedxSq = sigma1*stepLengthCm*e1*e1 + sigma2*stepLengthCm*e2*e2 + sigma3*stepLengthCm*e3*e3 + sigma3*stepLengthCm*sigmae3sq*(sigma3*stepLengthCm + 1.);
      dedxSq *= 1e-18;
      
      const double dedxsqlandau = 15.76*15.76*XI*XI*1e-12;
      const double dedxsqvavilov = XI * Emax * (1. - (beta * beta / 2.)) * 1.E-12;
//       
      std::cout << "dedxsq:  vavilov = " << dedxsqvavilov << " landau = " << dedxsqlandau << " sublandau = " << dedxSq << " Emax = " << Emax << " ealpha = " << ealpha << std::endl;;
//       std::cout << "sublandau\n";
      
    }
    
  
  }
#endif
  
  G4double Emaxmev = Emax*1e-3;
  G4double stepLengthmm = stepLengthCm*10;
  const double ePre = aTrack->GetStep()->GetPreStepPoint()->GetKineticEnergy() / GeV;
  const double ePost = aTrack->GetStep()->GetPostStepPoint()->GetKineticEnergy() / GeV;
  G4double eav = (ePre-ePost)*1e3;
//   G4double epremev = ePre*1e3;
  G4double ekinmev = 0.5*(ePre + ePost)*1e3;

//   std::cout << "eav = " << eav << std::endl;

  G4double dedxsqurban = 1e-6*fluct->SampleFluctuations(mate, aTrack->GetDynamicParticle(), Emaxmev, stepLengthmm, ekinmev);
//   G4double dispurban = fluct->Dispersion(mate, aTrack->GetDynamicParticle(), Emaxmev, stepLengthmm);
//   G4double dispurbansq = 1e-6*dispurban*dispurban;
  
  

//   std::cout << "dedxSq = " << dedxSq << " dedxsqurban = " << dedxsqurban << std::endl;
  
  G4double dedxSq = dedxsqurban;


//   dedxSq = 1.7*1.7*XI*XI*1e-12;
//   dedxSq = 15.76*15.76*XI*XI*1e-12; //corresponds to 0.996 quantile

//   if (true) {
//     if(XI / Emax < 0.01)
//       dedxSq *=
//         XI / Emax * 100;  // Quench for low Elos, see above: newVar=odVar *k/0.01
//   }
  
//   std::cout << "Geant4ePropagator xi = " << XI << " emax = " << Emax << " dedxSq = " << dedxSq << std::endl;

#ifdef G4EVERBOSE
  if(iverbose >= 2)
    G4cout << "G4EP:IONI: DEDX^2(GeV^2) " << dedxSq << " emass/GeV: " << eMass
           << " Emax/keV: " << Emax << "  k=Xi/Emax=" << XI / Emax << G4endl;

#endif
           
  Etot  = aTrack->GetTotalEnergy() / GeV;

  G4double pPre6 =
    (aTrack->GetStep()->GetPreStepPoint()->GetMomentum() / GeV).mag();
  pPre6 = std::pow(pPre6, 6);
  // Apply it to error
  const double res = Etot * Etot * dedxSq / pPre6;
#ifdef G4EVERBOSE
  if(iverbose >= 2)
    G4cout << "G4:IONI Etot/GeV: " << Etot << " err_dedx^2/GeV^2: " << dedxSq
           << " p^6: " << pPre6 << G4endl;
  if(iverbose >= 2)
    G4cout << "G4EP:IONI: error2_from_ionisation "
           << (Etot * Etot * dedxSq) / pPre6 << G4endl;
#endif

  return res;
}

//------------------------------------------------------------------------
void Geant4ePropagator::CalculateEffectiveZandA(const G4Material* mate,
                                                   G4double& effZ,
                                                   G4double& effA)
{
  effZ = 0.;
  effA = 0.;
  G4int ii, nelem = mate->GetNumberOfElements();
  const G4double* fracVec = mate->GetFractionVector();
  for(ii = 0; ii < nelem; ii++)
  {
    effZ += mate->GetElement(ii)->GetZ() * fracVec[ii];
    effA += mate->GetElement(ii)->GetA() * fracVec[ii] / g * mole;
  }
}

Eigen::Matrix<double, 5, 7> Geant4ePropagator::transportJacobianBzD(const Eigen::Matrix<double, 7, 1> &start, double s, double dEdx, double mass, double dBz) const {
  
  if (s==0.) {
    Eigen::Matrix<double, 5, 7> res;
    res.leftCols<5>() = Eigen::Matrix<double, 5, 5>::Identity();
    res.rightCols<2>() = Eigen::Matrix<double, 5, 2>::Zero();
    return res;
  }
  
  const GlobalPoint pos(start[0], start[1], start[2]);  
  const GlobalVector &bfield = theField->inInverseGeV(pos);
  
  const double Bx = bfield.x();
  const double By = bfield.y();
  const double Bz = bfield.z() + 2.99792458e-3*dBz;
  
  const double M0x = start[0];
  const double M0y = start[1];
  const double M0z = start[2];
  
  const Eigen::Matrix<double, 3, 1> W0 = start.segment<3>(3).normalized();
  const double W0x = W0[0];
  const double W0y = W0[1];
  const double W0z = W0[2];
  
  const double q = start[6];
  const double qop0 = q/start.segment<3>(3).norm();

  const double x0 = std::pow(mass, 2);
  const double x1 = std::pow(qop0, -2);
  const double x2 = std::sqrt(std::pow(q, 2)*x1 + x0);
  const double x3 = dEdx*s;
  const double x4 = x2 + x3;
  const double x5 = std::pow(-x0 + std::pow(x4, 2), -3.0/2.0);
  const double x6 = std::pow(W0z, 2);
  const double x7 = std::pow(W0x, 2);
  const double x8 = std::pow(W0y, 2);
  const double x9 = x7 + x8;
  const double x10 = 1.0/x9;
  const double x11 = std::pow(x10*x6 + 1, -1.0/2.0);
  const double x12 = std::pow(Bz, 2);
  const double x13 = std::pow(Bx, 2) + std::pow(By, 2) + x12;
  const double x14 = std::pow(x13, 5.0/2.0);
  const double x15 = std::sqrt(x13);
  const double x16 = s*x15;
  const double x17 = qop0*x16;
  const double x18 = std::sin(x17);
  const double x19 = x14*x18;
  const double x20 = x11*x19;
  const double x21 = std::sqrt(x9);
  const double x22 = 1.0/x21;
  const double x23 = W0z*x22;
  const double x24 = Bx*W0y;
  const double x25 = By*W0x;
  const double x26 = x22*x24 - x22*x25;
  const double x27 = std::pow(x13, 2);
  const double x28 = std::cos(x17);
  const double x29 = x28 - 1;
  const double x30 = x27*x29;
  const double x31 = x11*x30;
  const double x32 = x26*x31;
  const double x33 = std::pow(x13, 3);
  const double x34 = M0x*W0x;
  const double x35 = M0y*W0y;
  const double x36 = M0z*W0z + x34 + x35;
  const double x37 = M0z*(x22*x7 + x22*x8) - x23*x34 - x23*x35;
  const double x38 = W0z*x36 + x21*x37;
  const double x39 = x33*x38;
  const double x40 = x17 - x18;
  const double x41 = std::pow(x13, 3.0/2.0);
  const double x42 = Bx*x11;
  const double x43 = W0x*x22;
  const double x44 = By*x11;
  const double x45 = W0y*x22;
  const double x46 = Bz*x11;
  const double x47 = x23*x46;
  const double x48 = x42*x43 + x44*x45 + x47;
  const double x49 = x41*x48;
  const double x50 = x40*x49;
  const double x51 = Bz*x50 + qop0*x39 + x20*x23 + x32;
  const double x52 = 1.0/x33;
  const double x53 = x1*x52;
  const double x54 = x51*x53;
  const double x55 = x28*x33;
  const double x56 = x11*x55;
  const double x57 = x23*x56;
  const double x58 = x20*x26;
  const double x59 = s*x15 - x16*x28;
  const double x60 = Bz*x49;
  const double x61 = 1.0/qop0;
  const double x62 = x52*x61;
  const double x63 = x62*(s*x57 - s*x58 + x39 + x59*x60);
  const double x64 = x30*x48;
  const double x65 = -Bz*x64 + W0z*x11*x22*x28*x33 - x58;
  const double x66 = x52*x65;
  const double x67 = x21*x36;
  const double x68 = x33*x67;
  const double x69 = W0y*x68;
  const double x70 = W0y*x20;
  const double x71 = x23*x42;
  const double x72 = x43*x46;
  const double x73 = x71 - x72;
  const double x74 = x21*x30;
  const double x75 = -M0x*x45 + M0y*W0x*x22;
  const double x76 = W0z*x37;
  const double x77 = W0x*x75 - W0y*x76;
  const double x78 = x33*x77;
  const double x79 = By*x21;
  const double x80 = qop0*x69 + qop0*x78 + x50*x79 + x70 - x73*x74;
  const double x81 = x22*x53;
  const double x82 = s*x56;
  const double x83 = x19*x73;
  const double x84 = s*x21;
  const double x85 = x49*x59;
  const double x86 = W0y*x82 + x69 + x78 + x79*x85 + x83*x84;
  const double x87 = x22*x62;
  const double x88 = x45*x56;
  const double x89 = -By*x64 + x83 + x88;
  const double x90 = x52*x89;
  const double x91 = W0x*x68;
  const double x92 = W0x*x20;
  const double x93 = x23*x44;
  const double x94 = x45*x46;
  const double x95 = x93 - x94;
  const double x96 = W0x*x76 + W0y*x75;
  const double x97 = x33*x96;
  const double x98 = Bx*x21;
  const double x99 = qop0*x91 - qop0*x97 + x50*x98 + x74*x95 + x92;
  const double x100 = x19*x95;
  const double x101 = W0x*x82 - x100*x84 + x85*x98 + x91 - x97;
  const double x102 = -Bx*x64 + W0x*x11*x22*x28*x33 - x100;
  const double x103 = x102*x52;
  const double x104 = -x103*(x101*x87 - x81*x99) - x66*(-x54 + x63) - x90*(-x80*x81 + x86*x87);
  const double x105 = q*x4*x5;
  const double x106 = dEdx*x105;
  const double x107 = W0z*x10;
  const double x108 = W0x*x107;
  const double x109 = W0y*x107;
  const double x110 = Bz*x11 - x108*x42 - x109*x44;
  const double x111 = Bz*x110;
  const double x112 = x40*x41;
  const double x113 = x111*x112 + x20 - x23*x32;
  const double x114 = std::pow(x13, -6);
  const double x115 = x114*x61;
  const double x116 = x115*x65;
  const double x117 = W0z*x43;
  const double x118 = x109*x46 + x44;
  const double x119 = x110*x112*x98 - x117*x20 + x118*x74;
  const double x120 = x115*x22;
  const double x121 = x102*x120;
  const double x122 = W0z*x45;
  const double x123 = x108*x46 + x42;
  const double x124 = By*x110*x21*x40*x41 - x122*x20 - x123*x74;
  const double x125 = x120*x89;
  const double x126 = -x113*x116 - x119*x121 - x124*x125;
  const double x127 = Bx*W0x;
  const double x128 = By*W0y;
  const double x129 = x127*x22 + x128*x22;
  const double x130 = By*W0x*x11*x22 - x42*x45;
  const double x131 = Bz*x130;
  const double x132 = x112*x131 + x129*x31;
  const double x133 = x30*x46;
  const double x134 = -W0y*x133 + x112*x130*x79 + x92;
  const double x135 = Bx*x130*x21*x40*x41 - W0x*x133 - x70;
  const double x136 = -x116*x132 - x121*x135 - x125*x134;
  const double x137 = W0y*x102*x22*x52 - x43*x90;
  const double x138 = x103*x117 + x122*x90 - x21*x66;
  const double x139 = 6*Bz;
  const double x140 = x139/std::pow(x13, 4);
  const double x141 = x140*x61;
  const double x142 = x141*x51;
  const double x143 = x18*x41;
  const double x144 = 5*x143;
  const double x145 = x27*x28;
  const double x146 = qop0*s;
  const double x147 = x145*x146;
  const double x148 = x13*x29;
  const double x149 = 4*x148;
  const double x150 = x26*x46;
  const double x151 = x146*x150;
  const double x152 = Bz*qop0;
  const double x153 = 6*x27;
  const double x154 = x152*x153;
  const double x155 = x12*x48;
  const double x156 = 3*x15*x40;
  const double x157 = 1.0/x15;
  const double x158 = s*x152;
  const double x159 = Bz*qop0*s*x157 - x157*x158*x28;
  const double x160 = x112*x47 - x143*x151 + x144*x47 + x147*x47 + x149*x150 + x154*x38 + x155*x156 + x159*x60 + x50;
  const double x161 = x141*x22;
  const double x162 = Bz*W0y;
  const double x163 = qop0*x153*x67;
  const double x164 = x144*x46;
  const double x165 = x147*x46;
  const double x166 = W0z*x112;
  const double x167 = Bz*x73;
  const double x168 = x149*x21;
  const double x169 = x158*x73;
  const double x170 = x143*x21;
  const double x171 = Bz*x48;
  const double x172 = x156*x171;
  const double x173 = x159*x49;
  const double x174 = W0x*x31 + W0y*x164 + W0y*x165 + x154*x77 + x162*x163 + x166*x44 - x167*x168 + x169*x170 + x172*x79 + x173*x79;
  const double x175 = Bz*W0x;
  const double x176 = Bz*x95;
  const double x177 = x158*x95;
  const double x178 = W0x*x164 + W0x*x165 - W0y*x31 - x154*x96 + x163*x175 + x166*x42 + x168*x176 - x170*x177 + x172*x98 + x173*x98;
  const double x179 = -x103*(-x161*x99 + x178*x22*x52*x61) - x66*(-x142 + x160*x52*x61) - x90*(-x161*x80 + x174*x22*x52*x61);
  const double x180 = std::pow(x89, 2);
  const double x181 = std::pow(x102, 2);
  const double x182 = x114*x180 + x114*x181;
  const double x183 = 1.0/x182;
  const double x184 = x114*x183;
  const double x185 = 1.0/(x184*std::pow(x65, 2) + 1);
  const double x186 = s*x18;
  const double x187 = x11*std::pow(x13, 7.0/2.0);
  const double x188 = x187*x23;
  const double x189 = std::pow(x182, -1.0/2.0);
  const double x190 = x189*x52;
  const double x191 = x186*x187;
  const double x192 = x191*x45;
  const double x193 = s*x55;
  const double x194 = x193*x73;
  const double x195 = x19*x48;
  const double x196 = By*s*x195;
  const double x197 = (1.0/2.0)*x114;
  const double x198 = x197*x89;
  const double x199 = x191*x43;
  const double x200 = x193*x95;
  const double x201 = x102*x197;
  const double x202 = x66/std::pow(x182, 3.0/2.0);
  const double x203 = qop0*x18;
  const double x204 = qop0*x56;
  const double x205 = x187*x203;
  const double x206 = x205*x45;
  const double x207 = qop0*x55;
  const double x208 = x207*x73;
  const double x209 = By*qop0*x195;
  const double x210 = x205*x43;
  const double x211 = x207*x95;
  const double x212 = x185*(x190*(Bz*qop0*x14*x18*x48 - x188*x203 - x204*x26) + x202*(-x198*(-2*x206 + 2*x208 + 2*x209) - x201*(2*Bx*qop0*x14*x18*x48 - 2*x210 - 2*x211)));
  const double x213 = x109*x56;
  const double x214 = x110*x30;
  const double x215 = By*x214;
  const double x216 = x108*x56;
  const double x217 = x118*x19;
  const double x218 = Bx*x214;
  const double x219 = x43*x56;
  const double x220 = x19*x94;
  const double x221 = 2*x220;
  const double x222 = x130*x30;
  const double x223 = By*x222;
  const double x224 = Bx*x222;
  const double x225 = x139/std::pow(x13, 7);
  const double x226 = x20*x43;
  const double x227 = x30*x93;
  const double x228 = 10*x143;
  const double x229 = x145*x169;
  const double x230 = By*x171;
  const double x231 = 8*x148;
  const double x232 = By*qop0*x186*x60;
  const double x233 = x30*x71;
  const double x234 = x19*x72;
  const double x235 = x145*x177;
  const double x236 = Bx*x171;
  const double x237 = x102*x184;
  const double x238 = x184*x89;
  const double x239 = x237*(-x206 + x208 + x209) - x238*(Bx*qop0*x14*x18*x48 - x210 - x211);
  const double x240 = x30*(Bz*W0z + x127 + x128);
  const double x241 = -By*x240 + W0y*x55 + x19*(Bx*W0z - x175);
  const double x242 = x6 + x9;
  const double x243 = 1.0/x242;
  const double x244 = x114*x243;
  const double x245 = std::pow(x241, 2)*x244;
  const double x246 = -Bx*x240 + W0x*x28*x33 - x19*(By*W0z - x162);
  const double x247 = x244*std::pow(x246, 2);
  const double x248 = std::pow(x245 + x247, -1.0/2.0);
  const double x249 = x241*x248;
  const double x250 = x1*x99;
  const double x251 = std::pow(x10*x242, -1.0/2.0);
  const double x252 = x10*x251;
  const double x253 = x114*x252;
  const double x254 = x246*x248;
  const double x255 = x1*x80;
  const double x256 = x115*x252;
  const double x257 = x254*x256;
  const double x258 = x249*x256;
  const double x259 = qop0*x21;
  const double x260 = -qop0*x15*x28 + qop0*x15;
  const double x261 = x260*x49;
  const double x262 = W0y*x204 + x259*x83 + x261*x79;
  const double x263 = W0x*x204 - x100*x259 + x261*x98;
  const double x264 = x257*x262 - x258*x263;
  const double x265 = x252*x52;
  const double x266 = x251*x52;
  const double x267 = x225*x252;
  const double x268 = x249*x61;
  const double x269 = x254*x61;
  const double x270 = -Bz*x240 + W0z*x28*x33 - x19*(x24 - x25);
  const double x271 = x22*x243*x270/std::pow(x13, 9);
  const double x272 = x249*x271;
  const double x273 = x272*x61;
  const double x274 = x269*x271;
  const double x275 = x245*x248 + x247*x248;
  const double x276 = -x262*x273 - x263*x274 + x275*x52*x61*(qop0*x57 - qop0*x58 + x260*x60);
  const double x277 = x275*x62;
  const double x278 = x244*x270;
  const double x279 = x249*x278;
  const double x280 = x254*x278;
  const double x281 = x139*x22*x243*x270/std::pow(x13, 10);
  const double dqopdqop0 = std::pow(q, 3)*x4*x5/(std::pow(qop0, 3)*x2) - x104*x106;
  const double dqopdlam0 = -x106*x126;
  const double dqopdphi0 = -x106*x136;
  const double dqopdxt0 = -x106*x137;
  const double dqopdyt0 = -x106*x138;
  const double dqopdBz = -x106*x179;
  const double dqopdxi = -x105*x3;
  const double dlamdqop0 = x104*x212 + x185*(x190*(Bz*s*x14*x18*x48 - x186*x188 - x26*x82) + x202*(-x198*(-2*x192 + 2*x194 + 2*x196) - x201*(2*Bx*s*x14*x18*x48 - 2*x199 - 2*x200)));
  const double dlamdlam0 = x126*x212 + x185*(x190*(-x111*x30 + x23*x58 + x56) + x202*(-x198*(2*x123*x14*x18 - 2*x213 - 2*x215) - x201*(-2*x216 - 2*x217 - 2*x218)));
  const double dlamdphi0 = x136*x212 + x185*(x190*(-x129*x20 - x131*x30) + x202*(-x198*(2*x219 + x221 - 2*x223) - x201*(2*Bz*W0x*x11*x14*x18*x22 - 2*x224 - 2*x88)));
  const double dlamdxt0 = x137*x212;
  const double dlamdyt0 = x138*x212;
  const double dlamdBz = x179*x212 + x185*(-x140*x189*x65 + x190*(6*Bz*W0z*x11*x22*x27*x28 + qop0*s*x12*x18*x41*x48 - x144*x150 - x145*x151 - x146*x19*x47 - x149*x155 - x30*x47 - x64) + x202*(x180*x225 + x181*x225 - x198*(12*x145*x94 - x146*x221 + x167*x228 - 2*x226 - 2*x227 + 2*x229 - x230*x231 + 2*x232) - x201*(2*Bx*Bz*qop0*s*x18*x41*x48 + 12*Bz*W0x*x11*x22*x27*x28 + 2*W0y*x11*x14*x18*x22 - 2*x146*x234 - x176*x228 - x231*x236 - 2*x233 - 2*x235)));
  const double dlamdxi = 0;
  const double dphidqop0 = x104*x239 + x237*(-x192 + x194 + x196) - x238*(Bx*s*x14*x18*x48 - x199 - x200);
  const double dphidlam0 = x126*x239 + x237*(x123*x14*x18 - x213 - x215) - x238*(-x216 - x217 - x218);
  const double dphidphi0 = x136*x239 + x237*(x219 + x220 - x223) - x238*(Bz*W0x*x11*x14*x18*x22 - x224 - x88);
  const double dphidxt0 = x137*x239;
  const double dphidyt0 = x138*x239;
  const double dphidBz = x103*x183*(-x140*x89 + x52*(x144*x167 - x146*x220 - x149*x230 + x153*x28*x94 - x226 - x227 + x229 + x232)) + x179*x239 - x183*x90*(-x102*x140 + x52*(Bx*Bz*qop0*s*x18*x41*x48 + 6*Bz*W0x*x11*x22*x27*x28 + W0y*x11*x14*x18*x22 - x144*x176 - x146*x234 - x149*x236 - x233 - x235));
  const double dphidxi = 0;
  const double dxtdqop0 = -x101*x258 + x104*x264 + x249*x250*x253 - x253*x254*x255 + x257*x86;
  const double dxtdlam0 = -x119*x258 + x124*x257 + x126*x264;
  const double dxtdphi0 = x134*x257 - x135*x258 + x136*x264;
  const double dxtdxt0 = W0x*x254*x265 + W0y*x249*x265 + x137*x264;
  const double dxtdyt0 = x108*x249*x266 - x109*x254*x266 + x138*x264;
  const double dxtdBz = x174*x257 - x178*x258 + x179*x264 + x267*x268*x99 - x267*x269*x80;
  const double dxtdxi = 0;
  const double dytdqop0 = -x101*x274 + x104*x276 + x250*x254*x271 + x255*x272 - x273*x86 - x275*x54 + x275*x63;
  const double dytdlam0 = x113*x277 - x119*x274 - x124*x273 + x126*x276;
  const double dytdphi0 = x132*x277 - x134*x273 - x135*x274 + x136*x276;
  const double dytdxt0 = x137*x276 - x279*x43 + x280*x45;
  const double dytdyt0 = x117*x280 + x122*x279 + x138*x276 + x21*x275;
  const double dytdBz = -x142*x275 + x160*x275*x62 - x174*x273 - x178*x274 + x179*x276 + x268*x281*x80 + x269*x281*x99;
  const double dytdxi = 0;
  Eigen::Matrix<double, 5, 7> res;
  res(0,0) = dqopdqop0;
  res(0,1) = dqopdlam0;
  res(0,2) = dqopdphi0;
  res(0,3) = dqopdxt0;
  res(0,4) = dqopdyt0;
  res(0,5) = dqopdBz;
  res(0,6) = dqopdxi;
  res(1,0) = dlamdqop0;
  res(1,1) = dlamdlam0;
  res(1,2) = dlamdphi0;
  res(1,3) = dlamdxt0;
  res(1,4) = dlamdyt0;
  res(1,5) = dlamdBz;
  res(1,6) = dlamdxi;
  res(2,0) = dphidqop0;
  res(2,1) = dphidlam0;
  res(2,2) = dphidphi0;
  res(2,3) = dphidxt0;
  res(2,4) = dphidyt0;
  res(2,5) = dphidBz;
  res(2,6) = dphidxi;
  res(3,0) = dxtdqop0;
  res(3,1) = dxtdlam0;
  res(3,2) = dxtdphi0;
  res(3,3) = dxtdxt0;
  res(3,4) = dxtdyt0;
  res(3,5) = dxtdBz;
  res(3,6) = dxtdxi;
  res(4,0) = dytdqop0;
  res(4,1) = dytdlam0;
  res(4,2) = dytdphi0;
  res(4,3) = dytdxt0;
  res(4,4) = dytdyt0;
  res(4,5) = dytdBz;
  res(4,6) = dytdxi;

  res.col(5) *= 2.99792458e-3;
  
  return res;

  
}

Eigen::Matrix<double, 5, 7> Geant4ePropagator::transportJacobianBzMidpointD(const Eigen::Matrix<double, 7, 1> &start, const Eigen::Matrix<double, 7, 1> &end, double s, double dEdx, double mass, double dBz) const {

  if (s==0.) {
    Eigen::Matrix<double, 5, 7> res;
    res.leftCols<5>() = Eigen::Matrix<double, 5, 5>::Identity();
    res.rightCols<2>() = Eigen::Matrix<double, 5, 2>::Zero();
    return res;
  }

  const GlobalPoint pos(start[0], start[1], start[2]);
  const GlobalVector &bfield = theField->inInverseGeV(pos);

  const GlobalPoint posend(end[0], end[1], end[2]);
  const GlobalVector &bfieldend = theField->inInverseGeV(posend);

  const double Bx = 0.5*(bfield.x() + bfieldend.x());
  const double By = 0.5*(bfield.y() + bfieldend.y());
  const double Bz = 0.5*(bfield.z() + bfieldend.z()) + 2.99792458e-3*dBz;

  const double M0x = start[0];
  const double M0y = start[1];
  const double M0z = start[2];

  const Eigen::Matrix<double, 3, 1> W0 = start.segment<3>(3).normalized();
  const double W0x = W0[0];
  const double W0y = W0[1];
  const double W0z = W0[2];

  const double q = start[6];
  const double qop0 = q/start.segment<3>(3).norm();

  const double x0 = std::pow(mass, 2);
  const double x1 = std::pow(q, 2);
  const double x2 = std::pow(qop0, 2);
  const double x3 = std::sqrt(x0 + x1/x2);
  const double x4 = dEdx*s;
  const double x5 = x3 + x4;
  const double x6 = std::pow(-x0 + std::pow(x5, 2), -3.0/2.0);
  const double x7 = std::pow(W0x, 2);
  const double x8 = std::pow(W0y, 2);
  const double x9 = x7 + x8;
  const double x10 = std::sqrt(x9);
  const double x11 = 1.0/x10;
  const double x12 = W0z*x11;
  const double x13 = x0*x2;
  const double x14 = std::fabs(qop0);
  const double x15 = x14*x4;
  const double x16 = std::sqrt(x1 + x13);
  const double x17 = x15 + x16;
  const double x18 = std::sqrt(-x13 + std::pow(x17, 2));
  const double x19 = std::pow(Bx, 2);
  const double x20 = std::pow(By, 2);
  const double x21 = std::pow(Bz, 2);
  const double x22 = x19 + x20 + x21;
  const double x23 = std::sqrt(x22);
  const double x24 = s*x23;
  const double x25 = std::pow(s, 2);
  const double x26 = std::pow(dEdx, 2)*x25;
  const double x27 = x2*x26;
  const double x28 = x15*x16;
  const double x29 = x1 + x27 + 2*x28;
  const double x30 = std::pow(x29, -1.0/2.0);
  const double x31 = q*x14;
  const double x32 = qop0 + x30*x31;
  const double x33 = (1.0/2.0)*x32;
  const double x34 = x24*x33;
  const double x35 = std::sin(x34);
  const double x36 = 2*x35;
  const double x37 = x18*x36;
  const double x38 = std::pow(W0z, 2);
  const double x39 = 1.0/x9;
  const double x40 = std::pow(x38*x39 + 1, -1.0/2.0);
  const double x41 = std::pow(x22, 5.0/2.0);
  const double x42 = x40*x41;
  const double x43 = x37*x42;
  const double x44 = qop0*x18 + x31;
  const double x45 = std::pow(x22, 3);
  const double x46 = M0x*W0x;
  const double x47 = M0y*W0y;
  const double x48 = M0z*W0z + x46 + x47;
  const double x49 = M0z*(x11*x7 + x11*x8) - x12*x46 - x12*x47;
  const double x50 = W0z*x48 + x10*x49;
  const double x51 = x45*x50;
  const double x52 = W0y*x11;
  const double x53 = Bx*x52;
  const double x54 = W0x*x11;
  const double x55 = By*x54;
  const double x56 = x53 - x55;
  const double x57 = x18*x56;
  const double x58 = std::pow(x22, 2);
  const double x59 = std::cos(x34);
  const double x60 = x59 - 1;
  const double x61 = 2*x58*x60;
  const double x62 = x40*x61;
  const double x63 = x57*x62;
  const double x64 = s*x23*x44 - x37;
  const double x65 = std::pow(x22, 3.0/2.0);
  const double x66 = Bz*x40;
  const double x67 = x12*x66;
  const double x68 = Bx*x40;
  const double x69 = x54*x68;
  const double x70 = By*x40;
  const double x71 = x52*x70;
  const double x72 = x69 + x71;
  const double x73 = x67 + x72;
  const double x74 = x65*x73;
  const double x75 = x64*x74;
  const double x76 = Bz*x75 + x12*x43 + x44*x51 + x63;
  const double x77 = 1.0/x45;
  const double x78 = std::pow(x44, -2);
  const double x79 = (((qop0) > 0) - ((qop0) < 0));
  const double x80 = q*x79;
  const double x81 = 1.0/x18;
  const double x82 = qop0*x0;
  const double x83 = x4*x79;
  const double x84 = x82/x16;
  const double x85 = x81*((1.0/2.0)*x17*(2*x83 + 2*x84) - x82);
  const double x86 = qop0*x85 + x18 + x80;
  const double x87 = -x86;
  const double x88 = x78*x87;
  const double x89 = x77*x88;
  const double x90 = x36*x85;
  const double x91 = x42*x90;
  const double x92 = x56*x62;
  const double x93 = x40*x59;
  const double x94 = x18*x93;
  const double x95 = std::pow(x29, -3.0/2.0);
  const double x96 = x31*x95;
  const double x97 = x30*x80 + x96*(-qop0*x26 - x15*x84 - x16*x83) + 1;
  const double x98 = s*x97;
  const double x99 = x94*x98;
  const double x100 = x12*x45;
  const double x101 = x42*x57;
  const double x102 = x35*x98;
  const double x103 = x18*x59;
  const double x104 = x103*x24;
  const double x105 = s*x23*x86 - x104*x97 - x90;
  const double x106 = Bz*x74;
  const double x107 = x100*x99 - x101*x102 + x105*x106 + x12*x91 + x51*x86 + x85*x92;
  const double x108 = 1.0/x44;
  const double x109 = x108*x77;
  const double x110 = 1.0/x22;
  const double x111 = Bx*x54;
  const double x112 = x111*x66;
  const double x113 = By*x52;
  const double x114 = x113*x66;
  const double x115 = x21*x40;
  const double x116 = x19*x93;
  const double x117 = x20*x93;
  const double x118 = x23*x35;
  const double x119 = x52*x68;
  const double x120 = -x112*x59 + x112 - x114*x59 + x114 + x115*x12 + x116*x12 + x117*x12 - x118*x119 + x118*x54*x70;
  const double x121 = x110*x120;
  const double x122 = W0y*x45;
  const double x123 = x10*x48;
  const double x124 = x122*x123;
  const double x125 = W0y*x43;
  const double x126 = -M0x*x52 + M0y*W0x*x11;
  const double x127 = W0z*x49;
  const double x128 = W0x*x126 - W0y*x127;
  const double x129 = x128*x45;
  const double x130 = x12*x68;
  const double x131 = x54*x66;
  const double x132 = x130 - x131;
  const double x133 = x10*x18;
  const double x134 = x133*x61;
  const double x135 = By*x10;
  const double x136 = x124*x44 + x125 + x129*x44 - x132*x134 + x135*x75;
  const double x137 = x11*x89;
  const double x138 = x10*x61;
  const double x139 = x138*x85;
  const double x140 = x132*x133;
  const double x141 = x102*x41;
  const double x142 = x105*x74;
  const double x143 = W0y*x91 + x122*x99 + x124*x86 + x129*x86 - x132*x139 + x135*x142 + x140*x141;
  const double x144 = x109*x11;
  const double x145 = W0y*x20;
  const double x146 = x11*x40;
  const double x147 = x111*x70;
  const double x148 = By*x66;
  const double x149 = x12*x148;
  const double x150 = x147*x59;
  const double x151 = W0z*x59;
  const double x152 = x11*x151;
  const double x153 = x148*x152;
  const double x154 = x115*x59;
  const double x155 = x154*x52;
  const double x156 = -x118*x131 + x155;
  const double x157 = x116*x52 + x118*x130 + x145*x146 + x147 + x149 - x150 - x153 + x156;
  const double x158 = x110*x157;
  const double x159 = W0x*x45;
  const double x160 = x123*x159;
  const double x161 = W0x*x43;
  const double x162 = W0x*x127 + W0y*x126;
  const double x163 = x162*x45;
  const double x164 = x12*x70;
  const double x165 = x52*x66;
  const double x166 = x164 - x165;
  const double x167 = Bx*x10;
  const double x168 = x134*x166 + x160*x44 + x161 - x163*x44 + x167*x75;
  const double x169 = x133*x166;
  const double x170 = W0x*x91 + x139*x166 - x141*x169 + x142*x167 + x159*x99 + x160*x86 - x163*x86;
  const double x171 = W0x*x19;
  const double x172 = x53*x70;
  const double x173 = Bx*x66;
  const double x174 = x12*x173;
  const double x175 = x172*x59;
  const double x176 = x154*x54;
  const double x177 = x118*x165 + x176;
  const double x178 = x117*x54 - x118*x164 + x146*x171 - x152*x173 + x172 + x174 - x175 + x177;
  const double x179 = x110*x178;
  const double x180 = -x121*(x107*x109 + x76*x89) - x158*(x136*x137 + x143*x144) - x179*(x137*x168 + x144*x170);
  const double x181 = q*x5*x6;
  const double x182 = dEdx*x181;
  const double x183 = W0z*x39;
  const double x184 = W0x*x183;
  const double x185 = W0y*x183;
  const double x186 = Bz*x40 - x184*x68 - x185*x70;
  const double x187 = x64*x65;
  const double x188 = Bz*x187;
  const double x189 = -x12*x63 + x186*x188 + x43;
  const double x190 = std::pow(x22, -4);
  const double x191 = x108*x190;
  const double x192 = x120*x191;
  const double x193 = x40*x54;
  const double x194 = W0z*x37*x41;
  const double x195 = x185*x66;
  const double x196 = x134*(x195 + x70) + x167*x186*x187 - x193*x194;
  const double x197 = x11*x191;
  const double x198 = x178*x197;
  const double x199 = x40*x52;
  const double x200 = By*x10*x186*x64*x65 - x134*(x184*x66 + x68) - x194*x199;
  const double x201 = x157*x197;
  const double x202 = -x189*x192 - x196*x198 - x200*x201;
  const double x203 = x18*x62;
  const double x204 = By*W0x*x11*x40 - x119;
  const double x205 = x188*x204 + x203*(x111 + x113);
  const double x206 = x18*x66;
  const double x207 = x206*x61;
  const double x208 = Bx*x10*x204*x64*x65 - W0x*x207 - x125;
  const double x209 = -W0y*x207 + x135*x187*x204 + x161;
  const double x210 = -x192*x205 - x198*x208 - x201*x209;
  const double x211 = W0y*x11*x110*x178 - x158*x54;
  const double x212 = W0z*x54;
  const double x213 = W0z*x52;
  const double x214 = -x10*x121 + x158*x213 + x179*x212;
  const double x215 = 6*Bz;
  const double x216 = x191*x215*x76;
  const double x217 = 10*x35*x65;
  const double x218 = 6*x58;
  const double x219 = Bz*x218*x44;
  const double x220 = x57*x66;
  const double x221 = 8*x22*x60;
  const double x222 = s*x32;
  const double x223 = x206*x222*x58;
  const double x224 = x32*x35;
  const double x225 = s*x224*x65;
  const double x226 = 3*x23*x64*x73;
  const double x227 = 1.0/x23;
  const double x228 = s*x227;
  const double x229 = x103*x32;
  const double x230 = Bz*s*x227*x44 - Bz*x228*x229;
  const double x231 = x106*x230 + x152*x223 + x18*x217*x67 + x187*x67 + x21*x226 + x219*x50 + x220*x221 - x220*x225 + x75;
  const double x232 = x197*x215;
  const double x233 = Bz*W0x;
  const double x234 = x123*x218*x44;
  const double x235 = x206*x217;
  const double x236 = W0x*x59;
  const double x237 = Bz*x221;
  const double x238 = W0z*x187;
  const double x239 = Bz*x225;
  const double x240 = Bx*Bz;
  const double x241 = x10*x226;
  const double x242 = x230*x74;
  const double x243 = W0x*x235 - W0y*x203 - x162*x219 + x167*x242 + x169*x237 - x169*x239 + x223*x236 + x233*x234 + x238*x68 + x240*x241;
  const double x244 = Bz*W0y;
  const double x245 = W0y*x59;
  const double x246 = By*Bz;
  const double x247 = W0x*x203 + W0y*x235 + x128*x219 + x135*x242 - x140*x237 + x140*x239 + x223*x245 + x234*x244 + x238*x70 + x241*x246;
  const double x248 = -x121*(x108*x231*x77 - x216) - x158*(x108*x11*x247*x77 - x136*x232) - x179*(x108*x11*x243*x77 - x168*x232);
  const double x249 = x15*x17*x81;
  const double x250 = qop0*x249;
  const double x251 = x250*x78;
  const double x252 = x251*x77;
  const double x253 = x252*x76;
  const double x254 = x249*x36;
  const double x255 = x254*x42;
  const double x256 = -x27 - x28;
  const double x257 = x256*x96;
  const double x258 = s*x257;
  const double x259 = x258*x94;
  const double x260 = x258*x35;
  const double x261 = dEdx*qop0*x14*x17*x23*x25*x81 - x104*x257 - x254;
  const double x262 = x100*x259 - x101*x260 + x106*x261 + x12*x255 + x249*x92 + x250*x51;
  const double x263 = x11*x252;
  const double x264 = x138*x249;
  const double x265 = x260*x41;
  const double x266 = x261*x74;
  const double x267 = W0y*x255 + x122*x259 + x124*x250 + x129*x250 - x132*x264 + x135*x266 + x140*x265;
  const double x268 = W0x*x255 + x159*x259 + x160*x250 - x163*x250 + x166*x264 + x167*x266 - x169*x265;
  const double x269 = -x121*(x108*x262*x77 - x253) - x158*(x108*x11*x267*x77 - x136*x263) - x179*(x108*x11*x268*x77 - x168*x263);
  const double x270 = 1.0/x58;
  const double x271 = std::pow(x178, 2);
  const double x272 = std::pow(x157, 2);
  const double x273 = x270*x271 + x270*x272;
  const double x274 = 1.0/x273;
  const double x275 = x270*x274;
  const double x276 = 1.0/(std::pow(x120, 2)*x275 + 1);
  const double x277 = x22*x98;
  const double x278 = x277*x59;
  const double x279 = x12*x40;
  const double x280 = x19*x279;
  const double x281 = x24*x35;
  const double x282 = x281*x97;
  const double x283 = (1.0/2.0)*x282;
  const double x284 = x20*x279;
  const double x285 = std::pow(x273, -1.0/2.0);
  const double x286 = x110*x285;
  const double x287 = x131*x278;
  const double x288 = x19*x199;
  const double x289 = x282*x288;
  const double x290 = x115*x52;
  const double x291 = x282*x290;
  const double x292 = (1.0/2.0)*x270;
  const double x293 = x157*x292;
  const double x294 = x165*x278;
  const double x295 = x152*x70;
  const double x296 = x277*x295;
  const double x297 = x193*x20;
  const double x298 = x282*x297;
  const double x299 = x115*x54;
  const double x300 = x282*x299;
  const double x301 = x178*x292;
  const double x302 = x121/std::pow(x273, 3.0/2.0);
  const double x303 = x22*x33;
  const double x304 = x303*x59;
  const double x305 = x118*x33;
  const double x306 = x22*x32;
  const double x307 = x306*x59;
  const double x308 = x118*x32;
  const double x309 = x118*x70;
  const double x310 = x276*(x286*((1.0/2.0)*Bx*Bz*W0x*x11*x23*x32*x35*x40 + (1.0/2.0)*By*Bz*W0y*x11*x23*x32*x35*x40 + (1.0/2.0)*By*W0x*x11*x22*x32*x40*x59 - x119*x304 - x280*x305 - x284*x305) + x302*(-x293*(Bx*By*W0x*x11*x23*x32*x35*x40 + Bx*W0z*x11*x22*x32*x40*x59 + By*Bz*W0z*x11*x23*x32*x35*x40 - x131*x307 - x288*x308 - x290*x308) - x301*(x165*x307 + x174*x308 - x295*x306 - x297*x308 - x299*x308 + x309*x32*x53)));
  const double x311 = Bx*W0y;
  const double x312 = x183*x311*x70;
  const double x313 = x183*x40;
  const double x314 = x171*x313;
  const double x315 = x23*x36;
  const double x316 = x117*x184;
  const double x317 = x115*x183;
  const double x318 = x236*x317;
  const double x319 = Bx*x184*x70;
  const double x320 = x145*x313;
  const double x321 = x116*x185;
  const double x322 = x245*x317;
  const double x323 = x53*x66;
  const double x324 = x55*x66;
  const double x325 = x324*x59;
  const double x326 = x11*x93;
  const double x327 = x145*x326;
  const double x328 = x171*x326;
  const double x329 = 2*Bz;
  const double x330 = x270*x329;
  const double x331 = x227*x35;
  const double x332 = s*x33;
  const double x333 = x228*x33*x35;
  const double x334 = x21*x333;
  const double x335 = x333*x67;
  const double x336 = x329*x77;
  const double x337 = x152*x68;
  const double x338 = x227*x36;
  const double x339 = x224*x228;
  const double x340 = x148*x53;
  const double x341 = std::pow(Bz, 3);
  const double x342 = x339*x341;
  const double x343 = W0z*x21;
  const double x344 = x11*x343;
  const double x345 = x131*x20;
  const double x346 = x165*x19;
  const double x347 = x22*x258;
  const double x348 = x347*x59;
  const double x349 = x257*x281;
  const double x350 = (1.0/2.0)*x349;
  const double x351 = x131*x348;
  const double x352 = x288*x349;
  const double x353 = x290*x349;
  const double x354 = x165*x348;
  const double x355 = x295*x347;
  const double x356 = x297*x349;
  const double x357 = x299*x349;
  const double x358 = x178*x275;
  const double x359 = x157*x275;
  const double x360 = x358*((1.0/2.0)*Bx*By*W0x*x11*x23*x32*x35*x40 + (1.0/2.0)*Bx*W0z*x11*x22*x32*x40*x59 + (1.0/2.0)*By*Bz*W0z*x11*x23*x32*x35*x40 - x131*x304 - x288*x305 - x290*x305) - x359*(x165*x304 + x174*x305 - x295*x303 - x297*x305 - x299*x305 + x309*x33*x53);
  const double x361 = W0y*x39;
  const double x362 = W0z*x118;
  const double x363 = x333*x341;
  const double x364 = x136*x88;
  const double x365 = Bx*By;
  const double x366 = W0y*x365;
  const double x367 = Bz*W0z;
  const double x368 = Bx*x367;
  const double x369 = -By*x362 + x118*x244 + x171 + x20*x236 + x21*x236 - x366*x59 + x366 - x368*x59 + x368;
  const double x370 = x38 + x9;
  const double x371 = 1.0/x370;
  const double x372 = x270*x371;
  const double x373 = std::pow(x369, 2)*x372;
  const double x374 = By*x367;
  const double x375 = Bx*x362 + W0x*x365 - x118*x233 + x145 + x19*x245 + x21*x245 - x236*x365 - x374*x59 + x374;
  const double x376 = x372*std::pow(x375, 2);
  const double x377 = std::pow(x373 + x376, -1.0/2.0);
  const double x378 = x369*x377;
  const double x379 = std::pow(x370*x39, -1.0/2.0);
  const double x380 = x379*x39;
  const double x381 = x190*x380;
  const double x382 = x378*x381;
  const double x383 = x375*x377;
  const double x384 = x168*x383;
  const double x385 = x381*x384;
  const double x386 = x191*x380;
  const double x387 = x378*x386;
  const double x388 = x383*x386;
  const double x389 = x32*x94;
  const double x390 = x224*x41;
  const double x391 = -x229*x23 + x23*x44;
  const double x392 = x391*x74;
  const double x393 = x159*x389 + x167*x392 - x169*x390;
  const double x394 = x122*x389 + x135*x392 + x140*x390;
  const double x395 = x108*x190*x369*x377*x379*x39*x394 - x388*x393;
  const double x396 = x110*x379;
  const double x397 = x378*x396;
  const double x398 = x383*x396;
  const double x399 = std::pow(x22, -5);
  const double x400 = x108*x399;
  const double x401 = x215*x380;
  const double x402 = x378*x400;
  const double x403 = x136*x251;
  const double x404 = Bx*x233 + By*W0x*x118 + By*x244 - x118*x311 + x151*x19 + x151*x20 - x236*x240 - x245*x246 + x343;
  const double x405 = x11*x371*x404;
  const double x406 = x168*x378*x405;
  const double x407 = x399*x406;
  const double x408 = x383*x405;
  const double x409 = x399*x408;
  const double x410 = x373*x377 + x376*x377;
  const double x411 = x400*x408;
  const double x412 = x402*x405;
  const double x413 = x108*x410*x77*(x100*x389 - x101*x224 + x106*x391) - x393*x412 - x394*x411;
  const double x414 = x109*x410;
  const double x415 = x372*x404;
  const double x416 = x378*x415;
  const double x417 = x383*x415;
  const double x418 = x108*x215/std::pow(x22, 6);
  const double dqopdqop0 = std::pow(q, 3)*x5*x6/(std::pow(qop0, 3)*x3) - x180*x182;
  const double dqopdlam0 = -x182*x202;
  const double dqopdphi0 = -x182*x210;
  const double dqopdxt0 = -x182*x211;
  const double dqopdyt0 = -x182*x214;
  const double dqopdBz = -x182*x248;
  const double dqopdxi = -x181*x4 - x182*x269;
  const double dlamdqop0 = x180*x310 + x276*(x286*((1.0/2.0)*Bx*Bz*W0x*s*x11*x23*x35*x40*x97 + (1.0/2.0)*By*Bz*W0y*s*x11*x23*x35*x40*x97 + (1.0/2.0)*By*W0x*s*x11*x22*x40*x59*x97 - 1.0/2.0*x119*x278 - x280*x283 - x283*x284) + x302*(-x293*(Bx*By*W0x*s*x11*x23*x35*x40*x97 + Bx*W0z*s*x11*x22*x40*x59*x97 + By*Bz*W0z*s*x11*x23*x35*x40*x97 - x287 - x289 - x291) - x301*(x172*x282 + x174*x282 + x294 - x296 - x298 - x300)));
  const double dlamdlam0 = x202*x310 + x276*(x286*(x115 + x116 + x117 + x118*x185*x68 + x148*x183*x245 - x148*x185 + x173*x183*x236 - x173*x184 - x184*x309) + x302*(-x293*(2*Bx*By*W0x*W0z*x39*x40*x59 + 2*Bx*x23*x35*x40 + 2*By*Bz*x40 + 2*Bz*W0x*W0z*x23*x35*x39*x40 - 2*x148*x59 - 2*x319 - 2*x320 - 2*x321 - 2*x322) - x301*(2*Bx*By*W0y*W0z*x39*x40*x59 + 2*Bx*Bz*x40 - 2*x173*x59 - x195*x315 - 2*x312 - 2*x314 - x315*x70 - 2*x316 - 2*x318)));
  const double dlamdphi0 = x210*x310 + x276*(x286*(Bx*Bz*W0y*x11*x40*x59 + By*Bz*W0x*x11*x40 - x118*x69 - x118*x71 - x323 - x325) + x302*(-x293*(x165*x315 - 2*x172 + 2*x175 + 2*x176 + 2*x297 + 2*x328) - x301*(2*Bx*By*W0x*x11*x40 + 2*Bz*W0x*x11*x23*x35*x40 - 2*x150 - 2*x155 - 2*x288 - 2*x327)));
  const double dlamdxt0 = x211*x310;
  const double dlamdyt0 = x214*x310;
  const double dlamdBz = x248*x310 + x276*(-x120*x285*x330 + x286*(-x19*x335 - x20*x335 - x323*x331 - x323*x332*x59 + x324*x331 + x325*x332 + x334*x69 + x334*x71 - x59*x69 - x59*x71 + 2*x67 + x72) + x302*(x271*x336 + x272*x336 - x293*(Bx*By*Bz*W0x*s*x11*x227*x32*x35*x40 + Bx*Bz*W0z*s*x11*x32*x40*x59 + 2*Bx*Bz*W0z*x11*x227*x35*x40 + By*W0z*s*x11*x21*x227*x32*x35*x40 + 2*By*W0z*x11*x40 + 4*Bz*W0y*x11*x40*x59 - x176*x222 - x193*x315 - x199*x342 - 2*x295 - x299*x338 - x339*x346) - x301*(2*x130 + 4*x131*x59 - x149*x338 - x153*x222 + x155*x222 - x193*x342 + x199*x315 + x290*x338 - 2*x337 + x339*x340 + x339*x344*x68 - x339*x345)));
  const double dlamdxi = x269*x310 + x276*(x286*((1.0/2.0)*Bx*Bz*W0x*q*s*x11*x14*x23*x256*x35*x40*x95 + (1.0/2.0)*By*Bz*W0y*q*s*x11*x14*x23*x256*x35*x40*x95 + (1.0/2.0)*By*W0x*q*s*x11*x14*x22*x256*x40*x59*x95 - 1.0/2.0*x119*x348 - x280*x350 - x284*x350) + x302*(-x293*(Bx*By*W0x*q*s*x11*x14*x23*x256*x35*x40*x95 + Bx*W0z*q*s*x11*x14*x22*x256*x40*x59*x95 + By*Bz*W0z*q*s*x11*x14*x23*x256*x35*x40*x95 - x351 - x352 - x353) - x301*(x172*x349 + x174*x349 + x354 - x355 - x356 - x357)));
  const double dphidqop0 = x180*x360 + x358*((1.0/2.0)*Bx*By*W0x*s*x11*x23*x35*x40*x97 + (1.0/2.0)*Bx*W0z*s*x11*x22*x40*x59*x97 + (1.0/2.0)*By*Bz*W0z*s*x11*x23*x35*x40*x97 - 1.0/2.0*x287 - 1.0/2.0*x289 - 1.0/2.0*x291) - x359*(x172*x283 + x174*x283 + (1.0/2.0)*x294 - 1.0/2.0*x296 - 1.0/2.0*x298 - 1.0/2.0*x300);
  const double dphidlam0 = x202*x360 + x358*(Bx*By*W0x*W0z*x39*x40*x59 + Bx*x23*x35*x40 + By*Bz*x40 + Bz*W0x*W0z*x23*x35*x39*x40 - x148*x59 - x319 - x320 - x321 - x322) - x359*(Bx*By*W0y*W0z*x39*x40*x59 + Bx*Bz*x40 - x173*x59 - x309 - x312 - x314 - x316 - x318 - x361*x362*x66);
  const double dphidphi0 = x210*x360 + x358*(-x172 + x175 + x177 + x297 + x328) - x359*(Bx*By*W0x*x11*x40 - x150 - x156 - x288 - x327);
  const double dphidxt0 = x211*x360;
  const double dphidyt0 = x214*x360;
  const double dphidBz = -x158*x274*(x110*(x118*x199 + x130 + 2*x131*x59 - x149*x331 - x153*x332 + x155*x332 - x193*x363 + x290*x331 + x333*x340 + x333*x344*x68 - x333*x345 - x337) - x178*x330) + x179*x274*(x110*((1.0/2.0)*Bx*By*Bz*W0x*s*x11*x227*x32*x35*x40 + (1.0/2.0)*Bx*Bz*W0z*s*x11*x32*x40*x59 + Bx*Bz*W0z*x11*x227*x35*x40 + (1.0/2.0)*By*W0z*s*x11*x21*x227*x32*x35*x40 + By*W0z*x11*x40 + 2*Bz*W0y*x11*x40*x59 - x118*x193 - x176*x332 - x199*x363 - x295 - x299*x331 - x333*x346) - x157*x330) + x248*x360;
  const double dphidxi = x269*x360 + x358*((1.0/2.0)*Bx*By*W0x*q*s*x11*x14*x23*x256*x35*x40*x95 + (1.0/2.0)*Bx*W0z*q*s*x11*x14*x22*x256*x40*x59*x95 + (1.0/2.0)*By*Bz*W0z*q*s*x11*x14*x23*x256*x35*x40*x95 - 1.0/2.0*x351 - 1.0/2.0*x352 - 1.0/2.0*x353) - x359*(x172*x350 + x174*x350 + (1.0/2.0)*x354 - 1.0/2.0*x355 - 1.0/2.0*x356 - 1.0/2.0*x357);
  const double dxtdqop0 = x143*x387 - x170*x388 + x180*x395 + x364*x382 - x385*x88;
  const double dxtdlam0 = -x196*x388 + x200*x387 + x202*x395;
  const double dxtdphi0 = -x208*x388 + x209*x387 + x210*x395;
  const double dxtdxt0 = W0x*x39*x397 + x211*x395 + x361*x398;
  const double dxtdyt0 = x184*x398 - x185*x397 + x214*x395;
  const double dxtdBz = -x136*x401*x402 - x243*x388 + x247*x387 + x248*x395 + x384*x400*x401;
  const double dxtdxi = x251*x385 + x267*x387 - x268*x388 + x269*x395 - x382*x403;
  const double dytdqop0 = x107*x108*x410*x77 - x143*x411 - x170*x412 + x180*x413 - x364*x409 - x407*x88 + x410*x76*x77*x78*x87;
  const double dytdlam0 = x189*x414 - x196*x412 - x200*x411 + x202*x413;
  const double dytdphi0 = x205*x414 - x208*x412 - x209*x411 + x210*x413;
  const double dytdxt0 = x211*x413 + x416*x52 - x417*x54;
  const double dytdyt0 = x10*x410 + x212*x416 + x213*x417 + x214*x413;
  const double dytdBz = x109*x231*x410 + x136*x408*x418 - x216*x410 - x243*x412 - x247*x411 + x248*x413 + x406*x418;
  const double dytdxi = x109*x262*x410 + x251*x407 - x253*x410 - x267*x411 - x268*x412 + x269*x413 + x403*x409;
  Eigen::Matrix<double, 5, 7> res;
  res(0,0) = dqopdqop0;
  res(0,1) = dqopdlam0;
  res(0,2) = dqopdphi0;
  res(0,3) = dqopdxt0;
  res(0,4) = dqopdyt0;
  res(0,5) = dqopdBz;
  res(0,6) = dqopdxi;
  res(1,0) = dlamdqop0;
  res(1,1) = dlamdlam0;
  res(1,2) = dlamdphi0;
  res(1,3) = dlamdxt0;
  res(1,4) = dlamdyt0;
  res(1,5) = dlamdBz;
  res(1,6) = dlamdxi;
  res(2,0) = dphidqop0;
  res(2,1) = dphidlam0;
  res(2,2) = dphidphi0;
  res(2,3) = dphidxt0;
  res(2,4) = dphidyt0;
  res(2,5) = dphidBz;
  res(2,6) = dphidxi;
  res(3,0) = dxtdqop0;
  res(3,1) = dxtdlam0;
  res(3,2) = dxtdphi0;
  res(3,3) = dxtdxt0;
  res(3,4) = dxtdyt0;
  res(3,5) = dxtdBz;
  res(3,6) = dxtdxi;
  res(4,0) = dytdqop0;
  res(4,1) = dytdlam0;
  res(4,2) = dytdphi0;
  res(4,3) = dytdxt0;
  res(4,4) = dytdyt0;
  res(4,5) = dytdBz;
  res(4,6) = dytdxi;



  res.col(5) *= 2.99792458e-3;

  return res;


}
