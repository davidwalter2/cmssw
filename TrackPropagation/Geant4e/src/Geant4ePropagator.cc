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

  bool continuePropagation = true;
  while (continuePropagation) {
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
    
    //set the error matrix to null to disentangle MS and ionization contributions
    g4eTrajState.SetError(errnull);
    const int ierr = theG4eManager->PropagateOneStep(&g4eTrajState, mode);

    if (ierr != 0) {
      // propagation failed, return invalid track state
      return retDefault();
    }

    const double thisPathLength = TrackPropagation::g4doubleToCmsDouble(g4eTrajState.GetG4Track()->GetStepLength());

    const double pPre = g4eTrajState.GetMomentum().mag()/GeV;
    const double ePre = g4eTrajState.GetG4Track()->GetStep()->GetPreStepPoint()->GetTotalEnergy() / GeV;
    const double ePost = g4eTrajState.GetG4Track()->GetStep()->GetPostStepPoint()->GetTotalEnergy() / GeV;
    const double dEdx = (ePost - ePre)/thisPathLength;
    const double mass  = g4eTrajState.GetG4Track()->GetDynamicParticle()->GetMass() / GeV;

    if (std::abs(thisPathLength) > 0.) {
      dEdxlast = dEdx;
//       masslast = mass;
    }
    
//     std::cout << "thisPathLength = " << thisPathLength << std::endl;
//     const Matrix<double, 5, 7> transportJac = transportJacobian(statepre, thisPathLength, dEdx, mass);
    const Matrix<double, 5, 7> transportJac = transportJacobianBzD(statepre, thisPathLength, dEdx, mass, dBz);
    
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

//   std::cout << "propgationDirection = " << propagationDirection() << " finalPathLength = " << finalPathLength << std::endl;
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

Eigen::Matrix<double, 6, 1> Geant4ePropagator::transportResultD(const Eigen::Matrix<double, 7, 1> &start, double s, double dEdx, double mass, double bzoffset) const {
  
//   const GlobalTrajectoryParameters &globalSource = start.parameters();
  
  if (s==0.) {
    return start.head<6>();
  }
  
  const GlobalPoint pos(start[0], start[1], start[2]);  
  const GlobalVector &bfield = theField->inInverseGeV(pos);

  const double M0x = start[0];
  const double M0y = start[1];
  const double M0z = start[2];
  
  const Eigen::Matrix<double, 3, 1> W0 = start.segment<3>(3).normalized();
  const double W0x = W0[0];
  const double W0y = W0[1];
  const double W0z = W0[2];
    
  const double Bx = bfield.x();
  const double By = bfield.y();
  const double Bz = bfield.z() + 2.99792458e-3*bzoffset;
  
  const double q = start[6];
  const double qop0 = q/start.segment<3>(3).norm();

  const double xc0 = std::pow(Bx, 2) + std::pow(By, 2) + std::pow(Bz, 2);
  const double xc1 = qop0*s*std::sqrt(xc0);
  const double xc2 = std::sin(xc1);
  const double xc3 = std::pow(xc0, 5.0/2.0)*xc2;
  const double xc4 = By*W0z - Bz*W0y;
  const double xc5 = std::cos(xc1);
  const double xc6 = std::pow(xc0, 2)*(xc5 - 1);
  const double xc7 = Bx*W0x + By*W0y + Bz*W0z;
  const double xc8 = Bx*xc7;
  const double xc9 = std::pow(xc0, 3.0/2.0)*(xc1 - xc2);
  const double xc10 = std::pow(W0x, 2) + std::pow(W0y, 2);
  const double xc11 = std::pow(W0z, 2) + xc10;
  const double xc12 = xc11/xc10;
  const double xc13 = std::pow(xc0, 3);
  const double xc14 = 1.0/xc13;
  const double xc15 = xc11*xc14/(qop0*std::pow(xc10, 3.0/2.0)*std::pow(xc12, 3.0/2.0));
  const double xc16 = Bx*W0z - Bz*W0x;
  const double xc17 = xc7*xc9;
  const double xc18 = Bx*W0y - By*W0x;
  const double xc19 = xc13*xc5;
  const double xc20 = std::fabs(qop0);
  const double xc21 = std::pow(q, 2);
  const double xc22 = std::pow(qop0, 2);
  const double xc23 = xc14*std::sqrt(std::pow(dEdx, 2)*std::pow(s, 2)*xc20*xc22 + 2*dEdx*s*xc22*std::sqrt(std::pow(mass, 2)*xc22 + xc21) + xc20*xc21)/(std::sqrt(xc10)*std::sqrt(xc12)*std::pow(xc20, 3.0/2.0));
  const double xc24 = xc6*xc7;
  const double x = xc15*(W0x*xc3 + xc4*xc6 + xc8*xc9);
  const double y = xc15*(By*xc17 + W0y*xc3 - xc16*xc6);
  const double z = xc15*(Bz*xc17 + W0z*xc3 + xc18*xc6);
  const double px = xc23*(W0x*xc19 - xc3*xc4 - xc6*xc8);
  const double py = xc23*(-By*xc24 + W0y*xc19 + xc16*xc3);
  const double pz = xc23*(-Bz*xc24 + W0z*xc19 - xc18*xc3);
  Eigen::Matrix<double, 6, 1> rescart;
  rescart[0] = x;
  rescart[1] = y;
  rescart[2] = z;
  rescart[3] = px;
  rescart[4] = py;
  rescart[5] = pz;



  return rescart;
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
                
  const double x0 = std::pow(q, 2);
  const double x1 = x0/std::pow(qop0, 3);
  const double x2 = std::pow(qop0, -2);
  const double x3 = x0*x2;
  const double x4 = std::sqrt(std::pow(mass, 2) + x3);
  const double x5 = std::pow(dEdx, 2);
  const double x6 = std::pow(s, 2)*x5;
  const double x7 = dEdx*x4;
  const double x8 = s*x7;
  const double x9 = q/std::pow(x3 + x6 + 2*x8, 3.0/2.0);
  const double x10 = std::pow(W0z, 2);
  const double x11 = std::pow(W0x, 2);
  const double x12 = std::pow(W0y, 2);
  const double x13 = x11 + x12;
  const double x14 = 1.0/x13;
  const double x15 = std::pow(x10*x14 + 1, -1.0/2.0);
  const double x16 = std::pow(Bz, 2);
  const double x17 = std::pow(Bx, 2) + std::pow(By, 2) + x16;
  const double x18 = std::sqrt(x17);
  const double x19 = s*x18;
  const double x20 = qop0*x19;
  const double x21 = std::sin(x20);
  const double x22 = std::pow(x17, 5.0/2.0)*x21;
  const double x23 = x15*x22;
  const double x24 = std::sqrt(x13);
  const double x25 = 1.0/x24;
  const double x26 = W0z*x25;
  const double x27 = Bx*W0y;
  const double x28 = By*W0x;
  const double x29 = x25*x27 - x25*x28;
  const double x30 = std::pow(x17, 2);
  const double x31 = std::cos(x20);
  const double x32 = x31 - 1;
  const double x33 = x30*x32;
  const double x34 = x15*x33;
  const double x35 = x29*x34;
  const double x36 = std::pow(x17, 3);
  const double x37 = W0x*x25;
  const double x38 = W0z*x37;
  const double x39 = W0y*x25;
  const double x40 = W0z*x39;
  const double x41 = -M0x*x38 - M0y*x40 + M0z*(x11*x25 + x12*x25);
  const double x42 = x24*x41;
  const double x43 = x36*x42;
  const double x44 = x20 - x21;
  const double x45 = std::pow(x17, 3.0/2.0);
  const double x46 = Bx*x15;
  const double x47 = By*x15;
  const double x48 = Bz*x15;
  const double x49 = x26*x48;
  const double x50 = x37*x46 + x39*x47 + x49;
  const double x51 = x45*x50;
  const double x52 = x44*x51;
  const double x53 = Bz*x52 + qop0*x43 + x23*x26 + x35;
  const double x54 = 1.0/x36;
  const double x55 = x2*x54;
  const double x56 = x53*x55;
  const double x57 = x31*x36;
  const double x58 = x15*x57;
  const double x59 = x26*x58;
  const double x60 = x23*x29;
  const double x61 = -x19*x31 + x19;
  const double x62 = Bz*x51;
  const double x63 = 1.0/qop0;
  const double x64 = x54*x63;
  const double x65 = x64*(s*x59 - s*x60 + x43 + x61*x62);
  const double x66 = x33*x50;
  const double x67 = -Bz*x66 + x59 - x60;
  const double x68 = x54*x67;
  const double x69 = W0y*x23;
  const double x70 = x26*x46;
  const double x71 = x37*x48;
  const double x72 = x70 - x71;
  const double x73 = x24*x33;
  const double x74 = -M0x*x39 + M0y*x37;
  const double x75 = W0z*x41;
  const double x76 = W0x*x74 - W0y*x75;
  const double x77 = x36*x76;
  const double x78 = By*x24;
  const double x79 = qop0*x77 + x52*x78 + x69 - x72*x73;
  const double x80 = x25*x55;
  const double x81 = s*x58;
  const double x82 = x22*x72;
  const double x83 = s*x24;
  const double x84 = x51*x61;
  const double x85 = W0y*x81 + x77 + x78*x84 + x82*x83;
  const double x86 = x25*x64;
  const double x87 = x39*x58;
  const double x88 = -By*x66 + x82 + x87;
  const double x89 = x54*x88;
  const double x90 = W0x*x23;
  const double x91 = x26*x47;
  const double x92 = x39*x48;
  const double x93 = x91 - x92;
  const double x94 = W0x*x75 + W0y*x74;
  const double x95 = x36*x94;
  const double x96 = Bx*x24;
  const double x97 = -qop0*x95 + x52*x96 + x73*x93 + x90;
  const double x98 = x22*x93;
  const double x99 = W0x*x81 - x83*x98 + x84*x96 - x95;
  const double x100 = x37*x58;
  const double x101 = -Bx*x66 + x100 - x98;
  const double x102 = x101*x54;
  const double x103 = -x102*(-x80*x97 + x86*x99) - x68*(-x56 + x65) - x89*(-x79*x80 + x85*x86);
  const double x104 = x9*(-s*x5 - x7);
  const double x105 = W0z*x14;
  const double x106 = W0x*x105;
  const double x107 = W0y*x105;
  const double x108 = -x106*x46 - x107*x47 + x48;
  const double x109 = Bz*x108;
  const double x110 = x44*x45;
  const double x111 = x109*x110 + x23 - x26*x35;
  const double x112 = std::pow(x17, -6);
  const double x113 = x112*x63;
  const double x114 = x113*x67;
  const double x115 = x107*x48 + x47;
  const double x116 = x108*x110;
  const double x117 = x115*x73 + x116*x96 - x23*x38;
  const double x118 = x113*x25;
  const double x119 = x101*x118;
  const double x120 = x106*x48 + x46;
  const double x121 = x116*x78 - x120*x73 - x23*x40;
  const double x122 = x118*x88;
  const double x123 = -x111*x114 - x117*x119 - x121*x122;
  const double x124 = Bx*W0x;
  const double x125 = By*W0y;
  const double x126 = x124*x25 + x125*x25;
  const double x127 = x37*x47 - x39*x46;
  const double x128 = Bz*x127;
  const double x129 = x110*x128 + x126*x34;
  const double x130 = x33*x48;
  const double x131 = x110*x127;
  const double x132 = -W0y*x130 + x131*x78 + x90;
  const double x133 = -W0x*x130 + x131*x96 - x69;
  const double x134 = -x114*x129 - x119*x133 - x122*x132;
  const double x135 = x102*x39 - x37*x89;
  const double x136 = x102*x38 - x24*x68 + x40*x89;
  const double x137 = 6*Bz;
  const double x138 = x137/std::pow(x17, 4);
  const double x139 = x138*x63;
  const double x140 = x139*x53;
  const double x141 = x21*x45;
  const double x142 = 5*x141;
  const double x143 = x30*x31;
  const double x144 = qop0*s;
  const double x145 = x143*x144;
  const double x146 = x17*x32;
  const double x147 = 4*x146;
  const double x148 = x29*x48;
  const double x149 = x144*x148;
  const double x150 = Bz*qop0;
  const double x151 = 6*x30;
  const double x152 = x150*x151;
  const double x153 = x16*x50;
  const double x154 = 3*x18*x44;
  const double x155 = s*x150;
  const double x156 = x155/x18;
  const double x157 = -x156*x31 + x156;
  const double x158 = x64*(x110*x49 - x141*x149 + x142*x49 + x145*x49 + x147*x148 + x152*x42 + x153*x154 + x157*x62 + x52);
  const double x159 = x139*x25;
  const double x160 = x142*x48;
  const double x161 = x145*x48;
  const double x162 = W0z*x110;
  const double x163 = Bz*x72;
  const double x164 = x147*x24;
  const double x165 = x155*x72;
  const double x166 = x141*x24;
  const double x167 = Bz*x50;
  const double x168 = x154*x167;
  const double x169 = x157*x51;
  const double x170 = W0x*x34 + W0y*x160 + W0y*x161 + x152*x76 + x162*x47 - x163*x164 + x165*x166 + x168*x78 + x169*x78;
  const double x171 = Bz*x93;
  const double x172 = x155*x93;
  const double x173 = W0x*x160 + W0x*x161 - W0y*x34 - x152*x94 + x162*x46 + x164*x171 - x166*x172 + x168*x96 + x169*x96;
  const double x174 = -x102*(-x159*x97 + x173*x86) - x68*(-x140 + x158) - x89*(-x159*x79 + x170*x86);
  const double x175 = std::pow(x88, 2);
  const double x176 = std::pow(x101, 2);
  const double x177 = x112*x175 + x112*x176;
  const double x178 = 1.0/x177;
  const double x179 = x112*x178;
  const double x180 = 1.0/(x179*std::pow(x67, 2) + 1);
  const double x181 = s*x21;
  const double x182 = x15*std::pow(x17, 7.0/2.0);
  const double x183 = x182*x26;
  const double x184 = x22*x50;
  const double x185 = s*x184;
  const double x186 = std::pow(x177, -1.0/2.0);
  const double x187 = x186*x54;
  const double x188 = x181*x182;
  const double x189 = x188*x39;
  const double x190 = s*x57;
  const double x191 = x190*x72;
  const double x192 = By*x185;
  const double x193 = (1.0/2.0)*x112;
  const double x194 = x193*x88;
  const double x195 = x188*x37;
  const double x196 = x190*x93;
  const double x197 = Bx*x185;
  const double x198 = x101*x193;
  const double x199 = x68/std::pow(x177, 3.0/2.0);
  const double x200 = qop0*x21;
  const double x201 = qop0*x58;
  const double x202 = x182*x200;
  const double x203 = x202*x39;
  const double x204 = qop0*x57;
  const double x205 = x204*x72;
  const double x206 = qop0*x184;
  const double x207 = By*x206;
  const double x208 = x202*x37;
  const double x209 = x204*x93;
  const double x210 = Bx*x206;
  const double x211 = x180*(x187*(x150*x184 - x183*x200 - x201*x29) + x199*(-x194*(-2*x203 + 2*x205 + 2*x207) - x198*(-2*x208 - 2*x209 + 2*x210)));
  const double x212 = x107*x58;
  const double x213 = x120*x22;
  const double x214 = x108*x33;
  const double x215 = By*x214;
  const double x216 = x106*x58;
  const double x217 = x115*x22;
  const double x218 = Bx*x214;
  const double x219 = x22*x92;
  const double x220 = 2*x219;
  const double x221 = x127*x33;
  const double x222 = By*x221;
  const double x223 = x22*x71;
  const double x224 = 2*x223;
  const double x225 = Bx*x221;
  const double x226 = x151*x31;
  const double x227 = qop0*x181;
  const double x228 = x137/std::pow(x17, 7);
  const double x229 = x23*x37;
  const double x230 = 12*x143;
  const double x231 = x33*x91;
  const double x232 = 10*x141;
  const double x233 = x143*x165;
  const double x234 = By*x167;
  const double x235 = 8*x146;
  const double x236 = x227*x62;
  const double x237 = By*x236;
  const double x238 = x23*x39;
  const double x239 = x33*x70;
  const double x240 = x143*x172;
  const double x241 = Bx*x167;
  const double x242 = Bx*x236;
  const double x243 = x101*x179;
  const double x244 = x179*x88;
  const double x245 = x243*(-x203 + x205 + x207) - x244*(-x208 - x209 + x210);
  const double x246 = x33*(Bz*W0z + x124 + x125);
  const double x247 = -By*x246 + W0y*x57 + x22*(Bx*W0z - Bz*W0x);
  const double x248 = x10 + x13;
  const double x249 = 1.0/x248;
  const double x250 = x112*x249;
  const double x251 = std::pow(x247, 2)*x250;
  const double x252 = -Bx*x246 + W0x*x57 - x22*(By*W0z - Bz*W0y);
  const double x253 = x250*std::pow(x252, 2);
  const double x254 = std::pow(x251 + x253, -1.0/2.0);
  const double x255 = x247*x254;
  const double x256 = x2*x97;
  const double x257 = std::pow(x14*x248, -1.0/2.0);
  const double x258 = x14*x257;
  const double x259 = x112*x258;
  const double x260 = x252*x254;
  const double x261 = x2*x79;
  const double x262 = x113*x258;
  const double x263 = x260*x262;
  const double x264 = x255*x262;
  const double x265 = qop0*x24;
  const double x266 = qop0*x18;
  const double x267 = -x266*x31 + x266;
  const double x268 = x267*x51;
  const double x269 = W0y*x201 + x265*x82 + x268*x78;
  const double x270 = W0x*x201 - x265*x98 + x268*x96;
  const double x271 = x263*x269 - x264*x270;
  const double x272 = x258*x54;
  const double x273 = x257*x54;
  const double x274 = x228*x258;
  const double x275 = x255*x63;
  const double x276 = x260*x63;
  const double x277 = -Bz*x246 + W0z*x57 - x22*(x27 - x28);
  const double x278 = x249*x25*x277/std::pow(x17, 9);
  const double x279 = x255*x278;
  const double x280 = x279*x63;
  const double x281 = x276*x278;
  const double x282 = x251*x254 + x253*x254;
  const double x283 = x282*x64;
  const double x284 = -x269*x280 - x270*x281 + x283*(qop0*x59 - qop0*x60 + x267*x62);
  const double x285 = x250*x277;
  const double x286 = x255*x285;
  const double x287 = x260*x285;
  const double x288 = x137*x249*x25*x277/std::pow(x17, 10);
  const double dqopdqop0 = x103*x104 + x9*(dEdx*s*x1/x4 + x1);
  const double dqopdlam0 = x104*x123;
  const double dqopdphi0 = x104*x134;
  const double dqopdxt0 = x104*x135;
  const double dqopdyt0 = x104*x136;
  const double dqopdBz = x104*x174;
  const double dqopdxi = x9*(-x6 - x8);
  const double dlamdqop0 = x103*x211 + x180*(x187*(Bz*x185 - x181*x183 - x29*x81) + x199*(-x194*(-2*x189 + 2*x191 + 2*x192) - x198*(-2*x195 - 2*x196 + 2*x197)));
  const double dlamdlam0 = x123*x211 + x180*(x187*(-x109*x33 + x26*x60 + x58) + x199*(-x194*(-2*x212 + 2*x213 - 2*x215) - x198*(-2*x216 - 2*x217 - 2*x218)));
  const double dlamdphi0 = x134*x211 + x180*(x187*(-x126*x23 - x128*x33) + x199*(-x194*(2*x100 + x220 - 2*x222) - x198*(x224 - 2*x225 - 2*x87)));
  const double dlamdxt0 = x135*x211;
  const double dlamdyt0 = x136*x211;
  const double dlamdBz = x174*x211 + x180*(-x138*x186*x67 + x187*(-x142*x148 - x143*x149 - x144*x22*x49 - x147*x153 + x16*x227*x51 + x226*x49 - x33*x49 - x66) + x199*(x175*x228 + x176*x228 - x194*(-x144*x220 + x163*x232 - 2*x229 + x230*x92 - 2*x231 + 2*x233 - x234*x235 + 2*x237) - x198*(-x144*x224 - x171*x232 + x230*x71 - x235*x241 + 2*x238 - 2*x239 - 2*x240 + 2*x242)));
  const double dlamdxi = 0;
  const double dphidqop0 = x103*x245 + x243*(-x189 + x191 + x192) - x244*(-x195 - x196 + x197);
  const double dphidlam0 = x123*x245 + x243*(-x212 + x213 - x215) - x244*(-x216 - x217 - x218);
  const double dphidphi0 = x134*x245 + x243*(x100 + x219 - x222) - x244*(x223 - x225 - x87);
  const double dphidxt0 = x135*x245;
  const double dphidyt0 = x136*x245;
  const double dphidBz = x102*x178*(-x138*x88 + x54*(x142*x163 - x144*x219 - x147*x234 + x226*x92 - x229 - x231 + x233 + x237)) + x174*x245 - x178*x89*(-x101*x138 + x54*(-x142*x171 - x144*x223 - x147*x241 + x226*x71 + x238 - x239 - x240 + x242));
  const double dphidxi = 0;
  const double dxtdqop0 = x103*x271 + x255*x256*x259 - x259*x260*x261 + x263*x85 - x264*x99;
  const double dxtdlam0 = -x117*x264 + x121*x263 + x123*x271;
  const double dxtdphi0 = x132*x263 - x133*x264 + x134*x271;
  const double dxtdxt0 = W0x*x260*x272 + W0y*x255*x272 + x135*x271;
  const double dxtdyt0 = x106*x255*x273 - x107*x260*x273 + x136*x271;
  const double dxtdBz = x170*x263 - x173*x264 + x174*x271 + x274*x275*x97 - x274*x276*x79;
  const double dxtdxi = 0;
  const double dytdqop0 = x103*x284 + x256*x260*x278 + x261*x279 - x280*x85 - x281*x99 - x282*x56 + x282*x65;
  const double dytdlam0 = x111*x283 - x117*x281 - x121*x280 + x123*x284;
  const double dytdphi0 = x129*x283 - x132*x280 - x133*x281 + x134*x284;
  const double dytdxt0 = x135*x284 - x286*x37 + x287*x39;
  const double dytdyt0 = x136*x284 + x24*x282 + x286*x40 + x287*x38;
  const double dytdBz = -x140*x282 + x158*x282 - x170*x280 - x173*x281 + x174*x284 + x275*x288*x79 + x276*x288*x97;
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

// Eigen::Matrix<double, 5, 7> Geant4ePropagator::transportJacobianD(const Eigen::Matrix<double, 7, 1> &start, double s, double dEdx, double mass) const {
//
//   if (s==0.) {
//     Eigen::Matrix<double, 5, 7> res;
//     res.leftCols<5>() = Eigen::Matrix<double, 5, 5>::Identity();
//     res.rightCols<2>() = Eigen::Matrix<double, 5, 2>::Zero();
//     return res;
//   }
//
//   const GlobalPoint pos(start[0], start[1], start[2]);
//   const GlobalVector &bfield = theField->inInverseGeV(pos);
//
//   const Eigen::Matrix<double, 3, 1> Bv(bfield.x(), bfield.y(), bfield.z());
//   const double B = Bv.norm();
//   const Eigen::Matrix<double, 3, 1> Hv = Bv.normalized();
//   const double hx = Hv[0];
//   const double hy = Hv[1];
//   const double hz = Hv[2];
//
//   const double M0x = start[0];
//   const double M0y = start[1];
//   const double M0z = start[2];
//
//   const Eigen::Matrix<double, 3, 1> W0 = start.segment<3>(3).normalized();
//   const double W0x = W0[0];
//   const double W0y = W0[1];
//   const double W0z = W0[2];
//
//   const double q = start[6];
//   const double qop0 = q/start.segment<3>(3).norm();
//
//   const double x0 = std::pow(q, 2);
//   const double x1 = x0/std::pow(qop0, 3);
//   const double x2 = std::pow(qop0, -2);
//   const double x3 = x0*x2;
//   const double x4 = std::sqrt(std::pow(mass, 2) + x3);
//   const double x5 = std::pow(dEdx, 2);
//   const double x6 = std::pow(s, 2)*x5;
//   const double x7 = dEdx*x4;
//   const double x8 = s*x7;
//   const double x9 = q/std::pow(x3 + x6 + 2*x8, 3.0/2.0);
//   const double x10 = B*s;
//   const double x11 = qop0*x10;
//   const double x12 = std::cos(x11);
//   const double x13 = std::pow(W0z, 2);
//   const double x14 = std::pow(W0x, 2);
//   const double x15 = std::pow(W0y, 2);
//   const double x16 = x14 + x15;
//   const double x17 = 1.0/x16;
//   const double x18 = std::pow(x13*x17 + 1, -1.0/2.0);
//   const double x19 = x12*x18;
//   const double x20 = std::sqrt(x16);
//   const double x21 = 1.0/x20;
//   const double x22 = W0z*x21;
//   const double x23 = x19*x22;
//   const double x24 = W0y*hx;
//   const double x25 = x21*x24;
//   const double x26 = W0x*hy;
//   const double x27 = x21*x26;
//   const double x28 = x25 - x27;
//   const double x29 = std::sin(x11);
//   const double x30 = x18*x29;
//   const double x31 = x28*x30;
//   const double x32 = x12 - 1;
//   const double x33 = hx*x18;
//   const double x34 = W0x*x21;
//   const double x35 = hy*x18;
//   const double x36 = W0y*x21;
//   const double x37 = hz*x18;
//   const double x38 = x22*x37 + x33*x34 + x35*x36;
//   const double x39 = hz*x38;
//   const double x40 = x23 - x31 - x32*x39;
//   const double x41 = x22*x30;
//   const double x42 = x18*x32;
//   const double x43 = x28*x42;
//   const double x44 = B*qop0;
//   const double x45 = W0z*x34;
//   const double x46 = W0z*x36;
//   const double x47 = -M0x*x45 - M0y*x46 + M0z*(x14*x21 + x15*x21);
//   const double x48 = x20*x47;
//   const double x49 = x11 - x29;
//   const double x50 = x39*x49 + x41 + x43 + x44*x48;
//   const double x51 = 1.0/B;
//   const double x52 = x2*x51;
//   const double x53 = x50*x52;
//   const double x54 = x10*x12;
//   const double x55 = x10 - x54;
//   const double x56 = 1.0/qop0;
//   const double x57 = x51*x56;
//   const double x58 = x57*(B*x48 + x10*x23 - x10*x31 + x39*x55);
//   const double x59 = x19*x36;
//   const double x60 = x34*x37;
//   const double x61 = x22*x33 - x60;
//   const double x62 = x29*x61;
//   const double x63 = hy*x38;
//   const double x64 = x32*x63;
//   const double x65 = x59 + x62 - x64;
//   const double x66 = -M0x*x36 + M0y*x34;
//   const double x67 = W0z*x47;
//   const double x68 = W0x*x66 - W0y*x67;
//   const double x69 = B*x68;
//   const double x70 = x30*x36;
//   const double x71 = qop0*x69 + x20*(-x32*x61 + x49*x63 + x70);
//   const double x72 = x21*x52;
//   const double x73 = x71*x72;
//   const double x74 = x20*(x10*x59 + x10*x62 + x55*x63) + x69;
//   const double x75 = x21*x57;
//   const double x76 = x74*x75;
//   const double x77 = x19*x34;
//   const double x78 = x22*x35;
//   const double x79 = x36*x37;
//   const double x80 = x78 - x79;
//   const double x81 = x29*x80;
//   const double x82 = hx*x38;
//   const double x83 = -x32*x82 + x77 - x81;
//   const double x84 = W0x*x67 + W0y*x66;
//   const double x85 = B*x84;
//   const double x86 = x30*x34;
//   const double x87 = -qop0*x85 + x20*(x32*x80 + x49*x82 + x86);
//   const double x88 = x72*x87;
//   const double x89 = x20*(x10*x77 - x10*x81 + x55*x82) - x85;
//   const double x90 = x75*x89;
//   const double x91 = -x40*(-x53 + x58) - x65*(-x73 + x76) - x83*(-x88 + x90);
//   const double x92 = x9*(-s*x5 - x7);
//   const double x93 = W0z*x17;
//   const double x94 = W0x*x93;
//   const double x95 = W0y*x93;
//   const double x96 = -x33*x94 - x35*x95 + x37;
//   const double x97 = x49*x96;
//   const double x98 = hz*x97 - x22*x43 + x30;
//   const double x99 = x40*x57;
//   const double x100 = x37*x95;
//   const double x101 = hx*x97 - x30*x94 + x32*(x100 + x35);
//   const double x102 = x57*x83;
//   const double x103 = 1 - x12;
//   const double x104 = x33 + x37*x94;
//   const double x105 = hy*x97 + x103*x104 - x30*x95;
//   const double x106 = x57*x65;
//   const double x107 = -x101*x102 - x105*x106 - x98*x99;
//   const double x108 = W0x*hx;
//   const double x109 = x108*x21;
//   const double x110 = W0y*hy;
//   const double x111 = x110*x21;
//   const double x112 = -x33*x36 + x34*x35;
//   const double x113 = x112*x49;
//   const double x114 = hz*x113 + x42*(x109 + x111);
//   const double x115 = hy*x113 + x103*x79 + x86;
//   const double x116 = hx*x113 - x32*x60 - x70;
//   const double x117 = -x102*x116 - x106*x115 - x114*x99;
//   const double x118 = -x34*x65 + x36*x83;
//   const double x119 = -x20*x40 + x45*x83 + x46*x65;
//   const double x120 = x56/std::pow(B, 2);
//   const double x121 = x120*x50;
//   const double x122 = qop0*s;
//   const double x123 = x12*x122;
//   const double x124 = x122 - x123;
//   const double x125 = x57*(qop0*x48 + x122*x23 - x122*x31 + x124*x39);
//   const double x126 = x120*x21;
//   const double x127 = x126*x71;
//   const double x128 = qop0*x68 + x20*(x122*x59 + x122*x62 + x124*x63);
//   const double x129 = x128*x75;
//   const double x130 = x126*x87;
//   const double x131 = -qop0*x84 + x20*(x122*x77 - x122*x81 + x124*x82);
//   const double x132 = x131*x75;
//   const double x133 = -x40*(-x121 + x125) - x65*(-x127 + x129) - x83*(-x130 + x132);
//   const double x134 = std::pow(x65, 2) + std::pow(x83, 2);
//   const double x135 = 1.0/x134;
//   const double x136 = 1.0/(x135*std::pow(x40, 2) + 1);
//   const double x137 = -x25 + x27;
//   const double x138 = x137*x19;
//   const double x139 = x29*x39;
//   const double x140 = std::pow(x134, -1.0/2.0);
//   const double x141 = x10*x70;
//   const double x142 = x54*x61;
//   const double x143 = x10*x29;
//   const double x144 = x143*x63;
//   const double x145 = (1.0/2.0)*x65;
//   const double x146 = x10*x86;
//   const double x147 = -x78 + x79;
//   const double x148 = x147*x54;
//   const double x149 = x143*x82;
//   const double x150 = (1.0/2.0)*x83;
//   const double x151 = x40/std::pow(x134, 3.0/2.0);
//   const double x152 = x44*x70;
//   const double x153 = x12*x44;
//   const double x154 = x153*x61;
//   const double x155 = x29*x44;
//   const double x156 = x155*x63;
//   const double x157 = x44*x86;
//   const double x158 = x147*x153;
//   const double x159 = x155*x82;
//   const double x160 = x136*(x140*(x138*x44 + x139*x44 - x41*x44) + x151*(-x145*(-2*x152 + 2*x154 + 2*x156) - x150*(-2*x157 + 2*x158 + 2*x159)));
//   const double x161 = hz*x32;
//   const double x162 = x19*x95;
//   const double x163 = x104*x29;
//   const double x164 = x32*x96;
//   const double x165 = hy*x164;
//   const double x166 = x19*x94;
//   const double x167 = x29*(-x100 - x35);
//   const double x168 = hx*x164;
//   const double x169 = W0y*hz;
//   const double x170 = x21*x30;
//   const double x171 = x169*x170;
//   const double x172 = x112*x32;
//   const double x173 = hy*x172;
//   const double x174 = W0x*hz;
//   const double x175 = x170*x174;
//   const double x176 = hx*x172;
//   const double x177 = x122*x70;
//   const double x178 = x123*x61;
//   const double x179 = x122*x29;
//   const double x180 = x179*x63;
//   const double x181 = x122*x86;
//   const double x182 = x123*x147;
//   const double x183 = x179*x82;
//   const double x184 = x135*x83;
//   const double x185 = -x59;
//   const double x186 = x135*(x185 - x62 + x64);
//   const double x187 = x184*(-x152 + x154 + x156) + x186*(-x157 + x158 + x159);
//   const double x188 = W0z*hz + x108 + x110;
//   const double x189 = x188*x32;
//   const double x190 = W0y*x12 - hy*x189 + x29*(W0z*hx - x174);
//   const double x191 = x13 + x16;
//   const double x192 = 1.0/x191;
//   const double x193 = W0x*x12 - hx*x189 + x29*(-W0z*hy + x169);
//   const double x194 = x192*std::pow(x193, 2);
//   const double x195 = std::pow(x190, 2)*x192;
//   const double x196 = std::pow(x194 + x195, -1.0/2.0);
//   const double x197 = x196/std::sqrt(x17*x191);
//   const double x198 = x17*x197;
//   const double x199 = x190*x198;
//   const double x200 = x199*x87;
//   const double x201 = x193*x198;
//   const double x202 = x201*x71;
//   const double x203 = x201*x57;
//   const double x204 = x199*x57;
//   const double x205 = -x153 + x44;
//   const double x206 = x205*x63 + x44*x59 + x44*x62;
//   const double x207 = x193*x197;
//   const double x208 = x207*x75;
//   const double x209 = x205*x82 + x44*x77 - x44*x81;
//   const double x210 = x190*x197;
//   const double x211 = x210*x75;
//   const double x212 = x206*x208 - x209*x211;
//   const double x213 = x192*x196*(W0z*x12 - x161*x188 + x29*(-x24 + x26));
//   const double x214 = x193*x213;
//   const double x215 = x190*x213;
//   const double x216 = x194*x196 + x195*x196;
//   const double x217 = x215*x57;
//   const double x218 = x214*x57;
//   const double x219 = x216*x57;
//   const double x220 = -x206*x217 - x209*x218 + x219*(x205*x39 + x23*x44 - x31*x44);
//   const double dqopdqop0 = x9*(dEdx*s*x1/x4 + x1) + x91*x92;
//   const double dqopdlam0 = x107*x92;
//   const double dqopdphi0 = x117*x92;
//   const double dqopdxt0 = x118*x92;
//   const double dqopdyt0 = x119*x92;
//   const double dqopdB = x133*x92;
//   const double dqopdxi = x9*(-x6 - x8);
//   const double dlamdqop0 = x136*(x140*(x10*x138 + x10*x139 - x10*x41) + x151*(-x145*(-2*x141 + 2*x142 + 2*x144) - x150*(-2*x146 + 2*x148 + 2*x149))) + x160*x91;
//   const double dlamdlam0 = x107*x160 + x136*(x140*(-x137*x41 - x161*x96 + x19) + x151*(-x145*(-2*x162 + 2*x163 - 2*x165) - x150*(-2*x166 + 2*x167 - 2*x168)));
//   const double dlamdphi0 = x117*x160 + x136*(x140*(-x112*x161 + x30*(-x109 - x111)) + x151*(-x145*(2*x171 - 2*x173 + 2*x77) - x150*(2*x175 - 2*x176 - 2*x59)));
//   const double dlamdxt0 = x118*x160;
//   const double dlamdyt0 = x119*x160;
//   const double dlamdB = x133*x160 + x136*(x140*(x122*x138 + x122*x139 - x122*x41) + x151*(-x145*(-2*x177 + 2*x178 + 2*x180) - x150*(-2*x181 + 2*x182 + 2*x183)));
//   const double dlamdxi = 0;
//   const double dphidqop0 = x184*(-x141 + x142 + x144) + x186*(-x146 + x148 + x149) + x187*x91;
//   const double dphidlam0 = x107*x187 + x184*(-x162 + x163 - x165) + x186*(-x166 + x167 - x168);
//   const double dphidphi0 = x117*x187 + x184*(x171 - x173 + x77) + x186*(x175 - x176 + x185);
//   const double dphidxt0 = x118*x187;
//   const double dphidyt0 = x119*x187;
//   const double dphidB = x133*x187 + x184*(-x177 + x178 + x180) + x186*(-x181 + x182 + x183);
//   const double dphidxi = 0;
//   const double dxtdqop0 = x200*x52 - x202*x52 + x203*x74 - x204*x89 + x212*x91;
//   const double dxtdlam0 = -x101*x211 + x105*x208 + x107*x212;
//   const double dxtdphi0 = x115*x208 - x116*x211 + x117*x212;
//   const double dxtdxt0 = W0x*x201 + W0y*x199 + x118*x212;
//   const double dxtdyt0 = x119*x212 - x207*x95 + x210*x94;
//   const double dxtdB = x120*x200 - x120*x202 + x128*x203 - x131*x204 + x133*x212;
//   const double dxtdxi = 0;
//   const double dytdqop0 = x214*x88 - x214*x90 + x215*x73 - x215*x76 - x216*x53 + x216*x58 + x220*x91;
//   const double dytdlam0 = -x101*x218 - x105*x217 + x107*x220 + x219*x98;
//   const double dytdphi0 = x114*x219 - x115*x217 - x116*x218 + x117*x220;
//   const double dytdxt0 = x118*x220 + x214*x36 - x215*x34;
//   const double dytdyt0 = x119*x220 + x20*x216 + x214*x45 + x215*x46;
//   const double dytdB = -x121*x216 + x125*x216 + x127*x215 - x129*x215 + x130*x214 - x132*x214 + x133*x220;
//   const double dytdxi = 0;
//   Eigen::Matrix<double, 5, 7> res;
//   res(0,0) = dqopdqop0;
//   res(0,1) = dqopdlam0;
//   res(0,2) = dqopdphi0;
//   res(0,3) = dqopdxt0;
//   res(0,4) = dqopdyt0;
//   res(0,5) = dqopdB;
//   res(0,6) = dqopdxi;
//   res(1,0) = dlamdqop0;
//   res(1,1) = dlamdlam0;
//   res(1,2) = dlamdphi0;
//   res(1,3) = dlamdxt0;
//   res(1,4) = dlamdyt0;
//   res(1,5) = dlamdB;
//   res(1,6) = dlamdxi;
//   res(2,0) = dphidqop0;
//   res(2,1) = dphidlam0;
//   res(2,2) = dphidphi0;
//   res(2,3) = dphidxt0;
//   res(2,4) = dphidyt0;
//   res(2,5) = dphidB;
//   res(2,6) = dphidxi;
//   res(3,0) = dxtdqop0;
//   res(3,1) = dxtdlam0;
//   res(3,2) = dxtdphi0;
//   res(3,3) = dxtdxt0;
//   res(3,4) = dxtdyt0;
//   res(3,5) = dxtdB;
//   res(3,6) = dxtdxi;
//   res(4,0) = dytdqop0;
//   res(4,1) = dytdlam0;
//   res(4,2) = dytdphi0;
//   res(4,3) = dytdxt0;
//   res(4,4) = dytdyt0;
//   res(4,5) = dytdB;
//   res(4,6) = dytdxi;
//
//   res.col(5) *= 2.99792458e-3;
//
//   return res;
//
//
// }


