#include <sstream>

// Geant4e
#include "TrackPropagation/Geant4e/interface/ConvertFromToCLHEP.h"
#include "TrackPropagation/Geant4e/interface/Geant4ePropagator.h"
// #include "TrackPropagation/Geant4e/interface/G4ErrorPhysicsListCustom.h"

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
}

/** Destructor.
 */
Geant4ePropagator::~Geant4ePropagator() {
  LogDebug("Geant4e") << "Geant4ePropagator::~Geant4ePropagator()" << std::endl;

  // don't close the g4 Geometry here, because the propagator might have been
  // cloned
  // but there is only one, globally-shared Geometry
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
    //     theG4eManager->SetUserInitialization(new G4ErrorPhysicsListCustom());
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

std::tuple<TrajectoryStateOnSurface, Geant4ePropagator::AlgebraicMatrix57, AlgebraicMatrix55> Geant4ePropagator::propagateGenericWithJacobian(const FreeTrajectoryState &ftsStart,
                                                                                const Plane &pDest) const {
                                                                                  
  auto retDefault = []() {
    return std::tuple<TrajectoryStateOnSurface, AlgebraicMatrix57, AlgebraicMatrix55>();
  };
  
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
    return retDefault();

  // re-check propagation direction chosen in case of AnyDirection
  if (mode == G4ErrorMode_PropBackwards && !flipped)
    cmsInitMom = -cmsInitMom;

  CLHEP::Hep3Vector g4InitPos = TrackPropagation::globalPointToHep3Vector(cmsInitPos);
  CLHEP::Hep3Vector g4InitMom = TrackPropagation::globalVectorToHep3Vector(cmsInitMom * GeV);

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
    return retDefault();
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
  
  // initial jacobian is the identity matrix for the state components,
  // and 0 for b-field and material variations
  G4ErrorMatrix jac(5, 7);
  for (unsigned int i = 0; i < 5; ++i) {
    for (unsigned int j = 0; j < 7; ++j) {
      jac(i+1, j+1) = i == j ? 1. : 0.;
    }
  }
  
  const G4ErrorTrajErr errnull(5, 0);
  
  G4ErrorTrajErr g4errorEnd = g4eTrajState.GetError();
  
  G4ErrorTrajErr dQ(5, 0);


  bool continuePropagation = true;
  while (continuePropagation) {
    iterations++;
    LogDebug("Geant4e") << std::endl << "step count " << iterations << " step length " << finalPathLength;

    //set the error matrix to null to disentangle MS and ionization contributions
    g4eTrajState.SetError(errnull);
    const int ierr = theG4eManager->PropagateOneStep(&g4eTrajState, mode);

    if (ierr != 0) {
      // propagation failed, return invalid track state
      return retDefault();
    }

    const float thisPathLength = TrackPropagation::g4doubleToCmsDouble(g4eTrajState.GetG4Track()->GetStepLength());

    const G4ErrorMatrix &transferMatrix = g4eTrajState.GetTransfMat();
    
    // transport contribution to error
    g4errorEnd = g4errorEnd.similarity(transferMatrix).T();
    dQ = dQ.similarity(transferMatrix).T();
    
    // MS and ionization contributions to error
    
    G4ErrorTrajErr errMSI = g4eTrajState.GetError();
    
    // recompute ionization contribution (without truncation for the moment)
    errMSI(0+1, 0+1) = computeErrorIoni(g4eTrajState.GetG4Track());
    
    g4errorEnd += errMSI;
    
    // transport  + nominal energy loss contribution to jacobian
    jac = transferMatrix*jac;
    
    //TODO assess relevance of approximations (does the order matter? position/momentum before or after step?)
    // local b-field contribution to jacobian 
    const GlobalPoint pos = TrackPropagation::hepPoint3DToGlobalPoint(g4eTrajState.GetPosition());
    const GlobalVector mom = TrackPropagation::hep3VectorToGlobalVector(g4eTrajState.GetMomentum() / GeV);
    const GlobalVector momPre = TrackPropagation::hep3VectorToGlobalVector(g4eTrajState.GetG4Track()->GetStep()->GetPreStepPoint()->GetMomentum() / GeV);
    const double B = theField->inTesla(pos).mag();
    const double pinv = 1./mom.mag();
    
    for (unsigned int i = 1; i < 5; ++i) {
      jac(i+1, 5+1) += transferMatrix(i+1, 0+1)*pinv/B;
    }
    
    // local material contribution to jacobian
    const double ePre = g4eTrajState.GetG4Track()->GetStep()->GetPreStepPoint()->GetTotalEnergy() / GeV;
    const double ePost = g4eTrajState.GetG4Track()->GetStep()->GetPostStepPoint()->GetTotalEnergy() / GeV;
    const double dE = ePost - ePre;
        
//     std::cout << "ePost = " << ePost << " p = " << mom.mag() << std::endl;
    
    jac(0+1, 6+1) += -ePost*std::pow(pinv, 3)*dE;
    
    
    const double pPre = momPre.mag();
    const double betaPre = pPre/ePre;
    const double mass  = g4eTrajState.GetG4Track()->GetDynamicParticle()->GetMass() / GeV;
    
    G4ErrorTrajErr errMS = errMSI;
    errMS(0+1, 0+1) = 0.;
//     std::cout << "ppre = " << pPre << " betaPre = " << betaPre << " mass = " << mass << " ePre  = " << ePre << std::endl;
    dQ += 2.*pPre*(betaPre*betaPre + 2.*mass*mass/ePre/ePre)*errMS;
//     dQ += 2.*pPre*(betaPre*betaPre + 2.*mass*mass/ePre/ePre)*errMS;
    
//     std::cout << "dE = " << dE << std::endl;
    
    
//     std::cout << "jac" << std::endl;
//     std::cout << jac << std::endl;
    
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
  GlobalVector momEndGV = TrackPropagation::hep3VectorToGlobalVector(momEnd) / GeV;

  debugReportTrackState("final", posEndGV, finalRecoPos, momEndGV, momEnd, pDest);

  // Get the error covariance matrix from Geant4e. It comes in curvilinear
  // coordinates so use the appropiate CMS class
//   G4ErrorTrajErr g4errorEnd = g4eTrajState.GetError();

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

//   GlobalTrajectoryParameters tParsDest(posEndGV, momEndGV, ftsStart.charge(), theField);
  
//   auto const localf = pDest.toLocal(posEndGV);
//   auto const locald = pDest.toLocal<double>(posEndGV);
  
//   std::cout << "localf = " << localf << std::endl;
//   std::cout << "locald = " << locald << std::endl;
//   std::cout << "diff = " << localf-locald << std::endl;
  
  auto const localpos = pDest.toLocal<double>(posEndGV);
  auto const localmom = pDest.toLocal<double>(momEndGV);
  
  const LocalPoint localposcor(localpos.x() - localpos.z()*localmom.x()/localmom.z(), localpos.y() - localpos.z()*localmom.y()/localmom.z(), 0.);
  
  const LocalTrajectoryParameters localparms(localposcor, localmom, ftsStart.charge());
  
  const GlobalPoint posEndcor = pDest.toGlobal<double>(localposcor);
  
  GlobalTrajectoryParameters tParsDest(posEndcor, momEndGV, ftsStart.charge(), theField);
  

  SurfaceSideDefinition::SurfaceSide side;

  side = propagationDirection() == alongMomentum ? SurfaceSideDefinition::afterSurface
                                                 : SurfaceSideDefinition::beforeSurface;

  AlgebraicMatrix57 jacfinal;
  // convert from 1/p to q/p
  for (unsigned int i = 0; i < 5; ++i) {
    for (unsigned int j = 0; j < 7; ++j) {
      jacfinal(i, j) = jac(i+1, j+1);
      if (i==0) {
        jacfinal(i, j) *= ftsStart.charge();
      }
      if (j==0) {
        jacfinal(i, j) *= ftsStart.charge();
      }
    }
  }
  
  AlgebraicMatrix55 dQfinal;
  for (unsigned int i = 0; i < 5; ++i) {
    for (unsigned int j = 0; j < 5; ++j) {
      dQfinal(i, j) = ftsStart.charge()*dQ(i+1, j+1);
    }
  }
  
//   JacobianCurvilinearToLocal curv2Loc(pDest, localparms, tParsDest, *theField);
  const AlgebraicMatrix55 curv2Loc = curv2localJacobianAlt(tParsDest, pDest);
  auto const localerr = ROOT::Math::Similarity(curv2Loc, curvError.matrix());


//   return std::tuple<TrajectoryStateOnSurface, AlgebraicMatrix57, AlgebraicMatrix55>(TrajectoryStateOnSurface(tParsDest, curvError, pDest, side), jacfinal, dQfinal);
  
  return std::tuple<TrajectoryStateOnSurface, AlgebraicMatrix57, AlgebraicMatrix55>(TrajectoryStateOnSurface(localparms, localerr, pDest, theField, side), jacfinal, dQfinal);
}

std::tuple<TrajectoryStateOnSurface, Geant4ePropagator::AlgebraicMatrix57, AlgebraicMatrix55, double> Geant4ePropagator::propagateGenericWithJacobianAlt(const FreeTrajectoryState &ftsStart,
                                                                                const Plane &pDest) const {
                          
  using namespace Eigen;

                                                                                  
  auto retDefault = []() {
    return std::tuple<TrajectoryStateOnSurface, AlgebraicMatrix57, AlgebraicMatrix55, double>();
  };
  
//   std::cout << "start propagation" << std::endl;
  
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
    return retDefault();

  // re-check propagation direction chosen in case of AnyDirection
  if (mode == G4ErrorMode_PropBackwards && !flipped)
    cmsInitMom = -cmsInitMom;

  CLHEP::Hep3Vector g4InitPos = TrackPropagation::globalPointToHep3Vector(cmsInitPos);
  CLHEP::Hep3Vector g4InitMom = TrackPropagation::globalVectorToHep3Vector(cmsInitMom * GeV);

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
    return retDefault();
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
  const AlgebraicMatrix55 errstart = ftsStart.curvilinearError().matrix();
  Matrix<double, 5, 5> g4errorEnd = Map<const Matrix<double, 5, 5, RowMajor>>(errstart.Array());
  
  

  
  G4ErrorTrajErr dQ(5, 0);

  double dEdxlast = 0.;
  double masslast = 0.;

  bool continuePropagation = true;
  while (continuePropagation) {
    iterations++;
    LogDebug("Geant4e") << std::endl << "step count " << iterations << " step length " << finalPathLength;

    const GlobalPoint pos = TrackPropagation::hepPoint3DToGlobalPoint(g4eTrajState.GetPosition());
    const GlobalVector mom = TrackPropagation::hep3VectorToGlobalVector(g4eTrajState.GetMomentum() / GeV);
//     constexpr double mass = 0.1056583745;
    
    const FreeTrajectoryState statepre(pos, mom, ftsStart.charge(), &ftsStart.parameters().magneticField());
    
    //set the error matrix to null to disentangle MS and ionization contributions
    g4eTrajState.SetError(errnull);
    const int ierr = theG4eManager->PropagateOneStep(&g4eTrajState, mode);

    if (ierr != 0) {
      // propagation failed, return invalid track state
      return retDefault();
    }

    const double thisPathLength = TrackPropagation::g4doubleToCmsDouble(g4eTrajState.GetG4Track()->GetStepLength());

    const double ePre = g4eTrajState.GetG4Track()->GetStep()->GetPreStepPoint()->GetTotalEnergy() / GeV;
    const double ePost = g4eTrajState.GetG4Track()->GetStep()->GetPostStepPoint()->GetTotalEnergy() / GeV;
    const double dEdx = (ePost - ePre)/thisPathLength;
    const double mass  = g4eTrajState.GetG4Track()->GetDynamicParticle()->GetMass() / GeV;

    if (thisPathLength > 0.) {
      dEdxlast = dEdx;
      masslast = mass;
    }
    
//     std::cout << "thisPathLength = " << thisPathLength << std::endl;
//     const Matrix<double, 5, 7> transportJac = transportJacobian(statepre, thisPathLength, dEdx, mass);
    const Matrix<double, 5, 7> transportJac = transportJacobianBz(statepre, thisPathLength, dEdx, mass);
//     const Matrix<double, 5, 7> transportJac = transportJacobianBzAdvanced(statepre, thisPathLength, dEdx, mass);
    
//     const Matrix<double, 5, 7> transportJacAlt = transportJacobianBz(statepre, thisPathLength, dEdx, mass);

//     std::cout << "transportJac" << std::endl;
//     std::cout << transportJac << std::endl;
//     std::cout << "transportJacAlt" << std::endl;
//     std::cout << transportJacAlt << std::endl;
    
    
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
    
    
//     const G4ErrorMatrix &transferMatrix = g4eTrajState.GetTransfMat();
    
    // transport contribution to error
//     g4errorEnd = g4errorEnd.similarity(transferMatrix).T();
//     dQ = dQ.similarity(transferMatrix).T();
    
    // MS and ionization contributions to error
    
    G4ErrorTrajErr errMSI = g4eTrajState.GetError();
    
    // recompute ionization contribution (without truncation for the moment)
    errMSI(0+1, 0+1) = computeErrorIoni(g4eTrajState.GetG4Track());
    
//     std::cout << "errMSI" << std::endl;
//     std::cout << errMSI << std::endl;
    
    for (unsigned int i = 0; i < 5; i++) {
      for (unsigned int j = 0; j < 5; j++) {
        double w = 1.;
        if (i==0) {
          w *= ftsStart.charge();
        }
        if (j==0) {
          w *= ftsStart.charge();
        }
        g4errorEnd(i, j) += w*errMSI(i+1, j+1);
      }
    }
      
      
    
//     g4errorEnd += errMSI;
    
    // transport  + nominal energy loss contribution to jacobian
//     jac = transferMatrix*jac;
    jac = (transportJac.leftCols<5>()*jac).eval();

    //TODO assess relevance of approximations (does the order matter? position/momentum before or after step?)
    //b-field and material contributions to jacobian
    jac.rightCols<2>() += transportJac.rightCols<2>();
        
    
    //TODO assess relevance of approximations (does the order matter? position/momentum before or after step?)
//     // local b-field contribution to jacobian 
//     const GlobalPoint pos = TrackPropagation::hepPoint3DToGlobalPoint(g4eTrajState.GetPosition());
//     const GlobalVector mom = TrackPropagation::hep3VectorToGlobalVector(g4eTrajState.GetMomentum() / GeV);
//     const GlobalVector momPre = TrackPropagation::hep3VectorToGlobalVector(g4eTrajState.GetG4Track()->GetStep()->GetPreStepPoint()->GetMomentum() / GeV);
//     const double B = theField->inTesla(pos).mag();
//     const double pinv = 1./mom.mag();
//     
//     for (unsigned int i = 1; i < 5; ++i) {
//       jac(i+1, 5+1) += transferMatrix(i+1, 0+1)*pinv/B;
//     }
    
    // local material contribution to jacobian
//     const double ePre = g4eTrajState.GetG4Track()->GetStep()->GetPreStepPoint()->GetTotalEnergy() / GeV;
//     const double ePost = g4eTrajState.GetG4Track()->GetStep()->GetPostStepPoint()->GetTotalEnergy() / GeV;
//     const double dE = ePost - ePre;
        
//     std::cout << "ePost = " << ePost << " p = " << mom.mag() << std::endl;
    
//     jac(0+1, 6+1) += -ePost*std::pow(pinv, 3)*dE;
    
    
//     const double pPre = momPre.mag();
//     const double betaPre = pPre/ePre;
//     const double mass  = g4eTrajState.GetG4Track()->GetDynamicParticle()->GetMass() / GeV;
    
//     G4ErrorTrajErr errMS = errMSI;
//     errMS(0+1, 0+1) = 0.;
//     std::cout << "ppre = " << pPre << " betaPre = " << betaPre << " mass = " << mass << " ePre  = " << ePre << std::endl;
//     dQ += 2.*pPre*(betaPre*betaPre + 2.*mass*mass/ePre/ePre)*errMS;
//     dQ += 2.*pPre*(betaPre*betaPre + 2.*mass*mass/ePre/ePre)*errMS;
    
//     std::cout << "dE = " << dE << std::endl;
    
    
//     std::cout << "jac" << std::endl;
//     std::cout << jac << std::endl;
    
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
  const HepGeom::Vector3D<double> momEnd = g4eTrajState.GetMomentum();

  // use the hit on the the RECO plane as the final position to be d'accor with
  // the RecHit measurements
  const GlobalPoint posEndGV = TrackPropagation::hepPoint3DToGlobalPoint(finalRecoPos);
  GlobalVector momEndGV = TrackPropagation::hep3VectorToGlobalVector(momEnd) / GeV;

  debugReportTrackState("final", posEndGV, finalRecoPos, momEndGV, momEnd, pDest);

  // Get the error covariance matrix from Geant4e. It comes in curvilinear
  // coordinates so use the appropiate CMS class
//   G4ErrorTrajErr g4errorEnd = g4eTrajState.GetError();

//   AlgebraicMatrix55 endErrorMatrix;
//   Map<Matrix<double, 5, 5, RowMajor>>(endErrorMatrix.Array()) = g4errorEnd;
  
  AlgebraicSymMatrix55 endErrorMatrixSym;
  for (unsigned int i = 0; i < 5; i++) {
    for (unsigned int j = 0; j < 5; j++) {
      endErrorMatrixSym(i, j) = g4errorEnd(i, j);
    }
  } 
  
  CurvilinearTrajectoryError curvError(endErrorMatrixSym);
  
//   CurvilinearTrajectoryError curvError(
//       TrackPropagation::g4ErrorTrajErrToAlgebraicSymMatrix55(g4errorEnd, ftsStart.charge()));

  if (mode == G4ErrorMode_PropBackwards) {
    GlobalTrajectoryParameters endParm(
        posEndGV, momEndGV, ftsStart.parameters().charge(), &ftsStart.parameters().magneticField());

    // flip the momentum direction because it has been flipped before running
    // G4's backwards prop
    momEndGV = -momEndGV;
  }

  LogDebug("Geant4e") << "G4e -  Error matrix after propagation: " << std::endl << g4errorEnd;

  LogDebug("Geant4e") << "CMS -  Error matrix after propagation: " << std::endl << curvError.matrix();

//   GlobalTrajectoryParameters tParsDest(posEndGV, momEndGV, ftsStart.charge(), theField);
  
//   auto const localf = pDest.toLocal(posEndGV);
//   auto const locald = pDest.toLocal<double>(posEndGV);
  
//   std::cout << "localf = " << localf << std::endl;
//   std::cout << "locald = " << locald << std::endl;
//   std::cout << "diff = " << localf-locald << std::endl;
  
  auto const localpos = pDest.toLocal<double>(posEndGV);
  auto const localmom = pDest.toLocal<double>(momEndGV);
  
//   std::cout << "local z = " << localpos.z() << " local pz = " << localmom.z() << std::endl;
  
  const LocalPoint localposcor(localpos.x() - localpos.z()*localmom.x()/localmom.z(), localpos.y() - localpos.z()*localmom.y()/localmom.z(), 0.);
  
  const LocalTrajectoryParameters localparms(localposcor, localmom, ftsStart.charge());
  
  const GlobalPoint posEndcor = pDest.toGlobal<double>(localposcor);
  
  GlobalTrajectoryParameters tParsDest(posEndcor, momEndGV, ftsStart.charge(), theField);
  

  SurfaceSideDefinition::SurfaceSide side;

  side = propagationDirection() == alongMomentum ? SurfaceSideDefinition::afterSurface
                                                 : SurfaceSideDefinition::beforeSurface;

  AlgebraicMatrix57 jacfinal;
  Map<Matrix<double, 5, 7, RowMajor>>(jacfinal.Array()) = jac;
  
  
  // convert from 1/p to q/p
//   for (unsigned int i = 0; i < 5; ++i) {
//     for (unsigned int j = 0; j < 7; ++j) {
//       jacfinal(i, j) = jac(i+1, j+1);
//       if (i==0) {
//         jacfinal(i, j) *= ftsStart.charge();
//       }
//       if (j==0) {
//         jacfinal(i, j) *= ftsStart.charge();
//       }
//     }
//   }
  
  AlgebraicMatrix55 dQfinal;
//   for (unsigned int i = 0; i < 5; ++i) {
//     for (unsigned int j = 0; j < 5; ++j) {
//       dQfinal(i, j) = ftsStart.charge()*dQ(i+1, j+1);
//     }
//   }
  
//   JacobianCurvilinearToLocal curv2Loc(pDest, localparms, tParsDest, *theField);
//   const AlgebraicMatrix55 curv2Loc = curv2localJacobianAlt(tParsDest, pDest);
//   const AlgebraicMatrix55 curv2Loceloss = curv2localJacobianAlteloss(tParsDest, pDest, dEdxlast, masslast);
//   auto const localerr = ROOT::Math::Similarity(curv2Loc, curvError.matrix());
//   auto const localerreloss = ROOT::Math::Similarity(curv2Loceloss, curvError.matrix());
//   
//   std::cout << "curv2Loc" << std::endl;
//   std::cout << curv2Loc << std::endl;
//   std::cout << "curv2Loceloss" << std::endl;
//   std::cout << curv2Loceloss << std::endl;
//   std::cout << "localerr" << std::endl;
//   std::cout << localerr << std::endl;
//   std::cout << "localerreloss" << std::endl;
//   std::cout << localerreloss << std::endl;
  
  const AlgebraicMatrix55 curv2Loc = curv2localJacobianAlteloss(tParsDest, pDest, dEdxlast, masslast);
  auto const localerr = ROOT::Math::Similarity(curv2Loc, curvError.matrix());

//   return std::tuple<TrajectoryStateOnSurface, AlgebraicMatrix57, AlgebraicMatrix55>(TrajectoryStateOnSurface(tParsDest, curvError, pDest, side), jacfinal, dQfinal);
  
  return std::tuple<TrajectoryStateOnSurface, AlgebraicMatrix57, AlgebraicMatrix55, double>(TrajectoryStateOnSurface(localparms, localerr, pDest, theField, side), jacfinal, dQfinal, dEdxlast);
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
double Geant4ePropagator::computeErrorIoni(const G4Track* aTrack) const
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

  G4double Etot  = aTrack->GetTotalEnergy() / GeV;
  G4double beta  = aTrack->GetMomentum().mag() / GeV / Etot;
  G4double mass  = aTrack->GetDynamicParticle()->GetMass() / GeV;
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
  G4double dedxSq =
    XI * Emax * (1. - (beta * beta / 2.)) * 1.E-12;  // now in GeV^2
  /*The above  formula for var(1/p) good for dens scatterers. However, for MIPS
    passing through a gas it leads to overestimation. Further more for incident
    electrons the Emax is almost equal to incident energy. This leads  to
    k=Xi/Emax  as small as e-6  and gradually the cov matrix  explodes.
    http://www2.pv.infn.it/~rotondi/kalman_1.pdf
    Since I do not have enough info at the moment to implement Landau &
    sub-Landau models for k=Xi/Emax <0.01 I'll saturate k at this value for now
  */
  
//   dedxSq *= 0.;
  
  // Implementation based on PANDA Report PV/01-07 Section 7
  // TODO understand implications for thick scatterer split into small steps
  
//   G4double dedxSq;
  
//   const double kappa = XI/Emax;
//   std::cout << "steplength = " << stepLengthCm << " xi = " << XI*1e-6 << " xi/dx = " << XI*1e-6/stepLengthCm << " emax = " << Emax*1e-6 << " kappa = " << kappa << " effZ = " << effZ << std::endl;
  
//   if (kappa > 0.005) {
//     //vavilov distribution (or gaussian limit which is equivalent for the variance)
//   
//     dedxSq = XI * Emax * (1. - (beta * beta / 2.)) * 1.E-12;  // now in GeV^2
//   }
//   else {
//   
//   
//     const double I = 16.*pow(effZ,0.9);
//     const double f2 = effZ <= 2. ? 0. : 2./effZ;
//     const double f1 = 1. - f2;
//     const double e2 = 10.*effZ*effZ;
//     const double e1 = pow(I/pow(e2,f2),1./f1);
//     const double r = 0.4;
//     const double emaxev = Emax*1e6;
//     const double massev = mass*1e9;
//     
//     const double ePre = aTrack->GetStep()->GetPreStepPoint()->GetTotalEnergy() / GeV;
//     const double ePost = aTrack->GetStep()->GetPostStepPoint()->GetTotalEnergy() / GeV;
//     const double C = (ePre-ePost)/stepLengthCm*1e9;
//     
//     const double sigma1 = C*f1*(log(2.*massev*beta*beta*gamma*gamma/e1) - beta*beta)/e1/(log(2.*massev*beta*beta*gamma*gamma/I) - beta*beta)*(1.-r);
//     
//     const double sigma2 = C*f1*(log(2.*massev*beta*beta*gamma*gamma/e2) - beta*beta)/e2/(log(2.*massev*beta*beta*gamma*gamma/I) - beta*beta)*(1.-r);
//     
//     const double sigma3 = C*emaxev/I/(emaxev+I)/log((emaxev+I)/I)*r;
//     
//     const double Nc = (sigma1 + sigma2 + sigma3)*stepLengthCm;
//     
// //     std::cout << "Nc = " << Nc << std::endl;
//     
//     if (Nc > 50.) {
// //     if (false) {
//       //landau
//       dedxSq = 15.76*15.76*XI*XI*1e-12; //corresponds to 0.996 quantile
//       
//       const double dedxsqvavilov = XI * Emax * (1. - (beta * beta / 2.)) * 1.E-12;
//       
//       std::cout << "dedxsq: vavilov = " << dedxsqvavilov << " landau = " << dedxSq << std::endl;
//     }
//     else {
//       //sub-landau
//       const double alpha = 0.996;
// //       const double alpha = 1.;
//       const double ealpha = I/(1. - alpha*emaxev/(emaxev + I));
//       const double e3 = I*(emaxev +I)*log(ealpha/I)/emaxev;
//       const double e3sq = I*(emaxev + I)*(ealpha - I)/emaxev;
//       const double sigmae3sq = e3sq  - e3*e3;
//       
//       dedxSq = sigma1*stepLengthCm*e1*e1 + sigma2*stepLengthCm*e2*e2 + sigma3*stepLengthCm*e3*e3 + sigma3*stepLengthCm*sigmae3sq*(sigma3*stepLengthCm + 1.);
//       dedxSq *= 1e-18;
//       
//       const double dedxsqlandau = 15.76*15.76*XI*XI*1e-12;
//       const double dedxsqvavilov = XI * Emax * (1. - (beta * beta / 2.)) * 1.E-12;
// //       
//       std::cout << "dedxsq:  vavilov = " << dedxsqvavilov << " landau = " << dedxsqlandau << " sublandau = " << dedxSq << std::endl;;
//       
//     }
//     
//   
//   }
  
  
//   dedxSq = 1.7*1.7*XI*XI*1e-12;

//   if (false) {
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
                                                   G4double& effA) const
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

AlgebraicMatrix55 Geant4ePropagator::curv2localJacobianAlt(const GlobalTrajectoryParameters &globalSource, const Surface &surface) const {
  
  const GlobalVector &bfield = globalSource.magneticFieldInInverseGeV();
  
  CurvilinearTrajectoryParameters curvparms(globalSource.position(), globalSource.momentum(), globalSource.charge());

  const double qop0 = curvparms.Qbp();
  const double lam0 = curvparms.lambda();
  const double phi0 = curvparms.phi();
  const double xt0 = curvparms.xT();
  const double yt0 = curvparms.yT();
  
  const double p0 = globalSource.momentum().mag();
  const GlobalVector W0 = globalSource.momentum()/p0;
  const double W0x = W0.x();
  const double W0y = W0.y();
  const double W0z = W0.z();
    
  const LocalVector lx(1.,0.,0.);
  const LocalVector ly(0.,1.,0.);
  const LocalVector lz(0.,0.,1.);
  const LocalPoint l0(0., 0.);
  
  const GlobalVector I1 = surface.toGlobal<double>(lz);
  const GlobalVector J1 = surface.toGlobal<double>(lx);
  const GlobalVector K1 = surface.toGlobal<double>(ly);  
  const GlobalPoint r1 = surface.toGlobal<double>(l0);
  
  const double Ix1 = I1.x();
  const double Iy1 = I1.y();
  const double Iz1 = I1.z();
  
  const double Jx1 = J1.x();
  const double Jy1 = J1.y();
  const double Jz1 = J1.z();
  
  const double Kx1 = K1.x();
  const double Ky1 = K1.y();
  const double Kz1 = K1.z();
  
  const double rx1 = r1.x();
  const double ry1 = r1.y();
  const double rz1 = r1.z();
    
  const double B = bfield.mag();
  const GlobalVector H = bfield/B;
  const double hx = H.x();
  const double hy = H.y();
  const double hz = H.z();
      
  const double x0 = std::sin(lam0);
  const double x1 = Iz1*x0;
  const double x2 = std::cos(lam0);
  const double x3 = std::cos(phi0);
  const double x4 = Ix1*x3;
  const double x5 = std::sin(phi0);
  const double x6 = Iy1*x5;
  const double x7 = x1 + x2*x4 + x2*x6;
  const double x8 = std::pow(x7, -2);
  const double x9 = Iz1*Jx1;
  const double x10 = Iz1*Jy1;
  const double x11 = lam0 + phi0;
  const double x12 = -phi0;
  const double x13 = lam0 + x12;
  const double x14 = std::pow(Ix1*std::cos(x11) + Ix1*std::cos(x13) + Iy1*std::sin(x11) - Iy1*std::sin(x13) + 2*x1, -2);
  const double x15 = 2*Ix1;
  const double x16 = Jy1*x15;
  const double x17 = 2*Iy1;
  const double x18 = Jx1*x17;
  const double x19 = 2*lam0;
  const double x20 = std::cos(x19);
  const double x21 = phi0 + x19;
  const double x22 = std::cos(x21);
  const double x23 = std::sin(x21);
  const double x24 = Ix1*Jz1;
  const double x25 = Iy1*Jz1;
  const double x26 = x12 + x19;
  const double x27 = std::cos(x26);
  const double x28 = std::sin(x26);
  const double x29 = Ix1*W0y - Iy1*W0x;
  const double x30 = hz*x2;
  const double x31 = hx*x0 - x3*x30;
  const double x32 = hy*x0 - x30*x5;
  const double x33 = x2*(hx*x5 - hy*x3);
  const double x34 = x2*x3;
  const double x35 = x2*x5;
  const double x36 = Jx1*x34 + Jy1*x35 + Jz1*x0;
  const double x37 = -Ix1*x32 + Iy1*x31 - Iz1*x33;
  const double x38 = std::pow(W0x, 2) + std::pow(W0y, 2);
  const double x39 = std::pow(x38, -1.0/2.0);
  const double x40 = B*qop0*x39/std::pow(x7, 3);
  const double x41 = x40*(-x36*x37 + x7*(-Jx1*x32 + Jy1*x31 - Jz1*x33));
  const double x42 = Ix1*W0x*W0z + Iy1*W0y*W0z - Iz1*x38;
  const double x43 = Iz1*Kx1;
  const double x44 = Iz1*Ky1;
  const double x45 = Ky1*x15;
  const double x46 = Kx1*x17;
  const double x47 = Ix1*Kz1;
  const double x48 = Iy1*Kz1;
  const double x49 = Kx1*x34 + Ky1*x35 + Kz1*x0;
  const double x50 = x40*(-x37*x49 + x7*(-Kx1*x32 + Ky1*x31 - Kz1*x33));
  const double x51 = x39/x7;
  const double x52 = x38*x7;
  const double x53 = W0z*x7;
  const double dqopdqop0 = 1;
  const double dqopdlam0 = 0;
  const double dqopdphi0 = 0;
  const double dqopdxt0 = 0;
  const double dqopdyt0 = 0;
  const double ddxdzdqop0 = 0;
  const double ddxdzdlam0 = x8*(Jz1*x4 + Jz1*x6 - x10*x5 - x3*x9);
  const double ddxdzdphi0 = x14*(x10*x23 + x10*x28 + x16*x20 + x16 - x18*x20 - x18 - x22*x24 + x22*x9 - x23*x25 + x24*x27 - x25*x28 - x27*x9);
  const double ddxdzdxt0 = x29*x41;
  const double ddxdzdyt0 = x41*x42;
  const double ddydzdqop0 = 0;
  const double ddydzdlam0 = x8*(Kz1*x4 + Kz1*x6 - x3*x43 - x44*x5);
  const double ddydzdphi0 = x14*(x20*x45 - x20*x46 + x22*x43 - x22*x47 + x23*x44 - x23*x48 - x27*x43 + x27*x47 + x28*x44 - x28*x48 + x45 - x46);
  const double ddydzdxt0 = x29*x50;
  const double ddydzdyt0 = x42*x50;
  const double dxdqop0 = 0;
  const double dxdlam0 = 0;
  const double dxdphi0 = 0;
  const double dxdxt0 = x51*(x29*x36 + x7*(-Jx1*W0y + Jy1*W0x));
  const double dxdyt0 = x51*(Jz1*x52 + x36*x42 - x53*(Jx1*W0x + Jy1*W0y));
  const double dydqop0 = 0;
  const double dydlam0 = 0;
  const double dydphi0 = 0;
  const double dydxt0 = x51*(x29*x49 + x7*(-Kx1*W0y + Ky1*W0x));
  const double dydyt0 = x51*(Kz1*x52 + x42*x49 - x53*(Kx1*W0x + Ky1*W0y));
  AlgebraicMatrix55 res;
  res(0,0) = dqopdqop0;
  res(0,1) = dqopdlam0;
  res(0,2) = dqopdphi0;
  res(0,3) = dqopdxt0;
  res(0,4) = dqopdyt0;
  res(1,0) = ddxdzdqop0;
  res(1,1) = ddxdzdlam0;
  res(1,2) = ddxdzdphi0;
  res(1,3) = ddxdzdxt0;
  res(1,4) = ddxdzdyt0;
  res(2,0) = ddydzdqop0;
  res(2,1) = ddydzdlam0;
  res(2,2) = ddydzdphi0;
  res(2,3) = ddydzdxt0;
  res(2,4) = ddydzdyt0;
  res(3,0) = dxdqop0;
  res(3,1) = dxdlam0;
  res(3,2) = dxdphi0;
  res(3,3) = dxdxt0;
  res(3,4) = dxdyt0;
  res(4,0) = dydqop0;
  res(4,1) = dydlam0;
  res(4,2) = dydphi0;
  res(4,3) = dydxt0;
  res(4,4) = dydyt0;

  
  return res;
                                              
}

AlgebraicMatrix55 Geant4ePropagator::curv2localJacobianAlteloss(const GlobalTrajectoryParameters &globalSource, const Surface &surface, double dEdx, double mass) const {
  
  const GlobalVector &bfield = globalSource.magneticFieldInInverseGeV();
  
  CurvilinearTrajectoryParameters curvparms(globalSource.position(), globalSource.momentum(), globalSource.charge());

  const double qop0 = curvparms.Qbp();
  const double lam0 = curvparms.lambda();
  const double phi0 = curvparms.phi();
  const double xt0 = curvparms.xT();
  const double yt0 = curvparms.yT();
  
  const double q = globalSource.charge();
  
  const double p0 = globalSource.momentum().mag();
  const GlobalVector W0 = globalSource.momentum()/p0;
  const double W0x = W0.x();
  const double W0y = W0.y();
  const double W0z = W0.z();
    
  const LocalVector lx(1.,0.,0.);
  const LocalVector ly(0.,1.,0.);
  const LocalVector lz(0.,0.,1.);
  const LocalPoint l0(0., 0.);
  
  const GlobalVector I1 = surface.toGlobal<double>(lz);
  const GlobalVector J1 = surface.toGlobal<double>(lx);
  const GlobalVector K1 = surface.toGlobal<double>(ly);  
  const GlobalPoint r1 = surface.toGlobal<double>(l0);
  
  const double Ix1 = I1.x();
  const double Iy1 = I1.y();
  const double Iz1 = I1.z();
  
  const double Jx1 = J1.x();
  const double Jy1 = J1.y();
  const double Jz1 = J1.z();
  
  const double Kx1 = K1.x();
  const double Ky1 = K1.y();
  const double Kz1 = K1.z();
  
  const double rx1 = r1.x();
  const double ry1 = r1.y();
  const double rz1 = r1.z();
    
  const double B = bfield.mag();
  const GlobalVector H = bfield/B;
  const double hx = H.x();
  const double hy = H.y();
  const double hz = H.z();
        
  const double x0 = std::pow(q, 2);
  const double x1 = std::pow(x0, -3.0/2.0);
  const double x2 = Iy1*W0x;
  const double x3 = Ix1*W0y;
  const double x4 = std::pow(qop0, 2);
  const double x5 = std::pow(W0x, 2) + std::pow(W0y, 2);
  const double x6 = std::pow(x5, -1.0/2.0);
  const double x7 = std::sin(lam0);
  const double x8 = Iz1*x7;
  const double x9 = std::cos(lam0);
  const double x10 = std::cos(phi0);
  const double x11 = Ix1*x10;
  const double x12 = std::sin(phi0);
  const double x13 = Iy1*x12;
  const double x14 = x11*x9 + x13*x9 + x8;
  const double x15 = x6/x14;
  const double x16 = dEdx*q*x1*x15*x4*std::sqrt(std::pow(mass, 2)*x4 + x0);
  const double x17 = Ix1*W0x*W0z + Iy1*W0y*W0z - Iz1*x5;
  const double x18 = std::pow(x14, -2);
  const double x19 = Iz1*Jx1;
  const double x20 = Iz1*Jy1;
  const double x21 = lam0 + phi0;
  const double x22 = -phi0;
  const double x23 = lam0 + x22;
  const double x24 = std::pow(Ix1*std::cos(x21) + Ix1*std::cos(x23) + Iy1*std::sin(x21) - Iy1*std::sin(x23) + 2*x8, -2);
  const double x25 = 2*Ix1;
  const double x26 = Jy1*x25;
  const double x27 = 2*Iy1;
  const double x28 = Jx1*x27;
  const double x29 = 2*lam0;
  const double x30 = std::cos(x29);
  const double x31 = phi0 + x29;
  const double x32 = std::cos(x31);
  const double x33 = std::sin(x31);
  const double x34 = Ix1*Jz1;
  const double x35 = Iy1*Jz1;
  const double x36 = x22 + x29;
  const double x37 = std::cos(x36);
  const double x38 = std::sin(x36);
  const double x39 = -x2 + x3;
  const double x40 = hz*x9;
  const double x41 = hx*x7 - x10*x40;
  const double x42 = hy*x7 - x12*x40;
  const double x43 = x9*(hx*x12 - hy*x10);
  const double x44 = x10*x9;
  const double x45 = x12*x9;
  const double x46 = Jx1*x44 + Jy1*x45 + Jz1*x7;
  const double x47 = -Ix1*x42 + Iy1*x41 - Iz1*x43;
  const double x48 = B*qop0*x6/std::pow(x14, 3);
  const double x49 = x48*(x14*(-Jx1*x42 + Jy1*x41 - Jz1*x43) - x46*x47);
  const double x50 = Iz1*Kx1;
  const double x51 = Iz1*Ky1;
  const double x52 = Ky1*x25;
  const double x53 = Kx1*x27;
  const double x54 = Ix1*Kz1;
  const double x55 = Iy1*Kz1;
  const double x56 = Kx1*x44 + Ky1*x45 + Kz1*x7;
  const double x57 = x48*(x14*(-Kx1*x42 + Ky1*x41 - Kz1*x43) - x47*x56);
  const double x58 = x14*x5;
  const double x59 = W0z*x14;
  const double dqopdqop0 = std::pow(q, 3)*x1*std::fabs(qop0)/qop0;
  const double dqopdlam0 = 0;
  const double dqopdphi0 = 0;
  const double dqopdxt0 = x16*(x2 - x3);
  const double dqopdyt0 = -x16*x17;
  const double ddxdzdqop0 = 0;
  const double ddxdzdlam0 = x18*(Jz1*x11 + Jz1*x13 - x10*x19 - x12*x20);
  const double ddxdzdphi0 = x24*(x19*x32 - x19*x37 + x20*x33 + x20*x38 + x26*x30 + x26 - x28*x30 - x28 - x32*x34 - x33*x35 + x34*x37 - x35*x38);
  const double ddxdzdxt0 = x39*x49;
  const double ddxdzdyt0 = x17*x49;
  const double ddydzdqop0 = 0;
  const double ddydzdlam0 = x18*(Kz1*x11 + Kz1*x13 - x10*x50 - x12*x51);
  const double ddydzdphi0 = x24*(x30*x52 - x30*x53 + x32*x50 - x32*x54 + x33*x51 - x33*x55 - x37*x50 + x37*x54 + x38*x51 - x38*x55 + x52 - x53);
  const double ddydzdxt0 = x39*x57;
  const double ddydzdyt0 = x17*x57;
  const double dxdqop0 = 0;
  const double dxdlam0 = 0;
  const double dxdphi0 = 0;
  const double dxdxt0 = x15*(x14*(-Jx1*W0y + Jy1*W0x) + x39*x46);
  const double dxdyt0 = x15*(Jz1*x58 + x17*x46 - x59*(Jx1*W0x + Jy1*W0y));
  const double dydqop0 = 0;
  const double dydlam0 = 0;
  const double dydphi0 = 0;
  const double dydxt0 = x15*(x14*(-Kx1*W0y + Ky1*W0x) + x39*x56);
  const double dydyt0 = x15*(Kz1*x58 + x17*x56 - x59*(Kx1*W0x + Ky1*W0y));
  AlgebraicMatrix55 res;
  res(0,0) = dqopdqop0;
  res(0,1) = dqopdlam0;
  res(0,2) = dqopdphi0;
  res(0,3) = dqopdxt0;
  res(0,4) = dqopdyt0;
  res(1,0) = ddxdzdqop0;
  res(1,1) = ddxdzdlam0;
  res(1,2) = ddxdzdphi0;
  res(1,3) = ddxdzdxt0;
  res(1,4) = ddxdzdyt0;
  res(2,0) = ddydzdqop0;
  res(2,1) = ddydzdlam0;
  res(2,2) = ddydzdphi0;
  res(2,3) = ddydzdxt0;
  res(2,4) = ddydzdyt0;
  res(3,0) = dxdqop0;
  res(3,1) = dxdlam0;
  res(3,2) = dxdphi0;
  res(3,3) = dxdxt0;
  res(3,4) = dxdyt0;
  res(4,0) = dydqop0;
  res(4,1) = dydlam0;
  res(4,2) = dydphi0;
  res(4,3) = dydxt0;
  res(4,4) = dydyt0;

  
  return res;
                                              
}

Eigen::Matrix<double, 5, 7> Geant4ePropagator::transportJacobian(const FreeTrajectoryState &start, double s, double dEdx, double mass) const {
  
  if (s==0.) {
    Eigen::Matrix<double, 5, 7> res;
    res.leftCols<5>() = Eigen::Matrix<double, 5, 5>::Identity();
    res.rightCols<2>() = Eigen::Matrix<double, 5, 2>::Zero();
    return res;
  }
  
  const GlobalTrajectoryParameters &globalSource = start.parameters();
  const GlobalVector& bfield = start.parameters().magneticFieldInInverseGeV();

  const double M0x = globalSource.position().x();
  const double M0y = globalSource.position().y();
  const double M0z = globalSource.position().z();
  
  const double p0 = globalSource.momentum().mag();
  const GlobalVector W0 = globalSource.momentum()/p0;
  const double W0x = W0.x();
  const double W0y = W0.y();
  const double W0z = W0.z();
  
  const double B = bfield.mag();
  const GlobalVector H = bfield/B;
  const double hx = H.x();
  const double hy = H.y();
  const double hz = H.z();
  
  const double qop0 = globalSource.signedInverseMomentum();
  const double q = start.charge();
  
      
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
  const double x10 = B*s;
  const double x11 = qop0*x10;
  const double x12 = std::cos(x11);
  const double x13 = std::pow(W0z, 2);
  const double x14 = std::pow(W0x, 2);
  const double x15 = std::pow(W0y, 2);
  const double x16 = x14 + x15;
  const double x17 = 1.0/x16;
  const double x18 = std::pow(x13*x17 + 1, -1.0/2.0);
  const double x19 = x12*x18;
  const double x20 = std::sqrt(x16);
  const double x21 = 1.0/x20;
  const double x22 = W0z*x21;
  const double x23 = x19*x22;
  const double x24 = W0y*hx;
  const double x25 = x21*x24;
  const double x26 = W0x*hy;
  const double x27 = x21*x26;
  const double x28 = x25 - x27;
  const double x29 = std::sin(x11);
  const double x30 = x18*x29;
  const double x31 = x28*x30;
  const double x32 = x12 - 1;
  const double x33 = hx*x18;
  const double x34 = W0x*x21;
  const double x35 = hy*x18;
  const double x36 = W0y*x21;
  const double x37 = hz*x18;
  const double x38 = x22*x37 + x33*x34 + x35*x36;
  const double x39 = hz*x38;
  const double x40 = x23 - x31 - x32*x39;
  const double x41 = x22*x30;
  const double x42 = x18*x32;
  const double x43 = x28*x42;
  const double x44 = B*qop0;
  const double x45 = W0z*x34;
  const double x46 = W0z*x36;
  const double x47 = -M0x*x45 - M0y*x46 + M0z*(x14*x21 + x15*x21);
  const double x48 = x20*x47;
  const double x49 = x11 - x29;
  const double x50 = x39*x49 + x41 + x43 + x44*x48;
  const double x51 = 1.0/B;
  const double x52 = x2*x51;
  const double x53 = x50*x52;
  const double x54 = x10*x12;
  const double x55 = x10 - x54;
  const double x56 = 1.0/qop0;
  const double x57 = x51*x56;
  const double x58 = x57*(B*x48 + x10*x23 - x10*x31 + x39*x55);
  const double x59 = x19*x36;
  const double x60 = x34*x37;
  const double x61 = x22*x33 - x60;
  const double x62 = x29*x61;
  const double x63 = hy*x38;
  const double x64 = x32*x63;
  const double x65 = x59 + x62 - x64;
  const double x66 = -M0x*x36 + M0y*x34;
  const double x67 = W0z*x47;
  const double x68 = W0x*x66 - W0y*x67;
  const double x69 = B*x68;
  const double x70 = x30*x36;
  const double x71 = qop0*x69 + x20*(-x32*x61 + x49*x63 + x70);
  const double x72 = x21*x52;
  const double x73 = x71*x72;
  const double x74 = x20*(x10*x59 + x10*x62 + x55*x63) + x69;
  const double x75 = x21*x57;
  const double x76 = x74*x75;
  const double x77 = x19*x34;
  const double x78 = x22*x35;
  const double x79 = x36*x37;
  const double x80 = x78 - x79;
  const double x81 = x29*x80;
  const double x82 = hx*x38;
  const double x83 = -x32*x82 + x77 - x81;
  const double x84 = W0x*x67 + W0y*x66;
  const double x85 = B*x84;
  const double x86 = x30*x34;
  const double x87 = -qop0*x85 + x20*(x32*x80 + x49*x82 + x86);
  const double x88 = x72*x87;
  const double x89 = x20*(x10*x77 - x10*x81 + x55*x82) - x85;
  const double x90 = x75*x89;
  const double x91 = -x40*(-x53 + x58) - x65*(-x73 + x76) - x83*(-x88 + x90);
  const double x92 = x9*(-s*x5 - x7);
  const double x93 = W0z*x17;
  const double x94 = W0x*x93;
  const double x95 = W0y*x93;
  const double x96 = -x33*x94 - x35*x95 + x37;
  const double x97 = x49*x96;
  const double x98 = hz*x97 - x22*x43 + x30;
  const double x99 = x40*x57;
  const double x100 = x37*x95;
  const double x101 = hx*x97 - x30*x94 + x32*(x100 + x35);
  const double x102 = x57*x83;
  const double x103 = 1 - x12;
  const double x104 = x33 + x37*x94;
  const double x105 = hy*x97 + x103*x104 - x30*x95;
  const double x106 = x57*x65;
  const double x107 = -x101*x102 - x105*x106 - x98*x99;
  const double x108 = W0x*hx;
  const double x109 = x108*x21;
  const double x110 = W0y*hy;
  const double x111 = x110*x21;
  const double x112 = -x33*x36 + x34*x35;
  const double x113 = x112*x49;
  const double x114 = hz*x113 + x42*(x109 + x111);
  const double x115 = hy*x113 + x103*x79 + x86;
  const double x116 = hx*x113 - x32*x60 - x70;
  const double x117 = -x102*x116 - x106*x115 - x114*x99;
  const double x118 = -x34*x65 + x36*x83;
  const double x119 = -x20*x40 + x45*x83 + x46*x65;
  const double x120 = x56/std::pow(B, 2);
  const double x121 = x120*x50;
  const double x122 = qop0*s;
  const double x123 = x12*x122;
  const double x124 = x122 - x123;
  const double x125 = x57*(qop0*x48 + x122*x23 - x122*x31 + x124*x39);
  const double x126 = x120*x21;
  const double x127 = x126*x71;
  const double x128 = qop0*x68 + x20*(x122*x59 + x122*x62 + x124*x63);
  const double x129 = x128*x75;
  const double x130 = x126*x87;
  const double x131 = -qop0*x84 + x20*(x122*x77 - x122*x81 + x124*x82);
  const double x132 = x131*x75;
  const double x133 = -x40*(-x121 + x125) - x65*(-x127 + x129) - x83*(-x130 + x132);
  const double x134 = std::pow(x65, 2) + std::pow(x83, 2);
  const double x135 = 1.0/x134;
  const double x136 = 1.0/(x135*std::pow(x40, 2) + 1);
  const double x137 = -x25 + x27;
  const double x138 = x137*x19;
  const double x139 = x29*x39;
  const double x140 = std::pow(x134, -1.0/2.0);
  const double x141 = x10*x70;
  const double x142 = x54*x61;
  const double x143 = x10*x29;
  const double x144 = x143*x63;
  const double x145 = (1.0/2.0)*x65;
  const double x146 = x10*x86;
  const double x147 = -x78 + x79;
  const double x148 = x147*x54;
  const double x149 = x143*x82;
  const double x150 = (1.0/2.0)*x83;
  const double x151 = x40/std::pow(x134, 3.0/2.0);
  const double x152 = x44*x70;
  const double x153 = x12*x44;
  const double x154 = x153*x61;
  const double x155 = x29*x44;
  const double x156 = x155*x63;
  const double x157 = x44*x86;
  const double x158 = x147*x153;
  const double x159 = x155*x82;
  const double x160 = x136*(x140*(x138*x44 + x139*x44 - x41*x44) + x151*(-x145*(-2*x152 + 2*x154 + 2*x156) - x150*(-2*x157 + 2*x158 + 2*x159)));
  const double x161 = hz*x32;
  const double x162 = x19*x95;
  const double x163 = x104*x29;
  const double x164 = x32*x96;
  const double x165 = hy*x164;
  const double x166 = x19*x94;
  const double x167 = x29*(-x100 - x35);
  const double x168 = hx*x164;
  const double x169 = W0y*hz;
  const double x170 = x21*x30;
  const double x171 = x169*x170;
  const double x172 = x112*x32;
  const double x173 = hy*x172;
  const double x174 = W0x*hz;
  const double x175 = x170*x174;
  const double x176 = hx*x172;
  const double x177 = x122*x70;
  const double x178 = x123*x61;
  const double x179 = x122*x29;
  const double x180 = x179*x63;
  const double x181 = x122*x86;
  const double x182 = x123*x147;
  const double x183 = x179*x82;
  const double x184 = x135*x83;
  const double x185 = -x59;
  const double x186 = x135*(x185 - x62 + x64);
  const double x187 = x184*(-x152 + x154 + x156) + x186*(-x157 + x158 + x159);
  const double x188 = W0z*hz + x108 + x110;
  const double x189 = x188*x32;
  const double x190 = W0y*x12 - hy*x189 + x29*(W0z*hx - x174);
  const double x191 = x13 + x16;
  const double x192 = 1.0/x191;
  const double x193 = W0x*x12 - hx*x189 + x29*(-W0z*hy + x169);
  const double x194 = x192*std::pow(x193, 2);
  const double x195 = std::pow(x190, 2)*x192;
  const double x196 = std::pow(x194 + x195, -1.0/2.0);
  const double x197 = x196/std::sqrt(x17*x191);
  const double x198 = x17*x197;
  const double x199 = x190*x198;
  const double x200 = x199*x87;
  const double x201 = x193*x198;
  const double x202 = x201*x71;
  const double x203 = x201*x57;
  const double x204 = x199*x57;
  const double x205 = -x153 + x44;
  const double x206 = x205*x63 + x44*x59 + x44*x62;
  const double x207 = x193*x197;
  const double x208 = x207*x75;
  const double x209 = x205*x82 + x44*x77 - x44*x81;
  const double x210 = x190*x197;
  const double x211 = x210*x75;
  const double x212 = x206*x208 - x209*x211;
  const double x213 = x192*x196*(W0z*x12 - x161*x188 + x29*(-x24 + x26));
  const double x214 = x193*x213;
  const double x215 = x190*x213;
  const double x216 = x194*x196 + x195*x196;
  const double x217 = x215*x57;
  const double x218 = x214*x57;
  const double x219 = x216*x57;
  const double x220 = -x206*x217 - x209*x218 + x219*(x205*x39 + x23*x44 - x31*x44);
  const double dqopdqop0 = x9*(dEdx*s*x1/x4 + x1) + x91*x92;
  const double dqopdlam0 = x107*x92;
  const double dqopdphi0 = x117*x92;
  const double dqopdxt0 = x118*x92;
  const double dqopdyt0 = x119*x92;
  const double dqopdB = x133*x92;
  const double dqopdxi = x9*(-x6 - x8);
  const double dlamdqop0 = x136*(x140*(x10*x138 + x10*x139 - x10*x41) + x151*(-x145*(-2*x141 + 2*x142 + 2*x144) - x150*(-2*x146 + 2*x148 + 2*x149))) + x160*x91;
  const double dlamdlam0 = x107*x160 + x136*(x140*(-x137*x41 - x161*x96 + x19) + x151*(-x145*(-2*x162 + 2*x163 - 2*x165) - x150*(-2*x166 + 2*x167 - 2*x168)));
  const double dlamdphi0 = x117*x160 + x136*(x140*(-x112*x161 + x30*(-x109 - x111)) + x151*(-x145*(2*x171 - 2*x173 + 2*x77) - x150*(2*x175 - 2*x176 - 2*x59)));
  const double dlamdxt0 = x118*x160;
  const double dlamdyt0 = x119*x160;
  const double dlamdB = x133*x160 + x136*(x140*(x122*x138 + x122*x139 - x122*x41) + x151*(-x145*(-2*x177 + 2*x178 + 2*x180) - x150*(-2*x181 + 2*x182 + 2*x183)));
  const double dlamdxi = 0;
  const double dphidqop0 = x184*(-x141 + x142 + x144) + x186*(-x146 + x148 + x149) + x187*x91;
  const double dphidlam0 = x107*x187 + x184*(-x162 + x163 - x165) + x186*(-x166 + x167 - x168);
  const double dphidphi0 = x117*x187 + x184*(x171 - x173 + x77) + x186*(x175 - x176 + x185);
  const double dphidxt0 = x118*x187;
  const double dphidyt0 = x119*x187;
  const double dphidB = x133*x187 + x184*(-x177 + x178 + x180) + x186*(-x181 + x182 + x183);
  const double dphidxi = 0;
  const double dxtdqop0 = x200*x52 - x202*x52 + x203*x74 - x204*x89 + x212*x91;
  const double dxtdlam0 = -x101*x211 + x105*x208 + x107*x212;
  const double dxtdphi0 = x115*x208 - x116*x211 + x117*x212;
  const double dxtdxt0 = W0x*x201 + W0y*x199 + x118*x212;
  const double dxtdyt0 = x119*x212 - x207*x95 + x210*x94;
  const double dxtdB = x120*x200 - x120*x202 + x128*x203 - x131*x204 + x133*x212;
  const double dxtdxi = 0;
  const double dytdqop0 = x214*x88 - x214*x90 + x215*x73 - x215*x76 - x216*x53 + x216*x58 + x220*x91;
  const double dytdlam0 = -x101*x218 - x105*x217 + x107*x220 + x219*x98;
  const double dytdphi0 = x114*x219 - x115*x217 - x116*x218 + x117*x220;
  const double dytdxt0 = x118*x220 + x214*x36 - x215*x34;
  const double dytdyt0 = x119*x220 + x20*x216 + x214*x45 + x215*x46;
  const double dytdB = -x121*x216 + x125*x216 + x127*x215 - x129*x215 + x130*x214 - x132*x214 + x133*x220;
  const double dytdxi = 0;
  Eigen::Matrix<double, 5, 7> res;
  res(0,0) = dqopdqop0;
  res(0,1) = dqopdlam0;
  res(0,2) = dqopdphi0;
  res(0,3) = dqopdxt0;
  res(0,4) = dqopdyt0;
  res(0,5) = dqopdB;
  res(0,6) = dqopdxi;
  res(1,0) = dlamdqop0;
  res(1,1) = dlamdlam0;
  res(1,2) = dlamdphi0;
  res(1,3) = dlamdxt0;
  res(1,4) = dlamdyt0;
  res(1,5) = dlamdB;
  res(1,6) = dlamdxi;
  res(2,0) = dphidqop0;
  res(2,1) = dphidlam0;
  res(2,2) = dphidphi0;
  res(2,3) = dphidxt0;
  res(2,4) = dphidyt0;
  res(2,5) = dphidB;
  res(2,6) = dphidxi;
  res(3,0) = dxtdqop0;
  res(3,1) = dxtdlam0;
  res(3,2) = dxtdphi0;
  res(3,3) = dxtdxt0;
  res(3,4) = dxtdyt0;
  res(3,5) = dxtdB;
  res(3,6) = dxtdxi;
  res(4,0) = dytdqop0;
  res(4,1) = dytdlam0;
  res(4,2) = dytdphi0;
  res(4,3) = dytdxt0;
  res(4,4) = dytdyt0;
  res(4,5) = dytdB;
  res(4,6) = dytdxi;

  res.col(5) *= 2.99792458e-3;
  
  return res;

  
}


Eigen::Matrix<double, 5, 7> Geant4ePropagator::transportJacobianBz(const FreeTrajectoryState &start, double s, double dEdx, double mass) const {
  
  if (s==0.) {
    Eigen::Matrix<double, 5, 7> res;
    res.leftCols<5>() = Eigen::Matrix<double, 5, 5>::Identity();
    res.rightCols<2>() = Eigen::Matrix<double, 5, 2>::Zero();
    return res;
  }
  
  const GlobalTrajectoryParameters &globalSource = start.parameters();
  const GlobalVector& bfield = start.parameters().magneticFieldInInverseGeV();

  const double M0x = globalSource.position().x();
  const double M0y = globalSource.position().y();
  const double M0z = globalSource.position().z();
  
  const double p0 = globalSource.momentum().mag();
  const GlobalVector W0 = globalSource.momentum()/p0;
  const double W0x = W0.x();
  const double W0y = W0.y();
  const double W0z = W0.z();
  
  const double Bx = bfield.x();
  const double By = bfield.y();
  const double Bz = bfield.z();
  
  
  const double qop0 = globalSource.signedInverseMomentum();
  const double q = start.charge();
              
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


Eigen::Matrix<double, 5, 7> Geant4ePropagator::transportJacobianBzAdvanced(const FreeTrajectoryState &start, double s, double dEdx, double mass) const {
  
  using namespace Eigen;
  
  if (s==0.) {
    Eigen::Matrix<double, 5, 7> res;
    res.leftCols<5>() = Eigen::Matrix<double, 5, 5>::Identity();
    res.rightCols<2>() = Eigen::Matrix<double, 5, 2>::Zero();
    return res;
  }
  
  const GlobalTrajectoryParameters &globalSource = start.parameters();
  const GlobalVector& bfield = start.parameters().magneticFieldInInverseGeV();

  const double M0x = globalSource.position().x();
  const double M0y = globalSource.position().y();
  const double M0z = globalSource.position().z();
  
  const double p0 = globalSource.momentum().mag();
  const GlobalVector W0 = globalSource.momentum()/p0;
  const double W0x = W0.x();
  const double W0y = W0.y();
  const double W0z = W0.z();
  
  const double Bx = bfield.x();
  const double By = bfield.y();
  const double Bz = bfield.z();
  
  //compute gradient of magnetic field perpendicular to track
  
  const Matrix<double, 6, 5> curv2cartjac = curv2cartJacobianAlt(start);
  
  const double dx = 10e-4;
  
  Matrix<double, 5, 1> dxt0 = Matrix<double, 5, 1>::Zero();
  dxt0[3] = dx;
  
  Matrix<double, 5, 1> dyt0 = Matrix<double, 5, 1>::Zero();
  dyt0[4] = dx;
  
  const Vector3d dxt0cart = (curv2cartjac*dxt0).head<3>();
  const Vector3d dyt0cart = (curv2cartjac*dyt0).head<3>();
  
  const GlobalPoint xt0pos(M0x + dxt0cart[0], M0y + dxt0cart[1], M0z + dxt0cart[2]);
  const GlobalPoint yt0pos(M0x + dyt0cart[0], M0y + dyt0cart[1], M0z + dyt0cart[2]);
  
  
  
  const GlobalVector Bxt0 = start.parameters().magneticField().inInverseGeV(xt0pos);
  const GlobalVector Byt0 = start.parameters().magneticField().inInverseGeV(yt0pos);
  
  const GlobalVector dBdxt0 = (bfield - Bxt0)/dx;
  const GlobalVector dBdyt0 = (bfield - Byt0)/dx;
  
  std::cout << "dBdxt0" << std::endl;
  std::cout << dBdxt0 << std::endl;
  std::cout << "dBdyt0" << std::endl;
  std::cout << dBdyt0 << std::endl;
  
  const GlobalPoint spos(M0x + W0x*dx, M0y + W0y*dx, M0z + W0z*dx);
  const GlobalVector Bs = start.parameters().magneticField().inInverseGeV(spos);
  const GlobalVector dBds = (bfield - Bs)/dx;
  std::cout << "dBds" << std::endl;
  std::cout << dBds << std::endl;
  
  const double dBdxt0x = dBdxt0.x();
  const double dBdxt0y = dBdxt0.y();
  const double dBdxt0z = dBdxt0.z();
  
  const double dBdyt0x = dBdyt0.x();
  const double dBdyt0y = dBdyt0.y();
  const double dBdyt0z = dBdyt0.z();
  
  const double qop0 = globalSource.signedInverseMomentum();
  const double q = start.charge();
            
  const double x0 = std::pow(q, 2);
  const double x1 = x0/std::pow(qop0, 3);
  const double x2 = std::pow(qop0, 2);
  const double x3 = 1.0/x2;
  const double x4 = x0*x3;
  const double x5 = std::sqrt(std::pow(mass, 2) + x4);
  const double x6 = std::pow(dEdx, 2);
  const double x7 = std::pow(s, 2)*x6;
  const double x8 = dEdx*x5;
  const double x9 = s*x8;
  const double x10 = q/std::pow(x4 + x7 + 2*x9, 3.0/2.0);
  const double x11 = std::pow(W0x, 2);
  const double x12 = std::pow(W0y, 2);
  const double x13 = x11 + x12;
  const double x14 = std::pow(x13, -1.0/2.0);
  const double x15 = W0x*x14;
  const double x16 = std::pow(Bz, 2);
  const double x17 = std::pow(Bx, 2) + std::pow(By, 2) + x16;
  const double x18 = std::sqrt(x17);
  const double x19 = qop0*x18;
  const double x20 = s*x19;
  const double x21 = std::cos(x20);
  const double x22 = std::pow(W0z, 2);
  const double x23 = 1.0/x13;
  const double x24 = std::pow(x22*x23 + 1, -1.0/2.0);
  const double x25 = x21*x24;
  const double x26 = x15*x25;
  const double x27 = std::sin(x20);
  const double x28 = By*W0z;
  const double x29 = 1.0/x18;
  const double x30 = x24*x29;
  const double x31 = x14*x30;
  const double x32 = x28*x31;
  const double x33 = W0y*x14;
  const double x34 = Bz*x30;
  const double x35 = x33*x34;
  const double x36 = x32 - x35;
  const double x37 = x27*x36;
  const double x38 = x19*x21;
  const double x39 = -x19 + x38;
  const double x40 = 1.0/qop0;
  const double x41 = 1.0/x17;
  const double x42 = Bz*W0z;
  const double x43 = Bx*x30;
  const double x44 = x15*x43;
  const double x45 = By*x30;
  const double x46 = x33*x45;
  const double x47 = x44 + x46;
  const double x48 = x31*x42 + x47;
  const double x49 = x41*x48;
  const double x50 = x40*x49;
  const double x51 = x39*x50;
  const double x52 = -Bx*x51 + x26 - x37;
  const double x53 = s*x40;
  const double x54 = x27*x30;
  const double x55 = x3*x54;
  const double x56 = 1 - x21;
  const double x57 = x29*x56;
  const double x58 = x3*x57;
  const double x59 = -x20 + x27;
  const double x60 = x3*x49;
  const double x61 = x59*x60;
  const double x62 = s*x18;
  const double x63 = x21*x62;
  const double x64 = x50*(-x62 + x63);
  const double x65 = Bx*x61 - Bx*x64 - x15*x55 + x26*x53 + x36*x58 - x37*x53;
  const double x66 = x25*x33;
  const double x67 = x15*x34;
  const double x68 = Bx*W0z;
  const double x69 = x31*x68;
  const double x70 = x67 - x69;
  const double x71 = x27*x70;
  const double x72 = By*x51;
  const double x73 = x66 - x71 - x72;
  const double x74 = By*x61 - By*x64 - x33*x55 + x53*x66 - x53*x71 + x58*x70;
  const double x75 = W0z*x14;
  const double x76 = x25*x75;
  const double x77 = x33*x43;
  const double x78 = x15*x45;
  const double x79 = x77 - x78;
  const double x80 = x27*x79;
  const double x81 = -Bz*x51 + x76 - x80;
  const double x82 = Bz*x59;
  const double x83 = -Bz*x64 + x53*x76 - x53*x80 - x55*x75 + x58*x79 + x60*x82;
  const double x84 = -x52*x65 - x73*x74 - x81*x83;
  const double x85 = x10*(-s*x6 - x8);
  const double x86 = W0x*x23;
  const double x87 = x40*x54;
  const double x88 = W0z*x87;
  const double x89 = Bz*W0y;
  const double x90 = W0z*x23*x30;
  const double x91 = x89*x90;
  const double x92 = x40*x57;
  const double x93 = Bx*W0x;
  const double x94 = W0y*x23;
  const double x95 = x28*x30;
  const double x96 = x40*x41;
  const double x97 = x96*(x34 - x90*x93 - x94*x95);
  const double x98 = x59*x97;
  const double x99 = -Bx*x98 - x86*x88 - x92*(x45 + x91);
  const double x100 = x86*x95;
  const double x101 = x30*x68*x94;
  const double x102 = -x82*x97 + x87 - x92*(x100 - x101);
  const double x103 = x30*x42*x86;
  const double x104 = -By*x98 - x88*x94 - x92*(-x103 - x43);
  const double x105 = -x102*x81 - x104*x73 - x52*x99;
  const double x106 = -x77 + x78;
  const double x107 = x106*x96;
  const double x108 = -x107*x82 - x47*x92;
  const double x109 = Bz*x24*x56*x96;
  const double x110 = x107*x59;
  const double x111 = -By*x110 + x109*x33 + x15*x87;
  const double x112 = -Bx*x110 + x109*x15 - x33*x87;
  const double x113 = -x108*x81 - x111*x73 - x112*x52;
  const double x114 = Bx*dBdxt0x;
  const double x115 = By*dBdxt0y;
  const double x116 = Bz*dBdxt0z;
  const double x117 = x114 + x115 + x116;
  const double x118 = s*x117;
  const double x119 = x118*x41;
  const double x120 = std::pow(x17, 3.0/2.0);
  const double x121 = 1.0/x120;
  const double x122 = x121*(-x114 - x115 - x116);
  const double x123 = x122*x24;
  const double x124 = x27*x40;
  const double x125 = x124*x75;
  const double x126 = x40*x56;
  const double x127 = x122*x126;
  const double x128 = x50*x59;
  const double x129 = std::pow(x17, 2);
  const double x130 = 1.0/x129;
  const double x131 = x130*x40*x48;
  const double x132 = x131*(-2*x114 - 2*x115 - 2*x116);
  const double x133 = qop0*x29;
  const double x134 = x117*x133;
  const double x135 = s*x134;
  const double x136 = x135*x21;
  const double x137 = x50*(-x135 + x136);
  const double x138 = x30*x33;
  const double x139 = dBdxt0x*x138;
  const double x140 = x15*x30;
  const double x141 = dBdxt0y*x140;
  const double x142 = Bx*x33;
  const double x143 = x123*x142;
  const double x144 = x123*x15;
  const double x145 = By*x144;
  const double x146 = x30*x75;
  const double x147 = x123*x33;
  const double x148 = x123*x14;
  const double x149 = x96*(Bx*x144 + By*x147 + dBdxt0x*x140 + dBdxt0y*x138 + dBdxt0z*x146 + x148*x42);
  const double x150 = -Bz*x137 - dBdxt0z*x128 + x119*x76 - x119*x80 + x123*x125 - x127*x79 - x132*x82 - x149*x82 - x92*(x139 - x141 + x143 - x145);
  const double x151 = x132*x59;
  const double x152 = dBdxt0z*x140;
  const double x153 = dBdxt0x*x146;
  const double x154 = Bz*x144;
  const double x155 = x148*x68;
  const double x156 = x149*x59;
  const double x157 = -By*x137 - By*x151 - By*x156 - dBdxt0y*x128 + x119*x66 - x119*x71 + x124*x147 - x127*x70 + x15 - x92*(x152 - x153 + x154 - x155);
  const double x158 = dBdxt0y*x146;
  const double x159 = dBdxt0z*x138;
  const double x160 = x148*x28;
  const double x161 = Bz*x147;
  const double x162 = -Bx*x137 - Bx*x151 - Bx*x156 - dBdxt0x*x128 + x119*x26 - x119*x37 + x124*x144 - x127*x36 - x33 - x92*(x158 - x159 + x160 - x161);
  const double x163 = -x150*x81 - x157*x73 - x162*x52;
  const double x164 = Bx*dBdyt0x;
  const double x165 = By*dBdyt0y;
  const double x166 = Bz*dBdyt0z;
  const double x167 = x164 + x165 + x166;
  const double x168 = s*x41;
  const double x169 = x167*x168;
  const double x170 = -x164 - x165 - x166;
  const double x171 = x121*x24;
  const double x172 = x170*x171;
  const double x173 = x15*x172;
  const double x174 = x121*x126;
  const double x175 = x170*x174;
  const double x176 = x131*(-2*x164 - 2*x165 - 2*x166);
  const double x177 = x176*x59;
  const double x178 = x133*x167;
  const double x179 = s*x178;
  const double x180 = x179*x21;
  const double x181 = x50*(-x179 + x180);
  const double x182 = dBdyt0y*x146;
  const double x183 = dBdyt0z*x138;
  const double x184 = x14*x172;
  const double x185 = x184*x28;
  const double x186 = x172*x33;
  const double x187 = Bz*x186;
  const double x188 = x96*(Bx*x173 + By*x186 + dBdyt0x*x140 + dBdyt0y*x138 + dBdyt0z*x146 + x184*x42);
  const double x189 = x188*x59;
  const double x190 = -Bx*x177 - Bx*x181 - Bx*x189 - W0z*x15 - dBdyt0x*x128 + x124*x173 + x169*x26 - x169*x37 - x175*x36 - x92*(x182 - x183 + x185 - x187);
  const double x191 = dBdyt0z*x140;
  const double x192 = dBdyt0x*x146;
  const double x193 = Bz*x173;
  const double x194 = x184*x68;
  const double x195 = -By*x177 - By*x181 - By*x189 - W0z*x33 - dBdyt0y*x128 + x124*x186 + x169*x66 - x169*x71 - x175*x70 - x92*(x191 - x192 + x193 - x194);
  const double x196 = dBdyt0x*x138;
  const double x197 = dBdyt0y*x140;
  const double x198 = x142*x172;
  const double x199 = By*x173;
  const double x200 = -Bz*x181 - dBdyt0z*x128 + x11*x14 + x12*x14 + x125*x172 + x169*x76 - x169*x80 - x175*x79 - x176*x82 - x188*x82 - x92*(x196 - x197 + x198 - x199);
  const double x201 = -x190*x52 - x195*x73 - x200*x81;
  const double x202 = Bz*x168;
  const double x203 = Bz*x171;
  const double x204 = x203*x33;
  const double x205 = Bz*x174;
  const double x206 = x14*x42;
  const double x207 = x171*x206;
  const double x208 = Bx*x207;
  const double x209 = x16*x171;
  const double x210 = x15*x209;
  const double x211 = 2*x131;
  const double x212 = x211*x82;
  const double x213 = Bz*x133;
  const double x214 = s*x213;
  const double x215 = x21*x214;
  const double x216 = x50*(-x214 + x215);
  const double x217 = x15*x203;
  const double x218 = x96*(-Bx*x217 - By*x204 + x146 - x209*x75);
  const double x219 = x218*x59;
  const double x220 = By*x212 - By*x216 - By*x219 - x124*x204 + x202*x66 - x202*x71 + x205*x70 - x92*(x140 + x208 - x210);
  const double x221 = x209*x33;
  const double x222 = x14*x203*x28;
  const double x223 = Bx*x212 - Bx*x216 - Bx*x219 - x124*x217 + x202*x26 - x202*x37 + x205*x36 - x92*(-x138 + x221 - x222);
  const double x224 = By*x217;
  const double x225 = x142*x203;
  const double x226 = x16*x211;
  const double x227 = -Bz*x216 - x124*x207 - x128 + x168*x206*x25 - x202*x80 + x205*x79 - x218*x82 + x226*x59 - x92*(x224 - x225);
  const double x228 = -x220*x73 - x223*x52 - x227*x81;
  const double x229 = std::pow(x52, 2) + std::pow(x73, 2);
  const double x230 = 1.0/x229;
  const double x231 = 1.0/(x230*std::pow(x81, 2) + 1);
  const double x232 = x24*x27;
  const double x233 = x232*x75;
  const double x234 = Bz*x39;
  const double x235 = qop0*x27;
  const double x236 = s*x235;
  const double x237 = -x17*x236 + x18*x21 - x18;
  const double x238 = Bz*x50;
  const double x239 = std::pow(x229, -1.0/2.0);
  const double x240 = x232*x62;
  const double x241 = x15*x240;
  const double x242 = -x32 + x35;
  const double x243 = x242*x63;
  const double x244 = x39*x60;
  const double x245 = Bx*x244;
  const double x246 = x237*x50;
  const double x247 = Bx*x246;
  const double x248 = (1.0/2.0)*x52;
  const double x249 = x240*x33;
  const double x250 = -x67 + x69;
  const double x251 = x250*x63;
  const double x252 = By*x244;
  const double x253 = By*x246;
  const double x254 = (1.0/2.0)*x73;
  const double x255 = x81/std::pow(x229, 3.0/2.0);
  const double x256 = x235*x48;
  const double x257 = x19*x232;
  const double x258 = x15*x257;
  const double x259 = x242*x38;
  const double x260 = Bx*x256;
  const double x261 = x257*x33;
  const double x262 = x250*x38;
  const double x263 = By*x256;
  const double x264 = x231*(x239*(Bz*x256 + x106*x38 - x19*x233) + x255*(-x248*(-2*x258 + 2*x259 + 2*x260) - x254*(-2*x261 + 2*x262 + 2*x263)));
  const double x265 = W0z*x25;
  const double x266 = x265*x94;
  const double x267 = x27*(x103 + x43);
  const double x268 = x39*x97;
  const double x269 = By*x268;
  const double x270 = x265*x86;
  const double x271 = x27*(-x45 - x91);
  const double x272 = Bx*x268;
  const double x273 = x27*x67;
  const double x274 = x107*x39;
  const double x275 = Bx*x274;
  const double x276 = x27*x35;
  const double x277 = By*x274;
  const double x278 = qop0*x54;
  const double x279 = x278*x75;
  const double x280 = x2*x27;
  const double x281 = -x118*x280 + x134*x21 - x134;
  const double x282 = x118*x278;
  const double x283 = x15*x282;
  const double x284 = x136*x242;
  const double x285 = dBdxt0x*x51;
  const double x286 = x132*x39;
  const double x287 = Bx*x286;
  const double x288 = x27*(-x158 + x159 - x160 + x161);
  const double x289 = x281*x50;
  const double x290 = Bx*x289;
  const double x291 = x149*x39;
  const double x292 = Bx*x291;
  const double x293 = x282*x33;
  const double x294 = x136*x250;
  const double x295 = dBdxt0y*x51;
  const double x296 = By*x286;
  const double x297 = x27*(-x152 + x153 - x154 + x155);
  const double x298 = By*x289;
  const double x299 = By*x291;
  const double x300 = s*x167;
  const double x301 = s*x280;
  const double x302 = -x167*x301 + x178*x21 - x178;
  const double x303 = x278*x300;
  const double x304 = x15*x303;
  const double x305 = x180*x242;
  const double x306 = dBdyt0x*x51;
  const double x307 = x176*x39;
  const double x308 = Bx*x307;
  const double x309 = x27*(-x182 + x183 - x185 + x187);
  const double x310 = x302*x50;
  const double x311 = Bx*x310;
  const double x312 = x188*x39;
  const double x313 = Bx*x312;
  const double x314 = x303*x33;
  const double x315 = x180*x250;
  const double x316 = dBdyt0y*x51;
  const double x317 = By*x307;
  const double x318 = x27*(-x191 + x192 - x193 + x194);
  const double x319 = By*x310;
  const double x320 = By*x312;
  const double x321 = -Bz*x301 + x21*x213 - x213;
  const double x322 = x236*x67;
  const double x323 = x215*x242;
  const double x324 = x27*(x138 - x221 + x222);
  const double x325 = Bx*x234;
  const double x326 = 4*x131;
  const double x327 = x321*x50;
  const double x328 = Bx*x327;
  const double x329 = x218*x39;
  const double x330 = Bx*x329;
  const double x331 = x236*x35;
  const double x332 = x215*x250;
  const double x333 = x27*(-x140 - x208 + x210);
  const double x334 = By*x234;
  const double x335 = By*x327;
  const double x336 = By*x329;
  const double x337 = -x66;
  const double x338 = x230*(x337 + x71 + x72);
  const double x339 = x230*x52;
  const double x340 = x338*(-x258 + x259 + x260) + x339*(-x261 + x262 + x263);
  const double x341 = x129*x21;
  const double x342 = x120*x27;
  const double x343 = x17*(x21 - 1)*(By*W0y + x42 + x93);
  const double x344 = -By*x343 + W0y*x341 + x342*(-Bz*W0x + x68);
  const double x345 = x13 + x22;
  const double x346 = 1/(std::pow(x17, 4)*x345);
  const double x347 = std::pow(x344, 2)*x346;
  const double x348 = -Bx*x343 + W0x*x341 - x342*(x28 - x89);
  const double x349 = x346*std::pow(x348, 2);
  const double x350 = std::pow(x347 + x349, -1.0/2.0);
  const double x351 = x348*x350;
  const double x352 = x130*x345/(std::pow(x13, 3.0/2.0)*std::pow(x23*x345, 3.0/2.0));
  const double x353 = x351*x352;
  const double x354 = x344*x350;
  const double x355 = x352*x354;
  const double x356 = x353*x73 - x355*x52;
  const double x357 = x346*(-Bz*x343 + W0z*x341 - x342*(Bx*W0y - By*W0x));
  const double x358 = x354*x357;
  const double x359 = x351*x357;
  const double x360 = x347*x350 + x349*x350;
  const double x361 = -x358*x73 - x359*x52 + x360*x81;
  const double dqopdqop0 = x10*(dEdx*s*x1/x5 + x1) + x84*x85;
  const double dqopdlam0 = x105*x85;
  const double dqopdphi0 = x113*x85;
  const double dqopdxt0 = x163*x85;
  const double dqopdyt0 = x201*x85;
  const double dqopdBz = x228*x85;
  const double dqopdxi = x10*(-x7 - x9);
  const double dlamdqop0 = x231*(x239*(x106*x63 - x233*x62 + x234*x60 - x237*x238) + x255*(-x248*(-2*x241 + 2*x243 + 2*x245 - 2*x247) - x254*(-2*x249 + 2*x251 + 2*x252 - 2*x253))) + x264*x84;
  const double dlamdlam0 = x105*x264 + x231*(x239*(-x234*x97 + x25 + x27*(-x100 + x101)) + x255*(-x248*(-2*x270 + 2*x271 - 2*x272) - x254*(-2*x266 + 2*x267 - 2*x269)));
  const double dlamdphi0 = x113*x264 + x231*(x239*(-x107*x234 + x27*(-x44 - x46)) + x255*(-x248*(2*x273 - 2*x275 - 2*x66) - x254*(2*x26 + 2*x276 - 2*x277)));
  const double dlamdxt0 = x163*x264 + x231*(x239*(-dBdxt0z*x51 + x106*x136 - x118*x279 - x132*x234 - x149*x234 - x238*x281 + x27*(-x139 + x141 - x143 + x145)) + x255*(-x248*(-2*x283 + 2*x284 - 2*x285 - 2*x287 + 2*x288 - 2*x290 - 2*x292) - x254*(-2*x293 + 2*x294 - 2*x295 - 2*x296 + 2*x297 - 2*x298 - 2*x299)));
  const double dlamdyt0 = x201*x264 + x231*(x239*(-dBdyt0z*x51 + x106*x180 - x176*x234 - x188*x234 - x238*x302 + x27*(-x196 + x197 - x198 + x199) - x279*x300) + x255*(-x248*(-2*x304 + 2*x305 - 2*x306 - 2*x308 + 2*x309 - 2*x311 - 2*x313) - x254*(-2*x314 + 2*x315 - 2*x316 - 2*x317 + 2*x318 - 2*x319 - 2*x320)));
  const double dlamdBz = x228*x264 + x231*(x239*(-s*x206*x278 + x106*x215 - x218*x234 + x226*x39 - x238*x321 + x27*(-x224 + x225) - x51) + x255*(-x248*(-2*x322 + 2*x323 + 2*x324 + x325*x326 - 2*x328 - 2*x330) - x254*(x326*x334 - 2*x331 + 2*x332 + 2*x333 - 2*x335 - 2*x336)));
  const double dlamdxi = 0;
  const double dphidqop0 = x338*(-x241 + x243 + x245 - x247) + x339*(-x249 + x251 + x252 - x253) + x340*x84;
  const double dphidlam0 = x105*x340 + x338*(-x270 + x271 - x272) + x339*(-x266 + x267 - x269);
  const double dphidphi0 = x113*x340 + x338*(x273 - x275 + x337) + x339*(x26 + x276 - x277);
  const double dphidxt0 = x163*x340 + x338*(-x283 + x284 - x285 - x287 + x288 - x290 - x292) + x339*(-x293 + x294 - x295 - x296 + x297 - x298 - x299);
  const double dphidyt0 = x201*x340 + x338*(-x304 + x305 - x306 - x308 + x309 - x311 - x313) + x339*(-x314 + x315 - x316 - x317 + x318 - x319 - x320);
  const double dphidBz = x228*x340 + x338*(x211*x325 - x322 + x323 + x324 - x328 - x330) + x339*(x211*x334 - x331 + x332 + x333 - x335 - x336);
  const double dphidxi = 0;
  const double dxtdqop0 = x353*x74 - x355*x65 + x356*x84;
  const double dxtdlam0 = x104*x353 + x105*x356 - x355*x99;
  const double dxtdphi0 = x111*x353 - x112*x355 + x113*x356;
  const double dxtdxt0 = x157*x353 - x162*x355 + x163*x356;
  const double dxtdyt0 = -x190*x355 + x195*x353 + x201*x356;
  const double dxtdBz = x220*x353 - x223*x355 + x228*x356;
  const double dxtdxi = 0;
  const double dytdqop0 = -x358*x74 - x359*x65 + x360*x83 + x361*x84;
  const double dytdlam0 = x102*x360 - x104*x358 + x105*x361 - x359*x99;
  const double dytdphi0 = x108*x360 - x111*x358 - x112*x359 + x113*x361;
  const double dytdxt0 = x150*x360 - x157*x358 - x162*x359 + x163*x361;
  const double dytdyt0 = -x190*x359 - x195*x358 + x200*x360 + x201*x361;
  const double dytdBz = -x220*x358 - x223*x359 + x227*x360 + x228*x361;
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


Eigen::Matrix<double, 6, 5> Geant4ePropagator::curv2cartJacobianAlt(const FreeTrajectoryState &state) const {
  
  const GlobalTrajectoryParameters &globalSource = state.parameters();
  
  CurvilinearTrajectoryParameters curvparms(globalSource.position(), globalSource.momentum(), globalSource.charge());

  const double qop0 = curvparms.Qbp();
  const double lam0 = curvparms.lambda();
  const double phi0 = curvparms.phi();
  const double xt0 = curvparms.xT();
  const double yt0 = curvparms.yT();
  
  const double p0 = globalSource.momentum().mag();
  const GlobalVector W0 = globalSource.momentum()/p0;
  const double W0x = W0.x();
  const double W0y = W0.y();
  const double W0z = W0.z();
    
  const double x0 = std::sqrt(std::pow(W0x, 2) + std::pow(W0y, 2));
  const double x1 = 1.0/x0;
  const double x2 = W0y*x1;
  const double x3 = W0x*x1;
  const double x4 = 1.0/qop0;
  const double x5 = std::cos(phi0);
  const double x6 = 1.0/std::fabs(qop0);
  const double x7 = x6*std::cos(lam0);
  const double x8 = x5*x7;
  const double x9 = x6*std::sin(lam0);
  const double x10 = std::sin(phi0);
  const double x11 = x10*x7;
  const double dxdqop0 = 0;
  const double dxdlam0 = 0;
  const double dxdphi0 = 0;
  const double dxdxt0 = -x2;
  const double dxdyt0 = -W0z*x3;
  const double dydqop0 = 0;
  const double dydlam0 = 0;
  const double dydphi0 = 0;
  const double dydxt0 = x3;
  const double dydyt0 = -W0z*x2;
  const double dzdqop0 = 0;
  const double dzdlam0 = 0;
  const double dzdphi0 = 0;
  const double dzdxt0 = 0;
  const double dzdyt0 = x0;
  const double dpxdqop0 = -x4*x8;
  const double dpxdlam0 = -x5*x9;
  const double dpxdphi0 = -x11;
  const double dpxdxt0 = 0;
  const double dpxdyt0 = 0;
  const double dpydqop0 = -x11*x4;
  const double dpydlam0 = -x10*x9;
  const double dpydphi0 = x8;
  const double dpydxt0 = 0;
  const double dpydyt0 = 0;
  const double dpzdqop0 = -x4*x9;
  const double dpzdlam0 = x7;
  const double dpzdphi0 = 0;
  const double dpzdxt0 = 0;
  const double dpzdyt0 = 0;
  Eigen::Matrix<double, 6, 5> res;
  res(0,0) = dxdqop0;
  res(0,1) = dxdlam0;
  res(0,2) = dxdphi0;
  res(0,3) = dxdxt0;
  res(0,4) = dxdyt0;
  res(1,0) = dydqop0;
  res(1,1) = dydlam0;
  res(1,2) = dydphi0;
  res(1,3) = dydxt0;
  res(1,4) = dydyt0;
  res(2,0) = dzdqop0;
  res(2,1) = dzdlam0;
  res(2,2) = dzdphi0;
  res(2,3) = dzdxt0;
  res(2,4) = dzdyt0;
  res(3,0) = dpxdqop0;
  res(3,1) = dpxdlam0;
  res(3,2) = dpxdphi0;
  res(3,3) = dpxdxt0;
  res(3,4) = dpxdyt0;
  res(4,0) = dpydqop0;
  res(4,1) = dpydlam0;
  res(4,2) = dpydphi0;
  res(4,3) = dpydxt0;
  res(4,4) = dpydyt0;
  res(5,0) = dpzdqop0;
  res(5,1) = dpzdlam0;
  res(5,2) = dpzdphi0;
  res(5,3) = dpzdxt0;
  res(5,4) = dpzdyt0;
  
  return res;

  
}
