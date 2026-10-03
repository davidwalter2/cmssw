//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//
//
// Class Description:
//
//  Default physics list for GEANT4e (should not be overridden, unless by
//  experts). No multiple scattering and no production of secondaries.
//  The energy loss process is G4eMuIonisation or G4EnergyLossForExtrapolator
//  (depending on the value of the enviromental variable G4EELOSSEXTRAP)
//  It also defines the geant4e processes to limit the step:
//  G4eMagneticFieldLimitProcess, G4eStepLimitProcess.

// History:
// - Created:   P. Arce
// ---------------------------------------------------------------------

#ifndef TrackPropagation_G4ErrorPhysicsListForCVH_h
#define TrackPropagation_G4ErrorPhysicsListForCVH_h

#include "globals.hh"
#include "G4VUserPhysicsList.hh"

#include <string>
#include <vector>

class G4ErrorPhysicsListForCVH : public G4VUserPhysicsList {
public:  // with description
  // Default ctor: registers the canonical full CVH particle set
  // (gamma, e+-, mu+-, pi+-, K+-, p, anti_proton). Kept for back-compat
  // with any direct caller; new code should pass the narrowed list it
  // actually needs through the vector ctor.
  G4ErrorPhysicsListForCVH();

  // Construct with an explicit list of G4 particle names (e.g.
  // {"mu+", "mu-"} for a J/psi-only refit). Recognised names:
  //   "gamma", "e+", "e-", "mu+", "mu-", "pi+", "pi-",
  //   "kaon+", "kaon-", "proton", "anti_proton".
  // An unknown name aborts construction with a G4Exception.
  explicit G4ErrorPhysicsListForCVH(const std::vector<std::string>& particleNames);

  ~G4ErrorPhysicsListForCVH() override;

protected:
  void ConstructParticle() override;
  // constructs the particles listed in particleNames_

  void ConstructProcess() override;
  // construct physical processes

  void SetCuts() override;
  // SetCutsWithDefault

  virtual void ConstructEM();
  // constructs electromagnetic processes

private:
  // G4 particle names to register in ConstructParticle. Captured at
  // construction; const after.
  const std::vector<std::string> particleNames_;
};

#endif
