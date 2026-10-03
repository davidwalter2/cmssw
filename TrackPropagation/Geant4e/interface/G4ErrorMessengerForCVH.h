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
// Class Description:
//
// Messenger class for Geant4e processes limiting the step.

// History:
// - Created:   P. Arce
// --------------------------------------------------------------------

#ifndef TrackPropagation_G4ErrorMessengerForCVH_h
#define TrackPropagation_G4ErrorMessengerForCVH_h

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4ErrorStepLengthLimitProcess;
class G4ErrorMagFieldLimitProcess;
class G4ErrorEnergyLossForCVH;

//-----------------------------------------------------------------

class G4ErrorMessengerForCVH : public G4UImessenger {
public:  // with description
  G4ErrorMessengerForCVH(G4ErrorStepLengthLimitProcess* lengthAct,
                         G4ErrorMagFieldLimitProcess* magAct,
                         G4ErrorEnergyLossForCVH* elossAct);
  ~G4ErrorMessengerForCVH() override;

  void SetNewValue(G4UIcommand*, G4String) override;

private:
  G4ErrorStepLengthLimitProcess* StepLengthAction;
  G4ErrorMagFieldLimitProcess* MagFieldAction;
  G4ErrorEnergyLossForCVH* EnergyLossAction;

  G4UIdirectory* myDir;
  G4UIdirectory* myDirLimits;

  G4UIcmdWithADoubleAndUnit* StepLengthLimitCmd;
  G4UIcmdWithADouble* MagFieldLimitCmd;
  G4UIcmdWithADouble* EnergyLossCmd;
};

#endif
