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
// Class description:
//
// Description: Continuous Process to calculate energy loss without
//              fluctuations through G4EnergyLossForExtrapolator

// History:
//   Created:     2007-05-12
//   Author:      Pedro Arce
//
//   Modified:
//

#ifndef TrackPropagation_G4ErrorEnergyLossForCVH_h
#define TrackPropagation_G4ErrorEnergyLossForCVH_h

#include "globals.hh"
#include "G4VContinuousProcess.hh"
#include "G4ProcessType.hh"

class G4EnergyLossForExtrapolatorForCVH;

class G4ErrorEnergyLossForCVH : public G4VContinuousProcess {
public:
  explicit G4ErrorEnergyLossForCVH(const G4String& processName = "G4ErrorEnergyLossForCVH",
                                   G4ProcessType type = fElectromagnetic);

  ~G4ErrorEnergyLossForCVH() override;

  G4bool IsApplicable(const G4ParticleDefinition& aParticleType) override;
  // Returns true -> 'is applicable', for all charged particles.

  G4double GetContinuousStepLimit(const G4Track& aTrack, G4double, G4double currentMinimumStep, G4double&) override;

  G4VParticleChange* AlongStepDoIt(const G4Track& aTrack, const G4Step& aStep) override;
  // This is the method implementing the energy loss process.

  // Get and Set methods
  inline G4double GetStepLimit() const { return theStepLimit; }
  inline void SetStepLimit(G4double val) { theStepLimit = val; }

  // copy constructor and hide assignment operator
  G4ErrorEnergyLossForCVH(G4ErrorEnergyLossForCVH&) = delete;
  G4ErrorEnergyLossForCVH& operator=(const G4ErrorEnergyLossForCVH& right) = delete;

private:
  void InstantiateEforExtrapolator();

  G4EnergyLossForExtrapolatorForCVH* theELossForExtrapolator;

  G4double theStepLimit;
  G4double theFractionLimit = 0.2;
};

#endif
