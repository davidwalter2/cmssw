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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4UniversalFluctuationForExtrapolator
//
// Author:        V.Ivanchenko make a class with the Laszlo Urban model
//
// Creation date: 03.01.2002
//
// Modifications:
//
//
// Class Description:
//
// Implementation of energy loss fluctuations made by L.Urban in 2021

// -------------------------------------------------------------------
//

#ifndef G4UniversalFluctuationForExtrapolator_h
#define G4UniversalFluctuationForExtrapolator_h 1

#include "G4VEmFluctuationModel.hh"
#include "G4ParticleDefinition.hh"
#include "G4Poisson.hh"
#include <CLHEP/Random/RandomEngine.h>
#include "TrackPropagation/Geant4e/interface/G4TablesForExtrapolatorForCVH.h"

class G4UniversalFluctuationForExtrapolator : public G4VEmFluctuationModel {
public:
  explicit G4UniversalFluctuationForExtrapolator(const G4String& nam = "UniFluc");

  ~G4UniversalFluctuationForExtrapolator() override;

  G4double SampleFluctuations(const G4MaterialCutsCouple*,
                              const G4DynamicParticle*,
                              const G4double,
                              const G4double,
                              const G4double,
                              const G4double) override {
    return 0;
  };

  virtual G4double SampleFluctuations(const G4Material*, const G4DynamicParticle*, G4double, G4double, G4double);

  virtual G4double SampleFluctuations2(
      const G4Material*, const G4DynamicParticle*, G4double, G4double, G4double, G4double);

  G4double Dispersion(
      const G4Material*, const G4DynamicParticle*, const G4double, const G4double, const G4double) override {
    return 0;
  };

  virtual G4double Dispersion(const G4Material*, const G4DynamicParticle*, G4double, G4double);

  // Initialisation for a new particle type
  void InitialiseMe(const G4ParticleDefinition*) override;

  // Initialisation prestep
  void SetParticleAndCharge(const G4ParticleDefinition*, G4double q2) override;

  // hide assignment operator
  G4UniversalFluctuationForExtrapolator& operator=(const G4UniversalFluctuationForExtrapolator& right) = delete;
  G4UniversalFluctuationForExtrapolator(const G4UniversalFluctuationForExtrapolator&) = delete;

protected:
  virtual G4double SampleGlandz(CLHEP::HepRandomEngine* rndm, const G4Material*, const G4double tcut);

  inline void AddExcitation(G4double a, G4double e, G4double& eav, G4double& eloss, G4double& esig2, G4double& esig2tot);

  inline void SampleGauss(G4double eav, G4double esig2, G4double& eloss, G4double& esig2tot);

  inline void AddExcitation2(
      CLHEP::HepRandomEngine* rndm, G4double a, G4double e, G4double& eav, G4double& eloss, G4double& esig2);

  inline void SampleGauss2(CLHEP::HepRandomEngine* rndm, G4double eav, G4double esig2, G4double& eloss);

  // particle properties
  G4double particleMass = 0.0;
  G4double m_Inv_particleMass = DBL_MAX;
  G4double m_massrate = DBL_MAX;
  G4double chargeSquare = 1.0;

  // material properties
  G4double ipotFluct = 0.0;
  G4double ipotLogFluct = 0.0;
  G4double e0 = 0.0;
  G4double electronDensity;
  G4double f1Fluct;
  G4double f2Fluct;
  G4double e1Fluct;
  G4double e2Fluct;
  G4double e1LogFluct;
  G4double e2LogFluct;
  G4double esmall;
  G4double e1, e2;

  // model parameters
  G4double minNumberInteractionsBohr = 10.0;
  G4double minLoss;
  G4double nmaxCont = 8.0;
  G4double rate = 0.56;
  G4double fw = 4.0;
  G4double a0 = 42.0;
  G4double w2 = 0.0;
  G4double meanLoss = 0.0;

  const G4ParticleDefinition* particle = nullptr;
  const G4Material* lastMaterial;
  G4double* rndmarray = nullptr;
  G4int sizearray = 30;

  G4TablesForExtrapolatorForCVH* tables = nullptr;
  const G4PhysicsTable* table = nullptr;
  G4double massratio = 1.;
  G4double charge2ratio = 1.;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4UniversalFluctuationForExtrapolator::AddExcitation(
    G4double ax, G4double ex, G4double& eav, G4double& eloss, G4double& esig2, G4double& esig2tot) {
  if (ax > nmaxCont) {
    eav += ax * ex;
    esig2 += ax * ex * ex;
  } else {
    G4double p = ax;
    if (p > 0) {
      eloss += ((p + 1) - 1.) * ex;
      esig2tot += (ax + 1. / 3.) * ex * ex;
    }
  }
}

inline void G4UniversalFluctuationForExtrapolator::SampleGauss(G4double eav,
                                                               G4double esig2,
                                                               G4double& eloss,
                                                               G4double& esig2tot) {
  G4double x = eav;
  G4double sig = std::sqrt(esig2);
  if (eav < 0.25 * sig) {
    x += (1. - 1.) * eav;
    esig2tot += eav * eav / 3.;
  } else {
    x = eav;
    const double alpha = -eav / sig;
    const double beta = eav / sig;
    const double z = 0.5 * (std::erf(beta / std::sqrt(2.)) - std::erf(alpha / std::sqrt(2.)));
    const double phialpha = 1. / std::sqrt(2. * M_PI) * std::exp(-0.5 * alpha * alpha);
    const double phibeta = 1. / std::sqrt(2. * M_PI) * std::exp(-0.5 * beta * beta);
    esig2tot += sig * sig *
                (1. + (alpha * phialpha - beta * phibeta) / z - (phialpha - phibeta) * (phialpha - phibeta) / z / z);

    // Loop checking, 23-Feb-2016, Vladimir Ivanchenko
  }
  eloss += x;
}

inline void G4UniversalFluctuationForExtrapolator::AddExcitation2(CLHEP::HepRandomEngine* rndm,
                                                                  const G4double ax,
                                                                  const G4double ex,
                                                                  G4double& eav,
                                                                  G4double& eloss,
                                                                  G4double& esig2) {
  if (ax > nmaxCont) {
    eav += ax * ex;
    esig2 += ax * ex * ex;
  } else {
    const G4int p = (G4int)G4Poisson(ax);
    if (p > 0) {
      eloss += ((p + 1) - 2. * rndm->flat()) * ex;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4UniversalFluctuationForExtrapolator::SampleGauss2(CLHEP::HepRandomEngine* rndm,
                                                                const G4double eav,
                                                                const G4double esig2,
                                                                G4double& eloss) {
  G4double x = eav;
  const G4double sig = std::sqrt(esig2);
  if (eav < 0.25 * sig) {
    x += (2. * rndm->flat() - 1.) * eav;
  } else {
    do {
      x = G4RandGauss::shoot(rndm, eav, sig);
    } while (x < 0.0 || x > 2 * eav);
    // Loop checking, 23-Feb-2016, Vladimir Ivanchenko
  }
  eloss += x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
