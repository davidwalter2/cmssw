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
// $Id$
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
// Implementation of energy loss fluctuations

// -------------------------------------------------------------------
//

#ifndef G4UniversalFluctuationForExtrapolator_h
#define G4UniversalFluctuationForExtrapolator_h 1


#include "G4VEmFluctuationModel.hh"
#include "G4ParticleDefinition.hh"
#include "G4Poisson.hh"
#include <CLHEP/Random/RandomEngine.h>
#include "TrackPropagation/Geant4e/interface/G4TablesForExtrapolatorCustom.h"

class G4UniversalFluctuationForExtrapolator : public G4VEmFluctuationModel
{

public:

  explicit G4UniversalFluctuationForExtrapolator(const G4String& nam = "UniFluc");

  virtual ~G4UniversalFluctuationForExtrapolator();

  virtual G4double SampleFluctuations(const G4Material*,
                                      const G4DynamicParticle*,
                                      G4double,
                                      G4double,
                                      G4double);

  virtual G4double SampleFluctuations2(const G4Material*,
                                      const G4DynamicParticle*,
                                      G4double,
                                      G4double,
                                      G4double,
                                      G4double);

  virtual G4double SampleFluctuations(const G4MaterialCutsCouple*,
                                      const G4DynamicParticle*,
                                      G4double,
                                      G4double,
                                      G4double) override { return 0.; }

  virtual G4double Dispersion(const G4Material*,
                              const G4DynamicParticle*,
                              G4double,
                              G4double) override;

  virtual void InitialiseMe(const G4ParticleDefinition*) final;

  // Initialisation prestep
  virtual void SetParticleAndCharge(const G4ParticleDefinition*, 
                                    G4double q2) final;

private:

  inline void AddExcitation(G4double a, G4double e, G4double& eav,
                            G4double& eloss, G4double& esig2, G4double& esig2tot);

  inline void SampleGauss(G4double eav, G4double esig2,
                          G4double& eloss, G4double& esig2tot);

  inline void AddExcitation2(CLHEP::HepRandomEngine* rndm,
                            G4double a, G4double e, G4double& eav,
                            G4double& eloss, G4double& esig2);

  inline void SampleGauss2(CLHEP::HepRandomEngine* rndm,
                          G4double eav, G4double esig2,
                          G4double& eloss);

  // hide assignment operator
  G4UniversalFluctuationForExtrapolator & operator=(const  G4UniversalFluctuationForExtrapolator &right) = delete;
  G4UniversalFluctuationForExtrapolator(const  G4UniversalFluctuationForExtrapolator&) = delete;

  const G4ParticleDefinition* particle;
  const G4Material* lastMaterial;

  G4double particleMass;

  // Derived quantities
  G4double m_Inv_particleMass;
  G4double m_massrate;
  G4double chargeSquare;

  // data members to speed up the fluctuation calculation
  G4double ipotFluct;
  G4double electronDensity;
  
  G4double f1Fluct;
  G4double f2Fluct;
  G4double e1Fluct;
  G4double e2Fluct;
  G4double e1LogFluct;
  G4double e2LogFluct;
  G4double ipotLogFluct;
  G4double e0;
  G4double esmall;

  G4double e1,e2;

  G4double minNumberInteractionsBohr;
  G4double minLoss;
  G4double nmaxCont;
  G4double rate,a0,fw; 

  G4int     sizearray;
  G4double* rndmarray;

  G4TablesForExtrapolatorCustom *tables = nullptr;
  const G4PhysicsTable *table = nullptr;
  G4double massratio = 1.;
  G4double charge2ratio = 1.;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


inline void 
G4UniversalFluctuationForExtrapolator::AddExcitation(G4double ax, G4double ex, G4double& eav,
                                      G4double& eloss, G4double& esig2, G4double& esig2tot)
{
  if(ax > nmaxCont) {
    eav  += ax*ex;
    esig2 += ax*ex*ex;
  } else {
//     G4int p = G4Poisson(ax);
    G4double p = ax;
    if(p > 0) {
//       eloss += ((p + 1) - 2.*rndm->flat())*ex;
      eloss += ((p + 1) - 1.)*ex;
      esig2tot += (ax + 1./3.)*ex*ex;
    }
  }
}

inline void 
G4UniversalFluctuationForExtrapolator::SampleGauss(G4double eav, G4double esig2,
                                    G4double& eloss, G4double &esig2tot)
{
  G4double x = eav;
  G4double sig = std::sqrt(esig2);
  if(eav < 0.25*sig) {
//     x += (2.*rndm->flat() - 1.)*eav;
    x += (1. - 1.)*eav;
    esig2tot += eav*eav/3.;
  } else {
//     do {
// //       x = G4RandGauss::shoot(rndm, eav, sig);
//     } while (x < 0.0 || x > 2*eav);
    x = eav;
    // variance of truncated gaussian
//     esig2tot += sig*sig;
    const double alpha = -eav/sig;
    const double beta = eav/sig;
    const double z = 0.5*(std::erf(beta/std::sqrt(2.)) - std::erf(alpha/std::sqrt(2.)));
    const double phialpha = 1./std::sqrt(2.*M_PI)*std::exp(-0.5*alpha*alpha);
    const double phibeta = 1./std::sqrt(2.*M_PI)*std::exp(-0.5*beta*beta);
    esig2tot += sig*sig*(1. + (alpha*phialpha - beta*phibeta)/z - (phialpha - phibeta)*(phialpha - phibeta)/z/z);

    // Loop checking, 23-Feb-2016, Vladimir Ivanchenko
  }
  eloss += x;
} 

inline void
G4UniversalFluctuationForExtrapolator::AddExcitation2(CLHEP::HepRandomEngine* rndm,
                                      G4double ax, G4double ex, G4double& eav,
                                      G4double& eloss, G4double& esig2)
{
  if(ax > nmaxCont) {
    eav  += ax*ex;
    esig2 += ax*ex*ex;
  } else {
    G4int p = G4Poisson(ax);
    if(p > 0) { eloss += ((p + 1) - 2.*rndm->flat())*ex; }
  }
}

inline void
G4UniversalFluctuationForExtrapolator::SampleGauss2(CLHEP::HepRandomEngine* rndm,
                                    G4double eav, G4double esig2,
                                    G4double& eloss)
{
  G4double x = eav;
  G4double sig = std::sqrt(esig2);
  if(eav < 0.25*sig) {
    x += (2.*rndm->flat() - 1.)*eav;
  } else {
    do {
      x = G4RandGauss::shoot(rndm, eav, sig);
    } while (x < 0.0 || x > 2*eav);
    // Loop checking, 23-Feb-2016, Vladimir Ivanchenko
  }
  eloss += x;
}

#endif

