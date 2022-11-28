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
// GEANT4 Class file
//
//
// File name:     G4UniversalFluctuationForExtrapolator
//
// Author:        V. Ivanchenko for Laszlo Urban
// 
// Creation date: 03.01.2002
//
// Modifications: 
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackPropagation/Geant4e/interface/G4UniversalFluctuationForExtrapolator.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"
#include "G4Step.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4Log.hh"
#include "G4Exp.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4UniversalFluctuationForExtrapolator::G4UniversalFluctuationForExtrapolator(const G4String& nam)
 :G4VEmFluctuationModel(nam),
  particle(nullptr),
  minNumberInteractionsBohr(10.0),
  minLoss(10.*eV),
  nmaxCont(16.),
  rate(0.56),
  a0(50.),
  fw(4.00)
{
  lastMaterial = nullptr;
  particleMass = chargeSquare = ipotFluct = electronDensity = f1Fluct = f2Fluct 
    = e1Fluct = e2Fluct = e1LogFluct = e2LogFluct = ipotLogFluct = e0 = esmall 
    = e1 = e2 = 0.0;
  m_Inv_particleMass = m_massrate = DBL_MAX;
  sizearray = 30;
  rndmarray = new G4double[30];
  tables = new G4TablesForExtrapolatorCustom(0, 70, 1.*MeV, 10.*TeV, true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4UniversalFluctuationForExtrapolator::~G4UniversalFluctuationForExtrapolator()
{
  delete [] rndmarray;
  delete tables;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4UniversalFluctuationForExtrapolator::InitialiseMe(const G4ParticleDefinition* part)
{
  particle       = part;
  particleMass   = part->GetPDGMass();
  G4double q     = part->GetPDGCharge()/eplus;

  // Derived quantities
  m_Inv_particleMass = 1.0 / particleMass;
  m_massrate = electron_mass_c2 * m_Inv_particleMass ;
  chargeSquare   = q*q;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4double
G4UniversalFluctuationForExtrapolator::SampleFluctuations2(const G4Material* material,
                                           const G4DynamicParticle* dp,
                                           G4double tmax,
                                           G4double length,
                                           G4double ekin,
                                           G4double eloss)
{

  // Calculate actual loss from the mean loss.
  // The model used to get the fluctuations is essentially the same
  // as in Glandz in Geant3 (Cern program library W5013, phys332).
  // L. Urban et al. NIM A362, p.416 (1995) and Geant4 Physics Reference Manual

  // shortcut for very small loss or from a step nearly equal to the range
  // (out of validity of the model)
  //
//   G4double meanLoss = averageLoss;
//   const G4PhysicsTable *table = tables->GetPhysicsTable(fDedxMuon);
  size_t idx = 0;
  G4double dedx = ((*table)[material->GetIndex()])->Value(massratio*ekin, idx)*charge2ratio;
  G4double meanLoss = length*dedx;
//   std::cout << "meanLoss = " << meanLoss << std::endl;
//   G4double tkin  = dp->GetKineticEnergy();
  G4double tkin  = ekin;

  // std::cout << "eloss = " << eloss << " meanLoss = " << meanLoss << std::endl;

  // const double extraloss = eloss - meanLoss;
  const double extraloss = 0.;

  if (meanLoss < minLoss) { return meanLoss + extraloss; }

  if(dp->GetDefinition() != particle) { InitialiseMe(dp->GetDefinition()); }

  CLHEP::HepRandomEngine* rndmEngineF = G4Random::getTheEngine();

  G4double tau   = tkin * m_Inv_particleMass;
  G4double gam   = tau + 1.0;
  G4double gam2  = gam*gam;
  G4double beta2 = tau*(tau + 2.0)/gam2;

  G4double loss(0.), siga(0.);

  // const G4Material* material = couple->GetMaterial();

  // Gaussian regime
  // for heavy particles only and conditions
  // for Gauusian fluct. has been changed
  //
  if ((particleMass > electron_mass_c2) &&
      (meanLoss >= minNumberInteractionsBohr*tmax))
  {
    G4double tmaxkine = 2.*electron_mass_c2*beta2*gam2/
                        (1.+m_massrate*(2.*gam+m_massrate)) ;
    if (tmaxkine <= 2.*tmax)
    {
      electronDensity = material->GetElectronDensity();
      siga = sqrt((1.0/beta2 - 0.5) * twopi_mc2_rcl2 * tmax * length
                  * electronDensity * chargeSquare);

      G4double sn = meanLoss/siga;

      // thick target case
      if (sn >= 2.0) {

        G4double twomeanLoss = meanLoss + meanLoss;
        do {
          loss = G4RandGauss::shoot(rndmEngineF,meanLoss,siga);
          // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
        } while  (0.0 > loss || twomeanLoss < loss);

        // Gamma distribution
      } else {

        G4double neff = sn*sn;
        loss = meanLoss*G4RandGamma::shoot(rndmEngineF,neff,1.0)/neff;
      }
      //G4cout << "Gauss: " << loss << G4endl;
      return loss + extraloss;
    }
  }

  // Glandz regime : initialisation
  //
  if (material != lastMaterial) {
    f1Fluct      = material->GetIonisation()->GetF1fluct();
    f2Fluct      = material->GetIonisation()->GetF2fluct();
    e1Fluct      = material->GetIonisation()->GetEnergy1fluct();
    e2Fluct      = material->GetIonisation()->GetEnergy2fluct();
    e1LogFluct   = material->GetIonisation()->GetLogEnergy1fluct();
    e2LogFluct   = material->GetIonisation()->GetLogEnergy2fluct();
    ipotFluct    = material->GetIonisation()->GetMeanExcitationEnergy();
    ipotLogFluct = material->GetIonisation()->GetLogMeanExcEnergy();
    e0 = material->GetIonisation()->GetEnergy0fluct();
    esmall = 0.5*sqrt(e0*ipotFluct);
    lastMaterial = material;
  }

  // very small step or low-density material
  if(tmax <= e0) { return meanLoss + extraloss; }

  // width correction for small cuts
  G4double scaling = std::min(1.+0.5*CLHEP::keV/tmax,1.50);
  meanLoss /= scaling;

  G4double a1(0.0), a2(0.0), a3(0.0);

  loss = 0.0;

  e1 = e1Fluct;
  e2 = e2Fluct;

  if(tmax > ipotFluct) {
    G4double w2 = G4Log(2.*electron_mass_c2*beta2*gam2)-beta2;

    if(w2 > ipotLogFluct)  {
      if(w2 > e2LogFluct) {
	G4double C = meanLoss*(1.-rate)/(w2-ipotLogFluct);
	a1 = C*f1Fluct*(w2-e1LogFluct)/e1Fluct;
	a2 = C*f2Fluct*(w2-e2LogFluct)/e2Fluct;
      } else {
	a1 = meanLoss*(1.-rate)/e1;
      }
      if(a1 < a0) {
        G4double fwnow = 0.5+(fw-0.5)*sqrt(a1/a0);
        a1 /= fwnow;
        e1 *= fwnow;
      } else {
        a1 /= fw;
        e1 = fw*e1Fluct;
      }
    }
  }

  G4double w1 = tmax/e0;
  if(tmax > e0) {
    a3 = rate*meanLoss*(tmax-e0)/(e0*tmax*G4Log(w1));
    if(a1+a2 <= 0.) {
      a3 /= rate;
    }
  }
  //'nearly' Gaussian fluctuation if a1>nmaxCont&&a2>nmaxCont&&a3>nmaxCont
  G4double emean = 0.;
  G4double sig2e = 0.;

  // excitation of type 1
  if(a1 > 0.0) { AddExcitation2(rndmEngineF, a1, e1, emean, loss, sig2e); }

  // excitation of type 2
  if(a2 > 0.0) { AddExcitation2(rndmEngineF, a2, e2, emean, loss, sig2e); }

  if(sig2e > 0.0) { SampleGauss2(rndmEngineF, emean, sig2e, loss); }

  // ionisation
  if(a3 > 0.) {
    emean = 0.;
    sig2e = 0.;
    G4double p3 = a3;
    G4double alfa = 1.;
    if(a3 > nmaxCont)
      {
        alfa            = w1*(nmaxCont+a3)/(w1*nmaxCont+a3);
        G4double alfa1  = alfa*G4Log(alfa)/(alfa-1.);
        G4double namean = a3*w1*(alfa-1.)/((w1-1.)*alfa);
        emean          += namean*e0*alfa1;
        sig2e          += e0*e0*namean*(alfa-alfa1*alfa1);
        p3              = a3-namean;
      }

    G4double w2 = alfa*e0;
    if(tmax > w2) {
      G4double w  = (tmax-w2)/tmax;
      G4int nnb = G4Poisson(p3);
      if(nnb > 0) {
        if(nnb > sizearray) {
          sizearray = nnb;
          delete [] rndmarray;
          rndmarray = new G4double[nnb];
        }
        rndmEngineF->flatArray(nnb, rndmarray);
        for (G4int k=0; k<nnb; ++k) { loss += w2/(1.-w*rndmarray[k]); }
      }
    }
    if(sig2e > 0.0) { SampleGauss2(rndmEngineF, emean, sig2e, loss); }
  }

  loss *= scaling;

  return loss + extraloss;

}

G4double 
G4UniversalFluctuationForExtrapolator::SampleFluctuations(const G4Material* material,
                                           const G4DynamicParticle* dp,
                                           G4double tmax,
                                           G4double length,
                                           G4double ekin)
{
  // Calculate actual loss from the mean loss.
  // The model used to get the fluctuations is essentially the same
  // as in Glandz in Geant3 (Cern program library W5013, phys332).
  // L. Urban et al. NIM A362, p.416 (1995) and Geant4 Physics Reference Manual

  // shortcut for very small loss or from a step nearly equal to the range
  // (out of validity of the model)
  //
//   G4double meanLoss = averageLoss;
//   const G4PhysicsTable *table = tables->GetPhysicsTable(fDedxMuon);
  size_t idx = 0;
  G4double dedx = ((*table)[material->GetIndex()])->Value(massratio*ekin, idx)*charge2ratio;
  G4double meanLoss = length*dedx;
//   std::cout << "meanLoss = " << meanLoss << std::endl;
//   G4double tkin  = dp->GetKineticEnergy();
  G4double tkin  = ekin;
  //G4cout<< "Emean= "<< meanLoss<< " tmax= "<< tmax<< " L= "<<length<<G4endl;
  if (meanLoss < minLoss) { return 0.; }

  if(dp->GetDefinition() != particle) { InitialiseMe(dp->GetDefinition()); }

//   CLHEP::HepRandomEngine* rndmEngineF = G4Random::getTheEngine();
  
  G4double tau   = tkin * m_Inv_particleMass;            
  G4double gam   = tau + 1.0;
  G4double gam2  = gam*gam;
  G4double beta2 = tau*(tau + 2.0)/gam2;

  G4double loss(0.), siga(0.);

//   const G4Material* material = couple->GetMaterial();
  
  // Gaussian regime
  // for heavy particles only and conditions
  // for Gauusian fluct. has been changed 
  //
  if ((particleMass > electron_mass_c2) &&
      (meanLoss >= minNumberInteractionsBohr*tmax))
  {
    G4double tmaxkine = 2.*electron_mass_c2*beta2*gam2/
                        (1.+m_massrate*(2.*gam+m_massrate)) ;
    if (tmaxkine <= 2.*tmax)   
    {
      electronDensity = material->GetElectronDensity();
      siga = sqrt((1.0/beta2 - 0.5) * twopi_mc2_rcl2 * tmax * length
                  * electronDensity * chargeSquare);

      G4double sn = meanLoss/siga;
  
      // thick target case 
      if (sn >= 2.0) {

//         std::cout << "gaussian case\n";
        return siga*siga;

        G4double twomeanLoss = meanLoss + meanLoss;
        do {
//           loss = G4RandGauss::shoot(rndmEngineF,meanLoss,siga);
          loss  = meanLoss;
          // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
        } while  (0.0 > loss || twomeanLoss < loss);

        // Gamma distribution
      } else {

//         std::cout << "gamma case\n";
        G4double neff = sn*sn;
        return meanLoss*meanLoss/neff;
        loss = meanLoss;
//         loss = meanLoss*G4RandGamm:a:shoot(rndmEngineF,neff,1.0)/neff;
      }
      //G4cout << "Gauss: " << loss << G4endl;
      return loss;
    }
  }

  // Glandz regime : initialisation
  //
  if (material != lastMaterial) {
    f1Fluct      = material->GetIonisation()->GetF1fluct();
    f2Fluct      = material->GetIonisation()->GetF2fluct();
    e1Fluct      = material->GetIonisation()->GetEnergy1fluct();
    e2Fluct      = material->GetIonisation()->GetEnergy2fluct();
    e1LogFluct   = material->GetIonisation()->GetLogEnergy1fluct();
    e2LogFluct   = material->GetIonisation()->GetLogEnergy2fluct();
    ipotFluct    = material->GetIonisation()->GetMeanExcitationEnergy();
    ipotLogFluct = material->GetIonisation()->GetLogMeanExcEnergy();
    e0 = material->GetIonisation()->GetEnergy0fluct();
    esmall = 0.5*sqrt(e0*ipotFluct);  
    lastMaterial = material;   
  }

  // very small step or low-density material
  if(tmax <= e0) { return 0.; }

  // width correction for small cuts
  G4double scaling = std::min(1.+0.5*CLHEP::keV/tmax,1.50);
  meanLoss /= scaling;
  
//   std::cout << "urban model: scaling = " << scaling << std::endl;

  G4double a1(0.0), a2(0.0), a3(0.0);
    
  loss = 0.0;

  e1 = e1Fluct;
  e2 = e2Fluct;

  if(tmax > ipotFluct) {
    G4double w2 = G4Log(2.*electron_mass_c2*beta2*gam2)-beta2;

    if(w2 > ipotLogFluct)  {
      if(w2 > e2LogFluct) {
	G4double C = meanLoss*(1.-rate)/(w2-ipotLogFluct);
	a1 = C*f1Fluct*(w2-e1LogFluct)/e1Fluct;
	a2 = C*f2Fluct*(w2-e2LogFluct)/e2Fluct;
      } else {
	a1 = meanLoss*(1.-rate)/e1;
      }
      if(a1 < a0) { 
        G4double fwnow = 0.5+(fw-0.5)*sqrt(a1/a0);
        a1 /= fwnow;
        e1 *= fwnow;
      } else {
        a1 /= fw;
        e1 = fw*e1Fluct;
      }
    }   
  }

  G4double w1 = tmax/e0;
  if(tmax > e0) {
    a3 = rate*meanLoss*(tmax-e0)/(e0*tmax*G4Log(w1));
    if(a1+a2 <= 0.) { 
      a3 /= rate;
    }
  }
  //'nearly' Gaussian fluctuation if a1>nmaxCont&&a2>nmaxCont&&a3>nmaxCont  
  G4double emean = 0.;
  G4double sig2e = 0.;
  G4double esig2tot = 0.;

  // excitation of type 1
  if(a1 > 0.0) { AddExcitation(a1, e1, emean, loss, sig2e, esig2tot); }

  // excitation of type 2
  if(a2 > 0.0) { AddExcitation(a2, e2, emean, loss, sig2e, esig2tot); }

  if(sig2e > 0.0) { SampleGauss(emean, sig2e, loss, esig2tot); }

  // ionisation 
  if(a3 > 0.) {
    emean = 0.;
    sig2e = 0.;
    G4double p3 = a3;
    G4double alfa = 1.;
    if(a3 > nmaxCont)
      {
        alfa            = w1*(nmaxCont+a3)/(w1*nmaxCont+a3);
        G4double alfa1  = alfa*G4Log(alfa)/(alfa-1.);
        G4double namean = a3*w1*(alfa-1.)/((w1-1.)*alfa);
        emean          += namean*e0*alfa1;
        sig2e          += e0*e0*namean*(alfa-alfa1*alfa1);
        p3              = a3-namean;
      }

    G4double w2 = alfa*e0;
    if(tmax > w2) {
      G4double w  = (tmax-w2)/tmax;
//       G4int nnb = G4Poisson(p3);
//       if(nnb > 0) {
//         if(nnb > sizearray) {
//           sizearray = nnb;
//           delete [] rndmarray;
//           rndmarray = new G4double[nnb];
//         }
//         rndmEngineF->flatArray(nnb, rndmarray);
//         for (G4int k=0; k<nnb; ++k) { loss += w2/(1.-w*rndmarray[k]); }
//       }
//       const double f = -std::log(1.-w)*w2/w;
//       const double f2 = w2*w2/(1.-w);
      // const double alpha = 0.996;
//       const double alpha = 0.9999;
      // const double alpha = 1.;
      // const double alpha = 1. - 1e-1;
      const double alpha = 0.999;
//       const double ualpha = (1. - std::pow(1. - w, alpha))/w;
//       std::cout << "w = " << w << " ualpha = " << ualpha << std::endl;
      const double ualpha = alpha;
      const double f = -std::log(1.-ualpha*w)*w2/w;
      const double f2 = ualpha*w2*w2/(1.-ualpha*w);
      const double sigf2 = f2 - f*f;
      
      //w2/(1-uf*w) = f
      //w2 = f - f*w*uf
      //uf = (f - w2)/(fw)
//       const double uf = (f - w2)/(f*w);
      
      //sigf2 = w*w*w2*w2/(1. - uf*w)^4 * siguf^2
      // siguf^2 = sigf2*(1-uf*w)^4/(w*w*w2*w2)
//       const double siguf2 = sigf2*std::pow(1. - uf*w, 4)/(w*w*w2*w2);
      
//       std::cout << "uf = " << uf << std::endl;
      loss += p3*f;
      esig2tot += f*f*p3 + p3*sigf2;
      
//       esig2tot += f*f*p3 + p3*p3*sigf2;
//       loss += p3*w2/(1. - 0.5*w);
//       esig2tot += p3*w2*w2/std::pow(1. - 0.5*w, 2) + p3*p3*w*w*w2*w2/std::pow(1. - 0.5*w, 4)/12.;
//       esig2tot += p3*w2*w2/std::pow(1. - 0.5*w, 2) + p3*w*w*w2*w2/std::pow(1. - 0.5*w, 4)/12.;
//       esig2tot += p3*w2*w2/std::pow(1. - uf*w, 2) + p3*w*w*w2*w2/std::pow(1. - uf*w, 4)/12.;
    }
    if(sig2e > 0.0) { SampleGauss(emean, sig2e, loss, esig2tot); }
  }

  loss *= scaling;
  esig2tot *= scaling*scaling;
  
//   const double lossratio = loss/meanLoss;
//   std::cout << "lossratio = " << lossratio << std::endl;

//   return 0.1*esig2tot;
  return esig2tot;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4double G4UniversalFluctuationForExtrapolator::Dispersion(
                          const G4Material* material,
                          const G4DynamicParticle* dp,
                                G4double tmax,
                                G4double length)
{
  if(dp->GetDefinition() != particle) { InitialiseMe(dp->GetDefinition()); }

  electronDensity = material->GetElectronDensity();

  G4double gam   = (dp->GetKineticEnergy())*m_Inv_particleMass + 1.0;
  G4double beta2 = 1.0 - 1.0/(gam*gam);

  G4double siga  = (1.0/beta2 - 0.5) * twopi_mc2_rcl2 * tmax * length
                 * electronDensity * chargeSquare;

  return siga;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4UniversalFluctuationForExtrapolator::SetParticleAndCharge(const G4ParticleDefinition* part,
                                             G4double q2)
{
  if(part != particle) {
    particle       = part;
    particleMass   = part->GetPDGMass();

    // Derived quantities
    if( particleMass != 0.0 ){
      m_Inv_particleMass = 1.0 / particleMass;
      m_massrate = electron_mass_c2 * m_Inv_particleMass ;
    }else{
      m_Inv_particleMass = DBL_MAX;
      m_massrate = DBL_MAX;
    }
  }
  chargeSquare = q2;
  
  if (part == G4Electron::Electron()) {
    table = tables->GetPhysicsTable(fDedxElectron);
    massratio = 1.;
    charge2ratio = 1.;
  }
  else if (part == G4Positron::Positron()) {
    table = tables->GetPhysicsTable(fDedxPositron);
    massratio = 1.;
    charge2ratio = 1.;
  }
  else if (part == G4MuonPlus::MuonPlus() || part == G4MuonMinus::MuonMinus()) {
//     std::cout << "G4UniversalFluctuationForExtrapolator::SetParticleAndCharge: muon tables\n";
    table = tables->GetPhysicsTable(fDedxMuon);
    massratio = 1.;
    charge2ratio = 1.;
  }
  else {
    // scaled energy loss from proton tables
    table = tables->GetPhysicsTable(fDedxProton);
    massratio = proton_mass_c2/particleMass;
    charge2ratio = part->GetPDGCharge()*part->GetPDGCharge();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
