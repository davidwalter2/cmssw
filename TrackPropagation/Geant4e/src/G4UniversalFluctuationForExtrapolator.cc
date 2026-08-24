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
#include "TrackPropagation/Geant4e/interface/CGFQoPBlock.h"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4Log.hh"
#include "FWCore/Utilities/interface/Exception.h"

#include <cstdlib>
#include <cmath>

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Process-wide shared dE/dx / range / inv-range tables. Built once on the
// first ctor call (under extrMutex in MT builds); subsequent instances --
// including the per-stream Geant4ePropagator clones in the CVH refit --
// reuse the same pointer. Tables are read-only after construction. Leaked
// at process exit to keep ownership trivial.
G4TablesForExtrapolatorForCVH* G4UniversalFluctuationForExtrapolator::tables = nullptr;
#ifdef G4MULTITHREADED
G4Mutex G4UniversalFluctuationForExtrapolator::extrMutex = G4MUTEX_INITIALIZER;
#endif

G4UniversalFluctuationForExtrapolator::G4UniversalFluctuationForExtrapolator(const G4String& nam)
    : G4VEmFluctuationModel(nam), minLoss(10. * CLHEP::eV) {
  rndmarray = new G4double[sizearray];

  if (nullptr == tables) {
#ifdef G4MULTITHREADED
    G4MUTEXLOCK(&extrMutex);
    if (nullptr == tables) {
#endif
      // Grid from G4TablesForExtrapolatorForCVH, not restated here (it used to
      // be 70 bins to 10 TeV against the reference's 80 to 100 TeV).
      // `iononly = true` stays: this instance supplies the IONIZATION mean
      // loss, and the radiative fluctuation is a separate channel.
      tables = new G4TablesForExtrapolatorForCVH(0,
                                                G4TablesForExtrapolatorForCVH::kNbins,
                                                G4TablesForExtrapolatorForCVH::kEminMeV * CLHEP::MeV,
                                                G4TablesForExtrapolatorForCVH::kEmaxMeV * CLHEP::MeV,
                                                true);
#ifdef G4MULTITHREADED
    }
    G4MUTEXUNLOCK(&extrMutex);
#endif
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool G4UniversalFluctuationForExtrapolator::urban2021Enabled() {
  // DEFAULT OFF. NOTES_DELTASPEC s10.4 falsified the specific prediction this
  // harmonization was made for; it stays available as a diagnostic and is not
  // one of the four corrections that went default-on 2026-08-16.
  return cvhcgf::switches().ioniUrban2021;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool G4UniversalFluctuationForExtrapolator::exactDeltaEnabled() {
  // DEFAULT ON, so the model represents what the simulation does: Geant4
  // samples the exact knock-on cross section, and the 1/E^2 delta channel is an
  // approximation to it.  See the block comment on the defaults in
  // CGFQoPBlock.cc, including the attribution gate that is still open.
  // CVH_IONI_EXACTDELTA=1 selects the
  // exact knock-on cross section (regime 2/3, 13-wide `ioniurbanv` stride);
  // unset or =0 is the 1/E^2 delta channel with the 11-wide stride.
  return cvhcgf::switches().ioniExactDelta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// THE KOKOULIN VARIANCE INTEGRAL WAS HERE and went with the truncation
// (2026-08-24). It integrated (f_K - 1) over the ALPHA-TRUNCATED range to
// correct the variance this class returned. The variance is the record's own
// cumulant now, and `cvhcgf::kokoulinKappa2` does the same integral over the
// FULL range as part of it -- one formula, and it is the one the block CGF
// uses too. `cvhcgf::kokoulinFactor` carries f_K itself.

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

cvhcgf::IoniStep G4UniversalFluctuationForExtrapolator::toIoniStep(const UrbanFluctRecord& rec, int kokNbin) {
  // ONE mapping from the record to a block step, shared by this class (for the
  // variance it returns) and by Geant4ePropagator (for the block CGF). Two
  // copies is how the two would come to disagree about whether `a3` is a count
  // or an energy, or whether `scaling` has been applied -- both silent, both
  // ~1e-5. `gs` stays 1: the caller supplies the transport weight and the
  // standardization.
  const double gam = rec.scaling;
  cvhcgf::IoniStep s;
  s.regime = rec.regime;
  s.gsig2 = rec.gsig2;
  s.a1 = rec.a1;
  s.e1 = rec.e1 * gam;
  s.a2 = rec.a2;
  s.e2 = rec.e2 * gam;
  s.a3 = (rec.regime >= 2) ? rec.a3 * gam : rec.a3;  // count in regime 1, xi (an energy) in 2/3
  s.e0 = rec.e0r * gam;
  s.tmax = rec.tmaxr * gam;
  s.beta2 = rec.beta2;
  s.etot = rec.etot;
  s.kokNbin = (rec.regime >= 2) ? kokNbin : 0;
  s.gs = 1.;
  return s;
}

namespace {
  double blockKappa2Of(const G4UniversalFluctuationForExtrapolator::UrbanFluctRecord& rec, int kokNbin) {
    cvhcgf::Block b;
    b.ioni.push_back(G4UniversalFluctuationForExtrapolator::toIoniStep(rec, kokNbin));
    return cvhcgf::blockKappa2(b);
  }
}  // namespace

G4UniversalFluctuationForExtrapolator::~G4UniversalFluctuationForExtrapolator() {
  delete[] rndmarray;
  // tables is a process-wide shared static; do not delete here. The OS
  // reclaims memory on process exit.
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4UniversalFluctuationForExtrapolator::InitialiseMe(const G4ParticleDefinition* part) {
  particle = part;
  particleMass = part->GetPDGMass();
  const G4double q = part->GetPDGCharge() / CLHEP::eplus;

  // Derived quantities
  m_Inv_particleMass = 1.0 / particleMass;
  m_massrate = CLHEP::electron_mass_c2 * m_Inv_particleMass;
  chargeSquare = q * q;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4UniversalFluctuationForExtrapolator::SampleFluctuations2(const G4Material* material,
                                                                    const G4DynamicParticle* dp,
                                                                    G4double tmax,
                                                                    G4double length,
                                                                    G4double ekin,
                                                                    G4double eloss) {
  // Calculate actual loss from the mean loss.
  // The model used to get the fluctuations is essentially the same
  // as in Glandz in Geant3 (Cern program library W5013, phys332).
  // L. Urban et al. NIM A362, p.416 (1995) and Geant4 Physics Reference Manual

  // shortcut for very small loss or from a step nearly equal to the range
  // (out of validity of the model)
  //

  // Full rebind (not just InitialiseMe) BEFORE the dedx lookup -- see
  // SampleFluctuations for the rationale.
  if (dp->GetDefinition() != particle) {
    const G4double q = dp->GetDefinition()->GetPDGCharge() / CLHEP::eplus;
    SetParticleAndCharge(dp->GetDefinition(), q * q);
  }

  size_t idx = 0;
  G4double dedx = ((*table)[material->GetIndex()])->Value(massratio * ekin, idx) * charge2ratio;
  // SPECIES-DEPENDENT Tmax (CVH_REF_SPECIESDEDX). The SAME cvhcgf helper the
  // reference's ComputeDEDX uses, so the two cannot drift apart; exactly +0.0
  // for a proton/antiproton and never reached by a muon.
  if (speciesDedx && usesScaledProtonTable) {
    dedx += cvhcgf::speciesTmaxDedx(ekin, particleMass, speciesRefMass, chargeSquare, material->GetElectronDensity());
  }
  G4double meanLoss = length * dedx;

  G4double tkin = ekin;

  const double extraloss = 0.;

  if (meanLoss < minLoss) {
    return meanLoss + extraloss;
  }

  CLHEP::HepRandomEngine* rndmEngineF = G4Random::getTheEngine();

  G4double tau = tkin * m_Inv_particleMass;
  G4double gam = tau + 1.0;
  G4double gam2 = gam * gam;
  G4double beta2 = tau * (tau + 2.0) / gam2;

  G4double loss(0.), siga(0.);

  // Gaussian regime
  // for heavy particles only and conditions
  // for Gauusian fluct. has been changed
  //
  if ((particleMass > electron_mass_c2) && (meanLoss >= minNumberInteractionsBohr * tmax)) {
    G4double tmaxkine = 2. * electron_mass_c2 * beta2 * gam2 / (1. + m_massrate * (2. * gam + m_massrate));
    if (tmaxkine <= 2. * tmax) {
      electronDensity = material->GetElectronDensity();
      siga = sqrt((1.0 / beta2 - 0.5) * twopi_mc2_rcl2 * tmax * length * electronDensity * chargeSquare);

      G4double sn = meanLoss / siga;

      // thick target case
      if (sn >= 2.0) {
        G4double twomeanLoss = meanLoss + meanLoss;
        do {
          loss = G4RandGauss::shoot(rndmEngineF, meanLoss, siga);
          // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
        } while (0.0 > loss || twomeanLoss < loss);

        // Gamma distribution
      } else {
        G4double neff = sn * sn;
        loss = meanLoss * G4RandGamma::shoot(rndmEngineF, neff, 1.0) / neff;
      }
      return loss + extraloss;
    }
  }

  // Glandz regime : initialisation
  //
  if (material != lastMaterial) {
    f1Fluct = material->GetIonisation()->GetF1fluct();
    f2Fluct = material->GetIonisation()->GetF2fluct();
    e1Fluct = material->GetIonisation()->GetEnergy1fluct();
    e2Fluct = material->GetIonisation()->GetEnergy2fluct();
    e1LogFluct = material->GetIonisation()->GetLogEnergy1fluct();
    e2LogFluct = material->GetIonisation()->GetLogEnergy2fluct();
    ipotFluct = material->GetIonisation()->GetMeanExcitationEnergy();
    ipotLogFluct = material->GetIonisation()->GetLogMeanExcEnergy();
    e0 = material->GetIonisation()->GetEnergy0fluct();
    esmall = 0.5 * sqrt(e0 * ipotFluct);
    lastMaterial = material;
  }

  // very small step or low-density material
  if (tmax <= e0) {
    return meanLoss + extraloss;
  }

  // width correction for small cuts
  G4double scaling = std::min(1. + 0.5 * CLHEP::keV / tmax, 1.50);
  meanLoss /= scaling;

  G4double a1(0.0), a2(0.0), a3(0.0);

  loss = 0.0;

  e1 = e1Fluct;
  e2 = e2Fluct;

  // ------------------------------------------------------------------------
  // CVH_IONI_URBAN2021 -- harmonize the EXCITATION channels with the model the
  // simulation actually samples.
  //
  // This class is the PRE-2021 Urban model: two excitation channels at
  // e1Fluct and e2Fluct = 10 Zeff^2 eV with weights f1 = 1 - 2/Zeff,
  // f2 = 2/Zeff, and the fwnow floor 0.5.  Stock Geant4 11.2.2 -- which is
  // what `g4SimHits` runs, i.e. what the clean-propagation "data" is drawn
  // from -- is the 2021 model: ONE excitation channel at the mean excitation
  // energy ipotFluct, and the fwnow floor 0.1 (G4UniversalFluctuation.hh in
  // 11.2.2 has no f1Fluct/f2Fluct/e1Fluct/e2Fluct members at all).
  // NOTES_URBANSAMPLING measured the gap: +28 % excitation kappa2, +1.1 % of a
  // plane step's core variance, validated against 2.5e7 real stock samples.
  //
  // Two implementations of the same physics differing by 28 % in a channel is
  // a defect whether or not it moves the closure, so this brings the
  // extrapolator FORWARD to the simulation (stock Geant4 cannot be changed).
  // The body below is the 2021 Glandz initialisation verbatim -- it is the
  // same code as this class's own (dead) SampleGlandz, which is a copy of the
  // 2021 model.
  //
  // a2 is set to ZERO and e2 is LEFT AT e2Fluct: the offline consumers all
  // guard on `aj > 0`, so a zero-weight channel is inert, while keeping e2
  // preserves the material identification (e2 = 10 Zeff^2 eV) that the
  // downstream analysis uses.  The record stride is unchanged.
  //
  // DEFAULT OFF; composes with CVH_IONI_EXACTDELTA, which rescales a1 and a2
  // by a common factor and therefore does the right thing when a2 = 0.
  if (urban2021Enabled()) {
    e2 = e2Fluct;
    a2 = 0.0;
    e1 = ipotFluct;
    if (tmax > e1) {
      a1 = meanLoss * (1. - rate) / e1;
      if (a1 < a0) {
        const G4double fwnow = 0.1 + (fw - 0.1) * std::sqrt(a1 / a0);
        a1 /= fwnow;
        e1 *= fwnow;
      } else {
        a1 /= fw;
        e1 *= fw;
      }
    }
  } else if (tmax > ipotFluct) {
    G4double w2 = G4Log(2. * electron_mass_c2 * beta2 * gam2) - beta2;

    if (w2 > ipotLogFluct) {
      if (w2 > e2LogFluct) {
        G4double C = meanLoss * (1. - rate) / (w2 - ipotLogFluct);
        a1 = C * f1Fluct * (w2 - e1LogFluct) / e1Fluct;
        a2 = C * f2Fluct * (w2 - e2LogFluct) / e2Fluct;
      } else {
        a1 = meanLoss * (1. - rate) / e1;
      }
      if (a1 < a0) {
        G4double fwnow = 0.5 + (fw - 0.5) * sqrt(a1 / a0);
        a1 /= fwnow;
        e1 *= fwnow;
      } else {
        a1 /= fw;
        e1 = fw * e1Fluct;
      }
    }
  }
  // ------------------------------------------------------------------------

  G4double w1 = tmax / e0;
  if (tmax > e0) {
    a3 = rate * meanLoss * (tmax - e0) / (e0 * tmax * G4Log(w1));
    if (a1 + a2 <= 0.) {
      a3 /= rate;
    }
  }
  //'nearly' Gaussian fluctuation if a1>nmaxCont&&a2>nmaxCont&&a3>nmaxCont
  G4double emean = 0.;
  G4double sig2e = 0.;

  // excitation of type 1
  if (a1 > 0.0) {
    AddExcitation2(rndmEngineF, a1, e1, emean, loss, sig2e);
  }

  // excitation of type 2
  if (a2 > 0.0) {
    AddExcitation2(rndmEngineF, a2, e2, emean, loss, sig2e);
  }

  if (sig2e > 0.0) {
    SampleGauss2(rndmEngineF, emean, sig2e, loss);
  }

  // ionisation
  if (a3 > 0.) {
    emean = 0.;
    sig2e = 0.;
    G4double p3 = a3;
    G4double alfa = 1.;
    if (a3 > nmaxCont) {
      alfa = w1 * (nmaxCont + a3) / (w1 * nmaxCont + a3);
      G4double alfa1 = alfa * G4Log(alfa) / (alfa - 1.);
      G4double namean = a3 * w1 * (alfa - 1.) / ((w1 - 1.) * alfa);
      emean += namean * e0 * alfa1;
      sig2e += e0 * e0 * namean * (alfa - alfa1 * alfa1);
      p3 = a3 - namean;
    }

    G4double w2 = alfa * e0;
    if (tmax > w2) {
      G4double w = (tmax - w2) / tmax;
      G4int nnb = G4Poisson(p3);
      if (nnb > 0) {
        if (nnb > sizearray) {
          sizearray = nnb;
          delete[] rndmarray;
          rndmarray = new G4double[nnb];
        }
        rndmEngineF->flatArray(nnb, rndmarray);
        for (G4int k = 0; k < nnb; ++k) {
          loss += w2 / (1. - w * rndmarray[k]);
        }
      }
    }
    if (sig2e > 0.0) {
      SampleGauss2(rndmEngineF, emean, sig2e, loss);
    }
  }

  loss *= scaling;

  return loss + extraloss;
}

G4double G4UniversalFluctuationForExtrapolator::SampleFluctuations(
    const G4Material* material, const G4DynamicParticle* dp, G4double tmax, G4double length, G4double ekin) {
  recordValid_ = false;
  // Calculate actual loss from the mean loss.
  // The model used to get the fluctuations is essentially the same
  // as in Glandz in Geant3 (Cern program library W5013, phys332).
  // L. Urban et al. NIM A362, p.416 (1995) and Geant4 Physics Reference Manual

  // shortcut for very small loss or from a step nearly equal to the range
  // (out of validity of the model)
  //

  // Full rebind (not just InitialiseMe) BEFORE the dedx lookup: table /
  // massratio / charge2ratio must match the particle, otherwise a
  // particleNameOverride (kaon/pion/proton) would compute the variance
  // with the muon dE/dx table while the mean loss uses proton scaling.
  if (dp->GetDefinition() != particle) {
    const G4double q = dp->GetDefinition()->GetPDGCharge() / CLHEP::eplus;
    SetParticleAndCharge(dp->GetDefinition(), q * q);
  }

  size_t idx = 0;
  G4double dedx = ((*table)[material->GetIndex()])->Value(massratio * ekin, idx) * charge2ratio;
  // SPECIES-DEPENDENT Tmax (CVH_REF_SPECIESDEDX). The SAME cvhcgf helper the
  // reference's ComputeDEDX uses, so the two cannot drift apart; exactly +0.0
  // for a proton/antiproton and never reached by a muon.
  if (speciesDedx && usesScaledProtonTable) {
    dedx += cvhcgf::speciesTmaxDedx(ekin, particleMass, speciesRefMass, chargeSquare, material->GetElectronDensity());
  }
  G4double meanLoss = length * dedx;

  G4double tkin = ekin;
  if (meanLoss < minLoss) {
    return 0.;
  }

  G4double tau = tkin * m_Inv_particleMass;
  G4double gam = tau + 1.0;
  G4double gam2 = gam * gam;
  G4double beta2 = tau * (tau + 2.0) / gam2;

  G4double loss(0.), siga(0.);

  // Gaussian regime
  // for heavy particles only and conditions
  // for Gauusian fluct. has been changed
  //
  if ((particleMass > electron_mass_c2) && (meanLoss >= minNumberInteractionsBohr * tmax)) {
    G4double tmaxkine = 2. * electron_mass_c2 * beta2 * gam2 / (1. + m_massrate * (2. * gam + m_massrate));
    if (tmaxkine <= 2. * tmax) {
      electronDensity = material->GetElectronDensity();
      siga = sqrt((1.0 / beta2 - 0.5) * twopi_mc2_rcl2 * tmax * length * electronDensity * chargeSquare);

      G4double sn = meanLoss / siga;

      // thick target case
      if (sn >= 2.0) {
        record_ = UrbanFluctRecord{};
        record_.regime = 0;
        record_.gsig2 = siga * siga;
        recordValid_ = true;
        return siga * siga;

        G4double twomeanLoss = meanLoss + meanLoss;
        do {
          loss = meanLoss;
          // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
        } while (0.0 > loss || twomeanLoss < loss);

        // Gamma distribution
      } else {
        G4double neff = sn * sn;
        record_ = UrbanFluctRecord{};
        record_.regime = 0;
        record_.gsig2 = meanLoss * meanLoss / neff;
        recordValid_ = true;
        return meanLoss * meanLoss / neff;
        loss = meanLoss;
      }
      return loss;
    }
  }

  // Glandz regime : initialisation
  //
  if (material != lastMaterial) {
    f1Fluct = material->GetIonisation()->GetF1fluct();
    f2Fluct = material->GetIonisation()->GetF2fluct();
    e1Fluct = material->GetIonisation()->GetEnergy1fluct();
    e2Fluct = material->GetIonisation()->GetEnergy2fluct();
    e1LogFluct = material->GetIonisation()->GetLogEnergy1fluct();
    e2LogFluct = material->GetIonisation()->GetLogEnergy2fluct();
    ipotFluct = material->GetIonisation()->GetMeanExcitationEnergy();
    ipotLogFluct = material->GetIonisation()->GetLogMeanExcEnergy();
    e0 = material->GetIonisation()->GetEnergy0fluct();
    esmall = 0.5 * sqrt(e0 * ipotFluct);
    lastMaterial = material;
  }

  // very small step or low-density material
  if (tmax <= e0) {
    return 0.;
  }

  // width correction for small cuts
  G4double scaling = std::min(1. + 0.5 * CLHEP::keV / tmax, 1.50);
  meanLoss /= scaling;

  G4double a1(0.0), a2(0.0), a3(0.0);

  loss = 0.0;

  e1 = e1Fluct;
  e2 = e2Fluct;

  // ------------------------------------------------------------------------
  // CVH_IONI_URBAN2021 -- harmonize the EXCITATION channels with the model the
  // simulation actually samples.
  //
  // This class is the PRE-2021 Urban model: two excitation channels at
  // e1Fluct and e2Fluct = 10 Zeff^2 eV with weights f1 = 1 - 2/Zeff,
  // f2 = 2/Zeff, and the fwnow floor 0.5.  Stock Geant4 11.2.2 -- which is
  // what `g4SimHits` runs, i.e. what the clean-propagation "data" is drawn
  // from -- is the 2021 model: ONE excitation channel at the mean excitation
  // energy ipotFluct, and the fwnow floor 0.1 (G4UniversalFluctuation.hh in
  // 11.2.2 has no f1Fluct/f2Fluct/e1Fluct/e2Fluct members at all).
  // NOTES_URBANSAMPLING measured the gap: +28 % excitation kappa2, +1.1 % of a
  // plane step's core variance, validated against 2.5e7 real stock samples.
  //
  // Two implementations of the same physics differing by 28 % in a channel is
  // a defect whether or not it moves the closure, so this brings the
  // extrapolator FORWARD to the simulation (stock Geant4 cannot be changed).
  // The body below is the 2021 Glandz initialisation verbatim -- it is the
  // same code as this class's own (dead) SampleGlandz, which is a copy of the
  // 2021 model.
  //
  // a2 is set to ZERO and e2 is LEFT AT e2Fluct: the offline consumers all
  // guard on `aj > 0`, so a zero-weight channel is inert, while keeping e2
  // preserves the material identification (e2 = 10 Zeff^2 eV) that the
  // downstream analysis uses.  The record stride is unchanged.
  //
  // DEFAULT OFF; composes with CVH_IONI_EXACTDELTA, which rescales a1 and a2
  // by a common factor and therefore does the right thing when a2 = 0.
  if (urban2021Enabled()) {
    e2 = e2Fluct;
    a2 = 0.0;
    e1 = ipotFluct;
    if (tmax > e1) {
      a1 = meanLoss * (1. - rate) / e1;
      if (a1 < a0) {
        const G4double fwnow = 0.1 + (fw - 0.1) * std::sqrt(a1 / a0);
        a1 /= fwnow;
        e1 *= fwnow;
      } else {
        a1 /= fw;
        e1 *= fw;
      }
    }
  } else if (tmax > ipotFluct) {
    G4double w2 = G4Log(2. * electron_mass_c2 * beta2 * gam2) - beta2;

    if (w2 > ipotLogFluct) {
      if (w2 > e2LogFluct) {
        G4double C = meanLoss * (1. - rate) / (w2 - ipotLogFluct);
        a1 = C * f1Fluct * (w2 - e1LogFluct) / e1Fluct;
        a2 = C * f2Fluct * (w2 - e2LogFluct) / e2Fluct;
      } else {
        a1 = meanLoss * (1. - rate) / e1;
      }
      if (a1 < a0) {
        G4double fwnow = 0.5 + (fw - 0.5) * sqrt(a1 / a0);
        a1 /= fwnow;
        e1 *= fwnow;
      } else {
        a1 /= fw;
        e1 = fw * e1Fluct;
      }
    }
  }
  // ------------------------------------------------------------------------

  G4double w1 = tmax / e0;
  if (tmax > e0) {
    a3 = rate * meanLoss * (tmax - e0) / (e0 * tmax * G4Log(w1));
    if (a1 + a2 <= 0.) {
      a3 /= rate;
    }
  }

  // ------------------------------------------------------------------------
  // CVH_IONI_EXACTDELTA -- replace the a3 channel's pure 1/E^2 delta spectrum
  // by the exact spin-1/2 knock-on cross section, with its own (parameter
  // free) normalization.
  //
  // WHY.  The Urban a3 channel is the model's ENTIRE representation of hard
  // delta rays, and it is wrong in two independent ways:
  //
  //  (1) NORMALIZATION.  a3 C = rate * meanLoss / ln(tmax/e0) with rate = 0.56
  //      is an Urban parameterization of the excitation/ionization split.  The
  //      physical normalization is xi = 2 pi re^2 me c^2 n_e z^2 L / beta^2,
  //      the same quantity this class already uses in its regime-0 variance.
  //      Their ratio is ln(tmax/e0) / (rate * L_B) with L_B = meanLoss/xi the
  //      Bethe stopping number -- 1.285 for the layered toy at pT = 3 and
  //      1.329 at pT = 40, i.e. the a3 channel is ~28 % short and the shortfall
  //      is MOMENTUM DEPENDENT.
  //
  //  (2) SHAPE.  The exact cross section for a spin-1/2 projectile on a free
  //      electron (PDG "Passage of particles through matter" eq. 34.7; Geant4
  //      PRM "Muon ionisation" = G4MuBetheBlochModel) is
  //
  //          d sigma/dT ~ (1/T^2) [ 1 - beta^2 T/Tmax + T^2/(2 E^2) ]
  //
  //      i.e. pure 1/T^2 times a suppression that reaches 1 - beta^2 +
  //      Tmax^2/2E^2 = 0.025 at T = Tmax for a 3 GeV muon and 0.31 for a
  //      42 GeV one.  That is why the model's measured ENERGY deficit
  //      (1.19-1.26) is smaller than its RATE deficit (1.28-1.33), and why a
  //      single a3 scale cannot describe both momenta: the end-point behaviour
  //      follows Tmax/E, which runs 0.22 -> 0.79 between them.
  //
  // Both are fixed here.  The excitation channels are then rescaled (a1, a2
  // multiplied by a common factor, so their split is untouched) to hold the
  // block's MEAN loss at the dE/dx table value -- the mean is not in question
  // and must not move.
  //
  // WHY HERE AND NOT IN THE OFFLINE CF.  xi needs the material's electron
  // density and the step length.  The exported record carries neither, and it
  // cannot be reconstructed from the record: xi = meanLoss/L_B and L_B needs
  // the density-effect correction delta, which is a per-material Geant4 table.
  // The offline side can therefore apply the SHAPE correction but not the
  // NORMALIZATION one, and the normalization is the larger of the two.  Doing
  // it here also fixes the variance this function RETURNS, i.e. the track
  // fit's own Q matrix, which the offline route leaves wrong.
  //
  // With the switch off not one line below runs and the record is bit-identical
  // to the 1/E^2 delta channel's.
  //
  // `IoniExactDeltaT0` (MeV, 0 = use e0) raises the bottom of the exact
  // channel.  The free-electron cross section is not valid below ~1 keV, where
  // atomic binding matters and the Urban excitation channels are the stand-in;
  // moving T0 from e0 = 10 eV up to 1 keV moves ~12 % of the mean between the
  // two channels and is the systematic on this construction.  It is a scan
  // knob, not a tune: nothing is fitted to it.
  const bool exactDelta = exactDeltaEnabled();
  // ONE reader, shared with the offline consumer -- cvhcgf::ioniKokoulinEnabled.
  const bool kokoulinOn = cvhcgf::ioniKokoulinEnabled();
  const double exactT0 = cvhcgf::switches().ioniExactDeltaT0;
  // Geant4 branches the spectrum on the PDG SPIN, not on the particle id:
  // G4BetheBlochModel::SampleSecondaries adds 0.5 T^2/Etot^2 to the rejection
  // function `if (0.5 == spin)` and G4MuBetheBlochModel (spin 1/2 by
  // construction) always has it.  Mirror that, so the model matches the
  // simulation by construction rather than by our own physics judgement.
  const G4double pdgSpin = (particle != nullptr) ? particle->GetPDGSpin() : 0.5;
  const G4bool spinHalf = (pdgSpin == 0.5);
  G4bool useExact = false;
  G4double xiEx = 0., t0Ex = 0.;
  if (exactDelta && a3 > 0. && (a1 + a2) > 0.) {
    // pre-`scaling` units throughout, exactly like meanLoss / e1 / e2 / tmax:
    // the sampled loss is multiplied by `scaling` at the end, so xi -- which
    // scales like an energy -- is divided by it here.
    const G4double t0 = (exactT0 > e0) ? exactT0 : e0;
    if (tmax > t0) {
      const G4double etot = ekin + particleMass;
      const G4double xi =
          twopi_mc2_rcl2 * material->GetElectronDensity() * length * chargeSquare / (beta2 * scaling);
      // mean of the exact channel on [t0, tmax]:
      //   int T dN/dT = xi [ ln(tmax/t0) - beta^2 (tmax-t0)/tmax
      //                      + (tmax^2 - t0^2)/(4 E^2) ]
      G4double i1 = G4Log(tmax / t0) - beta2 * (tmax - t0) / tmax;
      if (spinHalf) {
        i1 += (tmax * tmax - t0 * t0) / (4. * etot * etot);
      }
      const G4double eDelta = xi * i1;
      const G4double eExc = a1 * e1 + a2 * e2;
      if (eDelta > 0. && eDelta < meanLoss && eExc > 0.) {
        const G4double fexc = (meanLoss - eDelta) / eExc;
        a1 *= fexc;
        a2 *= fexc;
        useExact = true;
        xiEx = xi;
        t0Ex = t0;
      }
    }
  }
  // ------------------------------------------------------------------------

  // EVERYTHING THAT USED TO BE HERE IS GONE (2026-08-24), WITH THE TRUNCATION.
  //
  // This entry point returns a VARIANCE and nothing else: `SampleFluctuations2`
  // keeps its own copy of the Urban logic for the sampling path, and the mean
  // computed here was never read. The variance is now the recorded model's own
  // second cumulant, formed once at the return statement, so the accumulation
  // that used to stand here -- `AddExcitation` into emean/sig2e, `SampleGauss`,
  // and the two ALPHA-TRUNCATED ionization branches (the 1/E^2 channel's
  // `ualpha` split and the exact channel's `ta` quantile, plus its Kokoulin
  // integral over the truncated range) -- computes a number nothing reads.
  //
  // Deleting it is not a tidy-up. It removes the SECOND formula for this
  // block's variance, and the thing that used to reconcile the two formulas
  // was the truncation quantile itself.


  // record the underlying (untruncated) Urban model of this step: Poisson
  // excitations (a1,e1), (a2,e2) plus a3 delta collisions on the 1/E^2
  // spectrum [e0, tmax]. The alfa/namean split and the ualpha truncation
  // above are numerical evaluation devices for the returned variance, NOT
  // part of the model, so they are deliberately not recorded.
  record_ = UrbanFluctRecord{};
  record_.regime = useExact ? (spinHalf ? 2 : 3) : 1;
  // gsig2 in the Glandz regime = the RETURNED (alpha-truncated) variance:
  // the offline CF fit uses it to reproduce the exact standardization the
  // track fit applied, so sigma-replica errors cannot leak into the fitted
  // material scale.
  // Filled from the cumulant just below: this is the column an offline
  // consumer uses to reproduce the variance the fit saw, so it must BE that
  // variance.
  record_.gsig2 = 0.;
  record_.a1 = a1;
  record_.e1 = e1;
  record_.a2 = a2;
  record_.e2 = e2;
  // regime 2 reuses the a3 / e0r slots with the meanings documented on
  // UrbanFluctRecord: a3 -> xi (an ENERGY, scaled by `scaling` downstream like
  // every other energy in the record), e0r -> the bottom of the exact channel.
  record_.a3 = useExact ? xiEx : a3;
  record_.e0r = useExact ? t0Ex : e0;
  record_.tmaxr = tmax;
  record_.scaling = scaling;
  // beta^2 and the total energy, PHYSICAL and unscaled. They are exported only
  // in regime 2/3 and they are not optional there: beta^2 and E cannot be
  // recovered from tmaxr alone without the particle mass, and assuming a muon
  // would put a 3 GeV kaon's beta^2 at 0.99999 instead of 0.976.
  record_.beta2 = beta2;
  record_.etot = ekin + particleMass;
  recordValid_ = true;

  // THE RETURNED VARIANCE IS THE RECORD'S OWN SECOND CUMULANT.
  //
  // There used to be two formulas for this block's variance -- the truncated
  // integral built above, and `cvhcgf::blockKappa2`, which is what the fit
  // normalizes the block CGF by -- and an alpha to reconcile them. Now there
  // is one, it is untruncated, and the two cannot disagree.
  //
  // The Kokoulin excess uses the FLUCTUATION model's bucket count
  // (`ioniKokoulinNbin`, default 96), not the block CF's: this is a
  // once-per-step closed-form quadrature and there is nothing to economize.
  //
  // IT IS NOT A GOOD GAUSSIAN WEIGHT and is not meant to be: for a 1/E^2
  // spectrum the second moment is dominated by the hard edge. That is exactly
  // why the fit's weight is the block's Fisher information (CgfQoPMode = 1,
  // the default) and why mode 0 is a diagnostic limit rather than a
  // configuration.
  const G4double k2 = blockKappa2Of(record_, cvhcgf::ioniKokoulinEnabled() ? cvhcgf::ioniKokoulinNbin() : 0);
  record_.gsig2 = k2;
  return k2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4UniversalFluctuationForExtrapolator::SampleGlandz(CLHEP::HepRandomEngine* rndmEngineF,
                                                             const G4Material*,
                                                             const G4double tcut) {
  G4double a1(0.0), a3(0.0);
  G4double loss = 0.0;
  G4double e1 = ipotFluct;

  if (tcut > e1) {
    a1 = meanLoss * (1. - rate) / e1;
    if (a1 < a0) {
      const G4double fwnow = 0.1 + (fw - 0.1) * std::sqrt(a1 / a0);
      a1 /= fwnow;
      e1 *= fwnow;
    } else {
      a1 /= fw;
      e1 *= fw;
    }
  }

  const G4double w1 = tcut / e0;
  a3 = rate * meanLoss * (tcut - e0) / (e0 * tcut * G4Log(w1));
  if (a1 <= 0.) {
    a3 /= rate;
  }

  //'nearly' Gaussian fluctuation if a1>nmaxCont&&a2>nmaxCont&&a3>nmaxCont
  G4double emean = 0.;
  G4double sig2e = 0.;

  // excitation of type 1
  if (a1 > 0.0) {
    AddExcitation2(rndmEngineF, a1, e1, emean, loss, sig2e);
  }

  if (sig2e > 0.0) {
    SampleGauss2(rndmEngineF, emean, sig2e, loss);
  }

  // ionisation
  if (a3 > 0.) {
    emean = 0.;
    sig2e = 0.;
    G4double p3 = a3;
    G4double alfa = 1.;
    if (a3 > nmaxCont) {
      alfa = w1 * (nmaxCont + a3) / (w1 * nmaxCont + a3);
      const G4double alfa1 = alfa * G4Log(alfa) / (alfa - 1.);
      const G4double namean = a3 * w1 * (alfa - 1.) / ((w1 - 1.) * alfa);
      emean += namean * e0 * alfa1;
      sig2e += e0 * e0 * namean * (alfa - alfa1 * alfa1);
      p3 = a3 - namean;
    }

    const G4double w3 = alfa * e0;
    if (tcut > w3) {
      const G4double w = (tcut - w3) / tcut;
      const G4int nnb = (G4int)G4Poisson(p3);
      if (nnb > 0) {
        if (nnb > sizearray) {
          sizearray = nnb;
          delete[] rndmarray;
          rndmarray = new G4double[nnb];
        }
        rndmEngineF->flatArray(nnb, rndmarray);
        for (G4int k = 0; k < nnb; ++k) {
          loss += w3 / (1. - w * rndmarray[k]);
        }
      }
    }
    if (sig2e > 0.0) {
      SampleGauss2(rndmEngineF, emean, sig2e, loss);
    }
  }
  //G4cout << "### loss=" << loss << G4endl;
  return loss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4UniversalFluctuationForExtrapolator::Dispersion(const G4Material* material,
                                                           const G4DynamicParticle* dp,
                                                           G4double tmax,
                                                           G4double length) {
  if (dp->GetDefinition() != particle) {
    InitialiseMe(dp->GetDefinition());
  }

  electronDensity = material->GetElectronDensity();

  G4double gam = (dp->GetKineticEnergy()) * m_Inv_particleMass + 1.0;
  G4double beta2 = 1.0 - 1.0 / (gam * gam);

  G4double siga = (1.0 / beta2 - 0.5) * twopi_mc2_rcl2 * tmax * length * electronDensity * chargeSquare;

  return siga;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4UniversalFluctuationForExtrapolator::SetParticleAndCharge(const G4ParticleDefinition* part, G4double q2) {
  if (part != particle) {
    particle = part;
    particleMass = part->GetPDGMass();

    // Derived quantities
    m_Inv_particleMass = 1.0 / particleMass;
    m_massrate = CLHEP::electron_mass_c2 * m_Inv_particleMass;
  }
  chargeSquare = q2;

  // Latched here rather than read per step; same single reader as the
  // extrapolator's, so the reference's mean and this mean cannot disagree
  // about whether the correction is on.
  speciesDedx = cvhcgf::referenceIsSpeciesDedx();
  // CHARGE-AWARE mean loss (cvhcgf::referenceIsChargeAware, CVH_REF_CHARGEAWARE).
  //
  // Until 2026-08-16 this function selected fDedxMuon / fDedxProton
  // UNCONDITIONALLY, with no isNegative dispatch, while
  // G4EnergyLossForExtrapolatorForCVH::ComputeDEDX -- which builds the
  // REFERENCE trajectory -- did dispatch. So with the switch on, the
  // reference's mean and this class's `meanLoss = length * dedx` disagreed by
  // exactly the charge-odd 3.2e-3 of Geant4's high-order stopping-power block
  // (2*(Barkas + Bloch) + Mott) on every NEGATIVE track. Diagnosed in
  // NOTES_SPECIESDEDX s2.1/s8, closed here.
  //
  // It is the same class of inconsistency CVH_REF_SPECIESDEDX was careful to
  // avoid, and it matters for the same reason: `meanLoss` sets the Urban
  // channel weights (a1, a2, a3), the excitation rescale that pins the block
  // to the mean, and in the Gaussian regime the returned loss itself. A noise
  // model centred on a different mean than the trajectory it is the noise of
  // is not a smaller error for being second order -- it is a different error.
  chargeAware = cvhcgf::referenceIsChargeAware();
  const bool negative = chargeAware && part->GetPDGCharge() < 0.0;
  usesScaledProtonTable = false;
  speciesRefMass = 0.;

  if (part == G4Electron::Electron()) {
    // electron and positron already have their OWN tables, built from their
    // own particle -- charge-awareness was never foreign to this class, the
    // muon and hadron branches were the omission.
    table = tables->GetPhysicsTable(fDedxElectron);
    massratio = 1.;
    charge2ratio = 1.;
  } else if (part == G4Positron::Positron()) {
    table = tables->GetPhysicsTable(fDedxPositron);
    massratio = 1.;
    charge2ratio = 1.;
  } else if (part == G4MuonPlus::MuonPlus() || part == G4MuonMinus::MuonMinus()) {
    table = tables->GetPhysicsTable(negative ? fDedxMuonMinus : fDedxMuon);
    massratio = 1.;
    charge2ratio = 1.;
  } else {
    // scaled energy loss from proton tables
    table = tables->GetPhysicsTable(negative ? fDedxAntiProton : fDedxProton);
    massratio = proton_mass_c2 / particleMass;
    charge2ratio = part->GetPDGCharge() * part->GetPDGCharge();
    usesScaledProtonTable = true;
    // The mass of the particle this table was BUILT from, not
    // CLHEP::proton_mass_c2 -- see cvhcgf::speciesTmaxDedx. It has to follow
    // the table actually selected one line above, exactly as
    // G4EnergyLossForExtrapolatorForCVH::tableParticleMass does, or the
    // antiproton would lose its exact null under the two switches together.
    speciesRefMass =
        negative ? G4AntiProton::AntiProton()->GetPDGMass() : G4Proton::Proton()->GetPDGMass();
  }

  // A NULL TABLE IS NOT ALLOWED TO BE QUIET, on either branch.
  //
  // ComputeValue-style lookups on a null table give 0, i.e. a noise model with
  // no mean loss at all -- the silent-weight failure mode this study has hit
  // before (NOTES_QVALID s13). With chargeAware on, a half-built table set
  // would do exactly that for every negative track and the fit would still
  // converge. The extrapolator throws for the same reason in its own
  // Initialisation; this is the fluctuation's copy of that guard.
  if ((chargeAware || speciesDedx) && (nullptr == table || table->empty())) {
    throw cms::Exception("G4UniversalFluctuationForExtrapolator")
        << "the dE/dx table for " << part->GetParticleName() << " is "
        << (nullptr == table ? "missing" : "empty") << " (CVH_REF_CHARGEAWARE="
        << (chargeAware ? 1 : 0) << ", CVH_REF_SPECIESDEDX=" << (speciesDedx ? 1 : 0) << ")";
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
