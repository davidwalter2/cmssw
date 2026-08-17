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
//---------------------------------------------------------------------------
//
// ClassName:    G4EnergyLossForExtrapolatorForCVH
//
// Description:  This class provide calculation of energy loss, fluctuation,
//               and msc angle
//
// Author:       09.12.04 V.Ivanchenko
//
// Modification:
// 08-04-05 Rename Propogator -> Extrapolator (V.Ivanchenko)
// 16-03-06 Add muon tables and fix bug in units (V.Ivanchenko)
// 21-03-06 Add verbosity defined in the constructor and Initialisation
//          start only when first public method is called (V.Ivanchenko)
// 03-05-06 Remove unused pointer G4Material* from number of methods (VI)
// 12-05-06 SEt linLossLimit=0.001 (VI)
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "TrackPropagation/Geant4e/interface/G4EnergyLossForExtrapolatorForCVH.h"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4ParticleTable.hh"
#include "TrackPropagation/Geant4e/interface/CGFQoPBlock.h"
#include "FWCore/Utilities/interface/Exception.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifdef G4MULTITHREADED
G4Mutex G4EnergyLossForExtrapolatorForCVH::extrMutex = G4MUTEX_INITIALIZER;
#endif

G4TablesForExtrapolatorForCVH* G4EnergyLossForExtrapolatorForCVH::tables = nullptr;

G4EnergyLossForExtrapolatorForCVH::G4EnergyLossForExtrapolatorForCVH(G4int verb)
    : maxEnergyTransfer(DBL_MAX), verbose(verb) {
  emin = 1. * CLHEP::MeV;
  emax = 100. * CLHEP::TeV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EnergyLossForExtrapolatorForCVH::~G4EnergyLossForExtrapolatorForCVH() {
  if (isMaster) {
    delete tables;
    tables = nullptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolatorForCVH::EnergyAfterStep(G4double kinEnergy,
                                                            G4double stepLength,
                                                            const G4Material* mat,
                                                            const G4ParticleDefinition* part) {
  G4double kinEnergyFinal = kinEnergy;
  if (SetupKinematics(part, mat, kinEnergy)) {
    G4double step = TrueStepLength(kinEnergy, stepLength, mat, part);
    G4double r = ComputeRange(kinEnergy, part, mat);
    if (r <= step) {
      kinEnergyFinal = 0.0;
    } else if (step < linLossLimit * r) {
      kinEnergyFinal -= step * ComputeDEDX(kinEnergy, part, mat);
    } else {
      G4double r1 = r - step;
      kinEnergyFinal = ComputeEnergy(r1, part, mat);
    }
  }
  return kinEnergyFinal;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolatorForCVH::EnergyBeforeStep(G4double kinEnergy,
                                                             G4double stepLength,
                                                             const G4Material* mat,
                                                             const G4ParticleDefinition* part) {
  G4double kinEnergyFinal = kinEnergy;

  if (SetupKinematics(part, mat, kinEnergy)) {
    G4double step = TrueStepLength(kinEnergy, stepLength, mat, part);
    G4double r = ComputeRange(kinEnergy, part, mat);

    if (step < linLossLimit * r) {
      kinEnergyFinal += step * ComputeDEDX(kinEnergy, part, mat);
    } else {
      G4double r1 = r + step;
      kinEnergyFinal = ComputeEnergy(r1, part, mat);
    }
  }
  return kinEnergyFinal;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolatorForCVH::TrueStepLength(G4double kinEnergy,
                                                           G4double stepLength,
                                                           const G4Material* mat,
                                                           const G4ParticleDefinition* part) {
  G4double res = stepLength;
  if (SetupKinematics(part, mat, kinEnergy)) {
    if (part == electron || part == positron) {
      const G4double x = stepLength * ComputeValue(kinEnergy, GetPhysicsTable(fMscElectron), mat->GetIndex());
      if (x < 0.2) {
        res *= (1.0 + 0.5 * x + x * x / 3.0);
      } else if (x < 0.9999) {
        res = -G4Log(1.0 - x) * stepLength / x;
      } else {
        res = ComputeRange(kinEnergy, part, mat);
      }
    } else {
      res = ComputeTrueStep(mat, part, kinEnergy, stepLength);
    }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4EnergyLossForExtrapolatorForCVH::SetupKinematics(const G4ParticleDefinition* part,
                                                          const G4Material* mat,
                                                          G4double kinEnergy) {
  if (mat->GetNumberOfMaterials() != nmat) {
    Initialisation();
  }
  if (nullptr == part || nullptr == mat || kinEnergy < CLHEP::keV) {
    return false;
  }
  if (part != currentParticle) {
    currentParticle = part;
    G4double q = part->GetPDGCharge() / eplus;
    charge2 = q * q;
  }
  if (mat != currentMaterial) {
    size_t i = mat->GetIndex();
    if (i >= nmat) {
      G4cout << "### G4EnergyLossForExtrapolatorForCVH WARNING: material index i= " << i
             << " above number of materials " << nmat << G4endl;
      return false;
    } else {
      currentMaterial = mat;
      electronDensity = mat->GetElectronDensity();
      radLength = mat->GetRadlen();
    }
  }
  if (kinEnergy != kineticEnergy) {
    kineticEnergy = kinEnergy;
    G4double mass = part->GetPDGMass();
    G4double tau = kinEnergy / mass;

    gam = tau + 1.0;
    bg2 = tau * (tau + 2.0);
    beta2 = bg2 / (gam * gam);
    tmax = kinEnergy;
    if (part == electron)
      tmax *= 0.5;
    else if (part != positron) {
      G4double r = CLHEP::electron_mass_c2 / mass;
      tmax = 2.0 * bg2 * CLHEP::electron_mass_c2 / (1.0 + 2.0 * gam * r + r * r);
    }
    tmax = std::min(tmax, maxEnergyTransfer);
  }
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4ParticleDefinition* G4EnergyLossForExtrapolatorForCVH::FindParticle(const G4String& name) {
  currentParticle = G4ParticleTable::GetParticleTable()->FindParticle(name);
  if (nullptr == currentParticle) {
    G4cout << "### G4EnergyLossForExtrapolatorForCVH WARNING: "
           << "FindParticle() fails to find " << name << G4endl;
  }
  return currentParticle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolatorForCVH::ComputeDEDX(G4double ekin,
                                                        const G4ParticleDefinition* part,
                                                        const G4Material* mat) {
  if (mat->GetNumberOfMaterials() != nmat) {
    Initialisation();
  }
  G4double x = 0.0;
  if (part == electron) {
    x = ComputeValue(ekin, GetPhysicsTable(fDedxElectron), mat->GetIndex());
  } else if (part == positron) {
    x = ComputeValue(ekin, GetPhysicsTable(fDedxPositron), mat->GetIndex());
  } else if (part == muonPlus || part == muonMinus) {
    x = ComputeValue(ekin, GetPhysicsTable(isNegative(part) ? fDedxMuonMinus : fDedxMuon), mat->GetIndex());
  } else {
    G4double e = ekin * CLHEP::proton_mass_c2 / part->GetPDGMass();
    G4double q = part->GetPDGCharge() / CLHEP::eplus;
    x = ComputeValue(e, GetPhysicsTable(isNegative(part) ? fDedxAntiProton : fDedxProton), mat->GetIndex()) * q * q;
    // SPECIES-DEPENDENT Tmax (CVH_REF_SPECIESDEDX). Exactly +0.0 for a proton
    // or an antiproton -- the table IS theirs -- and this branch is not
    // reached at all by a muon, so both are bit-for-bit nulls.
    x += totalDedxDelta(ekin, part, mat);
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolatorForCVH::speciesDedxDelta(G4double ekin,
                                                             const G4ParticleDefinition* part,
                                                             const G4Material* mat) const {
  if (!speciesDedx) {
    return 0.0;
  }
  const G4double q = part->GetPDGCharge() / CLHEP::eplus;
  const G4double d =
      cvhcgf::speciesTmaxDedx(ekin, part->GetPDGMass(), tableParticleMass(part), q * q, mat->GetElectronDensity());
  // Fail loud. A non-finite correction would propagate silently into the
  // reference trajectory as a NaN momentum and the fit would report a
  // convergence failure a long way from here.
  if (!std::isfinite(d)) {
    throw cms::Exception("G4EnergyLossForExtrapolatorForCVH")
        << "CVH_REF_SPECIESDEDX: non-finite dE/dx correction for " << part->GetParticleName() << " at ekin " << ekin
        << " MeV in " << mat->GetName();
  }
  return d;
}

G4double G4EnergyLossForExtrapolatorForCVH::radDedxDelta(G4double ekin,
                                                        const G4ParticleDefinition* part,
                                                        const G4Material* mat) {
  if (!hadronRad) {
    return 0.0;
  }
  const G4PhysicsTable* t = tables->GetHadronRadiativeTable(part);
  if (t == nullptr) {
    return 0.0;
  }
  // The table is built for THIS particle at its own kinetic energy, so ekin is
  // used directly -- no proton mass scaling. Radiative loss goes as q^2, like
  // ionization.
  const G4double q = part->GetPDGCharge() / CLHEP::eplus;
  const G4double d = ComputeValue(ekin, t, mat->GetIndex()) * q * q;
  if (!std::isfinite(d)) {
    throw cms::Exception("G4EnergyLossForExtrapolatorForCVH")
        << "CVH_REF_HADRAD: non-finite radiative dE/dx for " << part->GetParticleName() << " at ekin " << ekin
        << " MeV in " << mat->GetName();
  }
  // Radiative dE/dx is non-negative by definition and every tabulated value is
  // >= 0, so a negative here is log-spline undershoot, not physics. It is
  // routine in vacuum, where the true value is ~1e-27 MeV/mm and the
  // interpolation noise is the same order (measured: -3.0e-27 for a pi+ at
  // 30 MeV). Clamping is the correct reading; a genuine sign error would show
  // up as the closure moving the wrong way, which the physics validation
  // catches and a threshold here could not without also rejecting vacuum.
  return (d > 0.0) ? d : 0.0;
}

G4double G4EnergyLossForExtrapolatorForCVH::totalDedxDelta(G4double ekin,
                                                          const G4ParticleDefinition* part,
                                                          const G4Material* mat) {
  return speciesDedxDelta(ekin, part, mat) + radDedxDelta(ekin, part, mat);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolatorForCVH::speciesRangeDefect(G4double ekin,
                                                               const G4ParticleDefinition* part,
                                                               const G4Material* mat) {
  if (!speciesDedx || !(ekin > 0.0)) {
    return 0.0;
  }
  if (part == rdPart && mat == rdMat && ekin == rdEkin) {
    return rdValue;
  }
  const G4double m = part->GetPDGMass();
  const G4double massratio = CLHEP::proton_mass_c2 / m;
  const G4double q = part->GetPDGCharge() / CLHEP::eplus;
  const G4double q2 = q * q;
  const G4PhysicsTable* tab = GetPhysicsTable(isNegative(part) ? fDedxAntiProton : fDedxProton);
  const size_t idxMat = mat->GetIndex();

  // Composite Simpson in y = ln(1 + E'/m):  E' = m (e^y - 1),  dE' = (E'+m) dy.
  const G4int nb = cvhcgf::speciesDedxNbin();
  const G4double y1 = std::log1p(ekin / m);
  const G4double dy = y1 / nb;
  G4double acc = 0.0;
  for (G4int i = 0; i <= nb; ++i) {
    G4double f = 0.0;
    if (i > 0) {
      const G4double ep = m * std::expm1(i * dy);
      const G4double u = ComputeValue(ep * massratio, tab, idxMat) * q2;
      const G4double d = totalDedxDelta(ep, part, mat);
      const G4double ud = u + d;
      // u <= 0 means the table is missing or unpopulated at this energy; ud
      // <= 0 would mean the "correction" has eaten the whole stopping power.
      // Neither can happen for a real material (validated in Initialisation),
      // and contributing 0 rather than an infinity is the safe reading if the
      // very bottom of the grid ever produces one.
      if (u > 0.0 && ud > 0.0) {
        f = d / (u * ud) * (ep + m);
      }
    }
    const G4double c = (i == 0 || i == nb) ? 1.0 : ((i % 2) ? 4.0 : 2.0);
    acc += c * f;
  }
  const G4double res = acc * dy / 3.0;
  if (!std::isfinite(res)) {
    throw cms::Exception("G4EnergyLossForExtrapolatorForCVH")
        << "CVH_REF_SPECIESDEDX: non-finite range defect for " << part->GetParticleName() << " at ekin " << ekin
        << " MeV in " << mat->GetName();
  }
  rdPart = part;
  rdMat = mat;
  rdEkin = ekin;
  rdValue = res;
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolatorForCVH::ComputeRange(G4double ekin,
                                                         const G4ParticleDefinition* part,
                                                         const G4Material* mat) {
  if (mat->GetNumberOfMaterials() != nmat) {
    Initialisation();
  }
  G4double x = 0.0;
  if (part == electron) {
    x = ComputeValue(ekin, GetPhysicsTable(fRangeElectron), mat->GetIndex());
  } else if (part == positron) {
    x = ComputeValue(ekin, GetPhysicsTable(fRangePositron), mat->GetIndex());
  } else if (part == muonPlus || part == muonMinus) {
    x = ComputeValue(ekin, GetPhysicsTable(isNegative(part) ? fRangeMuonMinus : fRangeMuon), mat->GetIndex());
  } else {
    G4double massratio = CLHEP::proton_mass_c2 / part->GetPDGMass();
    G4double e = ekin * massratio;
    G4double q = part->GetPDGCharge() / CLHEP::eplus;
    x = ComputeValue(e, GetPhysicsTable(isNegative(part) ? fRangeAntiProton : fRangeProton), mat->GetIndex()) /
        (q * q * massratio);
    // R_corrected = R_uncorrected - D. Exactly +0.0 for a proton/antiproton
    // (every node's `d` is +0.0), so they are bit-for-bit nulls here too.
    if (speciesDedx) {
      const G4double d = speciesRangeDefect(ekin, part, mat);
      if (std::fabs(d) > 0.5 * x) {
        throw cms::Exception("G4EnergyLossForExtrapolatorForCVH")
            << "CVH_REF_SPECIESDEDX: range defect " << d << " mm is more than half the range " << x << " mm for "
            << part->GetParticleName() << " at ekin " << ekin << " MeV in " << mat->GetName()
            << " -- the perturbation treatment is not valid here";
      }
      x -= d;
    }
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolatorForCVH::ComputeEnergy(G4double range,
                                                          const G4ParticleDefinition* part,
                                                          const G4Material* mat) {
  if (mat->GetNumberOfMaterials() != nmat) {
    Initialisation();
  }
  G4double x = 0.0;
  if (part == electron) {
    x = ComputeValue(range, GetPhysicsTable(fInvRangeElectron), mat->GetIndex());
  } else if (part == positron) {
    x = ComputeValue(range, GetPhysicsTable(fInvRangePositron), mat->GetIndex());
  } else if (part == muonPlus || part == muonMinus) {
    x = ComputeValue(range, GetPhysicsTable(isNegative(part) ? fInvRangeMuonMinus : fInvRangeMuon), mat->GetIndex());
  } else {
    G4double massratio = CLHEP::proton_mass_c2 / part->GetPDGMass();
    G4double q = part->GetPDGCharge() / CLHEP::eplus;
    G4double r = range * massratio * q * q;
    const G4PhysicsTable* tab = GetPhysicsTable(isNegative(part) ? fInvRangeAntiProton : fInvRangeProton);
    x = ComputeValue(r, tab, mat->GetIndex()) / massratio;
    // ComputeEnergy must remain the NUMERICAL INVERSE of ComputeRange, or the
    // range branch of EnergyAfterStep loses a different energy from the dE/dx
    // branch by exactly the term this switch exists to remove. ComputeRange
    // returns R_0(E) - D(E), so solve R_0(E) = range + D(E) by fixed point.
    // The map contracts by |D'| dE/dr = |d|/(u+d) ~ 5e-3, so two iterations
    // leave 3e-5 of the correction, i.e. 1e-7 of the range.
    if (speciesDedx) {
      for (G4int it = 0; it < 2; ++it) {
        const G4double rr = (range + speciesRangeDefect(x, part, mat)) * massratio * q * q;
        if (!(rr > 0.0)) {
          break;
        }
        x = ComputeValue(rr, tab, mat->GetIndex()) / massratio;
      }
    }
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolatorForCVH::EnergyDispersion(G4double kinEnergy,
                                                             G4double stepLength,
                                                             const G4Material* mat,
                                                             const G4ParticleDefinition* part) {
  G4double sig2 = 0.0;
  if (SetupKinematics(part, mat, kinEnergy)) {
    G4double step = ComputeTrueStep(mat, part, kinEnergy, stepLength);
    sig2 = (1.0 / beta2 - 0.5) * CLHEP::twopi_mc2_rcl2 * tmax * step * electronDensity * charge2;
  }
  return sig2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolatorForCVH::AverageScatteringAngle(G4double kinEnergy,
                                                                   G4double stepLength,
                                                                   const G4Material* mat,
                                                                   const G4ParticleDefinition* part) {
  G4double theta = 0.0;
  if (SetupKinematics(part, mat, kinEnergy)) {
    G4double t = stepLength / radLength;
    G4double y = std::max(0.001, t);
    theta = 19.23 * CLHEP::MeV * std::sqrt(charge2 * t) * (1.0 + 0.038 * G4Log(y)) / (beta2 * gam * part->GetPDGMass());
  }
  return theta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EnergyLossForExtrapolatorForCVH::Initialisation() {
  if (verbose > 0) {
    G4cout << "### G4EnergyLossForExtrapolatorForCVH::Initialisation" << tables << G4endl;
  }
  electron = G4Electron::Electron();
  positron = G4Positron::Positron();
  proton = G4Proton::Proton();
  antiProton = G4AntiProton::AntiProton();
  muonPlus = G4MuonPlus::MuonPlus();
  muonMinus = G4MuonMinus::MuonMinus();
  // Same single reader the table build uses, so the dispatch below and the
  // set of tables that exist can never disagree.
  chargeAware = cvhcgf::referenceIsChargeAware();
  // Same pattern. The species correction builds no table of its own -- it is
  // added on top of the proton table's value -- so the latch here is the only
  // place it is read, and the memo below has to be dropped with it.
  speciesDedx = cvhcgf::referenceIsSpeciesDedx();
  hadronRad = cvhcgf::referenceHasHadronRadiative();
  rdPart = nullptr;
  rdMat = nullptr;
  rdEkin = -1.0;
  rdValue = 0.0;

  // initialisation for the 1st run
  if (nullptr == tables) {
#ifdef G4MULTITHREADED
    G4MUTEXLOCK(&extrMutex);
    if (nullptr == tables) {
#endif
      isMaster = true;
      // CVH_IONONLY=1 drops the RADIATIVE (brems + pair) mean from the
      // reference dE/dx table. Diagnostic for the mean-vs-mode energy-loss
      // bias found 2026-08-08: the reference subtracts the MEAN loss while the
      // typical muon loses the MODE, so the fit adds back energy never lost.
      // Radiative dE/dx ~ b*E is NEGLIGIBLE at low p (~0.2 MeV over the
      // tracker at 3 GeV) and only matters at high p, so this isolates the
      // flat/high-p term from the ionisation 1/p one.
      // Single reader, shared with the block CGF model, so the reference and
      // the noise model cannot disagree about the convention (see
      // cvhcgf::referenceIsIonOnly). Same value as the getenv it replaces.
      const bool _ionOnly = cvhcgf::referenceIsIonOnly();
      if (_ionOnly) {
        G4cout << "### G4EnergyLossForExtrapolatorForCVH: CVH_IONONLY set -- "
               << "radiative mean EXCLUDED from the dE/dx table" << G4endl;
      }
      tables = new G4TablesForExtrapolatorForCVH(verbose, nbins, emin, emax, _ionOnly);
      tables->Initialisation();
      if (cvhcgf::referenceIsChargeAware()) {
        // leading G4endl: the tables are built LAZILY, on the first
        // ComputeDEDX, so without it this banner lands in the middle of
        // whatever the caller was printing (it glued itself onto a
        // barkas_g4driver SPECIES line and broke its parser).
        G4cout << G4endl
               << "### G4EnergyLossForExtrapolatorForCVH: CVH_REF_CHARGEAWARE set -- "
               << "dE/dx, range and inverse-range tables built for G4MuonMinus and G4AntiProton "
               << "as well; the reference trajectory is no longer pinned to the positive particle" << G4endl;
      }
      if (speciesDedx) {
        // leading G4endl for the same reason as above: the tables are built
        // lazily, inside the first ComputeDEDX.
        G4cout << G4endl
               << "### G4EnergyLossForExtrapolatorForCVH: CVH_REF_SPECIESDEDX set -- "
               << "the hadron branch's proton-table dE/dx, range and inverse range are corrected by "
               << "xi ln(Tmax_species/Tmax_table); Simpson intervals for the range defect = "
               << cvhcgf::speciesDedxNbin() << G4endl;
      }
      nmat = G4Material::GetNumberOfMaterials();
      if (verbose > 0) {
        G4cout << "### G4EnergyLossForExtrapolator::BuildTables for " << nmat << " materials Nbins= " << nbins
               << " Emin(MeV)= " << emin << "  Emax(MeV)= " << emax << G4endl;
      }
#ifdef G4MULTITHREADED
    }
    G4MUTEXUNLOCK(&extrMutex);
#endif
  }

  // initialisation for the next run
  if (isMaster && G4Material::GetNumberOfMaterials() != nmat) {
#ifdef G4MULTITHREADED
    G4MUTEXLOCK(&extrMutex);
#endif
    tables->Initialisation();
#ifdef G4MULTITHREADED
    G4MUTEXUNLOCK(&extrMutex);
#endif
  }
  nmat = G4Material::GetNumberOfMaterials();

  // A missing or short table is NOT allowed to be quiet.
  //
  // ComputeValue returns 0.0 for a null table. With the switch OFF that gives
  // dE/dx = 0 and a reference that loses no energy -- visibly wrong. With it
  // ON the correction is ADDED to that zero, so the reference would integrate
  // the CORRECTION ALONE: a small, finite, plausible-looking stopping power,
  // and the fit would converge on it. The range defect would additionally
  // divide by it. This is the silent-weight failure mode this study has hit
  // before, so the switch refuses to run without the table it corrects.
  if (speciesDedx) {
    const ExtTableType want[] = {fDedxProton, fRangeProton, fInvRangeProton};
    for (ExtTableType t : want) {
      const G4PhysicsTable* p = GetPhysicsTable(t);
      if (nullptr == p || (G4int)p->length() < (G4int)nmat) {
        throw cms::Exception("G4EnergyLossForExtrapolatorForCVH")
            << "CVH_REF_SPECIESDEDX is set but the proton table it corrects is missing or short ("
            << (nullptr == p ? -1 : (G4int)p->length()) << " of " << nmat << " materials, table id " << (int)t << ")";
      }
    }
    if (chargeAware) {
      const ExtTableType wantNeg[] = {fDedxAntiProton, fRangeAntiProton, fInvRangeAntiProton};
      for (ExtTableType t : wantNeg) {
        const G4PhysicsTable* p = GetPhysicsTable(t);
        if (nullptr == p || (G4int)p->length() < (G4int)nmat) {
          throw cms::Exception("G4EnergyLossForExtrapolatorForCVH")
              << "CVH_REF_SPECIESDEDX with CVH_REF_CHARGEAWARE, but the antiproton table is missing or short ("
              << (nullptr == p ? -1 : (G4int)p->length()) << " of " << nmat << " materials, table id " << (int)t << ")";
        }
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
