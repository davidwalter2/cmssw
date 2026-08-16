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
// 08-04-05 Rename Propogator -> Extrapolator
// 16-03-06 Add muon tables
// 21-03-06 Add verbosity defined in the constructor and Initialisation
//          start only when first public method is called (V.Ivanchenko)
// 03-05-06 Remove unused pointer G4Material* from number of methods (VI)
// 28-07-07 Add maxEnergyTransfer for computation of energy loss (VI)
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef TrackPropagation_G4EnergyLossForExtrapolatorForCVH_h
#define TrackPropagation_G4EnergyLossForExtrapolatorForCVH_h 1

#include <vector>
#include <CLHEP/Units/PhysicalConstants.h>

#include "globals.hh"
#include "G4PhysicsTable.hh"
#include "TrackPropagation/Geant4e/interface/G4TablesForExtrapolatorForCVH.h"
#include "G4Log.hh"
#include "G4Threading.hh"
// needed complete, not forward-declared: isNegative() calls GetPDGCharge()
#include "G4ParticleDefinition.hh"

class G4Material;
class G4MaterialCutsCouple;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4EnergyLossForExtrapolatorForCVH {
public:
  explicit G4EnergyLossForExtrapolatorForCVH(G4int verb = 1);

  ~G4EnergyLossForExtrapolatorForCVH();

  void Initialisation();

  G4double ComputeDEDX(G4double kinEnergy, const G4ParticleDefinition*, const G4Material*);

  G4double ComputeRange(G4double kinEnergy, const G4ParticleDefinition*, const G4Material*);

  G4double ComputeEnergy(G4double range, const G4ParticleDefinition*, const G4Material*);

  G4double EnergyAfterStep(G4double kinEnergy, G4double step, const G4Material*, const G4ParticleDefinition*);

  G4double EnergyBeforeStep(G4double kinEnergy, G4double step, const G4Material*, const G4ParticleDefinition*);

  G4double TrueStepLength(G4double kinEnergy, G4double step, const G4Material*, const G4ParticleDefinition* part);

  inline G4double EnergyAfterStep(G4double kinEnergy, G4double step, const G4Material*, const G4String& particleName);

  inline G4double EnergyBeforeStep(G4double kinEnergy, G4double step, const G4Material*, const G4String& particleName);

  G4double AverageScatteringAngle(G4double kinEnergy,
                                  G4double step,
                                  const G4Material*,
                                  const G4ParticleDefinition* part);

  inline G4double AverageScatteringAngle(G4double kinEnergy,
                                         G4double step,
                                         const G4Material*,
                                         const G4String& particleName);

  inline G4double ComputeTrueStep(const G4Material*,
                                  const G4ParticleDefinition* part,
                                  G4double kinEnergy,
                                  G4double stepLength);

  G4double EnergyDispersion(G4double kinEnergy, G4double step, const G4Material*, const G4ParticleDefinition*);

  inline G4double EnergyDispersion(G4double kinEnergy, G4double step, const G4Material*, const G4String& particleName);

  inline void SetVerbose(G4int val);

  inline void SetMinKinEnergy(G4double);

  inline void SetMaxKinEnergy(G4double);

  inline void SetMaxEnergyTransfer(G4double);

  // hide assignment operator
  G4EnergyLossForExtrapolatorForCVH& operator=(const G4EnergyLossForExtrapolatorForCVH& right) = delete;
  G4EnergyLossForExtrapolatorForCVH(const G4EnergyLossForExtrapolatorForCVH&) = delete;

private:
  G4bool SetupKinematics(const G4ParticleDefinition*, const G4Material*, G4double kinEnergy);

  const G4ParticleDefinition* FindParticle(const G4String& name);

  inline G4double ComputeValue(G4double x, const G4PhysicsTable* table, size_t idxMat);

  inline const G4PhysicsTable* GetPhysicsTable(ExtTableType type) const;

  // "Should this track use the NEGATIVE particle's tables?"
  //
  // The charge comes from the track's own G4ParticleDefinition -- the
  // propagator builds it from the reconstructed charge
  // (Geant4ePropagator::generateParticleName, or the V0 ntuplizer's explicit
  // override) -- and never from a hard-coded PDG list, so mu/pi/K/p and every
  // one of their antiparticles, and any other singly-charged hadron the
  // propagator is ever asked for, are covered by the same two lines.
  //
  // Returns false unconditionally when the switch is off, which is what makes
  // the OFF path bit-identical rather than merely equivalent.
  inline G4bool isNegative(const G4ParticleDefinition* part) const;

  // SPECIES-DEPENDENT dE/dx (cvhcgf::referenceIsSpeciesDedx, CVH_REF_SPECIESDEDX).
  //
  // The PDG mass of the particle whose table the hadron branch is about to
  // read -- G4Proton, or G4AntiProton when the reference is also charge-aware.
  // The two corrections compose here and nowhere else: `isNegative` picks the
  // table, this picks the mass of the particle that table was built from, and
  // the Tmax correction is charge-EVEN so it is the same either way.
  //
  // It is deliberately NOT CLHEP::proton_mass_c2 (which is what the ENERGY
  // scaling above uses): G4Proton's PDG mass is 938.272013 MeV against CLHEP's
  // 938.27208816, and using the constant instead of the table particle's own
  // mass would give the proton a 1.4e-11 relative correction where the right
  // answer is an exact zero.
  inline G4double tableParticleMass(const G4ParticleDefinition* part) const;

  // The additive dE/dx correction, MeV/mm. Zero unless the switch is on AND
  // the particle is on the scaled-proton-table branch.
  G4double speciesDedxDelta(G4double ekin, const G4ParticleDefinition* part, const G4Material* mat) const;

  // D(E) = R_uncorrected(E) - R_corrected(E), the range the Tmax defect costs.
  //
  //     R(E) = INT_0^E dE'/dedx(E')   =>   D(E) = INT_0^E  d / (u (u + d)) dE'
  //
  // with u = the proton-table dE/dx and d = speciesDedxDelta. This is the
  // EXACT perturbation (no expansion in d/u): it is what makes ComputeRange
  // and ComputeEnergy describe the SAME corrected stopping power the dE/dx
  // branch uses, so that the two branches of EnergyAfterStep lose the same
  // energy over the same step.
  //
  // Composite Simpson in y = ln(1 + E'/m), which maps [0, E] onto [0, y1] with
  // the integrand's only structure -- the rise from 0 like beta^2 at low
  // beta*gamma -- resolved uniformly at any E. Not const: ComputeValue writes
  // the table-lookup hint `index`.
  G4double speciesRangeDefect(G4double ekin, const G4ParticleDefinition* part, const G4Material* mat);

#ifdef G4MULTITHREADED
  static G4Mutex extrMutex;
#endif
  static G4TablesForExtrapolatorForCVH* tables;

  const G4ParticleDefinition* currentParticle = nullptr;
  const G4ParticleDefinition* electron = nullptr;
  const G4ParticleDefinition* positron = nullptr;
  const G4ParticleDefinition* muonPlus = nullptr;
  const G4ParticleDefinition* muonMinus = nullptr;
  const G4ParticleDefinition* proton = nullptr;
  const G4ParticleDefinition* antiProton = nullptr;
  const G4Material* currentMaterial = nullptr;

  // True when the dE/dx / range / inverse-range tables also exist for the
  // NEGATIVE muon and the antiproton, i.e. when the reference trajectory is
  // charge-aware. Latched from cvhcgf::referenceIsChargeAware() at the same
  // moment the tables are built, so the dispatch below cannot ask for a table
  // that was not made. See CGFQoPBlock.h for the physics.
  G4bool chargeAware = false;

  // True when the hadron branch's proton-table lookup is corrected for the
  // species' own Tmax. Latched from cvhcgf::referenceIsSpeciesDedx() in
  // Initialisation, in the same place and for the same reason as chargeAware.
  G4bool speciesDedx = false;

  // One-entry memo of speciesRangeDefect. EnergyAfterStep and
  // GetContinuousStepLimit call ComputeRange with the SAME (particle,
  // material, ekin) within a step, so this halves the quadrature cost; it is
  // a pure memo of a pure function and holds no state of its own.
  const G4ParticleDefinition* rdPart = nullptr;
  const G4Material* rdMat = nullptr;
  G4double rdEkin = -1.0;
  G4double rdValue = 0.0;

  G4double electronDensity = 0.0;
  G4double radLength = 0.0;
  G4double charge2 = 0.0;
  G4double kineticEnergy = 0.0;
  G4double gam = 1.0;
  G4double bg2 = 0.0;
  G4double beta2 = 0.0;
  G4double tmax = 0.0;

  G4double linLossLimit = 0.01;
  G4double emin = 0.0;
  G4double emax = 0.0;
  G4double maxEnergyTransfer = 0.0;

  size_t index = 0;
  size_t nmat = 0;
  G4int nbins = 80;
  G4int verbose = 0;

  G4bool isMaster = false;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline const G4PhysicsTable* G4EnergyLossForExtrapolatorForCVH::GetPhysicsTable(ExtTableType type) const {
  return tables->GetPhysicsTable(type);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4EnergyLossForExtrapolatorForCVH::isNegative(const G4ParticleDefinition* part) const {
  return chargeAware && part->GetPDGCharge() < 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4EnergyLossForExtrapolatorForCVH::tableParticleMass(const G4ParticleDefinition* part) const {
  const G4ParticleDefinition* ref = isNegative(part) ? antiProton : proton;
  return (nullptr != ref) ? ref->GetPDGMass() : CLHEP::proton_mass_c2;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4EnergyLossForExtrapolatorForCVH::EnergyAfterStep(G4double kinEnergy,
                                                                   G4double step,
                                                                   const G4Material* mat,
                                                                   const G4String& name) {
  return EnergyAfterStep(kinEnergy, step, mat, FindParticle(name));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4EnergyLossForExtrapolatorForCVH::EnergyBeforeStep(G4double kinEnergy,
                                                                    G4double step,
                                                                    const G4Material* mat,
                                                                    const G4String& name) {
  return EnergyBeforeStep(kinEnergy, step, mat, FindParticle(name));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4EnergyLossForExtrapolatorForCVH::AverageScatteringAngle(G4double kinEnergy,
                                                                          G4double step,
                                                                          const G4Material* mat,
                                                                          const G4String& name) {
  return AverageScatteringAngle(kinEnergy, step, mat, FindParticle(name));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4EnergyLossForExtrapolatorForCVH::EnergyDispersion(G4double kinEnergy,
                                                                    G4double step,
                                                                    const G4Material* mat,
                                                                    const G4String& name) {
  return EnergyDispersion(kinEnergy, step, mat, FindParticle(name));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4EnergyLossForExtrapolatorForCVH::ComputeTrueStep(const G4Material* mat,
                                                                   const G4ParticleDefinition* part,
                                                                   G4double kinEnergy,
                                                                   G4double stepLength) {
  G4double theta = AverageScatteringAngle(kinEnergy, stepLength, mat, part);
  return stepLength * std::sqrt(1.0 + 0.625 * theta * theta);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4EnergyLossForExtrapolatorForCVH::ComputeValue(G4double x,
                                                                const G4PhysicsTable* table,
                                                                size_t idxMat) {
  return (nullptr != table) ? ((*table)[idxMat])->Value(x, index) : 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4EnergyLossForExtrapolatorForCVH::SetVerbose(G4int val) { verbose = val; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4EnergyLossForExtrapolatorForCVH::SetMinKinEnergy(G4double val) { emin = val; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4EnergyLossForExtrapolatorForCVH::SetMaxKinEnergy(G4double val) { emax = val; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4EnergyLossForExtrapolatorForCVH::SetMaxEnergyTransfer(G4double val) { maxEnergyTransfer = val; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
