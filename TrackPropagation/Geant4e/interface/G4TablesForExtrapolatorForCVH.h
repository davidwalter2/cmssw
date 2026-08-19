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
//---------------------------------------------------------------------------
//
// ClassName:    G4TablesForExtrapolatorForCVH
//
// Description:  This class keep dedx, range, inverse range tables
//               for extrapolator
//
// Author:       24.10.14 V.Ivanchenko
//
// Modification:
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef TrackPropagation_G4TablesForExtrapolatorForCVH_h
#define TrackPropagation_G4TablesForExtrapolatorForCVH_h 1

#include <map>
#include "globals.hh"
#include "G4PhysicsTable.hh"
#include "G4DataVector.hh"
#include <vector>

class G4ParticleDefinition;
class G4ProductionCuts;
class G4MaterialCutsCouple;
class G4LossTableBuilder;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

enum ExtTableType {
  fDedxElectron = 0,
  fDedxPositron,
  fDedxProton,
  fDedxMuon,
  fRangeElectron,
  fRangePositron,
  fRangeProton,
  fRangeMuon,
  fInvRangeElectron,
  fInvRangePositron,
  fInvRangeProton,
  fInvRangeMuon,
  fMscElectron,
  // CHARGE-AWARE reference (cvhcgf::referenceIsChargeAware, CVH_REF_CHARGEAWARE).
  //
  // The six tables below are the NEGATIVE partners of fDedx/fRange/fInvRange
  // Muon and Proton. They are built only when the switch is set, by the SAME
  // ComputeMuonDEDX / ComputeProtonDEDX with G4MuonMinus and G4AntiProton in
  // place of G4MuonPlus and G4Proton, and are nullptr otherwise.
  //
  // Appended at the END of the enum deliberately: the existing values are
  // stored in nothing, but they are read by a switch with no default and by
  // several call sites, and renumbering them buys nothing.
  fDedxMuonMinus,
  fRangeMuonMinus,
  fInvRangeMuonMinus,
  fDedxAntiProton,
  fRangeAntiProton,
  fInvRangeAntiProton
};

class G4TablesForExtrapolatorForCVH {
public:
  // THE GRID, in one place, because the two instantiation sites used to state
  // it separately and disagreed: the reference trajectory
  // (G4EnergyLossForExtrapolatorForCVH) built 80 bins over 1 MeV - 100 TeV
  // while the fluctuation model (G4UniversalFluctuationForExtrapolator) built
  // 70 bins over 1 MeV - 10 TeV. Both take it from here now, so the mean loss
  // the noise model reads and the mean loss the reference integrates come off
  // the same nodes by construction.
  //
  // The two grids happened to be ALIGNED -- (Emax/Emin)^(1/bins) is 10^0.1 for
  // both, i.e. 10 nodes per decade from the same 1 MeV -- so the short one was
  // the long one truncated, the node VALUES are whatever the same G4 models say
  // at the same energies, and the only thing that can move is the spline whose
  // second derivatives are solved over all nodes.
  //
  // MEASURED, not argued: `calibration_studies/resolution/gridharm_g4driver.cc`
  // fills both vectors with one analytic dE/dx of realistic curvature, calls
  // Geant4's own FillSecondDerivatives on each, and compares Value(E).
  //
  //   the 71 shared nodes coincide to    max |dE|/E = 4.1e-15
  //   long/short - 1 at 0.5 / 1 / 3.136 / 10 / 40 / 100 GeV and 1 TeV:
  //                                      |.| <= 2.2e-16   (i.e. rounding)
  //   long/short - 1 at 5 TeV            -1.7e-08
  //   long/short - 1 at 9 TeV             3.4e-05   <- the OLD grid's end
  //                                                    condition, and there the
  //                                                    long grid is the more
  //                                                    accurate of the two
  //                                                    (2.4e-06 against f)
  //
  // So this is inert at every energy the fit sees and an improvement in the
  // last decade of the grid it replaces. The one real behaviour change is above
  // 10 TeV, where G4PhysicsVector used to clamp and now does not.
  //
  // Deliberately NOT harmonised: `iononly`. The fluctuation model needs the
  // IONIZATION mean loss alone (radiative fluctuation is its own channel),
  // the reference needs the total. That difference is physics, not drift.
  static constexpr G4int kNbins = 80;
  static constexpr G4double kEminMeV = 1.;         // CLHEP::MeV
  static constexpr G4double kEmaxMeV = 1.e8;       // 100 * CLHEP::TeV

  explicit G4TablesForExtrapolatorForCVH(G4int verb, G4int bins, G4double e1, G4double e2, G4bool iononly = false);

  ~G4TablesForExtrapolatorForCVH();

  const G4PhysicsTable* GetPhysicsTable(ExtTableType type) const;

  // Lazily-built radiative (brems + pair) dE/dx for ONE hadron species, on the
  // shared grid but at the particle's OWN kinetic energy -- no proton mass
  // scaling, because radiative loss is not a function of beta*gamma. Returns
  // nullptr when CVH_REF_HADRAD is off or the particle is not a hadron the
  // proton table serves. Built on first use per particle and cached; a model
  // call per lookup would be far too slow (the range-defect integral alone
  // asks 16 times per step).
  const G4PhysicsTable* GetHadronRadiativeTable(const G4ParticleDefinition* part);

  void Initialisation();

  // hide assignment operator
  G4TablesForExtrapolatorForCVH& operator=(const G4TablesForExtrapolatorForCVH& right) = delete;
  G4TablesForExtrapolatorForCVH(const G4TablesForExtrapolatorForCVH&) = delete;

private:
  G4PhysicsTable* PrepareTable(G4PhysicsTable*);

  void ComputeElectronDEDX(const G4ParticleDefinition* part, G4PhysicsTable* table);

  void ComputeMuonDEDX(const G4ParticleDefinition* part, G4PhysicsTable* table);

  void ComputeProtonDEDX(const G4ParticleDefinition* part, G4PhysicsTable* table);

  void ComputeHadronRadiativeDEDX(const G4ParticleDefinition* part, G4PhysicsTable* table);

  void ComputeTrasportXS(const G4ParticleDefinition* part, G4PhysicsTable* table);

  std::vector<const G4MaterialCutsCouple*> couples;
  G4DataVector cuts;

  const G4ParticleDefinition* electron;
  const G4ParticleDefinition* positron;
  const G4ParticleDefinition* muonPlus;
  const G4ParticleDefinition* muonMinus;
  const G4ParticleDefinition* proton;
  const G4ParticleDefinition* antiProton;
  const G4ParticleDefinition* currentParticle = nullptr;

  G4LossTableBuilder* builder = nullptr;
  G4ProductionCuts* pcuts = nullptr;

  G4PhysicsTable* dedxElectron = nullptr;
  G4PhysicsTable* dedxPositron = nullptr;
  G4PhysicsTable* dedxMuon = nullptr;
  G4PhysicsTable* dedxProton = nullptr;
  // one radiative table per hadron species actually seen (CVH_REF_HADRAD)
  std::map<const G4ParticleDefinition*, G4PhysicsTable*> dedxHadRad;
  G4PhysicsTable* rangeElectron = nullptr;
  G4PhysicsTable* rangePositron = nullptr;
  G4PhysicsTable* rangeMuon = nullptr;
  G4PhysicsTable* rangeProton = nullptr;
  G4PhysicsTable* invRangeElectron = nullptr;
  G4PhysicsTable* invRangePositron = nullptr;
  G4PhysicsTable* invRangeMuon = nullptr;
  G4PhysicsTable* invRangeProton = nullptr;
  G4PhysicsTable* mscElectron = nullptr;

  // charge-aware partners; nullptr unless chargeAware
  G4PhysicsTable* dedxMuonMinus = nullptr;
  G4PhysicsTable* rangeMuonMinus = nullptr;
  G4PhysicsTable* invRangeMuonMinus = nullptr;
  G4PhysicsTable* dedxAntiProton = nullptr;
  G4PhysicsTable* rangeAntiProton = nullptr;
  G4PhysicsTable* invRangeAntiProton = nullptr;

  G4double emin;
  G4double emax;
  G4double mass = 0.0;
  G4double charge2 = 0.0;

  G4int verbose;
  G4int nbins;
  G4int nmat = 0;

  // Restore G4 10.4 behaviour: cubic-spline interpolation of the dE/dx /
  // range / inverse-range tables used by the CVH refit. In G4 10.4 the WMass
  // fork queried G4EmParameters::Instance()->Spline() which defaulted to
  // true; that getter was removed in G4 11.x with a class-level default of
  // false. The W mass calibration was tuned against spline-interpolated
  // tables, so keep that behaviour here.
  G4bool splineFlag = true;
  G4bool ionOnly;

  // cvhcgf::referenceIsChargeAware(), latched once in the constructor so that
  // the table build and the dispatch cannot see different values.
  G4bool chargeAware;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
