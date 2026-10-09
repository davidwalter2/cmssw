#ifndef TrackPropagation_Geant4e_RadiativeDEDXTable_h
#define TrackPropagation_Geant4e_RadiativeDEDXTable_h

// The radiative (bremsstrahlung + pair production) dE/dx of a CVH propagation
// step, from the Geant4 models the simulation and the mean-loss table use:
//
//   mu+-      G4MuBremsstrahlungModel, G4MuPairProductionModel built on mu+
//             (both charges, as G4TablesForExtrapolatorForCVH::ComputeMuonDEDX)
//   e+-       G4SeltzerBergerModel below 1 GeV, G4eBremsstrahlungRelModel from
//             1 GeV (G4eBremsstrahlung's split); no pair production
//   hadrons   G4hBremsstrahlungModel, G4hPairProductionModel built on the
//             particle itself
//
// all unrestricted (cut = kinetic energy).  The models integrate their
// differential cross sections numerically on every call, ~50 us per step for a
// muon, which made this evaluation the largest single cost of the CVH refit.
// The tabulated path reads the same function from tables instead -- per
// element for muons and hadrons, per material for e+- -- of ln(dE/dx) on a
// uniform grid in ln E, interpolated by local cubics, their nodes evaluated by
// the very models above.  Accuracy and construction are documented in
// RadiativeDEDXTable.cc.

class G4Material;
class G4ParticleDefinition;

namespace cvhrad {

  // Unrestricted bremsstrahlung and pair-production dE/dx of `part` at kinetic
  // energy `ekin` in `mat`, in Geant4 units (MeV/mm).  `tabulated` reads the
  // tables (built on first use, shared by all threads), otherwise the models
  // are evaluated directly; outside the tables' domain the tabulated path
  // evaluates directly as well.
  void radiativeDEDX(
      const G4ParticleDefinition *part, const G4Material *mat, double ekin, bool tabulated, double &brem, double &pair);

}  // namespace cvhrad

#endif
