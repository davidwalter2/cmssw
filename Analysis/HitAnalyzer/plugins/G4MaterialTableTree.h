#ifndef Analysis_HitAnalyzer_G4MaterialTableTree_h
#define Analysis_HitAnalyzer_G4MaterialTableTree_h

// The process-wide Geant4 material table as a TTree, one entry per G4Material:
//     index    G4Material::GetIndex() -- the value MoliereMsStep::materialIndex,
//              and therefore every `msmatv` entry, carries
//     name     G4Material name
//     density  g/cm3
//     elemZ    element Z
//     elemA    element atomic mass, g/mole (the element's own, e.g. H 1.00794)
//     elemW    element mass fraction (G4Material::GetFractionVector)
// These are the numbers the simulation's per-element cross sections use, so a
// reader resolves a step's composition exactly: the mass-averaged effZ/effA and
// the Moliere sum zzp1OverA of the `msmoliv` records are reproducible from it.
//
// Shared by the clean-propagation export (G4ePropagationExport) and the fit
// makers (ResidualGlobalCorrectionMakerBase), which write the same layout; the
// makers also resolve a step's composition through g4MaterialComposition.

#include "G4Material.hh"
#include "G4Element.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "TTree.h"

#include <string>
#include <vector>

namespace cvh {

  // Books the branches on `tree` (which must have none) and fills one entry per
  // material currently in G4Material::GetMaterialTable(). Call it once the
  // Geant4 geometry exists; the branch addresses are reset before returning.
  inline void fillG4MaterialTree(TTree &tree) {
    int index = -1;
    std::string name;
    double density = 0.;
    std::vector<double> elemZ, elemA, elemW;
    tree.Branch("index", &index);
    tree.Branch("name", &name);
    tree.Branch("density", &density);
    tree.Branch("elemZ", &elemZ);
    tree.Branch("elemA", &elemA);
    tree.Branch("elemW", &elemW);
    for (const G4Material *mat : *G4Material::GetMaterialTable()) {
      index = static_cast<int>(mat->GetIndex());
      name = mat->GetName();
      density = mat->GetDensity() / (CLHEP::g / CLHEP::cm3);
      elemZ.clear();
      elemA.clear();
      elemW.clear();
      const G4double *frac = mat->GetFractionVector();
      for (size_t i = 0; i < mat->GetNumberOfElements(); ++i) {
        const G4Element *el = mat->GetElement(i);
        elemZ.push_back(el->GetZ());
        elemA.push_back(el->GetA() / (CLHEP::g / CLHEP::mole));
        elemW.push_back(frac[i]);
      }
      tree.Fill();
    }
    tree.ResetBranchAddresses();
  }

  // The composition of material `index` in the same numbers the tree carries
  // (element Z, element atomic mass in g/mole, mass fraction): what the
  // in-maker nuclear-elastic family builds its per-element mixtures from.
  // False if the index is not in the table.
  inline bool g4MaterialComposition(int index, std::vector<double> &Z, std::vector<double> &A, std::vector<double> &W) {
    const G4MaterialTable *tab = G4Material::GetMaterialTable();
    if (tab == nullptr || index < 0 || static_cast<std::size_t>(index) >= tab->size())
      return false;
    const G4Material *mat = (*tab)[index];
    if (mat == nullptr || static_cast<int>(mat->GetIndex()) != index)
      return false;
    Z.clear();
    A.clear();
    W.clear();
    const G4double *frac = mat->GetFractionVector();
    for (size_t i = 0; i < mat->GetNumberOfElements(); ++i) {
      const G4Element *el = mat->GetElement(i);
      Z.push_back(el->GetZ());
      A.push_back(el->GetA() / (CLHEP::g / CLHEP::mole));
      W.push_back(frac[i]);
    }
    return true;
  }

}  // namespace cvh

#endif
