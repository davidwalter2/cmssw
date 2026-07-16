#ifndef TrackPropagation_Geant4e_MaterialGroupModel_h
#define TrackPropagation_Geant4e_MaterialGroupModel_h

// Global material model (Phase A): runtime classifier mapping a Geant4
// step (logical volume, midpoint r, z) to a material group, plus the
// per-group scaling values k_g. Groups are defined by ordered rules
// (name regex + r/z windows + z-side) in a materialGroups tier file
// produced by Analysis/HitAnalyzer/test/makeMaterialGroups.py; group 0
// is the catch-all "other".
//
// Per-step classification (rather than a volume->group map) is required
// because logical volumes are shared between layers (e.g. one
// PixelBarrelActive* LV placed at all three BPIX layers). The regex
// match per logical volume is memoized: each LV resolves once to the
// shortlist of rules whose name pattern matches; per step only the
// r/z/zside windows of that shortlist are evaluated.
//
// Thread-safety: intended to be owned per stream (one instance per
// residual-maker stream instance). The LV shortlist cache is built
// lazily on first use, after the G4 world exists.
//
// See Analysis/HitAnalyzer/doc/global-material-model-plan.md.

#include "SimG4Core/MagneticField/interface/Field.h"

#include <regex>
#include <string>
#include <unordered_map>
#include <vector>

class G4LogicalVolume;

class MaterialGroupModel : public sim::MaterialOffsetProvider {
public:
  // parameter type of the global material groups in the residual-makers'
  // (parmtype, DetId) registry; sentinel DetId(groupIdx), mirroring the
  // parmtype-14 scalar-potential field modes
  static constexpr int ParmTypeMaterialGlobal = 15;

  explicit MaterialGroupModel(const std::string &rulesFile);

  // groupId of a step; 0 = catch-all
  int classify(const G4LogicalVolume *lv, double r_cm, double z_cm) const;

  // k of the step's group plus any FD-closure injection
  double materialOffset(const G4LogicalVolume *lv, double r_cm, double z_cm) const override;

  int nGroups() const { return nGroups_; }
  const std::string &groupName(int g) const { return gnames_[g]; }
  double kValue(int g) const { return kval_[g]; }
  void setKValue(int g, double k) { kval_[g] = k; }
  double priorSigma(int g) const { return prior_[g]; }

  // FD-closure injection: adds eps to group g (g < 0 clears)
  void setInjection(int g, double eps) {
    injGroup_ = g;
    injEps_ = eps;
  }

private:
  struct Rule {
    int groupId;
    std::regex rex;
    double rmin, rmax, zmin, zmax;  // cm; <0 sentinel = unbounded
    int zside;                      // 0 = both
  };

  const std::vector<unsigned short> &shortlist(const G4LogicalVolume *lv) const;

  std::vector<Rule> rules_;
  int nGroups_ = 1;
  std::vector<double> kval_, prior_;
  std::vector<std::string> gnames_;
  int injGroup_ = -1;
  double injEps_ = 0.;

  // lazy per-LV shortlist memo (single-stream use; no locking)
  mutable std::unordered_map<const G4LogicalVolume *, std::vector<unsigned short>> lvCache_;
};

#endif
