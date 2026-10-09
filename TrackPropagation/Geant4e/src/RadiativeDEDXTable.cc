#include "TrackPropagation/Geant4e/interface/RadiativeDEDXTable.h"

#include "G4DataVector.hh"
#include "G4Element.hh"
#include "G4EmParameters.hh"
#include "G4Material.hh"
#include "G4MollerBhabhaModel.hh"
#include "G4MuBremsstrahlungModel.hh"
#include "G4MuPairProductionModel.hh"
#include "G4MuonPlus.hh"
#include "G4NistManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProductionCutsTable.hh"
#include "G4SeltzerBergerModel.hh"
#include "G4eBremsstrahlungRelModel.hh"
#include "G4hBremsstrahlungModel.hh"
#include "G4hPairProductionModel.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <CLHEP/Units/SystemOfUnits.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cfloat>
#include <cmath>
#include <fstream>
#include <map>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>

// THE TABLES.  ln(dE/dx) is tabulated on a uniform grid in ln E and read back
// by the cubic through the four nearest nodes.  The nodes are the models' own
// values, so a table can only add interpolation error, and that error is
// controlled where it can arise:
//
//   * WHAT IS TABULATED.  The muon and hadron models sum independent per-atom
//     losses over the material's elements (ComputeDEDXPerVolume: dedx +=
//     n_i loss(Z_i, T)), so they are tabulated per ELEMENT -- the per-atom
//     loss of their own ComputMuBremLoss / ComputMuPairLoss -- and summed with
//     the material's atom densities as the model does.  An element's table
//     serves every material containing it, and each node costs one element,
//     not a whole mixture.  The e+- models are not separable (the dielectric
//     suppression and the LPM effect depend on the material), so they are
//     tabulated per material.
//   * BREAKPOINTS.  A curve is split at the energies where the model's code
//     changes form, each computed from the model's own formula: thresholds
//     below which it returns zero, the e+- switch from Seltzer-Berger to the
//     relativistic model at 1 GeV, the material's LPM threshold (where
//     G4eBremsstrahlungRelModel changes its cross-section formula), the energy
//     nodes of the Seltzer-Berger data (bilinear in ln T, so the slope jumps at
//     each), the changes of the quadrature interval counts, and the energy at
//     which G4MuPairProductionModel's screening term zeta switches on.  Every
//     piece has its own grid (at least four intervals), ending a relative 1e-9
//     inside its breakpoints, and is never interpolated across; the 1e-9
//     slivers are evaluated directly.
//     The muon bremsstrahlung's electron term is clamped at zero, which each
//     of the 48 Gauss nodes of its energy integral crosses at its own energy
//     (a slope kink of up to 2e-3 in ln(dE/dx), 4e-2 for hydrogen): those
//     are breakpoints too.
//   * A CHECK OF EVERY INTERVAL.  The pair production and the LPM form carry
//     more such structure, clamps and switches inside their nested
//     quadratures and the LPM functions' own linear table, and interpolation
//     error at a kink falls only like the grid step, not its fourth power.  So
//     every interval is checked when its table is built: by the fourth
//     differences of the nodes, which bound the error of a kink or jump
//     anywhere in the cubic's stencil, and by the model itself at the
//     midpoint.  An interval that fails is re-tabulated on 16 sub-intervals,
//     each checked the same way, and a sub-interval that still fails (a jump
//     in the model) is evaluated directly.  The tolerance, kCheckTol, is
//     relative to the radiative loss and -- for e+-, whose radiative loss is
//     the larger part of the total above the critical energy -- to the
//     continuous loss the propagator takes as the difference.  It keeps the
//     bound local to every interval of every element and material, whatever
//     the composition or the Geant4 version.
//
// Outside the tabulated domain, at and below a model's threshold and in the
// intervals the check rejects, the models are evaluated directly, i.e.
// exactly as without the table.
//
// Measured against the models over the 467 materials of the CMS 2016 Geant4
// geometry, 1 MeV - 10 TeV: the largest relative error is 3.9e-7 (muon pair
// production just above its threshold), the rms 1e-9; relative to the step's
// total dE/dx, 7e-9 for muons.  Building every element's muon table takes
// 1.6 s, a material's e+- table 0.25 s.
//
// Tables are built on first use and shared by all threads: published through
// atomic pointers, built by whichever thread first needs them on that thread's
// own model instances.  Two threads racing on the same new table both build
// it and one copy is discarded; the copies are identical.

namespace cvhrad {

  namespace {

    // Table domain (kinetic energy): every track the CVH fit propagates.
    constexpr double kEminMeV = 1.;
    constexpr double kEmaxMeV = 1.e6;
    // grid densities, per decade of kinetic energy
    constexpr int kNodesPerDecadeElement = 256;
    constexpr int kNodesPerDecadeMaterial = 128;
    // error accepted in an interval, relative to the radiative loss and to the
    // continuous loss the propagator derives from it (see Curve::build)
    constexpr double kCheckTol = 3e-7;
    // relative distance (in ln E) of a piece's end nodes from its breakpoints,
    // and the shortest piece that is tabulated
    constexpr double kNudge = 1e-9;
    constexpr double kMinPiece = 1e-6;
    // elements are tabulated by integer Z
    constexpr int kMaxZ = 120;

    enum class Family { muon, electron, hadron };

    Family familyOf(const G4ParticleDefinition *part) {
      const int apdg = std::abs(part->GetPDGEncoding());
      return apdg == 11 ? Family::electron : (apdg == 13 ? Family::muon : Family::hadron);
    }

    // ---- the models: one set per thread, never deleted (G4 model statics) ----

    // The per-atom losses the muon/hadron ComputeDEDXPerVolume sums, exposed
    // (they are protected).  The subclasses add no state, so their
    // ComputeDEDXPerVolume is the base model's.
    template <class Base>
    struct BremAtoms : public Base {
      explicit BremAtoms(const G4ParticleDefinition *p) : Base(p) {}
      // unrestricted: ComputeDEDXPerVolume's cut is min(max(T, 0.9 keV), T) = T
      double atomLoss(double Z, double T) { return this->ComputMuBremLoss(Z, T, T); }

      // G4MuBremsstrahlungModel's electron term,
      //     fe = ln( rab2 m / ((1 + delta rmass / (m_e sqrt(e))) (m_e + delta sqrt(e) rab2)) ),
      // is clamped at zero, which it reaches at one value delta* of
      // delta = m^2 v / (2 (E - k)) (before its kinematic limit, where
      // delta = m_e).  The loss integral's 48 Gauss nodes sit at fixed
      // k/T = c, so each crosses delta* where
      //     delta* (1 - c) T^2 + (delta* (2 - c) m - m^2 c / 2) T + delta* m^2 = 0,
      // a slope kink of up to 2e-3 in ln(loss), ~4e-2 for hydrogen.  Valid
      // where the quadrature has its full 8 intervals, T / (T + m) >= 0.15.
      // (G4hBremsstrahlungModel has no electron term.)
      void electronTermBreaks(double Z, std::vector<double> &b) const {
        const int iz = std::min(std::max(G4lrint(Z), 1), 92);
        const double z13 = 1. / this->nist->GetZ13(iz);
        const double rab2 = ((iz == 1) ? this->bh1 : this->btf1) * z13 * z13;
        const double me = CLHEP::electron_mass_c2, m = this->mass;
        const double ka = this->rmass / (me * this->sqrte), kb = this->sqrte * rab2;
        // (1 + ka d)(m_e + kb d) = rab2 m, the positive root
        const double qa = ka * kb, qb = kb + ka * me, qc = me - rab2 * m;
        const double dstar = (-qb + std::sqrt(qb * qb - 4. * qa * qc)) / (2. * qa);
        for (int l = 0; l < 8; ++l) {
          for (int i = 0; i < 6; ++i) {
            const double c = (l + this->xgi[i]) / 8.;
            const double A = dstar * (1. - c), B = dstar * (2. - c) * m - 0.5 * m * m * c, C = dstar * m * m;
            const double disc = B * B - 4. * A * C;
            if (disc > 0.) {
              const double q = -0.5 * (B + std::copysign(std::sqrt(disc), B));
              for (double T : {q / A, C / q}) {
                if (T > 0. && T >= 0.15 / 0.85 * m) {
                  b.push_back(T);
                }
              }
            }
          }
        }
      }
    };
    template <class Base>
    struct PairAtoms : public Base {
      explicit PairAtoms(const G4ParticleDefinition *p) : Base(p) {}
      // MaxSecondaryEnergyForElement also sets the element's screening cache,
      // which ComputMuPairLoss reads
      double atomLoss(double Z, double T) {
        const double tmax = this->MaxSecondaryEnergyForElement(T, Z);
        return this->ComputMuPairLoss(Z, T, T, tmax);
      }
    };

    // muons: built on mu+ with unrestricted cuts, exactly as
    // G4TablesForExtrapolatorForCVH::ComputeMuonDEDX (whose table is built for
    // mu+ and used for both charges).
    struct MuonModels {
      const G4ParticleDefinition *muPlus = nullptr;
      PairAtoms<G4MuPairProductionModel> *pair = nullptr;
      BremAtoms<G4MuBremsstrahlungModel> *brem = nullptr;
    };
    const MuonModels &muonModels() {
      static thread_local MuonModels m;
      if (m.pair == nullptr) {
        m.muPlus = G4MuonPlus::MuonPlus();
        G4DataVector cuts(std::max<size_t>(G4Material::GetNumberOfMaterials(), 1), DBL_MAX);
        m.pair = new PairAtoms<G4MuPairProductionModel>(m.muPlus);
        m.brem = new BremAtoms<G4MuBremsstrahlungModel>(m.muPlus);
        m.pair->Initialise(m.muPlus, cuts);
        m.brem->Initialise(m.muPlus, cuts);
        m.pair->SetUseBaseMaterials(false);
        m.brem->SetUseBaseMaterials(false);
      }
      return m;
    }

    // hadrons: the mass-aware subclasses of the muon models the simulation runs
    // (hBrems / hPairProd), built on the particle itself.
    struct HadronModels {
      BremAtoms<G4hBremsstrahlungModel> *brem = nullptr;
      PairAtoms<G4hPairProductionModel> *pair = nullptr;
    };
    const HadronModels &hadronModels(const G4ParticleDefinition *part) {
      static thread_local std::map<const G4ParticleDefinition *, HadronModels> models;
      auto it = models.find(part);
      if (it == models.end()) {
        G4DataVector cuts(std::max<size_t>(G4Material::GetNumberOfMaterials(), 1), DBL_MAX);
        auto *hb = new BremAtoms<G4hBremsstrahlungModel>(part);
        auto *hp = new PairAtoms<G4hPairProductionModel>(part);
        hb->Initialise(part, cuts);
        hp->Initialise(part, cuts);
        hb->SetUseBaseMaterials(false);
        hp->SetUseBaseMaterials(false);
        it = models.emplace(part, HadronModels{hb, hp}).first;
      }
      return it->second;
    }

    // e+-: G4eBremsstrahlung's models, Seltzer-Berger below 1 GeV and the
    // relativistic model above (G4eBremsstrahlung::InitialiseEnergyLossProcess),
    // the split G4TablesForExtrapolatorForCVH::ComputeElectronDEDX builds the
    // mean from.  One pair per particle: Seltzer-Berger carries the positron
    // correction.
    struct ElectronModels {
      G4VEmModel *low = nullptr;
      G4VEmModel *high = nullptr;
      // the ionisation loss of G4TablesForExtrapolatorForCVH::ComputeElectronDEDX,
      // only to scale the table's tolerance (buildElectron)
      G4VEmModel *ioni = nullptr;
    };
    const ElectronModels &electronModels(const G4ParticleDefinition *part) {
      static thread_local std::map<const G4ParticleDefinition *, ElectronModels> models;
      auto it = models.find(part);
      if (it == models.end()) {
        const size_t n = std::max<size_t>({G4Material::GetNumberOfMaterials(),
                                           G4ProductionCutsTable::GetProductionCutsTable()->GetTableSize(),
                                           size_t(1)});
        G4DataVector cuts(n, DBL_MAX);
        auto *low = new G4SeltzerBergerModel();
        auto *high = new G4eBremsstrahlungRelModel();
        auto *ioni = new G4MollerBhabhaModel();
        low->Initialise(part, cuts);
        high->Initialise(part, cuts);
        ioni->Initialise(part, cuts);
        low->SetUseBaseMaterials(false);
        high->SetUseBaseMaterials(false);
        ioni->SetUseBaseMaterials(false);
        it = models.emplace(part, ElectronModels{low, high, ioni}).first;
      }
      return it->second;
    }

    // ---- direct evaluation ----

    // the material's dE/dx [0] bremsstrahlung, [1] pair production
    double directDEDX(Family fam, int ch, const G4ParticleDefinition *part, const G4Material *mat, double e) {
      switch (fam) {
        case Family::muon: {
          const MuonModels &m = muonModels();
          return ch == 0 ? m.brem->ComputeDEDXPerVolume(mat, m.muPlus, e, e)
                         : m.pair->ComputeDEDXPerVolume(mat, m.muPlus, e, e);
        }
        case Family::electron: {
          if (ch != 0) {
            return 0.;
          }
          const ElectronModels &m = electronModels(part);
          return (e < CLHEP::GeV ? m.low : m.high)->ComputeDEDXPerVolume(mat, part, e, e);
        }
        case Family::hadron: {
          const HadronModels &m = hadronModels(part);
          return ch == 0 ? m.brem->ComputeDEDXPerVolume(mat, part, e, e)
                         : m.pair->ComputeDEDXPerVolume(mat, part, e, e);
        }
      }
      return 0.;
    }

    // the per-atom loss of element Z (muons and hadrons)
    double atomLoss(Family fam, int ch, const G4ParticleDefinition *part, double Z, double e) {
      if (fam == Family::muon) {
        const MuonModels &m = muonModels();
        return ch == 0 ? m.brem->atomLoss(Z, e) : m.pair->atomLoss(Z, e);
      }
      const HadronModels &m = hadronModels(part);
      return ch == 0 ? m.brem->atomLoss(Z, e) : m.pair->atomLoss(Z, e);
    }

    // the model's threshold: ComputeDEDXPerVolume is zero at and below it
    double threshold(Family fam, int ch, const G4ParticleDefinition *part, const G4Material *mat) {
      switch (fam) {
        case Family::muon: {
          const MuonModels &m = muonModels();
          return ch == 0 ? m.brem->MinPrimaryEnergy(mat, m.muPlus, 0.) : m.pair->MinPrimaryEnergy(mat, m.muPlus, 0.);
        }
        case Family::electron:
          return electronModels(part).low->MinPrimaryEnergy(mat, part, 0.);
        case Family::hadron: {
          const HadronModels &m = hadronModels(part);
          return ch == 0 ? m.brem->MinPrimaryEnergy(mat, part, 0.) : m.pair->MinPrimaryEnergy(mat, part, 0.);
        }
      }
      return 0.;
    }

    // ---- breakpoints (kinetic energies), from the models' own formulas ----

    // G4MuBremsstrahlungModel::ComputMuBremLoss (also hadrons): the number of
    // quadrature intervals is int(vcut / 0.05) + 5, capped at 8, with
    // vcut = T / (T + m) for the unrestricted cut; it changes at
    // vcut = 0.05, 0.10, 0.15.
    void bremQuadratureBreaks(double mass, std::vector<double> &b) {
      for (int j = 1; j <= 3; ++j) {
        const double v = 0.05 * j;
        b.push_back(mass * v / (1. - v));
      }
    }

    // G4MuPairProductionModel::ComputeDMicroscopicCrossSection (also hadrons):
    // the screening term zeta is zero until E / (m + g1 Z^(2/3) E) exceeds
    // 35.221047195922 (the root of zeta's numerator), then grows like ln E --
    // a slope kink of up to ~0.1 in ln(loss) at that total energy E.
    void pairZetaBreak(double Z, double mass, std::vector<double> &b) {
      constexpr double kZ1Root = 35.221047195922;
      const double z13 = G4NistManager::Instance()->GetZ13(G4lrint(Z));
      const double g1 = (Z < 1.5) ? 4.4e-5 : 1.95e-5;
      const double den = 1. - kZ1Root * g1 * z13 * z13;
      if (den > 0.) {
        b.push_back(kZ1Root * mass / den - mass);
      }
    }

    // The Seltzer-Berger cross section is bilinear in (k/T, ln T) between the
    // nodes of its data files; the unrestricted loss samples fixed k/T, so its
    // slope jumps exactly at the ln T nodes.  The grid is read from the file
    // the model reads (G4SeltzerBergerModel::ReadData), parsed the same way
    // (G4Physics2DVector::Retrieve).
    const std::vector<double> &seltzerBergerLnT(int iz) {
      static std::mutex mtx;
      static std::map<int, std::vector<double>> grids;
      std::lock_guard<std::mutex> lk(mtx);
      auto it = grids.find(iz);
      if (it == grids.end()) {
        std::vector<double> y;
        std::ostringstream path;
        path << G4EmParameters::Instance()->GetDirLEDATA() << "/brem_SB/br" << iz;
        std::ifstream in(path.str());
        int k = 0, nx = 0, ny = 0;
        if (in >> k >> nx >> ny) {
          double v;
          for (int i = 0; i < nx && in >> v; ++i) {
            // the k/T grid
          }
          for (int j = 0; j < ny && in >> v; ++j) {
            y.push_back(v);
          }
          if (static_cast<int>(y.size()) != ny) {
            y.clear();
          }
        }
        it = grids.emplace(iz, std::move(y)).first;
      }
      return it->second;
    }

    // e+- bremsstrahlung: the model switch at 1 GeV; below it the Seltzer-
    // Berger data nodes and the changes of its quadrature interval count,
    // int(20 T / (T + m_e)) + 3; above it the material's LPM threshold,
    // G4eBremsstrahlungRelModel::SetupForMaterial's
    //     E_tot > sqrt(Migdal n_e) * LPMconstant * X0,
    // where the cross section switches to the LPM form.  False if the
    // Seltzer-Berger grid cannot be read: the curve is then not tabulated.
    bool electronBreaks(const G4Material *mat, const G4ParticleDefinition *part, std::vector<double> &b) {
      const double me = part->GetPDGMass();
      b.push_back(CLHEP::GeV);
      for (int j = 1; j < 20; ++j) {
        b.push_back(me * j / (20. - j));
      }
      for (const G4Element *el : *mat->GetElementVector()) {
        const std::vector<double> &y = seltzerBergerLnT(std::max(std::min(el->GetZasInt(), 100), 1));
        if (y.empty()) {
          return false;
        }
        for (double v : y) {
          b.push_back(std::exp(v) * CLHEP::MeV);
        }
      }
      if (G4EmParameters::Instance()->LPM()) {
        const double migdal = 4. * CLHEP::pi * CLHEP::classic_electr_radius * CLHEP::electron_Compton_length *
                              CLHEP::electron_Compton_length;
        const double lpm = CLHEP::fine_structure_const * CLHEP::electron_mass_c2 * CLHEP::electron_mass_c2 /
                           (4. * CLHEP::pi * CLHEP::hbarc);
        b.push_back(std::sqrt(migdal * mat->GetElectronDensity()) * (lpm * mat->GetRadlen()) - me);
      }
      return true;
    }

    // ---- one tabulated function of kinetic energy ----

    class Curve {
    public:
      // Tabulates f on (lo, hi], split at `breaks`, `nodesPerDecade` per
      // decade.  Every interval is checked (see check) against the relative
      // tolerance tol(E, f(E)) at its midpoint.
      template <class F, class T>
      void build(F f, T tol, double lo, double hi, std::vector<double> breaks, int nodesPerDecade);

      // false: ekin is not tabulated, evaluate directly; lnE = ln(ekin)
      bool value(double ekin, double lnE, double &out) const {
        if (segs_.empty() || !(ekin >= segs_.front().elo && ekin <= segs_.back().ehi)) {
          return false;
        }
        auto it =
            std::upper_bound(segs_.begin(), segs_.end(), ekin, [](double e, const Segment &s) { return e < s.elo; });
        const Segment &s = *(it - 1);
        if (ekin > s.ehi) {
          return false;
        }
        const double t = (lnE - s.x0) * s.invh;
        const int j = std::min(std::max(static_cast<int>(t), 0), s.n - 2);
        switch (flag_[s.off + j]) {
          case kCubic:
            out = std::exp(cubic(&y_[s.off], s.n, j, t));
            return true;
          case kRefined: {
            const std::size_t o = sub_[s.off + j];
            const double ts = (t - j) * kSub;
            const int k = std::min(std::max(static_cast<int>(ts), 0), kSub - 1);
            if (subBad_[o + k]) {
              return false;
            }
            out = std::exp(cubic(&subY_[o], kSub + 1, k, ts));
            return true;
          }
          default:
            return false;
        }
      }

      int nodes() const { return static_cast<int>(y_.size() + subY_.size()); }
      // base intervals evaluated directly, refined ones counted by their share
      double rejected() const {
        return std::count(flag_.begin(), flag_.end(), kDirect) +
               std::count(subBad_.begin(), subBad_.end(), 1) / double(kSub);
      }

    private:
      struct Segment {
        double elo, ehi;  // energies of the first and last node
        double x0, invh;  // ln E of the first node, 1 / grid step
        int n;            // nodes
        std::size_t off;  // first node in y_ (and its interval's entry in flag_, sub_)
      };

      // an interval that fails the check is re-tabulated on kSub sub-intervals
      static constexpr int kSub = 16;
      enum : unsigned char { kCubic = 0, kDirect = 1, kRefined = 2 };

      // cubic through the 4 of the nodes y[0..n) nearest to interval j
      // (one-sided at the ends), at t in units of the node spacing
      static double cubic(const double *y, int n, int j, double t) {
        const int k = std::min(std::max(j - 1, 0), n - 4);
        const double u = t - k, u1 = u - 1., u2 = u - 2., u3 = u - 3.;
        const double *p = y + k;
        return (-p[0] * u1 * u2 * u3 + 3. * p[1] * u * u2 * u3 - 3. * p[2] * u * u1 * u3 + p[3] * u * u1 * u2) / 6.;
      }

      // The check of the intervals of nodes y[0..n) (ln f at x0 + i h): fail[j]
      // unless both estimates of interval j's error pass.  The fourth
      // differences: a single slope kink or jump anywhere in the stencil gives
      // an error of at most 0.375 max|d4| over the fourth differences around
      // the interval (worst case over its position), even where it leaves the
      // midpoint untouched; on a smooth stretch this bounds the error 16x
      // over.  And the model itself at the midpoint (returned in vmid), which
      // sees structure narrower than the grid that the nodes cannot.
      template <class F, class T>
      static void check(F &f,
                        T &tol,
                        const double *y,
                        int n,
                        double x0,
                        double h,
                        std::vector<unsigned char> &fail,
                        std::vector<double> &vmid) {
        const int m = n - 1;
        std::vector<double> d4(n, 0.);
        for (int i = 2; i + 2 <= m; ++i) {
          d4[i] = std::abs(y[i - 2] - 4. * y[i - 1] + 6. * y[i] - 4. * y[i + 1] + y[i + 2]);
        }
        fail.assign(m, 0);
        vmid.assign(m, 0.);
        for (int j = 0; j < m; ++j) {
          double est = 0.;
          for (int i = std::max(j - 1, 2); i <= std::min(j + 2, m - 2); ++i) {
            est = std::max(est, 0.375 * d4[i]);
          }
          const double e = std::exp(x0 + (j + 0.5) * h);
          const double v = f(e);
          const double t = tol(e, v);
          vmid[j] = v;
          fail[j] = !(est <= t && std::abs(std::exp(cubic(y, n, j, j + 0.5)) - v) <= t * v);
        }
      }

      std::vector<Segment> segs_;
      std::vector<double> y_;              // ln f at the nodes
      std::vector<unsigned char> flag_;    // per interval (indexed as its left node)
      std::vector<std::size_t> sub_;       // refined interval: its first node in subY_
      std::vector<double> subY_;           // ln f at the sub-nodes, kSub + 1 per refined interval
      std::vector<unsigned char> subBad_;  // per sub-interval (indexed as subY_): evaluate directly
    };

    template <class F, class T>
    void Curve::build(F f, T tol, double lo, double hi, std::vector<double> breaks, int nodesPerDecade) {
      breaks.push_back(lo);
      breaks.push_back(hi);
      std::sort(breaks.begin(), breaks.end());
      breaks.erase(std::unique(breaks.begin(), breaks.end()), breaks.end());
      breaks.erase(std::remove_if(breaks.begin(), breaks.end(), [&](double b) { return !(b >= lo && b <= hi); }),
                   breaks.end());
      const double hmax = std::log(10.) / nodesPerDecade;
      std::vector<double> y, z, vmid, zmid;
      std::vector<unsigned char> fail, zfail;
      for (std::size_t ib = 0; ib + 1 < breaks.size(); ++ib) {
        const double x0 = std::log(breaks[ib]) + kNudge;
        const double x1 = std::log(breaks[ib + 1]) - kNudge;
        if (!(x1 - x0 > kMinPiece)) {
          continue;  // evaluated directly
        }
        // at least four intervals, so that every interval has a fourth difference
        const int m = std::max(4, static_cast<int>(std::ceil((x1 - x0) / hmax)));
        const double h = (x1 - x0) / m;
        bool ok = true;
        y.resize(m + 1);
        for (int j = 0; j <= m && ok; ++j) {
          const double v = f(std::exp(x0 + j * h));
          ok = (v > 0. && std::isfinite(v));
          y[j] = ok ? std::log(v) : 0.;
        }
        if (!ok) {
          continue;  // a zero or non-finite node: evaluated directly
        }
        const std::size_t off = y_.size();
        segs_.push_back({std::exp(x0), std::exp(x0 + m * h), x0, 1. / h, m + 1, off});
        y_.insert(y_.end(), y.begin(), y.end());
        flag_.insert(flag_.end(), m + 1, kCubic);
        sub_.insert(sub_.end(), m + 1, 0);
        check(f, tol, y.data(), m + 1, x0, h, fail, vmid);
        for (int j = 0; j < m; ++j) {
          if (!fail[j]) {
            continue;
          }
          // re-tabulate the interval on its own kSub sub-intervals, the
          // midpoint already known; a sub-interval that fails is evaluated
          // directly
          flag_[off + j] = kDirect;
          const double xj = x0 + j * h, hs = h / kSub;
          z.resize(kSub + 1);
          z[0] = y[j];
          z[kSub] = y[j + 1];
          bool zok = (vmid[j] > 0. && std::isfinite(vmid[j]));
          z[kSub / 2] = zok ? std::log(vmid[j]) : 0.;
          for (int k = 1; k < kSub && zok; ++k) {
            if (k != kSub / 2) {
              const double v = f(std::exp(xj + k * hs));
              zok = (v > 0. && std::isfinite(v));
              z[k] = zok ? std::log(v) : 0.;
            }
          }
          if (!zok) {
            continue;
          }
          check(f, tol, z.data(), kSub + 1, xj, hs, zfail, zmid);
          flag_[off + j] = kRefined;
          sub_[off + j] = subY_.size();
          subY_.insert(subY_.end(), z.begin(), z.end());
          subBad_.insert(subBad_.end(), zfail.begin(), zfail.end());
          subBad_.push_back(0);
        }
      }
    }

    // ---- the tables of one species ----

    // muons and hadrons: one element's per-atom losses
    struct ElementCurves {
      Curve ch[2];  // [0] bremsstrahlung, [1] pair production
    };

    ElementCurves *buildElement(Family fam,
                                const G4ParticleDefinition *part,
                                const G4Material *mat,
                                double Z,
                                int nodesPerDecade,
                                double checkTol) {
      auto *ec = new ElementCurves;
      const double mass = part->GetPDGMass();
      for (int c = 0; c < 2; ++c) {
        std::vector<double> b;
        if (c == 0) {
          bremQuadratureBreaks(mass, b);
          if (fam == Family::muon) {
            muonModels().brem->electronTermBreaks(Z, b);
          }
        } else {
          pairZetaBreak(Z, mass, b);
        }
        // The thresholds do not depend on the material.  The tolerance is
        // relative to the loss itself: below 1 TeV the radiative loss of a
        // muon or hadron is at most ~7x its ionisation loss (lead, tungsten),
        // which bounds the error passed on to the continuous loss the
        // propagator takes as the difference.
        ec->ch[c].build([&](double e) { return atomLoss(fam, c, part, Z, e); },
                        [&](double, double) { return checkTol; },
                        std::max(kEminMeV, threshold(fam, c, part, mat)),
                        kEmaxMeV,
                        std::move(b),
                        nodesPerDecade);
      }
      return ec;
    }

    // e+-: the material's bremsstrahlung loss.  The propagator takes the
    // continuous loss as the step's total minus this, and above the critical
    // energy that difference is the small one (1/1000 of the radiative loss
    // for 50 GeV in silicon), so the tolerance is relative to the smaller of
    // the two: the ionisation loss of the mean-loss table's model.
    Curve *buildElectron(const G4ParticleDefinition *part, const G4Material *mat, int nodesPerDecade, double checkTol) {
      auto *c = new Curve;
      std::vector<double> b;
      if (electronBreaks(mat, part, b)) {
        G4VEmModel *ioni = electronModels(part).ioni;
        c->build([&](double e) { return directDEDX(Family::electron, 0, part, mat, e); },
                 [&](double e, double v) {
                   return checkTol * std::min(1., ioni->ComputeDEDXPerVolume(mat, part, e, e) / v);
                 },
                 std::max(kEminMeV, threshold(Family::electron, 0, part, mat)),
                 kEmaxMeV,
                 std::move(b),
                 nodesPerDecade);
      }
      return c;
    }

    // Publishes a lazily built object through an atomic slot: the first
    // thread to need it builds it, a racing duplicate is discarded.
    template <class T, class B>
    const T *getOrBuild(std::atomic<const T *> &slot, B build) {
      const T *p = slot.load(std::memory_order_acquire);
      if (p == nullptr) {
        T *built = build();
        const T *expected = nullptr;
        if (slot.compare_exchange_strong(expected, built, std::memory_order_acq_rel)) {
          p = built;
        } else {
          delete built;
          p = expected;
        }
      }
      return p;
    }

    // A material's view of the tables: its elements' curves and the energy
    // range they cover, per process (muons, hadrons), or its own curve (e+-).
    struct MaterialEntry {
      std::vector<const ElementCurves *> elements;  // nullptr: element not tabulated
      double lo[2] = {0., 0.}, hi[2] = {-1., -1.};
      std::unique_ptr<const Curve> electron;
    };

    class SpeciesTables {
    public:
      SpeciesTables(const G4ParticleDefinition *part, Family fam)
          : part_(part),
            fam_(fam),
            // materials created after the first lookup (none in a CVH job,
            // whose geometry is complete before it propagates) get no table
            nmat_(G4Material::GetNumberOfMaterials()),
            materials_(new std::atomic<const MaterialEntry *>[nmat_]) {
        for (std::size_t i = 0; i < nmat_; ++i) {
          materials_[i].store(nullptr, std::memory_order_relaxed);
        }
        for (auto &e : elements_) {
          e.store(nullptr, std::memory_order_relaxed);
        }
      }

      const G4ParticleDefinition *particle() const { return part_; }
      Family family() const { return fam_; }

      // nullptr: no table for this material, evaluate directly
      const MaterialEntry *material(const G4Material *mat) {
        const std::size_t i = mat->GetIndex();
        if (i >= nmat_) {
          return nullptr;
        }
        return getOrBuild(materials_[i], [&] {
          auto *me = new MaterialEntry;
          if (fam_ == Family::electron) {
            me->electron.reset(buildElectron(part_, mat, kNodesPerDecadeMaterial, kCheckTol));
          } else {
            // every element's curve starts and ends at these nodes
            for (int c = 0; c < 2; ++c) {
              me->lo[c] = std::exp(std::log(std::max(kEminMeV, threshold(fam_, c, part_, mat))) + kNudge);
              me->hi[c] = std::exp(std::log(kEmaxMeV) - kNudge);
            }
            for (const G4Element *el : *mat->GetElementVector()) {
              const double Z = el->GetZ();
              const int iz = G4lrint(Z);
              me->elements.push_back(
                  (iz >= 1 && iz <= kMaxZ && Z == iz)
                      ? getOrBuild(elements_[iz],
                                   [&] { return buildElement(fam_, part_, mat, Z, kNodesPerDecadeElement, kCheckTol); })
                      : nullptr);
            }
          }
          return me;
        });
      }

    private:
      const G4ParticleDefinition *part_;
      Family fam_;
      std::size_t nmat_;
      std::unique_ptr<std::atomic<const MaterialEntry *>[]> materials_;
      std::array<std::atomic<const ElementCurves *>, kMaxZ + 1> elements_;
    };

    // The species' tables, keyed by the particle the models are built on
    // (mu+ for both muon charges).  Process-wide and never destroyed, like the
    // models; each thread keeps the pointers it has looked up.
    SpeciesTables &speciesTables(const G4ParticleDefinition *key, Family fam) {
      static thread_local std::vector<SpeciesTables *> seen;
      for (SpeciesTables *t : seen) {
        if (t->particle() == key) {
          return *t;
        }
      }
      static std::mutex mtx;
      static auto *all = new std::vector<SpeciesTables *>;
      SpeciesTables *found = nullptr;
      {
        std::lock_guard<std::mutex> lk(mtx);
        for (SpeciesTables *t : *all) {
          if (t->particle() == key) {
            found = t;
          }
        }
        if (found == nullptr) {
          found = new SpeciesTables(key, fam);
          all->push_back(found);
        }
      }
      seen.push_back(found);
      return *found;
    }

    // One process of a muon or hadron in `mat` from its elements' tables, as
    // ComputeDEDXPerVolume sums them.  Outside the tabulated range (where the
    // threshold lies) the material is evaluated directly; inside it an
    // element whose table rejects the energy is evaluated directly alone.
    double fromElements(Family fam,
                        int ch,
                        const G4ParticleDefinition *part,
                        const G4Material *mat,
                        const MaterialEntry &me,
                        double ekin) {
      if (!(ekin >= me.lo[ch] && ekin <= me.hi[ch])) {
        return directDEDX(fam, ch, part, mat, ekin);
      }
      const G4ElementVector *els = mat->GetElementVector();
      const double *natoms = mat->GetAtomicNumDensityVector();
      const double lnE = std::log(ekin);
      double dedx = 0.;
      for (std::size_t i = 0; i < me.elements.size(); ++i) {
        double loss;
        if (me.elements[i] == nullptr || !me.elements[i]->ch[ch].value(ekin, lnE, loss)) {
          loss = atomLoss(fam, ch, part, (*els)[i]->GetZ(), ekin);
        }
        dedx += loss * natoms[i];
      }
      return std::max(dedx, 0.);
    }

  }  // namespace

  void radiativeDEDX(const G4ParticleDefinition *part,
                     const G4Material *mat,
                     double ekin,
                     bool tabulated,
                     double &brem,
                     double &pair) {
    const Family fam = familyOf(part);
    const G4ParticleDefinition *key = (fam == Family::muon) ? muonModels().muPlus : part;
    const MaterialEntry *me = tabulated ? speciesTables(key, fam).material(mat) : nullptr;
    if (me == nullptr) {
      brem = directDEDX(fam, 0, key, mat, ekin);
      pair = directDEDX(fam, 1, key, mat, ekin);
      return;
    }
    if (fam == Family::electron) {
      if (me->electron == nullptr || !me->electron->value(ekin, std::log(ekin), brem)) {
        brem = directDEDX(fam, 0, key, mat, ekin);
      }
      pair = 0.;
      return;
    }
    brem = fromElements(fam, 0, key, mat, *me, ekin);
    pair = fromElements(fam, 1, key, mat, *me, ekin);
  }

}  // namespace cvhrad
