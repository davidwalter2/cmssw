// MaterialAuditAnalyzer
//
// Phase 0 of the global material model (see
// Analysis/HitAnalyzer/doc/global-material-model-plan.md): geometry-only
// audit of the Geant4 tracker material. Builds the same DDDWorld the CVH
// refit propagates through, then traces straight rays (geantino-style,
// no field, no physics) from the beamline across the tracker volume and
// tallies, per (logical volume, z-side):
//   - path length [cm]
//   - radiation-length fraction  sum(ds/X0)
//   - mass path  sum(ds * rho) [g/cm^2]
//   - approximate ionization energy loss for a 10 GeV muon [MeV]
//     (Bethe-Bloch from the material's electron density and mean
//      excitation energy, no density correction -- adequate for
//      grouping weights, ~10% high in dense materials)
//   - path-weighted mean r and |z| (for classifying the volume)
//   - number of rays crossing
// plus a per-eta profile of the total x/X0 (sanity check against the
// published tracker material-budget plots).
//
// The per-volume table is the input that defines the materialGroups50 /
// materialGroups100 grouping tiers.

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "SimG4Core/Geometry/interface/DDDWorld.h"
#include "SimG4Core/Geometry/interface/SensitiveDetectorCatalog.h"

#include "G4Navigator.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4GeometryManager.hh"
#include "G4ThreeVector.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include <cmath>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

class MaterialAuditAnalyzer : public edm::one::EDAnalyzer<edm::one::WatchRuns> {
public:
  explicit MaterialAuditAnalyzer(const edm::ParameterSet &);

private:
  void beginRun(edm::Run const &, edm::EventSetup const &) override;
  void endRun(edm::Run const &, edm::EventSetup const &) override {}
  void analyze(edm::Event const &, edm::EventSetup const &) override {}
  void endJob() override;

  void traceRays();
  static double muonDedx(const G4Material *mat, double pMuGeV);

  edm::ESGetToken<DDCompactView, IdealGeometryRecord> ddToken_;

  // ray grid + tracking volume bounds (cm)
  double rMax_, zMax_;
  double etaMax_;
  int nEta_, nPhi_;
  std::vector<double> vertexZ_;  // ray origins along the beamline (cm)
  double pMuGeV_;
  std::string outFile_;

  std::unique_ptr<DDDWorld> world_;
  std::unique_ptr<SensitiveDetectorCatalog> catalog_;
  bool done_ = false;

  struct Tally {
    double path = 0.;      // cm
    double xOverX0 = 0.;   // dimensionless
    double massPath = 0.;  // g/cm^2
    double dE = 0.;        // MeV
    double rPath = 0.;     // cm * cm (path-weighted r)
    double zPath = 0.;     // cm * cm (path-weighted |z|)
    unsigned long long nRays = 0;
    std::string material;
  };
  // key: (logical volume name, z-side of step midpoint, coarse r bin
  // [2 cm], coarse |z| bin [10 cm]). The r/z binning keeps shared logical
  // volumes placed at several layers (e.g. PixelBarrelActive* at all
  // three BPIX layers) separable by the downstream per-step classifier.
  struct Key {
    std::string name;
    int side, rbin, zbin;
    bool operator<(const Key &o) const {
      return std::tie(name, side, rbin, zbin) < std::tie(o.name, o.side, o.rbin, o.zbin);
    }
  };
  std::map<Key, Tally> tallies_;
  std::vector<double> etaProfileX0_;  // sum x/X0 per eta bin (averaged over phi, z0)
  std::vector<unsigned long long> etaProfileN_;
};

MaterialAuditAnalyzer::MaterialAuditAnalyzer(const edm::ParameterSet &cfg)
    : ddToken_(esConsumes<DDCompactView, IdealGeometryRecord, edm::Transition::BeginRun>()),
      rMax_(cfg.getParameter<double>("rMax")),
      zMax_(cfg.getParameter<double>("zMax")),
      etaMax_(cfg.getParameter<double>("etaMax")),
      nEta_(cfg.getParameter<int>("nEta")),
      nPhi_(cfg.getParameter<int>("nPhi")),
      vertexZ_(cfg.getParameter<std::vector<double>>("vertexZ")),
      pMuGeV_(cfg.getParameter<double>("muonMomentum")),
      outFile_(cfg.getParameter<std::string>("outFile")) {
  etaProfileX0_.assign(nEta_, 0.);
  etaProfileN_.assign(nEta_, 0ULL);
}

// Bethe-Bloch mean ionization dE/dx [MeV/cm] for a muon of momentum
// pMuGeV in the given material, from electron density and mean
// excitation energy. No density (delta) correction.
double MaterialAuditAnalyzer::muonDedx(const G4Material *mat, double pMuGeV) {
  const double mMu = 105.6583755 * CLHEP::MeV;
  const double me = CLHEP::electron_mass_c2;
  const double p = pMuGeV * CLHEP::GeV;
  const double E = std::sqrt(p * p + mMu * mMu);
  const double beta2 = p * p / (E * E);
  const double gamma = E / mMu;
  const double gamma2 = gamma * gamma;

  const double ne = mat->GetElectronDensity();  // 1/mm3 in G4 units
  const double I = mat->GetIonisation()->GetMeanExcitationEnergy();
  if (ne <= 0. || I <= 0.) {
    return 0.;
  }
  const double tmax =
      2. * me * beta2 * gamma2 / (1. + 2. * gamma * me / mMu + (me / mMu) * (me / mMu));
  const double re = CLHEP::classic_electr_radius;
  const double coef = CLHEP::twopi * re * re * me * ne / beta2;
  const double arg = 2. * me * beta2 * gamma2 * tmax / (I * I);
  if (arg <= 1.) {
    return 0.;
  }
  const double dedx = coef * (std::log(arg) - 2. * beta2);  // MeV/mm (G4 units)
  return dedx / CLHEP::MeV * CLHEP::cm;                     // -> MeV/cm
}

void MaterialAuditAnalyzer::beginRun(edm::Run const &, edm::EventSetup const &iSetup) {
  if (done_) {
    return;
  }
  done_ = true;

  edm::ESTransientHandle<DDCompactView> pDD = iSetup.getTransientHandle(ddToken_);
  catalog_ = std::make_unique<SensitiveDetectorCatalog>();
  world_ = std::make_unique<DDDWorld>(&(*pDD), nullptr, *catalog_, /*verb=*/0,
                                      /*cuts=*/false, /*pcut=*/false);
  // Close the geometry so navigation uses the voxel optimization.
  G4GeometryManager::GetInstance()->CloseGeometry(/*pOptimise=*/true);

  traceRays();
}

void MaterialAuditAnalyzer::traceRays() {
  G4Navigator nav;
  nav.SetWorldVolume(world_->GetWorldVolume());

  const double rMaxMM = rMax_ * CLHEP::cm;
  const double zMaxMM = zMax_ * CLHEP::cm;
  // step epsilon to cross boundaries cleanly
  const double eps = 1.e-4 * CLHEP::mm;

  unsigned long long nrays = 0, nsteps = 0;
  for (double z0cm : vertexZ_) {
    for (int ie = 0; ie < nEta_; ++ie) {
      const double eta = -etaMax_ + (ie + 0.5) * (2. * etaMax_ / nEta_);
      const double theta = 2. * std::atan(std::exp(-eta));
      for (int ip = 0; ip < nPhi_; ++ip) {
        const double phi = -M_PI + (ip + 0.5) * (2. * M_PI / nPhi_);
        const G4ThreeVector dir(std::sin(theta) * std::cos(phi),
                                std::sin(theta) * std::sin(phi), std::cos(theta));
        G4ThreeVector pos(0., 0., z0cm * CLHEP::cm);
        ++nrays;
        double rayX0 = 0.;

        G4VPhysicalVolume *vol = nav.LocateGlobalPointAndSetup(pos, &dir, false, false);
        std::string lastKeyName;
        int lastKeySide = 0;
        while (vol != nullptr) {
          G4double safety = 0.;
          G4double step = nav.ComputeStep(pos, dir, kInfinity, safety);
          if (step == kInfinity || step < 0.) {
            break;
          }
          const G4ThreeVector mid = pos + 0.5 * step * dir;
          const double rmid = mid.perp();
          const double zmid = mid.z();

          const G4LogicalVolume *lv = vol->GetLogicalVolume();
          const G4Material *mat = lv->GetMaterial();
          if (mat != nullptr && step > 0.) {
            const double dscm = step / CLHEP::cm;
            const double x0cm = mat->GetRadlen() / CLHEP::cm;
            const double rho = mat->GetDensity() / (CLHEP::g / CLHEP::cm3);
            const int side = zmid < 0. ? -1 : +1;
            const int rbin = static_cast<int>((rmid / CLHEP::cm) / 2.);
            const int zbin = static_cast<int>((std::abs(zmid) / CLHEP::cm) / 10.);
            auto &t = tallies_[{lv->GetName(), side, rbin, zbin}];
            if (t.material.empty()) {
              t.material = mat->GetName();
            }
            t.path += dscm;
            if (x0cm > 0.) {
              t.xOverX0 += dscm / x0cm;
              rayX0 += dscm / x0cm;
            }
            t.massPath += dscm * rho;
            t.dE += dscm * muonDedx(mat, pMuGeV_);
            t.rPath += dscm * (rmid / CLHEP::cm);
            t.zPath += dscm * (std::abs(zmid) / CLHEP::cm);
            if (lv->GetName() != lastKeyName || side != lastKeySide) {
              ++t.nRays;
              lastKeyName = lv->GetName();
              lastKeySide = side;
            }
            ++nsteps;
          }

          pos += (step + eps) * dir;
          if (pos.perp() > rMaxMM || std::abs(pos.z()) > zMaxMM) {
            break;
          }
          nav.SetGeometricallyLimitedStep();
          vol = nav.LocateGlobalPointAndSetup(pos, &dir, true, false);
        }
        etaProfileX0_[ie] += rayX0;
        ++etaProfileN_[ie];
      }
    }
  }
  edm::LogVerbatim("MaterialAudit") << "MaterialAuditAnalyzer: traced " << nrays << " rays, "
                                    << nsteps << " material steps, " << tallies_.size()
                                    << " (volume, side) tallies";
}

void MaterialAuditAnalyzer::endJob() {
  std::ofstream out(outFile_);
  out << "# MaterialAuditAnalyzer per-volume tally\n"
      << "# rays: eta grid " << nEta_ << " x phi " << nPhi_ << " x vertexZ " << vertexZ_.size()
      << ", |eta| < " << etaMax_ << ", bounds r < " << rMax_ << " cm, |z| < " << zMax_
      << " cm, muon p = " << pMuGeV_ << " GeV\n"
      << "# columns: volume  zside  nRays  path_cm  xOverX0  massPath_gcm2  dE_MeV  meanR_cm  "
         "meanAbsZ_cm  material\n";
  for (auto const &kv : tallies_) {
    auto const &t = kv.second;
    out << kv.first.name << "\t" << kv.first.side << "\t" << t.nRays << "\t" << t.path << "\t"
        << t.xOverX0 << "\t" << t.massPath << "\t" << t.dE << "\t"
        << (t.path > 0. ? t.rPath / t.path : 0.) << "\t" << (t.path > 0. ? t.zPath / t.path : 0.)
        << "\t" << t.material << "\n";
  }
  out << "# eta profile: etaCenter  meanXoverX0\n";
  for (int ie = 0; ie < nEta_; ++ie) {
    const double eta = -etaMax_ + (ie + 0.5) * (2. * etaMax_ / nEta_);
    out << "#ETA\t" << eta << "\t"
        << (etaProfileN_[ie] > 0 ? etaProfileX0_[ie] / etaProfileN_[ie] : 0.) << "\n";
  }
  edm::LogVerbatim("MaterialAudit") << "MaterialAuditAnalyzer: wrote " << outFile_;
}

DEFINE_FWK_MODULE(MaterialAuditAnalyzer);
