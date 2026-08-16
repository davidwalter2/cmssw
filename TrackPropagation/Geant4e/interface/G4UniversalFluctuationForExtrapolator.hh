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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4UniversalFluctuationForExtrapolator
//
// Author:        V.Ivanchenko make a class with the Laszlo Urban model
//
// Creation date: 03.01.2002
//
// Modifications:
//
//
// Class Description:
//
// Implementation of energy loss fluctuations made by L.Urban in 2021

// -------------------------------------------------------------------
//

#ifndef G4UniversalFluctuationForExtrapolator_h
#define G4UniversalFluctuationForExtrapolator_h 1

#include "G4VEmFluctuationModel.hh"
#include "G4ParticleDefinition.hh"
#include "G4Poisson.hh"
#include "G4Threading.hh"
#include <CLHEP/Random/RandomEngine.h>
#include "TrackPropagation/Geant4e/interface/G4TablesForExtrapolatorForCVH.h"

class G4UniversalFluctuationForExtrapolator : public G4VEmFluctuationModel {
public:
  explicit G4UniversalFluctuationForExtrapolator(const G4String& nam = "UniFluc");

  ~G4UniversalFluctuationForExtrapolator() override;

  // CDF fraction of the delta-electron (1/E^2) spectrum kept when computing
  // the truncated mean/variance of the ionization straggling (PANDA
  // PV/01-07 eq. 68-69). The fitted "resolution" is convention-dependent
  // through this cutoff (truncated sigma grows ~3x from alpha=0.99 to
  // 0.999), so the resolution-closure diagnostic scans it.
  void SetIoniTruncationAlpha(double a) { ioniTruncAlpha_ = a; }

  // Per-step Urban-model parameters recorded by the last SampleFluctuations
  // call. These are the exact ingredients of the analytic compound-Poisson
  // characteristic function of the step's ionization straggling:
  //   regime 0: near-Gaussian thick-target regime, variance gsig2 [MeV^2]
  //   regime 1: Glandz -- excitations as Poisson(a1) at e1 and Poisson(a2)
  //             at e2, plus a3 delta-electron collisions with spectrum
  //             ~1/E^2 on [e0r, tmaxr] (all energies MeV; the small-cut
  //             width correction `scaling` multiplies all energies).
  //             gsig2 then holds the RETURNED alpha-truncated variance,
  //             i.e. exactly what enters the track-fit Q matrix -- the
  //             offline CF fit standardizes with it.
  // Exported per resolution block by the CVH maker (doRes) so the offline
  // fit can use the untruncated non-Gaussian model -- unlike the truncated
  // variance, this has no convention dependence.
  //   regime 2/3: as regime 1, EXCEPT that the delta-ray channel is the exact
  //             knock-on cross section rather than a pure 1/E^2.  regime 2 is
  //             the SPIN-1/2 form (mu, p), regime 3 the SPIN-0 form (pi, K),
  //             mirroring G4BetheBlochModel's own `0.5 == spin` branch:
  //                 spin 1/2:  1 - beta^2 T/tmaxr + T^2/(2 E^2)
  //                 spin 0  :  1 - beta^2 T/tmaxr
  //             The spin is taken from G4ParticleDefinition::GetPDGSpin(), not
  //             from a PDG-id table.
  //             Both regimes carry two EXTRA record fields (beta2, etot) and,
  //             when the switch is on, two extra exported columns -- because
  //             beta^2 and E cannot be recovered from tmaxr without knowing the
  //             particle mass, and getting them wrong for a slow hadron is not
  //             a small error (beta^2 is 0.976 for a 3 GeV kaon against
  //             0.99999 for a 3 GeV muon, and it multiplies the whole
  //             suppression term).  They are PHYSICAL values and carry no
  //             `scaling` factor.
  //             Enabled by the environment switch CVH_IONI_EXACTDELTA; see
  //             the long comment in SampleFluctuations.  In this regime `a3`
  //             is NOT a collision count but the spectrum's normalization
  //             xi = 2 pi re^2 me c^2 n_e z^2 L / beta^2 (MeV), so that
  //                 dN/dT = (xi/T^2) [1 - beta^2 T/tmaxr + T^2/(2 E^2)]
  //             on [e0r, tmaxr], with E implied by tmaxr through the two-body
  //             kinematics.  It scales like an energy, i.e. the consumer
  //             multiplies it by `scaling` exactly as it does e1/e2/e0r/tmaxr.
  //             The stride and the field list are UNCHANGED, so a consumer
  //             that switches on `regime` needs no new columns; one that does
  //             not switch on it must refuse regime 2 rather than read `a3`
  //             as a count.
  struct UrbanFluctRecord {
    int regime = -1;
    double gsig2 = 0.;
    double a1 = 0., e1 = 0.;
    double a2 = 0., e2 = 0.;
    double a3 = 0., e0r = 0., tmaxr = 0.;
    double scaling = 1.;
    double beta2 = 0., etot = 0.;   // regime 2/3 only; physical, unscaled
  };
  const UrbanFluctRecord& lastRecord() const { return record_; }
  bool lastRecordValid() const { return recordValid_; }

  // Is CVH_IONI_EXACTDELTA set?  The exporter branches on this to decide the
  // ioniurbanv stride (11 when off -- byte-identical to the historical
  // record -- 13 when on, with beta2 and etot appended AFTER cs so that every
  // existing column index is unchanged).
  static bool exactDeltaEnabled();

  // CVH_IONI_KOKOULIN (read in exactly one place, cvhcgf::ioniKokoulinEnabled)
  // adds Geant4's own Kokoulin radiative correction to the regime-2/3 delta
  // channel's SECOND MOMENT, i.e. to the variance this class returns and hence
  // to `gsig2` and to the track fit's Q(0,0).  It changes NO other record
  // field: a1, e1, a2, e2, a3, e0r, tmaxr, scaling, beta2 and etot are
  // bit-identical with it on, and no new column is needed offline because f_K
  // depends only on T and E, both already carried.  The MEAN is deliberately
  // untouched -- Geant4's dE/dx table already contains the correction, so
  // `meanLoss` and the excitation rescale that pins the block to it are
  // already right and must not move.
  //
  // Is CVH_IONI_URBAN2021 set?  Harmonizes the EXCITATION channels with stock
  // Geant4 11.2.2's G4UniversalFluctuation (the 2021 Urban model, ONE channel
  // at the mean excitation energy) instead of the pre-2021 two-channel form
  // this class was forked from.  See the block comment in SampleFluctuations.
  static bool urban2021Enabled();

  G4double SampleFluctuations(const G4MaterialCutsCouple*,
                              const G4DynamicParticle*,
                              const G4double,
                              const G4double,
                              const G4double,
                              const G4double) override {
    return 0;
  };

  virtual G4double SampleFluctuations(const G4Material*, const G4DynamicParticle*, G4double, G4double, G4double);

  virtual G4double SampleFluctuations2(
      const G4Material*, const G4DynamicParticle*, G4double, G4double, G4double, G4double);

  G4double Dispersion(
      const G4Material*, const G4DynamicParticle*, const G4double, const G4double, const G4double) override {
    return 0;
  };

  virtual G4double Dispersion(const G4Material*, const G4DynamicParticle*, G4double, G4double);

  // Initialisation for a new particle type
  void InitialiseMe(const G4ParticleDefinition*) override;

  // Initialisation prestep
  void SetParticleAndCharge(const G4ParticleDefinition*, G4double q2) override;

  // hide assignment operator
  G4UniversalFluctuationForExtrapolator& operator=(const G4UniversalFluctuationForExtrapolator& right) = delete;
  G4UniversalFluctuationForExtrapolator(const G4UniversalFluctuationForExtrapolator&) = delete;

protected:
  virtual G4double SampleGlandz(CLHEP::HepRandomEngine* rndm, const G4Material*, const G4double tcut);

  inline void AddExcitation(G4double a, G4double e, G4double& eav, G4double& eloss, G4double& esig2, G4double& esig2tot);

  inline void SampleGauss(G4double eav, G4double esig2, G4double& eloss, G4double& esig2tot);

  inline void AddExcitation2(
      CLHEP::HepRandomEngine* rndm, G4double a, G4double e, G4double& eav, G4double& eloss, G4double& esig2);

  inline void SampleGauss2(CLHEP::HepRandomEngine* rndm, G4double eav, G4double esig2, G4double& eloss);

  // delta-electron tail truncation of the ionization variance (see setter)
  G4double ioniTruncAlpha_ = 0.999;

  // last-call Urban record (see lastRecord())
  UrbanFluctRecord record_;
  G4bool recordValid_ = false;

  // particle properties
  G4double particleMass = 0.0;
  G4double m_Inv_particleMass = DBL_MAX;
  G4double m_massrate = DBL_MAX;
  G4double chargeSquare = 1.0;

  // material properties
  G4double ipotFluct = 0.0;
  G4double ipotLogFluct = 0.0;
  G4double e0 = 0.0;
  G4double electronDensity;
  G4double f1Fluct;
  G4double f2Fluct;
  G4double e1Fluct;
  G4double e2Fluct;
  G4double e1LogFluct;
  G4double e2LogFluct;
  G4double esmall;
  G4double e1, e2;

  // model parameters
  G4double minNumberInteractionsBohr = 10.0;
  G4double minLoss;
  G4double nmaxCont = 8.0;
  G4double rate = 0.56;
  G4double fw = 4.0;
  G4double a0 = 42.0;
  G4double w2 = 0.0;
  G4double meanLoss = 0.0;

  const G4ParticleDefinition* particle = nullptr;
  const G4Material* lastMaterial = nullptr;
  G4double* rndmarray = nullptr;
  G4int sizearray = 30;

  // The dE/dx + range + inv-range tables (G4TablesForExtrapolatorForCVH)
  // are read-only after construction and identical for all instances
  // (always built with the same (0, 70, 1 MeV, 10 TeV, ionOnly=true)
  // parameters). Share them as a process-wide static so that multi-stream /
  // multi-thread CVH refits do not allocate ~tens of MB of duplicate
  // material × particle dE/dx tables per stream. Mirrors the existing
  // static pattern in G4EnergyLossForExtrapolatorForCVH::tables.
  static G4TablesForExtrapolatorForCVH* tables;
#ifdef G4MULTITHREADED
  static G4Mutex extrMutex;
#endif

  const G4PhysicsTable* table = nullptr;
  G4double massratio = 1.;
  G4double charge2ratio = 1.;

  // SPECIES-DEPENDENT Tmax (cvhcgf::referenceIsSpeciesDedx, CVH_REF_SPECIESDEDX).
  //
  // `meanLoss = length * dedx` above reads the SAME scaled proton table as
  // G4EnergyLossForExtrapolatorForCVH::ComputeDEDX and carries the same Tmax
  // defect. It sets the Urban channel weights (a1, a2, a3) and, in the
  // Gaussian regime, the returned loss itself, so the reference's mean and
  // the noise model's mean have to move together -- a fix that left them
  // disagreeing would be a new inconsistency in place of the old one.
  //
  // The effect on the VARIANCE is second order and tiny (the excitation
  // channels carry 5.4e-6 of the block's ionisation variance, so a 1.6e-2
  // change of e_exc moves the total by 8.6e-8 for a pion). It is done because
  // "the two must read the same table", not because it is numerically large.
  //
  // `usesScaledProtonTable` is set in SetParticleAndCharge and is true exactly
  // on the branch that does the mass scaling -- so electrons, positrons and
  // both muons are untouched by construction, as they are in the extrapolator.
  G4bool speciesDedx = false;
  G4bool usesScaledProtonTable = false;
  G4double speciesRefMass = 0.;

  // CHARGE-AWARE mean loss (cvhcgf::referenceIsChargeAware, CVH_REF_CHARGEAWARE).
  //
  // Latched in SetParticleAndCharge from the SAME single reader the
  // extrapolator's table build and dispatch use, so `meanLoss` here and the
  // reference trajectory's mean cannot disagree about whether the correction
  // is on -- which is precisely what they did, by the charge-odd 3.2e-3, from
  // the day CVH_REF_CHARGEAWARE was written until 2026-08-16.
  //
  // It selects fDedxMuonMinus / fDedxAntiProton for a negative track, i.e. the
  // tables G4TablesForExtrapolatorForCVH builds by calling ComputeMuonDEDX /
  // ComputeProtonDEDX with the other particle. No physics term is transcribed
  // here either.
  G4bool chargeAware = false;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4UniversalFluctuationForExtrapolator::AddExcitation(
    G4double ax, G4double ex, G4double& eav, G4double& eloss, G4double& esig2, G4double& esig2tot) {
  if (ax > nmaxCont) {
    eav += ax * ex;
    esig2 += ax * ex * ex;
  } else {
    G4double p = ax;
    if (p > 0) {
      eloss += ((p + 1) - 1.) * ex;
      esig2tot += (ax + 1. / 3.) * ex * ex;
    }
  }
}

inline void G4UniversalFluctuationForExtrapolator::SampleGauss(G4double eav,
                                                               G4double esig2,
                                                               G4double& eloss,
                                                               G4double& esig2tot) {
  G4double x = eav;
  G4double sig = std::sqrt(esig2);
  if (eav < 0.25 * sig) {
    x += (1. - 1.) * eav;
    esig2tot += eav * eav / 3.;
  } else {
    x = eav;
    const double alpha = -eav / sig;
    const double beta = eav / sig;
    const double z = 0.5 * (std::erf(beta / std::sqrt(2.)) - std::erf(alpha / std::sqrt(2.)));
    const double phialpha = 1. / std::sqrt(2. * M_PI) * std::exp(-0.5 * alpha * alpha);
    const double phibeta = 1. / std::sqrt(2. * M_PI) * std::exp(-0.5 * beta * beta);
    esig2tot += sig * sig *
                (1. + (alpha * phialpha - beta * phibeta) / z - (phialpha - phibeta) * (phialpha - phibeta) / z / z);

    // Loop checking, 23-Feb-2016, Vladimir Ivanchenko
  }
  eloss += x;
}

inline void G4UniversalFluctuationForExtrapolator::AddExcitation2(CLHEP::HepRandomEngine* rndm,
                                                                  const G4double ax,
                                                                  const G4double ex,
                                                                  G4double& eav,
                                                                  G4double& eloss,
                                                                  G4double& esig2) {
  if (ax > nmaxCont) {
    eav += ax * ex;
    esig2 += ax * ex * ex;
  } else {
    const G4int p = (G4int)G4Poisson(ax);
    if (p > 0) {
      eloss += ((p + 1) - 2. * rndm->flat()) * ex;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4UniversalFluctuationForExtrapolator::SampleGauss2(CLHEP::HepRandomEngine* rndm,
                                                                const G4double eav,
                                                                const G4double esig2,
                                                                G4double& eloss) {
  G4double x = eav;
  const G4double sig = std::sqrt(esig2);
  if (eav < 0.25 * sig) {
    x += (2. * rndm->flat() - 1.) * eav;
  } else {
    do {
      x = G4RandGauss::shoot(rndm, eav, sig);
    } while (x < 0.0 || x > 2 * eav);
    // Loop checking, 23-Feb-2016, Vladimir Ivanchenko
  }
  eloss += x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
