#ifndef TrackPropagation_Geant4ePropagator_h
#define TrackPropagation_Geant4ePropagator_h

#include <memory>
#include <utility>
#include <vector>

// CMS includes
// - Propagator
#include "TrackingTools/GeomPropagators/interface/Propagator.h"

// - Geant4e
#include "G4ErrorPropagatorData.hh"
#include "G4ErrorPropagatorManager.hh"
#include "G4ErrorSurfaceTarget.hh"

#include <Eigen/Core>

#include "TrackPropagation/Geant4e/interface/G4UniversalFluctuationForExtrapolator.hh"
#include "SimG4Core/MagneticField/interface/Field.h"

class MaterialGroupModel;

/** Propagator based on the Geant4e package. Uses the Propagator class
 *  in the TrackingTools/GeomPropagators package to define the interface.
 *  See that class for more details.
 */

class Geant4ePropagator : public Propagator {
public:
  typedef ROOT::Math::SMatrix<double, 5, 7, ROOT::Math::MatRepStd<double, 5, 7>> AlgebraicMatrix57;
  /** Constructor. Takes as arguments:
   *  * The magnetic field
   *  * The particle name whose properties will be used in the propagation.
   * Without the charge, i.e. "mu", "pi", ...
   *  * The propagation direction. It may be: alongMomentum, oppositeToMomentum
   */
  Geant4ePropagator(const MagneticField *field = nullptr,
                    std::string particleName = "mu",
                    PropagationDirection dir = alongMomentum,
                    double plimit = 1.0,
                    bool forCVH = false,
                    double ioniTruncAlpha = 0.999);

  ~Geant4ePropagator() override;

  /** The methods propagateWithPath() are identical to the corresponding
   *  methods propagate() in what concerns the resulting
   *  TrajectoryStateOnSurface, but they provide in addition the
   *  exact path length along the trajectory.
   *  All of these method calls are internally mapped to
   */

  std::pair<TrajectoryStateOnSurface, double> propagateWithPath(const FreeTrajectoryState &,
                                                                const Plane &) const override;

  std::pair<TrajectoryStateOnSurface, double> propagateWithPath(const FreeTrajectoryState &,
                                                                const Cylinder &) const override;

  std::pair<TrajectoryStateOnSurface, double> propagateWithPath(const TrajectoryStateOnSurface &,
                                                                const Plane &) const override;

  std::pair<TrajectoryStateOnSurface, double> propagateWithPath(const TrajectoryStateOnSurface &,
                                                                const Cylinder &) const override;

  // Deep-copy ctor: allocates an independent fluct so each Propagator
  // instance has its own per-event state. Required for Tier-3 MT: per-stream
  // residual-makers clone the ES-supplied propagator in produce() so no
  // mutable state is shared across streams.
  Geant4ePropagator(const Geant4ePropagator &other);
  Geant4ePropagator &operator=(const Geant4ePropagator &) = delete;

  Geant4ePropagator *clone() const override { return new Geant4ePropagator(*this); }

  const MagneticField *magneticField() const override { return theField; }

  // Optional `particleNameOverride` selects the Geant4 particle hypothesis
  // for this propagation step (e.g. "pi+", "pi-", "kaon+", "kaon-",
  // "proton", "anti_proton"). Empty string falls back to the propagator's
  // default `theParticleName` (set in the constructor; "mu" → "mu+"/"mu-").
  // Used by the V0 CVH ntuplizer so pion / proton / kaon daughters are
  // propagated with the correct dE/dx and energy-loss fluctuation rather
  // than as muons.
  // dB: 3-vector field offset (dBx, dBy, dBz) in Tesla applied to the field
  // along the propagation step. Replaces the older scalar dBz arg as part of
  // the scalar-potential B-field correction. The transport Jacobian exposes
  // dBx, dBy, dBz columns (cols 5,6,7); chain-rule scaling for the
  // scalar-potential modes is applied by the caller.
  std::tuple<bool,
             Eigen::Matrix<double, 7, 1>,
             Eigen::Matrix<double, 5, 5>,
             Eigen::Matrix<double, 5, 9>,
             double,
             Eigen::Matrix<double, 5, 5>,
             Eigen::Matrix<double, 5, 5>,
             double,
             double>
  propagateGenericWithJacobianAltD(const Eigen::Matrix<double, 7, 1> &ftsStart,
                                   const GloballyPositioned<double> &pDest,
                                   const Eigen::Vector3d &dB = Eigen::Vector3d::Zero(),
                                   double dxi = 0.,
                                   double dms = 0.,
                                   double dioni = 0.,
                                   double pforced = -1.,
                                   const std::string &particleNameOverride = std::string(),
                                   // Global material model (Phase A): when matGroups is set,
                                   // the energy loss of every step is scaled by the step
                                   // volume's group value k_g (on top of dxi), and, when
                                   // groupJacOut is also set, the per-group transported
                                   // d(state)/dk_g columns are accumulated into it (the sum
                                   // over groups equals the integrated dxi column of the 5x9
                                   // Jacobian exactly). See doc/global-material-model-plan.md.
                                   const MaterialGroupModel *matGroups = nullptr,
                                   std::vector<std::pair<int, Eigen::Matrix<double, 5, 1>>>
                                       *groupJacOut = nullptr,
                                   // Per-step field modes (leg-structure-free attribution):
                                   // when fieldModes is set, the correction field is applied
                                   // per step (evaluated at the step start, on top of the
                                   // constant dB argument) and, when modeJacOut is also set,
                                   // the per-mode transported d(state)/dc_i columns are
                                   // accumulated with the basis evaluated at each step
                                   // midpoint -- replacing the maker-side per-leg chain rule.
                                   const sim::FieldModeProvider *fieldModes = nullptr,
                                   std::vector<Eigen::Matrix<double, 5, 1>> *modeJacOut =
                                       nullptr) const;

  static void CalculateEffectiveZandA(const G4Material *mate, G4double &effZ, G4double &effA);

  bool GetForCVH() const { return forCVH_; }

  // Per-step Urban-model log for the physics-CF export (doRes): one entry
  // per Geant4 step that produced a nonzero ionization-fluctuation
  // contribution during the LAST propagateGenericWithJacobianAltD call
  // (cleared at entry when logging is enabled). cs = Etot/p^3 [GeV^-2] is
  // the linear map from the step's energy loss (record energies in MeV) to
  // the local q/p noise: delta(q/p) = cs * 1e-3 * dE[MeV].
  struct UrbanIoniStep {
    G4UniversalFluctuationForExtrapolator::UrbanFluctRecord rec;
    double cs = 0.;
  };
  void setIoniStepLogging(bool on) { ioniStepLogging_ = on; }
  const std::vector<UrbanIoniStep> &ioniStepLog() const { return ioniStepLog_; }

  // Phase B: per-step raw material/kinematic data for the offline Moliere
  // (screened-Rutherford compound-Poisson) model of the multiple-scattering
  // tail. Deliberately raw -- chi_c^2 / screening-angle formulas and the CF
  // live offline where they can be iterated without rebuilds.
  //   effZ, effA : effective Z, A of the step material
  //   xg         : areal density rho*d of the step [g/cm^2]
  //   pGeV, beta : momentum [GeV] and velocity at the step
  //   thp2       : projected-angle variance ACTUALLY entering Q for this
  //                step (incl. msfact), i.e. errMSIout(lambda,lambda)
  //   dOverX0    : step thickness in radiation lengths
  // Gated by the same setIoniStepLogging flag and cleared per propagate
  // call together with the Urban log.
  struct MoliereMsStep {
    double effZ = 0., effA = 0., xg = 0.;
    double pGeV = 0., beta = 0.;
    double thp2 = 0., dOverX0 = 0.;
    // material-group id of the step (MaterialGroupModel::classify; -1 when
    // no global material model is active) -- lets the offline fit tie the
    // MS scale to the parmtype-15 material groups
    int stepGroup = -1;
  };
  const std::vector<MoliereMsStep> &msStepLog() const { return msStepLog_; }

  // Per-step cumulative transport Jacobian, for the clean-propagation-test
  // model (Analysis/HitAnalyzer/plugins/G4ePropagationExport.cc).
  //
  // The step loop transports the accumulated noise and THEN adds the step's
  // own contribution, so the noise generated at step s is subsequently
  // transported by steps s+1..N only. Writing Jacc_s for the cumulative
  // transport from the leg start through the end of step s (the point where
  // n_s is injected), the exact transport of that noise to the end of the
  // leg is
  //     A_s = Jacc_N * Jacc_s^{-1}
  // which offline turns the per-step Urban/Moliere records into EXACT
  // per-step weights for any linear functional of the final state -- as
  // opposed to the RMS-matched scalar weight per pooled block that the fit
  // exports (resinfv) have to use. That distinction matters precisely in the
  // tails, which is what the clean test exists to probe.
  //
  // One entry per Geant4 step, unconditionally, so that nMs/nIoni (the sizes
  // of the two physics logs AFTER this step) give an unambiguous mapping
  // from a physics-log index back to its step: physics entry j belongs to the
  // first step whose n exceeds j. Separate flag and separate storage from the
  // Urban/Moliere logs so that existing consumers of those records see no
  // change in layout.
  struct StepTransport {
    double jacc[25] = {0.};  // cumulative 5x5 (row-major), curvilinear
    int nMs = 0;
    int nIoni = 0;
  };
  void setStepTransportLogging(bool on) { stepTransportLogging_ = on; }
  const std::vector<StepTransport> &stepTransportLog() const { return stepTransportLog_; }

private:
  typedef std::pair<TrajectoryStateOnSurface, double> TsosPP;
  typedef std::pair<bool, std::shared_ptr<G4ErrorTarget>> ErrorTargetPair;

  // Magnetic field
  const MagneticField *theField;

  // Name of the particle whose properties will be used in the propagation
  std::string theParticleName;

  // Per-exit-point failure counters for propagateGenericWithJacobianAltD.
  // Indices: 0 = configurePropagation (p < plimit), 1 = Geant4 step ierr != 0,
  // 2 = max path length / max iterations, 3 = in-flight momentum drain
  // below plimit (runaway leg ground down by material before reaching the
  // target plane), 4 = propagation reported success but the final state is
  // off the destination plane (e.g. pinned at the G4 world boundary),
  // 5 = state left the field model's validity region (all tracker modules
  // lie inside it, so exiting is unambiguous proof of a runaway leg --
  // fires earliest of the runaway guards). Dumped from the destructor.
  mutable std::array<unsigned long long, 6> propFailCounts_{{0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL}};
  mutable unsigned long long propTotalCalls_{0ULL};
  // Successful backward legs (anyDirection mode): counted so the frequency
  // of the momentum-flipped frame conversion stays observable per job.
  mutable unsigned long long propBackwardLegs_{0ULL};

  // Geant4 11.1 made G4ErrorPropagatorManager / G4ErrorPropagatorData
  // singletons G4ThreadLocal. Fetch them per-call via the static accessors
  // instead of caching pointers at construction (which would alias the
  // construction-thread's TLS instance from any other thread that uses
  // this propagator).
  double plimit_;

  // Transform a CMS Reco detector surface into a Geant4 Target for the error
  // propagation
  template <class SurfaceType>
  ErrorTargetPair transformToG4SurfaceTarget(const SurfaceType &pDest, bool moveTargetToEndOfSurface) const;

  template <class SurfaceType>
  ErrorTargetPair transformToG4SurfaceTargetD(const SurfaceType &pDest, bool moveTargetToEndOfSurface) const;

  // generates the Geant4 name for a particle from the
  // string stored in theParticleName ( set via constructor )
  // and the particle charge.
  // 'mu' as a basis for muon becomes 'mu+' or 'mu-', depening on the charge
  // This method only supports neutral and +/- 1e charges so far
  //
  // returns the generated string
  std::string generateParticleName(int charge) const;

  // flexible method which performs the actual propagation either for a plane or
  // cylinder surface type
  //
  // returns TSOS after the propagation and the path length
  template <class SurfaceType>
  std::pair<TrajectoryStateOnSurface, double> propagateGeneric(const FreeTrajectoryState &ftsStart,
                                                               const SurfaceType &pDest) const;

  // saves the Geant4 propagation direction (Forward or Backward) in the
  // provided variable reference mode and returns true if the propagation
  // direction could be set
  template <class SurfaceType>
  bool configurePropagation(G4ErrorMode &mode,
                            SurfaceType const &pDest,
                            GlobalPoint const &cmsInitPos,
                            GlobalVector const &cmsInitMom) const;

  // special case to determine the propagation direction if the CMS propagation
  // direction 'anyDirection' was set. This method is called by
  // configurePropagation and provides specific implementations for Plane and
  // Cylinder classes
  template <class SurfaceType>
  bool configureAnyPropagation(G4ErrorMode &mode,
                               SurfaceType const &pDest,
                               GlobalPoint const &cmsInitPos,
                               GlobalVector const &cmsInitMom) const;

  // Ensure Geant4 Error propagation is initialized, if not done so, yet
  // if the forceInit parameter is set to true, the initialization is performed,
  // even if already done before.
  // This can be necessary, when Geant4 needs to read in a new MagneticField
  // object, which changed during lumi section crossing
  void ensureGeant4eIsInitilized(bool forceInit) const;
  void ensureGeant4eIsInitilizedForCVH(bool forceInit) const;

  // returns the name of the SurfaceType. Mostly for debug outputs
  template <class SurfaceType>
  std::string getSurfaceType(SurfaceType const &surface) const;

  void debugReportPlaneSetup(GlobalPoint const &posPlane,
                             HepGeom::Point3D<double> const &surfPos,
                             GlobalVector const &normalPlane,
                             HepGeom::Normal3D<double> const &surfNorm,
                             const Plane &pDest) const;

  template <class SurfaceType>
  void debugReportTrackState(std::string const &currentContext,
                             GlobalPoint const &cmsInitPos,
                             CLHEP::Hep3Vector const &g4InitPos,
                             GlobalVector const &cmsInitMom,
                             CLHEP::Hep3Vector const &g4InitMom,
                             const SurfaceType &pDest) const;

  Eigen::Matrix<double, 5, 5> PropagateErrorMSC(const G4Track *aTrack, double pforced = -1.) const;

  std::pair<double, double> computeLandau(const G4Track *aTrack) const;

  double computeErrorIoni(const G4Track *aTrack, double pforced = -1.) const;

  Eigen::Matrix<double, 5, 9> transportJacobianBxByBzD(
      const Eigen::Matrix<double, 7, 1> &start, double s, double dEdx, double mass, const Eigen::Vector3d &dB) const;

  // mutable: allocation deferred from ctors to the first-call init block
  // in propagateGeneric / propagateGenericWithJacobianAltD (both const
  // methods). Allocating in the ctors aborted on MT worker threads where
  // the G4 navigator's world had not yet been set up by CvhWorker.
  mutable G4UniversalFluctuationForExtrapolator *fluct = nullptr;
  bool forCVH_ = false;

  // ionization-variance truncation passed through to fluct (see
  // G4UniversalFluctuationForExtrapolator::SetIoniTruncationAlpha)
  double ioniTruncAlpha_ = 0.999;

  // Urban step logging (see setIoniStepLogging); log is mutable because
  // computeErrorIoni is const.
  bool ioniStepLogging_ = false;
  mutable std::vector<UrbanIoniStep> ioniStepLog_;
  mutable std::vector<MoliereMsStep> msStepLog_;

  // per-step cumulative transport Jacobian log (see setStepTransportLogging)
  bool stepTransportLogging_ = false;
  mutable std::vector<StepTransport> stepTransportLog_;
};

#endif
