#include "ResidualGlobalCorrectionMakerBase.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Analysis/HitAnalyzer/interface/ParticleProperties.h"

// Sparse GBL design-matrix formulation (ported from the single-track
// maker). base.h provides Eigen/Core + Eigen/Eigenvalues + `using
// namespace Eigen`; the sparse solver path needs Eigen/Sparse too.
#include <algorithm>
#include <array>
#include <Eigen/Sparse>
#include <Eigen/Cholesky>

// required for Transient Tracks
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
// required for vtx fitting
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MultiTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MultiTrackPointingKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/CombinedKinematicConstraint.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "RecoLocalTracker/SiStripRecHitConverter/interface/StripCPE.h"
#include "RecoLocalTracker/Records/interface/TkStripCPERecord.h"
#include "RecoLocalTracker/ClusterParameterEstimator/interface/StripClusterParameterEstimator.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"

#include "Math/Vector4Dfwd.h"

#include "TrackPropagation/Geant4e/interface/Geant4ePropagator.h"
#include "TrackPropagation/Geant4e/interface/G4UniversalFluctuationForExtrapolator.hh"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"

#include <iomanip>
#include <map>
#include <iostream>
#include <limits>
#include <sstream>

namespace {

// Generic VertexCompositeCandidate decomposition.
//
// The candidate tree already encodes the fit structure, so no per-channel
// adapter is needed: a LEAF daughter (RecoChargedCandidate) is a single track,
// a COMPOSITE daughter is a two-track joint-fit subsystem. These two helpers
// are the whole of that rule.

// The track behind a leaf daughter, or nullptr if the node is not a leaf
// charged candidate (i.e. it is composite) or carries a null TrackRef.
const reco::Track* leafTrack(const reco::Candidate* c) {
  const auto* rcc = dynamic_cast<const reco::RecoChargedCandidate*>(c);
  if (!rcc || rcc->track().isNull()) return nullptr;
  return &*rcc->track();
}

// The two leaf tracks of a node's first two daughters. Fails (returns false,
// leaving `out` untouched) when the node has fewer than two daughters or when
// either of them is composite / trackless.
bool subsystemPair(const reco::Candidate* c, std::array<const reco::Track*, 2>& out) {
  if (!c || c->numberOfDaughters() < 2) return false;
  const reco::Track* t0 = leafTrack(c->daughter(0));
  const reco::Track* t1 = leafTrack(c->daughter(1));
  if (!t0 || !t1) return false;
  out = {{t0, t1}};
  return true;
}

}  // namespace


class ResidualGlobalCorrectionMakerTwoTrackG4e : public ResidualGlobalCorrectionMakerBase
{
public:
  ResidualGlobalCorrectionMakerTwoTrackG4e(const edm::ParameterSet &);
  ~ResidualGlobalCorrectionMakerTwoTrackG4e() {
    // fix-cvh-displaced-starting-state: report per-job counters for the
    // midPropagated mode. Total = number of (icons-phase × candidate) tries
    // when useStartingState='midPropagated'; Fallback = subset where the
    // analytical extrapolation failed on either daughter and the producer
    // fell back to perigee (Kalman-daughter) per-event.
    if (midPropagatedTotalCount_ > 0ULL) {
      const double fallbackPct = 100.0 * static_cast<double>(midPropagatedFallbackCount_)
                                       / static_cast<double>(midPropagatedTotalCount_);
      edm::LogInfo("ResidualGlobalCorrectionMakerTwoTrackG4e")
          << "midPropagated summary: total=" << midPropagatedTotalCount_
          << "  fallback=" << midPropagatedFallbackCount_
          << " (" << fallbackPct << "%)";
    }
    // Per-stream fit-outcome accounting (see counters below): attempted =
    // track pairs entering the icons/iteration loops; succeeded = pairs
    // whose fit survived all phases (tree filled).
    if (fitAttempted_ > 0ULL) {
      const unsigned long long failTotal = fitAttempted_ - fitSucceeded_;
      std::cout << "ResidualGlobalCorrectionMakerTwoTrackG4e fit summary"
                << "  attempted=" << fitAttempted_
                << "  succeeded=" << fitSucceeded_
                << "  failed=" << failTotal
                << " (" << (100. * failTotal / fitAttempted_) << "%)"
                << "  fail[kinfit]=" << fitFailKinFit_
                << "  fail[prop]=" << fitFailProp_
                << "  fail[hitupdate]=" << fitFailHitUpdate_
                << "  fail[chargeflip]=" << fitFailChargeFlip_
                << "  fail[nan]=" << fitFailNaN_
                << "  fail[ndof]=" << fitFailNdof_
                << "  skipped[samesign]=" << fitSkippedSameSign_
                << "  skipped[ndof<" << minNdof_ << "]=" << fitSkippedNdof_
                << "  skipped[hits<" << minPairHits_ << "]=" << fitSkippedHits_
                << "  skipped[leghits<" << minLegHits_ << "]=" << fitSkippedLegHits_
                << "  clamped[step]=" << fitStepClamped_
                << "  clampevents[step]=" << stepClampEvents_
                << "  backtracked[step]=" << fitStepBacktracked_
                << "  btchi2[step]=" << fitStepBtChi2_
                << "  btchi2events[step]=" << stepBtChi2Events_
                << "  inflated[seed]=" << fitSeedInflated_
                << std::endl;
      if (pixHitsSeen_ > 0ULL) {
        std::cout << "ResidualGlobalCorrectionMakerTwoTrackG4e pixel hit-quality summary"
                  << "  seen=" << pixHitsSeen_
                  << "  onEdge=" << pixHitsEdge_
                  << " (" << (100. * pixHitsEdge_ / pixHitsSeen_) << "%)"
                  << "  sizeX1=" << pixHitsSizeX1_
                  << " (" << (100. * pixHitsSizeX1_ / pixHitsSeen_) << "%)"
                  << "  demoted=" << pixHitsDemoted_
                  << " (" << (100. * pixHitsDemoted_ / pixHitsSeen_) << "%)"
                  << "  keepPixelEdgeHits=" << keepPixelEdgeHits_
                  << "  pixelMinSizeX=" << pixelMinSizeX_
                  << "  pixelMinSizeY=" << pixelMinSizeY_
                  << std::endl;
        // Pathology-class combination table (only non-empty bins). Class 0
        // is the healthy bulk; every other combination is a candidate for a
        // dedicated local-x/local-y correction parameter.
        for (unsigned int isub = 0; isub < 2; ++isub) {
          for (unsigned int icls = 0; icls < 16; ++icls) {
            const unsigned long long n = pixHitsClass_[isub][icls];
            if (n == 0ULL) continue;
            std::string label;
            if (icls & 1) label += "|edgeX";
            if (icls & 2) label += "|edgeY";
            if (icls & 4) label += "|sizeX1";
            if (icls & 8) label += "|sizeY1";
            if (label.empty()) label = "|clean";
            std::cout << "pixHitClass " << (isub == 0 ? "BPix" : "FPix")
                      << " " << label.substr(1)
                      << "  n=" << n
                      << " (" << (100. * n / pixHitsSeen_) << "%)"
                      << std::endl;
          }
        }
      }
    }
  }

// static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:

  virtual void beginStream(edm::StreamID) override;
  virtual void produce(edm::Event &, const edm::EventSetup &) override;
  
  bool doVtxConstraint_;
  bool doMassConstraint_;
  double massConstraint_;
  double massConstraintWidth_;
  double daughterMass1_;
  double daughterMass2_;
  double daughterMass1Err_;
  double daughterMass2Err_;

  // Geant4-particle base names per daughter slot (used to build the per-call
  // particle-name override passed to the Geant4e propagator). Empty string
  // = fall back to the propagator's default (muon). Recognised bases:
  // "mu", "pi", "kaon" --> "<base>+" / "<base>-"
  // "proton" --> "proton" / "anti_proton"
  // "e" --> "e+" / "e-"
  std::string daughterParticleName1_;
  std::string daughterParticleName2_;

  // Optional 2D-transverse pointing-angle constraint on the V0 (KS, Lambda):
  // requires the V0 momentum direction in xy to coincide with the flight
  // vector from the beamspot to the secondary vertex. Default off, so the
  // legacy behaviour for J/psi/Upsilon/D* runners (which never set these)
  // is preserved exactly.
  bool   doPointingConstraint_;
  double pointingSigma_;        // angular pointing resolution, radians

  // §1 of openspec/improve-cvh-refit-convergence: configurable knobs for
  // the CVH joint-refit convergence study. Defaults reproduce the
  // published Run2016H baseline bit-identically (nIters=10,
  // edmConvergence=1e-5, useStartingState="perigee",
  // debugPerIterDump=false). Read with existsAs-guards so legacy cfis
  // that don't set them continue to work.
  unsigned int nIters_;             // Gauss-Newton iteration cap per icons phase
  double       edmConvergence_;     // edmval < this -> break (was hard-coded 1e-5)
  std::string  useStartingState_;   // "perigee" (legacy) | "midPropagated" (lands iter-0 at perigee-midpoint, fix-cvh-displaced-starting-state)
  // Per-job counter for midPropagated mode: number of events where the
  // AnalyticalImpactPointExtrapolator returned an invalid TSOS for either
  // daughter, forcing a per-event fallback to the perigee (Kalman-daughter)
  // path. Logged from EndJob.
  mutable unsigned long long midPropagatedFallbackCount_ = 0ULL;
  mutable unsigned long long midPropagatedTotalCount_    = 0ULL;

  // Per-stream fit-outcome accounting, printed from the destructor.
  mutable unsigned long long fitAttempted_ = 0ULL;
  mutable unsigned long long fitSucceeded_ = 0ULL;
  mutable unsigned long long fitFailKinFit_ = 0ULL;      // seed kinematic vertex fit empty/inconsistent
  mutable unsigned long long fitFailProp_ = 0ULL;        // Geant4e propagation failed
  mutable unsigned long long fitFailHitUpdate_ = 0ULL;   // CPE re-evaluation (cloner) invalid
  mutable unsigned long long fitFailChargeFlip_ = 0ULL;  // q/p sign flip in parameter update
  mutable unsigned long long fitStepClamped_ = 0ULL;     // fits with >=1 momentum-floor-clamped GN step
  mutable unsigned long long stepClampEvents_ = 0ULL;    // individual clamped GN steps
  mutable unsigned long long fitStepBacktracked_ = 0ULL; // step halvings after a failed-leg iteration (retries)
  mutable unsigned long long fitStepBtChi2_ = 0ULL;      // fits with >=1 chi2 (Armijo) backtrack
  mutable unsigned long long stepBtChi2Events_ = 0ULL;   // individual chi2 step halvings
  mutable unsigned long long stepPrints_ = 0ULL;         // printouts emitted (rate limit)
  mutable unsigned long long fitSeedInflated_ = 0ULL;    // iteration-0 seed-momentum inflations (retries)
  // Momentum floor for the Gauss-Newton step clamp (GeV). Its ONLY job is to
  // keep the state out of the propagator's refusal region
  // (Geant4ePropagator.PropagationPtotLimit), so it must sit just above that
  // limit and below the soft-daughter spectrum -- the drivers derive it from
  // the limit, and must re-derive it whenever the limit moves. The 2.0 GeV
  // built-in default sits far above the current 0.2 GeV limit, and there it
  // PINS every daughter below 2 GeV at 2 GeV (momentum-high, chi2/ndof>>1,
  // and, when p_ref is already under the floor, a scale-to-zero frozen step).
  // On the flat-pT J/psi gun that is 12 % of candidates, carrying the whole
  // +0.21e-3 mass-scale offset.
  double clampMomentumFloor_ = 2.0;
  // Relative Gauss-Newton step damping. Per iteration a daughter's
  // momentum may change by at most this factor (default 2: p may at most halve
  // or double). Implemented as the effective floor max(clampMomentumFloor_,
  // p_ref/f) plus the symmetric upward cap p_ref*f. The lower bound is ALWAYS
  // strictly below p_ref, which removes both pathologies of the fixed floor:
  // no daughter is pinned at a fixed momentum, and the scale can no longer be
  // exactly zero (which froze the whole coupled step at its seed). <= 1
  // disables the relative cap and leaves the absolute floor as the only clamp.
  double maxMomentumStepFactor_ = 2.0;
  // chi2-based (Armijo) retroactive backtracking. The chi2 assembled in
  // iteration k is the REALIZED chi2 of the step taken at k-1; if it fails the
  // sufficient-decrease test, the k-1 linearization is restored (the same
  // snapshot the propagation-failure retry uses), the step is halved and the
  // iteration redone. No extra propagation on the accept path.
  bool stepBacktracking_ = true;
  unsigned int maxChi2Backtrack_ = 4;   // halvings per accepted step
  // First iteration at which the Armijo test may fire. Default 2, and NOT 1:
  // at iiter == 0 the GBL propagation/kink residuals are identically zero by
  // construction (the layer states ARE the propagated states, dx0 = 0), so the
  // chi2 at iteration 0 is a DIFFERENT objective from the one at every later
  // iteration and comparing across that boundary would backtrack every
  // candidate. From iteration 1 on the objective is the same function of the
  // state, so the first meaningful comparison is chi2(2) against chi2(1).
  unsigned int stepBacktrackFromIter_ = 2;
  double armijoC_ = 1.e-4;              // sufficient-decrease coefficient
  // Relative chi2 slack in the Armijo test. NOT a textbook line-search
  // tolerance: the CVH/GBL iteration does NOT
  // monotonically decrease r^T Vinv r -- the realized chi2 drifts UP by
  // ~0.3-0.5 per iteration even at 1/16 of the step (the model's predicted
  // decrease is never realized because every iteration re-propagates and
  // re-linearizes). A tolerance of 1e-3 therefore turns the test into a
  // permanent step-halver: 57 % of gun candidates backtracked, 15.6 halvings
  // each, 6x the propagation cost, for no change in the result. At 1.0
  // (the chi2 may not more than DOUBLE in one iteration) the test becomes a
  // pure DIVERGENCE TRAP: +0.5 % propagation on the gun ditrack smoke,
  // +0.2 % single track, fit output at the noise level.
  double armijoSlack_ = 1.0;
  unsigned int stepPrintLimit_ = 200;   // per-job cap on step-control printouts
  // Per-candidate leg-failure retry budgets (see the recovery block).
  unsigned int maxBacktracks_ = 4;
  unsigned int maxSeedInflations_ = 2;

  // Which two-track subsystem of each srcCandidates entry this instance fits.
  //
  // A VertexCompositeCandidate decomposes generically: a COMPOSITE daughter is
  // a joint-fit two-track subsystem (J/psi, K*0, phi, Ks, Lambda, pipi, D0), a
  // LEAF RecoChargedCandidate daughter is a single track (handled by the
  // single-track maker). This index selects which daughter to descend into, so
  // channels with two composite daughters (e.g. B0 -> J/psi K*0) are served by
  // two instances of this module rather than a per-channel splitter.
  //
  //   -1 (default) : legacy/auto -- fit the candidate directly when its first
  //                  two daughters are leaves, otherwise descend into
  //                  daughter(0). Keeps every existing driver bit-identical.
  //   >= 0         : descend into that daughter index and fit its two leaves.
  int subsystemDaughter_ = -1;

  // Optional per-candidate EDM ValueMap output for the NanoAOD path, keyed to
  // the srcCandidates collection. Off by default so the ALCARECO TTree drivers
  // are unaffected. Kinematics are emitted whenever produceValueMaps_ is set;
  // the (large) global-fit payload (globalIdxs / jacRef / jacMass / factored
  // Hessian) is emitted additionally only when fillGradsFactored_ is set.
  bool produceValueMaps_ = false;
  edm::EDPutTokenT<edm::ValueMap<float>> vmCorMass_, vmCorMassErr_, vmCorPt_, vmCorEta_, vmCorPhi_;
  edm::EDPutTokenT<edm::ValueMap<float>> vmMuPlusPt_, vmMuPlusEta_, vmMuPlusPhi_;
  edm::EDPutTokenT<edm::ValueMap<float>> vmMuMinusPt_, vmMuMinusEta_, vmMuMinusPhi_;
  edm::EDPutTokenT<edm::ValueMap<float>> vmEdmval_;
  edm::EDPutTokenT<edm::ValueMap<std::vector<int>>> vmGlobalIdxs_;
  edm::EDPutTokenT<edm::ValueMap<std::vector<float>>> vmJacRefMuPlus_, vmJacRefMuMinus_, vmJacMass_, vmHessFactor_;

  mutable unsigned long long fitFailNaN_ = 0ULL;         // NaN/inf parameter update
  mutable unsigned long long fitFailNdof_ = 0ULL;        // fit with no degrees of freedom (ndof <= 0)
  mutable unsigned long long fitSkippedSameSign_ = 0ULL; // same-sign pairs skipped pre-fit (not failures)
  mutable unsigned long long fitSkippedNdof_ = 0ULL;     // pre-fit ndof < minNdof_ (not failures)
  mutable unsigned long long fitSkippedHits_ = 0ULL;     // pre-fit valid-hit count < minPairHits_
  mutable unsigned long long fitSkippedLegHits_ = 0ULL;  // pre-fit weaker leg below minLegHits_

  // ---- THE MINIMUM-SIZE REQUIREMENT ON A PAIR --------------------------
  //
  // The two-track fit spends TEN state parameters on the common vertex, so its
  // degrees of freedom are
  //
  //     ndof = nvalid + nvalidpixel - 10 + (3 bs) + (1 pointing)
  //            + (1 vertex constraint) + (1 mass constraint)
  //
  // -- one measurement coordinate per strip hit, two per pixel hit. With the
  // vertex constraint on and nothing else that is `n_meas - 9`, so a pair
  // needs MORE THAN NINE measurement coordinates to have any degrees of
  // freedom at all (more than ten with the constraint off). At ndof == 0 the
  // system is exactly determined: chi2 is identically zero, chi2/ndof is 0/0,
  // and the factored-Hessian export indexes one past the end of the
  // eigenvalue vector and aborts the PROCESS -- which is what killed 28 % of
  // the first dymc_8p5M_260905 tasks before the signed-ndof gate below.
  //
  // Those pairs carry no information for the global fit and their vertex
  // residual is meaningless (at ndof == 1 the DCA is determined by the data
  // with one constraint left over, so its pull is not a resolution
  // measurement either). They are cut here, BEFORE the fit, on both readings
  // of "hits":
  //
  //   minNdof      minimum of the expression above evaluated for the
  //                unconstrained-mass pass. Default 1, i.e. n_meas >= 10 with
  //                the vertex constraint on and >= 11 with it off, exactly
  //                the requirement above. 0 restores the legacy behaviour
  //                (only the ndof <= 0 abort).
  //   minPairHits  minimum number of VALID HITS summed over the two legs
  //                (pixel hits counted once, not twice). Default -1 = auto =
  //                10 with the vertex constraint on, 11 with it off. 0
  //                disables.
  //   minLegHits   minimum valid hits on the WEAKER leg. Default 8, because
  //                a thin leg is background: by gen truth on DY, 82-93 % of
  //                the candidates a weaker leg below 8 removes are `dup` or
  //                `unmatched` (0.885 +- 0.026), for a signal efficiency of
  //                0.9983 +- 0.0004 -- and every candidate with a non-finite
  //                mass-resolution export has one. Measured on the 10 654
  //                candidates of `dy_vtxon`: every one of the five with a
  //                non-finite `Jpsi_sigmamass` has a leg of one or two valid
  //                hits -- (2,11), (1,15), (21,2), (14,1), (14,1) -- while
  //                their PAIR totals, 13 to 23 hits, sail through any
  //                pair-level cut. A pair sum cannot see a one-hit leg; this
  //                can. 0 disables. (The COMPANION cut on the vertex
  //                residual, |z_v| < 5, is deliberately NOT made here: it is
  //                a downstream selection, `resolution/selection.py`, so the
  //                tail stays in the trees and the terms that cut on it can
  //                normalise over the window they cut to.)
  //
  // Both are pre-fit, so a skipped pair costs nothing and every candidate that
  // survives is bit-identical to what the previous build wrote.
  int minNdof_ = 1;
  int minPairHits_ = -1;
  int minLegHits_ = 8;

  // Global material model: per-leg per-group dxi columns from the
  // propagator (reused buffer; see the global-material-model design note (kept outside the repository)).
  mutable std::vector<std::pair<int, Eigen::Matrix<double, 5, 1>>> groupJacs_;
  // Per-group PROCESS NOISE of the last propagation (the width counterpart of
  // groupJacs_'s mean). Filled only under doRes + the global material model;
  // see the parmtype-15 dV registration.
  mutable std::vector<std::pair<int, Eigen::Matrix<double, 5, 5>>> groupQs_;
  // per-step field-mode columns from the propagator (reused buffer)
  mutable std::vector<Eigen::Matrix<double, 5, 1>> modeJacs_;

  // Pixel hit-quality accounting (valid pixel hits entering the quality cut).
  mutable unsigned long long pixHitsSeen_ = 0ULL;
  mutable unsigned long long pixHitsEdge_ = 0ULL;       // cluster on the sensor boundary (isOnEdge)
  mutable unsigned long long pixHitsSizeX1_ = 0ULL;     // cluster sizeX == 1
  mutable unsigned long long pixHitsDemoted_ = 0ULL;    // demoted to inactive by the quality cut
  // Pathology-class table: bit0 = edge in local x (cluster touches a
  // sensor-boundary row), bit1 = edge in local y (boundary column),
  // bit2 = sizeX == 1 (single row -> x from one pixel), bit3 = sizeY == 1
  // (single column -> y from one pixel). Outer index 0 = BPix, 1 = FPix.
  mutable std::array<std::array<unsigned long long, 16>, 2> pixHitsClass_ {{{{0ULL}}, {{0ULL}}}};

  // Per-hit pixel diagnostics (gated by fillHitDiagnostics): last-iteration
  // local residuals dy0 + side-resolved pathology class of every valid
  // pixel hit on the two tracks. Class bits: 0=edge at -x boundary,
  // 1=edge at +x, 2=edge at -y, 3=edge at +y, 4=sizeX==1, 5=sizeY==1.
  bool fillHitDiagnostics_ = false;
  // Deweight pathological pixel hits (any class bit set) by scaling their
  // Vinv with 1e-6: the hit keeps its surface and trajectory state but
  // exerts no pull, so its dy0 is an (almost) unbiased residual w.r.t.
  // the surrounding fit -- the measurement mode for the per-side bias
  // attribution study.
  bool deweightPathoHits_ = false;
  std::vector<unsigned int> hitdiag_detid;
  std::vector<int> hitdiag_trk;     // 0/1 = position in the track pair
  std::vector<int> hitdiag_charge;  // charge of the track the hit is on
  std::vector<int> hitdiag_class;
  std::vector<float> hitdiag_dx, hitdiag_dy;    // dy0: hit - predicted, local x/y
  std::vector<float> hitdiag_exx, hitdiag_eyy;  // CPE local position variance
  std::vector<float> hitdiag_lx, hitdiag_ly;    // (deformation-corrected) hit local pos
  bool         debugPerIterDump_;   // emit per-iter vector branches when true

  // Per-iteration debug vectors (filled only when debugPerIterDump_=true).
  // Reset at the top of each event; push_back inside the iter loop.
  // Spans both icons=0 and icons=1 phases (concatenated in iteration order).
  std::vector<double> chisqval_iter;
  std::vector<double> edmval_iter;
  std::vector<double> edmvalref_iter;   // reference-block EDM = the convergence criterion
  std::vector<double> deltachisqval_iter;
  std::vector<double> mu_qoverp_iter;       // 2 entries per iter (muplus, muminus)
  std::vector<double> Jpsi_mass_iter;       // 1 entry per iter

  bool doL1Trigger_;
  edm::EDGetTokenT<L1GlobalTriggerReadoutRecord> inputL1ReadoutRecord_;
  edm::InputTag inputL1ReadoutRecordTag_;
  std::vector<std::string> l1Triggers_;
  std::vector<unsigned int> l1TriggerBitNumbers_;
  std::vector<int> l1TriggerDecisions_;
  
  float Jpsi_d;
  float Jpsi_x;
  float Jpsi_y;
  float Jpsi_z;
  float Jpsi_pt;
  float Jpsi_eta;
  float Jpsi_phi;
  float Jpsi_mass;
  
  float Jpsi_sigmamass;
  
  float Muplus_pt;
  float Muplus_eta;
  float Muplus_phi;
  
  float Muminus_pt;
  float Muminus_eta;
  float Muminus_phi;
  
  float Jpsikin_x;
  float Jpsikin_y;
  float Jpsikin_z;
  float Jpsikin_pt;
  float Jpsikin_eta;
  float Jpsikin_phi;
  float Jpsikin_mass;
  
  float Mupluskin_pt;
  float Mupluskin_eta;
  float Mupluskin_phi;
  
  float Muminuskin_pt;
  float Muminuskin_eta;
  float Muminuskin_phi;
  
  float Jpsitrk_pt;
  float Jpsitrk_eta;
  float Jpsitrk_phi;
  float Jpsitrk_mass;
  
  float Muplustrk_pt;
  float Muplustrk_eta;
  float Muplustrk_phi;
  
  float Muminustrk_pt;
  float Muminustrk_eta;
  float Muminustrk_phi;
  
  float Jpsicons_d = -99.f;
  float Jpsicons_x = -99.f;
  float Jpsicons_y = -99.f;
  float Jpsicons_z = -99.f;
  float Jpsicons_pt = -99.f;
  float Jpsicons_eta = -99.f;
  float Jpsicons_phi = -99.f;
  float Jpsicons_mass = -99.f;
  
  float Mupluscons_pt = -99.f;
  float Mupluscons_eta = -99.f;
  float Mupluscons_phi = -99.f;
  
  float Muminuscons_pt = -99.f;
  float Muminuscons_eta = -99.f;
  float Muminuscons_phi = -99.f;
  
  float Jpsikincons_x = -99.f;
  float Jpsikincons_y = -99.f;
  float Jpsikincons_z = -99.f;
  float Jpsikincons_pt = -99.f;
  float Jpsikincons_eta = -99.f;
  float Jpsikincons_phi = -99.f;
  float Jpsikincons_mass = -99.f;
  
  float Mupluskincons_pt = -99.f;
  float Mupluskincons_eta = -99.f;
  float Mupluskincons_phi = -99.f;
  
  float Muminuskincons_pt = -99.f;
  float Muminuskincons_eta = -99.f;
  float Muminuskincons_phi = -99.f;
  
  float Jpsigen_x;
  float Jpsigen_y;
  float Jpsigen_z;
  float Jpsigen_pt;
  float Jpsigen_eta;
  float Jpsigen_phi;
  float Jpsigen_mass;
  // ---- THE PRE-FSR RESONANCE ------------------------------------------
  //
  // `Jpsigen_mass` above is the POST-FSR pair: two status-1 muons matched in
  // dR < 0.1, i.e. what the tracker sees. The Z channel's FSR kernel is the
  // ratio of the two, so it needs the PRE-FSR mass as well, and today that
  // costs a separate FWLite pass over the MiniAOD (`zchannel/dump_gen_fsr.py`)
  // whose selection has to be kept in step with this one by hand.
  //
  // In the DY UL16 sample (powheg-MiNNLO + pythia8 + photos) NO muon carries
  // `fromHardProcessBeforeFSR`; what is there is the hard-process Z at
  // status 22 (first copy) and 62 (last copy), whose masses agree in
  // 4000/4000 events, and -- only in the 59 % of events that radiated -- the
  // pre-Photos muon pair at status 746. `m(mumu, 746) == m(Z, 62)` to an RMS
  // of 4e-6 GeV, so the status-62 resonance IS the pre-FSR mass and it exists
  // in every event (`zchannel/README.md`).
  //
  //   Jpsigenpre_mass    the |pdgId| in `genResonancePdgIds_` entry at
  //                      status 62, falling back to 22; -99 if absent (a
  //                      J/psi from a B decay has neither).
  //   Jpsigenpre_status  which copy was used, so a file says so itself.
  //   Jpsigenpre_masslep the status-746 lepton pair, the independent
  //                      cross-check; -99 when the event did not radiate.
  //   Jpsigen_massdressed the matched bare pair with every prompt status-1
  //                      photon within dR < 0.1 of either muon added back.
  float Jpsigenpre_mass;
  float Jpsigenpre_masslep;
  float Jpsigen_massdressed;
  int Jpsigenpre_status;
  std::vector<int> genResonancePdgIds_;
  
  float Muplusgen_pt;
  float Muplusgen_eta;
  float Muplusgen_phi;
  
  float Muminusgen_pt;
  float Muminusgen_eta;
  float Muminusgen_phi;

  float Muplusgen_dr;
  float Muminusgen_dr;

  // ---- GEN PROVENANCE OF THE TWO LEGS ----------------------------------
  //
  // `Mu*gen_dr` says only that SOME status-1 gen muon of the right charge sits
  // within dR < 0.1 of the leg. It does NOT say the candidate is a real
  // resonance decay: two reco tracks of the SAME muon, or a real muon paired
  // with a track from a heavy-flavour decay, both leave two "matched" legs.
  // Separating signal from combinatorial background needs the IDENTITY and the
  // ANCESTRY of the matched particle, which is what these carry. All are
  // filled only under doGen_, with the sentinels below otherwise.
  //
  //   Mu*gen_pdgId       pdgId of the matched status-1 gen particle (0 = none)
  //   Mu*gen_idx         its index in the gen collection (-1 = none). Two
  //                      candidates of one event sharing an index are two
  //                      reconstructions of ONE muon -- the combinatorial
  //                      background the inclusive vertex tail is made of.
  //   Mu*gen_motherPdgId pdgId of the first ancestor that is not itself a
  //                      muon copy: 23 for a Z daughter, 443/553 for a
  //                      quarkonium, a hadron id for a heavy-flavour decay
  //                      (0 = none found)
  //   Mu*gen_motherIdx   that ancestor's index (-1 = none)
  //   Mu*gen_isPrompt    GenStatusFlags::isPrompt(), asked through a
  //   Mu*gen_fromHardProcess  dynamic_cast (the collection is an
  //                      edm::View<reco::Candidate>, which does not expose the
  //                      flags) and left false when the cast fails
  //   Jpsigen_sameDecay  both legs matched, to DIFFERENT gen particles, whose
  //                      first non-muon ancestor is the SAME particle -- the
  //                      candidate is the two daughters of one decay
  int Muplusgen_pdgId;
  int Muminusgen_pdgId;
  int Muplusgen_idx;
  int Muminusgen_idx;
  int Muplusgen_motherPdgId;
  int Muminusgen_motherPdgId;
  int Muplusgen_motherIdx;
  int Muminusgen_motherIdx;
  bool Muplusgen_isPrompt;
  bool Muminusgen_isPrompt;
  bool Muplusgen_fromHardProcess;
  bool Muminusgen_fromHardProcess;
  bool Jpsigen_sameDecay;

  std::array<float, 3> Muplus_refParms;
  std::array<float, 3> Muminus_refParms;
  
  std::vector<float> Muplus_jacRef;
  std::vector<float> Muminus_jacRef;
  std::vector<float> Jpsi_jacMass;

  // ---- THE TWO LEGS' REFERENCE MOMENTUM COVARIANCE ----------------------
  //
  // `covrefmom` -- the (q/p, lambda, phi) covariance of BOTH legs at the
  // reference, the very matrix `Jpsi_sigmamass` is contracted out of -- is
  // exported so that the two second-order corrections the mass likelihood
  // needs are measured PER CANDIDATE rather than taken from MC:
  //
  //     A = sigma_rel1^2 + sigma_rel2^2 ,  B = 2 rho sigma_rel1 sigma_rel2
  //
  // is the Jensen term's whole content, and `f_ang`, the share of the mass
  // variance carried by the ANGLES rather than the two curvatures, is the
  // only thing the closed form `1.5 (sigma_m/m)^2` misses. On data there is
  // no MC to take them from, which is why the full symmetric 6x6 is written
  // out together with the two derived pieces that make it self-contained.
  //
  // ORDER IS (PLUS, MINUS), NOT the internal leg order: entries 0-2 are the
  // mu+ (q/p, lambda, phi) and 3-5 the mu-, permuted here by idxplus/idxminus
  // so that no consumer has to know `muchargearr`.
  //
  // `Jpsi_covrefmom` is the UPPER TRIANGLE, row-major: 21 floats, (0,0),
  // (0,1)...(0,5),(1,1)...(5,5). `Jpsi_jacrefmom` is dm/d(state) in the same
  // order, so `J C J^T` must reproduce `Jpsi_sigmamass^2` -- exported so the
  // matrix can be CHECKED rather than trusted.
  std::vector<float> Jpsi_covrefmom;
  std::vector<float> Jpsi_jacrefmom;
  // q/p of each leg at the reference, so sigma_rel = sqrt(C_ll)/|q/p| and the
  // sign conventions are fixed by the file rather than by a convention.
  float Jpsi_qoprefplus;
  float Jpsi_qoprefminus;
  // Derived, and the numbers the gates quote. sigma_rel is the RELATIVE
  // momentum resolution of the leg; `Jpsi_rhomom` is the leg-leg correlation
  // of d ln p (NOT of q/p: the two differ by sign(q+ q-), and it is d ln p
  // that enters B); `Jpsi_fang` = 1 - (J_kappa C J_kappa^T)/sigma_m^2 with
  // J_kappa the mass Jacobian restricted to the two q/p components.
  float Jpsi_sigmarelplus;
  float Jpsi_sigmarelminus;
  float Jpsi_rhomom;
  float Jpsi_fang;
  
  unsigned int Muplus_nhits;
  unsigned int Muplus_nvalid;
  unsigned int Muplus_nvalidpixel;
  unsigned int Muplus_nmatchedvalid;
  unsigned int Muplus_nambiguousmatchedvalid;
  // Per-muon counts of hits that survive the alignment-pass quality check
  // (`morehitquality` block at the per-hit loop). Symmetric with `nvalid`;
  // diverges from `nvalid` only if a future change introduces hit rejection
  // in `morehitquality` (currently always true).
  unsigned int Muplus_nvalidFinal;
  unsigned int Muplus_nvalidpixelFinal;

  unsigned int Muminus_nhits;
  unsigned int Muminus_nvalid;
  unsigned int Muminus_nvalidpixel;
  unsigned int Muminus_nmatchedvalid;
  unsigned int Muminus_nambiguousmatchedvalid;
  unsigned int Muminus_nvalidFinal;
  unsigned int Muminus_nvalidpixelFinal;

  // TRANSMISSION PROBE. Total MEAN energy loss that the fit's
  // reference trajectory actually applies between the reference point and the
  // outermost hit, summed over propagation steps of the LAST iteration of the
  // UNCONSTRAINED (icons==0) pass, in GeV. This is the denominator of the
  // "systematic transmission"
  //     T_sys = d(p_fit at reference) / d(assumed total energy loss)
  // measured by re-running with CVH_DEDX_SCALE != 1 and pairing per track.
  // It is the model's assumed loss, NOT the true (sim) loss: it is exactly
  // the quantity the coherent dE/dx re-centring moves.
  float Muplus_dEref;
  float Muminus_dEref;
  // THE WORST SINGLE PROPAGATION STEP OF THE LEG, as a fraction of its
  // momentum: max over surface-to-surface propagations of (E_in - E_out)/p_in.
  //
  // It is a QUALITY variable for the quadratic (hit-chi2) term's material
  // information, not a physics observable. The mean-loss bias of a group
  // grows with the step's fractional loss and is universal across groups
  // above dE/p ~ 0.1, while below 0.01 the pulls close (max/rms 0.96/0.19):
  // thick steps crossed by curlers with pT < 1 GeV carry 83 % of a group's
  // information and essentially all of its bias. With this and `_dEref` the
  // offline accumulation can impose a `dE_ref/p < 0.01` requirement WITHOUT
  // the step records, which is the whole point (they are 430 kB/candidate).
  float Muplus_maxfracloss;
  float Muminus_maxfracloss;

  // Per-muon pixel pathology-class counts of the hits USED in the fit
  // (16 combination bins, bit0=edgeX bit1=edgeY bit2=sizeX1 bit3=sizeY1;
  // bin 0 = clean), plus the count of demoted-to-inactive pixel hits.
  // Lets the analysis bin candidates by pathology content when the
  // veto is relaxed (keepPixelEdgeHits / pixelMinSizeX).
  std::vector<int> Muplus_pixClass;
  std::vector<int> Muminus_pixClass;
  unsigned int Muplus_npixDemoted;
  unsigned int Muminus_npixDemoted;

  bool Muplus_highpurity;
  bool Muminus_highpurity;

  // Track-level charges of the two daughters at idxplus/idxminus. The
  // names "plus"/"minus" reflect the candidate-builder ordering convention
  // (positive-charge first; daughter[0] = +q, daughter[1] = -q for opposite-
  // sign pairs); these branches expose the actual track charges so
  // downstream plotting can split by charge if needed.
  int Muplus_charge;
  int Muminus_charge;

  // Per-track dE/dx scalar estimators (Harmonic2 strip / pixel / joint),
  // projected onto the ALCARECO selected-track collection by
  // DeDxValueMapProjector in the skim and looked up via the candidate's
  // daughter TrackRef (which already points at the persisted ALCARECO
  // collection). Filled only when all four dE/dx parameters are configured
  // with non-empty InputTags (`readDeDx_`).
  bool readDeDx_;
  edm::EDGetTokenT<reco::TrackCollection> dedxSourceTracksToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> dedxHarmonic2Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> dedxPixelHarmonic2Token_;
  edm::EDGetTokenT<edm::ValueMap<float>> dedxAllHarmonic2Token_;
  float Muplus_dedxHarmonic2;
  float Muplus_dedxPixelHarmonic2;
  float Muplus_dedxAllHarmonic2;
  float Muminus_dedxHarmonic2;
  float Muminus_dedxPixelHarmonic2;
  float Muminus_dedxAllHarmonic2;

  bool Muplus_isMuon;
  bool Muplus_muonLoose;
  bool Muplus_muonMedium;
  bool Muplus_muonTight;
  bool Muplus_muonIsPF;
  bool Muplus_muonIsTracker;
  bool Muplus_muonIsGlobal;
  bool Muplus_muonIsStandalone;
  bool Muplus_muonInnerTrackBest;

  bool Muminus_isMuon;
  bool Muminus_muonLoose;
  bool Muminus_muonMedium;
  bool Muminus_muonTight;
  bool Muminus_muonIsPF;
  bool Muminus_muonIsTracker;
  bool Muminus_muonIsGlobal;
  bool Muminus_muonIsStandalone;
  bool Muminus_muonInnerTrackBest;

  float edmval_cons0;
  int niter_cons0;

  float dmassconvval = 0.;
  float dinvmasssqconvval = 0.;
  
  float dmassconvval_cons0 = 0.;
  float dinvmasssqconvval_cons0 = 0.;

// std::vector<float> hessv;

  edm::ESGetToken<TransientTrackingRecHitBuilder, TransientRecHitRecord> ttrhToken_;
  edm::ESGetToken<Propagator, TrackingComponentsRecord> g4ePropToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transTrackBuilderToken_;
  edm::ESGetToken<L1GtTriggerMenu, L1GtTriggerMenuRcd> l1MenuToken_;
  // The StripCPEfromTrackAngle instance the cloner uses, queried only for its
  // AlgoParam so the hit-resolution blocks can be labelled with the CPE's own
  // uProj -- the strip class variable (`hitres_classes.class_of`).
  edm::ESGetToken<StripClusterParameterEstimator, TkStripCPERecord> stripCPEToken_;

};


ResidualGlobalCorrectionMakerTwoTrackG4e::ResidualGlobalCorrectionMakerTwoTrackG4e(
    const edm::ParameterSet &iConfig)
    : ResidualGlobalCorrectionMakerBase(iConfig),
      ttrhToken_(esConsumes(edm::ESInputTag("", "WithAngleAndTemplate"))),
      g4ePropToken_(esConsumes(edm::ESInputTag("", "Geant4ePropagator"))),
      transTrackBuilderToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
      l1MenuToken_(esConsumes()),
      stripCPEToken_(esConsumes(edm::ESInputTag("", "StripCPEfromTrackAngle")))
{
  // The resolution-CF exponents this maker exports are the CANDIDATE-MASS
  // functional, not the single-track q/p one: different standardization,
  // different ionization sign. They must not share a branch name with it.
  cfprefix_ = "cfmass";
  cfGroupDelta_ = false;
  // Which |pdgId| counts as "the resonance" for the pre-FSR gen mass. The Z
  // is 23, prompt charmonium 443/100443, bottomonium 553/100553/200553.
  genResonancePdgIds_ = iConfig.existsAs<std::vector<int>>("genResonancePdgIds")
      ? iConfig.getParameter<std::vector<int>>("genResonancePdgIds")
      : std::vector<int>{23, 443, 100443, 553, 100553, 200553};
  doVtxConstraint_ = iConfig.getParameter<bool>("doVtxConstraint");
  doMassConstraint_ = iConfig.getParameter<bool>("doMassConstraint");
  // minimum size of a pair (see the member docs). `existsAs` so that every
  // existing cfi keeps working; the DEFAULTS are the cut, not a no-op.
  minNdof_ = iConfig.existsAs<int>("minNdof") ? iConfig.getParameter<int>("minNdof") : 1;
  minPairHits_ = iConfig.existsAs<int>("minPairHits")
                     ? iConfig.getParameter<int>("minPairHits") : -1;
  if (minPairHits_ < 0) {
    minPairHits_ = doVtxConstraint_ ? 10 : 11;
  }
  minLegHits_ = iConfig.existsAs<int>("minLegHits")
                    ? iConfig.getParameter<int>("minLegHits") : 8;
  // THE VERTEX-CONSTRAINT RESIDUAL, off by default so that no existing
  // configuration changes its output by a byte.
  exportVtxResidual_ = iConfig.existsAs<bool>("exportVtxResidual")
                           ? iConfig.getParameter<bool>("exportVtxResidual") : false;
  // Say out loud what the log-det term is doing, and to which families: it
  // changes the meaning of `gradv`/`hesspackedv`/`hessfactorv` and, for the
  // per-module families, the LAYOUT of `globalidxv`.
  if (exportVarianceGrads_) {
    std::string fams;
    if (varianceGradFamilies_.empty()) {
      fams = "8,9,10,11,15 (default)";
    } else {
      for (unsigned int f : varianceGradFamilies_) {
        fams += (fams.empty() ? "" : ",") + std::to_string(f);
      }
    }
    const bool wantHit = varianceFamilyWanted(8) || varianceFamilyWanted(9);
    const bool wantMat = varianceFamilyWanted(10) || varianceFamilyWanted(11);
    std::cout << "ResidualGlobalCorrectionMakerTwoTrackG4e: exportVarianceGrads ON, "
              << "families " << fams << ".  gradv/hess now differentiate the "
              << "COVARIANCE as well as the mean; the Hessian is the EXPECTED "
              << "(Fisher) one, so the mean-variance cross term is zero by "
              << "construction." << std::endl;
    if (wantHit || wantMat) {
      std::cout << "  NOTE: families 8/9/10/11 are per-module and are NOT columns "
                << "of this maker's parameter vector, so globalidxv/gradv/"
                << "hesspackedv/hessfactorv/jacrefv/Jpsi_jacMass GROW.  Use "
                << "varianceGradFamilies=15 to leave the layout unchanged."
                << std::endl;
    }
    if (wantHit) {
      std::cout << "  NOTE: this maker does not apply exp(corparms) to the hit "
                << "covariance (the single-track one does), so the parmtype-8/9 "
                << "derivatives are evaluated at k = 0 whatever corFiles says."
                << std::endl;
    }
    if (varianceFamilyWanted(15) && !exportMaterialNoise_) {
      std::cout << "  WARNING: family 15 requested but exportMaterialNoise=False, "
                << "so no parmtype-15 dV block is registered and k_g keeps its "
                << "mean-loss-only derivative." << std::endl;
    }
  }
  massConstraint_ = iConfig.getParameter<double>("massConstraint");
  massConstraintWidth_ = iConfig.getParameter<double>("massConstraintWidth");
  // Per-daughter Geant4 particle-name base. Channel cfis specify this
  // (e.g. "pi", "proton", "kaon") and the corresponding mass + mass
  // uncertainty are looked up from a single PDG table in
  // Analysis/HitAnalyzer/interface/ParticleProperties.h. Default empty
  // -> falls back to muon (legacy J/psi/Upsilon configuration).
  daughterParticleName1_ = iConfig.existsAs<std::string>("daughterParticleName1")
      ? iConfig.getParameter<std::string>("daughterParticleName1") : std::string("mu");
  daughterParticleName2_ = iConfig.existsAs<std::string>("daughterParticleName2")
      ? iConfig.getParameter<std::string>("daughterParticleName2") : std::string("mu");
  const auto props1 = ana_hitanalyzer::getParticleProperties(daughterParticleName1_);
  const auto props2 = ana_hitanalyzer::getParticleProperties(daughterParticleName2_);
  daughterMass1_    = props1.mass;
  daughterMass2_    = props2.mass;
  daughterMass1Err_ = props1.massErr;
  daughterMass2Err_ = props2.massErr;
  // 2D-transverse pointing-angle constraint (V0 channels). Default off.
  doPointingConstraint_ = iConfig.existsAs<bool>("doPointingConstraint")
      ? iConfig.getParameter<bool>("doPointingConstraint") : false;
  // Angular pointing resolution width, radians; ~1 mrad matches V0Producer's
  // cosThetaXY > 0.998 cut converted to a 1-sigma soft constraint.
  pointingSigma_ = iConfig.existsAs<double>("pointingSigma")
      ? iConfig.getParameter<double>("pointingSigma") : 1.e-3;

  // CVH-refit convergence knobs (openspec/improve-cvh-refit-convergence).
  // Defaults reproduce the published Run2016H baseline bit-identically.
  nIters_ = iConfig.existsAs<unsigned int>("nIters")
      ? iConfig.getParameter<unsigned int>("nIters") : 10u;
  edmConvergence_ = iConfig.existsAs<double>("edmConvergence")
      ? iConfig.getParameter<double>("edmConvergence") : 1.e-5;
  clampMomentumFloor_ = iConfig.existsAs<double>("clampMomentumFloor")
      ? iConfig.getParameter<double>("clampMomentumFloor") : 2.0;
  // Echo it once per maker instance: the floor silently decides whether soft
  // daughters are fitted or pinned at it, and a job log must record which
  // value was in force. Same line as the single-track maker.
  maxMomentumStepFactor_ = iConfig.existsAs<double>("maxMomentumStepFactor")
      ? iConfig.getParameter<double>("maxMomentumStepFactor") : 2.0;
  stepBacktracking_ = iConfig.existsAs<bool>("stepBacktracking")
      ? iConfig.getParameter<bool>("stepBacktracking") : true;
  maxChi2Backtrack_ = iConfig.existsAs<unsigned int>("maxChi2Backtrack")
      ? iConfig.getParameter<unsigned int>("maxChi2Backtrack") : 4u;
  stepBacktrackFromIter_ = iConfig.existsAs<unsigned int>("stepBacktrackFromIter")
      ? iConfig.getParameter<unsigned int>("stepBacktrackFromIter") : 2u;
  armijoC_ = iConfig.existsAs<double>("armijoC")
      ? iConfig.getParameter<double>("armijoC") : 1.e-4;
  armijoSlack_ = iConfig.existsAs<double>("armijoSlack")
      ? iConfig.getParameter<double>("armijoSlack") : 1.0;
  stepPrintLimit_ = iConfig.existsAs<unsigned int>("stepPrintLimit")
      ? iConfig.getParameter<unsigned int>("stepPrintLimit") : 200u;
  edm::LogPrint("ResidualGlobalCorrectionMakerTwoTrackG4e")
      << "[cvh] effective: clampMomentumFloor=" << clampMomentumFloor_ << " GeV"
      << ", maxMomentumStepFactor=" << maxMomentumStepFactor_
      << ", stepBacktracking=" << stepBacktracking_
      << " (fromIter=" << stepBacktrackFromIter_
      << ", maxChi2Backtrack=" << maxChi2Backtrack_
      << ", armijoC=" << armijoC_ << ", armijoSlack=" << armijoSlack_ << ")";
  maxBacktracks_ = iConfig.existsAs<unsigned int>("maxBacktracks")
      ? iConfig.getParameter<unsigned int>("maxBacktracks") : 4u;
  maxSeedInflations_ = iConfig.existsAs<unsigned int>("maxSeedInflations")
      ? iConfig.getParameter<unsigned int>("maxSeedInflations") : 2u;

  // Which two-track subsystem of each candidate to fit (see member comment).
  // Default -1 keeps the legacy flat-candidate behaviour.
  subsystemDaughter_ = iConfig.existsAs<int>("subsystemDaughter")
      ? iConfig.getParameter<int>("subsystemDaughter") : -1;
  
  // Pixel-pathology bias attribution (see member docs above).
  fillHitDiagnostics_ = iConfig.existsAs<bool>("fillHitDiagnostics")
      ? iConfig.getParameter<bool>("fillHitDiagnostics") : false;
  deweightPathoHits_ = iConfig.existsAs<bool>("deweightPathoHits")
      ? iConfig.getParameter<bool>("deweightPathoHits") : false;

  // NanoAOD path: emit per-candidate ValueMaps keyed to srcCandidates. Only
  // meaningful in the candidate-driven mode. Kinematics always; the global-fit
  // payload additionally when fillGradsFactored_ is set.
  produceValueMaps_ = iConfig.existsAs<bool>("produceValueMaps")
      ? iConfig.getParameter<bool>("produceValueMaps") : false;
  if (produceValueMaps_) {
    vmCorMass_    = produces<edm::ValueMap<float>>("corMass");
    vmCorMassErr_ = produces<edm::ValueMap<float>>("corMassErr");
    vmCorPt_      = produces<edm::ValueMap<float>>("corPt");
    vmCorEta_     = produces<edm::ValueMap<float>>("corEta");
    vmCorPhi_     = produces<edm::ValueMap<float>>("corPhi");
    vmMuPlusPt_   = produces<edm::ValueMap<float>>("muPlusPt");
    vmMuPlusEta_  = produces<edm::ValueMap<float>>("muPlusEta");
    vmMuPlusPhi_  = produces<edm::ValueMap<float>>("muPlusPhi");
    vmMuMinusPt_  = produces<edm::ValueMap<float>>("muMinusPt");
    vmMuMinusEta_ = produces<edm::ValueMap<float>>("muMinusEta");
    vmMuMinusPhi_ = produces<edm::ValueMap<float>>("muMinusPhi");
    vmEdmval_     = produces<edm::ValueMap<float>>("edmval");
    if (fillGradsFactored_) {
      vmGlobalIdxs_    = produces<edm::ValueMap<std::vector<int>>>("globalIdxs");
      vmJacRefMuPlus_  = produces<edm::ValueMap<std::vector<float>>>("jacRefMuPlus");
      vmJacRefMuMinus_ = produces<edm::ValueMap<std::vector<float>>>("jacRefMuMinus");
      vmJacMass_       = produces<edm::ValueMap<std::vector<float>>>("jacMass");
      vmHessFactor_    = produces<edm::ValueMap<std::vector<float>>>("hessFactor");
    }
  }

  useStartingState_ = iConfig.existsAs<std::string>("useStartingState")
      ? iConfig.getParameter<std::string>("useStartingState") : std::string("perigee");
  if (useStartingState_ != "perigee" && useStartingState_ != "midPropagated") {
    throw cms::Exception("Configuration")
        << "ResidualGlobalCorrectionMakerTwoTrackG4e: useStartingState='"
        << useStartingState_ << "' not supported. "
        << "Valid values are: 'perigee', 'midPropagated'.";
  }
  if (useStartingState_ == "midPropagated") {
    // Implementation per openspec/fix-cvh-displaced-starting-state:
    // each daughter's perigee FreeTrajectoryState is extrapolated (analytical
    // helix, AnalyticalImpactPointExtrapolator) to the 3D closest approach
    // to the midpoint of the two daughter perigees. The resulting TSOS is
    // used as the iter-0 reference state instead of the in-maker Kalman
    // fit's daughter-state-at-vertex. Motivation: in-maker Kalman fit on
    // prompt-tracking-perigee inputs can land mm-biased for displaced
    // J/psi-from-B vertices, putting iter-0 outside the GN convergence
    // basin (see the K_S vs J/psi residual diagnostic). The geometric
    // perigee-midpoint is a cheap displacement-aware seed that K_S
    // gets implicitly from V0Producer; midPropagated grants the same
    // to J/psi without an upstream Kalman fit.
    edm::LogInfo("ResidualGlobalCorrectionMakerTwoTrackG4e")
        << "useStartingState='midPropagated' active. Iter-0 reference state "
        << "for the joint refit will be the AnalyticalImpactPointExtrapolator "
        << "output at the perigee-midpoint. On extrapolation failure the "
        << "producer falls back to 'perigee' per-event and increments "
        << "midPropagatedFallbackCount_.";
  }
  debugPerIterDump_ = iConfig.existsAs<bool>("debugPerIterDump")
      ? iConfig.getParameter<bool>("debugPerIterDump") : false;
  edm::LogInfo("ResidualGlobalCorrectionMakerTwoTrackG4e")
      << "CVH convergence knobs: nIters=" << nIters_
      << ", edmConvergence=" << edmConvergence_
      << ", useStartingState=" << useStartingState_
      << ", debugPerIterDump=" << (debugPerIterDump_ ? "true" : "false");

  // Optional per-track dE/dx ValueMaps (Harmonic2 strip + pixel-only),
  // projected onto the ALCARECO selected-track collection by
  // DeDxValueMapProjector. All three of `dedxSourceTracks` (the ALCARECO
  // collection the ValueMaps are keyed on), `dedxHarmonic2`, and
  // `dedxPixelHarmonic2` must be set to enable the lookup; legacy
  // J/psi / Upsilon runners that don't set them get default-off.
  // Require all four to be present AND have non-empty labels. An empty
  // InputTag passes existsAs<> but fails at getByToken with ProductNotFound,
  // so use the label as the gate -- letting dimuon configs declare the
  // parameters as empty InputTag('') without enabling the lookup.
  auto hasLabel = [&](const char* name) {
    return iConfig.existsAs<edm::InputTag>(name) &&
           !iConfig.getParameter<edm::InputTag>(name).label().empty();
  };
  readDeDx_ = hasLabel("dedxSourceTracks") && hasLabel("dedxHarmonic2") &&
              hasLabel("dedxPixelHarmonic2") && hasLabel("dedxAllHarmonic2");
  if (readDeDx_) {
    dedxSourceTracksToken_   = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("dedxSourceTracks"));
    dedxHarmonic2Token_      = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("dedxHarmonic2"));
    dedxPixelHarmonic2Token_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("dedxPixelHarmonic2"));
    dedxAllHarmonic2Token_   = mayConsume<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("dedxAllHarmonic2"));
  }

  doL1Trigger_ = iConfig.existsAs<bool>("doL1Trigger") ? iConfig.getParameter<bool>("doL1Trigger") : false;
  if (doL1Trigger_) {
    const edm::InputTag l1Results = iConfig.existsAs<edm::InputTag>("l1Results")
                                        ? iConfig.getParameter<edm::InputTag>("l1Results")
                                        : edm::InputTag("gtDigis", "", "RECO");
    inputL1ReadoutRecordTag_ = l1Results;
    inputL1ReadoutRecord_ = consumes<L1GlobalTriggerReadoutRecord>(l1Results);
    l1Triggers_ = iConfig.getParameter<std::vector<std::string>>("l1Triggers");
    l1TriggerBitNumbers_.assign(l1Triggers_.size(), std::numeric_limits<unsigned int>::max());
    l1TriggerDecisions_.assign(l1Triggers_.size(), 0);
  }
}

void ResidualGlobalCorrectionMakerTwoTrackG4e::beginStream(edm::StreamID streamid)
{
  ResidualGlobalCorrectionMakerBase::beginStream(streamid);
  
  if (fillTrackTree_) {
    // Stage-2 per-row B+ candidate index, branched only when the cfi
    // configured bCandIdxSrc (additive, no-op for legacy J/psi/Upsilon/Z).
    if (!bCandIdxSrcTag_.label().empty()) {
      tree->Branch("bCandIdx", &bCandIdx);
    }
    // From CVH refit (no mass constraint): per-track parameters at the
    // joint two-track PCA from the GBL/Geant4e fit. Dimuon kinematics are
    // the sum of the per-track 4-vectors. Vertex (x, y, z) and signed
    // displacement d come from the joint-PCA state vector statepcaupd
    // (sign of d is tied to the charge of track[0]).
    tree->Branch("Jpsi_d", &Jpsi_d);
    tree->Branch("Jpsi_x", &Jpsi_x);
    tree->Branch("Jpsi_y", &Jpsi_y);
    tree->Branch("Jpsi_z", &Jpsi_z);
    tree->Branch("Jpsi_pt", &Jpsi_pt);
    tree->Branch("Jpsi_eta", &Jpsi_eta);
    tree->Branch("Jpsi_phi", &Jpsi_phi);
    tree->Branch("Jpsi_mass", &Jpsi_mass);
    // Per-event dimuon-mass uncertainty propagated from the CVH covariance.
    tree->Branch("Jpsi_sigmamass", &Jpsi_sigmamass);
    // The two legs' reference momentum covariance and the mass Jacobian it is
    // contracted with -- 33 floats, 132 B/candidate, 0.16 % of the slim
    // record. Always on: without them the Jensen and self-consistent-sigma
    // corrections have no truth-free inputs on data.
    tree->Branch("Jpsi_covrefmom", &Jpsi_covrefmom);
    tree->Branch("Jpsi_jacrefmom", &Jpsi_jacrefmom);
    tree->Branch("Jpsi_qoprefplus", &Jpsi_qoprefplus);
    tree->Branch("Jpsi_qoprefminus", &Jpsi_qoprefminus);
    tree->Branch("Jpsi_sigmarelplus", &Jpsi_sigmarelplus);
    tree->Branch("Jpsi_sigmarelminus", &Jpsi_sigmarelminus);
    tree->Branch("Jpsi_rhomom", &Jpsi_rhomom);
    tree->Branch("Jpsi_fang", &Jpsi_fang);

    tree->Branch("Muplus_pt", &Muplus_pt);
    tree->Branch("Muplus_eta", &Muplus_eta);
    tree->Branch("Muplus_phi", &Muplus_phi);

    tree->Branch("Muminus_pt", &Muminus_pt);
    tree->Branch("Muminus_eta", &Muminus_eta);
    tree->Branch("Muminus_phi", &Muminus_phi);

    // From CMSSW KinematicParticleVertexFitter (no mass constraint),
    // run on the post-CVH per-track states. This is the standard CMSSW
    // kinematic vertex fit -- distinct from the CVH/GBL refit above --
    // and provides a vertex-constrained per-particle 4-momentum at the
    // fitted decay vertex. Used as a cross-check of the CVH state.
    tree->Branch("Jpsikin_x", &Jpsikin_x);
    tree->Branch("Jpsikin_y", &Jpsikin_y);
    tree->Branch("Jpsikin_z", &Jpsikin_z);
    tree->Branch("Jpsikin_pt", &Jpsikin_pt);
    tree->Branch("Jpsikin_eta", &Jpsikin_eta);
    tree->Branch("Jpsikin_phi", &Jpsikin_phi);
    tree->Branch("Jpsikin_mass", &Jpsikin_mass);
    
    tree->Branch("Mupluskin_pt", &Mupluskin_pt);
    tree->Branch("Mupluskin_eta", &Mupluskin_eta);
    tree->Branch("Mupluskin_phi", &Mupluskin_phi);
    
    tree->Branch("Muminuskin_pt", &Muminuskin_pt);
    tree->Branch("Muminuskin_eta", &Muminuskin_eta);
    tree->Branch("Muminuskin_phi", &Muminuskin_phi);
    
    // From the raw input tracks, before any refit: per-track 4-momenta
    // are built from the input track px/py/pz and the daughter mass
    // hypothesis (trackMass[0], trackMass[1] from the channel cfi).
    // The dimu mass is just the 4-vector sum -- no vertex or mass
    // constraint applied. Useful as a no-fit reference / sanity check.
    tree->Branch("Jpsitrk_pt", &Jpsitrk_pt);
    tree->Branch("Jpsitrk_eta", &Jpsitrk_eta);
    tree->Branch("Jpsitrk_phi", &Jpsitrk_phi);
    tree->Branch("Jpsitrk_mass", &Jpsitrk_mass);
    
    tree->Branch("Muplustrk_pt", &Muplustrk_pt);
    tree->Branch("Muplustrk_eta", &Muplustrk_eta);
    tree->Branch("Muplustrk_phi", &Muplustrk_phi);
    
    tree->Branch("Muminustrk_pt", &Muminustrk_pt);
    tree->Branch("Muminustrk_eta", &Muminustrk_eta);
    tree->Branch("Muminustrk_phi", &Muminustrk_phi);
    
    // From CVH refit WITH mass constraint: same content as the no-suffix
    // Jpsi/Muplus/Muminus block above, but for the icons==1 pass where a
    // dimuon-mass-constraint chi^2 term is added to the GBL fit. Only
    // filled when doMassConstraint=True is set on the channel cfi
    // (currently False for KS, Lambda, D0; True for Jpsi).
    tree->Branch("Jpsicons_d", &Jpsicons_d);
    tree->Branch("Jpsicons_x", &Jpsicons_x);
    tree->Branch("Jpsicons_y", &Jpsicons_y);
    tree->Branch("Jpsicons_z", &Jpsicons_z);
    tree->Branch("Jpsicons_pt", &Jpsicons_pt);
    tree->Branch("Jpsicons_eta", &Jpsicons_eta);
    tree->Branch("Jpsicons_phi", &Jpsicons_phi);
    tree->Branch("Jpsicons_mass", &Jpsicons_mass);

    tree->Branch("Mupluscons_pt", &Mupluscons_pt);
    tree->Branch("Mupluscons_eta", &Mupluscons_eta);
    tree->Branch("Mupluscons_phi", &Mupluscons_phi);

    tree->Branch("Muminuscons_pt", &Muminuscons_pt);
    tree->Branch("Muminuscons_eta", &Muminuscons_eta);
    tree->Branch("Muminuscons_phi", &Muminuscons_phi);
    
    // From CMSSW KinematicConstrainedVertexFitter (with a TwoTrackMass
    // constraint) on the post-CVH per-track states. As above, distinct
    // from the CVH cons block: this is the standard CMSSW kinematic
    // vertex fit with the dimuon mass constrained to the channel's
    // expectedMass. Only filled when doMassConstraint=True.
    tree->Branch("Jpsikincons_x", &Jpsikincons_x);
    tree->Branch("Jpsikincons_y", &Jpsikincons_y);
    tree->Branch("Jpsikincons_z", &Jpsikincons_z);
    tree->Branch("Jpsikincons_pt", &Jpsikincons_pt);
    tree->Branch("Jpsikincons_eta", &Jpsikincons_eta);
    tree->Branch("Jpsikincons_phi", &Jpsikincons_phi);
    tree->Branch("Jpsikincons_mass", &Jpsikincons_mass);
    
    tree->Branch("Mupluskincons_pt", &Mupluskincons_pt);
    tree->Branch("Mupluskincons_eta", &Mupluskincons_eta);
    tree->Branch("Mupluskincons_phi", &Mupluskincons_phi);
    
    tree->Branch("Muminuskincons_pt", &Muminuskincons_pt);
    tree->Branch("Muminuskincons_eta", &Muminuskincons_eta);
    tree->Branch("Muminuskincons_phi", &Muminuskincons_phi);
    
    // Generator-level (MC truth, only filled when doGen_=True).
    // Mu*gen_dr = deltaR matching distance between the gen-muon and the
    // reco track (cut at 0.1 in the lookup; -1 if no match).
    tree->Branch("Jpsigen_x", &Jpsigen_x);
    tree->Branch("Jpsigen_y", &Jpsigen_y);
    tree->Branch("Jpsigen_z", &Jpsigen_z);
    tree->Branch("Jpsigen_pt", &Jpsigen_pt);
    tree->Branch("Jpsigen_eta", &Jpsigen_eta);
    tree->Branch("Jpsigen_phi", &Jpsigen_phi);
    tree->Branch("Jpsigen_mass", &Jpsigen_mass);
    tree->Branch("Jpsigenpre_mass", &Jpsigenpre_mass);
    tree->Branch("Jpsigenpre_masslep", &Jpsigenpre_masslep);
    tree->Branch("Jpsigenpre_status", &Jpsigenpre_status);
    tree->Branch("Jpsigen_massdressed", &Jpsigen_massdressed);

    tree->Branch("Muplusgen_pt", &Muplusgen_pt);
    tree->Branch("Muplusgen_eta", &Muplusgen_eta);
    tree->Branch("Muplusgen_phi", &Muplusgen_phi);

    tree->Branch("Muminusgen_pt", &Muminusgen_pt);
    tree->Branch("Muminusgen_eta", &Muminusgen_eta);
    tree->Branch("Muminusgen_phi", &Muminusgen_phi);

    tree->Branch("Muplusgen_dr", &Muplusgen_dr);
    tree->Branch("Muminusgen_dr", &Muminusgen_dr);

    // gen PROVENANCE (see the member docs): what the leg was matched TO, and
    // whether the two legs are the two daughters of one decay.
    tree->Branch("Muplusgen_pdgId", &Muplusgen_pdgId);
    tree->Branch("Muminusgen_pdgId", &Muminusgen_pdgId);
    tree->Branch("Muplusgen_idx", &Muplusgen_idx);
    tree->Branch("Muminusgen_idx", &Muminusgen_idx);
    tree->Branch("Muplusgen_motherPdgId", &Muplusgen_motherPdgId);
    tree->Branch("Muminusgen_motherPdgId", &Muminusgen_motherPdgId);
    tree->Branch("Muplusgen_motherIdx", &Muplusgen_motherIdx);
    tree->Branch("Muminusgen_motherIdx", &Muminusgen_motherIdx);
    tree->Branch("Muplusgen_isPrompt", &Muplusgen_isPrompt);
    tree->Branch("Muminusgen_isPrompt", &Muminusgen_isPrompt);
    tree->Branch("Muplusgen_fromHardProcess", &Muplusgen_fromHardProcess);
    tree->Branch("Muminusgen_fromHardProcess", &Muminusgen_fromHardProcess);
    tree->Branch("Jpsigen_sameDecay", &Jpsigen_sameDecay);
    
    // Per-track reference parameters at PCA (3-vector: q/pT, lambda, phi)
    // and Jacobians of those parameters and of the dimuon mass with
    // respect to the global correction parameters (alignment, B-field,
    // material). Used to re-apply updated calibrations downstream
    // without re-running the CVH fit. fillJac_ controls the (large)
    // Jacobian branches.
    tree->Branch("Muplus_refParms", Muplus_refParms.data(), "Muplus_refParms[3]/F");
    tree->Branch("Muminus_refParms", Muminus_refParms.data(), "Muminus_refParms[3]/F");

    if (fillJac_) {
      tree->Branch("Muplus_jacRef", &Muplus_jacRef);
      tree->Branch("Muminus_jacRef", &Muminus_jacRef);
      tree->Branch("Jpsi_jacMass", &Jpsi_jacMass);
    }
    
    // Per-track hit-content counters (from the original reco::Track)
    // and quality flag.
    // nhits = total rec-hits on track
    // nvalid = valid (non-rejected) hits
    // nvalidpixel = subset on pixel detectors
    // nmatchedvalid = valid hits compatible with the CVH fit
    // nambiguousmatchedvalid = valid hits with multiple compatible matches
    // highpurity = passes the standard high-purity track selection
    tree->Branch("Muplus_nhits", &Muplus_nhits);
    tree->Branch("Muplus_nvalid", &Muplus_nvalid);
    tree->Branch("Muplus_nvalidpixel", &Muplus_nvalidpixel);
    tree->Branch("Muplus_nvalidFinal", &Muplus_nvalidFinal);
    tree->Branch("Muplus_nvalidpixelFinal", &Muplus_nvalidpixelFinal);
    tree->Branch("Muplus_nmatchedvalid", &Muplus_nmatchedvalid);
    tree->Branch("Muplus_nambiguousmatchedvalid", &Muplus_nambiguousmatchedvalid);

    tree->Branch("Muminus_nhits", &Muminus_nhits);
    tree->Branch("Muminus_nvalid", &Muminus_nvalid);
    tree->Branch("Muminus_nvalidpixel", &Muminus_nvalidpixel);
    tree->Branch("Muminus_nvalidFinal", &Muminus_nvalidFinal);
    tree->Branch("Muminus_nvalidpixelFinal", &Muminus_nvalidpixelFinal);
    tree->Branch("Muminus_nmatchedvalid", &Muminus_nmatchedvalid);
    tree->Branch("Muminus_nambiguousmatchedvalid", &Muminus_nambiguousmatchedvalid);

    // transmission probe: assumed total mean eloss along the reference
    // trajectory (GeV), icons==0, last iteration. See member declaration.
    tree->Branch("Muplus_dEref", &Muplus_dEref);
    tree->Branch("Muminus_dEref", &Muminus_dEref);
    tree->Branch("Muplus_maxfracloss", &Muplus_maxfracloss);
    tree->Branch("Muminus_maxfracloss", &Muminus_maxfracloss);

    tree->Branch("Muplus_pixClass", &Muplus_pixClass);
    tree->Branch("Muminus_pixClass", &Muminus_pixClass);
    tree->Branch("Muplus_npixDemoted", &Muplus_npixDemoted);
    tree->Branch("Muminus_npixDemoted", &Muminus_npixDemoted);

    if (fillHitDiagnostics_) {
      tree->Branch("hitdiag_detid", &hitdiag_detid);
      tree->Branch("hitdiag_trk", &hitdiag_trk);
      tree->Branch("hitdiag_charge", &hitdiag_charge);
      tree->Branch("hitdiag_class", &hitdiag_class);
      tree->Branch("hitdiag_dx", &hitdiag_dx);
      tree->Branch("hitdiag_dy", &hitdiag_dy);
      tree->Branch("hitdiag_exx", &hitdiag_exx);
      tree->Branch("hitdiag_eyy", &hitdiag_eyy);
      tree->Branch("hitdiag_lx", &hitdiag_lx);
      tree->Branch("hitdiag_ly", &hitdiag_ly);
    }

    tree->Branch("Muplus_highpurity", &Muplus_highpurity);
    tree->Branch("Muminus_highpurity", &Muminus_highpurity);

    // Per-track charge.
    tree->Branch("Muplus_charge", &Muplus_charge);
    tree->Branch("Muminus_charge", &Muminus_charge);

    // Per-track dE/dx Harmonic-2 truncated estimators projected from
    // ALCARECO via DeDxValueMapProjector (re-keyed onto the selected
    // track collection). Three flavours: strip-only, pixel-only, and
    // joint strip+pixel. Optional via the readDeDx_ knob; values are
    // NaN if the corresponding ValueMap is not present in the input.
    if (readDeDx_) {
      tree->Branch("Muplus_dedxHarmonic2",       &Muplus_dedxHarmonic2);
      tree->Branch("Muplus_dedxPixelHarmonic2",  &Muplus_dedxPixelHarmonic2);
      tree->Branch("Muplus_dedxAllHarmonic2",    &Muplus_dedxAllHarmonic2);
      tree->Branch("Muminus_dedxHarmonic2",      &Muminus_dedxHarmonic2);
      tree->Branch("Muminus_dedxPixelHarmonic2", &Muminus_dedxPixelHarmonic2);
      tree->Branch("Muminus_dedxAllHarmonic2",   &Muminus_dedxAllHarmonic2);
    }

    // Muon-ID flags from a reco::Muon match against the track:
    // isMuon = a reco::Muon was found pointing at this track
    // muonLoose/Medium/Tight = standard CMS muon-ID working points
    // muonIsPF = particle-flow identified
    // muonIsTracker = has a tracker-only segment
    // muonIsGlobal = has a global (tracker+muon-system) fit
    // muonIsStandalone = has a stand-alone muon-system fit
    // muonInnerTrackBest = the matched reco::Muon's innerTrack is
    // the same as the input track (best-track
    // pointer comparison)
    tree->Branch("Muplus_isMuon", &Muplus_isMuon);
    tree->Branch("Muplus_muonLoose", &Muplus_muonLoose);
    tree->Branch("Muplus_muonMedium", &Muplus_muonMedium);
    tree->Branch("Muplus_muonTight", &Muplus_muonTight);
    tree->Branch("Muplus_muonIsPF", &Muplus_muonIsPF);
    tree->Branch("Muplus_muonIsTracker", &Muplus_muonIsTracker);
    tree->Branch("Muplus_muonIsGlobal", &Muplus_muonIsGlobal);
    tree->Branch("Muplus_muonIsStandalone", &Muplus_muonIsStandalone);
    tree->Branch("Muplus_muonInnerTrackBest", &Muplus_muonInnerTrackBest);

    tree->Branch("Muminus_isMuon", &Muminus_isMuon);
    tree->Branch("Muminus_muonLoose", &Muminus_muonLoose);
    tree->Branch("Muminus_muonMedium", &Muminus_muonMedium);
    tree->Branch("Muminus_muonTight", &Muminus_muonTight);
    tree->Branch("Muminus_muonIsPF", &Muminus_muonIsPF);
    tree->Branch("Muminus_muonIsTracker", &Muminus_muonIsTracker);
    tree->Branch("Muminus_muonIsGlobal", &Muminus_muonIsGlobal);
    tree->Branch("Muminus_muonIsStandalone", &Muminus_muonIsStandalone);
    tree->Branch("Muminus_muonInnerTrackBest", &Muminus_muonInnerTrackBest);

    // CVH fit-convergence diagnostics.
    // edmval_cons0 = -delta(chi^2) at the last GBL iteration of the
    // icons==0 (no-mass-constraint) pass; estimated
    // distance to the chi^2 minimum
    // niter_cons0 = number of iterations to convergence in that pass
    // dmassconvval, dinvmasssqconvval = bias on the dimuon mass and
    // on 1/mass^2 induced by the GBL convergence at
    // finite EDM (used downstream as a kernel-based
    // correction of the lineshape mean)
    // *_cons0 = same quantities recorded at the icons==0 pass
    // (no-mass-constraint), versus the unsuffixed
    // value from the final pass (icons==1 with mass
    // constraint, or icons==0 if no constraint)
    tree->Branch("edmval_cons0", &edmval_cons0);
    tree->Branch("niter_cons0", &niter_cons0);

    tree->Branch("dmassconvval", &dmassconvval);
    tree->Branch("dinvmasssqconvval", &dinvmasssqconvval);
    tree->Branch("dmassconvval_cons0", &dmassconvval_cons0);
    tree->Branch("dinvmasssqconvval_cons0", &dinvmasssqconvval_cons0);

    // openspec/improve-cvh-refit-convergence §1.5: per-iteration debug
    // dump. Branches present only when debugPerIterDump_=true. Each
    // vector spans both icons=0 and icons=1 phases, concatenated in
    // iteration order. mu_qoverp_iter has 2 entries per iteration
    // (Muplus, then Muminus); other vectors have 1 entry per iteration.
    if (debugPerIterDump_) {
      tree->Branch("chisqval_iter",       &chisqval_iter);
      tree->Branch("edmval_iter",         &edmval_iter);
      tree->Branch("edmvalref_iter",      &edmvalref_iter);
      tree->Branch("deltachisqval_iter",  &deltachisqval_iter);
      tree->Branch("mu_qoverp_iter",      &mu_qoverp_iter);
      tree->Branch("Jpsi_mass_iter",      &Jpsi_mass_iter);
    }

    // L1 trigger decisions: one boolean per configured l1Triggers_
    // path (from the channel cfi). Used downstream to re-weight or
    // categorise events by trigger.
    for (std::size_t itrig = 0; itrig < l1Triggers_.size(); ++itrig) {
      tree->Branch(l1Triggers_[itrig].c_str(), &l1TriggerDecisions_[itrig]);
    }

// tree->Branch("hessv", &hessv);

  }
}


// ------------ method called for each event ------------
void ResidualGlobalCorrectionMakerTwoTrackG4e::produce(edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  // Sync the material-group k values from corparms_ into the model so the
  // propagator's per-step provider applies the current calibration.
  if (globalMaterialModel_) {
    for (unsigned int g = 0; g < matGroupGlobalIdx_.size(); ++g) {
      matModel_->setKValue(g, corparms_[matGroupGlobalIdx_[g]]);
    }
  }

  const bool dogen = fitFromGenParms_;
 
  // CVH_LOCAL_UPDATE=1 -- EXPERIMENT (2026-08-10). Default false reproduces the
  // baseline bit-identically.
  //
  // false = the AN's scheme: between Gauss-Newton iterations only the REFERENCE
  // STATE is updated and each track is re-propagated UNSCATTERED with the
  // nominal energy loss, so material and Jacobians are evaluated on a
  // trajectory the particle is increasingly not on. The error is second order
  // in the scattering deviation (~1/p^2), which matches the observed low-p
  // breakdown (spike falls as ~p^-3.4, threshold-like below 3 GeV) far better
  // than the energy-loss 1/p.
  //
  // true = the standard GBL iteration: the per-layer states are updated from
  // the fitted kinks and carried forward, so the next propagation starts from
  // the DEFORMED trajectory and re-samples the material along it; Hp is
  // recomputed at the updated state and Q/dQMS/dQI are transformed to local
  // coordinates. dx0 then holds only the residual w.r.t. the propagated state,
  // so the kink prior is not double-counted.
  //
  // Josh: "all attempts to relax this constraint in the past introduced other
  // problems ... in principle you need a quasi-continuously deformed
  // propagation to the next layer". This switch is to find out WHICH problems.
  static const bool dolocalupdate = (getenv("CVH_LOCAL_UPDATE") != nullptr);
  
  using namespace edm;

  Handle<reco::TrackCollection> trackOrigH;
  iEvent.getByToken(inputTrackOrig_, trackOrigH);

  // Optional Stage-2 B+ candidate index, parallel to inputCandidates_
  // (or to the fallback legacy track-pair loop). Empty handle => bCandIdx
  // stays at its -1 sentinel and the branch was not added.
  Handle<std::vector<int>> bCandIdxH;
  if (!bCandIdxSrcTag_.label().empty()) {
    iEvent.getByToken(bCandIdxToken_, bCandIdxH);
  }

  Handle<reco::TrackCollection>   dedxSourceTracksH;
  Handle<edm::ValueMap<float>> dedxHarmonic2H;
  Handle<edm::ValueMap<float>> dedxPixelHarmonic2H;
  Handle<edm::ValueMap<float>> dedxAllHarmonic2H;
  if (readDeDx_) {
    iEvent.getByToken(dedxSourceTracksToken_,   dedxSourceTracksH);
    iEvent.getByToken(dedxHarmonic2Token_,      dedxHarmonic2H);
    iEvent.getByToken(dedxPixelHarmonic2Token_, dedxPixelHarmonic2H);
    iEvent.getByToken(dedxAllHarmonic2Token_,   dedxAllHarmonic2H);
  }


  // loop over gen particles

  auto globalGeometry = iSetup.getHandle(globalGeometryEventToken_);
  auto trackerTopology = iSetup.getHandle(trackerTopologyEventToken_);
  auto ttrh = iSetup.getHandle(ttrhToken_);

  // MT: bootstrap this TBB worker thread's G4 environment from the master
  // (world + per-thread navigator + per-thread magnetic field). Idempotent
  // per thread.
  worker_->ensureInitialized(iSetup.getData(cvhMasterToken_).cvhMaster());
  setG4RandomEngineForStream(iEvent.streamID());
  // Lazy-init the per-stream propagator clone. Safe AFTER ensureInitialized
  // has put the world in place on this thread.
  if (!streamPropagator_) {
    auto thePropagator = iSetup.getHandle(g4ePropToken_);
    const Geant4ePropagator *templateProp =
        dynamic_cast<const Geant4ePropagator*>(thePropagator.product());
    if (!templateProp) {
      throw cms::Exception("Configuration")
          << "ESProducer for label 'Geant4ePropagator' did not deliver a Geant4ePropagator";
    }
    streamPropagator_.reset(templateProp->clone());
    // Urban/Moliere step logs feed the per-candidate mass-CF export; only
    // pay the bookkeeping when the resolution machinery is on (same gating
    // as the single-track maker).
    streamPropagator_->setIoniStepLogging(doRes_ && (fillGrads_ || fillGradsFactored_));
  }
  const Geant4ePropagator *g4prop = streamPropagator_.get();
  const MagneticField* field = g4prop->magneticField();
  
// Handle<std::vector<reco::GenParticle>> genPartCollection;
  Handle<edm::View<reco::Candidate>> genPartCollection;
  Handle<math::XYZPointF> genXyz0;
  Handle<GenEventInfoProduct> genEventInfo;
  Handle<std::vector<int>> genPartBarcodes;
  Handle<std::vector<PileupSummaryInfo>> pileupSummary;
  if (doGen_) {
    iEvent.getByToken(GenParticlesToken_, genPartCollection);
    iEvent.getByToken(genXyz0Token_, genXyz0);
    iEvent.getByToken(genEventInfoToken_, genEventInfo);
    iEvent.getByToken(pileupSummaryToken_, pileupSummary);
  }
  
  std::vector<Handle<std::vector<PSimHit>>> simHits(inputSimHits_.size());
  edm::Handle<std::vector<SimTrack>> simTracks;
  if (doSim_) {
    iEvent.getByToken(genParticlesBarcodeToken_, genPartBarcodes);
    for (unsigned int isimhit = 0; isimhit<inputSimHits_.size(); ++isimhit) {
      iEvent.getByToken(inputSimHits_[isimhit], simHits[isimhit]);
    }
    iEvent.getByToken(inputSimTracks_, simTracks);
  }
  
  Handle<edm::View<reco::Muon> > muons;
  if (doMuons_) {
    iEvent.getByToken(inputMuons_, muons);
  }
  
  Handle<edm::TriggerResults> triggerResults;
  if (doTrigger_)  {
    iEvent.getByToken(inputTriggerResults_, triggerResults);
  }

  Handle<L1GlobalTriggerReadoutRecord> l1ReadoutRecord;
  if (doL1Trigger_) {
    iEvent.getByToken(inputL1ReadoutRecord_, l1ReadoutRecord);
  }

  KFUpdator updator;
  TkClonerImpl const& cloner = static_cast<TkTransientTrackingRecHitBuilder const *>(ttrh.product())->cloner();
  
  auto TTBuilder = iSetup.getHandle(transTrackBuilderToken_);
  KinematicParticleFactoryFromTransientTrack pFactory;

  Handle<reco::BeamSpot> bsH;
  iEvent.getByToken(inputBs_, bsH);
  
  const std::array<double, 2> trackMass = {{daughterMass1_, daughterMass2_}};
  const std::array<double, 2> trackMassErr = {{daughterMass1Err_, daughterMass2Err_}};
  const double massForConstraintHelpers = 0.5 * (daughterMass1_ + daughterMass2_);

  // Sparse GBL design-matrix formulation (ported from the single-track
  // maker ResidualGlobalCorrectionMakerG4e.cc). Replaces the dense
  // hessfull = MatrixXd::Zero(nparmsfull, nparmsfull) (~1.35 GB at 360
  // field modes) with rfull/Ffull/Jfull/Vinvfull (peak ~tens of MB,
  // sparseView'd before the heavy ops). Mathematically identical to the
  // old dense path: dxfull = -(2 FtVinvF)^-1 (2 FtVinvr) collapses to
  // dxfree = -(FtVinvF)^-1 FtVinvr (the factor of 2 cancels).
  VectorXd rfull;
  MatrixXd Ffull;
  MatrixXd Jfull;
  MatrixXd Vinvfull;
  VectorXd dxfree;

  VectorXd dxfull;
  MatrixXd dxdparms;
  VectorXd grad;
  MatrixXd hess;

  SparseMatrix<double> Fsparse;
  SparseMatrix<double> Vinvsparse;
  SparseMatrix<double> VinvF;
  SimplicialLDLT<SparseMatrix<double>> Cinvd;

  MatrixXd covfull;
  Matrix<double, 6, 6> covrefmom;
  // THE VERTEX FUNCTIONAL's influence row `w_v` (ncons) and its variance,
  // carried from the doRes block (where Vinv/F/C are still alive) to the
  // `Jpsi_jacVtx` emitter, which runs after the global-index remap. Members
  // for the same reason `covrefmom` is one: the two live in sibling scopes.
  Eigen::VectorXd wvtxinf;

  // doRes port (per-candidate mass-CF export): the material process-noise
  // derivative blocks dV_b (MS and ionization parts of Q), their row
  // ranges, and the per-block global-parameter labels. Registered fresh
  // each iteration of the unconstrained (icons==0) pass; the converged
  // iteration's content feeds the mass-projected influence export.
  const bool dores = doRes_;
  // The CPE's own uProj is the strip hit-class variable; the handle is taken
  // once per event, as in the single-track maker.
  const StripCPE *stripCPEForExport =
      doRes_ ? dynamic_cast<const StripCPE *>(iSetup.getHandle(stripCPEToken_).product()) : nullptr;

  std::vector<SparseMatrix<double>> dVs;
  std::vector<std::array<unsigned int, 2>> resblockrng;
  std::vector<unsigned int> resglobidx;
  // number of DISTINCT final columns that received a log-det contribution on
  // the exported pass; bounds the rank the variance block adds to `hess`
  // (rank(A+B) <= rank A + rank B), which is what the factored storage needs.
  unsigned int nvarcols = 0;
  // hit class per registered block, -1 for material (see reshitcls)
  std::vector<int> rescls_;
// FullPivLU<MatrixXd> Cinvd;
// ColPivHouseholderQR<MatrixXd> Cinvd;
  
  std::array<MatrixXd, 2> jacarr;
  
  run = iEvent.run();
  lumi = iEvent.luminosityBlock();
  event = iEvent.id().event();

  genweight = 1.;
  if (doGen_) {
    genweight = genEventInfo->weight();
    
    Pileup_nPU = pileupSummary->front().getPU_NumInteractions();
    Pileup_nTrueInt = pileupSummary->front().getTrueNumInteractions();
  }
  
  // trigger bits
  if (doTrigger_) {
    auto const &triggerNames = iEvent.triggerNames(*triggerResults);
    
    if (triggerNames.parameterSetID() != triggerNamesId_) {
      //trigger menu changed, update list of trigger path idxs
      
      triggerIdxs_.clear();
      
      for (auto const &trigger : triggers_) {
        const std::string basename = trigger + "_v";
// std::cout << "basename = " << basename << std::endl;
        std::size_t idx = triggerNames.size();
        for (std::size_t itrig = 0; itrig < triggerNames.size(); ++itrig) {
// auto findres = triggerNames.triggerName(itrig).find(basename);
// std::cout << triggerNames.triggerName(itrig) << " findres = " << findres << std::endl;
          if (triggerNames.triggerName(itrig).find(basename) == 0) {
            idx = itrig;
            break;
          }
        }
        triggerIdxs_.push_back(idx);
      }
      
// for (std::size_t itrig = 0; itrig < triggerIdxs_.size(); ++itrig) {
// std::cout << "itrig = " << itrig << ", idx = " << triggerIdxs_[itrig] << std::endl;
// }
      
      triggerNamesId_ = triggerNames.parameterSetID();
    }
    
    
    
    
    
    // set trigger decision bits
    for (std::size_t itrig = 0; itrig < triggerIdxs_.size(); ++itrig) {
      const std::size_t idx = triggerIdxs_[itrig];
// if (idx < triggerResults->size()) {
// std::cout << "itrig = " << itrig << " idx = " << idx << " accept = " << triggerResults->accept(idx) << std::endl;
// }
      triggerDecisions_[itrig] = idx < triggerResults->size() ? triggerResults->accept(idx) : false;
    }
    
    
// for (unsigned int i = 0; i < triggerNames.size(); ++i) {
// std::cout << i << " " << triggerNames.triggerName(i) << " accept = " << triggerResults->accept(i) << std::endl;
// }
// auto const idx = iEvent.triggerNames(*triggerResults).triggerIndex("HLT_Mu7p5_Track2_Jpsi");
// auto const idx = triggerNames.triggerIndex("HLT_Mu7p5_Track3p5_Jpsi_v4");
// auto const idx2 = triggerNames.triggerIndex("HLT_eawgawe");
// std::cout << "trigger names size = " << triggerNames.size() << std::endl;
// std::cout << "trigger index = " << idx << std::endl;
// std::cout << "trigger index2 = " << idx2 << std::endl;
  
  }

  if (doL1Trigger_) {
    auto l1Menu = iSetup.getHandle(l1MenuToken_);

    auto const &algorithmMap = l1Menu->gtAlgorithmMap();
    std::vector<std::string> missingTriggers;
    for (std::size_t itrig = 0; itrig < l1Triggers_.size(); ++itrig) {
      auto const menuIt = algorithmMap.find(l1Triggers_[itrig]);
      l1TriggerBitNumbers_[itrig] =
          menuIt != algorithmMap.end() ? static_cast<unsigned int>(menuIt->second.algoBitNumber())
                                       : std::numeric_limits<unsigned int>::max();
      if (menuIt == algorithmMap.end()) {
        missingTriggers.push_back(l1Triggers_[itrig]);
      }
    }

    std::fill(l1TriggerDecisions_.begin(), l1TriggerDecisions_.end(), 0);
    bool validDecisionWord = false;
    bool finalOr = false;
    unsigned int nFinalBitsSet = 0;
    unsigned int decisionWordSize = 0;
    if (l1ReadoutRecord.isValid()) {
      auto const &decisionWord = l1ReadoutRecord->decisionWord();
      validDecisionWord = !decisionWord.empty();
      finalOr = l1ReadoutRecord->decision();
      decisionWordSize = decisionWord.size();
      for (std::size_t ibit = 0; ibit < decisionWord.size(); ++ibit) {
        if (decisionWord.at(ibit)) {
          ++nFinalBitsSet;
        }
      }
      for (std::size_t itrig = 0; itrig < l1TriggerBitNumbers_.size(); ++itrig) {
        const unsigned int bit = l1TriggerBitNumbers_[itrig];
        l1TriggerDecisions_[itrig] =
            bit < decisionWord.size() ? static_cast<int>(decisionWord.at(bit)) : 0;
      }
    }

    std::ostringstream l1DecisionStream;
    l1DecisionStream << "run:lumi:event = " << run << ":" << lumi << ":" << event
                     << " L1 source=" << inputL1ReadoutRecordTag_.encode()
                     << " handleValid=" << (l1ReadoutRecord.isValid() ? 1 : 0)
                     << " decisionWordValid=" << (validDecisionWord ? 1 : 0)
                     << " decisionWordSize=" << decisionWordSize
                     << " finalOr=" << (finalOr ? 1 : 0)
                     << " nFinalBitsSet=" << nFinalBitsSet;
    if (!missingTriggers.empty()) {
      l1DecisionStream << " missingMenuNames=";
      for (std::size_t itrig = 0; itrig < missingTriggers.size(); ++itrig) {
        if (itrig != 0) {
          l1DecisionStream << ",";
        }
        l1DecisionStream << missingTriggers[itrig];
      }
    }
    l1DecisionStream << " L1 bits:";
    for (std::size_t itrig = 0; itrig < l1Triggers_.size(); ++itrig) {
      l1DecisionStream << " " << l1Triggers_[itrig] << "[";
      if (l1TriggerBitNumbers_[itrig] == std::numeric_limits<unsigned int>::max()) {
        l1DecisionStream << "missing";
      } else {
        l1DecisionStream << l1TriggerBitNumbers_[itrig];
      }
      l1DecisionStream << "]=" << l1TriggerDecisions_[itrig];
    }
    std::cout << l1DecisionStream.str() << std::endl;
  }
  
  // Build the list of (i, j) Track* pairs the CVH refit will iterate.
  // Fast path: when srcCandidates is configured, take each persisted
  // VertexCompositeCandidate's two RecoChargedCandidate daughters (stage-1
  // selection is trusted; no re-pairing or vertex pre-fit). Fallback: the
  // legacy j>i outer-product over the input TrackCollection.
  std::vector<std::array<const reco::Track*, 2>> trackPairs;
  // Stage-2 per-pair index, accumulated alongside trackPairs so the
  // per-pair loop below can fill bCandIdx positionally. Empty when no
  // bCandIdxSrc was configured -- per-row Fill() sees -1.
  std::vector<int> bCandIdxPerPair;
  // Index into the srcCandidates collection for each pair (candidate-driven
  // mode), so per-candidate ValueMap outputs can be written back positionally
  // (one entry per input candidate, sentinel for skipped/failed). -1 in the
  // legacy fallback (no candidate collection to key to).
  std::vector<int> candCollIdxPerPair;
  // The srcCandidates handle, kept at produce() scope so the ValueMaps can be
  // sized and keyed to it after the pair loop.
  Handle<reco::VertexCompositeCandidateCollection> vmCandH;
  if (!inputCandidatesTag_.label().empty()) {
    Handle<reco::VertexCompositeCandidateCollection> candH;
    iEvent.getByToken(inputCandidates_, candH);
    vmCandH = candH;
    trackPairs.reserve(candH->size());
    bCandIdxPerPair.reserve(candH->size());
    candCollIdxPerPair.reserve(candH->size());
    for (std::size_t ic = 0; ic < candH->size(); ++ic) {
      const auto& cand = (*candH)[ic];
      if (cand.numberOfDaughters() < 2) continue;

      // Resolve this instance's two-track subsystem from the candidate tree.
      std::array<const reco::Track*, 2> pair{{nullptr, nullptr}};
      bool ok = false;
      if (subsystemDaughter_ >= 0) {
        // Explicit selection: descend into the named composite daughter. Used
        // to give a second instance the X subsystem of a J/psi + X candidate.
        if (subsystemDaughter_ < static_cast<int>(cand.numberOfDaughters()))
          ok = subsystemPair(cand.daughter(subsystemDaughter_), pair);
      } else {
        // Auto: a flat candidate (both daughters leaves) is fit directly --
        // this is the legacy path and stays bit-identical. Otherwise the
        // candidate is nested (e.g. B+ -> [J/psi] K+), so descend into
        // daughter(0), which the stage-1 layout fixes as the composite.
        ok = subsystemPair(&cand, pair);
        if (!ok) ok = subsystemPair(cand.daughter(0), pair);
      }
      if (!ok) continue;

      trackPairs.push_back(pair);
      bCandIdxPerPair.push_back(
          (bCandIdxH.isValid() && ic < bCandIdxH->size()) ? (*bCandIdxH)[ic] : -1);
      candCollIdxPerPair.push_back(static_cast<int>(ic));
    }
  } else {
    if (trackOrigH->size() >= 2) {
      trackPairs.reserve(trackOrigH->size() * (trackOrigH->size() - 1) / 2);
    }
    for (auto itrack = trackOrigH->begin(); itrack != trackOrigH->end(); ++itrack)
      for (auto jtrack = itrack + 1; jtrack != trackOrigH->end(); ++jtrack) {
        trackPairs.push_back({{&*itrack, &*jtrack}});
        bCandIdxPerPair.push_back(-1);  // no candidate-level index in the legacy outer-product
        candCollIdxPerPair.push_back(-1);
      }
  }

  // Per-candidate ValueMap accumulators (sentinel-initialised to one entry per
  // input candidate; filled at the point the per-candidate tree row is written).
  const bool doVM = produceValueMaps_ && vmCandH.isValid();
  const std::size_t nVMCand = doVM ? vmCandH->size() : 0;
  std::vector<float> vmCorMassV(nVMCand, -99.f), vmCorMassErrV(nVMCand, -99.f),
      vmCorPtV(nVMCand, -99.f), vmCorEtaV(nVMCand, -99.f), vmCorPhiV(nVMCand, -99.f),
      vmMuPlusPtV(nVMCand, -99.f), vmMuPlusEtaV(nVMCand, -99.f), vmMuPlusPhiV(nVMCand, -99.f),
      vmMuMinusPtV(nVMCand, -99.f), vmMuMinusEtaV(nVMCand, -99.f), vmMuMinusPhiV(nVMCand, -99.f),
      vmEdmvalV(nVMCand, -99.f);
  std::vector<std::vector<int>> vmGlobalIdxsV(nVMCand);
  std::vector<std::vector<float>> vmJacRefMuPlusV(nVMCand), vmJacRefMuMinusV(nVMCand),
      vmJacMassV(nVMCand), vmHessFactorV(nVMCand);

  for (std::size_t ipair = 0; ipair < trackPairs.size(); ++ipair) {
    auto& trackPair = trackPairs[ipair];
    bCandIdx = bCandIdxPerPair[ipair];
    const int candCollIdx = candCollIdxPerPair[ipair];
    const reco::Track* itrack = trackPair[0];
    const reco::Track* jtrack = trackPair[1];
    if (itrack->isLooper() || jtrack->isLooper()) {
      continue;
    }

    // All two-track channels fit a neutral parent; a same-sign pair can
    // never satisfy the charge-sum requirement enforced after the update,
    // so it would waste a kinematic fit plus a full GN iteration and then
    // abort deterministically. Skip it up front (counted separately -- these
    // are not fit failures).
    if (itrack->charge() + jtrack->charge() != 0) {
      ++fitSkippedSameSign_;
      continue;
    }
    
    const reco::Candidate *mu0gen = nullptr;
    double drmin0 = 0.1;
    if (doGen_ && (!doSim_ || fitSimHitPositions_)) {
      for (auto const &genpart : *genPartCollection) {
        if (genpart.status() != 1) {
          continue;
        }
        if (std::abs(genpart.pdgId()) != 13) {
          continue;
        }
        
        const double dR0 = deltaR(genpart, *itrack);
        if (dR0 < drmin0 && genpart.charge() == itrack->charge()) {
          mu0gen = &genpart;
          drmin0 = dR0;
        }
      }
    }
    
    if (requireGen_ && mu0gen == nullptr) {
      continue;
    }

    const reco::TransientTrack itt = TTBuilder->build(*itrack);


    const reco::Muon *matchedmuon0 = nullptr;
    if (doMuons_) {
      for (auto const &muon : *muons) {
        if (muon.bestTrack()->algo() == itrack->algo()) {
          if ( (muon.bestTrack()->momentum() - itrack->momentum()).mag2() < 1e-3 ) {
            matchedmuon0 = &muon;
          }
        }
        else if (muon.innerTrack().isNonnull() && muon.innerTrack()->algo() == itrack->algo()) {
          if ( (muon.innerTrack()->momentum() - itrack->momentum()).mag2() < 1e-3 ) {
            matchedmuon0 = &muon;
          }
        }
      }
    }

    {
      std::array<ROOT::Math::PxPyPzMVector, 2> mutrkarr;
      mutrkarr[0] = ROOT::Math::PxPyPzMVector(itrack->px(), itrack->py(), itrack->pz(), trackMass[0]);
      mutrkarr[1] = ROOT::Math::PxPyPzMVector(jtrack->px(), jtrack->py(), jtrack->pz(), trackMass[1]);

      const reco::Candidate *mu1gen = nullptr;
      double drmin1 = 0.1;
      
      double massconstraintval = massConstraint_;
      if (doGen_ && (!doSim_ || fitSimHitPositions_)) {
        for (auto const &genpart : *genPartCollection) {
          if (genpart.status() != 1) {
            continue;
          }
          if (std::abs(genpart.pdgId()) != 13) {
            continue;
          }
          
          const double dR1 = deltaR(genpart, *jtrack);
          if (dR1 < drmin1 && genpart.charge() == jtrack->charge()) {
            mu1gen = &genpart;
            drmin1 = dR1;
          }
        }
        
      }
      
      if (requireGen_ && mu1gen == nullptr) {
        continue;
      }
      
      
      const reco::TransientTrack jtt = TTBuilder->build(*jtrack);

      const reco::Muon *matchedmuon1 = nullptr;
      if (doMuons_) {
        for (auto const &muon : *muons) {
          if (muon.bestTrack()->algo() == jtrack->algo()) {
            if ( (muon.bestTrack()->momentum() - jtrack->momentum()).mag2() < 1e-3 ) {
              matchedmuon1 = &muon;
            }
          }
          else if (muon.innerTrack().isNonnull() && muon.innerTrack()->algo() == jtrack->algo()) {
            if ( (muon.innerTrack()->momentum() - jtrack->momentum()).mag2() < 1e-3 ) {
              matchedmuon1 = &muon;
            }
          }
        }
      }
      
// std::cout << "massconstraintval = " << massconstraintval << std::endl;
    
      std::array<TransientTrackingRecHit::RecHitContainer, 2> hitsarr;

      // Per-muon pathology-class counts of pixel hits admitted to the fit
      // (16 combination bins, see Mu{plus,minus}_pixClass) + demoted count.
      std::array<std::array<int, 16>, 2> pixclassarr = {{{{0}}, {{0}}}};
      std::array<unsigned int, 2> npixdemotedarr = {{0u, 0u}};

      // prepare hits
      for (unsigned int id = 0; id < 2; ++id) {
        const reco::Track &track = id == 0 ? *itrack : *jtrack;
        auto &hits = hitsarr[id];
        hits.reserve(track.recHitsSize());
      
      
        for (auto it = track.recHitsBegin(); it != track.recHitsEnd(); ++it) {
          if ((*it)->geographicalId().det() != DetId::Tracker) {
            continue;
          }

          // Leg-structure-free mode: drop hitless surfaces (see the
          // single-track maker for the rationale).
          if (skipHitlessSurfaces_ && !(*it)->isValid()) {
            continue;
          }

          // hits on garbage-shifted modules: dropped (drop policy) or
          // re-inserted at the repaired-surface path position (reorder
          // policy); see the single-track producer for the rationale
          if (!garbageShiftReorderHits_ && garbageShiftModules_.count((*it)->geographicalId().rawId())) {
            continue;
          }

          const GeomDet* detectorG = globalGeometry->idToDet((*it)->geographicalId());
          const GluedGeomDet* detglued = dynamic_cast<const GluedGeomDet*>(detectorG);
          
          // split matched invalid hits
          if (detglued != nullptr && !(*it)->isValid()) {
// bool order = detglued->stereoDet()->surface().position().mag() > detglued->monoDet()->surface().position().mag();
            
            const auto stereopos = detglued->stereoDet()->surface().position();
            const auto monopos = detglued->monoDet()->surface().position();
            
            const Eigen::Vector3d stereoposv(stereopos.x(), stereopos.y(), stereopos.z());
            const Eigen::Vector3d monoposv(monopos.x(), monopos.y(), monopos.z());
            const Eigen::Vector3d trackmomv(track.momentum().x(), track.momentum().y(), track.momentum().z());
            
            bool order = (stereoposv - monoposv).dot(trackmomv) > 0.;
            
            const GeomDetUnit* detinner = order ? detglued->monoDet() : detglued->stereoDet();
            const GeomDetUnit* detouter = order ? detglued->stereoDet() : detglued->monoDet();
            
            if (garbageShiftReorderHits_ || !garbageShiftModules_.count(detinner->geographicalId().rawId())) {
              hits.push_back(TrackingRecHit::RecHitPointer(new InvalidTrackingRecHit(*detinner, (*it)->type())));
            }
            if (garbageShiftReorderHits_ || !garbageShiftModules_.count(detouter->geographicalId().rawId())) {
              hits.push_back(TrackingRecHit::RecHitPointer(new InvalidTrackingRecHit(*detouter, (*it)->type())));
            }
          }
          else {
            // apply hit quality criteria
            const bool ispixel = GeomDetEnumerators::isTrackerPixel(detectorG->subDetector());
            bool hitquality = false;
            if (applyHitQuality_ && (*it)->isValid()) {
              const TrackerSingleRecHit* tkhit = dynamic_cast<const TrackerSingleRecHit*>(*it);
              assert(tkhit != nullptr);
              
              if (ispixel) {
                const SiPixelRecHit *pixhit = dynamic_cast<const SiPixelRecHit*>(tkhit);
                const SiPixelCluster& cluster = *tkhit->cluster_pixel();
                assert(pixhit != nullptr);
                
// std::cout << "getSplitClusterErrorX = " << cluster.getSplitClusterErrorX() << std::endl;
                
// const double jpsieta = dimu_kinfit->currentState().freeTrajectoryState().momentum().eta();
// const double jpsipt = dimu_kinfit->currentState().freeTrajectoryState().momentum().perp();
// if (std::abs(jpsieta)>2. && jpsipt>20. && it == track.recHitsBegin()) {
// std::cout << "id = " << id << " detid = " << (*it)->geographicalId().rawId() << " minPixelRow = " << cluster.minPixelRow() << " maxPixelRow = " << cluster.maxPixelRow() << " minPixelCol = " << cluster.minPixelCol() << " maxPixelCol = " << cluster.maxPixelCol() << std::endl;
// }
                
                ++pixHitsSeen_;
                const bool onEdge = pixhit->isOnEdge();
                if (onEdge) ++pixHitsEdge_;
                if (cluster.sizeX() <= 1) ++pixHitsSizeX1_;
                // Direction-resolved pathology classes: the CPE edge flag
                // does not say WHICH boundary is touched, but the bias is
                // along the truncated coordinate, so resolve it from the
                // cluster extent vs the sensor edge rows/columns.
                const PixelTopology* pixtopo =
                    dynamic_cast<const PixelTopology*>(&detectorG->topology());
                assert(pixtopo != nullptr);
                const bool edgeX = pixtopo->isItEdgePixelInX(cluster.minPixelRow()) ||
                                   pixtopo->isItEdgePixelInX(cluster.maxPixelRow());
                const bool edgeY = pixtopo->isItEdgePixelInY(cluster.minPixelCol()) ||
                                   pixtopo->isItEdgePixelInY(cluster.maxPixelCol());
                const unsigned int icls = (edgeX ? 1u : 0u) |
                                          (edgeY ? 2u : 0u) |
                                          (cluster.sizeX() <= 1 ? 4u : 0u) |
                                          (cluster.sizeY() <= 1 ? 8u : 0u);
                const unsigned int isub =
                    GeomDetEnumerators::isBarrel(detectorG->subDetector()) ? 0 : 1;
                ++pixHitsClass_[isub][icls];
                // Boundary veto configurable via keepPixelEdgeHits; sizeX
                // threshold configurable via pixelMinSizeX (default 2 = legacy).
                hitquality = (keepPixelEdgeHits_ || !onEdge) && cluster.sizeX() >= pixelMinSizeX_
                          && cluster.sizeY() >= pixelMinSizeY_;
                if (!hitquality) ++pixHitsDemoted_;
                if (hitquality) ++pixclassarr[id][icls];
                else ++npixdemotedarr[id];
// hitquality = !pixhit->isOnEdge() && cluster.sizeX() > 1 && pixhit->qBin() < 2;
// hitquality = !pixhit->isOnEdge() && cluster.sizeX() > 1 && cluster.sizeY() > 1;
              }
              else {
                assert(tkhit->cluster_strip().isNonnull());
                const SiStripCluster& cluster = *tkhit->cluster_strip();
                const StripTopology* striptopology = dynamic_cast<const StripTopology*>(&(detectorG->topology()));
                assert(striptopology);
                
                const uint16_t firstStrip = cluster.firstStrip();
                const uint16_t lastStrip = cluster.firstStrip() + cluster.amplitudes().size() - 1;
                const bool isOnEdge = firstStrip == 0 || lastStrip == (striptopology->nstrips() - 1);
                
    // if (isOnEdge) {
    // std::cout << "strip hit isOnEdge" << std::endl;
    // }
                
// hitquality = !isOnEdge;
                hitquality = true;
              }
              
            }
            else {
              hitquality = true;
            }

            
            if (hitquality) {
              hits.push_back((*it)->cloneForFit(*detectorG));
            }
            else if (!skipHitlessSurfaces_) {
              hits.push_back(TrackingRecHit::RecHitPointer(new InvalidTrackingRecHit(*detectorG, TrackingRecHit::inactive)));
            }
          }
        }
        if (garbageShiftReorderHits_ && !garbageShiftModules_.empty()) {
          reorderGarbageShiftHits(hits, track.momentum());
        }
      }

      unsigned int nhits = 0;
      unsigned int nvalid = 0;
      unsigned int nvalidpixel = 0;
      unsigned int nvalidalign2d = 0;
      // Extra alignment-block columns from the pixel pathological-hit
      // class corrections (parmtypes 16-21): 2 per edge-x hit (mean+diff),
      // 2 per edge-y, 1 per sizeX1, 1 per sizeY1. Counted here so the
      // fillAlignGrads appends match nparsAlignment exactly.
      unsigned int nparsPixClass = 0;
      
      
      std::array<unsigned int, 2> nhitsarr = {{ 0, 0 }};
      std::array<unsigned int, 2> nvalidarr = {{ 0, 0 }};
      std::array<unsigned int, 2> nvalidpixelarr = {{ 0, 0 }};
      // Per-track Final counters parallel to nvalidarr / nvalidpixelarr,
      // incremented in the per-hit loop only when `morehitquality` passes
      // (currently always true). Filled to Mu{plus,minus}_nvalidFinal at
      // the per-event branch-write block. Symmetric with the existing
      // nvalid arrays so the per-muon final-hit counts are observable.
      std::array<unsigned int, 2> nvalidFinalarr = {{ 0, 0 }};
      std::array<unsigned int, 2> nvalidpixelFinalarr = {{ 0, 0 }};
      std::array<unsigned int, 2> nmatchedvalidarr = {{ 0, 0 }};
      std::array<unsigned int, 2> nambiguousmatchedvalidarr = {{ 0, 0 }};

      // transmission probe: per-leg sum of the propagator's applied mean
      // energy loss (GeV). Reset per leg at the top of each iteration's hit
      // loop, so after the loop it holds the last (converged) iteration.
      std::array<double, 2> dErefarr = {{ 0., 0. }};
      std::array<double, 2> maxFracLossArr = {{ 0., 0. }};
      
      const std::array<bool, 2> highpurityarr = {{ itrack->quality(reco::TrackBase::highPurity),
                                                  jtrack->quality(reco::TrackBase::highPurity) }};

      // second loop to count hits
      for (unsigned int id = 0; id < 2; ++id) {
        std::unordered_map<unsigned int, unsigned int> trackidmap;
        auto const &hits = hitsarr[id];
// layerStatesarr[id].reserve(hits.size());
        for (auto const &hit : hits) {

          const uint32_t gluedid = trackerTopology->glued(hit->geographicalId());
          const bool isglued = gluedid != 0;
          const DetId parmdetid = isglued ? DetId(gluedid) : hit->geographicalId();

          const DetId aligndetid = alignGlued_ ? parmdetid : hit->geographicalId();

          ++nhits;
          ++nhitsarr[id];
          if (hit->isValid()) {
            ++nvalid;
            ++nvalidarr[id];
            
            const bool ispixel = GeomDetEnumerators::isTrackerPixel(hit->det()->subDetector());
            if (ispixel) {
              ++nvalidpixel;
              ++nvalidpixelarr[id];
              if (pixelHitClassCorrections_) {
                const TrackerSingleRecHit* tkhit = dynamic_cast<const TrackerSingleRecHit*>(&*hit);
                const PixelTopology* pixtopo =
                    dynamic_cast<const PixelTopology*>(&hit->det()->topology());
                if (tkhit != nullptr && pixtopo != nullptr && tkhit->cluster_pixel().isNonnull()) {
                  const SiPixelCluster& cl = *tkhit->cluster_pixel();
                  const bool edgeX = cl.minPixelRow() == 0 ||
                                     cl.maxPixelRow() == pixtopo->nrows() - 1;
                  const bool edgeY = cl.minPixelCol() == 0 ||
                                     cl.maxPixelCol() == pixtopo->ncolumns() - 1;
                  if (pixelLorentzParam_) {
                    // dtanLA column on every valid pixel hit + the
                    // remaining empirical columns (17 edge-x-diff,
                    // 18/19 edge-y, 21 sizeY1); 16/20 replaced.
                    nparsPixClass += 1u + (edgeX ? 1u : 0u) + (edgeY ? 2u : 0u) +
                                     (cl.sizeY() <= 1 ? 1u : 0u);
                  } else {
                    nparsPixClass += (edgeX ? 2u : 0u) + (edgeY ? 2u : 0u) +
                                     (cl.sizeX() <= 1 ? 1u : 0u) + (cl.sizeY() <= 1 ? 1u : 0u);
                  }
                }
              }
            }
            
            
            const bool align2d = detidparms.count(std::make_pair(1, aligndetid));
            if (align2d) {
              ++nvalidalign2d;
            }
            
// std::cout << hit->localPosition() << std::endl;
            
// std::cout << "one hit:" << std:: endl;
            
            // count matching hits from sim tracks
            std::unordered_set<unsigned int> trackidset;
            if (doSim_) {
              for (auto const& simhith : simHits) {
                for (const PSimHit& simHit : *simhith) {
// if (simHit.detUnitId() == hit->geographicalId()) {
// // std::cout << "trackId = " << simHit.trackId() << " particleType = " << simHit.particleType() << "localpos = " << simHit.localPosition() << std::endl;
// }
                  
// if (simHit.detUnitId() == hit->geographicalId()) {
                  if (simHit.detUnitId() == hit->geographicalId() && std::abs(simHit.particleType()) == 13) {
                    //only count each simtrack once on a given detid
                    if (trackidset.count(simHit.trackId())) {
                      continue;
                    }
                    else {
                      trackidset.insert(simHit.trackId());
                    }
                    if (trackidmap.count(simHit.trackId())) {
                      ++trackidmap[simHit.trackId()];
                    }
                    else {
                      trackidmap[simHit.trackId()] = 1;
                    }
                  }
                }
              }
              if (trackidset.size() > 1) {
                ++nambiguousmatchedvalidarr[id];
              }
              
            }
            
          } 
        }
        
        // check for 50%+1 match
        for (auto const &pair : trackidmap) {
          const unsigned int trackid = pair.first;
          const unsigned int nmatch = pair.second;
          if (nmatch > nvalidarr[id]/2) {
            nmatchedvalidarr[id] = nmatch;
            //50% + 1 match found, find the sim track
            for (auto const& simTrack : *simTracks) {
              if (simTrack.trackId() == trackid) {
                // now find corresponding gen particle
                for (auto g = genPartCollection->begin(); g != genPartCollection->end(); ++g) {
                  const int genBarcode = (*genPartBarcodes)[g - genPartCollection->begin()];
                  if (genBarcode == simTrack.genpartIndex()) {
                    if (id == 0) {
                      mu0gen = &(*g);
                    }
                    else if (id == 1) {
                      mu1gen = &(*g);
                    }
                    break;
                  }
                }
                break;
              }
            }
            break;
          }
        }
        
      }

      if (nhitsarr[0] == 0 || nhitsarr[1] == 0) {
        continue;
      }

      // THE MINIMUM-SIZE REQUIREMENT (see the minNdof_ / minPairHits_ docs).
      // Evaluated for the unconstrained-mass pass, which is the one whose
      // vertex residual and chi2/ndof are exported; the mass-constrained pass
      // has one degree of freedom more and can never be the binding one.
      {
        const long long ndofpre = (long long)nvalid + (long long)nvalidpixel - 10LL
                                  + (bsConstraint_ ? 3LL : 0LL)
                                  + (doPointingConstraint_ ? 1LL : 0LL)
                                  + (doVtxConstraint_ ? 1LL : 0LL);
        if (minNdof_ > 0 && ndofpre < (long long)minNdof_) {
          ++fitSkippedNdof_;
          continue;
        }
        if (minPairHits_ > 0 && (int)nvalid < minPairHits_) {
          ++fitSkippedHits_;
          continue;
        }
        if (minLegHits_ > 0 && ((int)nvalidarr[0] < minLegHits_ ||
                                (int)nvalidarr[1] < minLegHits_)) {
          ++fitSkippedLegHits_;
          continue;
        }
      }
      
// if (mu0gen == nullptr || mu1gen == nullptr || mu0gen->eta()<2.2 || mu1gen->eta()<2.2) {
// continue;
// }
      
      
      AlgebraicSymMatrix55 null55;
      const CurvilinearTrajectoryError nullerr(null55);

      
// const unsigned int nparsAlignment = 2*nvalid + nvalidalign2d;
// const unsigned int nparsAlignment = 6*nvalid;
      const unsigned int nparsAlignment = 5*nvalid + nvalidalign2d + nparsPixClass;
      const unsigned int nFieldModes = fieldCorrection_->nModes();
      const unsigned int nparsBfield = nhits * nFieldModes;
      // Global material model: one slot per group per hit (uncrossed groups
      // contribute zero columns; shared global indices collapse like the
      // field-mode block). Legacy: one per-module eloss slot per hit.
      const unsigned int nMatGroups = globalMaterialModel_ ? matModel_->nGroups() : 0;
      const unsigned int nparsEloss = globalMaterialModel_ ? nhits * nMatGroups : nhits;
      const unsigned int npars = nparsAlignment + nparsBfield + nparsEloss;
      
      const unsigned int nstateparms = 10 + 5*nhits;
      const unsigned int nparmsfull = nstateparms + npars;

      // Sparse GBL: the state params that are actually solved for.
      // fitFromGenParms_ freezes the 10-dim vertex PCA (idx 0..9);
      // doVtxConstraint_ freezes the track-PCA distance (idx 6). A frozen
      // index is simply excluded from freestateidxs instead of being
      // deweighted with a 1e6 diagonal.
      using VectorXb = Matrix<bool, Dynamic, 1>;
      VectorXb freestatemask = VectorXb::Ones(nstateparms);
      if (fitFromGenParms_) {
        freestatemask.head<10>() = Matrix<bool, 10, 1>::Zero();
      }
      if (fitFromSimParms_) {
        freestatemask = VectorXb::Zero(nstateparms);
      }
      if (doVtxConstraint_) {
        freestatemask[6] = false;
      }
      std::vector<Eigen::Index> freestateidxs;
      freestateidxs.reserve(nstateparms);
      for (Eigen::Index istate = 0; istate < (Eigen::Index)nstateparms; ++istate) {
        if (freestatemask[istate]) {
          freestateidxs.push_back(istate);
        }
      }
      const unsigned int nstatefree = freestateidxs.size();

      bool valid = true;
      bool stepClampedThisFit = false;
      // Leg-failure recovery budgets (per candidate): step halvings for
      // failures at iiter > 0, seed-momentum inflations at iiter == 0.
      unsigned int nBacktracks = 0;
      unsigned int nSeedInflations = 0;
      // ---- chi2-based (Armijo) backtracking bookkeeping, per candidate ----
      // chisq0valPrev is the total chi2 at the previous iteration's
      // linearization point; predDecrPrev is the quadratic model's predicted
      // chi2 change of the step actually applied there. stepScaleApplied
      // accumulates the clamp scale within an iteration. nChi2Bt is the
      // halving count since the last accepted step, nChi2BtTotal a per-
      // candidate budget guaranteeing termination (a backtrack redoes the same
      // iteration index and so does not consume the niters budget).
      bool stepBtChi2ThisFit = false;
      unsigned int nChi2Bt = 0;
      unsigned int nChi2BtTotal = 0;
      double chisq0valPrev = std::numeric_limits<double>::quiet_NaN();
      double predDecrPrev = 0.;
      double stepScaleApplied = 1.;
      ++fitAttempted_;
      
      
      if (false) {
        const GlobalPoint fieldrefpoint(itrack->vertex().x(), itrack->vertex().y(), itrack->vertex().z());
        auto const fieldvalref = field->inTesla(fieldrefpoint);
        std::cout << "refpos: " << fieldrefpoint << " bfield = " << fieldvalref << std::endl;
      }

      const unsigned int nicons = doMassConstraint_ ? 2 : 1;
// const unsigned int nicons = doMassConstraint_ ? 3 : 1;

      // openspec/improve-cvh-refit-convergence §1.5: clear per-iter
      // debug vectors at the start of this candidate. push_back happens
      // inside the iter loop below; the vectors span both icons phases
      // concatenated in order.
      if (debugPerIterDump_) {
        chisqval_iter.clear();
        edmval_iter.clear();
        edmvalref_iter.clear();
        deltachisqval_iter.clear();
        mu_qoverp_iter.clear();
        Jpsi_mass_iter.clear();
      }


      // The *cons_* kinematics are written only by the icons != 0 (mass-
      // constrained) pass. With doMassConstraint_ off, nicons == 1 and that
      // pass never runs, so without this reset the branches would carry
      // whatever was in the member from the previous candidate -- or, on the
      // first candidate, uninitialised memory. Reset per candidate (one
      // tree->Fill() per candidate, after the icons loop) so the value is
      // always the -99 "not filled" sentinel used elsewhere in this maker.
      Jpsicons_d = -99.f;
      Jpsicons_x = -99.f;
      Jpsicons_y = -99.f;
      Jpsicons_z = -99.f;
      Jpsicons_pt = -99.f;
      Jpsicons_eta = -99.f;
      Jpsicons_phi = -99.f;
      Jpsicons_mass = -99.f;
      Jpsi_mass_unc = -99.f;
      Jpsi_covmassvtx = 0.f;
      Mupluscons_pt = -99.f;
      Mupluscons_eta = -99.f;
      Mupluscons_phi = -99.f;
      Muminuscons_pt = -99.f;
      Muminuscons_eta = -99.f;
      Muminuscons_phi = -99.f;
      Jpsikincons_x = -99.f;
      Jpsikincons_y = -99.f;
      Jpsikincons_z = -99.f;
      Jpsikincons_pt = -99.f;
      Jpsikincons_eta = -99.f;
      Jpsikincons_phi = -99.f;
      Jpsikincons_mass = -99.f;
      Mupluskincons_pt = -99.f;
      Mupluskincons_eta = -99.f;
      Mupluskincons_phi = -99.f;
      Muminuskincons_pt = -99.f;
      Muminuskincons_eta = -99.f;
      Muminuskincons_phi = -99.f;

      for (unsigned int icons = 0; icons < nicons; ++icons) {

        // Sparse GBL constraint-row count for this icons pass:
        //  - 5 propagation/MS rows per hit (both tracks; nhits is the sum)
        //  - 2 measurement rows per valid hit (the two-track dense maker
        //    uses a uniform 2-dim Fhit/dy0 for every valid hit; strip hits
        //    are handled by a near-singular Vinv on the unmeasured coord,
        //    NOT by emitting fewer rows -- so it is 2*nvalid, not the
        //    single-track maker's nvalid+nvalidpixel)
        //  - 3 beamspot rows ONCE PER PAIR when bsConstraint_ (off by
        //    default). It USED to be 3*2 -- the emission sat inside the
        //    per-leg `for (id...)` loop with no `id == 0` guard, so the same
        //    three rows (same residual, same Jacobian on the same three
        //    state indices, same weight) went in twice and `chisq0val`
        //    counted the beam chi2 twice, HALVING the effective luminous-
        //    region covariance -- while `ndof` below already counted +3.
        //    The luminous region is ONE Gaussian noise block shared by the
        //    two legs, like one extra hit; it enters once.
        //  - 1 pointing row when doPointingConstraint_ (off by default)
        //  - 1 J/psi-mass row on the constrained pass (icons==1)
        const unsigned int nbscons = bsConstraint_ ? 3u : 0u;
        const unsigned int npointcons = doPointingConstraint_ ? 1u : 0u;
        const unsigned int nmasscons = (icons == 1) ? 1u : 0u;
        const unsigned int ncons =
            5u * nhits + 2u * nvalid + nbscons + npointcons + nmasscons;

        // common vertex fit
        std::vector<RefCountedKinematicParticle> parts;
        
        float daughterMass1Err = trackMassErr[0];
        float daughterMass2Err = trackMassErr[1];
        float chisq = 0.;
        float ndf = 0.;
        parts.push_back(pFactory.particle(itt, trackMass[0], chisq, ndf, daughterMass1Err));
        parts.push_back(pFactory.particle(jtt, trackMass[1], chisq, ndf, daughterMass2Err));
        
        RefCountedKinematicTree kinTree;
        if (icons > 0) {
// double kinconstraintval = 1./std::sqrt(massconstraintval);
          double kinconstraintval = massconstraintval;
          TwoTrackMassKinematicConstraint constraint(kinconstraintval);
          KinematicConstrainedVertexFitter vtxFitter;
          kinTree = vtxFitter.fit(parts, &constraint);
        }
        else {
          KinematicParticleVertexFitter vtxFitter;
          kinTree = vtxFitter.fit(parts);
        }
        
        if (kinTree->isEmpty() || !kinTree->isConsistent()) {
// continue;
          std::cout << "Abort: invalid kinematic fit!\n";
          ++fitFailKinFit_;
          valid = false;
          break;
        }
        
        kinTree->movePointerToTheTop();
        RefCountedKinematicParticle dimu_kinfit = kinTree->currentParticle();
        const double m0 = dimu_kinfit->currentState().mass();

        if (debugPerIterDump_) {
          RefCountedKinematicVertex dbgvtx = kinTree->currentDecayVertex();
          std::cout << "dbgSeed: icons=" << icons
                    << " kinvtx=" << dbgvtx->position()
                    << " kinmass=" << m0
                    << " seed0(q,pt,eta)=(" << itrack->charge() << "," << itrack->pt() << "," << itrack->eta() << ")"
                    << " seed1(q,pt,eta)=(" << jtrack->charge() << "," << jtrack->pt() << "," << jtrack->eta() << ")"
                    << std::endl;
        }
        
        if (false) {
          // debug output
  // kinTree->movePointerToTheTop();
          
  // RefCountedKinematicParticle dimu_kinfit = kinTree->currentParticle();
          RefCountedKinematicVertex dimu_vertex = kinTree->currentDecayVertex();
          
          std::cout << dimu_kinfit->currentState().mass() << std::endl;
          std::cout << dimu_vertex->position() << std::endl;
        }
        
        const std::vector<RefCountedKinematicParticle> outparts = kinTree->finalStateParticles();
// std::array<Matrix<double, 7, 1>, 2> refftsarr = {{ outparts[0]->currentState().freeTrajectoryState(),
// outparts[1]->currentState().freeTrajectoryState() }};
        // reference FreeTrajectoryState array (refftsarr)
        std::array<Matrix<double, 7, 1>, 2> refftsarr;
        
        if (fitFromGenParms_ && mu0gen != nullptr && mu1gen != nullptr) {
          refftsarr[0][0] = mu0gen->vertex().x();
          refftsarr[0][1] = mu0gen->vertex().y();
          refftsarr[0][2] = mu0gen->vertex().z();
          refftsarr[0][3] = mu0gen->momentum().x();
          refftsarr[0][4] = mu0gen->momentum().y();
          refftsarr[0][5] = mu0gen->momentum().z();
          refftsarr[0][6] = mu0gen->charge();
          
          refftsarr[1][0] = mu1gen->vertex().x();
          refftsarr[1][1] = mu1gen->vertex().y();
          refftsarr[1][2] = mu1gen->vertex().z();
          refftsarr[1][3] = mu1gen->momentum().x();
          refftsarr[1][4] = mu1gen->momentum().y();
          refftsarr[1][5] = mu1gen->momentum().z();
          refftsarr[1][6] = mu1gen->charge();
        }
        else {
          // useStartingState_ == "midPropagated": analytical extrapolation of
          // each daughter's perigee FTS to the midpoint of the two perigees.
          // Falls back to the perigee (Kalman-daughter) path per-event when
          // either extrapolation returns an invalid TSOS.
          bool midPropOk = false;
          if (useStartingState_ == "midPropagated") {
            ++midPropagatedTotalCount_;
            const FreeTrajectoryState fts0 = itt.initialFreeState();
            const FreeTrajectoryState fts1 = jtt.initialFreeState();
            const GlobalPoint p0 = fts0.position();
            const GlobalPoint p1 = fts1.position();
            const GlobalPoint mid(0.5 * (p0.x() + p1.x()),
                                  0.5 * (p0.y() + p1.y()),
                                  0.5 * (p0.z() + p1.z()));
            AnalyticalImpactPointExtrapolator extrap(field);
            const TrajectoryStateOnSurface tsos0 = extrap.extrapolate(fts0, mid);
            const TrajectoryStateOnSurface tsos1 = extrap.extrapolate(fts1, mid);
            if (tsos0.isValid() && tsos1.isValid()) {
              for (unsigned int id = 0; id < 2; ++id) {
                const TrajectoryStateOnSurface& tsos = (id == 0 ? tsos0 : tsos1);
                const GlobalPoint  pos = tsos.globalPosition();
                const GlobalVector mom = tsos.globalMomentum();
                refftsarr[id][0] = pos.x();
                refftsarr[id][1] = pos.y();
                refftsarr[id][2] = pos.z();
                refftsarr[id][3] = mom.x();
                refftsarr[id][4] = mom.y();
                refftsarr[id][5] = mom.z();
                refftsarr[id][6] = static_cast<double>(tsos.charge());
              }
              midPropOk = true;
            } else {
              ++midPropagatedFallbackCount_;
            }
          }
          if (!midPropOk) {
            // perigee (Kalman-daughter) path -- the legacy / default code,
            // also the fallback when midPropagated extrapolation fails.
            for (unsigned int id = 0; id < 2; ++id) {
              const GlobalPoint pos = outparts[id]->currentState().freeTrajectoryState().position();
              const GlobalVector mom = outparts[id]->currentState().freeTrajectoryState().momentum();

              refftsarr[id][0] = pos.x();
              refftsarr[id][1] = pos.y();
              refftsarr[id][2] = pos.z();
              refftsarr[id][3] = mom.x();
              refftsarr[id][4] = mom.y();
              refftsarr[id][5] = mom.z();
              refftsarr[id][6] = outparts[id]->currentState().freeTrajectoryState().charge();
            }
          }
        }
        
        // THE COMMON-VERTEX CONSTRAINT FREEZES INDEX 6 AT ITS SEED VALUE, so
        // the seed must have d = 0 for the constraint to mean "one vertex".
        // The seeds above are each track's own PCA to the midpoint of the two
        // perigees (or the two perigees themselves), whose separation along
        // `n_hat` is the seed DCA -- O(10-100 um) for a J/psi, not zero.
        // Re-expressing the seed with d = 0 moves each reference point by
        // d_seed/2 along `n_hat` and changes nothing else (`twoTrackPca2cart`
        // is the exact inverse of `twoTrackCart2pca` at fixed momenta).
        if (doVtxConstraint_) {
          Matrix<double, 10, 1> statepcaseed = twoTrackCart2pca(refftsarr[0], refftsarr[1]);
          statepcaseed[6] = 0.;
          refftsarr = twoTrackPca2cart(statepcaseed);
        }

        std::array<std::vector<Matrix<double, 7, 1>>, 2> layerStatesarr;
        for (unsigned int id = 0; id < 2; ++id) {
          auto const &hits = hitsarr[id];
          layerStatesarr[id].reserve(hits.size());
        }
        

        double chisqvalold = std::numeric_limits<double>::max();
        // The objective CHANGES between icons passes (icons 1 adds the mass
        // constraint row), so the Armijo history must not cross that boundary.
        chisq0valPrev = std::numeric_limits<double>::quiet_NaN();
        predDecrPrev = 0.;
        nChi2Bt = 0;

        std::array<unsigned int, 2> trackstateidxarr;
        std::array<int, 2> muchargearr;
        
  // constexpr unsigned int niters = 1;
// constexpr unsigned int niters = 3;
// constexpr unsigned int niters = 5;
// constexpr unsigned int niters = 10;

// constexpr unsigned int niters = 1;
        // openspec/improve-cvh-refit-convergence §1.2: cap configurable via
        // nIters_ (cfi default 10 reproduces the baseline).
        const unsigned int niters = (dogen && !dolocalupdate) ? 1 : nIters_;


// const unsigned int niters = icons == 0 ? 10 : 1;
        

        for (unsigned int iiter=0; iiter<niters; ++iiter) {

          // Per-iter Final-counter reset (openspec §1.6). The per-hit loop
          // below runs once per iter; without resetting here, nvalidFinalarr
          // accumulates as nhits × niter. Resetting at the top of each iter
          // means the final value (read after the iter loop) is the LAST
          // iter's pass count -- exactly what "Final" semantically denotes.
          nvalidFinalarr = {{ 0, 0 }};
          nvalidpixelFinalarr = {{ 0, 0 }};

          // Backtracking snapshot: the linearization state at iteration
          // entry, BEFORE the reference update below applies dxfull. On a
          // failed propagation leg the iteration is redone from this state
          // with a halved step (iiter > 0) or an inflated seed momentum for
          // the failing daughter (iiter == 0: end-of-range protons whose
          // modeled dE/dx drains the seed trajectory).
          const std::array<Matrix<double, 7, 1>, 2> refftsarrSnap = refftsarr;
          const std::array<std::vector<Matrix<double, 7, 1>>, 2> layerStatesSnap = layerStatesarr;
          bool retryIter = false;
          int retryFailId = -1;
          stepScaleApplied = 1.;

          // Sparse GBL assembly buffers (replaces dense gradfull/hessfull).
          // Ffull = d(residual)/d(state)  [ncons x nstateparms]
          // Jfull = d(residual)/d(globalparm)  [ncons x npars]
          // Vinvfull = inverse covariance of the residual rows  [ncons x ncons]
          // rfull = residual vector  [ncons]
          rfull = VectorXd::Zero(ncons);
          Ffull = MatrixXd::Zero(ncons, nstateparms);
          Jfull = MatrixXd::Zero(ncons, npars);
          Vinvfull = MatrixXd::Zero(ncons, ncons);

          if (dores) {
            // fresh registration each iteration; the converged iteration's
            // content is what the export reads
            dVs.clear();
            resblockrng.clear();
            resglobidx.clear();
            resfamily_.clear();
            resvalidhit_.clear();
            rescls_.clear();
            // The raw step records are drained on the icons == 0 pass only
            // (see the `icons == 0` guard at the push sites), so they are
            // cleared on that pass only: cleared on the constrained pass too,
            // the tree would carry EMPTY step records whenever
            // doMassConstraint is on.
            if (icons == 0) {
              ioniurbanidx.clear();
              ioniurbanv.clear();
              ioniqscaleidx.clear();
              ioniqscalev.clear();
              radstepidx.clear();
              radstepv.clear();
              radstepspecv.clear();
              msmoliidx.clear();
              msmoliv.clear();
            }
          }

          // Per-hit diagnostics: cleared every iteration so the vectors
          // hold the LAST iteration's residuals when the tree row is
          // written for this candidate.
          if (fillHitDiagnostics_) {
            hitdiag_detid.clear();
            hitdiag_trk.clear();
            hitdiag_charge.clear();
            hitdiag_class.clear();
            hitdiag_dx.clear();
            hitdiag_dy.clear();
            hitdiag_exx.clear();
            hitdiag_eyy.clear();
            hitdiag_lx.clear();
            hitdiag_ly.clear();
          }

          // Running constraint-row cursor (single-track maker calls this
          // `icons`; renamed `irow` here because `icons` is the outer
          // constrained/unconstrained pass index in this two-track maker).
          unsigned int irow = 0;

          // The beam-line block's row offset and the numbers the export
          // needs, refreshed every iteration so the converged iteration's
          // values are what is written. -1 = the rows are off.
          int bsrow = -1;
          Matrix<double, 3, 3> bscovBS = Matrix<double, 3, 3>::Zero();
          Matrix<double, 2, 1> bswidtherr = Matrix<double, 2, 1>::Zero();
          Matrix<double, 3, 1> bsspot = Matrix<double, 3, 1>::Zero();
          Matrix<double, 2, 1> bsslope = Matrix<double, 2, 1>::Zero();
          Matrix<double, 3, 1> bswidth = Matrix<double, 3, 1>::Zero();
          double bschisq0 = 0.;

          globalidxv.clear();
          globalidxv.resize(npars, 0);
          
// nParms = npars;
// if (fillTrackTree_) {
// tree->SetBranchAddress("globalidxv", globalidxv.data());
// }
          
          std::array<Matrix<double, 5, 9>, 2> FdFmrefarr;
// std::array<unsigned int, 2> trackstateidxarr;
          std::array<unsigned int, 2> trackparmidxarr;
          
          unsigned int trackstateidx = 10;
          unsigned int parmidx = 0;
          unsigned int alignmentparmidx = 0;
          
          double chisq0val = 0.;
          
          if (iiter > 0) {
            //update current state from reference point state
            const Matrix<double, 10, 1> statepca = twoTrackCart2pca(refftsarr[0], refftsarr[1]);
            const Matrix<double, 10, 1> statepcaupd = statepca + dxfull.head<10>();

            refftsarr = twoTrackPca2cart(statepcaupd);
          }

          
  // const bool firsthitshared = hitsarr[0][0]->sharesInput(&(*hitsarr[1][0]), TrackingRecHit::some);
          
  // std::cout << "firsthitshared = " << firsthitshared << std::endl;
          
          // Per-track 3D field correction at each track's PCA reference point,
          // from the scalar-potential expansion. Replaces the old per-module
          // dBz lookup at the first hit's parmdetid.
          std::array<Eigen::Vector3d, 2> dBrefarr;
          for (unsigned int id = 0; id < 2; ++id) {
              const GlobalPoint refPos(refftsarr[id][0], refftsarr[id][1], refftsarr[id][2]);
              dBrefarr[id] = fieldCorrection_->getCorrectionAt(refPos, corparms_);
          }

          const Matrix<double, 10, 10> twotrackpca2curvref = twoTrackPca2curvJacobianD(refftsarr[0], refftsarr[1], field, dBrefarr[0], dBrefarr[1]);

          
          for (unsigned int id = 0; id < 2; ++id) {
    // FreeTrajectoryState refFts = outparts[id]->currentState().freeTrajectoryState();
// FreeTrajectoryState &refFts = refftsarr[id];
            
            // First resolution block belonging to THIS leg. The vertex
            // functional's ionization sign carries the leg's CHARGE (the two
            // legs have opposite ones, unlike the mass functional whose sign
            // is -1 for both), so each block has to know which leg it is on.
            resLegStart_[id] = resblockrng.size();

            Matrix<double, 7, 1> &refFts = refftsarr[id];
            auto &hits = hitsarr[id];

            // Build the per-track Geant4 particle name once: the propagator
            // uses the right particle hypothesis (pion / kaon / proton /
            // ...) instead of the muon default. Naming logic shared with
            // the single-track ntuplizer via the base-class helper.
            const std::string& baseName = (id == 0 ? daughterParticleName1_ : daughterParticleName2_);
            const std::string g4PartName = ana_hitanalyzer::g4ParticleName(baseName, static_cast<int>(refftsarr[id][6]));

            std::vector<Matrix<double, 7, 1>> &layerStates = layerStatesarr[id];
                      
            trackstateidxarr[id] = trackstateidx;
            trackparmidxarr[id] = parmidx;
            
            const unsigned int tracknhits = hits.size();

            Matrix<double, 7, 1> updtsos = refFts;

            // transmission probe: restart the assumed-eloss accumulator for
            // this leg on every iteration (and on every backtrack retry,
            // which re-enters here), so the value surviving the loop belongs
            // to the converged reference trajectory.
            dErefarr[id] = 0.;
            maxFracLossArr[id] = 0.;


            // THE BEAM-LINE (LUMINOUS-REGION) CONSTRAINT.
            //
            // ONE Gaussian noise block shared by the two legs -- three rows
            // on the vertex-PCA position -- so it is emitted ONCE PER PAIR,
            // `id == 0`, exactly like the pointing constraint below.  It
            // used to be emitted once per LEG: same residual, same Jacobian,
            // same weight, twice, which halved the effective covariance and
            // double-counted the beam chi2 while `ndof` counted +3.  See the
            // beam-line block in the base-class header for the full account
            // and for the leave-one-out residual the export builds from it.
            if (bsConstraint_ && id == 0) {
              constexpr unsigned int nlocalvtx = 3;
              constexpr unsigned int nlocal = nlocalvtx;
              constexpr unsigned int localvtxidx = 0;
              constexpr unsigned int fullvtxidx = 7;

              // `beamWidthScale_` multiplies the widths (covariance by its
              // square).  At 1e6 the rows are weightless: the gate that the
              // whole block reduces to the rows-OFF fit.
              const double sigb1 = beamWidthScale_*bsH->BeamWidthX();
              const double sigb2 = beamWidthScale_*bsH->BeamWidthY();
              const double sigb3 = beamWidthScale_*bsH->sigmaZ();
              const double dxdz = bsH->dxdz();
              const double dydz = bsH->dydz();
              const double x0 = bsH->x0();
              const double y0 = bsH->y0();
              const double z0 = bsH->z0();
              
              
              // covariance matrix of luminous region in global coordinates
              // taken from https://github.com/cms-sw/cmssw/blob/abc1f17b230effd629c9565fb5d95e527abcb294/RecoVertex/BeamSpotProducer/src/FcnBeamSpotFitPV.cc#L63-L90

              // FIXME xy correlation is not stored and assumed to be zero
              const double corrb12 = 0.;
              
              const double varb1 = sigb1*sigb1;
              const double varb2 = sigb2*sigb2;
              const double varb3 = sigb3*sigb3;
              
              Matrix<double, 3, 3> covBS = Matrix<double, 3, 3>::Zero();
              // parametrisation: rotation (dx/dz, dy/dz); covxy
              covBS(0,0) = varb1;
              covBS(1,0) = covBS(0,1) = corrb12*sigb1*sigb2;
              covBS(1,1) = varb2;
              covBS(2,0) = covBS(0,2) = dxdz*(varb3-varb1)-dydz*covBS(1,0);
              covBS(2,1) = covBS(1,2) = dydz*(varb3-varb2)-dxdz*covBS(1,0);
              covBS(2,2) = varb3;

              // THE COMMON VERTEX, not the leg's reference point.
              // `twoTrackPca2cart` puts the two reference points at
              // `x_v -+ (d/2) n_hat`, so `refftsarr[id].head<3>()` is the
              // vertex only when index 6 is frozen at zero.  The MIDPOINT is
              // `x_v` identically in both regimes, and its Jacobian is
              // exactly `I` on state indices 7,8,9 and ZERO on index 6 and
              // on the momenta -- which is what `Fbs = Identity` asserts.
              const Matrix<double, 3, 1> xv =
                  0.5*(refftsarr[0].head<3>() + refftsarr[1].head<3>());

              Matrix<double, 3, 1> dbs0;
              dbs0[0] = xv[0] - x0;
              dbs0[1] = xv[1] - y0;
              dbs0[2] = xv[2] - z0;

              const Matrix<double, 3, nlocal> Fbs = Matrix<double, 3, 3>::Identity();
              const Matrix<double, 3, 3> covBSinv = covBS.inverse();

              const double bschisq = dbs0.transpose()*covBSinv*dbs0;
              chisq0val += bschisq;

              // Sparse GBL row write: 3 beamspot rows constraining the
              // vertex-PCA position (state idx 7,8,9).  Fbs = Identity(3,3),
              // residual dbs0 = vertex - beamspot, weight covBSinv.
              bsrow = int(irow);
              bscovBS = covBS;
              bsspot << x0, y0, z0;
              bsslope << dxdz, dydz;
              bswidth << sigb1, sigb2, sigb3;
              // the RECORD's own errors on the two transverse widths -- the
              // prior width for the two `beamwidth_*` variance scales.  They
              // are quoted on sigma, so the prior on a VARIANCE scale
              // `k = (sigma'/sigma)^2` is `2 * BeamWidthXError/BeamWidthX`.
              bswidtherr << beamWidthScale_*bsH->BeamWidthXError(),
                            beamWidthScale_*bsH->BeamWidthYError();
              bschisq0 = bschisq;
              rfull.segment<3>(irow) = dbs0;
              Ffull.block(irow, fullvtxidx, 3, nlocalvtx) =
                  Fbs.leftCols<nlocalvtx>();
              Vinvfull.block<3, 3>(irow, irow) = covBSinv;

              // Registered as a RESOLUTION BLOCK, family 16, dV = covBS --
              // the same convention as the hit families (dV_b/dk = V_b for
              // `V_b -> e^k V_b`).  Without it `sum_b |a_b|^2 == sigma^2`
              // would miss the luminous region's own share for EVERY
              // functional.  `kBeamSpotGlobIdx` is a sentinel: family 16 has
              // no `detidparms` entry and adds no column to the quadratic
              // term.
              if (dores) {
                std::vector<Triplet<double>> bscoeffs;
                for (unsigned int a = 0; a < 3; ++a) {
                  for (unsigned int b = 0; b < 3; ++b) {
                    if (covBS(a, b) != 0.) {
                      bscoeffs.emplace_back(irow + a, irow + b, covBS(a, b));
                    }
                  }
                }
                SparseMatrix<double> &dVbs = dVs.emplace_back(ncons, ncons);
                dVbs.setFromTriplets(bscoeffs.begin(), bscoeffs.end());
                resblockrng.push_back({{irow, 3u}});
                resglobidx.push_back(kBeamSpotGlobIdx);
                resfamily_.push_back(kBeamSpotFamily);
                resvalidhit_.push_back(-1);
                rescls_.push_back(-1);
              }

              irow += 3;

            }

            // 2D-transverse pointing-angle constraint on the V0:
            // g(x) = (xv - xBS) * pVy - (yv - yBS) * pVx = 0
            // Soft-Gaussian chi^2 = g^2 / sigma_g^2 with
            // sigma_g = pointingSigma_ * Lxy * |p_xy| (linearised, small angle)
            // Touches vertex-state indices 0-5 (daughter qop/lambda/phi for both tracks)
            // and 7-8 (vertex xy). Index 6 (d0) and 9 (vertex z) are not coupled.
            // Applied once per iteration (guarded with id == 0) since it is intrinsically
            // a 2-track constraint, not per-daughter.
            if (doPointingConstraint_ && id == 0) {
              constexpr unsigned int nlocal = 10;          // full vertex-state block
              constexpr unsigned int fullvtxidx = 0;       // vertex state lives at 0..9

              // Linearisation-point momentum of each daughter (Cartesian, GeV)
              const double p1x = refftsarr[0][3];
              const double p1y = refftsarr[0][4];
              const double p1z = refftsarr[0][5];
              const double p2x = refftsarr[1][3];
              const double p2y = refftsarr[1][4];
              const double p2z = refftsarr[1][5];
              const double p1mag = std::sqrt(p1x*p1x + p1y*p1y + p1z*p1z);
              const double p2mag = std::sqrt(p2x*p2x + p2y*p2y + p2z*p2z);

              // Helix-state coordinates derived from the Cartesian linearisation point
              const double lam1 = (p1mag > 0.) ? std::asin(p1z / p1mag) : 0.;
              const double phi1 = std::atan2(p1y, p1x);
              const double lam2 = (p2mag > 0.) ? std::asin(p2z / p2mag) : 0.;
              const double phi2 = std::atan2(p2y, p2x);
              const double q1 = (refftsarr[0][6] >= 0.) ? +1. : -1.;
              const double q2 = (refftsarr[1][6] >= 0.) ? +1. : -1.;

              // V0 transverse momentum at SV
              const double pVx = p1x + p2x;
              const double pVy = p1y + p2y;
              const double pVxy = std::hypot(pVx, pVy);

              // Beamspot reference and flight vector in xy
              const double xBSp = bsH->x0();
              const double yBSp = bsH->y0();
              const double fx = refftsarr[0][0] - xBSp;
              const double fy = refftsarr[0][1] - yBSp;
              const double Lxy = std::hypot(fx, fy);

              // Constraint value (cm * GeV) and its width sigma_g.
              const double g = fx * pVy - fy * pVx;
              const double sigma_g = std::max(pointingSigma_ * Lxy * pVxy, 1.e-9);
              const double inv_var = 1. / (sigma_g * sigma_g);

              // Build 1x10 Jacobian J_g of g w.r.t. the local vertex-state perturbations.
              // Index layout: 0..2 = (qop,lam,phi)_track0, 3..5 = (qop,lam,phi)_track1,
              // 6 = d0 (no contribution), 7..9 = vertex (x,y,z); z no contribution.
              Matrix<double, 1, nlocal> Jg = Matrix<double, 1, nlocal>::Zero();

              // Per-track block: dpx/d(qop,lam,phi) and dpy/d(qop,lam,phi)
              {
                const double dpx_dqop = -p1x * p1mag * q1;
                const double dpy_dqop = -p1y * p1mag * q1;
                const double dpx_dlam = -p1mag * std::sin(lam1) * std::cos(phi1);
                const double dpy_dlam = -p1mag * std::sin(lam1) * std::sin(phi1);
                const double dpx_dphi = -p1y;
                const double dpy_dphi = +p1x;
                Jg(0, 0) = fx * dpy_dqop - fy * dpx_dqop;
                Jg(0, 1) = fx * dpy_dlam - fy * dpx_dlam;
                Jg(0, 2) = fx * dpy_dphi - fy * dpx_dphi;
              }
              {
                const double dpx_dqop = -p2x * p2mag * q2;
                const double dpy_dqop = -p2y * p2mag * q2;
                const double dpx_dlam = -p2mag * std::sin(lam2) * std::cos(phi2);
                const double dpy_dlam = -p2mag * std::sin(lam2) * std::sin(phi2);
                const double dpx_dphi = -p2y;
                const double dpy_dphi = +p2x;
                Jg(0, 3) = fx * dpy_dqop - fy * dpx_dqop;
                Jg(0, 4) = fx * dpy_dlam - fy * dpx_dlam;
                Jg(0, 5) = fx * dpy_dphi - fy * dpx_dphi;
              }
              Jg(0, 7) = +pVy;          // dg/dxv
              Jg(0, 8) = -pVx;          // dg/dyv
              // Jg(0, 6) and Jg(0, 9) are zero by construction.

              const double pointingChisq = g * g * inv_var;
              chisq0val += pointingChisq;

              // Sparse GBL row write: 1 pointing row. residual = g,
              // design = Jg (1 x 10 over the full vertex-PCA state at
              // fullvtxidx=0), weight = inv_var (scalar).
              rfull(irow) = g;
              Ffull.block(irow, fullvtxidx, 1, nlocal) = Jg;
              Vinvfull(irow, irow) = inv_var;
              irow += 1;
            }



            for (unsigned int ihit = 0; ihit < hits.size(); ++ihit) {
      // std::cout << "ihit " << ihit << std::endl;
              auto const& hit = hits[ihit];
              
              const uint32_t gluedid = trackerTopology->glued(hit->det()->geographicalId());
              const bool isglued = gluedid != 0;
              const DetId parmdetid = isglued ? DetId(gluedid) : hit->geographicalId();
              const GeomDet* parmDet = isglued ? globalGeometry->idToDet(parmdetid) : hit->det();

              const DetId aligndetid = alignGlued_ ? parmdetid : hit->geographicalId();

              const unsigned int elossglobalidx =
                  globalMaterialModel_ ? 0 : detidparms.at(std::make_pair(7, parmdetid));

              // 3D field correction at the propagation start. The per-mode
              // (Bx, By, Bz) basis values feed the chain-rule scaling of the
              // transport-Jacobian dBx/dBy/dBz columns (cols 5,6,7 of the
              // 5x9 transportJacobianBxByBzD).
              const GlobalPoint propStartPos(updtsos[0], updtsos[1], updtsos[2]);
              // Per-step mode: the provider applies the correction inside
              // the propagator; per-leg basis samples not needed.
              const Eigen::Vector3d dB = perStepFieldModes_
                  ? Eigen::Vector3d::Zero()
                  : fieldCorrection_->getCorrectionAt(propStartPos, corparms_);
              std::vector<double> dBxPerMode, dByPerMode, dBzPerMode;
              if (!perStepFieldModes_) {
                fieldCorrection_->getBxBasisAt(propStartPos, dBxPerMode);
                fieldCorrection_->getByBasisAt(propStartPos, dByPerMode);
                fieldCorrection_->getBzBasisAt(propStartPos, dBzPerMode);
              }

              // Global material model: leg-constant dxi is zero; per-step
              // group values are applied by the propagator's provider path
              // (k values synced from corparms_ at the top of produce).
              const double dxival = globalMaterialModel_ ? 0. : corparms_[elossglobalidx];

              const GloballyPositioned<double> &surface = surfacemapD_.at(hit->geographicalId());

              // Save input state so the FD-closure can re-propagate
              // from the same point with a perturbed dB.
              const Eigen::Matrix<double, 7, 1> propInputState = updtsos;

              auto propresult = g4prop->propagateGenericWithJacobianAltD(updtsos, surface, dB, dxival,
                                                                          0., 0., -1., g4PartName,
                                                                          matModel_.get(),
                                                                          matModel_ ? &groupJacs_ : nullptr,
                                                                          (matModel_ && doRes_ && exportMaterialNoise_) ? &groupQs_ : nullptr,
                                                                          fieldModeProvider_.get(),
                                                                          fieldModeProvider_ ? &modeJacs_ : nullptr);
              if (!std::get<0>(propresult)) {
                std::cout << "ResidualGlobalCorrectionMakerTwoTrackG4e ### Abort: Propagation Failed!"
                          << " icons = " << icons << " iiter = " << iiter
                          << " id = " << id << " ihit = " << ihit
                          << " seed0: q=" << itrack->charge() << " pt=" << itrack->pt() << " eta=" << itrack->eta()
                          << " seed1: q=" << jtrack->charge() << " pt=" << jtrack->pt() << " eta=" << jtrack->eta()
                          << std::endl;
                // Leg-failure recovery: redo the iteration with a halved
                // step (iiter > 0: the previous iterate propagated fine, so
                // the too-large update is the culprit) or an inflated seed
                // momentum for the failing daughter (iiter == 0). Only when
                // the retry budget is exhausted is the candidate lost.
                if ((iiter > 0 && nBacktracks < maxBacktracks_) || (iiter == 0 && nSeedInflations < maxSeedInflations_)) {
                  retryIter = true;
                  retryFailId = id;
                } else {
                  ++fitFailProp_;
                  valid = false;
                }
                break;
              }


              updtsos = std::get<1>(propresult);

              // transmission probe: accumulate the mean energy loss this
              // propagation step actually applied. Taken from the propagator
              // input/output states rather than from any dE/dx model call, so
              // it stays correct whatever scales the loss (CVH_DEDX_SCALE,
              // the global material model's k_g, the per-module dxi). The
              // difference is formed BEFORE any local state update below, so
              // only the propagation (not the fit) contributes.
              {
                const double m2 = trackMass[id] * trackMass[id];
                const double pIn = propInputState.segment<3>(3).norm();
                const double eIn = std::sqrt(pIn * pIn + m2);
                const double eOut = std::sqrt(updtsos.segment<3>(3).squaredNorm() + m2);
                dErefarr[id] += eIn - eOut;
                if (pIn > 0.) {
                  maxFracLossArr[id] = std::max(maxFracLossArr[id], (eIn - eOut) / pIn);
                }
              }

              const Matrix<double, 5, 5> Qcurv = std::get<2>(propresult);
              const Matrix<double, 5, 5> dQMScurv = std::get<5>(propresult);
              const Matrix<double, 5, 5> dQIcurv = std::get<6>(propresult);

              // doRes port: material-block global labels (glued detid
              // convention as in the single-track maker) + the propagator
              // step-record drains for the physics-CF export. The logs are
              // cleared at each propagate call, so the drain must happen
              // here, tagged with this block's labels. icons==0 only (the
              // unconstrained pass feeds the mass-CF export).
              unsigned int msglobalidx = 0;
              unsigned int ioniglobalidx = 0;
              if (dores && icons == 0) {
                const uint32_t gluedidprop = trackerTopology->glued(hit->geographicalId());
                const DetId propdetid = gluedidprop ? DetId(gluedidprop) : hit->geographicalId();
                msglobalidx = detidparms.at(std::make_pair(10, propdetid));
                ioniglobalidx = detidparms.at(std::make_pair(11, propdetid));
                for (auto const &ms : g4prop->msStepLog()) {
                  pushMsMoliStep(msglobalidx, ms);
                }
                for (auto const &us : g4prop->ioniStepLog()) {
                  ioniurbanidx.push_back(ioniglobalidx);
                  ioniurbanv.push_back(us.rec.regime);
                  ioniurbanv.push_back(us.rec.gsig2);
                  ioniurbanv.push_back(us.rec.a1);
                  ioniurbanv.push_back(us.rec.e1);
                  ioniurbanv.push_back(us.rec.a2);
                  ioniurbanv.push_back(us.rec.e2);
                  ioniurbanv.push_back(us.rec.a3);
                  ioniurbanv.push_back(us.rec.e0r);
                  ioniurbanv.push_back(us.rec.tmaxr);
                  ioniurbanv.push_back(us.rec.scaling);
                  ioniurbanv.push_back(us.cs);
                  // regime 2/3 (CVH_IONI_EXACTDELTA): two extra columns AFTER
                  // cs, so every existing column index is unchanged. The
                  // stride is `ioniurbanstride` (12, or 14 with the switch on).
                  if (G4UniversalFluctuationForExtrapolator::exactDeltaEnabled()) {
                    ioniurbanv.push_back(us.rec.beta2);
                    ioniurbanv.push_back(us.rec.etot);
                  }
                  // material group of the step, ALWAYS last (see
                  // UrbanIoniStep::stepGroup)
                  ioniurbanv.push_back(us.stepGroup);
                }

                // Ionization-block scale of this leg, [sc, nstep]. THIS MAKER
                // HAS NO CGF OVERRIDE HOOKS (it never calls setCgfOverride),
                // but that alone does NOT make the factor 1: with
                // CgfQoPMode >= 1 the propagator still takes its UNCACHED
                // branch and substitutes there. So the value is exported the
                // same way rather than hard-coded -- it is 1.0 because the
                // driver pins CgfQoPMode=0, and it will be right if that ever
                // changes.
                ioniqscaleidx.push_back(ioniglobalidx);
                ioniqscalev.push_back(g4prop->cgfQScale());
                ioniqscalev.push_back(static_cast<float>(g4prop->ioniStepLog().size()));

                // radiative steps of the same leg (see the single-track maker)
                for (auto const &rs : g4prop->radStepLog()) {
                  radstepidx.push_back(ioniglobalidx);
                  radstepv.push_back(rs.effZ);
                  radstepv.push_back(rs.effA);
                  radstepv.push_back(rs.xg);
                  radstepv.push_back(rs.etotGeV);
                  radstepv.push_back(rs.pGeV);
                  radstepv.push_back(rs.dOverX0);
                  radstepv.push_back(rs.stepCm);
                  radstepv.push_back(rs.dedxRad);
                  radstepv.push_back(rs.dedxBrem);
                  radstepv.push_back(rs.dedxPair);
                  radstepv.push_back(rs.cs);
                  radstepv.push_back(rs.stepGroup);   // column 11: the material group
                  for (int iv = 0; iv < RADSTEP_NV; ++iv) {
                    radstepspecv.push_back(rs.dNdvBrem[iv]);
                  }
                  for (int iv = 0; iv < RADSTEP_NV; ++iv) {
                    radstepspecv.push_back(rs.dNdvPair[iv]);
                  }
                }
              }

              if (debugPerIterDump_) {
                const auto& sp = surface.position();
                std::cout << "dbgHit: icons=" << icons << " iiter=" << iiter
                          << " id=" << id << " ihit=" << ihit
                          << " detid=" << hit->geographicalId().rawId()
                          << " valid=" << hit->isValid()
                          << " surf(r,z)=(" << std::hypot(sp.x(), sp.y()) << "," << sp.z() << ")"
                          << " in(r,z,p)=(" << std::hypot(propInputState[0], propInputState[1])
                          << "," << propInputState[2] << "," << propInputState.segment<3>(3).norm() << ")"
                          << " out(r,z,p)=(" << std::hypot(updtsos[0], updtsos[1])
                          << "," << updtsos[2] << "," << updtsos.segment<3>(3).norm() << ")"
                          << " dEdxlast=" << std::get<4>(propresult)
                          << std::endl;
              }
              const Matrix<double, 5, 9> FdFm = std::get<3>(propresult);
              const double dEdxlast = std::get<4>(propresult);

              const Matrix<double, 5, 5> Hm =
                  curv2localJacobianAltelossD(updtsos, field, surface, dEdxlast, trackMass[id], dB);
              
              const Matrix<double, 6, 1> localparmsprop = globalToLocal(updtsos, surface);

              Matrix<double, 6, 1> localparms = localparmsprop;
              
              Matrix<double, 5, 1> dx0 = Matrix<double, 5, 1>::Zero();
              if (dolocalupdate) {
                if (iiter==0) {
                  layerStates.push_back(updtsos);
                }
                else {
                  //current state from previous state on this layer
                  //save current parameters 
                  
                  Matrix<double, 7, 1>& oldtsos = layerStates[ihit];
                  const Matrix<double, 5, 5> Hold =
                      curv2localJacobianAltelossD(oldtsos, field, surface, dEdxlast, trackMass[id], dB);
                  const Matrix<double, 5, 1> dxlocal = Hold*dxfull.segment<5>(trackstateidx + 5*ihit);

                  localparms = globalToLocal(oldtsos, surface);

                  localparms.head<5>() += dxlocal;

                  oldtsos = localToGlobal(localparms, surface);

                  updtsos = oldtsos;

                  dx0 = (localparms - localparmsprop).head<5>();
                }
              }

              // curvilinear to local jacobian
              // const Matrix<double, 5, 5> &Hp = dolocalupdate ? curv2localJacobianAltelossD(updtsos, field, surface, dEdxlast, mmu, dB) : Hm;
              Matrix<double, 5, 5> Hp = Hm;
              if (dolocalupdate) {
                Hp = curv2localJacobianAltelossD(updtsos, field, surface, dEdxlast, trackMass[id], dB);
              }
              
              // const Matrix<double, 5, 5> &Q = dolocalupdate ? Hm*Qcurv*Hm.transpose() : Qcurv;

              Matrix<double, 5, 5> Q = Qcurv;
              if (dolocalupdate) {
                Q = Hm*Qcurv*Hm.transpose();
              }

              Matrix<double, 5, 5> dQMS = dQMScurv;
              Matrix<double, 5, 5> dQI = dQIcurv;
              if (dolocalupdate) {
                dQMS = Hm*dQMScurv*Hm.transpose();
                dQI = Hm*dQIcurv*Hm.transpose();
              }

              // ---- PARMTYPE-15: THE MATERIAL GROUP'S OWN PROCESS NOISE ----
              //
              // `k_g` scales the step's MEAN loss AND, coherently, its MS
              // covariance and ionization variance (`matStepFact` in the
              // propagator's M1 block).  Only the mean was ever
              // differentiated: the parmtype-15 column of
              // `transportJacobianBxByBzD` is the `dxi` column, whose one
              // non-zero row is `dqopdxi`.  Registering the group's own
              // dV here is the two-track counterpart of the single-track
              // maker's block (`ResidualGlobalCorrectionMakerG4e.cc`
              // "PARMTYPE-15" comment); together with the log-det assembly
              // further down it gives `k_g` its WIDTH term.
              //
              // `sum_g dQ_g == dQMS + dQI` exactly (the propagator sums the
              // same per-step `errMS + errI` into `groupQs_`, transported at
              // the same point and by the same Jacobian), so these blocks
              // are a RE-PARTITION of the parmtype-10/11 noise, not an
              // addition to it -- which is why `resinfcovgrp` is kept apart
              // from `resinfcov` below.
              const auto registerMatGroupNoise = [&](unsigned int r0) {
                if (!(dores && exportMaterialNoise_ && globalMaterialModel_)) {
                  return;
                }
                for (auto const &gq : groupQs_) {
                  if (gq.first < 0 || unsigned(gq.first) >= nMatGroups) {
                    continue;
                  }
                  Matrix<double, 5, 5> dQG = gq.second;
                  if (dolocalupdate) {
                    dQG = Hm * gq.second * Hm.transpose();
                  }
                  if (!(dQG.cwiseAbs().maxCoeff() > 0.)) {
                    continue;
                  }
                  std::vector<Triplet<double>> coeffs;
                  coeffs.reserve(25);
                  for (unsigned int ir = 0; ir < 5; ++ir) {
                    for (unsigned int ic = 0; ic < 5; ++ic) {
                      coeffs.emplace_back(r0 + ir, r0 + ic, dQG(ir, ic));
                    }
                  }
                  SparseMatrix<double> &dV = dVs.emplace_back(ncons, ncons);
                  dV.setFromTriplets(coeffs.begin(), coeffs.end());
                  resblockrng.push_back({{r0, 5}});
                  resglobidx.push_back(matGroupGlobalIdx_[gq.first]);
                  resfamily_.push_back(15);          // global material group
                  resvalidhit_.push_back(-1);
                  rescls_.push_back(-1);
                }
              };

              // Guarded inversion of the process noise: a (near-)zero-length
              // leg -- a displaced V0 vertex sitting on the first-hit layer --
              // has Q ~ 0, and a plain inverse poisons the solve with inf
              // (observed as the fail[nan] class: finite r and F, non-finite
              // Vinv, daughters with 1-3 valid hits at refR ~ 5 cm). Floor
              // the eigenvalues so the leg becomes an extremely stiff, rather
              // than exact, constraint.
              Matrix<double, 5, 5> Qinv;
              {
                const SelfAdjointEigenSolver<Matrix<double, 5, 5>> esq(Q);
                const double lmax = esq.eigenvalues()(4);
                const double lfloor = std::max(1e-10 * std::max(lmax, 0.), 1e-16);
                Matrix<double, 5, 1> linv;
                for (int k = 0; k < 5; ++k) {
                  linv(k) = 1. / std::max(esq.eigenvalues()(k), lfloor);
                }
                Qinv = esq.eigenvectors() * linv.asDiagonal() *
                       esq.eigenvectors().transpose();
              }

              // Build the per-hit field+eloss Jacobian: per-mode columns sum
              // the dBx, dBy, dBz transport-Jacobian columns (cols 5,6,7 of
              // the 5x9 transportJacobianBxByBzD) scaled by each mode's
              // Bx/By/Bz basis values at the propagation start; the last
              // column is d/dxi (FdFm.col(8) unchanged).
              const unsigned int nlocalbfield = nFieldModes;
              // Eloss block: one column per material group (global model) or
              // the single per-module dxi column (legacy).
              const unsigned int nlocaleloss = globalMaterialModel_ ? nMatGroups : 1;
              const unsigned int nlocalparms = nlocalbfield + nlocaleloss;

              Matrix<double, 5, Dynamic> dStateDparams(5, nlocalparms);
              if (perStepFieldModes_) {
                for (unsigned int imode = 0; imode < nlocalbfield; ++imode) {
                  dStateDparams.col(imode) = modeJacs_[imode];
                }
              } else {
                for (unsigned int imode = 0; imode < nlocalbfield; ++imode) {
                  dStateDparams.col(imode) = FdFm.col(5) * dBxPerMode[imode]
                                           + FdFm.col(6) * dByPerMode[imode]
                                           + FdFm.col(7) * dBzPerMode[imode];
                }
              }
              if (globalMaterialModel_) {
                // per-group dxi columns from the propagator (zero for groups
                // this leg did not cross)
                dStateDparams.rightCols(nlocaleloss).setZero();
                for (auto const &gc : groupJacs_) {
                  dStateDparams.col(nlocalbfield + gc.first) = gc.second;
                }
              } else {
                dStateDparams.col(nlocalbfield) = FdFm.col(8);
              }

              // ----- Numerical-FD closure (debug) ---------------
              // FDs only the basis-invariant curvilinear components
              // (qop, lambda, phi) -- those can be derived directly from
              // the global-cartesian 7-vector (px, py, pz, q) without
              // knowing the surface's local frame; (xt, yt) need a
              // surface-dependent transformation we don't reconstruct here.
              //
              // For each mode i we use a per-mode eps_i scaled so that
              // eps_i * ||dB_basis_i|| equals epsilonFDClosure_ (treat as a
              // target dB-perturbation magnitude in Tesla). Modes with
              // basis amplitude below 1e-15 at the test point are skipped
              // (analytic prediction is FP-noise; FD can't resolve).
              // Tests the 10 modes with the largest basis amplitude at this
              // point so the FD signal is well-conditioned across mode
              // counts. Runs once per job.
              if (runFDClosure_ && !didFDClosure_ && !perStepFieldModes_ && nlocalbfield > 0) {
                auto qopLamPhi = [](const Eigen::Matrix<double, 7, 1>& s) {
                  const double px = s(3), py = s(4), pz = s(5), q = s(6);
                  const double pT = std::sqrt(px * px + py * py);
                  const double pmag = std::sqrt(pT * pT + pz * pz);
                  Eigen::Vector3d v;
                  v(0) = q / pmag;             // qop
                  v(1) = std::atan2(pz, pT);   // lambda
                  v(2) = std::atan2(py, px);   // phi
                  return v;
                };
                const Eigen::Vector3d cNom = qopLamPhi(updtsos);
                const double dBtarget = epsilonFDClosure_;  // target dB |Tesla|
                // Rank modes by basis amplitude (sqrt(dBx^2+dBy^2+dBz^2)).
                std::vector<std::pair<double, unsigned int>> sorted;
                sorted.reserve(nlocalbfield);
                for (unsigned int i = 0; i < nlocalbfield; ++i) {
                  const double a = std::sqrt(dBxPerMode[i] * dBxPerMode[i]
                                             + dByPerMode[i] * dByPerMode[i]
                                             + dBzPerMode[i] * dBzPerMode[i]);
                  sorted.emplace_back(a, i);
                }
                std::sort(sorted.begin(), sorted.end(),
                          std::greater<std::pair<double, unsigned int>>());
                const unsigned int nTest = std::min<unsigned int>(10u, nlocalbfield);
                std::cout << "===== Numerical-FD closure ====="
                          << "  nFieldModes=" << nlocalbfield
                          << "  testing top-" << nTest << " modes by basis amplitude"
                          << "  dB_target=" << dBtarget << " T"
                          << "  comparing (qop, lambda, phi)" << std::endl;
                std::cout << std::scientific << std::setprecision(4);
                double worstRel = 0.0;
                for (unsigned int j = 0; j < nTest; ++j) {
                  const double basisAmp = sorted[j].first;
                  const unsigned int imode = sorted[j].second;
                  if (basisAmp < 1e-15) {
                    std::cout << "  mode " << imode
                              << " (rank " << j << ", basis=" << basisAmp
                              << "): below floor, skip" << std::endl;
                    continue;
                  }
                  const double eps = dBtarget / basisAmp;
                  const Eigen::Vector3d dBpert(
                      dB(0) + eps * dBxPerMode[imode],
                      dB(1) + eps * dByPerMode[imode],
                      dB(2) + eps * dBzPerMode[imode]);
                  auto pertResult = g4prop->propagateGenericWithJacobianAltD(
                      propInputState, surface, dBpert, dxival,
                      0., 0., -1., g4PartName);
                  if (!std::get<0>(pertResult)) {
                    std::cout << "  mode " << imode << ": pert prop failed" << std::endl;
                    continue;
                  }
                  const Eigen::Vector3d cPert = qopLamPhi(std::get<1>(pertResult));
                  const Eigen::Vector3d dCFD = (cPert - cNom) / eps;
                  // dStateDparams cols 0..2 are (qop, lambda, phi) of the
                  // 5-component curvilinear state.
                  const Eigen::Vector3d dCAn = dStateDparams.col(imode).head<3>();
                  Eigen::Vector3d rel;
                  for (int k = 0; k < 3; ++k) {
                    // qop is conserved under uniform-dB perturbation in
                    // helix transport (no eloss coupling here), so the
                    // analytic dCAn(0) is identically zero up to FP-noise
                    // (~1e-11). FD also lands at FP noise. Both small =
                    // skip the comparison; relative diff is meaningless.
                    constexpr double FP_FLOOR = 1e-10;
                    if (std::abs(dCAn(k)) < FP_FLOOR &&
                        std::abs(dCFD(k)) < FP_FLOOR) {
                      rel(k) = 0.0;
                      continue;
                    }
                    const double scale = std::max(std::abs(dCAn(k)), 1e-30);
                    rel(k) = std::abs(dCFD(k) - dCAn(k)) / scale;
                    if (rel(k) > worstRel) worstRel = rel(k);
                  }
                  std::cout << "  mode " << imode << " (rank " << j
                            << ", basis=" << basisAmp << ", eps=" << eps << ")"
                            << "  rel(qop,lam,phi)=[" << rel.transpose() << "]"
                            << std::endl;
                  std::cout << "      FD=[" << dCFD.transpose() << "]"
                            << "  an=[" << dCAn.transpose() << "]" << std::endl;
                }
                std::cout << "===== FD closure: worst rel = " << worstRel
                          << " over top-" << nTest << " modes ====="
                          << std::endl;
                std::cout.unsetf(std::ios_base::floatfield);
                didFDClosure_ = true;
              }
              // ------------------------------------------------------------

              if (ihit == 0) {
                constexpr unsigned int nvtxstate = 10;
                constexpr unsigned int nlocalstate = 5;
                const unsigned int nlocal = nvtxstate + nlocalstate + nlocalparms;

                constexpr unsigned int localvtxidx = 0;
                constexpr unsigned int localstateidx = localvtxidx + nvtxstate;
                constexpr unsigned int localparmidx = localstateidx + nlocalstate;

                constexpr unsigned int fullvtxidx = 0;
                const unsigned int fullstateidx = trackstateidx;
                const unsigned int fullparmidx = nstateparms + parmidx;

                const unsigned int vtxjacidx = 5*id;

                Matrix<double, 5, Dynamic> Fprop(5, nlocal);
                if (dolocalupdate) {
                  Fprop.middleCols<nvtxstate>(localvtxidx) = -Hm*FdFm.leftCols<5>()*twotrackpca2curvref.middleRows<5>(vtxjacidx);
                  Fprop.middleCols<nlocalstate>(localstateidx) = Hp;
                  Fprop.middleCols(localparmidx, nlocalparms) = -Hm * dStateDparams;
                }
                else {
                  Fprop.middleCols<nvtxstate>(localvtxidx) = -FdFm.leftCols<5>()*twotrackpca2curvref.middleRows<5>(vtxjacidx);
                  Fprop.middleCols<nlocalstate>(localstateidx) = Matrix<double, nlocalstate, nlocalstate>::Identity();
                  Fprop.middleCols(localparmidx, nlocalparms) = -dStateDparams;
                }

                const double propchisq = dx0.transpose()*Qinv*dx0;
                chisq0val += propchisq;

                // Sparse GBL row write: residual dx0, weight Qinv, design
                // = Fprop split into vertex-PCA + this-hit-state (Ffull)
                // and field+eloss (Jfull). (void)fullparmidx — Jfull is
                // npars-standalone so its column is parmidx, not the dense
                // nstateparms+parmidx.
                (void)fullparmidx;
                rfull.segment<5>(irow) = dx0;
                Ffull.block(irow, fullvtxidx, 5, nvtxstate) =
                    Fprop.middleCols<nvtxstate>(localvtxidx);
                Ffull.block(irow, fullstateidx, 5, nlocalstate) =
                    Fprop.middleCols<nlocalstate>(localstateidx);
                Jfull.block(irow, parmidx, 5, nlocalparms) =
                    Fprop.middleCols(localparmidx, nlocalparms);
                Vinvfull.block<5, 5>(irow, irow) = Qinv;
                if (dores && (icons == 0 || exportVarianceGrads_)) {
                  for (auto const *dQpart : {&dQMS, &dQI}) {
                    std::vector<Triplet<double>> coeffs;
                    for (unsigned int ir = 0; ir < 5; ++ir) {
                      for (unsigned int ic = 0; ic < 5; ++ic) {
                        coeffs.emplace_back(irow + ir, irow + ic, (*dQpart)(ir, ic));
                      }
                    }
                    SparseMatrix<double> &dV = dVs.emplace_back(ncons, ncons);
                    dV.setFromTriplets(coeffs.begin(), coeffs.end());
                    resblockrng.push_back({{irow, 5}});
                    resglobidx.push_back(dQpart == &dQMS ? msglobalidx : ioniglobalidx);
                    // family (parmtype) of the entry, for the cvhcf pooling
                    resfamily_.push_back(dQpart == &dQMS ? 10 : 11);
                    resvalidhit_.push_back(-1);   // material block, not a hit
                    rescls_.push_back(-1);
                  }
                  registerMatGroupNoise(irow);
                }
                irow += 5;
              }
              else {
                constexpr unsigned int nlocalstate = 10;
                const unsigned int nlocal = nlocalstate + nlocalparms;

                constexpr unsigned int localstateidx = 0;
                constexpr unsigned int localparmidx = localstateidx + nlocalstate;

                const unsigned int fullstateidx = trackstateidx + 5*(ihit - 1);
                const unsigned int fullparmidx = nstateparms + parmidx;

                Matrix<double, 5, Dynamic> Fprop(5, nlocal);
                if (dolocalupdate) {
                  Fprop.leftCols<5>() = -Hm*FdFm.leftCols<5>();
                  Fprop.middleCols<5>(5) = Hp;
                  Fprop.middleCols(localparmidx, nlocalparms) = -Hm * dStateDparams;
                }
                else {
                  Fprop.leftCols<5>() = -FdFm.leftCols<5>();
                  Fprop.middleCols<5>(5) = Matrix<double, 5, 5>::Identity();
                  Fprop.middleCols(localparmidx, nlocalparms) = -dStateDparams;
                }

                const double propchisq = dx0.transpose()*Qinv*dx0;
                chisq0val += propchisq;

                // Sparse GBL row write. nlocalstate=10 here spans the
                // previous + current hit state (Fprop.leftCols<10>);
                // localparmidx=10 is the field+eloss column group.
                (void)fullparmidx;
                rfull.segment<5>(irow) = dx0;
                Ffull.block(irow, fullstateidx, 5, nlocalstate) =
                    Fprop.leftCols(nlocalstate);
                Jfull.block(irow, parmidx, 5, nlocalparms) =
                    Fprop.middleCols(localparmidx, nlocalparms);
                Vinvfull.block<5, 5>(irow, irow) = Qinv;
                if (dores && (icons == 0 || exportVarianceGrads_)) {
                  for (auto const *dQpart : {&dQMS, &dQI}) {
                    std::vector<Triplet<double>> coeffs;
                    for (unsigned int ir = 0; ir < 5; ++ir) {
                      for (unsigned int ic = 0; ic < 5; ++ic) {
                        coeffs.emplace_back(irow + ir, irow + ic, (*dQpart)(ir, ic));
                      }
                    }
                    SparseMatrix<double> &dV = dVs.emplace_back(ncons, ncons);
                    dV.setFromTriplets(coeffs.begin(), coeffs.end());
                    resblockrng.push_back({{irow, 5}});
                    resglobidx.push_back(dQpart == &dQMS ? msglobalidx : ioniglobalidx);
                    // family (parmtype) of the entry, for the cvhcf pooling
                    resfamily_.push_back(dQpart == &dQMS ? 10 : 11);
                    resvalidhit_.push_back(-1);   // material block, not a hit
                    rescls_.push_back(-1);
                  }
                  registerMatGroupNoise(irow);
                }
                irow += 5;

              }

              for (unsigned int imode = 0; imode < nlocalbfield; ++imode) {
                globalidxv[parmidx++] = fieldCorrection_->basisGlobalIdx(imode);
              }
              if (globalMaterialModel_) {
                // One slot per material group, shared global indices
                for (unsigned int g = 0; g < nMatGroups; ++g) {
                  globalidxv[parmidx++] = matGroupGlobalIdx_[g];
                }
              } else {
                globalidxv[parmidx++] = elossglobalidx;
              }

              if (hit->isValid()) {

                //apply measurement update if applicable
                LocalTrajectoryParameters locparm(localparms[0],
                                                  localparms[1],
                                                  localparms[2],
                                                  localparms[3],
                                                  localparms[4],
                                                  localparms[5]);
                const TrajectoryStateOnSurface tsostmp(locparm, *hit->surface(), field);

                auto const& preciseHit = cloner.makeShared(hit, tsostmp);

                if (!preciseHit->isValid()) {
                  std::cout << "Abort: Failed updating hit" << std::endl;
                  ++fitFailHitUpdate_;
                  valid = false;
                  break;
                }

                const bool align2d = detidparms.count(std::make_pair(1, aligndetid));

                const Matrix<double, 2, 2> &Rglued = rgluemap_.at(preciseHit->geographicalId());
                const GloballyPositioned<double> &surfaceglued = surfacemapD_.at(parmdetid);

                auto fillAlignGrads = [&](auto Nalign) {
                  constexpr unsigned int nlocalstate = 2;
                  constexpr unsigned int localstateidx = 0;
                  constexpr unsigned int localalignmentidx = nlocalstate;
                  constexpr unsigned int localparmidx = localalignmentidx;

                  // abusing implicit template argument to pass
                  // a template value via std::integral_constant
                  constexpr unsigned int nlocalalignment = Nalign();
                  constexpr unsigned int nlocalparms = nlocalalignment;
                  constexpr unsigned int nlocal = nlocalstate + nlocalparms;
                  
                  const unsigned int fullstateidx = trackstateidx + 5*ihit + 3;
                  const unsigned int fullparmidx = nstateparms + nparsBfield + nparsEloss + alignmentparmidx;

                  const bool ispixel = GeomDetEnumerators::isTrackerPixel(preciseHit->det()->subDetector());

                  const bool hit1d = preciseHit->dimension() == 1;

                  const Matrix<double, 2, 2> Hu = Hp.bottomRightCorner<2,2>();

                  Matrix<double, 2, 1> dy0;
                  Matrix<double, 2, 2> Vinv;
                  // rotation from module to strip coordinates
                  Matrix2d R;

                  // Side-resolved pathology class of this hit (pixels only;
                  // stays 0 for strips and clean pixel hits). Bits:
                  // 0 = -x edge, 1 = +x edge, 2 = -y edge, 3 = +y edge,
                  // 4 = sizeX==1, 5 = sizeY==1. Used by the hit-diagnostic
                  // branches, the deweight mode, and the class-correction
                  // Jacobian columns appended after the alignment block.
                  // pixclsValid marks that the cluster classification
                  // succeeded (needed to tell a genuinely clean pixel hit
                  // from a failed cast, since both leave pixcls == 0).
                  int pixcls = 0;
                  bool pixclsValid = false;

                  const double lxcor = localparms[3];
                  const double lycor = localparms[4];


                  const Topology &topology = preciseHit->det()->topology();

                  // undo deformation correction
                  const LocalPoint lpnull(0., 0.);
                  const MeasurementPoint mpnull = topology.measurementPosition(lpnull);
                  const Topology::LocalTrackPred pred(tsostmp.localParameters().vector());

                  auto const defcorr = topology.localPosition(mpnull, pred) - topology.localPosition(mpnull);

                  // Rung-E closure mode (port of the single-track logic):
                  // fit SIMULATED hit positions with unchanged covariances.
                  // Per-leg sim-hit match by detid + signed particle type
                  // (mu- = +13 <-> charge -1 at updtsos[6]); only measured
                  // coordinates substituted (strip-y stays at the reco
                  // convention).
                  const PSimHit *simhit = nullptr;
                  if (fitSimHitPositions_ && doSim_) {
                    const int wanttype = updtsos[6] < 0. ? 13 : -13;
                    for (auto const &simhith : simHits) {
                      for (const PSimHit &sh : *simhith) {
                        if (sh.detUnitId() == preciseHit->geographicalId().rawId()
                            && sh.particleType() == wanttype) {
                          simhit = &sh;
                          break;
                        }
                      }
                      if (simhit != nullptr) {
                        break;
                      }
                    }
                  }
                  const bool usesimpos = simhit != nullptr;
                  const double hitxreco = preciseHit->localPosition().x() - defcorr.x();
                  const double hityreco = preciseHit->localPosition().y() - defcorr.y();
                  const double hitx = usesimpos ? simhit->localPosition().x() : hitxreco;
                  const double hity = (usesimpos && ispixel && !hit1d)
                                      ? simhit->localPosition().y() : hityreco;

                  double lyoffset = 0.;
                  double hitphival = -99.;
                  double localphival = -99.;

                  if (hit1d) {
                    const ProxyStripTopology *proxytopology = dynamic_cast<const ProxyStripTopology*>(&(preciseHit->det()->topology()));

                    dy0[0] = hitx - lxcor;
                    dy0[1] = hity - lycor;

                    const double striplength = proxytopology->stripLength();
                    const double yerr2 = striplength*striplength/12.;

                    Vinv = Matrix<double, 2, 2>::Zero();
                    Vinv(0,0) = 1./preciseHit->localPositionError().xx();
      // Vinv(1,1) = 1./yerr2;

                    R = Matrix2d::Identity();


      // std::cout << "1d hit, original x = " << preciseHit->localPosition().x() << " y = " << preciseHit->localPosition().y() << " corrected x = " << hitx << " y = " << hity << std::endl;
                  }
                  else {
                    // 2d hit
      // assert(align2d);

                    Matrix2d iV;
                    iV << preciseHit->localPositionError().xx(), preciseHit->localPositionError().xy(),
                          preciseHit->localPositionError().xy(), preciseHit->localPositionError().yy();
                    if (ispixel) {

                      dy0[0] = hitx - lxcor;
                      dy0[1] = hity - lycor;

                      Vinv = iV.inverse();

                      R = Matrix2d::Identity();

                      // Side-resolved pathology classification of this hit
                      // (see pixcls declaration above for the bit layout).
                      if (fillHitDiagnostics_ || deweightPathoHits_ ||
                          pixelHitClassCorrections_) {
                        const TrackerSingleRecHit* diagtkhit =
                            dynamic_cast<const TrackerSingleRecHit*>(&*preciseHit);
                        const PixelTopology* diagtopo =
                            dynamic_cast<const PixelTopology*>(&topology);
                        if (diagtkhit != nullptr && diagtopo != nullptr &&
                            diagtkhit->cluster_pixel().isNonnull()) {
                          const SiPixelCluster& cl = *diagtkhit->cluster_pixel();
                          if (cl.minPixelRow() == 0) pixcls |= 1 << 0;
                          if (cl.maxPixelRow() == diagtopo->nrows() - 1) pixcls |= 1 << 1;
                          if (cl.minPixelCol() == 0) pixcls |= 1 << 2;
                          if (cl.maxPixelCol() == diagtopo->ncolumns() - 1) pixcls |= 1 << 3;
                          if (cl.sizeX() <= 1) pixcls |= 1 << 4;
                          if (cl.sizeY() <= 1) pixcls |= 1 << 5;
                          pixclsValid = true;
                        }
                        // Any set bit = pathological: remove the pull but
                        // keep the surface/state so dy0 stays defined.
                        if (deweightPathoHits_ && pixcls != 0) {
                          Vinv *= 1e-6;
                        }
                        // dtanLA response weight of this hit (physics mode):
                        // size-1 = 1, x-edge = lorentzWedge, regular = wclean.
                        const double lorentzW = (pixcls & 0x10) ? lorentzWsize1_
                            : ((pixcls & 0x3) ? lorentzWedge_ : lorentzWclean_);
                        const double lorentzScale = pixclsValid
                            ? 0.5 * preciseHit->det()->surface().bounds().thickness()
                            : 0.;

                        // Injection test: simulate a true Lorentz-angle
                        // mismatch by shifting every valid pixel hit's
                        // local-x with the injected response weights.
                        if (injectLorentzTan_ != 0. && pixclsValid) {
                          const double winj = (pixcls & 0x10) ? lorentzWsize1_
                              : ((pixcls & 0x3) ? lorentzWedge_
                                 : (injectLorentzWclean_ > -900. ? injectLorentzWclean_
                                                                 : lorentzWclean_));
                          dy0[0] += lorentzScale * winj * injectLorentzTan_;
                        }

                        // Apply the current class-correction values (seeded
                        // from corFiles) to the residual: dy0 += J*theta
                        // with the same columns as appended to Jfull below,
                        // so a fitted theta zeroes the class-param gradients.
                        if (pixelHitClassCorrections_ &&
                            (pixcls != 0 || (pixelLorentzParam_ && pixclsValid))) {
                          const DetId pixdetid = preciseHit->geographicalId();
                          auto corval = [&](unsigned int pt) {
                            return corparms_[detidparms.at(std::make_pair(pt, pixdetid))];
                          };
                          if (pixelLorentzParam_ && pixclsValid) {
                            dy0[0] += lorentzScale * lorentzW * corval(22);
                          }
                          if (pixcls & 0x3) {
                            const double s = (pixcls & 0x2) ? 1. : -1.;
                            if (!pixelLorentzParam_) dy0[0] += corval(16);
                            dy0[0] += s * corval(17);
                          }
                          if (pixcls & 0xc) {
                            const double s = (pixcls & 0x8) ? 1. : -1.;
                            dy0[1] += corval(18) + s * corval(19);
                          }
                          if ((pixcls & 0x10) && !pixelLorentzParam_) dy0[0] += corval(20);
                          if (pixcls & 0x20) dy0[1] += corval(21);
                        }
                        if (fillHitDiagnostics_) {
                          hitdiag_detid.push_back(preciseHit->geographicalId().rawId());
                          hitdiag_trk.push_back(static_cast<int>(id));
                          hitdiag_charge.push_back(trackPair[id]->charge());
                          hitdiag_class.push_back(pixcls);
                          hitdiag_dx.push_back(dy0[0]);
                          hitdiag_dy.push_back(dy0[1]);
                          hitdiag_exx.push_back(iV(0, 0));
                          hitdiag_eyy.push_back(iV(1, 1));
                          hitdiag_lx.push_back(hitx);
                          hitdiag_ly.push_back(hity);
                        }
                      }
                    }
                    else {
                      // transform to polar coordinates to end the madness
                      //TODO handle the module deformations consistently here (currently equivalent to dropping/undoing deformation correction)

      // std::cout << "wedge\n" << std::endl;

                      const ProxyStripTopology *proxytopology = dynamic_cast<const ProxyStripTopology*>(&(preciseHit->det()->topology()));

                      const TkRadialStripTopology *radialtopology = dynamic_cast<const TkRadialStripTopology*>(&proxytopology->specificTopology());

                      const double rdir = radialtopology->yAxisOrientation();
                      const double radius = radialtopology->originToIntersection();

                      // Wedge modules measure LOCAL PHI: in sim-position mode
                      // the residual must use the FULL sim phi (sim x AND
                      // sim y; mixing sim-x with strip-center-y gives phi
                      // errors ~ x*dy/r, mm of arc).
                      double phihit = rdir*std::atan2(hitx, rdir*hity + radius);
                      double rhohit = std::sqrt(hitx*hitx + std::pow(rdir*hity + radius, 2));
                      if (usesimpos) {
                        const double lxs = simhit->localPosition().x();
                        const double lys = simhit->localPosition().y();
                        phihit = rdir*std::atan2(lxs, rdir*lys + radius);
                        rhohit = std::sqrt(lxs*lxs + std::pow(rdir*lys + radius, 2));
                      }

                      // invert original calculation of covariance matrix to
                      // extract variance on polar angle. MUST use the RECO
                      // phi: the CPE built xx around the reco position, so
                      // tt > 0 is only guaranteed there (sim-based tan(phi)
                      // can drive it negative -> indefinite system).
                      const double phihitreco = rdir*std::atan2(hitxreco, rdir*hityreco + radius);
                      const double detHeight = radialtopology->detHeight();
                      const double radsigma = detHeight*detHeight/12.;

                      const double t1 = std::tan(phihitreco);
                      const double t2 = t1*t1;

                      const double tt = preciseHit->localPositionError().xx() - t2*radsigma;

                      const double phierr2 = tt / std::pow(radialtopology->centreToIntersection(), 2);

                      const double striplength = detHeight * std::sqrt(1. + std::pow( hitx/(rdir*hity + radius), 2) );

                      const double rhoerr2 = striplength*striplength/12.;


      // std::cout << "rhohit = " << rhohit << " rhobar = " << rhobar << " rhoerr2lin = " << rhoerr2lin << " rhoerr2 = " << rhoerr2 << std::endl;

                      // TODO apply (inverse) corrections for module deformations here? (take into account for jacobian?)
                      const double phistate = rdir*std::atan2(lxcor, rdir*lycor + radius);
                      const double rhostate = std::sqrt(lxcor*lxcor + std::pow(rdir*lycor + radius, 2));

                      Vinv = Matrix<double, 2, 2>::Zero();
                      Vinv(0, 0) = 1./phierr2;
      // Vinv(1, 1) = 1./rhoerr2lin;

                      // jacobian from localx-localy to localphi-localrho
                      R = Matrix2d::Zero();

                      const double yp = rdir*lycor + radius;
                      const double invden = 1./(lxcor*lxcor + yp*yp);

                      // dphi / dx
                      R(0, 0) = rdir*yp*invden;
                      // dphi / dy
                      R(0, 1) = -lxcor*invden;
                      // drho / dx
                      R(1, 0) = lxcor/rhostate;
                      // drho / dy
                      R(1, 1) = rdir*(rdir*lycor + radius)/rhostate;


                      dy0[0] = phihit - phistate;
                      dy0[1] = rhohit - rhostate;

      // std::cout << "wedge hit, original x = " << preciseHit->localPosition().x() << " y = " << preciseHit->localPosition().y() << " corrected x = " << hitx << " y = " << hity << std::endl;

                    }
                  }

                  // alignment jacobian
                  Matrix<double, 2, 6> Aval = Matrix<double, 2, 6>::Zero();

                  // const Matrix<double, 6, 1> &localparmsalign = alignGlued_ ? globalToLocal(updtsos, surfaceglued) : localparms;

                  Matrix<double, 6, 1> localparmsalign = localparms;
                  if (alignGlued_) {
                    localparmsalign = globalToLocal(updtsos, surfaceglued);
                  }

                  const double localqopval = localparmsalign[0];
                  const double localdxdzval = localparmsalign[1];
                  const double localdydzval = localparmsalign[2];
                  const double localxval = localparmsalign[3];
                  const double localyval = localparmsalign[4];

                  const double localxvalorig = localparms[3];
                  const double localyvalorig = localparms[4];

                  //standard case

                  // dx/dx
                  Aval(0,0) = -1.;
                  // dy/dy
                  Aval(1,1) = -1.;
                  // dx/dz
                  Aval(0,2) = localdxdzval;
                  // dy/dz
                  Aval(1,2) = localdydzval;
                  // dx/dtheta_x
                  Aval(0,3) = localyval*localdxdzval;
                  // dy/dtheta_x
                  Aval(1,3) = localyval*localdydzval;
                  // dx/dtheta_y
                  Aval(0,4) = -localxval*localdxdzval;
                  // dy/dtheta_y
                  Aval(1,4) = -localxval*localdydzval;
                  // dx/dtheta_z
                  Aval(0,5) = localyvalorig;
                  // dy/dtheta_z
                  Aval(1,5) = -localxvalorig;

                  // const Matrix<double, 2, 6> &A = alignGlued_ ? Rglued*Aval : Aval;

                  Matrix<double, 2, 6> A = Aval;
                  if (alignGlued_) {
                    // glued alignment dofs only for out-of-plance
                    A.middleCols<3>(2) = Rglued*Aval.middleCols<3>(2);
                  }

                  double thetaincidence = std::asin(1./std::sqrt(std::pow(localdxdzval,2) + std::pow(localdydzval,2) + 1.));

      // bool morehitquality = applyHitQuality_ ? thetaincidence > 0.25 : true;
                  bool morehitquality = true;

                  if (morehitquality) {
                    nvalidFinalarr[id]++;
                    if (ispixel) {
                      nvalidpixelFinalarr[id]++;
                    }
                  }
                  else {
                    Vinv = Matrix<double, 2, 2>::Zero();
                  }

                  constexpr std::array<unsigned int, 6> alphaidxs = {{0, 2, 3, 4, 5, 1}};

                  Matrix<double, 2, nlocal> Fhit;
                  //TODO figure out why templated version doesn't work here (gcc bug?)
                  Fhit.leftCols(2) = -R*Hu;

                  for (unsigned int ialign = 0; ialign < nlocalalignment; ++ialign) {
                    Fhit.col(ialign + 2) = -R*A.col(alphaidxs[ialign]);
                  }

                  const double hitchisq = dy0.transpose()*Vinv*dy0;
                  chisq0val += hitchisq;

                  // Sparse GBL row write: 2-dim measurement residual dy0,
                  // weight Vinv; design = Fhit split into this-hit local
                  // x/y state (Ffull, cols [fullstateidx, +2)) and the
                  // alignment dofs (Jfull, col group at
                  // nparsBfield+nparsEloss+alignmentparmidx). The dense
                  // fullparmidx = nstateparms + that group.
                  (void)fullparmidx;
                  rfull.segment<2>(irow) = dy0;
                  Ffull.block(irow, fullstateidx, 2, nlocalstate) =
                      Fhit.middleCols(localstateidx, nlocalstate);
                  Jfull.block(irow,
                              nparsBfield + nparsEloss + alignmentparmidx,
                              2, nlocalparms) =
                      Fhit.middleCols(localparmidx, nlocalparms);
                  Vinvfull.block<2, 2>(irow, irow) = Vinv;

                  // ---- HIT-RESOLUTION BLOCKS (parmtype 8/9) --------------
                  //
                  // Registering the hit blocks here is what gives the mass
                  // functional hit shares to weigh, so the per-hit-class
                  // resolution parameters can be fitted from two-track
                  // candidates and not only from single tracks.
                  //
                  // `dVs` in THIS maker feeds nothing but the influence
                  // export, so registering a block cannot move the fit; and
                  // the parmtype-8/9 corrections are deliberately NOT applied
                  // to the covariance here (the single-track maker scales
                  // `iV` by exp(corparms)) -- that WOULD be a fit change, and
                  // this export is a derivative evaluated at the covariance
                  // the fit actually used.
                  //
                  // dV is the derivative of the MEASUREMENT-frame covariance
                  // w.r.t. the log-resolution parameter, split between the
                  // local-x and local-y families so the two sum to the full
                  // covariance -- the same decomposition the single-track
                  // maker makes, read here off `Vinv` because the rows are
                  // already in the measurement frame (local x/y on pixels,
                  // local phi on wedges, local x on rectangular strips).
                  // `icons == 0` is the influence export's pass; with the
                  // log-det term on, the exported grad/hess come from the
                  // LAST pass (icons == 1 under doMassConstraint), so the
                  // blocks have to exist there too or the term would be
                  // silently dropped.
                  if (dores && exportHitResBlocks_ && (icons == 0 || exportVarianceGrads_)) {
                    const bool pix2d = ispixel && !hit1d;
                    Matrix2d hitcov = Matrix2d::Zero();
                    bool hitcovok = false;
                    if (pix2d) {
                      const double det = Vinv(0, 0) * Vinv(1, 1) - Vinv(0, 1) * Vinv(1, 0);
                      if (std::abs(det) > 0.) {
                        hitcov = Vinv.inverse();
                        hitcovok = true;
                      }
                    } else if (Vinv(0, 0) > 0.) {
                      hitcov(0, 0) = 1. / Vinv(0, 0);
                      hitcovok = true;
                    }
                    if (hitcovok) {
                      const DetId hitdetid = preciseHit->geographicalId();
                      const unsigned int nb = pix2d ? 2u : 1u;
                      // the class variables, exactly `hitres_classes.class_of`
                      int clsSizeX = 1;
                      int clsQBin = -99;
                      float clsUProj = -99.f;
                      const TrackerSingleRecHit *clstkhit =
                          dynamic_cast<const TrackerSingleRecHit *>(&*preciseHit);
                      if (clstkhit != nullptr) {
                        if (ispixel) {
                          if (clstkhit->cluster_pixel().isNonnull()) {
                            clsSizeX = clstkhit->cluster_pixel()->sizeX();
                          }
                          const SiPixelRecHit *pxh = dynamic_cast<const SiPixelRecHit *>(clstkhit);
                          if (pxh != nullptr) {
                            clsQBin = pxh->qBin();
                          }
                        } else if (clstkhit->cluster_strip().isNonnull()) {
                          clsSizeX = clstkhit->cluster_strip()->amplitudes().size();
                          const StripGeomDetUnit *stripdu =
                              dynamic_cast<const StripGeomDetUnit *>(clstkhit->det());
                          if (stripCPEForExport != nullptr && stripdu != nullptr) {
                            clsUProj = stripCPEForExport->getAlgoParam(*stripdu, locparm).afullProjection;
                          }
                        }
                      }
                      const int subdet = hitdetid.subdetId();

                      // local x / phi
                      {
                        std::vector<Triplet<double>> coeffs;
                        coeffs.emplace_back(irow, irow, hitcov(0, 0));
                        if (pix2d) {
                          coeffs.emplace_back(irow, irow + 1, 0.5 * hitcov(0, 1));
                          coeffs.emplace_back(irow + 1, irow, 0.5 * hitcov(0, 1));
                        }
                        SparseMatrix<double> &dVx = dVs.emplace_back(ncons, ncons);
                        dVx.setFromTriplets(coeffs.begin(), coeffs.end());
                        resblockrng.push_back({{irow, nb}});
                        resglobidx.push_back(detidparms.at(std::make_pair(8, hitdetid)));
                        resfamily_.push_back(8);
                        resvalidhit_.push_back(int(id));
                        rescls_.push_back(hitResClassIndex(subdet, clsSizeX, clsUProj, clsQBin, false));
                      }
                      // local y (pixels only)
                      if (pix2d) {
                        std::vector<Triplet<double>> coeffs;
                        coeffs.emplace_back(irow, irow + 1, 0.5 * hitcov(0, 1));
                        coeffs.emplace_back(irow + 1, irow, 0.5 * hitcov(0, 1));
                        coeffs.emplace_back(irow + 1, irow + 1, hitcov(1, 1));
                        SparseMatrix<double> &dVy = dVs.emplace_back(ncons, ncons);
                        dVy.setFromTriplets(coeffs.begin(), coeffs.end());
                        resblockrng.push_back({{irow, 2u}});
                        resglobidx.push_back(detidparms.at(std::make_pair(9, hitdetid)));
                        resfamily_.push_back(9);
                        resvalidhit_.push_back(int(id));
                        rescls_.push_back(hitResClassIndex(subdet, clsSizeX, clsUProj, clsQBin, true));
                      }
                    }
                  }

                  irow += 2;
                  
                  for (unsigned int idim=0; idim<nlocalalignment; ++idim) {
                    const unsigned int iidx = alphaidxs[idim];
                    const DetId ialigndetid = iidx > 1 && iidx < 5 ? aligndetid : preciseHit->geographicalId();
                    const unsigned int xglobalidx = detidparms.at(std::make_pair(iidx, ialigndetid));
                    globalidxv[nparsBfield + nparsEloss + alignmentparmidx] = xglobalidx;
                    alignmentparmidx++;
                    if (alphaidxs[idim]==0) {
                      hitidxv.push_back(xglobalidx);
                    }
                  }

                  // Pixel pathological-hit class corrections (parmtypes
                  // 16-21): extra local-translation Jacobian columns gated
                  // on this hit's class, appended after the standard
                  // alignment dofs. Column = d(residual)/d(param): (1,0)
                  // for local-x-type params, (0,1) for local-y-type (same
                  // sign convention as the parmtype-0/1 columns: Fhit.col =
                  // -R*A.col with R = identity on pixels), and the diff
                  // params additionally carry the edge-side sign s (+1 at
                  // the hi boundary, -1 at lo). Counted in nparsPixClass.
                  if (pixelHitClassCorrections_ && ispixel &&
                      (pixcls != 0 || (pixelLorentzParam_ && pixclsValid))) {
                    const DetId pixdetid = preciseHit->geographicalId();
                    auto appendClassCol = [&](unsigned int parmtype, int coord, double sign) {
                      const unsigned int xglobalidx =
                          detidparms.at(std::make_pair(parmtype, pixdetid));
                      Jfull(irow - 2 + coord,
                            nparsBfield + nparsEloss + alignmentparmidx) = sign;
                      globalidxv[nparsBfield + nparsEloss + alignmentparmidx] = xglobalidx;
                      alignmentparmidx++;
                    };
                    if (pixelLorentzParam_ && pixclsValid) {
                      // dtanLA column on every valid pixel hit: J =
                      // (t/2) * w(class) on the local-x residual row.
                      const double w = (pixcls & 0x10) ? lorentzWsize1_
                          : ((pixcls & 0x3) ? lorentzWedge_ : lorentzWclean_);
                      appendClassCol(22, 0,
                          0.5 * preciseHit->det()->surface().bounds().thickness() * w);
                    }
                    if (pixcls & 0x3) {                    // edge in x
                      const double s = (pixcls & 0x2) ? 1. : -1.;
                      if (!pixelLorentzParam_) appendClassCol(16, 0, 1.);  // edge-x-mean
                      appendClassCol(17, 0, s);            // edge-x-diff
                    }
                    if (pixcls & 0xc) {                    // edge in y
                      const double s = (pixcls & 0x8) ? 1. : -1.;
                      appendClassCol(18, 1, 1.);           // edge-y-mean
                      appendClassCol(19, 1, s);            // edge-y-diff
                    }
                    if ((pixcls & 0x10) && !pixelLorentzParam_)
                      appendClassCol(20, 0, 1.);           // sizeX1
                    if (pixcls & 0x20) appendClassCol(21, 1, 1.);   // sizeY1
                  }
                };

                if (align2d) {
                  fillAlignGrads(std::integral_constant<unsigned int, 6>());
                }
                else {
                  fillAlignGrads(std::integral_constant<unsigned int, 5>());
                }
              }
            }

            if (!valid || retryIter) {
              break;
            }

            trackstateidx += 5*tracknhits;
          }

          if (retryIter) {
            // Restore the iteration-entry linearization state and redo this
            // iteration with the adjusted input.
            refftsarr = refftsarrSnap;
            layerStatesarr = layerStatesSnap;
            if (iiter > 0) {
              dxfull *= 0.5;
              ++nBacktracks;
              ++fitStepBacktracked_;
              iiter -= 1;  // loop ++ redoes the same iteration
            } else {
              // end-of-range daughter: inflate its seed momentum and let the
              // hits pull it back down
              refftsarr[retryFailId].segment<3>(3) *= 1.25;
              ++nSeedInflations;
              ++fitSeedInflated_;
              // unsigned wrap: ++ brings iiter back to 0
              iiter = std::numeric_limits<unsigned int>::max();
            }
            continue;
          }

          if (!valid) {
            break;
          }
          
    // MatrixXd massjac;

          const Matrix<double, 7, 1> &refFts0 = refftsarr[0];
          const Matrix<double, 7, 1> &refFts1 = refftsarr[1];    
          
          const Matrix<double, 6, 6> mhess = massHessianAltD(refFts0, refFts1, massForConstraintHelpers);
          const Matrix<double, 6, 6> mhessinvsq = massinvsqHessianAltD(refFts0, refFts1, massForConstraintHelpers);
          
          const double dmassconv = iiter > 0 ? 0.5*(mhess*covrefmom).trace() : 0.;
          const double dmassconvinvsq = iiter > 0 ? 0.5*(mhessinvsq*covrefmom).trace() : 0.;
          
          
          dmassconvval = dmassconv;
          dinvmasssqconvval = dmassconvinvsq;
          
          if (icons == 0) {
            dmassconvval_cons0 = dmassconv;
            dinvmasssqconvval_cons0 = dmassconvinvsq;
          }
          
          // add mass constraint to General Broken Lines (GBL) fit
          if (icons > 0) {
            //TODO simplify this to treat the 6 parameters contiguously (now that they are contiguous in the original vector)

            constexpr unsigned int nlocalstate0 = 3;
            constexpr unsigned int nlocalstate1 = 3;
            
            constexpr unsigned int nlocal = nlocalstate0 + nlocalstate1;
            
            constexpr unsigned int localstateidx0 = 0;
            constexpr unsigned int localstateidx1 = localstateidx0 + nlocalstate0;
            
            const unsigned int fullstateidx0 = 0;
            const unsigned int fullstateidx1 = 3;
            
            using MScalar = AANT<double, nlocal>;
            
            
            const ROOT::Math::PxPyPzMVector mom0(refFts0[3],
                                                    refFts0[4],
                                                    refFts0[5],
                                                    trackMass[0]);
            
            const ROOT::Math::PxPyPzMVector mom1(refFts1[3],
                                                    refFts1[4],
                                                    refFts1[5],
                                                    trackMass[1]);
            
            const double massval = (mom0 + mom1).mass();
            
            const Matrix<double, 1, 6> mjacalt = massJacobianAltD(refFts0, refFts1, massForConstraintHelpers);

// const double dmsq0 = massval - massconstraintval;
            const double dmsq0 = massval - massconstraintval - dmassconv;

            const Matrix<double, 1, 6> &Fmass = mjacalt;

            const double invSigmaMsq = 1./massConstraintWidth_/massConstraintWidth_;

            const double masschisq = dmsq0*dmsq0*invSigmaMsq;
            chisq0val += masschisq;

            // Sparse GBL row write: 1 J/psi-mass row (icons==1 pass only).
            // residual = dmsq0, weight = invSigmaMsq (scalar). Fmass (1x6)
            // splits into the two daughters' qop/lam/phi state blocks at
            // fullstateidx0=0 and fullstateidx1=3.
            rfull(irow) = dmsq0;
            Ffull.block(irow, fullstateidx0, 1, nlocalstate0) =
                Fmass.block(0, localstateidx0, 1, nlocalstate0);
            Ffull.block(irow, fullstateidx1, 1, nlocalstate1) =
                Fmass.block(0, localstateidx1, 1, nlocalstate1);
            Vinvfull(irow, irow) = invSigmaMsq;
            irow += 1;

          }

         

          
          
    // std::cout << nhits << std::endl;
    // std::cout << nvalid << std::endl;
    // std::cout << nvalidalign2d << std::endl;
    // std::cout << nparsAlignment << std::endl;
    // std::cout << alignmentparmidx << std::endl;
    // 
    // std::cout << nparsBfield << std::endl;
    // std::cout << nparsEloss << std::endl;
    // std::cout << parmidx << std::endl;
          
          assert(trackstateidx == nstateparms);
          assert(parmidx == (nparsBfield + nparsEloss));
          assert(alignmentparmidx == nparsAlignment);
          // Sparse GBL: the running constraint-row cursor must land
          // exactly on the pre-computed ncons, else rfull/Ffull/Jfull/
          // Vinvfull rows are misaligned. Cheap permanent invariant.
          assert(irow == ncons);
          
    // if (nhits != nvalid) {
    // continue;
    // }

          // Sparse GBL: the dense freezeparm() (deweight via 1e6 diagonal)
          // is replaced by excluding the index from freestateidxs --
          // already applied at the freestatemask construction near the
          // top (fitFromGenParms_ -> vtx 0..9, doVtxConstraint_ -> idx 6,
          // fitFromSimParms_ -> all). Nothing to do here.

// if (fitFromGenParms_) {
// freezeparm(2);
// for (unsigned int id = 0; id < 2; ++id) {
// freezeparm(trackstateidxarr[id] + 1);
// }
// }
          
  // if (fitFromGenParms_) {
  // for (unsigned int id = 0; id < 2; ++id) {
  // for (unsigned int i=1; i<3; ++i) {
  // freezeparm(trackstateidxarr[id] + i);
  // }
  // }
  // }
          
          //now do the expensive calculations and fill outputs
          
          //symmetrize the matrix (previous block operations do not guarantee that the needed blocks are filled)
          //TODO handle this more efficiently?
    // hessfull.triangularView<StrictlyLower>() = hessfull.triangularView<StrictlyUpper>().transpose();
          
    // for (unsigned int i=0; i<3; ++i) {
    // gradfull[i] = 0.;
    // hessfull.row(i) *= 0.;
    // hessfull.col(i) *= 0.;
    // hessfull(i,i) = 1e6;
    // }
          
    // for (auto trackstateidx : trackstateidxarr) {
    // for (unsigned int i = trackstateidx; i < (trackstateidx + 1); ++i) {
    // gradfull[i] = 0.;
    // hessfull.row(i) *= 0.;
    // hessfull.col(i) *= 0.;
    // hessfull(i,i) = 1e6;
    // }
    // }
          
    // {
    // unsigned int i = trackstateidxarr[1];
    // gradfull[i] = 0.;
    // hessfull.row(i) *= 0.;
    // hessfull.col(i) *= 0.;
    // hessfull(i,i) = 1e6; 
    // }
    // 
          
    // std::cout << "gradfull:" << std::endl;
    // std::cout << gradfull << std::endl;
    // 
    // std::cout << "gradfull.head(nstateparms):" << std::endl;
    // std::cout << gradfull.head(nstateparms) << std::endl;
    // 
    // std::cout << "gradfull.tail(npars):" << std::endl;
    // std::cout << gradfull.tail(npars) << std::endl;
    // 
    // std::cout << "hessfull.diagonal():" << std::endl;
    // std::cout << hessfull.diagonal() << std::endl;
          
          // ---- chi2-based (Armijo) retroactive backtracking ---------------
          // chisq0val is now complete (beamspot + pointing + propagation +
          // hits + mass) at the CURRENT linearization point, i.e. it is the
          // REALIZED chi2 of the step taken at the end of the previous
          // iteration. If it fails the sufficient-decrease test, restore that
          // linearization (the same snapshot the failed-leg retry uses), halve
          // the step and redo the iteration. Costs no extra propagation on the
          // accept path. The slack absorbs the chi2 wobble from relinearization
          // (the propagation/material model and the mass convolution term are
          // re-evaluated at the new state), so only a genuine blow-up fires.
          if (stepBacktracking_ && iiter >= stepBacktrackFromIter_ && std::isfinite(chisq0valPrev) &&
              nChi2Bt < maxChi2Backtrack_ && nChi2BtTotal < maxChi2Backtrack_ * niters) {
            const double thresh = chisq0valPrev + armijoC_ * predDecrPrev +
                                  armijoSlack_ * std::max(1., std::abs(chisq0valPrev));
            if (!(chisq0val <= thresh)) {
              refftsarr = refftsarrSnap;
              layerStatesarr = layerStatesSnap;
              dxfull *= 0.5;
              predDecrPrev *= 0.5;
              ++nChi2Bt;
              ++nChi2BtTotal;
              ++stepBtChi2Events_;
              if (!stepBtChi2ThisFit) {
                stepBtChi2ThisFit = true;
                ++fitStepBtChi2_;
              }
              if (stepPrints_ < stepPrintLimit_) {
                ++stepPrints_;
                std::cout << "GN step backtracked (two-track): icons = " << icons
                          << " iiter = " << iiter
                          << " chisq " << chisq0valPrev << " -> " << chisq0val
                          << " (thresh " << thresh << ") nbt = " << nChi2Bt
                          << " seed0(q,pt,eta)=(" << itrack->charge() << "," << itrack->pt()
                          << "," << itrack->eta() << ")"
                          << " seed1(q,pt,eta)=(" << jtrack->charge() << "," << jtrack->pt()
                          << "," << jtrack->eta() << ")" << std::endl;
              }
              iiter -= 1;  // loop ++ redoes the same iteration from the snapshot
              continue;
            }
          }
          // step accepted
          nChi2Bt = 0;
          chisq0valPrev = chisq0val;

          // Sparse GBL solve (replaces dense Cinvd=LDLT(2 Fs^T Vinv Fs);
          // dxfull=-Cinvd.solve(2 Fs^T Vinv r)). Mathematically identical:
          // the factor of 2 cancels in -(Fs^T Vinv Fs)^-1 Fs^T Vinv r.
          Fsparse = Ffull(Eigen::placeholders::all, freestateidxs).sparseView();
          Vinvsparse = Vinvfull.sparseView();
          VinvF = Vinvsparse * Fsparse;

          Cinvd.compute(Fsparse.transpose() * VinvF);

          dxfree = -Cinvd.solve(VinvF.transpose() * rfull);

          // Fail fast on a non-finite update: a drained/runaway propagation
          // leg (or a singular normal matrix) yields NaN/inf here, which
          // previously leaked into the charge-sum check and was miscounted
          // as a charge flip.
          if (!dxfree.allFinite()) {
            // Localize the poison: which solve input went non-finite, and
            // the displaced-vertex geometry of the candidate.
            const bool rOk = rfull.allFinite();
            const bool fOk = Eigen::MatrixXd(Ffull).allFinite();
            const bool vOk = Eigen::MatrixXd(Vinvfull).allFinite();
            const double r0 = std::hypot(refftsarr[0][0], refftsarr[0][1]);
            const double r1 = std::hypot(refftsarr[1][0], refftsarr[1][1]);
            std::cout << "Abort: non-finite parameter update from solve!"
                      << " icons = " << icons << " iiter = " << iiter
                      << " finite(r,F,Vinv)=(" << rOk << "," << fOk << "," << vOk << ")"
                      << " refR=(" << r0 << "," << r1 << ")"
                      << " nhits=(" << nhitsarr[0] << "," << nhitsarr[1] << ")"
                      << " nvalid=(" << nvalidarr[0] << "," << nvalidarr[1] << ")"
                      << " seed0(q,pt,eta)=(" << itrack->charge() << "," << itrack->pt() << "," << itrack->eta() << ")"
                      << " seed1(q,pt,eta)=(" << jtrack->charge() << "," << jtrack->pt() << "," << jtrack->eta() << ")"
                      << std::endl;
            ++fitFailNaN_;
            valid = false;
            break;
          }

          dxfull = VectorXd::Zero(nstateparms);
          dxfull(freestateidxs) = dxfree;

          // Momentum-floor safeguard on the Gauss-Newton step (port of the
          // single-track clamp). Only fatal update outcomes are prevented --
          // a daughter's q/p sign flipping, or its momentum dropping below
          // the configurable floor -- by scaling the WHOLE joint step vector
          // (both tracks + vertex are one coupled system: one common scale,
          // direction preserved). Converts most hard aborts at the
          // charge-sum check below into recoverable (or cleanly at-cap)
          // fits; essential for soft V0 daughters.
          {
            const Matrix<double, 10, 1> statepcaref =
                twoTrackCart2pca(refftsarr[0], refftsarr[1]);
            double stepscale = 1.;
            for (unsigned int id = 0; id < 2; ++id) {
              const double qopref = statepcaref[3 * id];
              const double dqop = dxfull[3 * id];
              if (qopref == 0. || dqop == 0.) {
                continue;
              }
              double s = 1.;
              if (maxMomentumStepFactor_ > 1.) {
                // Relative trust region in q/p, see the member comment. The
                // bound is always strictly inside p_ref, so a soft daughter is
                // neither pinned at the floor nor frozen at its seed. No charge
                // flip is permitted here (the two-track fit has no
                // ambiguous-charge use case), which the qopFlipAllow = 0
                // argument expresses.
                s = cvhstep::legStepScaleRel(qopref, dqop, clampMomentumFloor_,
                                             maxMomentumStepFactor_, 0., nullptr);
              } else {
                // LEGACY absolute-floor-only clamp
                const double qopupd = qopref + dqop;
                if (qopupd * qopref <= 0.) {
                  // sign flip: stop half-way toward q/p = 0
                  s = -0.5 * qopref / dqop;
                } else if (std::abs(qopupd) > 1. / clampMomentumFloor_) {
                  // p_upd below the floor: land exactly on p = floor, same charge
                  s = (std::copysign(1. / clampMomentumFloor_, qopref) - qopref) / dqop;
                }
              }
              if (s < stepscale) {
                stepscale = s;
              }
            }
            if (stepscale < 1.) {
              stepscale = std::max(stepscale, 0.);
              dxfree *= stepscale;
              dxfull *= stepscale;
              stepScaleApplied *= stepscale;
              ++stepClampEvents_;
              if (!stepClampedThisFit) {
                stepClampedThisFit = true;
                ++fitStepClamped_;
              }
              if (stepPrints_ < stepPrintLimit_) {
                ++stepPrints_;
                std::cout << "GN step clamped (two-track): icons = " << icons
                          << " iiter = " << iiter << " scale = " << stepscale
                          << " seed0(q,pt,eta)=(" << itrack->charge() << "," << itrack->pt()
                          << "," << itrack->eta() << ")"
                          << " seed1(q,pt,eta)=(" << jtrack->charge() << "," << jtrack->pt()
                          << "," << jtrack->eta() << ")" << std::endl;
              }
            }
          }

// std::cout << "dxfull vtx: " << dxfull.head<3>() << std::endl;
          
  // dxdparms = -Cinvd.solve(d2chisqdxdparms).transpose();
          
      // if (debugprintout_) {
      // std::cout << "dxrefdparms" << std::endl;
      // std::cout << dxdparms.leftCols<5>() << std::endl;
      // }
          
  // grad = dchisqdparms + dxdparms*dchisqdx;
          //TODO check the simplification
      // hess = d2chisqdparms2 + 2.*dxdparms*d2chisqdxdparms + dxdparms*d2chisqdx2*dxdparms.transpose();
  // hess = d2chisqdparms2 + dxdparms*d2chisqdxdparms;
          
          // Sparse GBL deltachi2: dense g^T dx + 0.5 dx^T H dx with
          // g=2Fs^TVinvr, H=2Fs^TVinvFs, dx=-H^-1 g reduces algebraically
          // to -(Fs^TVinvr)^T(Fs^TVinvFs)^-1(Fs^TVinvr) = rfull^T VinvF dxfree.
          const double deltachisq = (rfull.transpose() * VinvF * dxfree)(0, 0);

// std::cout << "iiter = " << iiter << ", deltachisq = " << deltachisq << std::endl;
// // 
// SelfAdjointEigenSolver<MatrixXd> es(d2chisqdx2, EigenvaluesOnly);
// const double condition = es.eigenvalues()[nstateparms-1]/es.eigenvalues()[0];
// std::cout << "eigenvalues:" << std::endl;
// std::cout << es.eigenvalues().transpose() << std::endl;
// std::cout << "condition: " << condition << std::endl;
          
          // Quadratic-model chi2 change of the step ACTUALLY applied. dxfree
          // has already been scaled by stepScaleApplied, so deltachisq = t*d
          // with d the full-step value and the model change is
          // d*(2t - t^2) = deltachisq*(2 - t).
          predDecrPrev = deltachisq * (2. - stepScaleApplied);

          chisqval = chisq0val + deltachisq;

          deltachisqval = chisq0val + deltachisq - chisqvalold;

          chisqvalold = chisq0val + deltachisq;
          
// ndof = 5*nhits + nvalid + nvalidalign2d - nstateparms;
          // COMPUTED SIGNED, THEN CLAMPED.  `ndof` is an unsigned member and
          // `nstateparms` is 10 + 5*nhits, so this expression is really
          // nvalid + nvalidpixel - 10 (+ the constraint rows): evaluated in
          // unsigned arithmetic it UNDERFLOWS to ~4e9 for a pair whose two
          // legs together carry fewer than ten valid-hit-equivalents, and
          // that value then reads as "full rank" in the factored-Hessian
          // block below (`nrank = min(ndof, nparsfinal)`).  Keeping the
          // arithmetic signed makes the degenerate case visible instead of
          // wrapping it.  For every candidate with a positive ndof the value
          // stored is bit-identical to the old expression.
          long long ndofsigned = 5LL*(long long)nhits + (long long)nvalid
                                 + (long long)nvalidpixel - (long long)nstateparms;

          if (bsConstraint_) {
            ndofsigned += 3;
          }

          if (doPointingConstraint_) {
            // 2D-transverse pointing adds 1 scalar constraint
            ++ndofsigned;
          }

          if (doVtxConstraint_) {
            ++ndofsigned;
          }

          if (icons == 1) {
            ++ndofsigned;
          }

          ndof = ndofsigned > 0 ? (unsigned int)ndofsigned : 0u;

          // A FIT WITH NO DEGREES OF FREEDOM IS A FAILED FIT, not a fit with
          // an empty rank.  The factored-Hessian export takes exactly the top
          // `nrank = min(ndof, nparsfinal)` eigenmodes of the mean Hessian and
          // reads `eigvals(nparsfinal - nrank)` to report the truncation gap;
          // at ndof == 0 that indexes ONE PAST THE END of the length-
          // nparsfinal eigenvalue vector and aborts inside Eigen
          // (DenseCoeffsBase::operator()'s `index < size()`), taking the whole
          // job with it -- a crashed cmsRun output has NO KEYS, so the whole
          // task is lost, not the one candidate.  That is what killed 37 of
          // the first 133 finished dymc_8p5M_260905 tasks (28 %) -- i.e.
          // 0.033 aborts per 1000 Z candidates at ~9850 candidates a chunk,
          // consistent with the neighbouring bins.  The ndof == 0 bin is
          // EXACTLY EMPTY across 197 417 candidates of 20 COMPLETED tasks
          // while ndof == 1 and 2 hold 4 and 5 -- the bin is empty because
          // landing in it kills the job, not because it is forbidden.  Two
          // examples, both real: nhits=(8,1) and nhits=(2,4).  The two-track
          // fit spends ten state parameters on the common vertex, so a
          // MiniAOD leg whose stored hit list is a handful is enough; the
          // J/psi ALCARECO tracks (full RECO hits) essentially never are.
          //
          // Such a candidate carries no information for the global fit --
          // rank-0 Hessian, undefined chi2/ndof -- so it is counted and
          // dropped through the same `valid` path as the other fit failures
          // rather than being written with a degenerate factorization.
          if (ndofsigned <= 0) {
            if (stepPrints_ < stepPrintLimit_) {
              ++stepPrints_;
              std::cout << "Abort: fit has no degrees of freedom!"
                        << " icons = " << icons << " iiter = " << iiter
                        << " ndof = " << ndofsigned
                        << " nhits=(" << nhitsarr[0] << "," << nhitsarr[1] << ")"
                        << " nvalid=(" << nvalidarr[0] << "," << nvalidarr[1] << ")"
                        << " nvalidpixel=(" << nvalidpixelarr[0] << "," << nvalidpixelarr[1] << ")"
                        << " seed0(q,pt,eta)=(" << itrack->charge() << "," << itrack->pt() << "," << itrack->eta() << ")"
                        << " seed1(q,pt,eta)=(" << jtrack->charge() << "," << jtrack->pt() << "," << jtrack->eta() << ")"
                        << std::endl;
            }
            ++fitFailNdof_;
            valid = false;
            break;
          }

// std::cout << "icons = " << icons << " iiter =" << iiter << " dx = " << refftsarr[0].position() - refftsarr[1].position() << std::endl;
          
    // std::cout << "dchisqdparms.head<6>()" << std::endl;
    // std::cout << dchisqdparms.head<6>() << std::endl;
    // 
    // std::cout << "grad.head<6>()" << std::endl;
    // std::cout << grad.head<6>() << std::endl;
    // 
    // std::cout << "d2chisqdparms2.topLeftCorner<6, 6>():" << std::endl;
    // std::cout << d2chisqdparms2.topLeftCorner<6, 6>() << std::endl;
    // std::cout << "hess.topLeftCorner<6, 6>():" << std::endl;
    // std::cout << hess.topLeftCorner<6, 6>() << std::endl;
    // 
    // std::cout << "dchisqdparms.segment<6>(nparsBfield+nparsEloss)" << std::endl;
    // std::cout << dchisqdparms.segment<6>(nparsBfield+nparsEloss) << std::endl;
    // 
    // std::cout << "grad.segment<6>(nparsBfield+nparsEloss)" << std::endl;
    // std::cout << grad.segment<6>(nparsBfield+nparsEloss) << std::endl;
    // 
    // std::cout << "d2chisqdparms2.block<6, 6>(nparsBfield+nparsEloss, nparsBfield+nparsEloss):" << std::endl;
    // std::cout << d2chisqdparms2.block<6, 6>(nparsBfield+nparsEloss, nparsBfield+nparsEloss) << std::endl;
    // std::cout << "hess.block<6, 6>(nparsBfield+nparsEloss, nparsBfield+nparsEloss):" << std::endl;
    // std::cout << hess.block<6, 6>(nparsBfield+nparsEloss, nparsBfield+nparsEloss) << std::endl;
    // // 
    // 
    // std::cout << "d2chisqdparms2.block<6, 6>(trackparmidxarr[1], trackparmidxarr[1]):" << std::endl;
    // std::cout << d2chisqdparms2.block<6, 6>(trackparmidxarr[1], trackparmidxarr[1]) << std::endl;
    // std::cout << "hess.block<6, 6>(trackparmidxarr[1], trackparmidxarr[1]):" << std::endl;
    // std::cout << hess.block<6, 6>(trackparmidxarr[1], trackparmidxarr[1]) << std::endl;
    // 
    // std::cout << "d2chisqdparms2.bottomRightCorner<6, 6>():" << std::endl;
    // std::cout << d2chisqdparms2.bottomRightCorner<6, 6>() << std::endl;
    // std::cout << "hess.bottomRightCorner<6, 6>():" << std::endl;
    // std::cout << hess.bottomRightCorner<6, 6>() << std::endl;

    // const double 
    // // const double corxi0plusminus = hess(1, trackparmidxarr[1] + 1)/std::sqrt(hess(1,1)*hess(trackparmidxarr[1] + 1, trackparmidxarr[1] + 1));
    // // const double corxi1plusminus = hess(3, trackparmidxarr[1] + 3)/std::sqrt(hess(3,3)*hess(trackparmidxarr[1] + 3, trackparmidxarr[1] + 3));
    // 
    // const double cor01plus = hess(1, 3)/std::sqrt(hess(1, 1)*hess(3, 3));
    // // const double cor01minus = hess(trackparmidxarr[1] + 1, trackparmidxarr[1] + 3)/std::sqrt(hess(trackparmidxarr[1] + 1, trackparmidxarr[1] + 1)*hess(trackparmidxarr[1] + 3, trackparmidxarr[1] + 3));
    // 
    // const double cor12plus = hess(3, 5)/std::sqrt(hess(3, 3)*hess(5, 5));
    // // const double cor12minus = hess(trackparmidxarr[1] + 3, trackparmidxarr[1] + 5)/std::sqrt(hess(trackparmidxarr[1] + 3, trackparmidxarr[1] + 3)*hess(trackparmidxarr[1] + 5, trackparmidxarr[1] + 5));
    // 
    // // std::cout << "corxi0plusminus = " << corxi0plusminus << std::endl;
    // // std::cout << "corxi1plusminus = " << corxi1plusminus << std::endl;
    // std::cout << "cor01plus = " << cor01plus << std::endl;
    // // std::cout << "cor01minus = " << cor01minus << std::endl;
    // std::cout << "cor12plus = " << cor12plus << std::endl;
    // // std::cout << "cor12minus = " << cor12minus << std::endl;
          
    // std::cout << "hess(1, 1)" << std::endl;
    // std::cout << hess(1, 1) << std::endl;
    // std::cout << "hess(trackparmidxarr[1] + 1, trackparmidxarr[1] + 1)" << std::endl;
    // std::cout << hess(trackparmidxarr[1] + 1, trackparmidxarr[1] + 1) << std::endl;
    // std::cout << "hess(1, trackparmidxarr[1] + 1)" << std::endl;
    // std::cout << hess(1, trackparmidxarr[1] + 1) << std::endl;
          
          // compute final kinematics
          
          kinTree->movePointerToTheTop();
          RefCountedKinematicVertex dimu_vertex = kinTree->currentDecayVertex();
          
          if (icons == 0) {
            Jpsikin_x = dimu_vertex->position().x();
            Jpsikin_y = dimu_vertex->position().y();
            Jpsikin_z = dimu_vertex->position().z();
          }
          else {
            Jpsikincons_x = dimu_vertex->position().x();
            Jpsikincons_y = dimu_vertex->position().y();
            Jpsikincons_z = dimu_vertex->position().z(); 
          }
          
          // apply the GBL fit results to the vertex position
          const Matrix<double, 10, 1> statepca = twoTrackCart2pca(refftsarr[0], refftsarr[1]);
          const Matrix<double, 10, 1> statepcaupd = statepca + dxfull.head<10>();

// std::cout << "statepcaupd d = " << statepcaupd[6] << std::endl;

          const bool firstplus = statepcaupd[0] > 0.;

          if (icons == 0) {
            // THE SIGNED TRACK-TRACK PCA DISTANCE, raw:
            //     d = n_hat . (x_b - x_a),   n_hat = (p_a x p_b)^
            // (`twoTrackCart2pca`). Swapping the two legs flips BOTH `n_hat`
            // and `x_b - x_a`, so `d` is invariant under the leg ordering and
            // is already a well-defined signed quantity; it must NOT be
            // multiplied by the charge of leg 0, which would make its sign
            // depend on an ordering it does not depend on.
            // Under `doVtxConstraint` index 6 is frozen at zero and `d` is
            // identically zero: `Jpsi_vtxres` is then the DCA the
            // unconstrained fit would have reported.
            Jpsi_d = statepcaupd[6];
            Jpsi_x = statepcaupd[7];
            Jpsi_y = statepcaupd[8];
            Jpsi_z = statepcaupd[9];
          }
          else {
            // same convention as `Jpsi_d` above
            Jpsicons_d = statepcaupd[6];
            Jpsicons_x = statepcaupd[7];
            Jpsicons_y = statepcaupd[8];
            Jpsicons_z = statepcaupd[9];
          }

          std::array<ROOT::Math::PxPyPzMVector, 2> muarr;
          std::array<Vector3d, 2> mucurvarr;
// std::array<int, 2> muchargearr;
          
    // std::cout << dimu_vertex->position() << std::endl;
          
          // apply the GBL fit results to the muon kinematics
          for (unsigned int id = 0; id < 2; ++id) {
            const Matrix<double, 7, 1> &refFts = refftsarr[id];
            const unsigned int trackstateidx = 3*id;

            const double qbpupd = statepcaupd[trackstateidx];
            const double lamupd = statepcaupd[trackstateidx + 1];
            const double phiupd = statepcaupd[trackstateidx + 2];
            
            const double charge = std::copysign(1., qbpupd);
            const double pupd = std::abs(1./qbpupd);
            
            const double pxupd = pupd*std::cos(lamupd)*std::cos(phiupd);
            const double pyupd = pupd*std::cos(lamupd)*std::sin(phiupd);
            const double pzupd = pupd*std::sin(lamupd);
            
            muarr[id] = ROOT::Math::PxPyPzMVector(pxupd, pyupd, pzupd, trackMass[id]);
            muchargearr[id] = charge;
            
            auto &refParms = mucurvarr[id];
            refParms << qbpupd, lamupd, phiupd;
            
    // auto const &refFts = outparts[id]->currentState().freeTrajectoryState();
// auto const &refFts = refftsarr[id];
// // auto const &jac = jacarr[id];
// unsigned int trackstateidx = trackstateidxarr[id];
// 
// // JacobianCurvilinearToCartesian curv2cart(refFts.parameters());
// // const AlgebraicMatrix65& jac = curv2cart.jacobian();
// // const Matrix<double, 6, 5> jac = curv2cartJacobianAlt(refFts);
// const AlgebraicVector6 glob = refFts.parameters().vector();
// 
// const Matrix<double, 3, 1> posupd = Map<const Matrix<double, 6, 1>>(glob.Array()).head<3>() + dxfull.head<3>();
// 
// // const Matrix<double, 3, 1> momupd = Map<const Matrix<double, 6, 1>>(glob.Array()).tail<3>() + Map<const Matrix<double, 6, 5, RowMajor>>(jac.Array()).bottomLeftCorner<3, 3>()*dxfull.segment<3>(trackstateidx);
// // const Matrix<double, 3, 1> momupd = Map<const Matrix<double, 6, 1>>(glob.Array()).tail<3>() + jac.bottomLeftCorner<3, 3>()*dxfull.segment<3>(trackstateidx);
// 
// const GlobalPoint pos(posupd[0], posupd[1], posupd[2]);
// // const GlobalVector mom(momupd[0], momupd[1], momupd[2]);
// // const double charge = std::copysign(1., refFts.charge()/refFts.momentum().mag() + dxfull[trackstateidx]);
// // std::cout << "before update: reffts:" << std::endl;
// // std::cout << refFts.parameters().vector() << std::endl;
// // std::cout << "charge " << refFts.charge() << std::endl;
// // updFts = FreeTrajectoryState(pos, mom, charge, field);
// 
// 
// const CurvilinearTrajectoryParameters curv(refFts.position(), refFts.momentum(), refFts.charge());
// 
// const double qbpupd = curv.Qbp() + dxfull(trackstateidx);
// const double lamupd = curv.lambda() + dxfull(trackstateidx + 1);
// const double phiupd = curv.phi() + dxfull(trackstateidx + 2);
// 
// const double charge = std::copysign(1., qbpupd);
// const double pupd = std::abs(1./qbpupd);
// 
// const double pxupd = pupd*std::cos(lamupd)*std::cos(phiupd);
// const double pyupd = pupd*std::cos(lamupd)*std::sin(phiupd);
// const double pzupd = pupd*std::sin(lamupd);
// 
// const GlobalVector mom(pxupd, pyupd, pzupd);
// 
// muarr[id] = ROOT::Math::PxPyPzMVector(pxupd, pyupd, pzupd, mmu);
// muchargearr[id] = charge;
// 
// // std::cout << "delta eta final = " << muarr[id].eta() - refFts.momentum().eta() << std::endl;
// 
// auto &refParms = mucurvarr[id];
// // CurvilinearTrajectoryParameters curvparms(refFts.position(), refFts.momentum(), refFts.charge());
// CurvilinearTrajectoryParameters curvparms(pos, mom, charge);
// // refParms << curvparms.Qbp(), curvparms.lambda(), curvparms.phi(), curvparms.xT(), curvparms.yT();
// refParms << curvparms.Qbp(), curvparms.lambda(), curvparms.phi();
// // refParms += dxcurv;

          }
          
          // *TODO* better handling of this case?
          if ( (muchargearr[0] + muchargearr[1]) != 0) {
            std::cout << "Abort: charge flip in parameter update!"
                      << " qbp0 = " << mucurvarr[0][0]
                      << " qbp1 = " << mucurvarr[1][0]
                      << " icons = " << icons << " iiter = " << iiter
                      << " seedq0 = " << itrack->charge() << " seedq1 = " << jtrack->charge()
                      << " seedpt0 = " << itrack->pt() << " seedpt1 = " << jtrack->pt()
                      << std::endl;
            ++fitFailChargeFlip_;
            valid = false;
            break;
          }

          // charge-ordered (was hardcoded 0/1, which made the charge-
          // required gen matching fail for the ~50% of candidates whose
          // first leg is the mu-; the constrained block below always
          // ordered by charge)
          const unsigned int idxplus = muchargearr[0] > 0 ? 0 : 1;
          const unsigned int idxminus = muchargearr[0] > 0 ? 1 : 0;
          
          const ROOT::Math::PxPyPzMVector jpsitrkmom = mutrkarr[0] + mutrkarr[1];
          
          Muplustrk_pt = mutrkarr[idxplus].pt();
          Muplustrk_eta = mutrkarr[idxplus].eta();
          Muplustrk_phi = mutrkarr[idxplus].phi();
          
          Muminustrk_pt = mutrkarr[idxminus].pt();
          Muminustrk_eta = mutrkarr[idxminus].eta();
          Muminustrk_phi = mutrkarr[idxminus].phi();
          
          Jpsitrk_pt = jpsitrkmom.pt();
          Jpsitrk_eta = jpsitrkmom.eta();
          Jpsitrk_phi = jpsitrkmom.phi();
          Jpsitrk_mass = jpsitrkmom.mass();
          
          
          
          
          if (icons == 0) {
          
            Muplus_pt = muarr[idxplus].pt();
            Muplus_eta = muarr[idxplus].eta();
            Muplus_phi = muarr[idxplus].phi();
            
            Muminus_pt = muarr[idxminus].pt();
            Muminus_eta = muarr[idxminus].eta();
            Muminus_phi = muarr[idxminus].phi();
            
            Mupluskin_pt = outparts[idxplus]->currentState().globalMomentum().perp();
            Mupluskin_eta = outparts[idxplus]->currentState().globalMomentum().eta();
            Mupluskin_phi = outparts[idxplus]->currentState().globalMomentum().phi();
            
            Muminuskin_pt = outparts[idxminus]->currentState().globalMomentum().perp();
            Muminuskin_eta = outparts[idxminus]->currentState().globalMomentum().eta();
            Muminuskin_phi = outparts[idxminus]->currentState().globalMomentum().phi();
          }
          else {
            Mupluscons_pt = muarr[idxplus].pt();
            Mupluscons_eta = muarr[idxplus].eta();
            Mupluscons_phi = muarr[idxplus].phi();
            
            Muminuscons_pt = muarr[idxminus].pt();
            Muminuscons_eta = muarr[idxminus].eta();
            Muminuscons_phi = muarr[idxminus].phi();
            
            Mupluskincons_pt = outparts[idxplus]->currentState().globalMomentum().perp();
            Mupluskincons_eta = outparts[idxplus]->currentState().globalMomentum().eta();
            Mupluskincons_phi = outparts[idxplus]->currentState().globalMomentum().phi();
            
            Muminuskincons_pt = outparts[idxminus]->currentState().globalMomentum().perp();
            Muminuskincons_eta = outparts[idxminus]->currentState().globalMomentum().eta();
            Muminuskincons_phi = outparts[idxminus]->currentState().globalMomentum().phi();
          }
          
// std::cout << "Muplus pt, eta, phi = " << Muplus_pt << ", " << Muplus_eta << ", " << Muplus_phi << std::endl;
// std::cout << "Muminus pt, eta, phi = " << Muminus_pt << ", " << Muminus_eta << ", " << Muminus_phi << std::endl;
          
          
          // Sparse GBL state covariance: dense 2*Cinvd_dense^-1 with
          // Cinvd_dense=2Fs^TVinvFs equals (Fs^TVinvFs)^-1 = the sparse
          // Cinvd.solve(I) (no factor 2). Scatter free->full so the
          // topLeftCorner<6,6>/<10,10> vertex-PCA views still work.
          MatrixXd covstate = MatrixXd::Zero(nstateparms, nstateparms);
          covstate(freestateidxs, freestateidxs) =
              Cinvd.solve(MatrixXd::Identity(nstatefree, nstatefree)).eval();

          covrefmom  = covstate.topLeftCorner<6, 6>();

          
// Matrix<double, 6, 6> covrefmom;
// covrefmom = Matrix<double, 6, 6>::Zero();
//
// constexpr std::array<unsigned int, 2> localidxs = {{ 0, 3 }};
// const std::array<unsigned int, 2> globalidxs = {{ trackstateidxarr[0], trackstateidxarr[1] }};
//
// for (unsigned int iidx = 0; iidx < localidxs.size(); ++iidx) {
// for (unsigned int jidx = 0; jidx < localidxs.size(); ++jidx) {
// covrefmom.block<3, 3>(localidxs[iidx], localidxs[jidx]) = covstate.block<3, 3>(globalidxs[iidx], globalidxs[jidx]);
// }
// }
          
          if (icons == 0) {
          
            Map<Matrix<float, 3, 1>>(Muplus_refParms.data()) = mucurvarr[idxplus].cast<float>();
            Map<Matrix<float, 3, 1>>(Muminus_refParms.data()) = mucurvarr[idxminus].cast<float>();
            
    // std::cout << "nstateparms = " << nstateparms << std::endl;
    // std::cout << "dxdparms " << dxdparms.rows() << " " << dxdparms.cols() << std::endl;
            
// Muplus_jacRef.resize(3*npars);
// Map<Matrix<float, 3, Dynamic, RowMajor>>(Muplus_jacRef.data(), 3, npars) = dxdparms.block(0, trackstateidxarr[idxplus], npars, 3).transpose().cast<float>();
// 
// Muminus_jacRef.resize(3*npars);
// Map<Matrix<float, 3, Dynamic, RowMajor>>(Muminus_jacRef.data(), 3, npars) = dxdparms.block(0, trackstateidxarr[idxminus], npars, 3).transpose().cast<float>();
            
            
// MatrixXd covstate = 2.*Cinvd.solve(MatrixXd::Identity(nstateparms,nstateparms));
            
// Matrix<double, 6, 6> covrefmom;
// 
// constexpr std::array<unsigned int, 2> localidxs = {{ 0, 3 }};
// const std::array<unsigned int, 2> globalidxs = {{ trackstateidxarr[0], trackstateidxarr[1] }};
// 
// for (unsigned int iidx = 0; iidx < localidxs.size(); ++iidx) {
// for (unsigned int jidx = 0; jidx < localidxs.size(); ++jidx) {
// covrefmom.block<3, 3>(localidxs[iidx], localidxs[jidx]) = covstate.block<3, 3>(globalidxs[iidx], globalidxs[jidx]);
// } 
// }
            
            
            
            const Matrix<double, 1, 6> mjacalt =
                massJacobianAltD(refftsarr[0], refftsarr[1], massForConstraintHelpers);


            Jpsi_sigmamass = std::sqrt((mjacalt*covrefmom*mjacalt.transpose())[0]);

            // ---- the covariance itself, permuted into (plus, minus) ------
            {
              const std::array<unsigned int, 6> perm = {{3 * idxplus, 3 * idxplus + 1, 3 * idxplus + 2,
                                                         3 * idxminus, 3 * idxminus + 1, 3 * idxminus + 2}};
              Jpsi_covrefmom.clear();
              Jpsi_covrefmom.reserve(21);
              for (unsigned int i = 0; i < 6; ++i) {
                for (unsigned int j = i; j < 6; ++j) {
                  Jpsi_covrefmom.push_back(float(covrefmom(perm[i], perm[j])));
                }
              }
              Jpsi_jacrefmom.assign(6, 0.f);
              for (unsigned int i = 0; i < 6; ++i) {
                Jpsi_jacrefmom[i] = float(mjacalt(0, perm[i]));
              }

              // q/p at the reference: state[6] is the charge, segment<3>(3)
              // the momentum -- the same expression massJacobianAltD uses.
              const double qopp = refftsarr[idxplus][6] / refftsarr[idxplus].segment<3>(3).norm();
              const double qopm = refftsarr[idxminus][6] / refftsarr[idxminus].segment<3>(3).norm();
              Jpsi_qoprefplus = float(qopp);
              Jpsi_qoprefminus = float(qopm);

              const double cpp = covrefmom(perm[0], perm[0]);
              const double cmm = covrefmom(perm[3], perm[3]);
              const double cpm = covrefmom(perm[0], perm[3]);
              const double srp = (std::abs(qopp) > 0.) ? std::sqrt(std::max(cpp, 0.)) / std::abs(qopp) : 0.;
              const double srm = (std::abs(qopm) > 0.) ? std::sqrt(std::max(cmm, 0.)) / std::abs(qopm) : 0.;
              Jpsi_sigmarelplus = float(srp);
              Jpsi_sigmarelminus = float(srm);
              // d ln p = -d(q/p)/(q/p), so the momentum correlation carries
              // sign(q+ q-) relative to the q/p one -- i.e. it flips for an
              // opposite-sign pair, which is every candidate here. Written
              // out rather than assumed, so a same-sign control sample is
              // still right.
              Jpsi_rhomom = (srp > 0. && srm > 0. && qopp != 0. && qopm != 0.)
                                ? float(cpm / (qopp * qopm) / (srp * srm))
                                : 0.f;
              // The ANGULAR share of the mass variance: everything the two
              // curvatures do not carry. `1.5 - f_ang` is the Jensen
              // coefficient the spec asks for.
              Matrix<double, 1, 6> jkappa = Matrix<double, 1, 6>::Zero();
              jkappa(0, perm[0]) = mjacalt(0, perm[0]);
              jkappa(0, perm[3]) = mjacalt(0, perm[3]);
              const double sm2 = double(Jpsi_sigmamass) * double(Jpsi_sigmamass);
              const double skap2 = (jkappa * covrefmom * jkappa.transpose())[0];
              Jpsi_fang = (sm2 > 0.) ? float(1. - skap2 / sm2) : 0.f;
            }

// std::cout << "covrefmom" << std::endl;
// std::cout << covrefmom << std::endl;
// std::cout << "Jpsi_sigmamass = " << Jpsi_sigmamass << std::endl;

            // Mass-projected influence export (per-candidate mass-CF
            // ingredients, doRes port). The candidate-mass error responds
            // linearly to the noise vector: delta m = w^T n with
            // w = (C a)^T F^T Vinv, a = mass Jacobian embedded at the
            // joint-state momentum block (entries 0..5). Per registered
            // material block b: resinfv = signed dof weights in
            // dV^{1/2}-standardized units (5 floats, sign carries the
            // Landau skew), resinfvarv = |u_b|^2 = the block's variance
            // contribution to sigma_m^2, resinfcov = their sum. The
            // GAUSSIAN remainder (hits + beamspot + pointing) is
            // sigma_m^2 - resinfcov by construction -- unlike the
            // single-track tree, resinfcov here does NOT include the hit
            // share. reseigidx labels the entries (MS/ioni global params,
            // matching msmoliidx/ioniurbanidx records).
            resinfv.clear();
            resinfvarv.clear();
            resinfcov = 0.;
            resinfcovhit = 0.f;
            resinfcovgrp = 0.f;
            reshitidx.clear();
            reshitcls.clear();
            cfhitclsv.clear();
            cfhitvv.clear();
            reseigidx.clear();
            reseigv.clear();
            resinfbv.clear();
            cfmsv.clear();
            cfdelv.clear();
            cfiorev.clear();
            cfioimv.clear();
            cfradrev.clear();
            cfradimv.clear();
            cfvgf = 0.f;
            cfok = false;
            cfnblock = 0;
            cfnpooled = 0;
            cfgrpv.clear();
            cfgrpmsv.clear();
            cfgrpdelv.clear();
            cfgrpiorev.clear();
            cfgrpioimv.clear();
            cfgrpradrev.clear();
            cfgrpradimv.clear();
            cfgrpclosure = 0.f;
            if (dores && !dVs.empty()) {
              // ================= THE sqrt(dV_b) CACHE ======================
              // `dV_b^{1/2}` is a property of the BLOCK, not of the
              // functional: the same matrix serves the MASS influence, the
              // VERTEX influence and the two BEAM influences.  It used to be
              // recomputed from scratch in each of those loops -- four
              // `SelfAdjointEigenSolver` calls per block per candidate, which
              // is where the beam functionals' +63 % of maker time went.  It
              // is now computed ONCE here and read everywhere.
              //
              // The ionization SIGN rule needs the dominant eigenVECTOR of
              // the same decomposition (sign-fixed by `udir(0) >= 0`), so
              // that is cached alongside, for family 11 only; every other
              // family leaves an empty vector and the consumers fall back to
              // sign +1 exactly as before.
              //
              // Bit-identity is not an approximation here: it is the same
              // solver on the same matrix, so every consumer sees the same
              // bits it computed for itself before.
              std::vector<MatrixXd> ressqrtdV(dVs.size());
              std::vector<VectorXd> resionidir(dVs.size());
              for (unsigned int ires = 0; ires < dVs.size(); ++ires) {
                const unsigned int r0c = resblockrng[ires][0];
                const unsigned int nbc = resblockrng[ires][1];
                const MatrixXd dVb = MatrixXd(dVs[ires]).block(r0c, r0c, nbc, nbc);
                const SelfAdjointEigenSolver<MatrixXd> eigv(dVb);
                ressqrtdV[ires] = eigv.eigenvectors() *
                                  eigv.eigenvalues().cwiseMax(0.).cwiseSqrt().asDiagonal() *
                                  eigv.eigenvectors().transpose();
                const int famc = ires < resfamily_.size() ? resfamily_[ires] : -1;
                if (famc == 11 && nbc > 0) {
                  int imax = 0;
                  for (int q = 1; q < int(nbc); ++q) {
                    if (std::abs(eigv.eigenvalues()(q)) > std::abs(eigv.eigenvalues()(imax))) {
                      imax = q;
                    }
                  }
                  VectorXd udirc = eigv.eigenvectors().col(imax);
                  if (udirc(0) < 0.) {
                    udirc = -udirc;
                  }
                  resionidir[ires] = udirc;
                }
              }
              // ============================================================
              VectorXd afull = VectorXd::Zero(nstateparms);
              afull.head<6>() = mjacalt.transpose();
              VectorXd afree = VectorXd::Zero(nstatefree);
              for (unsigned int i = 0; i < nstatefree; ++i) {
                afree(i) = afull(freestateidxs[i]);
              }
              const VectorXd wmass = VinvF*Cinvd.solve(afree);
              // the BEAM block's share of sigma_m^2, kept out of `resinfcov`
              // (material) and out of `resinfcovhit` so `cfvgf` keeps its
              // meaning: the Gaussian remainder = hits + beamspot + pointing.
              double vbsmass = 0., vbsmassx = 0., vbsmassy = 0.;
              for (unsigned int ires = 0; ires < dVs.size(); ++ires) {
                const unsigned int r0 = resblockrng[ires][0];
                const unsigned int nb = resblockrng[ires][1];
                const VectorXd ub = ressqrtdV[ires] * wmass.segment(r0, nb);
                reseigidx.push_back(resglobidx[ires]);
                reshitidx.push_back(ires < resvalidhit_.size() ? resvalidhit_[ires] : -1);
                reshitcls.push_back(static_cast<short>(ires < rescls_.size() ? rescls_[ires] : -1));
                for (unsigned int j = 0; j < 5; ++j) {
                  resinfv.push_back(j < nb ? ub(j) : 0.f);
                }
                const double vb = ub.squaredNorm();
                resinfvarv.push_back(vb);
                // MATERIAL and HIT shares are accumulated SEPARATELY.
                // `resinfcov` is the MATERIAL share alone, so
                // `cfmass_vgf = (sigma_m^2 - resinfcov)/sigma_m^2` is the
                // TOTAL Gaussian share (hits + beamspot + pointing), which is
                // what the self-consistent-sigma correction's
                // `a_i = (1 + f_hit) sigma_i/m_i` needs. Folding the hit
                // blocks into `resinfcov` would change both.
                const int fam = ires < resfamily_.size() ? resfamily_[ires] : -1;
                if (fam == 8 || fam == 9) {
                  resinfcovhit += float(vb);
                } else if (fam == kBeamSpotFamily) {
                  // the luminous region is a Gaussian block like a hit, but
                  // it is neither a hit nor material: its own accumulator.
                  vbsmass += vb;
                  // THE BEAM BLOCK'S VARIANCE SHARE, SPLIT BY DIRECTION.
                  // The luminous-region WIDTHS are physical parameters, so
                  // they float like a hit class: `sigma_x -> sqrt(k_x)
                  // sigma_x` sends `covBS -> D covBS D` with
                  // `D = diag(sqrt(k_x), sqrt(k_y), 1)`, and
                  //   v_b(k) = sum_ij w_i w_j C_ij d_i d_j ,  d = (rt kx, rt ky, 1)
                  // whose derivative at k = 1 is, EXACTLY,
                  //   dv/dk_x = w_x^2 C_xx + w_x w_y C_xy + w_x w_z C_xz
                  // -- each cross term split in half between the two
                  // directions.  So the half-split shares below ARE the
                  // first derivatives, they sum to `v_b` identically, and a
                  // LINEAR variance scale `(1 + eps)` on them is the same
                  // parameterisation the hit classes use (`--hit-mode
                  // linear`, where the card value IS `eps`).
                  {
                    const Matrix<double, 3, 1> wbb = wmass.segment<3>(r0);
                    const double cxx = bscovBS(0, 0), cyy = bscovBS(1, 1), czz = bscovBS(2, 2);
                    const double cxy = bscovBS(0, 1), cxz = bscovBS(0, 2), cyz = bscovBS(1, 2);
                    vbsmassx = wbb(0)*wbb(0)*cxx + wbb(0)*wbb(1)*cxy + wbb(0)*wbb(2)*cxz;
                    vbsmassy = wbb(1)*wbb(1)*cyy + wbb(0)*wbb(1)*cxy + wbb(1)*wbb(2)*cyz;
                    (void)czz;
                  }
                } else if (fam == 15) {
                  // The parmtype-15 blocks are a RE-PARTITION of the
                  // parmtype-10/11 noise (`sum_g dQ_g == dQMS + dQI`), not an
                  // addition to it, so they must NOT enter `resinfcov`: that
                  // would double-count the material share and, here, silently
                  // move `cfmass_vgf = (sigma_m^2 - resinfcov)/sigma_m^2`.
                  // Same split as the single-track maker.
                  resinfcovgrp += float(vb);
                } else {
                  resinfcov += vb;
                }
              }

              {
                const double sm2b = double(Jpsi_sigmamass) * double(Jpsi_sigmamass);
                Jpsi_massvbs = sm2b > 0. ? float(vbsmass / sm2b) : 0.f;
                Jpsi_massvbsx = sm2b > 0. ? float(vbsmassx / sm2b) : 0.f;
                Jpsi_massvbsy = sm2b > 0. ? float(vbsmassy / sm2b) : 0.f;
              }

              // Per-hit-class Gaussian shares, ascending in class index.
              {
                const double sm2c = double(Jpsi_sigmamass) * double(Jpsi_sigmamass);
                std::array<double, kNHitResClasses> vcls{};
                bool anycls = false;
                for (std::size_t ires = 0; ires < reshitcls.size(); ++ires) {
                  const int c = reshitcls[ires];
                  if (c < 0 || c >= kNHitResClasses) {
                    continue;
                  }
                  vcls[c] += resinfvarv[ires];
                  anycls = true;
                }
                if (anycls && sm2c > 0.) {
                  for (int c = 0; c < kNHitResClasses; ++c) {
                    if (vcls[c] == 0.) {
                      continue;
                    }
                    cfhitclsv.push_back(static_cast<short>(c));
                    cfhitvv.push_back(float(vcls[c] / sm2c));
                  }
                }
              }

              // ---- THE FUNCTIONALS THAT SHARE ONE EVALUATOR PASS --------
              //
              // The candidate MASS, the vertex DCA and the two whitened
              // BEAM-LINE pulls are four linear functionals of the SAME
              // converged fit: the same blocks, the same step records, and
              // different per-block weights.  Every `cvhcf` exponent
              // primitive depends on (weight, tau) ONLY through the product
              // `w tau`, so the four are ONE `cvhcf::trackExponents` pass on
              // the concatenated argument list { w_{b,k} tau_j } instead of
              // four passes over the same records -- and the pooling by
              // global index, the Moliere step parameters, the `gshape_elec`
              // rows and the radiative spectra are then built once rather
              // than four times.  It is EXACT: `phi_{aU}(tau) = phi_U(a tau)`
              // is an identity and the products are formed by the same
              // expression a single call forms them with, so each
              // functional's arrays are bitwise what its own call wrote.
              //
              // Each site below REGISTERS its weights here, in the order the
              // exported branches expect, and reads its results back from the
              // single pass that follows the beam block; nothing else about
              // them changes.  The weights are COPIED because `vabs` and the
              // two `vabsbs` live only inside their own blocks.
              std::vector<cvhcf::TrackInput> cfins;
              std::vector<std::vector<float>> cfvarw;
              std::vector<std::vector<float>> cfsgnw;
              std::vector<cvhcf::TrackResult> cfress;
              cfins.reserve(4);
              cfvarw.reserve(4);
              cfsgnw.reserve(4);
              int cfslotmass = -1;
              int cfslotvtx = -1;
              std::array<int, 2> cfslotbs = {{-1, -1}};
              auto cfRegister = [&](const float *var, int nvar, const float *sgn, double sigma,
                                    double ioniSign) -> int {
                cfvarw.emplace_back(var, var + nvar);
                cfsgnw.emplace_back();
                if (sgn != nullptr) {
                  cfsgnw.back().assign(sgn, sgn + nvar);
                }
                cvhcf::TrackInput ci;
                ci.resglobidx = resglobidx.data();
                ci.resfamily = resfamily_.data();
                ci.resvarv = cfvarw.back().data();
                ci.nres = int(std::min({resglobidx.size(), resfamily_.size(), cfvarw.back().size()}));
                ci.ressgn = cfsgnw.back().empty() ? nullptr : cfsgnw.back().data();
                ci.ms = {msmoliidx.data(), msmoliv.data(), int(msmoliidx.size()),
                         msmoliidx.empty() ? 0 : int(msmoliv.size() / msmoliidx.size())};
                ci.ioni = {ioniurbanidx.data(), ioniurbanv.data(), int(ioniurbanidx.size()),
                           ioniurbanidx.empty() ? 0 : int(ioniurbanv.size() / ioniurbanidx.size())};
                ci.qsc = {ioniqscaleidx.data(), ioniqscalev.data(), int(ioniqscaleidx.size()), 2};
                ci.rad = {radstepidx.data(), radstepv.data(), int(radstepidx.size()), RADSTEP_STRIDE};
                // material-group column of each record; see the single-track
                // maker for why these are set explicitly
                ci.ms.groupCol = ci.ms.stride >= 10 ? 9 : -1;
                ci.ioni.groupCol = ci.ioni.stride - 1;
                ci.rad.groupCol = RADSTEP_STRIDE - 1;
                ci.radspec = radstepspecv.data();
                ci.radvgrid = radvgrid.data();
                ci.radnv = int(radvgrid.size());
                ci.sigma = sigma;
                ci.ioniSign = ioniSign;
                ci.wantDelta = true;
                ci.wantGroups = exportCfGroupExponents_;
                // No functional's reference model splits the DELTA-RAY family
                // per group (`cf_mass_likelihood.build_pairs_tt` has no
                // `Sdel`), so the flat delta family is exported for comparison
                // but is not split.
                ci.wantGroupDelta = false;
                cfins.push_back(ci);
                return int(cfins.size()) - 1;
              };

              // ---- THE RESOLUTION-CF EXPONENTS, for the MASS functional ---
              //
              // Same blocks, same step records, DIFFERENT functional: the
              // standardization is sigma_m rather than sqrt(refCov(0,0)), and
              // the ionization (and radiative) weight carries the sign -1 for
              // BOTH legs and BOTH charges. That sign is not a convention:
              // dm = (dm/d(q/p)) q cs dE = -p^2 cs (dm/dp) dE < 0, i.e. q^2
              // removes the charge and an energy loss on EITHER muon can only
              // LOWER the pair mass, so the two legs' Landau skews ADD.
              // (cf_mass_likelihood.IONI_SGN; using sign(sum u_b) instead --
              // an arbitrary noise-eigenvector convention, +1 on 50.5 % of
              // candidates -- cancelled the skew and moved the unbinned scale
              // by 0.1e-3.)
              //
              // The DELTA-RAY family is computed but is not part of the
              // reference mass model (`build_pairs_tt` has no `Sdel`); it is
              // exported so the two can be compared, and the reader leaves it
              // out of the pairs cache unless asked.
              if (exportCfExponents_) {
                cfslotmass = cfRegister(resinfvarv.data(), int(resinfvarv.size()), nullptr, Jpsi_sigmamass, -1.);
                // The GAUSSIAN REMAINDER: hits + beamspot + pointing, i.e.
                // sigma_m^2 minus the material share, BY CONSTRUCTION --
                // unlike the single-track tree, resinfcov here does not carry
                // the hit blocks at all (they are not registered).
                const double sm2 = double(Jpsi_sigmamass) * double(Jpsi_sigmamass);
                cfvgf = (sm2 > 0.) ? float((sm2 - resinfcov) / sm2) : 0.f;
              }

              // ================= THE VERTEX-CONSTRAINT RESIDUAL ===========
              //
              // State index 6 is the SIGNED track-track PCA distance
              //     theta_6 = n_hat . (x_b - x_a),  n_hat = (p_a x p_b)^
              // (`twoTrackCart2pca`). Two free helices have 10 vertex
              // parameters, two through a common point have 9: a common-
              // vertex constraint removes exactly ONE dof and leaves exactly
              // ONE residual per candidate, the DCA the unconstrained fit
              // finds. Both muons come from one gen point, so on ideal
              // geometry its mean is ZERO BY CONSTRUCTION -- this is the
              // mass term with a delta kernel at zero: no kernel, no theory,
              // no PDG input, a pure resolution term.
              //
              // Both configurations of index 6 are supported, and the same
              // three numbers come out of each:
              //   * index 6 FREE (`doVtxConstraint == False`): the fit itself
              //     reports the DCA, so
              //       sigma_v^2 = C_66,  w_v = Vinv F C e_6,  r_v = theta_6
              //   * index 6 FROZEN: b is zero on every free index at
              //     convergence and nonzero only on index 6, and with
              //     F_6 = Ffull.col(6), h_f6 = F_f^T Vinv F_6, Cs = C h_f6,
              //       sigma_v^2 = 1/(h_66 - h_f6^T Cs)
              //       b_6 = -(Vinv F_6).r          (minus the half-gradient)
              //       r_v = theta_6^frozen + sigma_v^2 b_6   (one Newton step)
              //       w_v = sigma_v^2 (Vinv F_6 - Vinv F Cs)
              // In both, `sum_b |dV_b^{1/2} w_v,b|^2 == sigma_v^2` exactly
              // (one line: w_v^T V w_v), which is the closure gate
              // `Jpsi_vtxvchk`.
              if (exportVtxResidual_) {
                const VectorXd F6 = Ffull.col(6);
                const VectorXd VinvF6 = Vinvsparse * F6;
                // THE POST-STEP RESIDUAL.  `rfull` is the residual at the
                // LINEARISATION point, and in a GBL-type fit the reference
                // trajectory is built BY PROPAGATING, so every process-noise
                // row of `rfull` is identically zero there.  The vertex state
                // enters ONLY through the ihit == 0 propagation rows, so
                // `F_6^T Vinv rfull` vanishes identically on every candidate
                // and says nothing.  The gradient that matters is the one at
                // the OPTIMUM,
                //     rho = rfull + F_f dxfree
                // the same vector the exported `Rr` is built from: there the
                // free-index gradient is zero and the index-6 one is not.
                const VectorXd rho = rfull + Fsparse * dxfree;
                const double g6 = VinvF6.dot(rho);   // the half-gradient
                // THE CONVERGENCE GATE, DIMENSIONLESS.  `|g_i| sqrt(C_ii)` is
                // the remaining Newton step of free parameter i in units of
                // that parameter's own error, so it is directly comparable
                // with `|z_v| = |b_6| sigma_v`; the raw half-gradient is in
                // 1/(cm, rad, GeV^-1) units and says nothing on its own.
                // `r_v` below is built from the SAME pair (reference + one
                // step), so the two configurations of index 6 are compared at
                // the same linearisation point.
                const VectorXd gfree = VinvF.transpose() * rho;
                double bfree = 0.;
                for (unsigned int i = 0; i < nstatefree; ++i) {
                  const Eigen::Index is = freestateidxs[i];
                  const double ci = std::sqrt(std::max(covstate(is, is), 0.));
                  bfree = std::max(bfree, std::abs(gfree(i)) * ci);
                }
                // the fitted vertex covariance, packed upper-triangular --
                // what the rows-OFF run hands the leave-one-out construction
                {
                  const Matrix<double, 3, 3> Cv = covstate.block<3, 3>(7, 7);
                  Jpsi_covvtx = {{float(Cv(0, 0)), float(Cv(0, 1)), float(Cv(0, 2)),
                                  float(Cv(1, 1)), float(Cv(1, 2)), float(Cv(2, 2))}};
                }
                Jpsi_vtxbfree = float(bfree);
                Jpsi_vtxb6 = float(-g6);
                Jpsi_vtxfree = !doVtxConstraint_;

                // NO CHARGE RE-SIGN. theta_6 is INVARIANT under swapping the
                // two legs (`twoTrackCart2pca` flips both `n_hat` and
                // `x_b - x_a`), so it is already a well-defined signed DCA and
                // multiplying it by the charge of a leg would make its sign
                // depend on the leg ordering. The CF exponents below are built
                // from `w_v` in this same convention, so re-signing the
                // residual and not the weights would put the Landau skew on
                // the wrong side. `Jpsi_vtxfirstplus` records which leg is the
                // positive one.
                Jpsi_vtxfirstplus = firstplus;
                double sig2 = 0.;
                double rv = 0.;
                VectorXd wv;
                // THE UNCONSTRAINED MASS.  Freezing theta_6 at zero is a
                // CONDITIONING of the unconstrained solution,
                //     x_c = x_u - C_u e_6 sigma_v^-2 r_v ,
                // so the mass functional a = dm/dx transforms as
                //     m_u = m_c + cov(m, theta_6) sigma_v^-2 r_v ,
                //     cov(m, theta_6) = a^T C_u e_6 .
                // `a` lives on the momentum block alone (`afull.head<6>()`),
                // so its index-6 entry is zero and, with the frozen-index
                // block inverse C_u,f6 = -C h_f6 sigma_v^2,
                //     cov(m, theta_6) = -sigma_v^2 (C a_f).h_f6
                // and sigma_v^2 cancels out of `m_u` itself.  `dmdv` below is
                // that slope; with index 6 free the fit already reports the
                // unconstrained mass and the slope is zero, while the
                // covariance element is the plain a_f^T C e_6.
                double dmdv = 0.;
                double covmv = 0.;
                if (Jpsi_vtxfree) {
                  VectorXd afree6 = VectorXd::Zero(nstatefree);
                  for (unsigned int i = 0; i < nstatefree; ++i) {
                    if (freestateidxs[i] == 6) {
                      afree6(i) = 1.;
                    }
                  }
                  const VectorXd Ce6 = Cinvd.solve(afree6);
                  wv = VinvF * Ce6;
                  sig2 = covstate(6, 6);
                  rv = statepcaupd[6];
                  covmv = afree.dot(Ce6);
                } else {
                  const VectorXd hf6 = VinvF.transpose() * F6;
                  const VectorXd Cs = Cinvd.solve(hf6);
                  const double denom = F6.dot(VinvF6) - hf6.dot(Cs);
                  sig2 = denom > 0. ? 1. / denom : 0.;
                  wv = sig2 * (VinvF6 - VinvF * Cs);
                  rv = statepcaupd[6] - sig2 * g6;
                  dmdv = -Cinvd.solve(afree).dot(hf6);
                  covmv = sig2 * dmdv;
                }
                wvtxinf = wv;
                Jpsi_vtxres = float(rv);
                Jpsi_vtxsig = float(std::sqrt(std::max(sig2, 0.)));
                Jpsi_vtxz = Jpsi_vtxsig > 0.f ? Jpsi_vtxres / Jpsi_vtxsig : 0.f;
                Jpsi_vtxdchi2 = Jpsi_vtxz * Jpsi_vtxz;
                Jpsi_covmassvtx = float(covmv);
                Jpsi_mass_unc = float((muarr[0] + muarr[1]).M() + dmdv * rv);
                // ---- the per-block influence, variance shares and signs ---
                //
                // THE IONIZATION SIGN. The ionization CF is not even in its
                // weight (an energy loss goes one way), so a functional whose
                // influence changes sign from block to block needs the sign
                // to travel with the block -- `cvhcf::TrackInput::ressgn`.
                // `dQI` is rank one: in the curvilinear frame it is
                // `e_0 e_0^T sigma^2` and the local rotation spreads it over
                // all five rows, so the signed coefficient is `w_b . u` with
                // `u` the leading eigenvector of `dQI` oriented by its qop
                // component, and NOT the single row `w[r0]`, which the
                // rotation leaves no longer aligned with q/p. The leg CHARGE
                // multiplies it: the two legs have opposite charges, which is
                // exactly what the MASS functional's single `ioniSign = -1`
                // hides. `Jpsi_vtxsgnchk` applies the identical rule to the
                // MASS influence `wmass`, where the answer must be -1 on every
                // ionization block, so the convention is checked per candidate
                // rather than asserted.
                const unsigned int nb_all = dVs.size();
                resinfvtxv.clear();
                resinfvtxv.reserve(5 * nb_all);
                vtxvarv.clear();
                vtxvarv.reserve(nb_all);
                vtxsgnv.clear();
                vtxsgnv.reserve(nb_all);
                std::vector<float> vabs(nb_all, 0.f);
                double vsum = 0., vmat = 0., vhit = 0., vms = 0., vioni = 0., vbs = 0.;
                double vbsx = 0., vbsy = 0.;
                double mms = 0., mioni = 0.;
                double nioni = 0., nionimass = 0.;
                for (unsigned int ires = 0; ires < nb_all; ++ires) {
                  const unsigned int r0 = resblockrng[ires][0];
                  const unsigned int nb = resblockrng[ires][1];
                  const VectorXd wb = wv.segment(r0, nb);
                  const VectorXd ab = ressqrtdV[ires] * wb;
                  for (unsigned int j = 0; j < 5; ++j) {
                    resinfvtxv.push_back(j < nb ? float(ab(j)) : 0.f);
                  }
                  const double vb = ab.squaredNorm();
                  vabs[ires] = float(vb);
                  const int fam = ires < resfamily_.size() ? resfamily_[ires] : -1;
                  if (fam != 15) {
                    vsum += vb;
                  }
                  if (fam == 10 || fam == 11) {
                    vmat += vb;
                    (fam == 10 ? vms : vioni) += vb;
                    (fam == 10 ? mms : mioni) += double(resinfvarv[ires]);
                  } else if (fam == 8 || fam == 9) {
                    vhit += vb;
                  } else if (fam == kBeamSpotFamily) {
                    vbs += vb;
                    {
                      const Matrix<double, 3, 1> wbb = wv.segment<3>(r0);
                      const double cxx = bscovBS(0, 0), cyy = bscovBS(1, 1), czz = bscovBS(2, 2);
                      const double cxy = bscovBS(0, 1), cxz = bscovBS(0, 2), cyz = bscovBS(1, 2);
                      vbsx = wbb(0)*wbb(0)*cxx + wbb(0)*wbb(1)*cxy + wbb(0)*wbb(2)*cxz;
                      vbsy = wbb(1)*wbb(1)*cyy + wbb(0)*wbb(1)*cxy + wbb(1)*wbb(2)*cyz;
                      (void)czz;
                    }
                  }
                  float sg = 1.f;
                  if (fam == 11 && nb > 0 && resionidir[ires].size() == Eigen::Index(nb)) {
                    const VectorXd &udir = resionidir[ires];
                    const unsigned int leg = (ires >= resLegStart_[1]) ? 1u : 0u;
                    const double qleg = refftsarr[leg][6] >= 0. ? 1. : -1.;
                    sg = float(qleg * (wb.dot(udir) > 0. ? -1. : 1.));
                    // The same rule on the MASS influence, whose answer is
                    // known to be -1.  VARIANCE-WEIGHTED, because `cvhcf`
                    // pools blocks by global index and takes
                    // `sign(sum_i s_i v_i)`: a block carrying no variance
                    // cannot flip a pooled sign, and its own `w_b . u` is a
                    // ratio of two small numbers.
                    const VectorXd wmb = wmass.segment(r0, nb);
                    const double vm = double(resinfvarv[ires]);
                    nioni += vm;
                    if (qleg * (wmb.dot(udir) > 0. ? -1. : 1.) < 0.) {
                      nionimass += vm;
                    }
                  }
                  vtxsgnv.push_back(sg);
                }
                for (unsigned int ires = 0; ires < nb_all; ++ires) {
                  vtxvarv.push_back(sig2 > 0. ? float(vabs[ires] / sig2) : 0.f);
                }
                Jpsi_vtxvchk = sig2 > 0. ? float(std::abs(vsum / sig2 - 1.)) : 0.f;
                Jpsi_vtxvgf = sig2 > 0. ? float((sig2 - vmat) / sig2) : 0.f;
                Jpsi_vtxvhit = sig2 > 0. ? float(vhit / sig2) : 0.f;
                Jpsi_vtxvbs = sig2 > 0. ? float(vbs / sig2) : 0.f;
                Jpsi_vtxvbsx = sig2 > 0. ? float(vbsx / sig2) : 0.f;
                Jpsi_vtxvbsy = sig2 > 0. ? float(vbsy / sig2) : 0.f;
                Jpsi_vtxvms = sig2 > 0. ? float(vms / sig2) : 0.f;
                Jpsi_vtxvioni = sig2 > 0. ? float(vioni / sig2) : 0.f;
                const double sm2v = double(Jpsi_sigmamass) * double(Jpsi_sigmamass);
                Jpsi_massvms = sm2v > 0. ? float(mms / sm2v) : 0.f;
                Jpsi_massvioni = sm2v > 0. ? float(mioni / sm2v) : 0.f;
                Jpsi_vtxsgnchk = nioni > 0. ? float(nionimass / nioni) : 0.f;

                // per-hit-class Gaussian shares of sigma_v^2
                vtxhitclsv.clear();
                vtxhitvv.clear();
                {
                  std::array<double, kNHitResClasses> vcls{};
                  bool anycls = false;
                  for (std::size_t ires = 0; ires < reshitcls.size() && ires < nb_all; ++ires) {
                    const int c = reshitcls[ires];
                    if (c < 0 || c >= kNHitResClasses) {
                      continue;
                    }
                    vcls[c] += vabs[ires];
                    anycls = true;
                  }
                  if (anycls && sig2 > 0.) {
                    for (int c = 0; c < kNHitResClasses; ++c) {
                      if (vcls[c] == 0.) {
                        continue;
                      }
                      vtxhitclsv.push_back(static_cast<short>(c));
                      vtxhitvv.push_back(float(vcls[c] / sig2));
                    }
                  }
                }

                // ---- the CF exponents AT THE VERTEX WEIGHTS ---------------
                // Same blocks, same step records, third functional: the
                // standardization is sigma_v and the ionization sign travels
                // per block, so `ioniSign` is 1 and `ressgn` carries
                // everything.
                cfvtxmsv.clear();
                cfvtxdelv.clear();
                cfvtxiorev.clear();
                cfvtxioimv.clear();
                cfvtxradrev.clear();
                cfvtxradimv.clear();
                cfvtxgrpv.clear();
                cfvtxgrpmsv.clear();
                cfvtxgrpdelv.clear();
                cfvtxgrpiorev.clear();
                cfvtxgrpioimv.clear();
                cfvtxgrpradrev.clear();
                cfvtxgrpradimv.clear();
                cfvtxgrpvqmsv.clear();
                cfvtxgrpvqiov.clear();
                cfvtxgrpclosure = 0.f;
                Jpsi_vtxok = false;
                if (exportCfExponents_ && Jpsi_vtxsig > 0.f) {
                  cfslotvtx = cfRegister(vabs.data(), int(vabs.size()), vtxsgnv.data(), Jpsi_vtxsig, 1.);
                }
              }

              // ============= THE BEAM-LINE (LUMINOUS-REGION) RESIDUALS ====
              //
              // With the beam rows ON they are MEASUREMENT rows, so the
              // unbiased residual is the LEAVE-ONE-OUT (innovation) one: the
              // vertex the fit finds WITHOUT the beam rows, minus the beam
              // line at that vertex's z.  With
              //     A    = C[7:10, 7:10]           the fitted vertex covariance
              //     M    = covBS - A = Cov(rho_B)  the FITTED residual's covariance
              //     rho_B = x_v^fit - b0           the post-step residual on those rows
              // Sherman-Morrison gives, exactly,
              //     e_B      = covBS M^-1 rho_B            = x_v^{-B} - b0
              //     Cov(e_B) = covBS M^-1 covBS            = C_{-B} + covBS
              // and the CONDITIONAL-MEAN projector to the transverse plane
              //     P = [ I_2 , -s ] ,   s = covBS[0:2,2] / covBS[2,2]
              // is exactly "the beam line evaluated at that vertex's z",
              // because `P covBS P^T` IS the conditional covariance
              // Sigma_{xy|z}.  So the task's two readings coincide:
              //     Cov(r_bs) = P Cov(e_B) P^T = P C_{-B} P^T + Sigma_{xy|z}.
              // z is left out: its row is weightless (sigmaZ ~ 3.5 cm
              // against a ~100 um vertex error), and conditioning on it is
              // what makes the transverse pair well defined.
              //
              // The influence weights come from `rho = (I - F C F^T Vinv) n`:
              // with `G = P covBS M^-1` and `u_i` the vector that is `G_i^T`
              // on the three beam rows and zero elsewhere,
              //     w_i = u_i - Vinv F C (F^T u_i)
              // and `w_i^T V w_i = G_i M G_i^T = Cov(r_bs)_{ii}` identically
              // (the algebra: `(I - Q Vinv) V (I - Vinv Q) = V - Q` because
              // `Q Vinv Q = Q`).  `Jpsi_bsvchk` is that closure measured
              // block by block.
              //
              // A FREE INTERNAL GATE: `w_i` restricted to the beam rows is
              //     (I - covBS^-1 A) M^-1 covBS P_i^T = P_i^T
              // exactly, so `Jpsi_bsmeanbs` (= -w_i on those rows) MUST equal
              // -P to machine precision.  It is written out rather than
              // asserted.
              if (exportBsResidual_ && bsConstraint_ && bsrow >= 0) {
                Jpsi_bsres = {{0.f, 0.f}};
                Jpsi_bscov = {{0.f, 0.f, 0.f}};
                Jpsi_bsz = {{0.f, 0.f}};
                Jpsi_bschi2 = 0.f;
                Jpsi_bschi2fit = 0.f;
                Jpsi_bschi20 = float(bschisq0);
                Jpsi_bsvchk = 0.f;
                Jpsi_bsok = false;
                Jpsi_bscovlo = {{0.f, 0.f, 0.f, 0.f, 0.f, 0.f}};
                Jpsi_bslinv = {{0.f, 0.f, 0.f}};
                Jpsi_bsmeig = 0.f;
                Jpsi_bsvbs = {{0.f, 0.f}};
                Jpsi_bsvbsx = {{0.f, 0.f}};
                Jpsi_bsvbsy = {{0.f, 0.f}};
                Jpsi_bsvhit = {{0.f, 0.f}};
                Jpsi_bsvms = {{0.f, 0.f}};
                Jpsi_bsvioni = {{0.f, 0.f}};
                Jpsi_bssgnchk = {{0.f, 0.f}};
                Jpsi_bsmeanmass = {{0.f, 0.f, 0.f}};
                Jpsi_bsmeanvtx = {{0.f, 0.f, 0.f}};
                Jpsi_bsmeanbs = {{0.f, 0.f, 0.f, 0.f, 0.f, 0.f}};
                resinfbsv.clear();
                bsvarv.clear();
                cfbsmsv.clear();
                cfbsdelv.clear();
                cfbsiorev.clear();
                cfbsioimv.clear();
                cfbsradrev.clear();
                cfbsradimv.clear();
                cfbshitclsv.clear();
                cfbshitcompv.clear();
                cfbshitvv.clear();
                cfbsgrpv.clear();
                cfbsgrpcompv.clear();
                cfbsgrpmsv.clear();
                cfbsgrpdelv.clear();
                cfbsgrpiorev.clear();
                cfbsgrpioimv.clear();
                cfbsgrpradrev.clear();
                cfbsgrpradimv.clear();
                cfbsgrpvqmsv.clear();
                cfbsgrpvqiov.clear();
                cfbsgrpclosure = {{0.f, 0.f}};

                const unsigned int bsr0 = static_cast<unsigned int>(bsrow);
                const VectorXd rhobs = rfull + Fsparse*dxfree;
                const Matrix<double, 3, 1> rhoB = rhobs.segment<3>(bsr0);
                const Matrix<double, 3, 3> Abs = covstate.block<3, 3>(7, 7);
                const Matrix<double, 3, 3> Mbs = bscovBS - Abs;

                const Matrix<double, 3, 3> covBSinv = bscovBS.inverse();
                Jpsi_bschi2fit = float((rhoB.transpose()*covBSinv*rhoB)(0, 0));
                Jpsi_bsvtx = {{float(bsspot[0] + rhoB[0]),
                               float(bsspot[1] + rhoB[1]),
                               float(bsspot[2] + rhoB[2])}};
                Jpsi_bsspot = {{float(bsspot[0]), float(bsspot[1]), float(bsspot[2])}};
                Jpsi_bsslope = {{float(bsslope[0]), float(bsslope[1])}};
                Jpsi_bswidth = {{float(bswidth[0]), float(bswidth[1]), float(bswidth[2])}};
                Jpsi_bswidtherr = {{float(bswidtherr[0]), float(bswidtherr[1])}};

                // THE CONDITIONING FIX.  `M = covBS - A` is a DIFFERENCE of
                // two nearly equal covariances whenever the two tracks barely
                // constrain the vertex in some direction (A -> covBS), and
                // `A` comes out of a large sparse solve with its own
                // numerical error.  Building the innovation as
                // `covBS M^-1 covBS` then evaluates the WHOLE answer through
                // that amplification and can even return a non-positive
                // matrix -- which is what produced the 5 sigma tail.
                //
                // The "+" form puts the leading order in explicitly and
                // leaves only the CORRECTION to be amplified:
                //     C_{-B} = A + A M^-1 A            (the rows-OFF vertex
                //                                       covariance)
                //     Cov(e) = C_{-B} + covBS          (manifestly >= covBS,
                //                                       so positive definite)
                //     e_B    = rho_B + A M^-1 rho_B    (= covBS M^-1 rho_B)
                // Algebraically identical, numerically not: the answer is
                // never smaller than its own leading term, and `Cov(e)` can
                // no longer come back indefinite.  `Jpsi_bscovlo` exports
                // `C_{-B}`, which the rows-OFF run's `Jpsi_covvtx` must
                // reproduce (gate G7a).
                //
                // The guard below is a DEFINEDNESS test, not a clip: when the
                // vertex is determined by the beam rows alone (`A -> covBS`,
                // `M -> 0`) the leave-one-out residual does not exist, so it
                // is flagged (`Jpsi_bsok = false`) and not written as a
                // number.  It is tested on the smallest eigenvalue of `M`
                // relative to `covBS`'s own scale, which is dimensionally
                // clean because both are covariances of the same three
                // coordinates, and that ratio is EXPORTED (`Jpsi_bsmeig`) so
                // the degenerate corner can be studied rather than guessed.
                SelfAdjointEigenSolver<Matrix<double, 3, 3>> esM(Mbs);
                const double mmin = esM.eigenvalues().minCoeff();
                const double bsscale = bscovBS.diagonal().minCoeff();
                Jpsi_bsmeig = float(bsscale > 0. ? mmin/bsscale : 0.);
                if (mmin > 1.e-9*bsscale && bsscale > 0.) {
                  const LLT<Matrix<double, 3, 3>> lltM(Mbs);
                  const Matrix<double, 3, 3> AMinv =
                      lltM.info() == Eigen::Success
                          ? Matrix<double, 3, 3>(lltM.solve(Abs).transpose())
                          : Matrix<double, 3, 3>(Abs*Mbs.inverse());
                  const Matrix<double, 3, 1> eB = rhoB + AMinv*rhoB;
                  const Matrix<double, 3, 3> Clo = Abs + AMinv*Abs;
                  const Matrix<double, 3, 3> Se = Clo + bscovBS;
                  Jpsi_bscovlo = {{float(Clo(0, 0)), float(Clo(0, 1)), float(Clo(0, 2)),
                                   float(Clo(1, 1)), float(Clo(1, 2)), float(Clo(2, 2))}};

                  Matrix<double, 2, 3> P = Matrix<double, 2, 3>::Zero();
                  P(0, 0) = 1.;
                  P(1, 1) = 1.;
                  P(0, 2) = -bscovBS(0, 2)/bscovBS(2, 2);
                  P(1, 2) = -bscovBS(1, 2)/bscovBS(2, 2);

                  const Matrix<double, 2, 1> rbs = P*eB;
                  const Matrix<double, 2, 2> Cbs = P*Se*P.transpose();

                  Jpsi_bsres = {{float(rbs[0]), float(rbs[1])}};
                  Jpsi_bscov = {{float(Cbs(0, 0)), float(Cbs(0, 1)), float(Cbs(1, 1))}};

                  // LOWER-CHOLESKY whitening in the order (x, y): Cbs = L L^T
                  // with L lower triangular, z = L^-1 r.  z[0] is the x pull
                  // and z[1] the y pull GIVEN x.  The basis is a CHOICE and
                  // this is the one stated in the branch docs; the chi2
                  // below is basis independent.
                  const LLT<Matrix<double, 2, 2>> llt(Cbs);
                  const bool lok = llt.info() == Eigen::Success;
                  Matrix<double, 2, 2> Linv = Matrix<double, 2, 2>::Identity();
                  if (lok) {
                    const Matrix<double, 2, 1> zbs = llt.matrixL().solve(rbs);
                    Jpsi_bsz = {{float(zbs[0]), float(zbs[1])}};
                    Jpsi_bschi2 = float(zbs.squaredNorm());
                    Jpsi_bsok = true;
                    Linv = llt.matrixL().solve(Matrix<double, 2, 2>::Identity());
                    Jpsi_bslinv = {{float(Linv(0, 0)), float(Linv(1, 0)), float(Linv(1, 1))}};
                  }

                  // ---- the two influence rows, and everything at them ----
                  // THE TWO FUNCTIONALS ARE THE WHITENED PAIR, not the global
                  // x and y components.  `r = P e_B = Gbs rho_B` with
                  // `Gbs = P (I + A M^-1)` (the "+" form again), so the
                  // whitened pulls are `z = L^-1 r = (L^-1 Gbs) rho_B` and
                  // the functional's row is that product.  By construction
                  // `Cov(z) = I`: each pull has unit variance -- so the
                  // closure below is `sum_b |a_b|^2 == 1` -- and the two are
                  // UNCORRELATED, which is what stops their product from
                  // over-counting (the global pair's sandwich/quoted was
                  // 1.385 because `Cov_xy` is not zero).  The price, stated:
                  // neither pull is the response to a single beam-spot
                  // parameter any more; `Jpsi_bslinv` is exported so the
                  // global-basis response is one 2x2 multiply away.
                  const Matrix<double, 2, 3> Gbs =
                      Linv*(P*(Matrix<double, 3, 3>::Identity() + AMinv));
                  const unsigned int nb_bs = dVs.size();
                  resinfbsv.assign(2*5*nb_bs, 0.f);
                  bsvarv.assign(2*nb_bs, 0.f);
                  std::array<std::vector<float>, 2> vabsbs;
                  std::array<std::vector<float>, 2> sgnbs;
                  double vchkmax = 0.;
                  for (unsigned int k = 0; k < 2; ++k) {
                    VectorXd ubs = VectorXd::Zero(ncons);
                    ubs.segment<3>(bsr0) = Gbs.row(k).transpose();
                    VectorXd absfree = VectorXd::Zero(nstatefree);
                    for (unsigned int i = 0; i < nstatefree; ++i) {
                      const Eigen::Index is = freestateidxs[i];
                      if (is >= 7 && is <= 9) {
                        absfree(i) = Gbs(k, is - 7);
                      }
                    }
                    const VectorXd wbs = ubs - VinvF*Cinvd.solve(absfree);

                    for (unsigned int j = 0; j < 3; ++j) {
                      Jpsi_bsmeanbs[k*3 + j] = float(-wbs(bsr0 + j));
                    }

                    vabsbs[k].assign(nb_bs, 0.f);
                    sgnbs[k].assign(nb_bs, 1.f);
                    double vsumk = 0., vbsk = 0., vhitk = 0., vmsk = 0., vionik = 0.;
                    double vbsxk = 0., vbsyk = 0.;
                    for (unsigned int ires = 0; ires < nb_bs; ++ires) {
                      const unsigned int r0b = resblockrng[ires][0];
                      const unsigned int nbb = resblockrng[ires][1];
                      const VectorXd wb = wbs.segment(r0b, nbb);
                      const VectorXd ab = ressqrtdV[ires]*wb;
                      for (unsigned int j = 0; j < 5; ++j) {
                        resinfbsv[(k*nb_bs + ires)*5 + j] = j < nbb ? float(ab(j)) : 0.f;
                      }
                      const double vb = ab.squaredNorm();
                      vabsbs[k][ires] = float(vb);
                      const int fam = ires < resfamily_.size() ? resfamily_[ires] : -1;
                      if (fam != 15) {
                        vsumk += vb;
                      }
                      if (fam == 10) {
                        vmsk += vb;
                      } else if (fam == 11) {
                        vionik += vb;
                      } else if (fam == 8 || fam == 9) {
                        vhitk += vb;
                      } else if (fam == kBeamSpotFamily) {
                        vbsk += vb;
                        {
                          const Matrix<double, 3, 1> wbb = wbs.segment<3>(r0b);
                          const double cxx = bscovBS(0, 0), cyy = bscovBS(1, 1), czz = bscovBS(2, 2);
                          const double cxy = bscovBS(0, 1), cxz = bscovBS(0, 2), cyz = bscovBS(1, 2);
                          vbsxk = wbb(0)*wbb(0)*cxx + wbb(0)*wbb(1)*cxy + wbb(0)*wbb(2)*cxz;
                          vbsyk = wbb(1)*wbb(1)*cyy + wbb(0)*wbb(1)*cxy + wbb(1)*wbb(2)*cyz;
                          (void)czz;
                        }
                      }
                      // the per-block ionization sign, the same rule as the
                      // vertex functional's (the two legs carry opposite
                      // charges, so the sign travels with the block)
                      if (fam == 11 && nbb > 0 && resionidir[ires].size() == Eigen::Index(nbb)) {
                        const VectorXd &udir = resionidir[ires];
                        const unsigned int leg = (ires >= resLegStart_[1]) ? 1u : 0u;
                        const double qleg = refftsarr[leg][6] >= 0. ? 1. : -1.;
                        sgnbs[k][ires] = float(qleg*(wb.dot(udir) > 0. ? -1. : 1.));
                      }
                    }
                    // the whitened functional has UNIT variance by
                    // construction; `ckk` stays the normaliser so the closure
                    // reads as `sum_b |a_b|^2 / Cov_kk - 1` exactly as before
                    const double ckk = lok ? 1. : Cbs(k, k);
                    if (ckk > 0.) {
                      vchkmax = std::max(vchkmax, std::abs(vsumk/ckk - 1.));
                      Jpsi_bsvbs[k] = float(vbsk/ckk);
                      Jpsi_bsvbsx[k] = float(vbsxk/ckk);
                      Jpsi_bsvbsy[k] = float(vbsyk/ckk);
                      Jpsi_bsvhit[k] = float(vhitk/ckk);
                      Jpsi_bsvms[k] = float(vmsk/ckk);
                      Jpsi_bsvioni[k] = float(vionik/ckk);
                      for (unsigned int ires = 0; ires < nb_bs; ++ires) {
                        bsvarv[k*nb_bs + ires] = float(vabsbs[k][ires]/ckk);
                      }
                    }

                    // per-hit-class Gaussian shares of Cov_kk
                    if (ckk > 0.) {
                      std::array<double, kNHitResClasses> vcls{};
                      for (std::size_t ires = 0; ires < reshitcls.size() && ires < nb_bs; ++ires) {
                        const int c = reshitcls[ires];
                        if (c < 0 || c >= kNHitResClasses) {
                          continue;
                        }
                        vcls[c] += vabsbs[k][ires];
                      }
                      for (int c = 0; c < kNHitResClasses; ++c) {
                        if (vcls[c] == 0.) {
                          continue;
                        }
                        cfbshitcompv.push_back(static_cast<short>(k));
                        cfbshitclsv.push_back(static_cast<short>(c));
                        cfbshitvv.push_back(float(vcls[c]/ckk));
                      }
                    }

                    // ---- the CF exponents at THIS beam weight ------------
                    // The third and fourth functional of the same fit. They
                    // no longer cost a `trackExponents` call each: the four
                    // weight sets are evaluated together on the concatenated
                    // argument list after this block (see the registration
                    // above), so the pooling and the step-record work they
                    // used to repeat is paid once.
                    if (exportCfExponents_ && ckk > 0.) {
                      cfslotbs[k] = cfRegister(vabsbs[k].data(), int(vabsbs[k].size()), sgnbs[k].data(),
                                               std::sqrt(ckk), 1.);
                    }
                  }
                  Jpsi_bsvchk = float(vchkmax);
                }

                // THE MEAN TERM. A shift `delta` of the centroid moves the
                // beam residual row by `-delta`, so for ANY functional
                // `theta = w^T n`, `d theta / d b0 = -w[bs rows]`. The slope
                // response is that same entry times `(z_v - z0)`, which is
                // `Jpsi_bsvtx[2] - Jpsi_bsspot[2]`. No new global parameter
                // is added to the quadratic term by writing these.
                for (unsigned int j = 0; j < 3; ++j) {
                  Jpsi_bsmeanmass[j] = float(-wmass(bsr0 + j));
                }
                if (exportVtxResidual_ && wvtxinf.size() == Eigen::Index(ncons)) {
                  for (unsigned int j = 0; j < 3; ++j) {
                    Jpsi_bsmeanvtx[j] = float(-wvtxinf(bsr0 + j));
                  }
                }
              }

              // ========= ONE EVALUATOR PASS FOR EVERY FUNCTIONAL ==========
              // Every functional registered above is evaluated here together,
              // on the concatenated argument list { w_{b,k} tau_j }, and its
              // results are written to exactly the branches its own call used
              // to write, in the same order and layout.
              if (!cfins.empty()) {
                cfress.resize(cfins.size());
                cvhcf::trackExponents(cfins.data(), int(cfins.size()), cfress.data());

                auto storecf = [](const std::array<double, cvhcf::kNTau> &a, std::vector<float> &v) {
                  v.resize(cvhcf::kNTau);
                  for (int j = 0; j < cvhcf::kNTau; ++j)
                    v[j] = float(a[j]);
                };

                if (cfslotmass >= 0) {
                  const cvhcf::TrackResult &cfres = cfress[cfslotmass];
                  cfok = cfres.ok;
                  cfnblock = cfres.nblockms + cfres.nblockioni;
                  cfnpooled = cfres.npooled;
                  storecf(cfres.S.ms, cfmsv);
                  storecf(cfres.S.del, cfdelv);
                  storecf(cfres.S.ioRe, cfiorev);
                  storecf(cfres.S.ioIm, cfioimv);
                  storecf(cfres.S.radRe, cfradrev);
                  storecf(cfres.S.radIm, cfradimv);
                  storeCfGroups(cfres);
                }

                if (cfslotvtx >= 0) {
                  const cvhcf::TrackResult &cfresv = cfress[cfslotvtx];
                  Jpsi_vtxok = cfresv.ok;
                  storecf(cfresv.S.ms, cfvtxmsv);
                  storecf(cfresv.S.del, cfvtxdelv);
                  storecf(cfresv.S.ioRe, cfvtxiorev);
                  storecf(cfresv.S.ioIm, cfvtxioimv);
                  storecf(cfresv.S.radRe, cfvtxradrev);
                  storecf(cfresv.S.radIm, cfvtxradimv);
                  if (exportCfGroupExponents_) {
                    storeCfGroupsTo(cfresv,
                                    cfvtxgrpv,
                                    cfvtxgrpmsv,
                                    cfvtxgrpdelv,
                                    cfvtxgrpiorev,
                                    cfvtxgrpioimv,
                                    cfvtxgrpradrev,
                                    cfvtxgrpradimv,
                                    cfvtxgrpclosure,
                                    false,
                                    &cfvtxgrpvqmsv,
                                    &cfvtxgrpvqiov);
                  }
                }

                // the two beam components APPEND, so component 0 must be
                // written before component 1 exactly as the loop did
                for (unsigned int k = 0; k < 2; ++k) {
                  if (cfslotbs[k] < 0) {
                    continue;
                  }
                  const cvhcf::TrackResult &cfresb = cfress[cfslotbs[k]];
                  auto appendcf = [](const std::array<double, cvhcf::kNTau> &a, std::vector<float> &v) {
                    for (int j = 0; j < cvhcf::kNTau; ++j)
                      v.push_back(float(a[j]));
                  };
                  appendcf(cfresb.S.ms, cfbsmsv);
                  appendcf(cfresb.S.del, cfbsdelv);
                  appendcf(cfresb.S.ioRe, cfbsiorev);
                  appendcf(cfresb.S.ioIm, cfbsioimv);
                  appendcf(cfresb.S.radRe, cfbsradrev);
                  appendcf(cfresb.S.radIm, cfbsradimv);
                  if (exportCfGroupExponents_) {
                    std::vector<short> gtmpv;
                    std::vector<float> gms, gdel, giore, gioim, gradre, gradim, gvqms, gvqio;
                    float gclos = 0.f;
                    storeCfGroupsTo(cfresb, gtmpv, gms, gdel, giore, gioim, gradre, gradim,
                                    gclos, false, &gvqms, &gvqio);
                    cfbsgrpclosure[k] = gclos;
                    for (std::size_t ig = 0; ig < gtmpv.size(); ++ig) {
                      cfbsgrpv.push_back(gtmpv[ig]);
                      cfbsgrpcompv.push_back(static_cast<short>(k));
                      cfbsgrpvqmsv.push_back(gvqms[ig]);
                      cfbsgrpvqiov.push_back(gvqio[ig]);
                      for (int j = 0; j < cvhcf::kNTau; ++j) {
                        cfbsgrpmsv.push_back(gms[ig*cvhcf::kNTau + j]);
                        cfbsgrpiorev.push_back(giore[ig*cvhcf::kNTau + j]);
                        cfbsgrpioimv.push_back(gioim[ig*cvhcf::kNTau + j]);
                        cfbsgrpradrev.push_back(gradre[ig*cvhcf::kNTau + j]);
                        cfbsgrpradimv.push_back(gradim[ig*cvhcf::kNTau + j]);
                      }
                    }
                  }
                  Jpsi_bssgnchk[k] = cfresb.ok ? 1.f : 0.f;
                }
              }
            }

          }
  // 
  // (jacarr[idxplus].topLeftCorner(5, nstateparms)*dxdparms.transpose() + jacarr[idxplus].topRightCorner(5, npars)).cast<float>();
  // 
  // Muminus_jacRef.resize(3*npars);
  // Map<Matrix<float, 3, Dynamic, RowMajor>>(Muminus_jacRef.data(), 3, npars) = (jacarr[idxminus].topLeftCorner(5, nstateparms)*dxdparms.transpose() + jacarr[idxminus].topRightCorner(5, npars)).cast<float>();
          
          //TODO fix this
  // Muplus_jacRef.resize(5*npars);
  // Map<Matrix<float, 5, Dynamic, RowMajor>>(Muplus_jacRef.data(), 5, npars) = (jacarr[idxplus].topLeftCorner(5, nstateparms)*dxdparms.transpose() + jacarr[idxplus].topRightCorner(5, npars)).cast<float>();
  // 
  // Muminus_jacRef.resize(5*npars);
  // Map<Matrix<float, 5, Dynamic, RowMajor>>(Muminus_jacRef.data(), 5, npars) = (jacarr[idxminus].topLeftCorner(5, nstateparms)*dxdparms.transpose() + jacarr[idxminus].topRightCorner(5, npars)).cast<float>();
          
          auto const jpsimom = muarr[0] + muarr[1];
          
          if (icons == 0) {
            Jpsi_pt = jpsimom.pt();
            Jpsi_eta = jpsimom.eta();
            Jpsi_phi = jpsimom.phi();
            Jpsi_mass = jpsimom.mass();
          }
          else {
            Jpsicons_pt = jpsimom.pt();
            Jpsicons_eta = jpsimom.eta();
            Jpsicons_phi = jpsimom.phi();
            Jpsicons_mass = jpsimom.mass();
          }
          
          Muplus_nhits = nhitsarr[idxplus];
          Muplus_nvalid = nvalidarr[idxplus];
          Muplus_nvalidpixel = nvalidpixelarr[idxplus];
          Muplus_nvalidFinal = nvalidFinalarr[idxplus];
          Muplus_nvalidpixelFinal = nvalidpixelFinalarr[idxplus];
          Muplus_nmatchedvalid = nmatchedvalidarr[idxplus];
          Muplus_nambiguousmatchedvalid = nambiguousmatchedvalidarr[idxplus];
          
          Muminus_nhits = nhitsarr[idxminus];
          Muminus_nvalid = nvalidarr[idxminus];
          Muminus_nvalidpixel = nvalidpixelarr[idxminus];
          Muminus_nvalidFinal = nvalidFinalarr[idxminus];
          Muminus_nvalidpixelFinal = nvalidpixelFinalarr[idxminus];
          Muminus_nmatchedvalid = nmatchedvalidarr[idxminus];
          Muminus_nambiguousmatchedvalid = nambiguousmatchedvalidarr[idxminus];

          // transmission probe: the UNCONSTRAINED pass is the one whose
          // Mu{plus,minus}_pt the dE/dx scan uses, so pin the accumulator to
          // icons==0 (the mass-constrained pass re-propagates and would
          // otherwise overwrite it).
          if (icons == 0) {
            Muplus_dEref = dErefarr[idxplus];
            Muminus_dEref = dErefarr[idxminus];
            Muplus_maxfracloss = maxFracLossArr[idxplus];
            Muminus_maxfracloss = maxFracLossArr[idxminus];
          }

          Muplus_pixClass.assign(pixclassarr[idxplus].begin(), pixclassarr[idxplus].end());
          Muminus_pixClass.assign(pixclassarr[idxminus].begin(), pixclassarr[idxminus].end());
          Muplus_npixDemoted = npixdemotedarr[idxplus];
          Muminus_npixDemoted = npixdemotedarr[idxminus];
          
          Muplus_highpurity = highpurityarr[idxplus];
          Muminus_highpurity = highpurityarr[idxminus];

          Muplus_charge = muchargearr[idxplus];
          Muminus_charge = muchargearr[idxminus];

          if (readDeDx_) {
            // Find each fit track's index in the ALCARECO source collection by
            // matching the surviving TrackExtraRef.key() (preserved from the
            // original generalTracks-extras through the ALCAReco cloning).
            // NaN if not found.
            auto lookup = [&](const reco::Track& tk, float& outH, float& outP, float& outA) {
              const auto key = tk.extra().key();
              for (size_t k = 0; k < dedxSourceTracksH->size(); ++k) {
                if ((*dedxSourceTracksH)[k].extra().key() == key) {
                  const edm::Ref<reco::TrackCollection> r(dedxSourceTracksH, k);
                  outH = (*dedxHarmonic2H)[r];
                  outP = (*dedxPixelHarmonic2H)[r];
                  outA = (*dedxAllHarmonic2H)[r];
                  return;
                }
              }
              outH = std::numeric_limits<float>::quiet_NaN();
              outP = std::numeric_limits<float>::quiet_NaN();
              outA = std::numeric_limits<float>::quiet_NaN();
            };
            const std::array<const reco::Track*, 2> tkIts = {{ itrack, jtrack }};
            lookup(*tkIts[idxplus],  Muplus_dedxHarmonic2,  Muplus_dedxPixelHarmonic2,  Muplus_dedxAllHarmonic2);
            lookup(*tkIts[idxminus], Muminus_dedxHarmonic2, Muminus_dedxPixelHarmonic2, Muminus_dedxAllHarmonic2);
          }

          const ROOT::Math::PxPyPzMVector mompluskin(outparts[idxplus]->currentState().globalMomentum().x(),
                                                            outparts[idxplus]->currentState().globalMomentum().y(),
                                                            outparts[idxplus]->currentState().globalMomentum().z(),
                                                            trackMass[idxplus]);
          
          const ROOT::Math::PxPyPzMVector momminuskin(outparts[idxminus]->currentState().globalMomentum().x(),
                                                            outparts[idxminus]->currentState().globalMomentum().y(),
                                                            outparts[idxminus]->currentState().globalMomentum().z(),
                                                            trackMass[idxminus]);
          
          auto const jpsimomkin = mompluskin + momminuskin;
          
          if (icons == 0) {
            Jpsikin_pt = jpsimomkin.pt();
            Jpsikin_eta = jpsimomkin.eta();
            Jpsikin_phi = jpsimomkin.phi();
            Jpsikin_mass = jpsimomkin.mass();
          }
          else {
            Jpsikincons_pt = jpsimomkin.pt();
            Jpsikincons_eta = jpsimomkin.eta();
            Jpsikincons_phi = jpsimomkin.phi();
            Jpsikincons_mass = jpsimomkin.mass();            
          }
          
          const reco::Candidate *muplusgen = nullptr;
          const reco::Candidate *muminusgen = nullptr;

          // pre-FSR resonance + the photons the dressed mass needs, gathered
          // in the same pass as the muon match (see the Jpsigenpre_ member
          // docs). The container is an edm::View<reco::Candidate>, which does
          // NOT expose GenStatusFlags, so `isPrompt` is asked for through a
          // dynamic_cast and simply not required when the cast fails.
          const reco::Candidate *resgen = nullptr;
          int resstatus = -1;
          std::vector<const reco::Candidate *> pre746;
          std::vector<const reco::Candidate *> fsrphot;

          Muplusgen_dr = -99.;
          Muminusgen_dr = -99.;
          Muplusgen_pdgId = 0;
          Muminusgen_pdgId = 0;
          Muplusgen_idx = -1;
          Muminusgen_idx = -1;
          Muplusgen_motherPdgId = 0;
          Muminusgen_motherPdgId = 0;
          Muplusgen_motherIdx = -1;
          Muminusgen_motherIdx = -1;
          Muplusgen_isPrompt = false;
          Muminusgen_isPrompt = false;
          Muplusgen_fromHardProcess = false;
          Muminusgen_fromHardProcess = false;
          Jpsigen_sameDecay = false;
          Jpsigenpre_mass = -99.;
          Jpsigenpre_masslep = -99.;
          Jpsigenpre_status = -1;
          Jpsigen_massdressed = -99.;

          if (doGen_) {
            double drminplus = 0.1;
            double drminminus = 0.1;
            // address -> position in the collection, so that the matched
            // particle's mother can be reported as an INDEX (mother() hands
            // back a pointer into this same product) and two candidates that
            // matched the same muon can be recognised downstream.
            std::unordered_map<const reco::Candidate *, int> genidx;
            int igen = -1;

            for (auto const &genpart : *genPartCollection) {
              ++igen;
              genidx[&genpart] = igen;
              const int apid = std::abs(genpart.pdgId());
              const int st = genpart.status();
              // the hard-process resonance: last copy (62) wins over first
              // copy (22); their masses agree to 4e-6 GeV where both exist
              if ((st == 62 || st == 22) &&
                  std::find(genResonancePdgIds_.begin(), genResonancePdgIds_.end(), apid) !=
                      genResonancePdgIds_.end()) {
                if (resgen == nullptr || (st == 62 && resstatus != 62)) {
                  resgen = &genpart;
                  resstatus = st;
                }
              }
              // the pre-Photos lepton copies, present only in radiating events
              if (st == 746 && apid == 13) {
                pre746.push_back(&genpart);
              }
              if (st == 1 && apid == 22) {
                const reco::GenParticle *gp = dynamic_cast<const reco::GenParticle *>(&genpart);
                if (gp == nullptr || gp->statusFlags().isPrompt()) {
                  fsrphot.push_back(&genpart);
                }
              }
              if (st != 1) {
                continue;
              }
              if (apid != 13) {
                continue;
              }

// float dRplus = deltaR(genpart.phi(), muarr[idxplus].phi(), genpart.eta(), muarr[idxplus].eta());
              const double dRplus = deltaR(genpart, muarr[idxplus]);
              if (dRplus < drminplus && genpart.charge() > 0) {
                muplusgen = &genpart;
                Muplusgen_idx = igen;
                drminplus = dRplus;
              }

// float dRminus = deltaR(genpart.phi(), muarr[idxminus].phi(), genpart.eta(), muarr[idxminus].eta());
              const double dRminus = deltaR(genpart, muarr[idxminus]);
              if (dRminus < drminminus && genpart.charge() < 0) {
                muminusgen = &genpart;
                Muminusgen_idx = igen;
                drminminus = dRminus;
              }
            }

            if (muplusgen != nullptr) {
              Muplusgen_dr = drminplus;
            }

            if (muminusgen != nullptr) {
              Muminusgen_dr = drminminus;
            }

            // PROVENANCE. Walk up from the matched particle through its own
            // copies (pdgId == +-13 at any status) to the first ancestor that
            // is something else: the Z for a signal leg, a B/D/K/pi for a
            // decay leg. The guard bounds a pathological self-referential
            // chain; 50 is far above any real decay chain's depth.
            auto firstNonMuonMother = [](const reco::Candidate *p) -> const reco::Candidate * {
              const reco::Candidate *m = p;
              for (int guard = 0; guard < 50 && m != nullptr; ++guard) {
                if (m->numberOfMothers() == 0) {
                  return nullptr;
                }
                m = m->mother(0);
                if (m == nullptr) {
                  return nullptr;
                }
                if (std::abs(m->pdgId()) != 13) {
                  return m;
                }
              }
              return nullptr;
            };
            auto fillProv = [&](const reco::Candidate *mu, int &pdgOut, int &mpdgOut,
                                int &midxOut, bool &promptOut, bool &hardOut) {
              if (mu == nullptr) {
                return;
              }
              pdgOut = mu->pdgId();
              const reco::GenParticle *gp = dynamic_cast<const reco::GenParticle *>(mu);
              if (gp != nullptr) {
                promptOut = gp->statusFlags().isPrompt();
                hardOut = gp->statusFlags().fromHardProcess();
              }
              const reco::Candidate *mother = firstNonMuonMother(mu);
              if (mother != nullptr) {
                mpdgOut = mother->pdgId();
                auto it = genidx.find(mother);
                midxOut = it != genidx.end() ? it->second : -1;
              }
            };
            fillProv(muplusgen, Muplusgen_pdgId, Muplusgen_motherPdgId,
                     Muplusgen_motherIdx, Muplusgen_isPrompt, Muplusgen_fromHardProcess);
            fillProv(muminusgen, Muminusgen_pdgId, Muminusgen_motherPdgId,
                     Muminusgen_motherIdx, Muminusgen_isPrompt, Muminusgen_fromHardProcess);
            // The two legs are one decay only if they are two DIFFERENT gen
            // particles (an index collision is the same muon reconstructed
            // twice) sharing one non-muon ancestor.
            Jpsigen_sameDecay = muplusgen != nullptr && muminusgen != nullptr &&
                                Muplusgen_idx >= 0 && Muminusgen_idx >= 0 &&
                                Muplusgen_idx != Muminusgen_idx &&
                                Muplusgen_motherIdx >= 0 &&
                                Muplusgen_motherIdx == Muminusgen_motherIdx;

            if (resgen != nullptr) {
              Jpsigenpre_mass = resgen->mass();
              Jpsigenpre_status = resstatus;
            }
            if (pre746.size() == 2) {
              Jpsigenpre_masslep =
                  (ROOT::Math::PtEtaPhiMVector(pre746[0]->pt(), pre746[0]->eta(), pre746[0]->phi(), trackMass[0]) +
                   ROOT::Math::PtEtaPhiMVector(pre746[1]->pt(), pre746[1]->eta(), pre746[1]->phi(), trackMass[1]))
                      .mass();
            }

          }
          
          if (muplusgen != nullptr) {
            Muplusgen_pt = muplusgen->pt();
            Muplusgen_eta = muplusgen->eta();
            Muplusgen_phi = muplusgen->phi();
          }
          else {
            Muplusgen_pt = -99.;
            Muplusgen_eta = -99.;
            Muplusgen_phi = -99.;
          }
          
          if (muminusgen != nullptr) {
            Muminusgen_pt = muminusgen->pt();
            Muminusgen_eta = muminusgen->eta();
            Muminusgen_phi = muminusgen->phi();
          }
          else {
            Muminusgen_pt = -99.;
            Muminusgen_eta = -99.;
            Muminusgen_phi = -99.;
          }
          
          if (muplusgen != nullptr && muminusgen != nullptr) {
            auto const jpsigen = ROOT::Math::PtEtaPhiMVector(muplusgen->pt(), muplusgen->eta(), muplusgen->phi(), trackMass[0]) +
                                ROOT::Math::PtEtaPhiMVector(muminusgen->pt(), muminusgen->eta(), muminusgen->phi(), trackMass[1]);
            
            Jpsigen_pt = jpsigen.pt();
            Jpsigen_eta = jpsigen.eta();
            Jpsigen_phi = jpsigen.phi();
            Jpsigen_mass = jpsigen.mass();

            // DRESSED: the bare pair plus every prompt final-state photon
            // within dR < 0.1 of EITHER muon -- the same cone and the same
            // photon collection `zchannel/dump_gen_fsr.py` uses, so the two
            // definitions cannot drift.
            {
              ROOT::Math::PtEtaPhiMVector dressed = jpsigen;
              for (auto const *ph : fsrphot) {
                if (deltaR(*ph, *muplusgen) < 0.1 || deltaR(*ph, *muminusgen) < 0.1) {
                  dressed += ROOT::Math::PtEtaPhiMVector(ph->pt(), ph->eta(), ph->phi(), 0.);
                }
              }
              Jpsigen_massdressed = dressed.mass();
            }

            Jpsigen_x = muplusgen->vx();
            Jpsigen_y = muplusgen->vy();
            Jpsigen_z = muplusgen->vz();

            // genParticles:xyz0 (gen PV) is not kept in every ALCARECO
            // (the B->J/psi+X MC keeps only the recoGenParticles branch);
            // fall back to the -99 sentinel rather than throwing.
            genl3d = genXyz0.isValid()
                ? std::sqrt((muplusgen->vertex() - *genXyz0).mag2()) : -99.;
          }
          else {
            Jpsigen_pt = -99.;
            Jpsigen_eta = -99.;
            Jpsigen_phi = -99.;
            Jpsigen_mass = -99.;
            
            Jpsigen_x = -99.;
            Jpsigen_y = -99.;
            Jpsigen_z = -99.;

            genl3d = -99.;
          }


          Muplus_isMuon = false;
          Muplus_muonLoose = false;
          Muplus_muonMedium = false;
          Muplus_muonTight = false;
          Muplus_muonIsPF = false;
          Muplus_muonIsTracker = false;
          Muplus_muonIsGlobal = false;
          Muplus_muonIsStandalone = false;
          Muplus_muonInnerTrackBest = false;

          Muminus_isMuon = false;
          Muminus_muonLoose = false;
          Muminus_muonMedium = false;
          Muminus_muonTight = false;
          Muminus_muonIsPF = false;
          Muminus_muonIsTracker = false;
          Muminus_muonIsGlobal = false;
          Muminus_muonIsStandalone = false;
          Muminus_muonInnerTrackBest = false;

          if (doMuons_) {
            const reco::Muon *matchedmuonplus = idxplus == 0 ? matchedmuon0 : matchedmuon1;
            const reco::Muon *matchedmuonminus = idxminus == 0 ? matchedmuon0 : matchedmuon1;

            if (matchedmuonplus != nullptr) {
              Muplus_isMuon = true;
              Muplus_muonLoose = matchedmuonplus->passed(reco::Muon::CutBasedIdLoose);
              Muplus_muonMedium = matchedmuonplus->passed(reco::Muon::CutBasedIdMedium);
              Muplus_muonTight = matchedmuonplus->passed(reco::Muon::CutBasedIdTight);
              Muplus_muonIsPF = matchedmuonplus->isPFMuon();
              Muplus_muonIsTracker = matchedmuonplus->isTrackerMuon();
              Muplus_muonIsGlobal = matchedmuonplus->isGlobalMuon();
              Muplus_muonIsStandalone = matchedmuonplus->isStandAloneMuon();
              Muplus_muonInnerTrackBest = matchedmuonplus->muonBestTrackType() == reco::Muon::InnerTrack;
            }

            if (matchedmuonminus != nullptr) {
              Muminus_isMuon = true;
              Muminus_muonLoose = matchedmuonminus->passed(reco::Muon::CutBasedIdLoose);
              Muminus_muonMedium = matchedmuonminus->passed(reco::Muon::CutBasedIdMedium);
              Muminus_muonTight = matchedmuonminus->passed(reco::Muon::CutBasedIdTight);
              Muminus_muonIsPF = matchedmuonminus->isPFMuon();
              Muminus_muonIsTracker = matchedmuonminus->isTrackerMuon();
              Muminus_muonIsGlobal = matchedmuonminus->isGlobalMuon();
              Muminus_muonIsStandalone = matchedmuonminus->isStandAloneMuon();
              Muminus_muonInnerTrackBest = matchedmuonminus->muonBestTrackType() == reco::Muon::InnerTrack;
            }

          }
          
          
      // const Vector5d dxRef = dxfull.head<5>();
      // const Matrix5d Cinner = Cinvd.solve(MatrixXd::Identity(nstateparms,nstateparms)).topLeftCorner<5,5>();


          niter = iiter + 1;
          edmval = -deltachisq;

          // The reference-block EDM, and with it the convergence break, is
          // taken on the reference indices ACTUALLY SOLVED FOR. `covstate` is
          // the free-subspace covariance scattered into the full index space,
          // so a frozen reference index leaves a zero row and column: under
          // doVtxConstraint (index 6) the 10x10 block is singular, its inverse
          // is NaN, `edmvalref < edmConvergence_` is never true and every
          // candidate runs to `nIters`. Under fitFromGenParms the whole block
          // is frozen and the same happens.
          const Matrix<double, 10, 1> dxRef = dxfull.head<10>();
          std::vector<Eigen::Index> reffreeidxs;
          reffreeidxs.reserve(10);
          for (auto const idx : freestateidxs) {
            if (idx < 10) {
              reffreeidxs.push_back(idx);
            }
          }
          double deltachisqref = 0.;
          if (!reffreeidxs.empty()) {
            const MatrixXd covref = covstate(reffreeidxs, reffreeidxs);
            const MatrixXd hessref = covref.inverse();
            const VectorXd dxReffree = dxRef(reffreeidxs);
            deltachisqref = -0.5*dxReffree.transpose()*hessref*dxReffree;
          }

          edmvalref = -deltachisqref;

          if (std::isnan(edmval) || std::isinf(edmval)) {
            std::cout << "WARNING: invalid parameter update!!!" << " edmval = " << edmval << " deltachisqval = " << deltachisqval << std::endl;
            ++fitFailNaN_;
            valid = false;
            break;
          }
          
          if (icons == 0) {
            edmval_cons0 = edmval;
            niter_cons0 = niter;
          }

          // openspec/improve-cvh-refit-convergence §1.5: per-iter debug
          // dump. Records the full trace (chisqval, edmval, deltachisqval,
          // per-muon q/p, dimuon mass) across both icons phases for this
          // candidate. Pushed before the convergence-break check so the
          // break-triggering iteration is included.
          if (debugPerIterDump_) {
            // stdout mirror of the per-iter trace: unlike the tree branches
            // this also survives for candidates whose fit later aborts.
            std::cout << "dbgIter: icons=" << icons << " iiter=" << iiter
                      << " chisq=" << chisqval
                      << " deltachisq=" << deltachisqval
                      << " edmval=" << edmval
                      << " dxvtx10max=" << dxfull.head<10>().cwiseAbs().maxCoeff()
                      << " ref0(r,z,p)=(" << std::hypot(refftsarr[0][0], refftsarr[0][1])
                      << "," << refftsarr[0][2] << "," << refftsarr[0].segment<3>(3).norm() << ")"
                      << " ref1(r,z,p)=(" << std::hypot(refftsarr[1][0], refftsarr[1][1])
                      << "," << refftsarr[1][2] << "," << refftsarr[1].segment<3>(3).norm() << ")"
                      << std::endl;
            chisqval_iter.push_back(static_cast<double>(chisqval));
            edmval_iter.push_back(static_cast<double>(edmval));
            edmvalref_iter.push_back(static_cast<double>(edmvalref));
            deltachisqval_iter.push_back(static_cast<double>(deltachisqval));
            for (unsigned int id = 0; id < 2; ++id) {
              const double px = refftsarr[id][3];
              const double py = refftsarr[id][4];
              const double pz = refftsarr[id][5];
              const double pmag = std::sqrt(px*px + py*py + pz*pz);
              const double qop = (pmag > 0.0)
                                     ? refftsarr[id][6] / pmag
                                     : 0.0;
              mu_qoverp_iter.push_back(qop);
            }
            const double m0 = daughterMass1_;
            const double m1 = daughterMass2_;
            const double E0 = std::sqrt(refftsarr[0][3] * refftsarr[0][3] +
                                        refftsarr[0][4] * refftsarr[0][4] +
                                        refftsarr[0][5] * refftsarr[0][5] +
                                        m0 * m0);
            const double E1 = std::sqrt(refftsarr[1][3] * refftsarr[1][3] +
                                        refftsarr[1][4] * refftsarr[1][4] +
                                        refftsarr[1][5] * refftsarr[1][5] +
                                        m1 * m1);
            const double Px = refftsarr[0][3] + refftsarr[1][3];
            const double Py = refftsarr[0][4] + refftsarr[1][4];
            const double Pz = refftsarr[0][5] + refftsarr[1][5];
            const double m2 = (E0 + E1) * (E0 + E1) - Px * Px - Py * Py - Pz * Pz;
            Jpsi_mass_iter.push_back(m2 > 0.0 ? std::sqrt(m2) : 0.0);
          }

// std::cout << "icons = " << icons << " iiter = " << iiter << " edmval = " << edmval << " deltachisqval = " << deltachisqval << " chisqval = " << chisqval << std::endl;
// std::cout << "dxvtx" << std::endl;
          
// std::cout << "pt0 = " << refftsarr[0].momentum().perp() << " eta0 = " << refftsarr[0].momentum().eta() << " charge0 = " << refftsarr[0].charge() << " pt1 = " << refftsarr[1].momentum().perp() << " eta1 = " << refftsarr[1].momentum().eta() << " charge1 = " << refftsarr[1].charge() << std::endl;
// std::cout << "dxref" << std::endl;
          
// std::cout << "dxvtx:" << std::endl;
// std::cout << dxfull.head<3>() << std::endl;
// std::cout << "dxmom0" << std::endl;
// std::cout << dxfull.segment<3>(trackstateidxarr[0]) << std::endl;
// std::cout << "dxmom1" << std::endl;
// std::cout << dxfull.segment<3>(trackstateidxarr[1]) << std::endl;
// std::cout << "qop0 = " << refftsarr[0].signedInverseMomentum() + dxfull[trackstateidxarr[0]] << std::endl;
// std::cout << "qop1 = " << refftsarr[1].signedInverseMomentum() + dxfull[trackstateidxarr[1]] << std::endl;
          
          
// std::cout << "dx0" << std::endl;
// std::cout << dxfull.segment(trackstateidxarr[0], trackstateidxarr[1]-trackstateidxarr[0]) << std::endl;
// std::cout << "dx1" << std::endl;
// std::cout << dxfull.segment(trackstateidxarr[1], nstateparms - trackstateidxarr[1]) << std::endl;
// std::cout << dxfull.segment<3>(trackstateidxarr[0]) << std::endl;
// std::cout << dxfull.segment<3>(trackstateidxarr[1]) << std::endl;
         
// if (std::abs(deltachisqval)<1e-2) {
// break;
// }
          
// if (iiter > 0 && edmval < 1e-5) {
// break;
// }
          
          // openspec/improve-cvh-refit-convergence §1.3: threshold configurable
          // via edmConvergence_ (cfi default 1e-5 reproduces the baseline).
          if (iiter > 0 && dolocalupdate && edmval < edmConvergence_) {
            break;
          }
          else if (iiter > 0 && !dolocalupdate && edmvalref < edmConvergence_) {
            break;
          }
          
// if (iiter > 1 && std::abs(deltachisq[0])<1e-3) {
// break;
// }
      
        }
      
        if (!valid) {
          break;
        }

        std::unordered_map<unsigned int, unsigned int> idxmap;

        globalidxvfinal.clear();
        globalidxvfinal.reserve(npars);
        idxmap.reserve(npars);

        for (unsigned int idx : globalidxv) {
          if (!idxmap.count(idx)) {
            idxmap[idx] = globalidxvfinal.size();
            globalidxvfinal.push_back(idx);
          }
        }

        // ---- VARIANCE-ONLY COLUMNS (exportVarianceGrads) ------------------
        //
        // A parmtype-8/9/10/11 parameter enters this maker's -2lnL ONLY
        // through the covariance: unlike the single-track maker, the
        // two-track parameter vector is alignment + field + eloss/material
        // and carries no resolution slots (`npars = nparsAlignment +
        // nparsBfield + nparsEloss`).  So those families need columns
        // APPENDED here; their `Jfinal` columns stay exactly zero and only
        // the log-det block below writes to them.
        //
        // Parmtype 15 needs NOTHING appended -- `matGroupGlobalIdx_[g]` is
        // already a column of `globalidxv`, one slot per group per hit --
        // which is why `varianceGradFamilies = {15}` preserves the layout
        // and can be pooled with a production that ran without the switch.
        const std::size_t nparsmean = globalidxvfinal.size();
        if (dores && exportVarianceGrads_) {
          for (unsigned int ires = 0; ires < resglobidx.size(); ++ires) {
            const int fam = ires < resfamily_.size() ? resfamily_[ires] : -1;
            if (!varianceFamilyWanted(fam)) {
              continue;
            }
            const unsigned int idx = resglobidx[ires];
            if (!idxmap.count(idx)) {
              idxmap[idx] = globalidxvfinal.size();
              globalidxvfinal.push_back(idx);
            }
          }
        }
        (void)nparsmean;

        const unsigned int nparsfinal = globalidxvfinal.size();

        // Sparse GBL reduction (mirrors single-track maker 2536-2651).
        // Column-gather Jfull (ncons x npars, per-hit-replicated field
        // modes) down to Jfinal (ncons x nparsfinal, deduped global
        // params) via the same idxmap the dense path used. Replaces the
        // O(npars^2) dense d2chisqdparms2 gather (was ~164M iters at 360
        // modes) with an O(ncons*npars) sparse column-add.
        MatrixXd Jfinal = MatrixXd::Zero(ncons, nparsfinal);
        for (unsigned int i = 0; i < npars; ++i) {
          Jfinal.col(idxmap.at(globalidxv[i])) += Jfull.col(i);
        }
        const SparseMatrix<double> Jsparse = Jfinal.sparseView();

        // Residual projector R = Vinv - VinvF Cinvd^-1 (VinvF)^T and the
        // cheaper Rr = Vinv (r + Fs dxfree) == R r. grad/hess are the
        // GBL reduced gradient/Hessian wrt the global params (stored to
        // gradv/hesspackedv when fillGrads_). Algebraically identical to
        // the dense  grad = dchisqdparmsfinal + d2chisqdxdparmsfinal^T dxfull,
        // hess = d2chisqdparms2final + dxdparms d2chisqdxdparmsfinal
        // (the factor-of-2 / Schur-complement bookkeeping cancels).
        const SparseMatrix<double> FtVinv = VinvF.transpose();
        const SparseMatrix<double> Rsparse = Vinvsparse - VinvF*Cinvd.solve(FtVinv);
        const MatrixXd R = Rsparse;
        const VectorXd Rr = Vinvsparse*(rfull + Fsparse*dxfree);

        grad = 2.*Jsparse.transpose()*Rr;
        hess = 2.*Jsparse.transpose()*R*Jsparse;

        // dxdparms (nparsfinal x nstateparms): sensitivity of the fitted
        // state to the global params. Free-state columns filled from the
        // sparse solve; frozen-state columns stay zero. The Muplus/Muminus
        // _jacRef and Jpsi_jacMass emitters below consume this unchanged.
        dxdparms = MatrixXd::Zero(nparsfinal, nstateparms);
        dxdparms(Eigen::placeholders::all, freestateidxs) =
            -Cinvd.solve(VinvF.transpose()*Jsparse).transpose();

        // ================= THE VARIANCE (log-det) TERM =====================
        //
        // Up to here the exported gradient/Hessian are those of the QUADRATIC
        // form alone, `chi2(theta) = r^T R r` with `r -> r + J theta`, i.e.
        // every parameter enters only through the MEAN of the residuals.  For
        // a parameter that also moves the COVARIANCE -- parmtype 15 (a
        // material group scales its steps' MS covariance and ionization
        // variance by the same `exp(k_g)` that scales their mean loss), and
        // parmtypes 8/9/10/11 (which move nothing else) -- that is not the
        // derivative of the likelihood.
        //
        // The objective is the MARGINAL (REML) one, exactly as in the
        // single-track maker (`ResidualGlobalCorrectionMakerG4e.cc`, the
        // `gradll` loop):
        //
        //     -2 lnL(theta) = r^T R r + ln|V| + ln|C|,
        //         R = V^-1 - V^-1 F C^-1 F^T V^-1,   C = F^T V^-1 F,
        //
        // whose derivatives use  dR/dtheta = -R (dV/dtheta) R  and
        // d(ln|V| + ln|C|)/dtheta = tr(R dV/dtheta):
        //
        //     dG_i   = -r^T R dV_i R r  +  tr(dV_i R)
        //     dH_ij  =  tr(dV_i R dV_j R)                       (EXPECTED)
        //
        // Three properties of this, all deliberate and all shared with the
        // single-track code:
        //
        //  * THE LOCAL TRACK PARAMETERS.  `R` already carries the projection
        //    -V^-1 F C^-1 F^T V^-1, so the fitted state is profiled out and
        //    its implicit derivative vanishes by the envelope theorem; the
        //    `ln|C|` piece -- the difference between profiling and
        //    marginalizing -- is supplied automatically by using `R` rather
        //    than `V^-1` inside the trace.  The vertex / beamspot / pointing
        //    / mass constraint rows have theta-independent weights, so they
        //    register no `dV` and enter only through `R` and `C`.
        //
        //  * THE HESSIAN IS THE EXPECTED (FISHER) ONE.  The observed pieces
        //    `2 r^T R dV_i R dV_j R r` and the mean-variance cross term
        //    `-2 J^T R dV_i R r` are DROPPED -- the single-track maker has
        //    them behind `if (false)` for the same reason.  For a Gaussian
        //    the mean and variance blocks of the Fisher matrix are exactly
        //    orthogonal, so the cross term is zero IN EXPECTATION and keeping
        //    the observed one only adds noise; and the expected form is a
        //    Gram matrix, `tr(dV_i R dV_j R) = <R^1/2 dV_i R^1/2,
        //    R^1/2 dV_j R^1/2>_F`, hence positive semi-definite by
        //    construction, which `2 J^T R J` also is.  The exported `hess` is
        //    therefore PSD whatever the candidate does.
        //
        //  * WHAT dV IS EVALUATED AT.  `dV_i` is the derivative at the point
        //    the FIT USED.  For parmtypes 10/11/15 that is the covariance the
        //    propagator actually built (`k_g` included).  For parmtypes 8/9
        //    this maker does not apply `exp(corparms)` to the hit covariance
        //    at all (the single-track maker does), so those derivatives are
        //    at theta = 0 regardless of any `corFiles` -- see the doc.
        //
        // The blocks are evaluated in their own row range rather than as full
        // ncons x ncons sparse products: every `dV_i` is supported on
        // `resblockrng[i]`, so `tr(dV_i R) = tr(D_i R_ii)` and
        // `tr(dV_i R dV_j R) = tr(D_i R_ij D_j R_ji)` with `D` the small
        // dense block.  Identical algebra, ~5 x 5 matrices instead of
        // 400 x 400 sparse ones.
        nvarcols = 0;
        if (dores && exportVarianceGrads_) {
          // per candidate, so a candidate with no registered block writes an
          // EMPTY gradchisqv/gradllv rather than the previous candidate's
          gradchisqv.clear();
          gradllv.clear();
          hessvaridxv.clear();
          hessvarpackedv.clear();
          nHessVar = 0;
        }
        if (dores && exportVarianceGrads_ && !dVs.empty()) {
          std::vector<unsigned int> vcol, vr0, vnb;
          std::vector<MatrixXd> vD;
          vcol.reserve(dVs.size());
          vr0.reserve(dVs.size());
          vnb.reserve(dVs.size());
          vD.reserve(dVs.size());
          for (unsigned int ires = 0; ires < dVs.size(); ++ires) {
            const int fam = ires < resfamily_.size() ? resfamily_[ires] : -1;
            if (!varianceFamilyWanted(fam)) {
              continue;
            }
            const unsigned int r0 = resblockrng[ires][0];
            const unsigned int nb = resblockrng[ires][1];
            MatrixXd Db(nb, nb);
            for (unsigned int a = 0; a < nb; ++a) {
              for (unsigned int b = 0; b < nb; ++b) {
                Db(a, b) = dVs[ires].coeff(r0 + a, r0 + b);
              }
            }
            if (!(Db.cwiseAbs().maxCoeff() > 0.)) {
              continue;
            }
            vcol.push_back(idxmap.at(resglobidx[ires]));
            vr0.push_back(r0);
            vnb.push_back(nb);
            vD.emplace_back(std::move(Db));
          }

          const unsigned int nv = vcol.size();
          if (nv > 0) {
            {
              std::vector<unsigned int> uniq(vcol);
              std::sort(uniq.begin(), uniq.end());
              uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
              nvarcols = uniq.size();
            }

            VectorXd gradll = VectorXd::Zero(nparsfinal);
            // dVRr(:, p) = sum over the blocks of parameter p of dV_i R r.
            // Dense ncons x nparsfinal, but only the blocks' own rows are
            // ever written.
            MatrixXd dVRr = MatrixXd::Zero(ncons, nparsfinal);

            std::vector<VectorXd> vDRr(nv);
            for (unsigned int i = 0; i < nv; ++i) {
              const auto Rri = Rr.segment(vr0[i], vnb[i]);
              vDRr[i] = vD[i] * Rri;
              dVRr.block(vr0[i], vcol[i], vnb[i], 1) += vDRr[i];
              // log-det gradient: tr(dV_i R)
              gradll(vcol[i]) +=
                  (vD[i] * R.block(vr0[i], vr0[i], vnb[i], vnb[i])).trace();
            }

            // expected Hessian: tr(dV_i R dV_j R) = tr(D_i R_ij D_j R_ji).
            // Accumulated BOTH into `hess` (so `hesspackedv` stays complete)
            // and into its own small dense block over the variance columns,
            // which is what is shipped alongside `hessfactorv` -- see
            // `hessvaridxv` in the base class for why.
            std::vector<unsigned int> varcols(vcol);
            std::sort(varcols.begin(), varcols.end());
            varcols.erase(std::unique(varcols.begin(), varcols.end()),
                          varcols.end());
            std::unordered_map<unsigned int, unsigned int> varpos;
            varpos.reserve(varcols.size());
            for (unsigned int i = 0; i < varcols.size(); ++i) {
              varpos[varcols[i]] = i;
            }
            MatrixXd hessvar = MatrixXd::Zero(varcols.size(), varcols.size());
            for (unsigned int i = 0; i < nv; ++i) {
              const MatrixXd DiRij_base = vD[i];
              const unsigned int pi = varpos.at(vcol[i]);
              for (unsigned int j = 0; j <= i; ++j) {
                const MatrixXd T =
                    DiRij_base * R.block(vr0[i], vr0[j], vnb[i], vnb[j]);   // nb_i x nb_j
                const MatrixXd U =
                    vD[j] * R.block(vr0[j], vr0[i], vnb[j], vnb[i]);        // nb_j x nb_i
                const double hessres = (T.array() * U.transpose().array()).sum();
                const unsigned int pj = varpos.at(vcol[j]);
                hess(vcol[i], vcol[j]) += hessres;
                hessvar(pi, pj) += hessres;
                if (i != j) {
                  hess(vcol[j], vcol[i]) += hessres;
                  hessvar(pj, pi) += hessres;
                }
              }
            }
            // pack the block, upper triangle row-major, columns ascending
            nHessVar = varcols.size();
            hessvaridxv.assign(varcols.begin(), varcols.end());
            hessvarpackedv.clear();
            hessvarpackedv.reserve(nHessVar * (nHessVar + 1) / 2);
            for (unsigned int i = 0; i < nHessVar; ++i) {
              for (unsigned int j = i; j < nHessVar; ++j) {
                hessvarpackedv.push_back(float(hessvar(i, j)));
              }
            }

            // chi2 part of the variance gradient: -r^T R dV_i R r
            const SparseMatrix<double> dVRrsparse = dVRr.sparseView();
            grad += -dVRrsparse.transpose()*Rr;

            gradchisqv.clear();
            gradchisqv.resize(nparsfinal, 0.);
            Map<VectorXf>(gradchisqv.data(), nparsfinal) = grad.cast<float>();

            gradllv.clear();
            gradllv.resize(nparsfinal, 0.);
            Map<VectorXf>(gradllv.data(), nparsfinal) = gradll.cast<float>();

            grad += gradll;

            // The fitted state moves with a variance parameter too:
            //   dxhat/dtheta_i = C^-1 F^T V^-1 dV_i R r
            // (same expression the single-track maker adds).  This is what
            // puts the variance families into `Jpsi_jacMass` and the two
            // `_jacRef`, i.e. what lets the MASS term see them.
            dxdparms(Eigen::placeholders::all, freestateidxs) +=
                Cinvd.solve(FtVinv*dVRrsparse).transpose();
          }
        }

        // ---- IN-MAKER FINITE DIFFERENCE, at FIXED linearization -----------
        //
        // The propagator-level FD (perturb `k_g`, re-fit, difference `objval`)
        // cannot be sharp: the reference trajectory moves with the parameter,
        // so `V` moves with it through paths the analytic `dV` deliberately
        // omits, and the Gauss-Newton stopping tolerance leaves a
        // delta-INDEPENDENT residual that grows as 1/delta.  This one holds
        // `r`, `F` and `J` fixed and perturbs only the covariance,
        // `V -> V + s dV_i`, then re-does the profile from scratch:
        //
        //   obj(s) = r^T R(s) r + ln|V(s)| + ln|C(s)|,  C(s) = F^T V(s)^-1 F
        //
        // so (obj(+s) - obj(-s))/2s must reproduce the exported column to
        // O(s^2).  It tests the traces, the projector, the sign and the ln|C|
        // term -- everything except whether `dV` really is dV/dtheta.
        //
        // Only families 10/11/15 are done: their dV blocks span the 5x5
        // process-noise rows, whose V block is exactly the inverse of the
        // Vinv block there.  A hit block's rows are 2 wide with a
        // (near-)singular Vinv on the unmeasured strip coordinate, so the
        // same trick does not apply and the parmtype-8/9 columns rely on the
        // structural argument instead.
        if (exportVarianceGrads_ && varianceFDGlobalIdx_ != -1 && !dVs.empty()
            && nvarcols > 0) {
          // global index -> its blocks, keyed by row offset (one block per
          // propagation for a given global index, so the ranges are disjoint)
          std::map<unsigned int, std::map<unsigned int, MatrixXd>> perGlob;
          for (unsigned int ires = 0; ires < dVs.size(); ++ires) {
            const int fam = ires < resfamily_.size() ? resfamily_[ires] : -1;
            if (!(fam == 10 || fam == 11 || fam == 15)) {
              continue;
            }
            if (!varianceFamilyWanted(fam)) {
              continue;
            }
            const unsigned int gidx = resglobidx[ires];
            if (varianceFDGlobalIdx_ >= 0 && gidx != unsigned(varianceFDGlobalIdx_)) {
              continue;
            }
            const unsigned int r0 = resblockrng[ires][0];
            const unsigned int nb = resblockrng[ires][1];
            MatrixXd Db = MatrixXd::Zero(nb, nb);
            for (unsigned int a = 0; a < nb; ++a) {
              for (unsigned int b = 0; b < nb; ++b) {
                Db(a, b) = dVs[ires].coeff(r0 + a, r0 + b);
              }
            }
            auto &m = perGlob[gidx];
            auto it = m.find(r0);
            if (it == m.end()) {
              m.emplace(r0, std::move(Db));
            } else {
              it->second += Db;
            }
          }
          const double eps = varianceFDEps_;
          for (auto const &g : perGlob) {
            const unsigned int gidx = g.first;
            auto itcol = idxmap.find(gidx);
            if (itcol == idxmap.end()) {
              continue;
            }
            const unsigned int col = itcol->second;
            std::array<double, 2> objs{{0., 0.}};
            std::array<double, 2> chis{{0., 0.}};
            std::array<double, 2> lls{{0., 0.}};
            bool ok = true;
            for (int is = 0; is < 2 && ok; ++is) {
              const double sgn = is == 0 ? eps : -eps;
              MatrixXd Vinvs = Vinvfull;
              double dlogdetV = 0.;
              for (auto const &blk : g.second) {
                const unsigned int r0 = blk.first;
                const unsigned int nb = blk.second.rows();
                const MatrixXd Vinvb = Vinvfull.block(r0, r0, nb, nb);
                const double dib = Vinvb.determinant();
                if (!(std::abs(dib) > 0.)) {
                  ok = false;
                  break;
                }
                const MatrixXd Vb = Vinvb.inverse();
                const MatrixXd Vbs = Vb + sgn * blk.second;
                const double d0 = Vb.determinant();
                const double d1 = Vbs.determinant();
                if (!(d0 > 0. && d1 > 0.)) {
                  ok = false;
                  break;
                }
                dlogdetV += std::log(d1) - std::log(d0);
                Vinvs.block(r0, r0, nb, nb) = Vbs.inverse();
              }
              if (!ok) {
                break;
              }
              const SparseMatrix<double> Vinvss = Vinvs.sparseView();
              const SparseMatrix<double> VinvFs = Vinvss * Fsparse;
              SimplicialLDLT<SparseMatrix<double>> Cs(SparseMatrix<double>(
                  Fsparse.transpose() * VinvFs));
              if (Cs.info() != Eigen::Success) {
                ok = false;
                break;
              }
              const VectorXd dxs = -Cs.solve(VinvFs.transpose() * rfull);
              chis[is] = rfull.dot(VectorXd(Vinvss * (rfull + Fsparse * dxs)));
              // ln|V_0| is common to the two arms and cancels in the
              // difference, so it is not needed here (and `objlogdetv` is
              // filled below, not above).
              lls[is] = dlogdetV + Cs.vectorD().array().abs().log().sum();
              objs[is] = chis[is] + lls[is];
            }
            if (!ok) {
              continue;
            }
            std::cout << "VARFD glob=" << gidx
                      << " col=" << col
                      << " eps=" << eps
                      << " an=" << std::setprecision(12) << grad(col)
                      << " anchi=" << (col < gradchisqv.size() ? double(gradchisqv[col]) : 0.)
                      << " anll=" << (col < gradllv.size() ? double(gradllv[col]) : 0.)
                      << " fd=" << (objs[0] - objs[1]) / (2. * eps)
                      << " fdchi=" << (chis[0] - chis[1]) / (2. * eps)
                      << " fdll=" << (lls[0] - lls[1]) / (2. * eps)
                      << " ev=" << run << ":" << lumi << ":" << event
                      << std::setprecision(6) << std::endl;
          }
        }

        // The marginal objective itself, for finite-differencing the above.
        if (exportObjective_) {
          objchisq = rfull.dot(Rr);
          // ln|V| = -ln|Vinv|.  Vinv is block diagonal and positive definite
          // (the deweighted strip coordinate carries a tiny but non-zero
          // weight), so an LDLT is enough.
          // ln|V| = -ln|Vinv|, but Vinv is RANK DEFICIENT by construction:
          // the deweighted second coordinate of every 1-D strip hit carries
          // exactly zero weight, so a plain determinant is 0 and its log is
          // -inf.  What is wanted -- and what cancels in a finite difference,
          // because the null space is structural and delta-independent -- is
          // the PSEUDO-determinant over the non-null modes.  `objnullv`
          // records how many were dropped so a FD can assert that the two
          // arms dropped the same number.
          const SelfAdjointEigenSolver<MatrixXd> esV(Vinvfull, EigenvaluesOnly);
          const double lmaxV = esV.eigenvalues().maxCoeff();
          const double cutV = 1e-12 * std::max(lmaxV, 0.);
          double lsum = 0.;
          objnullv = 0;
          for (int k = 0; k < esV.eigenvalues().size(); ++k) {
            const double ev = esV.eigenvalues()(k);
            if (ev > cutV) {
              lsum += std::log(ev);
            } else {
              ++objnullv;
            }
          }
          objlogdetv = -lsum;
          objlogdetc = Cinvd.vectorD().array().abs().log().sum();
          objval = objchisq + objlogdetv + objlogdetc;
        }
        // ===================================================================

        if (icons == 0) {
          const unsigned int idxplus = muchargearr[0] > 0 ? 0 : 1;
          const unsigned int idxminus = muchargearr[0] > 0 ? 1 : 0;

          const unsigned int trackstateidxplus = 3*idxplus;
          const unsigned int trackstateidxminus = 3*idxminus;

          Muplus_jacRef.resize(3*nparsfinal);
          Map<Matrix<float, 3, Dynamic, RowMajor>>(Muplus_jacRef.data(), 3, nparsfinal) = dxdparms.block(0, trackstateidxplus, nparsfinal, 3).transpose().cast<float>();

          Muminus_jacRef.resize(3*nparsfinal);
          Map<Matrix<float, 3, Dynamic, RowMajor>>(Muminus_jacRef.data(), 3, nparsfinal) = dxdparms.block(0, trackstateidxminus, nparsfinal, 3).transpose().cast<float>();

          const Matrix<double, 1, 6> mjacalt =
              massJacobianAltD(refftsarr[0], refftsarr[1], massForConstraintHelpers);

          Jpsi_jacMass.resize(nparsfinal);
          Map<Matrix<float, 1, Dynamic, RowMajor>>(Jpsi_jacMass.data(), 1, nparsfinal) = (mjacalt*dxdparms.leftCols<6>().transpose()).cast<float>();

          // `Jpsi_jacVtx` = d theta_6^unconstrained / d(global params),
          // aligned with `globalidxv` exactly as `Jpsi_jacMass` is. It is the
          // D row a DATA fit needs for the MEAN (alignment / field) part of
          // the vertex term -- cheap, one row.
          //   index 6 free  : it is simply `dxdparms.col(6)`.
          //   index 6 frozen: that column is zero by construction, and the
          //     unconstrained response is `-w_v^T J`, which reduces to
          //     `dxdparms.col(6)` in the free case (same algebra), so the two
          //     branches agree where they overlap.
          if (exportVtxResidual_) {
            Jpsi_jacVtx.resize(nparsfinal);
            if (Jpsi_vtxfree) {
              Map<Matrix<float, 1, Dynamic, RowMajor>>(Jpsi_jacVtx.data(), 1, nparsfinal) =
                  dxdparms.col(6).transpose().cast<float>();
            } else if (wvtxinf.size() == Jfinal.rows()) {
              Map<Matrix<float, 1, Dynamic, RowMajor>>(Jpsi_jacVtx.data(), 1, nparsfinal) =
                  (-(wvtxinf.transpose() * Jfinal)).cast<float>();
            } else {
              std::fill(Jpsi_jacVtx.begin(), Jpsi_jacVtx.end(), 0.f);
            }
          }


          if (false) {
            std::cout << "Muplus pt eta phi = " << Muplus_pt << "  " << Muplus_eta << " " << Muplus_phi << std::endl;
            std::cout << "Muminus pt eta phi = " << Muminus_pt << "  " << Muminus_eta << " " << Muminus_phi << std::endl;

            std::cout << "Muplus_jacref:\n" << Map<Matrix<float, 3, Dynamic, RowMajor>>(Muplus_jacRef.data(), 3, nparsfinal) << std::endl;
            std::cout << "Muminus_jacref:\n" << Map<Matrix<float, 3, Dynamic, RowMajor>>(Muminus_jacRef.data(), 3, nparsfinal) << std::endl;
          }
        }

      }
      
      if (!valid) {
        continue;
      }
      ++fitSucceeded_;

// std::cout << "gradfull rows cols " << gradfull.rows() << " " << gradfull.cols() << "nstateparms = " << nstateparms << std::endl;
    
// auto const& dchisqdx = gradfull.head(nstateparms);
// auto const& dchisqdparms = gradfull.tail(npars);
//
// auto const& d2chisqdx2 = hessfull.topLeftCorner(nstateparms, nstateparms);
// auto const& d2chisqdxdparms = hessfull.topRightCorner(nstateparms, npars);
// auto const& d2chisqdparms2 = hessfull.bottomRightCorner(npars, npars);
//
// std::unordered_map<unsigned int, unsigned int> idxmap;
//
// globalidxvfinal.clear();
// globalidxvfinal.reserve(npars);
// idxmap.reserve(npars);
//
// for (unsigned int idx : globalidxv) {
// if (!idxmap.count(idx)) {
// idxmap[idx] = globalidxvfinal.size();
// globalidxvfinal.push_back(idx);
// }
// }
//
      const unsigned int nparsfinal = globalidxvfinal.size();
//
// VectorXd dchisqdparmsfinal = VectorXd::Zero(nparsfinal);
// MatrixXd d2chisqdxdparmsfinal = MatrixXd::Zero(nstateparms, nparsfinal);
// MatrixXd d2chisqdparms2final = MatrixXd::Zero(nparsfinal, nparsfinal);
//
// for (unsigned int i = 0; i < npars; ++i) {
// const unsigned int iidx = idxmap.at(globalidxv[i]);
// dchisqdparmsfinal[iidx] += dchisqdparms[i];
// d2chisqdxdparmsfinal.col(iidx) += d2chisqdxdparms.col(i);
// for (unsigned int j = 0; j < npars; ++j) {
// const unsigned int jidx = idxmap.at(globalidxv[j]);
// d2chisqdparms2final(iidx, jidx) += d2chisqdparms2(i, j);
// }
// }
//
// dxdparms = -Cinvd.solve(d2chisqdxdparmsfinal).transpose();
//
// // grad = dchisqdparmsfinal + dxdparms*dchisqdx;
// grad = dchisqdparmsfinal + d2chisqdxdparmsfinal.transpose()*dxfull;
// hess = d2chisqdparms2final + dxdparms*d2chisqdxdparmsfinal;
      
  // if (debugprintout_) {
  // std::cout << "dxrefdparms" << std::endl;
  // std::cout << dxdparms.leftCols<5>() << std::endl;
  // }
      
// grad = dchisqdparms + dxdparms*dchisqdx;
      
// std::cout << "dchisqdparms" << std::endl;
// std::cout << dchisqdparms.transpose() << std::endl;
// std::cout << "dxdparms*dchisqdx" << std::endl;
// std::cout << (dxdparms*dchisqdx).transpose() << std::endl;
// std::cout << "grad" << std::endl;
// std::cout << grad.transpose() << std::endl;
      //TODO check the simplification
  // hess = d2chisqdparms2 + 2.*dxdparms*d2chisqdxdparms + dxdparms*d2chisqdx2*dxdparms.transpose();
// hess = d2chisqdparms2 + dxdparms*d2chisqdxdparms;
  
// for (unsigned int iparm = 0; iparm < npars; ++iparm) {
// if (detidparmsrev[globalidxv[iparm]].first != 7) {
// hess.row(iparm) *= 0.;
// hess.col(iparm) *= 0.;
// hess(iparm, iparm) = 1e6;
// }
// }
      
// SelfAdjointEigenSolver<MatrixXd> es(hess, EigenvaluesOnly);
// const double condition = es.eigenvalues()[nstateparms-1]/es.eigenvalues()[0];
// std::cout << "hess eigenvalues:" << std::endl;
// std::cout << es.eigenvalues().transpose() << std::endl;
// std::cout << "condition: " << condition << std::endl;
      
// std::cout << "hess diagonal:" << std::endl;
// std::cout << hess.diagonal().transpose() << std::endl;
// 
// assert(es.eigenvalues()[0] > -1e-5);
// assert(hess.diagonal().minCoeff() > 0.);
      
      nParms = nparsfinal;

      gradv.clear();
      gradv.resize(nparsfinal,0.);

      if (fillTrackTree_ && (fillGrads_ || fillGradsFactored_)) {
        tree->SetBranchAddress("gradv", gradv.data());
      }
      
      //eigen representation of the underlying vector storage
      Map<VectorXf> gradout(gradv.data(), nparsfinal);

      gradout = grad.cast<float>();
      
        

      gradmax = 0.;
      for (unsigned int i=0; i<nparsfinal; ++i) {
        const float absval = std::abs(grad[i]);
        if (absval>gradmax) {
          gradmax = absval;
        }      
      }
      
      
      hessmax = 0.;
      for (unsigned int i=0; i<nparsfinal; ++i) {
        for (unsigned int j=i; j<nparsfinal; ++j) {
          const unsigned int iidx = globalidxvfinal[i];
          const unsigned int jidx = globalidxvfinal[j];
          
          const float absval = std::abs(hess(i,j));
          if (absval>hessmax) {
            hessmax = absval;
          }
          
        }
        
      }
        
      //fill packed hessian and indices
      const unsigned int nsym = nparsfinal*(1+nparsfinal)/2;
      hesspackedv.clear();
      hesspackedv.resize(nsym, 0.);

      nSym = nsym;
      if (fillTrackTree_ && fillGrads_) {
        tree->SetBranchAddress("hesspackedv", hesspackedv.data());
      }

      Map<VectorXf> hesspacked(hesspackedv.data(), nsym);
      const Map<const VectorXu> globalidx(globalidxvfinal.data(), nparsfinal);

      unsigned int packedidx = 0;
      for (unsigned int ipar = 0; ipar < nparsfinal; ++ipar) {
        const unsigned int segmentsize = nparsfinal - ipar;
        hesspacked.segment(packedidx, segmentsize) = hess.block<1, Dynamic>(ipar, ipar, 1, segmentsize).cast<float>();
        packedidx += segmentsize;
      }

      // Factored (low-rank) Hessian storage: hess = 2 J^T R J with
      // rank(R) = rank(Vinv) - nstatefree = ndof exactly -- the ndof
      // member counts the genuinely weighted residual rows (the strip
      // second coordinates carry exactly zero weight) minus the free
      // state parameters, plus the bs/pointing/vtx constraint rows and
      // the mass row on the constrained pass, matching the fit
      // configuration by construction. Keeping exactly the top-ndof
      // eigenmodes (rather than thresholding on the numerics) stores
      // nRank*nParms floats instead of nParms*(nParms+1)/2; everything
      // beyond index ndof is numerical noise of the double-precision
      // profiling, whose tail (observed up to ~1e-4 relative on rare
      // candidates) overlaps the smallest genuine modes (down to
      // ~4e-6 relative), so no eigenvalue cut separates the two -- the
      // count does. hessrankgap monitors the truncation boundary.
      // Convention: H = B^T B, B row-major (nRank x nParms), row k =
      // sqrt(lambda_k) * v_k^T.
      if (fillGradsFactored_) {
        MatrixXd hessmean = hess;
        if (nHessVar > 0 && hessvarpackedv.size()
                                == std::size_t(nHessVar) * (nHessVar + 1) / 2) {
          std::size_t ip = 0;
          for (unsigned int i = 0; i < nHessVar; ++i) {
            for (unsigned int j = i; j < nHessVar; ++j, ++ip) {
              const double v = hessvarpackedv[ip];
              hessmean(hessvaridxv[i], hessvaridxv[j]) -= v;
              if (i != j) {
                hessmean(hessvaridxv[j], hessvaridxv[i]) -= v;
              }
            }
          }
        }
        const SelfAdjointEigenSolver<MatrixXd> eshess(hessmean);
        const VectorXd& eigvals = eshess.eigenvalues();  // ascending

        // WHAT THE FACTORED FORM CANNOT DO.  `B^T B` is exact for the MEAN
        // Hessian `2 J^T R J`, whose rank is exactly ndof.  The log-det block
        // `tr(dV_i R dV_j R)` is a Gram matrix of the ncons x ncons objects
        // `R^1/2 dV_i R^1/2`, so its rank is the NUMBER of variance
        // parameters: it does not compress, and carrying it as extra rows of
        // B costs nvar*nParms floats -- measured at +57 % of a production
        // candidate with family 15 alone and +429 % with 8-11.  Its own
        // packed triangle costs nvar*(nvar+1)/2, ~20x less.  So `hessfactorv`
        // factors the MEAN block only, at `nRank = ndof`, and
        // the variance block rides in `hessvaridxv`/`hessvarpackedv`:
        //
        //     hess = B^T B + scatter(hessvarpackedv on hessvaridxv)
        //
        // `hesspackedv`, when written, is COMPLETE and needs no addition.
        const unsigned int nrank = std::min(ndof, nparsfinal);

        double keptmass = 0.;
        double droppedmass = 0.;
        for (unsigned int ieig = 0; ieig < nparsfinal; ++ieig) {
          const double lambda = eigvals(ieig);
          if (ieig >= nparsfinal - nrank) {
            keptmass += lambda;
          }
          else if (lambda > 0.) {
            droppedmass += lambda;
          }
        }

        // DEFENCE IN DEPTH.  `nrank == 0` (ndof == 0, or an empty parameter
        // list) makes `nparsfinal - nrank` one past the end of `eigvals` and
        // `nparsfinal - nrank - 1` wrap; the candidate-level gate at the ndof
        // computation already rejects that fit, so this branch is unreachable
        // in a healthy job -- but an out-of-range Eigen index aborts the
        // PROCESS, and a production is not the place to find out that some
        // other path can reach it. Bit-identical for every nrank > 0.
        const double lambdakeptmin = nrank > 0 ? eigvals(nparsfinal - nrank) : 0.;
        hessrankgap = (nrank > 0 && nrank < nparsfinal && lambdakeptmin > 0.)
                          ? std::max(0., eigvals(nparsfinal - nrank - 1)) / lambdakeptmin
                          : 0.;

        nRank = nrank;
        nFactor = nrank*nparsfinal;
        hessdroppedmass = keptmass > 0. ? droppedmass/keptmass : 0.;

        hessfactorv.clear();
        hessfactorv.resize(nFactor, 0.);
        if (fillTrackTree_) {
          tree->SetBranchAddress("hessfactorv", hessfactorv.data());
        }

        // rows ordered by decreasing eigenvalue; the max(0,.) guards the
        // fixed-count boundary against a numerically negative eigenvalue
        // in a (hypothetical) rank-deficient candidate
        Map<Matrix<float, Dynamic, Dynamic, RowMajor>> hessfactor(hessfactorv.data(), nrank, nparsfinal);
        for (unsigned int irank = 0; irank < nrank; ++irank) {
          const unsigned int ieig = nparsfinal - 1 - irank;
          hessfactor.row(irank) = (std::sqrt(std::max(0., eigvals(ieig)))*eshess.eigenvectors().col(ieig)).transpose().cast<float>();
        }
      }

      
// assert(globalidxvfinal.size() == (2*Muplus_nhits + 2*Muminus_nhits + 2*Muplus_nvalid + 2*Muminus_nvalid + Muplus_nvalidpixel + Muminus_nvalidpixel));

// hessv.resize(nparsfinal*nparsfinal);
// Map<Matrix<float, Dynamic, Dynamic, RowMajor>>(hessv.data(), nparsfinal, nparsfinal) = hess.cast<float>();

      // NanoAOD path: capture this candidate's refit result into the ValueMap
      // accumulators, positionally by srcCandidates index. This point is after
      // the icons loop (once per candidate); the Jpsi_*/Mu* members hold the
      // unconstrained (icons==0) result. Reached only when the candidate
      // completed the fit -- skipped/failed candidates keep their sentinel.
      // Kinematics always; the global-fit payload only when fillGradsFactored_
      // populated it.
      if (doVM && candCollIdx >= 0) {
        vmCorMassV[candCollIdx]    = Jpsi_mass;
        vmCorMassErrV[candCollIdx] = Jpsi_sigmamass;
        vmCorPtV[candCollIdx]      = Jpsi_pt;
        vmCorEtaV[candCollIdx]     = Jpsi_eta;
        vmCorPhiV[candCollIdx]     = Jpsi_phi;
        vmMuPlusPtV[candCollIdx]   = Muplus_pt;
        vmMuPlusEtaV[candCollIdx]  = Muplus_eta;
        vmMuPlusPhiV[candCollIdx]  = Muplus_phi;
        vmMuMinusPtV[candCollIdx]  = Muminus_pt;
        vmMuMinusEtaV[candCollIdx] = Muminus_eta;
        vmMuMinusPhiV[candCollIdx] = Muminus_phi;
        vmEdmvalV[candCollIdx]     = edmval;
        if (fillGradsFactored_) {
          vmGlobalIdxsV[candCollIdx].assign(globalidxvfinal.begin(), globalidxvfinal.end());
          vmJacRefMuPlusV[candCollIdx]  = Muplus_jacRef;
          vmJacRefMuMinusV[candCollIdx] = Muminus_jacRef;
          vmJacMassV[candCollIdx]       = Jpsi_jacMass;
          vmHessFactorV[candCollIdx]    = hessfactorv;
        }
      }

      if (fillTrackTree_) {
        tree->Fill();
      }
      

      

      
// const Matrix3d covvtx = Cinvd.solve(MatrixXd::Identity(nstateparms,nstateparms)).topLeftCorner<3,3>();
// 
// const double covqop0 = Cinvd.solve(MatrixXd::Identity(nstateparms,nstateparms))(trackstateidxarr[0], trackstateidxarr[0]);
// const double covqop1 = Cinvd.solve(MatrixXd::Identity(nstateparms,nstateparms))(trackstateidxarr[1], trackstateidxarr[1]);
// 
// const double covqop0kin = outparts[0]->currentState().freeTrajectoryState().curvilinearError().matrix()(0,0);
// const double covqop1kin = outparts[1]->currentState().freeTrajectoryState().curvilinearError().matrix()(0,0);
// 
// // Matrix<double, 1, 1> covmass = 2.*massjac*Cinvd.solve(MatrixXd::Identity(nstateparms,nstateparms))*massjac.transpose();
// 
// // const VectorXd cinvrow0 = Cinvd.solve(MatrixXd::Identity(nstateparms,nstateparms)).row(0).head(nstateparms);
// // 
// std::cout << "kinfit covariance:" << std::endl;
// std::cout << dimu_vertex->error().matrix() << std::endl;
// 
// std::cout << "GBL covariance:" << std::endl;
// std::cout << 2.*covvtx << std::endl;
// 
// std::cout << "kinfit qop0 covariance:" << std::endl;
// std::cout << covqop0kin << std::endl;
// 
// std::cout << "GBL qop0 covariance:" << std::endl;
// std::cout << 2.*covqop0 << std::endl;
// 
// std::cout << "kinfit qop1 covariance:" << std::endl;
// std::cout << covqop1kin << std::endl;
// 
// std::cout << "GBL qop1 covariance:" << std::endl;
// std::cout << 2.*covqop1 << std::endl;
      
// std::cout << "dqop0 beamline" << std::endl;
// std::cout << dxfull[trackstateidxarr[0]] << std::endl;
// std::cout << "dqop0 first layer" << std::endl;
// std::cout << dxfull[trackstateidxarr[0]+3] << std::endl;
// std::cout << "dqop0 second layer" << std::endl;
// std::cout << dxfull[trackstateidxarr[0]+6] << std::endl;
// 
// std::cout << "dqop1 beamline" << std::endl;
// std::cout << dxfull[trackstateidxarr[1]] << std::endl;
// std::cout << "dqop1 first layer" << std::endl;
// std::cout << dxfull[trackstateidxarr[1]+3] << std::endl;
// std::cout << "dqop1 second layer" << std::endl;
// std::cout << dxfull[trackstateidxarr[1]+6] << std::endl;
// 
// std::cout << "sigmam" << std::endl;
// std::cout << std::sqrt(covmass[0]) << std::endl;

// 
// std::cout << "cinvrow0" << std::endl;
// std::cout << cinvrow0 << std::endl;
      
      //TODO restore statejac stuff
// dxstate = statejac*dxfull;
// const Vector5d dxRef = dxstate.head<5>();
// const Matrix5d Cinner = (statejac*Cinvd.solve(MatrixXd::Identity(nstateparms,nstateparms))*statejac.transpose()).topLeftCorner<5,5>();
      
      //TODO fill outputs
      
    }
  }

  // NanoAOD path: put the per-candidate ValueMaps, keyed to the srcCandidates
  // collection (one entry per input candidate). Products are always put when
  // produceValueMaps_ is set (empty maps if there were no candidates).
  if (produceValueMaps_) {
    auto putF = [&](edm::EDPutTokenT<edm::ValueMap<float>>& tok, std::vector<float>& v) {
      edm::ValueMap<float> m;
      if (doVM) {
        edm::ValueMap<float>::Filler f(m);
        f.insert(vmCandH, v.begin(), v.end());
        f.fill();
      }
      iEvent.emplace(tok, std::move(m));
    };
    putF(vmCorMass_, vmCorMassV);
    putF(vmCorMassErr_, vmCorMassErrV);
    putF(vmCorPt_, vmCorPtV);
    putF(vmCorEta_, vmCorEtaV);
    putF(vmCorPhi_, vmCorPhiV);
    putF(vmMuPlusPt_, vmMuPlusPtV);
    putF(vmMuPlusEta_, vmMuPlusEtaV);
    putF(vmMuPlusPhi_, vmMuPlusPhiV);
    putF(vmMuMinusPt_, vmMuMinusPtV);
    putF(vmMuMinusEta_, vmMuMinusEtaV);
    putF(vmMuMinusPhi_, vmMuMinusPhiV);
    putF(vmEdmval_, vmEdmvalV);
    if (fillGradsFactored_) {
      {
        edm::ValueMap<std::vector<int>> m;
        if (doVM) {
          edm::ValueMap<std::vector<int>>::Filler f(m);
          f.insert(vmCandH, vmGlobalIdxsV.begin(), vmGlobalIdxsV.end());
          f.fill();
        }
        iEvent.emplace(vmGlobalIdxs_, std::move(m));
      }
      auto putVF = [&](edm::EDPutTokenT<edm::ValueMap<std::vector<float>>>& tok,
                       std::vector<std::vector<float>>& v) {
        edm::ValueMap<std::vector<float>> m;
        if (doVM) {
          edm::ValueMap<std::vector<float>>::Filler f(m);
          f.insert(vmCandH, v.begin(), v.end());
          f.fill();
        }
        iEvent.emplace(tok, std::move(m));
      };
      putVF(vmJacRefMuPlus_, vmJacRefMuPlusV);
      putVF(vmJacRefMuMinus_, vmJacRefMuMinusV);
      putVF(vmJacMass_, vmJacMassV);
      putVF(vmHessFactor_, vmHessFactorV);
    }
  }
}


DEFINE_FWK_MODULE(ResidualGlobalCorrectionMakerTwoTrackG4e);
