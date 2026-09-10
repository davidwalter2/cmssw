#include "ResidualGlobalCorrectionMakerBase.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "TrackPropagation/Geant4e/interface/Geant4ePropagator.h"
#include "TrackPropagation/Geant4e/interface/G4UniversalFluctuationForExtrapolator.hh"
#include "TrackPropagation/Geant4e/interface/MaterialGroupModel.h"
#include "Analysis/HitAnalyzer/interface/ParticleProperties.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


#include "DataFormats/PatCandidates/interface/Muon.h"

#include "Math/Vector4Dfwd.h"
#include "Math/Vector4D.h"

#include <algorithm>
#include <array>
#include <Eigen/Sparse>
#include <Eigen/Cholesky>

#include <iomanip>
#include <limits>
#include <chrono>
#include <iostream>

#include "TRandom.h"

#include "Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"

// Hit-resolution study: the strip CPE is queried DIRECTLY (in addition to
// going through the cloner) so its own independent variable -- the
// drift-included projected path in strip-pitch units -- can be exported per
// hit. Reconstructing it offline would need the per-module Lorentz drift.
#include "RecoLocalTracker/SiStripRecHitConverter/interface/StripCPE.h"
#include "RecoLocalTracker/Records/interface/TkStripCPERecord.h"
#include "RecoLocalTracker/ClusterParameterEstimator/interface/StripClusterParameterEstimator.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"



class ResidualGlobalCorrectionMakerG4e : public ResidualGlobalCorrectionMakerBase
{
public:
  ResidualGlobalCorrectionMakerG4e(const edm::ParameterSet &);
  ~ResidualGlobalCorrectionMakerG4e();

// static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  
  virtual void beginStream(edm::StreamID) override;
  virtual void produce(edm::Event &, const edm::EventSetup &) override;

  edm::EDGetTokenT<edm::Association<reco::TrackExtraCollection>> inputAssoc_;
  
  bool trackHighPurity = false;

  float muonPt;
  bool muonLoose;
  bool muonMedium;
  bool muonTight;
  bool muonIsPF;
  bool muonIsTracker;
  bool muonIsGlobal;
  bool muonIsStandalone;
  bool muonInnerTrackBest;
  
  bool trackExtraAssoc;
  
  std::vector<float> dEpred;
  std::vector<float> dE;
  std::vector<float> sigmadE;
  std::vector<float> Epred;
  std::vector<float> E;
  
  float outPt;
  float outEta;
  float outPhi;
  
  float outPtStart;
  float outEtaStart;
  float outPhiStart;

  float trackQopErr = 0.;

  edm::EDPutTokenT<edm::ValueMap<float>> outputCorPt_;
  edm::EDPutTokenT<edm::ValueMap<float>> outputCorEta_;
  edm::EDPutTokenT<edm::ValueMap<float>> outputCorPhi_;
  edm::EDPutTokenT<edm::ValueMap<int>> outputCorCharge_;
  edm::EDPutTokenT<edm::ValueMap<float>> outputCorDxy_;
  edm::EDPutTokenT<edm::ValueMap<float>> outputCorDz_;
  edm::EDPutTokenT<edm::ValueMap<float>> outputEdmval_;
  edm::EDPutTokenT<edm::ValueMap<int>> outputNValidHits_;
  edm::EDPutTokenT<edm::ValueMap<int>> outputNValidPixelHits_;

  edm::EDPutTokenT<edm::ValueMap<std::vector<int>>> outputGlobalIdxs_;

  edm::EDPutTokenT<edm::ValueMap<std::vector<float>>> outputJacRef_;
  edm::EDPutTokenT<edm::ValueMap<std::vector<float>>> outputMomCov_;

  // Optional refit-track output. The ValueMaps above publish only a 3x3
  // momentum block, and only when a muon association is configured, so they
  // cannot seed a TransientTrack. A downstream constrained B-vertex fit
  // (KinematicConstrainedVertexFitter, as Bmm5 does it) needs full 5-parameter
  // tracks with their 5x5 covariance, which is what this emits. Off by
  // default, so nominal behaviour is untouched.
  //
  // The collection is positionally aligned with the input tracks: entry i is
  // the refit of input track i. Entries whose fit failed keep a copy of the
  // input track and are flagged refitOk = 0, so a leg that was never refit
  // stays visible instead of silently masquerading as corrected.
  bool emitRefitTracks_ = false;
  edm::EDPutTokenT<reco::TrackCollection> outputRefitTracks_;
  edm::EDPutTokenT<edm::ValueMap<int>> outputRefitOk_;

  // Optional sanity cut on the refit momentum uncertainty. The NaN guard
  // below still lets through a finite-but-absurd tail (worst seen: ptErr of
  // 92 GeV on a ~1 GeV track), which would destabilise a downstream vertex
  // fit. A track above this relative threshold is emitted as the input copy
  // with refitOk = 0. Negative (the default) disables the cut entirely, so
  // nominal behaviour is unchanged and the tail stays available for study.
  double refitMaxRelPtErr_ = -1.;

  edm::ESGetToken<TransientTrackingRecHitBuilder, TransientRecHitRecord> ttrhToken_;
  edm::ESGetToken<Propagator, TrackingComponentsRecord> g4ePropToken_;
  edm::ESGetToken<StripClusterParameterEstimator, TkStripCPERecord> stripCPEToken_;

  SiStripClusterInfo siStripClusterInfo_;

  // Track mass + Geant4 particle base name used in the refit. Defaults
  // preserve the legacy J/psi/Upsilon (muon) configuration; channel cfis
  // can override to study e.g. pions, kaons, protons through the
  // single-track refit.
  double trackMass_;
  std::string trackParticleName_;

  // Per-stream fit-outcome accounting, printed from the destructor.
  // attempted = tracks entering the iterative refit (post-selection);
  // succeeded = fits that survived all iterations (tree/valuemap filled).
  mutable unsigned long long fitAttempted_ = 0ULL;
  mutable unsigned long long fitSucceeded_ = 0ULL;
  mutable unsigned long long fitFailProp_ = 0ULL;       // Geant4e propagation failed
  mutable unsigned long long fitFailHitUpdate_ = 0ULL;  // CPE re-evaluation (cloner) invalid
  mutable unsigned long long fitFailNaN_ = 0ULL;        // NaN/inf parameter update
  mutable unsigned long long fitStepClamped_ = 0ULL;    // fits with >=1 trust-region-clamped GN step
  mutable unsigned long long stepClampEvents_ = 0ULL;   // individual clamped GN steps
  mutable unsigned long long fitStepBacktracked_ = 0ULL;// fits with >=1 chi2 backtrack
  mutable unsigned long long stepBacktrackEvents_ = 0ULL;// individual chi2 step halvings
  mutable unsigned long long stepPrints_ = 0ULL;        // printouts emitted (rate limit)
  mutable unsigned long long fitChargeFlipAllowed_ = 0ULL;  // GN steps permitted to cross q/p=0 (high-p)
  mutable unsigned long long fitChargeHypTaken_ = 0ULL;  // tracks where the two-hyp fit kept the opposite charge

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

  // Convergence knobs (see ctor). Defaults reproduce the baseline.
  unsigned int nIters_ = 10;
  double edmConvergence_ = 1.e-5;
  unsigned int gnDampAfter_ = 0;   // 0 = damping off
  double gnDampFactor_ = 0.5;
  // Momentum above which the Gauss-Newton step is ALLOWED to cross q/p = 0
  // (a genuine charge flip). Crossing q/p = 0 is p -> inf, harmless at high p
  // where the seed charge is ambiguous; the sign-flip clamp is only a
  // divergence guard for stiffer (lower-p) tracks. Default 1e9 GeV = never
  // allow = legacy behaviour. The momentum floor below always applies.
  double allowChargeFlipAboveP_ = 1.e9;
  // Gauss-Newton momentum floor [GeV]: the updated |p| a step is allowed to
  // reach. Its ONLY job is to keep the state out of the propagator's
  // refusal region (Geant4ePropagator.PropagationPtotLimit), so it must
  // track that limit and nothing else. It was a hard-coded 2.0 GeV until
  // 2026-09-04, chosen when the propagation limit was 1.0 GeV; the drivers
  // have since lowered the limit to 0.2 GeV, and a 2 GeV clamp then PINS
  // every genuinely soft track at 2 GeV (momentum-high, chi2/ndof >> 1).
  // Configurable via the `clampMomentumFloor` cfi parameter, whose driver
  // default is derived from the propagation limit. Same name/semantics as
  // the two-track and N-track makers.
  double clampMomentumFloor_ = 2.0;
  // Relative Gauss-Newton step damping (2026-09-05). Per iteration a track's
  // momentum may change by at most this factor (default 2: p may at most halve
  // or double). Implemented as the effective floor max(clampMomentumFloor_,
  // p_ref/f) plus the symmetric upward cap p_ref*f, so the bound is ALWAYS
  // strictly below/above p_ref and no track is ever pinned at a fixed
  // momentum, whatever its true momentum is. <= 1 restores the legacy
  // absolute-floor-only clamp bit-identically.
  double maxMomentumStepFactor_ = 2.0;
  // chi2-based (Armijo) retroactive backtracking. The chi2 assembled at the
  // top of iteration k is the REALIZED chi2 of the step taken at k-1; if it
  // fails the sufficient-decrease test the previous linearization is restored,
  // that step is halved and the iteration is redone. Costs nothing when it
  // does not fire (the chi2 is assembled anyway).
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
  // tolerance: measured 2026-09-05, the CVH/GBL iteration does NOT
  // monotonically decrease r^T Vinv r -- the realized chi2 drifts UP by
  // ~0.3-0.5 per iteration even at 1/16 of the step (the model's predicted
  // decrease is never realized because every iteration re-propagates and
  // re-linearizes). A tolerance of 1e-3 therefore turns the test into a
  // permanent step-halver: 57 % of gun candidates backtracked, 15.6 halvings
  // each, 6x the propagation cost, for no change in the result. At 1.0
  // (the chi2 may not more than DOUBLE in one iteration) the test becomes a
  // pure DIVERGENCE TRAP: +0.5 % propagation on the gun ditrack smoke,
  // +0.2 % single track, fit output at the noise level. Scan in NOTES.md.
  double armijoSlack_ = 1.0;
  // Cap on the per-job "GN step clamped/backtracked" printouts. The relative
  // damping legitimately fires far more often than the absolute floor did, and
  // an uncapped print would bury the log (and the counters are exact anyway).
  unsigned int stepPrintLimit_ = 200;
  // Two-hypothesis charge resolution in one pass: when the nominal fit's clamp
  // catches a q/p sign crossing (nChargeFlipProtect>0), re-fit the opposite
  // charge and keep the lower-chi2 result. Inert for well-measured (low/mod p)
  // tracks whose fit never tries to cross q/p=0. Default on.
  bool twoHypothesisCharge_ = true;
  // Cosmic seeding: default (false) seeds from the inner (near-beamline) state
  // and propagates OUTWARD, so the upper half of a down-going cosmic propagates
  // BACKWARD (against the muon). When true, seed from the muon-ENTRY end (higher
  // y) and iterate hits along the muon (downward) -> the whole leg propagates
  // FORWARD. Moves the report point to the entry end; validates the backward
  // energy-loss handling (forward-only vs mixed on the same track).
  bool cosmicSeedFromEntry_ = false;
  // Force the outside-in cosmic handling regardless of track algo (needed for
  // split-track legs from cosmicTrackSplitter, whose algo 'cosmic' is not
  // recognised as ctf/cosmics). Default false.
  bool forceCosmic_ = false;
  // Two-hypothesis charge test: seed the fit with sign(q) * seedChargeSign_
  // (applied consistently to the G4 particle and the initial state), so the
  // fit converges to the OPPOSITE-charge minimum (the clamp keeps it there).
  // Run +1 and -1 and compare chi2 per track to identify genuine mis-ID
  // (opposite hypothesis fits better) vs correct seed (opposite diverges).
  // Default +1 = nominal seed.
  int seedChargeSign_ = 1;

  // Global material model validation knobs (Phase A; the model itself and
  // the parmtype-15 registration live in the base class). With a model
  // loaded, per-leg per-group dxi columns are accumulated by the
  // propagator and checked inline against the integrated dxi column (V2
  // sum identity); materialFDGroup >= 0 additionally runs a one-shot FD
  // closure (V1): the first leg crossing that group is re-propagated with
  // k_g += eps and d(q/p)/dk compared against the analytic column.
  int materialFDGroup_ = -1;
  double materialFDEps_ = 1e-3;
  mutable std::vector<std::pair<int, Eigen::Matrix<double, 5, 1>>> groupJacs_;
  // Per-group PROCESS NOISE of the last propagation (the width counterpart of
  // groupJacs_'s mean). Filled only under doRes + the global material model;
  // see the parmtype-15 dV registration.
  mutable std::vector<std::pair<int, Eigen::Matrix<double, 5, 5>>> groupQs_;
  // per-step field-mode columns from the propagator (reused buffer)
  mutable std::vector<Eigen::Matrix<double, 5, 1>> modeJacs_;
  mutable double v2MaxRelDiff_ = 0.;
  mutable unsigned long long v2Checks_ = 0ULL;
  mutable bool fdMatDone_ = false;

  // Per-iteration debug dump (mirrors the two-track maker): vector branches
  // recording the full chi2/EDM trajectory of each fit. Only booked/filled
  // when debugPerIterDump=true.
  bool debugPerIterDump_ = false;
  std::vector<double> chisqval_iter;
  std::vector<double> edmval_iter;
  std::vector<double> edmvalref_iter;
  std::vector<double> deltachisqval_iter;

  // Kink finder (decay-in-flight score test). For each material step the
  // alternative hypothesis of an unconstrained offset (q/p, dx/dz, dy/dz)
  // of the propagated state at the step end -- a decay kink plus momentum
  // step -- is scored against the converged fit without refitting (see the
  // scan block after the iteration loop). Components are in the local frame
  // of the step-end surface when dolocalupdate, curvilinear otherwise.
  bool doKinkFinder_ = false;
  // Synthetic-kink injection for the closure test: add the given offset to
  // the material residual dx0 at material step kinkInjectLayer (mimicking a
  // decay of that size there); the scan must recover it. -1 = off.
  int kinkInjectLayer_ = -1;
  double kinkInjectDqop_ = 0.;
  double kinkInjectDxdz_ = 0.;
  double kinkInjectDydz_ = 0.;
  // branch buffers (per material step; hold the winning charge hypothesis)
  std::vector<float> kinkDchisq;
  std::vector<float> kinkDchisqAngle;
  std::vector<float> kinkDchisqQop;
  std::vector<float> kinkDqop;
  std::vector<float> kinkDxdz;
  std::vector<float> kinkDydz;
  std::vector<float> kinkGlobalR;
  std::vector<float> kinkGlobalZ;
  float kinkMax = -1.f;
  int kinkMaxLayer = -1;
  // per-iteration cache: row offset of each material block in the constraint
  // vector, and the step-end global position (rebuilt every iteration; the
  // final iteration's layout is what the scan consumes)
  std::vector<unsigned int> kinkConsIdx;
  std::vector<float> kinkStepR;
  std::vector<float> kinkStepZ;

  // Geant4 decay/interaction truth for the gen-matched particle
  // (doSimDecayTruth_). Raw sim quantities only -- the physics
  // classification (decay vs nuclear interaction vs delta ray, kink angle,
  // daughter momentum fraction) is done offline, where the parent direction
  // at the vertex can be obtained by helix propagation.
  int simTrkFound = 0;
  int simTrkPdgId = 0;
  float simTrkP = -99.f;
  float simTrkPt = -99.f;
  float simTrkEta = -99.f;
  float simTrkPhi = -99.f;
  float simTrkVtxX = -99.f;
  float simTrkVtxY = -99.f;
  float simTrkVtxZ = -99.f;
  // one entry per SimVertex whose parent is the matched SimTrack
  std::vector<float> simVtxX;
  std::vector<float> simVtxY;
  std::vector<float> simVtxZ;
  std::vector<int> simVtxProcType;
  std::vector<int> simVtxNDau;
  // daughters, flattened; simDauVtx indexes into the simVtx* vectors
  std::vector<int> simDauVtx;
  std::vector<int> simDauPdgId;
  std::vector<float> simDauPx;
  std::vector<float> simDauPy;
  std::vector<float> simDauPz;
};

ResidualGlobalCorrectionMakerG4e::~ResidualGlobalCorrectionMakerG4e() {
  if (fitAttempted_ > 0ULL) {
    const unsigned long long failTotal = fitAttempted_ - fitSucceeded_;
    std::cout << "ResidualGlobalCorrectionMakerG4e fit summary"
              << "  attempted=" << fitAttempted_
              << "  succeeded=" << fitSucceeded_
              << "  failed=" << failTotal
              << " (" << (100. * failTotal / fitAttempted_) << "%)"
              << "  fail[prop]=" << fitFailProp_
              << "  fail[hitupdate]=" << fitFailHitUpdate_
              << "  fail[nan]=" << fitFailNaN_
              << "  clamped[step]=" << fitStepClamped_
              << "  clampevents[step]=" << stepClampEvents_
              << "  backtracked[step]=" << fitStepBacktracked_
              << "  backtrackevents[step]=" << stepBacktrackEvents_
              << "  chargeflip[allowed]=" << fitChargeFlipAllowed_
              << "  chargehyp[flipped]=" << fitChargeHypTaken_
              << std::endl;
    if (pixHitsSeen_ > 0ULL) {
      std::cout << "ResidualGlobalCorrectionMakerG4e pixel hit-quality summary"
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
      // Pathology-class combination table (only non-empty bins).
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
    if (v2Checks_ > 0ULL) {
      std::cout << "ResidualGlobalCorrectionMakerG4e material-group V2 (sum == dxi column)"
                << "  checks=" << v2Checks_
                << "  maxRelDiff=" << v2MaxRelDiff_
                << std::endl;
    }
  }
}


ResidualGlobalCorrectionMakerG4e::ResidualGlobalCorrectionMakerG4e(const edm::ParameterSet &iConfig)
    : ResidualGlobalCorrectionMakerBase(iConfig),
      ttrhToken_(esConsumes(edm::ESInputTag("", "WithAngleAndTemplate"))),
      g4ePropToken_(esConsumes(edm::ESInputTag("", "Geant4ePropagator"))),
      stripCPEToken_(esConsumes(edm::ESInputTag("", "StripCPEfromTrackAngle"))),
      siStripClusterInfo_(consumesCollector())
{

// inputAssoc_ = consumes<edm::Association<reco::TrackExtraCollection>>(edm::InputTag("muonReducedTrackExtras"));

  // Geant4 particle base name. Default = muon (legacy J/psi/Upsilon
  // configuration). Override via the cfi to e.g. "pi", "kaon", "proton"
  // to study other species. Mass is looked up from a single PDG table
  // in Analysis/HitAnalyzer/interface/ParticleProperties.h .
  trackParticleName_ = iConfig.existsAs<std::string>("trackParticleName")
      ? iConfig.getParameter<std::string>("trackParticleName") : std::string("mu");
  trackMass_ = ana_hitanalyzer::getParticleProperties(trackParticleName_).mass;

  // Convergence knobs, mirroring the two-track maker (existsAs-guarded so
  // legacy cfis keep the baseline behaviour: 10 iterations, EDM < 1e-5).
  nIters_ = iConfig.existsAs<unsigned int>("nIters")
      ? iConfig.getParameter<unsigned int>("nIters") : 10u;
  edmConvergence_ = iConfig.existsAs<double>("edmConvergence")
      ? iConfig.getParameter<double>("edmConvergence") : 1.e-5;
  debugPerIterDump_ = iConfig.existsAs<bool>("debugPerIterDump")
      ? iConfig.getParameter<bool>("debugPerIterDump") : false;
  // Gauss-Newton step damping from iteration gnDampAfter_ onward (0 = off).
  // Remedy for limit cycles between nearly-chi2-degenerate states: a
  // half-step from one cycle endpoint lands on the midpoint, collapsing
  // 2-cycles onto the true minimum. Applied only after the early
  // iterations so normal convergence (2-4 full steps) is untouched.
  gnDampAfter_ = iConfig.existsAs<unsigned int>("gnDampAfter")
      ? iConfig.getParameter<unsigned int>("gnDampAfter") : 0u;
  gnDampFactor_ = iConfig.existsAs<double>("gnDampFactor")
      ? iConfig.getParameter<double>("gnDampFactor") : 0.5;
  // Allow the GN step to cross q/p = 0 (charge flip) for tracks stiffer than
  // this momentum (default 1e9 GeV = never = legacy). Motivated by high-p
  // cosmics whose seed charge is ambiguous: freezing the seed sign can pin a
  // mis-identified track in a wrong-sign chi2 minimum.
  allowChargeFlipAboveP_ = iConfig.existsAs<double>("allowChargeFlipAboveP")
      ? iConfig.getParameter<double>("allowChargeFlipAboveP") : 1.e9;
  // Gauss-Newton momentum floor (see member comment). existsAs-guarded so a
  // cfi that does not set it keeps the historical 2.0 GeV.
  clampMomentumFloor_ = iConfig.existsAs<double>("clampMomentumFloor")
      ? iConfig.getParameter<double>("clampMomentumFloor") : 2.0;
  // Relative step damping and chi2 backtracking (see member comments).
  // maxMomentumStepFactor <= 1 => legacy absolute-floor clamp, bit-identical.
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
  // Echo it once per maker instance: the step control silently decides whether
  // soft tracks are fitted or pinned, and a job log must record what was in force.
  edm::LogPrint("ResidualGlobalCorrectionMakerG4e")
      << "[cvh] effective: clampMomentumFloor=" << clampMomentumFloor_
      << " GeV, allowChargeFlipAboveP=" << allowChargeFlipAboveP_ << " GeV"
      << ", maxMomentumStepFactor=" << maxMomentumStepFactor_
      << ", stepBacktracking=" << stepBacktracking_
      << " (fromIter=" << stepBacktrackFromIter_
      << ", maxChi2Backtrack=" << maxChi2Backtrack_
      << ", armijoC=" << armijoC_ << ", armijoSlack=" << armijoSlack_ << ")";
  // Two-hypothesis charge test: +1 = nominal seed, -1 = opposite-charge seed.
  seedChargeSign_ = iConfig.existsAs<int>("seedChargeSign")
      ? iConfig.getParameter<int>("seedChargeSign") : 1;
  // One-pass two-hypothesis charge resolution (default on).
  twoHypothesisCharge_ = iConfig.existsAs<bool>("twoHypothesisCharge")
      ? iConfig.getParameter<bool>("twoHypothesisCharge") : true;
  // Seed cosmic legs from the muon-entry end (fully forward propagation).
  cosmicSeedFromEntry_ = iConfig.existsAs<bool>("cosmicSeedFromEntry")
      ? iConfig.getParameter<bool>("cosmicSeedFromEntry") : false;
  // Force cosmic (outside-in) handling regardless of track algo.
  forceCosmic_ = iConfig.existsAs<bool>("forceCosmic")
      ? iConfig.getParameter<bool>("forceCosmic") : false;

  // Global material model validation knobs (the model itself is loaded by
  // the base class from materialGroupsFile).
  materialFDGroup_ = iConfig.existsAs<int>("materialFDGroup")
      ? iConfig.getParameter<int>("materialFDGroup") : -1;
  materialFDEps_ = iConfig.existsAs<double>("materialFDEps")
      ? iConfig.getParameter<double>("materialFDEps") : 1e-3;

  // Kink finder (existsAs-guarded so legacy cfis are untouched).
  doKinkFinder_ = iConfig.existsAs<bool>("doKinkFinder")
      ? iConfig.getParameter<bool>("doKinkFinder") : false;
  kinkInjectLayer_ = iConfig.existsAs<int>("kinkInjectLayer")
      ? iConfig.getParameter<int>("kinkInjectLayer") : -1;
  kinkInjectDqop_ = iConfig.existsAs<double>("kinkInjectDqop")
      ? iConfig.getParameter<double>("kinkInjectDqop") : 0.;
  kinkInjectDxdz_ = iConfig.existsAs<double>("kinkInjectDxdz")
      ? iConfig.getParameter<double>("kinkInjectDxdz") : 0.;
  kinkInjectDydz_ = iConfig.existsAs<double>("kinkInjectDydz")
      ? iConfig.getParameter<double>("kinkInjectDydz") : 0.;

  outputCorPt_ = produces<edm::ValueMap<float>>("corPt");
  outputCorEta_ = produces<edm::ValueMap<float>>("corEta");
  outputCorPhi_ = produces<edm::ValueMap<float>>("corPhi");
  outputCorDxy_ = produces<edm::ValueMap<float>>("corDxy");
  outputCorDz_ = produces<edm::ValueMap<float>>("corDz");
  outputCorCharge_ = produces<edm::ValueMap<int>>("corCharge");
  outputEdmval_ = produces<edm::ValueMap<float>>("edmval");
  outputNValidHits_ = produces<edm::ValueMap<int>>("nValidHits");
  outputNValidPixelHits_ = produces<edm::ValueMap<int>>("nValidPixelHits");

  outputGlobalIdxs_ = produces<edm::ValueMap<std::vector<int>>>("globalIdxs");

  outputJacRef_ = produces<edm::ValueMap<std::vector<float>>>("jacRef");
  outputMomCov_ = produces<edm::ValueMap<std::vector<float>>>("momCov");

  emitRefitTracks_ = iConfig.existsAs<bool>("emitRefitTracks")
      ? iConfig.getParameter<bool>("emitRefitTracks") : false;
  refitMaxRelPtErr_ = iConfig.existsAs<double>("refitMaxRelPtErr")
      ? iConfig.getParameter<double>("refitMaxRelPtErr") : -1.;
  if (emitRefitTracks_) {
    outputRefitTracks_ = produces<reco::TrackCollection>("refit");
    outputRefitOk_ = produces<edm::ValueMap<int>>("refitOk");
  }
}

void ResidualGlobalCorrectionMakerG4e::beginStream(edm::StreamID streamid)
{
  ResidualGlobalCorrectionMakerBase::beginStream(streamid);
  
  if (fillTrackTree_) {
    const int basketSize = 4*1024*1024;
    
    tree->Branch("trackPt", &trackPt, basketSize);
    tree->Branch("trackPtErr", &trackPtErr, basketSize);
    tree->Branch("trackEta", &trackEta, basketSize);
    tree->Branch("trackPhi", &trackPhi, basketSize);
    tree->Branch("trackCharge", &trackCharge, basketSize);
    tree->Branch("trackQopErr", &trackQopErr);
    // Per-iteration trajectory dump (debug): one entry per Gauss-Newton
    // iteration. edmvalref_iter is the reference-block EDM, i.e. the
    // actual convergence criterion.
    if (debugPerIterDump_) {
      tree->Branch("chisqval_iter",      &chisqval_iter);
      tree->Branch("edmval_iter",        &edmval_iter);
      tree->Branch("edmvalref_iter",     &edmvalref_iter);
      tree->Branch("deltachisqval_iter", &deltachisqval_iter);
    }
    // Kink-finder score-test outputs: one entry per material step (step 0 =
    // beamline->first hit, including the beampipe). kinkDchisq is the 3-dof
    // score-test Delta-chi2; the Angle/Qop variants test only the direction /
    // only the momentum-step subspace. kinkDqop/Dxdz/Dydz are the best-fit
    // offset at each step (the estimated kink).
    if (doKinkFinder_) {
      tree->Branch("kinkDchisq",      &kinkDchisq);
      tree->Branch("kinkDchisqAngle", &kinkDchisqAngle);
      tree->Branch("kinkDchisqQop",   &kinkDchisqQop);
      tree->Branch("kinkDqop",        &kinkDqop);
      tree->Branch("kinkDxdz",        &kinkDxdz);
      tree->Branch("kinkDydz",        &kinkDydz);
      tree->Branch("kinkGlobalR",     &kinkGlobalR);
      tree->Branch("kinkGlobalZ",     &kinkGlobalZ);
      tree->Branch("kinkMax",      &kinkMax,      basketSize);
      tree->Branch("kinkMaxLayer", &kinkMaxLayer, basketSize);
    }
    // Stage-2 per-row B+ candidate index, branched only when the cfi
    // configured bCandIdxSrc (additive, no-op for legacy J/psi/Upsilon/Z).
    if (!bCandIdxSrcTag_.label().empty()) {
      tree->Branch("bCandIdx", &bCandIdx, basketSize);
    }
    //workaround for older ROOT version inability to store std::array automatically
  // tree->Branch("trackOrigParms", trackOrigParms.data(), "trackOrigParms[5]/F", basketSize);
  // tree->Branch("trackOrigCov", trackOrigCov.data(), "trackOrigCov[25]/F", basketSize);
    tree->Branch("trackParms", trackParms.data(), "trackParms[5]/F", basketSize);
    tree->Branch("trackCov", trackCov.data(), "trackCov[25]/F", basketSize);
    
    tree->Branch("trackHighPurity", &trackHighPurity);

    tree->Branch("refParms_iter0", refParms_iter0.data(), "refParms_iter0[5]/F", basketSize);
    tree->Branch("refCov_iter0", refCov_iter0.data(), "refCov_iter0[25]/F", basketSize);
  // tree->Branch("refParms_iter2", refParms_iter2.data(), "refParms_iter2[5]/F", basketSize);
  // tree->Branch("refCov_iter2", refCov_iter2.data(), "refCov_iter2[25]/F", basketSize); 
    
    tree->Branch("refParms", refParms.data(), "refParms[5]/F", basketSize);
    tree->Branch("refCov", refCov.data(), "refCov[25]/F", basketSize);
    tree->Branch("genParms", genParms.data(), "genParms[5]/F", basketSize);

    tree->Branch("genPt", &genPt, basketSize);
    tree->Branch("genEta", &genEta, basketSize);
    tree->Branch("genPhi", &genPhi, basketSize);
    tree->Branch("genCharge", &genCharge, basketSize);
    
    tree->Branch("genX", &genX, basketSize);
    tree->Branch("genY", &genY, basketSize);
    tree->Branch("genZ", &genZ, basketSize);
    tree->Branch("genPdgId", &genPdgId, basketSize);
    tree->Branch("genDR", &genDR, basketSize);

    // Geant4 decay/interaction truth of the gen-matched particle. simVtx*
    // holds every SimVertex the matched SimTrack produced (decay, nuclear
    // interaction, delta ray, ...) with its G4 process subtype
    // (201 = Decay, 121 = hadronic inelastic, 111 = hadronic elastic,
    // 2 = ionisation/delta ray); simDau* are the flattened daughters,
    // simDauVtx pointing back into simVtx*.
    if (doSimDecayTruth_) {
      tree->Branch("simTrkFound", &simTrkFound, basketSize);
      tree->Branch("simTrkPdgId", &simTrkPdgId, basketSize);
      tree->Branch("simTrkP", &simTrkP, basketSize);
      tree->Branch("simTrkPt", &simTrkPt, basketSize);
      tree->Branch("simTrkEta", &simTrkEta, basketSize);
      tree->Branch("simTrkPhi", &simTrkPhi, basketSize);
      tree->Branch("simTrkVtxX", &simTrkVtxX, basketSize);
      tree->Branch("simTrkVtxY", &simTrkVtxY, basketSize);
      tree->Branch("simTrkVtxZ", &simTrkVtxZ, basketSize);
      tree->Branch("simVtxX", &simVtxX);
      tree->Branch("simVtxY", &simVtxY);
      tree->Branch("simVtxZ", &simVtxZ);
      tree->Branch("simVtxProcType", &simVtxProcType);
      tree->Branch("simVtxNDau", &simVtxNDau);
      tree->Branch("simDauVtx", &simDauVtx);
      tree->Branch("simDauPdgId", &simDauPdgId);
      tree->Branch("simDauPx", &simDauPx);
      tree->Branch("simDauPy", &simDauPy);
      tree->Branch("simDauPz", &simDauPz);
    }


    tree->Branch("normalizedChi2", &normalizedChi2, basketSize);
    
    tree->Branch("nHits", &nHits, basketSize);
    tree->Branch("nValidHits", &nValidHits, basketSize);
    tree->Branch("nValidPixelHits", &nValidPixelHits, basketSize);

    // openspec/improve-cvh-refit-convergence §2: the `nValidHitsFinal` and
    // `nValidPixelHitsFinal` branches previously emitted here were declared,
    // initialised to 0, and never incremented in this single-track producer
    // (the per-hit loop has no `morehitquality` quality gate). The branches
    // therefore wrote literal 0 for every event in the published Run2016H
    // sample, falsely suggesting that 100% of kaon hits had been dropped by
    // the refit (`Kbach_nValidHitsFinal=0` in the joined tree). Removed
    // entirely until a real per-hit rejection mechanism lands and the
    // counters can carry meaningful information.

    if (fillJac_) {
      tree->Branch("nJacRef", &nJacRef, basketSize);
      tree->Branch("jacrefv",jacrefv.data(),"jacrefv[nJacRef]/F", basketSize);
    }
    
    tree->Branch("dEpred", &dEpred);
    tree->Branch("dE", &dE);
    tree->Branch("sigmadE", &sigmadE);
    tree->Branch("E", &E);
    tree->Branch("Epred", &Epred);
    
// tree->Branch("outPt", &outPt);
// tree->Branch("outEta", &outEta);
// tree->Branch("outPhi", &outPhi);
// 
// tree->Branch("outPtStart", &outPtStart);
// tree->Branch("outEtaStart", &outEtaStart);
// tree->Branch("outPhiStart", &outPhiStart);

    tree->Branch("muonPt", &muonPt);
    tree->Branch("muonLoose", &muonLoose);
    tree->Branch("muonMedium", &muonMedium);
    tree->Branch("muonTight", &muonTight);

    tree->Branch("muonIsPF", &muonIsPF);
    tree->Branch("muonIsTracker", &muonIsTracker);
    tree->Branch("muonIsGlobal", &muonIsGlobal);
    tree->Branch("muonIsStandalone", &muonIsStandalone);

    tree->Branch("muonInnerTrackBest", &muonInnerTrackBest);

    tree->Branch("trackExtraAssoc", &trackExtraAssoc);
    
    // Per-hit CLASS variables, booked in BOTH fit modes. The offline CF
    // (cf_track_resolution.py) reads NOMINAL-fit productions
    // (fitFromGenParms=False) and needs to know which class each parmtype-8/9
    // resolution block belongs to, otherwise the hit families can only be
    // summed into one Gaussian -- which is exactly the term the hit-resolution
    // model is meant to replace. The gen-anchored pull study uses the same
    // branches, so there is one definition rather than two.
    if (fillTrackTree_) {
      tree->Branch("dxerr", &dxerr);
      tree->Branch("dyerr", &dyerr);
      tree->Branch("clusterSizeX", &clusterSizeX);
      tree->Branch("clusterSizeY", &clusterSizeY);
      tree->Branch("clusterCharge", &clusterCharge);
      tree->Branch("clusterChargeBin", &clusterChargeBin);
      tree->Branch("clusterOnEdge", &clusterOnEdge);
      tree->Branch("hitDetId", &hitDetId);
      tree->Branch("hitUProj", &hitUProj);
      tree->Branch("hitStripRec", &hitStripRec);
      tree->Branch("hitStripSim", &hitStripSim);
      tree->Branch("hitFirstStrip", &hitFirstStrip);
      tree->Branch("hitPitch", &hitPitch);
      tree->Branch("hitThickness", &hitThickness);
      tree->Branch("localdxdz", &localdxdz);
      tree->Branch("localdydz", &localdydz);
    }

    if (fitFromGenParms_) {
      tree->Branch("hitidxv", &hitidxv);
      tree->Branch("dxrecgen", &dxrecgen);
      tree->Branch("dyrecgen", &dyrecgen);
      tree->Branch("dxsimgen", &dxsimgen);
      tree->Branch("dysimgen", &dysimgen);
      tree->Branch("dxsimgenconv", &dxsimgenconv);
      tree->Branch("dysimgenconv", &dysimgenconv);
      tree->Branch("dxsimgenlocal", &dxsimgenlocal);
      tree->Branch("dysimgenlocal", &dysimgenlocal);
      tree->Branch("dxrecsim", &dxrecsim);
      tree->Branch("dyrecsim", &dyrecsim);
      
      tree->Branch("clusterSize", &clusterSize);
      
      tree->Branch("clusterProbXY", &clusterProbXY);
      tree->Branch("clusterSN", &clusterSN);

      tree->Branch("stripsToEdge", &stripsToEdge);

      tree->Branch("simHitNCand", &simHitNCand);
      
      tree->Branch("dxreccluster", &dxreccluster);
      tree->Branch("dyreccluster", &dyreccluster);

      tree->Branch("simlocalqop", &simlocalqop);
      tree->Branch("simlocaldxdz", &simlocaldxdz);
      tree->Branch("simlocaldydz", &simlocaldydz);
      tree->Branch("simlocalx", &simlocalx);
      tree->Branch("simlocaly", &simlocaly);
      
      tree->Branch("simlocalqopprop", &simlocalqopprop);
      tree->Branch("simlocaldxdzprop", &simlocaldxdzprop);
      tree->Branch("simlocaldydzprop", &simlocaldydzprop);
      tree->Branch("simlocalxprop", &simlocalxprop);
      tree->Branch("simlocalyprop", &simlocalyprop);

      tree->Branch("landauDelta", &landauDelta);
      tree->Branch("landauW", &landauW);
      
      tree->Branch("localqop", &localqop);
      tree->Branch("localx", &localx);
      tree->Branch("localy", &localy);
      
      tree->Branch("localqoperr", &localqoperr);
      tree->Branch("localdxdzerr", &localdxdzerr);
      tree->Branch("localdydzerr", &localdydzerr);
      tree->Branch("localxerr", &localxerr);
      tree->Branch("localyerr", &localyerr);

      tree->Branch("hitlocalx", &hitlocalx);
      tree->Branch("hitlocaly", &hitlocaly);
      
      tree->Branch("localqop_iter", &localqop_iter);
      tree->Branch("localdxdz_iter", &localdxdz_iter);
      tree->Branch("localdydz_iter", &localdydz_iter);
      tree->Branch("localx_iter", &localx_iter);
      tree->Branch("localy_iter", &localy_iter);

      tree->Branch("dxrecgen_iter", &dxrecgen_iter);
      tree->Branch("dyrecgen_iter", &dyrecgen_iter);
      
      tree->Branch("localqoperralt", &localqoperralt);
      
      tree->Branch("localphi", &localphi);
      tree->Branch("hitphi", &hitphi);

      
      tree->Branch("simtestz", &simtestz);
      tree->Branch("simtestvz", &simtestvz);
      tree->Branch("simtestrho", &simtestrho);
      tree->Branch("simtestzlocalref", &simtestzlocalref);
      tree->Branch("simtestdx", &simtestdx);
      tree->Branch("simtestdxrec", &simtestdxrec);
      tree->Branch("simtestdy", &simtestdy);
      tree->Branch("simtestdyrec", &simtestdyrec);
      tree->Branch("simtestdxprop", &simtestdxprop);
      tree->Branch("simtestdyprop", &simtestdyprop);
      tree->Branch("simtestdetid", &simtestdetid);
      
      tree->Branch("rx", &rx);
      tree->Branch("ry", &ry);
      
      tree->Branch("deigx", &deigx);
      tree->Branch("deigy", &deigy);
    }
    
    
    nJacRef = 0.;
  }
}


// ------------ method called for each event ------------
void ResidualGlobalCorrectionMakerG4e::produce(edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  siStripClusterInfo_.initEvent(iSetup);

  // Sync the material-group k values from corparms_ into the model so the
  // propagator's per-step provider applies the current calibration.
  if (globalMaterialModel_) {
    for (unsigned int g = 0; g < matGroupGlobalIdx_.size(); ++g) {
      matModel_->setKValue(g, corparms_[matGroupGlobalIdx_[g]]);
    }
  }

  const bool dogen = fitFromGenParms_;
  const bool dolocalupdate = fitFromSimParms_;
  // const bool dolocalupdate = true;
  
  const bool dores = doRes_;

  using namespace edm;

  Handle<reco::TrackCollection> trackOrigH;
  iEvent.getByToken(inputTrackOrig_, trackOrigH);

  // Optional Stage-2 B+ candidate index, parallel to trackOrigH. Empty
  // handle when bCandIdxSrc is not configured -- per-row bCandIdx stays at
  // its -1 sentinel and the branch was not added in beginStream().
  Handle<std::vector<int>> bCandIdxH;
  if (!bCandIdxSrcTag_.label().empty()) {
    iEvent.getByToken(bCandIdxToken_, bCandIdxH);
  }

  auto globalGeometry = iSetup.getHandle(globalGeometryEventToken_);
  auto trackerTopology = iSetup.getHandle(trackerTopologyEventToken_);
  auto ttrh = iSetup.getHandle(ttrhToken_);

  // MT: bootstrap this TBB worker thread's G4 environment from the master
  // (world + per-thread navigator + per-thread magnetic field). Idempotent
  // per thread; runs on every produce() but only does work on the first.
  worker_->ensureInitialized(iSetup.getData(cvhMasterToken_).cvhMaster());
  // Bind this thread's G4 RNG to the stream's CLHEP engine for this event.
  setG4RandomEngineForStream(iEvent.streamID());
  // Lazy-init the per-stream propagator clone on first call. Clone now is
  // safe -- the deep-copy ctor no longer eagerly allocates fluct (that
  // happens lazily in propagateGenericWithJacobianAltD's first-call init,
  // AFTER ensureGeant4eIsInitilizedForCVH has registered the G4Error
  // physics on the newly-bootstrapped thread).
  if (!streamPropagator_) {
    auto thePropagator = iSetup.getHandle(g4ePropToken_);
    const Geant4ePropagator *templateProp =
        dynamic_cast<const Geant4ePropagator*>(thePropagator.product());
    if (!templateProp) {
      throw cms::Exception("Configuration")
          << "ESProducer for label 'Geant4ePropagator' did not deliver a Geant4ePropagator";
    }
    streamPropagator_.reset(templateProp->clone());
    // Urban step log feeds the ioniurbanv physics-CF export; only pay the
    // bookkeeping when the resolution machinery is on and grads are kept.
    streamPropagator_->setIoniStepLogging(doRes_ && (fillGrads_ || fillGradsFactored_));
  }
  const Geant4ePropagator *g4prop = streamPropagator_.get();
  const MagneticField* field = g4prop->magneticField();
  
  // Track mass for energy / Jacobian calculations -- read from cfi
  // (default = muon mass for J/psi/Upsilon back-compat). Renamed from
  // the previous local `constexpr double trackmass = 0.1056583745;`.
  const double trackmass = trackMass_;

  Handle<reco::BeamSpot> bsH;
  iEvent.getByToken(inputBs_, bsH);

  
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

  edm::Handle<std::vector<SimVertex>> simVertices;
  if (doSimDecayTruth_) {
    if (!doSim_) {
      iEvent.getByToken(genParticlesBarcodeToken_, genPartBarcodes);
      iEvent.getByToken(inputSimTracks_, simTracks);
    }
    iEvent.getByToken(inputSimVertices_, simVertices);
  }

  Handle<edm::View<reco::Muon> > muons;
  if (doMuons_ || doMuonTrackAssoc_) {
    iEvent.getByToken(inputMuons_, muons);
  }
  
  Handle<edm::Association<std::vector<pat::Muon>>> muonAssoc;
  if (doMuonAssoc_) {
    iEvent.getByToken(inputMuonAssoc_, muonAssoc);
  }
  
  TkClonerImpl const& cloner = static_cast<TkTransientTrackingRecHitBuilder const *>(ttrh.product())->cloner();

  // Same StripCPEfromTrackAngle instance the cloner uses (single ES
  // product for the label), queried only for its AlgoParam so the
  // hit-resolution export can record the CPE's own uProj.
  // NOT gated on fitFromGenParms_: hitUProj is exactly the class variable the
  // offline CF needs on the NOMINAL-fit productions, and gating it there left
  // every strip hit at -99 (found by the conditioning study, which then put
  // the whole strip tracker in one uProj bin).
  const StripCPE *stripCPEForExport = nullptr;
  if (fillTrackTree_) {
    stripCPEForExport = dynamic_cast<const StripCPE*>(iSetup.getHandle(stripCPEToken_).product());
  }

  run = iEvent.run();
  lumi = iEvent.luminosityBlock();
  event = iEvent.id().event();

  genweight = 1.;
  if (doGen_) {
    genweight = genEventInfo->weight();

    Pileup_nPU = pileupSummary->front().getPU_NumInteractions();
    Pileup_nTrueInt = pileupSummary->front().getTrueNumInteractions();
  }

  std::vector<float> corPtV;
  std::vector<float> corEtaV;
  std::vector<float> corPhiV;
  std::vector<int> corChargeV;
  std::vector<float> corDxyV;
  std::vector<float> corDzV;
  std::vector<float> edmvalV;
  std::vector<int> nValidHitsV;
  std::vector<int> nValidPixelHitsV;

  std::vector<std::vector<int>> globalidxsV;
  std::vector<std::vector<float>> jacRefV;
  std::vector<std::vector<float>> momCovV;

  std::array<double, 3> refParmsMomD;

  // Refit-track output, positionally aligned with the input tracks. Seeded
  // with copies of the inputs so a failed fit leaves a valid (uncorrected)
  // entry that refitOk marks as such.
  std::vector<reco::Track> refitTracksV;
  std::vector<int> refitOkV;
  if (emitRefitTracks_) {
    refitTracksV.assign(trackOrigH->begin(), trackOrigH->end());
    refitOkV.assign(trackOrigH->size(), 0);
  }

  if (doMuonAssoc_) {
    corPtV.assign(muonAssoc->ref()->size(), -99.);
    corEtaV.assign(muonAssoc->ref()->size(), -99.);
    corPhiV.assign(muonAssoc->ref()->size(), -99.);
    corChargeV.assign(muonAssoc->ref()->size(), -99);
    corDxyV.assign(muonAssoc->ref()->size(), -99.);
    corDzV.assign(muonAssoc->ref()->size(), -99.);
    edmvalV.assign(muonAssoc->ref()->size(), -99.);
    nValidHitsV.assign(muonAssoc->ref()->size(), -99);
    nValidPixelHitsV.assign(muonAssoc->ref()->size(), -99);
    globalidxsV.assign(muonAssoc->ref()->size(), std::vector<int>());
    jacRefV.assign(muonAssoc->ref()->size(), std::vector<float>());
    momCovV.assign(muonAssoc->ref()->size(), std::vector<float>());
  }

// for (const reco::Track &track : *trackOrigH) {
  for (unsigned int itrack = 0; itrack < trackOrigH->size(); ++itrack) {
    const reco::Track &track = (*trackOrigH)[itrack];
    const reco::TrackRef trackref(trackOrigH, itrack);

    // Stage-2: positionally read the bCandIdx for this row, defaulting to
    // -1 if the optional input is absent or shorter than expected.
    bCandIdx = (bCandIdxH.isValid() && itrack < bCandIdxH->size())
        ? (*bCandIdxH)[itrack] : -1;

    // Note: the Geant4 particle name (charge-dependent) is built per fit
    // attempt inside the runHypothesis lambda below, so the two-hypothesis
    // fit can flip the charge consistently with refFts[6].

    const edm::Ref<std::vector<pat::Muon>> muonref = doMuonAssoc_ ? (*muonAssoc)[trackref] : edm::Ref<std::vector<pat::Muon>>();

    // forceCosmic_ makes split-leg tracks (algo 'cosmic', not matched below) use
    // the outside-in cosmic handling instead of the beamline-PCA parameterisation.
    const bool iscosmic = forceCosmic_ || track.algo() == reco::TrackBase::ctf || track.algo() == reco::TrackBase::cosmics;
    
    const bool dopca = !dogen && !iscosmic;
    
    if (track.isLooper()) {
      continue;
    }
    
    trackPt = track.pt();
    trackEta = track.eta();
    trackPhi = track.phi();
    trackCharge = track.charge();
    trackPtErr = track.ptError();
    trackQopErr = track.qoverpError();

// if (abs(trackEta) > 0.1) {
// continue;
// }

    trackHighPurity = track.quality(reco::TrackBase::highPurity);
    
    normalizedChi2 = track.normalizedChi2();
    
    auto const& tkparms = track.parameters();
    auto const& tkcov = track.covariance();
    trackParms.fill(0.);
    trackCov.fill(0.);
    //use eigen to fill raw memory
    Map<Vector5f>(trackParms.data()) = Map<const Vector5d>(tkparms.Array()).cast<float>();
    Map<Matrix<float, 5, 5, RowMajor> >(trackCov.data()).triangularView<Upper>() = Map<const Matrix<double, 5, 5, RowMajor> >(tkcov.Array()).cast<float>().triangularView<Upper>();
    
// std::cout << "track charge: " << track.charge() << " trackorig charge " << trackOrig.charge() << "inner state charge " << tms.back().updatedState().charge() << std::endl;
    
    const reco::Candidate* genpart = nullptr;
    
    genPt = -99.;
    genEta = -99.;
    genPhi = -99.;
    genCharge = -99;
    genX = -99.;
    genY = -99.;
    genZ = -99.;
    genParms.fill(0.);
    genl3d = -99.;
    genPdgId = 0;
    genDR = -99.f;
    simPabsFirst = -99.;
    simPabsLast = -99.;

    if (doSimDecayTruth_) {
      simTrkFound = 0;
      simTrkPdgId = 0;
      simTrkP = -99.f;
      simTrkPt = -99.f;
      simTrkEta = -99.f;
      simTrkPhi = -99.f;
      simTrkVtxX = -99.f;
      simTrkVtxY = -99.f;
      simTrkVtxZ = -99.f;
      simVtxX.clear();
      simVtxY.clear();
      simVtxZ.clear();
      simVtxProcType.clear();
      simVtxNDau.clear();
      simDauVtx.clear();
      simDauPdgId.clear();
      simDauPx.clear();
      simDauPy.clear();
      simDauPz.clear();
    }
    
    int genBarcode = -99;
    
    
    if (doGen_) {
      
      float drmin = genMatchDR_;

      for (auto g = genPartCollection->begin(); g != genPartCollection->end(); ++g)
      {
        if (g->status() != 1) {
          continue;
        }
        // Species allowed into the dR competition (see genMatchPdgIds_).
        const int abspdg = std::abs(g->pdgId());
        if (genMatchPdgIds_.empty()) {
          if (abspdg != genMatchPdgId_) {
            continue;
          }
        }
        else if (std::find(genMatchPdgIds_.begin(), genMatchPdgIds_.end(), abspdg)
                 == genMatchPdgIds_.end()) {
          continue;
        }
        // Same-charge requirement as in the two-track maker. Pure-dR
        // matching is unsafe in busy MC (B->J/psi+X): soft junk/hadron
        // tracks collinear with a muon (dR<0.1 is common for collimated
        // J/psi daughters) get anchored to the muon gen state, and with
        // fitFromGenParms that produces chi2 ~ 1e6 outliers that dominate
        // any variance-sensitive fit.
        if (g->charge() != track.charge()) {
          continue;
        }
        // Loose momentum-compatibility window (gen-matched muon tracks
        // agree to ~1%; hadron/junk mismatches are off by 10x). Configurable
        // because decay-in-flight studies need it loosened: the reconstructed
        // pT of a decayed kaon/pion follows the daughter, far below gen.
        if (std::abs(g->pt() - track.pt()) > genMatchPtWindow_ * g->pt()) {
          continue;
        }

        float dR = deltaR(*g, track);
        
        if (dR < drmin)
        {
          drmin = dR;
          
          genpart = &(*g);
          
          if (doSim_ || doSimDecayTruth_) {
            genBarcode = (*genPartBarcodes)[g - genPartCollection->begin()];
          }
          
          genPt = g->pt();
          genEta = g->eta();
          genPhi = g->phi();
          genCharge = g->charge();
          genPdgId = g->pdgId();
          genDR = dR;
          
          genX = g->vertex().x();
          genY = g->vertex().y();
          genZ = g->vertex().z();

          // genParticles:xyz0 (gen PV) is not kept in every ALCARECO
          // (the B->J/psi+X MC keeps only the recoGenParticles branch);
          // fall back to the -99 sentinel rather than throwing. (This
          // guard was originally a pixel-session working-tree fix that
          // was lost in the 2026-07-25 session disentangling.)
          genl3d = genXyz0.isValid()
              ? std::sqrt((g->vertex() - *genXyz0).mag2()) : -99.;

          auto const& vtx = g->vertex();

          // GEN REFERENCE PARAMETERS IN THE FIT'S OWN CONVENTION.
          //
          // They used to be hand-coded in the reco::TrackBase perigee
          // convention, which agrees with `refParms` on (q/p, lambda, phi,
          // d0) and NOT on the fifth:
          //   * `refParms[4]` is `cart2pca(...)[4]`, the ABSOLUTE z of the
          //     point of closest approach to the beamline, whereas the old
          //     line computed `dsz` -- and computed it wrongly, because
          //     `BeamSpot::position(z)` returns `Point(x(z), y(z), z)`, so
          //     the `vtx.z() - myBeamSpot.z()` term was IDENTICALLY ZERO and
          //     `genParms[4]` was the transverse remainder alone (~1e-4 cm
          //     against a z0 of centimetres).  Measured on 2000 mu-gun
          //     tracks: `Var((refParms[4]-genParms[4])/sigma) = 4.0e6` with a
          //     median pull of +375, against ~1.0 for the other four
          //     (NOTES 2026-09-09).
          //   * `genParms[3]` was `dxy`, which IS `cart2pca`'s `d0` for a
          //     beamline along z-hat -- no change in value, but now it is the
          //     same formula rather than two that happen to agree.
          // Calling `cart2pca` on the gen state makes gen and fitted
          // DEFINITIONALLY identical, which is what a residual needs.
          {
            Matrix<double, 7, 1> genstate;
            genstate << vtx.x(), vtx.y(), vtx.z(), g->px(), g->py(), g->pz(),
                double(g->charge());
            const Matrix<double, 5, 1> genpca = cart2pca(genstate, *bsH);
            for (unsigned int i = 0; i < 5; ++i) {
              genParms[i] = float(genpca[i]);
            }
          }
        }
        else {
          continue;
        }
      }

      // Cross-species arbitration: the winner of the dR competition must be
      // the species this pass is fitting, otherwise the track is treated as
      // unmatched. This removes the mismatch background without a pT cut --
      // a pT window would preferentially discard the genuine hard decays,
      // which are the signal.
      if (genpart != nullptr && !genMatchPdgIds_.empty()
          && std::abs(genPdgId) != genMatchPdgId_) {
        genpart = nullptr;
        genBarcode = -99;
        genPt = -99.;
        genEta = -99.;
        genPhi = -99.;
        genCharge = -99;
        genX = -99.;
        genY = -99.;
        genZ = -99.;
        genParms.fill(0.);
        genl3d = -99.;
      }
    }
    
// std::cout << "genPt = " << genPt << " genEta = " << genEta << " genPhi = " << genPhi << " genCharge = " << genCharge << " genX = " << genX << " genY = " << genY << " genZ = " << genZ << std::endl;
    
    if (requireGen_ && genpart == nullptr) {
      continue;
    }
    
    int simtrackid = -99;
    if (genpart != nullptr && (doSim_ || doSimDecayTruth_)) {
      for (auto const& simTrack : *simTracks) {
        if (simTrack.genpartIndex() == genBarcode) {
          simtrackid = simTrack.trackId();
          if (doSimDecayTruth_) {
            simTrkFound = 1;
            simTrkPdgId = simTrack.type();
            simTrkP = simTrack.momentum().P();
            simTrkPt = simTrack.momentum().pt();
            simTrkEta = simTrack.momentum().eta();
            simTrkPhi = simTrack.momentum().phi();
            if (!simTrack.noVertex()) {
              auto const& pv = (*simVertices)[simTrack.vertIndex()].position();
              simTrkVtxX = pv.x();
              simTrkVtxY = pv.y();
              simTrkVtxZ = pv.z();
            }
          }
          break;
        }
      }
    }

    // Every Geant4 vertex produced by the matched SimTrack: the decay
    // (processType 201) that this study targets, but also nuclear
    // interactions (111/121) and delta rays (2), which are the physics
    // backgrounds to the kink tag. Daughters are flattened with a back
    // pointer so the offline classification can use the full final state.
    if (doSimDecayTruth_ && simtrackid >= 0) {
      for (unsigned int iv = 0; iv < simVertices->size(); ++iv) {
        auto const& simVertex = (*simVertices)[iv];
        if (simVertex.noParent() || int(simVertex.parentIndex()) != simtrackid) {
          continue;
        }
        const int islot = int(simVtxX.size());
        auto const& vpos = simVertex.position();
        simVtxX.push_back(vpos.x());
        simVtxY.push_back(vpos.y());
        simVtxZ.push_back(vpos.z());
        simVtxProcType.push_back(int(simVertex.processType()));
        int ndau = 0;
        for (auto const& simTrack : *simTracks) {
          if (simTrack.noVertex() || simTrack.vertIndex() != int(iv)) {
            continue;
          }
          ++ndau;
          simDauVtx.push_back(islot);
          simDauPdgId.push_back(simTrack.type());
          simDauPx.push_back(simTrack.momentum().x());
          simDauPy.push_back(simTrack.momentum().y());
          simDauPz.push_back(simTrack.momentum().z());
        }
        simVtxNDau.push_back(ndau);
      }
    }


    
    muonPt = -99.;
    muonLoose = false;
    muonMedium = false;
    muonTight = false;
    muonIsTracker = false;
    muonIsGlobal = false;
    // muonIsPF was missing from this reset, so it carried over from the
    // previous track in the event: on the B->J/psi+X ALCARECO that made it
    // true for ~47% of hadron tracks that were never matched to a muon.
    muonIsPF = false;
    muonIsStandalone = false;
    muonInnerTrackBest = false;
    trackExtraAssoc = false;

    const reco::Muon *matchedmuon = nullptr;

    if (doMuons_) {
      for (auto const &muon : *muons) {
        if (muon.bestTrack()->algo() == track.algo()) {
          if ( (muon.bestTrack()->momentum() - track.momentum()).mag2() < 1e-3 ) {
            matchedmuon = &muon;
          }
        }
        else if (muon.innerTrack().isNonnull() && muon.innerTrack()->algo() == track.algo()) {
          if ( (muon.innerTrack()->momentum() - track.momentum()).mag2() < 1e-3 ) {
            matchedmuon = &muon;
          }
        }
      }
    }

    // Ref-safe muon match for ALCARECO inputs. The muon's bestTrack()/
    // innerTrack() refs point into generalTracks, which the ALCARECO drops,
    // so dereferencing them throws; and the shipped track->muon association
    // is keyed on the unselected track collection, not the preselected copy
    // this module runs on. Match on the muon's OWN four-momentum instead,
    // and read only values stored directly on the muon.
    if (doMuonTrackAssoc_ && muons.isValid()) {
      const reco::Muon *best = nullptr;
      double bestdr = 0.01;
      for (auto const &muon : *muons) {
        if (muon.charge() != track.charge()) {
          continue;
        }
        if (std::abs(muon.pt() - track.pt()) > 0.05*track.pt()) {
          continue;
        }
        const double dr = deltaR(muon, track);
        if (dr < bestdr) {
          bestdr = dr;
          best = &muon;
        }
      }
      if (best != nullptr) {
        muonPt = best->pt();
        muonLoose = best->passed(reco::Muon::CutBasedIdLoose);
        muonMedium = best->passed(reco::Muon::CutBasedIdMedium);
        muonTight = best->passed(reco::Muon::CutBasedIdTight);
        muonIsPF = best->isPFMuon();
        muonIsTracker = best->isTrackerMuon();
        muonIsGlobal = best->isGlobalMuon();
        muonIsStandalone = best->isStandAloneMuon();
        muonInnerTrackBest = best->muonBestTrackType() == reco::Muon::InnerTrack;
      }
    }

    if (matchedmuon != nullptr) {
      muonPt = matchedmuon->pt();
      muonLoose = matchedmuon->passed(reco::Muon::CutBasedIdLoose);
      muonMedium = matchedmuon->passed(reco::Muon::CutBasedIdMedium);
      muonTight = matchedmuon->passed(reco::Muon::CutBasedIdTight);
      muonIsPF = matchedmuon->isPFMuon();
      muonIsTracker = matchedmuon->isTrackerMuon();
      muonIsGlobal = matchedmuon->isGlobalMuon();
      muonIsStandalone = matchedmuon->isStandAloneMuon();
      muonInnerTrackBest = matchedmuon->muonBestTrackType() == reco::Muon::InnerTrack;
    }
    
    //prepare hits
    TransientTrackingRecHit::RecHitContainer hits;
    hits.reserve(track.recHitsSize());

// std::cout << "track: algo = " << track.algo() << " originalAlgo = " << track.originalAlgo() << std::endl;
//
//
// if (track.seedDirection() == oppositeToMomentum) {
// // std::cout << "track with oppositeToMomentum: algo = " << track.algo() << " originalAlgo = " << track.originalAlgo() << std::endl;
// std::cout << "track with oppositeToMomentum:" << std::endl;
// }


    
    for (auto it = track.recHitsBegin(); it != track.recHitsEnd(); ++it) {
// if (track.seedDirection() == oppositeToMomentum) {
// std::cout << "det = " << (*it)->geographicalId().det() << std::endl;
// }


      if ((*it)->geographicalId().det() != DetId::Tracker) {
        continue;
      }

      // Leg-structure-free mode: hitless surfaces (dead-module placeholders)
      // carry no measurement and, with the global models, no parameter
      // attribution either -- drop them so propagation goes hit to hit.
      // The material of the skipped modules is still crossed by the G4
      // steps of the longer leg.
      if (skipHitlessSurfaces_ && !(*it)->isValid()) {
        continue;
      }

      // hits on garbage-shifted modules: either dropped here (drop policy) or
      // re-inserted at the repaired-surface path position after collection
      // (reorder policy) -- the stored hit order comes from the garbage
      // constants and can imply backward propagation steps that abort the fit
      if (!garbageShiftReorderHits_ && garbageShiftModules_.count((*it)->geographicalId().rawId())) {
        continue;
      }

      const GeomDet* detectorG = globalGeometry->idToDet((*it)->geographicalId());
      const GluedGeomDet* detglued = dynamic_cast<const GluedGeomDet*>(detectorG);
      
// if (track.seedDirection() == oppositeToMomentum) {
// std::cout << "position mag = " << detectorG->surface().position().mag() << std::endl;
// }

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

// std::cout << "splitting matched hit, inner = " << detinner->geographicalId().rawId() <<" outer = " << detouter->geographicalId().rawId() << std::endl;
        if (hits.size() > 0) {
// std::cout << "previous = " << hits.back()->geographicalId().rawId() << std::endl;
          const bool duplicate = detinner->geographicalId() == hits.back()->geographicalId() || detouter->geographicalId() == hits.back()->geographicalId();

          if (duplicate) {
            std::cout << "WARNING: Duplicate hits from glued module splitting!\n";
          }
        }
        
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
// bool hitquality = true;
        bool hitquality = false;
        if (applyHitQuality_ && (*it)->isValid()) {
          const TrackerSingleRecHit* tkhit = dynamic_cast<const TrackerSingleRecHit*>(*it);
          assert(tkhit != nullptr);
          
          if (ispixel) {
            const SiPixelRecHit *pixhit = dynamic_cast<const SiPixelRecHit*>(tkhit);
            const SiPixelCluster& cluster = *tkhit->cluster_pixel();
            assert(pixhit != nullptr);
            
            ++pixHitsSeen_;
            const bool onEdge = pixhit->isOnEdge();
            if (onEdge) ++pixHitsEdge_;
            if (cluster.sizeX() <= 1) ++pixHitsSizeX1_;
            // Direction-resolved pathology classes (see the two-track maker).
            {
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
            }
            // Boundary veto configurable via keepPixelEdgeHits; sizeX
            // threshold configurable via pixelMinSizeX (default 2 = legacy).
            hitquality = (keepPixelEdgeHits_ || !onEdge) && cluster.sizeX() >= pixelMinSizeX_
                          && cluster.sizeY() >= pixelMinSizeY_;
            if (!hitquality) ++pixHitsDemoted_;
// hitquality = false;
          }
          else {
            assert(tkhit->cluster_strip().isNonnull());
            const SiStripCluster& cluster = *tkhit->cluster_strip();
            const StripTopology* striptopology = dynamic_cast<const StripTopology*>(&(detectorG->topology()));
            assert(striptopology);

            const uint16_t firstStrip = cluster.firstStrip();
            const uint16_t lastStrip = cluster.firstStrip() + cluster.amplitudes().size() - 1;
            const bool isOnEdge = firstStrip == 0 || lastStrip == (striptopology->nstrips() - 1);

            const bool isstereo = trackerTopology->isStereo((*it)->geographicalId());
            
            hitquality = true;
// hitquality = false;
            // hitquality = !isstereo;
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

    const unsigned int nhits = hits.size();

    if (nhits == 0) {
      continue;
    }

    nHits = nhits;

    // Cosmic forward-seeding option: if the muon enters at the OUTER end (higher
    // y than the inner end), reverse the hit order so the fit propagates from
    // the entry end along the muon (downward) = fully forward. The seed state
    // is chosen accordingly below. Default off (inner seed, upper half backward).
    bool cosmicFromOuter = false;
    if (cosmicSeedFromEntry_ && iscosmic && track.extra().isNonnull()) {
      // seed from whichever end is higher (the down-going muon's entry), so the
      // propagation follows the muon downward = forward. For split legs the
      // inner end is already the higher one (no-op); for a full cosmic the
      // outer extremity may be on top, in which case seed there + reverse hits.
      cosmicFromOuter =
          track.extra()->outerPosition().y() > track.extra()->innerPosition().y();
      if (cosmicFromOuter) {
        std::reverse(hits.begin(), hits.end());
      }
    }

    unsigned int nvalid = 0;
    unsigned int nvalidpixel = 0;
    unsigned int nvalidalign2d = 0;
    
    
    // count valid hits since this is needed to size the arrays
    for (auto const& hit : hits) {
      
      assert(hit->dimension()<=2);
      if (hit->isValid()) {
        nvalid += 1;
        
        const uint32_t gluedid = trackerTopology->glued(hit->geographicalId());
        const bool isglued = gluedid != 0;
        // const DetId parmdetid = isglued ? DetId(gluedid) : hit->geographicalId();

        // const DetId aligndetid = alignGlued_ ? parmdetid : hit->geographicalId();

        
        const bool align2d = detidparms.count(std::make_pair(1, hit->geographicalId()));
        // const bool align2d = detidparms.count(std::make_pair(1, aligndetid));
// 
        if (align2d) {
          nvalidalign2d += 1;
        }
        if (GeomDetEnumerators::isTrackerPixel(hit->det()->subDetector())) {
          nvalidpixel += 1;
        }
      }
    }

    if (nvalid == 0) {
      continue;
    }
    
    nValidHits = nvalid;
    nValidPixelHits = nvalidpixel;
    
    const unsigned int nparsAlignment = 5*nvalid + nvalidalign2d;
    const unsigned int nFieldModes = fieldCorrection_->nModes();
    const unsigned int nparsBfield = nhits * nFieldModes;
    // Global material model: one slot per group per hit (uncrossed groups
    // contribute zero columns; the idxmap collapses the shared global
    // indices exactly like the field-mode block). Legacy: one per-module
    // eloss slot per hit.
    const unsigned int nMatGroups = globalMaterialModel_ ? matModel_->nGroups() : 0;
    const unsigned int nparsEloss = globalMaterialModel_ ? nhits * nMatGroups : nhits;
// const unsigned int nparsRes = nhits + nvalid + nvalidpixel;
    const unsigned int nparsRes = dores ? 2*nhits + nvalid + nvalidpixel : 0;
    const unsigned int npars = nparsAlignment + nparsBfield + nparsEloss + nparsRes;
    
    const unsigned int nstateparms = 5*(nhits+1);
    const unsigned int nparmsfull = nstateparms + npars;
    
    const int nbscons = bsConstraint_ ? 3 : 0;
    const int ncons = 5*nhits + nvalid + nvalidpixel + nbscons;
    
    using VectorXb = Matrix<bool, Dynamic, 1>;
    VectorXb freestatemask = VectorXb::Ones(nstateparms);
    
    if (dogen) {
      freestatemask.head<5>() = Matrix<bool, 5, 1>::Zero();
    }
    if (fitFromSimParms_) {
      freestatemask = VectorXb::Zero(nstateparms);
    }

    // if (dogen) {
    // for (unsigned int istate = 0; istate < nstateparms; ++istate) {
    // if (istate % 5 == 0) {
    // freestatemask[istate] = false;
    // }
    // }
    // }
    
    std::vector<Eigen::Index> freestateidxs;
    freestateidxs.reserve(nstateparms);
    
    for (Eigen::Index istate = 0; istate < nstateparms; ++istate) {
      if (freestatemask[istate]) {
        freestateidxs.push_back(istate);
      }
    }
    
    const unsigned int nstatefree = freestateidxs.size();
    
    ndof = ncons - nstatefree;
    
    if (bsConstraint_) {
      ndof += 2;
    }
    
    MatrixXd validdxeigjac;
    VectorXd validdxeig;
    
    Matrix<float, Dynamic, 2> rxfull(nvalid, 2); 
    Matrix<float, Dynamic, 2> ryfull(nvalid, 2);
    
    globalidxv.clear();
    globalidxv.resize(npars, 0);
    
    VectorXd dxfull;
    MatrixXd dxdparms;
    VectorXd grad;
    MatrixXd hess;
// LDLT<MatrixXd> Cinvd;
    MatrixXd covfull = MatrixXd::Zero(nstateparms, nstateparms);
    
    VectorXd rfull;
    MatrixXd Ffull;
    MatrixXd Jfull;
    MatrixXd Vinvfull;
    VectorXd dxfree;
    
    MatrixXd Vinvfullalt;
    
    SparseMatrix<double> Fsparse;
    SparseMatrix<double> Vinvsparse;
    SparseMatrix<double> VinvF;
    
    SimplicialLDLT<SparseMatrix<double>> Cinvd;
    
    std::vector<SparseMatrix<double>> dVs;
    dVs.reserve(nhits + nvalid);

    // per dVs entry: [row offset, block size] within the constraint vector
    // and the entry's global parameter index -- inputs to the exact
    // block-eigenvalue export (reseigidx/reseigv)
    std::vector<std::array<unsigned int, 2>> resblockrng;
    resblockrng.reserve(nhits + nvalid);
    std::vector<unsigned int> resglobidx;
    resglobidx.reserve(nhits + nvalid);
    
    std::vector<unsigned int> residxs;
    residxs.reserve(nhits + nvalid);
    
    if (dogen && genpart==nullptr) {
      std::cout << "no gen part, skipping track\n";
      continue;
    }
    

    if (debugprintout_) {
      std::cout << "initial reference point parameters:" << std::endl;
      std::cout << track.parameters() << std::endl;
    }

    Matrix<double, 7, 1> refFts;
    
    if (dogen) {
      if (genpart==nullptr) {
        continue;
      }
      //init from gen state
      auto const& refpoint = genpart->vertex();
      auto const& trackmom = genpart->momentum();
      
      refFts[0] = refpoint.x();
      refFts[1] = refpoint.y();
      refFts[2] = refpoint.z();
      
      refFts[3] = trackmom.x();
      refFts[4] = trackmom.y();
      refFts[5] = trackmom.z();
      
      refFts[6] = genpart->charge();
    }
    else {
      //init from track state
      auto const& refpoint = track.referencePoint();
      auto const& trackmom = track.momentum();
      
      refFts[0] = refpoint.x();
      refFts[1] = refpoint.y();
      refFts[2] = refpoint.z();
      
      refFts[3] = trackmom.x();
      refFts[4] = trackmom.y();
      refFts[5] = trackmom.z();

      // seedChargeSign_ (default +1) flips the charge hypothesis (two-hypothesis test)
      refFts[6] = seedChargeSign_ * track.charge();

      // special case for cosmics, start from "inner" state with straight line extrapolation back 1cm to preserve the propagation logic
      // (the reference point is instead the PCA to the beamline and is therefore in the middle of the trajectory and not compatible
      // with the fitter/propagation logic
      // TODO change this so that cosmics directly use a local parameterization on the first hit?
      if (iscosmic) {
        // seed from the outer (entry) end when cosmicFromOuter, else the inner end
        auto const& startpoint = cosmicFromOuter ? track.extra()->outerPosition()
                                                 : track.extra()->innerPosition();
        auto const& startmom = cosmicFromOuter ? track.extra()->outerMomentum()
                                               : track.extra()->innerMomentum();

        const Eigen::Vector3d startpointv(startpoint.x(), startpoint.y(), startpoint.z());
        const Eigen::Vector3d startmomv(startmom.x(), startmom.y(), startmom.z());

        refFts.head<3>() = startpointv - 1.0*startmomv.normalized();
        refFts.segment<3>(3) = startmomv;

      }
    }

    // [openspec §9.4.a] One-shot kaon-q/p diagnostic. Prints the input
    // track perigee, the cartesian `refFts` state we just built, and the
    // implied input q/p. Repeated at end-of-iteration (around line ~2500
    // post-`qbpupd`) to compare against the refit's q/p. First N kaon
    // tracks only to keep logs short.
    static std::atomic<int> kaonQpDumpCount{0};
    constexpr int kKaonQpDumpMax = 5;
    // Gate on the B->J/psiK bachelor path (bCandIdxSrc configured) rather
    // than on the particle name: this keeps the diagnostic for the §9.4.c
    // mass-hypothesis A/B (kaonAsMuon=True => trackParticleName="mu" on the
    // kaon collection) WITHOUT firing in plain muon single-track jobs
    // (runCvhSingleTrack.py), which use this same plugin. fetch_add first
    // so the cap cannot be exceeded by concurrent streams.
    const bool dumpKaonQp = !bCandIdxSrcTag_.label().empty()
                          && (kaonQpDumpCount.fetch_add(1) < kKaonQpDumpMax);
    if (dumpKaonQp) {
      const double inP = std::sqrt(track.momentum().mag2());
      const double refP = refFts.segment<3>(3).norm();
      const double qOverPin = track.charge() / std::max(inP, 1e-12);
      const double qOverPref0 = refFts[6] / std::max(refP, 1e-12);
      std::cout << "===== [§9.4.a kaon q/p ENTRY] track #" << kaonQpDumpCount.load() << " =====\n"
                << "  input track : charge=" << track.charge()
                << "  pt=" << track.pt() << "  eta=" << track.eta()
                << "  phi=" << track.phi() << "  |p|=" << inP
                << "  q/p=" << qOverPin << "\n"
                << "  refFts(in)  : charge[6]=" << refFts[6]
                << "  mom=("  << refFts[3] << "," << refFts[4] << "," << refFts[5] << ")"
                << "  |p|="    << refP
                << "  q/p="    << qOverPref0 << "\n"
                << "  track.parameters() = " << track.parameters() << std::endl;
    }

    // ---- Two-hypothesis charge fit in one pass --------------------------
    // Snapshot the seed state (spatial parts fixed; refFts[6] set per attempt)
    // and run the whole fit as a lambda parameterised by the seed-charge sign.
    // The nominal attempt runs first; if its clamp caught a q/p sign crossing
    // (nChargeFlipProtect>0 -- the fit "wanted" the opposite charge, the
    // signature of a high-p charge mis-ID) the opposite hypothesis is fit and
    // the lower-chi2 result is kept. All fit products are written to member
    // branch buffers, so the members hold the winning attempt at tree->Fill.
    const Matrix<double, 7, 1> refFtsSeed = refFts;
    ++fitAttempted_;
    auto runHypothesis = [&](int effSeedSign) -> std::pair<bool, double> {
    // reset the mutable reference state to the seed and apply this attempt's charge
    refFts = refFtsSeed;
    refFts[6] = effSeedSign * track.charge();
    const std::string g4PartName = ana_hitanalyzer::g4ParticleName(trackParticleName_, effSeedSign * track.charge());

    if (dopca) {
      // enforce that reference state is really a consistent PCA to the beamline (by adjusting the reference position if needed)
      const Matrix<double, 5, 1> statepca = cart2pca(refFts, *bsH);
      refFts = pca2cart(statepca, *bsH);
    }
    if (dumpKaonQp) {
      const double refP = refFts.segment<3>(3).norm();
      std::cout << "  refFts(pca) : charge[6]=" << refFts[6]
                << "  |p|=" << refP
                << "  q/p=" << refFts[6]/std::max(refP, 1e-12) << std::endl;
    }

    std::vector<Matrix<double, 7, 1>> layerStates;
    layerStates.reserve(nhits);
    
    bool valid = true;
    bool stepClampedThisFit = false;
    // ---- chi2-based (Armijo) backtracking bookkeeping, per fit attempt ----
    // chisq0valPrev is the chi2 at the linearization point of the PREVIOUS
    // iteration; predDecrPrev is the quadratic model's predicted chi2 change
    // for the step that was actually applied there (already including any
    // clamp/damping scale). stepScaleApplied accumulates that scale within an
    // iteration. nChi2Bt is the halving count since the last accepted step;
    // nChi2BtTotal is a per-fit budget that guarantees termination (a backtrack
    // redoes the same iteration index and so does not consume the niters budget).
    bool stepBacktrackedThisFit = false;
    unsigned int nChi2Bt = 0;
    unsigned int nChi2BtTotal = 0;
    double chisq0valPrev = std::numeric_limits<double>::quiet_NaN();
    double predDecrPrev = 0.;
    double stepScaleApplied = 1.;
    nChargeFlipProtect = 0;
    if (debugPerIterDump_) {
      chisqval_iter.clear();
      edmval_iter.clear();
      edmvalref_iter.clear();
      deltachisqval_iter.clear();
    }

    const bool islikelihood = false;

    

    std::vector<double> localxsmearedsim(hits.size(), 0.);
    
    
    double chisqvalold = std::numeric_limits<double>::max();
    
    const bool anomDebug = false;
    
    const unsigned int niters = (dogen && !dolocalupdate) || (dogen && fitFromSimParms_) ? 1 : nIters_;
    
    // CGF IRLS (CVH_CGF_QOP=3): the realised q/p noise of each block from
    // the PREVIOUS solve, r_b^(t) = (Ffull dxfull)_b. Persists across
    // Gauss-Newton iterations; empty (hence r = 0) on the first, which is
    // exactly the t = 0 case of the derivation. Sized on first use.
    //
    // Carrying it across iterations is legitimate even though the reference
    // is regenerated every iteration: the reference update propagates a seed
    // shift through all layers, and a noiseless shift has (F delta)_b = 0 by
    // construction, so (F dx)_b is INVARIANT under it. It is precisely the
    // part of dx the reference cannot absorb.
    std::vector<double> cgfRprev;
    std::vector<unsigned int> cgfBlkRow;
    // Cached block weights and score tables (CVH_CGF_QOP_REFRESH). `I` and the
    // score table are properties of the block's DISTRIBUTION, so freezing them
    // after the first sweep is standard IRLS -- fixed weights, moving centre.
    // At the fixed point the stationarity condition is
    //     hit_grad + F^T psi(r) = 0
    // which contains no `I` at all, so `I` is a PRECONDITIONER: freezing it
    // changes the path and cannot change the fixed point. That makes
    // cached-vs-uncached a genuine path-independence test rather than a
    // model change (section 33).
    std::vector<double> cgfCacheQ, cgfCacheSig;
    std::vector<cvhcgf::Result> cgfCacheRes;

    // per-iteration accumulators for the reference energy loss; the CONVERGED
    // iteration's values are the ones exported (reset at the top of each)
    double dErefIter = 0.;
    double maxFracLossIter = 0.;

    for (unsigned int iiter=0; iiter<niters; ++iiter) {
      dErefIter = 0.;
      maxFracLossIter = 0.;
      // Linearization snapshot taken BEFORE the reference/layer-state update
      // below applies dxfull, so a failed sufficient-decrease test can restore
      // it and redo the iteration with a halved step.
      const Matrix<double, 7, 1> refFtsSnap = refFts;
      std::vector<Matrix<double, 7, 1>> layerStatesSnap;
      if (stepBacktracking_) {
        layerStatesSnap = layerStates;
      }
      stepScaleApplied = 1.;
      if (debugprintout_) {
        std::cout<< "iter " << iiter << std::endl;
      }

      // Iteration marker for the CGF diagnostic stream (CVH_CGF_QOP >= 2).
      // The propagator emits one `### CVHCGF` line per leg with no idea which
      // track or iteration it belongs to; without this marker the only way to
      // group them is to guess from the ordering, which breaks the moment a
      // propagation fails or two tracks cross similar radii.
      {
        if (cvhcgf::cgfQoPMode() >= 2) {
          std::cout << "### CVHITER iiter=" << iiter << " itrack=" << itrack
                    << " pt=" << trackPt << " eta=" << trackEta << " q=" << track.charge()
                    << std::endl;
        }
      }

// std::cout<< "iter " << iiter << std::endl;
            
      hitidxv.clear();
      hitidxv.reserve(nvalid);

      if (doKinkFinder_) {
        kinkConsIdx.clear();
        kinkConsIdx.reserve(nhits);
        kinkStepR.clear();
        kinkStepR.reserve(nhits);
        kinkStepZ.clear();
        kinkStepZ.reserve(nhits);
      }

      if (iiter == 0) {
        dxrecgen.clear();
        dxrecgen.reserve(nvalid);
        
        dyrecgen.clear();
        dyrecgen.reserve(nvalid);
        
        dxsimgen.clear();
        dxsimgen.reserve(nvalid);
        
        dysimgen.clear();
        dysimgen.reserve(nvalid);

        dxsimgenconv.clear();
        dxsimgenconv.reserve(nvalid);

        dysimgenconv.clear();
        dysimgenconv.reserve(nvalid);
        
        dxsimgenlocal.clear();
        dxsimgenlocal.reserve(nvalid);
        
        dysimgenlocal.clear();
        dysimgenlocal.reserve(nvalid);
        
        dxrecsim.clear();
        dxrecsim.reserve(nvalid);
        
        dyrecsim.clear();
        dyrecsim.reserve(nvalid);
        
        dxreccluster.clear();
        dxreccluster.reserve(nvalid);
        
        dyreccluster.clear();
        dyreccluster.reserve(nvalid);
        
        dxerr.clear();
        dxerr.reserve(nvalid);
        
        dyerr.clear();
        dyerr.reserve(nvalid);
        
        clusterSize.clear();
        clusterSize.reserve(nvalid);
        
        clusterSizeX.clear();
        clusterSizeX.reserve(nvalid);
        
        clusterSizeY.clear();
        clusterSizeY.reserve(nvalid);
        
        clusterCharge.clear();
        clusterCharge.reserve(nvalid);
        
        clusterChargeBin.clear();
        clusterChargeBin.reserve(nvalid);
        
        clusterOnEdge.clear();
        clusterOnEdge.reserve(nvalid);
        
        clusterProbXY.clear();
        clusterProbXY.reserve(nvalid);
        
        clusterSN.clear();
        clusterSN.reserve(nvalid);
        
        simlocalqop.clear();
        simlocaldxdz.clear();
        simlocaldydz.clear();
        simlocalx.clear();
        simlocaly.clear();

        simlocalqop.reserve(nvalid);
        simlocaldxdz.reserve(nvalid);
        simlocaldydz.reserve(nvalid);
        simlocalx.reserve(nvalid);
        simlocaly.reserve(nvalid);
        
        simlocalqopprop.clear();
        simlocaldxdzprop.clear();
        simlocaldydzprop.clear();
        simlocalxprop.clear();
        simlocalyprop.clear();

        simlocalqopprop.reserve(nvalid);
        simlocaldxdzprop.reserve(nvalid);
        simlocaldydzprop.reserve(nvalid);
        simlocalxprop.reserve(nvalid);
        simlocalyprop.reserve(nvalid);
        
        landauDelta.clear();
        landauW.clear();

        landauDelta.reserve(nvalid);
        landauW.reserve(nvalid);

        localqop.clear();
        localdxdz.clear();
        localdydz.clear();
        localx.clear();
        localy.clear();

        localqop.reserve(nvalid);
        localdxdz.reserve(nvalid);
        localdydz.reserve(nvalid);
        localx.reserve(nvalid);
        localy.reserve(nvalid);
        
        localqoperr.clear();
        localdxdzerr.clear();
        localdydzerr.clear();
        localxerr.clear();
        localyerr.clear();

        localqoperr.reserve(nvalid);
        localdxdzerr.reserve(nvalid);
        localdydzerr.reserve(nvalid);
        localxerr.reserve(nvalid);
        localyerr.reserve(nvalid);
        
        hitlocalx.clear();
        hitlocaly.clear();
        
        hitlocalx.reserve(nvalid);
        hitlocaly.reserve(nvalid);
        
        localqop_iter.assign(nvalid, -99.);
        localdxdz_iter.assign(nvalid, -99.);
        localdydz_iter.assign(nvalid, -99.);
        localx_iter.assign(nvalid, -99.);
        localy_iter.assign(nvalid, -99.);

        dxrecgen_iter.assign(nvalid, -99.);
        dyrecgen_iter.assign(nvalid, -99.);
        
        localqoperralt.assign(nvalid, -99.);

        localphi.clear();
        localphi.reserve(nvalid);

        hitphi.clear();
        hitphi.reserve(nvalid);
        
        dEpred.clear();
        dEpred.reserve(nvalid);
        
        dE.clear();
        dE.reserve(nvalid);
        
        sigmadE.clear();
        sigmadE.reserve(nvalid);
        
        E.clear();
        E.reserve(nvalid);
        
        Epred.clear();
        Epred.reserve(nvalid);

        stripsToEdge.clear();
        stripsToEdge.reserve(nvalid);

        simHitNCand.clear();
        simHitNCand.reserve(nvalid);
        hitDetId.clear();
        hitDetId.reserve(nvalid);
        hitUProj.clear();
        hitUProj.reserve(nvalid);
        hitStripRec.clear();
        hitStripRec.reserve(nvalid);
        hitStripSim.clear();
        hitStripSim.reserve(nvalid);
        hitFirstStrip.clear();
        hitFirstStrip.reserve(nvalid);
        hitPitch.clear();
        hitPitch.reserve(nvalid);
        hitThickness.clear();
        hitThickness.reserve(nvalid);

      }      
      
      rfull = VectorXd::Zero(ncons);
      Ffull = MatrixXd::Zero(ncons, nstateparms);
      cgfBlkRow.clear();
      Jfull = MatrixXd::Zero(ncons, npars);
      Vinvfull = MatrixXd::Zero(ncons, ncons);
      Vinvfullalt = MatrixXd::Zero(ncons, ncons);
      
      dVs.clear();
      resvalidhit_.clear();
      resfamily_.clear();
      residxs.clear();
      ioniurbanidx.clear();
      ioniurbanv.clear();
      ioniqscaleidx.clear();
      ioniqscalev.clear();
      radstepidx.clear();
      radstepv.clear();
      radstepspecv.clear();
      msmoliidx.clear();
      msmoliv.clear();
      reseigidx.clear();
      reseigv.clear();
      resinfv.clear();
      resinfvarv.clear();
      reshitidx.clear();
      reshitcls.clear();
      cfhitclsv.clear();
      cfhitvv.clear();
      resinfcov = 0.;
      resinfcovhit = 0.f;
      resinfcovgrp = 0.f;
      resinfbv.clear();
      resblockrng.clear();
      resglobidx.clear();
      
      validdxeigjac = MatrixXd::Zero(2*nvalid, nstateparms);
      
      
      double chisq0val = 0.;
      
      dxpxb1 = -99.;
      dypxb1 = -99.;
      dxttec9rphi = -99.;
      dxttec9stereo = -99.;
      dxttec4rphi = -99.;
      dxttec4stereo = -99.;
      
      dxttec4rphisimgen = -99.;
      dyttec4rphisimgen = -99.;
      dxttec4rphirecsim = -99.;
      
      dxttec9rphisimgen = -99.;
      dyttec9rphisimgen = -99.;
      
      simlocalxref = -99.;
      simlocalyref = -99.;
      
      unsigned int ivalidhit = 0;
      
      unsigned int icons = 0;
      unsigned int iparm = 0;

      if (iiter > 0) {
        //update current state from reference point state (errors not needed beyond first iteration)
        
        auto const& dxlocal = dxfull.head<5>();

        
        if (dopca) {
          // position transformation from track-beamline pca to cartesian is non-linear so do it exactly
          const Matrix<double, 5, 1> statepca = cart2pca(refFts, *bsH);
          const Matrix<double, 5, 1> statepcaupd = statepca + dxlocal;

          refFts = pca2cart(statepcaupd, *bsH);

        }
        else {
          // position transformation from curvilinear to cartesian is linear so just use the jacobian

          const Matrix<double, 6, 5> jac = curv2cartJacobianAltD(refFts);
          const Matrix<double, 6, 1> globupd = refFts.head<6>() + jac*dxlocal;

          const double qbp = refFts[6]/refFts.segment<3>(3).norm();
          const double lam = std::atan(refFts[5]/std::sqrt(refFts[3]*refFts[3] + refFts[4]*refFts[4]));
          const double phi = std::atan2(refFts[4], refFts[3]);

          const double qbpupd = qbp + dxlocal(0);
          const double lamupd = lam + dxlocal(1);
          const double phiupd = phi + dxlocal(2);

          const double charge = std::copysign(1., qbpupd);
          const double pupd = std::abs(1./qbpupd);

          const double pxupd = pupd*std::cos(lamupd)*std::cos(phiupd);
          const double pyupd = pupd*std::cos(lamupd)*std::sin(phiupd);
          const double pzupd = pupd*std::sin(lamupd);

          refFts.head<3>() = globupd.head<3>();
          refFts[3] = pxupd;
          refFts[4] = pyupd;
          refFts[5] = pzupd;
          refFts[6] = charge;
        }
      }
      
      Matrix<double, 5, 5> ref2curvjac = dopca ? pca2curvJacobianD(refFts, field, *bsH) : Matrix<double, 5, 5>::Identity();

      Matrix<double, 7, 1> updtsos = refFts;

      Matrix<double, 7, 1> propfromtsos = refFts;

      Matrix<double, 5, 5> Qtot = Matrix<double, 5, 5>::Zero();

      float e = genpart == nullptr ? -99. : std::sqrt(genpart->momentum().mag2() + trackmass*trackmass);
      float epred = std::sqrt(refFts.segment<3>(3).squaredNorm() + trackmass*trackmass);
      
      
      if (bsConstraint_) {
        // apply beamspot constraint
        // TODO add residual corrections for beamspot parameters?
        // TODO make this a 2D constraint?
        
        //TODO convert this to Ffull/Vinvfull parameterization
        
        constexpr unsigned int nlocalcons = 3;
        constexpr unsigned int nlocalstate = 5;
        constexpr unsigned int nlocal = nlocalstate;
        constexpr unsigned int localstateidx = 0;
        constexpr unsigned int fullstateidx = 0;
        
        const Matrix<double, 6, 5> jac = dopca ? pca2cartJacobianD(refFts, *bsH) : curv2cartJacobianAltD(refFts);
        
        const double sigb1 = bsH->BeamWidthX();
        const double sigb2 = bsH->BeamWidthY();
        const double sigb3 = bsH->sigmaZ();
        const double dxdz = bsH->dxdz();
        const double dydz = bsH->dydz();
        const double x0 = bsH->x0();
        const double y0 = bsH->y0();
        const double z0 = bsH->z0();
        
        
        // covariance matrix of luminous region in global coordinates
        // taken from https://github.com/cms-sw/cmssw/blob/abc1f17b230effd629c9565fb5d95e527abcb294/RecoVertex/BeamSpotProducer/src/FcnBeamSpotFitPV.cc#L63-L90

        // TODO check consistency of this parameterization

        // FIXME xy correlation is not stored and assumed to be zero
        const double corrb12 = 0.;
        
        const double varb1 = sigb1*sigb1;
        const double varb2 = sigb2*sigb2;
        const double varb3 = sigb3*sigb3;
        
        Matrix<double, nlocalcons, nlocalcons> covBS = Matrix<double, nlocalcons, nlocalcons>::Zero();
        // parametrisation: rotation (dx/dz, dy/dz); covxy
        covBS(0,0) = varb1;
        covBS(1,0) = covBS(0,1) = corrb12*sigb1*sigb2;
        covBS(1,1) = varb2;
        covBS(2,0) = covBS(0,2) = dxdz*(varb3-varb1)-dydz*covBS(1,0);
        covBS(2,1) = covBS(1,2) = dydz*(varb3-varb2)-dxdz*covBS(1,0);
        covBS(2,2) = varb3;
        
        Matrix<double, nlocalcons, 1> dbs0;
        dbs0[0] = refFts[0] - x0;
        dbs0[1] = refFts[1] - y0;
        dbs0[2] = refFts[2] - z0;
        
        const Matrix<double, nlocalcons, nlocalstate> Fbs = jac.topRows<3>();
        const Matrix<double, nlocalcons, nlocalcons> covBSinv = covBS.inverse();
        
        rfull.segment<nlocalcons>(icons) = dbs0;
        Ffull.block<nlocalcons, nlocalstate>(icons, fullstateidx) = Fbs;
        Vinvfull.block<nlocalcons, nlocalcons>(icons, icons) = covBSinv;
        Vinvfullalt.block<nlocalcons, nlocalcons>(icons, icons) = covBSinv;
        
        icons += nlocalcons;
      }
      
      
      for (unsigned int ihit = 0; ihit < hits.size(); ++ihit) {
// std::cout << "iiter = " << iiter << " ihit " << ihit << std::endl;

        auto const& hit = hits[ihit];
        
        // auto const &surface = surfacemapD_.at(hit->geographicalId());
        GloballyPositioned<double> surface = surfacemapD_.at(hit->geographicalId());

        // check propagation direction to choose correct bfield and material parameters
        // when propagating inside-out the parameters correspond to the target module
        // when propagating outside-in (e.g. for cosmics) they correspond to the source module
        // For the first hit always use the target module
        // (only relevant for cosmics)

        bool sourceParms = false;
        if (iscosmic && ihit > 0) {
          auto const &lzhat = surface.rotation().z();
          auto const &pos = surface.position();

          const double zdotpos = lzhat.x()*pos.x() + lzhat.y()*pos.y() + lzhat.z()*pos.z();

          const Point3DBase<double, GlobalTag> globalpos(updtsos[0], updtsos[1], updtsos[2]);
          const Point3DBase<double, LocalTag> localpos = surface.toLocal(globalpos);

          sourceParms = zdotpos*localpos.z() > 0.;
        }

// std::cout << "ihit = " << ihit << " sourceParms = " << sourceParms << std::endl;

        auto const& prophit = sourceParms ? hits[ihit - 1] : hit;
        const uint32_t gluedidprop = trackerTopology->glued(prophit->geographicalId());
        const bool isgluedprop = gluedidprop != 0;
        const DetId propdetid = isgluedprop ? DetId(gluedidprop) : prophit->geographicalId();

        const uint32_t gluedid = trackerTopology->glued(hit->geographicalId());
        const bool isglued = gluedid != 0;
        const DetId parmdetid = isglued ? DetId(gluedid) : hit->geographicalId();
        const DetId aligndetid = alignGlued_ ? parmdetid : hit->geographicalId();
        
        const unsigned int elossglobalidx =
            globalMaterialModel_ ? 0 : detidparms.at(std::make_pair(7, propdetid));
        const unsigned int msglobalidx = dores ? detidparms.at(std::make_pair(10, propdetid)) : 0;
        const unsigned int ioniglobalidx = dores ? detidparms.at(std::make_pair(11, propdetid)) : 0;

        // 3D field correction at the start of the propagation step. The
        // per-mode (Bx, By, Bz) basis values at the same point feed the
        // chain-rule scaling for the transport-Jacobian dBx/dBy/dBz columns
        // (cols 5,6,7 of the 5x9 transportJacobianBxByBzD).
        const GlobalPoint propStartPos(updtsos[0], updtsos[1], updtsos[2]);
        // Per-step mode: the provider applies the correction inside the
        // propagator (dB argument stays zero for FD perturbations) and the
        // per-leg basis samples are not needed.
        const Eigen::Vector3d dB = perStepFieldModes_
            ? Eigen::Vector3d::Zero()
            : fieldCorrection_->getCorrectionAt(propStartPos, corparms_);
        std::vector<double> dBxPerMode, dByPerMode, dBzPerMode;
        if (!perStepFieldModes_) {
          fieldCorrection_->getBxBasisAt(propStartPos, dBxPerMode);
          fieldCorrection_->getByBasisAt(propStartPos, dByPerMode);
          fieldCorrection_->getBzBasisAt(propStartPos, dBzPerMode);
        }

        // Global material model: the leg-constant dxi is zero; the per-step
        // group values k_g (synced from corparms_ into the model at the top
        // of the fit) are applied by the propagator's provider path.
        const double dxival = globalMaterialModel_ ? 0. : corparms_[elossglobalidx];
        const double dmsval = dores ? corparms_[msglobalidx] : dxival;
        const double dionival = dores ? corparms_[ioniglobalidx] : dxival;
        
        const PSimHit *simhit = nullptr;
        int simhitNCandVal = 0;
        std::vector<const PSimHit*> simhitCands;
        
        if (doSim_) {
          // The track can leave MORE THAN ONE PSimHit on a module -- a curling
          // low-momentum hadron, or a re-entry after a large-angle nuclear
          // elastic scatter. Taking the first candidate produced a 0.126 %
          // population of proton strip hits with |(rec-sim)/sigma| > 10 and a
          // median |rec - sim| of 1.35 mm, which would be read as a hit
          // resolution tail; the kaon gun shows 40x less of it, so it is
          // proton kinematics and not a defect of the hit. The CHOICE among
          // candidates is deferred to the point where the propagated local
          // position exists (see the simhitCands re-selection below); the
          // first candidate is kept here so the fitFromSimParms / simhitdebug
          // paths behave exactly as before.
          //
          // Species from the CONFIGURED gen-match hypothesis, not a hardcoded
          // muon. With |particleType|==13 the sim-hit machinery silently
          // no-ops on the kaon/pion/proton guns: simhit stays null, so
          // dxrecsim/dyrecsim are -99 and fitSimHitPositions quietly falls
          // back to reco positions (usesimpos is && simhit != nullptr).
          for (auto const& simhith : simHits) {
            for (const PSimHit& simHit : *simhith) {
              if (simHit.detUnitId() == hit->geographicalId() && int(simHit.trackId()) == simtrackid && std::abs(simHit.particleType()) == genMatchPdgId_) {
                simhitCands.push_back(&simHit);
              }
            }
          }
          simhitNCandVal = int(simhitCands.size());
          if (!simhitCands.empty()) {
            simhit = simhitCands.front();
          }
        }

        // const bool simhitdebug = true;
        const bool simhitdebug = false;

        if (simhitdebug && !simhit) {
          std::cout << "WARNING: no simhit for debugging scenario, abort!\n";
          valid = false;
          break;
        }

        if ((simhitdebug || fitFromSimParms_) && simhit) {
          //move surface to plane of sim entry point for consistency in debugging vs sim state

          const double lz = simhit->entryPoint().z();
          const Vector3DBase<double, LocalTag> dlocal(0., 0., lz);
          const Vector3DBase<double, GlobalTag> dglobal = surface.toGlobal(dlocal);
          surface.move(dglobal);
        }

// std::cout << "iiter = " << iiter << " ihit = " << ihit << " updtsos:\n" << updtsos << std::endl;

        // auto const &propresult = g4prop->propagateGenericWithJacobianAltD(updtsos, surface, dbetaval, dxival, dradval);
        // auto const &propresult = g4prop->propagateGenericWithJacobianAltD(propfromtsos, surface, dbetaval, dxival, dradval);

        // Save the input state so the FD closure block can re-run
        // the propagation with a perturbed dB starting from the same point.
        const Eigen::Matrix<double, 7, 1> propInputState =
            simhitdebug ? propfromtsos : updtsos;

        // cached block weight for this leg, if we have one and are not
        // refreshing on this iteration
        {
          // `cvhcgf` is the single reader; the maker must not carry its own
          // copy of the schedule, or the propagator and the fit could disagree
          // about which estimator is running.
          const int cgfRefresh = cvhcgf::cgfQoPRefresh();
          const bool refreshNow = (iiter == 0) ||
                                  (cgfRefresh > 0 && (iiter % cgfRefresh) == 0);
          if (!refreshNow && ihit < cgfCacheQ.size() && cgfCacheQ[ihit] > 0.) {
            g4prop->setCgfOverride(cgfCacheQ[ihit]);
          }
        }
        auto const &propresult = simhitdebug
            ? g4prop->propagateGenericWithJacobianAltD(propfromtsos, surface, dB, dxival, dmsval, dionival, -1., g4PartName,
                                                       matModel_.get(), matModel_ ? &groupJacs_ : nullptr,
                                                       (matModel_ && doRes_ && exportMaterialNoise_) ? &groupQs_ : nullptr,
                                                       fieldModeProvider_.get(),
                                                       fieldModeProvider_ ? &modeJacs_ : nullptr)
            : g4prop->propagateGenericWithJacobianAltD(updtsos,      surface, dB, dxival, dmsval, dionival, -1., g4PartName,
                                                       matModel_.get(), matModel_ ? &groupJacs_ : nullptr,
                                                       (matModel_ && doRes_ && exportMaterialNoise_) ? &groupQs_ : nullptr,
                                                       fieldModeProvider_.get(),
                                                       fieldModeProvider_ ? &modeJacs_ : nullptr);


        if (simhitdebug && simhit) {
          // reset propagation start point to sim hit for debugging purposes
          auto const &surfaceideal = surfacemapIdealD_.at(hit->geographicalId());
          const Point3DBase<double, GlobalTag> pos = surfaceideal.toGlobal(simhit->entryPoint());
          const Vector3DBase<double, GlobalTag> mom = surfaceideal.toGlobal(simhit->momentumAtEntry());

          propfromtsos[0] = pos.x();
          propfromtsos[1] = pos.y();
          propfromtsos[2] = pos.z();
          propfromtsos[3] = mom.x();
          propfromtsos[4] = mom.y();
          propfromtsos[5] = mom.z();
          propfromtsos[6] = genpart->charge();
        }

        if (!std::get<0>(propresult)) {
          std::cout << "Abort: Propagation Failed!"
                    << " iiter = " << iiter << " ihit = " << ihit
                    << " seed: q=" << track.charge() << " pt=" << track.pt()
                    << " eta=" << track.eta() << std::endl;
          ++fitFailProp_;
          {
            const DetId fdet = hit->geographicalId();
            std::cout << "CVHDIAG reason=PROP ihit=" << ihit << " nhit=" << hits.size()
                      << " detid=" << fdet.rawId() << " subdet=" << fdet.subdetId()
                      << " layer=" << trackerTopology->layer(fdet)
                      << " sx=" << surface.position().x() << " sy=" << surface.position().y()
                      << " sz=" << surface.position().z()
                      << " trkpt=" << trackPt << " trketa=" << trackEta << " trkphi=" << trackPhi
                      << " iiter=" << iiter << std::endl;
          }
          valid = false;
          break;
        }
        
        updtsos = std::get<1>(propresult);

        // Reference energy loss of this propagation, and the running maximum
        // of its fractional size. Taken from the propagator's input/output
        // STATES rather than from any dE/dx model call, so it stays correct
        // whatever scales the loss (CVH_DEDX_SCALE, the material model's k_g,
        // the per-module dxi), and formed before any local state update, so
        // only the propagation contributes. Same construction as the
        // two-track maker's `dErefarr`.
        {
          const double m2 = trackmass * trackmass;
          const double pIn = propInputState.segment<3>(3).norm();
          const double eIn = std::sqrt(pIn * pIn + m2);
          const double eOut = std::sqrt(updtsos.segment<3>(3).squaredNorm() + m2);
          dErefIter += eIn - eOut;
          if (pIn > 0.) {
            maxFracLossIter = std::max(maxFracLossIter, (eIn - eOut) / pIn);
          }
        }

        const Matrix<double, 5, 5> Qcurv = std::get<2>(propresult);
        // record the (possibly just-computed) block weight and score table
        if (g4prop != nullptr && g4prop->cgfBlockValid()) {
          if (cgfCacheQ.size() <= ihit) {
            cgfCacheQ.resize(ihit + 1, -1.);
            cgfCacheSig.resize(ihit + 1, 0.);
            cgfCacheRes.resize(ihit + 1);
          }
          cgfCacheQ[ihit] = Qcurv(0, 0);
          cgfCacheSig[ihit] = g4prop->cgfBlockSigma();
          cgfCacheRes[ihit] = g4prop->cgfBlock();
        }
        const Matrix<double, 5, 9> FdFmcurv = std::get<3>(propresult);
        const double dEdxlast = std::get<4>(propresult);
        const Matrix<double, 5, 5> dQMScurv = std::get<5>(propresult);
        const Matrix<double, 5, 5> dQIcurv = std::get<6>(propresult);
// const Matrix<double, 5, 5> dQcurv = Qcurv;
        const double deltaTotal = std::get<7>(propresult);
        const double wTotal = std::get<8>(propresult);

        // Global material model Phase A validations (see member docs).
        if (matModel_) {
          // V2: the per-group columns must sum to the integrated dxi column
          Eigen::Matrix<double, 5, 1> gsum = Eigen::Matrix<double, 5, 1>::Zero();
          for (auto const &gc : groupJacs_) {
            gsum += gc.second;
          }
          const double ref = FdFmcurv.col(8).norm();
          const double rel = ref > 0. ? (gsum - FdFmcurv.col(8)).norm() / ref
                                      : (gsum - FdFmcurv.col(8)).norm();
          if (rel > v2MaxRelDiff_) {
            v2MaxRelDiff_ = rel;
          }
          ++v2Checks_;

          // V1: one-shot FD closure of one group's column on the first leg
          // that actually crosses the injected group
          if (materialFDGroup_ >= 0 && !fdMatDone_) {
            auto itg = std::find_if(groupJacs_.begin(), groupJacs_.end(),
                                    [this](auto const &e) { return e.first == materialFDGroup_; });
            if (itg != groupJacs_.end() && itg->second.cwiseAbs().maxCoeff() > 0.) {
              fdMatDone_ = true;
              const Eigen::Matrix<double, 5, 1> anacol = itg->second;
              matModel_->setInjection(materialFDGroup_, materialFDEps_);
              auto const &pertres = g4prop->propagateGenericWithJacobianAltD(
                  propInputState, surface, dB, dxival, dmsval, dionival, -1., g4PartName,
                  matModel_.get(), nullptr);
              matModel_->setInjection(-1, 0.);
              if (std::get<0>(pertres)) {
                const auto &ftsNom = std::get<1>(propresult);
                const auto &ftsPert = std::get<1>(pertres);
                const double qopNom = ftsNom[6] / ftsNom.segment<3>(3).norm();
                const double qopPert = ftsPert[6] / ftsPert.segment<3>(3).norm();
                const double fd = (qopPert - qopNom) / materialFDEps_;
                std::cout << "material FD closure (V1): group " << materialFDGroup_
                          << " (" << matModel_->groupName(materialFDGroup_) << ")"
                          << " eps = " << materialFDEps_
                          << " dqop/dk analytic = " << anacol[0]
                          << " fd = " << fd
                          << " ratio fd/analytic = "
                          << (anacol[0] != 0. ? fd / anacol[0] : 0.) << std::endl;

                // V3: THE PER-GROUP PROCESS NOISE, the width counterpart of
                // V1's mean.  `setInjection(g, eps)` multiplies group g's
                // `matStepFact` by e^eps, and that factor multiplies the
                // step's MS covariance and ionization variance as well as its
                // mean loss, so
                //     Q(eps) - Q(0) = (e^eps - 1) dQ_g + O(eps^2)
                // with dQ_g exactly the block this maker now registers as the
                // parmtype-15 resolution family.  Comparing at the Q level
                // rather than end-to-end on the chi2 is deliberate: the
                // gradient/Hessian assembly downstream is the SAME code the
                // parmtype-10/11 families already go through and is validated
                // by them; what is new is only this matrix.
                auto itq = std::find_if(groupQs_.begin(), groupQs_.end(),
                                        [this](auto const &e) { return e.first == materialFDGroup_; });
                if (itq != groupQs_.end()) {
                  const Matrix<double, 5, 5> &anaQ = itq->second;
                  const Matrix<double, 5, 5> qNom = std::get<2>(propresult);
                  const Matrix<double, 5, 5> qPert = std::get<2>(pertres);
                  const double scale = std::expm1(materialFDEps_);
                  const Matrix<double, 5, 5> fdQ = (qPert - qNom) / scale;
                  const double refn = anaQ.cwiseAbs().maxCoeff();
                  const double relQ = refn > 0. ? (fdQ - anaQ).cwiseAbs().maxCoeff() / refn : -1.;
                  // and the sum rule the split has to obey exactly
                  Matrix<double, 5, 5> qsum = Matrix<double, 5, 5>::Zero();
                  for (auto const &gq : groupQs_) {
                    qsum += gq.second;
                  }
                  const Matrix<double, 5, 5> qmsi = dQMScurv + dQIcurv;
                  const double refs = qmsi.cwiseAbs().maxCoeff();
                  const double relS = refs > 0. ? (qsum - qmsi).cwiseAbs().maxCoeff() / refs : -1.;
                  std::cout << "material FD closure (V3, process noise): group " << materialFDGroup_
                            << " max|dQ_g| = " << refn
                            << " max|fd - analytic|/max|dQ_g| = " << relQ
                            << " ; sum rule max|sum_g dQ_g - (dQMS+dQI)|/max = " << relS << std::endl;
                }
              }
            }
          }
        }
        
        Matrix<double, 5, 9> FdFm = FdFmcurv;
        if (ihit == 0) {
          // extra jacobian from reference state to curvilinear potentially needed
          FdFm.leftCols<5>() = FdFmcurv.leftCols<5>()*ref2curvjac;
        }
        

        // propfromtsos = updtsos;

// //zero energy loss contribution
// FdFm.rightCols<1>() *= 0.;


        Qtot = (FdFmcurv.leftCols<5>()*Qtot*FdFmcurv.leftCols<5>().transpose()).eval();
        Qtot += Qcurv;

        if (simhitdebug) {
          Qtot = Qcurv;
        }

        const Matrix<double, 5, 5> Hm = curv2localJacobianAltelossD(updtsos, field, surface, dEdxlast, trackmass, dB);
        
        
        const float enext = simhit == nullptr ? -99. : std::sqrt(std::pow(simhit->pabs(), 2) + trackmass*trackmass) - 0.5*simhit->energyLoss();        
        const float eprednext = std::sqrt(updtsos.segment<3>(3).squaredNorm() + trackmass*trackmass);
        
        const float dEval = e > 0. && enext > 0. ? enext - e : -99.;
        const float dEpredval = eprednext - epred;
        
        e = enext;
        epred = eprednext;
        
        const float sigmadEval = std::pow(updtsos.segment<3>(3).norm(), 3)/e*std::sqrt(Qcurv(0, 0));
        

        const Matrix<double, 6, 1> localparmsprop = globalToLocal(updtsos, surface);

        Matrix<double, 6, 1> localparms = localparmsprop;
        
        // update state from previous iteration
        //momentum kink residual
        Matrix<double, 5, 1> dx0 = Matrix<double, 5, 1>::Zero();
        
        if (iiter == 0 && fitFromSimParms_) {
          if (simhit == nullptr) {
            std::cout << "ABORT: Propagating from sim parameters, but sim hit is not matched!\n";
            valid = false;
            break;
          }

          auto const &surfaceideal = surfacemapIdealD_.at(hit->geographicalId());


          // alternate version with propagation from entry state

          const Point3DBase<double, LocalTag> simlocalpos = simhit->entryPoint();
          const Vector3DBase<double, LocalTag> simlocalmom = simhit->momentumAtEntry();

          const Point3DBase<double, GlobalTag> simglobalpos = surfaceideal.toGlobal(simlocalpos);
          const Vector3DBase<double, GlobalTag> simglobalmom = surfaceideal.toGlobal(simlocalmom);

          updtsos[0] = simglobalpos.x();
          updtsos[1] = simglobalpos.y();
          updtsos[2] = simglobalpos.z();
          updtsos[3] = simglobalmom.x();
          updtsos[4] = simglobalmom.y();
          updtsos[5] = simglobalmom.z();
          updtsos[6] = genpart->charge();

          if (false) {
            auto const &propresultsim = g4prop->propagateGenericWithJacobianAltD(updtsos, surface, dB, dxival, dmsval, dionival, -1., g4PartName);

            if (!std::get<0>(propresultsim)) {
              std::cout << "Abort: Sim state Propagation Failed!" << std::endl;
              valid = false;
              break;
            }

            updtsos = std::get<1>(propresultsim);
          }

          localparms = globalToLocal(updtsos, surface);

          dx0 = (localparms - localparmsprop).head<5>();

        }
        
        
        if (dolocalupdate) {
          if (iiter==0) {
            layerStates.push_back(updtsos);
          }
          else {
            //current state from previous state on this layer
            //save current parameters

            Matrix<double, 7, 1>& oldtsos = layerStates[ihit];
            const Matrix<double, 5, 5> Hold = curv2localJacobianAltelossD(oldtsos, field, surface, dEdxlast, trackmass, dB);
            const Matrix<double, 5, 1> dxlocal = Hold*dxfull.segment<5>(5*(ihit+1));

            localparms = globalToLocal(oldtsos, surface);

            localparms.head<5>() += dxlocal;

            oldtsos = localToGlobal(localparms, surface);

            updtsos = oldtsos;

            dx0 = (localparms - localparmsprop).head<5>();
          }
        }
        
        Matrix<double, 6, 1> localparmsalignprop = localparms;
        Matrix<double, 7, 1> updtsosalign = updtsos;
        
        if (fitFromSimParms_) {
        // re-propagate to unmodified surface

          
          auto const &surfaceideal = surfacemapIdealD_.at(hit->geographicalId());


          // alternate version with propagation from entry state

          const Point3DBase<double, LocalTag> simlocalpos = simhit->entryPoint();
          const Vector3DBase<double, LocalTag> simlocalmom = simhit->momentumAtEntry();

          const Point3DBase<double, GlobalTag> simglobalpos = surfaceideal.toGlobal(simlocalpos);
          const Vector3DBase<double, GlobalTag> simglobalmom = surfaceideal.toGlobal(simlocalmom);

          Matrix<double, 7, 1> simtsos;
          simtsos[0] = simglobalpos.x();
          simtsos[1] = simglobalpos.y();
          simtsos[2] = simglobalpos.z();
          simtsos[3] = simglobalmom.x();
          simtsos[4] = simglobalmom.y();
          simtsos[5] = simglobalmom.z();
          simtsos[6] = genpart->charge();          
          
          auto const &surfacealign = surfacemapD_.at(hit->geographicalId());
          auto const &propresultsalign = g4prop->propagateGenericWithJacobianAltD(simtsos, surfacealign, dB, dxival, dmsval, dionival, -1., g4PartName);

          if (!std::get<0>(propresultsalign)) {
            std::cout << "WARNING propagation for alignment failed!\n";
            valid = false;
            break;
          }
          
          updtsosalign = std::get<1>(propresultsalign);
          
          localparmsalignprop = globalToLocal(updtsosalign, surfacealign);
        }

        // Now that the PROPAGATED local position exists, choose among the
        // candidate sim hits on this module the one closest to it. The
        // tie-break is measurement-independent, so it cannot pull rec - sim
        // toward zero the way breaking it on the reco position would.
        if (simhitCands.size() > 1) {
          double bestd2 = std::numeric_limits<double>::max();
          for (const PSimHit *cand : simhitCands) {
            const double dxp = cand->localPosition().x() - localparmsalignprop[3];
            const double dyp = cand->localPosition().y() - localparmsalignprop[4];
            const double d2 = dxp*dxp + dyp*dyp;
            if (d2 < bestd2) {
              bestd2 = d2;
              simhit = cand;
            }
          }
        }

        //TODO optimize this without ternary functions
        
        Matrix<double, 5, 5> Hp = Hm;
        
        if (dolocalupdate) {
          Hp = curv2localJacobianAltelossD(updtsos, field, surface, dEdxlast, trackmass, dB);
        }

        Matrix<double, 5, 5> Q = Qcurv;
        if (dolocalupdate) {
          Q = Hm*Qcurv*Hm.transpose();
        }
        
        const Matrix<double, 5, 5> Qinv = Q.inverse();
        
        Matrix<double, 5, 5> dQMS = dQMScurv;
        Matrix<double, 5, 5> dQI = dQIcurv;
        if (dolocalupdate) {
          dQMS = Hm*dQMScurv*Hm.transpose();
          dQI = Hm*dQIcurv*Hm.transpose();
        }
        
        {
          constexpr unsigned int nlocalcons = 5;
          constexpr unsigned int nlocalstateparms = 5;
          // Field block expands to nFieldModes columns (one per scalar-potential
          // mode), each scaled by its Bz basis value at the propagation start.
          const unsigned int nlocalbfield = nFieldModes;
          // Eloss block: one column per material group (global model) or the
          // single per-module dxi column (legacy).
          const unsigned int nlocaleloss = globalMaterialModel_ ? nMatGroups : 1;
          const unsigned int nlocalparms = nlocalbfield + nlocaleloss;

          const unsigned int fullstateidx = 5*ihit;
          const unsigned int fullparmidx = iparm;

          if (doKinkFinder_) {
            kinkConsIdx.push_back(icons);
            kinkStepR.push_back(std::sqrt(updtsos[0]*updtsos[0] + updtsos[1]*updtsos[1]));
            kinkStepZ.push_back(updtsos[2]);
            if (kinkInjectLayer_ >= 0 && int(ihit) == kinkInjectLayer_) {
              // closure test: mimic a decay at this step by offsetting the
              // material residual (state minus propagated) by the injected
              // kink; the score-test scan must recover deltahat = injected
              dx0[0] += kinkInjectDqop_;
              dx0[1] += kinkInjectDxdz_;
              dx0[2] += kinkInjectDydz_;
            }
          }

          // ------------------------------------------------------------------
          // IRLS RE-CENTRING of the q/p process-noise row (CVH_CGF_QOP=3).
          // OFF unless the mode is exactly 3; modes 0/1/2 leave dx0 untouched.
          //
          // The surrogate replaces the block's -ln p by
          //     1/2 (r - mu)^T I (r - mu),   mu = r - psi(r)/I
          // so the residual the least-squares sees is (r - mu) = psi(r)/I,
          // evaluated at the current iterate. With the outer Gauss-Newton loop
          // already iterating, taking `current` = this iterate is exactly
          // FISHER SCORING (the increment solves I delta = psi), and at the
          // fixed point the surrogate score equals the true score EXACTLY.
          //
          // In the propagator's standardized units z = r/sigma with
          // sigma^2 = the leg's nominal Gaussian q/p variance:
          //     r_used = sigma * (1/I_z) * psi_z(r / sigma)
          // and the weight Q(0,0) = sigma^2 * (1/I_z) is already substituted.
          //
          // GAUSSIAN-LIMIT IDENTITY, and it is the wiring control: for a
          // Gaussian block psi_z(z) = z / (1/I_z), so r_used = r identically.
          // CVH_CGF_QOP_GAUSSPSI=1 substitutes exactly that psi and MUST
          // reproduce mode 1 bit for bit -- which separates "are the units
          // right" from "is the score right".
          //
          // Only component 0 is touched. dx0 is the local 5-parameter
          // residual whose first component is q/p in both the local and the
          // curvilinear bases, which is why the standardization by the
          // CURVILINEAR sigma is the matching one.
          const int cgfMode = cvhcgf::cgfQoPMode();
          static const bool cgfGaussPsi = (getenv("CVH_CGF_QOP_GAUSSPSI") != nullptr);
          if (cgfMode == 3) {
            // Block index within this track, in layer order -- the same order
            // every iteration, which is what makes cgfRprev addressable.
            const unsigned int iblk = cgfBlkRow.size();
            cgfBlkRow.push_back(icons);
            const double rprev = (iblk < cgfRprev.size()) ? cgfRprev[iblk] : 0.;
            if (ihit < cgfCacheRes.size() && cgfCacheRes[ihit].ok) {
              const cvhcgf::Result &cgfb = cgfCacheRes[ihit];
              const double csig = cgfCacheSig[ihit];
              if (csig > 0. && cgfb.invFisher > 0.) {
                // rfull_b = -mu_b = -( r^(t) - psi_r(r^(t))/I_r ), and in the
                // propagator's standardized units
                //     psi_r/I_r = sigma * (1/I_z) * psi_z(r/sigma)
                // so           rfull_b = -r^(t) + sigma (1/I_z) psi_z(r/sigma).
                //
                // GAUSSIAN LIMIT: psi_z(z) = z/(1/I_z) gives
                // sigma (1/I_z) (r/sigma)/(1/I_z) = r, hence rfull_b = 0
                // IDENTICALLY, at every iteration and every r -- not merely at
                // r = 0. That is the check that this matches the derivation,
                // and CVH_CGF_QOP_GAUSSPSI=1 runs it.
                const double z = rprev / csig;
                double psiz;
                if (cgfGaussPsi) {
                  psiz = z / cgfb.invFisher;
                } else {
                  bool clamped = false;
                  psiz = cvhcgf::scoreAt(cgfb, z, &clamped);
                  if (clamped) {
                    ++nCgfClamp_;
                  }
                }
                // Written as the NON-GAUSSIAN part of the score, not as
                // `-r + sigma (1/I) psi`. The two are algebraically identical,
                // but the latter is a catastrophic cancellation: the Gaussian
                // piece of psi reconstructs r, so the answer is the difference
                // of two O(r) numbers and comes out at r * O(eps) instead of
                // exactly zero. Measured: that residual is 1e-17 of r -- and
                // the FIT AMPLIFIES IT TO 2.7e-5 on the fitted q/p, so the
                // difference is not cosmetic (section 25.2).
                //
                // With d = psi_z(z) - z/(1/I_z) the Gaussian limit is
                // bit-exactly zero, and the t = 0 case reduces to
                // sigma (1/I_z) psi_z(0) as before.
                const double dscore = psiz - z / cgfb.invFisher;
                // SHRINKAGE. The re-centring is measured to be anti-correlated
                // with the truth residual and to add scatter of a comparable
                // size, so its variance-optimal scale is not 1 (NOTES_CGFFIT
                // s85). `CgfRecentreDamping` is that scale; 1.0 reproduces
                // every earlier stage bit for bit.
                dx0[0] = cvhcgf::cgfRecentreDamping() * csig * cgfb.invFisher * dscore;
                ++nCgfRecentre_;
                if (getenv("CVH_CGF_QOP_DEBUG") != nullptr) {
                  static int ndbg = 0;
                  if (ndbg < 40 && rprev != 0.) {
                    ++ndbg;
                    std::cout << "### CGFDBG iiter=" << iiter << " iblk=" << iblk
                              << " icons=" << icons
                              << " rprev=" << rprev << " csig=" << csig
                              << " invI=" << cgfb.invFisher << " psiz=" << psiz
                              << " dx0out=" << dx0[0] << std::endl;
                  }
                }
              }
            }
          }
          // ------------------------------------------------------------------

          rfull.segment<nlocalcons>(icons) = dx0;

          // Build the 5 x nlocalparms field+eloss Jacobian: per-mode columns
          // sum the dBx, dBy, dBz transport-Jacobian columns scaled by each
          // mode's Bx/By/Bz basis values at the propagation start; the last
          // column is the unchanged d/dxi column FdFm.col(8).
          Matrix<double, 5, Dynamic> dStateDparams(5, nlocalparms);
          if (perStepFieldModes_) {
            // per-step columns from the propagator, already leg-integrated
            // with the basis sampled at every step midpoint
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
            // this leg did not cross); ihit == 0 reference-state Jacobian
            // correction does not apply to the perturbation columns, but the
            // columns were accumulated on the curvilinear trajectory like
            // col 8, so they get the same treatment as the legacy column
            // (none beyond FdFm assembly).
            dStateDparams.rightCols(nlocaleloss).setZero();
            for (auto const &gc : groupJacs_) {
              dStateDparams.col(nlocalbfield + gc.first) = gc.second;
            }
          } else {
            dStateDparams.col(nlocalbfield) = FdFm.col(8);
          }

          // ----- Numerical-FD closure (debug only) ----------------
          // Validates the analytic Bx/By/Bz chain rule by perturbing dB at the
          // propagation start by epsilon * (dBxPerMode[i], dByPerMode[i],
          // dBzPerMode[i]) for the first few modes, re-propagating, and
          // finite-differencing the 5-component endpoint state. Compares
          // against dStateDparams.col(imode). Runs once per job at the first
          // chain-rule site that has nlocalbfield > 0.
          // Per-step variant: perturb the provider's applied field by
          // eps * basis_i and re-propagate with the provider active --
          // the FD then probes the same per-step application the analytic
          // columns describe.
          // Compares the basis-invariant curvilinear components (qop,
          // lambda, phi), per-mode eps scaled to a target |dB| at the leg
          // start, top modes by basis amplitude -- mirroring the per-leg
          // closure in the two-track maker.
          if (runFDClosure_ && !didFDClosure_ && perStepFieldModes_ && nlocalbfield > 0) {
            didFDClosure_ = true;
            auto qopLamPhi = [](const Eigen::Matrix<double, 7, 1> &s) {
              const double px = s(3), py = s(4), pz = s(5), q = s(6);
              const double pT = std::sqrt(px * px + py * py);
              const double pmag = std::sqrt(pT * pT + pz * pz);
              return Eigen::Vector3d(q / pmag, std::atan2(pz, pT), std::atan2(py, px));
            };
            const Eigen::Vector3d cNom = qopLamPhi(updtsos);
            const double dBtarget = epsilonFDClosure_;
            const double *bxs = nullptr;
            const double *bys = nullptr;
            const double *bzs = nullptr;
            double bdummy[3];
            fieldModeProvider_->sampleAt(propInputState[0], propInputState[1], propInputState[2],
                                         bdummy, bxs, bys, bzs);
            std::vector<std::pair<double, unsigned int>> sorted;
            sorted.reserve(nlocalbfield);
            for (unsigned int i = 0; i < nlocalbfield; ++i) {
              sorted.emplace_back(
                  std::sqrt(bxs[i] * bxs[i] + bys[i] * bys[i] + bzs[i] * bzs[i]), i);
            }
            std::sort(sorted.begin(), sorted.end(),
                      std::greater<std::pair<double, unsigned int>>());
            const unsigned int nTest = std::min<unsigned int>(10u, nlocalbfield);
            std::cout << "===== Per-step field FD closure =====  nModes=" << nlocalbfield
                      << "  testing top-" << nTest << " by basis amplitude  dB_target="
                      << dBtarget << " T  comparing (qop, lambda, phi)" << std::endl;
            std::cout << std::scientific << std::setprecision(4);
            double worstRel = 0.0;
            for (unsigned int j = 0; j < nTest; ++j) {
              const double basisAmp = sorted[j].first;
              const unsigned int imode = sorted[j].second;
              if (basisAmp < 1e-15) {
                continue;
              }
              const double eps = dBtarget / basisAmp;
              const Eigen::Vector3d anac = modeJacs_[imode].head<3>();
              fieldModeProvider_->setInjection(imode, eps);
              auto const &pertResult = g4prop->propagateGenericWithJacobianAltD(
                  propInputState, surface, dB, dxival, dmsval, dionival, -1., g4PartName,
                  matModel_.get(), nullptr, nullptr, fieldModeProvider_.get(), nullptr);
              fieldModeProvider_->setInjection(-1, 0.);
              if (!std::get<0>(pertResult)) {
                std::cout << "  mode " << imode << ": perturbed propagation failed" << std::endl;
                continue;
              }
              const Eigen::Vector3d dCFD = (qopLamPhi(std::get<1>(pertResult)) - cNom) / eps;
              Eigen::Vector3d rel;
              for (int k = 0; k < 3; ++k) {
                rel(k) = std::abs(dCFD(k) - anac(k)) / std::max(std::abs(anac(k)), 1e-30);
                if (rel(k) > worstRel) {
                  worstRel = rel(k);
                }
              }
              std::cout << "  mode " << imode << " (basis=" << basisAmp
                        << ")  rel(qop,lam,phi) = " << rel.transpose() << std::endl;
            }
            std::cout << "===== per-step FD closure: worst rel = " << worstRel << " ====="
                      << std::endl;
            std::cout.unsetf(std::ios_base::floatfield);
          }

          if (runFDClosure_ && !didFDClosure_ && !perStepFieldModes_ && nlocalbfield > 0) {
            const Matrix<double, 5, 1> stateNom =
                Eigen::Matrix<double, 5, 1>(updtsos.head<5>());
            const double eps = epsilonFDClosure_;
            const unsigned int nTest = std::min<unsigned int>(10u, nlocalbfield);
            std::cout << "===== Numerical-FD closure ====="
                      << "  nFieldModes=" << nlocalbfield
                      << "  testing " << nTest << " modes"
                      << "  eps=" << eps << std::endl;
            std::cout << std::scientific << std::setprecision(4);
            double worstRel = 0.0;
            for (unsigned int imode = 0; imode < nTest; ++imode) {
              const Eigen::Vector3d dBpert(
                  dB(0) + eps * dBxPerMode[imode],
                  dB(1) + eps * dByPerMode[imode],
                  dB(2) + eps * dBzPerMode[imode]);
              auto pertResult = g4prop->propagateGenericWithJacobianAltD(
                  propInputState, surface, dBpert, dxival,
                  dmsval, dionival, -1., g4PartName);
              if (!std::get<0>(pertResult)) {
                std::cout << "  mode " << imode
                          << ": perturbed propagation failed" << std::endl;
                continue;
              }
              const Matrix<double, 5, 1> statePert =
                  Eigen::Matrix<double, 5, 1>(std::get<1>(pertResult).head<5>());
              const Matrix<double, 5, 1> dStateFD = (statePert - stateNom) / eps;
              const Matrix<double, 5, 1> dStateAn = dStateDparams.col(imode);
              Matrix<double, 5, 1> rel;
              for (int k = 0; k < 5; ++k) {
                const double scale = std::max(std::abs(dStateAn(k)), 1e-30);
                rel(k) = std::abs(dStateFD(k) - dStateAn(k)) / scale;
                if (rel(k) > worstRel) worstRel = rel(k);
              }
              std::cout << "  mode " << imode
                        << "  max|FD-an|/|an| = " << rel.maxCoeff()
                        << "  FD=[" << dStateFD.transpose() << "]"
                        << "  an=[" << dStateAn.transpose() << "]"
                        << std::endl;
            }
            std::cout << "===== FD closure: worst rel = " << worstRel
                      << " over " << nTest << " modes ====="
                      << std::endl;
            std::cout.unsetf(std::ios_base::floatfield);
            didFDClosure_ = true;
          }
          // ------------------------------------------------------------------

          if (dolocalupdate) {
            Ffull.block<nlocalcons, nlocalstateparms>(icons, fullstateidx) = -Hm*FdFm.leftCols<nlocalstateparms>();
            Ffull.block<nlocalcons, nlocalstateparms>(icons, fullstateidx + nlocalstateparms) = Hp;
            Jfull.block(icons, fullparmidx, nlocalcons, nlocalparms) = -Hm * dStateDparams;
          }
          else {
            Ffull.block<nlocalcons, nlocalstateparms>(icons, fullstateidx) = -FdFm.leftCols<nlocalstateparms>();
            Ffull.block<nlocalcons, nlocalstateparms>(icons, fullstateidx + nlocalstateparms) = Matrix<double, nlocalcons, nlocalstateparms>::Identity();
            Jfull.block(icons, fullparmidx, nlocalcons, nlocalparms) = -dStateDparams;
          }

          Vinvfull.block<nlocalcons, nlocalcons>(icons, icons) = Qinv;
// Vinvfullalt.block<nlocalcons, nlocalcons>(icons, icons) = Qinv;
          Vinvfullalt.block<nlocalcons, nlocalcons>(icons, icons) = (Q - 0.1*dQMS).inverse();

          if (dores) {
            std::vector<Triplet<double>> coeffs;
            for (unsigned int irow = 0; irow < nlocalcons; ++irow) {
              for (unsigned int icol = 0; icol < nlocalcons; ++icol) {
                coeffs.emplace_back(icons + irow, icons + icol, dQMS(irow, icol));
              }
            }
            SparseMatrix<double> &dV = dVs.emplace_back(ncons, ncons);
            dV.setFromTriplets(coeffs.begin(), coeffs.end());
            // MS resolution parameter slot sits right after the bfield block
            // and the eloss slot(s).
            residxs.push_back(iparm + nlocalbfield + nlocaleloss);
            resblockrng.push_back({{icons, nlocalcons}});
            resglobidx.push_back(msglobalidx);
            resvalidhit_.push_back(-1);          // material block, not a hit
            resfamily_.push_back(10);            // multiple scattering

            // Phase B export: Moliere raw step data of the same leg (log
            // sync argument as for the Urban export below).
            for (auto const &ms : g4prop->msStepLog()) {
              msmoliidx.push_back(msglobalidx);
              msmoliv.push_back(ms.effZ);
              msmoliv.push_back(ms.effA);
              msmoliv.push_back(ms.xg);
              msmoliv.push_back(ms.pGeV);
              msmoliv.push_back(ms.beta);
              msmoliv.push_back(ms.thp2);
              msmoliv.push_back(ms.dOverX0);
              // per-element Moliere sums (2026-08-08): effZ/effA are mass
              // averages and both parameters are non-linear in Z
              msmoliv.push_back(ms.zzp1OverA);
              msmoliv.push_back(ms.lnScreenW);
              msmoliv.push_back(ms.stepGroup);
            }
          }

          if (dores) {
            std::vector<Triplet<double>> coeffs;
            for (unsigned int irow = 0; irow < nlocalcons; ++irow) {
              for (unsigned int icol = 0; icol < nlocalcons; ++icol) {
                coeffs.emplace_back(icons + irow, icons + icol, dQI(irow, icol));
              }
            }
            SparseMatrix<double> &dV = dVs.emplace_back(ncons, ncons);
            dV.setFromTriplets(coeffs.begin(), coeffs.end());
            // Ionization resolution slot is the second-after-eloss; offset is
            // bfield block + eloss slot(s) + msres.
            residxs.push_back(iparm + nlocalbfield + nlocaleloss + 1);
            resblockrng.push_back({{icons, nlocalcons}});
            resglobidx.push_back(ioniglobalidx);
            resvalidhit_.push_back(-1);          // material block, not a hit
            resfamily_.push_back(11);            // ionization

            // Physics-CF export: the propagator's Urban step log corresponds
            // to the leg propagation whose dQI was stored above (the log is
            // cleared at each propagate call; the FD-closure / sim re-runs
            // that would overwrite it are debug-gated off in production).
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
              // regime 2/3 (CVH_IONI_EXACTDELTA): two extra columns AFTER cs,
              // so every existing column index is unchanged and the stride is
              // 11 exactly when the switch is off.
              if (G4UniversalFluctuationForExtrapolator::exactDeltaEnabled()) {
                ioniurbanv.push_back(us.rec.beta2);
                ioniurbanv.push_back(us.rec.etot);
              }
              // material group of the step, ALWAYS last so that neither the
              // base nor the exact-delta column indices move (see
              // UrbanIoniStep::stepGroup)
              ioniurbanv.push_back(us.stepGroup);
            }

            // The scale the CGF substitution applied to THIS leg's ionization
            // block, with the number of step rows just pushed so that pooled
            // blocks (several legs under one global index) can be split back
            // apart exactly. 1.0 under CgfQoPMode=0. See the member docs.
            ioniqscaleidx.push_back(ioniglobalidx);
            ioniqscalev.push_back(g4prop->cgfQScale());
            ioniqscalev.push_back(static_cast<float>(g4prop->ioniStepLog().size()));

            // Radiative steps of the same leg, one row per Geant4 step (the
            // radiative log has an entry for every step, unlike the Urban log)
            // plus the two per-process dN/dv shapes. Same log-sync argument as
            // above: the propagator clears these at each propagate call.
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
              radstepv.push_back(rs.stepGroup);   // column 11, appended 2026-09-06
              for (int iv = 0; iv < RADSTEP_NV; ++iv) {
                radstepspecv.push_back(rs.dNdvBrem[iv]);
              }
              for (int iv = 0; iv < RADSTEP_NV; ++iv) {
                radstepspecv.push_back(rs.dNdvPair[iv]);
              }
            }
          }


          // ---- PARMTYPE-15: THE MATERIAL GROUP'S OWN PROCESS NOISE ------
          //
          // `k_g` scales the step's MEAN loss AND, coherently, its MS
          // covariance and ionization variance (`matStepFact` in the
          // propagator's M1 block).  The mean dependence has always been
          // differentiated -- it is the parmtype-15 column of
          // `transportJacobianBxByBzD`, whose only non-zero row is `dqopdxi`
          // -- but the WIDTH dependence never was.  So the quadratic term
          // measured a group's mean loss only while the mass CF measured its
          // width: two functionals of one parameter, one of them blind, which
          // is exactly the configuration in which a -37 % `tec_services` pull
          // can sit unexplained (NOTES 2026-09-06, `resolution/qmsmodel/`).
          //
          // dV/dk_g is the same object the parmtype-10/11 blocks use, split by
          // group: the propagator accumulates (errMS + errI) per group in
          // `groupQs_`, transported by the same Jacobian as dQ/dQ2 and
          // localized by the same `Hm`, so `sum_g dQ_g == dQMS + dQI` exactly.
          if (dores && exportMaterialNoise_ && globalMaterialModel_) {
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
              coeffs.reserve(nlocalcons * nlocalcons);
              for (unsigned int irow = 0; irow < nlocalcons; ++irow) {
                for (unsigned int icol = 0; icol < nlocalcons; ++icol) {
                  coeffs.emplace_back(icons + irow, icons + icol, dQG(irow, icol));
                }
              }
              SparseMatrix<double> &dV = dVs.emplace_back(ncons, ncons);
              dV.setFromTriplets(coeffs.begin(), coeffs.end());
              // this group's own column of THIS propagation's parameter block
              residxs.push_back(iparm + nlocalbfield + gq.first);
              resblockrng.push_back({{icons, nlocalcons}});
              resglobidx.push_back(matGroupGlobalIdx_[gq.first]);
              resvalidhit_.push_back(-1);
              resfamily_.push_back(15);          // global material group
            }
          }

          icons += nlocalcons;

          // One slot per scalar-potential mode, sharing the same global index
          // across all hits (the idxmap collapses these in the final Jacobian).
          for (unsigned int imode = 0; imode < nlocalbfield; ++imode) {
            globalidxv[iparm++] = fieldCorrection_->basisGlobalIdx(imode);
          }

          if (globalMaterialModel_) {
            // One slot per material group, shared global indices as above.
            for (unsigned int g = 0; g < nMatGroups; ++g) {
              globalidxv[iparm++] = matGroupGlobalIdx_[g];
            }
          } else {
            globalidxv[iparm++] = elossglobalidx;
          }

          if (dores) {
            globalidxv[iparm++] = msglobalidx;
            globalidxv[iparm++] = ioniglobalidx;
          }
        }

        if (hit->isValid()) {

          //apply measurement update if applicable
          LocalTrajectoryParameters locparm(localparmsalignprop[0],
                                            localparmsalignprop[1],
                                            localparmsalignprop[2],
                                            localparmsalignprop[3],
                                            localparmsalignprop[4],
                                            localparmsalignprop[5]);
          const TrajectoryStateOnSurface tsostmp(locparm, *hit->surface(), field);

          auto const& preciseHit = cloner.makeShared(hit, tsostmp);
          if (!preciseHit->isValid()) {
            std::cout << "Abort: Failed updating hit" << std::endl;
            ++fitFailHitUpdate_;
            {
              const DetId fdet = hit->geographicalId();
              std::cout << "CVHDIAG reason=HITUPD ihit=" << ihit << " nhit=" << hits.size()
                        << " detid=" << fdet.rawId() << " subdet=" << fdet.subdetId()
                        << " layer=" << trackerTopology->layer(fdet)
                        << " trkpt=" << trackPt << " trketa=" << trackEta << " trkphi=" << trackPhi
                        << " iiter=" << iiter << std::endl;
            }
            valid = false;
            break;
          }

          const bool align2d = detidparms.count(std::make_pair(1, preciseHit->geographicalId()));
          // const bool align2d = detidparms.count(std::make_pair(1, aligndetid));

          const Matrix<double, 2, 2> &Rglued = rgluemap_.at(preciseHit->geographicalId());
          const GloballyPositioned<double> &surfaceglued = surfacemapD_.at(parmdetid);

          const Matrix<double, 5, 5> curvcov = covfull.block<5, 5>(5*(ihit+1), 5*(ihit+1));

          const Matrix<double, 2, 1> localconv = localPositionConvolutionD(updtsos, curvcov, surface);

          {
            constexpr unsigned int nlocalstate = 2;

            const unsigned int nlocalalignment = align2d ? 6 : 5;
            const unsigned int nlocalparms = nlocalalignment;
            
            const unsigned int fullstateidx = 5*(ihit+1) + 3;
            const unsigned int fullparmidx = iparm;

            const bool ispixel = GeomDetEnumerators::isTrackerPixel(preciseHit->det()->subDetector());
            
            const bool hit1d = preciseHit->dimension() == 1;
            
            const Matrix<double, 2, 2> Hu = Hp.bottomRightCorner<2,2>();

            Matrix<double, 2, 1> dy0;
            Matrix<double, 2, 2> iV;
            // rotation from module to strip coordinates
            Matrix2d R;
            
            const double lxcor = localparmsalignprop[3];
            const double lycor = localparmsalignprop[4];


            const Topology &topology = preciseHit->det()->topology();

            // undo deformation correction
            const LocalPoint lpnull(0., 0.);
            const MeasurementPoint mpnull = topology.measurementPosition(lpnull);
            const Topology::LocalTrackPred pred(tsostmp.localParameters().vector());

            auto const defcorr = topology.localPosition(mpnull, pred) - topology.localPosition(mpnull);

            // Rung-E closure mode: fit the SIMULATED hit positions instead
            // of the reconstructed cluster positions (assigned covariances
            // unchanged, so the fit weighting is identical) -- removes the
            // cluster-position (CPE) layer from the response, leaving only
            // the propagator/fit itself. ONLY the measured coordinates are
            // substituted: x for strips, x and y for pixels. The strip-y is
            // unmeasured and reco puts it at the strip center by convention;
            // substituting the true (cm-scale different) sim-y feeds the
            // wedge-module local-polar variance conversion and de-syncs V
            // from the registered dV blocks (found via the identity guard:
            // 33% violations up to 9e3 before this restriction).
            const bool usesimpos = fitSimHitPositions_ && simhit != nullptr;
            const double hitxreco = preciseHit->localPosition().x() - defcorr.x();
            const double hityreco = preciseHit->localPosition().y() - defcorr.y();
            const double hitx = usesimpos ? simhit->localPosition().x() : hitxreco;
            const double hity = (usesimpos && ispixel && !hit1d)
                                          ? simhit->localPosition().y() : hityreco;

            double dxrecsimval = -99.;
            double dyrecsimval = -99.;

            if (simhit != nullptr) {
              // fit-transmission export: pabs at the first/last matched
              // sim hit = the track's true in-tracker momentum profile
              if (simPabsFirst < 0.) {
                simPabsFirst = simhit->pabs();
              }
              simPabsLast = simhit->pabs();
              dxrecsimval = hitx - simhit->localPosition().x();
              dyrecsimval = hity - simhit->localPosition().y();
            }

            double lyoffset = 0.;
            double hitphival = -99.;
            double localphival = -99.;

            if (false) {
              std::cout << "ihit = " << ihit << " orig x y = " << preciseHit->localPosition().x() << " " << preciseHit->localPosition().y() << " defcor x y = " << defcorr.x() << " " << defcorr.y() << " hit x y = " << hitx << " " << hity << std::endl;
            }

            if (hit1d) {
              const ProxyStripTopology *proxytopology = dynamic_cast<const ProxyStripTopology*>(&(preciseHit->det()->topology()));

              dy0[0] = hitx - lxcor;
              dy0[1] = hity - lycor;

              const double striplength = proxytopology->stripLength();
              const double yerr2 = striplength*striplength/12.;
              
              iV = Matrix<double, 2, 2>::Zero();
              iV(0, 0) = preciseHit->localPositionError().xx();
              //iV(1, 1) = yerr2;
              
              R = Matrix2d::Identity();


// std::cout << "1d hit, original x = " << preciseHit->localPosition().x() << " y = " << preciseHit->localPosition().y() << " corrected x = " << hitx << " y = " << hity << std::endl;
            }
            else {
              // 2d hit
// assert(align2d);
              
              
              if (ispixel) {

                dy0[0] = hitx - lxcor;
                dy0[1] = hity - lycor;

                iV << preciseHit->localPositionError().xx(), preciseHit->localPositionError().xy(),
                      preciseHit->localPositionError().xy(), preciseHit->localPositionError().yy();

                R = Matrix2d::Identity();
              }
              else {
                // transform to polar coordinates to end the madness
                //TODO handle the module deformations consistently here (currently equivalent to dropping/undoing deformation correction)

// std::cout << "wedge\n" << std::endl;

                const ProxyStripTopology *proxytopology = dynamic_cast<const ProxyStripTopology*>(&(preciseHit->det()->topology()));

                const TkRadialStripTopology *radialtopology = dynamic_cast<const TkRadialStripTopology*>(&proxytopology->specificTopology());

                const double rdir = radialtopology->yAxisOrientation();
                const double radius = radialtopology->originToIntersection();

                // Wedge modules measure LOCAL PHI. In sim-position mode the
                // residual must use the FULL sim phi (sim x AND sim y):
                // mixing sim-x with the strip-center reco-y gives phi errors
                // of order x*dy/r -- millimeters of arc (found as 5.8-sigma-
                // wide mass pulls with a -40 MeV mean).
                double phihit = rdir*std::atan2(hitx, rdir*hity + radius);
                double rhohit = std::sqrt(hitx*hitx + std::pow(rdir*hity + radius, 2));
                if (usesimpos) {
                  const double lxs = simhit->localPosition().x();
                  const double lys = simhit->localPosition().y();
                  phihit = rdir*std::atan2(lxs, rdir*lys + radius);
                  rhohit = std::sqrt(lxs*lxs + std::pow(rdir*lys + radius, 2));
                }

                // invert original calculation of covariance matrix to extract
                // variance on polar angle. The inversion MUST use the RECO
                // hit phi: the CPE built xx = phierr2*c2i^2 + tan^2*radsigma
                // around the reco position, so tt > 0 is only guaranteed
                // there. With substituted sim positions (fitSimHitPositions)
                // a sim-based tan(phi) can drive tt negative -> negative
                // variance -> indefinite normal equations (found via the
                // export identity guard + per-iteration dumps).
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

                if (simhit != nullptr) {
                  const double lxsim = simhit->localPosition().x();
                  const double lysim = simhit->localPosition().y();

                  const double phisim = rdir*std::atan2(lxsim, rdir*lysim + radius);
                  const double rhosim = std::sqrt(lxsim*lxsim + std::pow(rdir*lysim + radius, 2));

                  dxrecsimval = phihit - phisim;
                  dyrecsimval = rhohit - rhosim;
                }

                iV = Matrix<double, 2, 2>::Zero();
                iV(0, 0) = phierr2;
                //iV(1, 1) = rhoerr2;

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
            
            // Configurable multiplier on the ASSIGNED hit covariance. Applied
            // here, where iV is finished and before anything reads it, so
            // Vinvfull, the resolution dV blocks and the dxerr export are all
            // consistent with one another. Defaults are 1.0 = no change.
            // This is the dead `scalecov = ispixel ? 0.8 : 1.2` of the
            // Vinvfullalt branch turned into something measurable.
            {
              const double covscale = ispixel ? hitCovScalePixel_ : hitCovScaleStrip_;
              if (covscale != 1.0) {
                iV *= covscale;
              }
            }

            rxfull.row(ivalidhit) = R.row(0).cast<float>();
            ryfull.row(ivalidhit) = R.row(1).cast<float>();
            
            validdxeigjac.block<2,2>(2*ivalidhit, 3*(ihit+1)) = R*Hp.bottomRightCorner<2,2>();
            
            // alignment jacobian
            Matrix<double, 2, 6> Aval = Matrix<double, 2, 6>::Zero();

// const Matrix<double, 6, 1> &localparmsalign = alignGlued_ ? globalToLocal(updtsos, surfaceglued) : localparms;
            
// Matrix<double, 6, 1> localparmsalign = localparms;
            Matrix<double, 6, 1> localparmsalign = localparmsalignprop;
            if (alignGlued_) {
              localparmsalign = globalToLocal(updtsosalign, surfaceglued);
            }

            const double localqopval = localparmsalign[0];
            const double localdxdzval = localparmsalign[1];
            const double localdydzval = localparmsalign[2];
            const double localxval = localparmsalign[3];
            const double localyval = localparmsalign[4];

            const double localxvalorig = localparmsalignprop[3];
            const double localyvalorig = localparmsalignprop[4];
                        
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
            
            Matrix<double, 2, 6> A = R*Aval;
            if (alignGlued_) {
              // glued alignment dofs only for out-of-plance
              A.middleCols<3>(2) = R*Rglued*Aval.middleCols<3>(2);
            }
            
            const unsigned int xresglobalidx = dores ? detidparms.at(std::make_pair(8, hit->geographicalId())) : 0;
            const unsigned int yresglobalidx = dores && ispixel ? detidparms.at(std::make_pair(9, hit->geographicalId())) : 0;

            // apply correction to hit covariance from previous iterations if applicable
            const double xxresfact = dores ? std::exp(corparms_[xresglobalidx]) : 1.;
            const double yyresfact = dores && ispixel ? std::exp(corparms_[yresglobalidx]) : 1.;
            const double xyresfact = std::sqrt(xxresfact*yyresfact);
            
// std::cout << "xxresfact = " << xxresfact << " yyresfact = " << yyresfact << " xyresfact = " << xyresfact << std::endl;

            iV(0, 0) *= xxresfact;
            iV(0, 1) *= xyresfact;
            iV(1, 0) *= xyresfact;
            iV(1, 1) *= yyresfact;
            
            // local x/phi resolution variation
            if (dores) {
              // scale diagonal element by exp(parm), preserve correlations
              const double dVxx = iV(0, 0);
              
              std::vector<Triplet<double>> coeffs;
              coeffs.reserve(3);
              coeffs.emplace_back(icons, icons, dVxx);
              
              // cross terms only relevant/exist for pixel hits
              if (ispixel) {
                const double dVxy = 0.5*iV(0, 1);
                coeffs.emplace_back(icons, icons + 1, dVxy);
                coeffs.emplace_back(icons + 1, icons, dVxy);
              }

              SparseMatrix<double> &dV = dVs.emplace_back(ncons, ncons);
              dV.setFromTriplets(coeffs.begin(), coeffs.end());
              residxs.push_back(iparm + nlocalalignment);
              resblockrng.push_back({{icons, ispixel ? 2u : 1u}});
              resglobidx.push_back(xresglobalidx);
              resvalidhit_.push_back(int(ivalidhit));
              resfamily_.push_back(8);           // hit resolution, local x
            }
            
            // local y resolution variation
            if (dores && ispixel) {
              // scale diagonal element by exp(parm), preserve correlations
              const double dVxy = 0.5*iV(0, 1);
              const double dVyy = iV(1, 1);
              
              std::vector<Triplet<double>> coeffs;
              coeffs.reserve(3);
              coeffs.emplace_back(icons, icons + 1, dVxy);
              coeffs.emplace_back(icons + 1, icons, dVxy);
              coeffs.emplace_back(icons + 1, icons + 1, dVyy);

              SparseMatrix<double> &dV = dVs.emplace_back(ncons, ncons);
              dV.setFromTriplets(coeffs.begin(), coeffs.end());
              residxs.push_back(iparm + nlocalalignment + 1);
              resblockrng.push_back({{icons, 2u}});
              resglobidx.push_back(yresglobalidx);
              resvalidhit_.push_back(int(ivalidhit));
              resfamily_.push_back(9);           // hit resolution, local y
            }

            constexpr std::array<unsigned int, 6> alphaidxs = {{0, 2, 3, 4, 5, 1}};

// const double scalecov = hit1d ? 1.2 : 1.0;
            const double scalecov = ispixel ? 0.8 : 1.2;
            
            if (ispixel) {
              constexpr unsigned int nlocalcons = 2;
              
              rfull.segment<nlocalcons>(icons) = dy0;
              Ffull.block<nlocalcons, nlocalstate>(icons, fullstateidx) = -R*Hu;
              
              for (unsigned int ialign = 0; ialign < nlocalalignment; ++ialign) {
                Jfull.block<nlocalcons, 1>(icons, fullparmidx + ialign) = -A.col(alphaidxs[ialign]);
              }
              
              Vinvfull.block<nlocalcons, nlocalcons>(icons, icons) = iV.inverse();
              Vinvfullalt.block<nlocalcons, nlocalcons>(icons, icons) = (scalecov*iV).inverse();
              
              icons += nlocalcons;
            }
            else {
              constexpr unsigned int nlocalcons = 1;
              
              rfull(icons) = dy0(0);
              Ffull.block<1, nlocalstate>(icons, fullstateidx) = (-R*Hu).row(0);
              
              for (unsigned int ialign = 0; ialign < nlocalalignment; ++ialign) {
                Jfull(icons, fullparmidx + ialign) = -A(0, alphaidxs[ialign]);
              }
              
              Vinvfull(icons, icons) = 1./iV(0, 0);
              Vinvfullalt(icons, icons) = 1./(scalecov*iV)(0, 0);
              
              icons += nlocalcons;
            }
            
            for (unsigned int idim=0; idim<nlocalalignment; ++idim) {
              const unsigned int iidx = alphaidxs[idim];
              const DetId ialigndetid = iidx > 1 && iidx < 5 ? aligndetid : preciseHit->geographicalId();
              const unsigned int xglobalidx = detidparms.at(std::make_pair(iidx, ialigndetid));
              globalidxv[iparm] = xglobalidx;
              iparm++;
              if (alphaidxs[idim]==0) {
                hitidxv.push_back(xglobalidx);
              }
            }
            
            if (dores) {
              globalidxv[iparm] = xresglobalidx;
              ++iparm;
            }
            
            if (dores && ispixel) {
              globalidxv[iparm] = yresglobalidx;
              ++iparm;
            }
            
            localqop_iter[ivalidhit] = localqopval;
            localdxdz_iter[ivalidhit] = localdxdzval;
            localdydz_iter[ivalidhit] = localdydzval;
            localx_iter[ivalidhit] = localxval;
            localy_iter[ivalidhit] = localyval;

            dxrecgen_iter[ivalidhit] = dy0[0];
            dyrecgen_iter[ivalidhit] = dy0[1];

            const double localqopvar = covfull(5*(ihit+1), 5*(ihit+1));
            localqoperralt[ivalidhit] = std::sqrt(localqopvar);
            
            if (iiter == 0) {
              
              // fill hit validation information
// Vector2d dyrecgenlocal;
// dyrecgenlocal << dy0[0].value().value(), dy0[1].value().value();
// const Vector2d dyrecgeneig = R*dyrecgenlocal;
// dxrecgen.push_back(dyrecgeneig[0]);
// dyrecgen.push_back(dyrecgeneig[1]);
              dxrecgen.push_back(dy0[0]);
              dyrecgen.push_back(dy0[1]);
              
              dxerr.push_back(std::sqrt(iV(0,0)));
              dyerr.push_back(std::sqrt(iV(1,1)));

              const Matrix<double, 6, 1> localstatedebug = globalToLocal(updtsos, surface);

              localqop.push_back(localstatedebug[0]);
              localdxdz.push_back(localstatedebug[1]);
              localdydz.push_back(localstatedebug[2]);
              localx.push_back(localstatedebug[3]);
              localy.push_back(localstatedebug[4]);

// localqop.push_back(localqopval);
// localdxdz.push_back(localdxdzval);
// localdydz.push_back(localdydzval);
// localx.push_back(localxval);
// localy.push_back(localyval);
              
              const Matrix<double, 5, 5> Qtotlocal = Hp*Qtot*Hp.transpose();
              
              localqoperr.push_back(std::sqrt(Qtotlocal(0, 0)));
              localdxdzerr.push_back(std::sqrt(Qtotlocal(1, 1)));
              localdydzerr.push_back(std::sqrt(Qtotlocal(2, 2)));
              localxerr.push_back(std::sqrt(Qtotlocal(3, 3)));
              localyerr.push_back(std::sqrt(Qtotlocal(4, 4)));
              
              hitlocalx.push_back(hitx);
              hitlocaly.push_back(hity);
              
              dEpred.push_back(dEpredval);
              sigmadE.push_back(sigmadEval);
              Epred.push_back(epred);
              E.push_back(e);
              
              hitphi.push_back(hitphival);
              localphi.push_back(localphival);

              landauDelta.push_back(deltaTotal);
              landauW.push_back(wTotal);

              const TrackerSingleRecHit* tkhit = dynamic_cast<const TrackerSingleRecHit*>(preciseHit.get());
              assert(tkhit != nullptr);

              hitDetId.push_back(preciseHit->geographicalId().rawId());
              hitThickness.push_back(tkhit->det()->surface().bounds().thickness());
              
              if (ispixel) {
                const SiPixelCluster& cluster = *tkhit->cluster_pixel();
                const SiPixelRecHit *pixhit = dynamic_cast<const SiPixelRecHit*>(tkhit);
                assert(pixhit != nullptr);
  // std::cout << "pixel cluster sizeX = " << cluster.sizeX() <<" sizeY = " << cluster.sizeY() << std::endl;
                clusterSize.push_back(cluster.size());
                clusterSizeX.push_back(cluster.sizeX());
                clusterSizeY.push_back(cluster.sizeY());
                clusterCharge.push_back(cluster.charge());
                clusterChargeBin.push_back(pixhit->qBin());
                clusterOnEdge.push_back(pixhit->isOnEdge());
                if (pixhit->hasFilledProb()) {
                  clusterProbXY.push_back(pixhit->clusterProbability(0));
                }
                else {
                  clusterProbXY.push_back(2.);
                }
                
                clusterSN.push_back(-99.);
                stripsToEdge.push_back(-99);

                // uProj is a STRIP CPE concept; the pixel template is
                // indexed by (angle, qbin) instead, and those are already
                // exported as localdxdz/localdydz and clusterChargeBin.
                hitUProj.push_back(-99.f);
                hitStripRec.push_back(-99.f);
                hitStripSim.push_back(-99.f);
                hitFirstStrip.push_back(-99);
                const PixelTopology *pixtopology =
                    dynamic_cast<const PixelTopology*>(&(tkhit->det()->topology()));
                hitPitch.push_back(pixtopology != nullptr ? pixtopology->pitch().first : -99.f);
              }
              else {
                const StripTopology* striptopology = dynamic_cast<const StripTopology*>(&(tkhit->det()->topology()));
                const SiStripCluster& cluster = *tkhit->cluster_strip();
                siStripClusterInfo_.setCluster(cluster, preciseHit->geographicalId().rawId());
                clusterSize.push_back(cluster.amplitudes().size());
                clusterSizeX.push_back(cluster.amplitudes().size());
                clusterSizeY.push_back(1);
                clusterCharge.push_back(cluster.charge());
                clusterChargeBin.push_back(-99);
                clusterOnEdge.push_back(-99);
                clusterProbXY.push_back(-99.);
                clusterSN.push_back(siStripClusterInfo_.signalOverNoise());


                const uint16_t firstStrip = cluster.firstStrip();
                const uint16_t lastStrip = cluster.firstStrip() + cluster.amplitudes().size() - 1;
                stripsToEdge.push_back(std::min<int>(firstStrip, striptopology->nstrips() - 1 - lastStrip));

                // The CPE's own independent variable. getAlgoParam is given
                // the SAME LocalTrajectoryParameters the cloner passed, so
                // this is the uProj that produced localPositionError() above,
                // not a reconstruction of it.
                float uProjVal = -99.f;
                const StripGeomDetUnit *stripdu =
                    dynamic_cast<const StripGeomDetUnit*>(tkhit->det());
                if (stripCPEForExport != nullptr && stripdu != nullptr) {
                  uProjVal = stripCPEForExport->getAlgoParam(*stripdu, locparm).afullProjection;
                }
                hitUProj.push_back(uProjVal);
                hitPitch.push_back(striptopology->localPitch(preciseHit->localPosition()));

                // Strip COORDINATES, so the true impact point can be expressed
                // in the lattice frame of the cluster that measured it. The
                // eta / S-curve is E[rec - sim | position within the cluster],
                // and referring the phase to the ABSOLUTE lattice instead
                // mixes the two possible strip assignments near a boundary:
                // measured on the muon gun, that compresses the N=1 ramp from
                // the geometric +-0.5 to +-0.29 and its slope from -1 to -0.49.
                hitStripRec.push_back(striptopology->strip(preciseHit->localPosition()));
                hitStripSim.push_back(simhit != nullptr
                    ? striptopology->strip(simhit->localPosition()) : -99.f);
                hitFirstStrip.push_back(int(cluster.firstStrip()));
              }
              
  // if (ispixel) {
  // const SiPixelCluster& cluster = *tkhit->cluster_pixel();
  // dxreccluster.push_back(cluster.x() - preciseHit->localPosition().x());
  // dyreccluster.push_back(cluster.y() - preciseHit->localPosition().y());
  // }
  // else {
  // dxreccluster.push_back(-99.);
  // dyreccluster.push_back(-99.);
  // }
              
              if (doSim_) {
                if (simhit != nullptr) {
                  
// std::cout << "ihit = " << ihit << " eloss = " << simhit->energyLoss() << std::endl;
                  
                  Vector2d dy0simgenlocal;
                  dy0simgenlocal << simhit->localPosition().x() - localxval,
                                  simhit->localPosition().y() - localyval;
                  const Vector2d dysimgeneig = R*dy0simgenlocal;
                  dxsimgen.push_back(dysimgeneig[0]);
                  dysimgen.push_back(dysimgeneig[1]);
  // dxsimgen.push_back(simhit->localPosition().x() - updtsos.localPosition().x());
  // dysimgen.push_back(simhit->localPosition().y() - updtsos.localPosition().y());
                  
                  Vector2d dy0simgenlocalconv;
                  dy0simgenlocalconv << simhit->localPosition().x() - localxval - localconv[0],
                                  simhit->localPosition().y() - localyval - localconv[1];
                  const Vector2d dysimgeneigconv = R*dy0simgenlocalconv;
                  dxsimgenconv.push_back(dysimgeneigconv[0]);
                  dysimgenconv.push_back(dysimgeneigconv[1]);
                  
                  dxsimgenlocal.push_back(dy0simgenlocal[0]);
                  dysimgenlocal.push_back(dy0simgenlocal[1]);
                  
                  // Vector2d dyrecsimlocal;
                  // dyrecsimlocal << preciseHit->localPosition().x() - simhit->localPosition().x(),
                  // preciseHit->localPosition().y() - simhit->localPosition().y();
                  // const Vector2d dyrecsimeig = R*dyrecsimlocal;
                  // dxrecsim.push_back(dyrecsimeig[0]);
                  // dyrecsim.push_back(dyrecsimeig[1]);

                  dxrecsim.push_back(dxrecsimval);
                  dyrecsim.push_back(dyrecsimval);
                  simHitNCand.push_back(simhitNCandVal);


                  
                  dE.push_back(dEval);
                  
                  auto const momsim = simhit->momentumAtEntry();

                  const double eentry = std::sqrt(std::pow(simhit->pabs(), 2) + trackmass*trackmass);
                  const double emid = eentry - 0.5*simhit->energyLoss();
                  const double simqopval = genpart->charge()/std::sqrt(emid*emid - trackmass*trackmass);
// std::cout << "eloss = " << simhit->energyLoss() << std::endl;
                  
                  // "hybrid state" trying to adjust for entry point -> midpoint
                  // simlocalqop.push_back(simqopval);
                  // simlocaldxdz.push_back(momsim.x()/momsim.z());
                  // simlocaldydz.push_back(momsim.y()/momsim.z());
                  // simlocalx.push_back(simhit->localPosition().x());
                  // simlocaly.push_back(simhit->localPosition().y());

                  // just fill entry state for debugging
                  simlocalqop.push_back(genpart->charge()/simhit->pabs());
                  simlocaldxdz.push_back(momsim.x()/momsim.z());
                  simlocaldydz.push_back(momsim.y()/momsim.z());
                  simlocalx.push_back(simhit->entryPoint().x());
                  simlocaly.push_back(simhit->entryPoint().y());
                    
                  if (false) {

                    const Point3DBase<double, LocalTag> simlocalpos = simhit->entryPoint();
                    const Vector3DBase<double, LocalTag> simlocalmom = simhit->momentumAtEntry();

          // std::cout << "simlocalpos" << simlocalpos << std::endl;

                    const Point3DBase<double, GlobalTag> simglobalpos = surface.toGlobal(simlocalpos);
                    const Vector3DBase<double, GlobalTag> simglobalmom = surface.toGlobal(simlocalmom);

                    Matrix<double, 7, 1> simtsos;
                    simtsos[0] = simglobalpos.x();
                    simtsos[1] = simglobalpos.y();
                    simtsos[2] = simglobalpos.z();
                    simtsos[3] = simglobalmom.x();
                    simtsos[4] = simglobalmom.y();
                    simtsos[5] = simglobalmom.z();
                    simtsos[6] = genpart->charge();

                    auto propresultsim = g4prop->propagateGenericWithJacobianAltD(simtsos, surface, dB, dxival, 0., 0., -1., g4PartName);

                    if (std::get<0>(propresultsim)) {
                      simtsos = std::get<1>(propresultsim);

                      const Point3DBase<double, GlobalTag> simglobalposprop(simtsos[0], simtsos[1], simtsos[2]);
                      const Vector3DBase<double, GlobalTag> simglobalmomprop(simtsos[3], simtsos[4], simtsos[5]);

                      const Point3DBase<double, LocalTag> simlocalposprop = surface.toLocal(simglobalposprop);
                      const Vector3DBase<double, LocalTag> simlocalmomprop = surface.toLocal(simglobalmomprop);

                      Matrix<double, 5, 1> simlocalparms;
                      simlocalparms[0] = simtsos[6]/simtsos.segment<3>(3).norm();
                      simlocalparms[1] = simlocalmomprop.x()/simlocalmomprop.z();
                      simlocalparms[2] = simlocalmomprop.y()/simlocalmomprop.z();
                      simlocalparms[3] = simlocalposprop.x();
                      simlocalparms[4] = simlocalposprop.y(); 
                      
                      simlocalqopprop.push_back(simlocalparms[0]);
                      simlocaldxdzprop.push_back(simlocalparms[1]);
                      simlocaldydzprop.push_back(simlocalparms[2]);
                      simlocalxprop.push_back(simlocalparms[3]);
                      simlocalyprop.push_back(simlocalparms[4]);
                    }
                    else {
                      simlocalqopprop.push_back(-99.);
                      simlocaldxdzprop.push_back(-99.);
                      simlocaldydzprop.push_back(-99.);
                      simlocalxprop.push_back(-99.);
                      simlocalyprop.push_back(-99.); 
                    }
                  }
                  else {
                    simlocalqopprop.push_back(-99.);
                    simlocaldxdzprop.push_back(-99.);
                    simlocaldydzprop.push_back(-99.);
                    simlocalxprop.push_back(-99.);
                    simlocalyprop.push_back(-99.);
                  }
                  
                }
                else {
                  dxsimgen.push_back(-99.);
                  dysimgen.push_back(-99.);
                  dxsimgenconv.push_back(-99.);
                  dysimgenconv.push_back(-99.);
                  dxsimgenlocal.push_back(-99.);
                  dysimgenlocal.push_back(-99.);
                  dxrecsim.push_back(-99.);
                  dyrecsim.push_back(-99.);
                  simHitNCand.push_back(simhitNCandVal);
                  dE.push_back(-99.);
                  
                  simlocalqop.push_back(-99.);
                  simlocaldxdz.push_back(-99.);
                  simlocaldydz.push_back(-99.);
                  simlocalx.push_back(-99.);
                  simlocaly.push_back(-99.);
                  
                  simlocalqopprop.push_back(-99.);
                  simlocaldxdzprop.push_back(-99.);
                  simlocaldydzprop.push_back(-99.);
                  simlocalxprop.push_back(-99.);
                  simlocalyprop.push_back(-99.);
                }
              }

            }
            
          };
          
          ivalidhit++;
            
        }
        
      }
            
      if (!valid) {
        break;
      }

      
      assert(iparm == npars);
      
// Fsparse = Ffull.rightCols(nstatefree).sparseView();
      Fsparse = Ffull(Eigen::placeholders::all, freestateidxs).sparseView();
      Vinvsparse = Vinvfull.sparseView();

      // ---- chi2-based (Armijo) retroactive backtracking -------------------
      // rfull has just been assembled at the CURRENT linearization point, i.e.
      // this is the REALIZED chi2 of the step taken at the end of the previous
      // iteration. If it fails the sufficient-decrease test, restore the
      // previous linearization, halve that step and redo the iteration. This
      // costs no extra propagation on the accept path (the chi2 is assembled
      // anyway) and one repeated iteration per halving on the reject path.
      // The slack term absorbs the chi2 wobble from relinearization (the
      // propagation/material model is re-evaluated at the new state), so only
      // a genuine chi2 blow-up triggers a halving.
      const double chisq0valNow = rfull.dot(Vinvsparse * rfull);
      if (stepBacktracking_ && iiter >= stepBacktrackFromIter_ && std::isfinite(chisq0valPrev) &&
          nChi2Bt < maxChi2Backtrack_ && nChi2BtTotal < maxChi2Backtrack_ * niters) {
        const double thresh = chisq0valPrev + armijoC_ * predDecrPrev +
                              armijoSlack_ * std::max(1., std::abs(chisq0valPrev));
        if (!(chisq0valNow <= thresh)) {
          refFts = refFtsSnap;
          layerStates = layerStatesSnap;
          dxfull *= 0.5;
          predDecrPrev *= 0.5;  // conservative: the model decrease shrinks with t
          ++nChi2Bt;
          ++nChi2BtTotal;
          ++stepBacktrackEvents_;
          if (!stepBacktrackedThisFit) {
            stepBacktrackedThisFit = true;
            ++fitStepBacktracked_;
          }
          if (stepPrints_ < stepPrintLimit_) {
            ++stepPrints_;
            std::cout << "GN step backtracked: iiter = " << iiter
                      << " chisq " << chisq0valPrev << " -> " << chisq0valNow
                      << " (thresh " << thresh << ")"
                      << " nbt = " << nChi2Bt
                      << " seed(q,pt,eta)=(" << track.charge() << "," << track.pt() << "," << track.eta() << ")"
                      << std::endl;
          }
          iiter -= 1;  // loop ++ redoes the same iteration from the snapshot
          continue;
        }
      }
      // step accepted
      nChi2Bt = 0;
      chisq0valPrev = chisq0valNow;

      VinvF = Vinvsparse*Fsparse;


      
      Cinvd.compute(Fsparse.transpose()*VinvF);


      // randomize residuals given covariance matrix
      if (false) {
        const SparseMatrix<double> Vinvsparsealt = Vinvfullalt.sparseView();
        const SparseMatrix<double> VinvFalt = Vinvsparsealt*Fsparse;

        const SimplicialLDLT<SparseMatrix<double>> Cinvdalt(Fsparse.transpose()*VinvFalt);
        
        const SparseMatrix<double> FtVinvalt = VinvFalt.transpose();
        const SparseMatrix<double> Rsparsealt = Vinvsparsealt - VinvFalt*Cinvdalt.solve(FtVinvalt);
        const MatrixXd Ralt = Rsparsealt;

        const Eigen::SelfAdjointEigenSolver<MatrixXd> Reig(Ralt);

        rfull = VectorXd::Zero(ncons);
        for (int ieig = ncons - ndof; ieig < ncons; ++ieig) {
          const double rnd = gRandom->Gaus();
          rfull += (rnd/std::sqrt(Reig.eigenvalues()(ieig)))*Reig.eigenvectors().col(ieig);
        }

      }
      
// SimplicialLLT<SparseMatrix<double>> lltsparse;
// lltsparse.compute(Fsparse.transpose()*VinvF);
// 
// const SparseMatrix<double> lfact = lltsparse.matrixL();
// std::cout << "lfact size = " << lfact.rows()*lfact.cols() << " nonzeros = " << lfact.nonZeros() << std::endl;
      
// const SparseMatrix<double> testsparse = Fsparse.transpose()*VinvF;
// const MatrixXd testdense = testsparse;
// SimplicialLDLT<SparseMatrix<double>> ldltsparse;
// LDLT<MatrixXd> ldltdense(testdense.rows());
// 
// constexpr unsigned int ntest = 10000;
// 
// double sparsetime;
// double densetime;
// 
// {
// auto start = std::chrono::high_resolution_clock::now();
// for (unsigned int itest = 0; itest < ntest; ++itest) {
// ldltsparse.compute(testsparse);
// }
// auto stop = std::chrono::high_resolution_clock::now();
// 
// auto duration = stop - start;
// sparsetime = duration.count();
// std::cout << "sparse decomp: " << sparsetime << std::endl;
// }
// 
// {
// auto start = std::chrono::high_resolution_clock::now();
// for (unsigned int itest = 0; itest < ntest; ++itest) {
// ldltdense.compute(testdense);
// }
// auto stop = std::chrono::high_resolution_clock::now();
// 
// auto duration = stop - start;
// densetime = duration.count();
// std::cout << "dense decomp: " << densetime << std::endl;
// }
// 
// const double timeratio = sparsetime/densetime;
// std::cout << "size = " << testdense.rows() << " fillfactor = " << double(testsparse.nonZeros())/double(testsparse.rows()*testsparse.cols()) << " sparse/dense = " << timeratio << std::endl;
      
      dxfree = -Cinvd.solve(VinvF.transpose()*rfull);

      // Fail fast on a non-finite update (drained/runaway propagation leg
      // or singular normal matrix) instead of leaking NaN into the state.
      if (!dxfree.allFinite()) {
        std::cout << "Abort: non-finite parameter update from solve!"
                  << " iiter = " << iiter
                  << " seed(q,pt,eta)=(" << track.charge() << "," << track.pt() << "," << track.eta() << ")"
                  << std::endl;
        ++fitFailNaN_;
        valid = false;
        break;
      }

      dxfull = VectorXd::Zero(nstateparms);
      dxfull(freestateidxs) = dxfree;

      // Momentum-floor safeguard on the Gauss-Newton step. Large relative
      // q/p steps are legitimate (forward tracks routinely repull q/p by
      // O(100%) in the first iterations), so a generic step cap would bind
      // on healthy fits. Two update outcomes are guarded, by scaling the
      // whole step vector (direction preserved):
      //  - the updated momentum drops below clampMomentumFloor_ (which must
      //    stay safely above the propagator's PropagationPtotLimit refusal
      //    that would otherwise kill such fits) -- ALWAYS clamped;
      //  - q/p changes sign (charge flip). Crossing q/p = 0 is p -> inf, so a
      //    flip is only dangerous for a stiff track (where it signals a
      //    diverging step). For a high-p track (p > allowChargeFlipAboveP_)
      //    the seed charge is genuinely ambiguous and the flip is a legitimate
      //    correction, so it is PERMITTED (subject to the momentum floor on the
      //    far side). Default allowChargeFlipAboveP_ = 1e9 GeV -> never allowed.
      {
        const double pref = refFts.segment<3>(3).norm();
        const double qopref = pref > 0. ? refFts[6] / pref : 0.;
        const double dqop = dxfull[0];
        const double qopupd = qopref + dqop;
        const double pFloor = clampMomentumFloor_;  // GeV
        // |q/p| threshold below which a sign flip is allowed (high p, stiff track)
        const double qopFlipAllow = 1. / std::max(allowChargeFlipAboveP_, pFloor);
        double stepscale = 1.;
        if (maxMomentumStepFactor_ > 1.) {
          // NEW (2026-09-05): relative trust region in q/p. The lower bound
          // max(pFloor, p_ref/f) is always strictly below p_ref, so a track
          // whose true momentum is below pFloor is no longer pinned there and
          // the scale can no longer be zero. The upward cap p_ref*f guards the
          // opposite runaway (p -> inf / sign flip of a stiff track) and, at
          // f = 2, reproduces the legacy half-way-to-zero flip rule exactly.
          bool flipProtect = false;
          stepscale = cvhstep::legStepScaleRel(qopref, dqop, pFloor,
                                               maxMomentumStepFactor_,
                                               qopFlipAllow, &flipProtect);
          if (flipProtect) {
            ++nChargeFlipProtect;
          }
        } else if (qopref != 0. && dqop != 0.) {
          // LEGACY absolute-floor-only clamp (maxMomentumStepFactor <= 1)
          if (std::abs(qopupd) > 1. / pFloor) {
            // far-side momentum below the floor: land on p = pFloor, same
            // charge (catches diverging steps, incl. sign flips toward low p)
            stepscale = (std::copysign(1. / pFloor, qopref) - qopref) / dqop;
          } else if (qopupd * qopref <= 0. && std::abs(qopref) > qopFlipAllow) {
            // sign flip of a track too stiff for the flip to be ambiguous:
            // treat as divergence, stop half-way toward q/p = 0. This track
            // "wanted" the opposite charge -> flag for the two-hypothesis fit.
            stepscale = -0.5 * qopref / dqop;
            ++nChargeFlipProtect;
          }
          // else: benign same-sign step, or a high-p (ambiguous) charge flip
          // with far-side p >= pFloor -> allowed unchanged
        }
        if (qopref != 0. && dqop != 0. && qopupd * qopref <= 0. && stepscale >= 1.) {
          ++fitChargeFlipAllowed_;
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
            std::cout << "GN step clamped: iiter = " << iiter
                      << " qopref = " << qopref << " dqop = " << dqop
                      << " scale = " << stepscale
                      << " seed(q,pt,eta)=(" << track.charge() << "," << track.pt() << "," << track.eta() << ")"
                      << std::endl;
          }
        }
      }

      // Limit-cycle damping (see ctor comment). Applied after the momentum
      // floor and before the EDM bookkeeping so the recorded edmval refers
      // to the step actually taken.
      if (gnDampAfter_ > 0 && iiter >= gnDampAfter_) {
        dxfree *= gnDampFactor_;
        dxfull *= gnDampFactor_;
        stepScaleApplied *= gnDampFactor_;
      }

      // CGF IRLS: the realised block residual of the step just taken, for the
      // next iteration's re-centring. Taken AFTER the momentum-floor clamp and
      // the damping, so it refers to the step actually applied. One dot
      // product per block; only the q/p row (icons + 0) is needed.
      if (!cgfBlkRow.empty()) {
        cgfRprev.assign(cgfBlkRow.size(), 0.);
        for (unsigned int ib = 0; ib < cgfBlkRow.size(); ++ib) {
          cgfRprev[ib] = Ffull.row(cgfBlkRow[ib]).dot(dxfull);
        }
      }

      const double deltachisq = rfull.transpose()*VinvF*dxfree;
      edmval = -deltachisq;

      // Quadratic-model chi2 change of the step ACTUALLY applied (dxfree has
      // already been scaled by stepScaleApplied above, so deltachisq = t*d with
      // d the full-step value; the model change is d*(2t - t^2) =
      // deltachisq*(2 - t)). Consumed by the Armijo test at the top of the
      // NEXT iteration, where the realized chi2 becomes available.
      predDecrPrev = deltachisq * (2. - stepScaleApplied);

      // Realized chi2 tracking (this maker historically never filled the
      // chisqval/deltachisqval members in-loop -- they were stale storage;
      // the always-stored chisqval branch is still overwritten with the
      // final residual-projector value after the loop as before). chi2 at
      // the current linearization point plus the predicted change of this
      // step, mirroring the two-track maker; deltachisqval is the realized
      // iteration-to-iteration change (per-iteration debug dump).
      {
        const double chisq0val = chisq0valNow;
        const double chisqcur = chisq0val + deltachisq;
        deltachisqval = chisqcur - chisqvalold;
        chisqvalold = chisqcur;
        chisqval = chisqcur;
      }

      const Vector5d dxref = dxfull.head<5>();

// std::cout << "iiter = " << iiter << std::endl;
// std::cout << "dxfree.head<5>():\n" << dxfree.head<5>() << std::endl;
// std::cout << "dxfull.head<5>():\n" << dxfull.head<5>() << std::endl;
// std::cout << "dxref:\n" <<dxref << std::endl;

      covfull = MatrixXd::Zero(nstateparms, nstateparms);

      //TODO temporary needed to work around Eigen fancy indexing limitation (maybe fixed in newer Eigen versions?)
      covfull(freestateidxs, freestateidxs) = Cinvd.solve(MatrixXd::Identity(nstatefree, nstatefree)).eval();

      const Matrix<double, 5, 5> Cinner = covfull.topLeftCorner<5,5>();
      
      if (dogen) {
        edmvalref = 0.;
      }
      else {
        const Matrix<double, 5, 5> hessref = Cinner.inverse();
        const double deltachisqref = -0.5*dxref.transpose()*hessref*dxref;
        edmvalref = -deltachisqref;
      }

      
      //fill output with corrected state and covariance at reference point
      refParms.fill(0.);
      refCov.fill(0.);
      
      const double qbp = refFts[6]/refFts.segment<3>(3).norm();
      const double lam = std::atan(refFts[5]/std::sqrt(refFts[3]*refFts[3] + refFts[4]*refFts[4]));
      const double phi = std::atan2(refFts[4], refFts[3]);
      
      const double qbpupd = qbp + dxref[0];
      const double lamupd = lam + dxref[1];
      const double phiupd = phi + dxref[2];

      // [openspec §9.4.a] One-shot kaon-q/p diagnostic, exit side. Compares
      // the post-iteration refit q/p (qbpupd) to the input track q/p, plus
      // the iteration-end |p|.
      if (dumpKaonQp) {
        const double refP = refFts.segment<3>(3).norm();
        const double inP = std::sqrt(track.momentum().mag2());
        std::cout << "===== [§9.4.a kaon q/p EXIT] track itrack=" << itrack << " =====\n"
                  << "  refFts(end) : charge[6]=" << refFts[6]
                  << "  |p|=" << refP
                  << "  qbp_pre=" << qbp << "  dxref[0]=" << dxref[0]
                  << "  qbpupd="  << qbpupd << "\n"
                  << "  input     : q/p=" << track.charge()/std::max(inP, 1e-12)
                  << "  |p|=" << inP << "\n"
                  << "  ratio refit/input  q/p=" << qbpupd/(track.charge()/std::max(inP, 1e-12))
                  << "  |p|=" << refP/std::max(inP, 1e-12)
                  << "  lam ref/in=" << lamupd << "/" << std::atan(track.momentum().z()/track.pt())
                  << std::endl;
      }

      refParms[0] = qbpupd;
      refParms[1] = lamupd;
      refParms[2] = phiupd;
      // Position (impact) parameters at the PCA to the beamspot: d0 (=dxy) and
      // z0 (=dz). Reference value from the converged state (cart2pca of refFts)
      // plus the reference-block update, mirroring the momentum fill above.
      // Resolves the longstanding "fill position parameters" TODO.
      const Matrix<double, 5, 1> statepcaRef = cart2pca(refFts, *bsH);
      refParms[3] = statepcaRef[3] + dxref[3];
      refParms[4] = statepcaRef[4] + dxref[4];

      // PHI IS AN ANGLE AND THE RESIDUAL IS NOT.  `genParms[2]` comes out of
      // `cart2pca` in (-pi, pi]; `refParms[2] = phi + dxref[2]` is the fitted
      // value and is NOT wrapped, so a track sitting on the branch cut has
      // `refParms[2] - genParms[2] = +-2pi`, i.e. ~4e4 sigma. Measured 8
      // tracks in 20 000 (0.04 %), enough to make `Var(z_phi) = 2.9e4`
      // (NOTES 2026-09-09).  Put the gen value on the same branch as the
      // fitted one HERE, once, rather than leaving every reader to remember.
      if (genpart != nullptr) {
        const double dphi = double(refParms[2]) - double(genParms[2]);
        if (std::abs(dphi) > M_PI) {
          genParms[2] = float(double(genParms[2]) + 2. * M_PI * std::round(dphi / (2. * M_PI)));
        }
      }

      refParmsMomD[0] = qbpupd;
      refParmsMomD[1] = lamupd;
      refParmsMomD[2] = phiupd;
      
      Map<Matrix<float, 5, 5, RowMajor> >(refCov.data()).triangularView<Upper>() = (Cinner).cast<float>().triangularView<Upper>();
      
      if (iiter==0) {
        refParms_iter0 = refParms;
        refCov_iter0 = refCov;
      }


        
      
      niter = iiter + 1;
      // the CONVERGED iteration's reference energy loss is the exported one,
      // for the same reason `niter` is written here
      dEref = float(dErefIter);
      maxfracloss = float(maxFracLossIter);

      // Per-iteration trajectory record (pushed before the convergence
      // break so the break-triggering iteration is included).
      if (debugPerIterDump_) {
        chisqval_iter.push_back(static_cast<double>(chisqval));
        edmval_iter.push_back(static_cast<double>(edmval));
        edmvalref_iter.push_back(static_cast<double>(edmvalref));
        deltachisqval_iter.push_back(static_cast<double>(deltachisqval));
      }

      if (std::isnan(edmval) || std::isinf(edmval)) {
        std::cout << "WARNING: invalid parameter update!!!" << " edmval = " << edmval << " lamupd = " << lamupd << " deltachisqval = " << deltachisqval << std::endl;
        ++fitFailNaN_;
        std::cout << "CVHDIAG reason=PARUPD ihit=-1 nhit=" << nValidHits
                  << " detid=0 subdet=-1 layer=-1"
                  << " trkpt=" << trackPt << " trketa=" << trackEta << " trkphi=" << trackPhi
                  << " iiter=" << iiter << " edmval=" << edmval << " lamupd=" << lamupd << std::endl;
        valid = false;
        break;
      }
      
      
// std::cout << "iiter = " << iiter << " edmval = " << edmval << " edmvalref " << edmvalref << " deltachisqval = " << deltachisqval << " chisqval = " << chisqval << std::endl;

      if (anomDebug) {
        std::cout << "anomDebug: iiter = " << iiter << " edmval = " << edmval << " deltachisqval = " << deltachisqval << " chisqval = " << chisqval << std::endl;
      }

      
      if (iiter > 0 && dolocalupdate && edmval < edmConvergence_) {
        break;
      }
      else if (iiter > 0 && !dolocalupdate && edmvalref < edmConvergence_) {
        break;
      }
    
    }
    
    if (!valid) {
      // this hypothesis failed to converge; caller decides what to keep
      return {false, std::numeric_limits<double>::max()};
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
    
    const unsigned int nparsfinal = globalidxvfinal.size();
    
    MatrixXd Jfinal = MatrixXd::Zero(ncons, nparsfinal);
    
    for (unsigned int i = 0; i < npars; ++i) {
      const unsigned int iidx = idxmap.at(globalidxv[i]);
      Jfinal.col(iidx) += Jfull.col(i);
    }
    
    const SparseMatrix<double> Jsparse = Jfinal.sparseView();
        
    std::vector<unsigned int> residxsfinal;
    residxsfinal.reserve(residxs.size());
    for (auto idx : residxs) {
      const unsigned int iidx = idxmap.at(globalidxv[idx]);
      residxsfinal.push_back(iidx);
    }
    
// dxdparms = -Cinvd.solve(d2chisqdxdparmsfinal).transpose();
    
// grad = dchisqdparmsfinal + dxdparms*dchisqdx; 
// grad = dchisqdparmsfinal + d2chisqdxdparmsfinal.transpose()*dxfull;
// hess = d2chisqdparms2final + dxdparms*d2chisqdxdparmsfinal;
    
    nParms = nparsfinal;

    //TODO avoid explicitly storing the transpose?
    const SparseMatrix<double> FtVinv = VinvF.transpose();
// const MatrixXd R = Vinvsparse - VinvF*Cinvd.solve(FtVinv);
    const SparseMatrix<double> Rsparse = Vinvsparse - VinvF*Cinvd.solve(FtVinv);
    const MatrixXd R = Rsparse;
    
    
// constexpr unsigned int ntest = 1000;
// 
// auto time0 = std::chrono::high_resolution_clock::now();
// for (unsigned int itest = 0; itest < ntest; ++itest) {
// const SparseMatrix<double> R1 = Vinvsparse - VinvF*Cinvd.solve(FtVinv);
// const MatrixXd R2 = R1;
// }
// 
// auto time1 = std::chrono::high_resolution_clock::now();
// for (unsigned int itest = 0; itest < ntest; ++itest) {
// const MatrixXd R1 = Vinvsparse - VinvF*Cinvd.solve(FtVinv);
// const SparseMatrix<double> R2 = R1.sparseView();
// }
// 
// auto time2 = std::chrono::high_resolution_clock::now();
// 
// 
// const double d0 = (time1-time0).count();
// const double d1 = (time2-time1).count();
// // const double d2 = (time3-time2).count();
// // const double d3 = (time4-time3).count();
// // 
// std::cout << "d0 = " << d0 << " d1 = " << d1 << std::endl;
    
    
// std::vector<Triplet<Matrix<double, 3, 3>>> triplets;
// // SparseMatrix<Matrix<double, Dynamic, Dynamic, 0, 5, 5>> testsparse(100, 100);
// SparseMatrix<Matrix<double, 3, 3>> testsparse(100, 100);
// triplets.emplace_back(0, 0, Matrix<double, 3, 3>::Zero());
// 
// testsparse.setFromTriplets(triplets.begin(), triplets.end());
// testsparse.insert(0, 0) = Matrix<double, 3, 3>::Identity();
    
// std::cout << "R total size = " << R.rows()*R.cols() << " nonzeros = " << R.nonZeros() << std::endl;
    
// std::cout << "VinvF total size = " << VinvF.rows()*VinvF.cols() << " nonzeros = " << VinvF.nonZeros() << std::endl;


    //TODO make R dense (and find the best solution for sparse = sparse*dense case)
    
    
// chisqval = rfull.transpose()*R*rfull;
// 
// grad = 2.*Jsparse.transpose()*R*rfull;
// hess = 2.*Jsparse.transpose()*R*Jsparse;
    
    //TODO check this against explicit version
    const VectorXd Rr = Vinvsparse*(rfull + Fsparse*dxfree);
    
// const VectorXd Rrtest = R*rfull;
// const double Rrdiff = (Rr-Rrtest).array().square().sum();
// std::cout << "Rrdiff = " << Rrdiff << std::endl;
    
// constexpr unsigned int ntest = 1000;
// 
// auto time0 = std::chrono::high_resolution_clock::now();
// for (unsigned int itest = 0; itest < ntest; ++itest) {
// const VectorXd Rrtest = Vinvsparse*(rfull + Fsparse*dxfull);
// }
// 
// auto time1 = std::chrono::high_resolution_clock::now();
// for (unsigned int itest = 0; itest < ntest; ++itest) {
// // VectorXd Rrtest = Vinvsparse*rfull;
// // Rrtest += Vinvsparse*Fsparse*dxfull;
// const VectorXd Rrtest = R*rfull;
// }
// 
// auto time2 = std::chrono::high_resolution_clock::now();
// 
// 
// const double d0 = (time1-time0).count();
// const double d1 = (time2-time1).count();
// 
// std::cout << "d0 = " << d0 << " d1 = " << d1 << std::endl;
    
    
    chisqval = rfull.transpose()*Rr;
    
    grad = 2.*Jsparse.transpose()*Rr;
    hess = 2.*Jsparse.transpose()*R*Jsparse;
    
    
    dxdparms = MatrixXd::Zero(nparsfinal, nstateparms);

    dxdparms(Eigen::placeholders::all, freestateidxs) = -Cinvd.solve(VinvF.transpose()*Jsparse).transpose();

    // ---- Kink finder: decay-in-flight score test -------------------------
    // For each material step i, score the alternative hypothesis that the
    // propagated state acquires an unconstrained offset delta =
    // (dqop, ddxdz, ddydz) at the step end -- a decay kink plus momentum
    // step. With the residual model r_i(delta) = dx0_i - S*delta (S selects
    // components 0..2 of material block i), the chi2 gradient and
    // state-marginalized Hessian w.r.t. delta at the converged fit are
    //   g = -2 * Rr[block_i](0:3),   H = 2 * R[block_i, block_i](0:3, 0:3)
    // giving dchisq = g^T H^+ g / 2 (~ chi2(3) under the null) and best-fit
    // kink deltahat = -H^+ g, with no refit needed. The angle-only (2x2) and
    // qop-only (1x1) sub-tests separate hard elastic scatters (angle, no
    // momentum step) from decays (correlated angle + momentum step; for a
    // decay deltahat_qop * q > 0, the daughter is softer). H is singular in
    // directions the downstream hits cannot constrain (notably the last
    // step), handled by an eigenvalue-thresholded pseudo-inverse: fully
    // unconstrained directions contribute zero.
    // These 3 pseudo-parameters per step are deliberately NOT registered in
    // parmset/detidparms, so nothing leaks into the alignment-fit outputs
    // (jacrefv/gradv/hesspacked).
    if (doKinkFinder_) {
      const unsigned int nsteps = kinkConsIdx.size();
      kinkDchisq.assign(nsteps, 0.f);
      kinkDchisqAngle.assign(nsteps, 0.f);
      kinkDchisqQop.assign(nsteps, 0.f);
      kinkDqop.assign(nsteps, 0.f);
      kinkDxdz.assign(nsteps, 0.f);
      kinkDydz.assign(nsteps, 0.f);
      kinkGlobalR.assign(kinkStepR.begin(), kinkStepR.end());
      kinkGlobalZ.assign(kinkStepZ.begin(), kinkStepZ.end());
      kinkMax = -1.f;
      kinkMaxLayer = -1;

      // score test with pseudo-inverse: returns dchisq, fills deltahat
      const auto kinkScore = [](const MatrixXd &H, const VectorXd &g, VectorXd &deltahat) -> double {
        deltahat = VectorXd::Zero(g.size());
        SelfAdjointEigenSolver<MatrixXd> eig(H);
        const double lmax = eig.eigenvalues().cwiseMax(0.).maxCoeff();
        if (!(lmax > 0.)) {
          return 0.;
        }
        // Relative eigenvalue floor: directions this weakly constrained by
        // the downstream hits carry no usable kink information, and keeping
        // them amplifies noise as proj^2/lambda (observed as dchisq > total
        // chi2 on ~3e-3 of ideal-geometry tracks at 1e-8).
        const double thresh = 1e-6*lmax;
        double dchisq = 0.;
        for (Eigen::Index k = 0; k < eig.eigenvalues().size(); ++k) {
          const double lambda = eig.eigenvalues()(k);
          if (lambda > thresh) {
            const double proj = eig.eigenvectors().col(k).dot(g);
            dchisq += 0.5*proj*proj/lambda;
            deltahat -= (proj/lambda)*eig.eigenvectors().col(k);
          }
        }
        return dchisq;
      };

      for (unsigned int istep = 0; istep < nsteps; ++istep) {
        const unsigned int r0 = kinkConsIdx[istep];

        const MatrixXd H3 = 2.*R.block(r0, r0, 3, 3);
        const VectorXd g3 = -2.*Rr.segment(r0, 3);
        VectorXd d3;
        kinkDchisq[istep] = kinkScore(H3, g3, d3);
        kinkDqop[istep] = d3(0);
        kinkDxdz[istep] = d3(1);
        kinkDydz[istep] = d3(2);

        const MatrixXd H2 = 2.*R.block(r0 + 1, r0 + 1, 2, 2);
        const VectorXd g2 = -2.*Rr.segment(r0 + 1, 2);
        VectorXd d2;
        kinkDchisqAngle[istep] = kinkScore(H2, g2, d2);

        const MatrixXd H1 = 2.*R.block(r0, r0, 1, 1);
        const VectorXd g1 = -2.*Rr.segment(r0, 1);
        VectorXd d1;
        kinkDchisqQop[istep] = kinkScore(H1, g1, d1);

        if (kinkDchisq[istep] > kinkMax) {
          kinkMax = kinkDchisq[istep];
          kinkMaxLayer = istep;
        }
      }

      if (kinkInjectLayer_ >= 0) {
        std::cout << "kink injection closure: injected layer " << kinkInjectLayer_
                  << " delta = (" << kinkInjectDqop_ << ", " << kinkInjectDxdz_
                  << ", " << kinkInjectDydz_ << ")";
        if (kinkInjectLayer_ < int(nsteps)) {
          std::cout << "  recovered deltahat = (" << kinkDqop[kinkInjectLayer_]
                    << ", " << kinkDxdz[kinkInjectLayer_]
                    << ", " << kinkDydz[kinkInjectLayer_] << ")"
                    << "  dchisq(injected) = " << kinkDchisq[kinkInjectLayer_];
        }
        std::cout << "  kinkMax = " << kinkMax << " at layer " << kinkMaxLayer
                  << " of " << nsteps << std::endl;
      }
    }
    // ---- end kink finder -------------------------------------------------

    //additional contributions from resolution variations

    // exact block-eigenvalue export: eigenvalues of dV_b^{1/2} R_bb dV_b^{1/2}
    // per resolution entry (leg). Small dense blocks (<= 5x5); descending,
    // zero-padded to 5 floats. Offline validation: sum over the legs of a
    // parameter of sum(lambda) reproduces its gradllv entry.
    // The per-hit (complement) residual branches are reset for EVERY row,
    // not only inside the block that fills them: a track that never reaches
    // the export (a failed fit, a gen-frozen pass) must write empty arrays
    // rather than the previous track's.
    phresd = 0;
    phresnmeas = 0;
    phresnfree = 0;
    phreschi2 = 0.f;
    phresvchk = 0.f;
    phresrankgap = 0.f;
    phresgchk = 0.f;
    phresqrank = 0;
    phresnref = 0;
    phresok = false;
    phcfnok = 0;
    phcfms = 0.f;
    phcfgrpclosure = 0.f;
    phresz.clear();
    phresraw.clear();
    phresrow.clear();
    phreshit.clear();
    phresdim.clear();
    phrescls.clear();
    phrespiv.clear();
    phresinflat.clear();
    phresvarv.clear();
    phresbv.clear();
    phcfmsv.clear();
    phcfdelv.clear();
    phcfiorev.clear();
    phcfioimv.clear();
    phcfradrev.clear();
    phcfradimv.clear();
    phcfvgf.clear();
    phcfgrpcomp.clear();
    phcfgrpv.clear();
    phcfgrpmsv.clear();
    phcfgrpdelv.clear();
    phcfgrpiorev.clear();
    phcfgrpioimv.clear();
    phcfgrpradrev.clear();
    phcfgrpradimv.clear();
    phcfgrpvqms.clear();
    phcfgrpvqio.clear();
    phcfhitcomp.clear();
    phcfhitcls.clear();
    phcfhitv.clear();

    if (dores && fillTrackTree_ && (fillGrads_ || fillGradsFactored_)) {
      // Influence of the noise on the 5 reference parameters:
      // W5 = Vinv F C E5 (one solve with 5 RHS + one sparse matmul); zero
      // in gen-frozen fits, where the reference state is not free and the
      // noise never propagates to it.
      MatrixXd W5 = MatrixXd::Zero(ncons, 5);
      {
        MatrixXd E5 = MatrixXd::Zero(nstatefree, 5);
        bool anyfree = false;
        for (unsigned int i = 0; i < nstatefree; ++i) {
          if (freestateidxs[i] < 5) {
            E5(i, freestateidxs[i]) = 1.;
            anyfree = true;
          }
        }
        if (anyfree) {
          W5 = VinvF*Cinvd.solve(E5);
        }
      }
      const VectorXd wqop = W5.col(0);

      for (unsigned int ires = 0; ires < dVs.size(); ++ires) {
        const unsigned int r0 = resblockrng[ires][0];
        const unsigned int nb = resblockrng[ires][1];
        const MatrixXd dVb = MatrixXd(dVs[ires]).block(r0, r0, nb, nb);
        SelfAdjointEigenSolver<MatrixXd> eigv(dVb);
        const MatrixXd sqrtdV = eigv.eigenvectors() *
                                eigv.eigenvalues().cwiseMax(0.).cwiseSqrt().asDiagonal() *
                                eigv.eigenvectors().transpose();
        const MatrixXd B = sqrtdV * R.block(r0, r0, nb, nb) * sqrtdV;
        SelfAdjointEigenSolver<MatrixXd> eigB(B);
        reseigidx.push_back(resglobidx[ires]);
        // which valid hit this block belongs to (-1 = material). Without it
        // the parmtype-8/9 blocks can only be matched to hits by guessing the
        // ordering, which breaks the moment a propagation fails.
        reshitidx.push_back(resvalidhit_[ires]);
        // The hit's CLASS, from the per-hit variables this maker already
        // exports (`hitres_classes.class_of`, index into the canonical
        // 18-entry list). It is redundant here -- an offline reader can form
        // it from `reshitidx` + hitDetId/clusterSizeX/hitUProj/
        // clusterChargeBin -- but writing it means the two makers' trees
        // carry the SAME quantity under the same name, and the two-track one
        // has no per-hit variables to form it from.
        {
          const int fam = ires < resfamily_.size() ? resfamily_[ires] : -1;
          const int ih = resvalidhit_[ires];
          int cls = -1;
          if ((fam == 8 || fam == 9) && ih >= 0 && std::size_t(ih) < hitDetId.size() &&
              std::size_t(ih) < clusterSizeX.size() && std::size_t(ih) < hitUProj.size() &&
              std::size_t(ih) < clusterChargeBin.size()) {
            cls = hitResClassIndex(int((hitDetId[ih] >> 25) & 0x7), clusterSizeX[ih], hitUProj[ih],
                                   clusterChargeBin[ih], fam == 9);
          }
          reshitcls.push_back(static_cast<short>(cls));
        }
        for (unsigned int j = 0; j < 5; ++j) {
          reseigv.push_back(j < nb ? std::max(eigB.eigenvalues()(nb - 1 - j), 0.) : 0.f);
        }
        // influence-weight export, aligned with reseigidx: raw dof weights
        // and the block's variance contribution to the fitted q/p
        const double vb = wqop.segment(r0, nb).transpose() * dVb * wqop.segment(r0, nb);
        resinfvarv.push_back(vb);
        {
          const int fam = ires < resfamily_.size() ? resfamily_[ires] : -1;
          // The parmtype-15 blocks are a RE-PARTITION of the parmtype-10/11
          // noise, not an addition to it, so they must NOT enter `resinfcov`:
          // that would double-count the material share and break the offline
          // coverage cut `|resinfcov/refCov(0,0) - 1| < 5e-3`.
          if (fam == 15) {
            resinfcovgrp += float(vb);
          } else {
            resinfcov += vb;
          }
          // The hit share on its own, so the two trees expose the same
          // decomposition. Unlike the two-track maker, `resinfcov` HERE has
          // always included the hit blocks, and that is left alone.
          if (fam == 8 || fam == 9) {
            resinfcovhit += float(vb);
          }
        }
        for (unsigned int j = 0; j < 5; ++j) {
          resinfv.push_back(j < nb ? wqop(r0 + j) : 0.f);
        }
        // generalized functional export: B_b = M_b dV_b^{1/2} (5 x nb,
        // row-major, dof-padded to 5); reuses sqrtdV from the eigenvalue
        // block above
        const MatrixXd Bb = W5.block(r0, 0, nb, 5).transpose() * sqrtdV;
        for (unsigned int p = 0; p < 5; ++p) {
          for (unsigned int j = 0; j < 5; ++j) {
            resinfbv.push_back(j < nb ? Bb(p, j) : 0.f);
          }
        }
      }

      // ---- THE RESOLUTION-CF EXPONENTS, for the q/p functional -----------
      //
      // Here and not offline because the block weights `sqrt(v_b/sq2)/sigma`
      // are the FIT's own influence coefficients and do not exist until it has
      // converged -- which is also why the raw step records had to be exported
      // at all. Everything the offline extractor reads is in scope now, so the
      // 6 x 64 floats it would have spent 2.2 s and 430 kB producing cost a
      // few ms here.
      //
      // The inputs are the EXPORT ARRAYS, not the propagator's logs, so the
      // pooling is identical to `cf_track_resolution.extract`'s join by
      // construction and cannot drift from it.
      cvhcf::TrackInput cfin;
      // sigma from the FLOAT the tree carries, so the in-maker weight is
      // exactly the one a reader of the same file would have recovered.
      // Hoisted out of the fill block below because both consumers need it.
      const double c00 = refCov[0];
      if (exportCfExponents_ || exportPerHitResidual_) {
        cfin.resglobidx = resglobidx.data();
        cfin.resfamily = resfamily_.data();
        cfin.resvarv = resinfvarv.data();
        cfin.nres = int(std::min({resglobidx.size(), resfamily_.size(), resinfvarv.size()}));
        cfin.ms = {msmoliidx.data(), msmoliv.data(), int(msmoliidx.size()),
                   msmoliidx.empty() ? 0 : int(msmoliv.size() / msmoliidx.size())};
        cfin.ioni = {ioniurbanidx.data(), ioniurbanv.data(), int(ioniurbanidx.size()),
                     ioniurbanidx.empty() ? 0 : int(ioniurbanv.size() / ioniurbanidx.size())};
        cfin.qsc = {ioniqscaleidx.data(), ioniqscalev.data(), int(ioniqscaleidx.size()), 2};
        cfin.rad = {radstepidx.data(), radstepv.data(), int(radstepidx.size()), RADSTEP_STRIDE};
        // The material-group column of each record: `msmoliv` has carried it
        // at column 9 since the global material model landed; `ioniurbanv`
        // and `radstepv` carry it as their LAST column (2026-09-06). Set
        // explicitly rather than inferred so that a stride change cannot
        // silently relabel a group as a physics quantity.
        cfin.ms.groupCol = cfin.ms.stride >= 10 ? 9 : -1;
        cfin.ioni.groupCol = cfin.ioni.stride - 1;
        cfin.rad.groupCol = RADSTEP_STRIDE - 1;
        cfin.radspec = radstepspecv.data();
        cfin.radvgrid = radvgrid.data();
        cfin.radnv = int(radvgrid.size());
        cfin.sigma = c00 > 0. ? std::sqrt(c00) : 0.;
        // THE CHARGE. `ioniurbanv`'s cs = E/p^3 is positive for every track
        // and the physical map is d(q/p) = q cs dE, so the ionization (and
        // radiative) step weight is charge-signed -- the same factor
        // `Geant4ePropagator` puts into the in-fit CGF block's `gs`.
        cfin.ioniSign = refParms[0] >= 0.f ? 1. : -1.;
        cfin.wantDelta = true;
        cfin.wantGroups = exportCfGroupExponents_;
        cfin.wantGroupDelta = true;   // the q/p functional's model uses S_del
      }
      if (exportCfExponents_) {
        cvhcf::TrackResult cfres;
        cvhcf::trackExponents(cfin, cfres);
        cfok = cfres.ok;
        cfnblock = cfres.nblockms + cfres.nblockioni;
        cfnpooled = cfres.npooled;
        cfvgf = (c00 > 0.) ? float(cfres.vgauss / c00) : 0.f;
        auto storecf = [](const std::array<double, cvhcf::kNTau> &a, std::vector<float> &v) {
          v.resize(cvhcf::kNTau);
          for (int j = 0; j < cvhcf::kNTau; ++j)
            v[j] = float(a[j]);
        };
        storecf(cfres.S.ms, cfmsv);
        storecf(cfres.S.del, cfdelv);
        storecf(cfres.S.ioRe, cfiorev);
        storecf(cfres.S.ioIm, cfioimv);
        storecf(cfres.S.radRe, cfradrev);
        storecf(cfres.S.radIm, cfradimv);
        storeCfGroups(cfres);
        // Per-hit-class Gaussian shares of the q/p variance, ascending in
        // class. Here `vgauss` (and hence `cfqop_vgf`) IS the sum over the
        // parmtype-8/9 blocks, so `sum_c cfqop_hitv == cfqop_vgf` exactly.
        {
          std::array<double, kNHitResClasses> vcls{};
          bool anycls = false;
          for (std::size_t i = 0; i < reshitcls.size() && i < resinfvarv.size(); ++i) {
            const int c = reshitcls[i];
            if (c < 0 || c >= kNHitResClasses) {
              continue;
            }
            vcls[c] += resinfvarv[i];
            anycls = true;
          }
          if (anycls && c00 > 0.) {
            for (int c = 0; c < kNHitResClasses; ++c) {
              if (vcls[c] == 0.) {
                continue;
              }
              cfhitclsv.push_back(static_cast<short>(c));
              cfhitvv.push_back(float(vcls[c] / c00));
            }
          }
        }
      }

      // ---- THE PER-HIT (COMPLEMENT) RESIDUAL VECTOR ----------------------
      //
      // The q/p functional above, and the whole truth-referenced prototype it
      // feeds, need `genParms`.  On DATA there is none, and what is left is
      // the part of the constraint residual the fit has NOT absorbed:
      //
      //     rho = V R r = r + F dxfree ,  Cov(rho) = V R V = V - F C F^T
      //
      // of rank `d = ncons - nstatefree`, which for this fit
      // (`ncons = 5 nhits + nvalid + nvalidpixel`, `nstateparms = 5(nhits+1)`)
      // is exactly `n_meas - 5`.  The kink rows of `rho` are a deterministic
      // function of the measurement rows, and the Mahalanobis form of a
      // Gaussian is invariant under a bijective map of its support, so
      // restricting to the MEASUREMENT rows loses nothing: with
      // `G = V_mm - F_m C F_m^T` (rank d) and any `Cw Cw^T = G`,
      // `z = Cw^+ rho_m` has `Cov(z) = I_d` and `sum_k z_k^2 = r^T R r`,
      // the fit's own chi2.  Both identities are exported as gates
      // (`phres_d`, `phres_chi2`, `phres_vchk`).
      //
      // THE BASIS is the LDL^T of `G` in MEASUREMENT-ROW ORDER, i.e. hit
      // order inner to outer (and, within a pixel, the first local
      // coordinate then the second), with the `n_meas - d = 5` null pivots
      // skipped.  Component k is then that row's post-fit residual
      // CONDITIONED ON the inner rows' -- local to a hit, which is what a
      // term whose parameters are per-hit-class needs.  It is deliberately
      // NOT the Kalman filter innovation sequence: that is a different
      // orthonormal basis of the same d-space (it conditions on the RAW
      // inner measurements and drops the FIRST five components rather than
      // the last five), it needs a sequential filter pass that nothing here
      // has, and the product-of-marginals likelihood the offline term forms
      // is basis dependent, so the choice is stated and its cross-dependence
      // measured rather than assumed away.
      if (exportPerHitResidual_) {
        // (1) the measurement rows in hit order, from the parmtype-8 blocks,
        //     with the hit index, the local coordinate and the hit class.
        std::vector<unsigned int> mrow;
        std::vector<short> mhit, mdim, mcls;
        mrow.reserve(nvalid + nvalidpixel);
        for (unsigned int ires = 0; ires < resfamily_.size(); ++ires) {
          if (resfamily_[ires] != 8) {
            continue;
          }
          const unsigned int r0 = resblockrng[ires][0];
          const unsigned int nb = resblockrng[ires][1];
          for (unsigned int j = 0; j < nb; ++j) {
            mrow.push_back(r0 + j);
            mhit.push_back(static_cast<short>(resvalidhit_[ires]));
            mdim.push_back(static_cast<short>(j));
            mcls.push_back(-1);
          }
        }
        // the class of the SECOND pixel coordinate lives on the parmtype-9
        // entry of the same hit, so the join is on (hit, coordinate).
        for (unsigned int ires = 0; ires < resfamily_.size(); ++ires) {
          const int fam = resfamily_[ires];
          if (fam != 8 && fam != 9) {
            continue;
          }
          const short want = (fam == 8) ? 0 : 1;
          const int ih = resvalidhit_[ires];
          const short cls = ires < reshitcls.size() ? reshitcls[ires] : short(-1);
          for (std::size_t i = 0; i < mrow.size(); ++i) {
            if (mhit[i] == ih && mdim[i] == want) {
              mcls[i] = cls;
            }
          }
        }
        const int nm = static_cast<int>(mrow.size());
        phresnmeas = nm;
        phresnfree = static_cast<int>(nstatefree);

        if (nm > 5 && nstatefree > 0 && rfull.size() == Eigen::Index(ncons)) {
          const int nfree = static_cast<int>(nstatefree);

          // (2) THE WHITENING.
          //
          // The obvious route -- assemble `G = V_mm - F_m C F_m^T` and
          // factorize it -- does not work, and it is worth saying why,
          // because it looks like it should.  `G` is a Schur complement: its
          // five null directions are a DIFFERENCE of two nearly equal
          // numbers, and the solve they come from (`C = (F^T V^-1 F)^-1`) is
          // badly conditioned, a thin layer giving its kink block very small
          // process noise.  MEASURED on a 200-track smoke: the null
          // eigenvalues of the assembled `G` come out at ~1e-9 of the
          // diagonal instead of 1e-16, `max_k |(Psi^T G Psi)_kk - 1|` reaches
          // 1e9, and the LAST whitened component of a track -- the one that
          // closes the d-dimensional space -- misses `sum_b v^(k)_b = 1` by
          // up to a factor of 19.  One step of iterative refinement made it
          // worse, not better.
          //
          // So the construction is done in the STANDARDIZED space with
          // ORTHOGONAL operations only.  With `Fw = V^-1/2 F`,
          //     R = V^-1/2 (I - P) V^-1/2 ,   P = Fw (Fw^T Fw)^-1 Fw^T
          //                                     = Q1 Q1^T
          // is an identity, so the projector is a QR of `Fw` and never an
          // inverse:
          //   * `s = (I - Q1 Q1^T) V^-1/2 r` is the standardized post-fit
          //     residual and `chi2 = |s|^2 = r^T R r`;
          //   * `Gs = I - Q1_m Q1_m^T` is its covariance on the measurement
          //     rows, and its null eigenvalues now sit at 1e-16 of the unit
          //     diagonal because `Q1_m Q1_m^T` is a submatrix of an
          //     orthogonal projector and not a cancelling product of solves.
          // `V` is block diagonal over exactly the blocks `resblockrng`
          // enumerates (5 process-noise rows per propagation, 1 or 2 per
          // valid hit), so `V^-1/2` is a per-block symmetric square root of
          // `Vinvfull` and costs nothing.
          std::vector<std::array<unsigned int, 2>> vblk;
          {
            std::vector<char> cov(ncons, 0);
            for (unsigned int ires = 0; ires < resfamily_.size(); ++ires) {
              const int fam = resfamily_[ires];
              if (fam != 8 && fam != 10) {
                continue;
              }
              const unsigned int r0 = resblockrng[ires][0];
              const unsigned int nb = resblockrng[ires][1];
              vblk.push_back({{r0, nb}});
              for (unsigned int j = 0; j < nb; ++j) {
                cov[r0 + j] = 1;
              }
            }
            // whatever the resolution families do not cover (the beamspot
            // constraint rows) becomes its own contiguous block.
            unsigned int i = 0;
            while (i < unsigned(ncons)) {
              if (cov[i]) {
                ++i;
                continue;
              }
              unsigned int j = i;
              while (j < unsigned(ncons) && !cov[j]) {
                ++j;
              }
              vblk.push_back({{i, j - i}});
              i = j;
            }
          }
          MatrixXd Vih = MatrixXd::Zero(ncons, ncons);
          for (auto const &b : vblk) {
            const unsigned int r0 = b[0];
            const unsigned int nb = b[1];
            const MatrixXd Vib = Vinvfull.block(r0, r0, nb, nb);
            SelfAdjointEigenSolver<MatrixXd> esb(0.5*(Vib + Vib.transpose()));
            Vih.block(r0, r0, nb, nb) =
                esb.eigenvectors()*esb.eigenvalues().cwiseMax(0.).cwiseSqrt().asDiagonal()*
                esb.eigenvectors().transpose();
          }

          MatrixXd Ffree(ncons, nfree);
          for (int i = 0; i < int(ncons); ++i) {
            for (int j = 0; j < nfree; ++j) {
              Ffree(i, j) = Ffull(i, freestateidxs[j]);
            }
          }
          const MatrixXd Fw = Vih*Ffree;
          const VectorXd sfull = Vih*rfull;
          ColPivHouseholderQR<MatrixXd> qrFw(Fw);
          qrFw.setThreshold(1e-12);
          const int qrank = int(qrFw.rank());
          phresqrank = qrank;
          const MatrixXd Q1 = qrFw.householderQ()*MatrixXd::Identity(ncons, qrank);
          const VectorXd qts = Q1.transpose()*sfull;
          const VectorXd sstd = sfull - Q1*qts;

          MatrixXd Q1m(nm, qrank);
          VectorXd sm(nm);
          for (int i = 0; i < nm; ++i) {
            Q1m.row(i) = Q1.row(mrow[i]);
            sm(i) = sstd(mrow[i]);
          }
          MatrixXd G = -(Q1m*Q1m.transpose());
          G.diagonal().array() += 1.;

          // (3) hit-order factorization of `Gs`, rank imposed at
          //     `dexp = ncons - rank(Fw)`.  `S = U sqrt(Lambda)` truncated to
          //     `dexp`, then modified Gram-Schmidt (with one
          //     re-orthogonalization pass) on the ROWS of `S` in
          //     measurement-row order: row k contributes a component iff it
          //     is not already in the span of the inner rows, which is a test
          //     on a NORM.  What comes out is the LDL^T's own triangular
          //     whitener, `z_k = (s_k - sum_{j<k} c_kj z_j)/|e_k|`, computed
          //     stably.  `phres_rankgap` is lambda_dexp/lambda_(dexp+1), so
          //     the rank decision is auditable rather than trusted.
          int dexp = int(ncons) - qrank;
          if (dexp > nm) {
            dexp = nm;
          }
          std::vector<double> piv(nm, 0.);
          std::vector<int> keep;
          MatrixXd Psi;
          VectorXd zvec;
          if (dexp > 0) {
            SelfAdjointEigenSolver<MatrixXd> esG(G);
            const VectorXd ev = esG.eigenvalues();  // ascending
            MatrixXd S(nm, dexp);
            for (int c = 0; c < dexp; ++c) {
              const int idx = nm - 1 - c;
              S.col(c) = esG.eigenvectors().col(idx)*std::sqrt(std::max(ev(idx), 0.));
            }
            phresrankgap = (dexp < nm && ev(nm - dexp) > 0.)
                               ? float(ev(nm - dexp)/std::max(std::abs(ev(nm - dexp - 1)), 1e-300))
                               : 0.f;

            MatrixXd Qm = MatrixXd::Zero(dexp, dexp);
            MatrixXd T = MatrixXd::Zero(dexp, nm);
            int nkeep = 0;
            const double gtol = 1e-12;
            for (int k = 0; k < nm && nkeep < dexp; ++k) {
              VectorXd e = S.row(k).transpose();
              VectorXd trow = VectorXd::Zero(nm);
              trow(k) = 1.;
              const double s0 = e.squaredNorm();
              for (int pass = 0; pass < 2; ++pass) {
                for (int j = 0; j < nkeep; ++j) {
                  const double c = Qm.row(j).dot(e);
                  e -= c*Qm.row(j).transpose();
                  trow -= c*T.row(j).transpose();
                }
              }
              const double n2 = e.squaredNorm();
              if (!(n2 > gtol*std::max(s0, 1e-300))) {
                continue;
              }
              const double nn = std::sqrt(n2);
              Qm.row(nkeep) = e.transpose()/nn;
              T.row(nkeep) = trow.transpose()/nn;
              piv[k] = n2;
              keep.push_back(k);
              ++nkeep;
            }
            if (nkeep > 0) {
              Psi = T.topRows(nkeep).transpose();
              zvec = T.topRows(nkeep)*sm;
              // the whitener's own closure, against the SAME `Gs` it was
              // built from: this separates a defect in the factorization from
              // a defect in the influence `W`, which additionally goes
              // through `V^-1/2` and `Q1`.
              const MatrixXd PGP = Psi.transpose()*G*Psi;
              double gc = 0.;
              for (int kk = 0; kk < nkeep; ++kk) {
                gc = std::max(gc, std::abs(PGP(kk, kk) - 1.));
              }
              phresgchk = float(gc);
            }
          }
          const int nd = static_cast<int>(keep.size());
          phresd = nd;

          if (nd > 0) {
            // (4) the per-component influence.  `z_k = Wstd[:,k]^T (V^-1/2 n)`
            //     with `Wstd = E_m Psi - Q1 Q1_m^T Psi`, so in the noise's own
            //     units `W = V^-1/2 Wstd`, and
            //     `sum_b W_b^T dV_b W_b = |Wstd_k|^2 = 1` by construction of
            //     the whitener -- now entirely out of an orthogonal projector.
            MatrixXd Wstd = -(Q1*(Q1m.transpose()*Psi));
            for (int i = 0; i < nm; ++i) {
              Wstd.row(mrow[i]) += Psi.row(i);
            }
            const MatrixXd W = Vih*Wstd;

            double chi2z = 0.;
            for (int kk = 0; kk < nd; ++kk) {
              const int k = keep[kk];
              const double zk = zvec(kk);
              chi2z += zk*zk;
              phresz.push_back(float(zk));
              phresrow.push_back(static_cast<short>(k));
              phreshit.push_back(mhit[k]);
              phresdim.push_back(mdim[k]);
              phrescls.push_back(mcls[k]);
              phrespiv.push_back(float(piv[k]));
              phresinflat.push_back(float(G(k, k)/piv[k]));
            }
            phreschi2 = float(chi2z);
            // the STANDARDIZED post-fit residual before whitening, in
            // measurement-row order: dimensionless, and what a
            // density-vs-model plot wants.
            phresraw.reserve(nm);
            for (int i = 0; i < nm; ++i) {
              phresraw.push_back(float(sm(i)));
            }

            // (3b) THE TRUTH-REFERENCED COMPONENTS, appended to the same
            //      arrays.  `r = refParms - genParms` whitened by the lower
            //      Cholesky factor of `refCov` is exactly the prototype of
            //      `calibration_studies/resolution/hitlik`, and its influence
            //      is `W5 L^-T` -- the SAME `W5` the q/p export above already
            //      built.  Carrying it here rather than in a second
            //      production means the two terms are read out of one file
            //      with one convention, which is what the complementarity
            //      test (`F^T R = 0` -> the two are uncorrelated, so their
            //      information must add) and the joint fit need.
            //
            //      Component 0 of this set IS the established q/p functional:
            //      column 0 of `L^-T` is `e_0/sigma_qp`, so its per-block
            //      weight `sqrt(v^(0)_b/sq2)` equals the `cfqop_*` weight
            //      `sqrt(v_b/sq2)/sigma` term by term.  That makes the
            //      agreement of `phcf_*` component `d` with `cfqop_*` a
            //      closure test of the whole per-component route, and it is
            //      why nothing here is a second implementation of the model.
            //
            //      On DATA there is no `genParms`, `nref` is 0, and only the
            //      per-hit components survive -- which is the point of the
            //      whole exercise.
            MatrixXd Wall = W;
            VectorXd zall(nd);
            for (int kk = 0; kk < nd; ++kk) {
              zall(kk) = zvec(kk);
            }
            int nref = 0;
            if (perHitRefComponents_ && genpart != nullptr) {
              Matrix<double, 5, 5> Cref;
              for (unsigned int a = 0; a < 5; ++a) {
                for (unsigned int b = 0; b < 5; ++b) {
                  Cref(a, b) = double(refCov[5*std::min(a, b) + std::max(a, b)]);
                }
              }
              const LLT<Matrix<double, 5, 5>> llt(Cref);
              if (llt.info() == Eigen::Success) {
                const Matrix<double, 5, 5> Lref = llt.matrixL();
                Matrix<double, 5, 1> rres;
                for (unsigned int a = 0; a < 5; ++a) {
                  rres(a) = double(refParms[a]) - double(genParms[a]);
                }
                const Matrix<double, 5, 1> zref = Lref.triangularView<Eigen::Lower>().solve(rres);
                // W5 L^-T : the whitened influence, column j of which is the
                // influence of the noise on `zref(j)`.
                const Matrix<double, 5, 5> LinvT =
                    Lref.transpose().triangularView<Eigen::Upper>().solve(
                        Matrix<double, 5, 5>::Identity());
                const MatrixXd Wref = W5*LinvT;
                nref = 5;
                Wall.conservativeResize(ncons, nd + nref);
                Wall.rightCols(nref) = Wref;
                zall.conservativeResize(nd + nref);
                for (int j = 0; j < nref; ++j) {
                  zall(nd + j) = zref(j);
                  phresz.push_back(float(zref(j)));
                  phresrow.push_back(static_cast<short>(-1 - j));
                  phreshit.push_back(-1);
                  phresdim.push_back(-1);
                  phrescls.push_back(-1);
                  phrespiv.push_back(1.f);
                  phresinflat.push_back(float(Cref(j, j)/std::max(Lref(j, j)*Lref(j, j), 1e-300)));
                }
              }
            }
            phresnref = nref;
            const int ntot = nd + nref;

            // (4) the per-(block, component) variance shares, signed by the
            //     qop-row influence.  `sum_b v^(k)_b == 1` over the
            //     non-parmtype-15 blocks (15 is a RE-PARTITION of 10/11, not
            //     an addition to it) -- the export's own closure test.
            const std::size_t nres = dVs.size();
            std::vector<std::vector<float>> shares(ntot, std::vector<float>(nres, 0.f));
            std::vector<std::vector<float>> signs(ntot, std::vector<float>(nres, 1.f));
            std::vector<double> vsum(ntot, 0.);
            phresvarv.assign(static_cast<std::size_t>(ntot)*nres, 0.f);
            for (std::size_t ires = 0; ires < nres; ++ires) {
              const unsigned int r0 = resblockrng[ires][0];
              const unsigned int nb = resblockrng[ires][1];
              MatrixXd dVb = MatrixXd::Zero(nb, nb);
              for (int c = 0; c < dVs[ires].outerSize(); ++c) {
                for (SparseMatrix<double>::InnerIterator it(dVs[ires], c); it; ++it) {
                  const int rr = int(it.row()) - int(r0);
                  const int cc = int(it.col()) - int(r0);
                  if (rr >= 0 && rr < int(nb) && cc >= 0 && cc < int(nb)) {
                    dVb(rr, cc) = it.value();
                  }
                }
              }
              const int fam = ires < resfamily_.size() ? resfamily_[ires] : -1;
              // THE PER-COMPONENT INFLUENCE VECTORS `A_b = W_b^T dV_b^{1/2}`
              // (ntot x nb, dof-padded to 5), the exact analogue of the q/p
              // export's `resinfbv`.  `|A_b[k]|^2` is the variance share
              // below, so this adds nothing to the likelihood -- what it adds
              // is `A_b[j] . A_b[k]`, and with it the FOURTH CROSS CUMULANT
              // between components, which is the size of what the
              // product-of-marginals likelihood drops.  For a multiple-
              // scattering block the two projected angles are iso-Gaussian
              // plus a common radial tail, so
              //   kappa(z_j,z_j,z_k,z_k) ~ (A_j.A_k)^2 + |A_j|^2|A_k|^2/2
              // against a diagonal 3/2 |A_j|^4, and neither term is available
              // from the variance shares alone.
              MatrixXd sqrtdVb;
              if (perHitInfluenceBlocks_) {
                SelfAdjointEigenSolver<MatrixXd> esv(dVb);
                // dV of a pixel's parmtype-8 block is INDEFINITE (zero on the
                // yy diagonal, the xy correlation off it), so the negative
                // eigenvalue is clamped exactly as the `resinfbv` export
                // already does. |Vxy| is small, hence the 1e-7 closure of
                // `sum_b B_b B_b^T` against `refCov`.
                sqrtdVb = esv.eigenvectors()*
                          esv.eigenvalues().cwiseMax(0.).cwiseSqrt().asDiagonal()*
                          esv.eigenvectors().transpose();
              }
              // THE BLOCK'S NOISE DIRECTION, and with it the SIGN of its
              // weight.  The ionization (and radiative) CF is not even in
              // its weight -- an energy loss goes one way -- so each
              // component needs its own sign.  `dQI` is rank one: in the
              // curvilinear frame it is `e_0 e_0^T sigma^2`, and the
              // local-frame rotation `Hm` turns that into `u u^T sigma^2`
              // with `u = Hm e_0` spread over all five rows, so the signed
              // coefficient is `W_b . u` and NOT `W[r0, k]`.  `u` is oriented
              // by its qop component, and the product carries one further
              // MINUS because the process-noise constraint row is built as
              // (propagated - state) rather than the other way round.
              //
              // None of that is asserted: reference component 0 IS the q/p
              // functional (see (3b)), so requiring it to reproduce the
              // validated `cfqop_ioni_im` / `cfqop_rad_im` fixes the
              // convention, and the three candidate conventions are
              // distinguished by their worst-track relative difference --
              // `W[r0,k]` with no minus: 2.0 (i.e. the opposite sign, every
              // track); `W[r0,k]` with the minus: 5.3e-2; `W_b . u` with the
              // minus: 1.2e-7, which is the float32 storage of the
              // reference.  The real families (`ms`, `del`, `ioni_re`,
              // `rad_re`) agree at 1.2e-7 under all three, as they must --
              // they are even.
              VectorXd udir = VectorXd::Zero(nb);
              if (fam == 11 && nb > 0) {
                SelfAdjointEigenSolver<MatrixXd> esb(dVb);
                int imax = 0;
                for (int q = 1; q < int(nb); ++q) {
                  if (std::abs(esb.eigenvalues()(q)) > std::abs(esb.eigenvalues()(imax))) {
                    imax = q;
                  }
                }
                udir = esb.eigenvectors().col(imax);
                if (udir(0) < 0.) {
                  udir = -udir;
                }
              }
              for (int kk = 0; kk < ntot; ++kk) {
                const VectorXd wk = Wall.block(r0, kk, nb, 1);
                double v = wk.transpose()*dVb*wk;
                if (!(v > 0.)) {
                  v = 0.;
                }
                const float sg = (fam == 11 && wk.dot(udir) > 0.) ? -1.f : 1.f;
                shares[kk][ires] = float(v);
                signs[kk][ires] = sg;
                phresvarv[static_cast<std::size_t>(kk)*nres + ires] = float(sg*v);
                if (perHitInfluenceBlocks_) {
                  const VectorXd ab = sqrtdVb*wk;
                  for (unsigned int j = 0; j < 5; ++j) {
                    phresbv.push_back(j < nb ? float(ab(j)) : 0.f);
                  }
                }
                if (fam != 15) {
                  vsum[kk] += v;
                }
              }
            }
            double vchk = 0.;
            for (int kk = 0; kk < ntot; ++kk) {
              vchk = std::max(vchk, std::abs(vsum[kk] - 1.));
            }
            phresvchk = float(vchk);

            // (5) the CF exponents, one `cvhcf` pass per component.  The
            //     standardization is 1 (z already has unit variance), so the
            //     block weight is `sqrt(v^(k)_b / sq2)` -- exactly what
            //     `extract_res5.py` evaluates offline.  `phcf_msec` is the
            //     wall clock of this loop, to be read against the ~2 s the
            //     fit itself costs.
            const auto tcf0 = std::chrono::steady_clock::now();
            const int nresi = int(std::min(resglobidx.size(), resfamily_.size()));
            double grpclos = 0.;
            for (int kk = 0; kk < ntot; ++kk) {
              cvhcf::TrackInput ci = cfin;
              ci.sigma = 1.;
              ci.ioniSign = refParms[0] >= 0.f ? 1. : -1.;
              if (perHitShareMin_ > 0.) {
                for (auto &x : shares[kk]) {
                  if (double(x) < perHitShareMin_) {
                    x = 0.f;
                  }
                }
              }
              ci.resvarv = shares[kk].data();
              ci.ressgn = signs[kk].data();
              ci.nres = nresi;
              ci.wantDelta = true;
              ci.wantGroups = perHitCfGroups_;
              ci.wantGroupDelta = true;
              cvhcf::TrackResult cr;
              cvhcf::trackExponents(ci, cr);
              if (cr.ok) {
                ++phcfnok;
              }
              for (int j = 0; j < cvhcf::kNTau; ++j) {
                phcfmsv.push_back(float(cr.S.ms[j]));
                phcfdelv.push_back(float(cr.S.del[j]));
                phcfiorev.push_back(float(cr.S.ioRe[j]));
                phcfioimv.push_back(float(cr.S.ioIm[j]));
                phcfradrev.push_back(float(cr.S.radRe[j]));
                phcfradimv.push_back(float(cr.S.radIm[j]));
              }
              phcfvgf.push_back(float(cr.vgauss));
              if (perHitCfGroups_) {
                std::array<double, cvhcf::kNTau> sms{}, sdel{}, sre{}, sim{}, srre{}, srim{};
                for (auto const &gg : cr.groups) {
                  phcfgrpcomp.push_back(static_cast<short>(kk));
                  phcfgrpv.push_back(static_cast<short>(gg.group));
                  phcfgrpvqms.push_back(float(gg.vqms));
                  phcfgrpvqio.push_back(float(gg.vqio));
                  for (int j = 0; j < cvhcf::kNTau; ++j) {
                    phcfgrpmsv.push_back(float(gg.S.ms[j]));
                    phcfgrpdelv.push_back(float(gg.S.del[j]));
                    phcfgrpiorev.push_back(float(gg.S.ioRe[j]));
                    phcfgrpioimv.push_back(float(gg.S.ioIm[j]));
                    phcfgrpradrev.push_back(float(gg.S.radRe[j]));
                    phcfgrpradimv.push_back(float(gg.S.radIm[j]));
                    sms[j] += gg.S.ms[j];
                    sdel[j] += gg.S.del[j];
                    sre[j] += gg.S.ioRe[j];
                    sim[j] += gg.S.ioIm[j];
                    srre[j] += gg.S.radRe[j];
                    srim[j] += gg.S.radIm[j];
                  }
                }
                double dmax = 0., smax = 0.;
                auto cmp = [&](const std::array<double, cvhcf::kNTau> &a,
                               const std::array<double, cvhcf::kNTau> &b) {
                  for (int j = 0; j < cvhcf::kNTau; ++j) {
                    dmax = std::max(dmax, std::abs(a[j] - b[j]));
                    smax = std::max(smax, std::abs(b[j]));
                  }
                };
                cmp(sms, cr.S.ms);
                cmp(sdel, cr.S.del);
                cmp(sre, cr.S.ioRe);
                cmp(sim, cr.S.ioIm);
                cmp(srre, cr.S.radRe);
                cmp(srim, cr.S.radIm);
                if (smax > 0.) {
                  grpclos = std::max(grpclos, dmax/smax);
                }
              }
              // the component's Gaussian variance by hit class; sums to
              // `phcf_vgf` of the same component by construction.
              {
                std::array<double, kNHitResClasses> vcls{};
                for (std::size_t i = 0; i < reshitcls.size() && i < shares[kk].size(); ++i) {
                  const int c = reshitcls[i];
                  if (c < 0 || c >= kNHitResClasses) {
                    continue;
                  }
                  vcls[c] += shares[kk][i];
                }
                for (int c = 0; c < kNHitResClasses; ++c) {
                  if (vcls[c] == 0.) {
                    continue;
                  }
                  phcfhitcomp.push_back(static_cast<short>(kk));
                  phcfhitcls.push_back(static_cast<short>(c));
                  phcfhitv.push_back(float(vcls[c]));
                }
              }
            }
            phcfgrpclosure = float(grpclos);
            phcfms = float(std::chrono::duration<double, std::milli>(
                               std::chrono::steady_clock::now() - tcf0).count());
            phresok = (nd == dexp) && (nd == nm - 5) && (qrank == nfree) && (phcfnok == ntot);
          }
        }
      }
    }

    std::vector<SparseMatrix<double>> dVRs;
    dVRs.reserve(dVs.size());
    
    MatrixXd dVRr = MatrixXd::Zero(ncons, nparsfinal);
    
    VectorXd gradll = VectorXd::Zero(nparsfinal);
    
    for (unsigned int ires = 0; ires < dVs.size(); ++ires) {
      auto const &dVi = dVs[ires];
      const unsigned int residxi = residxsfinal[ires];
      
      SparseMatrix<double> &dViR = dVRs.emplace_back();
      dViR = dVi*Rsparse;
      
// const MatrixXd Rdense = R;
      
// const SparseMatrix<double> dViRtest = dVi*(Rdense.sparseView().transpose().transpose());
// const SparseMatrix<double> dViRtest = Rdense.sparseView().transpose();

// 
// MatrixXd tmpdense = MatrixXd::Zero(ncons, ncons);
// 
// constexpr unsigned int ntest = 1000;
// 
// auto time0 = std::chrono::high_resolution_clock::now();
// for (unsigned int itest = 0; itest < ntest; ++itest) {
// const SparseMatrix<double> dViRtest = (dVi*Rdense).sparseView();
// }
// auto time1 = std::chrono::high_resolution_clock::now();
// 
// for (unsigned int itest = 0; itest < ntest; ++itest) {
// const SparseMatrix<double> Rsparse = Rdense.sparseView();
// const SparseMatrix<double> dViRtest = dVi*Rsparse;
// }
// 
// auto time2 = std::chrono::high_resolution_clock::now();
// 
// for (unsigned int itest = 0; itest < ntest; ++itest) {
// const SparseMatrix<double> dViRtest = dVi*R;
// }
// 
// auto time3 = std::chrono::high_resolution_clock::now();
// 
// 
// for (unsigned int itest = 0; itest < ntest; ++itest) {
// tmpdense = dVi*Rdense;
// // const SparseMatrix<double> dViRtest = tmpdense.sparseView();
// }
// 
// auto time4 = std::chrono::high_resolution_clock::now();
// 
// const double d0 = (time1-time0).count();
// const double d1 = (time2-time1).count();
// const double d2 = (time3-time2).count();
// const double d3 = (time4-time3).count();
// 
// std::cout << "d0 = " << d0 << " d1 = " << d1 << " d2 = " << d2 << " d3 = " << d3 << std::endl;
      
// const double testfail = (dVi*MatrixXd::Identity(ncons, ncons)).sparseView();
      
// const int testcount = (dVi*MatrixXd(R)).nonZeros();
// const Matrix<double, 10, 10> testdense = Matrix<double, 10, 10>::Identity();
// const int testcount2 = testdense.nonZeros();
// 
// std::cout << "testcount2 = " << testcount2 << std::endl;
// 
// std::cout << "dViR total size = " << dViR.rows()*dViR.cols() << " testcount = " << testcount << std::endl;

      
// const MatrixXd ident = MatrixXd::Identity(ncons, ncons);
      
// SparseMatrix<double> testsparse = SparseMatrix<double>(dVi*ident);
      
// std::cout << "dViR total size = " << dViR.rows()*dViR.cols() << " nonzeros = " << dViR.nonZeros() << std::endl;
      
// dVRr.col(residxi) += dViR*rfull;
      dVRr.col(residxi) += dVi*Rr;
      
      //TODO optimize traces for sparse matrices
      
// const double gradres = MatrixXd(dViR).trace();
      
      //trace is not implemented directly for sparse matrices so use sum of diagonals explicitly
      const double gradres = dViR.diagonal().sum();
      
// const double gradres = -rfull.transpose()*R*dViR*rfull + MatrixXd(dViR).trace();
      
// grad(residxi) += gradres;
      gradll(residxi) += gradres;
            
// VectorXd hesscross = -2.*Jsparse.transpose()*R*dViR*rfull;
      // adjust diagonal element since it will be added twice
// hesscross(residxi) *= 0.5;
      
// hess.col(residxi) += hesscross;
// hess.row(residxi) += hesscross.transpose();
      
// if (!dogen) {
// dxdparms.row(residxi) += (M*dViR*rfull).transpose();
// }
      
      const SparseMatrix<double> dViRT = dViR.transpose();
      
      for (unsigned int jres = 0; jres <= ires; ++jres) {
        auto const &dVjR = dVRs[jres];
        const unsigned int residxj = residxsfinal[jres];
        
// const SparseMatrix<double> dViRdVjR = dViR*dVjR;
        
// const double hessres = -MatrixXd(dViR*dVjR).trace();
// const double hessres = 2.*rfull.transpose()*R*dViRdVjR*rfull - MatrixXd(dViRdVjR).trace();
        
        // below is equivalent to -(dViR*dVjR).trace()
// const double hessres = -dViRT.cwiseProduct(dVjR).sum();
        const double hessres = dViRT.cwiseProduct(dVjR).sum();
// const double hessresalt = -MatrixXd(dViR*dVjR).trace();
        
// std::cout << "ires = " << ires << " jres = " << jres << " hessres = " << hessres << " hessresalt = " << hessresalt << std::endl;
        
        hess(residxi, residxj) += hessres;
        if (ires != jres) {
          hess(residxj, residxi) += hessres;
        }
        
      }
      

    
      
      
// // const SparseMatrix<double> dViR = (dVi*R).sparseView();
// const SparseMatrix<double> dViR = dVi*R;
// // const MatrixXd dViR = dVi*R;
// 
// //TODO optimized trace for sparse matrices
// //TODO optimize use of sparse vs dense matrices (and avoid redundant conversions)
// 
// grad(residxi) += -rfull.transpose()*R*dViR*rfull + MatrixXd(dViR).trace();
// 
// VectorXd hesscross = -2.*Jsparse.transpose()*R*dViR*rfull;
// 
// // adjust diagonal element since it will be added twice
// hesscross(residxi) *= 0.5;
// 
// hess.col(residxi) += hesscross;
// hess.row(residxi) += hesscross.transpose();
// 
// if (!dogen) {
// dxdparms.row(residxi) += (M*dViR*rfull).transpose();
// }
// 
// for (unsigned int jres = 0; jres < dVs.size(); ++jres) {
// auto const &dVj = dVs[jres];
// const unsigned int residxj = residxsfinal[jres];
// 
// // const SparseMatrix<double> dViRdVjR = (dViR*dVj*R).sparseView();
// const SparseMatrix<double> dViRdVjR = dViR*dVj*R;
// // const MatrixXd dViRdVjR = dViR*dVj*R;
// 
// const double hessres = 2.*rfull.transpose()*R*dViRdVjR*rfull - MatrixXd(dViRdVjR).trace();
// 
// hess(residxi, residxj) += hessres;
// if (residxi != residxj) {
// hess(residxj, residxi) += hessres;
// }
// }
      
    }
    
    const SparseMatrix<double> dVRrsparse = dVRr.sparseView();
    
// std::cout << "dVRrsparse size = " << dVRrsparse.rows()*dVRrsparse.cols() << " nonzeros = " << dVRrsparse.nonZeros() << std::endl;
    
// grad += -(rfull.transpose()*R*dVRrsparse).transpose();
// grad += -(Rr.transpose()*dVRrsparse).transpose();
    grad += -dVRrsparse.transpose()*Rr;
    
    if (false) {
      hess += 2.*dVRrsparse.transpose()*R*dVRrsparse;
      
      const SparseMatrix<double> hesscross = -2.*Jsparse.transpose()*Rsparse*dVRrsparse;    
      hess += hesscross;
      hess += hesscross.transpose();
    }
    
// hess += hesscross + hesscross.transpose();

    //TODO deduplicate with above
    dxdparms(Eigen::placeholders::all, freestateidxs) += Cinvd.solve(FtVinv*dVRrsparse).transpose();
    
// if (!dogen) {
// //TODO deduplicate with above
// // const MatrixXd M = Cinvd.solve(FtVinv);
// // dxdparms += (M*dVRrsparse).transpose();
// dxdparms += Cinvd.solve(FtVinv*dVRrsparse).transpose();
// }
// 
    
    //TODO check this against element-wise version
// grad += -(rfull.transpose()*R*dVRrsparse).transpose();
// hess += 2.*dVRrsparse.transpose()*R*dVRrsparse;
// 
// MatrixXd hesstest = 2.*dVRrsparse.transpose()*R*dVRrsparse;
// MatrixXd hesstest2 = MatrixXd::Zero(nparsfinal, nparsfinal);
// for (unsigned int ires = 0; ires < dVs.size(); ++ires) {
// auto const &dVi = dVs[ires];
// const unsigned int residxi = residxsfinal[ires];
// 
// const SparseMatrix<double> &dViR = dVRs[ires];
// 
// for (unsigned int jres = 0; jres <= ires; ++jres) {
// auto const &dVjR = dVRs[jres];
// const unsigned int residxj = residxsfinal[jres];
// 
// const double hesstestval = 2.*rfull.transpose()*R*dViR*dVjR*rfull;
// // const double hesstestval = 2.*rfull.transpose()*dViR*dVjR*rfull;
// 
// // const double hesstestval2 = 2.*rfull.transpose()*dViR.transpose()*R*dVjR*rfull;
// 
// hesstest2(residxi, residxj) += hesstestval;
// // if (residxi != residxj) {
// if (ires != jres) {
// hesstest2(residxj, residxi) += hesstestval;
// }
// 
// // std::cout << "ires = " << ires << " jres = " << jres << " residxi = " << residxi << " residxj = " << residxj << " hesstestval = " << hesstestval << " hesstest(residxi, residxj) = " << hesstest(residxi, residxj) << std::endl;
// 
// }
// }
// 
// for (unsigned int ipar = 0; ipar < nparsfinal; ++ipar) {
// for (unsigned int jpar = 0; jpar < nparsfinal; ++jpar) {
// const double diff = hesstest2(ipar, jpar) - hesstest(ipar, jpar);
// 
// if (std::fabs(diff) > 1e-16 && std::fabs(diff)/hesstest(ipar, jpar) > 1e-3) {
// std::cout << "ipar = " << ipar << " jpar = " << jpar << " hesstest(ipar, jpar) = " << hesstest(ipar, jpar) << " hessest2(ipar, jpar) = " << hesstest2(ipar, jpar) << std::endl;
// }
// 
// }
// }
// 
// const double diffsq = (hesstest2-hesstest).array().square().sum();
// std::cout << "diffsq = " << diffsq << std::endl;
    
    
    
// MatrixXd hesscross = -2.*Jsparse.transpose()*R*dVRr;
// //correct diagonal since it gets added twice
// hesscross.diagonal() *= 0.5;
    
// hess += hesscross + hesscross.transpose();
    
// std::cout << "hesscross.bottomRightCorner<10,10>():\n" << hesscross.bottomRightCorner<10,10>() << std::endl;
    
// if (!dogen) {
// dxdparms += (M*dVRrsparse).transpose();
// }
    
    // The marginal objective the log-det gradient differentiates, in double.
    // Same definition and same code as the two-track maker's, so the two
    // trees' finite-difference checks are the same check.
    if (exportObjective_) {
      objchisq = rfull.dot(Rr);
      // ln|V| = -ln|Vinv| as a PSEUDO-determinant: Vinv is rank deficient by
      // construction (the deweighted strip coordinates carry exactly zero
      // weight).  The null space is structural, so it cancels in a finite
      // difference; `objnullv` records its size so that can be asserted.
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

    gradchisqv.clear();
    gradchisqv.resize(nparsfinal, 0.);

    Map<VectorXf>(gradchisqv.data(), nparsfinal) = grad.cast<float>();

    gradllv.clear();
    gradllv.resize(nparsfinal, 0.);
    Map<VectorXf>(gradllv.data(), nparsfinal) = gradll.cast<float>();
    
    grad += gradll;
    
    gradv.clear();
    jacrefv.clear();

    gradv.resize(nparsfinal,0.);
    jacrefv.resize(5*nparsfinal, 0.);
    
    nJacRef = 5*nparsfinal;
    if (fillTrackTree_ && fillGrads_) {
      tree->SetBranchAddress("gradv", gradv.data());
    }
    if (fillTrackTree_ && fillJac_) {
      tree->SetBranchAddress("jacrefv", jacrefv.data());
    }
    
    //eigen representation of the underlying vector storage
    Map<VectorXf> gradout(gradv.data(), nparsfinal);
    Map<Matrix<float, 5, Dynamic, RowMajor> > jacrefout(jacrefv.data(), 5, nparsfinal);
    
// jacrefout = dxdparms.leftCols<5>().transpose().cast<float>(); 
    jacrefout = ( (dxdparms).leftCols<5>().transpose() ).cast<float>();  


    if (false) {
      const double refpt = std::fabs(1./refParms[0])*std::sin(M_PI_2 - refParms[1]);
      const double refphi = refParms[2];
      const double reftheta = M_PI_2 - refParms[1];
      const double refeta = -std::log(std::tan(0.5*reftheta));

      std::cout << "ref pt eta phi: " << refpt << " " << refeta << " " << refphi << std::endl;
      std::cout << "jacref qop lam phi:\n" << jacrefout.topRows<3>() << std::endl;

    }

    

    gradout = grad.cast<float>();
    
    
    rx.resize(2*nvalid);
    Map<Matrix<float, Dynamic, 2, RowMajor> > rxout(rx.data(), nvalid, 2);
    rxout = rxfull;
// std::cout << "rx:" << std::endl;
// for (auto elem : rx) {
// std::cout << elem << " ";
// }
// std::cout << std::endl;
// std::cout << rx << std::endl;
    
    ry.resize(2*nvalid);
    Map<Matrix<float, Dynamic, 2, RowMajor> > ryout(ry.data(), nvalid, 2);
    ryout = ryfull;
    
// deigx.resize(nvalid);
// deigy.resize(nvalid);
// 
// validdxeig = validdxeigjac*dxfull;
// 
// for (unsigned int ivalid = 0; ivalid < nvalid; ++ivalid) {
// deigx[ivalid] = validdxeig[2*ivalid];
// deigy[ivalid] = validdxeig[2*ivalid + 1];
// }
    
    float refPt = dogen ? genpart->pt() : std::abs(1./refParms[0])*std::sin(M_PI_2 - refParms[1]);

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

    if (debugprintout_) {
      const Matrix5d Cinner = (Cinvd.solve(MatrixXd::Identity(nstateparms,nstateparms))).topLeftCorner<5,5>();
      std::cout << "hess debug" << std::endl;
      std::cout << "track parms" << std::endl;
      std::cout << tkparms << std::endl;
  // std::cout << "dxRef" << std::endl;
  // std::cout << dxRef << std::endl;
      std::cout << "original cov" << std::endl;
      std::cout << track.covariance() << std::endl;
      std::cout << "recomputed cov" << std::endl;
      std::cout << 2.*Cinner << std::endl;
    }

// std::cout << "dxinner/dparms" << std::endl;
// std::cout << dxdparms.bottomRows<5>() << std::endl;
// std::cout << "grad" << std::endl;
// std::cout << grad << std::endl;
// std::cout << "hess diagonal" << std::endl;
// std::cout << hess.diagonal() << std::endl;
// std::cout << "hess0 diagonal" << std::endl;
// std::cout << d2chisqdparms2.diagonal() << std::endl;
// std::cout << "hess1 diagonal" << std::endl;
// std::cout << 2.*(dxdparms.transpose()*d2chisqdxdparms).diagonal() << std::endl;
// std::cout << "hess2 diagonal" << std::endl;
// std::cout << (dxdparms.transpose()*d2chisqdx2*dxdparms).diagonal() << std::endl;
    
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

// std::cout << "refParms[0]: " << refParms[0] << std::endl;

    // Refit track, rebuilt from the converged reference state. Deliberately
    // NOT inside the muonref gate below: the bachelor-track use case has no
    // muon association at all, which is exactly the case the ValueMaps miss.
    // covfull's top-left 5x5 is the reference-state covariance in the same
    // (qoverp, lambda, phi, dxy, dz) basis reco::TrackBase uses.
    if (emitRefitTracks_) {
      const double qbp = refParmsMomD[0];
      const double lam = refParmsMomD[1];
      const double phi = refParmsMomD[2];
      const double pmag = 1. / std::max(std::abs(qbp), 1e-12);
      const math::XYZVector mom(pmag * std::cos(lam) * std::cos(phi),
                                pmag * std::cos(lam) * std::sin(phi),
                                pmag * std::sin(lam));
      const double dxyv = refParms[3];
      const double dzv = refParms[4];
      const math::XYZPoint refpt(-dxyv * std::sin(phi), dxyv * std::cos(phi), dzv);

      // Guard the covariance. Clamped / marginally-converged fits can leave a
      // non-finite or non-positive reference block, and a NaN covariance would
      // silently poison any downstream vertex fit. Such tracks are emitted as
      // the input copy with refitOk = 0 rather than as a corrected track.
      const Matrix<double, 5, 5> c5 = covfull.topLeftCorner<5, 5>();
      bool covok = c5.allFinite();
      for (int i = 0; i < 5 && covok; ++i)
        covok = c5(i, i) > 0.;

      if (covok) {
        reco::TrackBase::CovarianceMatrix cov;
        for (int i = 0; i < 5; ++i)
          for (int j = i; j < 5; ++j)
            cov(i, j) = c5(i, j);

        // NB: chisqval is the CVH fit's full chi2 (all hits plus the
        // scattering/material terms), and ndof here is only the nominal
        // tracker convention -- so normalizedChi2() on a refit track is NOT
        // a standard track-quality measure and must not be cut on as one.
        const int ndof = std::max(1, 2 * static_cast<int>(nvalid) - 5);
        const reco::Track cand(chisqval, ndof, refpt, mom,
                               static_cast<int>(std::copysign(1., qbp)),
                               cov, track.algo());

        // Optional relative-ptErr sanity cut (disabled when negative).
        const bool pterrok =
            refitMaxRelPtErr_ < 0. ||
            (cand.pt() > 0. && cand.ptError() / cand.pt() < refitMaxRelPtErr_);
        if (pterrok) {
          refitTracksV[itrack] = cand;
          refitOkV[itrack] = 1;
        }
      }
    }

    // Muon-association outputs (written per attempt; the winning attempt's
    // values persist since the fit is re-run to restore the nominal if the
    // opposite loses). Inert for cosmics/single-track (muonref null).
    if (muonref.isNonnull()) {
      const double qbpupd = refParmsMomD[0];
      const double lamupd = refParmsMomD[1];
      const double phiupd = refParmsMomD[2];

      const double ptupd = std::cos(lamupd)/std::abs(qbpupd);
      const double chargeupd = std::copysign(1., qbpupd);

      const double thetaupd = M_PI_2 - lamupd;
      const double etaupd = -std::log(std::tan(0.5*thetaupd));

      corPtV[muonref.index()] = ptupd;
      corEtaV[muonref.index()] = etaupd;
      corPhiV[muonref.index()] = phiupd;
      corChargeV[muonref.index()] = chargeupd;
      corDxyV[muonref.index()] = refParms[3];
      corDzV[muonref.index()] = refParms[4];
      edmvalV[muonref.index()] = edmval;
      nValidHitsV[muonref.index()] = nvalid;
      nValidPixelHitsV[muonref.index()] = nvalidpixel;

      auto &iglobalidxv = globalidxsV[muonref.key()];
      iglobalidxv.clear();
      iglobalidxv.reserve(globalidxvfinal.size());
      for (auto val : globalidxvfinal) {
        iglobalidxv.push_back(val);
      }

      auto &ijacrefv = jacRefV[muonref.key()];
      ijacrefv.assign(3*nparsfinal, 0.);
      //eigen representation of the underlying vector storage
      Map<Matrix<float, 3, Dynamic, RowMajor> > momjacrefout(ijacrefv.data(), 3, nparsfinal);
      momjacrefout = ( (dxdparms).leftCols<3>().transpose() ).cast<float>();

      auto &imomCov = momCovV[muonref.key()];
      imomCov.assign(3*3, 0.);
      Map<Matrix<float, 3, 3, RowMajor>> momCovout(imomCov.data(), 3, 3);
      momCovout.triangularView<Upper>() = covfull.topLeftCorner<3,3>().cast<float>();
    }

    return {true, static_cast<double>(chisqval)};
    };  // ---- end runHypothesis lambda --------------------------------------

    // Attempt 1: nominal seed charge.
    std::pair<bool, double> fitres = runHypothesis(seedChargeSign_);
    const unsigned int flipProtectNom = nChargeFlipProtect;
    chargeHypFlipped = 0u;
    // Attempt 2: only if the nominal fit converged but its clamp caught a q/p
    // sign crossing (it wanted the opposite charge -- the high-p mis-ID
    // signature). Keep whichever hypothesis has the lower chi2.
    if (twoHypothesisCharge_ && fitres.first && flipProtectNom > 0u) {
      const double chisqNom = fitres.second;
      const std::pair<bool, double> fitresOpp = runHypothesis(-seedChargeSign_);
      if (fitresOpp.first && fitresOpp.second < chisqNom) {
        // opposite charge fits better: member buffers already hold it
        fitres = fitresOpp;
        chargeHypFlipped = 1u;
        ++fitChargeHypTaken_;
      } else {
        // opposite worse or failed: re-run the nominal to restore its outputs
        fitres = runHypothesis(seedChargeSign_);
      }
    }
    // restore the diagnostic flag to the nominal attempt's value
    nChargeFlipProtect = flipProtectNom;

    if (!fitres.first) {
      continue;
    }
    ++fitSucceeded_;

    if (fillTrackTree_) {
      tree->Fill();
    }

  }

  edm::ValueMap<float> corPtMap;
  edm::ValueMap<float> corEtaMap;
  edm::ValueMap<float> corPhiMap;
  edm::ValueMap<int> corChargeMap;
  edm::ValueMap<float> corDxyMap;
  edm::ValueMap<float> corDzMap;
  edm::ValueMap<float> edmvalMap;
  edm::ValueMap<int> nValidHitsMap;
  edm::ValueMap<int> nValidPixelHitsMap;

  edm::ValueMap<std::vector<int>> globalidxsMap;

  edm::ValueMap<std::vector<float>> jacRefMap;
  edm::ValueMap<std::vector<float>> momCovMap;

  if (doMuonAssoc_) {
    edm::ValueMap<float>::Filler corPtMapFiller(corPtMap);
    corPtMapFiller.insert(muonAssoc->ref(), std::make_move_iterator(corPtV.begin()), std::make_move_iterator(corPtV.end()));
    corPtMapFiller.fill();

    edm::ValueMap<float>::Filler corEtaMapFiller(corEtaMap);
    corEtaMapFiller.insert(muonAssoc->ref(), std::make_move_iterator(corEtaV.begin()), std::make_move_iterator(corEtaV.end()));
    corEtaMapFiller.fill();

    edm::ValueMap<float>::Filler corPhiMapFiller(corPhiMap);
    corPhiMapFiller.insert(muonAssoc->ref(), std::make_move_iterator(corPhiV.begin()), std::make_move_iterator(corPhiV.end()));
    corPhiMapFiller.fill();

    edm::ValueMap<int>::Filler corChargeMapFiller(corChargeMap);
    corChargeMapFiller.insert(muonAssoc->ref(), std::make_move_iterator(corChargeV.begin()), std::make_move_iterator(corChargeV.end()));
    corChargeMapFiller.fill();

    edm::ValueMap<float>::Filler corDxyMapFiller(corDxyMap);
    corDxyMapFiller.insert(muonAssoc->ref(), std::make_move_iterator(corDxyV.begin()), std::make_move_iterator(corDxyV.end()));
    corDxyMapFiller.fill();

    edm::ValueMap<float>::Filler corDzMapFiller(corDzMap);
    corDzMapFiller.insert(muonAssoc->ref(), std::make_move_iterator(corDzV.begin()), std::make_move_iterator(corDzV.end()));
    corDzMapFiller.fill();

    edm::ValueMap<float>::Filler edmvalMapFiller(edmvalMap);
    edmvalMapFiller.insert(muonAssoc->ref(), std::make_move_iterator(edmvalV.begin()), std::make_move_iterator(edmvalV.end()));
    edmvalMapFiller.fill();

    edm::ValueMap<int>::Filler nValidHitsMapFiller(nValidHitsMap);
    nValidHitsMapFiller.insert(muonAssoc->ref(), std::make_move_iterator(nValidHitsV.begin()), std::make_move_iterator(nValidHitsV.end()));
    nValidHitsMapFiller.fill();

    edm::ValueMap<int>::Filler nValidPixelHitsMapFiller(nValidPixelHitsMap);
    nValidPixelHitsMapFiller.insert(muonAssoc->ref(), std::make_move_iterator(nValidPixelHitsV.begin()), std::make_move_iterator(nValidPixelHitsV.end()));
    nValidPixelHitsMapFiller.fill();

    edm::ValueMap<std::vector<int>>::Filler globalidxsMapFiller(globalidxsMap);
    globalidxsMapFiller.insert(muonAssoc->ref(), std::make_move_iterator(globalidxsV.begin()), std::make_move_iterator(globalidxsV.end()));
    globalidxsMapFiller.fill();

    edm::ValueMap<std::vector<float>>::Filler jacRefMapFiller(jacRefMap);
    jacRefMapFiller.insert(muonAssoc->ref(), std::make_move_iterator(jacRefV.begin()), std::make_move_iterator(jacRefV.end()));
    jacRefMapFiller.fill();

    edm::ValueMap<std::vector<float>>::Filler momCovMapFiller(momCovMap);
    momCovMapFiller.insert(muonAssoc->ref(), std::make_move_iterator(momCovV.begin()), std::make_move_iterator(momCovV.end()));
    momCovMapFiller.fill();
  }

  if (emitRefitTracks_) {
    auto refitOut = std::make_unique<reco::TrackCollection>(std::move(refitTracksV));
    const edm::OrphanHandle<reco::TrackCollection> refitH =
        iEvent.put(outputRefitTracks_, std::move(refitOut));
    // Keyed to the emitted collection, so refitOk[i] describes refit track i.
    auto okOut = std::make_unique<edm::ValueMap<int>>();
    edm::ValueMap<int>::Filler okFiller(*okOut);
    okFiller.insert(refitH, refitOkV.begin(), refitOkV.end());
    okFiller.fill();
    iEvent.put(outputRefitOk_, std::move(okOut));
  }

  iEvent.emplace(outputCorPt_, std::move(corPtMap));
  iEvent.emplace(outputCorEta_, std::move(corEtaMap));
  iEvent.emplace(outputCorPhi_, std::move(corPhiMap));
  iEvent.emplace(outputCorCharge_, std::move(corChargeMap));
  iEvent.emplace(outputCorDxy_, std::move(corDxyMap));
  iEvent.emplace(outputCorDz_, std::move(corDzMap));
  iEvent.emplace(outputEdmval_, std::move(edmvalMap));
  iEvent.emplace(outputNValidHits_, std::move(nValidHitsMap));
  iEvent.emplace(outputNValidPixelHits_, std::move(nValidPixelHitsMap));

  iEvent.emplace(outputGlobalIdxs_, std::move(globalidxsMap));

  iEvent.emplace(outputJacRef_, std::move(jacRefMap));
  iEvent.emplace(outputMomCov_, std::move(momCovMap));

}

DEFINE_FWK_MODULE(ResidualGlobalCorrectionMakerG4e);
