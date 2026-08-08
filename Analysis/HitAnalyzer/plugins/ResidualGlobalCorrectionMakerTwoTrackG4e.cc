#include "ResidualGlobalCorrectionMakerBase.h"
#include "Analysis/HitAnalyzer/interface/ParticleProperties.h"

// Sparse GBL design-matrix formulation (ported from the single-track
// maker). base.h provides Eigen/Core + Eigen/Eigenvalues + `using
// namespace Eigen`; the sparse solver path needs Eigen/Sparse too.
#include <Eigen/Sparse>

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

#include "Math/Vector4Dfwd.h"

#include "TrackPropagation/Geant4e/interface/Geant4ePropagator.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"

#include <iomanip>
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
                << "  skipped[samesign]=" << fitSkippedSameSign_
                << "  clamped[step]=" << fitStepClamped_
                << "  backtracked[step]=" << fitStepBacktracked_
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
  mutable unsigned long long fitStepBacktracked_ = 0ULL; // step halvings after a failed-leg iteration (retries)
  mutable unsigned long long fitSeedInflated_ = 0ULL;    // iteration-0 seed-momentum inflations (retries)
  // Momentum floor for the Gauss-Newton step clamp (GeV). 2 GeV suits
  // J/psi muons (as in the single-track maker); V0 drivers lower it to
  // sit above the propagation floor but below the soft-daughter spectrum.
  double clampMomentumFloor_ = 2.0;
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
  mutable unsigned long long fitSkippedSameSign_ = 0ULL; // same-sign pairs skipped pre-fit (not failures)

  // Global material model: per-leg per-group dxi columns from the
  // propagator (reused buffer; see doc/global-material-model-plan.md).
  mutable std::vector<std::pair<int, Eigen::Matrix<double, 5, 1>>> groupJacs_;
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
  
  float Jpsicons_d;
  float Jpsicons_x;
  float Jpsicons_y;
  float Jpsicons_z;
  float Jpsicons_pt;
  float Jpsicons_eta;
  float Jpsicons_phi;
  float Jpsicons_mass;
  
  float Mupluscons_pt;
  float Mupluscons_eta;
  float Mupluscons_phi;
  
  float Muminuscons_pt;
  float Muminuscons_eta;
  float Muminuscons_phi;
  
  float Jpsikincons_x;
  float Jpsikincons_y;
  float Jpsikincons_z;
  float Jpsikincons_pt;
  float Jpsikincons_eta;
  float Jpsikincons_phi;
  float Jpsikincons_mass;
  
  float Mupluskincons_pt;
  float Mupluskincons_eta;
  float Mupluskincons_phi;
  
  float Muminuskincons_pt;
  float Muminuskincons_eta;
  float Muminuskincons_phi;
  
  float Jpsigen_x;
  float Jpsigen_y;
  float Jpsigen_z;
  float Jpsigen_pt;
  float Jpsigen_eta;
  float Jpsigen_phi;
  float Jpsigen_mass;
  
  float Muplusgen_pt;
  float Muplusgen_eta;
  float Muplusgen_phi;
  
  float Muminusgen_pt;
  float Muminusgen_eta;
  float Muminusgen_phi;

  float Muplusgen_dr;
  float Muminusgen_dr;
  
  std::array<float, 3> Muplus_refParms;
  std::array<float, 3> Muminus_refParms;
  
  std::vector<float> Muplus_jacRef;
  std::vector<float> Muminus_jacRef;
  std::vector<float> Jpsi_jacMass;
  
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

};


ResidualGlobalCorrectionMakerTwoTrackG4e::ResidualGlobalCorrectionMakerTwoTrackG4e(
    const edm::ParameterSet &iConfig)
    : ResidualGlobalCorrectionMakerBase(iConfig),
      ttrhToken_(esConsumes(edm::ESInputTag("", "WithAngleAndTemplate"))),
      g4ePropToken_(esConsumes(edm::ESInputTag("", "Geant4ePropagator"))),
      transTrackBuilderToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
      l1MenuToken_(esConsumes())
{
  doVtxConstraint_ = iConfig.getParameter<bool>("doVtxConstraint");
  doMassConstraint_ = iConfig.getParameter<bool>("doMassConstraint");
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

    tree->Branch("Muplusgen_pt", &Muplusgen_pt);
    tree->Branch("Muplusgen_eta", &Muplusgen_eta);
    tree->Branch("Muplusgen_phi", &Muplusgen_phi);

    tree->Branch("Muminusgen_pt", &Muminusgen_pt);
    tree->Branch("Muminusgen_eta", &Muminusgen_eta);
    tree->Branch("Muminusgen_phi", &Muminusgen_phi);

    tree->Branch("Muplusgen_dr", &Muplusgen_dr);
    tree->Branch("Muminusgen_dr", &Muminusgen_dr);
    
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
 
  constexpr bool dolocalupdate = false;
  
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

  // doRes port (per-candidate mass-CF export): the material process-noise
  // derivative blocks dV_b (MS and ionization parts of Q), their row
  // ranges, and the per-block global-parameter labels. Registered fresh
  // each iteration of the unconstrained (icons==0) pass; the converged
  // iteration's content feeds the mass-projected influence export.
  const bool dores = doRes_;
  std::vector<SparseMatrix<double>> dVs;
  std::vector<std::array<unsigned int, 2>> resblockrng;
  std::vector<unsigned int> resglobidx;
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

      // Sparse GBL: the state params that are actually solved for. Mirrors
      // the dense freezeparm() logic (old lines ~2282-2298): fitFromGenParms_
      // freezes the 10-dim vertex PCA (idx 0..9); doVtxConstraint_ freezes
      // the track-PCA distance (idx 6). A frozen index is simply excluded
      // from freestateidxs instead of being deweighted with a 1e6 diagonal.
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

      for (unsigned int icons = 0; icons < nicons; ++icons) {

        // Sparse GBL constraint-row count for this icons pass:
        //  - 5 propagation/MS rows per hit (both tracks; nhits is the sum)
        //  - 2 measurement rows per valid hit (the two-track dense maker
        //    uses a uniform 2-dim Fhit/dy0 for every valid hit; strip hits
        //    are handled by a near-singular Vinv on the unmeasured coord,
        //    NOT by emitting fewer rows -- so it is 2*nvalid, not the
        //    single-track maker's nvalid+nvalidpixel)
        //  - 3 beamspot rows per track when bsConstraint_ (off by default)
        //  - 1 pointing row when doPointingConstraint_ (off by default)
        //  - 1 J/psi-mass row on the constrained pass (icons==1)
        const unsigned int nbscons = bsConstraint_ ? 3u * 2u : 0u;
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
        
        std::array<std::vector<Matrix<double, 7, 1>>, 2> layerStatesarr;
        for (unsigned int id = 0; id < 2; ++id) {
          auto const &hits = hitsarr[id];
          layerStatesarr[id].reserve(hits.size());
        }
        

        double chisqvalold = std::numeric_limits<double>::max();

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
            ioniurbanidx.clear();
            ioniurbanv.clear();
            msmoliidx.clear();
            msmoliv.clear();
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
            
            
            if (bsConstraint_) {
              // apply beamspot constraint
              // TODO add residual corrections for beamspot parameters?
              // TODO impose constraints on individual tracks when not applying common vertex constraint?
              
              constexpr unsigned int nlocalvtx = 3;
              constexpr unsigned int nlocal = nlocalvtx;
              constexpr unsigned int localvtxidx = 0;
              constexpr unsigned int fullvtxidx = 7;

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

              Matrix<double, 3, 1> dbs0;
              dbs0[0] = refFts[0] - x0;
              dbs0[1] = refFts[1] - y0;
              dbs0[2] = refFts[2] - z0;

              const Matrix<double, 3, nlocal> Fbs = Matrix<double, 3, 3>::Identity();
              const Matrix<double, 3, 3> covBSinv = covBS.inverse();

              const double bschisq = dbs0.transpose()*covBSinv*dbs0;
              chisq0val += bschisq;

              // Sparse GBL row write: 3 beamspot rows constraining the
              // vertex-PCA position (state idx 7,8,9). Emitted once per
              // track id (matching the dense path's per-id accumulation;
              // ncons accounts for 3*2). Fbs = Identity(3,3), residual
              // dbs0 = vertex - beamspot, weight covBSinv.
              rfull.segment<3>(irow) = dbs0;
              Ffull.block(irow, fullvtxidx, 3, nlocalvtx) =
                  Fbs.leftCols<nlocalvtx>();
              Vinvfull.block<3, 3>(irow, irow) = covBSinv;
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
                if (dores && icons == 0) {
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
                  }
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
                if (dores && icons == 0) {
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
                  }
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
              const double qopupd = qopref + dqop;
              double s = 1.;
              if (qopupd * qopref <= 0.) {
                // sign flip: stop half-way toward q/p = 0
                s = -0.5 * qopref / dqop;
              } else if (std::abs(qopupd) > 1. / clampMomentumFloor_) {
                // p_upd below the floor: land exactly on p = floor, same charge
                s = (std::copysign(1. / clampMomentumFloor_, qopref) - qopref) / dqop;
              }
              if (s < stepscale) {
                stepscale = s;
              }
            }
            if (stepscale < 1.) {
              stepscale = std::max(stepscale, 0.);
              dxfree *= stepscale;
              dxfull *= stepscale;
              if (!stepClampedThisFit) {
                stepClampedThisFit = true;
                ++fitStepClamped_;
              }
              std::cout << "GN step clamped (two-track): icons = " << icons
                        << " iiter = " << iiter << " scale = " << stepscale
                        << " seed0(q,pt,eta)=(" << itrack->charge() << "," << itrack->pt()
                        << "," << itrack->eta() << ")"
                        << " seed1(q,pt,eta)=(" << jtrack->charge() << "," << jtrack->pt()
                        << "," << jtrack->eta() << ")" << std::endl;
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
          
          chisqval = chisq0val + deltachisq;

          deltachisqval = chisq0val + deltachisq - chisqvalold;

          chisqvalold = chisq0val + deltachisq;
          
// ndof = 5*nhits + nvalid + nvalidalign2d - nstateparms;
          ndof = 5*nhits + nvalid + nvalidpixel - nstateparms;
          
          if (bsConstraint_) {
            ndof += 3;
          }

          if (doPointingConstraint_) {
            // 2D-transverse pointing adds 1 scalar constraint
            ++ndof;
          }

          if (doVtxConstraint_) {
            ++ndof;
          }
          
          if (icons == 1) {
            ++ndof;
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
            // define sign of d wrt charge of tracks
            Jpsi_d = firstplus ? statepcaupd[6] : -statepcaupd[6];
            Jpsi_x = statepcaupd[7];
            Jpsi_y = statepcaupd[8];
            Jpsi_z = statepcaupd[9];
          }
          else {
            // define sign of d wrt charge of tracks
            Jpsicons_d = firstplus ? statepcaupd[6] : -statepcaupd[6];
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
            reseigidx.clear();
            reseigv.clear();
            resinfbv.clear();
            if (dores && !dVs.empty()) {
              VectorXd afull = VectorXd::Zero(nstateparms);
              afull.head<6>() = mjacalt.transpose();
              VectorXd afree = VectorXd::Zero(nstatefree);
              for (unsigned int i = 0; i < nstatefree; ++i) {
                afree(i) = afull(freestateidxs[i]);
              }
              const VectorXd wmass = VinvF*Cinvd.solve(afree);
              for (unsigned int ires = 0; ires < dVs.size(); ++ires) {
                const unsigned int r0 = resblockrng[ires][0];
                const unsigned int nb = resblockrng[ires][1];
                const MatrixXd dVb = MatrixXd(dVs[ires]).block(r0, r0, nb, nb);
                SelfAdjointEigenSolver<MatrixXd> eigv(dVb);
                const MatrixXd sqrtdV = eigv.eigenvectors() *
                                        eigv.eigenvalues().cwiseMax(0.).cwiseSqrt().asDiagonal() *
                                        eigv.eigenvectors().transpose();
                const VectorXd ub = sqrtdV * wmass.segment(r0, nb);
                reseigidx.push_back(resglobidx[ires]);
                for (unsigned int j = 0; j < 5; ++j) {
                  resinfv.push_back(j < nb ? ub(j) : 0.f);
                }
                const double vb = ub.squaredNorm();
                resinfvarv.push_back(vb);
                resinfcov += vb;
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
          
          Muplusgen_dr = -99.;
          Muminusgen_dr = -99.;

          if (doGen_) {
            double drminplus = 0.1;
            double drminminus = 0.1;

            for (auto const &genpart : *genPartCollection) {
              if (genpart.status() != 1) {
                continue;
              }
              if (std::abs(genpart.pdgId()) != 13) {
                continue;
              }

// float dRplus = deltaR(genpart.phi(), muarr[idxplus].phi(), genpart.eta(), muarr[idxplus].eta());
              const double dRplus = deltaR(genpart, muarr[idxplus]);
              if (dRplus < drminplus && genpart.charge() > 0) {
                muplusgen = &genpart;
                drminplus = dRplus;
              }

// float dRminus = deltaR(genpart.phi(), muarr[idxminus].phi(), genpart.eta(), muarr[idxminus].eta());
              const double dRminus = deltaR(genpart, muarr[idxminus]);
              if (dRminus < drminminus && genpart.charge() < 0) {
                muminusgen = &genpart;
                drminminus = dRminus;
              }
            }

            if (muplusgen != nullptr) {
              Muplusgen_dr = drminplus;
            }

            if (muminusgen != nullptr) {
              Muminusgen_dr = drminminus;
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

          const Matrix<double, 10, 1> dxRef = dxfull.head<10>();
          const Matrix<double, 10, 10> covref = covstate.topLeftCorner<10, 10>();

          const Matrix<double, 10, 10> hessref = covref.inverse();

          const double deltachisqref = -0.5*dxRef.transpose()*hessref*dxRef;

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
        const SelfAdjointEigenSolver<MatrixXd> eshess(hess);
        const VectorXd& eigvals = eshess.eigenvalues();  // ascending

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

        const double lambdakeptmin = eigvals(nparsfinal - nrank);
        hessrankgap = (nrank < nparsfinal && lambdakeptmin > 0.)
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
