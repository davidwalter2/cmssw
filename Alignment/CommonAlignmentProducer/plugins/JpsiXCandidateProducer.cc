/// \class JpsiXCandidateProducer
///
/// EDProducer that builds nested reco::VertexCompositeCandidateCollection
/// objects for B-meson and quarkonium decays of the form mother -> J/psi + X,
/// where X is either a single bachelor track (B+, Bc channels) or a second
/// VertexCompositeCandidate representing an intermediate resonance (B0, Bs,
/// Lb, psi(2S) channels). Operates in two modes controlled by the xMode
/// parameter:
///
///   "track" mode (B+, Bc):
///     Iterates over generalTracks, skipping tracks with pT < minBachelorPt.
///     The bachelor TrackRef stored in the output RecoChargedCandidate points
///     directly into generalTracks so that the downstream chain
///     CompositeDaughterTrackProducer -> extra().key() -> remapper is intact.
///     A separate AlignmentTrackSelectorModule-based pre-filter MUST NOT be
///     used here because TrackCollectionStoreManager rewrites TrackExtraRefs,
///     breaking the remapper lookup.
///
///   "vcc" mode (B0->K*0, B0->Ks, Bs->phi, Lb->Lambda, psi(2S)->Ks):
///     Combines each J/psi candidate with each candidate from a second
///     VertexCompositeCandidateCollection (K*0, Ks, phi, Lambda).
///
/// In both modes the output candidate has:
///   daughter(0) = J/psi (cloned VertexCompositeCandidate, daughters intact)
///   daughter(1) = X (RecoChargedCandidate or nested VertexCompositeCandidate)
///
/// The J/psi four-momentum used in the mother mass calculation is the stored
/// candidate p4() -- raw track-based, consistent with
/// TwoBodyDecayCandidateProducer convention. The mass window is wide enough
/// to accommodate J/psi mass resolution spread without constraint.
///
/// Optional quality parameters (all backward-compatible, default = no cut):
///   minJpsiPt             early pT cut on J/psi before combination loop
///   minMotherPt           mother pT cut after mass window (both modes)
///   maxBachelorEta        |eta| cut on bachelor (track mode only)
///   maxJpsiAlphaBS        [deprecated] J/psi pointing angle to beamspot
///   maxBachelorIPToJpsiVertex  3D DCA of bachelor straight-line to J/psi
///                         vertex midpoint (track mode only)
///   maxBachelorMuTrackDOCA  per-(bachelor-leaf, J/psi-muon) track-track DOCA
///                         cut (cm). In track mode the bachelor is a single
///                         track; in vcc mode the producer loops over the
///                         VCC's leaf tracks and rejects the J/psi-X
///                         candidate if any (leaf, muon) pair exceeds the
///                         threshold. Straight-line approximation; matches
///                         Phase-1 bkmm_kaon_mu{1,2}_doca physics scale.
///
/// Phase 3 parameters (B-level Kalman vertex fit, all default = no cut):
///   maxMotherAlphaBS      B pointing angle (2D transverse cosalpha) using
///                         Kalman-fitted B vertex AND the closest-in-z primary
///                         vertex as origin. Despite the historical name,
///                         the cut is NOT beamspot-based - the beamspot has
///                         ~5-15 cm z-spread that dominated the original 3D
///                         angle and made it random vs B pointing. Requires
///                         offlinePrimaryVertices.
///   minBVtxProb           B Kalman vertex fit probability (chi2/ndof)
///   minBLxyOverSigma      B transverse flight significance from the closest-
///                         in-z offline PV (requires offlinePrimaryVertices)

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalError.h"
#include "TMath.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

// Helix propagation primitives used to evaluate each daughter's 4-momentum
// at a candidate-level common point (not at the track's PCA-to-beamline).
// Removes the ~0.1-0.5 mrad per-track angle bias at displaced vertices.
// See openspec/changes/finalize-jpsi-x-preset-b-production/design.md
// Section "Decision 3a".
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

// Kinematic fitter used to (i) apply a J/psi mass constraint to the dimuon
// when computing the B candidate 4-momentum under any preset (the constraint
// improves B mass resolution from ~30 MeV to ~15 MeV, essential to making
// real-data peaks visible), and (ii) under preset C, perform a multi-track
// vertex fit with the same J/psi mass constraint applied to the muon pair
// (replaces the old plain-Kalman path).
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"

#include <cmath>
#include <limits>
#include <set>
#include <string>
#include <utility>
#include <vector>

class JpsiXCandidateProducer : public edm::stream::EDProducer<> {
public:
  explicit JpsiXCandidateProducer(const edm::ParameterSet& cfg)
      : jpsiToken_(consumes<reco::VertexCompositeCandidateCollection>(
            cfg.getParameter<edm::InputTag>("jpsiSrc"))),
        trackMode_(cfg.getParameter<std::string>("xMode") == "track"),
        minMotherMass_(cfg.getParameter<double>("minMotherMass")),
        maxMotherMass_(cfg.getParameter<double>("maxMotherMass")),
        motherPdgId_(cfg.getParameter<int>("motherPdgId")),
        minJpsiPt_(cfg.existsAs<double>("minJpsiPt")
                       ? cfg.getParameter<double>("minJpsiPt") : 0.0),
        minMotherPt_(cfg.existsAs<double>("minMotherPt")
                         ? cfg.getParameter<double>("minMotherPt") : 0.0),
        maxJpsiAlphaBS_(cfg.existsAs<double>("maxJpsiAlphaBS")
                            ? cfg.getParameter<double>("maxJpsiAlphaBS")
                            : std::numeric_limits<double>::max()),
        applyAlphaBS_(cfg.existsAs<double>("maxJpsiAlphaBS")) {
    // Mother pdgId is used as the species tag on every emitted candidate
    // (VertexCompositeCandidate) and must be non-zero for the invariant that
    // downstream reads daughter->pdgId() as the mass-hypothesis source of
    // truth. openspec change add-jpsi-x-muons-and-preprod-refinements.
    if (motherPdgId_ == 0) {
      throw cms::Exception("Configuration")
          << "JpsiXCandidateProducer requires motherPdgId != 0; "
          << "downstream consumers read daughter.pdgId() as the species tag.";
    }
    if (trackMode_) {
      trackToken_ = consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("trackSrc"));
      minBachelorPt_ = cfg.getParameter<double>("minBachelorPt");
      bachelorMass_ = cfg.getParameter<double>("bachelorMass");
      bachelorPdgId_ = cfg.getParameter<int>("bachelorPdgId");
      if (bachelorPdgId_ == 0) {
        throw cms::Exception("Configuration")
            << "JpsiXCandidateProducer (track mode) requires bachelorPdgId != 0.";
      }
      maxBachelorEta_ = cfg.existsAs<double>("maxBachelorEta")
                            ? cfg.getParameter<double>("maxBachelorEta")
                            : std::numeric_limits<double>::max();
      maxBachelorIPToJpsiVertex_ = cfg.existsAs<double>("maxBachelorIPToJpsiVertex")
                                       ? cfg.getParameter<double>("maxBachelorIPToJpsiVertex")
                                       : std::numeric_limits<double>::max();
    }
    // Phase-1-tuned: track-track DOCA between each bachelor leaf track and each
    // J/psi-daughter muon track. Available to both track-mode (bachelor =
    // generalTracks track) and vcc-mode (bachelor = VCC, looped over daughters).
    maxBachelorMuTrackDOCA_ = cfg.existsAs<double>("maxBachelorMuTrackDOCA")
                                  ? cfg.getParameter<double>("maxBachelorMuTrackDOCA")
                                  : std::numeric_limits<double>::max();
    if (!trackMode_) {
      xToken_ = consumes<reco::VertexCompositeCandidateCollection>(
          cfg.getParameter<edm::InputTag>("xSrc"));
    }

    // J/psi mass constraint on the dimuon (default off). Under preset B we
    // want raw track-sum dimuon kinematics so that the alignment fit and any
    // Stage-2 CVH refit downstream see unconstrained inputs. Under preset C
    // the multi-track Kalman fit applies its own TwoTrackMassKinematicConstraint
    // internally; this flag, when true, also activates the dimuon-only
    // KinematicParticleFitter as the fallback dimuon p4 if the multi-track
    // fit fails. See openspec/changes/finalize-jpsi-x-preset-b-production.
    applyJpsiMassConstraint_ = cfg.existsAs<bool>("applyJpsiMassConstraint")
                                   ? cfg.getParameter<bool>("applyJpsiMassConstraint")
                                   : false;

    // Beamspot: needed by deprecated J/psi-level alphaBS or new B-level alphaBS.
    bool needBS = applyAlphaBS_ || cfg.existsAs<double>("maxMotherAlphaBS");
    if (needBS) {
      bsToken_ = consumes<reco::BeamSpot>(
          cfg.existsAs<edm::InputTag>("beamSpotSrc")
              ? cfg.getParameter<edm::InputTag>("beamSpotSrc")
              : edm::InputTag("offlineBeamSpot"));
    }

    // Phase 3: B-level Kalman vertex fit parameters.
    if (cfg.existsAs<double>("maxMotherAlphaBS")) {
      maxMotherAlphaBS_ = cfg.getParameter<double>("maxMotherAlphaBS");
      applyMotherAlphaBS_ = true;
    }
    if (cfg.existsAs<double>("minBVtxProb"))
      minBVtxProb_ = cfg.getParameter<double>("minBVtxProb");
    if (cfg.existsAs<double>("minBLxyOverSigma")) {
      minBLxyOverSigma_ = cfg.getParameter<double>("minBLxyOverSigma");
      applyBLxyOverSigma_ = (minBLxyOverSigma_ > 0.0);
    }
    applyBVtxFit_ = (minBVtxProb_ > 0.0 || applyBLxyOverSigma_ || applyMotherAlphaBS_);
    // PV collection needed for both the Lxy/sigma cut AND the (now PV-based)
    // cosalpha cut. Either of those being active triggers the consume.
    if (applyBLxyOverSigma_ || applyMotherAlphaBS_) {
      pvToken_ = consumes<reco::VertexCollection>(
          cfg.existsAs<edm::InputTag>("primaryVerticesSrc")
              ? cfg.getParameter<edm::InputTag>("primaryVerticesSrc")
              : edm::InputTag("offlinePrimaryVertices"));
    }

    produces<reco::VertexCompositeCandidateCollection>();
  }

  ~JpsiXCandidateProducer() override {
    edm::LogInfo("JpsiXCandidateProducer")
        << "Producer summary: "
        << "J/psi-constraint attempts=" << n_jpsi_constraint_attempted_
        << " fallbacks=" << n_jpsi_constraint_fallback_
        << "; B-vertex-fit attempts=" << n_b_vertex_fit_attempted_
        << " fallbacks=" << n_b_vertex_fit_fallback_
        << "; pair-prop attempts=" << n_pair_propagation_attempted_
        << " fallbacks=" << n_pair_propagation_fallback_
        << "; single-prop attempts=" << n_single_propagation_attempted_
        << " fallbacks=" << n_single_propagation_fallback_;
  }

  void produce(edm::Event& evt, const edm::EventSetup& iSetup) override {
    edm::Handle<reco::VertexCompositeCandidateCollection> jpsiH;
    evt.getByToken(jpsiToken_, jpsiH);

    const reco::BeamSpot* bs = nullptr;
    if (applyAlphaBS_ || applyMotherAlphaBS_) {
      edm::Handle<reco::BeamSpot> bsH;
      evt.getByToken(bsToken_, bsH);
      bs = bsH.product();
    }

    const reco::VertexCollection* pvs = nullptr;
    if (applyBLxyOverSigma_ || applyMotherAlphaBS_) {
      edm::Handle<reco::VertexCollection> pvH;
      evt.getByToken(pvToken_, pvH);
      pvs = pvH.product();
    }

    // TransientTrackBuilder is always needed: the J/psi mass constraint runs
    // under every preset (it's a mass-resolution correctness fix, not a
    // preset-selectable feature), and the constrained-vertex fit under
    // preset C needs the same builder.
    const TransientTrackBuilder* ttb = nullptr;
    edm::ESHandle<TransientTrackBuilder> ttbH;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttbH);
    ttb = ttbH.product();

    auto out = std::make_unique<reco::VertexCompositeCandidateCollection>();

    if (trackMode_) {
      edm::Handle<reco::TrackCollection> trackH;
      evt.getByToken(trackToken_, trackH);
      produceTrackMode(*jpsiH, *trackH, trackH, bs, pvs, ttb, *out);
    } else {
      edm::Handle<reco::VertexCompositeCandidateCollection> xH;
      evt.getByToken(xToken_, xH);
      produceVccMode(*jpsiH, *xH, bs, pvs, ttb, *out);
    }

    LogDebug("JpsiXCandidateProducer") << "Built " << out->size() << " candidates.";
    evt.put(std::move(out));
  }

private:
  // alphaBS using the J/psi candidate vertex (deprecated: for B->J/psiX the
  // correct observable is B-level alphaBS, see computeAlphaBSFromPoint below).
  static double computeAlphaBS(const reco::VertexCompositeCandidate& jpsi,
                                const reco::BeamSpot& bs) {
    const double vx = jpsi.vx() - bs.x(jpsi.vz());
    const double vy = jpsi.vy() - bs.y(jpsi.vz());
    const double vz = jpsi.vz() - bs.z0();
    const double r = std::sqrt(vx * vx + vy * vy + vz * vz);
    if (r < 1e-10) return 0.0;
    const double cosAlpha =
        (jpsi.px() * vx + jpsi.py() * vy + jpsi.pz() * vz) / (jpsi.p() * r);
    return std::acos(std::min(1.0, std::max(-1.0, cosAlpha)));
  }

  // B-level pointing angle, 2D transverse, with the matched PV as origin.
  // This is the correct observable: the transverse plane is dominated by the
  // real B flight (cT*beta*gamma projected on transverse, ~500 um for B+),
  // not by the PV-vs-bunch-center z offset that dominated the previous
  // 3D-from-beamspot version (see fix-jpsi-x-cosalpha-bug-and-jpsi-mass-
  // constraint for the bug history).
  static double computeAlphaPV(double vx, double vy,
                                double px, double py,
                                const reco::Vertex& pv) {
    const double dvx = vx - pv.x();
    const double dvy = vy - pv.y();
    const double rT = std::sqrt(dvx * dvx + dvy * dvy);
    if (rT < 1e-10) return 0.0;
    const double pT = std::sqrt(px * px + py * py);
    if (pT < 1e-10) return 0.0;
    const double cosAlpha = (px * dvx + py * dvy) / (pT * rT);
    return std::acos(std::min(1.0, std::max(-1.0, cosAlpha)));
  }

  // Pick the PV closest in z to the candidate's B vertex; fall back to
  // pvs.front() (the highest-sumPt^2 PV) if no PV is within 1 cm in z.
  // Standard B-physics convention for matching B candidates to their
  // production vertex.
  static const reco::Vertex& bestPVForCandidate(double vz_B,
                                                  const reco::VertexCollection& pvs) {
    size_t best_idx = 0;
    double best_dz = std::abs(vz_B - pvs[0].z());
    for (size_t i = 1; i < pvs.size(); ++i) {
      const double dz = std::abs(vz_B - pvs[i].z());
      if (dz < best_dz) { best_dz = dz; best_idx = i; }
    }
    if (best_dz > 1.0) {
      edm::LogWarning("JpsiXCandidateProducer")
          << "No PV within 1 cm of B vertex z=" << vz_B
          << " (closest dz=" << best_dz << "); falling back to pvs.front()";
      return pvs.front();
    }
    return pvs[best_idx];
  }

  // Straight-line approximation to the 3D DCA between two tracks. Matches the
  // Phase-1 quantity bkmm_kaon_mu{1,2}_doca used to tune the cut: same physics
  // scale (cm), straight-line within ~10% of helix-based DOCA for the short
  // distances probed by the cut (~0.03 cm). For parallel tracks falls back to
  // the point-on-second-track distance.
  static double trackTrackDCA(const reco::Track& t1, const reco::Track& t2) {
    const double inv1 = 1.0 / t1.p();
    const double inv2 = 1.0 / t2.p();
    const double u1x = t1.px() * inv1, u1y = t1.py() * inv1, u1z = t1.pz() * inv1;
    const double u2x = t2.px() * inv2, u2y = t2.py() * inv2, u2z = t2.pz() * inv2;
    const double cx = u1y * u2z - u1z * u2y;
    const double cy = u1z * u2x - u1x * u2z;
    const double cz = u1x * u2y - u1y * u2x;
    const double cN = std::sqrt(cx * cx + cy * cy + cz * cz);
    if (cN < 1e-12) {
      return trackDCAToPoint(t1, t2.vx(), t2.vy(), t2.vz());
    }
    const double dx = t2.vx() - t1.vx();
    const double dy = t2.vy() - t1.vy();
    const double dz = t2.vz() - t1.vz();
    return std::fabs(dx * cx + dy * cy + dz * cz) / cN;
  }

  // 3D DCA of a track (straight-line approximation) to a given 3D point.
  static double trackDCAToPoint(const reco::Track& tr, double px, double py, double pz) {
    const double invP = 1.0 / tr.p();
    const double ux = tr.px() * invP;
    const double uy = tr.py() * invP;
    const double uz = tr.pz() * invP;
    const double dx = px - tr.vx();
    const double dy = py - tr.vy();
    const double dz = pz - tr.vz();
    const double cx = dy * uz - dz * uy;
    const double cy = dz * ux - dx * uz;
    const double cz = dx * uy - dy * ux;
    return std::sqrt(cx * cx + cy * cy + cz * cz);
  }

  // Signed pdgId for a charged daughter given the (positive) species code.
  // Charged leptons have pdgId sign opposite to their charge (mu- = +13);
  // hadrons carry the sign of their charge (K+ = +321, pi+ = +211).
  static int signedPdgId(int species, int charge) {
    const int a = std::abs(species);
    const bool lepton = (a == 11 || a == 13 || a == 15);
    return lepton ? -charge * a : charge * a;
  }

  // Recursively collect TrackRefs of all leaf RecoChargedCandidates.
  static void collectLeafTrackRefs(const reco::Candidate& cand,
                                    std::vector<reco::TrackRef>& refs) {
    if (cand.numberOfDaughters() == 0) {
      const auto* rcc = dynamic_cast<const reco::RecoChargedCandidate*>(&cand);
      if (rcc) {
        const auto& ref = rcc->track();
        if (ref.isNonnull()) refs.push_back(ref);
      }
    } else {
      for (size_t i = 0; i < cand.numberOfDaughters(); ++i)
        collectLeafTrackRefs(*cand.daughter(i), refs);
    }
  }

  // Recursively collect track keys for self-combination avoidance.
  void collectLeafTrackKeys(const reco::Candidate& cand,
                             const edm::ProductID& srcId,
                             std::set<size_t>& keys) const {
    if (cand.numberOfDaughters() == 0) {
      const auto* rcc = dynamic_cast<const reco::RecoChargedCandidate*>(&cand);
      if (rcc) {
        const auto& ref = rcc->track();
        if (ref.isNonnull() && ref.id() == srcId) keys.insert(ref.key());
      }
    } else {
      for (size_t i = 0; i < cand.numberOfDaughters(); ++i)
        collectLeafTrackKeys(*cand.daughter(i), srcId, keys);
    }
  }

  // Apply a J/psi mass constraint to the dimuon and return the constrained
  // J/psi 4-momentum. Used under any preset that doesn't run the full multi-
  // track constrained vertex fit (i.e. preset B). On fit failure returns
  // false; the caller falls back to the unconstrained J/psi 4-momentum and
  // increments the fallback counter.
  bool constrainJpsi4Momentum(const std::vector<reco::TrackRef>& muRefs,
                                 const TransientTrackBuilder& ttb,
                                 reco::Particle::LorentzVector& lvJpsiOut) const {
    if (muRefs.size() != 2 || !muRefs[0].isNonnull() || !muRefs[1].isNonnull())
      return false;

    // Local non-const copies: CMSSW's KinematicParticleFactory and
    // MassKinematicConstraint take their mass/sigma args by non-const ref.
    ParticleMass muonMass = kMuonMass_;
    float muonMassSigma   = kMuonMassSigma_;
    ParticleMass jpsiMass = kJpsiMass_;
    float jpsiMassSigma   = kJpsiMassSigma_;

    KinematicParticleFactoryFromTransientTrack pFactory;
    std::vector<RefCountedKinematicParticle> dimuon;
    dimuon.reserve(2);
    for (const auto& muRef : muRefs) {
      reco::TransientTrack tt = ttb.build(muRef);
      if (!tt.isValid()) return false;
      float chi2 = 0.0f, ndf = 0.0f;
      dimuon.push_back(pFactory.particle(tt, muonMass, chi2, ndf, muonMassSigma));
    }

    KinematicParticleVertexFitter pvFitter;
    RefCountedKinematicTree tree;
    try {
      tree = pvFitter.fit(dimuon);
    } catch (...) {
      return false;
    }
    if (!tree.get() || tree->isEmpty()) return false;

    KinematicParticleFitter cFitter;
    auto jpsiConstraint = std::unique_ptr<KinematicConstraint>(
        new MassKinematicConstraint(jpsiMass, jpsiMassSigma));
    try {
      tree = cFitter.fit(jpsiConstraint.get(), tree);
    } catch (...) {
      return false;
    }
    if (!tree.get() || tree->isEmpty()) return false;

    tree->movePointerToTheTop();
    RefCountedKinematicParticle jpsiP = tree->currentParticle();
    if (!jpsiP || !jpsiP->currentState().isValid()) return false;

    auto state = jpsiP->currentState();
    GlobalVector p3 = state.globalMomentum();
    double m = state.mass();
    double E = std::sqrt(p3.mag2() + m * m);
    lvJpsiOut = reco::Particle::LorentzVector(p3.x(), p3.y(), p3.z(), E);
    return true;
  }

  // ------------------------------------------------------------------
  // Helix propagation of daughter tracks to a candidate-level common
  // point. Runs under EVERY preset (no preset gating).
  //
  // Why: tracks are 5-parameter helices reported at their PCA-to-beamline.
  // For displaced decays (B+ ctau ~ 491 um, similar for Bc/B0/Bs/Lb) each
  // track's momentum *direction* differs by ~0.1-0.5 mrad between the
  // PCA-to-beamline and the actual decay vertex (rotation = arc/rho in a
  // B field). Summing daughter 4-vectors at different points along
  // different helices biases the invariant mass by ~1-3 MeV per candidate,
  // scaling linearly with the mother ctau. Resolution impact is invisible
  // (<<30 MeV detector floor) but the kinematics-correlated bias matters
  // for the alignment use case where Stage-2 CVH ingests raw preset-B
  // tracks. See openspec/changes/finalize-jpsi-x-preset-b-production.
  //
  // Implementation uses two existing CMSSW classes:
  //   propagatePair        -> ClosestApproachInRPhi on two helices
  //                            (gives both PCA points + GlobalTrajectoryParameters
  //                             at the PCAs + crossing point in one call)
  //   propagateSingleToPoint -> AnalyticalImpactPointExtrapolator on one
  //                              helix to a given GlobalPoint.

  // Two-helix propagation: returns true on success. On success, lvA_out /
  // lvB_out are the 4-vectors built from each daughter's propagated
  // (px, py, pz) and species mass; crossing_out is the helix-helix
  // crossing point used as the pair's vertex estimate.
  bool propagatePair(const reco::TrackRef& a,
                     const reco::TrackRef& b,
                     double mA, double mB,
                     const TransientTrackBuilder& ttb,
                     reco::Particle::LorentzVector& lvA_out,
                     reco::Particle::LorentzVector& lvB_out,
                     GlobalPoint& crossing_out) const {
    if (!a.isNonnull() || !b.isNonnull()) return false;
    reco::TransientTrack tta = ttb.build(a);
    reco::TransientTrack ttb_ = ttb.build(b);
    if (!tta.isValid() || !ttb_.isValid()) return false;
    FreeTrajectoryState ftsA = tta.initialFreeState();
    FreeTrajectoryState ftsB = ttb_.initialFreeState();

    ClosestApproachInRPhi cap;
    bool ok = false;
    try {
      ok = cap.calculate(ftsA, ftsB);
    } catch (...) {
      return false;
    }
    if (!ok || !cap.status()) return false;

    auto pars = cap.trajectoryParameters();
    GlobalVector pA = pars.first.momentum();
    GlobalVector pB = pars.second.momentum();
    double EA = std::sqrt(pA.mag2() + mA * mA);
    double EB = std::sqrt(pB.mag2() + mB * mB);
    lvA_out = reco::Particle::LorentzVector(pA.x(), pA.y(), pA.z(), EA);
    lvB_out = reco::Particle::LorentzVector(pB.x(), pB.y(), pB.z(), EB);
    crossing_out = cap.crossingPoint();
    return true;
  }

  // Single-track propagation to a fixed GlobalPoint via
  // AnalyticalImpactPointExtrapolator. Used for the bachelor in track mode,
  // propagated to the dimuon crossing point.
  bool propagateSingleToPoint(const reco::TrackRef& t,
                              double m,
                              const GlobalPoint& target,
                              const TransientTrackBuilder& ttb,
                              reco::Particle::LorentzVector& lv_out) const {
    if (!t.isNonnull()) return false;
    reco::TransientTrack tt = ttb.build(t);
    if (!tt.isValid()) return false;
    FreeTrajectoryState fts = tt.initialFreeState();

    AnalyticalImpactPointExtrapolator extrap(ttb.field());
    TrajectoryStateOnSurface tsos;
    try {
      tsos = extrap.extrapolate(fts, target);
    } catch (...) {
      return false;
    }
    if (!tsos.isValid()) return false;

    GlobalVector p = tsos.globalMomentum();
    double E = std::sqrt(p.mag2() + m * m);
    lv_out = reco::Particle::LorentzVector(p.x(), p.y(), p.z(), E);
    return true;
  }

  // Result of the multi-track constrained vertex fit (preset C path).
  struct ConstrainedBFitResult {
    reco::Particle::LorentzVector lvM;
    GlobalPoint vertex;
    GlobalError vertexError;
    double chi2;
    int ndof;
  };

  // Build kinematic particles from (track ref, mass-hypothesis) pairs.
  bool buildKinematicParticles(
      const std::vector<std::pair<reco::TrackRef, double>>& refsAndMass,
      const TransientTrackBuilder& ttb,
      std::vector<RefCountedKinematicParticle>& particlesOut) const {
    KinematicParticleFactoryFromTransientTrack pFactory;
    particlesOut.reserve(refsAndMass.size());
    for (const auto& [ref, mass] : refsAndMass) {
      if (!ref.isNonnull()) return false;
      reco::TransientTrack tt = ttb.build(ref);
      if (!tt.isValid()) return false;
      float chi2 = 0.0f, ndf = 0.0f;
      // Local non-const copies: factory takes non-const refs for these.
      ParticleMass particleMass = mass;
      float massSigma = (std::fabs(mass - static_cast<double>(kMuonMass_)) < 1e-6)
                            ? kMuonMassSigma_ : 1.0e-4f;
      particlesOut.push_back(pFactory.particle(tt, particleMass, chi2, ndf, massSigma));
    }
    return true;
  }

  // Full multi-track constrained vertex fit. The first two entries in
  // refsAndMass MUST be the two muons (used by TwoTrackMassKinematicConstraint
  // which constrains the invariant mass of the first two particles).
  bool constrainedBVertexFit(
      const std::vector<std::pair<reco::TrackRef, double>>& refsAndMass,
      const TransientTrackBuilder& ttb,
      ConstrainedBFitResult& resultOut) const {
    if (refsAndMass.size() < 3) return false;
    std::vector<RefCountedKinematicParticle> particles;
    if (!buildKinematicParticles(refsAndMass, ttb, particles)) return false;

    // Local non-const copy: TwoTrackMassKinematicConstraint takes
    // ParticleMass& (non-const reference) by API design.
    ParticleMass jpsiMass = kJpsiMass_;
    KinematicConstrainedVertexFitter kcvFitter;
    auto jpsiC = std::unique_ptr<MultiTrackKinematicConstraint>(
        new TwoTrackMassKinematicConstraint(jpsiMass));
    RefCountedKinematicTree fitTree;
    try {
      fitTree = kcvFitter.fit(particles, jpsiC.get());
    } catch (...) {
      return false;
    }
    if (!fitTree.get() || fitTree->isEmpty()) return false;

    fitTree->movePointerToTheTop();
    RefCountedKinematicParticle bCand = fitTree->currentParticle();
    RefCountedKinematicVertex bVertex = fitTree->currentDecayVertex();
    if (!bCand || !bCand->currentState().isValid()) return false;
    if (!bVertex || !bVertex->vertexIsValid()) return false;

    auto state = bCand->currentState();
    GlobalVector p3 = state.globalMomentum();
    double m = state.mass();
    double E = std::sqrt(p3.mag2() + m * m);
    resultOut.lvM = reco::Particle::LorentzVector(p3.x(), p3.y(), p3.z(), E);
    resultOut.vertex = bVertex->position();
    resultOut.vertexError = bVertex->error();
    resultOut.chi2 = bVertex->chiSquared();
    resultOut.ndof = static_cast<int>(bVertex->degreesOfFreedom());
    return true;
  }

  // Evaluate the three preset-C cuts on a constrained-vertex-fit result.
  // Returns false if the candidate should be dropped.
  bool passesBVertexCuts(const ConstrainedBFitResult& r,
                          const reco::VertexCollection* pvs) const {
    if (minBVtxProb_ > 0.0) {
      double prob = TMath::Prob(r.chi2, r.ndof);
      if (prob < minBVtxProb_) return false;
    }
    const reco::Vertex* matchedPV = nullptr;
    if ((applyBLxyOverSigma_ || applyMotherAlphaBS_) && pvs && !pvs->empty()) {
      matchedPV = &bestPVForCandidate(r.vertex.z(), *pvs);
    }
    if (applyBLxyOverSigma_ && matchedPV) {
      double dx = r.vertex.x() - matchedPV->position().x();
      double dy = r.vertex.y() - matchedPV->position().y();
      double lxy = std::sqrt(dx * dx + dy * dy);
      if (lxy < 1e-10) return false;
      const reco::Vertex::Error& pvErr = matchedPV->error();
      double ux = dx / lxy, uy = dy / lxy;
      double sig2 = ux * ux * (r.vertexError.cxx() + pvErr(0, 0))
                  + 2 * ux * uy * (r.vertexError.cyx() + pvErr(0, 1))
                  + uy * uy * (r.vertexError.cyy() + pvErr(1, 1));
      if (lxy / std::sqrt(std::max(sig2, 1e-20)) < minBLxyOverSigma_) return false;
    }
    if (applyMotherAlphaBS_ && matchedPV) {
      if (computeAlphaPV(r.vertex.x(), r.vertex.y(),
                         r.lvM.px(), r.lvM.py(), *matchedPV)
          > maxMotherAlphaBS_)
        return false;
    }
    return true;
  }

  void produceTrackMode(const reco::VertexCompositeCandidateCollection& jpsiCands,
                        const reco::TrackCollection& tracks,
                        const edm::Handle<reco::TrackCollection>& trackH,
                        const reco::BeamSpot* bs,
                        const reco::VertexCollection* pvs,
                        const TransientTrackBuilder* ttb,
                        reco::VertexCompositeCandidateCollection& out) {
    for (const auto& jpsi : jpsiCands) {
      if (jpsi.pt() < minJpsiPt_) continue;
      if (applyAlphaBS_ && computeAlphaBS(jpsi, *bs) > maxJpsiAlphaBS_) continue;

      // Collect J/psi daughter track keys to avoid self-combinations.
      std::set<size_t> jpsiDaughterKeys;
      collectLeafTrackKeys(jpsi, trackH.id(), jpsiDaughterKeys);

      // Collect the J/psi muon refs once per J/psi (used by mass constraint
      // AND, under preset C, the multi-track constrained vertex fit).
      std::vector<reco::TrackRef> muRefs;
      collectLeafTrackRefs(jpsi, muRefs);

      // Helix propagation of the dimuon to the dimuon helix-helix PCA.
      // Computed once per J/psi (vJpsiPCA reused as the propagation target
      // for every bachelor in the inner loop). Falls back to raw jpsi.p4()
      // if the propagation fails.
      reco::Particle::LorentzVector lvJpsi = jpsi.p4();
      GlobalPoint vJpsiPCA(jpsi.vx(), jpsi.vy(), jpsi.vz());
      bool propPair = false;
      if (ttb && muRefs.size() == 2) {
        ++n_pair_propagation_attempted_;
        reco::Particle::LorentzVector lvMu0, lvMu1;
        if (propagatePair(muRefs[0], muRefs[1], kMuonMass_, kMuonMass_,
                          *ttb, lvMu0, lvMu1, vJpsiPCA)) {
          lvJpsi = lvMu0 + lvMu1;
          propPair = true;
        } else {
          ++n_pair_propagation_fallback_;
        }
      }

      // Optional J/psi mass constraint on top of the propagated dimuon.
      // Off by default (preset B); on under preset C, where it is
      // structurally required by the multi-track Kalman fit and used as
      // the fallback dimuon p4 when the Kalman fit fails.
      if (applyJpsiMassConstraint_ && ttb && muRefs.size() == 2) {
        ++n_jpsi_constraint_attempted_;
        reco::Particle::LorentzVector lvJpsiConstrained;
        if (constrainJpsi4Momentum(muRefs, *ttb, lvJpsiConstrained))
          lvJpsi = lvJpsiConstrained;
        else
          ++n_jpsi_constraint_fallback_;
      }

      for (size_t iT = 0; iT < tracks.size(); ++iT) {
        const reco::Track& tr = tracks[iT];
        if (tr.pt() < minBachelorPt_) continue;
        if (std::fabs(tr.eta()) > maxBachelorEta_) continue;
        if (jpsiDaughterKeys.count(iT)) continue;
        if (trackDCAToPoint(tr, jpsi.vx(), jpsi.vy(), jpsi.vz()) > maxBachelorIPToJpsiVertex_)
          continue;

        // Bachelor track to each J/psi-daughter muon track-track DOCA cut.
        if (maxBachelorMuTrackDOCA_ < std::numeric_limits<double>::max()) {
          bool docaFail = false;
          for (const auto& muRef : muRefs) {
            if (!muRef.isNonnull()) continue;
            if (trackTrackDCA(tr, *muRef) > maxBachelorMuTrackDOCA_) { docaFail = true; break; }
          }
          if (docaFail) continue;
        }

        double eBach = std::sqrt(tr.p() * tr.p() + bachelorMass_ * bachelorMass_);
        reco::Particle::LorentzVector lvBach(tr.px(), tr.py(), tr.pz(), eBach);
        // Propagate the bachelor to the dimuon crossing point (when the
        // dimuon propagation succeeded; otherwise we have no good target,
        // and the raw bachelor 4-vector at its native PCA is the best we
        // can do).
        if (propPair && ttb) {
          ++n_single_propagation_attempted_;
          reco::Particle::LorentzVector lvBachProp;
          if (propagateSingleToPoint(reco::TrackRef(trackH, iT), bachelorMass_,
                                     vJpsiPCA, *ttb, lvBachProp))
            lvBach = lvBachProp;
          else
            ++n_single_propagation_fallback_;
        }
        reco::Particle::LorentzVector lvM = lvJpsi + lvBach;
        // Note: lvM/mass below are computed from the J/psi-constrained 4-mom
        // (preset B path) or from the multi-track-constrained fit (preset C,
        // overwritten further down) - both are tighter than the raw sum.
        const double mass = lvM.M();
        if (mass < minMotherMass_ || mass > maxMotherMass_) continue;
        if (lvM.Pt() < minMotherPt_) continue;

        reco::Particle::Point vtxM(
            0.5 * (tr.vx() + jpsi.vx()),
            0.5 * (tr.vy() + jpsi.vy()),
            0.5 * (tr.vz() + jpsi.vz()));

        if (applyBVtxFit_) {
          // Multi-track constrained vertex fit. First two entries MUST be
          // the muons (TwoTrackMassKinematicConstraint constrains the first
          // two particles' invariant mass).
          ++n_b_vertex_fit_attempted_;
          std::vector<std::pair<reco::TrackRef, double>> refsAndMass;
          refsAndMass.reserve(3);
          for (const auto& muRef : muRefs)
            refsAndMass.emplace_back(muRef, kMuonMass_);
          refsAndMass.emplace_back(reco::TrackRef(trackH, iT), bachelorMass_);
          ConstrainedBFitResult fr;
          if (constrainedBVertexFit(refsAndMass, *ttb, fr)) {
            if (!passesBVertexCuts(fr, pvs)) continue;
            lvM = fr.lvM;
            vtxM = reco::Particle::Point(fr.vertex.x(), fr.vertex.y(), fr.vertex.z());
          } else {
            // Fit failed - keep candidate with J/psi-constrained sum + midpoint
            // vertex; preset-C cuts cannot be evaluated, so don't reject just
            // because the fit failed. Counter logged at end of job.
            ++n_b_vertex_fit_fallback_;
          }
        }

        int pdgBach = signedPdgId(bachelorPdgId_, tr.charge());
        reco::RecoChargedCandidate dBach(tr.charge(), lvBach,
                                         {tr.vx(), tr.vy(), tr.vz()}, pdgBach);
        dBach.setTrack(reco::TrackRef(trackH, iT));

        const int chargeM = jpsi.charge() + tr.charge();
        reco::VertexCompositeCandidate cand(chargeM, lvM, vtxM, motherPdgId_);
        cand.addDaughter(jpsi);   // daughter(0): J/psi with its daughters intact
        cand.addDaughter(dBach);  // daughter(1): bachelor RecoChargedCandidate
        out.push_back(std::move(cand));
      }
    }
  }

  void produceVccMode(const reco::VertexCompositeCandidateCollection& jpsiCands,
                      const reco::VertexCompositeCandidateCollection& xCands,
                      const reco::BeamSpot* bs,
                      const reco::VertexCollection* pvs,
                      const TransientTrackBuilder* ttb,
                      reco::VertexCompositeCandidateCollection& out) {
    for (const auto& jpsi : jpsiCands) {
      if (jpsi.pt() < minJpsiPt_) continue;
      if (applyAlphaBS_ && computeAlphaBS(jpsi, *bs) > maxJpsiAlphaBS_) continue;

      // Collect J/psi muon refs once per J/psi (used by the dimuon helix
      // propagation, the J/psi mass constraint if enabled, and the
      // multi-track Kalman fit under preset C).
      std::vector<reco::TrackRef> muRefs;
      collectLeafTrackRefs(jpsi, muRefs);

      // Helix propagation of the dimuon to its helix-helix PCA. The
      // crossing point is unused in vcc mode (the X-side gets its own
      // pair propagation below), so we only need the propagated lvJpsi.
      reco::Particle::LorentzVector lvJpsi = jpsi.p4();
      if (ttb && muRefs.size() == 2) {
        ++n_pair_propagation_attempted_;
        reco::Particle::LorentzVector lvMu0, lvMu1;
        GlobalPoint vJpsiPCA_unused;
        if (propagatePair(muRefs[0], muRefs[1], kMuonMass_, kMuonMass_,
                          *ttb, lvMu0, lvMu1, vJpsiPCA_unused)) {
          lvJpsi = lvMu0 + lvMu1;
        } else {
          ++n_pair_propagation_fallback_;
        }
      }

      // Optional J/psi mass constraint (preset C only).
      if (applyJpsiMassConstraint_ && ttb && muRefs.size() == 2) {
        ++n_jpsi_constraint_attempted_;
        reco::Particle::LorentzVector lvJpsiConstrained;
        if (constrainJpsi4Momentum(muRefs, *ttb, lvJpsiConstrained))
          lvJpsi = lvJpsiConstrained;
        else
          ++n_jpsi_constraint_fallback_;
      }

      for (const auto& xCand : xCands) {
        // Helix propagation of the X sub-resonance to its helix-helix PCA.
        // Each daughter's mass hypothesis is read from its stored
        // RecoChargedCandidate::mass() (set upstream by
        // TwoBodyDecayCandidateProducer for K*0/phi, V0Producer for Ks/Lambda).
        std::vector<reco::TrackRef> xDaughterRefs;
        collectLeafTrackRefs(xCand, xDaughterRefs);

        // Veto X candidates that reuse one of the J/psi muon tracks. The
        // K*0/phi/pipi producers run on generalTracks with no muon veto, so
        // a J/psi muon can be re-paired into the X candidate; combining it
        // with its own J/psi would double-count the track (and, under
        // preset C, feed the same track twice to the multi-track fit).
        // Mirrors the jpsiDaughterKeys veto in track mode.
        bool sharesTrackWithJpsi = false;
        for (const auto& xRef : xDaughterRefs) {
          for (const auto& muRef : muRefs) {
            if (xRef.id() == muRef.id() && xRef.key() == muRef.key()) {
              sharesTrackWithJpsi = true;
              break;
            }
          }
          if (sharesTrackWithJpsi) break;
        }
        if (sharesTrackWithJpsi) continue;

        reco::Particle::LorentzVector lvX = xCand.p4();
        if (ttb && xDaughterRefs.size() == 2 &&
            xCand.numberOfDaughters() == 2) {
          ++n_pair_propagation_attempted_;
          const auto* rcc0 = dynamic_cast<const reco::RecoChargedCandidate*>(xCand.daughter(0));
          const auto* rcc1 = dynamic_cast<const reco::RecoChargedCandidate*>(xCand.daughter(1));
          if (rcc0 && rcc1) {
            reco::Particle::LorentzVector lvD0, lvD1;
            GlobalPoint vXPCA;
            if (propagatePair(xDaughterRefs[0], xDaughterRefs[1],
                              rcc0->mass(), rcc1->mass(),
                              *ttb, lvD0, lvD1, vXPCA)) {
              lvX = lvD0 + lvD1;
            } else {
              ++n_pair_propagation_fallback_;
            }
          } else {
            ++n_pair_propagation_fallback_;
          }
        }

        reco::Particle::LorentzVector lvM = lvJpsi + lvX;
        const double mass = lvM.M();
        if (mass < minMotherMass_ || mass > maxMotherMass_) continue;
        if (lvM.Pt() < minMotherPt_) continue;

        // VCC daughter-track to each J/psi-daughter muon track-track DOCA cut.
        // Loop over the VCC's leaf tracks (K, pi for K*0; K, K for phi; pi, pi
        // for Ks; p, pi for Lambda) and require the DOCA against each muon
        // below the configured threshold. Reject the J/psi-X candidate if any
        // (daughter, muon) pair exceeds the cut.
        if (maxBachelorMuTrackDOCA_ < std::numeric_limits<double>::max()) {
          bool docaFail = false;
          for (const auto& muRef : muRefs) {
            if (!muRef.isNonnull()) continue;
            for (const auto& xRef : xDaughterRefs) {
              if (!xRef.isNonnull()) continue;
              if (trackTrackDCA(*muRef, *xRef) > maxBachelorMuTrackDOCA_) {
                docaFail = true; break;
              }
            }
            if (docaFail) break;
          }
          if (docaFail) continue;
        }

        reco::Particle::Point vtxM(
            0.5 * (jpsi.vx() + xCand.vx()),
            0.5 * (jpsi.vy() + xCand.vy()),
            0.5 * (jpsi.vz() + xCand.vz()));

        if (applyBVtxFit_) {
          // Multi-track constrained vertex fit over (mu, mu, x1, x2).
          // First two MUST be muons. The X-side daughter masses come from
          // each VCC daughter's stored .mass() (built by TwoBodyDecayCandidate-
          // Producer with the channel-specific PDG mass hypothesis).
          ++n_b_vertex_fit_attempted_;
          std::vector<std::pair<reco::TrackRef, double>> refsAndMass;
          refsAndMass.reserve(2 + xDaughterRefs.size());
          for (const auto& muRef : muRefs)
            refsAndMass.emplace_back(muRef, kMuonMass_);
          for (size_t iD = 0; iD < xCand.numberOfDaughters(); ++iD) {
            const auto* rcc = dynamic_cast<const reco::RecoChargedCandidate*>(xCand.daughter(iD));
            if (!rcc || !rcc->track().isNonnull()) {
              refsAndMass.clear();
              break;
            }
            refsAndMass.emplace_back(rcc->track(), rcc->mass());
          }
          if (refsAndMass.size() >= 4) {
            ConstrainedBFitResult fr;
            if (constrainedBVertexFit(refsAndMass, *ttb, fr)) {
              if (!passesBVertexCuts(fr, pvs)) continue;
              lvM = fr.lvM;
              vtxM = reco::Particle::Point(fr.vertex.x(), fr.vertex.y(), fr.vertex.z());
            } else {
              ++n_b_vertex_fit_fallback_;
            }
          } else {
            ++n_b_vertex_fit_fallback_;
          }
        }

        const int chargeM = jpsi.charge() + xCand.charge();
        reco::VertexCompositeCandidate cand(chargeM, lvM, vtxM, motherPdgId_);
        cand.addDaughter(jpsi);   // daughter(0): J/psi
        cand.addDaughter(xCand);  // daughter(1): intermediate resonance
        out.push_back(std::move(cand));
      }
    }
  }

  edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> jpsiToken_;
  edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> xToken_;
  edm::EDGetTokenT<reco::BeamSpot> bsToken_;
  edm::EDGetTokenT<reco::VertexCollection> pvToken_;

  bool trackMode_;
  double minMotherMass_;
  double maxMotherMass_;
  int motherPdgId_;

  double minJpsiPt_ = 0.;
  double minMotherPt_ = 0.;
  double maxJpsiAlphaBS_ = std::numeric_limits<double>::max();
  bool applyAlphaBS_;

  double minBachelorPt_ = 0.;
  double bachelorMass_ = 0.;
  int bachelorPdgId_ = 0;
  double maxBachelorEta_ = std::numeric_limits<double>::max();
  double maxBachelorIPToJpsiVertex_ = std::numeric_limits<double>::max();
  double maxBachelorMuTrackDOCA_ = std::numeric_limits<double>::max();

  // J/psi mass constraint on the dimuon: false under preset B (default),
  // true under preset C (where the multi-track Kalman fit needs it).
  bool applyJpsiMassConstraint_ = false;

  // Phase 3: B-level Kalman vertex fit
  double maxMotherAlphaBS_ = std::numeric_limits<double>::max();
  bool applyMotherAlphaBS_ = false;
  double minBVtxProb_ = 0.0;
  double minBLxyOverSigma_ = 0.0;
  bool applyBVtxFit_ = false;
  bool applyBLxyOverSigma_ = false;

  // J/psi mass constraint (applied to dimuon under any preset).
  // Constants as static-const float because ParticleMass is a typedef of float.
  static constexpr float kMuonMass_      = 0.10565837f;
  static constexpr float kMuonMassSigma_ = 1.0e-6f;
  static constexpr float kJpsiMass_      = 3.0969f;
  static constexpr float kJpsiMassSigma_ = 1.0e-6f;

  // Fallback counters for the end-of-job summary log.
  mutable std::size_t n_jpsi_constraint_attempted_ = 0;
  mutable std::size_t n_jpsi_constraint_fallback_  = 0;
  mutable std::size_t n_b_vertex_fit_attempted_    = 0;
  mutable std::size_t n_b_vertex_fit_fallback_     = 0;
  mutable std::size_t n_pair_propagation_attempted_   = 0;
  mutable std::size_t n_pair_propagation_fallback_    = 0;
  mutable std::size_t n_single_propagation_attempted_ = 0;
  mutable std::size_t n_single_propagation_fallback_  = 0;
};

DEFINE_FWK_MODULE(JpsiXCandidateProducer);
