// JpsiXKinematicFitProducer
//
// Builds a properly FITTED mother candidate from a stage-1 J/psi + X
// candidate, using the RecoVertex/KinematicFit machinery that Bmm5 uses.
//
// Why this exists: the stage-1 candidate's p4() is a raw four-vector SUM
// ("raw track-based", by design, so the CVH refit sees unconstrained inputs),
// and the two-track CVH maker's corMass is the refit DIMUON mass -- not a
// mother mass. So nothing upstream produces a fitted m(mu mu X).
//
// DUAL ARM (add-cvh-refit-tracks-into-bfit). The fit is run TWICE per candidate
// with an identical configuration, differing only in the input tracks:
//   RAW   arm ("rawFit*", "rawDimuon*")  -- the stage-1 tracks.
//   REFIT arm ("refFit*", "refDimuon*")  -- the CVH single-track refit of each
//                                           leg, falling back to the raw track
//                                           for any leg with refitOk = 0.
// Both arms apply the J/psi mass + common-vertex constraints exactly ONCE, in
// the fit (`inFit`, the Bmm5 convention), so the two are a clean A/B: the only
// difference is the per-track parameters and covariances. The driver aliases the
// unsuffixed physics columns (cvhFit*, dimuon*) to the REFIT arm.
//
// Refit-track mapping is LEAF-KEYED. Each single-track maker emits a `refit`
// TrackCollection (1:1 with its input tracks) + a `refitOk` map; its input-track
// producer (CandidateLeafTrackProducer) emits parallel `candIdx` / `leafIdx`
// vectors. So refit track j belongs to candidate candIdx[j], leaf leafIdx[j] --
// and we assemble, per candidate, the refit track for each leaf position in the
// SAME canonical leaf ordering `collectLeaves` produces here. A leg with
// refitOk = 0 (or no refit at all) uses its raw track; `nLegsRefit` counts how
// many legs were substituted.
//
// Modes (`jpsiConstraint`), kept for A/B on the constraint placement:
//   "inFit"    -- Bmm5: apply the J/psi mass constraint inside the mother fit.
//   "upstream" -- muons already J/psi-constrained upstream; no second constraint.
//   "cascade"  -- fit+mass-constrain the dimuon, collapse it, then fit vs the
//                 bachelor with the vertex floating.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"

#include "TMath.h"

namespace ana_hitanalyzer {

namespace {
constexpr double kMuonMass = 0.1056583745;
constexpr double kMuonSigma = 1.6e-9;
constexpr double kJpsiMass = 3.0969;

// PDG mass for a leaf daughter, from the signed pdgId stage-1 stamped on it.
double massForPdgId(int pdgid) {
  switch (std::abs(pdgid)) {
    case 13: return kMuonMass;
    case 211: return 0.13957039;   // pi
    case 321: return 0.493677;     // K
    case 2212: return 0.9382720813;  // p
    default: return kMuonMass;
  }
}

// Leaf tracks of a candidate node, depth-first. Mirrors the decomposition the
// CVH maker (and CandidateLeafTrackProducer's leafIdx) use: a composite daughter
// recurses, a leaf yields its track.
void collectLeaves(const reco::Candidate* c,
                   std::vector<const reco::RecoChargedCandidate*>& out) {
  if (c == nullptr) return;
  const auto* rcc = dynamic_cast<const reco::RecoChargedCandidate*>(c);
  if (rcc != nullptr) {
    out.push_back(rcc);
    return;
  }
  for (size_t i = 0; i < c->numberOfDaughters(); ++i)
    collectLeaves(c->daughter(i), out);
}

// A gen particle whose PDG id carries a b (anti)quark: for mesons the b sits
// in the hundreds digit (B+ 521, B0 511, Bs 531, Bc 541), for baryons in the
// thousands digit (Lambda_b 5122). Used to gen-match the reco B candidate.
inline bool isBHadron(int pdgId) {
  const int p = std::abs(pdgId);
  return (p / 100) % 10 == 5 || (p / 1000) % 10 == 5;
}

// ---------------------------------------------------------------------------
// Per-leg generator matching and the common-ancestor search.
//
// Ported from Bmm5/NanoAOD, which solved this for the same channel:
//   dr_match             <- GenBmmProducer.cc:74
//   is_acceptable        <- src/CommonTools.cc:105
//   get_mother           <- src/CommonTools.cc:117
//   find_common_ancestor <- src/CommonTools.cc:132
//
// Two details are load-bearing and are the reason this is ported rather than
// written fresh:
//
//  * the ancestor search permutes the generation DEPTH assigned to each leg.
//    In B+ -> J/psi(-> mu mu) K+ the muons sit two generations below the B and
//    the kaon one, so any fixed-depth comparison never finds the B.
//  * is_acceptable skips quarks, gluons and protons. Without it the search
//    always succeeds -- every pair of particles in the event shares the beam
//    proton as an ancestor -- and the whole truth test silently becomes
//    vacuous while continuing to look like it works.
// ---------------------------------------------------------------------------

// A match needs BOTH angular and momentum agreement. dR alone picks the wrong
// track routinely inside a b jet.
inline bool genDrMatch(const reco::Candidate& reco, const reco::Candidate& gen,
                       double maxDR, double maxRelDPt) {
  if (gen.pt() <= 0.) return false;
  return std::abs(reco.pt() - gen.pt()) / gen.pt() < maxRelDPt &&
         reco::deltaR(reco, gen) < maxDR;
}

inline bool isAcceptableAncestor(const reco::Candidate* c) {
  if (c == nullptr) return false;
  const int p = std::abs(c->pdgId());
  if (p < 10) return false;     // quarks and leptons-as-partons
  if (p == 2212) return false;  // beam protons
  if (p == 21) return false;    // gluons
  return true;
}

// depth 0 = the immediate mother
inline const reco::Candidate* genMotherAtDepth(const reco::Candidate* c, unsigned depth) {
  if (c == nullptr) return nullptr;
  const reco::Candidate* m = c->mother();
  unsigned i = 0;
  while (isAcceptableAncestor(m) && i < depth) {
    ++i;
    m = m->mother();
  }
  return isAcceptableAncestor(m) ? m : nullptr;
}

inline std::vector<unsigned> depthsFromPermutation(const std::vector<unsigned>& elements) {
  std::vector<unsigned> out;
  unsigned counter = 0;
  for (unsigned e : elements) {
    if (e == 0) ++counter;
    else { out.push_back(counter); counter = 0; }
  }
  out.push_back(counter);
  return out;
}

const reco::Candidate* findCommonAncestor(const std::vector<const reco::Candidate*>& parts,
                                          unsigned maxDepth = 10) {
  const size_t n = parts.size();
  if (n == 0) return nullptr;
  for (unsigned depth = 0; depth < maxDepth; ++depth) {
    std::vector<unsigned> elements;
    for (unsigned i = 0; i < depth; ++i) elements.push_back(0);
    for (size_t i = 0; i + 1 < n; ++i) elements.push_back(1);
    do {
      const auto depths = depthsFromPermutation(elements);
      const reco::Candidate* common = nullptr;
      for (size_t i = 0; i < n; ++i) {
        const reco::Candidate* m = genMotherAtDepth(parts[i], depths[i]);
        if (m == nullptr) { common = nullptr; break; }
        if (common == nullptr) common = m;
        if (common != m) { common = nullptr; break; }
      }
      if (common != nullptr) return common;
    } while (std::next_permutation(elements.begin(), elements.end()));
  }
  return nullptr;
}

// One refit "leg" source: a single-track maker's refit collection + refitOk map,
// keyed back to (candidate, leaf) via the input-track producer's index vectors.
struct RefitLegTokens {
  edm::EDGetTokenT<reco::TrackCollection> tracks;
  edm::EDGetTokenT<edm::ValueMap<int>> refitOk;
  edm::EDGetTokenT<std::vector<int>> candIdx;
  edm::EDGetTokenT<std::vector<int>> leafIdx;
};

// Per-candidate fit result for one arm.
struct ArmOut {
  float mass = -99.f, massErr = -99.f, pt = -99.f, eta = -99.f, phi = -99.f;
  float vchi2 = -99.f, vndof = -99.f, vprob = -99.f;
  int ok = 0;
  float mmVtxProb = -99.f, mmAlphaBS = -99.f, mmSl3d = -99.f, mmSl3dPV = -99.f;
  // UNCONSTRAINED dimuon mass from the same fit (gap G11). The dimuon fit
  // below applies no mass constraint, so this is a genuine free J/psi mass
  // -- the natural cross-check on the muon momentum scale. It was being
  // computed and thrown away; the only J/psi mass reaching the NanoAOD was
  // the mother fit's CONSTRAINED parameter, pinned to the PDG value.
  float mmMass = -99.f, mmMassErr = -99.f;
  // Vertex block, consumed by CandidateVertexGeometryProducer. The same
  // fixed set every fitter emits, so CandidateVertexGeometryProducer can be
  // instantiated per arm instead of each fitter deriving its own geometry.
  // motherPt/Eta/Phi are the FITTED momentum -- the candidate's own p4() is
  // the raw four-vector sum and must not be used for a pointing angle.
  float vtxX = -99.f, vtxY = -99.f, vtxZ = -99.f;
  float cXX = -99.f, cXY = -99.f, cXZ = -99.f, cYY = -99.f, cYZ = -99.f, cZZ = -99.f;
};
// Per-event storage for one arm's vertex block. A struct rather than another
// dozen positional arguments on store()/putArm().
struct BlockVecs {
  std::vector<float> x, y, z, cxx, cxy, cxz, cyy, cyz, czz, pt, eta, phi;
  explicit BlockVecs(std::size_t n)
      : x(n, -99.f), y(n, -99.f), z(n, -99.f), cxx(n, -99.f), cxy(n, -99.f),
        cxz(n, -99.f), cyy(n, -99.f), cyz(n, -99.f), czz(n, -99.f),
        pt(n, -99.f), eta(n, -99.f), phi(n, -99.f) {}
  void store(std::size_t i, const ArmOut& a) {
    x[i] = a.vtxX; y[i] = a.vtxY; z[i] = a.vtxZ;
    cxx[i] = a.cXX; cxy[i] = a.cXY; cxz[i] = a.cXZ;
    cyy[i] = a.cYY; cyz[i] = a.cYZ; czz[i] = a.cZZ;
    pt[i] = a.pt; eta[i] = a.eta; phi[i] = a.phi;
  }
};
}  // namespace

class JpsiXKinematicFitProducer : public edm::stream::EDProducer<> {
public:
  explicit JpsiXKinematicFitProducer(const edm::ParameterSet&);
  ~JpsiXKinematicFitProducer() override = default;

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  const edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> candToken_;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbToken_;
  std::vector<RefitLegTokens> refitLegs_;
  edm::EDGetTokenT<reco::BeamSpot> bsToken_;
  edm::EDGetTokenT<reco::VertexCollection> pvToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle>> genToken_;
  const bool useRefitTracks_;
  const bool useBeamSpot_;
  const bool usePV_;
  const bool haveGen_;
  const std::string mode_;
  const double jpsiMass_;
  const double maxChi2_;

  // Both arms emit the same set of fit + dimuon handles, prefixed raw*/ref*.
  struct ArmTokens {
    edm::EDPutTokenT<edm::ValueMap<float>> mass, massErr, pt, eta, phi, vchi2, vndof, vprob;
    edm::EDPutTokenT<edm::ValueMap<int>> ok;
    edm::EDPutTokenT<edm::ValueMap<float>> mmVtxProb, mmAlphaBS, mmSl3d, mmSl3dPV;
    edm::EDPutTokenT<edm::ValueMap<float>> mmMass, mmMassErr;
    // vertex block
    edm::EDPutTokenT<edm::ValueMap<float>> vtxX, vtxY, vtxZ,
        cXX, cXY, cXZ, cYY, cYZ, cZZ, mPt, mEta, mPhi;
  } raw_, ref_;
  edm::EDPutTokenT<edm::ValueMap<int>> outNLegsRefit_;
  // Generator matching (MC only; sentinels on data). Arm-independent.
  edm::EDPutTokenT<edm::ValueMap<float>> outGenBMass_, outGenBPt_, outGenBEta_,
      outGenBPhi_, outGenBDR_;
  edm::EDPutTokenT<edm::ValueMap<int>> outGenBPdgId_, outGenBIdx_;
  // Per-leg generator match (leg 0/1 = muons, leg 2 = bachelor) and the common
  // ancestor of the matched legs. The ancestor's species IS the truth category.
  std::array<edm::EDPutTokenT<edm::ValueMap<int>>, 3> outLegGenIdx_, outLegGenPdgId_,
      outLegGenMotherPdgId_;
  std::array<edm::EDPutTokenT<edm::ValueMap<float>>, 3> outLegGenPt_, outLegGenDR_;
  edm::EDPutTokenT<edm::ValueMap<int>> outAncPdgId_, outAncIdx_, outNLegsMatched_;
  edm::EDPutTokenT<edm::ValueMap<float>> outAncPt_, outAncMass_;
  const double genLegMaxDR_, genLegMaxRelDPt_;
  const unsigned genAncestorMaxDepth_;

  ArmTokens declareArm(const std::string& prefix);
};

JpsiXKinematicFitProducer::ArmTokens
JpsiXKinematicFitProducer::declareArm(const std::string& p) {
  ArmTokens t;
  t.mass = produces<edm::ValueMap<float>>(p + "FitMass");
  t.massErr = produces<edm::ValueMap<float>>(p + "FitMassErr");
  t.pt = produces<edm::ValueMap<float>>(p + "FitPt");
  t.eta = produces<edm::ValueMap<float>>(p + "FitEta");
  t.phi = produces<edm::ValueMap<float>>(p + "FitPhi");
  t.vchi2 = produces<edm::ValueMap<float>>(p + "FitVtxChi2");
  t.vndof = produces<edm::ValueMap<float>>(p + "FitVtxNdof");
  t.vprob = produces<edm::ValueMap<float>>(p + "FitVtxProb");
  t.ok = produces<edm::ValueMap<int>>(p + "FitOk");
  t.mmVtxProb = produces<edm::ValueMap<float>>(p + "DimuonVtxProb");
  t.mmAlphaBS = produces<edm::ValueMap<float>>(p + "DimuonAlphaBS");
  t.mmSl3d = produces<edm::ValueMap<float>>(p + "DimuonSxy");
  t.mmSl3dPV = produces<edm::ValueMap<float>>(p + "DimuonSl3d");
  t.mmMass = produces<edm::ValueMap<float>>(p + "DimuonMass");
  t.mmMassErr = produces<edm::ValueMap<float>>(p + "DimuonMassErr");
  // Vertex block, prefixed per arm (rawVtxX / refVtxX ...). The geometry
  // module takes the prefix as configuration, so one implementation serves
  // every arm and every channel.
  t.vtxX = produces<edm::ValueMap<float>>(p + "VtxX");
  t.vtxY = produces<edm::ValueMap<float>>(p + "VtxY");
  t.vtxZ = produces<edm::ValueMap<float>>(p + "VtxZ");
  t.cXX = produces<edm::ValueMap<float>>(p + "VtxCovXX");
  t.cXY = produces<edm::ValueMap<float>>(p + "VtxCovXY");
  t.cXZ = produces<edm::ValueMap<float>>(p + "VtxCovXZ");
  t.cYY = produces<edm::ValueMap<float>>(p + "VtxCovYY");
  t.cYZ = produces<edm::ValueMap<float>>(p + "VtxCovYZ");
  t.cZZ = produces<edm::ValueMap<float>>(p + "VtxCovZZ");
  t.mPt = produces<edm::ValueMap<float>>(p + "MotherPt");
  t.mEta = produces<edm::ValueMap<float>>(p + "MotherEta");
  t.mPhi = produces<edm::ValueMap<float>>(p + "MotherPhi");
  return t;
}

JpsiXKinematicFitProducer::JpsiXKinematicFitProducer(const edm::ParameterSet& cfg)
    : candToken_(consumes<reco::VertexCompositeCandidateCollection>(
          cfg.getParameter<edm::InputTag>("src"))),
      ttbToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
      useRefitTracks_(!cfg.getParameter<std::vector<edm::ParameterSet>>("refitLegs").empty()),
      useBeamSpot_(!cfg.getParameter<edm::InputTag>("beamSpot").label().empty()),
      usePV_(!cfg.getParameter<edm::InputTag>("primaryVertices").label().empty()),
      haveGen_(!cfg.getParameter<edm::InputTag>("genParticles").label().empty()),
      mode_(cfg.getParameter<std::string>("jpsiConstraint")),
      jpsiMass_(cfg.getParameter<double>("jpsiMass")),
      maxChi2_(cfg.getParameter<double>("maxChi2")),
      genLegMaxDR_(cfg.getUntrackedParameter<double>("genLegMaxDR", 0.02)),
      genLegMaxRelDPt_(cfg.getUntrackedParameter<double>("genLegMaxRelDPt", 0.1)),
      genAncestorMaxDepth_(
          cfg.getUntrackedParameter<unsigned>("genAncestorMaxDepth", 10u)) {
  if (mode_ != "inFit" && mode_ != "upstream" && mode_ != "cascade")
    throw cms::Exception("Configuration")
        << "jpsiConstraint must be inFit|upstream|cascade, got '" << mode_ << "'";
  for (const auto& leg : cfg.getParameter<std::vector<edm::ParameterSet>>("refitLegs")) {
    RefitLegTokens t;
    t.tracks = consumes<reco::TrackCollection>(leg.getParameter<edm::InputTag>("tracks"));
    t.refitOk = consumes<edm::ValueMap<int>>(leg.getParameter<edm::InputTag>("refitOk"));
    t.candIdx = consumes<std::vector<int>>(leg.getParameter<edm::InputTag>("candIdx"));
    t.leafIdx = consumes<std::vector<int>>(leg.getParameter<edm::InputTag>("leafIdx"));
    refitLegs_.push_back(t);
  }
  if (useBeamSpot_)
    bsToken_ = consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamSpot"));
  if (usePV_)
    pvToken_ = consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("primaryVertices"));
  if (haveGen_)
    genToken_ = consumes<edm::View<reco::GenParticle>>(cfg.getParameter<edm::InputTag>("genParticles"));

  outGenBMass_ = produces<edm::ValueMap<float>>("genBMass");
  outGenBPt_ = produces<edm::ValueMap<float>>("genBPt");
  outGenBEta_ = produces<edm::ValueMap<float>>("genBEta");
  outGenBPhi_ = produces<edm::ValueMap<float>>("genBPhi");
  outGenBDR_ = produces<edm::ValueMap<float>>("genBDR");
  outGenBPdgId_ = produces<edm::ValueMap<int>>("genBPdgId");
  outGenBIdx_ = produces<edm::ValueMap<int>>("genBIdx");
  for (int i = 0; i < 3; ++i) {
    const std::string L = "leg" + std::to_string(i);
    outLegGenIdx_[i] = produces<edm::ValueMap<int>>(L + "GenIdx");
    outLegGenPdgId_[i] = produces<edm::ValueMap<int>>(L + "GenPdgId");
    outLegGenMotherPdgId_[i] = produces<edm::ValueMap<int>>(L + "GenMotherPdgId");
    outLegGenPt_[i] = produces<edm::ValueMap<float>>(L + "GenPt");
    outLegGenDR_[i] = produces<edm::ValueMap<float>>(L + "GenDR");
  }
  outAncPdgId_ = produces<edm::ValueMap<int>>("genAncestorPdgId");
  outAncIdx_ = produces<edm::ValueMap<int>>("genAncestorIdx");
  outNLegsMatched_ = produces<edm::ValueMap<int>>("nLegsGenMatched");
  outAncPt_ = produces<edm::ValueMap<float>>("genAncestorPt");
  outAncMass_ = produces<edm::ValueMap<float>>("genAncestorMass");
  outNLegsRefit_ = produces<edm::ValueMap<int>>("nLegsRefit");
  raw_ = declareArm("raw");
  ref_ = declareArm("ref");
}

void JpsiXKinematicFitProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<reco::VertexCompositeCandidateCollection> candH;
  iEvent.getByToken(candToken_, candH);
  const auto& ttb = iSetup.getData(ttbToken_);

  edm::Handle<reco::BeamSpot> bsH;
  if (useBeamSpot_) iEvent.getByToken(bsToken_, bsH);
  const bool haveBS = useBeamSpot_ && bsH.isValid();

  edm::Handle<reco::VertexCollection> pvH;
  if (usePV_) iEvent.getByToken(pvToken_, pvH);
  const bool havePV = usePV_ && pvH.isValid() && !pvH->empty();

  edm::Handle<edm::View<reco::GenParticle>> genH;
  if (haveGen_) iEvent.getByToken(genToken_, genH);
  const bool haveGen = haveGen_ && genH.isValid();

  const size_t n = candH->size();

  // Assemble the per-candidate leaf -> refit-track map from every leg source.
  // refitByCand[ic][leafIdx] = refit track (only for legs with refitOk = 1).
  std::vector<std::unordered_map<int, const reco::Track*>> refitByCand(n);
  for (const auto& leg : refitLegs_) {
    edm::Handle<reco::TrackCollection> trkH;
    edm::Handle<edm::ValueMap<int>> okH;
    edm::Handle<std::vector<int>> candIdxH;
    edm::Handle<std::vector<int>> leafIdxH;
    iEvent.getByToken(leg.tracks, trkH);
    iEvent.getByToken(leg.refitOk, okH);
    iEvent.getByToken(leg.candIdx, candIdxH);
    iEvent.getByToken(leg.leafIdx, leafIdxH);
    if (!trkH.isValid() || !okH.isValid() || !candIdxH.isValid() || !leafIdxH.isValid())
      continue;
    const size_t m = trkH->size();
    if (candIdxH->size() != m || leafIdxH->size() != m) continue;  // desync guard
    for (size_t j = 0; j < m; ++j) {
      const reco::TrackRef rref(trkH, j);
      if ((*okH)[rref] != 1) continue;  // leg not refit -> raw fallback downstream
      const int ic = (*candIdxH)[j];
      const int g = (*leafIdxH)[j];
      if (ic >= 0 && static_cast<size_t>(ic) < n)
        refitByCand[ic][g] = &(*trkH)[j];
    }
  }

  // Sentinel-filled: one entry per input candidate, so a failed fit stays
  // visible as fitOk = 0 rather than vanishing from the collection.
  std::vector<float> rMass(n, -99.f), rMassErr(n, -99.f), rPt(n, -99.f), rEta(n, -99.f),
      rPhi(n, -99.f), rVchi2(n, -99.f), rVndof(n, -99.f), rVprob(n, -99.f);
  std::vector<int> rOk(n, 0);
  std::vector<float> rMmVtxProb(n, -99.f), rMmAlphaBS(n, -99.f), rMmSl3d(n, -99.f),
      rMmSl3dPV(n, -99.f), rMmMass(n, -99.f), rMmMassErr(n, -99.f);
  std::vector<float> fMass(n, -99.f), fMassErr(n, -99.f), fPt(n, -99.f), fEta(n, -99.f),
      fPhi(n, -99.f), fVchi2(n, -99.f), fVndof(n, -99.f), fVprob(n, -99.f);
  std::vector<int> fOk(n, 0);
  std::vector<float> fMmVtxProb(n, -99.f), fMmAlphaBS(n, -99.f), fMmSl3d(n, -99.f),
      fMmSl3dPV(n, -99.f), fMmMass(n, -99.f), fMmMassErr(n, -99.f);
  BlockVecs rBlock(n), fBlock(n);
  std::vector<int> nLegsRefit(n, 0);
  std::vector<float> genBMass(n, -99.f), genBPt(n, -99.f), genBEta(n, -99.f),
      genBPhi(n, -99.f), genBDR(n, 9.9f);
  std::vector<int> genBPdgId(n, 0), genBIdx(n, -1);
  // Per-leg generator match. Leg order is the canonical collectLeaves order,
  // i.e. the same ordering that fills mu0/mu1/kaon TrackIdx downstream.
  std::array<std::vector<int>, 3> legGenIdx, legGenPdgId, legGenMotherPdgId;
  std::array<std::vector<float>, 3> legGenPt, legGenDR;
  for (int i = 0; i < 3; ++i) {
    legGenIdx[i].assign(n, -1);
    legGenPdgId[i].assign(n, 0);
    legGenMotherPdgId[i].assign(n, 0);
    legGenPt[i].assign(n, -99.f);
    legGenDR[i].assign(n, -99.f);
  }
  std::vector<int> ancPdgId(n, 0), ancIdx(n, -1), nLegsMatched(n, 0);
  std::vector<float> ancPt(n, -99.f), ancMass(n, -99.f);

  KinematicParticleFactoryFromTransientTrack factory;

  // One arm's worth of work for a candidate: build KinematicParticles from the
  // given per-leaf tracks, do the dimuon-quality fit, then the constrained
  // mother fit. Pure function of (tracks, masses) + the event-level handles.
  auto runArm = [&](const std::vector<const reco::Track*>& trks,
                    const std::vector<double>& masses) -> ArmOut {
    ArmOut a;
    std::vector<reco::TransientTrack> tts;
    tts.reserve(trks.size());
    for (const auto* t : trks) tts.push_back(ttb.build(*t));

    // Dimuon (J/psi) fit-quality handles the analysis path cuts on. An
    // independent unconstrained vertex fit of the two muon legs -> vertex
    // probability, XY pointing angle wrt the beamspot (alphaBS), 2D transverse
    // significance wrt the beamspot (Sxy) and true 3D significance wrt the PV.
    try {
      std::vector<RefCountedKinematicParticle> mm;
      for (size_t il = 0; il < 2; ++il) {
        float m = kMuonMass, s = kMuonSigma, c2 = 0., nd = 0.;
        mm.push_back(factory.particle(tts[il], m, c2, nd, s));
      }
      KinematicParticleVertexFitter mmFit;
      RefCountedKinematicTree mmTree = mmFit.fit(mm);
      if (!mmTree->isEmpty() && mmTree->isValid()) {
        mmTree->movePointerToTheTop();
        auto mmPart = mmTree->currentParticle();
        auto mmVtx = mmTree->currentDecayVertex();
        if (mmPart->currentState().isValid() && mmVtx->vertexIsValid()) {
          // Free dimuon mass -- no constraint was applied in this fit.
          a.mmMass = static_cast<float>(mmPart->currentState().mass());
          {
            const double mErr2 = mmPart->currentState().kinematicParametersError()
                                     .matrix()(6, 6);
            if (mErr2 > 0.) a.mmMassErr = static_cast<float>(std::sqrt(mErr2));
          }
          const double c2 = mmVtx->chiSquared();
          const double nd = mmVtx->degreesOfFreedom();
          if (nd > 0.)
            a.mmVtxProb = TMath::Prob(c2, static_cast<int>(std::lround(nd)));
          const auto sv = mmVtx->position();
          const auto mom = mmPart->currentState().globalMomentum();
          if (haveBS) {
            const double dx = sv.x() - bsH->x(sv.z());
            const double dy = sv.y() - bsH->y(sv.z());
            const double lxy = std::hypot(dx, dy);
            const double ptv = mom.perp();
            if (lxy > 0. && ptv > 0.) {
              const double cosA = (dx * mom.x() + dy * mom.y()) / (lxy * ptv);
              a.mmAlphaBS = std::acos(std::max(-1.0, std::min(1.0, cosA)));
            }
            if (lxy > 0.) {
              const auto e = mmVtx->error().matrix();
              const double ux = dx / lxy, uy = dy / lxy;
              const double bw2 =
                  bsH->BeamWidthX() * bsH->BeamWidthX() * ux * ux +
                  bsH->BeamWidthY() * bsH->BeamWidthY() * uy * uy;
              const double var = ux * ux * e(0, 0) + uy * uy * e(1, 1) +
                                 2 * ux * uy * e(0, 1) + bw2;
              if (var > 0.) a.mmSl3d = lxy / std::sqrt(var);
            }
          }
          // True 3D flight-length significance of the dimuon vertex wrt the
          // production PV (chosen the Bmm5 way: the good PV that MINIMIZES the
          // candidate trajectory's 3D impact parameter). Significance uses the
          // full 3D covariance of both vertices projected on the flight dir.
          if (havePV) {
            const reco::TransientTrack mmTT = mmPart->refittedTransientTrack();
            const reco::Vertex* bestPV = nullptr;
            double bestIP = 1e18;
            for (const auto& pv : *pvH) {
              if (!pv.isValid() || pv.isFake() || pv.ndof() < 4.) continue;
              const auto ip3d = IPTools::absoluteImpactParameter3D(mmTT, pv);
              if (ip3d.first && ip3d.second.value() < bestIP) {
                bestIP = ip3d.second.value();
                bestPV = &pv;
              }
            }
            if (bestPV) {
              const double dx = sv.x() - bestPV->x();
              const double dy = sv.y() - bestPV->y();
              const double dz = sv.z() - bestPV->z();
              const double l3d = std::sqrt(dx * dx + dy * dy + dz * dz);
              if (l3d > 0.) {
                const auto e = mmVtx->error().matrix();
                const double ux = dx / l3d, uy = dy / l3d, uz = dz / l3d;
                const double c00 = e(0, 0) + bestPV->covariance(0, 0);
                const double c11 = e(1, 1) + bestPV->covariance(1, 1);
                const double c22 = e(2, 2) + bestPV->covariance(2, 2);
                const double c01 = e(0, 1) + bestPV->covariance(0, 1);
                const double c02 = e(0, 2) + bestPV->covariance(0, 2);
                const double c12 = e(1, 2) + bestPV->covariance(1, 2);
                const double var = ux * ux * c00 + uy * uy * c11 + uz * uz * c22 +
                                   2 * (ux * uy * c01 + ux * uz * c02 + uy * uz * c12);
                if (var > 0.) a.mmSl3dPV = l3d / std::sqrt(var);
              }
            }
          }
        }
      }
    } catch (const std::exception&) {
      // leave the dimuon sentinels
    }

    try {
      RefCountedKinematicTree tree;
      if (mode_ == "cascade") {
        std::vector<RefCountedKinematicParticle> mumu;
        for (size_t il = 0; il < 2; ++il) {
          float m = kMuonMass, s = kMuonSigma, c2 = 0., nd = 0.;
          mumu.push_back(factory.particle(tts[il], m, c2, nd, s));
        }
        KinematicParticleVertexFitter vfit;
        RefCountedKinematicTree jtree = vfit.fit(mumu);
        if (jtree->isEmpty() || !jtree->isValid()) return a;
        KinematicParticleFitter csFitter;
        auto jpsiC = std::make_unique<MassKinematicConstraint>(
            static_cast<float>(jpsiMass_), static_cast<float>(kMuonSigma));
        jtree = csFitter.fit(jpsiC.get(), jtree);
        if (jtree->isEmpty() || !jtree->isValid()) return a;
        jtree->movePointerToTheTop();
        RefCountedKinematicParticle jpsiPart = jtree->currentParticle();

        std::vector<RefCountedKinematicParticle> parts;
        parts.push_back(jpsiPart);
        for (size_t il = 2; il < tts.size(); ++il) {
          float m = masses[il], s = 1.e-6, c2 = 0., nd = 0.;
          parts.push_back(factory.particle(tts[il], m, c2, nd, s));
        }
        KinematicParticleVertexFitter bfit;
        tree = bfit.fit(parts);
      } else {
        std::vector<RefCountedKinematicParticle> parts;
        for (size_t il = 0; il < tts.size(); ++il) {
          float m = masses[il], s = (il < 2 ? kMuonSigma : 1.e-6), c2 = 0., nd = 0.;
          parts.push_back(factory.particle(tts[il], m, c2, nd, s));
        }
        KinematicConstrainedVertexFitter kcv;
        if (mode_ == "inFit") {
          ParticleMass jm = jpsiMass_;  // ctor takes a non-const ref
          TwoTrackMassKinematicConstraint jpsiC(jm);
          tree = kcv.fit(parts, &jpsiC);
        } else {
          tree = kcv.fit(parts);
        }
      }

      if (tree->isEmpty() || !tree->isValid()) return a;
      tree->movePointerToTheTop();
      RefCountedKinematicParticle top = tree->currentParticle();
      RefCountedKinematicVertex vtx = tree->currentDecayVertex();
      if (!top->currentState().isValid() || !vtx->vertexIsValid()) return a;

      const double chi2 = vtx->chiSquared();
      const double ndof = vtx->degreesOfFreedom();
      if (maxChi2_ > 0. && chi2 > maxChi2_) return a;

      const auto st = top->currentState();
      const auto p4 = st.globalMomentum();
      a.mass = st.mass();
      a.massErr = std::sqrt(std::max(0., st.kinematicParametersError().matrix()(6, 6)));
      a.pt = p4.perp();
      a.eta = p4.eta();
      a.phi = p4.phi();
      a.vchi2 = chi2;
      a.vndof = ndof;
      // Vertex block. KinematicVertex gives a genuine fitted position and its
      // 3x3 error, which is exactly what the geometry module needs.
      {
        const GlobalPoint vp = vtx->position();
        a.vtxX = vp.x(); a.vtxY = vp.y(); a.vtxZ = vp.z();
        const auto ve = vtx->error().matrix();
        a.cXX = ve(0, 0); a.cXY = ve(0, 1); a.cXZ = ve(0, 2);
        a.cYY = ve(1, 1); a.cYZ = ve(1, 2); a.cZZ = ve(2, 2);
      }
      a.vprob = (ndof > 0.) ? TMath::Prob(chi2, static_cast<int>(std::lround(ndof))) : -1.;
      a.ok = 1;
    } catch (const std::exception& e) {
      LogDebug("JpsiXKinematicFitProducer") << "fit failed: " << e.what();
    }
    return a;
  };

  for (size_t ic = 0; ic < n; ++ic) {
    const auto& cand = (*candH)[ic];

    // Generator match (MC only): the closest last-copy b-hadron in dR to the
    // raw candidate direction. Arm-independent; done before any skip.
    if (haveGen) {
      double bestDR = 9.9;
      for (size_t ig = 0; ig < genH->size(); ++ig) {
        const auto& g = (*genH)[ig];
        if (!isBHadron(g.pdgId()) || !g.statusFlags().isLastCopy()) continue;
        const double dr = reco::deltaR(cand, g);
        if (dr < bestDR) {
          bestDR = dr;
          genBDR[ic] = static_cast<float>(dr);
          genBMass[ic] = g.mass();
          genBPt[ic] = g.pt();
          genBEta[ic] = g.eta();
          genBPhi[ic] = g.phi();
          genBPdgId[ic] = g.pdgId();
          genBIdx[ic] = static_cast<int>(ig);
        }
      }
    }

    std::vector<const reco::RecoChargedCandidate*> leaves;
    collectLeaves(&cand, leaves);
    if (leaves.size() < 3) continue;  // need a dimuon plus at least one bachelor

    // Raw per-leaf tracks (stage-1) and the refit substitution (leaf-keyed).
    std::vector<const reco::Track*> rawTrk, refTrk;
    std::vector<double> masses;
    bool badleg = false;
    int nRefit = 0;
    const auto& refMap = refitByCand[ic];
    for (size_t il = 0; il < leaves.size(); ++il) {
      if (leaves[il]->track().isNull()) { badleg = true; break; }
      const reco::Track* raw = &(*leaves[il]->track());
      rawTrk.push_back(raw);
      auto it = refMap.find(static_cast<int>(il));
      if (it != refMap.end()) { refTrk.push_back(it->second); ++nRefit; }
      else refTrk.push_back(raw);
      masses.push_back(massForPdgId(leaves[il]->pdgId()));
    }
    if (badleg) continue;
    nLegsRefit[ic] = nRefit;

    // ---- per-leg generator match + common ancestor (MC only) --------------
    if (haveGen && genH.isValid()) {
      std::vector<const reco::Candidate*> matched;
      matched.reserve(3);
      const size_t nLeg = std::min<size_t>(3, leaves.size());
      for (size_t il = 0; il < nLeg; ++il) {
        const reco::RecoChargedCandidate& leg = *leaves[il];
        const reco::GenParticle* best = nullptr;
        int bestIdx = -1;
        double bestDR = 1e9;
        for (size_t ig = 0; ig < genH->size(); ++ig) {
          const auto& g = (*genH)[ig];
          if (g.status() != 1) continue;             // final state only
          if (g.charge() != leg.charge()) continue;  // charge must agree
          if (!genDrMatch(leg, g, genLegMaxDR_, genLegMaxRelDPt_)) continue;
          const double dr = reco::deltaR(leg, g);
          if (dr < bestDR) { bestDR = dr; best = &g; bestIdx = static_cast<int>(ig); }
        }
        if (best == nullptr) continue;  // stays -1; never a forced nearest neighbour
        legGenIdx[il][ic] = bestIdx;
        legGenPdgId[il][ic] = best->pdgId();
        legGenPt[il][ic] = best->pt();
        legGenDR[il][ic] = static_cast<float>(bestDR);
        if (best->mother() != nullptr)
          legGenMotherPdgId[il][ic] = best->mother()->pdgId();
        matched.push_back(best);
      }
      nLegsMatched[ic] = static_cast<int>(matched.size());
      // Only ask for an ancestor when every leg matched; a common ancestor of
      // two legs out of three says nothing about what the candidate is.
      if (matched.size() == nLeg && nLeg >= 3) {
        const reco::Candidate* anc = findCommonAncestor(matched, genAncestorMaxDepth_);
        if (anc != nullptr) {
          ancPdgId[ic] = anc->pdgId();
          ancPt[ic] = anc->pt();
          ancMass[ic] = anc->mass();
          // Row in the Gen table, which is 1:1 with genParticles (the table is
          // produced with cut=''), found by identity against the same view.
          for (size_t ig = 0; ig < genH->size(); ++ig) {
            if (&(*genH)[ig] == anc) { ancIdx[ic] = static_cast<int>(ig); break; }
          }
        }
      }
    }

    ArmOut raw = runArm(rawTrk, masses);
    ArmOut ref = useRefitTracks_ ? runArm(refTrk, masses) : raw;

    auto store = [&](const ArmOut& a, size_t i,
                     std::vector<float>& M, std::vector<float>& ME, std::vector<float>& P,
                     std::vector<float>& E, std::vector<float>& PH, std::vector<float>& VC,
                     std::vector<float>& VN, std::vector<float>& VP, std::vector<int>& OK,
                     std::vector<float>& MV, std::vector<float>& MA, std::vector<float>& MS,
                     std::vector<float>& MSP, std::vector<float>& MM,
                     std::vector<float>& MME) {
      M[i] = a.mass; ME[i] = a.massErr; P[i] = a.pt; E[i] = a.eta; PH[i] = a.phi;
      VC[i] = a.vchi2; VN[i] = a.vndof; VP[i] = a.vprob; OK[i] = a.ok;
      MV[i] = a.mmVtxProb; MA[i] = a.mmAlphaBS; MS[i] = a.mmSl3d; MSP[i] = a.mmSl3dPV;
      MM[i] = a.mmMass; MME[i] = a.mmMassErr;
    };
    store(raw, ic, rMass, rMassErr, rPt, rEta, rPhi, rVchi2, rVndof, rVprob, rOk,
          rMmVtxProb, rMmAlphaBS, rMmSl3d, rMmSl3dPV, rMmMass, rMmMassErr);
    store(ref, ic, fMass, fMassErr, fPt, fEta, fPhi, fVchi2, fVndof, fVprob, fOk,
          fMmVtxProb, fMmAlphaBS, fMmSl3d, fMmSl3dPV, fMmMass, fMmMassErr);
    rBlock.store(ic, raw);
    fBlock.store(ic, ref);
  }

  auto put = [&](auto token, const auto& vals) {
    using T = typename std::decay_t<decltype(vals)>::value_type;
    auto out = std::make_unique<edm::ValueMap<T>>();
    typename edm::ValueMap<T>::Filler f(*out);
    f.insert(candH, vals.begin(), vals.end());
    f.fill();
    iEvent.put(token, std::move(out));
  };
  auto putArm = [&](const ArmTokens& t,
                    const std::vector<float>& M, const std::vector<float>& ME,
                    const std::vector<float>& P, const std::vector<float>& E,
                    const std::vector<float>& PH, const std::vector<float>& VC,
                    const std::vector<float>& VN, const std::vector<float>& VP,
                    const std::vector<int>& OK, const std::vector<float>& MV,
                    const std::vector<float>& MA, const std::vector<float>& MS,
                    const std::vector<float>& MSP, const std::vector<float>& MM,
                    const std::vector<float>& MME) {
    put(t.mass, M); put(t.massErr, ME); put(t.pt, P); put(t.eta, E); put(t.phi, PH);
    put(t.vchi2, VC); put(t.vndof, VN); put(t.vprob, VP); put(t.ok, OK);
    put(t.mmVtxProb, MV); put(t.mmAlphaBS, MA); put(t.mmSl3d, MS); put(t.mmSl3dPV, MSP);
    put(t.mmMass, MM); put(t.mmMassErr, MME);
  };
  putArm(raw_, rMass, rMassErr, rPt, rEta, rPhi, rVchi2, rVndof, rVprob, rOk,
         rMmVtxProb, rMmAlphaBS, rMmSl3d, rMmSl3dPV, rMmMass, rMmMassErr);
  putArm(ref_, fMass, fMassErr, fPt, fEta, fPhi, fVchi2, fVndof, fVprob, fOk,
         fMmVtxProb, fMmAlphaBS, fMmSl3d, fMmSl3dPV, fMmMass, fMmMassErr);
  auto putBlock = [&](const ArmTokens& t, const BlockVecs& b) {
    put(t.vtxX, b.x); put(t.vtxY, b.y); put(t.vtxZ, b.z);
    put(t.cXX, b.cxx); put(t.cXY, b.cxy); put(t.cXZ, b.cxz);
    put(t.cYY, b.cyy); put(t.cYZ, b.cyz); put(t.cZZ, b.czz);
    put(t.mPt, b.pt); put(t.mEta, b.eta); put(t.mPhi, b.phi);
  };
  putBlock(raw_, rBlock);
  putBlock(ref_, fBlock);
  put(outNLegsRefit_, nLegsRefit);
  put(outGenBMass_, genBMass);
  put(outGenBPt_, genBPt);
  put(outGenBEta_, genBEta);
  put(outGenBPhi_, genBPhi);
  put(outGenBDR_, genBDR);
  put(outGenBPdgId_, genBPdgId);
  put(outGenBIdx_, genBIdx);
  for (int i = 0; i < 3; ++i) {
    put(outLegGenIdx_[i], legGenIdx[i]);
    put(outLegGenPdgId_[i], legGenPdgId[i]);
    put(outLegGenMotherPdgId_[i], legGenMotherPdgId[i]);
    put(outLegGenPt_[i], legGenPt[i]);
    put(outLegGenDR_[i], legGenDR[i]);
  }
  put(outAncPdgId_, ancPdgId);
  put(outAncIdx_, ancIdx);
  put(outNLegsMatched_, nLegsMatched);
  put(outAncPt_, ancPt);
  put(outAncMass_, ancMass);
}

}  // namespace ana_hitanalyzer

using ana_hitanalyzer::JpsiXKinematicFitProducer;
DEFINE_FWK_MODULE(JpsiXKinematicFitProducer);
