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
// Modes (`jpsiConstraint`), selectable because the two are not equivalent and
// should be compared rather than argued about:
//
//   "inFit"   -- the Bmm5 convention. Take UNCONSTRAINED muon tracks and apply
//                the J/psi mass constraint inside the mother fit, via
//                KinematicConstrainedVertexFitter + TwoTrackMassKinematicConstraint.
//   "upstream"-- take muons already J/psi-mass-constrained by the CVH refit and
//                run the mother fit with NO further mass constraint. Applying
//                it twice would double-count: the J/psi mass error shrinks
//                artificially and the mother mass error is biased low.
//   "cascade" -- fit the dimuon first, collapse it into a single
//                KinematicParticle carrying its full covariance, then fit that
//                against the bachelor with the vertex FLOATING. Equivalent to
//                the simultaneous fit in the Gaussian limit only because the
//                J/psi decays at the mother vertex; it would not be for a
//                long-lived intermediate.
//
// Track source: `srcTracks` is optional. When set, a leaf daughter's track is
// replaced by the positionally-matched entry of that collection (the CVH
// refit output), gated on its `refitOk` map. Unset -> the stage-1 tracks.

#include <memory>
#include <vector>

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

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
// CVH maker uses: a composite daughter recurses, a leaf yields its track.
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
}  // namespace

class JpsiXKinematicFitProducer : public edm::stream::EDProducer<> {
public:
  explicit JpsiXKinematicFitProducer(const edm::ParameterSet&);
  ~JpsiXKinematicFitProducer() override = default;

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  const edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> candToken_;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbToken_;
  edm::EDGetTokenT<reco::TrackCollection> refitTrackToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> refitOkToken_;
  const bool useRefitTracks_;
  const std::string mode_;
  const double jpsiMass_;
  const double maxChi2_;

  edm::EDPutTokenT<edm::ValueMap<float>> outMass_, outMassErr_, outPt_, outEta_,
      outPhi_, outVtxChi2_, outVtxNdof_, outVtxProb_;
  edm::EDPutTokenT<edm::ValueMap<int>> outOk_;
};

JpsiXKinematicFitProducer::JpsiXKinematicFitProducer(const edm::ParameterSet& cfg)
    : candToken_(consumes<reco::VertexCompositeCandidateCollection>(
          cfg.getParameter<edm::InputTag>("src"))),
      ttbToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
      useRefitTracks_(!cfg.getParameter<edm::InputTag>("srcTracks").label().empty()),
      mode_(cfg.getParameter<std::string>("jpsiConstraint")),
      jpsiMass_(cfg.getParameter<double>("jpsiMass")),
      maxChi2_(cfg.getParameter<double>("maxChi2")) {
  if (mode_ != "inFit" && mode_ != "upstream" && mode_ != "cascade")
    throw cms::Exception("Configuration")
        << "jpsiConstraint must be inFit|upstream|cascade, got '" << mode_ << "'";
  if (useRefitTracks_) {
    refitTrackToken_ = consumes<reco::TrackCollection>(
        cfg.getParameter<edm::InputTag>("srcTracks"));
    refitOkToken_ = consumes<edm::ValueMap<int>>(
        cfg.getParameter<edm::InputTag>("srcRefitOk"));
  }
  outMass_ = produces<edm::ValueMap<float>>("fitMass");
  outMassErr_ = produces<edm::ValueMap<float>>("fitMassErr");
  outPt_ = produces<edm::ValueMap<float>>("fitPt");
  outEta_ = produces<edm::ValueMap<float>>("fitEta");
  outPhi_ = produces<edm::ValueMap<float>>("fitPhi");
  outVtxChi2_ = produces<edm::ValueMap<float>>("fitVtxChi2");
  outVtxNdof_ = produces<edm::ValueMap<float>>("fitVtxNdof");
  outVtxProb_ = produces<edm::ValueMap<float>>("fitVtxProb");
  outOk_ = produces<edm::ValueMap<int>>("fitOk");
}

void JpsiXKinematicFitProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<reco::VertexCompositeCandidateCollection> candH;
  iEvent.getByToken(candToken_, candH);
  const auto& ttb = iSetup.getData(ttbToken_);

  edm::Handle<reco::TrackCollection> refitH;
  edm::Handle<edm::ValueMap<int>> refitOkH;
  if (useRefitTracks_) {
    iEvent.getByToken(refitTrackToken_, refitH);
    iEvent.getByToken(refitOkToken_, refitOkH);
  }

  const size_t n = candH->size();
  // Sentinel-filled: one entry per input candidate, so a failed fit stays
  // visible as fitOk = 0 rather than vanishing from the collection.
  std::vector<float> mass(n, -99.f), massErr(n, -99.f), pt(n, -99.f), eta(n, -99.f),
      phi(n, -99.f), vchi2(n, -99.f), vndof(n, -99.f), vprob(n, -99.f);
  std::vector<int> ok(n, 0);

  KinematicParticleFactoryFromTransientTrack factory;

  for (size_t ic = 0; ic < n; ++ic) {
    const auto& cand = (*candH)[ic];

    std::vector<const reco::RecoChargedCandidate*> leaves;
    collectLeaves(&cand, leaves);
    if (leaves.size() < 3) continue;  // need a dimuon plus at least one bachelor

    // Stage-1 orders daughter(0) = the J/psi composite, so its two leaves come
    // first; the remainder are bachelors.
    std::vector<reco::TransientTrack> tts;
    std::vector<double> masses;
    bool badleg = false;
    for (size_t il = 0; il < leaves.size(); ++il) {
      const reco::Track* trk = nullptr;
      if (leaves[il]->track().isNull()) { badleg = true; break; }
      const reco::TrackRef tref = leaves[il]->track();
      // Refit-track substitution requires a refit collection keyed to the SAME
      // track collection the daughters reference. The single-track maker's
      // output is keyed to the per-candidate BACHELOR collection instead, so a
      // tref.key() lookup would cross collections and grab the wrong track --
      // the product-id equality below is what prevents that. It is therefore
      // inert for the current single-track wiring (different product) and will
      // self-enable only when fed a refit collection built from the daughters'
      // own track product (the N-body maker, or a full alignment-track refit).
      if (useRefitTracks_ && refitH.isValid() && tref.id() == refitH.id() &&
          tref.key() < refitH->size()) {
        const reco::TrackRef rref(refitH, tref.key());
        if ((*refitOkH)[rref] == 1) trk = &(*refitH)[tref.key()];
      }
      if (trk == nullptr) trk = &(*tref);
      tts.push_back(ttb.build(*trk));
      masses.push_back(massForPdgId(leaves[il]->pdgId()));
    }
    if (badleg || tts.size() != leaves.size()) continue;

    try {
      RefCountedKinematicTree tree;

      if (mode_ == "cascade") {
        // Fit the dimuon on its own, mass-constrain it, then combine the
        // resulting single particle with the bachelors, vertex floating.
        std::vector<RefCountedKinematicParticle> mumu;
        for (size_t il = 0; il < 2; ++il) {
          float m = kMuonMass, s = kMuonSigma;
          float c2 = 0., nd = 0.;
          mumu.push_back(factory.particle(tts[il], m, c2, nd, s));
        }
        KinematicParticleVertexFitter vfit;
        RefCountedKinematicTree jtree = vfit.fit(mumu);
        if (jtree->isEmpty() || !jtree->isValid()) continue;
        KinematicParticleFitter csFitter;
        auto jpsiC = std::make_unique<MassKinematicConstraint>(
            static_cast<float>(jpsiMass_), static_cast<float>(kMuonSigma));
        jtree = csFitter.fit(jpsiC.get(), jtree);
        if (jtree->isEmpty() || !jtree->isValid()) continue;
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
          // Bmm5: constrain the first two (the J/psi legs) inside the fit.
          ParticleMass jm = jpsiMass_;  // ctor takes a non-const ref
          TwoTrackMassKinematicConstraint jpsiC(jm);
          tree = kcv.fit(parts, &jpsiC);
        } else {
          // "upstream": the muons already carry the CVH J/psi constraint;
          // applying it again here would double-count it.
          tree = kcv.fit(parts);
        }
      }

      if (tree->isEmpty() || !tree->isValid()) continue;
      tree->movePointerToTheTop();
      RefCountedKinematicParticle top = tree->currentParticle();
      RefCountedKinematicVertex vtx = tree->currentDecayVertex();
      if (!top->currentState().isValid() || !vtx->vertexIsValid()) continue;

      const double chi2 = vtx->chiSquared();
      const double ndof = vtx->degreesOfFreedom();
      if (maxChi2_ > 0. && chi2 > maxChi2_) continue;

      const auto st = top->currentState();
      const auto p4 = st.globalMomentum();
      mass[ic] = st.mass();
      massErr[ic] = std::sqrt(std::max(0., st.kinematicParametersError().matrix()(6, 6)));
      pt[ic] = p4.perp();
      eta[ic] = p4.eta();
      phi[ic] = p4.phi();
      vchi2[ic] = chi2;
      vndof[ic] = ndof;
      vprob[ic] = (ndof > 0.) ? TMath::Prob(chi2, static_cast<int>(std::lround(ndof))) : -1.;
      ok[ic] = 1;
    } catch (const std::exception& e) {
      // A failed kinematic fit is expected combinatorics, not an error:
      // leave the sentinels and the fitOk = 0 flag in place.
      LogDebug("JpsiXKinematicFitProducer") << "fit failed: " << e.what();
    }
  }

  auto put = [&](auto token, const auto& vals) {
    using T = typename std::decay_t<decltype(vals)>::value_type;
    auto out = std::make_unique<edm::ValueMap<T>>();
    typename edm::ValueMap<T>::Filler f(*out);
    f.insert(candH, vals.begin(), vals.end());
    f.fill();
    iEvent.put(token, std::move(out));
  };
  put(outMass_, mass);
  put(outMassErr_, massErr);
  put(outPt_, pt);
  put(outEta_, eta);
  put(outPhi_, phi);
  put(outVtxChi2_, vchi2);
  put(outVtxNdof_, vndof);
  put(outVtxProb_, vprob);
  put(outOk_, ok);
}

}  // namespace ana_hitanalyzer

using ana_hitanalyzer::JpsiXKinematicFitProducer;
DEFINE_FWK_MODULE(JpsiXKinematicFitProducer);
