/// \class TwoBodyDecayCandidateProducer
///
/// EDProducer that builds a reco::VertexCompositeCandidateCollection from
/// pairs of reco::Tracks passing a configurable mass window + charge +
/// acoplanarity selection. Used as the upstream candidate-builder for the
/// dimuon ALCAREco streams (TkAlJpsiMuMu, TkAlUpsilonMuMu, TkAlZMuMu): pair
/// finding is moved here so the candidate object survives downstream and the
/// AlignmentTrackSelector (via V0DaughterTrackProducer) only has to apply
/// per-track quality cuts -- mirroring the V0 ALCAREco pattern.
///
/// Each output candidate carries two reco::RecoChargedCandidate daughters
/// whose TrackRefs point into the input track collection. After the
/// AlignmentTrackSelector clones the daughter tracks into the ALCAREco
/// output, VertexCompositeCandidateRemapper rewrites these TrackRefs to
/// point at the cloned collection so downstream consumers can navigate
/// candidate -> daughter -> track without going back to generalTracks.
///
/// Unlike AlignmentTwoBodyDecayTrackSelector, all passing pairs are emitted
/// (multiplicity > 1 allowed); ordering by mother pT is preserved as the
/// natural iteration order.

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include <Math/Vector4D.h>
#include <algorithm>
#include <cmath>
#include <set>
#include <vector>

class TwoBodyDecayCandidateProducer : public edm::stream::EDProducer<> {
public:
  explicit TwoBodyDecayCandidateProducer(const edm::ParameterSet& cfg)
      : trackToken_(consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("src"))),
        muonTag_(cfg.getParameter<edm::InputTag>("muonSrc")),
        useMuonFilter_(!muonTag_.label().empty()),
        minMass_(cfg.getParameter<double>("minMass")),
        maxMass_(cfg.getParameter<double>("maxMass")),
        daughterMass_(cfg.getParameter<double>("daughterMass")),
        daughterPdgId_(cfg.getParameter<int>("daughterPdgId")),
        motherPdgId_(cfg.getParameter<int>("motherPdgId")),
        applyChargeFilter_(cfg.getParameter<bool>("applyChargeFilter")),
        charge_(cfg.getParameter<int>("charge")),
        useUnsignedCharge_(cfg.getParameter<bool>("useUnsignedCharge")),
        applyAcoplanarityFilter_(cfg.getParameter<bool>("applyAcoplanarityFilter")),
        acoplanarDistance_(cfg.getParameter<double>("acoplanarDistance")) {
    if (useMuonFilter_) {
      muonToken_ = consumes<edm::View<reco::Muon>>(muonTag_);
    }
    produces<reco::VertexCompositeCandidateCollection>();
  }

  void produce(edm::Event& evt, const edm::EventSetup&) override {
    edm::Handle<reco::TrackCollection> trackH;
    evt.getByToken(trackToken_, trackH);

    auto out = std::make_unique<reco::VertexCompositeCandidateCollection>();

    // If requested, restrict to tracks that are the innerTrack of a muon in
    // muonSrc. Replicates AlignmentGlobalTrackSelector::applyGlobalMuonFilter
    // upstream of the candidate combinatorics.
    std::vector<size_t> goodIdx;
    goodIdx.reserve(trackH->size());
    if (useMuonFilter_) {
      edm::Handle<edm::View<reco::Muon>> muonH;
      evt.getByToken(muonToken_, muonH);
      std::set<size_t> muTrackKeys;
      for (const auto& mu : *muonH) {
        const auto& tr = mu.innerTrack();
        if (tr.isNonnull() && tr.id() == trackH.id())
          muTrackKeys.insert(tr.key());
      }
      for (size_t i = 0; i < trackH->size(); ++i)
        if (muTrackKeys.count(i)) goodIdx.push_back(i);
    } else {
      for (size_t i = 0; i < trackH->size(); ++i) goodIdx.push_back(i);
    }

    if (goodIdx.size() < 2) {
      evt.put(std::move(out));
      return;
    }

    // Score each passing pair by mother pT to give a deterministic, useful
    // iteration order in the output collection.
    using P4 = ROOT::Math::PxPyPzMVector;
    using PairItem = std::tuple<double, size_t, size_t>;  // (motherPt, idxA, idxB)
    std::vector<PairItem> ranked;
    ranked.reserve(goodIdx.size() * (goodIdx.size() - 1) / 2);

    const double mD = daughterMass_;
    for (size_t a = 0; a < goodIdx.size(); ++a) {
      const reco::Track& trA = (*trackH)[goodIdx[a]];
      P4 p4A(trA.px(), trA.py(), trA.pz(), mD);
      for (size_t b = a + 1; b < goodIdx.size(); ++b) {
        const reco::Track& trB = (*trackH)[goodIdx[b]];

        if (applyChargeFilter_) {
          int sumQ = trA.charge() + trB.charge();
          if (useUnsignedCharge_) sumQ = std::abs(sumQ);
          if (sumQ != charge_) continue;
        }

        if (applyAcoplanarityFilter_) {
          if (std::fabs(reco::deltaPhi(trA.phi(), trB.phi() - M_PI)) >= acoplanarDistance_)
            continue;
        }

        P4 p4B(trB.px(), trB.py(), trB.pz(), mD);
        const auto p4M = p4A + p4B;
        const double mass = p4M.M();
        if (mass <= minMass_ || mass >= maxMass_) continue;

        ranked.emplace_back(p4M.Pt(), goodIdx[a], goodIdx[b]);
      }
    }

    if (ranked.empty()) {
      evt.put(std::move(out));
      return;
    }

    std::sort(ranked.begin(), ranked.end(),
              [](const PairItem& x, const PairItem& y) { return std::get<0>(x) > std::get<0>(y); });

    out->reserve(ranked.size());
    for (const auto& item : ranked) {
      size_t iA = std::get<1>(item);
      size_t iB = std::get<2>(item);
      // Match V0Producer's positive-charge-first convention: when the pair
      // is opposite-sign, daughter(0) carries the positive track. For same-
      // sign pairs (allowed when applyChargeFilter is off, e.g. J/psi
      // QC/control), preserve track-index order.
      if ((*trackH)[iB].charge() > (*trackH)[iA].charge())
        std::swap(iA, iB);

      const reco::Track& trA = (*trackH)[iA];
      const reco::Track& trB = (*trackH)[iB];

      P4 p4A(trA.px(), trA.py(), trA.pz(), mD);
      P4 p4B(trB.px(), trB.py(), trB.pz(), mD);
      const auto p4M = p4A + p4B;

      // Daughter PDG id sign follows track charge; pdgId convention here is
      // such that daughterPdgId_ > 0 represents the negative-charge daughter
      // (e.g. mu- is +13). Sign with -charge() to match.
      const int pdgA = (trA.charge() < 0) ? +daughterPdgId_ : -daughterPdgId_;
      const int pdgB = (trB.charge() < 0) ? +daughterPdgId_ : -daughterPdgId_;

      reco::Particle::LorentzVector lvA(p4A.Px(), p4A.Py(), p4A.Pz(), p4A.E());
      reco::Particle::LorentzVector lvB(p4B.Px(), p4B.Py(), p4B.Pz(), p4B.E());
      const reco::Particle::Point vtxA(trA.vx(), trA.vy(), trA.vz());
      const reco::Particle::Point vtxB(trB.vx(), trB.vy(), trB.vz());

      reco::RecoChargedCandidate dA(trA.charge(), lvA, vtxA, pdgA);
      dA.setTrack(reco::TrackRef(trackH, iA));
      reco::RecoChargedCandidate dB(trB.charge(), lvB, vtxB, pdgB);
      dB.setTrack(reco::TrackRef(trackH, iB));

      // Mother position: midpoint of daughter reference points (no vertex
      // fit at this stage; the post-CVH refit reassigns vertex/kinematics).
      const reco::Particle::Point vtxM(0.5 * (trA.vx() + trB.vx()),
                                       0.5 * (trA.vy() + trB.vy()),
                                       0.5 * (trA.vz() + trB.vz()));
      reco::Particle::LorentzVector lvM(p4M.Px(), p4M.Py(), p4M.Pz(), p4M.E());
      const int chargeM = trA.charge() + trB.charge();

      reco::VertexCompositeCandidate cand(chargeM, lvM, vtxM, motherPdgId_);
      cand.addDaughter(dA);
      cand.addDaughter(dB);
      out->push_back(std::move(cand));
    }

    LogDebug("TwoBodyDecayCandidateProducer")
        << "Built " << out->size() << " candidates from " << goodIdx.size() << " input tracks.";
    evt.put(std::move(out));
  }

private:
  edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  edm::InputTag muonTag_;
  bool useMuonFilter_;
  edm::EDGetTokenT<edm::View<reco::Muon>> muonToken_;

  double minMass_;
  double maxMass_;
  double daughterMass_;
  int daughterPdgId_;
  int motherPdgId_;

  bool applyChargeFilter_;
  int charge_;
  bool useUnsignedCharge_;

  bool applyAcoplanarityFilter_;
  double acoplanarDistance_;
};

DEFINE_FWK_MODULE(TwoBodyDecayCandidateProducer);
