/// \class CompositeDaughterTrackProducer
///
/// EDProducer that accepts one or more VertexCompositeCandidateCollection
/// inputs (each potentially nested, e.g. B -> J/psi(mu mu) K*0(K pi)),
/// recursively extracts all leaf RecoChargedCandidate tracks, deduplicates
/// by (ProductID, track index), and emits one merged TrackCollection.
///
/// Each output Track is a shallow copy of the original generalTracks entry:
/// the TrackExtraRef is preserved intact so that extra().key() still equals
/// the generalTracks index. This invariant is required by both:
///
///   AlignmentTrackSelectorWithIndexMapModule -- which uses
///     TrackCollectionStoreManager to deep-clone hits/extras; the extra()
///     reference guides the hit copy.
///
///   VertexCompositeCandidateRemapper -- which builds its origToSel lookup
///     as: origToSel[intermediate[originalIndex[j]].extra().key()] = j, and
///     then resolves candidate daughter.track().key() (a generalTracks index)
///     to a cloned-collection index j.
///
/// All inputs are assumed to carry TrackRefs pointing into the same
/// generalTracks collection. Deduplication uses (ProductID, key) pairs so
/// it also handles the (rare) case where inputs reference different upstream
/// collections.

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

#include <set>
#include <utility>
#include <vector>

class CompositeDaughterTrackProducer : public edm::stream::EDProducer<> {
public:
  explicit CompositeDaughterTrackProducer(const edm::ParameterSet& cfg) {
    for (const auto& tag : cfg.getParameter<std::vector<edm::InputTag>>("srcs"))
      srcTokens_.push_back(
          consumes<reco::VertexCompositeCandidateCollection>(tag));
    produces<reco::TrackCollection>();
  }

  void produce(edm::Event& evt, const edm::EventSetup&) override {
    auto outTracks = std::make_unique<reco::TrackCollection>();
    // Deduplicate by (ProductID, key) to handle inputs from different
    // upstream collections while avoiding double-counting within one.
    std::set<std::pair<edm::ProductID, unsigned int>> seen;

    for (const auto& tok : srcTokens_) {
      edm::Handle<reco::VertexCompositeCandidateCollection> vccH;
      evt.getByToken(tok, vccH);
      for (const auto& cand : *vccH)
        extractLeaves(cand, *outTracks, seen);
    }

    LogDebug("CompositeDaughterTrackProducer")
        << "Extracted " << outTracks->size() << " unique leaf tracks.";
    evt.put(std::move(outTracks));
  }

private:
  // Recurse into the candidate tree; push shallow Track copies for each
  // unique leaf RecoChargedCandidate encountered.
  static void extractLeaves(const reco::Candidate& cand,
                             reco::TrackCollection& out,
                             std::set<std::pair<edm::ProductID, unsigned int>>& seen) {
    if (cand.numberOfDaughters() == 0) {
      const auto* rcc = dynamic_cast<const reco::RecoChargedCandidate*>(&cand);
      if (!rcc) return;
      const reco::TrackRef& ref = rcc->track();
      if (!ref.isNonnull()) return;
      // Deduplication key: ProductID identifies the source collection,
      // key() is the track index within that collection.
      auto dedupKey = std::make_pair(ref.id(), ref.key());
      if (!seen.insert(dedupKey).second) return;
      // Shallow copy: the Track object is copied but its internal
      // TrackExtraRef still points into the original (generalTracks)
      // extras collection. extra().key() == generalTracks track index.
      out.push_back(*ref);
    } else {
      for (size_t i = 0; i < cand.numberOfDaughters(); ++i)
        extractLeaves(*cand.daughter(i), out, seen);
    }
  }

  std::vector<edm::EDGetTokenT<reco::VertexCompositeCandidateCollection>> srcTokens_;
};

DEFINE_FWK_MODULE(CompositeDaughterTrackProducer);
