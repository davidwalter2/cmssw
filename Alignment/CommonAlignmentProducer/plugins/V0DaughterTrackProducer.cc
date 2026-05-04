/// \class V0DaughterTrackProducer
///
/// EDProducer that emits a small reco::TrackCollection containing the unique
/// daughter tracks of V0 candidates from generalV0Candidates. Designed as a
/// thin upstream producer for the standard AlignmentTrackSelectorModule
/// chain: feed its output as the `src` of AlignmentTrackSelector and the
/// existing TrackCollectionStoreManager machinery clones the tracks +
/// extras + hits + clusters into the ALCARECO output, exactly as it does
/// when reading directly from generalTracks.
///
/// The trick: when reco::Track is copied into the new collection, its
/// internal reco::TrackExtraRef is preserved -- it still points back to
/// generalTracks::TrackExtras. So the cloning chain inside
/// helper::TrackCollectionStoreManager works unchanged.
///
/// Used by the ALCARECOTkAlKsToPiPi and ALCARECOTkAlLambdaToProtonPi paths
/// to restrict the AlignmentTrackSelector input to V0 daughters only,
/// avoiding the combinatorial track-pair search done by
/// AlignmentTwoBodyDecayTrackSelector when run on the full generalTracks
/// collection.

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include <set>

class V0DaughterTrackProducer : public edm::stream::EDProducer<> {
public:
  explicit V0DaughterTrackProducer(const edm::ParameterSet& cfg)
      : v0Token_(consumes<reco::VertexCompositeCandidateCollection>(
            cfg.getParameter<edm::InputTag>("src"))) {
    produces<reco::TrackCollection>();
  }

  void produce(edm::Event& evt, const edm::EventSetup&) override {
    edm::Handle<reco::VertexCompositeCandidateCollection> v0Handle;
    evt.getByToken(v0Token_, v0Handle);

    auto out = std::make_unique<reco::TrackCollection>();
    std::set<size_t> seen;
    for (const auto& cand : *v0Handle) {
      for (size_t i = 0; i < cand.numberOfDaughters(); ++i) {
        const auto* d = dynamic_cast<const reco::RecoChargedCandidate*>(cand.daughter(i));
        if (!d || d->track().isNull()) continue;
        if (seen.insert(d->track().key()).second) {
          out->push_back(*d->track());  // shallow copy: TrackExtraRef preserved
        }
      }
    }
    evt.put(std::move(out));
  }

private:
  edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> v0Token_;
};

DEFINE_FWK_MODULE(V0DaughterTrackProducer);
