/// \class AlignmentTracksFromV0Selector
///
/// EDFilter that produces a self-contained track collection (with track
/// extras, rec hits, and strip + pixel clusters cloned via
/// helper::TrackCollectionStoreManager) containing only the unique daughter
/// tracks of V0 candidates from generalV0Candidates.
///
/// Implementation pattern parallels AlignmentTrackSelectorModule: a small
/// "selector struct" with a select() method is wrapped via
/// ObjectSelectorStream<> into an edm::stream::EDFilter, inheriting all five
/// produces<>() calls from helper::TrackSelectorBase and the cloneAndStore +
/// put() machinery from helper::TrackCollectionStoreManager. This keeps the
/// output byte-for-byte aligned with AlignmentTrackSelectorModule's, and any
/// future additions to the base / store-manager are picked up automatically.
///
/// Used by the ALCARECOTkAlKsToPiPi and ALCARECOTkAlLambdaToProtonPi paths
/// to store only V0-relevant tracks for the alignment fit, avoiding the
/// combinatorial track-pair search done by AlignmentTwoBodyDecayTrackSelector
/// when run on the full generalTracks collection.
///
/// Parameters:
///   src   : InputTag - generalTracks (the framework's main input; required)
///   v0src : InputTag - generalV0Candidates:Kshort or :Lambda
///   filter: bool     - if True, drop the event if zero daughter tracks survive
///                       (default False; in our path the V0 CandViewCountFilter
///                       upstream already guarantees >=1 V0 candidate)

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "CommonTools/UtilAlgos/interface/ObjectSelectorStream.h"
#include "CommonTools/RecoAlgos/interface/TrackSelector.h"   // TrackCollectionStoreManager + TrackSelectorBase + StoreManagerTrait specialization

#include <set>
#include <vector>

namespace {

  struct V0DaughterTrackConfigSelector {
    typedef std::vector<const reco::Track*> container;
    typedef container::const_iterator       const_iterator;
    typedef reco::TrackCollection           collection;

    V0DaughterTrackConfigSelector(const edm::ParameterSet& cfg, edm::ConsumesCollector&& iC)
        : v0Token_(iC.consumes<reco::VertexCompositeCandidateCollection>(
              cfg.getParameter<edm::InputTag>("v0src"))) {}

    const_iterator begin() const { return selected_.begin(); }
    const_iterator end()   const { return selected_.end(); }
    size_t         size()  const { return selected_.size(); }

    void select(const edm::Handle<reco::TrackCollection>& trackHandle,
                const edm::Event& evt, const edm::EventSetup&) {
      selected_.clear();
      edm::Handle<reco::VertexCompositeCandidateCollection> v0Handle;
      evt.getByToken(v0Token_, v0Handle);

      // Collect unique daughter track keys. V0Producer's daughters are
      // RecoChargedCandidate wrapping a TrackRef into the same generalTracks
      // collection that the framework just handed us as `trackHandle`.
      std::set<size_t> daughterKeys;
      for (const auto& cand : *v0Handle) {
        for (size_t i = 0; i < cand.numberOfDaughters(); ++i) {
          const auto* d = dynamic_cast<const reco::RecoChargedCandidate*>(cand.daughter(i));
          if (d && d->track().isNonnull()) {
            daughterKeys.insert(d->track().key());
          }
        }
      }
      selected_.reserve(daughterKeys.size());
      for (auto key : daughterKeys) {
        if (key < trackHandle->size()) {
          selected_.push_back(&(*trackHandle)[key]);
        }
      }
    }

  private:
    container selected_;
    edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> v0Token_;
  };

}  // namespace

// ObjectSelectorStream<> resolves to ObjectSelector<> with Base =
// helper::TrackSelectorBase and StoreManager = helper::TrackCollectionStoreManager
// (via the StoreManagerTrait<reco::TrackCollection, ...> specialization in
// CommonTools/RecoAlgos/interface/TrackSelector.h). This gives us the five
// produces<>() calls and the cloneAndStore + put() flow for free.
typedef ObjectSelectorStream<V0DaughterTrackConfigSelector> AlignmentTracksFromV0Selector;
DEFINE_FWK_MODULE(AlignmentTracksFromV0Selector);
