/// \class AlignmentTrackSelectorWithIndexMapModule
///
/// Drop-in replacement for AlignmentTrackSelectorModule that, in addition
/// to the standard cloned outputs (tracks + extras + hits + clusters via
/// helper::TrackCollectionStoreManager), emits an
/// `edm::ValueMap<unsigned int>` keyed on the cloned-track collection.
/// The value at each entry is the index of the corresponding track in the
/// *source* TrackCollection (the module's `src` input).
///
/// Why this exists: TrackCollectionStoreManager rewrites every cloned
/// track's TrackExtraRef to a fresh sequential key into the new selected-
/// extras collection (see CommonTools/RecoAlgos/src/TrackSelector.cc:32),
/// so `selectedTrack.extra().key()` does NOT recover the source-track
/// index. Downstream consumers that need to navigate from a cloned track
/// back into the source collection (DeDxValueMapProjector, the candidate
/// remapper, ...) are broken without the side-channel ValueMap this
/// module provides.
///
/// Implementation: the index is computed by pointer arithmetic on the
/// selector chain's output -- TrackConfigSelector::select() returns
/// `vector<const reco::Track*>` pointing INTO the source collection's
/// contiguous storage, so `(p - &source[0])` is the exact source index.
/// No track-parameter fingerprinting is involved.
///
/// Configuration is identical to AlignmentTrackSelectorModule; just the
/// module label changes. The source-index ValueMap is published with the
/// label "originalIndex".

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "Alignment/CommonAlignmentProducer/interface/AlignmentTrackSelector.h"
#include "Alignment/CommonAlignmentProducer/interface/AlignmentGlobalTrackSelector.h"
#include "Alignment/CommonAlignmentProducer/interface/AlignmentTwoBodyDecayTrackSelector.h"
#include "Alignment/CommonAlignmentProducer/interface/AlignmentThreeBodyDecayTrackSelector.h"

#include "CommonTools/RecoAlgos/interface/TrackSelector.h"

#include "DataFormats/Common/interface/OrphanHandle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"

#include <memory>
#include <vector>

namespace {

// Mirror of the TrackConfigSelector struct from AlignmentTrackSelectorModule.cc;
// kept local to avoid sharing internal-linkage state across plugins.
struct TrackConfigSelector {
  typedef std::vector<const reco::Track*> container;
  typedef container::const_iterator const_iterator;
  typedef reco::TrackCollection collection;

  TrackConfigSelector(const edm::ParameterSet& cfg, edm::ConsumesCollector&& iC)
      : theBaseSelector(cfg, iC),
        theGlobalSelector(cfg.getParameter<edm::ParameterSet>("GlobalSelector"), iC),
        theTwoBodyDecaySelector(cfg.getParameter<edm::ParameterSet>("TwoBodyDecaySelector"), iC),
        theThreeBodyDecaySelector(cfg.getParameter<edm::ParameterSet>("ThreeBodyDecaySelector"), iC),
        theBaseSwitch(theBaseSelector.useThisFilter()),
        theGlobalSwitch(theGlobalSelector.useThisFilter()),
        theTwoBodyDecaySwitch(theTwoBodyDecaySelector.useThisFilter()),
        theThreeBodyDecaySwitch(theThreeBodyDecaySelector.useThisFilter()) {}

  const_iterator begin() const { return theSelectedTracks.begin(); }
  const_iterator end() const { return theSelectedTracks.end(); }
  size_t size() const { return theSelectedTracks.size(); }

  void select(const edm::Handle<reco::TrackCollection>& c, const edm::Event& evt, const edm::EventSetup& es) {
    theSelectedTracks.clear();
    theSelectedTracks.reserve(c->size());
    for (const auto& trk : *c) theSelectedTracks.push_back(&trk);
    if (theBaseSwitch)        theSelectedTracks = theBaseSelector.select(theSelectedTracks, evt, es);
    if (theGlobalSwitch)      theSelectedTracks = theGlobalSelector.select(theSelectedTracks, evt, es);
    if (theTwoBodyDecaySwitch)   theSelectedTracks = theTwoBodyDecaySelector.select(theSelectedTracks, evt, es);
    if (theThreeBodyDecaySwitch) theSelectedTracks = theThreeBodyDecaySelector.select(theSelectedTracks, evt, es);
  }

private:
  container theSelectedTracks;
  AlignmentTrackSelector theBaseSelector;
  AlignmentGlobalTrackSelector theGlobalSelector;
  AlignmentTwoBodyDecayTrackSelector theTwoBodyDecaySelector;
  AlignmentThreeBodyDecayTrackSelector theThreeBodyDecaySelector;
  bool theBaseSwitch, theGlobalSwitch, theTwoBodyDecaySwitch, theThreeBodyDecaySwitch;
};

}  // namespace

class AlignmentTrackSelectorWithIndexMapModule : public edm::stream::EDFilter<> {
public:
  explicit AlignmentTrackSelectorWithIndexMapModule(const edm::ParameterSet& cfg)
      : srcToken_(consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("src"))),
        selector_(cfg, consumesCollector()),
        filterEmpty_(cfg.getParameter<bool>("filter")) {
    const std::string alias = cfg.getParameter<std::string>("@module_label");
    produces<reco::TrackCollection>().setBranchAlias(alias + "Tracks");
    produces<reco::TrackExtraCollection>().setBranchAlias(alias + "TrackExtras");
    produces<TrackingRecHitCollection>().setBranchAlias(alias + "RecHits");
    produces<edmNew::DetSetVector<SiPixelCluster>>().setBranchAlias(alias + "PixelClusters");
    produces<edmNew::DetSetVector<SiStripCluster>>().setBranchAlias(alias + "StripClusters");
    produces<edm::ValueMap<unsigned int>>("originalIndex");
  }

  bool filter(edm::Event& evt, const edm::EventSetup& iSetup) override {
    edm::Handle<reco::TrackCollection> srcH;
    evt.getByToken(srcToken_, srcH);

    selector_.select(srcH, evt, iSetup);

    if (filterEmpty_ && selector_.size() == 0) {
      // Still need to put the (empty) products to satisfy the framework,
      // because this is an EDFilter that DECLARES these products via produces<>.
      auto emptyTracks = std::make_unique<reco::TrackCollection>();
      auto emptyExtras = std::make_unique<reco::TrackExtraCollection>();
      auto emptyHits = std::make_unique<TrackingRecHitCollection>();
      auto emptyPx = std::make_unique<edmNew::DetSetVector<SiPixelCluster>>();
      auto emptySt = std::make_unique<edmNew::DetSetVector<SiStripCluster>>();
      auto emptyVm = std::make_unique<edm::ValueMap<unsigned int>>();
      evt.put(std::move(emptyTracks));
      evt.put(std::move(emptyExtras));
      evt.put(std::move(emptyHits));
      evt.put(std::move(emptyPx));
      evt.put(std::move(emptySt));
      evt.put(std::move(emptyVm), "originalIndex");
      return false;
    }

    // Compute each surviving track's source-collection index by pointer
    // arithmetic against the contiguous source storage.
    std::vector<unsigned int> srcIdx;
    srcIdx.reserve(selector_.size());
    const reco::Track* base = (srcH->empty() ? nullptr : &(*srcH)[0]);
    for (auto it = selector_.begin(); it != selector_.end(); ++it) {
      srcIdx.push_back(static_cast<unsigned int>(*it - base));
    }

    // Standard track + extras + hits + clusters cloning.
    helper::TrackCollectionStoreManager mgr(srcH);
    mgr.cloneAndStore(selector_.begin(), selector_.end(), evt);
    edm::OrphanHandle<reco::TrackCollection> selH = mgr.put(evt);

    auto vm = std::make_unique<edm::ValueMap<unsigned int>>();
    edm::ValueMap<unsigned int>::Filler filler(*vm);
    filler.insert(selH, srcIdx.begin(), srcIdx.end());
    filler.fill();
    evt.put(std::move(vm), "originalIndex");

    return true;
  }

private:
  edm::EDGetTokenT<reco::TrackCollection> srcToken_;
  TrackConfigSelector selector_;
  bool filterEmpty_;
};

DEFINE_FWK_MODULE(AlignmentTrackSelectorWithIndexMapModule);
