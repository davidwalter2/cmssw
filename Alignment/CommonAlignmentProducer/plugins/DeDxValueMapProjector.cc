/// \class DeDxValueMapProjector
///
/// EDProducer that re-keys a ValueMap<DeDxData> from the original
/// generalTracks collection onto a deep-copied selected-track collection
/// (e.g. the ALCARECOTkAl{Ks,Lambda,Dst}* outputs), via an explicit
/// reference chain (no kinematic fingerprinting).
///
/// Reference chain per cloned track at index j of `selectedTracks`:
///
///   selH[j] --(originalIndexMap[selH[j]])--> intermediateTracks[ii]
///          --(intermediateTracks[ii].extra().key())--> generalTracks[gt]
///          --Ref<TrackCollection>(sourceTracks, gt)--> ValueMap lookup
///
/// For the V0 streams (Ks, Lambda), `intermediateTracks` is the
/// V0DaughterTrackProducer output: it stores Tracks BY VALUE preserving
/// the original TrackExtraRef whose key is the generalTracks index. For
/// the D* stream (where AlignmentTrackSelector takes generalTracks
/// directly), `intermediateTracks == sourceTracks == generalTracks` and
/// `extra().key()` is just the self index, so the same chain trivially
/// applies.
///
/// Why this exists at all: helper::TrackCollectionStoreManager (used by
/// AlignmentTrackSelector) rewrites every cloned track's TrackExtraRef
/// to a fresh sequential key into the new selected-extras collection,
/// destroying the link back to generalTracks. AlignmentTrackSelectorWithIndexMapModule
/// emits a side-channel ValueMap<unsigned int> ("originalIndex") that
/// records each cloned track's source-collection index by pointer
/// arithmetic on the selector chain's output -- this is the
/// `originalIndexMap` consumed here.

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

#include <limits>
#include <memory>
#include <vector>

class DeDxValueMapProjector : public edm::global::EDProducer<> {
public:
  explicit DeDxValueMapProjector(const edm::ParameterSet& cfg)
      : selTracksToken_(consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("selectedTracks"))),
        intermediateTracksToken_(
            consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("intermediateTracks"))),
        sourceTracksToken_(consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("sourceTracks"))),
        sourceValueMapToken_(consumes<edm::ValueMap<reco::DeDxData>>(cfg.getParameter<edm::InputTag>("sourceValueMap"))),
        origIndexMapToken_(consumes<edm::ValueMap<unsigned int>>(cfg.getParameter<edm::InputTag>("originalIndexMap"))) {
    produces<edm::ValueMap<float>>();
  }

  void produce(edm::StreamID, edm::Event& evt, const edm::EventSetup&) const override {
    edm::Handle<reco::TrackCollection> selH;
    evt.getByToken(selTracksToken_, selH);
    edm::Handle<reco::TrackCollection> intermH;
    evt.getByToken(intermediateTracksToken_, intermH);
    edm::Handle<reco::TrackCollection> srcH;
    evt.getByToken(sourceTracksToken_, srcH);
    edm::Handle<edm::ValueMap<reco::DeDxData>> vmH;
    evt.getByToken(sourceValueMapToken_, vmH);
    edm::Handle<edm::ValueMap<unsigned int>> origIdxH;
    evt.getByToken(origIndexMapToken_, origIdxH);

    constexpr unsigned int kMissing = std::numeric_limits<unsigned int>::max();
    std::vector<float> values;
    values.reserve(selH->size());
    for (size_t j = 0; j < selH->size(); ++j) {
      const edm::Ref<reco::TrackCollection> selRef(selH, j);
      const unsigned int intermIdx = (*origIdxH)[selRef];
      if (intermIdx == kMissing || intermIdx >= intermH->size()) {
        values.push_back(std::numeric_limits<float>::quiet_NaN());
        continue;
      }
      const auto& imTk = (*intermH)[intermIdx];
      if (imTk.extra().isNull()) {
        values.push_back(std::numeric_limits<float>::quiet_NaN());
        continue;
      }
      const size_t gtIdx = imTk.extra().key();
      if (gtIdx >= srcH->size()) {
        values.push_back(std::numeric_limits<float>::quiet_NaN());
        continue;
      }
      const edm::Ref<reco::TrackCollection> srcRef(srcH, gtIdx);
      values.push_back(static_cast<float>((*vmH)[srcRef].dEdx()));
    }

    auto out = std::make_unique<edm::ValueMap<float>>();
    edm::ValueMap<float>::Filler filler(*out);
    filler.insert(selH, values.begin(), values.end());
    filler.fill();
    evt.put(std::move(out));
  }

private:
  const edm::EDGetTokenT<reco::TrackCollection> selTracksToken_;
  const edm::EDGetTokenT<reco::TrackCollection> intermediateTracksToken_;
  const edm::EDGetTokenT<reco::TrackCollection> sourceTracksToken_;
  const edm::EDGetTokenT<edm::ValueMap<reco::DeDxData>> sourceValueMapToken_;
  const edm::EDGetTokenT<edm::ValueMap<unsigned int>> origIndexMapToken_;
};

DEFINE_FWK_MODULE(DeDxValueMapProjector);
