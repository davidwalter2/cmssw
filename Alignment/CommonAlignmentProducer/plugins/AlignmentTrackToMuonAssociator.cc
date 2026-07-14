/// \class AlignmentTrackToMuonAssociator
///
/// EDProducer emitting an edm::Association<reco::MuonCollection> keyed on a
/// deduplicated cloned track collection (e.g. `ALCARECOTkAlJpsiX`) and
/// valued into the persisted muon collection (e.g.
/// `ALCARECOTkAlJpsiXLooseMuons`). For each cloned track:
///   - if the track is the innerTrack of some muon in the persisted muon
///     collection, the association returns a MuonRef to that muon;
///   - otherwise the association returns a null MuonRef.
///
/// The chain from a cloned track back to a `generalTracks` index follows
/// the same pattern as `DeDxValueMapProjector`:
///
///   selH[j] --(originalIndexMap[selH[j]])--> intermediateTracks[ii]
///          --(intermediateTracks[ii].extra().key())--> generalTracks index
///
/// The muons' `innerTrack()` refs point into `generalTracks`, so we build a
/// reverse map (generalTracks index -> muon index) once per event and
/// look up each cloned track through it.
///
/// openspec change add-jpsi-x-muons-and-preprod-refinements.

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include <limits>
#include <memory>
#include <unordered_map>
#include <vector>

class AlignmentTrackToMuonAssociator : public edm::global::EDProducer<> {
public:
  explicit AlignmentTrackToMuonAssociator(const edm::ParameterSet& cfg)
      : selTracksToken_(consumes<reco::TrackCollection>(
            cfg.getParameter<edm::InputTag>("selectedTracks"))),
        intermediateTracksToken_(consumes<reco::TrackCollection>(
            cfg.getParameter<edm::InputTag>("intermediateTracks"))),
        origIndexMapToken_(consumes<edm::ValueMap<unsigned int>>(
            cfg.getParameter<edm::InputTag>("originalIndexMap"))),
        muonsToken_(consumes<reco::MuonCollection>(
            cfg.getParameter<edm::InputTag>("muons"))) {
    produces<edm::Association<reco::MuonCollection>>();
  }

  void produce(edm::StreamID, edm::Event& evt, const edm::EventSetup&) const override {
    edm::Handle<reco::TrackCollection> selH;
    evt.getByToken(selTracksToken_, selH);
    edm::Handle<reco::TrackCollection> intermH;
    evt.getByToken(intermediateTracksToken_, intermH);
    edm::Handle<edm::ValueMap<unsigned int>> origIdxH;
    evt.getByToken(origIndexMapToken_, origIdxH);
    edm::Handle<reco::MuonCollection> muonsH;
    evt.getByToken(muonsToken_, muonsH);

    // Reverse lookup: generalTracks index -> muon index.
    std::unordered_map<std::size_t, std::size_t> gtIdxToMuonIdx;
    gtIdxToMuonIdx.reserve(muonsH->size() * 2);
    for (std::size_t iMu = 0; iMu < muonsH->size(); ++iMu) {
      const auto& mu = (*muonsH)[iMu];
      const auto& tr = mu.innerTrack();
      if (tr.isNonnull()) {
        gtIdxToMuonIdx.emplace(tr.key(), iMu);
      }
    }

    constexpr unsigned int kMissing = std::numeric_limits<unsigned int>::max();
    std::vector<int> muonIndex(selH->size(), -1);
    for (std::size_t j = 0; j < selH->size(); ++j) {
      const edm::Ref<reco::TrackCollection> selRef(selH, j);
      const unsigned int intermIdx = (*origIdxH)[selRef];
      if (intermIdx == kMissing || intermIdx >= intermH->size()) continue;
      const auto& imTk = (*intermH)[intermIdx];
      if (imTk.extra().isNull()) continue;
      const std::size_t gtIdx = imTk.extra().key();
      auto it = gtIdxToMuonIdx.find(gtIdx);
      if (it != gtIdxToMuonIdx.end()) {
        muonIndex[j] = static_cast<int>(it->second);
      }
    }

    auto out = std::make_unique<edm::Association<reco::MuonCollection>>(muonsH);
    edm::Association<reco::MuonCollection>::Filler filler(*out);
    filler.insert(selH, muonIndex.begin(), muonIndex.end());
    filler.fill();
    evt.put(std::move(out));
  }

private:
  const edm::EDGetTokenT<reco::TrackCollection> selTracksToken_;
  const edm::EDGetTokenT<reco::TrackCollection> intermediateTracksToken_;
  const edm::EDGetTokenT<edm::ValueMap<unsigned int>> origIndexMapToken_;
  const edm::EDGetTokenT<reco::MuonCollection> muonsToken_;
};

DEFINE_FWK_MODULE(AlignmentTrackToMuonAssociator);
