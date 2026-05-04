/// \class DeDxValueMapProjector
///
/// EDProducer that re-keys a ValueMap<float> from the original generalTracks
/// collection onto a deep-copied selected-track collection (e.g. the
/// ALCARECOTkAl{Ks,Lambda,Dst}* outputs).
///
/// The selected tracks are deep copies of generalTracks (made by
/// V0DaughterTrackProducer + AlignmentTrackSelectorModule), but each Track
/// retains its TrackExtraRef, whose key() is the index into
/// generalTracks::TrackExtras. By the standard CMSSW reco convention the
/// extras and the tracks are produced in parallel, so that key is also the
/// index into generalTracks itself. We use it to look up the original
/// ValueMap value and emit a new ValueMap keyed on the selected tracks.
///
/// Used by the V0 ALCARECO skims to persist per-track dE/dx (dedxHarmonic2,
/// dedxPixelHarmonic2) for downstream particle-ID without keeping the full
/// generalTracks collection.

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

#include <memory>
#include <vector>

/// Reads the standard `ValueMap<reco::DeDxData>` (as produced by
/// `dedxHarmonic2` / `dedxPixelHarmonic2`) keyed onto generalTracks, and
/// emits a `ValueMap<float>` (just the `.dEdx()` scalar) keyed onto the
/// selected/cloned track collection. The keying trick relies on the cloned
/// Track preserving its TrackExtraRef.key() into the source collection.
class DeDxValueMapProjector : public edm::global::EDProducer<> {
public:
  explicit DeDxValueMapProjector(const edm::ParameterSet& cfg)
      : selTracksToken_(consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("selectedTracks"))),
        srcTracksToken_(consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("sourceTracks"))),
        srcValueMapToken_(consumes<edm::ValueMap<reco::DeDxData>>(cfg.getParameter<edm::InputTag>("sourceValueMap"))) {
    produces<edm::ValueMap<float>>();
  }

  void produce(edm::StreamID, edm::Event& evt, const edm::EventSetup&) const override {
    edm::Handle<reco::TrackCollection> selH;
    evt.getByToken(selTracksToken_, selH);
    edm::Handle<reco::TrackCollection> srcH;
    evt.getByToken(srcTracksToken_, srcH);
    edm::Handle<edm::ValueMap<reco::DeDxData>> vmH;
    evt.getByToken(srcValueMapToken_, vmH);

    std::vector<float> values;
    values.reserve(selH->size());
    for (const auto& tk : *selH) {
      const size_t origIdx = tk.extra().key();
      const edm::Ref<reco::TrackCollection> srcRef(srcH, origIdx);
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
  const edm::EDGetTokenT<reco::TrackCollection> srcTracksToken_;
  const edm::EDGetTokenT<edm::ValueMap<reco::DeDxData>> srcValueMapToken_;
};

DEFINE_FWK_MODULE(DeDxValueMapProjector);
