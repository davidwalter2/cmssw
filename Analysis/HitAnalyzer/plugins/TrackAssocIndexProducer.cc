// TrackAssocIndexProducer
//
// Turns a persisted edm::Association<C> (track -> muon, or track -> PV) into a
// flat-tree-friendly edm::ValueMap<int> giving, per track, the row index of the
// associated object in C (or -1 when there is no association). This is the
// Track->Muon / Track->PV analogue of CandidateLeafTrackIndexProducer: a flat
// NanoAOD carries no EDM refs, so the association is only preservable as an
// integer column, e.g. Muon_<x>[Track_muonIdx[i]].
//
// Templated on the association target so one implementation serves both links;
// concrete typedefs + DEFINE_FWK_MODULE are at the bottom.
//
// originalIndex bridge (Track->PV, task 4.6c)
// -------------------------------------------
// The persisted track->PV Association is keyed to `generalTracks`, NOT to the
// cloned ALCARECO alignment tracks that feed the NanoAOD Track table. So the
// direct `assocH->contains(trackH.id())` lookup below fails and every index
// stays -1. The AlCaReco keeps a `ValueMap<uint> originalIndex` mapping each
// alignment track to its `generalTracks` row precisely to close this gap. When
// `originalIndex` and `pvSrc` are configured, any track left at -1 by the direct
// lookup is resolved by: originalIndex[i] -> the generalTracks key -> the PV
// whose track list contains that key (built once per event from the persisted
// vertices' own track refs). Left empty (the Track->Muon case, whose
// association IS keyed to the alignment tracks), the bridge is inert.

#include <memory>
#include <unordered_map>
#include <vector>

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

namespace ana_hitanalyzer {

template <typename TColl>
class TrackAssocIndexProducer : public edm::stream::EDProducer<> {
public:
  explicit TrackAssocIndexProducer(const edm::ParameterSet& cfg)
      : trackToken_(consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("trackSrc"))),
        assocToken_(consumes<edm::Association<TColl>>(cfg.getParameter<edm::InputTag>("association"))),
        origIndexTag_(cfg.getParameter<edm::InputTag>("originalIndex")),
        pvSrcTag_(cfg.getParameter<edm::InputTag>("pvSrc")),
        useBridge_(!origIndexTag_.label().empty() && !pvSrcTag_.label().empty()) {
    if (useBridge_) {
      origIndexToken_ = consumes<edm::ValueMap<unsigned int>>(origIndexTag_);
      pvToken_ = consumes<reco::VertexCollection>(pvSrcTag_);
    }
    put_ = produces<edm::ValueMap<int>>();
  }

  void produce(edm::Event& iEvent, const edm::EventSetup&) override {
    edm::Handle<reco::TrackCollection> trackH;
    iEvent.getByToken(trackToken_, trackH);
    edm::Handle<edm::Association<TColl>> assocH;
    iEvent.getByToken(assocToken_, assocH);

    std::vector<int> idx(trackH->size(), -1);
    // The association may be keyed to a different track product than trackSrc
    // (e.g. the offlinePrimaryVertices track->PV map is keyed to generalTracks,
    // not the persisted alignment tracks). operator[] would throw an
    // InvalidReference in that case, so only look up when the association
    // actually contains our product; otherwise the originalIndex bridge below
    // (Track->PV) closes the gap and, absent that, every index stays -1.
    if (assocH->contains(trackH.id())) {
      for (size_t i = 0; i < trackH->size(); ++i) {
        const reco::TrackRef tref(trackH, i);
        const edm::Ref<TColl> target = (*assocH)[tref];
        if (target.isNonnull()) idx[i] = static_cast<int>(target.key());
      }
    } else if (useBridge_) {
      bridgeViaOriginalIndex(iEvent, trackH, idx);
    }

    auto out = std::make_unique<edm::ValueMap<int>>();
    typename edm::ValueMap<int>::Filler f(*out);
    f.insert(trackH, idx.begin(), idx.end());
    f.fill();
    iEvent.put(put_, std::move(out));
  }

private:
  // Track->PV via the originalIndex bridge. Compiled for every TColl but only
  // meaningful (and only wired) for the vertex case; the map from a
  // generalTracks key to its PV row is read from the persisted vertices' own
  // track references, which carry the generalTracks key even though the
  // generalTracks collection itself is not in the AlCaReco.
  void bridgeViaOriginalIndex(edm::Event& iEvent, const edm::Handle<reco::TrackCollection>& trackH,
                              std::vector<int>& idx) {
    edm::Handle<edm::ValueMap<unsigned int>> origH;
    iEvent.getByToken(origIndexToken_, origH);
    edm::Handle<reco::VertexCollection> pvH;
    iEvent.getByToken(pvToken_, pvH);
    if (!origH.isValid() || !pvH.isValid()) return;

    // generalTracks key -> PV row, from each vertex's stored track refs.
    std::unordered_map<unsigned int, int> keyToPV;
    for (size_t iv = 0; iv < pvH->size(); ++iv) {
      const reco::Vertex& v = (*pvH)[iv];
      for (auto it = v.tracks_begin(); it != v.tracks_end(); ++it) {
        // TrackBaseRef into generalTracks; key() is the generalTracks row.
        keyToPV.emplace(static_cast<unsigned int>(it->key()), static_cast<int>(iv));
      }
    }
    if (keyToPV.empty()) return;  // vertices carry no track refs -> nothing to bridge

    for (size_t i = 0; i < trackH->size(); ++i) {
      if (idx[i] >= 0) continue;
      const unsigned int gkey = (*origH)[reco::TrackRef(trackH, i)];
      auto found = keyToPV.find(gkey);
      if (found != keyToPV.end()) idx[i] = found->second;
    }
  }

  const edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  const edm::EDGetTokenT<edm::Association<TColl>> assocToken_;
  const edm::InputTag origIndexTag_;
  const edm::InputTag pvSrcTag_;
  const bool useBridge_;
  edm::EDGetTokenT<edm::ValueMap<unsigned int>> origIndexToken_;
  edm::EDGetTokenT<reco::VertexCollection> pvToken_;
  edm::EDPutTokenT<edm::ValueMap<int>> put_;
};

using TrackToMuonIndexProducer = TrackAssocIndexProducer<reco::MuonCollection>;
using TrackToVertexIndexProducer = TrackAssocIndexProducer<reco::VertexCollection>;

}  // namespace ana_hitanalyzer

using ana_hitanalyzer::TrackToMuonIndexProducer;
using ana_hitanalyzer::TrackToVertexIndexProducer;
DEFINE_FWK_MODULE(TrackToMuonIndexProducer);
DEFINE_FWK_MODULE(TrackToVertexIndexProducer);
