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

#include <memory>
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
        assocToken_(consumes<edm::Association<TColl>>(cfg.getParameter<edm::InputTag>("association"))) {
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
    // actually contains our product; otherwise every index stays -1.
    if (assocH->contains(trackH.id())) {
      for (size_t i = 0; i < trackH->size(); ++i) {
        const reco::TrackRef tref(trackH, i);
        const edm::Ref<TColl> target = (*assocH)[tref];
        if (target.isNonnull()) idx[i] = static_cast<int>(target.key());
      }
    }

    auto out = std::make_unique<edm::ValueMap<int>>();
    typename edm::ValueMap<int>::Filler f(*out);
    f.insert(trackH, idx.begin(), idx.end());
    f.fill();
    iEvent.put(put_, std::move(out));
  }

private:
  const edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  const edm::EDGetTokenT<edm::Association<TColl>> assocToken_;
  edm::EDPutTokenT<edm::ValueMap<int>> put_;
};

using TrackToMuonIndexProducer = TrackAssocIndexProducer<reco::MuonCollection>;
using TrackToVertexIndexProducer = TrackAssocIndexProducer<reco::VertexCollection>;

}  // namespace ana_hitanalyzer

using ana_hitanalyzer::TrackToMuonIndexProducer;
using ana_hitanalyzer::TrackToVertexIndexProducer;
DEFINE_FWK_MODULE(TrackToMuonIndexProducer);
DEFINE_FWK_MODULE(TrackToVertexIndexProducer);
