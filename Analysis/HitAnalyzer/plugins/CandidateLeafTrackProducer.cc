// CandidateLeafTrackProducer
//
// Generic replacement for JpsiKCandidateSplitter's bachelor output. Emits the
// tracks of the DIRECT leaf daughters of each candidate (the bachelors), plus a
// parallel candidate-index vector so the single-track CVH maker can key its
// per-candidate outputs the same way it did off the splitter's bCandIdx.
//
// The decomposition rule (shared with the two-track maker): a candidate's
// COMPOSITE daughters are two-track subsystems (fit by the two-track maker); its
// DIRECT LEAF daughters are single tracks (fit by the single-track maker off
// this collection). It is generic across channels:
//   B+ -> [J/psi] K            : composite {J/psi}, leaf {K}        -> 1 bachelor
//   Bc -> [J/psi] pi           : composite {J/psi}, leaf {pi}       -> 1 bachelor
//   B0 -> [J/psi] [K*0]        : composite {J/psi, K*0}, leaf {}    -> 0 bachelors
// so a two-composite channel simply yields no single tracks, which is correct.
//
// Track value-copies preserve their TrackExtraRef, so the single-track maker
// reaches hits exactly as it did with the splitter's output.

#include <memory>
#include <vector>

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

namespace ana_hitanalyzer {

class CandidateLeafTrackProducer : public edm::stream::EDProducer<> {
public:
  explicit CandidateLeafTrackProducer(const edm::ParameterSet& cfg)
      : candToken_(consumes<reco::VertexCompositeCandidateCollection>(
            cfg.getParameter<edm::InputTag>("src"))) {
    tracksPut_ = produces<reco::TrackCollection>();
    idxPut_ = produces<std::vector<int>>("candIdx");
  }

  void produce(edm::Event& iEvent, const edm::EventSetup&) override {
    edm::Handle<reco::VertexCompositeCandidateCollection> candH;
    iEvent.getByToken(candToken_, candH);

    auto tracks = std::make_unique<reco::TrackCollection>();
    auto candIdx = std::make_unique<std::vector<int>>();

    for (size_t ic = 0; ic < candH->size(); ++ic) {
      const auto& cand = (*candH)[ic];
      for (size_t id = 0; id < cand.numberOfDaughters(); ++id) {
        const auto* rcc = dynamic_cast<const reco::RecoChargedCandidate*>(cand.daughter(id));
        // Only DIRECT leaf daughters -- composite daughters are the two-track
        // maker's job. A null track ref is skipped.
        if (rcc == nullptr || rcc->track().isNull()) continue;
        tracks->push_back(*rcc->track());  // value copy keeps the TrackExtraRef
        candIdx->push_back(static_cast<int>(ic));
      }
    }

    iEvent.put(tracksPut_, std::move(tracks));
    iEvent.put(idxPut_, std::move(candIdx));
  }

private:
  const edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> candToken_;
  edm::EDPutTokenT<reco::TrackCollection> tracksPut_;
  edm::EDPutTokenT<std::vector<int>> idxPut_;
};

}  // namespace ana_hitanalyzer

using ana_hitanalyzer::CandidateLeafTrackProducer;
DEFINE_FWK_MODULE(CandidateLeafTrackProducer);
