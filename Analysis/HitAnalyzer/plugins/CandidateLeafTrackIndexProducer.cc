// CandidateLeafTrackIndexProducer
//
// A flat NanoAOD tree has no EDM references, so the fact that a candidate's
// kaon *is* a particular row of the Track table is lost unless we write it as
// an integer column. This producer emits, per candidate, the row index of each
// leaf track in a target track collection, so a downstream reader can do
//   Track_dedxHarmonic2[BuJpsiK_kaonTrackIdx[i]]
// and recover the kaon's dE/dx, hits, refit quantities, etc.
//
// The index is the TrackRef key. It is only meaningful when the daughter
// TrackRefs point into the same product as `trackSrc`; the producer checks the
// ProductID and writes -1 (a sentinel) otherwise, so a cross-collection
// mismatch is visible rather than silently wrong.
//
// Leaf order follows the stage-1 layout: daughter(0) is the J/psi composite,
// so its two muons are leaves 0 and 1; the remaining leaves are bachelors.
// Emitted maps: mu0TrackIdx, mu1TrackIdx, bach0TrackIdx (and bach1TrackIdx for
// the two-bachelor channels).

#include <memory>
#include <vector>

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/ValueMap.h"
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

namespace {
void collectLeaves(const reco::Candidate* c,
                   std::vector<const reco::RecoChargedCandidate*>& out) {
  if (c == nullptr) return;
  const auto* rcc = dynamic_cast<const reco::RecoChargedCandidate*>(c);
  if (rcc != nullptr) {
    out.push_back(rcc);
    return;
  }
  for (size_t i = 0; i < c->numberOfDaughters(); ++i)
    collectLeaves(c->daughter(i), out);
}
}  // namespace

class CandidateLeafTrackIndexProducer : public edm::stream::EDProducer<> {
public:
  explicit CandidateLeafTrackIndexProducer(const edm::ParameterSet& cfg)
      : candToken_(consumes<reco::VertexCompositeCandidateCollection>(
            cfg.getParameter<edm::InputTag>("src"))),
        trackToken_(consumes<reco::TrackCollection>(
            cfg.getParameter<edm::InputTag>("trackSrc"))) {
    for (const char* n : {"mu0TrackIdx", "mu1TrackIdx", "bach0TrackIdx", "bach1TrackIdx"})
      puts_[n] = produces<edm::ValueMap<int>>(n);
  }

  void produce(edm::Event& iEvent, const edm::EventSetup&) override {
    edm::Handle<reco::VertexCompositeCandidateCollection> candH;
    iEvent.getByToken(candToken_, candH);
    edm::Handle<reco::TrackCollection> trackH;
    iEvent.getByToken(trackToken_, trackH);
    const edm::ProductID trackId = trackH.id();

    const size_t n = candH->size();
    std::vector<int> mu0(n, -1), mu1(n, -1), b0(n, -1), b1(n, -1);

    for (size_t ic = 0; ic < n; ++ic) {
      std::vector<const reco::RecoChargedCandidate*> leaves;
      collectLeaves(&(*candH)[ic], leaves);
      auto key = [&](size_t il) -> int {
        if (il >= leaves.size() || leaves[il]->track().isNull()) return -1;
        const reco::TrackRef r = leaves[il]->track();
        return (r.id() == trackId) ? static_cast<int>(r.key()) : -1;
      };
      mu0[ic] = key(0);
      mu1[ic] = key(1);
      b0[ic] = key(2);
      b1[ic] = key(3);
    }

    auto put = [&](const char* name, const std::vector<int>& v) {
      auto out = std::make_unique<edm::ValueMap<int>>();
      edm::ValueMap<int>::Filler f(*out);
      f.insert(candH, v.begin(), v.end());
      f.fill();
      iEvent.put(puts_[name], std::move(out));
    };
    put("mu0TrackIdx", mu0);
    put("mu1TrackIdx", mu1);
    put("bach0TrackIdx", b0);
    put("bach1TrackIdx", b1);
  }

private:
  const edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> candToken_;
  const edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  std::map<std::string, edm::EDPutTokenT<edm::ValueMap<int>>> puts_;
};

}  // namespace ana_hitanalyzer

using ana_hitanalyzer::CandidateLeafTrackIndexProducer;
DEFINE_FWK_MODULE(CandidateLeafTrackIndexProducer);
