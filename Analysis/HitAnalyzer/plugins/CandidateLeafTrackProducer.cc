// CandidateLeafTrackProducer
//
// Emits, per candidate, the tracks of a chosen subset of leaf daughters, plus
// two parallel index vectors so a single-track CVH maker's per-track outputs
// can be keyed back to the candidate AND to the leaf position the downstream
// kinematic fit uses.
//
// `mode` selects the subset:
//   "bachelor" (default) : the DIRECT leaf daughters of the candidate (the
//                          bachelors -- e.g. the B+ -> [J/psi] K kaon). Fit with
//                          the bachelor mass hypothesis by the single-track maker.
//   "jpsi"               : the leaf tracks of daughter(0) (the J/psi composite,
//                          which stage-1 always puts first) -- i.e. the two
//                          muons. Fit with the muon hypothesis.
//
// The decomposition rule (shared with the two-track maker): a candidate's
// COMPOSITE daughters are two-track subsystems; its DIRECT LEAF daughters are
// single tracks. Generic across channels:
//   B+ -> [J/psi] K            : composite {J/psi}, leaf {K}     -> bachelor {K}, jpsi {mu,mu}
//   Bc -> [J/psi] pi           : composite {J/psi}, leaf {pi}    -> bachelor {pi}
//   B0 -> [J/psi] [K*0]        : composite {J/psi, K*0}, leaf {} -> bachelor {}
//
// Emitted products:
//   reco::TrackCollection (unlabeled) : the selected leaf tracks (value copies,
//                                       preserving the TrackExtraRef so the
//                                       single-track maker reaches the hits).
//   std::vector<int> "candIdx"        : parent candidate row of each track.
//   std::vector<int> "leafIdx"        : the track's position in the candidate's
//                                       canonical leaf ordering (depth-first,
//                                       daughter(0) first) -- the SAME ordering
//                                       JpsiXKinematicFitProducer uses, so the
//                                       fit can map (candIdx, leafIdx) -> refit
//                                       track with no positional assumptions.

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
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"

namespace ana_hitanalyzer {

namespace {
// Leaf tracks of a candidate node, depth-first. Identical to the ordering
// JpsiXKinematicFitProducer / CandidateLeafTrackIndexProducer use, so leafIdx
// is a shared key across the three producers.
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

class CandidateLeafTrackProducer : public edm::stream::EDProducer<> {
public:
  explicit CandidateLeafTrackProducer(const edm::ParameterSet& cfg)
      : candToken_(consumes<reco::VertexCompositeCandidateCollection>(
            cfg.getParameter<edm::InputTag>("src"))),
        mode_(cfg.existsAs<std::string>("mode") ? cfg.getParameter<std::string>("mode")
                                                : std::string("bachelor")) {
    if (mode_ != "bachelor" && mode_ != "jpsi")
      throw cms::Exception("Configuration")
          << "CandidateLeafTrackProducer mode must be bachelor|jpsi, got '" << mode_ << "'";
    tracksPut_ = produces<reco::TrackCollection>();
    candIdxPut_ = produces<std::vector<int>>("candIdx");
    leafIdxPut_ = produces<std::vector<int>>("leafIdx");
  }

  void produce(edm::Event& iEvent, const edm::EventSetup&) override {
    edm::Handle<reco::VertexCompositeCandidateCollection> candH;
    iEvent.getByToken(candToken_, candH);

    auto tracks = std::make_unique<reco::TrackCollection>();
    auto candIdx = std::make_unique<std::vector<int>>();
    auto leafIdx = std::make_unique<std::vector<int>>();

    const bool wantJpsi = (mode_ == "jpsi");

    for (size_t ic = 0; ic < candH->size(); ++ic) {
      const auto& cand = (*candH)[ic];

      // Canonical leaf ordering -> the global leafIdx of each selected leaf.
      std::vector<const reco::RecoChargedCandidate*> full;
      collectLeaves(&cand, full);

      // Membership sets (pointer identity into `full`).
      //  jpsi     : leaves under daughter(0) (the J/psi composite).
      //  bachelor : the DIRECT leaf daughters of the candidate.
      std::vector<const reco::RecoChargedCandidate*> jpsiLeaves;
      if (cand.numberOfDaughters() > 0)
        collectLeaves(cand.daughter(0), jpsiLeaves);

      for (size_t g = 0; g < full.size(); ++g) {
        const reco::RecoChargedCandidate* leaf = full[g];
        if (leaf->track().isNull()) continue;

        bool selected = false;
        if (wantJpsi) {
          for (const auto* j : jpsiLeaves)
            if (j == leaf) { selected = true; break; }
        } else {
          // bachelor: leaf is a DIRECT daughter of the candidate.
          for (size_t id = 0; id < cand.numberOfDaughters(); ++id)
            if (cand.daughter(id) == static_cast<const reco::Candidate*>(leaf)) {
              selected = true;
              break;
            }
        }
        if (!selected) continue;

        tracks->push_back(*leaf->track());  // value copy keeps the TrackExtraRef
        candIdx->push_back(static_cast<int>(ic));
        leafIdx->push_back(static_cast<int>(g));
      }
    }

    iEvent.put(tracksPut_, std::move(tracks));
    iEvent.put(candIdxPut_, std::move(candIdx));
    iEvent.put(leafIdxPut_, std::move(leafIdx));
  }

private:
  const edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> candToken_;
  const std::string mode_;
  edm::EDPutTokenT<reco::TrackCollection> tracksPut_;
  edm::EDPutTokenT<std::vector<int>> candIdxPut_;
  edm::EDPutTokenT<std::vector<int>> leafIdxPut_;
};

}  // namespace ana_hitanalyzer

using ana_hitanalyzer::CandidateLeafTrackProducer;
DEFINE_FWK_MODULE(CandidateLeafTrackProducer);
