/// \class VertexCompositeCandidateRemapper
///
/// EDProducer that emits a copy of an input
/// reco::VertexCompositeCandidateCollection where every leaf
/// reco::RecoChargedCandidate daughter has its TrackRef rewritten to point
/// at a *filtered* TrackCollection (typically the ALCAREco output of an
/// AlignmentTrackSelector). Candidates whose daughters did not all survive
/// the upstream track selection are dropped.
///
/// Reference chain to map a candidate daughter (`Ref<TrackCollection>` into
/// the original generalTracks-style collection) to a selected-track index:
///
///   daughter.track().key()  ==  intermediateTracks[i].extra().key()
///                           |     for some i in intermediateTracks
///                           |
///                           ==  originalIndexMap[selH[j]]  for some j
///
/// where:
///  * `originalIndexMap` is the side-channel ValueMap<unsigned int> emitted
///    by AlignmentTrackSelectorWithIndexMapModule, recording each cloned
///    track's index in its source collection;
///  * `intermediateTracks` is that source collection -- typically the
///    V0DaughterTrackProducer output (whose tracks preserve the original
///    generalTracks key in `extra().key()`), or generalTracks itself when
///    AlignmentTrackSelector takes generalTracks directly.
///
/// We invert the chain once per event into
///   `genTracksKey -> selH index`
/// and walk the candidates, rebuilding leaves with TrackRefs into the
/// filtered collection.
///
/// Recurses into nested composite daughters (e.g. for B -> J/psi(->mumu) K),
/// so the same module handles flat 2/3-body candidates and nested
/// topologies. No kinematic fingerprinting anywhere -- everything is
/// reference-based.

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include <limits>
#include <memory>
#include <unordered_map>

class VertexCompositeCandidateRemapper : public edm::stream::EDProducer<> {
public:
  explicit VertexCompositeCandidateRemapper(const edm::ParameterSet& cfg)
      : candToken_(consumes<reco::VertexCompositeCandidateCollection>(
            cfg.getParameter<edm::InputTag>("srcCandidates"))),
        selectedTracksToken_(
            consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("selectedTracks"))),
        intermediateTracksToken_(
            consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("intermediateTracks"))),
        origIndexMapToken_(
            consumes<edm::ValueMap<unsigned int>>(cfg.getParameter<edm::InputTag>("originalIndexMap"))) {
    produces<reco::VertexCompositeCandidateCollection>();
  }

  void produce(edm::Event& evt, const edm::EventSetup&) override {
    edm::Handle<reco::VertexCompositeCandidateCollection> candH;
    evt.getByToken(candToken_, candH);

    edm::Handle<reco::TrackCollection> selH;
    evt.getByToken(selectedTracksToken_, selH);

    edm::Handle<reco::TrackCollection> intermH;
    evt.getByToken(intermediateTracksToken_, intermH);

    edm::Handle<edm::ValueMap<unsigned int>> mapH;
    evt.getByToken(origIndexMapToken_, mapH);

    // Build inverse map: generalTracks-style key -> selH index, by chaining
    //   selH[j] -> intermH[origIdx] -> intermH[origIdx].extra().key()
    constexpr unsigned int kMissing = std::numeric_limits<unsigned int>::max();
    std::unordered_map<unsigned int, size_t> origToSel;
    origToSel.reserve(selH->size() * 2);
    for (size_t j = 0; j < selH->size(); ++j) {
      const edm::Ref<reco::TrackCollection> ref(selH, j);
      const unsigned int intermIdx = (*mapH)[ref];
      if (intermIdx == kMissing || intermIdx >= intermH->size()) continue;
      const auto& imTk = (*intermH)[intermIdx];
      if (imTk.extra().isNull()) continue;
      origToSel.emplace(static_cast<unsigned int>(imTk.extra().key()), j);
    }

    auto out = std::make_unique<reco::VertexCompositeCandidateCollection>();
    out->reserve(candH->size());
    for (const auto& cand : *candH) {
      reco::VertexCompositeCandidate remapped = makeShell(cand);
      if (rebuildDaughters(cand, remapped, origToSel, selH))
        out->push_back(std::move(remapped));
    }

    LogDebug("VertexCompositeCandidateRemapper")
        << "Remapped " << out->size() << " / " << candH->size()
        << " candidates onto " << selH->size() << " filtered tracks.";
    evt.put(std::move(out));
  }

private:
  // Vertex-only copy of a candidate (daughters rebuilt separately). When the
  // source is itself a VertexCompositeCandidate, carry the vertex covariance
  // and fit chi2/ndof over -- V0Producer / vertex-fitted candidates store
  // physics content there that would otherwise be silently dropped.
  static reco::VertexCompositeCandidate makeShell(const reco::Candidate& c) {
    if (const auto* vcc = dynamic_cast<const reco::VertexCompositeCandidate*>(&c)) {
      reco::Candidate::CovarianceMatrix cov;
      vcc->fillVertexCovariance(cov);
      return reco::VertexCompositeCandidate(
          c.charge(), c.p4(), c.vertex(), cov, vcc->vertexChi2(), vcc->vertexNdof(), c.pdgId());
    }
    return reco::VertexCompositeCandidate(c.charge(), c.p4(), c.vertex(), c.pdgId());
  }

  // Copy daughters of `src` into `dst`, rewriting TrackRefs of leaf
  // RecoChargedCandidates via origToSel. Returns false if any leaf is
  // missing from the filtered collection (-> caller drops the whole tree).
  static bool rebuildDaughters(const reco::Candidate& src,
                               reco::CompositeCandidate& dst,
                               const std::unordered_map<unsigned int, size_t>& origToSel,
                               const edm::Handle<reco::TrackCollection>& selH) {
    for (size_t i = 0; i < src.numberOfDaughters(); ++i) {
      const reco::Candidate* d = src.daughter(i);
      if (!d) return false;

      if (const auto* ch = dynamic_cast<const reco::RecoChargedCandidate*>(d)) {
        const auto& origRef = ch->track();
        if (origRef.isNull()) return false;
        const auto it = origToSel.find(static_cast<unsigned int>(origRef.key()));
        if (it == origToSel.end()) return false;
        reco::RecoChargedCandidate newDau(*ch);
        newDau.setTrack(reco::TrackRef(selH, it->second));
        dst.addDaughter(newDau);
      } else if (d->numberOfDaughters() > 0) {
        // Nested composite: rebuild recursively.
        reco::VertexCompositeCandidate nested = makeShell(*d);
        if (!rebuildDaughters(*d, nested, origToSel, selH)) return false;
        dst.addDaughter(nested);
      } else {
        // Leaf that is neither charged-track-bearing nor composite -- we
        // can't remap it; drop the whole candidate to stay consistent.
        return false;
      }
    }
    return true;
  }

  edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> candToken_;
  edm::EDGetTokenT<reco::TrackCollection> selectedTracksToken_;
  edm::EDGetTokenT<reco::TrackCollection> intermediateTracksToken_;
  edm::EDGetTokenT<edm::ValueMap<unsigned int>> origIndexMapToken_;
};

DEFINE_FWK_MODULE(VertexCompositeCandidateRemapper);
