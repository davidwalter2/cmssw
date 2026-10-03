// JpsiKCandidateSplitter
//
// Stage-2 Stage-1-output adapter for B+ -> J/psi K+ (and other one-track-
// bachelor channels). Consumes a `VertexCompositeCandidateCollection` of B+
// candidates (each with two daughters: a J/psi VCC and a bachelor
// RecoChargedCandidate) and emits, in the same per-B+ order:
//
//   - a `reco::VertexCompositeCandidateCollection` of dimuon (J/psi) sub-
//     candidates -- shallow copies of daughter(0). Their internal muon-RCC
//     daughters keep pointing at the persisted track collection, so the
//     downstream `ResidualGlobalCorrectionMakerTwoTrackG4e` (`srcCandidates=
//     <this>`, `src=ALCARECOTkAlJpsiX`) resolves tracks via the candidate's
//     fast path (RecoChargedCandidate::track()).
//   - a `reco::TrackCollection` of bachelor tracks -- value-copies of
//     daughter(1)::track() dereferenced. The Track value copy preserves its
//     internal TrackExtraRef into the persisted vector<reco::TrackExtra>, so
//     the single-track maker (`src=<this>`) reaches hits the usual way.
//   - two parallel `std::vector<int>` index branches, in the SAME order as the
//     emitted dimuon / bachelor collections, holding the source B+ candidate
//     index. Per-candidate joining downstream uses these.
//
// Species disambiguation is by the daughter `pdgId()` set at Stage 1 by
// `JpsiXCandidateProducer` (bachelor) / `TwoBodyDecayCandidateProducer`
// (J/psi muons and K*0/phi daughters). This splitter reads the bachelor
// pdgId and the two J/psi-daughter muon pdgIds off the candidate and
// propagates them onto its output as parallel `std::vector<int>` branches
// (`bachelorPdgId`, `muon0PdgId`, `muon1PdgId`). Downstream Stage-2
// mass-hypothesis lookup consumes those branches instead of hard-coding
// per-channel masses.
//
// Collection LAYOUT still defines the daughter ORDER
// (daughter(0) = J/psi VCC, daughter(1) = bachelor RCC): only the species
// TAG moves to the pdgId. See openspec change
// add-jpsi-x-muons-and-preprod-refinements.
//
// No fit code, no candidate-level cuts. Stage-1 already selected.

#include <memory>
#include <vector>

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

namespace ana_hitanalyzer {

class JpsiKCandidateSplitter : public edm::stream::EDProducer<> {
public:
  explicit JpsiKCandidateSplitter(const edm::ParameterSet& iConfig);
  ~JpsiKCandidateSplitter() override = default;

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> bplusToken_;
  // Diagnostic counters (one per produce() call). Total across the job is
  // assembled by the outer harness when needed; we just emit LogDebug here.
  // Not thread-shared; one instance per stream.
};

JpsiKCandidateSplitter::JpsiKCandidateSplitter(const edm::ParameterSet& iConfig)
    : bplusToken_(consumes<reco::VertexCompositeCandidateCollection>(
          iConfig.getParameter<edm::InputTag>("src"))) {
  produces<reco::VertexCompositeCandidateCollection>("dimuon");
  produces<reco::TrackCollection>("bachelor");
  produces<std::vector<int>>("dimuonBCandIdx");
  produces<std::vector<int>>("bachelorBCandIdx");
  // Species tags read from daughter->pdgId() at Stage 1 (openspec change
  // add-jpsi-x-muons-and-preprod-refinements). Downstream Stage-2 lookup
  // consumes these instead of hard-coding per-channel masses.
  produces<std::vector<int>>("bachelorPdgId");
  produces<std::vector<int>>("muon0PdgId");
  produces<std::vector<int>>("muon1PdgId");
}

void JpsiKCandidateSplitter::produce(edm::Event& iEvent, const edm::EventSetup&) {
  edm::Handle<reco::VertexCompositeCandidateCollection> bplusH;
  iEvent.getByToken(bplusToken_, bplusH);

  auto outDimuon = std::make_unique<reco::VertexCompositeCandidateCollection>();
  auto outBachelor = std::make_unique<reco::TrackCollection>();
  auto outDimuonIdx = std::make_unique<std::vector<int>>();
  auto outBachelorIdx = std::make_unique<std::vector<int>>();
  auto outBachelorPdg = std::make_unique<std::vector<int>>();
  auto outMuon0Pdg = std::make_unique<std::vector<int>>();
  auto outMuon1Pdg = std::make_unique<std::vector<int>>();

  int nSkip = 0;
  outDimuon->reserve(bplusH->size());
  outBachelor->reserve(bplusH->size());
  outDimuonIdx->reserve(bplusH->size());
  outBachelorIdx->reserve(bplusH->size());
  outBachelorPdg->reserve(bplusH->size());
  outMuon0Pdg->reserve(bplusH->size());
  outMuon1Pdg->reserve(bplusH->size());

  for (std::size_t ib = 0; ib < bplusH->size(); ++ib) {
    const reco::VertexCompositeCandidate& bplus = (*bplusH)[ib];
    if (bplus.numberOfDaughters() < 2) {
      ++nSkip;
      continue;
    }
    // daughter(0) is the J/psi VCC, daughter(1) is the bachelor RCC, per
    // `JpsiXCandidateProducer.cc`. LAYOUT is trusted; species/mass hypothesis
    // is looked up from daughter->pdgId() (Stage-1 non-zero invariant,
    // openspec add-jpsi-x-muons-and-preprod-refinements).
    const auto* dJpsi = dynamic_cast<const reco::VertexCompositeCandidate*>(
        bplus.daughter(0));
    const auto* dBach = dynamic_cast<const reco::RecoChargedCandidate*>(
        bplus.daughter(1));
    if (!dJpsi || !dBach) {
      ++nSkip;
      continue;
    }
    if (dJpsi->numberOfDaughters() != 2) {
      ++nSkip;
      continue;
    }
    if (dBach->track().isNull()) {
      ++nSkip;
      continue;
    }

    outDimuon->push_back(*dJpsi);
    outDimuonIdx->push_back(static_cast<int>(ib));

    outBachelor->push_back(*dBach->track());
    outBachelorIdx->push_back(static_cast<int>(ib));

    outBachelorPdg->push_back(dBach->pdgId());
    outMuon0Pdg->push_back(dJpsi->daughter(0)->pdgId());
    outMuon1Pdg->push_back(dJpsi->daughter(1)->pdgId());
  }

  if (nSkip > 0) {
    LogDebug("JpsiKCandidateSplitter")
        << "skipped " << nSkip << " B+ candidates with unexpected daughter layout"
        << " out of " << bplusH->size() << " total";
  }

  iEvent.put(std::move(outDimuon), "dimuon");
  iEvent.put(std::move(outBachelor), "bachelor");
  iEvent.put(std::move(outDimuonIdx), "dimuonBCandIdx");
  iEvent.put(std::move(outBachelorIdx), "bachelorBCandIdx");
  iEvent.put(std::move(outBachelorPdg), "bachelorPdgId");
  iEvent.put(std::move(outMuon0Pdg), "muon0PdgId");
  iEvent.put(std::move(outMuon1Pdg), "muon1PdgId");
}

}  // namespace ana_hitanalyzer

using ana_hitanalyzer::JpsiKCandidateSplitter;
DEFINE_FWK_MODULE(JpsiKCandidateSplitter);
