// DiMuonTrackVertexCandidateProducer
//
// Builds opposite-sign muon-track pairs from a reco::TrackCollection (the
// muon inner tracks made by TrackProducerFromPatMuons) into a
// reco::VertexCompositeCandidateCollection, in a configurable invariant-mass
// window. Each candidate carries two reco::RecoChargedCandidate daughters that
// wrap the paired reco::Tracks (via TrackRef), which is exactly the form the
// CVH two-track refit (ResidualGlobalCorrectionMakerTwoTrackG4e) consumes on
// its srcCandidates path -- it re-does the vertex/mass fit, so the candidate
// here is only a pairing + selection vehicle (the vertex/mass stored are
// placeholders from a trivial 4-vector sum).
//
// Resonance-agnostic: pick the mass window to select Z / J/psi / Upsilon; the
// downstream refit and NanoAOD table are identical for all of them.

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

#include <memory>
#include <vector>

class DiMuonTrackVertexCandidateProducer : public edm::global::EDProducer<> {
public:
  explicit DiMuonTrackVertexCandidateProducer(const edm::ParameterSet& iConfig)
      : tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("src"))),
        massMin_(iConfig.getParameter<double>("massMin")),
        massMax_(iConfig.getParameter<double>("massMax")),
        oppositeSign_(iConfig.getParameter<bool>("oppositeSign")) {
    produces<reco::VertexCompositeCandidateCollection>();
  }

  void produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup&) const override {
    edm::Handle<reco::TrackCollection> tracksH;
    iEvent.getByToken(tracksToken_, tracksH);

    auto out = std::make_unique<reco::VertexCompositeCandidateCollection>();
    AddFourMomenta addp4;

    const std::size_t n = tracksH->size();
    for (std::size_t i = 0; i < n; ++i) {
      const reco::Track& ti = (*tracksH)[i];
      for (std::size_t j = i + 1; j < n; ++j) {
        const reco::Track& tj = (*tracksH)[j];
        if (oppositeSign_ && (ti.charge() + tj.charge() != 0))
          continue;

        const math::XYZTLorentzVector p4i(ti.px(), ti.py(), ti.pz(),
                                          std::sqrt(ti.p() * ti.p() + kMuonMass2));
        const math::XYZTLorentzVector p4j(tj.px(), tj.py(), tj.pz(),
                                          std::sqrt(tj.p() * tj.p() + kMuonMass2));
        const double mass = (p4i + p4j).mass();
        if (mass < massMin_ || mass > massMax_)
          continue;

        reco::RecoChargedCandidate di;
        di.setCharge(ti.charge());
        di.setP4(p4i);
        di.setVertex(ti.vertex());
        di.setTrack(reco::TrackRef(tracksH, i));

        reco::RecoChargedCandidate dj;
        dj.setCharge(tj.charge());
        dj.setP4(p4j);
        dj.setVertex(tj.vertex());
        dj.setTrack(reco::TrackRef(tracksH, j));

        reco::VertexCompositeCandidate cand;
        cand.addDaughter(di);
        cand.addDaughter(dj);
        addp4.set(cand);  // sets p4 + charge from the daughters
        out->push_back(cand);
      }
    }

    iEvent.put(std::move(out));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("tracksfrommuons"))
        ->setComment("muon inner-track collection (e.g. TrackProducerFromPatMuons output)");
    desc.add<double>("massMin", 50.0)->setComment("min mu-mu invariant mass (GeV)");
    desc.add<double>("massMax", 150.0)->setComment("max mu-mu invariant mass (GeV)");
    desc.add<bool>("oppositeSign", true);
    descriptions.addWithDefaultLabel(desc);
  }

private:
  static constexpr double kMuonMass = 0.1056583745;
  static constexpr double kMuonMass2 = kMuonMass * kMuonMass;

  const edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  const double massMin_;
  const double massMax_;
  const bool oppositeSign_;
};

DEFINE_FWK_MODULE(DiMuonTrackVertexCandidateProducer);
