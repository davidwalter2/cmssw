// TrackBeamspotConstraintValueMapProducer
//
// MuonBeamspotConstraintValueMapProducer (RecoMuon/GlobalTrackingTools, the
// NanoAOD Muon_bsConstrained* columns) for a bare reco::Track collection: the
// same beam-spot constraint of each track (the primary vertex when the beam
// spot widths are badly measured, the unconstrained track when both fail), the
// same ValueMap names. The AlCaReco NanoAOD builds its Muon table on the
// AlCaReco muon tracks, which play the role of muonBestTrack().

#include <memory>
#include <vector>

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "RecoVertex/KalmanVertexFit/interface/SingleTrackVertexConstraint.h"
#include "RecoVertex/VertexPrimitives/interface/VertexException.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

class TrackBeamspotConstraintValueMapProducer : public edm::global::EDProducer<> {
public:
  explicit TrackBeamspotConstraintValueMapProducer(const edm::ParameterSet& config)
      : trackToken_(consumes<reco::TrackCollection>(config.getParameter<edm::InputTag>("src"))),
        beamSpotToken_(consumes<reco::BeamSpot>(config.getParameter<edm::InputTag>("beamspot"))),
        primaryVertexToken_(consumes<reco::VertexCollection>(config.getParameter<edm::InputTag>("vertices"))),
        ttbToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))) {
    produces<edm::ValueMap<float>>("muonBSConstrainedPt");
    produces<edm::ValueMap<float>>("muonBSConstrainedPtErr");
    produces<edm::ValueMap<float>>("muonBSConstrainedChi2");
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src")->setComment("reco::Track collection");
    desc.add<edm::InputTag>("beamspot", edm::InputTag("offlineBeamSpot"))->setComment("Beam spot collection");
    desc.add<edm::InputTag>("vertices", edm::InputTag("offlinePrimaryVertices"))
        ->setComment("Primary vertex collection");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  void produce(edm::StreamID, edm::Event& event, const edm::EventSetup& setup) const override {
    const auto tracks = event.getHandle(trackToken_);
    const auto beamSpotHandle = event.getHandle(beamSpotToken_);
    const auto& ttkb = setup.getData(ttbToken_);

    std::vector<float> pts, ptErrs, chi2s;
    pts.reserve(tracks->size());
    ptErrs.reserve(tracks->size());
    chi2s.reserve(tracks->size());

    for (std::size_t i = 0; i < tracks->size(); ++i) {
      const reco::TrackRef track(tracks, i);
      bool tbd = true;
      if (beamSpotHandle.isValid()) {
        // Protect for mis-reconstructed beamspots (SingleTrackVertexConstraint
        // uses the width for the constraint, not the error)
        if ((beamSpotHandle->BeamWidthXError() / beamSpotHandle->BeamWidthX() < 0.3) &&
            (beamSpotHandle->BeamWidthYError() / beamSpotHandle->BeamWidthY() < 0.3)) {
          try {
            SingleTrackVertexConstraint::BTFtuple btft = stvc_.constrain(ttkb.build(track), *beamSpotHandle);
            if (std::get<0>(btft)) {
              const reco::Track& trkBS = std::get<1>(btft).track();
              pts.push_back(trkBS.pt());
              ptErrs.push_back(trkBS.ptError());
              chi2s.push_back(std::get<2>(btft));
              tbd = false;
            }
          } catch (const VertexException& exc) {
            // Update failed; give up.
          }
        }
      }
      if (tbd) {
        // Invalid BS; use PV instead
        const auto pvHandle = event.getHandle(primaryVertexToken_);
        if (pvHandle.isValid() && !pvHandle->empty()) {
          const auto& pv = pvHandle->at(0);
          VertexState pvs = VertexState(GlobalPoint(Basic3DVector<float>(pv.position())), GlobalError(pv.covariance()));
          try {
            SingleTrackVertexConstraint::BTFtuple btft = stvc_.constrain(ttkb.build(track), pvs);
            if (std::get<0>(btft)) {
              const reco::Track& trkBS = std::get<1>(btft).track();
              pts.push_back(trkBS.pt());
              ptErrs.push_back(trkBS.ptError());
              chi2s.push_back(std::get<2>(btft));
              tbd = false;
            }
          } catch (const VertexException& exc) {
            // Update failed; give up.
          }
        }
      }
      if (tbd) {
        // Fall-back case, keep the unconstrained values
        pts.push_back(track->pt());
        ptErrs.push_back(track->ptError());
        chi2s.push_back(-1.f);
      }
    }

    put(event, tracks, pts, "muonBSConstrainedPt");
    put(event, tracks, ptErrs, "muonBSConstrainedPtErr");
    put(event, tracks, chi2s, "muonBSConstrainedChi2");
  }

  static void put(edm::Event& event,
                  const edm::Handle<reco::TrackCollection>& tracks,
                  const std::vector<float>& vals,
                  const char* label) {
    auto valueMap = std::make_unique<edm::ValueMap<float>>();
    edm::ValueMap<float>::Filler filler(*valueMap);
    filler.insert(tracks, vals.begin(), vals.end());
    filler.fill();
    event.put(std::move(valueMap), label);
  }

  const edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  const edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  const edm::EDGetTokenT<reco::VertexCollection> primaryVertexToken_;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbToken_;
  const SingleTrackVertexConstraint stvc_;
};

DEFINE_FWK_MODULE(TrackBeamspotConstraintValueMapProducer);
