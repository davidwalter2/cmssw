// TrackPatImpactParameterProducer
//
// The impact parameters a pat::Muon carries (PATMuonProducer::embedHighLevel),
// computed for a bare reco::Track collection and emitted as ValueMap<float>,
// one per NanoAOD Muon column they fill:
//
//   dxy, dxyErr      transverse IP w.r.t. the leading primary vertex (PV2D)
//   dz, dzErr        longitudinal IP w.r.t. the leading primary vertex (PVDZ)
//   dxybs, dxybsErr  transverse IP w.r.t. the beam spot (BS2D)
//   ip3d, sip3d      |3D IP| w.r.t. the leading primary vertex and its
//                    significance (PV3D)
//
// The AlCaReco NanoAOD builds its Muon table on the AlCaReco muon tracks (the
// AlCaReco carries no reco::Muon); this producer gives those columns the
// definitions of the MiniAOD NanoAOD. The formulas, the choice of the leading
// vertex and the behaviour without a vertex (a default-constructed one, with
// an ip3d error of -1) are PATMuonProducer's.

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

class TrackPatImpactParameterProducer : public edm::stream::EDProducer<> {
public:
  explicit TrackPatImpactParameterProducer(const edm::ParameterSet& cfg)
      : tracksToken_(consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("src"))),
        pvToken_(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("pvSrc"))),
        bsToken_(consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamSpot"))),
        ttbToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))) {
    for (const char* name : kNames)
      produces<edm::ValueMap<float>>(name);
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src")->setComment("reco::Track collection");
    desc.add<edm::InputTag>("pvSrc", edm::InputTag("offlinePrimaryVertices"))
        ->setComment("primary vertices; the leading one is the reference, as for pat::Muon");
    desc.add<edm::InputTag>("beamSpot", edm::InputTag("offlineBeamSpot"));
    descriptions.addWithDefaultLabel(desc);
  }

  void produce(edm::Event& evt, const edm::EventSetup& es) override {
    const auto tracks = evt.getHandle(tracksToken_);
    const auto& vertices = evt.get(pvToken_);
    const auto& beamSpot = evt.get(bsToken_);
    const auto& builder = es.getData(ttbToken_);

    const bool pvValid = !vertices.empty();
    const reco::Vertex pv = pvValid ? vertices.front() : reco::Vertex();

    const std::size_t n = tracks->size();
    std::vector<std::vector<float>> vals(std::size(kNames), std::vector<float>(n));
    for (std::size_t i = 0; i < n; ++i) {
      const reco::TrackRef track(tracks, i);
      const reco::TransientTrack tt = builder.build(track);
      const GlobalVector dir(track->px(), track->py(), track->pz());
      const Measurement1D ip3d = IPTools::signedImpactParameter3D(tt, dir, pv).second;
      const double ip3dErr = pvValid ? ip3d.error() : -1.;

      vals[0][i] = track->dxy(pv.position());
      vals[1][i] = track->dxyError(pv.position(), pv.covariance());
      vals[2][i] = track->dz(pv.position());
      vals[3][i] = std::hypot(track->dzError(), pv.zError());
      vals[4][i] = track->dxy(beamSpot);
      vals[5][i] = track->dxyError(beamSpot);
      vals[6][i] = std::abs(ip3d.value());
      vals[7][i] = std::abs(ip3d.value() / ip3dErr);
    }

    for (std::size_t k = 0; k < std::size(kNames); ++k) {
      auto out = std::make_unique<edm::ValueMap<float>>();
      edm::ValueMap<float>::Filler filler(*out);
      filler.insert(tracks, vals[k].begin(), vals[k].end());
      filler.fill();
      evt.put(std::move(out), kNames[k]);
    }
  }

private:
  static constexpr const char* kNames[] = {"dxy", "dxyErr", "dz", "dzErr", "dxybs", "dxybsErr", "ip3d", "sip3d"};

  const edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  const edm::EDGetTokenT<reco::VertexCollection> pvToken_;
  const edm::EDGetTokenT<reco::BeamSpot> bsToken_;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbToken_;
};

DEFINE_FWK_MODULE(TrackPatImpactParameterProducer);
