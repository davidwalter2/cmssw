// TrackImpactParameterProducer
//
// Emits, per track, the impact parameters measured with respect to a primary
// vertex, as ValueMap<float>: d0 (transverse), dzPV (longitudinal), and an
// uncertainty for each.
//
// Why this exists
// ---------------
// The NanoAOD Track table's `dxy`/`dz` come from reco::Track's own accessors,
// which are measured from the ORIGIN (0,0,0). Profiling `dxy` against track phi
// gives a clean sinusoid whose amplitude equals the median transverse beam
// offset -- 0.117 cm in 2016H data, 0.190 cm in this simulation -- so the
// stored quantity is unusable as an impact parameter, and the apparent data/MC
// disagreement is only the two beamspots sitting at different positions.
//
// reco::Track::dxy(const Point&) does the right thing, but a
// SimpleTrackFlatTableProducer `Var` string evaluates member functions on the
// track alone and cannot reach the vertex collection, and a search of
// CommonTools and PhysicsTools found no existing ValueMap producer for this.
// Hence a producer.
//
// Vertex choice
// -------------
// Each track's own associated PV is used where one exists, via the same
// ValueMap<int> that fills `Track_pvIdx` (produced by TrackAssocIndexProducer,
// which already handles the originalIndex bridge onto the cloned alignment
// tracks). Tracks with no association fall back to the leading vertex; this is
// reported per event at LogDebug so a systematic failure of the association
// cannot hide as a slightly-wrong impact parameter.
//
// Uncertainty
// -----------
// d0Err combines the track's own dxyError with the vertex's transverse position
// uncertainty projected onto the track's impact-parameter direction, which for a
// track of azimuth phi is (-sin phi, cos phi):
//
//   var(d0) = dxyError^2 + sin^2(phi) Vxx - 2 sin(phi) cos(phi) Vxy + cos^2(phi) Vyy
//
// dzPVErr adds dzError and the vertex z uncertainty in quadrature. Both are
// approximations in that the track is not refit without the vertex, so for a
// track used in the vertex fit the uncertainty is mildly optimistic. That is the
// standard treatment and is stated here rather than hidden.

#include <cmath>
#include <memory>
#include <vector>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

namespace ana_hitanalyzer {

class TrackImpactParameterProducer : public edm::stream::EDProducer<> {
public:
  explicit TrackImpactParameterProducer(const edm::ParameterSet& cfg)
      : trackToken_(consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("trackSrc"))),
        pvToken_(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("pvSrc"))),
        pvIdxTag_(cfg.getParameter<edm::InputTag>("pvIdx")),
        usePvIdx_(!pvIdxTag_.label().empty()),
        sentinel_(cfg.getParameter<double>("sentinel")) {
    if (usePvIdx_)
      pvIdxToken_ = consumes<edm::ValueMap<int>>(pvIdxTag_);
    produces<edm::ValueMap<float>>("d0");
    produces<edm::ValueMap<float>>("dzPV");
    produces<edm::ValueMap<float>>("d0Err");
    produces<edm::ValueMap<float>>("dzPVErr");
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("trackSrc");
    desc.add<edm::InputTag>("pvSrc", edm::InputTag("offlinePrimaryVertices"));
    desc.add<edm::InputTag>("pvIdx", edm::InputTag(""))
        ->setComment("ValueMap<int> of per-track PV row (-1 = none); empty = always use the leading vertex");
    desc.add<double>("sentinel", -99.)
        ->setComment("value written when no vertex is available at all");
    descriptions.addDefault(desc);
  }

  void produce(edm::Event& evt, const edm::EventSetup&) override {
    edm::Handle<reco::TrackCollection> tracks;
    evt.getByToken(trackToken_, tracks);
    edm::Handle<reco::VertexCollection> vertices;
    evt.getByToken(pvToken_, vertices);

    edm::Handle<edm::ValueMap<int>> pvIdx;
    const bool havePvIdx = usePvIdx_ && evt.getByToken(pvIdxToken_, pvIdx) && pvIdx.isValid();

    const size_t n = tracks->size();
    std::vector<float> d0(n, static_cast<float>(sentinel_));
    std::vector<float> dzPV(n, static_cast<float>(sentinel_));
    std::vector<float> d0Err(n, static_cast<float>(sentinel_));
    std::vector<float> dzPVErr(n, static_cast<float>(sentinel_));

    unsigned int nFallback = 0;
    if (!vertices->empty()) {
      for (size_t i = 0; i < n; ++i) {
        const reco::Track& trk = (*tracks)[i];

        int iv = -1;
        if (havePvIdx) {
          const edm::Ref<reco::TrackCollection> ref(tracks, i);
          iv = (*pvIdx)[ref];
        }
        if (iv < 0 || static_cast<size_t>(iv) >= vertices->size()) {
          iv = 0;  // leading vertex
          ++nFallback;
        }
        const reco::Vertex& pv = (*vertices)[iv];

        d0[i] = static_cast<float>(trk.dxy(pv.position()));
        dzPV[i] = static_cast<float>(trk.dz(pv.position()));

        // Project the vertex transverse covariance onto the impact-parameter
        // direction (-sin phi, cos phi) and add the track term in quadrature.
        const double s = std::sin(trk.phi());
        const double c = std::cos(trk.phi());
        const double vtxVarXY =
            s * s * pv.covariance(0, 0) - 2. * s * c * pv.covariance(0, 1) + c * c * pv.covariance(1, 1);
        const double varD0 = trk.dxyError() * trk.dxyError() + std::max(0., vtxVarXY);
        const double varDz = trk.dzError() * trk.dzError() + std::max(0., pv.covariance(2, 2));
        d0Err[i] = static_cast<float>(std::sqrt(varD0));
        dzPVErr[i] = static_cast<float>(std::sqrt(varDz));
      }
    }

    if (nFallback > 0)
      LogDebug("TrackImpactParameterProducer")
          << nFallback << " of " << n << " tracks had no PV association and used the leading vertex";

    put(evt, tracks, d0, "d0");
    put(evt, tracks, dzPV, "dzPV");
    put(evt, tracks, d0Err, "d0Err");
    put(evt, tracks, dzPVErr, "dzPVErr");
  }

private:
  static void put(edm::Event& evt,
                  const edm::Handle<reco::TrackCollection>& tracks,
                  const std::vector<float>& vals,
                  const std::string& label) {
    auto out = std::make_unique<edm::ValueMap<float>>();
    edm::ValueMap<float>::Filler filler(*out);
    filler.insert(tracks, vals.begin(), vals.end());
    filler.fill();
    evt.put(std::move(out), label);
  }

  const edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  const edm::EDGetTokenT<reco::VertexCollection> pvToken_;
  const edm::InputTag pvIdxTag_;
  edm::EDGetTokenT<edm::ValueMap<int>> pvIdxToken_;
  const bool usePvIdx_;
  const double sentinel_;
};

}  // namespace ana_hitanalyzer

using ana_hitanalyzer::TrackImpactParameterProducer;
DEFINE_FWK_MODULE(TrackImpactParameterProducer);
