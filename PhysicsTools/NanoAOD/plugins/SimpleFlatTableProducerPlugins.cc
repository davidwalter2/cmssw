#include "DataFormats/Candidate/interface/Candidate.h"
#include "PhysicsTools/NanoAOD/interface/SimpleFlatTableProducer.h"
typedef SimpleFlatTableProducer<reco::Candidate> SimpleCandidateFlatTableProducer;

typedef SimpleCollectionFlatTableProducer<reco::Candidate> SimpleCandidateCollectionFlatTableProducer;

#include "DataFormats/TrackReco/interface/Track.h"
typedef SimpleFlatTableProducer<reco::Track> SimpleTrackFlatTableProducer;

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
typedef SimpleFlatTableProducer<reco::SuperCluster> SimpleSuperclusterFlatTableProducer;

#include "DataFormats/JetReco/interface/PFJet.h"
typedef SimpleFlatTableProducer<reco::PFJet> SimplePFJetFlatTableProducer;

#include "DataFormats/JetReco/interface/GenJet.h"
typedef SimpleFlatTableProducer<reco::GenJet> SimpleGenJetFlatTableProducer;

#include "DataFormats/VertexReco/interface/Vertex.h"
typedef SimpleFlatTableProducer<reco::Vertex> SimpleVertexFlatTableProducer;

#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
typedef SimpleFlatTableProducer<reco::VertexCompositePtrCandidate> SimpleSecondaryVertexFlatTableProducer;

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
typedef SimpleFlatTableProducer<reco::GenParticle> SimpleGenParticleFlatTableProducer;

typedef SimpleTypedExternalFlatTableProducer<reco::Candidate, reco::Candidate>
    SimpleCandidate2CandidateFlatTableProducer;

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
typedef EventSingletonSimpleFlatTableProducer<GenEventInfoProduct> SimpleGenEventFlatTableProducer;

#include "SimDataFormats/GeneratorProducts/interface/GenFilterInfo.h"
typedef LumiSingletonSimpleFlatTableProducer<GenFilterInfo> SimpleGenFilterFlatTableProducerLumi;

#include "SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h"
typedef EventSingletonSimpleFlatTableProducer<HTXS::HiggsClassification> SimpleHTXSFlatTableProducer;

#include "DataFormats/ProtonReco/interface/ForwardProton.h"
typedef SimpleFlatTableProducer<reco::ForwardProton> SimpleProtonTrackFlatTableProducer;

#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"
typedef SimpleFlatTableProducer<CTPPSLocalTrackLite> SimpleLocalTrackFlatTableProducer;

#include "DataFormats/Math/interface/Point3D.h"
typedef EventSingletonSimpleFlatTableProducer<math::XYZPointF> SimpleXYZPointFlatTableProducer;

#include "DataFormats/OnlineMetaData/interface/OnlineLuminosityRecord.h"
typedef EventSingletonSimpleFlatTableProducer<OnlineLuminosityRecord> SimpleOnlineLuminosityFlatTableProducer;

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
typedef EventSingletonSimpleFlatTableProducer<reco::BeamSpot> SimpleBeamspotFlatTableProducer;

#include "DataFormats/TrackReco/interface/Track.h"
typedef SimpleFlatTableProducer<TrajectorySeed> SimpleTrajectorySeedFlatTableProducer;

#include "DataFormats/MuonSeed/interface/L2MuonTrajectorySeed.h"
typedef SimpleFlatTableProducer<L2MuonTrajectorySeed> SimpleL2MuonTrajectorySeedFlatTableProducer;

#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
typedef SimpleFlatTableProducer<reco::Track> SimpleTriggerTrackFlatTableProducer;

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
typedef SimpleFlatTableProducer<reco::GsfTrack> SimpleGsfTrackFlatTableProducer;

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
typedef SimpleFlatTableProducer<pat::CompositeCandidate> SimpleCompositeCandidateFlatTableProducer;

// RECO-tier muons. The stock Muon table is built from pat::Muon, but AlCaReco
// streams persist plain reco::Muon (no PAT embedding), and the generic
// reco::Candidate producer reaches only the kinematics, not the muon-specific
// accessors (isGlobalMuon / isTrackerMuon / numberOfMatches).
#include "DataFormats/MuonReco/interface/Muon.h"
typedef SimpleFlatTableProducer<reco::Muon> SimpleMuonFlatTableProducer;

// Detector conditions carried by the AlCaReco (scalersRawToDigi): the magnet
// current is what the CVH calibration cares about. One row per DcsStatus.
#include "DataFormats/Scalers/interface/DcsStatus.h"
typedef SimpleFlatTableProducer<DcsStatus> SimpleDcsStatusFlatTableProducer;

// RECO-tier composite candidates (V0s, quarkonia, B candidates). Needed
// because string expressions on a reco::Candidate view cannot reach daughters
// ("method daughter returned void"); the concrete type can.
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
typedef SimpleFlatTableProducer<reco::VertexCompositeCandidate>
    SimpleVertexCompositeCandidateFlatTableProducer;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SimpleCandidateFlatTableProducer);
DEFINE_FWK_MODULE(SimpleCandidateCollectionFlatTableProducer);
DEFINE_FWK_MODULE(SimpleTrackFlatTableProducer);
DEFINE_FWK_MODULE(SimpleSuperclusterFlatTableProducer);
DEFINE_FWK_MODULE(SimplePFJetFlatTableProducer);
DEFINE_FWK_MODULE(SimpleGenJetFlatTableProducer);
DEFINE_FWK_MODULE(SimpleVertexFlatTableProducer);
DEFINE_FWK_MODULE(SimpleSecondaryVertexFlatTableProducer);
DEFINE_FWK_MODULE(SimpleGenParticleFlatTableProducer);
DEFINE_FWK_MODULE(SimpleCandidate2CandidateFlatTableProducer);
DEFINE_FWK_MODULE(SimpleGenEventFlatTableProducer);
DEFINE_FWK_MODULE(SimpleGenFilterFlatTableProducerLumi);
DEFINE_FWK_MODULE(SimpleHTXSFlatTableProducer);
DEFINE_FWK_MODULE(SimpleProtonTrackFlatTableProducer);
DEFINE_FWK_MODULE(SimpleLocalTrackFlatTableProducer);
DEFINE_FWK_MODULE(SimpleXYZPointFlatTableProducer);
DEFINE_FWK_MODULE(SimpleOnlineLuminosityFlatTableProducer);
DEFINE_FWK_MODULE(SimpleBeamspotFlatTableProducer);
DEFINE_FWK_MODULE(SimpleTrajectorySeedFlatTableProducer);
DEFINE_FWK_MODULE(SimpleL2MuonTrajectorySeedFlatTableProducer);
DEFINE_FWK_MODULE(SimpleTriggerTrackFlatTableProducer);
DEFINE_FWK_MODULE(SimpleGsfTrackFlatTableProducer);
DEFINE_FWK_MODULE(SimpleCompositeCandidateFlatTableProducer);
DEFINE_FWK_MODULE(SimpleMuonFlatTableProducer);
DEFINE_FWK_MODULE(SimpleDcsStatusFlatTableProducer);
DEFINE_FWK_MODULE(SimpleVertexCompositeCandidateFlatTableProducer);
