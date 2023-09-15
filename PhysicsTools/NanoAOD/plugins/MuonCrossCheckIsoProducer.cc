#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include "ROOT/RVec.hxx"
#include <Math/Vector4D.h>
#include "TVector2.h"
#include "TLorentzVector.h"

int chargedHadronVertex( const reco::VertexCollection& vertices, const pat::PackedCandidate& pfcand ) {

  //reco::PFCandidate cand = pfcand;
  //auto const & track = cand.trackRef();  
  size_t  iVertex = 0;
  unsigned int index=0;
  // no vertex found with this track. 

  // optional: as a secondary solution, associate the closest vertex in z
  //if ( checkClosestZVertex_ ) {

    double dzmin = 10000;
    double ztrack = pfcand.vertex().z();
    bool foundVertex = false;
    index = 0;
    for(auto iv=vertices.begin(); iv!=vertices.end(); ++iv, ++index) {

      double dz = fabs(ztrack - iv->z());
      if(dz<dzmin) {
	dzmin = dz; 
	iVertex = index;
	foundVertex = true;
      }
    }

    if( foundVertex ) 
      return iVertex;  

  //}


  return -1 ;
}

namespace Deltas {
  double deltaPhi(float phi1, float phi2) {
    double result = phi1 - phi2;
    while (result > M_PI) result -= 2.0*M_PI;
    while (result <= -1.0*M_PI) result += 2.0*M_PI;
    return result;
  }

  double deltaR2(float eta1, float phi1, float eta2, float phi2) {
    double deta = eta1-eta2;
    double dphi = deltaPhi(phi1,phi2);
    return deta*deta + dphi*dphi;
  }
}

class MuonCrossCheckIsoProducer : public edm::stream::EDProducer<>
{
public:
  explicit MuonCrossCheckIsoProducer(const edm::ParameterSet &);
  ~MuonCrossCheckIsoProducer() {}

private:

  virtual void produce(edm::Event &, const edm::EventSetup &) override;

  using map_t = edm::ValueMap<float>;

  float maxdr_;
  float mindr_;
  float maxdeltaz_;
  float minptchg_;
  float minptneu_;
  float minptpho_;
  float minptpu_;

  int neutrpho_;

  edm::EDGetTokenT<std::vector<pat::Muon> > muonsToken_;
  edm::EDGetTokenT<std::vector<pat::PackedCandidate> > packedpfcandidatesToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex> > verticesToken_;
  edm::EDPutTokenT<map_t> outputTotalIso_;
  edm::EDPutTokenT<map_t> outputChargedHadronIso_;
  edm::EDPutTokenT<map_t> outputNeutralHadronIso_;
  edm::EDPutTokenT<map_t> outputPhotonIso_;
  edm::EDPutTokenT<map_t> outputPUIso_;
};


MuonCrossCheckIsoProducer::MuonCrossCheckIsoProducer(const edm::ParameterSet &iConfig)

{
  maxdr_ = iConfig.getParameter<double>("maxdr");
  mindr_ = iConfig.getParameter<double>("mindr");
  maxdeltaz_ = iConfig.getParameter<double>("maxdeltaz");
  minptchg_ = iConfig.getParameter<double>("minptchg");
  minptneu_ = iConfig.getParameter<double>("minptneu");
  minptpho_ = iConfig.getParameter<double>("minptpho");
  minptpu_ = iConfig.getParameter<double>("minptpu");
  neutrpho_ = iConfig.getParameter<int>("calculateNeutralPhoton");
  muonsToken_ = consumes<std::vector<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muonInputTag"));
  packedpfcandidatesToken_ = consumes<std::vector<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("pfCandidateInputTag"));
  verticesToken_ = consumes<std::vector<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("pvInputTag"));

  outputTotalIso_ = produces<map_t>("vtxAgnosticTotalIso");
  outputChargedHadronIso_ = produces<map_t>("vtxAgnosticChargedHadronIso");
  outputNeutralHadronIso_ = produces<map_t>("vtxAgnosticNeutralHadronIso");
  outputPhotonIso_ = produces<map_t>("vtxAgnosticPhotonIso");
  outputPUIso_ = produces<map_t>("vtxAgnosticPUIso");
}

void MuonCrossCheckIsoProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  edm::Handle<std::vector<pat::Muon> > muons;
  iEvent.getByToken(muonsToken_, muons);

  edm::Handle<std::vector<pat::PackedCandidate> > packedpfcandidates;
  iEvent.getByToken(packedpfcandidatesToken_, packedpfcandidates);

  edm::Handle<std::vector<reco::Vertex> > vertices;
  iEvent.getByToken(verticesToken_, vertices);

  std::vector<float> charged, neutral, photon, pu, total;
  charged.reserve(muons->size());
  neutral.reserve(muons->size());
  photon.reserve(muons->size());
  pu.reserve(muons->size());
  total.reserve(muons->size());
  std::vector<pat::PackedCandidate> PUCands, ChargedCands, NeutralCands, PhotonCands;
  for (auto const &packedpfcandidate : *packedpfcandidates) {
    auto pdgid=packedpfcandidate.pdgId();
    if ((abs(pdgid)==11)||(abs(pdgid)==13)) continue;
    if (abs(packedpfcandidate.charge())>0) {
      int ivertex = chargedHadronVertex( *vertices, packedpfcandidate );
      if ((ivertex == -1)||(ivertex == 0)) {if (packedpfcandidate.pt()>minptchg_) ChargedCands.push_back(packedpfcandidate);}
      else {if (packedpfcandidate.pt()>minptpu_) PUCands.push_back(packedpfcandidate);}
    }
    else {
      if (pdgid==22) {if (packedpfcandidate.et()>minptpho_) PhotonCands.push_back(packedpfcandidate);}
      else {if (packedpfcandidate.et()>minptneu_) NeutralCands.push_back(packedpfcandidate);}
    }
  }
  for (auto const &muon : *muons) {
    if (maxdr_>0.35) if (muon.pt()>10.) std::cout<<"Muon pt: "<<muon.pt()<<"\teta: "<<muon.eta()<<"\tphi: "<<muon.phi()<<"\n";
    float Charged=0., Neutral=0., Photon=0., PU=0., Total=0.;
    if (maxdr_>0.35) if (muon.pt()>10.) std::cout<<"reeval PU:\n";
    for (auto const &cand : PUCands) {
      if (Deltas::deltaR2(muon.eta(),muon.phi(),cand.eta(),cand.phi())>(maxdr_*maxdr_)) continue;
      if (Deltas::deltaR2(muon.eta(),muon.phi(),cand.eta(),cand.phi())<(mindr_*mindr_)) continue;
      PU+=cand.pt();
      if (maxdr_>0.35) if (muon.pt()>10.) std::cout<<"pt: "<<cand.pt()<<"\teta: "<<cand.eta()<<"\tphi: "<<cand.phi()<<"\tdeltar: "<<Deltas::deltaR2(muon.eta(),muon.phi(),cand.eta(),cand.phi())<<"\n";
    }
    if (maxdr_>0.35) if (muon.pt()>10.) std::cout<<"reeval Charged:\n";
    for (auto const &cand : ChargedCands) {
      if (Deltas::deltaR2(muon.eta(),muon.phi(),cand.eta(),cand.phi())>(maxdr_*maxdr_)) continue;
      if (Deltas::deltaR2(muon.eta(),muon.phi(),cand.eta(),cand.phi())<(mindr_*mindr_)) continue;
      Charged+=cand.pt();
      if (maxdr_>0.35) if (muon.pt()>10.) std::cout<<"pt: "<<cand.pt()<<"\teta: "<<cand.eta()<<"\tphi: "<<cand.phi()<<"\tdeltar: "<<Deltas::deltaR2(muon.eta(),muon.phi(),cand.eta(),cand.phi())<<"\n";
    }
    Charged/=muon.pt();
    PU/=muon.pt();
    if (neutrpho_!=0) {
      for (auto const &cand : NeutralCands) {
        if (Deltas::deltaR2(muon.eta(),muon.phi(),cand.eta(),cand.phi())>(maxdr_*maxdr_)) continue;
        if (Deltas::deltaR2(muon.eta(),muon.phi(),cand.eta(),cand.phi())<(mindr_*mindr_)) continue;
        Neutral+=cand.et();
      }
      for (auto const &cand : PhotonCands) {
        if (Deltas::deltaR2(muon.eta(),muon.phi(),cand.eta(),cand.phi())>(maxdr_*maxdr_)) continue;
        if (Deltas::deltaR2(muon.eta(),muon.phi(),cand.eta(),cand.phi())<(mindr_*mindr_)) continue;
        Photon+=cand.et();
      }
      Neutral/=muon.pt();
      Photon/=muon.pt();
    }
    else {
      Neutral=(maxdr_>0.35) ? muon.pfIsolationR04().sumNeutralHadronEt/muon.pt() : muon.pfIsolationR03().sumNeutralHadronEt/muon.pt();
      Photon=(maxdr_>0.35) ? muon.pfIsolationR04().sumPhotonEt/muon.pt() : muon.pfIsolationR03().sumPhotonEt/muon.pt();
    }
    Total=Charged+std::max(0.,Neutral+Photon-0.5*PU);
    charged.push_back(Charged);
    neutral.push_back(Neutral);
    photon.push_back(Photon);
    pu.push_back(PU);
    total.push_back(Total);
  }
  map_t totalMap, chargedMap, neutralMap, photonMap, puMap;
  map_t::Filler totalMapFiller(totalMap), chargedMapFiller(chargedMap), neutralMapFiller(neutralMap), photonMapFiller(photonMap), puMapFiller(puMap);
  totalMapFiller.insert(muons, std::make_move_iterator(total.begin()), std::make_move_iterator(total.end()));
  totalMapFiller.fill();
  chargedMapFiller.insert(muons, std::make_move_iterator(charged.begin()), std::make_move_iterator(charged.end()));
  chargedMapFiller.fill();
  neutralMapFiller.insert(muons, std::make_move_iterator(neutral.begin()), std::make_move_iterator(neutral.end()));
  neutralMapFiller.fill();
  photonMapFiller.insert(muons, std::make_move_iterator(photon.begin()), std::make_move_iterator(photon.end()));
  photonMapFiller.fill();
  puMapFiller.insert(muons, std::make_move_iterator(pu.begin()), std::make_move_iterator(pu.end()));
  puMapFiller.fill();
  iEvent.emplace(outputTotalIso_, std::move(totalMap));
  iEvent.emplace(outputChargedHadronIso_, std::move(chargedMap));
  iEvent.emplace(outputNeutralHadronIso_, std::move(neutralMap));
  iEvent.emplace(outputPhotonIso_, std::move(photonMap));
  iEvent.emplace(outputPUIso_, std::move(puMap));
}

DEFINE_FWK_MODULE(MuonCrossCheckIsoProducer);
