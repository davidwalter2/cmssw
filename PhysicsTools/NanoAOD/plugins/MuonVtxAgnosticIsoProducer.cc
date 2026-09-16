#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

namespace deltas {
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

class MuonVtxAgnosticIsoProducer : public edm::stream::EDProducer<>
{
public:
  explicit MuonVtxAgnosticIsoProducer(const edm::ParameterSet &);
  ~MuonVtxAgnosticIsoProducer() {}

private:

  virtual void produce(edm::Event &, const edm::EventSetup &) override;

  using map_t = edm::ValueMap<float>;

  float maxdr_;
  float mindrchg_;
  float mindrrest_;
  float maxdeltaz_;
  float minptchg_;
  float minptneu_;
  float minptpho_;
  float minptpu_;

  int neutrpho_;

  edm::EDGetTokenT<std::vector<pat::Muon> > muonsToken_;
  edm::EDGetTokenT<std::vector<pat::PackedCandidate> > packedpfcandidatesToken_;
  edm::EDPutTokenT<map_t> outputTotalIso_;
  edm::EDPutTokenT<map_t> outputChargedHadronIso_;
  edm::EDPutTokenT<map_t> outputNeutralHadronIso_;
  edm::EDPutTokenT<map_t> outputPhotonIso_;
  edm::EDPutTokenT<map_t> outputPUIso_;
};


MuonVtxAgnosticIsoProducer::MuonVtxAgnosticIsoProducer(const edm::ParameterSet &iConfig)

{
  maxdr_ = iConfig.getParameter<double>("maxdr");
  mindrchg_ = iConfig.getParameter<double>("mindrchg");
  mindrrest_ = iConfig.getParameter<double>("mindrrest");
  maxdeltaz_ = iConfig.getParameter<double>("maxdeltaz");
  minptchg_ = iConfig.getParameter<double>("minptchg");
  minptneu_ = iConfig.getParameter<double>("minptneu");
  minptpho_ = iConfig.getParameter<double>("minptpho");
  minptpu_ = iConfig.getParameter<double>("minptpu");
  neutrpho_ = iConfig.getParameter<int>("calculateNeutralPhoton");
  muonsToken_ = consumes<std::vector<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muonInputTag"));
  packedpfcandidatesToken_ = consumes<std::vector<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("pfCandidateInputTag"));

  outputTotalIso_ = produces<map_t>("vtxAgnosticTotalIso");
  outputChargedHadronIso_ = produces<map_t>("vtxAgnosticChargedHadronIso");
  outputNeutralHadronIso_ = produces<map_t>("vtxAgnosticNeutralHadronIso");
  outputPhotonIso_ = produces<map_t>("vtxAgnosticPhotonIso");
  outputPUIso_ = produces<map_t>("vtxAgnosticPUIso");
}

void MuonVtxAgnosticIsoProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  edm::Handle<std::vector<pat::Muon> > muons;
  iEvent.getByToken(muonsToken_, muons);

  edm::Handle<std::vector<pat::PackedCandidate> > packedpfcandidates;
  iEvent.getByToken(packedpfcandidatesToken_, packedpfcandidates);

  std::vector<float> charged, neutral, photon, pu, total;
  charged.reserve(muons->size());
  neutral.reserve(muons->size());
  photon.reserve(muons->size());
  pu.reserve(muons->size());
  total.reserve(muons->size());
  for (auto const &muon : *muons) {
    float Charged=0., Neutral=0., Photon=0., PU=0., Total=0.;
    for (auto const &packedpfcandidate : *packedpfcandidates) {
      auto pdgid=packedpfcandidate.pdgId();
      if ((abs(pdgid)==11)||(abs(pdgid)==13)) continue;
      if (deltas::deltaR2(muon.eta(),muon.phi(),packedpfcandidate.eta(),packedpfcandidate.phi())>(maxdr_*maxdr_)) continue;
      if (abs(packedpfcandidate.charge())>0) {
        if (abs(muon.vz()-packedpfcandidate.vz()) < maxdeltaz_) {
          if (deltas::deltaR2(muon.eta(),muon.phi(),packedpfcandidate.eta(),packedpfcandidate.phi())<(mindrchg_*mindrchg_)) continue;
          if (packedpfcandidate.pt()>minptchg_) Charged+=packedpfcandidate.pt();
        }
        else {
          if (deltas::deltaR2(muon.eta(),muon.phi(),packedpfcandidate.eta(),packedpfcandidate.phi())<(mindrrest_*mindrrest_)) continue;
          if (packedpfcandidate.pt()>minptpu_) PU+=packedpfcandidate.pt();
        }
      }
      else if (neutrpho_!=0) {
        if (deltas::deltaR2(muon.eta(),muon.phi(),packedpfcandidate.eta(),packedpfcandidate.phi())<(mindrrest_*mindrrest_)) continue;
        if (pdgid==22) {if (packedpfcandidate.et()>minptpho_) Photon+=packedpfcandidate.et();}
        else {if (packedpfcandidate.et()>minptneu_) Neutral+=packedpfcandidate.et();}
      }
    }
    Charged/=muon.pt();
    PU/=muon.pt();
    if (neutrpho_!=0) {
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

DEFINE_FWK_MODULE(MuonVtxAgnosticIsoProducer);
