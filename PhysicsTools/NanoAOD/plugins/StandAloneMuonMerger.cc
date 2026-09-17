#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

class StandAloneMuonMerger : public edm::stream::EDProducer<>
{
public:
  explicit StandAloneMuonMerger(const edm::ParameterSet &);
  ~StandAloneMuonMerger() {}

//   static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  
  virtual void produce(edm::Event &, const edm::EventSetup &) override;
  
  edm::EDGetTokenT<reco::TrackCollection> inputTrackStandAlone_;
  edm::EDGetTokenT<reco::TrackCollection> inputTrackStandAloneUpdatedAtVtx_;
  edm::EDGetTokenT<std::vector<pat::Muon>> inputMuons_;
  edm::EDPutTokenT<reco::TrackCollection> outputTrack_;
  edm::EDPutTokenT<edm::ValueMap<bool>> outputUpdatedAtVtxMap_;
  edm::EDPutTokenT<edm::ValueMap<bool>> outputMuonUpdatedAtVtxMap_;

};


StandAloneMuonMerger::StandAloneMuonMerger(const edm::ParameterSet &iConfig)

{
  inputTrackStandAlone_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("standAlone"));
  inputTrackStandAloneUpdatedAtVtx_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("standAloneUpdatedAtVtx"));
  inputMuons_ = consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"));
  
  outputTrack_ = produces<reco::TrackCollection>();
  outputUpdatedAtVtxMap_ = produces<edm::ValueMap<bool>>("updatedAtVtx");
  outputMuonUpdatedAtVtxMap_ = produces<edm::ValueMap<bool>>("muonUpdatedAtVtx");
}

// ------------ method called for each event  ------------
void StandAloneMuonMerger::produce(edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  using namespace edm;
  
  Handle<reco::TrackCollection> tracksStandAlone;
  iEvent.getByToken(inputTrackStandAlone_, tracksStandAlone);

  Handle<reco::TrackCollection> tracksStandAloneUpdatedAtVtx;
  iEvent.getByToken(inputTrackStandAloneUpdatedAtVtx_, tracksStandAloneUpdatedAtVtx);
  
  Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(inputMuons_, muons);
  
  reco::TrackCollection tracksOut;
  
  // fill collection of disambiguated standalone muon tracks plus value map indicating whether
  // the selected track is updated at vertex or not
  
  std::vector<bool> updatedAtVtxV;
  updatedAtVtxV.reserve(tracksStandAlone->size());
  
  for (auto const &standAlone : *tracksStandAlone) {
    auto const extraIdx = standAlone.extra().key();
    
    const reco::Track *track = &standAlone;
    bool updatedAtVtx = false;
    
    // replicate logic in GlobalMuonProducer, take the updatedAtVtx track if it exists and has
    // the same eta sign as the original, otherwise take the original
    for (auto const &standAloneUpdatedAtVtx : *tracksStandAloneUpdatedAtVtx) {
      if (standAloneUpdatedAtVtx.extra().key() == extraIdx) {
        const int etaFlip1 = (standAloneUpdatedAtVtx.eta() * standAlone.eta()) < 0 ? -1 : 1;
        if (etaFlip1==1) {
          track = &standAloneUpdatedAtVtx;
          updatedAtVtx = true;
        }
        break;
      }
    }
    
    tracksOut.emplace_back(*track);
    updatedAtVtxV.push_back(updatedAtVtx);
  }
    
  // fill value map for the muon objects indicating whether standalone muon track is updated at vertex
  
  std::vector<bool> muonUpdatedAtVtxV;
  muonUpdatedAtVtxV.reserve(muons->size());
  
  for (auto const &muon : *muons) {
    bool isUpdatedAtVertex = false;
      
    // can't compare references because of track embedding in pat::Muon, so match by track extra index and check kinematics instead
    if (muon.standAloneMuon().isNonnull()) {
      auto const &muonStandAlone = *muon.standAloneMuon();
      
      for (auto const &updatedStandAlone : *tracksStandAloneUpdatedAtVtx) {
        if (updatedStandAlone.extra().key() == muonStandAlone.extra().key()) {
          if (muonStandAlone.pt() == updatedStandAlone.pt() && muonStandAlone.eta() == updatedStandAlone.eta() && muonStandAlone.phi() == updatedStandAlone.phi()) {
            isUpdatedAtVertex = true;
          }
          break;
        }
      }
    }
    
    muonUpdatedAtVtxV.push_back(isUpdatedAtVertex);
  }
  
  edm::ValueMap<bool> muonUpdatedAtVtxMap;
  
  edm::ValueMap<bool>::Filler muonUpdatedAtVtxMapFiller(muonUpdatedAtVtxMap);
  muonUpdatedAtVtxMapFiller.insert(muons, std::make_move_iterator(muonUpdatedAtVtxV.begin()), std::make_move_iterator(muonUpdatedAtVtxV.end()));
  muonUpdatedAtVtxMapFiller.fill();
  
  
  auto outputTrackHandle = iEvent.emplace(outputTrack_, std::move(tracksOut));
  
  // finish filling valuemap for output tracks
  edm::ValueMap<bool> updatedAtVtxMap;
  
  edm::ValueMap<bool>::Filler updatedAtVtxMapFiller(updatedAtVtxMap);
  updatedAtVtxMapFiller.insert(outputTrackHandle, std::make_move_iterator(updatedAtVtxV.begin()), std::make_move_iterator(updatedAtVtxV.end()));
  updatedAtVtxMapFiller.fill();
  
  iEvent.emplace(outputUpdatedAtVtxMap_, std::move(updatedAtVtxMap));
  
  iEvent.emplace(outputMuonUpdatedAtVtxMap_, std::move(muonUpdatedAtVtxMap));

  
}

DEFINE_FWK_MODULE(StandAloneMuonMerger);
