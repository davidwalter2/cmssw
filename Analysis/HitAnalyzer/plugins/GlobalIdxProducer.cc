#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"

class GlobalIdxProducer : public edm::stream::EDProducer<>
{
public:
  explicit GlobalIdxProducer(const edm::ParameterSet &);
  ~GlobalIdxProducer() {}

//   static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:

  virtual void produce(edm::Event &, const edm::EventSetup &) override;

  using map_t = edm::ValueMap<std::vector<int>>;

  edm::EDGetTokenT<map_t> inputGlobalIdxs0_;
  edm::EDGetTokenT<map_t> inputGlobalIdxs1_;
  edm::EDPutTokenT<map_t> outputGlobalIdxs_;


};


GlobalIdxProducer::GlobalIdxProducer(const edm::ParameterSet &iConfig)

{
  inputGlobalIdxs0_ = consumes<map_t>(iConfig.getParameter<edm::InputTag>("src0"));
  inputGlobalIdxs1_ = consumes<map_t>(iConfig.getParameter<edm::InputTag>("src1"));

  outputGlobalIdxs_ = produces<map_t>();
}

// ------------ method called for each event  ------------
void GlobalIdxProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  using namespace edm;

  Handle<map_t> idxs0;
  iEvent.getByToken(inputGlobalIdxs0_, idxs0);

  Handle<map_t> idxs1;
  iEvent.getByToken(inputGlobalIdxs1_, idxs1);


  // take idxs from idxs0, unless idxs0 is empty, then take them from idxs1
  map_t idxsout = *idxs0;

  for (auto it = idxsout.begin(); it != idxsout.end(); ++it) {
    auto const id = it.id();
    auto const size = it.size();
    if (!idxs1->contains(id)) {
      throw std::runtime_error("product ids should always match");
    }
    for (std::size_t i = 0; i < size; ++i) {
      auto &valout = idxsout.get(id, i);
      auto const &val1 = idxs1->get(id, i);

      if (!valout.empty() && ! val1.empty() && valout != val1) {
        throw std::runtime_error("idxs should always be either empty or equal;");
      }

      if (valout.empty() && !val1.empty()) {
        valout = val1;
      }
    }
  }

  iEvent.emplace(outputGlobalIdxs_, std::move(idxsout));

}

DEFINE_FWK_MODULE(GlobalIdxProducer);
