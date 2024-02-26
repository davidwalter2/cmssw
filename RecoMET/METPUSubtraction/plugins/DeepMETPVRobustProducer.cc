#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"

struct DeepMETPVRobustCache {
  std::atomic<tensorflow::GraphDef*> graph_def;
};

class DeepMETPVRobustProducer : public edm::stream::EDProducer<edm::GlobalCache<DeepMETPVRobustCache>> {
public:
  explicit DeepMETPVRobustProducer(const edm::ParameterSet&, const DeepMETPVRobustCache*);
  void produce(edm::Event& event, const edm::EventSetup& setup) override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  // static methods for handling the global cache
  static std::unique_ptr<DeepMETPVRobustCache> initializeGlobalCache(const edm::ParameterSet&);
  static void globalEndJob(DeepMETPVRobustCache*);

private:
  const edm::EDGetTokenT<std::vector<pat::PackedCandidate>> pf_token_;
  const edm::EDGetTokenT<edm::ValueMap<double>> PFPVRobustDxy_;
  const edm::EDGetTokenT<edm::ValueMap<double>> PFPVRobustDz_;
  const edm::EDGetTokenT<edm::ValueMap<double>> PFPVRobustPuppiWeight_;

  const bool usePUPPI_;
  const float norm_;
  const bool ignore_leptons_;
  const unsigned int max_n_pf_;
  const bool do_print_;

  tensorflow::Session* session_;

  tensorflow::Tensor input_;
  tensorflow::Tensor input_cat0_;
  tensorflow::Tensor input_cat1_;

  inline static const std::unordered_map<int, int32_t> charge_embedding_{{-1, 0}, {0, 1}, {1, 2}};
  inline static const std::unordered_map<int, int32_t> pdg_id_embedding_{
      {-211, 0}, {-13, 1}, {-11, 2}, {0, 3}, {1, 4}, {2, 5}, {11, 6}, {13, 7}, {22, 8}, {130, 9}, {211, 10}};
};

namespace {
  float scale_and_rm_outlier(float val, float scale) {
    float ret_val = val * scale;
    if (ret_val > 1e6 || ret_val < -1e6)
      return 0.;
    return ret_val;
  }
}  // namespace

DeepMETPVRobustProducer::DeepMETPVRobustProducer(const edm::ParameterSet& cfg, const DeepMETPVRobustCache* cache)
    : pf_token_(consumes<std::vector<pat::PackedCandidate>>(cfg.getParameter<edm::InputTag>("pf_src"))),
      PFPVRobustDxy_(consumes<edm::ValueMap<double>>(cfg.getParameter<edm::InputTag>("PFPVRobustDxy"))),
      PFPVRobustDz_(consumes<edm::ValueMap<double>>(cfg.getParameter<edm::InputTag>("PFPVRobustDz"))),
      PFPVRobustPuppiWeight_(consumes<edm::ValueMap<double>>(cfg.getParameter<edm::InputTag>("PFPVRobustPuppiWeight"))),
      usePUPPI_(cfg.getParameter<bool>("usePUPPI")),
      norm_(cfg.getParameter<double>("norm_factor")),
      ignore_leptons_(cfg.getParameter<bool>("ignore_leptons")),
      max_n_pf_(cfg.getParameter<unsigned int>("max_n_pf")),
      do_print_(cfg.getParameter<bool>("do_print")),
      session_(tensorflow::createSession(cache->graph_def)) {
  produces<pat::METCollection>();

  const int n_input = usePUPPI_ ? 7 : 6;
  const tensorflow::TensorShape shape({1, max_n_pf_, n_input});
  const tensorflow::TensorShape cat_shape({1, max_n_pf_, 1});

  input_ = tensorflow::Tensor(tensorflow::DT_FLOAT, shape);
  input_cat0_ = tensorflow::Tensor(tensorflow::DT_FLOAT, cat_shape);
  input_cat1_ = tensorflow::Tensor(tensorflow::DT_FLOAT, cat_shape);
}

void DeepMETPVRobustProducer::produce(edm::Event& event, const edm::EventSetup& setup) {
  edm::Handle<std::vector<pat::PackedCandidate>> hPFProduct;
  event.getByToken(pf_token_, hPFProduct);
  const std::vector<pat::PackedCandidate>* pfCol = hPFProduct.product();

  edm::Handle<edm::ValueMap<double>> PFPVRobustDxy;
  edm::Handle<edm::ValueMap<double>> PFPVRobustDz;
  edm::Handle<edm::ValueMap<double>> PFPVRobustPuppiWeight;

  event.getByToken(PFPVRobustDxy_, PFPVRobustDxy);
  event.getByToken(PFPVRobustDz_, PFPVRobustDz);
  if (usePUPPI_)
    event.getByToken(PFPVRobustPuppiWeight_, PFPVRobustPuppiWeight);

  const tensorflow::NamedTensorList input_list = {
      {"input", input_}, {"input_cat0", input_cat0_}, {"input_cat1", input_cat1_}};

  // Set all inputs to zero
  input_.flat<float>().setZero();
  input_cat0_.flat<float>().setZero();
  input_cat1_.flat<float>().setZero();

  size_t i_pf = 0;
  float px_leptons = 0.;
  float py_leptons = 0.;
  const float scale = 1. / norm_;
  for (unsigned ipf = 0; ipf < pfCol->size(); ipf++) {
    const auto& pf = (*pfCol)[ipf];
    edm::Ptr<pat::PackedCandidate> pf_ptr(hPFProduct, ipf);

    if (ignore_leptons_) {
      int pdg_id = std::abs(pf.pdgId());
      if (pdg_id == 11 || pdg_id == 13) {
        px_leptons += pf.px();
        py_leptons += pf.py();
        continue;
      }
    }

    // fill the tensor
    // dz, eta, mass, pt, puppi, px, py
    float* ptr = &input_.tensor<float, 3>()(0, i_pf, 0);
    //*ptr = pf.dz();
    *ptr = (*PFPVRobustDz)[pf_ptr];
    *(++ptr) = pf.eta();
    *(++ptr) = pf.mass();
    *(++ptr) = scale_and_rm_outlier(pf.pt(), scale);
    if (usePUPPI_) {
      //*(++ptr) = pf.puppiWeight();
      *(++ptr) = (*PFPVRobustPuppiWeight)[pf_ptr];
    }
    *(++ptr) = scale_and_rm_outlier(pf.px(), scale);
    *(++ptr) = scale_and_rm_outlier(pf.py(), scale);
    input_cat0_.tensor<float, 3>()(0, i_pf, 0) = charge_embedding_.at(pf.charge());
    input_cat1_.tensor<float, 3>()(0, i_pf, 0) = pdg_id_embedding_.at(pf.pdgId());

    ++i_pf;
    if (i_pf == max_n_pf_) {
      break;  // output a warning?
    }
  }

  std::vector<tensorflow::Tensor> outputs;
  const std::vector<std::string> output_names = {"output/BiasAdd"};

  // run the inference and return met
  tensorflow::run(session_, input_list, output_names, &outputs);

  // The DNN directly estimates the missing px and py, not the recoil
  float px = outputs[0].tensor<float, 2>()(0, 0) * norm_;
  float py = outputs[0].tensor<float, 2>()(0, 1) * norm_;

  px -= px_leptons;
  py -= py_leptons;

  if (do_print_) {
    std::cout << "DeepMETPVRobustProducer: MET = " << std::hypot(px, py) << " (px, py) = (" << px << ", " << py << ")"
              << std::endl;
  }

  auto pf_mets = std::make_unique<pat::METCollection>();
  const reco::Candidate::LorentzVector p4(px, py, 0., std::hypot(px, py));
  pf_mets->emplace_back(reco::MET(p4, {}));
  event.put(std::move(pf_mets));
}

std::unique_ptr<DeepMETPVRobustCache> DeepMETPVRobustProducer::initializeGlobalCache(const edm::ParameterSet& params) {
  // this method is supposed to create, initialize and return a DeepMETPVRobustCache instance
  std::unique_ptr<DeepMETPVRobustCache> cache = std::make_unique<DeepMETPVRobustCache>();

  // load the graph def and save it
  std::string graphPath = params.getParameter<std::string>("graph_path");
  if (!graphPath.empty()) {
    graphPath = edm::FileInPath(graphPath).fullPath();
    cache->graph_def = tensorflow::loadGraphDef(graphPath);
  }

  return cache;
}

void DeepMETPVRobustProducer::globalEndJob(DeepMETPVRobustCache* cache) { delete cache->graph_def; }

void DeepMETPVRobustProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("pf_src", edm::InputTag("packedPFCandidates"));
  desc.add<edm::InputTag>("PFPVRobustDxy", edm::InputTag("puppiPVRobust", "PFPVRobustDxy"));
  desc.add<edm::InputTag>("PFPVRobustDz", edm::InputTag("puppiPVRobust", "PFPVRobustDz"));
  desc.add<edm::InputTag>("PFPVRobustPuppiWeight", edm::InputTag("puppiPVRobust", "PFPVRobustPuppiWeight"));
  desc.add<bool>("usePUPPI", true);
  desc.add<bool>("ignore_leptons", false);
  desc.add<double>("norm_factor", 50.);
  desc.add<unsigned int>("max_n_pf", 4500);
  desc.add<std::string>("graph_path", "RecoMET/METPUSubtraction/data/deepmet_pvrobust/deepmet_pvrobust.pb");
  desc.add<bool>("do_print", false);
  descriptions.add("deepMETPVRobustProducer", desc);
}

DEFINE_FWK_MODULE(DeepMETPVRobustProducer);
