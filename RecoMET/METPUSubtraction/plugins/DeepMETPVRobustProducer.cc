// WMass: DeepMET evaluated on the PV-robust inputs of PuppiPVRobustProducer --
// the packed candidates' dz and PUPPI weight with respect to the vertex closest
// to the leading muon instead of the leading PV -- with the dedicated models of
// WmassNanoProd_10_6_26 (PR #33): 7 inputs (dz, eta, mass, pt, puppi weight, px,
// py) or 6 without PUPPI. Structure of the 15_0 DeepMETProducer (SessionCache).
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "RecoMET/METPUSubtraction/interface/DeepMETHelp.h"

using namespace deepmet_helper;

class DeepMETPVRobustProducer : public edm::stream::EDProducer<edm::GlobalCache<tensorflow::SessionCache>> {
public:
  explicit DeepMETPVRobustProducer(const edm::ParameterSet&, const tensorflow::SessionCache*);
  void produce(edm::Event& event, const edm::EventSetup& setup) override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  static std::unique_ptr<tensorflow::SessionCache> initializeGlobalCache(const edm::ParameterSet&);
  static void globalEndJob(tensorflow::SessionCache*) {}

private:
  const edm::EDGetTokenT<std::vector<pat::PackedCandidate>> pf_token_;
  const edm::EDGetTokenT<edm::ValueMap<double>> dxy_token_;
  const edm::EDGetTokenT<edm::ValueMap<double>> dz_token_;
  const edm::EDGetTokenT<edm::ValueMap<double>> puppi_token_;
  const bool usePUPPI_;
  const float norm_;
  const bool ignore_leptons_;
  const unsigned int max_n_pf_;

  const tensorflow::Session* session_;

  tensorflow::Tensor input_;
  tensorflow::Tensor input_cat0_;
  tensorflow::Tensor input_cat1_;
};

DeepMETPVRobustProducer::DeepMETPVRobustProducer(const edm::ParameterSet& cfg, const tensorflow::SessionCache* cache)
    : pf_token_(consumes<std::vector<pat::PackedCandidate>>(cfg.getParameter<edm::InputTag>("pf_src"))),
      dxy_token_(consumes<edm::ValueMap<double>>(cfg.getParameter<edm::InputTag>("PFPVRobustDxy"))),
      dz_token_(consumes<edm::ValueMap<double>>(cfg.getParameter<edm::InputTag>("PFPVRobustDz"))),
      puppi_token_(consumes<edm::ValueMap<double>>(cfg.getParameter<edm::InputTag>("PFPVRobustPuppiWeight"))),
      usePUPPI_(cfg.getParameter<bool>("usePUPPI")),
      norm_(cfg.getParameter<double>("norm_factor")),
      ignore_leptons_(cfg.getParameter<bool>("ignore_leptons")),
      max_n_pf_(cfg.getParameter<unsigned int>("max_n_pf")),
      session_(cache->getSession()) {
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
  const edm::ValueMap<double>& dz = event.get(dz_token_);
  const edm::ValueMap<double>* puppi = usePUPPI_ ? &event.get(puppi_token_) : nullptr;

  const tensorflow::NamedTensorList input_list = {
      {"input", input_}, {"input_cat0", input_cat0_}, {"input_cat1", input_cat1_}};

  input_.flat<float>().setZero();
  input_cat0_.flat<float>().setZero();
  input_cat1_.flat<float>().setZero();

  size_t i_pf = 0;
  float px_leptons = 0.;
  float py_leptons = 0.;
  const float scale = 1. / norm_;
  for (unsigned ipf = 0; ipf < hPFProduct->size(); ipf++) {
    const auto& pf = (*hPFProduct)[ipf];
    const edm::Ptr<pat::PackedCandidate> pf_ptr(hPFProduct, ipf);

    if (ignore_leptons_) {
      int pdg_id = std::abs(pf.pdgId());
      if (pdg_id == 11 || pdg_id == 13) {
        px_leptons += pf.px();
        py_leptons += pf.py();
        continue;
      }
    }

    // dz, eta, mass, pt, [puppi], px, py -- the training layout of the PV-robust models
    float* ptr = &input_.tensor<float, 3>()(0, i_pf, 0);
    *ptr = dz[pf_ptr];
    *(++ptr) = pf.eta();
    *(++ptr) = pf.mass();
    *(++ptr) = scale_and_rm_outlier(pf.pt(), scale);
    if (usePUPPI_)
      *(++ptr) = (*puppi)[pf_ptr];
    *(++ptr) = scale_and_rm_outlier(pf.px(), scale);
    *(++ptr) = scale_and_rm_outlier(pf.py(), scale);
    input_cat0_.tensor<float, 3>()(0, i_pf, 0) = charge_embedding.at(pf.charge());
    input_cat1_.tensor<float, 3>()(0, i_pf, 0) = pdg_id_embedding.at(pf.pdgId());

    ++i_pf;
    if (i_pf == max_n_pf_) {
      break;
    }
  }

  std::vector<tensorflow::Tensor> outputs;
  const std::vector<std::string> output_names = {"output/BiasAdd"};
  tensorflow::run(session_, input_list, output_names, &outputs);

  // the DNN estimates the missing px and py directly
  float px = outputs[0].tensor<float, 2>()(0, 0) * norm_;
  float py = outputs[0].tensor<float, 2>()(0, 1) * norm_;
  px -= px_leptons;
  py -= py_leptons;

  auto pf_mets = std::make_unique<pat::METCollection>();
  const reco::Candidate::LorentzVector p4(px, py, 0., std::hypot(px, py));
  pf_mets->emplace_back(reco::MET(p4, {}));
  event.put(std::move(pf_mets));
}

std::unique_ptr<tensorflow::SessionCache> DeepMETPVRobustProducer::initializeGlobalCache(
    const edm::ParameterSet& params) {
  std::string graphPath = edm::FileInPath(params.getParameter<std::string>("graph_path")).fullPath();
  return std::make_unique<tensorflow::SessionCache>(graphPath);
}

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
  descriptions.add("deepMETPVRobustProducer", desc);
}

DEFINE_FWK_MODULE(DeepMETPVRobustProducer);
