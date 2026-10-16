// RecoVertexTableProducer
//
// The PV and OtherPV tables of the stock NanoAOD VertexTableProducer for an
// input without PF candidates and secondary vertices (the AlCaReco): the same
// columns, names, documentation and precisions from the reco::Vertex
// collection and its score ValueMap, except PV_sumpt2, PV_sumpx and PV_sumpy,
// which are sums over the PF candidates.

#include <memory>
#include <string>
#include <vector>

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

class RecoVertexTableProducer : public edm::global::EDProducer<> {
public:
  explicit RecoVertexTableProducer(const edm::ParameterSet& cfg)
      : pvs_(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("pvSrc"))),
        pvsScore_(consumes<edm::ValueMap<float>>(cfg.getParameter<edm::InputTag>("pvSrc"))),
        goodPvCut_(cfg.getParameter<std::string>("goodPvCut"), true),
        goodPvCutString_(cfg.getParameter<std::string>("goodPvCut")),
        pvName_(cfg.getParameter<std::string>("pvName")) {
    produces<nanoaod::FlatTable>("pv");
    produces<nanoaod::FlatTable>("otherPVs");
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("pvSrc", edm::InputTag("offlinePrimaryVertices"))
        ->setComment("primary vertices and, under the same tag, their score ValueMap<float>");
    desc.add<std::string>("goodPvCut", "!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2")
        ->setComment("selection of the good vertices counted in npvsGood");
    desc.add<std::string>("pvName", "PV");
    descriptions.addWithDefaultLabel(desc);
  }

  void produce(edm::StreamID, edm::Event& evt, const edm::EventSetup&) const override {
    const auto pvsIn = evt.getHandle(pvs_);
    const auto& pvsScoreProd = evt.get(pvsScore_);
    const auto& pvs = *pvsIn;

    auto pvTable = std::make_unique<nanoaod::FlatTable>(1, pvName_, true);
    pvTable->addColumnValue<float>("ndof", pvs[0].ndof(), "main primary vertex number of degree of freedom", 8);
    pvTable->addColumnValue<float>("x", pvs[0].position().x(), "main primary vertex position x coordinate", 10);
    pvTable->addColumnValue<float>("y", pvs[0].position().y(), "main primary vertex position y coordinate", 10);
    pvTable->addColumnValue<float>("z", pvs[0].position().z(), "main primary vertex position z coordinate", 16);
    pvTable->addColumnValue<float>("chi2", pvs[0].normalizedChi2(), "main primary vertex reduced chi2", 8);
    int goodPVs = 0;
    for (const auto& pv : pvs)
      if (goodPvCut_(pv))
        goodPVs++;
    pvTable->addColumnValue<uint8_t>("npvs", pvs.size(), "total number of reconstructed primary vertices");
    pvTable->addColumnValue<uint8_t>(
        "npvsGood", goodPVs, "number of good reconstructed primary vertices. selection:" + goodPvCutString_);
    pvTable->addColumnValue<float>(
        "score", pvsScoreProd.get(pvsIn.id(), 0), "main primary vertex score, i.e. sum pt2 of clustered objects", 8);

    auto otherPVsTable =
        std::make_unique<nanoaod::FlatTable>(pvs.size() > 4 ? 3 : pvs.size() - 1, "Other" + pvName_, false);
    std::vector<float> pvsz;
    std::vector<float> pvscores;
    for (std::size_t i = 1; i < pvs.size() && i < 4; i++) {
      pvsz.push_back(pvs[i].position().z());
      pvscores.push_back(pvsScoreProd.get(pvsIn.id(), i));
    }
    otherPVsTable->addColumn<float>("z", pvsz, "Z position of other primary vertices, excluding the main PV", 8);
    otherPVsTable->addColumn<float>("score", pvscores, "scores of other primary vertices, excluding the main PV", 8);

    evt.put(std::move(pvTable), "pv");
    evt.put(std::move(otherPVsTable), "otherPVs");
  }

private:
  const edm::EDGetTokenT<reco::VertexCollection> pvs_;
  const edm::EDGetTokenT<edm::ValueMap<float>> pvsScore_;
  const StringCutObjectSelector<reco::Vertex> goodPvCut_;
  const std::string goodPvCutString_;
  const std::string pvName_;
};

DEFINE_FWK_MODULE(RecoVertexTableProducer);
