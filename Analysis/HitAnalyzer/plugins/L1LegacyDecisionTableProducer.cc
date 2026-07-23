// L1LegacyDecisionTableProducer
//
// The 2016 AlCaReco carries the LEGACY L1 record (L1GlobalTriggerReadoutRecord,
// "gtDigis"), not the Run3 gtStage2Digis that stock NanoAOD L1 tooling reads.
// This emits a singleton NanoAOD FlatTable "L1" with the finalOR and the full
// 128-bit algo + 64-bit technical decision words packed into uint32 columns,
// so any algo/tech bit is recoverable offline once the L1 menu is known:
//   fired(bit N) = (L1_algo{N/32} >> (N%32)) & 1
//
// Packing (rather than 192 individual bool branches) keeps the event size tiny
// and lossless. The L1 menu itself is intentionally not applied here.

#include <memory>
#include <vector>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

namespace ana_hitanalyzer {

namespace {
// Pack a vector<bool> into ceil(n/32) little-endian uint32 words.
std::vector<uint32_t> packBits(const std::vector<bool>& bits, unsigned nwords) {
  std::vector<uint32_t> words(nwords, 0u);
  for (size_t i = 0; i < bits.size(); ++i)
    if (bits[i]) words[i / 32] |= (1u << (i % 32));
  return words;
}
}  // namespace

class L1LegacyDecisionTableProducer : public edm::stream::EDProducer<> {
public:
  explicit L1LegacyDecisionTableProducer(const edm::ParameterSet& cfg)
      : token_(consumes<L1GlobalTriggerReadoutRecord>(cfg.getParameter<edm::InputTag>("src"))),
        name_(cfg.getParameter<std::string>("name")) {
    produces<nanoaod::FlatTable>();
  }

  void produce(edm::Event& iEvent, const edm::EventSetup&) override {
    edm::Handle<L1GlobalTriggerReadoutRecord> h;
    iEvent.getByToken(token_, h);

    auto tab = std::make_unique<nanoaod::FlatTable>(1, name_, /*singleton=*/true);
    if (h.isValid()) {
      const std::vector<uint32_t> algo = packBits(h->decisionWord(), 4);          // 128 bits
      const std::vector<uint32_t> tech = packBits(h->technicalTriggerWord(), 2);  // 64 bits
      tab->template addColumnValue<bool>("finalOR", h->finalOR() != 0,
                                         "legacy L1 global final-OR (event fired L1)");
      for (unsigned w = 0; w < algo.size(); ++w)
        tab->template addColumnValue<uint32_t>(
            "algo" + std::to_string(w), algo[w],
            "algo decision bits " + std::to_string(32 * w) + ".." + std::to_string(32 * w + 31));
      for (unsigned w = 0; w < tech.size(); ++w)
        tab->template addColumnValue<uint32_t>(
            "tech" + std::to_string(w), tech[w],
            "technical trigger bits " + std::to_string(32 * w) + ".." + std::to_string(32 * w + 31));
    } else {
      tab->template addColumnValue<bool>("finalOR", false, "record absent");
    }
    iEvent.put(std::move(tab));
  }

private:
  const edm::EDGetTokenT<L1GlobalTriggerReadoutRecord> token_;
  const std::string name_;
};

}  // namespace ana_hitanalyzer

using ana_hitanalyzer::L1LegacyDecisionTableProducer;
DEFINE_FWK_MODULE(L1LegacyDecisionTableProducer);
