// JpsiXCandidatePreselectorForCorrections
//
// Candidate preselection for the CVH *grads* path -- the gradient/Hessian
// production that feeds the global-correction fit, not the physics refit. The
// N-body maker exposes no candidate-level cuts beyond minValidHitsPerLeg, and
// every quantity the selection needs already exists in the same job as a
// JpsiXKinematicFitProducer ValueMap, so the filtering belongs here: read those
// maps, emit a reduced VertexCompositeCandidateCollection, and let the maker run
// over the survivors.
//
// TWO RULES that are easy to get wrong and expensive if you do.
//
// 1. EVERY CUT READS THE RAW ARM. The `ref*` ValueMaps are computed from
//    CVH-refit tracks, i.e. from the very corrections the calibration is
//    fitting, so selecting on them would make the selection a function of the
//    answer. The driver aliases the unsuffixed physics columns (cvhFit*,
//    dimuon*) to the REFIT arm, so "the obvious tag" is the wrong one here.
//    This module therefore takes each map as an explicit InputTag and the
//    driver wires the `raw*` instances.
//
// 2. AT MOST ONE CANDIDATE PER EVENT is a correctness condition, not an
//    optimisation. At ~1.7 candidates/event the same dimuon's hits would
//    otherwise enter the summed global Hessian several times, and the sum would
//    no longer be the information of an independent sample -- unlike the J/psi
//    and K0s reference channels, whose AlCaRecos give effectively unique track
//    pairs. When several candidates survive, the one with the largest mother
//    vertex probability wins.
//
// Selection philosophy: loosen the KINEMATICS, keep the DISPLACEMENT.
// Displacement and vertex-quality cuts suppress combinatorial background without
// touching bachelor kinematics; the analysis-level kinematic cuts sculpt exactly
// the soft/forward bachelor phase space the channel exists to probe (the
// bachelor spectrum has median pT ~0.53 GeV with 82% below 1 GeV, so the Bmm5
// analysis cut pT>1 keeps 18% of them). Defaults here are the loose nominal;
// the driver's selTight=True switch supplies the analysis values as a labelled
// systematic variation.
//
// Leaf order follows the stage-1 layout, as in CandidateLeafTrackIndexProducer:
// daughter(0) is the J/psi composite, so leaves 0 and 1 are the muons and the
// remaining leaves are bachelors.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

namespace ana_hitanalyzer {

namespace {
void collectLeaves(const reco::Candidate* c,
                   std::vector<const reco::RecoChargedCandidate*>& out) {
  if (c == nullptr) return;
  const auto* rcc = dynamic_cast<const reco::RecoChargedCandidate*>(c);
  if (rcc != nullptr) {
    out.push_back(rcc);
    return;
  }
  for (size_t i = 0; i < c->numberOfDaughters(); ++i)
    collectLeaves(c->daughter(i), out);
}

// Cut names in the order they are applied, for the end-of-job cutflow. Keep in
// sync with the `fail(k)` indices in produce().
const char* const kCutNames[] = {"input",     "fitOk",      "massWindow",
                                 "bVtxProb",  "mmVtxProb",  "mmAlphaBS",
                                 "mmSl3d",    "muonKin",    "bachelorKin",
                                 "dimuonPt",  "bestPerEvent"};
constexpr unsigned kNCuts = sizeof(kCutNames) / sizeof(kCutNames[0]);
}  // namespace

class JpsiXCandidatePreselectorForCorrections : public edm::stream::EDProducer<> {
public:
  explicit JpsiXCandidatePreselectorForCorrections(const edm::ParameterSet& cfg)
      : candToken_(consumes<reco::VertexCompositeCandidateCollection>(
            cfg.getParameter<edm::InputTag>("src"))),
        fitOkToken_(consumes<edm::ValueMap<int>>(
            cfg.getParameter<edm::InputTag>("fitOk"))),
        massToken_(consumes<edm::ValueMap<float>>(
            cfg.getParameter<edm::InputTag>("fitMass"))),
        bVtxProbToken_(consumes<edm::ValueMap<float>>(
            cfg.getParameter<edm::InputTag>("fitVtxProb"))),
        mmVtxProbToken_(consumes<edm::ValueMap<float>>(
            cfg.getParameter<edm::InputTag>("dimuonVtxProb"))),
        mmAlphaToken_(consumes<edm::ValueMap<float>>(
            cfg.getParameter<edm::InputTag>("dimuonAlphaBS"))),
        mmSl3dToken_(consumes<edm::ValueMap<float>>(
            cfg.getParameter<edm::InputTag>("dimuonSl3d"))),
        massCentre_(cfg.getParameter<double>("massCentre")),
        massHalfWindow_(cfg.getParameter<double>("massHalfWindow")),
        bVtxProbMin_(cfg.getParameter<double>("fitVtxProbMin")),
        mmVtxProbMin_(cfg.getParameter<double>("dimuonVtxProbMin")),
        mmAlphaMax_(cfg.getParameter<double>("dimuonAlphaBSMax")),
        mmSl3dMin_(cfg.getParameter<double>("dimuonSl3dMin")),
        muonPtMin_(cfg.getParameter<double>("muonPtMin")),
        muonEtaMax_(cfg.getParameter<double>("muonEtaMax")),
        bachPtMin_(cfg.getParameter<double>("bachelorPtMin")),
        bachPtMax_(cfg.getParameter<double>("bachelorPtMax")),
        bachEtaMax_(cfg.getParameter<double>("bachelorEtaMax")),
        dimuonPtMin_(cfg.getParameter<double>("dimuonPtMin")),
        bestPerEvent_(cfg.getParameter<bool>("bestPerEvent")),
        counts_(kNCuts, 0ULL) {
    produces<reco::VertexCompositeCandidateCollection>();
    // Provenance: index of each survivor in the input collection, so a
    // downstream consumer (or a debugging session) can go back to the parent.
    produces<std::vector<int>>("srcIdx");
  }

  void produce(edm::Event& iEvent, const edm::EventSetup&) override {
    edm::Handle<reco::VertexCompositeCandidateCollection> candH;
    iEvent.getByToken(candToken_, candH);
    edm::Handle<edm::ValueMap<int>> fitOkH;
    iEvent.getByToken(fitOkToken_, fitOkH);
    edm::Handle<edm::ValueMap<float>> massH, bVtxProbH, mmVtxProbH, mmAlphaH, mmSl3dH;
    iEvent.getByToken(massToken_, massH);
    iEvent.getByToken(bVtxProbToken_, bVtxProbH);
    iEvent.getByToken(mmVtxProbToken_, mmVtxProbH);
    iEvent.getByToken(mmAlphaToken_, mmAlphaH);
    iEvent.getByToken(mmSl3dToken_, mmSl3dH);

    auto out = std::make_unique<reco::VertexCompositeCandidateCollection>();
    auto outIdx = std::make_unique<std::vector<int>>();

    const size_t n = candH->size();
    counts_[0] += n;

    // Surviving candidates and their ranking quantity (mother vertex prob).
    std::vector<std::pair<float, size_t>> kept;
    kept.reserve(n);

    for (size_t ic = 0; ic < n; ++ic) {
      const edm::Ref<reco::VertexCompositeCandidateCollection> ref(candH, ic);
      unsigned stage = 1;
      auto pass = [&](bool ok) {
        if (ok) { ++counts_[stage]; ++stage; }
        return ok;
      };

      if (!pass((*fitOkH)[ref] != 0)) continue;
      const float m = (*massH)[ref];
      if (!pass(std::abs(m - massCentre_) < massHalfWindow_)) continue;
      const float bvp = (*bVtxProbH)[ref];
      if (!pass(bvp > bVtxProbMin_)) continue;
      if (!pass((*mmVtxProbH)[ref] > mmVtxProbMin_)) continue;
      if (!pass((*mmAlphaH)[ref] < mmAlphaMax_)) continue;
      if (!pass((*mmSl3dH)[ref] > mmSl3dMin_)) continue;

      // Kinematics from the stage-1 leaves (raw track-based four-vectors; the
      // candidate's own p4 is a raw sum by design).
      std::vector<const reco::RecoChargedCandidate*> leaves;
      collectLeaves(&(*candH)[ic], leaves);

      // A candidate with fewer than three leaves is not a J/psi + bachelor and
      // fails here rather than dropping out uncounted -- so the cutflow's
      // muonKin row accounts for every candidate that reaches it.
      bool muonsOk = leaves.size() >= 3;
      for (size_t il = 0; muonsOk && il < 2; ++il)
        if (leaves[il]->pt() < muonPtMin_ ||
            std::abs(leaves[il]->eta()) > muonEtaMax_) muonsOk = false;
      if (!pass(muonsOk)) continue;

      bool bachOk = true;
      for (size_t il = 2; il < leaves.size(); ++il)
        if (leaves[il]->pt() < bachPtMin_ || leaves[il]->pt() > bachPtMax_ ||
            std::abs(leaves[il]->eta()) > bachEtaMax_) bachOk = false;
      if (!pass(bachOk)) continue;

      const double mmPt = (leaves[0]->p4() + leaves[1]->p4()).pt();
      if (!pass(mmPt > dimuonPtMin_)) continue;

      kept.emplace_back(bvp, ic);
    }

    if (bestPerEvent_ && kept.size() > 1) {
      // Largest mother vertex probability wins. max_element on the pair
      // compares the probability first, which is what we want.
      const auto best = std::max_element(kept.begin(), kept.end());
      const auto chosen = *best;
      kept.assign(1, chosen);
    }
    counts_[kNCuts - 1] += kept.size();

    for (const auto& kv : kept) {
      out->push_back((*candH)[kv.second]);
      outIdx->push_back(static_cast<int>(kv.second));
    }

    iEvent.put(std::move(out));
    iEvent.put(std::move(outIdx), "srcIdx");
  }

  void endStream() override {
    const double denom = counts_[0] > 0 ? static_cast<double>(counts_[0]) : 1.;
    std::string msg = "JpsiXCandidatePreselectorForCorrections cutflow (candidates):";
    for (unsigned k = 0; k < kNCuts; ++k) {
      char buf[128];
      snprintf(buf, sizeof(buf), "\n  %-14s %10llu  cum.eff %7.4f", kCutNames[k],
               static_cast<unsigned long long>(counts_[k]),
               static_cast<double>(counts_[k]) / denom);
      msg += buf;
    }
    edm::LogPrint("JpsiXCandidatePreselectorForCorrections") << msg;
  }

private:
  const edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> candToken_;
  const edm::EDGetTokenT<edm::ValueMap<int>> fitOkToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> massToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> bVtxProbToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> mmVtxProbToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> mmAlphaToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> mmSl3dToken_;

  const double massCentre_;
  const double massHalfWindow_;
  const double bVtxProbMin_;
  const double mmVtxProbMin_;
  const double mmAlphaMax_;
  const double mmSl3dMin_;
  const double muonPtMin_;
  const double muonEtaMax_;
  const double bachPtMin_;
  const double bachPtMax_;
  const double bachEtaMax_;
  const double dimuonPtMin_;
  const bool bestPerEvent_;

  std::vector<unsigned long long> counts_;
};

}  // namespace ana_hitanalyzer

using ana_hitanalyzer::JpsiXCandidatePreselectorForCorrections;
DEFINE_FWK_MODULE(JpsiXCandidatePreselectorForCorrections);
