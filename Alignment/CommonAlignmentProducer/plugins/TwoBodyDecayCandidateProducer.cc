/// \class TwoBodyDecayCandidateProducer
///
/// EDProducer that builds a reco::VertexCompositeCandidateCollection from
/// pairs of reco::Tracks passing a configurable mass window + charge +
/// acoplanarity selection. Used as the upstream candidate-builder for the
/// dimuon ALCAREco streams (TkAlJpsiMuMu, TkAlUpsilonMuMu, TkAlZMuMu,
/// TkAlJpsiX): pair finding is moved here so the candidate object survives
/// downstream and the AlignmentTrackSelector only has to apply per-track
/// quality cuts -- mirroring the V0 ALCAREco pattern.
///
/// Each output candidate carries two reco::RecoChargedCandidate daughters
/// whose TrackRefs point into the input track collection. After the
/// AlignmentTrackSelector clones the daughter tracks into the ALCAREco
/// output, VertexCompositeCandidateRemapper rewrites these TrackRefs to
/// point at the cloned collection so downstream consumers can navigate
/// candidate -> daughter -> track without going back to generalTracks.
///
/// Unlike AlignmentTwoBodyDecayTrackSelector, all passing pairs are emitted
/// (multiplicity > 1 allowed); ordering by mother pT is preserved as the
/// natural iteration order.
///
/// Optional asymmetric-mass mode: set firstDaughterMass / secondDaughterMass
/// to assign different mass hypotheses to the positive- and negative-charge
/// daughters respectively (e.g. K+pi- for K*0 reconstruction). Set
/// tryBothChargeAssignments = True to also try the swapped assignment (K-pi+
/// hypothesis) per pair; all mass-window-passing assignments are emitted as
/// separate candidates. All existing configs that pass only daughterMass are
/// unaffected.
///
/// Optional vertex fit: set applyVertexFit = True to run a KalmanVertexFitter
/// on each candidate pair and apply minVtxProb (chi2 probability). Used for
/// K*0 and phi sub-resonance producers in TkAlJpsiX. The fitted vertex
/// position replaces the track midpoint in the output VCC.

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "TMath.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <vector>

class TwoBodyDecayCandidateProducer : public edm::stream::EDProducer<> {
public:
  explicit TwoBodyDecayCandidateProducer(const edm::ParameterSet& cfg)
      : trackToken_(consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("src"))),
        muonTag_(cfg.getParameter<edm::InputTag>("muonSrc")),
        useMuonFilter_(!muonTag_.label().empty()),
        minMass_(cfg.getParameter<double>("minMass")),
        maxMass_(cfg.getParameter<double>("maxMass")),
        daughterMass_(cfg.getParameter<double>("daughterMass")),
        // Optional asymmetric-mass parameters; default to the symmetric daughterMass value.
        firstDaughterMass_(cfg.existsAs<double>("firstDaughterMass")
                               ? cfg.getParameter<double>("firstDaughterMass")
                               : cfg.getParameter<double>("daughterMass")),
        secondDaughterMass_(cfg.existsAs<double>("secondDaughterMass")
                                ? cfg.getParameter<double>("secondDaughterMass")
                                : cfg.getParameter<double>("daughterMass")),
        daughterPdgId_(cfg.getParameter<int>("daughterPdgId")),
        firstDaughterPdgId_(cfg.existsAs<int>("firstDaughterPdgId")
                                ? cfg.getParameter<int>("firstDaughterPdgId")
                                : cfg.getParameter<int>("daughterPdgId")),
        secondDaughterPdgId_(cfg.existsAs<int>("secondDaughterPdgId")
                                 ? cfg.getParameter<int>("secondDaughterPdgId")
                                 : cfg.getParameter<int>("daughterPdgId")),
        motherPdgId_(cfg.getParameter<int>("motherPdgId")),
        applyChargeFilter_(cfg.getParameter<bool>("applyChargeFilter")),
        charge_(cfg.getParameter<int>("charge")),
        useUnsignedCharge_(cfg.getParameter<bool>("useUnsignedCharge")),
        applyAcoplanarityFilter_(cfg.getParameter<bool>("applyAcoplanarityFilter")),
        acoplanarDistance_(cfg.getParameter<double>("acoplanarDistance")),
        tryBothChargeAssignments_(cfg.existsAs<bool>("tryBothChargeAssignments")
                                      ? cfg.getParameter<bool>("tryBothChargeAssignments")
                                      : false),
        maxDaughterEta_(cfg.existsAs<double>("maxDaughterEta")
                            ? cfg.getParameter<double>("maxDaughterEta")
                            : std::numeric_limits<double>::max()),
        minDaughterPt_(cfg.existsAs<double>("minDaughterPt")
                            ? cfg.getParameter<double>("minDaughterPt")
                            : 0.0),
        applyVertexFit_(cfg.existsAs<bool>("applyVertexFit")
                            ? cfg.getParameter<bool>("applyVertexFit")
                            : false),
        minVtxProb_(cfg.existsAs<double>("minVtxProb")
                        ? cfg.getParameter<double>("minVtxProb")
                        : 0.0) {
    if (useMuonFilter_) {
      muonToken_ = consumes<edm::View<reco::Muon>>(muonTag_);
    }
    produces<reco::VertexCompositeCandidateCollection>();
  }

  void produce(edm::Event& evt, const edm::EventSetup& iSetup) override {
    edm::Handle<reco::TrackCollection> trackH;
    evt.getByToken(trackToken_, trackH);

    auto out = std::make_unique<reco::VertexCompositeCandidateCollection>();

    // If requested, restrict to tracks that are the innerTrack of a muon in
    // muonSrc. Replicates AlignmentGlobalTrackSelector::applyGlobalMuonFilter
    // upstream of the candidate combinatorics.
    std::vector<size_t> goodIdx;
    goodIdx.reserve(trackH->size());
    if (useMuonFilter_) {
      edm::Handle<edm::View<reco::Muon>> muonH;
      evt.getByToken(muonToken_, muonH);
      std::set<size_t> muTrackKeys;
      for (const auto& mu : *muonH) {
        const auto& tr = mu.innerTrack();
        if (tr.isNonnull() && tr.id() == trackH.id())
          muTrackKeys.insert(tr.key());
      }
      for (size_t i = 0; i < trackH->size(); ++i)
        if (muTrackKeys.count(i)) goodIdx.push_back(i);
    } else {
      for (size_t i = 0; i < trackH->size(); ++i) goodIdx.push_back(i);
    }

    if (goodIdx.size() < 2) {
      evt.put(std::move(out));
      return;
    }

    // Each entry: (motherPt, iPosTrack, iNegTrack, massSwapped).
    // massSwapped=false: positive track gets firstDaughterMass, negative gets secondDaughterMass.
    // massSwapped=true:  positive track gets secondDaughterMass, negative gets firstDaughterMass.
    using PairItem = std::tuple<double, size_t, size_t, bool>;
    std::vector<PairItem> ranked;
    ranked.reserve(goodIdx.size() * (goodIdx.size() - 1) / 2);

    for (size_t a = 0; a < goodIdx.size(); ++a) {
      const reco::Track& trA = (*trackH)[goodIdx[a]];
      for (size_t b = a + 1; b < goodIdx.size(); ++b) {
        const reco::Track& trB = (*trackH)[goodIdx[b]];

        if (applyChargeFilter_) {
          int sumQ = trA.charge() + trB.charge();
          if (useUnsignedCharge_) sumQ = std::abs(sumQ);
          if (sumQ != charge_) continue;
        }

        if (std::fabs(trA.eta()) > maxDaughterEta_) continue;
        if (std::fabs(trB.eta()) > maxDaughterEta_) continue;
        if (trA.pt() < minDaughterPt_) continue;
        if (trB.pt() < minDaughterPt_) continue;

        if (applyAcoplanarityFilter_) {
          if (std::fabs(reco::deltaPhi(trA.phi(), trB.phi() - M_PI)) >= acoplanarDistance_)
            continue;
        }

        // Assign pos/neg following the positive-first convention so that
        // daughter(0) = positive-charge track in all output candidates.
        size_t iPos = goodIdx[a], iNeg = goodIdx[b];
        if (trB.charge() > trA.charge()) { iPos = goodIdx[b]; iNeg = goodIdx[a]; }

        const reco::Track& trPos = (*trackH)[iPos];
        const reco::Track& trNeg = (*trackH)[iNeg];

        // Assignment 1: positive -> firstDaughterMass, negative -> secondDaughterMass.
        {
          double ePos = std::sqrt(trPos.p() * trPos.p() + firstDaughterMass_ * firstDaughterMass_);
          double eNeg = std::sqrt(trNeg.p() * trNeg.p() + secondDaughterMass_ * secondDaughterMass_);
          reco::Particle::LorentzVector p4M(
              trPos.px() + trNeg.px(), trPos.py() + trNeg.py(),
              trPos.pz() + trNeg.pz(), ePos + eNeg);
          const double mass = p4M.M();
          if (mass > minMass_ && mass < maxMass_)
            ranked.emplace_back(p4M.Pt(), iPos, iNeg, false);
        }

        // Assignment 2 (swapped): positive -> secondDaughterMass, negative -> firstDaughterMass.
        if (tryBothChargeAssignments_) {
          double ePos = std::sqrt(trPos.p() * trPos.p() + secondDaughterMass_ * secondDaughterMass_);
          double eNeg = std::sqrt(trNeg.p() * trNeg.p() + firstDaughterMass_ * firstDaughterMass_);
          reco::Particle::LorentzVector p4M(
              trPos.px() + trNeg.px(), trPos.py() + trNeg.py(),
              trPos.pz() + trNeg.pz(), ePos + eNeg);
          const double mass = p4M.M();
          if (mass > minMass_ && mass < maxMass_)
            ranked.emplace_back(p4M.Pt(), iPos, iNeg, true);
        }
      }
    }

    if (ranked.empty()) {
      evt.put(std::move(out));
      return;
    }

    std::sort(ranked.begin(), ranked.end(),
              [](const PairItem& x, const PairItem& y) { return std::get<0>(x) > std::get<0>(y); });

    // Get TransientTrackBuilder once if vertex fit is requested.
    const TransientTrackBuilder* ttb = nullptr;
    edm::ESHandle<TransientTrackBuilder> ttbH;
    if (applyVertexFit_) {
      iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttbH);
      ttb = ttbH.product();
    }

    unsigned nDroppedByVtx = 0;
    out->reserve(ranked.size());
    for (const auto& item : ranked) {
      size_t iPos = std::get<1>(item);
      size_t iNeg = std::get<2>(item);
      bool swapped = std::get<3>(item);

      const reco::Track& trPos = (*trackH)[iPos];
      const reco::Track& trNeg = (*trackH)[iNeg];

      double mPos = swapped ? secondDaughterMass_ : firstDaughterMass_;
      double mNeg = swapped ? firstDaughterMass_ : secondDaughterMass_;
      // PDG sign convention (matching existing code): positive PDG for the
      // negative-charge daughter (particle), negative PDG for the positive-charge
      // daughter (anti-particle). firstDaughterPdgId_ describes the positive-side
      // species in the standard assignment; secondDaughterPdgId_ the negative-side.
      int pdgPos = swapped ? -secondDaughterPdgId_ : -firstDaughterPdgId_;
      int pdgNeg = swapped ? +firstDaughterPdgId_ : +secondDaughterPdgId_;

      double ePos = std::sqrt(trPos.p() * trPos.p() + mPos * mPos);
      double eNeg = std::sqrt(trNeg.p() * trNeg.p() + mNeg * mNeg);

      reco::Particle::LorentzVector lvPos(trPos.px(), trPos.py(), trPos.pz(), ePos);
      reco::Particle::LorentzVector lvNeg(trNeg.px(), trNeg.py(), trNeg.pz(), eNeg);

      reco::RecoChargedCandidate dPos(trPos.charge(), lvPos,
                                      {trPos.vx(), trPos.vy(), trPos.vz()}, pdgPos);
      dPos.setTrack(reco::TrackRef(trackH, iPos));
      reco::RecoChargedCandidate dNeg(trNeg.charge(), lvNeg,
                                      {trNeg.vx(), trNeg.vy(), trNeg.vz()}, pdgNeg);
      dNeg.setTrack(reco::TrackRef(trackH, iNeg));

      reco::Particle::LorentzVector lvM = lvPos + lvNeg;
      reco::Particle::Point vtxM(0.5 * (trPos.vx() + trNeg.vx()),
                                  0.5 * (trPos.vy() + trNeg.vy()),
                                  0.5 * (trPos.vz() + trNeg.vz()));
      const int chargeM = trPos.charge() + trNeg.charge();

      if (applyVertexFit_) {
        std::vector<reco::TransientTrack> tts;
        tts.push_back(ttb->build(reco::TrackRef(trackH, iPos)));
        tts.push_back(ttb->build(reco::TrackRef(trackH, iNeg)));
        KalmanVertexFitter kvf;
        TransientVertex fv = kvf.vertex(tts);
        if (!fv.isValid()) { ++nDroppedByVtx; continue; }
        double prob = TMath::Prob(fv.totalChiSquared(), (int)fv.degreesOfFreedom());
        if (prob < minVtxProb_) { ++nDroppedByVtx; continue; }
        vtxM = reco::Particle::Point(fv.position().x(), fv.position().y(), fv.position().z());
      }

      reco::VertexCompositeCandidate cand(chargeM, lvM, vtxM, motherPdgId_);
      cand.addDaughter(dPos);
      cand.addDaughter(dNeg);
      out->push_back(std::move(cand));
    }

    LogDebug("TwoBodyDecayCandidateProducer")
        << "Built " << out->size() << " candidates from " << goodIdx.size()
        << " input tracks; " << nDroppedByVtx << " dropped by vertex fit.";
    evt.put(std::move(out));
  }

private:
  edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  edm::InputTag muonTag_;
  bool useMuonFilter_;
  edm::EDGetTokenT<edm::View<reco::Muon>> muonToken_;

  double minMass_;
  double maxMass_;
  double daughterMass_;
  double firstDaughterMass_;
  double secondDaughterMass_;
  int daughterPdgId_;
  int firstDaughterPdgId_;
  int secondDaughterPdgId_;
  int motherPdgId_;

  bool applyChargeFilter_;
  int charge_;
  bool useUnsignedCharge_;

  bool applyAcoplanarityFilter_;
  double acoplanarDistance_;

  bool tryBothChargeAssignments_;

  double maxDaughterEta_;
  double minDaughterPt_;

  bool applyVertexFit_;
  double minVtxProb_;
};

DEFINE_FWK_MODULE(TwoBodyDecayCandidateProducer);
