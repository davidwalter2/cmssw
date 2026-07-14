#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/Event.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "TLorentzVector.h"

#include "Alignment/CommonAlignmentProducer/interface/AlignmentThreeBodyDecayTrackSelector.h"

using namespace edm;
using namespace std;

namespace {
  struct ThreeBodyCandidate {
    double pt;
    const reco::Track* kaon;
    const reco::Track* pion;
    const reco::Track* softPion;
  };

  void addUniqueTrack(AlignmentThreeBodyDecayTrackSelector::Tracks& tracks, const reco::Track* track) {
    if (std::find(tracks.begin(), tracks.end(), track) == tracks.end()) {
      tracks.push_back(track);
    }
  }
}

AlignmentThreeBodyDecayTrackSelector::AlignmentThreeBodyDecayTrackSelector(const edm::ParameterSet& cfg,
                                                                           edm::ConsumesCollector& iC) {
  (void)iC;
  LogDebug("Alignment") << "> applying three body decay Trackfilter ...";

  theMassrangeSwitch = cfg.getParameter<bool>("applyMassrangeFilter");
  theIntermediateMassSwitch = cfg.getParameter<bool>("applyIntermediateMassrangeFilter");
  theMassDifferenceSwitch = cfg.getParameter<bool>("applyMassDifferenceFilter");
  theChargeSwitch = cfg.getParameter<bool>("applyChargeFilter");

  theMinMass = cfg.getParameter<double>("minXMass");
  theMaxMass = cfg.getParameter<double>("maxXMass");
  theMinIntermediateMass = cfg.getParameter<double>("minIntermediateMass");
  theMaxIntermediateMass = cfg.getParameter<double>("maxIntermediateMass");
  theMinMassDifference = cfg.getParameter<double>("minMassDifference");
  theMaxMassDifference = cfg.getParameter<double>("maxMassDifference");
  theFirstDaughterMass = cfg.getParameter<double>("firstDaughterMass");
  theSecondDaughterMass = cfg.getParameter<double>("secondDaughterMass");
  theThirdDaughterMass = cfg.getParameter<double>("thirdDaughterMass");
  theFirstDaughterPtMin =
      cfg.existsAs<double>("firstDaughterPtMin") ? cfg.getParameter<double>("firstDaughterPtMin") : -1.;
  theSecondDaughterPtMin =
      cfg.existsAs<double>("secondDaughterPtMin") ? cfg.getParameter<double>("secondDaughterPtMin") : -1.;
  theThirdDaughterPtMin =
      cfg.existsAs<double>("thirdDaughterPtMin") ? cfg.getParameter<double>("thirdDaughterPtMin") : -1.;
  theCandNumber = cfg.getParameter<unsigned int>("numberOfCandidates");

  theCharge = cfg.getParameter<int>("charge");
  theUnsignedSwitch = cfg.getParameter<bool>("useUnsignedCharge");
  eventsChecked_ = 0;
  eventsWithCandidates_ = 0;
  totalPassingCandidates_ = 0;
  totalSelectedFlatTracks_ = 0;
}

AlignmentThreeBodyDecayTrackSelector::~AlignmentThreeBodyDecayTrackSelector() {
  edm::LogInfo("Alignment") << "AlignmentThreeBodyDecayTrackSelector D* summary"
                            << " eventsChecked=" << eventsChecked_
                            << " eventsWithCandidates=" << eventsWithCandidates_
                            << " totalPassingCandidates=" << totalPassingCandidates_
                            << " totalSelectedFlatTracks=" << totalSelectedFlatTracks_;
}

bool AlignmentThreeBodyDecayTrackSelector::useThisFilter() {
  return theMassrangeSwitch || theIntermediateMassSwitch || theMassDifferenceSwitch || theChargeSwitch;
}

AlignmentThreeBodyDecayTrackSelector::Tracks AlignmentThreeBodyDecayTrackSelector::select(const Tracks& tracks,
                                                                                          const edm::Event& iEvent,
                                                                                          const edm::EventSetup& iSetup) {
  (void)iEvent;
  (void)iSetup;
  Tracks result = tracks;

  if (useThisFilter())
    result = checkMass(result);

  LogDebug("Alignment") << ">  ThreeBodyDecay tracks all,kept: " << tracks.size() << "," << result.size();
  return result;
}

AlignmentThreeBodyDecayTrackSelector::Tracks AlignmentThreeBodyDecayTrackSelector::checkMass(const Tracks& cands) const {
  Tracks result;
  ++eventsChecked_;

  if (cands.size() < 3)
    return result;

  vector<ThreeBodyCandidate> candCollection;

  for (unsigned int iK = 0; iK < cands.size(); ++iK) {
    const reco::Track* kaon = cands[iK];
    if (theFirstDaughterPtMin >= 0. && kaon->pt() < theFirstDaughterPtMin)
      continue;
    const int kaonCharge = kaon->charge();

    TLorentzVector kaonP4;
    kaonP4.SetXYZT(kaon->px(),
                   kaon->py(),
                   kaon->pz(),
                   std::sqrt(kaon->p() * kaon->p() + theFirstDaughterMass * theFirstDaughterMass));

    for (unsigned int iPi = 0; iPi < cands.size(); ++iPi) {
      if (iPi == iK)
        continue;

      const reco::Track* pion = cands[iPi];
      if (theSecondDaughterPtMin >= 0. && pion->pt() < theSecondDaughterPtMin)
        continue;
      if (pion->charge() != -kaonCharge)
        continue;

      TLorentzVector pionP4;
      pionP4.SetXYZT(pion->px(),
                     pion->py(),
                     pion->pz(),
                     std::sqrt(pion->p() * pion->p() + theSecondDaughterMass * theSecondDaughterMass));

      const TLorentzVector intermediate = kaonP4 + pionP4;
      if (theIntermediateMassSwitch &&
          !(intermediate.M() > theMinIntermediateMass && intermediate.M() < theMaxIntermediateMass))
        continue;

      for (unsigned int iPis = 0; iPis < cands.size(); ++iPis) {
        if (iPis == iK || iPis == iPi)
          continue;

        const reco::Track* softPion = cands[iPis];
        if (theThirdDaughterPtMin >= 0. && softPion->pt() < theThirdDaughterPtMin)
          continue;
        if (softPion->charge() != pion->charge())
          continue;

        if (theChargeSwitch && !checkCharge(kaon, pion, softPion))
          continue;

        TLorentzVector softPionP4;
        softPionP4.SetXYZT(softPion->px(),
                           softPion->py(),
                           softPion->pz(),
                           std::sqrt(softPion->p() * softPion->p() + theThirdDaughterMass * theThirdDaughterMass));

        const TLorentzVector mother = intermediate + softPionP4;
        const double massDifference = mother.M() - intermediate.M();

        if (theMassrangeSwitch && !(mother.M() > theMinMass && mother.M() < theMaxMass))
          continue;

        if (theMassDifferenceSwitch &&
            !(massDifference > theMinMassDifference && massDifference < theMaxMassDifference))
          continue;

        candCollection.push_back({mother.Pt(), kaon, pion, softPion});
      }
    }
  }

  if (candCollection.empty())
    return result;
  ++eventsWithCandidates_;
  totalPassingCandidates_ += candCollection.size();

  sort(candCollection.begin(), candCollection.end(), [](const ThreeBodyCandidate& a, const ThreeBodyCandidate& b) {
    return a.pt > b.pt;
  });

  const unsigned int maxCandidates = theCandNumber == 0 ? candCollection.size() : theCandNumber;
  for (unsigned int i = 0; i < candCollection.size() && i < maxCandidates; ++i) {
    const ThreeBodyCandidate& candidate = candCollection[i];
    addUniqueTrack(result, candidate.kaon);
    addUniqueTrack(result, candidate.pion);
    addUniqueTrack(result, candidate.softPion);
  }
  totalSelectedFlatTracks_ += result.size();

  return result;
}

bool AlignmentThreeBodyDecayTrackSelector::checkCharge(const reco::Track* kaon,
                                                       const reco::Track* pion,
                                                       const reco::Track* softPion) const {
  const bool validD0Charges = (kaon->charge() + pion->charge() == 0);
  const bool validSoftCharge = (pion->charge() == softPion->charge());
  if (!validD0Charges || !validSoftCharge)
    return false;

  int totalCharge = kaon->charge() + pion->charge() + softPion->charge();
  if (theUnsignedSwitch)
    totalCharge = std::abs(totalCharge);

  return totalCharge == theCharge;
}
