#ifndef Alignment_CommonAlignmentAlgorithm_AlignmentThreeBodyDecayTrackSelector_h
#define Alignment_CommonAlignmentAlgorithm_AlignmentThreeBodyDecayTrackSelector_h

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <vector>

#include <DataFormats/TrackReco/interface/TrackFwd.h>

namespace edm {
  class Event;
  class EventSetup;
}  // namespace edm

class AlignmentThreeBodyDecayTrackSelector {
public:
  typedef std::vector<const reco::Track*> Tracks;

  AlignmentThreeBodyDecayTrackSelector(const edm::ParameterSet& cfg, edm::ConsumesCollector& iC);
  ~AlignmentThreeBodyDecayTrackSelector();

  Tracks select(const Tracks& tracks, const edm::Event& iEvent, const edm::EventSetup& iSetup);
  bool useThisFilter();

private:
  Tracks checkMass(const Tracks& cands) const;
  bool checkCharge(const reco::Track* kaon, const reco::Track* pion, const reco::Track* softPion) const;

  bool theMassrangeSwitch;
  bool theIntermediateMassSwitch;
  bool theMassDifferenceSwitch;
  bool theChargeSwitch;

  double theMinMass;
  double theMaxMass;
  double theMinIntermediateMass;
  double theMaxIntermediateMass;
  double theMinMassDifference;
  double theMaxMassDifference;
  double theFirstDaughterMass;
  double theSecondDaughterMass;
  double theThirdDaughterMass;
  double theFirstDaughterPtMin;
  double theSecondDaughterPtMin;
  double theThirdDaughterPtMin;

  unsigned int theCandNumber;

  int theCharge;
  bool theUnsignedSwitch;

  mutable unsigned long long eventsChecked_;
  mutable unsigned long long eventsWithCandidates_;
  mutable unsigned long long totalPassingCandidates_;
  mutable unsigned long long totalSelectedFlatTracks_;
};

#endif
