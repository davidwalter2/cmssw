#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"

#include "Analysis/HitAnalyzer/interface/ParticleProperties.h"
#include "TMath.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

namespace {
struct FitQuality {
  bool isValid = false;
  double chi2 = -1.;
  double dof = -1.;
  double pval = -1.;
  bool isGood = false;
};

struct CandidateSummary {
  float dstPt;
  reco::Track kaon;
  reco::Track pion;
};
}  // namespace

class DstToD0PiCandidateProducer : public edm::stream::EDProducer<> {
public:
  explicit DstToD0PiCandidateProducer(const edm::ParameterSet&);
  ~DstToD0PiCandidateProducer() override = default;

private:
  void produce(edm::Event&, const edm::EventSetup&) override;
  bool checkCharge(const reco::Track* kaon, const reco::Track* pion, const reco::Track* softPion) const;
  FitQuality fitQuality(RefCountedKinematicTree tree) const;
  RefCountedKinematicTree fitD0(const edm::EventSetup&, const reco::Track&, const reco::Track&) const;
  RefCountedKinematicTree fitDst(const edm::EventSetup&, const reco::Track&, const RefCountedKinematicParticle&) const;

  edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transTrackBuilderToken_;

  double kaonMass_;
  double pionMass_;
  double softPionMass_;
  double kaonMassErr_;
  double pionMassErr_;
  double minD0Mass_;
  double maxD0Mass_;
  double minDstMass_;
  double maxDstMass_;
  double minDeltaMass_;
  double maxDeltaMass_;
  double pvalMin_;
  bool applyChargeFilter_;
  bool useUnsignedCharge_;
  int targetCharge_;
};

DstToD0PiCandidateProducer::DstToD0PiCandidateProducer(const edm::ParameterSet& iConfig)
    : trackToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      transTrackBuilderToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
      // Kaon / pion masses + mass uncertainties looked up from the PDG
      // table (Analysis/HitAnalyzer/interface/ParticleProperties.h).
      // Single source of truth for daughter-particle properties.
      kaonMass_   (ana_hitanalyzer::getParticleProperties("kaon").mass),
      pionMass_   (ana_hitanalyzer::getParticleProperties("pi").mass),
      softPionMass_(ana_hitanalyzer::getParticleProperties("pi").mass),
      kaonMassErr_(ana_hitanalyzer::getParticleProperties("kaon").massErr),
      pionMassErr_(ana_hitanalyzer::getParticleProperties("pi").massErr),
      minD0Mass_(iConfig.getParameter<double>("minD0Mass")),
      maxD0Mass_(iConfig.getParameter<double>("maxD0Mass")),
      minDstMass_(iConfig.getParameter<double>("minDstMass")),
      maxDstMass_(iConfig.getParameter<double>("maxDstMass")),
      minDeltaMass_(iConfig.getParameter<double>("minDeltaMass")),
      maxDeltaMass_(iConfig.getParameter<double>("maxDeltaMass")),
      pvalMin_(iConfig.getParameter<double>("pvalMin")),
      applyChargeFilter_(iConfig.getParameter<bool>("applyChargeFilter")),
      useUnsignedCharge_(iConfig.getParameter<bool>("useUnsignedCharge")),
      targetCharge_(iConfig.getParameter<int>("charge")) {
  produces<reco::TrackCollection>();
}

bool DstToD0PiCandidateProducer::checkCharge(const reco::Track* kaon,
                                             const reco::Track* pion,
                                             const reco::Track* softPion) const {
  if (!applyChargeFilter_) {
    return true;
  }

  const bool validD0Charges = (kaon->charge() + pion->charge() == 0);
  const bool validSoftCharge = (pion->charge() == softPion->charge());
  if (!validD0Charges || !validSoftCharge) {
    return false;
  }

  int totalCharge = kaon->charge() + pion->charge() + softPion->charge();
  if (useUnsignedCharge_) {
    totalCharge = std::abs(totalCharge);
  }

  return totalCharge == targetCharge_;
}

FitQuality DstToD0PiCandidateProducer::fitQuality(RefCountedKinematicTree tree) const {
  FitQuality quality;
  if (!tree || !tree->isValid()) {
    return quality;
  }

  tree->movePointerToTheTop();
  RefCountedKinematicVertex vertex = tree->currentDecayVertex();
  if (!vertex || !vertex->vertexIsValid()) {
    return quality;
  }

  quality.isValid = true;
  quality.chi2 = vertex->chiSquared();
  quality.dof = vertex->degreesOfFreedom();
  quality.pval = TMath::Prob(quality.chi2, quality.dof);
  quality.isGood = quality.pval > pvalMin_;
  return quality;
}

RefCountedKinematicTree DstToD0PiCandidateProducer::fitD0(const edm::EventSetup& iSetup,
                                                          const reco::Track& kaon,
                                                          const reco::Track& pion) const {
  auto ttBuilder = iSetup.getHandle(transTrackBuilderToken_);

  reco::TransientTrack kaonTT = ttBuilder->build(kaon);
  reco::TransientTrack pionTT = ttBuilder->build(pion);

  KinematicParticleFactoryFromTransientTrack particleFactory;
  std::vector<RefCountedKinematicParticle> parts;
  float chi = 0.f;
  float ndf = 0.f;
  float kaonMassErr = kaonMassErr_;
  float pionMassErr = pionMassErr_;
  parts.push_back(particleFactory.particle(kaonTT, kaonMass_, chi, ndf, kaonMassErr));
  parts.push_back(particleFactory.particle(pionTT, pionMass_, chi, ndf, pionMassErr));

  KinematicParticleVertexFitter fitter;
  return fitter.fit(parts);
}

RefCountedKinematicTree DstToD0PiCandidateProducer::fitDst(const edm::EventSetup& iSetup,
                                                           const reco::Track& softPion,
                                                           const RefCountedKinematicParticle& d0) const {
  auto ttBuilder = iSetup.getHandle(transTrackBuilderToken_);

  reco::TransientTrack softPionTT = ttBuilder->build(softPion);
  reco::TransientTrack d0TT = d0->refittedTransientTrack();

  KinematicParticleFactoryFromTransientTrack particleFactory;
  std::vector<RefCountedKinematicParticle> parts;
  float chi = 0.f;
  float ndf = 0.f;
  float d0Mass = d0->currentState().mass();
  float d0MassErr = std::sqrt(std::abs(d0->currentState().kinematicParametersError().matrix()(6, 6)));
  float pionMassErr = pionMassErr_;
  parts.push_back(particleFactory.particle(softPionTT, softPionMass_, chi, ndf, pionMassErr));
  parts.push_back(particleFactory.particle(d0TT, d0Mass, chi, ndf, d0MassErr));

  KinematicParticleVertexFitter fitter;
  return fitter.fit(parts);
}

void DstToD0PiCandidateProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(trackToken_, tracks);

  auto selectedTracks = std::make_unique<reco::TrackCollection>();
  std::vector<CandidateSummary> candidates;

  if (tracks.isValid() && tracks->size() >= 3) {
    for (unsigned int iK = 0; iK < tracks->size(); ++iK) {
      const reco::Track* kaon = &(*tracks)[iK];
      // K and D0 pion are written to the output collection and refit
      // downstream by the CVH; loopers are dropped at the input stage so the
      // CVH propagator never has to handle them. The soft pion is only used
      // here for the D* mass / charge tag and is not in the output, so its
      // looper status is left as-is (the dM cut downstream filters mis-
      // measured soft pions on its own).
      if (kaon->isLooper()) {
        continue;
      }
      const int kaonCharge = kaon->charge();

      for (unsigned int iPi = 0; iPi < tracks->size(); ++iPi) {
        if (iPi == iK) {
          continue;
        }

        const reco::Track* pion = &(*tracks)[iPi];
        if (pion->isLooper()) {
          continue;
        }
        if (pion->charge() != -kaonCharge) {
          continue;
        }

        RefCountedKinematicTree d0Tree = fitD0(iSetup, *kaon, *pion);
        FitQuality d0Quality = fitQuality(d0Tree);
        if (!d0Quality.isGood) {
          continue;
        }

        d0Tree->movePointerToTheTop();
        RefCountedKinematicParticle d0 = d0Tree->currentParticle();
        const double d0Mass = d0->currentState().mass();
        if (d0Mass < minD0Mass_ || d0Mass > maxD0Mass_) {
          continue;
        }

        for (unsigned int iPis = 0; iPis < tracks->size(); ++iPis) {
          if (iPis == iK || iPis == iPi) {
            continue;
          }

          const reco::Track* softPion = &(*tracks)[iPis];
          if (softPion->charge() != pion->charge()) {
            continue;
          }
          if (!checkCharge(kaon, pion, softPion)) {
            continue;
          }

          RefCountedKinematicTree dstTree = fitDst(iSetup, *softPion, d0);
          FitQuality dstQuality = fitQuality(dstTree);
          if (!dstQuality.isGood) {
            continue;
          }

          dstTree->movePointerToTheTop();
          RefCountedKinematicParticle dst = dstTree->currentParticle();
          const double dstMass = dst->currentState().mass();
          const double deltaMass = dstMass - d0Mass;

          if (dstMass < minDstMass_ || dstMass > maxDstMass_) {
            continue;
          }
          if (deltaMass < minDeltaMass_ || deltaMass > maxDeltaMass_) {
            continue;
          }

          candidates.push_back(CandidateSummary{static_cast<float>(dst->currentState().globalMomentum().perp()),
                                                *kaon,
                                                *pion});
        }
      }
    }
  }

  std::sort(candidates.begin(), candidates.end(), [](const CandidateSummary& a, const CandidateSummary& b) {
    return a.dstPt > b.dstPt;
  });

  if (!candidates.empty()) {
    selectedTracks->push_back(candidates.front().kaon);
    selectedTracks->push_back(candidates.front().pion);
  }

  iEvent.put(std::move(selectedTracks));
}

DEFINE_FWK_MODULE(DstToD0PiCandidateProducer);
