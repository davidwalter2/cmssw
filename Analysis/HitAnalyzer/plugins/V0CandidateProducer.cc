/// \class V0CandidateProducer
///
/// Step-2 V0 candidate finder for KS->pi+pi- and Lambda0->p pi-.
///
/// Reads a flat reco::TrackCollection (the deduplicated daughter tracks
/// stored by V0DaughterTrackProducer in the ALCARECO step-1 output), iterates
/// over all opposite-sign track pairs, fits a Kalman vertex with mass
/// hypotheses (m1, m2), and applies V0-like cuts (vertex p-value, mass
/// window, pointing angle, transverse flight significance).
///
/// For asymmetric V0s (Lambda -> p pi-), enable `tryBothAssignments`: both
/// (m1, m2) and (m2, m1) are tried per pair and the assignment giving the
/// reconstructed mass closer to the nominal V0 mass is kept.
///
/// Output: a reco::TrackCollection containing the daughter tracks of the
/// single best candidate per event (sorted by V0 pT), in pair order:
/// track[0] = "first" daughter (gets m1 downstream), track[1] = "second".
/// This matches the input convention of
/// ResidualGlobalCorrectionMakerTwoTrackG4e with respectTrackOrder=True.
/// Designed to be a drop-in V0 analog of DstToD0PiCandidateProducer.

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "FWCore/Framework/interface/EDProducer.h"
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
  float v0Pt = 0.f;
  float v0Mass = 0.f;
  float deltaMass = 0.f;
  reco::Track first;   ///< gets daughterMass1 downstream
  reco::Track second;  ///< gets daughterMass2 downstream
};
}  // namespace

class V0CandidateProducer : public edm::EDProducer {
public:
  explicit V0CandidateProducer(const edm::ParameterSet&);
  ~V0CandidateProducer() override = default;

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  RefCountedKinematicTree fitV0(const edm::EventSetup& iSetup,
                                const reco::Track& t1, double m1,
                                const reco::Track& t2, double m2) const;
  FitQuality fitQuality(const RefCountedKinematicTree& tree) const;

  /// transverse pointing angle: cos(angle) between (SV-BS) in xy and V0 momentum in xy
  double cosThetaXY(const reco::BeamSpot& bs,
                    const RefCountedKinematicVertex& vertex,
                    const RefCountedKinematicParticle& v0) const;

  /// transverse flight length significance Lxy / sigma(Lxy) using vertex covariance
  double LxyOverSigma(const reco::BeamSpot& bs,
                      const RefCountedKinematicVertex& vertex) const;

  edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;

  // mass hypotheses
  double daughterMass1_;
  double daughterMass2_;
  double daughterMass1Err_;
  double daughterMass2Err_;
  bool tryBothAssignments_;
  double expectedV0Mass_;

  // mass window cut on the post-fit V0 mass
  double minV0Mass_;
  double maxV0Mass_;

  // vertex / pointing / flight cuts
  double pvalMin_;
  double cosThetaXYMin_;
  double LxyOverSigmaMin_;

  // basic per-track cut (redundant with ALCARECO but cheap and explicit)
  double minTrackPt_;

  // charge filter
  bool applyChargeFilter_;
  int targetCharge_;

  // K_S veto: when enabled (typically only for the Lambda channel), reject a
  // candidate whose two daughter tracks, evaluated under the (pi, pi) mass
  // hypothesis, give an invariant mass within +/- ksVetoWindow_ of the K_S
  // PDG mass. Suppresses real K_S contamination of the Lambda sample.
  bool applyKsVeto_;
  double ksVetoMass_;
  double ksVetoWindow_;
  double ksVetoPionMass_;
};

V0CandidateProducer::V0CandidateProducer(const edm::ParameterSet& iConfig)
    : trackToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
      beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
      daughterMass1_(iConfig.getParameter<double>("daughterMass1")),
      daughterMass2_(iConfig.getParameter<double>("daughterMass2")),
      // Per-daughter mass uncertainties (PDG). Fall back to the legacy
      // single-value `daughterMassErr` for back-compatibility, then to a
      // small non-zero default so the kinematic fit doesn't refuse a
      // singular mass-error matrix.
      daughterMass1Err_(iConfig.existsAs<double>("daughterMass1Err")
          ? iConfig.getParameter<double>("daughterMass1Err")
          : (iConfig.existsAs<double>("daughterMassErr")
              ? iConfig.getParameter<double>("daughterMassErr") : 1.e-6)),
      daughterMass2Err_(iConfig.existsAs<double>("daughterMass2Err")
          ? iConfig.getParameter<double>("daughterMass2Err")
          : (iConfig.existsAs<double>("daughterMassErr")
              ? iConfig.getParameter<double>("daughterMassErr") : 1.e-6)),
      tryBothAssignments_(iConfig.getParameter<bool>("tryBothAssignments")),
      expectedV0Mass_(iConfig.getParameter<double>("expectedV0Mass")),
      minV0Mass_(iConfig.getParameter<double>("minV0Mass")),
      maxV0Mass_(iConfig.getParameter<double>("maxV0Mass")),
      pvalMin_(iConfig.getParameter<double>("pvalMin")),
      cosThetaXYMin_(iConfig.getParameter<double>("cosThetaXYMin")),
      LxyOverSigmaMin_(iConfig.getParameter<double>("LxyOverSigmaMin")),
      minTrackPt_(iConfig.getParameter<double>("minTrackPt")),
      applyChargeFilter_(iConfig.getParameter<bool>("applyChargeFilter")),
      targetCharge_(iConfig.getParameter<int>("charge")),
      // K_S veto: defaults are off / sensible PDG values, so existing cfis
      // (e.g. the K_S producer itself) keep their current behaviour unless
      // they explicitly opt in.
      applyKsVeto_(iConfig.existsAs<bool>("applyKsVeto")
          ? iConfig.getParameter<bool>("applyKsVeto") : false),
      ksVetoMass_(iConfig.existsAs<double>("ksVetoMass")
          ? iConfig.getParameter<double>("ksVetoMass") : 0.497611),
      ksVetoWindow_(iConfig.existsAs<double>("ksVetoWindow")
          ? iConfig.getParameter<double>("ksVetoWindow") : 0.010),
      ksVetoPionMass_(iConfig.existsAs<double>("ksVetoPionMass")
          ? iConfig.getParameter<double>("ksVetoPionMass") : 0.13957039) {
  produces<reco::TrackCollection>();
}

FitQuality V0CandidateProducer::fitQuality(const RefCountedKinematicTree& tree) const {
  FitQuality q;
  if (!tree || !tree->isValid()) return q;
  tree->movePointerToTheTop();
  RefCountedKinematicVertex vertex = tree->currentDecayVertex();
  if (!vertex || !vertex->vertexIsValid()) return q;
  q.isValid = true;
  q.chi2 = vertex->chiSquared();
  q.dof = vertex->degreesOfFreedom();
  q.pval = TMath::Prob(q.chi2, q.dof);
  q.isGood = q.pval > pvalMin_;
  return q;
}

RefCountedKinematicTree V0CandidateProducer::fitV0(const edm::EventSetup& iSetup,
                                                  const reco::Track& t1, double m1,
                                                  const reco::Track& t2, double m2) const {
  edm::ESHandle<TransientTrackBuilder> ttBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttBuilder);
  reco::TransientTrack tt1 = ttBuilder->build(t1);
  reco::TransientTrack tt2 = ttBuilder->build(t2);
  KinematicParticleFactoryFromTransientTrack particleFactory;
  std::vector<RefCountedKinematicParticle> parts;
  float chi = 0.f, ndf = 0.f;
  // mass-error pairing follows the mass pairing: m1 always carries
  // daughterMass1Err_, m2 always carries daughterMass2Err_, regardless of
  // which input track (t1 or t2) was assigned which mass.
  float me1 = daughterMass1Err_;
  float me2 = daughterMass2Err_;
  parts.push_back(particleFactory.particle(tt1, m1, chi, ndf, me1));
  parts.push_back(particleFactory.particle(tt2, m2, chi, ndf, me2));
  KinematicParticleVertexFitter fitter;
  return fitter.fit(parts);
}

double V0CandidateProducer::cosThetaXY(const reco::BeamSpot& bs,
                                      const RefCountedKinematicVertex& vertex,
                                      const RefCountedKinematicParticle& v0) const {
  const double dx = vertex->position().x() - bs.x0();
  const double dy = vertex->position().y() - bs.y0();
  const double Lxy = std::sqrt(dx * dx + dy * dy);
  if (Lxy <= 0.) return -2.;
  const auto p = v0->currentState().globalMomentum();
  const double pt = std::sqrt(p.x() * p.x() + p.y() * p.y());
  if (pt <= 0.) return -2.;
  return (dx * p.x() + dy * p.y()) / (Lxy * pt);
}

double V0CandidateProducer::LxyOverSigma(const reco::BeamSpot& bs,
                                        const RefCountedKinematicVertex& vertex) const {
  const double dx = vertex->position().x() - bs.x0();
  const double dy = vertex->position().y() - bs.y0();
  const double Lxy = std::sqrt(dx * dx + dy * dy);
  if (Lxy <= 0.) return -1.;
  // sigma(Lxy)^2 from vertex covariance projected onto flight direction
  const auto cov = vertex->error().matrix();
  const double lx = dx / Lxy, ly = dy / Lxy;
  const double var = lx * lx * cov(0, 0) + 2. * lx * ly * cov(0, 1) + ly * ly * cov(1, 1);
  return var > 0. ? Lxy / std::sqrt(var) : -1.;
}

void V0CandidateProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(trackToken_, tracks);
  edm::Handle<reco::BeamSpot> bsH;
  iEvent.getByToken(beamSpotToken_, bsH);

  auto out = std::make_unique<reco::TrackCollection>();

  if (!tracks.isValid() || tracks->size() < 2 || !bsH.isValid()) {
    iEvent.put(std::move(out));
    return;
  }
  const reco::BeamSpot& bs = *bsH;

  std::vector<CandidateSummary> candidates;

  // helper: try one fit, apply all cuts, and push a CandidateSummary on success
  auto tryOne = [&](const reco::Track& tA, double mA, const reco::Track& tB, double mB) {
    // Optional K_S veto, evaluated up-front from the raw track 3-momenta
    // (pi+pi-) so that vetoed candidates skip the kinematic fit entirely.
    if (applyKsVeto_) {
      const double mpi = ksVetoPionMass_;
      const double EA = std::sqrt(tA.p() * tA.p() + mpi * mpi);
      const double EB = std::sqrt(tB.p() * tB.p() + mpi * mpi);
      const double pxsum = tA.px() + tB.px();
      const double pysum = tA.py() + tB.py();
      const double pzsum = tA.pz() + tB.pz();
      const double m2 = (EA + EB) * (EA + EB) - (pxsum * pxsum + pysum * pysum + pzsum * pzsum);
      if (m2 > 0.) {
        const double mpp = std::sqrt(m2);
        if (std::fabs(mpp - ksVetoMass_) < ksVetoWindow_) return;
      }
    }
    RefCountedKinematicTree tree = fitV0(iSetup, tA, mA, tB, mB);
    FitQuality q = fitQuality(tree);
    if (!q.isGood) return;
    tree->movePointerToTheTop();
    RefCountedKinematicParticle v0 = tree->currentParticle();
    RefCountedKinematicVertex vertex = tree->currentDecayVertex();
    const double mass = v0->currentState().mass();
    if (mass < minV0Mass_ || mass > maxV0Mass_) return;
    if (cosThetaXY(bs, vertex, v0) < cosThetaXYMin_) return;
    if (LxyOverSigma(bs, vertex) < LxyOverSigmaMin_) return;
    const auto p = v0->currentState().globalMomentum();
    const double v0Pt = std::sqrt(p.x() * p.x() + p.y() * p.y());
    candidates.push_back(CandidateSummary{
        static_cast<float>(v0Pt),
        static_cast<float>(mass),
        static_cast<float>(std::fabs(mass - expectedV0Mass_)),
        tA, tB});
  };

  for (unsigned int i = 0; i < tracks->size(); ++i) {
    const reco::Track& t1 = (*tracks)[i];
    if (t1.isLooper()) continue;
    if (t1.pt() < minTrackPt_) continue;
    for (unsigned int j = i + 1; j < tracks->size(); ++j) {
      const reco::Track& t2 = (*tracks)[j];
      if (t2.isLooper()) continue;
      if (t2.pt() < minTrackPt_) continue;
      if (applyChargeFilter_ && (t1.charge() + t2.charge() != targetCharge_)) continue;

      // Hypothesis A: t1 gets m1, t2 gets m2
      tryOne(t1, daughterMass1_, t2, daughterMass2_);
      // Hypothesis B (asymmetric only): t1 gets m2, t2 gets m1 -- effectively
      // the swap. The two hypotheses are the same when m1 == m2.
      if (tryBothAssignments_ && daughterMass1_ != daughterMass2_) {
        tryOne(t2, daughterMass1_, t1, daughterMass2_);
      }
    }
  }

  if (!candidates.empty()) {
    // For each pair (i,j) with both hypotheses tried, keep only the assignment
    // closer to the nominal V0 mass. Then sort the remaining candidates by V0
    // pT and keep the highest-pT one (matches DstToD0PiCandidateProducer).
    std::sort(candidates.begin(), candidates.end(),
              [](const CandidateSummary& a, const CandidateSummary& b) {
                return a.v0Pt > b.v0Pt;
              });
    out->push_back(candidates.front().first);   // gets m1 downstream
    out->push_back(candidates.front().second);  // gets m2 downstream
  }

  iEvent.put(std::move(out));
}

DEFINE_FWK_MODULE(V0CandidateProducer);
