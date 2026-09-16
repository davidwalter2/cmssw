// WMass: PUPPI weights and impact parameters of the packed PF candidates with
// respect to a "robust" primary vertex -- the reconstructed vertex closest in dz
// to the leading loose muon (pt > muonPtMin), or, when no vertex is within
// muonVertexDzMax, a pseudo-vertex at the beam spot with the muon's z. The
// hard-scatter vertex of a W -> mu nu event is then never lost to a mis-ranked
// PV. Ported from WmassNanoProd_10_6_26 (PR #33, PV-robust DeepMET): the 10_6
// plugin was a copy of the PuppiProducer of its time with this choice added; this
// one is the packed-candidate path of the 15_0 PuppiProducer with the same
// additions, but it runs the 10_6 PuppiContainer (private copy, puppi106): the
// PV-robust DeepMET models were trained on those weights, and the release
// container computes different ones for the same inputs. Products:
//   PVRobustIndex (int)                      the chosen vertex, -1 = beam-spot pseudo-vertex
//   PVMuonIndex (int)                        the muon in the input collection, -1 = none
//                                            (10_6 wrote the last index when no muon passed)
//   ValueMap<float>                          PUPPI weight per input candidate
//   pat::PackedCandidateCollection           the inputs with the new PUPPI weight (p4 untouched)
//   PFPVRobustDxy, PFPVRobustDz (ValueMap<double>)  w.r.t. the robust vertex, 0 for neutrals
//   PFPVRobustPuppiWeight (ValueMap<double>)
#include "PuppiAlgo106.h"
#include "PuppiContainer106.h"
#include "RecoObj106.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <memory>

class PuppiPVRobustProducer : public edm::stream::EDProducer<> {
public:
  explicit PuppiPVRobustProducer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  typedef reco::VertexCollection VertexCollection;
  typedef edm::View<reco::Candidate> CandidateView;
  typedef std::vector<pat::Muon> MuonCollection;

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  edm::EDGetTokenT<CandidateView> tokenPFCandidates_;
  edm::EDGetTokenT<VertexCollection> tokenVertices_;
  edm::EDGetTokenT<MuonCollection> tokenMuons_;
  edm::EDGetTokenT<reco::BeamSpot> tokenBeamSpot_;
  edm::EDPutTokenT<int> ptokenPVIndex_;
  edm::EDPutTokenT<int> ptokenMuonIndex_;
  edm::EDPutTokenT<edm::ValueMap<float>> ptokenPupOut_;
  edm::EDPutTokenT<pat::PackedCandidateCollection> ptokenPackedPuppiCandidates_;
  edm::EDPutTokenT<edm::ValueMap<double>> ptokenDxy_;
  edm::EDPutTokenT<edm::ValueMap<double>> ptokenDz_;
  edm::EDPutTokenT<edm::ValueMap<double>> ptokenPuppiWeight_;

  bool fPuppiNoLep;
  bool fUseFromPVLooseTight;
  bool fUseDZ;
  double fDZCut;
  double fEtaMinUseDZ;
  double fPtMaxCharged;
  double fEtaMaxCharged;
  uint fNumOfPUVtxsForCharged;
  double fDZCutForChargedFromPUVtxs;
  double fMuonPtMin;
  double fMuonVertexDzMax;
  int fVtxNdofCut;
  double fVtxZCut;
  std::unique_ptr<puppi106::PuppiContainer> fPuppiContainer;
};

// ------------------------------------------------------------------------------------------
PuppiPVRobustProducer::PuppiPVRobustProducer(const edm::ParameterSet& iConfig)
    : fPuppiContainer(std::make_unique<puppi106::PuppiContainer>(iConfig)) {
  fPuppiNoLep = iConfig.getParameter<bool>("puppiNoLep");
  fUseFromPVLooseTight = iConfig.getParameter<bool>("UseFromPVLooseTight");
  fUseDZ = iConfig.getParameter<bool>("UseDeltaZCut");
  fDZCut = iConfig.getParameter<double>("DeltaZCut");
  fEtaMinUseDZ = iConfig.getParameter<double>("EtaMinUseDeltaZ");
  fPtMaxCharged = iConfig.getParameter<double>("PtMaxCharged");
  fEtaMaxCharged = iConfig.getParameter<double>("EtaMaxCharged");
  fNumOfPUVtxsForCharged = iConfig.getParameter<uint>("NumOfPUVtxsForCharged");
  fDZCutForChargedFromPUVtxs = iConfig.getParameter<double>("DeltaZCutForChargedFromPUVtxs");
  fMuonPtMin = iConfig.getParameter<double>("muonPtMin");
  fMuonVertexDzMax = iConfig.getParameter<double>("muonVertexDzMax");
  fVtxNdofCut = iConfig.getParameter<int>("vtxNdofCut");
  fVtxZCut = iConfig.getParameter<double>("vtxZCut");

  tokenPFCandidates_ = consumes<CandidateView>(iConfig.getParameter<edm::InputTag>("candName"));
  tokenVertices_ = consumes<VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexName"));
  tokenMuons_ = consumes<MuonCollection>(iConfig.getParameter<edm::InputTag>("muonName"));
  tokenBeamSpot_ = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotName"));

  ptokenPVIndex_ = produces<int>("PVRobustIndex");
  ptokenMuonIndex_ = produces<int>("PVMuonIndex");
  ptokenPupOut_ = produces<edm::ValueMap<float>>();
  ptokenPackedPuppiCandidates_ = produces<pat::PackedCandidateCollection>();
  ptokenDxy_ = produces<edm::ValueMap<double>>("PFPVRobustDxy");
  ptokenDz_ = produces<edm::ValueMap<double>>("PFPVRobustDz");
  ptokenPuppiWeight_ = produces<edm::ValueMap<double>>("PFPVRobustPuppiWeight");
}

// ------------------------------------------------------------------------------------------
void PuppiPVRobustProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<CandidateView> hPFProduct;
  iEvent.getByToken(tokenPFCandidates_, hPFProduct);
  const reco::VertexCollection& pvCol = iEvent.get(tokenVertices_);
  const MuonCollection& muonCol = iEvent.get(tokenMuons_);
  const reco::BeamSpot& beamSpot = iEvent.get(tokenBeamSpot_);
  const math::XYZPoint beamPoint(beamSpot.x0(), beamSpot.y0(), beamSpot.z0());

  // the robust vertex: the one closest in dz to the leading loose muon;
  // without such a muon the leading PV, as in the standard PUPPI
  int iLV = 0;
  int iMuon = -1;
  double muonz = 0;
  for (size_t im = 0; im < muonCol.size(); ++im) {
    const pat::Muon& muon = muonCol[im];
    if (muon.pt() < fMuonPtMin || !muon.isLooseMuon())
      continue;
    int iVClosest = -1;
    double dzClosest = fMuonVertexDzMax;
    for (size_t iv = 0; iv < pvCol.size(); ++iv) {
      const double dz = std::abs(muon.muonBestTrack()->dz(pvCol[iv].position()));
      if (dz < dzClosest) {
        dzClosest = dz;
        iVClosest = iv;
      }
    }
    iLV = iVClosest;
    iMuon = im;
    // if the muon is not close to any vertex, a pseudo-vertex at the beam spot with the muon's z
    muonz = muon.muonBestTrack()->dz(beamPoint) + beamPoint.z();
    break;
  }
  const math::XYZPoint pvPoint =
      iLV >= 0 ? math::XYZPoint(pvCol[iLV].position()) : math::XYZPoint(beamSpot.x0(), beamSpot.y0(), muonz);
  if (iLV != 0)
    LogDebug("PuppiPVRobustProducer") << "run " << iEvent.run() << " lumi " << iEvent.luminosityBlock() << " event "
                                      << iEvent.id().event() << ": robust PV " << iLV << " (muon " << iMuon << ", z "
                                      << pvPoint.z() << ", PV0 z " << (pvCol.empty() ? 0. : pvCol[0].z()) << ")";

  int npv = 0;
  for (auto const& vtx : pvCol) {
    if (!vtx.isFake() && vtx.ndof() >= fVtxNdofCut && std::abs(vtx.z()) <= fVtxZCut)
      ++npv;
  }

  const size_t nCand = hPFProduct->size();
  std::vector<puppi106::RecoObj> recoObjCollection;
  recoObjCollection.reserve(nCand);
  std::vector<double> dxyVals(nCand, 0.);
  std::vector<double> dzVals(nCand, 0.);
  for (size_t ic = 0; ic < nCand; ++ic) {
    const reco::Candidate& aPF = (*hPFProduct)[ic];
    const pat::PackedCandidate* lPack = dynamic_cast<const pat::PackedCandidate*>(&aPF);
    if (lPack == nullptr)
      throw edm::Exception(edm::errors::LogicError, "PuppiPVRobustProducer: inputs are not PackedCandidates");
    puppi106::RecoObj pReco;
    pReco.pt = aPF.pt();
    pReco.eta = aPF.eta();
    pReco.phi = aPF.phi();
    pReco.m = aPF.mass();
    pReco.rapidity = aPF.rapidity();
    pReco.charge = aPF.charge();
    pReco.id = 0;  // 0: to be calculated, 1: from PV, 2: from PU, 3: lepton left out
    pReco.dZ = 0;
    pReco.d0 = 0;
    const bool isLepton = ((std::abs(aPF.pdgId()) == 11) || (std::abs(aPF.pdgId()) == 13));

    if (lPack->vertexRef().isNonnull() && std::abs(pReco.charge) > 0) {
      // impact parameters w.r.t. the robust vertex, not the candidate's own vertex
      const double pDZ = lPack->dz(pvPoint);
      const double pD0 = lPack->dxy(pvPoint);
      pReco.dZ = pDZ;
      pReco.d0 = pD0;
      dzVals[ic] = pDZ;
      dxyVals[ic] = pD0;

      if (fPuppiNoLep && isLepton) {
        pReco.id = 3;
      } else if (iLV >= 0) {
        // the standard assignment, with the association to the robust vertex
        const int fromPV = lPack->fromPV(iLV);
        if (fromPV == 0) {
          pReco.id = 2;
          if ((fNumOfPUVtxsForCharged > 0) and (std::abs(pDZ) < fDZCutForChargedFromPUVtxs)) {
            // vertex-splitting recovery; from vertex 0, the robust vertex need not be the leading one
            for (size_t puVtx_idx = 0; puVtx_idx <= (fNumOfPUVtxsForCharged + 1) && puVtx_idx < pvCol.size();
                 ++puVtx_idx) {
              if (lPack->fromPV(puVtx_idx) >= 2) {
                pReco.id = 1;
                break;
              }
            }
          }
        } else if (fromPV == pat::PackedCandidate::PVUsedInFit) {
          pReco.id = 1;
        } else if (fromPV == pat::PackedCandidate::PVTight || fromPV == pat::PackedCandidate::PVLoose) {
          pReco.id = 0;
          if ((fPtMaxCharged > 0) and (pReco.pt > fPtMaxCharged))
            pReco.id = 1;
          else if (std::abs(pReco.eta) > fEtaMaxCharged)
            pReco.id = 1;
          else if ((fUseDZ) && (std::abs(pReco.eta) >= fEtaMinUseDZ))
            pReco.id = (std::abs(pDZ) < fDZCut) ? 1 : 2;
          else if (fUseFromPVLooseTight && fromPV == pat::PackedCandidate::PVLoose)
            pReco.id = 2;
          else if (fUseFromPVLooseTight && fromPV == pat::PackedCandidate::PVTight)
            pReco.id = 1;
        }
      } else {
        // pseudo-vertex: no vertex association to use, dz alone decides
        if ((fPtMaxCharged > 0) and (pReco.pt > fPtMaxCharged))
          pReco.id = 1;
        else if (std::abs(pReco.eta) > fEtaMaxCharged)
          pReco.id = 1;
        else
          pReco.id = (std::abs(pDZ) < fDZCut) ? 1 : 2;
      }
    }
    recoObjCollection.push_back(pReco);
  }

  fPuppiContainer->initialize(recoObjCollection);
  fPuppiContainer->setNPV(npv);
  const std::vector<double> lWeights = fPuppiContainer->puppiWeights();  // one per input candidate
  const std::vector<int>& recoToPup = fPuppiContainer->recoToPup();      // -1: not a puppi particle

  edm::ValueMap<float> lPupOut;
  edm::ValueMap<float>::Filler lPupFiller(lPupOut);
  lPupFiller.insert(hPFProduct, lWeights.begin(), lWeights.end());
  lPupFiller.fill();

  // the inputs with the new PUPPI weight; the four-momenta are left as they are
  pat::PackedCandidateCollection packedPuppiCandidates;
  packedPuppiCandidates.reserve(nCand);
  for (size_t ic = 0; ic < nCand; ++ic) {
    const reco::Candidate& aCand = (*hPFProduct)[ic];
    pat::PackedCandidate pCand(*dynamic_cast<const pat::PackedCandidate*>(&aCand));
    if (fPuppiNoLep)
      pCand.setPuppiWeight(pCand.puppiWeight(), lWeights[ic]);
    else
      pCand.setPuppiWeight(lWeights[ic], pCand.puppiWeightNoLep());
    pCand.setSourceCandidatePtr(aCand.sourceCandidatePtr(0));
    packedPuppiCandidates.push_back(pCand);
  }

  auto fillDouble = [&](const std::vector<double>& vals) {
    edm::ValueMap<double> map;
    edm::ValueMap<double>::Filler filler(map);
    filler.insert(hPFProduct, vals.begin(), vals.end());
    filler.fill();
    return map;
  };
  // as in 10_6: the DeepMET input map carries 0 for candidates the container dropped
  std::vector<double> puppiWeights(nCand, 0.);
  for (size_t ic = 0; ic < nCand; ++ic)
    puppiWeights[ic] = recoToPup[ic] >= 0 ? lWeights[ic] : 0.;

  iEvent.emplace(ptokenPVIndex_, iLV);
  iEvent.emplace(ptokenMuonIndex_, iMuon);
  iEvent.emplace(ptokenPupOut_, std::move(lPupOut));
  iEvent.emplace(ptokenPackedPuppiCandidates_, std::move(packedPuppiCandidates));
  iEvent.emplace(ptokenDxy_, fillDouble(dxyVals));
  iEvent.emplace(ptokenDz_, fillDouble(dzVals));
  iEvent.emplace(ptokenPuppiWeight_, fillDouble(puppiWeights));
}

// ------------------------------------------------------------------------------------------
void PuppiPVRobustProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  // the PuppiProducer / PuppiContainer parameters
  desc.add<bool>("puppiDiagnostics", false);
  desc.add<bool>("puppiNoLep", false);
  desc.add<bool>("UseFromPVLooseTight", false);
  desc.add<bool>("UseDeltaZCut", true);
  desc.add<double>("DeltaZCut", 0.3);
  desc.add<double>("EtaMinUseDeltaZ", 0.);
  desc.add<double>("PtMaxCharged", -1.);
  desc.add<double>("EtaMaxCharged", 99999.);
  desc.add<double>("PtMaxPhotons", -1.);
  desc.add<double>("EtaMaxPhotons", 2.5);
  desc.add<double>("PtMaxNeutrals", 200.);
  desc.add<double>("PtMaxNeutralsStartSlope", 0.);
  desc.add<uint>("NumOfPUVtxsForCharged", 0);
  desc.add<double>("DeltaZCutForChargedFromPUVtxs", 0.2);
  desc.add<int>("vtxNdofCut", 4);
  desc.add<double>("vtxZCut", 24);
  desc.add<edm::InputTag>("candName", edm::InputTag("packedPFCandidates"));
  desc.add<edm::InputTag>("vertexName", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<bool>("applyCHS", true);
  desc.add<bool>("invertPuppi", false);
  desc.add<bool>("useExp", false);
  desc.add<double>("MinPuppiWeight", .01);
  // the robust vertex
  desc.add<edm::InputTag>("muonName", edm::InputTag("slimmedMuons"));
  desc.add<edm::InputTag>("beamSpotName", edm::InputTag("offlineBeamSpot"));
  desc.add<double>("muonPtMin", 10.)->setComment("the leading loose muon above this pt chooses the vertex");
  desc.add<double>("muonVertexDzMax", 0.2)
      ->setComment("max |dz(muon, vertex)| in cm; beyond it a pseudo-vertex at the beam spot with the muon's z");

  puppi106::PuppiAlgo::fillDescriptionsPuppiAlgo(desc);

  descriptions.add("PuppiPVRobustProducer", desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(PuppiPVRobustProducer);
