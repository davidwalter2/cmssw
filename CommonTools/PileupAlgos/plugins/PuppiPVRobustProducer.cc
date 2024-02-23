// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/Common/interface/Association.h"
//Main File
#include "CommonTools/PileupAlgos/plugins/PuppiPVRobustProducer.h"
#include "CommonTools/PileupAlgos/interface/PuppiCandidate.h"

// ------------------------------------------------------------------------------------------
PuppiPVRobustProducer::PuppiPVRobustProducer(const edm::ParameterSet& iConfig) {
  fPuppiDiagnostics = iConfig.getParameter<bool>("puppiDiagnostics");
  fPuppiForLeptons = iConfig.getParameter<bool>("puppiForLeptons");
  fUseFromPVLooseTight = iConfig.getParameter<bool>("UseFromPVLooseTight");
  fUseDZ = iConfig.getParameter<bool>("UseDeltaZCut");
  fDZCut = iConfig.getParameter<double>("DeltaZCut");
  fEtaMinUseDZ = iConfig.getParameter<double>("EtaMinUseDeltaZ");
  fPtMaxCharged = iConfig.getParameter<double>("PtMaxCharged");
  fEtaMaxCharged = iConfig.getParameter<double>("EtaMaxCharged");
  fNumOfPUVtxsForCharged = iConfig.getParameter<uint>("NumOfPUVtxsForCharged");
  fDZCutForChargedFromPUVtxs = iConfig.getParameter<double>("DeltaZCutForChargedFromPUVtxs");
  fUseExistingWeights = iConfig.getParameter<bool>("useExistingWeights");
  fUseWeightsNoLep = iConfig.getParameter<bool>("useWeightsNoLep");
  fClonePackedCands = iConfig.getParameter<bool>("clonePackedCands");
  fVtxNdofCut = iConfig.getParameter<int>("vtxNdofCut");
  fVtxZCut = iConfig.getParameter<double>("vtxZCut");
  fPuppiContainer = std::unique_ptr<PuppiContainer>(new PuppiContainer(iConfig));

  tokenPFCandidates_ = consumes<CandidateView>(iConfig.getParameter<edm::InputTag>("candName"));
  tokenVertices_ = consumes<VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexName"));
  tokenMuons_ = consumes<MuonCollection>(iConfig.getParameter<edm::InputTag>("muonName"));
  offlinebeamSpot_ = consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));

  produces<edm::ValueMap<float>>();
  produces<edm::ValueMap<LorentzVector>>();
  produces<edm::ValueMap<reco::CandidatePtr>>();

  if (fUseExistingWeights || fClonePackedCands)
    produces<pat::PackedCandidateCollection>();
  else
    produces<reco::PFCandidateCollection>();

  if (fPuppiDiagnostics) {
    produces<double>("PuppiNAlgos");
    produces<std::vector<double>>("PuppiRawAlphas");
    produces<std::vector<double>>("PuppiAlphas");
    produces<std::vector<double>>("PuppiAlphasMed");
    produces<std::vector<double>>("PuppiAlphasRms");
  }
}
// ------------------------------------------------------------------------------------------
PuppiPVRobustProducer::~PuppiPVRobustProducer() {}
// ------------------------------------------------------------------------------------------
void PuppiPVRobustProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Get PFCandidate Collection
  edm::Handle<CandidateView> hPFProduct;
  iEvent.getByToken(tokenPFCandidates_, hPFProduct);
  const CandidateView* pfCol = hPFProduct.product();

  // Get vertex collection w/PV as the first entry?
  edm::Handle<reco::VertexCollection> hVertexProduct;
  iEvent.getByToken(tokenVertices_, hVertexProduct);
  const reco::VertexCollection* pvCol = hVertexProduct.product();

  // Get muon collection
  edm::Handle<MuonCollection> hMuonProduct;
  iEvent.getByToken(tokenMuons_, hMuonProduct);
  const MuonCollection* muonCol = hMuonProduct.product();

  // Get BeamSpot
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(offlinebeamSpot_, beamSpotHandle);

  if (beamSpotHandle.isValid()) {
    beamSpot = *beamSpotHandle;
  } else {
    edm::LogError("PuppiPVRobustProducer") << "No beam spot available from EventSetup"
                                           << "\n...skip event";
    return;
  }
  const math::XYZPoint beamPoint(beamSpot.x0(), beamSpot.y0(), beamSpot.z0());

  // leading vertex closest to the "leading" muon passing the loose ID
  int iLV = 0;
  double muonz = 0;
  //double pv_z = 0;
  for (auto const& muon : *muonCol) {
    if (muon.pt() < 5)
      continue;
    if (!muon.isLooseMuon())
      continue;

    //std::cout << "muonz " << muon.vz() << " BestTrack_z " << muon.muonBestTrack()->vz() << std::endl;
    // is this the best way to get the z position of the muon?
    muonz = muon.vz();

    int iVClosest = -1;
    int iV = 0;
    // minimum dz requirement for the closest vertex
    double dzClosest = 0.2;
    for (auto const& aV : *pvCol) {
      // is this the best way to get the dz between the vertex and the leading muon?
      if (fabs(muon.muonBestTrack()->dz(aV.position())) < dzClosest) {
        dzClosest = fabs(muon.muonBestTrack()->dz(aV.position()));
        iVClosest = iV;

        //pv_z = aV.z();
      }
      iV++;
    }
    iLV = iVClosest;
    break;
  }
  //std::cout << "Leading vertex closest to the leading muon passing the loose ID: " << iLV << " muonz " << muonz
  //          << " pv_z " << pv_z << std::endl;

  math::XYZPoint pvPoint;
  if (iLV >= 0) {
    auto const& aLV = pvCol->at(iLV);
    pvPoint = math::XYZPoint(aLV.x(), aLV.y(), aLV.z());
  } else {
    // failed to find a vertex close to the leading muon within 0.2 cm
    // use the beamspot and muon z position to make a "fake" vertex
    pvPoint = math::XYZPoint(beamSpot.x0(), beamSpot.y0(), muonz);
  }
  //auto const& aLV = *pvCol->begin();
  //const math::XYZPoint pvPoint(aLV.x(), aLV.y(), aLV.z());

  int npv = 0;
  const reco::VertexCollection::const_iterator vtxEnd = pvCol->end();
  for (reco::VertexCollection::const_iterator vtxIter = pvCol->begin(); vtxEnd != vtxIter; ++vtxIter) {
    if (!vtxIter->isFake() && vtxIter->ndof() >= fVtxNdofCut && std::abs(vtxIter->z()) <= fVtxZCut)
      npv++;
  }

  //Fill the reco objects
  fRecoObjCollection.clear();
  fRecoObjCollection.reserve(pfCol->size());
  for (auto const& aPF : *pfCol) {
    RecoObj pReco;
    pReco.pt = aPF.pt();
    pReco.eta = aPF.eta();
    pReco.phi = aPF.phi();
    pReco.m = aPF.mass();
    pReco.rapidity = aPF.rapidity();
    pReco.charge = aPF.charge();

    // not sure why we need this "if" here
    if (aPF.vertexRef().isNonnull()) {
      // id: 0 to be calculated; 1: from PV; 2: from PU
      pReco.id = 0;

      if (std::abs(pReco.charge) == 0) {
        // for neutral particles, d0 and dZ are always 0
        pReco.dZ = 0;
        pReco.d0 = 0;
      } else {
        double pDZ = aPF.dz(pvPoint);
        double pD0 = aPF.dxy(pvPoint);

        pReco.dZ = pDZ;
        // i dont think d0 is used anywhere though?
        pReco.d0 = pD0;

        if (iLV >= 0) {
          // leading vertex closest to the leading muon exists/reconstructed
          if (aPF.fromPV(iLV) == 0) {
            pReco.id = 2;
            if ((fNumOfPUVtxsForCharged > 0) and (std::abs(pDZ) < fDZCutForChargedFromPUVtxs)) {
              // for vertex splitting case
              for (size_t puVtx_idx = 0; puVtx_idx <= (fNumOfPUVtxsForCharged + 1) && puVtx_idx < pvCol->size();
                   ++puVtx_idx) {
                // loop from 0th now, since the iLV is not necessarily the 0th vertex
                if (aPF.fromPV(puVtx_idx) >= 2) {
                  pReco.id = 1;
                  break;
                }
              }
            }
          } else if (aPF.fromPV(iLV) == (pat::PackedCandidate::PVUsedInFit)) {
            pReco.id = 1;
          } else if (aPF.fromPV(iLV) == (pat::PackedCandidate::PVTight) ||
                     aPF.fromPV(iLV) == (pat::PackedCandidate::PVLoose)) {
            pReco.id = 0;
            if ((fPtMaxCharged > 0) and (pReco.pt > fPtMaxCharged))
              pReco.id = 1;
            else if (std::abs(pReco.eta) > fEtaMaxCharged)
              pReco.id = 1;
            else if ((fUseDZ) && (std::abs(pReco.eta) >= fEtaMinUseDZ))
              pReco.id = (std::abs(pDZ) < fDZCut) ? 1 : 2;
            else if (fUseFromPVLooseTight && aPF.fromPV(iLV) == (pat::PackedCandidate::PVLoose))
              pReco.id = 2;
            else if (fUseFromPVLooseTight && aPF.fromPV(iLV) == (pat::PackedCandidate::PVTight))
              pReco.id = 1;
          }
        } else {
          // no PV closest to the leading muon within 0.2cm
          // use dZ as the discriminator for determining LV or PU
          if ((fPtMaxCharged > 0) and (pReco.pt > fPtMaxCharged))
            pReco.id = 1;
          else if (std::abs(pReco.eta) > fEtaMaxCharged)
            pReco.id = 1;
          else
            pReco.id = (std::abs(pDZ) < fDZCut) ? 1 : 2;
        }
      }

      //std::cout << "PF pt " << aPF.pt() << " eta " << aPF.eta() << " phi " << aPF.phi() << " charge " << aPF.charge()
      //          << " pdgId " << aPF.pdgId() << " d0 " << aPF.dxy() << " d0_BS " << aPF.dxy(beamPoint) << " d0_PV "
      //          << aPF.dxy(pvPoint) << " dz " << aPF.dz() << " dz_BS " << aPF.dz(beamPoint) << " dz_PV "
      //          << aPF.dz(pvPoint) << " id " << pReco.id << std::endl;
    }

    fRecoObjCollection.push_back(pReco);
  }

  assert(fRecoObjCollection.size() == pfCol->size());

  fPuppiContainer->initialize(fRecoObjCollection);
  fPuppiContainer->setNPV(npv);

  //Compute the weights and get the particles
  std::vector<double> lWeights = fPuppiContainer->puppiWeights();
  std::vector<PuppiCandidate> lCandidates = fPuppiContainer->puppiParticles();

  //Fill it into the event
  std::unique_ptr<edm::ValueMap<float>> lPupOut(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler lPupFiller(*lPupOut);
  lPupFiller.insert(hPFProduct, lWeights.begin(), lWeights.end());
  lPupFiller.fill();

  // This is a dummy to access the "translate" method which is a
  // non-static member function even though it doesn't need to be.
  // Will fix in the future.
  static const reco::PFCandidate dummySinceTranslateIsNotStatic;

  // Fill a new PF/Packed Candidate Collection and write out the ValueMap of the new p4s.
  // Since the size of the ValueMap must be equal to the input collection, we need
  // to search the "puppi" particles to find a match for each input. If none is found,
  // the input is set to have a four-vector of 0,0,0,0
  fPuppiCandidates.reset(new PFOutputCollection);
  fPackedPuppiCandidates.reset(new PackedOutputCollection);
  std::unique_ptr<edm::ValueMap<LorentzVector>> p4PupOut(new edm::ValueMap<LorentzVector>());
  LorentzVectorCollection puppiP4s;
  std::vector<reco::CandidatePtr> values(hPFProduct->size());

  int val = -1;
  for (auto const& aCand : *hPFProduct) {
    val++;
    std::unique_ptr<pat::PackedCandidate> pCand;
    std::unique_ptr<reco::PFCandidate> pfCand;
    if (fUseExistingWeights || fClonePackedCands) {
      const pat::PackedCandidate* cand = dynamic_cast<const pat::PackedCandidate*>(&aCand);
      if (!cand)
        throw edm::Exception(edm::errors::LogicError, "PuppiPVRobustProducer: inputs are not PackedCandidates");
      pCand.reset(new pat::PackedCandidate(*cand));
    } else {
      auto id = dummySinceTranslateIsNotStatic.translatePdgIdToType(aCand.pdgId());
      const reco::PFCandidate* cand = dynamic_cast<const reco::PFCandidate*>(&aCand);
      pfCand.reset(new reco::PFCandidate(cand ? *cand : reco::PFCandidate(aCand.charge(), aCand.p4(), id)));
    }
    LorentzVector pVec;

    //get an index to a pup in lCandidates: either fUseExistingWeights with no skips or get from fPuppiContainer
    int iPuppiMatched = fUseExistingWeights ? val : fPuppiContainer->recoToPup()[val];
    if (iPuppiMatched >= 0) {
      auto const& puppiMatched = lCandidates[iPuppiMatched];
      pVec.SetPxPyPzE(puppiMatched.px, puppiMatched.py, puppiMatched.pz, puppiMatched.e);
      if (fClonePackedCands && (!fUseExistingWeights)) {
        if (fPuppiForLeptons)
          pCand->setPuppiWeight(pCand->puppiWeight(), lWeights[val]);
        else
          pCand->setPuppiWeight(lWeights[val], pCand->puppiWeightNoLep());
      }
    } else {
      pVec.SetPxPyPzE(0, 0, 0, 0);
      if (fClonePackedCands && (!fUseExistingWeights)) {
        pCand->setPuppiWeight(0, 0);
      }
    }
    puppiP4s.push_back(pVec);

    if (fUseExistingWeights || fClonePackedCands) {
      // already have PUPPI weight, do not reset kinematic
      //pCand->setP4(pVec);
      pCand->setSourceCandidatePtr(aCand.sourceCandidatePtr(0));
      fPackedPuppiCandidates->push_back(*pCand);
    } else {
      pfCand->setP4(pVec);
      pfCand->setSourceCandidatePtr(aCand.sourceCandidatePtr(0));
      fPuppiCandidates->push_back(*pfCand);
    }
  }

  //Compute the modified p4s
  edm::ValueMap<LorentzVector>::Filler p4PupFiller(*p4PupOut);
  p4PupFiller.insert(hPFProduct, puppiP4s.begin(), puppiP4s.end());
  p4PupFiller.fill();

  iEvent.put(std::move(lPupOut));
  iEvent.put(std::move(p4PupOut));
  if (fUseExistingWeights || fClonePackedCands) {
    edm::OrphanHandle<pat::PackedCandidateCollection> oh = iEvent.put(std::move(fPackedPuppiCandidates));
    for (unsigned int ic = 0, nc = oh->size(); ic < nc; ++ic) {
      reco::CandidatePtr pkref(oh, ic);
      values[ic] = pkref;
    }
  } else {
    edm::OrphanHandle<reco::PFCandidateCollection> oh = iEvent.put(std::move(fPuppiCandidates));
    for (unsigned int ic = 0, nc = oh->size(); ic < nc; ++ic) {
      reco::CandidatePtr pkref(oh, ic);
      values[ic] = pkref;
    }
  }
  std::unique_ptr<edm::ValueMap<reco::CandidatePtr>> pfMap_p(new edm::ValueMap<reco::CandidatePtr>());
  edm::ValueMap<reco::CandidatePtr>::Filler filler(*pfMap_p);
  filler.insert(hPFProduct, values.begin(), values.end());
  filler.fill();
  iEvent.put(std::move(pfMap_p));

  //////////////////////////////////////////////
  if (fPuppiDiagnostics && !fUseExistingWeights) {
    // all the different alphas per particle
    // THE alpha per particle
    std::unique_ptr<std::vector<double>> theAlphas(new std::vector<double>(fPuppiContainer->puppiAlphas()));
    std::unique_ptr<std::vector<double>> theAlphasMed(new std::vector<double>(fPuppiContainer->puppiAlphasMed()));
    std::unique_ptr<std::vector<double>> theAlphasRms(new std::vector<double>(fPuppiContainer->puppiAlphasRMS()));
    std::unique_ptr<std::vector<double>> alphas(new std::vector<double>(fPuppiContainer->puppiRawAlphas()));
    std::unique_ptr<double> nalgos(new double(fPuppiContainer->puppiNAlgos()));

    iEvent.put(std::move(alphas), "PuppiRawAlphas");
    iEvent.put(std::move(nalgos), "PuppiNAlgos");
    iEvent.put(std::move(theAlphas), "PuppiAlphas");
    iEvent.put(std::move(theAlphasMed), "PuppiAlphasMed");
    iEvent.put(std::move(theAlphasRms), "PuppiAlphasRms");
  }
}

// ------------------------------------------------------------------------------------------
void PuppiPVRobustProducer::beginJob() {}
// ------------------------------------------------------------------------------------------
void PuppiPVRobustProducer::endJob() {}
// ------------------------------------------------------------------------------------------
void PuppiPVRobustProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<bool>("puppiDiagnostics", false);
  desc.add<bool>("puppiForLeptons", false);
  desc.add<bool>("UseFromPVLooseTight", false);
  desc.add<bool>("UseDeltaZCut", true);
  desc.add<double>("DeltaZCut", 0.3);
  desc.add<double>("EtaMinUseDeltaZ", 0.);
  desc.add<double>("PtMaxCharged", 0.);
  desc.add<double>("EtaMaxCharged", 99999.);
  desc.add<double>("PtMaxNeutrals", 200.);
  desc.add<double>("PtMaxNeutralsStartSlope", 0.);
  desc.add<uint>("NumOfPUVtxsForCharged", 0);
  desc.add<double>("DeltaZCutForChargedFromPUVtxs", 0.2);
  desc.add<bool>("useExistingWeights", false);
  desc.add<bool>("useWeightsNoLep", false);
  desc.add<bool>("clonePackedCands", false);
  desc.add<int>("vtxNdofCut", 4);
  desc.add<double>("vtxZCut", 24);
  desc.add<edm::InputTag>("candName", edm::InputTag("particleFlow"));
  desc.add<edm::InputTag>("vertexName", edm::InputTag("offlinePrimaryVertices"));
  desc.add<edm::InputTag>("muonName", edm::InputTag("slimmedMuons"));
  desc.add<bool>("applyCHS", true);
  desc.add<bool>("invertPuppi", false);
  desc.add<bool>("useExp", false);
  desc.add<double>("MinPuppiWeight", .01);

  PuppiAlgo::fillDescriptionsPuppiAlgo(desc);

  descriptions.add("PuppiPVRobustProducer", desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(PuppiPVRobustProducer);
