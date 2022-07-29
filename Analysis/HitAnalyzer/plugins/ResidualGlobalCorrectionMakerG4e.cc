#include "ResidualGlobalCorrectionMakerBase.h"
#include "MagneticFieldOffset.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "TrackPropagation/Geant4e/interface/Geant4ePropagator.h"


#include "DataFormats/PatCandidates/interface/Muon.h"

#include "Math/Vector4Dfwd.h"
#include "Math/Vector4D.h"

#include <Eigen/Sparse>

#include "TRandom.h"

#include "Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h"



class ResidualGlobalCorrectionMakerG4e : public ResidualGlobalCorrectionMakerBase
{
public:
  explicit ResidualGlobalCorrectionMakerG4e(const edm::ParameterSet &);
  ~ResidualGlobalCorrectionMakerG4e() {}

//   static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  
  virtual void beginStream(edm::StreamID) override;
  virtual void produce(edm::Event &, const edm::EventSetup &) override;

  edm::EDGetTokenT<edm::Association<reco::TrackExtraCollection>> inputAssoc_;
  
  float muonPt;
  bool muonLoose;
  bool muonMedium;
  bool muonTight;
  bool muonIsPF;
  bool muonIsTracker;
  bool muonIsGlobal;
  bool muonIsStandalone;
  bool muonInnerTrackBest;
  
  bool trackExtraAssoc;
  
  std::vector<float> dEpred;
  std::vector<float> dE;
  std::vector<float> sigmadE;
  std::vector<float> Epred;
  std::vector<float> E;
  
  float outPt;
  float outEta;
  float outPhi;
  
  float outPtStart;
  float outEtaStart;
  float outPhiStart;

  edm::EDPutTokenT<edm::ValueMap<float>> outputCorPt_;
  edm::EDPutTokenT<edm::ValueMap<float>> outputCorEta_;
  edm::EDPutTokenT<edm::ValueMap<float>> outputCorPhi_;
  edm::EDPutTokenT<edm::ValueMap<int>> outputCorCharge_;
  edm::EDPutTokenT<edm::ValueMap<float>> outputEdmval_;

  edm::EDPutTokenT<edm::ValueMap<std::vector<int>>> outputGlobalIdxs_;

  edm::EDPutTokenT<edm::ValueMap<std::vector<float>>> outputJacRef_;
  edm::EDPutTokenT<edm::ValueMap<std::vector<float>>> outputMomCov_;

};


ResidualGlobalCorrectionMakerG4e::ResidualGlobalCorrectionMakerG4e(const edm::ParameterSet &iConfig) : ResidualGlobalCorrectionMakerBase(iConfig) 
{
  
//   inputAssoc_ = consumes<edm::Association<reco::TrackExtraCollection>>(edm::InputTag("muonReducedTrackExtras"));

  outputCorPt_ = produces<edm::ValueMap<float>>("corPt");
  outputCorEta_ = produces<edm::ValueMap<float>>("corEta");
  outputCorPhi_ = produces<edm::ValueMap<float>>("corPhi");
  outputCorCharge_ = produces<edm::ValueMap<int>>("corCharge");
  outputEdmval_ = produces<edm::ValueMap<float>>("edmval");

  outputGlobalIdxs_ = produces<edm::ValueMap<std::vector<int>>>("globalIdxs");

  outputJacRef_ = produces<edm::ValueMap<std::vector<float>>>("jacRef");
  outputMomCov_ = produces<edm::ValueMap<std::vector<float>>>("momCov");
  
}

void ResidualGlobalCorrectionMakerG4e::beginStream(edm::StreamID streamid)
{
  ResidualGlobalCorrectionMakerBase::beginStream(streamid);
  
  if (fillTrackTree_) {
    const int basketSize = 4*1024*1024;
    
    tree->Branch("trackPt", &trackPt, basketSize);
    tree->Branch("trackPtErr", &trackPtErr, basketSize);
    tree->Branch("trackEta", &trackEta, basketSize);
    tree->Branch("trackPhi", &trackPhi, basketSize);
    tree->Branch("trackCharge", &trackCharge, basketSize);
    //workaround for older ROOT version inability to store std::array automatically
  //   tree->Branch("trackOrigParms", trackOrigParms.data(), "trackOrigParms[5]/F", basketSize);
  //   tree->Branch("trackOrigCov", trackOrigCov.data(), "trackOrigCov[25]/F", basketSize);
    tree->Branch("trackParms", trackParms.data(), "trackParms[5]/F", basketSize);
    tree->Branch("trackCov", trackCov.data(), "trackCov[25]/F", basketSize);
    
    tree->Branch("refParms_iter0", refParms_iter0.data(), "refParms_iter0[5]/F", basketSize);
    tree->Branch("refCov_iter0", refCov_iter0.data(), "refCov_iter0[25]/F", basketSize);
  //   tree->Branch("refParms_iter2", refParms_iter2.data(), "refParms_iter2[5]/F", basketSize);
  //   tree->Branch("refCov_iter2", refCov_iter2.data(), "refCov_iter2[25]/F", basketSize);  
    
    tree->Branch("refParms", refParms.data(), "refParms[5]/F", basketSize);
    tree->Branch("refCov", refCov.data(), "refCov[25]/F", basketSize);
    tree->Branch("genParms", genParms.data(), "genParms[5]/F", basketSize);

    tree->Branch("genPt", &genPt, basketSize);
    tree->Branch("genEta", &genEta, basketSize);
    tree->Branch("genPhi", &genPhi, basketSize);
    tree->Branch("genCharge", &genCharge, basketSize);
    
    tree->Branch("genX", &genX, basketSize);
    tree->Branch("genY", &genY, basketSize);
    tree->Branch("genZ", &genZ, basketSize);
    
    tree->Branch("normalizedChi2", &normalizedChi2, basketSize);
    
    tree->Branch("nHits", &nHits, basketSize);
    tree->Branch("nValidHits", &nValidHits, basketSize);
    tree->Branch("nValidPixelHits", &nValidPixelHits, basketSize);
    
    tree->Branch("nValidHitsFinal", &nValidHitsFinal);
    tree->Branch("nValidPixelHitsFinal", &nValidPixelHitsFinal);

    if (fillJac_) {
      tree->Branch("nJacRef", &nJacRef, basketSize);
      tree->Branch("jacrefv",jacrefv.data(),"jacrefv[nJacRef]/F", basketSize);
    }
    
    tree->Branch("dEpred", &dEpred);
    tree->Branch("dE", &dE);
    tree->Branch("sigmadE", &sigmadE);
    tree->Branch("E", &E);
    tree->Branch("Epred", &Epred);
    
//     tree->Branch("outPt", &outPt);
//     tree->Branch("outEta", &outEta);
//     tree->Branch("outPhi", &outPhi);
//     
//     tree->Branch("outPtStart", &outPtStart);
//     tree->Branch("outEtaStart", &outEtaStart);
//     tree->Branch("outPhiStart", &outPhiStart);

    tree->Branch("muonPt", &muonPt);
    tree->Branch("muonLoose", &muonLoose);
    tree->Branch("muonMedium", &muonMedium);
    tree->Branch("muonTight", &muonTight);

    tree->Branch("muonIsPF", &muonIsPF);
    tree->Branch("muonIsTracker", &muonIsTracker);
    tree->Branch("muonIsGlobal", &muonIsGlobal);
    tree->Branch("muonIsStandalone", &muonIsStandalone);

    tree->Branch("muonInnerTrackBest", &muonInnerTrackBest);

    tree->Branch("trackExtraAssoc", &trackExtraAssoc);
    
    if (fitFromGenParms_) {

  //     tree->Branch("dxpxb1", &dxpxb1);
  //     tree->Branch("dypxb1", &dypxb1);
  //     
  //     tree->Branch("dxttec9rphi", &dxttec9rphi);
  //     tree->Branch("dxttec9stereo", &dxttec9stereo);
  //     
  //     tree->Branch("dxttec4rphi", &dxttec4rphi);
  //     tree->Branch("dxttec4stereo", &dxttec4stereo);
  //     
  //     tree->Branch("dxttec4rphisimgen", &dxttec4rphisimgen);
  //     tree->Branch("dyttec4rphisimgen", &dyttec4rphisimgen);
  //     tree->Branch("dxttec4rphirecsim", &dxttec4rphirecsim);
  //     
  //     tree->Branch("dxttec9rphisimgen", &dxttec9rphisimgen);
  //     tree->Branch("dyttec9rphisimgen", &dyttec9rphisimgen);
  //     
  //     tree->Branch("simlocalxref", &simlocalxref);
  //     tree->Branch("simlocalyref", &simlocalyref);
      
      tree->Branch("hitidxv", &hitidxv);
      tree->Branch("dxrecgen", &dxrecgen);
      tree->Branch("dyrecgen", &dyrecgen);
      tree->Branch("dxsimgen", &dxsimgen);
      tree->Branch("dysimgen", &dysimgen);
      tree->Branch("dxsimgenconv", &dxsimgenconv);
      tree->Branch("dysimgenconv", &dysimgenconv);
      tree->Branch("dxsimgenlocal", &dxsimgenlocal);
      tree->Branch("dysimgenlocal", &dysimgenlocal);
      tree->Branch("dxrecsim", &dxrecsim);
      tree->Branch("dyrecsim", &dyrecsim);
      tree->Branch("dxerr", &dxerr);
      tree->Branch("dyerr", &dyerr);
      
      tree->Branch("clusterSize", &clusterSize);
      tree->Branch("clusterSizeX", &clusterSizeX);
      tree->Branch("clusterSizeY", &clusterSizeY);
      tree->Branch("clusterCharge", &clusterCharge);
      tree->Branch("clusterChargeBin", &clusterChargeBin);
      tree->Branch("clusterOnEdge", &clusterOnEdge);
      
      tree->Branch("clusterProbXY", &clusterProbXY);
      tree->Branch("clusterSN", &clusterSN);

      tree->Branch("stripsToEdge", &stripsToEdge);
      
      tree->Branch("dxreccluster", &dxreccluster);
      tree->Branch("dyreccluster", &dyreccluster);

      tree->Branch("simlocalqop", &simlocalqop);
      tree->Branch("simlocaldxdz", &simlocaldxdz);
      tree->Branch("simlocaldydz", &simlocaldydz);
      tree->Branch("simlocalx", &simlocalx);
      tree->Branch("simlocaly", &simlocaly);
      
      tree->Branch("localqop", &localqop);
      tree->Branch("localdxdz", &localdxdz);
      tree->Branch("localdydz", &localdydz);
      tree->Branch("localx", &localx);
      tree->Branch("localy", &localy);

      tree->Branch("localqop_iter", &localqop_iter);
      tree->Branch("localdxdz_iter", &localdxdz_iter);
      tree->Branch("localdydz_iter", &localdydz_iter);
      tree->Branch("localx_iter", &localx_iter);
      tree->Branch("localy_iter", &localy_iter);
      
      tree->Branch("localqoperr", &localqoperr);
      
      tree->Branch("localphi", &localphi);
      tree->Branch("hitphi", &hitphi);

      
      tree->Branch("simtestz", &simtestz);
      tree->Branch("simtestvz", &simtestvz);
      tree->Branch("simtestrho", &simtestrho);
      tree->Branch("simtestzlocalref", &simtestzlocalref);
      tree->Branch("simtestdx", &simtestdx);
      tree->Branch("simtestdxrec", &simtestdxrec);
      tree->Branch("simtestdy", &simtestdy);
      tree->Branch("simtestdyrec", &simtestdyrec);
      tree->Branch("simtestdxprop", &simtestdxprop);
      tree->Branch("simtestdyprop", &simtestdyprop);
      tree->Branch("simtestdetid", &simtestdetid);
      
      tree->Branch("rx", &rx);
      tree->Branch("ry", &ry);
      
      tree->Branch("deigx", &deigx);
      tree->Branch("deigy", &deigy);
    }
    
    
    nJacRef = 0.;
  }
}


// ------------ method called for each event  ------------
void ResidualGlobalCorrectionMakerG4e::produce(edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  
  const bool dogen = fitFromGenParms_;

  constexpr bool dolocalupdate = false;

  using namespace edm;

  Handle<reco::TrackCollection> trackOrigH;
  iEvent.getByToken(inputTrackOrig_, trackOrigH);

  edm::ESHandle<GlobalTrackingGeometry> globalGeometry;
  iSetup.get<GlobalTrackingGeometryRecord>().get(globalGeometry);

  edm::ESHandle<TrackerTopology> trackerTopology;
  iSetup.get<TrackerTopologyRcd>().get(trackerTopology);
  
  edm::ESHandle<TransientTrackingRecHitBuilder> ttrh;
  iSetup.get<TransientRecHitRecord>().get("WithAngleAndTemplate",ttrh);
  
  ESHandle<Propagator> thePropagator;
  iSetup.get<TrackingComponentsRecord>().get("Geant4ePropagator", thePropagator);
  
  const Geant4ePropagator *g4prop = dynamic_cast<const Geant4ePropagator*>(thePropagator.product());
  const MagneticField* field = thePropagator->magneticField();
  
  ESHandle<MagneticField> fieldh;
  iSetup.get<IdealMagneticFieldRecord>().get("", fieldh);
  std::unique_ptr<MagneticFieldOffset> fieldOffset = std::make_unique<MagneticFieldOffset>(&(*fieldh));
  
  std::unique_ptr<PropagatorWithMaterial> fPropagator = std::make_unique<PropagatorWithMaterial>(alongMomentum, 0.105, fieldOffset.get(), 1.6, true,  -1., true);
  
  constexpr double mmu = 0.1056583745;

  Handle<reco::BeamSpot> bsH;
  iEvent.getByToken(inputBs_, bsH);

  
  Handle<edm::View<reco::Candidate>> genPartCollection;
  Handle<math::XYZPointF> genXyz0;
  Handle<GenEventInfoProduct> genEventInfo;
  Handle<std::vector<int>> genPartBarcodes;
  Handle<std::vector<PileupSummaryInfo>> pileupSummary;
  if (doGen_) {
    iEvent.getByToken(GenParticlesToken_, genPartCollection);
    iEvent.getByToken(genXyz0Token_, genXyz0);
    iEvent.getByToken(genEventInfoToken_, genEventInfo);
    iEvent.getByToken(pileupSummaryToken_, pileupSummary);
  }
  
  std::vector<Handle<std::vector<PSimHit>>> simHits(inputSimHits_.size());
  edm::Handle<std::vector<SimTrack>> simTracks;
  if (doSim_) {
    iEvent.getByToken(genParticlesBarcodeToken_, genPartBarcodes);
    for (unsigned int isimhit = 0; isimhit<inputSimHits_.size(); ++isimhit) {
      iEvent.getByToken(inputSimHits_[isimhit], simHits[isimhit]);
    }
    iEvent.getByToken(inputSimTracks_, simTracks);
  }
  
  Handle<edm::View<reco::Muon> > muons;
  if (doMuons_) {
    iEvent.getByToken(inputMuons_, muons);
  }
  
  Handle<edm::Association<std::vector<pat::Muon>>> muonAssoc;
  if (doMuonAssoc_) {
    iEvent.getByToken(inputMuonAssoc_, muonAssoc);
  }
  
  TkClonerImpl const& cloner = static_cast<TkTransientTrackingRecHitBuilder const *>(ttrh.product())->cloner();

  run = iEvent.run();
  lumi = iEvent.luminosityBlock();
  event = iEvent.id().event();

  genweight = 1.;
  if (doGen_) {
    genweight = genEventInfo->weight();

    Pileup_nPU = pileupSummary->front().getPU_NumInteractions();
    Pileup_nTrueInt = pileupSummary->front().getTrueNumInteractions();
  }

  std::vector<float> corPtV;
  std::vector<float> corEtaV;
  std::vector<float> corPhiV;
  std::vector<int> corChargeV;
  std::vector<float> edmvalV;

  std::vector<std::vector<int>> globalidxsV;
  std::vector<std::vector<float>> jacRefV;
  std::vector<std::vector<float>> momCovV;

  std::array<double, 3> refParmsMomD;

  if (doMuonAssoc_) {
    corPtV.assign(muonAssoc->ref()->size(), -99.);
    corEtaV.assign(muonAssoc->ref()->size(), -99.);
    corPhiV.assign(muonAssoc->ref()->size(), -99.);
    corChargeV.assign(muonAssoc->ref()->size(), -99);
    edmvalV.assign(muonAssoc->ref()->size(), -99.);
    globalidxsV.assign(muonAssoc->ref()->size(), std::vector<int>());
    jacRefV.assign(muonAssoc->ref()->size(), std::vector<float>());
    momCovV.assign(muonAssoc->ref()->size(), std::vector<float>());
  }

//   for (const reco::Track &track : *trackOrigH) {
  for (unsigned int itrack = 0; itrack < trackOrigH->size(); ++itrack) {
    const reco::Track &track = (*trackOrigH)[itrack];
    const reco::TrackRef trackref(trackOrigH, itrack);

    const edm::Ref<std::vector<pat::Muon>> muonref = doMuonAssoc_ ? (*muonAssoc)[trackref] : edm::Ref<std::vector<pat::Muon>>();

    const bool iscosmic = track.algo() == reco::TrackBase::ctf || track.algo() == reco::TrackBase::cosmics;
    
    const bool dopca = !dogen && !iscosmic;
    
    if (track.isLooper()) {
      continue;
    }
    
    trackPt = track.pt();
    trackEta = track.eta();
    trackPhi = track.phi();
    trackCharge = track.charge();
    trackPtErr = track.ptError();
    
    normalizedChi2 = track.normalizedChi2();
    
    auto const& tkparms = track.parameters();
    auto const& tkcov = track.covariance();
    trackParms.fill(0.);
    trackCov.fill(0.);
    //use eigen to fill raw memory
    Map<Vector5f>(trackParms.data()) = Map<const Vector5d>(tkparms.Array()).cast<float>();
    Map<Matrix<float, 5, 5, RowMajor> >(trackCov.data()).triangularView<Upper>() = Map<const Matrix<double, 5, 5, RowMajor> >(tkcov.Array()).cast<float>().triangularView<Upper>();
    
//     std::cout << "track charge: " << track.charge() << " trackorig charge " << trackOrig.charge() << "inner state charge " << tms.back().updatedState().charge() << std::endl;
    
    const reco::Candidate* genpart = nullptr;
    
    genPt = -99.;
    genEta = -99.;
    genPhi = -99.;
    genCharge = -99;
    genX = -99.;
    genY = -99.;
    genZ = -99.;
    genParms.fill(0.);
    genl3d = -99.;
    
    int genBarcode = -99;
    
    
    if (doGen_) {
      
      float drmin = 0.1;
      
      for (auto g = genPartCollection->begin(); g != genPartCollection->end(); ++g)
      {
        if (g->status() != 1) {
          continue;
        }
        if (std::abs(g->pdgId()) != 13) {
          continue;
        }
        
        float dR = deltaR(*g, track);
        
        if (dR < drmin)
        {
          drmin = dR;
          
          genpart = &(*g);
          
          if (doSim_) {
            genBarcode = (*genPartBarcodes)[g - genPartCollection->begin()];
          }
          
          genPt = g->pt();
          genEta = g->eta();
          genPhi = g->phi();
          genCharge = g->charge();
          
          genX = g->vertex().x();
          genY = g->vertex().y();
          genZ = g->vertex().z();

          genl3d = std::sqrt((g->vertex() - *genXyz0).mag2());

          auto const& vtx = g->vertex();
          auto const& myBeamSpot = bsH->position(vtx.z());
          
          //q/|p|
          genParms[0] = g->charge()/g->p();
          //lambda
          genParms[1] = M_PI_2 - g->momentum().theta();
          //phi
          genParms[2] = g->phi();
          //dxy
          genParms[3] = (-(vtx.x() - myBeamSpot.x()) * g->py() + (vtx.y() - myBeamSpot.y()) * g->px()) / g->pt();
          //dsz
          genParms[4] = (vtx.z() - myBeamSpot.z()) * g->pt() / g->p() -
            ((vtx.x() - myBeamSpot.x()) * g->px() + (vtx.y() - myBeamSpot.y()) * g->py()) / g->pt() * g->pz() / g->p();
        }
        else {
          continue;
        }
      }
    }
    
//     std::cout << "genPt = " << genPt << " genEta = " << genEta << " genPhi = " << genPhi << " genCharge = " << genCharge << " genX = " << genX << " genY = " << genY << " genZ = " << genZ << std::endl;
    
    int simtrackid = -99;
    if (genpart != nullptr && doSim_) {
      for (auto const& simTrack : *simTracks) {
        if (simTrack.genpartIndex() == genBarcode) {
          simtrackid = simTrack.trackId();
          break;
        }
      }
    }
    
    
    muonPt = -99.;
    muonLoose = false;
    muonMedium = false;
    muonTight = false;
    muonIsTracker = false;
    muonIsGlobal = false;
    muonIsStandalone = false;
    muonInnerTrackBest = false;
    trackExtraAssoc = false;

    const reco::Muon *matchedmuon = nullptr;

    if (doMuons_) {
      for (auto const &muon : *muons) {
        if (muon.bestTrack()->algo() == track.algo()) {
          if ( (muon.bestTrack()->momentum() - track.momentum()).mag2() < 1e-3 ) {
            matchedmuon = &muon;
          }
        }
        else if (muon.innerTrack().isNonnull() && muon.innerTrack()->algo() == track.algo()) {
          if ( (muon.innerTrack()->momentum() - track.momentum()).mag2() < 1e-3 ) {
            matchedmuon = &muon;
          }
        }
      }
    }

    if (matchedmuon != nullptr) {
      muonPt = matchedmuon->pt();
      muonLoose = matchedmuon->passed(reco::Muon::CutBasedIdLoose);
      muonMedium = matchedmuon->passed(reco::Muon::CutBasedIdMedium);
      muonTight = matchedmuon->passed(reco::Muon::CutBasedIdTight);
      muonIsPF = matchedmuon->isPFMuon();
      muonIsTracker = matchedmuon->isTrackerMuon();
      muonIsGlobal = matchedmuon->isGlobalMuon();
      muonIsStandalone = matchedmuon->isStandAloneMuon();
      muonInnerTrackBest = matchedmuon->muonBestTrackType() == reco::Muon::InnerTrack;
    }
    
    //prepare hits
    TransientTrackingRecHit::RecHitContainer hits;
    hits.reserve(track.recHitsSize());

//     std::cout << "track: algo = " << track.algo() << " originalAlgo = " << track.originalAlgo() << std::endl;
//
//
//     if (track.seedDirection() == oppositeToMomentum) {
// //       std::cout << "track with oppositeToMomentum: algo = " << track.algo() << " originalAlgo = " << track.originalAlgo() << std::endl;
//       std::cout << "track with oppositeToMomentum:" << std::endl;
//     }


    
    for (auto it = track.recHitsBegin(); it != track.recHitsEnd(); ++it) {
//       if (track.seedDirection() == oppositeToMomentum) {
//         std::cout << "det = " << (*it)->geographicalId().det() << std::endl;
//       }


      if ((*it)->geographicalId().det() != DetId::Tracker) {
        continue;
      }

      const GeomDet* detectorG = globalGeometry->idToDet((*it)->geographicalId());
      const GluedGeomDet* detglued = dynamic_cast<const GluedGeomDet*>(detectorG);
      
//       if (track.seedDirection() == oppositeToMomentum) {
//         std::cout << "position mag = " << detectorG->surface().position().mag() << std::endl;
//       }

      // split matched invalid hits
      if (detglued != nullptr && !(*it)->isValid()) {
//         bool order = detglued->stereoDet()->surface().position().mag() > detglued->monoDet()->surface().position().mag();
        
        const auto stereopos = detglued->stereoDet()->surface().position();
        const auto monopos = detglued->monoDet()->surface().position();
        
        const Eigen::Vector3d stereoposv(stereopos.x(), stereopos.y(), stereopos.z());
        const Eigen::Vector3d monoposv(monopos.x(), monopos.y(), monopos.z());
        const Eigen::Vector3d trackmomv(track.momentum().x(), track.momentum().y(), track.momentum().z());
        
        bool order = (stereoposv - monoposv).dot(trackmomv) > 0.;
        

        const GeomDetUnit* detinner = order ? detglued->monoDet() : detglued->stereoDet();
        const GeomDetUnit* detouter = order ? detglued->stereoDet() : detglued->monoDet();
        
        hits.push_back(TrackingRecHit::RecHitPointer(new InvalidTrackingRecHit(*detinner, (*it)->type())));
        hits.push_back(TrackingRecHit::RecHitPointer(new InvalidTrackingRecHit(*detouter, (*it)->type())));
        
      }
      else {
        // apply hit quality criteria
        const bool ispixel = GeomDetEnumerators::isTrackerPixel(detectorG->subDetector());
//         bool hitquality = true;
        bool hitquality = false;
        if (applyHitQuality_ && (*it)->isValid()) {
          const TrackerSingleRecHit* tkhit = dynamic_cast<const TrackerSingleRecHit*>(*it);
          assert(tkhit != nullptr);
          
          if (ispixel) {
            const SiPixelRecHit *pixhit = dynamic_cast<const SiPixelRecHit*>(tkhit);
            const SiPixelCluster& cluster = *tkhit->cluster_pixel();
            assert(pixhit != nullptr);
            
            hitquality = !pixhit->isOnEdge() && cluster.sizeX() > 1;
          }
          else {
            assert(tkhit->cluster_strip().isNonnull());
            const SiStripCluster& cluster = *tkhit->cluster_strip();
            const StripTopology* striptopology = dynamic_cast<const StripTopology*>(&(detectorG->topology()));
            assert(striptopology);
            
            const uint16_t firstStrip = cluster.firstStrip();
            const uint16_t lastStrip = cluster.firstStrip() + cluster.amplitudes().size() - 1;
            const bool isOnEdge = firstStrip == 0 || lastStrip == (striptopology->nstrips() - 1);
            
            hitquality = true;
          }
          
        }
        else {
          hitquality = true;
        }
        
        if (hitquality) {
          hits.push_back((*it)->cloneForFit(*detectorG));
        }
        else {
          hits.push_back(TrackingRecHit::RecHitPointer(new InvalidTrackingRecHit(*detectorG, TrackingRecHit::inactive)));
        }
      }
    }

    const unsigned int nhits = hits.size();

    if (nhits == 0) {
      continue;
    }

    nHits = nhits;

    unsigned int nvalid = 0;
    unsigned int nvalidpixel = 0;
    unsigned int nvalidalign2d = 0;
    
    
    // count valid hits since this is needed to size the arrays
    for (auto const& hit : hits) {
      
      assert(hit->dimension()<=2);
      if (hit->isValid()) {
        nvalid += 1;
        
        const uint32_t gluedid = trackerTopology->glued(hit->geographicalId());
        const bool isglued = gluedid != 0;
        const DetId parmdetid = isglued ? DetId(gluedid) : hit->geographicalId();
        
//         const bool align2d = detidparms.count(std::make_pair(1, hit->geographicalId()));
        const bool align2d = detidparms.count(std::make_pair(1, parmdetid));
//         
        if (align2d) {
          nvalidalign2d += 1;
        }
        if (GeomDetEnumerators::isTrackerPixel(hit->det()->subDetector())) {
          nvalidpixel += 1;
        }
      }
    }

    if (nvalid == 0) {
      continue;
    }
    
    nValidHits = nvalid;
    nValidPixelHits = nvalidpixel;
    
    nValidHitsFinal = 0;
    nValidPixelHitsFinal = 0;
    
    const unsigned int nparsAlignment = 5*nvalid + nvalidalign2d;
    const unsigned int nparsBfield = nhits;
    const unsigned int nparsEloss = nhits;
    const unsigned int npars = nparsAlignment + nparsBfield + nparsEloss;
    
    const unsigned int nstateparms = 5*(nhits+1);
    const unsigned int nparmsfull = nstateparms + npars;
    
    VectorXd gradfull;
    MatrixXd hessfull;
    
    std::vector<MatrixXd> dhessv;

    
    MatrixXd validdxeigjac;
    VectorXd validdxeig;
    
    Matrix<float, Dynamic, 2> rxfull(nvalid, 2); 
    Matrix<float, Dynamic, 2> ryfull(nvalid, 2);
    
    globalidxv.clear();
    globalidxv.resize(npars, 0);
    
    VectorXd dxfull;
    MatrixXd dxdparms;
    VectorXd grad;
    MatrixXd hess;
    LDLT<MatrixXd> Cinvd;
    MatrixXd covfull = MatrixXd::Zero(nstateparms, nstateparms);
    
    if (dogen && genpart==nullptr) {
      std::cout << "no gen part, skipping track\n";
      continue;
    }
    

    if (debugprintout_) {
      std::cout << "initial reference point parameters:" << std::endl;
      std::cout << track.parameters() << std::endl;
    }

    Matrix<double, 7, 1> refFts;
    
    if (dogen) {
      if (genpart==nullptr) {
        continue;
      }
      //init from gen state
      auto const& refpoint = genpart->vertex();
      auto const& trackmom = genpart->momentum();
      
      refFts[0] = refpoint.x();
      refFts[1] = refpoint.y();
      refFts[2] = refpoint.z();
      
      refFts[3] = trackmom.x();
      refFts[4] = trackmom.y();
      refFts[5] = trackmom.z();
      
      refFts[6] = genpart->charge();
    }
    else {
      //init from track state
      auto const& refpoint = track.referencePoint();
      auto const& trackmom = track.momentum();
      
      refFts[0] = refpoint.x();
      refFts[1] = refpoint.y();
      refFts[2] = refpoint.z();
      
      refFts[3] = trackmom.x();
      refFts[4] = trackmom.y();
      refFts[5] = trackmom.z();
      
      refFts[6] = track.charge();
      
      // special case for cosmics, start from "inner" state with straight line extrapolation back 1cm to preserve the propagation logic
      // (the reference point is instead the PCA to the beamline and is therefore in the middle of the trajectory and not compatible
      // with the fitter/propagation logic
      // TODO change this so that cosmics directly use a local parameterization on the first hit?
      if (iscosmic) {
        auto const& startpoint = track.extra()->innerPosition();
        auto const& startmom = track.extra()->innerMomentum();
        
        const Eigen::Vector3d startpointv(startpoint.x(), startpoint.y(), startpoint.z());
        const Eigen::Vector3d startmomv(startmom.x(), startmom.y(), startmom.z());
        
        refFts.head<3>() = startpointv - 1.0*startmomv.normalized();
        refFts.segment<3>(3) = startmomv;

      }
    }

    if (dopca) {
      // enforce that reference state is really a consistent PCA to the beamline (by adjusting the reference position if needed)
      const Matrix<double, 5, 1> statepca = cart2pca(refFts, *bsH);
      refFts = pca2cart(statepca, *bsH);
    }

    std::vector<Matrix<double, 7, 1>> layerStates;
    layerStates.reserve(nhits);
    
    bool valid = true;

    const bool islikelihood = false;

    

    std::vector<double> localxsmearedsim(hits.size(), 0.);
    
    
    double chisqvalold = std::numeric_limits<double>::max();
    
    const bool anomDebug = false;
    
    const unsigned int niters = (dogen && !dolocalupdate) || (dogen && fitFromSimParms_) ? 1 : 10;
    
    for (unsigned int iiter=0; iiter<niters; ++iiter) {
      if (debugprintout_) {
        std::cout<< "iter " << iiter << std::endl;
      }

//       std::cout<< "iter " << iiter << std::endl;
            
      hitidxv.clear();
      hitidxv.reserve(nvalid);
      
      if (iiter == 0) {
        dxrecgen.clear();
        dxrecgen.reserve(nvalid);
        
        dyrecgen.clear();
        dyrecgen.reserve(nvalid);
        
        dxsimgen.clear();
        dxsimgen.reserve(nvalid);
        
        dysimgen.clear();
        dysimgen.reserve(nvalid);

        dxsimgenconv.clear();
        dxsimgenconv.reserve(nvalid);

        dysimgenconv.clear();
        dysimgenconv.reserve(nvalid);
        
        dxsimgenlocal.clear();
        dxsimgenlocal.reserve(nvalid);
        
        dysimgenlocal.clear();
        dysimgenlocal.reserve(nvalid);
        
        dxrecsim.clear();
        dxrecsim.reserve(nvalid);
        
        dyrecsim.clear();
        dyrecsim.reserve(nvalid);
        
        dxreccluster.clear();
        dxreccluster.reserve(nvalid);
        
        dyreccluster.clear();
        dyreccluster.reserve(nvalid);
        
        dxerr.clear();
        dxerr.reserve(nvalid);
        
        dyerr.clear();
        dyerr.reserve(nvalid);
        
        clusterSize.clear();
        clusterSize.reserve(nvalid);
        
        clusterSizeX.clear();
        clusterSizeX.reserve(nvalid);
        
        clusterSizeY.clear();
        clusterSizeY.reserve(nvalid);
        
        clusterCharge.clear();
        clusterCharge.reserve(nvalid);
        
        clusterChargeBin.clear();
        clusterChargeBin.reserve(nvalid);
        
        clusterOnEdge.clear();
        clusterOnEdge.reserve(nvalid);
        
        clusterProbXY.clear();
        clusterProbXY.reserve(nvalid);
        
        clusterSN.clear();
        clusterSN.reserve(nvalid);
        
        simlocalqop.clear();
        simlocaldxdz.clear();
        simlocaldydz.clear();
        simlocalx.clear();
        simlocaly.clear();

        simlocalqop.reserve(nvalid);
        simlocaldxdz.reserve(nvalid);
        simlocaldydz.reserve(nvalid);
        simlocalx.reserve(nvalid);
        simlocaly.reserve(nvalid);
        
        localqop.clear();
        localdxdz.clear();
        localdydz.clear();
        localx.clear();
        localy.clear();

        localqop.reserve(nvalid);
        localdxdz.reserve(nvalid);
        localdydz.reserve(nvalid);
        localx.reserve(nvalid);
        localy.reserve(nvalid);
        
        localqop_iter.assign(nvalid, -99.);
        localdxdz_iter.assign(nvalid, -99.);
        localdydz_iter.assign(nvalid, -99.);
        localx_iter.assign(nvalid, -99.);
        localy_iter.assign(nvalid, -99.);
        
        localqoperr.assign(nvalid, -99.);

        localphi.clear();
        localphi.reserve(nvalid);

        hitphi.clear();
        hitphi.reserve(nvalid);
        
        dEpred.clear();
        dEpred.reserve(nvalid);
        
        dE.clear();
        dE.reserve(nvalid);
        
        sigmadE.clear();
        sigmadE.reserve(nvalid);
        
        E.clear();
        E.reserve(nvalid);
        
        Epred.clear();
        Epred.reserve(nvalid);

        stripsToEdge.clear();
        stripsToEdge.reserve(nvalid);

      }      
      
        
//       const bool islikelihood = true;
//       const bool islikelihood = iiter > 2;
//       const bool islikelihood = false;
      
//       const bool islikelihood = iiter > 0;
//       const bool islikelihood = true;
      
      gradfull = VectorXd::Zero(nparmsfull);
      hessfull = MatrixXd::Zero(nparmsfull, nparmsfull);
      
      if (islikelihood) {
        dhessv.assign(nhits, MatrixXd::Zero(nstateparms, nstateparms));
      }
      
      
      validdxeigjac = MatrixXd::Zero(2*nvalid, nstateparms);
      
      
      double chisq0val = 0.;
      
      dxpxb1 = -99.;
      dypxb1 = -99.;
      dxttec9rphi = -99.;
      dxttec9stereo = -99.;
      dxttec4rphi = -99.;
      dxttec4stereo = -99.;
      
      dxttec4rphisimgen = -99.;
      dyttec4rphisimgen = -99.;
      dxttec4rphirecsim = -99.;
      
      dxttec9rphisimgen = -99.;
      dyttec9rphisimgen = -99.;
      
      simlocalxref = -99.;
      simlocalyref = -99.;
      
      
      
      unsigned int parmidx = 0;
      unsigned int alignmentparmidx = 0;
      unsigned int ivalidhit = 0;

      if (iiter > 0) {
        //update current state from reference point state (errors not needed beyond first iteration)
//         JacobianCurvilinearToCartesian curv2cart(refFts.parameters());
//         const AlgebraicMatrix65& jac = curv2cart.jacobian();
//         const AlgebraicVector6 glob = refFts.parameters().vector();
        
        auto const& dxlocal = dxfull.head<5>();
//         const Matrix<double, 6, 1> globupd = Map<const Matrix<double, 6, 1>>(glob.Array()) + Map<const Matrix<double, 6, 5, RowMajor>>(jac.Array())*dxlocal;
//         const Matrix<double, 6, 1> globupd = Map<const Matrix<double, 6, 1>>(glob.Array()) + jac*dxlocal;

        
        

        
        if (dopca) {
          // position transformation from track-beamline pca to cartesian is non-linear so do it exactly
          const Matrix<double, 5, 1> statepca = cart2pca(refFts, *bsH);
          const Matrix<double, 5, 1> statepcaupd = statepca + dxlocal;

          refFts = pca2cart(statepcaupd, *bsH);

        }
        else {
          // position transformation from curvilinear to cartesian is linear so just use the jacobian

          const Matrix<double, 6, 5> jac = curv2cartJacobianAltD(refFts);
          const Matrix<double, 6, 1> globupd = refFts.head<6>() + jac*dxlocal;

          const double qbp = refFts[6]/refFts.segment<3>(3).norm();
          const double lam = std::atan(refFts[5]/std::sqrt(refFts[3]*refFts[3] + refFts[4]*refFts[4]));
          const double phi = std::atan2(refFts[4], refFts[3]);

          const double qbpupd = qbp + dxlocal(0);
          const double lamupd = lam + dxlocal(1);
          const double phiupd = phi + dxlocal(2);

          const double charge = std::copysign(1., qbpupd);
          const double pupd = std::abs(1./qbpupd);

          const double pxupd = pupd*std::cos(lamupd)*std::cos(phiupd);
          const double pyupd = pupd*std::cos(lamupd)*std::sin(phiupd);
          const double pzupd = pupd*std::sin(lamupd);

          refFts.head<3>() = globupd.head<3>();
          refFts[3] = pxupd;
          refFts[4] = pyupd;
          refFts[5] = pzupd;
          refFts[6] = charge;
        }
      }
      
      Matrix<double, 5, 5> ref2curvjac = dopca ? pca2curvJacobianD(refFts, field, *bsH) : Matrix<double, 5, 5>::Identity();

      Matrix<double, 7, 1> updtsos = refFts;

      Matrix<double, 5, 5> Qtot = Matrix<double, 5, 5>::Zero();
      
      float e = genpart == nullptr ? -99. : std::sqrt(genpart->momentum().mag2() + mmu*mmu);
      float epred = std::sqrt(refFts.segment<3>(3).squaredNorm() + mmu*mmu);
      
      
      if (bsConstraint_) {
        // apply beamspot constraint
        // TODO add residual corrections for beamspot parameters?
        
        constexpr unsigned int nlocalstate = 5;
        
        constexpr unsigned int nlocal = nlocalstate;
        
        constexpr unsigned int localstateidx = 0;
        
        constexpr unsigned int fullstateidx = 0;

        using BSScalar = AANT<double, nlocal>;
        
//         JacobianCurvilinearToCartesian curv2cart(refFts.parameters());
//         const AlgebraicMatrix65& jac = curv2cart.jacobian();
        const Matrix<double, 6, 5> jac = dopca ? pca2cartJacobianD(refFts, *bsH) : curv2cartJacobianAltD(refFts);
        
        const double sigb1 = bsH->BeamWidthX();
        const double sigb2 = bsH->BeamWidthY();
        const double sigb3 = bsH->sigmaZ();
        const double dxdz = bsH->dxdz();
        const double dydz = bsH->dydz();
        const double x0 = bsH->x0();
        const double y0 = bsH->y0();
        const double z0 = bsH->z0();
        
        
        // covariance matrix of luminous region in global coordinates
        // taken from https://github.com/cms-sw/cmssw/blob/abc1f17b230effd629c9565fb5d95e527abcb294/RecoVertex/BeamSpotProducer/src/FcnBeamSpotFitPV.cc#L63-L90

        // FIXME xy correlation is not stored and assumed to be zero
        const double corrb12 = 0.;
        
        const double varb1 = sigb1*sigb1;
        const double varb2 = sigb2*sigb2;
        const double varb3 = sigb3*sigb3;
        
        Matrix<double, 3, 3> covBS = Matrix<double, 3, 3>::Zero();
        // parametrisation: rotation (dx/dz, dy/dz); covxy
        covBS(0,0) = varb1;
        covBS(1,0) = covBS(0,1) = corrb12*sigb1*sigb2;
        covBS(1,1) = varb2;
        covBS(2,0) = covBS(0,2) = dxdz*(varb3-varb1)-dydz*covBS(1,0);
        covBS(2,1) = covBS(1,2) = dydz*(varb3-varb2)-dxdz*covBS(1,0);
        covBS(2,2) = varb3;

//         std::cout << "covBS:" << std::endl;
//         std::cout << covBS << std::endl;
        
        Matrix<BSScalar, 5, 1> du = Matrix<BSScalar, 5, 1>::Zero();
        for (unsigned int j=0; j<du.size(); ++j) {
          init_twice_active_var(du[j], nlocal, localstateidx + j);
        }
        
        Matrix<BSScalar, 3, 1> dbs0;
        dbs0[0] = BSScalar(refFts[0] - x0);
        dbs0[1] = BSScalar(refFts[1] - y0);
        dbs0[2] = BSScalar(refFts[2] - z0);
        
//         std::cout << "dposition / d(qop, lambda, phi) (should be 0?):" << std::endl;
//         std::cout << Map<const Matrix<double, 6, 5, RowMajor>>(jac.Array()).topLeftCorner<3,3>() << std::endl;
        
//         const Matrix<BSScalar, 3, 2> jacpos = Map<const Matrix<double, 6, 5, RowMajor>>(jac.Array()).topRightCorner<3,2>().cast<BSScalar>();
        const Matrix<BSScalar, 3, 5> jacpos = jac.topRows<3>().cast<BSScalar>();
        const Matrix<BSScalar, 3, 3> covBSinv = covBS.inverse().cast<BSScalar>();
        
        const Matrix<BSScalar, 3, 1> dbs = dbs0 + jacpos*du;
        const BSScalar chisq = dbs.transpose()*covBSinv*dbs;

        chisq0val += chisq.value().value();
        
        auto const& gradlocal = chisq.value().derivatives();
        //fill local hessian
        Matrix<double, nlocal, nlocal> hesslocal;
        for (unsigned int j=0; j<nlocal; ++j) {
          hesslocal.row(j) = chisq.derivatives()[j].derivatives();
        }
        
        //fill global gradient
        gradfull.segment<nlocalstate>(fullstateidx) += gradlocal.head<nlocalstate>();

        //fill global hessian (upper triangular blocks only)
        hessfull.block<nlocalstate,nlocalstate>(fullstateidx, fullstateidx) += hesslocal.topLeftCorner<nlocalstate,nlocalstate>();        
        
      }
      
      
      for (unsigned int ihit = 0; ihit < hits.size(); ++ihit) {
//         std::cout << "ihit " << ihit << std::endl;
        auto const& hit = hits[ihit];
        
        // check propagation direction to choose correct bfield and material parameters
        // when propagating inside-out the parameters correspond to the target module
        // when propagating outside-in (e.g. for cosmics) they correspond to the source module
        // For the first hit always use the target module

        bool sourceParms = false;
//         if (ihit > 0) {
        if (false) {
          auto const& lasthit = hits[ihit - 1];
          const GloballyPositioned<double> &lastsurface = surfacemapD_.at(lasthit->geographicalId());

          const Vector3DBase<double, LocalTag> localmomz(0., 0., 1.);
          const Vector3DBase<double, GlobalTag> globalmomz = lastsurface.toGlobal(localmomz);

          const Vector3DBase<double, GlobalTag> globalmomtmp(updtsos[3], updtsos[4], updtsos[5]);
          const Vector3DBase<double, LocalTag> localmomtmp = lastsurface.toLocal(globalmomtmp);

          // is local z dir facing out or in
          const double localzoutin = std::copysign(1.0, lastsurface.position().x()*globalmomz.x() + lastsurface.position().y()*globalmomz.y() + lastsurface.position().z()*globalmomz.z());
          const double localzdir = std::copysign(1.0, localmomtmp.z());

          sourceParms = localzoutin*localzdir < 0.;
        }

//         std::cout << "ihit = " << ihit << " sourceParms = " << sourceParms << std::endl;

        auto const& prophit = sourceParms ? hits[ihit - 1] : hit;
        const uint32_t gluedid = trackerTopology->glued(prophit->det()->geographicalId());
        const bool isglued = gluedid != 0;
        const DetId parmdetid = isglued ? DetId(gluedid) : prophit->geographicalId();
        const GeomDet* parmDet = isglued ? globalGeometry->idToDet(parmdetid) : prophit->det();
        const double xifraction = isglued ? prophit->det()->surface().mediumProperties().xi()/parmDet->surface().mediumProperties().xi() : 1.;
        
        const unsigned int bfieldglobalidx = detidparms.at(std::make_pair(6, parmdetid));                
        const unsigned int elossglobalidx = detidparms.at(std::make_pair(7, parmdetid));
        
        const double dbetaval = corparms_[bfieldglobalidx];
        const double dxival = corparms_[elossglobalidx];
        
        const PSimHit *simhit = nullptr;
        
        if (doSim_) {
          for (auto const& simhith : simHits) {
            for (const PSimHit& simHit : *simhith) {
              if (simHit.detUnitId() == hit->geographicalId() && int(simHit.trackId()) == simtrackid && std::abs(simHit.particleType()) == 13) {
                simhit = &simHit;
                break;
              }
            }
            if (simhit != nullptr) {
              break;
            }
          }
        }
        
                

//         const DetId partnerid = isglued ? trackerTopology->partnerDetId(hit->det()->geographicalId()) : DetId();
//         
//         const bool isfront = ihit != (hits.size() - 1) && isglued && hits[ihit+1]->det()->geographicalId() == partnerid;
//         const bool isback = ihit !=0 && isglued && hits[ihit-1]->det()->geographicalId() == partnerid;
//         
//         const double xifraction = isfront ? 0. : 1.;
                 
      
        const ROOT::Math::PxPyPzMVector momtmp(updtsos[3], updtsos[4], updtsos[5], mmu);
      
//         if (std::abs(momtmp.eta()) > 4.0) {
//           std::cout << "WARNING:  Invalid state!!!" << std::endl;
//           valid = false;
//           break;
//         }
        
        
        
//         const Matrix<double, 5, 5> curvcov = covfull.block<5, 5>(5*ihit, 5*ihit);


//         const double sigmaqop = std::sqrt(std::max(curvcov(0,0), 0.));
//
//         const double oldqop = updtsos[6]/updtsos.segment<3>(3).norm();
//
//         const Matrix<double, 6, 5> curv2cartforconv = curv2cartJacobianAltD(updtsos);
//
//         Matrix<double, 5, 1> curvvar = Matrix<double, 5, 1>::Zero();
//         curvvar[0] = 0.5*sigmaqop;
//
//         auto tsosqopup = updtsos;
//         tsosqopup.head<6>() += curv2cartforconv*curvvar;
//
//         auto tsosqopdown = updtsos;
//         tsosqopdown.head<6>() += -curv2cartforconv*curvvar;

        
//        auto const &surfaceip1 = *hits[ihit+1]->surface();
//           auto const &surface = *hit->surface();
//           const Plane &surface = *hit->surface();
//           const Surface &surface = *hit->surface();
        const GloballyPositioned<double> &surface = surfacemapD_.at(hit->geographicalId());
//         auto const &surface = *surfacemap_.at(hit->geographicalId());
//           auto propresult = thePropagator->propagateWithPath(updtsos, surface);
//           auto propresult = g4prop->propagateGenericWithJacobian(*updtsos.freeState(), surface);
        
//         std::cout << "propagation for iiter = " << iiter << " ihit = " << ihit << std::endl;
        
//         const Point3DBase<double, GlobalTag> crosspostmp(updtsos[0], updtsos[1], updtsos[2]);
//         const Vector3DBase<double, GlobalTag> crossmomtmp(updtsos[3], updtsos[4], updtsos[5]);
//         
//         if (surface.toLocal(crosspostmp).z() * surface.toLocal(crossmomtmp).z() > 0) {
//           std::cout << "Abort: wrong propagation direction!\n";
//           valid = false;
//           break;
//         }
        
        auto propresult = g4prop->propagateGenericWithJacobianAltD(updtsos, surface, dbetaval, dxival);

//           propresult = fPropagator->geometricalPropagator().propagateWithPath(updtsos, *hits[ihit+1]->surface());
        if (!std::get<0>(propresult)) {
          std::cout << "Abort: Propagation Failed!" << std::endl;
          valid = false;
          break;
        }
        
//           auto propresultorig = g4prop->propagateGenericWithJacobian(*updtsos.freeState(), surface);
// //           
//           std::cout << "jac" << std::endl;
//           std::cout << std::get<1>(propresult) << std::endl;
//           std::cout << "jacorig" << std::endl;
//           std::cout << std::get<1>(propresultorig) << std::endl;
//           
//           std::cout << "err" << std::endl;
//           std::cout << std::get<0>(propresult).localError().matrix() << std::endl;
//           std::cout << "errorig" << std::endl;
//           std::cout << std::get<0>(propresultorig).localError().matrix() << std::endl;
        
//           std::cout << "pPre = " << updtsos.globalMomentum().mag() << " pPost = " << std::get<0>(propresult).globalMomentum().mag() << std::endl;
        
//           assert(std::get<0>(propresult).globalMomentum().mag() <= updtsos.globalMomentum().mag());
        
//         const Matrix<double, 5, 7> FdFm = Map<const Matrix<double, 5, 7, RowMajor>>(std::get<1>(propresult).Array());
//           FdFm = localTransportJacobian(updtsos, propresult, false);
        
//         const Matrix<double, 5, 5> dQ = Map<const Matrix<double, 5, 5, RowMajor>>(std::get<2>(propresult).Array());
//         const Matrix<double, 5, 5> dQ = Matrix<double, 5, 5>::Zero();
        
        updtsos = std::get<1>(propresult);
        const Matrix<double, 5, 5> Qcurv = std::get<2>(propresult);
        const Matrix<double, 5, 7> FdFmcurv = std::get<3>(propresult);
        const double dEdxlast = std::get<4>(propresult);
//         const Matrix<double, 5, 5> dQcurv = std::get<5>(propresult);
        
        Matrix<double, 5, 7> FdFm = FdFmcurv;
        if (ihit == 0) {
          // extra jacobian from reference state to curvilinear potentially needed
          FdFm.leftCols<5>() = (FdFm.leftCols<5>()*ref2curvjac).eval();
        }

//         if (ihit == (hits.size() - 1)) {
//           dEdxout = dEdxlast;
//         }
        Qtot = (FdFm.leftCols<5>()*Qtot*FdFm.leftCols<5>().transpose()).eval();
        Qtot += Qcurv;






//         auto propresultqopup = g4prop->propagateGenericWithJacobianAltD(tsosqopup, surface, dbetaval, dxival);
//         const Matrix<double, 5, 7> FdFmqopup = std::get<3>(propresultqopup);
//
//         auto propresultqopdown = g4prop->propagateGenericWithJacobianAltD(tsosqopdown, surface, dbetaval, dxival);
//         const Matrix<double, 5, 7> FdFmqopdown = std::get<3>(propresultqopdown);
//
//         const double qopden = sigmaqop > 0. ? sigmaqop : 1.;
//         const Matrix<double, 5, 1> d2xdqop2 = (FdFmqopup.col(0) - FdFmqopdown.col(0))/qopden;
//
//         const Matrix<double, 5, 1> curvconv = 0.5*sigmaqop*sigmaqop*d2xdqop2;
//
//         std::cout << "iiter = " << iiter << " ihit = " << ihit << " oldqop = " << oldqop << " sigmaqop = " << sigmaqop << " curvconv:\n" << curvconv << std::endl;

        
//         const GlobalPoint postmp(updtsos[0], updtsos[1], updtsos[2]);
//         const GlobalVector bvtmp = field->inTesla(postmp);
//         
//         std::cout << "postmp" << std::endl;
//         std::cout << postmp << std::endl;
//         std::cout << "bvtmp" << std::endl;
//         std::cout << bvtmp << std::endl;
        
//         if (trackEta > 2.2) {
//           std::cout << "nominal bfield" << std::endl;
//           std::cout << updtsos.magneticField()->inTesla(updtsos.globalPosition()) << std::endl;
//           
//           const double dx = 10e-4;
//           
//           GlobalPoint altpos(updtsos.globalPosition().x(), updtsos.globalPosition().y(), updtsos.globalPosition().z() + dx);
//           
//           const GlobalVector bgrad = (updtsos.magneticField()->inTesla(altpos) - updtsos.magneticField()->inTesla(updtsos.globalPosition()))/dx;
//           
//           std::cout << "bgrad" << std::endl;
//           std::cout << bgrad << std::endl;
//         }
        
        
//         std::cout << "updtsos localposition = " << updtsos.localPosition() << std::endl;
        
        


        
//         std::cout << "FdFm" << std::endl;
//         std::cout << FdFm << std::endl;
        
//         std::cout << "ihit = " << ihit << " p = " << updtsos.globalMomentum().mag() << std::endl;
        
        


        
        // curvilinear to local jacobian
//         JacobianCurvilinearToLocal curv2localm(updtsos.surface(), updtsos.localParameters(), *updtsos.magneticField());
//         const AlgebraicMatrix55& curv2localjacm = curv2localm.jacobian();
//         const Matrix<double, 5, 5> Hm = Map<const Matrix<double, 5, 5, RowMajor>>(curv2localjacm.Array()); 
//         const Matrix<double, 5, 5> Hm = curv2localJacobianAlt(updtsos);

//         const Matrix<double, 5, 5> Hm = curv2localJacobianAltelossD(updtsos, field, surface, dEdxlast, mmu, dbetaval);
        const Matrix<double, 5, 5> Hm = curv2localhybridJacobianAltelossD(updtsos, field, surface, dEdxlast, mmu, dbetaval);

//         const Matrix<double, 5, 5> Hmalt = curv2localJacobianAltelossD(updtsos, field, surface, dEdxlast, mmu, dbetaval);

//         std::cout << "Hm:\n" << Hm << std::endl;
//         std::cout << "Hmalt:\n" << Hmalt << std::endl;

//         const Matrix<double, 5, 11> Hmalign = curv2localJacobianAlignmentD(updtsos, field, surface, dEdxlast, mmu, dbetaval);

//         std::cout << "Hm:\n" << Hm << std::endl;
//         std::cout << "Hmalign:\n" << Hmalign << std::endl;

        // compute convolution correction in local coordinates (BEFORE material effects are applied)
//         const Matrix<double, 2, 1> dxlocalconv = localPositionConvolution(updtsos);

        //get the process noise matrix
//         AlgebraicMatrix55 const Qmat = updtsos.localError().matrix();
//         const Map<const Matrix<double, 5, 5, RowMajor>>Q(Qmat.Array());
        const Matrix<double, 5, 5> Q = dolocalupdate ? Hm*Qcurv*Hm.transpose() : Qcurv;
        
//         const Matrix<double, 5, 5> dQ = Hm*dQcurv*Hm.transpose();
        
//         const Map<const Matrix<double, 5, 5, RowMajor>>Qorig(Qmat.Array());
//         Matrix<double, 5, 5> Q = Qorig;
//         Q(0,0) = 100.*dqop*dqop;
        
//         std::cout<< "Q" << std::endl;
//         std::cout<< "ihit = " << ihit << " Q" << std::endl;
//         std::cout<< Q << std::endl;
//         
//         std::cout<< "dQ" << std::endl;
//         std::cout<< dQ << std::endl;
        
        
        const float enext = simhit == nullptr ? -99. : std::sqrt(std::pow(simhit->pabs(), 2) + mmu*mmu) - 0.5*simhit->energyLoss();        
        const float eprednext = std::sqrt(updtsos.segment<3>(3).squaredNorm() + mmu*mmu);
        
        const float dEval = e > 0. && enext > 0. ? enext - e : -99.;
        const float dEpredval = eprednext - epred;
        
        e = enext;
        epred = eprednext;
        
        const float sigmadEval = std::pow(updtsos.segment<3>(3).norm(), 3)/e*std::sqrt(Q(0, 0));
        

        Matrix<double, 5, 1> localparmsprop;
        const Point3DBase<double, GlobalTag> posprop(updtsos[0], updtsos[1], updtsos[2]);
//         const Vector3DBase<double, GlobalTag> momprop(updtsos[3], updtsos[4], updtsos[5]);
        
        const Point3DBase<double, LocalTag> localpos = surface.toLocal(posprop);
//         const Vector3DBase<double, LocalTag> localmom = surface.toLocal(momprop);
        
//         const Point3DBase<double, LocalTag> localpos = toLocal(surface, posprop);
//         const Vector3DBase<double, LocalTag> localmom = toLocal(surface, momprop);
        
        localparmsprop[0] = updtsos[6]/updtsos.segment<3>(3).norm();
//         localparmsprop[1] = localmom.x()/localmom.z();
//         localparmsprop[2] = localmom.y()/localmom.z();
        localparmsprop[1] = std::atan(updtsos[5]/std::sqrt(updtsos[3]*updtsos[3] + updtsos[4]*updtsos[4]));
        localparmsprop[2] = std::atan2(updtsos[4], updtsos[3]);
        localparmsprop[3] = localpos.x();
        localparmsprop[4] = localpos.y();        
        
        Matrix<double, 5, 1> localparms = localparmsprop;
        
        // update state from previous iteration
        //momentum kink residual
//         AlgebraicVector5 idx0(0., 0., 0., 0., 0.);
        Matrix<double, 5, 1> idx0 = Matrix<double, 5, 1>::Zero();
        
        if (iiter == 0 && fitFromSimParms_) {
          if (simhit == nullptr) {
            valid = false;
            break;
          }

          // alternate version with propagation from entry state

          const Point3DBase<double, LocalTag> simlocalpos = simhit->entryPoint();
          const Vector3DBase<double, LocalTag> simlocalmom = simhit->momentumAtEntry();

//           std::cout << "simlocalpos" << simlocalpos << std::endl;

          const Point3DBase<double, GlobalTag> simglobalpos = surface.toGlobal(simlocalpos);
          const Vector3DBase<double, GlobalTag> simglobalmom = surface.toGlobal(simlocalmom);

          updtsos[0] = simglobalpos.x();
          updtsos[1] = simglobalpos.y();
          updtsos[2] = simglobalpos.z();
          updtsos[3] = simglobalmom.x();
          updtsos[4] = simglobalmom.y();
          updtsos[5] = simglobalmom.z();
          updtsos[6] = genpart->charge();

          auto propresultsim = g4prop->propagateGenericWithJacobianAltD(updtsos, surface, dbetaval, dxival);

          if (!std::get<0>(propresultsim)) {
            std::cout << "Abort: Sim state Propagation Failed!" << std::endl;
            valid = false;
            break;
          }

          updtsos = std::get<1>(propresultsim);

          const Point3DBase<double, GlobalTag> simglobalposprop(updtsos[0], updtsos[1], updtsos[2]);
          const Vector3DBase<double, GlobalTag> simglobalmomprop(updtsos[3], updtsos[4], updtsos[5]);

          const Point3DBase<double, LocalTag> simlocalposprop = surface.toLocal(simglobalposprop);
          const Vector3DBase<double, LocalTag> simlocalmomprop = surface.toLocal(simglobalmomprop);

          localparms[0] = updtsos[6]/updtsos.segment<3>(3).norm();
          localparms[1] = simlocalmomprop.x()/simlocalmomprop.z();
          localparms[2] = simlocalmomprop.y()/simlocalmomprop.z();
          localparms[3] = simlocalposprop.x();
          localparms[4] = simlocalposprop.y();

          idx0 = localparms - localparmsprop;

        }
        
        
        if (dolocalupdate) {
          if (iiter==0) {
  //         if (true) {
            layerStates.push_back(updtsos);
  //           layerStatesStart.push_back(updtsos);


          }
          else {
            //current state from previous state on this layer
            //save current parameters

            Matrix<double, 7, 1>& oldtsos = layerStates[ihit];
//             const Matrix<double, 5, 5> Hold = curv2localJacobianAltelossD(oldtsos, field, surface, dEdxlast, mmu, dbetaval);
            const Matrix<double, 5, 5> Hold = curv2localhybridJacobianAltelossD(oldtsos, field, surface, dEdxlast, mmu, dbetaval);
            const Matrix<double, 5, 1> dxlocal = Hold*dxfull.segment<5>(5*(ihit+1));

            const Point3DBase<double, GlobalTag> pos(oldtsos[0], oldtsos[1], oldtsos[2]);
            const Point3DBase<double, LocalTag> localpos = surface.toLocal(pos);

//             const Point3DBase<double, LocalTag> localposupd(localpos.x() + dxlocal[3], localpos.y() + dxlocal[4], localpos.z());
            const Point3DBase<double, LocalTag> localposupd(localpos.x() + dxlocal[3], localpos.y() + dxlocal[4], 0.);
            const Point3DBase<double, GlobalTag> posupd = surface.toGlobal(localposupd);


//             const Vector3DBase<double, GlobalTag> mom(oldtsos[3], oldtsos[4], oldtsos[5]);
//             const Vector3DBase<double, LocalTag> localmom = surface.toLocal(mom);

//             const double dxdz = localmom.x()/localmom.z();
//             const double dydz = localmom.y()/localmom.z();



//             const double dxdzupd = dxdz + dxlocal[1];
//             const double dydzupd = dydz + dxlocal[2];

            const double qop = oldtsos[6]/oldtsos.segment<3>(3).norm();
            const double lam = std::atan(oldtsos[5]/std::sqrt(oldtsos[3]*oldtsos[3] + oldtsos[4]*oldtsos[4]));
            const double phi = std::atan2(oldtsos[4], oldtsos[3]);

            const double qopupd = qop + dxlocal[0];
            const double lamupd = lam + dxlocal[1];
            const double phiupd = phi + dxlocal[2];

            const double pupd = std::abs(1./qopupd);
            const double charge = std::copysign(1., qopupd);

            const double pxupd = pupd*std::cos(lamupd)*std::cos(phiupd);
            const double pyupd = pupd*std::cos(lamupd)*std::sin(phiupd);
            const double pzupd = pupd*std::sin(lamupd);

//             const double signpz = std::copysign(1., localmom.z());
//             const double localmomfact = signpz/std::sqrt(1. + dxdzupd*dxdzupd + dydzupd*dydzupd);
//             const Vector3DBase<double, LocalTag> localmomupd(pupd*dxdzupd*localmomfact, pupd*dydzupd*localmomfact, pupd*localmomfact);
//             const Vector3DBase<double, GlobalTag> momupd = surface.toGlobal(localmomupd);

            oldtsos[0] = posupd.x();
            oldtsos[1] = posupd.y();
            oldtsos[2] = posupd.z();
            oldtsos[3] = pxupd;
            oldtsos[4] = pyupd;
            oldtsos[5] = pzupd;
//             oldtsos[3] = momupd.x();
//             oldtsos[4] = momupd.y();
//             oldtsos[5] = momupd.z();
            oldtsos[6] = charge;

            updtsos = oldtsos;

            localparms[0] = qopupd;
//             localparms[1] = dxdzupd;
//             localparms[2] = dydzupd;
//             localparms[1] = std::atan(updtsos[5]/std::sqrt(updtsos[3]*updtsos[3] + updtsos[4]*updtsos[4]));
//             localparms[2] = std::atan2(updtsos[4], updtsos[3]);
            localparms[1] = lamupd;
            localparms[2] = phiupd;
            localparms[3] = localposupd.x();
            localparms[4] = localposupd.y();

            idx0 = localparms - localparmsprop;

          }
        }
        
//         if (false) {
//           //current state from previous state on this layer
//           //save current parameters    
//           //debug version
//           
//           std::cout << "update ihit = " << ihit << std::endl;
//           
//           Matrix<double, 7, 1>& oldtsos = layerStates[ihit];
//           
// //           std::cout << "proptsos" << std::endl;
// //           std::cout << updtsos.transpose() << std::endl;
// //           std::cout << "oldtsos" << std::endl;
// //           std::cout << oldtsos.transpose() << std::endl;
//           
// //           JacobianCurvilinearToLocal curv2localold(oldtsos.surface(), oldtsos.localParameters(), *oldtsos.magneticField());
// //           const AlgebraicMatrix55& curv2localjacold = curv2localold.jacobian();
// //           const Matrix<double, 5, 5> Hold = Map<const Matrix<double, 5, 5, RowMajor>>(curv2localjacold.Array()); 
// //           const Matrix<double, 5, 5> Hold = curv2localJacobianAlt(oldtsos);
//           const Matrix<double, 5, 5> Hold = curv2localJacobianAltelossD(oldtsos, field, surface, dEdxlast, mmu);
// //           const Matrix<double, 6, 5> jacold = curv2cartJacobianAltD(oldtsos);
// 
//           
//           
// //           auto const& dxlocal = Hold*dxfull.segment<5>(5*(ihit+1));
//           
//           const Matrix<double, 5, 1> dxlocal = Hold*dxfull.segment<5>(5*(ihit+1));
//           
// //           const Matrix<double, 6, 1> dxcart = jacold*dxfull.segment<5>(5*(ihit+1));
//           
//           
//           
//           std::cout << "dxfull.segment<5>(5*(ihit+1))" << std::endl;
//           std::cout << dxfull.segment<5>(5*(ihit+1)).transpose() << std::endl;
//           std::cout << "dxlocal" << std::endl;
//           std::cout << dxlocal << std::endl;
// //           std::cout << "dxcart" <<std::endl;
// //           std::cout << dxcart << std::endl;
//           
//           const Vector3DBase<double, LocalTag> deltaposlocal(dxlocal[3], dxlocal[4], 0.);
//           const Vector3DBase<double, GlobalTag> deltaposglobal = toGlobal(surface, deltaposlocal);
//           
//           std::cout << "deltaposlocal" << std::endl;
//           std::cout << deltaposlocal << std::endl;
//           std::cout << "deltaposglobal" << std::endl;
//           std::cout << deltaposglobal << std::endl;
//           
//           std::cout << "deltaposlocal.mag()" << std::endl;
//           std::cout << deltaposlocal.mag() << std::endl;
//           std::cout << "deltaposglobal.mag()" << std::endl;
//           std::cout << deltaposglobal.mag() << std::endl;
//           
// //           Vector3DBase<double, LocalTag> localx(1., 0., 0.);
// //           Vector3DBase<double, LocalTag> localy(0., 1., 0.);
//           
//           const Point3DBase<double, GlobalTag> pos(oldtsos[0], oldtsos[1], oldtsos[2]);
// //           const Point3DBase<double, LocalTag> localpos = surface.toLocal<double>(pos);
//           const Point3DBase<double, LocalTag> localpos = toLocal(surface, pos);
//           
//           const Point3DBase<double, LocalTag> localposupd(localpos.x() + dxlocal[3], localpos.y() + dxlocal[4], localpos.z());
// //           const Point3DBase<double, GlobalTag> posupd = surface.toGlobal<double>(localposupd);
//           const Point3DBase<double, GlobalTag> posupd = toGlobal(surface, localposupd);
//           const Point3DBase<double, GlobalTag> postest = toGlobal(surface, localpos);
//          
//           
//           const Vector3DBase<double, GlobalTag> mom(oldtsos[3], oldtsos[4], oldtsos[5]);
// //           const Vector3DBase<double, LocalTag> localmom = surface.toLocal<double>(mom);
//           const Vector3DBase<double, LocalTag> localmom = toLocal(surface, mom);
//           
//           
//           std::cout << "localposupd - localpos" << std::endl;
//           std::cout << localposupd - localpos << std::endl;
//           std::cout << "posupd - pos" << std::endl;
//           std::cout << posupd - pos << std::endl;
//           std::cout << "postest - pos" << std::endl;
//           std::cout << postest - pos << std::endl;
//           
//           Matrix<double, 3, 3> rot;
//           rot << surface.rotation().xx(), surface.rotation().xy(), surface.rotation().xz(),
//                   surface.rotation().yx(), surface.rotation().yy(), surface.rotation().yz(),
//                   surface.rotation().zx(), surface.rotation().zy(), surface.rotation().zz();
//                   
//           std::cout << "rotpre" << std::endl;
//           std::cout << rot << std::endl;
//                   
//           for (unsigned int i = 0; i < 3; ++i) {
//             for (unsigned int j = 0; j < i; ++j) {
//               rot.row(i) = (rot.row(i) - (rot.row(i)*rot.row(j).transpose())[0]/rot.row(j).squaredNorm()*rot.row(j)).eval();
//             }
//           }
//           for (unsigned int i = 0; i < 3; ++i) {
//               rot.row(i) = (rot.row(i)/rot.row(i).norm()).eval();
//           }
//           
//                   
//           const Matrix<double, 3, 3> Itest = rot*rot.transpose();
//           
//           std::cout << "rot" << std::endl;
//           std::cout << rot << std::endl;
//           std::cout << "Itest" << std::endl;
//           std::cout << Itest << std::endl;
//           
//           const Matrix<double, 3, 1> deltposlocaleig(dxlocal[3], dxlocal[4], 0.);
//           const Matrix<double, 3, 1> deltposglobaleig = rot.transpose()*deltposlocaleig;
//           
//           std::cout << "deltposlocaleig" << std::endl;
//           std::cout << deltposlocaleig.transpose() << std::endl;
//           std::cout << "deltposglobaleig" << std::endl;
//           std::cout << deltposglobaleig.transpose() << std::endl;
//           
//           std::cout << "deltposlocaleig.norm()" << std::endl;
//           std::cout << deltposlocaleig.norm() << std::endl;
//           std::cout << "deltposglobaleig.norm()" << std::endl;
//           std::cout << deltposglobaleig.norm() << std::endl;
//           
//           const Matrix<double, 3, 1> sposeig(surface.position().x(), surface.position().y(), surface.position().z());
//           
//           const Matrix<double, 3, 1> poseig(oldtsos[0], oldtsos[1], oldtsos[2]);
//           const Matrix<double, 3, 1> localposeig = rot*(poseig - sposeig);
//           const Matrix<double, 3, 1> postesteig = rot.transpose()*localposeig + sposeig;
//           
//           std::cout << "postesteig - poseig" << std::endl;
//           std::cout << (postesteig - poseig).transpose() << std::endl;
//           
//           
//           
//           
//           const double dxdz = localmom.x()/localmom.z();
//           const double dydz = localmom.y()/localmom.z();
//           
//           
//           
//           const double dxdzupd = dxdz + dxlocal[1];
//           const double dydzupd = dydz + dxlocal[2];
//           
//           const double qop = oldtsos[6]/oldtsos.segment<3>(3).norm();
//           const double qopupd = qop + dxlocal[0];
//           
//           const double pupd = std::abs(1./qopupd);
//           const double charge = std::copysign(1., qopupd);
//           
//           const double signpz = std::copysign(1., localmom.z());
//           const double localmomfact = signpz/std::sqrt(1. + dxdzupd*dxdzupd + dydzupd*dydzupd);
//           const Vector3DBase<double, LocalTag> localmomupd(pupd*dxdzupd*localmomfact, pupd*dydzupd*localmomfact, pupd*localmomfact);
// //           const Vector3DBase<double, GlobalTag> momupd = surface.toGlobal<double>(localmomupd);
//           const Vector3DBase<double, GlobalTag> momupd = toGlobal(surface, localmomupd);
//           
//           const Matrix<double, 7, 1> oldtsosorig = oldtsos;
//           
//           oldtsos[0] = posupd.x();
//           oldtsos[1] = posupd.y();
//           oldtsos[2] = posupd.z();
//           oldtsos[3] = momupd.x();
//           oldtsos[4] = momupd.y();
//           oldtsos[5] = momupd.z();
//           oldtsos[6] = charge;          
//           
//           updtsos = oldtsos;
//           
//           std::cout << "oldtsos:" << std::endl;
//           std::cout << oldtsosorig.transpose() << std::endl;
//           std::cout << "updtsos:" << std::endl;
//           std::cout << updtsos.transpose() << std::endl;
//           std::cout << "diff:" << std::endl;
//           std::cout << (updtsos - oldtsosorig).transpose() << std::endl;
//           
// //           std::cout << "updtsos" << std::endl;
// //           std::cout << updtsos.transpose() << std::endl;
// //           
// //           std::cout << "updtsos p = " << updtsos.segment<3>(3).norm() << " target p = " << pupd << std::endl;
//           
//           localparms[0] = qopupd;
//           localparms[1] = dxdzupd;
//           localparms[2] = dydzupd;
//           localparms[3] = localposupd.x();
//           localparms[4] = localposupd.y();
//           
//           idx0 = localparms - localparmsprop;
//           
//           
//         }
        
        //apply measurement update if applicable
//         std::cout << "constructing preciseHit" << std::endl;
//         const Matrix<double, 7, 1> starttsos = layerStatesStart[ihit];
//         const GlobalPoint pos(starttsos[0], starttsos[1], starttsos[2]);
//         const GlobalVector mom(starttsos[3], starttsos[4], starttsos[5]);
//         const GlobalTrajectoryParameters glob(pos, mom, starttsos[6], field);
        const GlobalPoint pos(updtsos[0], updtsos[1], updtsos[2]);
        const GlobalVector mom(updtsos[3], updtsos[4], updtsos[5]);
        const GlobalTrajectoryParameters glob(pos, mom, updtsos[6], field);
//         const TrajectoryStateOnSurface tsostmp(glob, surface);
        const TrajectoryStateOnSurface tsostmp(glob, *hit->surface());

        auto const& preciseHit = hit->isValid() ? cloner.makeShared(hit, tsostmp) : hit;
        if (hit->isValid() && !preciseHit->isValid()) {
          std::cout << "Abort: Failed updating hit" << std::endl;
          valid = false;
          break;
        }
        
//         const uint32_t gluedid = trackerTopology->glued(preciseHit->det()->geographicalId());
//         const bool isglued = gluedid != 0;
//         const DetId parmdetid = isglued ? DetId(gluedid) : preciseHit->geographicalId();
//         const bool align2d = detidparms.count(std::make_pair(1, parmdetid));
//         const GeomDet* parmDet = isglued ? globalGeometry->idToDet(parmdetid) : preciseHit->det();
        
//         const bool align2d = detidparms.count(std::make_pair(1, preciseHit->geographicalId()));
        const bool align2d = detidparms.count(std::make_pair(1, parmdetid));

        const Matrix<double, 2, 2> &Rglued = rgluemap_.at(preciseHit->geographicalId());
        const GloballyPositioned<double> &surfaceglued = surfacemapD_.at(parmdetid);
        
        // curvilinear to local jacobian
//         JacobianCurvilinearToLocal curv2localp(updtsos.surface(), updtsos.localParameters(), *updtsos.magneticField());
//         const AlgebraicMatrix55& curv2localjacp = curv2localp.jacobian();
//         const Matrix<double, 5, 5> Hp = Map<const Matrix<double, 5, 5, RowMajor>>(curv2localjacp.Array()); 
//         const Matrix<double, 5, 5> Hp = curv2localJacobianAlt(updtsos);

//         const Matrix<double, 5, 5> Hp = curv2localJacobianAltelossD(updtsos, field, surface, dEdxlast, mmu, dbetaval);
        const Matrix<double, 5, 5> Hp = curv2localhybridJacobianAltelossD(updtsos, field, surface, dEdxlast, mmu, dbetaval);
//         const Matrix<double, 5, 11> Hpalign = curv2localJacobianAlignmentD(updtsos, field, surface, dEdxlast, mmu, dbetaval);

//         const Matrix<double, 5, 11> dHalign = Hpalign - Hmalign;

//         const Matrix<double, 6, 1> aligngradprop = 2.*(Hpalign.rightCols<6>() - Hmalign.rightCols<6>()).transpose()*Q.inverse()*idx0;

//         const Matrix<double, 6, 6> alignhessprop = 2.*(Hpalign.rightCols<6>() - Hmalign.rightCols<6>()).transpose()*Q.inverse()*(Hpalign.rightCols<6>() - Hmalign.rightCols<6>());

//         std::cout << "genEta = " << genEta << " genPt = " << genPt << " iiter = " << iiter << " ihit = " << ihit << " valid = " << hit->isValid() << std::endl;
//         std::cout << "idx0:\n" << idx0 << std::endl;
//         std::cout << "Hpalign:\n" << Hpalign << std::endl;
//         std::cout << "Hmalign:\n" << Hmalign << std::endl;
//         std::cout << "dHalign:\n" << dHalign << std::endl;
//         std::cout << "aligngradprop:\n" << aligngradprop << std::endl;
//         std::cout << "alignhessprop:\n" << alignhessprop << std::endl;
        
//         const Matrix<double, 2, 8> Hpalign = curv2localJacobianAltelossDalign(updtsos, field, surface, dEdxlast, mmu, dbetaval);

        
        const Matrix<double, 5, 5> curvcov = covfull.block<5, 5>(5*(ihit+1), 5*(ihit+1));
//         
//         const Matrix<double, 2, 1> localconv = localPositionConvolutionD(updtsos, Qtot, surface);
        const Matrix<double, 2, 1> localconv = localPositionConvolutionD(updtsos, curvcov, surface);
//         
//         std::cout << "iiter = " << iiter << " ihit = " << ihit << " localconv\n" << localconv << std::endl;
        
//         if (abs(localconv[0]) > 5e-3 || abs(localconv[1]) > 5e-3) {
//           std::cout << "localconv" << std::endl;
//           std::cout << localconv << std::endl;
//         }
        

        //FIXME take care of this elsewhere for the moment
//         const bool genconstraint = dogen && ihit==0;
//         const bool genconstraint = false;
        
        if (true) {

          //momentum kink residual
//           const Vector5d dx0 = Map<const Vector5d>(idx0.Array());
          const Vector5d dx0 = idx0;
          

//           const uint32_t gluedidip1 = trackerTopology->glued(hits[ihit + 1]->det()->geographicalId());
//           const bool isgluedip1 = gluedidip1 != 0;
//           const DetId parmdetidip1 = isgluedip1 ? DetId(gluedidip1) : hits[ihit + 1]->geographicalId();
//           const unsigned int bfieldidx = detidparms.at(std::make_pair(6, parmdetidip1));
//           fieldOffset->setOffset(corparms_[bfieldidx]);
          
//           if (ihit==0) {
//             FreeTrajectoryState tmpfts(updtsospos, updtsos.globalParameters().momentum(), updtsos.charge(), field);
//             propresult = fPropagator->geometricalPropagator().propagateWithPath(tmpfts, *hits[ihit+1]->surface());
//           }
//           else {
//             propresult = fPropagator->geometricalPropagator().propagateWithPath(updtsos, *hits[ihit+1]->surface());
//           }
          

// //           auto const &surfaceip1 = *hits[ihit+1]->surface();
//           auto const &surfaceip1 = *surfacemap_.at(hits[ihit+1]->geographicalId());
//           propresult = thePropagator->propagateWithPath(updtsos, surfaceip1);
// //           propresult = fPropagator->geometricalPropagator().propagateWithPath(updtsos, *hits[ihit+1]->surface());
//           if (!propresult.first.isValid()) {
//             std::cout << "Abort: Propagation Failed!" << std::endl;
//             valid = false;
//             break;
//           }
          
//           const double dqopnext = propresult.first.signedInverseMomentum() - updtsos.signedInverseMomentum();
          
          if (true) {

            constexpr unsigned int nlocalstate = 10;
            constexpr unsigned int nlocalbfield = 1;
            constexpr unsigned int nlocaleloss = 1;
            constexpr unsigned int nlocalparms = nlocalbfield + nlocaleloss;
            
            constexpr unsigned int nlocal = nlocalstate + nlocalbfield + nlocaleloss;
            
            using MSScalar = AANT<double, nlocal>;
            
            constexpr unsigned int localstateidx = 0;

            constexpr unsigned int localparmidx = localstateidx + nlocalstate;
            
            const unsigned int fullstateidx = 5*ihit;
            const unsigned int fullparmidx = nstateparms + parmidx;
            
            
            Matrix<MSScalar, 5, 5> Fstate = FdFm.leftCols<5>().cast<MSScalar>();
            Matrix<MSScalar, 5, 1> Fb = FdFm.col(5).cast<MSScalar>();
            Matrix<MSScalar, 5, 1> Fxi = FdFm.col(6).cast<MSScalar>();
            
            Matrix<MSScalar, 5, 5> Hmstate = Hm.cast<MSScalar>();
            Matrix<MSScalar, 5, 5> Hpstate = Hp.cast<MSScalar>();
            
            Matrix<MSScalar, 5, 5> Qinv = Q.inverse().cast<MSScalar>();
                                    
            // initialize active scalars for state parameters
            Matrix<MSScalar, 5, 1> dum = Matrix<MSScalar, 5, 1>::Zero();
            //suppress gradients of reference point parameters when fitting with gen constraint
            for (unsigned int j=0; j<dum.size(); ++j) {
              init_twice_active_var(dum[j], nlocal, localstateidx + j);
            }

            Matrix<MSScalar, 5, 1> du = Matrix<MSScalar, 5, 1>::Zero();
            for (unsigned int j=0; j<du.size(); ++j) {
              init_twice_active_var(du[j], nlocal, localstateidx + 5 + j);
            }
            
            MSScalar dbeta(0.);
            init_twice_active_var(dbeta, nlocal, localparmidx);

            MSScalar dxi(0.);
            init_twice_active_var(dxi, nlocal, localparmidx + 1);
                        
//             const Matrix<MSScalar, 5, 1> dprop = dolocalupdate ? dx0.cast<MSScalar>() + Hpstate*du - Hmstate*Fstate*dum - Hmstate*Fb*dbeta - Hmstate*Fxi*dxi : du - Fstate*dum - Fb*dbeta - Fxi*dxi;
//             const Matrix<MSScalar, 5, 1> dprop = dx0.cast<MSScalar>() + du - Fstate*dum - Fb*dbeta;
            
            Matrix<MSScalar, 5, 1> dprop;
            if (dolocalupdate) {
              dprop = dx0.cast<MSScalar>() + Hpstate*du - Hmstate*Fstate*dum - Hmstate*Fb*dbeta - Hmstate*Fxi*dxi;
            }
            else {
              dprop = du - Fstate*dum - Fb*dbeta - Fxi*dxi;
            }
            
            
//             const MSScalar chisq = dprop.transpose()*Qinv*dprop;
            
            MSScalar chisq;
            
            if (!islikelihood) {
              chisq = dprop.transpose()*Qinv*dprop;
            }
            else {
//               const Matrix<MSScalar, 5, 5> dQdqop = dQ.cast<MSScalar>();
//               const Matrix<MSScalar, 5, 5> dQdxi = dQ.cast<MSScalar>();
//               const Matrix<MSScalar, 5, 5> Qinvfull = (1. - dxi)*Qinv;
//               const Matrix<MSScalar, 5, 5> Qinvfull = Qinv/(1. + dxi);
              const Matrix<MSScalar, 5, 5> Qinvfull = Qinv*(1. - dxi);
//               const Matrix<MSScalar, 5, 5> Qinvfull = Qinv - dxi*Qinv*dQdxi*Qinv;
              chisq = dprop.transpose()*Qinvfull*dprop;
//               chisq = dprop.transpose()*Qinvfull*dprop + 5.*internal::log(1. + dxi);
                            
              const MSScalar dchisq = -dprop.transpose()*Qinv*dprop;
//               const MSScalar dchisq = -dprop.transpose()*Qinv*dQdxi*Qinv*dprop;
              
              //TODO should this be 11x11 instead?
              //TODO check additional factor of 2
              
//               std::cout << "extra xi grad = " << dchisq.value().value() << std::endl;
              
              
              Matrix<double, nlocal, nlocal> dhesslocal;
              for (unsigned int j=0; j<nlocal; ++j) {
                dhesslocal.row(j) = dchisq.derivatives()[j].derivatives();
              }
              
//               std::cout << "ihit = " << ihit << " dhesslocal:" << std::endl;
//               std::cout << dhesslocal << std::endl;

//               const unsigned int fullxiidx = fullparmidx + 1;
              
              dhessv[ihit].block<nlocalstate, nlocalstate>(fullstateidx, fullstateidx) += dhesslocal.topLeftCorner<nlocalstate, nlocalstate>();
                            
//               for (unsigned int i=1; i<nlocal; ++i) {
//                 const unsigned int ifull = i < nlocalstate ? fullstateidx + i : (fullparmidx + i) - nlocalstate;
//                 for (unsigned int j=0; j<nlocalstate; ++j) {
//                   const unsigned int jfull = fullstateidx + j;
//                   for (unsigned int k=0; k<nlocalstate; ++k) {
//                     const unsigned int kfull = fullstateidx + k;
//                     if (j==0) {
//                       dhessv[ifull](jfull, kfull) += dhesslocal(i, k);
//                     }
//                     else if (k==0) {
//                       dhessv[ifull](jfull, kfull) += dhesslocal(i, j);
//                     }
//                   }
//                 }  
//               }
              
//               for (unsigned int j=i; j<nlocal; ++j) {
//                 const unsigned int jfull = j < nlocal ? fullstateidx + j : (fullparmidx + j) - nlocalstate;
//                 for (unsigned int k=j; k<nlocal; ++k) {
//                   const unsigned int kfull = k < nlocal ? fullstateidx + k : (fullparmidx + k) - nlocalstate;
//                   const double d3val = dchisq.derivatives()[j].derivatives()[k];
//                   
//                   
//                   if (j < nlocalstate && k < nlocalstate) {
//                     dhessv[ifull](jfull, kfull) += d3val;
//                     if (jfull != kfull) {
//                       dhessv[ifull](kfull, jfull) += d3val;
//                     }
//                   }
// //                   if (j < nlocalstate) {
// //                     if (ifull != kfull) {
// //                       dhessv[kfull](ifull, jfull) += d3val;
// //                       if (ifull != jfull) {
// //                         dhessv[kfull](jfull, ifull) += d3val;
// //                       }
// //                     }
// //                   }
// //                   if (k < nlocalstate) {
// //                     if (jfull != ifull) {
// //                       dhessv[jfull](ifull, kfull) += d3val;
// //                       if (ifull != kfull) {
// //                         dhessv[jfull](kfull, ifull) += d3val;
// //                       }
// //                     }
// //                   }
//                 }
//               }
              
//               std::cout << "dhessv" << std::endl;
//               std::cout << dhessv[ihit] << std::endl;
                
              
              
            }
            
            chisq0val += chisq.value().value();
            
//             std::cout << "ihit = " << ihit << " MS chisq = " << chisq.value().value() << std::endl;
            
//             auto const& gradlocal = chisq.value().derivatives();
            Matrix<double, nlocal, 1> gradlocal = chisq.value().derivatives();
            //fill local hessian
            Matrix<double, nlocal, nlocal> hesslocal;
            for (unsigned int j=0; j<nlocal; ++j) {
              hesslocal.row(j) = chisq.derivatives()[j].derivatives();
            }
            
            if (islikelihood) {
              const double fact = 1./(1. - dxi.value().value());
              gradlocal[localparmidx + 1] += 5.*fact;
              hesslocal(localparmidx + 1, localparmidx + 1) += 5.*fact*fact;
            }
            
//             std::cout << "material: ihit = " << ihit << std::endl;
//             std::cout << "dx0" << std::endl;
//             std::cout << dx0.transpose() << std::endl;
//             std::cout << "Q" << std::endl;
//             std::cout << Q << std::endl;
//             std::cout << "gradlocal" << std::endl;
//             std::cout << gradlocal.transpose() << std::endl;
            
//             std::cout << "gradlocal" << std::endl;
//             std::cout << gradlocal << std::endl;
//             std::cout << "hesslocal" << std::endl;
//             std::cout << hesslocal << std::endl;
            
            //fill global gradient
            gradfull.segment<nlocalstate>(fullstateidx) += gradlocal.head<nlocalstate>();
            gradfull.segment<nlocalparms>(fullparmidx) += gradlocal.segment<nlocalparms>(localparmidx);

            //fill global hessian (upper triangular blocks only)
            hessfull.block<nlocalstate,nlocalstate>(fullstateidx, fullstateidx) += hesslocal.topLeftCorner<nlocalstate,nlocalstate>();
            hessfull.block<nlocalstate,nlocalparms>(fullstateidx, fullparmidx) += hesslocal.topRightCorner<nlocalstate, nlocalparms>();
            hessfull.block<nlocalparms, nlocalparms>(fullparmidx, fullparmidx) += hesslocal.bottomRightCorner<nlocalparms, nlocalparms>();
            

          }
          
//           const unsigned int bfieldglobalidx = detidparms.at(std::make_pair(6, parmdetid));
          globalidxv[parmidx] = bfieldglobalidx;
          parmidx++;
          
//           const unsigned int elossglobalidx = detidparms.at(std::make_pair(7, parmdetid));
          globalidxv[parmidx] = elossglobalidx;
          parmidx++;
          
          //backwards propagation jacobian (local to local) to be used at the next layer
//           FdFm = curv2curvTransportJacobian(*updtsos.freeState(), propresult, false);
//           FdFm = localTransportJacobian(updtsos, propresult, false);
//           dqop = dqopnext;
          
        }

        if (preciseHit->isValid()) {

          auto fillAlignGrads = [&](auto Nalign) {
            constexpr unsigned int nlocalstate = 2;
            constexpr unsigned int localstateidx = 0;
            constexpr unsigned int localalignmentidx = nlocalstate;
            constexpr unsigned int localparmidx = localalignmentidx;

            // abusing implicit template argument to pass
            // a template value via std::integral_constant
            constexpr unsigned int nlocalalignment = Nalign();
            constexpr unsigned int nlocalparms = nlocalalignment;
            constexpr unsigned int nlocal = nlocalstate + nlocalparms;

            using AlignScalar = AANT<double, nlocal>;
            
            const unsigned int fullstateidx = 5*(ihit+1) + 3;
            const unsigned int fullparmidx = nstateparms + nparsBfield + nparsEloss + alignmentparmidx;

            const bool ispixel = GeomDetEnumerators::isTrackerPixel(preciseHit->det()->subDetector());

            //TODO add hit validation stuff
            //TODO add simhit stuff
            
//             const PSimHit *matchedsim = nullptr;
//             for (auto const& simhith : simHits) {
//               for (const PSimHit& simHit : *simhith) {
//                 if (std::abs(simHit.particleType()) == 13 && simHit.detUnitId() == preciseHit->geographicalId()) {
//                   matchedsim = &simHit;
//                   break;
//                 }
//               }
//               if (matchedsim) {
//                 break;
//               }
//             }
            
//             const bool hit1d = preciseHit->dimension() == 1 || (detidlayermap.at(preciseHit->geographicalId())[0] == 0 && detidlayermap.at(preciseHit->geographicalId())[1] == 1);
//             const bool hit1d = preciseHit->dimension() == 1 || ispixel;
            
            const bool hit1d = preciseHit->dimension() == 1;
            
            
            Matrix<AlignScalar, 2, 2> Hu = Hp.bottomRightCorner<2,2>().cast<AlignScalar>();

            Matrix<AlignScalar, 2, 1> dy0;
            Matrix<AlignScalar, 2, 2> Vinv;
            // rotation from module to strip coordinates
//             Matrix<AlignScalar, 2, 2> R;
            Matrix2d R;
            
            const double lxcor = localparms[3];
            const double lycor = localparms[4];

//             const double lxcor = localparms[3] - localconv[0];
//             const double lycor = localparms[4] - localconv[1];

            const Topology &topology = preciseHit->det()->topology();

            // undo deformation correction
            const LocalPoint lpnull(0., 0.);
            const MeasurementPoint mpnull = topology.measurementPosition(lpnull);
            const Topology::LocalTrackPred pred(tsostmp.localParameters().vector());

            auto const defcorr = topology.localPosition(mpnull, pred) - topology.localPosition(mpnull);

            const double hitx = preciseHit->localPosition().x() - defcorr.x();
            const double hity = preciseHit->localPosition().y() - defcorr.y();

            double lyoffset = 0.;
            double hitphival = -99.;
            double localphival = -99.;

//             if (preciseHit->dimension() == 1) {
            if (hit1d) {
              const ProxyStripTopology *proxytopology = dynamic_cast<const ProxyStripTopology*>(&(preciseHit->det()->topology()));

//               std::cout << "1d hit" << std::endl;
//               assert(!align2d);
//               dy0[0] = AlignScalar(matchedsim->localPosition().x() - updtsos.localPosition().x());
              
//               dy0[0] = AlignScalar(preciseHit->localPosition().x() - localparms[3]);
//               dy0[0] = AlignScalar(preciseHit->localPosition().x() - lxcor);
              dy0[0] = AlignScalar(hitx - lxcor);
//               dy0[0] = AlignScalar(preciseHit->localPosition().x() - localparms[3] + localconv[0]);
              dy0[1] = 0. - lycor;
//               dy0[1] = AlignScalar(0.);
              
//               bool simvalid = false;
//               for (auto const& simhith : simHits) {
//                 for (const PSimHit& simHit : *simhith) {
//                   if (simHit.detUnitId() == preciseHit->geographicalId()) {                      
//                     
//                     dy0[0] = AlignScalar(simHit.localPosition().x() - updtsos.localPosition().x());
//                     
//                     simvalid = true;
//                     break;
//                   }
//                 }
//                 if (simvalid) {
//                   break;
//                 }
//               }

              const double striplength = proxytopology->stripLength();
              const double yerr2 = striplength*striplength/12.;
              
              Vinv = Matrix<AlignScalar, 2, 2>::Zero();
              Vinv(0,0) = 1./preciseHit->localPositionError().xx();
//               Vinv(1,1) = 1./yerr2;
              
//               R = Matrix<AlignScalar, 2, 2>::Identity();
              R = Matrix2d::Identity();
            }
            else {
              // 2d hit
//               assert(align2d);
              
              Matrix2d iV;
              iV << preciseHit->localPositionError().xx(), preciseHit->localPositionError().xy(),
                    preciseHit->localPositionError().xy(), preciseHit->localPositionError().yy();
              if (ispixel) {
//               if (true) {
//                 std::cout << "2d pixel hit, subdet = " << preciseHit->det()->subDetector() << " track eta = " << trackEta << std::endl;
//                 std::cout << iV << std::endl;
                
                //take 2d hit as-is for pixels
//                 dy0[0] = AlignScalar(matchedsim->localPosition().x() - updtsos.localPosition().x());
//                 dy0[1] = AlignScalar(matchedsim->localPosition().y() - updtsos.localPosition().y());
                
//                 dy0[0] = AlignScalar(preciseHit->localPosition().x() - localparms[3]);
//                 dy0[1] = AlignScalar(preciseHit->localPosition().y() - localparms[4]);

//                 dy0[0] = AlignScalar(preciseHit->localPosition().x() - lxcor);
//                 dy0[1] = AlignScalar(preciseHit->localPosition().y() - lycor);

                dy0[0] = AlignScalar(hitx - lxcor);
                dy0[1] = AlignScalar(hity - lycor);
                
//                 dy0[0] = AlignScalar(preciseHit->localPosition().x() - localparms[3] + localconv[0]);
//                 dy0[1] = AlignScalar(preciseHit->localPosition().y() - localparms[4] + localconv[1]);
              
                Vinv = iV.inverse().cast<AlignScalar>();
                //FIXME various temporary hacks;
                
//                 dy0[1] = AlignScalar(0.);
//                 Vinv = Matrix<AlignScalar, 2, 2>::Zero();
//                 Vinv(0,0) = 1./preciseHit->localPositionError().xx();
                
//                 if (GeomDetEnumerators::isEndcap(preciseHit->det()->subDetector())) {
//                 if (GeomDetEnumerators::isBarrel(preciseHit->det()->subDetector())) {
//                   PXBDetId detidtest(preciseHit->det()->geographicalId());
//                   int layertest = detidtest.layer();
//                   
//                   if (layertest > 1) {
//                     Vinv = Matrix<AlignScalar, 2, 2>::Zero();
//                   }
//                   
// //                   Vinv = Matrix<AlignScalar, 2, 2>::Zero();
// //                   dy0[0] = AlignScalar(0.);
// //                   dy0[1] = AlignScalar(0.);
//                 }
                
//                 bool simvalid = false;
//                 for (auto const& simhith : simHits) {
//                   for (const PSimHit& simHit : *simhith) {
//                     if (simHit.detUnitId() == preciseHit->geographicalId()) {                      
//                       
//                       if (GeomDetEnumerators::isBarrel(preciseHit->det()->subDetector())) {
//                         dy0[0] = AlignScalar(simHit.localPosition().x() - updtsos.localPosition().x());
//                         dy0[1] = AlignScalar(simHit.localPosition().y() - updtsos.localPosition().y());
//                       }
//                       
// //                       dy0[1] = AlignScalar(0.);
//                       
//                       simvalid = true;
//                       break;
//                     }
//                   }
//                   if (simvalid) {
//                     break;
//                   }
//                 }
                
                
//                 R = Matrix<AlignScalar, 2, 2>::Identity();
                R = Matrix2d::Identity();
              }
              else {
                constexpr bool dopolar = true;
//                 constexpr bool dopolar = false;

                if (dopolar) {
                  // transform to polar coordinates to end the madness
                  //TODO handle the module deformations consistently here (currently equivalent to dropping/undoing deformation correction)

//                   std::cout << "wedge\n" << std::endl;

                  const ProxyStripTopology *proxytopology = dynamic_cast<const ProxyStripTopology*>(&(preciseHit->det()->topology()));

//                   // undo deformation correction
//                   const Topology::LocalTrackPred pred(tsostmp.localParameters().vector());
//
//                   auto const defcorr = proxytopology->localPosition(0., pred) - proxytopology->localPosition(0.);
//
//                   const double hitx = preciseHit->localPosition().x() - defcorr.x();
//                   const double hity = preciseHit->localPosition().y() - defcorr.y();

                  const TkRadialStripTopology *radialtopology = dynamic_cast<const TkRadialStripTopology*>(&proxytopology->specificTopology());

//                   if (preciseHit->det()->geographicalId() == 470438309 || preciseHit->det()->geographicalId() == 470438310) {
//                     std::cout << "detid = " << preciseHit->det()->geographicalId().rawId() << " yAxisOrientation = " << radialtopology->yAxisOrientation() << " yCentreOfStripPlane = " << radialtopology->yCentreOfStripPlane() << std::endl;
//                   }

                  const double rdir = radialtopology->yAxisOrientation();
                  const double radius = radialtopology->originToIntersection();

                  const double phihit = rdir*std::atan2(hitx, rdir*hity + radius);
                  const double rhohit = std::sqrt(hitx*hitx + std::pow(rdir*hity + radius, 2));

                  // invert original calculation of covariance matrix to extract variance on polar angle
                  const double detHeight = radialtopology->detHeight();
                  const double radsigma = detHeight*detHeight/12.;

                  const double t1 = std::tan(phihit);
                  const double t2 = t1*t1;

                  const double tt = preciseHit->localPositionError().xx() - t2*radsigma;

                  const double phierr2 = tt / std::pow(radialtopology->centreToIntersection(), 2);

                  const double striplength = detHeight * std::sqrt(1. + std::pow( hitx/(rdir*hity + radius), 2) );

                  const double rhoerr2lin = striplength*striplength/12.;

                  // TODO confirm that rhohit is really at midpoint
                  const double rho0 = rhohit - 0.5*striplength;
                  const double rho1 = rhohit + 0.5*striplength;

                  const double rhobar = (2./3.)*(std::pow(rho1, 3) - std::pow(rho0, 3))/(rho1*rho1 - rho0*rho0);

                  const double rhosqbar = 0.5*(rho1*rho1 + rho0*rho0);

                  const double rhoerr2 = rhosqbar - rhobar*rhobar;

//                   std::cout << "rhohit = " << rhohit << " rhobar = " << rhobar << " rhoerr2lin = " << rhoerr2lin << " rhoerr2 = " << rhoerr2 << std::endl;

                  // TODO apply (inverse) corrections for module deformations here? (take into account for jacobian?)
                  const double phistate = rdir*std::atan2(lxcor, rdir*lycor + radius);
                  const double rhostate = std::sqrt(lxcor*lxcor + std::pow(rdir*lycor + radius, 2));

                  Vinv = Matrix<AlignScalar, 2, 2>::Zero();
                  Vinv(0, 0) = 1./phierr2;
//                   Vinv(1, 1) = 1./rhoerr2;
//                   Vinv(1, 1) = 1./rhoerr2lin;

                  // jacobian from localx-localy to localphi-localy (module bounds are also rectangular in this coordinate system)
                  R = Matrix2d::Zero();

                  const double yp = rdir*lycor + radius;
                  const double invden = 1./(lxcor*lxcor + yp*yp);

                  // dphi / dx
                  R(0, 0) = rdir*yp*invden;
                  // dphi / dy
                  R(0, 1) = -lxcor*invden;
//                   //dy / dy
//                   R(1, 1) = 1.;
                  // drho / dx
                  R(1, 0) = lxcor/rhostate;
                  // drho / dy
                  R(1, 1) = rdir*(rdir*lycor + radius)/rhostate;


                  dy0[0] = phihit - phistate;
//                   dy0[1] = hity - lycor;
//                   dy0[1] = 0.;
//                   dy0[1] = rhobar - rhostate;
                  dy0[1] = rhohit - rhostate;


                  const Matrix<double, 2, 2> Huval = Hp.bottomRightCorner<2,2>();
                  const Matrix<double, 2, 2> statecov = R*Huval*curvcov.bottomRightCorner<2, 2>()*Huval.transpose()*R.transpose();

                  const double phibound = std::abs(radialtopology->phiOfOneEdge());

                  const double sigmaphi = std::sqrt(std::max(statecov(0, 0), 0.));

                  hitphival = phihit;
                  localphival = phistate;

//                   const double phithres = 2e-2;
// //                   const double nsigmaveto = 3.;
// //                   const double nsigmaveto = 10.;
//                   // zero the hit if the state is sufficiently close to the boundary
// //                   if (iiter > 0 && (phistate < (-phibound + phithres) || phistate > (phibound - phithres)  ) ) {
//                   if ((phistate < (-phibound + phithres) || phistate > (phibound - phithres)  ) ) {
// //                   if (iiter > 0 &&  ((phistate - nsigmaveto*sigmaphi) < (-phibound) || (phistate + nsigmaveto*sigmaphi) > phibound) ) {
// //                     std::cout << "veto:  phistate = " << phistate << " sigmaphi = " << sigmaphi << " phibound = " << phibound << std::endl;
//                     Vinv(0, 0) = 0.;
//                     dy0[0] = 0.;
//                     dy0[1] = 0.;
//                   }


                }
                else {
                  // standard treatment

                  // diagonalize and take only smallest eigenvalue for 2d hits in strip wedge modules,
                  // since the constraint parallel to the strip is spurious
                  SelfAdjointEigenSolver<Matrix2d> eigensolver(iV);

                  R = eigensolver.eigenvectors().transpose();
                  if (R(0,0) < 0.) {
                    R.row(0) *= -1.;
                  }
                  if (R(1,1) <0.) {
                    R.row(1) *= -1.;
                  }

                  Matrix<double, 2, 1> dy0local;

                  dy0local[0] = preciseHit->localPosition().x() - lxcor;
                  dy0local[1] = preciseHit->localPosition().y() - lycor;


                  const Matrix<double, 2, 1> dy0eig = R*dy0local;


                  //TODO deal properly with rotations (rotate back to module local coords?)
                  dy0[0] = AlignScalar(dy0eig[0]);
                  dy0[1] = AlignScalar(0.);

                  Vinv = Matrix<AlignScalar, 2, 2>::Zero();
                  Vinv(0,0) = AlignScalar(1./eigensolver.eigenvalues()[0]);
                }
              }
            }
            
            rxfull.row(ivalidhit) = R.row(0).cast<float>();
            ryfull.row(ivalidhit) = R.row(1).cast<float>();
            
            validdxeigjac.block<2,2>(2*ivalidhit, 3*(ihit+1)) = R*Hp.bottomRightCorner<2,2>();
            
            const Matrix<AlignScalar, 2, 2> Ralign = R.cast<AlignScalar>();
            
            Matrix<AlignScalar, 2, 1> dx = Matrix<AlignScalar, 2, 1>::Zero();
            for (unsigned int j=0; j<dx.size(); ++j) {
              init_twice_active_var(dx[j], nlocal, localstateidx + j);
            }
            
            Matrix<AlignScalar, 6, 1> dalpha = Matrix<AlignScalar, 6, 1>::Zero();
            // order in which to use parameters, especially relevant in case nlocalalignment < 6
            constexpr std::array<unsigned int, 6> alphaidxs = {{0, 2, 3, 4, 5, 1}};
            for (unsigned int idim=0; idim<nlocalalignment; ++idim) {
//               init_twice_active_var(dalpha[idim], nlocal, localalignmentidx+idim);
              init_twice_active_var(dalpha[alphaidxs[idim]], nlocal, localalignmentidx+idim);
            }
            
            // alignment jacobian
            Matrix<double, 2, 6> Aval = Matrix<double, 2, 6>::Zero();

            const Point3DBase<double, GlobalTag> posprop(updtsos[0], updtsos[1], updtsos[2]);
            const Point3DBase<double, LocalTag> localposglued = surfaceglued.toLocal(posprop);

            const Vector3DBase<double, GlobalTag> momprop(updtsos[3], updtsos[4], updtsos[5]);
            const Vector3DBase<double, LocalTag> localmomglued = surfaceglued.toLocal(momprop);

            const double localqopval = localparms[0];
//             const double localdxdzval = localparms[1];
//             const double localdydzval = localparms[2];
            const double localdxdzval = localmomglued.x()/localmomglued.z();
            const double localdydzval = localmomglued.y()/localmomglued.z();
//             const double localxval = localparms[3];
//             const double localyval = localparms[4];
//             const double localxval = lxcor;
//             const double localyval = lycor;
            const double localxval = localposglued.x();
            const double localyval = localposglued.y();
//             const double localyval = localparms[4] + lyoffset;
                        

            //standard case

            // dx/dx
            Aval(0,0) = 1.;
            // dy/dy
            Aval(1,1) = 1.;
            // dx/dz
            Aval(0,2) = localdxdzval;
            // dy/dz
            Aval(1,2) = localdydzval;
            // dx/dtheta_x
            Aval(0,3) = -localyval*localdxdzval;
            // dy/dtheta_x
            Aval(1,3) = -localyval*localdydzval;
            // dx/dtheta_y
            Aval(0,4) = -localxval*localdxdzval;
            // dy/dtheta_y
            Aval(1,4) = -localxval*localdydzval;
            // dx/dtheta_z
            Aval(0,5) = -localyval;
            // dy/dtheta_z
            Aval(1,5) = localxval;

            const Matrix<double, 2, 6> Avalglued = Rglued*Aval;

//             Matrix<AlignScalar, 2, 6> A = Aval.cast<AlignScalar>();
            Matrix<AlignScalar, 2, 6> A = Avalglued.cast<AlignScalar>();


//             const Matrix<double, 2, 6> Rlt = R*Aval;
//             std::cout << "Aval:\n" << Aval << std::endl;
//             std::cout << "Rlt:\n" << Rlt << std::endl;



            
            
//             Matrix<double, 2, 6> Atest = Matrix<double, 2, 6>::Zero();
//             // dx/dx
//             Atest(0,0) = 1.;
//             // dy/dy
//             Atest(1,1) = 1.;
//             // dx/dz
//             Atest(0,2) = localdxdzval;
//             // dy/dz
//             Atest(1,2) = localdydzval;
//             // dx/dtheta_x
//             Atest(0,3) = -localyval*localdxdzval;
//             // dy/dtheta_x
//             Atest(1,3) = -localyval*localdydzval;
//             // dx/dtheta_y
//             Atest(0,4) = -localxval*localdxdzval;
//             // dy/dtheta_y
//             Atest(1,4) = -localxval*localdydzval;
//             // dx/dtheta_z
//             Atest(0,5) = -localyval;
//             // dy/dtheta_z
//             Atest(1,5) = localxval;
//             
//             const Matrix<double, 2, 2> Hutest = Hp.bottomRightCorner<2,2>();
//             
//             const Matrix<double, 2, 2> Hutestalt = Hpalign.leftCols<2>();
//             const Matrix<double, 2, 6> Atestalt = Hpalign.rightCols<6>();
//             
//             std::cout << "Hutest\n" << Hutest << std::endl;
//             std::cout << "Hutestalt\n" << Hutestalt << std::endl;
//             std::cout << "Atest\n" << Atest << std::endl;
//             std::cout << "Atestalt\n" << Atestalt << std::endl;


            
            
//             std::cout << "strip local z shift gradient: " << (Ralign*A.col(2))[0].value().value() << std::endl;
            
            // rotation from alignment basis to module local coordinates
//             Matrix<AlignScalar, 2, 2> A;
//             if (isglued) {
//               const GlobalVector modx = preciseHit->det()->surface().toGlobal(LocalVector(1.,0.,0.));
//               const GlobalVector mody = preciseHit->det()->surface().toGlobal(LocalVector(0.,1.,0.));
//               
//               const GlobalVector gluedx = parmDet->surface().toGlobal(LocalVector(1.,0.,0.));
//               const GlobalVector gluedy = parmDet->surface().toGlobal(LocalVector(0.,1.,0.));
//               
//               A(0,0) = AlignScalar(modx.dot(gluedx));
//               A(0,1) = AlignScalar(modx.dot(gluedy));
//               A(1,0) = AlignScalar(mody.dot(gluedx));
//               A(1,1) = AlignScalar(mody.dot(gluedy));
//             }
//             else {
//               A = Matrix<AlignScalar, 2, 2>::Identity();
//             }
// 
//             Matrix<AlignScalar, 2, 1> dh = dy0 - R*Hu*dx - R*A*dalpha;
            
                      
            double thetaincidence = std::asin(1./std::sqrt(std::pow(localdxdzval,2) + std::pow(localdydzval,2) + 1.));
            
//             bool morehitquality = applyHitQuality_ ? thetaincidence > 0.25 : true;
            bool morehitquality = true;
            
            if (morehitquality) {
              nValidHitsFinal++;
              if (ispixel) {
                nValidPixelHitsFinal++;
              }
            }
            else {
              Vinv = Matrix<AlignScalar, 2, 2>::Zero();
            }

            Matrix<AlignScalar, 2, 1> dh = dy0 - Ralign*Hu*dx - Ralign*A*dalpha;
//             Matrix<AlignScalar, 2, 1> dh = dy0 - Ralign*dx - Ralign*A*dalpha;
            AlignScalar chisq = dh.transpose()*Vinv*dh;
            
            
            chisq0val += chisq.value().value();
            
//             std::cout << "ihit = " << ihit << " meas chisq = " << chisq.value().value() << std::endl;

            auto const& gradlocal = chisq.value().derivatives();
            //fill local hessian
            Matrix<double, nlocal, nlocal> hesslocal;
            for (unsigned int j=0; j<nlocal; ++j) {
              hesslocal.row(j) = chisq.derivatives()[j].derivatives();
            }
            
//             std::cout << "hit: ihit = " << ihit << " dy0 = " << dy0[0].value().value() << " " << dy0[1].value().value() << std::endl;
//             std::cout << "gradlocal" << std::endl;
//             std::cout << gradlocal.transpose() << std::endl;
            
            
//             std::cout << "updtsos.localPosition()" << std::endl;
//             std::cout << updtsos.localPosition() << std::endl;
//             std::cout << "A" << std::endl;
//             std::cout << A << std::endl;
//             
//             std::cout << "gradlocal" << std::endl;
//             std::cout << gradlocal.transpose() << std::endl;
//             
//             std::cout << "hesslocal" << std::endl;
//             std::cout << hesslocal << std::endl;
            
//             Matrix<double, nlocal, 1> gradloctest0;
//             Matrix<double, 1, 1> gradloctest1;
//             Matrix<double, 2, 1> gradloctest2;
            
//             std::cout << "nlocalalignment: " << nlocalalignment << " nlocal: " << nlocal << std::endl;
//             std::cout << "gradlocal type: " << typeid(gradlocal).name() << std::endl;
//             std::cout << "gradloctest0 type: " << typeid(gradloctest0).name() << std::endl;
//             std::cout << "gradloctest1 type: " << typeid(gradloctest1).name() << std::endl;
//             std::cout << "gradloctest2 type: " << typeid(gradloctest2).name() << std::endl;
//             
//             std::cout << "nhits: " << nhits << " nvalid: " << nvalid << " nvalidalign2d: " << nvalidalign2d << " ihit: " << ihit << std::endl;
//             std::cout << "gradfull.size(): " << gradfull.size() << " nlocalstate: " << nlocalstate << " fullstateidx: " << fullstateidx << " nlocalparms: " << nlocalparms << " fullparmidx: " << fullparmidx << std::endl;
            
            // FIXME the templated block functions don't work here for some reason
            //fill global gradient
            gradfull.segment<nlocalstate>(fullstateidx) += gradlocal.head(nlocalstate);
            gradfull.segment<nlocalparms>(fullparmidx) += gradlocal.segment(localparmidx, nlocalparms);

            //fill global hessian (upper triangular blocks only)
            hessfull.block<nlocalstate,nlocalstate>(fullstateidx, fullstateidx) += hesslocal.topLeftCorner(nlocalstate,nlocalstate);
            hessfull.block<nlocalstate,nlocalparms>(fullstateidx, fullparmidx) += hesslocal.topRightCorner(nlocalstate, nlocalparms);
            hessfull.block<nlocalparms, nlocalparms>(fullparmidx, fullparmidx) += hesslocal.bottomRightCorner(nlocalparms, nlocalparms);
            
            for (unsigned int idim=0; idim<nlocalalignment; ++idim) {
//               const unsigned int xglobalidx = detidparms.at(std::make_pair(alphaidxs[idim], preciseHit->geographicalId()));
              const unsigned int xglobalidx = detidparms.at(std::make_pair(alphaidxs[idim], parmdetid));
              globalidxv[nparsBfield + nparsEloss + alignmentparmidx] = xglobalidx;
              alignmentparmidx++;
              if (alphaidxs[idim]==0) {
                hitidxv.push_back(xglobalidx);
              }
            }
            
            localqop_iter[ivalidhit] = localqopval;
            localdxdz_iter[ivalidhit] = localdxdzval;
            localdydz_iter[ivalidhit] = localdydzval;
            localx_iter[ivalidhit] = localxval;
            localy_iter[ivalidhit] = localyval;
            
            const double localqopvar = covfull(5*(ihit+1), 5*(ihit+1));
            localqoperr[ivalidhit] = std::sqrt(localqopvar);
            
            if (iiter == 0) {
              
              // fill hit validation information
//               Vector2d dyrecgenlocal;
//               dyrecgenlocal << dy0[0].value().value(), dy0[1].value().value();
//               const Vector2d dyrecgeneig = R*dyrecgenlocal;
//               dxrecgen.push_back(dyrecgeneig[0]);
//               dyrecgen.push_back(dyrecgeneig[1]);
              dxrecgen.push_back(dy0[0].value().value());
              dyrecgen.push_back(dy0[1].value().value());
              
              dxerr.push_back(1./std::sqrt(Vinv(0,0).value().value()));
              dyerr.push_back(1./std::sqrt(Vinv(1,1).value().value()));

              localqop.push_back(localqopval);
              localdxdz.push_back(localdxdzval);
              localdydz.push_back(localdydzval);
              localx.push_back(localxval);
              localy.push_back(localyval);
              
              dEpred.push_back(dEpredval);
              sigmadE.push_back(sigmadEval);
              Epred.push_back(epred);
              E.push_back(e);
              
              hitphi.push_back(hitphival);
              localphi.push_back(localphival);

              const TrackerSingleRecHit* tkhit = dynamic_cast<const TrackerSingleRecHit*>(preciseHit.get());
              assert(tkhit != nullptr);
              
              if (ispixel) {
                const SiPixelCluster& cluster = *tkhit->cluster_pixel();
                const SiPixelRecHit *pixhit = dynamic_cast<const SiPixelRecHit*>(tkhit);
                assert(pixhit != nullptr);
  //               std::cout << "pixel cluster sizeX = " << cluster.sizeX() <<" sizeY = " << cluster.sizeY() << std::endl;
                clusterSize.push_back(cluster.size());
                clusterSizeX.push_back(cluster.sizeX());
                clusterSizeY.push_back(cluster.sizeY());
                clusterCharge.push_back(cluster.charge());
                clusterChargeBin.push_back(pixhit->qBin());
                clusterOnEdge.push_back(pixhit->isOnEdge());
                if (pixhit->hasFilledProb()) {
                  clusterProbXY.push_back(pixhit->clusterProbability(0));
                }
                else {
                  clusterProbXY.push_back(2.);
                }
                
                clusterSN.push_back(-99.);
                stripsToEdge.push_back(-99);
              }
              else {
                const StripTopology* striptopology = dynamic_cast<const StripTopology*>(&(tkhit->det()->topology()));
                const SiStripCluster& cluster = *tkhit->cluster_strip();
  //               siStripClusterInfo_.setCluster(cluster, preciseHit->geographicalId().rawId());
                SiStripClusterInfo clusterInfo = SiStripClusterInfo(cluster, iSetup, preciseHit->geographicalId().rawId());
                clusterSize.push_back(cluster.amplitudes().size());
                clusterSizeX.push_back(cluster.amplitudes().size());
                clusterSizeY.push_back(1);
                clusterCharge.push_back(cluster.charge());
                clusterChargeBin.push_back(-99);
                clusterOnEdge.push_back(-99);
                clusterProbXY.push_back(-99.);
                clusterSN.push_back(clusterInfo.signalOverNoise());


                const uint16_t firstStrip = cluster.firstStrip();
                const uint16_t lastStrip = cluster.firstStrip() + cluster.amplitudes().size() - 1;
                stripsToEdge.push_back(std::min<int>(firstStrip, striptopology->nstrips() - 1 - lastStrip));
              }
              
  //             if (ispixel) {
  //               const SiPixelCluster& cluster = *tkhit->cluster_pixel();
  //               dxreccluster.push_back(cluster.x() - preciseHit->localPosition().x());
  //               dyreccluster.push_back(cluster.y() - preciseHit->localPosition().y());
  //             }
  //             else {
  //               dxreccluster.push_back(-99.);
  //               dyreccluster.push_back(-99.);
  //             }
              
              if (doSim_) {
                if (simhit != nullptr) {
                  
//                   std::cout << "ihit = " << ihit << " eloss = " << simhit->energyLoss() << std::endl;
                  
                  Vector2d dy0simgenlocal;
                  dy0simgenlocal << simhit->localPosition().x() - localxval,
                                  simhit->localPosition().y() - localyval;
                  const Vector2d dysimgeneig = R*dy0simgenlocal;
                  dxsimgen.push_back(dysimgeneig[0]);
                  dysimgen.push_back(dysimgeneig[1]);
  //                     dxsimgen.push_back(simhit->localPosition().x() - updtsos.localPosition().x());
  //                     dysimgen.push_back(simhit->localPosition().y() - updtsos.localPosition().y());
                  
                  Vector2d dy0simgenlocalconv;
                  dy0simgenlocalconv << simhit->localPosition().x() - localxval - localconv[0],
                                  simhit->localPosition().y() - localyval - localconv[1];
                  const Vector2d dysimgeneigconv = R*dy0simgenlocalconv;
                  dxsimgenconv.push_back(dysimgeneigconv[0]);
                  dysimgenconv.push_back(dysimgeneigconv[1]);
                  
                  dxsimgenlocal.push_back(dy0simgenlocal[0]);
                  dysimgenlocal.push_back(dy0simgenlocal[1]);
                  
                  Vector2d dyrecsimlocal;
                  dyrecsimlocal << preciseHit->localPosition().x() - simhit->localPosition().x(),
                                  preciseHit->localPosition().y() - simhit->localPosition().y();
                  const Vector2d dyrecsimeig = R*dyrecsimlocal;
                  dxrecsim.push_back(dyrecsimeig[0]);
                  dyrecsim.push_back(dyrecsimeig[1]);
                  
                  dE.push_back(dEval);
                  
                  auto const momsim = simhit->momentumAtEntry();

                  const double eentry = std::sqrt(std::pow(simhit->pabs(), 2) + mmu*mmu);
                  const double emid = eentry - 0.5*simhit->energyLoss();
                  const double simqopval = genpart->charge()/std::sqrt(emid*emid - mmu*mmu);
//                   std::cout << "eloss = " << simhit->energyLoss() << std::endl;
                  
                  simlocalqop.push_back(simqopval);
                  simlocaldxdz.push_back(momsim.x()/momsim.z());
                  simlocaldydz.push_back(momsim.y()/momsim.z());
                  simlocalx.push_back(simhit->localPosition().x());
                  simlocaly.push_back(simhit->localPosition().y());
                  
                }
                else {
                  dxsimgen.push_back(-99.);
                  dysimgen.push_back(-99.);
                  dxsimgenconv.push_back(-99.);
                  dysimgenconv.push_back(-99.);
                  dxsimgenlocal.push_back(-99.);
                  dysimgenlocal.push_back(-99.);
                  dxrecsim.push_back(-99.);
                  dyrecsim.push_back(-99.);
                  dE.push_back(-99.);
                  
                  simlocalqop.push_back(-99.);
                  simlocaldxdz.push_back(-99.);
                  simlocaldydz.push_back(-99.);
                  simlocalx.push_back(-99.);
                  simlocaly.push_back(-99.);
                }
              }

            }
            
          };
                    
//           if (align2d) {
//             fillAlignGrads(std::integral_constant<unsigned int, 3>());
//           }
//           else {
//             fillAlignGrads(std::integral_constant<unsigned int, 2>());
//           }
          if (align2d) {
            fillAlignGrads(std::integral_constant<unsigned int, 6>());
          }
          else {
            fillAlignGrads(std::integral_constant<unsigned int, 5>());
          }
//           fillAlignGrads(std::integral_constant<unsigned int, 6>());
          
          ivalidhit++;
            
        }
        
//         std::cout << "hit " << ihit << " isvalid " << preciseHit->isValid() << std::endl;
//         std::cout << "global position: " << updtsos.globalParameters().position() << std::endl;
        //hit information
        //FIXME consolidate this special cases into templated function(s)
//         if (preciseHit->isValid()) {
      }
            
      if (!valid) {
        break;
      }
      
      assert(parmidx == (nparsBfield + nparsEloss));
      assert(alignmentparmidx == nparsAlignment);
      

      auto freezeparm = [&](unsigned int idx) {
        gradfull[idx] = 0.;
        hessfull.row(idx) *= 0.;
        hessfull.col(idx) *= 0.;
        hessfull(idx,idx) = 1e6;
      };

      //fake constraint on reference point parameters
      if (dogen) {
        for (unsigned int i=0; i<5; ++i) {
          freezeparm(i);
        }
      }

      if (fitFromSimParms_) {
        for (unsigned int i = 5; i < nstateparms; ++i) {
          freezeparm(i);
        }
      }
      

//       // brute force constraint on reference point parameters
//       if (dogen) {
//         const double sigmaqop = 1e-6*1./genpart->p();
//         hessfull(0,0) += 1./sigmaqop/sigmaqop;
//         hessfull(1,1) += 1./1e-8/1e-8;
//         hessfull(2,2) += 1./1e-8/1e-8;
//         hessfull(3,3) += 1./1e-7/1e-7;
//         hessfull(4,4) += 1./1e-7/1e-7;
//       }
      

      
//       std::cout << "dchisqdx" << std::endl;
//       std::cout << dchisqdx << std::endl;
//       std::cout << "d2chisqdx2 diagonal" << std::endl;
//       std::cout << d2chisqdx2.diagonal() << std::endl;
//       std::cout << "d2chisqdx2" << std::endl;
//       std::cout << d2chisqdx2 << std::endl;
//       
//       auto const& eigenvalues = d2chisqdx2.eigenvalues();
//       std::cout << "iiter = " << iiter << std::endl;
//       std::cout << "d2chisqdx2 eigenvalues" << std::endl;
//       std::cout << eigenvalues << std::endl;
//       std::cout << "condition = " << eigenvalues[nstateparms-1]/eigenvalues[0] << std::endl;
      
//       auto const& Cinvd = d2chisqdx2.ldlt();
      
      
      //now do the expensive calculations and fill outputs
      auto const& dchisqdx = gradfull.head(nstateparms);
      auto const& dchisqdparms = gradfull.tail(npars);
      
      auto const& d2chisqdx2 = hessfull.topLeftCorner(nstateparms, nstateparms);
      auto const& d2chisqdxdparms = hessfull.topRightCorner(nstateparms, npars);
      auto const& d2chisqdparms2 = hessfull.bottomRightCorner(npars, npars);
      
      Cinvd.compute(d2chisqdx2);
      

      
      dxfull = -Cinvd.solve(dchisqdx);
      
      
//       const Matrix<double, 1, 1> deltachisq = dchisqdx.transpose()*dxfull + 0.5*dxfull.transpose()*d2chisqdx2*dxfull;
      
      const Matrix<double, 1, 1> deltachisq = 0.5*dchisqdx.transpose()*dxfull;
      
//       if (islikelihood) {
//         dxfull *= 0.5;
//       }
      
//       const Matrix<double, 1, 1> deltachisqalt = dchisqdx.transpose()*dxfull + 0.5*dxfull.transpose()*d2chisqdx2*dxfull;
// //       
//       const VectorXd deltachisqv = 0.5*dchisqdx.array()*dxfull.array();
// //       const VectorXd deltachisqvabs = deltachisqv.array().abs();
// // //       const unsigned int maxidx = deltachisqvabs.maxCoeff();
// //       
//       unsigned int maxidx = 0;
//       double maxval = 0.;
//       for (unsigned int i=0; i < deltachisqv.size(); ++i) {
//         const double absval = std::abs(deltachisqv[i]);
//         if (absval > maxval) {
//           maxval = absval;
//           maxidx = i;
//         }
//       }
// //       
//       std::cout << "deltachisq = " << deltachisq[0] << " deltachisqalt = " << deltachisqalt << " size = " << deltachisqv.size() <<  " maxidx = " << maxidx << " maxval = " << maxval << std::endl;
// //       
// //       
//       std::cout << "i dchisqdx dxfull deltachisqv" << std::endl;
//       for (unsigned int i=0; i < deltachisqv.size(); ++i) {
//         std::cout << i << " " << dchisqdx[i] << " " << dxfull[i] << " " << deltachisqv[i] << std::endl;
//       }

      
//       std::cout << "deltachisqv" << std::endl;
//       std::cout << deltachisqv << std::endl;
      
      
      
      chisqval = chisq0val + deltachisq[0];
      
      deltachisqval = chisq0val + deltachisq[0] - chisqvalold;
      
      chisqvalold = chisq0val + deltachisq[0];
        
//       ndof = 5*nhits + nvalid + nvalidalign2d - nstateparms;
      ndof = 5*nhits + nvalid + nvalidpixel - nstateparms;
//       ndof = 5*nhits + 2.*nvalid - nstateparms;
      
      if (bsConstraint_) {
        ndof += 2;
      }
      
      if (dogen) {
        ndof += 5;
      }
      
      if (fitFromSimParms_) {
        ndof += nstateparms - 5;
      }
      
      covfull = 2.*Cinvd.solve(MatrixXd::Identity(nstateparms,nstateparms));
      

      
//       std::cout << "iiter = " << iiter << ", deltachisq = " << deltachisq[0] << std::endl;
      
//       const Vector5d dxRef = dx.head<5>();
// //       const Vector5d dxRef = -Cinvd.solve(dchisqdx).head<5>();
//       const Matrix5d Cinner = Cinvd.solve(MatrixXd::Identity(nstateparms,nstateparms)).topLeftCorner<5,5>();
      
//       dxdparms = -Cinvd.solve(d2chisqdxdparms).transpose();
//       
//       grad = dchisqdparms + dxdparms*dchisqdx;
//       hess = d2chisqdparms2 + 2.*dxdparms*d2chisqdxdparms + dxdparms*d2chisqdx2*dxdparms.transpose();
//       
//       std::cout << "dxfull" << std::endl;
//       std::cout << dxfull << std::endl;
//       std::cout << "errsq" << std::endl;
//       std::cout << Cinvd.solve(MatrixXd::Identity(nstateparms,nstateparms)).diagonal() << std::endl;
      
      const Vector5d dxRef = dxfull.head<5>();
//       const Matrix5d Cinner = 2.*(Cinvd.solve(MatrixXd::Identity(nstateparms,nstateparms))).topLeftCorner<5,5>();
      const Matrix5d Cinner = covfull.topLeftCorner<5,5>();
      
//       const Matrix5d Cinner = Cinvd.solve(MatrixXd::Identity(nstateparms,nstateparms)).topLeftCorner<5,5>();

      if (debugprintout_) {
        std::cout<< "dxRef" << std::endl;
        std::cout<< dxRef << std::endl;
      }
      
      const Matrix<double, 5, 5> hessref = Cinner.inverse();
      
      const double deltachisqref = -0.5*dxRef.transpose()*hessref*dxRef;
      
      edmvalref = -deltachisqref;
      
      
      //fill output with corrected state and covariance at reference point
      refParms.fill(0.);
      refCov.fill(0.);
//       const AlgebraicVector5& refVec = track.parameters();
//       CurvilinearTrajectoryParameters curvparms(refFts.position(), refFts.momentum(), refFts.charge());
//       const AlgebraicVector5& refVec = curvparms.vector();
//       Map<Vector5f>(refParms.data()) = (Map<const Vector5d>(refVec.Array()) + dxRef).cast<float>();
      
      const double qbp = refFts[6]/refFts.segment<3>(3).norm();
      const double lam = std::atan(refFts[5]/std::sqrt(refFts[3]*refFts[3] + refFts[4]*refFts[4]));
      const double phi = std::atan2(refFts[4], refFts[3]);
      
      const double qbpupd = qbp + dxRef[0];
      const double lamupd = lam + dxRef[1];
      const double phiupd = phi + dxRef[2];
      
      refParms[0] = qbpupd;
      refParms[1] = lamupd;
      refParms[2] = phiupd;
      //TODO (longstanding) fix filling of position parameters

      refParmsMomD[0] = qbpupd;
      refParmsMomD[1] = lamupd;
      refParmsMomD[2] = phiupd;
      
      Map<Matrix<float, 5, 5, RowMajor> >(refCov.data()).triangularView<Upper>() = (Cinner).cast<float>().triangularView<Upper>();
      
      if (iiter==0) {
        refParms_iter0 = refParms;
        refCov_iter0 = refCov;
      }

//       // fill values for valuemap
//       if (muonref.isNonnull()) {
//         const double ptupd = std::cos(lamupd)/std::abs(qbpupd);
//         const double chargeupd = std::copysign(1., qbpupd);
//
//         const double thetaupd = M_PI_2 - lamupd;
//         const double etaupd = -std::log(std::tan(0.5*thetaupd));
//
//         corPtV[muonref.index()] = ptupd;
//         corEtaV[muonref.index()] = etaupd;
//         corPhiV[muonref.index()] = phiupd;
//         corChargeV[muonref.index()] = chargeupd;
//       }

      
      
         
      //current state from previous state on this layer
      //save current parameters          
//       TrajectoryStateOnSurface& oldtsosout = layerStates.back();
//       
//       if (iiter == 0) {
//         outPtStart = oldtsosout.globalMomentum().perp();
//         outEtaStart = oldtsosout.globalMomentum().eta();
//         outPhiStart = oldtsosout.globalMomentum().phi();
//       }
//       
// //           JacobianCurvilinearToLocal curv2localold(oldtsos.surface(), oldtsos.localParameters(), *oldtsos.magneticField());
// //           const AlgebraicMatrix55& curv2localjacold = curv2localold.jacobian();
// //           const Matrix<double, 5, 5> Hold = Map<const Matrix<double, 5, 5, RowMajor>>(curv2localjacold.Array()); 
// //           const Matrix<double, 5, 5> Hold = curv2localJacobianAlt(oldtsos);
//       const Matrix<double, 5, 5> Hold = curv2localJacobianAlteloss(oldtsosout, dEdxlast, mmu);
//       
//       const AlgebraicVector5 local = oldtsosout.localParameters().vector();
// //           auto const& dxlocal = dxfull.segment<5>(5*(ihit+1));
//       auto const& dxlocal = Hold*dxfull.tail<5>();
//       const Matrix<double, 5, 1> localupd = Map<const Matrix<double, 5, 1>>(local.Array()) + dxlocal;
//       AlgebraicVector5 localvecupd(localupd[0],localupd[1],localupd[2],localupd[3],localupd[4]);
//       
//       
//       const LocalTrajectoryParameters localparms(localvecupd, oldtsosout.localParameters().pzSign());
//       
//       const TrajectoryStateOnSurface updtsosout(localparms, oldtsosout.surface(), field, oldtsosout.surfaceSide());
//       
//       outPt = updtsosout.globalMomentum().perp();
//       outEta = updtsosout.globalMomentum().eta();
//       outPhi = updtsosout.globalMomentum().phi();

        
      
      niter = iiter + 1;
      edmval = -deltachisq[0];
      
      const double threshparam = dolocalupdate ? edmval : std::fabs(deltachisqval);

//       if (iiter>0 && threshparam > 1e4) {
//         anomDebug = true;
//       }

//       if ( !iscosmic && (std::isnan(edmval) || std::isinf(edmval) || std::abs(lamupd) > M_PI_2 || (iiter>0 && threshparam > 1e5) || (iiter>1 && threshparam > 1e4) )) {
      if (std::isnan(edmval) || std::isinf(edmval)) {
        std::cout << "WARNING: invalid parameter update!!!" << " edmval = " << edmval << " lamupd = " << lamupd << " deltachisqval = " << deltachisqval << std::endl;
        valid = false;
        break;
      }
      

      
//       const double nll0val = 0.5*chisq0val - 0.5*logdetH;
//       const double nllval = nll0val + deltachisq[0];
      
      const double ptupd = (1./std::abs(qbpupd))*std::cos(lamupd);

//       std::cout << "iiter = " << iiter << " edmval = " << edmval << " chisqval = " << chisqval << " chisq0val = " << chisq0val << " qbpupd = " << qbpupd << " lamupd = " << lamupd << " phiupd = " << phiupd << " ptupd = " << ptupd << std::endl;
      
      
//       std::cout << "iiter = " << iiter << " edmval = " << edmval << " chisqval = " << chisqval << " chisq0val = " << chisq0val << " logdetH = " << logdetH << " nll0val = " << nll0val << " nllval = " << nllval << std::endl;

      
      
//       std::cout << "iiter = " << iiter << " edmval = " << edmval << " chisqval = " << chisqval << std::endl;

//       std::cout << "iiter = " << iiter << " edmval = " << edmval << " deltachisqval = " << deltachisqval << " chisqval = " << chisqval << std::endl;
      
//       std::cout << "iiter = " << iiter << " edmval = " << edmval << " edmvalref " << edmvalref << " deltachisqval = " << deltachisqval << " chisqval = " << chisqval << std::endl;

      if (anomDebug) {
        std::cout << "anomDebug: iiter = " << iiter << " edmval = " << edmval << " deltachisqval = " << deltachisqval << " chisqval = " << chisqval << std::endl;
      }

//       std::cout <<"dxRef" << std::endl;
//       std::cout << dxRef << std::endl;
//       std::cout <<"dxOut" << std::endl;
//       std::cout << dxfull.tail<5>() << std::endl;
      
//       std::cout << "dxfull" << std::endl;
//       std::cout << dxfull << std::endl;
      
//       if (iiter > 1 && std::abs(deltachisq[0])<1e-3) {
//         break;
//       }
      
//       if (std::abs(deltachisqval)<1e-2) {
//         break;
//       }
      
      if (iiter > 0 && dolocalupdate && edmval < 1e-5) {
        break;
      }
//       else if (iiter > 0 && !dolocalupdate && std::fabs(deltachisqval)<1e-5) {
      else if (iiter > 0 && !dolocalupdate && edmvalref < 1e-5) {
        break;
      }
//       else if (iiter==2) {
//         refParms_iter2 = refParms;
//         refCov_iter2 = refCov;        
//       }
      
//       std::cout << "refParms" << std::endl;
//       std::cout << Map<const Vector5f>(refParms.data()) << std::endl;
      
//   //     gradv.clear();
//       jacrefv.clear();
// 
//   //     gradv.resize(npars,0.);
//       jacrefv.resize(5*npars, 0.);
//       
//       nJacRef = 5*npars;
//   //     tree->SetBranchAddress("gradv", gradv.data());
//       tree->SetBranchAddress("jacrefv", jacrefv.data());
//       
//       //eigen representation of the underlying vector storage
//   //     Map<VectorXf> gradout(gradv.data(), npars);
//       Map<Matrix<float, 5, Dynamic, RowMajor> > jacrefout(jacrefv.data(), 5, npars);
//       
//       jacrefout = dxdparms.leftCols<5>().transpose().cast<float>();
    
    }
    
    if (!valid) {
      continue;
    }
    
    if (!nValidHitsFinal) {
      continue;
    }
    
//     if (islikelihood) {
    if (false) {
      for (unsigned int ihit = 0; ihit < nhits; ++ihit) {
        
        if (dogen) {
  //       if (false) {
          for (unsigned int i=0; i<5; ++i) {
            dhessv[ihit].row(i) *= 0.;
            dhessv[ihit].col(i) *= 0.;
          }
        }
        
        dhessv[ihit] = (Cinvd.solve(dhessv[ihit])).eval();
//         const double sumfull = dhessv[ihit].array().abs().sum();
//         const double sumblock = dhessv[ihit].block<10, 10>(5*ihit, 5*ihit).array().abs().sum();
//         const double diffsum = sumfull - sumblock;
//         std::cout << "ihit = " << ihit << " diffsum:" << diffsum << std::endl;
//         std::cout << "ihit = " << ihit << " dhessv:" << std::endl;
//         std::cout << dhessv[ihit] << std::endl;
//         Matrix<int, Dynamic, Dynamic> testmat = Matrix<int, Dynamic, Dynamic>::Zero(nstateparms, nstateparms);
//         for (unsigned int j = 0; j < nstateparms; ++j) {
//           for (unsigned int k = 0; k < nstateparms; ++k) {
//             if (std::abs(dhessv[ihit](j,k)) > 0.) {
//               testmat(j,k) = 1;
//             }
//           }
//         }
//         std::cout << "ihit = " << ihit << " testmat:" << std::endl;
//         std::cout << testmat << std::endl;
        
      }
      
      for (unsigned int ihit = 0; ihit < nhits; ++ihit) {
        const double gradval = -dhessv[ihit].trace();
        const unsigned int ifull = nstateparms + 2*ihit + 1;
        gradfull[ifull] += gradval;
        std::cout << "ihit = " << ihit << " gradval = " << gradval << std::endl;
        for (unsigned int jhit = ihit; jhit < nhits; ++jhit) {
          const unsigned int jfull = nstateparms + 2*jhit + 1;
          const double hessval = dhessv[ihit].transpose().cwiseProduct(dhessv[jhit]).sum();
          std::cout << "ihit = " << ihit << " jhit = " << jhit << " hessval = " << hessval << std::endl;
          
          hessfull(ifull, jfull) += hessval;
          if (ihit != jhit) {
            hessfull(jfull, ifull) += hessval;
          }
        }
      }
      
    }
    
    
    auto const& dchisqdx = gradfull.head(nstateparms);
    auto const& dchisqdparms = gradfull.tail(npars);
    
    auto const& d2chisqdx2 = hessfull.topLeftCorner(nstateparms, nstateparms);
    auto const& d2chisqdxdparms = hessfull.topRightCorner(nstateparms, npars);
    auto const& d2chisqdparms2 = hessfull.bottomRightCorner(npars, npars);
    
    std::unordered_map<unsigned int, unsigned int> idxmap;
    
    globalidxvfinal.clear();
    globalidxvfinal.reserve(npars);
    idxmap.reserve(npars);
    
    for (unsigned int idx : globalidxv) {
      if (!idxmap.count(idx)) {
        idxmap[idx] = globalidxvfinal.size();
        globalidxvfinal.push_back(idx);
      }
    }
    
    const unsigned int nparsfinal = globalidxvfinal.size();
    
    VectorXd dchisqdparmsfinal = VectorXd::Zero(nparsfinal);
    MatrixXd d2chisqdxdparmsfinal = MatrixXd::Zero(nstateparms, nparsfinal);
    MatrixXd d2chisqdparms2final = MatrixXd::Zero(nparsfinal, nparsfinal);
    
    for (unsigned int i = 0; i < npars; ++i) {
      const unsigned int iidx = idxmap.at(globalidxv[i]);
      dchisqdparmsfinal[iidx] += dchisqdparms[i];
      d2chisqdxdparmsfinal.col(iidx) += d2chisqdxdparms.col(i);
      for (unsigned int j = 0; j < npars; ++j) {
        const unsigned int jidx = idxmap.at(globalidxv[j]);
        d2chisqdparms2final(iidx, jidx) += d2chisqdparms2(i, j);
      }
    }
    
    dxdparms = -Cinvd.solve(d2chisqdxdparmsfinal).transpose();
    
//     grad = dchisqdparmsfinal + dxdparms*dchisqdx; 
    grad = dchisqdparmsfinal + d2chisqdxdparmsfinal.transpose()*dxfull;
    hess = d2chisqdparms2final + dxdparms*d2chisqdxdparmsfinal;


//     std::cout <<"dchisqdparmsfinal:\n" << dchisqdparmsfinal << std::endl;
//     std::cout << "grad\n" << grad << std::endl;
//     std::cout << "dxfull\n" << dxfull << std::endl;

//     std::cout << "gradalt\n" << gradalt << std::endl;

    
//     dxdparms = -Cinvd.solve(d2chisqdxdparms).transpose();
    
    
//     if (islikelihood) {
//     
//       for (unsigned int ihit = 0; ihit < nhits; ++ihit) {
//         const unsigned int iparmfull = 2*ihit + 1;
//         dxdparms.row(iparmfull) += -(dhessv[ihit]*dxfull).transpose();
//       }
//     
//     }
    
//     if (debugprintout_) {
//       std::cout << "dxrefdparms" << std::endl;
//       std::cout << dxdparms.leftCols<5>() << std::endl;
//     }
    
//     grad = dchisqdparms + dxdparms*dchisqdx;
    //TODO check the simplification
//     hess = d2chisqdparms2 + 2.*dxdparms*d2chisqdxdparms + dxdparms*d2chisqdx2*dxdparms.transpose();
//     hess = d2chisqdparms2 + dxdparms*d2chisqdxdparms;
    
//     std::cout << "grad:" << std::endl;
//     std::cout << grad.transpose() << std::endl;
//     std::cout << "hess diagonal:" << std::endl;
//     std::cout << hess.diagonal().transpose() << std::endl;
    
//     const Vector5d dxRef = dxfull.head<5>();
//     const Matrix5d Cinner = Cinvd.solve(MatrixXd::Identity(nstateparms,nstateparms)).topLeftCorner<5,5>();

//     std::cout << "dchisqdparms.head<6>()" << std::endl;
//     std::cout << dchisqdparms.head<6>() << std::endl;
//     
//     std::cout << "grad.head<6>()" << std::endl;
//     std::cout << grad.head<6>() << std::endl;
//     
//     std::cout << "d2chisqdparms2.topLeftCorner<6, 6>():" << std::endl;
//     std::cout << d2chisqdparms2.topLeftCorner<6, 6>() << std::endl;
//     std::cout << "hess.topLeftCorner<6, 6>():" << std::endl;
//     std::cout << hess.topLeftCorner<6, 6>() << std::endl;
//     
//     std::cout << "dchisqdparms.segment<6>(nparsBfield+nparsEloss)" << std::endl;
//     std::cout << dchisqdparms.segment<6>(nparsBfield+nparsEloss) << std::endl;
//     
//     std::cout << "grad.segment<6>(nparsBfield+nparsEloss)" << std::endl;
//     std::cout << grad.segment<6>(nparsBfield+nparsEloss) << std::endl;
//     
//     std::cout << "d2chisqdparms2.block<6, 6>(nparsBfield+nparsEloss, nparsBfield+nparsEloss):" << std::endl;
//     std::cout << d2chisqdparms2.block<6, 6>(nparsBfield+nparsEloss, nparsBfield+nparsEloss) << std::endl;
//     std::cout << "hess.block<6, 6>(nparsBfield+nparsEloss, nparsBfield+nparsEloss):" << std::endl;
//     std::cout << hess.block<6, 6>(nparsBfield+nparsEloss, nparsBfield+nparsEloss) << std::endl;
//     
//     std::cout << "d2chisqdparms2.bottomRightCorner<6, 6>():" << std::endl;
//     std::cout << d2chisqdparms2.bottomRightCorner<6, 6>() << std::endl;
//     std::cout << "hess.bottomRightCorner<6, 6>():" << std::endl;
//     std::cout << hess.bottomRightCorner<6, 6>() << std::endl;
    
    nParms = nparsfinal;

    
    gradv.clear();
    jacrefv.clear();

    gradv.resize(nparsfinal,0.);
    jacrefv.resize(5*nparsfinal, 0.);
    
    nJacRef = 5*nparsfinal;
    if (fillTrackTree_ && fillGrads_) {
      tree->SetBranchAddress("gradv", gradv.data());
    }
    if (fillTrackTree_ && fillJac_) {
      tree->SetBranchAddress("jacrefv", jacrefv.data());
    }
    
    //eigen representation of the underlying vector storage
    Map<VectorXf> gradout(gradv.data(), nparsfinal);
    Map<Matrix<float, 5, Dynamic, RowMajor> > jacrefout(jacrefv.data(), 5, nparsfinal);
    
//     jacrefout = dxdparms.leftCols<5>().transpose().cast<float>();    
    jacrefout = ( (dxdparms).leftCols<5>().transpose() ).cast<float>();  
    

    gradout = grad.cast<float>();
    
    
    rx.resize(2*nvalid);
    Map<Matrix<float, Dynamic, 2, RowMajor> > rxout(rx.data(), nvalid, 2);
    rxout = rxfull;
//     std::cout << "rx:" << std::endl;
//     for (auto elem : rx) {
//       std::cout << elem << " ";
//     }
//     std::cout << std::endl;
//     std::cout << rx << std::endl;
    
    ry.resize(2*nvalid);
    Map<Matrix<float, Dynamic, 2, RowMajor> > ryout(ry.data(), nvalid, 2);
    ryout = ryfull;
    
    deigx.resize(nvalid);
    deigy.resize(nvalid);
    
    validdxeig = validdxeigjac*dxfull;
    
    for (unsigned int ivalid = 0; ivalid < nvalid; ++ivalid) {
      deigx[ivalid] = validdxeig[2*ivalid];
      deigy[ivalid] = validdxeig[2*ivalid + 1];
    }
    
//     unsigned int ivalid = 0;
//     for (unsigned int ihit = 0; ihit < nhits; ++ihit) {
//       auto const& hit = hits[ihit];
//       if (hit->isValid()) {
// //         if (ihit < (nhits-1)) {
// //           std::cout << "ihit = " << ihit << ", ivalid = " << ivalid << std::endl;
// //           std::cout << "dxfull.segment<3>(3*(ihit+1)):" << std::endl;
// //           std::cout << dxfull.segment<3>(3*(ihit+1)) << std::endl;
// //           std::cout << "dxstate.segment<5>(5*(ihit+1))" << std::endl;
// //           std::cout << dxstate.segment<5>(5*(ihit+1)) << std::endl;
// //         }
//         const unsigned int dxidx = 3*(ihit + 1);
//         dlocalx[ivalid] = dxfull[dxidx];
//         dlocaly[ivalid] = dxfull[dxidx + 1];
//         ++ivalid;
//       }
//     }
    
    
    float refPt = dogen ? genpart->pt() : std::abs(1./refParms[0])*std::sin(M_PI_2 - refParms[1]);

    gradmax = 0.;
    for (unsigned int i=0; i<nparsfinal; ++i) {
      const float absval = std::abs(grad[i]);
      if (absval>gradmax) {
        gradmax = absval;
      }      
    }
    
    
    hessmax = 0.;
    for (unsigned int i=0; i<nparsfinal; ++i) {
      for (unsigned int j=i; j<nparsfinal; ++j) {
        const unsigned int iidx = globalidxvfinal[i];
        const unsigned int jidx = globalidxvfinal[j];
        
        const float absval = std::abs(hess(i,j));
        if (absval>hessmax) {
          hessmax = absval;
        }
        
      }
      
    }
    
//     if (gradmax < 1e5 && refPt > 5.5) {
//       //fill aggregrate gradient and hessian
//       for (unsigned int i=0; i<nparsfinal; ++i) {
//         gradagg[globalidxvfinal[i]] += grad[i];
//       }
//       
//       hessmax = 0.;
//       for (unsigned int i=0; i<nparsfinal; ++i) {
//         for (unsigned int j=i; j<nparsfinal; ++j) {
//           const unsigned int iidx = globalidxvfinal[i];
//           const unsigned int jidx = globalidxvfinal[j];
//           
//           const float absval = std::abs(hess(i,j));
//           if (absval>hessmax) {
//             hessmax = absval;
//           }
//           
//           const std::pair<unsigned int, unsigned int> key = std::make_pair(std::min(iidx,jidx), std::max(iidx,jidx));
//           
//           auto it = hessaggsparse.find(key);
//           if (it==hessaggsparse.end()) {
//             hessaggsparse[key] = hess(i,j);
//           }
//           else {
//             it->second += hess(i,j);
//           }
//         }
//       }
//     }
    
    
    if (debugprintout_) {
      const Matrix5d Cinner = (Cinvd.solve(MatrixXd::Identity(nstateparms,nstateparms))).topLeftCorner<5,5>();
      std::cout << "hess debug" << std::endl;
      std::cout << "track parms" << std::endl;
      std::cout << tkparms << std::endl;
  //     std::cout << "dxRef" << std::endl;
  //     std::cout << dxRef << std::endl;
      std::cout << "original cov" << std::endl;
      std::cout << track.covariance() << std::endl;
      std::cout << "recomputed cov" << std::endl;
      std::cout << 2.*Cinner << std::endl;
    }

//     std::cout << "dxinner/dparms" << std::endl;
//     std::cout << dxdparms.bottomRows<5>() << std::endl;
//     std::cout << "grad" << std::endl;
//     std::cout << grad << std::endl;
//     std::cout << "hess diagonal" << std::endl;
//     std::cout << hess.diagonal() << std::endl;
//     std::cout << "hess0 diagonal" << std::endl;
//     std::cout << d2chisqdparms2.diagonal() << std::endl;
//     std::cout << "hess1 diagonal" << std::endl;
//     std::cout << 2.*(dxdparms.transpose()*d2chisqdxdparms).diagonal() << std::endl;
//     std::cout << "hess2 diagonal" << std::endl;
//     std::cout << (dxdparms.transpose()*d2chisqdx2*dxdparms).diagonal() << std::endl;
    
    //fill packed hessian and indices
    const unsigned int nsym = nparsfinal*(1+nparsfinal)/2;
    hesspackedv.clear();    
    hesspackedv.resize(nsym, 0.);
    
    nSym = nsym;
    if (fillTrackTree_ && fillGrads_) {
      tree->SetBranchAddress("hesspackedv", hesspackedv.data());
    }
    
    Map<VectorXf> hesspacked(hesspackedv.data(), nsym);
    const Map<const VectorXu> globalidx(globalidxvfinal.data(), nparsfinal);

    unsigned int packedidx = 0;
    for (unsigned int ipar = 0; ipar < nparsfinal; ++ipar) {
      const unsigned int segmentsize = nparsfinal - ipar;
      hesspacked.segment(packedidx, segmentsize) = hess.block<1, Dynamic>(ipar, ipar, 1, segmentsize).cast<float>();
      packedidx += segmentsize;
    }

//     std::cout << "refParms[0]: " << refParms[0] << std::endl;

    if (fillTrackTree_) {
      tree->Fill();
    }



    if (muonref.isNonnull()) {
      const double qbpupd = refParmsMomD[0];
      const double lamupd = refParmsMomD[1];
      const double phiupd = refParmsMomD[2];

      const double ptupd = std::cos(lamupd)/std::abs(qbpupd);
      const double chargeupd = std::copysign(1., qbpupd);

      const double thetaupd = M_PI_2 - lamupd;
      const double etaupd = -std::log(std::tan(0.5*thetaupd));

      corPtV[muonref.index()] = ptupd;
      corEtaV[muonref.index()] = etaupd;
      corPhiV[muonref.index()] = phiupd;
      corChargeV[muonref.index()] = chargeupd;
      edmvalV[muonref.index()] = edmval;

      auto &iglobalidxv = globalidxsV[muonref.key()];
      iglobalidxv.clear();
      iglobalidxv.reserve(globalidxvfinal.size());
      for (auto val : globalidxvfinal) {
        iglobalidxv.push_back(val);
      }

      auto &ijacrefv = jacRefV[muonref.key()];
      ijacrefv.assign(3*nparsfinal, 0.);
      //eigen representation of the underlying vector storage
      Map<Matrix<float, 3, Dynamic, RowMajor> > momjacrefout(ijacrefv.data(), 3, nparsfinal);
      momjacrefout = ( (dxdparms).leftCols<3>().transpose() ).cast<float>();

      auto &imomCov = momCovV[muonref.key()];
      imomCov.assign(3*3, 0.);
      Map<Matrix<float, 3, 3, RowMajor>> momCovout(imomCov.data(), 3, 3);
      momCovout.triangularView<Upper>() = covfull.topLeftCorner<3,3>().cast<float>();
//       std::cout << "covfull.topLeftCorner<3,3>()" << std::endl;
//       std::cout << covfull.topLeftCorner<3,3>() << std::endl;
//       std::cout << "momcovout" << std::endl;
//       std::cout << momCovout << std::endl;
    }



  }

  edm::ValueMap<float> corPtMap;
  edm::ValueMap<float> corEtaMap;
  edm::ValueMap<float> corPhiMap;
  edm::ValueMap<int> corChargeMap;
  edm::ValueMap<float> edmvalMap;

  edm::ValueMap<std::vector<int>> globalidxsMap;

  edm::ValueMap<std::vector<float>> jacRefMap;
  edm::ValueMap<std::vector<float>> momCovMap;

  if (doMuonAssoc_) {
    edm::ValueMap<float>::Filler corPtMapFiller(corPtMap);
    corPtMapFiller.insert(muonAssoc->ref(), std::make_move_iterator(corPtV.begin()), std::make_move_iterator(corPtV.end()));
    corPtMapFiller.fill();

    edm::ValueMap<float>::Filler corEtaMapFiller(corEtaMap);
    corEtaMapFiller.insert(muonAssoc->ref(), std::make_move_iterator(corEtaV.begin()), std::make_move_iterator(corEtaV.end()));
    corEtaMapFiller.fill();

    edm::ValueMap<float>::Filler corPhiMapFiller(corPhiMap);
    corPhiMapFiller.insert(muonAssoc->ref(), std::make_move_iterator(corPhiV.begin()), std::make_move_iterator(corPhiV.end()));
    corPhiMapFiller.fill();

    edm::ValueMap<int>::Filler corChargeMapFiller(corChargeMap);
    corChargeMapFiller.insert(muonAssoc->ref(), std::make_move_iterator(corChargeV.begin()), std::make_move_iterator(corChargeV.end()));
    corChargeMapFiller.fill();

    edm::ValueMap<float>::Filler edmvalMapFiller(edmvalMap);
    edmvalMapFiller.insert(muonAssoc->ref(), std::make_move_iterator(edmvalV.begin()), std::make_move_iterator(edmvalV.end()));
    edmvalMapFiller.fill();

    edm::ValueMap<std::vector<int>>::Filler globalidxsMapFiller(globalidxsMap);
    globalidxsMapFiller.insert(muonAssoc->ref(), std::make_move_iterator(globalidxsV.begin()), std::make_move_iterator(globalidxsV.end()));
    globalidxsMapFiller.fill();

    edm::ValueMap<std::vector<float>>::Filler jacRefMapFiller(jacRefMap);
    jacRefMapFiller.insert(muonAssoc->ref(), std::make_move_iterator(jacRefV.begin()), std::make_move_iterator(jacRefV.end()));
    jacRefMapFiller.fill();

    edm::ValueMap<std::vector<float>>::Filler momCovMapFiller(momCovMap);
    momCovMapFiller.insert(muonAssoc->ref(), std::make_move_iterator(momCovV.begin()), std::make_move_iterator(momCovV.end()));
    momCovMapFiller.fill();
  }

  iEvent.emplace(outputCorPt_, std::move(corPtMap));
  iEvent.emplace(outputCorEta_, std::move(corEtaMap));
  iEvent.emplace(outputCorPhi_, std::move(corPhiMap));
  iEvent.emplace(outputCorCharge_, std::move(corChargeMap));
  iEvent.emplace(outputEdmval_, std::move(edmvalMap));

  iEvent.emplace(outputGlobalIdxs_, std::move(globalidxsMap));

  iEvent.emplace(outputJacRef_, std::move(jacRefMap));
  iEvent.emplace(outputMomCov_, std::move(momCovMap));

}

DEFINE_FWK_MODULE(ResidualGlobalCorrectionMakerG4e);
