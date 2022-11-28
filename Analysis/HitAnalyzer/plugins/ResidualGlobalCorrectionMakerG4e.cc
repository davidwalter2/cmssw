#include "ResidualGlobalCorrectionMakerBase.h"
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
  
  bool trackHighPurity = false;

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

  float trackQopErr = 0.;

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
    tree->Branch("trackQopErr", &trackQopErr);
    //workaround for older ROOT version inability to store std::array automatically
  //   tree->Branch("trackOrigParms", trackOrigParms.data(), "trackOrigParms[5]/F", basketSize);
  //   tree->Branch("trackOrigCov", trackOrigCov.data(), "trackOrigCov[25]/F", basketSize);
    tree->Branch("trackParms", trackParms.data(), "trackParms[5]/F", basketSize);
    tree->Branch("trackCov", trackCov.data(), "trackCov[25]/F", basketSize);
    
    tree->Branch("trackHighPurity", &trackHighPurity);

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
      
      tree->Branch("simlocalqopprop", &simlocalqopprop);
      tree->Branch("simlocaldxdzprop", &simlocaldxdzprop);
      tree->Branch("simlocaldydzprop", &simlocaldydzprop);
      tree->Branch("simlocalxprop", &simlocalxprop);
      tree->Branch("simlocalyprop", &simlocalyprop);

      tree->Branch("landauDelta", &landauDelta);
      tree->Branch("landauW", &landauW);
      
      tree->Branch("localqop", &localqop);
      tree->Branch("localdxdz", &localdxdz);
      tree->Branch("localdydz", &localdydz);
      tree->Branch("localx", &localx);
      tree->Branch("localy", &localy);
      
      tree->Branch("localqoperr", &localqoperr);
      tree->Branch("localdxdzerr", &localdxdzerr);
      tree->Branch("localdydzerr", &localdydzerr);
      tree->Branch("localxerr", &localxerr);
      tree->Branch("localyerr", &localyerr);

      tree->Branch("hitlocalx", &hitlocalx);
      tree->Branch("hitlocaly", &hitlocaly);
      
      tree->Branch("localqop_iter", &localqop_iter);
      tree->Branch("localdxdz_iter", &localdxdz_iter);
      tree->Branch("localdydz_iter", &localdydz_iter);
      tree->Branch("localx_iter", &localx_iter);
      tree->Branch("localy_iter", &localy_iter);

      tree->Branch("dxrecgen_iter", &dxrecgen_iter);
      tree->Branch("dyrecgen_iter", &dyrecgen_iter);
      
      tree->Branch("localqoperralt", &localqoperralt);
      
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
  const bool dolocalupdate = fitFromSimParms_;
  // const bool dolocalupdate = true;
  
  const bool dores = doRes_;

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
    trackQopErr = track.qoverpError();

//     if (abs(trackEta) > 0.1) {
//       continue;
//     }

    trackHighPurity = track.quality(reco::TrackBase::highPurity);
    
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

//         std::cout << "splitting matched hit, inner = " << detinner->geographicalId().rawId() <<" outer = " << detouter->geographicalId().rawId() << std::endl;
        if (hits.size() > 0) {
//           std::cout << "previous = " << hits.back()->geographicalId().rawId() << std::endl;
          const bool duplicate = detinner->geographicalId() == hits.back()->geographicalId() || detouter->geographicalId() == hits.back()->geographicalId();

          if (duplicate) {
            std::cout << "WARNING: Duplicate hits from glued module splitting!\n";
          }
        }
        
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
//             hitquality = false;
          }
          else {
            assert(tkhit->cluster_strip().isNonnull());
            const SiStripCluster& cluster = *tkhit->cluster_strip();
            const StripTopology* striptopology = dynamic_cast<const StripTopology*>(&(detectorG->topology()));
            assert(striptopology);

            const uint16_t firstStrip = cluster.firstStrip();
            const uint16_t lastStrip = cluster.firstStrip() + cluster.amplitudes().size() - 1;
            const bool isOnEdge = firstStrip == 0 || lastStrip == (striptopology->nstrips() - 1);

            const bool isstereo = trackerTopology->isStereo((*it)->geographicalId());
            
            hitquality = true;
//             hitquality = false;
            // hitquality = !isstereo;
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

        const DetId aligndetid = alignGlued_ ? parmdetid : hit->geographicalId();

        
//         const bool align2d = detidparms.count(std::make_pair(1, hit->geographicalId()));
        const bool align2d = detidparms.count(std::make_pair(1, aligndetid));
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
//     const unsigned int nparsRes = nhits + nvalid + nvalidpixel;
    const unsigned int nparsRes = dores ? 2*nhits + nvalid + nvalidpixel : 0;
    const unsigned int npars = nparsAlignment + nparsBfield + nparsEloss + nparsRes;
    
    const unsigned int nstateparms = 5*(nhits+1);
    const unsigned int nparmsfull = nstateparms + npars;
    
    const int nbscons = bsConstraint_ ? 3 : 0;
    const int ncons = 5*nhits + nvalid + nvalidpixel + nbscons;
    
    using VectorXb = Matrix<bool, Dynamic, 1>;
    VectorXb freestatemask = VectorXb::Ones(nstateparms);
    
    if (dogen) {
      freestatemask.head<5>() = Matrix<bool, 5, 1>::Zero();
    }
    if (fitFromSimParms_) {
      freestatemask = VectorXb::Zero(nstateparms);
    }

    // if (dogen) {
    //   for (unsigned int istate = 0; istate < nstateparms; ++istate) {
    //     if (istate % 5 == 0) {
    //       freestatemask[istate] = false;
    //     }
    //   }
    // }
    
    std::vector<Eigen::Index> freestateidxs;
    freestateidxs.reserve(nstateparms);
    
    for (Eigen::Index istate = 0; istate < nstateparms; ++istate) {
      if (freestatemask[istate]) {
        freestateidxs.push_back(istate);
      }
    }
    
    const unsigned int nstatefree = freestateidxs.size();
    
    ndof = ncons - nstatefree;
    
    if (bsConstraint_) {
      ndof += 2;
    }
    
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
//     LDLT<MatrixXd> Cinvd;
    MatrixXd covfull = MatrixXd::Zero(nstateparms, nstateparms);
    
    VectorXd rfull;
    MatrixXd Ffull;
    MatrixXd Jfull;
    MatrixXd Vinvfull;
    VectorXd dxfree;
    
    MatrixXd Vinvfullalt;
    
    SparseMatrix<double> Fsparse;
    SparseMatrix<double> Vinvsparse;
    SparseMatrix<double> VinvF;
    
    SimplicialLDLT<SparseMatrix<double>> Cinvd;
    
    std::vector<SparseMatrix<double>> dVs;
    dVs.reserve(nhits + nvalid);
    
    std::vector<unsigned int> residxs;
    residxs.reserve(nhits + nvalid);
    
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
        
        simlocalqopprop.clear();
        simlocaldxdzprop.clear();
        simlocaldydzprop.clear();
        simlocalxprop.clear();
        simlocalyprop.clear();

        simlocalqopprop.reserve(nvalid);
        simlocaldxdzprop.reserve(nvalid);
        simlocaldydzprop.reserve(nvalid);
        simlocalxprop.reserve(nvalid);
        simlocalyprop.reserve(nvalid);
        
        landauDelta.clear();
        landauW.clear();

        landauDelta.reserve(nvalid);
        landauW.reserve(nvalid);

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
        
        localqoperr.clear();
        localdxdzerr.clear();
        localdydzerr.clear();
        localxerr.clear();
        localyerr.clear();

        localqoperr.reserve(nvalid);
        localdxdzerr.reserve(nvalid);
        localdydzerr.reserve(nvalid);
        localxerr.reserve(nvalid);
        localyerr.reserve(nvalid);
        
        hitlocalx.clear();
        hitlocaly.clear();
        
        hitlocalx.reserve(nvalid);
        hitlocaly.reserve(nvalid);
        
        localqop_iter.assign(nvalid, -99.);
        localdxdz_iter.assign(nvalid, -99.);
        localdydz_iter.assign(nvalid, -99.);
        localx_iter.assign(nvalid, -99.);
        localy_iter.assign(nvalid, -99.);

        dxrecgen_iter.assign(nvalid, -99.);
        dyrecgen_iter.assign(nvalid, -99.);
        
        localqoperralt.assign(nvalid, -99.);

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
      
      rfull = VectorXd::Zero(ncons);
      Ffull = MatrixXd::Zero(ncons, nstateparms);
      Jfull = MatrixXd::Zero(ncons, npars);
      Vinvfull = MatrixXd::Zero(ncons, ncons);
      Vinvfullalt = MatrixXd::Zero(ncons, ncons);
      
      dVs.clear();
      residxs.clear();
      
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
      
      unsigned int ivalidhit = 0;
      
      unsigned int icons = 0;
      unsigned int iparm = 0;

      if (iiter > 0) {
        //update current state from reference point state (errors not needed beyond first iteration)
        
        auto const& dxlocal = dxfull.head<5>();

        
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

      Matrix<double, 7, 1> propfromtsos = refFts;

      Matrix<double, 5, 5> Qtot = Matrix<double, 5, 5>::Zero();

      float e = genpart == nullptr ? -99. : std::sqrt(genpart->momentum().mag2() + mmu*mmu);
      float epred = std::sqrt(refFts.segment<3>(3).squaredNorm() + mmu*mmu);
      
      
      if (bsConstraint_) {
        // apply beamspot constraint
        // TODO add residual corrections for beamspot parameters?
        // TODO make this a 2D constraint?
        
        //TODO convert this to Ffull/Vinvfull parameterization
        
        constexpr unsigned int nlocalcons = 3;
        constexpr unsigned int nlocalstate = 5;
        constexpr unsigned int nlocal = nlocalstate;
        constexpr unsigned int localstateidx = 0;
        constexpr unsigned int fullstateidx = 0;
        
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

        // TODO check consistency of this parameterization

        // FIXME xy correlation is not stored and assumed to be zero
        const double corrb12 = 0.;
        
        const double varb1 = sigb1*sigb1;
        const double varb2 = sigb2*sigb2;
        const double varb3 = sigb3*sigb3;
        
        Matrix<double, nlocalcons, nlocalcons> covBS = Matrix<double, nlocalcons, nlocalcons>::Zero();
        // parametrisation: rotation (dx/dz, dy/dz); covxy
        covBS(0,0) = varb1;
        covBS(1,0) = covBS(0,1) = corrb12*sigb1*sigb2;
        covBS(1,1) = varb2;
        covBS(2,0) = covBS(0,2) = dxdz*(varb3-varb1)-dydz*covBS(1,0);
        covBS(2,1) = covBS(1,2) = dydz*(varb3-varb2)-dxdz*covBS(1,0);
        covBS(2,2) = varb3;
        
        Matrix<double, nlocalcons, 1> dbs0;
        dbs0[0] = refFts[0] - x0;
        dbs0[1] = refFts[1] - y0;
        dbs0[2] = refFts[2] - z0;
        
        const Matrix<double, nlocalcons, nlocalstate> Fbs = jac.topRows<3>();
        const Matrix<double, nlocalcons, nlocalcons> covBSinv = covBS.inverse();
        
        rfull.segment<nlocalcons>(icons) = dbs0;
        Ffull.block<nlocalcons, nlocalstate>(icons, fullstateidx) = Fbs;
        Vinvfull.block<nlocalcons, nlocalcons>(icons, icons) = covBSinv;
        Vinvfullalt.block<nlocalcons, nlocalcons>(icons, icons) = covBSinv;
        
        icons += nlocalcons;
      }
      
      
      for (unsigned int ihit = 0; ihit < hits.size(); ++ihit) {
//         std::cout << "iiter = " << iiter << " ihit " << ihit << std::endl;

        auto const& hit = hits[ihit];
        
        // auto const &surface = surfacemapD_.at(hit->geographicalId());
        GloballyPositioned<double> surface = surfacemapD_.at(hit->geographicalId());

        // check propagation direction to choose correct bfield and material parameters
        // when propagating inside-out the parameters correspond to the target module
        // when propagating outside-in (e.g. for cosmics) they correspond to the source module
        // For the first hit always use the target module
        // (only relevant for cosmics)

        bool sourceParms = false;
        if (iscosmic && ihit > 0) {
          auto const &lzhat = surface.rotation().z();
          auto const &pos = surface.position();

          const double zdotpos = lzhat.x()*pos.x() + lzhat.y()*pos.y() + lzhat.z()*pos.z();

          const Point3DBase<double, GlobalTag> globalpos(updtsos[0], updtsos[1], updtsos[2]);
          const Point3DBase<double, LocalTag> localpos = surface.toLocal(globalpos);

          sourceParms = zdotpos*localpos.z() > 0.;
        }

//         std::cout << "ihit = " << ihit << " sourceParms = " << sourceParms << std::endl;

        auto const& prophit = sourceParms ? hits[ihit - 1] : hit;
        const uint32_t gluedidprop = trackerTopology->glued(prophit->geographicalId());
        const bool isgluedprop = gluedidprop != 0;
        const DetId propdetid = isgluedprop ? DetId(gluedidprop) : prophit->geographicalId();

        const uint32_t gluedid = trackerTopology->glued(hit->geographicalId());
        const bool isglued = gluedid != 0;
        const DetId parmdetid = isglued ? DetId(gluedid) : hit->geographicalId();
        const DetId aligndetid = alignGlued_ ? parmdetid : hit->geographicalId();
        
        const unsigned int bfieldglobalidx = detidparms.at(std::make_pair(6, propdetid));
        const unsigned int elossglobalidx = detidparms.at(std::make_pair(7, propdetid));
        const unsigned int msglobalidx = dores ? detidparms.at(std::make_pair(10, propdetid)) : 0;
        const unsigned int ioniglobalidx = dores ? detidparms.at(std::make_pair(11, propdetid)) : 0;
        
        const double dbetaval = corparms_[bfieldglobalidx];
        const double dxival = corparms_[elossglobalidx];
        const double dmsval = dores ? corparms_[msglobalidx] : dxival;
        const double dionival = dores ? corparms_[ioniglobalidx] : dxival;
        
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

        // const bool simhitdebug = true;
        const bool simhitdebug = false;

        if (simhitdebug && !simhit) {
          std::cout << "WARNING: no simhit for debugging scenario, abort!\n";
          valid = false;
          break;
        }

        if ((simhitdebug || fitFromSimParms_) && simhit) {
          //move surface to plane of sim entry point for consistency in debugging vs sim state

          const double lz = simhit->entryPoint().z();
          const Vector3DBase<double, LocalTag> dlocal(0., 0., lz);
          const Vector3DBase<double, GlobalTag> dglobal = surface.toGlobal(dlocal);
          surface.move(dglobal);
        }

//         std::cout << "iiter = " << iiter << " ihit = " << ihit << " updtsos:\n" << updtsos << std::endl;
        
        // auto const &propresult = g4prop->propagateGenericWithJacobianAltD(updtsos, surface, dbetaval, dxival, dradval);
        // auto const &propresult = g4prop->propagateGenericWithJacobianAltD(propfromtsos, surface, dbetaval, dxival, dradval);
        
        auto const &propresult = simhitdebug ? g4prop->propagateGenericWithJacobianAltD(propfromtsos, surface, dbetaval, dxival, dmsval, dionival) : g4prop->propagateGenericWithJacobianAltD(updtsos, surface, dbetaval, dxival, dmsval, dionival);


        if (simhitdebug && simhit) {
          // reset propagation start point to sim hit for debugging purposes
          auto const &surfaceideal = surfacemapIdealD_.at(hit->geographicalId());
          const Point3DBase<double, GlobalTag> pos = surfaceideal.toGlobal(simhit->entryPoint());
          const Vector3DBase<double, GlobalTag> mom = surfaceideal.toGlobal(simhit->momentumAtEntry());

          propfromtsos[0] = pos.x();
          propfromtsos[1] = pos.y();
          propfromtsos[2] = pos.z();
          propfromtsos[3] = mom.x();
          propfromtsos[4] = mom.y();
          propfromtsos[5] = mom.z();
          propfromtsos[6] = genpart->charge();
        }

        if (!std::get<0>(propresult)) {
          std::cout << "Abort: Propagation Failed!" << std::endl;
          valid = false;
          break;
        }
        
        updtsos = std::get<1>(propresult);
        const Matrix<double, 5, 5> Qcurv = std::get<2>(propresult);
        const Matrix<double, 5, 7> FdFmcurv = std::get<3>(propresult);
        const double dEdxlast = std::get<4>(propresult);
        const Matrix<double, 5, 5> dQMScurv = std::get<5>(propresult);
        const Matrix<double, 5, 5> dQIcurv = std::get<6>(propresult);
//         const Matrix<double, 5, 5> dQcurv = Qcurv;
        const double deltaTotal = std::get<7>(propresult);
        const double wTotal = std::get<8>(propresult);
        
        Matrix<double, 5, 7> FdFm = FdFmcurv;
        if (ihit == 0) {
          // extra jacobian from reference state to curvilinear potentially needed
          FdFm.leftCols<5>() = FdFmcurv.leftCols<5>()*ref2curvjac;
        }
        

        // propfromtsos = updtsos;

//         //zero energy loss contribution
//         FdFm.rightCols<1>() *= 0.;


        Qtot = (FdFmcurv.leftCols<5>()*Qtot*FdFmcurv.leftCols<5>().transpose()).eval();
        Qtot += Qcurv;

        if (simhitdebug) {
          Qtot = Qcurv;
        }

        const Matrix<double, 5, 5> Hm = curv2localJacobianAltelossD(updtsos, field, surface, dEdxlast, mmu, dbetaval);
        
        
        const float enext = simhit == nullptr ? -99. : std::sqrt(std::pow(simhit->pabs(), 2) + mmu*mmu) - 0.5*simhit->energyLoss();        
        const float eprednext = std::sqrt(updtsos.segment<3>(3).squaredNorm() + mmu*mmu);
        
        const float dEval = e > 0. && enext > 0. ? enext - e : -99.;
        const float dEpredval = eprednext - epred;
        
        e = enext;
        epred = eprednext;
        
        const float sigmadEval = std::pow(updtsos.segment<3>(3).norm(), 3)/e*std::sqrt(Qcurv(0, 0));
        

        const Matrix<double, 6, 1> localparmsprop = globalToLocal(updtsos, surface);

        Matrix<double, 6, 1> localparms = localparmsprop;
        
        // update state from previous iteration
        //momentum kink residual
        Matrix<double, 5, 1> dx0 = Matrix<double, 5, 1>::Zero();
        
        if (iiter == 0 && fitFromSimParms_) {
          if (simhit == nullptr) {
            std::cout << "ABORT: Propagating from sim parameters, but sim hit is not matched!\n";
            valid = false;
            break;
          }

          auto const &surfaceideal = surfacemapIdealD_.at(hit->geographicalId());


          // alternate version with propagation from entry state

          const Point3DBase<double, LocalTag> simlocalpos = simhit->entryPoint();
          const Vector3DBase<double, LocalTag> simlocalmom = simhit->momentumAtEntry();

          const Point3DBase<double, GlobalTag> simglobalpos = surfaceideal.toGlobal(simlocalpos);
          const Vector3DBase<double, GlobalTag> simglobalmom = surfaceideal.toGlobal(simlocalmom);

          updtsos[0] = simglobalpos.x();
          updtsos[1] = simglobalpos.y();
          updtsos[2] = simglobalpos.z();
          updtsos[3] = simglobalmom.x();
          updtsos[4] = simglobalmom.y();
          updtsos[5] = simglobalmom.z();
          updtsos[6] = genpart->charge();

          if (false) {
            auto const &propresultsim = g4prop->propagateGenericWithJacobianAltD(updtsos, surface, dbetaval, dxival, dmsval, dionival);

            if (!std::get<0>(propresultsim)) {
              std::cout << "Abort: Sim state Propagation Failed!" << std::endl;
              valid = false;
              break;
            }

            updtsos = std::get<1>(propresultsim);
          }

          localparms = globalToLocal(updtsos, surface);

          dx0 = (localparms - localparmsprop).head<5>();

        }
        
        
        if (dolocalupdate) {
          if (iiter==0) {
            layerStates.push_back(updtsos);
          }
          else {
            //current state from previous state on this layer
            //save current parameters

            Matrix<double, 7, 1>& oldtsos = layerStates[ihit];
            const Matrix<double, 5, 5> Hold = curv2localJacobianAltelossD(oldtsos, field, surface, dEdxlast, mmu, dbetaval);
            const Matrix<double, 5, 1> dxlocal = Hold*dxfull.segment<5>(5*(ihit+1));

            localparms = globalToLocal(oldtsos, surface);

            localparms.head<5>() += dxlocal;

            oldtsos = localToGlobal(localparms, surface);

            updtsos = oldtsos;

            dx0 = (localparms - localparmsprop).head<5>();
          }
        }
        
        Matrix<double, 6, 1> localparmsalignprop = localparms;
        Matrix<double, 7, 1> updtsosalign = updtsos;
        
        if (fitFromSimParms_) {
        // re-propagate to unmodified surface

          
          auto const &surfaceideal = surfacemapIdealD_.at(hit->geographicalId());


          // alternate version with propagation from entry state

          const Point3DBase<double, LocalTag> simlocalpos = simhit->entryPoint();
          const Vector3DBase<double, LocalTag> simlocalmom = simhit->momentumAtEntry();

          const Point3DBase<double, GlobalTag> simglobalpos = surfaceideal.toGlobal(simlocalpos);
          const Vector3DBase<double, GlobalTag> simglobalmom = surfaceideal.toGlobal(simlocalmom);

          Matrix<double, 7, 1> simtsos;
          simtsos[0] = simglobalpos.x();
          simtsos[1] = simglobalpos.y();
          simtsos[2] = simglobalpos.z();
          simtsos[3] = simglobalmom.x();
          simtsos[4] = simglobalmom.y();
          simtsos[5] = simglobalmom.z();
          simtsos[6] = genpart->charge();          
          
          auto const &surfacealign = surfacemapD_.at(hit->geographicalId());
          auto const &propresultsalign = g4prop->propagateGenericWithJacobianAltD(simtsos, surfacealign, dbetaval, dxival, dmsval, dionival);

          if (!std::get<0>(propresultsalign)) {
            std::cout << "WARNING propagation for alignment failed!\n";
            valid = false;
            break;
          }
          
          updtsosalign = std::get<1>(propresultsalign);
          
          localparmsalignprop = globalToLocal(updtsosalign, surfacealign);
        }

        //TODO optimize this without ternary functions
        
        Matrix<double, 5, 5> Hp = Hm;
        
        if (dolocalupdate) {
          Hp = curv2localJacobianAltelossD(updtsos, field, surface, dEdxlast, mmu, dbetaval);
        }

        Matrix<double, 5, 5> Q = Qcurv;
        if (dolocalupdate) {
          Q = Hm*Qcurv*Hm.transpose();
        }
        
        const Matrix<double, 5, 5> Qinv = Q.inverse();
        
        Matrix<double, 5, 5> dQMS = dQMScurv;
        Matrix<double, 5, 5> dQI = dQIcurv;
        if (dolocalupdate) {
          dQMS = Hm*dQMScurv*Hm.transpose();
          dQI = Hm*dQIcurv*Hm.transpose();
        }
        
        {
          constexpr unsigned int nlocalcons = 5;
          constexpr unsigned int nlocalstateparms = 5;
          constexpr unsigned int nlocalbfield = 1;
          constexpr unsigned int nlocaleloss = 1;
          constexpr unsigned int nlocalparms = nlocalbfield + nlocaleloss;

          const unsigned int fullstateidx = 5*ihit;
          const unsigned int fullparmidx = iparm;

          rfull.segment<nlocalcons>(icons) = dx0;
          
          if (dolocalupdate) {
            Ffull.block<nlocalcons, nlocalstateparms>(icons, fullstateidx) = -Hm*FdFm.leftCols<nlocalstateparms>();
            Ffull.block<nlocalcons, nlocalstateparms>(icons, fullstateidx + nlocalstateparms) = Hp;
            Jfull.block<nlocalcons, nlocalparms>(icons, fullparmidx) = -Hm*FdFm.rightCols<nlocalparms>();
          }
          else {
            Ffull.block<nlocalcons, nlocalstateparms>(icons, fullstateidx) = -FdFm.leftCols<nlocalstateparms>();
            Ffull.block<nlocalcons, nlocalstateparms>(icons, fullstateidx + nlocalstateparms) = Matrix<double, nlocalcons, nlocalstateparms>::Identity();
            Jfull.block<nlocalcons, nlocalparms>(icons, fullparmidx) = -FdFm.rightCols<nlocalparms>();
          }
          
          Vinvfull.block<nlocalcons, nlocalcons>(icons, icons) = Qinv;
//           Vinvfullalt.block<nlocalcons, nlocalcons>(icons, icons) = Qinv;
          Vinvfullalt.block<nlocalcons, nlocalcons>(icons, icons) = (Q - 0.1*dQMS).inverse();
          
          if (dores) {
            std::vector<Triplet<double>> coeffs;
            for (unsigned int irow = 0; irow < nlocalcons; ++irow) {
              for (unsigned int icol = 0; icol < nlocalcons; ++icol) {
                coeffs.emplace_back(icons + irow, icons + icol, dQMS(irow, icol));
              }
            }
            SparseMatrix<double> &dV = dVs.emplace_back(ncons, ncons);
            dV.setFromTriplets(coeffs.begin(), coeffs.end());
            // MS resolution parameter is decoupled from energy loss
            residxs.push_back(iparm + 2);
          }

          if (dores) {
            std::vector<Triplet<double>> coeffs;
            for (unsigned int irow = 0; irow < nlocalcons; ++irow) {
              for (unsigned int icol = 0; icol < nlocalcons; ++icol) {
                coeffs.emplace_back(icons + irow, icons + icol, dQI(irow, icol));
              }
            }
            SparseMatrix<double> &dV = dVs.emplace_back(ncons, ncons);
            dV.setFromTriplets(coeffs.begin(), coeffs.end());
            // energy loss resolution parameter is decoupled from the energy loss mean value
            residxs.push_back(iparm + 3);
          }

          
          icons += nlocalcons;

          globalidxv[iparm] = bfieldglobalidx;
          iparm++;

          globalidxv[iparm] = elossglobalidx;
          iparm++;

          if (dores) {
            globalidxv[iparm] = msglobalidx;
            iparm++;
            
            globalidxv[iparm] = ioniglobalidx;
            iparm++;
          }
        }

        if (hit->isValid()) {

          //apply measurement update if applicable
          LocalTrajectoryParameters locparm(localparmsalignprop[0],
                                            localparmsalignprop[1],
                                            localparmsalignprop[2],
                                            localparmsalignprop[3],
                                            localparmsalignprop[4],
                                            localparmsalignprop[5]);
          const TrajectoryStateOnSurface tsostmp(locparm, *hit->surface(), field);

          auto const& preciseHit = cloner.makeShared(hit, tsostmp);
          if (!preciseHit->isValid()) {
            std::cout << "Abort: Failed updating hit" << std::endl;
            valid = false;
            break;
          }

          const bool align2d = detidparms.count(std::make_pair(1, aligndetid));

          const Matrix<double, 2, 2> &Rglued = rgluemap_.at(preciseHit->geographicalId());
          const GloballyPositioned<double> &surfaceglued = surfacemapD_.at(parmdetid);

          const Matrix<double, 5, 5> curvcov = covfull.block<5, 5>(5*(ihit+1), 5*(ihit+1));

          const Matrix<double, 2, 1> localconv = localPositionConvolutionD(updtsos, curvcov, surface);

          {
            constexpr unsigned int nlocalstate = 2;

            const unsigned int nlocalalignment = align2d ? 6 : 5;
            const unsigned int nlocalparms = nlocalalignment;
            
            const unsigned int fullstateidx = 5*(ihit+1) + 3;
            const unsigned int fullparmidx = iparm;

            const bool ispixel = GeomDetEnumerators::isTrackerPixel(preciseHit->det()->subDetector());
            
            const bool hit1d = preciseHit->dimension() == 1;
            
            const Matrix<double, 2, 2> Hu = Hp.bottomRightCorner<2,2>();

            Matrix<double, 2, 1> dy0;
            Matrix<double, 2, 2> iV;
            // rotation from module to strip coordinates
            Matrix2d R;
            
            const double lxcor = localparmsalignprop[3];
            const double lycor = localparmsalignprop[4];


            const Topology &topology = preciseHit->det()->topology();

            // undo deformation correction
            const LocalPoint lpnull(0., 0.);
            const MeasurementPoint mpnull = topology.measurementPosition(lpnull);
            const Topology::LocalTrackPred pred(tsostmp.localParameters().vector());

            auto const defcorr = topology.localPosition(mpnull, pred) - topology.localPosition(mpnull);

            const double hitx = preciseHit->localPosition().x() - defcorr.x();
            const double hity = preciseHit->localPosition().y() - defcorr.y();

            double dxrecsimval = -99.;
            double dyrecsimval = -99.;

            if (simhit != nullptr) {
              dxrecsimval = hitx - simhit->localPosition().x();
              dyrecsimval = hity - simhit->localPosition().y();
            }

            double lyoffset = 0.;
            double hitphival = -99.;
            double localphival = -99.;
            

            if (hit1d) {
              const ProxyStripTopology *proxytopology = dynamic_cast<const ProxyStripTopology*>(&(preciseHit->det()->topology()));

              dy0[0] = hitx - lxcor;
              dy0[1] = hity - lycor;

              const double striplength = proxytopology->stripLength();
              const double yerr2 = striplength*striplength/12.;
              
              iV = Matrix<double, 2, 2>::Zero();
              iV(0, 0) = preciseHit->localPositionError().xx();
              //iV(1, 1) = yerr2;
              
              R = Matrix2d::Identity();


//               std::cout << "1d hit, original x = " << preciseHit->localPosition().x() << " y = " << preciseHit->localPosition().y() << " corrected x = " << hitx << " y = " << hity << std::endl;
            }
            else {
              // 2d hit
//               assert(align2d);
              
              
              if (ispixel) {

                dy0[0] = hitx - lxcor;
                dy0[1] = hity - lycor;

                iV << preciseHit->localPositionError().xx(), preciseHit->localPositionError().xy(),
                      preciseHit->localPositionError().xy(), preciseHit->localPositionError().yy();

                R = Matrix2d::Identity();
              }
              else {
                // transform to polar coordinates to end the madness
                //TODO handle the module deformations consistently here (currently equivalent to dropping/undoing deformation correction)

//                   std::cout << "wedge\n" << std::endl;

                const ProxyStripTopology *proxytopology = dynamic_cast<const ProxyStripTopology*>(&(preciseHit->det()->topology()));

                const TkRadialStripTopology *radialtopology = dynamic_cast<const TkRadialStripTopology*>(&proxytopology->specificTopology());

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

                const double rhoerr2 = striplength*striplength/12.;


//                   std::cout << "rhohit = " << rhohit << " rhobar = " << rhobar << " rhoerr2lin = " << rhoerr2lin << " rhoerr2 = " << rhoerr2 << std::endl;

                // TODO apply (inverse) corrections for module deformations here? (take into account for jacobian?)
                const double phistate = rdir*std::atan2(lxcor, rdir*lycor + radius);
                const double rhostate = std::sqrt(lxcor*lxcor + std::pow(rdir*lycor + radius, 2));

                if (simhit != nullptr) {
                  const double lxsim = simhit->localPosition().x();
                  const double lysim = simhit->localPosition().y();

                  const double phisim = rdir*std::atan2(lxsim, rdir*lysim + radius);
                  const double rhosim = std::sqrt(lxsim*lxsim + std::pow(rdir*lysim + radius, 2));

                  dxrecsimval = phihit - phisim;
                  dyrecsimval = rhohit - rhosim;
                }

                iV = Matrix<double, 2, 2>::Zero();
                iV(0, 0) = phierr2;
                //iV(1, 1) = rhoerr2;

                // jacobian from localx-localy to localphi-localrho
                R = Matrix2d::Zero();

                const double yp = rdir*lycor + radius;
                const double invden = 1./(lxcor*lxcor + yp*yp);

                // dphi / dx
                R(0, 0) = rdir*yp*invden;
                // dphi / dy
                R(0, 1) = -lxcor*invden;
                // drho / dx
                R(1, 0) = lxcor/rhostate;
                // drho / dy
                R(1, 1) = rdir*(rdir*lycor + radius)/rhostate;


                dy0[0] = phihit - phistate;
                dy0[1] = rhohit - rhostate;

//                 std::cout << "wedge hit, original x = " << preciseHit->localPosition().x() << " y = " << preciseHit->localPosition().y() << " corrected x = " << hitx << " y = " << hity << std::endl;

              }
            }
            
            rxfull.row(ivalidhit) = R.row(0).cast<float>();
            ryfull.row(ivalidhit) = R.row(1).cast<float>();
            
            validdxeigjac.block<2,2>(2*ivalidhit, 3*(ihit+1)) = R*Hp.bottomRightCorner<2,2>();
            
            // alignment jacobian
            Matrix<double, 2, 6> Aval = Matrix<double, 2, 6>::Zero();

//             const Matrix<double, 6, 1> &localparmsalign = alignGlued_ ? globalToLocal(updtsos, surfaceglued) : localparms;
            
//             Matrix<double, 6, 1> localparmsalign = localparms;
            Matrix<double, 6, 1> localparmsalign = localparmsalignprop;
            if (alignGlued_) {
              localparmsalign = globalToLocal(updtsosalign, surfaceglued);
            }

            const double localqopval = localparmsalign[0];
            const double localdxdzval = localparmsalign[1];
            const double localdydzval = localparmsalign[2];
            const double localxval = localparmsalign[3];
            const double localyval = localparmsalign[4];
                        
            //standard case

            // dx/dx
            Aval(0,0) = -1.;
            // dy/dy
            Aval(1,1) = -1.;
            // dx/dz
            Aval(0,2) = localdxdzval;
            // dy/dz
            Aval(1,2) = localdydzval;
            // dx/dtheta_x
            Aval(0,3) = localyval*localdxdzval;
            // dy/dtheta_x
            Aval(1,3) = localyval*localdydzval;
            // dx/dtheta_y
            Aval(0,4) = -localxval*localdxdzval;
            // dy/dtheta_y
            Aval(1,4) = -localxval*localdydzval;
            // dx/dtheta_z
            Aval(0,5) = localyval;
            // dy/dtheta_z
            Aval(1,5) = -localxval;

//             const Matrix<double, 2, 6> &A = alignGlued_ ? Rglued*Aval : Aval;
            
            Matrix<double, 2, 6> A = R*Aval;
            if (alignGlued_) {
              A = R*Rglued*Aval;
            }
            
            const unsigned int xresglobalidx = dores ? detidparms.at(std::make_pair(8, hit->geographicalId())) : 0;
            const unsigned int yresglobalidx = dores && ispixel ? detidparms.at(std::make_pair(9, hit->geographicalId())) : 0;

            // apply correction to hit covariance from previous iterations if applicable
            const double xxresfact = dores ? std::exp(corparms_[xresglobalidx]) : 1.;
            const double yyresfact = dores && ispixel ? std::exp(corparms_[yresglobalidx]) : 1.;
            const double xyresfact = std::sqrt(xxresfact*yyresfact);
            
//             std::cout << "xxresfact = " << xxresfact << " yyresfact = " << yyresfact << " xyresfact = " << xyresfact << std::endl;

            iV(0, 0) *= xxresfact;
            iV(0, 1) *= xyresfact;
            iV(1, 0) *= xyresfact;
            iV(1, 1) *= yyresfact;
            
            // local x/phi resolution variation
            if (dores) {
              // scale diagonal element by exp(parm), preserve correlations
              const double dVxx = iV(0, 0);
              
              std::vector<Triplet<double>> coeffs;
              coeffs.reserve(3);
              coeffs.emplace_back(icons, icons, dVxx);
              
              // cross terms only relevant/exist for pixel hits
              if (ispixel) {
                const double dVxy = 0.5*iV(0, 1);
                coeffs.emplace_back(icons, icons + 1, dVxy);
                coeffs.emplace_back(icons + 1, icons, dVxy);
              }

              SparseMatrix<double> &dV = dVs.emplace_back(ncons, ncons);
              dV.setFromTriplets(coeffs.begin(), coeffs.end());
              residxs.push_back(iparm + nlocalalignment);
            }
            
            // local y resolution variation
            if (dores && ispixel) {
              // scale diagonal element by exp(parm), preserve correlations
              const double dVxy = 0.5*iV(0, 1);
              const double dVyy = iV(1, 1);
              
              std::vector<Triplet<double>> coeffs;
              coeffs.reserve(3);
              coeffs.emplace_back(icons, icons + 1, dVxy);
              coeffs.emplace_back(icons + 1, icons, dVxy);
              coeffs.emplace_back(icons + 1, icons + 1, dVyy);

              SparseMatrix<double> &dV = dVs.emplace_back(ncons, ncons);
              dV.setFromTriplets(coeffs.begin(), coeffs.end());
              residxs.push_back(iparm + nlocalalignment + 1);
            }

            constexpr std::array<unsigned int, 6> alphaidxs = {{0, 2, 3, 4, 5, 1}};

//             const double scalecov = hit1d ? 1.2 : 1.0;
            const double scalecov = ispixel ? 0.8 : 1.2;
            
            if (ispixel) {
              constexpr unsigned int nlocalcons = 2;
              
              rfull.segment<nlocalcons>(icons) = dy0;
              Ffull.block<nlocalcons, nlocalstate>(icons, fullstateidx) = -R*Hu;
              
              for (unsigned int ialign = 0; ialign < nlocalalignment; ++ialign) {
                Jfull.block<nlocalcons, 1>(icons, fullparmidx + ialign) = -A.col(alphaidxs[ialign]);
              }
              
              Vinvfull.block<nlocalcons, nlocalcons>(icons, icons) = iV.inverse();
              Vinvfullalt.block<nlocalcons, nlocalcons>(icons, icons) = (scalecov*iV).inverse();
              
              icons += nlocalcons;
            }
            else {
              constexpr unsigned int nlocalcons = 1;
              
              rfull(icons) = dy0(0);
              Ffull.block<1, nlocalstate>(icons, fullstateidx) = (-R*Hu).row(0);
              
              for (unsigned int ialign = 0; ialign < nlocalalignment; ++ialign) {
                Jfull(icons, fullparmidx + ialign) = -A(0, alphaidxs[ialign]);
              }
              
              Vinvfull(icons, icons) = 1./iV(0, 0);
              Vinvfullalt(icons, icons) = 1./(scalecov*iV)(0, 0);
              
              icons += nlocalcons;
            }
            
            for (unsigned int idim=0; idim<nlocalalignment; ++idim) {
              const unsigned int xglobalidx = detidparms.at(std::make_pair(alphaidxs[idim], aligndetid));
              globalidxv[iparm] = xglobalidx;
              iparm++;
              if (alphaidxs[idim]==0) {
                hitidxv.push_back(xglobalidx);
              }
            }
            
            if (dores) {
              globalidxv[iparm] = xresglobalidx;
              ++iparm;
            }
            
            if (dores && ispixel) {
              globalidxv[iparm] = yresglobalidx;
              ++iparm;
            }
            
            localqop_iter[ivalidhit] = localqopval;
            localdxdz_iter[ivalidhit] = localdxdzval;
            localdydz_iter[ivalidhit] = localdydzval;
            localx_iter[ivalidhit] = localxval;
            localy_iter[ivalidhit] = localyval;

            dxrecgen_iter[ivalidhit] = dy0[0];
            dyrecgen_iter[ivalidhit] = dy0[1];

            const double localqopvar = covfull(5*(ihit+1), 5*(ihit+1));
            localqoperralt[ivalidhit] = std::sqrt(localqopvar);
            
            if (iiter == 0) {
              
              // fill hit validation information
//               Vector2d dyrecgenlocal;
//               dyrecgenlocal << dy0[0].value().value(), dy0[1].value().value();
//               const Vector2d dyrecgeneig = R*dyrecgenlocal;
//               dxrecgen.push_back(dyrecgeneig[0]);
//               dyrecgen.push_back(dyrecgeneig[1]);
              dxrecgen.push_back(dy0[0]);
              dyrecgen.push_back(dy0[1]);
              
              dxerr.push_back(std::sqrt(iV(0,0)));
              dyerr.push_back(std::sqrt(iV(1,1)));

              const Matrix<double, 6, 1> localstatedebug = globalToLocal(updtsos, surface);

              localqop.push_back(localstatedebug[0]);
              localdxdz.push_back(localstatedebug[1]);
              localdydz.push_back(localstatedebug[2]);
              localx.push_back(localstatedebug[3]);
              localy.push_back(localstatedebug[4]);

//               localqop.push_back(localqopval);
//               localdxdz.push_back(localdxdzval);
//               localdydz.push_back(localdydzval);
//               localx.push_back(localxval);
//               localy.push_back(localyval);
              
              const Matrix<double, 5, 5> Qtotlocal = Hp*Qtot*Hp.transpose();
              
              localqoperr.push_back(std::sqrt(Qtotlocal(0, 0)));
              localdxdzerr.push_back(std::sqrt(Qtotlocal(1, 1)));
              localdydzerr.push_back(std::sqrt(Qtotlocal(2, 2)));
              localxerr.push_back(std::sqrt(Qtotlocal(3, 3)));
              localyerr.push_back(std::sqrt(Qtotlocal(4, 4)));
              
              hitlocalx.push_back(hitx);
              hitlocaly.push_back(hity);
              
              dEpred.push_back(dEpredval);
              sigmadE.push_back(sigmadEval);
              Epred.push_back(epred);
              E.push_back(e);
              
              hitphi.push_back(hitphival);
              localphi.push_back(localphival);

              landauDelta.push_back(deltaTotal);
              landauW.push_back(wTotal);

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
                  
                  // Vector2d dyrecsimlocal;
                  // dyrecsimlocal << preciseHit->localPosition().x() - simhit->localPosition().x(),
                  //                 preciseHit->localPosition().y() - simhit->localPosition().y();
                  // const Vector2d dyrecsimeig = R*dyrecsimlocal;
                  // dxrecsim.push_back(dyrecsimeig[0]);
                  // dyrecsim.push_back(dyrecsimeig[1]);

                  dxrecsim.push_back(dxrecsimval);
                  dyrecsim.push_back(dyrecsimval);


                  
                  dE.push_back(dEval);
                  
                  auto const momsim = simhit->momentumAtEntry();

                  const double eentry = std::sqrt(std::pow(simhit->pabs(), 2) + mmu*mmu);
                  const double emid = eentry - 0.5*simhit->energyLoss();
                  const double simqopval = genpart->charge()/std::sqrt(emid*emid - mmu*mmu);
//                   std::cout << "eloss = " << simhit->energyLoss() << std::endl;
                  
                  // "hybrid state" trying to adjust for entry point -> midpoint
                  // simlocalqop.push_back(simqopval);
                  // simlocaldxdz.push_back(momsim.x()/momsim.z());
                  // simlocaldydz.push_back(momsim.y()/momsim.z());
                  // simlocalx.push_back(simhit->localPosition().x());
                  // simlocaly.push_back(simhit->localPosition().y());

                  // just fill entry state for debugging
                  simlocalqop.push_back(genpart->charge()/simhit->pabs());
                  simlocaldxdz.push_back(momsim.x()/momsim.z());
                  simlocaldydz.push_back(momsim.y()/momsim.z());
                  simlocalx.push_back(simhit->entryPoint().x());
                  simlocaly.push_back(simhit->entryPoint().y());
                    
                  if (false) {

                    const Point3DBase<double, LocalTag> simlocalpos = simhit->entryPoint();
                    const Vector3DBase<double, LocalTag> simlocalmom = simhit->momentumAtEntry();

          //           std::cout << "simlocalpos" << simlocalpos << std::endl;

                    const Point3DBase<double, GlobalTag> simglobalpos = surface.toGlobal(simlocalpos);
                    const Vector3DBase<double, GlobalTag> simglobalmom = surface.toGlobal(simlocalmom);

                    Matrix<double, 7, 1> simtsos;
                    simtsos[0] = simglobalpos.x();
                    simtsos[1] = simglobalpos.y();
                    simtsos[2] = simglobalpos.z();
                    simtsos[3] = simglobalmom.x();
                    simtsos[4] = simglobalmom.y();
                    simtsos[5] = simglobalmom.z();
                    simtsos[6] = genpart->charge();

                    auto propresultsim = g4prop->propagateGenericWithJacobianAltD(simtsos, surface, dbetaval, dxival);

                    if (std::get<0>(propresultsim)) {
                      simtsos = std::get<1>(propresultsim);

                      const Point3DBase<double, GlobalTag> simglobalposprop(simtsos[0], simtsos[1], simtsos[2]);
                      const Vector3DBase<double, GlobalTag> simglobalmomprop(simtsos[3], simtsos[4], simtsos[5]);

                      const Point3DBase<double, LocalTag> simlocalposprop = surface.toLocal(simglobalposprop);
                      const Vector3DBase<double, LocalTag> simlocalmomprop = surface.toLocal(simglobalmomprop);

                      Matrix<double, 5, 1> simlocalparms;
                      simlocalparms[0] = simtsos[6]/simtsos.segment<3>(3).norm();
                      simlocalparms[1] = simlocalmomprop.x()/simlocalmomprop.z();
                      simlocalparms[2] = simlocalmomprop.y()/simlocalmomprop.z();
                      simlocalparms[3] = simlocalposprop.x();
                      simlocalparms[4] = simlocalposprop.y(); 
                      
                      simlocalqopprop.push_back(simlocalparms[0]);
                      simlocaldxdzprop.push_back(simlocalparms[1]);
                      simlocaldydzprop.push_back(simlocalparms[2]);
                      simlocalxprop.push_back(simlocalparms[3]);
                      simlocalyprop.push_back(simlocalparms[4]);
                    }
                    else {
                      simlocalqopprop.push_back(-99.);
                      simlocaldxdzprop.push_back(-99.);
                      simlocaldydzprop.push_back(-99.);
                      simlocalxprop.push_back(-99.);
                      simlocalyprop.push_back(-99.); 
                    }
                  }
                  else {
                    simlocalqopprop.push_back(-99.);
                    simlocaldxdzprop.push_back(-99.);
                    simlocaldydzprop.push_back(-99.);
                    simlocalxprop.push_back(-99.);
                    simlocalyprop.push_back(-99.);
                  }
                  
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
                  
                  simlocalqopprop.push_back(-99.);
                  simlocaldxdzprop.push_back(-99.);
                  simlocaldydzprop.push_back(-99.);
                  simlocalxprop.push_back(-99.);
                  simlocalyprop.push_back(-99.);
                }
              }

            }
            
          };
          
          ivalidhit++;
            
        }
        
      }
            
      if (!valid) {
        break;
      }

      
      assert(iparm == npars);
      
//       Fsparse = Ffull.rightCols(nstatefree).sparseView();
      Fsparse = Ffull(Eigen::all, freestateidxs).sparseView();
      Vinvsparse = Vinvfull.sparseView();
      
      VinvF = Vinvsparse*Fsparse;


      
      Cinvd.compute(Fsparse.transpose()*VinvF);


      // randomize residuals given covariance matrix
      if (false) {
        const SparseMatrix<double> Vinvsparsealt = Vinvfullalt.sparseView();
        const SparseMatrix<double> VinvFalt = Vinvsparsealt*Fsparse;

        const SimplicialLDLT<SparseMatrix<double>> Cinvdalt(Fsparse.transpose()*VinvFalt);
        
        const SparseMatrix<double> FtVinvalt = VinvFalt.transpose();
        const SparseMatrix<double> Rsparsealt = Vinvsparsealt - VinvFalt*Cinvdalt.solve(FtVinvalt);
        const MatrixXd Ralt = Rsparsealt;

        const Eigen::SelfAdjointEigenSolver<MatrixXd> Reig(Ralt);

        rfull = VectorXd::Zero(ncons);
        for (int ieig = ncons - ndof; ieig < ncons; ++ieig) {
          const double rnd = gRandom->Gaus();
          rfull += (rnd/std::sqrt(Reig.eigenvalues()(ieig)))*Reig.eigenvectors().col(ieig);
        }

      }
      
//       SimplicialLLT<SparseMatrix<double>> lltsparse;
//       lltsparse.compute(Fsparse.transpose()*VinvF);
//       
//       const SparseMatrix<double> lfact = lltsparse.matrixL();
//       std::cout << "lfact size = " << lfact.rows()*lfact.cols() << " nonzeros = " << lfact.nonZeros() << std::endl;
      
//       const SparseMatrix<double> testsparse = Fsparse.transpose()*VinvF;
//       const MatrixXd testdense = testsparse;
//       SimplicialLDLT<SparseMatrix<double>> ldltsparse;
//       LDLT<MatrixXd> ldltdense(testdense.rows());
//       
//       constexpr unsigned int ntest = 10000;
// 
//       double sparsetime;
//       double densetime;
//       
//       {
//         auto start = std::chrono::high_resolution_clock::now();
//         for (unsigned int itest = 0; itest < ntest; ++itest) {
//           ldltsparse.compute(testsparse);
//         }
//         auto stop = std::chrono::high_resolution_clock::now();
//         
//         auto duration = stop - start;
//         sparsetime = duration.count();
//         std::cout << "sparse decomp: " << sparsetime << std::endl;
//       }
// 
//       {
//         auto start = std::chrono::high_resolution_clock::now();
//         for (unsigned int itest = 0; itest < ntest; ++itest) {
//           ldltdense.compute(testdense);
//         }
//         auto stop = std::chrono::high_resolution_clock::now();
//         
//         auto duration = stop - start;
//         densetime = duration.count();
//         std::cout << "dense decomp: " << densetime << std::endl;
//       }
//       
//       const double timeratio = sparsetime/densetime;
//       std::cout << "size = " << testdense.rows() << " fillfactor = " << double(testsparse.nonZeros())/double(testsparse.rows()*testsparse.cols()) << " sparse/dense = " << timeratio << std::endl;
      
      dxfree = -Cinvd.solve(VinvF.transpose()*rfull);

      const double deltachisq = rfull.transpose()*VinvF*dxfree;
      edmval = -deltachisq;
      
      dxfull = VectorXd::Zero(nstateparms);
      dxfull(freestateidxs) = dxfree;
      
      const Vector5d dxref = dxfull.head<5>();

//       std::cout << "iiter = " << iiter << std::endl;
//       std::cout << "dxfree.head<5>():\n" << dxfree.head<5>() << std::endl;
//       std::cout << "dxfull.head<5>():\n" << dxfull.head<5>() << std::endl;
//       std::cout << "dxref:\n" <<dxref << std::endl;

      covfull = MatrixXd::Zero(nstateparms, nstateparms);

      //TODO temporary needed to work around Eigen fancy indexing limitation (maybe fixed in newer Eigen versions?)
      covfull(freestateidxs, freestateidxs) = Cinvd.solve(MatrixXd::Identity(nstatefree, nstatefree)).eval();

      const Matrix<double, 5, 5> Cinner = covfull.topLeftCorner<5,5>();
      
      if (dogen) {
        edmvalref = 0.;
      }
      else {
        const Matrix<double, 5, 5> hessref = Cinner.inverse();
        const double deltachisqref = -0.5*dxref.transpose()*hessref*dxref;
        edmvalref = -deltachisqref;
      }

      
      //fill output with corrected state and covariance at reference point
      refParms.fill(0.);
      refCov.fill(0.);
      
      const double qbp = refFts[6]/refFts.segment<3>(3).norm();
      const double lam = std::atan(refFts[5]/std::sqrt(refFts[3]*refFts[3] + refFts[4]*refFts[4]));
      const double phi = std::atan2(refFts[4], refFts[3]);
      
      const double qbpupd = qbp + dxref[0];
      const double lamupd = lam + dxref[1];
      const double phiupd = phi + dxref[2];
      
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


        
      
      niter = iiter + 1;
      
      if (std::isnan(edmval) || std::isinf(edmval)) {
        std::cout << "WARNING: invalid parameter update!!!" << " edmval = " << edmval << " lamupd = " << lamupd << " deltachisqval = " << deltachisqval << std::endl;
        valid = false;
        break;
      }
      
      
//       std::cout << "iiter = " << iiter << " edmval = " << edmval << " edmvalref " << edmvalref << " deltachisqval = " << deltachisqval << " chisqval = " << chisqval << std::endl;

      if (anomDebug) {
        std::cout << "anomDebug: iiter = " << iiter << " edmval = " << edmval << " deltachisqval = " << deltachisqval << " chisqval = " << chisqval << std::endl;
      }

      
      if (iiter > 0 && dolocalupdate && edmval < 1e-5) {
        break;
      }
      else if (iiter > 0 && !dolocalupdate && edmvalref < 1e-5) {
        break;
      }
    
    }
    
    if (!valid) {
      // chisqval = -99.;
      // tree->Fill();
      continue;
    }
        
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
    
    MatrixXd Jfinal = MatrixXd::Zero(ncons, nparsfinal);
    
    for (unsigned int i = 0; i < npars; ++i) {
      const unsigned int iidx = idxmap.at(globalidxv[i]);
      Jfinal.col(iidx) += Jfull.col(i);
    }
    
    const SparseMatrix<double> Jsparse = Jfinal.sparseView();
        
    std::vector<unsigned int> residxsfinal;
    residxsfinal.reserve(residxs.size());
    for (auto idx : residxs) {
      const unsigned int iidx = idxmap.at(globalidxv[idx]);
      residxsfinal.push_back(iidx);
    }
    
//     dxdparms = -Cinvd.solve(d2chisqdxdparmsfinal).transpose();
    
//     grad = dchisqdparmsfinal + dxdparms*dchisqdx; 
//     grad = dchisqdparmsfinal + d2chisqdxdparmsfinal.transpose()*dxfull;
//     hess = d2chisqdparms2final + dxdparms*d2chisqdxdparmsfinal;
    
    nParms = nparsfinal;

    //TODO avoid explicitly storing the transpose?
    const SparseMatrix<double> FtVinv = VinvF.transpose();
//     const MatrixXd R = Vinvsparse - VinvF*Cinvd.solve(FtVinv);
    const SparseMatrix<double> Rsparse = Vinvsparse - VinvF*Cinvd.solve(FtVinv);
    const MatrixXd R = Rsparse;
    
    
//     constexpr unsigned int ntest = 1000;
//     
//     auto time0 = std::chrono::high_resolution_clock::now();
//     for (unsigned int itest = 0; itest < ntest; ++itest) {
//       const SparseMatrix<double> R1 = Vinvsparse - VinvF*Cinvd.solve(FtVinv);
//       const MatrixXd R2 = R1;
//     }
//     
//     auto time1 = std::chrono::high_resolution_clock::now();
//     for (unsigned int itest = 0; itest < ntest; ++itest) {
//       const MatrixXd R1 = Vinvsparse - VinvF*Cinvd.solve(FtVinv);
//       const SparseMatrix<double> R2 = R1.sparseView();
//     }
//     
//     auto time2 = std::chrono::high_resolution_clock::now();
// 
//     
//     const double d0 = (time1-time0).count();
//     const double d1 = (time2-time1).count();
// //       const double d2 = (time3-time2).count();
// //       const double d3 = (time4-time3).count();
// //       
//       std::cout << "d0 = " << d0 << " d1 = " << d1 << std::endl;
    
    
//     std::vector<Triplet<Matrix<double, 3, 3>>> triplets;
// //     SparseMatrix<Matrix<double, Dynamic, Dynamic, 0, 5, 5>> testsparse(100, 100);
//     SparseMatrix<Matrix<double, 3, 3>> testsparse(100, 100);
//     triplets.emplace_back(0, 0, Matrix<double, 3, 3>::Zero());
//     
//     testsparse.setFromTriplets(triplets.begin(), triplets.end());
//     testsparse.insert(0, 0) = Matrix<double, 3, 3>::Identity();
    
//     std::cout << "R total size = " << R.rows()*R.cols() << " nonzeros = " << R.nonZeros() << std::endl;
    
//     std::cout << "VinvF total size = " << VinvF.rows()*VinvF.cols() << " nonzeros = " << VinvF.nonZeros() << std::endl;


    //TODO make R dense (and find the best solution for sparse = sparse*dense case)
    
    
//     chisqval = rfull.transpose()*R*rfull;
//     
//     grad = 2.*Jsparse.transpose()*R*rfull;
//     hess = 2.*Jsparse.transpose()*R*Jsparse;
    
    //TODO check this against explicit version
    const VectorXd Rr = Vinvsparse*(rfull + Fsparse*dxfree);
    
//     const VectorXd Rrtest = R*rfull;
//     const double Rrdiff = (Rr-Rrtest).array().square().sum();
//     std::cout << "Rrdiff = " << Rrdiff << std::endl;
    
//     constexpr unsigned int ntest = 1000;
//     
//     auto time0 = std::chrono::high_resolution_clock::now();
//     for (unsigned int itest = 0; itest < ntest; ++itest) {
//       const VectorXd Rrtest = Vinvsparse*(rfull + Fsparse*dxfull);
//     }
//     
//     auto time1 = std::chrono::high_resolution_clock::now();
//     for (unsigned int itest = 0; itest < ntest; ++itest) {
// //       VectorXd Rrtest = Vinvsparse*rfull;
// //       Rrtest += Vinvsparse*Fsparse*dxfull;
//       const VectorXd Rrtest = R*rfull;
//     }
//     
//     auto time2 = std::chrono::high_resolution_clock::now();
// 
//     
//     const double d0 = (time1-time0).count();
//     const double d1 = (time2-time1).count();
// 
//     std::cout << "d0 = " << d0 << " d1 = " << d1 << std::endl;
    
    
    chisqval = rfull.transpose()*Rr;
    
    grad = 2.*Jsparse.transpose()*Rr;
    hess = 2.*Jsparse.transpose()*R*Jsparse;
    
    
    dxdparms = MatrixXd::Zero(nparsfinal, nstateparms);
    
    dxdparms(Eigen::all, freestateidxs) = -Cinvd.solve(VinvF.transpose()*Jsparse).transpose();
    
    
    //additional contributions from resolution variations
    
    std::vector<SparseMatrix<double>> dVRs;
    dVRs.reserve(dVs.size());
    
    MatrixXd dVRr = MatrixXd::Zero(ncons, nparsfinal);
    
    VectorXd gradll = VectorXd::Zero(nparsfinal);
    
    for (unsigned int ires = 0; ires < dVs.size(); ++ires) {
      auto const &dVi = dVs[ires];
      const unsigned int residxi = residxsfinal[ires];
      
      SparseMatrix<double> &dViR = dVRs.emplace_back();
      dViR = dVi*Rsparse;
      
//       const MatrixXd Rdense = R;
      
//       const SparseMatrix<double> dViRtest = dVi*(Rdense.sparseView().transpose().transpose());
//       const SparseMatrix<double> dViRtest = Rdense.sparseView().transpose();

//       
//       MatrixXd tmpdense = MatrixXd::Zero(ncons, ncons);
//       
//       constexpr unsigned int ntest = 1000;
//       
//       auto time0 = std::chrono::high_resolution_clock::now();
//       for (unsigned int itest = 0; itest < ntest; ++itest) {
//         const SparseMatrix<double> dViRtest = (dVi*Rdense).sparseView();
//       }
//       auto time1 = std::chrono::high_resolution_clock::now();
// 
//       for (unsigned int itest = 0; itest < ntest; ++itest) {
//         const SparseMatrix<double> Rsparse = Rdense.sparseView();
//         const SparseMatrix<double> dViRtest = dVi*Rsparse;
//       }
//       
//       auto time2 = std::chrono::high_resolution_clock::now();
//       
//       for (unsigned int itest = 0; itest < ntest; ++itest) {
//         const SparseMatrix<double> dViRtest = dVi*R;
//       }
//       
//       auto time3 = std::chrono::high_resolution_clock::now();
//       
//       
//       for (unsigned int itest = 0; itest < ntest; ++itest) {
//         tmpdense = dVi*Rdense;
// //         const SparseMatrix<double> dViRtest = tmpdense.sparseView();
//       }
//       
//       auto time4 = std::chrono::high_resolution_clock::now();
//       
//       const double d0 = (time1-time0).count();
//       const double d1 = (time2-time1).count();
//       const double d2 = (time3-time2).count();
//       const double d3 = (time4-time3).count();
//       
//       std::cout << "d0 = " << d0 << " d1 = " << d1 << " d2 = " << d2 << " d3 = " << d3 << std::endl;
      
//       const double testfail = (dVi*MatrixXd::Identity(ncons, ncons)).sparseView();
      
//       const int testcount = (dVi*MatrixXd(R)).nonZeros();
//       const Matrix<double, 10, 10> testdense = Matrix<double, 10, 10>::Identity();
//       const int testcount2 = testdense.nonZeros();
//       
//       std::cout << "testcount2 = " << testcount2 << std::endl;
//       
//       std::cout << "dViR total size = " << dViR.rows()*dViR.cols() << " testcount = " << testcount << std::endl;

      
//       const MatrixXd ident = MatrixXd::Identity(ncons, ncons);
      
//       SparseMatrix<double> testsparse = SparseMatrix<double>(dVi*ident);
      
//       std::cout << "dViR total size = " << dViR.rows()*dViR.cols() << " nonzeros = " << dViR.nonZeros() << std::endl;
      
//       dVRr.col(residxi) += dViR*rfull;
      dVRr.col(residxi) += dVi*Rr;
      
      //TODO optimize traces for sparse matrices
      
//       const double gradres = MatrixXd(dViR).trace();
      
      //trace is not implemented directly for sparse matrices so use sum of diagonals explicitly
      const double gradres = dViR.diagonal().sum();
      
//       const double gradres = -rfull.transpose()*R*dViR*rfull + MatrixXd(dViR).trace();
      
//       grad(residxi) += gradres;
      gradll(residxi) += gradres;
            
//       VectorXd hesscross = -2.*Jsparse.transpose()*R*dViR*rfull;
      // adjust diagonal element since it will be added twice
//       hesscross(residxi) *= 0.5;
      
//       hess.col(residxi) += hesscross;
//       hess.row(residxi) += hesscross.transpose();
      
//       if (!dogen) {
//         dxdparms.row(residxi) += (M*dViR*rfull).transpose();
//       }
      
      const SparseMatrix<double> dViRT = dViR.transpose();
      
      for (unsigned int jres = 0; jres <= ires; ++jres) {
        auto const &dVjR = dVRs[jres];
        const unsigned int residxj = residxsfinal[jres];
        
//         const SparseMatrix<double> dViRdVjR = dViR*dVjR;
        
//         const double hessres = -MatrixXd(dViR*dVjR).trace();
//         const double hessres = 2.*rfull.transpose()*R*dViRdVjR*rfull - MatrixXd(dViRdVjR).trace();
        
        // below is equivalent to -(dViR*dVjR).trace()
//         const double hessres = -dViRT.cwiseProduct(dVjR).sum();
        const double hessres = dViRT.cwiseProduct(dVjR).sum();
//         const double hessresalt = -MatrixXd(dViR*dVjR).trace();
        
//         std::cout << "ires = " << ires << " jres = " << jres << " hessres = " << hessres << " hessresalt = " << hessresalt << std::endl;
        
        hess(residxi, residxj) += hessres;
        if (ires != jres) {
          hess(residxj, residxi) += hessres;
        }
        
      }
      

    
      
      
// //       const SparseMatrix<double> dViR = (dVi*R).sparseView();
//       const SparseMatrix<double> dViR = dVi*R;
// //       const MatrixXd dViR = dVi*R;
//       
//       //TODO optimized trace for sparse matrices
//       //TODO optimize use of sparse vs dense matrices (and avoid redundant conversions)
//       
//       grad(residxi) += -rfull.transpose()*R*dViR*rfull + MatrixXd(dViR).trace();
//       
//       VectorXd hesscross = -2.*Jsparse.transpose()*R*dViR*rfull;
//       
//       // adjust diagonal element since it will be added twice
//       hesscross(residxi) *= 0.5;
//       
//       hess.col(residxi) += hesscross;
//       hess.row(residxi) += hesscross.transpose();
//       
//       if (!dogen) {
//         dxdparms.row(residxi) += (M*dViR*rfull).transpose();
//       }
//       
//       for (unsigned int jres = 0; jres < dVs.size(); ++jres) {
//         auto const &dVj = dVs[jres];
//         const unsigned int residxj = residxsfinal[jres];
//         
// //         const SparseMatrix<double> dViRdVjR = (dViR*dVj*R).sparseView();
//         const SparseMatrix<double> dViRdVjR = dViR*dVj*R;
// //         const MatrixXd dViRdVjR = dViR*dVj*R;
//         
//         const double hessres = 2.*rfull.transpose()*R*dViRdVjR*rfull - MatrixXd(dViRdVjR).trace();
//         
//         hess(residxi, residxj) += hessres;
//         if (residxi != residxj) {
//           hess(residxj, residxi) += hessres;
//         }
//       }
      
    }
    
    const SparseMatrix<double> dVRrsparse = dVRr.sparseView();
    
//     std::cout << "dVRrsparse size = " << dVRrsparse.rows()*dVRrsparse.cols() << " nonzeros = " << dVRrsparse.nonZeros() << std::endl;
    
//     grad += -(rfull.transpose()*R*dVRrsparse).transpose();
//     grad += -(Rr.transpose()*dVRrsparse).transpose();
    grad += -dVRrsparse.transpose()*Rr;
    
    if (false) {
      hess += 2.*dVRrsparse.transpose()*R*dVRrsparse;
      
      const SparseMatrix<double> hesscross = -2.*Jsparse.transpose()*Rsparse*dVRrsparse;    
      hess += hesscross;
      hess += hesscross.transpose();
    }
    
//     hess += hesscross + hesscross.transpose();

    //TODO deduplicate with above
    dxdparms(Eigen::all, freestateidxs) += Cinvd.solve(FtVinv*dVRrsparse).transpose();
    
//     if (!dogen) {
//     //TODO deduplicate with above
// //       const MatrixXd M = Cinvd.solve(FtVinv);
// //       dxdparms += (M*dVRrsparse).transpose();
//       dxdparms += Cinvd.solve(FtVinv*dVRrsparse).transpose();
//     }
//     
    
    //TODO check this against element-wise version
//     grad += -(rfull.transpose()*R*dVRrsparse).transpose();
//     hess += 2.*dVRrsparse.transpose()*R*dVRrsparse;
//     
//     MatrixXd hesstest = 2.*dVRrsparse.transpose()*R*dVRrsparse;
//     MatrixXd hesstest2 = MatrixXd::Zero(nparsfinal, nparsfinal);
//     for (unsigned int ires = 0; ires < dVs.size(); ++ires) {
//       auto const &dVi = dVs[ires];
//       const unsigned int residxi = residxsfinal[ires];
//       
//       const SparseMatrix<double> &dViR = dVRs[ires];
//       
//       for (unsigned int jres = 0; jres <= ires; ++jres) {
//         auto const &dVjR = dVRs[jres];
//         const unsigned int residxj = residxsfinal[jres];
//         
//         const double hesstestval = 2.*rfull.transpose()*R*dViR*dVjR*rfull;
// //         const double hesstestval = 2.*rfull.transpose()*dViR*dVjR*rfull;
//         
// //         const double hesstestval2 = 2.*rfull.transpose()*dViR.transpose()*R*dVjR*rfull;
//         
//         hesstest2(residxi, residxj) += hesstestval;
// //         if (residxi != residxj) {
//         if (ires != jres) {
//           hesstest2(residxj, residxi) += hesstestval;
//         }
//         
// //         std::cout << "ires = " << ires << " jres = " << jres << " residxi = " << residxi << " residxj = " << residxj << " hesstestval = " << hesstestval << " hesstest(residxi, residxj) = " << hesstest(residxi, residxj) << std::endl;
//         
//       }
//     }
//     
//     for (unsigned int ipar = 0; ipar < nparsfinal; ++ipar) {
//       for (unsigned int jpar = 0; jpar < nparsfinal; ++jpar) {
//         const double diff = hesstest2(ipar, jpar) - hesstest(ipar, jpar);
//         
//         if (std::fabs(diff) > 1e-16 && std::fabs(diff)/hesstest(ipar, jpar) > 1e-3) {
//           std::cout << "ipar = " << ipar << " jpar = " << jpar << " hesstest(ipar, jpar) = " << hesstest(ipar, jpar) << " hessest2(ipar, jpar) = " << hesstest2(ipar, jpar) << std::endl;
//         }
//         
//       }
//     }
//     
//     const double diffsq = (hesstest2-hesstest).array().square().sum();
//     std::cout << "diffsq = " << diffsq << std::endl;
    
    
    
//     MatrixXd hesscross = -2.*Jsparse.transpose()*R*dVRr;
//     //correct diagonal since it gets added twice
//     hesscross.diagonal() *= 0.5;
    
//     hess += hesscross + hesscross.transpose();
    
//     std::cout << "hesscross.bottomRightCorner<10,10>():\n" << hesscross.bottomRightCorner<10,10>() << std::endl;
    
//     if (!dogen) {
//       dxdparms += (M*dVRrsparse).transpose();
//     }
    
    gradchisqv.clear();
    gradchisqv.resize(nparsfinal, 0.);
    
    Map<VectorXf>(gradchisqv.data(), nparsfinal) = grad.cast<float>();
    
    grad += gradll;
    
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
    
//     deigx.resize(nvalid);
//     deigy.resize(nvalid);
//     
//     validdxeig = validdxeigjac*dxfull;
//     
//     for (unsigned int ivalid = 0; ivalid < nvalid; ++ivalid) {
//       deigx[ivalid] = validdxeig[2*ivalid];
//       deigy[ivalid] = validdxeig[2*ivalid + 1];
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
