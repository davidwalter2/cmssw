#include "ResidualGlobalCorrectionMakerBase.h"
#include "MagneticFieldOffset.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "DataFormats/Math/interface/deltaR.h"


class ResidualGlobalCorrectionMaker : public ResidualGlobalCorrectionMakerBase
{
public:
  explicit ResidualGlobalCorrectionMaker(const edm::ParameterSet &);
  ~ResidualGlobalCorrectionMaker() {}

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
  
};


ResidualGlobalCorrectionMaker::ResidualGlobalCorrectionMaker(const edm::ParameterSet &iConfig) : ResidualGlobalCorrectionMakerBase(iConfig) 
{
  
  inputAssoc_ = consumes<edm::Association<reco::TrackExtraCollection>>(edm::InputTag("muonReducedTrackExtras"));
  
}

void ResidualGlobalCorrectionMaker::beginStream(edm::StreamID streamid)
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
    tree->Branch("nJacRef", &nJacRef, basketSize);
    
    tree->Branch("nValidHitsFinal", &nValidHitsFinal);
    tree->Branch("nValidPixelHitsFinal", &nValidPixelHitsFinal);
    
    tree->Branch("jacrefv",jacrefv.data(),"jacrefv[nJacRef]/F", basketSize);
    

    
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
    
    tree->Branch("dxreccluster", &dxreccluster);
    tree->Branch("dyreccluster", &dyreccluster);
    
    tree->Branch("localqop", &localqop);
    tree->Branch("localdxdz", &localdxdz);
    tree->Branch("localdydz", &localdydz);
    tree->Branch("localx", &localx);
    tree->Branch("localy", &localy);
    
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
    
    
    nJacRef = 0.;
  }
}


// ------------ method called for each event  ------------
void ResidualGlobalCorrectionMaker::produce(edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  
  const bool dogen = fitFromGenParms_;
  
  using namespace edm;

  Handle<reco::TrackCollection> trackOrigH;
  iEvent.getByToken(inputTrackOrig_, trackOrigH);
  
  
//   bool foundmodule = false;
//   for (const reco::Track &track : *trackOrigH) {
//     if (track.innerDetId() == 302055944) {
//       foundmodule = true;
//       break;
//     }
// //     for (auto it = track.recHitsBegin(); it != track.recHitsEnd(); ++it) {
// //       if ((*it)->geographicalId().rawId() == 302055944) {
// //         foundmodule = true;
// //         break;
// //       }
// //     }
// //     if (foundmodule) {
// //       break;
// //     }
//   }
//   if (!foundmodule) {
// //     printf("not found, returning\n");
//     return;
//   }
  


  // loop over gen particles

  edm::ESHandle<GlobalTrackingGeometry> globalGeometry;
  iSetup.get<GlobalTrackingGeometryRecord>().get(globalGeometry);
  
//   edm::ESHandle<TrackerGeometry> globalGeometry;
//   iSetup.get<TrackerDigiGeometryRecord>().get("idealForDigi", globalGeometry);
  
  edm::ESHandle<TrackerTopology> trackerTopology;
  iSetup.get<TrackerTopologyRcd>().get(trackerTopology);
  
//   ESHandle<MagneticField> magfield;
//   iSetup.get<IdealMagneticFieldRecord>().get(magfield);
//   auto field = magfield.product();
  
  edm::ESHandle<TransientTrackingRecHitBuilder> ttrh;
  iSetup.get<TransientRecHitRecord>().get("WithAngleAndTemplate",ttrh);
  
//   ESHandle<Propagator> thePropagator;
//   iSetup.get<TrackingComponentsRecord>().get("RungeKuttaTrackerPropagator", thePropagator);
//   iSetup.get<TrackingComponentsRecord>().get("PropagatorWithMaterial", thePropagator);
//   iSetup.get<TrackingComponentsRecord>().get("PropagatorWithMaterialParabolicMf", thePropagator);
//   iSetup.get<TrackingComponentsRecord>().get("Geant4ePropagator", thePropagator);
//   const MagneticField* field = thePropagator->magneticField();
  
//   edm::ESHandle<TrajectoryFitter> fit;
//   iSetup.get<TrajectoryFitter::Record>().get("G4eFitterSmoother", fit);
//   const KFTrajectoryFitter *kffit = dynamic_cast<const KFTrajectoryFitter*>(fit.product());
//   const Propagator *thePropagator = kffit->propagator();
  
  
//   ESHandle<Propagator> theAnalyticPropagator;
//   iSetup.get<TrackingComponentsRecord>().get("PropagatorWithMaterial", theAnalyticPropagator);
  
  ESHandle<MagneticField> fieldh;
  iSetup.get<IdealMagneticFieldRecord>().get("", fieldh);
  std::unique_ptr<MagneticFieldOffset> fieldOffset = std::make_unique<MagneticFieldOffset>(&(*fieldh));
  
  std::unique_ptr<PropagatorWithMaterial> fPropagator = std::make_unique<PropagatorWithMaterial>(alongMomentum, 0.105, fieldOffset.get(), 1.6, true,  -1., true);
  
  const MagneticField* field = fPropagator->magneticField();
  
//   Handle<TrajTrackAssociationCollection> trackH;
//   Handle<reco::TrackCollection> trackH;
//   iEvent.getByToken(inputTrack_, trackH);
  

  
//   Handle<std::vector<int> > indicesH;
//   iEvent.getByToken(inputIndices_, indicesH);
  
//   Handle<std::vector<Trajectory> > trajH;
//   iEvent.getByToken(inputTraj_, trajH);
  
  

  
  Handle<reco::BeamSpot> bsH;
  iEvent.getByToken(inputBs_, bsH);

  
  Handle<std::vector<reco::GenParticle>> genPartCollection;
  Handle<std::vector<int>> genPartBarcodes;
  if (doGen_) {
    iEvent.getByToken(GenParticlesToken_, genPartCollection);
    iEvent.getByToken(genParticlesBarcodeToken_, genPartBarcodes);
  }
  
//   Handle<std::vector<PSimHit>> tecSimHits;
  std::vector<Handle<std::vector<PSimHit>>> simHits(inputSimHits_.size());
  edm::Handle<std::vector<SimTrack>> simTracks;
  if (doSim_) {
    for (unsigned int isimhit = 0; isimhit<inputSimHits_.size(); ++isimhit) {
      iEvent.getByToken(inputSimHits_[isimhit], simHits[isimhit]);
    }
    iEvent.getByToken(inputSimTracks_, simTracks);
  }
  
  Handle<reco::MuonCollection> muons;
  if (doMuons_) {
    iEvent.getByToken(inputMuons_, muons);
  }
  
  Handle<edm::Association<reco::TrackExtraCollection>> assoc;
  if (doMuons_) {
    iEvent.getByToken(inputAssoc_, assoc);
  }
  
//   if (doSim_) {
//     iEvent.getByToken(inputSimHits_, tecSimHits);
//   }
  
//   const float mass = 0.105;
//   const float maxDPhi = 1.6;
//   PropagatorWithMaterial rPropagator(oppositeToMomentum, mass, field, maxDPhi, true, -1., false);
//   PropagatorWithMaterial fPropagator(alongMomentum, mass, field, maxDPhi, true, -1., false);
  
//   std::unique_ptr<PropagatorWithMaterial> fPropagator(static_cast<PropagatorWithMaterial*>(thePropagator->clone()));
//   fPropagator->setPropagationDirection(alongMomentum);
//   
//   std::unique_ptr<PropagatorWithMaterial> fAnalyticPropagator(static_cast<PropagatorWithMaterial*>(theAnalyticPropagator->clone()));
//   fAnalyticPropagator->setPropagationDirection(alongMomentum);
  
  KFUpdator updator;
  TkClonerImpl const& cloner = static_cast<TkTransientTrackingRecHitBuilder const *>(ttrh.product())->cloner();

  
//   siStripClusterInfo_.initEvent(iSetup);
  
//   edm::ESHandle<Alignments> globalPositionRcd;
//   iSetup.get<GlobalPositionRcd>().get(globalPositionRcd);
//   
//   printf("globalPositionRcd translation = %e, %e, %e\n", globalPositionRcd->m_align.front().translation().x(),
//                                                         globalPositionRcd->m_align.front().translation().y(),
//                                                         globalPositionRcd->m_align.front().translation().z());
//   std::cout << "globalPositionRcd rotation" << globalPositionRcd->m_align.front().rotation() << std::endl;
  
  // set up cylindrical surface for beam pipe
//   const double ABe = 9.0121831;
//   const double ZBe = 4.;
//   const double K =  0.307075*1e-3;
//   const double dr = 0.08;
// //   const double xibeampipe = 0.5*K*dr*ZBe/ABe;
//   const double xibeampipe = 0.*0.5*K*dr*ZBe/ABe;
  
  
  
  
//   auto beampipe = Cylinder::build(Surface::PositionType(0.,0.,0.), Surface::RotationType(), 2.94);
//   beampipe->setMediumProperties(MediumProperties(0., xibeampipe));
  
//   std::cout << "xi beampipe: " << xibeampipe << std::endl;
  
//   const GeomDet *testdet = nullptr;
//   //debugging
//   for (const GeomDet* det : globalGeometry->detUnits()) {
//     if (!det) {
//       continue;
//     }
//     
//     if (det->subDetector() == GeomDetEnumerators::TEC) {
//       const DetId& detid = det->geographicalId();
// //       TECDetId detid(det->geographicalId());
// //       layer = -1 * (detid.side() == 1) * detid.wheel() + (detid.side() == 2) * detid.wheel();
//       unsigned int side = trackerTopology->tecSide(detid);
//       unsigned int wheel = trackerTopology->tecWheel(detid);
//       int layer = -1 * (side == 1) * wheel + (side == 2) * wheel;
//       bool stereo = trackerTopology->isStereo(det->geographicalId());
//       
//       if (layer == -9) {
//         testdet = det;
//         break;
//         
//       }
//     }
//     
//     
//     
//   }
//   
//   if (testdet) {
//     const GlobalPoint center = testdet->surface().toGlobal(LocalPoint(1.,0.));
//     
//     const GlobalVector centerv(center.x(), center.y(), center.z());
//     const GlobalVector dir = centerv/centerv.mag();
//     const double sintheta = dir.perp();
//     const GlobalVector mom = (100000./sintheta)*dir;
//     const GlobalPoint pos(0.,0.,0.);
//     
//     FreeTrajectoryState ftsplus(pos, mom, 1., field);
//     FreeTrajectoryState ftsminus(pos, mom, -1., field);
//     
//     const TrajectoryStateOnSurface tsosplus = fPropagator->propagate(ftsplus, testdet->surface());
//     const TrajectoryStateOnSurface tsosminus = fPropagator->propagate(ftsminus, testdet->surface());
//     
//     std::cout << "global target" << std::endl;
//     std::cout << center << std::endl;
//     
//     std::cout << "momentum" << std::endl;
//     std::cout << mom << std::endl;
//     
//     std::cout << "tsosplus local:" << std::endl;
//     std::cout << tsosplus.localPosition() << std::endl;
//     std::cout << "tsosminus local:" << std::endl;
//     std::cout << tsosminus.localPosition() << std::endl;
//     
//     std::cout << "delta local" << std::endl;
//     std::cout << tsosplus.localPosition() - tsosminus.localPosition() << std::endl;
//     
//   }
  
  
  simtestz = -99.;
  simtestvz = -99.;
  simtestrho = -99.;
  simtestzlocalref = -99.;
  simtestdx = -99.;
  simtestdxrec = -99.;
  simtestdy = -99.;
  simtestdyrec = -99.;
  simtestdxprop = -99.;
  simtestdyprop = -99.;
  simtestdetid = 0;
  
  if (false) {
    
    //sim hit debugging
    const reco::GenParticle* genmuon = nullptr;
    for (const reco::GenParticle& genPart : *genPartCollection) {
      if (genPart.status()==1 && std::abs(genPart.pdgId()) == 13) {
        genmuon = &genPart;
        break;
      }
    }
    
    if (genmuon) {
      genPt = genmuon->pt();
      genCharge = genmuon->charge();
      genEta = genmuon->eta();
      genPhi = genmuon->phi();
      
      simtestvz = genmuon->vertex().z();
      
      auto const& refpoint = genmuon->vertex();
      auto const& trackmom = genmuon->momentum();
      const GlobalPoint refpos(refpoint.x(), refpoint.y(), refpoint.z());
      const GlobalVector refmom(trackmom.x(), trackmom.y(), trackmom.z()); 
  //     const GlobalTrajectoryParameters refglobal(refpos, refmom, genmuon->charge(), field);
      
  //       std::cout << "gen ref state" << std::endl;
  //       std::cout << refpos << std::endl;
  //       std::cout << refmom << std::endl;
  //       std::cout << genpart->charge() << std::endl;
      
      //zero uncertainty on generated parameters
  //       AlgebraicSymMatrix55 nullerr;
  //       const CurvilinearTrajectoryError referr(nullerr);
      
      const FreeTrajectoryState fts = FreeTrajectoryState(refpos, refmom, genmuon->charge(), field);
      
      std::cout << "gen muon charge: " << genmuon->charge() << std::endl;
      
      TrajectoryStateOnSurface tsos;

      std::vector<const PSimHit*> simhitsflat;
      for (auto const& simhith : simHits) {
        for (const PSimHit& simHit : *simhith) {
          simhitsflat.push_back(&simHit);
        }
      }
      
      auto simhitcompare = [&](const PSimHit *hit0, const PSimHit *hit1) {
        return globalGeometry->idToDet(hit0->detUnitId())->surface().toGlobal(hit0->localPosition()).mag() < globalGeometry->idToDet(hit1->detUnitId())->surface().toGlobal(hit1->localPosition()).mag();
      };
      
      std::sort(simhitsflat.begin(), simhitsflat.end(), simhitcompare);
      
      unsigned int ihit = 0;
      for (auto const& simHitp : simhitsflat) {
          auto const &simHit = *simHitp;
          if (std::abs(simHit.particleType()) != 13) {
            continue;
          }
          
//           if (std::abs(simHit.localPosition().z()) > 1e-9) {
//             continue;
//           }
          
//       for (auto const& simhith : simHits) {
//         for (const PSimHit& simHit : *simhith) {
          
          const GeomDet *detectorG = globalGeometry->idToDet(simHit.detUnitId());
          
          bool isbarrel = detectorG->subDetector() == GeomDetEnumerators::PixelBarrel || detectorG->subDetector() == GeomDetEnumerators::TIB || detectorG->subDetector() == GeomDetEnumerators::TOB;
          
          float absz = std::abs(detectorG->surface().toGlobal(LocalVector(0.,0.,1.)).z());
          
          bool idealdisk = absz == 1.;
          bool idealbarrel = absz<1e-9;
          
  //         idealdisk = false;
  //         idealbarrel = false;
          
          bool isstereo = trackerTopology->isStereo(simHit.detUnitId());

          std::cout << "isbarrel: " << isbarrel << " idealbarrel: " << idealbarrel << " idealdisk: " << idealdisk << "stereo: " << isstereo << " globalpos: " << detectorG->surface().position() << std::endl;
          
          LocalPoint proplocal(0.,0.,0.);
          
//           auto const propresult = fPropagator->geometricalPropagator().propagateWithPath(fts, detectorG->surface());
//           if (propresult.first.isValid()) {
//             proplocal = propresult.first.localPosition();
//           }
          
          if (!tsos.isValid()) {
            tsos = fPropagator->geometricalPropagator().propagate(fts, detectorG->surface());
          }
          else {
            tsos = fPropagator->geometricalPropagator().propagate(tsos, detectorG->surface());
          }
          
          if (tsos.isValid()) {
            proplocal = tsos.localPosition();
          }
          
          
          
          
  //         Vector3d refprop;
          
  //         LocalTrajectoryParameters
          Point3DBase<double, LocalTag> reflocal(0, 0., 0.);

          simtestz = detectorG->surface().position().z();
          simtestrho = detectorG->surface().position().perp();
          

          GlobalPoint refglobal;
          
          auto const simhitglobal = detectorG->surface().toGlobal(Point3DBase<double, LocalTag>(simHit.localPosition().x(),
                                                                                                simHit.localPosition().y(),
                                                                                                simHit.localPosition().z()));
          
          const Vector3d Msim(simhitglobal.x(), simhitglobal.y(), simhitglobal.z());
          
          auto const propglobal = detectorG->surface().toGlobal(Point3DBase<double, LocalTag>(proplocal.x(),
                                                                                                proplocal.y(),
                                                                                                proplocal.z()));
          
          
          const Vector3d Mprop(propglobal.x(), propglobal.y(), propglobal.z());
          
          Vector3d M(genmuon->vertex().x(),
                                  genmuon->vertex().y(),
                                  genmuon->vertex().z());
          
          Vector3d P(genmuon->momentum().x(),
                                  genmuon->momentum().y(),
                                  genmuon->momentum().z());
          
          
          
          
  //         if (true) {
          for (unsigned int iref=0; iref<1; ++iref) {
            const double zs = detectorG->surface().position().z();
            
            const Vector3d T0 = P.normalized();
            
  //           const Vector3d T0 = P.normalized();
            
            const Vector3d H(0.,0.,1.);
            
            const double rho = fts.transverseCurvature();
            
            double s;
            
            if (idealdisk) {
              s = (zs - M[2])/T0[2];
            }
            else if (idealbarrel) {
              HelixBarrelPlaneCrossingByCircle crossing(GlobalPoint(M[0],M[1],M[2]), GlobalVector(P[0],P[1],P[2]), rho);
              s = crossing.pathLength(detectorG->surface()).second;
            }
            else {
              HelixArbitraryPlaneCrossing crossing(Basic3DVector<float>(M[0],M[1],M[2]), Basic3DVector<float>(P[0],P[1],P[2]), rho);
              s = crossing.pathLength(detectorG->surface()).second;
  //             s = propresult.second;
            }
            
            const Vector3d HcrossT = H.cross(T0);
            const double alpha = HcrossT.norm();
            const Vector3d N0 = HcrossT.normalized();
            
            const double gamma = T0[2];
            const double q = genmuon->charge();
            const double Q = -3.8*2.99792458e-3*q/P.norm();
            const double theta = Q*s;
            
            const Vector3d dM = gamma*(theta-std::sin(theta))/Q*H + std::sin(theta)/Q*T0 + alpha*(1.-std::cos(theta))/Q*N0;
            M = M + dM;
            const Vector3d dT = gamma*(1.-std::cos(theta))*H + std::cos(theta)*T0 + alpha*std::sin(theta)*N0;
            const Vector3d T = T0 + dT;
            const double pmag = P.norm();
            P = pmag*T;
            
            refglobal = GlobalPoint(M[0], M[1], M[2]);
            reflocal = detectorG->surface().toLocal(Point3DBase<double, GlobalTag>(M[0], M[1], M[2]));
            simtestzlocalref = reflocal.z();
            
            const Vector3d xhat = Vector3d(0.,0.,1.).cross(M).normalized();
            
//             const LocalVector localxhat = LocalVector(1.,0.,0.);
//             const GlobalVector globalxhat = detectorG->surface().toGlobal(localxhat);
//             const Vector3d xhat(globalxhat.x(), globalxhat.y(), globalxhat.z());
            
            const double dx = xhat.dot(Msim-M);
            const double dxrec = xhat.dot(Mprop - Msim);
            
            const Vector3d yhat = Vector3d(0.,0.,1.).cross(xhat).normalized();
            
            const double dy = yhat.dot(Msim-M);
            const double dyrec = yhat.dot(Mprop - Msim);
            
            simtestdx = dx;
//             simtestdxrec = dxrec;
            simtestdxrec = proplocal.x() - simHit.localPosition().x();
            simtestdxprop = xhat.dot(Mprop-M);
            
            simtestdy = dy;
//             simtestdyrec = dyrec;
            simtestdyrec = proplocal.y() - simHit.localPosition().y();
            simtestdyprop = yhat.dot(Mprop-M);
            
            simtestdetid = detectorG->geographicalId().rawId();
            
            if (idealdisk) {
              break;
            }
            
  //           refprop = M;
            
  //           const Vector3d Mprop(updtsosnomat.globalPosition().x(),
  //                                 updtsosnomat.globalPosition().y(),
  //                                 updtsosnomat.globalPosition().z()); 
          }
          
          tree->Fill();
  //         else {
  //           const TrajectoryStateOnSurface propresult = fAnalyticPropagator->geometricalPropagator().propagate(fts, detectorG->surface());
  //           if (propresult.isValid()) {
  //             reflocal = propresult.localPosition();
  // //             refprop << propresult.globalPosition().x(), propresult.globalPosition().y(), propresult.globalPosition().z();
  //           }
  // //           else {
  // //             refprop << 0., 0., 0.;
  // //           }
  //         }
          
          
  //         const LocalPoint reflocal = detectorG->surface().toLocal(GlobalPoint(refprop[0], refprop[1], refprop[2]));
          

          
          const LocalPoint simlocal = simHit.localPosition();
          
//           std::cout << "isbarrel: " << isbarrel << " idealbarrel: " << idealbarrel << " idealdisk: " << idealdisk << "stereo: " << isstereo << " globalpos: " << detectorG->surface().position() << std::endl;
          std::cout << "detid: " << simHit.detUnitId() << std::endl;
          std::cout << "local z to global: " << detectorG->surface().toGlobal(LocalVector(0.,0.,1.)) << std::endl;
          std::cout << "ref      : " << reflocal << std::endl;
          std::cout << "proplocal: " << proplocal << std::endl;
          std::cout << "simlocal: " << simlocal << std::endl;
          std::cout << "sim entry point: " << simHit.entryPoint() << std::endl;
          std::cout << "sim exit point: " << simHit.exitPoint() << std::endl;
          std::cout << "refglobal: " << refglobal << std::endl;
          std::cout << "propglobal: " << propglobal << std::endl;
          std::cout << "simglobal: " << simhitglobal << std::endl;
          std::cout << "sim-ref : " << simlocal - reflocal << std::endl;
          std::cout << "sim-ref (global): " << simhitglobal - refglobal << std::endl;
          std::cout << "prop-ref: " << proplocal - reflocal << std::endl;
          std::cout << "sim-prop: " << simlocal - proplocal << std::endl;
          std::cout << "sim-prop (global): " << simhitglobal - propglobal << std::endl;
          
          
          std::cout << "simtestdx: " << simtestdx << std::endl;
          std::cout << "simtestdy: " << simtestdy << std::endl;
          
          std::cout << "simtestdxrec: " << simtestdxrec << std::endl;
          std::cout << "simtestdyrec: " << simtestdyrec << std::endl;
          
          std::cout << "simtestdxprop: " << simtestdxprop << std::endl;
          std::cout << "simtestdyprop: " << simtestdyprop << std::endl;
          
//           assert(std::abs(simtestdy)<0.1e-4);
          
//           assert(simHit.entryPoint().z()*simHit.exitPoint().z() < 0.);
          
//           assert(std::abs(simlocal.z())<1e-4);
          
//           assert(std::abs(simtestdxrec) < 40e-4);
//           assert(simtestdxprop > -950e-4);
          
          ++ihit;
          
  //         if (simHit.detUnitId() == preciseHit->geographicalId()) {                      
  //           dxsimgen.push_back(simHit.localPosition().x() - updtsos.localPosition().x());
  //           dysimgen.push_back(simHit.localPosition().y() - updtsos.localPosition().y());
  //           
  //           dxrecsim.push_back(preciseHit->localPosition().x() - simHit.localPosition().x());
  //           dyrecsim.push_back(-99.);
  //           
  //           simvalid = true;
  //           break;
  //         }
//         }
      }
      
      
    }
    return;
    
  }
  
//   TkClonerImpl hitCloner;
//   TKCloner const* cloner = static_cast<TkTransientTrackingRecHitBuilder const *>(builder)->cloner()
//   TrajectoryStateCombiner combiner;
  
  run = iEvent.run();
  lumi = iEvent.luminosityBlock();
  event = iEvent.id().event();
  
//   for (const reco::Track &track : *trackOrigH) {
  for (unsigned int itrack = 0; itrack < trackOrigH->size(); ++itrack) {
    const reco::Track &track = (*trackOrigH)[itrack];
    const reco::TrackRef trackref(trackOrigH, itrack);
//     const Trajectory& traj = (*trajH)[itraj];
    
//     const edm::Ref<std::vector<Trajectory> > trajref(trajH, j);
//     const reco::Track& track = *(*trackH)[trajref];
//     const reco::Track& track = (*trackH)[itraj];
//     const reco::Track& trackOrig = (*trackOrigH)[(*indicesH)[j]];

//     std::cout << "j " << j << " (*indicesH)[j] " << (*indicesH)[j] <<std::endl;
    
    if (track.isLooper()) {
      continue;
    }
    
    
    trackPt = track.pt();
    trackEta = track.eta();
    trackPhi = track.phi();
    trackCharge = track.charge();
    trackPtErr = track.ptError();
    
    normalizedChi2 = track.normalizedChi2();
    
//     std::cout << "track pt: " << trackPt << " track eta: " << trackEta << " trackCharge: " << trackCharge << " qop: " << track.parameters()[0] << std::endl;
    
    auto const& tkparms = track.parameters();
    auto const& tkcov = track.covariance();
    trackParms.fill(0.);
    trackCov.fill(0.);
    //use eigen to fill raw memory
    Map<Vector5f>(trackParms.data()) = Map<const Vector5d>(tkparms.Array()).cast<float>();
    Map<Matrix<float, 5, 5, RowMajor> >(trackCov.data()).triangularView<Upper>() = Map<const Matrix<double, 5, 5, RowMajor> >(tkcov.Array()).cast<float>().triangularView<Upper>();
    
//     std::cout << "track charge: " << track.charge() << " trackorig charge " << trackOrig.charge() << "inner state charge " << tms.back().updatedState().charge() << std::endl;
    
    const reco::GenParticle* genpart = nullptr;
    
    genPt = -99.;
    genEta = -99.;
    genPhi = -99.;
    genCharge = -99;
    genX = -99.;
    genY = -99.;
    genZ = -99.;
    genParms.fill(0.);
    
    int genBarcode = -99;
    
    
    if (doGen_) {
//       bool isjpsi = false;
//       for (std::vector<reco::GenParticle>::const_iterator g = genPartCollection->begin(); g != genPartCollection->end(); ++g)
//       {
//         if (std::abs(g->pdgId()) == 443) {
//           isjpsi = true;
//           break;
//         }
//       }
      
      float drmin = 0.1;
      
      for (std::vector<reco::GenParticle>::const_iterator g = genPartCollection->begin(); g != genPartCollection->end(); ++g)
      {
        if (g->status() != 1) {
          continue;
        }
        if (std::abs(g->pdgId()) != 13) {
          continue;
        }
        
//         if (isjpsi && g->charge() != track.charge()) {
//           continue;
//         }
        
        float dR = deltaR(*g, track);
        
        if (dR < drmin)
        {
          drmin = dR;
          
          genpart = &(*g);
          
          genBarcode = (*genPartBarcodes)[g - genPartCollection->begin()];
          
          genPt = g->pt();
          genEta = g->eta();
          genPhi = g->phi();
          genCharge = g->charge();
          
          genX = g->vertex().x();
          genY = g->vertex().y();
          genZ = g->vertex().z();
          
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

    if (doMuons_) {
      for (auto const &muon : *muons) {
        if (muon.innerTrack() == trackref) {
          muonPt = muon.pt();
          muonLoose = muon.passed(reco::Muon::CutBasedIdLoose);
          muonMedium = muon.passed(reco::Muon::CutBasedIdMedium);
          muonTight = muon.passed(reco::Muon::CutBasedIdTight);
          muonIsPF = muon.isPFMuon();
          muonIsTracker = muon.isTrackerMuon();
          muonIsGlobal = muon.isGlobalMuon();
          muonIsStandalone = muon.isStandAloneMuon();
          if (muon.muonBestTrack() == trackref) {
            muonInnerTrackBest = true;
          }
        } 
      }
      if (assoc->contains(track.extra().id())) {
        const reco::TrackExtraRef &trackextraref = (*assoc)[track.extra()];
        if (trackextraref.isNonnull()) {
          trackExtraAssoc = true;
        }
      }
    }
    
    

//     PropagationDirection rpropdir = traj.direction();
//     PropagationDirection fpropdir = rpropdir == alongMomentum ? oppositeToMomentum : alongMomentum;
    
    //TODO properly handle the outside-in case
    assert(track.seedDirection() == alongMomentum);
    
    //prepare hits
    TransientTrackingRecHit::RecHitContainer hits;
    hits.reserve(track.recHitsSize());
//     hits.reserve(track.recHitsSize()+1);
//     hits.push_back(RecHitPointer(InvalidTrackingRecHit());
    
//     printf("building hits\n");
    
    std::set<std::array<int, 3>> hitlayers;
    
    for (auto it = track.recHitsBegin(); it != track.recHitsEnd(); ++it) {
      const GeomDet* detectorG = globalGeometry->idToDet((*it)->geographicalId());
      const GluedGeomDet* detglued = dynamic_cast<const GluedGeomDet*>(detectorG);
      
      // split matched invalid hits
      if (detglued != nullptr && !(*it)->isValid()) {
        bool order = detglued->stereoDet()->surface().position().mag() > detglued->monoDet()->surface().position().mag();
        const GeomDetUnit* detinner = order ? detglued->monoDet() : detglued->stereoDet();
        const GeomDetUnit* detouter = order ? detglued->stereoDet() : detglued->monoDet();
        
        hits.push_back(TrackingRecHit::RecHitPointer(new InvalidTrackingRecHit(*detinner, (*it)->type())));
        hits.push_back(TrackingRecHit::RecHitPointer(new InvalidTrackingRecHit(*detouter, (*it)->type())));
        
//         auto const& layerinner = detidlayermap.at(detinner->geographicalId());
//         auto const& layerouter = detidlayermap.at(detouter->geographicalId());
        
//         if (!hitlayers.count(layerinner))
//           hits.push_back(TrackingRecHit::RecHitPointer(new InvalidTrackingRecHit(*detinner, (*it)->type())));
//         hitlayers.insert(layerinner);
//         if (!hitlayers.count(layerouter))
//           hits.push_back(TrackingRecHit::RecHitPointer(new InvalidTrackingRecHit(*detouter, (*it)->type())));
//         hitlayers.insert(layerouter);
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
//             hitquality = !pixhit->isOnEdge();
//             hitquality = cluster.sizeX() > 1;
            
//             hitquality = pixhit->hasFilledProb() && pixhit->clusterProbability(0) > 0.000125 && pixhit->qBin()>0 && pixhit->qBin()<4;
            
//             if (pixhit->hasFilledProb() && pixhit->clusterProbability(0) > 0.000125 && pixhit->qBin()>0 && pixhit->qBin()<4) {
//             if (pixhit->hasFilledProb() && pixhit->clusterProbability(0) > 0.000125 && pixhit->qBin()>0 && pixhit->qBin()<3) {
//               hitquality = true;
//             }
          }
          else {
            assert(tkhit->cluster_strip().isNonnull());
            const SiStripCluster& cluster = *tkhit->cluster_strip();
            const StripTopology* striptopology = dynamic_cast<const StripTopology*>(&(detectorG->topology()));
            assert(striptopology);
            
            const uint16_t firstStrip = cluster.firstStrip();
            const uint16_t lastStrip = cluster.firstStrip() + cluster.amplitudes().size() - 1;
            const bool isOnEdge = firstStrip == 0 || lastStrip == (striptopology->nstrips() - 1);
            
//             if (isOnEdge) {
//               std::cout << "strip hit isOnEdge" << std::endl;
//             }
            
//             hitquality = !isOnEdge;
            hitquality = true;
            
            
            
            
//             SiStripClusterInfo clusterInfo = SiStripClusterInfo(cluster, iSetup, (*it)->geographicalId().rawId());
//             if (clusterInfo.signalOverNoise() > 12.) {
//               hitquality = true;
//             }
//             hitquality = true;
          }
          
        }
        else {
          hitquality = true;
        }
        
//         std::cout << "detid: " << (*it)->geographicalId().rawId() << std::endl;
        
//         bool foundsim = false;
//         for (auto const& simhith : simHits) {
//           for (const PSimHit& simHit : *simhith) {
//             if (std::abs(simHit.particleType()) == 13 && simHit.detUnitId() == (*it)->geographicalId()) {
//               foundsim = true;
//               break;
//             }
//           }
//           if (foundsim) {
//             break;
//           }
//         }
//         
//         if (!foundsim) {
//           hitquality = false;
//         }
        
//         for (auto const& simhith : simHits) {
//           for (const PSimHit& simHit : *simhith) {
//             if (simHit.detUnitId() == (*it)->geographicalId()) {
// //               std::cout << "particle type: " << simHit.particleType() << std::endl;
// //               std::cout << "track id: " << simHit.trackId() << std::endl;
// //               std::cout << "entry point: " << simHit.entryPoint() << std::endl;
// //               std::cout << "exit point: " << simHit.exitPoint() << std::endl;
// //               std::cout << "local position: " << simHit.localPosition() << std::endl;
//               if (std::abs(simHit.particleType()) == 13 && std::abs(simHit.localPosition().z()) > 1e-4) {
//                 std::cout << "anomalous simhit!" << std::endl;
//                 hitquality = false;
//               }
//             }
//           }
//         }


        if (hitquality) {
          hits.push_back((*it)->cloneForFit(*detectorG));
        }
        else {
          hits.push_back(TrackingRecHit::RecHitPointer(new InvalidTrackingRecHit(*detectorG, TrackingRecHit::inactive)));
        }

        
//         auto const &layer = detidlayermap.at((*it)->geographicalId());
//         
//         if (!hitlayers.count(layer)) {
//           if (hitquality) {
//             hits.push_back((*it)->cloneForFit(*detectorG));
//           }
//           else {
//             hits.push_back(TrackingRecHit::RecHitPointer(new InvalidTrackingRecHit(*detectorG, TrackingRecHit::inactive)));
//           }
//         }
//         hitlayers.insert(layer);
      }
    }
    
//     printf("done building hits\n");
    
//     const unsigned int nhits = track.recHitsSize();
    const unsigned int nhits = hits.size();
    nHits = nhits;
//     unsigned int npixhits = 0;

    unsigned int nvalid = 0;
    unsigned int nvalidpixel = 0;
    unsigned int nvalidalign2d = 0;
    
//     std::set<std::array<int, 3>> hitlayers;
    
    // count valid hits since this is needed to size the arrays
    for (auto const& hit : hits) {
//       auto const &layers = detidlayermap.at(hit->geographicalId());
//       if (hitlayers.count(layers)) {
//         std::cout << "WARNING: multiple hits on the same layer!!!" << std::endl;
//         std::cout << layers[0] << " " << layers[1] << " " << layers[2] << std::endl;
//       }
//       hitlayers.insert(layers);
      
      assert(hit->dimension()<=2);
      if (hit->isValid()) {
        nvalid += 1;
        
//         const uint32_t gluedid = trackerTopology->glued(hit->geographicalId());
//         const bool isglued = gluedid != 0;
//         const DetId parmdetid = isglued ? DetId(gluedid) : hit->geographicalId();
//         const bool align2d = detidparms.count(std::make_pair(1, parmdetid));
//         const bool align2d = detidparms.count(std::make_pair(2, hit->geographicalId()));
        
        const bool align2d = detidparms.count(std::make_pair(1, hit->geographicalId()));
//         
        if (align2d) {
          nvalidalign2d += 1;
        }
        if (GeomDetEnumerators::isTrackerPixel(hit->det()->subDetector())) {
          nvalidpixel += 1;
        }
      }
    }
    
//     //count valid hits since this is needed to size the arrays
//     auto const& hitsbegin = track.recHitsBegin();
//     for (unsigned int ihit = 0; ihit < track.recHitsSize(); ++ihit) {
//       auto const& hit = *(hitsbegin + ihit);
//       if (hit->isValid() && hit->dimension()<=2) {
//         nvalid += 1;
//         
//         const GeomDet *detectorG = globalGeometry->idToDet(hit->geographicalId());
//         if (hit->dimension()==2 && GeomDetEnumerators::isTrackerPixel(detectorG->subDetector())) {
// //         if (hit->dimension()==2) {
//           nvalidpixel += 1;
//         }
//       }
//     }
    
    nValidHits = nvalid;
    nValidPixelHits = nvalidpixel;
    
    nValidHitsFinal = 0;
    nValidPixelHitsFinal = 0;
    
//     const unsigned int nstriphits = nhits-npixhits;
//     const unsigned int nparsAlignment = nstriphits + 2*npixhits;
//     const unsigned int nvalidstrip = nvalid - nvalidpixel;
//     const unsigned int nparsAlignment = nvalidstrip + 2*nvalidpixel;
    const unsigned int nparsAlignment = 2*nvalid + nvalidalign2d;
//     const unsigned int nparsAlignment = 6*nvalid;
    const unsigned int nparsBfield = nhits;
    const unsigned int nparsEloss = nhits - 1;
    const unsigned int npars = nparsAlignment + nparsBfield + nparsEloss;
    
//     const unsigned int nstateparms = 5*(nhits+1);
    const unsigned int nstateparms = 3*(nhits+1) - 1;
//     const unsigned int nstateparms = 3*nhits - 1;
    const unsigned int nparmsfull = nstateparms + npars;
    
    
    const unsigned int nstateparmspost = 5*(nhits+1);
    
//     std::cout << "nhits " << nhits << std::endl;
//     std::cout << "nstateparms " << nstateparms << std::endl;
//     std::cout << "nparmsfull " << nparmsfull << std::endl;
//     std::cout << "nparmsfull " << nparmsfull << std::endl;
//     std::cout << "nparmsfull " << nparmsfull << std::endl;
//     std::cout << "nparmsfull " << nparmsfull << std::endl;
//     std::cout << "nparmsfull " << nparmsfull << std::endl;
    
//     const unsigned int npropparms = 5*(nhits-1);
//     const unsigned int nhitparms = 2*nhits;
//     const unsigned int nmomparms = 3*(nhits-1);
//     const unsigned int nposparms = 2*(nhits-1);
//     constexpr unsigned int nrefparms = 5;
    

    
    //active double for autodiff gradients
//     using Adouble = AutoDiffScalar<VectorXd>;
//     using AVectorXd = Matrix<Adouble, Dynamic, 1>;
//     //double double for autodiff hessians
//     using AAdouble = AutoDiffScalar<AVectorXd>;
    

    
//     using AAXd = AANT<double, Dynamic>;
//     using AAdouble = AAXd;
//     
//     using AA2d = AANT<double, 2>;
//     using AA3d = AANT<double, 3>;
//     using AA4d = AANT<double, 4>;
//     using AA12d = AANT<double, 12>;
//     
//     using ScalarConst = AANT<double, 0>;
    
//     using AConstd = AutoDiffScalar<VectorXd>;
//     using AConstd = AutoDiffScalar<Matrix<double, 0, 0>>;
    
    
//     using VectorXAd = Matrix<AScalar, Dynamic, 1>;
//     using MatrixXAd = Matrix<AScalar, Dynamic, Dynamic>;
    
    //two position parameters and and one alignment parameter
    using StripHitScalar = AANT<double, 3>;;
    
    using StripHit1DJacobian = Matrix<StripHitScalar, 1, 2>;
    
    using StripHitVector = Matrix<StripHitScalar, 2, 1>;
    using StripHit2DCovariance = Matrix<StripHitScalar, 2, 2>;
    using StripHit2DJacobian = Matrix<StripHitScalar, 2, 2>;

    
    
    //two hit dimensions and two alignment parameters
    using PixelHit2DScalar = AANT<double, 4>;
    using PixelHit2DVector = Matrix<PixelHit2DScalar, 2, 1>;
    using PixelHit2DCovariance = Matrix<PixelHit2DScalar, 2, 2>;
    using PixelHit2DJacobian = Matrix<PixelHit2DScalar, 2, 2>;
    
    
    //2x5 state parameters, one bfield parameter, and one material parameter
//     using MSScalar = AANT<double, 11>;;
    using MSScalar = AANT<double, 13>;
    using MSVector = Matrix<MSScalar, 5, 1>;
    using MSProjection = Matrix<MSScalar, 5, 5>;
    using MSJacobian = Matrix<MSScalar, 5, 5>;
    using MSCovariance = Matrix<MSScalar, 5, 5>;

    using BSScalar = AANT<double, 2>;
    
//     using HitProjection = Matrix<AAdouble, 2, 5>;
//     using HitCovariance = Matrix<AAdouble, 2, 2>;
//     using HitVector = Matrix<AAdouble, 2, 1>;
    
//     evector<HitCovarianceMatrix> Vinv(nhits, HitCovarianceMatrix::Zero());
//     evector<HitProjection> Hh(nhits, HitProjection::Zero());
//     evector<HitVector> dy0(nhits, HitVector::Zero());
//     evector<StateVector> dx(nhits, StateVector::Zero());
//     //initialize backpropagation indexing
//     for (unsigned int i=0; i<nhits; ++i) {
// //       StateVector& dxi = dx[i];
// //       dxi.derivatives().resize(nstateparms);
//       for (unsigned int j=0; j<5; ++j) {
// //         dx[i][j].derivatives() = VectorXd::Unit(nstateparms, 5*i + j);
//         init_twice_active_var(dx[i][j], nstateparms, 5*i +j);
//       }
//     }
    
    
    VectorXd gradfull;
    MatrixXd hessfull;
    
    MatrixXd statejac;
    VectorXd dxstate;
    
    MatrixXd validdxeigjac;
    VectorXd validdxeig;
    
    Matrix<float, Dynamic, 2> rxfull(nvalid, 2); 
    Matrix<float, Dynamic, 2> ryfull(nvalid, 2);
    
//     VectorXd gradfull = chisq.value().derivatives();
//     MatrixXd hessfull = MatrixXd::Zero(nparmsfull, nparmsfull);
//     for (unsigned int i=0; i<nstateparms; ++i) {
//       hessfull.row(i) = chisq.derivatives()[i].derivatives();
//     }
    
    
    globalidxv.clear();
    globalidxv.resize(npars, 0);
    
    nParms = npars;
    if (fillTrackTree_) {
      tree->SetBranchAddress("globalidxv", globalidxv.data());
    }
    
//     TrajectoryStateOnSurface currtsos;
    
    
    
    VectorXd dxfull;
    MatrixXd dxdparms;
    VectorXd grad;
    MatrixXd hess;
    LDLT<MatrixXd> Cinvd;
    
    if (dogen && genpart==nullptr) {
      continue;
    }
    
//     if (dogen && genpart->eta()>-2.3) {
//       continue;
//     }

//     if (genpart==nullptr) {
//       continue;
//     }
//     if (genpart->pt()>10.) {
//       continue;
//     }
//     if (genpart->pt()<100.) {
//       continue;
//     }
//     if (genpart->eta()>-2.3) {
//       continue;
//     }
    
    if (debugprintout_) {
      std::cout << "initial reference point parameters:" << std::endl;
      std::cout << track.parameters() << std::endl;
    }

//     //prepare hits
//     TransientTrackingRecHit::RecHitContainer hits;
//     hits.reserve(track.recHitsSize());
//     for (auto it = track.recHitsBegin(); it != track.recHitsEnd(); ++it) {
//       const GeomDet *detectorG = globalGeometry->idToDet((*it)->geographicalId());
//       hits.push_back((*it)->cloneForFit(*detectorG));
//     }
    
//     // fix mixed up clusters?
//     for (unsigned int ihit=0; ihit<(hits.size()-1); ++ihit) {
//       TrackerSingleRecHit* hit = const_cast<TrackerSingleRecHit*>(dynamic_cast<const TrackerSingleRecHit*>(hits[ihit].get()));
//       TrackerSingleRecHit* nexthit = const_cast<TrackerSingleRecHit*>(dynamic_cast<const TrackerSingleRecHit*>(hits[ihit+1].get()));
//       
// //      const TrackingRecHitSingle* nexthit = hits[ihit+1];
//       
//       if (!hit || !nexthit) {
//         continue;
//       }
//       
//       const DetId partnerid = trackerTopology->partnerDetId(hit->geographicalId());
// //       
//       if (partnerid == nexthit->geographicalId()) {
// //         std::cout << "swapping clusters" << std::endl;
//         const OmniClusterRef::ClusterStripRef cluster = hit->cluster_strip();
//         const OmniClusterRef::ClusterStripRef nextcluster = nexthit->cluster_strip();
//         
//         hit->setClusterStripRef(nextcluster);
//         nexthit->setClusterStripRef(cluster);
//       }
// 
//     }

//     if (genpart==nullptr) {
//       continue;
//     }
//     
//     if (genpart->eta()<-2.4 || genpart->eta()>-2.3) {
//       continue;
//     }
    
    
    FreeTrajectoryState refFts;
    
    if (dogen) {
//     if (true) {
      
      if (genpart==nullptr) {
        continue;
      }
      //init from gen state
      auto const& refpoint = genpart->vertex();
      auto const& trackmom = genpart->momentum();
      const GlobalPoint refpos(refpoint.x(), refpoint.y(), refpoint.z());
      const GlobalVector refmom(trackmom.x(), trackmom.y(), trackmom.z()); 
      const GlobalTrajectoryParameters refglobal(refpos, refmom, genpart->charge(), field);
      
//       std::cout << "gen ref state" << std::endl;
//       std::cout << refpos << std::endl;
//       std::cout << refmom << std::endl;
//       std::cout << genpart->charge() << std::endl;
      
      //zero uncertainty on generated parameters
//       AlgebraicSymMatrix55 nullerr;
//       const CurvilinearTrajectoryError referr(nullerr);
      
      refFts = FreeTrajectoryState(refpos, refmom, genpart->charge(), field);
//       refFts = FreeTrajectoryState(refglobal, referr);
    }
    else {
      //init from track state
      auto const& refpoint = track.referencePoint();
      auto const& trackmom = track.momentum();
      const GlobalPoint refpos(refpoint.x(), refpoint.y(), refpoint.z());
      const GlobalVector refmom(trackmom.x(), trackmom.y(), trackmom.z()); 
//       const GlobalTrajectoryParameters refglobal(refpos, refmom, track.charge(), field);
//       const CurvilinearTrajectoryError referr(track.covariance());
      
      //null uncertainty (tracking process noise sum only)
//       AlgebraicSymMatrix55 nullerr;
//       const CurvilinearTrajectoryError referr(nullerr);
      
      refFts = FreeTrajectoryState(refpos, refmom, track.charge(), field);
//       refFts = FreeTrajectoryState(refglobal, referr);
    }

//     std::vector<std::pair<TrajectoryStateOnSurface, double>> layerStates;
    std::vector<TrajectoryStateOnSurface> layerStates;
    layerStates.reserve(nhits);
    
    bool valid = true;

    
//     //do propagation and prepare states
//     auto propresult = fPropagator->propagateWithPath(refFts, *hits.front()->surface());
//     if (!propresult.first.isValid()) {
//       std::cout << "Abort: Propagation from reference point failed" << std::endl;
//       continue;
//     }
//     layerStates.push_back(propresult);
//     
//     for (auto const& hit : hits) {
//       propresult = fPropagator->propagateWithPath(layerStates.back().first, *hit->surface());
//       if (!propresult.first.isValid()) {
//         std::cout << "Abort: Propagation failed" << std::endl;
//         valid = false;
//         break;
//       }
//       layerStates.push_back(propresult);
//     }
//     
//     if (!valid) {
//       continue;
//     }
    

    
    //inflate errors
//     refFts.rescaleError(100.);
    
    
//     unsigned int ntotalhitdim = 0;
//     unsigned int alignmentidx = 0;
//     unsigned int bfieldidx = 0;
//     unsigned int elossidx = 0;
    
//     constexpr unsigned int niters = 1;
    constexpr unsigned int niters = 10;
    
    for (unsigned int iiter=0; iiter<niters; ++iiter) {
      if (debugprintout_) {
        std::cout<< "iter " << iiter << std::endl;
      }
      
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

      }      
      
        
      const bool islikelihood = false;
//       const bool islikelihood = iiter > 0;
//       const bool islikelihood = true;
      
      gradfull = VectorXd::Zero(nparmsfull);
      hessfull = MatrixXd::Zero(nparmsfull, nparmsfull);
      statejac = MatrixXd::Zero(nstateparmspost, nparmsfull);
      
      validdxeigjac = MatrixXd::Zero(2*nvalid, nstateparms);
      
      evector<std::array<Matrix<double, 8, 8>, 11> > dhessv;
      if (islikelihood) {
        dhessv.resize(nhits-1);
      }
      
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
      
      
      
      const uint32_t gluedid0 = trackerTopology->glued(hits[0]->det()->geographicalId());
      const bool isglued0 = gluedid0 != 0;
      const DetId parmdetid0 = isglued0 ? DetId(gluedid0) : hits[0]->geographicalId();
      const unsigned int bfieldidx = detidparms.at(std::make_pair(6, parmdetid0));
      fieldOffset->setOffset(corparms_[bfieldidx]);
      
//       dhessv.reserve(nhits-1);
      
      unsigned int parmidx = 0;
      unsigned int alignmentparmidx = 0;
      unsigned int ivalidhit = 0;

      if (iiter > 0) {
        //update current state from reference point state (errors not needed beyond first iteration)
        JacobianCurvilinearToCartesian curv2cart(refFts.parameters());
        const AlgebraicMatrix65& jac = curv2cart.jacobian();
        const AlgebraicVector6 glob = refFts.parameters().vector();
        
        auto const& dxlocal = dxstate.head<5>();
        const Matrix<double, 6, 1> globupd = Map<const Matrix<double, 6, 1>>(glob.Array()) + Map<const Matrix<double, 6, 5, RowMajor>>(jac.Array())*dxlocal;
        
//         const Vector6d dglob = Map<const Matrix<double, 6, 5, RowMajor>>(jac.Array())*dxlocal;
//         
//         std::cout << "iiter = " << iiter << std::endl;
//         std::cout << "dxlocal "  << dxlocal << std::endl;
//         std::cout << "glob " << glob << std::endl;
//         std::cout << "dglob " << dglob << std::endl;
        
        const GlobalPoint pos(globupd[0], globupd[1], globupd[2]);
        const GlobalVector mom(globupd[3], globupd[4], globupd[5]);
        const double charge = std::copysign(1., refFts.charge()/refFts.momentum().mag() + dxlocal[0]);
//         std::cout << "before update: reffts:" << std::endl;
//         std::cout << refFts.parameters().vector() << std::endl;
//         std::cout << "charge " << refFts.charge() << std::endl;
        refFts = FreeTrajectoryState(pos, mom, charge, field);
//         std::cout << "after update: reffts:" << std::endl;
//         std::cout << refFts.parameters().vector() << std::endl;
//         std::cout << "charge " << refFts.charge() << std::endl;
//         currentFts = refFts;
      }
      
//       Matrix5d Hlm = Matrix5d::Identity();
//       currentFts = refFts;
//       TrajectoryStateOnSurface currentTsos;

//       ;
      
//       auto const &surface0orig = *hits[0]->surface();
      auto const &surface0 = *surfacemap_.at(hits[0]->geographicalId());
//       printf("detid = %u, parmdetid = %u, old = %p, new  = %p\n", hits[0]->geographicalId().rawId(), parmdetid0.rawId(), &surface0orig, &surface0);
//       std::cout << "old xi = " << surface0orig.mediumProperties().xi() << " new xi = " << surface0.mediumProperties().xi() << " dxi = " << surface0.mediumProperties().xi() - surface0orig.mediumProperties().xi() << std::endl;
      auto propresult = fPropagator->geometricalPropagator().propagateWithPath(refFts, surface0);
//       auto propresult = fPropagator->geometricalPropagator().propagateWithPath(refFts, *hits[0]->surface());
//       auto propresult = fPropagator->geometricalPropagator().propagateWithPath(refFts, *beampipe);
      if (!propresult.first.isValid()) {
        std::cout << "Abort: Propagation of reference state Failed!" << std::endl;
        valid = false;
        break;
      }
      
//       std::cout << "position on beampipe " << propresult.first.globalParameters().position() << std::endl;
      
      const Matrix<double, 5, 6> FdFp = curv2localTransportJacobian(refFts, propresult, false);

      Matrix<double, 2, 2> J = FdFp.block<2, 2>(3, 3);
      // (du/dalphap)^-1
      Matrix<double, 2, 2> Sinv = FdFp.block<2, 2>(3, 1).inverse();
      // du/dqopp
      Matrix<double, 2, 1> D = FdFp.block<2, 1>(3, 0);
      // du/dBp
      Matrix<double, 2, 1> Bpref = FdFp.block<2, 1>(3, 5);

      constexpr unsigned int jacstateidxout = 0;
      constexpr unsigned int jacstateidxin = 0;
      
      // qop_i
      statejac(jacstateidxout, jacstateidxin + 2) = 1.;
      // d(lambda, phi)_i/dqop_i
      statejac.block<2, 1>(jacstateidxout + 1, jacstateidxin + 2) = -Sinv*D;
      // d(lambda, phi)_i/(dxy, dsz)
      statejac.block<2, 2>(jacstateidxout + 1, jacstateidxin) = -Sinv*J;
      // d(lambda, phi)_i/du_(i+1)
      statejac.block<2, 2>(jacstateidxout + 1, jacstateidxin + 3) = Sinv;
      // d(lambda, phi) / dbeta
      statejac.block<2, 1>(jacstateidxout + 1, nstateparms + parmidx) = -Sinv*Bpref;
      // dxy
      statejac(jacstateidxout + 3, jacstateidxin) = 1.;
      // dsz
      statejac(jacstateidxout + 4, jacstateidxin + 1) = 1.;
      
      if (bsConstraint_) {
        // apply beamspot constraint
        // TODO add residual corrections for beamspot parameters?
        
        constexpr unsigned int nlocalstate = 2;
        constexpr unsigned int nlocalbs = 0;
        constexpr unsigned int nlocalparms = nlocalbs;
        
        constexpr unsigned int nlocal = nlocalstate + nlocalparms;
        
        constexpr unsigned int localstateidx = 0;
        constexpr unsigned int localparmidx = localstateidx + nlocalstate;
        
        constexpr unsigned int fullstateidx = 0;
        const unsigned int fullparmidx = nstateparms + parmidx;
        
        JacobianCurvilinearToCartesian curv2cart(refFts.parameters());
        const AlgebraicMatrix65& jac = curv2cart.jacobian();
        
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
        
        Matrix<BSScalar, 2, 1> du = Matrix<BSScalar, 2, 1>::Zero();
        for (unsigned int j=0; j<du.size(); ++j) {
          init_twice_active_var(du[j], nlocal, localstateidx + j);
        }
        
        Matrix<BSScalar, 3, 1> dbs0;
        dbs0[0] = BSScalar(refFts.position().x() - x0);
        dbs0[1] = BSScalar(refFts.position().y() - y0);
        dbs0[2] = BSScalar(refFts.position().z() - z0);
        
//         std::cout << "dposition / d(qop, lambda, phi) (should be 0?):" << std::endl;
//         std::cout << Map<const Matrix<double, 6, 5, RowMajor>>(jac.Array()).topLeftCorner<3,3>() << std::endl;
        
        const Matrix<BSScalar, 3, 2> jacpos = Map<const Matrix<double, 6, 5, RowMajor>>(jac.Array()).topRightCorner<3,2>().cast<BSScalar>();
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
        gradfull.segment<nlocalparms>(fullparmidx) += gradlocal.segment<nlocalparms>(localparmidx);

        //fill global hessian (upper triangular blocks only)
        hessfull.block<nlocalstate,nlocalstate>(fullstateidx, fullstateidx) += hesslocal.topLeftCorner<nlocalstate,nlocalstate>();
        hessfull.block<nlocalstate,nlocalparms>(fullstateidx, fullparmidx) += hesslocal.topRightCorner<nlocalstate, nlocalparms>();
        hessfull.block<nlocalparms, nlocalparms>(fullparmidx, fullparmidx) += hesslocal.bottomRightCorner<nlocalparms, nlocalparms>();
        
        
      }
      
      
      Matrix<double, 5, 6> FdFm = curv2localTransportJacobian(refFts, propresult, true);
      
      for (unsigned int ihit = 0; ihit < hits.size(); ++ihit) {
//         std::cout << "ihit " << ihit << std::endl;
        auto const& hit = hits[ihit];
        
        const uint32_t gluedid = trackerTopology->glued(hit->det()->geographicalId());
        const bool isglued = gluedid != 0;
        const DetId parmdetid = isglued ? DetId(gluedid) : hit->geographicalId();
        const GeomDet* parmDet = isglued ? globalGeometry->idToDet(parmdetid) : hit->det();
        const double xifraction = isglued ? hit->det()->surface().mediumProperties().xi()/parmDet->surface().mediumProperties().xi() : 1.;
        
//         const DetId partnerid = isglued ? trackerTopology->partnerDetId(hit->det()->geographicalId()) : DetId();
//         
//         const bool isfront = ihit != (hits.size() - 1) && isglued && hits[ihit+1]->det()->geographicalId() == partnerid;
//         const bool isback = ihit !=0 && isglued && hits[ihit-1]->det()->geographicalId() == partnerid;
//         
//         const double xifraction = isfront ? 0. : 1.;
                
        TrajectoryStateOnSurface updtsos = propresult.first;
        
        //apply measurement update if applicable
//         std::cout << "constructing preciseHit" << std::endl;
        auto const& preciseHit = hit->isValid() ? cloner.makeShared(hit, updtsos) : hit;
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
        
        const bool align2d = detidparms.count(std::make_pair(1, preciseHit->geographicalId()));

        
        // compute convolution correction in local coordinates (BEFORE material effects are applied)
//         const Matrix<double, 2, 1> dxlocalconv = localPositionConvolution(updtsos);
         
        // curvilinear to local jacobian
        JacobianCurvilinearToLocal curv2localm(updtsos.surface(), updtsos.localParameters(), *updtsos.magneticField());
        const AlgebraicMatrix55& curv2localjacm = curv2localm.jacobian();
//         const Matrix<double, 5, 5> Hm = Map<const Matrix<double, 5, 5, RowMajor>>(curv2localjacm.Array()); 
        
        //energy loss jacobian
        const Matrix<double, 5, 6> EdE = materialEffectsJacobian(updtsos, fPropagator->materialEffectsUpdator());
//         const Matrix<double, 5, 6> EdE = materialEffectsJacobianVar(updtsos, fPropagator->materialEffectsUpdator());
        
//         const Matrix<double, 5, 6> EdEVar = materialEffectsJacobianVar(updtsos, fPropagator->materialEffectsUpdator());
//         
//         std::cout << "EdE:" << std::endl;
//         std::cout << EdE << std::endl;
//         std::cout << "EdEVar:" << std::endl;
//         std::cout << EdEVar << std::endl;
       
        //process noise jacobians
        const std::array<Matrix<double, 5, 5>, 5> dQs = processNoiseJacobians(updtsos, fPropagator->materialEffectsUpdator());
        
        //TODO update code to allow doing this in one step with nominal update
        //temporary tsos to extract process noise without loss of precision
        TrajectoryStateOnSurface tmptsos(updtsos);
        tmptsos.update(tmptsos.localParameters(),
                        LocalTrajectoryError(0.,0.,0.,0.,0.),
                        tmptsos.surface(),
                        tmptsos.magneticField(),
                        tmptsos.surfaceSide());
        
        //apply the state update from the material effects
        bool ok = fPropagator->materialEffectsUpdator().updateStateInPlace(tmptsos, alongMomentum);
        if (!ok) {
          std::cout << "Abort: material update failed" << std::endl;
          valid = false;
          break;
        }
        
        const AlgebraicVector5 dxeloss = tmptsos.localParameters().vector() - updtsos.localParameters().vector();
        
        // compute convolution effects
//         const AlgebraicVector5 dlocalconv = localMSConvolution(updtsos, fPropagator->materialEffectsUpdator());
        
//         const GlobalPoint updtsospos = updtsos.globalParameters().position();
//         std::cout << "before material update: " << updtsos.globalParameters().position() << " " << updtsos.globalParameters().momentum() << std::endl;
        ok = fPropagator->materialEffectsUpdator().updateStateInPlace(updtsos, alongMomentum);
        if (!ok) {
          std::cout << "Abort: material update failed" << std::endl;
          valid = false;
          break;
        }
//         std::cout << "after material update: " << updtsos.globalParameters().position() << " " << updtsos.globalParameters().momentum() << std::endl;
        
        
//         std::cout << "local parameters" << std::endl;
//         std::cout << updtsos.localParameters().vector() << std::endl;
//         std::cout << "dlocalconv" << std::endl;
//         std::cout << dlocalconv << std::endl;
//         
//         // apply convolution effects
//         const LocalTrajectoryParameters localupd(updtsos.localParameters().vector() + dlocalconv,
//                                                  updtsos.localParameters().pzSign());
//         updtsos.update(localupd,
// //                        updtsos.localError(),
//                        updtsos.surface(),
//                        updtsos.magneticField(),
//                        updtsos.surfaceSide());
        
        //get the process noise matrix
        AlgebraicMatrix55 const Qmat = tmptsos.localError().matrix();
        const Map<const Matrix<double, 5, 5, RowMajor>>Q(Qmat.Array());
        std::cout<< "Q" << std::endl;
        std::cout<< Q << std::endl;
        

        // update state from previous iteration
        //momentum kink residual
        AlgebraicVector5 idx0(0., 0., 0., 0., 0.);
        if (iiter==0) {
          layerStates.push_back(updtsos);
        }
        else {          
          //current state from previous state on this layer
          //save current parameters          
          TrajectoryStateOnSurface& oldtsos = layerStates[ihit];
          
//           JacobianCurvilinearToLocal curv2localold(oldtsos.surface(), oldtsos.localParameters(), *oldtsos.magneticField());
//           const AlgebraicMatrix55& curv2localjacold = curv2localold.jacobian();
//           const Matrix<double, 5, 5> Hold = Map<const Matrix<double, 5, 5, RowMajor>>(curv2localjacold.Array()); 
          
          const AlgebraicVector5 local = oldtsos.localParameters().vector();
          auto const& dxlocal = dxstate.segment<5>(5*(ihit+1));
          const Matrix<double, 5, 1> localupd = Map<const Matrix<double, 5, 1>>(local.Array()) + dxlocal;
          AlgebraicVector5 localvecupd(localupd[0],localupd[1],localupd[2],localupd[3],localupd[4]);
          
          idx0 = localvecupd - updtsos.localParameters().vector();
          
          const LocalTrajectoryParameters localparms(localvecupd, oldtsos.localParameters().pzSign());
          
//           std::cout << "before update: oldtsos:" << std::endl;
//           std::cout << oldtsos.localParameters().vector() << std::endl;
          oldtsos.update(localparms, oldtsos.surface(), field, oldtsos.surfaceSide());
//           std::cout << "after update: oldtsos:" << std::endl;
//           std::cout << oldtsos.localParameters().vector() << std::endl;
          updtsos = oldtsos;

        }
        
        // curvilinear to local jacobian
        JacobianCurvilinearToLocal curv2localp(updtsos.surface(), updtsos.localParameters(), *updtsos.magneticField());
        const AlgebraicMatrix55& curv2localjacp = curv2localp.jacobian();
        const Matrix<double, 5, 5> Hp = Map<const Matrix<double, 5, 5, RowMajor>>(curv2localjacp.Array()); 
        

        //FIXME take care of this elsewhere for the moment
        const bool genconstraint = dogen && ihit==0;
//         const bool genconstraint = false;
        
        if (ihit < (nhits-1)) {

          //momentum kink residual
          const Vector5d dx0 = Map<const Vector5d>(idx0.Array());
          

          const uint32_t gluedidip1 = trackerTopology->glued(hits[ihit + 1]->det()->geographicalId());
          const bool isgluedip1 = gluedidip1 != 0;
          const DetId parmdetidip1 = isgluedip1 ? DetId(gluedidip1) : hits[ihit + 1]->geographicalId();
          const unsigned int bfieldidx = detidparms.at(std::make_pair(6, parmdetidip1));
          fieldOffset->setOffset(corparms_[bfieldidx]);
          
//           if (ihit==0) {
//             FreeTrajectoryState tmpfts(updtsospos, updtsos.globalParameters().momentum(), updtsos.charge(), field);
//             propresult = fPropagator->geometricalPropagator().propagateWithPath(tmpfts, *hits[ihit+1]->surface());
//           }
//           else {
//             propresult = fPropagator->geometricalPropagator().propagateWithPath(updtsos, *hits[ihit+1]->surface());
//           }
          

//           auto const &surfaceip1 = *hits[ihit+1]->surface();
          auto const &surfaceip1 = *surfacemap_.at(hits[ihit+1]->geographicalId());
          propresult = fPropagator->geometricalPropagator().propagateWithPath(updtsos, surfaceip1);
//           propresult = fPropagator->geometricalPropagator().propagateWithPath(updtsos, *hits[ihit+1]->surface());
          if (!propresult.first.isValid()) {
            std::cout << "Abort: Propagation Failed!" << std::endl;
            valid = false;
            break;
          }
          
//           auto const propresultalt = thePropagator->propagateWithPath(updtsos, surfaceip1);
//           
//           if (propresultalt.first.isValid()) {
//             std::cout << "nominal" << std::endl;
// //             std::cout << propresult.first.localPosition() << std::endl;
//             std::cout << propresult.first.globalMomentum().mag() << std::endl;
//             std::cout << "g4e" << std::endl;
//   //           std::cout << propresultalt.first.isValid() << std::endl;
//             std::cout << propresultalt.first.globalMomentum().mag() << std::endl;
//           }
          
          
          if (true) {
//           if (false) {
            //forward propagation jacobian (local to local)
            const Matrix<double, 5, 6> FdFp = localTransportJacobian(updtsos, propresult, false);
            
//             const Matrix<double, 5, 6> FdFpAlt = localTransportJacobianAlt(updtsos, propresult, false);
//             
//             std::cout << "FdFp" << std::endl;
//             std::cout << FdFp << std::endl;
//             std::cout << "FdFpAlt" << std::endl;
//             std::cout << FdFpAlt << std::endl;
            
//             const float finitedB = 1e-2;
//             fieldOffset->setOffset(finitedB);
//             
//             auto const propresultalt = fPropagator->geometricalPropagator().propagateWithPath(updtsos, surfaceip1);
//             
//             fieldOffset->setOffset(0.);
// 
// //             JacobianCurvilinearToCartesian curv2cartdebug(propresult.first.parameters());
// //             const AlgebraicMatrix65& jac = curv2cart.jacobian();
//             
//             JacobianCurvilinearToLocal curv2localdebug(propresult.first.surface(), propresult.first.localParameters(), *propresult.first.magneticField());
//             const AlgebraicMatrix55& curv2localjacdebug = curv2localdebug.jacobian();
//             const Matrix<double, 5, 5> Hdebug = Map<const Matrix<double, 5, 5, RowMajor>>(curv2localjacdebug.Array()); 
//             
//             const AlgebraicVector5 dparms = (propresultalt.first.localParameters().vector() - propresult.first.localParameters().vector())/finitedB;
//             
//             const Vector5d jacnom = Hdebug*FdFp.block<5,1>(0, 5);
//             
//             std::cout << "qop0 = " << updtsos.signedInverseMomentum() << std::endl;
//             std::cout << "nominal b field jacobian:" << std::endl;
//             std::cout << jacnom << std::endl;
//             std::cout << "finite diff b field jacobian:" << std::endl;
//             std::cout << dparms << std::endl;
            

            Matrix<double, 2, 2> J = FdFp.block<2, 2>(3, 3);
            // (du/dalphap)^-1
            Matrix<double, 2, 2> Sinv = FdFp.block<2, 2>(3, 1).inverse();
            // du/dqopp
            Matrix<double, 2, 1> D = FdFp.block<2, 1>(3, 0);
            // du/dBp
            Matrix<double, 2, 1> Bstate = FdFp.block<2, 1>(3, 5);
            
            const unsigned int jacstateidxout = 5*(ihit+1);
            const unsigned int jacstateidxin = 3*(ihit+1);
            
            // qop_i
            statejac(jacstateidxout, jacstateidxin + 2) = 1.;
            // dalpha_i/dqop_i
            statejac.block<2, 1>(jacstateidxout + 1, jacstateidxin + 2) = -Sinv*D;
            // dalpha_i/du_i
            statejac.block<2, 2>(jacstateidxout + 1, jacstateidxin) = -Sinv*J;
            // dalpha_i/du_(i+1)
            statejac.block<2, 2>(jacstateidxout + 1, jacstateidxin + 3) = Sinv;
            // d(lambda, phi) / dbeta
            statejac.block<2, 1>(jacstateidxout + 1, nstateparms + parmidx) = -Sinv*Bstate;
            // xlocal_i
            statejac(jacstateidxout + 3, jacstateidxin) = 1.;
            // ylocal_i
            statejac(jacstateidxout + 4, jacstateidxin + 1) = 1.;

//             std::cout << "FdFm" << std::endl;
//             std::cout << FdFm << std::endl;
//             std::cout << "FdFp" << std::endl;
//             std::cout << FdFp << std::endl;
            
            constexpr unsigned int nlocalstate = 8;
            constexpr unsigned int nlocalbfield = 3;
            constexpr unsigned int nlocaleloss = 2;
            constexpr unsigned int nlocalparms = nlocalbfield + nlocaleloss;
            
            constexpr unsigned int nlocal = nlocalstate + nlocalbfield + nlocaleloss;
            
            constexpr unsigned int localstateidx = 0;
  //           constexpr unsigned int localbfieldidx = localstateidx + nlocalstate;
  //           constexpr unsigned int localelossidx = localbfieldidx + nlocalbfield;
            constexpr unsigned int localparmidx = localstateidx + nlocalstate;
            
            const unsigned int fullstateidx = 3*ihit;
//             const unsigned int fullstateidx = 3*(ihit-1);
            const unsigned int fullparmidx = (nstateparms + parmidx) - 2;
            
//             std::cout << "ihit = " << ihit << " nstateparms = " << nstateparms << " parmidx = " << parmidx << " fullparmidx = " << fullparmidx << std::endl;
             
            // individual pieces, now starting to cast to active scalars for autograd,
            // as in eq (3) of https://doi.org/10.1016/j.cpc.2011.03.017
            // du/dum
            Matrix<MSScalar, 2, 2> Jm = FdFm.block<2, 2>(3, 3).cast<MSScalar>();
            // (du/dalpham)^-1
            Matrix<MSScalar, 2, 2> Sinvm = FdFm.block<2, 2>(3, 1).inverse().cast<MSScalar>();
            // du/dqopm
            Matrix<MSScalar, 2, 1> Dm = FdFm.block<2, 1>(3, 0).cast<MSScalar>();
            // du/dBm
            Matrix<MSScalar, 2, 1> Bm = FdFm.block<2, 1>(3, 5).cast<MSScalar>();

            // du/dup
            Matrix<MSScalar, 2, 2> Jp = FdFp.block<2, 2>(3, 3).cast<MSScalar>();
            // (du/dalphap)^-1
            Matrix<MSScalar, 2, 2> Sinvp = FdFp.block<2, 2>(3, 1).inverse().cast<MSScalar>();
            // du/dqopp
            Matrix<MSScalar, 2, 1> Dp = FdFp.block<2, 1>(3, 0).cast<MSScalar>();
            // du/dBp
            Matrix<MSScalar, 2, 1> Bp = FdFp.block<2, 1>(3, 5).cast<MSScalar>();
            
//             std::cout << "Jm" << std::endl;
//             std::cout << Jm << std::endl;
//             std::cout << "Sinvm" << std::endl;
//             std::cout << Sinvm << std::endl;
//             std::cout << "Dm" << std::endl;
//             std::cout << Dm << std::endl;
//             std::cout << "Bm" << std::endl;
//             std::cout << Bm << std::endl;
//             
//             std::cout << "Jp" << std::endl;
//             std::cout << Jp << std::endl;
//             std::cout << "Sinvp" << std::endl;
//             std::cout << Sinvp << std::endl;
//             std::cout << "Dp" << std::endl;
//             std::cout << Dp << std::endl;
//             std::cout << "Bp" << std::endl;
//             std::cout << Bp << std::endl;
            
            // energy loss jacobians
  //           const MSJacobian E = EdE.leftCols<5>().cast<MSScalar>();
  //           const MSVector dE = EdE.rightCols<1>().cast<MSScalar>();
            
            // fraction of material on this layer compared to glued layer if relevant
//             double xifraction = isglued ? preciseHit->det()->surface().mediumProperties().xi()/parmDet->surface().mediumProperties().xi() : 1.;
            
//             std::cout << "xifraction: " << xifraction << std::endl;
            
            const MSScalar Eqop(EdE(0,0));
            const Matrix<MSScalar, 1, 2> Ealpha = EdE.block<1, 2>(0, 1).cast<MSScalar>();
//             const MSScalar dE(EdE(0,5));
            const MSScalar dE(xifraction*EdE(0,5));
//             (void)EdE;
            
            const MSScalar muE(dxeloss[0]);
            
//             std::cout<<"EdE" << std::endl;
//             std::cout << EdE << std::endl;
            
            //energy loss inverse variance
            MSScalar invSigmaE(1./Q(0,0));
            
            // multiple scattering inverse covariance
            Matrix<MSScalar, 2, 2> Qinvms = Q.block<2,2>(1,1).inverse().cast<MSScalar>();
                        
            // initialize active scalars for state parameters
            Matrix<MSScalar, 2, 1> dum = Matrix<MSScalar, 2, 1>::Zero();
            //suppress gradients of reference point parameters when fitting with gen constraint
            for (unsigned int j=0; j<dum.size(); ++j) {
              init_twice_active_var(dum[j], nlocal, localstateidx + j);
              //FIXME this would be the correct condition if we were using it
//               if (dogen && ihit < 2) {
//               if (genconstraint) {
//                 init_twice_active_null(dum[j], nlocal);
//               }
//               else {
//                 init_twice_active_var(dum[j], nlocal, localstateidx + j);
//               }
            }

            
            MSScalar dqopm(0.);
            init_twice_active_var(dqopm, nlocal, localstateidx + 2);
            
//             //suppress gradients of reference point parameters when fitting with gen constraint
//             if (genconstraint) {
//               init_twice_active_null(dqopm, nlocal);
//             }
//             else {
//               init_twice_active_var(dqopm, nlocal, localstateidx + 2);
//             }

            Matrix<MSScalar, 2, 1> du = Matrix<MSScalar, 2, 1>::Zero();
            for (unsigned int j=0; j<du.size(); ++j) {
              init_twice_active_var(du[j], nlocal, localstateidx + 3 + j);
//               if (genconstraint) {
//                 init_twice_active_null(du[j], nlocal);
//               }
//               else {
//                 init_twice_active_var(du[j], nlocal, localstateidx + 3 + j);
//               }
            }
            
            MSScalar dqop(0.);
            init_twice_active_var(dqop, nlocal, localstateidx + 5);

            Matrix<MSScalar, 2, 1> dup = Matrix<MSScalar, 2, 1>::Zero();
            for (unsigned int j=0; j<dup.size(); ++j) {
              init_twice_active_var(dup[j], nlocal, localstateidx + 6 + j);
            }
  
            // initialize active scalars for correction parameters
            
            // only used for gen constraint
            MSScalar dbetam(0.);
            
            MSScalar dbeta(0.);
            init_twice_active_var(dbeta, nlocal, localparmidx + 2);
//             if (!isback) {
//               init_twice_active_var(dbeta, nlocal, localparmidx + 2);
//             }
            
            MSScalar dxi(0.);
            init_twice_active_var(dxi, nlocal, localparmidx + 3);
            
            MSScalar dbetap(0.);
            init_twice_active_var(dbetap, nlocal, localparmidx + 4);
//             if (!isfront) {
//               init_twice_active_var(dbetap, nlocal, localparmidx + 4);
//             }
            
            if (dogen && ihit==0) {
              du = Bpref.cast<MSScalar>()*dbeta;
            }
            else if (dogen && ihit==1) {
              init_twice_active_var(dbetam, nlocal, localparmidx);
              dum = Bpref.cast<MSScalar>()*dbetam;
            }
            
            //multiple scattering kink term
            
//             Matrix<MSScalar, 2, 2> Halphalamphim = Hm.block<2,2>(1, 1).cast<MSScalar>();
//             Matrix<MSScalar, 2, 2> Halphaum = Hm.block<2,2>(1, 3).cast<MSScalar>();
            
//             Matrix<MSScalar, 2, 2> Halphalamphip = Hp.block<2,2>(1, 1).cast<MSScalar>();
//             Matrix<MSScalar, 2, 2> Halphaup = Hp.block<2,2>(1, 3).cast<MSScalar>();
            
            const Matrix<MSScalar, 2, 1> dalpha0 = dx0.segment<2>(1).cast<MSScalar>();
   
            const Matrix<MSScalar, 2, 1> dlamphim = Sinvm*(dum - Jm*du - Dm*dqopm - Bm*dbeta);
            const Matrix<MSScalar, 2, 1> dlamphip = Sinvp*(dup - Jp*du - Dp*dqop - Bp*dbetap);
            
            const Matrix<MSScalar, 2, 1> dalpham = dlamphim;
            const Matrix<MSScalar, 2, 1> dalphap = dlamphip;
            
            
//             const Matrix<MSScalar, 2, 1> dalpham = Sinvm*(dum - Jm*du - Dm*dqopm - Bm*dbeta);
//             const Matrix<MSScalar, 2, 1> dalphap = Sinvp*(dup - Jp*du - Dp*dqop - Bp*dbetap);
//             const Matrix<MSScalar, 2, 1> dalpham = Sinvm*(dum - Jm*du - Dm*dqopm);
//             const Matrix<MSScalar, 2, 1> dalphap = Sinvp*(dup - Jp*du - Dp*dqop);
            
            
            const MSScalar deloss0(dx0[0]);

            
            
            
//             const Matrix<MSScalar, 2, 1> dms = dalpha0 + dalphap - dalpham;
//             const MSScalar chisqms = dms.transpose()*Qinvms*dms;
//             //energy loss term
//             
//             
//             const MSScalar deloss = deloss0 + dqop - Eqop*dqopm - (Ealpha*dalpham)[0] - dE*dxi;
//             const MSScalar chisqeloss = deloss*deloss*invSigmaE;
//             
//             const MSScalar chisq = chisqms + chisqeloss;
            
            

            
//             const bool dolikelihood = false;
//           
            MSScalar chisq;
            
            if (!islikelihood) {
              //standard chisquared contribution
              
              const Matrix<MSScalar, 2, 1> dms = dalpha0 + dalphap - dalpham;
              const MSScalar chisqms = dms.transpose()*Qinvms*dms;
              //energy loss term
              
              
              const MSScalar deloss = deloss0 + dqop - Eqop*dqopm - (Ealpha*dalpham)[0] - dE*dxi;
//               const MSScalar deloss = deloss0 + dqop - dqopm - dE*dxi;
              const MSScalar chisqeloss = deloss*invSigmaE*deloss;
              
              chisq = chisqms + chisqeloss;
              
              chisq0val += chisq.value().value();
            }
            else {
//               islikelihood = true;
              //maximum likelihood contribution 
              const MSCovariance dQdqop = dQs[0].cast<MSScalar>();
//               const MSCovariance dQddxdz = dQs[1].cast<MSScalar>();
//               const MSCovariance dQddydz = dQs[2].cast<MSScalar>();
//               const MSCovariance dQdxi = dQs[3].cast<MSScalar>();
              
//               const MSCovariance dQ = dqopm*dQdqop + dalpham[0]*dQddxdz + dalpham[1]*dQddydz + dxi*dQdxi;
//               const MSCovariance dQ = 0.5*(dqop+dqopm)*dQdqop;
              const MSCovariance dQ = dqopm*dQdqop;
//               const MSCovariance dQ = 0.5*(dqop+dqopm)*dQdqop + 0.5*(dalpham[0] + dalphap[0])*dQddxdz + 0.5*(dalpham[1]+dalphap[1])*dQddydz + dxi*dQdxi;
              
              const Matrix<MSScalar, 2, 2> Qmsnom = Q.block<2,2>(1,1).cast<MSScalar>();
              const Matrix<MSScalar, 2, 2> Qmsnominv = Qmsnom.inverse();
              const Matrix<MSScalar, 2, 2> Qmsinv = Qmsnominv - Qmsnominv*dQ.block<2,2>(1,1)*Qmsnominv;
              
              
//               const Matrix<MSScalar, 2, 2> Qms = Q.block<2,2>(1,1).cast<MSScalar>() + dQ.block<2,2>(1,1);
//               const Matrix<MSScalar, 2, 2> Qmsinv = Qms.inverse();
//               const MSScalar logdetQms = Eigen::log(Qms.determinant());
              
              const Matrix<MSScalar, 2, 1> dms = dalpha0 + dalphap - dalpham;
              MSScalar chisqms = dms.transpose()*Qmsinv*dms;
//               chisqms = chisqms + logdetQms;
              
              //energy loss term
//               const MSScalar sigmaE = MSScalar(Q(0,0)) + dQ(0,0);
//               const MSScalar sigmaEinv = 1./sigmaE;
              
              const MSScalar sigmaEnom = MSScalar(Q(0,0));
              const MSScalar sigmaEnominv = 1./sigmaEnom;
              
              const MSScalar sigmaEinv = sigmaEnominv - sigmaEnominv*dQ(0,0)*sigmaEnominv;
              
//               const MSScalar logsigmaE = Eigen::log(sigmaE);
              
              const MSScalar deloss = deloss0 + dqop - Eqop*dqopm - (Ealpha*dalpham)[0] - dE*dxi;
              MSScalar chisqeloss = deloss*sigmaEinv*deloss;
//               chisqeloss = chisqeloss + logsigmaE;
              
              chisq = chisqms + chisqeloss;
              
              //compute contributions to hessian matrix-vector derivative
              for (unsigned int i=0; i<nlocal; ++i) {
                MSScalar x(0.);
                init_twice_active_var(x, nlocal, i);
                
                Matrix<MSScalar, 2, 2> dQmsinv;
                for (unsigned int j=0; j<2; ++j) {
                  for (unsigned int k=0; k<2; ++k) {
                    dQmsinv(j,k) = MSScalar(Qmsinv(j,k).value().derivatives()[i]);
                  }
                }
                const MSScalar dSigmaEinv(sigmaEinv.value().derivatives()[i]);
                
                MSScalar dchisqms = dms.transpose()*dQmsinv*x*dms;
//                 dchisqms = 3.*dchisqms;
                MSScalar dchisqeloss = deloss*deloss*dSigmaEinv*x;
//                 dchisqeloss = 3.*dchisqeloss;
                const MSScalar dchisq = dchisqms + dchisqeloss;
                
                //TODO should this be 11x11 instead?
                //TODO check additional factor of 2
                for (unsigned int j=0; j<8; ++j) {
                  for (unsigned int k=0; k<8; ++k) {
                    dhessv[ihit][i](j,k) = dchisq.derivatives()[j].derivatives()[k];
                  }
                }
                
              }
              
            }
            
          
            
//             const MSScalar chisq = chisqms;

  //           std::cout << "chisq.value()" << std::endl;
  //           std::cout << chisq.value() << std::endl;
  //           std::cout << "chisq.value().derivatives()" << std::endl;
  //           std::cout << chisq.value().derivatives() << std::endl;
  //           std::cout << "chisq.derivatives()[0].derivatives()" << std::endl;
  //           std::cout << chisq.derivatives()[0].derivatives() << std::endl;
            
            
            //           const MSVector dms = dx0 + H*dx - E*Hprop*F*dxprev - E*Hprop*dF*dbeta - dE*dxi;
            
            
            
            
  //           MSScalar chisq;
  //           
  //           if (ihit==0 || ihit == (nhits-1)) {
  //             //standard fit
  //             const MSVector dms = dx0 + H*dx - E*Hprop*F*dxprev - E*Hprop*dF*dbeta - dE*dxi;
  //             chisq = dms.transpose()*Qinv*dms;            
  //           }
  //           else {
  //             //maximum likelihood fit
  //             const MSVector dxprop = Hprop*F*dxprev;
  //             const MSCovariance dQdxprop0 = dQs[0].cast<MSScalar>();
  //             const MSCovariance dQdxprop1 = dQs[1].cast<MSScalar>();
  //             const MSCovariance dQdxprop2 = dQs[2].cast<MSScalar>();
  //             const MSCovariance dQdxi = dQs[3].cast<MSScalar>();
  //             
  //             const MSCovariance dQ = dxprop[0]*dQdxprop0 + dxprop[1]*dQdxprop1 + dxprop[2]*dQdxprop2 + dxi*dQdxi;
  //             
  // //             const MSCovariance dQdxprop0 = dQs[0].cast<MSScalar>();
  // //             const MSCovariance d2Qdxprop02 = dQs[1].cast<MSScalar>();
  // //   //         
  // //             const MSCovariance dQ = dxprop[0]*dQdxprop0 + 0.5*dxprop[0]*dxprop[0]*d2Qdxprop02;
  //             
  // //             const Matrix<MSScalar, 3, 3> Qms = iQ.topLeftCorner<3,3>().cast<MSScalar>() + dQ.topLeftCorner<3,3>();
  // //             Qinv.topLeftCorner<3,3>() = Qms.inverse();
  //             const Matrix<MSScalar, 2, 2> Qms = iQ.block<2,2>(1,1).cast<MSScalar>() + dQ.block<2,2>(1,1);
  //             Qinv.block<2,2>(1,1) = Qms.inverse();
  //             
  //             const MSScalar logdetQ = Eigen::log(Qms.determinant());
  // 
  //             const MSVector dms = dx0 + H*dx - E*dxprop - E*Hprop*dF*dbeta - dE*dxi;
  //             chisq = dms.transpose()*Qinv*dms;
  //             chisq = chisq + logdetQ;            
  //             
  //           }
            
    //         MSCovariance Q = iQ.cast<MSScalar>();
            
    //         const MSVector dxprop = Hprop*F*dxprev;
    //         const MSCovariance dQdxprop0 = dQs[0].cast<MSScalar>();
    //         const MSCovariance dQdxprop1 = dQs[1].cast<MSScalar>();
    //         const MSCovariance dQdxprop2 = dQs[2].cast<MSScalar>();
    //         const MSCovariance dQdxi = dQs[3].cast<MSScalar>();
    // //         
    //         const MSCovariance dQ = dxprop[0]*dQdxprop0 + dxprop[1]*dQdxprop1 + dxprop[2]*dQdxprop2 + dxi*dQdxi;
      
    //         const MSVector dxprop = Hprop*F*dxprev;
    //         const MSCovariance dQdxprop0 = dQs[0].cast<MSScalar>();
    //         const MSCovariance d2Qdxprop02 = dQs[1].cast<MSScalar>();
    //         
    //         const MSCovariance dQ = dxprop[0]*dQdxprop0 + 0.5*dxprop[0]*dxprop[0]*d2Qdxprop02;
            
    //         MSCovariance Qinv = MSCovariance::Zero();
    //         Qinv(3,3) = MSScalar(1./epsxy/epsxy);
    //         Qinv(4,4) = MSScalar(1./epsxy/epsxy);
    //         Qinv.block<2,2>(1,1) = iQ.block<2,2>(1,1).inverse().cast<MSScalar>();
    //         const MSScalar Qelos = MSScalar(iQ(0,0)) + dQ(0,0);
    //         Qinv(0,0) = 1./Qelos;
    // //         const Matrix<MSScalar, 3, 3> Qms = iQ.topLeftCorner<3,3>().cast<MSScalar>() + dQ.topLeftCorner<3,3>();
    // //         Qinv.topLeftCorner<3,3>() = Qms.inverse();
    // //         const MSScalar logdetQ = Eigen::log(Qms.determinant());
    //         const MSScalar logdetQ = Eigen::log(Qelos);
    // //         
    //         const MSVector dms = dx0 + H*dx - E*dxprop - E*Hprop*dF*dbeta - dE*dxi;
    //         MSScalar chisq = dms.transpose()*Qinv*dms;
    //         chisq = chisq + logdetQ;
            
    //         const MSCovariance Qinvmod = Qinv - Qinv*dQ*Qinv;
    //         const MSScalar dlogdetQ = Eigen::log(1. + (Qinv*dQ).trace());
    //         
    //         const MSVector dms = dx0 + H*dx - E*dxprop - E*Hprop*dF*dbeta - dE*dxi;
    //         MSScalar chisq = dms.transpose()*Qinvmod*dms;
    //         chisq = chisq + dlogdetQ;

            
            auto const& gradlocal = chisq.value().derivatives();
            //fill local hessian
            Matrix<double, nlocal, nlocal> hesslocal;
            for (unsigned int j=0; j<nlocal; ++j) {
              hesslocal.row(j) = chisq.derivatives()[j].derivatives();
            }
            
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
            
            const unsigned int bfieldglobalidx = detidparms.at(std::make_pair(6, parmdetid));
            globalidxv[parmidx] = bfieldglobalidx;
            parmidx++;
            
            const unsigned int elossglobalidx = detidparms.at(std::make_pair(7, parmdetid));
            globalidxv[parmidx] = elossglobalidx;
            parmidx++;
          }
                    
          //backwards propagation jacobian (local to local) to be used at the next layer
          FdFm = localTransportJacobian(updtsos, propresult, true);
          
        }
        else {
          
          // special case for last hit, assume zero scattering angle and nominal energy loss on last layer
          Matrix<double, 2, 2> J = FdFm.block<2, 2>(3, 3);
          // (du/dalphap)^-1
          Matrix<double, 2, 2> Sinv = FdFm.block<2, 2>(3, 1).inverse();
          // du/dqopp
          Matrix<double, 2, 1> D = FdFm.block<2, 1>(3, 0);
          // du/dBp
          Matrix<double, 2, 1> Bstate = FdFm.block<2, 1>(3, 5);
          
          const unsigned int jacstateidxout = 5*(ihit+1);
          const unsigned int jacstateidxin = 3*(ihit+1);
          
          // qop_i
          //FIXME this should account for the energy loss, but only needed for outer momentum estimate which is not currently used
          statejac(jacstateidxout, jacstateidxin - 1) = 1.;
          
          // dalpha_i/dqop_i
          statejac.block<2, 1>(jacstateidxout + 1, jacstateidxin - 1) = -Sinv*D;
          // dalpha_i/du_i
          statejac.block<2, 2>(jacstateidxout + 1, jacstateidxin) = -Sinv*J;
          // dalpha_i/du_(i-1)
          statejac.block<2, 2>(jacstateidxout + 1, jacstateidxin - 3) = Sinv;
          // d(lambda, phi) / dbeta
          statejac.block<2, 1>(jacstateidxout + 1, nstateparms + parmidx) = -Sinv*Bstate;
          // xlocal_i
          statejac(jacstateidxout + 3, jacstateidxin) = 1.;
          // ylocal_i
          statejac(jacstateidxout + 4, jacstateidxin + 1) = 1.;

          
          const unsigned int bfieldglobalidx = detidparms.at(std::make_pair(6, parmdetid));
          globalidxv[parmidx] = bfieldglobalidx;
          parmidx++; 
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
            
            //FIXME dirty hack to abuse state idx for reference point magnetic field
            const unsigned int fullstateidx = genconstraint ? nstateparms : 3*(ihit+1);
  //           const unsigned int fullstateidx = 3*ihit;
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
            
            
//             Matrix<AlignScalar, 2, 2> Hu = Hp.bottomRightCorner<2,2>().cast<AlignScalar>();

            Matrix<AlignScalar, 2, 1> dy0;
            Matrix<AlignScalar, 2, 2> Vinv;
            // rotation from module to strip coordinates
//             Matrix<AlignScalar, 2, 2> R;
            Matrix2d R;
            if (preciseHit->dimension() == 1) {
//               dy0[0] = AlignScalar(matchedsim->localPosition().x() - updtsos.localPosition().x());
              dy0[0] = AlignScalar(preciseHit->localPosition().x() - updtsos.localPosition().x());
              dy0[1] = AlignScalar(0.);
              
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
              
              Vinv = Matrix<AlignScalar, 2, 2>::Zero();
              Vinv(0,0) = 1./preciseHit->localPositionError().xx();
              
//               R = Matrix<AlignScalar, 2, 2>::Identity();
              R = Matrix2d::Identity();
            }
            else {
              // 2d hit
              Matrix2d iV;
              iV << preciseHit->localPositionError().xx(), preciseHit->localPositionError().xy(),
                    preciseHit->localPositionError().xy(), preciseHit->localPositionError().yy();
              if (ispixel) {
                //take 2d hit as-is for pixels
//                 dy0[0] = AlignScalar(matchedsim->localPosition().x() - updtsos.localPosition().x());
//                 dy0[1] = AlignScalar(matchedsim->localPosition().y() - updtsos.localPosition().y());
                dy0[0] = AlignScalar(preciseHit->localPosition().x() - updtsos.localPosition().x());
                dy0[1] = AlignScalar(preciseHit->localPosition().y() - updtsos.localPosition().y());
              
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
                // diagonalize and take only smallest eigenvalue for 2d hits in strip wedge modules,
                // since the constraint parallel to the strip is spurious
                SelfAdjointEigenSolver<Matrix2d> eigensolver(iV);
//                 const Matrix2d& v = eigensolver.eigenvectors();
                R = eigensolver.eigenvectors().transpose();
                if (R(0,0) < 0.) {
                  R.row(0) *= -1.;
                }
                if (R(1,1) <0.) {
                  R.row(1) *= -1.;
                }
                
                Matrix<double, 2, 1> dy0local;
//                 dy0local[0] = matchedsim->localPosition().x() - updtsos.localPosition().x();
//                 dy0local[1] = matchedsim->localPosition().y() - updtsos.localPosition().y();
                dy0local[0] = preciseHit->localPosition().x() - updtsos.localPosition().x();
                dy0local[1] = preciseHit->localPosition().y() - updtsos.localPosition().y();
                
//                 bool simvalid = false;
//                 for (auto const& simhith : simHits) {
//                   for (const PSimHit& simHit : *simhith) {
//                     if (simHit.detUnitId() == preciseHit->geographicalId()) {                      
//                       
//                       dy0local[0] = simHit.localPosition().x() - updtsos.localPosition().x();
//                       dy0local[1] = simHit.localPosition().y() - updtsos.localPosition().y();
//                       
//                       simvalid = true;
//                       break;
//                     }
//                   }
//                   if (simvalid) {
//                     break;
//                   }
//                 }
                
                const Matrix<double, 2, 1> dy0eig = R*dy0local;
                
                //TODO deal properly with rotations (rotate back to module local coords?)
                dy0[0] = AlignScalar(dy0eig[0]);
                dy0[1] = AlignScalar(0.);
                
                Vinv = Matrix<AlignScalar, 2, 2>::Zero();
                Vinv(0,0) = AlignScalar(1./eigensolver.eigenvalues()[0]);      
                
//                 R = v.transpose().cast<AlignScalar>();
                
              }
            }
            
            rxfull.row(ivalidhit) = R.row(0).cast<float>();
            ryfull.row(ivalidhit) = R.row(1).cast<float>();
            
            validdxeigjac.block<2,2>(2*ivalidhit, 3*(ihit+1)) = R*Hp.bottomRightCorner<2,2>();
            
            const Matrix<AlignScalar, 2, 2> Ralign = R.cast<AlignScalar>();
            
            Matrix<AlignScalar, 2, 1> dx = Matrix<AlignScalar, 2, 1>::Zero();
            AlignScalar dbeta(0.);
            if (!genconstraint) {
              for (unsigned int j=0; j<dx.size(); ++j) {
                init_twice_active_var(dx[j], nlocal, localstateidx + j);
              }
            }
            else {
              init_twice_active_var(dbeta, nlocal, localstateidx);
              dx = Bpref.cast<AlignScalar>()*dbeta;
            }

            Matrix<AlignScalar, 6, 1> dalpha = Matrix<AlignScalar, 6, 1>::Zero();
            // order in which to use parameters, especially relevant in case nlocalalignment < 6
            constexpr std::array<unsigned int, 6> alphaidxs = {{5, 0, 1, 2, 3, 4}};
            for (unsigned int idim=0; idim<nlocalalignment; ++idim) {
//               init_twice_active_var(dalpha[idim], nlocal, localalignmentidx+idim);
              init_twice_active_var(dalpha[alphaidxs[idim]], nlocal, localalignmentidx+idim);
            }
            
            // alignment jacobian
            Matrix<AlignScalar, 2, 6> A = Matrix<AlignScalar, 2, 6>::Zero();

                        
            // dx/dx
            A(0,0) = AlignScalar(1.);
            // dy/dy
            A(1,1) = AlignScalar(1.);
            // dx/dz
            A(0,2) = updtsos.localParameters().dxdz();
            // dy/dz
            A(1,2) = updtsos.localParameters().dydz();
            // dx/dtheta_x
            A(0,3) = -updtsos.localPosition().y()*updtsos.localParameters().dxdz();
            // dy/dtheta_x
            A(1,3) = -updtsos.localPosition().y()*updtsos.localParameters().dydz();
            // dx/dtheta_y
            A(0,4) = -updtsos.localPosition().x()*updtsos.localParameters().dxdz();
            // dy/dtheta_y
            A(1,4) = -updtsos.localPosition().x()*updtsos.localParameters().dydz();
            // dx/dtheta_z
            A(0,5) = -updtsos.localPosition().y();
            // dy/dtheta_z
            A(1,5) = updtsos.localPosition().x();
            
            
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
            
                      
            double thetaincidence = std::asin(1./std::sqrt(std::pow(updtsos.localParameters().dxdz(),2) + std::pow(updtsos.localParameters().dydz(),2) + 1.));
            
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

            Matrix<AlignScalar, 2, 1> dh = dy0 - Ralign*dx - Ralign*A*dalpha;
            AlignScalar chisq = dh.transpose()*Vinv*dh;
            
            chisq0val += chisq.value().value();

            auto const& gradlocal = chisq.value().derivatives();
            //fill local hessian
            Matrix<double, nlocal, nlocal> hesslocal;
            for (unsigned int j=0; j<nlocal; ++j) {
              hesslocal.row(j) = chisq.derivatives()[j].derivatives();
            }
            
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
              const unsigned int xglobalidx = detidparms.at(std::make_pair(alphaidxs[idim], preciseHit->geographicalId()));
              globalidxv[nparsBfield + nparsEloss + alignmentparmidx] = xglobalidx;
              alignmentparmidx++;
              if (alphaidxs[idim]==0) {
                hitidxv.push_back(xglobalidx);
              }
            }
            
            if (iiter == 0) {
              
              // fill hit validation information
              Vector2d dyrecgenlocal;
              dyrecgenlocal << dy0[0].value().value(), dy0[1].value().value();
              const Vector2d dyrecgeneig = R*dyrecgenlocal;
              dxrecgen.push_back(dyrecgeneig[0]);
              dyrecgen.push_back(dyrecgeneig[1]);
              
              dxerr.push_back(1./std::sqrt(Vinv(0,0).value().value()));
              dyerr.push_back(1./std::sqrt(Vinv(1,1).value().value()));

              localqop.push_back(updtsos.localParameters().qbp());
              localdxdz.push_back(updtsos.localParameters().dxdz());
              localdydz.push_back(updtsos.localParameters().dydz());
              localx.push_back(updtsos.localPosition().x());
              localy.push_back(updtsos.localPosition().y());
              
              
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
              }
              else {
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
                bool simvalid = false;
                for (auto const& simhith : simHits) {
                  for (const PSimHit& simHit : *simhith) {
                    if (simHit.detUnitId() == preciseHit->geographicalId()) {
                      
                      if (int(simHit.trackId()) != simtrackid) {
                        continue;
                      }
                      
  //                     std::cout << "entry point: " << simHit.entryPoint() << std::endl;
  //                     std::cout << "exit point: " << simHit.exitPoint() << std::endl;
  //                     std::cout << "local position: " << simHit.localPosition() << std::endl;
  //                     std::cout << "particle type: " << simHit.particleType() << std::endl;
  //                     std::cout << "trackId: " << simHit.trackId() << std::endl;
  //                     
  // //                     if (simHit.entryPoint().z() * simHit.exitPoint().z() >=0.) {
  //                     if (std::abs(simHit.localPosition().z()) > 1e-4) {
  //                       std::cout << "anomalous simhit!" << std::endl;
  //                     }
                      
  //                     assert(simHit.entryPoint().z() * simHit.exitPoint().z() < 0.);
                      
                      Vector2d dysimgenlocal;
                      dysimgenlocal << simHit.localPosition().x() - updtsos.localPosition().x(),
                                      simHit.localPosition().y() - updtsos.localPosition().y();
                      const Vector2d dysimgeneig = R*dysimgenlocal;
                      dxsimgen.push_back(dysimgeneig[0]);
                      dysimgen.push_back(dysimgeneig[1]);
  //                     dxsimgen.push_back(simHit.localPosition().x() - updtsos.localPosition().x());
  //                     dysimgen.push_back(simHit.localPosition().y() - updtsos.localPosition().y());
                      
                      
                      Vector2d dyrecsimlocal;
                      dyrecsimlocal << preciseHit->localPosition().x() - simHit.localPosition().x(),
                                      preciseHit->localPosition().y() - simHit.localPosition().y();
                      const Vector2d dyrecsimeig = R*dyrecsimlocal;
                      dxrecsim.push_back(dyrecsimeig[0]);
                      dyrecsim.push_back(dyrecsimeig[1]);
                                      
  //                     dxrecsim.push_back(preciseHit->localPosition().x() - simHit.localPosition().x());
  //                     dyrecsim.push_back(preciseHit->localPosition().y() - simHit.localPosition().y());
                      
                      simvalid = true;
                      break;
                    }
                  }
                  if (simvalid) {
                    break;
                  }
                }
                if (!simvalid) {
                  dxsimgen.push_back(-99.);
                  dysimgen.push_back(-99.);
                  dxrecsim.push_back(-99.);
                  dyrecsim.push_back(-99.);
                }
              }
            }
            
          };
                    
          if (align2d) {
            fillAlignGrads(std::integral_constant<unsigned int, 3>());
          }
          else {
            fillAlignGrads(std::integral_constant<unsigned int, 2>());
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
      
      //fake constraint on reference point parameters
      if (dogen) {
//       if (false) {
        for (unsigned int i=0; i<5; ++i) {
          gradfull[i] = 0.;
          hessfull.row(i) *= 0.;
          hessfull.col(i) *= 0.;
          hessfull(i,i) = 1e6;
        }
//         //b field from reference point not consistently used in this case
//         gradfull[nstateparms] = 0.;
//         hessfull.row(nstateparms) *= 0.;
//         hessfull.col(nstateparms) *= 0.;
      }
      
      //now do the expensive calculations and fill outputs
      auto const& dchisqdx = gradfull.head(nstateparms);
      auto const& dchisqdparms = gradfull.tail(npars);
      
      auto const& d2chisqdx2 = hessfull.topLeftCorner(nstateparms, nstateparms);
      auto const& d2chisqdxdparms = hessfull.topRightCorner(nstateparms, npars);
      auto const& d2chisqdparms2 = hessfull.bottomRightCorner(npars, npars);
      
//       std::cout << "dchisqdx" << std::endl;
//       std::cout << dchisqdx << std::endl;
//       std::cout << "d2chisqdx2 diagonal" << std::endl;
//       std::cout << d2chisqdx2.diagonal() << std::endl;
//       std::cout << "d2chisqdx2" << std::endl;
//       std::cout << d2chisqdx2 << std::endl;
//       
//       auto const& eigenvalues = d2chisqdx2.eigenvalues();
//       std::cout << "d2chisqdx2 eigenvalues" << std::endl;
//       std::cout << eigenvalues << std::endl;
      
//       auto const& Cinvd = d2chisqdx2.ldlt();
      Cinvd.compute(d2chisqdx2);
      
      
      if (islikelihood) {
        const MatrixXd Cfull = Cinvd.solve(MatrixXd::Identity(nstateparms,nstateparms));
        
        // add ln det terms to gradient and hessian
  //       MatrixXd dhessfulli;
  //       MatrixXd dhessfullj;
  //       VectorXd dgradfull;
        // TODO should this cover correction parameter part of the matrix as well?
        for (unsigned int ihit=0; ihit<(nhits-1); ++ihit) {
          constexpr unsigned int localstateidx = 0;
          const unsigned int fullstateidx = 3*ihit;
          
          auto const& Cblock = Cfull.block<8,8>(fullstateidx, fullstateidx);
          
  //         dhessfulli = MatrixXd::Zero(nstateparms, nstateparms);
  //         dhessfullj = MatrixXd::Zero(nstateparms, nstateparms);
          
          //TODO fill correction parameter block as well
          for (unsigned int i=0; i<8; ++i) {
            gradfull[fullstateidx + i] += (Cblock*dhessv[ihit][i]).trace();
            for (unsigned int j=0; j<8; ++j) {
              hessfull(fullstateidx + i, fullstateidx + j) += (-Cblock*dhessv[ihit][j]*Cblock*dhessv[ihit][i]).trace();
            }
          }
          
        }
        
        Cinvd.compute(d2chisqdx2);
      
      }
      
      dxfull = -Cinvd.solve(dchisqdx);
      
      dxstate = statejac.leftCols(nstateparms)*dxfull;
      
      const Matrix<double, 1, 1> deltachisq = dchisqdx.transpose()*dxfull + 0.5*dxfull.transpose()*d2chisqdx2*dxfull;
      
      chisqval = chisq0val + deltachisq[0];
        
      ndof = 3*(nhits - 1) + nvalid + nvalidalign2d - nstateparms;
      
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
      
      const Vector5d dxRef = dxstate.head<5>();
      const Matrix5d Cinner = (statejac.leftCols(nstateparms)*Cinvd.solve(MatrixXd::Identity(nstateparms,nstateparms))*statejac.leftCols(nstateparms).transpose()).topLeftCorner<5,5>();
      
//       const Matrix5d Cinner = Cinvd.solve(MatrixXd::Identity(nstateparms,nstateparms)).topLeftCorner<5,5>();

      if (debugprintout_) {
        std::cout<< "dxRef" << std::endl;
        std::cout<< dxRef << std::endl;
      }
      
      //fill output with corrected state and covariance at reference point
      refParms.fill(0.);
      refCov.fill(0.);
//       const AlgebraicVector5& refVec = track.parameters();
      CurvilinearTrajectoryParameters curvparms(refFts.position(), refFts.momentum(), refFts.charge());
      const AlgebraicVector5& refVec = curvparms.vector();
      Map<Vector5f>(refParms.data()) = (Map<const Vector5d>(refVec.Array()) + dxRef).cast<float>();
      Map<Matrix<float, 5, 5, RowMajor> >(refCov.data()).triangularView<Upper>() = (2.*Cinner).cast<float>().triangularView<Upper>();
      
      if (iiter==0) {
        refParms_iter0 = refParms;
        refCov_iter0 = refCov;
      }
      
      niter = iiter + 1;
      edmval = -deltachisq[0];
      
//       std::cout << "iiter = " << iiter << " edmval = " << edmval << std::endl;
      
      if (iiter > 1 && std::abs(deltachisq[0])<1e-3) {
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
    
    if (!nValidPixelHitsFinal) {
      continue;
    }
    
    auto const& dchisqdx = gradfull.head(nstateparms);
    auto const& dchisqdparms = gradfull.tail(npars);
    
    auto const& d2chisqdx2 = hessfull.topLeftCorner(nstateparms, nstateparms);
    auto const& d2chisqdxdparms = hessfull.topRightCorner(nstateparms, npars);
    auto const& d2chisqdparms2 = hessfull.bottomRightCorner(npars, npars);
    
    dxdparms = -Cinvd.solve(d2chisqdxdparms).transpose();
    
//     if (debugprintout_) {
//       std::cout << "dxrefdparms" << std::endl;
//       std::cout << dxdparms.leftCols<5>() << std::endl;
//     }
    
    grad = dchisqdparms + dxdparms*dchisqdx;
    //TODO check the simplification
//     hess = d2chisqdparms2 + 2.*dxdparms*d2chisqdxdparms + dxdparms*d2chisqdx2*dxdparms.transpose();
    hess = d2chisqdparms2 + dxdparms*d2chisqdxdparms;
    
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
    
    gradv.clear();
    jacrefv.clear();

    gradv.resize(npars,0.);
    jacrefv.resize(5*npars, 0.);
    
    nJacRef = 5*npars;
    if (fillTrackTree_ && fillGrads_) {
      tree->SetBranchAddress("gradv", gradv.data());
    }
    if (fillTrackTree_) {
      tree->SetBranchAddress("jacrefv", jacrefv.data());
    }
    
    //eigen representation of the underlying vector storage
    Map<VectorXf> gradout(gradv.data(), npars);
    Map<Matrix<float, 5, Dynamic, RowMajor> > jacrefout(jacrefv.data(), 5, npars);
    
//     jacrefout = dxdparms.leftCols<5>().transpose().cast<float>();    
    jacrefout = ( (dxdparms*statejac.leftCols(nstateparms).transpose()).leftCols<5>().transpose() + statejac.topRightCorner(5, npars) ).cast<float>();  
    
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
    for (unsigned int i=0; i<npars; ++i) {
      const float absval = std::abs(grad[i]);
      if (absval>gradmax) {
        gradmax = absval;
      }      
    }
    
    
    hessmax = 0.;
    for (unsigned int i=0; i<npars; ++i) {
      for (unsigned int j=i; j<npars; ++j) {
        const unsigned int iidx = globalidxv[i];
        const unsigned int jidx = globalidxv[j];
        
        const float absval = std::abs(hess(i,j));
        if (absval>hessmax) {
          hessmax = absval;
        }
        
      }
      
    }
    
//     if (gradmax < 1e5 && refPt > 5.5) {
//       //fill aggregrate gradient and hessian
//       for (unsigned int i=0; i<npars; ++i) {
//         gradagg[globalidxv[i]] += grad[i];
//       }
//       
//       hessmax = 0.;
//       for (unsigned int i=0; i<npars; ++i) {
//         for (unsigned int j=i; j<npars; ++j) {
//           const unsigned int iidx = globalidxv[i];
//           const unsigned int jidx = globalidxv[j];
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
      const Matrix5d Cinner = (statejac*Cinvd.solve(MatrixXd::Identity(nstateparms,nstateparms))*statejac.transpose()).topLeftCorner<5,5>();
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
    const unsigned int nsym = npars*(1+npars)/2;
    hesspackedv.clear();    
    hesspackedv.resize(nsym, 0.);
    
    nSym = nsym;
    if (fillTrackTree_ && fillGrads_) {
      tree->SetBranchAddress("hesspackedv", hesspackedv.data());
    }
    
    Map<VectorXf> hesspacked(hesspackedv.data(), nsym);
    const Map<const VectorXu> globalidx(globalidxv.data(), npars);

    unsigned int packedidx = 0;
    for (unsigned int ipar = 0; ipar < npars; ++ipar) {
      const unsigned int segmentsize = npars - ipar;
      hesspacked.segment(packedidx, segmentsize) = hess.block<1, Dynamic>(ipar, ipar, 1, segmentsize).cast<float>();
      packedidx += segmentsize;
    }

    if (fillTrackTree_) {
      tree->Fill();
    }
  }
}

DEFINE_FWK_MODULE(ResidualGlobalCorrectionMaker);
