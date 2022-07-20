// system include files
#include <memory>

#include "ResidualGlobalCorrectionMakerBase.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
// #include "DataFormats/SiStripDetId/interface/TIDDetId.h"
// #include "DataFormats/SiStripDetId/interface/TIBDetId.h"
// #include "DataFormats/SiStripDetId/interface/TOBDetId.h"
// #include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "Geometry/CommonTopologies/interface/TkRadialStripTopology.h"
#include "DataFormats/Math/interface/approx_log.h"
#include "DataFormats/Math/interface/AlgebraicROOTObjects.h"
#include "TrackingTools/AnalyticalJacobians/interface/AnalyticalCurvilinearJacobian.h"
#include "TrackingTools/AnalyticalJacobians/interface/JacobianLocalToCurvilinear.h"
#include "TrackingTools/AnalyticalJacobians/interface/JacobianCurvilinearToLocal.h"
#include "TrackingTools/AnalyticalJacobians/interface/JacobianCurvilinearToCartesian.h"
#include "TrackingTools/AnalyticalJacobians/interface/JacobianCartesianToCurvilinear.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/ErrorFrameTransformer.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TrajectoryParametrization/interface/CurvilinearTrajectoryParameters.h"
#include "TrackingTools/KalmanUpdators/interface/KFSwitching1DUpdator.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h" 
#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryStateCombiner.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "DataFormats/TrackerRecHit2D/interface/TkCloner.h"
#include "DataFormats/TrackingRecHit/interface/KfComponentsHolder.h"
#include "DataFormats/Math/interface/invertPosDefMatrix.h"
#include "DataFormats/Math/interface/ProjectMatrix.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h" 

#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/CommonDetUnit/interface/GluedGeomDet.h"
#include "DataFormats/TrackerRecHit2D/interface/TrackerSingleRecHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderWithPropagator.h"
#include "DataFormats/TrackingRecHit/interface/InvalidTrackingRecHit.h"
#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "TrackingTools/GeomPropagators/interface/HelixBarrelPlaneCrossingByCircle.h"
#include "TrackingTools/GeomPropagators/interface/HelixArbitraryPlaneCrossing.h"
#include "CondFormats/Alignment/interface/Alignments.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "RecoLocalTracker/SiStripClusterizer/interface/SiStripClusterInfo.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "Alignment/CommonAlignment/interface/Utilities.h"




// #include "../interface/OffsetMagneticField.h"
// #include "../interface/ParmInfo.h"


#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "functions.h"
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/AutoDiff>
#include<Eigen/StdVector>
#include <iostream>
#include <functional>


//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ResidualGlobalCorrectionMakerBase::ResidualGlobalCorrectionMakerBase(const edm::ParameterSet &iConfig)

{
  //now do what ever initialization is needed
//   inputTraj_ = consumes<std::vector<Trajectory>>(edm::InputTag("TrackRefitter"));
//   inputTrack_ = consumes<TrajTrackAssociationCollection>(edm::InputTag("TrackRefitter"));
//   inputTrack_ = consumes<reco::TrackCollection>(edm::InputTag("TrackRefitter"));
//   inputIndices_ = consumes<std::vector<int> >(edm::InputTag("TrackRefitter"));
  
  
  inputTrackOrig_ = consumes<reco::TrackCollection>(edm::InputTag(iConfig.getParameter<edm::InputTag>("src")));

  
  fitFromGenParms_ = iConfig.getParameter<bool>("fitFromGenParms");
  fitFromSimParms_ = iConfig.getParameter<bool>("fitFromSimParms");
  fillTrackTree_ = iConfig.getParameter<bool>("fillTrackTree");
  fillGrads_ = iConfig.getParameter<bool>("fillGrads");
  fillJac_ = iConfig.getParameter<bool>("fillJac");
  fillRunTree_ = iConfig.getParameter<bool>("fillRunTree");
  doGen_ = iConfig.getParameter<bool>("doGen");
  doSim_ = iConfig.getParameter<bool>("doSim");
  bsConstraint_ = iConfig.getParameter<bool>("bsConstraint");
  applyHitQuality_ = iConfig.getParameter<bool>("applyHitQuality");
  doMuons_ = iConfig.getParameter<bool>("doMuons");
  doTrigger_ = iConfig.getParameter<bool>("doTrigger");
  corFile_ = iConfig.getParameter<std::string>("corFile");

  inputBs_ = consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));

  
  if (doGen_) {
//     GenParticlesToken_ = consumes<std::vector<reco::GenParticle>>(edm::InputTag("genParticles"));
//     GenParticlesToken_ = consumes<edm::View<reco::Candidate>>(edm::InputTag("genParticles"));
    GenParticlesToken_ = consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("genParticles"));
    genXyz0Token_ = consumes<math::XYZPointF>(edm::InputTag("genParticles","xyz0"));
    genEventInfoToken_ = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
    pileupSummaryToken_ = consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileupInfo"));
  }
  
  if (doSim_) {
    genParticlesBarcodeToken_ = consumes<std::vector<int>>(edm::InputTag("genParticles"));

//     inputSimHits_ = consumes<std::vector<PSimHit>>(edm::InputTag("g4SimHits","TrackerHitsTECLowTof"));
    std::vector<std::string> labels;
    labels.push_back("TrackerHitsPixelBarrelLowTof");
    labels.push_back("TrackerHitsPixelEndcapLowTof");
    labels.push_back("TrackerHitsTECLowTof");
    labels.push_back("TrackerHitsTIBLowTof");
    labels.push_back("TrackerHitsTIDLowTof");
    labels.push_back("TrackerHitsTOBLowTof");
    
    for (const std::string& label : labels) {
      inputSimHits_.push_back(consumes<std::vector<PSimHit>>(edm::InputTag("g4SimHits", label)));
    }
    
    inputSimTracks_ = consumes<std::vector<SimTrack>>(edm::InputTag("g4SimHits"));
  }
  
  if (doMuons_) {
//     inputMuons_ = consumes<reco::MuonCollection>(edm::InputTag(iConfig.getParameter<edm::InputTag>("muons")));
    inputMuons_ = consumes<edm::View<reco::Muon>>(edm::InputTag(iConfig.getParameter<edm::InputTag>("muons")));
  }
  
  if (doTrigger_) {
    inputTriggerResults_ = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT"));
    triggers_ = iConfig.getParameter<std::vector<std::string>>("triggers");
    triggerDecisions_.assign(triggers_.size(), 0);
  }
  
  inputGeometry_ = consumes<int>(edm::InputTag("geopro"));
  
  debugprintout_ = false;

  doMuonAssoc_ = iConfig.getParameter<bool>("doMuonAssoc");
  if (doMuonAssoc_) {
    inputMuonAssoc_ = consumes<edm::Association<std::vector<pat::Muon>>>(iConfig.getParameter<edm::InputTag>("src"));
  }


//   fout = new TFile("trackTreeGrads.root", "RECREATE");
//   fout = new TFile("trackTreeGradsdebug.root", "RECREATE");
//   fout = new TFile("trackTreeGrads.root", "RECREATE");
  //TODO this needs a newer root version
//   fout->SetCompressionAlgorithm(ROOT::kLZ4);
//   fout->SetCompressionLevel(3);
  
//   edm::Service<TgFileService> fs;
  
//   tree = new TTree("tree", "tree");
  

  outprefix = iConfig.getUntrackedParameter<std::string>("outprefix", "globalcor");

}

ResidualGlobalCorrectionMakerBase::~ResidualGlobalCorrectionMakerBase()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//


// ------------ method called once each job just before starting event loop  ------------
void ResidualGlobalCorrectionMakerBase::beginStream(edm::StreamID streamid)
{
  if (fillTrackTree_ || fillRunTree_) {
    std::stringstream filenamestream;
//     filenamestream << "globalcor_" << streamid.value() << ".root";
    filenamestream << outprefix << "_" << streamid.value() << ".root";
    fout = new TFile(filenamestream.str().c_str(), "RECREATE");
  }
  
//   runtree = new TTree("runtree","");
//   gradtree = fs->make<TTree>("gradtree","");
//   hesstree = fs->make<TTree>("hesstree","");
  
  
  if (fillTrackTree_) {
    tree = new TTree("tree","");
    const int basketSize = 4*1024*1024;
    tree->SetAutoFlush(0);
    
    tree->Branch("nParms", &nParms, basketSize);
//     tree->Branch("globalidxv", globalidxv.data(), "globalidxv[nParms]/i", basketSize);
    tree->Branch("globalidxv", &globalidxvfinal, basketSize);
    
    if (fillGrads_) {
      tree->Branch("gradv", gradv.data(), "gradv[nParms]/F", basketSize);
      tree->Branch("nSym", &nSym, basketSize);
      tree->Branch("hesspackedv", hesspackedv.data(), "hesspackedv[nSym]/F", basketSize);
      
      tree->Branch("gradmax", &gradmax);
      tree->Branch("hessmax", &hessmax);
    }
    
    tree->Branch("run", &run);
    tree->Branch("lumi", &lumi);
    tree->Branch("event", &event);
    
    tree->Branch("edmval", &edmval);
    tree->Branch("deltachisqval", &deltachisqval);
    tree->Branch("niter", &niter);
    
    tree->Branch("chisqval", &chisqval);
    tree->Branch("ndof", &ndof);

    tree->Branch("genweight", &genweight);
    
    tree->Branch("Pileup_nPU", &Pileup_nPU);
    tree->Branch("Pileup_nTrueInt", &Pileup_nTrueInt);

    tree->Branch("genl3d", &genl3d);
    
    nParms = 0.;

    
    for (std::size_t itrig = 0; itrig < triggers_.size(); ++itrig) {
      tree->Branch(triggers_[itrig].c_str(), &triggerDecisions_[itrig]);
    }
    
  }

}

// ------------ method called once each job just after ending the event loop  ------------
void ResidualGlobalCorrectionMakerBase::endStream()
{
  if (fout != nullptr) {

    fout->cd();

  //   TTree *gradtree = new TTree("gradtree","");
  //   unsigned int idx;
  //   double gradval;
  //   gradtree->Branch("idx",&idx);
  //   gradtree->Branch("gradval",&gradval);
  //   for (unsigned int i=0; i<gradagg.size(); ++i) {
  //     idx = i;
  //     gradval = gradagg[i];
  //     gradtree->Fill();
  //   }
  //
  //   TTree *hesstree = new TTree("hesstree","");
  //   unsigned int iidx;
  //   unsigned int jidx;
  //   double hessval;
  //   hesstree->Branch("iidx",&iidx);
  //   hesstree->Branch("jidx",&jidx);
  //   hesstree->Branch("hessval",&hessval);
  //
  //   for (auto const& item : hessaggsparse) {
  //     iidx = item.first.first;
  //     jidx = item.first.second;
  //     hessval = item.second;
  //     hesstree->Fill();
  //   }

    fout->Write();
    fout->Close();
  }
}

// ------------ method called when starting to processes a run  ------------

void 
ResidualGlobalCorrectionMakerBase::beginRun(edm::Run const& run, edm::EventSetup const& es)
{
  
  edm::ESHandle<GlobalTrackingGeometry> globalGeometry;
  es.get<GlobalTrackingGeometryRecord>().get(globalGeometry);
  
  edm::ESHandle<TrackerGeometry> globalGeometryIdeal;
  es.get<TrackerDigiGeometryRecord>().get("idealForDigi", globalGeometryIdeal);
  
  edm::ESHandle<TrackerTopology> trackerTopology;
  es.get<TrackerTopologyRcd>().get(trackerTopology);
  
//   edm::ESHandle<Propagator> thePropagator;
//   es.get<TrackingComponentsRecord>().get("RungeKuttaTrackerPropagator", thePropagator);
//   const MagneticField* field = thePropagator->magneticField();

  edm::ESHandle<MagneticField> magfield;
  es.get<IdealMagneticFieldRecord>().get(magfield);
  auto field = magfield.product();
  
  std::set<std::pair<int, DetId> > parmset;
  
  

  
  for (const GeomDet* det : globalGeometry->detUnits()) {
    if (!det) {
      continue;
    }
    if (GeomDetEnumerators::isTracker(det->subDetector())) {
      
//       std::cout << "detid: " << det->geographicalId().rawId() << std::endl;
      
//       std::cout << "detid: " << det->geographicalId().rawId() << " subdet: " << det->subDetector() << " isStereo: " << trackerTopology->isStereo(det->geographicalId()) << " isRphi: " << trackerTopology->isRPhi(det->geographicalId()) << " glued: " << trackerTopology->glued(det->geographicalId()) << " stack: " << trackerTopology->stack(det->geographicalId()) << " upper: " << trackerTopology->upper(det->geographicalId()) << " lower: " << trackerTopology->lower(det->geographicalId()) << " partner: " << trackerTopology->partnerDetId(det->geographicalId()).rawId() <<" xi: " << det->surface().mediumProperties().xi() << std::endl;

      const bool ispixel = GeomDetEnumerators::isTrackerPixel(det->subDetector());
      const bool isendcap = GeomDetEnumerators::isEndcap(det->subDetector());

      const uint32_t gluedid = trackerTopology->glued(det->geographicalId());
      const bool isglued = gluedid != 0;
      const DetId parmdetid = isglued ? DetId(gluedid) : det->geographicalId();
      
//       const uint32_t gluedid = trackerTopology->glued(det->geographicalId());
//       const bool isglued = gluedid != 0;
// //       const bool align2d = ispixel || isglued || isendcap;      
//       const DetId parmdetid = isglued ? DetId(gluedid) : det->geographicalId();
      
      const bool align2d = ispixel || isendcap;
//       const bool align2d = true;
//       const bool align2d = ispixel;
//       const bool align2d = false;
//       const bool align2d = isendcap && !ispixel;

      
      //always have parameters for local x alignment, in-plane rotation, bfield, and e-loss
      parmset.emplace(0, det->geographicalId());
//       parmset.emplace(1, det->geographicalId());
      parmset.emplace(2, det->geographicalId());
      parmset.emplace(3, det->geographicalId());
      parmset.emplace(4, det->geographicalId());
      parmset.emplace(5, det->geographicalId());
      
      if (align2d) {
        //local y alignment parameters only for pixels and wedge modules
        parmset.emplace(1, det->geographicalId());
      }
      // bfield and material parameters are associated to glued detids where applicable
      parmset.emplace(6, parmdetid);
      parmset.emplace(7, parmdetid);
    }
  }
  
//   const unsigned int netabins = 48;
//   const unsigned int nphibins = 36;
//   if (hetaphi == nullptr) {
//     hetaphi = new TH2D("hetaphi", "", netabins, -2.4, 2.4, nphibins, -M_PI, M_PI);
//   }
//   
//   for (unsigned int ibin = 0; ibin < (netabins + 2); ++ibin) {
//     for (unsigned int jbin = 1; jbin < (nphibins + 1); ++jbin) {
//       const unsigned int bin = hetaphi->GetBin(ibin, jbin);
//       parmset.emplace(8, bin);
//     }
//   }
  
  if (detidparms.empty()) {

    
    detidparms.clear();
    detidparmsrev.clear();
    
    TTree *runtree = nullptr;

    unsigned int iidx;
    int parmtype;
    unsigned int rawdetid;
    int subdet;
    int layer;
    int stereo;
    double x;
    double y;
    double z;
    double eta;
    double phi;
    double rho;
    double xi;
    double bx;
    double by;
    double bz;
    double bradial;
    double baxial;
    double b0;
    double b0trivial;
    double dx;
    double dy;
    double dz;
    double dtheta;
    double nx;
    double ny;
    double nz;
    double lxx;
    double lxy;
    double lxz;
    double lyx;
    double lyy;
    double lyz;



  //   assert(0);
    if (fillRunTree_) {
    //   TFile *runfout = new TFile("trackTreeGradsParmInfo.root", "RECREATE");
      fout->cd();
      runtree = new TTree("runtree", "");

      runtree->Branch("iidx", &iidx);
      runtree->Branch("parmtype", &parmtype);
      runtree->Branch("rawdetid", &rawdetid);
      runtree->Branch("subdet", &subdet);
      runtree->Branch("layer", &layer);
      runtree->Branch("stereo", &stereo);
      runtree->Branch("x", &x);
      runtree->Branch("y", &y);
      runtree->Branch("z", &z);
      runtree->Branch("eta", &eta);
      runtree->Branch("phi", &phi);
      runtree->Branch("rho", &rho);
      runtree->Branch("xi", &xi);
      runtree->Branch("bx", &bx);
      runtree->Branch("by", &by);
      runtree->Branch("bz", &bz);
      runtree->Branch("bradial", &bradial);
      runtree->Branch("baxial", &baxial);
      runtree->Branch("b0", &b0);
      runtree->Branch("b0trivial", &b0trivial);
      runtree->Branch("dx", &dx);
      runtree->Branch("dy", &dy);
      runtree->Branch("dz", &dz);
      runtree->Branch("dtheta", &dtheta);
      runtree->Branch("nx", &nx);
      runtree->Branch("ny", &ny);
      runtree->Branch("nz", &nz);
      runtree->Branch("lxx", &lxx);
      runtree->Branch("lxy", &lxy);
      runtree->Branch("lxz", &lxz);
      runtree->Branch("lyx", &lyx);
      runtree->Branch("lyy", &lyy);
      runtree->Branch("lyz", &lyz);

    }
    
    unsigned int globalidx = 0;
    for (const auto& key: parmset) {
//       std::cout << "parmtype = " << key.first << " detid = " << key.second.rawId() << std::endl;
      
      //fill info
      iidx = globalidx;
      parmtype = key.first;
      
      if (parmtype < 8) {
      
        const DetId& detid = key.second;
        const GeomDet* det = globalGeometry->idToDet(detid);
        
        const Surface &surface = det->surface();
        
        const Surface &surfaceIdeal = globalGeometryIdeal->idToDet(detid)->surface();
        
        const GloballyPositioned<double> surfaceD = surfaceToDouble(surface);
        const GloballyPositioned<double> surfaceDIdeal = surfaceToDouble(surfaceIdeal);
        
        const Vector3DBase<double, GlobalTag> dpos = surfaceDIdeal.position() - surfaceD.position();
        
        dx = dpos.x();
        dy = dpos.y();
        dz = dpos.z();
        
        const Vector3DBase<double, LocalTag> localx(1., 0., 0.);
        const Vector3DBase<double, GlobalTag> xglob = surfaceD.toGlobal(localx);
        const Vector3DBase<double, GlobalTag> xglobIdeal = surfaceDIdeal.toGlobal(localx);
        
        const double dtheta = std::acos(xglob.dot(xglobIdeal));
        
        
    //     if (detid.rawId() == 302122272) {
    //       std::cout << "width: " << det->surface().bounds().width() << std::endl;
    //       std::cout << "length: " << det->surface().bounds().length() << std::endl;
    //       std::cout << "thickness: " << det->surface().bounds().thickness() << std::endl;
    //     }
        
        layer = 0;
        stereo = 0;
    //     int subdet = det->subDetector();
    //     float eta = det->surface().position().eta();

        if (det->subDetector() == GeomDetEnumerators::PixelBarrel)
        {
          PXBDetId detid(det->geographicalId());
          layer = detid.layer();
        }
        else if (det->subDetector() == GeomDetEnumerators::PixelEndcap)
        {
          PXFDetId detid(det->geographicalId());
          layer = -1 * (detid.side() == 1) * detid.disk() + (detid.side() == 2) * detid.disk();
        }
        else if (det->subDetector() == GeomDetEnumerators::TIB)
        {
    //       TIBDetId detid(det->geographicalId());
    //       layer = detid.layer();
          layer = trackerTopology->tibLayer(det->geographicalId());
          stereo = trackerTopology->isStereo(det->geographicalId());
        }
        else if (det->subDetector() == GeomDetEnumerators::TOB)
        {
    //       TOBDetId detid(det->geographicalId());
    //       layer = detid.layer();
          layer = trackerTopology->tobLayer(det->geographicalId());
          stereo = trackerTopology->isStereo(det->geographicalId());
        }
        else if (det->subDetector() == GeomDetEnumerators::TID)
        {
          unsigned int side = trackerTopology->tidSide(detid);
          unsigned int wheel = trackerTopology->tidWheel(detid);
          layer = -1 * (side == 1) * wheel + (side == 2) * wheel;
          stereo = trackerTopology->isStereo(det->geographicalId());

        }
        else if (det->subDetector() == GeomDetEnumerators::TEC)
        {
    //       TECDetId detid(det->geographicalId());
    //       layer = -1 * (detid.side() == 1) * detid.wheel() + (detid.side() == 2) * detid.wheel();
          unsigned int side = trackerTopology->tecSide(detid);
          unsigned int wheel = trackerTopology->tecWheel(detid);
          layer = -1 * (side == 1) * wheel + (side == 2) * wheel;
          stereo = trackerTopology->isStereo(det->geographicalId());
        }
        
    //     ParmInfo parminfo;
    //     parminfo.parmtype = key.first;
    //     parminfo.subdet = det->subDetector();
    //     parminfo.layer = layer;
    //     parminfo.x = det->surface().position().x();
    //     parminfo.y = det->surface().position().y();
    //     parminfo.z = det->surface().position().z();
    //     parminfo.eta = det->surface().position().eta();
    //     parminfo.phi = det->surface().position().phi();
    //     parminfo.rho = det->surface().position().perp();


        rawdetid = detid;
        subdet = det->subDetector();
        //layer already set above
        x = det->surface().position().x();
        y = det->surface().position().y();
        z = det->surface().position().z();
        eta = det->surface().position().eta();
        phi = det->surface().position().phi();
        rho = det->surface().position().perp();
        xi = det->surface().mediumProperties().xi();
        
//         nx = det->surface().normalVector().x();
//         ny = det->surface().normalVector().y();
//         nz = det->surface().normalVector().z();

        const LocalVector lx(1.,0.,0.);
        const LocalVector ly(0.,1.,0.);
        const LocalVector lz(0.,0.,1.);

        const GlobalVector J1 = det->surface().toGlobal<double>(lx);
        const GlobalVector K1 = det->surface().toGlobal<double>(ly);
        const GlobalVector I1 = det->surface().toGlobal<double>(lz);

        lxx = J1.x();
        lxy = J1.y();
        lxz = J1.z();

        lyx = K1.x();
        lyy = K1.y();
        lyz = K1.z();

        nx = I1.x();
        ny = I1.y();
        nz = I1.z();
        
        auto const bfieldval = field->inTesla(det->surface().position());
        bx = bfieldval.x();
        by = bfieldval.y();
        bz = bfieldval.z();
        
        const GlobalVector posv(det->surface().position().x(), det->surface().position().y(), det->surface().position().z());
        
        const GlobalVector posxy(det->surface().position().x(), det->surface().position().y(), 0.);
        bradial = bfieldval.dot(posxy)/posxy.mag();
        
        const GlobalVector rphi = posxy.cross(GlobalVector(0.,0.,1.));
        baxial = bfieldval.dot(rphi)/rphi.mag();
        
        const GlobalVector b0dir = posv.cross(rphi)/posv.mag()/rphi.mag();
        
        b0 = bfieldval.dot(b0dir)/b0dir.mag();
        
        const GlobalVector btrivial(0., 0., 3.8);
        b0trivial = btrivial.dot(b0dir)/b0dir.mag();
        
        detidlayermap[detid] = {{ subdet, layer, stereo }};
      }
      else {
        rawdetid = -99;
        subdet = -99;
        layer = -99;
        stereo = -99;
        x = -99.;
        y = -99.;
        z = -99.;
        rho = -99.;
        xi = -99.;
        bx = -99.;
        by = -99.;
        bz = -99.;
        
        const unsigned int etaphiidx = key.second;
        int etabin;
        int phibin;
        int kbin;
        
        hetaphi->GetBinXYZ(etaphiidx, etabin, phibin, kbin);
        eta = hetaphi->GetXaxis()->GetBinCenter(etabin);
        phi = hetaphi->GetYaxis()->GetBinCenter(phibin);
        
      }

      //fill map
      detidparms.emplace(key, globalidx);
      detidparmsrev.emplace_back(key);
      
      globalidx++;
      
      if (fillRunTree_) {
        runtree->Fill();
      }
    }
    
  //   runfout->Write();
  //   runfout->Close();
    
    unsigned int nglobal = detidparms.size();
  //   std::sort(detidparms.begin(), detidparms.end());
    std::cout << "nglobalparms = " << detidparms.size() << std::endl;
    
    //initialize gradient
    if (!gradagg.size()) {
      gradagg.resize(nglobal, 0.);
    }
    
    // load corrections from previous iteration if applicable
    corparms_.assign(parmset.size(), 0.);
    
    if (!corFile_.empty()) {
      TFile *corfile = TFile::Open(corFile_.c_str());
      TTree *cortree = (TTree*)corfile->Get("parmtree");
      
      unsigned int idx;
      float val;
      
//       cortree->SetBranchAddress("idx", &idx);
      cortree->SetBranchAddress("x", &val);
      
      const unsigned int nparms = cortree->GetEntries();
      assert(nparms == parmset.size());
      for (unsigned int iparm = 0; iparm < nparms; ++iparm) {
        cortree->GetEntry(iparm);
        corparms_[iparm] = val;
//         corparms_[idx] = val;
      }
    }
  
  }
  
  surfacemap_.clear();
  surfacemapD_.clear();
  // fill map of modified surfaces with results of previous iteration if applicable
  for (const GeomDet* det : globalGeometry->detUnits()) {
//   for (const GeomDet* det : globalGeometryIdeal->detUnits()) {
    if (!det) {
      continue;
    }
    if (GeomDetEnumerators::isTracker(det->subDetector())) {
      
      const bool ispixel = GeomDetEnumerators::isTrackerPixel(det->subDetector());
      const bool isendcap = GeomDetEnumerators::isEndcap(det->subDetector());
      
      const uint32_t gluedid = trackerTopology->glued(det->geographicalId());
      const bool isglued = gluedid != 0;
      const DetId parmdetid = isglued ? DetId(gluedid) : det->geographicalId();
      const GeomDet* parmDet = isglued ? globalGeometry->idToDet(parmdetid) : det;
      
      const double xifraction = isglued ? det->surface().mediumProperties().xi()/parmDet->surface().mediumProperties().xi() : 1.;
      
//       const bool align2d = ispixel;
//       const bool align2d = false;

      const bool align2d = detidparms.count(std::make_pair(1, det->geographicalId()));

      

      
      const unsigned int dxidx = detidparms.at(std::make_pair(0, det->geographicalId()));
      const unsigned int dthetaidx = detidparms.at(std::make_pair(5, det->geographicalId()));
//       const unsigned int dbidx = detidparms.at(std::make_pair(6, parmdetid));
      const unsigned int dxiidx = detidparms.at(std::make_pair(7, parmdetid));
      
      // n.b. sign is flipped for alignment terms
      const float dx = -corparms_[dxidx];
//       const float dy = corparms_[dyidx];
      const float dtheta = -corparms_[dthetaidx];
//       const float db = corparms_[dbdx];
      const float dxi = corparms_[dxiidx];
      
      float dy = 0.;
      if (align2d) {
        const unsigned int dyidx = detidparms.at(std::make_pair(1, det->geographicalId())); 
        dy = -corparms_[dyidx];
      }
      
      const Surface &surface = det->surface();
//       Point3DBase<double, GlobalTag> pos = surface.position();
//       
//       //re-orthogonalize
//       Matrix<double, 3, 3> rot;
//       rot << surface.rotation().xx(), surface.rotation().xy(), surface.rotation().xz(),
//               surface.rotation().yx(), surface.rotation().yy(), surface.rotation().yz(),
//               surface.rotation().zx(), surface.rotation().zy(), surface.rotation().zz();
//               
// //       std::cout << "rot check pre:" << std::endl;
// //       std::cout << rot.transpose()*rot << std::endl;
//             
//       for (unsigned int i = 0; i < 3; ++i) {
//         for (unsigned int j = 0; j < i; ++j) {
//           rot.row(i) = (rot.row(i) - (rot.row(i)*rot.row(j).transpose())[0]/rot.row(j).squaredNorm()*rot.row(j)).eval();
//         }
//       }
//       for (unsigned int i = 0; i < 3; ++i) {
//           rot.row(i) = (rot.row(i)/rot.row(i).norm()).eval();
//       }
//       
// //       std::cout << "rot check post:" << std::endl;
// //       std::cout << rot.transpose()*rot << std::endl;
//       
//       const TkRotation<double> tkrot(rot(0,0), rot(0,1), rot(0,2),
//                                      rot(1,0), rot(1,1), rot(1,2),
//                                      rot(2,0), rot(2,1), rot(2,2));
//       
//       
//       GloballyPositioned<double> surfaceD(pos, tkrot);
      
      GloballyPositioned<double> surfaceD = surfaceToDouble(surface);
      
      
      ReferenceCountingPointer<Plane> plane = Plane::build(det->surface());
//       std::shared_ptr<Plane> plane = std::make_shared<Plane>(det->surface());
      
      //move, rotate and modify material
      plane->move(plane->toGlobal(LocalVector(dx, dy, 0.)));
      plane->rotate(Surface::RotationType(Surface::RotationType::BasicVector(plane->toGlobal(LocalVector(0.,0.,1.))), dtheta));
      
      surfaceD.move(surfaceD.toGlobal(Vector3DBase<double, LocalTag>(dx, dy, 0.)));
      surfaceD.rotate(TkRotation<double>(TkRotation<double>::BasicVector(surfaceD.toGlobal(Vector3DBase<double, LocalTag>(0.,0.,1.))), dtheta));
      
//       Matrix<double, 3, 3> rot2;
//       rot2 << surfaceD.rotation().xx(), surfaceD.rotation().xy(), surfaceD.rotation().xz(),
//               surfaceD.rotation().yx(), surfaceD.rotation().yy(), surfaceD.rotation().yz(),
//               surfaceD.rotation().zx(), surfaceD.rotation().zy(), surfaceD.rotation().zz();
//       
//       std::cout << "rot2 check:" << std::endl;
//       std::cout << rot2.transpose()*rot2 << std::endl;
              
//       if (plane->mediumProperties().xi() + xifraction*dxi < 0.05*plane->mediumProperties().xi()) {
//         std::cout << "xi value clipped!" << std::endl;
//       }
      
      const float radLen = plane->mediumProperties().radLen();
      const float xi = std::max(plane->mediumProperties().xi() + xifraction*dxi, 0.05*plane->mediumProperties().xi());
//       const float xi = 0.05*plane->mediumProperties().xi();
      const MediumProperties mp(radLen, xi);
//       plane->setMediumProperties(mp);
      
//       printf("in beginRun, detid = %u, parmdetid = %u, oldp = %p, newp = %p\n", det->geographicalId().rawId(), parmdetid.rawId(), &det->surface(), &(*plane));
//       std::cout << "oldxi = " << det->surface().mediumProperties().xi() << " newxi = " << plane->mediumProperties().xi() << " dxi = " << dxi << " xifraction = " << xifraction << std::endl;
      
//       surfacemap_[det->geographicalId()] = plane;
      
      surfacemapD_[det->geographicalId()] = surfaceD;
      

    }
    
  }
    
  
}


GloballyPositioned<double> ResidualGlobalCorrectionMakerBase::surfaceToDouble(const Surface &surface) const {
  Point3DBase<double, GlobalTag> pos = surface.position();
  //re-orthogonalize
  Matrix<double, 3, 3> rot;
  rot << surface.rotation().xx(), surface.rotation().xy(), surface.rotation().xz(),
          surface.rotation().yx(), surface.rotation().yy(), surface.rotation().yz(),
          surface.rotation().zx(), surface.rotation().zy(), surface.rotation().zz();
          
//       std::cout << "rot check pre:" << std::endl;
//       std::cout << rot.transpose()*rot << std::endl;
        
  for (unsigned int i = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < i; ++j) {
      rot.row(i) = (rot.row(i) - (rot.row(i)*rot.row(j).transpose())[0]/rot.row(j).squaredNorm()*rot.row(j)).eval();
    }
  }
  for (unsigned int i = 0; i < 3; ++i) {
      rot.row(i) = (rot.row(i)/rot.row(i).norm()).eval();
  }
  
//       std::cout << "rot check post:" << std::endl;
//       std::cout << rot.transpose()*rot << std::endl;
  
  const TkRotation<double> tkrot(rot(0,0), rot(0,1), rot(0,2),
                                  rot(1,0), rot(1,1), rot(1,2),
                                  rot(2,0), rot(2,1), rot(2,2));
  
  
  return GloballyPositioned<double>(pos, tkrot);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ResidualGlobalCorrectionMakerBase::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


Matrix<double, 5, 6> ResidualGlobalCorrectionMakerBase::hybrid2curvJacobianD(const Matrix<double, 7, 1> &state, const MagneticField *field, double dBz) const {
  
  const GlobalPoint pos(state[0], state[1], state[2]);  
  const GlobalVector &bfield = field->inInverseGeV(pos);
  const Matrix<double, 3, 1> Bv(bfield.x(), bfield.y(), double(bfield.z()) + 2.99792458e-3*dBz);
  
  const double q = state[6];
  
  const double qop0 = q/state.segment<3>(3).norm();
  const double lam0 = std::atan(state[5]/std::sqrt(state[3]*state[3] + state[4]*state[4]));
  const double phi0 = std::atan2(state[4], state[3]);
  
  const Matrix<double, 3, 1> W0 = state.segment<3>(3).normalized();
  const double W0x = W0[0];
  const double W0y = W0[1];
  const double W0z = W0[2];  

  const double x0 = state[0];
  const double y0 = state[1];
  const double z0 = state[2];
  
  const double B = Bv.norm();
  const Matrix<double, 3, 1> H = Bv.normalized();
  const double hx = H[0];
  const double hy = H[1];
  const double hz = H[2];
  
  const double xf0 = std::cos(lam0);
  const double xf1 = std::pow(xf0, 2);
  const double xf2 = std::sqrt(xf1);
  const double xf3 = 2*phi0;
  const double xf4 = std::sin(xf3);
  const double xf5 = std::cos(xf3);
  const double xf6 = B*qop0;
  const double xf7 = (1.0/4.0)*xf6*std::sqrt(2*std::cos(2*lam0) + 2);
  const double xf8 = std::tan(lam0);
  const double xf9 = std::sin(phi0);
  const double xf10 = std::cos(phi0);
  const double xf11 = std::sin(lam0);
  const double xf12 = xf6*(hx*xf10*xf11 + hy*xf11*xf9 - hz*xf0);
  const double xf13 = std::pow(W0x, 2) + std::pow(W0y, 2);
  const double xf14 = std::pow(xf13, -1.0/2.0);
  const double xf15 = W0x*xf9;
  const double xf16 = W0y*xf10;
  const double xf17 = xf1*(xf15 - xf16);
  const double xf18 = W0x*W0z;
  const double xf19 = W0y*W0z;
  const double xf20 = xf0*xf10*xf18 + xf0*xf19*xf9 - xf11*xf13;
  const double xf21 = xf0*xf20;
  const double dqopdqop0 = 1;
  const double dqopdlam0 = 0;
  const double dqopdphi0 = 0;
  const double dqopdx0 = 0;
  const double dqopdy0 = 0;
  const double dqopdz0 = 0;
  const double dlamdqop0 = 0;
  const double dlamdlam0 = xf2/xf0;
  const double dlamdphi0 = 0;
  const double dlamdx0 = xf7*(hx*xf4 - hy*xf5 - hy);
  const double dlamdy0 = xf7*(-hx*xf5 + hx - hy*xf4);
  const double dlamdz0 = xf2*xf6*xf8*(hx*xf9 - hy*xf10);
  const double dphidqop0 = 0;
  const double dphidlam0 = 0;
  const double dphidphi0 = 1;
  const double dphidx0 = -xf10*xf12;
  const double dphidy0 = -xf12*xf9;
  const double dphidz0 = -xf12*xf8;
  const double dxtdqop0 = 0;
  const double dxtdlam0 = 0;
  const double dxtdphi0 = 0;
  const double dxtdx0 = -xf14*(W0y + xf10*xf17);
  const double dxtdy0 = xf14*(W0x - xf17*xf9);
  const double dxtdz0 = xf0*xf11*xf14*(-xf15 + xf16);
  const double dytdqop0 = 0;
  const double dytdlam0 = 0;
  const double dytdphi0 = 0;
  const double dytdx0 = xf14*(xf10*xf21 - xf18);
  const double dytdy0 = xf14*(-xf19 + xf21*xf9);
  const double dytdz0 = xf14*(xf11*xf20 + xf13);
  Matrix<double, 5, 6> res;
  res(0,0) = dqopdqop0;
  res(0,1) = dqopdlam0;
  res(0,2) = dqopdphi0;
  res(0,3) = dqopdx0;
  res(0,4) = dqopdy0;
  res(0,5) = dqopdz0;
  res(1,0) = dlamdqop0;
  res(1,1) = dlamdlam0;
  res(1,2) = dlamdphi0;
  res(1,3) = dlamdx0;
  res(1,4) = dlamdy0;
  res(1,5) = dlamdz0;
  res(2,0) = dphidqop0;
  res(2,1) = dphidlam0;
  res(2,2) = dphidphi0;
  res(2,3) = dphidx0;
  res(2,4) = dphidy0;
  res(2,5) = dphidz0;
  res(3,0) = dxtdqop0;
  res(3,1) = dxtdlam0;
  res(3,2) = dxtdphi0;
  res(3,3) = dxtdx0;
  res(3,4) = dxtdy0;
  res(3,5) = dxtdz0;
  res(4,0) = dytdqop0;
  res(4,1) = dytdlam0;
  res(4,2) = dytdphi0;
  res(4,3) = dytdx0;
  res(4,4) = dytdy0;
  res(4,5) = dytdz0;
  
  return res;

}

Matrix<double, 7, 1> ResidualGlobalCorrectionMakerBase::pca2cart(const Matrix<double, 5, 1> &statepca, const reco::BeamSpot &bs) const {
  const double x0bs = bs.x0();
  const double y0bs = bs.y0();
  const double z0bs = bs.z0();
  const double dxdzbs = bs.dxdz();
  const double dydzbs = bs.dydz();

  const Matrix<double, 3, 1> bv(dxdzbs, dydzbs, 1.);

  const double qop0 = statepca[0];
  const double lam0 = statepca[1];
  const double phi0 = statepca[2];
  const double d0 = statepca[3];
  const double z0 = statepca[4];

  const double p0 = std::abs(1./qop0);
  const double q = std::copysign(1., qop0);

  const double px = p0*std::cos(lam0)*std::cos(phi0);
  const double py = p0*std::cos(lam0)*std::sin(phi0);
  const double pz = p0*std::sin(lam0);

  const Matrix<double, 3, 1> mom(px, py, pz);

  // pca on beamline
  Matrix<double, 3, 1> pcabs;
  pcabs[0] = x0bs + (z0 - z0bs)*dxdzbs;
  pcabs[1] = y0bs + (z0 - z0bs)*dydzbs;
  pcabs[2] = z0;

  const Matrix<double, 3, 1> dv = bv.cross(mom).normalized()*d0;

  const Matrix<double, 3, 1> xyz = pcabs + dv;

  Matrix<double, 7, 1> res;
  res.head<3>() = xyz;
  res.segment<3>(3) = mom;
  res[6] = q;

  return res;

}

Matrix<double, 5, 1> ResidualGlobalCorrectionMakerBase::cart2pca(const Matrix<double, 7, 1> &state, const reco::BeamSpot &bs) const {

  const double x0bs = bs.x0();
  const double y0bs = bs.y0();
  const double z0bs = bs.z0();
  const double dxdzbs = bs.dxdz();
  const double dydzbs = bs.dydz();

  const Matrix<double, 3, 1> xyz0bs(x0bs, y0bs, z0bs);
  const Matrix<double, 3, 1> bv(dxdzbs, dydzbs, 1.);
  const Matrix<double, 3, 1> bhat = bv.normalized();

  const Matrix<double, 3, 1> xyz0 = state.head<3>();
  const Matrix<double, 3, 1> mom0 = state.segment<3>(3);

  // track reference point in beamspot local coordinates
  const double vpll = (xyz0 - xyz0bs).dot(bhat);

  // point of closest approach on beamline
  const Matrix<double, 3, 1> pcabs = xyz0bs + vpll*bhat;
  const double z0 = pcabs[2];

  // vector from pca on beamline to track reference point
  const Matrix<double, 3, 1> dv = xyz0 - pcabs;
  const double d0 = bhat.cross(mom0).normalized().dot(dv);

  const double q = state[6];

  const double qop0 = q/state.segment<3>(3).norm();
  const double lam0 = std::atan(state[5]/std::sqrt(state[3]*state[3] + state[4]*state[4]));
  const double phi0 = std::atan2(state[4], state[3]);

  Matrix<double, 5, 1> res;
  res[0] = qop0;
  res[1] = lam0;
  res[2] = phi0;
  res[3] = d0;
  res[4] = z0;

  return res;

}

Matrix<double, 5, 5> ResidualGlobalCorrectionMakerBase::pca2curvJacobianD(const Matrix<double, 7, 1> &state, const MagneticField *field, const reco::BeamSpot &bs,  double dBz) const {

  const GlobalPoint pos(state[0], state[1], state[2]);
  const GlobalVector &bfield = field->inInverseGeV(pos);
  const Matrix<double, 3, 1> Bv(bfield.x(), bfield.y(), double(bfield.z()) + 2.99792458e-3*dBz);

  const double x0bs = bs.x0();
  const double y0bs = bs.y0();
  const double z0bs = bs.z0();
  const double dxdzbs = bs.dxdz();
  const double dydzbs = bs.dydz();

  const Matrix<double, 3, 1> W0 = state.segment<3>(3).normalized();
  const double W0x = W0[0];
  const double W0y = W0[1];
  const double W0z = W0[2];

  const Matrix<double, 5, 1> statepca = cart2pca(state, bs);
  const double qop0 = statepca[0];
  const double lam0 = statepca[1];
  const double phi0 = statepca[2];
  const double d0 = statepca[3];
  const double z0 = statepca[4];

  const double B = Bv.norm();
  const Matrix<double, 3, 1> H = Bv.normalized();
  const double hx = H[0];
  const double hy = H[1];
  const double hz = H[2];

  const double xf0 = std::sin(lam0);
  const double xf1 = std::sin(phi0);
  const double xf2 = std::pow(xf1, 2);
  const double xf3 = std::cos(lam0);
  const double xf4 = std::pow(xf3, 2);
  const double xf5 = xf2*xf4;
  const double xf6 = std::cos(phi0);
  const double xf7 = std::pow(xf6, 2);
  const double xf8 = xf4*xf7;
  const double xf9 = xf5 + xf8;
  const double xf10 = 1.0/xf9;
  const double xf11 = xf0*xf3;
  const double xf12 = std::sqrt(xf4)*(-hx*xf1 + hy*xf6);
  const double xf13 = xf12/xf3;
  const double xf14 = B*d0;
  const double xf15 = qop0*xf14;
  const double xf16 = dxdzbs*xf1;
  const double xf17 = dydzbs*xf6;
  const double xf18 = xf16 - xf17;
  const double xf19 = std::pow(dxdzbs, 2);
  const double xf20 = std::pow(dydzbs, 2);
  const double xf21 = dxdzbs*xf0;
  const double xf22 = xf3*xf6;
  const double xf23 = dydzbs*xf0;
  const double xf24 = xf1*xf3;
  const double xf25 = std::pow(-2*xf16*xf17*xf4 - xf19*xf8 + xf19 + xf20*xf4*xf7 - xf20*xf4 + xf20 - 2*xf21*xf22 - 2*xf23*xf24 + xf4, -1.0/2.0);
  const double xf26 = xf18*xf25;
  const double xf27 = xf15*xf26;
  const double xf28 = -xf1*xf23 - xf21*xf6 + xf3;
  const double xf29 = xf25*xf28;
  const double xf30 = -dxdzbs*xf22 - dydzbs*xf24 - xf0;
  const double xf31 = B*qop0;
  const double xf32 = xf30*xf31;
  const double xf33 = std::tan(lam0);
  const double xf34 = hx*xf33*xf6 + hy*xf1*xf33 - hz;
  const double xf35 = xf15*xf3;
  const double xf36 = xf21 - xf22;
  const double xf37 = xf14*xf36;
  const double xf38 = z0 - z0bs;
  const double xf39 = dydzbs*xf38 + y0bs;
  const double xf40 = xf23 - xf24;
  const double xf41 = std::pow(xf18, 2);
  const double xf42 = std::pow(xf36, 2) + xf4*xf41 + std::pow(xf40, 2);
  const double xf43 = std::sqrt(xf42);
  const double xf44 = B*xf39*xf43 - xf37;
  const double xf45 = 1.0/qop0;
  const double xf46 = std::pow(W0x, 2);
  const double xf47 = std::pow(W0y, 2);
  const double xf48 = xf46 + xf47;
  const double xf49 = std::pow(xf48, -1.0/2.0);
  const double xf50 = W0x*xf49;
  const double xf51 = 1.0/xf43;
  const double xf52 = 1.0/B;
  const double xf53 = xf51*xf52;
  const double xf54 = xf50*xf53;
  const double xf55 = xf45*xf54;
  const double xf56 = xf14*xf40;
  const double xf57 = dxdzbs*xf38 + x0bs;
  const double xf58 = B*xf43;
  const double xf59 = xf57*xf58;
  const double xf60 = xf56 + xf59;
  const double xf61 = W0y*xf49;
  const double xf62 = xf53*xf61;
  const double xf63 = xf45*xf62;
  const double xf64 = std::pow(qop0, -2);
  const double xf65 = qop0*xf56 + qop0*xf59;
  const double xf66 = B*qop0*xf39*xf43 - qop0*xf37;
  const double xf67 = W0x*xf1 - W0y*xf6;
  const double xf68 = xf3*xf67;
  const double xf69 = d0*xf49;
  const double xf70 = dxdzbs*xf3;
  const double xf71 = xf0*xf6;
  const double xf72 = dydzbs*xf3;
  const double xf73 = xf0*xf1;
  const double xf74 = -xf11*xf41 + (1.0/2.0)*xf36*(2*xf70 + 2*xf71) + (1.0/2.0)*xf40*(2*xf72 + 2*xf73);
  const double xf75 = B*qop0*xf39*xf51*xf74 - xf15*(xf70 + xf71);
  const double xf76 = xf31*xf51;
  const double xf77 = xf15*(xf72 + xf73) + xf57*xf74*xf76;
  const double xf78 = -xf74;
  const double xf79 = std::pow(xf42, -3.0/2.0);
  const double xf80 = xf45*xf52*xf79;
  const double xf81 = xf78*xf80;
  const double xf82 = xf50*xf66;
  const double xf83 = xf61*xf65;
  const double xf84 = dxdzbs*xf6;
  const double xf85 = dydzbs*xf1;
  const double xf86 = (1.0/2.0)*xf18*xf4*(2*xf84 + 2*xf85) - xf22*xf40 + xf24*xf36;
  const double xf87 = B*qop0*xf39*xf51*xf86 - xf15*xf24;
  const double xf88 = B*qop0*xf51*xf57*xf86 - xf15*xf22;
  const double xf89 = -xf86;
  const double xf90 = xf80*xf89;
  const double xf91 = xf36*xf51;
  const double xf92 = xf40*xf51;
  const double xf93 = xf30*xf49;
  const double xf94 = W0z*xf50;
  const double xf95 = xf45*xf53*xf94;
  const double xf96 = W0z*xf61;
  const double xf97 = xf45*xf53*xf96;
  const double xf98 = xf46*xf49 + xf47*xf49;
  const double xf99 = xf18*xf3;
  const double xf100 = xf58*z0;
  const double xf101 = qop0*xf100 + xf15*xf99;
  const double xf102 = xf51*xf98;
  const double xf103 = -W0x*W0z*xf22 - W0y*W0z*xf24 + xf0*xf48;
  const double xf104 = xf65*xf94;
  const double xf105 = xf66*xf96;
  const double dqopdqop0 = 1;
  const double dqopdlam0 = 0;
  const double dqopdphi0 = 0;
  const double dqopdd0 = 0;
  const double dqopdz0 = 0;
  const double dlamdqop0 = 0;
  const double dlamdlam0 = xf13*xf27 + (xf0*(xf11*xf2 + xf11*xf7)/std::pow(xf9, 3.0/2.0) + xf3/std::sqrt(xf9))/(std::pow(xf0, 2)*xf10 + 1);
  const double dlamdphi0 = xf12*xf15*xf29;
  const double dlamdd0 = 0;
  const double dlamdz0 = xf13*xf32;
  const double dphidqop0 = 0;
  const double dphidlam0 = xf27*xf34;
  const double dphidphi0 = xf10*xf5 + xf10*xf8 + xf29*xf34*xf35;
  const double dphidd0 = 0;
  const double dphidz0 = xf32*xf34;
  const double dxtdqop0 = xf44*xf55 - xf54*xf64*xf66 - xf60*xf63 + xf62*xf64*xf65;
  const double dxtdlam0 = xf26*xf68*xf69 + xf55*xf75 - xf63*xf77 + xf81*xf82 - xf81*xf83;
  const double dxtdphi0 = xf29*xf4*xf67*xf69 + xf55*xf87 - xf63*xf88 + xf82*xf90 - xf83*xf90;
  const double dxtdd0 = -xf50*xf91 - xf61*xf92;
  const double dxtdz0 = -dxdzbs*xf61 + dydzbs*xf50 + xf68*xf93;
  const double dytdqop0 = W0x*W0z*xf49*xf51*xf52*xf64*xf65 + W0y*W0z*xf49*xf51*xf52*xf64*xf66 - xf101*xf102*xf52*xf64 - xf44*xf97 + xf45*xf51*xf52*xf98*(xf100 + xf14*xf99) - xf60*xf95;
  const double dytdlam0 = d0*xf103*xf18*xf25*xf49 + xf101*xf45*xf52*xf78*xf79*xf98 - xf104*xf81 - xf105*xf81 + xf45*xf51*xf52*xf98*(B*qop0*xf51*xf74*z0 - xf0*xf15*xf18) - xf75*xf97 - xf77*xf95;
  const double dytdphi0 = d0*xf103*xf25*xf28*xf3*xf49 + xf101*xf45*xf52*xf79*xf89*xf98 - xf104*xf90 - xf105*xf90 + xf45*xf51*xf52*xf98*(xf35*(xf84 + xf85) + xf76*xf86*z0) - xf87*xf97 - xf88*xf95;
  const double dytdd0 = xf102*xf99 + xf91*xf96 - xf92*xf94;
  const double dytdz0 = -dxdzbs*xf94 - dydzbs*xf96 + xf103*xf93 + xf98;
  Matrix<double, 5, 5> res;
  res(0,0) = dqopdqop0;
  res(0,1) = dqopdlam0;
  res(0,2) = dqopdphi0;
  res(0,3) = dqopdd0;
  res(0,4) = dqopdz0;
  res(1,0) = dlamdqop0;
  res(1,1) = dlamdlam0;
  res(1,2) = dlamdphi0;
  res(1,3) = dlamdd0;
  res(1,4) = dlamdz0;
  res(2,0) = dphidqop0;
  res(2,1) = dphidlam0;
  res(2,2) = dphidphi0;
  res(2,3) = dphidd0;
  res(2,4) = dphidz0;
  res(3,0) = dxtdqop0;
  res(3,1) = dxtdlam0;
  res(3,2) = dxtdphi0;
  res(3,3) = dxtdd0;
  res(3,4) = dxtdz0;
  res(4,0) = dytdqop0;
  res(4,1) = dytdlam0;
  res(4,2) = dytdphi0;
  res(4,3) = dytdd0;
  res(4,4) = dytdz0;

  return res;

}


Matrix<double, 6, 5> ResidualGlobalCorrectionMakerBase::pca2cartJacobianD(const Matrix<double, 7, 1> &state, const reco::BeamSpot &bs) const {

  const double x0bs = bs.x0();
  const double y0bs = bs.y0();
  const double z0bs = bs.z0();
  const double dxdzbs = bs.dxdz();
  const double dydzbs = bs.dydz();

  const Matrix<double, 5, 1> statepca = cart2pca(state, bs);
  const double qop0 = statepca[0];
  const double lam0 = statepca[1];
  const double phi0 = statepca[2];
  const double d0 = statepca[3];
  const double z0 = statepca[4];

  const double xf0 = std::sin(lam0);
  const double xf1 = dxdzbs*xf0;
  const double xf2 = std::cos(phi0);
  const double xf3 = dydzbs*xf2;
  const double xf4 = std::cos(lam0);
  const double xf5 = std::pow(dxdzbs, 2);
  const double xf6 = xf2*xf4;
  const double xf7 = std::sin(phi0);
  const double xf8 = dxdzbs*xf7;
  const double xf9 = xf6*xf8;
  const double xf10 = std::pow(xf2, 2);
  const double xf11 = std::pow(dydzbs, 2);
  const double xf12 = xf1 - xf6;
  const double xf13 = xf4*xf7;
  const double xf14 = dydzbs*xf0 - xf13;
  const double xf15 = -xf3 + xf8;
  const double xf16 = std::pow(xf12, 2) + std::pow(xf14, 2) + std::pow(xf15, 2)*std::pow(xf4, 2);
  const double xf17 = d0/std::pow(xf16, 3.0/2.0);
  const double xf18 = dxdzbs*xf2 + dydzbs*xf7;
  const double xf19 = xf15*xf4;
  const double xf20 = xf12*xf7 - xf14*xf2 + xf18*xf19;
  const double xf21 = xf17*xf4;
  const double xf22 = std::pow(xf16, -1.0/2.0);
  const double xf23 = dxdzbs*xf4;
  const double xf24 = std::pow(xf7, 2);
  const double xf25 = xf0*xf11;
  const double xf26 = 1.0/qop0;
  const double xf27 = 1.0/std::fabs(qop0);
  const double xf28 = xf27*xf6;
  const double xf29 = xf0*xf27;
  const double xf30 = xf13*xf27;
  const double dxdqop0 = 0;
  const double dxdlam0 = xf17*(std::pow(dydzbs, 3)*xf10*xf4 - dydzbs*xf10*xf4*xf5 + dydzbs*xf10*xf4 + dydzbs*xf4*xf5 + xf0*xf5*xf7 - xf1*xf3 - 2*xf11*xf9 - xf9);
  const double dxdphi0 = xf21*(-xf14*xf20 - xf16*xf2);
  const double dxdd0 = xf14*xf22;
  const double dxdz0 = dxdzbs;
  const double dydqop0 = 0;
  const double dydlam0 = xf17*(-std::pow(dxdzbs, 3)*xf24*xf4 + dxdzbs*dydzbs*xf0*xf7 + dxdzbs*xf11*xf24*xf4 + 2*dydzbs*xf2*xf4*xf5*xf7 + dydzbs*xf2*xf4*xf7 - xf11*xf23 - xf2*xf25 - xf23*xf24);
  const double dydphi0 = xf21*(xf12*xf20 - xf16*xf7);
  const double dydd0 = -xf12*xf22;
  const double dydz0 = dydzbs;
  const double dzdqop0 = 0;
  const double dzdlam0 = xf15*xf17*(dxdzbs*xf2*xf4 + dydzbs*xf4*xf7 - xf0*xf5 - xf25);
  const double dzdphi0 = xf21*(-std::pow(xf15, 2)*xf4*(dxdzbs*xf6 + dydzbs*xf13 + xf0) + xf16*xf18);
  const double dzdd0 = xf19*xf22;
  const double dzdz0 = 1;
  const double dpxdqop0 = -xf26*xf28;
  const double dpxdlam0 = -xf2*xf29;
  const double dpxdphi0 = -xf30;
  const double dpxdd0 = 0;
  const double dpxdz0 = 0;
  const double dpydqop0 = -xf26*xf30;
  const double dpydlam0 = -xf29*xf7;
  const double dpydphi0 = xf28;
  const double dpydd0 = 0;
  const double dpydz0 = 0;
  const double dpzdqop0 = -xf26*xf29;
  const double dpzdlam0 = xf27*xf4;
  const double dpzdphi0 = 0;
  const double dpzdd0 = 0;
  const double dpzdz0 = 0;
  Matrix<double, 6, 5> res;
  res(0,0) = dxdqop0;
  res(0,1) = dxdlam0;
  res(0,2) = dxdphi0;
  res(0,3) = dxdd0;
  res(0,4) = dxdz0;
  res(1,0) = dydqop0;
  res(1,1) = dydlam0;
  res(1,2) = dydphi0;
  res(1,3) = dydd0;
  res(1,4) = dydz0;
  res(2,0) = dzdqop0;
  res(2,1) = dzdlam0;
  res(2,2) = dzdphi0;
  res(2,3) = dzdd0;
  res(2,4) = dzdz0;
  res(3,0) = dpxdqop0;
  res(3,1) = dpxdlam0;
  res(3,2) = dpxdphi0;
  res(3,3) = dpxdd0;
  res(3,4) = dpxdz0;
  res(4,0) = dpydqop0;
  res(4,1) = dpydlam0;
  res(4,2) = dpydphi0;
  res(4,3) = dpydd0;
  res(4,4) = dpydz0;
  res(5,0) = dpzdqop0;
  res(5,1) = dpzdlam0;
  res(5,2) = dpzdphi0;
  res(5,3) = dpzdd0;
  res(5,4) = dpzdz0;

  return res;
}

std::array<Matrix<double, 7, 1>, 2> ResidualGlobalCorrectionMakerBase::twoTrackPca2cart(const Matrix<double, 10, 1> &statepca) const {

  const double qopa0 = statepca[0];
  const double lama0 = statepca[1];
  const double phia0 = statepca[2];
  const double qopb0 = statepca[3];
  const double lamb0 = statepca[4];
  const double phib0 = statepca[5];

  const double d0 = statepca[6];

  const Matrix<double, 3, 1> xyz0 = statepca.tail<3>();

  const double pa0 = std::abs(1./qopa0);
  const double qa = std::copysign(1., qopa0);

  const double pax = pa0*std::cos(lama0)*std::cos(phia0);
  const double pay = pa0*std::cos(lama0)*std::sin(phia0);
  const double paz = pa0*std::sin(lama0);

  const double pb0 = std::abs(1./qopb0);
  const double qb = std::copysign(1., qopb0);

  const double pbx = pb0*std::cos(lamb0)*std::cos(phib0);
  const double pby = pb0*std::cos(lamb0)*std::sin(phib0);
  const double pbz = pb0*std::sin(lamb0);

  const Matrix<double, 3, 1> moma0(pax, pay, paz);
  const Matrix<double, 3, 1> momb0(pbx, pby, pbz);

  const Matrix<double, 3, 1> dvhat = moma0.cross(momb0).normalized();
  const Matrix<double, 3, 1> dv = d0*dvhat;

  const Matrix<double, 3, 1> xyza0 = xyz0 - 0.5*dv;
  const Matrix<double, 3, 1> xyzb0 = xyz0 + 0.5*dv;

  std::array<Matrix<double, 7, 1>, 2> res;

  Matrix<double, 7, 1> &resa = res[0];
  Matrix<double, 7, 1> &resb = res[1];

  resa.head<3>() = xyza0;
  resa.segment<3>(3) = moma0;
  resa[6] = qa;

  resb.head<3>() = xyzb0;
  resb.segment<3>(3) = momb0;
  resb[6] = qb;

  return res;
}

Matrix<double, 10, 1> ResidualGlobalCorrectionMakerBase::twoTrackCart2pca(const Matrix<double, 7, 1> &state0, const Matrix<double, 7, 1> &state1) const {

  const Matrix<double, 3, 1> xyza0 = state0.head<3>();
  const Matrix<double, 3, 1> moma0 = state0.segment<3>(3);

  const Matrix<double, 3, 1> xyzb0 = state1.head<3>();
  const Matrix<double, 3, 1> momb0 = state1.segment<3>(3);

  const Matrix<double, 3, 1> xyz0 = 0.5*(xyza0 + xyzb0);
  const Matrix<double, 3, 1> dv = xyzb0 - xyza0;

  const Matrix<double, 3, 1> dvhat = moma0.cross(momb0).normalized();

  const double d0 = dvhat.dot(dv);

  const double qa = state0[6];

  const double qopa0 = qa/state0.segment<3>(3).norm();
  const double lama0 = std::atan(state0[5]/std::sqrt(state0[3]*state0[3] + state0[4]*state0[4]));
  const double phia0 = std::atan2(state0[4], state0[3]);

  const double qb = state1[6];

  const double qopb0 = qb/state1.segment<3>(3).norm();
  const double lamb0 = std::atan(state1[5]/std::sqrt(state1[3]*state1[3] + state1[4]*state1[4]));
  const double phib0 = std::atan2(state1[4], state1[3]);

  Matrix<double, 10, 1> res;
  res[0] = qopa0;
  res[1] = lama0;
  res[2] = phia0;
  res[3] = qopb0;
  res[4] = lamb0;
  res[5] = phib0;
  res[6] = d0;
  res.tail<3>() = xyz0;

  return res;
}

Matrix<double, 10, 10> ResidualGlobalCorrectionMakerBase::twoTrackPca2curvJacobianD(const Matrix<double, 7, 1> &state0, const Matrix<double, 7, 1> &state1, const MagneticField *field, double dBz0, double dBz1) const {

  const GlobalPoint posa(state0[0], state0[1], state0[2]);
  const GlobalVector &bfielda = field->inInverseGeV(posa);
  const Matrix<double, 3, 1> Bva(bfielda.x(), bfielda.y(), double(bfielda.z()) + 2.99792458e-3*dBz0);

  const GlobalPoint posb(state1[0], state1[1], state1[2]);
  const GlobalVector &bfieldb = field->inInverseGeV(posb);
  const Matrix<double, 3, 1> Bvb(bfieldb.x(), bfieldb.y(), double(bfieldb.z()) + 2.99792458e-3*dBz1);

  const Matrix<double, 3, 1> Wa0 = state0.segment<3>(3).normalized();
  const double Wa0x = Wa0[0];
  const double Wa0y = Wa0[1];
  const double Wa0z = Wa0[2];

  const Matrix<double, 3, 1> Wb0 = state1.segment<3>(3).normalized();
  const double Wb0x = Wb0[0];
  const double Wb0y = Wb0[1];
  const double Wb0z = Wb0[2];

  const Matrix<double, 10, 1> statepca = twoTrackCart2pca(state0, state1);
  const double qopa0 = statepca[0];
  const double lama0 = statepca[1];
  const double phia0 = statepca[2];
  const double qopb0 = statepca[3];
  const double lamb0 = statepca[4];
  const double phib0 = statepca[5];
  const double d0 = statepca[6];
  const double x0 = statepca[7];
  const double y0 = statepca[8];
  const double z0 = statepca[9];

  const double Ba = Bva.norm();
  const Matrix<double, 3, 1> Ha = Bva.normalized();
  const double hax = Ha[0];
  const double hay = Ha[1];
  const double haz = Ha[2];

  const double Bb = Bvb.norm();
  const Matrix<double, 3, 1> Hb = Bvb.normalized();
  const double hbx = Hb[0];
  const double hby = Hb[1];
  const double hbz = Hb[2];

  const double xf0 = std::sin(phib0);
  const double xf1 = std::sin(lama0);
  const double xf2 = std::cos(lamb0);
  const double xf3 = xf1*xf2;
  const double xf4 = xf0*xf3;
  const double xf5 = std::sin(phia0);
  const double xf6 = std::cos(lama0);
  const double xf7 = std::sin(lamb0);
  const double xf8 = xf6*xf7;
  const double xf9 = xf5*xf8;
  const double xf10 = xf4 - xf9;
  const double xf11 = Ba*d0;
  const double xf12 = xf10*xf11;
  const double xf13 = std::pow(xf6, 2);
  const double xf14 = std::pow(xf2, 2);
  const double xf15 = phia0 - phib0;
  const double xf16 = std::sin(xf15);
  const double xf17 = std::pow(xf16, 2);
  const double xf18 = std::cos(phib0);
  const double xf19 = xf18*xf3;
  const double xf20 = std::cos(phia0);
  const double xf21 = xf20*xf8;
  const double xf22 = xf19 - xf21;
  const double xf23 = std::pow(xf10, 2) + xf13*xf14*xf17 + std::pow(xf22, 2);
  const double xf24 = std::sqrt(xf23);
  const double xf25 = 2*xf24;
  const double xf26 = Ba*x0*xf25;
  const double xf27 = 1.0/qopa0;
  const double xf28 = 1.0/Ba;
  const double xf29 = 1.0/xf24;
  const double xf30 = (1.0/2.0)*xf29;
  const double xf31 = xf28*xf30;
  const double xf32 = xf27*xf31;
  const double xf33 = xf32*(xf12 + xf26);
  const double xf34 = qopa0*xf12 + qopa0*xf26;
  const double xf35 = xf31/std::pow(qopa0, 2);
  const double xf36 = xf34*xf35;
  const double xf37 = xf20*xf6;
  const double xf38 = xf11*xf22;
  const double xf39 = xf32*(2*Ba*xf24*y0 - xf38);
  const double xf40 = 2*Ba*qopa0*xf24*y0 - qopa0*xf38;
  const double xf41 = xf35*xf40;
  const double xf42 = xf5*xf6;
  const double xf43 = -xf37*(xf33 - xf36) - xf42*(xf39 - xf41);
  const double xf44 = Ba*qopa0;
  const double xf45 = std::sqrt(xf13)*xf44*(-hax*xf5 + hay*xf20);
  const double xf46 = xf45/xf6;
  const double xf47 = std::pow(xf5, 2);
  const double xf48 = xf13*xf47;
  const double xf49 = std::pow(xf20, 2);
  const double xf50 = xf13*xf49;
  const double xf51 = xf48 + xf50;
  const double xf52 = 1.0/xf51;
  const double xf53 = xf1*xf6;
  const double xf54 = d0*xf30;
  const double xf55 = xf16*xf54;
  const double xf56 = xf1*xf7;
  const double xf57 = xf5*xf56;
  const double xf58 = xf2*xf6;
  const double xf59 = xf0*xf58;
  const double xf60 = (1.0/2.0)*xf10;
  const double xf61 = xf20*xf56;
  const double xf62 = xf18*xf58;
  const double xf63 = (1.0/2.0)*xf22;
  const double xf64 = -xf1*xf14*xf17*xf6 + xf60*(2*xf57 + 2*xf59) + xf63*(2*xf61 + 2*xf62);
  const double xf65 = (1.0/2.0)/std::pow(xf23, 3.0/2.0);
  const double xf66 = -xf64*xf65;
  const double xf67 = d0*xf16*xf58;
  const double xf68 = -xf3*xf55 + xf66*xf67;
  const double xf69 = xf57 + xf59;
  const double xf70 = qopa0*xf11;
  const double xf71 = 2*xf29;
  const double xf72 = xf44*xf71;
  const double xf73 = xf32*(x0*xf64*xf72 + xf69*xf70);
  const double xf74 = xf27*xf28;
  const double xf75 = xf66*xf74;
  const double xf76 = xf34*xf75;
  const double xf77 = xf61 + xf62;
  const double xf78 = xf32*(2*Ba*qopa0*xf29*xf64*y0 - xf70*xf77);
  const double xf79 = xf40*xf75;
  const double xf80 = -xf1*xf68 - xf37*(xf73 + xf76) - xf42*(xf78 + xf79);
  const double xf81 = std::cos(xf15);
  const double xf82 = xf54*xf58*xf81;
  const double xf83 = xf13*xf14*xf16*xf81;
  const double xf84 = -xf10*xf20*xf6*xf7 + xf22*xf9 + xf83;
  const double xf85 = -xf84;
  const double xf86 = xf65*xf67;
  const double xf87 = xf82 + xf85*xf86;
  const double xf88 = xf32*(2*Ba*qopa0*x0*xf29*xf84 - xf21*xf70);
  const double xf89 = xf65*xf74;
  const double xf90 = xf85*xf89;
  const double xf91 = xf34*xf90;
  const double xf92 = xf32*(2*Ba*qopa0*xf29*xf84*y0 - xf70*xf9);
  const double xf93 = xf40*xf90;
  const double xf94 = -xf1*xf87 - xf37*(xf88 + xf91) - xf42*(xf92 + xf93);
  const double xf95 = xf0*xf56;
  const double xf96 = xf2*xf42;
  const double xf97 = xf18*xf56;
  const double xf98 = xf2*xf37;
  const double xf99 = -xf13*xf17*xf2*xf7 + xf60*(-2*xf95 - 2*xf96) + xf63*(-2*xf97 - 2*xf98);
  const double xf100 = -xf99;
  const double xf101 = xf100*xf86 - xf55*xf8;
  const double xf102 = xf100*xf89;
  const double xf103 = xf102*xf34;
  const double xf104 = -xf95 - xf96;
  const double xf105 = xf32*(x0*xf72*xf99 + xf104*xf70);
  const double xf106 = xf102*xf40;
  const double xf107 = -xf97 - xf98;
  const double xf108 = xf32*(2*Ba*qopa0*xf29*xf99*y0 - xf107*xf70);
  const double xf109 = -xf1*xf101 - xf37*(xf103 + xf105) - xf42*(xf106 + xf108);
  const double xf110 = -xf10*xf19 + xf22*xf4 + xf83;
  const double xf111 = xf110*xf86 - xf82;
  const double xf112 = -xf110;
  const double xf113 = xf112*xf72;
  const double xf114 = xf32*(x0*xf113 + xf19*xf70);
  const double xf115 = xf110*xf89;
  const double xf116 = xf115*xf34;
  const double xf117 = xf32*(xf113*y0 + xf4*xf70);
  const double xf118 = xf115*xf40;
  const double xf119 = -xf1*xf111 - xf37*(xf114 + xf116) - xf42*(xf117 + xf118);
  const double xf120 = xf16*xf30;
  const double xf121 = xf29*xf60;
  const double xf122 = -xf120*xf3*xf6 - xf121*xf37 + (1.0/2.0)*xf22*xf29*xf5*xf6;
  const double xf123 = std::tan(lama0);
  const double xf124 = xf44*(hax*xf123*xf20 + hay*xf123*xf5 - haz);
  const double xf125 = std::pow(Wa0x, 2);
  const double xf126 = std::pow(Wa0y, 2);
  const double xf127 = xf125 + xf126;
  const double xf128 = std::pow(xf127, -1.0/2.0);
  const double xf129 = Wa0y*xf128;
  const double xf130 = Wa0x*xf128;
  const double xf131 = Wa0x*xf5 - Wa0y*xf20;
  const double xf132 = xf128*xf131;
  const double xf133 = xf132*xf6;
  const double xf134 = xf29*xf63;
  const double xf135 = xf13*xf132;
  const double xf136 = Wa0z*xf130;
  const double xf137 = Wa0z*xf129;
  const double xf138 = -Wa0x*Wa0z*xf37 - Wa0y*Wa0z*xf42 + xf1*xf127;
  const double xf139 = xf128*xf138;
  const double xf140 = xf125*xf128 + xf126*xf128;
  const double xf141 = xf120*xf58;
  const double xf142 = -xf68;
  const double xf143 = Bb*d0;
  const double xf144 = qopb0*xf143;
  const double xf145 = Bb*qopb0;
  const double xf146 = xf145*xf71;
  const double xf147 = 1.0/qopb0;
  const double xf148 = 1.0/Bb;
  const double xf149 = xf148*xf30;
  const double xf150 = xf147*xf149;
  const double xf151 = xf150*(xf144*xf77 + xf146*xf64*y0);
  const double xf152 = xf143*xf22;
  const double xf153 = Bb*xf25*y0;
  const double xf154 = qopb0*xf152 + qopb0*xf153;
  const double xf155 = xf147*xf148;
  const double xf156 = xf155*xf66;
  const double xf157 = xf154*xf156;
  const double xf158 = xf0*xf2;
  const double xf159 = xf150*(2*Bb*qopb0*x0*xf29*xf64 - xf144*xf69);
  const double xf160 = xf10*xf143;
  const double xf161 = 2*Bb*qopb0*x0*xf24 - qopb0*xf160;
  const double xf162 = xf156*xf161;
  const double xf163 = xf18*xf2;
  const double xf164 = -xf142*xf7 - xf158*(xf151 + xf157) - xf163*(xf159 + xf162);
  const double xf165 = std::sqrt(xf14)*xf145*(-hbx*xf0 + hby*xf18);
  const double xf166 = xf165/xf2;
  const double xf167 = -xf87;
  const double xf168 = xf146*xf84;
  const double xf169 = xf150*(xf144*xf9 + xf168*y0);
  const double xf170 = xf155*xf65;
  const double xf171 = xf170*xf85;
  const double xf172 = xf154*xf171;
  const double xf173 = xf150*(x0*xf168 + xf144*xf21);
  const double xf174 = xf161*xf171;
  const double xf175 = -xf158*(xf169 + xf172) - xf163*(xf173 + xf174) - xf167*xf7;
  const double xf176 = xf150*(xf152 + xf153);
  const double xf177 = xf149/std::pow(qopb0, 2);
  const double xf178 = xf154*xf177;
  const double xf179 = xf150*(2*Bb*x0*xf24 - xf160);
  const double xf180 = xf161*xf177;
  const double xf181 = -xf158*(xf176 - xf178) - xf163*(xf179 - xf180);
  const double xf182 = std::pow(xf0, 2);
  const double xf183 = xf14*xf182;
  const double xf184 = std::pow(xf18, 2);
  const double xf185 = xf14*xf184;
  const double xf186 = xf183 + xf185;
  const double xf187 = 1.0/xf186;
  const double xf188 = xf2*xf7;
  const double xf189 = -xf101;
  const double xf190 = xf100*xf170;
  const double xf191 = xf154*xf190;
  const double xf192 = xf150*(xf107*xf144 + xf146*xf99*y0);
  const double xf193 = xf161*xf190;
  const double xf194 = xf150*(2*Bb*qopb0*x0*xf29*xf99 - xf104*xf144);
  const double xf195 = -xf158*(xf191 + xf192) - xf163*(xf193 + xf194) - xf189*xf7;
  const double xf196 = -xf111;
  const double xf197 = xf150*(2*Bb*qopb0*xf112*xf29*y0 - xf144*xf4);
  const double xf198 = xf110*xf170;
  const double xf199 = xf154*xf198;
  const double xf200 = xf150*(2*Bb*qopb0*x0*xf112*xf29 - xf144*xf19);
  const double xf201 = xf161*xf198;
  const double xf202 = -xf158*(xf197 + xf199) - xf163*(xf200 + xf201) - xf196*xf7;
  const double xf203 = xf120*xf2*xf8 + xf121*xf163 - xf134*xf158;
  const double xf204 = std::tan(lamb0);
  const double xf205 = xf145*(hbx*xf18*xf204 + hby*xf0*xf204 - hbz);
  const double xf206 = std::pow(Wb0x, 2);
  const double xf207 = std::pow(Wb0y, 2);
  const double xf208 = xf206 + xf207;
  const double xf209 = std::pow(xf208, -1.0/2.0);
  const double xf210 = Wb0x*xf209;
  const double xf211 = Wb0y*xf209;
  const double xf212 = xf209*(Wb0x*xf0 - Wb0y*xf18);
  const double xf213 = xf2*xf212;
  const double xf214 = xf14*xf212;
  const double xf215 = Wb0z*xf211;
  const double xf216 = Wb0z*xf210;
  const double xf217 = xf206*xf209 + xf207*xf209;
  const double xf218 = -Wb0x*Wb0z*xf163 - Wb0y*Wb0z*xf158 + xf208*xf7;
  const double xf219 = xf209*xf218;
  const double dqopadqopa0 = 1;
  const double dqopadlama0 = 0;
  const double dqopadphia0 = 0;
  const double dqopadqopb0 = 0;
  const double dqopadlamb0 = 0;
  const double dqopadphib0 = 0;
  const double dqopadd0 = 0;
  const double dqopadx0 = 0;
  const double dqopady0 = 0;
  const double dqopadz0 = 0;
  const double dlamadqopa0 = xf43*xf46;
  const double dlamadlama0 = xf46*xf80 + (xf1*(xf47*xf53 + xf49*xf53)/std::pow(xf51, 3.0/2.0) + xf6/std::sqrt(xf51))/(std::pow(xf1, 2)*xf52 + 1);
  const double dlamadphia0 = xf46*xf94;
  const double dlamadqopb0 = 0;
  const double dlamadlamb0 = xf109*xf46;
  const double dlamadphib0 = xf119*xf46;
  const double dlamadd0 = xf122*xf46;
  const double dlamadx0 = -xf20*xf45;
  const double dlamady0 = -xf45*xf5;
  const double dlamadz0 = -xf1*xf46;
  const double dphiadqopa0 = xf124*xf43;
  const double dphiadlama0 = xf124*xf80;
  const double dphiadphia0 = xf124*xf94 + xf48*xf52 + xf50*xf52;
  const double dphiadqopb0 = 0;
  const double dphiadlamb0 = xf109*xf124;
  const double dphiadphib0 = xf119*xf124;
  const double dphiadd0 = xf122*xf124;
  const double dphiadx0 = -xf124*xf37;
  const double dphiady0 = -xf124*xf42;
  const double dphiadz0 = -xf1*xf124;
  const double dxtadqopa0 = -xf129*xf33 + xf129*xf36 + xf130*xf39 - xf130*xf41 + xf133*xf43;
  const double dxtadlama0 = -xf129*xf73 - xf129*xf76 + xf130*xf78 + xf130*xf79 + xf133*xf80;
  const double dxtadphia0 = -xf129*xf88 - xf129*xf91 + xf130*xf92 + xf130*xf93 + xf133*xf94;
  const double dxtadqopb0 = 0;
  const double dxtadlamb0 = -xf103*xf129 - xf105*xf129 + xf106*xf130 + xf108*xf130 + xf109*xf133;
  const double dxtadphib0 = -xf114*xf129 - xf116*xf129 + xf117*xf130 + xf118*xf130 + xf119*xf133;
  const double dxtadd0 = -xf121*xf129 + xf122*xf128*xf131*xf6 - xf130*xf134;
  const double dxtadx0 = -xf129 - xf135*xf20;
  const double dxtady0 = xf130 - xf135*xf5;
  const double dxtadz0 = -xf132*xf53;
  const double dytadqopa0 = -xf136*xf33 + xf136*xf36 - xf137*xf39 + xf137*xf41 + xf139*xf43;
  const double dytadlama0 = xf128*xf138*xf80 - xf136*xf73 - xf136*xf76 - xf137*xf78 - xf137*xf79 + xf140*xf68;
  const double dytadphia0 = xf128*xf138*xf94 - xf136*xf88 - xf136*xf91 - xf137*xf92 - xf137*xf93 + xf140*xf87;
  const double dytadqopb0 = 0;
  const double dytadlamb0 = xf101*xf140 - xf103*xf136 - xf105*xf136 - xf106*xf137 - xf108*xf137 + xf109*xf128*xf138;
  const double dytadphib0 = xf111*xf140 - xf114*xf136 - xf116*xf136 - xf117*xf137 - xf118*xf137 + xf119*xf128*xf138;
  const double dytadd0 = -xf121*xf136 + xf122*xf139 + xf134*xf137 + xf140*xf141;
  const double dytadx0 = -xf136 - xf139*xf37;
  const double dytady0 = -xf137 - xf139*xf42;
  const double dytadz0 = -xf1*xf139 + xf140;
  const double dqopbdqopa0 = 0;
  const double dqopbdlama0 = 0;
  const double dqopbdphia0 = 0;
  const double dqopbdqopb0 = 1;
  const double dqopbdlamb0 = 0;
  const double dqopbdphib0 = 0;
  const double dqopbdd0 = 0;
  const double dqopbdx0 = 0;
  const double dqopbdy0 = 0;
  const double dqopbdz0 = 0;
  const double dlambdqopa0 = 0;
  const double dlambdlama0 = xf164*xf166;
  const double dlambdphia0 = xf166*xf175;
  const double dlambdqopb0 = xf166*xf181;
  const double dlambdlamb0 = xf166*xf195 + (xf2/std::sqrt(xf186) + xf7*(xf182*xf188 + xf184*xf188)/std::pow(xf186, 3.0/2.0))/(xf187*std::pow(xf7, 2) + 1);
  const double dlambdphib0 = xf166*xf202;
  const double dlambdd0 = xf166*xf203;
  const double dlambdx0 = -xf165*xf18;
  const double dlambdy0 = -xf0*xf165;
  const double dlambdz0 = -xf166*xf7;
  const double dphibdqopa0 = 0;
  const double dphibdlama0 = xf164*xf205;
  const double dphibdphia0 = xf175*xf205;
  const double dphibdqopb0 = xf181*xf205;
  const double dphibdlamb0 = xf195*xf205;
  const double dphibdphib0 = xf183*xf187 + xf185*xf187 + xf202*xf205;
  const double dphibdd0 = xf203*xf205;
  const double dphibdx0 = -xf163*xf205;
  const double dphibdy0 = -xf158*xf205;
  const double dphibdz0 = -xf205*xf7;
  const double dxtbdqopa0 = 0;
  const double dxtbdlama0 = xf151*xf210 + xf157*xf210 - xf159*xf211 - xf162*xf211 + xf164*xf213;
  const double dxtbdphia0 = xf169*xf210 + xf172*xf210 - xf173*xf211 - xf174*xf211 + xf175*xf213;
  const double dxtbdqopb0 = xf176*xf210 - xf178*xf210 - xf179*xf211 + xf180*xf211 + xf181*xf213;
  const double dxtbdlamb0 = xf191*xf210 + xf192*xf210 - xf193*xf211 - xf194*xf211 + xf195*xf213;
  const double dxtbdphib0 = xf197*xf210 + xf199*xf210 - xf200*xf211 - xf201*xf211 + xf202*xf213;
  const double dxtbdd0 = xf121*xf211 + xf134*xf210 + xf203*xf213;
  const double dxtbdx0 = -xf18*xf214 - xf211;
  const double dxtbdy0 = -xf0*xf214 + xf210;
  const double dxtbdz0 = -xf188*xf212;
  const double dytbdqopa0 = 0;
  const double dytbdlama0 = xf142*xf217 - xf151*xf215 - xf157*xf215 - xf159*xf216 - xf162*xf216 + xf164*xf209*xf218;
  const double dytbdphia0 = xf167*xf217 - xf169*xf215 - xf172*xf215 - xf173*xf216 - xf174*xf216 + xf175*xf209*xf218;
  const double dytbdqopb0 = -xf176*xf215 + xf178*xf215 - xf179*xf216 + xf180*xf216 + xf181*xf219;
  const double dytbdlamb0 = xf189*xf217 - xf191*xf215 - xf192*xf215 - xf193*xf216 - xf194*xf216 + xf195*xf209*xf218;
  const double dytbdphib0 = xf196*xf217 - xf197*xf215 - xf199*xf215 - xf200*xf216 - xf201*xf216 + xf202*xf209*xf218;
  const double dytbdd0 = xf121*xf216 - xf134*xf215 - xf141*xf217 + xf203*xf219;
  const double dytbdx0 = -xf163*xf219 - xf216;
  const double dytbdy0 = -xf158*xf219 - xf215;
  const double dytbdz0 = xf217 - xf219*xf7;
  Matrix<double, 10, 10> res;
  res(0,0) = dqopadqopa0;
  res(0,1) = dqopadlama0;
  res(0,2) = dqopadphia0;
  res(0,3) = dqopadqopb0;
  res(0,4) = dqopadlamb0;
  res(0,5) = dqopadphib0;
  res(0,6) = dqopadd0;
  res(0,7) = dqopadx0;
  res(0,8) = dqopady0;
  res(0,9) = dqopadz0;
  res(1,0) = dlamadqopa0;
  res(1,1) = dlamadlama0;
  res(1,2) = dlamadphia0;
  res(1,3) = dlamadqopb0;
  res(1,4) = dlamadlamb0;
  res(1,5) = dlamadphib0;
  res(1,6) = dlamadd0;
  res(1,7) = dlamadx0;
  res(1,8) = dlamady0;
  res(1,9) = dlamadz0;
  res(2,0) = dphiadqopa0;
  res(2,1) = dphiadlama0;
  res(2,2) = dphiadphia0;
  res(2,3) = dphiadqopb0;
  res(2,4) = dphiadlamb0;
  res(2,5) = dphiadphib0;
  res(2,6) = dphiadd0;
  res(2,7) = dphiadx0;
  res(2,8) = dphiady0;
  res(2,9) = dphiadz0;
  res(3,0) = dxtadqopa0;
  res(3,1) = dxtadlama0;
  res(3,2) = dxtadphia0;
  res(3,3) = dxtadqopb0;
  res(3,4) = dxtadlamb0;
  res(3,5) = dxtadphib0;
  res(3,6) = dxtadd0;
  res(3,7) = dxtadx0;
  res(3,8) = dxtady0;
  res(3,9) = dxtadz0;
  res(4,0) = dytadqopa0;
  res(4,1) = dytadlama0;
  res(4,2) = dytadphia0;
  res(4,3) = dytadqopb0;
  res(4,4) = dytadlamb0;
  res(4,5) = dytadphib0;
  res(4,6) = dytadd0;
  res(4,7) = dytadx0;
  res(4,8) = dytady0;
  res(4,9) = dytadz0;
  res(5,0) = dqopbdqopa0;
  res(5,1) = dqopbdlama0;
  res(5,2) = dqopbdphia0;
  res(5,3) = dqopbdqopb0;
  res(5,4) = dqopbdlamb0;
  res(5,5) = dqopbdphib0;
  res(5,6) = dqopbdd0;
  res(5,7) = dqopbdx0;
  res(5,8) = dqopbdy0;
  res(5,9) = dqopbdz0;
  res(6,0) = dlambdqopa0;
  res(6,1) = dlambdlama0;
  res(6,2) = dlambdphia0;
  res(6,3) = dlambdqopb0;
  res(6,4) = dlambdlamb0;
  res(6,5) = dlambdphib0;
  res(6,6) = dlambdd0;
  res(6,7) = dlambdx0;
  res(6,8) = dlambdy0;
  res(6,9) = dlambdz0;
  res(7,0) = dphibdqopa0;
  res(7,1) = dphibdlama0;
  res(7,2) = dphibdphia0;
  res(7,3) = dphibdqopb0;
  res(7,4) = dphibdlamb0;
  res(7,5) = dphibdphib0;
  res(7,6) = dphibdd0;
  res(7,7) = dphibdx0;
  res(7,8) = dphibdy0;
  res(7,9) = dphibdz0;
  res(8,0) = dxtbdqopa0;
  res(8,1) = dxtbdlama0;
  res(8,2) = dxtbdphia0;
  res(8,3) = dxtbdqopb0;
  res(8,4) = dxtbdlamb0;
  res(8,5) = dxtbdphib0;
  res(8,6) = dxtbdd0;
  res(8,7) = dxtbdx0;
  res(8,8) = dxtbdy0;
  res(8,9) = dxtbdz0;
  res(9,0) = dytbdqopa0;
  res(9,1) = dytbdlama0;
  res(9,2) = dytbdphia0;
  res(9,3) = dytbdqopb0;
  res(9,4) = dytbdlamb0;
  res(9,5) = dytbdphib0;
  res(9,6) = dytbdd0;
  res(9,7) = dytbdx0;
  res(9,8) = dytbdy0;
  res(9,9) = dytbdz0;

  return res;
}

Matrix<double, 5, 5> ResidualGlobalCorrectionMakerBase::curv2localJacobianAltelossD(const Matrix<double, 7, 1> &state, const MagneticField *field, const GloballyPositioned<double> &surface, double dEdx, double mass, double dBz) const {
  
  
  const GlobalPoint pos(state[0], state[1], state[2]);  
  const GlobalVector &bfield = field->inInverseGeV(pos);
  const Matrix<double, 3, 1> Bv(bfield.x(), bfield.y(), double(bfield.z()) + 2.99792458e-3*dBz);
  
  const double q = state[6];
  
  const double qop0 = q/state.segment<3>(3).norm();
  const double lam0 = std::atan(state[5]/std::sqrt(state[3]*state[3] + state[4]*state[4]));
  const double phi0 = std::atan2(state[4], state[3]);
  
  const Matrix<double, 3, 1> W0 = state.segment<3>(3).normalized();
  const double W0x = W0[0];
  const double W0y = W0[1];
  const double W0z = W0[2];
    
  const Vector3DBase<double, LocalTag> lx(1.,0.,0.);
  const Vector3DBase<double, LocalTag> ly(0.,1.,0.);
  const Vector3DBase<double, LocalTag> lz(0.,0.,1.);
  const Point3DBase<double, LocalTag> l0(0., 0., 0.);
  
//   const Vector3DBase<double, GlobalTag> I1 = surface.toGlobal<double>(lz);
//   const Vector3DBase<double, GlobalTag> J1 = surface.toGlobal<double>(lx);
//   const Vector3DBase<double, GlobalTag> K1 = surface.toGlobal<double>(ly);  
//   const Point3DBase<double, GlobalTag> r1 = surface.toGlobal<double>(l0);
  
  const Vector3DBase<double, GlobalTag> I1 = surface.toGlobal(lz);
  const Vector3DBase<double, GlobalTag> J1 = surface.toGlobal(lx);
  const Vector3DBase<double, GlobalTag> K1 = surface.toGlobal(ly);  
  const Point3DBase<double, GlobalTag> r1 = surface.toGlobal(l0);
  
//   const Vector3DBase<double, GlobalTag> I1 = toGlobal(surface, lz);
//   const Vector3DBase<double, GlobalTag> J1 = toGlobal(surface, lx);
//   const Vector3DBase<double, GlobalTag> K1 = toGlobal(surface, ly);  
//   const Point3DBase<double, GlobalTag> r1 = toGlobal(surface, l0);
  
  const double Ix1 = I1.x();
  const double Iy1 = I1.y();
  const double Iz1 = I1.z();
  
  const double Jx1 = J1.x();
  const double Jy1 = J1.y();
  const double Jz1 = J1.z();
  
  const double Kx1 = K1.x();
  const double Ky1 = K1.y();
  const double Kz1 = K1.z();
  
  const double rx1 = r1.x();
  const double ry1 = r1.y();
  const double rz1 = r1.z();
    
  const double B = Bv.norm();
  const Matrix<double, 3, 1> H = Bv.normalized();
  const double hx = H[0];
  const double hy = H[1];
  const double hz = H[2];

  const double x0 = std::pow(q, 2);
  const double x1 = std::pow(x0, 3.0/2.0);
  const double x2 = Ix1*W0y - Iy1*W0x;
  const double x3 = std::pow(qop0, 2);
  const double x4 = std::pow(W0x, 2) + std::pow(W0y, 2);
  const double x5 = std::pow(x4, -1.0/2.0);
  const double x6 = std::sin(lam0);
  const double x7 = Iz1*x6;
  const double x8 = std::cos(lam0);
  const double x9 = std::cos(phi0);
  const double x10 = Ix1*x9;
  const double x11 = std::sin(phi0);
  const double x12 = Iy1*x11;
  const double x13 = x10*x8 + x12*x8 + x7;
  const double x14 = x5/x13;
  const double x15 = dEdx*q*x14*x3*std::sqrt(std::pow(mass, 2)*x3 + x0)/x1;
  const double x16 = Ix1*W0x*W0z + Iy1*W0y*W0z - Iz1*x4;
  const double x17 = std::pow(x13, -2);
  const double x18 = Iz1*Jx1;
  const double x19 = Iz1*Jy1;
  const double x20 = lam0 + phi0;
  const double x21 = -phi0;
  const double x22 = lam0 + x21;
  const double x23 = std::pow(Ix1*std::cos(x20) + Ix1*std::cos(x22) + Iy1*std::sin(x20) - Iy1*std::sin(x22) + 2*x7, -2);
  const double x24 = 2*Ix1;
  const double x25 = Jy1*x24;
  const double x26 = 2*Iy1;
  const double x27 = Jx1*x26;
  const double x28 = 2*lam0;
  const double x29 = std::cos(x28);
  const double x30 = phi0 + x28;
  const double x31 = std::cos(x30);
  const double x32 = std::sin(x30);
  const double x33 = Ix1*Jz1;
  const double x34 = Iy1*Jz1;
  const double x35 = x21 + x28;
  const double x36 = std::cos(x35);
  const double x37 = std::sin(x35);
  const double x38 = x8*x9;
  const double x39 = x11*x8;
  const double x40 = Jx1*x38 + Jy1*x39 + Jz1*x6;
  const double x41 = hz*x8;
  const double x42 = hy*x6 - x11*x41;
  const double x43 = x8*(hx*x11 - hy*x9);
  const double x44 = hx*x6 - x41*x9;
  const double x45 = -Ix1*x42 + Iy1*x44 - Iz1*x43;
  const double x46 = B*qop0*x5/std::pow(x13, 3);
  const double x47 = x46*(x13*(-Jx1*x42 + Jy1*x44 - Jz1*x43) - x40*x45);
  const double x48 = Iz1*Kx1;
  const double x49 = Iz1*Ky1;
  const double x50 = Ky1*x24;
  const double x51 = Kx1*x26;
  const double x52 = Ix1*Kz1;
  const double x53 = Iy1*Kz1;
  const double x54 = Kx1*x38 + Ky1*x39 + Kz1*x6;
  const double x55 = x46*(x13*(-Kx1*x42 + Ky1*x44 - Kz1*x43) - x45*x54);
  const double x56 = x13*x4;
  const double x57 = W0z*x13;
  const double dqopdqop0 = x1*std::fabs(qop0)/(std::pow(q, 3)*qop0);
  const double dqopdlam0 = 0;
  const double dqopdphi0 = 0;
  const double dqopdxt0 = -x15*x2;
  const double dqopdyt0 = -x15*x16;
  const double ddxdzdqop0 = 0;
  const double ddxdzdlam0 = x17*(Jz1*x10 + Jz1*x12 - x11*x19 - x18*x9);
  const double ddxdzdphi0 = x23*(x18*x31 - x18*x36 + x19*x32 + x19*x37 + x25*x29 + x25 - x27*x29 - x27 - x31*x33 - x32*x34 + x33*x36 - x34*x37);
  const double ddxdzdxt0 = x2*x47;
  const double ddxdzdyt0 = x16*x47;
  const double ddydzdqop0 = 0;
  const double ddydzdlam0 = x17*(Kz1*x10 + Kz1*x12 - x11*x49 - x48*x9);
  const double ddydzdphi0 = x23*(x29*x50 - x29*x51 + x31*x48 - x31*x52 + x32*x49 - x32*x53 - x36*x48 + x36*x52 + x37*x49 - x37*x53 + x50 - x51);
  const double ddydzdxt0 = x2*x55;
  const double ddydzdyt0 = x16*x55;
  const double dxdqop0 = 0;
  const double dxdlam0 = 0;
  const double dxdphi0 = 0;
  const double dxdxt0 = x14*(x13*(-Jx1*W0y + Jy1*W0x) + x2*x40);
  const double dxdyt0 = x14*(Jz1*x56 + x16*x40 - x57*(Jx1*W0x + Jy1*W0y));
  const double dydqop0 = 0;
  const double dydlam0 = 0;
  const double dydphi0 = 0;
  const double dydxt0 = x14*(x13*(-Kx1*W0y + Ky1*W0x) + x2*x54);
  const double dydyt0 = x14*(Kz1*x56 + x16*x54 - x57*(Kx1*W0x + Ky1*W0y));
  Matrix<double, 5, 5> res;
  res(0,0) = dqopdqop0;
  res(0,1) = dqopdlam0;
  res(0,2) = dqopdphi0;
  res(0,3) = dqopdxt0;
  res(0,4) = dqopdyt0;
  res(1,0) = ddxdzdqop0;
  res(1,1) = ddxdzdlam0;
  res(1,2) = ddxdzdphi0;
  res(1,3) = ddxdzdxt0;
  res(1,4) = ddxdzdyt0;
  res(2,0) = ddydzdqop0;
  res(2,1) = ddydzdlam0;
  res(2,2) = ddydzdphi0;
  res(2,3) = ddydzdxt0;
  res(2,4) = ddydzdyt0;
  res(3,0) = dxdqop0;
  res(3,1) = dxdlam0;
  res(3,2) = dxdphi0;
  res(3,3) = dxdxt0;
  res(3,4) = dxdyt0;
  res(4,0) = dydqop0;
  res(4,1) = dydlam0;
  res(4,2) = dydphi0;
  res(4,3) = dydxt0;
  res(4,4) = dydyt0;
  
  return res;
                                              
}

Matrix<double, 5, 5> ResidualGlobalCorrectionMakerBase::curv2localJacobianSimpleD(const Matrix<double, 7, 1> &state, const MagneticField *field, const GloballyPositioned<double> &surface, double dEdx, double mass, double dBz) const {
  
  
  const GlobalPoint pos(state[0], state[1], state[2]);  
  const GlobalVector &bfield = field->inInverseGeV(pos);
  const Matrix<double, 3, 1> Bv(bfield.x(), bfield.y(), double(bfield.z()) + 2.99792458e-3*dBz);
  
  const double q = state[6];
  
  const double qop0 = q/state.segment<3>(3).norm();
  const double lam0 = std::atan(state[5]/std::sqrt(state[3]*state[3] + state[4]*state[4]));
  const double phi0 = std::atan2(state[4], state[3]);
  
  const Matrix<double, 3, 1> W0 = state.segment<3>(3).normalized();
  const double W0x = W0[0];
  const double W0y = W0[1];
  const double W0z = W0[2];
    
  const Vector3DBase<double, LocalTag> lx(1.,0.,0.);
  const Vector3DBase<double, LocalTag> ly(0.,1.,0.);
  const Vector3DBase<double, LocalTag> lz(0.,0.,1.);
  const Point3DBase<double, LocalTag> l0(0., 0., 0.);
  
//   const Vector3DBase<double, GlobalTag> I1 = surface.toGlobal<double>(lz);
//   const Vector3DBase<double, GlobalTag> J1 = surface.toGlobal<double>(lx);
//   const Vector3DBase<double, GlobalTag> K1 = surface.toGlobal<double>(ly);  
//   const Point3DBase<double, GlobalTag> r1 = surface.toGlobal<double>(l0);
  
  const Vector3DBase<double, GlobalTag> I1 = surface.toGlobal(lz);
  const Vector3DBase<double, GlobalTag> J1 = surface.toGlobal(lx);
  const Vector3DBase<double, GlobalTag> K1 = surface.toGlobal(ly);  
  const Point3DBase<double, GlobalTag> r1 = surface.toGlobal(l0);
  
//   const Vector3DBase<double, GlobalTag> I1 = toGlobal(surface, lz);
//   const Vector3DBase<double, GlobalTag> J1 = toGlobal(surface, lx);
//   const Vector3DBase<double, GlobalTag> K1 = toGlobal(surface, ly);  
//   const Point3DBase<double, GlobalTag> r1 = toGlobal(surface, l0);
  
  const double Ix1 = I1.x();
  const double Iy1 = I1.y();
  const double Iz1 = I1.z();
  
  const double Jx1 = J1.x();
  const double Jy1 = J1.y();
  const double Jz1 = J1.z();
  
  const double Kx1 = K1.x();
  const double Ky1 = K1.y();
  const double Kz1 = K1.z();
  
  const double rx1 = r1.x();
  const double ry1 = r1.y();
  const double rz1 = r1.z();
    
  const double B = Bv.norm();
  const Matrix<double, 3, 1> H = Bv.normalized();
  const double hx = H[0];
  const double hy = H[1];
  const double hz = H[2];

  const double x0 = std::pow(q, 2);
  const double x1 = std::pow(x0, 3.0/2.0);
  const double x2 = Ix1*W0y - Iy1*W0x;
  const double x3 = std::pow(qop0, 2);
  const double x4 = std::pow(W0x, 2) + std::pow(W0y, 2);
  const double x5 = std::pow(x4, -1.0/2.0);
  const double x6 = std::sin(lam0);
  const double x7 = Iz1*x6;
  const double x8 = std::cos(lam0);
  const double x9 = std::cos(phi0);
  const double x10 = Ix1*x9;
  const double x11 = std::sin(phi0);
  const double x12 = Iy1*x11;
  const double x13 = x10*x8 + x12*x8 + x7;
  const double x14 = x5/x13;
  const double x15 = dEdx*q*x14*x3*std::sqrt(std::pow(mass, 2)*x3 + x0)/x1;
  const double x16 = Ix1*W0x*W0z + Iy1*W0y*W0z - Iz1*x4;
  const double x17 = std::pow(x13, -2);
  const double x18 = Iz1*Jx1;
  const double x19 = Iz1*Jy1;
  const double x20 = lam0 + phi0;
  const double x21 = -phi0;
  const double x22 = lam0 + x21;
  const double x23 = std::pow(Ix1*std::cos(x20) + Ix1*std::cos(x22) + Iy1*std::sin(x20) - Iy1*std::sin(x22) + 2*x7, -2);
  const double x24 = 2*Ix1;
  const double x25 = Jy1*x24;
  const double x26 = 2*Iy1;
  const double x27 = Jx1*x26;
  const double x28 = 2*lam0;
  const double x29 = std::cos(x28);
  const double x30 = phi0 + x28;
  const double x31 = std::cos(x30);
  const double x32 = std::sin(x30);
  const double x33 = Ix1*Jz1;
  const double x34 = Iy1*Jz1;
  const double x35 = x21 + x28;
  const double x36 = std::cos(x35);
  const double x37 = std::sin(x35);
  const double x38 = x8*x9;
  const double x39 = x11*x8;
  const double x40 = Jx1*x38 + Jy1*x39 + Jz1*x6;
  const double x41 = hz*x8;
  const double x42 = hy*x6 - x11*x41;
  const double x43 = x8*(hx*x11 - hy*x9);
  const double x44 = hx*x6 - x41*x9;
  const double x45 = -Ix1*x42 + Iy1*x44 - Iz1*x43;
  const double x46 = B*qop0*x5/std::pow(x13, 3);
  const double x47 = x46*(x13*(-Jx1*x42 + Jy1*x44 - Jz1*x43) - x40*x45);
  const double x48 = Iz1*Kx1;
  const double x49 = Iz1*Ky1;
  const double x50 = Ky1*x24;
  const double x51 = Kx1*x26;
  const double x52 = Ix1*Kz1;
  const double x53 = Iy1*Kz1;
  const double x54 = Kx1*x38 + Ky1*x39 + Kz1*x6;
  const double x55 = x46*(x13*(-Kx1*x42 + Ky1*x44 - Kz1*x43) - x45*x54);
  const double x56 = x13*x4;
  const double x57 = W0z*x13;
  const double dqopdqop0 = x1*std::fabs(qop0)/(std::pow(q, 3)*qop0);
  const double dqopdlam0 = 0;
  const double dqopdphi0 = 0;
  const double dqopdxt0 = -x15*x2;
  const double dqopdyt0 = -x15*x16;
  const double ddxdzdqop0 = 0;
  const double ddxdzdlam0 = x17*(Jz1*x10 + Jz1*x12 - x11*x19 - x18*x9);
  const double ddxdzdphi0 = x23*(x18*x31 - x18*x36 + x19*x32 + x19*x37 + x25*x29 + x25 - x27*x29 - x27 - x31*x33 - x32*x34 + x33*x36 - x34*x37);
  const double ddxdzdxt0 = x2*x47;
  const double ddxdzdyt0 = x16*x47;
  const double ddydzdqop0 = 0;
  const double ddydzdlam0 = x17*(Kz1*x10 + Kz1*x12 - x11*x49 - x48*x9);
  const double ddydzdphi0 = x23*(x29*x50 - x29*x51 + x31*x48 - x31*x52 + x32*x49 - x32*x53 - x36*x48 + x36*x52 + x37*x49 - x37*x53 + x50 - x51);
  const double ddydzdxt0 = x2*x55;
  const double ddydzdyt0 = x16*x55;
  const double dxdqop0 = 0;
  const double dxdlam0 = 0;
  const double dxdphi0 = 0;
  const double dxdxt0 = x14*(x13*(-Jx1*W0y + Jy1*W0x) + x2*x40);
  const double dxdyt0 = x14*(Jz1*x56 + x16*x40 - x57*(Jx1*W0x + Jy1*W0y));
  const double dydqop0 = 0;
  const double dydlam0 = 0;
  const double dydphi0 = 0;
  const double dydxt0 = x14*(x13*(-Kx1*W0y + Ky1*W0x) + x2*x54);
  const double dydyt0 = x14*(Kz1*x56 + x16*x54 - x57*(Kx1*W0x + Ky1*W0y));
  Matrix<double, 5, 5> res;
  res(0,0) = dqopdqop0;
  res(0,1) = dqopdlam0;
  res(0,2) = dqopdphi0;
  res(0,3) = dqopdxt0;
  res(0,4) = dqopdyt0;
  res(1,0) = ddxdzdqop0;
  res(1,1) = ddxdzdlam0;
  res(1,2) = ddxdzdphi0;
  res(1,3) = ddxdzdxt0;
  res(1,4) = ddxdzdyt0;
  res(2,0) = ddydzdqop0;
  res(2,1) = ddydzdlam0;
  res(2,2) = ddydzdphi0;
  res(2,3) = ddydzdxt0;
  res(2,4) = ddydzdyt0;
  res(3,0) = dxdqop0;
  res(3,1) = dxdlam0;
  res(3,2) = dxdphi0;
  res(3,3) = dxdxt0;
  res(3,4) = dxdyt0;
  res(4,0) = dydqop0;
  res(4,1) = dydlam0;
  res(4,2) = dydphi0;
  res(4,3) = dydxt0;
  res(4,4) = dydyt0;
  
  return res;
                                              
}

Matrix<double, 6, 5> ResidualGlobalCorrectionMakerBase::curv2cartJacobianAltD(const Matrix<double, 7, 1> &state) const {
  
  const double charge = state[6];
  
  const double qop0 = charge/state.segment<3>(3).norm();
  const double lam0 = std::atan(state[5]/std::sqrt(state[3]*state[3] + state[4]*state[4]));
  const double phi0 = std::atan2(state[4], state[3]);
  
  const Matrix<double, 3, 1> W0 = state.segment<3>(3).normalized();
  const double W0x = W0[0];
  const double W0y = W0[1];
  const double W0z = W0[2];
    
  const double x0 = std::sqrt(std::pow(W0x, 2) + std::pow(W0y, 2));
  const double x1 = 1.0/x0;
  const double x2 = W0y*x1;
  const double x3 = W0x*x1;
  const double x4 = 1.0/qop0;
  const double x5 = std::cos(phi0);
  const double x6 = 1.0/std::fabs(qop0);
  const double x7 = x6*std::cos(lam0);
  const double x8 = x5*x7;
  const double x9 = x6*std::sin(lam0);
  const double x10 = std::sin(phi0);
  const double x11 = x10*x7;
  const double dxdqop0 = 0;
  const double dxdlam0 = 0;
  const double dxdphi0 = 0;
  const double dxdxt0 = -x2;
  const double dxdyt0 = -W0z*x3;
  const double dydqop0 = 0;
  const double dydlam0 = 0;
  const double dydphi0 = 0;
  const double dydxt0 = x3;
  const double dydyt0 = -W0z*x2;
  const double dzdqop0 = 0;
  const double dzdlam0 = 0;
  const double dzdphi0 = 0;
  const double dzdxt0 = 0;
  const double dzdyt0 = x0;
  const double dpxdqop0 = -x4*x8;
  const double dpxdlam0 = -x5*x9;
  const double dpxdphi0 = -x11;
  const double dpxdxt0 = 0;
  const double dpxdyt0 = 0;
  const double dpydqop0 = -x11*x4;
  const double dpydlam0 = -x10*x9;
  const double dpydphi0 = x8;
  const double dpydxt0 = 0;
  const double dpydyt0 = 0;
  const double dpzdqop0 = -x4*x9;
  const double dpzdlam0 = x7;
  const double dpzdphi0 = 0;
  const double dpzdxt0 = 0;
  const double dpzdyt0 = 0;
  Matrix<double, 6, 5> res;
  res(0,0) = dxdqop0;
  res(0,1) = dxdlam0;
  res(0,2) = dxdphi0;
  res(0,3) = dxdxt0;
  res(0,4) = dxdyt0;
  res(1,0) = dydqop0;
  res(1,1) = dydlam0;
  res(1,2) = dydphi0;
  res(1,3) = dydxt0;
  res(1,4) = dydyt0;
  res(2,0) = dzdqop0;
  res(2,1) = dzdlam0;
  res(2,2) = dzdphi0;
  res(2,3) = dzdxt0;
  res(2,4) = dzdyt0;
  res(3,0) = dpxdqop0;
  res(3,1) = dpxdlam0;
  res(3,2) = dpxdphi0;
  res(3,3) = dpxdxt0;
  res(3,4) = dpxdyt0;
  res(4,0) = dpydqop0;
  res(4,1) = dpydlam0;
  res(4,2) = dpydphi0;
  res(4,3) = dpydxt0;
  res(4,4) = dpydyt0;
  res(5,0) = dpzdqop0;
  res(5,1) = dpzdlam0;
  res(5,2) = dpzdphi0;
  res(5,3) = dpzdxt0;
  res(5,4) = dpzdyt0;
  
  return res;

  
}



Matrix<double, 2, 1> ResidualGlobalCorrectionMakerBase::localPositionConvolutionD(const Matrix<double, 7, 1>& state, const Matrix<double, 5, 5> &curvcov, const GloballyPositioned<double> &surface) const {

  const double q = state[6];

  const double qop = q/state.segment<3>(3).norm();
  const double lam = std::atan(state[5]/std::sqrt(state[3]*state[3] + state[4]*state[4]));
  const double phi = std::atan2(state[4], state[3]);
  const double xt = 0.;
  const double yt = 0.;

  // curvilinear parameters
//   const CurvilinearTrajectoryParameters curv(tsos.globalPosition(), tsos.globalMomentum(), tsos.charge());
//   const double qop = curv.Qbp();
//   const double lam = curv.lambda();
//   const double phi = curv.phi();
//   const double xt = curv.xT();
//   const double yt = curv.yT();
//   const double xt = 0.;
//   const double yt = 0.;

  const Matrix<double, 3, 1> p0(state[3], state[4], state[5]);
  const Matrix<double, 3, 1> W0 = p0.normalized();
  const Matrix<double, 3, 1> zhat(0., 0., 1.);

  const Matrix<double, 3, 1> U0 = zhat.cross(W0).normalized();
  const Matrix<double, 3, 1> V0 = W0.cross(U0);

//   std::cout << "U0" << std::endl;
//   std::cout << U0 << std::endl;
//   std::cout << "V0" << std::endl;
//   std::cout << V0 << std::endl;

//   const Matrix<double, 3, 1> x0alt = xt*U0 + yt*V0;
//   std::cout << "global pos" << std::endl;
//   std::cout << tsos.globalPosition() << std::endl;
//   std::cout << "x0alt" << std::endl;
//   std::cout << x0alt << std::endl;
//   std::cout << "xt" << std::endl;
//   std::cout << xt << std::endl;
//   std::cout << "yt" << std::endl;
//   std::cout << yt << std::endl;

  const Vector3DBase<double, LocalTag> lx(1.,0.,0.);
  const Vector3DBase<double, LocalTag> ly(0.,1.,0.);
  const Vector3DBase<double, LocalTag> lz(0.,0.,1.);
  const Vector3DBase<double, GlobalTag> I = surface.toGlobal(lz);
  const Vector3DBase<double, GlobalTag> J = surface.toGlobal(lx);
  const Vector3DBase<double, GlobalTag> K = surface.toGlobal(ly);

  const Point3DBase<double, LocalTag> l0(0., 0.);
  const Point3DBase<double, GlobalTag> r = surface.toGlobal(l0);

  const double Ux = U0[0];
  const double Uy = U0[1];
//   const double Uz = U0[2];

  const double Vx = V0[0];
  const double Vy = V0[1];
  const double Vz = V0[2];

  const double Ix = I.x();
  const double Iy = I.y();
  const double Iz = I.z();

  const double Jx = J.x();
  const double Jy = J.y();
  const double Jz = J.z();

  const double Kx = K.x();
  const double Ky = K.y();
  const double Kz = K.z();

  const double rx = r.x();
  const double ry = r.y();
  const double rz = r.z();

  const double pos0x = state[0];
  const double pos0y = state[1];
  const double pos0z = state[2];

  //sympy stuff goes here
  const double x0 = std::sin(lam);
  const double x1 = Iz*x0;
  const double x2 = std::cos(lam);
  const double x3 = std::cos(phi);
  const double x4 = Ix*x3;
  const double x5 = x2*x4;
  const double x6 = std::sin(phi);
  const double x7 = Iy*x6;
  const double x8 = x2*x7;
  const double x9 = x5 + x8;
  const double x10 = x1 + x9;
  const double x11 = 1.0/x10;
  const double x12 = Ix*rx;
  const double x13 = Iy*ry;
  const double x14 = Iz*rz;
  const double x15 = Ix*pos0x;
  const double x16 = Iy*pos0y;
  const double x17 = Iz*pos0z;
  const double x18 = Ux*xt;
  const double x19 = Ix*x18;
  const double x20 = Vx*yt;
  const double x21 = Ix*x20;
  const double x22 = Uy*xt;
  const double x23 = Iy*x22;
  const double x24 = Vy*yt;
  const double x25 = Iy*x24;
  const double x26 = Vz*yt;
  const double x27 = Iz*x26;
  const double x28 = x12 + x13 + x14 - x15 - x16 - x17 - x19 - x21 - x23 - x25 - x27;
  const double x29 = pos0z + x26;
  const double x30 = Iz*x2;
  const double x31 = x0*x4;
  const double x32 = x0*x7;
  const double x33 = x30 - x31 - x32;
  const double x34 = x2*x28;
  const double x35 = x29*x33 + x34;
  const double x36 = Jz*x11;
  const double x37 = pos0x + x18 + x20;
  const double x38 = x0*x28;
  const double x39 = -x3*x38 + x33*x37;
  const double x40 = Jx*x11;
  const double x41 = pos0y + x22 + x24;
  const double x42 = x33*x41 - x38*x6;
  const double x43 = Jy*x11;
  const double x44 = std::pow(x10, -2);
  const double x45 = -x30 + x31 + x32;
  const double x46 = x44*x45;
  const double x47 = x10*x29 + x38;
  const double x48 = Jz*x47;
  const double x49 = x3*x34;
  const double x50 = x10*x37 + x49;
  const double x51 = Jx*x50;
  const double x52 = x34*x6;
  const double x53 = x10*x41 + x52;
  const double x54 = Jy*x46;
  const double x55 = -x5 - x8;
  const double x56 = -x1 + x55;
  const double x57 = -x12 - x13 - x14 + x15 + x16 + x17 + x19 + x21 + x23 + x25 + x27;
  const double x58 = x2*x57;
  const double x59 = x3*x58;
  const double x60 = x37*x56 + x59;
  const double x61 = x41*x56 + x58*x6;
  const double x62 = x29*x56 - x38;
  const double x63 = Jz*x46;
  const double x64 = 2*x35;
  const double x65 = Jx*x46;
  const double x66 = 2*x39;
  const double x67 = 2*x42;
  const double x68 = std::pow(x10, -3);
  const double x69 = x45*x68;
  const double x70 = x69*(-2*x30 + 2*x31 + 2*x32);
  const double x71 = Jy*x53;
  const double x72 = x0*x6;
  const double x73 = Ix*x72;
  const double x74 = x0*x3;
  const double x75 = Iy*x74;
  const double x76 = x73 - x75;
  const double x77 = x29*x36;
  const double x78 = x2*x3;
  const double x79 = Iy*x78;
  const double x80 = x2*x6;
  const double x81 = Ix*x80;
  const double x82 = x79 - x81;
  const double x83 = x29*x82;
  const double x84 = x41*x76 + x57*x74;
  const double x85 = x37*x76 - x57*x72;
  const double x86 = x44*(-x73 + x75);
  const double x87 = -x79 + x81;
  const double x88 = x44*x87;
  const double x89 = Jz*x88;
  const double x90 = x41*x82 + x49;
  const double x91 = -x52;
  const double x92 = x37*x82 + x91;
  const double x93 = Jx*x88;
  const double x94 = Jy*x88;
  const double x95 = -2*x79 + 2*x81;
  const double x96 = x69*x95;
  const double x97 = Ix*Ux;
  const double x98 = Iy*Uy;
  const double x99 = -x97 - x98;
  const double x100 = x36*x99;
  const double x101 = x0*x99;
  const double x102 = x97 + x98;
  const double x103 = Ux*x33 + x102*x74;
  const double x104 = Uy*x33 + x102*x72;
  const double x105 = x78*x99;
  const double x106 = Ux*x10 + x105;
  const double x107 = Uy*x10 + x80*x99;
  const double x108 = Ix*Vx;
  const double x109 = Iy*Vy;
  const double x110 = Iz*Vz;
  const double x111 = x108 + x109 + x110;
  const double x112 = Vx*x33 + x111*x74;
  const double x113 = Vy*x33 + x111*x72;
  const double x114 = -x108 - x109 - x110;
  const double x115 = x114*x2;
  const double x116 = Vz*x33 + x115;
  const double x117 = Vz*x10 + x0*x114;
  const double x118 = x117*x46;
  const double x119 = x115*x3;
  const double x120 = Vx*x10 + x119;
  const double x121 = Vy*x10 + x115*x6;
  const double x122 = 2*x83;
  const double x123 = x37*x55 + x59;
  const double x124 = x41*x55 + x91;
  const double x125 = x44*x9;
  const double x126 = 2*x90;
  const double x127 = 2*x92;
  const double x128 = x68*x87*x95;
  const double x129 = Ux*x82 + x102*x80;
  const double x130 = Uy*x82 + x105;
  const double x131 = Vz*x82;
  const double x132 = Vx*x82 + x111*x80;
  const double x133 = Vy*x82 + x119;
  const double x134 = Kz*x11;
  const double x135 = Kx*x11;
  const double x136 = Ky*x11;
  const double x137 = Kz*x47;
  const double x138 = Kx*x46;
  const double x139 = Ky*x46;
  const double x140 = Kz*x46;
  const double x141 = Kx*x50;
  const double x142 = Ky*x53;
  const double x143 = x134*x29;
  const double x144 = Kz*x88;
  const double x145 = Kx*x88;
  const double x146 = Ky*x88;
  const double x147 = x134*x99;
  const double shat = x11*x28;
  const double dvdqop = 0;
  const double d2vdqopdqop = 0;
  const double d2vdqopdlam = 0;
  const double d2vdqopdphi = 0;
  const double d2vdqopdxt = 0;
  const double d2vdqopdyt = 0;
  const double dvdlam = x35*x36 + x39*x40 + x42*x43 + x46*x48 + x46*x51 + x53*x54;
  const double d2vdlamdlam = x36*x47 + x36*x62 + x40*x50 + x40*x60 + x43*x53 + x43*x61 + x48*x70 + x51*x70 + x54*x67 + x63*x64 + x65*x66 + x70*x71;
  const double d2vdlamdphi = x35*x89 + x39*x93 + x40*x85 + x42*x94 + x43*x84 + x48*x86 + x48*x96 + x51*x86 + x51*x96 + x54*x90 + x63*x83 + x65*x92 + x71*x86 + x71*x96 + x76*x77;
  const double d2vdlamdxt = x100*x2 + x101*x63 + x103*x40 + x104*x43 + x106*x65 + x107*x54;
  const double d2vdlamdyt = Jz*x118 + x112*x40 + x113*x43 + x116*x36 + x120*x65 + x121*x54;
  const double dvdphi = x40*x92 + x43*x90 + x48*x88 + x51*x88 + x53*x94 + x77*x82;
  const double d2vdphidphi = x122*x89 + x123*x40 + x124*x43 + x125*x48 + x125*x51 + x125*x71 + x126*x94 + x127*x93 + x128*x48 + x128*x51 + x128*x71 + x55*x77;
  const double d2vdphidxt = x101*x89 + x106*x93 + x107*x94 + x129*x40 + x130*x43;
  const double d2vdphidyt = x117*x89 + x120*x93 + x121*x94 + x131*x36 + x132*x40 + x133*x43;
  const double dvdxt = x0*x100 + x106*x40 + x107*x43;
  const double d2vdxtdxt = 0;
  const double d2vdxtdyt = 0;
  const double dvdyt = x117*x36 + x120*x40 + x121*x43;
  const double d2vdytdyt = 0;
  const double dwdqop = 0;
  const double d2wdqopdqop = 0;
  const double d2wdqopdlam = 0;
  const double d2wdqopdphi = 0;
  const double d2wdqopdxt = 0;
  const double d2wdqopdyt = 0;
  const double dwdlam = x134*x35 + x135*x39 + x136*x42 + x137*x46 + x138*x50 + x139*x53;
  const double d2wdlamdlam = x134*x47 + x134*x62 + x135*x50 + x135*x60 + x136*x53 + x136*x61 + x137*x70 + x138*x66 + x139*x67 + x140*x64 + x141*x70 + x142*x70;
  const double d2wdlamdphi = x135*x85 + x136*x84 + x137*x86 + x137*x96 + x138*x92 + x139*x90 + x140*x83 + x141*x86 + x141*x96 + x142*x86 + x142*x96 + x143*x76 + x144*x35 + x145*x39 + x146*x42;
  const double d2wdlamdxt = x101*x140 + x103*x135 + x104*x136 + x106*x138 + x107*x139 + x147*x2;
  const double d2wdlamdyt = Kz*x118 + x112*x135 + x113*x136 + x116*x134 + x120*x138 + x121*x139;
  const double dwdphi = x135*x92 + x136*x90 + x137*x88 + x143*x82 + x145*x50 + x146*x53;
  const double d2wdphidphi = x122*x144 + x123*x135 + x124*x136 + x125*x137 + x125*x141 + x125*x142 + x126*x146 + x127*x145 + x128*x137 + x128*x141 + x128*x142 + x143*x55;
  const double d2wdphidxt = x101*x144 + x106*x145 + x107*x146 + x129*x135 + x130*x136;
  const double d2wdphidyt = x117*x144 + x120*x145 + x121*x146 + x131*x134 + x132*x135 + x133*x136;
  const double dwdxt = x0*x147 + x106*x135 + x107*x136;
  const double d2wdxtdxt = 0;
  const double d2wdxtdyt = 0;
  const double dwdyt = x117*x134 + x120*x135 + x121*x136;
  const double d2wdytdyt = 0;
  Matrix<double, 5, 5> d2vdx2;
  d2vdx2(0, 0) = d2vdqopdqop;
  d2vdx2(0, 1) = d2vdqopdlam;
  d2vdx2(0, 2) = d2vdqopdphi;
  d2vdx2(0, 3) = d2vdqopdxt;
  d2vdx2(0, 4) = d2vdqopdyt;
  d2vdx2(1, 0) = d2vdqopdlam;
  d2vdx2(1, 1) = d2vdlamdlam;
  d2vdx2(1, 2) = d2vdlamdphi;
  d2vdx2(1, 3) = d2vdlamdxt;
  d2vdx2(1, 4) = d2vdlamdyt;
  d2vdx2(2, 0) = d2vdqopdphi;
  d2vdx2(2, 1) = d2vdlamdphi;
  d2vdx2(2, 2) = d2vdphidphi;
  d2vdx2(2, 3) = d2vdphidxt;
  d2vdx2(2, 4) = d2vdphidyt;
  d2vdx2(3, 0) = d2vdqopdxt;
  d2vdx2(3, 1) = d2vdlamdxt;
  d2vdx2(3, 2) = d2vdphidxt;
  d2vdx2(3, 3) = d2vdxtdxt;
  d2vdx2(3, 4) = d2vdxtdyt;
  d2vdx2(4, 0) = d2vdqopdyt;
  d2vdx2(4, 1) = d2vdlamdyt;
  d2vdx2(4, 2) = d2vdphidyt;
  d2vdx2(4, 3) = d2vdxtdyt;
  d2vdx2(4, 4) = d2vdytdyt;
  Matrix<double, 5, 5> d2wdx2;
  d2wdx2(0, 0) = d2wdqopdqop;
  d2wdx2(0, 1) = d2wdqopdlam;
  d2wdx2(0, 2) = d2wdqopdphi;
  d2wdx2(0, 3) = d2wdqopdxt;
  d2wdx2(0, 4) = d2wdqopdyt;
  d2wdx2(1, 0) = d2wdqopdlam;
  d2wdx2(1, 1) = d2wdlamdlam;
  d2wdx2(1, 2) = d2wdlamdphi;
  d2wdx2(1, 3) = d2wdlamdxt;
  d2wdx2(1, 4) = d2wdlamdyt;
  d2wdx2(2, 0) = d2wdqopdphi;
  d2wdx2(2, 1) = d2wdlamdphi;
  d2wdx2(2, 2) = d2wdphidphi;
  d2wdx2(2, 3) = d2wdphidxt;
  d2wdx2(2, 4) = d2wdphidyt;
  d2wdx2(3, 0) = d2wdqopdxt;
  d2wdx2(3, 1) = d2wdlamdxt;
  d2wdx2(3, 2) = d2wdphidxt;
  d2wdx2(3, 3) = d2wdxtdxt;
  d2wdx2(3, 4) = d2wdxtdyt;
  d2wdx2(4, 0) = d2wdqopdyt;
  d2wdx2(4, 1) = d2wdlamdyt;
  d2wdx2(4, 2) = d2wdphidyt;
  d2wdx2(4, 3) = d2wdxtdyt;
  d2wdx2(4, 4) = d2wdytdyt;
  Matrix<double, 5, 1> dvdx;
  dvdx[0] = dvdqop;
  dvdx[1] = dvdlam;
  dvdx[2] = dvdphi;
  dvdx[3] = dvdxt;
  dvdx[4] = dvdyt;
  Matrix<double, 5, 1> dwdx;
  dwdx[0] = dwdqop;
  dwdx[1] = dwdlam;
  dwdx[2] = dwdphi;
  dwdx[3] = dwdxt;
  dwdx[4] = dwdyt;


  Matrix<double, 2, 1> res;
  res[0] = 0.5*(d2vdx2*curvcov).trace();
  res[1] = 0.5*(d2wdx2*curvcov).trace();


  return res;

}

Matrix<double, 1, 6> ResidualGlobalCorrectionMakerBase::massJacobianAltD(const Matrix<double, 7, 1> &state0, const Matrix<double, 7, 1> &state1, double dmass) const {
  
  const double mp = dmass;
  
  const double qop0 = state0[6]/state0.segment<3>(3).norm();
  const double lam0 = std::atan(state0[5]/std::sqrt(state0[3]*state0[3] + state0[4]*state0[4]));
  const double phi0 = std::atan2(state0[4], state0[3]);
  
  const double qop1 = state1[6]/state1.segment<3>(3).norm();
  const double lam1 = std::atan(state1[5]/std::sqrt(state1[3]*state1[3] + state1[4]*state1[4]));
  const double phi1 = std::atan2(state1[4], state1[3]);

  const double xf0 = std::pow(mp, 2);
  const double xf1 = std::sin(lam0);
  const double xf2 = std::sin(lam1);
  const double xf3 = 1.0/std::fabs(qop0);
  const double xf4 = 1.0/std::fabs(qop1);
  const double xf5 = xf3*xf4;
  const double xf6 = xf2*xf5;
  const double xf7 = std::sin(phi0);
  const double xf8 = std::sin(phi1);
  const double xf9 = xf7*xf8;
  const double xf10 = std::cos(lam0);
  const double xf11 = std::cos(lam1);
  const double xf12 = xf11*xf5;
  const double xf13 = xf10*xf12;
  const double xf14 = 2*xf13;
  const double xf15 = std::cos(phi0);
  const double xf16 = std::cos(phi1);
  const double xf17 = xf15*xf16;
  const double xf18 = std::pow(qop0, -2);
  const double xf19 = std::sqrt(xf0 + xf18);
  const double xf20 = std::pow(qop1, -2);
  const double xf21 = std::sqrt(xf0 + xf20);
  const double xf22 = std::sqrt(2*xf0 - 2*xf1*xf6 - xf14*xf17 - xf14*xf9 + 2*xf19*xf21);
  const double xf23 = xf1*xf2;
  const double xf24 = xf18*xf4*(((qop0) > 0) - ((qop0) < 0));
  const double xf25 = xf10*xf11;
  const double xf26 = xf24*xf25;
  const double xf27 = 1.0/xf22;
  const double xf28 = xf10*xf6;
  const double xf29 = xf1*xf12;
  const double xf30 = xf13*xf16*xf7;
  const double xf31 = xf13*xf15*xf8;
  const double xf32 = xf20*xf3*(((qop1) > 0) - ((qop1) < 0));
  const double xf33 = xf25*xf32;
  const double m = xf22;
  const double dmdqop0 = xf27*(xf17*xf26 + xf23*xf24 + xf26*xf9 - xf21/(std::pow(qop0, 3)*xf19));
  const double dmdlam0 = xf27*(xf17*xf29 - xf28 + xf29*xf9);
  const double dmdphi0 = xf27*(xf30 - xf31);
  const double dmdqop1 = xf27*(xf17*xf33 + xf23*xf32 + xf33*xf9 - xf19/(std::pow(qop1, 3)*xf21));
  const double dmdlam1 = xf27*(xf17*xf28 + xf28*xf9 - xf29);
  const double dmdphi1 = xf27*(-xf30 + xf31);
  Matrix<double, 1, 6> res;
  res(0,0) = dmdqop0;
  res(0,1) = dmdlam0;
  res(0,2) = dmdphi0;
  res(0,3) = dmdqop1;
  res(0,4) = dmdlam1;
  res(0,5) = dmdphi1;
  
//   std::cout << "massJacobianAlt m = " << m << std::endl;




  
  return res;
}

Matrix<double, 6, 6> ResidualGlobalCorrectionMakerBase::massHessianAltD(const Matrix<double, 7, 1> &state0, const Matrix<double, 7, 1> &state1, double dmass) const {
  
  const double mp = dmass;
  
  const double qop0 = state0[6]/state0.segment<3>(3).norm();
  const double lam0 = std::atan(state0[5]/std::sqrt(state0[3]*state0[3] + state0[4]*state0[4]));
  const double phi0 = std::atan2(state0[4], state0[3]);
  
  const double qop1 = state1[6]/state1.segment<3>(3).norm();
  const double lam1 = std::atan(state1[5]/std::sqrt(state1[3]*state1[3] + state1[4]*state1[4]));
  const double phi1 = std::atan2(state1[4], state1[3]);

  const double xf0 = std::pow(mp, 2);
  const double xf1 = std::sin(lam0);
  const double xf2 = std::sin(lam1);
  const double xf3 = 1.0/std::fabs(qop0);
  const double xf4 = 1.0/std::fabs(qop1);
  const double xf5 = xf3*xf4;
  const double xf6 = xf2*xf5;
  const double xf7 = xf1*xf6;
  const double xf8 = std::cos(lam0);
  const double xf9 = std::cos(lam1);
  const double xf10 = xf5*xf9;
  const double xf11 = xf10*xf8;
  const double xf12 = std::sin(phi0);
  const double xf13 = std::sin(phi1);
  const double xf14 = xf12*xf13;
  const double xf15 = xf11*xf14;
  const double xf16 = std::cos(phi0);
  const double xf17 = std::cos(phi1);
  const double xf18 = xf16*xf17;
  const double xf19 = xf11*xf18;
  const double xf20 = std::pow(qop0, -2);
  const double xf21 = xf0 + xf20;
  const double xf22 = std::sqrt(xf21);
  const double xf23 = std::pow(qop1, -2);
  const double xf24 = xf0 + xf23;
  const double xf25 = std::sqrt(xf24);
  const double xf26 = 2*xf0 - 2*xf15 - 2*xf19 + 2*xf22*xf25 - 2*xf7;
  const double xf27 = std::pow(xf26, -1.0/2.0);
  const double xf28 = xf1*xf2;
  const double xf29 = 2*xf28;
  const double xf30 = std::pow(qop0, -3);
  const double xf31 = (((qop0) > 0) - ((qop0) < 0));
  const double xf32 = xf31*xf4;
  const double xf33 = xf30*xf32;
  const double xf34 = 2*xf33;
  const double xf35 = xf8*xf9;
  const double xf36 = xf14*xf35;
  const double xf37 = xf18*xf35;
  const double xf38 = 1.0/xf22;
  const double xf39 = xf20*xf32;
  const double xf40 = xf35*xf39;
  const double xf41 = xf14*xf40 + xf18*xf40 - xf25*xf30*xf38 + xf28*xf39;
  const double xf42 = -xf41;
  const double xf43 = std::pow(xf26, -3.0/2.0);
  const double xf44 = xf41*xf43;
  const double xf45 = xf1*xf9;
  const double xf46 = xf39*xf45;
  const double xf47 = xf27*(-xf14*xf46 - xf18*xf46 + xf2*xf20*xf31*xf4*xf8);
  const double xf48 = xf1*xf10;
  const double xf49 = xf14*xf48 + xf18*xf48 - xf2*xf3*xf4*xf8;
  const double xf50 = -xf49;
  const double xf51 = xf12*xf17;
  const double xf52 = -xf13*xf16*xf20*xf31*xf4*xf8*xf9 + xf40*xf51;
  const double xf53 = -xf27*xf52;
  const double xf54 = xf11*xf51 - xf13*xf16*xf3*xf4*xf8*xf9;
  const double xf55 = -xf54;
  const double xf56 = (((qop1) > 0) - ((qop1) < 0));
  const double xf57 = xf23*xf56;
  const double xf58 = xf20*xf31*xf57;
  const double xf59 = std::pow(qop1, -3);
  const double xf60 = 1.0/xf25;
  const double xf61 = xf27*(-xf28*xf58 + xf30*xf38*xf59*xf60 - xf36*xf58 - xf37*xf58);
  const double xf62 = xf3*xf57;
  const double xf63 = xf35*xf62;
  const double xf64 = xf14*xf63 + xf18*xf63 - xf22*xf59*xf60 + xf28*xf62;
  const double xf65 = -xf64;
  const double xf66 = xf2*xf8;
  const double xf67 = xf39*xf66;
  const double xf68 = xf27*(xf1*xf20*xf31*xf4*xf9 - xf14*xf67 - xf18*xf67);
  const double xf69 = xf6*xf8;
  const double xf70 = -xf1*xf3*xf4*xf9 + xf14*xf69 + xf18*xf69;
  const double xf71 = -xf70;
  const double xf72 = xf27*xf52;
  const double xf73 = xf43*xf49;
  const double xf74 = xf15 + xf19;
  const double xf75 = xf27*(xf7 + xf74);
  const double xf76 = -xf1*xf13*xf16*xf3*xf4*xf9 + xf48*xf51;
  const double xf77 = -xf27*xf76;
  const double xf78 = xf45*xf62;
  const double xf79 = xf27*(-xf14*xf78 - xf18*xf78 + xf2*xf23*xf3*xf56*xf8);
  const double xf80 = xf27*(-xf11 - xf14*xf7 - xf18*xf7);
  const double xf81 = xf27*xf76;
  const double xf82 = xf43*xf54;
  const double xf83 = xf27*xf74 + xf55*xf82;
  const double xf84 = -xf13*xf16*xf23*xf3*xf56*xf8*xf9 + xf51*xf63;
  const double xf85 = -xf27*xf84;
  const double xf86 = -xf13*xf16*xf2*xf3*xf4*xf8 + xf51*xf69;
  const double xf87 = -xf27*xf86;
  const double xf88 = -xf27*xf74;
  const double xf89 = xf43*xf64;
  const double xf90 = xf3*xf56*xf59;
  const double xf91 = 2*xf90;
  const double xf92 = xf62*xf66;
  const double xf93 = xf27*(xf1*xf23*xf3*xf56*xf9 - xf14*xf92 - xf18*xf92);
  const double xf94 = xf27*xf84;
  const double xf95 = xf43*xf70;
  const double xf96 = xf27*xf86;
  const double xf97 = xf43*xf55;
  const double d2mdqop0dqop0 = xf27*(-xf29*xf33 - xf34*xf36 - xf34*xf37 + 3*xf25*xf38/std::pow(qop0, 4) - xf25/(std::pow(qop0, 6)*std::pow(xf21, 3.0/2.0))) + xf42*xf44;
  const double d2mdqop0dlam0 = xf44*xf50 + xf47;
  const double d2mdqop0dphi0 = xf44*xf55 + xf53;
  const double d2mdqop0dqop1 = xf44*xf65 + xf61;
  const double d2mdqop0dlam1 = xf44*xf71 + xf68;
  const double d2mdqop0dphi1 = xf44*xf54 + xf72;
  const double d2mdlam0dqop0 = xf42*xf73 + xf47;
  const double d2mdlam0dlam0 = xf50*xf73 + xf75;
  const double d2mdlam0dphi0 = xf55*xf73 + xf77;
  const double d2mdlam0dqop1 = xf65*xf73 + xf79;
  const double d2mdlam0dlam1 = xf71*xf73 + xf80;
  const double d2mdlam0dphi1 = xf54*xf73 + xf81;
  const double d2mdphi0dqop0 = xf42*xf82 + xf53;
  const double d2mdphi0dlam0 = xf50*xf82 + xf77;
  const double d2mdphi0dphi0 = xf83;
  const double d2mdphi0dqop1 = xf65*xf82 + xf85;
  const double d2mdphi0dlam1 = xf71*xf82 + xf87;
  const double d2mdphi0dphi1 = xf43*std::pow(xf54, 2) + xf88;
  const double d2mdqop1dqop0 = xf42*xf89 + xf61;
  const double d2mdqop1dlam0 = xf50*xf89 + xf79;
  const double d2mdqop1dphi0 = xf55*xf89 + xf85;
  const double d2mdqop1dqop1 = xf27*(-xf29*xf90 - xf36*xf91 - xf37*xf91 + 3*xf22*xf60/std::pow(qop1, 4) - xf22/(std::pow(qop1, 6)*std::pow(xf24, 3.0/2.0))) + xf65*xf89;
  const double d2mdqop1dlam1 = xf71*xf89 + xf93;
  const double d2mdqop1dphi1 = xf64*xf82 + xf94;
  const double d2mdlam1dqop0 = xf42*xf95 + xf68;
  const double d2mdlam1dlam0 = xf50*xf95 + xf80;
  const double d2mdlam1dphi0 = xf55*xf95 + xf87;
  const double d2mdlam1dqop1 = xf65*xf95 + xf93;
  const double d2mdlam1dlam1 = xf71*xf95 + xf75;
  const double d2mdlam1dphi1 = xf70*xf82 + xf96;
  const double d2mdphi1dqop0 = xf42*xf97 + xf72;
  const double d2mdphi1dlam0 = xf50*xf97 + xf81;
  const double d2mdphi1dphi0 = xf43*std::pow(xf55, 2) + xf88;
  const double d2mdphi1dqop1 = xf65*xf97 + xf94;
  const double d2mdphi1dlam1 = xf71*xf97 + xf96;
  const double d2mdphi1dphi1 = xf83;
  Matrix<double, 6, 6> res;
  res(0,0) = d2mdqop0dqop0;
  res(0,1) = d2mdqop0dlam0;
  res(0,2) = d2mdqop0dphi0;
  res(0,3) = d2mdqop0dqop1;
  res(0,4) = d2mdqop0dlam1;
  res(0,5) = d2mdqop0dphi1;
  res(1,0) = d2mdlam0dqop0;
  res(1,1) = d2mdlam0dlam0;
  res(1,2) = d2mdlam0dphi0;
  res(1,3) = d2mdlam0dqop1;
  res(1,4) = d2mdlam0dlam1;
  res(1,5) = d2mdlam0dphi1;
  res(2,0) = d2mdphi0dqop0;
  res(2,1) = d2mdphi0dlam0;
  res(2,2) = d2mdphi0dphi0;
  res(2,3) = d2mdphi0dqop1;
  res(2,4) = d2mdphi0dlam1;
  res(2,5) = d2mdphi0dphi1;
  res(3,0) = d2mdqop1dqop0;
  res(3,1) = d2mdqop1dlam0;
  res(3,2) = d2mdqop1dphi0;
  res(3,3) = d2mdqop1dqop1;
  res(3,4) = d2mdqop1dlam1;
  res(3,5) = d2mdqop1dphi1;
  res(4,0) = d2mdlam1dqop0;
  res(4,1) = d2mdlam1dlam0;
  res(4,2) = d2mdlam1dphi0;
  res(4,3) = d2mdlam1dqop1;
  res(4,4) = d2mdlam1dlam1;
  res(4,5) = d2mdlam1dphi1;
  res(5,0) = d2mdphi1dqop0;
  res(5,1) = d2mdphi1dlam0;
  res(5,2) = d2mdphi1dphi0;
  res(5,3) = d2mdphi1dqop1;
  res(5,4) = d2mdphi1dlam1;
  res(5,5) = d2mdphi1dphi1;


  
  return res;
}

Matrix<double, 1, 6> ResidualGlobalCorrectionMakerBase::massinvsqJacobianAltD(const Matrix<double, 7, 1> &state0, const Matrix<double, 7, 1> &state1, double dmass) const {
  
  const double mp = dmass;
  
  const double qop0 = state0[6]/state0.segment<3>(3).norm();
  const double lam0 = std::atan(state0[5]/std::sqrt(state0[3]*state0[3] + state0[4]*state0[4]));
  const double phi0 = std::atan2(state0[4], state0[3]);
  
  const double qop1 = state1[6]/state1.segment<3>(3).norm();
  const double lam1 = std::atan(state1[5]/std::sqrt(state1[3]*state1[3] + state1[4]*state1[4]));
  const double phi1 = std::atan2(state1[4], state1[3]);

  const double xf0 = std::sin(lam0);
  const double xf1 = std::sin(lam1);
  const double xf2 = xf0*xf1;
  const double xf3 = std::pow(qop0, -2);
  const double xf4 = 1.0/std::fabs(qop1);
  const double xf5 = 2*xf4;
  const double xf6 = xf3*xf5*(((qop0) > 0) - ((qop0) < 0));
  const double xf7 = std::sin(phi0);
  const double xf8 = std::sin(phi1);
  const double xf9 = xf7*xf8;
  const double xf10 = std::cos(lam0);
  const double xf11 = std::cos(lam1);
  const double xf12 = xf10*xf11;
  const double xf13 = xf12*xf6;
  const double xf14 = std::cos(phi0);
  const double xf15 = std::cos(phi1);
  const double xf16 = xf14*xf15;
  const double xf17 = std::pow(mp, 2);
  const double xf18 = std::sqrt(xf17 + xf3);
  const double xf19 = std::pow(qop1, -2);
  const double xf20 = std::sqrt(xf17 + xf19);
  const double xf21 = 1.0/std::fabs(qop0);
  const double xf22 = xf21*xf5;
  const double xf23 = xf1*xf22;
  const double xf24 = xf11*xf22;
  const double xf25 = xf10*xf24;
  const double xf26 = 1.0/std::pow(-xf0*xf23 - xf16*xf25 + 2*xf17 + 2*xf18*xf20 - xf25*xf9, 2);
  const double xf27 = xf0*xf24;
  const double xf28 = -2*xf10*xf11*xf14*xf21*xf4*xf8 + xf15*xf25*xf7;
  const double xf29 = 2*xf19*xf21*(((qop1) > 0) - ((qop1) < 0));
  const double xf30 = xf12*xf29;
  const double xf31 = xf10*xf23;
  const double dminvsqdqop0 = xf26*(-xf13*xf16 - xf13*xf9 - xf2*xf6 + 2*xf20/(std::pow(qop0, 3)*xf18));
  const double dminvsqdlam0 = xf26*(2*xf1*xf10*xf21*xf4 - xf16*xf27 - xf27*xf9);
  const double dminvsqdphi0 = -xf26*xf28;
  const double dminvsqdqop1 = xf26*(-xf16*xf30 - xf2*xf29 - xf30*xf9 + 2*xf18/(std::pow(qop1, 3)*xf20));
  const double dminvsqdlam1 = xf26*(2*xf0*xf11*xf21*xf4 - xf16*xf31 - xf31*xf9);
  const double dminvsqdphi1 = xf26*xf28;
  Matrix<double, 1, 6> res;
  res(0,0) = dminvsqdqop0;
  res(0,1) = dminvsqdlam0;
  res(0,2) = dminvsqdphi0;
  res(0,3) = dminvsqdqop1;
  res(0,4) = dminvsqdlam1;
  res(0,5) = dminvsqdphi1;
  
//   std::cout << "massJacobianAlt m = " << m << std::endl;




  
  return res;
}

Matrix<double, 6, 6> ResidualGlobalCorrectionMakerBase::massinvsqHessianAltD(const Matrix<double, 7, 1> &state0, const Matrix<double, 7, 1> &state1, double dmass) const {
  
  const double mp = dmass;
  
  const double qop0 = state0[6]/state0.segment<3>(3).norm();
  const double lam0 = std::atan(state0[5]/std::sqrt(state0[3]*state0[3] + state0[4]*state0[4]));
  const double phi0 = std::atan2(state0[4], state0[3]);
  
  const double qop1 = state1[6]/state1.segment<3>(3).norm();
  const double lam1 = std::atan(state1[5]/std::sqrt(state1[3]*state1[3] + state1[4]*state1[4]));
  const double phi1 = std::atan2(state1[4], state1[3]);

  const double xf0 = std::sin(lam0);
  const double xf1 = std::sin(lam1);
  const double xf2 = xf0*xf1;
  const double xf3 = std::pow(qop0, -3);
  const double xf4 = 1.0/std::fabs(qop1);
  const double xf5 = (((qop0) > 0) - ((qop0) < 0));
  const double xf6 = xf4*xf5;
  const double xf7 = 4*xf3*xf6;
  const double xf8 = std::cos(lam0);
  const double xf9 = xf7*xf8;
  const double xf10 = std::cos(lam1);
  const double xf11 = std::sin(phi0);
  const double xf12 = std::sin(phi1);
  const double xf13 = xf11*xf12;
  const double xf14 = xf10*xf13;
  const double xf15 = std::cos(phi0);
  const double xf16 = std::cos(phi1);
  const double xf17 = xf15*xf16;
  const double xf18 = xf10*xf17;
  const double xf19 = std::pow(mp, 2);
  const double xf20 = std::pow(qop0, -2);
  const double xf21 = xf19 + xf20;
  const double xf22 = std::pow(qop1, -2);
  const double xf23 = xf19 + xf22;
  const double xf24 = std::sqrt(xf23);
  const double xf25 = std::sqrt(xf21);
  const double xf26 = 1.0/xf25;
  const double xf27 = 1.0/std::fabs(qop0);
  const double xf28 = xf27*xf4;
  const double xf29 = xf1*xf28;
  const double xf30 = 2*xf0;
  const double xf31 = xf29*xf30;
  const double xf32 = xf10*xf28;
  const double xf33 = 2*xf8;
  const double xf34 = xf32*xf33;
  const double xf35 = xf13*xf34 + xf17*xf34;
  const double xf36 = xf31 + xf35;
  const double xf37 = 2*xf19 + 2*xf24*xf25 - xf36;
  const double xf38 = 1.0/std::pow(xf37, 2);
  const double xf39 = 4*xf0;
  const double xf40 = xf20*xf6;
  const double xf41 = xf1*xf40;
  const double xf42 = 4*xf8;
  const double xf43 = xf40*xf42;
  const double xf44 = -xf14*xf43 - xf18*xf43 + 4*xf24*xf26*xf3 - xf39*xf41;
  const double xf45 = xf1*xf30;
  const double xf46 = xf33*xf40;
  const double xf47 = 1.0/std::pow(xf37, 3);
  const double xf48 = xf47*(-xf14*xf46 - xf18*xf46 + 2*xf24*xf26*xf3 - xf40*xf45);
  const double xf49 = xf33*xf41;
  const double xf50 = xf10*xf30;
  const double xf51 = xf40*xf50;
  const double xf52 = xf38*(xf13*xf51 + xf17*xf51 - xf49);
  const double xf53 = xf32*xf39;
  const double xf54 = 4*xf1*xf27*xf4*xf8 - xf13*xf53 - xf17*xf53;
  const double xf55 = xf11*xf16;
  const double xf56 = xf10*xf46;
  const double xf57 = xf12*xf15;
  const double xf58 = xf55*xf56 - xf56*xf57;
  const double xf59 = xf38*xf58;
  const double xf60 = -4*xf10*xf12*xf15*xf27*xf4*xf8 + xf32*xf42*xf55;
  const double xf61 = -xf60;
  const double xf62 = (((qop1) > 0) - ((qop1) < 0));
  const double xf63 = xf20*xf22*xf5*xf62;
  const double xf64 = xf33*xf63;
  const double xf65 = 1.0/xf24;
  const double xf66 = std::pow(qop1, -3);
  const double xf67 = xf38*(xf14*xf64 + xf18*xf64 - 2*xf26*xf3*xf65*xf66 + xf45*xf63);
  const double xf68 = xf27*xf62;
  const double xf69 = xf22*xf68;
  const double xf70 = xf1*xf69;
  const double xf71 = xf42*xf69;
  const double xf72 = -xf14*xf71 - xf18*xf71 + 4*xf25*xf65*xf66 - xf39*xf70;
  const double xf73 = xf38*(xf13*xf49 + xf17*xf49 - xf51);
  const double xf74 = xf29*xf8;
  const double xf75 = 4*xf74;
  const double xf76 = 4*xf0*xf10*xf27*xf4 - xf13*xf75 - xf17*xf75;
  const double xf77 = -xf38*xf58;
  const double xf78 = xf30*xf32;
  const double xf79 = xf47*(2*xf1*xf27*xf4*xf8 - xf13*xf78 - xf17*xf78);
  const double xf80 = -xf36*xf38;
  const double xf81 = xf55*xf78 - xf57*xf78;
  const double xf82 = xf38*xf81;
  const double xf83 = xf33*xf70;
  const double xf84 = xf50*xf69;
  const double xf85 = xf38*(xf13*xf84 + xf17*xf84 - xf83);
  const double xf86 = xf38*(xf13*xf31 + xf17*xf31 + xf34);
  const double xf87 = -xf38*xf81;
  const double xf88 = -2*xf10*xf12*xf15*xf27*xf4*xf8 + xf34*xf55;
  const double xf89 = -xf47*xf88;
  const double xf90 = -xf35*xf38;
  const double xf91 = xf33*xf69;
  const double xf92 = xf10*xf91;
  const double xf93 = xf55*xf92 - xf57*xf92;
  const double xf94 = xf38*xf93;
  const double xf95 = 2*xf74;
  const double xf96 = xf55*xf95 - xf57*xf95;
  const double xf97 = xf38*xf96;
  const double xf98 = xf35*xf38;
  const double xf99 = xf47*(-xf14*xf91 - xf18*xf91 + 2*xf25*xf65*xf66 - xf45*xf69);
  const double xf100 = 4*xf66*xf68;
  const double xf101 = xf100*xf8;
  const double xf102 = xf38*(xf13*xf83 + xf17*xf83 - xf84);
  const double xf103 = -xf38*xf93;
  const double xf104 = xf47*(2*xf0*xf10*xf27*xf4 - xf13*xf95 - xf17*xf95);
  const double xf105 = -xf38*xf96;
  const double xf106 = xf47*xf88;
  const double d2minvsqdqop0dqop0 = xf38*(xf14*xf9 + xf18*xf9 + xf2*xf7 - 6*xf24*xf26/std::pow(qop0, 4) + 2*xf24/(std::pow(qop0, 6)*std::pow(xf21, 3.0/2.0))) + xf44*xf48;
  const double d2minvsqdqop0dlam0 = xf48*xf54 + xf52;
  const double d2minvsqdqop0dphi0 = xf48*xf61 + xf59;
  const double d2minvsqdqop0dqop1 = xf48*xf72 + xf67;
  const double d2minvsqdqop0dlam1 = xf48*xf76 + xf73;
  const double d2minvsqdqop0dphi1 = xf48*xf60 + xf77;
  const double d2minvsqdlam0dqop0 = xf44*xf79 + xf52;
  const double d2minvsqdlam0dlam0 = xf54*xf79 + xf80;
  const double d2minvsqdlam0dphi0 = xf61*xf79 + xf82;
  const double d2minvsqdlam0dqop1 = xf72*xf79 + xf85;
  const double d2minvsqdlam0dlam1 = xf76*xf79 + xf86;
  const double d2minvsqdlam0dphi1 = xf60*xf79 + xf87;
  const double d2minvsqdphi0dqop0 = xf44*xf89 + xf59;
  const double d2minvsqdphi0dlam0 = xf54*xf89 + xf82;
  const double d2minvsqdphi0dphi0 = xf61*xf89 + xf90;
  const double d2minvsqdphi0dqop1 = xf72*xf89 + xf94;
  const double d2minvsqdphi0dlam1 = xf76*xf89 + xf97;
  const double d2minvsqdphi0dphi1 = xf60*xf89 + xf98;
  const double d2minvsqdqop1dqop0 = xf44*xf99 + xf67;
  const double d2minvsqdqop1dlam0 = xf54*xf99 + xf85;
  const double d2minvsqdqop1dphi0 = xf61*xf99 + xf94;
  const double d2minvsqdqop1dqop1 = xf38*(xf100*xf2 + xf101*xf14 + xf101*xf18 - 6*xf25*xf65/std::pow(qop1, 4) + 2*xf25/(std::pow(qop1, 6)*std::pow(xf23, 3.0/2.0))) + xf72*xf99;
  const double d2minvsqdqop1dlam1 = xf102 + xf76*xf99;
  const double d2minvsqdqop1dphi1 = xf103 + xf60*xf99;
  const double d2minvsqdlam1dqop0 = xf104*xf44 + xf73;
  const double d2minvsqdlam1dlam0 = xf104*xf54 + xf86;
  const double d2minvsqdlam1dphi0 = xf104*xf61 + xf97;
  const double d2minvsqdlam1dqop1 = xf102 + xf104*xf72;
  const double d2minvsqdlam1dlam1 = xf104*xf76 + xf80;
  const double d2minvsqdlam1dphi1 = xf104*xf60 + xf105;
  const double d2minvsqdphi1dqop0 = xf106*xf44 + xf77;
  const double d2minvsqdphi1dlam0 = xf106*xf54 + xf87;
  const double d2minvsqdphi1dphi0 = xf106*xf61 + xf98;
  const double d2minvsqdphi1dqop1 = xf103 + xf106*xf72;
  const double d2minvsqdphi1dlam1 = xf105 + xf106*xf76;
  const double d2minvsqdphi1dphi1 = xf106*xf60 + xf90;
  Matrix<double, 6, 6> res;
  res(0,0) = d2minvsqdqop0dqop0;
  res(0,1) = d2minvsqdqop0dlam0;
  res(0,2) = d2minvsqdqop0dphi0;
  res(0,3) = d2minvsqdqop0dqop1;
  res(0,4) = d2minvsqdqop0dlam1;
  res(0,5) = d2minvsqdqop0dphi1;
  res(1,0) = d2minvsqdlam0dqop0;
  res(1,1) = d2minvsqdlam0dlam0;
  res(1,2) = d2minvsqdlam0dphi0;
  res(1,3) = d2minvsqdlam0dqop1;
  res(1,4) = d2minvsqdlam0dlam1;
  res(1,5) = d2minvsqdlam0dphi1;
  res(2,0) = d2minvsqdphi0dqop0;
  res(2,1) = d2minvsqdphi0dlam0;
  res(2,2) = d2minvsqdphi0dphi0;
  res(2,3) = d2minvsqdphi0dqop1;
  res(2,4) = d2minvsqdphi0dlam1;
  res(2,5) = d2minvsqdphi0dphi1;
  res(3,0) = d2minvsqdqop1dqop0;
  res(3,1) = d2minvsqdqop1dlam0;
  res(3,2) = d2minvsqdqop1dphi0;
  res(3,3) = d2minvsqdqop1dqop1;
  res(3,4) = d2minvsqdqop1dlam1;
  res(3,5) = d2minvsqdqop1dphi1;
  res(4,0) = d2minvsqdlam1dqop0;
  res(4,1) = d2minvsqdlam1dlam0;
  res(4,2) = d2minvsqdlam1dphi0;
  res(4,3) = d2minvsqdlam1dqop1;
  res(4,4) = d2minvsqdlam1dlam1;
  res(4,5) = d2minvsqdlam1dphi1;
  res(5,0) = d2minvsqdphi1dqop0;
  res(5,1) = d2minvsqdphi1dlam0;
  res(5,2) = d2minvsqdphi1dphi0;
  res(5,3) = d2minvsqdphi1dqop1;
  res(5,4) = d2minvsqdphi1dlam1;
  res(5,5) = d2minvsqdphi1dphi1;
  
  return res;
}

//define this as a plug-in
// DEFINE_FWK_MODULE(ResidualGlobalCorrectionMakerBase);
