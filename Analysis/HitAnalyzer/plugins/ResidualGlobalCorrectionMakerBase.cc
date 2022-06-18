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
    genEventInfoToken_ = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
    pileupSummaryToken_ = consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("addPileupInfo"));
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
  const double x1 = std::pow(x0, -3.0/2.0);
  const double x2 = Iy1*W0x;
  const double x3 = Ix1*W0y;
  const double x4 = std::pow(qop0, 2);
  const double x5 = std::pow(W0x, 2) + std::pow(W0y, 2);
  const double x6 = std::pow(x5, -1.0/2.0);
  const double x7 = std::sin(lam0);
  const double x8 = Iz1*x7;
  const double x9 = std::cos(lam0);
  const double x10 = std::cos(phi0);
  const double x11 = Ix1*x10;
  const double x12 = std::sin(phi0);
  const double x13 = Iy1*x12;
  const double x14 = x11*x9 + x13*x9 + x8;
  const double x15 = x6/x14;
  const double x16 = dEdx*q*x1*x15*x4*std::sqrt(std::pow(mass, 2)*x4 + x0);
  const double x17 = Ix1*W0x*W0z + Iy1*W0y*W0z - Iz1*x5;
  const double x18 = std::pow(x14, -2);
  const double x19 = Iz1*Jx1;
  const double x20 = Iz1*Jy1;
  const double x21 = lam0 + phi0;
  const double x22 = -phi0;
  const double x23 = lam0 + x22;
  const double x24 = std::pow(Ix1*std::cos(x21) + Ix1*std::cos(x23) + Iy1*std::sin(x21) - Iy1*std::sin(x23) + 2*x8, -2);
  const double x25 = 2*Ix1;
  const double x26 = Jy1*x25;
  const double x27 = 2*Iy1;
  const double x28 = Jx1*x27;
  const double x29 = 2*lam0;
  const double x30 = std::cos(x29);
  const double x31 = phi0 + x29;
  const double x32 = std::cos(x31);
  const double x33 = std::sin(x31);
  const double x34 = Ix1*Jz1;
  const double x35 = Iy1*Jz1;
  const double x36 = x22 + x29;
  const double x37 = std::cos(x36);
  const double x38 = std::sin(x36);
  const double x39 = -x2 + x3;
  const double x40 = hz*x9;
  const double x41 = hx*x7 - x10*x40;
  const double x42 = hy*x7 - x12*x40;
  const double x43 = x9*(hx*x12 - hy*x10);
  const double x44 = x10*x9;
  const double x45 = x12*x9;
  const double x46 = Jx1*x44 + Jy1*x45 + Jz1*x7;
  const double x47 = -Ix1*x42 + Iy1*x41 - Iz1*x43;
  const double x48 = B*qop0*x6/std::pow(x14, 3);
  const double x49 = x48*(x14*(-Jx1*x42 + Jy1*x41 - Jz1*x43) - x46*x47);
  const double x50 = Iz1*Kx1;
  const double x51 = Iz1*Ky1;
  const double x52 = Ky1*x25;
  const double x53 = Kx1*x27;
  const double x54 = Ix1*Kz1;
  const double x55 = Iy1*Kz1;
  const double x56 = Kx1*x44 + Ky1*x45 + Kz1*x7;
  const double x57 = x48*(x14*(-Kx1*x42 + Ky1*x41 - Kz1*x43) - x47*x56);
  const double x58 = x14*x5;
  const double x59 = W0z*x14;
  const double dqopdqop0 = std::pow(q, 3)*x1*std::fabs(qop0)/qop0;
  const double dqopdlam0 = 0;
  const double dqopdphi0 = 0;
  const double dqopdxt0 = x16*(x2 - x3);
  const double dqopdyt0 = -x16*x17;
  const double ddxdzdqop0 = 0;
  const double ddxdzdlam0 = x18*(Jz1*x11 + Jz1*x13 - x10*x19 - x12*x20);
  const double ddxdzdphi0 = x24*(x19*x32 - x19*x37 + x20*x33 + x20*x38 + x26*x30 + x26 - x28*x30 - x28 - x32*x34 - x33*x35 + x34*x37 - x35*x38);
  const double ddxdzdxt0 = x39*x49;
  const double ddxdzdyt0 = x17*x49;
  const double ddydzdqop0 = 0;
  const double ddydzdlam0 = x18*(Kz1*x11 + Kz1*x13 - x10*x50 - x12*x51);
  const double ddydzdphi0 = x24*(x30*x52 - x30*x53 + x32*x50 - x32*x54 + x33*x51 - x33*x55 - x37*x50 + x37*x54 + x38*x51 - x38*x55 + x52 - x53);
  const double ddydzdxt0 = x39*x57;
  const double ddydzdyt0 = x17*x57;
  const double dxdqop0 = 0;
  const double dxdlam0 = 0;
  const double dxdphi0 = 0;
  const double dxdxt0 = x15*(x14*(-Jx1*W0y + Jy1*W0x) + x39*x46);
  const double dxdyt0 = x15*(Jz1*x58 + x17*x46 - x59*(Jx1*W0x + Jy1*W0y));
  const double dydqop0 = 0;
  const double dydlam0 = 0;
  const double dydphi0 = 0;
  const double dydxt0 = x15*(x14*(-Kx1*W0y + Ky1*W0x) + x39*x56);
  const double dydyt0 = x15*(Kz1*x58 + x17*x56 - x59*(Kx1*W0x + Ky1*W0y));
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
