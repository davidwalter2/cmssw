// -*- C++ -*-
//
// Package:    TrackAnalysis/HitAnalyzer
// Class:      HitAnalyzer
//
/**\class HitAnalyzer HitAnalyzer.cc TrackAnalysis/HitAnalyzer/plugins/HitAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Michail Bachtis
//         Created:  Mon, 21 Mar 2016 14:17:37 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

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


#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "functions.h"
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <iostream>

using namespace Eigen;

constexpr unsigned int max_n = 25; //!< In order to avoid use of dynamic memory

typedef Matrix<double, Dynamic, Dynamic, 0, max_n, max_n> MatrixNd;
typedef Array<double, Dynamic, Dynamic, 0, max_n, max_n> ArrayNd;
typedef Matrix<double, Dynamic, Dynamic, 0, 2 * max_n, 2 * max_n> Matrix2Nd;
typedef Matrix<double, Dynamic, Dynamic, 0, 3 * max_n, 3 * max_n> Matrix3Nd;
typedef Matrix<double, 2, Dynamic, 0, 2, max_n> Matrix2xNd;
typedef Array<double, 2, Dynamic, 0, 2, max_n> Array2xNd;
typedef Matrix<double, 3, Dynamic, 0, 3, max_n> Matrix3xNd;
typedef Matrix<double, Dynamic, 3, 0, max_n, 3> MatrixNx3d;
typedef Matrix<double, Dynamic, 5, 0, max_n, 5> MatrixNx5d;
typedef Matrix<double, Dynamic, 1, 0, max_n, 1> VectorNd;
typedef Matrix<double, Dynamic, 1, 0, 2 * max_n, 1> Vector2Nd;
typedef Matrix<double, Dynamic, 1, 0, 3 * max_n, 1> Vector3Nd;
typedef Matrix<double, 1, Dynamic, 1, 1, max_n> RowVectorNd;
typedef Matrix<double, 1, Dynamic, 1, 1, 2 * max_n> RowVector2Nd;
typedef Matrix<double, 2, 2> Matrix2d;
typedef Matrix<double, 5, 5> Matrix5d;
typedef Matrix<double, 6, 6> Matrix6d;
typedef Matrix<double, 5, 6> Matrix56d;
typedef Matrix<double, 5, 1> Vector5d;
typedef Matrix<double, 6, 1> Vector6d;

typedef ROOT::Math::SMatrix<double, 2> SMatrix22;

//
// class declaration
//

class HitAnalyzer : public edm::EDAnalyzer
{
public:
  explicit HitAnalyzer(const edm::ParameterSet &);
  ~HitAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event &, const edm::EventSetup &) override;
  virtual void endJob() override;

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
  //  edm::EDGetTokenT<reco::TrackCollection>      inputTracks_;
  edm::EDGetTokenT<std::vector<Trajectory>> inputTraj_;
  edm::EDGetTokenT<std::vector<reco::GenParticle>> GenParticlesToken_;

  TFile *fout;
  TTree *tree;

  const int N = 25;

  int n;
  std::vector<float> z;
  std::vector<float> eta;
  std::vector<float> phi;
  std::vector<float> r;

  std::vector<float> xUnc;
  std::vector<float> yUnc;
  std::vector<float> zUnc;
  std::vector<float> etaUnc;
  std::vector<float> phiUnc;
  std::vector<float> rUnc;

  std::vector<float> pt;

  std::vector<int> detector;
  std::vector<int> stereo;
  std::vector<int> glued;
  std::vector<int> layer;

  //local positions and errors
  std::vector<float> localx;
  std::vector<float> localy;
  std::vector<float> localz;
  std::vector<float> localxErr;
  std::vector<float> localyErr;
  std::vector<float> localxyErr;

  std::vector<float> globalrErr;
  std::vector<float> globalzErr;
  std::vector<float> globalrphiErr;
  std::vector<float> localx_state;
  std::vector<float> localy_state;

  //material stuff
  std::vector<std::vector<double>> trackQf;
  std::vector<std::vector<double>> trackHf;
  std::vector<std::vector<double>> trackFf;
  std::vector<std::vector<double>> trackQr;
  std::vector<std::vector<double>> trackHr;
  std::vector<std::vector<double>> trackFr;
  std::vector<std::vector<double>> trackC;

  std::vector<std::vector<double>> updState;
  std::vector<std::vector<double>> backPropState;
  std::vector<std::vector<double>> forwardPropState;

  std::vector<std::vector<double>> updStateLocal;
  std::vector<std::vector<double>> fwdPredStateLocal;
  std::vector<std::vector<double>> bkgPredStateLocal;
  
  std::vector<std::vector<double>> curvpars;
  std::vector<int> hitDimension;

  float trackEta;
  float trackPhi;
  float trackPt;
  float trackPtErr;
  float trackZ0;
  float trackX0;
  float trackY0;
  float trackCharge;

  float genPt;
  float genCharge;
  int ninvalidHits;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HitAnalyzer::HitAnalyzer(const edm::ParameterSet &iConfig)

{
  //now do what ever initialization is needed
  //  inputTracks_ = consumes<reco::TrackCollection>(edm::InputTag("TrackRefitter"));
  inputTraj_ = consumes<std::vector<Trajectory>>(edm::InputTag("TrackRefitter"));
  GenParticlesToken_ = consumes<std::vector<reco::GenParticle>>(edm::InputTag("genParticles"));

  n = 0;

  fout = new TFile("trackTree.root", "RECREATE");
  tree = new TTree("tree", "tree");

  tree->Branch("n", &n, "n/I");
  tree->Branch("ninvalidHits", &ninvalidHits, "ninvalidHits/I");
  tree->Branch("z", &z);
  tree->Branch("eta", &eta);
  tree->Branch("phi", &phi);
  tree->Branch("r", &r);
  tree->Branch("pt", &pt);

  tree->Branch("xUnc", &xUnc);
  tree->Branch("yUnc", &yUnc);
  tree->Branch("zUnc", &zUnc);
  tree->Branch("etaUnc", &etaUnc);
  tree->Branch("phiUnc", &phiUnc);
  tree->Branch("rUnc", &rUnc);

  tree->Branch("stereo", &stereo);
  tree->Branch("glued", &glued);
  tree->Branch("detector", &detector);
  tree->Branch("layer", &layer);

  tree->Branch("trackPt", &trackPt, "trackPt/F");
  tree->Branch("trackPtErr", &trackPtErr, "trackPtErr/F");
  tree->Branch("trackEta", &trackEta, "trackEta/F");
  tree->Branch("trackPhi", &trackPhi, "trackPhi/F");
  tree->Branch("trackX0", &trackX0, "trackX0/F");
  tree->Branch("trackY0", &trackY0, "trackY0/F");
  tree->Branch("trackZ0", &trackZ0, "trackZ0/F");
  tree->Branch("trackCharge", &trackCharge, "trackCharge/F");

  tree->Branch("genPt", &genPt, "genPt/F");
  tree->Branch("genCharge", &genCharge, "genCharge/F");

  tree->Branch("localx", &localx);
  tree->Branch("localy", &localy);
  tree->Branch("localz", &localz);
  tree->Branch("localx_state", &localx_state);
  tree->Branch("localy_state", &localy_state);
  tree->Branch("localxErr", &localxErr);
  tree->Branch("localyErr", &localyErr);
  tree->Branch("localxyErr", &localxyErr);
  tree->Branch("globalrErr", &globalrErr);
  tree->Branch("globalzErr", &globalzErr);
  tree->Branch("globalrphiErr", &globalrphiErr);
  tree->Branch("trackQf", &trackQf);
  tree->Branch("trackHf", &trackHf);
  tree->Branch("trackFf", &trackFf);
  tree->Branch("trackQr", &trackQr);
  tree->Branch("trackHr", &trackHr);
  tree->Branch("trackFr", &trackFr);
  tree->Branch("trackC", &trackC);
  tree->Branch("updState", &updState);
  tree->Branch("backPropState", &backPropState);
  tree->Branch("forwardPropState", &forwardPropState);
  tree->Branch("updStateLocal", &updStateLocal);
  tree->Branch("fwdPredStateLocal", &fwdPredStateLocal);
  tree->Branch("bkgPredStateLocal", &bkgPredStateLocal);
  tree->Branch("curvpars", &curvpars);
  tree->Branch("hitDimension", &hitDimension);
}

HitAnalyzer::~HitAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void HitAnalyzer::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  using namespace edm;

  Handle<std::vector<reco::GenParticle>> genPartCollection;
  iEvent.getByToken(GenParticlesToken_, genPartCollection);

  auto genParticles = *genPartCollection.product();

  // loop over gen particles

  edm::ESHandle<GlobalTrackingGeometry> globalGeometry;
  iSetup.get<GlobalTrackingGeometryRecord>().get(globalGeometry);

  Handle<std::vector<Trajectory>> trajH;
  iEvent.getByToken(inputTraj_, trajH);

  ESHandle<MagneticField> magfield;
  iSetup.get<IdealMagneticFieldRecord>().get(magfield);
  auto field = magfield.product();
  
  edm::ESHandle<TransientTrackingRecHitBuilder> ttrh;
  iSetup.get<TransientRecHitRecord>().get("WithAngleAndTemplate",ttrh);

  for (unsigned int j = 0; j < trajH->size(); ++j)
  {

    //     const reco::Track& track = (*trackH)[i];
    //     if (track.lost()>0)
    //       continue;

    const std::vector<TrajectoryMeasurement> &tms = (*trajH)[j].measurements();

    ////

    if (((*trajH)[j].direction()) == alongMomentum)
    {

      TrajectoryStateOnSurface measurement = (*trajH)[j].firstMeasurement().updatedState();

      trackPt = measurement.globalMomentum().perp();
      //FIX BUG
      trackPtErr = sqrt(measurement.curvilinearError().matrix()[0][0]) * trackPt;

      trackEta = measurement.globalMomentum().eta();
      trackPhi = measurement.globalMomentum().phi();
      trackX0 = measurement.globalPosition().x();
      trackY0 = measurement.globalPosition().y();
      trackZ0 = measurement.globalPosition().z();
      trackCharge = measurement.charge();
    }
    else
    {

      TrajectoryStateOnSurface measurement = (*trajH)[j].lastMeasurement().updatedState();

      trackPt = measurement.globalMomentum().perp();
      trackPtErr = sqrt(measurement.curvilinearError().matrix()[0][0]) * trackPt;
      trackEta = measurement.globalMomentum().eta();
      trackPhi = measurement.globalMomentum().phi();
      trackX0 = measurement.globalPosition().x();
      trackY0 = measurement.globalPosition().y();
      trackZ0 = measurement.globalPosition().z();
      trackCharge = measurement.charge();
    }

    //     printf("First point %f %f %f  - %f %f %f\n",trackX0,trackY0,trackZ0,trackPt,trackEta,trackPhi);

    for (std::vector<reco::GenParticle>::const_iterator g = genParticles.begin(); g != genParticles.end(); ++g)
    {

      float dR = deltaR(g->phi(), trackPhi, g->eta(), trackEta);

      if (dR < 0.15)
      {
        genPt = g->pt();
        genCharge = g->charge();
      }
      else
        continue;
    }

    ////
    n = 0;
    ninvalidHits = 0;

    z.clear();
    eta.clear();
    phi.clear();
    r.clear();

    xUnc.clear();
    yUnc.clear();
    zUnc.clear();
    etaUnc.clear();
    phiUnc.clear();
    rUnc.clear();

    pt.clear();

    detector.clear();
    stereo.clear();
    glued.clear();
    layer.clear();

    //local positions and errors
    localx.clear();
    localy.clear();
    localz.clear();
    localxErr.clear();
    localyErr.clear();
    localxyErr.clear();

    localx_state.clear();
    localy_state.clear();

    globalrErr.clear();
    globalzErr.clear();
    globalrphiErr.clear();

    updState.clear();
    backPropState.clear();
    forwardPropState.clear();

    updStateLocal.clear();
    fwdPredStateLocal.clear();
    bkgPredStateLocal.clear();
    
    //material stuff
    trackQf.clear();
    trackHf.clear();
    trackFf.clear();
    trackQr.clear();
    trackHr.clear();
    trackFr.clear();
    trackC.clear();
    
    hitDimension.clear();

    const float mass = 0.1395703;
    const float maxDPhi = 1.6;
//     printf("traj propdir = %i\n", int((*trajH)[j].direction()));
    

    PropagationDirection rpropdir = (*trajH)[j].direction();
    PropagationDirection fpropdir = rpropdir == alongMomentum ? oppositeToMomentum : alongMomentum;
    
    PropagatorWithMaterial rPropagator(rpropdir, mass, field, maxDPhi, true, -1., false);
    PropagatorWithMaterial fPropagator(fpropdir, mass, field, maxDPhi, true, -1., false);
    
    KFSwitching1DUpdator updator;
    
    for (unsigned int i = 0; i < tms.size(); ++i)
    {
      if (!tms[i].recHit()->isValid()) {
        ++ninvalidHits;
      }
    }
    if (ninvalidHits>0) {
      continue;
    }
    
    for (unsigned int i = 0; i < tms.size(); ++i)
    {
      TrajectoryStateOnSurface updatedState = tms[i].updatedState();

//       if (!tms[i].recHit()->isValid()){
//         ninvalidHits++;
//         continue;
//       }
      
      if (!updatedState.isValid())
        continue;

      pt.push_back(updatedState.globalMomentum().perp());

      const GeomDet *detectorG = globalGeometry->idToDet(tms[i].recHit()->geographicalId());
      LocalPoint local(tms[i].recHit()->localPosition().x(), tms[i].recHit()->localPosition().y(), tms[i].recHit()->localPosition().z());

      if (detectorG->subDetector() == GeomDetEnumerators::PixelBarrel)
      {
        PXBDetId detid(tms[i].recHit()->rawId());
        layer.push_back(detid.layer());
        stereo.push_back(0);
        glued.push_back(0);
      }
      else if (detectorG->subDetector() == GeomDetEnumerators::PixelEndcap)
      {
        PXFDetId detid(tms[i].recHit()->rawId());
        layer.push_back(-1 * (detid.side() == 1) * detid.disk() + (detid.side() == 2) * detid.disk());
        stereo.push_back(0);
        glued.push_back(0);
      }
      else if (detectorG->subDetector() == GeomDetEnumerators::TIB)
      {
        SiStripDetId detid(tms[i].recHit()->rawId());
//         TIBDetId detid(tms[i].recHit()->rawId());
        layer.push_back(0);
        if (detid.stereo() != 0)
          stereo.push_back(1);
        if (detid.glued() != 0)
          glued.push_back(1);
      }
      else if (detectorG->subDetector() == GeomDetEnumerators::TOB)
      {
        SiStripDetId detid(tms[i].recHit()->rawId());
//         TOBDetId detid(tms[i].recHit()->rawId());
//         layer.push_back(detid.layer());
        layer.push_back(0);
        if (detid.stereo() != 0)
          stereo.push_back(1);
        if (detid.glued() != 0)
          glued.push_back(1);
      }
      else if (detectorG->subDetector() == GeomDetEnumerators::TID)
      {
        SiStripDetId detid(tms[i].recHit()->rawId());
//         TIDDetId detid(tms[i].recHit()->rawId());
//         layer.push_back(-1 * (detid.side() == 1) * detid.wheel() + (detid.side() == 2) * detid.wheel());
        layer.push_back(0);
        if (detid.stereo() != 0)
          stereo.push_back(1);
        if (detid.glued() != 0)
          glued.push_back(1);
      }
      else if (detectorG->subDetector() == GeomDetEnumerators::TEC)
      {
        SiStripDetId detid(tms[i].recHit()->rawId());
//         TECDetId detid(tms[i].recHit()->rawId());
//         layer.push_back(-1 * (detid.side() == 1) * detid.wheel() + (detid.side() == 2) * detid.wheel());
        layer.push_back(0);
        if (detid.stereo() != 0)
          stereo.push_back(1);
        if (detid.glued() != 0)
          glued.push_back(1);
      }

      // material info

      // Get surface
      const Surface &surface = updatedState.surface();


      TrajectoryStateOnSurface const& fwdtsos = tms[i].forwardPredictedState();
      TrajectoryStateOnSurface const& revtsos = tms[i].backwardPredictedState();
//       TrajectoryStateOnSurface const& revtsos = tms[i].updatedState();
      auto const& hit = tms[i].recHit();
      
//       TrajectoryStateCombiner combiner;
//       revtsos = combiner(revtsos,fwdtsos);
      
      
      JacobianCurvilinearToLocal hf(fwdtsos.surface(), fwdtsos.localParameters(), *fwdtsos.magneticField());
      const AlgebraicMatrix55 &jachf = hf.jacobian();
      std::vector<double> Hf(jachf.Array(), jachf.Array()+25);
      trackHf.push_back(Hf);
      
      JacobianCurvilinearToLocal hr(revtsos.surface(), revtsos.localParameters(), *revtsos.magneticField());
      const AlgebraicMatrix55 &jachr = hr.jacobian();
      std::vector<double> Hr(jachr.Array(), jachr.Array()+25);
      trackHr.push_back(Hr);

      const AlgebraicSymMatrix55 &covariance = updatedState.curvilinearError().matrix();
      AlgebraicMatrix55 covfull(covariance);
      std::vector<double> C(covfull.Array(), covfull.Array()+25);
      trackC.push_back(C);

// //       LocalTrajectoryParameters updPars = updatedState.localParameters();
//       CurvilinearTrajectoryParameters curvPrs(updatedState.globalPosition(), updatedState.globalMomentum(), updatedState.charge());
// //       AlgebraicVector5 curvPrsv = updPars.vector();
// 
//       Vector5d curvPrsv_copy;
//       curvPrsv_copy << curvPrsv[0], curvPrsv[1], curvPrsv[2], curvPrsv[3], curvPrsv[4];
//       std::vector<double> curvPrsv_vec;
//       curvPrsv_vec.resize(5);
//       Eigen::Map<Vector5d>(curvPrsv_vec.data(), 5, 1) = curvPrsv_copy;
//       updState.push_back(curvPrsv_vec);
      
      CurvilinearTrajectoryParameters curvPrs(updatedState.globalPosition(), updatedState.globalMomentum(), updatedState.charge());
      AlgebraicVector5 ustatevec = curvPrs.vector();
      std::vector<double> ustate(ustatevec.Array(), ustatevec.Array()+5);
      updState.push_back(ustate);
      
      CurvilinearTrajectoryParameters fcurvPrs(fwdtsos.globalPosition(), fwdtsos.globalMomentum(), fwdtsos.charge());
      AlgebraicVector5 fstatevec = fcurvPrs.vector();
      std::vector<double> fstate(fstatevec.Array(), fstatevec.Array()+5);
      forwardPropState.push_back(fstate);
      
      CurvilinearTrajectoryParameters rcurvPrs(revtsos.globalPosition(), revtsos.globalMomentum(), revtsos.charge());
      AlgebraicVector5 rstatevec = rcurvPrs.vector();
      std::vector<double> rstate(rstatevec.Array(), rstatevec.Array()+5);
      backPropState.push_back(rstate);
      
      AlgebraicVector5 ustateveclocal = updatedState.localParameters().vector();
      std::vector<double> ustatelocal(ustateveclocal.Array(), ustateveclocal.Array()+5);
      updStateLocal.push_back(ustatelocal);

      AlgebraicVector5 fstateveclocal = fwdtsos.localParameters().vector();
      std::vector<double> fstatelocal(fstateveclocal.Array(), fstateveclocal.Array()+5);
      fwdPredStateLocal.push_back(fstatelocal);
      
      AlgebraicVector5 bstateveclocal = revtsos.localParameters().vector();
      std::vector<double> bstatelocal(bstateveclocal.Array(), bstateveclocal.Array()+5);
      bkgPredStateLocal.push_back(bstatelocal);

      if (i== (tms.size()-1)) {
        std::vector<double> Ff(25,0.);
        trackFf.push_back(Ff);
        
        std::vector<double> Q(9,0.);
        trackQf.push_back(Q);
      }
      else {
        //compute material effects
        SurfaceSideDefinition::SurfaceSide side = fpropdir == alongMomentum ? SurfaceSideDefinition::afterSurface : SurfaceSideDefinition::beforeSurface;
        TrajectoryStateOnSurface tmptsos(fwdtsos.localParameters(),
                                         LocalTrajectoryError(),
                                         fwdtsos.surface(),
                                         fwdtsos.magneticField(),
                                         side);
        rPropagator.materialEffectsUpdator().updateStateInPlace(tmptsos,rpropdir);
        AlgebraicSymMatrix55 const &qmat = tmptsos.localError().matrix();
        std::vector<double> Q(9,0.);
        for (unsigned int ierr=0, ijerr=0; ierr<3; ++ierr) {
          for (unsigned int jerr=0; jerr<3; ++jerr, ++ijerr) {
            Q[ijerr] = qmat[ierr][jerr];
          }
        }
        trackQf.push_back(Q);
        
        //compute transport jacobian
        //recompute updated state on adjacent measurement
        TrajectoryStateOnSurface const& adjfwdtsos = tms[i+1].forwardPredictedState();
        auto preciseHit = static_cast<TkTransientTrackingRecHitBuilder const*>(ttrh.product())->cloner().makeShared(tms[i+1].recHit(), adjfwdtsos);  //TODO this should be updated with the 
//         auto preciseHit = tms[i+1].recHit();
        TrajectoryStateOnSurface const& updstate = updator.update(adjfwdtsos, *preciseHit);
        
        const auto &propresult = fPropagator.propagateWithPath(updstate, fwdtsos.surface());
        AnalyticalCurvilinearJacobian curvjac;
        GlobalVector h = updstate.globalParameters().magneticFieldInInverseGeV(updstate.globalParameters().position());
        curvjac.computeFullJacobian(updstate.globalParameters(), propresult.first.globalParameters().position(), propresult.first.globalParameters().momentum(), h, propresult.second);
        const AlgebraicMatrix55 &jacF = curvjac.jacobian();
        
        double eloss = fwdtsos.localParameters().signedInverseMomentum()/tmptsos.localParameters().signedInverseMomentum();
        AlgebraicMatrix55 elossmat;
        elossmat(0,0) = eloss;
        elossmat(1,1) = 1.;
        elossmat(2,2) = 1.;
        elossmat(3,3) = 1.;
        elossmat(4,4) = 1.;
        
        AlgebraicMatrix55 Fmat = fpropdir==alongMomentum ? AlgebraicMatrix55(eloss*jacF) : AlgebraicMatrix55(jacF*eloss);
        
        std::vector<double> F(Fmat.Array(), Fmat.Array()+25);
        trackFf.push_back(F);                              
        
      }
      
      if (i==0) {
        std::vector<double> Ff(25,0.);
        trackFr.push_back(Ff);
        
        std::vector<double> Q(9,0.);
        trackQr.push_back(Q);
      }
      else {
        //compute material effects
        SurfaceSideDefinition::SurfaceSide side = rpropdir == alongMomentum ? SurfaceSideDefinition::afterSurface : SurfaceSideDefinition::beforeSurface;
        TrajectoryStateOnSurface tmptsos(revtsos.localParameters(),
                                         LocalTrajectoryError(),
                                         revtsos.surface(),
                                         revtsos.magneticField(),
                                         side);
        fPropagator.materialEffectsUpdator().updateStateInPlace(tmptsos,fpropdir);
        AlgebraicSymMatrix55 const &qmat = tmptsos.localError().matrix();
        std::vector<double> Q(9,0.);
        for (unsigned int ierr=0, ijerr=0; ierr<3; ++ierr) {
          for (unsigned int jerr=0; jerr<3; ++jerr, ++ijerr) {
            Q[ijerr] = qmat[ierr][jerr];
          }
        }
        trackQr.push_back(Q);
        
        //compute transport jacobian
        //recompute updated state on adjacent measurement
        TrajectoryStateOnSurface const& adjrevtsos = tms[i-1].backwardPredictedState();
        auto preciseHit = tms[i-1].recHit(); //n.b. this is exactly the rechit used for the update and propagation during the smoothing
        
        TrajectoryStateOnSurface const& updstate = updator.update(adjrevtsos, *preciseHit);
        
        const auto &propresult = rPropagator.propagateWithPath(updstate, revtsos.surface());
        AnalyticalCurvilinearJacobian curvjac;
        GlobalVector h = updstate.globalParameters().magneticFieldInInverseGeV(updstate.globalParameters().position());
        curvjac.computeFullJacobian(updstate.globalParameters(), propresult.first.globalParameters().position(), propresult.first.globalParameters().momentum(), h, propresult.second);
        const AlgebraicMatrix55 &jacF = curvjac.jacobian();
        
        double eloss = revtsos.localParameters().signedInverseMomentum()/tmptsos.localParameters().signedInverseMomentum();
        AlgebraicMatrix55 elossmat;
        elossmat(0,0) = eloss;
        elossmat(1,1) = 1.;
        elossmat(2,2) = 1.;
        elossmat(3,3) = 1.;
        elossmat(4,4) = 1.;
        
        AlgebraicMatrix55 Fmat = rpropdir==alongMomentum ? AlgebraicMatrix55(eloss*jacF) : AlgebraicMatrix55(jacF*eloss);
        
        std::vector<double> F(Fmat.Array(), Fmat.Array()+25);
        trackFr.push_back(F);                              
        
      }

      detector.push_back(detectorG->subDetector());
      GlobalPoint corrected = detectorG->toGlobal(local);
      z.push_back(corrected.z());
      eta.push_back(corrected.eta());
      phi.push_back(corrected.phi());
      r.push_back(corrected.perp());

      float xhyb = r[n] * TMath::Cos(phi[n]);
      float yhyb = r[n] * TMath::Sin(phi[n]);

      LocalVector perp(0., 0., 1.);
      GlobalVector globperp = detectorG->toGlobal(perp);
      GlobalPoint globpoint(xhyb, yhyb, z[n]);

      CurvilinearTrajectoryParameters cp(globpoint, globperp, trackCharge);

      AlgebraicVector5 cpvec =cp.vector();

      Vector5d vec_curv;
      vec_curv << cpvec[0], cpvec[1], cpvec[2], cpvec[3], cpvec[4];
      std::vector<double> curv;
      
      curv.resize(5);
      Eigen::Map<Vector5d>(curv.data(), 5, 1) = vec_curv;

      curvpars.push_back(curv);

      xUnc.push_back(tms[i].recHit()->globalPosition().x());
      yUnc.push_back(tms[i].recHit()->globalPosition().y());
      zUnc.push_back(tms[i].recHit()->globalPosition().z());
      etaUnc.push_back(tms[i].recHit()->globalPosition().eta());
      phiUnc.push_back(tms[i].recHit()->globalPosition().phi());
      rUnc.push_back(tms[i].recHit()->globalPosition().perp());

      localx.push_back(tms[i].recHit()->localPosition().x());
      localy.push_back(tms[i].recHit()->localPosition().y());
      localz.push_back(tms[i].recHit()->localPosition().z());

      localx_state.push_back(updatedState.localPosition().x());
      localy_state.push_back(updatedState.localPosition().y());

      localxErr.push_back(tms[i].recHit()->localPositionError().xx());
      localyErr.push_back(tms[i].recHit()->localPositionError().yy());
      localxyErr.push_back(tms[i].recHit()->localPositionError().xy());

      // globalrErr.push_back(tms[i].recHit()->errorGlobalR();
      // globalzErr.push_back(tms[i].recHit()->errorGlobalZ();
      // globalrphiErr.push_back(tms[i].recHit()->errorGlobalRPhi();
      GlobalError globalError = ErrorFrameTransformer::transform(tms[i].recHit()->localPositionError(), surface);
      globalrphiErr.push_back(r[n] * std::sqrt(float(globalError.phierr(corrected))));
      globalrErr.push_back(std::sqrt(float(globalError.rerr(corrected))));
      globalzErr.push_back(std::sqrt(float(globalError.czz())));

      hitDimension.push_back(tms[i].recHit()->dimension());
      
      n = n + 1;
    }

    tree->Fill();
  }
}

// ------------ method called once each job just before starting event loop  ------------
void HitAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void HitAnalyzer::endJob()
{
  fout->cd();
  fout->Write();
  fout->Close();
}

// ------------ method called when starting to processes a run  ------------
/*
void 
HitAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
HitAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
HitAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
HitAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HitAnalyzer::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HitAnalyzer);
