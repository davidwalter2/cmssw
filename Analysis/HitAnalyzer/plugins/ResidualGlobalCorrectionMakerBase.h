#ifndef HitAnalyzer_ResidualGlobalCorrectionMakerBase_h
#define HitAnalyzer_ResidualGlobalCorrectionMakerBase_h


#include <memory>
#include <unordered_set>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
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

#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"


#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/Track/interface/SimTrack.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Provenance/interface/ParameterSetID.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// #include "../interface/ParmInfo.h"

#include "Analysis/HitAnalyzer/interface/ScalarPotentialFieldCorrection.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "TrackPropagation/Geant4e/interface/Geant4ePropagator.h"
#include "TrackPropagation/Geant4e/interface/CvhMasterThread.h"
#include "TrackPropagation/Geant4e/interface/CvhMasterRecord.h"
#include "TrackPropagation/Geant4e/interface/CvhWorker.h"


#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH2D.h"
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/AutoDiff>
#include<Eigen/StdVector>
#include <iostream>
#include <functional>

#include "TrackPropagation/Geant4e/interface/MaterialGroupModel.h"
#include "Analysis/HitAnalyzer/interface/ScalarPotFieldModeProvider.h"

using namespace Eigen;

constexpr unsigned int max_n = 25; //!< In order to avoid use of dynamic memory

//too big for stack :(
// typedef Matrix<double, Dynamic, Dynamic, 0, 5*max_n, 5*max_n> GlobalParameterMatrix;
// typedef Matrix<double, Dynamic, 1, 0, 5*max_n, 1> GlobalParameterVector;
// typedef Matrix<double, Dynamic, Dynamic, 0, 5*max_n, 2*max_n> AlignmentJacobianMatrix;
// typedef Matrix<double, Dynamic, Dynamic, 0, 5*max_n, 2*max_n> TransportJacobianMatrix;

typedef Matrix<double, 5, 5> Matrix5d;
typedef Matrix<double, 5, 1> Vector5d;
typedef Matrix<double, 6, 1> Vector6d;
typedef Matrix<float, 5, 5> Matrix5f;
typedef Matrix<float, 5, 1> Vector5f;
typedef Matrix<unsigned int, Dynamic, 1> VectorXu;

typedef MatrixXd GlobalParameterMatrix;
typedef VectorXd GlobalParameterVector;
typedef MatrixXd AlignmentJacobianMatrix;
typedef MatrixXd TransportJacobianMatrix;
typedef MatrixXd ELossJacobianMatrix;

//
// class declaration
//

template<typename T>
using evector = std::vector<T, Eigen::aligned_allocator<T>>;


namespace std{
  template <>
  struct hash<std::pair<unsigned int, unsigned int>> {
      size_t operator() (const std::pair<unsigned int, unsigned int>& s) const {
        unsigned long long rawkey = s.first;
        rawkey = rawkey << 32;
        rawkey += s.second;
        return std::hash<unsigned long long>()(rawkey); 
      }
  };
}

namespace pat {
  class Muon;
}

//double active scalar for autodiff grad+hessian
//template arguments are datatype, and size of gradient
//In this form the gradients are null length unless/until the variable
//is initialized with init_twice_active_var (though the corresponding
//memory is still used on the stack)
template<typename T, int N>
using AANT = AutoDiffScalar<Matrix<AutoDiffScalar<Matrix<T, N, 1>>, Dynamic, 1, 0, N, 1>>;

class ResidualGlobalCorrectionMakerBase
    : public edm::stream::EDProducer<>
{
public:
  // The shared CVH G4 master (CvhMasterThread) is now an EventSetup product on
  // CvhMasterRecord (built once per job by CvhMasterESProducer). Each maker
  // simply esConsumes it (cvhMasterToken_) and, on the first produce() of each
  // TBB worker thread, bootstraps that thread's G4 state from the shared
  // master via worker_->ensureInitialized. This replaces the old
  // edm::GlobalCache<CvhMasterThread> ownership (which was per-module-label and
  // therefore could not be shared across multiple CVH producers in one job).
  ResidualGlobalCorrectionMakerBase(const edm::ParameterSet &);
  ~ResidualGlobalCorrectionMakerBase();

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

protected:

  virtual void beginStream(edm::StreamID) override;
// virtual void analyze(const edm::Event &, const edm::EventSetup &) override;
  virtual void endStream() override;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;

  GloballyPositioned<double> surfaceToDouble(const Surface &surface) const;

  // re-place hits on garbage-shifted modules at the path position implied by
  // their REPAIRED surfaces: the stored hit order comes from the garbage
  // constants and can put the repaired surface behind its predecessor,
  // which aborts the forward-only Geant4e propagation. Sane hits keep their
  // stored order; each flagged hit is inserted where its projection along
  // the track momentum fits between them.
  void reorderGarbageShiftHits(TransientTrackingRecHit::RecHitContainer &hits,
                               const math::XYZVector &trackmom) const;

  GloballyPositioned<double> surfaceToDouble(const Surface &surface, const Basic3DVector<double> &gz) const;

  void applyAlignment(GloballyPositioned<double> &surface, const DetId &detid) const;

  Matrix<double, 6, 1> globalToLocal(const Matrix<double, 7, 1> &state, const GloballyPositioned<double> &surface) const;

  Matrix<double, 7, 1> localToGlobal(const Matrix<double, 6, 1> &localstate, const GloballyPositioned<double> &surface) const;

  Matrix<double, 5, 5> curv2localJacobianAltelossD(const Matrix<double, 7, 1> &state, const MagneticField *field, const GloballyPositioned<double> &surface, double dEdx, double mass, const Eigen::Vector3d &dB = Eigen::Vector3d::Zero()) const;

  Matrix<double, 5, 5> curv2localhybridJacobianAltelossD(const Matrix<double, 7, 1> &state, const MagneticField *field, const GloballyPositioned<double> &surface, double dEdx, double mass, const Eigen::Vector3d &dB = Eigen::Vector3d::Zero()) const;

  Matrix<double, 5, 11> curv2localJacobianAlignmentD(const Matrix<double, 7, 1> &state, const MagneticField *field, const GloballyPositioned<double> &surface, double dEdx, double mass, const Eigen::Vector3d &dB = Eigen::Vector3d::Zero()) const;

  Matrix<double, 6, 5> curv2cartJacobianAltD(const Matrix<double, 7, 1> &state) const;

  Matrix<double, 5, 6> hybrid2curvJacobianD(const Matrix<double, 7, 1> &state, const MagneticField *field, const Eigen::Vector3d &dB = Eigen::Vector3d::Zero()) const;

  Matrix<double, 7, 1> pca2cart(const Matrix<double, 5, 1> &statepca, const reco::BeamSpot &bs) const;

  Matrix<double, 5, 1> cart2pca(const Matrix<double, 7, 1> &state, const reco::BeamSpot &bs) const;

  Matrix<double, 5, 5> pca2curvJacobianD(const Matrix<double, 7, 1> &state, const MagneticField *field, const reco::BeamSpot &bs, const Eigen::Vector3d &dB = Eigen::Vector3d::Zero()) const;

  Matrix<double, 6, 5> pca2cartJacobianD(const Matrix<double, 7, 1> &state, const reco::BeamSpot &bs) const;

  std::array<Matrix<double, 7, 1>, 2> twoTrackPca2cart(const Matrix<double, 10, 1> &statepca) const;

  Matrix<double, 10, 1> twoTrackCart2pca(const Matrix<double, 7, 1> &state0, const Matrix<double, 7, 1> &state1) const;

  Matrix<double, 10, 10> twoTrackPca2curvJacobianD(const Matrix<double, 7, 1> &state0, const Matrix<double, 7, 1> &state1, const MagneticField *field, const Eigen::Vector3d &dB0 = Eigen::Vector3d::Zero(), const Eigen::Vector3d &dB1 = Eigen::Vector3d::Zero()) const;


// Matrix<double, 10, 1> twoTrackCart2pcaJacobianD(const Matrix<double, 7, 1> &state0, const Matrix<double, 7, 1> &state 1, const MagneticField *field, const reco::BeamSpot &bs, double dBz0 = 0., double dBz1 = 0.);

  Matrix<double, 2, 1> localPositionConvolutionD(const Matrix<double, 7, 1>& state, const Matrix<double, 5, 5> &curvcov, const GloballyPositioned<double> &surface) const;

  template <typename T>
  void init_twice_active_var(T &ad, const unsigned int d_num, const unsigned int idx) const;
  
  template <typename T>
  void init_twice_active_null(T &ad, const unsigned int d_num) const;
  
  
  Matrix<double, 1, 6> massJacobianAltD(const Matrix<double, 7, 1> &state0, const Matrix<double, 7, 1> &state1, double dmass) const;
  
  Matrix<double, 6, 6> massHessianAltD(const Matrix<double, 7, 1> &state0, const Matrix<double, 7, 1> &state1, double dmass) const;

  Matrix<double, 1, 6> massinvsqJacobianAltD(const Matrix<double, 7, 1> &state0, const Matrix<double, 7, 1> &state1, double dmass) const;
  
  Matrix<double, 6, 6> massinvsqHessianAltD(const Matrix<double, 7, 1> &state0, const Matrix<double, 7, 1> &state1, double dmass) const;  

  // ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<Trajectory>> inputTraj_;
// edm::EDGetTokenT<std::vector<reco::GenParticle>> GenParticlesToken_;
  edm::EDGetTokenT<edm::View<reco::Candidate>> GenParticlesToken_;
  edm::EDGetTokenT<math::XYZPointF> genXyz0Token_;
  edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_;
  edm::EDGetTokenT<std::vector<int>> genParticlesBarcodeToken_;
// edm::EDGetTokenT<TrajTrackAssociationCollection> inputTrack_;
  
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> pileupSummaryToken_;
  
  edm::EDGetTokenT<reco::TrackCollection> inputTrack_;
  edm::EDGetTokenT<reco::TrackCollection> inputTrackOrig_;
  // Optional fast-path: when srcCandidates is set, the per-pair iteration
  // walks this VertexCompositeCandidateCollection (each candidate's two
  // RecoChargedCandidate daughters) instead of doing a j>i outer-product
  // over inputTrackOrig_. inputCandidatesTag_.label().empty() is the
  // runtime switch.
  edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> inputCandidates_;
  edm::InputTag inputCandidatesTag_;
  // Optional Stage-2 per-row B+ candidate index. When `bCandIdxSrc` is
  // configured non-empty, a parallel std::vector<int> (same length + order
  // as the per-row source collection -- trackOrigH for the single-track
  // maker, inputCandidates_ for the two-track maker's fast path) is read
  // and an `int bCandIdx` branch is added to the per-row output tree.
  // Empty => no branch added; bCandIdx member stays at sentinel -1.
  edm::EDGetTokenT<std::vector<int>> bCandIdxToken_;
  edm::InputTag bCandIdxSrcTag_;
  edm::EDGetTokenT<std::vector<int> > inputIndices_;
  edm::EDGetTokenT<reco::BeamSpot> inputBs_;
// edm::EDGetTokenT<std::vector<PSimHit>> inputSimHits_;
  std::vector<edm::EDGetTokenT<std::vector<PSimHit>>> inputSimHits_;
  edm::EDGetTokenT<std::vector<SimTrack>> inputSimTracks_;
  
// edm::EDGetTokenT<reco::MuonCollection> inputMuons_;
  edm::EDGetTokenT<edm::View<reco::Muon>> inputMuons_;
  edm::EDGetTokenT<int> inputGeometry_;

  edm::EDGetTokenT<edm::Association<std::vector<pat::Muon>>> inputMuonAssoc_;
  bool doMuonAssoc_;
  
  edm::EDGetTokenT<edm::TriggerResults> inputTriggerResults_;
  bool doTrigger_;
  edm::ParameterSetID triggerNamesId_;
  std::vector<std::string> triggers_;
  std::vector<std::size_t> triggerIdxs_;
  std::vector<int> triggerDecisions_;

  edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> globalGeometryToken_;
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackerGeomIdealToken_;
  edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> trackerTopologyToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magfieldToken_;
  // Event-transition tokens used by derived classes' produce() (Phase 2
  // residual makers query the geometry / topology per event). The
  // BeginRun-scoped tokens above serve beginRun() in this base class.
  edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> globalGeometryEventToken_;
  edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> trackerTopologyEventToken_;

  // The shared CVH G4 master, consumed per-event; worker_->ensureInitialized
  // attaches this thread's G4 state to it on the first produce() per thread.
  edm::ESGetToken<CvhMasterThread, CvhMasterRecord> cvhMasterToken_;
  
  
  std::vector<std::string> corFiles_;
  
// SiStripClusterInfo siStripClusterInfo_;

  
  TFile *fout = nullptr;
  TTree *tree = nullptr;
// TTree *runtree;
// TTree *gradtree;
// TTree *hesstree;

  float trackEta;
  float trackPhi;
  float trackPt;
  float trackPtErr;
  float trackCharge;

  float normalizedChi2;
  
  float genPt;
  float genEta;
  float genPhi;
  float genCharge;
  
  float genX;
  float genY;
  float genZ;
  
  unsigned int nHits;
  unsigned int nValidHits;
  unsigned int nValidPixelHits;
  unsigned int nParms;
  unsigned int nJacRef;
  unsigned int nSym;

  // Factored Hessian storage (fillGradsFactored_): the reduced Hessian
  // wrt the global params has rank exactly ndof (the measurement content
  // incl. constraint rows; the deweighted strip coordinates carry exactly
  // zero weight) << nParms, so it is stored as B (nRank x nParms,
  // row-major, H = B^T B summed over rows) instead of the packed dense
  // upper triangle. nRank = min(ndof, nParms) by counting -- not by
  // eigenvalue threshold, since the numerical-noise tail of the
  // double-precision profiling overlaps the smallest genuine modes.
  // nFactor = nRank*nParms is the flat branch dimension.
  unsigned int nRank;
  unsigned int nFactor;
  // relative eigenvalue mass dropped by the rank truncation,
  // sum(dropped lambda)/sum(kept lambda) -- monitoring quantity
  float hessdroppedmass;
  // spectral-gap monitor: lambda_{first dropped}/lambda_{last kept}
  // (0 if nothing dropped). Healthy candidates sit at ~1e-9..1e-4;
  // a value approaching 1 would flag a defect in the rank counting
  // (a genuine eigenmode being dropped).
  float hessrankgap;

  // Stage-2 per-row B+ candidate index (filled per Fill() call when the
  // optional bCandIdxToken_ is configured; -1 sentinel otherwise).
  int bCandIdx = -1;
  
  float gradmax;
  float hessmax;
  
  std::array<float, 5> trackOrigParms;
  std::array<float, 25> trackOrigCov;
  
  
  std::array<float, 5> trackParms;
  std::array<float, 25> trackCov;
  
  std::array<float, 5> refParms_iter0;
  std::array<float, 25> refCov_iter0;

  std::array<float, 5> refParms_iter2;
  std::array<float, 25> refCov_iter2;
  
  std::array<float, 5> refParms;
  std::array<float, 25> refCov;
  
  std::array<float, 5> genParms;
  
  std::vector<float> gradv;
  std::vector<float> jacrefv;
  std::vector<unsigned int> globalidxv;
  std::vector<unsigned int> globalidxvfinal;
  
  std::vector<float> hesspackedv;
  std::vector<float> hessfactorv;

  std::vector<unsigned int> hitidxv;
  std::vector<float> dxrecgen;
  std::vector<float> dyrecgen;
  std::vector<float> dxsimgen;
  std::vector<float> dysimgen;
  std::vector<float> dxsimgenconv;
  std::vector<float> dysimgenconv;
  std::vector<float> dxsimgenlocal;
  std::vector<float> dysimgenlocal;
  std::vector<float> dxrecsim;
  std::vector<float> dyrecsim;
  std::vector<float> dxerr;
  std::vector<float> dyerr;
  
  std::vector<int> clusterSize;
  std::vector<int> clusterSizeX;
  std::vector<int> clusterSizeY;
  std::vector<int> clusterCharge;

  std::vector<int> stripsToEdge;
  
  std::vector<int> clusterChargeBin;
  std::vector<int> clusterOnEdge;
  
  std::vector<float> clusterProbXY;
  std::vector<float> clusterSN;
  
  std::vector<float> dxreccluster;
  std::vector<float> dyreccluster;
  
  std::vector<float> simlocalqop;
  std::vector<float> simlocaldxdz;
  std::vector<float> simlocaldydz;
  std::vector<float> simlocalx;
  std::vector<float> simlocaly;
  
  std::vector<float> simlocalqopprop;
  std::vector<float> simlocaldxdzprop;
  std::vector<float> simlocaldydzprop;
  std::vector<float> simlocalxprop;
  std::vector<float> simlocalyprop;

  std::vector<float> landauDelta;
  std::vector<float> landauW;
  
  std::vector<float> localqop;
  std::vector<float> localdxdz;
  std::vector<float> localdydz;
  std::vector<float> localx;
  std::vector<float> localy;
  
  std::vector<float> localqoperr;
  std::vector<float> localdxdzerr;
  std::vector<float> localdydzerr;
  std::vector<float> localxerr;
  std::vector<float> localyerr;
  
  std::vector<float> hitlocalx;
  std::vector<float> hitlocaly;
  
  std::vector<float> localqop_iter;
  std::vector<float> localdxdz_iter;
  std::vector<float> localdydz_iter;
  std::vector<float> localx_iter;
  std::vector<float> localy_iter;

  std::vector<float> dxrecgen_iter;
  std::vector<float> dyrecgen_iter;
  
  std::vector<float> localqoperralt;

  std::vector<float> localphi;
  std::vector<float> hitphi;
  
  std::map<std::pair<int, DetId>, unsigned int> detidparms;
  std::vector<std::pair<int, DetId>> detidparmsrev;
  std::map<DetId, std::array<int, 3>> detidlayermap;
  
  std::map<DetId, ReferenceCountingPointer<Plane>> surfacemap_;
  std::map<DetId, GloballyPositioned<double>> surfacemapD_;
  std::map<DetId, GloballyPositioned<double>> surfacemapIdealD_;
  std::map<DetId, Eigen::Matrix<double, 2, 2>> rgluemap_;

  std::vector<std::vector<double>> corparmsIncremental_;
  std::vector<double> corparms_;
  
  unsigned int run;
  unsigned int lumi;
  unsigned long long event;
  
  std::vector<double> gradagg;
  
  std::unordered_map<std::pair<unsigned int, unsigned int>, double> hessaggsparse;
  
  bool fitFromGenParms_;
  bool fitFromSimParms_;
  bool fillTrackTree_;
  bool fillGrads_;
  bool fillGradsFactored_;
  bool fillJac_;
  bool fillRunTree_;
  bool alignGlued_ = true;
  double gluedGarbageTiltThreshold_ = 0.05;
  double moduleGarbageShiftThreshold_ = 0.4;
  // modules failing the local-consensus displacement check; their hits are
  // either re-ordered to match the repaired surface positions (default) or
  // excluded from the fit entirely, per garbageShiftReorderHits_
  std::unordered_set<uint32_t> garbageShiftModules_;
  bool garbageShiftReorderHits_ = true;
  
  bool debugprintout_;
  
  bool doGen_;
  bool doSim_;
  // Rung-E closure: substitute simulated hit positions for cluster
  // positions in the fit (needs doSim=True for the sim-hit matching)
  bool fitSimHitPositions_ = false;
  bool doMuons_;
  bool requireGen_;
  
  bool bsConstraint_;
  
  bool applyHitQuality_;

  int genMatchPdgId_ = 13;
  double genMatchPtWindow_ = 0.5;

  // Keep pixel hits whose cluster touches the sensor boundary (isOnEdge) in
  // the fit instead of demoting them to inactive. The sizeX CPE-quality
  // requirement (below) is unaffected. Default false = legacy behaviour.
  bool keepPixelEdgeHits_ = false;

  // Minimum pixel cluster size in x for a hit to stay in the fit
  // (CPE x-resolution needs charge sharing between >=2 pixels).
  // Default 2 = legacy sizeX>1 cut; 1 admits all clusters.
  int pixelMinSizeX_ = 2;

  // Register the pixel pathological-hit class correction parameters
  // (parmtypes 16-21, per pixel module, mean/diff basis for the edge
  // classes) and emit their Jacobian columns in the fit. Meant to be
  // used together with keepPixelEdgeHits=True pixelMinSizeX=1 so the
  // pathological hits are actually in the fit. Default false = catalog
  // unchanged.
  bool pixelHitClassCorrections_ = false;
  int genMatchPdgId_ = 13;

  // Physics parameterization of the local-x drift effects: replace the
  // empirical parmtypes 16 (edge-x-mean) and 20 (sizeX1) with a single
  // delta-tan(thetaL) parameter per pixel module (parmtype 22) whose
  // Jacobian column is attached to EVERY valid pixel hit with the
  // class-dependent response weight w (size-1: 1, x-edge: lorentzWedge,
  // regular: lorentzWclean) and per-hit scale t/2 (t = sensor thickness).
  // Coupling the regular hits makes dtanLA separable from the local-x
  // alignment translation (uniform response) in a simultaneous fit.
  // Requires pixelHitClassCorrections.
  bool pixelLorentzParam_ = false;
  double lorentzWclean_ = 0.69;   // measured, 100k twins: +0.686 +- 0.013
  double lorentzWsize1_ = 0.07;   // measured: +0.066 +- 0.033 (pixel-center quantisation)
  double lorentzWedge_ = 0.75;    // measured: lo +0.77+-0.19, hi +0.74+-0.18 (side-symmetric)
  // Injection test: shift every valid pixel hit's local-x position by
  // 0.5*t*w_inj(class)*injectLorentzTan, simulating a true Lorentz-angle
  // mismatch of the given size. injectLorentzWclean lets the injected
  // clean-hit response differ from the model one (response-model error
  // study); < -900 = use lorentzWclean.
  double injectLorentzTan_ = 0.;
  double injectLorentzWclean_ = -999.;

  bool doRes_ = false;
  bool useIdealGeometry_ = false;

  // Label of the magnetic field ESProducer the CVH refit should consume
  // (e.g. "ScalarPot3DMf"). Routed to the residual-makers via the cfi
  // MagneticFieldLabel parameter.
  std::string fieldlabel_;

  // Scalar-potential B-field correction. Replaces the per-module dBz block
  // with the spherical-harmonic coefficients of the magnetic scalar potential
  // loaded from a coefficient dump file (mfs/dump_coeffs_for_cmssw.py).
  std::string scalarPotentialInitFile_;
  std::unique_ptr<ana_hitanalyzer::ScalarPotentialFieldCorrection> fieldCorrection_;

  // Global material model (doc/global-material-model-plan.md).
  // materialGroupsFile loads a grouping-tier rules file (Phase A
  // validation hook usable on its own); globalMaterialModel=true
  // additionally REPLACES the per-module material parameters (parmtype 7)
  // with parmtype-15 sentinel entries, one per group (exclusive switch).
  bool globalMaterialModel_ = false;
  std::string materialGroupsFile_;
  std::unique_ptr<MaterialGroupModel> matModel_;
  std::vector<unsigned int> matGroupGlobalIdx_;  // groupId -> corparms_ index

  // Per-step field modes (leg-structure-free attribution): when true, the
  // scalar-potential correction is applied per Geant4 step via the
  // provider, and the per-mode derivative columns come from the propagator
  // instead of the per-leg chain rule.
  bool perStepFieldModes_ = false;
  std::unique_ptr<ana_hitanalyzer::ScalarPotFieldModeProvider> fieldModeProvider_;

  // Skip hitless module surfaces in the fit hit list (dead-module
  // placeholders and quality-demoted hits): valid only with the global
  // material model (per-module leg attribution would otherwise break).
  bool skipHitlessSurfaces_ = false;

  // Numerical-FD closure (debug only; one-shot per job).
  bool runFDClosure_ = false;
  double epsilonFDClosure_ = 1e-4;
  // mutable so the per-job "did the test" flag can be set inside
  // analyze() without making the whole maker non-const.
  mutable bool didFDClosure_ = false;
  
  float dxpxb1;
  float dypxb1;
  
  float dxttec9rphisimgen;
  float dyttec9rphisimgen;
  float dxttec9rphi;
  float dxttec9stereo;

  float dxttec4rphisimgen;
  float dyttec4rphisimgen;
  float dxttec4rphirecsim;
  
  float dxttec4rphi;
  float dxttec4stereo;
  
  float simlocalxref;
  float simlocalyref;
  
  float simtestz;
  float simtestvz;
  float simtestzlocalref;
  float simtestrho;
  float simtestdx;
  float simtestdxrec;
  float simtestdy;
  float simtestdyrec;
  float simtestdxprop;
  float simtestdyprop;
  unsigned int simtestdetid;
  
  std::vector<float> rx;
  std::vector<float> ry;
  
  std::vector<float> deigx;
  std::vector<float> deigy;
  
  float edmval;
  float edmvalref;
  float deltachisqval;
  unsigned int niter;
  // per-track count of GN iterations where the clamp caught a q/p sign
  // crossing (charge-flip protection through p->inf). >0 flags a track that
  // "wanted" the opposite charge -> trigger the two-hypothesis second fit.
  unsigned int nChargeFlipProtect;
  // 1 if the two-hypothesis fit kept the OPPOSITE charge (opposite converged
  // with lower chi2 than the nominal seed), 0 otherwise.
  unsigned int chargeHypFlipped;
  
  float chisqval;
  unsigned int ndof;

  float genweight;
  
  int Pileup_nPU = 0;
  float Pileup_nTrueInt = 0.;

  float genl3d = -99.;
  // fit-transmission study: true momentum (pabs, GeV) at the FIRST and
  // LAST sim hit matched along the trajectory (doSim); the difference is
  // the track's true in-tracker energy loss
  float simPabsFirst = -99.;
  float simPabsLast = -99.;
  
  std::vector<float> gradchisqv;
  // log-det trace term nu = tr(dV_i R) per parameter, stored DIRECTLY so
  // small values (ionization: nu ~ 1e-7) keep full relative precision --
  // the offline reconstruction gradv - gradchisqv loses them to float32
  // cancellation. Nonzero only for resolution parameters (doRes).
  std::vector<float> gradllv;

  // Physics-CF export (doRes + grads): per Geant4 step with a nonzero
  // ionization-fluctuation contribution, the Urban-model compound-Poisson
  // parameters + the q/p mapping coefficient, tagged by the parmtype-11
  // global parameter index of the leg. 11 floats per step:
  // [regime, gsig2, a1, e1, a2, e2, a3, e0, tmax, scaling, cs]
  // (energies MeV; cs = Etot/p^3 in GeV^-2, see Geant4ePropagator).
  std::vector<unsigned int> ioniurbanidx;
  std::vector<float> ioniurbanv;

  // Phase B analogue for multiple scattering: per Geant4 step, raw
  // material/kinematic data for the offline Moliere compound-Poisson tail
  // model, tagged by the parmtype-10 global parameter index of the leg.
  // 8 floats per step: [effZ, effA, x(g/cm2), p(GeV), beta,
  // thp2-as-in-Q, d/X0, materialGroup] (see Geant4ePropagator::MoliereMsStep).
  std::vector<unsigned int> msmoliidx;
  std::vector<float> msmoliv;

  // Per resolution entry (leg), the exact eigenvalues of the block
  // quadratic form dV_b^{1/2} R_bb dV_b^{1/2} (descending, zero-padded to
  // 5), tagged by the entry's global parameter index. Replaces the
  // two-moment r_eff approximation in the offline CF fits and breaks the
  // k/s_est degeneracy. Validation: sum over legs of sum(lambda) equals
  // the parameter's gradllv entry.
  std::vector<unsigned int> reseigidx;
  std::vector<float> reseigv;

  // q/p influence weights: the linearized-fit response to the noise vector
  // n is delta x_ref = C F^T Vinv n, so the q/p row w = (C e_qop)^T F^T Vinv
  // decomposes the fitted q/p error into independent per-block noise
  // contributions delta(q/p) = sum_b w_b^T n_b. Per resolution entry
  // (aligned with reseigidx): resinfv = the raw dof weights w_b in fit
  // units (zero-padded to 5), resinfvarv = the variance contribution
  // w_b^T dV_b w_b. resinfcov = sum_b resinfvarv, which equals refCov(0,0)
  // exactly when every noise block carries a resolution entry (coverage
  // check for the per-track CF-product resolution prediction).
  std::vector<float> resinfv;
  std::vector<float> resinfvarv;
  float resinfcov = 0.;

  // Generalized influence export: per resolution entry the 5x5 (row-major,
  // dof-padded) matrix B_b = M_b dV_b^{1/2}, where M_b (5 x nb) is the
  // response of the 5 reference parameters to the block's noise dofs.
  // For ANY linear functional a of the reference state (mass Jacobian,
  // pT, ...), the block's signed noise weights are a^T B_b (in
  // dV^{1/2}-standardized units) and sum_b |a^T B_b|^2 = a^T C a exactly
  // (per-candidate identity against refCov). Signs preserve the Landau
  // skew of the ionization contribution.
  std::vector<float> resinfbv;


  TH2D *hetaphi = nullptr;

  std::string outprefix;

// bool filledRunTree_;

  // MT support: each stream owns its own Geant4ePropagator clone, plus a
  // CvhWorker that does per-thread G4 setup (world + magnetic field) on
  // first produce(). The clone is fetched / built lazily in derived
  // produce() AFTER worker_->ensureInitialized has installed the world on
  // this TBB worker thread; before then, allocating the propagator's fluct
  // would abort in G4WentzelVIModel::Initialise (NULL world). See
  // geant4e_multithreading_exploration.txt.
  std::unique_ptr<CvhWorker> worker_;
  std::unique_ptr<Geant4ePropagator> streamPropagator_;

  // Wire the stream's CLHEP engine into Geant4's thread-local engine at the
  // top of every produce() call (G4Random::setTheEngine is itself
  // thread-local; setting it once in beginStream wouldn't survive TBB
  // thread migration). Derived classes call this helper from produce().
  // Returns the engine reference so callers may seed any private RNG too.
  CLHEP::HepRandomEngine &setG4RandomEngineForStream(edm::StreamID) const;

};

template <typename T>
void ResidualGlobalCorrectionMakerBase::init_twice_active_var(T &ad, const unsigned int d_num, const unsigned int idx) const {
  // initialize derivative direction in value field of outer active variable
  ad.value().derivatives() = T::DerType::Scalar::DerType::Unit(d_num,idx);
  // initialize derivatives direction of the variable
  ad.derivatives() = T::DerType::Unit(d_num,idx);
  // initialize Hessian matrix of variable to zero
  for(unsigned int idx=0;idx<d_num;idx++){
    ad.derivatives()(idx).derivatives()  = T::DerType::Scalar::DerType::Zero(d_num);
  }
}

template <typename T>
void ResidualGlobalCorrectionMakerBase::init_twice_active_null(T &ad, const unsigned int d_num) const {
  // initialize derivative direction in value field of outer active variable
  ad.value().derivatives() = T::DerType::Scalar::DerType::Zero(d_num);
  // initialize derivatives direction of the variable
  ad.derivatives() = T::DerType::Zero(d_num);
  // initialize Hessian matrix of variable to zero
  for(unsigned int idx=0;idx<d_num;idx++){
    ad.derivatives()(idx).derivatives()  = T::DerType::Scalar::DerType::Zero(d_num);
  }
}

#endif
