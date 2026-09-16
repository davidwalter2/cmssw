#ifndef HitAnalyzer_ResidualGlobalCorrectionMakerBase_h
#define HitAnalyzer_ResidualGlobalCorrectionMakerBase_h


#include <cmath>
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
#include "SimDataFormats/Vertex/interface/SimVertex.h"

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
#include "TrackPropagation/Geant4e/interface/CvhCfExponents.h"
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
  // Geant4 vertex collection, consumed only for the decay-truth block
  // (doSimDecayTruth_). Kept separate from inputSimHits_/doSim_ because the
  // B->J/psi+X ALCARECO keeps SimTracks+SimVertices but NOT the PSimHit
  // collections, so the full doSim_ path cannot be switched on there.
  edm::EDGetTokenT<std::vector<SimVertex>> inputSimVertices_;

// edm::EDGetTokenT<reco::MuonCollection> inputMuons_;
  edm::EDGetTokenT<edm::View<reco::Muon>> inputMuons_;
  edm::EDGetTokenT<int> inputGeometry_;

  edm::EDGetTokenT<edm::Association<std::vector<pat::Muon>>> inputMuonAssoc_;
  bool doMuonAssoc_;
  // Ref-safe muon match against the ALCARECO loose-muon collection: the
  // muon's own bestTrack()/innerTrack() refs point into generalTracks, which
  // the ALCARECO does not keep, so dereferencing them throws. Only quantities
  // stored directly on the muon (selector bits, type bits, p4) may be read.
  // See the matching code in the G4e maker.
  bool doMuonTrackAssoc_ = false;
  
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
  // pdgId of the winner of the dR competition and its dR to the track:
  // the match-quality handles the decay studies need offline.
  int genPdgId = 0;
  float genDR = -99.f;
  
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
  // The VARIANCE (log-det) block of the Hessian, shipped SEPARATELY from
  // `hessfactorv` (see the export doc).  `hesspackedv` is complete on its own;
  // `hessfactorv` factors only the MEAN block `2 J^T R J`, whose rank is
  // exactly ndof, and this is what has to be added to it:
  //
  //     hess = B^T B + scatter(hessvarpackedv on hessvaridxv)
  //
  // The reason is volume.  The variance block is a Gram matrix of the
  // ncons x ncons objects R^1/2 dV_i R^1/2, so its rank is the NUMBER of
  // variance parameters and it does not compress: carrying it as extra rows
  // of B costs nvar*nParms floats (+57 % of a production candidate at family
  // 15 alone, +429 % with 8-11), while its own packed triangle costs
  // nvar*(nvar+1)/2 -- a factor ~20 less.
  //
  // The presence of `hessvaridxv` IS the flag: in a file without it
  // `hessfactorv` is complete, in a file with it this addition is required.
  // `globalfit/extract.py` does it.
  unsigned int nHessVar = 0;
  std::vector<unsigned int> hessvaridxv;   // columns, ascending, into nParms
  std::vector<float> hessvarpackedv;       // upper triangle, row-major
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
  // transient, per fit iteration: block -> valid-hit index (see reshitidx)
  std::vector<int> resvalidhit_;
  // Family (parmtype) of each resolution entry, parallel to `resvalidhit_`
  // and to the makers' local `resglobidx`: 8/9 = hit, 10 = multiple
  // scattering, 11 = ionization. Pushed at the same four sites as
  // `resvalidhit_`, because `cvhcf` has to tell an MS block from an
  // ionization one to pool them the way the offline `extract()` does, and the
  // valid-hit index only separates hits from material.
  std::vector<int> resfamily_;
  // TWO-TRACK ONLY: index of the first resolution block belonging to leg 1,
  // so a block's LEG (and with it its charge) is `ires >= resLegStart_[1]`.
  // The vertex functional's ionization sign is charge-signed per block --
  // the two legs have opposite charges, which the mass functional's single
  // `ioniSign = -1` hides.
  std::array<std::size_t, 2> resLegStart_ = {{0, 0}};

  // Multipliers on the assigned hit covariance (1.0 = the CPE value as-is).
  double hitCovScalePixel_ = 1.0;
  double hitCovScaleStrip_ = 1.0;

  // Hit-resolution study exports (per valid hit, filled on iteration 0
  // alongside dxerr/dxrecsim). hitUProj is the CPE's OWN independent
  // variable -- the track path across the sensor projected onto the
  // measurement direction, in strip-pitch units, INCLUDING the Lorentz
  // drift -- so a pull binned in it is binned in the parametrisation being
  // tested rather than in a proxy. -99 for pixels.
  // Number of PSimHits on this module belonging to the fitted sim track. > 1
  // means the closest-to-the-propagated-state tie-break above actually chose.
  std::vector<int> simHitNCand;
  std::vector<unsigned int> hitDetId;
  std::vector<float> hitUProj;
  // strip coordinates: reco, truth, and the cluster's first strip, so the
  // true impact point can be referred to the cluster's own lattice
  std::vector<float> hitStripRec;
  std::vector<float> hitStripSim;
  std::vector<int> hitFirstStrip;
  std::vector<float> hitPitch;
  std::vector<float> hitThickness;
  
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
  // Geant4 decay/interaction truth for the gen-matched particle: walks the
  // SimVertex collection for vertices whose parent is the matched SimTrack.
  // Needs doGen_ (the gen match provides the barcode -> SimTrack link) but
  // NOT doSim_ (no PSimHits required).
  bool doSimDecayTruth_ = false;
  // Rung-E closure: substitute simulated hit positions for cluster
  // positions in the fit (needs doSim=True for the sim-hit matching)
  bool fitSimHitPositions_ = false;
  bool doMuons_;
  bool requireGen_;
  
  bool bsConstraint_;
  
  bool applyHitQuality_;

  int genMatchPdgId_ = 13;
  double genMatchPtWindow_ = 0.5;
  // Species allowed to compete for the dR match. When non-empty every listed
  // species is matched and the winner must be genMatchPdgId_, otherwise the
  // track counts as unmatched. Without this a soft gen hadron can win the
  // match to an unrelated (usually J/psi muon) track whenever the pT window
  // is opened up, which is exactly what decay studies have to do.
  std::vector<int> genMatchPdgIds_;
  double genMatchDR_ = 0.1;

  // Keep pixel hits whose cluster touches the sensor boundary (isOnEdge) in
  // the fit instead of demoting them to inactive. The sizeX CPE-quality
  // requirement (below) is unaffected. Default false = legacy behaviour.
  bool keepPixelEdgeHits_ = false;

  // Minimum pixel cluster size in x for a hit to stay in the fit
  // (CPE x-resolution needs charge sharing between >=2 pixels).
  // Default 2 = legacy sizeX>1 cut; 1 admits all clusters.
  int pixelMinSizeX_ = 2;
  // Default 1 = NO cut. sizeY==1 is the GEOMETRIC outcome for a
  // perpendicular track (sizeY ~ 1 + 1.9|cot theta|), unlike
  // sizeX==1 which means the expected Lorentz sharing failed.
  // Provided as an A/B DIAGNOSTIC for the BPix local-x residual,
  // not as a recommended selection: it would drop 12.6% of BPix
  // and 34.6% of FPix hits, concentrated at central eta.
  int pixelMinSizeY_ = 1;

  // Register the pixel pathological-hit class correction parameters
  // (parmtypes 16-21, per pixel module, mean/diff basis for the edge
  // classes) and emit their Jacobian columns in the fit. Meant to be
  // used together with keepPixelEdgeHits=True pixelMinSizeX=1 so the
  // pathological hits are actually in the fit. Default false = catalog
  // unchanged.
  bool pixelHitClassCorrections_ = false;

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

  // CGF q/p IRLS re-centring counters (CVH_CGF_QOP=3). Diagnostics only; the
  // clamp count is the one that matters -- it says how often a residual fell
  // outside the block's own support, where the score is continued rather than
  // evaluated.
  mutable unsigned long long nCgfRecentre_ = 0;
  mutable unsigned long long nCgfClamp_ = 0;
  
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
  
  float simtestz = -99.f;
  float simtestvz = -99.f;
  float simtestzlocalref = -99.f;
  float simtestrho = -99.f;
  float simtestdx = -99.f;
  float simtestdxrec = -99.f;
  float simtestdy = -99.f;
  float simtestdyrec = -99.f;
  float simtestdxprop = -99.f;
  float simtestdyprop = -99.f;
  unsigned int simtestdetid = 0;
  
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
  unsigned int nChargeFlipProtect = 0;
  // 1 if the two-hypothesis fit kept the OPPOSITE charge (opposite converged
  // with lower chi2 than the nominal seed), 0 otherwise.
  unsigned int chargeHypFlipped = 0;
  
  float chisqval;
  // The MARGINAL objective  r^T R r + ln|V| + ln|C|  in double precision,
  // written only under `exportObjective_`.  `chisqval` is the first term
  // alone and is a float, which is 3-4 digits short of what a
  // finite-difference of the log-det gradient needs.
  double objval = 0.;
  // its three pieces, for diagnosing which one a FD mismatch is in
  double objchisq = 0.;
  double objlogdetv = 0.;
  double objlogdetc = 0.;
  // number of exactly-null modes of Vinv dropped from ln|V| (structural:
  // one per deweighted strip coordinate).  Must match between the two
  // arms of a finite difference.
  int objnullv = 0;
  // REFERENCE ENERGY LOSS of the track, and the worst single propagation
  // step's fractional loss. See the two-track maker's `Mu*_dEref` /
  // `Mu*_maxfracloss` for what they are for: a `dE_ref/p < 0.01` quality
  // requirement on the quadratic term's material information, imposable
  // WITHOUT the 430 kB/candidate step records. Two floats.
  float dEref = 0.f;
  float maxfracloss = 0.f;
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

  // THE SCALE APPLIED TO EACH IONIZATION BLOCK.
  //
  // `ioniurbanv`'s gsig2 column is the RECORD's variance; the matrix that
  // entered the fit (and hence `resinfvarv` = w_b^T dV_b w_b) is
  // Q_applied = sc * dQ2_record whenever the CGF substitution ran
  // (CgfQoPMode 1/3; see Geant4ePropagator::cgfQScale). sc == 1.0 exactly
  // otherwise -- CgfQoPMode=0, every two-track fit that pins mode 0, and any
  // leg with no block. Without it the offline `var` normalisation
  // sqrt(v_b / sum_steps gsig2 cs^2) is wrong by sqrt(sc) (~400x in the
  // exponent on mode-1 files).
  //
  // ONE ENTRY PER LEG, tagged with the SAME parmtype-11 global index as that
  // leg's `ioniurbanidx` rows, and 2 floats per entry:
  //
  //     [ sc , nIoniSteps ]
  //
  // `nIoniSteps` is how many `ioniurbanv` rows this leg contributed. It is
  // there because several legs can share one global index (a track crossing
  // the same module twice; ~1/3 of blocks), and the pooled step sum is then
  // sum_legs sc_l * (step sum of leg l), NOT a single scalar times the pooled
  // sum. The legs' step rows are contiguous and in drain order for a given
  // index, so the counts split them exactly -- and the split self-checks
  // (the counts must add up to the number of pooled rows).
  std::vector<unsigned int> ioniqscaleidx;
  std::vector<float> ioniqscalev;

  // RADIATIVE (brems + pair) STEP EXPORT, for the offline
  // radiative CF term (cf_brems_exact.py). One entry per Geant4 step of the
  // leg -- i.e. aligned 1:1 with `msmoliv`, NOT with `ioniurbanv`, which only
  // has rows for steps that produced a fluctuation record -- tagged with the
  // leg's parmtype-11 global index (the q/p block, since this is a q/p noise
  // channel and its offline weight comes from the ionization block).
  //
  // `radstepv`: RADSTEP_STRIDE = 12 floats per step, in exactly the order
  // G4ePropagationExport.cc writes its `radv` branch, so cf_brems_exact.py's
  // R_* column constants, step_spectrum() and rad_exponent() apply verbatim:
  //     [effZ, effA, xg(g/cm2), etot(GeV), p(GeV), d/X0, step(cm),
  //      dedxRad, dedxBrem, dedxPair (all GeV/cm), cs(GeV^-2), stepGroup]
  // Column 11 (`stepGroup`) is the step's material group, which the per-group
  // CF export keys on. Read the stride from `radstepstride`, never assume it.
  // `radstepspecv`: 2*RADSTEP_NV = 96 floats per step, dN/dv SHAPE for
  // bremsstrahlung (48) then pair production (48) on the shared v grid.
  // These are shapes only; each must be renormalized offline to its own
  // process mean dedxBrem/dedxPair (see Geant4ePropagator::RadiativeStep).
  // `radvgrid`: the shared v grid, RADSTEP_NV points, written on EVERY entry
  // (it is constant, so it compresses to nothing) so no reader has to
  // hard-code it. `radstepstride`/`radstepnv` are the two strides as scalar
  // branches for the same reason.
  static constexpr int RADSTEP_STRIDE = 12;
  static constexpr int RADSTEP_NV = Geant4ePropagator::kNRadV;
  std::vector<unsigned int> radstepidx;
  std::vector<float> radstepv;
  std::vector<float> radstepspecv;
  std::vector<float> radvgrid;
  int radstepstride = RADSTEP_STRIDE;
  int radstepnv = RADSTEP_NV;

  // `ioniurbanv`'s stride is not a compile-time constant -- it is 12, or 14
  // with CVH_IONI_EXACTDELTA on -- so it is written as its own scalar branch
  // rather than left to be inferred from `size(v)/size(idx)`.
  int ioniurbanstride = 0;

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
  // Valid-hit index of each resolution block, parallel to reseigidx; -1 for
  // the material (parmtype 10/11) blocks. Lets the offline CF give every
  // parmtype-8/9 block its own per-hit class instead of summing the hit
  // families into a single Gaussian.
  std::vector<int> reshitidx;
  // HIT CLASS of each resolution entry: -1 for a material block, else the
  // index into the canonical 18-class list `hitResClassIndex` implements
  // (`calibration_studies/resolution/hitres_classes.py:CLASSES`).  Filled by
  // the TWO-TRACK maker, whose tree carries no per-hit variables at all; the
  // single-track tree leaves it empty because it exports `hitDetId`,
  // `hitUProj`, `clusterSizeX` and `clusterChargeBin` per hit plus
  // `reshitidx`, which is strictly more information.
  std::vector<short> reshitcls;
  // Sum of v_b over the HIT families (parmtype 8/9).  A SEPARATE scalar
  // deliberately: `resinfcov` keeps its meaning of "the material share", so
  //     cfmass_vgf = (sigma_m^2 - resinfcov)/sigma_m^2
  // still means the TOTAL Gaussian share (hits + beamspot + pointing).  The
  // self-consistent-sigma correction's `a_i = (1 + f_hit) sigma_i/m_i` uses
  // exactly that total, so folding the hit share into `resinfcov` would be a
  // silent physics change.
  float resinfcovhit = 0.f;
  // Sum of v_b over the parmtype-15 MATERIAL-GROUP blocks. They are a
  // RE-PARTITION of the same process noise the parmtype-10/11 blocks carry
  // (`sum_g dQ_g == dQMS + dQI` exactly), so adding them to `resinfcov` would
  // double-count it and break the offline coverage check
  // `|resinfcov/refCov(0,0) - 1| < 5e-3`. They get their own scalar instead,
  // and `resinfcovgrp` should equal the parmtype-10 + parmtype-11 part of
  // `resinfcov` per candidate -- a free closure test of the split.
  float resinfcovgrp = 0.f;
  // Per-hit-class Gaussian variance shares of the candidate: `cf*_hitcls` is
  // the ascending class index and `cf*_hitv` the summed v_b of that class
  // DIVIDED by the functional's sigma^2, i.e. directly the `v_{c,i}` of
  //     v_i(eps) = v_other,i + sum_c H(eps_c) v_{c,i}
  // and `v_other = vgf - sum_c v_c` is formed offline.
  std::vector<short> cfhitclsv;
  std::vector<float> cfhitvv;
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


  // ======================================================================
  // IN-MAKER RESOLUTION-CF EXPONENTS (cvhcf).
  //
  // The four family exponents of this track / candidate on the 64-point
  // standardized tau grid, computed in the doRes pass from the SAME flat
  // arrays above -- see TrackPropagation/Geant4e/interface/CvhCfExponents.h
  // for why they are computed here and not offline.
  //
  // Branch names carry `cfprefix_`: `cfqop_*` in the single-track maker (the
  // q/p functional) and `cfmass_*` in the two-track one (the candidate-mass
  // functional). The two are DIFFERENT functionals of the same blocks -- they
  // differ in the standardization sigma and in the sign the ionization and
  // radiative weights carry -- so they must not share a name.
  std::string cfprefix_ = "cfqop";
  std::vector<float> cfmsv;     // S_ms   (real)
  std::vector<float> cfdelv;    // S_del  (real, delta-ray recoil minus carve)
  std::vector<float> cfiorev;   // Re S_ioni
  std::vector<float> cfioimv;   // Im S_ioni
  std::vector<float> cfradrev;  // Re S_rad
  std::vector<float> cfradimv;  // Im S_rad
  // The GAUSSIAN share of the functional's variance: sum_b v_b over the hit
  // families divided by refCov(0,0) in the single-track maker, and
  // (sigma_m^2 - resinfcov)/sigma_m^2 -- hits + beamspot + pointing -- in the
  // two-track one. It is what `vgf` means in the offline caches.
  float cfvgf = 0.f;
  // 0 iff a registered MS/ionization block had no step rows under its global
  // index, i.e. the model would be missing a block. The offline extractor
  // drops such a track; the flag lets the reader drop it identically.
  bool cfok = false;
  int cfnblock = 0;   // MS + ionization blocks entering the exponents
  int cfnpooled = 0;  // of those, how many pooled more than one leg

  // ---- THE PER-MATERIAL-GROUP SPLIT OF THE SAME EXPONENTS ---------------
  //
  // Every step-level exponent is LINEAR in the step's material amount at
  // fixed composition, so with the fit's influence weights held fixed
  //
  //     S_f(tau; k) = S_f^fixed(tau) + sum_g A(k_g) S_{f,g}(tau),
  //     A(k) = exp(k)  (the propagator's own `matStepFact` convention)
  //
  // is exact and `k = 0` reproduces the flat exponents above. It is what lets
  // the material amounts be floated as the parmtype-15 material-group
  // parameters the hit chi2 already carries, rather than through per-family
  // `k_hit/k_ms/k_ioni/k_rad` knobs.
  //
  // Layout, sparse over the groups the candidate actually touched (~22 of 42
  // on a J/psi gun candidate): `cf*_grp` is the ascending group id, and each
  // family array is `n_grp * kNTau` floats, ROW-MAJOR in (group, tau). The
  // per-candidate group count is `cf*_grp.size()`, so no pointer branch is
  // needed inside an entry. The group ids are parmtype-15 indices of the
  // material-group file named in the runtree.
  //
  // Cost: ~22 groups x 5 families x 64 x float32 = 27 kB/candidate against
  // 1.4 kB flat, hence the switch. A rank-16 PCA of the tau axis would cost
  // 6.8 kB at a relative error of 4e-7 (ms) to 1e-3 (rad); this version
  // exports the RAW rows so that the basis can be chosen from data.
  std::vector<short> cfgrpv;      // ascending material-group id
  // THE FIT'S OWN Q variance of each group, in units of the functional's
  // sigma^2 (Rossi's `thp2` for MS, `ioniSq2` for ionization, nothing for the
  // radiative and delta channels the fit's Q does not have). It is what the
  // Gaussian chi2 the whole exercise is measured against actually assumes,
  // and it cannot be recovered from the exponents -- those carry the MODEL's
  // full Moliere second moment, 14 % larger. `sum_g (vqms+vqio) + vgf == 1`.
  std::vector<float> cfgrpvqmsv, cfgrpvqiov;
  std::vector<float> cfgrpmsv;    // S_ms per group
  std::vector<float> cfgrpdelv;   // S_del per group (single-track only)
  std::vector<float> cfgrpiorev;  // Re S_ioni per group
  std::vector<float> cfgrpioimv;  // Im S_ioni per group
  std::vector<float> cfgrpradrev; // Re S_rad per group
  std::vector<float> cfgrpradimv; // Im S_rad per group
  // max_j |sum_g S_{f,g}(tau_j) - S_f(tau_j)| over the exported families,
  // relative to max_j |S_f|. It is float64 round-off (1e-14 or below) by
  // construction; exporting it means a file can be AUDITED rather than
  // trusted, at 4 bytes.
  float cfgrpclosure = 0.f;

  // Clear + fill the cfgrp* branches from one `cvhcf` result, and set
  // `cfgrpclosure`. Shared so that the single-track and two-track makers
  // cannot lay the arrays out differently.
  void storeCfGroups(const cvhcf::TrackResult &res);
  // The same, into caller-supplied arrays, so a SECOND functional of the same
  // candidate (the vertex one) gets an identical layout for free.
  void storeCfGroupsTo(const cvhcf::TrackResult &res,
                       std::vector<short> &grpv,
                       std::vector<float> &msv,
                       std::vector<float> &delv,
                       std::vector<float> &iorev,
                       std::vector<float> &ioimv,
                       std::vector<float> &radrev,
                       std::vector<float> &radimv,
                       float &closure,
                       bool wantDelta,
                       std::vector<float> *vqms = nullptr,
                       std::vector<float> *vqio = nullptr);

  // ======================================================================
  // THE VERTEX-CONSTRAINT RESIDUAL of the TWO-TRACK fit (`exportVtxResidual_`)
  //
  // State index 6 of the 10-dim vertex-PCA block is the SIGNED track-track
  // PCA distance
  //      theta_6 = n_hat . (x_b - x_a),   n_hat = (p_a x p_b).normalized()
  // (`twoTrackCart2pca`).  Two free helices have 10 vertex parameters, two
  // through a common point have 9, so a common-vertex constraint removes
  // exactly ONE degree of freedom and leaves ONE residual per candidate: the
  // DCA the unconstrained fit finds.  On ideal geometry both muons come from
  // one gen point, so its mean is ZERO BY CONSTRUCTION -- no kernel, no
  // theory, no PDG input.  It is the mass term with a delta kernel at zero.
  //
  //   index 6 FREE (`doVtxConstraint == False`):
  //      sigma_v^2 = C_66,  w_v = Vinv F C e_6,  r_v = statepcaupd[6]
  //   index 6 FROZEN:
  //      sigma_v^2 = 1/(h_66 - h_6f C h_6f^T),  b_6 = -(Vinv F_6).r,
  //      r_v = sigma_v^2 b_6,  w_v = sigma_v^2 (Vinv F_6 - Vinv F C h_6f^T)
  // Both satisfy `sum_b |dV_b^{1/2} w_v,b|^2 == sigma_v^2` exactly (= the
  // gate `vtxvchk`).
  bool exportVtxResidual_ = false;
  float Jpsi_vtxres = 0.f;    // r_v, cm, the raw swap-invariant signed DCA
  float Jpsi_vtxsig = 0.f;    // sigma_v, cm
  float Jpsi_vtxz = 0.f;      // the pull r_v/sigma_v
  float Jpsi_vtxb6 = 0.f;     // -(Vinv F_6).r at convergence (0 when 6 free)
  float Jpsi_vtxdchi2 = 0.f;  // z_v^2 = chi2(constrained) - chi2(unconstrained)
  float Jpsi_vtxvchk = 0.f;   // |sum_b v_b / sigma_v^2 - 1|, fam != 15
  float Jpsi_vtxvgf = 0.f;    // the GAUSSIAN (hit) share of sigma_v^2
  // max over the FREE indices of |g_i| sqrt(C_ii), i.e. the largest remaining
  // Newton step in units of that parameter's own error -- dimensionless, and
  // directly comparable with |z_v| (which is |b_6| sigma_v).
  float Jpsi_vtxbfree = 0.f;
  // Was leg 0 the POSITIVE muon?  theta_6 is invariant under swapping the two
  // legs (n_hat and x_b - x_a both flip), so it is already a well-defined
  // signed DCA and neither `Jpsi_vtxres` nor `Jpsi_d` re-signs it by charge;
  // this flag is here for anyone who wants a charge-ordered convention.
  bool Jpsi_vtxfirstplus = false;
  // The variance shares of sigma_v^2 by FAMILY: hits (parmtype 8/9), multiple
  // scattering (10) and ionization (11).  They sum to 1 (the fit's Q has no
  // radiative or delta channel, which is why the CF is wider than 1).
  float Jpsi_vtxvhit = 0.f, Jpsi_vtxvms = 0.f, Jpsi_vtxvioni = 0.f;
  // the same split for the MASS functional, for the side-by-side composition
  float Jpsi_massvms = 0.f, Jpsi_massvioni = 0.f;
  // THE UNCONSTRAINED MASS, and the covariance element that produces it.
  //
  // Freezing theta_6 at zero is a conditioning of the unconstrained solution,
  //     x_c = x_u - C_u e_6 sigma_v^-2 r_v ,
  // so for the (linearised) mass functional a = dm/dx,
  //     m_u = m_c + cov(m, theta_6) sigma_v^-2 r_v ,
  //     cov(m, theta_6) = a^T C_u e_6 = sigma_v^2 (a_6 - a_f^T C h_f6) ,
  // and a_6 = 0 because the mass depends on the two momenta only. So the
  // slope is just -(C a_f).h_f6 and sigma_v^2 cancels out of the mass:
  // the analysis can form EITHER mass from one constrained fit.
  // With index 6 free the fit already reports the unconstrained mass, and
  // `Jpsi_mass_unc` is `Jpsi_mass` while `cov(m, theta_6) = a_f^T C e_6`
  // is still written out.
  float Jpsi_mass_unc = -99.f;    // GeV
  float Jpsi_covmassvtx = 0.f;    // GeV cm
  bool Jpsi_vtxfree = false;  // was index 6 a free parameter?
  bool Jpsi_vtxok = false;
  // The convention gate: the general per-block ionization sign rule applied
  // to the MASS influence must return -1 on every ionization block (the
  // validated `cfmass` convention).  This is the fraction that does.
  float Jpsi_vtxsgnchk = 0.f;
  std::vector<float> resinfvtxv;  // a_{b,v} = dV_b^{1/2} w_v,b, 5/block padded
  std::vector<float> vtxvarv;     // v_b / sigma_v^2, aligned with reseigidx
  std::vector<float> vtxsgnv;     // the per-block ionization sign
  std::vector<short> vtxhitclsv;  // per-hit-class Gaussian shares of sigma_v^2
  std::vector<float> vtxhitvv;
  std::vector<float> cfvtxmsv, cfvtxdelv, cfvtxiorev, cfvtxioimv, cfvtxradrev, cfvtxradimv;
  std::vector<short> cfvtxgrpv;
  std::vector<float> cfvtxgrpmsv, cfvtxgrpdelv, cfvtxgrpiorev, cfvtxgrpioimv, cfvtxgrpradrev,
      cfvtxgrpradimv;
  std::vector<float> cfvtxgrpvqmsv, cfvtxgrpvqiov;
  float cfvtxgrpclosure = 0.f;
  // d theta_6^unconstrained / d(global params), aligned with `globalidxv`
  // exactly as `Jpsi_jacMass` is: the D row a DATA fit needs for the mean
  // (alignment / field) part of this term.
  std::vector<float> Jpsi_jacVtx;
  // The FITTED common vertex's 3x3 covariance (state indices 7,8,9), cm^2,
  // packed (xx, xy, xz, yy, yz, zz).  Exported in BOTH regimes and with the
  // beam rows on or off, because the leave-one-out beam residual is built
  // from it: the rows-OFF run's `Jpsi_x/y/z` and this matrix are what the ON
  // export must reproduce (gate G7).
  std::array<float, 6> Jpsi_covvtx = {{0.f, 0.f, 0.f, 0.f, 0.f, 0.f}};

  // ======================================================================
  // THE BEAM-LINE (LUMINOUS-REGION) CONSTRAINT, AND ITS TWO RESIDUALS
  // ======================================================================
  //
  // The luminous region is a GAUSSIAN NOISE BLOCK, not a prior on a global
  // parameter: constraining the common vertex to the beam line with the
  // beam-width covariance is mathematically ONE EXTRA HIT shared by the two
  // legs -- three rows with Jacobian `Fbs = I` on the vertex-PCA position
  // (state indices 7,8,9), residual `dbs0 = x_v - b0` and weight
  // `covBS^-1`.  It is right for a PROMPT resonance (Z) and WRONG for a
  // non-prompt one (the J/psi's B-decay fraction), so nothing here is
  // enabled by default.
  //
  // TWO DEFECTS of the original code, both fixed (branch `beamline-260913`):
  //
  //  1. THE ROWS ENTERED TWICE.  The emission sat inside the per-leg
  //     `for (id...)` loop with no `id == 0` guard (the pointing constraint
  //     right below it HAS one), so the same three rows -- same residual,
  //     same Jacobian on the same three state indices, same weight -- were
  //     written once per leg and `chisq0val += bschisq` ran twice.  The
  //     effective luminous-region covariance was therefore HALVED and the
  //     beam chi2 double-counted, while `ndof` already counted +3.  Now
  //     emitted once, `nbscons = 3`.
  //
  //  2. THE RESIDUAL WAS THE LEG'S REFERENCE POINT, not the vertex.
  //     `twoTrackPca2cart` puts the two reference points at
  //     `x_v -+ (d/2) n_hat`, so with index 6 free `refFts.head<3>()` is NOT
  //     `x_v` and `Fbs = I` on indices 7,8,9 is not its Jacobian.  The
  //     MIDPOINT `(x_a + x_b)/2` IS `x_v` identically, with Jacobian exactly
  //     `I` on 7,8,9 and zero on index 6 and on the momenta, so the residual
  //     is built from it and the Jacobian is exact in BOTH regimes.
  //
  // `beamWidthScale_` multiplies the three widths (covariance by its
  // square).  It exists for the gate: at 1e6 the rows are weightless and the
  // fit must reproduce the rows-OFF fit on every export.
  //
  // THE BLOCK IS REGISTERED AS A RESOLUTION BLOCK, family 16, `dV = covBS`
  // (the same convention as the hit families 8/9, whose `dV` is the hit
  // covariance: `V_b -> e^k V_b` has `dV_b/dk = V_b`).  Registering it is
  // what keeps `sum_b |a_b|^2 == sigma^2` EXACT for every functional with
  // the rows on -- the mass (`resinfv`), the vertex (`resinfvtxv`) and the
  // two beam residuals below.  Family 16 is a FAMILY LABEL ONLY: it has no
  // `detidparms` entry and adds no column to the quadratic term, so
  // `resglobidx` carries the sentinel `kBeamSpotGlobIdx` and
  // `varianceFamilyWanted(16)` is never true in any shipped configuration.
  // `cvhcf::trackExponents` counts it into `vgauss` alongside 8/9 (it IS a
  // Gaussian block), so the reader's group closure
  // `sum_g (vqms + vqio) + vgauss/sigma^2 == 1` stays exact.
  //
  // THE BEAM-SPOT GLOBAL PARAMETERS NEED NO SEPARATE EXPORT.  A shift `delta`
  // of the centroid x0 (y0) moves the beam residual row by `-delta`, and any
  // functional is `theta = w^T n`, so
  //     d theta / d x0 = -w[bs row x] ,     d theta / d y0 = -w[bs row y]
  // -- its own influence weight on the beam block's x (y) row.  For the
  // slopes, moving `dxdz` by `ds` moves the same row by `-ds (z_v - z0)`, so
  //     d theta / d dxdz = -w[bs row x] (z_v - z0)
  // (the tilt also enters `covBS`'s off-diagonal, a second-order effect on
  // the MEAN and ignored here).  `Jpsi_bsmean*` below are exactly those
  // weights, and `Jpsi_bsvtx` carries `z_v` so the slope rows follow.  No new
  // global parameters are added to the quadratic term by this study.
  //
  // THE TWO RESIDUALS.  With the rows ON the beam rows are MEASUREMENT rows,
  // so their unbiased residual is the LEAVE-ONE-OUT (innovation) residual:
  // the vertex the fit finds WITHOUT the beam rows, minus the beam line.
  // With `A = C[7:10,7:10]` (the fitted vertex covariance, rows in),
  // `M = covBS - A = Cov(rho_B)` and `rho_B = x_v^fit - b0` the post-step
  // residual on those rows,
  //     e_B    = rho_B + A M^-1 rho_B         = x_v^{-B} - b0
  //     C_{-B} = A + A M^-1 A                 (the rows-OFF vertex covariance)
  //     Cov(e) = C_{-B} + covBS
  // (Sherman-Morrison.)  THE "+" FORM IS DELIBERATE.  The same quantities are
  // `covBS M^-1 rho_B` and `covBS M^-1 covBS`, which is how they were first
  // written -- and that form evaluates the WHOLE answer through `M^-1`.  `M`
  // is a DIFFERENCE of two nearly equal covariances (it degenerates as the
  // tracks stop constraining the vertex, `A -> covBS`) and `A` carries the
  // numerical error of a large sparse solve, so the product could lose the
  // leading order entirely and could return an INDEFINITE covariance -- which
  // is what produced the 5 sigma tail in the first pass of this study.  In the
  // "+" form the leading order is explicit, only the correction is amplified,
  // and `Cov(e) >= covBS` holds by construction.  `A M^-1` comes from an `LLT`
  // solve of `M`, not an explicit inverse.  `Jpsi_bscovlo` exports `C_{-B}`,
  // which the rows-OFF run's `Jpsi_covvtx` must reproduce (gate G7a).
  // Projecting to the transverse plane with the
  // CONDITIONAL-MEAN projector
  //     P = [ I_2 , -s ] ,   s = covBS[0:2,2] / covBS[2,2]
  // is exactly "the beam line evaluated at that vertex's z": `P covBS P^T`
  // is the conditional covariance `Sigma_{xy|z}`, so
  //     r_bs      = P e_B                              (2, cm)
  //     Cov(r_bs) = P Cov(e) P^T = P C_{-B} P^T + Sigma_{xy|z}
  // -- the task's two readings are the same formula.  z is LEFT OUT: its row
  // is weightless (sigmaZ ~ 3.5 cm against a ~100 um vertex error).
  //
  // THE TWO CF FUNCTIONALS ARE THE WHITENED PAIR, not the global x and y
  // components.  `Cov(r_bs) = L L^T` (lower Cholesky, order x then y) and
  // `z = L^-1 r_bs`, so `Cov(z) = I`: each pull has UNIT variance and the two
  // are UNCORRELATED.  That is what stops the product of the two terms from
  // over-counting -- with the global pair `Cov_xy` is not zero (correlation
  // ~0.2) and multiplying the two likelihoods treated them as independent,
  // which showed up as sandwich/quoted = 1.385 for the `bs` channel alone.
  // The price, stated: neither pull is the response to a SINGLE beam-spot
  // parameter any more, so `Jpsi_bslinv` (the packed `L^-1`) is exported and
  // the global-basis reading of any weight is one 2x2 multiply away.
  //
  // The influence weights follow from `rho = (I - F C F^T Vinv) n`:
  //     w_k = u_k - Vinv F C F^T u_k ,   u_k = G_k^T on the 3 beam rows
  // with `G = L^-1 P (I + A M^-1)` (the "+" form again), and then
  // `sum_b |dV_b^{1/2} w_k,b|^2 == 1` exactly (gate `Jpsi_bsvchk`).
  //
  // `dV_b^{1/2}` is CACHED once per candidate (`ressqrtdV` in the two-track
  // maker) and shared by the mass, vertex and both beam influence loops: it
  // is a property of the block, not of the functional, and recomputing its
  // eigendecomposition four times per block was most of the beam
  // functionals' cost.
  bool exportBsResidual_ = false;
  double beamWidthScale_ = 1.0;
  static constexpr unsigned int kBeamSpotGlobIdx = 0xffffffffu;
  static constexpr int kBeamSpotFamily = 16;
  // r_bs, cm, (x, y) -- the leave-one-out vertex minus the beam line at z_v
  std::array<float, 2> Jpsi_bsres = {{0.f, 0.f}};
  // Cov(r_bs), cm^2, packed (xx, xy, yy)
  std::array<float, 3> Jpsi_bscov = {{0.f, 0.f, 0.f}};
  // the whitened pulls, LOWER-CHOLESKY basis in the order (x, y):
  // Cov = L L^T with L lower triangular, z = L^-1 r_bs.  So z[0] is the x
  // pull and z[1] the y pull GIVEN x -- stated because the basis is a
  // choice, and the eigenbasis would mix the two.
  std::array<float, 2> Jpsi_bsz = {{0.f, 0.f}};
  float Jpsi_bschi2 = 0.f;     // r_bs^T Cov^-1 r_bs == rho_B^T M^-1 rho_B (2 dof)
  float Jpsi_bschi2fit = 0.f;  // the FITTED beam rows' chi2, rho_B^T covBS^-1 rho_B
  float Jpsi_bschi20 = 0.f;    // the same at the LINEARISATION point, dbs0^T covBS^-1 dbs0
  float Jpsi_bsvchk = 0.f;     // max_i |sum_b |a_{b,i}|^2 / Cov_ii - 1|
  bool Jpsi_bsok = false;
  // THE LEAVE-ONE-OUT VERTEX COVARIANCE `C_{-B} = A + A M^-1 A`, cm^2, packed
  // (xx, xy, xz, yy, yz, zz) -- the vertex covariance the SAME fit would have
  // WITHOUT the beam rows, predicted from the rows-ON fit.  The rows-OFF
  // run's `Jpsi_covvtx` must reproduce it (gate G7a), and
  // `Cov(e_B) = C_{-B} + covBS` is the well-conditioned "+" form of the
  // innovation covariance: manifestly >= covBS, so it can never come back
  // indefinite the way `covBS M^-1 covBS` could.
  std::array<float, 6> Jpsi_bscovlo = {{0.f, 0.f, 0.f, 0.f, 0.f, 0.f}};
  // the lower-Cholesky inverse of `Cov(r_bs)`, packed (00, 10, 11): the map
  // from the RAW transverse pair to the WHITENED pair, `z = L^-1 r`.  The two
  // CF functionals are the whitened pulls, so this is what takes an influence
  // weight (or a beam-spot mean response) back to the global x/y basis.
  std::array<float, 3> Jpsi_bslinv = {{0.f, 0.f, 0.f}};
  // smallest eigenvalue of `M = covBS - C_vtx(ON)` in units of covBS's own
  // scale: the conditioning of the leave-one-out step, 1 when the tracks say
  // nothing about the vertex and 0 when the beam rows are redundant.  A
  // DIAGNOSTIC for the tail, not a cut.
  float Jpsi_bsmeig = 0.f;
  // the fitted common vertex (cm) and the beam spot it is measured against,
  // so the slope response and any re-derivation offline need no extra record
  std::array<float, 3> Jpsi_bsvtx = {{0.f, 0.f, 0.f}};
  std::array<float, 3> Jpsi_bsspot = {{0.f, 0.f, 0.f}};
  std::array<float, 2> Jpsi_bsslope = {{0.f, 0.f}};   // (dxdz, dydz) as used
  std::array<float, 3> Jpsi_bswidth = {{0.f, 0.f, 0.f}};  // the SCALED widths, cm
  // THE MEAN TERM: -w[bs rows] for each functional, i.e. d(functional)/d(x0,y0,z0).
  // Units: GeV/cm for the mass, dimensionless for the vertex (cm/cm) and the
  // beam residuals.  Slope response = the x (y) entry times (z_v - z0).
  std::array<float, 3> Jpsi_bsmeanmass = {{0.f, 0.f, 0.f}};
  std::array<float, 3> Jpsi_bsmeanvtx = {{0.f, 0.f, 0.f}};
  // WHITENED pull 1 then whitened pull 2 (it was the global x then y
  // functional before the pair was whitened): `w_k` on the beam rows is
  // `(L^-1 P)_k^T`, so these are `-(L^-1 P)`.
  std::array<float, 6> Jpsi_bsmeanbs = {{0.f, 0.f, 0.f, 0.f, 0.f, 0.f}};
  // per-block influence for the TWO beam functionals (the WHITENED pulls),
  // COMPONENT MAJOR:
  // entry (k, b, j) at ((k * nblock) + b) * 5 + j, k = 0 (x) / 1 (y)
  std::vector<float> resinfbsv;
  // v_{b,i} / Cov_ii, same component-major layout without the 5-padding
  std::vector<float> bsvarv;
  // the family shares of Cov_ii, i = x then y: beam block, hits, MS, ionization
  std::array<float, 2> Jpsi_bsvbs = {{0.f, 0.f}};
  std::array<float, 2> Jpsi_bsvhit = {{0.f, 0.f}};
  std::array<float, 2> Jpsi_bsvms = {{0.f, 0.f}};
  std::array<float, 2> Jpsi_bsvioni = {{0.f, 0.f}};
  // the beam block's share of the MASS and VERTEX functionals' variance
  float Jpsi_massvbs = 0.f, Jpsi_vtxvbs = 0.f;
  // THE SAME SHARE SPLIT BY BEAM DIRECTION -- the derivative of the block's
  // variance share with respect to a scale on `sigma_x^2` (`...x`) and on
  // `sigma_y^2` (`...y`), each cross term split half and half, so
  // `...x + ...y + (the z part) == the total share` identically.  These are
  // what makes the luminous-region WIDTHS float in the CF terms exactly as a
  // hit class does: a LINEAR variance scale `(1 + eps)` on them.
  float Jpsi_massvbsx = 0.f, Jpsi_massvbsy = 0.f;
  float Jpsi_vtxvbsx = 0.f, Jpsi_vtxvbsy = 0.f;
  std::array<float, 2> Jpsi_bsvbsx = {{0.f, 0.f}};
  std::array<float, 2> Jpsi_bsvbsy = {{0.f, 0.f}};
  // the RECORD's errors on the two transverse widths, cm (BeamWidthXError /
  // BeamWidthYError as read) -- the prior width for the two scales
  std::array<float, 2> Jpsi_bswidtherr = {{0.f, 0.f}};
  // CF exponents at the two beam weights, component major [2 * kNTau]
  std::vector<float> cfbsmsv, cfbsdelv, cfbsiorev, cfbsioimv, cfbsradrev, cfbsradimv;
  std::vector<short> cfbshitclsv;   // (component, class) pairs: comp in cfbshitcompv
  std::vector<short> cfbshitcompv;
  std::vector<float> cfbshitvv;
  // per-material-group exponents at the two beam weights; the group key is
  // (component, group) via `cfbsgrpcompv` / `cfbsgrpv`
  std::vector<short> cfbsgrpv, cfbsgrpcompv;
  std::vector<float> cfbsgrpmsv, cfbsgrpdelv, cfbsgrpiorev, cfbsgrpioimv, cfbsgrpradrev,
      cfbsgrpradimv;
  std::vector<float> cfbsgrpvqmsv, cfbsgrpvqiov;
  std::array<float, 2> cfbsgrpclosure = {{0.f, 0.f}};
  std::array<float, 2> Jpsi_bssgnchk = {{0.f, 0.f}};

  // ======================================================================
  // THE PER-HIT (COMPLEMENT) RESIDUAL EXPORT -- the DATA version of the
  // hit-residual CF likelihood.  See `exportPerHitResidual_`.
  //
  // A truth-referenced form whitens `refParms - genParms` and needs MC.
  // What exists on data is the part of the constraint residual the fit has
  // NOT absorbed:
  //
  //     rho = V R r ,   R = V^-1 - V^-1 F C F^T V^-1 ,   Cov(rho) = V R V
  //
  // of rank `d = ncons - nstatefree = n_meas - 5`.  Restricted to the
  // MEASUREMENT rows it loses nothing (the kink rows of `rho` are a
  // deterministic function of them, and the Mahalanobis form is invariant
  // under a bijection of the support), so the export is built on the
  // `n_meas = nvalid + nvalidpixel` measurement rows in HIT ORDER, whitened
  // by the LDL^T of `G = V_mm - F_m C F_m^T` with the `n_meas - d` null
  // pivots skipped.  Then `Cov(z) = I_d` and `sum_k z_k^2 = r^T R r`, the
  // fit's own chi2 -- both are exported as gates.
  // ======================================================================
  int phresd = 0;         // d, the number of whitened components kept
  int phresnmeas = 0;     // n_meas = nvalid + nvalidpixel
  int phresnfree = 0;     // nstatefree, so `d == n_meas - 5` can be audited
  float phreschi2 = 0.f;  // sum_k z_k^2   (GATE 1: == chisqval)
  float phresvchk = 0.f;  // max_k |sum_b v^(k)_b - 1|   (GATE 2)
  // lambda_dexp / lambda_(dexp+1) of G, descending: the separation the
  // IMPOSED rank `dexp = ncons - nstatefree` rests on.  Large = the rank
  // decision is unambiguous; near 1 = it is not, and the row should be cut.
  float phresrankgap = 0.f;
  // max_k |(Psi^T G Psi)_kk - 1|: the whitener against the SAME G it was
  // built from.  Small here with `phres_vchk` large means the factorization
  // is right and the influence W is not, and vice versa.
  float phresgchk = 0.f;
  // rank(Fw) from the column-pivoted QR; equals `nstatefree` on a healthy
  // fit, and `d = ncons - phres_qrank` by construction.
  int phresqrank = 0;
  // number of TRUTH-REFERENCED components appended after the `phres_d`
  // per-hit ones (5 on gen-matched MC, 0 on data).  Every per-component array
  // -- `phresz`, `phresvarv`, the `phcf_*` exponents, the hit-class shares --
  // has `phres_d + phres_nref` entries, the last `nref` of which are the
  // Cholesky-whitened `refParms - genParms`; their `phreshit`/`phresdim`/
  // `phrescls` are -1 and their `phresrow` is -1-j.
  int phresnref = 0;
  bool phresok = false;
  // [d] the whitened complement residual, the quantity the likelihood eats
  std::vector<float> phresz;
  // [n_meas] the post-fit residual itself, in hit order (pre-whitening)
  std::vector<float> phresraw;
  // [d], all parallel: which measurement row led component k, that row's
  // valid-hit index, its local coordinate (0 = first, 1 = the pixel's
  // second) and its hit-resolution class (the `reshitcls` code).
  std::vector<short> phresrow;
  std::vector<short> phreshit;
  std::vector<short> phresdim;
  std::vector<short> phrescls;
  // [d] the LDL pivot and the conditioning `G_kk / pivot_k`.  A cut on the
  // latter is a cut on the FIT'S COVARIANCE, not on the residual, so it
  // cannot bias the distribution being measured.
  std::vector<float> phrespiv;
  std::vector<float> phresinflat;
  // [d * nres], COMPONENT MAJOR (`k*nres + b`): the block's variance share of
  // component k, `v^(k)_b = W[b-rows, k]^T dV_b W[b-rows, k]`, carrying the
  // SIGN of `W[r0, k]` (the qop-row influence) so the ionization/radiative
  // weight can be signed offline exactly as in the maker.  `sum_b |v|` over
  // the non-parmtype-15 blocks is 1 for every k.
  std::vector<float> phresvarv;
  // [nres * ntot * 5] BLOCK MAJOR then component then dof:
  // `A_b = W_b^T dV_b^{1/2}`, dof-padded to 5.  `|A_b[k]|^2` is
  // `|phresvarv|`, so this is redundant for the likelihood; it exists for the
  // fourth CROSS cumulant between components, which needs `A_b[j] . A_b[k]`.
  std::vector<float> phresbv;
  // The per-component CF exponents, same six families and same `cftau` grid
  // as the `cf*` block, laid out [d * kNTau].
  std::vector<float> phcfmsv, phcfdelv, phcfiorev, phcfioimv, phcfradrev, phcfradimv;
  std::vector<float> phcfvgf;  // [d] the parmtype-8/9 (Gaussian) share
  // The per-(component, material group) split: `phcfgrpcomp` and `phcfgrpv`
  // are the (k, group) key of each slot, the arrays are [nslot * kNTau].
  std::vector<short> phcfgrpcomp, phcfgrpv;
  std::vector<float> phcfgrpmsv, phcfgrpdelv, phcfgrpiorev, phcfgrpioimv, phcfgrpradrev, phcfgrpradimv;
  // per slot: the group's share of the component's variance under the FIT'S
  // OWN Q (Rossi `thp2` for MS, `ioniSq2` for ionization, nothing radiative
  // or delta -- the fit's Q has neither).  This is what the Gaussian hit-chi2
  // assumes, and it is NOT recoverable from the exponents, which carry the
  // model's full Moliere second moment (14 % larger).  Closure:
  // `sum_slots(vqms + vqio) + phcf_vgf == 1` for every component.
  std::vector<float> phcfgrpvqms, phcfgrpvqio;
  float phcfgrpclosure = 0.f;
  // Per-(component, hit class) Gaussian variance shares; `sum_c` over a
  // component is that component's `phcfvgf`.
  std::vector<short> phcfhitcomp, phcfhitcls;
  std::vector<float> phcfhitv;
  int phcfnok = 0;      // components whose cvhcf call returned ok
  float phcfms = 0.f;   // wall-clock ms spent in the d cvhcf calls
  // Runtree: the tau grid and the model provenance, written once per global
  // parameter (constant, so ROOT compresses them away) rather than per event.
  std::vector<float> cftau;
  std::string cfmodel;

  // Export switches.
  //   exportCfExponents_  -- compute and write the cf* branches at all.
  //   exportStepRecords_  -- write the RAW per-step export the exponents are
  //                          built from. It is the whole 430 kB/candidate;
  //                          with the exponents validated it is only needed
  //                          to re-derive them with a different model.
  bool exportCfExponents_ = true;
  bool exportStepRecords_ = false;
  //   exportCfGroupExponents_ -- additionally split those exponents by
  //                          material group (see cfgrp*). +26 kB/candidate,
  //                          so it is opt-in and off by default.
  bool exportCfGroupExponents_ = false;
  //   exportPerHitResidual_ -- build and write the per-hit (complement)
  //                          residual block above.  OFF by default: it costs
  //                          `d` extra `cvhcf` evaluations per track.
  bool exportPerHitResidual_ = false;
  //   perHitCfGroups_    -- also split the per-component exponents by
  //                          material group.  On by default WHEN the block is
  //                          on: without it the material amounts cannot be
  //                          floated, which is the whole point.
  bool perHitCfGroups_ = true;
  //   perHitRefComponents_ -- append the 5 truth-referenced reference-state
  //                          components (`refParms - genParms` whitened by
  //                          `refCov`) to the same arrays, so one file
  //                          carries both terms in one convention.  MC only;
  //                          silently 0 components without a gen match.
  bool perHitRefComponents_ = true;
  //   perHitInfluenceBlocks_ -- also write `phresbv`, the per-(block,
  //                          component) influence VECTORS. ~23 kB/track, and
  //                          the only consumer is the cross-cumulant sizing
  //                          of the composite-likelihood approximation.
  bool perHitInfluenceBlocks_ = true;
  //   perHitShareMin_    -- zero a block's share of a component when it is
  //                          below this FRACTION of that component's unit
  //                          variance, so `cvhcf` skips it.  A controlled
  //                          approximation (the model variance moves by at
  //                          most the dropped share), 0 = keep everything,
  //                          which is the default and what the gates run at.
  double perHitShareMin_ = 0.;
  // Does this maker's functional use the delta-recoil family? The q/p one
  // does; the mass one's offline reference (`build_pairs_tt`) does not, so
  // the two-track maker does not pay for a sixth per-group array.
  bool cfGroupDelta_ = true;
  //   exportHitResBlocks_ -- register the parmtype-8/9 HIT-RESOLUTION dV
  //                          blocks in the two-track maker, which is what
  //                          makes the per-hit-class resolution parameters
  //                          fittable from two-track candidates (the
  //                          single-track maker always registers them).
  //                          It is EXPORT ONLY -- the two-track `dVs` feed
  //                          nothing but the influence export -- so it
  //                          cannot move the fit.
  bool exportHitResBlocks_ = true;
  //   exportMaterialNoise_ -- register the parmtype-15 MATERIAL-GROUP process
  //                          noise as a resolution family, so the quadratic
  //                          term's gradient and Hessian for k_g carry the
  //                          group's WIDTH as well as its mean loss. It
  //                          CHANGES the exported G and H of the parmtype-15
  //                          columns (not the track fit, which does not read
  //                          `dVs`), so it is opt-in; with False the group
  //                          amount enters only through the mean loss.
  bool exportMaterialNoise_ = false;

  //   exportVarianceGrads_ -- add the VARIANCE (log-det) part of the profiled
  //                          -2lnL to the exported global gradient and
  //                          Hessian of the TWO-TRACK maker; it is the
  //                          two-track counterpart of the single-track
  //                          maker's `gradll`.  Opt-in: with it off `k_g`
  //                          (and every parmtype-8..11 family) enters the
  //                          exported gradient and Hessian only through the
  //                          MEAN loss.
  //
  //                          The objective differentiated is the REML/marginal
  //                          one, the same one the single-track maker uses:
  //                              -2lnL = r^T R r + ln|V| + ln|C|,
  //                              R = V^-1 - V^-1 F C^-1 F^T V^-1,
  //                              C = F^T V^-1 F,
  //                          so the local track parameters are integrated out
  //                          rather than merely profiled, and the implicit
  //                          derivative through the fitted state vanishes by
  //                          the envelope theorem.
  bool exportVarianceGrads_ = false;
  //   varianceGradFamilies_ -- which parmtypes get the log-det treatment.
  //                          EMPTY means {8, 9, 10, 11, 15}.  {15} is the
  //                          LAYOUT-PRESERVING subset: the material-group
  //                          global indices are already columns of
  //                          `globalidxv` (one slot per group per hit), so
  //                          enabling only 15 adds no parameter and the
  //                          output can still be pooled with a production
  //                          that ran without the switch.  Families 8-11 are
  //                          per-module and are NOT columns of the two-track
  //                          parameter vector, so enabling them APPENDS
  //                          columns (and grows `jacrefv` / `Jpsi_jacMass` /
  //                          `hessfactorv` with them).
  std::vector<unsigned int> varianceGradFamilies_;
  // Is `fam` (a parmtype) one of the families whose log-det derivative is
  // exported?  Empty `varianceGradFamilies_` means the full default set.
  bool varianceFamilyWanted(int fam) const {
    if (!exportVarianceGrads_ || fam < 0) {
      return false;
    }
    if (varianceGradFamilies_.empty()) {
      return fam == 8 || fam == 9 || fam == 10 || fam == 11 || fam == 15;
    }
    for (unsigned int f : varianceGradFamilies_) {
      if (int(f) == fam) {
        return true;
      }
    }
    return false;
  }
  //   exportObjective_ -- write `objval`, the value of the marginal objective
  //                          above in DOUBLE precision.  Debug/validation
  //                          only: it costs an ncons x ncons LDLT per
  //                          candidate and exists so that the new gradient
  //                          can be finite-differenced against the thing it
  //                          claims to be the derivative of.
  bool exportObjective_ = false;
  //   varianceFDGlobalIdx_ / varianceFDEps_ -- IN-MAKER finite difference of
  //                          the variance gradient, at FIXED linearization.
  //                          `>= 0` does that one global index, `-2` does
  //                          every enabled variance column of families
  //                          10/11/15 (whose dV blocks are the 5x5 process
  //                          noise, the only ones whose V block is separately
  //                          invertible in its own row range).  It perturbs
  //                          V -> V + s dV_i with r, F and J held fixed and
  //                          re-does the profile, so it tests the ASSEMBLY --
  //                          the traces, the projector, the sign, the ln|C|
  //                          term -- to O(s^2), which the propagator-level FD
  //                          cannot do because there the reference trajectory
  //                          moves with the parameter.  Prints `VARFD ...`.
  int varianceFDGlobalIdx_ = -1;
  double varianceFDEps_ = 1e-3;

  // THE CANONICAL 18 HIT CLASSES, the C++ image of
  // `calibration_studies/resolution/hitres_classes.py:class_of`:
  //   0-3   pix_x_q0..q3    (pixel local x, by template charge bin)
  //   4-7   pix_y_q0..q3    (pixel local y)
  //   8-17  str_N1_lo, str_N1_hi, ... str_N5_lo, str_N5_hi
  //                         (strips, by cluster width N clipped to 1..5 and
  //                          uProj below/above 0.25)
  // The ORDER is part of the file format: a relabelling would silently
  // rename every hit of every file, so it is fixed here and in the python
  // by the same enumeration and must not be reordered.
  static constexpr int kNHitResClasses = 18;
  static int hitResClassIndex(int subdet, int sizeX, float uProj, int qBin, bool isY);


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


// ---------------------------------------------------------------------------
// Gauss-Newton step control shared by the single-track, two-track and N-track
// makers.
//
// The damping is a RELATIVE trust region in q/p, on top of the absolute
// momentum floor that keeps a leg above the propagator's refusal threshold
// (p < PropagationPtotLimit). A bare absolute floor does not work on its own:
// scaling a step so that the leg lands exactly on the floor is a hard
// non-linearity at a fixed momentum, so every track whose TRUE momentum is
// below the floor is pinned at it, and a leg whose reference momentum is
// ALREADY below the floor gets a negative scale, which max(s,0) turns into a
// frozen (zero) step.
//
// With the trust region, per iteration a leg's momentum may change by at most
// a factor f = maxMomentumStepFactor (default 2, i.e. p may at most halve or
// double). The effective lower bound is
//     p_lo = max(absFloor, p_ref/f)   if p_ref > absFloor
//     p_lo = p_ref/f                  otherwise
// which is ALWAYS strictly below p_ref -- so the scale is never zero, no leg
// is ever pinned at a fixed momentum, and a leg already below the absolute
// floor can still climb out. The upper bound p_hi = p_ref*f catches the
// opposite runaway (a stiff track pulled toward p -> inf / across q/p = 0).
namespace cvhstep {

  // Largest s in (0,1] such that qop_ref + s*dqop respects the momentum window
  // described above for ONE leg. Callers take the minimum over legs and scale
  // the whole (coupled) step vector by it, so the step direction is preserved.
  //
  //   absFloor     absolute momentum floor [GeV]; must stay above the
  //                propagator's PropagationPtotLimit refusal.
  //   f            max per-iteration momentum change factor; f <= 1 disables
  //                the relative window.
  //   qopFlipAllow |q/p| at or below which a genuine charge flip is permitted
  //                (i.e. p_ref >= allowChargeFlipAboveP). Such a flip is let
  //                through subject only to the momentum floor on the far side.
  //   flipProtect  set to true when a NON-permitted flip was capped (the
  //                nChargeFlipProtect bookkeeping).
  //
  // At the default f = 2 the cap on a non-permitted flip stops the step
  // half-way toward q/p = 0.
  inline double legStepScaleRel(double qopref, double dqop, double absFloor,
                                double f, double qopFlipAllow,
                                bool *flipProtect = nullptr) {
    if (qopref == 0. || dqop == 0. || !(f > 1.) || !std::isfinite(dqop)) {
      return 1.;
    }
    const double pref = std::abs(1. / qopref);
    double pLo = pref / f;
    if (pref > absFloor) {
      pLo = std::max(pLo, absFloor);
    }
    const double qopHi = 1. / pLo;               // |q/p| may not exceed this
    const double qopLo = std::abs(qopref) / f;   // |q/p| may not fall below this
    const double qopupd = qopref + dqop;
    const bool flip = qopupd * qopref <= 0.;
    const bool flipAllowed = std::abs(qopref) <= qopFlipAllow;
    double s = 1.;
    if (flip && flipAllowed) {
      // permitted (ambiguous-charge) flip: only the momentum floor applies
      if (std::abs(qopupd) > qopHi) {
        s = (std::copysign(qopHi, qopref) - qopref) / dqop;
      }
    } else if (flip) {
      // divergence-like flip: treat as the extreme upward step and cap it
      s = (std::copysign(qopLo, qopref) - qopref) / dqop;
      if (flipProtect != nullptr) {
        *flipProtect = true;
      }
    } else if (std::abs(qopupd) > qopHi) {
      s = (std::copysign(qopHi, qopref) - qopref) / dqop;
    } else if (std::abs(qopupd) < qopLo) {
      s = (std::copysign(qopLo, qopref) - qopref) / dqop;
    }
    if (!(s < 1.)) {
      return 1.;
    }
    // pLo < p_ref and qopLo < |qopref| strictly, so s > 0 by construction;
    // the clamp below only defends against a non-finite dqop.
    return std::isfinite(s) ? std::max(s, 0.) : 0.;
  }

}  // namespace cvhstep

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
