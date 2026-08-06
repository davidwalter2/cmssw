// CandidateVertexGeometryProducer
//
// Turns a fitted decay vertex into PV/BS-relative geometry. Does NO fitting.
//
// Why this exists: the CVH refit
// code emits only what comes out of the fit. The joint N-body maker has no
// pointing constraint, so it never reads a primary vertex -- giving it a PV
// token would put event content the fit does not use inside the fit code. This
// module is the stepping stone: it consumes the maker's vertex + covariance
// ValueMaps together with the PV and beam-spot collections, and produces the
// pointing and flight-significance variables the analysis selection needs.
//
// It is deliberately generic over the vertex source: it knows nothing about
// how many tracks went into the vertex, so the 2-body, 3-body and future
// 4-track / displaced-daughter configurations all reuse it unchanged.
//
// NO PREFILTER (proposal, 2026-07-28). Every selection lives at histmaker
// level. Outputs are 1:1 with the candidate collection and pre-filled with
// -99; a candidate whose upstream fit failed (fitOk = 0, or a sentinel vertex)
// keeps the sentinel rather than being dropped.
//
// Outputs (all ValueMap<float>, keyed to srcCandidates):
//   alphaBS  TRANSVERSE angle between the fitted momentum and the BS->vertex
//            flight direction. Transverse, not 3D: the beam spot constrains
//            x-y tightly and z not at all.
//   lxy      transverse flight distance from the beam spot
//   sxy      lxy / sigma(lxy)
//   l3d      3D flight distance from the closest-z primary vertex
//   sl3d     l3d / sigma(l3d)
//   pvIdx    index of the primary vertex used (-1 if none)
//
// NO FIT happens here (D27b): no DOCA, no midpoint, no iteration, no
// propagation. Every output is a dot product, a norm, an acos, or a covariance
// projected onto the flight direction. Vertex-FINDING happened upstream.
// Consequence: only a fitter producing a genuine fitted vertex may feed this.
// The two-track CVH maker with doVtxConstraint=False produces a MIDPOINT that
// neither track passes through -- it must not be wired in here.

#include <cctype>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"

namespace {
constexpr float kSentinel = -99.f;

// Symmetric 3x3, row-major upper storage: xx xy xz yy yz zz.
struct Sym3 {
  double m[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
  double project(const double u[3]) const {
    double s = 0.;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) s += u[i] * m[i][j] * u[j];
    return s;
  }
};
}  // namespace

class CandidateVertexGeometryProducer : public edm::global::EDProducer<> {
 public:
  explicit CandidateVertexGeometryProducer(const edm::ParameterSet&);
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

 private:
  edm::EDGetTokenT<edm::View<reco::Candidate>> candToken_;
  edm::EDGetTokenT<reco::BeamSpot> bsToken_;
  edm::EDGetTokenT<reco::VertexCollection> pvToken_;
  // vertex position + covariance, as emitted by the maker
  edm::EDGetTokenT<edm::ValueMap<float>> vxT_, vyT_, vzT_;
  edm::EDGetTokenT<edm::ValueMap<float>> cxxT_, cxyT_, cxzT_, cyyT_, cyzT_, czzT_;
  // FITTED momentum, part of the vertex block. NOT candidate.p4(),
  // which is the raw four-vector sum -- a pointing angle taken from that would
  // use uncorrected kinematics.
  edm::EDGetTokenT<edm::ValueMap<float>> mptT_, metaT_, mphiT_;
  edm::EDGetTokenT<edm::ValueMap<int>> okT_;
  const bool haveOk_;

  edm::EDPutTokenT<edm::ValueMap<float>> oAlphaBS_, oLxy_, oSxy_, oL3d_, oSl3d_;
  edm::EDPutTokenT<edm::ValueMap<int>> oPvIdx_;
};

CandidateVertexGeometryProducer::CandidateVertexGeometryProducer(
    const edm::ParameterSet& cfg)
    : candToken_(consumes<edm::View<reco::Candidate>>(
          cfg.getParameter<edm::InputTag>("srcCandidates"))),
      bsToken_(consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamSpot"))),
      pvToken_(consumes<reco::VertexCollection>(
          cfg.getParameter<edm::InputTag>("primaryVertices"))),
      haveOk_(!cfg.getParameter<edm::InputTag>("fitOk").label().empty()) {
  const std::string src = cfg.getParameter<edm::InputTag>("vertexSrc").label();
  // Instance-name prefix for the vertex block. Empty for a module that emits
  // one block (the N-body maker); "raw"/"ref" for JpsiXKinematicFitProducer,
  // which emits one per arm. Capitalisation follows the emitter: "" -> vtxX,
  // "raw" -> rawVtxX.
  const std::string pf = cfg.getParameter<std::string>("blockPrefix");
  auto cap = [&](const char* inst) {
    std::string n(inst);
    if (!pf.empty()) n[0] = static_cast<char>(std::toupper(n[0]));
    return pf + n;
  };
  auto vm = [&](const char* inst) {
    return consumes<edm::ValueMap<float>>(edm::InputTag(src, cap(inst)));
  };
  vxT_ = vm("vtxX");
  vyT_ = vm("vtxY");
  vzT_ = vm("vtxZ");
  cxxT_ = vm("vtxCovXX");
  cxyT_ = vm("vtxCovXY");
  cxzT_ = vm("vtxCovXZ");
  cyyT_ = vm("vtxCovYY");
  cyzT_ = vm("vtxCovYZ");
  czzT_ = vm("vtxCovZZ");
  mptT_ = vm("motherPt");
  metaT_ = vm("motherEta");
  mphiT_ = vm("motherPhi");
  if (haveOk_)
    okT_ = consumes<edm::ValueMap<int>>(cfg.getParameter<edm::InputTag>("fitOk"));

  oAlphaBS_ = produces<edm::ValueMap<float>>("alphaBS");
  oLxy_ = produces<edm::ValueMap<float>>("lxy");
  oSxy_ = produces<edm::ValueMap<float>>("sxy");
  oL3d_ = produces<edm::ValueMap<float>>("l3d");
  oSl3d_ = produces<edm::ValueMap<float>>("sl3d");
  oPvIdx_ = produces<edm::ValueMap<int>>("pvIdx");
}

void CandidateVertexGeometryProducer::produce(edm::StreamID, edm::Event& iEvent,
                                              const edm::EventSetup&) const {
  edm::Handle<edm::View<reco::Candidate>> candH;
  iEvent.getByToken(candToken_, candH);
  const std::size_t n = candH.isValid() ? candH->size() : 0;

  std::vector<float> alphaBS(n, kSentinel), lxy(n, kSentinel), sxy(n, kSentinel),
      l3d(n, kSentinel), sl3d(n, kSentinel);
  std::vector<int> pvIdx(n, -1);

  edm::Handle<reco::BeamSpot> bsH;
  iEvent.getByToken(bsToken_, bsH);
  edm::Handle<reco::VertexCollection> pvH;
  iEvent.getByToken(pvToken_, pvH);

  auto get = [&](const edm::EDGetTokenT<edm::ValueMap<float>>& t) {
    edm::Handle<edm::ValueMap<float>> h;
    iEvent.getByToken(t, h);
    return h;
  };
  const auto vx = get(vxT_), vy = get(vyT_), vz = get(vzT_);
  const auto cxx = get(cxxT_), cxy = get(cxyT_), cxz = get(cxzT_);
  const auto cyy = get(cyyT_), cyz = get(cyzT_), czz = get(czzT_);
  const auto mpt = get(mptT_), meta = get(metaT_), mphi = get(mphiT_);
  if (!(mpt.isValid() && meta.isValid() && mphi.isValid()))
    throw cms::Exception("VertexBlockIncomplete")
        << "CandidateVertexGeometryProducer requires the FITTED momentum "
           "(motherPt/motherEta/motherPhi) from the vertex block. Falling back "
           "to candidate.p4() would silently use the raw four-vector sum.";
  edm::Handle<edm::ValueMap<int>> okH;
  if (haveOk_) iEvent.getByToken(okT_, okH);

  const bool haveVtx = vx.isValid() && vy.isValid() && vz.isValid();
  const bool haveCov = cxx.isValid() && cxy.isValid() && cxz.isValid() &&
                       cyy.isValid() && cyz.isValid() && czz.isValid();

  for (std::size_t i = 0; haveVtx && i < n; ++i) {
    const edm::Ptr<reco::Candidate> ptr = candH->ptrAt(i);
    if (haveOk_ && okH.isValid() && (*okH)[ptr] != 1) continue;  // sentinel stays
    const double X = (*vx)[ptr], Y = (*vy)[ptr], Z = (*vz)[ptr];
    if (X == kSentinel && Y == kSentinel && Z == kSentinel) continue;

    Sym3 C;
    if (haveCov) {
      C.m[0][0] = (*cxx)[ptr];
      C.m[0][1] = C.m[1][0] = (*cxy)[ptr];
      C.m[0][2] = C.m[2][0] = (*cxz)[ptr];
      C.m[1][1] = (*cyy)[ptr];
      C.m[1][2] = C.m[2][1] = (*cyz)[ptr];
      C.m[2][2] = (*czz)[ptr];
    }

    // FITTED momentum from the block -- see the token comment above.
    const double pt = (*mpt)[ptr], eta = (*meta)[ptr], phi = (*mphi)[ptr];
    if (pt <= 0.) continue;
    const double px = pt * std::cos(phi), py = pt * std::sin(phi);
    const double pz = pt * std::sinh(eta);
    const double pmag = std::sqrt(px * px + py * py + pz * pz);
    if (pmag <= 0.) continue;

    // ---- transverse, relative to the beam spot -----------------------------
    if (bsH.isValid()) {
      // The BS is evaluated at the vertex z, so its x-y slope is respected.
      const double bx = bsH->x(Z), by = bsH->y(Z);
      const double dx = X - bx, dy = Y - by;
      const double L = std::sqrt(dx * dx + dy * dy);
      lxy[i] = L;
      if (L > 0.) {
        // Direction of flight; sigma along it, including the BS width. The BS
        // contribution is added, not neglected: it is comparable to the vertex
        // error in x-y for a well-measured B.
        const double u[3] = {dx / L, dy / L, 0.};
        double var = C.project(u);
        var += u[0] * u[0] * bsH->BeamWidthX() * bsH->BeamWidthX() +
               u[1] * u[1] * bsH->BeamWidthY() * bsH->BeamWidthY();
        if (var > 0.) sxy[i] = L / std::sqrt(var);
      }
      // TRANSVERSE pointing angle -- x-y only, matching route 1's
      // `dimuonAlphaBS` ("XY pointing angle wrt BS").
      //
      // It must NOT be done in 3D from the beam spot: the primary vertex is
      // spread over +-8 cm in z along the beam while the B flies ~0.05 cm, so
      // a 3D flight vector taken from the BS centre is dominated by
      // (z_PV - z_BS) and carries essentially no information about the B
      // direction. (Measured: 3D-from-BS gives a median angle of 2.20 rad,
      // i.e. 126 deg, which is nonsense for a pointing variable.) The beam
      // spot constrains x-y tightly and z not at all, so the transverse
      // projection is the only well-defined choice against it.
      const double ptmag = std::sqrt(px * px + py * py);
      if (ptmag > 0. && L > 0.) {
        double ct = (dx * px + dy * py) / (L * ptmag);
        ct = std::max(-1., std::min(1., ct));
        alphaBS[i] = std::acos(ct);
      }
    }

    // ---- 3D, relative to the closest-z primary vertex ----------------------
    if (pvH.isValid() && !pvH->empty()) {
      int best = -1;
      double bestdz = 1e9;
      for (std::size_t k = 0; k < pvH->size(); ++k) {
        const double dz = std::abs((*pvH)[k].z() - Z);
        if (dz < bestdz) {
          bestdz = dz;
          best = static_cast<int>(k);
        }
      }
      if (best >= 0) {
        pvIdx[i] = best;
        const reco::Vertex& pv = (*pvH)[best];
        const double dx = X - pv.x(), dy = Y - pv.y(), dz = Z - pv.z();
        const double L = std::sqrt(dx * dx + dy * dy + dz * dz);
        l3d[i] = L;
        if (L > 0.) {
          const double u[3] = {dx / L, dy / L, dz / L};
          double var = C.project(u);
          // PV error along the same direction -- the two vertices are
          // independent, so the variances add.
          for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) var += u[a] * pv.covariance(a, b) * u[b];
          if (var > 0.) sl3d[i] = L / std::sqrt(var);
        }
      }
    }
  }

  auto putF = [&](const edm::EDPutTokenT<edm::ValueMap<float>>& tok,
                  const std::vector<float>& v) {
    edm::ValueMap<float> m;
    if (candH.isValid()) {
      edm::ValueMap<float>::Filler f(m);
      f.insert(candH, v.begin(), v.end());
      f.fill();
    }
    iEvent.emplace(tok, std::move(m));
  };
  putF(oAlphaBS_, alphaBS);
  putF(oLxy_, lxy);
  putF(oSxy_, sxy);
  putF(oL3d_, l3d);
  putF(oSl3d_, sl3d);
  {
    edm::ValueMap<int> m;
    if (candH.isValid()) {
      edm::ValueMap<int>::Filler f(m);
      f.insert(candH, pvIdx.begin(), pvIdx.end());
      f.fill();
    }
    iEvent.emplace(oPvIdx_, std::move(m));
  }
}

DEFINE_FWK_MODULE(CandidateVertexGeometryProducer);
