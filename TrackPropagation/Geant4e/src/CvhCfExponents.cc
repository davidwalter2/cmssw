#include "TrackPropagation/Geant4e/interface/CvhCfExponents.h"

#include "TrackPropagation/Geant4e/interface/CGFQoPBlock.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <map>
#include <mutex>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

// THE ARITHMETIC IS THE REFERENCE'S.  The knock-on quadrature is the
// difference of two Filon sums that cancel to ~1e-6 on a thin step, so the
// last bits of every product show at the scale the validation measures.  The
// release compiles with -march=x86-64-v3, where GCC contracts a*b + c into a
// fused multiply-add; numpy in its baseline arithmetic and the standalone
// validation build do not.  Contraction is off for this file (not for the
// propagator, whose numerics it would move), so the maker, the shim and the
// reference form every product and sum alike.
#pragma GCC optimize("fp-contract=off")

namespace cvhcf {

  //==========================================================================
  // 0. THE TAU GRID
  //==========================================================================
  namespace {
    struct TauGrid {
      std::array<double, kNTau> t;
      TauGrid() {
        // The offline grid is np.linspace(0, 14, 448); numpy computes
        // linspace as `start + step*i` with step = (stop-start)/(N-1) and
        // then FORCES the last sample to `stop`. Only the stride-4 subset is
        // exported and its last element (index 252) is not the endpoint, so
        // the plain `i*step` form reproduces the reference exactly.
        const double step = kTauFullMax / static_cast<double>(kTauFullN - 1);
        for (int i = 0; i < kNTau; ++i)
          t[i] = static_cast<double>(i * kTauStride) * step;
      }
    };
    const TauGrid &tauGridObj() {
      static const TauGrid g;
      return g;
    }
  }  // namespace

  const double *tauGrid() { return tauGridObj().t.data(); }

  void Exponents::clear() {
    ms.fill(0.);
    ioRe.fill(0.);
    ioIm.fill(0.);
    radRe.fill(0.);
    radIm.fill(0.);
    kxRe.fill(0.);
    kxIm.fill(0.);
    kjRe.fill(0.);
    kjIm.fill(0.);
  }

  //==========================================================================
  // 1. THE SHAPE TABLES  (cf_ms_exact._GE)
  //==========================================================================
  namespace {

    // The shape tables' id, recorded in `modelTag`: the file name stem.
    const std::string kShapeTableId = "cvhcf_gshape_elec_v2";

    struct ShapeTables {
      int nTau = 0, nY = 0;
      std::vector<double> gtau, elecY, hard, dipole;
      // cached logs of the Y axis, for the row interpolation
      double lnY0 = 0., dlnY = 0.;
      std::string path;
      bool loaded = false;
    };

    ShapeTables &tables() {
      static ShapeTables t;
      return t;
    }
    std::once_flag &tablesOnce() {
      static std::once_flag f;
      return f;
    }

    void readTables(const std::string &pathIn) {
      ShapeTables &T = tables();
      std::string path = pathIn;
      if (path.empty())
        path = edm::FileInPath("TrackPropagation/Geant4e/data/" + kShapeTableId + ".bin").fullPath();
      std::ifstream f(path, std::ios::binary);
      if (!f)
        throw cms::Exception("CvhCfExponents") << "cannot open the Moliere shape table " << path;
      char magic[9] = {0};
      f.read(magic, 8);
      if (std::strncmp(magic, "CVHCFGSH", 8) != 0)
        throw cms::Exception("CvhCfExponents") << path << " is not a cvhcf shape table";
      int hdr[3] = {0, 0, 0};
      f.read(reinterpret_cast<char *>(hdr), sizeof(hdr));
      if (hdr[0] != 2)
        throw cms::Exception("CvhCfExponents") << path << ": unsupported table version " << hdr[0];
      T.nTau = hdr[1];
      T.nY = hdr[2];
      if (T.nTau < 2 || T.nY < 2)
        throw cms::Exception("CvhCfExponents") << path << ": nonsensical table dimensions";
      auto rd = [&](std::vector<double> &v, std::size_t n) {
        v.resize(n);
        f.read(reinterpret_cast<char *>(v.data()), static_cast<std::streamsize>(n * sizeof(double)));
        if (!f)
          throw cms::Exception("CvhCfExponents") << path << ": truncated table";
      };
      rd(T.gtau, T.nTau);
      rd(T.elecY, T.nY);
      rd(T.hard, static_cast<std::size_t>(T.nY) * T.nTau);
      rd(T.dipole, static_cast<std::size_t>(T.nY) * T.nTau);
      if (f.peek() != std::char_traits<char>::eof())
        throw cms::Exception("CvhCfExponents") << path << ": trailing bytes after the shape tables";
      T.lnY0 = std::log(T.elecY[0]);
      T.dlnY = std::log(T.elecY[1]) - std::log(T.elecY[0]);
      T.path = path;
      T.loaded = true;
    }

  }  // namespace

  void loadShapeTables(const std::string &path) {
    // `call_once` with an explicit path is a first-caller-wins contract: the
    // tables are process-global (they are the model), so two different tables
    // in one job would silently mean two different models.
    std::call_once(tablesOnce(), [&]() { readTables(path); });
    if (!path.empty() && tables().path != path)
      throw cms::Exception("CvhCfExponents")
          << "the cvhcf shape tables are already loaded from " << tables().path << "; " << path << " was requested";
  }

  namespace {
    inline const ShapeTables &T() {
      std::call_once(tablesOnce(), [&]() { readTables(std::string()); });
      return tables();
    }

    // --- np.interp on the (log-spaced but NOT assumed so) tau axis --------
    // `left = nan` in the reference is immediately overwritten by the tiny-tau
    // quadratic branch, so it is implemented directly; `right = row[-1]`.
    inline double interpTau(const ShapeTables &t, const double *row, double tau) {
      tau = std::fabs(tau);
      const double t0 = t.gtau[0];
      if (tau < t0) {
        const double r = tau / t0;
        return row[0] * r * r;
      }
      if (tau >= t.gtau[t.nTau - 1])
        return row[t.nTau - 1];
      // binary search: largest j with gtau[j] <= tau
      int lo = 0, hi = t.nTau - 1;
      while (hi - lo > 1) {
        const int mid = (lo + hi) >> 1;
        if (t.gtau[mid] <= tau)
          lo = mid;
        else
          hi = mid;
      }
      const double slope = (row[lo + 1] - row[lo]) / (t.gtau[lo + 1] - t.gtau[lo]);
      return slope * (tau - t.gtau[lo]) + row[lo];
    }

    // --- gshape_elec's row: linear in log(y2max) between two table rows ---
    void elecRow(const ShapeTables &t, double y2max, bool dipole, std::vector<double> &row) {
      double ly = std::log(y2max);
      ly = std::max(t.lnY0, std::min(ly, std::log(t.elecY[t.nY - 1])));
      const double fi = (ly - t.lnY0) / t.dlnY;
      int i0 = static_cast<int>(std::floor(fi));
      i0 = std::max(0, std::min(i0, t.nY - 2));
      const double fr = fi - i0;
      const std::vector<double> &tab = dipole ? t.dipole : t.hard;
      const double *r0 = tab.data() + static_cast<std::size_t>(i0) * t.nTau;
      const double *r1 = r0 + t.nTau;
      row.resize(t.nTau);
      for (int i = 0; i < t.nTau; ++i)
        row[i] = (1. - fr) * r0[i] + fr * r1[i];
    }

  }  // namespace

  //==========================================================================
  // 2. MOLIERE PARAMETERS  (cf_ms_exact.moliere_params, the default switches)
  //==========================================================================
  namespace {
    constexpr double kAlphaEm = 1.0 / 137.036;
    constexpr double kMeGeV = 0.51099895000e-3;
    // MS_CHI0_G4 = True: chi_0 = alpha m_e / 0.88534
    const double kChi0G4 = kAlphaEm * 0.51099895e-3 / 0.88534;
    constexpr double kG4ScreenF = 1.0;    // cf_ms_exact.G4_SCREEN_F
    constexpr double kFfConstnMeV2 = 6.937e-6;  // MS_FF_G4 = True

    struct MolPars {
      double chic2, chia2, thff2;
    };

    MolPars moliereParams(double effZ, double effA, double xg, double pGeV, double beta,
                          double zzp1OverA, double lnScreenW, bool haveSums) {
      MolPars p{0., 0., 0.};
      if (haveSums && zzp1OverA > 0.)
        p.chic2 = 0.157e-6 * zzp1OverA * xg / (pGeV * pGeV * beta * beta);
      else
        p.chic2 = 0.157e-6 * effZ * (effZ + 1.) / effA * xg / (pGeV * pGeV * beta * beta);
      const double az = kAlphaEm * effZ / beta;
      const double chi0 = kChi0G4 * std::pow(effZ, 1. / 3.) / pGeV;
      if (haveSums && std::isfinite(lnScreenW) && lnScreenW != 0.) {
        const double c = kChi0G4 / pGeV;
        p.chia2 = c * c * std::exp(lnScreenW);
      } else {
        const double g4screen = 1. + kG4ScreenF * std::exp(-effZ * effZ * 1.0e-3);
        p.chia2 = chi0 * chi0 * (1.13 + 3.76 * az * az) * g4screen;
      }
      const double pmev = pGeV * 1.0e3;
      p.thff2 = 2. / (kFfConstnMeV2 * std::pow(std::max(effA, 1.), 0.54) * pmev * pmev);
      return p;
    }

    // numpy's np.round(x, 3): multiply by 1000, rint (ties to even), divide.
    inline double round3(double x) { return std::rint(x * 1000.0) / 1000.0; }

    inline double clipd(double x, double lo, double hi) { return std::max(lo, std::min(x, hi)); }
  }  // namespace

  //==========================================================================
  // 3. THE MULTIPLE-SCATTERING FAMILY
  //==========================================================================
  namespace {

    // Per-step quantities the MS and delta-ray terms share.
    struct MsStep {
      double effZ = 0., effA = 0., xg = 0., pGeV = 0., beta = 0., thp2 = 0.;
      double chic2 = 0., chia2 = 0., thff2 = 0.;
      double lgRow = 0.;  // round3(log(clip(theta_FF/chi_a, 10, 1e7)))
      double lgEle = 0.;  // round3(log(clip(theta_e,max/chi_a, 10, 1e7)))
      double fN = 0., fE = 0.;
    };

    // One record's step, or false where the reference's `ok` (thp2 > 0) and
    // `act` (chi_c^2, chi_a^2 > 0) masks drop it.  `xgFrac` is an entry's
    // share of the record's material (cf_rows.ms_rows scales xg by it before
    // the Moliere parameters are formed, as here).
    bool msStepOf(const double *r, int stride, double xgFrac, const RowConfig &cfg, MsStep &s) {
      const bool haveSums = stride >= 10;
      if (!(r[5] > 0.))
        return false;
      s = MsStep();
      s.effZ = r[0];
      s.effA = r[1];
      s.xg = (xgFrac != 1.0) ? r[2] * xgFrac : r[2];
      s.pGeV = r[3];
      s.beta = r[4];
      s.thp2 = r[5];
      {
        const MolPars mp =
            moliereParams(s.effZ, s.effA, s.xg, s.pGeV, s.beta, haveSums ? r[7] : 0., haveSums ? r[8] : 0., haveSums);
        if (!(mp.chic2 > 0.) || !(mp.chia2 > 0.))
          return false;
        s.chic2 = mp.chic2;
        s.chia2 = mp.chia2;
        s.thff2 = mp.thff2;
        const double ym = std::sqrt(s.thff2 / s.chia2);
        s.lgRow = round3(std::log(clipd(ym, 1e1, 1e7)));
        // The atomic-electron ceiling, G4WentzelOKandVIxSection's own:
        // 1 - cos(theta_e,max) = Tmax m_e / p^2, theta^2 = 2 (1 - cos).
        // The mass comes from the record as m = p sqrt(1-beta^2)/beta --
        // `msmoliv` carries no PDG code, and the reference does the same.
        const double bt = clipd(s.beta, 1e-9, 1. - 1e-15);
        const double gam = 1. / std::sqrt(1. - bt * bt);
        const double mgev = s.pGeV / (bt * gam);
        const double bg = bt * gam;
        const double rat = kMeGeV / mgev;
        double tmx = 2. * kMeGeV * bg * bg / (1. + 2. * gam * rat + rat * rat);
        // Under the knock-on joint law the collisions above the e- production
        // threshold are the knock-on channel's, with their exact two-body
        // kinematics, and the electron term stops at that threshold, as
        // Geant4's msc does (cf_track_resolution._knockon_carries_delta).
        if (cfg.knockonJoint)
          tmx = std::min(tmx, cfg.knockonTcut * 1e-3);
        const double te = 2. * tmx * kMeGeV / (s.pGeV * s.pGeV);
        const double yme = std::min(std::sqrt(te / s.chia2), ym);
        // the electron ceiling reaches down to the production threshold's
        // angle, a few chi_a (the shape table starts at y = 10^-1/2)
        s.lgEle = round3(std::log(clipd(yme, std::pow(10., -0.5), 1e7)));
        s.fN = s.effZ / (s.effZ + 1.);
        s.fE = 1. / (s.effZ + 1.);
      }
      return true;
    }

    // The steps of `n` whole records (the vectorised path of the reference).
    void buildMsSteps(const double *rows, int stride, int n, const RowConfig &cfg, std::vector<MsStep> &out) {
      out.clear();
      out.reserve(static_cast<std::size_t>(n));
      MsStep s;
      for (int i = 0; i < n; ++i)
        if (msStepOf(rows + static_cast<std::size_t>(i) * stride, stride, 1.0, cfg, s))
          out.push_back(s);
    }

    // `ms_step_exponent` at MS_ELEC_TMAX = MS_ELEC_EDGE = MS_FINE_G = 1 and
    // MS_SNAP_YMAX = MS_WVI_SPLIT = 0. ADDS into S.
    //
    // The GROUPING and the ORDER of accumulation follow the reference exactly:
    // ascending rounded ln(ymax) for the nuclear piece, and inside each such
    // group ascending rounded ln(yme) for the electron piece. The grouping is
    // not an optimization the reference happens to make -- the rounding to
    // 1e-3 in ln CHANGES the ceiling that `gshape_elec` is called with, so it
    // is part of the model and has to be reproduced, not improved on.
    // `wstd` is now a LIST of `nk` weights and `S` is `nk * nt`, laid out
    // functional-major: the concatenated argument list { wstd[k] tau[j] }.
    // The per-row work above the k loop -- the `gshape_elec` row, the step
    // grouping, `sqrt(chi_a^2)` -- is the part that is now done once instead
    // of nk times; the argument itself is formed by the SAME expression, in
    // the same association, so each k is bitwise the single-weight answer.
    void msExponentImpl(
        const std::vector<MsStep> &st, const double *wstd, int nk, const double *tau, int nt, double *S) {
      if (st.empty() || nk <= 0)
        return;
      const ShapeTables &tb = T();
      std::vector<double> rows;
      rows.reserve(st.size());
      for (const MsStep &s : st)
        rows.push_back(s.lgRow);
      std::sort(rows.begin(), rows.end());
      rows.erase(std::unique(rows.begin(), rows.end()), rows.end());

      std::vector<double> row, lyes;
      for (double r : rows) {
        const double ymr = std::exp(r);
        elecRow(tb, ymr * ymr, /*dipole=*/true, row);
        for (const MsStep &s : st) {
          if (s.lgRow != r)
            continue;
          const double sq = std::sqrt(s.chia2);
          const double w = s.chic2 * s.fN / s.chia2;
          for (int k = 0; k < nk; ++k) {
            const double wk = wstd[k];
            double *Sk = S + static_cast<std::size_t>(k) * nt;
            for (int j = 0; j < nt; ++j)
              Sk[j] += w * interpTau(tb, row.data(), sq * (wk * tau[j]));
          }
        }
        lyes.clear();
        for (const MsStep &s : st)
          if (s.lgRow == r)
            lyes.push_back(s.lgEle);
        std::sort(lyes.begin(), lyes.end());
        lyes.erase(std::unique(lyes.begin(), lyes.end()), lyes.end());
        for (double v : lyes) {
          const double ymv = std::exp(v);
          elecRow(tb, ymv * ymv, /*dipole=*/false, row);
          for (const MsStep &s : st) {
            if (s.lgRow != r || s.lgEle != v)
              continue;
            const double sq = std::sqrt(s.chia2);
            const double w = s.chic2 * s.fE / s.chia2;
            for (int k = 0; k < nk; ++k) {
              const double wk = wstd[k];
              double *Sk = S + static_cast<std::size_t>(k) * nt;
              for (int j = 0; j < nt; ++j)
                Sk[j] += w * interpTau(tb, row.data(), sq * (wk * tau[j]));
            }
          }
        }
      }
    }

  }  // namespace

  //==========================================================================
  // 4. THE IONISATION BLOCK VARIANCE
  //==========================================================================
  double ioniSq2(const float *rows, int stride, int n, const float *qsc, int nqsc) {
    if (rows == nullptr || n <= 0)
      return 0.;
    // `steps[:,1] * gq * gq`, left to right: the association is the
    // reference's, so the no-scale path is bit-identical to the caches.
    std::vector<double> w2(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
      const float *r = rows + static_cast<std::size_t>(i) * stride;
      const double gq = static_cast<double>(r[10]) * 1e-3;
      w2[i] = static_cast<double>(r[1]) * gq * gq;
    }
    double plain = 0.;
    for (double x : w2)
      plain += x;
    if (qsc == nullptr || nqsc <= 0)
      return plain;
    // Several legs can share one global index (a track crossing the same
    // module twice, ~1/3 of blocks) and each carries its own scale, so the
    // pooled sum is sum_l sc_l * (step sum of leg l).
    long long tot = 0;
    for (int i = 0; i < nqsc; ++i)
      tot += static_cast<long long>(qsc[2 * i + 1]);
    if (tot != static_cast<long long>(n)) {
      double msc = 0.;
      for (int i = 0; i < nqsc; ++i)
        msc += qsc[2 * i];
      return plain * (msc / nqsc);
    }
    double s = 0.;
    int off = 0;
    for (int i = 0; i < nqsc; ++i) {
      const int ns = static_cast<int>(qsc[2 * i + 1]);
      double leg = 0.;
      for (int k = 0; k < ns; ++k)
        leg += w2[off + k];
      s += static_cast<double>(qsc[2 * i]) * leg;
      off += ns;
    }
    return s;
  }


  //==========================================================================
  // 5. THE NUCLEAR-ELASTIC FAMILY  (ksclosure/nucel/nucel_tables.Table.rows,
  //     through which cf_nucel_exact.rows_exponent and
  //     ks_nucel_cf.step_family evaluate; layout and rules in
  //     data/make_cvhcf_nucel_tables.py)
  //==========================================================================
  namespace {

    const std::string kNucelTableId = "cvhcf_nucel_v2";
    // elements below this share of a node's collisions are left out of the
    // material's joint bins (trace elements; the families keep them) --
    // nucel_tables.JSHARE
    constexpr double kNucelJShare = 1e-5;

    struct NucelTable {
      int nS = 0, nE = 0, nU = 0, nD = 0, nR = 0;
      std::vector<double> ug, dee, prate, lpr, binC, elemA;
      std::vector<int> pdg, elemZ;
      double lug0 = 0., dlug = 0.;
      struct Species {
        std::vector<double> edges, nodeP;
        std::vector<int> seg;
        // node range [segA[k], segB[k]] of model segment k (0 .. nEdge)
        std::vector<int> segA, segB;
      };
      std::vector<Species> sp;
      std::vector<double> lnMu;              // [s][e][r]
      std::vector<std::size_t> recBase;      // first record of species s
      std::vector<float> fixed;              // [rec][7 + nU]
      std::vector<int> b0, nb;               // nonzero recoil-bin range per record
      std::vector<std::size_t> binOff;       // offset into `bins` per record
      std::vector<float> bins;               // per record: P[nb], D1[nb], D2[nb]
      std::vector<int> nj;                   // joint bins per record
      std::vector<std::size_t> jOff;         // offset into `jbins` per record
      std::vector<float> jbins;              // per record: theta[nj], dE[nj], w[nj]
      std::string path;

      std::size_t rec(int s, int e, int n) const {
        return recBase[s] + static_cast<std::size_t>(e) * sp[s].nodeP.size() + n;
      }
      int speciesIndex(int code) const {
        for (int s = 0; s < nS; ++s)
          if (pdg[s] == code)
            return s;
        return -1;
      }
      int elementIndex(int z) const {
        for (int e = 0; e < nE; ++e)
          if (elemZ[e] == z)
            return e;
        return -1;
      }
    };

    NucelTable &nucTable() {
      static NucelTable t;
      return t;
    }
    std::once_flag &nucOnce() {
      static std::once_flag f;
      return f;
    }

    void readNucel(const std::string &pathIn) {
      NucelTable &T = nucTable();
      std::string path = pathIn;
      if (path.empty())
        path = edm::FileInPath("TrackPropagation/Geant4e/data/" + kNucelTableId + ".bin").fullPath();
      std::ifstream f(path, std::ios::binary);
      if (!f)
        throw cms::Exception("CvhCfExponents") << "cannot open the nuclear-elastic table " << path;
      auto rdraw = [&](void *dst, std::size_t nbytes) {
        f.read(reinterpret_cast<char *>(dst), static_cast<std::streamsize>(nbytes));
        if (!f)
          throw cms::Exception("CvhCfExponents") << path << ": truncated nuclear-elastic table";
      };
      auto rdd = [&](std::vector<double> &v, std::size_t n) {
        v.resize(n);
        rdraw(v.data(), n * sizeof(double));
      };
      auto rdi = [&](std::vector<int> &v, std::size_t n) {
        v.resize(n);
        rdraw(v.data(), n * sizeof(int));
      };
      static_assert(sizeof(int) == 4 && sizeof(float) == 4 && sizeof(double) == 8, "table layout");
      char magic[8];
      rdraw(magic, 8);
      if (std::strncmp(magic, "CVHNUCEL", 8) != 0)
        throw cms::Exception("CvhCfExponents") << path << " is not a cvhcf nuclear-elastic table";
      int hdr[6];
      rdraw(hdr, sizeof(hdr));
      if (hdr[0] != 2)
        throw cms::Exception("CvhCfExponents") << path << ": unsupported nuclear-elastic table version " << hdr[0];
      T.nS = hdr[1];
      T.nE = hdr[2];
      T.nU = hdr[3];
      T.nD = hdr[4];
      T.nR = hdr[5];
      if (T.nS < 1 || T.nE < 1 || T.nU < 4 || T.nD < 1 || T.nR < 2)
        throw cms::Exception("CvhCfExponents") << path << ": nonsensical nuclear-elastic table dimensions";
      rdd(T.ug, T.nU);
      rdd(T.dee, T.nD + 1);
      rdd(T.prate, T.nR);
      rdi(T.pdg, T.nS);
      rdi(T.elemZ, T.nE);
      rdd(T.elemA, T.nE);
      T.sp.resize(T.nS);
      std::size_t nrec = 0;
      T.recBase.resize(T.nS);
      for (int s = 0; s < T.nS; ++s) {
        int nn[2];
        rdraw(nn, sizeof(nn));
        NucelTable::Species &S = T.sp[s];
        rdd(S.edges, nn[0]);
        rdd(S.nodeP, nn[1]);
        rdi(S.seg, nn[1]);
        if (nn[1] < 1)
          throw cms::Exception("CvhCfExponents") << path << ": species " << T.pdg[s] << " has no kernel nodes";
        S.segA.assign(nn[0] + 1, -1);
        S.segB.assign(nn[0] + 1, -1);
        for (int n = 0; n < nn[1]; ++n) {
          const int k = S.seg[n];
          if (k < 0 || k > nn[0] || (n > 0 && (k < S.seg[n - 1] || !(S.nodeP[n] > S.nodeP[n - 1]))))
            throw cms::Exception("CvhCfExponents") << path << ": species " << T.pdg[s] << ": bad node segments";
          if (S.segA[k] < 0)
            S.segA[k] = n;
          S.segB[k] = n;
        }
        for (int k = 0; k <= nn[0]; ++k)
          if (S.segA[k] < 0)
            throw cms::Exception("CvhCfExponents")
                << path << ": species " << T.pdg[s] << ": model segment " << k << " has no kernel node";
        T.recBase[s] = nrec;
        nrec += static_cast<std::size_t>(T.nE) * nn[1];
      }
      rdd(T.lnMu, static_cast<std::size_t>(T.nS) * T.nE * T.nR);
      T.fixed.resize(nrec * (7 + T.nU));
      rdraw(T.fixed.data(), T.fixed.size() * sizeof(float));
      std::vector<int> idx;
      rdi(idx, 2 * nrec);
      T.b0.resize(nrec);
      T.nb.resize(nrec);
      T.binOff.resize(nrec);
      std::size_t tot = 0;
      for (std::size_t r = 0; r < nrec; ++r) {
        T.b0[r] = idx[2 * r];
        T.nb[r] = idx[2 * r + 1];
        if (T.b0[r] < 0 || T.nb[r] < 0 || T.b0[r] + T.nb[r] > T.nD)
          throw cms::Exception("CvhCfExponents") << path << ": recoil bin range out of bounds";
        T.binOff[r] = tot;
        tot += 3 * static_cast<std::size_t>(T.nb[r]);
      }
      T.bins.resize(tot);
      rdraw(T.bins.data(), tot * sizeof(float));
      rdi(T.nj, nrec);
      T.jOff.resize(nrec);
      std::size_t jtot = 0;
      for (std::size_t r = 0; r < nrec; ++r) {
        if (T.nj[r] < 0)
          throw cms::Exception("CvhCfExponents") << path << ": negative joint bin count";
        T.jOff[r] = jtot;
        jtot += 3 * static_cast<std::size_t>(T.nj[r]);
      }
      T.jbins.resize(jtot);
      rdraw(T.jbins.data(), jtot * sizeof(float));
      char tail[8];
      rdraw(tail, 8);
      if (std::strncmp(tail, "CVHNUEND", 8) != 0 || f.peek() != std::char_traits<char>::eof())
        throw cms::Exception("CvhCfExponents") << path << ": bad nuclear-elastic table trailer / size";
      T.lpr.resize(T.nR);
      for (int r = 0; r < T.nR; ++r)
        T.lpr[r] = std::log(T.prate[r]);
      T.binC.resize(T.nD);
      for (int b = 0; b < T.nD; ++b)
        T.binC[b] = 0.5 * (T.dee[b] + T.dee[b + 1]);
      T.lug0 = std::log(T.ug[0]);
      T.dlug = (std::log(T.ug[T.nU - 1]) - T.lug0) / (T.nU - 1);
      T.path = path;
    }

    inline const NucelTable &NT() {
      std::call_once(nucOnce(), [&]() { readNucel(std::string()); });
      return nucTable();
    }

    // log-log interpolation of a ln(mu) row on the rate grid: linear in ln p
    // between the bracketing nodes, constant beyond the ends.
    inline double logLog(const NucelTable &t, const double *lmu, double lp) {
      int j = static_cast<int>(std::lower_bound(t.lpr.begin(), t.lpr.end(), lp) - t.lpr.begin()) - 1;
      j = std::max(0, std::min(j, t.nR - 2));
      const double f = clipd((lp - t.lpr[j]) / (t.lpr[j + 1] - t.lpr[j]), 0., 1.);
      return (1. - f) * lmu[j] + f * lmu[j + 1];
    }

  }  // namespace

  // The mixture record of one (species, material): what `Table.mixture`
  // holds, with the recoil bins reduced to their uniform-bin form.
  struct NucelMix {
    int isp = -1;
    std::vector<double> lmu;   // ln mu_mat on the rate grid
    std::vector<double> mom;   // [node][7]: <theta>, <theta^2>, <dE>, <dE^2>, P_low, M1_low, M2_low
    std::vector<double> gm1;   // [node][nU]: g - 1
    std::vector<int> binStart; // [node + 1] into bP / bMu / bHw
    std::vector<double> bP, bMu, bHw;  // P_b, the bin mean, sqrt(3) x the bin's std
    // the joint (theta, dE) bins of every node: its elements' bins, each
    // weight times the element's collision share at the node
    std::vector<int> jStart;   // [node + 1] into jTh / jDe / jW
    std::vector<double> jTh, jDe, jW;
  };

  void NucelExponents::clear() {
    ang.fill(0.);
    recRe.fill(0.);
    recIm.fill(0.);
    jntRe.fill(0.);
    jntIm.fill(0.);
    N = 0.;
  }

  int nucelPdg(const std::string &n) {
    if (n == "pi+")
      return 211;
    if (n == "pi-")
      return -211;
    if (n == "kaon+")
      return 321;
    if (n == "kaon-")
      return -321;
    if (n == "proton")
      return 2212;
    if (n == "anti_proton")
      return -2212;
    return 0;
  }

  void loadNucelTables(const std::string &path) {
    std::call_once(nucOnce(), [&]() { readNucel(path); });
    if (!path.empty() && nucTable().path != path)
      throw cms::Exception("CvhCfExponents")
          << "the nuclear-elastic table is already loaded from " << nucTable().path << "; " << path << " was requested";
  }

  const std::string &nucelTableId() { return kNucelTableId; }

  NucelMixtures::NucelMixtures(NucelResolver resolver) : resolver_(std::move(resolver)) {}
  NucelMixtures::~NucelMixtures() = default;

  const NucelMix &NucelMixtures::get(int pdgCode, int matIndex) {
    const auto key = std::make_pair(pdgCode, matIndex);
    auto it = cache_.find(key);
    if (it != cache_.end())
      return *it->second;
    const NucelTable &t = NT();
    const int isp = t.speciesIndex(pdgCode);
    if (isp < 0)
      throw cms::Exception("CvhCfExponents") << "no nuclear-elastic table for species " << pdgCode;
    NucelComposition c;
    if (matIndex < 0 || !resolver_ || !resolver_(matIndex, c) || c.Z.size() != c.A.size() ||
        c.Z.size() != c.W.size())
      throw cms::Exception("CvhCfExponents") << "nuclear-elastic family: material index " << matIndex
                                             << " does not resolve to a composition";
    // elements with a mass share (Z >= 1, w > 0), as `Table.composition`
    std::vector<int> ei;
    std::vector<double> scale;
    for (std::size_t i = 0; i < c.Z.size(); ++i) {
      if (!(c.Z[i] >= 1.) || !(c.W[i] > 0.))
        continue;
      const int z = static_cast<int>(std::lround(c.Z[i]));
      const int e = t.elementIndex(z);
      if (e < 0)
        throw cms::Exception("CvhCfExponents") << "nuclear-elastic family: material index " << matIndex
                                               << " has element Z=" << z << ", which the table "
                                               << t.path << " does not carry";
      if (!(c.A[i] > 0.))
        throw cms::Exception("CvhCfExponents")
            << "nuclear-elastic family: material index " << matIndex << " has a non-positive atomic mass";
      ei.push_back(e);
      // w_i A_tab / A_i: at a fixed per-atom cross section the mass
      // attenuation scales as 1/A
      scale.push_back(c.W[i] * t.elemA[e] / c.A[i]);
    }
    if (ei.empty())
      throw cms::Exception("CvhCfExponents") << "nuclear-elastic family: material index " << matIndex
                                             << " has no element";
    const std::size_t nc = ei.size();
    const NucelTable::Species &S = t.sp[isp];
    const int nn = static_cast<int>(S.nodeP.size());
    auto m = std::make_unique<NucelMix>();
    m->isp = isp;
    // the rate on its own grid
    m->lmu.resize(t.nR);
    for (int r = 0; r < t.nR; ++r) {
      double s = 0.;
      for (std::size_t i = 0; i < nc; ++i)
        s += scale[i] * std::exp(t.lnMu[(static_cast<std::size_t>(isp) * t.nE + ei[i]) * t.nR + r]);
      m->lmu[r] = std::log(s);
    }
    // the kernel records at every node: rate-weighted mixtures of the
    // element records (linear in every stored entry, the recoil bins
    // included)
    const int nfix = 7 + t.nU;
    m->mom.assign(static_cast<std::size_t>(nn) * 7, 0.);
    m->gm1.assign(static_cast<std::size_t>(nn) * t.nU, 0.);
    m->binStart.assign(nn + 1, 0);
    m->jStart.assign(nn + 1, 0);
    std::vector<double> share(nc), P(t.nD), D1(t.nD), D2(t.nD);
    for (int n = 0; n < nn; ++n) {
      const double lp = std::log(S.nodeP[n]);
      double rs = 0.;
      for (std::size_t i = 0; i < nc; ++i) {
        share[i] = scale[i] * std::exp(logLog(t, &t.lnMu[(static_cast<std::size_t>(isp) * t.nE + ei[i]) * t.nR], lp));
        rs += share[i];
      }
      for (std::size_t i = 0; i < nc; ++i)
        share[i] /= rs;
      std::fill(P.begin(), P.end(), 0.);
      std::fill(D1.begin(), D1.end(), 0.);
      std::fill(D2.begin(), D2.end(), 0.);
      double *mo = &m->mom[static_cast<std::size_t>(n) * 7];
      double *gm = &m->gm1[static_cast<std::size_t>(n) * t.nU];
      for (std::size_t i = 0; i < nc; ++i) {
        const std::size_t r = t.rec(isp, ei[i], n);
        const float *fx = &t.fixed[r * nfix];
        for (int k = 0; k < 7; ++k)
          mo[k] += share[i] * static_cast<double>(fx[k]);
        for (int k = 0; k < t.nU; ++k)
          gm[k] += share[i] * static_cast<double>(fx[7 + k]);
        const int b0 = t.b0[r], nb = t.nb[r];
        const float *bb = &t.bins[t.binOff[r]];
        for (int b = 0; b < nb; ++b) {
          P[b0 + b] += share[i] * static_cast<double>(bb[b]);
          D1[b0 + b] += share[i] * static_cast<double>(bb[nb + b]);
          D2[b0 + b] += share[i] * static_cast<double>(bb[2 * nb + b]);
        }
      }
      // each bin uniform over its own mean -+ sqrt(3) x its own std
      for (int b = 0; b < t.nD; ++b) {
        if (!(P[b] > 0.))
          continue;
        const double d = D1[b] / P[b];
        m->bP.push_back(P[b]);
        m->bMu.push_back(t.binC[b] + d);
        m->bHw.push_back(std::sqrt(3. * std::max(D2[b] / P[b] - d * d, 0.)));
      }
      m->binStart[n + 1] = static_cast<int>(m->bP.size());
      for (std::size_t i = 0; i < nc; ++i) {
        if (share[i] < kNucelJShare)
          continue;
        const std::size_t r = t.rec(isp, ei[i], n);
        const int nj = t.nj[r];
        const float *jb = &t.jbins[t.jOff[r]];
        for (int k = 0; k < nj; ++k) {
          m->jTh.push_back(static_cast<double>(jb[k]));
          m->jDe.push_back(static_cast<double>(jb[nj + k]));
          m->jW.push_back(static_cast<double>(jb[2 * nj + k]) * share[i]);
        }
      }
      m->jStart[n + 1] = static_cast<int>(m->jTh.size());
    }
    const NucelMix &ref = *m;
    cache_.emplace(key, std::move(m));
    return ref;
  }

  namespace {

    // The node pair bracketing p inside p's model segment and the ln-p weight
    // of the upper one (`Table.nodes`); the end node alone outside the
    // segment's nodes.
    inline void nucelNodes(const NucelTable::Species &S, double p, int &j0, int &j1, double &f) {
      const int seg = static_cast<int>(std::upper_bound(S.edges.begin(), S.edges.end(), p) - S.edges.begin());
      const int a = S.segA[seg], b = S.segB[seg];
      const double *P = S.nodeP.data();
      if (p <= P[a]) {
        j0 = j1 = a;
        f = 0.;
        return;
      }
      if (p >= P[b]) {
        j0 = j1 = b;
        f = 0.;
        return;
      }
      int j = static_cast<int>(std::lower_bound(P + a, P + b + 1, p) - P) - 1;
      j = std::max(a, std::min(j, b));
      const int jj = std::min(j + 1, b);
      j0 = j;
      j1 = jj;
      f = (jj > j) ? clipd(std::log(p / P[j]) / std::log(P[jj] / P[j]), 0., 1.) : 0.;
    }

    // g - 1 of one node record at u >= 0 (`Table.gm1_node`).
    inline double nucelGm1(const NucelTable &t, const double *gm, double theta2, double u) {
      if (u < t.ug[0])
        return -0.25 * u * u * theta2;
      if (u > t.ug[t.nU - 1])
        return gm[t.nU - 1];
      const double tt = (std::log(u) - t.lug0) / t.dlug;
      int i = static_cast<int>(std::floor(tt));
      i = std::max(0, std::min(i, t.nU - 2));
      const double s = tt - i;
      const double p1 = gm[i], p2 = gm[i + 1];
      const double p0 = (i == 0) ? (gm[0] * 2 - gm[1]) : gm[i - 1];
      const double p3 = (i == t.nU - 2) ? (gm[t.nU - 1] * 2 - gm[t.nU - 2]) : gm[i + 2];
      return 0.5 * (2 * p1 + (p2 - p0) * s + (2 * p0 - 5 * p1 + 4 * p2 - p3) * s * s +
                    (3 * p1 - p0 - 3 * p2 + p3) * s * s * s);
    }

    // sin(y)/y - 1 without the cancellation near 0 (`nucel_tables.sincm1`):
    // the Taylor series, Horner in y^2 and complete to double precision,
    // below |y| = 1; the direct form above.
    inline double sincm1(double y) {
      static constexpr double c[9] = {-0.16666666666666666,
                                      0.008333333333333333,
                                      -0.0001984126984126984,
                                      2.7557319223985893e-06,
                                      -2.505210838544172e-08,
                                      1.6059043836821613e-10,
                                      -7.647163731819816e-13,
                                      2.8114572543455206e-15,
                                      -8.22063524662433e-18};
      if (!(std::fabs(y) < 1.0))
        return std::sin(y) / y - 1.0;
      const double y2 = y * y;
      double p = c[8];
      for (int k = 7; k >= 0; --k)
        p = p * y2 + c[k];
      return p * y2;
    }

    // ADDS c * (h(v_j) - 1) of node n at v_j = wq tau_j (`Table.hm1_node`),
    // formed without cancellation and normalised exactly (P_low + sum_b P_b
    // = 1 drops the constant):
    //   h - 1 = i v M1_low - v^2 M2_low / 2 + sum_b P_b [e^{i v mu_b} sinc(v a_b) - 1],
    //   e^{ix} s - 1 = -2 sin^2(x/2) s + (s - 1) + i sin(x) s.
    void nucelRecoilNodeAt(
        const NucelMix &m, int n, double wq, double c, const double *tau, int nt, double *re, double *im) {
      const double *mo = &m.mom[static_cast<std::size_t>(n) * 7];
      const int b0 = m.binStart[n], b1 = m.binStart[n + 1];
      for (int j = 0; j < nt; ++j) {
        const double v = wq * tau[j];
        double sr = 0., si = 0.;
        for (int b = b0; b < b1; ++b) {
          const double x = v * m.bMu[b], y = v * m.bHw[b];
          const double sc = (y != 0.) ? std::sin(y) / y : 1.;
          const double h = std::sin(0.5 * x);
          sr += m.bP[b] * ((-2.0 * h * h) * sc + sincm1(y));
          si += m.bP[b] * (std::sin(x) * sc);
        }
        re[j] += c * (-0.5 * v * v * mo[6] + sr);
        im[j] += c * (v * mo[5] + si);
      }
    }

    // The q/p-equivalent of a recoil energy loss T [MeV] at the step's E, p
    // [MeV] under the exact 1/p map (cf_knockon.t_eff).
    inline double nucelTEff(double T, double E, double p) {
      const double M2 = E * E - p * p;
      const double d = E - T;
      const double pp = std::sqrt(std::max(d * d - M2, 1e-300));
      return p * p * T * (2.0 * E - T) / (E * pp * (p + pp));
    }

    // ADDS c x the joint term of node n (`Table.joint_node`) at the arguments
    // b_j = bw tau_j (angular, >= 0) and a_j = aw tau_j (recoil, signed): over
    // the node's joint bins, the deflections rescaled by sc = p_node / p,
    //   sum_i W_i [(J0(b theta_i sc) - 1)(e^{i a X_i} - 1) + (e^{i a X_i} - e^{i a dE_i})],
    // X = T_eff(dE) at the row's (E, p) [GeV] under the exact map (the second
    // bracket is then the recoil's own map correction), else dE (where it
    // vanishes and is not formed).  The exact map ends where the knock-on
    // channel's does (`lawTop`): a recoil that would leave the primary with
    // less than pminFrac of its momentum keeps the linear map -- the upper
    // node's recoils near the kinematic end point can exceed the row's kinetic
    // energy, the recoil being carried at fixed momentum transfer.  Two
    // sin/cos pairs per bin and point: e^{iy} - 1 at the half angle, e^{iyx} -
    // e^{iyd} in sum-to-product form.
    void nucelJointNodeAt(const NucelMix &m,
                          int n,
                          double sc,
                          double bw,
                          double aw,
                          bool exact,
                          double pminFrac,
                          double E,
                          double p,
                          double c,
                          const double *tau,
                          int nt,
                          double *re,
                          double *im,
                          std::vector<double> &thp,
                          std::vector<double> &X) {
      const int k0 = m.jStart[n], nb = m.jStart[n + 1] - k0;
      if (nb <= 0)
        return;
      const double *de = &m.jDe[k0];
      const double *W = &m.jW[k0];
      thp.resize(nb);
      X.resize(nb);
      const double Em = 1e3 * E, pm = 1e3 * p, q = pminFrac * pm;
      const double tcap = Em - std::sqrt(Em * Em - pm * pm + q * q);
      for (int i = 0; i < nb; ++i) {
        thp[i] = m.jTh[k0 + i] * sc;
        X[i] = (exact && de[i] < tcap) ? nucelTEff(de[i], Em, pm) : de[i];
      }
      for (int j = 0; j < nt; ++j) {
        const double b = bw * tau[j], a = aw * tau[j];
        double sr = 0., si = 0.;
        for (int i = 0; i < nb; ++i) {
          const double jm1 = besselJ0(b * thp[i]) - 1.0;
          const double yx = a * X[i];
          // e^{iy} - 1 = (-2 s^2, 2 s c) at the half angle
          const double hs = std::sin(0.5 * yx), hc = std::cos(0.5 * yx);
          double r = jm1 * (-2.0 * hs * hs);
          double q = jm1 * (2.0 * hs * hc);
          if (exact) {
            // e^{iyx} - e^{iyd} by the sum-to-product forms
            const double yd = a * de[i];
            const double yp = 0.5 * (yx + yd), ym = std::sin(0.5 * (yx - yd));
            r = r + (-2.0 * std::sin(yp) * ym);
            q = q + (2.0 * std::cos(yp) * ym);
          }
          sr += r * W[i];
          si += q * W[i];
        }
        re[j] += c * sr;
        im[j] += c * si;
      }
    }

    // The species' mass [GeV] (cf_nucel_exact._SPECIES_MASS): the recoil's
    // exact map needs the step's energy.
    inline double nucelMassOf(int pdgCode) {
      switch (std::abs(pdgCode)) {
        case 211:
          return 0.13957039;
        case 321:
          return 0.493677;
        case 2212:
          return 0.93827209;
        default:
          break;
      }
      throw cms::Exception("CvhCfExponents") << "nuclear-elastic family: no mass for species " << pdgCode;
    }

    // THE ROW FUNCTION (`Table.rows`) over the rows `sel` (all n when null):
    // angular per entry, recoil per row (grouped over rows sharing a mixture,
    // a node and a weight, as the reference groups them), joint per row.
    // Returns sum_s N_s over the selected rows.
    template <typename Rec>
    double nucelRowsImpl(const double *tau,
                         int nt,
                         const Rec *rows,
                         int stride,
                         const int *sel,
                         int nsel,
                         const int *mat,
                         const int *pdgv,
                         const double *mass,
                         const int *rid,
                         const double *wb,
                         const double *frac,
                         int ne,
                         const double *wqRow,
                         const double *wbMid,
                         const RowConfig &cfg,
                         NucelMixtures &mix,
                         double *ang,
                         double *recRe,
                         double *recIm,
                         double *jntRe,
                         double *jntIm) {
      if (rows == nullptr || nsel <= 0 || stride < 4 || nt <= 0)
        return 0.;
      const NucelTable &t = NT();
      struct Row {
        const NucelMix *m = nullptr;
        double N = 0., p = 0.;
        int nd[2] = {0, 0};
        double w[2] = {0., 0.};
      };
      std::unordered_map<int, Row> live;
      std::vector<int> order;
      double Ntot = 0.;
      for (int q = 0; q < nsel; ++q) {
        const int s = sel ? sel[q] : q;
        if (pdgv[s] == 0)
          continue;
        const Rec *r = rows + static_cast<std::size_t>(s) * stride;
        const double xg = r[2], p = r[3];
        if (!(xg > 0.) || !(p > 0.))
          continue;
        const NucelMix &m = mix.get(pdgv[s], mat[s]);
        const double Ns = std::exp(logLog(t, m.lmu.data(), std::log(p))) * xg;
        Ntot += Ns;
        if (!(Ns > 0.))
          continue;
        Row &R = live[s];
        R.m = &m;
        R.N = Ns;
        R.p = p;
        double f;
        nucelNodes(t.sp[m.isp], p, R.nd[0], R.nd[1], f);
        R.w[0] = 1. - f;
        R.w[1] = f;
        order.push_back(s);
      }
      if (order.empty())
        return Ntot;
      // angular, per entry
      for (int e = 0; e < ne; ++e) {
        auto it = live.find(rid[e]);
        if (it == live.end() || wb[e] == 0.)
          continue;
        const Row &R = it->second;
        const double aw = std::fabs(wb[e]);
        for (int q = 0; q < 2; ++q) {
          if (R.w[q] == 0.)
            continue;
          const int n = R.nd[q];
          const double sc = t.sp[R.m->isp].nodeP[n] / R.p;  // fixed momentum transfer
          const double *gm = &R.m->gm1[static_cast<std::size_t>(n) * t.nU];
          const double th2 = R.m->mom[static_cast<std::size_t>(n) * 7 + 1];
          const double c = R.N * frac[e] * R.w[q];
          for (int j = 0; j < nt; ++j)
            ang[j] += c * nucelGm1(t, gm, th2, (aw * tau[j]) * sc);
        }
      }
      if (!cfg.nucelRecoil)
        return Ntot;
      // recoil, per row: once per (mixture, node, weight)
      struct Grp {
        const NucelMix *m;
        int n;
        double wq, c;
      };
      std::vector<Grp> grp;
      for (int s : order) {
        const Row &R = live[s];
        if (wqRow[s] == 0.)
          continue;
        for (int q = 0; q < 2; ++q) {
          if (R.w[q] == 0.)
            continue;
          auto it = std::find_if(grp.begin(), grp.end(), [&](const Grp &g) {
            return g.m == R.m && g.n == R.nd[q] && g.wq == wqRow[s];
          });
          if (it == grp.end())
            grp.push_back(Grp{R.m, R.nd[q], wqRow[s], R.N * R.w[q]});
          else
            it->c += R.N * R.w[q];
        }
      }
      for (const Grp &g : grp)
        nucelRecoilNodeAt(*g.m, g.n, g.wq, g.c, tau, nt, recRe, recIm);
      if (!cfg.nucelJoint)
        return Ntot;
      // joint, per row at (|wb_mid|, wq)
      std::vector<double> thp, X;
      for (int s : order) {
        const Row &R = live[s];
        if (wqRow[s] == 0.)
          continue;
        const double E = std::sqrt(R.p * R.p + mass[s] * mass[s]);
        const double bw = std::fabs(wbMid[s]);
        for (int q = 0; q < 2; ++q) {
          if (R.w[q] == 0.)
            continue;
          const int n = R.nd[q];
          nucelJointNodeAt(*R.m,
                           n,
                           t.sp[R.m->isp].nodeP[n] / R.p,
                           bw,
                           wqRow[s],
                           cfg.qopExact,
                           cfg.pminFrac,
                           E,
                           R.p,
                           R.N * R.w[q],
                           tau,
                           nt,
                           jntRe,
                           jntIm,
                           thp,
                           X);
        }
      }
      return Ntot;
    }

  }  // namespace

  double nucelRows(const double *tau,
                   int nt,
                   const double *rows,
                   int stride,
                   int n,
                   const int *mat,
                   const int *pdg,
                   const double *mass,
                   const int *rid,
                   const double *wb,
                   const double *frac,
                   int ne,
                   const double *wqRow,
                   const double *wbMid,
                   const RowConfig &cfg,
                   NucelMixtures &mix,
                   double *ang,
                   double *recRe,
                   double *recIm,
                   double *jntRe,
                   double *jntIm) {
    if (rows == nullptr || mat == nullptr || pdg == nullptr || n <= 0)
      return 0.;
    return nucelRowsImpl<double>(tau, nt, rows, stride, nullptr, n, mat, pdg, mass, rid, wb, frac, ne, wqRow,
                                 wbMid, cfg, mix, ang, recRe, recIm, jntRe, jntIm);
  }

  namespace {

    // The family of one functional from the MS rows and the block weights
    // the pooling pass formed: `wms` / `wio` are (global index, weight) of
    // the MS and (signed) ionisation blocks, ascending in index.  Each MS row
    // is one entry at its block's angular weight; its recoil weight is the
    // ionisation weight of the same leg at the same step (the radiative row
    // parallel to it: block index and cs), and the joint term takes both.
    void nucelFunctional(const TrackInput &sh,
                         const TrackInput &in,
                         const std::vector<std::pair<unsigned int, double>> &wms,
                         const std::vector<std::pair<unsigned int, double>> &wio,
                         const RowConfig &cfg,
                         TrackResult &out) {
      out.nucel = true;
      out.nuc.clear();
      out.nucGroups.clear();
      const StepRows &ms = sh.ms;
      if (ms.n <= 0)
        return;
      if (sh.msmat == nullptr || sh.mspdg == nullptr)
        throw cms::Exception("CvhCfExponents") << "nuclear-elastic family requested without msmat / mspdg";
      if (sh.rad.n != ms.n || sh.rad.v == nullptr || sh.rad.stride < 11)
        throw cms::Exception("CvhCfExponents")
            << "nuclear-elastic family: the radiative rows (" << sh.rad.n << ") are not parallel to the MS rows ("
            << ms.n << "); the recoil weight cannot be paired";
      auto look = [](const std::vector<std::pair<unsigned int, double>> &v, unsigned int g) {
        auto it = std::lower_bound(
            v.begin(), v.end(), g, [](const std::pair<unsigned int, double> &a, unsigned int b) { return a.first < b; });
        return (it != v.end() && it->first == g) ? it->second : 0.;
      };
      const int n = ms.n;
      std::vector<double> wA(n, 0.), wR(n, 0.), mass(n, 0.), fr(n, 1.);
      std::vector<int> rid(n);
      bool any = false;
      for (int i = 0; i < n; ++i) {
        rid[i] = i;
        if (sh.mspdg[i] == 0)
          continue;
        any = true;
        const float *r = ms.v + static_cast<std::size_t>(i) * ms.stride;
        const float *rr = sh.rad.v + static_cast<std::size_t>(i) * sh.rad.stride;
        // the pairing is exact by construction; check it on the momentum
        // both records carry (the same double, stored as float twice)
        if (rr[4] != r[3])
          throw cms::Exception("CvhCfExponents")
              << "nuclear-elastic family: radiative row " << i << " (p = " << rr[4]
              << ") is not the MS row's step (p = " << r[3] << ")";
        wA[i] = look(wms, ms.idx[i]);
        wR[i] = look(wio, sh.rad.idx[i]) * static_cast<double>(rr[10]) * 1e-3;
        mass[i] = nucelMassOf(sh.mspdg[i]);
      }
      if (!any)
        return;
      const double *tau = tauGrid();
      std::vector<double> wbE;
      auto run = [&](const int *sel, int nsel, NucelExponents &e) {
        // one entry per row, at the row's block weight
        wbE.resize(nsel);
        for (int q = 0; q < nsel; ++q)
          wbE[q] = wA[sel ? sel[q] : q];
        e.N += nucelRowsImpl<float>(tau,
                                    kNTau,
                                    ms.v,
                                    ms.stride,
                                    sel,
                                    nsel,
                                    sh.msmat,
                                    sh.mspdg,
                                    mass.data(),
                                    sel ? sel : rid.data(),
                                    wbE.data(),
                                    fr.data(),
                                    nsel,
                                    wR.data(),
                                    wA.data(),
                                    cfg,
                                    *in.nucel,
                                    e.ang.data(),
                                    e.recRe.data(),
                                    e.recIm.data(),
                                    e.jntRe.data(),
                                    e.jntIm.data());
      };
      const bool groups = in.wantGroups && ms.groupCol >= 0 && ms.groupCol < ms.stride;
      if (!groups) {
        run(nullptr, n, out.nuc);
        return;
      }
      // per material group, the flat family the sum of the groups
      std::map<int, std::vector<int>> byGroup;
      for (int i = 0; i < n; ++i)
        if (sh.mspdg[i] != 0)
          byGroup[static_cast<int>(ms.v[static_cast<std::size_t>(i) * ms.stride + ms.groupCol])].push_back(i);
      for (const auto &kv : byGroup) {
        NucelExponents ge;
        run(kv.second.data(), static_cast<int>(kv.second.size()), ge);
        if (ge.N == 0.)
          continue;
        out.nuc.N += ge.N;
        for (int j = 0; j < kNTau; ++j) {
          out.nuc.ang[j] += ge.ang[j];
          out.nuc.recRe[j] += ge.recRe[j];
          out.nuc.recIm[j] += ge.recIm[j];
          out.nuc.jntRe[j] += ge.jntRe[j];
          out.nuc.jntIm[j] += ge.jntIm[j];
        }
        out.nucGroups.emplace_back(kv.first, ge);
      }
    }

  }  // namespace

  //==========================================================================
  // 6. THE ROW FUNCTIONS  (cf_rows, cf_knockon, cf_brems_exact,
  //     cf_track_resolution.ioni_step_exponent)
  //
  // A port of the reference, statement by statement: the same node sets, the
  // same quadratures, and the reference's own association of every product
  // and sum, so the two agree to rounding and the closure that evaluates
  // through these functions (via the shim) tests the model the makers export.
  //==========================================================================
  const RowConfig &productionRowConfig() {
    static const RowConfig c;
    return c;
  }

  namespace {

    using cplx = std::complex<double>;

    //------------------------------------------------------------------------
    // NUMPY'S ARITHMETIC.  The knock-on map piece is the difference of two
    // Filon sums of the same spectrum -- in T_eff and in T -- which nearly
    // cancel, so the reference's last bits are visible at the scale the gate
    // measures.  The quadrature below therefore forms every complex product,
    // quotient and sum as numpy does: umath's complex division (Smith's
    // algorithm, with its reciprocal scale), umath's complex expm1, and the
    // pairwise summation of a contiguous reduction (np.sum along the last
    // axis, 8-way unrolled below 128 values).
    //------------------------------------------------------------------------
    double pwSum(const double *a, std::ptrdiff_t n) {
      if (n < 8) {
        double r = 0.;
        for (std::ptrdiff_t i = 0; i < n; ++i)
          r += a[i];
        return r;
      }
      if (n <= 128) {
        double r[8];
        for (int k = 0; k < 8; ++k)
          r[k] = a[k];
        std::ptrdiff_t i = 8;
        for (; i < n - (n % 8); i += 8)
          for (int k = 0; k < 8; ++k)
            r[k] += a[i + k];
        double res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        for (; i < n; ++i)
          res += a[i];
        return res;
      }
      std::ptrdiff_t n2 = n / 2;
      n2 -= n2 % 8;
      return pwSum(a, n2) + pwSum(a + n2, n - n2);
    }

    // The complex reduction counts DOUBLES, so its lanes are four complex
    // values wide and its split point is in doubles.
    cplx pwSumC(const cplx *a, std::ptrdiff_t m) {
      const std::ptrdiff_t n = 2 * m;
      if (n < 8) {
        double rr = 0., ri = 0.;
        for (std::ptrdiff_t i = 0; i < m; ++i) {
          rr += a[i].real();
          ri += a[i].imag();
        }
        return cplx(rr, ri);
      }
      if (n <= 128) {
        double r[8];
        for (int k = 0; k < 4; ++k) {
          r[2 * k] = a[k].real();
          r[2 * k + 1] = a[k].imag();
        }
        std::ptrdiff_t i = 8;
        for (; i < n - (n % 8); i += 8)
          for (int k = 0; k < 4; ++k) {
            r[2 * k] += a[i / 2 + k].real();
            r[2 * k + 1] += a[i / 2 + k].imag();
          }
        double rr = (r[0] + r[2]) + (r[4] + r[6]);
        double ri = (r[1] + r[3]) + (r[5] + r[7]);
        for (; i < n; i += 2) {
          rr += a[i / 2].real();
          ri += a[i / 2].imag();
        }
        return cplx(rr, ri);
      }
      std::ptrdiff_t n2 = n / 2;
      n2 -= n2 % 8;
      const cplx s1 = pwSumC(a, n2 / 2), s2 = pwSumC(a + n2 / 2, (n - n2) / 2);
      return cplx(s1.real() + s2.real(), s1.imag() + s2.imag());
    }

    // np.interp (numpy's compiled `interp`: the bracketing node found by
    // binary search, slope * (x - xp_j) + fp_j, the end values beyond).
    double npInterp(double x, const double *xp, const double *fp, int n) {
      if (x > xp[n - 1])
        return fp[n - 1];
      if (x < xp[0])
        return fp[0];
      int lo = 0, hi = n;
      while (hi - lo > 1) {
        const int mid = (lo + hi) >> 1;
        if (x >= xp[mid])
          lo = mid;
        else
          hi = mid;
      }
      const int j = lo;
      if (j == n - 1)
        return fp[j];
      if (xp[j] == x)
        return fp[j];
      const double slope = (fp[j + 1] - fp[j]) / (xp[j + 1] - xp[j]);
      return slope * (x - xp[j]) + fp[j];
    }

  }  // namespace

  //------------------------------------------------------------------------
  // SCIPY'S BESSEL FUNCTIONS.  J0 (the knock-on joint law) and K1 (the
  // muon's photon-angle law) are Cephes' rational/Chebyshev approximations,
  // the ones scipy.special evaluates, rather than the toolchain's, so the
  // model does not depend on which library the build links.
  //------------------------------------------------------------------------
  namespace {
    double polevl(double x, const double *c, int N) {
      double ans = c[0];
      for (int i = 1; i <= N; ++i)
        ans = ans * x + c[i];
      return ans;
    }
    double p1evl(double x, const double *c, int N) {
      double ans = x + c[0];
      for (int i = 1; i < N; ++i)
        ans = ans * x + c[i];
      return ans;
    }
    double chbevl(double x, const double *a, int n) {
      double b0 = a[0], b1 = 0., b2 = 0.;
      for (int i = 1; i < n; ++i) {
        b2 = b1;
        b1 = b0;
        b0 = x * b1 - b2 + a[i];
      }
      return 0.5 * (b0 - b2);
    }

    constexpr double kJ0PP[7] = {7.96936729297347051624E-4, 8.28352392107440799803E-2, 1.23953371646414299388E0,
                                 5.44725003058768775090E0,  8.74716500199817011941E0,  5.30324038235394892183E0,
                                 9.99999999999999997821E-1};
    constexpr double kJ0PQ[7] = {9.24408810558863637013E-4, 8.56288474354474431428E-2, 1.25352743901058953537E0,
                                 5.47097740330417105182E0,  8.76190883237069594232E0,  5.30605288235394617618E0,
                                 1.00000000000000000218E0};
    constexpr double kJ0QP[8] = {-1.13663838898469149931E-2, -1.28252718670509318512E0, -1.95539544257735972385E1,
                                 -9.32060152123768231369E1,  -1.77681167980488050595E2, -1.47077505154951170175E2,
                                 -5.14105326766599330220E1,  -6.05014350600728481186E0};
    constexpr double kJ0QQ[7] = {6.43178256118178023184E1, 8.56430025976980587198E2, 3.88240183605401609683E3,
                                 7.24046774195652478189E3, 5.93072701187316984827E3, 2.06209331660327847417E3,
                                 2.42005740240291393179E2};
    constexpr double kJ0DR1 = 5.78318596294678452118E0;
    constexpr double kJ0DR2 = 3.04712623436620863991E1;
    constexpr double kJ0RP[4] = {-4.79443220978201773821E9, 1.95617491946556577543E12, -2.49248344360967716204E14,
                                 9.70862251047306323952E15};
    constexpr double kJ0RQ[8] = {4.99563147152651017219E2,  1.73785401676374683123E5,  4.84409658339962045305E7,
                                 1.11855537045356834862E10, 2.11277520115489217587E12, 3.10518229857422583814E14,
                                 3.18121955943204943306E16, 1.71086294081043136091E18};
    constexpr double kSqrt2OPi = 7.9788456080286535587989E-1;
    constexpr double kPiO4 = 0.78539816339744830962;

    constexpr double kI1A[29] = {
        2.77791411276104639959E-18, -2.11142121435816608115E-17, 1.55363195773620046921E-16,
        -1.10559694773538630805E-15, 7.60068429473540693410E-15,  -5.04218550472791168711E-14,
        3.22379336594557470981E-13,  -1.98397439776494371520E-12, 1.17361862988909016308E-11,
        -6.66348972350202774223E-11, 3.62559028155211703701E-10,  -1.88724975172282928790E-9,
        9.38153738649577178388E-9,   -4.44505912879632808065E-8,  2.00329475355213526229E-7,
        -8.56872026469545474066E-7,  3.47025130813767847674E-6,   -1.32731636560394358279E-5,
        4.78156510755005422638E-5,   -1.61760815825896745588E-4,  5.12285956168575772895E-4,
        -1.51357245063125314899E-3,  4.15642294431288815669E-3,   -1.05640848946261981558E-2,
        2.47264490306265168283E-2,   -5.29459812080949914269E-2,  1.02643658689847095384E-1,
        -1.76416518357834055153E-1,  2.52587186443633654823E-1};
    constexpr double kI1B[25] = {
        7.51729631084210481353E-18,  4.41434832307170791151E-18,  -4.65030536848935832153E-17,
        -3.20952592199342395980E-17, 2.96262899764595013876E-16,  3.30820231092092828324E-16,
        -1.88035477551078244854E-15, -3.81440307243700780478E-15, 1.04202769841288027642E-14,
        4.27244001671195135429E-14,  -2.10154184277266431302E-14, -4.08355111109219731823E-13,
        -7.19855177624590851209E-13, 2.03562854414708950722E-12,  1.41258074366137813316E-11,
        3.25260358301548823856E-11,  -1.89749581235054123450E-11, -5.58974346219658380687E-10,
        -3.83538038596423702205E-9,  -2.63146884688951950684E-8,  -2.51223623787020892529E-7,
        -3.88256480887769039346E-6,  -1.10588938762623716291E-4,  -9.76109749136146840777E-3,
        7.78576235018280120474E-1};
    constexpr double kK1A[11] = {-7.02386347938628759343E-18, -2.42744985051936593393E-15, -6.66690169419932900609E-13,
                                 -1.41148839263352776110E-10, -2.21338763073472585583E-8,  -2.43340614156596823496E-6,
                                 -1.73028895751305206302E-4,  -6.97572385963986435018E-3,  -1.22611180822657148235E-1,
                                 -3.53155960776544875667E-1,  1.52530022733894777053E0};
    constexpr double kK1B[25] = {
        -5.75674448366501715755E-18, 1.79405087314755922667E-17,  -5.68946255844285935196E-17,
        1.83809354436663880070E-16,  -6.05704724837331885336E-16, 2.03870316562433424052E-15,
        -7.01983709041831346144E-15, 2.47715442448130437068E-14,  -8.97670518232499435011E-14,
        3.34841966607842919884E-13,  -1.28917396095102890680E-12, 5.13963967348173025100E-12,
        -2.12996783842756842877E-11, 9.21831518760500529508E-11,  -4.19035475934189648750E-10,
        2.01504975519703286596E-9,   -1.03457624656780970260E-8,  5.74108412545004946722E-8,
        -3.50196060308781257119E-7,  2.40648494783721712015E-6,   -1.93619797416608296024E-5,
        1.95215518471351631108E-4,   -2.85781685962277938680E-3,  1.03923736576817238437E-1,
        2.72062619048444266945E0};

    double besselI1(double x) {
      double z = std::fabs(x);
      if (z <= 8.0) {
        const double y = (z / 2.0) - 2.0;
        z = chbevl(y, kI1A, 29) * z * std::exp(z);
      } else {
        z = std::exp(z) * chbevl(32.0 / z - 2.0, kI1B, 25) / std::sqrt(z);
      }
      return (x < 0.0) ? -z : z;
    }

  }  // namespace

  double besselJ0(double x) {
    if (x < 0)
      x = -x;
    if (x <= 5.0) {
      const double z = x * x;
      if (x < 1.0e-5)
        return 1.0 - z / 4.0;
      double p = (z - kJ0DR1) * (z - kJ0DR2);
      p = p * polevl(z, kJ0RP, 3) / p1evl(z, kJ0RQ, 8);
      return p;
    }
    const double w = 5.0 / x;
    double q = 25.0 / (x * x);
    double p = polevl(q, kJ0PP, 6) / polevl(q, kJ0PQ, 6);
    q = polevl(q, kJ0QP, 7) / p1evl(q, kJ0QQ, 7);
    const double xn = x - kPiO4;
    p = p * std::cos(xn) - w * q * std::sin(xn);
    return p * kSqrt2OPi / std::sqrt(x);
  }

  double besselK1(double x) {
    if (x == 0.0)
      return std::numeric_limits<double>::infinity();
    if (x < 0.0)
      return std::numeric_limits<double>::quiet_NaN();
    const double z = 0.5 * x;
    if (x <= 2.0) {
      const double y = x * x - 2.0;
      return std::log(z) * besselI1(x) + chbevl(y, kK1A, 11) / x;
    }
    return std::exp(-x) * chbevl(8.0 / x - 2.0, kK1B, 25) / std::sqrt(x);
  }

  namespace {

    //------------------------------------------------------------------------
    // THE KNOCK-ON QUADRATURE (cf_knockon.nodes, _fmom, filon, simpson)
    //------------------------------------------------------------------------
    constexpr double kKoMeMeV = 0.51099895;
    constexpr double kKoMuMeV = 105.6583745;
    constexpr double kKoKokMuMin = 1000.0;  // G4MuBetheBlochModel::lowestKinEnergy [MeV]
    // the charged species a track can be: e, mu, pi, K, p [GeV]
    // (cf_brems_exact._SPECIES_MASSES)
    constexpr double kSpeciesMassesGeV[5] = {0.51099895e-3, 0.1056583745, 0.13957039, 0.493677, 0.93827208816};
    constexpr int kSpeciesMu = 1;

    // The species of an ionisation record (cf_track_resolution.
    // _record_species_mass), as an index into kSpeciesMassesGeV: the mass
    // E sqrt(1 - beta^2) [MeV] its float32 columns give -- to a few MeV at high
    // momentum -- snapped to the nearest species the record's regime admits:
    // 2 the spin-1/2 law (mu, p), 3 the spin-0 law (pi, K), 4/5 e+-, any other
    // regime every species.
    inline int recordSpecies(double b2, double E, int reg) {
      static constexpr int c2[2] = {1, 4}, c3[2] = {2, 3}, c45[1] = {0}, call[5] = {0, 1, 2, 3, 4};
      const int *c = call;
      int nc = 5;
      if (reg == 2) {
        c = c2;
        nc = 2;
      } else if (reg == 3) {
        c = c3;
        nc = 2;
      } else if (reg == 4 || reg == 5) {
        c = c45;
        nc = 1;
      }
      const double m = E * std::sqrt(std::max(1.0 - b2, 0.0));
      int best = c[0];
      double bd = std::fabs(1e3 * kSpeciesMassesGeV[c[0]] - m);
      for (int k = 1; k < nc; ++k) {
        const double d = std::fabs(1e3 * kSpeciesMassesGeV[c[k]] - m);
        if (d < bd) {
          bd = d;
          best = c[k];
        }
      }
      return best;
    }
    constexpr int kFmomNser = 14;           // cf_knockon._NSER
    // 1/m for the series' divisors m = 1 .. kFmomNser + 2, the same doubles
    // 1.0 / m gives, looked up rather than divided per term
    struct Recips {
      double v[kFmomNser + 3];
      constexpr Recips() : v() {
        for (int m = 1; m < kFmomNser + 3; ++m)
          v[m] = 1.0 / static_cast<double>(m);
      }
    };
    constexpr Recips kRecip;
    constexpr double kFmomZser = 0.2;       // cf_knockon._ZSER

    // np.geomspace(lo, hi, n), appended: 10**linspace(log10 lo, log10 hi)
    // with both ends forced.
    void geomspaceApp(double lo, double hi, int n, std::vector<double> &out) {
      const double a = std::log10(lo), b = std::log10(hi);
      const double step = (b - a) / static_cast<double>(n - 1);
      const std::size_t first = out.size();
      for (int i = 0; i < n; ++i)
        out.push_back(std::pow(10.0, static_cast<double>(i) * step + a));
      out[first] = lo;
      out.back() = hi;
    }

    // cf_knockon.nodes: nper per decade on [max(e0, tlo), 0.9 tmax], 30
    // log-spaced in (tmax - T) down to 1e-7 tmax, and tmax; odd count.
    void koNodes(double e0, double tmax, int nper, double tlo, std::vector<double> &x) {
      x.clear();
      const double lo = std::max(e0, tlo);
      const double top = 0.9 * tmax;
      const int n1 = std::max(static_cast<int>(std::ceil(std::log10(top / lo) * nper)), 2);
      geomspaceApp(lo, top, n1 + 1, x);
      std::vector<double> d;
      geomspaceApp(0.1, 1e-7, 31, d);
      for (int i = 1; i < 31; ++i)
        x.push_back(tmax - tmax * d[i]);
      x.push_back(tmax);
      if (x.size() % 2 == 0)
        x.insert(x.begin() + 1, std::sqrt(x[0] * x[1]));
    }

    // f_n(z) = INT_0^1 s^n e^{z s} ds, n = 0, 1, 2 (cf_knockon._fmom), at the
    // purely imaginary z = i y the Filon rule meets (z = 1j a D), in numpy's
    // complex arithmetic: the series below |z| = 0.2, the closed forms beyond.
    void fmomI(double y, cplx f[3]) {
      if (std::fabs(y) < kFmomZser) {
        double tr = 1.0, ti = 0.0;
        double ar[3] = {0., 0., 0.}, ai[3] = {0., 0., 0.};
        for (int k = 0; k < kFmomNser; ++k) {
          if (k) {
            // term * w / k: the product with w = i y, then umath's division
            // by the real k (a multiplication by its reciprocal)
            const double mr = -ti * y, mi = tr * y;
            const double scl = kRecip.v[k];
            tr = mr * scl;
            ti = mi * scl;
          }
          for (int n = 0; n < 3; ++n) {
            const double scl = kRecip.v[n + k + 1];
            ar[n] += tr * scl;
            ai[n] += ti * scl;
          }
        }
        for (int n = 0; n < 3; ++n)
          f[n] = cplx(ar[n], ai[n]);
        return;
      }
      const double c = std::cos(y), s = std::sin(y);  // np.exp(w)
      // f0 = expm1(w) / w
      const double sh = std::sin(y / 2);
      const double emr = -(2 * sh * sh), emi = s;
      {
        const double scl = 1.0 / y;
        f[0] = cplx(emi * scl, -emr * scl);
      }
      // f1 = (e^w (w - 1) + 1) / w^2
      {
        const double nr = (c * -1.0 - s * y) + 1.0;
        const double ni = c * y + s * -1.0;
        const double scl = 1.0 / (-(y * y));
        f[1] = cplx(nr * scl, ni * scl);
      }
      // f2 = (e^w (w^2 - 2 w + 2) - 2) / w^3
      {
        const double br = -(y * y) + 2.0, bi = -(2.0 * y);
        const double nr = (c * br - s * bi) - 2.0;
        const double ni = c * bi + s * br;
        const double scl = 1.0 / ((-(y * y)) * y);
        f[2] = cplx(ni * scl, -nr * scl);
      }
    }

    // The panels of one amplitude h on the odd node set x (Newton form of the
    // quadratic through each panel's three nodes).
    struct Panels {
      std::vector<double> x0, D, h0, c1, c2;
      int n() const { return static_cast<int>(D.size()); }
    };
    void makePanels(const double *x, const double *h, int nx, Panels &P) {
      const int np = (nx - 1) / 2;
      P.x0.resize(np);
      P.D.resize(np);
      P.h0.resize(np);
      P.c1.resize(np);
      P.c2.resize(np);
      for (int i = 0; i < np; ++i) {
        const double x0 = x[2 * i], x1 = x[2 * i + 1], x2 = x[2 * i + 2];
        const double h0 = h[2 * i], h1 = h[2 * i + 1], h2 = h[2 * i + 2];
        const double d01 = (h1 - h0) / (x1 - x0);
        const double d12 = (h2 - h1) / (x2 - x1);
        const double c2 = (d12 - d01) / (x2 - x0);
        P.c1[i] = d01 - c2 * (x1 - x0);
        P.c2[i] = c2;
        P.x0[i] = x0;
        P.D[i] = x2 - x0;
        P.h0[i] = h0;
      }
    }

    // cf_knockon.filon at one a: INT e^{i a x} h dx, phase exact per panel.
    cplx filonAt(const Panels &P, double a, std::vector<cplx> &buf) {
      const int np = P.n();
      buf.resize(np);
      cplx f[3];
      for (int i = 0; i < np; ++i) {
        const double D = P.D[i];
        fmomI(a * D, f);
        const double ph = a * P.x0[i];
        const double t1r = D * std::cos(ph), t1i = D * std::sin(ph);
        const double k1 = P.c1[i] * D, k2 = P.c2[i] * D * D;
        const double ur = (P.h0[i] * f[0].real() + k1 * f[1].real()) + k2 * f[2].real();
        const double ui = (P.h0[i] * f[0].imag() + k1 * f[1].imag()) + k2 * f[2].imag();
        buf[i] = cplx(t1r * ur - t1i * ui, t1r * ui + t1i * ur);
      }
      return pwSumC(buf.data(), np);
    }

    // cf_knockon._fmom_centred at z = i y: f_n (`fmomI`), g_n = INT s^n
    // (e^{zs} - 1) ds and q_n = INT s^n (e^{zs} - 1 - zs) ds, by the series
    // below |z| = _ZSER (no cancellation) and from the closed-form f_n above.
    void fmomCentredI(double y, cplx f[3], cplx g[3], cplx q[3]) {
      if (std::fabs(y) < kFmomZser) {
        double tr = 1.0, ti = 0.0;
        double fr[3] = {0., 0., 0.}, fi[3] = {0., 0., 0.};
        double gr[3] = {0., 0., 0.}, gi[3] = {0., 0., 0.};
        double qr[3] = {0., 0., 0.}, qi[3] = {0., 0., 0.};
        for (int k = 0; k < kFmomNser; ++k) {
          if (k) {
            const double mr = -ti * y, mi = tr * y;
            const double scl = kRecip.v[k];
            tr = mr * scl;
            ti = mi * scl;
          }
          for (int n = 0; n < 3; ++n) {
            const double scl = kRecip.v[n + k + 1];
            const double t_r = tr * scl, t_i = ti * scl;
            fr[n] += t_r;
            fi[n] += t_i;
            if (k >= 1) {
              gr[n] += t_r;
              gi[n] += t_i;
            }
            if (k >= 2) {
              qr[n] += t_r;
              qi[n] += t_i;
            }
          }
        }
        for (int n = 0; n < 3; ++n) {
          f[n] = cplx(fr[n], fi[n]);
          g[n] = cplx(gr[n], gi[n]);
          q[n] = cplx(qr[n], qi[n]);
        }
        return;
      }
      fmomI(y, f);
      for (int n = 0; n < 3; ++n) {
        g[n] = cplx(f[n].real() - 1.0 / (n + 1), f[n].imag());
        q[n] = cplx(g[n].real(), g[n].imag() - y * (1.0 / static_cast<double>(n + 2)));
      }
    }

    // sin(y) - y (cf_knockon._sinm): the Taylor series below |y| = 0.5
    inline double sinm(double y) {
      static constexpr double c[8] = {-1.0 / 6.0,
                                      1.0 / 120.0,
                                      -1.0 / 5040.0,
                                      1.0 / 362880.0,
                                      -1.0 / 39916800.0,
                                      1.0 / 6227020800.0,
                                      -1.0 / 1307674368000.0,
                                      1.0 / 355687428096000.0};
      if (!(std::fabs(y) < 0.5))
        return std::sin(y) - y;
      const double y2 = y * y;
      double p = c[7];
      for (int k = 6; k >= 0; --k)
        p = p * y2 + c[k];
      return (p * y2) * y;
    }

    // cf_knockon.filon_centred at one a: INT h (e^{i a x} - 1 - i a x) dx on
    // the panels of `filonAt`, the centring subtracted inside each panel:
    // E0 e^{zs} + (e^{zs} - 1 - zs) + i y (e^{zs} - 1), E0 = e^{iy} - 1 - iy.
    cplx filonCentredAt(const Panels &P, double a, std::vector<cplx> &buf) {
      const int np = P.n();
      buf.resize(np);
      cplx f[3], g[3], q[3];
      for (int i = 0; i < np; ++i) {
        const double D = P.D[i];
        const double y = a * P.x0[i];
        fmomCentredI(a * D, f, g, q);
        const double hy = std::sin(0.5 * y);
        const double er = -2.0 * hy * hy, ei = sinm(y);
        const double k1 = P.c1[i] * D, k2 = P.c2[i] * D * D, h0 = P.h0[i];
        const double fr = (h0 * f[0].real() + k1 * f[1].real()) + k2 * f[2].real();
        const double fi = (h0 * f[0].imag() + k1 * f[1].imag()) + k2 * f[2].imag();
        const double qr = (h0 * q[0].real() + k1 * q[1].real()) + k2 * q[2].real();
        const double qi = (h0 * q[0].imag() + k1 * q[1].imag()) + k2 * q[2].imag();
        const double gr = (h0 * g[0].real() + k1 * g[1].real()) + k2 * g[2].real();
        const double gi = (h0 * g[0].imag() + k1 * g[1].imag()) + k2 * g[2].imag();
        // E0 * F + Q + (i y) * G
        const double tr = ((er * fr - ei * fi) + qr) + (-y * gi);
        const double ti = ((er * fi + ei * fr) + qi) + y * gr;
        buf[i] = cplx(D * tr, D * ti);
      }
      return pwSumC(buf.data(), np);
    }

    // cf_knockon.simpson: the same panels at a = 0.
    double simpsonOf(const Panels &P, std::vector<double> &buf) {
      const int np = P.n();
      buf.resize(np);
      for (int i = 0; i < np; ++i) {
        const double D = P.D[i];
        buf[i] = D * ((P.h0[i] + P.c1[i] * D / 2.0) + P.c2[i] * D * D / 3.0);
      }
      return pwSum(buf.data(), np);
    }

    //------------------------------------------------------------------------
    // THE KNOCK-ON LAW (cf_knockon.rate, law_rate, law_top, t_eff, dteff_dt,
    // theta_kick) on a node set.
    //------------------------------------------------------------------------
    // Kokoulin's factor f_K(T) of G4MuBetheBlochModel, as the reference forms
    // it (cf_track_resolution.kokoulin_factor): its electron mass is the
    // reference's, 0.51099895 MeV, where the in-fit CGF's `cvhcgf::
    // kokoulinFactor` takes CLHEP's 0.510998910 (Geant4's own) -- a 8e-8
    // relative difference in ln(1 + 2T/m_e) that no physics sees and a
    // bit-level comparison does.
    constexpr double kKokAlphaPrime = 1.0 / (2.0 * 3.14159265358979323846 * 137.035999084);
    constexpr double kKokTMin = 0.1;  // G4MuBetheBlochModel::limitKinEnergy [MeV]
    // f_K - 1 on the correction's range (100 keV, E - m_mu), where the caller
    // keeps T (cf_track_resolution.kokoulin_excess)
    inline double kokExcess(double T, double E) {
      const double a1 = std::log(1.0 + 2.0 * T / kKoMeMeV);
      const double a3 = std::log(4.0 * E * (E - T) / (kKoMuMeV * kKoMuMeV));
      return kKokAlphaPrime * a1 * (a3 - a1);
    }
    inline double kokFactor(double T, double E) {
      if (!(T > kKokTMin) || !(T < E - kKoMuMeV))
        return 1.0;
      return 1.0 + kokExcess(T, E);
    }

    //------------------------------------------------------------------------
    // THE EXACT KNOCK-ON SPECTRUM IN CLOSED FORM (cf_track_resolution.
    // exact_delta_exponent, _delta_terms_exact, _ein_neg) in numpy's
    // arithmetic -- e^{i alpha} - 1 as (-2 sin^2(alpha/2), sin alpha), not
    // cos alpha - 1, which loses the real part below alpha ~ 1e-8 where an
    // e+- record's xi/e0 ~ 1e3 multiplies it -- on a BATCH of (step, tau)
    // pairs: the reference evaluates one call's pairs together and stops its
    // series when the whole call has converged, so a batch is one reference
    // call.  Ein(-i u) beyond |u| = 2 is Cin(|u|) - i sgn(u) Si(|u|)
    // (`cvhcgf::siCin`) where the reference takes gamma + log(-s) + E1(-s).
    //------------------------------------------------------------------------
    constexpr double kDtSer = 2.0;  // cf_track_resolution._DT_SER
    constexpr int kDtNser = 80;     // cf_track_resolution._DT_NSER

    // Ein(-s) at s = i u for a batch of u
    void einNegBatch(const std::vector<double> &u, std::vector<double> &er, std::vector<double> &ei) {
      const std::size_t n = u.size();
      er.assign(n, 0.);
      ei.assign(n, 0.);
      std::vector<std::size_t> sm;
      for (std::size_t i = 0; i < n; ++i) {
        if (std::fabs(u[i]) <= kDtSer) {
          sm.push_back(i);
        } else {
          double si, cin;
          cvhcgf::siCin(std::fabs(u[i]), si, cin);
          er[i] = cin;
          ei[i] = (u[i] >= 0.) ? -si : si;
        }
      }
      if (sm.empty())
        return;
      const std::size_t m = sm.size();
      std::vector<double> tr(m, 1.), ti(m, 0.), ar(m, 0.), ai(m, 0.);
      for (int k = 1; k < kDtNser; ++k) {
        const double sk = 1.0 / static_cast<double>(k);
        double mt = 0., ma = 0.;
        for (std::size_t q = 0; q < m; ++q) {
          const double uq = u[sm[q]];
          // term = term * s / k, s = (+-0, u)
          const double mr = -ti[q] * uq, mi = tr[q] * uq;
          tr[q] = mr * sk;
          ti[q] = mi * sk;
          // acc = acc - term / k
          ar[q] = ar[q] - tr[q] * sk;
          ai[q] = ai[q] - ti[q] * sk;
          mt = std::max(mt, std::hypot(tr[q], ti[q]));
          ma = std::max(ma, std::hypot(ar[q], ai[q]));
        }
        if (k > 3 && mt < 1e-19 * std::max(ma, 1e-300))
          break;
      }
      for (std::size_t q = 0; q < m; ++q) {
        er[sm[q]] = ar[q];
        ei[sm[q]] = ai[q];
      }
    }

    // (J0, J1, J2) at alpha = a e0, w = tmax/e0, for a batch
    void deltaTermsBatch(const std::vector<double> &al,
                         const std::vector<double> &w,
                         std::vector<cplx> &J0,
                         std::vector<cplx> &J1,
                         std::vector<cplx> &J2) {
      const std::size_t n = al.size();
      J0.assign(n, cplx(0., 0.));
      J1.assign(n, cplx(0., 0.));
      J2.assign(n, cplx(0., 0.));
      std::vector<std::size_t> ser, bg;
      for (std::size_t i = 0; i < n; ++i)
        (std::fabs(al[i] * w[i]) <= kDtSer ? ser : bg).push_back(i);
      if (!ser.empty()) {
        const std::size_t m = ser.size();
        std::vector<double> txr(m, 1.), txi(m, 0.), tyr(m, 1.), tyi(m, 0.);
        std::vector<double> a0r(m, 0.), a0i(m, 0.), a1r(m, 0.), a1i(m, 0.), a2r(m, 0.), a2i(m, 0.);
        double mw = 0.;
        for (std::size_t q = 0; q < m; ++q)
          mw = (q == 0) ? w[ser[q]] : std::max(mw, w[ser[q]]);
        for (int k = 1; k < kDtNser; ++k) {
          const double sk = 1.0 / static_cast<double>(k);
          double mt = 0., m2 = 0.;
          for (std::size_t q = 0; q < m; ++q) {
            const std::size_t i = ser[q];
            const double x = al[i] * w[i], y = al[i];
            // tx = tx * x / k, ty = ty * y / k (x, y purely imaginary)
            double r = -txi[q] * x, im = txr[q] * x;
            txr[q] = r * sk;
            txi[q] = im * sk;
            r = -tyi[q] * y;
            im = tyr[q] * y;
            tyr[q] = r * sk;
            tyi[q] = im * sk;
            if (k >= 2) {
              const double iw = 1.0 / w[i];
              const double s0 = 1.0 / (k - 1.0), s2 = 1.0 / (k + 1.0);
              a0r[q] = a0r[q] + (txr[q] * iw - tyr[q]) * s0;
              a0i[q] = a0i[q] + (txi[q] * iw - tyi[q]) * s0;
              a1r[q] = a1r[q] + (txr[q] - tyr[q]) * sk;
              a1i[q] = a1i[q] + (txi[q] - tyi[q]) * sk;
              a2r[q] = a2r[q] + (txr[q] * w[i] - tyr[q]) * s2;
              a2i[q] = a2i[q] + (txi[q] * w[i] - tyi[q]) * s2;
              mt = std::max(mt, std::hypot(txr[q], txi[q]));
              m2 = std::max(m2, std::hypot(a2r[q], a2i[q]));
            }
          }
          if (k >= 2 && mt * mw < 1e-19 * std::max(m2, 1e-300))
            break;
        }
        for (std::size_t q = 0; q < m; ++q) {
          J0[ser[q]] = cplx(a0r[q], a0i[q]);
          J1[ser[q]] = cplx(a1r[q], a1i[q]);
          J2[ser[q]] = cplx(a2r[q], a2i[q]);
        }
      }
      if (!bg.empty()) {
        const std::size_t m = bg.size();
        std::vector<double> uy(m), ux(m), eyr, eyi, exr, exi;
        for (std::size_t q = 0; q < m; ++q) {
          uy[q] = al[bg[q]];
          ux[q] = al[bg[q]] * w[bg[q]];
        }
        einNegBatch(uy, eyr, eyi);
        einNegBatch(ux, exr, exi);
        for (std::size_t q = 0; q < m; ++q) {
          const std::size_t i = bg[q];
          const double a = al[i], ww = w[i], aw = ux[q];
          const double er = eyr[q] - exr[q], ei = eyi[q] - exi[q];
          const double hy = std::sin(a / 2), hx = std::sin(aw / 2);
          const double my_r = -(2 * hy * hy), my_i = std::sin(a);
          const double mx_r = -(2 * hx * hx), mx_i = std::sin(aw);
          const double iw = 1.0 / ww;
          J0[i] = cplx((my_r - mx_r * iw) + (-(a * ei)), (my_i - mx_i * iw) + a * er);
          J1[i] = cplx(er, ei - a * (ww - 1.0));
          const double nr = std::cos(aw) - std::cos(a), ni = std::sin(aw) - std::sin(a);
          const double ia = 1.0 / a;
          J2[i] = cplx(ni * ia - (ww - 1.0), -nr * ia - (a * (ww * ww - 1.0)) * 0.5);
        }
      }
    }

    // One reference call of exact_delta_exponent: sum over the steps of
    // xi/e0 [J0 - (beta^2 e0/tmax) J1 (+ spin 1/2: e0^2/(2E^2) J2)] at
    // a = gs tau, ACCUMULATED into S.
    struct ExStep {
      double xi, e0, tmax, beta2, etot, gs;
    };
    void exactDeltaBatch(const std::vector<ExStep> &st, const double *tau, int nt, bool half, double *Sre, double *Sim) {
      if (st.empty())
        return;
      std::vector<double> al(st.size() * nt), w(st.size() * nt);
      for (std::size_t s = 0; s < st.size(); ++s)
        for (int j = 0; j < nt; ++j) {
          al[s * nt + j] = (st[s].gs * tau[j]) * st[s].e0;
          w[s * nt + j] = st[s].tmax / st[s].e0;
        }
      std::vector<cplx> J0, J1, J2;
      deltaTermsBatch(al, w, J0, J1, J2);
      std::vector<double> accr(nt, 0.), acci(nt, 0.);
      for (std::size_t s = 0; s < st.size(); ++s) {
        const ExStep &e = st[s];
        const double pref = e.xi / e.e0, c1 = e.beta2 * e.e0 / e.tmax;
        const double coef = pref * ((e.e0 * e.e0) / (2.0 * (e.etot * e.etot)));
        for (int j = 0; j < nt; ++j) {
          const std::size_t i = s * nt + j;
          double vr = pref * (J0[i].real() - c1 * J1[i].real());
          double vi = pref * (J0[i].imag() - c1 * J1[i].imag());
          if (half) {
            vr = vr + coef * J2[i].real();
            vi = vi + coef * J2[i].imag();
          }
          accr[j] += vr;
          acci[j] += vi;
        }
      }
      for (int j = 0; j < nt; ++j) {
        Sre[j] += accr[j];
        Sim[j] += acci[j];
      }
    }

    void lawRate(const std::vector<double> &T,
                 double xi,
                 double tmax,
                 double b2,
                 double E,
                 int reg,
                 bool kok,
                 const RowConfig &cfg,
                 std::vector<double> &r) {
      const std::size_t n = T.size();
      r.resize(n);
      if (reg == 2 || reg == 3) {
        for (std::size_t i = 0; i < n; ++i) {
          double v = (xi / (T[i] * T[i])) * (1.0 - b2 * T[i] / tmax);
          if (reg == 2)
            v = v + xi / (2.0 * E * E);
          if (kok) {
            double f = kokFactor(T[i], E);
            if (!(T[i] > cfg.ioniKokoulinTcut))
              f = 1.0;
            v = v * (1.0 + cfg.ioniKokoulin * (f - 1.0));
          }
          r[i] = v;
        }
        return;
      }
      const double gam = E / kKoMeMeV;
      if (reg == 4) {
        const double gg = (2.0 * gam - 1.0) / (gam * gam);
        for (std::size_t i = 0; i < n; ++i) {
          const double x = T[i] / tmax;
          const double y = x / (1.0 - x);
          const double v = ((1.0 + x * x * (1.0 - gg)) + y * y) - gg * y;
          r[i] = (x <= 0.5) ? (xi / (T[i] * T[i])) * v : 0.0;
        }
        return;
      }
      if (reg == 5) {
        const double yy = 1.0 / (1.0 + gam);
        const double y12 = 1.0 - 2.0 * yy;
        const double c1 = 2.0 - yy * yy;
        const double c2 = y12 * (3.0 + yy * yy);
        const double c4 = std::pow(y12, 3.0);
        const double c3 = c4 + y12 * y12;
        for (std::size_t i = 0; i < n; ++i) {
          const double x = T[i] / tmax;
          const double v =
              1.0 + b2 * (((-c1 * x + c2 * x * x) - c3 * std::pow(x, 3.0)) + c4 * std::pow(x, 4.0));
          r[i] = (x <= 1.0) ? (xi / (T[i] * T[i])) * v : 0.0;
        }
        return;
      }
      throw cms::Exception("CvhCfExponents") << "knock-on law: no regime " << reg;
    }

    // where the law's channel ends: the record's tmax, T0/2 for Moller, and
    // never past the PMIN_FRAC floor on the primary's outgoing momentum
    inline double lawTop(int reg, double tmax, double E, double p, const RowConfig &cfg) {
      const double top = (reg == 4) ? 0.5 * tmax : tmax;
      const double q = cfg.pminFrac * p;
      const double tcap = E - std::sqrt(E * E - p * p + q * q);
      return (tcap < top) ? tcap : top;
    }

    inline double ppOfT(double T, double E, double p) {
      const double M2 = E * E - p * p;
      const double d = E - T;
      return std::sqrt(std::max(d * d - M2, 1e-300));
    }

    //------------------------------------------------------------------------
    // THE IONISATION CHANNEL'S REGIME-4/5 LAW (cf_knockon.law_base_exponent):
    // INT_{e0}^{top} dN (e^{i a T} - 1 - i a T) dT, the 1/T^2 part in closed
    // form and the law's O(T/T0) remainder by Filon-Simpson.
    //------------------------------------------------------------------------
    void lawBaseExponent(double xi,
                         double e0,
                         double tmax,
                         double b2,
                         double E,
                         int reg,
                         double gs,
                         const double *tau,
                         int nt,
                         const RowConfig &cfg,
                         double *Sre,
                         double *Sim) {
      const double top = (reg == 4) ? 0.5 * tmax : tmax;
      std::vector<double> T, r, h, hT, bufR;
      std::vector<cplx> bufC;
      koNodes(e0, top, cfg.knockonNPerDec, e0, T);
      lawRate(T, xi, tmax, b2, E, reg, false, cfg, r);
      h.resize(T.size());
      hT.resize(T.size());
      for (std::size_t i = 0; i < T.size(); ++i) {
        h[i] = r[i] - xi / (T[i] * T[i]);
        hT[i] = h[i] * T[i];
      }
      Panels Ph, PhT;
      makePanels(T.data(), h.data(), static_cast<int>(T.size()), Ph);
      makePanels(T.data(), hT.data(), static_cast<int>(T.size()), PhT);
      const double sh = simpsonOf(Ph, bufR), shT = simpsonOf(PhT, bufR);
      // the 1/T^2 part: exact_delta_exponent with beta^2 = 0, no spin term
      std::vector<double> cr(nt, 0.), ci(nt, 0.);
      exactDeltaBatch({ExStep{xi, e0, top, 0.0, E, gs}}, tau, nt, false, cr.data(), ci.data());
      for (int j = 0; j < nt; ++j) {
        const double a = gs * tau[j];
        double sr = cr[j], si = ci[j];
        const cplx F = filonAt(Ph, a, bufC);
        sr = sr + F.real();
        si = si + F.imag();
        sr = sr - sh;
        si = si - a * shT;
        Sre[j] += sr;
        Sim[j] += si;
      }
    }

  }  // namespace


  //--------------------------------------------------------------------------
  // cf_rows.ioni_rows -> cf_track_resolution.ioni_step_exponent(steps, 1, tau)
  // with the record's q/p-per-MeV column times the row's weight.
  //--------------------------------------------------------------------------
  void ioniRows(const double *tau,
                int nt,
                const double *rows,
                int stride,
                int n,
                const double *wq,
                const RowConfig &cfg,
                double *Sre,
                double *Sim) {
    if (rows == nullptr || n <= 0 || nt <= 0)
      return;
    if (stride < 11)
      throw cms::Exception("CvhCfExponents") << "ioniurbanv stride " << stride << " < 11";
    struct St {
      int reg;
      double gsig2, a1, e1, a2, e2, a3, e0, tmax, gam, gs, b2, et;
    };
    const bool gauge = (cfg.ioniA3Scale != 1.0 || cfg.ioniExcScale != 1.0);
    std::vector<St> st(n);
    for (int i = 0; i < n; ++i) {
      const double *r = rows + static_cast<std::size_t>(i) * stride;
      St &s = st[i];
      s.reg = static_cast<int>(r[0]);
      s.gsig2 = r[1];
      s.gam = r[9];
      s.a1 = gauge ? r[2] * cfg.ioniExcScale : r[2];
      s.a2 = gauge ? r[4] * cfg.ioniExcScale : r[4];
      s.a3 = gauge ? r[6] * cfg.ioniA3Scale : r[6];
      s.e1 = r[3] * s.gam;
      s.e2 = r[5] * s.gam;
      s.e0 = r[7] * s.gam;
      s.tmax = r[8] * s.gam;
      s.gs = (r[10] * wq[i]) * 1e-3;
      s.b2 = (stride >= 13) ? r[11] : 0.;
      s.et = (stride >= 13) ? r[12] : 0.;
    }
    std::vector<double> acc(nt), accI(nt);
    auto flush = [&]() {
      for (int j = 0; j < nt; ++j) {
        Sre[j] += acc[j];
        Sim[j] += accI[j];
      }
    };
    // regime 0: the Gaussian channel, -tau^2/2 sum gsig2 gs^2
    {
      std::vector<double> v;
      for (const St &s : st)
        if (s.reg == 0)
          v.push_back(s.gsig2 * (s.gs * s.gs));
      if (!v.empty()) {
        const double sum0 = pwSum(v.data(), static_cast<std::ptrdiff_t>(v.size()));
        for (int j = 0; j < nt; ++j)
          Sre[j] += (-0.5 * (tau[j] * tau[j])) * sum0;
      }
    }
    // the excitations, Poisson at fixed energy, centred; channel by channel
    for (int ch = 0; ch < 2; ++ch) {
      std::fill(acc.begin(), acc.end(), 0.);
      std::fill(accI.begin(), accI.end(), 0.);
      bool any = false;
      for (const St &s : st) {
        if (s.reg == 0)
          continue;
        const double aj = ch ? s.a2 : s.a1, ej = ch ? s.e2 : s.e1;
        if (!(aj > 0.) || !(ej > 0.))
          continue;
        any = true;
        const double ge = s.gs * ej;
        for (int j = 0; j < nt; ++j) {
          const double th = ge * tau[j];
          acc[j] += aj * (std::cos(th) - 1.);
          accI[j] += aj * (std::sin(th) - th);
        }
      }
      if (any)
        flush();
    }
    auto actDelta = [](const St &s) { return s.a3 > 0. && s.tmax > s.e0 && s.e0 > 0.; };
    // regime 1: a3 collisions from 1/E^2 on [e0, tmax]
    {
      std::fill(acc.begin(), acc.end(), 0.);
      std::fill(accI.begin(), accI.end(), 0.);
      bool any = false;
      for (const St &s : st) {
        if (s.reg == 0 || s.reg == 2 || s.reg == 3 || s.reg == 4 || s.reg == 5 || !actDelta(s))
          continue;
        any = true;
        const double ge0 = s.gs * s.e0, w = s.tmax / s.e0;
        for (int j = 0; j < nt; ++j) {
          const std::complex<double> d = cvhcgf::deltaTerm(ge0 * tau[j], w);
          acc[j] += s.a3 * d.real();
          accI[j] += s.a3 * d.imag();
        }
      }
      if (any)
        flush();
    }
    // regime 2/3: the exact Bethe-Bloch knock-on spectrum, xi = a3 gam, and
    // Geant4's Kokoulin correction for muons
    for (int spin = 2; spin <= 3; ++spin) {
      bool any = false;
      for (const St &s : st)
        if (s.reg == spin && actDelta(s))
          any = true;
      if (!any)
        continue;
      if (stride < 13)
        throw cms::Exception("CvhCfExponents")
            << "ioniurbanv regime " << spin << " record with stride " << stride
            << ": beta^2 and E are not in the record, and cannot be recovered without the particle mass";
      const bool half = (spin == 2);
      {
        std::vector<ExStep> ex;
        for (const St &s : st) {
          if (s.reg != spin || !actDelta(s))
            continue;
          const double tmx = (cfg.ioniTmaxScale != 1.0) ? s.tmax * cfg.ioniTmaxScale : s.tmax;
          ex.push_back(ExStep{s.a3 * s.gam, s.e0, tmx, s.b2, s.et, s.gs});
        }
        exactDeltaBatch(ex, tau, nt, half, Sre, Sim);
      }
      if (cfg.ioniKokoulin != 0.0) {
        // cf_track_resolution._kokoulin_exponent: INT_{lo}^{tmax} (f_K - 1) dN
        // (e^{iaT} - 1 - iaT) dT, lo = max(e0, tcut, 100 keV), by the centred
        // Filon-Simpson on the knock-on channel's nodes, for the muon steps
        // above 1 GeV
        std::vector<double> T, r, h, kr(nt, 0.), ki(nt, 0.);
        std::vector<cplx> bufC;
        Panels Ph;
        for (const St &s : st) {
          if (s.reg != spin || !actDelta(s))
            continue;
          const double t0 = std::max(s.e0, std::max(cfg.ioniKokoulinTcut, kKokTMin));
          const double lo0 = std::max(t0, kKokTMin);
          const double tmx = (cfg.ioniTmaxScale != 1.0) ? s.tmax * cfg.ioniTmaxScale : s.tmax;
          if (recordSpecies(s.b2, s.et, spin) != kSpeciesMu || !(lo0 < tmx) || !(s.et - kKoMuMeV > kKoKokMuMin))
            continue;
          const double xi = s.a3 * s.gam;
          koNodes(lo0, tmx, cfg.knockonNPerDec, lo0, T);
          lawRate(T, xi, tmx, s.b2, s.et, spin, false, cfg, r);
          const std::size_t nn = T.size();
          h.resize(nn);
          for (std::size_t i = 0; i < nn; ++i)
            h[i] = kokExcess(T[i], s.et) * r[i];
          makePanels(T.data(), h.data(), static_cast<int>(nn), Ph);
          for (int j = 0; j < nt; ++j) {
            const cplx F = filonCentredAt(Ph, s.gs * tau[j], bufC);
            kr[j] += F.real();
            ki[j] += F.imag();
          }
        }
        for (int j = 0; j < nt; ++j) {
          Sre[j] += cfg.ioniKokoulin * kr[j];
          Sim[j] += cfg.ioniKokoulin * ki[j];
        }
      }
    }
    // regime 4/5: the e+- Moller / Bhabha law
    for (const St &s : st) {
      if ((s.reg != 4 && s.reg != 5) || !actDelta(s))
        continue;
      if (stride < 13)
        throw cms::Exception("CvhCfExponents")
            << "ioniurbanv regime " << s.reg << " record with stride " << stride << ": beta^2 and E are not in the record";
      const double tmx = s.tmax * cfg.ioniTmaxScale;
      lawBaseExponent(s.a3 * s.gam, s.e0, tmx, s.b2, s.et, s.reg, s.gs, tau, nt, cfg, Sre, Sim);
    }
  }

  //--------------------------------------------------------------------------
  // cf_rows.ms_rows
  //--------------------------------------------------------------------------
  void msRows(const double *tau,
              int nt,
              const double *rows,
              int stride,
              int n,
              const int *rid,
              const double *wb,
              const double *frac,
              int ne,
              double scale,
              const RowConfig &cfg,
              double *S) {
    if (rows == nullptr || n <= 0 || ne <= 0 || nt <= 0)
      return;
    if (stride < 8)
      throw cms::Exception("CvhCfExponents") << "msmoliv stride " << stride << " < 8";
    std::vector<double> loc(nt, 0.);
    bool whole = (ne == n);
    for (int e = 0; whole && e < ne; ++e)
      whole = (frac[e] == 1.0 && wb[e] == wb[0] && rid[e] == e);
    if (whole) {
      // one weight for the whole record set (a fit block): the vectorised
      // path of the reference
      if (!(wb[0] > 0.))
        return;
      std::vector<MsStep> st;
      buildMsSteps(rows, stride, n, cfg, st);
      msExponentImpl(st, wb, 1, tau, nt, loc.data());
      for (int j = 0; j < nt; ++j)
        S[j] += scale * loc[j];
      return;
    }
    std::vector<MsStep> one(1);
    for (int e = 0; e < ne; ++e) {
      if (wb[e] <= 0.)
        continue;
      if (!msStepOf(rows + static_cast<std::size_t>(rid[e]) * stride, stride, frac[e], cfg, one[0]))
        continue;
      std::fill(loc.begin(), loc.end(), 0.);
      msExponentImpl(one, wb + e, 1, tau, nt, loc.data());
      for (int j = 0; j < nt; ++j)
        S[j] += scale * loc[j];
    }
  }

  //--------------------------------------------------------------------------
  // cf_brems_exact.refine_spectra
  //--------------------------------------------------------------------------
  void refineSpectra(const double *recs,
                     int rstride,
                     int n,
                     const double *spec,
                     const double *vg,
                     int nv,
                     int nsub,
                     std::vector<double> &vf,
                     std::vector<double> &specf) {
    const int nf = (nv - 1) * nsub + 1;
    std::vector<double> lv(nv), idx(nv);
    for (int i = 0; i < nv; ++i) {
      lv[i] = std::log(vg[i]);
      idx[i] = i;
    }
    // np.linspace(0, nv - 1, nf): i * step + 0 with the end forced
    vf.resize(nf);
    std::vector<double> lvf(nf);
    {
      const double step = (nv - 1.0) / static_cast<double>(nf - 1);
      for (int i = 0; i < nf; ++i) {
        const double u = (i == nf - 1) ? (nv - 1.0) : static_cast<double>(i) * step + 0.0;
        vf[i] = std::exp(npInterp(u, idx.data(), lv.data(), nv));
      }
      // the exported nodes exactly: exp(log(v)) is not always v, and the
      // spectrum's support below is decided by comparing against them
      for (int i = 0; i < nv; ++i)
        vf[static_cast<std::size_t>(i) * nsub] = vg[i];
      for (int i = 0; i < nf; ++i)
        lvf[i] = std::log(vf[i]);
    }
    specf.assign(static_cast<std::size_t>(n) * 2 * nf, 0.);
    std::vector<double> xp, fp, f(nf);
    std::vector<int> pos;
    for (int s = 0; s < n; ++s) {
      const double *rec = recs + static_cast<std::size_t>(s) * rstride;
      const double *sp = spec + static_cast<std::size_t>(s) * 2 * nv;
      const double E = rec[3], p = rec[4];
      const double vmax = (E > 0.0) ? (E - speciesMass(E, p)) / E : 0.0;
      for (int h = 0; h < 2; ++h) {
        const double *y = sp + h * nv;
        pos.clear();
        for (int i = 0; i < nv; ++i)
          if (y[i] > 0.0)
            pos.push_back(i);
        if (pos.size() < 2)
          continue;
        const int i0 = pos.front(), iL = pos.back(), i1 = pos[pos.size() - 2];
        xp.resize(pos.size());
        fp.resize(pos.size());
        for (std::size_t k = 0; k < pos.size(); ++k) {
          xp[k] = lv[pos[k]];
          fp[k] = std::log(y[pos[k]]);
        }
        std::fill(f.begin(), f.end(), 0.);
        for (int i = 0; i < nf; ++i)
          if (vf[i] >= vg[i0] && vf[i] <= vg[iL])
            f[i] = std::exp(npInterp(lvf[i], xp.data(), fp.data(), static_cast<int>(xp.size())));
        if (vmax > vg[iL]) {
          const double sl = std::log(y[iL] / y[i1]) / (lv[iL] - lv[i1]);
          for (int i = 0; i < nf; ++i)
            if (vf[i] > vg[iL] && vf[i] <= vmax)
              f[i] = y[iL] * std::pow(vf[i] / vg[iL], sl);
          int last = -1;
          for (int i = 0; i < nf; ++i)
            if (vf[i] <= vmax)
              last = i;
          if (last + 1 < nf) {
            const double a_ = vf[last], b_ = vf[last + 1];
            const double fa = f[last];
            double A;
            if (std::fabs(sl + 1.0) > 1e-9)
              A = fa * a_ / (sl + 1.0) * (std::pow(vmax / a_, sl + 1.0) - 1.0);
            else
              A = fa * a_ * std::log(vmax / a_);
            f[last + 1] = std::max(2.0 * A / (b_ - a_) - fa, 0.0);
          }
        }
        std::copy(f.begin(), f.end(), specf.begin() + (static_cast<std::size_t>(s) * 2 + h) * nf);
      }
    }
  }

  //--------------------------------------------------------------------------
  // cf_rows.rad_rows -> cf_brems_exact.rad_exponent
  //--------------------------------------------------------------------------
  namespace {
    constexpr double kRadMeGeV = 0.51099895e-3;
    constexpr double kRadMuGeV = 0.1056583745;

    // e+- or not (cf_brems_exact.is_epm), from a mass rebuilt from float32
    // records: good to a few MeV, and the lightest heavy species is at 106 MeV
    inline bool radIsEpm(double M) { return M < 0.01; }

    // E[J0(c w)] - 1 of the photon-angle law (cf_brems_exact.rad_angle_cfm1),
    // c = (projected weight) x theta / w: G4ModifiedTsai for e+-,
    // G4ModifiedMephi (c K1(c)) for every heavier species (muBrems and hBrems
    // share it).
    inline double radAngleCfm1(double c, bool epm) {
      if (epm) {
        const double a1 = 1.6, a2 = 1.6 / 3.0;
        const double u1 = a1 * c, u2 = a2 * c;
        return 0.25 * std::expm1(-1.5 * std::log1p(u1 * u1)) + 0.75 * std::expm1(-1.5 * std::log1p(u2 * u2));
      }
      return ((c > 0.0) ? c * besselK1(std::max(c, 1e-300)) : 1.0) - 1.0;
    }
  }  // namespace

  double speciesMass(double E, double p) {
    const double M = std::sqrt(std::max(E * E - p * p, 0.0));
    int best = 0;
    double bd = std::fabs(kSpeciesMassesGeV[0] - M);
    for (int i = 1; i < 5; ++i) {
      const double d = std::fabs(kSpeciesMassesGeV[i] - M);
      if (d < bd) {
        bd = d;
        best = i;
      }
    }
    return kSpeciesMassesGeV[best];
  }

  void radRows(const double *tau,
               int nt,
               const double *recs,
               int rstride,
               int n,
               const double *spec,
               const double *vg,
               int nv,
               const int *rid,
               const double *wq,
               const double *wb,
               const double *frac,
               int ne,
               bool exactQop,
               double *Sre,
               double *Sim) {
    if (recs == nullptr || spec == nullptr || vg == nullptr || n <= 0 || ne <= 0 || nv < 2 || nt <= 0)
      return;
    if (rstride < 11)
      throw cms::Exception("CvhCfExponents") << "radiative record stride " << rstride << " < 11";
    std::vector<double> wtrap(nv, 0.), dv(nv - 1);
    for (int i = 0; i + 1 < nv; ++i)
      dv[i] = vg[i + 1] - vg[i];
    for (int i = 0; i + 1 < nv; ++i)
      wtrap[i] += 0.5 * dv[i];
    for (int i = 0; i + 1 < nv; ++i)
      wtrap[i + 1] += 0.5 * dv[i];
    std::vector<double> dNb(nv), dNp(nv), dNdv(nv), tz(nv - 1), T(nv), X(nv), c(nv), qb(nv), q(nv), a(nt), re(nv),
        im(nv);
    for (int e = 0; e < ne; ++e) {
      const double *rec = recs + static_cast<std::size_t>(rid[e]) * rstride;
      const double *sp = spec + static_cast<std::size_t>(rid[e]) * 2 * nv;
      const double w = wq[e];
      const double E = rec[3];
      const double L = rec[6] * frac[e];
      // step_spectrum_parts: each process normalised to its own mean loss
      for (int proc = 0; proc < 2; ++proc) {
        double *out = proc ? dNp.data() : dNb.data();
        std::fill(out, out + nv, 0.);
        const double dE = rec[8 + proc] * L;
        if (dE <= 0.)
          continue;
        const double *shape = sp + proc * nv;
        for (int i = 0; i + 1 < nv; ++i) {
          const double y0 = (shape[i] * vg[i]) * E, y1 = (shape[i + 1] * vg[i + 1]) * E;
          tz[i] = (dv[i] * (y1 + y0)) / 2.0;
        }
        const double norm = pwSum(tz.data(), nv - 1);
        if (norm > 0.)
          for (int i = 0; i < nv; ++i)
            out[i] = shape[i] * (dE / norm);
      }
      bool any = false, anyB = false;
      for (int i = 0; i < nv; ++i) {
        dNdv[i] = dNb[i] + dNp[i];
        any = any || dNdv[i] > 0.;
        anyB = anyB || dNb[i] > 0.;
      }
      if (!any)
        continue;
      const double bw = wb ? wb[e] : 0.0;
      const double p = rec[4];
      const double cs = rec[10];
      for (int j = 0; j < nt; ++j)
        a[j] = tau[j] * cs * w;
      if (bw != 0.0 && anyB) {
        // the recoil: theta = (k/p') theta_gamma, theta_gamma = u M/E
        const double M = speciesMass(E, p);
        const bool epm = radIsEpm(M);
        for (int i = 0; i < nv; ++i) {
          const double Ti = vg[i] * E;
          const double d = E - Ti;
          const double pp = std::sqrt(std::max(d * d - M * M, (1e-3 * p) * (1e-3 * p)));
          c[i] = Ti / pp * M / E;
          X[i] = exactQop ? p * p * Ti * (2.0 * E - Ti) / (E * pp * (p + pp)) : Ti;
          qb[i] = wtrap[i] * dNb[i];
        }
        const double *tb = tau;
        for (int j = 0; j < nt; ++j) {
          const double tbw = tb[j] * bw;
          double sr = 0., si = 0.;
          for (int i = 0; i < nv; ++i) {
            const double g = radAngleCfm1(tbw * c[i], epm);
            const double ph = a[j] * X[i];
            sr += (std::cos(ph) * g) * qb[i];
            si += (std::sin(ph) * g) * qb[i];
          }
          Sre[j] += sr;
          Sim[j] += si;
        }
      }
      for (int i = 0; i < nv; ++i) {
        T[i] = vg[i] * E;
        q[i] = wtrap[i] * dNdv[i];
      }
      if (exactQop) {
        // p' floored at 1e-3 p: a loss that leaves the primary below that
        // never reaches a plane, and the floor keeps T_eff finite
        const double M = speciesMass(E, p);
        for (int i = 0; i < nv; ++i) {
          const double d = E - T[i];
          const double pp = std::sqrt(std::max(d * d - M * M, (1e-3 * p) * (1e-3 * p)));
          X[i] = p * p * T[i] * (2.0 * E - T[i]) / (E * pp * (p + pp));
        }
        for (int j = 0; j < nt; ++j) {
          double sr = 0., si = 0.;
          for (int i = 0; i < nv; ++i) {
            const double x = a[j] * T[i], xe = a[j] * X[i];
            double r, m;
            if (std::fabs(xe) < 1e-4) {
              r = -0.5 * (xe * xe);
              m = (xe - x) - std::pow(xe, 3.0) / 6.;
            } else {
              const double sh = std::sin(0.5 * xe);
              r = -2.0 * sh * sh;
              m = std::sin(xe) - x;
            }
            sr += r * q[i];
            si += m * q[i];
          }
          Sre[j] += sr;
          Sim[j] += si;
        }
      } else {
        for (int j = 0; j < nt; ++j) {
          double sr = 0., si = 0.;
          for (int i = 0; i < nv; ++i) {
            const double x = a[j] * T[i];
            double r, m;
            if (std::fabs(x) < 1e-4) {
              r = -0.5 * (x * x);
              m = std::pow(x, 3.0) / 6.;
            } else {
              const double sh = std::sin(0.5 * x);
              r = -2.0 * sh * sh;
              m = std::sin(x) - x;
            }
            sr += r * q[i];
            si += m * q[i];
          }
          Sre[j] += sr;
          Sim[j] += si;
        }
      }
    }
  }

  //--------------------------------------------------------------------------
  // cf_rows.knockon_rows -> cf_knockon.knockon_rows / step_correction.
  // `mask` (null: every row) selects the rows evaluated; the thin-step node
  // density is set by the largest xi of ALL `n` rows, so a block split by
  // material group evaluates each group's rows exactly as the whole block
  // would.
  //--------------------------------------------------------------------------
  namespace {
    void knockonRowsMasked(const double *tau,
                           int nt,
                           const double *rows,
                           int stride,
                           int n,
                           const int *rid,
                           const double *wq,
                           const double *wb,
                           const double *frac,
                           int ne,
                           KnockonPart part,
                           const RowConfig &cfg,
                           const char *mask,
                           double *Sre,
                           double *Sim) {
    if (!cfg.knockonActive() || rows == nullptr || n <= 0 || ne <= 0 || nt <= 0)
      return;
    bool anyReg = false;
    for (int i = 0; i < n; ++i) {
      const int reg = static_cast<int>(rows[static_cast<std::size_t>(i) * stride]);
      anyReg = anyReg || (reg >= 2 && reg <= 5);
    }
    // regime-1 records carry no exact knock-on spectrum (and no beta^2, E)
    if (!anyReg)
      return;
    if (stride < 13)
      throw cms::Exception("CvhCfExponents") << "knock-on channel needs the stride-13 exact-delta record";
    const bool joint = cfg.knockonJoint, exact = cfg.qopExact;
    if (part == KnockonPart::kMap && !exact)
      return;  // INT dN [e^{iaT} - e^{iaT}] is identically zero
    if (part == KnockonPart::kJoint && !joint)
      return;
    const bool kokOn = cfg.ioniKokoulin != 0.0;
    struct R {
      int reg;
      double xi, e0, tmax, b2, E, g, p;
      bool mu, act;
    };
    std::vector<R> rr(n);
    double xmax = 0.;
    bool anyAct = false;
    for (int i = 0; i < n; ++i) {
      const double *r = rows + static_cast<std::size_t>(i) * stride;
      R &s = rr[i];
      s.reg = static_cast<int>(r[0]);
      const double gam = r[9];
      s.xi = r[6] * gam * cfg.ioniA3Scale;
      s.e0 = r[7] * gam;
      s.tmax = r[8] * gam * cfg.ioniTmaxScale;
      s.b2 = r[11];
      s.E = r[12];
      s.g = r[10] * 1e-3;
      s.p = s.E * std::sqrt(s.b2);
      s.mu = recordSpecies(s.b2, s.E, s.reg) == kSpeciesMu && (s.E - kKoMuMeV > kKoKokMuMin);
      s.act = (s.reg >= 2 && s.reg <= 5) && s.xi > 0. && s.tmax > s.e0 && s.e0 > 0.;
      if (s.act) {
        xmax = anyAct ? std::max(xmax, s.xi) : s.xi;
        anyAct = true;
      }
    }
    // the entries of each record, in their own order
    std::vector<std::vector<int>> ents(n);
    for (int e = 0; e < ne; ++e)
      if (rid[e] >= 0 && rid[e] < n)
        ents[rid[e]].push_back(e);
    std::vector<double> T, dN, X, dX, th, amp, ampX;
    std::vector<cplx> bufC;
    Panels PB, PA, Pjx;
    for (int s = 0; s < n; ++s) {
      const R &k = rr[s];
      if (!k.act || ents[s].empty() || (mask != nullptr && !mask[s]))
        continue;
      const bool thin = k.xi < cfg.knockonThin * xmax;
      koNodes(k.e0,
              lawTop(k.reg, k.tmax, k.E, k.p, cfg),
              thin ? cfg.knockonNPerDecThin : cfg.knockonNPerDec,
              cfg.knockonTcut,
              T);
      const std::size_t nn = T.size();
      lawRate(T, k.xi, k.tmax, k.b2, k.E, k.reg, kokOn && k.mu, cfg, dN);
      X.resize(nn);
      dX.resize(nn);
      th.resize(nn);
      for (std::size_t i = 0; i < nn; ++i) {
        const double pp = ppOfT(T[i], k.E, k.p);
        if (exact) {
          X[i] = k.p * k.p * T[i] * (2.0 * k.E - T[i]) / (k.E * pp * (k.p + pp));
          dX[i] = (std::pow(k.p, 3.0) / k.E) * (k.E - T[i]) / std::pow(pp, 3.0);
        } else {
          X[i] = T[i];
          dX[i] = 1.0;
        }
        if (joint)
          th[i] = std::sqrt(std::max(2.0 * kKoMeMeV * T[i] * (1.0 - T[i] / k.tmax), 0.0)) / pp;
      }
      // the tau-independent amplitudes: dN in T, dN / dX in X (map, all)
      makePanels(T.data(), dN.data(), static_cast<int>(nn), PB);
      if (exact && part != KnockonPart::kJoint) {
        ampX.resize(nn);
        for (std::size_t i = 0; i < nn; ++i)
          ampX[i] = dN[i] / dX[i];
        makePanels(X.data(), ampX.data(), static_cast<int>(nn), PA);
      }
      amp.resize(nn);
      ampX.resize(nn);
      for (int e : ents[s]) {
        const double al = wq[e] * k.g;
        const double be = wb[e];
        if (al == 0.0 && (be == 0.0 || !joint))
          continue;
        const double fr = frac[e];
        for (int j = 0; j < nt; ++j) {
          const double ta = tau[j] * al, tb = tau[j] * be;
          double vr = 0., vi = 0.;
          if (part == KnockonPart::kMap) {
            const cplx A = filonAt(PA, ta, bufC), B = filonAt(PB, ta, bufC);
            vr = A.real() - B.real();
            vi = A.imag() - B.imag();
          } else if (part == KnockonPart::kJoint) {
            // INT dN e^{iaX} (J - 1): the collision's deflection, jointly
            // with its loss
            for (std::size_t i = 0; i < nn; ++i) {
              amp[i] = dN[i] * (besselJ0(tb * th[i]) - 1.0);
              ampX[i] = amp[i] / dX[i];
            }
            makePanels(exact ? X.data() : T.data(), exact ? ampX.data() : amp.data(), static_cast<int>(nn), Pjx);
            const cplx A = filonAt(Pjx, ta, bufC);
            vr = A.real();
            vi = A.imag();
          } else {
            // INT dN [e^{iaX} J - e^{iaT}]
            const cplx B = filonAt(PB, ta, bufC);
            for (std::size_t i = 0; i < nn; ++i) {
              const double J = joint ? besselJ0(tb * th[i]) : 1.0;
              ampX[i] = exact ? (dN[i] * J) / dX[i] : dN[i] * J;
            }
            makePanels(exact ? X.data() : T.data(), ampX.data(), static_cast<int>(nn), Pjx);
            const cplx A = filonAt(Pjx, ta, bufC);
            vr = A.real() - B.real();
            vi = A.imag() - B.imag();
          }
          Sre[j] += vr * fr;
          Sim[j] += vi * fr;
        }
      }
    }
  }

  }  // namespace

  void knockonRows(const double *tau,
                   int nt,
                   const double *rows,
                   int stride,
                   int n,
                   const int *rid,
                   const double *wq,
                   const double *wb,
                   const double *frac,
                   int ne,
                   KnockonPart part,
                   const RowConfig &cfg,
                   double *Sre,
                   double *Sim) {
    knockonRowsMasked(tau, nt, rows, stride, n, rid, wq, wb, frac, ne, part, cfg, nullptr, Sre, Sim);
  }

  //==========================================================================
  // 7. THE FIT-LEVEL FAMILIES  (cf_rows.fit_families, with the pooling of
  //    cf_track_resolution.extract / cf_mass_likelihood.build_pairs_tt)
  //
  // THE POOLING: for each material family in (10, 11), for each DISTINCT
  // global parameter index among the resolution entries of that family, pool
  // v_b over the entries, gather the step rows whose index VALUE matches, and
  // form the block's scalar standardized weight -- sqrt(v_b / sq2) / sigma,
  // sq2 the fit's own variance of the block (Rossi's thp2 sum for scattering,
  // `ioniSq2` for ionisation), signed for ionisation by the functional's
  // `ioniSign` and the block's variance-weighted influence sign.  Ascending
  // index order is `np.unique`'s.
  //
  // THE FAMILIES, from the ROW FUNCTIONS at those weights, one entry per
  // record:
  //   Sms   every scattering block at |w|;
  //   Sio   every ionisation block at w;
  //   Skx   the knock-on exact map and
  //   Skj   the knock-on joint piece, on the same ionisation rows at w, with
  //         the angular weight beta_g of the block: the xg-weighted rms of the
  //         scattering weights of the Geant4 steps it pairs with -- the
  //         radiative rows are one per step, PARALLEL to `msmoliv`, and carry
  //         their leg's ionisation block index (cf_rows.fit_pairing);
  //   Srad  the radiative rows of each ionisation block (joined on the index
  //         VALUE) at w, on the spectra refined once per track, each row's
  //         recoil at the scattering weight of its own step.
  // plus the nuclear-elastic family of hadrons (section 5).
  //
  // THE PER-GROUP SPLIT.  Every family is a sum over rows, so a family split
  // by the rows' material group is exact; when a functional asks for the
  // split, its flat families are formed AS the sum of the group parts (in
  // ascending group order), so `sum_g groups[i].S == S` bit for bit.  The
  // knock-on's thin-step node density is a property of the whole block (its
  // largest xi), so the split evaluates each group's rows with the block's.
  //
  // THE MULTI-FUNCTIONAL PASS.  The makers form several linear functionals of
  // one fit (the candidate mass, the vertex DCA, the beam-line pulls); they
  // share every block and every record and differ only in the per-block
  // weights.  The weight-independent work -- the pooling, the widened
  // records, the refined spectra, the Moliere step parameters -- is done once;
  // the families are then evaluated per functional.
  //==========================================================================
  namespace {

    // The group entry of `v` (kept ascending in `group`), created on first use.
    GroupExponents &groupEntry(std::vector<GroupExponents> &v, int g) {
      auto it = std::lower_bound(
          v.begin(), v.end(), g, [](const GroupExponents &a, int b) { return a.group < b; });
      if (it != v.end() && it->group == g)
        return *it;
      GroupExponents e;
      e.group = g;
      return *v.insert(it, e);
    }

    // Per-STEP ionisation variance contribution with the per-leg CGF scale
    // applied, so that a material group's share of `ioniSq2` is a plain sum
    // over its own steps.  Same association, left to right, as `ioniSq2`.
    void ioniStepVar(const float *rows, int stride, int n, const float *qsc, int nqsc, std::vector<double> &w2) {
      w2.assign(static_cast<std::size_t>(std::max(n, 0)), 0.);
      if (rows == nullptr || n <= 0)
        return;
      for (int i = 0; i < n; ++i) {
        const float *r = rows + static_cast<std::size_t>(i) * stride;
        const double gq = static_cast<double>(r[10]) * 1e-3;
        w2[i] = static_cast<double>(r[1]) * gq * gq;
      }
      if (qsc == nullptr || nqsc <= 0)
        return;
      long long tot = 0;
      for (int i = 0; i < nqsc; ++i)
        tot += static_cast<long long>(qsc[2 * i + 1]);
      if (tot != static_cast<long long>(n)) {
        double msc = 0.;
        for (int i = 0; i < nqsc; ++i)
          msc += qsc[2 * i];
        const double f = msc / nqsc;
        for (double &x : w2)
          x *= f;
        return;
      }
      int off = 0;
      for (int i = 0; i < nqsc; ++i) {
        const int ns = static_cast<int>(qsc[2 * i + 1]);
        for (int k = 0; k < ns; ++k)
          w2[off + k] *= static_cast<double>(qsc[2 * i]);
        off += ns;
      }
    }

    // The material group of record `i` (-1 when the records carry none).
    inline int groupOf(const StepRows &r, int i) {
      return (r.groupCol >= 0 && r.groupCol < r.stride)
                 ? static_cast<int>(r.v[static_cast<std::size_t>(i) * r.stride + r.groupCol])
                 : -1;
    }

    // The rows `rows` of a record set partitioned by material group,
    // ascending; one part (group -1 or the rows' single group) when `split`
    // is false or the rows carry no group.
    void partition(const StepRows &r,
                   const std::vector<int> &rows,
                   bool split,
                   std::vector<std::pair<int, std::vector<int>>> &parts) {
      parts.clear();
      if (!split) {
        parts.emplace_back(rows.empty() ? -1 : groupOf(r, rows.front()), rows);
        return;
      }
      std::map<int, std::vector<int>> m;
      for (int i : rows)
        m[groupOf(r, i)].push_back(i);
      for (auto &kv : m)
        parts.emplace_back(kv.first, std::move(kv.second));
    }

    // The contiguous double copy of records `idx` of a widened record set.
    void gather(const std::vector<double> &src, int stride, const std::vector<int> &idx, std::vector<double> &dst) {
      dst.resize(idx.size() * static_cast<std::size_t>(stride));
      for (std::size_t i = 0; i < idx.size(); ++i)
        std::memcpy(&dst[i * stride], &src[static_cast<std::size_t>(idx[i]) * stride], stride * sizeof(double));
    }

    // One pooled block: its global index, its record rows, and per
    // functional its weight (standardized; signed for ionisation), its share
    // of the standardized variance under the fit's Q, and whether the
    // functional carries it.
    struct Block {
      unsigned int g = 0;
      std::vector<int> rows;
      std::vector<double> w, vq;
      std::vector<char> act;
    };

    // S += src over kNTau (one family of one group or the flat result)
    inline void addTo(std::array<double, kNTau> &dst, const double *src) {
      for (int j = 0; j < kNTau; ++j)
        dst[j] += src[j];
    }

    void trackExponentsImpl(const TrackInput *in, int nf, TrackResult *out) {
      for (int k = 0; k < nf; ++k) {
        out[k].ok = true;
        out[k].vgauss = 0.;
        out[k].nblockms = out[k].nblockioni = out[k].npooled = 0;
        out[k].S.clear();
        out[k].groups.clear();
        out[k].nucel = false;
        out[k].nuc.clear();
        out[k].nucGroups.clear();
      }
      if (nf <= 0)
        return;

      // A functional with no standardization or no resolution entries is
      // dropped; the others carry on.
      std::vector<char> alive(nf, 1);
      int nalive = 0, kref = -1;
      for (int k = 0; k < nf; ++k) {
        if (!(in[k].sigma > 0.) || in[k].nres <= 0) {
          out[k].ok = false;
          alive[k] = 0;
        } else {
          if (kref < 0)
            kref = k;
          ++nalive;
        }
      }
      if (nalive == 0)
        return;
      // The step records, the (global index, family) arrays and the model
      // configuration are SHARED by contract; read from the first usable
      // functional.
      const TrackInput &sh = in[kref];
      const RowConfig &cfg = sh.rowConfig ? *sh.rowConfig : productionRowConfig();
      const double *tau = tauGrid();

      // The GAUSSIAN families: 8/9 the hit blocks, 16 the two-track maker's
      // beam-line block (a fixed covariance, no Landau channel), so that
      // `sum_g (vqms + vqio) + vgauss/sigma^2 == 1` stays exact.
      int nresmax = 0;
      for (int k = 0; k < nf; ++k) {
        if (!alive[k])
          continue;
        nresmax = std::max(nresmax, in[k].nres);
        for (int i = 0; i < in[k].nres; ++i)
          if (in[k].resfamily[i] == 8 || in[k].resfamily[i] == 9 || in[k].resfamily[i] == 16)
            out[k].vgauss += in[k].resvarv[i];
      }

      //---------------------------------------------------------------- pool
      std::vector<Block> msb, iob;
      {
        std::vector<unsigned int> gs;
        std::vector<double> vpool(nf), vsig(nf);
        std::vector<int> nent(nf);
        std::vector<float> blk, qs;
        for (int fam = 10; fam <= 11; ++fam) {
          const StepRows &rows = (fam == 10) ? sh.ms : sh.ioni;
          gs.clear();
          for (int i = 0; i < nresmax; ++i)
            if (sh.resfamily[i] == fam)
              gs.push_back(sh.resglobidx[i]);
          std::sort(gs.begin(), gs.end());
          gs.erase(std::unique(gs.begin(), gs.end()), gs.end());
          for (unsigned int g : gs) {
            std::vector<int> kact;
            for (int k = 0; k < nf; ++k) {
              vpool[k] = vsig[k] = 0.;
              nent[k] = 0;
              if (!alive[k])
                continue;
              const TrackInput &I = in[k];
              for (int i = 0; i < I.nres; ++i)
                if (I.resfamily[i] == fam && I.resglobidx[i] == g) {
                  vpool[k] += I.resvarv[i];
                  vsig[k] += (I.ressgn != nullptr ? double(I.ressgn[i]) : 1.) * I.resvarv[i];
                  ++nent[k];
                }
              if (vpool[k] > 0.)
                kact.push_back(k);
            }
            if (kact.empty())
              continue;
            Block b;
            b.g = g;
            for (int i = 0; i < rows.n; ++i)
              if (rows.idx[i] == g)
                b.rows.push_back(i);
            if (b.rows.empty()) {
              // A registered block with no step rows: the offline extractor
              // drops such a track, so the functionals that carry the block
              // are dropped (`ok` false) rather than kept with a model that
              // misses it.
              for (int k : kact) {
                out[k].ok = false;
                alive[k] = 0;
                --nalive;
              }
              if (nalive == 0)
                return;
              continue;
            }
            for (int k : kact)
              if (nent[k] > 1)
                ++out[k].npooled;
            b.w.assign(nf, 0.);
            b.vq.assign(nf, 0.);
            b.act.assign(nf, 0);
            if (fam == 10) {
              for (int k : kact)
                ++out[k].nblockms;
              double sq2 = 0.;
              for (int i : b.rows)
                sq2 += rows.v[static_cast<std::size_t>(i) * rows.stride + 5];
              if (!(sq2 > 0.))
                continue;
              for (int k : kact) {
                b.w[k] = std::sqrt(vpool[k] / sq2) / in[k].sigma;
                b.vq[k] = vpool[k] / (in[k].sigma * in[k].sigma);
                b.act[k] = 1;
              }
              msb.push_back(std::move(b));
            } else {
              for (int k : kact)
                ++out[k].nblockioni;
              blk.resize(b.rows.size() * static_cast<std::size_t>(rows.stride));
              for (std::size_t i = 0; i < b.rows.size(); ++i)
                std::memcpy(&blk[i * rows.stride],
                            rows.v + static_cast<std::size_t>(b.rows[i]) * rows.stride,
                            rows.stride * sizeof(float));
              qs.clear();
              for (int i = 0; i < sh.qsc.n; ++i)
                if (sh.qsc.idx[i] == g) {
                  qs.push_back(sh.qsc.v[static_cast<std::size_t>(i) * sh.qsc.stride]);
                  qs.push_back(sh.qsc.v[static_cast<std::size_t>(i) * sh.qsc.stride + 1]);
                }
              const double sq2 = ioniSq2(blk.data(),
                                         rows.stride,
                                         static_cast<int>(b.rows.size()),
                                         qs.empty() ? nullptr : qs.data(),
                                         static_cast<int>(qs.size() / 2));
              if (!(sq2 > 0.))
                continue;
              for (int k : kact) {
                const double sgnblk = (vsig[k] < 0.) ? -1. : 1.;
                b.w[k] = in[k].ioniSign * sgnblk * (std::sqrt(vpool[k] / sq2) / in[k].sigma);
                b.vq[k] = vpool[k] / (in[k].sigma * in[k].sigma);
                b.act[k] = 1;
              }
              iob.push_back(std::move(b));
            }
          }
        }
      }

      bool anyGroups = false, anyNucel = false, anyKnock = false;
      for (int k = 0; k < nf; ++k) {
        if (!alive[k])
          continue;
        anyGroups = anyGroups || in[k].wantGroups;
        anyNucel = anyNucel || (in[k].wantNucel && in[k].nucel != nullptr);
        anyKnock = anyKnock || in[k].wantKnockon;
      }
      const bool knockMap = anyKnock && cfg.knockonActive() && cfg.qopExact;
      const bool knockJoint = anyKnock && cfg.knockonActive() && cfg.knockonJoint;

      //---------------------------------------------- the widened records
      const StepRows &MS = sh.ms, &IO = sh.ioni, &RD = sh.rad;
      std::vector<double> msd(MS.v, MS.v + static_cast<std::size_t>(std::max(MS.n, 0)) * MS.stride);
      std::vector<double> iod(IO.v, IO.v + static_cast<std::size_t>(std::max(IO.n, 0)) * IO.stride);
      const bool haveRad = RD.n > 0 && RD.v != nullptr && sh.radspec != nullptr && sh.radvgrid != nullptr &&
                           sh.radnv > 1;
      std::vector<double> radd, vf, specf;
      int nvf = 0;
      if (haveRad) {
        if (RD.stride < 11)
          throw cms::Exception("CvhCfExponents") << "radiative record stride " << RD.stride << " < 11";
        // the radiative rows pair each Geant4 step with its scattering row
        // (fit_pairing): one per step, PARALLEL to `msmoliv`
        if (RD.n != MS.n)
          throw cms::Exception("CvhCfExponents") << "the radiative rows (" << RD.n
                                                 << ") are not parallel to the scattering rows (" << MS.n << ")";
        for (int i = 0; i < RD.n; ++i)
          if (RD.v[static_cast<std::size_t>(i) * RD.stride + 4] != MS.v[static_cast<std::size_t>(i) * MS.stride + 3])
            throw cms::Exception("CvhCfExponents")
                << "radiative row " << i << " (p = " << RD.v[static_cast<std::size_t>(i) * RD.stride + 4]
                << ") is not the scattering row's step (p = " << MS.v[static_cast<std::size_t>(i) * MS.stride + 3]
                << ")";
        radd.assign(RD.v, RD.v + static_cast<std::size_t>(RD.n) * RD.stride);
        std::vector<double> specd(sh.radspec, sh.radspec + static_cast<std::size_t>(RD.n) * 2 * sh.radnv);
        std::vector<double> vgd(sh.radvgrid, sh.radvgrid + sh.radnv);
        refineSpectra(radd.data(), RD.stride, RD.n, specd.data(), vgd.data(), sh.radnv, cfg.radNsub, vf, specf);
        nvf = static_cast<int>(vf.size());
      }

      std::vector<std::pair<int, std::vector<int>>> parts;
      std::vector<double> rowsd, specsub, buf, re(kNTau), im(kNTau), wq, wb, ones;
      std::vector<int> rid;

      //------------------------------------------------------------ Sms
      // All functionals of a block at once: the Moliere step parameters and
      // the shape rows they interpolate are weight independent.
      {
        std::vector<MsStep> st;
        std::vector<int> kk;
        std::vector<double> wk;
        for (const Block &b : msb) {
          kk.clear();
          wk.clear();
          for (int k = 0; k < nf; ++k)
            if (alive[k] && b.act[k]) {
              kk.push_back(k);
              wk.push_back(std::fabs(b.w[k]));
            }
          if (kk.empty())
            continue;
          const int na = static_cast<int>(kk.size());
          partition(MS, b.rows, anyGroups, parts);
          double sq2 = 0.;
          for (int i : b.rows)
            sq2 += MS.v[static_cast<std::size_t>(i) * MS.stride + 5];
          for (const auto &pt : parts) {
            gather(msd, MS.stride, pt.second, rowsd);
            buildMsSteps(rowsd.data(), MS.stride, static_cast<int>(pt.second.size()), cfg, st);
            buf.assign(static_cast<std::size_t>(na) * kNTau, 0.);
            msExponentImpl(st, wk.data(), na, tau, kNTau, buf.data());
            double sq2g = 0.;
            for (int i : pt.second)
              sq2g += MS.v[static_cast<std::size_t>(i) * MS.stride + 5];
            for (int a = 0; a < na; ++a) {
              const int k = kk[a];
              const double *Sa = buf.data() + static_cast<std::size_t>(a) * kNTau;
              addTo(out[k].S.ms, Sa);
              if (in[k].wantGroups) {
                GroupExponents &ge = groupEntry(out[k].groups, pt.first);
                addTo(ge.S.ms, Sa);
                ge.vqms += (parts.size() == 1) ? b.vq[k] : b.vq[k] * (sq2g / sq2);
              }
            }
          }
        }
      }

      //------------------------------------------ per functional: the rest
      for (int k = 0; k < nf; ++k) {
        if (!alive[k])
          continue;
        const TrackInput &I = in[k];
        TrackResult &R = out[k];
        const bool grp = I.wantGroups;
        // the scattering and ionisation block weights of this functional,
        // ascending in global index (the pairing, the nuclear-elastic family)
        std::vector<std::pair<unsigned int, double>> wms, wio;
        for (const Block &b : msb)
          if (b.act[k])
            wms.emplace_back(b.g, std::fabs(b.w[k]));
        for (const Block &b : iob)
          if (b.act[k])
            wio.emplace_back(b.g, b.w[k]);
        auto look = [](const std::vector<std::pair<unsigned int, double>> &v, unsigned int g) {
          auto it = std::lower_bound(
              v.begin(), v.end(), g, [](const std::pair<unsigned int, double> &a, unsigned int b) { return a.first < b; });
          return (it != v.end() && it->first == g) ? it->second : 0.;
        };

        // cf_rows.fit_pairing: each radiative row's angular weight (its
        // step's scattering block) and each ionisation block's beta, the
        // xg-weighted rms of those over the block's steps
        std::vector<double> wbRad;
        std::map<unsigned int, double> beta;
        if (haveRad) {
          wbRad.resize(RD.n);
          for (int s = 0; s < RD.n; ++s)
            wbRad[s] = look(wms, MS.idx[s]);
          std::map<unsigned int, std::pair<std::vector<double>, std::vector<double>>> acc;
          for (int s = 0; s < RD.n; ++s) {
            auto &a = acc[RD.idx[s]];
            const double xg = msd[static_cast<std::size_t>(s) * MS.stride + 2];
            a.first.push_back(xg * (wbRad[s] * wbRad[s]));
            a.second.push_back(xg);
          }
          for (const auto &kv : acc) {
            const double den = pwSum(kv.second.second.data(), static_cast<std::ptrdiff_t>(kv.second.second.size()));
            beta[kv.first] =
                (den > 0.) ? std::sqrt(pwSum(kv.second.first.data(), static_cast<std::ptrdiff_t>(kv.second.first.size())) / den)
                           : 0.;
          }
        }

        //------------------------------------------ Sio, Skx, Skj
        std::vector<double> w2all;
        std::vector<char> mask;
        for (const Block &b : iob) {
          if (!b.act[k])
            continue;
          const int ns = static_cast<int>(b.rows.size());
          gather(iod, IO.stride, b.rows, rowsd);
          partition(IO, b.rows, grp, parts);
          wq.assign(ns, b.w[k]);
          // the group shares of the block's fit-Q variance
          double w2tot = 0.;
          if (grp && parts.size() > 1) {
            std::vector<float> blk(static_cast<std::size_t>(ns) * IO.stride);
            for (int i = 0; i < ns; ++i)
              std::memcpy(&blk[static_cast<std::size_t>(i) * IO.stride],
                          IO.v + static_cast<std::size_t>(b.rows[i]) * IO.stride,
                          IO.stride * sizeof(float));
            std::vector<float> qs;
            for (int i = 0; i < sh.qsc.n; ++i)
              if (sh.qsc.idx[i] == b.g) {
                qs.push_back(sh.qsc.v[static_cast<std::size_t>(i) * sh.qsc.stride]);
                qs.push_back(sh.qsc.v[static_cast<std::size_t>(i) * sh.qsc.stride + 1]);
              }
            ioniStepVar(blk.data(), IO.stride, ns, qs.empty() ? nullptr : qs.data(), static_cast<int>(qs.size() / 2),
                        w2all);
            for (double x : w2all)
              w2tot += x;
          }
          // position of each block row in `b.rows` (the partition holds
          // record indices)
          std::map<int, int> pos;
          for (int i = 0; i < ns; ++i)
            pos[b.rows[i]] = i;
          const double bt = beta.count(b.g) ? beta[b.g] : 0.;
          for (const auto &pt : parts) {
            const int np = static_cast<int>(pt.second.size());
            // ionisation: the group's own rows
            std::vector<double> sub;
            gather(iod, IO.stride, pt.second, sub);
            std::fill(re.begin(), re.end(), 0.);
            std::fill(im.begin(), im.end(), 0.);
            ioniRows(tau, kNTau, sub.data(), IO.stride, np, wq.data(), cfg, re.data(), im.data());
            addTo(R.S.ioRe, re.data());
            addTo(R.S.ioIm, im.data());
            GroupExponents *ge = nullptr;
            if (grp) {
              ge = &groupEntry(R.groups, pt.first);
              addTo(ge->S.ioRe, re.data());
              addTo(ge->S.ioIm, im.data());
              if (parts.size() == 1) {
                ge->vqio += b.vq[k];
              } else {
                double w2g = 0.;
                for (int i : pt.second)
                  w2g += w2all[pos[i]];
                ge->vqio += (w2tot > 0.) ? b.vq[k] * (w2g / w2tot) : 0.;
              }
            }
            // the knock-on, on the WHOLE block's rows (its largest xi sets
            // the thin-step node density), evaluated on this group's
            if (I.wantKnockon && (knockMap || knockJoint)) {
              mask.assign(ns, 0);
              for (int i : pt.second)
                mask[pos[i]] = 1;
              rid.resize(ns);
              for (int i = 0; i < ns; ++i)
                rid[i] = i;
              wb.assign(ns, bt);
              ones.assign(ns, 1.);
              if (knockMap) {
                std::fill(re.begin(), re.end(), 0.);
                std::fill(im.begin(), im.end(), 0.);
                knockonRowsMasked(tau, kNTau, rowsd.data(), IO.stride, ns, rid.data(), wq.data(), wb.data(),
                                  ones.data(), ns, KnockonPart::kMap, cfg, mask.data(), re.data(), im.data());
                addTo(R.S.kxRe, re.data());
                addTo(R.S.kxIm, im.data());
                if (ge) {
                  addTo(ge->S.kxRe, re.data());
                  addTo(ge->S.kxIm, im.data());
                }
              }
              if (knockJoint) {
                std::fill(re.begin(), re.end(), 0.);
                std::fill(im.begin(), im.end(), 0.);
                knockonRowsMasked(tau, kNTau, rowsd.data(), IO.stride, ns, rid.data(), wq.data(), wb.data(),
                                  ones.data(), ns, KnockonPart::kJoint, cfg, mask.data(), re.data(), im.data());
                addTo(R.S.kjRe, re.data());
                addTo(R.S.kjIm, im.data());
                if (ge) {
                  addTo(ge->S.kjRe, re.data());
                  addTo(ge->S.kjIm, im.data());
                }
              }
            }
          }
        }

        //------------------------------------------ Srad
        if (haveRad) {
          std::vector<int> rrows;
          for (const Block &b : iob) {
            if (!b.act[k])
              continue;
            rrows.clear();
            for (int s = 0; s < RD.n; ++s)
              if (RD.idx[s] == b.g)
                rrows.push_back(s);
            if (rrows.empty())
              continue;
            partition(RD, rrows, grp, parts);
            for (const auto &pt : parts) {
              const int np = static_cast<int>(pt.second.size());
              gather(radd, RD.stride, pt.second, rowsd);
              specsub.resize(static_cast<std::size_t>(np) * 2 * nvf);
              for (int i = 0; i < np; ++i)
                std::memcpy(&specsub[static_cast<std::size_t>(i) * 2 * nvf],
                            &specf[static_cast<std::size_t>(pt.second[i]) * 2 * nvf],
                            2 * nvf * sizeof(double));
              rid.resize(np);
              wb.resize(np);
              for (int i = 0; i < np; ++i) {
                rid[i] = i;
                wb[i] = wbRad[pt.second[i]];
              }
              wq.assign(np, b.w[k]);
              ones.assign(np, 1.);
              std::fill(re.begin(), re.end(), 0.);
              std::fill(im.begin(), im.end(), 0.);
              radRows(tau, kNTau, rowsd.data(), RD.stride, np, specsub.data(), vf.data(), nvf, rid.data(), wq.data(),
                      wb.data(), ones.data(), np, cfg.qopExact, re.data(), im.data());
              addTo(R.S.radRe, re.data());
              addTo(R.S.radIm, im.data());
              if (grp) {
                GroupExponents &ge = groupEntry(R.groups, pt.first);
                addTo(ge.S.radRe, re.data());
                addTo(ge.S.radIm, im.data());
              }
            }
          }
        }

        //------------------------------------------ nuclear elastic
        if (anyNucel && I.wantNucel && I.nucel != nullptr)
          nucelFunctional(sh, I, wms, wio, cfg, R);
      }

      // A functional that asked for the split carries its flat families AS
      // the sum of the group parts, ascending in group, so the closure
      // `sum_g S_g == S` is exact.
      for (int k = 0; k < nf; ++k) {
        if (!alive[k] || !in[k].wantGroups)
          continue;
        Exponents flat;
        flat.clear();
        for (const GroupExponents &g : out[k].groups) {
          addTo(flat.ms, g.S.ms.data());
          addTo(flat.ioRe, g.S.ioRe.data());
          addTo(flat.ioIm, g.S.ioIm.data());
          addTo(flat.radRe, g.S.radRe.data());
          addTo(flat.radIm, g.S.radIm.data());
          addTo(flat.kxRe, g.S.kxRe.data());
          addTo(flat.kxIm, g.S.kxIm.data());
          addTo(flat.kjRe, g.S.kjRe.data());
          addTo(flat.kjIm, g.S.kjIm.data());
        }
        out[k].S = flat;
      }
    }

  }  // namespace

  void trackExponents(const TrackInput &in, TrackResult &out) { trackExponentsImpl(&in, 1, &out); }

  void trackExponents(const TrackInput *in, int nfunc, TrackResult *out) { trackExponentsImpl(in, nfunc, out); }

  //==========================================================================
  // 8. PROVENANCE
  //==========================================================================
  std::string modelTag(const RowConfig &c, bool withNucel) {
    char b[512];
    std::snprintf(b,
                  sizeof(b),
                  "cvhcf/3 rows=cf_rows tau=stride4of448<=8 tab=%s "
                  "ms:elecTmax=1,elecEdge=1,fineG=1,snapYmax=0,wviSplit=0 "
                  "ioni:kokoulin=%g,kokTcut=%g,kokQuad=filon,species=regimeSnap,a3=%g,exc=%g,tmaxScale=%g "
                  "knockon:joint=%d,qopExact=%d,nperdec=%d,nperdecThin=%d,thin=%g,tcut=%gMeV,pminFrac=%g "
                  "rad:nsub=%d,recoil=1,mass=species",
                  kShapeTableId.c_str(),
                  c.ioniKokoulin,
                  c.ioniKokoulinTcut,
                  c.ioniA3Scale,
                  c.ioniExcScale,
                  c.ioniTmaxScale,
                  int(c.knockonJoint),
                  int(c.qopExact),
                  c.knockonNPerDec,
                  c.knockonNPerDecThin,
                  c.knockonThin,
                  c.knockonTcut,
                  c.pminFrac,
                  c.radNsub);
    std::string t(b);
    if (withNucel)
      t += " nuc=" + nucelTableId() + ":recoil=" + std::to_string(int(c.nucelRecoil)) +
           ",joint=" + std::to_string(int(c.nucelJoint));
    return t;
  }

}  // namespace cvhcf
