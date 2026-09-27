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
#include <mutex>
#include <utility>
#include <vector>

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
    del.fill(0.);
    ioRe.fill(0.);
    ioIm.fill(0.);
    radRe.fill(0.);
    radIm.fill(0.);
  }

  //==========================================================================
  // 1. THE SHAPE TABLES  (cf_ms_exact._GE, cf_delta_ray._PHI)
  //==========================================================================
  namespace {

    struct ShapeTables {
      int nTau = 0, nY = 0, nPhi = 0;
      std::vector<double> gtau, elecY, hard, dipole, philx, phi;
      // cached logs of the Y axis, for the row interpolation
      double lnY0 = 0., dlnY = 0.;
      double lnPhi0 = 0., dlnPhi = 0.;
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
        path = edm::FileInPath("TrackPropagation/Geant4e/data/cvhcf_gshape_elec_v1.bin").fullPath();
      std::ifstream f(path, std::ios::binary);
      if (!f)
        throw cms::Exception("CvhCfExponents") << "cannot open the Moliere shape table " << path;
      char magic[9] = {0};
      f.read(magic, 8);
      if (std::strncmp(magic, "CVHCFGSH", 8) != 0)
        throw cms::Exception("CvhCfExponents") << path << " is not a cvhcf shape table";
      int hdr[4] = {0, 0, 0, 0};
      f.read(reinterpret_cast<char *>(hdr), sizeof(hdr));
      if (hdr[0] != 1)
        throw cms::Exception("CvhCfExponents") << path << ": unsupported table version " << hdr[0];
      T.nTau = hdr[1];
      T.nY = hdr[2];
      T.nPhi = hdr[3];
      if (T.nTau < 2 || T.nY < 2 || T.nPhi < 2)
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
      rd(T.philx, T.nPhi);
      rd(T.phi, T.nPhi);
      T.lnY0 = std::log(T.elecY[0]);
      T.dlnY = std::log(T.elecY[1]) - std::log(T.elecY[0]);
      T.lnPhi0 = T.philx[0];
      T.dlnPhi = T.philx[1] - T.philx[0];
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

    // --- cf_delta_ray.phi_tab --------------------------------------------
    inline double phiTab(const ShapeTables &t, double x) {
      const double lx = std::log(std::max(x, 1e-300));
      if (lx < t.philx[0])
        return t.phi[0] + 0.5 * (lx - t.philx[0]);
      if (lx > t.philx[t.nPhi - 1])
        return -std::exp(-2. * lx);
      // uniform grid in ln x
      double fi = (lx - t.lnPhi0) / t.dlnPhi;
      int j = static_cast<int>(fi);
      j = std::max(0, std::min(j, t.nPhi - 2));
      const double slope = (t.phi[j + 1] - t.phi[j]) / (t.philx[j + 1] - t.philx[j]);
      return slope * (lx - t.philx[j]) + t.phi[j];
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
  // 3. THE MULTIPLE-SCATTERING AND DELTA-RAY FAMILIES
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

    // The reference's `ok` (thp2 > 0) and `act` (chi_c^2, chi_a^2 > 0) masks,
    // applied while the Moliere parameters are formed.
    void buildMsSteps(const float *rows, int stride, int n, std::vector<MsStep> &out) {
      out.clear();
      out.reserve(static_cast<std::size_t>(n));
      const bool haveSums = stride >= 10;
      for (int i = 0; i < n; ++i) {
        const float *r = rows + static_cast<std::size_t>(i) * stride;
        if (!(r[5] > 0.f))
          continue;
        MsStep s;
        s.effZ = r[0];
        s.effA = r[1];
        s.xg = r[2];
        s.pGeV = r[3];
        s.beta = r[4];
        s.thp2 = r[5];
        const MolPars mp =
            moliereParams(s.effZ, s.effA, s.xg, s.pGeV, s.beta, haveSums ? r[7] : 0., haveSums ? r[8] : 0., haveSums);
        if (!(mp.chic2 > 0.) || !(mp.chia2 > 0.))
          continue;
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
        const double tmx = 2. * kMeGeV * bg * bg / (1. + 2. * gam * rat + rat * rat);
        const double te = 2. * tmx * kMeGeV / (s.pGeV * s.pGeV);
        const double yme = std::min(std::sqrt(te / s.chia2), ym);
        s.lgEle = round3(std::log(clipd(yme, 1e1, 1e7)));
        s.fN = s.effZ / (s.effZ + 1.);
        s.fE = 1. / (s.effZ + 1.);
        out.push_back(s);
      }
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

    //------------------------------------------------------------------------
    // THE DISCRETE DELTA-RAY RECOIL  (cf_delta_ray)
    //
    // Above the e- production cut Geant4 makes a REAL knock-on electron and
    // conserves 4-momentum, so the muon is deflected by p_eT = sqrt(2 m_e T).
    // Moliere's Z(Z+1) already carries that scattering CONTINUOUSLY, so the
    // discrete block REPLACES a variance-matched share of it (`carve_factor`)
    // rather than adding to it -- which is why this needs the block's own MS
    // exponent as an input.
    //------------------------------------------------------------------------
    constexpr double kDelME = 0.51099895e-3;  // GeV
    constexpr double kDelHalfK = 0.1535;      // MeV cm^2/g
    constexpr double kDelTcut = 0.35e-3;      // GeV  (CF_DELTA_TCUT)
    constexpr double kDelTmaxCap = 0.05;      // GeV  (CF_DELTA_TMAXCAP)

    // cf_delta_ray._tmx: note the beta clip here is 1 - 1e-12, NOT the
    // 1 - 1e-15 the xi and spin factors use. The reference has both and they
    // are not interchangeable at beta -> 1.
    inline double delTmx(double pGeV, double beta) {
      const double bt = clipd(beta, 1e-9, 1. - 1e-12);
      const double g = 1. / std::sqrt(1. - bt * bt);
      const double m = std::max(pGeV / (bt * g), 1e-6);
      const double bg = bt * g;
      const double r = kDelME / std::max(m, 1e-6);
      const double t = 2. * kDelME * bg * bg / (1. + 2. * g * r + r * r);
      return std::min(t, kDelTmaxCap);
    }

    // The two halves of the delta family, so that a per-group split can reuse
    // the block's carve. `acc` may be null, in which case only vd/vms are
    // formed (the tau loop is the expensive part).
    void delAccumImpl(const float *rows, int stride, int n, const double *wstd, int nk, const double *tau, int nt,
                      double *acc, double &vd, double &vms) {
      const ShapeTables &tb = T();
      vd = 0.;
      vms = 0.;
      for (int i = 0; i < n; ++i)
        vms += rows[static_cast<std::size_t>(i) * stride + 5];
      for (int i = 0; i < n; ++i) {
        const float *r = rows + static_cast<std::size_t>(i) * stride;
        const double effZ = r[0], effA = r[1], xg = r[2], p = r[3];
        const double bt = clipd(r[4], 1e-9, 1. - 1e-15);
        const double tmx = delTmx(p, r[4]);
        if (!(xg > 0.) || !(tmx > kDelTcut) || !(p > 0.))
          continue;
        const double xi = kDelHalfK * (effZ / std::max(effA, 1.)) * xg / std::max(bt * bt, 1e-9) * 1e-3;
        vd += xi * kDelME * std::log(tmx / kDelTcut) / (p * p);
        if (acc == nullptr)
          continue;
        const double sc = std::sqrt(kDelTcut), sh = std::sqrt(tmx);
        // the (1 - beta^2 T/Tmax) spin-0 term of the PDG delta spectrum, as a
        // first-order variance correction
        const double spin = 1. - 0.5 * bt * bt / std::log(std::max(tmx / kDelTcut, 1.0001));
        for (int k = 0; k < nk; ++k) {
          const double a0 = wstd[k] * std::sqrt(2. * kDelME) / p;
          double *ak = acc + static_cast<std::size_t>(k) * nt;
          for (int j = 0; j < nt; ++j) {
            const double a = a0 * tau[j];
            ak[j] += xi * a * a * (phiTab(tb, a * sc) - phiTab(tb, a * sh)) * spin;
          }
        }
      }
    }

    // `Sms` and `S` are `nk * nt`. The carve is a property of the BLOCK and
    // not of the weight (`vd` and `vms` do not depend on it), so it is formed
    // once and applied to every functional.
    void delExponentImpl(const float *rows, int stride, int n, const double *wstd, int nk, const double *tau, int nt,
                         const double *Sms, double *S) {
      if (n <= 0 || nk <= 0)
        return;
      const std::size_t m = static_cast<std::size_t>(nk) * nt;
      std::vector<double> acc(m, 0.);
      double vd = 0., vms = 0.;
      delAccumImpl(rows, stride, n, wstd, nk, tau, nt, acc.data(), vd, vms);
      const double carve = (vms > 0.) ? clipd(vd / vms, 0., 0.5) : 0.;
      for (std::size_t j = 0; j < m; ++j)
        S[j] += acc[j] - carve * Sms[j];
    }

  }  // namespace

  namespace {

    // THE BLOCK PRIMITIVES ON THE CONCATENATED ARGUMENT LIST.  `wstd` holds
    // `nk` weights; every output slab is `nk * kNTau`, functional-major. The
    // public single-weight primitives below are these at `nk == 1` -- not
    // merely equal to them, literally the same call.
    void msBlockMulti(const float *rows, int stride, int n, const double *wstd, int nk, double *Sms, double *Sdel) {
      if (rows == nullptr || n <= 0 || stride < 8 || nk <= 0)
        return;
      // `buildMsSteps` -- the Moliere parameters, the screening angles, the
      // electron ceiling and their roundings -- is the weight-independent
      // half of the block, and is now paid once for all `nk`.
      std::vector<MsStep> st;
      buildMsSteps(rows, stride, n, st);
      const std::size_t m = static_cast<std::size_t>(nk) * kNTau;
      std::vector<double> local(m, 0.);
      msExponentImpl(st, wstd, nk, tauGrid(), kNTau, local.data());
      if (Sms != nullptr)
        for (std::size_t j = 0; j < m; ++j)
          Sms[j] += local[j];
      if (Sdel != nullptr)
        delExponentImpl(rows, stride, n, wstd, nk, tauGrid(), kNTau, local.data(), Sdel);
    }

    void delBlockCarvedMulti(const float *rows, int stride, int n, const double *wstd, int nk, const double *Sms,
                             double carve, double *Sdel) {
      if (rows == nullptr || n <= 0 || stride < 8 || nk <= 0 || Sdel == nullptr || Sms == nullptr)
        return;
      const std::size_t m = static_cast<std::size_t>(nk) * kNTau;
      std::vector<double> acc(m, 0.);
      double vd = 0., vms = 0.;
      delAccumImpl(rows, stride, n, wstd, nk, tauGrid(), kNTau, acc.data(), vd, vms);
      for (std::size_t j = 0; j < m; ++j)
        Sdel[j] += acc[j] - carve * Sms[j];
    }

  }  // namespace

  double delCarveFactor(const float *rows, int stride, int n) {
    if (rows == nullptr || n <= 0 || stride < 8)
      return 0.;
    double vd = 0., vms = 0.;
    const double one = 1.;
    delAccumImpl(rows, stride, n, &one, 1, tauGrid(), kNTau, nullptr, vd, vms);
    return (vms > 0.) ? clipd(vd / vms, 0., 0.5) : 0.;
  }

  void delBlockCarved(
      const float *rows, int stride, int n, double wstd, const double *Sms, double carve, double *Sdel) {
    delBlockCarvedMulti(rows, stride, n, &wstd, 1, Sms, carve, Sdel);
  }

  void msBlock(const float *rows, int stride, int n, double wstd, double *Sms, double *Sdel) {
    msBlockMulti(rows, stride, n, &wstd, 1, Sms, Sdel);
  }

  //==========================================================================
  // 4. THE IONIZATION AND RADIATIVE FAMILIES
  //==========================================================================
  namespace {

    // One pooled ionization block, as cvhcgf sees it. The energies carry the
    // record's `scaling` factor -- and `a3` does too, but ONLY in regimes 2/3
    // where the slot holds an energy (xi) rather than a collision count. That
    // is the branch `cf_track_resolution.ioni_step_exponent` takes, and the
    // same one `Geant4ePropagator::toIoniStep` takes on the fit side.
    void buildIoniSteps(const float *rows, int stride, int n, double gsScale, std::vector<cvhcgf::IoniStep> &out) {
      out.clear();
      out.reserve(static_cast<std::size_t>(n));
      for (int i = 0; i < n; ++i) {
        const float *r = rows + static_cast<std::size_t>(i) * stride;
        cvhcgf::IoniStep s;
        s.regime = static_cast<int>(r[0]);
        s.gsig2 = r[1];
        const double gam = r[9];
        s.a1 = r[2];
        s.e1 = r[3] * gam;
        s.a2 = r[4];
        s.e2 = r[5] * gam;
        s.a3 = (s.regime >= 2) ? (static_cast<double>(r[6]) * gam) : static_cast<double>(r[6]);
        s.e0 = r[7] * gam;
        s.tmax = r[8] * gam;
        s.gs = gsScale * (static_cast<double>(r[10]) * 1e-3);
        if (stride >= 13) {
          s.beta2 = r[11];
          s.etot = r[12];
        } else if (s.regime >= 2) {
          throw cms::Exception("CvhCfExponents")
              << "ioniurbanv stride " << stride << " carries a regime-" << s.regime
              << " record but no beta^2/E columns; they cannot be recovered without the particle mass";
        }
        // Kokoulin is NOT part of the offline reference model (default OFF on
        // both sides -- cvhcgf::ioniKokoulinEnabled and
        // cf_track_resolution.IONI_KOKOULIN read the same switch), and it is
        // deliberately not wired in here: the exported exponent must be the
        // model the closure is quoted against.
        s.kokNbin = 0;
        out.push_back(s);
      }
    }

    struct RadRec {
      std::vector<double> dNdv;
      double etotGeV = 0.;
      double gs = 0.;  // (residual units) per GeV: cs * w, in the reference's order
    };

  }  // namespace

  namespace {

    // The ionization block on the concatenated argument list. `blockExponent`
    // reads the weight only as `gs * t`, so the steps are PARSED ONCE at unit
    // weight -- `gs` then holds the step's own dE -> residual map -- and each
    // functional's weight is folded in by the same multiplication
    // `buildIoniSteps` would have done (`1.0 * x == x` exactly, so the cached
    // factor is bit-for-bit the record's own).
    void ioniBlockMulti(
        const float *rows, int stride, int n, const double *wstdSigned, int nk, double *Sre, double *Sim) {
      if (rows == nullptr || n <= 0 || stride < 11 || stride > 14 || nk <= 0)
        return;
      cvhcgf::Block blk;
      buildIoniSteps(rows, stride, n, 1.0, blk.ioni);
      const std::size_t nst = blk.ioni.size();
      std::vector<double> gsraw(nst);
      for (std::size_t i = 0; i < nst; ++i)
        gsraw[i] = blk.ioni[i].gs;
      const double *tau = tauGrid();
      for (int k = 0; k < nk; ++k) {
        for (std::size_t i = 0; i < nst; ++i)
          blk.ioni[i].gs = wstdSigned[k] * gsraw[i];
        double *Rk = Sre + static_cast<std::size_t>(k) * kNTau;
        double *Ik = Sim + static_cast<std::size_t>(k) * kNTau;
        for (int j = 0; j < kNTau; ++j) {
          const std::complex<double> s = cvhcgf::blockExponent(blk, tau[j]);
          Rk[j] += s.real();
          Ik[j] += s.imag();
        }
      }
    }

    // The radiative channel of the same block. `makeRadSpectrum` -- the
    // expensive half, and the reason the radiative channel is two thirds of
    // the CF cost -- knows nothing about the functional, so the spectrum of a
    // step is built ONCE and every `k` reads it.
    void radBlockMulti(const float *rows, int stride, int n, const float *spec, const float *vgridf, int nv,
                       const double *wstdSigned, int nk, double *Sre, double *Sim) {
      if (rows == nullptr || spec == nullptr || vgridf == nullptr || n <= 0 || nv < 2 || nk <= 0)
        return;
      std::vector<double> vgrid(nv), shapeB(nv), shapeP(nv), dNdv(nv), wtrap(nv, 0.);
      for (int i = 0; i < nv; ++i)
        vgrid[i] = vgridf[i];
      for (int i = 0; i + 1 < nv; ++i) {
        const double dv = 0.5 * (vgrid[i + 1] - vgrid[i]);
        wtrap[i] += dv;
        wtrap[i + 1] += dv;
      }
      const double *tau = tauGrid();
      std::vector<double> a(kNTau), ve(nv);
      for (int ir = 0; ir < n; ++ir) {
        const float *r = rows + static_cast<std::size_t>(ir) * stride;
        const float *sp = spec + static_cast<std::size_t>(ir) * 2 * nv;
        for (int i = 0; i < nv; ++i) {
          shapeB[i] = sp[i];
          shapeP[i] = sp[nv + i];
        }
        const double etot = r[3];    // GeV
        const double stepCm = r[6];  // cm
        if (!(etot > 0.))
          continue;
        // GeV throughout, matching cf_brems_exact: dE/norm is unit-free, so
        // the reference's own units are used rather than the propagator's MeV
        // convention, and `gs` below carries no 1e-3 -- the radiative records
        // store GeV where the Urban ones store MeV.
        cvhcgf::makeRadSpectrum(vgrid.data(), shapeB.data(), shapeP.data(), r[8] * stepCm, r[9] * stepCm, etot, nv,
                                dNdv.data());
        bool any = false;
        for (int i = 0; i < nv; ++i)
          if (dNdv[i] > 0.) {
            any = true;
            break;
          }
        if (!any)
          continue;
        for (int i = 0; i < nv; ++i)
          ve[i] = vgrid[i] * etot;
        for (int k = 0; k < nk; ++k) {
          const double gs = static_cast<double>(r[10]) * wstdSigned[k];
          if (gs == 0.)
            continue;
          // The reference forms `a = tau * cs * w` and `x = outer(a, v * E)`.
          // Keeping that association is what makes the two agree to the last
          // bits rather than merely to 1e-16.
          for (int j = 0; j < kNTau; ++j)
            a[j] = tau[j] * gs;
          double *Rk = Sre + static_cast<std::size_t>(k) * kNTau;
          double *Ik = Sim + static_cast<std::size_t>(k) * kNTau;
          for (int i = 0; i < nv; ++i) {
            const double q = wtrap[i] * dNdv[i];
            if (q == 0.)
              continue;
            for (int j = 0; j < kNTau; ++j) {
              const double x = a[j] * ve[i];
              if (std::fabs(x) < 1e-4) {
                Rk[j] += q * (-0.5 * x * x);
                Ik[j] += q * (x * x * x / 6.);
              } else {
                // cos x - 1 as -2 sin^2(x/2), the cancellation-free form; the
                // same two sines cf_brems_exact and cvhcgf::blockExponent
                // spend.
                const double s2 = std::sin(0.5 * x);
                Rk[j] += q * (-2.0 * s2 * s2);
                Ik[j] += q * (std::sin(x) - x);
              }
            }
          }
        }
      }
    }

  }  // namespace

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

  void ioniBlock(const float *rows, int stride, int n, double wstdSigned, double *Sre, double *Sim) {
    // Accepted strides: 11 = the base record, 13 = with the exact-delta
    // beta^2/E columns, 12 / 14 = the same two with the material-group column
    // appended. Anything else is not an `ioniurbanv` row and is refused rather
    // than mis-parsed -- but note the guard SILENTLY returns, so a stride that
    // is not on this list zeroes the ionization exponent of every candidate;
    // keep the list in step with the writer.
    ioniBlockMulti(rows, stride, n, &wstdSigned, 1, Sre, Sim);
  }

  void radBlock(const float *rows, int stride, int n, const float *spec, const float *vgridf, int nv,
                double wstdSigned, double *Sre, double *Sim) {
    radBlockMulti(rows, stride, n, spec, vgridf, nv, &wstdSigned, 1, Sre, Sim);
  }

  //==========================================================================
  // 4b. THE NUCLEAR-ELASTIC FAMILY  (ksclosure/nucel/nucel_tables.Table and
  //     ks_nucel_cf.step_family; layout and rules in
  //     data/make_cvhcf_nucel_tables.py)
  //==========================================================================
  namespace {

    const std::string kNucelTableId = "cvhcf_nucel_v1";

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
      if (hdr[0] != 1)
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
  };

  void NucelExponents::clear() {
    ang.fill(0.);
    recRe.fill(0.);
    recIm.fill(0.);
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

    // g - 1 of one node record at u >= 0 (`Table.g_node` minus one).
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

    // ADDS c * (h(v_j) - 1) of node n over the tau grid, v_j = wq tau_j
    // (`Table.h_node`).  The tau grid is uniform (tau_j = j dt), so each bin's
    // phase exp(i j theta_b) and the sin(j phi_b) of its sinc are advanced by
    // complex recurrence -- two sincos per bin instead of two per bin and tau
    // point -- in four independent lanes (j, j+1, j+2, j+3 stepped by
    // exp(4 i theta)), so that the tau points of one bin carry no dependency
    // on each other.  The recurrence drifts by ~j ulp (<= 1e-14 here), far
    // below the float32 floor of the export.
    constexpr int kLanes = 4;
    constexpr int kNBlk = (kNTau - 1 + kLanes - 1) / kLanes;  // blocks covering j = 1 .. kNTau-1

    void nucelRecoilNode(const NucelMix &m, int n, double wq, double c, double *re, double *im) {
      const double *tau = tauGrid();
      const double *mo = &m.mom[static_cast<std::size_t>(n) * 7];
      // the recoil below the first bin, by its partial moments
      for (int j = 0; j < kNTau; ++j) {
        const double v = wq * tau[j];
        re[j] += c * (mo[4] - 1. - 0.5 * v * v * mo[6]);
        im[j] += c * (v * mo[5]);
      }
      const int b0 = m.binStart[n], b1 = m.binStart[n + 1];
      if (b1 <= b0)
        return;
      const double dt = tau[1];  // tau_j = j dt
      // accumulators for j = 1 .. kLanes * kNBlk (the tail beyond kNTau-1 is dropped)
      double ar[kLanes * kNBlk] = {0.}, ai[kLanes * kNBlk] = {0.};
      double ij[kLanes * kNBlk];
      for (int k = 0; k < kLanes * kNBlk; ++k)
        ij[k] = 1. / (k + 1);
      double p0 = 0.;
      for (int b = b0; b < b1; ++b) {
        const double Pb = c * m.bP[b];
        p0 += Pb;
        const double th = wq * dt * m.bMu[b];
        const double ph = wq * dt * m.bHw[b];
        const double c1 = std::cos(th), s1 = std::sin(th);
        const double d1 = std::cos(ph), e1 = std::sin(ph);
        // sinc(j ph) = Im exp(i j ph) / (j ph); a zero-width bin has sinc 1
        const double iph = (ph != 0.) ? 1. / ph : 0.;
        const double sz = (ph != 0.) ? 0. : 1.;
        // lanes at j = 1..4, and the step exp(4 i theta), exp(4 i phi)
        double zr[kLanes], zi[kLanes], yr[kLanes], yi[kLanes];
        zr[0] = c1;
        zi[0] = s1;
        yr[0] = d1;
        yi[0] = e1;
        for (int l = 1; l < kLanes; ++l) {
          zr[l] = zr[l - 1] * c1 - zi[l - 1] * s1;
          zi[l] = zr[l - 1] * s1 + zi[l - 1] * c1;
          yr[l] = yr[l - 1] * d1 - yi[l - 1] * e1;
          yi[l] = yr[l - 1] * e1 + yi[l - 1] * d1;
        }
        const double c4 = zr[kLanes - 1], s4 = zi[kLanes - 1];
        const double d4 = yr[kLanes - 1], e4 = yi[kLanes - 1];
        for (int k = 0; k < kNBlk; ++k) {
          double *__restrict pr = ar + k * kLanes;
          double *__restrict pi = ai + k * kLanes;
          const double *__restrict pj = ij + k * kLanes;
          for (int l = 0; l < kLanes; ++l) {
            const double q = Pb * (yi[l] * iph * pj[l] + sz);
            pr[l] += q * zr[l];
            pi[l] += q * zi[l];
          }
          for (int l = 0; l < kLanes; ++l) {
            const double zr2 = zr[l] * c4 - zi[l] * s4;
            zi[l] = zr[l] * s4 + zi[l] * c4;
            zr[l] = zr2;
            const double yr2 = yr[l] * d4 - yi[l] * e4;
            yi[l] = yr[l] * e4 + yi[l] * d4;
            yr[l] = yr2;
          }
        }
      }
      re[0] += p0;
      for (int j = 1; j < kNTau; ++j) {
        re[j] += ar[j - 1];
        im[j] += ai[j - 1];
      }
    }

    // One row's contribution: ADDS N_s (g - 1) and N_s (h - 1) into the
    // caller's arrays and returns N_s (0 when the row has no species).
    double nucelRow(const NucelTable &t,
                    const float *r,
                    int matIndex,
                    int pdgCode,
                    double wAng,
                    double wRec,
                    NucelMixtures &mix,
                    double *ang,
                    double *re,
                    double *im) {
      if (pdgCode == 0)
        return 0.;
      const double xg = r[2];
      const double p = r[3];
      if (!(xg > 0.) || !(p > 0.))
        return 0.;
      const NucelMix &m = mix.get(pdgCode, matIndex);
      const NucelTable::Species &S = t.sp[m.isp];
      const double Ns = std::exp(logLog(t, m.lmu.data(), std::log(p))) * xg;
      if (!(Ns > 0.))
        return 0.;
      int j0, j1;
      double f;
      nucelNodes(S, p, j0, j1, f);
      const int nodes[2] = {j0, j1};
      const double wts[2] = {1. - f, f};
      const double *tau = tauGrid();
      if (wAng != 0.) {
        for (int q = 0; q < 2; ++q) {
          if (wts[q] == 0.)
            continue;
          const int n = nodes[q];
          const double sc = S.nodeP[n] / p;  // fixed momentum transfer
          const double *gm = &m.gm1[static_cast<std::size_t>(n) * t.nU];
          const double th2 = m.mom[static_cast<std::size_t>(n) * 7 + 1];
          const double c = Ns * wts[q];
          for (int j = 0; j < kNTau; ++j)
            ang[j] += c * nucelGm1(t, gm, th2, (std::fabs(wAng) * tau[j]) * sc);
        }
      }
      if (wRec != 0.) {
        for (int q = 0; q < 2; ++q) {
          if (wts[q] == 0.)
            continue;
          nucelRecoilNode(m, nodes[q], wRec, Ns * wts[q], re, im);
        }
      }
      return Ns;
    }

  }  // namespace

  void nucelSteps(const float *rows,
                  int stride,
                  int n,
                  const int *mat,
                  const int *pdgv,
                  const double *wAng,
                  const double *wRec,
                  NucelMixtures &mix,
                  NucelExponents &out) {
    if (rows == nullptr || mat == nullptr || pdgv == nullptr || n <= 0 || stride < 4)
      return;
    const NucelTable &t = NT();
    for (int i = 0; i < n; ++i)
      out.N += nucelRow(t,
                        rows + static_cast<std::size_t>(i) * stride,
                        mat[i],
                        pdgv[i],
                        wAng ? wAng[i] : 0.,
                        wRec ? wRec[i] : 0.,
                        mix,
                        out.ang.data(),
                        out.recRe.data(),
                        out.recIm.data());
  }

  namespace {

    // The family of one functional from the MS rows and the block weights
    // the pooling pass formed: `wms` / `wio` are (global index, weight) of
    // the MS and (signed) ionisation blocks, ascending in index.
    void nucelFunctional(const TrackInput &sh,
                         const TrackInput &in,
                         const std::vector<std::pair<unsigned int, double>> &wms,
                         const std::vector<std::pair<unsigned int, double>> &wio,
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
      const NucelTable &t = NT();
      auto look = [](const std::vector<std::pair<unsigned int, double>> &v, unsigned int g) {
        auto it = std::lower_bound(
            v.begin(), v.end(), g, [](const std::pair<unsigned int, double> &a, unsigned int b) { return a.first < b; });
        return (it != v.end() && it->first == g) ? it->second : 0.;
      };
      const bool groups = in.wantGroups && ms.groupCol >= 0 && ms.groupCol < ms.stride;
      std::array<double, kNTau> a{}, re{}, im{};
      for (int i = 0; i < ms.n; ++i) {
        if (sh.mspdg[i] == 0)
          continue;
        const float *r = ms.v + static_cast<std::size_t>(i) * ms.stride;
        const float *rr = sh.rad.v + static_cast<std::size_t>(i) * sh.rad.stride;
        // the pairing is exact by construction; check it on the momentum
        // both records carry (the same double, stored as float twice)
        if (rr[4] != r[3])
          throw cms::Exception("CvhCfExponents")
              << "nuclear-elastic family: radiative row " << i << " (p = " << rr[4]
              << ") is not the MS row's step (p = " << r[3] << ")";
        const double wA = look(wms, ms.idx[i]);
        const double wR = look(wio, sh.rad.idx[i]) * static_cast<double>(rr[10]) * 1e-3;
        if (groups) {
          a.fill(0.);
          re.fill(0.);
          im.fill(0.);
          const double Ns = nucelRow(t, r, sh.msmat[i], sh.mspdg[i], wA, wR, *in.nucel, a.data(), re.data(), im.data());
          if (Ns == 0.)
            continue;
          const int g = static_cast<int>(r[ms.groupCol]);
          auto it = std::lower_bound(out.nucGroups.begin(),
                                     out.nucGroups.end(),
                                     g,
                                     [](const std::pair<int, NucelExponents> &x, int y) { return x.first < y; });
          if (it == out.nucGroups.end() || it->first != g)
            it = out.nucGroups.insert(it, std::make_pair(g, NucelExponents()));
          NucelExponents &ge = it->second;
          ge.N += Ns;
          out.nuc.N += Ns;
          for (int j = 0; j < kNTau; ++j) {
            ge.ang[j] += a[j];
            ge.recRe[j] += re[j];
            ge.recIm[j] += im[j];
            out.nuc.ang[j] += a[j];
            out.nuc.recRe[j] += re[j];
            out.nuc.recIm[j] += im[j];
          }
        } else {
          out.nuc.N += nucelRow(t,
                                r,
                                sh.msmat[i],
                                sh.mspdg[i],
                                wA,
                                wR,
                                *in.nucel,
                                out.nuc.ang.data(),
                                out.nuc.recRe.data(),
                                out.nuc.recIm.data());
        }
      }
    }

  }  // namespace

  //==========================================================================
  // 5. THE POOLING
  //
  // Identical to `cf_track_resolution.extract`: for each material family in
  // (10, 11), for each DISTINCT global parameter index among the resolution
  // entries of that family, pool v_b over the entries, gather the step rows
  // whose index VALUE matches, and form the block's scalar standardized
  // weight. Ascending index order is `np.unique`'s, and only the float
  // summation order depends on it -- which is exactly what a bit-level
  // comparison sees.
  //==========================================================================
  namespace {

    // Accumulate `src` into the (group -> exponents) list, creating the entry
    // on first use. The list is kept ASCENDING in `group` so that the export
    // order is a property of the candidate and not of the propagation.
    GroupExponents &groupEntry(std::vector<GroupExponents> &v, int g) {
      auto it = std::lower_bound(
          v.begin(), v.end(), g, [](const GroupExponents &a, int b) { return a.group < b; });
      if (it != v.end() && it->group == g)
        return *it;
      GroupExponents e;
      e.group = g;
      return *v.insert(it, e);
    }

    // Per-STEP ionization variance contribution, with the per-leg CGF scale
    // applied, so that a material group's share of `ioniSq2` is a plain sum
    // over its own steps.  Same association, left to right, as `ioniSq2`.
    void ioniStepVar(
        const float *rows, int stride, int n, const float *qsc, int nqsc, std::vector<double> &w2) {
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

    Exponents &groupSlot(std::vector<GroupExponents> &v, int g) {
      auto it = std::lower_bound(
          v.begin(), v.end(), g, [](const GroupExponents &a, int b) { return a.group < b; });
      if (it != v.end() && it->group == g)
        return it->S;
      GroupExponents e;
      e.group = g;
      return v.insert(it, e)->S;
    }

    void addInto(Exponents &dst, const Exponents &src) {
      for (int j = 0; j < kNTau; ++j) {
        dst.ms[j] += src.ms[j];
        dst.del[j] += src.del[j];
        dst.ioRe[j] += src.ioRe[j];
        dst.ioIm[j] += src.ioIm[j];
        dst.radRe[j] += src.radRe[j];
        dst.radIm[j] += src.radIm[j];
      }
    }

    // THE SIX-FAMILY SCRATCH SLAB of the multi-functional pass: family-major,
    // then functional, then tau, i.e. `buf[(f * nact + a) * kNTau + j]`, so
    // each family's `nact * kNTau` block is contiguous and can be handed to a
    // `*Multi` primitive directly.
    enum SlabFam { kSlabMs = 0, kSlabDel, kSlabIoRe, kSlabIoIm, kSlabRadRe, kSlabRadIm, kNSlabFam };

    inline double *slab(std::vector<double> &buf, int f, int nact) {
      return buf.data() + static_cast<std::size_t>(f) * nact * kNTau;
    }

    // `addInto` for functional `a` of such a slab. `useDel` false adds a
    // literal zero to the delta family, which is what the single-functional
    // path does when it hands `addInto` a cleared `Exponents` -- the families
    // a branch did not fill are added, not skipped.
    void addIntoSlice(Exponents &dst, const std::vector<double> &buf, int a, int nact, bool useDel) {
      const std::size_t o = static_cast<std::size_t>(a) * kNTau;
      const double *ms = buf.data() + static_cast<std::size_t>(kSlabMs) * nact * kNTau + o;
      const double *del = buf.data() + static_cast<std::size_t>(kSlabDel) * nact * kNTau + o;
      const double *iore = buf.data() + static_cast<std::size_t>(kSlabIoRe) * nact * kNTau + o;
      const double *ioim = buf.data() + static_cast<std::size_t>(kSlabIoIm) * nact * kNTau + o;
      const double *radre = buf.data() + static_cast<std::size_t>(kSlabRadRe) * nact * kNTau + o;
      const double *radim = buf.data() + static_cast<std::size_t>(kSlabRadIm) * nact * kNTau + o;
      for (int j = 0; j < kNTau; ++j) {
        dst.ms[j] += ms[j];
        dst.del[j] += useDel ? del[j] : 0.;
        dst.ioRe[j] += iore[j];
        dst.ioIm[j] += ioim[j];
        dst.radRe[j] += radre[j];
        dst.radIm[j] += radim[j];
      }
    }

    // The distinct material groups of `rows`, ascending. `col < 0` (rows that
    // carry no group) collapses to the single group -1, which is also what a
    // job with no global material model produces.
    void rowGroups(const float *rows, int stride, int n, int col, std::vector<int> &gs) {
      gs.clear();
      if (col < 0 || col >= stride) {
        gs.push_back(-1);
        return;
      }
      for (int i = 0; i < n; ++i)
        gs.push_back(static_cast<int>(rows[static_cast<std::size_t>(i) * stride + col]));
      std::sort(gs.begin(), gs.end());
      gs.erase(std::unique(gs.begin(), gs.end()), gs.end());
    }

    // Copy the rows of `rows` belonging to group `g` into `dst` (contiguous).
    // With `col < 0` every row belongs to the single group -1.
    int selectGroup(const float *rows, int stride, int n, int col, int g, std::vector<float> &dst) {
      dst.clear();
      if (col < 0 || col >= stride) {
        dst.assign(rows, rows + static_cast<std::size_t>(n) * stride);
        return n;
      }
      int m = 0;
      for (int i = 0; i < n; ++i) {
        const float *r = rows + static_cast<std::size_t>(i) * stride;
        if (static_cast<int>(r[col]) != g)
          continue;
        dst.insert(dst.end(), r, r + stride);
        ++m;
      }
      return m;
    }

  }  // namespace

  //==========================================================================
  // 6. THE EVALUATOR PASS, ON THE CONCATENATED ARGUMENT LIST
  //
  // `nf` functionals of ONE converged fit -- in the two-track maker the
  // candidate MASS, the vertex DCA and the two whitened BEAM-LINE pulls --
  // share every block and every step record and differ only in the scalar
  // weight they give a block,
  //
  //     w_{b,k} = s_{b,k} sqrt(v_{b,k}/sq2_b) / sigma_k ,
  //
  // and in the ionization sign they carry. Every primitive above reads the
  // weight only through the PRODUCT `w tau`, so the `nf` functionals are the
  // same primitive evaluated on the concatenated list { w_{b,k} tau_j }, and
  // the weight-independent half of the block -- the pooling by global index,
  // the row gather, `sq2`, the Moliere step parameters and the `gshape_elec`
  // row they interpolate, `makeRadSpectrum`, the per-group row selections --
  // is paid ONCE instead of `nf` times.
  //
  // EXACTNESS. `phi_{aU}(tau) = phi_U(a tau)` is an identity, and the products
  // `w_{b,k} tau_j` are formed by the same expression, in the same
  // association, that a single-functional call forms them with. Each
  // functional's exponents are therefore BITWISE what one call per functional
  // would have written -- there is no reference functional and no rescaled
  // grid, only a longer list of arguments.
  //
  // THE POOLING is unchanged and is still `cf_track_resolution.extract`'s: for
  // each material family in (10, 11), for each DISTINCT global parameter index
  // among the resolution entries of that family, pool v_b over the entries,
  // gather the step rows whose index VALUE matches, and form the block's
  // scalar standardized weight. What is now per functional is only what was
  // always per functional: the pooled variance and its sign, the
  // standardization, and the `want*` flags.
  //==========================================================================
  namespace {

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
      // dropped exactly as the single-functional path drops it; the others
      // carry on. `alive` is what keeps that per functional rather than per
      // call, here and at the missing-step-rows exit below.
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
      // The step records and the (global index, family) arrays are SHARED by
      // contract; they are read from the first usable functional.
      const TrackInput &sh = in[kref];

      // The GAUSSIAN families. 8/9 are the hit blocks; 16 is the beam-line
      // (luminous-region) block of the two-track maker, which is a Gaussian
      // noise block of exactly the same kind -- a fixed covariance, no Landau
      // channel -- so it belongs here and NOT with the material families.
      // Counting it keeps the reader's closure
      // `sum_g (vqms + vqio) + vgauss/sigma^2 == 1` exact with the beam rows
      // on. No maker that does not register family 16 is affected.
      for (int k = 0; k < nf; ++k) {
        if (!alive[k])
          continue;
        for (int i = 0; i < in[k].nres; ++i)
          if (in[k].resfamily[i] == 8 || in[k].resfamily[i] == 9 || in[k].resfamily[i] == 16)
            out[k].vgauss += in[k].resvarv[i];
      }

      // The `want*` flags are per functional. Shared work is done when ANY
      // functional asks for it and is scattered only to those that do; in the
      // maker all four agree, so nothing extra is computed in practice.
      bool anyGroups = false, anyDelta = false, anyGroupDelta = false;
      int nresmax = 0;
      for (int k = 0; k < nf; ++k) {
        if (!alive[k])
          continue;
        anyGroups = anyGroups || in[k].wantGroups;
        anyDelta = anyDelta || in[k].wantDelta;
        anyGroupDelta = anyGroupDelta || (in[k].wantDelta && in[k].wantGroupDelta);
        nresmax = std::max(nresmax, in[k].nres);
      }

      std::vector<unsigned int> gs;
      std::vector<int> sel;
      std::vector<float> blk;   // the block's step rows, contiguous
      std::vector<float> qs;    // the block's [scale, nsteps] pairs
      std::vector<float> rblk, rspec;
      // per-group scratch (untouched unless some functional wants groups)
      std::vector<int> gsteps;
      std::vector<float> gblk, grblk, grspec;
      std::vector<double> w2all;
      // per-functional scratch
      std::vector<double> vpool(nf, 0.), vsig(nf, 0.), vqblk(nf, 0.);
      std::vector<int> nent(nf, 0);
      std::vector<int> kact;     // the functionals ACTIVE on this block
      std::vector<double> wact;  // their weights, packed to match `kact`
      std::vector<double> bufA, bufB;
      // The nuclear-elastic family needs every block's weight before it can
      // pair a step with its leg's ionisation block, so the weights are
      // collected here (ascending in global index, as `gs` is) and the family
      // is formed after the pooling pass.
      bool anyNucel = false;
      for (int k = 0; k < nf; ++k)
        anyNucel = anyNucel || (alive[k] && in[k].wantNucel && in[k].nucel != nullptr);
      std::vector<std::vector<std::pair<unsigned int, double>>> nucWms(anyNucel ? nf : 0),
          nucWio(anyNucel ? nf : 0);

      for (int fam = 10; fam <= 11; ++fam) {
        const StepRows &rows = (fam == 10) ? sh.ms : sh.ioni;
        gs.clear();
        for (int i = 0; i < nresmax; ++i)
          if (sh.resfamily[i] == fam)
            gs.push_back(sh.resglobidx[i]);
        std::sort(gs.begin(), gs.end());
        gs.erase(std::unique(gs.begin(), gs.end()), gs.end());

        for (unsigned int g : gs) {
          kact.clear();
          for (int k = 0; k < nf; ++k) {
            vpool[k] = 0.;
            vsig[k] = 0.;
            nent[k] = 0;
            if (!alive[k])
              continue;
            const TrackInput &I = in[k];
            for (int i = 0; i < I.nres; ++i)
              if (I.resfamily[i] == fam && I.resglobidx[i] == g) {
                vpool[k] += I.resvarv[i];
                // variance-weighted sign of the block's influence; +1
                // everywhere when the caller supplies none, so
                // `vsig == vpool > 0`.
                vsig[k] += (I.ressgn != nullptr ? double(I.ressgn[i]) : 1.) * I.resvarv[i];
                ++nent[k];
              }
            if (vpool[k] > 0.)
              kact.push_back(k);
          }
          if (kact.empty())
            continue;
          sel.clear();
          for (int i = 0; i < rows.n; ++i)
            if (rows.idx[i] == g)
              sel.push_back(i);
          if (sel.empty()) {
            // A registered block with no step rows. The offline extractor
            // DROPS such a track (`ok = False`), so the flag is exported and
            // the reader drops it identically rather than silently keeping a
            // track whose model is missing a block. Only the functionals that
            // HAVE this block are dropped: one pooling no variance here would
            // never have seen it on its own.
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
          const int ns = static_cast<int>(sel.size());
          blk.resize(static_cast<std::size_t>(ns) * rows.stride);
          for (int i = 0; i < ns; ++i)
            std::memcpy(&blk[static_cast<std::size_t>(i) * rows.stride],
                        rows.v + static_cast<std::size_t>(sel[i]) * rows.stride, rows.stride * sizeof(float));
          const int nact = static_cast<int>(kact.size());
          const std::size_t nslab = static_cast<std::size_t>(kNSlabFam) * nact * kNTau;

          if (fam == 10) {
            for (int k : kact)
              ++out[k].nblockms;
            double sq2 = 0.;
            for (int i = 0; i < ns; ++i)
              sq2 += blk[static_cast<std::size_t>(i) * rows.stride + 5];
            if (!(sq2 > 0.))
              continue;
            wact.resize(nact);
            for (int a = 0; a < nact; ++a) {
              const int k = kact[a];
              wact[a] = std::sqrt(vpool[k] / sq2) / in[k].sigma;
              // the block's share of the STANDARDIZED variance under the
              // fit's Q
              vqblk[k] = vpool[k] / (in[k].sigma * in[k].sigma);
              if (anyNucel && in[k].wantNucel)
                nucWms[k].emplace_back(g, wact[a]);
            }

            rowGroups(blk.data(), rows.stride, ns, anyGroups ? sh.ms.groupCol : -1, gsteps);
            // THE WHOLE BLOCK, ONCE, FOR EVERY FUNCTIONAL. The flat exponent
            // is the same evaluation whether or not the block straddles two
            // materials, so it is hoisted out of the split; a block that does
            // not straddle then hands the SAME numbers to the group slot,
            // which is what keeps `sum_g == flat` exact to the last bit.
            bufA.assign(nslab, 0.);
            msBlockMulti(blk.data(), rows.stride, ns, wact.data(), nact, slab(bufA, kSlabMs, nact),
                         anyDelta ? slab(bufA, kSlabDel, nact) : nullptr);
            for (int a = 0; a < nact; ++a) {
              const int k = kact[a];
              const double *Sm = slab(bufA, kSlabMs, nact) + static_cast<std::size_t>(a) * kNTau;
              const double *Sd = slab(bufA, kSlabDel, nact) + static_cast<std::size_t>(a) * kNTau;
              for (int j = 0; j < kNTau; ++j) {
                out[k].S.ms[j] += Sm[j];
                if (in[k].wantDelta)
                  out[k].S.del[j] += Sd[j];
              }
            }
            if (anyGroups && gsteps.size() == 1) {
              for (int a = 0; a < nact; ++a) {
                const int k = kact[a];
                if (!in[k].wantGroups)
                  continue;
                const double *Sm = slab(bufA, kSlabMs, nact) + static_cast<std::size_t>(a) * kNTau;
                const double *Sd = slab(bufA, kSlabDel, nact) + static_cast<std::size_t>(a) * kNTau;
                GroupExponents &ge = groupEntry(out[k].groups, gsteps[0]);
                ge.vqms += vqblk[k];
                for (int j = 0; j < kNTau; ++j) {
                  ge.S.ms[j] += Sm[j];
                  if (in[k].wantDelta && in[k].wantGroupDelta)
                    ge.S.del[j] += Sd[j];
                }
              }
            } else if (anyGroups) {
              // A block that straddles two materials: the split is a second
              // pass. The carve is a RATIO over the WHOLE block (see
              // delCarveFactor) and does not depend on the weight, so it is
              // computed once and reused by every group and every functional.
              const double carve = anyGroupDelta ? delCarveFactor(blk.data(), rows.stride, ns) : 0.;
              for (int gg : gsteps) {
                const int mg = selectGroup(blk.data(), rows.stride, ns, sh.ms.groupCol, gg, gblk);
                if (mg <= 0)
                  continue;
                bufB.assign(nslab, 0.);
                msBlockMulti(gblk.data(), rows.stride, mg, wact.data(), nact, slab(bufB, kSlabMs, nact), nullptr);
                if (anyGroupDelta)
                  delBlockCarvedMulti(gblk.data(), rows.stride, mg, wact.data(), nact, slab(bufB, kSlabMs, nact),
                                      carve, slab(bufB, kSlabDel, nact));
                // the group's share of the block's fit-Q variance is its share
                // of `sq2`, which for multiple scattering is the `thp2` column.
                double sq2g = 0.;
                for (int i = 0; i < mg; ++i)
                  sq2g += gblk[static_cast<std::size_t>(i) * rows.stride + 5];
                for (int a = 0; a < nact; ++a) {
                  const int k = kact[a];
                  if (!in[k].wantGroups)
                    continue;
                  GroupExponents &ge = groupEntry(out[k].groups, gg);
                  ge.vqms += vqblk[k] * (sq2g / sq2);
                  addIntoSlice(ge.S, bufB, a, nact, in[k].wantDelta && in[k].wantGroupDelta);
                }
              }
            }
          } else {
            for (int k : kact)
              ++out[k].nblockioni;
            qs.clear();
            for (int i = 0; i < sh.qsc.n; ++i)
              if (sh.qsc.idx[i] == g) {
                qs.push_back(sh.qsc.v[static_cast<std::size_t>(i) * sh.qsc.stride]);
                qs.push_back(sh.qsc.v[static_cast<std::size_t>(i) * sh.qsc.stride + 1]);
              }
            const double sq2 =
                ioniSq2(blk.data(), rows.stride, ns, qs.empty() ? nullptr : qs.data(), static_cast<int>(qs.size() / 2));
            if (!(sq2 > 0.))
              continue;
            wact.resize(nact);
            for (int a = 0; a < nact; ++a) {
              const int k = kact[a];
              const double sgnblk = (vsig[k] < 0.) ? -1. : 1.;
              wact[a] = in[k].ioniSign * sgnblk * (std::sqrt(vpool[k] / sq2) / in[k].sigma);
              vqblk[k] = vpool[k] / (in[k].sigma * in[k].sigma);
              if (anyNucel && in[k].wantNucel)
                nucWio[k].emplace_back(g, wact[a]);
            }

            rowGroups(blk.data(), rows.stride, ns, anyGroups ? sh.ioni.groupCol : -1, gsteps);
            bufA.assign(nslab, 0.);
            ioniBlockMulti(blk.data(), rows.stride, ns, wact.data(), nact, slab(bufA, kSlabIoRe, nact),
                           slab(bufA, kSlabIoIm, nact));
            for (int a = 0; a < nact; ++a) {
              const int k = kact[a];
              const double *Re = slab(bufA, kSlabIoRe, nact) + static_cast<std::size_t>(a) * kNTau;
              const double *Im = slab(bufA, kSlabIoIm, nact) + static_cast<std::size_t>(a) * kNTau;
              for (int j = 0; j < kNTau; ++j) {
                out[k].S.ioRe[j] += Re[j];
                out[k].S.ioIm[j] += Im[j];
              }
            }
            if (anyGroups && gsteps.size() == 1) {
              for (int a = 0; a < nact; ++a) {
                const int k = kact[a];
                if (!in[k].wantGroups)
                  continue;
                const double *Re = slab(bufA, kSlabIoRe, nact) + static_cast<std::size_t>(a) * kNTau;
                const double *Im = slab(bufA, kSlabIoIm, nact) + static_cast<std::size_t>(a) * kNTau;
                GroupExponents &ge = groupEntry(out[k].groups, gsteps[0]);
                ge.vqio += vqblk[k];
                for (int j = 0; j < kNTau; ++j) {
                  ge.S.ioRe[j] += Re[j];
                  ge.S.ioIm[j] += Im[j];
                }
              }
            } else if (anyGroups) {
              // per-step variance, so a group's share of `sq2` is a plain sum
              ioniStepVar(blk.data(), rows.stride, ns, qs.empty() ? nullptr : qs.data(),
                          static_cast<int>(qs.size() / 2), w2all);
              double w2tot = 0.;
              for (double x : w2all)
                w2tot += x;
              for (int gg : gsteps) {
                const int mg = selectGroup(blk.data(), rows.stride, ns, sh.ioni.groupCol, gg, gblk);
                if (mg <= 0)
                  continue;
                bufB.assign(nslab, 0.);
                ioniBlockMulti(gblk.data(), rows.stride, mg, wact.data(), nact, slab(bufB, kSlabIoRe, nact),
                               slab(bufB, kSlabIoIm, nact));
                double w2g = 0.;
                if (sh.ioni.groupCol >= 0) {
                  for (int i = 0; i < ns; ++i) {
                    if (static_cast<int>(blk[static_cast<std::size_t>(i) * rows.stride + sh.ioni.groupCol]) == gg)
                      w2g += w2all[i];
                  }
                }
                for (int a = 0; a < nact; ++a) {
                  const int k = kact[a];
                  if (!in[k].wantGroups)
                    continue;
                  GroupExponents &ge = groupEntry(out[k].groups, gg);
                  ge.vqio += (w2tot > 0.) ? vqblk[k] * (w2g / w2tot) : 0.;
                  addIntoSlice(ge.S, bufB, a, nact, in[k].wantDelta && in[k].wantGroupDelta);
                }
              }
            }

            // The radiative channel of the SAME block: same weight, same sign.
            // The join is on the global index VALUE and never on the row
            // position -- the radiative rows are 1:1 with `msmoliv`, not with
            // `ioniurbanv`.
            if (sh.rad.n > 0 && sh.radspec != nullptr && sh.radvgrid != nullptr && sh.radnv > 1) {
              rblk.clear();
              rspec.clear();
              for (int i = 0; i < sh.rad.n; ++i) {
                if (sh.rad.idx[i] != g)
                  continue;
                const float *r = sh.rad.v + static_cast<std::size_t>(i) * sh.rad.stride;
                rblk.insert(rblk.end(), r, r + sh.rad.stride);
                const float *sp = sh.radspec + static_cast<std::size_t>(i) * 2 * sh.radnv;
                rspec.insert(rspec.end(), sp, sp + 2 * sh.radnv);
              }
              if (!rblk.empty()) {
                const int nr = static_cast<int>(rblk.size() / sh.rad.stride);
                // The radiative channel is two thirds of the whole CF cost
                // (`nsteps x nv x ntau` trigonometry on a spectrum that has to
                // be built first), so sharing the spectrum across the
                // functionals matters most here.
                rowGroups(rblk.data(), sh.rad.stride, nr, anyGroups ? sh.rad.groupCol : -1, gsteps);
                bufA.assign(nslab, 0.);
                radBlockMulti(rblk.data(), sh.rad.stride, nr, rspec.data(), sh.radvgrid, sh.radnv, wact.data(), nact,
                              slab(bufA, kSlabRadRe, nact), slab(bufA, kSlabRadIm, nact));
                for (int a = 0; a < nact; ++a) {
                  const int k = kact[a];
                  const double *Re = slab(bufA, kSlabRadRe, nact) + static_cast<std::size_t>(a) * kNTau;
                  const double *Im = slab(bufA, kSlabRadIm, nact) + static_cast<std::size_t>(a) * kNTau;
                  for (int j = 0; j < kNTau; ++j) {
                    out[k].S.radRe[j] += Re[j];
                    out[k].S.radIm[j] += Im[j];
                  }
                }
                if (anyGroups && gsteps.size() == 1) {
                  for (int a = 0; a < nact; ++a) {
                    const int k = kact[a];
                    if (!in[k].wantGroups)
                      continue;
                    const double *Re = slab(bufA, kSlabRadRe, nact) + static_cast<std::size_t>(a) * kNTau;
                    const double *Im = slab(bufA, kSlabRadIm, nact) + static_cast<std::size_t>(a) * kNTau;
                    Exponents &gsl = groupSlot(out[k].groups, gsteps[0]);
                    for (int j = 0; j < kNTau; ++j) {
                      gsl.radRe[j] += Re[j];
                      gsl.radIm[j] += Im[j];
                    }
                  }
                } else if (anyGroups) {
                  // The spectra ride along with their rows, so the subset has
                  // to be taken on BOTH arrays with one index walk.
                  for (int gg : gsteps) {
                    grblk.clear();
                    grspec.clear();
                    for (int i = 0; i < nr; ++i) {
                      const float *r = rblk.data() + static_cast<std::size_t>(i) * sh.rad.stride;
                      if (sh.rad.groupCol >= 0 && static_cast<int>(r[sh.rad.groupCol]) != gg)
                        continue;
                      grblk.insert(grblk.end(), r, r + sh.rad.stride);
                      const float *sp = rspec.data() + static_cast<std::size_t>(i) * 2 * sh.radnv;
                      grspec.insert(grspec.end(), sp, sp + 2 * sh.radnv);
                    }
                    if (grblk.empty())
                      continue;
                    bufB.assign(nslab, 0.);
                    radBlockMulti(grblk.data(), sh.rad.stride, static_cast<int>(grblk.size() / sh.rad.stride),
                                  grspec.data(), sh.radvgrid, sh.radnv, wact.data(), nact,
                                  slab(bufB, kSlabRadRe, nact), slab(bufB, kSlabRadIm, nact));
                    for (int a = 0; a < nact; ++a) {
                      const int k = kact[a];
                      if (!in[k].wantGroups)
                        continue;
                      addIntoSlice(groupSlot(out[k].groups, gg), bufB, a, nact,
                                   in[k].wantDelta && in[k].wantGroupDelta);
                    }
                  }
                }
              }
            }
          }
        }
      }

      if (anyNucel) {
        for (int k = 0; k < nf; ++k)
          if (alive[k] && in[k].wantNucel && in[k].nucel != nullptr)
            nucelFunctional(sh, in[k], nucWms[k], nucWio[k], out[k]);
      }
    }

  }  // namespace

  void trackExponents(const TrackInput &in, TrackResult &out) { trackExponentsImpl(&in, 1, &out); }

  void trackExponents(const TrackInput *in, int nfunc, TrackResult *out) { trackExponentsImpl(in, nfunc, out); }

  //==========================================================================
  const std::string &modelTag() {
    static const std::string tag =
        "cvhcf/1 tau=stride4of448<=8 ms:elecTmax=1,elecEdge=1,fineG=1,snapYmax=0,wviSplit=0 "
        "ioni:kokoulin=0,a3=1,exc=1,tmaxScale=1 rad:cf_brems_exact del:tcut=0.35MeV,tmaxcap=50MeV "
        "tab=cvhcf_gshape_elec_v1";
    return tag;
  }

  std::string modelTag(bool withNucel) {
    return withNucel ? modelTag() + " nuc=" + nucelTableId() : modelTag();
  }

}  // namespace cvhcf
