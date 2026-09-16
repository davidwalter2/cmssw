#ifndef TrackPropagation_Geant4e_CvhCfExponents_h
#define TrackPropagation_Geant4e_CvhCfExponents_h

// IN-FIT EVALUATION OF THE PER-CANDIDATE RESOLUTION-CF EXPONENTS.
//
// WHAT THIS IS FOR
// ----------------
// The offline resolution model describes a fitted track's q/p error -- and a
// two-track candidate's mass error -- by the CHARACTERISTIC FUNCTION of the
// sum of its independent process-noise blocks:
//
//     phi_z(t) = phi_hit(t) * exp( S_ms(t) + S_ioni(t) + S_rad(t) + S_del(t) )
//
// with z the residual standardized by the fit's own sigma. Each S is a
// compound-Poisson log-CF built from the Geant4 STEP RECORDS of the block and
// the block's scalar standardized weight `wstd = sqrt(v_b/sq2)/sigma`.
//
// Building the exponents offline instead would mean exporting every step
// record -- `ioniurbanv`, `msmoliv`, `radstepv`/`radstepspecv`, `reseigv`,
// `resinfv` -- at 430 kB and 2.2 s per candidate, i.e. 16 TB and 24k
// core-hours at the 40M candidates the full calibration needs, to carry what
// is in the end 6 x 64 floats.
//
// The weights are known only AFTER the fit converges (they are the fit's own
// influence coefficients), so this runs in the makers' doRes pass, from the
// same flat arrays the tree would have carried. Reading the EXPORT ARRAYS
// rather than the propagator's internal logs is deliberate: it makes the
// in-maker pooling identical BY CONSTRUCTION to the offline `extract()` join,
// including the ~1/3 of blocks that pool two legs under one global index.
//
// FIDELITY
// --------
// This is a port of, and is validated against, the offline reference:
//   * S_ms   -- `cf_track_resolution.ms_step_exponent` (Moliere compound
//               Poisson with the universal shape G(tau; ymax), the Z:Z^2
//               nuclear/electron split and the electron kinematic ceiling),
//               with `cf_ms_exact.moliere_params` / `gshape_elec`;
//   * S_ioni -- `cf_track_resolution.ioni_step_exponent` (Urban regimes 0/1
//               plus the exact spin-1/2 and spin-0 knock-on channels), which
//               is ALREADY implemented in `cvhcgf` for the in-fit Fisher
//               weight and is reused here rather than duplicated;
//   * S_rad  -- `cf_brems_exact.rad_exponent`, on the propagator's own
//               tabulated brems/pair spectra (`cvhcgf::makeRadSpectrum`);
//   * S_del  -- `cf_delta_ray.delta_step_exponent` minus `carve_factor`
//               times S_ms of the same block (the discrete delta-ray recoil
//               REPLACES part of Moliere's continuous Z(Z+1), it does not
//               add to it).
//
// THE PORTED SWITCH CONFIGURATION is the offline production default and is
// recorded in `modelTag()` (which the makers write into the runtree, so a file
// says which model produced its exponents):
//     MS_ELEC_TMAX = 1, MS_ELEC_EDGE = 1, MS_FINE_G = 1, MS_SNAP_YMAX = 0,
//     MS_WVI_SPLIT = 0, IONI_KOKOULIN = 0, IONI_A3_SCALE = IONI_EXC_SCALE =
//     IONI_TMAX_SCALE = 1, CF_DELTA = 1 (tcut 0.35 MeV, tmax cap 50 MeV).
// The diagnostic knobs the offline module carries (MS_WVI_SPLIT, MS_ELEC_EDGE
// 0/0.5/2/3, the a3/excitation gauges, Kokoulin) are NOT ported: they exist to
// SIZE a systematic on a small sample, which is exactly the case where the
// offline extractor is still affordable. Porting them would put four dead
// branches in the inner loop of a 40M-candidate production.
//
// THE MOLIERE SHAPE TABLES are not rebuilt here. `cf_ms_exact._build_elec_
// tables` is 5.8e8 evaluations of J0(x)-1 on a 141 x 1600 x 2560 grid; both
// the ~20 s of startup and the dependence on WHICH J0 the toolchain ships
// (libstdc++'s is not Cephes, which is what scipy uses) are avoided by loading
// the reference's own table from `data/cvhcf_gshape_elec_v1.bin`, written by
// `data/make_cvhcf_gshape_tables.py`. The C++ then does the identical linear
// interpolation on the identical numbers.

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace cvhcf {

  // THE EXPORT GRID.  64 points, tau in [0, 7.8926], the STRIDE-4 SUBSET of
  // the offline `TG = np.linspace(0, 14, 448)`.
  //
  // Why a subset and not `linspace(0, 8, 64)`: it makes the validation exact.
  // Every one of these tau is a point of the offline grid, so the in-maker
  // exponent can be compared against the reference at the SAME argument with
  // no interpolation standing between them.  The truncation it implies is
  // harmless: refitting the unbinned mass likelihood on this grid moves alpha
  // by -1.6e-8 (-1.6e-5 in the 1e-3 units the fit reports), i.e. 1/1000 of the
  // full-sample statistical error, and every other fitted parameter by less
  // than 0.001 sigma.
  //
  // The grid is written into the runtree next to the exponents, so no consumer
  // has to reconstruct it.
  constexpr int kNTau = 64;
  constexpr int kTauStride = 4;      // of the offline 448-point grid
  constexpr int kTauFullN = 448;
  constexpr double kTauFullMax = 14.0;
  const double *tauGrid();

  // One block's / one track's exponents, on `tauGrid()`.
  struct Exponents {
    std::array<double, kNTau> ms{};
    std::array<double, kNTau> del{};
    std::array<double, kNTau> ioRe{};
    std::array<double, kNTau> ioIm{};
    std::array<double, kNTau> radRe{};
    std::array<double, kNTau> radIm{};
    void clear();
  };

  // A jagged export array: `idx[i]` is the global parameter index of row `i`,
  // `v` holds `n * stride` floats. Exactly the (msmoliidx, msmoliv) /
  // (ioniurbanidx, ioniurbanv) / (radstepidx, radstepv) pairs of the makers.
  struct StepRows {
    const unsigned int *idx = nullptr;
    const float *v = nullptr;
    int n = 0;       // rows
    int stride = 0;  // floats per row
    // Column holding the step's MATERIAL GROUP (the parmtype-15 index of the
    // global material model), or -1 when the rows carry none. `msmoliv` has
    // it at column 9; `ioniurbanv` and `radstepv` carry it as their LAST
    // column.  Only read when `TrackInput::wantGroups`.
    int groupCol = -1;
  };

  // Everything one track (or one two-track candidate) needs.
  struct TrackInput {
    // Resolution entries, parallel arrays of length `nres`, in the order the
    // maker pushed them (i.e. aligned with reseigidx / resinfvarv).
    const unsigned int *resglobidx = nullptr;
    const int *resfamily = nullptr;  // 8/9 = hit, 10 = MS, 11 = ionization
    const float *resvarv = nullptr;  // resinfvarv, the block's v_b
    int nres = 0;

    StepRows ms;    // msmoliidx / msmoliv, stride 10 (>= 8 tolerated)
    StepRows ioni;  // ioniurbanidx / ioniurbanv, stride 11 or 13
    StepRows qsc;   // ioniqscaleidx / ioniqscalev, stride 2: [scale, nsteps]
    StepRows rad;   // radstepidx / radstepv, stride RADSTEP_STRIDE = 11
    // radstepspecv (2*nv floats per rad row) and the shared v grid.
    const float *radspec = nullptr;
    const float *radvgrid = nullptr;
    int radnv = 0;

    // The standardization: sqrt(refCov(0,0)) for the q/p functional,
    // Jpsi_sigmamass for the candidate-mass one. Taken from the FLOAT the
    // tree carries, so the in-maker weight is the one the offline reader
    // would have recovered.
    double sigma = 0.;

    // The sign the ionization (and radiative) step weight carries.
    //   q/p functional : the track CHARGE. `ioniurbanv`'s cs = E/p^3 is
    //                    positive for every track and the physical map is
    //                    d(q/p) = q cs dE, so the exponent's weight is
    //                    charge-signed (cf_track_resolution `chg * wstd`).
    //   mass functional: -1 for BOTH legs and BOTH charges -- an energy loss
    //                    on either muon can only LOWER the pair mass
    //                    (cf_mass_likelihood.IONI_SGN).
    double ioniSign = 1.;

    // OPTIONAL PER-BLOCK SIGN, parallel to `resvarv` (length `nres`).
    //
    // `ioniSign` is ONE number because the q/p (and the candidate-mass)
    // functional's influence has the same sign on every ionization block. A
    // general linear functional of the same fit does not: one whitened
    // component of the PER-HIT (complement) residual vector has its own
    // signed influence coefficient `W[r0, k]` on each block's qop row, and
    // the ionization CF is not even in its weight (an energy loss can only go
    // one way), so the sign has to travel with the block.
    //
    // When set, a pooled block's weight is
    //     ioniSign * sign(sum_i s_i v_i) * sqrt(vpool/sq2) / sigma
    // over the entries `i` sharing the block's global index -- the
    // variance-weighted sign, which for the (common) single-entry block is
    // just `s_i`. NULL reproduces the scalar behaviour bit for bit (all
    // `s_i = +1`, so the sum is `vpool > 0`).
    const float *ressgn = nullptr;

    bool wantDelta = true;  // the discrete delta-ray recoil family

    // SPLIT THE EXPONENTS BY MATERIAL GROUP as well as accumulating the flat
    // ones.  This is what lets the offline fit float the AMOUNT of material
    // per group instead of four per-family k knobs: every step-level exponent
    // is linear in the step's material amount at fixed composition, so
    //
    //     S_f(tau; k) = S_f^fixed(tau) + sum_g A(k_g) S_{f,g}(tau)
    //
    // is exact with the fit's influence weights held fixed.  Off by default:
    // ~22 live groups per candidate
    // multiply the 1.4 kB flat export by ~20.
    bool wantGroups = false;
    // Include the delta-recoil family in the per-group split.  The q/p
    // functional's model uses `S_del`; the mass functional's reference
    // (`cf_mass_likelihood.build_pairs_tt`) does not, so the two-track maker
    // leaves it out of the per-group arrays and saves a sixth of them.
    bool wantGroupDelta = false;
  };

  // One material group's share of a track's exponents.
  struct GroupExponents {
    int group = -1;
    Exponents S;
    // THE FIT'S OWN Q VARIANCE for this group, in units of `sigma^2` -- i.e.
    // the group's share of the standardized functional's variance under the
    // MODEL THE FIT USED (Rossi's `thp2` for multiple scattering, `ioniSq2`
    // for ionization, and nothing for the radiative or delta channels,
    // because the fit's `Q` has neither).  It is what the Gaussian chi2 the
    // whole exercise is measured against actually assumes, and it cannot be
    // recovered from the exponents: those carry the MODEL's (full Moliere)
    // second moment, which is 14 % larger.  `sum_g (vqms + vqio) + vgauss/
    // sigma^2 == 1` by construction, which is the closure a reader checks.
    double vqms = 0.;
    double vqio = 0.;
  };

  struct TrackResult {
    // False when a registered MS/ionization block has NO step rows under its
    // global index. The offline extractor DROPS such a track (`ok = False`),
    // so the flag is exported and the reader drops it identically rather than
    // silently keeping a track whose model is missing a block.
    bool ok = false;
    double vgauss = 0.;  // sum of v_b over the GAUSSIAN families (8, 9, and 16 = beam line)
    int nblockms = 0, nblockioni = 0, npooled = 0;
    Exponents S;
    // The per-group split of `S`, ascending in `group`; empty unless
    // `TrackInput::wantGroups`.  `sum_g groups[i].S == S` to float64
    // round-off -- they are the same per-step sums associated differently,
    // NOT two models -- which is the validation gate the makers report.
    std::vector<GroupExponents> groups;
  };

  // THE ENTRY POINT. Pools by global parameter index exactly as
  // `cf_track_resolution.extract` does and accumulates the four families.
  void trackExponents(const TrackInput &in, TrackResult &out);

  // THE MULTI-FUNCTIONAL ENTRY POINT: `nfunc` functionals of the SAME fit in
  // ONE pass over the step records.
  //
  // The two-track maker forms four linear functionals of one converged fit --
  // the candidate MASS, the vertex DCA, and the two whitened BEAM-LINE pulls.
  // They share every block and every step record and differ only in the
  // per-block scalar weight they give it,
  //
  //     w_{b,k} = s_{b,k} sqrt(v_{b,k}/sq2_b) / sigma_k ,
  //
  // and in the ionization sign they carry. EVERY exponent primitive here
  // depends on (weight, tau) ONLY through the product w tau -- the Moliere
  // shape is read at `sqrt(chi_a^2) w tau`, the delta and ionization channels
  // at `gs w tau`, the radiative one at `cs w tau` -- so the k functionals are
  // the SAME primitive evaluated on the CONCATENATED argument list
  // { w_{b,k} tau_j }_{k,j}. Everything that does not depend on the weight --
  // the pooling by global index, the block gather, `sq2`, the Moliere step
  // parameters and the `gshape_elec` row they interpolate, the radiative
  // spectra `makeRadSpectrum` builds, the per-group row selections -- is then
  // done ONCE instead of `nfunc` times.
  //
  // This is EXACT, not an approximation: phi_{aU}(tau) = phi_U(a tau) is an
  // identity, and the products `w_{b,k} tau_j` are formed by the same
  // expression, in the same association, that the single-functional path forms
  // them with. Each functional's exponents are therefore BITWISE what one call
  // per functional would have produced; the concatenation changes WHICH points
  // are evaluated, never how.
  //
  // `in[k]` may differ in `resvarv`, `ressgn`, `nres`, `sigma`, `ioniSign`,
  // `wantDelta`, `wantGroups` and `wantGroupDelta`. The STEP RECORDS (`ms`,
  // `ioni`, `qsc`, `rad`, `radspec`, `radvgrid`, `radnv`) and the
  // (`resglobidx`, `resfamily`) arrays are shared and are read from the first
  // usable entry; passing entries that disagree on them is a caller error.
  // `out` must have `nfunc` elements. `nfunc == 1` is bit-identical to the
  // single-functional entry point, which is implemented as exactly that call.
  void trackExponents(const TrackInput *in, int nfunc, TrackResult *out);

  //------------------------------------------------------------------------
  // THE PER-BLOCK PRIMITIVES.
  //
  // Public because the validation runs through them: the offline reference is
  // written per block (`ms_step_exponent(steps, wstd, tau)` and friends), so a
  // like-for-like comparison needs the same entry points, and going through
  // the full `trackExponents` would confound a formula error with a pooling
  // error. They are built into a ctypes shim for that comparison, exactly as
  // `cvhcgf`'s primitives already are.
  //
  // Every one of them ACCUMULATES into its output (length kNTau), and every
  // one takes the rows as the FLOATS the tree carries, so no conversion
  // convention can differ between the maker and the reader.

  // One pooled multiple-scattering block: `rows` are its `msmoliv` records.
  // `Sms` gets the Moliere exponent; `Sdel`, when non-null, gets the discrete
  // delta-ray recoil MINUS the carve (i.e. the family as the cache stores it).
  void msBlock(const float *rows, int stride, int n, double wstd, double *Sms, double *Sdel);

  // `cf_track_resolution.ioni_sq2`: the block variance the FIT used.
  // `qsc` are the [scale, nsteps] pairs of the legs sharing the block (may be
  // null), `nqsc` their count.
  double ioniSq2(const float *rows, int stride, int n, const float *qsc, int nqsc);

  // One pooled ionization block at an ALREADY SIGNED standardized weight.
  void ioniBlock(const float *rows, int stride, int n, double wstdSigned, double *Sre, double *Sim);

  // THE BLOCK'S DELTA-RECOIL CARVE FACTOR, `clip(v_delta/v_moliere, 0, 0.5)`.
  // `msBlock` computes and applies its own; a per-group split has to reuse
  // the WHOLE block's, because the carve is a RATIO: per-group ratios would
  // not sum back to the block's exponent, whereas
  //   S_del,g = delta(steps_g) - carve(steps_ALL) * S_ms,g
  // does, exactly.
  double delCarveFactor(const float *rows, int stride, int n);

  // The delta-recoil family at an externally supplied carve factor. `Sms`
  // must be the exponent of the SAME rows (the group's own, not the block's).
  void delBlockCarved(
      const float *rows, int stride, int n, double wstd, const double *Sms, double carve, double *Sdel);

  // The radiative channel of the same block: `rows` are its `radstepv`
  // records, `spec` its `radstepspecv` rows (2*nv floats each), `vgrid` the
  // shared v grid. Same weight and sign as the ionization block.
  void radBlock(const float *rows,
                int stride,
                int n,
                const float *spec,
                const float *vgrid,
                int nv,
                double wstdSigned,
                double *Sre,
                double *Sim);

  // Provenance: the ported switch configuration and the shape-table id, for
  // the runtree. Stable for a given model; changes only when the model does.
  const std::string &modelTag();

  // Load the Moliere shape tables. Called automatically on first use; exposed
  // so a job can fail at configuration time rather than mid-event, and so the
  // standalone validation driver can point at a table explicitly.
  // Thread-safe, idempotent; throws cms::Exception if the file is unusable.
  void loadShapeTables(const std::string &path = std::string());

}  // namespace cvhcf

#endif
