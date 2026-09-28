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
//     phi_z(t) = phi_hit(t) * exp( S_ms + S_ioni + S_rad + S_kx + S_kj (+ S_nuc) )
//
// with z the residual standardized by the fit's own sigma. Each S is a
// compound-Poisson log-CF built from the Geant4 STEP RECORDS of the blocks
// and each block's scalar standardized weight `wstd = sqrt(v_b/sq2)/sigma`.
//
// Building the exponents offline instead would mean exporting every step
// record -- `ioniurbanv`, `msmoliv`, `radstepv`/`radstepspecv`, `reseigv`,
// `resinfv` -- at 430 kB per candidate, i.e. 16 TB at the 40M candidates the
// full calibration needs, to carry what is in the end a few x 64 floats.
//
// The weights are known only AFTER the fit converges (they are the fit's own
// influence coefficients), so this runs in the makers' doRes pass, from the
// same flat arrays the tree would have carried. Reading the EXPORT ARRAYS
// rather than the propagator's internal logs is deliberate: it makes the
// in-maker pooling identical BY CONSTRUCTION to the offline `extract()` join,
// including the ~1/3 of blocks that pool two legs under one global index.
//
// ONE IMPLEMENTATION OF EVERY CHANNEL
// -----------------------------------
// The physics is the offline ROW API, calibration_studies/resolution/
// cf_rows.py, ported function by function ("THE ROW FUNCTIONS" below):
//   * scattering   -- `ms_rows` -> cf_track_resolution.ms_step_exponent
//                     (Moliere with G(tau; ymax), the Z:Z^2 split, the
//                     electron ceiling, which stops at the e- production
//                     threshold under the knock-on joint law);
//   * ionisation   -- `ioni_rows` -> ioni_step_exponent (Urban regimes 0/1,
//                     the exact Bethe-Bloch knock-on spectrum with Kokoulin's
//                     correction for muons, Moller and Bhabha for e+-);
//   * radiation    -- `rad_rows` -> cf_brems_exact.rad_exponent (spectra
//                     refined `refine_spectra`, the exact 1/p map, the
//                     primary's recoil against the photon);
//   * knock-on     -- `knockon_rows` -> cf_knockon (the hard collision as one
//                     event: exact map and joint deflection);
//   * nuclear elastic, hadrons only -- NucelExponents below;
// and the fit-level assembly is `cf_rows.fit_families` (`trackExponents`).
// The clean-propagation closure evaluates its model through these same
// functions (calibration_studies/resolution/cvhcf_rows.py), so it tests the
// model the makers export; cxx/gate_cvhcf_rows.py measures the port against
// the Python reference: bit for bit on the knock-on quadrature and the
// spectrum refinement, to rounding elsewhere.
//
// THE SWITCHES are the reference modules' globals, carried by `RowConfig`
// and recorded in `modelTag()` (which the makers write into the runtree, so a
// file says which model produced its exponents).  The scattering channel's
// diagnostic knobs (MS_WVI_SPLIT, MS_ELEC_EDGE 0/0.5/2/3, MS_SNAP_YMAX) are
// not ported: the port is MS_ELEC_TMAX = MS_ELEC_EDGE = MS_FINE_G = 1,
// MS_SNAP_YMAX = MS_WVI_SPLIT = 0, and the offline binding refuses anything
// else.
//
// THE MOLIERE SHAPE TABLES are not rebuilt here. `cf_ms_exact._build_elec_
// tables` is 5.8e8 evaluations of J0(x)-1 on a 141 x 1600 x 2560 grid; both
// the ~20 s of startup and the dependence on WHICH J0 the toolchain ships
// are avoided by loading the reference's own table from
// `data/cvhcf_gshape_elec_v2.bin`, written by `data/make_cvhcf_gshape_tables.py`.
// The C++ then does the identical linear interpolation on the identical
// numbers.  The J0 and K1 the other channels evaluate are Cephes' (scipy's),
// ported (`besselJ0`, `besselK1`).

#include <array>
#include <cstddef>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <utility>
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

  class NucelMixtures;  // the nuclear-elastic kernels, below

  //------------------------------------------------------------------------
  // THE ROW FUNCTIONS  (offline reference: calibration_studies/resolution/
  // cf_rows.py and what it calls).
  //
  // Every process-noise channel as a function of the exported Geant4 step
  // RECORDS -- `ioniurbanv` (ionisation, knock-on), `msmoliv` (scattering),
  // `radstepv` + spectrum (radiation) -- and a list of ENTRIES, each
  // (record index rid, q/p weight wq, angular weight wb, fraction frac):
  //
  //   wq    z per unit q/p, the track charge included (the records' q/p maps
  //         are for a positive charge);
  //   wb    z per radian, the length of the projected angular weight
  //         (w_lambda, w_phi / cos lambda): the kicks are isotropic in 2D, so
  //         a channel sees only that length;
  //   frac  the entry's share of the record's material / rate (1/NSUB for a
  //         sub-step, 1 for a whole record).
  //
  // The SAME functions serve the clean-propagation closure (rows at the exact
  // transport weights of each step, expanded into sub-steps) and the fit-level
  // CF of the makers (rows at the fit's block influence weights, one entry
  // per record, `trackExponents`), so the closure tests exactly the model the
  // makers export.  They take an ARBITRARY tau array (the makers use
  // `tauGrid()`) and the records as doubles (the makers widen the floats the
  // tree carries, which is exact), and ACCUMULATE into their outputs.
  //
  // The model's switches are the reference modules' globals, carried by
  // `RowConfig` so that the offline side can hand its module state over
  // (the shim) and the makers run `productionRowConfig()`.
  struct RowConfig {
    // cf_knockon: the joint law of loss and deflection (KNOCKON_JOINT), the
    // exact energy -> q/p map (QOP_EXACT; also the radiative channel's map),
    // the knock-on quadrature's nodes per decade for ordinary and thin steps
    // (KNOCKON_NPERDEC, KNOCKON_NPERDEC_THIN), the thin-step threshold as a
    // fraction of the row set's largest xi (KNOCKON_THIN), the joint law's
    // lower limit [MeV] (KNOCKON_TCUT) and the p' floor (PMIN_FRAC)
    bool knockonJoint = true;
    bool qopExact = true;
    // cf_knockon.QOP_LOG: the exact map's q/p variable on a logarithmic scale,
    // |q/p|_ref ln(|q/p|/|q/p|_ref), whose per-collision changes q ln(p/p')/p
    // add over successive collisions exactly (repeated radiation compounds
    // without error); with qopExact only
    bool qopLog = false;
    int knockonNPerDec = 40;
    int knockonNPerDecThin = 10;
    double knockonThin = 1e-2;
    double knockonTcut = 0.99e-3;
    double pminFrac = 1e-3;
    // cf_brems_exact: the refinement of the exported radiative spectra
    // (RAD_NSUB)
    int radNsub = 4;
    // cf_track_resolution: Geant4's Kokoulin correction to the muon knock-on
    // spectrum (IONI_KOKOULIN, its lower limit IONI_KOKOULIN_TCUT [MeV]; its
    // quadrature takes the knock-on nodes, knockonNPerDec) and the diagnostic
    // gauges IONI_A3_SCALE, IONI_EXC_SCALE, IONI_TMAX_SCALE
    double ioniKokoulin = 1.0;
    double ioniKokoulinTcut = 0.0;
    double ioniA3Scale = 1.0;
    double ioniExcScale = 1.0;
    double ioniTmaxScale = 1.0;
    // cf_nucel_exact: the nuclear-elastic recoil family (NUCEL_RECOIL) and
    // the joint law of one collision's deflection and recoil (NUCEL_JOINT,
    // with the recoil only)
    bool nucelRecoil = true;
    bool nucelJoint = true;
    bool knockonActive() const { return knockonJoint || qopExact; }
  };
  // The offline modules' defaults.
  const RowConfig &productionRowConfig();

  // scipy.special's j0 and k1 (Cephes), evaluated by the joint knock-on law
  // and the heavy species' photon-angle law.
  double besselJ0(double x);
  double besselK1(double x);

  // The projectile's mass [GeV] of a radiative record, sqrt(E^2 - p^2)
  // snapped to the nearest charged species (e, mu, pi, K, p) -- the float32
  // (E, p) of a maker's record carry it only to a few MeV
  // (cf_brems_exact.species_mass).
  double speciesMass(double E, double p);

  // The ionisation channel (cf_rows.ioni_rows): Urban excitations and each
  // record's knock-on law -- 1/E^2 (regime 1), Bethe-Bloch spin 1/2 and 0
  // with Kokoulin for muons (2/3), Moller (4) and Bhabha (5) -- mapped
  // LINEARLY into q/p, each row at its own weight wq[i] (n rows, stride
  // >= 11; >= 13 for regimes 2-5).  Complex, centred.
  void ioniRows(const double *tau,
                int nt,
                const double *rows,
                int stride,
                int n,
                const double *wq,
                const RowConfig &cfg,
                double *Sre,
                double *Sim);

  // Multiple scattering (cf_rows.ms_rows): one isotropic 2D Moliere kick per
  // entry, the entry's share of the record's material (`msmoliv`, stride >=
  // 8) at angular weight wb; entries with wb <= 0 carry nothing.  Under the
  // knock-on joint law (cfg.knockonJoint) the electron term stops at the e-
  // production threshold cfg.knockonTcut: the collisions above it are the
  // knock-on channel's.  Real.
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
              double *S);

  // The exported radiative spectra (48 points, brems then pair per row, on
  // `vg`) on a grid `nsub` times finer (cf_brems_exact.refine_spectra):
  // log-log between points, the last interval continued to the step's
  // kinematic end point (the mass `speciesMass` takes from the record) with
  // an area-preserving straddling node.  `vf` gets
  // the (nv-1)*nsub+1 fine nodes, `specf` n rows of 2*vf.size().
  void refineSpectra(const double *recs,
                     int rstride,
                     int n,
                     const double *spec,
                     const double *vg,
                     int nv,
                     int nsub,
                     std::vector<double> &vf,
                     std::vector<double> &specf);

  // Radiation (cf_rows.rad_rows / cf_brems_exact.rad_exponent) on REFINED
  // spectra (`spec`: 2*nv per record on `vg`): each emission the joint event
  // of its q/p change (the exact 1/p map under `exactQop`, linear otherwise)
  // and the primary's recoil against the photon at angular weight wb
  // (ModifiedTsai for e+-, ModifiedMephi for every heavier species); frac
  // scales the step length.  Complex, centred on the linear mean.
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
               double *Sim,
               bool qopLog = false);

  // The hard knock-on collision exactly (cf_rows.knockon_rows /
  // cf_knockon.knockon_rows) on the IONISATION rows: the collisions above the
  // e- production threshold as ONE event each -- the energy loss with the
  // exact 1/p map and the deflection jointly -- as the correction of the
  // ionisation channel's linear map, per record's law (regimes 2-5).  The
  // scattering channel's electron term stops at the threshold when the joint
  // law is on (msRows), so the deflection above it is carried here alone.
  //   kAll    INT dN [e^{i a X} J - e^{i a T}]
  //   kMap    INT dN [e^{i a X} - e^{i a T}]           (zero unless qopExact)
  //   kJoint  INT dN e^{i a X} (J - 1)                  (zero unless knockonJoint)
  // X = T_eff under qopExact, J = J0(b theta) under knockonJoint.  The
  // thin-step node density is set by the largest xi of the `n` rows passed.
  enum class KnockonPart { kAll = 0, kMap = 1, kJoint = 2 };
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
                   double *Sim);

  // One track's (or candidate's) families, on `tauGrid()`: cf_rows.
  // fit_families' Sms (real), Sio, Srad, Skx and Skj (complex).
  struct Exponents {
    std::array<double, kNTau> ms{};
    std::array<double, kNTau> ioRe{};
    std::array<double, kNTau> ioIm{};
    std::array<double, kNTau> radRe{};
    std::array<double, kNTau> radIm{};
    // the hard knock-on collision on the ionisation rows: the exact energy ->
    // q/p map (kx, a nonzero mean: the Jensen excess) and the joint law of
    // loss and deflection (kj); zero with the switches off or without
    // `TrackInput::wantKnockon`
    std::array<double, kNTau> kxRe{};
    std::array<double, kNTau> kxIm{};
    std::array<double, kNTau> kjRe{};
    std::array<double, kNTau> kjIm{};
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

    // SPLIT THE EXPONENTS BY MATERIAL GROUP as well as accumulating the flat
    // ones.  This is what lets the offline fit float the AMOUNT of material
    // per group instead of four per-family k knobs: every step-level exponent
    // is linear in the step's material amount at fixed composition, so
    //
    //     S_f(tau; k) = S_f^fixed(tau) + sum_g A(k_g) S_{f,g}(tau)
    //
    // is exact with the fit's influence weights held fixed.  Off by default:
    // ~22 live groups per candidate multiply the flat export by ~20.  With
    // it, the flat families are formed as the sum of the group parts, so the
    // two agree bit for bit.
    bool wantGroups = false;

    // THE NUCLEAR-ELASTIC FAMILY (see NucelExponents).  Computed only when
    // `wantNucel` and `nucel` are set; then `msmat` and `mspdg` must be
    // parallel to `ms` (one G4Material index and one species PDG code per
    // row, 0 = no elastic channel) and `rad` must be parallel to `ms` as the
    // makers write it (checked: a mismatch throws).  Split by the step's
    // material group (`ms.groupCol`) as well when `wantGroups`.
    bool wantNucel = false;
    const int *msmat = nullptr;
    const int *mspdg = nullptr;
    NucelMixtures *nucel = nullptr;

    // THE KNOCK-ON FAMILIES kx, kj (under the model's switches).  A
    // functional that does not use them -- the per-hit components, whose
    // offline reference carries them only on request -- skips their cost.
    bool wantKnockon = true;

    // The model's switches (`productionRowConfig()` when null).  Shared by
    // every functional of a multi-functional pass.
    const RowConfig *rowConfig = nullptr;
  };

  //------------------------------------------------------------------------
  // THE NUCLEAR-ELASTIC FAMILY (hadElastic), for hadron tracks.
  //
  // A hadron also takes `hadElastic` collisions: rare (N ~ 0.02-0.1 per
  // track), large (25-100 mrad at 3 GeV) isotropic kicks that no Gaussian
  // block and no Moliere exponent carries.  They enter ONLY the resolution
  // CF, as compound-Poisson families -- never the fit's process noise Q or its
  // weights, which stay frozen (a Gaussian variance for a rare large kick
  // would outweigh multiple scattering and inflate every hadron's errors).
  //
  //   angular:  S_ang(tau) = sum_s N_s ( g_{m(s)}(p_s; wb_s tau) - 1 )
  //   recoil:   S_rec(tau) = sum_s N_s ( h_{m(s)}(p_s; wq_s tau) - 1 )
  //   joint:    S_jnt(tau) = sum_s N_s J_{m(s)}(p_s; wb_s tau, wq_s tau)
  //
  // over the `msmoliv` rows s, with N_s = mu_{m(s)}(p_s) xg_s the expected
  // collisions of the step (rate mu in cm^2/g, xg in g/cm^2), g the projected
  // single-collision CF of the deflection and h the CF of the kinetic energy
  // given to the recoiling nucleus [MeV].  One collision's CF is
  // E[J0(b theta) e^{i a X(dE)}]; the two families carry E[J0(b theta)] +
  // E[e^{i a dE}] - 1, and J is the difference -- the recoil is a function of
  // the deflection (two-body kinematics; on hydrogen a pion loses ~15 % of its
  // momentum at ~37 deg), and X is the exact 1/p map of the loss under
  // `RowConfig::qopExact`.  m(s) is the step's material with ONE TARGET PER
  // ELEMENT: a compound's kernel is the rate-weighted mixture of its
  // elements' kernels (hydrogen in the tracker composites is ~3x wider in
  // angle and ~12x harder in recoil than carbon).
  //
  //   * wb_s: in the fit, the MS block's own weight sqrt(v_b / sum thp2) /
  //     sigma -- a collision is an isotropic 2D kick at a point in the step,
  //     like a Moliere one, so it rides on the same weight.  Real, not
  //     centred (the kick's mean projection is zero).
  //   * wq_s: the ionisation weight of the SAME leg at the SAME step,
  //     wq_s = w_io(block) * cs_s * 1e-3, with w_io the leg's ionisation block
  //     weight (it carries `ioniSign`: the charge for q/p, -1 for a mass),
  //     cs_s = E/p^3 of the step and 1e-3 for the MeV of the recoil against
  //     the GeV of cs.  The ionisation block of step s is `rad.idx[s]` and cs_s
  //     is `rad.v[s][10]`: the radiative records are written one per Geant4
  //     step, PARALLEL to `msmoliv` (the Urban ionisation records are not --
  //     a step without a valid Urban record has none -- so they can never be
  //     paired by row).  Not centred: the Geant4e reference energy loss
  //     carries no elastic recoil, so no mean was subtracted.
  //
  // THE TABLES are `data/cvhcf_nucel_v2.bin` (writer
  // `data/make_cvhcf_nucel_tables.py`, whose docstring is the byte layout and
  // the evaluation rules implemented here): per (species, element, momentum
  // node) from Geant4's own species-specific elastic model and cross section,
  // with the joint (theta, dE) law of the same collisions.  Offline
  // reference: calibration_studies resolution/ksclosure/nucel/nucel_tables.py
  // (`Table.rows`, read from the same file), through which both the
  // clean-propagation closure (`cf_nucel_exact.rows_exponent`) and the
  // fit-level family (`ks_nucel_cf.step_family`) evaluate.
  //
  // SPECIES come from the maker's configured particle and the track charge
  // (p and pbar share a mass, not a model).  Muons and electrons have no
  // hadElastic: species 0, no family, zero cost.

  // The element composition of one Geant4 material: element Z, the element's
  // atomic mass [g/mole] and its mass fraction, as G4Material holds them.
  struct NucelComposition {
    std::vector<double> Z, A, W;
  };
  // Resolves a G4Material index to its composition; false if unknown.
  using NucelResolver = std::function<bool(int, NucelComposition &)>;

  // The table's PDG code for a Geant4 particle name ("pi+", "kaon-",
  // "proton", "anti_proton", ...); 0 for a particle without a tabulated
  // elastic channel (mu+-, e+-, ...).
  int nucelPdg(const std::string &g4ParticleName);

  // Load the kernel tables (automatic on first use; thread-safe, idempotent,
  // first caller wins as for the shape tables).  Throws cms::Exception if the
  // file is unusable.
  void loadNucelTables(const std::string &path = std::string());
  // The table id recorded in `modelTag(cfg, true)`.
  const std::string &nucelTableId();

  // The per-(species, material) mixture kernels, built on first use from the
  // resolver's composition and kept for the life of the object.  NOT
  // thread-safe: one per stream (the makers own one each).  The table itself
  // is process-global and read-only.
  struct NucelMix;
  class NucelMixtures {
  public:
    explicit NucelMixtures(NucelResolver resolver);
    ~NucelMixtures();
    NucelMixtures(const NucelMixtures &) = delete;
    NucelMixtures &operator=(const NucelMixtures &) = delete;
    // The mixture of species `pdg` in material `matIndex`.  Throws
    // cms::Exception when the material cannot be resolved or has an element
    // the table does not carry: a silently missing target would bias the
    // family, and every element of the tracker geometry is tabulated.
    const NucelMix &get(int pdg, int matIndex);
    std::size_t size() const { return cache_.size(); }

  private:
    NucelResolver resolver_;
    std::map<std::pair<int, int>, std::unique_ptr<NucelMix>> cache_;
  };

  // One track's (or candidate's) nuclear-elastic exponents.
  struct NucelExponents {
    std::array<double, kNTau> ang{};
    std::array<double, kNTau> recRe{};
    std::array<double, kNTau> recIm{};
    std::array<double, kNTau> jntRe{};
    std::array<double, kNTau> jntIm{};
    double N = 0.;  // expected collisions, sum_s N_s over the hadron rows
    void clear();
  };

  // THE ROW FUNCTION (`nucel_tables.Table.rows`), in the entry convention of
  // the other row functions.  Over `n` `msmoliv` rows (xg column 2, p column
  // 3) with material `mat[s]`, species `pdg[s]` (0: no channel) and mass
  // `mass[s]` [GeV], ACCUMULATES
  //   ang       per entry e (row rid[e], angular weight wb[e], share frac[e]
  //             of the row's collisions): N_s frac (g(|wb| tau) - 1);
  //   rec       per row at the signed weight wqRow[s] [z per MeV] under
  //             cfg.nucelRecoil: N_s (h(wq tau) - 1);
  //   jnt       per row at (|wbMid[s]|, wqRow[s]) under cfg.nucelJoint (with
  //             the recoil): N_s J, X = T_eff(dE) under cfg.qopExact;
  // on an arbitrary tau array; a zero weight carries no term.  Returns
  // sum_s N_s over the rows with a species.
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
                   double *jntIm);

  // One material group's share of a track's exponents.
  struct GroupExponents {
    int group = -1;
    Exponents S;
    // THE FIT'S OWN Q VARIANCE for this group, in units of `sigma^2` -- i.e.
    // the group's share of the standardized functional's variance under the
    // MODEL THE FIT USED (Rossi's `thp2` for multiple scattering, `ioniSq2`
    // for ionization, and nothing for the radiative and knock-on families,
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
    // `TrackInput::wantGroups`, and then `sum_g groups[i].S == S` exactly
    // (the flat families are formed as that sum), which is the closure the
    // makers report.
    std::vector<GroupExponents> groups;
    // The nuclear-elastic family; filled only under `TrackInput::wantNucel`
    // (`nucel` false otherwise), and its per-material-group split, ascending
    // in group, under `wantGroups` as well.  The split is its own list: a
    // step can carry an elastic term in a group where no other family does.
    bool nucel = false;
    NucelExponents nuc;
    std::vector<std::pair<int, NucelExponents>> nucGroups;
  };

  // THE ENTRY POINT: cf_rows.fit_families of one track (or candidate) from
  // the export arrays, with the pooling of `cf_track_resolution.extract`
  // (section 7 of the source describes both).
  void trackExponents(const TrackInput &in, TrackResult &out);

  // THE MULTI-FUNCTIONAL ENTRY POINT: `nfunc` functionals of the SAME fit in
  // one pass.  The two-track maker forms four linear functionals of one
  // converged fit -- the candidate MASS, the vertex DCA, and the two whitened
  // BEAM-LINE pulls -- which share every block and every step record and
  // differ only in the per-block weight they give it,
  //     w_{b,k} = s_{b,k} sqrt(v_{b,k}/sq2_b) / sigma_k ,
  // and in the ionisation sign they carry.  The weight-independent work (the
  // pooling, the widened records, the refined radiative spectra, the Moliere
  // step parameters) is done once; each functional's families are what its
  // own call would have produced.
  //
  // `in[k]` may differ in `resvarv`, `ressgn`, `nres`, `sigma`, `ioniSign`,
  // `wantGroups`, `wantNucel` and `wantKnockon`.  The STEP RECORDS (`ms`,
  // `ioni`, `qsc`, `rad`, `radspec`, `radvgrid`, `radnv`), the (`resglobidx`,
  // `resfamily`) arrays and `rowConfig` are shared and are read from the
  // first usable entry; passing entries that disagree on them is a caller
  // error.  `out` must have `nfunc` elements.
  void trackExponents(const TrackInput *in, int nfunc, TrackResult *out);

  // `cf_track_resolution.ioni_sq2`: the block variance the FIT used.
  // `qsc` are the [scale, nsteps] pairs of the legs sharing the block (may be
  // null), `nqsc` their count.
  double ioniSq2(const float *rows, int stride, int n, const float *qsc, int nqsc);

  // Provenance for the runtree: the model's switches, the shape table, and
  // the nuclear-elastic table id when that family is exported.
  std::string modelTag(const RowConfig &cfg = productionRowConfig(), bool withNucel = false);

  // Load the Moliere shape tables. Called automatically on first use; exposed
  // so a job can fail at configuration time rather than mid-event, and so the
  // standalone validation driver can point at a table explicitly.
  // Thread-safe, idempotent; throws cms::Exception if the file is unusable.
  void loadShapeTables(const std::string &path = std::string());

}  // namespace cvhcf

#endif
