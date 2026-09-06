# The resolution-CF export of the CVH makers

*(2026-09-05. Companion to `TrackPropagation/Geant4e/interface/CvhCfExponents.h`,
which carries the physics documentation; this file is the maker-side contract:
what is written, what it replaces, and how to turn the old export off.)*

## What the makers write

The offline resolution model describes a fitted track's q/p error — and a
two-track candidate's mass error — by the characteristic function of the sum of
its independent process-noise blocks,

    phi_z(t) = phi_hit(t) * exp( S_ms(t) + S_ioni(t) + S_rad(t) + S_del(t) )

with `z` the residual standardized by the fit's own sigma. Since 2026-09-05 the
four exponents are computed **inside the makers**, in the `doRes` pass, and
written on a fixed 64-point tau grid:

| branch (single track) | branch (two track) | type | meaning |
|---|---|---|---|
| `cfqop_ms`      | `cfmass_ms`      | 64 float | Moliere multiple scattering, real |
| `cfqop_ioni_re` | `cfmass_ioni_re` | 64 float | Re of the Urban / exact-delta ionization exponent |
| `cfqop_ioni_im` | `cfmass_ioni_im` | 64 float | Im of the same |
| `cfqop_rad_re`  | `cfmass_rad_re`  | 64 float | Re of the brems + pair exponent |
| `cfqop_rad_im`  | `cfmass_rad_im`  | 64 float | Im of the same |
| `cfqop_del`     | `cfmass_del`     | 64 float | discrete delta-ray recoil, minus the Moliere carve |
| `cfqop_vgf`     | `cfmass_vgf`     | float | the GAUSSIAN share of the functional's variance |
| `cfqop_ok`      | `cfmass_ok`      | bool | 0 iff a registered block had no step rows |
| `cfqop_nblock`  | `cfmass_nblock`  | int | MS + ionization blocks entering the exponents |
| `cfqop_npooled` | `cfmass_npooled` | int | of those, how many pooled more than one leg |

plus, once, in the **runtree**:

| branch | meaning |
|---|---|
| `cftau`   | the 64 tau values the exponents are sampled at |
| `cfmodel` | the switch configuration and shape-table id that produced them |

### The per-material-group split (2026-09-06, `exportCfGroupExponents`)

Every step-level exponent is **linear in the step's material amount at fixed
composition** — Moliere's `chi_c^2` and `Omega_0` carry the step length, the
Urban `a_1, a_2, a_3` carry it, the radiative mean number of emissions carries
it, and `xi` of the delta channel carries it — so with the fit's influence
weights held fixed

    S_f(tau; k) = S_f^fixed(tau) + sum_g A(k_g) S_{f,g}(tau),   A(k) = exp(k)

is exact, and `k = 0` reproduces the flat exponents above. `k_g` is the
parmtype-15 material-group amount the propagator already applies to every step
(`matStepFact`) and the hit chi2 already floats. That is what replaces the four
unphysical `k_hit / k_ms / k_ioni / k_rad` knobs.

| branch (single track) | branch (two track) | type | meaning |
|---|---|---|---|
| `cfqop_grp`         | `cfmass_grp`         | n int16 | material-group ids, ASCENDING |
| `cfqop_grp_ms`      | `cfmass_grp_ms`      | n*64 float | `S_ms` of each group, row-major (group, tau) |
| `cfqop_grp_ioni_re` | `cfmass_grp_ioni_re` | n*64 float | |
| `cfqop_grp_ioni_im` | `cfmass_grp_ioni_im` | n*64 float | |
| `cfqop_grp_rad_re`  | `cfmass_grp_rad_re`  | n*64 float | |
| `cfqop_grp_rad_im`  | `cfmass_grp_rad_im`  | n*64 float | |
| `cfqop_grp_del`     | *(absent)*           | n*64 float | the q/p functional's model uses `S_del`; the mass functional's reference (`build_pairs_tt`) does not |
| `cfqop_grp_closure` | `cfmass_grp_closure` | float | `max_j |sum_g S_g - S|` over the exported families, normalized by `max_j |S|` |

The per-candidate group count is the length of `cf*_grp`; there is no pointer
branch, because an entry knows its own length. The group ids index the
material-group file the runtree names.

`cf*_grp_closure` is exported so a file can be AUDITED rather than trusted. It
is float64 round-off by construction — the flat and per-group paths are the
same per-step sums associated differently, and a block whose steps all belong
to one group (the common case) runs the primitives ONCE and adds the same
doubles to both, so `sum_g` is bitwise equal there.

The DELTA-RECOIL family needs one care: its carve factor
`clip(v_delta/v_moliere, 0, 0.5)` is a RATIO over the whole block, so the split
is

    S_del,g = delta_exponent(steps_g) - carve(steps_ALL) * S_ms,g

which sums back exactly; a per-group carve would not.

**Cost**: ~22 live groups of 42 on a J/psi-gun candidate, 5 families, 64 points,
float32 = **~27 kB/candidate** against 1.4 kB flat. Hence
`exportCfGroupExponents`, default **False**. The raw rows are written on
purpose: a rank-16 PCA of the tau axis costs 6.8 kB at a relative error of
4.2e-7 (ms) to 1.2e-3 (rad), but the basis should be fitted to the data rather
than guessed, and this is the export that lets that happen.

**The group is now on every step record.** `msmoliv` has carried it at column 9
since the global material model landed; since 2026-09-06 `ioniurbanv` and
`radstepv` carry it as their LAST column, so the offline
`groups.pair_ioni_rows` heuristic (take the `n_ioni` Moliere rows of largest
areal density) is no longer needed. Strides changed accordingly:

| array | before | after |
|---|---|---|
| `ioniurbanv` | 11, or 13 with `CVH_IONI_EXACTDELTA` | **12 / 14**, new scalar branch `ioniurbanstride` |
| `radstepv` | 11 (`radstepstride`) | **12** |

Every existing column index is unchanged. **A reader that hard-codes 11/13 will
mis-parse a new file** — read `ioniurbanstride` / `radstepstride`.

### The hit-class blocks and their Gaussian shares (2026-09-06)

| branch | type | meaning |
|---|---|---|
| `reshitcls` | n_res int16 | hit CLASS of each resolution entry, −1 for material |
| `resinfcovhit` | float | `sum_b v_b` over the HIT families (parmtype 8/9) alone |
| `cfqop_hitcls` / `cfmass_hitcls` | int16 | classes present, ascending |
| `cfqop_hitv` / `cfmass_hitv` | float | `v_c / sigma^2` of that class |

The 18 classes are the C++ image of
`calibration_studies/resolution/hitres_classes.py:class_of`, in the SAME index
order (0–3 `pix_x_q0..3`, 4–7 `pix_y_q0..3`, 8–17 `str_N1_lo` … `str_N5_hi`).
The index, not the name, is what a file stores, so the order is part of the
format and must not be reordered.

They feed the hit half of the same construction,

    v_i(eps) = v_other,i + sum_c H(eps_c) v_{c,i},   Re S += -0.5 v_i tau^2

with `v_other = vgf - sum_c v_c` formed offline.

**`cfmass_vgf` and `resinfcov` do NOT change.** The two-track maker now
registers parmtype-8/9 blocks (`exportHitResBlocks`, default True) — it never
did, which is why the per-hit-class parameters could not be fitted at all — but
their variances go into the NEW `resinfcovhit`, not into `resinfcov`. That is
deliberate: `cfmass_vgf = (sigma_m^2 - resinfcov)/sigma_m^2` must keep meaning
the TOTAL Gaussian share (hits + beamspot + pointing), because that is what the
self-consistent-sigma correction's `a_i = (1 + f_hit) sigma_i/m_i` uses and what
every cache built before the hit blocks existed assumes. Folding them in would
have been a silent physics change.

The two-track `dVs` feed nothing but the influence export, so registering a
block cannot move the fit; and the parmtype-8/9 corrections are deliberately
NOT applied to the two-track hit covariance (the single-track maker scales `iV`
by `exp(corparms)`), because that WOULD move it.

In the single-track tree `resinfcov` has always included the hit blocks and
still does; `resinfcovhit` is the same sum on its own, and `sum_c cfqop_hitv`
equals `cfqop_vgf` exactly.

### The two legs' reference momentum covariance (2026-09-06)

| branch | type | meaning |
|---|---|---|
| `Jpsi_covrefmom` | 21 float | UPPER TRIANGLE, row-major, of the symmetric 6x6 covariance of `(q/p, lambda, phi)` for mu+ (0–2) then mu− (3–5) |
| `Jpsi_jacrefmom` | 6 float | `dm/d(state)` in the same order |
| `Jpsi_qoprefplus` / `Jpsi_qoprefminus` | float | `q/p` at the reference |
| `Jpsi_sigmarelplus` / `Jpsi_sigmarelminus` | float | `sqrt(C_ll)/|q/p_l|` |
| `Jpsi_rhomom` | float | leg–leg correlation of `d ln p` (NOT of `q/p`: they differ by `sign(q+ q−)`) |
| `Jpsi_fang` | float | `1 - (J_kappa C J_kappa^T)/sigma_m^2`, the ANGULAR share of the mass variance |

33 floats, 132 B/candidate, 0.16 % of the slim record, always on. The order is
(plus, minus) rather than the internal leg order, so no consumer has to know
`muchargearr`, and `J C J^T` reproduces `Jpsi_sigmamass^2` — exported precisely
so the matrix can be checked rather than trusted.

This is what makes the second-order corrections of
`resolution/oddmoment/MASSCFTERM_SPEC.md` truth-free on DATA. The Jensen term
needs `A = sigma_rel1^2 + sigma_rel2^2` and `B = 2 rho sigma_rel1 sigma_rel2`,
and the closed form `1.5 (sigma_m/m)^2` misses exactly `f_ang`; until now all
three had to be taken from an MC measurement (rho = 0, f_ang = 0.11, the closed
form 4.6 % high).

### The pre-FSR gen mass (2026-09-06)

| branch | meaning |
|---|---|
| `Jpsigenpre_mass` | mass of the hard-process resonance: `|pdgId|` in `genResonancePdgIds` at status 62 (last copy), falling back to 22 (first copy); −99 if absent |
| `Jpsigenpre_status` | which copy was used |
| `Jpsigenpre_masslep` | the status-746 (pre-Photos) lepton pair, the independent cross-check; −99 when the event did not radiate |
| `Jpsigen_massdressed` | the matched bare pair plus every prompt status-1 photon within dR < 0.1 of either muon |

`Jpsigen_mass` is and stays the POST-FSR pair. In the DY UL16 sample no muon
carries `fromHardProcessBeforeFSR`; what exists is the Z at status 22/62 (masses
identical in 4000/4000 events) and, in the 59 % of events that radiated, the
status-746 pair, with `m(mumu, 746) == m(Z, 62)` to an RMS of 4e-6 GeV. So the
status-62 resonance IS the pre-FSR mass and it exists in every event
(`calibration_studies/zchannel/README.md`). With these branches
`zfsr_kernel.py` runs off the production's own pairs cache and inherits the
analysis selection exactly, instead of a separate FWLite pass whose selection
has to be kept in step by hand.

`genResonancePdgIds` defaults to `{23, 443, 100443, 553, 100553, 200553}`. A
J/psi from a B decay has no status-22/62 copy, so `Jpsigenpre_mass` is −99 there
— that is information, not a failure. The gen container is an
`edm::View<reco::Candidate>`, which does not expose `GenStatusFlags`, so the
photons' `isPrompt` is asked for through a `dynamic_cast` and simply not
required when the cast fails.

`cfqop_*` and `cfmass_*` are **different functionals of the same blocks** — they
differ in the standardization sigma (`sqrt(refCov(0,0))` against
`Jpsi_sigmamass`) and in the sign the ionization and radiative weights carry
(the track charge against −1 for both legs and both charges). They are named
apart so that nothing can read one as the other.

`vgf` also means different things in the two trees, and both are what the
offline cache called `vgf`:

* single track — `sum_b v_b` over the hit families (parmtype 8/9) divided by
  `refCov(0,0)`;
* two track — `(sigma_m^2 − resinfcov)/sigma_m^2`, i.e. hits + beamspot +
  pointing by construction, because the two-track maker does not register hit
  blocks at all.

## Why in the maker

The block weights are `sqrt(v_b/sq2)/sigma` with `v_b` the fit's own influence
coefficients. They do not exist until the fit has converged — which is *why* the
raw Geant4 step records had to be exported in the first place. Doing the
evaluation in the `doRes` pass, from the same flat arrays the tree would have
carried, removes both halves of the cost:

Measured (compressed TREE bytes per entry; the runtree is a fixed 4.54 MB per
file and is quoted separately):

| | single track | two track |
|---|---|---|
| the exponents | 1 453 | 1 419 |
| the raw step records they replace | 165 377 | 329 466 |
| `exportStepRecords=True` | 230 180 | 448 055 |
| `exportStepRecords=False` | **64 803** | **118 589** |
| per-candidate evaluation time | 41.7 ms | 95.8 ms |
| the offline extraction it replaces | 2.2 s | 2.2 s |

Projected over 20M J/psi + 8M Z + 10M Upsilon = 38M dimuon candidates: 17.0 TB
today, **4.51 TB** with the switch off, and ~0.84 TB once `fillGradsFactored`
also replaces the dense Hessian (which is 4.10 TB of the 4.51). The exponents
themselves are 54 GB.

i.e. the switch removes 72 % of the two-track tree, and the evaluation is 53x faster
than the offline extraction of the same object — inside a fit that costs O(1 s)
per candidate, so 4–10 % of it. Two thirds of that time is the radiative
channel (`nsteps x nv x ntau` trigonometric evaluations); a further factor ~2
is available there from the half-angle identity and was not taken, because the
exactness of the port is worth more at this size.

What is LEFT after the switch is dominated by `hesspackedv` (59 kB/track,
108 kB/candidate) — a separate problem with a separate solution already in the
tree (`fillGradsFactored`).

The inputs are the **export arrays**, not the propagator's internal logs, so
the in-maker pooling by global parameter index is identical to the offline
`cf_track_resolution.extract` join by construction, including the ~1/3 of blocks
that pool two legs (each with its own `ioniqscalev` entry) under one index.

## The tau grid

64 points, the **stride-4 subset of the offline `np.linspace(0, 14, 448)`**
truncated at 8, i.e. `tau_i = 4*i*14/447`, `tau_63 = 7.8926`.

A subset rather than a fresh `linspace(0, 8, 64)` because every point is then
also a point of the offline grid, so the in-maker exponent can be compared
against the reference *at the same argument* with no interpolation in between.
It is also the grid the 2026-09-04 compression study measured
(`calibration_studies/resolution/cfcompress/gridtest.py`, row `stride4/t<=8.0`):
refitting the unbinned mass likelihood on it moves alpha by −1.6e-8, about
1/1000 of the full-sample statistical error, and every other fitted parameter by
less than 0.001 sigma.

## The slimming switch

```
cmsRun runCvhResClosure.py  ... exportStepRecords=False
cmsRun runCvhJpsiGenMC.py   ... exportStepRecords=False
```

`exportStepRecords` (default **False** since 2026-09-06; it was True) governs
the RAW per-step export that the exponents are built from:

* dropped when False: `ioniurbanidx`/`ioniurbanv`, `msmoliidx`/`msmoliv`,
  `radstepidx`/`radstepv`/`radstepspecv`, `reseigv`, `resinfv`, `resinfbv`;
* **kept regardless**: `reseigidx`, `reshitidx`, `reshitcls`, `resinfvarv`,
  `resinfcov`, `resinfcovhit`, `ioniqscaleidx`/`ioniqscalev`, `radvgrid`,
  `radstepstride`, `radstepnv`, `ioniurbanstride`, and every `cf*` branch.

The default flipped because every production since 2026-09-05 set it to False
explicitly (`production/config_jpsimc20M.sh`, `config_dymc8p5M.sh`) and the
exponents it feeds are validated against the offline reference, so leaving the
default at True meant a new driver silently wrote 72 % of a two-track tree in
records nothing reads. `exportStepRecords=True` still reproduces the old output
(modulo the two appended group columns above).

The keep list is what a reader still needs and what costs nothing: `resinfvarv`
+ `reseigidx` + `reshitidx` are the per-block variance shares and their hit
association, which is what gives every parmtype-8/9 block its own measured
per-hit CF class (`--hitmode class`); `resinfcov` is the coverage check
(`|resinfcov/refCov(0,0) − 1| < 5e-3`) that the offline extractor drops tracks
on; `ioniqscalev` records the CGF substitution factor the fit applied. All of
them are O(100 B).

`resinfbv` (the 5×5 `B_b = M_b dV_b^{1/2}` per block, ~6 kB/candidate) goes with
the raw records rather than the keep list: its only consumer is
`cf_mass_likelihood.leg_exponents`, the single-track pairing route that the
two-track maker's own `cfmass_*` export supersedes.

`exportCfExponents` (default True) governs the exponents themselves. Turning
both off reproduces a pre-2026-09-05 tree exactly.

**Do not turn `exportStepRecords` off on a sample whose model has not been
validated yet.** Its purpose is to make a DIFFERENT model derivable from the
same files; once it is off, changing the model means re-running the fit.

## One switch that must agree on both sides

`cvhcf` exports the **Kokoulin-OFF** ionization model (`cfmodel` says
`ioni:kokoulin=0`), which is what the offline production builds: the shard
runner sets `CVH_IONI_KOKOULIN=0` explicitly, for a term that costs >7x and
moves the model by ~1e-3. The PROPAGATOR may still have `IoniKokoulin=True` --
that is the correction in the fit's Q matrix and in `ioniurbanv`'s `gsig2`, a
different object.

`cf_track_resolution.IONI_KOKOULIN` defaults to ON when the environment does
not say otherwise, so an offline extraction run without `CVH_IONI_KOKOULIN=0`
builds a DIFFERENT model from the one the maker exported. That is exactly what
`cfmodel` is in the file for: compare it before mixing caches.

## Reading it back

`calibration_studies/resolution/cf_inmaker.py` turns these branches into the
same npz caches `cf_track_resolution --extract` and
`cf_mass_likelihood --pairs-tt` used to write:

```bash
python3 cf_inmaker.py extract --files '<glob>/globalcor_resclosure_0.root' --cache runs/cf_trackres_<tag>.npz
python3 cf_inmaker.py pairs   --files '<glob>/globalcor_0.root'            --cache runs/cf_masspairs_<tag>.npz
python3 cf_inmaker.py closure --cache runs/cf_trackres_<tag>.npz           # the reference closure, on the 64-point grid
```

## Validating it

`calibration_studies/resolution/cvhcf_validate.py` compares `cvhcf` against the
offline reference on real records, at two levels — per pooled block (a failure
is a formula error) and per whole track (a failure that block level passes is a
pooling error):

```bash
source /work/submit/david_w/ZMass/mfs/.venv/bin/activate
export CMSSW_SRC=/work/submit/david_w/ZMass/CMSSW_15_0_19_patch2_dev/src
./cxx/build_cvhcf.sh                      # compiles the maker's own .cc into a ctypes shim
python3 cvhcf_validate.py --file <globalcor.root> --ntracks 40
```

The requirement is 1e-6 absolute on the exponent. Measured on the reference
smokes: **1.2e-11** (single track, 40 tracks / 515 blocks) and **1.6e-11**
(two track, 28 candidates / 734 blocks) per block and per track, and 9.5e-7 on
the floats the maker wrote — which IS the float32 storage floor at those |S|,
every other family being 100x below it.

End to end (`cvhcf_e2e_260905.sh`, 1891 gun ditrack candidates): the unbinned
mass likelihood moves by `d(alpha) = 2e-10` in its 1e-3 units, i.e. 2e-9 of its
own sigma, and the single-track even closure is identical to every printed
digit. See `Documents/Resolution/NOTES.md`, 2026-09-05.
