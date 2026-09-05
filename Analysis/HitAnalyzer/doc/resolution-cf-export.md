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

| | before | after |
|---|---|---|
| per-candidate export | ~430 kB | ~1.6 kB |
| per-candidate offline extraction | 2.2 s | 0 |
| at 40M candidates | 16 TB, 24k core-hours | ~64 GB, none |

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

`exportStepRecords` (default **True**) governs the RAW per-step export that the
exponents are built from:

* dropped when False: `ioniurbanidx`/`ioniurbanv`, `msmoliidx`/`msmoliv`,
  `radstepidx`/`radstepv`/`radstepspecv`, `reseigv`, `resinfv`, `resinfbv`;
* **kept regardless**: `reseigidx`, `reshitidx`, `resinfvarv`, `resinfcov`,
  `ioniqscaleidx`/`ioniqscalev`, `radvgrid`, `radstepstride`, `radstepnv`.

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

The requirement is 1e-6 absolute on the exponent; the measured agreement is
~1e-11 (see `Documents/Resolution/NOTES.md`, 2026-09-05).
