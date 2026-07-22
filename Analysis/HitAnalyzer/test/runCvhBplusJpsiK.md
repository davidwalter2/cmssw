# `runCvhBplusJpsiK.py` — Stage-2 CVH refit driver for B⁺ → J/ψ K⁺

This document explains what the `cmsRun` driver
[`runCvhBplusJpsiK.py`](runCvhBplusJpsiK.py) does, with the bulk devoted to the
**CVH refit** itself: how the Gauss–Newton iteration works, which knobs control
it, how to ask for *more steps*, and where the hard-coded **maximum propagation
length before the track is declared to "miss every module"** lives.

---

## 1. What the driver schedules

The driver runs one `cms.Process('CVHBPLUS')` that chains three modules on a
single path:

```
offlineBeamSpot  ->  jpsiKCandidateSplitter  ->  <one CVH maker>
```

1. **`offlineBeamSpot`** — standard `BeamSpotProducer`.
2. **`jpsiKCandidateSplitter`**
   ([`JpsiKCandidateSplitter_cfi.py`](../python/JpsiKCandidateSplitter_cfi.py))
   — takes the upstream B⁺ `VertexCompositeCandidate`s (daughter 0 = J/ψ VCC,
   daughter 1 = bachelor-kaon `RecoChargedCandidate`) and splits them into a
   `dimuon` collection, a `bachelor` kaon track collection, and the parallel
   `…BCandIdx` index vectors used to join the two refit rows back together
   offline on `(run, lumi, event, bCandIdx)`.
3. **CVH maker(s)**, chosen by `mode=`:
   - `mode=both`   → **default.** Both makers in ONE `cmsRun` job.
   - `mode=dimuon` → `globalCorJpsiK` (two-track maker, mass-constraint ON)
   - `mode=kaon`   → `globalCorJpsiKKaon` (single-track bachelor maker)

`mode=both` used to crash: each maker owned its own `CvhMasterThread` (a
per-module-label `edm::GlobalCache`), so the second `G4MTRunManagerKernel`
tripped Geant4's single-master singleton. **That is fixed** — `e33c8d`
made the master an EventSetup product on `CvhMasterRecord`, built once per
job by `cvhMasterESProducer` and `esConsumes`-ed by every maker. Verified:
both fit summaries print from a single process.

Because both makers now run in one job, the per-candidate results are paired
in memory and **no offline join is needed**.

The two-track maker reads the **nested** stage-1 candidate directly
(`srcCandidates=ALCARECOTkAlJpsiXBPlusResonances`, the default): it descends
composite daughters itself, so the splitter is not needed on the dimuon side.
Verified bit-identical to the legacy pre-split input — both give
`attempted=15 succeeded=15`, pixel `seen=110 sizeX1=29 demoted=29` on the same
5 events. The splitter is still scheduled only to feed the single-track kaon
maker its `bachelor` track collection.

This document focuses on the **dimuon two-track maker**
([`ResidualGlobalCorrectionMakerTwoTrackG4e.cc`](../plugins/ResidualGlobalCorrectionMakerTwoTrackG4e.cc)),
since that is where the mass-constrained CVH vertex fit lives.

### Typical invocations

```bash
# Both makers, one job (default)
cmsRun runCvhBplusJpsiK.py input=/path/to/alcareco_sl1.root nEvents=200

# Soft bachelor kaons need a lowered propagation floor: with the default
# plimit=1 GeV most bachelors abort with Geant4e fail[plimit].
cmsRun runCvhBplusJpsiK.py input=/path/to/alcareco_sl1.root nEvents=200 plimit=0.05

# Single-maker A/B runs
cmsRun runCvhBplusJpsiK.py input=… nEvents=200 mode=dimuon
cmsRun runCvhBplusJpsiK.py input=… nEvents=200 mode=kaon
```

---

## 2. The CVH refit, conceptually

CVH ("Continuous Vertex/Hit", the global-correction G4e refit) re-fits each
track by propagating it through the **real CMS material and magnetic field with
Geant4e**, layer by layer, and solving a sparse Gauss–Newton (GBL) linear system
that simultaneously corrects the trajectory at every layer. For the dimuon side
of B⁺ → J/ψ K⁺ the two muon tracks are fit **jointly** with a common vertex and,
on the constrained pass, a **J/ψ mass constraint**.

There are **two nested loops** you need to understand before touching any knob:

### 2.1 Outer loop — constraint phases (`icons`)

```cpp
const unsigned int nicons = doMassConstraint_ ? 2 : 1;     // line 1381
for (unsigned int icons = 0; icons < nicons; ++icons) { … } // line 1396
```

With `doMassConstraint=True` (the default for `globalCorJpsiK`,
[`…TwoTrackJpsiKMuMuG4e_cfi.py:59`](../python/ResidualGlobalCorrectionMakerTwoTrackJpsiKMuMuG4e_cfi.py)),
the maker runs **two phases**:

- **`icons=0`** — unconstrained common-vertex fit (`KinematicParticleVertexFitter`).
- **`icons=1`** — adds the J/ψ mass constraint
  (`TwoTrackMassKinematicConstraint`, `KinematicConstrainedVertexFitter`),
  pulling m(μμ) to `massConstraint=3.0969 GeV` with width
  `massConstraintWidth=9.29e-5 GeV`.

Each phase seeds the next; the mass row is only added on `icons==1`
([line 1410](../plugins/ResidualGlobalCorrectionMakerTwoTrackG4e.cc#L1410)).

### 2.2 Inner loop — Gauss–Newton iterations (`iiter`)

This is the loop you most likely mean by "number of steps":

```cpp
// line 1518
const unsigned int niters = (dogen && !dolocalupdate) ? 1 : nIters_;
for (unsigned int iiter = 0; iiter < niters; ++iiter) { … }   // line 1524
```

Each iteration:

1. Propagates both tracks through every layer with Geant4e (this is where the
   propagation-length limit of §5 bites).
2. Assembles the sparse GBL system (5 propagation/MS rows per hit + 2
   measurement rows per valid hit + constraint rows).
3. Solves for the parameter update `dxfull` and the χ² change.
4. Computes the **EDM** (estimated distance to minimum):

   ```cpp
   edmval = -deltachisq;   // line 3175
   ```

5. Checks convergence and breaks early if converged (§4).

So the *total* number of propagation/solve steps per candidate is
roughly `nicons × niters` minus whatever early-convergence saves.

---

## 3. Increasing the number of steps (`nIters`)

`nIters` is the **per-phase Gauss–Newton iteration cap**. It is a registered
VarParsing knob, so you set it straight on the command line:

```bash
cmsRun runCvhBplusJpsiK.py input=… mode=dimuon nIters=50
```

- Registered at [`runCvhBplusJpsiK.py:78`](runCvhBplusJpsiK.py#L78), default **10**.
- Passed to the **dimuon** maker as `nIters=cms.uint32(int(opts.nIters))`
  ([line 157](runCvhBplusJpsiK.py#L157)). Note the driver does **not** forward
  `nIters` / `edmConvergence` / `useStartingState` / `debug` to the kaon-side
  maker clone ([lines 177–187](runCvhBplusJpsiK.py#L177)), so the single-track
  bachelor refit always runs at the built-in defaults (10, 1e-5, `perigee`).
  If you want to sweep iterations on the kaon side too, add those clones to the
  `globalCorJpsiKKaon` block.
- Read in the maker with an `existsAs` guard so legacy cfis still default to 10:
  [`…TwoTrackG4e.cc:345`](../plugins/ResidualGlobalCorrectionMakerTwoTrackG4e.cc#L345).
- Consumed as the loop bound at
  [line 1518](../plugins/ResidualGlobalCorrectionMakerTwoTrackG4e.cc#L1518).

The default of **10 reproduces the published Run2016H baseline bit-identically**.
The openspec change that introduced this knob sweeps the matrix `{10, 20, 50,
100}`. Note the cap is **per `icons` phase**, so `nIters=50` with the mass
constraint on means up to `2 × 50 = 100` iterations per candidate (early
convergence usually stops it well short).

> If you raise `nIters` and *also* want to see whether the extra iterations are
> actually doing anything, turn on the per-iter debug dump with `debug=True`
> (§6) and inspect `edmval_iter` / `Jpsi_mass_iter` per candidate.

---

## 4. Convergence threshold (`edmConvergence`)

`nIters` is only the *ceiling*. The loop normally exits earlier when the EDM
drops below `edmConvergence`:

```cpp
// lines 3266–3271
if (iiter > 0 && dolocalupdate && edmval    < edmConvergence_)  break;
else if (iiter > 0 && !dolocalupdate && edmvalref < edmConvergence_) break;
```

- Registered at [`runCvhBplusJpsiK.py:82`](runCvhBplusJpsiK.py#L82), default
  **1e-5** (matches the baseline). Matrix sweeps `{1e-5, 1e-3, 1e-2}`.
- Read at [`…TwoTrackG4e.cc:347`](../plugins/ResidualGlobalCorrectionMakerTwoTrackG4e.cc#L347).

Interaction with `nIters`:

| Goal | What to change |
|------|----------------|
| Let the fit run longer / settle harder | raise `nIters`, keep `edmConvergence` tight (1e-5) |
| Stop sooner once "good enough" | loosen `edmConvergence` (1e-3) |
| Force a fixed iteration count (no early stop) | set `edmConvergence` to something tiny like `1e-30` and rely on `nIters` |

`useStartingState` (default `perigee`) selects the iteration-0 reference state.
`midPropagated` is wired but **throws `cms::Exception`** until the propagation
helper lands ([`…TwoTrackG4e.cc:357`](../plugins/ResidualGlobalCorrectionMakerTwoTrackG4e.cc#L357)),
so leave it on `perigee`.

---

## 5. Maximum propagation length before a track "misses every module"

This is the limit you asked about. It is **not** (currently) a VarParsing knob —
it is hard-coded inside the Geant4e propagator. Understanding it requires
knowing how the per-layer propagation actually terminates.

### 5.1 How a single propagation step terminates

For each layer, the CVH maker calls
`Geant4ePropagator::propagateGenericWithJacobianAltD`
([`Geant4ePropagator.cc:609`](../../../TrackPropagation/Geant4e/src/Geant4ePropagator.cc#L609)).
Inside, Geant4e steps the track toward the target surface in a `while` loop,
accumulating `finalPathLength` (in **cm** —
`g4doubleToCmsDouble` divides the G4 mm value by `cm`,
[`ConvertFromToCLHEP.h:46`](../../../TrackPropagation/Geant4e/interface/ConvertFromToCLHEP.h#L46)):

```cpp
finalPathLength += thisPathLength;                 // line 894

if (std::fabs(finalPathLength) > 10000.0f) {       // line 896  <-- THE LIMIT
  // reached maximum path length, bail out
  ++propFailCounts_[2];
  std::cout << "Geant4e fail[maxlen] …";
  return retDefault();                             // invalid state
}

if (theG4eManager->GetPropagator()->CheckIfLastStep(...)) {  // hit the target
  continuePropagation = false;
}
```

So the loop ends in exactly one of two ways:

- **Success:** `CheckIfLastStep` fires — the track reached the destination
  surface (the next module/layer).
- **`fail[maxlen]`:** the track has travelled **more than `10000.0` cm
  (= 100 m)** without ever intersecting the target surface. This is the
  *"maximum length allowed without hitting any module"*. Beyond it the
  propagator gives up and returns an invalid `TrajectoryStateOnSurface`, which
  the maker treats as a failed hit.

### 5.2 Where it is, and how to change it

The 100 m cap appears **twice** in the propagator:

| Line | Function | Used by |
|------|----------|---------|
| [`Geant4ePropagator.cc:526`](../../../TrackPropagation/Geant4e/src/Geant4ePropagator.cc#L526) | `propagateGeneric` (legacy, `!forCVH_` path) | non-CVH reco |
| [`Geant4ePropagator.cc:896`](../../../TrackPropagation/Geant4e/src/Geant4ePropagator.cc#L896) | `propagateGenericWithJacobianAltD` | **the CVH refit** |

For the CVH refit, **edit line 896** (the `AltD` variant). To allow longer
propagations (e.g. very low-pT helices that spiral a long way before reaching a
layer), raise the `10000.0f`:

```cpp
if (std::fabs(finalPathLength) > 50000.0f) {   // 500 m, e.g.
```

There is **no cfi/VarParsing override** for this value today — it is a literal,
so changing it requires recompiling `TrackPropagation/Geant4e`
(`scram b`). If you expect to sweep it, the clean change is to promote it to a
member read from a PSet (mirroring how `plimit_` is already plumbed, §5.4) and
wire it through `runCvhBplusJpsiK.py`; until then it is a source edit.

> **Symptom that you are hitting this limit:** lines like
> `Geant4e fail[maxlen]  iter=…  pathLen=…  pT=…  eta=…` on stdout, and a
> non-zero `exit3[maxlen]` in the destructor summary
> ([`Geant4ePropagator.cc:140`](../../../TrackPropagation/Geant4e/src/Geant4ePropagator.cc#L140)):
> ```
> Geant4ePropagator::propagateGenericWithJacobianAltD summary
>    exit1[plimit]=…   exit2[ierr]=…   exit3[maxlen]=…
> ```

### 5.3 Per-step length limit (distinct from the total cap)

Do not confuse the **total-path cap above** with the **per-step length limit**,
which is the maximum length of a *single* Geant4e step:

```cpp
G4UImanager::GetUIpointer()->ApplyCommand("/geant4e/limits/stepLength 10.0 mm");
```

Set on the master in
[`CvhMaster.cc:123`](../../../TrackPropagation/Geant4e/src/CvhMaster.cc#L123)
and re-applied per worker in
[`CvhWorker.cc:110`](../../../TrackPropagation/Geant4e/src/CvhWorker.cc#L110).
Smaller per-step length = more steps = finer material/field sampling but slower;
this is the granularity knob, whereas `10000.0f` is the give-up distance.

### 5.4 Low-momentum cut (`plimit`) — the *other* propagation exit

The third propagation failure mode is the momentum floor:

```cpp
if (!configurePropagation(mode, pDest, cmsInitPos, cmsInitMom)) {
  ++propFailCounts_[0];                 // line 712
  std::cout << "Geant4e fail[plimit] …";
  return retDefault();
}
```

This is governed by `Geant4ePropagator.PropagationPtotLimit`, which **is** a
VarParsing knob (`plimit`, default 1.0 GeV/c):

```bash
cmsRun runCvhBplusJpsiK.py input=… mode=kaon plimit=0.05
```

- Registered at [`runCvhBplusJpsiK.py:70`](runCvhBplusJpsiK.py#L70).
- Applied at [`runCvhBplusJpsiK.py:218`](runCvhBplusJpsiK.py#L218)
  (`process.Geant4ePropagator.PropagationPtotLimit`).
- Lowering it to 0.05 recovers the **soft-bachelor-kaon tail** (the §9.8 study).
  On the dimuon side it is a no-op — muons clear 1.0 GeV/c trivially.

**Summary of the three propagation exits** (`propFailCounts_[0..2]`):

| Counter | Tag | Cause | Knob |
|---------|-----|-------|------|
| `[0]` | `plimit` | momentum below floor | `plimit` (VarParsing) |
| `[1]` | `ierr`   | G4 `PropagateOneStep` returned an error | — |
| `[2]` | `maxlen` | travelled >100 m without hitting the target | hard-coded `10000.0f`, line 896 |

---

## 6. Diagnostics: per-iteration debug dump (`debug=True`)

To see the convergence trajectory of every candidate, enable the per-iter dump:

```bash
cmsRun runCvhBplusJpsiK.py input=… mode=dimuon nIters=50 debug=True
```

- Registered at [`runCvhBplusJpsiK.py:91`](runCvhBplusJpsiK.py#L91), passed as
  `debugPerIterDump` ([line 160](runCvhBplusJpsiK.py#L160)).
- Adds vector branches to the dimuon tree, filled at
  [`…TwoTrackG4e.cc:3202`](../plugins/ResidualGlobalCorrectionMakerTwoTrackG4e.cc#L3202):
  `chisqval_iter`, `edmval_iter`, `deltachisqval_iter`, `mu_qoverp_iter`,
  `Jpsi_mass_iter`. The vectors concatenate both `icons` phases in order.
- Costs ~80 B/event — use it only for the matrix deep-dive, not production.

Scalar convergence summaries are always present: `edmval_cons0` and
`niter_cons0` record the EDM and iteration count of the **unconstrained**
phase ([lines 3192–3195](../plugins/ResidualGlobalCorrectionMakerTwoTrackG4e.cc#L3192)).

---

## 7. Magnetic field and the FD closure

- **Field** — by default (`useScalarPot3D=True`) the driver loads the
  scalar-potential 3D field producer (`ScalarPot3DMf` label) and points both the
  `Geant4ePropagator` and the makers at it
  ([`runCvhBplusJpsiK.py:202`](runCvhBplusJpsiK.py#L202)). `useScalarPot3D=False`
  falls back to the stock CMSSW field for the §9.4.b A/B test.
- **FD closure** — `runFDClosure=True` runs a one-shot finite-difference
  Jacobian cross-check (`epsilonFDClosure`, default 1e-4), printing
  `===== Numerical-FD closure =====`. Use it to validate analytic Jacobians, not
  in bulk processing.
- **Mass-hypothesis A/B test** — `kaonAsMuon=True` swaps the *propagation*
  hypothesis on the kaon maker from `kaon` to `mu` while keeping the bachelor
  kaon tracks as input ([§9.4.c, line 184](runCvhBplusJpsiK.py#L184)).

---

## 8. Quick reference — knobs you set on the command line

| Knob | Default | Effect | Source |
|------|---------|--------|--------|
| `input` | — (required) | ALCARECO sl1 file | [L36](runCvhBplusJpsiK.py#L36) |
| `nEvents` | 200 | events to process (`-1` = all) | [L39](runCvhBplusJpsiK.py#L39) |
| `mode` | `dimuon` | which maker (`dimuon`/`kaon`; `both` crashes) | [L48](runCvhBplusJpsiK.py#L48) |
| **`nIters`** | **10** | **Gauss–Newton iteration cap per phase — raise for more steps** | [L78](runCvhBplusJpsiK.py#L78) |
| `edmConvergence` | 1e-5 | early-stop EDM threshold | [L82](runCvhBplusJpsiK.py#L82) |
| `useStartingState` | `perigee` | iter-0 reference (`midPropagated` throws) | [L86](runCvhBplusJpsiK.py#L86) |
| `plimit` | 1.0 | low-momentum propagation floor [GeV/c] | [L70](runCvhBplusJpsiK.py#L70) |
| `debug` | False | per-iter convergence-trace branches | [L91](runCvhBplusJpsiK.py#L91) |
| `useScalarPot3D` | True | scalar-pot 3D field vs stock CMSSW field | [L59](runCvhBplusJpsiK.py#L59) |
| `kaonAsMuon` | False | kaon-side mass-hypothesis A/B | [L65](runCvhBplusJpsiK.py#L65) |
| `runFDClosure` / `epsilonFDClosure` | False / 1e-4 | FD Jacobian closure test | [L52](runCvhBplusJpsiK.py#L52) |
| `fillJac` | True | store per-track Jacobians | [L41](runCvhBplusJpsiK.py#L41) |

**Not a command-line knob** (requires a source edit + recompile):

| Quantity | Value | Location |
|----------|-------|----------|
| Max propagation length before `fail[maxlen]` (CVH path) | `10000.0f` cm (100 m) | [`Geant4ePropagator.cc:896`](../../../TrackPropagation/Geant4e/src/Geant4ePropagator.cc#L896) |
| Same, legacy non-CVH path | `10000.0f` cm | [`Geant4ePropagator.cc:526`](../../../TrackPropagation/Geant4e/src/Geant4ePropagator.cc#L526) |
| Per-step length limit | `10.0 mm` | [`CvhMaster.cc:123`](../../../TrackPropagation/Geant4e/src/CvhMaster.cc#L123), [`CvhWorker.cc:110`](../../../TrackPropagation/Geant4e/src/CvhWorker.cc#L110) |
