# Global material model — design and implementation plan

Status: PLANNED (design v1, 2026-07-16). Owner: David Walter.
The material-sector analogue of the scalar-potential B-field reduction
(15k per-module dBz -> ~50 parmtype-14 modes); replaces ~25k per-module
energy-loss parameters with O(10-100) parameters attached to physical
material volumes, making the material correction independent of a track's
origin and direction.

## 1. Motivation

The legacy material correction assigns one free parameter per silicon
module, k_xi, scaling the energy loss of the *leg ending at that module*:
dE = e^{k_xi} (dE/dx)_0 ds. Because every prompt track from the nominal IP
traverses the same material sequence on the way to a given module, the
per-module parameter absorbs the *integrated path* material. This is
correct only for prompt tracks:

- **Displaced-vertex daughters** (K_S, Lambda, B -> J/psi K) originate
  beyond the beam pipe and inner services; the inner material folded into
  the per-module parameters was never crossed.
- **Cosmics** cross the detector top-down, in the opposite radial order,
  at incidence angles the prompt parameterization never sees.

These are exactly the channels needed to break the A–epsilon (field vs
material) degeneracy and to extend the calibration program. A model whose
parameters live on the *material itself* — evaluated along each track's
actual path — is origin- and direction-independent by construction.

Second driver: tractability. ~25k eloss + ~15k field parameters make the
unbinned calibration fit expensive; ~50 field modes (done, parmtype 14)
plus O(30-60) material groups make the combined fit small enough for
GPU-resident unbinned likelihoods.

## 2. Physics model

Assign each material *group* g (a set of Geant4 volumes, active or
passive) a scale s_g = e^{k_g} multiplying the material amount a step
sees. Staged scope:

- **M0 (baseline, this plan):** scale the mean ionization only —
  dE|_step = e^{k_g(step)} (dE/dx)_0 ds. Semantically identical to the
  legacy xi, but resolved per step and attached to the volume crossed
  rather than the leg destination. Uses the existing SetMaterialOffset /
  dxi machinery in the propagator.
- **M1 (extension):** coherent density scaling — the same k_g also scales
  the energy-loss fluctuation variance and the multiple-scattering
  strength (1/X0) in the Q matrix. Physically more correct (a missing
  cable is missing for MS too), but it changes fit *weights*, not just
  means. DECIDED: implemented as a two-step scheme — fitted M0 k_g values
  are propagated into the Q matrix at the next refit iteration, never
  fitted through the weights; adoption gated on a resolution-closure
  comparison before/after M0.

The correction enters the propagation exactly where the per-module xi
enters today; only the lookup changes from "per leg" to "per step, by
volume group".

## 3. Volume -> group mapping

- Source of truth: the Geant4 geometry (DDDWorld, built once by
  CvhMasterThread). During propagation each G4 step knows its physical
  volume/touchable; map G4VPhysicalVolume* -> groupId with a hash table
  built at master initialization (const afterwards, shared by workers —
  same pattern as the master field).
- Grouping defined by logical-volume name patterns (DDD names) + optional
  z-side/layer qualifiers, stored in a versioned text file
  (analogous to the scalar-potential init file):
  `pattern  zside  groupId  k_init  prior_sigma`.
- **Two grouping tiers** (DECIDED 2026-07-16), mirroring the B-field
  50/360-mode dual model — both defined from the same Phase 0 audit and
  interchangeable via the init file:
  - `materialGroups50.txt` — O(50) groups with symmetry assumptions:
    z-symmetric barrel supports and services, phi-symmetric everywhere,
    TID/TEC rings merged per disk; sensors per layer (own groups, see
    section 6).
  - `materialGroups100.txt` — O(100) groups relaxing the symmetries:
    z+/z- split services and endcap elements, finer ring/disk splits,
    barrel support half-shells where the audit shows material contrast.
  The tier comparison itself is a validation handle: fitted physics
  (J/psi closure) must agree between tiers; where the O(100) tier pulls
  a symmetry-split pair apart significantly, the O(50) symmetry
  assumption is measured to be wrong.
- Data-driven split/merge refinement within each tier later (split groups
  whose fitted k_g pulls differ by channel/eta; merge pairs with fit
  correlation > 0.95).
- **Phase 0 audit tool** defines the groups from numbers, not guesswork:
  a geantino-ray tally over (eta, phi) producing per-volume integral
  dE/dx and 1/X0 tables — both the grouping input and a standalone
  material-budget validation deliverable.

## 4. Propagator changes (Geant4ePropagator)

Current state: `propagateGenericWithJacobianAltD` composes per-step 5x9
transport Jacobians (5 state + dBx,dBy,dBz + dxi); the dxi column is
accumulated over the whole leg (`jac.rightCols<4>() += ...`), i.e. one
integrated material-scaling derivative per leg. `SetMaterialOffset(dxi)`
applies one xi value for the whole leg.

Changes:
1. Per step, look up groupId; maintain `map<groupId, Vector5d>` per leg:
   the step's dxi column transported to the leg end with the same
   composition as today, summed into its group's slot. A leg crosses
   ~1-5 groups, so this is a handful of 5-vectors.
2. The existing single dxi column equals the sum over groups — kept, and
   used as an exact validation identity (V2 below).
3. Applying the correction: SetMaterialOffset becomes per-step — the
   wrapper takes the groupId-resolved k_g (values held by the maker,
   passed per leg as today's dxi is).
4. Backward legs: each group column gets the same P-conjugation as the
   existing rightCols(4) block; dE sign flip unchanged.
5. Return: per-leg vector of (groupId, column, dE_g, ds_g) alongside the
   legacy outputs.

## 5. Maker integration (ResidualGlobalCorrectionMaker*)

- New parameter type `ParmTypeMaterialGlobal = 15`, sentinel
  DetId(groupIdx) — mirror the parmtype-14 registration, runtree
  bookkeeping, corFiles update path, and init-file seeding.
- Gradient assembly: chain rule leg-by-leg exactly as for the field
  modes; per-track global columns grow by the number of groups crossed
  (~10-30) while the ~17 per-module dxi columns per track are removed —
  net reduction in stored gradients.
- Config: `globalMaterialModel` flag + `materialGroupsFile`; when on, the
  per-module dxi parameters are not registered (exclusive by default; a
  both-on mode only for dedicated degeneracy studies).
- Factored Hessian (fillGradsFactored) unaffected structurally.
- Scope order: two-track maker first (calibration), single-track next,
  three-track (B -> J/psi K) when that lands.

## 6. Degeneracies and constraints

- Prompt-only data constrain only the *cumulative* material along
  IP-to-module rays — group combinations along a ray are degenerate
  (this is the legacy situation made explicit). Broken by: displaced
  daughters (start mid-sequence), cosmics (reverse order, off-axis
  paths), and species (kaon vs muon dE/dx at same pT reweights groups
  differently and separates material from field).
- Priors: k_g ~ N(0, sigma_g) with sigma_g from the material-budget
  uncertainty of that group (few %% for well-known elements like the beam
  pipe, looser for services). Keeps under-constrained groups bounded
  until the displaced/cosmic channels weigh in.
- The overall-scale combination correlates with the epsilon term of the
  residual curvature model — same bookkeeping as the field A term;
  document in the combined-fit spec.

## 7. Validation plan (each gate blocks the next phase)

- **V0** all k_g = 0: bit-identical to current code (smoke + 10k, both
  makers) — the accounting must be a pure spectator.
- **V1** FD closure per group: inject k_g = eps in propagation, compare
  Delta(state) against J_g * eps per leg, forward and backward legs
  (mirrors the Phase B.5 field closure).
- **V2** sum identity: Sigma_g J_g == legacy integrated dxi column,
  exact per leg.
- **V3** MC closure: refit with injected k_g != 0 applied in propagation;
  the fit must recover -k_g per group on J/psi MC.
- **V4** prompt data: fitted k_g vs legacy per-module solution projected
  onto groups; J/psi mass closure vs (pT, eta) at least as good as
  legacy.
- **V5** origin/direction independence — the point of the model: the
  same k_g set must close J/psi (prompt), K_S/Lambda (displaced,
  runGlobalCorRecKs*/Lambda* drivers exist in 10_6_dev and need porting),
  and cosmics (new driver + ALCARECO). Per-channel pulls of k_g are the
  headline plot. **POSTPONED (2026-07-16): the K_S/Lambda/cosmics samples
  are not yet produced — until they are, the program runs on J/psi
  alone.** Consequences while J/psi-only: (a) only along-ray cumulative
  combinations of k_g are constrained; priors carry the decomposition,
  so V4 must report the fitted-vs-prior pull per group honestly (groups
  at their prior are unmeasured, not confirmed); (b) J/psi-only
  substitutes for V5 in a weaker form: closure vs (pT, eta, phi) and the
  50-vs-100 tier comparison; (c) the model can and should still be
  *implemented* now — parameters attached to volumes lose nothing by
  being fitted first with prompt data, and the displaced/cosmic channels
  drop in later without any code change.
- **V6** combined field+material fit: correlation matrix, A–epsilon
  separation with B -> J/psi K kaons.

## 8. Performance and storage

- Per-step pointer-hash lookup: O(1), negligible against field evaluation
  (31%% of producer time).
- Track gradient payload: +O(10-30) global columns, -~17 per-module dxi
  columns -> net smaller NanoAOD.
- Fit side: material block shrinks ~25k -> O(50); combined with the
  50-mode field block the global part of the unbinned fit becomes
  trivially GPU-resident.

## 9. Phasing

- **Phase 0 — audit & grouping (standalone):** geantino tally tool,
  material-budget tables, materialGroups.txt v1. No fit code touched.
- **Phase A — propagator:** per-group columns + per-step offset
  application; gates V0/V1/V2.
- **Phase B — maker:** parmtype-15 registration, gradient storage,
  ALCARECO smoke; V0 repeated end-to-end.
- **Phase C — fit side:** calibration_studies integration, V3/V4.
- **Phase D — channels (DEFERRED until samples exist):** K_S/Lambda
  driver port, cosmics driver + samples; V5. Not on the critical path —
  Phases 0-C complete on J/psi alone.
- **Phase E — combined calibration:** joint field+material fit, V6;
  document in AN-XX-XXX.

Development in the 15_0 port (CMSSW_15_0_19_patch2_dev branch); PR after
Phase B, mirroring the scalar-potential development flow.

## 10. Decisions (2026-07-16) and remaining open points

Decided:
- **M0 first; M1 as two-step.** k_g scales the mean dE/dx only in the
  fit; after an M0 fit converges, fitted k_g values are propagated into
  the Q matrix (MS + fluctuations) for the next refit iteration rather
  than fitted through the weights. Adoption of M1 gated on a resolution
  closure comparison before/after M0.
- **Two grouping tiers** O(50) (symmetry assumptions) / O(100) (relaxed),
  see section 3.
- **Sensors in their own groups**, separate from same-layer supports, so
  the sensor-vs-service degeneracy shows up as a fit correlation rather
  than being hidden (section 6).
- **Priors from engineering budget uncertainties** as plain Gaussians;
  photon-conversion / nuclear-interaction maps used as a cross-check
  overlay on the Phase 0 audit, not imported as fit constraints.
- **One parameter set per data-taking era**; stability checked by fitting
  in run-blocks as a validation, not modeled.
- **Exclusive switch** vs legacy per-module dxi (globalMaterialModel
  flag). Both-registered is supported only in the frozen-legacy mode for
  the V4 comparison (registered but not floating, via the existing
  corFiles freeze notion).
- **J/psi-only for now**: displaced/cosmic validation (V5, Phase D)
  deferred until those samples are produced.

Open, with owner/trigger:
- **Application contract (linear analysis-level correction vs re-seeded
  refit; where the linearity of the stored Jacobians breaks): decision
  deferred to Phase B/C design time.** To be brought back to David with
  a concrete recommendation and the V3 numbers (injected-k_g size at
  which the linear application diverges from the refit) in hand.
- Geant4/DDD description as the accuracy ceiling: shape errors (missing
  or misplaced elements) are not representable by any k_g; diagnosed by
  localized, channel-dependent pulls that splitting does not cure, and
  by the Phase 0 audit vs conversion maps. Budget one round of
  geometry-description fixes.
- Split/merge refinement loop cadence within each tier (data-driven,
  starts after V4).
