## Single-track CVH gen-closure driver with resolution parameters (doRes).
##
## Re-runs the 2022 resolution-correction closure test with the current
## (2026) hit selection and fit infrastructure: refits gen-matched muon
## tracks of the inclusive B->J/psi+X MC ALCARECO (2016 postVFP, split=1,
## gen kept) with fitFromGenParms=True (reference state frozen to gen ->
## no weak modes, no mass constraint needed) and doRes=True, which
## registers the log-variance resolution parameter families
##   parmtype  8: local-x/phi hit resolution scale (per module)
##   parmtype  9: local-y hit resolution scale (per pixel module)
##   parmtype 10: multiple-scattering process-noise scale (per glued detid)
##   parmtype 11: ionization fluctuation scale (per glued detid)
## and emits their gradient/Hessian rows (incl. the log-det terms) into
## the fillGrads output. The closure criterion: solving the global system
## with any resolution family floated (calibration_studies/
## global_corrections/fit_global_grads.py --parmtypes 10) must return
## parameters consistent with zero. The 2022 attempt failed this because
## non-Gaussian straggling/MS tails pulled the second-moment estimator;
## see Documents/Resolution/NOTES.md.
##
## Assembled from runCvhSingleTrack.py (single-track maker block) and
## runCvhJpsiGenMC.py (MC conditions, source handling, gen anchoring).
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os

from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
from Configuration.AlCa.GlobalTag import GlobalTag

opts = VarParsing.VarParsing('analysis')
opts.register('tightG4eStepper', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'tighten the Geant4e field-integration tolerances in the FIT to '
              'match the simulation (see the block after geantRefit_cff). '
              'False reproduces every result before 2026-08-08.')
opts.register('input', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'comma-separated absolute paths or root:// URLs of ALCARECO files')
opts.register('inputFileList', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'text file with one input path per line (the MC has ~10 events '
              'per file, so runs typically need many files); combined with '
              'input= if both are given')
opts.register('nEvents', 500, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, 'number of events to process (-1 = all)')
opts.register('particle', 'mu', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'track hypothesis + gen-match species: mu (default), pi, kaon, '
              'proton. Sets trackParticleName, the gen-match |pdgId| and the '
              'CvhMaster G4 particle list. The TkAlJpsiX MC track collection '
              'contains the non-muon B daughters, so pi/kaon closures run on '
              'the same input files.')
opts.register('doRes', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'register the resolution parameter families (parmtypes 8-11) '
              'and emit their gradient/Hessian rows (default True; set False '
              'for a doRes-off reference run with identical selection)')
opts.register('exportStepRecords', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'write the RAW per-step resolution export (ioniurbanv, msmoliv, '
              'radstepv/radstepspecv, reseigv, resinfv, resinfbv). It is 430 kB '
              'per candidate -- 16 TB over the 40M candidates of the full '
              'calibration -- and its only consumer was the offline exponent '
              'extractor, which the in-maker cfqop_*/cfmass_* export replaces. '
              'Default FALSE since 2026-09-06 -- every production since '
              '2026-09-05 set it so explicitly and the exponents are '
              'validated; set True to re-derive them under a different model')
opts.register('materialFDGroup', -1, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'one-shot finite-difference closure of one material group: the '
              "group's mean-loss column (V1) and, with doRes on, its PROCESS "
              'NOISE block (V3, the parmtype-15 dV this maker registers). -1 = off')
opts.register('materialFDEps', 1e-3, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'the injected delta k_g of the material FD closure')
opts.register('exportMaterialNoise', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'register the parmtype-15 MATERIAL-GROUP process noise as a '
              'resolution family, so the quadratic term differentiates a '
              "group's WIDTH (its MS covariance and ionization variance, both "
              'scaled by exp(k_g)) as well as its mean loss. It CHANGES the '
              'exported gradient and Hessian of the parmtype-15 columns -- not '
              'the track fit -- so it is off by default')
opts.register('exportVarianceGrads', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'add the VARIANCE (log-det) part of the profiled -2lnL to the '
              'exported global gradient and Hessian: '
              '-r^T V^-1 dV V^-1 r + tr(V^-1 dV) and the expected curvature '
              'tr(dV R dV R). Without it a parameter that moves the '
              'covariance (parmtype 15 through exp(k_g), and 8/9/10/11 '
              'entirely) enters the quadratic hit-chi2 term only through the '
              'MEAN. Two-track maker only (the single-track one always had '
              'it); OFF reproduces the pre-2026-09-06 gradients bit for bit')
opts.register('varianceGradFamilies', [], VarParsing.VarParsing.multiplicity.list,
              VarParsing.VarParsing.varType.int,
              'which parmtypes exportVarianceGrads covers; empty = '
              '{8,9,10,11,15}. 15 alone is the LAYOUT-PRESERVING subset (the '
              'material-group globals are already columns of globalidxv), so '
              'its output can be pooled with a production that ran without '
              'the switch; 8/9/10/11 append per-module columns')
opts.register('exportObjective', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'write objval/objchisq/objlogdetv/objlogdetc, the marginal '
              'objective r^T R r + ln|V| + ln|C| in double precision. '
              'Validation only -- it costs an ncons x ncons LDLT per '
              'candidate -- and exists so the new gradient can be '
              'finite-differenced against what it claims to differentiate')
opts.register('exportHitResBlocks', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'register the parmtype-8/9 HIT-RESOLUTION dV blocks in the '
              'influence export (reseigidx/resinfvarv/reshitcls + the '
              'cf*_hitcls/cf*_hitv per-class shares). The single-track maker '
              'has always done it; the two-track one did not, which is why the '
              'per-hit-class resolution parameters have never been fitted. '
              'Export only -- it cannot move the fit. Set False to reproduce a '
              'pre-2026-09-06 two-track tree')
opts.register('exportCfGroupExponents', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'additionally split the CF exponents by parmtype-15 MATERIAL '
              'GROUP (cf*_grp, cf*_grp_ms, ...). This is what lets the fit '
              'float the material amount per group instead of four per-family '
              'k knobs. ~27 kB/candidate against 1.4 kB for the flat '
              'exponents, so it is off unless the output feeds the joint '
              'material+field fit')
opts.register('exportCfExponents', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'compute the resolution-CF exponents in the maker and write them '
              'on the 64-point tau grid (cfqop_* single-track, cfmass_* '
              'two-track). Default True')
opts.register('propagationPtotLimit', 0.2, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'G4e propagation momentum floor [GeV]; cfi default was 1.0')
opts.register('maxMomentumStepFactor', 2.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'RELATIVE Gauss-Newton step damping: the max factor by which a '
              'track momentum may change in one iteration (default 2 = p may '
              'at most halve or double). Implemented as the effective floor '
              'max(clampMomentumFloor, p_ref/f) and the symmetric cap p_ref*f, '
              'so the bound is ALWAYS strictly inside p_ref -- unlike the bare '
              'absolute floor it can neither pin a genuinely soft track at a '
              'fixed momentum nor scale the step to exactly zero. Set <=1 to '
              'switch it off and get the legacy absolute-floor-only clamp '
              '(bit-identical).')
opts.register('stepBacktracking', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'chi2-based (Armijo) retroactive step backtracking. The chi2 '
              'assembled in iteration k is the realized chi2 of the step taken '
              'at k-1; if it fails the sufficient-decrease test the previous '
              'linearization is restored, that step is halved and the iteration '
              'redone. Costs no extra propagation on the accept path.')
opts.register('stepBacktrackFromIter', 2, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'first Gauss-Newton iteration at which the chi2 backtracking test '
              'may fire (default 2). NOT 1: at iteration 0 the GBL '
              'propagation/kink residuals are identically zero by construction, '
              'so the iteration-0 chi2 is a different objective from every later '
              'one and comparing across that boundary would backtrack every '
              'candidate.')
opts.register('maxChi2Backtrack', 4, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'max chi2-backtracking halvings per accepted step (default 4). '
              'NOTE: distinct from the two/N-track maxBacktracks, which is the '
              'failed-propagation-leg retry budget.')
opts.register('armijoC', 1.e-4, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'Armijo sufficient-decrease coefficient c1 (default 1e-4).')
opts.register('armijoSlack', 1.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'relative chi2 slack in the Armijo test (default 1.0 = the chi2 may '
              'not more than DOUBLE in one iteration). This is a DIVERGENCE TRAP, '
              'not a line-search tolerance: the CVH/GBL iteration does NOT '
              'monotonically decrease r^T Vinv r (it drifts up ~0.3-0.5 per '
              'iteration even at 1/16 step), so a textbook 1e-4..1e-3 makes it '
              'halve the step forever -- 57 % of gun candidates, 15.6 halvings '
              'each, 6x the propagation cost, no change in the result. Scan in '
              'NOTES.md 2026-09-05.')
opts.register('clampMomentumFloor', -1.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'Gauss-Newton momentum floor [GeV] for the refit step clamp. '
              '<0 (default) = derive it from propagationPtotLimit as '
              '1.25*plimit. THE FLOOR MUST STAY ABOVE THE PROPAGATION LIMIT: '
              'the clamp exists only to keep the Gauss-Newton state out of '
              'the propagator refusal region, so it has to bracket that limit '
              'from above with a small margin. It must NOT be set any higher '
              'than that -- a floor above the physical momentum spectrum pins '
              'soft tracks at the floor (momentum-high, chi2/ndof >> 1) and, '
              'where p_ref is already below it, scales the step to zero and '
              'freezes the fit at its seed. The hard-coded 2.0 GeV floor '
              'against the 0.2 GeV limit did exactly that to 12 % of the '
              'flat-pT J/psi-gun candidates and carried the entire +0.21e-3 '
              'mass-scale offset (NOTES.md 2026-09-04).')
opts.register('doSimHits', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'read tracker PSimHits and match them to hits (input must keep them)')
opts.register('fitSimHitPositions', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'rung-E closure: fit simulated hit positions instead of cluster '
              'positions (covariances unchanged); requires doSimHits=True')
opts.register('trackSrc', 'ALCARECOTkAlJpsiX', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'input track collection (ALCARECOTkAlJpsiX for the custom '
              'B->JpsiX ALCARECO, ALCARECOTkAlJpsiMuMu for standard TkAl)')
opts.register('fitFromGenParms', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'freeze the reference parameters to the gen-muon values '
              '(gen-closure mode, default True)')
opts.register('fillJac', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'store per-track Jacobians (not needed for the grads solve)')
opts.register('fillGrads', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'store per-track gradient + packed Hessian (default True; '
              'input to fit_global_grads.py)')
opts.register('useIdealGeometry', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'use ideal (uncorrected) tracker geometry instead of the '
              'MC-production alignment from the GT (default False)')
opts.register('useScalarPot3D', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'use the spherical-harmonic scalar-potential field model (default)')
opts.register('useOpera3D', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'use the full 3D TOSCA volumetric grid (160812) as baseline field '
              '(takes precedence over useScalarPot3D)')
opts.register('useDefaultField', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'use the UNLABELLED default CMSSW field already loaded by '
              'MagneticField_cff: VolumeBasedMagneticField 160812 with '
              'useParametrizedTrackerField=True -> OAE_1103l_071212 inside the '
              'tracker. THIS IS THE CORRECT SETTING FOR THIS DRIVER on standard '
              'MC: the SIM propagates through OAE, so refitting with the full 3D '
              'grid (or ScalarPot3D) injects a SIM-vs-refit field difference that '
              'leaks an eta/phi-coherent shift into the gen-matched pull width -- '
              'measured 2026-08-07 as dp/p ~ 8e-4 across eta (0.155% of unit '
              'variance, a LOWER bound since OAE is phi-symmetric by construction). '
              'Takes precedence over useOpera3D and useScalarPot3D.')
opts.register('scalarPot3DInitFile', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'coefficient dump file produced by mfs/dump_coeffs_for_cmssw.py '
              '(always required: parmtype-14 registration needs it)')
opts.register('numberOfThreads', 1, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'framework numberOfThreads (numberOfStreams follows the same value)')
opts.register('propagationDirection', 'anyDirection', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'Geant4ePropagator PropagationDirection (anyDirection = per-leg '
              'forward/backward choice, default; alongMomentum = legacy '
              'forward-only)')
opts.register('keepPixelEdgeHits', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'keep pixel hits whose cluster touches the sensor boundary '
              '(default False = baseline exclusion; the 2022 attempt predates '
              'this cut, so the baseline is the interesting configuration)')
opts.register('pixelMinSizeX', 2, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'minimum pixel cluster size in x for a hit to stay in the fit '
              '(default 2 = baseline sizeX>1 cut; 1 admits all clusters)')
opts.register('globalTag', '106X_mcRun2_asymptotic_v17',
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'conditions GT. MUST match the GT the input MC was PRODUCED '
              'with: the fit re-evaluates hit positions with its own CPEs, '
              'so a mismatched SiPixelLorentzAngle / SiPixelTemplate payload '
              'shifts local-x per module and fakes a pixel hit-quality bias. '
              'Default 106X_mcRun2_asymptotic_v17 is right for the '
              'B->J/psi+X MC (produced in CMSSW_10_6_20_patch1). The private '
              'mu-gun simprod uses auto:run2_design -> 131X_mcRun2_design_v3, '
              'whose pixel templates are SiPixelTemplates38T_2010_2011_mc; '
              'pass globalTag=131X_mcRun2_design_v3 for those samples.')
opts.register('pixelMinSizeY', 1, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'minimum pixel cluster sizeY. Default 1 = NO cut, which is '
              'deliberate: sizeY==1 is the GEOMETRIC outcome for a track '
              'crossing perpendicular (sizeY ~ 1 + 1.9|cot theta|), unlike '
              'sizeX==1 which means the expected Lorentz sharing failed. '
              'Set to 2 only as an A/B DIAGNOSTIC -- it removes 12.6% of BPix '
              'and 34.6% of FPix hits, concentrated at central eta.')
opts.register('hitCovScalePixel', 1.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'multiplier on the assigned PIXEL hit covariance (1.0 = the CPE '
              'value as-is). Turns the dead scalecov=0.8 hypothesis into a '
              'measurable configuration.')
opts.register('hitCovScaleStrip', 1.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'multiplier on the assigned STRIP hit covariance (1.0 = as-is; '
              'the dead hypothesis was 1.2)')
opts.register('corFile', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'optional correction file (parmtree/x) applied via corparms_, '
              'e.g. fitted pixel-hit class corrections')
opts.register('nIters', 10, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'Gauss-Newton iteration cap (default 10; gen-anchored fits run 1)')
opts.register('edmConvergence', 1e-5, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'EDM convergence threshold on the reference-state block')
opts.register('gnDampAfter', 0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'damp Gauss-Newton steps from this iteration on (0 = off)')
opts.register('gnDampFactor', 0.5, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'damping factor applied to the GN step when gnDampAfter is active')
opts.register('debugPerIterDump', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'per-iteration debug trace; use with eventsToProcess on a few events')
opts.register('eventsToProcess', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'comma-separated run:event list to select specific events (empty = all)')
_defaultGroupsFile = os.path.join(os.environ.get('CMSSW_BASE', ''),
                                  'src/Analysis/HitAnalyzer/data/materialGroups50.txt')
opts.register('materialGroupsFile', _defaultGroupsFile, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'global material model grouping-tier rules file; empty = off')
opts.register('globalMaterialModel', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'parmtype-15 global material groups instead of per-module '
              'parmtype 7 (default True = current production model)')
opts.register('perStepFieldModes', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'per-Geant4-step field-mode attribution (default True)')
opts.register('skipHitlessSurfaces', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'drop hitless module surfaces from the fit (default True; '
              'effective only with globalMaterialModel=True)')
# The CVH energy-loss / process-noise switches, so this driver can select the
# ESTIMATOR as well as the sample. Without this hook the closure campaigns
# could only ever run the cfi default -- which since 9a7c692 is the CGF Fisher
# weight, making the legacy truncated-Q refit (CgfQoPMode=0) unreachable from
# the one driver every resolution study uses.
import TrackPropagation.Geant4e.cvhSwitches as cvhSwitches
cvhSwitches.register(opts)

opts.parseArguments()
if not opts.scalarPot3DInitFile:
    raise SystemExit(
        "scalarPot3DInitFile=<path> is required (coefficient dump file): "
        "the basis evaluator in globalCor needs it for parmtype-14 registration.")

process = cms.Process("BENCH", Run2_2016)

process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.EndOfProcess_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.GeometrySimDB_cff")

# Conditions the MC was produced with (CMSSW_10_6_20_patch1 production
# chain) -- alignment/CPE/beamspot consistent with the simulated detector,
# which is what a gen-closure fit must use.
process.GlobalTag = GlobalTag(process.GlobalTag, opts.globalTag, "")
process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(
        record=cms.string("GeometryFileRcd"),
        tag=cms.string("XMLFILE_Geometry_2016_81YV1_Extended2016_mc"),
        label=cms.untracked.string("Extended"),
    ),
)
process.XMLFromDBSource.label = cms.string("Extended")

process.load("TrackPropagation.Geant4e.geantRefit_cff")

# --- Geant4e field-integration precision in the FIT --------------------------
# Josh: "really really really important" for the CVH momentum scale -- and it
# was only ever applied to the SIMULATION. geantRefit_cff builds geopro with
#     MagneticField = _g4SimHits.MagneticField.clone()
# i.e. the CFI DEFAULTS, which nothing here overrode. Measured 2026-08-08:
#     DeltaOneStepTracker      1e-4   (our SIM: 1e-5)   10x looser
#     DeltaIntersectionTracker 1e-6   (our SIM: 1e-6)   same
#     DeltaOneStep             1e-3   (our SIM: 1e-5)  100x looser
#     DeltaIntersection        1e-4   (our SIM: 1e-6)  100x looser
# So every CVH refit propagated at 10-100x looser precision than the simulation
# it is compared against. A chord error is a COHERENT trajectory displacement,
# not a random one, so it biases the momentum rather than broadening it -- the
# signature of the unexplained ~1e-4 single-track scale offset and the +4.5e-4
# J/psi two-track mass bias.
# Default True: matching the SIM is the physically defensible choice, and the
# flag exists so the change can be measured rather than assumed.
if opts.tightG4eStepper:
    _fsp = process.geopro.MagneticField.ConfGlobalMFM.OCMS.StepperParam
    import os as _os
    if _os.environ.get("CVH_ABSURD_STEPPER"):
        # DIAGNOSTIC: deliberately absurd tolerances. If the output is still
        # bit-identical, geopro is not scheduled / these params do not reach
        # the CVH propagation path at all.
        _fsp.DeltaOneStepTracker = 1.0
        _fsp.DeltaIntersectionTracker = 1.0
        _fsp.DeltaOneStep = 1.0
        _fsp.DeltaIntersection = 1.0
        print(">>> ABSURD STEPPER SET:", _fsp.DeltaOneStep.value())
    else:
        _fsp.DeltaOneStepTracker = 1e-5
    _fsp.DeltaIntersectionTracker = 1e-6
    _fsp.DeltaOneStep = 1e-5
    _fsp.DeltaIntersection = 1e-6
from TrackPropagation.Geant4e.cvhMaster_cfi import CvhMasterPSet

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(opts.nEvents))

_paths = [p.strip() for p in opts.input.split(',') if p.strip()]
if opts.inputFileList:
    with open(opts.inputFileList) as _f:
        _paths += [l.strip() for l in _f if l.strip() and not l.startswith('#')]
assert _paths, "must set input=<paths> and/or inputFileList=<file> on the cmsRun command line"
_urls = [p if p.startswith(("root://", "file:")) else "file:" + p for p in _paths]
process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(*_urls),
    secondaryFileNames=cms.untracked.vstring(),
    # The condor MC production has a small tail of corrupt files; skip them
    # instead of aborting (pre-scan the filelist for the
    # FormatIncompatibility cases, see calibration_studies/pixelhits).
    skipBadFiles=cms.untracked.bool(True),
    # Every condor job numbers its events from the same (run=1, lumi=1)
    # range, so distinct physics events collide in (run, lumi, event) and
    # the default duplicate check silently drops most of the sample.
    duplicateCheckMode=cms.untracked.string('noDuplicateCheck'),
)

if opts.eventsToProcess:
    process.source.eventsToProcess = cms.untracked.VEventRange(
        *[s.strip() for s in opts.eventsToProcess.split(',') if s.strip()])

process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(int(opts.numberOfThreads)),
    numberOfStreams=cms.untracked.uint32(int(opts.numberOfThreads)),
    numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Per-stream CLHEP engine for the residual-maker (MT-safe CvhMaster path).
process.RandomNumberGeneratorService.globalCor = cms.PSet(
    initialSeed=cms.untracked.uint32(123456789),
    engineName=cms.untracked.string('HepJamesRandom'),
)

process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")

_PARTICLES = {
    'mu':     dict(pdg=13,   g4=["mu+", "mu-"]),
    'pi':     dict(pdg=211,  g4=["pi+", "pi-"]),
    'kaon':   dict(pdg=321,  g4=["kaon+", "kaon-"]),
    'proton': dict(pdg=2212, g4=["proton", "anti_proton"]),
}
if opts.particle not in _PARTICLES:
    raise SystemExit(f"unknown particle={opts.particle} (use mu/pi/kaon/proton)")
_pcfg = _PARTICLES[opts.particle]

# Single-track CVH refit of the gen-matched muon tracks. requireGen with
# the dR<0.1 status-1 |pdgId|==13 matching drops the other B daughters
# (kaon/pion tracks) of the TkAlJpsiX track collection.
process.globalCor = cms.EDProducer(
    "ResidualGlobalCorrectionMakerG4e",
    src=cms.InputTag(opts.trackSrc),
    fitFromGenParms=cms.bool(bool(opts.fitFromGenParms)),
    fitFromSimParms=cms.bool(False),
    fillTrackTree=cms.bool(True),
    fillGrads=cms.bool(bool(opts.fillGrads)),
    fillJac=cms.bool(bool(opts.fillJac)),
    fillRunTree=cms.bool(True),
    doGen=cms.bool(True),
    genParticles=cms.InputTag("genParticles"),
    pileupInfo=cms.InputTag("addPileupInfo"),
    doSim=cms.bool(bool(opts.doSimHits)),
    fitSimHitPositions=cms.untracked.bool(bool(opts.fitSimHitPositions)),
    requireGen=cms.bool(True),
    doMuons=cms.bool(False),
    doMuonAssoc=cms.bool(False),
    doTrigger=cms.bool(False),
    doRes=cms.bool(bool(opts.doRes)),
    exportStepRecords=cms.bool(bool(opts.exportStepRecords)),
    exportCfExponents=cms.bool(bool(opts.exportCfExponents)),
    exportCfGroupExponents=cms.bool(bool(opts.exportCfGroupExponents)),
    exportHitResBlocks=cms.bool(bool(opts.exportHitResBlocks)),
    exportMaterialNoise=cms.bool(bool(opts.exportMaterialNoise)),
    exportVarianceGrads=cms.bool(bool(opts.exportVarianceGrads)),
    varianceGradFamilies=cms.vuint32(*[int(x) for x in opts.varianceGradFamilies]),
    exportObjective=cms.bool(bool(opts.exportObjective)),

    useIdealGeometry=cms.bool(bool(opts.useIdealGeometry)),
    bsConstraint=cms.bool(False),
    applyHitQuality=cms.bool(True),
    keepPixelEdgeHits=cms.bool(bool(opts.keepPixelEdgeHits)),
    pixelMinSizeX=cms.int32(int(opts.pixelMinSizeX)),
    pixelMinSizeY=cms.int32(int(opts.pixelMinSizeY)),
    corFiles=cms.vstring(*( [opts.corFile] if opts.corFile else [] )),
    triggers=cms.vstring(),
    hitCovScalePixel=cms.double(float(opts.hitCovScalePixel)),
    hitCovScaleStrip=cms.double(float(opts.hitCovScaleStrip)),
    trackParticleName=cms.string(opts.particle),
    genMatchPdgId=cms.int32(_pcfg['pdg']),
    MagneticFieldLabel=cms.string(""),
    scalarPotentialInitFile=cms.string(opts.scalarPot3DInitFile),
    nIters=cms.uint32(int(opts.nIters)),
    edmConvergence=cms.double(float(opts.edmConvergence)),
    debugPerIterDump=cms.bool(bool(opts.debugPerIterDump)),
    gnDampAfter=cms.uint32(int(opts.gnDampAfter)),
    gnDampFactor=cms.double(float(opts.gnDampFactor)),
    runFDClosure=cms.bool(False),
    epsilonFDClosure=cms.double(1e-4),
    materialGroupsFile=cms.string(opts.materialGroupsFile),
    materialFDGroup=cms.int32(int(opts.materialFDGroup)),
    materialFDEps=cms.double(float(opts.materialFDEps)),
    globalMaterialModel=cms.bool(bool(opts.globalMaterialModel)),
    perStepFieldModes=cms.bool(bool(opts.perStepFieldModes)),
    skipHitlessSurfaces=cms.bool(bool(opts.skipHitlessSurfaces) and bool(opts.globalMaterialModel)),
    outprefix=cms.untracked.string("globalcor_resclosure"),
    # MT G4Error master (GlobalCache); particle set follows the hypothesis.
    CvhMaster=CvhMasterPSet.clone(Particles=cms.vstring(*_pcfg['g4'])),
)

# B-field model wiring (identical to the other 15_0 drivers).
if opts.useDefaultField:
    # Consume the unlabelled field MagneticField_cff already placed in the
    # EventSetup (OAE inside the tracker) -- the same field the SIM used.
    # Nothing to instantiate: the empty label IS the default producer's
    # label, and the CPEs keep their default (empty) label too, so the
    # Lorentz drift matches as well.
    fieldlabel = ""
elif opts.useOpera3D:
    from MagneticField.Engine.volumeBasedMagneticField_160812_cfi import \
        VolumeBasedMagneticFieldESProducer as Opera3DMagneticFieldProducer
    from MagneticField.Engine.volumeBasedMagneticField_160812_cfi import magfield as MagneticFieldGeometry
    process.magfield = MagneticFieldGeometry
    process.es_prefer_magfield_cvhrefit = cms.ESPrefer("XMLIdealGeometryESSource", "magfield")
    process.Opera3DMagneticFieldProducer = Opera3DMagneticFieldProducer
    fieldlabel = "grid_160812_3_8t"
    process.Opera3DMagneticFieldProducer.label = fieldlabel
    process.Opera3DMagneticFieldProducer.useParametrizedTrackerField = cms.bool(False)
    for _cpe in ("stripCPEESProducer", "StripCPEfromTrackAngleESProducer",
                 "siPixelTemplateDBObjectESProducer", "templates"):
        if hasattr(process, _cpe):
            getattr(process, _cpe).MagneticFieldLabel = fieldlabel
elif not opts.useScalarPot3D:
    raise RuntimeError(
        "useScalarPot3D=False is no longer supported; use the ScalarPot3D model.")
else:
    from MagneticField.ParametrizedEngine.parametrizedMagneticField_ScalarPot3D_cfi \
        import ParametrizedMagneticFieldProducer as ScalarPot3DMagneticFieldProducer
    process.ScalarPot3DMagneticFieldProducer = ScalarPot3DMagneticFieldProducer.clone()
    process.ScalarPot3DMagneticFieldProducer.parameters.InitFile = opts.scalarPot3DInitFile
    fieldlabel = "ScalarPot3DMf"
    process.ScalarPot3DMagneticFieldProducer.label = fieldlabel
process.geopro.MagneticFieldLabel = fieldlabel
process.Geant4ePropagator.MagneticFieldLabel = fieldlabel
from TrackPropagation.Geant4e.cvhMasterESProducer_cfi import cvhMasterESProducer
process.cvhMasterESProducer = cvhMasterESProducer.clone()
process.cvhMasterESProducer.MagneticFieldLabel = cms.string(fieldlabel)
# CVH-specific propagator path: instantiates the custom fluct
# (G4UniversalFluctuationForExtrapolator). Without this, computeErrorIoni
# dereferences a null fluct->table on the first event.
process.Geant4ePropagator.ForCVH = cms.bool(True)
process.Geant4ePropagator.PropagationDirection = cms.string(opts.propagationDirection)
# Momentum floor for the G4e propagation. The cfi default (1.0 GeV) drops
# ~2.5% of the flat-pT-gun daughters (fail[plimit], median seed p 0.65);
# 0.2 GeV matches the ditrack/V0 configuration and recovers them.
process.Geant4ePropagator.PropagationPtotLimit = cms.double(float(opts.propagationPtotLimit))
# Gauss-Newton momentum floor for the refit step clamp, derived from the
# propagation limit unless given explicitly (see the option's help). The pair
# is echoed because the two are only correct together.
_clampFloor = (float(opts.clampMomentumFloor) if float(opts.clampMomentumFloor) > 0.
               else 1.25 * float(opts.propagationPtotLimit))
if _clampFloor <= float(opts.propagationPtotLimit):
    raise RuntimeError(
        "clampMomentumFloor (%g GeV) must be ABOVE propagationPtotLimit (%g GeV): "
        "the Gauss-Newton clamp exists to keep the state out of the propagator's "
        "refusal region." % (_clampFloor, float(opts.propagationPtotLimit)))
process.globalCor.clampMomentumFloor = cms.double(_clampFloor)
print("[cvh] effective: PropagationPtotLimit=%g GeV, clampMomentumFloor=%g GeV"
      % (float(opts.propagationPtotLimit), _clampFloor))
# Relative step damping + chi2 backtracking (2026-09-05 replacement for the
# bare momentum floor; see the option help and NOTES.md). Defaults are the new
# behaviour; maxMomentumStepFactor<=1 + stepBacktracking=False reproduce the
# legacy clamp bit-identically.
process.globalCor.maxMomentumStepFactor = cms.double(float(opts.maxMomentumStepFactor))
process.globalCor.stepBacktracking = cms.bool(bool(opts.stepBacktracking))
process.globalCor.maxChi2Backtrack = cms.uint32(int(opts.maxChi2Backtrack))
process.globalCor.stepBacktrackFromIter = cms.uint32(int(opts.stepBacktrackFromIter))
process.globalCor.armijoC = cms.double(float(opts.armijoC))
process.globalCor.armijoSlack = cms.double(float(opts.armijoSlack))
print("[cvh] effective: maxMomentumStepFactor=%g, stepBacktracking=%s "
      "(fromIter=%d, maxChi2Backtrack=%d, armijoC=%g, armijoSlack=%g)"
      % (float(opts.maxMomentumStepFactor), bool(opts.stepBacktracking),
         int(opts.stepBacktrackFromIter), int(opts.maxChi2Backtrack),
         float(opts.armijoC), float(opts.armijoSlack)))
process.globalCor.MagneticFieldLabel = cms.string(fieldlabel)

# After every explicit propagator assignment above, so a command-line switch
# wins over the driver's own defaults and the effective state is echoed once.
cvhSwitches.apply(process, opts)

process.reconstruction_step = cms.Path(process.offlineBeamSpot * process.globalCor)
process.schedule = cms.Schedule(process.reconstruction_step)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
