## CVH two-track (dimuon) refit driven straight off UL16 MiniAODv2.
##
## Production driver for the Z -> mumu leg of the CVH full-scale feasibility
## test.  The track source is the pat::Muon tracker tracks of `slimmedMuons`
## (TrackProducerFromPatMuons), i.e. exactly the path the WMass custom NanoAOD
## uses (PhysicsTools/NanoAOD muons_cff: tracksfrommuons ->
## diMuonTrackVertexCandidates -> trackrefitdimuon).
##
## It is the MiniAOD counterpart of runCvhJpsiGenMC.py and exposes the same
## production switches -- the two drivers must be able to write outputs that
## are POOLABLE in one global fit, which means the same global parameter
## catalogue (same doRes / globalMaterialModel / materialGroupsFile / 50-mode
## scalar-potential dump), so those options live here with the same names and
## the same defaults.
##
## Grown from calibration_studies/transfer/readback_check/runCvhDimuonMiniAOD.py
## (a readback smoke) plus everything the profiling driver
## calibration_studies/production/profiling/runCvhProfile.py exposed.
##
## MC -> default CMSSW field (no ScalarPot3D override), per CLAUDE.md: the SIM
## propagated through the OAE-parametrised tracker field, so a closure/feasibility
## fit has to use the same field.
##
## Typical production invocation (see calibration_studies/production/
## config_dymc8p5M.sh, which is the authority):
##
##   cmsRun runCvhDimuonMiniAOD.py \
##     input=<one MiniAOD file> skipEvents=0 nEvents=25000 numberOfThreads=1 \
##     doRes=True exportCfExponents=True exportStepRecords=False \
##     fillJac=True fillGrads=False fillGradsFactored=True \
##     fitFromGenParms=False doSimHits=False doGen=True requireGen=False \
##     useIdealGeometry=False useDefaultField=True \
##     globalTag=106X_mcRun2_asymptotic_v17 \
##     doMassConstraint=False massMin=60 massMax=120 \
##     CgfQoPMode=0 tightG4eStepper=True \
##     propagationPtotLimit=0.2 maxMomentumStepFactor=2.0 stepBacktracking=True \
##     scalarPot3DInitFile=<...custom50.txt>
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os

from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
from Configuration.AlCa.GlobalTag import GlobalTag

opts = VarParsing.VarParsing('analysis')

# ---------------------------------------------------------------- input ----
opts.register('input', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'comma-separated absolute paths or root:// URLs of MiniAOD files')
opts.register('inputFileList', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'text file with one input path per line; combined with input= '
              'if both are given')
opts.register('nEvents', 500, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, 'number of events (-1 = all)')
opts.register('skipEvents', 0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'skip the first N events of the input (PoolSource skipEvents). '
              'Combined with nEvents this splits ONE input file across several '
              'batch tasks -- an 80k-event MiniAOD file is ~7 h of CVH in a '
              'single job, so production chunks it (skipEvents=k*C nEvents=C). '
              'NOTE: PoolSource counts the skip across the WHOLE fileNames '
              'list, so chunking is only well-defined with a single input '
              'file per task; this is asserted below.')
opts.register('eventsToProcess', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'comma-separated run:event list to select specific events; '
              'empty = all')
opts.register('muonSrc', 'slimmedMuons', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'pat::Muon collection feeding TrackProducerFromPatMuons')
opts.register('innerTrackOnly', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'TrackProducerFromPatMuons: use only the inner track')
opts.register('muonPtMin', -1.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'TrackProducerFromPatMuons pt threshold (-1 = no cut)')
opts.register('massMin', 60.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'min mu-mu mass for the candidate producer (Z window: 60)')
opts.register('massMax', 120.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'max mu-mu mass for the candidate producer (Z window: 120)')
opts.register('numberOfThreads', 1, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'framework numberOfThreads (numberOfStreams follows). 1 also '
              'makes the output exactly <outprefix>_0.root')
opts.register('outprefix', 'globalcor', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'output tree file prefix; the file is <outprefix>_<stream>.root')

# ---------------------------------------------------------- conditions ----
opts.register('globalTag', '106X_mcRun2_asymptotic_v17',
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'conditions. MUST match the sample: 106X_mcRun2_asymptotic_v17 '
              'is the UL16 MiniAODv2 chain. Getting this wrong is not '
              'cosmetic -- a mismatched GT supplies the wrong pixel templates '
              'and alignment for the simulated detector.')
opts.register('useIdealGeometry', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'ideal tracker geometry instead of the aligned one from the GT '
              '(default False = aligned, the real-geometry configuration)')

# ---------------------------------------------------------------- gen -----
opts.register('doGen', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'read gen particles / gen weight / pileup and fill the gen '
              'branches (MC only)')
opts.register('genParticles', 'prunedGenParticles',
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'gen particle collection. MiniAOD RETAG: the AOD-level '
              '"genParticles" does not exist in MiniAOD, the pruned copy '
              'does. Z muons (pt ~ 20-60 GeV, status 1) are always kept by '
              'the standard pruning.')
opts.register('pileupInfo', 'slimmedAddPileupInfo',
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'PileupSummaryInfo collection. MiniAOD RETAG of "addPileupInfo". '
              'Dereferenced UNGUARDED by the maker when doGen=True, so it has '
              'to be right.')
opts.register('requireGen', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'drop candidates whose tracks do not gen-match a status-1 muon '
              'within dR<0.1 and the same charge. Default False, UNLIKE the '
              'J/psi driver: there the all-pairs loop needs it to reject '
              'combinatorics, whereas an opposite-sign muon pair in a 60-120 '
              'GeV window on DY is the Z by construction. Keeping it False '
              'also removes any dependence of the yield on the MiniAOD gen '
              'pruning thresholds. The Mu{plus,minus}gen_* branches are '
              'filled either way, so the cut can be made offline.')
opts.register('fitFromGenParms', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'freeze the 10-dim vertex/kinematic reference block to the gen '
              'muon values (gen-closure mode). Default False = a '
              'reconstruction-level fit. NOTE the freeze is NOT conditioned '
              'on the gen match succeeding, so True is a physics choice, not '
              'a fallback.')
opts.register('doSimHits', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'read tracker PSimHits (MiniAOD has none -- keep False)')

# ------------------------------------------------------------ triggers ----
opts.register('doTrigger', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'store the per-path HLT decisions of the Z trigger list below')
opts.register('applyHltFilter', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'pre-filter events on the Z trigger list (default False)')

# ------------------------------------------------------------- outputs ----
opts.register('fillJac', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool, 'store per-track Jacobians')
opts.register('fillGrads', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'store per-event gradient + PACKED Hessian (hesspackedv). NOTE '
              'this does not save the packing CPU when off -- the packed '
              'Hessian is computed unconditionally, the flag only controls '
              'whether it is written.')
opts.register('fillGradsFactored', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'store per-event gradient + LOW-RANK factored Hessian '
              '(H = B^T B as nRank/nFactor/hessfactorv); ~9x smaller than '
              'hesspackedv at this mode count. Not mutually exclusive with '
              'fillGrads -- both would be written.')
opts.register('doRes', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'register the resolution families (parmtypes 8/9/10/11) and '
              'enable the CF export block. NOTE the in-maker CF export is '
              'gated on doRes AND (fillGrads OR fillGradsFactored): each is '
              'nearly free alone and the pair costs ~10x (0.06 -> 0.63 '
              's/candidate on Z), because only the pair turns the export on.')
opts.register('exportStepRecords', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'write the RAW per-step resolution export (ioniurbanv, msmoliv, '
              'radstep*, reseigv, resinf*). ~290 kB/candidate on Z (317 vs 28 '
              'kB/cand measured), while the exponents it feeds are written '
              'directly by the in-maker cfmass_* export, so it is off here; '
              'set True to re-derive the exponents under a different model. '
              'Always passed explicitly: the maker cfi defaults it True when '
              'the parameter is absent.')
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
              'tr(dV R dV R). Two-track maker only. With it OFF a parameter '
              'that moves the covariance (parmtype 15 through exp(k_g), and '
              '8/9/10/11 entirely) enters the quadratic hit-chi2 term only '
              'through the MEAN')
opts.register('varianceGradFamilies', [], VarParsing.VarParsing.multiplicity.list,
              VarParsing.VarParsing.varType.int,
              'which parmtypes exportVarianceGrads covers; empty = '
              '{8,9,10,11,15}. 15 alone is LAYOUT-PRESERVING: the '
              'material-group globals are already columns of globalidxv, so '
              'only their VALUES change. 8/9/10/11 are per-module and APPEND '
              'columns to globalidxv and to every array indexed by it')
opts.register('exportObjective', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'write objval/objchisq/objlogdetv/objlogdetc, the marginal '
              'objective r^T R r + ln|V| + ln|C| in double precision. '
              'Validation only -- it costs an ncons x ncons LDLT per '
              'candidate -- and exists so the gradient can be '
              'finite-differenced against what it claims to differentiate')
opts.register('varianceFDGlobalIdx', -1, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'IN-MAKER finite difference of the variance gradient at FIXED '
              'linearization: >=0 does that one global index, -2 does every '
              'enabled variance column of families 10/11/15. Perturbs '
              'V -> V + s dV_i with r/F/J held fixed and re-does the profile, '
              'so it tests the assembly (traces, projector, sign, ln|C|) to '
              'O(s^2). Prints VARFD lines')
opts.register('varianceFDEps', 1e-3, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'the s of varianceFDGlobalIdx')
opts.register('exportHitResBlocks', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'register the parmtype-8/9 HIT-RESOLUTION dV blocks in the '
              'influence export (reseigidx/resinfvarv/reshitcls + the '
              'cf*_hitcls/cf*_hitv per-class shares); they are what the '
              'per-hit-class resolution parameters are fitted from. Export '
              'only -- it cannot move the fit. Set False to leave them out of '
              'the tree')
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
              'compute the resolution-CF exponents in the maker and write '
              'them on the 64-point tau grid (cfmass_* for the two-track fit)')
opts.register('produceValueMaps', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'emit the per-candidate EDM ValueMaps. Default False: this '
              'driver has no output module, so nothing would consume them.')

# --------------------------------------------------- global parameters ----
opts.register('scalarPot3DInitFile', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'coefficient dump file produced by mfs/dump_coeffs_for_cmssw.py. '
              'ALWAYS REQUIRED: the maker uses it to register the parmtype-14 '
              'modes and their Jacobian columns even when another model '
              'supplies the field. Use the 50-mode ("custom50") dump -- the '
              '360-mode one costs 0.21 s/candidate in the per-step basis '
              'evaluation against 0.013 s at 50 modes, and would give a '
              'DIFFERENT global parameter catalogue from the J/psi leg.')
_defaultGroupsFile = os.path.join(os.environ.get('CMSSW_BASE', ''),
                                  'src/Analysis/HitAnalyzer/data/materialGroups50.txt')
opts.register('materialGroupsFile', _defaultGroupsFile,
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'global material model grouping-tier rules file (42 groups from '
              'materialGroups50.txt). MUST match the J/psi leg for the two '
              'outputs to be poolable; empty = off')
opts.register('globalMaterialModel', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'replace the per-module material parameters (parmtype 7) with '
              'the parmtype-15 global material groups (exclusive switch)')
opts.register('perStepFieldModes', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'apply the scalar-potential correction and attribute the '
              'per-mode derivatives per Geant4 step rather than piecewise-'
              'constant per leg')
opts.register('skipHitlessSurfaces', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'drop hitless module surfaces from the fit; effective only with '
              'globalMaterialModel=True (auto-disabled otherwise)')

# ------------------------------------------------------------- the fit ----
opts.register('doVtxConstraint', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'apply the common-vertex constraint in the two-track fit: state '
              'index 6, the signed track-track PCA distance, is frozen at zero '
              'and the fitted mass is the vertex-constrained one. Jpsi_mass_unc '
              'carries the unconstrained mass, so either can be formed offline. '
              'False leaves index 6 free (a plain two-track fit through the PCA)')
opts.register('minLegHits', 8, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'minimum valid hits on the WEAKER leg. Default 8, 0 = off: a '
              'thin leg is BACKGROUND -- by gen truth on DY, 82-93 % of what '
              'this removes is a duplicate or unmatched pairing '
              '(0.885 +- 0.026) at 0.9983 +- 0.0004 signal efficiency -- and '
              'every candidate with a non-finite exported mass resolution in '
              'dy_vtxon has a leg of one or two hits while its PAIR total is '
              '13-23, so a pair-sum cut cannot see them and this can.')
opts.register('minNdof', 1, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'minimum degrees of freedom of the two-track fit, required '
              'BEFORE the fit. ndof = nvalid + nvalidpixel - 10 (+3 beamspot, '
              '+1 pointing, +1 vertex constraint), i.e. one coordinate per '
              'strip hit and two per pixel hit against the ten state '
              'parameters the common vertex costs. The default 1 means more '
              'than nine measurement coordinates with the vertex constraint '
              'on and more than ten with it off: at ndof == 0 the fit is '
              'exactly determined (chi2 identically zero, chi2/ndof 0/0) and '
              'the factored-Hessian export aborts the process. 0 disables.')
opts.register('minPairHits', -1, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'minimum number of VALID HITS summed over the two legs (pixel '
              'hits counted once). -1 = auto = 10 with the vertex constraint '
              'on, 11 with it off -- the same requirement as minNdof read on '
              'hits rather than on measurement coordinates. 0 disables.')
opts.register('exportVtxResidual', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'export the VERTEX-CONSTRAINT RESIDUAL (state index 6, the '
              'track-track PCA distance) as a CF resolution term')
opts.register('doMassConstraint', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'apply the dimuon mass constraint (adds a second nicons pass)')
opts.register('massConstraint', 91.1876, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float, 'mass constraint value [GeV]')
opts.register('massConstraintWidth', 2.4952, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float, 'mass constraint width [GeV]')
opts.register('nIters', 10, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'Gauss-Newton iteration cap per constraint phase')
opts.register('edmConvergence', 1e-5, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'EDM convergence threshold on the reference-state block '
              '(0 disables early stopping)')
opts.register('keepPixelEdgeHits', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'keep pixel hits whose cluster touches the sensor boundary')
opts.register('pixelMinSizeX', 2, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'minimum pixel cluster size in x for a hit to stay in the fit')

# ---------------------------------------------------- propagation / G4 ----
opts.register('tightG4eStepper', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'tighten the Geant4e field-integration tolerances in the FIT to '
              'match the simulation. See the block after geantRefit_cff for '
              'WHICH PSet has to be patched in this release -- it is NOT '
              'geopro any more.')
opts.register('propagationDirection', 'anyDirection',
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'Geant4ePropagator PropagationDirection. "anyDirection" picks '
              'forward/backward per leg from the target-plane geometry, '
              'recovering legs whose target plane is marginally behind the '
              'state; "alongMomentum" is the legacy forward-only behaviour.')
opts.register('propagationPtotLimit', 0.2, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'G4e propagation momentum floor [GeV]; the cfi default is 1.0')
opts.register('clampMomentumFloor', -1.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'Gauss-Newton momentum floor [GeV] for the refit step clamp. '
              '<0 (default) = derive it from propagationPtotLimit as '
              '1.25*plimit. THE FLOOR MUST STAY ABOVE THE PROPAGATION LIMIT '
              'and no higher: a floor above the physical momentum spectrum '
              'pins soft tracks at the floor and, where p_ref is already '
              'below it, freezes the fit at its seed.')
opts.register('maxMomentumStepFactor', 2.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'RELATIVE Gauss-Newton step damping: the max factor by which a '
              'track momentum may change in one iteration (2 = p may at most '
              'halve or double). <=1 switches it off and leaves '
              'clampMomentumFloor as the only bound.')
opts.register('stepBacktracking', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'chi2-based (Armijo) retroactive step backtracking; costs no '
              'extra propagation on the accept path')
opts.register('stepBacktrackFromIter', 2, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'first Gauss-Newton iteration at which the chi2 backtracking '
              'test may fire. NOT 1: the iteration-0 chi2 is a different '
              'objective (the GBL propagation/kink residuals are identically '
              'zero there).')
opts.register('maxChi2Backtrack', 4, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'max chi2-backtracking halvings per accepted step')
opts.register('armijoC', 1.e-4, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'Armijo sufficient-decrease coefficient c1')
opts.register('armijoSlack', 1.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'relative chi2 slack in the Armijo test (1.0 = the chi2 may not '
              'more than DOUBLE in one iteration). This is a DIVERGENCE TRAP, '
              'not a line-search tolerance.')

# ------------------------------------------------------------- B field ----
opts.register('useDefaultField', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'use the UNLABELLED default CMSSW field, i.e. what '
              'MagneticField_cff already loaded: VolumeBasedMagneticField '
              '160812 with useParametrizedTrackerField=True -> the '
              'OAE_1103l_071212 tracker parametrization. DEFAULT TRUE HERE '
              'because this driver runs on standard MC: the SIM propagated '
              'through OAE, and refitting with the full 3D grid instead '
              'leaks an eta/phi-coherent ~8e-4 dp/p pattern into the pull '
              'width. Takes precedence over useOpera3D and useScalarPot3D.')
opts.register('useOpera3D', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'use the full 3D TOSCA volumetric grid (160812) as the baseline '
              'field; reserved as a deliberate injected-field systematic on a '
              'subset, not as the baseline for a closure test.')
opts.register('useScalarPot3D', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'use the spherical-harmonic scalar-potential field model in the '
              'refit (only used when useDefaultField and useOpera3D are both '
              'False)')

# The CVH physics/estimator switches (CgfQoPMode, IoniExactDelta, ...) as
# command-line options, so this driver PINS the estimator instead of silently
# inheriting the cfi default. This is the switch that decided the Z cost:
# under CgfQoPMode>=1 the propagator runs cvhcgf::inverseFisher on every
# propagate call (two 262144-point FFTs) at ~88 calls/candidate, i.e.
# 8.5 s/candidate against 0.63 s at mode 0 -- and the two-track maker has no
# setCgfOverride hook, so it pays for a block it never uses.
import TrackPropagation.Geant4e.cvhSwitches as cvhSwitches
cvhSwitches.register(opts)

opts.parseArguments()

# TWO-TRACK DEFAULT IS THE Q-MATRIX ESTIMATOR. Under CgfQoPMode >= 1 the
# fluctuation model returns the UNTRUNCATED ionization second cumulant as the
# leg's q/p variance on the assumption that the maker substitutes the Fisher
# weight through setCgfOverride. Only the single-track maker does; this maker
# has no hooks, so a non-zero mode would run the fit with an inflated
# ionization variance AND pay ~14x for a block it never uses.
if opts.CgfQoPMode < 0:
    opts.CgfQoPMode = 0
    print('[cvh] two-track driver: CgfQoPMode not given, defaulting to 0 '
          '(legacy truncated-Q); the two-track maker has no CGF override hooks')

assert opts.input or opts.inputFileList, \
    "must set input=<paths> and/or inputFileList=<file>"
if not opts.scalarPot3DInitFile:
    raise SystemExit(
        "scalarPot3DInitFile=<path> is required (coefficient dump file): "
        "the basis evaluator in globalCor needs it for chain-rule columns "
        "even when another model supplies the baseline field.")

# Single-muon paths of the 2016 menu; only read when doTrigger=True. Missing
# paths record False rather than throwing, so an era mismatch is not fatal.
Z_TRIGGERS = [
    "HLT_IsoMu24",
    "HLT_IsoTkMu24",
    "HLT_Mu50",
    "HLT_TkMu50",
    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
    "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",
]

process = cms.Process("CVHDIMUMINI", Run2_2016)
process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.GeometrySimDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

# Conditions the MC was produced with (106X UL2016 MC chain), plus the 2016
# DDD geometry payload the CVH G4 world is built from.
process.GlobalTag = GlobalTag(process.GlobalTag, opts.globalTag, "")
process.GlobalTag.toGet = cms.VPSet(cms.PSet(
    record=cms.string("GeometryFileRcd"),
    tag=cms.string("XMLFILE_Geometry_2016_81YV1_Extended2016_mc"),
    label=cms.untracked.string("Extended"),
))
process.XMLFromDBSource.label = cms.string("Extended")

process.load("TrackPropagation.Geant4e.geantRefit_cff")

# --- Geant4e field-integration precision in the FIT --------------------------
# Josh: "really really really important" for the CVH momentum scale -- it was
# only ever applied to the SIMULATION. The CMSSW cfi defaults are
#     DeltaOneStepTracker      1e-4   (our SIM: 1e-5)   10x looser
#     DeltaIntersectionTracker 1e-6   (our SIM: 1e-6)   same
#     DeltaOneStep             1e-3   (our SIM: 1e-5)  100x looser
#     DeltaIntersection        1e-4   (our SIM: 1e-6)  100x looser
# A chord error is a COHERENT trajectory displacement, not a random one, so it
# biases the momentum rather than broadening it.
#
# WHICH PSet: in this release the CVH G4 world and its field manager are built
# by CvhMaster/CvhWorker from `cvhMasterESProducer.MagneticField`
# (CvhMaster.cc: `m_pField(p.getParameter<edm::ParameterSet>("MagneticField"))`
# -> sim::FieldBuilder; CvhWorker.cc does the same per thread). `geopro` is
# neither scheduled nor consumed here, and geantRefit_cff gives it an
# INDEPENDENT clone of the same g4SimHits PSet, so patching geopro alone has
# no effect on the propagation. Both are patched below: cvhMasterESProducer
# because it is the one that acts, geopro so the two never disagree if it is
# ever put back on a path.
_TIGHT_STEPPER = dict(DeltaOneStepTracker=1e-5, DeltaIntersectionTracker=1e-6,
                      DeltaOneStep=1e-5, DeltaIntersection=1e-6)


def _tighten(fieldPSet):
    for _k, _v in _TIGHT_STEPPER.items():
        setattr(fieldPSet.ConfGlobalMFM.OCMS.StepperParam, _k, cms.double(_v))


process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(opts.nEvents))

_paths = [p.strip() for p in opts.input.split(',') if p.strip()]
if opts.inputFileList:
    with open(opts.inputFileList) as _f:
        _paths += [l.strip() for l in _f if l.strip() and not l.startswith('#')]
assert _paths, "must set input=<paths> and/or inputFileList=<file>"
# Accept local paths (prepend "file:") or xrootd URLs as-is.
_urls = [p if p.startswith(("root://", "file:")) else "file:" + p for p in _paths]
process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(*_urls),
    secondaryFileNames=cms.untracked.vstring(),
    # Tolerate an unreadable input rather than aborting a many-file job.
    # WARNING: a skipped input is SILENT -- the job then writes a valid but
    # EMPTY output. The batch wrapper guards this by failing on
    # "fit summary attempted=0".
    skipBadFiles=cms.untracked.bool(True),
    duplicateCheckMode=cms.untracked.string('noDuplicateCheck'),
)

if int(opts.skipEvents) > 0:
    assert len(_urls) == 1, (
        "skipEvents is only well-defined with exactly one input file per job "
        "(PoolSource skips across the concatenated fileNames list); got %d"
        % len(_urls))
    process.source.skipEvents = cms.untracked.uint32(int(opts.skipEvents))

if opts.eventsToProcess:
    process.source.eventsToProcess = cms.untracked.VEventRange(
        *[s.strip() for s in opts.eventsToProcess.split(',') if s.strip()])

process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(int(opts.numberOfThreads)),
    numberOfStreams=cms.untracked.uint32(int(opts.numberOfThreads)),
    numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# --- B field ----------------------------------------------------------------
if opts.useDefaultField:
    # Consume the unlabelled field MagneticField_cff already put in the
    # EventSetup: VolumeBasedMagneticField 160812 with
    # useParametrizedTrackerField=True, i.e. OAE_1103l_071212 inside the
    # tracker -- the field the SIM propagated through. Nothing to
    # instantiate: the empty label IS the default producer's label, and the
    # CPEs keep their default (empty) label so the Lorentz drift uses the
    # same field as everything else.
    fieldlabel = ""
elif opts.useOpera3D:
    from MagneticField.Engine.volumeBasedMagneticField_160812_cfi import \
        VolumeBasedMagneticFieldESProducer as Opera3DMagneticFieldProducer
    from MagneticField.Engine.volumeBasedMagneticField_160812_cfi import \
        magfield as MagneticFieldGeometry
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
        "useScalarPot3D=False is not supported; the legacy non-thread-safe "
        "wrapper class is not part of this port. Use the ScalarPot3D model.")
else:
    from MagneticField.ParametrizedEngine.parametrizedMagneticField_ScalarPot3D_cfi \
        import ParametrizedMagneticFieldProducer as ScalarPot3DMagneticFieldProducer
    process.ScalarPot3DMagneticFieldProducer = ScalarPot3DMagneticFieldProducer.clone()
    process.ScalarPot3DMagneticFieldProducer.parameters.InitFile = opts.scalarPot3DInitFile
    fieldlabel = "ScalarPot3DMf"
    process.ScalarPot3DMagneticFieldProducer.label = fieldlabel

process.geopro.MagneticFieldLabel = fieldlabel
process.Geant4ePropagator.MagneticFieldLabel = fieldlabel
# Activate the CVH-specific propagator path (custom fluct + its table pointer).
process.Geant4ePropagator.ForCVH = cms.bool(True)
process.Geant4ePropagator.PropagationDirection = cms.string(opts.propagationDirection)
process.Geant4ePropagator.PropagationPtotLimit = cms.double(float(opts.propagationPtotLimit))

# The Geant4 master is the shared EventSetup product (CvhMasterRecord),
# consumed by the maker via esConsumes; it builds its master G4 field with
# SimG4Core's FieldBuilder on top of the same labelled magnetic field the
# propagator consumes. `Particles` is left at the cfi default (the full
# canonical CVH set): the ES product is shared across every maker in a job,
# and the narrowing the J/psi driver applies to its own private CvhMaster PSet
# would not be safe to bake in here.
from TrackPropagation.Geant4e.cvhMasterESProducer_cfi import cvhMasterESProducer
process.cvhMasterESProducer = cvhMasterESProducer.clone()
process.cvhMasterESProducer.MagneticFieldLabel = cms.string(fieldlabel)
if opts.tightG4eStepper:
    _tighten(process.cvhMasterESProducer.MagneticField)   # the one that acts
    _tighten(process.geopro.MagneticField)                # kept in step
print("[cvh] effective: tightG4eStepper=%s -> cvhMasterESProducer.MagneticField"
      ".ConfGlobalMFM.OCMS.StepperParam %s"
      % (bool(opts.tightG4eStepper),
         {k: getattr(process.cvhMasterESProducer.MagneticField.ConfGlobalMFM
                     .OCMS.StepperParam, k).value()
          for k in sorted(_TIGHT_STEPPER)}))

# --- MiniAOD muon tracker tracks -> reco::TrackCollection (+ pat::Muon assoc)
process.tracksfrommuons = cms.EDProducer(
    "TrackProducerFromPatMuons",
    src=cms.InputTag(opts.muonSrc),
    innerTrackOnly=cms.bool(bool(opts.innerTrackOnly)),
    ptMin=cms.double(float(opts.muonPtMin)),
)

from Analysis.HitAnalyzer.diMuonTrackVertexCandidates_cfi import diMuonTrackVertexCandidates
from Analysis.HitAnalyzer.ResidualGlobalCorrectionMakerDiMuonG4e_cfi import \
    ResidualGlobalCorrectionMakerDiMuonG4e

process.diMuonTrackVertexCandidates = diMuonTrackVertexCandidates.clone(
    src="tracksfrommuons",
    massMin=float(opts.massMin),
    massMax=float(opts.massMax),
)

# Gauss-Newton momentum floor for the refit step clamp, derived from the
# propagation limit unless given explicitly. The pair is echoed by the maker
# itself ("[cvh] effective: clampMomentumFloor=...") because the two are only
# correct together.
_clampFloor = (float(opts.clampMomentumFloor) if float(opts.clampMomentumFloor) > 0.
               else 1.25 * float(opts.propagationPtotLimit))
if _clampFloor <= float(opts.propagationPtotLimit):
    raise RuntimeError(
        "clampMomentumFloor (%g GeV) must be ABOVE propagationPtotLimit (%g GeV): "
        "the Gauss-Newton clamp exists to keep the state out of the propagator's "
        "refusal region." % (_clampFloor, float(opts.propagationPtotLimit)))

process.trackrefitdimuon = ResidualGlobalCorrectionMakerDiMuonG4e.clone(
    src=cms.InputTag("tracksfrommuons"),
    srcCandidates=cms.InputTag("diMuonTrackVertexCandidates"),
    MagneticFieldLabel=cms.string(fieldlabel),
    scalarPotentialInitFile=cms.string(opts.scalarPot3DInitFile),
    materialGroupsFile=cms.string(opts.materialGroupsFile),
    globalMaterialModel=cms.bool(bool(opts.globalMaterialModel)),
    perStepFieldModes=cms.bool(bool(opts.perStepFieldModes)),
    skipHitlessSurfaces=cms.bool(bool(opts.skipHitlessSurfaces)
                                 and bool(opts.globalMaterialModel)),
    # --- outputs
    fillTrackTree=cms.bool(True),
    # fillRunTree carries the GLOBAL PARAMETER CATALOGUE (`runtree`); without
    # it the per-candidate globalidxv indices cannot be interpreted at all.
    fillRunTree=cms.bool(True),
    fillJac=cms.bool(bool(opts.fillJac)),
    fillGrads=cms.bool(bool(opts.fillGrads)),
    fillGradsFactored=cms.untracked.bool(bool(opts.fillGradsFactored)),
    doRes=cms.bool(bool(opts.doRes)),
    exportStepRecords=cms.bool(bool(opts.exportStepRecords)),
    exportCfExponents=cms.bool(bool(opts.exportCfExponents)),
    exportCfGroupExponents=cms.bool(bool(opts.exportCfGroupExponents)),
    exportHitResBlocks=cms.bool(bool(opts.exportHitResBlocks)),
    exportMaterialNoise=cms.bool(bool(opts.exportMaterialNoise)),
    exportVarianceGrads=cms.bool(bool(opts.exportVarianceGrads)),
    varianceGradFamilies=cms.vuint32(*[int(x) for x in opts.varianceGradFamilies]),
    exportObjective=cms.bool(bool(opts.exportObjective)),
    varianceFDGlobalIdx=cms.int32(int(opts.varianceFDGlobalIdx)),
    varianceFDEps=cms.double(float(opts.varianceFDEps)),
    produceValueMaps=cms.bool(bool(opts.produceValueMaps)),
    # --- gen (MiniAOD retags)
    doGen=cms.bool(bool(opts.doGen)),
    genParticles=cms.InputTag(opts.genParticles),
    pileupInfo=cms.InputTag(opts.pileupInfo),
    requireGen=cms.bool(bool(opts.requireGen)),
    fitFromGenParms=cms.bool(bool(opts.fitFromGenParms)),
    doSim=cms.bool(bool(opts.doSimHits)),
    # --- triggers
    doTrigger=cms.bool(bool(opts.doTrigger)),
    triggers=cms.vstring(*Z_TRIGGERS),
    # --- the fit
    useIdealGeometry=cms.bool(bool(opts.useIdealGeometry)),
    doVtxConstraint=cms.bool(bool(opts.doVtxConstraint)),
    minNdof=cms.int32(int(opts.minNdof)),
    minPairHits=cms.int32(int(opts.minPairHits)),
    minLegHits=cms.int32(int(opts.minLegHits)),
    exportVtxResidual=cms.bool(bool(opts.exportVtxResidual)),
    doMassConstraint=cms.bool(bool(opts.doMassConstraint)),
    massConstraint=cms.double(float(opts.massConstraint)),
    massConstraintWidth=cms.double(float(opts.massConstraintWidth)),
    keepPixelEdgeHits=cms.bool(bool(opts.keepPixelEdgeHits)),
    pixelMinSizeX=cms.int32(int(opts.pixelMinSizeX)),
    nIters=cms.uint32(int(opts.nIters)),
    edmConvergence=cms.double(float(opts.edmConvergence)),
    outprefix=cms.untracked.string(opts.outprefix),
)
# Step damping. Set after the clone so a command-line switch always wins.
process.trackrefitdimuon.clampMomentumFloor = cms.double(_clampFloor)
process.trackrefitdimuon.maxMomentumStepFactor = cms.double(float(opts.maxMomentumStepFactor))
process.trackrefitdimuon.stepBacktracking = cms.bool(bool(opts.stepBacktracking))
process.trackrefitdimuon.maxChi2Backtrack = cms.uint32(int(opts.maxChi2Backtrack))
process.trackrefitdimuon.stepBacktrackFromIter = cms.uint32(int(opts.stepBacktrackFromIter))
process.trackrefitdimuon.armijoC = cms.double(float(opts.armijoC))
process.trackrefitdimuon.armijoSlack = cms.double(float(opts.armijoSlack))
print("[cvh] effective: PropagationPtotLimit=%g GeV, clampMomentumFloor=%g GeV"
      % (float(opts.propagationPtotLimit), _clampFloor))

# Per-stream CLHEP engine for the residual maker (the MT path wires it into
# Geant4's thread-local RNG at the top of every produce()).
process.RandomNumberGeneratorService.trackrefitdimuon = cms.PSet(
    initialSeed=cms.untracked.uint32(423456789),
    engineName=cms.untracked.string('HepJamesRandom'),
)

# After every explicit propagator assignment above, so a command-line switch
# wins over the driver's own defaults and the effective state is echoed once.
cvhSwitches.apply(process, opts)

process.hltFilter = cms.EDFilter(
    "HLTHighLevel",
    HLTPaths=cms.vstring(*[t + "_v*" for t in Z_TRIGGERS]),
    eventSetupPathsKey=cms.string(""),
    andOr=cms.bool(True),
    throw=cms.bool(False),
    TriggerResultsTag=cms.InputTag("TriggerResults", "", "HLT"),
)

# NO BeamSpotProducer: MiniAOD already carries `offlineBeamSpot`, which is the
# hard-coded tag the maker consumes; producing a second one in this process
# would shadow it.
_seq = (process.tracksfrommuons *
        process.diMuonTrackVertexCandidates *
        process.trackrefitdimuon)
process.reconstruction_step = cms.Path((process.hltFilter * _seq)
                                       if opts.applyHltFilter else _seq)
process.schedule = cms.Schedule(process.reconstruction_step)

print("[cvh] effective: field=%s(label=%r) GT=%s geometry=2016 DDD "
      "massWindow=[%g,%g] doRes=%s exportCfExponents=%s exportStepRecords=%s "
      "fillJac=%s fillGrads=%s fillGradsFactored=%s doGen=%s requireGen=%s "
      "fitFromGenParms=%s doSim=%s globalMaterialModel=%s(%s) "
      "scalarPot3DInitFile=%s nThreads=%d skipEvents=%d nEvents=%d"
      % ("default" if opts.useDefaultField else
         ("opera3d" if opts.useOpera3D else "scalarpot3d"), fieldlabel,
         opts.globalTag, float(opts.massMin), float(opts.massMax),
         bool(opts.doRes), bool(opts.exportCfExponents),
         bool(opts.exportStepRecords), bool(opts.fillJac), bool(opts.fillGrads),
         bool(opts.fillGradsFactored), bool(opts.doGen), bool(opts.requireGen),
         bool(opts.fitFromGenParms), bool(opts.doSimHits),
         bool(opts.globalMaterialModel),
         os.path.basename(opts.materialGroupsFile) or "off",
         os.path.basename(opts.scalarPot3DInitFile),
         int(opts.numberOfThreads), int(opts.skipEvents), int(opts.nEvents)))
