# Archive of retired CMSSW_15_0_19_patch2_* test areas (removed 2026-11-03)

Files whose exact content was never committed anywhere, copied verbatim from the areas'
`src/` (paths below `archive/<area>/` mirror `src/`). Not for merging; dev has moved on.

| area | base | what it was | status |
| --- | --- | --- | --- |
| alignslope_lap | copy of CMSSW_15_0_19_patch2 @ nano15 5d3b689c16b (+ alignslope) | 10/9 Laplace / marginal-likelihood alignment-gradient term (`laplaceGrad`, `gradlapv`) fixing the localUpdate=True TEC z-expansion (704 -> 72 um) | **useful**: validated fix, not in dev; needed only if localUpdate=True returns |
| alignslope | same base | 10/7-10/8 alignSlopeMode 0-3, hitResPost/hitVinv/hitAlignJac debug branches, CVH_IONI_SCALE / CVH_FREEZE_CPE_ITER knobs | diagnostic, superseded by alignslope_lap |
| chordctl, chordval, chordver, chordver2 | chord-speedup development, 10/9 | control / validation / variant builds of the analytic layer chord | superseded by 1e2ad33 (PR #49) |
| cyclediag | 10/6 | limit-cycle diagnostics in the single-track maker | diagnostic |
| ebuild | 9/27 | private build of the nanoV15 merge; tracker.xml material test geometry, older resolution-cf-export.md | superseded |
| lucheck | 10/3 | localUpdate diagnostics (also as calibration_studies localupdate_check/lucheck_diagnostics.patch) | diagnostic |
| speed | branch speed @ 8fbd0f6af23 | 10/3 prototype of the radiative dE/dx table + within-volume relocation (uncommitted diff) | superseded by 4ac946b / 1e2ad33 (PR #49) |

Other backups: `backup/nano15-wip-260923` (stash + July test cfgs of CMSSW_15_0_19_patch2),
`backup/trklik-frozenmodel-260917` (two-track frozen-model export of the trklik area).
