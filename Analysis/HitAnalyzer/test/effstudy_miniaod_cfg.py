## eb96caef efficiency A/B on the WMass SingleMuon data MiniAOD.
##
## Thin wrapper over the production data nano config so the CVH refit runs
## EXACTLY as in production (muon chain -> tracksfrommuons -> trackrefit,
## GlobalTag 106X_dataRun2_v35, real data alignment via GeometryRecoDB). We:
##   * flip trackrefit.fillTrackTree on (production leaves it off -> ValueMaps),
##   * set the A/B toggle trackrefit.gluedGarbageTiltThreshold,
##   * (MINIMAL) read tracks straight off slimmedMuons and run only
##     geopro -> tracksfrommuons -> trackrefit, skipping the whole nano
##     (the WMass MiniAOD embeds the muonBestTrack recHits),
##   * drop the NANOAOD writer (the maker writes <outprefix>_<stream>.root).
##
## Parameters via cmsRun VarParsing (for the slurm harness) OR env vars:
##   input / EFF_FILES    input file (root:// URL or /store LFN); repeatable via comma
##   gluedTiltThr/EFF_THR  0.05 = fix ON, 1e9 = fix OFF (buggy)
##   nEvents / EFF_NEV     events (-1 = all)
##   minimal / EFF_MINIMAL 1 = fast path (default), 0 = full nanoSequence
import os
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

opts = VarParsing.VarParsing('analysis')
opts.register('input', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'input file (root:// or /store LFN)')
opts.register('gluedTiltThr', float(os.environ.get('EFF_THR', '0.05')),
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float, '0.05=fix, 1e9=bug')
opts.register('nEvents', int(os.environ.get('EFF_NEV', '50')),
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, 'events (-1=all)')
opts.register('minimal', int(os.environ.get('EFF_MINIMAL', '1')),
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, '1=fast path, 0=full nano')
opts.parseArguments()

_PROD = ("/work/submit/david_w/ZMass/CMSSW_10_6_26/src/Configuration/"
         "WMassNanoProduction/configs/NanoV9DataPostVFP_cfg.py")
exec(open(_PROD).read())  # defines `process`

THR = float(opts.gluedTiltThr)
TAG = "fix" if THR < 1.0 else "bug"

# input files: VarParsing 'input=' (single) or env EFF_FILES (comma list) else keep
FILES = opts.input or os.environ.get("EFF_FILES", "").strip()
if FILES:
    def _url(f):
        f = f.strip()
        return f if f.startswith(("root://", "file:", "/store/")) else "file:" + f
    process.source.fileNames = cms.untracked.vstring(
        *[_url(f) for f in FILES.split(",") if f.strip()])

process.maxEvents.input = cms.untracked.int32(int(opts.nEvents))

process.trackrefit.fillTrackTree = cms.bool(True)
process.trackrefit.gluedGarbageTiltThreshold = cms.untracked.double(THR)
process.trackrefit.outprefix = cms.untracked.string("effstudy_miniaod_" + TAG)

if int(opts.minimal):
    process.tracksfrommuons.src = cms.InputTag("slimmedMuons")
    process.effpath = cms.Path(process.geopro + process.tracksfrommuons + process.trackrefit)
    process.schedule = cms.Schedule(process.effpath)
    print("[effstudy] MINIMAL path: slimmedMuons -> tracksfrommuons -> trackrefit")
else:
    process.schedule = cms.Schedule(process.nanoAOD_step)

print("[effstudy] THR=%g TAG=%s nEvents=%s minimal=%s nFiles=%d outprefix=%s"
      % (THR, TAG, opts.nEvents, opts.minimal, len(process.source.fileNames),
         process.trackrefit.outprefix.value()))
