## Standalone unit test for JpsiKCandidateSplitter. Runs the splitter
## alone on one preset-B sl1 ALCARECO file and writes the four output
## collections to a small ROOT file so we can compare per-event counts.
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

opts = VarParsing.VarParsing('analysis')
opts.register('input', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'input ALCARECO file')
opts.register('nEvents', 200, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, 'events to process')
opts.parseArguments()

process = cms.Process("SPLIT")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

assert opts.input, "must set input=<path>"
_url = opts.input if opts.input.startswith(("root://", "file:")) else "file:" + opts.input
process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring(_url),
)
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(opts.nEvents))

from Analysis.HitAnalyzer.JpsiKCandidateSplitter_cfi import jpsiKCandidateSplitter
process.split = jpsiKCandidateSplitter.clone()

process.out = cms.OutputModule("PoolOutputModule",
    fileName=cms.untracked.string("splitter_unit_test.root"),
    outputCommands=cms.untracked.vstring(
        "drop *",
        # keep the new outputs
        "keep *_split_*_*",
        # also keep the original B+ collection for cross-check
        "keep *_ALCARECOTkAlJpsiX_*_*",
        "keep *_ALCARECOTkAlJpsiXBPlusResonances_*_*",
    ),
)

process.p = cms.Path(process.split)
process.e = cms.EndPath(process.out)
