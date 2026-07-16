## Phase 0 of the global material model: geometry-only geantino-ray audit
## of the Geant4 tracker material (see doc/global-material-model-plan.md).
## Uses the same DB XML geometry as the CVH refit drivers; no field, no
## physics, no input data (EmptySource pinned to a 2016 run).
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
from Configuration.AlCa.GlobalTag import GlobalTag

opts = VarParsing.VarParsing('analysis')
opts.register('outFile', 'material_audit.txt', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'output tally file')
opts.register('nEta', 100, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, 'eta bins')
opts.register('nPhi', 72, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, 'phi bins')
opts.parseArguments()

process = cms.Process("AUDIT", Run2_2016)

process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.GeometrySimDB_cff")

process.GlobalTag = GlobalTag(process.GlobalTag, "auto:run2_data", "")
process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(
        record=cms.string("GeometryFileRcd"),
        tag=cms.string("XMLFILE_Geometry_2016_81YV1_Extended2016_mc"),
        label=cms.untracked.string("Extended"),
    ),
)
process.XMLFromDBSource.label = cms.string("Extended")

process.source = cms.Source("EmptySource",
                            firstRun=cms.untracked.uint32(278769))
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))

process.audit = cms.EDAnalyzer(
    "MaterialAuditAnalyzer",
    rMax=cms.double(125.),   # cm: past the outer TOB + support tube wall
    zMax=cms.double(300.),   # cm: past TEC9 + bulkhead
    etaMax=cms.double(2.5),
    nEta=cms.int32(int(opts.nEta)),
    nPhi=cms.int32(int(opts.nPhi)),
    vertexZ=cms.vdouble(-5., 0., 5.),  # cm along the beamline
    muonMomentum=cms.double(10.),      # GeV, for the dE/dx tally
    outFile=cms.string(opts.outFile),
)

process.p = cms.Path(process.audit)
