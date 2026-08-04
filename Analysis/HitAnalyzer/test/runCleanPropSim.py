# runCleanPropSim.py -- ground-truth side of the clean propagation test.
#
# Simulates the SAME single particle, with a COMPLETELY FIXED initial state,
# N times through the full CMS Geant4 simulation, and writes the true 5D local
# state where it enters each silicon sensor (SimHitStateNtuplizer). The only
# source of event-to-event variation is the Geant4 random engine -- there is no
# vertex smearing, no kinematic smearing, no pileup, no FSR, no reconstruction,
# no selection. The spread of the output IS the propagation kernel.
#
#   cmsRun runCleanPropSim.py nEvents=20000 pt=10 eta=0 phi=0.35 seed=1 \
#          output=simstates.root
#
# Conditions deliberately match the rung-D' private sim recipe (design GT,
# ideal geometry, DB grid field, era Run2_2016) so that the model side
# (runCleanPropModel.py) propagates through identical geometry and field.

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
from Configuration.AlCa.GlobalTag import GlobalTag

opts = VarParsing.VarParsing('analysis')
opts.register('nEvents', 1000, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, 'number of simulated events')
opts.register('pt', 10.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float, 'fixed transverse momentum [GeV]')
opts.register('eta', 0.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float, 'fixed pseudorapidity')
opts.register('phi', 0.35, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float, 'fixed azimuth [rad]')
opts.register('partId', 13, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, 'PDG id of the gun particle (13 = mu-)')
opts.register('seed', 1, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, 'Geant4 random seed (the ONLY thing to vary across tasks)')
opts.register('output', 'simstates.root', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'output ntuple')
opts.parseArguments()

process = cms.Process('CLEANSIM', Run2_2016)

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Conditions/geometry pinned to exactly what the CVH refit drivers use
# (runCvhResClosure.py), so that the model side propagates through the same
# material. NOTE: 'auto:run2_design' + the 'Ideal' XML label does NOT work here
# -- it yields a geometry with no pixel sensitive volumes at all (verified:
# TrackerHitsPixelBarrelLowTof is empty for a muon that must cross BPIX), so
# the muon would only be seen from TIB outwards.
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mcRun2_asymptotic_v17', '')
process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(
        record=cms.string('GeometryFileRcd'),
        tag=cms.string('XMLFILE_Geometry_2016_81YV1_Extended2016_mc'),
        label=cms.untracked.string('Extended'),
    ),
)
process.XMLFromDBSource.label = cms.string('Extended')

process.source = cms.Source('EmptySource')
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(int(opts.nEvents)))
process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(1),
    numberOfStreams=cms.untracked.uint32(1),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# --- the fixed-state gun -----------------------------------------------------
# FlatRandomPtGunProducer shoots RandFlat(min, max); with min == max that is
# exactly the endpoint, bit for bit, and consumes the same number of randoms
# every event. The production vertex is hard-coded at the origin.
process.generator = cms.EDProducer(
    'FlatRandomPtGunProducer',
    PGunParameters=cms.PSet(
        PartID=cms.vint32(int(opts.partId)),
        MinPt=cms.double(float(opts.pt)), MaxPt=cms.double(float(opts.pt)),
        MinEta=cms.double(float(opts.eta)), MaxEta=cms.double(float(opts.eta)),
        MinPhi=cms.double(float(opts.phi)), MaxPhi=cms.double(float(opts.phi)),
    ),
    Verbosity=cms.untracked.int32(0),
    AddAntiParticle=cms.bool(False),
    firstRun=cms.untracked.uint32(1),
    psethack=cms.string('fixed-state single particle for the clean propagation test'),
)

# no vertex smearing at all (the standard VtxSmeared modules would reintroduce
# an event-by-event initial state, which is precisely what must not happen).
process.VtxSmeared = cms.EDProducer(
    'PassThroughEvtVtxGenerator',
    src=cms.InputTag('generator', 'unsmeared'),
    readDB=cms.bool(False),
)

# --- seeds -------------------------------------------------------------------
# Only Geant4 is allowed to vary: the generator and the vertex generator are
# deterministic by construction, but they still get distinct seeds so that
# nothing is silently shared.
process.RandomNumberGeneratorService.generator.initialSeed = 1
process.RandomNumberGeneratorService.VtxSmeared.initialSeed = 1
process.RandomNumberGeneratorService.g4SimHits.initialSeed = 100000 + int(opts.seed)

# --- ntuple ------------------------------------------------------------------
_SIMHIT_TAGS = [
    'TrackerHitsPixelBarrelLowTof', 'TrackerHitsPixelBarrelHighTof',
    'TrackerHitsPixelEndcapLowTof', 'TrackerHitsPixelEndcapHighTof',
    'TrackerHitsTIBLowTof', 'TrackerHitsTIBHighTof',
    'TrackerHitsTIDLowTof', 'TrackerHitsTIDHighTof',
    'TrackerHitsTOBLowTof', 'TrackerHitsTOBHighTof',
    'TrackerHitsTECLowTof', 'TrackerHitsTECHighTof',
]
process.simstates = cms.EDAnalyzer(
    'SimHitStateNtuplizer',
    simHitTags=cms.VInputTag(*[cms.InputTag('g4SimHits', t) for t in _SIMHIT_TAGS]),
    pdgId=cms.int32(int(opts.partId)),
    trackId=cms.uint32(1),
)
process.TFileService = cms.Service('TFileService', fileName=cms.string(opts.output))

process.p = cms.Path(
    process.generator
    * process.VtxSmeared
    * process.generatorSmeared
    * process.g4SimHits
    * process.simstates
)
process.schedule = cms.Schedule(process.p)
