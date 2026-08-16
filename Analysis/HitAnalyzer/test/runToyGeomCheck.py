# Milestone 1 of the toy-geometry clean-propagation test: does a minimal
# homogeneous DDD geometry load, and will Geant4 transport a muon through it?
#
# This deliberately does NOT produce hits yet. The plan (see NOTES 2026-08-12)
# is to cut the test loose from tracker DetIds entirely:
#
#   sim   : a stepping action records the true 5D state at chosen planes,
#           replacing PSimHits -- which is what currently forces the real
#           tracker geometry, since PSimHits only exist where the tracker
#           sensitive detector runs.
#   model : G4ePropagationExport gains a path taking explicit planes instead of
#           target DetIds (it currently throws if a DetId is not in the
#           TrackerGeometry).
#
# Both are additive; nothing existing changes behaviour. But neither is worth
# writing until we know CMSSW's Geant4 will run on a non-CMS geometry at all,
# which is the single biggest unknown in the whole plan. Hence this check.
#
# usage:
#   cmsRun runToyGeomCheck.py [events=10] [pt=3] [eta=0.30]

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.Eras.Era_Run2_2016_cff import Run2_2016

opts = VarParsing('analysis')
opts.register('events', 10, VarParsing.multiplicity.singleton,
              VarParsing.varType.int, 'number of muons to transport')
opts.register('pt', 3.0, VarParsing.multiplicity.singleton,
              VarParsing.varType.float, 'muon pT [GeV]')
opts.register('eta', 0.30, VarParsing.multiplicity.singleton,
              VarParsing.varType.float, 'muon eta')
opts.register('output', 'toystates.root', VarParsing.multiplicity.singleton,
              VarParsing.varType.string, 'output file')
opts.register('loosestepper', False, VarParsing.multiplicity.singleton,
              VarParsing.varType.bool,
              'use the CMSSW default (loose) stepper instead of the tight one')
opts.register('maxstep', -1.0, VarParsing.multiplicity.singleton,
              VarParsing.varType.float,
              'G4 field-propagation max step [cm]; <0 keeps the default 150')
opts.register('seed', 0, VarParsing.multiplicity.singleton,
              VarParsing.varType.int,
              'Geant4 random seed offset. 0 (default) keeps the historical '
              'seeding EXACTLY, so an unseeded rerun still reproduces the '
              'existing hs*.root samples bit-for-bit. Any other value gives '
              'an independent stream, which is what makes the sim splittable '
              'across jobs: the initial state is fixed and only the G4 seed '
              'varies between events, so N jobs at different seeds are N '
              'chunks of the same sample and the offline side just globs them.')
opts.parseArguments()

process = cms.Process('TOYGEOM', Run2_2016)

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
# The field map is delivered through conditions, so without a GlobalTag the
# IdealMagneticFieldRecord has no valid IOV. Same tag the real clean-propagation
# drivers pin, so the toy and the existing test see an identical field.
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mcRun2_asymptotic_v17', '')

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(opts.events))
process.source = cms.Source('EmptySource')

# ---------------------------------------------------------------- geometry
# The toy replaces the CMS geometry entirely. Only the pieces the toy actually
# needs are loaded: the standard material definitions (so Silicon/Air resolve),
# rotations, the world extent, and the toy volume itself. NOTE this deliberately
# does NOT load GeometryRecoDB/GeometrySimDB or the DB XML tag the real
# cleanprop drivers use -- there is no tracker here at all.
# A MINIMAL standalone geometry does not work: g4SimHits pulls
# CaloSimulationConstants, which is produced from HcalSimNumberingRecord, which
# needs the HCal geometry -- and OnlySDs does not gate that. Satisfying CMSSW's
# calo/muon EventSetup web is the expensive part, not the volumes.
#
# So keep the FULL standard 2016 list and substitute only the tracker VOLUME
# xmls with the toy. Every record stays satisfied; the tracking region -- the
# only part this test touches -- becomes homogeneous.
from Geometry.CMSCommonData.cmsExtendedGeometry2016aXML_cfi import \
    XMLIdealGeometryESSource as _stdGeom

_DROP = ('Geometry/TrackerCommonData/', 'Geometry/TrackerSimData/',
         'Geometry/TrackerRecoData/',
         # PLT/BCM sit just outside the pixel forward and position themselves
         # off pixfwd:ZPixelForward / pixfwd:RootStartZ, so they cannot survive
         # the tracker removal. They are luminosity detectors at |z| ~ 1.8 m,
         # irrelevant to a track that stops inside the tracking volume.
         'Geometry/ForwardCommonData/data/pltbcm.xml',
         'Geometry/ForwardCommonData/data/bcm1f.xml',
         'Geometry/ForwardSimData/data/bcm1fsens.xml',
         'Geometry/ForwardCommonData/data/plt.xml')
# trackermaterial.xml must survive the cull even though it lives under
# TrackerCommonData: it defines MATERIALS ONLY, no volumes, and
# Geometry/EcalCommonData/data/ectkcable.xml (ECAL, kept) references them.
# Without it DDG4Builder::convertMaterial throws "material is not valid from
# the Detector Description" on the dangling reference, on a G4 worker thread,
# which aborts the job with no readable message.
_KEEP = ('Geometry/TrackerCommonData/data/trackermaterial.xml',)
_files = [f for f in _stdGeom.geomXMLFiles
          if f in _KEEP or not f.startswith(_DROP)]
# the toy defines tracker:Tracker, which cmsTracker.xml (kept, it lives in
# CMSCommonData) positions into cms:CMSE. It must come after cms.xml, which
# defines the TrackBeam*/TrackCalorR constants the envelope references.
_files.append('Analysis/HitAnalyzer/data/tracker.xml')

# Load the STANDARD geometry sequence, then override only its file list. It is
# not enough to supply an XMLIdealGeometryESSource by hand: that provides the
# geometry but none of the ESProducers derived from it, so g4SimHits still fails
# on HcalSimNumberingRecord. This brings the ecal/hcal/muon constants producers
# along. The tracker numbering producer comes too but is never invoked, since
# ESProducers are lazy and nothing here consumes a tracker record.
_nstd = len(_stdGeom.geomXMLFiles)
process.load('Configuration.Geometry.GeometryExtended2016_cff')
process.XMLIdealGeometryESSource.geomXMLFiles = cms.vstring(*_files)

# Drop the tracker numbering producer. g4SimHits consumes only DDCompactView,
# but TrackerGeometricDetESModule PRODUCES INTO THE SAME IdealGeometryRecord, so
# making that record valid invokes it -- and with no tracker volumes it throws
# "The first child of the DDFilteredView is not what is expected". Nothing in
# this job consumes GeometricDet, so removing the producer is safe.
del process.trackerNumberingGeometry
print('[toy] geometry: %d of %d standard XML files kept, tracker volumes '
      'replaced by the homogeneous envelope'
      % (len(_files) - 1, _nstd))

# ---------------------------------------------------------------- generator
# fixed initial state, exactly as the real clean-propagation sim does: the only
# thing varying between events is the Geant4 seed, so the spread IS the
# propagation kernel
process.generator = cms.EDProducer(
    'FlatRandomPtGunProducer',
    PGunParameters=cms.PSet(
        PartID=cms.vint32(13),
        MinPt=cms.double(opts.pt), MaxPt=cms.double(opts.pt),
        MinEta=cms.double(opts.eta), MaxEta=cms.double(opts.eta),
        MinPhi=cms.double(0.70), MaxPhi=cms.double(0.70),
    ),
    AddAntiParticle=cms.bool(False),
    Verbosity=cms.untracked.int32(0),
    firstRun=cms.untracked.uint32(1),
)

# ---------------------------------------------------------------- Geant4
# Sensitive detectors. Do NOT empty OnlySDs: sensitiveDetectorMakers.cc does
#     if (chosenMakers.empty()) { ...create every registered maker... }
# so an empty list means "all", not "none" -- which is how a Hcal TEST BEAM
# detector appeared and demanded HcalTB06BeamParameters. Instead filter the
# default list, dropping the SDs whose volumes no longer exist:
#   TkAccumulating -> the tracker is now the homogeneous envelope, and it would
#                     also need the GeometricDet we deleted above
#   PLT / BCM1F    -> dropped from the geometry (they position off pixfwd)
# The calo and muon SDs are kept: those volumes are still there and their
# records come from the standard geometry sequence. TrackHits/CaloHits keep
# their defaults so the surviving SDs' produces() declarations still match.
_dropSD = ('TkAccumulatingSensitiveDetector', 'PLTSensitiveDetector',
           'BCM1FSensitiveDetector')
process.g4SimHits.OnlySDs = cms.vstring(
    *[s for s in process.g4SimHits.OnlySDs if s not in _dropSD])
# SURFACE-INTERSECTION PRECISION. runCleanPropSim.py tightens the Geant4
# stepper 100x over the CMSSW defaults, and the toy MUST match it or the two
# tests are not measuring the same truth. Defaults are DeltaIntersection =
# 1e-4 mm = 0.1 um, which is 2.4 % of sigma at the innermost plane (4.2 um) --
# not negligible, and largest exactly where the homogeneous toy was worst.
# Tightening the track-surface intersection in simulation is itself one of the
# documented CVH ingredients (AN-21-131), so it is not an optional detail.
_sp = process.g4SimHits.MagneticField.ConfGlobalMFM.OCMS.StepperParam
if opts.loosestepper:
    _sp.DeltaOneStepTracker = 1e-4
    _sp.DeltaIntersectionTracker = 1e-6
    _sp.DeltaOneStep = 1e-3
    _sp.DeltaIntersection = 1e-4
    print('[toy] LOOSE stepper (CMSSW defaults)')
else:
    _sp.DeltaOneStepTracker = 1e-5
    _sp.DeltaIntersectionTracker = 1e-6
    _sp.DeltaOneStep = 1e-5
    _sp.DeltaIntersection = 1e-6
    print('[toy] TIGHT stepper, matching runCleanPropSim.py')

process.g4SimHits.UseMagneticField = cms.bool(True)
# Step-length knob under test. Default is 150 cm, i.e. never limiting inside a
# 1 mm layer. Any interpretation of a closure change REQUIRES first showing the
# knob does something: StepperParam in this package was once found completely
# inert while looking like a valid control, so an unchanged result is ambiguous
# between "no step dependence" and "knob does nothing".
if opts.maxstep > 0:
    process.g4SimHits.MagneticField.ConfGlobalMFM.OCMS.MaxStep = cms.double(opts.maxstep)
    print('[toy] G4 MaxStep set to %g cm' % opts.maxstep)
process.g4SimHits.Physics.DefaultCutValue = cms.double(1.0)

# PassThroughEvtVtxGenerator, not a real VtxSmeared: the whole point of the
# clean-propagation test is ONE fixed initial state, with only the Geant4 seed
# varying between events. Standard vertex smearing would reintroduce an
# event-by-event initial state. Same construction as runCleanPropSim.py.
process.VtxSmeared = cms.EDProducer(
    'PassThroughEvtVtxGenerator',
    src=cms.InputTag('generator', 'unsmeared'),
    readDB=cms.bool(False),
)
process.RandomNumberGeneratorService.generator.initialSeed = 1
process.RandomNumberGeneratorService.VtxSmeared.initialSeed = 1
# g4SimHits carries the ONLY stream that matters here (generator and VtxSmeared
# emit a fixed state every event), so splitting the sample across jobs means
# giving each job its own g4SimHits seed and nothing else. seed=0 leaves the
# service default untouched, which is what every existing sample used.
if opts.seed:
    process.RandomNumberGeneratorService.g4SimHits.initialSeed = opts.seed
    print('[toy] g4SimHits initialSeed = %d' % opts.seed)

# The watcher replaces PSimHits: it records the true 5D state wherever the
# primary ends a step ON a shell boundary (fGeomBoundary), which is exact.
# Radii must match the shell edges in data/tracker.xml.
process.g4SimHits.Watchers = cms.VPSet(cms.PSet(
    type=cms.string('ToyStateNtuplizer'),
    radii=cms.vdouble(4.2000, 7.3000, 10.2000, 26.9000, 35.2000, 40.3000, 43.2000, 51.8000, 61.6000, 67.9000, 75.0000, 85.5000, 94.6000, 106.8000),
    tolerance=cms.double(1e-4),
    output=cms.string(opts.output),
))
# one thread / one stream, as runCleanPropSim.py does: the watcher owns its
# output file, and Geant4 MT would give one watcher instance per worker
process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(1),
    numberOfStreams=cms.untracked.uint32(1),
)

process.p = cms.Path(
    process.generator
    * process.VtxSmeared
    * process.generatorSmeared
    * process.g4SimHits
)

process.MessageLogger.cerr.FwkReport.reportEvery = 1
