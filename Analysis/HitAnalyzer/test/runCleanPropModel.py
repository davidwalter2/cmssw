# runCleanPropModel.py -- model side of the clean propagation test.
#
# Propagates ONE fixed initial state deterministically with the CVH Geant4e
# propagator through the same sensor surfaces the simulation crosses, and
# exports the reference state, the leg transport Jacobians and noise matrices,
# and the per-step Urban/Moliere physics records plus per-step cumulative
# transport. Everything the offline characteristic-function model needs.
#
#   cmsRun runCleanPropModel.py pt=10 eta=0.3 phi=0.35 \
#          targets=targets.txt output=model.root
#
# `targets` is the file written by cf_propagation_test.py --targets from the
# simulation output: one line per crossed module, "<detid> <localZ>", ordered
# along the trajectory. localZ is the sensor entry face in the DetUnit frame.
#
# Conditions are pinned to the same values as runCleanPropSim.py. The tracker
# geometry is the IDEAL one by default because PSimHit local coordinates are.

import math
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
from Configuration.AlCa.GlobalTag import GlobalTag

opts = VarParsing.VarParsing('analysis')
opts.register('pt', 10.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float, 'fixed transverse momentum [GeV]')
opts.register('eta', 0.3, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float, 'fixed pseudorapidity')
opts.register('phi', 0.35, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float, 'fixed azimuth [rad]')
opts.register('partId', 13, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int, 'PDG id (13 = mu-)')
opts.register('targets', '', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'target list file: "<detid> <localZ>" per line')
opts.register('useIdealGeometry', True, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool, 'use the ideal tracker geometry (must be True to match PSimHits)')
opts.register('useOpera3D', False, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool, 'use the raw 160812 grid field instead of the default')
opts.register('ioniTruncationAlpha', 0.999, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float, 'ionization variance truncation alpha')
opts.register('output', 'model.root', VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string, 'output ntuple')
opts.register('stepLength', 10.0, VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.float,
              'Geant4e max step [mm]; formerly hard-coded in Geant4ePropagator.cc')
opts.parseArguments()

assert opts.targets, 'must pass targets=<file>'
_detids, _zoff = [], []
with open(opts.targets) as fh:
    for line in fh:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        a, b = line.split()
        _detids.append(int(a))
        _zoff.append(float(b))
assert _detids, 'empty target list'

process = cms.Process('CLEANMODEL', Run2_2016)

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# provides the "idealForDigi" labelled TrackerGeometry (the frame PSimHits use)
process.load('Geometry.TrackerGeometryBuilder.idealForDigiTrackerGeometry_cff')

process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mcRun2_asymptotic_v17', '')
process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(
        record=cms.string('GeometryFileRcd'),
        tag=cms.string('XMLFILE_Geometry_2016_81YV1_Extended2016_mc'),
        label=cms.untracked.string('Extended'),
    ),
)
process.XMLFromDBSource.label = cms.string('Extended')

process.load('TrackPropagation.Geant4e.geantRefit_cff')
from TrackPropagation.Geant4e.cvhMaster_cfi import CvhMasterPSet

process.source = cms.Source('EmptySource')
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))
process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(1),
    numberOfStreams=cms.untracked.uint32(1),
)

process.RandomNumberGeneratorService.propExport = cms.PSet(
    initialSeed=cms.untracked.uint32(123456789),
    engineName=cms.untracked.string('HepJamesRandom'),
)

fieldlabel = ""
if opts.useOpera3D:
    from MagneticField.Engine.volumeBasedMagneticField_160812_cfi import \
        VolumeBasedMagneticFieldESProducer as Opera3DMagneticFieldProducer
    from MagneticField.Engine.volumeBasedMagneticField_160812_cfi import magfield as MagneticFieldGeometry
    process.magfield = MagneticFieldGeometry
    process.es_prefer_magfield_cvhrefit = cms.ESPrefer("XMLIdealGeometryESSource", "magfield")
    process.Opera3DMagneticFieldProducer = Opera3DMagneticFieldProducer
    fieldlabel = "grid_160812_3_8t"
    process.Opera3DMagneticFieldProducer.label = fieldlabel
    process.Opera3DMagneticFieldProducer.useParametrizedTrackerField = cms.bool(False)

process.geopro.MagneticFieldLabel = fieldlabel
process.Geant4ePropagator.MagneticFieldLabel = fieldlabel
process.Geant4ePropagator.ForCVH = cms.bool(True)
process.Geant4ePropagator.PropagationDirection = "anyDirection"
process.Geant4ePropagator.IoniTruncationAlpha = cms.double(float(opts.ioniTruncationAlpha))
process.Geant4ePropagator.StepLengthLimit = cms.double(float(opts.stepLength))

# PDG id -> (G4 particle name, charge). NOTE the sign conventions differ by
# species class and getting this wrong is silent: leptons carry charge
# -sign(pdg) (mu- is +13), while mesons and baryons carry +sign(pdg)
# (pi+ is +211, K+ is +321). The old lepton-only rule gave K-/pi- a POSITIVE
# charge, i.e. a reference trajectory bending the wrong way.
_SPECIES = {
      13: ('mu-', -1.),      -13: ('mu+', +1.),
     211: ('pi+', +1.),     -211: ('pi-', -1.),
     321: ('kaon+', +1.),   -321: ('kaon-', -1.),
    2212: ('proton', +1.), -2212: ('anti_proton', -1.),
}
_pdg = int(opts.partId)
if _pdg not in _SPECIES:
    raise ValueError(f'partId={_pdg} not supported; known: {sorted(_SPECIES)}')
_g4name, _charge = _SPECIES[_pdg]

from TrackPropagation.Geant4e.cvhMasterESProducer_cfi import cvhMasterESProducer
process.cvhMasterESProducer = cvhMasterESProducer.clone()
process.cvhMasterESProducer.MagneticFieldLabel = cms.string(fieldlabel)
# the species pair on top of the always-on mandatory set (gamma/e+/e-/mu+/mu-/
# proton are registered regardless; see cvhMaster_cfi)
_pair = {13: ('mu+', 'mu-'), 211: ('pi+', 'pi-'),
         321: ('kaon+', 'kaon-'), 2212: ('proton', 'anti_proton')}[abs(_pdg)]
process.cvhMasterESProducer.Particles = cms.vstring(
    *dict.fromkeys(('gamma', 'e+', 'e-', 'mu+', 'mu-', 'proton') + _pair))

_theta = 2. * math.atan(math.exp(-float(opts.eta)))
_pt = float(opts.pt)
_pz = _pt / math.tan(_theta)
_px = _pt * math.cos(float(opts.phi))
_py = _pt * math.sin(float(opts.phi))

process.propExport = cms.EDAnalyzer(
    'G4ePropagationExport',
    initialPosition=cms.vdouble(0., 0., 0.),
    initialMomentum=cms.vdouble(_px, _py, _pz),
    charge=cms.double(_charge),
    particleName=cms.string(_g4name),
    targetDetIds=cms.vuint32(*_detids),
    targetLocalZ=cms.vdouble(*_zoff),
    targetPlaneOrigin=cms.vdouble(),
    targetPlaneNormal=cms.vdouble(),
    targetPlaneU=cms.vdouble(),
    useIdealGeometry=cms.bool(bool(opts.useIdealGeometry)),
)
process.TFileService = cms.Service('TFileService', fileName=cms.string(opts.output))

process.p = cms.Path(process.propExport)
process.schedule = cms.Schedule(process.p)
