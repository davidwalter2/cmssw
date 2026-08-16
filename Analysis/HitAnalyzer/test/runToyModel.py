# Toy-geometry MODEL side of the clean-propagation test: the deterministic
# reference propagation through the HOMOGENEOUS volume, exporting the same
# per-step physics records and transport Jacobians as the real one.
#
# Differs from runCleanPropModel.py in exactly two ways:
#   * geometry: tracker volumes replaced by the homogeneous envelope, same
#     substitution as runToyGeomCheck.py
#   * targets: EXPLICIT PLANES instead of DetIds. The toy has no tracker and no
#     DetIds, and G4ePropagationExport throws on a DetId it cannot find in the
#     TrackerGeometry, so the surfaces are given directly. The planes are the
#     tangent planes at the reference crossings (toyPlanes_pt3.py), whose
#     azimuth was validated against the sim: predicted 0.9043 vs measured 0.906.
#
# usage: cmsRun runToyModel.py output=model_toy.root

import math
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
from Configuration.AlCa.GlobalTag import GlobalTag
import toyPlanes_pt3 as planes

opts = VarParsing('analysis')
opts.register('output', 'model_toy.root', VarParsing.multiplicity.singleton,
              VarParsing.varType.string, 'output file')
opts.register('pt', 3.0, VarParsing.multiplicity.singleton,
              VarParsing.varType.float, 'muon pT [GeV]')
opts.register('eta', 0.30, VarParsing.multiplicity.singleton,
              VarParsing.varType.float, 'muon eta')
opts.register('phi', 0.70, VarParsing.multiplicity.singleton,
              VarParsing.varType.float, 'muon phi')
opts.register('partId', 13, VarParsing.multiplicity.singleton,
              VarParsing.varType.int, 'PDG id (13 = mu-)')
opts.register('useIdealGeometry', True, VarParsing.multiplicity.singleton,
              VarParsing.varType.bool, 'use the ideal tracker geometry (must be True to match PSimHits)')
opts.register('useOpera3D', False, VarParsing.multiplicity.singleton,
              VarParsing.varType.bool, 'use the raw 160812 grid field instead of the default')
opts.register('ioniTruncationAlpha', 0.999, VarParsing.multiplicity.singleton,
              VarParsing.varType.float, 'ionization variance truncation alpha')
opts.register('stepLength', 10.0, VarParsing.multiplicity.singleton,
              VarParsing.varType.float,
              'Geant4e max step [mm]. Was hard-coded to 10.0 in '
              'Geant4ePropagator.cc; in a homogeneous medium the model sits '
              'exactly ON this ceiling (measured median step 1.0000 cm) while '
              'the SIM steps far more finely, and Moliere is not additive '
              'between the two.')
opts.parseArguments()

process = cms.Process('CLEANMODEL', Run2_2016)

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
from Geometry.CMSCommonData.cmsExtendedGeometry2016aXML_cfi import \
    XMLIdealGeometryESSource as _stdGeom
_DROP = ('Geometry/TrackerCommonData/', 'Geometry/TrackerSimData/',
         'Geometry/TrackerRecoData/',
         'Geometry/ForwardCommonData/data/pltbcm.xml',
         'Geometry/ForwardCommonData/data/bcm1f.xml',
         'Geometry/ForwardSimData/data/bcm1fsens.xml',
         'Geometry/ForwardCommonData/data/plt.xml')
_KEEP = ('Geometry/TrackerCommonData/data/trackermaterial.xml',)
_files = [f for f in _stdGeom.geomXMLFiles
          if f in _KEEP or not f.startswith(_DROP)]
_files.append('Analysis/HitAnalyzer/data/tracker.xml')
process.load('Configuration.Geometry.GeometryExtended2016_cff')
process.XMLIdealGeometryESSource.geomXMLFiles = cms.vstring(*_files)
del process.trackerNumberingGeometry
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# provides the "idealForDigi" labelled TrackerGeometry (the frame PSimHits use)

process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mcRun2_asymptotic_v17', '')

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
    targetDetIds=cms.vuint32(),
    targetLocalZ=cms.vdouble(),
    targetPlaneOrigin=cms.vdouble(*planes.origin),
    targetPlaneNormal=cms.vdouble(*planes.normal),
    targetPlaneU=cms.vdouble(*planes.uaxis),
    useIdealGeometry=cms.bool(bool(opts.useIdealGeometry)),
)
process.TFileService = cms.Service('TFileService', fileName=cms.string(opts.output))

process.p = cms.Path(process.propExport)
process.schedule = cms.Schedule(process.p)
