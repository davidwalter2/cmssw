"""Minimal cmsRun smoke test for the ScalarPot3DMagneticField ESProducer.

Brings up the labelled field as the default `IdealMagneticFieldRecord`
(label="") and prints the field at a few hard-coded points via the
testMagneticField analyzer.

Run inside the el7 container:

    cmsRun MagneticField/ParametrizedEngine/test/testScalarPot3DField_cfg.py \\
        initFile=/work/submit/david_w/ZMass/mfs/data/fitresults/polyfit3d_full_coeffs_lmax18_cmsswnorm.txt
"""

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

opts = VarParsing('analysis')
opts.register('initFile', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
              'coefficient dump file path (required)')
opts.parseArguments()
if not opts.initFile:
    raise SystemExit("initFile=<path> is required")

process = cms.Process("ScalarPot3DTest")

process.MessageLogger = cms.Service("MessageLogger",
    cerr = cms.untracked.PSet(threshold = cms.untracked.string('INFO')),
    destinations = cms.untracked.vstring('cerr'),
)

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))

# Load our cfi. Override label="" so the standard `IdealMagneticFieldRecord`
# accessor (no label) returns our field.
process.load('MagneticField.ParametrizedEngine.parametrizedMagneticField_ScalarPot3D_cfi')

# Standalone job: nothing else provides the IdealMagneticFieldRecord IOV
# (the cfi deliberately ships no ESSource -- real jobs get it from the
# GlobalTag / standard field config).
process.idealMagneticFieldRecordSource = cms.ESSource("EmptyESSource",
    recordName = cms.string('IdealMagneticFieldRecord'),
    iovIsRunNotTime = cms.bool(True),
    firstValid = cms.vuint32(1)
)
process.ParametrizedMagneticFieldProducer.label = ''
process.ParametrizedMagneticFieldProducer.parameters.InitFile = opts.initFile

# Use queryField via a small input file (a few hand-picked points
# spanning inside-validity and outside-validity).
import os, tempfile
points = [
    (0.0,    0.0,   0.0),     # origin
    (50.0,   30.0,  100.0),   # inside tracker, canonical
    (76.207, -69.251, 12.435),# inside tracker, row 1
    (200.0,  100.0,  250.0),  # near tracker edge
    (350.0,  0.0,    0.0),    # outside validity sphere -> expect (0,0,0)
]
fd, ptfile = tempfile.mkstemp(prefix='scalarpot3d_pts_', suffix='.txt')
with os.fdopen(fd, 'w') as f:
    for x, y, z in points:
        f.write("{0} {1} {2}\n".format(x, y, z))

process.test = cms.EDAnalyzer("testMagneticField",
    outputTable = cms.untracked.string(''),
    inputTable  = cms.untracked.string(ptfile),
    inputTableType = cms.untracked.string('xyz'),
    InnerRadius = cms.untracked.double(0.0),
    OuterRadius = cms.untracked.double(150.0),
    HalfLength  = cms.untracked.double(280.0),
    numberOfPoints = cms.untracked.int32(5),
    resolution  = cms.untracked.double(1.0),  # loose; we just want the values printed
)

process.p = cms.Path(process.test)
