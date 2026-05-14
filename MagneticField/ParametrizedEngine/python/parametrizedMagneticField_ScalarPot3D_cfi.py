"""ESSource + ESProducer config for the spherical-harmonic
scalar-potential B-field model .

Registers a labelled `IdealMagneticFieldRecord` so CVH-side consumers
can be redirected via `MagneticFieldLabel = 'ScalarPot3DMf'`.

Usage in nano_cff.py / runCvhJpsi.py:

    from MagneticField.ParametrizedEngine.parametrizedMagneticField_ScalarPot3D_cfi \\
        import ParametrizedMagneticFieldProducer as ScalarPot3DMagneticFieldProducer
    process.ScalarPot3DMagneticFieldProducer = ScalarPot3DMagneticFieldProducer
    process.ScalarPot3DMagneticFieldProducer.parameters.InitFile = '<path-to-dump.txt>'
    fieldlabel = 'ScalarPot3DMf'
    process.ScalarPot3DMagneticFieldProducer.label = fieldlabel
    # ... and the seven `process.<consumer>.MagneticFieldLabel = fieldlabel`
"""

import FWCore.ParameterSet.Config as cms

idealMagneticFieldRecordSource = cms.ESSource("EmptyESSource",
    recordName = cms.string('IdealMagneticFieldRecord'),
    iovIsRunNotTime = cms.bool(True),
    firstValid = cms.vuint32(1)
)

ParametrizedMagneticFieldProducer = cms.ESProducer("ParametrizedMagneticFieldProducer",
    version = cms.string('ScalarPot3D'),
    parameters = cms.PSet(
        # Path to a coefficient dump file produced by mfs/dump_coeffs_for_cmssw.py.
        # Required (cannot be empty).
        InitFile = cms.string(''),
        # Sphere radius (cm) outside which inTesla returns zero.
        # The default 320 cm matches the typical mfs fit domain
        # (sphere320: r<290, |z|<316, R<320). a future update will replace this
        # hard cutoff with a C^1 cosine blend to the underlying CMSSW
        # volume-based field.
        ValidityRadius = cms.double(320.0),
    ),
    label = cms.untracked.string('ScalarPot3DMf'),
)
