#

import FWCore.ParameterSet.Config as cms

process = cms.Process("MAGNETICFIELDTEST")

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1),
)

# Example configuration for the magnetic field

# Uncomment ONE of the following:

### Uniform field
#process.load("Configuration.StandardSequences.MagneticField_0T_cff")
#process.localUniform.ZFieldInTesla = 3.8


### Full field map, static configuration for each field value -> dynamic map?
#process.load("Configuration.StandardSequences.MagneticField_20T_cff")
#process.load("Configuration.StandardSequences.MagneticField_30T_cff")
#process.load("Configuration.StandardSequences.MagneticField_35T_cff")
# process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load("Configuration.StandardSequences.MagneticField_40T_cff")

# # instead of standard sequence: 
# # process.VolumeBasedMagneticFieldESProducer = cms.ESProducer("VolumeBasedMagneticFieldESProducerFromDB",
# process.VolumeBasedMagneticFieldESProducer = cms.ESProducer("VolumeBasedMagneticFieldESProducer",
#     label = cms.untracked.string(''),
#     debugBuilder = cms.untracked.bool(False),
#     valueOverride = cms.int32(-1), # Force value of current (in A); take the value from DB if < 0.
#     useParametrizedTrackerField = cms.untracked.bool(False),
#     version = cms.string('grid_160812_3_8t_Run1'),
#     geometryVersion = cms.int32(160812),
# )

# process.idealMagneticFieldRecordSource = cms.ESSource("EmptyESSource",
#     recordName = cms.string('IdealMagneticFieldRecord'),
#     iovIsRunNotTime = cms.bool(True),
#     firstValid = cms.vuint32(1)
# )

# # Parabolic parametrized magnetic field used for track building (scaled to nominal map closest to current from runInfo)
# process.ParabolicParametrizedMagneticFieldProducer = cms.ESProducer("AutoParametrizedMagneticFieldProducer",
#     version = cms.string('Parabolic'),
#     label = cms.untracked.string('ParabolicMf'),
#     valueOverride = cms.int32(-1)
# )

# process.VolumeBasedMagneticFieldESProducer.valueOverride = 18268
# process.ParabolicParametrizedMagneticFieldProducer.valueOverride = 18268


# -> static map?
# process.load("MagneticField.Engine.volumeBasedMagneticField_090322_2pi_scaled_cfi")

# process.load("MagneticField.Engine.volumeBasedMagneticField_1103l_cfi")

# process.load("MagneticField.Engine.volumeBasedMagneticField_120812_smallYE4_cfi")
# process.load("MagneticField.Engine.volumeBasedMagneticField_120812_largeYE4_cfi")

# process.load("MagneticField.Engine.volumeBasedMagneticField_130503_smallYE4_cfi")
# process.load("MagneticField.Engine.volumeBasedMagneticField_130503_largeYE4_cfi")

process.load("MagneticField.Engine.volumeBasedMagneticField_160812_cfi")
process.VolumeBasedMagneticFieldESProducer.version = "grid_160812_3_8t" # grid_160812_3_8t, grid_160812_3_8t_Run1

# process.load("MagneticField.Engine.volumeBasedMagneticField_170812_cfi")
# process.VolumeBasedMagneticFieldESProducer.version = "grid_170812_3_8t" # grid_170812_3_8t, grid_170812_3_8t_Run1, grid_170812_3_8t_SX5

# disable 2D parameterization in tracker and use slower but more accurate 3D splines of original data
process.VolumeBasedMagneticFieldESProducer.useParametrizedTrackerField = cms.bool(False) 

# activate MTCC field map data
# process.load("MagneticField.ParametrizedEngine.parametrizedMagneticField_PolyFit2D_cfi")
# process.load("MagneticField.ParametrizedEngine.parametrizedMagneticField_PolyFit3D_cfi")


### Configuration to select map based on recorded current in the DB
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')


#process.GlobalTag = GlobalTag(process.GlobalTag,'auto:phase1_2017_realistic', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
# process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
# process.load("MagneticField.Engine.autoMagneticFieldProducer_cfi")
# process.AutoMagneticFieldESProducer.valueOverride = 18000


### Set scaling factors
#process.VolumeBasedMagneticFieldESProducer.scalingVolumes = ( 802 , )
#process.VolumeBasedMagneticFieldESProducer.scalingFactors = ( 1.5 , )


# process.MessageLogger = cms.Service("MessageLogger",
#     categories   = cms.untracked.vstring("MagneticField"),
#     destinations = cms.untracked.vstring("cout"),
#     cout = cms.untracked.PSet(  
#     noLineBreaks = cms.untracked.bool(True),
#     threshold = cms.untracked.string("INFO"),
#     INFO = cms.untracked.PSet(
#       limit = cms.untracked.int32(0)
#     ),
#     WARNING = cms.untracked.PSet(
#       limit = cms.untracked.int32(0)
#     ),
#     MagneticField = cms.untracked.PSet(
#      limit = cms.untracked.int32(10000000)
#     )
#   )
# )

process.queryField  = cms.EDAnalyzer("queryField")
process.p1 = cms.Path(process.queryField)

