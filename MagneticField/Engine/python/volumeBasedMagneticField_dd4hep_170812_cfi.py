"""
MF map, model version 170812, MF geometry read from xml with DD4hep.

Backport to CMSSW_15_0 of the upstream cfi (cms-sw/cmssw#51526, 20_1_X). The
170812 model differs from 160812 only in the corrected ss400 steel
magnetization curves; the MF geometry is unchanged (geometryVersion 160812).

The 15_0 DD4hep builder reads the MERGED interpolation tables
MagneticField/Interpolation/data/grid_170812_3_8t/merged.{bin,index}
(cms-data/MagneticField-Interpolation#5; not in the 15_0 data package), which
must therefore be on CMSSW_SEARCH_PATH, e.g. in
$CMSSW_BASE/external/$SCRAM_ARCH/data/MagneticField/Interpolation/data/grid_170812_3_8t/
(shipped by CRAB with config.JobType.sendExternalFolder = True).
The DDD builder of 15_0 (volumeBasedMagneticField_160812_cfi style) cannot read
merged tables and would need the 11136 individual table files instead.
"""

import FWCore.ParameterSet.Config as cms

from MagneticField.Engine.volumeBasedMagneticField_160812_cfi import VBFConfig_160812
from MagneticField.Engine.volumeBasedMagneticField_160812_cfi import ParametrizedMagneticFieldProducer

DDDetectorESProducer = cms.ESSource("DDDetectorESProducer",
    confGeomXMLFiles = cms.FileInPath('MagneticField/GeomBuilder/data/cms-mf-geometry_160812.xml'),
    rootDDName = cms.string('cmsMagneticField:MAGF'),
    appendToDataLabel = cms.string('magfield')
)

DDCompactViewMFESProducer = cms.ESProducer("DDCompactViewMFESProducer",
    appendToDataLabel = cms.string('magfield')
)

VBFConfig_170812 = VBFConfig_160812.clone(
    version = cms.string('grid_170812_3_8t'),
    geometryVersion = cms.int32(160812),
)

VolumeBasedMagneticFieldESProducer = cms.ESProducer("DD4hep_VolumeBasedMagneticFieldESProducer",
    VBFConfig_170812,
    useMergeFileIfAvailable = cms.bool(True),
    appendToDataLabel = cms.string(''),
)
