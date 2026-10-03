"""OPERA/TOSCA volume-based magnetic field for the CVH refit (one place for the
NanoAOD customise and the standalone drivers).

The field is added as a LABELLED ES product with the full 3D interpolation
grid inside the tracker (useParametrizedTrackerField = False) and, optionally,
routed into the strip/pixel CPEs so the Lorentz drift of the hit
re-evaluation uses the same field. The caller routes the returned label into
the propagator, the refit makers and the shared G4 master.

version "170812" (default): the latest OPERA model (corrected ss400 steel BH
  curves; inside the tracker it agrees with 160812 to <= 1e-5 T). Built with
  the DD4hep builder, which in CMSSW_15_0 is the one that reads the MERGED
  tables grid_170812_3_8t/merged.{bin,index} (76 MB; byte-identical to
  cms-data/MagneticField-Interpolation#5, which is not in the 15_0 data
  package). The tables must be on CMSSW_SEARCH_PATH -- see TABLES_HELP.
version "160812": the release tables (XML/DDD builder), as in WMass/cmssw#42.
"""
import os

import FWCore.ParameterSet.Config as cms

OPERA_DEFAULT_VERSION = "170812"

_CPE_PRODUCERS = ("stripCPEESProducer", "StripCPEfromTrackAngleESProducer",
                  "siPixelTemplateDBObjectESProducer", "templates")

TABLES_HELP = """\
The OPERA 170812 interpolation tables are not in the CMSSW_15_0 data package.
Put merged.bin and merged.index (76 MB) where FileInPath finds them, e.g.

  D=$CMSSW_BASE/external/$SCRAM_ARCH/data/MagneticField/Interpolation/data/grid_170812_3_8t
  mkdir -p $D && cp <tables>/grid_170812_3_8t/merged.{bin,index} $D/

Sources: cms-data/MagneticField-Interpolation PR #5 (namapane), or the local
copy /work/submit/david_w/ZMass/MagneticField-Interpolation/grid_170812_3_8t.
CRAB ships that directory with config.JobType.sendExternalFolder = True
(strip the debug symbols of $CMSSW_BASE/lib first to stay under the 120 MB
sandbox limit), or point CMSSW_SEARCH_PATH to a CVMFS copy on the worker."""


def findTables(tableset):
    """Directory of `tableset` on CMSSW_SEARCH_PATH (merged file), or None."""
    rel = os.path.join("MagneticField", "Interpolation", "data", tableset)
    for base in os.environ.get("CMSSW_SEARCH_PATH", "").split(":"):
        if base and os.path.isfile(os.path.join(base, rel, "merged.bin")) \
                and os.path.isfile(os.path.join(base, rel, "merged.index")):
            return os.path.join(base, rel)
    return None


def setupOpera3DField(process, version=OPERA_DEFAULT_VERSION, routeCPEs=True, checkTables=True):
    """Add the OPERA field as process.Opera3DMagneticFieldProducer; return its label."""
    if version == "170812":
        import MagneticField.Engine.volumeBasedMagneticField_dd4hep_170812_cfi as _mf
        label = "grid_170812_3_8t"
        if checkTables and findTables(label) is None:
            raise RuntimeError("setupOpera3DField: %s tables not found on CMSSW_SEARCH_PATH.\n%s"
                               % (label, TABLES_HELP))
        # MF geometry (DD4hep) under the 'magfield' label, consumed only by
        # this producer (the default field of a Run 2 process comes from the
        # DB). DDDetectorESProducer labels its MF product with its MODULE label
        # (+ appendToDataLabel), so the source must be named 'magfield'.
        if hasattr(process, "magfield") and \
                getattr(process, "magfield").type_() != "DDDetectorESProducer":
            raise RuntimeError("setupOpera3DField: process.magfield already exists (%s)"
                               % getattr(process, "magfield").type_())
        process.magfield = _mf.DDDetectorESProducer.clone(appendToDataLabel=cms.string(""))
        process.cvhMFCompactView = _mf.DDCompactViewMFESProducer.clone()
        producer = _mf.VolumeBasedMagneticFieldESProducer.clone()
    elif version == "160812":
        from MagneticField.Engine.volumeBasedMagneticField_160812_cfi import \
            VolumeBasedMagneticFieldESProducer as _producer160812
        from MagneticField.Engine.volumeBasedMagneticField_160812_cfi import magfield as _magfield
        label = "grid_160812_3_8t"
        process.magfield = _magfield
        process.es_prefer_magfield_cvhrefit = cms.ESPrefer("XMLIdealGeometryESSource", "magfield")
        producer = _producer160812.clone()
    else:
        raise RuntimeError("setupOpera3DField: version must be '170812' or '160812'")
    producer.label = cms.untracked.string(label)
    producer.useParametrizedTrackerField = cms.bool(False)
    process.Opera3DMagneticFieldProducer = producer
    if routeCPEs:
        for cpe in _CPE_PRODUCERS:
            if hasattr(process, cpe):
                getattr(process, cpe).MagneticFieldLabel = cms.string(label)
    return label
