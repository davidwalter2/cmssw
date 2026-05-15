# Configuration template for the CVH G4 master (CvhMasterThread).
#
# Each consumer of the residual-maker (ResidualGlobalCorrectionMakerG4e,
# ResidualGlobalCorrectionMakerTwoTrackG4e) reads a `CvhMaster` PSet from
# its own cfi. The pset below is the canonical default; clone it into your
# residual-maker cfi as `CvhMaster = CvhMasterPSet.clone(...)`.
#
# The MagneticField subpset matches what geopro used to receive from
# SimG4Core/Application/g4SimHits_cfi; it parametrises the per-thread
# CMSFieldManager / stepper / chord-finder that CvhWorker builds on each
# TBB worker's first produce(). DD4hep / legacy DD geometry is selected by
# `g4GeometryDD4hepSource` (default False = legacy DD).

import FWCore.ParameterSet.Config as cms

from SimG4Core.Application.g4SimHits_cfi import g4SimHits as _g4SimHits

CvhMasterPSet = cms.PSet(
    g4GeometryDD4hepSource = cms.bool(False),
    UseMagneticField = cms.bool(True),
    # Label of the MagneticField ESProducer the master should consume. Empty
    # = default unlabelled MagneticField (used by tests with the stock map);
    # override to e.g. "ScalarPot3DMf" to wire in a labelled producer (the
    # CVH refit's actual production path).
    MagneticFieldLabel = cms.string(""),
    MagneticField = _g4SimHits.MagneticField.clone(),
)

# DD4hep era hookup: same toModify pattern the legacy geopro cfi used.
from Configuration.ProcessModifiers.dd4hep_cff import dd4hep
dd4hep.toModify(CvhMasterPSet, g4GeometryDD4hepSource=True)
