# The shared CVH Geant4 master as an EventSetup product (CvhMasterRecord).
#
# One instance per job, consumed via esConsumes by every CVH residual maker,
# so any number/mix of makers (single-track nominal/ideal/bs, two-track, V0,
# B->J/psi K) share the single G4 master. Replaces the old per-producer
# `CvhMaster = CvhMasterPSet.clone(...)` block.
#
# Configuration is the same CvhMasterPSet the producers used to carry. Because
# the master is now shared, `Particles` should be the UNION of every channel's
# needs -- the CvhMasterPSet default is already the full canonical set, so no
# per-job narrowing is applied here.
#
# Set the field label before use, e.g.:
#     process.cvhMasterESProducer.MagneticFieldLabel = "ScalarPot3DMf"

import FWCore.ParameterSet.Config as cms

from TrackPropagation.Geant4e.cvhMaster_cfi import CvhMasterPSet

cvhMasterESProducer = cms.ESProducer(
    "CvhMasterESProducer",
    **CvhMasterPSet.clone().parameters_()
)
