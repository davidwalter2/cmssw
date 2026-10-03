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
    # G4 particle names to register in G4ErrorPhysicsListForCVH. Default is
    # the full canonical CVH set; narrow it per job to skip per-thread
    # process registration for particles you won't propagate.
    # For J/psi -> mu mu set Particles = cms.vstring("mu+","mu-").
    # For K_s -> pi pi: ("pi+","pi-"). For B -> J/psi K: ("mu+","mu-","kaon+","kaon-").
    # For Lambda -> p pi: ("anti_proton","pi+","pi-").
    # An unknown name aborts master-thread initG4 with a clear G4Exception.
    #
    # Always-on mandatory particles (registered regardless of this list):
    #   gamma, e+, e-  -- G4PhysicsListHelper::CheckParticleList aborts
    #                     otherwise (Run0101 "Missing EM basic particle").
    #   mu+, mu-, proton -- referenced by G4TablesForExtrapolatorForCVH's
    #                       dE/dx + range table generation. Race condition
    #                       at high thread count (Run0271 / PART10116) if
    #                       they aren't pre-registered on master.
    # Cost of the mandatory set is small (just particle definitions); EM
    # physics processes are only attached for entries explicitly listed
    # here (the gamma Compton/conv/photoelectric block needs "gamma" in
    # the list to attach).
    Particles = cms.vstring(
        "gamma",
        "e+", "e-",
        "mu+", "mu-",
        "pi+", "pi-",
        "kaon+", "kaon-",
        "proton", "anti_proton",
    ),
)

# DD4hep era hookup: same toModify pattern the legacy geopro cfi used.
from Configuration.ProcessModifiers.dd4hep_cff import dd4hep
dd4hep.toModify(CvhMasterPSet, g4GeometryDD4hepSource=True)
