import FWCore.ParameterSet.Config as cms
from CommonTools.PileupAlgos.puppiPVRobust_cff import puppiPVRobust

# WMass: DeepMET on the PV-robust inputs (plugins/DeepMETPVRobustProducer.cc),
# the two models of the 10_6 production, leptons removed from the inputs and
# added back (ignore_leptons) as in the 10_6 nano_cff
deepMETsPVRobust = cms.EDProducer("DeepMETPVRobustProducer",
    pf_src = cms.InputTag("packedPFCandidates"),
    PFPVRobustDxy = cms.InputTag("puppiPVRobust", "PFPVRobustDxy"),
    PFPVRobustDz = cms.InputTag("puppiPVRobust", "PFPVRobustDz"),
    PFPVRobustPuppiWeight = cms.InputTag("puppiPVRobust", "PFPVRobustPuppiWeight"),
    usePUPPI = cms.bool(True),
    ignore_leptons = cms.bool(True),
    norm_factor = cms.double(50.),
    max_n_pf = cms.uint32(4500),
    graph_path = cms.string("RecoMET/METPUSubtraction/data/deepmet_pvrobust/deepmet_pvrobust.pb"),
)
deepMETsPVRobustNoPUPPI = deepMETsPVRobust.clone(
    usePUPPI = False,
    graph_path = "RecoMET/METPUSubtraction/data/deepmet_pvrobust/deepmet_pvrobust_nopuppi.pb",
)

deepMETPVRobustTask = cms.Task(puppiPVRobust, deepMETsPVRobust, deepMETsPVRobustNoPUPPI)
