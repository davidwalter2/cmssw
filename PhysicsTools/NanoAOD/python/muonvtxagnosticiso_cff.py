import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import ExtVar

# WMass (from WmassNanoProd_10_6_26, PR #31): PF isolation of the muon computed
# from the packed PF candidates without any primary-vertex association -- the
# charged component takes every charged candidate within |dz(muon, cand)| <
# maxdeltaz, the PU component the rest; neutral/photon components are taken from
# the standard pfIsolationR0{3,4} unless calculateNeutralPhoton != 0.

muonvtxagniso04 = cms.EDProducer("MuonVtxAgnosticIsoProducer",
    muonInputTag = cms.InputTag("linkedObjects","muons"),
    pfCandidateInputTag = cms.InputTag("packedPFCandidates"),
    maxdr = cms.double(0.4),
    mindrchg = cms.double(1e-4),
    mindrrest = cms.double(0.01),
    maxdeltaz = cms.double(0.2),
    minptchg = cms.double(0.0),
    minptneu = cms.double(0.5),
    minptpho = cms.double(0.5),
    minptpu = cms.double(0.5),
    calculateNeutralPhoton = cms.int32(0)
)

muonvtxagniso03 = muonvtxagniso04.clone()
muonvtxagniso03.maxdr = cms.double(0.3)


vtxAgnIsoVariables = cms.PSet(
    vtxAgnPfRelIso04_chg = ExtVar(cms.InputTag("muonvtxagniso04:vtxAgnosticChargedHadronIso"), float, doc="Manually computed PF relative isolation to avoid vertex selection (DR=0.4, charged component)"),
    vtxAgnPfRelIso04_neu = ExtVar(cms.InputTag("muonvtxagniso04:vtxAgnosticNeutralHadronIso"), float, doc="Manually computed PF relative isolation to avoid vertex selection (DR=0.4, neutral component)"),
    vtxAgnPfRelIso04_pho = ExtVar(cms.InputTag("muonvtxagniso04:vtxAgnosticPhotonIso"), float, doc="Manually computed PF relative isolation to avoid vertex selection (DR=0.4, photon component)"),
    vtxAgnPfRelIso04_pu = ExtVar(cms.InputTag("muonvtxagniso04:vtxAgnosticPUIso"), float, doc="Manually computed PF relative isolation to avoid vertex selection (DR=0.4, PU component)"),
    vtxAgnPfRelIso04_all = ExtVar(cms.InputTag("muonvtxagniso04:vtxAgnosticTotalIso"), float, doc="Manually computed PF relative isolation to avoid vertex selection (DR=0.4, combination with delta beta corrections)"),
    vtxAgnPfRelIso03_chg = ExtVar(cms.InputTag("muonvtxagniso03:vtxAgnosticChargedHadronIso"), float, doc="Manually computed PF relative isolation to avoid vertex selection (DR=0.3, charged component)"),
    vtxAgnPfRelIso03_neu = ExtVar(cms.InputTag("muonvtxagniso03:vtxAgnosticNeutralHadronIso"), float, doc="Manually computed PF relative isolation to avoid vertex selection (DR=0.3, neutral component)"),
    vtxAgnPfRelIso03_pho = ExtVar(cms.InputTag("muonvtxagniso03:vtxAgnosticPhotonIso"), float, doc="Manually computed PF relative isolation to avoid vertex selection (DR=0.3, photon component)"),
    vtxAgnPfRelIso03_pu = ExtVar(cms.InputTag("muonvtxagniso03:vtxAgnosticPUIso"), float, doc="Manually computed PF relative isolation to avoid vertex selection (DR=0.3, PU component)"),
    vtxAgnPfRelIso03_all = ExtVar(cms.InputTag("muonvtxagniso03:vtxAgnosticTotalIso"), float, doc="Manually computed PF relative isolation to avoid vertex selection (DR=0.3, combination with delta beta corrections)"),
)

muonVtxAgnosticIsoTask = cms.Task(muonvtxagniso04, muonvtxagniso03)


def nanoAOD_addVtxAgnosticIso(process):
    """Run the two isolation producers on the muons of the Muon table and add
    the ten vtxAgnPfRelIso0{3,4}_* columns (data and MC)."""
    process.muonvtxagniso04 = muonvtxagniso04
    process.muonvtxagniso03 = muonvtxagniso03
    process.muonvtxagniso04.muonInputTag = process.muonTable.src
    process.muonvtxagniso03.muonInputTag = process.muonTable.src
    process.muonTablesTask.add(process.muonvtxagniso04, process.muonvtxagniso03)
    for name in vtxAgnIsoVariables.parameterNames_():
        setattr(process.muonTable.externalVariables, name, getattr(vtxAgnIsoVariables, name))
    return process
