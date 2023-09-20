import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import ExtVar

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
