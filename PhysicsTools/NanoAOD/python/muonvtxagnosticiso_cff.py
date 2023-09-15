import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import ExtVar

muonvtxagniso04 = cms.EDProducer("MuonVtxAgnosticIsoProducer",
    muonInputTag = cms.InputTag("linkedObjects","muons"),
    pfCandidateInputTag = cms.InputTag("packedPFCandidates"),
    maxdr = cms.double(0.4),
    mindr = cms.double(0.005),
    maxdeltaz = cms.double(0.2),
    minptchg = cms.double(0.0),
    minptneu = cms.double(0.0),
    minptpho = cms.double(0.0),
    minptpu = cms.double(0.0),
    calculateNeutralPhoton = cms.int32(1)
)

muonvtxagniso03 = muonvtxagniso04.clone()
muonvtxagniso03.maxdr = cms.double(0.3)

muoncrosscheckiso04 = cms.EDProducer("MuonCrossCheckIsoProducer",
    muonInputTag = muonvtxagniso04.muonInputTag,
    pfCandidateInputTag = muonvtxagniso04.pfCandidateInputTag,
    pvInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
    maxdr = muonvtxagniso04.maxdr,
    mindr = muonvtxagniso04.mindr,
    maxdeltaz = muonvtxagniso04.maxdeltaz,
    minptchg = muonvtxagniso04.minptchg,
    minptneu = muonvtxagniso04.minptneu,
    minptpho = muonvtxagniso04.minptpho,
    minptpu = muonvtxagniso04.minptpu,
    calculateNeutralPhoton = muonvtxagniso04.calculateNeutralPhoton
)

muoncrosscheckiso03 = muoncrosscheckiso04.clone()
muoncrosscheckiso03.maxdr = muonvtxagniso03.maxdr


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

crossCheckIsoVariables = cms.PSet(
    crossCheckPfRelIso04_chg = ExtVar(cms.InputTag("muoncrosscheckiso04:vtxAgnosticChargedHadronIso"), float, doc="Manually computed PF relative isolation to avoid vertex selection (DR=0.4, charged component)"),
    crossCheckPfRelIso04_neu = ExtVar(cms.InputTag("muoncrosscheckiso04:vtxAgnosticNeutralHadronIso"), float, doc="Manually computed PF relative isolation to avoid vertex selection (DR=0.4, neutral component)"),
    crossCheckPfRelIso04_pho = ExtVar(cms.InputTag("muoncrosscheckiso04:vtxAgnosticPhotonIso"), float, doc="Manually computed PF relative isolation to avoid vertex selection (DR=0.4, photon component)"),
    crossCheckPfRelIso04_pu = ExtVar(cms.InputTag("muoncrosscheckiso04:vtxAgnosticPUIso"), float, doc="Manually computed PF relative isolation to avoid vertex selection (DR=0.4, PU component)"),
    crossCheckPfRelIso04_all = ExtVar(cms.InputTag("muoncrosscheckiso04:vtxAgnosticTotalIso"), float, doc="Manually computed PF relative isolation to avoid vertex selection (DR=0.4, combination with delta beta corrections)"),
    crossCheckPfRelIso03_chg = ExtVar(cms.InputTag("muoncrosscheckiso03:vtxAgnosticChargedHadronIso"), float, doc="Manually computed PF relative isolation to avoid vertex selection (DR=0.3, charged component)"),
    crossCheckPfRelIso03_neu = ExtVar(cms.InputTag("muoncrosscheckiso03:vtxAgnosticNeutralHadronIso"), float, doc="Manually computed PF relative isolation to avoid vertex selection (DR=0.3, neutral component)"),
    crossCheckPfRelIso03_pho = ExtVar(cms.InputTag("muoncrosscheckiso03:vtxAgnosticPhotonIso"), float, doc="Manually computed PF relative isolation to avoid vertex selection (DR=0.3, photon component)"),
    crossCheckPfRelIso03_pu = ExtVar(cms.InputTag("muoncrosscheckiso03:vtxAgnosticPUIso"), float, doc="Manually computed PF relative isolation to avoid vertex selection (DR=0.3, PU component)"),
    crossCheckPfRelIso03_all = ExtVar(cms.InputTag("muoncrosscheckiso03:vtxAgnosticTotalIso"), float, doc="Manually computed PF relative isolation to avoid vertex selection (DR=0.3, combination with delta beta corrections)"),
)
