import FWCore.ParameterSet.Config as cms
from CommonTools.PileupAlgos.Puppi_cff import puppiCentral, puppiForward

# WMass: PUPPI on MiniAOD with the primary vertex chosen by the leading loose
# muon (see plugins/PuppiPVRobustProducer.cc). The parameters are the ones of
# the 10_6 PRODUCTION (provenance of the NanoV9MCPostVFP_TrackFitV722_NanoProdv6
# files), not of the 10_6 Puppi_cff.py: there the MET recalibration's PUPPI v15
# tune (UpdatePuppiTuneV15: PtMaxCharged 20, EtaMinUseDeltaZ 2.4,
# PtMaxNeutralsStartSlope 20, NumOfPUVtxsForCharged 2, central etaMin -0.01)
# also reached puppiPVRobust, and the PV-robust DeepMET values of the
# production nano are only reproduced with these values.
puppiPVRobust = cms.EDProducer("PuppiPVRobustProducer",
    puppiDiagnostics = cms.bool(False),
    puppiNoLep = cms.bool(False),
    UseFromPVLooseTight = cms.bool(False),
    UseDeltaZCut = cms.bool(True),
    EtaMinUseDeltaZ = cms.double(2.4),
    DeltaZCut = cms.double(0.3),
    NumOfPUVtxsForCharged = cms.uint32(2),
    DeltaZCutForChargedFromPUVtxs = cms.double(0.2),
    PtMaxCharged = cms.double(20.),
    EtaMaxCharged = cms.double(99999.),
    PtMaxPhotons = cms.double(-1.),
    EtaMaxPhotons = cms.double(2.5),
    PtMaxNeutrals = cms.double(200.),
    PtMaxNeutralsStartSlope = cms.double(20.),
    candName = cms.InputTag('packedPFCandidates'),
    vertexName = cms.InputTag('offlineSlimmedPrimaryVertices'),
    muonName = cms.InputTag('slimmedMuons'),
    beamSpotName = cms.InputTag('offlineBeamSpot'),
    muonPtMin = cms.double(10.),
    muonVertexDzMax = cms.double(0.2),
    applyCHS = cms.bool(True),
    invertPuppi = cms.bool(False),
    useExp = cms.bool(False),
    MinPuppiWeight = cms.double(0.01),
    vtxNdofCut = cms.int32(4),
    vtxZCut = cms.double(24),
    algos = cms.VPSet(
        cms.PSet(
            etaMin = cms.vdouble(-0.01),
            etaMax = cms.vdouble(2.5),
            ptMin = cms.vdouble(0.),
            MinNeutralPt = cms.vdouble(0.2),
            MinNeutralPtSlope = cms.vdouble(0.015),
            RMSEtaSF = cms.vdouble(1.0),
            MedEtaSF = cms.vdouble(1.0),
            EtaMaxExtrap = cms.double(2.0),
            puppiAlgos = puppiCentral
        ),
        cms.PSet(
            etaMin = cms.vdouble( 2.5,  3.0),
            etaMax = cms.vdouble( 3.0, 10.0),
            ptMin = cms.vdouble( 0.0,  0.0),
            MinNeutralPt = cms.vdouble( 1.7,  2.0),
            MinNeutralPtSlope = cms.vdouble(0.08, 0.08),
            RMSEtaSF = cms.vdouble(1.20, 0.95),
            MedEtaSF = cms.vdouble(0.90, 0.75),
            EtaMaxExtrap = cms.double( 2.0),
            puppiAlgos = puppiForward
        ),
    )
)
