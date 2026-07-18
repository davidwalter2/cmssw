import FWCore.ParameterSet.Config as cms

from PhysicsTools.NanoAOD.nano_eras_cff import *
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.simplePATMuonFlatTableProducer_cfi import simplePATMuonFlatTableProducer

import PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi

# this below is used only in some eras
slimmedMuonsUpdated = cms.EDProducer("PATMuonUpdater",
    src = cms.InputTag("slimmedMuons"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    computeMiniIso = cms.bool(False),
    fixDxySign = cms.bool(True),
    pfCandsForMiniIso = cms.InputTag("packedPFCandidates"),
    miniIsoParams = PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi.patMuons.miniIsoParams, # so they're in sync
    recomputeMuonBasicSelectors = cms.bool(False),
    recomputeSoftMuonMvaRun3 = cms.bool(False),
    softMvaRun3Model = PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi.patMuons.softMvaRun3Model,
)

(run2_nanoAOD_106Xv2 | run3_nanoAOD_pre142X).toModify(
    slimmedMuonsUpdated, recomputeMuonBasicSelectors=True, recomputeSoftMuonMvaRun3=True,
)

isoForMu = cms.EDProducer("MuonIsoValueMapProducer",
    src = cms.InputTag("slimmedMuonsUpdated"),
    relative = cms.bool(False),
    rho_MiniIso = cms.InputTag("fixedGridRhoFastjetAll"),
    EAFile_MiniIso = cms.FileInPath("PhysicsTools/NanoAOD/data/effAreaMuons_cone03_pfNeuHadronsAndPhotons_94X.txt"),
)

ptRatioRelForMu = cms.EDProducer("MuonJetVarProducer",
    srcJet = cms.InputTag("updatedJetsPuppi"),
    srcLep = cms.InputTag("slimmedMuonsUpdated"),
    srcVtx = cms.InputTag("offlineSlimmedPrimaryVertices"),
)

muonMVAID = cms.EDProducer("EvaluateMuonMVAID",
    src = cms.InputTag("slimmedMuonsUpdated"),
    weightFile =  cms.FileInPath("RecoMuon/MuonIdentification/data/mvaID.onnx"),
    backend = cms.string('ONNX'),
    name = cms.string("muonMVAID"),
    outputTensorName= cms.string("probabilities"),
    inputTensorName= cms.string("float_input"),
    outputNames = cms.vstring(["probGOOD", "wpMedium", "wpTight"]),
    batch_eval =cms.bool(True),
    outputFormulas = cms.vstring(["at(1)", "? at(1) > 0.08 ? 1 : 0", "? at(1) > 0.20 ? 1 : 0"]),
    variables = cms.VPSet(
        cms.PSet( name = cms.string("LepGood_global_muon"), expr = cms.string("isGlobalMuon")),
        cms.PSet( name = cms.string("LepGood_validFraction"), expr = cms.string("?innerTrack.isNonnull?innerTrack().validFraction:-99")),
        cms.PSet( name = cms.string("Muon_norm_chi2_extended")),
        cms.PSet( name = cms.string("LepGood_local_chi2"), expr = cms.string("combinedQuality().chi2LocalPosition")),
        cms.PSet( name = cms.string("LepGood_kink"), expr = cms.string("combinedQuality().trkKink")),
        cms.PSet( name = cms.string("LepGood_segmentComp"), expr = cms.string("segmentCompatibility")),
        cms.PSet( name = cms.string("Muon_n_Valid_hits_extended")),
        cms.PSet( name = cms.string("LepGood_n_MatchedStations"), expr = cms.string("numberOfMatchedStations()")),
        cms.PSet( name = cms.string("LepGood_Valid_pixel"), expr = cms.string("?innerTrack.isNonnull()?innerTrack().hitPattern().numberOfValidPixelHits():-99")),
        cms.PSet( name = cms.string("LepGood_tracker_layers"), expr = cms.string("?innerTrack.isNonnull()?innerTrack().hitPattern().trackerLayersWithMeasurement():-99")),
        cms.PSet( name = cms.string("LepGood_pt"), expr = cms.string("pt")),
        cms.PSet( name = cms.string("LepGood_eta"), expr = cms.string("eta")),
    )
)






slimmedMuonsWithUserData = cms.EDProducer("PATMuonUserDataEmbedder",
     src = cms.InputTag("slimmedMuonsUpdated"),
     userFloats = cms.PSet(
        miniIsoChg = cms.InputTag("isoForMu:miniIsoChg"),
        miniIsoAll = cms.InputTag("isoForMu:miniIsoAll"),
        ptRatio = cms.InputTag("ptRatioRelForMu:ptRatio"),
        ptRel = cms.InputTag("ptRatioRelForMu:ptRel"),
        jetNDauChargedMVASel = cms.InputTag("ptRatioRelForMu:jetNDauChargedMVASel"),
        mvaIDMuon_wpMedium = cms.InputTag("muonMVAID:wpMedium"),
        mvaIDMuon_wpTight = cms.InputTag("muonMVAID:wpTight"),
        mvaIDMuon = cms.InputTag("muonMVAID:probGOOD")
     ),
     userCands = cms.PSet(
        jetForLepJetVar = cms.InputTag("ptRatioRelForMu:jetForLepJetVar") # warning: Ptr is null if no match is found
     ),
)


finalMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("slimmedMuonsWithUserData"),
    cut = cms.string("pt > 15 || (pt > 3 && (passed('CutBasedIdLoose') || passed('SoftCutBasedId') || passed('SoftMvaId') || passed('CutBasedIdGlobalHighPt') || passed('CutBasedIdTrkHighPt')))")
)

# lower the muon pt threshold to 2 GeV
(run3_nanoAOD_2025 | run3_nanoAOD_devel).toModify(
    finalMuons,
    cut = cms.string("pt > 15 || (pt > 2 && (passed('CutBasedIdLoose') || passed('SoftCutBasedId') || passed('SoftMvaId') || passed('CutBasedIdGlobalHighPt') || passed('CutBasedIdTrkHighPt')))")
)


finalLooseMuons = cms.EDFilter("PATMuonRefSelector", # for isotrack cleaning
    src = cms.InputTag("slimmedMuonsWithUserData"),
    cut = cms.string("pt > 3 && track.isNonnull && isLooseMuon")
)

muonPROMPTMVA= cms.EDProducer("MuonBaseMVAValueMapProducer",
    src = cms.InputTag("linkedObjects","muons"),
    weightFile =  cms.FileInPath("PhysicsTools/NanoAOD/data/mu_BDTG_2022.weights.xml"),
    backend = cms.string("TMVA"),
    name = cms.string("muonPROMPTMVA"),
    isClassifier = cms.bool(True),
    variables = cms.VPSet(
        cms.PSet( name = cms.string("LepGood_pt"), expr = cms.string("pt")),
        cms.PSet( name = cms.string("LepGood_eta"), expr = cms.string("eta")),
        cms.PSet( name = cms.string("LepGood_pfRelIso03_all"), expr = cms.string("(pfIsolationR03().sumChargedHadronPt + max(pfIsolationR03().sumNeutralHadronEt + pfIsolationR03().sumPhotonEt - pfIsolationR03().sumPUPt/2,0.0))/pt")),
        cms.PSet( name = cms.string("LepGood_miniRelIsoCharged"), expr = cms.string("userFloat('miniIsoChg')/pt")),
        cms.PSet( name = cms.string("LepGood_miniRelIsoNeutral"), expr = cms.string("(userFloat('miniIsoAll')-userFloat('miniIsoChg'))/pt")),
        cms.PSet( name = cms.string("LepGood_jetNDauChargedMVASel"), expr = cms.string("?userCand('jetForLepJetVar').isNonnull()?userFloat('jetNDauChargedMVASel'):0")),
        cms.PSet( name = cms.string("LepGood_jetPtRelv2"), expr = cms.string("?userCand('jetForLepJetVar').isNonnull()?userFloat('ptRel'):0")),
        cms.PSet( name = cms.string("LepGood_jetDF"), expr = cms.string("?userCand('jetForLepJetVar').isNonnull()?max(userCand('jetForLepJetVar').bDiscriminator('pfDeepFlavourJetTags:probbb')+userCand('jetForLepJetVar').bDiscriminator('pfDeepFlavourJetTags:probb')+userCand('jetForLepJetVar').bDiscriminator('pfDeepFlavourJetTags:problepb'),0.0):0.0")),
        cms.PSet( name = cms.string("LepGood_jetPtRatio"), expr = cms.string("?userCand('jetForLepJetVar').isNonnull()?min(userFloat('ptRatio'),1.5):1.0/(1.0+(pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - pfIsolationR04().sumPUPt/2,0.0))/pt)")),
        cms.PSet( name = cms.string("LepGood_sip3d"), expr = cms.string("abs(dB('PV3D')/edB('PV3D'))")),
        cms.PSet( name = cms.string("LepGood_dxy"), expr = cms.string("log(abs(dB('PV2D')))")),
        cms.PSet( name = cms.string("LepGood_dz"), expr = cms.string("log(abs(dB('PVDZ')))")),
        cms.PSet( name = cms.string("LepGood_segmentComp"), expr = cms.string("segmentCompatibility")),
        )
)

_legacy_muon_BDT_variable = cms.VPSet(
    cms.PSet( name = cms.string("LepGood_pt"), expr = cms.string("pt")),
    cms.PSet( name = cms.string("LepGood_eta"), expr = cms.string("eta")),
    cms.PSet( name = cms.string("LepGood_jetNDauChargedMVASel"), expr = cms.string("?userCand('jetForLepJetVar').isNonnull()?userFloat('jetNDauChargedMVASel'):0")),
    cms.PSet( name = cms.string("LepGood_miniRelIsoCharged"), expr = cms.string("userFloat('miniIsoChg')/pt")),
    cms.PSet( name = cms.string("LepGood_miniRelIsoNeutral"), expr = cms.string("(userFloat('miniIsoAll')-userFloat('miniIsoChg'))/pt")),
    cms.PSet( name = cms.string("LepGood_jetPtRelv2"), expr = cms.string("?userCand('jetForLepJetVar').isNonnull()?userFloat('ptRel'):0")),
    cms.PSet( name = cms.string("LepGood_jetDF"), expr = cms.string("?userCand('jetForLepJetVar').isNonnull()?max(userCand('jetForLepJetVar').bDiscriminator('pfDeepFlavourJetTags:probbb')+userCand('jetForLepJetVar').bDiscriminator('pfDeepFlavourJetTags:probb')+userCand('jetForLepJetVar').bDiscriminator('pfDeepFlavourJetTags:problepb'),0.0):0.0")),
    cms.PSet( name = cms.string("LepGood_jetPtRatio"), expr = cms.string("?userCand('jetForLepJetVar').isNonnull()?min(userFloat('ptRatio'),1.5):1.0/(1.0+(pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - pfIsolationR04().sumPUPt/2,0.0))/pt)")),
    cms.PSet( name = cms.string("LepGood_dxy"), expr = cms.string("log(abs(dB('PV2D')))")),
    cms.PSet( name = cms.string("LepGood_sip3d"), expr = cms.string("abs(dB('PV3D')/edB('PV3D'))")),
    cms.PSet( name = cms.string("LepGood_dz"), expr = cms.string("log(abs(dB('PVDZ')))")),
    cms.PSet( name = cms.string("LepGood_segmentComp"), expr = cms.string("segmentCompatibility")),
)

muonMVALowPt = muonPROMPTMVA.clone(
    weightFile =  cms.FileInPath("PhysicsTools/NanoAOD/data/mu_BDTG_lowpt.weights.xml"),
    name = cms.string("muonMVALowPt"),
    variables = _legacy_muon_BDT_variable
)

run2_muon_2016.toModify(
    muonPROMPTMVA,
    weightFile = "PhysicsTools/NanoAOD/data/mu_BDTG_2016.weights.xml",
    variables =	_legacy_muon_BDT_variable
)

(run2_muon_2017 | run2_muon_2018).toModify(
    muonPROMPTMVA,
    weightFile =  cms.FileInPath("PhysicsTools/NanoAOD/data/mu_BDTG_2017.weights.xml"),
    variables = _legacy_muon_BDT_variable
)

from PhysicsTools.PatAlgos.muonTagInfos_cfi import muonTagInfos as _muonTagInfos
muonPNetVariables = _muonTagInfos.clone(
    src = cms.InputTag("linkedObjects","muons"),
    leptonVars = cms.PSet(
        MuonSelected_LepGood_pt = cms.string("pt"),
        MuonSelected_LepGood_eta = cms.string("eta"),
        MuonSelected_LepGood_jetNDauChargedMVASel = cms.string("?userCand('jetForLepJetVar').isNonnull()?userFloat('jetNDauChargedMVASel'):0"),
        MuonSelected_LepGood_miniRelIsoCharged = cms.string("userFloat('miniIsoChg')/pt"),
        MuonSelected_LepGood_miniRelIsoNeutral = cms.string("(userFloat('miniIsoAll')-userFloat('miniIsoChg'))/pt"),
        MuonSelected_LepGood_jetPtRelv2 = cms.string("?userCand('jetForLepJetVar').isNonnull()?userFloat('ptRel'):0"),
        MuonSelected_LepGood_jetDF = cms.string("?userCand('jetForLepJetVar').isNonnull()?max(userCand('jetForLepJetVar').bDiscriminator('pfDeepFlavourJetTags:probbb')+userCand('jetForLepJetVar').bDiscriminator('pfDeepFlavourJetTags:probb')+userCand('jetForLepJetVar').bDiscriminator('pfDeepFlavourJetTags:problepb'),0.0):0.0"),
        MuonSelected_LepGood_jetPtRatio = cms.string("?userCand('jetForLepJetVar').isNonnull()?min(userFloat('ptRatio'),1.5):1.0/(1.0+(pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - pfIsolationR04().sumPUPt/2,0.0))/pt)"),
        MuonSelected_dxy = cms.string("log(abs(dB('PV2D')))"),
        MuonSelected_sip3d = cms.string("abs(dB('PV3D')/edB('PV3D'))"),
        MuonSelected_dz = cms.string("log(abs(dB('PVDZ')))"),
        MuonSelected_LepGood_dz = cms.string("log(abs(dB('PVDZ')))"),
        MuonSelected_segmentComp = cms.string("segmentCompatibility"),
        MuonSelected_global_muon = cms.string("isGlobalMuon"),
        MuonSelected_validFraction = cms.string("?innerTrack.isNonnull?innerTrack().validFraction:-99"),
        MuonSelected_local_chi2 = cms.string("combinedQuality().chi2LocalPosition"),
        MuonSelected_kink = cms.string("combinedQuality().trkKink"),
        MuonSelected_n_MatchedStations = cms.string("numberOfMatchedStations()"),
        MuonSelected_Valid_pixel = cms.string("?innerTrack.isNonnull()?innerTrack().hitPattern().numberOfValidPixelHits():-99"),
        MuonSelected_tracker_layers = cms.string("?innerTrack.isNonnull()?innerTrack().hitPattern().trackerLayersWithMeasurement():-99"),
        MuonSelected_mvaId=cms.string("userFloat('mvaIDMuon')"),
    ),
    leptonVarsExt = cms.PSet(
        MuonSelected_mvaTTH=cms.InputTag("muonPROMPTMVA"),
    ),
    pfVars = cms.PSet(
        PF_pt=cms.string("pt"),
        PF_charge=cms.string("charge"),
        PF_isElectron=cms.string("?abs(pdgId)==11?1:0"),
        PF_isMuon=cms.string("?abs(pdgId)==13?1:0"),
        PF_isNeutralHadron=cms.string("?abs(pdgId)==130?1:0"),
        PF_isPhoton=cms.string("?abs(pdgId)==22?1:0"),
        PF_isChargedHadron=cms.string("?abs(pdgId)==211?1:0"),
        PF_puppiWeightNoLep=cms.string("puppiWeightNoLep"),
        PF_fromPV=cms.string("fromPV"),
        PF_numberOfPixelHits=cms.string("numberOfPixelHits"),
        PF_dzSig=cms.string("?hasTrackDetails?dz/max(dzError,1.e-6):0"),
        PF_dxySig=cms.string("?hasTrackDetails?dxy/max(dxyError,1.e-6):0"),
        PF_hcalFraction=cms.string("hcalFraction"),
        PF_trackerLayersWithMeasurement=cms.string("?hasTrackDetails?bestTrack().hitPattern().trackerLayersWithMeasurement:0"),
        PF_mask=cms.string("1"),
    ),
    svVars = cms.PSet(
        SV_eta=cms.string("eta"),
        SV_phi=cms.string("phi"),
        SV_pt=cms.string("pt"),
        SV_ndof=cms.string("vertexNdof"),
        SV_chi2=cms.string("vertexChi2"),
        SV_nTracks=cms.string("numberOfDaughters"),
        SV_mass=cms.string("mass"),
        SV_mask=cms.string("1"),
    ),
)

from PhysicsTools.PatAlgos.muonPNetTags_cfi import muonPNetTags as _muonPNetTags
muonPNetScores = _muonPNetTags.clone(
    src = cms.InputTag("muonPNetVariables"),
    srcLeps = cms.InputTag("linkedObjects", "muons"),
    model_path = 'PhysicsTools/NanoAOD/data/PNetMuonId/model.onnx',
    preprocess_json = 'PhysicsTools/NanoAOD/data/PNetMuonId/preprocess.json',
    flav_names = cms.vstring(["light", "prompt", "tau", "heavy"]),
)

from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *
muonBSConstrain = cms.EDProducer("MuonBeamspotConstraintValueMapProducer",
    src = cms.InputTag("linkedObjects","muons"),
)

muonTable = simplePATMuonFlatTableProducer.clone(
    src = cms.InputTag("linkedObjects","muons"),
    name = cms.string("Muon"),
    doc  = cms.string("slimmedMuons after basic selection (" + finalMuons.cut.value()+")"),
    variables = cms.PSet(CandVars,
        ptErr   = Var("bestTrack().ptError()", float, doc = "ptError of the muon track", precision=6),
        tunepRelPt = Var("tunePMuonBestTrack().pt/pt",float,doc="TuneP relative pt, tunePpt/pt",precision=6),
        tuneP_pterr = Var("tunePMuonBestTrack().ptError()", float, doc = "pTerr from tunePMuonBestTrack", precision=6),
        tuneP_charge = Var("? tunePMuonBestTrack().isNonnull() && tunePMuonBestTrack().isAvailable() ? tunePMuonBestTrack().charge(): -99", float, doc="tunePMuonBestTrack() charge",precision=6),
        dz = Var("dB('PVDZ')",float,doc="dz (with sign) wrt first PV, in cm",precision=10),
        dzErr = Var("abs(edB('PVDZ'))",float,doc="dz uncertainty, in cm",precision=6),
        dxybs = Var("dB('BS2D')",float,doc="dxy (with sign) wrt the beam spot, in cm",precision=10),
        dxybsErr = Var("edB('BS2D')",float,doc="dxy uncertainty wrt the beam spot, in cm", precision=6),
        dxy = Var("dB('PV2D')",float,doc="dxy (with sign) wrt first PV, in cm",precision=10),
        dxyErr = Var("edB('PV2D')",float,doc="dxy uncertainty, in cm",precision=6),
        ip3d = Var("abs(dB('PV3D'))",float,doc="3D impact parameter wrt first PV, in cm",precision=10),
        sip3d = Var("abs(dB('PV3D')/edB('PV3D'))",float,doc="3D impact parameter significance wrt first PV",precision=10),
        segmentComp   = Var("segmentCompatibility()", float, doc = "muon segment compatibility", precision=14), # keep higher precision since people have cuts with 3 digits on this
        nStations = Var("numberOfMatchedStations", "uint8", doc = "number of matched stations with default arbitration (segment & track)"),
        nTrackerLayers = Var("?track.isNonnull?innerTrack().hitPattern().trackerLayersWithMeasurement():0", "uint8", doc = "number of layers in the tracker"),
        bestTrackType = Var("muonBestTrackType()", "uint8", doc = "Type of track used (1=inner, 2=STA, 3=global, 4=TPFMS, 5=Picky, 6=DYT)"),
        highPurity = Var("?track.isNonnull?innerTrack().quality('highPurity'):0", bool, doc = "inner track is high purity"),
        jetIdx = Var("?hasUserCand('jet')?userCand('jet').key():-1", "int16", doc="index of the associated jet (-1 if none)"),
        svIdx = Var("?hasUserCand('vertex')?userCand('vertex').key():-1", "int16", doc="index of matching secondary vertex"),
        tkRelIso = Var("isolationR03().sumPt/pt",float,doc="Tracker-based relative isolation dR=0.3 for highPt, trkIso/pt",precision=6),
        miniPFRelIso_chg = Var("userFloat('miniIsoChg')/pt",float,doc="mini PF relative isolation, charged component"),
        miniPFRelIso_all = Var("userFloat('miniIsoAll')/pt",float,doc="mini PF relative isolation, total (with scaled rho*EA PU corrections)"),
        pfRelIso03_chg = Var("pfIsolationR03().sumChargedHadronPt/pt",float,doc="PF relative isolation dR=0.3, charged component"),
        pfRelIso03_all = Var("(pfIsolationR03().sumChargedHadronPt + max(pfIsolationR03().sumNeutralHadronEt + pfIsolationR03().sumPhotonEt - pfIsolationR03().sumPUPt/2,0.0))/pt",float,doc="PF relative isolation dR=0.3, total (deltaBeta corrections)"),
        pfRelIso04_all = Var("(pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - pfIsolationR04().sumPUPt/2,0.0))/pt",float,doc="PF relative isolation dR=0.4, total (deltaBeta corrections)"),
        jetRelIso = Var("?userCand('jetForLepJetVar').isNonnull()?(1./userFloat('ptRatio'))-1.:-1.",float,doc="Relative isolation in matched jet (1/ptRatio-1), -1 if none",precision=8),
        jetPtRelv2 = Var("?userCand('jetForLepJetVar').isNonnull()?userFloat('ptRel'):0",float,doc="Relative momentum of the lepton with respect to the closest jet after subtracting the lepton",precision=8),
        jetDF = Var("?userCand('jetForLepJetVar').isNonnull()?max(userCand('jetForLepJetVar').bDiscriminator('pfDeepFlavourJetTags:probbb')+userCand('jetForLepJetVar').bDiscriminator('pfDeepFlavourJetTags:probb')+userCand('jetForLepJetVar').bDiscriminator('pfDeepFlavourJetTags:problepb'),0.0):0.0",float,doc="value of the DEEPJET b tagging algorithm discriminator of the associated jet (0 if none)",precision=8,lazyEval=True),
        tightCharge = Var("?(muonBestTrack().ptError()/muonBestTrack().pt() < 0.2)?2:0", "uint8", doc="Tight charge criterion using pterr/pt of muonBestTrack (0:fail, 2:pass)"),
        looseId  = Var("passed('CutBasedIdLoose')",bool, doc="muon is loose muon"),
        isPFcand = Var("isPFMuon",bool,doc="muon is PF candidate"),
        isGlobal = Var("isGlobalMuon",bool,doc="muon is global muon"),
        isTracker = Var("isTrackerMuon",bool,doc="muon is tracker muon"),
        isStandalone = Var("isStandAloneMuon",bool,doc="muon is a standalone muon"),
        mediumId = Var("passed('CutBasedIdMedium')",bool,doc="cut-based ID, medium WP"),
        mediumPromptId = Var("passed('CutBasedIdMediumPrompt')",bool,doc="cut-based ID, medium prompt WP"),
        tightId = Var("passed('CutBasedIdTight')",bool,doc="cut-based ID, tight WP"),
        softId = Var("passed('SoftCutBasedId')",bool,doc="soft cut-based ID"),
        softMvaId = Var("passed('SoftMvaId')",bool,doc="soft MVA ID"),
        softMva = Var("softMvaValue()",float,doc="soft MVA ID score",precision=6),
        softMvaRun3 = Var("softMvaRun3Value()",float,doc="soft MVA Run3 ID score",precision=6),
        highPtId = Var("?passed('CutBasedIdGlobalHighPt')?2:passed('CutBasedIdTrkHighPt')","uint8",doc="high-pT cut-based ID (1 = tracker high pT, 2 = global high pT, which includes tracker high pT)"),
        pfIsoId = Var("passed('PFIsoVeryLoose')+passed('PFIsoLoose')+passed('PFIsoMedium')+passed('PFIsoTight')+passed('PFIsoVeryTight')+passed('PFIsoVeryVeryTight')","uint8",doc="PFIso ID from miniAOD selector (1=PFIsoVeryLoose, 2=PFIsoLoose, 3=PFIsoMedium, 4=PFIsoTight, 5=PFIsoVeryTight, 6=PFIsoVeryVeryTight)"),
        tkIsoId = Var("?passed('TkIsoTight')?2:passed('TkIsoLoose')","uint8",doc="TkIso ID (1=TkIsoLoose, 2=TkIsoTight)"),
        miniIsoId = Var("passed('MiniIsoLoose')+passed('MiniIsoMedium')+passed('MiniIsoTight')+passed('MiniIsoVeryTight')","uint8",doc="MiniIso ID from miniAOD selector (1=MiniIsoLoose, 2=MiniIsoMedium, 3=MiniIsoTight, 4=MiniIsoVeryTight)"),
        mvaMuID = Var("userFloat('mvaIDMuon')", float, doc="MVA-based ID score",precision=6),
        mvaMuID_WP = Var("userFloat('mvaIDMuon_wpMedium') + userFloat('mvaIDMuon_wpTight')","uint8",doc="MVA-based ID selector WPs (1=MVAIDwpMedium,2=MVAIDwpTight)"),
        multiIsoId = Var("?passed('MultiIsoMedium')?2:passed('MultiIsoLoose')","uint8",doc="MultiIsoId from miniAOD selector (1=MultiIsoLoose, 2=MultiIsoMedium)"),
        puppiIsoId = Var("passed('PuppiIsoLoose')+passed('PuppiIsoMedium')+passed('PuppiIsoTight')", "uint8", doc="PuppiIsoId from miniAOD selector (1=Loose, 2=Medium, 3=Tight)"),
        triggerIdLoose = Var("passed('TriggerIdLoose')",bool,doc="TriggerIdLoose ID"),
        inTimeMuon = Var("passed('InTimeMuon')",bool,doc="inTimeMuon ID"),
        jetNDauCharged = Var("?userCand('jetForLepJetVar').isNonnull()?userFloat('jetNDauChargedMVASel'):0", "uint8", doc="number of charged daughters of the closest jet"),
        VXBS_Cov00 = Var("? tunePMuonBestTrack().isNonnull() && tunePMuonBestTrack().isAvailable() ? tunePMuonBestTrack().covariance(0,0) : -999",float,doc="0, 0 element of the VXBS Covariance matrix", precision=16),
        VXBS_Cov03 = Var("? tunePMuonBestTrack().isNonnull() && tunePMuonBestTrack().isAvailable() ? tunePMuonBestTrack().covariance(0,3) : -999",float,doc="0, 3 element of the VXBS Covariance matrix", precision=16),
        VXBS_Cov33 = Var("? tunePMuonBestTrack().isNonnull() && tunePMuonBestTrack().isAvailable() ? tunePMuonBestTrack().covariance(3,3) : -999",float,doc="3, 3 element of the VXBS Covariance matrix", precision=16),
        ),
    externalVariables = cms.PSet(
        promptMVA = ExtVar(cms.InputTag("muonPROMPTMVA"),float, doc="Prompt MVA lepton ID score. Corresponds to the previous mvaTTH",precision=14),
        mvaLowPt = ExtVar(cms.InputTag("muonMVALowPt"),float, doc="Low pt muon ID score",precision=14),
        pnScore_prompt = ExtVar(cms.InputTag("muonPNetScores:prompt"),float, doc="PNet muon ID score for lepton from W/Z/H bosons", precision=14),
        pnScore_heavy = ExtVar(cms.InputTag("muonPNetScores:heavy"),float, doc="PNet muon ID score for lepton from B or D hadrons", precision=14),
        pnScore_light = ExtVar(cms.InputTag("muonPNetScores:light"),float, doc="PNet muon ID score for lepton from hadrons w/o b or c quarks OR w/o generator matching", precision=14),
        pnScore_tau = ExtVar(cms.InputTag("muonPNetScores:tau"),float, doc="PNet muon ID score for decay of tau to light leptons (mu)", precision=14),
        fsrPhotonIdx = ExtVar(cms.InputTag("leptonFSRphotons:muFsrIndex"), "int16", doc="Index of the lowest-dR/ET2 among associated FSR photons"),
        bsConstrainedPt = ExtVar(cms.InputTag("muonBSConstrain:muonBSConstrainedPt"),float, doc="pT with beamspot constraint",precision=-1),
        bsConstrainedPtErr = ExtVar(cms.InputTag("muonBSConstrain:muonBSConstrainedPtErr"),float, doc="pT error with beamspot constraint ",precision=6),
        bsConstrainedChi2 = ExtVar(cms.InputTag("muonBSConstrain:muonBSConstrainedChi2"),float, doc="chi2 of beamspot constraint",precision=6),
    ),
)

# Increase precision of eta and phi
muonTable.variables.eta.precision = 16
muonTable.variables.phi.precision = 16

# --- CVH single-muon-track refit (WMass custom NanoAOD) --------------------
# Turns muon inner tracks into reco::Track + pat::Muon association, then runs
# the CVH refit. `trackrefit` (nominal, real geometry) runs for data and MC;
# the MC-only `trackrefitideal` (ideal geometry) and `trackrefitbs` (beamspot
# constraint) variants run too when isMC=True. All share the one EventSetup G4
# master (cvhMasterESProducer / CvhMasterRecord); coexistence is what Stage A2
# unblocked. See nanoAOD_addCvhMuonBranches.
from Analysis.HitAnalyzer.ResidualGlobalCorrectionMakerMuonG4e_cfi import ResidualGlobalCorrectionMakerMuonG4e  # noqa: E402

tracksfrommuons = cms.EDProducer("TrackProducerFromPatMuons",
    src = cms.InputTag("linkedObjects", "muons"),
    innerTrackOnly = cms.bool(False),
    ptMin = cms.double(-1.),
)

trackrefit = ResidualGlobalCorrectionMakerMuonG4e.clone()
# MC-only variants (added to the process only when isMC=True).
trackrefitideal = ResidualGlobalCorrectionMakerMuonG4e.clone(useIdealGeometry = cms.bool(True))
trackrefitbs = ResidualGlobalCorrectionMakerMuonG4e.clone(bsConstraint = cms.bool(True))
# Union of the correction-parameter indices touched by the nominal + ideal
# refits, so the downstream fit has one consistent index vector (matches the
# 10_6-tip scheme). GlobalIdxProducer merges exactly two inputs.
mergedGlobalIdxs = cms.EDProducer("GlobalIdxProducer",
    src0 = cms.InputTag("trackrefit", "globalIdxs"),
    src1 = cms.InputTag("trackrefitideal", "globalIdxs"),
)

# --- Dimuon (two-track) CVH refit ------------------------------------------
# Opposite-sign muon-track pairs in a mass window -> generic two-track CVH refit
# (no mass constraint, resonance-agnostic) -> per-pair "Dimuon" NanoAOD table.
# Shares the one G4 master with the single-track refit (Stage A2). Data + MC.
from Analysis.HitAnalyzer.diMuonTrackVertexCandidates_cfi import diMuonTrackVertexCandidates  # noqa: E402
from Analysis.HitAnalyzer.ResidualGlobalCorrectionMakerDiMuonG4e_cfi import ResidualGlobalCorrectionMakerDiMuonG4e  # noqa: E402

trackrefitdimuon = ResidualGlobalCorrectionMakerDiMuonG4e.clone()

dimuonTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("diMuonTrackVertexCandidates"),
    cut = cms.string(""),
    name = cms.string("Dimuon"),
    doc = cms.string("Opposite-sign dimuon pairs with the CVH two-track refit"),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(CandVars),
    externalVariables = cms.PSet(
        cvhMass = ExtVar(cms.InputTag("trackrefitdimuon:corMass"), float, doc="CVH two-track refit dimuon mass", precision=-1),
        cvhMassErr = ExtVar(cms.InputTag("trackrefitdimuon:corMassErr"), float, doc="CVH refit dimuon mass uncertainty", precision=10),
        cvhPt = ExtVar(cms.InputTag("trackrefitdimuon:corPt"), float, doc="CVH refit dimuon pt", precision=-1),
        cvhEta = ExtVar(cms.InputTag("trackrefitdimuon:corEta"), float, doc="CVH refit dimuon eta", precision=12),
        cvhPhi = ExtVar(cms.InputTag("trackrefitdimuon:corPhi"), float, doc="CVH refit dimuon phi", precision=12),
        cvhMuPlusPt = ExtVar(cms.InputTag("trackrefitdimuon:muPlusPt"), float, doc="CVH refit mu+ pt", precision=-1),
        cvhMuPlusEta = ExtVar(cms.InputTag("trackrefitdimuon:muPlusEta"), float, doc="CVH refit mu+ eta", precision=12),
        cvhMuPlusPhi = ExtVar(cms.InputTag("trackrefitdimuon:muPlusPhi"), float, doc="CVH refit mu+ phi", precision=12),
        cvhMuMinusPt = ExtVar(cms.InputTag("trackrefitdimuon:muMinusPt"), float, doc="CVH refit mu- pt", precision=-1),
        cvhMuMinusEta = ExtVar(cms.InputTag("trackrefitdimuon:muMinusEta"), float, doc="CVH refit mu- eta", precision=12),
        cvhMuMinusPhi = ExtVar(cms.InputTag("trackrefitdimuon:muMinusPhi"), float, doc="CVH refit mu- phi", precision=12),
        cvhEdmval = ExtVar(cms.InputTag("trackrefitdimuon:edmval"), float, doc="CVH refit estimated distance to minimum", precision=10),
    ),
)

# Global-fit payload (large): attached only when the refit's fillGradsFactored
# flag is on. Full precision (the downstream global fit needs it).
dimuonVecVarsTable = cms.EDProducer("FlattenedCandValueMapVectorTableProducer",
    name = cms.string("Dimuon"),
    src = cms.InputTag("diMuonTrackVertexCandidates"),
    cut = cms.string(""),
    doc = cms.string("CVH two-track refit global-fit payload"),
    variables = cms.PSet(
        cvhGlobalIdxs = ExtVar(cms.InputTag("trackrefitdimuon:globalIdxs"), "std::vector<int>", doc="global correction-parameter indices", precision=16),
        cvhJacRefMuPlus = ExtVar(cms.InputTag("trackrefitdimuon:jacRefMuPlus"), "std::vector<float>", doc="d(mu+ refParms)/d(globalparms)", precision=-1),
        cvhJacRefMuMinus = ExtVar(cms.InputTag("trackrefitdimuon:jacRefMuMinus"), "std::vector<float>", doc="d(mu- refParms)/d(globalparms)", precision=-1),
        cvhJacMass = ExtVar(cms.InputTag("trackrefitdimuon:jacMass"), "std::vector<float>", doc="d(mass)/d(globalparms)", precision=-1),
        cvhHessFactor = ExtVar(cms.InputTag("trackrefitdimuon:hessFactor"), "std::vector<float>", doc="factored Hessian B (row-major nRank x nParms)", precision=-1),
    )
)

muonExternalVecVarsTable = cms.EDProducer("FlattenedCandValueMapVectorTableProducer",
    name = cms.string(muonTable.name.value()),
    src = muonTable.src,
    cut = muonTable.cut,
    doc = muonTable.doc,
    variables = cms.PSet(
        # can you declare a max number of bits here?  technically for the moment this needs 16 bits, but might eventually need 17 or 18
        cvhmergedGlobalIdxs = ExtVar(cms.InputTag("trackrefit:globalIdxs"), "std::vector<int>", doc="Indices for correction parameters", precision = 16),
        # optimal precision tbd, but presumably can work the same way as for scalar floats
        cvhJacRef = ExtVar(cms.InputTag("trackrefit:jacRef"), "std::vector<float>", doc="jacobian for corrections", precision = 12),
        # optimal precision tbd, but presumably can work the same way as for scalar floats
        cvhMomCov = ExtVar(cms.InputTag("trackrefit:momCov"), "std::vector<float>", doc="covariance matrix for qop, lambda, phi", precision = 12),
    )
)

def _cvhScalarVars(extVars, tag, suffix, note):
    """Add the 7 scalar cvh<suffix>* ExtVars from a refit producer `tag`."""
    setattr(extVars, "cvh%sPt" % suffix, ExtVar(cms.InputTag(tag + ":corPt"), float, doc="Refitted track pt" + note, precision=-1))
    setattr(extVars, "cvh%sEta" % suffix, ExtVar(cms.InputTag(tag + ":corEta"), float, doc="Refitted track eta" + note, precision=12))
    setattr(extVars, "cvh%sPhi" % suffix, ExtVar(cms.InputTag(tag + ":corPhi"), float, doc="Refitted track phi" + note, precision=12))
    setattr(extVars, "cvh%sCharge" % suffix, ExtVar(cms.InputTag(tag + ":corCharge"), int, doc="Refitted track charge" + note))
    setattr(extVars, "cvh%sDxy" % suffix, ExtVar(cms.InputTag(tag + ":corDxy"), float, doc="Refitted track dxy (d0) wrt beamspot" + note, precision=12))
    setattr(extVars, "cvh%sDz" % suffix, ExtVar(cms.InputTag(tag + ":corDz"), float, doc="Refitted track dz (z0) wrt beamspot" + note, precision=12))
    setattr(extVars, "cvh%sEdmval" % suffix, ExtVar(cms.InputTag(tag + ":edmval"), float, doc="Refitted estimated distance to minimum" + note, precision=10))
    setattr(extVars, "cvh%sNValidHits" % suffix, ExtVar(cms.InputTag(tag + ":nValidHits"), int, doc="Number of valid hits in refit" + note))
    setattr(extVars, "cvh%sNValidPixelHits" % suffix, ExtVar(cms.InputTag(tag + ":nValidPixelHits"), int, doc="Number of valid pixel hits in refit" + note))


def nanoAOD_addCvhMuonBranches(process, initFile=None, isMC=False, useScalarPot3D=None):
    """Attach the CVH-refit muon branches (Muon_cvh*) to the muon table.

    Wires the single-muon-track CVH refit into the NanoAOD muon tables:
    defines process.tracksfrommuons (inner tracks + pat::Muon association) and
    process.trackrefit (nominal, real geometry, data + MC), adds them and the
    vector-value-map table to muonTablesTask, and sets up the 3D scalar-
    potential field + Geant4e propagator the refit needs
    (nano_cff.setup3DFieldForRefit).

    isMC=True also wires the MC-only variants -- process.trackrefitideal
    (ideal geometry) and process.trackrefitbs (beamspot constraint) -- plus
    process.mergedGlobalIdxs (union of the nominal+ideal correction-parameter
    indices). All three refits share the one EventSetup G4 master; running them
    together in one job is what Stage A2 unblocked. The nominal+ideal branch
    set reproduces the validated 10_6-tip contract (cvh*/cvhideal* +
    cvhmergedGlobalIdxs from mergedGlobalIdxs); the beamspot set (cvhbs*) is
    added as a self-contained parallel set with its own cvhbsGlobalIdxs.

    initFile: scalar-potential coefficient dump (mfs/dump_coeffs_for_cmssw.py
    output). If None, setup3DFieldForRefit falls back to CVH_SCALARPOT_INITFILE
    / its built-in default.

    useScalarPot3D: baseline field for the refit. None (default) means data uses
    the accurate ScalarPot3D map and MC keeps the DEFAULT CMSSW field (consistent
    with the field the simulation used); pass True/False to override.
    """
    from PhysicsTools.NanoAOD.nano_cff import setup3DFieldForRefit
    if useScalarPot3D is None:
        useScalarPot3D = not isMC
    extVars = process.muonTable.externalVariables
    vecVars = muonExternalVecVarsTable.variables

    # Nominal (data + MC).
    _cvhScalarVars(extVars, "trackrefit", "", "")
    process.tracksfrommuons = tracksfrommuons
    process.trackrefit = trackrefit
    _producers = [process.tracksfrommuons, process.trackrefit]

    if isMC:
        # MC-only ideal + beamspot refits + merged nominal/ideal indices.
        process.trackrefitideal = trackrefitideal
        process.trackrefitbs = trackrefitbs
        process.mergedGlobalIdxs = mergedGlobalIdxs
        _producers += [process.trackrefitideal, process.trackrefitbs, process.mergedGlobalIdxs]
        _cvhScalarVars(extVars, "trackrefitideal", "ideal", " (ideal geometry)")
        _cvhScalarVars(extVars, "trackrefitbs", "bs", " (beamspot constraint)")
        # Merged nominal+ideal index vector (10_6-tip contract).
        vecVars.cvhmergedGlobalIdxs = ExtVar(cms.InputTag("mergedGlobalIdxs"), "std::vector<int>", doc="Indices for correction parameters (merged nominal+ideal)", precision=16)
        vecVars.cvhidealJacRef = ExtVar(cms.InputTag("trackrefitideal:jacRef"), "std::vector<float>", doc="jacobian for corrections (ideal geometry)", precision=12)
        vecVars.cvhidealMomCov = ExtVar(cms.InputTag("trackrefitideal:momCov"), "std::vector<float>", doc="covariance matrix for qop, lambda, phi (ideal geometry)", precision=12)
        # Beamspot-constrained set (self-contained: own index vector).
        vecVars.cvhbsGlobalIdxs = ExtVar(cms.InputTag("trackrefitbs:globalIdxs"), "std::vector<int>", doc="Indices for correction parameters (beamspot constraint)", precision=16)
        vecVars.cvhbsJacRef = ExtVar(cms.InputTag("trackrefitbs:jacRef"), "std::vector<float>", doc="jacobian for corrections (beamspot constraint)", precision=12)
        vecVars.cvhbsMomCov = ExtVar(cms.InputTag("trackrefitbs:momCov"), "std::vector<float>", doc="covariance matrix for qop, lambda, phi (beamspot constraint)", precision=12)

    process.muonExternalVecVarsTable = muonExternalVecVarsTable
    _producers.append(process.muonExternalVecVarsTable)

    # Dimuon (two-track) CVH refit + Dimuon table (data + MC). Shares the one
    # G4 master with the single-track refit.
    process.diMuonTrackVertexCandidates = diMuonTrackVertexCandidates
    process.trackrefitdimuon = trackrefitdimuon
    process.dimuonTable = dimuonTable
    _producers += [process.diMuonTrackVertexCandidates, process.trackrefitdimuon, process.dimuonTable]
    if trackrefitdimuon.fillGradsFactored.value():
        process.dimuonVecVarsTable = dimuonVecVarsTable
        _producers.append(process.dimuonVecVarsTable)

    process.muonTablesTask.add(*_producers)

    # Field + Geant4e propagator + shared G4 master. Data -> ScalarPot3D map;
    # MC -> default (sim-consistent) field.
    setup3DFieldForRefit(process, initFile=initFile, useScalarPot3D=useScalarPot3D)
    return process



# Revert back to AK4 CHS jets for Run 2
run2_muon.toModify(
    ptRatioRelForMu,srcJet="updatedJets"
)


muonsMCMatchForTable = cms.EDProducer("MCMatcher",       # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = muonTable.src,                         # final reco collection
    matched     = cms.InputTag("finalGenParticles"),     # final mc-truth particle collection
    mcPdgId     = cms.vint32(13),               # one or more PDG ID (13 = mu); absolute values (see below)
    checkCharge = cms.bool(False),              # True = require RECO and MC objects to have the same charge
    mcStatus    = cms.vint32(1),                # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.3),              # Minimum deltaR for the match
    maxDPtRel   = cms.double(0.5),              # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),     # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True),    # False = just match input in order; True = pick lowest deltaR pair first
)

muonMCTable = cms.EDProducer("CandMCMatchTableProducer",
    src     = muonTable.src,
    mcMap   = cms.InputTag("muonsMCMatchForTable"),
    objName = muonTable.name,
    objType = muonTable.name, #cms.string("Muon"),
    branchName = cms.string("genPart"),
    docString = cms.string("MC matching to status==1 muons"),
)

muonTask = cms.Task(slimmedMuonsUpdated,isoForMu,ptRatioRelForMu,slimmedMuonsWithUserData,finalMuons,finalLooseMuons)
muonMCTask = cms.Task(muonsMCMatchForTable,muonMCTable)
muonTablesTask = cms.Task(muonPROMPTMVA,muonMVALowPt,muonBSConstrain,muonTable,muonMVAID,muonPNetVariables,muonPNetScores)

