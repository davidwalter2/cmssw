# Muon tag-and-probe NanoAOD (NANOTP) for CMSSW 15_0, ported from
# WmassNanoProd_10_6_26 (nanoTP_cff.py + generalTracks_cff.py +
# StandAloneMuonMerger). One PAT+NANO job from AOD:
#
#   cmsDriver.py NANO --step PAT,NANO --era Run2_2016 --procModifiers run2_miniAOD_UL \
#       --customise PhysicsTools/NanoAOD/nanoTP_cff.customizeNANOTP[,PhysicsTools/NanoAOD/nano_cff.nanoGenWmassCustomize]
#
# Deliberately WITHOUT the run2_nanoAOD_106Xv2 modifier: that modifier
# re-clusters PUPPI "from MiniAOD" and re-wires the tau/PNet chain, which is
# circular when PAT runs in the same process. It only touches jets, taus, MET,
# low-pT electrons and the muon selector recomputation (done below), none of
# which this nano keeps.
#
# Layout follows the 10_6 T&P sequence rather than the stock nano: the stock
# PATObjectCrossLinker consumes jets, electrons, photons and taus and would drag
# the whole 15_0 PAT jet/tau/egamma machinery over the 10_6 AOD. Here the Muon
# table is keyed to a plain PATMuonSelector, and only the tables listed in
# _keep* stay scheduled -- producers in the nano tasks run only when a scheduled
# table consumes them.
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, ExtVar
from PhysicsTools.NanoAOD.generalTracks_cff import generalTrackTable, standaloneMuonTable, standaloneMuonUpdatedAtVtxTable, mergedStandaloneMuonTable
from PhysicsTools.NanoAOD.standAloneMuonMerger_cfi import mergedStandAloneMuons
from PhysicsTools.NanoAOD.muonvtxagnosticiso_cff import nanoAOD_addVtxAgnosticIso
from PhysicsTools.NanoAOD.nano_cff import nanoAOD_wmassMuonVariables

# table tasks the T&P nano never had (10_6: muons, vertices, isolated tracks,
# tracks, standalone muons, trigger objects, L1 bits; MC: gen particles,
# weights, gen vertex, PU, muon gen match)
_dropTablesCommon = ("jetPuppiTablesTask", "jetPuppiForMETTask", "jetAK8TablesTask", "jetConstituentsTablesTask",
                     "tauTablesTask", "boostedTauTablesTask", "electronTablesTask",
                     "lowPtElectronTablesTask", "photonTablesTask", "metTablesTask",
                     "fsrTablesTask", "softActivityTablesTask", "extraFlagsTableTask")
_dropTablesMC = ("jetMCTask", "electronMCTask", "lowPtElectronMCTask", "photonMCTask",
                 "tauMCTask", "boostedTauMCTask", "metMCTable", "ttbarCatMCProducersTask",
                 "ttbarCategoryTableTask", "tauSpinnerTableTask", "genProtonTablesTask")
_dropTablesData = ("protonTable", "multiRPTable", "singleRPTable")
# Muon columns that only make sense with the jet cross-link / jet user data
_dropMuonVars = ("jetIdx", "jetRelIso", "jetPtRelv2", "jetDF", "jetNDauCharged")
_dropMuonExtVars = ("promptMVA", "mvaLowPt", "pnScore_prompt", "pnScore_heavy", "pnScore_light", "pnScore_tau", "fsrPhotonIdx")
_dropMuonTableModules = ("muonPROMPTMVA", "muonMVALowPt", "muonPNetVariables", "muonPNetScores")

def _removeFromTask(process, taskName, members):
    task = getattr(process, taskName, None)
    if task is None:
        return
    for m in members:
        if hasattr(process, m) and task.contains(getattr(process, m)):
            task.remove(getattr(process, m))

def _retargetLinkedObjects(process):
    """Every consumer of linkedObjects:muons -> linkedMuons, :vertices -> slimmedSecondaryVertices."""
    def fix(tag):
        if isinstance(tag, cms.InputTag) and tag.getModuleLabel() == "linkedObjects":
            if tag.getProductInstanceLabel() == "muons":
                return cms.InputTag("linkedMuons")
            if tag.getProductInstanceLabel() == "vertices":
                return cms.InputTag("slimmedSecondaryVertices")
        return None
    def walk(pset):
        for name in pset.parameterNames_():
            par = getattr(pset, name)
            new = fix(par)
            if new is not None:
                setattr(pset, name, new)
            elif isinstance(par, cms.PSet):
                walk(par)
            elif isinstance(par, cms.VInputTag):
                tags = [fix(t) or t for t in par]
                setattr(pset, name, cms.VInputTag(*tags))
    for coll in (process.producers_(), process.filters_(), process.analyzers_()):
        for mod in coll.values():
            walk(mod)

def customizeNANOTP(process):
    # --- muon selection: keep standalone-only muons down to standalone pT > 15 GeV
    # through both the PAT selector and the nano final selection (10_6 logic)
    passStandalone = "(standAloneMuon().isNonnull() && standAloneMuon().pt() > 15)"
    process.selectedPatMuons.cut = cms.string("||".join([passStandalone, process.selectedPatMuons.cut.value()]))
    process.finalMuons.cut = cms.string("||".join([passStandalone, process.finalMuons.cut.value()]))

    # --- no rekeying of TrackExtra refs to the reduced collections: the
    # *ExtraIdx columns index the AOD collections, so muons, standalone tracks
    # and general tracks can be cross-referenced
    process.slimmedMuons.trackExtraAssocs = cms.VInputTag()

    # --- run2_nanoAOD_106Xv2 (not used, see above) recomputes the muon
    # selectors and the Run 3 soft MVA on 106X input; keep that
    process.slimmedMuonsUpdated.recomputeMuonBasicSelectors = True
    process.slimmedMuonsUpdated.recomputeSoftMuonMvaRun3 = True

    # --- no jets: cut the lepton-jet user data out of the muon chain
    process.muonTask.remove(process.ptRatioRelForMu)
    for f in ("ptRatio", "ptRel", "jetNDauChargedMVASel"):
        delattr(process.slimmedMuonsWithUserData.userFloats, f)
    delattr(process.slimmedMuonsWithUserData.userCands, "jetForLepJetVar")
    for v in _dropMuonVars:
        if hasattr(process.muonTable.variables, v):
            delattr(process.muonTable.variables, v)
    for v in _dropMuonExtVars:
        if hasattr(process.muonTable.externalVariables, v):
            delattr(process.muonTable.externalVariables, v)
    _removeFromTask(process, "muonTablesTask", _dropMuonTableModules)

    # --- muon collection for the tables: plain selector instead of the
    # cross-linker (which would consume jets, electrons, photons and taus)
    process.linkedMuons = cms.EDFilter("PATMuonSelector", src = process.finalMuons.src, cut = process.finalMuons.cut)
    process.muonTask.add(process.linkedMuons)
    _retargetLinkedObjects(process)

    # --- isolated tracks cleaned against loose muons only (no electrons here),
    # and without the impact-parameter requirement (10_6 1230c724004, as in
    # nanoAOD_wmassContent)
    process.finalIsolatedTracks.finalLeptons = cms.VInputTag("finalLooseMuons")
    process.finalIsolatedTracks.cut = cms.string(
        "((pt>5 && (abs(pdgId) == 11 || abs(pdgId) == 13)) || pt > 10) && (abs(pdgId) < 15 || abs(eta) < 2.5) && "
        "((pfIsolationDR03().chargedHadronIso < 5 && pt < 25) || pfIsolationDR03().chargedHadronIso/pt < 0.2)")

    # --- Muon table: 10_6 W-mass columns + T&P cross-reference columns
    nanoAOD_wmassMuonVariables(process)
    v = process.muonTable.variables
    v.standaloneExtraIdx = Var('? standAloneMuon().isNonnull() ? standAloneMuon().extra().key() : -99', 'int', precision=-1, doc='Index of the standalone-muon TrackExtra in the original collection')
    v.innerTrackExtraIdx = Var('? innerTrack().isNonnull() ? innerTrack().extra().key() : -99', 'int', precision=-1, doc='Index of the innerTrack TrackExtra in the original collection')
    v.vx = Var('vx', 'float', precision=-1, doc='Muon X position')
    v.vy = Var('vy', 'float', precision=-1, doc='Muon Y position')
    v.vz = Var('vz', 'float', precision=-1, doc='Muon Z position')

    # --- standalone-muon disambiguation (GlobalMuonProducer logic) keyed to the
    # nano muon collection, plus the flag on the Muon table
    process.mergedStandAloneMuons = mergedStandAloneMuons.clone(muons = process.muonTable.src)
    process.muonTable.externalVariables = cms.PSet(process.muonTable.externalVariables,
        isStandAloneUpdatedAtVtx = ExtVar(cms.InputTag("mergedStandAloneMuons:muonUpdatedAtVtx"), bool, doc="is standalone muon track updated at vertex"),
    )
    nanoAOD_addVtxAgnosticIso(process)

    # --- track tables: general tracks (pT > 8), the two standalone collections
    # and the merged one
    process.generalTrackTable = generalTrackTable.clone()
    process.standaloneMuonTable = standaloneMuonTable.clone()
    process.standaloneMuonUpdatedAtVtxTable = standaloneMuonUpdatedAtVtxTable.clone()
    process.mergedStandaloneMuonTable = mergedStandaloneMuonTable.clone()
    process.tnpTrackTablesTask = cms.Task(process.mergedStandAloneMuons, process.generalTrackTable,
                                          process.standaloneMuonTable, process.standaloneMuonUpdatedAtVtxTable,
                                          process.mergedStandaloneMuonTable)
    process.nanoTableTaskCommon.add(process.tnpTrackTablesTask)

    # --- unschedule what the T&P nano never had
    _removeFromTask(process, "nanoTableTaskCommon", _dropTablesCommon)
    _removeFromTask(process, "nanoTableTaskFS", _dropTablesMC)
    _removeFromTask(process, "protonTablesTask", _dropTablesData)
    return process

def customizeNANOTPLowPU(process):
    """The 2017 low-PU run (2017H) T&P nano on the UL re-reco AOD: the T&P
    content plus the muon trigger objects of the run's HI-style menu (the
    Mu17 filter, as in nanoAOD_wmassLowPU); electrons are not part of the
    T&P nano, so the Ele20/Ele17HI selection is not carried."""
    from PhysicsTools.NanoAOD.triggerObjects_cff import mksel
    process = customizeNANOTP(process)
    process.triggerObjectTable.selections.Muon = cms.PSet(
        id = cms.int32(13),
        sel = cms.string("type(83) && pt > 5 && (coll('hltIterL3MuonCandidates') || (pt > 45 && coll('hltHighPtTkMuonCands')) || (pt > 95 && coll('hltOldL3MuonCandidates')))"),
        l1seed = cms.string("type(-81)"), l1deltaR = cms.double(0.5),
        l2seed = cms.string("type(83) && coll('hltL2MuonCandidates')"), l2deltaR = cms.double(0.3),
        skipObjectsNotPassingQualityBits = cms.bool(True),
        qualityBits = cms.VPSet(
            mksel("filter('hltL3fL1sMu10lqL1f0L2f10L3Filtered17')", "Mu17"),
        ),
    )
    return process
