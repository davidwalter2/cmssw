import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.nano_eras_cff import *
from PhysicsTools.NanoAOD.jetsAK4_CHS_cff import *
from PhysicsTools.NanoAOD.jetsAK4_Puppi_cff import *
from PhysicsTools.NanoAOD.jetsAK8_cff import *
from PhysicsTools.NanoAOD.jetMC_cff import *
from PhysicsTools.NanoAOD.jetConstituents_cff import *
from PhysicsTools.NanoAOD.muons_cff import *
from PhysicsTools.NanoAOD.taus_cff import *
from PhysicsTools.NanoAOD.boostedTaus_cff import *
from PhysicsTools.NanoAOD.electrons_cff import *
from PhysicsTools.NanoAOD.lowPtElectrons_cff import *
from PhysicsTools.NanoAOD.photons_cff import *
from PhysicsTools.NanoAOD.globals_cff import *
from PhysicsTools.NanoAOD.extraflags_cff import *
from PhysicsTools.NanoAOD.ttbarCategorization_cff import *
from PhysicsTools.NanoAOD.genparticles_cff import *
from PhysicsTools.NanoAOD.particlelevel_cff import *
from PhysicsTools.NanoAOD.genWeightsTable_cfi import *
from PhysicsTools.NanoAOD.tauSpinnerTable_cfi import *
from PhysicsTools.NanoAOD.genVertex_cff import *
from PhysicsTools.NanoAOD.vertices_cff import *
from PhysicsTools.NanoAOD.met_cff import *
from PhysicsTools.NanoAOD.triggerObjects_cff import *
from PhysicsTools.NanoAOD.isotracks_cff import *
from PhysicsTools.NanoAOD.protons_cff import *
from PhysicsTools.NanoAOD.NanoAODEDMEventContent_cff import *
from PhysicsTools.NanoAOD.fsrPhotons_cff import *
from PhysicsTools.NanoAOD.softActivity_cff import *

nanoMetadata = cms.EDProducer("UniqueStringProducer",
    strings = cms.PSet(
        tag = cms.string("untagged"),
    )
)

linkedObjects = cms.EDProducer("PATObjectCrossLinker",
   jets=cms.InputTag("finalJetsPuppi"),
   muons=cms.InputTag("finalMuons"),
   electrons=cms.InputTag("finalElectrons"),
   lowPtElectrons=cms.InputTag("finalLowPtElectrons"),
   taus=cms.InputTag("finalTaus"),
   boostedTaus=cms.InputTag("finalBoostedTaus"),
   photons=cms.InputTag("finalPhotons"),
   vertices=cms.InputTag("slimmedSecondaryVertices")
)

from PhysicsTools.NanoAOD.lhcInfoProducer_cfi import lhcInfoProducer
lhcInfoTable = lhcInfoProducer.clone()
(~run3_common).toModify(
    lhcInfoTable, useNewLHCInfo=False
)

nanoTableTaskCommon = cms.Task(
    cms.Task(nanoMetadata),
    jetPuppiTask, jetPuppiForMETTask, jetAK8Task, jetConstituentsTask,
    extraFlagsProducersTask, muonTask, tauTask, boostedTauTask,
    electronTask , lowPtElectronTask, photonTask,
    vertexTask, isoTrackTask, jetAK8LepTask,  # must be after all the leptons
    softActivityTask,
    cms.Task(linkedObjects),
    jetPuppiTablesTask, jetAK8TablesTask, jetConstituentsTablesTask,
    muonTablesTask, fsrTablesTask, tauTablesTask, boostedTauTablesTask,
    electronTablesTask, lowPtElectronTablesTask, photonTablesTask,
    globalTablesTask, vertexTablesTask, metTablesTask, extraFlagsTableTask,
    isoTrackTablesTask,softActivityTablesTask
)

(run2_muon | run2_egamma).toReplaceWith(
    nanoTableTaskCommon,
    nanoTableTaskCommon.copyAndAdd(chsJetUpdateTask)
)

nanoSequenceCommon = cms.Sequence(nanoTableTaskCommon)

nanoSequenceOnlyFullSim = cms.Sequence(triggerObjectTablesTask)
nanoSequenceOnlyData = cms.Sequence(cms.Sequence(protonTablesTask) + lhcInfoTable)

nanoSequence = cms.Sequence(nanoSequenceCommon + nanoSequenceOnlyData + nanoSequenceOnlyFullSim)

nanoTableTaskFS = cms.Task(
    genParticleTask, particleLevelTask, jetMCTask, muonMCTask, electronMCTask, lowPtElectronMCTask, photonMCTask,
    tauMCTask, boostedTauMCTask,
    metMCTable, ttbarCatMCProducersTask, globalTablesMCTask, ttbarCategoryTableTask,
    genWeightsTableTask, genVertexTablesTask, genParticleTablesTask, genProtonTablesTask, particleLevelTablesTask, tauSpinnerTableTask
)

nanoSequenceFS = cms.Sequence(nanoSequenceCommon + cms.Sequence(nanoTableTaskFS))

# GenVertex only stored in newer MiniAOD
nanoSequenceMC = nanoSequenceFS.copy()
nanoSequenceMC.insert(nanoSequenceFS.index(nanoSequenceCommon)+1,nanoSequenceOnlyFullSim)


def _fixPNetInputCollection(process):
    # fix circular module dependency in ParticleNetFromMiniAOD TagInfos when slimmedTaus is updated
    if hasattr(process, 'slimmedTaus'):
        for mod in process.producers.keys():
            if 'ParticleNetFromMiniAOD' in mod and 'TagInfos' in mod:
                getattr(process, mod).taus = 'slimmedTaus::@skipCurrentProcess'


# modifier which adds new tauIDs
import RecoTauTag.RecoTau.tools.runTauIdMVA as tauIdConfig
def nanoAOD_addTauIds(process, idsToRun=[], addPNetCHS=False, addUParTPuppi=False):
    originalTauName = 'slimmedTaus::@skipCurrentProcess'
    updatedTauName = None

    if idsToRun:  # no-empty list of tauIDs to run
        updatedTauName = 'slimmedTausUpdated'
        tauIdEmbedder = tauIdConfig.TauIDEmbedder(process, debug=False,
                                                  originalTauName=originalTauName,
                                                  updatedTauName=updatedTauName,
                                                  postfix="ForNano",
                                                  toKeep=idsToRun)
        tauIdEmbedder.runTauID()
        process.tauTask.add(process.rerunMvaIsolationTaskForNano, getattr(process, updatedTauName))
        originalTauName = updatedTauName

    from PhysicsTools.PatAlgos.patTauHybridProducer_cfi import patTauHybridProducer
    if addPNetCHS:
        jetCollection = "updatedJets"
        TagName = "pfParticleNetFromMiniAODAK4CHSCentralJetTags"
        tag_prefix = "byUTagCHS"
        updatedTauName = originalTauName.split(':')[0] + 'WithPNetCHS'
        # PNet tagger used for CHS jets
        from RecoBTag.ONNXRuntime.pfParticleNetFromMiniAODAK4_cff import pfParticleNetFromMiniAODAK4CHSCentralJetTags
        Discriminators = [TagName + ":" + tag for tag in pfParticleNetFromMiniAODAK4CHSCentralJetTags.flav_names.value()]

        # Define "hybridTau" producer
        setattr(process, updatedTauName, patTauHybridProducer.clone(
            src=originalTauName,
            jetSource=jetCollection,
            dRMax=0.4,
            jetPtMin=15,
            jetEtaMax=2.5,
            UTagLabel=TagName,
            UTagScoreNames=Discriminators,
            tagPrefix=tag_prefix,
            tauScoreMin=-1,
            vsJetMin=0.05,
            checkTauScoreIsBest=False,
            chargeAssignmentProbMin=0.2,
            addGenJetMatch=False,
            genJetMatch=""
        ))
        process.tauTask.add(process.chsJetUpdateTask, getattr(process, updatedTauName))
        originalTauName = updatedTauName

    if addUParTPuppi:
        jetCollection = "updatedJetsPuppi"
        TagName = "pfUnifiedParticleTransformerAK4JetTags"
        tag_prefix = "byUTagPUPPI"
        updatedTauName = originalTauName.split(':')[0] + 'WithUParTPuppi'
        # Unified ParT Tagger used for PUPPI jets
        from RecoBTag.ONNXRuntime.pfUnifiedParticleTransformerAK4JetTags_cfi import pfUnifiedParticleTransformerAK4JetTags
        Discriminators = [TagName + ":" + tag for tag in pfUnifiedParticleTransformerAK4JetTags.flav_names.value()]

        # Define "hybridTau" producer
        setattr(process, updatedTauName, patTauHybridProducer.clone(
            src=originalTauName,
            jetSource=jetCollection,
            dRMax=0.4,
            jetPtMin=15,
            jetEtaMax=2.5,
            UTagLabel=TagName,
            UTagScoreNames=Discriminators,
            tagPrefix=tag_prefix,
            tauScoreMin=-1,
            vsJetMin=0.05,
            checkTauScoreIsBest=False,
            chargeAssignmentProbMin=0.2,
            addGenJetMatch=False,
            genJetMatch=""
        ))
        process.tauTask.add(getattr(process, updatedTauName))
        originalTauName = updatedTauName

    if updatedTauName is not None:
        process.slimmedTaus = getattr(process, updatedTauName).clone()
        process.tauTask.replace(getattr(process, updatedTauName), process.slimmedTaus)
        delattr(process, updatedTauName)
        _fixPNetInputCollection(process)

    return process


def nanoAOD_addBoostedTauIds(process, idsToRun=[]):
    if idsToRun:  # no-empty list of tauIDs to run
        boostedTauIdEmbedder = tauIdConfig.TauIDEmbedder(process, debug=False,
                                                         originalTauName="slimmedTausBoosted::@skipCurrentProcess",
                                                         updatedTauName="slimmedTausBoostedNewID",
                                                         postfix="BoostedForNano",
                                                         toKeep=idsToRun)
        boostedTauIdEmbedder.runTauID()

        process.slimmedTausBoosted = process.slimmedTausBoostedNewID.clone()
        del process.slimmedTausBoostedNewID
        process.boostedTauTask.add(process.rerunMvaIsolationTaskBoostedForNano, process.slimmedTausBoosted)

    return process


from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
def nanoAOD_activateVID(process):

    switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD, electronTask)
    for modname in electron_id_modules_WorkingPoints_nanoAOD.modules:
        setupAllVIDIdsInModule(process, modname, setupVIDElectronSelection)

    process.electronTask.add(process.egmGsfElectronIDTask)

    # do not call this to avoid resetting photon IDs in VID, if called before inside makePuppiesFromMiniAOD
    switchOnVIDPhotonIdProducer(process, DataFormat.MiniAOD, photonTask)
    for modname in photon_id_modules_WorkingPoints_nanoAOD.modules:
        setupAllVIDIdsInModule(process, modname, setupVIDPhotonSelection)

    process.photonTask.add(process.egmPhotonIDTask)

    return process


def nanoAOD_customizeCommon(process):

    process = nanoAOD_activateVID(process)

    nanoAOD_rePuppi_switch = cms.PSet(
        useExistingWeights=cms.bool(False),
        reclusterAK4MET=cms.bool(False),
        reclusterAK8=cms.bool(False),
    )

    # recompute Puppi weights, and remake AK4, AK8 Puppi jets and PuppiMET
    (run2_nanoAOD_106Xv2 | run3_nanoAOD_pre142X | nanoAOD_rePuppi).toModify(
        nanoAOD_rePuppi_switch, useExistingWeights=False, reclusterAK4MET=True, reclusterAK8=True
    )

    # PAT running in the same job (NANO from AOD, e.g. the muon tag-and-probe
    # nano): PAT builds the PUPPI jets and MET in this release already, and the
    # "from MiniAOD" re-clustering replaces slimmedJetsAK8 before the PAT
    # customisation (applySubstructure) runs on it -> skip it.
    if hasattr(process, "packedPFCandidates"):
        nanoAOD_rePuppi_switch.useExistingWeights = True
        nanoAOD_rePuppi_switch.reclusterAK4MET = False
        nanoAOD_rePuppi_switch.reclusterAK8 = False

    runOnMC = True
    if hasattr(process, "NANOEDMAODoutput") or hasattr(process, "NANOAODoutput"):
        runOnMC = False
    from PhysicsTools.PatAlgos.tools.puppiJetMETReclusteringFromMiniAOD_cff import puppiJetMETReclusterFromMiniAOD
    puppiJetMETReclusterFromMiniAOD(process,
                                    runOnMC=runOnMC,
                                    useExistingWeights=nanoAOD_rePuppi_switch.useExistingWeights.value(),
                                    reclusterAK4MET=nanoAOD_rePuppi_switch.reclusterAK4MET.value(),
                                    reclusterAK8=nanoAOD_rePuppi_switch.reclusterAK8.value(),
                                    )

    if not(nanoAOD_rePuppi_switch.useExistingWeights) and (nanoAOD_rePuppi_switch.reclusterAK4MET or nanoAOD_rePuppi_switch.reclusterAK8):
        process = UsePuppiWeightFromValueMapForPFCandTable(process)

    # This function is defined in jetsAK4_Puppi_cff.py
    process = nanoAOD_addDeepInfoAK4(process,
                                     addParticleNet=nanoAOD_addDeepInfoAK4_switch.nanoAOD_addParticleNet_switch,
                                     addRobustParTAK4=nanoAOD_addDeepInfoAK4_switch.nanoAOD_addRobustParTAK4Tag_switch,
                                     addUnifiedParTAK4=nanoAOD_addDeepInfoAK4_switch.nanoAOD_addUnifiedParTAK4Tag_switch
                                     )

    # Needs to run PNet on CHS jets to update the tau collections
    run2_nanoAOD_106Xv2.toModify(
        nanoAOD_addDeepInfoAK4CHS_switch, nanoAOD_addParticleNet_switch=True,
    )
    # This function is defined in jetsAK4_CHS_cff.py
    process = nanoAOD_addDeepInfoAK4CHS(
        process, addDeepBTag=nanoAOD_addDeepInfoAK4CHS_switch.nanoAOD_addDeepBTag_switch,
        addDeepFlavour=nanoAOD_addDeepInfoAK4CHS_switch.nanoAOD_addDeepFlavourTag_switch,
        addParticleNet=nanoAOD_addDeepInfoAK4CHS_switch.nanoAOD_addParticleNet_switch,
        addRobustParTAK4=nanoAOD_addDeepInfoAK4CHS_switch.nanoAOD_addRobustParTAK4Tag_switch,
        addUnifiedParTAK4=nanoAOD_addDeepInfoAK4CHS_switch.nanoAOD_addUnifiedParTAK4Tag_switch)

    # This function is defined in jetsAK8_cff.py
    process = nanoAOD_addDeepInfoAK8(
        process, addDeepBTag=nanoAOD_addDeepInfoAK8_switch.nanoAOD_addDeepBTag_switch,
        addDeepBoostedJet=nanoAOD_addDeepInfoAK8_switch.nanoAOD_addDeepBoostedJet_switch,
        addDeepDoubleX=nanoAOD_addDeepInfoAK8_switch.nanoAOD_addDeepDoubleX_switch,
        addDeepDoubleXV2=nanoAOD_addDeepInfoAK8_switch.nanoAOD_addDeepDoubleXV2_switch,
        addParticleNetMassLegacy=nanoAOD_addDeepInfoAK8_switch.nanoAOD_addParticleNetMassLegacy_switch,
        addParticleNetLegacy=nanoAOD_addDeepInfoAK8_switch.nanoAOD_addParticleNetLegacy_switch,
        addParticleNet=nanoAOD_addDeepInfoAK8_switch.nanoAOD_addParticleNet_switch,
        addGlobalParT=nanoAOD_addDeepInfoAK8_switch.nanoAOD_addGlobalParT_switch,
        jecPayload=nanoAOD_addDeepInfoAK8_switch.jecPayload)

    nanoAOD_tau_switch = cms.PSet(
        idsToAdd=cms.vstring(),
        addPNetCHS=cms.bool(False),
        addUParTPuppi=cms.bool(False)
    )
    (run2_nanoAOD_106Xv2).toModify(
        nanoAOD_tau_switch, idsToAdd=["deepTau2018v2p5"]
    )
    (run2_nanoAOD_106Xv2 | run3_nanoAOD_pre142X).toModify(
        nanoAOD_tau_switch, addPNetCHS=True, addUParTPuppi=True,
    )
    nanoAOD_addTauIds(process,
                      idsToRun=nanoAOD_tau_switch.idsToAdd.value(),
                      addPNetCHS=nanoAOD_tau_switch.addPNetCHS.value(),
                      addUParTPuppi=nanoAOD_tau_switch.addUParTPuppi.value(),
                      )

    nanoAOD_boostedTau_switch = cms.PSet(
        idsToAdd=cms.vstring()
    )
    run2_nanoAOD_106Xv2.toModify(
        nanoAOD_boostedTau_switch,
        idsToAdd=["mvaIso", "mvaIsoNewDM", "mvaIsoDR0p3", "againstEle", "boostedDeepTauRunIIv2p0"])
    run3_nanoAOD_pre142X.toModify(
        nanoAOD_boostedTau_switch, idsToAdd=["boostedDeepTauRunIIv2p0"]
    )
    nanoAOD_addBoostedTauIds(process, nanoAOD_boostedTau_switch.idsToAdd.value())
    
    from PhysicsTools.NanoAOD.leptonTimeLifeInfo_common_cff import addTimeLifeInfoBase
    process = addTimeLifeInfoBase(process)
    
    process = nanoAOD_refineFastSim_puppiJet(process)
    process = nanoAOD_refineFastSim_bTagDeepFlav(process, nanoAOD_addDeepInfoAK4CHS_switch.nanoAOD_addDeepFlavourTag_switch)

    return process

###increasing the precision of selected GenParticles.
import os

# Default scalar-potential coefficient dump for the CVH refit field model.
# The scalar-potential coefficient dump for the CVH refit. Shipped in the CMSSW
# area as MagneticField/ParametrizedEngine/data/... so that a CRAB sandbox
# carries it (CRAB packs the data/ directories); resolved through CMSSW_BASE,
# then CMSSW_RELEASE_BASE, the way edm::FileInPath does. Override with the
# CVH_SCALARPOT_INITFILE env var (cmsDriver --customise cannot pass function
# arguments).
_SCALARPOT_INITFILE_REL = "MagneticField/ParametrizedEngine/data/polyfit3d_full_coeffs_lmax18_cmsswnorm.txt"

def _resolveInPath(rel):
    for base in (os.environ.get("CMSSW_BASE"), os.environ.get("CMSSW_RELEASE_BASE")):
        if base and os.path.isfile(os.path.join(base, "src", rel)):
            return os.path.join(base, "src", rel)
    return rel

_DEFAULT_SCALARPOT_INITFILE = _resolveInPath(_SCALARPOT_INITFILE_REL)


# Correction model of the CVH refits in the NanoAOD, i.e. which global
# parameters the stored Jacobians (Muon_cvhJacRef, Dimuon payload) refer to:
#
#   "module"  one Bz offset (parmtype 6) and one material scale (parmtype 7) per
#             tracker module, dead modules included (their hitless surfaces stay
#             in the fit) -- the 10_6 / W-mass scheme. DATA baseline field: the
#             OPERA/TOSCA finite-element map CVH_OPERA_VERSION (default 170812,
#             the latest model) as the full 3D grid, without the tracker
#             parametrization (Analysis/HitAnalyzer/python/cvhOperaField.py; the
#             170812 tables must be on CMSSW_SEARCH_PATH, see TABLES_HELP there).
#             No scalar-potential or material-group Jacobians.
#   "global"  the scalar-potential field modes (parmtype 14) and the global
#             material groups (parmtype 15); DATA baseline field: ScalarPot3D.
#
# MC always keeps the default field (the one the simulation used).
CVH_CORRECTION_MODEL = "module"
CVH_OPERA_VERSION = "170812"
_CVH_REFITS = ("trackrefit", "trackrefitideal", "trackrefitbs", "trackrefitdimuon")


def _cvhSimGeometry(process):
    """The sim geometry (DDD) of the detector era of the process, for the G4
    world of the CVH refit: ("db", <GeometryFileRcd tag>) for Run 2 -- the
    XML blob the MC global tags of the era carry as "Extended" -- or
    ("xml", <XMLIdealGeometryESSource cfi>) for Run 3."""
    from Configuration.Eras.Modifier_phase1Pixel_cff import phase1Pixel
    from Configuration.Eras.Modifier_run2_HCAL_2018_cff import run2_HCAL_2018
    from Configuration.Eras.Modifier_run3_common_cff import run3_common
    from Configuration.Eras.Modifier_run3_egamma_2023_cff import run3_egamma_2023
    from Configuration.Eras.Modifier_stage2L1Trigger_2024_cff import stage2L1Trigger_2024
    from Configuration.Eras.Modifier_run3_SiPixel_2025_cff import run3_SiPixel_2025
    uses = process.isUsingModifier
    if not uses(phase1Pixel):
        return "db", "XMLFILE_Geometry_2016_81YV1_Extended2016_mc"
    if not uses(run3_common):
        if uses(run2_HCAL_2018):
            return "db", "XMLFILE_Geometry_101YV4_Extended2018_mc"
        return "db", "XMLFILE_Geometry_92YV5_Extended2017Plan1_mc"
    if uses(run3_SiPixel_2025):
        return "xml", "Geometry.CMSCommonData.cmsExtendedGeometry2025XML_cfi"
    if uses(stage2L1Trigger_2024):
        return "xml", "Geometry.CMSCommonData.cmsExtendedGeometry2024XML_cfi"
    if uses(run3_egamma_2023):
        return "xml", "Geometry.CMSCommonData.cmsExtendedGeometry2023XML_cfi"
    return "xml", "Geometry.CMSCommonData.cmsExtendedGeometry2021XML_cfi"


def setup3DFieldForRefit(process, initFile=None, useScalarPot3D=True, correctionModel=None):
    """Set up the Geant4e propagator + shared G4 master the CVH muon refit needs,
    choose the correction model of every refit maker, and choose the baseline
    magnetic field.

    correctionModel: "module" or "global" (default CVH_CORRECTION_MODEL), see
    the comment above CVH_CORRECTION_MODEL.

    useScalarPot3D=True (DATA): replace the baseline field with the data field of
    the model -- the OPERA 160812 3D grid for "module", ScalarPot3D for "global"
    -- routed into every CVH-side consumer (the propagator, the strip/pixel CPEs,
    the refit makers and the shared G4 master). The name is historical.
    useScalarPot3D=False (MC): keep the DEFAULT CMSSW field -- it must stay
    consistent with the field the simulation used (see CLAUDE.md, "custom field
    map for data only"); only the geometry / propagator / master are set up.

    geopro is loaded (for the propagator) but is NOT scheduled -- the shared
    CvhMasterThread ES product owns the G4 world.

    initFile ("global" only): scalar-potential coefficient dump (basis +
    coefficients). Falls back to CVH_SCALARPOT_INITFILE then
    _DEFAULT_SCALARPOT_INITFILE.
    """
    model = correctionModel or CVH_CORRECTION_MODEL
    if model not in ("module", "global"):
        raise RuntimeError("setup3DFieldForRefit: correctionModel must be 'module' or 'global'")
    if model == "global":
        if initFile is None:
            initFile = os.environ.get("CVH_SCALARPOT_INITFILE",
                                      _DEFAULT_SCALARPOT_INITFILE)
        if not initFile:
            raise RuntimeError(
                "setup3DFieldForRefit: no scalar-potential coefficient file "
                "(pass initFile= or set CVH_SCALARPOT_INITFILE)")

    # DDD sim geometry (DDCompactView on IdealGeometryRecord) needed by the
    # Geant4e propagator / CvhMaster G4 world. The stock NANO process only
    # provides the DB reco tracker geometry, not the DDD sim geometry, so load
    # it here (matches the standalone drivers). Additive: DDCompactView is a
    # different data type in IdealGeometryRecord than the DB GeometricDet.
    # The detector of the era: the data GTs carry no "Extended" sim geometry,
    # and the 2016 one must not be used for the Phase-1 tracker (2017 on).
    _geoKind, _geoSource = _cvhSimGeometry(process)
    if _geoKind == "db":
        process.load("Configuration.StandardSequences.GeometrySimDB_cff")
        process.GlobalTag.toGet.append(
            cms.PSet(
                record=cms.string("GeometryFileRcd"),
                tag=cms.string(_geoSource),
                label=cms.untracked.string("Extended"),
            )
        )
        if hasattr(process, "XMLFromDBSource"):
            process.XMLFromDBSource.label = cms.string("Extended")
    else:
        # Run 3: the DB only has DD4hep-format blobs from 2023 on; the release's
        # DDD XML of the same detector feeds the DDD G4 world
        process.load(_geoSource)

    process.load("TrackPropagation.Geant4e.geantRefit_cff")

    _refits = _CVH_REFITS
    # Correction model of every refit maker (data and MC alike: it defines the
    # Jacobians stored in the NanoAOD and is independent of the baseline field).
    for _refit in _refits:
        if not hasattr(process, _refit):
            continue
        _m = getattr(process, _refit)
        if model == "module":
            _m.perModuleBfield = cms.bool(True)
            _m.globalMaterialModel = cms.bool(False)
            _m.perStepFieldModes = cms.bool(False)
            _m.skipHitlessSurfaces = cms.bool(False)
            # the group file only feeds the global model (all k_init = 0, so
            # loading it would change nothing but the cost)
            _m.materialGroupsFile = cms.string("")
            _m.scalarPotentialInitFile = cms.string("")
        else:
            _m.perModuleBfield = cms.bool(False)
            _m.scalarPotentialInitFile = cms.string(initFile)

    # Baseline field. DATA: the model's data field, routed into every CVH-side
    # consumer. MC: keep the default CMSSW field (label "") so the refit stays
    # consistent with the field the simulation used -- no override.
    fieldlabel = ""
    if useScalarPot3D:
        if model == "module":
            # OPERA/TOSCA 3D grid, no tracker parametrization (cvhOperaField)
            from Analysis.HitAnalyzer.cvhOperaField import setupOpera3DField
            fieldlabel = setupOpera3DField(process, version=CVH_OPERA_VERSION, routeCPEs=False)
        else:
            from MagneticField.ParametrizedEngine.parametrizedMagneticField_ScalarPot3D_cfi \
                import ParametrizedMagneticFieldProducer as ScalarPot3DMagneticFieldProducer
            process.ScalarPot3DMagneticFieldProducer = ScalarPot3DMagneticFieldProducer.clone()
            process.ScalarPot3DMagneticFieldProducer.parameters.InitFile = initFile
            fieldlabel = "ScalarPot3DMf"
            process.ScalarPot3DMagneticFieldProducer.label = fieldlabel
        for _consumer in ("geopro", "Geant4ePropagator",
                          "stripCPEESProducer", "StripCPEfromTrackAngleESProducer",
                          "siPixelTemplateDBObjectESProducer", "templates") + _refits:
            if hasattr(process, _consumer):
                getattr(process, _consumer).MagneticFieldLabel = cms.string(fieldlabel)

    # Shared CVH Geant4 master (CvhMasterRecord), consumed by every CVH maker
    # via esConsumes. Uses the same baseline field as the refit -- the model's
    # data field for data, the default field (label "") for MC.
    from TrackPropagation.Geant4e.cvhMasterESProducer_cfi import cvhMasterESProducer
    process.cvhMasterESProducer = cvhMasterESProducer.clone()
    process.cvhMasterESProducer.MagneticFieldLabel = cms.string(fieldlabel)
    # The G4 world is built from the DDD sim geometry loaded above, also in the
    # Run 3 eras, whose dd4hep process modifier would switch the master and the
    # propagator to a cms::DDCompactView this job does not provide.
    process.cvhMasterESProducer.g4GeometryDD4hepSource = cms.bool(False)
    if hasattr(process, "geopro"):
        process.geopro.GeoFromDD4hep = cms.bool(False)

    # Activate the CVH-specific propagator path (custom fluctuation table) and
    # the per-leg forward/backward propagation choice.
    process.Geant4ePropagator.ForCVH = cms.bool(True)
    process.Geant4ePropagator.PropagationDirection = cms.string("anyDirection")

    # Per-stream CLHEP engine for each refit (Geant4 thread-local RNG). Seeds
    # kept < 9e8 (CLHEP HepJamesRandom range) and distinct per instance.
    if not hasattr(process, "RandomNumberGeneratorService"):
        process.load("Configuration.StandardSequences.Services_cff")
    for _refit, _seed in (("trackrefit", 123456789),
                          ("trackrefitideal", 223456789),
                          ("trackrefitbs", 323456789),
                          ("trackrefitdimuon", 423456789)):
        if hasattr(process, _refit):
            setattr(process.RandomNumberGeneratorService, _refit, cms.PSet(
                initialSeed=cms.untracked.uint32(_seed),
                engineName=cms.untracked.string("HepJamesRandom"),
            ))
    return process


def nanoAOD_addCvhMuon(process, initFile=None):
    """cmsDriver --customise entry point (DATA): nominal CVH muon refit only.

    Usage:
      cmsDriver.py ... --customise PhysicsTools/NanoAOD/nano_cff.nanoAOD_addCvhMuon
    (set the coefficient dump via the CVH_SCALARPOT_INITFILE env var).
    """
    from PhysicsTools.NanoAOD.muons_cff import nanoAOD_addCvhMuonBranches
    return nanoAOD_addCvhMuonBranches(process, initFile=initFile, isMC=False)


def nanoAOD_addCvhMuonMC(process, initFile=None):
    """cmsDriver --customise entry point (MC): nominal + ideal + beamspot refits.

    Adds the MC-only trackrefitideal / trackrefitbs variants and mergedGlobalIdxs
    on top of the nominal refit. All three share the one EventSetup G4 master.

    Usage:
      cmsDriver.py ... --customise PhysicsTools/NanoAOD/nano_cff.nanoAOD_addCvhMuonMC
    """
    from PhysicsTools.NanoAOD.muons_cff import nanoAOD_addCvhMuonBranches
    return nanoAOD_addCvhMuonBranches(process, initFile=initFile, isMC=True)


def nanoAOD_cvhPixelClassHits(process):
    """cmsDriver add-on, chained AFTER nanoAOD_addCvhMuon[MC]: re-admit the
    pixel edge / single-pixel hits in every CVH refit and emit their class-
    correction columns (parmtypes 16-21). The calibration applied to the
    resulting nano must be derived with the same setting (catalog +8640).

        --customise PhysicsTools/NanoAOD/nano_cff.nanoAOD_addCvhMuon,PhysicsTools/NanoAOD/nano_cff.nanoAOD_cvhPixelClassHits
    """
    from PhysicsTools.NanoAOD.muons_cff import cvhPixelClassHits
    return cvhPixelClassHits(process)


def nanoWmassGenCustomize(process):
    pdgSelection="?(abs(pdgId) == 11|| abs(pdgId)==13 || abs(pdgId)==15 ||abs(pdgId)== 12 || abs(pdgId)== 14 || abs(pdgId)== 16|| abs(pdgId)== 24|| pdgId== 23)"
    # Keep precision same as default RECO for selected particles
    ptPrecision="{}?{}:{}".format(pdgSelection, CandVars.pt.precision.value(),genParticleTable.variables.pt.precision.value())
    process.genParticleTable.variables.pt.precision=cms.string(ptPrecision)
    phiPrecision="{} ? {} : {}".format(pdgSelection, CandVars.phi.precision.value(), genParticleTable.variables.phi.precision.value())
    process.genParticleTable.variables.phi.precision=cms.string(phiPrecision)
    etaPrecision="{} ? {} : {}".format(pdgSelection, CandVars.eta.precision.value(), genParticleTable.variables.eta.precision.value())
    process.genParticleTable.variables.eta.precision=cms.string(etaPrecision)
    return process


### WMass customizations (ported from WmassNanoProd_10_6_26)
### Full precision for the GenParticles the W/Z analyses fit on.
def customizeGenLeptonPrecision(process):
    """Store pt/eta/phi at full precision (mantissa bits = -1) for leptons,
    neutrinos, top, W, Z and Higgs; everything else keeps the genParticleTable
    default. The stock nanoWmassGenCustomize above is NOT equivalent: it uses
    the CandVars precision (12 bits for eta/phi) and omits top and Higgs.
    """
    pdgSelection="?(abs(pdgId) == 11|| abs(pdgId)==13 || abs(pdgId)==15 ||abs(pdgId)== 12 || abs(pdgId)== 14 || abs(pdgId)== 16|| abs(pdgId)== 6|| abs(pdgId)== 24|| pdgId== 23|| pdgId== 25)"

    # Keep full precision for selected particles
    ptPrecision="{}?{}:{}".format(pdgSelection, -1, genParticleTable.variables.pt.precision.value())
    process.genParticleTable.variables.pt.precision=cms.string(ptPrecision)
    phiPrecision="{} ? {} : {}".format(pdgSelection, -1, genParticleTable.variables.phi.precision.value())
    process.genParticleTable.variables.phi.precision=cms.string(phiPrecision)
    etaPrecision="{} ? {} : {}".format(pdgSelection, -1, genParticleTable.variables.eta.precision.value())
    process.genParticleTable.variables.eta.precision=cms.string(etaPrecision)

    return process

def nanoGenWmassCustomize(process):
    """The 10_6 production entry point for MC (WMassNanoProduction makeNanoV9MC*.sh):
      cmsDriver.py ... --customise PhysicsTools/NanoAOD/nano_cff.nanoGenWmassCustomize

    Gen precision, the full LHE record, and the LHE weight tables in the grouped
    10_6 layout (LHEPdfWeightAltSetN, MEParamWeight[AltSetN], ...) that the
    analysis code reads; the stock LHEScaleWeight/LHEPdfWeight tables are
    replaced by the grouped ones, genWeight / PSWeight / LHEReweightingWeight
    stay from the stock producer.
    """
    process = customizeGenLeptonPrecision(process)

    process.lheInfoTable.storeAllLHEInfo = cms.bool(True)

    from PhysicsTools.NanoAOD.lheWeightGroups_cff import nanoAOD_addLHEWeightGroups
    process = nanoAOD_addLHEWeightGroups(process)

    return process

def nanoAOD_wmassContent(process):
    """WMass content of the 15_0 NanoAOD on UL2016 MiniAOD (data and MC), the
    non-CVH customisations of WmassNanoProd_10_6_26 that are not upstream:

    * drop the boosted-tau chain: boosted taus play no role in the W/Z
      analyses and the 10_6 nano did not have them. (With the UL 106X global
      tags the stock rerun of their MVA isolation would also fail -- the
      RecoTauTag GBRForest payloads only exist in the 150X tags the production
      uses since 2026-09-17.) Only unschedules the tasks; PATObjectCrossLinker
      accepts an empty boostedTaus tag.
    * the vertex-agnostic muon isolation (muonvtxagnosticiso_cff, 10_6 PR #31)
    * the extra Muon columns of the 10_6 table: inner-track algo, kink finder,
      tracker/pixel hit counts, the pfRelIso04 components, the standalone track
    * IsoTrack: no impact-parameter requirement (10_6 1230c724004)
    * the PV-robust PUPPI + DeepMET (deepMETPVRobust_cff, 10_6 PR #33):
      DeepMETPVRobust[NoPUPPI]_pt/phi, PVRobustIndex, PVMuonIndex

    Already upstream in 15_0, nothing to do: the pt > 15 muon pass-through,
    Muon_isStandalone, Muon_svIdx, the SV matching of jets/taus, the string
    precision fix of SimpleFlatTableProducer, GenVtx_*.
    Safe to call before nanoAOD_customizeCommon (appended by cmsDriver last).
    """
    for taskName, members in (("nanoTableTaskCommon", ("boostedTauTask", "boostedTauTablesTask")),
                              ("nanoTableTaskFS", ("boostedTauMCTask",))):
        task = getattr(process, taskName, None)
        if task is None:
            continue
        for m in members:
            if hasattr(process, m) and task.contains(getattr(process, m)):
                task.remove(getattr(process, m))
    process.linkedObjects.boostedTaus = cms.InputTag("")

    from PhysicsTools.NanoAOD.muonvtxagnosticiso_cff import nanoAOD_addVtxAgnosticIso
    process = nanoAOD_addVtxAgnosticIso(process)

    process = nanoAOD_wmassMuonVariables(process)

    from PhysicsTools.NanoAOD.deepMETPVRobust_cff import nanoAOD_addDeepMETPVRobust
    process = nanoAOD_addDeepMETPVRobust(process)

    # 10_6: no |dxy| < 0.2 && |dz| < 0.1 requirement on the isolated tracks
    # (stock 15_0 keeps it below pt 15)
    process.finalIsolatedTracks.cut = cms.string(
        "((pt>5 && (abs(pdgId) == 11 || abs(pdgId) == 13)) || pt > 10) && (abs(pdgId) < 15 || abs(eta) < 2.5) && "
        "((pfIsolationDR03().chargedHadronIso < 5 && pt < 25) || pfIsolationDR03().chargedHadronIso/pt < 0.2)")
    return process


def nanoAOD_wmassLowPU(process):
    """The 2017 low-PU run (2017H, 13 TeV, ~0 PU) on the UL re-reconstruction
    (RunIILowPUSummer20UL17MiniAODv2 / Run2017H-UL2017_MiniAODv2). The 10_6
    production ran on the 94X low-PU MiniAOD behind a run2_nanoAOD_LowPU era
    modifier; on UL input the standard run2_nanoAOD_106Xv2 path applies and
    only the low-PU-specific content of that era is left, as a customise
    (Category H of the migration):

    * trigger objects of the HI-style menu of the run (nanoAOD_wmassLowPUTriggers)
    * DeepMET re-run with the low-PU models (deepmet_lowPU[_Resp].pb, leptons
      removed from the inputs and added back), replacing the stock
      DeepMETResolutionTune / DeepMETResponseTune tables, which would carry
      the standard-PU models' values stored in the MiniAOD

    Left behind with the 94X era: the VID / scale-smearing source rewiring,
    the puppiIsoId / softMva removal, the Run2017_LowPU_v2 electron
    scale/smearing file and the ecalCorr column (a workaround of the 94X
    MiniAODv2 E/p bug; the UL2017 file of the era applies), the lhcInfoTable
    removal. Unlike 10_6, the production adds the CVH refit
    (nanoAOD_addCvhMuon[MC]) to the low-PU runs as to every other campaign.
    """
    process = nanoAOD_wmassLowPUTriggers(process)

    from RecoMET.METPUSubtraction.deepMETProducer_cfi import deepMETProducer
    process.deepMETsResolutionTuneLowPU = deepMETProducer.clone(
        graph_path = "PhysicsTools/NanoAOD/data/deepmetmodel/deepmet_lowPU.pb",
        ignore_leptons = True,
    )
    process.deepMETsResponseTuneLowPU = deepMETProducer.clone(
        graph_path = "PhysicsTools/NanoAOD/data/deepmetmodel/deepmet_lowPU_Resp.pb",
        ignore_leptons = True,
    )
    for table, producer in (("deepMetResolutionTuneTable", "deepMETsResolutionTuneLowPU"),
                            ("deepMetResponseTuneTable", "deepMETsResponseTuneLowPU")):
        mod = getattr(process, table)
        mod.src = cms.InputTag(producer)
        mod.variables.pt = Var("pt", float, doc=mod.variables.pt.doc.value() + " (low-PU model, leptons excluded from the inputs)", precision=-1)
        mod.variables.phi = Var("phi", float, doc=mod.variables.phi.doc.value() + " (low-PU model, leptons excluded from the inputs)", precision=12)
    process.metTablesTask.add(process.deepMETsResolutionTuneLowPU, process.deepMETsResponseTuneLowPU)
    return process


def nanoAOD_wmassLowPU5TeV(process):
    """The 2017 pp reference run at 5.02 TeV (2017G) on the UL re-reconstruction
    (RunIISummer20UL17pp5TeVMiniAODv2 / Run2017G-UL2017_MiniAODv2): the
    trigger objects of its HI-style menu (the same paths as 2017H), stock
    DeepMET."""
    return nanoAOD_wmassLowPUTriggers(process)


def nanoAOD_beamSpotFromTag(process, tag):
    """Re-make offlineBeamSpot in this process from the BeamSpotObjectsRcd
    payload of `tag` (a GlobalTag override), for MC whose input files carry a
    beamspot that does not describe the simulated collisions.

    Every module that reads 'offlineBeamSpot' without a process name then gets
    this one instead of the input file's: the BeamSpot and PVBS tables, the
    beamspot-constrained muon pt, the lepton/tau impact points (refittedPV),
    the dxybs sign of the lepton updaters, the electron ID conversion veto,
    the low-pt electrons, the secondary vertices of the re-run DeepJet, the
    PV-robust PUPPI fallback point and the CVH refits (reference point of the
    single-track states, beamspot constraint of trackrefitbs and of the dimuon
    common vertex). Quantities computed upstream with the old beamspot and only
    read here are NOT corrected: the primary vertices of the input file and, in
    the nano from MiniAOD, the magnitudes of the pat::Muon / pat::Electron
    dB(BS2D) (dxybs; the tag-and-probe nano runs PAT in the same job and
    recomputes them with the new beamspot)."""
    process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")
    process.GlobalTag.toGet.append(cms.PSet(
        record = cms.string("BeamSpotObjectsRcd"),
        tag = cms.string(tag),
    ))
    process.nanoBeamSpotOverrideTask = cms.Task(process.offlineBeamSpot)
    process.schedule.associate(process.nanoBeamSpotOverrideTask)
    return process


def nanoAOD_beamSpotEarly2018MC(process):
    """MC generated with the 2018 vertex smearing (Realistic25ns13TeVEarly2018Collision:
    x, y = +108, +417 um, 7 um wide) but reconstructed with the 2017 MC beamspot
    (BeamSpotObjects_Realistic25ns_13TeVCollisions_Early2017_v1_mc of
    106X_mc2017_realistic_v9For2017H_v1: -248, +693 um), i.e. the POWHEG-MiNNLO
    W/Z, ttbar and single-top samples of RunIILowPUSummer20UL17 (McM wmLHEGS
    requests; the pomflux samples of the campaign were generated with the 2017
    smearing and are consistent). Their MiniAOD offlineBeamSpot is ~450 um off
    the collisions, which biases every beamspot-constrained quantity (e.g. the
    CVH dimuon mass by -1.9%). The 2018 MC beamspot tag describes the generated
    vertices (checked against genParticles:xyz0)."""
    return nanoAOD_beamSpotFromTag(process, "BeamSpotObjects_Realistic25ns_13TeVCollisions_Early2018_v1_mc")


def nanoAOD_wmassLowPUTriggers(process):
    """Trigger objects of the HI-style menu of the 2017 low-PU runs (2017G at
    5.02 TeV, 2017H at 13 TeV): the Electron and Muon selections carry the
    Ele20 / Ele17HI and Mu17 filters (the 10_6 selections_lowPU) instead of
    the standard-menu bits."""
    from PhysicsTools.NanoAOD.triggerObjects_cff import mksel
    process.triggerObjectTable.selections.Electron = cms.PSet(
        doc = cms.string("PixelMatched e/gamma, low-PU 2017 menu"),  # this may also select photons!
        id = cms.int32(11),
        sel = cms.string("type(92) && pt > 7 && coll('hltEgammaCandidates') && filter('*PixelMatchFilter')"),
        l1seed = cms.string("type(-98)"), l1deltaR = cms.double(0.3),
        skipObjectsNotPassingQualityBits = cms.bool(True),
        qualityBits = cms.VPSet(
            mksel("filter('hltEle20WPLoose1GsfTrackIsoFilter*')", "Ele20"),
            mksel("filter('hltEle17WPLoose1GsfTrackIsoFilterForHI')", "Ele17HI"),
        ),
    )
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


def nanoAOD_wmassMuonVariables(process):
    """The Muon columns the 10_6 custom NanoAOD had on top of the stock table
    (muons_cff.py of WmassNanoProd_10_6_26; standalone* also exist in the MUO
    POG custom_muon_cff with the same definitions)."""
    v = process.muonTable.variables
    v.innerTrackAlgo = Var('? innerTrack().isNonnull() ? innerTrack().algo() : -99', 'int', precision=-1, doc='Track algo enum, check DataFormats/TrackReco/interface/TrackBase.h for details.')
    v.innerTrackOriginalAlgo = Var('? innerTrack().isNonnull() ? innerTrack().originalAlgo() : -99', 'int', precision=-1, doc='Track original algo enum')
    v.trkKink = Var("combinedQuality().trkKink", float, doc="kink finder output")
    v.nTrackerValidHits = Var("?track.isNonnull?innerTrack().hitPattern().numberOfValidHits():0", int, doc="number of valid hits in the tracker")
    v.nPixelValidHits = Var("?track.isNonnull?innerTrack().hitPattern().numberOfValidPixelHits():0", int, doc="number of valid hits in the pixel detector")
    v.pfRelIso04_chg = Var("pfIsolationR04().sumChargedHadronPt/pt", float, doc="PF relative isolation dR=0.4, charged component")
    v.pfRelIso04_neu = Var("pfIsolationR04().sumNeutralHadronEt/pt", float, doc="PF relative isolation dR=0.4, neutral component")
    v.pfRelIso04_pho = Var("pfIsolationR04().sumPhotonEt/pt", float, doc="PF relative isolation dR=0.4, photon component")
    v.pfRelIso04_pu = Var("pfIsolationR04().sumPUPt/pt", float, doc="PF relative isolation dR=0.4, PU component")
    v.standalonePt = Var("? standAloneMuon().isNonnull() ? standAloneMuon().pt() : -1", float, doc="pt of the standalone muon", precision=14)
    v.standaloneEta = Var("? standAloneMuon().isNonnull() ? standAloneMuon().eta() : -99", float, doc="eta of the standalone muon", precision=14)
    v.standalonePhi = Var("? standAloneMuon().isNonnull() ? standAloneMuon().phi() : -99", float, doc="phi of the standalone muon", precision=14)
    v.standaloneCharge = Var("? standAloneMuon().isNonnull() ? standAloneMuon().charge() : -99", float, doc="charge of the standalone muon", precision=14)
    v.standaloneNumberOfValidHits = Var('? standAloneMuon().isNonnull() ? standAloneMuon().numberOfValidHits() : -1', 'int', precision=-1, doc='Number of valid hits in track (standalone chambers)')
    return process

