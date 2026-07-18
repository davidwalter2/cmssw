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
# Override with the CVH_SCALARPOT_INITFILE env var (cmsDriver --customise
# cannot pass function arguments).
_DEFAULT_SCALARPOT_INITFILE = (
    "/work/submit/david_w/ZMass/mfs/data/fitresults/"
    "polyfit3d_full_coeffs_lmax18_cmsswnorm.txt")


def setup3DFieldForRefit(process, initFile=None, useScalarPot3D=True):
    """Set up the Geant4e propagator + shared G4 master the CVH muon refit needs,
    and choose the baseline magnetic field.

    useScalarPot3D=True (DATA): replace the baseline field with the accurate 3D
    scalar-potential map (ScalarPot3D), routed into every CVH-side consumer (the
    propagator, the strip/pixel CPEs, the refit makers, and the shared G4
    master). useScalarPot3D=False (MC): keep the DEFAULT CMSSW field -- it must
    stay consistent with the field the simulation used (see CLAUDE.md, "custom
    field map for data only"); only the geometry / propagator / master are set
    up, with no field override.

    The parmtype-14 correction-mode basis (scalarPotentialInitFile on the
    makers) is set in BOTH cases -- it defines the Jacobian modes stored in
    NanoAOD and is independent of the baseline field.

    geopro is loaded (for the propagator) but is NOT scheduled -- the shared
    CvhMasterThread ES product owns the G4 world.

    initFile: scalar-potential coefficient dump (basis + coefficients). Falls
    back to CVH_SCALARPOT_INITFILE then _DEFAULT_SCALARPOT_INITFILE.
    """
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
    process.load("Configuration.StandardSequences.GeometrySimDB_cff")
    process.GlobalTag.toGet.append(
        cms.PSet(
            record=cms.string("GeometryFileRcd"),
            tag=cms.string("XMLFILE_Geometry_2016_81YV1_Extended2016_mc"),
            label=cms.untracked.string("Extended"),
        )
    )
    if hasattr(process, "XMLFromDBSource"):
        process.XMLFromDBSource.label = cms.string("Extended")

    process.load("TrackPropagation.Geant4e.geantRefit_cff")

    _refits = ("trackrefit", "trackrefitideal", "trackrefitbs")
    # The parmtype-14 correction-mode basis is needed by the makers in both
    # data and MC (it defines the Jacobians stored in NanoAOD) and is
    # independent of the baseline field.
    for _refit in _refits:
        if hasattr(process, _refit):
            getattr(process, _refit).scalarPotentialInitFile = cms.string(initFile)

    # Baseline field. DATA: the accurate ScalarPot3D map, routed into every
    # CVH-side consumer. MC: keep the default CMSSW field (label "") so the
    # refit stays consistent with the field the simulation used -- no override.
    fieldlabel = ""
    if useScalarPot3D:
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
    # via esConsumes. Uses the same baseline field as the refit -- ScalarPot3D
    # for data, the default field (label "") for MC.
    from TrackPropagation.Geant4e.cvhMasterESProducer_cfi import cvhMasterESProducer
    process.cvhMasterESProducer = cvhMasterESProducer.clone()
    process.cvhMasterESProducer.MagneticFieldLabel = cms.string(fieldlabel)

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
                          ("trackrefitbs", 323456789)):
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

