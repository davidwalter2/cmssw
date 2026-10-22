# AlCaReco NanoAOD: the CVH single-track refit of the custom NanoAOD and its
# Muon columns, for the muon tracks of a tracker-alignment AlCaReco
# (TkAlJpsiMuMu, TkAlUpsilonMuMu).
#
# The AlCaReco keeps the selected muon tracks with their hits, the primary
# vertices, the trigger results and, in MC, the generator record; it has no
# reco::Muon, PF candidates, jets or trigger objects. Hence:
#   * the Muon table has one row per AlCaReco track (all muon tracks by the
#     skim's selection) and the columns of the MiniAOD NanoAOD's Muon table
#     that the track determines: kinematics (muon mass and pdgId), ptErr,
#     tightCharge, impact parameters (pat::Muon definitions), the beam-spot
#     constrained pT, hit counts, quality, algorithm, the track covariance
#     elements, and all CVH columns (cvh*, in MC also cvhideal*, and the
#     Jacobian payload), with the same names, definitions and precisions;
#     muonBestTrack() and tunePMuonBestTrack() are the AlCaReco track;
#   * PV/OtherPV (without the PF-candidate sums), BeamSpot (from the
#     conditions, the AlCaReco keeps none), the HLT bits, the LHC info (data)
#     and in MC GenPart (without the PF-based iso), Muon_genPart*, GenVtx,
#     Pileup, Generator, genWeight and GenFilter, as in the MiniAOD NanoAOD;
#   * the muon ID, isolation and standalone-track columns, and the jets, MET,
#     taus, electrons, photons, isolated tracks, secondary vertices, trigger
#     objects, L1 bits and MET filters are not in the AlCaReco and not written.
#
#
# Two versions, chosen by the flag `calib` of nanoAOD_alcareco:
#   * calib=False: the content above, from which the module-level corrections
#     are APPLIED (Muon_cvhJacRef and the indices), as in the MiniAOD NanoAOD;
#   * calib=True (J/psi only): additionally the content from which they are
#     DERIVED: the Dimuon table of the two-track CVH refit of the opposite-sign
#     track pairs, with the common-vertex and the J/psi mass constraint, and its
#     global-fit payload (the gradient and the factored Hessian of the fit's
#     objective wrt the global parameters, their indices, d(mass)/d(parms) and
#     the per-leg Jacobians, what runCvhJpsi.py writes for the standalone
#     calibration), plus the parameter catalog in the Runs tree (CvhParm_*).
#     Both versions have the same correction model, so the corrections derived
#     from one apply to the Jacobians of the other.
#
# cmsDriver loads this cff in place of nano_cff, so its default NANO
# customisation is nanoAOD_customizeCommon below (the CVH field setup):
#   cmsDriver.py NANO -s NANO:PhysicsTools/NanoAOD/alcareco_nano_cff.alcarecoNanoSequence[MC] \
#       --customise PhysicsTools/NanoAOD/alcareco_nano_cff.nanoAOD_alcareco<skim>[Calib] ...
# with <skim> TkAlJpsiMuMu (the default track source) or TkAlUpsilonMuMu, and
# the Calib suffix (TkAlJpsiMuMu only) for calib=True.
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, ExtVar
from PhysicsTools.NanoAOD.nano_cff import nanoMetadata, lhcInfoTable
from PhysicsTools.NanoAOD.globals_cff import beamSpotTable, genTable, genFilterTable
from PhysicsTools.NanoAOD.globals_cff import puTable as _puTable
from PhysicsTools.NanoAOD.vertices_cff import vertexTable as _vertexTable
from PhysicsTools.NanoAOD.genWeightsTable_cfi import genWeightsTable
from PhysicsTools.NanoAOD.genVertex_cff import genVertexTable, genVertexT0Table
from PhysicsTools.NanoAOD.genparticles_cff import finalGenParticles as _finalGenParticles
from PhysicsTools.NanoAOD.genparticles_cff import genParticleTable as _genParticleTable
from PhysicsTools.NanoAOD import muons_cff as _muons
from RecoVertex.BeamSpotProducer.BeamSpot_cfi import offlineBeamSpot
from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import TransientTrackBuilderESProducer

_DEFAULT_SKIM = "TkAlJpsiMuMu"
_tracks = cms.InputTag("ALCARECO" + _DEFAULT_SKIM)

# --- CVH refits, keyed to the AlCaReco tracks (no muon association) ----------
trackrefit = _muons.trackrefit.clone(src = _tracks, doMuonAssoc = False, trackKeyedOutputs = cms.bool(True))
trackrefitideal = _muons.trackrefitideal.clone(src = _tracks, doMuonAssoc = False, trackKeyedOutputs = cms.bool(True))
mergedGlobalIdxs = _muons.mergedGlobalIdxs.clone()

# --- Muon table on the AlCaReco tracks ----------------------------------------
muonIP = cms.EDProducer("TrackPatImpactParameterProducer",
    src = _tracks,
    pvSrc = cms.InputTag("offlinePrimaryVertices"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
)
muonBSConstrain = cms.EDProducer("TrackBeamspotConstraintValueMapProducer",
    src = _tracks,
    beamspot = cms.InputTag("offlineBeamSpot"),
    vertices = cms.InputTag("offlinePrimaryVertices"),
)

def _muonVar(name, expr):
    """The MiniAOD NanoAOD's Muon column `name`, evaluated as `expr` on the track."""
    v = getattr(_muons.muonTable.variables, name).clone()
    v.expr = expr
    return v

def _muonExtVar(name, src):
    """The MiniAOD NanoAOD's Muon external column `name`, read from `src`."""
    v = getattr(_muons.muonTable.externalVariables, name).clone()
    v.src = src
    return v

_v = _muons.muonTable.variables
_e = _muons.muonTable.externalVariables
muonTable = cms.EDProducer("SimpleTrackFlatTableProducer",
    src = _tracks,
    cut = cms.string(""),
    name = cms.string("Muon"),
    doc = cms.string("AlCaReco muon tracks (ALCARECO%s)" % _DEFAULT_SKIM),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(
        pt = _v.pt.clone(),
        eta = _v.eta.clone(),
        phi = _v.phi.clone(),
        mass = _muonVar("mass", "0.1056583755"),
        pdgId = _muonVar("pdgId", "-13*charge"),
        charge = _v.charge.clone(),
        ptErr = _muonVar("ptErr", "ptError()"),
        tightCharge = _muonVar("tightCharge", "?(ptError()/pt() < 0.2)?2:0"),
        nTrackerLayers = _muonVar("nTrackerLayers", "hitPattern().trackerLayersWithMeasurement()"),
        highPurity = _muonVar("highPurity", "quality('highPurity')"),
        innerTrackAlgo = Var("algo()", "int", precision=-1, doc="Track algo enum, check DataFormats/TrackReco/interface/TrackBase.h for details."),
        innerTrackOriginalAlgo = Var("originalAlgo()", "int", precision=-1, doc="Track original algo enum"),
        nTrackerValidHits = Var("hitPattern().numberOfValidHits()", int, doc="number of valid hits in the tracker"),
        nPixelValidHits = Var("hitPattern().numberOfValidPixelHits()", int, doc="number of valid hits in the pixel detector"),
        VXBS_Cov00 = _muonVar("VXBS_Cov00", "covariance(0,0)"),
        VXBS_Cov03 = _muonVar("VXBS_Cov03", "covariance(0,3)"),
        VXBS_Cov33 = _muonVar("VXBS_Cov33", "covariance(3,3)"),
    ),
    externalVariables = cms.PSet(
        dxy = ExtVar(cms.InputTag("muonIP:dxy"), float, doc=_v.dxy.doc.value(), precision=_v.dxy.precision.value()),
        dxyErr = ExtVar(cms.InputTag("muonIP:dxyErr"), float, doc=_v.dxyErr.doc.value(), precision=_v.dxyErr.precision.value()),
        dz = ExtVar(cms.InputTag("muonIP:dz"), float, doc=_v.dz.doc.value(), precision=_v.dz.precision.value()),
        dzErr = ExtVar(cms.InputTag("muonIP:dzErr"), float, doc=_v.dzErr.doc.value(), precision=_v.dzErr.precision.value()),
        dxybs = ExtVar(cms.InputTag("muonIP:dxybs"), float, doc=_v.dxybs.doc.value(), precision=_v.dxybs.precision.value()),
        dxybsErr = ExtVar(cms.InputTag("muonIP:dxybsErr"), float, doc=_v.dxybsErr.doc.value(), precision=_v.dxybsErr.precision.value()),
        ip3d = ExtVar(cms.InputTag("muonIP:ip3d"), float, doc=_v.ip3d.doc.value(), precision=_v.ip3d.precision.value()),
        sip3d = ExtVar(cms.InputTag("muonIP:sip3d"), float, doc=_v.sip3d.doc.value(), precision=_v.sip3d.precision.value()),
        bsConstrainedPt = _e.bsConstrainedPt.clone(),
        bsConstrainedPtErr = _e.bsConstrainedPtErr.clone(),
        bsConstrainedChi2 = _e.bsConstrainedChi2.clone(),
    ),
)
_muons._cvhScalarVars(muonTable.externalVariables, "trackrefit", "", "")

# MC: the ideal-geometry refit, an extension of the Muon table
muonCvhIdealTable = cms.EDProducer("SimpleTrackFlatTableProducer",
    src = _tracks,
    cut = cms.string(""),
    name = cms.string("Muon"),
    doc = cms.string("CVH refit with the ideal geometry"),
    singleton = cms.bool(False),
    extension = cms.bool(True),
    variables = cms.PSet(),
    externalVariables = cms.PSet(),
)
_muons._cvhScalarVars(muonCvhIdealTable.externalVariables, "trackrefitideal", "ideal", " (ideal geometry)")

# The Jacobian payload, as muonExternalVecVarsTable of the MiniAOD NanoAOD:
# data: the nominal refit; MC: the indices merged with the ideal refit's, and
# the ideal refit's Jacobian and covariance
_vecVars = _muons.muonExternalVecVarsTable.variables
muonExternalVecVarsTable = cms.EDProducer("FlattenedTrackValueMapVectorTableProducer",
    name = cms.string("Muon"),
    src = _tracks,
    cut = cms.string(""),
    doc = muonTable.doc,
    variables = cms.PSet(
        cvhmergedGlobalIdxs = _vecVars.cvhmergedGlobalIdxs.clone(),
        cvhJacRef = _vecVars.cvhJacRef.clone(),
        cvhMomCov = _vecVars.cvhMomCov.clone(),
    ),
)
muonExternalVecVarsMCTable = muonExternalVecVarsTable.clone()
muonExternalVecVarsMCTable.variables.cvhmergedGlobalIdxs = ExtVar(cms.InputTag("mergedGlobalIdxs"), "std::vector<int>", doc="Indices for correction parameters (merged nominal+ideal)", precision=16)
muonExternalVecVarsMCTable.variables.cvhidealJacRef = ExtVar(cms.InputTag("trackrefitideal:jacRef"), "std::vector<float>", doc="jacobian for corrections (ideal geometry)", precision=12)
muonExternalVecVarsMCTable.variables.cvhidealMomCov = ExtVar(cms.InputTag("trackrefitideal:momCov"), "std::vector<float>", doc="covariance matrix for qop, lambda, phi (ideal geometry)", precision=12)

# MC: Muon_genPartIdx / genPartFlav, matching candidates made from the tracks
# (same order, one per track) to finalGenParticles
muonCands = cms.EDProducer("ConcreteChargedCandidateProducer",
    src = _tracks,
    particleType = cms.string("mu+"),
)
muonsMCMatchForTable = _muons.muonsMCMatchForTable.clone(src = cms.InputTag("muonCands"))
muonMCTable = _muons.muonMCTable.clone(src = cms.InputTag("muonCands"), objName = "Muon", objType = "Muon")

# --- event tables ----------------------------------------------------------------
vertexTable = cms.EDProducer("RecoVertexTableProducer",
    pvSrc = cms.InputTag("offlinePrimaryVertices"),
    goodPvCut = _vertexTable.goodPvCut,
    pvName = _vertexTable.pvName,
)
puTable = _puTable.clone(src = "addPileupInfo", pvsrc = "offlinePrimaryVertices")
# the AlCaReco MC keeps the full generator record: prune it as the MiniAOD's
finalGenParticles = _finalGenParticles.clone(src = "genParticles")
genParticleTable = _genParticleTable.clone(externalVariables = cms.PSet())  # iso: from the packed gen particles

# --- calibration content (calib=True, J/psi) -----------------------------------
# Underscored, so that loading the cff does not put them into every process.
# The opposite-sign pairs of the AlCaReco tracks in a mass window wider than the
# skim's (2.7-3.4 GeV, tracks pt > 0.8 GeV, |eta| < 3.5); no muon pt or |eta|
# cut is added.
from Analysis.HitAnalyzer.diMuonTrackVertexCandidates_cfi import diMuonTrackVertexCandidates as _diMuonCands
_calibDimuonCandidates = _diMuonCands.clone(src = _tracks, massMin = 2.6, massMax = 3.6)
# The two-track refit as for the standalone J/psi calibration (runCvhJpsi.py):
# common vertex, J/psi mass constraint, no luminous-region constraint (the
# B -> J/psi X fraction is displaced), the factored-Hessian payload. Its
# correction model and field are those of trackrefit (setup3DFieldForRefit).
_calibTrackrefitDimuon = _muons.trackrefitdimuon.clone(
    src = _tracks,
    srcCandidates = cms.InputTag("dimuonCandidates"),
    doVtxConstraint = True,
    bsConstraint = False,
    doMassConstraint = True,
    massConstraint = 3.0969,
    massConstraintWidth = 1e-5,
    fillGradsFactored = cms.untracked.bool(True),
    produceValueMaps = True,
    produceCatalog = cms.bool(True),
)
_calibDimuonTable = _muons.dimuonTable.clone(
    src = cms.InputTag("dimuonCandidates"),
    doc = cms.string("Opposite-sign AlCaReco muon-track pairs with the CVH two-track refit "
                     "(common vertex, J/psi mass constraint); the refit's mass and kinematics are "
                     "those without the mass constraint"),
)
_ev = _calibDimuonTable.externalVariables
_ev.cvhEdmvalRef = ExtVar(cms.InputTag("trackrefitdimuon:edmvalRef"), float, doc="CVH refit estimated distance to minimum in the reference parametrization", precision=-1)
_ev.cvhChi2 = ExtVar(cms.InputTag("trackrefitdimuon:chi2"), float, doc="CVH refit chi2", precision=-1)
_ev.cvhNdof = ExtVar(cms.InputTag("trackrefitdimuon:ndof"), int, doc="CVH refit number of degrees of freedom")
_ev.cvhMuPlusConsPt = ExtVar(cms.InputTag("trackrefitdimuon:muPlusConsPt"), float, doc="CVH refit mu+ pt with the mass constraint", precision=-1)
_ev.cvhMuMinusConsPt = ExtVar(cms.InputTag("trackrefitdimuon:muMinusConsPt"), float, doc="CVH refit mu- pt with the mass constraint", precision=-1)
_ev.cvhMuPlusNValid = ExtVar(cms.InputTag("trackrefitdimuon:muPlusNValid"), int, doc="valid hits of the mu+ leg in the refit")
_ev.cvhMuPlusNValidPixel = ExtVar(cms.InputTag("trackrefitdimuon:muPlusNValidPixel"), int, doc="valid pixel hits of the mu+ leg in the refit")
_ev.cvhMuMinusNValid = ExtVar(cms.InputTag("trackrefitdimuon:muMinusNValid"), int, doc="valid hits of the mu- leg in the refit")
_ev.cvhMuMinusNValidPixel = ExtVar(cms.InputTag("trackrefitdimuon:muMinusNValidPixel"), int, doc="valid pixel hits of the mu- leg in the refit")
_ev.cvhGradMax = ExtVar(cms.InputTag("trackrefitdimuon:gradMax"), float, doc="largest |gradient| entry", precision=-1)
_calibDimuonVecVarsTable = _muons.dimuonVecVarsTable.clone(src = cms.InputTag("dimuonCandidates"))
_calibDimuonVecVarsTable.variables.cvhGrad = ExtVar(cms.InputTag("trackrefitdimuon:grad"), "std::vector<float>", doc="d(objective)/d(globalparms) at cvhGlobalIdxs (the gradient of the mass-constrained fit)", precision=-1)
_calibParmCatalogTable = cms.EDProducer("CvhParmCatalogRunTableProducer",
    src = cms.InputTag("trackrefitdimuon", "catalog"),
    name = cms.string("CvhParm"),
)

# --- sequences ---------------------------------------------------------------------
alcarecoNanoTask = cms.Task(
    nanoMetadata, offlineBeamSpot, trackrefit, muonIP, muonBSConstrain,
    muonTable, muonExternalVecVarsTable, vertexTable, beamSpotTable,
)
alcarecoNanoSequence = cms.Sequence(lhcInfoTable, alcarecoNanoTask)

alcarecoNanoTaskMC = cms.Task(
    nanoMetadata, offlineBeamSpot, trackrefit, trackrefitideal, mergedGlobalIdxs,
    muonIP, muonBSConstrain, muonTable, muonCvhIdealTable, muonExternalVecVarsMCTable,
    muonCands, finalGenParticles, muonsMCMatchForTable, muonMCTable,
    vertexTable, beamSpotTable, puTable, genParticleTable, genTable, genWeightsTable,
    genFilterTable, genVertexTable, genVertexT0Table,
)
alcarecoNanoSequenceMC = cms.Sequence(alcarecoNanoTaskMC)

_trackConsumers = ("trackrefit", "trackrefitideal", "muonIP", "muonBSConstrain", "muonTable",
                   "muonCvhIdealTable", "muonExternalVecVarsTable", "muonExternalVecVarsMCTable", "muonCands")


def nanoAOD_alcareco(process, skim, calib=False):
    """Run on the tracks of the AlCaReco skim `skim` (ALCARECO<skim>).

    calib=True (TkAlJpsiMuMu only): add the content from which the module-level
    corrections are derived, the two-track J/psi refit with its Dimuon table,
    global-fit payload and parameter catalog (see the top of this file)."""
    for label in _trackConsumers:
        if hasattr(process, label):
            getattr(process, label).src = cms.InputTag("ALCARECO" + skim)
    for label in ("muonTable", "muonExternalVecVarsTable", "muonExternalVecVarsMCTable"):
        if hasattr(process, label):
            getattr(process, label).doc = "AlCaReco muon tracks (ALCARECO%s)" % skim
    if calib:
        if skim != "TkAlJpsiMuMu":
            raise RuntimeError("nanoAOD_alcareco: calib=True is for TkAlJpsiMuMu only (the J/psi mass constraint)")
        tracks = cms.InputTag("ALCARECO" + skim)
        process.dimuonCandidates = _calibDimuonCandidates.clone(src = tracks)
        process.trackrefitdimuon = _calibTrackrefitDimuon.clone(src = tracks)
        process.dimuonTable = _calibDimuonTable.clone()
        process.dimuonVecVarsTable = _calibDimuonVecVarsTable.clone()
        process.cvhParmCatalogTable = _calibParmCatalogTable.clone()
        isMC = process.nanoAOD_step.contains(process.trackrefitideal)
        task = process.alcarecoNanoTaskMC if isMC else process.alcarecoNanoTask
        task.add(process.dimuonCandidates, process.trackrefitdimuon, process.dimuonTable,
                 process.dimuonVecVarsTable)
        # on the path, not in the task: it has no event product, so unscheduled
        # it would never see the event that carries the catalog
        process.nanoAOD_step += process.cvhParmCatalogTable
        # cmsDriver runs this before nanoAOD_customizeCommon, whose
        # setup3DFieldForRefit gives trackrefitdimuon the correction model, field
        # and seed of the CVH refits; if that has run already, copy them now.
        if hasattr(process, "cvhMasterESProducer"):
            from PhysicsTools.NanoAOD.nano_cff import _cvhCopyRefitSetup
            _cvhCopyRefitSetup(process, "trackrefitdimuon", 423456789)
    return process


def nanoAOD_alcarecoTkAlJpsiMuMu(process):
    return nanoAOD_alcareco(process, "TkAlJpsiMuMu")


def nanoAOD_alcarecoTkAlJpsiMuMuCalib(process):
    return nanoAOD_alcareco(process, "TkAlJpsiMuMu", calib=True)


def nanoAOD_alcarecoTkAlUpsilonMuMu(process):
    return nanoAOD_alcareco(process, "TkAlUpsilonMuMu")


def nanoAOD_customizeCommon(process):
    """cmsDriver's default NANO customisation, which is this function when the
    NANO step loads this cff: the field, geometry and Geant4e setup of the CVH
    refits (nano_cff.setup3DFieldForRefit), with the OPERA field map for data
    and the default field for MC, the field of the simulation; in MC also the
    full-precision gen leptons of the MiniAOD NanoAOD (nanoGenWmassCustomize).
    MC is the sequence alcarecoNanoSequenceMC, which runs the ideal-geometry
    refit."""
    from PhysicsTools.NanoAOD.nano_cff import setup3DFieldForRefit, customizeGenLeptonPrecision
    isMC = process.nanoAOD_step.contains(process.trackrefitideal)
    if isMC:
        process = customizeGenLeptonPrecision(process)
    return setup3DFieldForRefit(process, useScalarPot3D=not isMC)
