import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, ExtVar
from PhysicsTools.NanoAOD.met_cff import simpleSingletonPATMETFlatTableProducer
from PhysicsTools.NanoAOD.globals_cff import globalVariablesTableProducer
from RecoMET.METPUSubtraction.deepMETPVRobust_cff import puppiPVRobust, deepMETsPVRobust, deepMETsPVRobustNoPUPPI

# WMass: the PV-robust DeepMET of the 10_6 custom NanoAOD (Category F of the
# migration). 10_6 routed the estimate through pat::MET corrections
# (RawDeepPVRobust) and the MET slimmer; the 15_0 nano does not re-run the
# slimmer, so the tables read the producers' pat::MET directly -- same numbers,
# same branch names.
deepMetPVRobustTable = simpleSingletonPATMETFlatTableProducer.clone(
    src = cms.InputTag("deepMETsPVRobust"),
    name = cms.string("DeepMETPVRobust"),
    doc = cms.string("DeepMET with PV-robust inputs (dz and PUPPI weight w.r.t. the vertex closest to the leading muon), for mW"),
    variables = cms.PSet(
        pt = Var("pt", float, doc="DeepMET PVRobust pt", precision=-1),
        phi = Var("phi", float, doc="DeepMET PVRobust phi", precision=12),
    ),
)
deepMetPVRobustNoPUPPITable = deepMetPVRobustTable.clone(
    src = "deepMETsPVRobustNoPUPPI",
    name = "DeepMETPVRobustNoPUPPI",
    doc = "DeepMET with PV-robust inputs, without PUPPI, for mW",
    variables = cms.PSet(
        pt = Var("pt", float, doc="DeepMET PVRobustNoPUPPI pt", precision=-1),
        phi = Var("phi", float, doc="DeepMET PVRobustNoPUPPI phi", precision=12),
    ),
)
pvRobustTable = globalVariablesTableProducer.clone(
    name = cms.string(""),
    variables = cms.PSet(
        PVRobustIndex = ExtVar(cms.InputTag("puppiPVRobust:PVRobustIndex"), "int", doc="index of the PV closest to the leading loose muon; -1 means failure to find any vertex close to the leading muon within 0.2 cm, use the beamSpot and muon vz in this case"),
        PVMuonIndex = ExtVar(cms.InputTag("puppiPVRobust:PVMuonIndex"), "int", doc="index of the slimmedMuons collection that is used to determine the robust PV index; -1 if no muon with pt > 10 passed the loose ID"),
    ),
)

deepMETPVRobustTablesTask = cms.Task(puppiPVRobust, deepMETsPVRobust, deepMETsPVRobustNoPUPPI,
                                     deepMetPVRobustTable, deepMetPVRobustNoPUPPITable, pvRobustTable)


def nanoAOD_addDeepMETPVRobust(process):
    """Schedule the PV-robust PUPPI + DeepMET producers and their tables
    (DeepMETPVRobust_pt/phi, DeepMETPVRobustNoPUPPI_pt/phi, PVRobustIndex,
    PVMuonIndex) with the MET tables; data and MC."""
    process.puppiPVRobust = puppiPVRobust
    process.deepMETsPVRobust = deepMETsPVRobust
    process.deepMETsPVRobustNoPUPPI = deepMETsPVRobustNoPUPPI
    process.deepMetPVRobustTable = deepMetPVRobustTable
    process.deepMetPVRobustNoPUPPITable = deepMetPVRobustNoPUPPITable
    process.pvRobustTable = pvRobustTable
    process.metTablesTask.add(process.puppiPVRobust, process.deepMETsPVRobust, process.deepMETsPVRobustNoPUPPI,
                              process.deepMetPVRobustTable, process.deepMetPVRobustNoPUPPITable, process.pvRobustTable)
    return process
