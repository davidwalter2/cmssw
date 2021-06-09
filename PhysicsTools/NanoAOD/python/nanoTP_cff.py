from __future__ import print_function
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.globals_cff import *
from PhysicsTools.NanoAOD.muons_cff import *
from PhysicsTools.NanoAOD.electrons_cff import *
from PhysicsTools.NanoAOD.genparticles_cff import *
from PhysicsTools.NanoAOD.particlelevel_cff import *
from PhysicsTools.NanoAOD.genWeightsTable_cfi import *
from PhysicsTools.NanoAOD.genVertex_cff import *
from PhysicsTools.NanoAOD.vertices_cff import *
from PhysicsTools.NanoAOD.met_cff import *
from PhysicsTools.NanoAOD.triggerObjects_cff import *
from PhysicsTools.NanoAOD.isotracks_cff import *
from PhysicsTools.NanoAOD.generalTracks_cff import *
from PhysicsTools.NanoAOD.NanoAODEDMEventContent_cff import *
from PhysicsTools.PatAlgos.slimming.miniAOD_tools import *

nanoMetadata = cms.EDProducer("UniqueStringProducer",
    strings = cms.PSet(
        tag = cms.string("untagged"),
    )
)

l1bits=cms.EDProducer("L1TriggerResultsConverter", src=cms.InputTag("gtStage2Digis"), legacyL1=cms.bool(False),
                      storeUnprefireableBit=cms.bool(True), src_ext=cms.InputTag("gtStage2Digis"))

finalIsolatedTracks.finalLeptons = ["finalLooseMuons"]

finalMuons.src = "slimmedMuonsUpdated"
finalLooseMuons.src = "slimmedMuonsUpdated"
muonTable.src = "finalMuons"
muonTable.externalVariables = cms.PSet()
muonTable.variables = cms.PSet(muonTable.variables,
        standalonePt = Var("? standAloneMuon().isNonnull() ? standAloneMuon().pt() : -1", float, doc = "pt of the standalone muon", precision=14),
        standaloneEta = Var("? standAloneMuon().isNonnull() ? standAloneMuon().eta() : -99", float, doc = "eta of the standalone muon", precision=14),
        standalonePhi = Var("? standAloneMuon().isNonnull() ? standAloneMuon().phi() : -99", float, doc = "phi of the standalone muon", precision=14),
        standaloneCharge = Var("? standAloneMuon().isNonnull() ? standAloneMuon().charge() : -99", float, doc = "phi of the standalone muon", precision=14),
)

muonSimpleSequence= cms.Sequence(slimmedMuonsUpdated+isoForMu + finalMuons + finalLooseMuons )

nanotpSequence = cms.Sequence(
        nanoMetadata + triggerObjectTables + l1bits +
        genParticleSequence + genParticleTables + genWeightsTable + genVertexTables + 
        muonSequence + vertexSequence+
        isoTrackSequence + # must be after all the leptons
		puTable + genTable +
        muonMC + muonTable + vertexTables+ isoTrackTables + generalTrackTable
        )

def customizeMuonPassThrough(process):
	passStandalone = "(standAloneMuon().isNonnull() && standAloneMuon().pt() > 15)"
	process.selectedPatMuons.cut = cms.string("||".join([passStandalone, process.selectedPatMuons.cut.value()]))
	process.finalMuons.cut = cms.string("||".join([passStandalone, process.finalMuons.cut.value()]))
	return process

