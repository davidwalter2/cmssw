from __future__ import print_function
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
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
from PhysicsTools.NanoAOD.NanoAODEDMEventContent_cff import *

nanoMetadata = cms.EDProducer("UniqueStringProducer",
    strings = cms.PSet(
        tag = cms.string("untagged"),
    )
)

finalIsolatedTracks.finalLeptons = ["finalLooseMuons"]

finalMuons.src = "slimmedMuonsUpdated"
finalLooseMuons.src = "slimmedMuonsUpdated"
muonTable.src = "finalMuons"
muonTable.externalVariables = cms.PSet()


muonSimpleSequence= cms.Sequence(slimmedMuonsUpdated+isoForMu + finalMuons + finalLooseMuons )

nanotpSequence = cms.Sequence(
        nanoMetadata + 
        genParticleSequence + genVertexTables + 
        muonSequence + vertexSequence+
        isoTrackSequence + # must be after all the leptons
        muonMC + muonTable + vertexTables+ isoTrackTables
        )
