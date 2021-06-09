import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

generalTrackTable = cms.EDProducer("SimpleTrackFlatTableProducer",
    src = cms.InputTag("generalTracks"),
    cut = cms.string("pt > 15"), # filtered already above
    name = cms.string("Track"),
    doc  = cms.string("General tracks with pt > 15 GeV"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(P3Vars,
        dz = Var("dz",float,doc="dz (with sign) wrt first PV, in cm",precision=10),
        dxy = Var("dxy",float,doc="dxy (with sign) wrt first PV, in cm",precision=10),
        charge = Var("charge", int, doc="electric charge"),
        normChiSq = Var("normalizedChi2", float, precision=14, doc="Chi^2/ndof"),
        numberOfValidHits = Var('numberOfValidHits()', 'int', precision=-1, doc='Number of valid hits in track'),
        numberOfLostHits = Var('numberOfLostHits()', 'int', precision=-1, doc='Number of lost hits in track'),
    ),
)


