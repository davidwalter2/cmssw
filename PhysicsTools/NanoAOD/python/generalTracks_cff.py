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
        ## next three for marc's TnP
        trackAlgo = Var('algo()', 'int', precision=-1, doc='Track algo enum, check DataFormats/TrackReco/interface/TrackBase.h for details.'),
        trackOriginalAlgo = Var('originalAlgo()', 'int', precision=-1, doc='Track original algo enum'),
        qualityMask = Var('qualityMask()', 'int', precision=-1, doc='Quality mask of the track.'),
        extraIdx = Var('extra().key()', 'int', precision=-1, doc='Index of the TrackExtra in the original collection'),
    ),
)

standaloneMuonTable = cms.EDProducer("SimpleTrackFlatTableProducer",
    src = cms.InputTag("standAloneMuons"),
    cut = cms.string(""),
    name = cms.string("StandAloneMuon"),
    doc  = cms.string("Extra standalone Muons from track collection"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(P3Vars,
        dz = Var("dz",float,doc="dz (with sign) wrt first PV, in cm",precision=10),
        dxy = Var("dxy",float,doc="dxy (with sign) wrt first PV, in cm",precision=10),
        charge = Var("charge", int, doc="electric charge"),
        normChiSq = Var("normalizedChi2", float, precision=14, doc="Chi^2/ndof"),
        numberOfValidHits = Var('numberOfValidHits()', 'int', precision=-1, doc='Number of valid hits in track'),
        numberOfLostHits = Var('numberOfLostHits()', 'int', precision=-1, doc='Number of lost hits in track'),
        ## next three for marc's TnP
        trackAlgo = Var('algo()', 'int', precision=-1, doc='Track algo enum, check DataFormats/TrackReco/interface/TrackBase.h for details.'),
        trackOriginalAlgo = Var('originalAlgo()', 'int', precision=-1, doc='Track original algo enum'),
        qualityMask = Var('qualityMask()', 'int', precision=-1, doc='Quality mask of the track.'),
        extraIdx = Var('extra().key()', 'int', precision=-1, doc='Index of the TrackExtra in the original collection'),
    ),
)

standaloneMuonUpdatedAtVtxTable = cms.EDProducer("SimpleTrackFlatTableProducer",
    src = cms.InputTag("standAloneMuons:UpdatedAtVtx"),
    cut = cms.string(""),
    name = cms.string("StandAloneMuonUpdatedAtVtx"),
    doc  = cms.string("Extra standalone Muons updated at vertex from track collection"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(P3Vars,
        dz = Var("dz",float,doc="dz (with sign) wrt first PV, in cm",precision=10),
        dxy = Var("dxy",float,doc="dxy (with sign) wrt first PV, in cm",precision=10),
        charge = Var("charge", int, doc="electric charge"),
        normChiSq = Var("normalizedChi2", float, precision=14, doc="Chi^2/ndof"),
        numberOfValidHits = Var('numberOfValidHits()', 'int', precision=-1, doc='Number of valid hits in track'),
        numberOfLostHits = Var('numberOfLostHits()', 'int', precision=-1, doc='Number of lost hits in track'),
        ## next three for marc's TnP
        trackAlgo = Var('algo()', 'int', precision=-1, doc='Track algo enum, check DataFormats/TrackReco/interface/TrackBase.h for details.'),
        trackOriginalAlgo = Var('originalAlgo()', 'int', precision=-1, doc='Track original algo enum'),
        qualityMask = Var('qualityMask()', 'int', precision=-1, doc='Quality mask of the track.'),
        extraIdx = Var('extra().key()', 'int', precision=-1, doc='Index of the TrackExtra in the original collection'),
    ),
)

mergedStandaloneMuonTable = cms.EDProducer("SimpleTrackFlatTableProducer",
    src = cms.InputTag("mergedStandAloneMuons"),
    cut = cms.string(""),
    name = cms.string("MergedStandAloneMuon"),
    doc  = cms.string("Disambiguated collection of standalone Muons following logic of GlobalMuonProducer"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(P3Vars,
        dz = Var("dz",float,doc="dz (with sign) wrt first PV, in cm",precision=10),
        dxy = Var("dxy",float,doc="dxy (with sign) wrt first PV, in cm",precision=10),
        charge = Var("charge", int, doc="electric charge"),
        normChiSq = Var("normalizedChi2", float, precision=14, doc="Chi^2/ndof"),
        numberOfValidHits = Var('numberOfValidHits()', 'int', precision=-1, doc='Number of valid hits in track'),
        numberOfLostHits = Var('numberOfLostHits()', 'int', precision=-1, doc='Number of lost hits in track'),
        ## next three for marc's TnP
        trackAlgo = Var('algo()', 'int', precision=-1, doc='Track algo enum, check DataFormats/TrackReco/interface/TrackBase.h for details.'),
        trackOriginalAlgo = Var('originalAlgo()', 'int', precision=-1, doc='Track original algo enum'),
        qualityMask = Var('qualityMask()', 'int', precision=-1, doc='Quality mask of the track.'),
        extraIdx = Var('extra().key()', 'int', precision=-1, doc='Index of the TrackExtra in the original collection'),
    ),
    externalVariables = cms.PSet(
        isUpdatedAtVtx = ExtVar(cms.InputTag("mergedStandAloneMuons:updatedAtVtx"),bool, doc="is standalone muon track updated at vertex"),
    )
)
