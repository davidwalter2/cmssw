# Flat-pT J/psi -> mumu gun for the rung-E closure sample: produced INSIDE
# CMSSW_15_0 so the simulation Geant4 matches the CVH refit propagator.
# Kinematics roughly cover the TkAlJpsiMuMu phase space (muon pT 3-20).
import FWCore.ParameterSet.Config as cms

generator = cms.EDFilter(
    "Pythia8PtGun",
    PGunParameters=cms.PSet(
        ParticleID=cms.vint32(443),
        AddAntiParticle=cms.bool(False),
        MinPt=cms.double(5.0),
        MaxPt=cms.double(30.0),
        MinEta=cms.double(-2.4),
        MaxEta=cms.double(2.4),
        MinPhi=cms.double(-3.14159265359),
        MaxPhi=cms.double(3.14159265359),
    ),
    PythiaParameters=cms.PSet(
        pythia8CommonSettingsBlock=cms.PSet(
            parameterSets=cms.vstring()
        ),
        jpsiDecay=cms.vstring(
            "443:onMode = off",
            "443:onIfMatch = 13 -13",
        ),
        parameterSets=cms.vstring("jpsiDecay"),
    ),
    Verbosity=cms.untracked.int32(0),
    firstRun=cms.untracked.uint32(1),
    psethack=cms.string("Jpsi mumu flat pT gun"),
)
