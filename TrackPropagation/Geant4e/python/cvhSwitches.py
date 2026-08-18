"""cmsRun options for the CVH energy-loss switches.

These were `CVH_*` ENVIRONMENT VARIABLES until 2026-08-18.  They are
ParameterSet parameters on `Geant4ePropagator` now (read by
`cvhcgf::configure` in GeantPropagatorESProducer's constructor), which is the
CMSSW convention and which puts them in the output file's provenance -- an
environment variable left no trace in the data, so a model file could not be
told apart from one exported with different physics.

This module is the ONE place that maps a command-line option onto the PSet
parameter, so a driver picks them all up with two calls:

    import TrackPropagation.Geant4e.cvhSwitches as cvhSwitches
    cvhSwitches.register(opts)          # before opts.parseArguments()
    ...
    cvhSwitches.apply(process, opts)    # after the propagator is in `process`

Every option defaults to -1 meaning "not specified": the driver then leaves the
cfi's own default alone, so passing nothing reproduces the default
configuration exactly and an arm only overrides what it names.
"""

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

# option name == PSet parameter name, deliberately: the cmsRun command line,
# the config dump and edmProvDump then all show the same string, so a run log
# can be matched to a file's provenance by grep.
BOOLS = (
    "IoniExactDelta",
    "IoniKokoulin",
    "ReferenceChargeAware",
    "ReferenceSpeciesDedx",
    "ReferenceHadronRadiative",
    "ReferenceIonizationOnly",
    "IoniUrban2021",
    "DumpEmParameters",
    "EmHarmonise",
)
INTS = (
    "ReferenceSpeciesDedxNbin",
    "IoniKokoulinNbin",
)
FLOATS = (
    "IoniExactDeltaT0",
)

# DedxScale defaults to 1.0, not -1, so it needs its own "unset" sentinel
# rather than the <0 rule the others use.
FLOATS_POS = (
    "DedxScale",
)


def register(opts):
    for n in BOOLS + INTS:
        opts.register(n, -1, VarParsing.multiplicity.singleton,
                      VarParsing.varType.int,
                      f"Geant4ePropagator.{n} (-1 = leave the cfi default)")
    for n in FLOATS + FLOATS_POS:
        opts.register(n, -1.0, VarParsing.multiplicity.singleton,
                      VarParsing.varType.float,
                      f"Geant4ePropagator.{n} (<0 = leave the cfi default)")
    return opts


def apply(process, opts, propagator="Geant4ePropagator"):
    """Override only what was explicitly passed, and say so in the log."""
    prop = getattr(process, propagator)
    named = []
    for n in BOOLS:
        v = getattr(opts, n)
        if v >= 0:
            setattr(prop, n, cms.bool(bool(v)))
            named.append(f"{n}={bool(v)}")
    for n in INTS:
        v = getattr(opts, n)
        if v >= 0:
            setattr(prop, n, cms.int32(int(v)))
            named.append(f"{n}={int(v)}")
    for n in FLOATS + FLOATS_POS:
        v = getattr(opts, n)
        if v >= 0.:
            setattr(prop, n, cms.double(float(v)))
            named.append(f"{n}={float(v)}")
    print("[cvh] switches overridden on the command line: "
          + (", ".join(named) if named else "NONE (cfi defaults)"))
    # Always echo the resulting state: the point of moving off the environment
    # is that what ran is recoverable, and the log is the first place anyone
    # looks before edmProvDump.
    print("[cvh] effective: " + ", ".join(
        f"{n}={getattr(prop, n).value()}" for n in BOOLS + INTS + FLOATS + FLOATS_POS))
    return process
