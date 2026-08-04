# Rung-E private production helper: CASTOR geometry no longer exists in
# 15_X, so the 2016-era chain must not mix/digitize/reconstruct it.
import FWCore.ParameterSet.Config as cms


def customise_nocastor(process):
    if hasattr(process, "mix") and hasattr(process.mix, "mixObjects") \
            and hasattr(process.mix.mixObjects, "mixCH"):
        mixCH = process.mix.mixObjects.mixCH
        # input (VInputTag) and subdets (vstring) are PARALLEL lists --
        # filter them together or the module constructor reads garbage
        keep = [i for i, x in enumerate(mixCH.input)
                if "Castor" not in x.getProductInstanceLabel()]
        mixCH.input = cms.VInputTag(*[mixCH.input[i] for i in keep])
        mixCH.subdets = cms.vstring(*[mixCH.subdets[i] for i in keep])
        if hasattr(mixCH, "crossingFrames"):
            mixCH.crossingFrames = cms.untracked.vstring(
                *[s for s in mixCH.crossingFrames if "Castor" not in s])
    # the castor DIGITIZER lives as a sub-PSet of the MixingModule --
    # deleting sequence modules does not remove its consumes
    for pset_name in ("digitizers", "theDigitizersValid"):
        if hasattr(process, "mix") and hasattr(process.mix, pset_name):
            pset = getattr(process.mix, pset_name)
            for dname in list(pset.parameterNames_()):
                if "castor" in dname.lower():
                    delattr(pset, dname)
    # castor EDAliases (e.g. simCastorDigis -> mix) reference digitizer
    # products that no longer exist
    for aname in list(process.aliases_()):
        if "astor" in aname.lower():
            delattr(process, aname)
    # drop every castor module from every sequence/path it appears in
    for name in list(process.producers_()) + list(process.filters_()) \
            + list(process.analyzers_()):
        if "astor" in name.lower():
            mod = getattr(process, name)
            for pname, path in list(process.paths_().items()) \
                    + list(process.endpaths_().items()):
                path.remove(mod)
            for sname in list(process.sequences_()):
                getattr(process, sname).remove(mod)
            delattr(process, name)
    # keep the output content free of castor products
    for out in ("RECOSIMoutput", "RAWSIMoutput", "FEVTDEBUGoutput"):
        if hasattr(process, out):
            getattr(process, out).outputCommands.append(
                "drop *_*astor*_*_*")
    return process
