import FWCore.ParameterSet.Config as cms


def setAlcaRecoSplitLevel(process):
    """Write ALCARECO PoolOutputModules with ROOT split level 1 instead of
    the PoolOutputModule default of 99.

    Fully-split (99) ALCARECO triggers a ROOT schema-evolution defect
    (I/O rules not applied for split branches, ROOT #19773): when the
    SiStripCluster collection -- whose persistent layout grew between
    CMSSW 10_6 (v11) and 15_0 (v14, CMSSW PR #47094) -- is read in
    CMSSW >= 15_0, the v11->v14 read rule does not fire and the clusters
    come back empty. Split level 1 keeps the collection items streamed
    together so the rule fires correctly. Size impact is ~1-2% for this
    skimmed tier. This mirrors MiniAOD, which already stores its
    slimmedMuonTrackExtras SiStripCluster at split level 1.

    Applied to every ALCARECO-tier PoolOutputModule so future re-RECO /
    ALCARECO production is forward-readable in 15_X without a repack.
    """
    for m in process.outputModules_().values():
        if (m.type_() == "PoolOutputModule"
                and hasattr(m, "dataset")
                and m.dataset.dataTier.value() == "ALCARECO"):
            m.splitLevel = cms.untracked.int32(1)
    return process
