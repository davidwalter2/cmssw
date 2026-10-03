#  DataFormats/SiStripCluster

## `SiStripApproximateClusterCollection`

The class `SiStripApproximateClusterCollection` is part of the RAW data, and any changes must be backwards compatible. In order to ensure it can be read by all future CMSSW releases, there is a `TestSiStripApproximateClusterCollection` unit test, which makes use of the `TestReadSiStripApproximateClusterCollection` analyzer and the `TestWriteSiStripApproximateClusterCollection` producer. The unit test checks that the object can be read properly from

* a file written by the same release
* files written by (some) earlier releases

If the persistent format of class `SiStripApproximateClusterCollection` gets changed in the future, please adjust the `TestReadSiStripApproximateClusterCollection` and `TestWriteSiStripApproximateClusterCollection` modules accordingly. It is important that every member container has some content in this test. Please also add new files to the [https://github.com/cms-data/DataFormats-SiStripCluster/](https://github.com/cms-data/DataFormats-SiStripCluster/) repository, and update the `TestSiStripApproximateClusterCollection` unit test to read the newly created files. The file name should contain the release or pre-release with which it was written and the split level.

## Reading SiStripCluster versions 10-13 (CMSSW <= 14_X) in 15_X (local W-mass patch)

The v10-v13 -> v14 read rule in `classes_def.xml` has no `source`: ROOT reads
`amplitudes_` and `firstStrip_` by automatic member-wise schema evolution (same
name and type on file in every version) and the rule only fills the
`barycenter_`/`charge_` cache that v14 keeps for every cluster (v12/v13 stored it
for approximate clusters only, v10/v11 not at all). The release rules sourced the
on-file members (`onfile.amplitudes_`, ...); ROOT does not serve those for
branches written at split level > 1 (root-project/root#19773), so every <= 14_X
file with the clusters at split 99 (AOD, RECO, ALCARECO, 2024 PromptReco MINIAOD)
was read as EMPTY clusters (firstStrip 0, no amplitudes, barycenter 0, charge 0) --
and the 15_0 re-MINI of 2024 (MiniAODv6 / MINIv6NANOv15) wrote those empty
clusters to disk. Not covered: the approximateMask flag of v12/v13 approximate
(rawprime) clusters is not set; their stored barycenter_/charge_ are kept.
