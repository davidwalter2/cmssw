#  DataFormats/SiStripCluster

## `SiStripApproximateClusterCollection`

The class `SiStripApproximateClusterCollection` is part of the RAW data, and any changes must be backwards compatible. In order to ensure it can be read by all future CMSSW releases, there is a `TestSiStripApproximateClusterCollection` unit test, which makes use of the `TestReadSiStripApproximateClusterCollection` analyzer and the `TestWriteSiStripApproximateClusterCollection` producer. The unit test checks that the object can be read properly from

* a file written by the same release
* files written by (some) earlier releases

If the persistent format of class `SiStripApproximateClusterCollection` gets changed in the future, please adjust the `TestReadSiStripApproximateClusterCollection` and `TestWriteSiStripApproximateClusterCollection` modules accordingly. It is important that every member container has some content in this test. Please also add new files to the [https://github.com/cms-data/DataFormats-SiStripCluster/](https://github.com/cms-data/DataFormats-SiStripCluster/) repository, and update the `TestSiStripApproximateClusterCollection` unit test to read the newly created files. The file name should contain the release or pre-release with which it was written and the split level.

## Reading SiStripCluster versions 10-13 (CMSSW <= 14_X) in 15_X (local W-mass patch)

The v10-v11 and v12-v13 -> v14 read rules in `classes_def.xml` are the release's, with
their `target` lists reordered so that the first target is the last source
(`firstStrip_` for v10-v11, `charge_` for v12-v13). ROOT attaches a rule to the
element of its first target and inserts it after the last source; when two or more
sources still exist under the same name and the first target is not the last
source, the rule is never executed for a class read member-wise from split
sub-branches. With the release order (`amplitudes_` first) every <= 14_X file with
the clusters at split level > 1 (AOD, RECO, ALCARECO, 2024 PromptReco MINIAOD) was
read as EMPTY clusters (firstStrip 0, no amplitudes, barycenter 0, charge 0) -- and
the 15_0 re-MINI of 2024 (MiniAODv6 / MINIv6NANOv15) wrote those empty clusters to
disk. Split level 1 was not affected. With the reordered targets, files written by
10_6, 13_0 and 14_0 at split level 99 read cluster for cluster the same as their
split-level-1 copies (detId, firstStrip, amplitudes, barycenter, charge). The rule
code is unchanged, so v12/v13 approximate (rawprime) clusters keep their stored
barycenter/charge and get the approximateMask flag.
