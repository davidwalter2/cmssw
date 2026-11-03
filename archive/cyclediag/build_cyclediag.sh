#!/bin/bash
cd /work/submit/david_w/ZMass/CMSSW_15_0_19_patch2_cyclediag
source /cvmfs/cms.cern.ch/cmsset_default.sh
scram b ProjectRename > build.log 2>&1
cd src
eval $(scramv1 runtime -sh)
scram b -j 32 >> ../build.log 2>&1
echo "EXIT $?" >> ../build.log
