#!/bin/sh
#
# Template of the shell script for submitting a CONDOR job
#
# Need to be substituted:
#    - C M S S W B A S E - local release base ($CMSSW_BASE)
#    - D I R E C T O R Y - condor work and output dir
#    - P R E F I X - some generic name for the set of jobs (like ttbar, wjets)
#    - J O B I D  - job number
#
#
#_____ setup the environment ____________________________________________
#

source /cvmfs/cms.cern.ch/cmsset_default.sh

#get setup coffea area
xrdcp root://cmseos.fnal.gov//store/user/drankin/dazsle_coffea.tgz ./dazsle_coffea.tgz
tar -zxf dazsle_coffea.tgz

tar -zxf coffeaenv.tar.gz
source coffeaenv/bin/activate
export PYTHONPATH=${PWD}:${PYTHONPATH}

mkdir test
cd test
xrdcp -f root://cmseos.fnal.gov//store/user/drankin/SCRIPTNAME .
python SCRIPTNAME --year YEAR --starti STARTNUM --endi ENDNUM --selsamples SAMPLE --outname htt_test

#move output to eos
xrdcp -f htt_test.coffea EOSOUT
