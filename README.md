# boostedhiggs

## Quickstart
```bash
# check your platform: CC7 shown below, for SL6 it would be "x86_64-slc6-gcc8-opt"
source /cvmfs/sft.cern.ch/lcg/views/LCG_96python3/x86_64-centos7-gcc8-opt/setup.sh  # or .csh, etc.
git clone git@github.com:drankincms/boostedhiggs.git
cd boostedhiggs
pip install --user --editable .
```

## For Condor:
# Only necessary once
```
bash envSetup.sh
```
# For submission
```
cd condor
python CoffeaSubmit.py Apr24 run_htt.py 50 1
```
The arguments are: 
1) label (and output directory)
2) script to use (searches in test/ directory)
3) number of files per job 
4) whether a re-tar of environment/processors is necessary (usually is)
