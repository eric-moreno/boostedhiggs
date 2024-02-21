#!/bin/bash

python3 -m pip install onnxruntime==1.12.0
python3 -m pip install correctionlib
#python3 -m pip install numpy==1.19.0
#python3 -m pip install --upgrade cython

# print the version numbers of everything in the environment
python -m pip freeze


# make dir for output
mkdir outfiles

# run code (e.g. for run.py)
python SCRIPTNAME --year YEAR --starti STARTNUM --endi ENDNUM --sample SAMPLE --processor PROCESSOR --condor --fileset metadata.json ADDOPTS

#move output to eos
xrdcp -f outfiles/*.hist EOSOUT
