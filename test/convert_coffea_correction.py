from __future__ import print_function, division
import gzip
import json
import os

import uproot
import matplotlib.pyplot as plt
import numpy as np

from coffea import hist
from coffea.util import load, save
from coffea.lookup_tools import extractor

import pickle
import gzip
import math

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--filename',     dest='filename',   default="../boostedhiggs/data/trig_sf_corr",          help="filename")
    args = parser.parse_args()

    ext = extractor()
    ext.add_weight_sets(["* * %s.json"%args.filename])
    ext.finalize()

    evaluator = ext.make_evaluator()
    for key in evaluator.keys():
        print("\t", key)
        print("\t\t",evaluator[key])
    
    save(evaluator,"%s.coffea"%args.filename)
