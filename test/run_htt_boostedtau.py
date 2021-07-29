import os
import numpy as np
from coffea import processor, util, hist
import json

from boostedhiggs import HtautauProcessor_BoostedTau, BTagEfficiency
from coffea.nanoaod import NanoEvents

import argparse

def run_processor(year,selsamples,starti,endi,outname):
    p = HtautauProcessor_BoostedTau(year=year)
    
    files = {}
    
    with open('../data/filesetBoostedTau.json', 'r') as f:
        newfiles = json.load(f)
        files.update(newfiles)
    
    selfiles = {k: files[k][starti:endi] for k in selsamples}

    args = {'nano': True, 'workers': 1, 'savemetrics': True}
    out, metrics = processor.run_uproot_job(selfiles, 'Events', p, processor.futures_executor, args)
    
    util.save(out, '%s.coffea'%outname)


if __name__ == "__main__":
    #python run_htt_boostedtau.py --year 2017 --selsamples boostedTau_GluGluHTauTau_boostedTaua_13TeV_user --starti 0 --endi 2 --outname htt_boostedtau
    parser = argparse.ArgumentParser()
    parser.add_argument('--year',       dest='year',       default='2017',       help="year",        type=str)
    parser.add_argument('--starti',     dest='starti',     default=0,            help="start index", type=int)
    parser.add_argument('--endi',       dest='endi',       default=-1,           help="end index",   type=int)
    parser.add_argument('--selsamples', dest='selsamples', default=[],           help='selsamples',  nargs='+')
    parser.add_argument('--outname',    dest='outname',    default='htt_test',   help='outname')
    args = parser.parse_args()

    run_processor(args.year,args.selsamples,args.starti,args.endi,args.outname)

    possible = [
      "boostedTau_GluGluHTauTau_boostedTaua_13TeV_user",
      "boostedTau_QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8",
      "boostedTau_WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8",
    ]
        #TTToHadronic_TuneCP5_13TeV-powheg-pythia8 TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8
        #WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8 WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8 WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8 WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8 WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8
        #QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8 QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8 QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8 QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8 QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8 QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8
#python run_htt_lepid.py --year 2017 --starti 0 --endi 100 --selsamples --outname htt_runtest_lepid_QCD

    #print(possible)
#    python run_htt.py --year 2017 --starti 0 --endi 100 --selsamples DYJetsToLL_M-50_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8 --outname htt_runtest
