from __future__ import print_function, division
import gzip
import json
import os

import uproot
import matplotlib.pyplot as plt
import numpy as np

from coffea import hist
from coffea.util import load

import pickle
import gzip
import math

import argparse
import processmap
from hists_map import *

def drawEfficiency(hnum,hdenom,axis_label,fig_lab='',title=''):
    print(fig_lab,title)
    data_err_opts = {
    'linestyle':'none',
    'marker': '.',
    'markersize': 10.,
    'color':'k',
    'elinewidth': 1,
    'emarker': '_'
    }
    axisDict = {axis_label}
    fig,ax,_ = hist.plotratio(num=hnum.sum(*[ax for ax in hnum.axes() if ax.name not in axisDict]),
                              denom=hdenom.sum(*[ax for ax in hdenom.axes() if ax.name not in axisDict]),
                              error_opts=data_err_opts,
                              guide_opts={},
                              unc='clopper-pearson')
    ax.set_ylim(0.2,1.2)
    ax.set_title(title)
    fig.savefig("triggereff_%s_%s.pdf"%(fig_lab,axis_label))

def getPlots(args):
    odir = 'plots/trigger/'
    os.system('mkdir -p %s'%odir)
    pwd = os.getcwd()

    # open hists
    hists_unmapped = load('%s.coffea'%args.hists)
    os.chdir(odir)
    print(hists_unmapped)

    # map to hists
    hists_mapped = {}
    for key, val in hists_unmapped.items():
        if isinstance(val, hist.Hist):
            hists_mapped[key] = processmap.apply(val)

    print(hists_mapped)
    # properties: denom and num
    trigger_plots = {'had_hadtrig': {'vars':['genH_pt','jet0_pt'],
                                     'denom':'hadseljet0Trignone',
                                     'num':'hadseljet0Trighad',
                                     'title':'HLT_PFHT1050 || HLT_AK8*_TrimMass50'},
                     'ele_hadtrig': {'vars':['genH_pt','genEle_pt'],
                                     'denom':'eleseljet0Trignone',
                                     'num':'eleseljet0Trighad',
                                     'title':'HLT_PFHT1050 || HLT_AK8*_TrimMass50'},
                     'ele_eletrig': {'vars':['genH_pt','genEle_pt'],
                                     'denom':'eleseljet0Trignone',
                                     'num':'eleseljet0Trigele',
                                     'title':'HLT_Ele115_CaloIdVT_GsfTrkIdT and others'},
                     'ele_vveletrig': {'vars':['genH_pt','genEle_pt'],
                                       'denom':'eleseljet0Trignone',
                                       'num':'eleseljet0Trigvvlele',
                                       'title':'HLT_Ele15_IsoVVVL_PFHT450',},
                     'mu_hadtrig': {'vars':['genH_pt','genMu_pt'],
                                     'denom':'museljet0Trignone',
                                     'num':'museljet0Trighad',
                                     'title':'HLT_PFHT1050 || HLT_AK8*_TrimMass50',},
                     'mu_mutrig': {'vars':['genH_pt','genMu_pt'],
                                   'denom':'museljet0Trignone',
                                   'num':'museljet0Trigmu',
                                   'title':'HLT_Mu50'},
                     'mu_vvmutrig': {'vars':['genH_pt','genMu_pt'],
                                     'denom':'museljet0Trignone',
                                     'num':'museljet0Trigvvlmu',
                                     'title':'HLT_Mu15_IsoVVVL_PFHT600'},
                     'mu_mettrig': {'vars':['genH_pt','MET_pt'],
                                    'denom':'museljet0Trignone',
                                    'num':'museljet0Trigmet',
                                    'title':'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight'},
                     'mu_btagtrig': {'vars':['genH_pt','genMu_pt'],
                                     'denom':'museljet0Trignone',
                                     'num':'museljet0Trigbtagmu',
                                     'title':'HLT_BTagMu_AK8Jet300_Mu5'},
                     }

    for key,item in trigger_plots.items():
        h_num = hists_mapped[item['num']]
        h_denom = hists_mapped[item['denom']]
        for var in item['vars']:
            drawEfficiency(h_num,h_denom,var,fig_lab=key+'_'+var,title=item['title'])

    os.chdir(pwd)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--hists',      dest='hists',    default="hists",  help="hists pickle name")
    parser.add_argument('--tag_num',    dest='tag_num',  default="",       help="tag numerator")
    parser.add_argument('--tag_denom',  dest='tag_denom',  default="",     help="tag denominator")
    parser.add_argument('--sel',        dest='sel',      default='',       help='')
    args = parser.parse_args()

    getPlots(args)
