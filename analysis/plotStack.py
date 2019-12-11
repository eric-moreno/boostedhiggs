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

plt.rcParams.update({
        'font.size': 14,
        'axes.titlesize': 18,
        'axes.labelsize': 18,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'text.usetex': False,
        })

fill_opts = {
    'edgecolor': (0,0,0,0.3),
    'alpha': 0.8
    }
err_opts = {
    'label':'Stat. Unc.',
    'hatch':'///',
    'facecolor':'none',
    'edgecolor':(0,0,0,.5),
    'linewidth': 0
    }

def drawStack(h,var_name,var_label,lumifb,vars_cut):
    exceptions = ['process', var_name]
    for var,val in vars_cut.items():
        exceptions.append(var)
    x = h.sum(*[ax for ax in h.axes() if ax.name not in exceptions])
    for var,val in vars_cut.items():
        if var!=var_name:
            print('integrating ',var,val[0],val[1])
            x = x.integrate(var,slice(val[0],val[1]))
    if var_name in vars_cut.keys():
        x = x[:, vars_cut[var_name][0]:vars_cut[var_name][1]]

    xaxis = var_name
    x.axis(xaxis).label = var_label
    for ih,hkey in enumerate(x.identifiers('process')):
        x.identifiers('process')[ih].label = process_latex[hkey.name]

    x.axis('process').sorting = 'integral'
    fig,ax = plt.subplots(1,1, figsize=(8,8))
    hist.plot1d(x[nosig],
                overlay='process',ax=ax,
                clear=False,
                stack=True,
                fill_opts=fill_opts,
                error_opts=err_opts,
                )
    all_bkg = 0
    for key,val in x[nosig].values().items():
        all_bkg+=val.sum()
    x_nobkg = x[nobkg]
    x_nobkg.scale(0.29)
    all_sig = 0
    for key,val in x_nobkg.values().items():
        all_sig +=val.sum()
    x_nobkg.scale(50)
    print('%.4f %.4f %.4f'%(all_bkg,all_sig,all_sig/math.sqrt(all_bkg)))
    hist.plot1d(x_nobkg,ax=ax,overlay='process',clear=False,error_opts={'color': 'aquamarine','linewidth':2})
    ax.autoscale(axis='x', tight=True)
    #ax.set_xlim(20, 200)
    ax.set_ylim(0, None)
    ax.ticklabel_format(axis='x', style='sci')
    old_handles, old_labels = ax.get_legend_handles_labels()
    leg = ax.legend(handles=old_handles,labels=old_labels,title='Hadronic trigger')
    lumi = plt.text(1., 1., r"%i fb$^{-1}$ (13 TeV)"%lumifb,fontsize=16,horizontalalignment='right',verticalalignment='bottom',transform=ax.transAxes)
    fig.savefig("stack_%s_lumi%i.pdf"%(var_name,lumifb))

def getPlots(args):
    lumifb = args.lumi
    tag = args.tag

    odir = 'plots/%s/'%tag
    os.system('mkdir -p %s'%odir)
    pwd = os.getcwd()

    # open hists
    hists_unmapped = load('%s.coffea'%args.hists)
    os.chdir(odir)

    # map to hists
    hists_mapped = {}
    for key, val in hists_unmapped.items():
        if isinstance(val, hist.Hist):
            hists_mapped[key] = processmap.apply(val)
    print(hists_mapped)
    # normalize to lumi
    for h in hists_mapped.values():
        h.scale({p: lumifb for p in h.identifiers('process')}, axis="process")

    # properties
    #hist_name = 'jet_museljet'
    hist_name = 'jet_museljet'
    var_name = 'jetMu_pt'
    var_label = r'Jet p$_{T}$ (GeV)'
    vars_cut =  {'jetMu_lsf3':[0.78,1],
                 'jetMu_msd':[20,200],
                 }
    h = hists_mapped[hist_name]
        
    drawStack(h,var_name,var_label,lumifb,vars_cut)

    os.chdir(pwd)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--hists',      dest='hists',    default="hists",  help="hists pickle name")
    parser.add_argument('--tag',        dest='tag',      default="",       help="tag")
    parser.add_argument('--lumi',       dest='lumi',     default=50,       help="lumi")
    parser.add_argument('--sel',        dest='sel',      default='',       help='')
    args = parser.parse_args()

    getPlots(args)
