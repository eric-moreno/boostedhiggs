from __future__ import print_function, division
import gzip
import json
import os

import uproot
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.styles.ROOT)
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
        #'text.usetex': False,
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

def drawSolo(h,sel,var_name,var_label,plottitle,lumifb,vars_cut,regionsel,savename):
    exceptions = ['process', var_name]
    for var,val in vars_cut.items():
        exceptions.append(var)
    if (regionsel is not ''):
        exceptions.append('region')
    print([ax.name for ax in h.axes()])
    x = h.sum(*[ax for ax in h.axes() if ax.name not in exceptions])
    for var,val in vars_cut.items():
        if var!=var_name:
            print('integrating ',var,val[0],val[1])
            x = x.integrate(var,slice(val[0],val[1]))
    for reg in regionsel:
        print('integrating ',reg)
        x = x.integrate('region',reg)
    if var_name in vars_cut.keys():
        x = x[:, vars_cut[var_name][0]:vars_cut[var_name][1]]

    xaxis = var_name
    x.axis(xaxis).label = var_label
    for ih,hkey in enumerate(x.identifiers('process')):
        x.identifiers('process')[ih].label = process_latex[hkey.name]

    x.axis('process').sorting = 'integral'
    fig,ax = plt.subplots()
    hist.plot1d(x,
                overlay='process',ax=ax,
                clear=False,
                fill_opts=fill_opts,
                error_opts=err_opts,
                )
    ax.autoscale(axis='x', tight=True)
    #ax.set_xlim(20, 200)
    ax.set_ylim(0, None)
    ax.ticklabel_format(axis='x', style='sci')
    old_handles, old_labels = ax.get_legend_handles_labels()
    leg = ax.legend(handles=old_handles,labels=old_labels,title=r'$%s$'%plottitle)
    lumi = plt.text(1., 1., r"%.1f fb$^{-1}$ (13 TeV)"%lumifb,fontsize=16,horizontalalignment='right',verticalalignment='bottom',transform=ax.transAxes)
    cmstext = plt.text(0., 1., "CMS",fontsize=20,horizontalalignment='left',verticalalignment='bottom',transform=ax.transAxes, fontweight='bold')
    addtext = plt.text(0.085, 1., "Simulation Preliminary",fontsize=16,horizontalalignment='left',verticalalignment='bottom',transform=ax.transAxes, style='italic')
    #hep.cms.cmslabel(ax, data=False, paper=False, year='2017')
    fig.savefig("solo_%s_%s_%s_lumi%i.pdf"%(sel,var_name,savename,lumifb))
    ax.semilogy()
    minvals = []
    for xd in x.values():
        minvals.append(min(np.trim_zeros(x.values()[xd])))
    decsplit = str(min(minvals)).split('.')
    if (int(decsplit[0])==0):
        logmin = 0.1**float(len(decsplit[1])-len(decsplit[1].lstrip('0'))+2)
    else:
        logmin = 10.**float(len(decsplit[0])-1)
    ax.set_ylim(logmin/10., None)
    fig.savefig("solo_%s_%s_%s_lumi%i_logy.pdf"%(sel,var_name,savename,lumifb))

def getPlots(args):
    print(args.lumi)
    lumifb = float(args.lumi)
    tag = args.tag
    savename = args.savetag

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
    # normalize to lumi
    for h in hists_mapped.values():
        h.scale({p: lumifb for p in h.identifiers('process')}, axis="process")
    

    # properties
    hist_name = args.hist
    var_name = args.var
    var_label = r"$%s$"%args.varlabel
    vars_cut =  {}
    #print(args.sel)
    if (len(args.sel)%3==0):
      for vi in range(int(len(args.sel)/3)):
        vars_cut[args.sel[vi*3]] = [float(args.sel[vi*3+1]), float(args.sel[vi*3+2])]
    print(vars_cut)
    h = hists_mapped[hist_name]
    print(h)
        
    drawSolo(h,args.hist,var_name,var_label,args.title,lumifb,vars_cut,args.regions,savename)

    os.chdir(pwd)

if __name__ == "__main__":
    #ex. python plot_solo.py --hists htt_test --tag test --var jet_pt --varlabel 'p_{T}(jet)' --title Test --lumi 41.5 --sel lep_pt 20. 200. --regions hadel_signal  --hist trigeff --savetag leppt_20
    parser = argparse.ArgumentParser()
    parser.add_argument('--hists',      dest='hists',    default="hists",      help="hists pickle name")
    parser.add_argument('--tag',        dest='tag',      default="",           help="tag")
    parser.add_argument('--savetag',    dest='savetag',  default="",           help="savetag")
    parser.add_argument('--var',        dest='var',      default="",           help="var")
    parser.add_argument('--varlabel',   dest='varlabel', default="",           help="varlabel")
    parser.add_argument('--title',      dest='title',    default="",           help="title")
    parser.add_argument('--lumi',       dest='lumi',     default=50.,          help="lumi",       type=float)
    parser.add_argument('--sel',        dest='sel',      default='',           help='selection',  nargs='+')
    parser.add_argument('--regions',    dest='regions',  default='',           help='regionsel',  nargs='+')
    parser.add_argument('--hist',       dest='hist',     default='',           help='histname')
    args = parser.parse_args()

    getPlots(args)
