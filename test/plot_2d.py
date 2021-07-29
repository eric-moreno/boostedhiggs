from __future__ import print_function, division
import gzip
import json
import os

import uproot
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.image as mplimg
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
line_opts = {
    'color': 'aquamarine',
    'linewidth':2,
    #'marker':'None',
    'linestyle':'dashed',
    'drawstyle':'steps'}

colormap = 'viridis'
overflow = 'all'

def drawStack(h,sel,var1_name,var1_label,var2_name,var2_label,plottitle,lumifb,vars_cut,regionsel,sample,savename,xlimits,ylimits,rebin1,rebin2):
    exceptions = ['process', var1_name, var2_name]
    for var,val in vars_cut.items():
        exceptions.append(var)
    if (regionsel is not ''):
        exceptions.append('region')
    print([ax.name for ax in h.axes()])
    x = h.sum(*[ax for ax in h.axes() if ax.name not in exceptions],overflow=overflow)
    for var,val in vars_cut.items():
        if var!=var1_name and var!=var2_name:
            print('integrating ',var,val[0],val[1])
            x = x.integrate(var,slice(val[0],val[1]),overflow=overflow)
    for reg in regionsel:
        x = x.integrate('region',reg)
    if var1_name in vars_cut.keys():
        x = x[:, vars_cut[var1_name][0]:vars_cut[var1_name][1]]
    if var2_name in vars_cut.keys():
        x = x[:, vars_cut[var2_name][0]:vars_cut[var2_name][1]]


    histo = x

    if (sample=='sig'):
        histo = x[nobkg]
    elif (sample=='bkg'):
        histo = x[nosig]
    elif (sample!='all'):
        histo = x[sample]

    histo = histo.sum(*[ax for ax in histo.axes() if ax.name is 'process'])

    xaxis = var1_name
    histo.axis(xaxis).label = var1_label
    yaxis = var2_name
    histo.axis(yaxis).label = var2_label

    print(histo.axis(xaxis).edges())
    print(histo.axis(yaxis).edges())

    if len(rebin1)==1:
        if (rebin1[0]>1): 
            histo = histo.rebin(xaxis,int(rebin1[0]))
    else:
        histo = histo.rebin(xaxis,hist.Bin(xaxis, var1_label, rebin1))

    if len(rebin2)==1:
        if (rebin2[0]>1): 
            histo = histo.rebin(yaxis,int(rebin2[0]))
    else:
        histo = histo.rebin(yaxis,hist.Bin(yaxis, var2_label, rebin2))

    fig,ax = plt.subplots()
    hist.plot2d(histo,
                ax=ax,
                clear=True,
                xaxis=var1_name,
                xoverflow=overflow,
                yoverflow=overflow,
                patch_opts={'cmap':colormap}
                )
    ax.autoscale(axis='x', tight=True)
    ax.autoscale(axis='y', tight=True)
    if len(xlimits)==2:
        try:
            ax.set_xlim(float(xlimits[0]), None)
        except:
            pass
        try:
            ax.set_xlim(None, float(xlimits[1]))
        except:
            pass
    if len(ylimits)==2:
        try:
            ax.set_ylim(float(ylimits[0]), None)
        except:
            pass
        try:
            ax.set_ylim(None, float(ylimits[1]))
        except:
            pass
    ax.ticklabel_format(axis='x', style='sci')
    ax.ticklabel_format(axis='y', style='sci')
    lumi = plt.text(1., 1., r"%.1f fb$^{-1}$ (13 TeV)"%lumifb,fontsize=16,horizontalalignment='right',verticalalignment='bottom',transform=ax.transAxes)
    cmstext = plt.text(0., 1., "CMS",fontsize=19,horizontalalignment='left',verticalalignment='bottom',transform=ax.transAxes, fontweight='bold')
    addtext = plt.text(0.105, 1., "Simulation Preliminary",fontsize=16,horizontalalignment='left',verticalalignment='bottom',transform=ax.transAxes, style='italic')
    #hep.cms.cmslabel(ax, data=False, paper=False, year='2017')
    fig.savefig("plot2d_%s_%s_vs_%s_%s_%s_lumi%i.pdf"%(sel,var1_name,var2_name,sample,savename,lumifb))
    fig,ax = plt.subplots()
    hist.plot2d(histo,
                ax=ax,
                clear=True,
                xaxis=var1_name,
                xoverflow=overflow,
                yoverflow=overflow,
                patch_opts={'cmap':colormap,'norm':colors.LogNorm()}
                )
    ax.autoscale(axis='x', tight=True)
    ax.autoscale(axis='y', tight=True)
    if len(xlimits)==2:
        try:
            ax.set_xlim(float(xlimits[0]), None)
        except:
            pass
        try:
            ax.set_xlim(None, float(xlimits[1]))
        except:
            pass
    if len(ylimits)==2:
        try:
            ax.set_ylim(float(ylimits[0]), None)
        except:
            pass
        try:
            ax.set_ylim(None, float(ylimits[1]))
        except:
            pass
    ax.ticklabel_format(axis='x', style='sci')
    ax.ticklabel_format(axis='y', style='sci')
    lumi = plt.text(1., 1., r"%.1f fb$^{-1}$ (13 TeV)"%lumifb,fontsize=16,horizontalalignment='right',verticalalignment='bottom',transform=ax.transAxes)
    cmstext = plt.text(0., 1., "CMS",fontsize=19,horizontalalignment='left',verticalalignment='bottom',transform=ax.transAxes, fontweight='bold')
    addtext = plt.text(0.105, 1., "Simulation Preliminary",fontsize=16,horizontalalignment='left',verticalalignment='bottom',transform=ax.transAxes, style='italic')
    fig.savefig("plot2d_%s_%s_vs_%s_%s_%s_lumi%i_logy.pdf"%(sel,var1_name,var2_name,sample,savename,lumifb))

    plt.close('all') 

def getPlots(args,returnHist=False):
    print(args.lumi)
    lumifb = float(args.lumi)
    tag = args.tag
    savename = args.savetag

    odir = 'plots/%s/'%tag
    os.system('mkdir -p %s'%odir)
    pwd = os.getcwd()

    hists_mapped = {}
    if args.dirname is not None:
        args.hists = [h[:-7] for h in glob.glob('%s/*.coffea'%args.dirname)]
    for h in args.hists:
        # open hists
        hists_unmapped = load('%s.coffea'%h)
        # map to hists
        for key, val in hists_unmapped.items():
            if isinstance(val, hist.Hist):
                if not args.noprocmap:
                    if key in hists_mapped:
                        hists_mapped[key] = hists_mapped[key].add(processmap.apply(val))
                    else:
                        hists_mapped[key] = processmap.apply(val)
                else:
                    hists_mapped[key] = val

    os.chdir(odir)
    # normalize to lumi
    for h in hists_mapped.values():
        if args.noprocmap:
            hp = hist.Cat("process", "process")
            h = h.group("dataset",hp,{d.name:[d.name] for d in h.axis("dataset").identifiers()})
        h.scale({p: lumifb for p in h.identifiers('process')}, axis="process")
    
    # properties
    hist_name = args.hist
    var1_name = args.var1
    var1_label = r"%s"%args.var1label
    var2_name = args.var2
    var2_label = r"%s"%args.var2label
    vars_cut =  {}
    #print(args.sel)
    if (len(args.sel)%3==0):
      for vi in range(int(len(args.sel)/3)):
        if (args.sel[vi*3+1]=='neginf'):
          vars_cut[args.sel[vi*3]] = [None, float(args.sel[vi*3+2])]
        elif (args.sel[vi*3+2]=='inf'):
          vars_cut[args.sel[vi*3]] = [float(args.sel[vi*3+1]), None]
        else:
          vars_cut[args.sel[vi*3]] = [float(args.sel[vi*3+1]), float(args.sel[vi*3+2])]
    print(vars_cut)
    h = hists_mapped[hist_name]
    print(h)
        
    if not returnHist:
        drawStack(h,args.hist,var1_name,var1_label,var2_name,var2_label,args.title,lumifb,vars_cut,args.regions,args.sample,savename,args.xlimits,args.ylimits,args.rebin1,args.rebin2)
        os.chdir(pwd)
        return None

    else:
        os.chdir(pwd)
        return h

if __name__ == "__main__":
    #ex. python plot_2d.py --hists ../condor/May09/hists_sum --tag May09 --var1 lep_pt --var1label '$p_{T}(\mu)$' --var2 lep_jet_dr --var2label '$\Delta R(j,\mu)$' --title '$\tau_{h}\mu$' --lumi 41.5 --regions hadmu_signal --hist lep_kin --sample sig --savetag hadmu
    parser = argparse.ArgumentParser()
    parser.add_argument('--hists',      dest='hists',     default="hists",      help="hists pickle name", nargs='+')
    parser.add_argument('--dirname',    dest='dirname',   default=None,         help="hists dir")
    parser.add_argument('--tag',        dest='tag',       default="",           help="tag")
    parser.add_argument('--savetag',    dest='savetag',   default="",           help="savetag")
    parser.add_argument('--var1',       dest='var1',      default="",           help="var1")
    parser.add_argument('--var1label',  dest='var1label', default="",           help="var1label")
    parser.add_argument('--rebin1',     dest='rebin1',    default=[1],          help='rebin1',     type=float,   nargs='+')
    parser.add_argument('--var2',       dest='var2',      default="",           help="var2")
    parser.add_argument('--var2label',  dest='var2label', default="",           help="var2label")
    parser.add_argument('--rebin2',     dest='rebin2',    default=[1],          help='rebin2',     type=float,   nargs='+')
    parser.add_argument('--title',      dest='title',     default="",           help="title")
    parser.add_argument('--lumi',       dest='lumi',      default=50.,          help="lumi",       type=float)
    parser.add_argument('--sel',        dest='sel',       default='',           help='selection',  nargs='+')
    parser.add_argument('--regions',    dest='regions',   default='',           help='regionsel',  nargs='+')
    parser.add_argument('--hist',       dest='hist',      default='',           help='histname')
    parser.add_argument('--sample',     dest='sample',    default='all',        help='sample')
    parser.add_argument('--xlimits',    dest='xlimits',   default='',           help='xlimits',    nargs='+')
    parser.add_argument('--ylimits',    dest='ylimits',   default='',           help='ylimits',    nargs='+')
    parser.add_argument('--noprocmap',  dest='noprocmap', action='store_true',  help='noprocmap')
    args = parser.parse_args()

    getPlots(args)
