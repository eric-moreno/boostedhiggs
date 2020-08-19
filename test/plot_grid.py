from __future__ import print_function, division
import gzip
import json
import os
import sys

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
line_opts = {
    'color': 'aquamarine',
    'linewidth':2,
    #'marker':'None',
    'linestyle':'dashed',
    'drawstyle':'steps'}

overflow = 'allnan'

def drawGrid(h,sel,var_name,var_label,plottitle,lumifb,vars_cut,regionsel,savename,gridsize,gridmod,grid_cut,grid_labels):
    print(gridsize,gridmod)
    exceptions = ['process', var_name]
    for var in vars_cut:
        exceptions.append(var)
    for var in grid_cut[0::3]:
        exceptions.append(var)
    if (regionsel is not ''):
        exceptions.append('region')
    print([ax.name for ax in h.axes()])
    x = h.sum(*[ax for ax in h.axes() if ax.name not in exceptions],overflow=overflow)
    for reg in regionsel:
        print('integrating ',reg)
        x = x.integrate('region',reg)
    for var,val in vars_cut.items():
        if var!=var_name:
            print('integrating ',var,val[0],val[1])
            x = x.integrate(var,slice(val[0],val[1]))
            #x = x.integrate(var,slice(val[0],val[1]),overflow=overflow)
    if var_name in vars_cut.keys():
        x = x[:, vars_cut[var_name][0]:vars_cut[var_name][1]]
    x_full = x

    fig,axs = plt.subplots(nrows=int((gridsize-1)/gridmod)+1, ncols=gridmod, squeeze=False, gridspec_kw={'wspace':0.4},figsize=(float(gridmod)*4.5, float(int((gridsize-1)/gridmod)+1)*6))
    print(axs)

    logmin = []
    print(grid_cut)
    for isel in range(gridsize):
        x = x_full.copy()
        if grid_cut[isel*3]!=var_name:
            print('integrating ',grid_cut[isel*3],grid_cut[isel*3+1],grid_cut[isel*3+2])
            x = x.integrate(grid_cut[isel*3],slice(grid_cut[isel*3+1],grid_cut[isel*3+2]))
        if var_name==grid_cut[isel*3]:
            x = x[:, grid_cut[isel*3+1]:grid_cut[isel*3+2]]
    
        xaxis = var_name
        x.axis(xaxis).label = var_label
        for ih,hkey in enumerate(x.identifiers('process')):
            x.identifiers('process')[ih].label = process_latex[hkey.name]
    
        x.axis('process').sorting = 'integral'
        the_axis = axs[int(isel/gridmod)][isel%gridmod]
        hist.plot1d(x[nosig],
                    overlay='process',ax=the_axis,
                    clear=False,
                    stack=True,
                    fill_opts=fill_opts,
                    error_opts=err_opts,
                    overflow=overflow
                    )
    
        all_bkg = 0
        for key,val in x[nosig].values().items():
            all_bkg+=val.sum()
        x_nobkg = x[nobkg]
        x_nobkg.scale(50)
    
        all_sig = 0
        for key,val in x_nobkg.values().items():
            all_sig +=val.sum()
        print('%.4f %.4f %.4f'%(all_bkg,all_sig,all_sig/math.sqrt(all_bkg)))
        hist.plot1d(x_nobkg,ax=the_axis,
                    overlay='process',
                    clear=False,
                    line_opts=line_opts,
                    overflow=overflow)
    
        the_axis.autoscale(axis='x', tight=True)
        #ax.set_xlim(20, 200)
        the_axis.set_ylim(0, None)
        the_axis.ticklabel_format(axis='x', style='sci')
        old_handles, old_labels = the_axis.get_legend_handles_labels()
        new_labels = []
        for xl in old_labels:
            if ('ggH' in xl): xl = xl + " (x 50)"
            new_labels.append(xl)
        try:
            title_add = ', %s'%grid_labels[isel]
        except:
            title_add = ''
        leg = the_axis.legend(handles=old_handles,labels=new_labels,title=r'%s%s'%(plottitle,title_add),frameon=True,framealpha=1.0,facecolor='white',fontsize=8)
        #lumi = plt.text(1., 1., r"%.1f fb$^{-1}$ (13 TeV)"%lumifb,fontsize=16,horizontalalignment='right',verticalalignment='bottom',transform=the_axis.transAxes)
        #cmstext = plt.text(0., 1., "CMS",fontsize=19,horizontalalignment='left',verticalalignment='bottom',transform=the_axis.transAxes, fontweight='bold')
        #addtext = plt.text(0.085, 1., "Simulation Preliminary",fontsize=16,horizontalalignment='left',verticalalignment='bottom',transform=the_axis.transAxes, style='italic')


        minvals = []
        for xd in x.values():
          if not np.trim_zeros(x.values()[xd]).any(): 
            continue
          minvals.append(min(np.trim_zeros(x.values()[xd])))
        if not minvals:
          decsplit = ['0','0']
        else:
          decsplit = str(min(minvals)).split('.')
        if (int(decsplit[0])==0):
          logmin.append(0.1**float(len(decsplit[1])-len(decsplit[1].lstrip('0'))+2))
        else:
          logmin.append(10.**float(len(decsplit[0])-1))

    #hep.cms.cmslabel(ax, data=False, paper=False, year='2017')
    fig.savefig("grid_%s_%s_%s_lumi%i.pdf"%(sel,var_name,savename,lumifb))
    for isel in range(gridsize):
        the_axis = axs[int(isel/gridmod)][isel%gridmod]
        the_axis.semilogy()
        the_axis.set_ylim(logmin[isel]/10., None)
    fig.savefig("grid_%s_%s_%s_lumi%i_logy.pdf"%(sel,var_name,savename,lumifb))

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
    var_label = r"%s"%args.varlabel
    vars_cut =  {}
    grid_cut = {}
    #print(args.sel)
    gridsize = 1
    if (len(args.sel)%3==0):
      for vi in range(int(len(args.sel)/3)):
        if (args.sel[vi*3+1]=='neginf'):
          vars_cut[args.sel[vi*3]] = [None, float(args.sel[vi*3+2])]
        elif (args.sel[vi*3+2]=='inf'):
          vars_cut[args.sel[vi*3]] = [float(args.sel[vi*3+1]), None]
        else:
          vars_cut[args.sel[vi*3]] = [float(args.sel[vi*3+1]), float(args.sel[vi*3+2])]
    grid_cut = []
    for vi in range(len(args.selgrid)):
      if (args.selgrid[vi]=='neginf' or args.selgrid[vi]=='inf'): 
          grid_cut.append(None)
      else:
          try: grid_cut.append(float(args.selgrid[vi]))
          except: grid_cut.append(str(args.selgrid[vi]))
    gridsize = int(len(grid_cut)/3)

    print(vars_cut)
    print(grid_cut)
    h = hists_mapped[hist_name]
    print(h)
    gridmod = args.gridmod
    if (args.gridmod==0):
        gridmod = gridsize
        
    if (gridsize<=1): print('plot_grid.py expects at least 2 grid selections, for a single selection use plot_stacks')
    elif (len(grid_cut)%3!=0): print('number of selgrid arguments must be divisible by 3: name lo hi')
    else: drawGrid(h,args.hist,var_name,var_label,args.title,lumifb,vars_cut,args.regions,savename,gridsize,gridmod,grid_cut,args.labelgrid)

    os.chdir(pwd)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--hists',      dest='hists',     default="hists",      help="hists pickle name")
    parser.add_argument('--tag',        dest='tag',       default="",           help="tag")
    parser.add_argument('--savetag',    dest='savetag',   default="",           help="savetag")
    parser.add_argument('--var',        dest='var',       default="",           help="var")
    parser.add_argument('--varlabel',   dest='varlabel',  default="",           help="varlabel")
    parser.add_argument('--title',      dest='title',     default="",           help="title")
    parser.add_argument('--lumi',       dest='lumi',      default=50.,          help="lumi",            type=float)
    parser.add_argument('--sel',        dest='sel',       default='',           help='selection',       nargs='+')
    parser.add_argument('--selgrid',    dest='selgrid',   default='',           help='grid selection',  nargs='+')
    parser.add_argument('--labelgrid',  dest='labelgrid', default='',           help='grid labels',     nargs='+')
    parser.add_argument('--regions',    dest='regions',   default='',           help='regionsel',       nargs='+')
    parser.add_argument('--hist',       dest='hist',      default='',           help='histname')
    parser.add_argument('--gridmod',    dest='gridmod',   default=0,            help='gridmod',         type=int)
    args = parser.parse_args()

    getPlots(args)
