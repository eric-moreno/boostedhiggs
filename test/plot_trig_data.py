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
    #'linestyle':'-',
    'marker': '.',
    'markersize': 10.,
    'color':'k',
    'elinewidth': 1,
    'emarker': '-'
    }

#overflow_behavior = 'all'
overflow_behavior = 'over'

def drawTrigEff(h,var_name,var_label,vars_cut,num_sel,plot_title,plot_label,samplename):
    print(h)
    #print(h.values())
    #exceptions = [var_name,'dataset']
    exceptions = [var_name,'process']
    for var,val in vars_cut.items():
        exceptions.append(var)
    for var,val in num_sel.items():
        exceptions.append(var)
    print(exceptions)
    x = h.sum(*[ax for ax in h.axes() if ax.name not in exceptions],overflow='all')
    if var_name in num_sel.keys():
        print("%s cannot be a variable in numerator selection"%var_name)
        return
    for var,val in vars_cut.items():
        if var!=var_name:
            print('integrating ',var,val[0],val[1])
            if (len(val)==2):
                x = x.integrate(var,slice(val[0],val[1]))#,overflow=overflow_behavior)
            elif(type(val) is str):
                x = x.integrate(var,val)#,overflow=overflow_behavior)
            elif(len(val)==1):
                x = x.integrate(var,val[0])#,overflow=overflow_behavior)
    x_num = x.copy()
    #x_num = x_num.sum(*[ax for ax in x_num.axes() if ax.name in num_sel],overflow='all') #same axes as numerator
    #x_num.clear()

    xlist = []
    for var,val in num_sel.items():
        if var!=var_name:
            print('integrating ',var,val)
            print(var,val)
            if (len(val)==2):
                #xlist.append(x.integrate(var,slice(val[0],val[1])))#,overflow=overflow_behavior))
                x_num = x_num.integrate(var,slice(val[0],val[1]))#,overflow=overflow_behavior))
            elif(len(val)==1):
                #xlist.append(x.integrate(var,val[0]))#,overflow=overflow_behavior))
                x_num = x_num.integrate(var,val[0])#,overflow=overflow_behavior))
    #for xadd in xlist:
    #    x_num.add(xadd)
    x = x.sum(*[ax for ax in x.axes() if ax.name in num_sel],overflow='all')

    print(x.values())
    print(x_num.values())

    xaxis = var_name
    #x_num.axis(xaxis).label = "Efficiency"

    fig,ax = plt.subplots(1,1, figsize=(8,8))

    dataset = samplename
    if (dataset == ""):
        #dataset = x.axis("process").identifiers()[0]
        dataset = "data"
    histl = hist.Cat("histl", dataset)
    x_scale = x.group("process",histl,{"Total":dataset})
    x_num_scale = x_num.group("process",histl,{"Pass":dataset})
    total = list(x.sum(var_name).values().values())[0]
    maxval = list(x.values().values())[0].max()
    x_scale.scale(1./maxval)
    x_num_scale.scale(1./maxval)
    print(x_scale.values())
    print(x_num_scale.values())
    hist.plot1d(x_scale,
              overlay='histl',
              ax=ax,
              clear=False,
              fill_opts={'edgecolor':'k'},#,'facecolor':'r'},
              )
    hist.plot1d(x_num_scale,
              overlay='histl',
              ax=ax,
              clear=False,
              fill_opts={'edgecolor':'k'},#,'facecolor':'b'},
              )

    x = x.sum(*["process"],overflow='allnan')
    x_num = x_num.sum(*["process"],overflow='allnan')
    hist.plotratio(x_num,x,
                #overlay='dataset',
                ax=ax,
                clear=False,
                #stack=True,
                #denom_fill_opts=fill_opts,
                error_opts=err_opts,
                unc='clopper-pearson'
                )
    ax.autoscale(axis='x', tight=True)
    #ax.set_xlim(20, 200)
    ax.set_ylim(0, 1.3)
    ax.ticklabel_format(axis='x', style='sci')
    #old_handles, old_labels = ax.get_legend_handles_labels()
    #for x in old_labels:
    #    if ('H(125)' in x): x = x + " (x 50)"
    #leg = ax.legend(handles=old_handles,labels=old_labels,title='Hadronic trigger')
    ax.set_ylabel("Efficiency")
    #com_sample = plt.text(1., 1., r"(13 TeV)",fontsize=16,horizontalalignment='right',verticalalignment='bottom',transform=ax.transAxes)
    trigtitle = plt.text(0., 1., r"%s"%plot_title,fontsize=16,horizontalalignment='left',verticalalignment='bottom',transform=ax.transAxes)
    fig.savefig("ratio_data_%s_%s.pdf"%(plot_label,var_name))

def getPlots(args):
    tag = args.tag

    odir = 'plots/%s/'%tag
    os.system('mkdir -p %s'%odir)
    pwd = os.getcwd()

    # open hists
    hists_unmapped = load('%s.coffea'%args.hists)

    os.chdir(odir)

    # map to hists
    h_trig = None
    #for key in hists_unmapped:
    #    if (key==args.histname):
    #        h_trig = hists_unmapped[key]
    for key, val in hists_unmapped.items():
        if isinstance(val, hist.Hist):
            if key==args.histname:
                h_trig = processmap.apply(val)
    

    #print(h_trig)
    vars_cuts = {}
    ic = 0
    while ic < len(args.varcuts):
        if (args.varcuts[ic]=='region'):
            vars_cuts[args.varcuts[ic]] = args.varcuts[ic+1]
            ic = ic + 2
        else:
            vars_cuts[args.varcuts[ic]] = [float(args.varcuts[ic+1]) if args.varcuts[ic+1]!='None' else None, float(args.varcuts[ic+2]) if args.varcuts[ic+2]!='None' else None]
            ic = ic + 3
 
    num_sels = {}
    ic = 0
    while ic < len(args.numsel):
        if (args.numsel[ic]=='region'):
            num_sels[args.numsel[ic]] = args.numsel[ic+1]
            ic = ic + 2
        else:
            num_sels[args.numsel[ic]] = [float(args.numsel[ic+1]) if args.numsel[ic+1]!='None' else None, float(args.numsel[ic+2]) if args.numsel[ic+2]!='None' else None]
            ic = ic + 3
 
    drawTrigEff(h_trig,args.varname,args.varlabel,vars_cuts,num_sels,args.title,args.label,args.sample)

    os.chdir(pwd)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--hists',      dest='hists',    default="hists",      help="hists pickle name")
    parser.add_argument('--histname',   dest='histname', default="trigeff",    help="hist name")
    parser.add_argument('--tag',        dest='tag',      default="",           help="tag")
    parser.add_argument('--varname',    dest='varname',  default="",           help="varname")
    parser.add_argument('--varlabel',   dest='varlabel', default="",           help="varlabel")
    parser.add_argument('--varcuts',    dest='varcuts',  default="",           help="varcuts",    nargs='+')
    parser.add_argument('--numsel',     dest='numsel',   default="",           help="numsel",     nargs='+')
    parser.add_argument('--title',      dest='title',    default="",           help="title")
    parser.add_argument('--label',      dest='label',    default="",           help="label")
    parser.add_argument('--sample',     dest='sample',   default="",           help="sample")
    args = parser.parse_args()

    getPlots(args)
