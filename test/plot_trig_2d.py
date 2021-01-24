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
#import processmap
#from hists_map import *

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

def drawTrigEff(h,var1_name,var1_label,var2_name,var2_label,vars_cut,num_sel,plot_title,plot_label,sample,debug=False):
    print(h)
    #print(h.values())
    exceptions = [var1_name,var2_name,'dataset']
    for var,val in vars_cut.items():
        exceptions.append(var)
    for var,val in num_sel.items():
        exceptions.append(var)
    print(exceptions)
    x = h.sum(*[ax for ax in h.axes() if ax.name not in exceptions],overflow='allnan')
    if var1_name in num_sel.keys() or var2_name in num_sel.keys():
        print("%s and %s cannot be a variable in numerator selection"%(var1_name,var2_name))
        return
    for var,val in vars_cut.items():
        if var!=var1_name and var!=var2_name:
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
        if var!=var1_name and var!=var2_name:
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
    x = x.sum(*[ax for ax in x.axes() if ax.name in num_sel],overflow='allnan')

    #print(x.values())
    #print(x_num.values())

    fig,ax = plt.subplots(1,1, figsize=(8,8))

    histl = hist.Cat("histl", sample)
    x = x.sum(*["dataset"],overflow='allnan')
    x_num = x_num.sum(*["dataset"],overflow='allnan')

    x_bins = x.axis(var1_name).edges()
    y_bins = x.axis(var2_name).edges()

    den_arr = np.array(x.values(overflow='all')[()])
    num_arr = np.array(x_num.values(overflow='all')[()])

    if ([ax.name for ax in x.axes()][0]==var1_name):
        den_arr = np.transpose(den_arr)
        num_arr = np.transpose(num_arr)

    den_arr[:][1] = den_arr[:][1] + den_arr[:][0]
    den_arr[:][-2] = den_arr[:][-2] + den_arr[:][-1]
    num_arr[:][1] = num_arr[:][1] + num_arr[:][0]
    num_arr[:][-2] = num_arr[:][-2] + num_arr[:][-1]

    den_arr[1][:] = den_arr[1][:] + den_arr[0][:]
    den_arr[-2][:] = den_arr[-2][:] + den_arr[-1][:]
    num_arr[1][:] = num_arr[1][:] + num_arr[0][:]
    num_arr[-2][:] = num_arr[-2][:] + num_arr[-1][:]

    den_arr = np.delete(den_arr,-1,0)
    den_arr = np.delete(den_arr,0,0)
    den_arr = np.delete(den_arr,-1,1)
    den_arr = np.delete(den_arr,0,1)
    num_arr = np.delete(num_arr,-1,0)
    num_arr = np.delete(num_arr,0,0)
    num_arr = np.delete(num_arr,-1,1)
    num_arr = np.delete(num_arr,0,1)

    eff_range_arr = hist.clopper_pearson_interval(num_arr, den_arr)

    #print(eff_range_arr[0])
    #print(num_arr/den_arr)
    #print(eff_range_arr[1])
    #print(x_bins)
    #print(y_bins)

    #print(num_arr)
    #print(den_arr)
    #print(x_bins)
    #print(y_bins)

    plt.hist2d([(x_bins[ix]+x_bins[ix+1])/2. for ix in range(len(x_bins)-1) for iy in range(len(y_bins)-1)], 
               [(y_bins[iy]+y_bins[iy+1])/2. for ix in range(len(x_bins)-1) for iy in range(len(y_bins)-1)], 
               bins=[x_bins,y_bins], weights=np.divide(num_arr.flatten('F'),den_arr.flatten('F'), out=np.zeros_like(num_arr.flatten('F')), where=den_arr.flatten('F')!=0))

    for i in range(len(y_bins)-1):
        for j in range(len(x_bins)-1):
            ax.text((x_bins[j]+x_bins[j+1])/2.,(y_bins[i]+y_bins[i+1])/2., "{:0.2f}".format(num_arr[i,j]/den_arr[i,j]) if den_arr[i,j]>0. else "", 
                color="k", ha="center", va="center")#, fontweight="bold")

    ax.autoscale(axis='x', tight=True)
    ax.autoscale(axis='y', tight=True)
    ax.ticklabel_format(axis='x', style='sci')
    #old_handles, old_labels = ax.get_legend_handles_labels()
    #for x in old_labels:
    #    if ('ggh' in x): x = x + " (x 50)"
    #leg = ax.legend(handles=old_handles,labels=old_labels,title='Hadronic trigger')
    ax.set_xlabel(var1_label)
    ax.set_ylabel(var2_label)
    com_sample = plt.text(1., 1., r"(13 TeV)",fontsize=16,horizontalalignment='right',verticalalignment='bottom',transform=ax.transAxes)
    trigtitle = plt.text(0., 1., r"%s"%plot_title,fontsize=16,horizontalalignment='left',verticalalignment='bottom',transform=ax.transAxes)
    fig.savefig("ratio2d_%s_%s_%s.pdf"%(plot_label,var1_name,var2_name))

    if (debug):
        plt.clf()
        fig,ax = plt.subplots(1,1, figsize=(8,8))
        plt.hist2d([x_bins[ix] for ix in range(len(x_bins)-1) for iy in range(len(y_bins)-1)], 
               [y_bins[iy] for ix in range(len(x_bins)-1) for iy in range(len(y_bins)-1)], 
               bins=[x_bins,y_bins], weights=num_arr.flatten('F'))
        ax.autoscale(axis='x', tight=True)
        ax.autoscale(axis='y', tight=True)
        ax.ticklabel_format(axis='x', style='sci')
        ax.set_xlabel(var1_label)
        ax.set_ylabel(var2_label)
        com_sample = plt.text(1., 1., r"(13 TeV)",fontsize=16,horizontalalignment='right',verticalalignment='bottom',transform=ax.transAxes)
        trigtitle = plt.text(0., 1., r"%s"%plot_title,fontsize=16,horizontalalignment='left',verticalalignment='bottom',transform=ax.transAxes)
        for i in range(len(y_bins)-1):
            for j in range(len(x_bins)-1):
                ax.text((x_bins[j]+x_bins[j+1])/2.,(y_bins[i]+y_bins[i+1])/2., str(int(num_arr[i,j])), 
                    color="k", ha="center", va="center")#, fontweight="bold")
        fig.savefig("ratio2d_debug_%s_%s_%s.pdf"%(plot_label,var1_name,var2_name))

def getPlots(args):
    tag = args.tag

    odir = 'plots/%s/'%tag
    os.system('mkdir -p %s'%odir)
    pwd = os.getcwd()

    h_trig = None
    for hs in args.hists:
        # open hists
        hists_unmapped = load('%s.coffea'%hs)

        # map to hists
        for key in hists_unmapped:
            if (key==args.histname):
                if not h_trig: h_trig = hists_unmapped[key]
                else: h_trig = h_trig + hists_unmapped[key]
    
    os.chdir(odir)

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
 
    drawTrigEff(h_trig,args.var1name,args.var1label,args.var2name,args.var2label,vars_cuts,num_sels,args.title,args.label,args.sample,args.debug)

    os.chdir(pwd)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--hists',      dest='hists',    default=["hists"],    help="hists pickle name",   nargs='+')
    parser.add_argument('--histname',   dest='histname',  default="trigeff",    help="hist name")
    parser.add_argument('--tag',        dest='tag',       default="",           help="tag")
    parser.add_argument('--var1name',   dest='var1name',  default="",           help="var1name")
    parser.add_argument('--var1label',  dest='var1label', default="",           help="var1label")
    parser.add_argument('--var2name',   dest='var2name',  default="",           help="var2name")
    parser.add_argument('--var2label',  dest='var2label', default="",           help="var2label")
    parser.add_argument('--varcuts',    dest='varcuts',   default="",           help="varcuts",    nargs='+')
    parser.add_argument('--numsel',     dest='numsel',    default="",           help="numsel",     nargs='+')
    parser.add_argument('--title',      dest='title',     default="",           help="title")
    parser.add_argument('--label',      dest='label',     default="",           help="label")
    parser.add_argument('--sample',     dest='sample',    default="",           help="sample")
    parser.add_argument('--debug',      dest='debug',     action='store_true',  help="debug")
    args = parser.parse_args()

    getPlots(args)
