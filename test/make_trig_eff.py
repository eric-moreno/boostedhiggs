from __future__ import print_function, division
import gzip
import json
import os

import uproot
import matplotlib.pyplot as plt
import numpy as np

from coffea import hist
from coffea.util import load, save

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

chanlist = ["hadhad", "hadel", "hadmu"]
histnames = {
    'hadhad': "trigeff_m",
    'hadel': "trigeff_m",
    'hadmu': "trigeff_m",
}
varcuts_data = {
    "hadhad": {"region": "hadhad_signal_150", "trig_pass_ref": [0.5, 1.5]},
    "hadel": {"region": "hadel_signal_150", "trig_pass_ref": [0.5, 1.5], "jet_msd": [40.,None]},
    "hadmu": {"region": "hadmu_signal_150", "trig_pass_ref": [0.5, 1.5], "jet_msd": [40.,None]},
}
varcuts_mc = {
    "hadhad": {"region": "hadhad_signal_150", "trig_pass_ref": [0.5, 1.5]},
    "hadel": {"region": "hadel_signal_150", "trig_pass_ref": [0.5, 1.5], "jet_msd": [40.,None]},
    "hadmu": {"region": "hadmu_signal_150", "trig_pass_ref": [0.5, 1.5], "jet_msd": [40.,None]},
}
var1names = {
    "hadhad": "jet_pt",
    "hadel": "jet_pt",
    "hadmu": "jet_pt",
}
var1labels = {
    "hadhad": "$p_{T}(jet)$",
    "hadel": "$p_{T}(jet)$",
    "hadmu": "$p_{T}(jet)$",
}
rebin1 = {
    "hadhad": [300.,350.,400.,450.,500.,550.,600.,650.,700.],
    "hadel": [300.,350.,400.,500.,700.],
    "hadmu": [300.,350.,400.,500.,700.],
}
var2names = {
    "hadhad": "jet_msd",
    "hadel": "lep_pt",
    "hadmu": "lep_pt",
}
var2labels = {
    "hadhad": "$m_{SD}(jet)$",
    "hadel": "$p_{T}(e)$",
    "hadmu": "$p_{T}(\mu)$",
}
rebin2 = {
    "hadhad": [0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.],
    "hadel": [20.,32.,44.,56.,68.,92.,116.,140.],
    "hadmu": [20.,32.,44.,56.,68.,92.,116.,140.],
}
numsels = {
    "hadhad": {"trig_pass_hadhad": [0.5, 1.5]},
    "hadel": {"trig_pass_hadel": [0.5, 1.5]},
    "hadmu": {"trig_pass_hadmu": [0.5, 1.5]},
}

#overflow_behavior = 'all'
overflow_behavior = 'over'

def getTrigEff(h,var1_name,var2_name,vars_cut,num_sel,rebins1,rebins2):
#def drawTrigEff(h,var1_name,var1_label,var2_name,var2_label,vars_cut,num_sel,plot_title,plot_label):
    print(h)
    #print(h.values())
    exceptions = [var1_name,var2_name,'dataset']
    for var,val in vars_cut.items():
        exceptions.append(var)
    for var,val in num_sel.items():
        exceptions.append(var)
    print(exceptions)
    x = h.sum(*[ax for ax in h.axes() if ax.name not in exceptions],overflow='all')
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
    x = x.sum(*[ax for ax in x.axes() if ax.name in num_sel],overflow='none')

    #print(x.values())
    #print(x_num.values())

    x = x.sum(*["dataset"],overflow='allnan')
    x_num = x_num.sum(*["dataset"],overflow='allnan')

    #x = x.rebin(var1_name, hist.Bin(var1_name+"_new", var1_name+"_new", rebins1))
    #x = x.rebin(var2_name, hist.Bin(var2_name+"_new", var2_name+"_new", rebins2))
    #x_num = x_num.rebin(var1_name, hist.Bin(var1_name+"_new", var1_name+"_new", rebins1))
    #x_num = x_num.rebin(var2_name, hist.Bin(var2_name+"_new", var2_name+"_new", rebins2))
    x = x.rebin(var1_name, hist.Bin(var1_name, var1_name, rebins1))
    x = x.rebin(var2_name, hist.Bin(var2_name, var2_name, rebins2))
    x_num = x_num.rebin(var1_name, hist.Bin(var1_name, var1_name, rebins1))
    x_num = x_num.rebin(var2_name, hist.Bin(var2_name, var2_name, rebins2))

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

    print(num_arr)
    print(den_arr)
    print(x_bins)
    print(y_bins)

    eff_range_arr = hist.clopper_pearson_interval(num_arr, den_arr)

    return np.transpose(np.divide(num_arr,den_arr, out=np.zeros_like(num_arr), where=den_arr!=0)),x_bins,y_bins,np.transpose(eff_range_arr,[0,2,1])

def getHists(filename_data,filename_mc,hadel_w,hadmu_w,hadhad_w):

    eff_hists_data = {}
    eff_hists_mc = {}
    eff_hists_data_int = {}
    eff_hists_mc_int = {}
    x_bins = {}
    y_bins = {}
    chan_w = {'hadhad':hadhad_w,'hadel':hadel_w,'hadmu':hadmu_w}

    for chan in chanlist:
        h_trig = None
        for f_d in filename_data:
            # open hists
            hists_unmapped_data = load('%s.coffea'%f_d)
    
            # map to hists
            for key in hists_unmapped_data:
                if (key==histnames[chan]):
                    if not h_trig:
                        h_trig = hists_unmapped_data[key]
                    else:
                        h_trig = h_trig + hists_unmapped_data[key]
            
            eff_hists_data[chan],_,_,eff_hists_data_int[chan] = getTrigEff(h_trig,var1names[chan],var2names[chan],varcuts_data[chan],numsels[chan],rebin1[chan],rebin2[chan])
            #drawTrigEff(h_trig,args.var1name,args.var1label,args.var2name,args.var2label,vars_cuts,num_sels,args.title,args.label)

        h_trig = None
        if (len(chan_w[chan]) != len(filename_mc)):
            chan_w[chan] = [1. for f in filename_mc]
        for i,f_m in enumerate(filename_mc):
            # open hists
            hists_unmapped_mc = load('%s.coffea'%f_m)

            # map to hists
            for key in hists_unmapped_mc:
                if (key==histnames[chan]):
                    if (chan_w[chan][i] != 1.): 
                        hists_unmapped_mc[key].scale(chan_w[chan][i])
                    if not h_trig:
                        h_trig = hists_unmapped_mc[key]
                    else:
                        h_trig = h_trig + hists_unmapped_mc[key]
            
            eff_hists_mc[chan],x_bins[chan],y_bins[chan],eff_hists_mc_int[chan] = getTrigEff(h_trig,var1names[chan],var2names[chan],varcuts_mc[chan],numsels[chan],rebin1[chan],rebin2[chan])
            #drawTrigEff(h_trig,args.var1name,args.var1label,args.var2name,args.var2label,vars_cuts,num_sels,args.title,args.label)

    return eff_hists_data,eff_hists_mc,x_bins,y_bins,eff_hists_data_int,eff_hists_mc_int
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--hists_data',      dest='hists_data',    default="hists_data",      help="hists pickle name (data)",    nargs='+')
    parser.add_argument('--hists_mc',        dest='hists_mc',      default="hists_mc",        help="hists pickle name (MC)",    nargs='+')
    parser.add_argument('--hadel_w',         dest='hadel_w',       default=[1.],         help="HadEl File Weights (MC)",   nargs='+', type=float)
    parser.add_argument('--hadmu_w',         dest='hadmu_w',       default=[1.],         help="HadMu File Weights (MC)",   nargs='+', type=float)
    parser.add_argument('--hadhad_w',        dest='hadhad_w',      default=[1.],        help="HadHad File Weights (MC)",  nargs='+', type=float)
    #parser.add_argument('--histname',   dest='histname', default="trigeff",    help="hist name")
    parser.add_argument('--tag',        dest='tag',      default="trig_sf_debug",             help="tag")
    parser.add_argument('--output',     dest='output',   default="../boostedhiggs/data/trig_sf_corr",          help="output")
    parser.add_argument('--year',       dest='year',     default="2017",             help="year")
    #parser.add_argument('--varname',    dest='varname',  default="",           help="varname")
    #parser.add_argument('--varlabel',   dest='varlabel', default="",           help="varlabel")
    #parser.add_argument('--varcuts',    dest='varcuts',  default="",           help="varcuts",    nargs='+')
    #parser.add_argument('--numsel',     dest='numsel',   default="",           help="numsel",     nargs='+')
    #parser.add_argument('--title',      dest='title',    default="",           help="title")
    #parser.add_argument('--label',      dest='label',    default="",           help="label")
    args = parser.parse_args()

    #python make_trig_eff.py --hists_data ../condor/May27_Trig/hists_trig_Run2017CDEF --hists_mc ../condor/May27_Trig/hists_trig_QCD

    eff_hists_data,eff_hists_mc,x_bins,y_bins,eff_hists_data_int,eff_hists_mc_int = getHists(args.hists_data,args.hists_mc,args.hadel_w,args.hadmu_w,args.hadhad_w)

    h_trig_sf = {}
    arr_sf = {}
    arr_sf_up = {}
    arr_sf_down = {}
    h_trig_eff_mc = {}
    h_trig_eff_data = {}
    for chan in chanlist:
         h_trig_sf[args.year+"_trigsf_"+chan+"_nom"] = hist.Hist("Trigger Scale Factor (%s) Nominal"%chan,
                  hist.Bin(var1names[chan], var1labels[chan], x_bins[chan]),
                  hist.Bin(var2names[chan], var2labels[chan], y_bins[chan]),
                  )
         h_trig_sf[args.year+"_trigsf_"+chan+"_up"] = hist.Hist("Trigger Scale Factor (%s) Up"%chan,
                  hist.Bin(var1names[chan], var1labels[chan], x_bins[chan]),
                  hist.Bin(var2names[chan], var2labels[chan], y_bins[chan]),
                  )
         h_trig_sf[args.year+"_trigsf_"+chan+"_down"] = hist.Hist("Trigger Scale Factor (%s) Down"%chan,
                  hist.Bin(var1names[chan], var1labels[chan], x_bins[chan]),
                  hist.Bin(var2names[chan], var2labels[chan], y_bins[chan]),
                  )
         h_trig_eff_mc[chan] = hist.Hist("Trigger Efficiency, MC (%s)"%chan,
                  hist.Bin(var1names[chan], var1labels[chan], x_bins[chan]),
                  hist.Bin(var2names[chan], var2labels[chan], y_bins[chan]),
                  )
         h_trig_eff_data[chan] = hist.Hist("Trigger Efficiency, Data (%s)"%chan,
                  hist.Bin(var1names[chan], var1labels[chan], x_bins[chan]),
                  hist.Bin(var2names[chan], var2labels[chan], y_bins[chan]),
                  )
         inputs = {}
         inputs_up = {}
         inputs_down = {}
         inputs_mc = {}
         inputs_data = {}
         inputs[var1names[chan]] = np.array([(x_bins[chan][ix]+x_bins[chan][ix+1])/2. for ix in range(len(x_bins[chan])-1) for iy in range(len(y_bins[chan])-1)])
         inputs[var2names[chan]]  = np.array([(y_bins[chan][iy]+y_bins[chan][iy+1])/2. for ix in range(len(x_bins[chan])-1) for iy in range(len(y_bins[chan])-1)])
         inputs_up[var1names[chan]] = np.array([(x_bins[chan][ix]+x_bins[chan][ix+1])/2. for ix in range(len(x_bins[chan])-1) for iy in range(len(y_bins[chan])-1)])
         inputs_up[var2names[chan]]  = np.array([(y_bins[chan][iy]+y_bins[chan][iy+1])/2. for ix in range(len(x_bins[chan])-1) for iy in range(len(y_bins[chan])-1)])
         inputs_down[var1names[chan]] = np.array([(x_bins[chan][ix]+x_bins[chan][ix+1])/2. for ix in range(len(x_bins[chan])-1) for iy in range(len(y_bins[chan])-1)])
         inputs_down[var2names[chan]]  = np.array([(y_bins[chan][iy]+y_bins[chan][iy+1])/2. for ix in range(len(x_bins[chan])-1) for iy in range(len(y_bins[chan])-1)])
         inputs_mc[var1names[chan]] = np.array([(x_bins[chan][ix]+x_bins[chan][ix+1])/2. for ix in range(len(x_bins[chan])-1) for iy in range(len(y_bins[chan])-1)])
         inputs_mc[var2names[chan]]  = np.array([(y_bins[chan][iy]+y_bins[chan][iy+1])/2. for ix in range(len(x_bins[chan])-1) for iy in range(len(y_bins[chan])-1)])
         inputs_data[var1names[chan]] = np.array([(x_bins[chan][ix]+x_bins[chan][ix+1])/2. for ix in range(len(x_bins[chan])-1) for iy in range(len(y_bins[chan])-1)])
         inputs_data[var2names[chan]]  = np.array([(y_bins[chan][iy]+y_bins[chan][iy+1])/2. for ix in range(len(x_bins[chan])-1) for iy in range(len(y_bins[chan])-1)])
         inputs["weight"]  = np.divide(eff_hists_data[chan],eff_hists_mc[chan],out=np.ones_like(eff_hists_data[chan]), where=eff_hists_mc[chan]!=0.).flatten()
         inputs_up["weight"]  = np.divide(eff_hists_data_int[chan][1],eff_hists_mc_int[chan][0],out=np.ones_like(eff_hists_data_int[chan][1]), where=eff_hists_mc_int[chan][0]!=0.).flatten()
         inputs_down["weight"]  = np.divide(eff_hists_data_int[chan][0],eff_hists_mc_int[chan][1],out=np.ones_like(eff_hists_data_int[chan][0]), where=eff_hists_mc_int[chan][1]!=0.).flatten()
         arr_sf[chan] = inputs["weight"]
         arr_sf_up[chan] = inputs_up["weight"]
         arr_sf_down[chan] = inputs_down["weight"]
         inputs_mc["weight"]  = eff_hists_mc[chan].flatten()
         inputs_data["weight"]  = eff_hists_data[chan].flatten()
         
         h_trig_sf[args.year+"_trigsf_"+chan+"_nom"].fill(**inputs)
         h_trig_sf[args.year+"_trigsf_"+chan+"_up"].fill(**inputs_up)
         h_trig_sf[args.year+"_trigsf_"+chan+"_down"].fill(**inputs_down)
         h_trig_eff_mc[chan].fill(**inputs_mc)
         h_trig_eff_data[chan].fill(**inputs_data)

    for chan in chanlist:
        fig,ax = plt.subplots(1,1, figsize=(8,8))
        hist.plot2d(h_trig_sf[args.year+"_trigsf_"+chan+"_nom"],
                ax=ax,
                clear=True,
                xaxis=var1names[chan],
                )
        for i in range(len(y_bins[chan])-1):
            for j in range(len(x_bins[chan])-1):
                ax.text((x_bins[chan][j]+x_bins[chan][j+1])/2.,(y_bins[chan][i]+y_bins[chan][i+1])/2., "{:0.2f}".format(np.reshape(arr_sf[chan],(len(x_bins[chan])-1,len(y_bins[chan])-1))[j,i]) if eff_hists_mc[chan][j,i]>0. else "",
                    color="k", ha="center", va="center")#, fontweight="bold")
        
        fig.savefig("%s/trig_sf_debug_%s.pdf"%(args.tag,chan))

    for chan in chanlist:
        fig,ax = plt.subplots(1,1, figsize=(8,8))
        hist.plot2d(h_trig_sf[args.year+"_trigsf_"+chan+"_up"],
                ax=ax,
                clear=True,
                xaxis=var1names[chan],
                )
        for i in range(len(y_bins[chan])-1):
            for j in range(len(x_bins[chan])-1):
                ax.text((x_bins[chan][j]+x_bins[chan][j+1])/2.,(y_bins[chan][i]+y_bins[chan][i+1])/2., "{:0.2f}".format(np.reshape(arr_sf_up[chan]-arr_sf[chan],(len(x_bins[chan])-1,len(y_bins[chan])-1))[j,i]) if eff_hists_mc_int[chan][0][j,i]>0. else "",
                    color="k", ha="center", va="center")#, fontweight="bold")
        
        fig.savefig("%s/trig_sf_debug_up_%s.pdf"%(args.tag,chan))

    for chan in chanlist:
        fig,ax = plt.subplots(1,1, figsize=(8,8))
        hist.plot2d(h_trig_sf[args.year+"_trigsf_"+chan+"_down"],
                ax=ax,
                clear=True,
                xaxis=var1names[chan],
                )
        for i in range(len(y_bins[chan])-1):
            for j in range(len(x_bins[chan])-1):
                ax.text((x_bins[chan][j]+x_bins[chan][j+1])/2.,(y_bins[chan][i]+y_bins[chan][i+1])/2., "{:0.2f}".format(np.reshape(arr_sf_down[chan]-arr_sf[chan],(len(x_bins[chan])-1,len(y_bins[chan])-1))[j,i]) if eff_hists_mc_int[chan][1][j,i]>0. else "",
                    color="k", ha="center", va="center")#, fontweight="bold")
        
        fig.savefig("%s/trig_sf_debug_down_%s.pdf"%(args.tag,chan))

    for chan in chanlist:
        fig,ax = plt.subplots(1,1, figsize=(8,8))
        hist.plot2d(h_trig_eff_data[chan],
                ax=ax,
                clear=True,
                xaxis=var1names[chan],
                )
        for i in range(len(y_bins[chan])-1):
            for j in range(len(x_bins[chan])-1):
                ax.text((x_bins[chan][j]+x_bins[chan][j+1])/2.,(y_bins[chan][i]+y_bins[chan][i+1])/2., "{:0.2f}".format(eff_hists_data[chan][j,i]) if eff_hists_data[chan][j,i]>0. else "",
                    color="k", ha="center", va="center")#, fontweight="bold")
        
        fig.savefig("%s/trig_eff_data_debug_%s.pdf"%(args.tag,chan))

    for chan in chanlist:
        fig,ax = plt.subplots(1,1, figsize=(8,8))
        hist.plot2d(h_trig_eff_mc[chan],
                ax=ax,
                clear=True,
                xaxis=var1names[chan],
                )
        for i in range(len(y_bins[chan])-1):
            for j in range(len(x_bins[chan])-1):
                ax.text((x_bins[chan][j]+x_bins[chan][j+1])/2.,(y_bins[chan][i]+y_bins[chan][i+1])/2., "{:0.2f}".format(eff_hists_mc[chan][j,i]) if eff_hists_mc[chan][j,i]>0. else "",
                    color="k", ha="center", va="center")#, fontweight="bold")
        
        fig.savefig("%s/trig_eff_mc_debug_%s.pdf"%(args.tag,chan))
        
    print(h_trig_sf)
    save(h_trig_sf,"%s_%s.coffea"%(args.output,args.year))
