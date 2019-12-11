from __future__ import print_function, division
import os

import uproot
import matplotlib.pyplot as plt

import numpy as np
import scipy.interpolate
import argparse

from coffea import hist
from coffea.util import load
from sklearn.metrics import roc_curve
import processmap

def roc(h, bkg, sig, direction=-1):
    ax = h.axes()[-1]
    #print(h)
    bkgvals = h.project("process", bkg).values()[()]
    #print(bkgvals)
    sigvals = h.project("process", sig).values()[()]
    bkgeff_cut = np.cumsum(bkgvals[::direction])
    bkgeff_cut = bkgeff_cut[::direction] / bkgeff_cut[-1]
    sigeff_cut = np.cumsum(sigvals[::direction])
    sigeff_cut = sigeff_cut[::direction] / sigeff_cut[-1]
    interp = scipy.interpolate.interp1d(ax.centers(), np.c_[bkgeff_cut, sigeff_cut], axis=0, )
    #scores = np.concatenate((sigvals, bkgvals), axis=0)
    return ax.centers(), interp 

def getRoc(args):
    
    # vars
    vars = args.vars.split(',')

    # open hists
    hists_unmapped = load('%s.coffea'%args.hists)
    print(hists_unmapped)

    # map to hists
    hists_mapped = {}
    for key, val in hists_unmapped.items():
        if isinstance(val, hist.Hist):
            hists_mapped[key] = processmap.apply(val)
    print('hists mapped ',hists_mapped)

    # build roc for all vars
    vals = {}
    rocs = {}
    labels = {}

    for lep in ['ele','mu']:
        for jet in ['jet0','jet1']:
            vals['%s_%s'%(lep,jet)] = {}
            rocs['%s_%s'%(lep,jet)] = {}

    for var in vars:
        for lep in ['ele','mu']:
            for jet in ['jet0','jet1']:
                print('getting roc for ',var,lep,jet)
                hist_name = 'roc_%ssel%s'%(lep,jet)
                if 'lsf' in var:
                    var_name = jet+'_'+var
                    var_cut_dir = -1
                else:
                    var_name = lep+'0_'+var
                    var_cut_dir = 1
        
                # get hist
                h = hists_mapped[hist_name]
                print(h)
                print([ax for ax in h.axes() if ax.name not in {'process', var_name}])
                x = h.sum(*[ax for ax in h.axes() if ax.name not in {'process', var_name}])

                bkg = 'qcd'
                sig = 'h125'

                #vals['%s_%s'%(lep,jet)][var], rocs['%s_%s'%(lep,jet)][var] = roc(x, bkg, sig, direction=var_cut_dir)
                labels[var] = var

                # plot variable
                fig,ax = plt.subplots(1,1, figsize=(8,8))
                print(x)
                hist.plot1d(x,ax=ax,overlay='process',clear=False,density=True)
                fig.savefig("plots/rocs/lsf_%s_%s_%s.png"%(var,lep,jet))

    return vals,rocs,labels

def drawRoc(args,vals,rocs,labels,newtag="",newvar=""):
    # plot roc
    import matplotlib.patheffects
    fig, ax = plt.subplots(1,1)

    for key,roc in rocs.items():
        ax.plot(*roc(vals[key]).T, label=labels[key])

    ax.plot([0,1], [0,1], linestyle='--', color='grey')
    x_eps = np.linspace(0,1,100)
    y_eps = np.sqrt(x_eps)
    ax.plot(x_eps, y_eps, linestyle='-.', color='silver', label=r'$\epsilon_S = \sqrt{\epsilon_B}$')

    axtitle = r"Lepton in jet"
    tag = args.tag + newtag
    axtitle += "\n p$_T >$ 300 GeV"
    axtitle += "\n lep p$_T > 27$ GeV, $|\eta| < 2.4$"
    axtitle += "\n SIP < 4, $|dz|$ < 0.1, $|d0| < 0.05$"
    ax.legend(title=axtitle,loc='lower right')
    ax.autoscale(tight=True)
    bkg = "QCD"
    sig = "ggH"
    ax.set_ylabel("Signal efficiency (%s)"%sig)
    ax.set_xlabel("Background efficiency (%s)"%bkg)
    output=tag
    fig.savefig("plots/rocs/"+output)
    new_handles, new_labels = ax.get_legend_handles_labels()
    return new_handles, new_labels

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--hists',      dest='hists',    default="hists",     help="hists pickle name")
    parser.add_argument('--tag',        dest='tag',      default="",          help="tag")
    parser.add_argument('--vars',       dest='vars',     default="lsf3,miso,iso",  help="variables to compare performance"),
    parser.add_argument('--all',        dest='all', action='store_true', default=False, help="draw all rocs in one plot"),
    args = parser.parse_args()

    # get rocs
    vals,rocs,labels = getRoc(args)

    '''
    for key,val in vals.items():
        vals_i = val
        rocs_i = rocs[key]
        drawRoc(args,vals_i,rocs_i,labels,key,"")

        '''
