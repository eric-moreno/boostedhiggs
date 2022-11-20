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

import pickle
import gzip
import math
from cycler import cycler

import argparse
import processmap
from hists_map import *


import rhalphalib as rl
import scipy.stats
from scipy.special import erfinv
import pickle
#rl.util.install_roofit_helpers()
rl.ParametericSample.PreferRooParametricHist = True

import plot_stack
import samplelists

overflow_sum = 'allnan'

def expo_sample(norm, scale, obs):
    cdf = scipy.stats.expon.cdf(scale=scale, x=obs.binning) * norm
    return (np.diff(cdf), obs.binning, obs.name)


def gaus_sample(norm, loc, scale, obs):
    cdf = scipy.stats.norm.cdf(loc=loc, scale=scale, x=obs.binning) * norm
    return (np.diff(cdf), obs.binning, obs.name)

singleBinCR = True
includeLowMass = False

lowmassbin = 2
highmassbin = -1

fig, axs = plt.subplots(4, 2, sharex=True)
axs = axs.reshape(4,2)

fig_shape, axs_shape = plt.subplots(4, 2, sharex=True)
axs_shape = axs_shape.reshape(4,2)

def createLepHad(sig_hist, top_cr_hist, wlnu_cr_hist, qcd_cr_hist, sig_faillep, top_cr_faillep, wlnu_cr_faillep, qcd_cr_faillep, sig_faildphi, top_cr_faildphi, wlnu_cr_faildphi, qcd_cr_faildphi, var_name, mttbins, ptbins, leptype, tmpdir, label, usingData, nnCut, nnCut_loose, metCut, lowMetCut, h_pt_min, masspoints, shaperegion, year, doHtt, noSyst):

    #could be made region-dependent

    samp_combinations = {
        'data_obs'   : ['data'],
        'top'        : ['tt-had', 'tt-semilep', 'tt-dilep','st'],
        'htt125'     : ['h125'],
        'multijet'   : ['qcd'],
        'dy'         : ['zem','ztt'],
        'wlnu'       : ['wjets'],
        'vvqq'       : ['vv','vqq'],
        'ignore'     : [],
    }
    for m in masspoints:
        if doHtt:
            samp_combinations['ignore'].append('phi%s'%m)
        else:
            samp_combinations['phitt%s'%m] = ['phi%s'%m]
        
    sig_hists = {}
    top_cr_hists = {}
    wlnu_cr_hists = {}
    qcd_cr_hists = {}

    sig_hists["lowmet_faildphi"] = sig_faildphi.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,metCut),'under')
    top_cr_hists["lowmet_faildphi"] = top_cr_faildphi.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,metCut),'under')
    wlnu_cr_hists["lowmet_faildphi"] = wlnu_cr_faildphi.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,metCut),'under')
    qcd_cr_hists["lowmet_faildphi"] = qcd_cr_faildphi.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,metCut),'under')

    sig_hists["faildphi"] = sig_faildphi.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(metCut,None),'over')
    top_cr_hists["faildphi"] = top_cr_faildphi.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(metCut,None),'over')
    wlnu_cr_hists["faildphi"] = wlnu_cr_faildphi.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(metCut,None),'over')
    qcd_cr_hists["faildphi"] = qcd_cr_faildphi.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(metCut,None),'over')

    sig_hists["lowmet_faillep"] = sig_faillep.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,metCut),'under')
    top_cr_hists["lowmet_faillep"] = top_cr_faillep.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,metCut),'under')
    wlnu_cr_hists["lowmet_faillep"] = wlnu_cr_faillep.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,metCut),'under')
    qcd_cr_hists["lowmet_faillep"] = qcd_cr_faillep.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,metCut),'under')

    sig_hists["faillep"] = sig_faillep.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(metCut,None),'over')
    top_cr_hists["faillep"] = top_cr_faillep.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(metCut,None),'over')
    wlnu_cr_hists["faillep"] = wlnu_cr_faillep.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(metCut,None),'over')
    qcd_cr_hists["faillep"] = qcd_cr_faillep.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(metCut,None),'over')

    sig_hists["lowmet"] = sig_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,metCut),'under')
    top_cr_hists["lowmet"] = top_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,metCut),'under')
    wlnu_cr_hists["lowmet"] = wlnu_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,metCut),'under')
    qcd_cr_hists["lowmet"] = qcd_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,metCut),'under')

    sig_hists["nom"] = sig_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(metCut,None),'over')
    top_cr_hists["nom"] = top_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(metCut,None),'over')
    wlnu_cr_hists["nom"] = wlnu_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(metCut,None),'over')
    qcd_cr_hists["nom"] = qcd_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(metCut,None),'over')

    #jec = rl.NuisanceParameter('CMS_jec', 'lnN')
    top_norm = rl.NuisanceParameter('CMS_top_norm', 'lnN')
    wlnu_norm = rl.NuisanceParameter('CMS_wlnu_norm', 'lnN')
    qcd_norm = rl.NuisanceParameter('CMS_qcd_norm', 'lnN')
    dy_norm = rl.NuisanceParameter('CMS_dy_norm', 'lnN')
    vvqq_norm = rl.NuisanceParameter('CMS_vvqq_norm', 'lnN')
    trig = rl.NuisanceParameter('CMS_trig_had%s'%leptype, 'lnN')
    lumi = rl.NuisanceParameter('CMS_lumi', 'lnN')

    qcd_fail = rl.NuisanceParameter('qcd_Rfail_had%s'%leptype, 'shape')
    qcd_loosepass = rl.NuisanceParameter('qcd_Rloosepass_had%s'%leptype, 'shape')
    qcd_pass = rl.NuisanceParameter('qcd_Rpass_had%s'%leptype, 'shape')
    #qcdnormSF = rl.IndependentParameter('qcdnormSF_had%s'%leptype, 1., 0, 10)
    qcdnormSF = rl.NuisanceParameter('qcdnormSF_had%s'%leptype, 'lnN')
    #qcdtop_fail = rl.NuisanceParameter('qcd_top_Rfail_had%s'%leptype, 'shape')
    #qcdtop_loosepass = rl.NuisanceParameter('qcd_top_Rloosepass_had%s'%leptype, 'shape')
    #qcdtop_pass = rl.NuisanceParameter('qcd_top_Rpass_had%s'%leptype, 'shape')
    #qcdnormSF_top = rl.IndependentParameter('qcdnormSF_top_had%s'%leptype, 1., 0, 10)
    qcdnormSF_top = rl.NuisanceParameter('qcdnormSF_top_had%s'%leptype, 'lnN')
    #qcdwlnu_fail = rl.NuisanceParameter('qcd_wlnu_Rfail_had%s'%leptype, 'shape')
    #qcdwlnu_loosepass = rl.NuisanceParameter('qcd_wlnu_Rloosepass_had%s'%leptype, 'shape')
    #qcdwlnu_pass = rl.NuisanceParameter('qcd_wlnu_Rpass_had%s'%leptype, 'shape')
    #qcdnormSF_wlnu = rl.IndependentParameter('qcdnormSF_wlnu_had%s'%leptype, 1., 0, 10)
    qcdnormSF_wlnu = rl.NuisanceParameter('qcdnormSF_wlnu_had%s'%leptype, 'lnN')

    #m_scale = rl.NuisanceParameter('massscale_had%s'%leptype, 'shape')

    topeffSF = rl.IndependentParameter('topeffSF_had%s'%leptype, 1., 0, 10)
    #topnormSF = rl.IndependentParameter('topnormSF_had%s'%leptype, 1., 0, 10)
    wlnueffSF = rl.IndependentParameter('wlnueffSF_had%s'%leptype, 1., 0, 10)
    #wlnunormSF = rl.IndependentParameter('wlnunormSF_had%s'%leptype, 1., 0, 10)

    topLeffSF = rl.IndependentParameter('topLeffSF_had%s'%leptype, 1., 0, 10)
    wlnuLeffSF = rl.IndependentParameter('wlnuLeffSF_had%s'%leptype, 1., 0, 10)

    rdy = rl.IndependentParameter('r_dy_had%s'%leptype, 1., 0, 10)
    #rdy = rl.NuisanceParameter('r_dy_had%s'%leptype, 'lnN')
    rh125 = rl.NuisanceParameter('r_h125_had%s'%leptype, 'lnN')
    dy_eff = rl.IndependentParameter('dy_eff_had%s'%leptype, 1., 0, 10)

    toppt = rl.NuisanceParameter('toppt', 'shape')
    uescale = rl.NuisanceParameter('uescale', 'shape')
    jescale = rl.NuisanceParameter('jescale', 'shape')
    jeresol = rl.NuisanceParameter('jeresol', 'shape')
    l1prefire = rl.NuisanceParameter('l1prefire', 'shape')
    syst_dict = {
        samp:{
            'uescale':[uescale,'UESDown','UESUp'],
            'jescale':[jescale,'JESDown','JESUp'],
            'jeresol':[jeresol,'JERDown','JERUp'],
        } for samp in samp_combinations if samp not in ['data_obs','ignore']
    }
    syst_dict['top']['toppt'] = [toppt,'nominal','TopPtReweightUp']
    if year not in ['2018']:
        for samp in syst_dict:
            syst_dict[samp]['l1prefire'] = [l1prefire,'L1PreFiringDown','L1PreFiringUp']
    if noSyst:
        syst_dict = {}

    npt = len(ptbins) - 1
    ptslices = [slice(ptbins[ipt],ptbins[ipt+1]) for ipt in range(npt)]
    mtt = rl.Observable(var_name, mttbins[lowmassbin:highmassbin])
    mttone = rl.Observable(var_name+'_one', np.array([mttbins[1],mttbins[-2]]))

    qcdedges = [sig_hists["nom"].axis('h_pt').edges()[ix] for ix in range(len(sig_hists["nom"].axis('h_pt').edges()))]
    qcdbins = [[ix for ix in range(len(qcdedges)-1) if qcdedges[ix]>=(ptbins[ipt] if ptbins[ipt] is not None else -999999.) and qcdedges[ix+1]<=(ptbins[ipt+1] if ptbins[ipt+1] is not None else 999999.)] for ipt in range(npt)]
    print('qcdbins',qcdbins)
    
    model = rl.Model("had%sModel"%leptype)

    nnCut_lep = nnCut
    nnCut_lep_loose = nnCut_loose

    def intRegion(inhist,theregion,systematic='nominal',samplelist=None,debug=False,hptslice=slice(h_pt_min,None),mslice=None,mrebin=mttbins):
        if theregion=='pass':
            theslice = slice(nnCut_lep,None)
            overflow_str = 'over'
        elif theregion=='loosepass':
            theslice = slice(nnCut_lep_loose,nnCut_lep)
            overflow_str = 'none'
        elif theregion=='fail':
            theslice = slice(None,nnCut_lep_loose)
            overflow_str = 'under'
        else:
            print("Unknown region",theregion)
            return

        the_int = inhist.integrate('nn_disc',theslice,overflow_str).integrate('systematic',systematic)
        if debug:
            print('\t',nnCut_lep_loose,nnCut_lep)
            print('\t',overflow_str,the_int.values())

        if hptslice is not None:
            if hptslice.start is not None and hptslice.stop is not None:
                overflow_str = 'none'
            elif hptslice.start is None and hptslice.stop is not None:
                overflow_str = 'under'
            elif hptslice.start is not None and hptslice.stop is None:
                overflow_str = 'over'
            else:
                overflow_str = 'allnan'
            the_int = the_int.integrate('h_pt',hptslice,overflow_str)
            if debug:
                print('\t',hptslice)
                print('\t',overflow_str,the_int.values())
        if mrebin is not None:
            the_int = the_int.rebin('massreg',hist.Bin('massreg','massreg',mrebin))
        if mslice is not None:
            if mslice.start is not None and mslice.stop is not None:
                overflow_str = 'none'
            elif mslice.start is None and mslice.stop is not None:
                overflow_str = 'under'
            elif mslice.start is not None and mslice.stop is None:
                overflow_str = 'over'
            else:
                overflow_str = 'allnan'
            the_int = the_int.integrate('massreg',mslice,overflow_str)
            if debug:
                print('\t',overflow_str,the_int.values())

        if samplelist is not None:
            the_int = the_int.integrate('sample',samplelist).values(sumw2=True)[()]
        else:
            the_int = the_int.sum('sample').values(sumw2=True)
            if () in the_int:
                the_int = the_int[()]
            else:
                the_int = np.zeros(len(mttbins)-1,dtype=np.float32)

        if debug:
            print('\tdebug',the_int)

        return the_int

    def getQCDfromData(inhist,region,default=1.,systematic='nominal',hptslice=slice(h_pt_min,None),mslice=None):
        qcd_data_full = intRegion(inhist['data_obs'],region,systematic=systematic,hptslice=hptslice,mslice=mslice)[0][lowmassbin:highmassbin]
        other_mc = intRegion(inhist,region,samplelist=[s.name for s in inhist.identifiers('sample') if s.name!='data_obs' and s.name!='multijet'],systematic=systematic,hptslice=hptslice,mslice=mslice)
        qcd_data = qcd_data_full - other_mc[0][lowmassbin:highmassbin]
        #qcd_data = np.array([v if v>default else 0. for v in qcd_data])
        qcd_data = np.clip(qcd_data,default,None)
        qcd_data[qcd_data==default] = 0.
        qcd_data_w2 = np.clip(qcd_data_full-other_mc[1][lowmassbin:highmassbin],default,None)
        qcd_data_w2[qcd_data_w2==default] = 0.
        qcd_data_int = hist.poisson_interval(qcd_data,qcd_data_w2)
        qcd_temp = qcd_data
        qcd_temp_dn = np.array([qcd_data_int[0][bi] if qcd_data_int[0][bi] >= 0. and not np.isnan(qcd_data_int[0][bi]) else 0. for bi in range(len(qcd_data))])
        qcd_temp_up = np.array([qcd_data_int[1][bi] if qcd_data_int[1][bi] >= 0. and not np.isnan(qcd_data_int[1][bi]) else 0. for bi in range(len(qcd_data))])
        return qcd_temp,qcd_temp_dn,qcd_temp_up

    qcd_data_hists = {}
    top_data_hists = {}
    wlnu_data_hists = {}
    sig_data_hists = {}

    qcdfrom_sig = {}
    qcdfrom_top = {}
    qcdfrom_wlnu = {}
    qcdfrom_qcd = {}
    for pfregion in ["fail","loosepass","pass"]:
        qcdfrom_sig[pfregion] = {}
        qcdfrom_top[pfregion] = {}
        qcdfrom_wlnu[pfregion] = {}
        qcdfrom_qcd[pfregion] = {}

    def getQCDRatios(region, category="fail", verbose=False, defval=0.):
        qcd_data_hists[region] = [getQCDfromData(qcd_cr_hists[region],category,hptslice=slice(qcd_cr_hists[region].axis('h_pt').edges()[ix-1] if ix>1 else None,qcd_cr_hists[region].axis('h_pt').edges()[ix] if ix<len(qcd_cr_hists[region].axis('h_pt').edges())-1 else None),default=defval) for ix in range(len(qcd_cr_hists[region].axis('h_pt').edges()))]
        top_data_hists[region] = [getQCDfromData(top_cr_hists[region],category,hptslice=slice(top_cr_hists[region].axis('h_pt').edges()[ix-1] if ix>1 else None,top_cr_hists[region].axis('h_pt').edges()[ix] if ix<len(top_cr_hists[region].axis('h_pt').edges())-1 else None),default=defval) for ix in range(len(top_cr_hists[region].axis('h_pt').edges()))]
        wlnu_data_hists[region] = [getQCDfromData(wlnu_cr_hists[region],category,hptslice=slice(wlnu_cr_hists[region].axis('h_pt').edges()[ix-1] if ix>1 else None,wlnu_cr_hists[region].axis('h_pt').edges()[ix] if ix<len(wlnu_cr_hists[region].axis('h_pt').edges())-1 else None),default=defval) for ix in range(len(wlnu_cr_hists[region].axis('h_pt').edges()))]
        sig_data_hists[region] = [getQCDfromData(sig_hists[region],category,hptslice=slice(sig_hists[region].axis('h_pt').edges()[ix-1] if ix>1 else None,sig_hists[region].axis('h_pt').edges()[ix] if ix<len(sig_hists[region].axis('h_pt').edges())-1 else None),default=defval) for ix in range(len(sig_hists[region].axis('h_pt').edges()))]
        # 0 = nom, 1 = dn, 2 = up
    
        #qcdratio_sig[region] = [[np.sum(np.stack([(sig_data_hists[region][ix][t]) for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else sig_data_hists[region][t]),axis=0) for ptbin in range(npt)] for t in range(3)]
        #qcdratio_top[region] = [[np.sum(np.stack([(top_data_hists[region][ix][t]) for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else top_data_hists[region][t]),axis=0) for ptbin in range(npt)] for t in range(3)]
        #qcdratio_wlnu[region] = [[np.sum(np.stack([(wlnu_data_hists[region][ix][t]) for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else wlnu_data_hists[region][t]),axis=0) for ptbin in range(npt)] for t in range(3)]
        #qcdratio_qcd[region] = [[np.sum(np.stack([(qcd_data_hists[region][ix][t]) for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else qcd_data_hists[region][t]),axis=0) for ptbin in range(npt)] for t in range(3)] #does this work for len(qcdbins[ptbin])=0?

        qcdfrom_sig[category][region] = [[np.sum(np.stack([(sig_data_hists[region][ix][t]) for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else sig_data_hists[region][t]),axis=0) for ptbin in range(npt)] for t in range(3)]
        qcdfrom_top[category][region] = [[np.sum(np.stack([(top_data_hists[region][ix][t]) for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else top_data_hists[region][t]),axis=0) for ptbin in range(npt)] for t in range(3)]
        qcdfrom_wlnu[category][region] = [[np.sum(np.stack([(wlnu_data_hists[region][ix][t]) for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else wlnu_data_hists[region][t]),axis=0) for ptbin in range(npt)] for t in range(3)]
        qcdfrom_qcd[category][region] = [[np.sum(np.stack([(qcd_data_hists[region][ix][t]) for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else qcd_data_hists[region][t]),axis=0) for ptbin in range(npt)] for t in range(3)] #does this work for len(qcdbins[ptbin])=0?

        if verbose:
            print(region,category,'ratio sig ',qcdfrom_sig[category][region])
            print(region,category,'ratio top ',qcdfrom_top[category][region])
            print(region,category,'ratio wlnu',qcdfrom_wlnu[category][region])
            print(region,category,'ratio qcd ',qcdfrom_qcd[category][region])
        
    for qcdregion in ["faildphi","lowmet_faildphi","faillep","lowmet_faillep","lowmet","nom"]:
        for pfregion in ["fail","loosepass","pass"]:
            getQCDRatios(qcdregion,category=pfregion)

    #qcdfrom_sig_err = [np.abs(qcdfrom_sig["lowmet"][ix]-qcdfrom_sig["faildphi"][ix])/2. for ix in range(len(qcdfrom_sig["lowmet"]))]
    #qcdfrom_sig = [(qcdfrom_sig["lowmet"][ix]+qcdfrom_sig["faildphi"][ix])/2. for ix in range(len(qcdfrom_sig["lowmet"]))]
    #qcdfrom_top_err = [np.abs(qcdfrom_top["lowmet"][ix]-qcdfrom_top["faildphi"][ix])/2. for ix in range(len(qcdfrom_top["lowmet"]))]
    #qcdfrom_top = [(qcdfrom_top["lowmet"][ix]+qcdfrom_top["faildphi"][ix])/2. for ix in range(len(qcdfrom_top["lowmet"]))]
    #qcdfrom_wlnu_err = [np.abs(qcdfrom_wlnu["lowmet"][ix]-qcdfrom_wlnu["faildphi"][ix])/2. for ix in range(len(qcdfrom_wlnu["lowmet"]))]
    #qcdfrom_wlnu = [(qcdfrom_wlnu["lowmet"][ix]+qcdfrom_wlnu["faildphi"][ix])/2. for ix in range(len(qcdfrom_wlnu["lowmet"]))]

    #print('qcdfrom_sig',qcdfrom_sig, qcdfrom_sig_err)
    #print('qcdfrom_top',qcdfrom_top, qcdfrom_top_err)
    #print('qcdfrom_wlnu',qcdfrom_wlnu, qcdfrom_wlnu_err)

    qcdfail_temp = {}
    qcdfail_temp_dn = {}
    qcdfail_temp_up = {}

    for region in ['pass','loosepass','fail']:
        theqcdregion = "lowmet"
        failtemp = [getQCDfromData(sig_hists[theqcdregion],region,hptslice=slice(sig_hists[theqcdregion].axis('h_pt').edges()[ix-1] if ix>1 else None,sig_hists[theqcdregion].axis('h_pt').edges()[ix] if ix<len(sig_hists[theqcdregion].axis('h_pt').edges())-1 else None),default=0.) for ix in range(len(sig_hists[theqcdregion].axis('h_pt').edges()))]
        qcdfail_temp[region] = [[np.sum(np.stack([(failtemp[ix][t]) for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else failtemp[t]),axis=0) for ptbin in range(npt)] for t in range(3)]
        qcdfail_temp_dn[region] = qcdfail_temp[region][:][1]
        qcdfail_temp_up[region] = qcdfail_temp[region][:][2]
        qcdfail_temp[region] = qcdfail_temp[region][:][0]

    qcdratio_F = {
        'pass':[qcdfail_temp['pass'][ipt]/qcdfail_temp['fail'][ipt] for ipt in range(npt)],
        'loosepass':[qcdfail_temp['loosepass'][ipt]/qcdfail_temp['fail'][ipt] for ipt in range(npt)],
        'fail':[np.ones_like(qcdfail_temp['fail'][ipt]) for ipt in range(npt)],
    }
    qcdratio_F_dn = {
        'pass':[qcdfail_temp_dn['pass'][ipt]/qcdfail_temp_dn['fail'][ipt] for ipt in range(npt)],
        'loosepass':[qcdfail_temp_dn['loosepass'][ipt]/qcdfail_temp_dn['fail'][ipt] for ipt in range(npt)],
        'fail':[np.ones_like(qcdfail_temp['fail'][ipt]) for ipt in range(npt)],
    }
    qcdratio_F_up = {
        'pass':[qcdfail_temp_up['pass'][ipt]/qcdfail_temp_up['fail'][ipt] for ipt in range(npt)],
        'loosepass':[qcdfail_temp_up['loosepass'][ipt]/qcdfail_temp_up['fail'][ipt] for ipt in range(npt)],
        'fail':[np.ones_like(qcdfail_temp['fail'][ipt]) for ipt in range(npt)],
    }
    #for region in ['pass', 'loosepass', 'fail']:
    #    qcdratio_F[region] = np.sum(qcdratio_F[region]*qcd_from_data)/np.sum(qcd_from_data)
    #    qcdratio_F_dn[region] = np.sum(qcdratio_F_dn[region]*qcd_from_data)/np.sum(qcd_from_data)
    #    qcdratio_F_up[region] = np.sum(qcdratio_F_up[region]*qcd_from_data)/np.sum(qcd_from_data)
    print('qcdratio_F',qcdratio_F)
    print('qcdratio_F_dn',qcdratio_F_dn)
    print('qcdratio_F_up',qcdratio_F_up)

    mbinend=-4
    for iregion,region in enumerate(["loosepass","pass"]):
        for ipt in range(npt):
            axs_shape[3,int(leptype=='el')].step((mttbins[lowmassbin:highmassbin-1]+mttbins[lowmassbin+1:highmassbin])[:mbinend]/2.,qcdratio_F[region][ipt][:mbinend], color='C%i'%(iregion*npt+ipt),label='%s, %i<pT<%i'%(region, ptbins[ipt] if ptbins[ipt] is not None else 0, ptbins[ipt+1] if ptbins[ipt+1] is not None else 9999),where="mid")
            axs_shape[3,int(leptype=='el')].fill_between((mttbins[lowmassbin:highmassbin-1]+mttbins[lowmassbin+1:highmassbin])[:mbinend]/2., qcdratio_F_dn[region][ipt][:mbinend], qcdratio_F_up[region][ipt][:mbinend], alpha=0.5, color='C%i'%(iregion*npt+ipt),step="mid")

    qcdfail_temp = {}
    qcdfail_temp_dn = {}
    qcdfail_temp_up = {}

    for qcdregion in ["lowmet","qcdlowmet","qcdnom"]:
        qcdfail_temp[qcdregion] = {}
        qcdfail_temp_dn[qcdregion] = {}
        qcdfail_temp_up[qcdregion] = {}
        for region in ['pass','loosepass','fail']:
            #qcdfail_temp[region], qcdfail_temp_dn[region], qcdfail_temp_up[region] = getQCDfromData(qcd_cr_hists["lowmet"],region,default=0.)
            qcdfail_temp[qcdregion][region], qcdfail_temp_dn[qcdregion][region], qcdfail_temp_up[qcdregion][region] = getQCDfromData(sig_hists[qcdregion] if "qcd" not in qcdregion else qcd_cr_hists[qcdregion.replace('qcd','')],region,default=0.)
            #qcdfail_temp[region] = np.array([np.sum(qcdfail_temp[region])])
            #qcdfail_temp_dn[region] = np.array([np.sum(qcdfail_temp_dn[region])])
            #qcdfail_temp_up[region] = np.array([np.sum(qcdfail_temp_up[region])])

    qcdratio_F = {}
    qcdratio_F_dn = {}
    qcdratio_F_up = {}

    for qcdregion in ["lowmet","qcdlowmet","qcdnom"]:
        qcdratio_F[qcdregion] = {}
        qcdratio_F_dn[qcdregion] = {}
        qcdratio_F_up[qcdregion] = {}
        for ireg,reg in enumerate(["fail","loosepass","pass"]):
            #qcdratio_F[qcdregion][reg] = np.ones_like(qcdfail_temp[qcdregion][reg]) if reg in shaperegion else np.divide(qcdfail_temp[qcdregion][reg],qcdfail_temp[qcdregion][shaperegion[ireg]], out=qcdfail_temp[qcdregion][reg], where=qcdfail_temp[qcdregion][shaperegion[ireg]]!=0.)
            #qcdratio_F_dn[qcdregion][reg] = np.ones_like(qcdfail_temp_dn[qcdregion][reg]) if reg in shaperegion else np.divide(qcdfail_temp_dn[qcdregion][reg],qcdfail_temp_dn[qcdregion][shaperegion[ireg]], out=qcdfail_temp_dn[qcdregion][reg], where=qcdfail_temp_dn[qcdregion][shaperegion[ireg]]!=0.)
            #qcdratio_F_up[qcdregion][reg] = np.ones_like(qcdfail_temp_up[qcdregion][reg]) if reg in shaperegion else np.divide(qcdfail_temp_up[qcdregion][reg],qcdfail_temp_up[qcdregion][shaperegion[ireg]], out=qcdfail_temp_up[qcdregion][reg], where=qcdfail_temp_up[qcdregion][shaperegion[ireg]]!=0.)
            qcdratio_F[qcdregion][reg] = np.ones_like(qcdfail_temp[qcdregion][reg])*(1. if reg in shaperegion else np.sum(qcdfail_temp[qcdregion][reg])/np.sum(qcdfail_temp[qcdregion][shaperegion[ireg]]))
            qcdratio_F_dn[qcdregion][reg] = np.ones_like(qcdfail_temp_dn[qcdregion][reg])*(1. if reg in shaperegion else np.sum(qcdfail_temp_dn[qcdregion][reg])/np.sum(qcdfail_temp_dn[qcdregion][shaperegion[ireg]]))
            qcdratio_F_up[qcdregion][reg] = np.ones_like(qcdfail_temp_up[qcdregion][reg])*(1. if reg in shaperegion else np.sum(qcdfail_temp_up[qcdregion][reg])/np.sum(qcdfail_temp_up[qcdregion][shaperegion[ireg]]))
    for qcdreg in qcdratio_F:
        for reg in qcdratio_F[qcdreg]:
            qcdratio_F[qcdreg][reg] = np.nan_to_num(qcdratio_F[qcdreg][reg])
            qcdratio_F_dn[qcdreg][reg] = np.nan_to_num(qcdratio_F_dn[qcdreg][reg])
            qcdratio_F_up[qcdreg][reg] = np.nan_to_num(qcdratio_F_up[qcdreg][reg])

    qcd_from_data_sig = {reg:{pf:[getQCDfromData(sig_hists[reg],shaperegion[ipf],hptslice=slice(sig_hists[reg].axis('h_pt').edges()[ix-1] if ix>1 else None,sig_hists[reg].axis('h_pt').edges()[ix] if ix<len(sig_hists[reg].axis('h_pt').edges())-1 else None)) for ix in range(len(sig_hists[reg].axis('h_pt').edges()))] for ipf,pf in enumerate(['fail', 'loosepass', 'pass'])} for reg in ["nom","faildphi","lowmet_faildphi","faillep","lowmet_faillep","lowmet"]}
    qcd_from_data_qcd = {reg:{pf:[getQCDfromData(qcd_cr_hists[reg],shaperegion[ipf],hptslice=slice(qcd_cr_hists[reg].axis('h_pt').edges()[ix-1] if ix>1 else None,qcd_cr_hists[reg].axis('h_pt').edges()[ix] if ix<len(qcd_cr_hists[reg].axis('h_pt').edges())-1 else None)) for ix in range(len(qcd_cr_hists[reg].axis('h_pt').edges()))] for ipf,pf in enumerate(['fail', 'loosepass', 'pass'])} for reg in ["nom","faildphi","lowmet_faildphi","faillep","lowmet_faillep","lowmet"]}
    qcd_from_data_top = {reg:{pf:[getQCDfromData(top_cr_hists[reg],shaperegion[ipf],hptslice=slice(top_cr_hists[reg].axis('h_pt').edges()[ix-1] if ix>1 else None,top_cr_hists[reg].axis('h_pt').edges()[ix] if ix<len(top_cr_hists[reg].axis('h_pt').edges())-1 else None)) for ix in range(len(top_cr_hists[reg].axis('h_pt').edges()))] for ipf,pf in enumerate(['fail', 'loosepass', 'pass'])} for reg in ["nom","faildphi","lowmet_faildphi","faillep","lowmet_faillep","lowmet"]}
    qcd_from_data_wlnu = {reg:{pf:[getQCDfromData(wlnu_cr_hists[reg],shaperegion[ipf],hptslice=slice(wlnu_cr_hists[reg].axis('h_pt').edges()[ix-1] if ix>1 else None,wlnu_cr_hists[reg].axis('h_pt').edges()[ix] if ix<len(wlnu_cr_hists[reg].axis('h_pt').edges())-1 else None)) for ix in range(len(wlnu_cr_hists[reg].axis('h_pt').edges()))] for ipf,pf in enumerate(['fail', 'loosepass', 'pass'])} for reg in ["nom","faildphi","lowmet_faildphi","faillep","lowmet_faillep","lowmet"]}
    #print('qcd_from_data_sig',qcd_from_data_sig)
    #print('qcd_from_data_top',qcd_from_data_top)
    #print('qcd_from_data_wlnu',qcd_from_data_wlnu)
    
#    qcdratio_sig_dn = [((qcdfrom_sig["lowmet"][ptbin])/(qcdfrom_qcd["lowmet"][ptbin])) for ptbin in range(npt)]
#    qcdratio_sig_up = [((qcdfrom_sig["faildphi"][ptbin])/(qcdfrom_qcd["faildphi"][ptbin])) for ptbin in range(npt)]
#    qcdratio_sig = [((qcdfrom_sig["lowmet"][ptbin]/qcdfrom_qcd["lowmet"][ptbin]) + (qcdfrom_sig["faildphi"][ptbin]/qcdfrom_qcd["faildphi"][ptbin]))/2. for ptbin in range(npt)]
#
#    qcdratio_top_dn = [((qcdfrom_top["lowmet"][ptbin])/(qcdfrom_qcd["lowmet"][ptbin])) for ptbin in range(npt)]
#    qcdratio_top_up = [((qcdfrom_top["faildphi"][ptbin])/(qcdfrom_qcd["faildphi"][ptbin])) for ptbin in range(npt)]
#    qcdratio_top = [((qcdfrom_top["lowmet"][ptbin]/qcdfrom_qcd["lowmet"][ptbin]) + (qcdfrom_top["faildphi"][ptbin]/qcdfrom_qcd["faildphi"][ptbin]))/2. for ptbin in range(npt)]
#
#    qcdratio_wlnu_dn = [((qcdfrom_wlnu["lowmet"][ptbin])/(qcdfrom_qcd["lowmet"][ptbin])) for ptbin in range(npt)]
#    qcdratio_wlnu_up = [((qcdfrom_wlnu["faildphi"][ptbin])/(qcdfrom_qcd["faildphi"][ptbin])) for ptbin in range(npt)]
#    qcdratio_wlnu = [((qcdfrom_wlnu["lowmet"][ptbin]/qcdfrom_qcd["lowmet"][ptbin]) + (qcdfrom_wlnu["faildphi"][ptbin]/qcdfrom_qcd["faildphi"][ptbin]))/2. for ptbin in range(npt)]

    qcdratio_sig_dn = {}
    qcdratio_sig_up = {}
    qcdratio_sig = {}

    qcdratio_qcd_dn = {}
    qcdratio_qcd_up = {}
    qcdratio_qcd = {}

    qcdratio_top_dn = {}
    qcdratio_top_up = {}
    qcdratio_top = {}

    qcdratio_wlnu_dn = {}
    qcdratio_wlnu_up = {}
    qcdratio_wlnu = {}

    for pfreg in ["fail","loosepass","pass"]:
        qcdratio_sig_dn[pfreg] = {}
        qcdratio_sig_up[pfreg] = {}
        qcdratio_sig[pfreg] = {}
        qcdratio_qcd_dn[pfreg] = {}
        qcdratio_qcd_up[pfreg] = {}
        qcdratio_qcd[pfreg] = {}
        qcdratio_top_dn[pfreg] = {}
        qcdratio_top_up[pfreg] = {}
        qcdratio_top[pfreg] = {}
        qcdratio_wlnu_dn[pfreg] = {}
        qcdratio_wlnu_up[pfreg] = {}
        qcdratio_wlnu[pfreg] = {}
        for qcdreg in ["faildphi","lowmet_faildphi","faillep","lowmet_faillep","lowmet","nom"]:
            qcdratio_sig_dn[pfreg][qcdreg] = []
            qcdratio_sig_up[pfreg][qcdreg] = []
            qcdratio_sig[pfreg][qcdreg] = []
            qcdratio_qcd_dn[pfreg][qcdreg] = []
            qcdratio_qcd_up[pfreg][qcdreg] = []
            qcdratio_qcd[pfreg][qcdreg] = []
            qcdratio_top_dn[pfreg][qcdreg] = []
            qcdratio_top_up[pfreg][qcdreg] = []
            qcdratio_top[pfreg][qcdreg] = []
            qcdratio_wlnu_dn[pfreg][qcdreg] = []
            qcdratio_wlnu_up[pfreg][qcdreg] = []
            qcdratio_wlnu[pfreg][qcdreg] = []
            for ptbin in range(npt):
                marr_dn = []
                marr_up = []
                marr = []
                for ix in range(len(qcdfrom_qcd[pfreg][qcdreg][0][ptbin])):
                    if qcdfrom_qcd[pfreg][qcdreg][0][ptbin][ix]>0.:
                        marr_dn.append(qcdfrom_qcd[pfreg]["nom"][1][ptbin][ix]/qcdfrom_qcd[pfreg][qcdreg][2][ptbin][ix])
                        marr_up.append(qcdfrom_qcd[pfreg]["nom"][2][ptbin][ix]/qcdfrom_qcd[pfreg][qcdreg][1][ptbin][ix])
                        marr.append(qcdfrom_qcd[pfreg]["nom"][0][ptbin][ix]/qcdfrom_qcd[pfreg][qcdreg][0][ptbin][ix])
                    else:
                        if ix==0:
                            marr_dn.append(1.)
                            marr_up.append(1.)
                            marr.append(1.)
                        else:
                            marr_dn.append(marr_dn[-1])
                            marr_up.append(marr_up[-1])
                            marr.append(marr[-1])
                qcdratio_sig_dn[pfreg][qcdreg].append(np.array(marr_dn))
                qcdratio_sig_up[pfreg][qcdreg].append(np.array(marr_up))
                qcdratio_sig[pfreg][qcdreg].append(np.array(marr))
    
                marr_dn = []
                marr_up = []
                marr = []
                for ix in range(len(qcdfrom_sig[pfreg][qcdreg][0][ptbin])):
                    if qcdfrom_qcd[pfreg][qcdreg][0][ptbin][ix]>0.:
                        marr_dn.append(qcdfrom_sig[pfreg][qcdreg][1][ptbin][ix]/qcdfrom_qcd[pfreg][qcdreg][2][ptbin][ix])
                        marr_up.append(qcdfrom_sig[pfreg][qcdreg][2][ptbin][ix]/qcdfrom_qcd[pfreg][qcdreg][1][ptbin][ix])
                        marr.append(qcdfrom_sig[pfreg][qcdreg][0][ptbin][ix]/qcdfrom_qcd[pfreg][qcdreg][0][ptbin][ix])
                    else:
                        if ix==0:
                            marr_dn.append(1.)
                            marr_up.append(1.)
                            marr.append(1.)
                        else:
                            marr_dn.append(marr_dn[-1])
                            marr_up.append(marr_up[-1])
                            marr.append(marr[-1])
                qcdratio_qcd_dn[pfreg][qcdreg].append(np.array(marr_dn))
                qcdratio_qcd_up[pfreg][qcdreg].append(np.array(marr_up))
                qcdratio_qcd[pfreg][qcdreg].append(np.array(marr))
    
                marr_dn = []
                marr_up = []
                marr = []
                for ix in range(len(qcdfrom_top[pfreg][qcdreg][0][ptbin])):
                    if qcdfrom_qcd[pfreg][qcdreg][0][ptbin][ix]>0.:
                        marr_dn.append(qcdfrom_top[pfreg][qcdreg][1][ptbin][ix]/qcdfrom_qcd[pfreg][qcdreg][2][ptbin][ix])
                        marr_up.append(qcdfrom_top[pfreg][qcdreg][2][ptbin][ix]/qcdfrom_qcd[pfreg][qcdreg][1][ptbin][ix])
                        marr.append(qcdfrom_top[pfreg][qcdreg][0][ptbin][ix]/qcdfrom_qcd[pfreg][qcdreg][0][ptbin][ix])
                    else:
                        if ix==0:
                            marr_dn.append(1.)
                            marr_up.append(1.)
                            marr.append(1.)
                        else:
                            marr_dn.append(marr_dn[-1])
                            marr_up.append(marr_up[-1])
                            marr.append(marr[-1])
                qcdratio_top_dn[pfreg][qcdreg].append(np.array(marr_dn))
                qcdratio_top_up[pfreg][qcdreg].append(np.array(marr_up))
                qcdratio_top[pfreg][qcdreg].append(np.array(marr))
    
                marr_dn = []
                marr_up = []
                marr = []
                for ix in range(len(qcdfrom_wlnu[pfreg][qcdreg][0][ptbin])):
                    if qcdfrom_qcd[pfreg][qcdreg][0][ptbin][ix]>0.:
                        marr_dn.append(qcdfrom_wlnu[pfreg][qcdreg][1][ptbin][ix]/qcdfrom_qcd[pfreg][qcdreg][2][ptbin][ix])
                        marr_up.append(qcdfrom_wlnu[pfreg][qcdreg][2][ptbin][ix]/qcdfrom_qcd[pfreg][qcdreg][1][ptbin][ix])
                        marr.append(qcdfrom_wlnu[pfreg][qcdreg][0][ptbin][ix]/qcdfrom_qcd[pfreg][qcdreg][0][ptbin][ix])
                    else:
                        if ix==0:
                            marr_dn.append(1.)
                            marr_up.append(1.)
                            marr.append(1.)
                        else:
                            marr_dn.append(marr_dn[-1])
                            marr_up.append(marr_up[-1])
                            marr.append(marr[-1])
                qcdratio_wlnu_dn[pfreg][qcdreg].append(np.array(marr_dn))
                qcdratio_wlnu_up[pfreg][qcdreg].append(np.array(marr_up))
                qcdratio_wlnu[pfreg][qcdreg].append(np.array(marr))
                

    #qcdratio_sig_up = [np.divide((qcdratio_qcd["nom"][2][ptbin]),(qcdratio_qcd["lowmet"][0][ptbin]), out=np.zeros_like(qcdratio_qcd["nom"][2][ptbin]), where=qcdratio_qcd["lowmet"][0][ptbin]!=0.) for ptbin in range(npt)]
    #qcdratio_sig_dn = [np.divide((qcdratio_qcd["nom"][1][ptbin]),(qcdratio_qcd["lowmet"][0][ptbin]), out=np.zeros_like(qcdratio_qcd["nom"][2][ptbin]), where=qcdratio_qcd["lowmet"][0][ptbin]!=0.) for ptbin in range(npt)]
    #qcdratio_sig = [np.divide((qcdratio_qcd["nom"][0][ptbin]),(qcdratio_qcd["lowmet"][0][ptbin]), out=np.zeros_like(qcdratio_qcd["nom"][2][ptbin]), where=qcdratio_qcd["lowmet"][0][ptbin]!=0.) for ptbin in range(npt)]
    
    for pfreg in ["fail","loosepass","pass"]:
        #for iq,qcdreg in enumerate(["faildphi","faillep","lowmet"]):
        for iq,qcdreg in enumerate(["lowmet","faillep","nom"]):
            for ptbin in range(npt):
                print(ptbins[ptbin],ptbins[ptbin+1])
                axs[0,int(leptype=='el')].step((mttbins[lowmassbin:highmassbin-1]+mttbins[lowmassbin+1:highmassbin])[:mbinend]/2.,qcdratio_sig[pfreg][qcdreg][ptbin][:mbinend], color='C%i'%(iq+2*ptbin),label='%s, %i<pT<%i'%(qcdreg,ptbins[ptbin] if not ptbins[ptbin] is None else 0,ptbins[ptbin+1] if not ptbins[ptbin+1] is None else 9999),where="mid")
                axs[0,int(leptype=='el')].fill_between((mttbins[lowmassbin:highmassbin-1]+mttbins[lowmassbin+1:highmassbin])[:mbinend]/2., qcdratio_sig_dn[pfreg][qcdreg][ptbin][:mbinend], qcdratio_sig_up[pfreg][qcdreg][ptbin][:mbinend], alpha=0.5, color='C%i'%(iq+2*ptbin),step="mid")
                axs[1,int(leptype=='el')].step((mttbins[lowmassbin:highmassbin-1]+mttbins[lowmassbin+1:highmassbin])[:mbinend]/2.,qcdratio_qcd[pfreg][qcdreg][ptbin][:mbinend], color='C%i'%(iq+2*ptbin),label='%s, %i<pT<%i'%(qcdreg,ptbins[ptbin] if not ptbins[ptbin] is None else 0,ptbins[ptbin+1] if not ptbins[ptbin+1] is None else 9999),where="mid")
                axs[1,int(leptype=='el')].fill_between((mttbins[lowmassbin:highmassbin-1]+mttbins[lowmassbin+1:highmassbin])[:mbinend]/2., qcdratio_qcd_dn[pfreg][qcdreg][ptbin][:mbinend], qcdratio_qcd_up[pfreg][qcdreg][ptbin][:mbinend], alpha=0.5, color='C%i'%(iq+2*ptbin),step="mid")
                axs[2,int(leptype=='el')].step((mttbins[lowmassbin:highmassbin-1]+mttbins[lowmassbin+1:highmassbin])[:mbinend]/2.,qcdratio_top[pfreg][qcdreg][ptbin][:mbinend], color='C%i'%(iq+2*ptbin),label='%s, %i<pT<%i'%(qcdreg,ptbins[ptbin] if not ptbins[ptbin] is None else 0,ptbins[ptbin+1] if not ptbins[ptbin+1] is None else 9999),where="mid")
                axs[2,int(leptype=='el')].fill_between((mttbins[lowmassbin:highmassbin-1]+mttbins[lowmassbin+1:highmassbin])[:mbinend]/2., qcdratio_top_dn[pfreg][qcdreg][ptbin][:mbinend], qcdratio_top_up[pfreg][qcdreg][ptbin][:mbinend], alpha=0.5, color='C%i'%(iq+2*ptbin),step="mid")
                axs[3,int(leptype=='el')].step((mttbins[lowmassbin:highmassbin-1]+mttbins[lowmassbin+1:highmassbin])[:mbinend]/2.,qcdratio_wlnu[pfreg][qcdreg][ptbin][:mbinend], color='C%i'%(iq+2*ptbin),label='%s, %i<pT<%i'%(qcdreg,ptbins[ptbin] if not ptbins[ptbin] is None else 0,ptbins[ptbin+1] if not ptbins[ptbin+1] is None else 9999),where="mid")
                axs[3,int(leptype=='el')].fill_between((mttbins[lowmassbin:highmassbin-1]+mttbins[lowmassbin+1:highmassbin])[:mbinend]/2., qcdratio_wlnu_dn[pfreg][qcdreg][ptbin][:mbinend], qcdratio_qcd_up[pfreg][qcdreg][ptbin][:mbinend], alpha=0.5, color='C%i'%(iq+2*ptbin),step="mid")

    for pfreg in ["fail","loosepass","pass"]:
        for qcdreg in ["faildphi","lowmet_faildphi","faillep","lowmet_faillep","lowmet","nom"]:
            for ptbin in range(npt):
                qcdratio_sig[pfreg][qcdreg][ptbin] = np.sum(qcdratio_sig[pfreg][qcdreg][ptbin]*qcdfrom_sig[pfreg][qcdreg][0][ptbin])/np.sum(qcdfrom_sig[pfreg][qcdreg][0][ptbin])
                qcdratio_qcd[pfreg][qcdreg][ptbin] = np.sum(qcdratio_qcd[pfreg][qcdreg][ptbin]*qcdfrom_qcd[pfreg][qcdreg][0][ptbin])/np.sum(qcdfrom_qcd[pfreg][qcdreg][0][ptbin])
                qcdratio_top[pfreg][qcdreg][ptbin] = np.sum(qcdratio_top[pfreg][qcdreg][ptbin]*qcdfrom_top[pfreg][qcdreg][0][ptbin])/np.sum(qcdfrom_top[pfreg][qcdreg][0][ptbin])
                qcdratio_wlnu[pfreg][qcdreg][ptbin] = np.sum(qcdratio_wlnu[pfreg][qcdreg][ptbin]*qcdfrom_wlnu[pfreg][qcdreg][0][ptbin])/np.sum(qcdfrom_wlnu[pfreg][qcdreg][0][ptbin])
                qcdratio_sig_dn[pfreg][qcdreg][ptbin] = np.sum(qcdratio_sig_dn[pfreg][qcdreg][ptbin]*qcdfrom_sig[pfreg][qcdreg][1][ptbin])/np.sum(qcdfrom_sig[pfreg][qcdreg][1][ptbin])
                qcdratio_qcd_dn[pfreg][qcdreg][ptbin] = np.sum(qcdratio_qcd_dn[pfreg][qcdreg][ptbin]*qcdfrom_qcd[pfreg][qcdreg][1][ptbin])/np.sum(qcdfrom_qcd[pfreg][qcdreg][1][ptbin])
                qcdratio_top_dn[pfreg][qcdreg][ptbin] = np.sum(qcdratio_top_dn[pfreg][qcdreg][ptbin]*qcdfrom_top[pfreg][qcdreg][1][ptbin])/np.sum(qcdfrom_top[pfreg][qcdreg][1][ptbin])
                qcdratio_wlnu_dn[pfreg][qcdreg][ptbin] = np.sum(qcdratio_wlnu_dn[pfreg][qcdreg][ptbin]*qcdfrom_wlnu[pfreg][qcdreg][1][ptbin])/np.sum(qcdfrom_wlnu[pfreg][qcdreg][1][ptbin])
                qcdratio_sig_up[pfreg][qcdreg][ptbin] = np.sum(qcdratio_sig_up[pfreg][qcdreg][ptbin]*qcdfrom_sig[pfreg][qcdreg][2][ptbin])/np.sum(qcdfrom_sig[pfreg][qcdreg][2][ptbin])
                qcdratio_qcd_up[pfreg][qcdreg][ptbin] = np.sum(qcdratio_qcd_up[pfreg][qcdreg][ptbin]*qcdfrom_qcd[pfreg][qcdreg][2][ptbin])/np.sum(qcdfrom_qcd[pfreg][qcdreg][2][ptbin])
                qcdratio_top_up[pfreg][qcdreg][ptbin] = np.sum(qcdratio_top_up[pfreg][qcdreg][ptbin]*qcdfrom_top[pfreg][qcdreg][2][ptbin])/np.sum(qcdfrom_top[pfreg][qcdreg][2][ptbin])
                qcdratio_wlnu_up[pfreg][qcdreg][ptbin] = np.sum(qcdratio_wlnu_up[pfreg][qcdreg][ptbin]*qcdfrom_wlnu[pfreg][qcdreg][2][ptbin])/np.sum(qcdfrom_wlnu[pfreg][qcdreg][2][ptbin])

    #print('qcdratio_sig',qcdratio_sig,qcdratio_sig_dn,qcdratio_sig_up)
    #print('qcdratio_qcd',qcdratio_qcd,qcdratio_qcd_dn,qcdratio_qcd_up)
    #print('qcdratio_top',qcdratio_top,qcdratio_top_dn,qcdratio_top_up)
    #print('qcdratio_wlnu',qcdratio_wlnu,qcdratio_wlnu_dn,qcdratio_wlnu_up)
    
    #SR               - [F] [L] [P]
    #SR (noAnti)      - [F]
    #mIso CR          - [F] [L] [P]
    #mIso CR (noAnti) - [F] [L] [P]
    #derive ratio from {SR (noAnti) [F]} / {mIso CR (noAnti) [F]} ==> qcdratio

    #CR - SR (R_miso)
    #F - L   (RQCD_FL) f(pT)
    #F - P   (RQCD_FP) f(pT)

    #QCD[F] * R_miso = SR[F]
    #QCD[F] * R_miso * RQCD_FL = SR[L]
    #QCD[F] * R_miso * RQCD_FP = SR[P]

    #check mIso CR [L] & [P]
    
    for ptbin in range(npt):
        for iregion,region in enumerate(['fail','loosepass','pass']):
            ch = rl.Channel("ptbin%d%s%s%s" % (ptbin, region, 'had%s'%leptype, year))
            model.addChannel(ch)

            isPass = region=='pass'
            isLoosePass = region=='loosepass'

            thehist = sig_hists["nom"]
    
            for sName in thehist.identifiers('sample'):
                if sName.name=='ignore':
                    continue
                if sName.name=='data_obs' and isPass and not usingData:
                    continue
                tempint = intRegion(thehist[sName],region,hptslice=ptslices[ptbin])
                templ = (tempint[0][lowmassbin:highmassbin], mtt.binning, mtt.name, tempint[1][lowmassbin:highmassbin])
                for iv,val in enumerate(templ[0]):
                    if val<0.:
                        templ[0][iv] = 0.
                if sName.name=='data_obs':
                    ch.setObservation(templ, read_sumw2=True)
                else:
                    if sName.name=='multijet':
                        qcdpred = np.stack([(qcd_from_data_sig["lowmet"][shaperegion[iregion]][ix][0]*qcdratio_sig[shaperegion[iregion]]["lowmet"][ptbin]*qcdratio_F["qcdnom"][region][ix]+qcd_from_data_qcd["nom"][shaperegion[iregion]][ix][0]*qcdratio_qcd[shaperegion[iregion]]["nom"][ptbin]*qcdratio_F["lowmet"][region][ix])/2. for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else [np.zeros_like(qcd_from_data_sig["lowmet"][region][0][0])])
                        qcdpred_dn = np.clip(np.sum(np.stack([qcd_from_data_sig["lowmet"][shaperegion[iregion]][ix][1] * qcdratio_sig[shaperegion[iregion]]["lowmet"][ptbin] * qcdratio_F["qcdnom"][region][ix] for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else [np.zeros_like(qcd_from_data_sig["lowmet"][region][0][0])]),axis=0),0.,None)
                        qcdpred_up = np.sum(np.stack([qcd_from_data_qcd["nom"][shaperegion[iregion]][ix][2] * qcdratio_qcd[shaperegion[iregion]]["nom"][ptbin] * qcdratio_F["lowmet"][region][ix] for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else [np.zeros_like(qcd_from_data_qcd["nom"][region][0][0])]),axis=0)
                        print('qcdpred',qcdpred)
                        qcdpred = np.sum(qcdpred,axis=0)
                        #print('qcdpred_dn',qcdpred_dn)
                        #print('qcdpred_up',qcdpred_up)
                        axs_shape[iregion,int(leptype=='el')].step((mttbins[lowmassbin:highmassbin-1]+mttbins[lowmassbin+1:highmassbin])[:mbinend]/2.,qcdpred[:mbinend], color='C%i'%(ptbin+0*iregion),label='%s, %i<pT<%i'%(region,ptbins[ptbin] if not ptbins[ptbin] is None else 0,ptbins[ptbin+1] if not ptbins[ptbin+1] is None else 9999),where="mid")
                        axs_shape[iregion,int(leptype=='el')].fill_between((mttbins[lowmassbin:highmassbin-1]+mttbins[lowmassbin+1:highmassbin])[:mbinend]/2., qcdpred_dn[:mbinend], qcdpred_up[:mbinend], alpha=0.5, color='C%i'%(ptbin+0*iregion),step="mid")
                        templ = (qcdpred, mtt.binning, mtt.name, np.ones_like(qcdpred)*1000000.)
                    if doHtt:
                        stype = rl.Sample.SIGNAL if 'htt125' in sName.name else rl.Sample.BACKGROUND
                    else:
                        stype = rl.Sample.SIGNAL if 'phi' in sName.name else rl.Sample.BACKGROUND
                    sample = rl.TemplateSample(ch.name + '_' + sName.name, stype, templ)
        
                    #sample.setParamEffect(jec, 1.05)
                    if sName.name in syst_dict:
                        nom = intRegion(thehist[sName],region,hptslice=ptslices[ptbin])[0][lowmassbin:highmassbin]
                        for syst in syst_dict[sName.name]:
                            syst_params = syst_dict[sName.name][syst]
                            dnmod = intRegion(thehist[sName],region,systematic=syst_params[1],hptslice=ptslices[ptbin])[0][lowmassbin:highmassbin]
                            upmod = intRegion(thehist[sName],region,systematic=syst_params[2],hptslice=ptslices[ptbin])[0][lowmassbin:highmassbin]
                            sample.setParamEffect(syst_params[0], np.divide(upmod, nom, out=np.ones_like(nom), where=nom>0.) if (upmod!=nom).all() else np.ones_like(nom)*1.001, np.divide(dnmod, nom, out=np.ones_like(nom), where=nom>0.) if (dnmod!=nom).all() else np.ones_like(nom)*0.999)

                    if sName.name=='top':
                        sample.setParamEffect(top_norm, 1.05)
                    if sName.name=='wlnu':
                        sample.setParamEffect(wlnu_norm, 1.10)
                    if sName.name=='multijet':
                        sample.setParamEffect(qcd_norm, 1.50)
                        qcd_shape_dn = np.divide(qcdpred_dn, qcdpred, out=np.ones_like(qcdpred), where=qcdpred>0.)
                        #qcd_shape_dn = qcd_shape_dn*np.sum(qcdpred)/np.sum(qcdpred_dn)
                        qcd_shape_up = np.divide(qcdpred_up, qcdpred, out=np.ones_like(qcdpred), where=qcdpred>0.)
                        #qcd_shape_up = qcd_shape_up*np.sum(qcdpred)/np.sum(qcdpred_up)
                        sample.setParamEffect(qcd_pass if isPass else qcd_loosepass if isLoosePass else qcd_fail, qcd_shape_dn, qcd_shape_up)
                    if sName.name=='vvqq':
                        sample.setParamEffect(vvqq_norm, 1.05)
                    if sName.name=='dy':
                        sample.setParamEffect(dy_norm, 1.05)
                    if sName.name=='htt125' or sName.name=='dy' or 'phi' in sName.name:
                        nom_full = intRegion(thehist[sName],region,hptslice=ptslices[ptbin])[0]
                        shift_dn = nom_full[lowmassbin+1:highmassbin+1 if highmassbin>0 else None]
                        shift_up = nom_full[lowmassbin-1:highmassbin-1]
                        #sample.setParamEffect(m_scale, np.divide(shift_dn, nom_full[lowmassbin:highmassbin], out=np.ones_like(nom_full[lowmassbin:highmassbin]), where=nom_full[lowmassbin:highmassbin]>0.), np.divide(shift_up, nom_full[lowmassbin:highmassbin], out=np.ones_like(nom_full[lowmassbin:highmassbin]), where=nom_full[lowmassbin:highmassbin]>0.))
                        #print(leptype,'htt125 ptbin',ptbin,region,nom_full[lowmassbin:highmassbin])
                    sample.setParamEffect(trig, 1.02)
                    sample.setParamEffect(lumi, 1.024)
            
                    ch.addSample(sample)
            if not usingData and isPass:
                ch.setObservation((np.zeros(len(mtt.binning)-1),mtt.binning, mtt.name))

    h125vals = [intRegion(thehist['htt125'],'pass',hptslice=ptslices[ptbin])[0][lowmassbin:highmassbin] for ptbin in range(npt)]
    if not doHtt:
        phi50vals = [intRegion(thehist['phitt50'],'pass',hptslice=ptslices[ptbin])[0][lowmassbin:highmassbin] for ptbin in range(npt)]
    dyvals = [intRegion(thehist['dy'],'pass',hptslice=ptslices[ptbin])[0][lowmassbin:highmassbin] for ptbin in range(npt)]
    wlnuvals = [intRegion(thehist['wlnu'],'pass',hptslice=ptslices[ptbin])[0][lowmassbin:highmassbin] for ptbin in range(npt)]
    if doHtt:
        probs = [np.prod([math.erfc(h125vals[ipt][ib]/(math.sqrt(abs(dyvals[ipt][ib]+wlnuvals[ipt][ib])) if dyvals[ipt][ib]+wlnuvals[ipt][ib]>2. else 2.)) for ib in range(len(h125vals[ipt]))]) for ipt in range(npt)]
    else:
        probs = [np.prod([math.erfc(phi50vals[ipt][ib]/(math.sqrt(abs(dyvals[ipt][ib]+wlnuvals[ipt][ib]+h125vals[ipt][ib])) if dyvals[ipt][ib]+wlnuvals[ipt][ib]+h125vals[ipt][ib]>2. else 2.)) for ib in range(len(phi50vals[ipt]))]) for ipt in range(npt)]
    print(leptype,'PROB',np.prod(probs),'->',erfinv(1.-np.prod(probs)),'\t\t',probs)

    for ptbin in range(npt):
        failCh = model['ptbin%dfail%s%s' % (ptbin,'had%s'%leptype,year)]
        passCh = model['ptbin%dpass%s%s' % (ptbin,'had%s'%leptype,year)]
        loosePassCh = model['ptbin%dloosepass%s%s' % (ptbin,'had%s'%leptype,year)]

        qcdpass = passCh['multijet']
        qcdloosepass = loosePassCh['multijet']
        qcdfail = failCh['multijet']
        qcdpass.setParamEffect(qcdnormSF, 1.50)
        qcdloosepass.setParamEffect(qcdnormSF, 1.50)
        qcdfail.setParamEffect(qcdnormSF, 1.50)
        #qcdpass.setParamEffect(qcdnormSF, 1*qcdnormSF)
        #qcdloosepass.setParamEffect(qcdnormSF, 1*qcdnormSF)
        #qcdfail.setParamEffect(qcdnormSF, 1*qcdnormSF)

        toppass = passCh['top']
        toploosepass = loosePassCh['top']
        topfail = failCh['top']
        topPF = (toploosepass.getExpectation(nominal=True).sum() + toppass.getExpectation(nominal=True).sum()) / topfail.getExpectation(nominal=True).sum()
        topLPF = toppass.getExpectation(nominal=True).sum() / toploosepass.getExpectation(nominal=True).sum()
        topRLPF = 1./(1.+(1./topLPF)) # P/(L+P)
        toppass.setParamEffect(topLeffSF, 1*topLeffSF)
        toploosepass.setParamEffect(topLeffSF, (1 - topLeffSF) * topLPF + 1)
        toploosepass.setParamEffect(topeffSF, 1*topeffSF)
        topfail.setParamEffect(topeffSF, (1 - topeffSF) * (1-topRLPF) * topPF + 1)
        #toppass.setParamEffect(topnormSF, 1*topnormSF)
        #toploosepass.setParamEffect(topnormSF, 1*topnormSF)
        #topfail.setParamEffect(topnormSF, 1*topnormSF)

        wlnupass = passCh['wlnu']
        wlnuloosepass = loosePassCh['wlnu']
        wlnufail = failCh['wlnu']
        wlnuPF = (wlnuloosepass.getExpectation(nominal=True).sum() + wlnupass.getExpectation(nominal=True).sum()) / wlnufail.getExpectation(nominal=True).sum()
        wlnuLPF = wlnupass.getExpectation(nominal=True).sum() / wlnuloosepass.getExpectation(nominal=True).sum()
        wlnuRLPF = 1./(1.+(1./wlnuLPF)) # P/(L+P)
        wlnupass.setParamEffect(wlnuLeffSF, 1*wlnuLeffSF)
        wlnuloosepass.setParamEffect(wlnuLeffSF, (1 - wlnuLeffSF) * wlnuLPF + 1)
        wlnuloosepass.setParamEffect(wlnueffSF, 1*wlnueffSF)
        wlnufail.setParamEffect(wlnueffSF, (1 - wlnueffSF) * (1-wlnuRLPF) * wlnuPF + 1)
        #wlnupass.setParamEffect(wlnunormSF, 1*wlnunormSF)
        #wlnuloosepass.setParamEffect(wlnunormSF, 1*wlnunormSF)
        #wlnufail.setParamEffect(wlnunormSF, 1*wlnunormSF)

        dypass = passCh['dy']
        dyloosepass = loosePassCh['dy']
        dyfail = failCh['dy']
        dyLP = dypass.getExpectation(nominal=True).sum() / dyloosepass.getExpectation(nominal=True).sum()
        httLP = passCh['htt125'].getExpectation(nominal=True).sum() / loosePassCh['htt125'].getExpectation(nominal=True).sum()
        if not doHtt:
            phittLP = {m:passCh['phitt%s'%m].getExpectation(nominal=True).sum() / loosePassCh['phitt%s'%m].getExpectation(nominal=True).sum() if loosePassCh['phitt%s'%m].getExpectation(nominal=True).sum() > 0. else 1. for m in masspoints}
        dypass.setParamEffect(dy_eff, 1*dy_eff)
        dyloosepass.setParamEffect(dy_eff, (1 - dy_eff) * dyLP + 1)
        passCh['htt125'].setParamEffect(dy_eff, 1*dy_eff)
        loosePassCh['htt125'].setParamEffect(dy_eff, (1 - dy_eff) * httLP + 1)
        if not doHtt:
            for m in masspoints:
                passCh['phitt%s'%m].setParamEffect(dy_eff, 1*dy_eff)
                loosePassCh['phitt%s'%m].setParamEffect(dy_eff, (1 - dy_eff) * phittLP[m] + 1)
        dypass.setParamEffect(rdy, 1*rdy)
        dyloosepass.setParamEffect(rdy, 1*rdy)
        dyfail.setParamEffect(rdy, 1*rdy)
        #dypass.setParamEffect(rdy, 1.05)
        #dyloosepass.setParamEffect(rdy, 1.05)
        #dyfail.setParamEffect(rdy, 1.05)
        if not doHtt:
            passCh['htt125'].setParamEffect(rh125, 1.10)
            loosePassCh['htt125'].setParamEffect(rh125, 1.10)
            failCh['htt125'].setParamEffect(rh125, 1.10)


    # Fill in top CR
    for iregion,region in enumerate(['pass', 'loosepass', 'fail']):
        ch = rl.Channel("topCR%s%s%s" % (region, 'had%s'%leptype, year))
        model.addChannel(ch)

        isPass = region=='pass'
        isLoosePass = region=='loosepass'

        thehist = top_cr_hists["nom"]

        for sName in thehist.identifiers('sample'):
            if sName.name=='ignore':
                continue
            tempint = intRegion(thehist[sName],region)
            if singleBinCR:
                templ = (np.array([np.sum(tempint[0][lowmassbin:highmassbin])]), mttone.binning, mttone.name, np.array([np.sum(tempint[1][lowmassbin:highmassbin])]))
            else:
                templ = (tempint[0][lowmassbin:highmassbin], mtt.binning, mtt.name, tempint[1][lowmassbin:highmassbin])
            for iv,val in enumerate(templ[0]):
                if val<0.:
                    templ[0][iv] = 0.
            if sName.name=='data_obs':
                ch.setObservation(templ, read_sumw2=True)
            else:
                if sName.name=='multijet':
                    qcdpred = np.sum(np.stack([(qcd_from_data_top["lowmet"][shaperegion[iregion]][ix][0]*qcdratio_top[shaperegion[iregion]]["lowmet"][ptbin]*qcdratio_F["qcdnom"][region][ix]+qcd_from_data_top["nom"][shaperegion[iregion]][ix][0]*qcdratio_top[shaperegion[iregion]]["nom"][ptbin]*qcdratio_F["lowmet"][region][ix])/2. for ix in range(len(qcd_from_data_top["lowmet"][region]))]),axis=0)
                    qcdpred_dn = np.clip(np.sum(np.stack([qcd_from_data_top["lowmet"][shaperegion[iregion]][ix][1] * qcdratio_top[shaperegion[iregion]]["lowmet"][ptbin] * qcdratio_F["qcdnom"][region][ix] for ix in range(len(qcd_from_data_top["lowmet"][region]))]),axis=0),0.,None)
                    qcdpred_up = np.sum(np.stack([qcd_from_data_qcd["nom"][shaperegion[iregion]][ix][2] * qcdratio_qcd[shaperegion[iregion]]["nom"][ptbin] * qcdratio_F["lowmet"][region][ix] for ix in range(len(qcd_from_data_qcd["nom"][region]))]),axis=0)
                    if singleBinCR:
                        templ = (np.array([np.sum(qcdpred)]), mttone.binning, mttone.name, np.array([1000000.]))
                    else:
                        templ = (qcdpred, mtt.binning, mtt.name, np.ones_like(qcdpred)*1000000.)
                if doHtt:
                    stype = rl.Sample.SIGNAL if 'htt125' in sName.name else rl.Sample.BACKGROUND
                else:
                    stype = rl.Sample.SIGNAL if 'phi' in sName.name else rl.Sample.BACKGROUND
                sample = rl.TemplateSample(ch.name + '_' + sName.name, stype, templ)
    
                #sample.setParamEffect(jec, 1.05)
                if sName.name in syst_dict:
                    nom_base = intRegion(thehist[sName],region)[0][lowmassbin:highmassbin]
                    for syst in syst_dict[sName.name]:
                        syst_params = syst_dict[sName.name][syst]
                        dnmod = intRegion(thehist[sName],region,systematic=syst_params[1])[0][lowmassbin:highmassbin]
                        upmod = intRegion(thehist[sName],region,systematic=syst_params[2])[0][lowmassbin:highmassbin]
                        if singleBinCR:
                            dnmod = np.array([np.sum(dnmod)])
                            nom = np.array([np.sum(nom_base)])
                            upmod = np.array([np.sum(upmod)])
                        else:
                            nom = nom_base
                        sample.setParamEffect(syst_params[0], np.divide(upmod, nom, out=np.ones_like(nom), where=nom>0.) if (upmod!=nom).all() else np.ones_like(nom)*1.001, np.divide(dnmod, nom, out=np.ones_like(nom), where=nom>0.) if (dnmod!=nom).all() else np.ones_like(nom)*0.999)

                if sName.name=='top':
                    sample.setParamEffect(top_norm, 1.05)
                if sName.name=='wlnu':
                    sample.setParamEffect(wlnu_norm, 1.10)
                if sName.name=='multijet':
                    sample.setParamEffect(qcd_norm, 1.50)
                    if singleBinCR:
                        err_dn = np.sqrt(np.sum((qcdpred - qcdpred_dn)*(qcdpred - qcdpred_dn)))
                        err_up = np.sqrt(np.sum((qcdpred_up - qcdpred)*(qcdpred_up - qcdpred)))
                        qcdpred = np.array([np.sum(qcdpred)])
                        qcdpred_dn = qcdpred - err_dn
                        qcdpred_up = qcdpred + err_up
                    #sample.setParamEffect(qcdtop_pass if isPass else qcdtop_loosepass if isLoosePass else qcdtop_fail, np.divide(qcdpred_dn, qcdpred, out=np.ones_like(qcdpred), where=qcdpred>0.), np.divide(qcdpred_up, qcdpred, out=np.ones_like(qcdpred), where=qcdpred>0.))
                    sample.setParamEffect(qcd_pass if isPass else qcd_loosepass if isLoosePass else qcd_fail, np.divide(qcdpred_dn, qcdpred, out=np.ones_like(qcdpred), where=qcdpred>0.), np.divide(qcdpred_up, qcdpred, out=np.ones_like(qcdpred), where=qcdpred>0.))
                if sName.name=='vvqq':
                    sample.setParamEffect(vvqq_norm, 1.05)
                if sName.name=='dy':
                    sample.setParamEffect(dy_norm, 1.05)
                if sName.name=='htt125' or sName.name=='dy' or 'phi' in sName.name:
                    nom_full = intRegion(thehist[sName],region)[0]
                    shift_dn = nom_full[lowmassbin+1:highmassbin+1 if highmassbin>0 else None]
                    shift_up = nom_full[lowmassbin-1:highmassbin-1]
                    #sample.setParamEffect(m_scale, np.divide(shift_dn, nom_full[lowmassbin:highmassbin], out=np.ones_like(nom_full[lowmassbin:highmassbin]), where=nom_full[lowmassbin:highmassbin]>0.), np.divide(shift_up, nom_full[lowmassbin:highmassbin], out=np.ones_like(nom_full[lowmassbin:highmassbin]), where=nom_full[lowmassbin:highmassbin]>0.))
                sample.setParamEffect(trig, 1.02)
                sample.setParamEffect(lumi, 1.024)
        
                ch.addSample(sample)

    qcdpass = model['topCRpasshad%s%s'%(leptype,year)]['multijet']
    qcdloosepass = model['topCRloosepasshad%s%s'%(leptype,year)]['multijet']
    qcdfail = model['topCRfailhad%s%s'%(leptype,year)]['multijet']
    #qcdpass.setParamEffect(qcdnormSF_top, 1.50)
    #qcdloosepass.setParamEffect(qcdnormSF_top, 1.50)
    #qcdfail.setParamEffect(qcdnormSF_top, 1.50)
    qcdpass.setParamEffect(qcdnormSF_top, 1.50)
    qcdloosepass.setParamEffect(qcdnormSF_top, 1.50)
    qcdfail.setParamEffect(qcdnormSF_top, 1.50)
    #qcdpass.setParamEffect(qcdnormSF_top, 1*qcdnormSF_top)
    #qcdloosepass.setParamEffect(qcdnormSF_top, 1*qcdnormSF_top)
    #qcdfail.setParamEffect(qcdnormSF_top, 1*qcdnormSF_top)

    toppass = model['topCRpasshad%s%s'%(leptype,year)]['top']
    toploosepass = model['topCRloosepasshad%s%s'%(leptype,year)]['top']
    topfail = model['topCRfailhad%s%s'%(leptype,year)]['top']
    topPF = (toploosepass.getExpectation(nominal=True).sum() + toppass.getExpectation(nominal=True).sum()) / topfail.getExpectation(nominal=True).sum()
    topLPF = toppass.getExpectation(nominal=True).sum() / toploosepass.getExpectation(nominal=True).sum()
    topRLPF = 1./(1.+(1./topLPF)) # P/(L+P)
    toppass.setParamEffect(topLeffSF, 1*topLeffSF)
    toploosepass.setParamEffect(topLeffSF, (1 - topLeffSF) * topLPF + 1)
    toploosepass.setParamEffect(topeffSF, 1*topeffSF)
    topfail.setParamEffect(topeffSF, (1 - topeffSF) * (1-topRLPF) * topPF + 1)
    #toppass.setParamEffect(topnormSF, 1*topnormSF)
    #toploosepass.setParamEffect(topnormSF, 1*topnormSF)
    #topfail.setParamEffect(topnormSF, 1*topnormSF)

    wlnupass = model['topCRpasshad%s%s'%(leptype,year)]['wlnu']
    wlnuloosepass = model['topCRloosepasshad%s%s'%(leptype,year)]['wlnu']
    wlnufail = model['topCRfailhad%s%s'%(leptype,year)]['wlnu']
    wlnuPF = (wlnuloosepass.getExpectation(nominal=True).sum() + wlnupass.getExpectation(nominal=True).sum()) / wlnufail.getExpectation(nominal=True).sum()
    wlnuLPF = wlnupass.getExpectation(nominal=True).sum() / wlnuloosepass.getExpectation(nominal=True).sum()
    wlnuRLPF = 1./(1.+(1./wlnuLPF)) # P/(L+P)
    wlnupass.setParamEffect(wlnuLeffSF, 1*wlnuLeffSF)
    wlnuloosepass.setParamEffect(wlnuLeffSF, (1 - wlnuLeffSF) * wlnuLPF + 1)
    wlnuloosepass.setParamEffect(wlnueffSF, 1*wlnueffSF)
    wlnufail.setParamEffect(wlnueffSF, (1 - wlnueffSF) * (1-wlnuRLPF) * wlnuPF + 1)
    #wlnupass.setParamEffect(wlnunormSF, 1*wlnunormSF)
    #wlnuloosepass.setParamEffect(wlnunormSF, 1*wlnunormSF)
    #wlnufail.setParamEffect(wlnunormSF, 1*wlnunormSF)

    dypass = model['topCRpasshad%s%s'%(leptype,year)]['dy']
    dyloosepass = model['topCRloosepasshad%s%s'%(leptype,year)]['dy']
    dyfail = model['topCRfailhad%s%s'%(leptype,year)]['dy']
    dyLP = dypass.getExpectation(nominal=True).sum() / dyloosepass.getExpectation(nominal=True).sum()
    httLP = model['topCRpasshad%s%s'%(leptype,year)]['htt125'].getExpectation(nominal=True).sum() / model['topCRloosepasshad%s%s'%(leptype,year)]['htt125'].getExpectation(nominal=True).sum()
    if not doHtt:
        phittLP = {m:model['topCRpasshad%s%s'%(leptype,year)]['phitt%s'%m].getExpectation(nominal=True).sum() / model['topCRloosepasshad%s%s'%(leptype,year)]['phitt%s'%m].getExpectation(nominal=True).sum() if model['topCRloosepasshad%s%s'%(leptype,year)]['phitt%s'%m].getExpectation(nominal=True).sum() > 0. else 1. for m in masspoints}
    dypass.setParamEffect(dy_eff, 1*dy_eff)
    dyloosepass.setParamEffect(dy_eff, (1 - dy_eff)* dyLP + 1)
    model['topCRpasshad%s%s'%(leptype,year)]['htt125'].setParamEffect(dy_eff, 1*dy_eff)
    model['topCRloosepasshad%s%s'%(leptype,year)]['htt125'].setParamEffect(dy_eff, (1 - dy_eff) * httLP + 1)
    if not doHtt:
        for m in masspoints:
            model['topCRpasshad%s%s'%(leptype,year)]['phitt%s'%m].setParamEffect(dy_eff, 1*dy_eff)
            model['topCRloosepasshad%s%s'%(leptype,year)]['phitt%s'%m].setParamEffect(dy_eff, (1 - dy_eff) * phittLP[m] + 1)
    dypass.setParamEffect(rdy, 1*rdy)
    dyloosepass.setParamEffect(rdy, 1*rdy)
    dyfail.setParamEffect(rdy, 1*rdy)
    #dypass.setParamEffect(rdy, 1.05)
    #dyloosepass.setParamEffect(rdy, 1.05)
    #dyfail.setParamEffect(rdy, 1.05)
    if not doHtt:
        model['topCRpasshad%s%s'%(leptype,year)]['htt125'].setParamEffect(rh125, 1.10)
        model['topCRloosepasshad%s%s'%(leptype,year)]['htt125'].setParamEffect(rh125, 1.10)
        model['topCRfailhad%s%s'%(leptype,year)]['htt125'].setParamEffect(rh125, 1.10)

    # Fill in wlnu CR
    for iregion,region in enumerate(['pass', 'fail', 'loosepass']):
        ch = rl.Channel("wlnuCR%s%s%s" % (region, 'had%s'%leptype, year))
        model.addChannel(ch)

        isPass = region=='pass'
        isLoosePass = region=='loosepass'

        thehist = wlnu_cr_hists["nom"]

        for sName in thehist.identifiers('sample'):
            if sName.name=='ignore':
                continue
            tempint = intRegion(thehist[sName],region)
            if singleBinCR:
                templ = (np.array([np.sum(tempint[0][lowmassbin:highmassbin])]), mttone.binning, mttone.name, np.array([np.sum(tempint[1][lowmassbin:highmassbin])]))
            else:
                templ = (tempint[0][lowmassbin:highmassbin], mtt.binning, mtt.name, tempint[1][lowmassbin:highmassbin])
            for iv,val in enumerate(templ[0]):
                if val<0.:
                    templ[0][iv] = 0.
            if sName.name=='data_obs':
                ch.setObservation(templ, read_sumw2=True)
            else:
                if sName.name=='multijet':
                    qcdpred = np.sum(np.stack([(qcd_from_data_wlnu["lowmet"][shaperegion[iregion]][ix][0]*qcdratio_wlnu[shaperegion[iregion]]["lowmet"][ptbin]*qcdratio_F["qcdnom"][region][ix]+qcd_from_data_qcd["nom"][region][ix][0]*qcdratio_qcd[shaperegion[iregion]]["nom"][ptbin]*qcdratio_F["lowmet"][region][ix])/2. for ix in range(len(qcd_from_data_wlnu["lowmet"][region]))]),axis=0)
                    qcdpred_dn = np.clip(np.sum(np.stack([qcd_from_data_wlnu["lowmet"][shaperegion[iregion]][ix][1] * qcdratio_wlnu[shaperegion[iregion]]["lowmet"][ptbin] * qcdratio_F["qcdnom"][region][ix] for ix in range(len(qcd_from_data_wlnu["lowmet"][region]))]),axis=0),0.,None)
                    qcdpred_up = np.sum(np.stack([qcd_from_data_qcd["nom"][shaperegion[iregion]][ix][2] * qcdratio_qcd[shaperegion[iregion]]["nom"][ptbin] * qcdratio_F["lowmet"][region][ix] for ix in range(len(qcd_from_data_qcd["nom"][region]))]),axis=0)
                    if singleBinCR:
                        templ = (np.array([np.sum(qcdpred)]), mttone.binning, mttone.name, np.array([1000000.]))
                    else:
                        templ = (qcdpred, mtt.binning, mtt.name, np.ones_like(qcdpred)*1000000.)
                if doHtt:
                    stype = rl.Sample.SIGNAL if 'htt125' in sName.name else rl.Sample.BACKGROUND
                else:
                    stype = rl.Sample.SIGNAL if 'phi' in sName.name else rl.Sample.BACKGROUND
                sample = rl.TemplateSample(ch.name + '_' + sName.name, stype, templ)
    
                #sample.setParamEffect(jec, 1.05)
                if sName.name in syst_dict:
                    nom_base = intRegion(thehist[sName],region)[0][lowmassbin:highmassbin]
                    for syst in syst_dict[sName.name]:
                        syst_params = syst_dict[sName.name][syst]
                        dnmod = intRegion(thehist[sName],region,systematic=syst_params[1])[0][lowmassbin:highmassbin]
                        upmod = intRegion(thehist[sName],region,systematic=syst_params[2])[0][lowmassbin:highmassbin]
                        if singleBinCR:
                            dnmod = np.array([np.sum(dnmod)])
                            nom = np.array([np.sum(nom_base)])
                            upmod = np.array([np.sum(upmod)])
                        else:
                            nom = nom_base
                        sample.setParamEffect(syst_params[0], np.divide(upmod, nom, out=np.ones_like(nom), where=nom>0.) if (upmod!=nom).all() else np.ones_like(nom)*1.001, np.divide(dnmod, nom, out=np.ones_like(nom), where=nom>0.) if (dnmod!=nom).all() else np.ones_like(nom)*0.999)

                if sName.name=='top':
                    sample.setParamEffect(top_norm, 1.05)
                if sName.name=='wlnu':
                    sample.setParamEffect(wlnu_norm, 1.10)
                if sName.name=='multijet':
                    sample.setParamEffect(qcd_norm, 1.50)
                    if singleBinCR:
                        err_dn = np.sqrt(np.sum((qcdpred - qcdpred_dn)*(qcdpred - qcdpred_dn)))
                        err_up = np.sqrt(np.sum((qcdpred_up - qcdpred)*(qcdpred_up - qcdpred)))
                        qcdpred = np.array([np.sum(qcdpred)])
                        qcdpred_dn = qcdpred - err_dn
                        qcdpred_up = qcdpred + err_up
                    #sample.setParamEffect(qcdwlnu_pass if isPass else qcdwlnu_loosepass if isLoosePass else qcdwlnu_fail, np.divide(qcdpred_dn, qcdpred, out=np.ones_like(qcdpred), where=qcdpred>0.), np.divide(qcdpred_up, qcdpred, out=np.ones_like(qcdpred), where=qcdpred>0.))
                    sample.setParamEffect(qcd_pass if isPass else qcd_loosepass if isLoosePass else qcd_fail, np.divide(qcdpred_dn, qcdpred, out=np.ones_like(qcdpred), where=qcdpred>0.), np.divide(qcdpred_up, qcdpred, out=np.ones_like(qcdpred), where=qcdpred>0.))
                if sName.name=='vvqq':
                    sample.setParamEffect(vvqq_norm, 1.05)
                if sName.name=='dy':
                    sample.setParamEffect(dy_norm, 1.05)
                if sName.name=='htt125' or sName.name=='dy' or 'phi' in sName.name:
                    nom_full = intRegion(thehist[sName],region)[0]
                    shift_dn = nom_full[lowmassbin+1:highmassbin+1 if highmassbin>0 else None]
                    shift_up = nom_full[lowmassbin-1:highmassbin-1]
                    #sample.setParamEffect(m_scale, np.divide(shift_dn, nom_full[lowmassbin:highmassbin], out=np.ones_like(nom_full[lowmassbin:highmassbin]), where=nom_full[lowmassbin:highmassbin]>0.), np.divide(shift_up, nom_full[lowmassbin:highmassbin], out=np.ones_like(nom_full[lowmassbin:highmassbin]), where=nom_full[lowmassbin:highmassbin]>0.))
                sample.setParamEffect(trig, 1.02)
                sample.setParamEffect(lumi, 1.024)
        
                ch.addSample(sample)

    qcdpass = model['wlnuCRpasshad%s%s'%(leptype,year)]['multijet']
    qcdloosepass = model['wlnuCRloosepasshad%s%s'%(leptype,year)]['multijet']
    qcdfail = model['wlnuCRfailhad%s%s'%(leptype,year)]['multijet']
    #qcdpass.setParamEffect(qcdnormSF_wlnu, 1.50)
    #qcdloosepass.setParamEffect(qcdnormSF_wlnu, 1.50)
    #qcdfail.setParamEffect(qcdnormSF_wlnu, 1.50)
    qcdpass.setParamEffect(qcdnormSF_wlnu, 1.50)
    qcdloosepass.setParamEffect(qcdnormSF_wlnu, 1.50)
    qcdfail.setParamEffect(qcdnormSF_wlnu, 1.50)
    #qcdpass.setParamEffect(qcdnormSF_wlnu, 1*qcdnormSF_wlnu)
    #qcdloosepass.setParamEffect(qcdnormSF_wlnu, 1*qcdnormSF_wlnu)
    #qcdfail.setParamEffect(qcdnormSF_wlnu, 1*qcdnormSF_wlnu)

    toppass = model['wlnuCRpasshad%s%s'%(leptype,year)]['top']
    toploosepass = model['wlnuCRloosepasshad%s%s'%(leptype,year)]['top']
    topfail = model['wlnuCRfailhad%s%s'%(leptype,year)]['top']
    topPF = (toploosepass.getExpectation(nominal=True).sum() + toppass.getExpectation(nominal=True).sum()) / topfail.getExpectation(nominal=True).sum()
    topLPF = toppass.getExpectation(nominal=True).sum() / toploosepass.getExpectation(nominal=True).sum()
    topRLPF = 1./(1.+(1./topLPF)) # P/(L+P)
    toppass.setParamEffect(topLeffSF, 1*topLeffSF)
    toploosepass.setParamEffect(topLeffSF, (1 - topLeffSF) * topLPF + 1)
    toploosepass.setParamEffect(topeffSF, 1*topeffSF)
    topfail.setParamEffect(topeffSF, (1 - topeffSF) * (1-topRLPF) * topPF + 1)
    #toppass.setParamEffect(topnormSF, 1*topnormSF)
    #toploosepass.setParamEffect(topnormSF, 1*topnormSF)
    #topfail.setParamEffect(topnormSF, 1*topnormSF)

    wlnupass = model['wlnuCRpasshad%s%s'%(leptype,year)]['wlnu']
    wlnuloosepass = model['wlnuCRloosepasshad%s%s'%(leptype,year)]['wlnu']
    wlnufail = model['wlnuCRfailhad%s%s'%(leptype,year)]['wlnu']
    wlnuPF = (wlnuloosepass.getExpectation(nominal=True).sum() + wlnupass.getExpectation(nominal=True).sum()) / wlnufail.getExpectation(nominal=True).sum()
    wlnuLPF = wlnupass.getExpectation(nominal=True).sum() / wlnuloosepass.getExpectation(nominal=True).sum()
    wlnuRLPF = 1./(1.+(1./wlnuLPF)) # P/(L+P)
    wlnupass.setParamEffect(wlnuLeffSF, 1*wlnuLeffSF)
    wlnuloosepass.setParamEffect(wlnuLeffSF, (1 - wlnuLeffSF) * wlnuLPF + 1)
    wlnuloosepass.setParamEffect(wlnueffSF, 1*wlnueffSF)
    wlnufail.setParamEffect(wlnueffSF, (1 - wlnueffSF) * (1-wlnuRLPF) * wlnuPF + 1)
    #wlnupass.setParamEffect(wlnunormSF, 1*wlnunormSF)
    #wlnuloosepass.setParamEffect(wlnunormSF, 1*wlnunormSF)
    #wlnufail.setParamEffect(wlnunormSF, 1*wlnunormSF)

    dypass = model['wlnuCRpasshad%s%s'%(leptype,year)]['dy']
    dyloosepass = model['wlnuCRloosepasshad%s%s'%(leptype,year)]['dy']
    dyfail = model['wlnuCRfailhad%s%s'%(leptype,year)]['dy']
    dyLP = dypass.getExpectation(nominal=True).sum() / dyloosepass.getExpectation(nominal=True).sum()
    httLP = model['wlnuCRpasshad%s%s'%(leptype,year)]['htt125'].getExpectation(nominal=True).sum() / model['wlnuCRloosepasshad%s%s'%(leptype,year)]['htt125'].getExpectation(nominal=True).sum()
    if not doHtt:
        phittLP = {m:model['wlnuCRpasshad%s%s'%(leptype,year)]['phitt%s'%m].getExpectation(nominal=True).sum() / model['wlnuCRloosepasshad%s%s'%(leptype,year)]['phitt%s'%m].getExpectation(nominal=True).sum() if model['wlnuCRloosepasshad%s%s'%(leptype,year)]['phitt%s'%m].getExpectation(nominal=True).sum() > 0. else 1. for m in masspoints}
    dypass.setParamEffect(dy_eff, 1*dy_eff)
    dyloosepass.setParamEffect(dy_eff, (1 - dy_eff) * dyLP + 1)
    model['wlnuCRpasshad%s%s'%(leptype,year)]['htt125'].setParamEffect(dy_eff, 1*dy_eff)
    model['wlnuCRloosepasshad%s%s'%(leptype,year)]['htt125'].setParamEffect(dy_eff, (1 - dy_eff) * httLP + 1)
    if not doHtt:
        for m in masspoints:
            model['wlnuCRpasshad%s%s'%(leptype,year)]['phitt%s'%m].setParamEffect(dy_eff, 1*dy_eff)
            model['wlnuCRloosepasshad%s%s'%(leptype,year)]['phitt%s'%m].setParamEffect(dy_eff, (1 - dy_eff) * phittLP[m] + 1)
    dypass.setParamEffect(rdy, 1*rdy)
    dyloosepass.setParamEffect(rdy, 1*rdy)
    dyfail.setParamEffect(rdy, 1*rdy)
    #dypass.setParamEffect(rdy, 1.05)
    #dyloosepass.setParamEffect(rdy, 1.05)
    #dyfail.setParamEffect(rdy, 1.05)
    if not doHtt:
        model['wlnuCRpasshad%s%s'%(leptype,year)]['htt125'].setParamEffect(rh125, 1.10)
        model['wlnuCRloosepasshad%s%s'%(leptype,year)]['htt125'].setParamEffect(rh125, 1.10)
        model['wlnuCRfailhad%s%s'%(leptype,year)]['htt125'].setParamEffect(rh125, 1.10)

    with open(os.path.join("%s/%s"%(str(tmpdir),label), 'had%sModel.pkl'%leptype), "wb") as fout:
        pickle.dump(model, fout, protocol=2)

    #model.renderCombine(os.path.join("%s/%s"%(str(tmpdir),label), 'had%sModel'%leptype))

def createHadHad(sig_hist, top_cr_hist, wlnu_cr_hist, qcd_cr_hist, sig_faillep, top_cr_faillep, wlnu_cr_faillep, qcd_cr_faillep, sig_faildphi, top_cr_faildphi, wlnu_cr_faildphi, qcd_cr_faildphi, var_name, mttbins, ptbins, tmpdir, label, usingData, nnCut_met, nnCut_met_loose, metCut, lowMetCut, h_pt_min, masspoints, shaperegion, year, doHtt, noSyst):

    #could be made region-dependent
    samp_combinations = {
        'data_obs'   : ['data'],
        'top'        : ['tt-had', 'tt-semilep', 'tt-dilep', 'st'],
        'htt125'     : ['h125'],
        'multijet'   : ['qcd'],
        'dy'         : ['zem','ztt'],
        'wlnu'       : ['wjets'],
        'vvqq'        : ['vv','vqq'],
        'ignore'     : [],
    }
    for m in masspoints:
        if doHtt:
            samp_combinations['ignore'].append('phi%s'%m)
        else:
            samp_combinations['phitt%s'%m] = ['phi%s'%m]

    for ident in sig_hist.identifiers("process"):
        if 'phi' in ident.name:
            ignoreSig = True
            for m in masspoints:
                if ident.name=='phi%s'%m: 
                    ignoreSig = False
            if ignoreSig:
                print('adding %s to ignore list'%ident.name)
                samp_combinations['ignore'].append(ident.name)

    sig_hists = {}
    top_cr_hists = {}
    wlnu_cr_hists = {}
    qcd_cr_hists = {}

    print(samp_combinations)
    
    sig_hists["lowmet_faildphi"] = sig_faildphi.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,metCut),'under')
    top_cr_hists["lowmet_faildphi"] = top_cr_faildphi.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,metCut),'under')
    wlnu_cr_hists["lowmet_faildphi"] = wlnu_cr_faildphi.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,metCut),'under')
    qcd_cr_hists["lowmet_faildphi"] = qcd_cr_faildphi.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,metCut),'under')

    sig_hists["faildphi"] = sig_faildphi.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(metCut,None),'over')
    top_cr_hists["faildphi"] = top_cr_faildphi.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,None),'over')
    wlnu_cr_hists["faildphi"] = wlnu_cr_faildphi.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,None),'over')
    qcd_cr_hists["faildphi"] = qcd_cr_faildphi.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,None),'over')

    sig_hists["lowmet_faillep"] = sig_faillep.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,metCut),'under')
    top_cr_hists["lowmet_faillep"] = top_cr_faillep.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,metCut),'under')
    wlnu_cr_hists["lowmet_faillep"] = wlnu_cr_faillep.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,metCut),'under')
    qcd_cr_hists["lowmet_faillep"] = qcd_cr_faillep.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,metCut),'under')

    sig_hists["faillep"] = sig_faillep.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(metCut,None),'over')
    top_cr_hists["faillep"] = top_cr_faillep.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,None),'over')
    wlnu_cr_hists["faillep"] = wlnu_cr_faillep.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,None),'over')
    qcd_cr_hists["faillep"] = qcd_cr_faillep.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,None),'over')

    sig_hists["lowmet"] = sig_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,metCut),'under')
    top_cr_hists["lowmet"] = top_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,metCut),'under')
    wlnu_cr_hists["lowmet"] = wlnu_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,metCut),'under')
    qcd_cr_hists["lowmet"] = qcd_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,metCut),'under')

    sig_hists["nom"] = sig_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(metCut,None),'over')
    top_cr_hists["nom"] = top_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,None),'over')
    wlnu_cr_hists["nom"] = wlnu_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,None),'over')
    qcd_cr_hists["nom"] = qcd_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('met_pt',slice(lowMetCut,None),'over')

    #jec = rl.NuisanceParameter('CMS_jec', 'lnN')
    top_norm = rl.NuisanceParameter('CMS_top_norm', 'lnN')
    wlnu_norm = rl.NuisanceParameter('CMS_wlnu_norm', 'lnN')
    qcd_norm = rl.NuisanceParameter('CMS_qcd_norm', 'lnN')
    dy_norm = rl.NuisanceParameter('CMS_dy_norm', 'lnN')
    vvqq_norm = rl.NuisanceParameter('CMS_vvqq_norm', 'lnN')
    trig = rl.NuisanceParameter('CMS_trig_hadhad', 'lnN')
    lumi = rl.NuisanceParameter('CMS_lumi', 'lnN')

    qcd_fail = rl.NuisanceParameter('qcd_Rfail_hadhad', 'shape')
    qcd_loosepass = rl.NuisanceParameter('qcd_Rloosepass_hadhad', 'shape')
    qcd_pass = rl.NuisanceParameter('qcd_Rpass_hadhad', 'shape')
    #qcdnormSF = rl.IndependentParameter('qcdnormSF_hadhad', 1., 0, 10)
    qcdnormSF = rl.NuisanceParameter('qcdnormSF_hadhad', 'lnN')
    qcdtop_fail = rl.NuisanceParameter('qcd_top_Rfail_hadhad', 'shape')
    qcdtop_loosepass = rl.NuisanceParameter('qcd_top_Rloosepass_hadhad', 'shape')
    qcdtop_pass = rl.NuisanceParameter('qcd_top_Rpass_hadhad', 'shape')
    #qcdnormSF_top = rl.IndependentParameter('qcdnormSF_top_hadhad', 1., 0, 10)
    #qcdnormSF_top = rl.NuisanceParameter('qcdnormSF_top_hadhad', 'lnN')
    qcdwlnu_fail = rl.NuisanceParameter('qcd_wlnu_Rfail_hadhad', 'shape')
    qcdwlnu_loosepass = rl.NuisanceParameter('qcd_wlnu_Rloosepass_hadhad', 'shape')
    qcdwlnu_pass = rl.NuisanceParameter('qcd_wlnu_Rpass_hadhad', 'shape')
    #qcdnormSF_wlnu = rl.IndependentParameter('qcdnormSF_wlnu_hadhad', 1., 0, 10)
    #qcdnormSF_wlnu = rl.NuisanceParameter('qcdnormSF_wlnu_hadhad', 'lnN')

    #m_scale = rl.NuisanceParameter('massscale_hadhad', 'shape')

    topeffSF = rl.IndependentParameter('topeffSF_hadhad', 1., 0, 10)
    #topnormSF = rl.IndependentParameter('topnormSF_hadhad', 1., 0, 10)
    wlnueffSF = rl.IndependentParameter('wlnueffSF_hadhad', 1., 0, 10)
    #wlnunormSF = rl.IndependentParameter('wlnunormSF_hadhad', 1., 0, 10)

    topLeffSF = rl.IndependentParameter('topLeffSF_hadhad', 1., 0, 10)
    wlnuLeffSF = rl.IndependentParameter('wlnuLeffSF_hadhad', 1., 0, 10)

    rdy = rl.IndependentParameter('r_dy_hadhad', 1., 0, 10)
    #rdy = rl.NuisanceParameter('r_dy_hadhad', 'lnN')
    rh125 = rl.NuisanceParameter('r_h125_hadhad', 'lnN')
    dy_eff = rl.IndependentParameter('dy_eff_hadhad', 1., 0, 10)

    toppt = rl.NuisanceParameter('toppt', 'shape')
    uescale = rl.NuisanceParameter('uescale', 'shape')
    jescale = rl.NuisanceParameter('jescale', 'shape')
    jeresol = rl.NuisanceParameter('jeresol', 'shape')
    l1prefire = rl.NuisanceParameter('l1prefire', 'shape')
    syst_dict = {
        samp:{
            'uescale':[uescale,'UESDown','UESUp'],
            'jescale':[jescale,'JESDown','JESUp'],
            'jeresol':[jeresol,'JERDown','JERUp'],
        } for samp in samp_combinations if samp not in ['data_obs','ignore']
    }
    syst_dict['top']['toppt'] = [toppt,'nominal','TopPtReweightUp']
    if year not in ['2018']:
        for samp in syst_dict:
            syst_dict[samp]['l1prefire'] = [l1prefire,'L1PreFiringDown','L1PreFiringUp']
    if noSyst:
        syst_dict = {}

    npt = len(ptbins) - 1
    ptslices = [slice(ptbins[ipt],ptbins[ipt+1]) for ipt in range(npt)]
    mtt = rl.Observable(var_name, mttbins[lowmassbin:highmassbin])
    mttone = rl.Observable(var_name+'_one', np.array([mttbins[1],mttbins[-2]]))

    qcdedges = [sig_hists["nom"].axis('h_pt').edges()[ix] for ix in range(len(sig_hists["nom"].axis('h_pt').edges()))]
    qcdbins = [[ix for ix in range(len(qcdedges)-1) if qcdedges[ix]>=(ptbins[ipt] if ptbins[ipt] is not None else -999999.) and qcdedges[ix+1]<=(ptbins[ipt+1] if ptbins[ipt+1] is not None else 999999.)] for ipt in range(npt)]
    print('qcdbins',qcdbins)
    
    model = rl.Model("hadhadModel")

    def intRegion(inhist,theregion,systematic='nominal',samplelist=None,debug=False,hptslice=slice(h_pt_min,None),mslice=None,mrebin=mttbins,nnCut=nnCut_met,nnCut_loose=nnCut_met_loose):
        if theregion=='pass':
            theslice = slice(nnCut,None)
            overflow_str = 'over'
        elif theregion=='loosepass':
            theslice = slice(nnCut_loose,nnCut)
            overflow_str = 'none'
        elif theregion=='fail':
            theslice = slice(None,nnCut_loose)
            overflow_str = 'under'
        else:
            print("Unknown region",theregion)
            return

        the_int = inhist.integrate('nn_disc',theslice,overflow_str).integrate('systematic',systematic)
        if debug:
            print('\t',nnCut_loose,nnCut)
            print('\t',overflow_str,the_int.values())

        if hptslice is not None:
            if hptslice.start is not None and hptslice.stop is not None:
                overflow_str = 'none'
            elif hptslice.start is None and hptslice.stop is not None:
                overflow_str = 'under'
            elif hptslice.start is not None and hptslice.stop is None:
                overflow_str = 'over'
            else:
                overflow_str = 'allnan'
            the_int = the_int.integrate('h_pt',hptslice,overflow_str)
            if debug:
                print('\t',hptslice)
                print('\t',overflow_str,the_int.values())
        if mrebin is not None:
            the_int = the_int.rebin('massreg',hist.Bin('massreg','massreg',mrebin))
        if mslice is not None:
            if mslice.start is not None and mslice.stop is not None:
                overflow_str = 'none'
            elif mslice.start is None and mslice.stop is not None:
                overflow_str = 'under'
            elif mslice.start is not None and mslice.stop is None:
                overflow_str = 'over'
            else:
                overflow_str = 'allnan'
            the_int = the_int.integrate('massreg',mslice,overflow_str)
            if debug:
                print('\t',overflow_str,the_int.values())

        if samplelist is not None:
            the_int = the_int.integrate('sample',samplelist).values(sumw2=True)[()]
        else:
            the_int = the_int.sum('sample').values(sumw2=True)
            if () in the_int:
                the_int = the_int[()]
            else:
                the_int = np.zeros(len(mttbins)-1,dtype=np.float32)

        if debug:
            print('\tdebug',the_int)

        return the_int

    def getQCDfromDataFullMass(inhist,region,default=1.,systematic='nominal',hptslice=slice(h_pt_min,None),mslice=None):
        qcd_data_full = intRegion(inhist['data_obs'],region,systematic=systematic,hptslice=hptslice,mslice=mslice)[0]
        other_mc = intRegion(inhist,region,samplelist=[s.name for s in inhist.identifiers('sample') if s.name!='data_obs' and s.name!='multijet'],systematic=systematic,hptslice=hptslice,mslice=mslice)
        qcd_data = qcd_data_full - other_mc[0]
        #qcd_data = np.array([v if v>0. else default for v in qcd_data])
        qcd_data = np.clip(qcd_data,default,None)
        qcd_data[qcd_data==default] = 0.
        #qcd_data_w2 = np.clip(qcd_data_full-other_mc[1],default,None)
        qcd_data_w2 = np.clip(qcd_data,default,None)
        qcd_data_w2[qcd_data_w2==default] = 0.
        qcd_data_int = hist.poisson_interval(qcd_data,qcd_data_w2)
        qcd_temp = qcd_data
        qcd_temp_dn = np.array([qcd_data_int[0][bi] if qcd_data_int[0][bi] > 0. and not np.isnan(qcd_data_int[0][bi]) else 0. for bi in range(len(qcd_data))])
        qcd_temp_up = np.array([qcd_data_int[1][bi] if qcd_data_int[1][bi] > 0. and not np.isnan(qcd_data_int[1][bi]) else 0. for bi in range(len(qcd_data))])
        return qcd_temp,qcd_temp_dn,qcd_temp_up

    def getQCDfromData(inhist,region,default=1.,systematic='nominal',hptslice=slice(h_pt_min,None),mslice=None):
        tmp = getQCDfromDataFullMass(inhist,region,default,systematic,hptslice,mslice)
        return tmp[0][lowmassbin:highmassbin],tmp[1][lowmassbin:highmassbin],tmp[2][lowmassbin:highmassbin]

    #         [F]    [L]    [P]
    # SR       A      B      C
    # SR(inv)  D      E      F
    # CR       G      H      I
    # CR(inv)  J      K      L
    #
    # shape = D/Ji
    # alt   = G/Ji
    #
    #  A    = D*Gi/Ji
    #  Aalt = G*Di/Ji
    #  B    = D*Hi/Ji
    #  Balt = G*Ei/Ji
    #  C    = D*Ii/Ji
    #  Calt = G*Fi/Ji

    qcd_data_hists = {}
    top_data_hists = {}
    wlnu_data_hists = {}
    sig_data_hists = {}

    qcdfrom_sig = {}
    qcdfrom_top = {}
    qcdfrom_wlnu = {}
    qcdfrom_qcd = {}
    for pfregion in ["fail","loosepass","pass"]:
        qcdfrom_sig[pfregion] = {}
        qcdfrom_top[pfregion] = {}
        qcdfrom_wlnu[pfregion] = {}
        qcdfrom_qcd[pfregion] = {}

    def getQCDRatios(region, category="fail", verbose=False, defval=0.):
        qcd_data_hists[region] = [getQCDfromData(qcd_cr_hists[region],category,hptslice=slice(qcd_cr_hists[region].axis('h_pt').edges()[ix-1] if ix>1 else None,qcd_cr_hists[region].axis('h_pt').edges()[ix] if ix<len(qcd_cr_hists[region].axis('h_pt').edges())-1 else None),default=defval) for ix in range(len(qcd_cr_hists[region].axis('h_pt').edges()))]
        top_data_hists[region] = [getQCDfromData(top_cr_hists[region],category,hptslice=slice(top_cr_hists[region].axis('h_pt').edges()[ix-1] if ix>1 else None,top_cr_hists[region].axis('h_pt').edges()[ix] if ix<len(top_cr_hists[region].axis('h_pt').edges())-1 else None),default=defval) for ix in range(len(top_cr_hists[region].axis('h_pt').edges()))]
        wlnu_data_hists[region] = [getQCDfromData(wlnu_cr_hists[region],category,hptslice=slice(wlnu_cr_hists[region].axis('h_pt').edges()[ix-1] if ix>1 else None,wlnu_cr_hists[region].axis('h_pt').edges()[ix] if ix<len(wlnu_cr_hists[region].axis('h_pt').edges())-1 else None),default=defval) for ix in range(len(wlnu_cr_hists[region].axis('h_pt').edges()))]
        sig_data_hists[region] = [getQCDfromData(sig_hists[region],category,hptslice=slice(sig_hists[region].axis('h_pt').edges()[ix-1] if ix>1 else None,sig_hists[region].axis('h_pt').edges()[ix] if ix<len(sig_hists[region].axis('h_pt').edges())-1 else None),default=defval) for ix in range(len(sig_hists[region].axis('h_pt').edges()))]
        # 0 = nom, 1 = dn, 2 = up

        #qcdratio_sig[region] = [[np.sum(np.stack([(sig_data_hists[region][ix][t]) for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else sig_data_hists[region][t]),axis=0) for ptbin in range(npt)] for t in range(3)]
        #qcdratio_top[region] = [[np.sum(np.stack([(top_data_hists[region][ix][t]) for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else top_data_hists[region][t]),axis=0) for ptbin in range(npt)] for t in range(3)]
        #qcdratio_wlnu[region] = [[np.sum(np.stack([(wlnu_data_hists[region][ix][t]) for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else wlnu_data_hists[region][t]),axis=0) for ptbin in range(npt)] for t in range(3)]
        #qcdratio_qcd[region] = [[np.sum(np.stack([(qcd_data_hists[region][ix][t]) for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else qcd_data_hists[region][t]),axis=0) for ptbin in range(npt)] for t in range(3)] #does this work for len(qcdbins[ptbin])=0?

        qcdfrom_sig[category][region] = [[np.sum(np.stack([(sig_data_hists[region][ix][t]) for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else sig_data_hists[region][t]),axis=0) for ptbin in range(npt)] for t in range(3)]
        qcdfrom_top[category][region] = [[np.sum(np.stack([(top_data_hists[region][ix][t]) for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else top_data_hists[region][t]),axis=0) for ptbin in range(npt)] for t in range(3)]
        qcdfrom_wlnu[category][region] = [[np.sum(np.stack([(wlnu_data_hists[region][ix][t]) for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else wlnu_data_hists[region][t]),axis=0) for ptbin in range(npt)] for t in range(3)]
        qcdfrom_qcd[category][region] = [[np.sum(np.stack([(qcd_data_hists[region][ix][t]) for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else qcd_data_hists[region][t]),axis=0) for ptbin in range(npt)] for t in range(3)] #does this work for len(qcdbins[ptbin])=0?

        if verbose:
            print(region,category,'ratio sig ',qcdfrom_sig[category][region])
            print(region,category,'ratio top ',qcdfrom_top[category][region])
            print(region,category,'ratio wlnu',qcdfrom_wlnu[category][region])
            print(region,category,'ratio qcd ',qcdfrom_qcd[category][region])

    for qcdregion in ["faildphi","lowmet_faildphi","faillep","lowmet_faillep","lowmet","nom"]:
        for pfregion in ["fail","loosepass","pass"]:
            getQCDRatios(qcdregion,category=pfregion,verbose=qcdregion=="faildphi" or qcdregion=="nom")
    qcdfail_temp = {}
    qcdfail_temp_dn = {}
    qcdfail_temp_up = {}

    for region in ['pass','loosepass','fail']:
        theqcdregion = "faillep"
        failtemp = [getQCDfromData(sig_hists[theqcdregion],region,hptslice=slice(sig_hists[theqcdregion].axis('h_pt').edges()[ix-1] if ix>1 else None,sig_hists[theqcdregion].axis('h_pt').edges()[ix] if ix<len(sig_hists[theqcdregion].axis('h_pt').edges())-1 else None),default=0.) for ix in range(len(sig_hists[theqcdregion].axis('h_pt').edges()))]
        qcdfail_temp[region] = [[np.sum(np.stack([(failtemp[ix][t]) for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else failtemp[t]),axis=0) for ptbin in range(npt)] for t in range(3)]
        qcdfail_temp_dn[region] = qcdfail_temp[region][:][1]
        qcdfail_temp_up[region] = qcdfail_temp[region][:][2]
        qcdfail_temp[region] = qcdfail_temp[region][:][0]

    qcdratio_F = {
        'pass':[qcdfail_temp['pass'][ipt]/qcdfail_temp['fail'][ipt] for ipt in range(npt)],
        'loosepass':[qcdfail_temp['loosepass'][ipt]/qcdfail_temp['fail'][ipt] for ipt in range(npt)],
        'fail':[np.ones_like(qcdfail_temp['fail'][ipt]) for ipt in range(npt)],
    }
    qcdratio_F_dn = {
        'pass':[qcdfail_temp_dn['pass'][ipt]/qcdfail_temp_dn['fail'][ipt] for ipt in range(npt)],
        'loosepass':[qcdfail_temp_dn['loosepass'][ipt]/qcdfail_temp_dn['fail'][ipt] for ipt in range(npt)],
        'fail':[np.ones_like(qcdfail_temp['fail'][ipt]) for ipt in range(npt)],
    }
    qcdratio_F_up = {
        'pass':[qcdfail_temp_up['pass'][ipt]/qcdfail_temp_up['fail'][ipt] for ipt in range(npt)],
        'loosepass':[qcdfail_temp_up['loosepass'][ipt]/qcdfail_temp_up['fail'][ipt] for ipt in range(npt)],
        'fail':[np.ones_like(qcdfail_temp['fail'][ipt]) for ipt in range(npt)],
    }
    #for region in ['pass', 'loosepass', 'fail']:
    #    qcdratio_F[region] = np.sum(qcdratio_F[region]*qcd_from_data)/np.sum(qcd_from_data)
    #    qcdratio_F_dn[region] = np.sum(qcdratio_F_dn[region]*qcd_from_data)/np.sum(qcd_from_data)
    #    qcdratio_F_up[region] = np.sum(qcdratio_F_up[region]*qcd_from_data)/np.sum(qcd_from_data)
    print('qcdratio_F',qcdratio_F)
    print('qcdratio_F_dn',qcdratio_F_dn)
    print('qcdratio_F_up',qcdratio_F_up)

    qcdfail_temp = {}
    qcdfail_temp_dn = {}
    qcdfail_temp_up = {}

    for qcdregion in ["faildphi","qcdfaildphi","qcdnom","qcdfaillep","faillep"]:
        qcdfail_temp[qcdregion] = {}
        qcdfail_temp_dn[qcdregion] = {}
        qcdfail_temp_up[qcdregion] = {}
        for region in ['pass','loosepass','fail']:
            #qcdfail_temp[region], qcdfail_temp_dn[region], qcdfail_temp_up[region] = getQCDfromData(qcd_cr_hists["lowmet"],region,default=0.)
            qcdfail_temp[qcdregion][region], qcdfail_temp_dn[qcdregion][region], qcdfail_temp_up[qcdregion][region] = getQCDfromData(sig_hists[qcdregion] if "qcd" not in qcdregion else qcd_cr_hists[qcdregion.replace('qcd','')],region,default=0.)
            #qcdfail_temp[region] = np.array([np.sum(qcdfail_temp[region])])
            #qcdfail_temp_dn[region] = np.array([np.sum(qcdfail_temp_dn[region])])
            #qcdfail_temp_up[region] = np.array([np.sum(qcdfail_temp_up[region])])

    qcdratio_F = {}
    qcdratio_F_dn = {}
    qcdratio_F_up = {}

    for qcdregion in ["faildphi","qcdfaildphi","qcdnom","qcdfaillep","faillep"]:
        qcdratio_F[qcdregion] = {}
        qcdratio_F_dn[qcdregion] = {}
        qcdratio_F_up[qcdregion] = {}
        for ireg,reg in enumerate(["fail","loosepass","pass"]):
            #qcdratio_F[qcdregion][reg] = np.ones_like(qcdfail_temp[qcdregion][reg]) if reg in shaperegion else np.divide(qcdfail_temp[qcdregion][reg],qcdfail_temp[qcdregion][shaperegion[ireg]], out=qcdfail_temp[qcdregion][reg], where=qcdfail_temp[qcdregion][shaperegion[ireg]]!=0.)
            #qcdratio_F_dn[qcdregion][reg] = np.ones_like(qcdfail_temp_dn[qcdregion][reg]) if reg in shaperegion else np.divide(qcdfail_temp_dn[qcdregion][reg],qcdfail_temp_dn[qcdregion][shaperegion[ireg]], out=qcdfail_temp_dn[qcdregion][reg], where=qcdfail_temp_dn[qcdregion][shaperegion[ireg]]!=0.)
            #qcdratio_F_up[qcdregion][reg] = np.ones_like(qcdfail_temp_up[qcdregion][reg]) if reg in shaperegion else np.divide(qcdfail_temp_up[qcdregion][reg],qcdfail_temp_up[qcdregion][shaperegion[ireg]], out=qcdfail_temp_up[qcdregion][reg], where=qcdfail_temp_up[qcdregion][shaperegion[ireg]]!=0.)
            qcdratio_F[qcdregion][reg] = np.ones_like(qcdfail_temp[qcdregion][reg])*(1. if reg in shaperegion else np.sum(qcdfail_temp[qcdregion][reg])/np.sum(qcdfail_temp[qcdregion][shaperegion[ireg]]))
            qcdratio_F_dn[qcdregion][reg] = np.ones_like(qcdfail_temp_dn[qcdregion][reg])*(1. if reg in shaperegion else np.sum(qcdfail_temp_dn[qcdregion][reg])/np.sum(qcdfail_temp_dn[qcdregion][shaperegion[ireg]]))
            qcdratio_F_up[qcdregion][reg] = np.ones_like(qcdfail_temp_up[qcdregion][reg])*(1. if reg in shaperegion else np.sum(qcdfail_temp_up[qcdregion][reg])/np.sum(qcdfail_temp_up[qcdregion][shaperegion[ireg]]))

    for qcdreg in qcdratio_F:
        for reg in qcdratio_F[qcdreg]:
            qcdratio_F[qcdreg][reg] = np.nan_to_num(qcdratio_F[qcdreg][reg])
            qcdratio_F_dn[qcdreg][reg] = np.nan_to_num(qcdratio_F_dn[qcdreg][reg])
            qcdratio_F_up[qcdreg][reg] = np.nan_to_num(qcdratio_F_up[qcdreg][reg])

    qcd_from_data_sig = {reg:{pf:[getQCDfromData(sig_hists[reg],shaperegion[ipf],hptslice=slice(sig_hists[reg].axis('h_pt').edges()[ix-1] if ix>1 else None,sig_hists[reg].axis('h_pt').edges()[ix] if ix<len(sig_hists[reg].axis('h_pt').edges())-1 else None)) for ix in range(len(sig_hists[reg].axis('h_pt').edges()))] for ipf,pf in enumerate(['fail', 'loosepass', 'pass'])} for reg in ["nom","faildphi","lowmet_faildphi","faillep","lowmet_faillep","lowmet"]}
    qcd_from_data_qcd = {reg:{pf:[getQCDfromData(qcd_cr_hists[reg],shaperegion[ipf],hptslice=slice(qcd_cr_hists[reg].axis('h_pt').edges()[ix-1] if ix>1 else None,qcd_cr_hists[reg].axis('h_pt').edges()[ix] if ix<len(qcd_cr_hists[reg].axis('h_pt').edges())-1 else None)) for ix in range(len(qcd_cr_hists[reg].axis('h_pt').edges()))] for ipf,pf in enumerate(['fail', 'loosepass', 'pass'])} for reg in ["nom","faildphi","lowmet_faildphi","faillep","lowmet_faillep","lowmet"]}
    qcd_from_data_top = {reg:{pf:[getQCDfromData(top_cr_hists[reg],shaperegion[ipf],hptslice=slice(top_cr_hists[reg].axis('h_pt').edges()[ix-1] if ix>1 else None,top_cr_hists[reg].axis('h_pt').edges()[ix] if ix<len(top_cr_hists[reg].axis('h_pt').edges())-1 else None)) for ix in range(len(top_cr_hists[reg].axis('h_pt').edges()))] for ipf,pf in enumerate(['fail', 'loosepass', 'pass'])} for reg in ["nom","faildphi","lowmet_faildphi","faillep","lowmet_faillep","lowmet"]}
    qcd_from_data_wlnu = {reg:{pf:[getQCDfromData(wlnu_cr_hists[reg],shaperegion[ipf],hptslice=slice(wlnu_cr_hists[reg].axis('h_pt').edges()[ix-1] if ix>1 else None,wlnu_cr_hists[reg].axis('h_pt').edges()[ix] if ix<len(wlnu_cr_hists[reg].axis('h_pt').edges())-1 else None)) for ix in range(len(wlnu_cr_hists[reg].axis('h_pt').edges()))] for ipf,pf in enumerate(['fail', 'loosepass', 'pass'])} for reg in ["nom","faildphi","lowmet_faildphi","faillep","lowmet_faillep","lowmet"]}
    
    print('WLNU QCD FAILLEP',qcd_from_data_wlnu['faillep'])
    print(' QCD QCD FAILLEP',qcd_from_data_qcd['faillep'])

    qcdratio_sig_dn = {}
    qcdratio_sig_up = {}
    qcdratio_sig = {}

    qcdratio_qcd_dn = {}
    qcdratio_qcd_up = {}
    qcdratio_qcd = {}

    qcdratio_top_dn = {}
    qcdratio_top_up = {}
    qcdratio_top = {}

    qcdratio_wlnu_dn = {}
    qcdratio_wlnu_up = {}
    qcdratio_wlnu = {}

    for pfreg in ["fail","loosepass","pass"]:
        qcdratio_sig_dn[pfreg] = {}
        qcdratio_sig_up[pfreg] = {}
        qcdratio_sig[pfreg] = {}
        qcdratio_qcd_dn[pfreg] = {}
        qcdratio_qcd_up[pfreg] = {}
        qcdratio_qcd[pfreg] = {}
        qcdratio_top_dn[pfreg] = {}
        qcdratio_top_up[pfreg] = {}
        qcdratio_top[pfreg] = {}
        qcdratio_wlnu_dn[pfreg] = {}
        qcdratio_wlnu_up[pfreg] = {}
        qcdratio_wlnu[pfreg] = {}
        for qcdreg in ["faildphi","lowmet_faildphi","faillep","lowmet_faillep","lowmet","nom"]:
            qcdratio_sig_dn[pfreg][qcdreg] = []
            qcdratio_sig_up[pfreg][qcdreg] = []
            qcdratio_sig[pfreg][qcdreg] = []
            qcdratio_qcd_dn[pfreg][qcdreg] = []
            qcdratio_qcd_up[pfreg][qcdreg] = []
            qcdratio_qcd[pfreg][qcdreg] = []
            qcdratio_top_dn[pfreg][qcdreg] = []
            qcdratio_top_up[pfreg][qcdreg] = []
            qcdratio_top[pfreg][qcdreg] = []
            qcdratio_wlnu_dn[pfreg][qcdreg] = []
            qcdratio_wlnu_up[pfreg][qcdreg] = []
            qcdratio_wlnu[pfreg][qcdreg] = []
            for ptbin in range(npt):
                marr_dn = []
                marr_up = []
                marr = []
                for ix in range(len(qcdfrom_qcd[pfreg][qcdreg][0][ptbin])):
                    if qcdfrom_qcd[pfreg][qcdreg][0][ptbin][ix]>0.:
                        marr_dn.append(qcdfrom_qcd[pfreg]["nom"][1][ptbin][ix]/qcdfrom_qcd[pfreg][qcdreg][2][ptbin][ix])
                        marr_up.append(qcdfrom_qcd[pfreg]["nom"][2][ptbin][ix]/qcdfrom_qcd[pfreg][qcdreg][1][ptbin][ix])
                        marr.append(qcdfrom_qcd[pfreg]["nom"][0][ptbin][ix]/qcdfrom_qcd[pfreg][qcdreg][0][ptbin][ix])
                    else:
                        if ix==0:
                            marr_dn.append(1.)
                            marr_up.append(1.)
                            marr.append(1.)
                        else:
                            marr_dn.append(marr_dn[-1])
                            marr_up.append(marr_up[-1])
                            marr.append(marr[-1])
                qcdratio_sig_dn[pfreg][qcdreg].append(np.array(marr_dn))
                qcdratio_sig_up[pfreg][qcdreg].append(np.array(marr_up))
                qcdratio_sig[pfreg][qcdreg].append(np.array(marr))

                marr_dn = []
                marr_up = []
                marr = []
                for ix in range(len(qcdfrom_sig[pfreg][qcdreg][0][ptbin])):
                    if qcdfrom_qcd[pfreg][qcdreg][0][ptbin][ix]>0.:
                        marr_dn.append(qcdfrom_sig[pfreg][qcdreg][1][ptbin][ix]/qcdfrom_qcd[pfreg][qcdreg][2][ptbin][ix])
                        marr_up.append(qcdfrom_sig[pfreg][qcdreg][2][ptbin][ix]/qcdfrom_qcd[pfreg][qcdreg][1][ptbin][ix])
                        marr.append(qcdfrom_sig[pfreg][qcdreg][0][ptbin][ix]/qcdfrom_qcd[pfreg][qcdreg][0][ptbin][ix])
                    else:
                        if ix==0:
                            marr_dn.append(1.)
                            marr_up.append(1.)
                            marr.append(1.)
                        else:
                            marr_dn.append(marr_dn[-1])
                            marr_up.append(marr_up[-1])
                            marr.append(marr[-1])
                qcdratio_qcd_dn[pfreg][qcdreg].append(np.array(marr_dn))
                qcdratio_qcd_up[pfreg][qcdreg].append(np.array(marr_up))
                qcdratio_qcd[pfreg][qcdreg].append(np.array(marr))

                marr_dn = []
                marr_up = []
                marr = []
                for ix in range(len(qcdfrom_top[pfreg][qcdreg][0][ptbin])):
                    if qcdfrom_qcd[pfreg][qcdreg][0][ptbin][ix]>0.:
                        marr_dn.append(qcdfrom_top[pfreg][qcdreg][1][ptbin][ix]/qcdfrom_qcd[pfreg][qcdreg][2][ptbin][ix])
                        marr_up.append(qcdfrom_top[pfreg][qcdreg][2][ptbin][ix]/qcdfrom_qcd[pfreg][qcdreg][1][ptbin][ix])
                        marr.append(qcdfrom_top[pfreg][qcdreg][0][ptbin][ix]/qcdfrom_qcd[pfreg][qcdreg][0][ptbin][ix])
                    else:
                        if ix==0:
                            marr_dn.append(1.)
                            marr_up.append(1.)
                            marr.append(1.)
                        else:
                            marr_dn.append(marr_dn[-1])
                            marr_up.append(marr_up[-1])
                            marr.append(marr[-1])
                qcdratio_top_dn[pfreg][qcdreg].append(np.array(marr_dn))
                qcdratio_top_up[pfreg][qcdreg].append(np.array(marr_up))
                qcdratio_top[pfreg][qcdreg].append(np.array(marr))

                marr_dn = []
                marr_up = []
                marr = []
                for ix in range(len(qcdfrom_wlnu[pfreg][qcdreg][0][ptbin])):
                    if qcdfrom_qcd[pfreg][qcdreg][0][ptbin][ix]>0.:
                        marr_dn.append(qcdfrom_wlnu[pfreg][qcdreg][1][ptbin][ix]/qcdfrom_qcd[pfreg][qcdreg][2][ptbin][ix])
                        marr_up.append(qcdfrom_wlnu[pfreg][qcdreg][2][ptbin][ix]/qcdfrom_qcd[pfreg][qcdreg][1][ptbin][ix])
                        marr.append(qcdfrom_wlnu[pfreg][qcdreg][0][ptbin][ix]/qcdfrom_qcd[pfreg][qcdreg][0][ptbin][ix])
                    else:
                        if ix==0:
                            marr_dn.append(1.)
                            marr_up.append(1.)
                            marr.append(1.)
                        else:
                            marr_dn.append(marr_dn[-1])
                            marr_up.append(marr_up[-1])
                            marr.append(marr[-1])
                qcdratio_wlnu_dn[pfreg][qcdreg].append(np.array(marr_dn))
                qcdratio_wlnu_up[pfreg][qcdreg].append(np.array(marr_up))
                qcdratio_wlnu[pfreg][qcdreg].append(np.array(marr))

    for pfreg in ["fail","loosepass","pass"]:
        for qcdreg in ["faildphi","lowmet_faildphi","faillep","lowmet_faillep","lowmet","nom"]:
            for ptbin in range(npt):
                qcdratio_sig[pfreg][qcdreg][ptbin] = np.sum(qcdratio_sig[pfreg][qcdreg][ptbin]*qcdfrom_sig[pfreg][qcdreg][0][ptbin])/np.sum(qcdfrom_sig[pfreg][qcdreg][0][ptbin])
                qcdratio_qcd[pfreg][qcdreg][ptbin] = np.sum(qcdratio_qcd[pfreg][qcdreg][ptbin]*qcdfrom_qcd[pfreg][qcdreg][0][ptbin])/np.sum(qcdfrom_qcd[pfreg][qcdreg][0][ptbin])
                qcdratio_top[pfreg][qcdreg][ptbin] = np.sum(qcdratio_top[pfreg][qcdreg][ptbin]*qcdfrom_top[pfreg][qcdreg][0][ptbin])/np.sum(qcdfrom_top[pfreg][qcdreg][0][ptbin])
                qcdratio_wlnu[pfreg][qcdreg][ptbin] = np.sum(qcdratio_wlnu[pfreg][qcdreg][ptbin]*qcdfrom_wlnu[pfreg][qcdreg][0][ptbin])/np.sum(qcdfrom_wlnu[pfreg][qcdreg][0][ptbin])
                qcdratio_sig_dn[pfreg][qcdreg][ptbin] = np.sum(qcdratio_sig_dn[pfreg][qcdreg][ptbin]*qcdfrom_sig[pfreg][qcdreg][1][ptbin])/np.sum(qcdfrom_sig[pfreg][qcdreg][1][ptbin])
                qcdratio_qcd_dn[pfreg][qcdreg][ptbin] = np.sum(qcdratio_qcd_dn[pfreg][qcdreg][ptbin]*qcdfrom_qcd[pfreg][qcdreg][1][ptbin])/np.sum(qcdfrom_qcd[pfreg][qcdreg][1][ptbin])
                qcdratio_top_dn[pfreg][qcdreg][ptbin] = np.sum(qcdratio_top_dn[pfreg][qcdreg][ptbin]*qcdfrom_top[pfreg][qcdreg][1][ptbin])/np.sum(qcdfrom_top[pfreg][qcdreg][1][ptbin])
                qcdratio_wlnu_dn[pfreg][qcdreg][ptbin] = np.sum(qcdratio_wlnu_dn[pfreg][qcdreg][ptbin]*qcdfrom_wlnu[pfreg][qcdreg][1][ptbin])/np.sum(qcdfrom_wlnu[pfreg][qcdreg][1][ptbin])
                qcdratio_sig_up[pfreg][qcdreg][ptbin] = np.sum(qcdratio_sig_up[pfreg][qcdreg][ptbin]*qcdfrom_sig[pfreg][qcdreg][2][ptbin])/np.sum(qcdfrom_sig[pfreg][qcdreg][2][ptbin])
                qcdratio_qcd_up[pfreg][qcdreg][ptbin] = np.sum(qcdratio_qcd_up[pfreg][qcdreg][ptbin]*qcdfrom_qcd[pfreg][qcdreg][2][ptbin])/np.sum(qcdfrom_qcd[pfreg][qcdreg][2][ptbin])
                qcdratio_top_up[pfreg][qcdreg][ptbin] = np.sum(qcdratio_top_up[pfreg][qcdreg][ptbin]*qcdfrom_top[pfreg][qcdreg][2][ptbin])/np.sum(qcdfrom_top[pfreg][qcdreg][2][ptbin])
                qcdratio_wlnu_up[pfreg][qcdreg][ptbin] = np.sum(qcdratio_wlnu_up[pfreg][qcdreg][ptbin]*qcdfrom_wlnu[pfreg][qcdreg][2][ptbin])/np.sum(qcdfrom_wlnu[pfreg][qcdreg][2][ptbin])


    for ptbin in range(npt):
        for iregion,region in enumerate(['fail','loosepass','pass']):
            ch = rl.Channel("ptbin%d%s%s%s" % (ptbin, region, 'hadhad', year))
            model.addChannel(ch)

            isPass = region=='pass'
            isLoosePass = region=='loosepass'

            thehist = sig_hists["nom"]
    
            for sName in thehist.identifiers('sample'):
                if sName.name=='ignore':
                    continue
                if sName.name=='data_obs' and isPass and not usingData:
                    continue
                tempint = intRegion(thehist[sName],region,hptslice=ptslices[ptbin])
                templ = (tempint[0][lowmassbin:highmassbin], mtt.binning, mtt.name, tempint[1][lowmassbin:highmassbin])
                #print('\tDEBUG',region,ptbin,sName.name,templ[0])
                for iv,val in enumerate(templ[0]):
                    if val<0.:
                        templ[0][iv] = 0.
                        #print('\t\t',templ[0])
                if sName.name=='data_obs':
                    ch.setObservation(templ, read_sumw2=True)
                else:
                    if sName.name=='multijet':
                        #qcdpred = np.stack([(qcd_from_data_sig["faildphi"][shaperegion[iregion]][ix][0]*qcdratio_sig[shaperegion[iregion]]["faildphi"][ptbin]*qcdratio_F["qcdnom"][region][ix]+qcd_from_data_qcd["nom"][shaperegion[iregion]][ix][0]*qcdratio_qcd[shaperegion[iregion]]["nom"][ptbin]*qcdratio_F["faildphi"][region][ix])/2. for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else [np.zeros_like(qcd_from_data_sig["faildphi"][region][0][0])])
                        #qcdpred_dn = np.clip(np.sum(np.stack([qcd_from_data_sig["faildphi"][shaperegion[iregion]][ix][1] * qcdratio_sig[shaperegion[iregion]]["faildphi"][ptbin] * qcdratio_F["qcdnom"][region][ix] for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else [np.zeros_like(qcd_from_data_sig["faildphi"][region][0][0])]),axis=0),0.,None)
                        #qcdpred_up = np.sum(np.stack([qcd_from_data_qcd["nom"][shaperegion[iregion]][ix][2] * qcdratio_qcd[shaperegion[iregion]]["nom"][ptbin] * qcdratio_F["faildphi"][region][ix] for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else [np.zeros_like(qcd_from_data_qcd["nom"][region][0][0])]),axis=0)
                        #qcdpred = np.stack([(qcd_from_data_sig["faillep"][shaperegion[iregion]][ix][0]*qcdratio_sig[shaperegion[iregion]]["faillep"][ptbin]*qcdratio_F["qcdnom"][region][ix]+qcd_from_data_qcd["faildphi"][shaperegion[iregion]][ix][0]*qcdratio_sig[shaperegion[iregion]]["faildphi"][ptbin]*qcdratio_F["qcdnom"][region][ix])/2. for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else [np.zeros_like(qcd_from_data_sig["faillep"][region][0][0])])
                        #qcdpred_dn = np.clip(np.sum(np.stack([qcd_from_data_sig["faillep"][shaperegion[iregion]][ix][1] * qcdratio_sig[shaperegion[iregion]]["faillep"][ptbin] * qcdratio_F["qcdnom"][region][ix] for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else [np.zeros_like(qcd_from_data_sig["faillep"][region][0][0])]),axis=0),0.,None)
                        #qcdpred_up = np.sum(np.stack([qcd_from_data_sig["faildphi"][shaperegion[iregion]][ix][2] * qcdratio_sig[shaperegion[iregion]]["faildphi"][ptbin] * qcdratio_F["qcdnom"][region][ix] for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else [np.zeros_like(qcd_from_data_sig["faildphi"][region][0][0])]),axis=0)
                        qcdpred = np.stack([qcd_from_data_qcd["nom"][shaperegion[iregion]][ix][0]*qcdratio_qcd[shaperegion[iregion]]["nom"][ptbin]*qcdratio_F["qcdnom"][region][ix] for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else [np.zeros_like(qcd_from_data_sig["nom"][region][0][0])])
                        qcdpred_dn = np.clip(np.sum(np.stack([qcd_from_data_qcd["nom"][shaperegion[iregion]][ix][1] * qcdratio_qcd[shaperegion[iregion]]["nom"][ptbin] * qcdratio_F["qcdnom"][region][ix] for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else [np.zeros_like(qcd_from_data_sig["nom"][region][0][0])]),axis=0),0.,None)
                        qcdpred_up = np.clip(np.sum(np.stack([qcd_from_data_qcd["nom"][shaperegion[iregion]][ix][2] * qcdratio_qcd[shaperegion[iregion]]["nom"][ptbin] * qcdratio_F["qcdnom"][region][ix] for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else [np.zeros_like(qcd_from_data_sig["nom"][region][0][0])]),axis=0),0.,None)
                        print('qcdpred',qcdpred)
                        qcdpred = np.sum(qcdpred,axis=0)

                        templ = (qcdpred, mtt.binning, mtt.name, np.ones_like(qcdpred)*1000000.)
                    if doHtt:
                        stype = rl.Sample.SIGNAL if 'htt125' in sName.name else rl.Sample.BACKGROUND
                    else:
                        stype = rl.Sample.SIGNAL if 'phi' in sName.name else rl.Sample.BACKGROUND
                    sample = rl.TemplateSample(ch.name + '_' + sName.name, stype, templ)
        
                    #sample.setParamEffect(jec, 1.05)
                    if sName.name in syst_dict:
                        nom = intRegion(thehist[sName],region,hptslice=ptslices[ptbin])[0][lowmassbin:highmassbin]
                        for syst in syst_dict[sName.name]:
                            syst_params = syst_dict[sName.name][syst]
                            dnmod = intRegion(thehist[sName],region,syst_params[1],hptslice=ptslices[ptbin])[0][lowmassbin:highmassbin]
                            upmod = intRegion(thehist[sName],region,syst_params[2],hptslice=ptslices[ptbin])[0][lowmassbin:highmassbin]
                            sample.setParamEffect(syst_params[0], np.divide(upmod, nom, out=np.ones_like(nom), where=nom>0.) if (upmod!=nom).all() else np.ones_like(nom)*1.001, np.divide(dnmod, nom, out=np.ones_like(nom), where=nom>0.) if (dnmod!=nom).all() else np.ones_like(nom)*0.999)

                    if sName.name=='top':
                        sample.setParamEffect(top_norm, 1.05)
                    if sName.name=='wlnu':
                        sample.setParamEffect(wlnu_norm, 1.10)
                    if sName.name=='multijet':
                        sample.setParamEffect(qcd_norm, 1.50)
                        qcd_shape_dn = np.divide(qcdpred_dn, qcdpred, out=np.ones_like(qcdpred), where=qcdpred>0.)
                        #qcd_shape_dn = qcd_shape_dn*np.sum(qcdpred)/np.sum(qcdpred_dn)
                        qcd_shape_up = np.divide(qcdpred_up, qcdpred, out=np.ones_like(qcdpred), where=qcdpred>0.)
                        #qcd_shape_up = qcd_shape_up*np.sum(qcdpred)/np.sum(qcdpred_up)
                        sample.setParamEffect(qcd_pass if isPass else qcd_loosepass if isLoosePass else qcd_fail, qcd_shape_dn, qcd_shape_up)
                    if sName.name=='vvqq':
                        sample.setParamEffect(vvqq_norm, 1.05)
                    if sName.name=='dy':
                        sample.setParamEffect(dy_norm, 1.05)
                    if sName.name=='htt125' or sName.name=='dy' or 'phi' in sName.name:
                        nom_full = intRegion(thehist[sName],region,hptslice=ptslices[ptbin])[0]
                        shift_dn = nom_full[lowmassbin+1:highmassbin+1 if highmassbin>0 else None]
                        shift_up = nom_full[lowmassbin-1:highmassbin-1]
                        #sample.setParamEffect(m_scale, np.divide(shift_dn, nom_full[lowmassbin:highmassbin], out=np.ones_like(nom_full[lowmassbin:highmassbin]), where=nom_full[lowmassbin:highmassbin]>0.), np.divide(shift_up, nom_full[lowmassbin:highmassbin], out=np.ones_like(nom_full[lowmassbin:highmassbin]), where=nom_full[lowmassbin:highmassbin]>0.))
                        #print('had htt125 ptbin',ptbin,region,nom_full[lowmassbin:highmassbin])
                    sample.setParamEffect(trig, 1.02)
                    sample.setParamEffect(lumi, 1.024)
            
                    ch.addSample(sample)
            if not usingData and isPass:
                ch.setObservation((np.zeros(len(mtt.binning)-1),mtt.binning, mtt.name))
    
 
    print('htt125')
    h125vals = [intRegion(thehist['htt125'],'pass',hptslice=ptslices[ptbin])[lowmassbin:highmassbin] for ptbin in range(npt)]
    if not doHtt:
        print('phitt50')
        phi50vals = [intRegion(thehist['phitt50'],'pass',hptslice=ptslices[ptbin])[lowmassbin:highmassbin] for ptbin in range(npt)]
    print('dy')
    dyvals = [intRegion(thehist['dy'],'pass',hptslice=ptslices[ptbin])[lowmassbin:highmassbin] for ptbin in range(npt)]
    print('wlnu')
    wlnuvals = [intRegion(thehist['wlnu'],'pass',hptslice=ptslices[ptbin])[lowmassbin:highmassbin] for ptbin in range(npt)]
    if doHtt:
        probs = [np.prod([math.erfc(h125vals[ipt][ib]/(math.sqrt(abs(dyvals[ipt][ib]+wlnuvals[ipt][ib])) if dyvals[ipt][ib]+wlnuvals[ipt][ib]>2. else 2.)) for ib in range(len(h125vals[ipt]))]) for ipt in range(npt)]
    else:
        probs = [np.prod([math.erfc(phi50vals[ipt][ib]/(math.sqrt(abs(dyvals[ipt][ib]+wlnuvals[ipt][ib]+h125vals[ipt][ib])) if dyvals[ipt][ib]+wlnuvals[ipt][ib]+h125vals[ipt][ib]>2. else 2.)) for ib in range(len(h125vals[ipt]))]) for ipt in range(npt)]
    print('had','PROB',np.prod(probs),'->',erfinv(1.-np.prod(probs)),'\t\t',probs)
    

    for ptbin in range(npt):
        failCh = model['ptbin%dfail%s%s' % (ptbin,'hadhad',year)]
        passCh = model['ptbin%dpass%s%s' % (ptbin,'hadhad',year)]
        loosePassCh = model['ptbin%dloosepass%s%s' % (ptbin,'hadhad',year)]

        qcdpass = passCh['multijet']
        qcdloosepass = loosePassCh['multijet']
        qcdfail = failCh['multijet']
        qcdpass.setParamEffect(qcdnormSF, 1.20)
        qcdloosepass.setParamEffect(qcdnormSF, 1.20)
        qcdfail.setParamEffect(qcdnormSF, 1.20)
        #qcdpass.setParamEffect(qcdnormSF, 1*qcdnormSF)
        #qcdloosepass.setParamEffect(qcdnormSF, 1*qcdnormSF)
        #qcdfail.setParamEffect(qcdnormSF, 1*qcdnormSF)

        toppass = passCh['top']
        toploosepass = loosePassCh['top']
        topfail = failCh['top']
        topPF = (toploosepass.getExpectation(nominal=True).sum() + toppass.getExpectation(nominal=True).sum()) / topfail.getExpectation(nominal=True).sum()
        topLPF = toppass.getExpectation(nominal=True).sum() / toploosepass.getExpectation(nominal=True).sum()
        topRLPF = 1./(1.+(1./topLPF)) # P/(L+P)
        toppass.setParamEffect(topLeffSF, 1*topLeffSF)
        toploosepass.setParamEffect(topLeffSF, (1 - topLeffSF) * topLPF + 1)
        toploosepass.setParamEffect(topeffSF, 1*topeffSF)
        topfail.setParamEffect(topeffSF, (1 - topeffSF) * (1-topRLPF) * topPF + 1)
        #toppass.setParamEffect(topnormSF, 1*topnormSF)
        #toploosepass.setParamEffect(topnormSF, 1*topnormSF)
        #topfail.setParamEffect(topnormSF, 1*topnormSF)

        wlnupass = passCh['wlnu']
        wlnuloosepass = loosePassCh['wlnu']
        wlnufail = failCh['wlnu']
        wlnuPF = (wlnuloosepass.getExpectation(nominal=True).sum() + wlnupass.getExpectation(nominal=True).sum()) / wlnufail.getExpectation(nominal=True).sum()
        wlnuLPF = wlnupass.getExpectation(nominal=True).sum() / wlnuloosepass.getExpectation(nominal=True).sum()
        wlnuRLPF = 1./(1.+(1./wlnuLPF)) # P/(L+P)
        wlnupass.setParamEffect(wlnuLeffSF, 1*wlnuLeffSF)
        wlnuloosepass.setParamEffect(wlnuLeffSF, (1 - wlnuLeffSF) * wlnuLPF + 1)
        wlnuloosepass.setParamEffect(wlnueffSF, 1*wlnueffSF)
        wlnufail.setParamEffect(wlnueffSF, (1 - wlnueffSF) * (1-wlnuRLPF) * wlnuPF + 1)
        #wlnupass.setParamEffect(wlnunormSF, 1*wlnunormSF)
        #wlnuloosepass.setParamEffect(wlnunormSF, 1*wlnunormSF)
        #wlnufail.setParamEffect(wlnunormSF, 1*wlnunormSF)

        dypass = passCh['dy']
        dyloosepass = loosePassCh['dy']
        dyfail = failCh['dy']
        dyLP = dypass.getExpectation(nominal=True).sum() / dyloosepass.getExpectation(nominal=True).sum()
        httLP = passCh['htt125'].getExpectation(nominal=True).sum() / loosePassCh['htt125'].getExpectation(nominal=True).sum()
        if not doHtt:
            phittLP = {m:passCh['phitt%s'%m].getExpectation(nominal=True).sum() / loosePassCh['phitt%s'%m].getExpectation(nominal=True).sum() if loosePassCh['phitt%s'%m].getExpectation(nominal=True).sum()>0. else 1. for m in masspoints}
        dypass.setParamEffect(dy_eff, 1*dy_eff)
        dyloosepass.setParamEffect(dy_eff, (1 - dy_eff) * dyLP + 1)
        passCh['htt125'].setParamEffect(dy_eff, 1*dy_eff)
        loosePassCh['htt125'].setParamEffect(dy_eff, (1 - dy_eff) * httLP + 1)
        if not doHtt:
            for m in masspoints:
                passCh['phitt%s'%m].setParamEffect(dy_eff, 1*dy_eff)
                loosePassCh['phitt%s'%m].setParamEffect(dy_eff, (1 - dy_eff) * phittLP[m] + 1)
        dypass.setParamEffect(rdy, 1*rdy)
        dyloosepass.setParamEffect(rdy, 1*rdy)
        dyfail.setParamEffect(rdy, 1*rdy)
        #dypass.setParamEffect(rdy, 1.05)
        #dyloosepass.setParamEffect(rdy, 1.05)
        #dyfail.setParamEffect(rdy, 1.05)
        if not doHtt:
            passCh['htt125'].setParamEffect(rh125, 1.10)
            loosePassCh['htt125'].setParamEffect(rh125, 1.10)
            failCh['htt125'].setParamEffect(rh125, 1.10)


    # Fill in top CR
    for region in ['pass', 'loosepass', 'fail']:
        ch = rl.Channel("topCR%s%s%s" % (region, 'hadhad', year))
        model.addChannel(ch)

        isPass = region=='pass'
        isLoosePass = region=='loosepass'

        thehist = top_cr_hists["nom"]

        for sName in thehist.identifiers('sample'):
            if sName.name=='ignore':
                continue
            tempint = intRegion(thehist[sName],region)
            if singleBinCR:
                templ = (np.array([np.sum(tempint[0][lowmassbin:highmassbin])]), mttone.binning, mttone.name, np.array([np.sum(tempint[1][lowmassbin:highmassbin])]))
                if includeLowMass:
                    templ = (np.array([np.sum(intRegion(thehist[sName],region,mrebin=None))]), mttone.binning, mttone.name)
            else:
                templ = (tempint[0][lowmassbin:highmassbin], mtt.binning, mtt.name, tempint[1][lowmassbin:highmassbin])
            for iv,val in enumerate(templ[0]):
                if val<0.:
                    templ[0][iv] = 0.
            if sName.name=='data_obs':
                ch.setObservation(templ, read_sumw2=True)
            else:
                if sName.name=='multijet':
                    #qcdpred = np.sum(np.stack([(qcd_from_data_top["faildphi"][shaperegion[iregion]][ix][0]*qcdratio_top[shaperegion[iregion]]["faildphi"][ptbin]*qcdratio_F["qcdnom"][region][ix]+qcd_from_data_top["nom"][shaperegion[iregion]][ix][0]*qcdratio_top[shaperegion[iregion]]["nom"][ptbin]*qcdratio_F["faildphi"][region][ix])/2. for ix in range(len(qcd_from_data_top["faildphi"][region]))]),axis=0)
                    #qcdpred_dn = np.clip(np.sum(np.stack([qcd_from_data_top["faildphi"][shaperegion[iregion]][ix][1] * qcdratio_top[shaperegion[iregion]]["faildphi"][ptbin] * qcdratio_F["qcdnom"][region][ix] for ix in range(len(qcd_from_data_top["faildphi"][region]))]),axis=0),0.,None)
                    #qcdpred_up = np.sum(np.stack([qcd_from_data_qcd["nom"][shaperegion[iregion]][ix][2] * qcdratio_qcd[shaperegion[iregion]]["nom"][ptbin] * qcdratio_F["faildphi"][region][ix] for ix in range(len(qcd_from_data_qcd["nom"][region]))]),axis=0)
                    qcdpred = np.sum(np.stack([qcd_from_data_top["faillep"][shaperegion[iregion]][ix][0]*qcdratio_top[shaperegion[iregion]]["faillep"][ptbin]*qcdratio_F["qcdnom"][region][ix] for ix in range(len(qcd_from_data_top["faillep"][region]))]),axis=0)
                    qcdpred_dn = np.clip(np.sum(np.stack([qcd_from_data_top["faillep"][shaperegion[iregion]][ix][1] * qcdratio_top[shaperegion[iregion]]["faillep"][ptbin] * qcdratio_F["qcdnom"][region][ix] for ix in range(len(qcd_from_data_top["faillep"][region]))]),axis=0),0.,None)
                    qcdpred_up = np.clip(np.sum(np.stack([qcd_from_data_top["faillep"][shaperegion[iregion]][ix][2] * qcdratio_top[shaperegion[iregion]]["faillep"][ptbin] * qcdratio_F["qcdnom"][region][ix] for ix in range(len(qcd_from_data_top["faillep"][region]))]),axis=0),0.,None)

                    if singleBinCR:
                        if includeLowMass:
                            qcdpred = np.sum(np.array(qcd_from_data_full_top[region]), axis=0)
                            qcdpred_dn = np.sum(np.array(qcd_from_data_full_top[region]) - np.array(qcd_from_data_full_top_err[region]), axis=0)
                            qcdpred_up = np.sum(np.array(qcd_from_data_full_top[region]) + np.array(qcd_from_data_full_top_err[region]), axis=0)
                        templ = (np.array([np.sum(qcdpred)]), mttone.binning, mttone.name, np.array([1000000.]))
                    else:
                        templ = (qcdpred, mtt.binning, mtt.name, np.ones_like(qcdpred)*1000000.)
                if doHtt:
                    stype = rl.Sample.SIGNAL if 'htt125' in sName.name else rl.Sample.BACKGROUND
                else:
                    stype = rl.Sample.SIGNAL if 'phi' in sName.name else rl.Sample.BACKGROUND
                sample = rl.TemplateSample(ch.name + '_' + sName.name, stype, templ)
    
                #sample.setParamEffect(jec, 1.05)
                if sName.name in syst_dict:
                    nom_base = intRegion(thehist[sName],region)[0][lowmassbin:highmassbin]
                    if singleBinCR:
                        if includeLowMass:
                            nom_base = intRegion(thehist[sName],region,mrebin=None)
                    for syst in syst_dict[sName.name]:
                        syst_params = syst_dict[sName.name][syst]
                        dnmod = intRegion(thehist[sName],region,systematic=syst_params[1])[0][lowmassbin:highmassbin]
                        upmod = intRegion(thehist[sName],region,systematic=syst_params[2])[0][lowmassbin:highmassbin]
                        if singleBinCR:
                            if includeLowMass:
                                dnmod = intRegion(thehist[sName],region,mrebin=None,systematic=syst_params[1])
                                upmod = intRegion(thehist[sName],region,mrebin=None,systematic=syst_params[2])
                            dnmod = np.array([np.sum(dnmod)])
                            nom = np.array([np.sum(nom_base)])
                            upmod = np.array([np.sum(upmod)])
                        else:
                            nom = nom_base
                        sample.setParamEffect(syst_params[0], np.divide(upmod, nom, out=np.ones_like(nom), where=nom>0.) if (upmod!=nom).all() else np.ones_like(nom)*1.001, np.divide(dnmod, nom, out=np.ones_like(nom), where=nom>0.) if (dnmod!=nom).all() else np.ones_like(nom)*0.999)

                if sName.name=='top':
                    sample.setParamEffect(top_norm, 1.05)
                if sName.name=='wlnu':
                    sample.setParamEffect(wlnu_norm, 1.10)
                if sName.name=='multijet':
                    sample.setParamEffect(qcd_norm, 1.20)
                    if singleBinCR:
                        err_dn = np.sqrt(np.sum((qcdpred - qcdpred_dn)*(qcdpred - qcdpred_dn)))
                        err_up = np.sqrt(np.sum((qcdpred_up - qcdpred)*(qcdpred_up - qcdpred)))
                        qcdpred = np.array([np.sum(qcdpred)])
                        qcdpred_dn = qcdpred - err_dn
                        qcdpred_up = qcdpred + err_up
                    sample.setParamEffect(qcdtop_pass if isPass else qcdtop_loosepass if isLoosePass else qcdtop_fail, np.divide(qcdpred_dn, qcdpred, out=np.ones_like(qcdpred), where=qcdpred>0.), np.divide(qcdpred_up, qcdpred, out=np.ones_like(qcdpred), where=qcdpred>0.))
                if sName.name=='vvqq':
                    sample.setParamEffect(vvqq_norm, 1.05)
                if sName.name=='dy':
                    sample.setParamEffect(dy_norm, 1.05)
                if sName.name=='htt125' or sName.name=='dy' or 'phi' in sName.name:
                    nom_full = intRegion(thehist[sName],region)
                    shift_dn = nom_full[lowmassbin+1:highmassbin+1 if highmassbin>0 else None]
                    shift_up = nom_full[lowmassbin-1:highmassbin-1]
                    #sample.setParamEffect(m_scale, np.divide(shift_dn, nom_full[lowmassbin:highmassbin], out=np.ones_like(nom_full[lowmassbin:highmassbin]), where=nom_full[lowmassbin:highmassbin]>0.), np.divide(shift_up, nom_full[lowmassbin:highmassbin], out=np.ones_like(nom_full[lowmassbin:highmassbin]), where=nom_full[lowmassbin:highmassbin]>0.))
                sample.setParamEffect(trig, 1.02)
                sample.setParamEffect(lumi, 1.024)
        
                ch.addSample(sample)

    qcdpass = model['topCRpasshadhad%s'%year]['multijet']
    qcdloosepass = model['topCRloosepasshadhad%s'%year]['multijet']
    qcdfail = model['topCRfailhadhad%s'%year]['multijet']
    #qcdpass.setParamEffect(qcdnormSF_top, 1.20)
    #qcdloosepass.setParamEffect(qcdnormSF_top, 1.20)
    #qcdfail.setParamEffect(qcdnormSF_top, 1.20)
    qcdpass.setParamEffect(qcdnormSF, 1.20)
    qcdloosepass.setParamEffect(qcdnormSF, 1.20)
    qcdfail.setParamEffect(qcdnormSF, 1.20)
    #qcdpass.setParamEffect(qcdnormSF_top, 1*qcdnormSF_top)
    #qcdloosepass.setParamEffect(qcdnormSF_top, 1*qcdnormSF_top)
    #qcdfail.setParamEffect(qcdnormSF_top, 1*qcdnormSF_top)

    toppass = model['topCRpasshadhad%s'%year]['top']
    toploosepass = model['topCRloosepasshadhad%s'%year]['top']
    topfail = model['topCRfailhadhad%s'%year]['top']
    topPF = (toploosepass.getExpectation(nominal=True).sum() + toppass.getExpectation(nominal=True).sum()) / topfail.getExpectation(nominal=True).sum()
    topLPF = toppass.getExpectation(nominal=True).sum() / toploosepass.getExpectation(nominal=True).sum()
    topRLPF = 1./(1.+(1./topLPF)) # P/(L+P)
    toppass.setParamEffect(topLeffSF, 1*topLeffSF)
    toploosepass.setParamEffect(topLeffSF, (1 - topLeffSF) * topLPF + 1)
    toploosepass.setParamEffect(topeffSF, 1*topeffSF)
    topfail.setParamEffect(topeffSF, (1 - topeffSF) * (1-topRLPF) * topPF + 1)
    #toppass.setParamEffect(topnormSF, 1*topnormSF)
    #toploosepass.setParamEffect(topnormSF, 1*topnormSF)
    #topfail.setParamEffect(topnormSF, 1*topnormSF)

    wlnupass = model['topCRpasshadhad%s'%year]['wlnu']
    wlnuloosepass = model['topCRloosepasshadhad%s'%year]['wlnu']
    wlnufail = model['topCRfailhadhad%s'%year]['wlnu']
    wlnuPF = (wlnuloosepass.getExpectation(nominal=True).sum() + wlnupass.getExpectation(nominal=True).sum()) / wlnufail.getExpectation(nominal=True).sum()
    wlnuLPF = wlnupass.getExpectation(nominal=True).sum() / wlnuloosepass.getExpectation(nominal=True).sum()
    wlnuRLPF = 1./(1.+(1./wlnuLPF)) # P/(L+P)
    wlnupass.setParamEffect(wlnuLeffSF, 1*wlnuLeffSF)
    wlnuloosepass.setParamEffect(wlnuLeffSF, (1 - wlnuLeffSF) * wlnuLPF + 1)
    wlnuloosepass.setParamEffect(wlnueffSF, 1*wlnueffSF)
    wlnufail.setParamEffect(wlnueffSF, (1 - wlnueffSF) * (1-wlnuRLPF) * wlnuPF + 1)
    #wlnupass.setParamEffect(wlnunormSF, 1*wlnunormSF)
    #wlnuloosepass.setParamEffect(wlnunormSF, 1*wlnunormSF)
    #wlnufail.setParamEffect(wlnunormSF, 1*wlnunormSF)

    dypass = model['topCRpasshadhad%s'%year]['dy']
    dyloosepass = model['topCRloosepasshadhad%s'%year]['dy']
    dyfail = model['topCRfailhadhad%s'%year]['dy']
    dyLP = dypass.getExpectation(nominal=True).sum() / dyloosepass.getExpectation(nominal=True).sum()
    httLP = model['topCRpasshadhad%s'%year]['htt125'].getExpectation(nominal=True).sum() / model['topCRloosepasshadhad%s'%year]['htt125'].getExpectation(nominal=True).sum()
    if not doHtt:
        phittLP = {m:model['topCRpasshadhad%s'%year]['phitt%s'%m].getExpectation(nominal=True).sum() / model['topCRloosepasshadhad%s'%year]['phitt%s'%m].getExpectation(nominal=True).sum() if model['topCRloosepasshadhad%s'%year]['phitt%s'%m].getExpectation(nominal=True).sum() > 0. else 1. for m in masspoints}
    dypass.setParamEffect(dy_eff, 1*dy_eff)
    dyloosepass.setParamEffect(dy_eff, (1 - dy_eff)* dyLP + 1)
    model['topCRpasshadhad%s'%year]['htt125'].setParamEffect(dy_eff, 1*dy_eff)
    model['topCRloosepasshadhad%s'%year]['htt125'].setParamEffect(dy_eff, (1 - dy_eff) * httLP + 1)
    if not doHtt:
        for m in masspoints:
            model['topCRpasshadhad%s'%year]['phitt%s'%m].setParamEffect(dy_eff, 1*dy_eff)
            model['topCRloosepasshadhad%s'%year]['phitt%s'%m].setParamEffect(dy_eff, (1 - dy_eff) * phittLP[m] + 1)
    dypass.setParamEffect(rdy, 1*rdy)
    dyloosepass.setParamEffect(rdy, 1*rdy)
    dyfail.setParamEffect(rdy, 1*rdy)
    #dypass.setParamEffect(rdy, 1.05)
    #dyloosepass.setParamEffect(rdy, 1.05)
    #dyfail.setParamEffect(rdy, 1.05)
    if not doHtt:
        model['topCRpasshadhad%s'%year]['htt125'].setParamEffect(rh125, 1.10)
        model['topCRloosepasshadhad%s'%year]['htt125'].setParamEffect(rh125, 1.10)
        model['topCRfailhadhad%s'%year]['htt125'].setParamEffect(rh125, 1.10)

    # Fill in wlnu CR
    for region in ['pass', 'fail', 'loosepass']:
        ch = rl.Channel("wlnuCR%s%s%s" % (region, 'hadhad', year))
        model.addChannel(ch)

        isPass = region=='pass'
        isLoosePass = region=='loosepass'

        thehist = wlnu_cr_hists["nom"]

        for sName in thehist.identifiers('sample'):
            if sName.name=='ignore':
                continue
            tempint = intRegion(thehist[sName],region)
            if singleBinCR:
                templ = (np.array([np.sum(tempint[0][lowmassbin:highmassbin])]), mttone.binning, mttone.name, np.array([np.sum(tempint[1][lowmassbin:highmassbin])]))
                if includeLowMass:
                    templ = (np.array([np.sum(intRegion(thehist[sName],region,mrebin=None))]), mttone.binning, mttone.name)
            else:
                templ = (tempint[0][lowmassbin:highmassbin], mtt.binning, mtt.name, tempint[1][lowmassbin:highmassbin])
            for iv,val in enumerate(templ[0]):
                if val<0.:
                    templ[0][iv] = 0.
            if sName.name=='data_obs':
                ch.setObservation(templ, read_sumw2=True)
            else:
                if sName.name=='multijet':
                    #qcdpred = np.sum(np.stack([(qcd_from_data_wlnu["faildphi"][shaperegion[iregion]][ix][0]*qcdratio_wlnu[shaperegion[iregion]]["faildphi"][ptbin]*qcdratio_F["qcdnom"][region][ix]+qcd_from_data_qcd["nom"][region][ix][0]*qcdratio_qcd[shaperegion[iregion]]["nom"][ptbin]*qcdratio_F["faildphi"][region][ix])/2. for ix in range(len(qcd_from_data_wlnu["faildphi"][region]))]),axis=0)
                    #qcdpred_dn = np.clip(np.sum(np.stack([qcd_from_data_wlnu["faildphi"][shaperegion[iregion]][ix][1] * qcdratio_wlnu[shaperegion[iregion]]["faildphi"][ptbin] * qcdratio_F["qcdnom"][region][ix] for ix in range(len(qcd_from_data_wlnu["faildphi"][region]))]),axis=0),0.,None)
                    #qcdpred_up = np.sum(np.stack([qcd_from_data_qcd["nom"][shaperegion[iregion]][ix][2] * qcdratio_qcd[shaperegion[iregion]]["nom"][ptbin] * qcdratio_F["faildphi"][region][ix] for ix in range(len(qcd_from_data_qcd["nom"][region]))]),axis=0)
                    qcdpred = np.sum(np.stack([qcd_from_data_wlnu["faillep"][shaperegion[iregion]][ix][0]*qcdratio_wlnu[shaperegion[iregion]]["faillep"][ptbin]*qcdratio_F["qcdnom"][region][ix] for ix in range(len(qcd_from_data_wlnu["faillep"][region]))]),axis=0)
                    qcdpred_dn = np.clip(np.sum(np.stack([qcd_from_data_wlnu["faillep"][shaperegion[iregion]][ix][1] * qcdratio_wlnu[shaperegion[iregion]]["faillep"][ptbin] * qcdratio_F["qcdnom"][region][ix] for ix in range(len(qcd_from_data_wlnu["faillep"][region]))]),axis=0),0.,None)
                    qcdpred_up = np.clip(np.sum(np.stack([qcd_from_data_wlnu["faillep"][shaperegion[iregion]][ix][2] * qcdratio_wlnu[shaperegion[iregion]]["faillep"][ptbin] * qcdratio_F["qcdnom"][region][ix] for ix in range(len(qcd_from_data_wlnu["faillep"][region]))]),axis=0),0.,None)

                    if singleBinCR:
                        if includeLowMass:
                            qcdpred = np.sum(np.array(qcd_from_data_full_wlnu[region]), axis=0)
                            qcdpred_dn = np.sum(np.array(qcd_from_data_full_wlnu[region]) - np.array(qcd_from_data_full_wlnu_err[region]), axis=0)
                            qcdpred_up = np.sum(np.array(qcd_from_data_full_wlnu[region]) + np.array(qcd_from_data_full_wlnu_err[region]), axis=0)
                        templ = (np.array([np.sum(qcdpred)]), mttone.binning, mttone.name, np.array([1000000.]))
                    else:
                        templ = (qcdpred, mtt.binning, mtt.name, np.ones_like(qcdpred)*1000000.)
                if doHtt:
                    stype = rl.Sample.SIGNAL if 'htt125' in sName.name else rl.Sample.BACKGROUND
                else:
                    stype = rl.Sample.SIGNAL if 'phi' in sName.name else rl.Sample.BACKGROUND
                sample = rl.TemplateSample(ch.name + '_' + sName.name, stype, templ)
    
                #sample.setParamEffect(jec, 1.05)
                if sName.name in syst_dict:
                    nom_base = intRegion(thehist[sName],region)[0][lowmassbin:highmassbin]
                    if singleBinCR:
                        if includeLowMass:
                            nom_base = intRegion(thehist[sName],region,mrebin=None)
                    for syst in syst_dict[sName.name]:
                        syst_params = syst_dict[sName.name][syst]
                        dnmod = intRegion(thehist[sName],region,systematic=syst_params[1])[0][lowmassbin:highmassbin]
                        upmod = intRegion(thehist[sName],region,systematic=syst_params[2])[0][lowmassbin:highmassbin]
                        if singleBinCR:
                            if includeLowMass:
                                dnmod = intRegion(thehist[sName],region,mrebin=None,systematic=syst_params[1])
                                upmod = intRegion(thehist[sName],region,mrebin=None,systematic=syst_params[2])
                            dnmod = np.array([np.sum(dnmod)])
                            nom = np.array([np.sum(nom_base)])
                            upmod = np.array([np.sum(upmod)])
                        else:
                            nom = nom_base
                        sample.setParamEffect(syst_params[0], np.divide(upmod, nom, out=np.ones_like(nom), where=nom>0.) if (upmod!=nom).all() else np.ones_like(nom)*1.001, np.divide(dnmod, nom, out=np.ones_like(nom), where=nom>0.) if (dnmod!=nom).all() else np.ones_like(nom)*0.999)

                if sName.name=='top':
                    sample.setParamEffect(top_norm, 1.05)
                if sName.name=='wlnu':
                    sample.setParamEffect(wlnu_norm, 1.10)
                if sName.name=='multijet':
                    sample.setParamEffect(qcd_norm, 1.20)
                    if singleBinCR:
                        print(qcdpred_dn)
                        print(qcdpred)
                        print(qcdpred_up)
                        err_dn = np.sqrt(np.sum((qcdpred - qcdpred_dn)*(qcdpred - qcdpred_dn)))
                        err_up = np.sqrt(np.sum((qcdpred_up - qcdpred)*(qcdpred_up - qcdpred)))
                        qcdpred = np.array([np.sum(qcdpred)])
                        qcdpred_dn = qcdpred - err_dn
                        qcdpred_up = qcdpred + err_up
                        print('QCDPRED_COMB ',qcdpred,'+',err_up,'-',err_dn)
                    sample.setParamEffect(qcdwlnu_pass if isPass else qcdwlnu_loosepass if isLoosePass else qcdwlnu_fail, np.divide(qcdpred_dn, qcdpred, out=np.ones_like(qcdpred), where=qcdpred>0.), np.divide(qcdpred_up, qcdpred, out=np.ones_like(qcdpred), where=qcdpred>0.))
                if sName.name=='vvqq':
                    sample.setParamEffect(vvqq_norm, 1.05)
                if sName.name=='dy':
                    sample.setParamEffect(dy_norm, 1.05)
                if sName.name=='htt125' or sName.name=='dy' or 'phi' in sName.name:
                    nom_full = intRegion(thehist[sName],region)
                    shift_dn = nom_full[lowmassbin+1:highmassbin+1 if highmassbin>0 else None]
                    shift_up = nom_full[lowmassbin-1:highmassbin-1]
                    #sample.setParamEffect(m_scale, np.divide(shift_dn, nom_full[lowmassbin:highmassbin], out=np.ones_like(nom_full[lowmassbin:highmassbin]), where=nom_full[lowmassbin:highmassbin]>0.), np.divide(shift_up, nom_full[lowmassbin:highmassbin], out=np.ones_like(nom_full[lowmassbin:highmassbin]), where=nom_full[lowmassbin:highmassbin]>0.))
                sample.setParamEffect(trig, 1.02)
                sample.setParamEffect(lumi, 1.024)
        
                ch.addSample(sample)

    qcdpass = model['wlnuCRpasshadhad%s'%year]['multijet']
    qcdloosepass = model['wlnuCRloosepasshadhad%s'%year]['multijet']
    qcdfail = model['wlnuCRfailhadhad%s'%year]['multijet']
    #qcdpass.setParamEffect(qcdnormSF_wlnu, 1.20)
    #qcdloosepass.setParamEffect(qcdnormSF_wlnu, 1.20)
    #qcdfail.setParamEffect(qcdnormSF_wlnu, 1.20)
    qcdpass.setParamEffect(qcdnormSF, 1.20)
    qcdloosepass.setParamEffect(qcdnormSF, 1.20)
    qcdfail.setParamEffect(qcdnormSF, 1.20)
    #qcdpass.setParamEffect(qcdnormSF_wlnu, 1*qcdnormSF_wlnu)
    #qcdloosepass.setParamEffect(qcdnormSF_wlnu, 1*qcdnormSF_wlnu)
    #qcdfail.setParamEffect(qcdnormSF_wlnu, 1*qcdnormSF_wlnu)

    toppass = model['wlnuCRpasshadhad%s'%year]['top']
    toploosepass = model['wlnuCRloosepasshadhad%s'%year]['top']
    topfail = model['wlnuCRfailhadhad%s'%year]['top']
    topPF = (toploosepass.getExpectation(nominal=True).sum() + toppass.getExpectation(nominal=True).sum()) / topfail.getExpectation(nominal=True).sum()
    topLPF = toppass.getExpectation(nominal=True).sum() / toploosepass.getExpectation(nominal=True).sum()
    topRLPF = 1./(1.+(1./topLPF)) # P/(L+P)
    toppass.setParamEffect(topLeffSF, 1*topLeffSF)
    toploosepass.setParamEffect(topLeffSF, (1 - topLeffSF) * topLPF + 1)
    toploosepass.setParamEffect(topeffSF, 1*topeffSF)
    topfail.setParamEffect(topeffSF, (1 - topeffSF) * (1-topRLPF) * topPF + 1)
    #toppass.setParamEffect(topnormSF, 1*topnormSF)
    #toploosepass.setParamEffect(topnormSF, 1*topnormSF)
    #topfail.setParamEffect(topnormSF, 1*topnormSF)

    wlnupass = model['wlnuCRpasshadhad%s'%year]['wlnu']
    wlnuloosepass = model['wlnuCRloosepasshadhad%s'%year]['wlnu']
    wlnufail = model['wlnuCRfailhadhad%s'%year]['wlnu']
    wlnuPF = (wlnuloosepass.getExpectation(nominal=True).sum() + wlnupass.getExpectation(nominal=True).sum()) / wlnufail.getExpectation(nominal=True).sum()
    wlnuLPF = wlnupass.getExpectation(nominal=True).sum() / wlnuloosepass.getExpectation(nominal=True).sum()
    wlnuRLPF = 1./(1.+(1./wlnuLPF)) # P/(L+P)
    wlnupass.setParamEffect(wlnuLeffSF, 1*wlnuLeffSF)
    wlnuloosepass.setParamEffect(wlnuLeffSF, (1 - wlnuLeffSF) * wlnuLPF + 1)
    wlnuloosepass.setParamEffect(wlnueffSF, 1*wlnueffSF)
    wlnufail.setParamEffect(wlnueffSF, (1 - wlnueffSF) * (1-wlnuRLPF) * wlnuPF + 1)
    #wlnupass.setParamEffect(wlnunormSF, 1*wlnunormSF)
    #wlnuloosepass.setParamEffect(wlnunormSF, 1*wlnunormSF)
    #wlnufail.setParamEffect(wlnunormSF, 1*wlnunormSF)

    dypass = model['wlnuCRpasshadhad%s'%year]['dy']
    dyloosepass = model['wlnuCRloosepasshadhad%s'%year]['dy']
    dyfail = model['wlnuCRfailhadhad%s'%year]['dy']
    dyLP = dypass.getExpectation(nominal=True).sum() / dyloosepass.getExpectation(nominal=True).sum()
    httLP = model['wlnuCRpasshadhad%s'%year]['htt125'].getExpectation(nominal=True).sum() / model['wlnuCRloosepasshadhad%s'%year]['htt125'].getExpectation(nominal=True).sum()
    if not doHtt:
        phittLP = {m:model['wlnuCRpasshadhad%s'%year]['phitt%s'%m].getExpectation(nominal=True).sum() / model['wlnuCRloosepasshadhad%s'%year]['phitt%s'%m].getExpectation(nominal=True).sum() if model['wlnuCRloosepasshadhad%s'%year]['phitt%s'%m].getExpectation(nominal=True).sum() > 0. else 1. for m in masspoints}
    dypass.setParamEffect(dy_eff, 1*dy_eff)
    dyloosepass.setParamEffect(dy_eff, (1 - dy_eff) * dyLP + 1)
    model['wlnuCRpasshadhad%s'%year]['htt125'].setParamEffect(dy_eff, 1*dy_eff)
    model['wlnuCRloosepasshadhad%s'%year]['htt125'].setParamEffect(dy_eff, (1 - dy_eff) * httLP + 1)
    if not doHtt:
        for m in masspoints:
            model['wlnuCRpasshadhad%s'%year]['phitt%s'%m].setParamEffect(dy_eff, 1*dy_eff)
            model['wlnuCRloosepasshadhad%s'%year]['phitt%s'%m].setParamEffect(dy_eff, (1 - dy_eff) * phittLP[m] + 1)
    dypass.setParamEffect(rdy, 1*rdy)
    dyloosepass.setParamEffect(rdy, 1*rdy)
    dyfail.setParamEffect(rdy, 1*rdy)
    #dypass.setParamEffect(rdy, 1.05)
    #dyloosepass.setParamEffect(rdy, 1.05)
    #dyfail.setParamEffect(rdy, 1.05)
    if not doHtt:
        model['wlnuCRpasshadhad%s'%year]['htt125'].setParamEffect(rh125, 1.10)
        model['wlnuCRloosepasshadhad%s'%year]['htt125'].setParamEffect(rh125, 1.10)
        model['wlnuCRfailhadhad%s'%year]['htt125'].setParamEffect(rh125, 1.10)

    with open(os.path.join("%s/%s"%(str(tmpdir),label), 'hadhadModel.pkl'), "wb") as fout:
        pickle.dump(model, fout, protocol=2)

    #model.renderCombine(os.path.join("%s/%s"%(str(tmpdir),label), 'hadhadModel'))

def getHist(h,var_name,lumifb,vars_cut,regionsel,blind,sigscale,rebin,debug=False):
    if debug: print(h)
    exceptions = ['process', var_name]
    for var in vars_cut:
        exceptions.append(var)
    if (regionsel!=''):
        exceptions.append('region')
    if debug: print([ax.name for ax in h.axes()])
    x = h.sum(*[ax for ax in h.axes() if ax.name not in exceptions],overflow=overflow_sum)
    for reg in regionsel:
        if debug: print('integrating ',reg)
        x = x.integrate('region',reg)
    for var,val in vars_cut.items():
        if len(val)==0:
            continue
        if var!=var_name:
            if debug: print('integrating ',var,val[0],val[1])
            if val[0] is not None and val[1] is not None:
                overflow_str = 'none'
            elif val[0] is None and val[1] is not None:
                overflow_str = 'under'
            elif val[0] is not None and val[1] is None:
                overflow_str = 'over'
            else:
                overflow_str = 'allnan'
            x = x.integrate(var,slice(val[0],val[1]),overflow=overflow_str)
            #x = x.integrate(var,slice(val[0],val[1]),overflow=overflow)
    if var_name in vars_cut.keys():
        if debug: print('integrating ',var_name,vars_cut[var_name][0],vars_cut[var_name][1])
        x = x.rebin(var_name, hist.Bin(var_name, var_name, [xe for xe in x.axis(var_name).edges() if (xe >= vars_cut[var_name][0] and xe <= vars_cut[var_name][1])]))
        #x = x.rebin(var_name, x.axis(var_name).reduced(slice(vars_cut[var_name][0],vars_cut[var_name][1])))

    xaxis = var_name
    for ih,hkey in enumerate(x.identifiers('process')):
        x.identifiers('process')[ih].label = process_latex[hkey.name]

    if len(rebin)==1:
        if (rebin[0]>1): 
            x = x.rebin(xaxis,int(rebin[0]))
    else:
        x = x.rebin(xaxis,hist.Bin(xaxis, var_name, rebin))
 

    x_nobkg = x[nobkg]
    x_nosig = x[nosig]

    # normalize to lumi
    x_nosig.scale({p: lumifb for p in x_nosig.identifiers('process')}, axis="process")
    x_nobkg.scale({p: lumifb*float(sigscale) for p in x_nobkg.identifiers('process')}, axis="process")

    x_data = x['data']

    return x_nobkg+x_nosig+x_data

def makeCards(args):

    mcSamples = samplelists.getSamplesMC(args.year)
    dataSamples = samplelists.getSamplesData(args.year)

    tag = args.tag

    odir = 'cards/%s/'%tag
    os.system('mkdir -p %s/%s'%(odir,args.label))
    pwd = os.getcwd()

    nnCut_hadhad = args.nnCutHH
    nnCut_hadhad_loose = args.nnCutHHL
    nnCut_hadhad_met = args.nnCutHHMET
    nnCut_hadhad_met_loose = args.nnCutHHMETL
    nnCut_hadel = args.nnCutEH
    nnCut_hadel_loose = args.nnCutEHL
    nnCut_hadmu = args.nnCutMH
    nnCut_hadmu_loose = args.nnCutMHL
    metCutLep = args.metCutLep
    lowMetCutLep = args.lowMetCutLep
    metCutHad = args.metCutHad
    lowMetCutHad = args.lowMetCutHad
    h_pt_min = args.hPtCut
    
    includeData = False
    signalmasses = ['10', '20', '30', '40', '50', '75', '100', '125', '150', '200', '250', '300']
    #mttbins = np.array([0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,200.,250.,300.,350.,400.,450.,500.])
    mttbins = np.array([0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,200.,250.,300.,350.,400.])
    # properties
    region_dict = {
        'hadhad_signal':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadhad_signal'],
        'hadhad_signal_met':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadhad_signal_met'],
        'hadhad_cr_mu':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[],'met_pt':[],'systematic':[]}, includeData, [1], 'hadhad_cr_mu'],
        'hadhad_cr_b_met':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[],'met_pt':[],'systematic':[]}, includeData, [1], 'hadhad_cr_b_met'],
        'hadhad_cr_b_mu_iso':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[],'met_pt':[],'systematic':[]}, includeData, [1], 'hadhad_cr_b_mu_iso'],
        'hadhad_cr_mu_iso':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[],'met_pt':[],'systematic':[]}, includeData, [1], 'hadhad_cr_mu_iso'],
        'hadhad_cr_dphi_inv':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadhad_cr_dphi_inv'],
        'hadhad_cr_b_met_dphi_inv':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[],'met_pt':[],'systematic':[]}, includeData, [1], 'hadhad_cr_b_met_dphi_inv'],
        'hadhad_cr_mu_dphi_inv':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[],'met_pt':[],'systematic':[]}, includeData, [1], 'hadhad_cr_mu_dphi_inv'],
        'hadhad_cr_b_mu_iso_dphi_inv':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[],'met_pt':[],'systematic':[]}, includeData, [1], 'hadhad_cr_b_mu_iso_dphi_inv'],
        'hadhad_cr_mu_iso_dphi_inv':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[],'met_pt':[],'systematic':[]}, includeData, [1], 'hadhad_cr_mu_iso_dphi_inv'],
        'hadhad_cr_anti_inv':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadhad_cr_anti_inv'],
        'hadhad_cr_b_met_anti_inv':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[],'met_pt':[],'systematic':[]}, includeData, [1], 'hadhad_cr_b_met_anti_inv'],
        'hadhad_cr_mu_anti_inv':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[],'met_pt':[],'systematic':[]}, includeData, [1], 'hadhad_cr_mu_anti_inv'],
        'hadhad_cr_b_mu_iso_anti_inv':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[],'met_pt':[],'systematic':[]}, includeData, [1], 'hadhad_cr_b_mu_iso_anti_inv'],
        'hadhad_cr_mu_iso_anti_inv':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[],'met_pt':[],'systematic':[]}, includeData, [1], 'hadhad_cr_mu_iso_anti_inv'],
        'hadel_signal':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadel_signal'],
        'hadmu_signal':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadmu_signal'],
        'hadel_cr_ztag_inv':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadel_cr_ztag_inv'],
        'hadmu_cr_ztag_inv':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadmu_cr_ztag_inv'],
        'hadel_cr_dphi_inv':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadel_cr_dphi_inv'],
        'hadmu_cr_dphi_inv':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadmu_cr_dphi_inv'],
        'hadel_cr_b':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadel_cr_b'],
        'hadmu_cr_b':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadmu_cr_b'],
        'hadel_cr_w':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadel_cr_w'],
        'hadmu_cr_w':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadmu_cr_w'],
        'hadel_cr_qcd':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadel_cr_qcd'],
        'hadmu_cr_qcd':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadmu_cr_qcd'],
        'hadel_cr_b_ztag_inv':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadel_cr_b_ztag_inv'],
        'hadmu_cr_b_ztag_inv':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadmu_cr_b_ztag_inv'],
        'hadel_cr_w_ztag_inv':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadel_cr_w_ztag_inv'],
        'hadmu_cr_w_ztag_inv':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadmu_cr_w_ztag_inv'],
        'hadel_cr_qcd_ztag_inv':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadel_cr_qcd_ztag_inv'],
        'hadmu_cr_qcd_ztag_inv':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadmu_cr_qcd_ztag_inv'],
        'hadel_cr_b_dphi_inv':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadel_cr_b_dphi_inv'],
        'hadmu_cr_b_dphi_inv':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadmu_cr_b_dphi_inv'],
        'hadel_cr_w_dphi_inv':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadel_cr_w_dphi_inv'],
        'hadmu_cr_w_dphi_inv':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadmu_cr_w_dphi_inv'],
        'hadel_cr_qcd_dphi_inv':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadel_cr_qcd_dphi_inv'],
        'hadmu_cr_qcd_dphi_inv':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadmu_cr_qcd_dphi_inv'],
    }

    xs = {}
    with open('fileset/xsecs.json', 'r') as f:
        xs = json.load(f)
    scale1 = {k: args.sigxs/xs[k] for k in xs if 'Spin0' in k and args.sigxs>0.}
    print(scale1)

    full_dict = {}
    h_dict = {reg:None for reg in region_dict}

    for hs in mcSamples + dataSamples:
        # open hists
        try:
            with open("%s%s.hist"%(args.hist,hs), 'rb') as f:
                hists_unmapped = pickle.load(f)
        except:
            print('ERROR: unable to open file %s%s.hist'%(args.hist,hs))

        for reg in region_dict:
            reg_opt = region_dict[reg]
     
            # map to hists
            for key, val in hists_unmapped.items():
                if isinstance(val, hist.Hist) and key==reg_opt[0]:
                    toadd = val
                    toadd.scale({p: scale1[p.name] if p.name in scale1 else 1. for p in toadd.identifiers('dataset')}, axis="dataset")
                    if h_dict[reg] is not None:
                        h_dict[reg] = h_dict[reg].add(processmap.apply(toadd).integrate('region',reg_opt[5]))
                    else:
                        h_dict[reg] = processmap.apply(toadd).integrate('region',reg_opt[5])

    for reg in region_dict:
        reg_opt = region_dict[reg]
        #print(h)
        xhist = getHist(h_dict[reg],reg_opt[1],args.lumi,reg_opt[2],'',reg_opt[3],args.sigscale,reg_opt[4],args.sigxs)
        #print(xhist)
        full_dict[reg] = xhist

    hadptbins = np.array([h_pt_min]+[float(b) if b!="None" else None for b in args.hPtBinsHad])
    lepptbins = np.array([h_pt_min]+[float(b) if b!="None" else None for b in args.hPtBinsLep])

    #createHadHad(full_dict['hadhad_signal_met'], full_dict['hadhad_cr_b_mu_iso'], full_dict['hadhad_cr_mu_iso'], full_dict['hadhad_cr_mu'], full_dict['hadhad_cr_anti_inv'], full_dict['hadhad_cr_b_mu_iso_anti_inv'], full_dict['hadhad_cr_mu_iso_anti_inv'], full_dict['hadhad_cr_mu_anti_inv'], full_dict['hadhad_cr_dphi_inv'], full_dict['hadhad_cr_b_mu_iso_dphi_inv'], full_dict['hadhad_cr_mu_iso_dphi_inv'], full_dict['hadhad_cr_mu_dphi_inv'], 'massreg', mttbins, hadptbins, odir, args.label, args.unblind, nnCut_hadhad_met, nnCut_hadhad_met_loose, metCutHad, lowMetCutHad, h_pt_min, signalmasses, args.shapeRegionsHad)
    if not args.noHad: 
        createHadHad(full_dict['hadhad_signal_met'], full_dict['hadhad_cr_b_met'], full_dict['hadhad_cr_mu_iso'], full_dict['hadhad_cr_mu'], full_dict['hadhad_cr_anti_inv'], full_dict['hadhad_cr_b_met_anti_inv'], full_dict['hadhad_cr_mu_iso_anti_inv'], full_dict['hadhad_cr_mu_anti_inv'], full_dict['hadhad_cr_dphi_inv'], full_dict['hadhad_cr_b_met_dphi_inv'], full_dict['hadhad_cr_mu_iso_dphi_inv'], full_dict['hadhad_cr_mu_dphi_inv'], 'massreg', mttbins, hadptbins, odir, args.label, args.unblind, nnCut_hadhad_met, nnCut_hadhad_met_loose, metCutHad, lowMetCutHad, h_pt_min, signalmasses, args.shapeRegionsHad, args.year, args.doHtt, args.noSyst)
    if not args.noLep:
        createLepHad(full_dict['hadel_signal'], full_dict['hadel_cr_b'], full_dict['hadel_cr_w'], full_dict['hadel_cr_qcd'], full_dict['hadel_cr_ztag_inv'], full_dict['hadel_cr_b_ztag_inv'], full_dict['hadel_cr_w_ztag_inv'], full_dict['hadel_cr_qcd_ztag_inv'], full_dict['hadel_cr_dphi_inv'], full_dict['hadel_cr_b_dphi_inv'], full_dict['hadel_cr_w_dphi_inv'], full_dict['hadel_cr_qcd_dphi_inv'], 'massreg', mttbins, lepptbins, "el", odir, args.label, args.unblind, nnCut_hadel, nnCut_hadel_loose, metCutLep, lowMetCutLep, h_pt_min, signalmasses, args.shapeRegionsLep, args.year, args.doHtt, args.noSyst)
        createLepHad(full_dict['hadmu_signal'], full_dict['hadmu_cr_b'], full_dict['hadmu_cr_w'], full_dict['hadmu_cr_qcd'], full_dict['hadmu_cr_ztag_inv'], full_dict['hadmu_cr_b_ztag_inv'], full_dict['hadmu_cr_w_ztag_inv'], full_dict['hadmu_cr_qcd_ztag_inv'], full_dict['hadmu_cr_dphi_inv'], full_dict['hadmu_cr_b_dphi_inv'], full_dict['hadmu_cr_w_dphi_inv'], full_dict['hadmu_cr_qcd_dphi_inv'], 'massreg', mttbins, lepptbins, "mu", odir, args.label, args.unblind, nnCut_hadmu, nnCut_hadmu_loose, metCutLep, lowMetCutLep, h_pt_min, signalmasses, args.shapeRegionsLep, args.year, args.doHtt, args.noSyst)

    for ax1 in axs:
        for ax in ax1:
            ax.legend()
            ax.set_ylim(0,3)
            ax.set_ylabel('QCD ratio')
            ax.set_xlabel('m [GeV]')

    fig.savefig("%s/%s/qcdratio.png"%(odir,args.label))
    fig.savefig("%s/%s/qcdratio.pdf"%(odir,args.label))

    for ax1 in axs_shape:
        for ax in ax1:
            ax.legend()
            #ax.set_ylim(0,3)
            ax.set_ylabel('QCD prediction')
            ax.set_xlabel('m [GeV]')
    axs_shape[3,0].set_ylim(0.,0.5)
    axs_shape[3,1].set_ylim(0.,0.5)

    fig_shape.savefig("%s/%s/qcdpred.png"%(odir,args.label))
    fig_shape.savefig("%s/%s/qcdpred.pdf"%(odir,args.label))

    import datetime

    f = open("%s/%s/info.txt"%(odir,args.label), "w")
    f.write("created: " + str(datetime.datetime.now()))
    f.write("\nhist: " + args.hist)
    f.write("\nyear: " + args.year)
    f.write("\nlumi: " + str(args.lumi))
    f.write("\ntag: " + args.tag)
    f.write("\nlabel: " + args.label)
    f.write("\nsigscale: " + str(args.sigscale))
    f.write("\nunblind: " + str(args.unblind))
    f.write("\nHadHad cuts: %s, %s"%(args.nnCutHH,args.nnCutHHL))
    f.write("\nHadHadMet cuts: %s, %s"%(args.nnCutHHMET,args.nnCutHHMETL))
    f.write("\nHadMu cuts: %s, %s"%(args.nnCutMH,args.nnCutMHL))
    f.write("\nHadEl cuts: %s, %s"%(args.nnCutEH,args.nnCutEHL))
    f.write("\nLow MET cut: " + str(args.lowMetCutLep))
    f.write("\nMET cut: " + str(args.metCutLep))
    f.write("\nLow MET cut: " + str(args.lowMetCutHad))
    f.write("\nMET cut: " + str(args.metCutHad))
    f.write("\nhPt cut: " + str(args.hPtCut))
    f.write("\nLepBins: " + ", ".join(args.hPtBinsLep))
    f.write("\nHadBins: " + ", ".join(args.hPtBinsHad))
    f.write("\nShape Regions (Lep): " + ", ".join(args.shapeRegionsLep))
    f.write("\nShape Regions (Had): " + ", ".join(args.shapeRegionsHad))
    f.close()

    if args.plots:
    
        pwd = os.getcwd()
        os.chdir('%s/%s'%(odir,args.label))
    
        for rn in ['hadhad_signal_met','hadhad_cr_b_mu_iso','hadhad_cr_mu_iso','hadhad_cr_mu']:
          if rn=='hadhad_signal':
            nncutl = nnCut_hadhad_loose
            nncutt = nnCut_hadhad
          elif rn=='hadhad_signal_met':
            nncutl = nnCut_hadhad_met_loose
            nncutt = nnCut_hadhad_met
          else:
            nncutl = nnCut_hadhad_met_loose
            nncutt = nnCut_hadhad_met
          plot_stack.drawStack(full_dict[rn],'HadhadModel','massreg',r"$m_{NN}$",'%s fail'%rn,'all',1.,{'nn_disc':[None,nncutl]},{},[],rn+'_fail')
          plot_stack.drawStack(full_dict[rn],'HadhadModel','massreg',r"$m_{NN}$",'%s loosepass'%rn,'all',1.,{'nn_disc':[nncutl,nncutt]},{},[],rn+'_loosepass')
          plot_stack.drawStack(full_dict[rn],'HadhadModel','massreg',r"$m_{NN}$",'%s pass'%rn,'all',1.,{'nn_disc':[nncutt,None]},{},[],rn+'_pass', blind='' if args.unblind or 'signal' not in rn else ['None',110.])
        for lep in ['el','mu']:
          if lep=='el':
            nncutl = nnCut_hadel_loose
            nncutt = nnCut_hadel
          else:
            nncutl = nnCut_hadmu_loose
            nncutt = nnCut_hadmu
          for rn in ['had%s_signal'%lep,'had%s_cr_b'%lep,'had%s_cr_w'%lep,'had%s_cr_qcd'%lep]:
            plot_stack.drawStack(full_dict[rn],'Had%sModel'%lep,'massreg',r"$m_{NN}$",'%s fail'%rn,'all',1.,{'nn_disc':[None,nncutl],'met_pt':[metCutLep,None]},{},[],rn+'_fail')
            plot_stack.drawStack(full_dict[rn],'Had%sModel'%lep,'massreg',r"$m_{NN}$",'%s loosepass'%rn,'all',1.,{'nn_disc':[nncutl,nncutt],'met_pt':[metCutLep,None]},{},[],rn+'_loosepass')
            plot_stack.drawStack(full_dict[rn],'Had%sModel'%lep,'massreg',r"$m_{NN}$",'%s pass'%rn,'all',1.,{'nn_disc':[nncutt,None],'met_pt':[metCutLep,None]},{},[],rn+'_pass', blind='' if args.unblind or 'signal' not in rn else ['None',110.])
    
        os.chdir(pwd)

if __name__ == "__main__":

    if not sys.warnoptions:
        import warnings
        warnings.simplefilter("ignore")
        os.environ["PYTHONWARNINGS"] = "ignore" # Also affect subprocesses

    #ex. python makeCards.py --hist ../../boostedhiggs/condor/Mar08_NN/hists_sum_ --year 2017 --lumi 36.7 --tag Mar08 --label 01
    parser = argparse.ArgumentParser()
    parser.add_argument('--hist',                 dest='hist',         default="hists_sum_",          help="hists pickle prefix")
    parser.add_argument('--year',                 dest='year',         default="2017",                help="year")
    parser.add_argument('--lumi',                 dest='lumi',         default=50.,                   help="lumi",       type=float)
    parser.add_argument('--tag',                  dest='tag',          default="",                    help="tag")
    parser.add_argument('--label',                dest='label',        default="",                    help="label")
    parser.add_argument('--sigscale',             dest='sigscale',     default=1.,                    help="sigscale",   type=float)
    parser.add_argument('--unblind',              dest='unblind',      action="store_true",           help="unblind")
    parser.add_argument('--plots',                dest='plots',        action="store_true",           help="plots")
    parser.add_argument('--nnCutHH',              dest='nnCutHH',      default=0.9999,                help="nncut hadhad",   type=float)
    parser.add_argument('--nnCutHHMET',           dest='nnCutHHMET',   default=0.9999,                help="nncut hadhad met",   type=float)
    parser.add_argument('--nnCutHHL',             dest='nnCutHHL',     default=0.995,                 help="nncut hadhad loose",   type=float)
    parser.add_argument('--nnCutHHMETL',          dest='nnCutHHMETL',  default=0.995,                 help="nncut hadhad met loose",   type=float)
    parser.add_argument('--nnCutEH',              dest='nnCutEH',      default=0.98,                  help="nncut elhad",   type=float)
    parser.add_argument('--nnCutMH',              dest='nnCutMH',      default=0.98,                  help="nncut muhad",   type=float)
    parser.add_argument('--nnCutEHL',             dest='nnCutEHL',     default=0.9,                   help="nncut elhad loose",   type=float)
    parser.add_argument('--nnCutMHL',             dest='nnCutMHL',     default=0.9,                   help="nncut muhad loose",   type=float)
    parser.add_argument('--metCutLep',            dest='metCutLep',    default=50.,                   help="met cut (lep)",   type=float)
    parser.add_argument('--lowMetCutLep',         dest='lowMetCutLep', default=20.,                   help="low met cut (lep)",   type=float)
    parser.add_argument('--metCutHad',            dest='metCutHad',    default=150.,                  help="met cut (had)",   type=float)
    parser.add_argument('--lowMetCutHad',         dest='lowMetCutHad', default=75.,                  help="low met cut (had)",   type=float)
    parser.add_argument('--hPtBinsHad',           dest='hPtBinsHad',   default=['300','None'],        help="h pt bins (had)",   type=str,  nargs='+')
    parser.add_argument('--hPtBinsLep',           dest='hPtBinsLep',   default=['300','None'],        help="h pt bins (lep)",   type=str,  nargs='+')
    parser.add_argument('--hPtCut',               dest='hPtCut',       default=250.,                  help="h pt cut",   type=float)
    parser.add_argument('--sigxs',                dest='sigxs',        default=0.,                    help="signal xs to scale to",   type=float)
    parser.add_argument('--shapeRegionsLep',      dest='shapeRegionsLep',      default=["fail","fail","fail"],                  help="shape regions (lep)",   type=str, nargs=3)
    parser.add_argument('--shapeRegionsHad',      dest='shapeRegionsHad',      default=["fail","fail","fail"],                  help="shape regions (had)",   type=str, nargs=3)
    parser.add_argument('--noHad',              dest='noHad',      action="store_true",           help="noHad")
    parser.add_argument('--noLep',              dest='noLep',      action="store_true",           help="noLep")
    parser.add_argument('--doHtt',              dest='doHtt',      action="store_true",           help="doHtt")
    parser.add_argument('--noSyst',             dest='noSyst',     action="store_true",           help="noSyst")

    args = parser.parse_args()

    makeCards(args)
