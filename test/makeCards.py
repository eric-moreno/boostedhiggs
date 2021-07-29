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
from cycler import cycler

import argparse
import processmap
from hists_map import *


import rhalphalib as rl
import scipy.stats
from scipy.special import erfinv
import pickle
import ROOT
rl.util.install_roofit_helpers()
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
includeLowMass = True

def createLepHad(sig_hist, top_cr_hist, wlnu_cr_hist, qcd_cr_hist, var_name, mttbins, ptbins, leptype, tmpdir, label, usingData, nnCut, nnCut_loose, metCut, h_pt_min, antilepcut):

    #could be made region-dependent

    samp_combinations = {
        'data_obs'   : ['data'],
        'top'        : ['tt-had', 'tt-semilep', 'tt-dilep','st'],
        'htt125'     : ['h125'],
        'multijet'   : ['qcd'],
        'ztt'        : ['zll'],
        'wlnu'       : ['wjets'],
        'vvqq'       : ['vv','vqq'],
        'ignore'     : [],
    }

    sig_faillep = sig_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('antilep',slice(None,antilepcut),'under').integrate('met_pt',slice(metCut,None),'over')
    top_cr_faillep = top_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('antilep',slice(None,antilepcut),'under').integrate('met_pt',slice(metCut,None),'over')
    wlnu_cr_faillep = wlnu_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('antilep',slice(None,antilepcut),'under').integrate('met_pt',slice(metCut,None),'over')
    qcd_cr_faillep = qcd_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('antilep',slice(None,antilepcut),'under').integrate('met_pt',slice(metCut,None),'over')

    sig_lowmet = sig_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('antilep',slice(antilepcut,None),'over').integrate('met_pt',slice(None,metCut),'under')
    top_cr_lowmet = top_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('antilep',slice(antilepcut,None),'over').integrate('met_pt',slice(None,metCut),'under')
    wlnu_cr_lowmet = wlnu_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('antilep',slice(antilepcut,None),'over').integrate('met_pt',slice(None,metCut),'under')
    qcd_cr_lowmet = qcd_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('antilep',slice(antilepcut,None),'over').integrate('met_pt',slice(None,metCut),'under')

    sig_hist = sig_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('antilep',slice(antilepcut,None),'over').integrate('met_pt',slice(metCut,None),'over')
    top_cr_hist = top_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('antilep',slice(antilepcut,None),'over').integrate('met_pt',slice(metCut,None),'over')
    wlnu_cr_hist = wlnu_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('antilep',slice(antilepcut,None),'over').integrate('met_pt',slice(metCut,None),'over')
    qcd_cr_hist = qcd_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('antilep',slice(antilepcut,None),'over').integrate('met_pt',slice(metCut,None),'over')

    #jec = rl.NuisanceParameter('CMS_jec', 'lnN')
    top_norm = rl.NuisanceParameter('CMS_top_norm', 'lnN')
    wlnu_norm = rl.NuisanceParameter('CMS_wlnu_norm', 'lnN')
    qcd_norm = rl.NuisanceParameter('CMS_qcd_norm', 'lnN')
    #ztt_norm = rl.NuisanceParameter('CMS_ztt_norm', 'lnN')
    vqq_norm = rl.NuisanceParameter('CMS_vqq_norm', 'lnN')
    trig = rl.NuisanceParameter('CMS_trig_had%s'%leptype, 'lnN')
    lumi = rl.NuisanceParameter('CMS_lumi', 'lnN')

    qcd_fail = rl.NuisanceParameter('qcd_Rfail_had%s'%leptype, 'shape')
    qcd_loosepass = rl.NuisanceParameter('qcd_Rloosepass_had%s'%leptype, 'shape')
    qcd_pass = rl.NuisanceParameter('qcd_Rpass_had%s'%leptype, 'shape')
    #qcdnormSF = rl.IndependentParameter('qcdnormSF_had%s'%leptype, 1., 0, 10)
    qcdnormSF = rl.NuisanceParameter('qcdnormSF_had%s'%leptype, 'lnN')
    qcdtop_fail = rl.NuisanceParameter('qcd_top_Rfail_had%s'%leptype, 'shape')
    qcdtop_loosepass = rl.NuisanceParameter('qcd_top_Rloosepass_had%s'%leptype, 'shape')
    qcdtop_pass = rl.NuisanceParameter('qcd_top_Rpass_had%s'%leptype, 'shape')
    #qcdnormSF_top = rl.IndependentParameter('qcdnormSF_top_had%s'%leptype, 1., 0, 10)
    qcdnormSF_top = rl.NuisanceParameter('qcdnormSF_top_had%s'%leptype, 'lnN')
    qcdwlnu_fail = rl.NuisanceParameter('qcd_wlnu_Rfail_had%s'%leptype, 'shape')
    qcdwlnu_loosepass = rl.NuisanceParameter('qcd_wlnu_Rloosepass_had%s'%leptype, 'shape')
    qcdwlnu_pass = rl.NuisanceParameter('qcd_wlnu_Rpass_had%s'%leptype, 'shape')
    #qcdnormSF_wlnu = rl.IndependentParameter('qcdnormSF_wlnu_had%s'%leptype, 1., 0, 10)
    qcdnormSF_wlnu = rl.NuisanceParameter('qcdnormSF_wlnu_had%s'%leptype, 'lnN')

    #m_scale = rl.NuisanceParameter('massscale_had%s'%leptype, 'shape')

    topeffSF = rl.IndependentParameter('topeffSF_had%s'%leptype, 1., 0, 10)
    topnormSF = rl.IndependentParameter('topnormSF_had%s'%leptype, 1., 0, 10)
    wlnueffSF = rl.IndependentParameter('wlnueffSF_had%s'%leptype, 1., 0, 10)
    wlnunormSF = rl.IndependentParameter('wlnunormSF_had%s'%leptype, 1., 0, 10)

    topLeffSF = rl.IndependentParameter('topLeffSF_had%s'%leptype, 1., 0, 10)
    wlnuLeffSF = rl.IndependentParameter('wlnuLeffSF_had%s'%leptype, 1., 0, 10)

    rztt = rl.IndependentParameter('r_ztt_had%s'%leptype, 1., 0, 10)
    #rztt = rl.NuisanceParameter('r_ztt_had%s'%leptype, 'lnN')
    ztt_eff = rl.IndependentParameter('ztt_eff_had%s'%leptype, 1., 0, 10)

    toppt = rl.NuisanceParameter('toppt', 'shape')
    syst_dict = {
        "top":{
            'toppt':[toppt,'nominal','TopPtReweightUp'],
        },
    }

    npt = len(ptbins) - 1
    ptslices = [slice(ptbins[ipt],ptbins[ipt+1]) for ipt in range(npt)]
    mtt = rl.Observable(var_name, mttbins[1:-1])
    mttone = rl.Observable(var_name+'_one', np.array([mttbins[1],mttbins[-2]]))

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
            the_int = the_int.integrate('sample',samplelist).values()[()]
        else:
            the_int = the_int.sum('sample').values()[()]

        if debug:
            print('\tdebug',the_int)

        return the_int

    def getQCDfromData(inhist,region,default=1.,systematic='nominal',hptslice=slice(h_pt_min,None),mslice=None):
        qcd_data = intRegion(inhist['data_obs'],region,systematic=systematic,hptslice=hptslice,mslice=mslice)[1:-1]
        other_mc = intRegion(inhist,region,samplelist=[s.name for s in inhist.identifiers('sample') if s.name!='data_obs' and s.name!='multijet'],systematic=systematic,hptslice=hptslice,mslice=mslice)[1:-1]
        qcd_data = qcd_data - other_mc
        qcd_data = np.array([v if v>0. else default for v in qcd_data])
        qcd_data_int = hist.poisson_interval(qcd_data,qcd_data)
        qcd_temp = qcd_data
        qcd_temp_dn = np.array([qcd_data_int[1][bi] if qcd_data_int[1][bi] > 0. and not np.isnan(qcd_data_int[1][bi]) else default for bi in range(len(qcd_data))])
        qcd_temp_up = np.array([qcd_data_int[0][bi] if qcd_data_int[0][bi] > 0. and not np.isnan(qcd_data_int[0][bi]) else default for bi in range(len(qcd_data))])
        return qcd_temp,qcd_temp_dn,qcd_temp_up

    qcd_data_faillep = [getQCDfromData(qcd_cr_faillep,"fail",hptslice=slice(qcd_cr_faillep.axis('h_pt').edges()[ix-1] if ix>1 else None,qcd_cr_faillep.axis('h_pt').edges()[ix] if ix<len(qcd_cr_faillep.axis('h_pt').edges())-1 else None))[0] for ix in range(len(qcd_cr_faillep.axis('h_pt').edges()))]
    top_data_faillep = [getQCDfromData(top_cr_faillep,"fail",hptslice=slice(top_cr_faillep.axis('h_pt').edges()[ix-1] if ix>1 else None,top_cr_faillep.axis('h_pt').edges()[ix] if ix<len(top_cr_faillep.axis('h_pt').edges())-1 else None))[0] for ix in range(len(top_cr_faillep.axis('h_pt').edges()))]
    wlnu_data_faillep = [getQCDfromData(wlnu_cr_faillep,"fail",hptslice=slice(wlnu_cr_faillep.axis('h_pt').edges()[ix-1] if ix>1 else None,wlnu_cr_faillep.axis('h_pt').edges()[ix] if ix<len(wlnu_cr_faillep.axis('h_pt').edges())-1 else None))[0] for ix in range(len(wlnu_cr_faillep.axis('h_pt').edges()))]
    sig_data_faillep = [getQCDfromData(sig_faillep,"fail",hptslice=slice(sig_faillep.axis('h_pt').edges()[ix-1] if ix>1 else None,sig_faillep.axis('h_pt').edges()[ix] if ix<len(sig_faillep.axis('h_pt').edges())-1 else None))[0] for ix in range(len(sig_faillep.axis('h_pt').edges()))]

    #qcdratio_sig_faillep = [sig_data_faillep[ix]/qcd_data_faillep[ix] for ix in range(len(qcd_data_faillep))]
    #qcdratio_top_faillep = [top_data_faillep[ix]/qcd_data_faillep[ix] for ix in range(len(qcd_data_faillep))]
    #qcdratio_wlnu_faillep = [wlnu_data_faillep[ix]/qcd_data_faillep[ix] for ix in range(len(qcd_data_faillep))]
    qcdratio_sig_faillep = [np.sum(sig_data_faillep[ix])/np.sum(qcd_data_faillep[ix]) for ix in range(len(qcd_data_faillep))]
    qcdratio_top_faillep = [np.sum(top_data_faillep[ix])/np.sum(qcd_data_faillep[ix]) for ix in range(len(qcd_data_faillep))]
    qcdratio_wlnu_faillep = [np.sum(wlnu_data_faillep[ix])/np.sum(qcd_data_faillep[ix]) for ix in range(len(qcd_data_faillep))]
    #print('faillep ratio sig ',qcdratio_sig_faillep)
    #print('faillep ratio top ',qcdratio_top_faillep)
    #print('faillep ratio wlnu',qcdratio_wlnu_faillep)

    qcd_data_lowmet = [getQCDfromData(qcd_cr_lowmet,"fail",hptslice=slice(qcd_cr_lowmet.axis('h_pt').edges()[ix-1] if ix>1 else None,qcd_cr_lowmet.axis('h_pt').edges()[ix] if ix<len(qcd_cr_lowmet.axis('h_pt').edges())-1 else None))[0] for ix in range(len(qcd_cr_lowmet.axis('h_pt').edges()))]
    top_data_lowmet = [getQCDfromData(top_cr_lowmet,"fail",hptslice=slice(top_cr_lowmet.axis('h_pt').edges()[ix-1] if ix>1 else None,top_cr_lowmet.axis('h_pt').edges()[ix] if ix<len(top_cr_lowmet.axis('h_pt').edges())-1 else None))[0] for ix in range(len(top_cr_lowmet.axis('h_pt').edges()))]
    wlnu_data_lowmet = [getQCDfromData(wlnu_cr_lowmet,"fail",hptslice=slice(wlnu_cr_lowmet.axis('h_pt').edges()[ix-1] if ix>1 else None,wlnu_cr_lowmet.axis('h_pt').edges()[ix] if ix<len(wlnu_cr_lowmet.axis('h_pt').edges())-1 else None))[0] for ix in range(len(wlnu_cr_lowmet.axis('h_pt').edges()))]
    sig_data_lowmet = [getQCDfromData(sig_lowmet,"fail",hptslice=slice(sig_lowmet.axis('h_pt').edges()[ix-1] if ix>1 else None,sig_lowmet.axis('h_pt').edges()[ix] if ix<len(sig_lowmet.axis('h_pt').edges())-1 else None))[0] for ix in range(len(sig_lowmet.axis('h_pt').edges()))]

    #qcdratio_sig_lowmet = [sig_data_lowmet[ix]/qcd_data_lowmet[ix] for ix in range(len(qcd_data_lowmet))]
    #qcdratio_top_lowmet = [top_data_lowmet[ix]/qcd_data_lowmet[ix] for ix in range(len(qcd_data_lowmet))]
    #qcdratio_wlnu_lowmet = [wlnu_data_lowmet[ix]/qcd_data_lowmet[ix] for ix in range(len(qcd_data_lowmet))]
    qcdratio_sig_lowmet = [np.sum(sig_data_lowmet[ix])/np.sum(qcd_data_lowmet[ix]) for ix in range(len(qcd_data_lowmet))]
    qcdratio_top_lowmet = [np.sum(top_data_lowmet[ix])/np.sum(qcd_data_lowmet[ix]) for ix in range(len(qcd_data_lowmet))]
    qcdratio_wlnu_lowmet = [np.sum(wlnu_data_lowmet[ix])/np.sum(qcd_data_lowmet[ix]) for ix in range(len(qcd_data_lowmet))]
    #print('lowmet ratio sig ',qcdratio_sig_lowmet)
    #print('lowmet ratio top ',qcdratio_top_lowmet)
    #print('lowmet ratio wlnu',qcdratio_wlnu_lowmet)

    qcdratio_sig = [(qcdratio_sig_lowmet[ix]+qcdratio_sig_faillep[ix])/2. for ix in range(len(qcdratio_sig_faillep))]
    qcdratio_sig_err = [np.abs(qcdratio_sig_lowmet[ix]-qcdratio_sig[ix])/2. for ix in range(len(qcdratio_sig_faillep))]
    qcdratio_top = [(qcdratio_top_lowmet[ix]+qcdratio_top_faillep[ix])/2. for ix in range(len(qcdratio_top_faillep))]
    qcdratio_top_err = [np.abs(qcdratio_top_lowmet[ix]-qcdratio_top[ix])/2. for ix in range(len(qcdratio_top_faillep))]
    qcdratio_wlnu = [(qcdratio_wlnu_lowmet[ix]+qcdratio_wlnu_faillep[ix])/2. for ix in range(len(qcdratio_wlnu_faillep))]
    qcdratio_wlnu_err = [np.abs(qcdratio_wlnu_lowmet[ix]-qcdratio_wlnu[ix])/2. for ix in range(len(qcdratio_wlnu_faillep))]

    #print('qcdratio_sig',qcdratio_sig, qcdratio_sig_err)
    #print('qcdratio_top',qcdratio_top, qcdratio_top_err)
    #print('qcdratio_wlnu',qcdratio_wlnu, qcdratio_wlnu_err)

    qcdfail_temp = {}
    qcdfail_temp_dn = {}
    qcdfail_temp_up = {}

    for region in ['pass','loosepass','fail']:
        qcdfail_temp[region], qcdfail_temp_dn[region], qcdfail_temp_up[region] = getQCDfromData(qcd_cr_faillep,region,default=0.)
        qcdfail_temp[region] = np.array([np.sum(qcdfail_temp[region])])
        qcdfail_temp_dn[region] = np.array([np.sum(qcdfail_temp_dn[region])])
        qcdfail_temp_up[region] = np.array([np.sum(qcdfail_temp_up[region])])

    qcdratio_F = {
        'pass':qcdfail_temp['pass']/qcdfail_temp['fail'],
        'loosepass':qcdfail_temp['loosepass']/qcdfail_temp['fail'],
        'fail':np.ones_like(qcdfail_temp['fail']),
    }
    qcdratio_F_dn = {
        'pass':qcdfail_temp_dn['pass']/qcdfail_temp_dn['fail'],
        'loosepass':qcdfail_temp_dn['loosepass']/qcdfail_temp_dn['fail'],
        'fail':np.ones_like(qcdfail_temp['fail']),
    }
    qcdratio_F_up = {
        'pass':qcdfail_temp_up['pass']/qcdfail_temp_up['fail'],
        'loosepass':qcdfail_temp_up['loosepass']/qcdfail_temp_up['fail'],
        'fail':np.ones_like(qcdfail_temp['fail']),
    }
    #for region in ['pass', 'loosepass', 'fail']:
    #    qcdratio_F[region] = np.sum(qcdratio_F[region]*qcd_from_data)/np.sum(qcd_from_data)
    #    qcdratio_F_dn[region] = np.sum(qcdratio_F_dn[region]*qcd_from_data)/np.sum(qcd_from_data)
    #    qcdratio_F_up[region] = np.sum(qcdratio_F_up[region]*qcd_from_data)/np.sum(qcd_from_data)
    #print('qcdratio_F',qcdratio_F)
    #print('qcdratio_F_dn',qcdratio_F_dn)
    #print('qcdratio_F_up',qcdratio_F_up)

    qcd_from_data = [getQCDfromData(qcd_cr_hist,"fail",hptslice=slice(qcd_cr_hist.axis('h_pt').edges()[ix-1] if ix>1 else None,qcd_cr_hist.axis('h_pt').edges()[ix] if ix<len(qcd_cr_hist.axis('h_pt').edges())-1 else None))[0] for ix in range(len(qcd_cr_hist.axis('h_pt').edges()))]
    #print('qcd_from_data',qcd_from_data)
    qcdedges = [qcd_cr_hist.axis('h_pt').edges()[ix] for ix in range(len(qcd_cr_hist.axis('h_pt').edges()))]
    qcdbins = [[ix for ix in range(len(qcdedges)-1) if qcdedges[ix]>=(ptbins[ipt] if ptbins[ipt] is not None else -999999.) and qcdedges[ix+1]<=(ptbins[ipt+1] if ptbins[ipt+1] is not None else 999999.)] for ipt in range(npt)]

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
        for region in ['pass','loosepass','fail']:
            ch = rl.Channel("ptbin%d%s" % (ptbin, region))
            model.addChannel(ch)

            isPass = region=='pass'
            isLoosePass = region=='loosepass'

            thehist = sig_hist
    
            for sName in thehist.identifiers('sample'):
                if sName.name=='ignore':
                    continue
                if sName.name=='data_obs' and isPass and not usingData:
                    continue
                templ = (intRegion(thehist[sName],region,hptslice=ptslices[ptbin])[1:-1], mtt.binning, mtt.name)
                for iv,val in enumerate(templ[0]):
                    if val<0.:
                        templ[0][iv] = 0.
                if sName.name=='data_obs':
                    ch.setObservation(templ)
                else:
                    if sName.name=='multijet':
                        qcdpred = np.sum(np.stack([qcd_from_data[ix] * qcdratio_sig[ix] * np.mean(qcdratio_F[region]) for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else [np.zeros_like(qcd_from_data[0])]),axis=0)
                        qcdpred_dn = np.sum(np.stack([qcd_from_data[ix] * (qcdratio_sig[ix] - qcdratio_sig_err[ix]) * np.mean(qcdratio_F[region]) for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else [np.zeros_like(qcd_from_data[0])]),axis=0)
                        qcdpred_up = np.sum(np.stack([qcd_from_data[ix] * (qcdratio_sig[ix] + qcdratio_sig_err[ix]) * np.mean(qcdratio_F[region]) for ix in qcdbins[ptbin]] if len(qcdbins[ptbin])>0 else [np.zeros_like(qcd_from_data[0])]),axis=0)
                        #print('qcdpred',qcdpred)
                        #print('qcdpred_dn',qcdpred_dn)
                        #print('qcdpred_up',qcdpred_up)
                        templ = (qcdpred, mtt.binning, mtt.name)
                    stype = rl.Sample.SIGNAL if sName.name == 'htt125' else rl.Sample.BACKGROUND
                    sample = rl.TemplateSample(ch.name + '_' + sName.name, stype, templ)
        
                    #sample.setParamEffect(jec, 1.05)
                    if sName.name in syst_dict:
                        nom = intRegion(thehist[sName],region,hptslice=ptslices[ptbin])[1:-1]
                        for syst in syst_dict[sName.name]:
                            syst_params = syst_dict[sName.name][syst]
                            dnmod = intRegion(thehist[sName],region,systematic=syst_params[1],hptslice=ptslices[ptbin])[1:-1]
                            upmod = intRegion(thehist[sName],region,systematic=syst_params[2],hptslice=ptslices[ptbin])[1:-1]
                            sample.setParamEffect(syst_params[0], np.divide(dnmod, nom, out=np.ones_like(nom), where=nom>0.) if (dnmod!=nom).all() else np.ones_like(nom)*0.999, np.divide(upmod, nom, out=np.ones_like(nom), where=nom>0.) if (upmod!=nom).all() else np.ones_like(nom)*1.001)

                    #if sName.name=='top':
                    #    sample.setParamEffect(top_norm, 1.05)
                    #if sName.name=='wlnu':
                    #    sample.setParamEffect(wlnu_norm, 1.10)
                    if sName.name=='multijet':
                        sample.setParamEffect(qcd_norm, 1.20)
                        sample.setParamEffect(qcd_pass if isPass else qcd_loosepass if isLoosePass else qcd_fail, np.divide(qcdpred_dn, qcdpred, out=np.ones_like(qcdpred), where=qcdpred>0.), np.divide(qcdpred_up, qcdpred, out=np.ones_like(qcdpred), where=qcdpred>0.))
                    if sName.name=='vvqq':
                        sample.setParamEffect(vqq_norm, 1.10)
                    #if sName.name=='ztt':
                    #    sample.setParamEffect(ztt_norm, 1.05)
                    if sName.name=='htt125' or sName.name=='ztt':
                        nom_full = intRegion(thehist[sName],region,hptslice=ptslices[ptbin])
                        shift_dn = nom_full[2:]
                        shift_up = nom_full[:-2]
                        #sample.setParamEffect(m_scale, np.divide(shift_dn, nom_full[1:-1], out=np.ones_like(nom_full[1:-1]), where=nom_full[1:-1]>0.), np.divide(shift_up, nom_full[1:-1], out=np.ones_like(nom_full[1:-1]), where=nom_full[1:-1]>0.))
                        #print(leptype,'htt125 ptbin',ptbin,region,nom_full[1:-1])
                    sample.setParamEffect(trig, 1.02)
                    sample.setParamEffect(lumi, 1.024)
            
                    ch.addSample(sample)
            if not usingData and isPass:
                ch.setObservation((np.zeros(len(mtt.binning)-1),mtt.binning, mtt.name))

    h125vals = [intRegion(thehist['htt125'],'pass',hptslice=ptslices[ptbin])[1:-1] for ptbin in range(npt)]
    zttvals = [intRegion(thehist['ztt'],'pass',hptslice=ptslices[ptbin])[1:-1] for ptbin in range(npt)]
    wlnuvals = [intRegion(thehist['wlnu'],'pass',hptslice=ptslices[ptbin])[1:-1] for ptbin in range(npt)]
    probs = [np.prod([math.erfc(h125vals[ipt][ib]/(math.sqrt(abs(zttvals[ipt][ib]+wlnuvals[ipt][ib])) if zttvals[ipt][ib]+wlnuvals[ipt][ib]>2. else 2.)) for ib in range(len(h125vals[ipt]))]) for ipt in range(npt)]
    print(leptype,'PROB',np.prod(probs),'->',erfinv(1.-np.prod(probs)),'\t\t',probs)

    for ptbin in range(npt):
        failCh = model['ptbin%dfail' % ptbin]
        passCh = model['ptbin%dpass' % ptbin]
        loosePassCh = model['ptbin%dloosepass' % ptbin]

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
        toppass.setParamEffect(topnormSF, 1*topnormSF)
        toploosepass.setParamEffect(topnormSF, 1*topnormSF)
        topfail.setParamEffect(topnormSF, 1*topnormSF)

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
        wlnupass.setParamEffect(wlnunormSF, 1*wlnunormSF)
        wlnuloosepass.setParamEffect(wlnunormSF, 1*wlnunormSF)
        wlnufail.setParamEffect(wlnunormSF, 1*wlnunormSF)

        zttpass = passCh['ztt']
        zttloosepass = loosePassCh['ztt']
        zttfail = failCh['ztt']
        zttLP = zttpass.getExpectation(nominal=True).sum() / zttloosepass.getExpectation(nominal=True).sum()
        httLP = passCh['htt125'].getExpectation(nominal=True).sum() / loosePassCh['htt125'].getExpectation(nominal=True).sum()
        zttpass.setParamEffect(ztt_eff, 1*ztt_eff)
        zttloosepass.setParamEffect(ztt_eff, (1 - ztt_eff) * zttLP + 1)
        passCh['htt125'].setParamEffect(ztt_eff, 1*ztt_eff)
        loosePassCh['htt125'].setParamEffect(ztt_eff, (1 - ztt_eff) * httLP + 1)
        zttpass.setParamEffect(rztt, 1*rztt)
        zttloosepass.setParamEffect(rztt, 1*rztt)
        zttfail.setParamEffect(rztt, 1*rztt)
        #zttpass.setParamEffect(rztt, 1.05)
        #zttloosepass.setParamEffect(rztt, 1.05)
        #zttfail.setParamEffect(rztt, 1.05)


    # Fill in top CR
    for region in ['pass', 'loosepass', 'fail']:
        ch = rl.Channel("topCR%s" % (region, ))
        model.addChannel(ch)

        isPass = region=='pass'
        isLoosePass = region=='loosepass'

        thehist = top_cr_hist

        for sName in thehist.identifiers('sample'):
            if sName.name=='ignore':
                continue
            if singleBinCR:
                templ = (np.array([np.sum(intRegion(thehist[sName],region)[1:-1])]), mttone.binning, mttone.name)
            else:
                templ = (intRegion(thehist[sName],region)[1:-1], mtt.binning, mtt.name)
            for iv,val in enumerate(templ[0]):
                if val<0.:
                    templ[0][iv] = 0.
            if sName.name=='data_obs':
                ch.setObservation(templ)
            else:
                if sName.name=='multijet':
                    qcdpred = np.sum(np.stack([qcd_from_data[ix] * qcdratio_top[ix] * np.mean(qcdratio_F[region]) for ix in range(len(qcd_from_data))]),axis=0)
                    qcdpred_dn = np.sum(np.stack([qcd_from_data[ix] * (qcdratio_top[ix] - qcdratio_top_err[ix]) * np.mean(qcdratio_F[region]) for ix in range(len(qcd_from_data))]),axis=0)
                    qcdpred_up = np.sum(np.stack([qcd_from_data[ix] * (qcdratio_top[ix] + qcdratio_top_err[ix]) * np.mean(qcdratio_F[region]) for ix in range(len(qcd_from_data))]),axis=0)
                    if singleBinCR:
                        templ = (np.array([np.sum(qcdpred)]), mttone.binning, mttone.name)
                    else:
                        templ = (qcdpred, mtt.binning, mtt.name)
                stype = rl.Sample.SIGNAL if sName.name == 'htt125' else rl.Sample.BACKGROUND
                sample = rl.TemplateSample(ch.name + '_' + sName.name, stype, templ)
    
                #sample.setParamEffect(jec, 1.05)
                if sName.name in syst_dict:
                    nom_base = intRegion(thehist[sName],region)[1:-1]
                    for syst in syst_dict[sName.name]:
                        syst_params = syst_dict[sName.name][syst]
                        dnmod = intRegion(thehist[sName],region,systematic=syst_params[1])[1:-1]
                        upmod = intRegion(thehist[sName],region,systematic=syst_params[2])[1:-1]
                        if singleBinCR:
                            dnmod = np.array([np.sum(dnmod)])
                            nom = np.array([np.sum(nom_base)])
                            upmod = np.array([np.sum(upmod)])
                        else:
                            nom = nom_base
                        sample.setParamEffect(syst_params[0], np.divide(dnmod, nom, out=np.ones_like(nom), where=nom>0.) if (dnmod!=nom).all() else np.ones_like(nom)*0.999, np.divide(upmod, nom, out=np.ones_like(nom), where=nom>0.) if (upmod!=nom).all() else np.ones_like(nom)*1.001)

                #if sName.name=='top':
                #    sample.setParamEffect(top_norm, 1.05)
                #if sName.name=='wlnu':
                    #sample.setParamEffect(wlnu_norm, 1.10)
                if sName.name=='multijet':
                    sample.setParamEffect(qcd_norm, 1.20)
                    if singleBinCR:
                        qcdpred_dn = np.array([np.sum(qcdpred_dn)])
                        qcdpred = np.array([np.sum(qcdpred)])
                        qcdpred_up = np.array([np.sum(qcdpred_up)])
                    sample.setParamEffect(qcdtop_pass if isPass else qcdtop_loosepass if isLoosePass else qcdtop_fail, np.divide(qcdpred_dn, qcdpred, out=np.ones_like(qcdpred), where=qcdpred>0.), np.divide(qcdpred_up, qcdpred, out=np.ones_like(qcdpred), where=qcdpred>0.))
                if sName.name=='vvqq':
                    sample.setParamEffect(vqq_norm, 1.10)
                #if sName.name=='ztt':
                #    sample.setParamEffect(ztt_norm, 1.10)
                if sName.name=='htt125' or sName.name=='ztt':
                    nom_full = intRegion(thehist[sName],region)
                    shift_dn = nom_full[2:]
                    shift_up = nom_full[:-2]
                    #sample.setParamEffect(m_scale, np.divide(shift_dn, nom_full[1:-1], out=np.ones_like(nom_full[1:-1]), where=nom_full[1:-1]>0.), np.divide(shift_up, nom_full[1:-1], out=np.ones_like(nom_full[1:-1]), where=nom_full[1:-1]>0.))
                sample.setParamEffect(trig, 1.02)
                sample.setParamEffect(lumi, 1.024)
        
                ch.addSample(sample)

    qcdpass = model['topCRpass']['multijet']
    qcdloosepass = model['topCRloosepass']['multijet']
    qcdfail = model['topCRfail']['multijet']
    qcdpass.setParamEffect(qcdnormSF_top, 1.20)
    qcdloosepass.setParamEffect(qcdnormSF_top, 1.20)
    qcdfail.setParamEffect(qcdnormSF_top, 1.20)
    #qcdpass.setParamEffect(qcdnormSF_top, 1*qcdnormSF_top)
    #qcdloosepass.setParamEffect(qcdnormSF_top, 1*qcdnormSF_top)
    #qcdfail.setParamEffect(qcdnormSF_top, 1*qcdnormSF_top)

    toppass = model['topCRpass']['top']
    toploosepass = model['topCRloosepass']['top']
    topfail = model['topCRfail']['top']
    topPF = (toploosepass.getExpectation(nominal=True).sum() + toppass.getExpectation(nominal=True).sum()) / topfail.getExpectation(nominal=True).sum()
    topLPF = toppass.getExpectation(nominal=True).sum() / toploosepass.getExpectation(nominal=True).sum()
    topRLPF = 1./(1.+(1./topLPF)) # P/(L+P)
    toppass.setParamEffect(topLeffSF, 1*topLeffSF)
    toploosepass.setParamEffect(topLeffSF, (1 - topLeffSF) * topLPF + 1)
    toploosepass.setParamEffect(topeffSF, 1*topeffSF)
    topfail.setParamEffect(topeffSF, (1 - topeffSF) * (1-topRLPF) * topPF + 1)
    toppass.setParamEffect(topnormSF, 1*topnormSF)
    toploosepass.setParamEffect(topnormSF, 1*topnormSF)
    topfail.setParamEffect(topnormSF, 1*topnormSF)

    wlnupass = model['topCRpass']['wlnu']
    wlnuloosepass = model['topCRloosepass']['wlnu']
    wlnufail = model['topCRfail']['wlnu']
    wlnuPF = (wlnuloosepass.getExpectation(nominal=True).sum() + wlnupass.getExpectation(nominal=True).sum()) / wlnufail.getExpectation(nominal=True).sum()
    wlnuLPF = wlnupass.getExpectation(nominal=True).sum() / wlnuloosepass.getExpectation(nominal=True).sum()
    wlnuRLPF = 1./(1.+(1./wlnuLPF)) # P/(L+P)
    wlnupass.setParamEffect(wlnuLeffSF, 1*wlnuLeffSF)
    wlnuloosepass.setParamEffect(wlnuLeffSF, (1 - wlnuLeffSF) * wlnuLPF + 1)
    wlnuloosepass.setParamEffect(wlnueffSF, 1*wlnueffSF)
    wlnufail.setParamEffect(wlnueffSF, (1 - wlnueffSF) * (1-wlnuRLPF) * wlnuPF + 1)
    wlnupass.setParamEffect(wlnunormSF, 1*wlnunormSF)
    wlnuloosepass.setParamEffect(wlnunormSF, 1*wlnunormSF)
    wlnufail.setParamEffect(wlnunormSF, 1*wlnunormSF)

    zttpass = model['topCRpass']['ztt']
    zttloosepass = model['topCRloosepass']['ztt']
    zttfail = model['topCRfail']['ztt']
    zttLP = zttpass.getExpectation(nominal=True).sum() / zttloosepass.getExpectation(nominal=True).sum()
    httLP = model['topCRpass']['htt125'].getExpectation(nominal=True).sum() / model['topCRloosepass']['htt125'].getExpectation(nominal=True).sum()
    zttpass.setParamEffect(ztt_eff, 1*ztt_eff)
    zttloosepass.setParamEffect(ztt_eff, (1 - ztt_eff)* zttLP + 1)
    model['topCRpass']['htt125'].setParamEffect(ztt_eff, 1*ztt_eff)
    model['topCRloosepass']['htt125'].setParamEffect(ztt_eff, (1 - ztt_eff) * httLP + 1)
    zttpass.setParamEffect(rztt, 1*rztt)
    zttloosepass.setParamEffect(rztt, 1*rztt)
    zttfail.setParamEffect(rztt, 1*rztt)
    #zttpass.setParamEffect(rztt, 1.05)
    #zttloosepass.setParamEffect(rztt, 1.05)
    #zttfail.setParamEffect(rztt, 1.05)

    # Fill in wlnu CR
    for region in ['pass', 'fail', 'loosepass']:
        ch = rl.Channel("wlnuCR%s" % (region, ))
        model.addChannel(ch)

        isPass = region=='pass'
        isLoosePass = region=='loosepass'

        thehist = wlnu_cr_hist

        for sName in thehist.identifiers('sample'):
            if sName.name=='ignore':
                continue
            if singleBinCR:
                templ = (np.array([np.sum(intRegion(thehist[sName],region)[1:-1])]), mttone.binning, mttone.name)
            else:
                templ = (intRegion(thehist[sName],region)[1:-1], mtt.binning, mtt.name)
            for iv,val in enumerate(templ[0]):
                if val<0.:
                    templ[0][iv] = 0.
            if sName.name=='data_obs':
                ch.setObservation(templ)
            else:
                if sName.name=='multijet':
                    qcdpred = np.sum(np.stack([qcd_from_data[ix] * qcdratio_wlnu[ix] * np.mean(qcdratio_F[region]) for ix in range(len(qcd_from_data))]),axis=0)
                    qcdpred_dn = np.sum(np.stack([qcd_from_data[ix] * (qcdratio_wlnu[ix] - qcdratio_wlnu_err[ix]) * np.mean(qcdratio_F[region]) for ix in range(len(qcd_from_data))]),axis=0)
                    qcdpred_up = np.sum(np.stack([qcd_from_data[ix] * (qcdratio_wlnu[ix] + qcdratio_wlnu_err[ix]) * np.mean(qcdratio_F[region]) for ix in range(len(qcd_from_data))]),axis=0)
                    if singleBinCR:
                        templ = (np.array([np.sum(qcdpred)]), mttone.binning, mttone.name)
                    else:
                        templ = (qcdpred, mtt.binning, mtt.name)
                stype = rl.Sample.SIGNAL if sName.name == 'htt125' else rl.Sample.BACKGROUND
                sample = rl.TemplateSample(ch.name + '_' + sName.name, stype, templ)
    
                #sample.setParamEffect(jec, 1.05)
                if sName.name in syst_dict:
                    nom_base = intRegion(thehist[sName],region)[1:-1]
                    for syst in syst_dict[sName.name]:
                        syst_params = syst_dict[sName.name][syst]
                        dnmod = intRegion(thehist[sName],region,systematic=syst_params[1])[1:-1]
                        upmod = intRegion(thehist[sName],region,systematic=syst_params[2])[1:-1]
                        if singleBinCR:
                            dnmod = np.array([np.sum(dnmod)])
                            nom = np.array([np.sum(nom_base)])
                            upmod = np.array([np.sum(upmod)])
                        else:
                            nom = nom_base
                        sample.setParamEffect(syst_params[0], np.divide(dnmod, nom, out=np.ones_like(nom), where=nom>0.) if (dnmod!=nom).all() else np.ones_like(nom)*0.999, np.divide(upmod, nom, out=np.ones_like(nom), where=nom>0.) if (upmod!=nom).all() else np.ones_like(nom)*1.001)

                #if sName.name=='top':
                #    sample.setParamEffect(top_norm, 1.05)
                #if sName.name=='wlnu':
                #    sample.setParamEffect(wlnu_norm, 1.10)
                if sName.name=='multijet':
                    sample.setParamEffect(qcd_norm, 1.20)
                    if singleBinCR:
                        qcdpred_dn = np.array([np.sum(qcdpred_dn)])
                        qcdpred = np.array([np.sum(qcdpred)])
                        qcdpred_up = np.array([np.sum(qcdpred_up)])
                    sample.setParamEffect(qcdwlnu_pass if isPass else qcdwlnu_loosepass if isLoosePass else qcdwlnu_fail, np.divide(qcdpred_dn, qcdpred, out=np.ones_like(qcdpred), where=qcdpred>0.), np.divide(qcdpred_up, qcdpred, out=np.ones_like(qcdpred), where=qcdpred>0.))
                if sName.name=='vvqq':
                    sample.setParamEffect(vqq_norm, 1.10)
                #if sName.name=='ztt':
                #    sample.setParamEffect(ztt_norm, 1.10)
                if sName.name=='htt125' or sName.name=='ztt':
                    nom_full = intRegion(thehist[sName],region)
                    shift_dn = nom_full[2:]
                    shift_up = nom_full[:-2]
                    #sample.setParamEffect(m_scale, np.divide(shift_dn, nom_full[1:-1], out=np.ones_like(nom_full[1:-1]), where=nom_full[1:-1]>0.), np.divide(shift_up, nom_full[1:-1], out=np.ones_like(nom_full[1:-1]), where=nom_full[1:-1]>0.))
                sample.setParamEffect(trig, 1.02)
                sample.setParamEffect(lumi, 1.024)
        
                ch.addSample(sample)

    qcdpass = model['wlnuCRpass']['multijet']
    qcdloosepass = model['wlnuCRloosepass']['multijet']
    qcdfail = model['wlnuCRfail']['multijet']
    qcdpass.setParamEffect(qcdnormSF_wlnu, 1.20)
    qcdloosepass.setParamEffect(qcdnormSF_wlnu, 1.20)
    qcdfail.setParamEffect(qcdnormSF_wlnu, 1.20)
    #qcdpass.setParamEffect(qcdnormSF_wlnu, 1*qcdnormSF_wlnu)
    #qcdloosepass.setParamEffect(qcdnormSF_wlnu, 1*qcdnormSF_wlnu)
    #qcdfail.setParamEffect(qcdnormSF_wlnu, 1*qcdnormSF_wlnu)

    toppass = model['wlnuCRpass']['top']
    toploosepass = model['wlnuCRloosepass']['top']
    topfail = model['wlnuCRfail']['top']
    topPF = (toploosepass.getExpectation(nominal=True).sum() + toppass.getExpectation(nominal=True).sum()) / topfail.getExpectation(nominal=True).sum()
    topLPF = toppass.getExpectation(nominal=True).sum() / toploosepass.getExpectation(nominal=True).sum()
    topRLPF = 1./(1.+(1./topLPF)) # P/(L+P)
    toppass.setParamEffect(topLeffSF, 1*topLeffSF)
    toploosepass.setParamEffect(topLeffSF, (1 - topLeffSF) * topLPF + 1)
    toploosepass.setParamEffect(topeffSF, 1*topeffSF)
    topfail.setParamEffect(topeffSF, (1 - topeffSF) * (1-topRLPF) * topPF + 1)
    toppass.setParamEffect(topnormSF, 1*topnormSF)
    toploosepass.setParamEffect(topnormSF, 1*topnormSF)
    topfail.setParamEffect(topnormSF, 1*topnormSF)

    wlnupass = model['wlnuCRpass']['wlnu']
    wlnuloosepass = model['wlnuCRloosepass']['wlnu']
    wlnufail = model['wlnuCRfail']['wlnu']
    wlnuPF = (wlnuloosepass.getExpectation(nominal=True).sum() + wlnupass.getExpectation(nominal=True).sum()) / wlnufail.getExpectation(nominal=True).sum()
    wlnuLPF = wlnupass.getExpectation(nominal=True).sum() / wlnuloosepass.getExpectation(nominal=True).sum()
    wlnuRLPF = 1./(1.+(1./wlnuLPF)) # P/(L+P)
    wlnupass.setParamEffect(wlnuLeffSF, 1*wlnuLeffSF)
    wlnuloosepass.setParamEffect(wlnuLeffSF, (1 - wlnuLeffSF) * wlnuLPF + 1)
    wlnuloosepass.setParamEffect(wlnueffSF, 1*wlnueffSF)
    wlnufail.setParamEffect(wlnueffSF, (1 - wlnueffSF) * (1-wlnuRLPF) * wlnuPF + 1)
    wlnupass.setParamEffect(wlnunormSF, 1*wlnunormSF)
    wlnuloosepass.setParamEffect(wlnunormSF, 1*wlnunormSF)
    wlnufail.setParamEffect(wlnunormSF, 1*wlnunormSF)

    zttpass = model['wlnuCRpass']['ztt']
    zttloosepass = model['wlnuCRloosepass']['ztt']
    zttfail = model['wlnuCRfail']['ztt']
    zttLP = zttpass.getExpectation(nominal=True).sum() / zttloosepass.getExpectation(nominal=True).sum()
    httLP = model['wlnuCRpass']['htt125'].getExpectation(nominal=True).sum() / model['wlnuCRloosepass']['htt125'].getExpectation(nominal=True).sum()
    zttpass.setParamEffect(ztt_eff, 1*ztt_eff)
    zttloosepass.setParamEffect(ztt_eff, (1 - ztt_eff) * zttLP + 1)
    model['wlnuCRpass']['htt125'].setParamEffect(ztt_eff, 1*ztt_eff)
    model['wlnuCRloosepass']['htt125'].setParamEffect(ztt_eff, (1 - ztt_eff) * httLP + 1)
    zttpass.setParamEffect(rztt, 1*rztt)
    zttloosepass.setParamEffect(rztt, 1*rztt)
    zttfail.setParamEffect(rztt, 1*rztt)
    #zttpass.setParamEffect(rztt, 1.05)
    #zttloosepass.setParamEffect(rztt, 1.05)
    #zttfail.setParamEffect(rztt, 1.05)

    with open(os.path.join("%s/%s"%(str(tmpdir),label), 'had%sModel.pkl'%leptype), "wb") as fout:
        pickle.dump(model, fout, protocol=2)

    model.renderCombine(os.path.join("%s/%s"%(str(tmpdir),label), 'had%sModel'%leptype))

def createHadHad(sig_hist, sig_met_hist, top_cr_hist, wlnu_cr_hist, qcd_cr_hist, var_name, mttbins, ptbins, tmpdir, label, usingData, nnCut_jet, nnCut_jet_loose, nnCut_met, nnCut_met_loose, metCut, h_pt_min):

    #could be made region-dependent
    samp_combinations = {
        'data_obs'   : ['data'],
        'top'        : ['tt-had', 'tt-semilep', 'tt-dilep', 'st'],
        'htt125'     : ['h125'],
        'multijet'   : ['qcd'],
        'ztt'        : ['zll'],
        'wlnu'       : ['wjets'],
        'vqq'        : ['vv'],
        'ignore'     : ['vqq'],
    }

    sig_jet_invdphi = sig_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('jetmet_dphi',slice(1.6,None),'over')
    sig_met_invdphi = sig_met_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('jetmet_dphi',slice(1.6,None),'over')
    top_cr_invdphi = top_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('jetmet_dphi',slice(1.6,None),'over')
    wlnu_cr_invdphi = wlnu_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('jetmet_dphi',slice(1.6,None),'over')
    qcd_cr_invdphi = qcd_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('jetmet_dphi',slice(1.6,None),'over')

    sig_jet_hist = sig_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('jetmet_dphi',slice(None,1.6),'under')
    sig_met_hist = sig_met_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('jetmet_dphi',slice(None,1.6),'under')
    top_cr_hist = top_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('jetmet_dphi',slice(None,1.6),'under')
    wlnu_cr_hist = wlnu_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('jetmet_dphi',slice(None,1.6),'under')
    qcd_cr_hist = qcd_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('jetmet_dphi',slice(None,1.6),'under')

    #jec = rl.NuisanceParameter('CMS_jec', 'lnN')
    top_norm = rl.NuisanceParameter('CMS_top_norm', 'lnN')
    wlnu_norm = rl.NuisanceParameter('CMS_wlnu_norm', 'lnN')
    qcd_norm = rl.NuisanceParameter('CMS_qcd_norm', 'lnN')
    qcd_norm_top = rl.NuisanceParameter('CMS_qcd_norm_top', 'lnN')
    qcd_norm_wlnu = rl.NuisanceParameter('CMS_qcd_norm_wlnu', 'lnN')
    #ztt_norm = rl.NuisanceParameter('CMS_ztt_norm', 'lnN')
    vqq_norm = rl.NuisanceParameter('CMS_vqq_norm', 'lnN')
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
    qcdnormSF_top = rl.NuisanceParameter('qcdnormSF_top_hadhad', 'lnN')
    qcdwlnu_fail = rl.NuisanceParameter('qcd_wlnu_Rfail_hadhad', 'shape')
    qcdwlnu_loosepass = rl.NuisanceParameter('qcd_wlnu_Rloosepass_hadhad', 'shape')
    qcdwlnu_pass = rl.NuisanceParameter('qcd_wlnu_Rpass_hadhad', 'shape')
    #qcdnormSF_wlnu = rl.IndependentParameter('qcdnormSF_wlnu_hadhad', 1., 0, 10)
    qcdnormSF_wlnu = rl.NuisanceParameter('qcdnormSF_wlnu_hadhad', 'lnN')

    #m_scale = rl.NuisanceParameter('massscale_hadhad', 'shape')

    topeffSF = rl.IndependentParameter('topeffSF_hadhad', 1., 0, 10)
    topnormSF = rl.IndependentParameter('topnormSF_hadhad', 1., 0, 10)
    wlnueffSF = rl.IndependentParameter('wlnueffSF_hadhad', 1., 0, 10)
    wlnunormSF = rl.IndependentParameter('wlnunormSF_hadhad', 1., 0, 10)

    topLeffSF = rl.IndependentParameter('topLeffSF_hadhad', 1., 0, 10)
    wlnuLeffSF = rl.IndependentParameter('wlnuLeffSF_hadhad', 1., 0, 10)

    rztt = rl.IndependentParameter('r_ztt_hadhad', 1., 0, 10)
    #rztt = rl.NuisanceParameter('r_ztt_hadhad', 'lnN')
    ztt_eff = rl.IndependentParameter('ztt_eff_hadhad', 1., 0, 10)

    toppt = rl.NuisanceParameter('toppt', 'shape')
    syst_dict = {
        "top":{
            'toppt':[toppt,'nominal','TopPtReweightUp'],
        },
    }

    npt = len(ptbins) - 1
    ptslices = [slice(ptbins[ipt],ptbins[ipt+1]) for ipt in range(npt)]
    mtt = rl.Observable(var_name, mttbins[1:-1])
    mttone = rl.Observable(var_name+'_one', np.array([mttbins[1],mttbins[-2]]))

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
            the_int = the_int.integrate('sample',samplelist).values()[()]
        else:
            the_int = the_int.sum('sample').values()[()]

        if debug:
            print('\tdebug',the_int)

        return the_int

    def getQCDfromData(inhist,region,default=1.,systematic='nominal',hptslice=slice(h_pt_min,None),mslice=None):
        qcd_data = intRegion(inhist['data_obs'],region,systematic=systematic,hptslice=hptslice,mslice=mslice)[1:-1]
        other_mc = intRegion(inhist,region,samplelist=[s.name for s in inhist.identifiers('sample') if s.name!='data_obs' and s.name!='multijet'],systematic=systematic,hptslice=hptslice,mslice=mslice)[1:-1]
        qcd_data = qcd_data - other_mc
        qcd_data = np.array([v if v>0. else default for v in qcd_data])
        qcd_data_int = hist.poisson_interval(qcd_data,qcd_data)
        qcd_temp = qcd_data
        qcd_temp_dn = np.array([qcd_data_int[1][bi] if qcd_data_int[1][bi] > 0. and not np.isnan(qcd_data_int[1][bi]) else default for bi in range(len(qcd_data))])
        qcd_temp_up = np.array([qcd_data_int[0][bi] if qcd_data_int[0][bi] > 0. and not np.isnan(qcd_data_int[0][bi]) else default for bi in range(len(qcd_data))])
        return qcd_temp,qcd_temp_dn,qcd_temp_up

    def getQCDfromDataFullMass(inhist,region,default=1.,systematic='nominal',hptslice=slice(h_pt_min,None),mslice=None):
        qcd_data = intRegion(inhist['data_obs'],region,systematic=systematic,hptslice=hptslice,mslice=mslice,mrebin=None)
        other_mc = intRegion(inhist,region,samplelist=[s.name for s in inhist.identifiers('sample') if s.name!='data_obs' and s.name!='multijet'],systematic=systematic,hptslice=hptslice,mslice=mslice,mrebin=None)
        qcd_data = qcd_data - other_mc
        qcd_data = np.array([v if v>0. else default for v in qcd_data])
        qcd_data_int = hist.poisson_interval(qcd_data,qcd_data)
        qcd_temp = qcd_data
        qcd_temp_dn = np.array([qcd_data_int[1][bi] if qcd_data_int[1][bi] > 0. and not np.isnan(qcd_data_int[1][bi]) else default for bi in range(len(qcd_data))])
        qcd_temp_up = np.array([qcd_data_int[0][bi] if qcd_data_int[0][bi] > 0. and not np.isnan(qcd_data_int[0][bi]) else default for bi in range(len(qcd_data))])
        return qcd_temp,qcd_temp_dn,qcd_temp_up


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


    qcd_data_sel = {}
    top_data_sel = {}
    wlnu_data_sel = {}
    sig_data_sel = {}

    qcd_data_invdphi = {}
    top_data_invdphi = {}
    wlnu_data_invdphi = {}
    sig_data_invdphi = {}

    for reg in ["fail","loosepass","pass"]:
        qcd_data_sel[reg]  = [getQCDfromData(qcd_cr_hist, reg, hptslice=ptslices[ipt])[0] for ipt in range(npt)]
        top_data_sel[reg]  = [getQCDfromData(top_cr_hist, reg, hptslice=ptslices[ipt])[0] for ipt in range(npt)]
        wlnu_data_sel[reg] = [getQCDfromData(wlnu_cr_hist, reg, hptslice=ptslices[ipt])[0] for ipt in range(npt)]
        sig_data_sel[reg]  = [getQCDfromData(sig_met_hist, reg, hptslice=ptslices[ipt])[0] for ipt in range(npt)]

        qcd_data_invdphi[reg]  = [getQCDfromData(qcd_cr_invdphi, reg, hptslice=ptslices[ipt])[0] for ipt in range(npt)]
        top_data_invdphi[reg]  = [getQCDfromData(top_cr_invdphi, reg, hptslice=ptslices[ipt])[0] for ipt in range(npt)]
        wlnu_data_invdphi[reg] = [getQCDfromData(wlnu_cr_invdphi, reg, hptslice=ptslices[ipt])[0] for ipt in range(npt)]
        sig_data_invdphi[reg]  = [getQCDfromData(sig_met_invdphi, reg, hptslice=ptslices[ipt])[0] for ipt in range(npt)]

    qcd_data_full_sel = {}
    top_data_full_sel = {}
    wlnu_data_full_sel = {}
    sig_data_full_sel = {}

    qcd_data_full_invdphi = {}
    top_data_full_invdphi = {}
    wlnu_data_full_invdphi = {}
    sig_data_full_invdphi = {}

    for reg in ["fail","loosepass","pass"]:
        qcd_data_full_sel[reg]  = [getQCDfromDataFullMass(qcd_cr_hist, reg, hptslice=ptslices[ipt])[0] for ipt in range(npt)]
        top_data_full_sel[reg]  = [getQCDfromDataFullMass(top_cr_hist, reg, hptslice=ptslices[ipt])[0] for ipt in range(npt)]
        wlnu_data_full_sel[reg] = [getQCDfromDataFullMass(wlnu_cr_hist, reg, hptslice=ptslices[ipt])[0] for ipt in range(npt)]
        sig_data_full_sel[reg]  = [getQCDfromDataFullMass(sig_met_hist, reg, hptslice=ptslices[ipt])[0] for ipt in range(npt)]

        qcd_data_full_invdphi[reg]  = [getQCDfromDataFullMass(qcd_cr_invdphi, reg, hptslice=ptslices[ipt])[0] for ipt in range(npt)]
        top_data_full_invdphi[reg]  = [getQCDfromDataFullMass(top_cr_invdphi, reg, hptslice=ptslices[ipt])[0] for ipt in range(npt)]
        wlnu_data_full_invdphi[reg] = [getQCDfromDataFullMass(wlnu_cr_invdphi, reg, hptslice=ptslices[ipt])[0] for ipt in range(npt)]
        sig_data_full_invdphi[reg]  = [getQCDfromDataFullMass(sig_met_invdphi, reg, hptslice=ptslices[ipt])[0] for ipt in range(npt)]

    #  A    = D*Gi/Ji
    #  Aalt = G*Di/Ji
    #  D = sig_inv_F
    #  GV = qcd_sel_F

    print('had')

    qcd_from_data_sig = {}
    qcd_from_data_sig_err = {}
    qcd_from_data_top = {}
    qcd_from_data_top_err = {}
    qcd_from_data_wlnu = {}
    qcd_from_data_wlnu_err = {}

    for reg in ["fail","loosepass","pass"]:
        qcd_sig = [((sig_data_invdphi["fail"][ipt]/np.sum(sig_data_invdphi["fail"][ipt]))+(qcd_data_sel["fail"][ipt]/np.sum(qcd_data_sel["fail"][ipt])))*np.sum(sig_data_invdphi[reg][ipt])*np.sum(qcd_data_sel[reg][ipt])/2. for ipt in range(npt)]
        qcd_from_data_sig_err[reg] = [np.abs((sig_data_invdphi["fail"][ipt]*np.sum(qcd_data_sel[reg][ipt])*np.sum(sig_data_invdphi[reg][ipt])/np.sum(sig_data_invdphi["fail"][ipt]))-qcd_sig[ipt])/(2.*np.sum(qcd_data_invdphi["fail"][ipt])) for ipt in range(npt)]
        qcd_from_data_sig[reg] = [qcd_sig[ipt]/np.sum(qcd_data_invdphi["fail"][ipt]) for ipt in range(npt)]
        qcd_top = [((top_data_invdphi["fail"][ipt]/np.sum(top_data_invdphi["fail"][ipt]))+(qcd_data_sel["fail"][ipt]/np.sum(qcd_data_sel["fail"][ipt])))*np.sum(top_data_invdphi[reg][ipt])*np.sum(qcd_data_sel[reg][ipt])/2. for ipt in range(npt)]
        qcd_from_data_top_err[reg] = [np.abs((top_data_invdphi["fail"][ipt]*np.sum(qcd_data_sel[reg][ipt])*np.sum(top_data_invdphi[reg][ipt])/np.sum(top_data_invdphi["fail"][ipt]))-qcd_top[ipt])/(2.*np.sum(qcd_data_invdphi["fail"][ipt])) for ipt in range(npt)]
        qcd_from_data_top[reg] = [qcd_top[ipt]/np.sum(qcd_data_invdphi["fail"][ipt]) for ipt in range(npt)]
        qcd_wlnu = [((wlnu_data_invdphi["fail"][ipt]/np.sum(wlnu_data_invdphi["fail"][ipt]))+(qcd_data_sel["fail"][ipt]/np.sum(qcd_data_sel["fail"][ipt])))*np.sum(wlnu_data_invdphi[reg][ipt])*np.sum(qcd_data_sel[reg][ipt])/2. for ipt in range(npt)]
        qcd_from_data_wlnu_err[reg] = [np.abs((wlnu_data_invdphi["fail"][ipt]*np.sum(qcd_data_sel[reg][ipt])*np.sum(wlnu_data_invdphi[reg][ipt])/np.sum(wlnu_data_invdphi["fail"][ipt]))-qcd_wlnu[ipt])/(2.*np.sum(qcd_data_invdphi["fail"][ipt])) for ipt in range(npt)]
        qcd_from_data_wlnu[reg] = [qcd_wlnu[ipt]/np.sum(qcd_data_invdphi["fail"][ipt]) for ipt in range(npt)]

    #print('nom',qcd_from_data_sig)
    #print('err',qcd_from_data_sig_err)

    qcd_from_data_full_sig = {}
    qcd_from_data_full_sig_err = {}
    qcd_from_data_full_top = {}
    qcd_from_data_full_top_err = {}
    qcd_from_data_full_wlnu = {}
    qcd_from_data_full_wlnu_err = {}

    for reg in ["fail","loosepass","pass"]:
        qcd_sig = [((sig_data_full_invdphi["fail"][ipt]/np.sum(sig_data_full_invdphi["fail"][ipt]))+(qcd_data_full_sel["fail"][ipt]/np.sum(qcd_data_full_sel["fail"][ipt])))*np.sum(sig_data_full_invdphi[reg][ipt])*np.sum(qcd_data_full_sel[reg][ipt])/2. for ipt in range(npt)]
        qcd_from_data_full_sig_err[reg] = [np.abs((sig_data_full_invdphi["fail"][ipt]*np.sum(qcd_data_full_sel[reg][ipt])*np.sum(sig_data_full_invdphi[reg][ipt])/np.sum(sig_data_full_invdphi["fail"][ipt]))-qcd_sig[ipt])/(2.*np.sum(qcd_data_full_invdphi["fail"][ipt])) for ipt in range(npt)]
        qcd_from_data_full_sig[reg] = [qcd_sig[ipt]/np.sum(qcd_data_full_invdphi["fail"][ipt]) for ipt in range(npt)]
        qcd_top = [((top_data_full_invdphi["fail"][ipt]/np.sum(top_data_full_invdphi["fail"][ipt]))+(qcd_data_full_sel["fail"][ipt]/np.sum(qcd_data_full_sel["fail"][ipt])))*np.sum(top_data_full_invdphi[reg][ipt])*np.sum(qcd_data_full_sel[reg][ipt])/2. for ipt in range(npt)]
        qcd_from_data_full_top_err[reg] = [np.abs((top_data_full_invdphi["fail"][ipt]*np.sum(qcd_data_full_sel[reg][ipt])*np.sum(top_data_full_invdphi[reg][ipt])/np.sum(top_data_full_invdphi["fail"][ipt]))-qcd_top[ipt])/(2.*np.sum(qcd_data_full_invdphi["fail"][ipt])) for ipt in range(npt)]
        qcd_from_data_full_top[reg] = [qcd_top[ipt]/np.sum(qcd_data_full_invdphi["fail"][ipt]) for ipt in range(npt)]
        qcd_wlnu = [((wlnu_data_full_invdphi["fail"][ipt]/np.sum(wlnu_data_full_invdphi["fail"][ipt]))+(qcd_data_full_sel["fail"][ipt]/np.sum(qcd_data_full_sel["fail"][ipt])))*np.sum(wlnu_data_full_invdphi[reg][ipt])*np.sum(qcd_data_full_sel[reg][ipt])/2. for ipt in range(npt)]
        qcd_from_data_full_wlnu_err[reg] = [np.abs((wlnu_data_full_invdphi["fail"][ipt]*np.sum(qcd_data_full_sel[reg][ipt])*np.sum(wlnu_data_full_invdphi[reg][ipt])/np.sum(wlnu_data_full_invdphi["fail"][ipt]))-qcd_wlnu[ipt])/(2.*np.sum(qcd_data_full_invdphi["fail"][ipt])) for ipt in range(npt)]
        qcd_from_data_full_wlnu[reg] = [qcd_wlnu[ipt]/np.sum(qcd_data_full_invdphi["fail"][ipt]) for ipt in range(npt)]

    #print('nom_full',qcd_from_data_full_sig)
    #print('err_full',qcd_from_data_full_sig_err)

    for ptbin in range(npt):
        for region in ['pass','loosepass','fail']:
            ch = rl.Channel("ptbin%d%s" % (ptbin, region))
            model.addChannel(ch)

            isPass = region=='pass'
            isLoosePass = region=='loosepass'

            thehist = sig_met_hist
    
            for sName in thehist.identifiers('sample'):
                if sName.name=='ignore':
                    continue
                if sName.name=='data_obs' and isPass and not usingData:
                    continue
                templ = (intRegion(thehist[sName],region,hptslice=ptslices[ptbin])[1:-1], mtt.binning, mtt.name)
                #print('\tDEBUG',region,ptbin,sName.name,templ[0])
                for iv,val in enumerate(templ[0]):
                    if val<0.:
                        templ[0][iv] = 0.
                        #print('\t\t',templ[0])
                if sName.name=='data_obs':
                    ch.setObservation(templ)
                else:
                    if sName.name=='multijet':
                        qcdpred = qcd_from_data_sig[region][ptbin]
                        qcdpred_dn = qcd_from_data_sig[region][ptbin] - qcd_from_data_sig_err[region][ptbin]
                        qcdpred_up = qcd_from_data_sig[region][ptbin] + qcd_from_data_sig_err[region][ptbin]
                        #print('qcdpred',qcdpred)
                        #print('qcdpred_dn',qcdpred_dn)
                        #print('qcdpred_up',qcdpred_up)
                        templ = (qcdpred, mtt.binning, mtt.name)
                    stype = rl.Sample.SIGNAL if sName.name == 'htt125' else rl.Sample.BACKGROUND
                    sample = rl.TemplateSample(ch.name + '_' + sName.name, stype, templ)
        
                    #sample.setParamEffect(jec, 1.05)
                    if sName.name in syst_dict:
                        nom = templ[0]
                        for syst in syst_dict[sName.name]:
                            syst_params = syst_dict[sName.name][syst]
                            dnmod = intRegion(thehist[sName],region,syst_params[1],hptslice=ptslices[ptbin])[1:-1]
                            upmod = intRegion(thehist[sName],region,syst_params[2],hptslice=ptslices[ptbin])[1:-1]
                            sample.setParamEffect(syst_params[0], np.divide(dnmod, nom, out=np.ones_like(nom), where=nom>0.) if (dnmod!=nom).all() else np.ones_like(nom)*0.999, np.divide(upmod, nom, out=np.ones_like(nom), where=nom>0.) if (upmod!=nom).all() else np.ones_like(nom)*1.001)

                    #if sName.name=='top':
                    #    sample.setParamEffect(top_norm, 1.05)
                    #if sName.name=='wlnu':
                    #    sample.setParamEffect(wlnu_norm, 1.10)
                    if sName.name=='multijet':
                        sample.setParamEffect(qcd_norm, 1.20)
                        sample.setParamEffect(qcd_pass if isPass else qcd_loosepass if isLoosePass else qcd_fail, np.divide(qcdpred_dn, qcdpred, out=np.ones_like(qcdpred), where=qcdpred>0.), np.divide(qcdpred_up, qcdpred, out=np.ones_like(qcdpred), where=qcdpred>0.))
                    if sName.name=='vqq':
                        sample.setParamEffect(vqq_norm, 1.10)
                    #if sName.name=='ztt':
                    #    sample.setParamEffect(ztt_norm, 1.10)
                    if sName.name=='htt125' or sName.name=='ztt':
                        nom_full = intRegion(thehist[sName],region)
                        shift_dn = nom_full[2:]
                        shift_up = nom_full[:-2]
                        #sample.setParamEffect(m_scale, np.divide(shift_dn, nom_full[1:-1], out=np.ones_like(nom_full[1:-1]), where=nom_full[1:-1]>0.), np.divide(shift_up, nom_full[1:-1], out=np.ones_like(nom_full[1:-1]), where=nom_full[1:-1]>0.))
                        #print('had htt125 ptbin',ptbin,region,nom_full[1:-1])
                    sample.setParamEffect(trig, 1.02)
                    sample.setParamEffect(lumi, 1.024)
            
                    ch.addSample(sample)
            if not usingData and isPass:
                ch.setObservation((np.zeros(len(mtt.binning)-1),mtt.binning, mtt.name))
    
    h125vals = [intRegion(thehist['htt125'],'pass',hptslice=ptslices[ptbin])[1:-1] for ptbin in range(npt)]
    zttvals = [intRegion(thehist['ztt'],'pass',hptslice=ptslices[ptbin])[1:-1] for ptbin in range(npt)]
    wlnuvals = [intRegion(thehist['wlnu'],'pass',hptslice=ptslices[ptbin])[1:-1] for ptbin in range(npt)]
    probs = [np.prod([math.erfc(h125vals[ipt][ib]/(math.sqrt(abs(zttvals[ipt][ib]+wlnuvals[ipt][ib])) if zttvals[ipt][ib]+wlnuvals[ipt][ib]>2. else 2.)) for ib in range(len(h125vals[ipt]))]) for ipt in range(npt)]
    print('had','PROB',np.prod(probs),'->',erfinv(1.-np.prod(probs)),'\t\t',probs)
    

    for ptbin in range(npt):
        failCh = model['ptbin%dfail' % ptbin]
        passCh = model['ptbin%dpass' % ptbin]
        loosePassCh = model['ptbin%dloosepass' % ptbin]

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
        toppass.setParamEffect(topnormSF, 1*topnormSF)
        toploosepass.setParamEffect(topnormSF, 1*topnormSF)
        topfail.setParamEffect(topnormSF, 1*topnormSF)

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
        wlnupass.setParamEffect(wlnunormSF, 1*wlnunormSF)
        wlnuloosepass.setParamEffect(wlnunormSF, 1*wlnunormSF)
        wlnufail.setParamEffect(wlnunormSF, 1*wlnunormSF)

        zttpass = passCh['ztt']
        zttloosepass = loosePassCh['ztt']
        zttfail = failCh['ztt']
        zttLP = zttpass.getExpectation(nominal=True).sum() / zttloosepass.getExpectation(nominal=True).sum()
        httLP = passCh['htt125'].getExpectation(nominal=True).sum() / loosePassCh['htt125'].getExpectation(nominal=True).sum()
        zttpass.setParamEffect(ztt_eff, 1*ztt_eff)
        zttloosepass.setParamEffect(ztt_eff, (1 - ztt_eff) * zttLP + 1)
        passCh['htt125'].setParamEffect(ztt_eff, 1*ztt_eff)
        loosePassCh['htt125'].setParamEffect(ztt_eff, (1 - ztt_eff) * httLP + 1)
        zttpass.setParamEffect(rztt, 1*rztt)
        zttloosepass.setParamEffect(rztt, 1*rztt)
        zttfail.setParamEffect(rztt, 1*rztt)
        #zttpass.setParamEffect(rztt, 1.05)
        #zttloosepass.setParamEffect(rztt, 1.05)
        #zttfail.setParamEffect(rztt, 1.05)


    # Fill in top CR
    for region in ['pass', 'loosepass', 'fail']:
        ch = rl.Channel("topCR%s" % (region, ))
        model.addChannel(ch)

        isPass = region=='pass'
        isLoosePass = region=='loosepass'

        thehist = top_cr_hist

        for sName in thehist.identifiers('sample'):
            if sName.name=='ignore':
                continue
            if singleBinCR:
                templ = (np.array([np.sum(intRegion(thehist[sName],region)[1:-1])]), mttone.binning, mttone.name)
                if includeLowMass:
                    templ = (np.array([np.sum(intRegion(thehist[sName],region,mrebin=None))]), mttone.binning, mttone.name)
            else:
                templ = (intRegion(thehist[sName],region)[1:-1], mtt.binning, mtt.name)
            for iv,val in enumerate(templ[0]):
                if val<0.:
                    templ[0][iv] = 0.
            if sName.name=='data_obs':
                ch.setObservation(templ)
            else:
                if sName.name=='multijet':
                    qcdpred = np.sum(np.array(qcd_from_data_top[region]), axis=0)
                    qcdpred_dn = np.sum(np.array(qcd_from_data_top[region]) - np.array(qcd_from_data_top_err[region]), axis=0)
                    qcdpred_up = np.sum(np.array(qcd_from_data_top[region]) + np.array(qcd_from_data_top_err[region]), axis=0)
                    if singleBinCR:
                        if includeLowMass:
                            qcdpred = np.sum(np.array(qcd_from_data_full_top[region]), axis=0)
                            qcdpred_dn = np.sum(np.array(qcd_from_data_full_top[region]) - np.array(qcd_from_data_full_top_err[region]), axis=0)
                            qcdpred_up = np.sum(np.array(qcd_from_data_full_top[region]) + np.array(qcd_from_data_full_top_err[region]), axis=0)
                        templ = (np.array([np.sum(qcdpred)]), mttone.binning, mttone.name)
                    else:
                        templ = (qcdpred, mtt.binning, mtt.name)
                stype = rl.Sample.SIGNAL if sName.name == 'htt125' else rl.Sample.BACKGROUND
                sample = rl.TemplateSample(ch.name + '_' + sName.name, stype, templ)
    
                #sample.setParamEffect(jec, 1.05)
                if sName.name in syst_dict:
                    nom_base = intRegion(thehist[sName],region)[1:-1]
                    if singleBinCR:
                        if includeLowMass:
                            nom_base = intRegion(thehist[sName],region,mrebin=None)
                    for syst in syst_dict[sName.name]:
                        syst_params = syst_dict[sName.name][syst]
                        dnmod = intRegion(thehist[sName],region,systematic=syst_params[1])[1:-1]
                        upmod = intRegion(thehist[sName],region,systematic=syst_params[2])[1:-1]
                        if singleBinCR:
                            if includeLowMass:
                                dnmod = intRegion(thehist[sName],region,mrebin=None,systematic=syst_params[1])
                                upmod = intRegion(thehist[sName],region,mrebin=None,systematic=syst_params[2])
                            dnmod = np.array([np.sum(dnmod)])
                            nom = np.array([np.sum(nom_base)])
                            upmod = np.array([np.sum(upmod)])
                        else:
                            nom = nom_base
                        sample.setParamEffect(syst_params[0], np.divide(dnmod, nom, out=np.ones_like(nom), where=nom>0.) if (dnmod!=nom).all() else np.ones_like(nom)*0.999, np.divide(upmod, nom, out=np.ones_like(nom), where=nom>0.) if (upmod!=nom).all() else np.ones_like(nom)*1.001)

                #if sName.name=='top':
                #    sample.setParamEffect(top_norm, 1.05)
                #if sName.name=='wlnu':
                #    sample.setParamEffect(wlnu_norm, 1.10)
                if sName.name=='multijet':
                    sample.setParamEffect(qcd_norm, 1.20)
                    if singleBinCR:
                        qcdpred_dn = np.array([np.sum(qcdpred_dn)])
                        qcdpred = np.array([np.sum(qcdpred)])
                        qcdpred_up = np.array([np.sum(qcdpred_up)])
                    sample.setParamEffect(qcdtop_pass if isPass else qcdtop_loosepass if isLoosePass else qcdtop_fail, np.divide(qcdpred_dn, qcdpred, out=np.ones_like(qcdpred), where=qcdpred>0.), np.divide(qcdpred_up, qcdpred, out=np.ones_like(qcdpred), where=qcdpred>0.))
                if sName.name=='vqq':
                    sample.setParamEffect(vqq_norm, 1.10)
                #if sName.name=='ztt':
                #    sample.setParamEffect(ztt_norm, 1.10)
                if sName.name=='htt125' or sName.name=='ztt':
                    nom_full = intRegion(thehist[sName],region)
                    shift_dn = nom_full[2:]
                    shift_up = nom_full[:-2]
                    #sample.setParamEffect(m_scale, np.divide(shift_dn, nom_full[1:-1], out=np.ones_like(nom_full[1:-1]), where=nom_full[1:-1]>0.), np.divide(shift_up, nom_full[1:-1], out=np.ones_like(nom_full[1:-1]), where=nom_full[1:-1]>0.))
                sample.setParamEffect(trig, 1.02)
                sample.setParamEffect(lumi, 1.024)
        
                ch.addSample(sample)

    qcdpass = model['topCRpass']['multijet']
    qcdloosepass = model['topCRloosepass']['multijet']
    qcdfail = model['topCRfail']['multijet']
    qcdpass.setParamEffect(qcdnormSF_top, 1.20)
    qcdloosepass.setParamEffect(qcdnormSF_top, 1.20)
    qcdfail.setParamEffect(qcdnormSF_top, 1.20)
    #qcdpass.setParamEffect(qcdnormSF_top, 1*qcdnormSF_top)
    #qcdloosepass.setParamEffect(qcdnormSF_top, 1*qcdnormSF_top)
    #qcdfail.setParamEffect(qcdnormSF_top, 1*qcdnormSF_top)

    toppass = model['topCRpass']['top']
    toploosepass = model['topCRloosepass']['top']
    topfail = model['topCRfail']['top']
    topPF = (toploosepass.getExpectation(nominal=True).sum() + toppass.getExpectation(nominal=True).sum()) / topfail.getExpectation(nominal=True).sum()
    topLPF = toppass.getExpectation(nominal=True).sum() / toploosepass.getExpectation(nominal=True).sum()
    topRLPF = 1./(1.+(1./topLPF)) # P/(L+P)
    toppass.setParamEffect(topLeffSF, 1*topLeffSF)
    toploosepass.setParamEffect(topLeffSF, (1 - topLeffSF) * topLPF + 1)
    toploosepass.setParamEffect(topeffSF, 1*topeffSF)
    topfail.setParamEffect(topeffSF, (1 - topeffSF) * (1-topRLPF) * topPF + 1)
    toppass.setParamEffect(topnormSF, 1*topnormSF)
    toploosepass.setParamEffect(topnormSF, 1*topnormSF)
    topfail.setParamEffect(topnormSF, 1*topnormSF)

    wlnupass = model['topCRpass']['wlnu']
    wlnuloosepass = model['topCRloosepass']['wlnu']
    wlnufail = model['topCRfail']['wlnu']
    wlnuPF = (wlnuloosepass.getExpectation(nominal=True).sum() + wlnupass.getExpectation(nominal=True).sum()) / wlnufail.getExpectation(nominal=True).sum()
    wlnuLPF = wlnupass.getExpectation(nominal=True).sum() / wlnuloosepass.getExpectation(nominal=True).sum()
    wlnuRLPF = 1./(1.+(1./wlnuLPF)) # P/(L+P)
    wlnupass.setParamEffect(wlnuLeffSF, 1*wlnuLeffSF)
    wlnuloosepass.setParamEffect(wlnuLeffSF, (1 - wlnuLeffSF) * wlnuLPF + 1)
    wlnuloosepass.setParamEffect(wlnueffSF, 1*wlnueffSF)
    wlnufail.setParamEffect(wlnueffSF, (1 - wlnueffSF) * (1-wlnuRLPF) * wlnuPF + 1)
    wlnupass.setParamEffect(wlnunormSF, 1*wlnunormSF)
    wlnuloosepass.setParamEffect(wlnunormSF, 1*wlnunormSF)
    wlnufail.setParamEffect(wlnunormSF, 1*wlnunormSF)

    zttpass = model['topCRpass']['ztt']
    zttloosepass = model['topCRloosepass']['ztt']
    zttfail = model['topCRfail']['ztt']
    zttLP = zttpass.getExpectation(nominal=True).sum() / zttloosepass.getExpectation(nominal=True).sum()
    httLP = model['topCRpass']['htt125'].getExpectation(nominal=True).sum() / model['topCRloosepass']['htt125'].getExpectation(nominal=True).sum()
    zttpass.setParamEffect(ztt_eff, 1*ztt_eff)
    zttloosepass.setParamEffect(ztt_eff, (1 - ztt_eff)* zttLP + 1)
    model['topCRpass']['htt125'].setParamEffect(ztt_eff, 1*ztt_eff)
    model['topCRloosepass']['htt125'].setParamEffect(ztt_eff, (1 - ztt_eff) * httLP + 1)
    zttpass.setParamEffect(rztt, 1*rztt)
    zttloosepass.setParamEffect(rztt, 1*rztt)
    zttfail.setParamEffect(rztt, 1*rztt)
    #zttpass.setParamEffect(rztt, 1.05)
    #zttloosepass.setParamEffect(rztt, 1.05)
    #zttfail.setParamEffect(rztt, 1.05)

    # Fill in wlnu CR
    for region in ['pass', 'fail', 'loosepass']:
        ch = rl.Channel("wlnuCR%s" % (region, ))
        model.addChannel(ch)

        isPass = region=='pass'
        isLoosePass = region=='loosepass'

        thehist = wlnu_cr_hist

        for sName in thehist.identifiers('sample'):
            if sName.name=='ignore':
                continue
            if singleBinCR:
                templ = (np.array([np.sum(intRegion(thehist[sName],region)[1:-1])]), mttone.binning, mttone.name)
                if includeLowMass:
                    templ = (np.array([np.sum(intRegion(thehist[sName],region,mrebin=None))]), mttone.binning, mttone.name)
            else:
                templ = (intRegion(thehist[sName],region)[1:-1], mtt.binning, mtt.name)
            for iv,val in enumerate(templ[0]):
                if val<0.:
                    templ[0][iv] = 0.
            if sName.name=='data_obs':
                ch.setObservation(templ)
            else:
                if sName.name=='multijet':
                    qcdpred = np.sum(np.array(qcd_from_data_wlnu[region]), axis=0)
                    qcdpred_dn = np.sum(np.array(qcd_from_data_wlnu[region]) - np.array(qcd_from_data_wlnu_err[region]), axis=0)
                    qcdpred_up = np.sum(np.array(qcd_from_data_wlnu[region]) + np.array(qcd_from_data_wlnu_err[region]), axis=0)
                    if singleBinCR:
                        if includeLowMass:
                            qcdpred = np.sum(np.array(qcd_from_data_full_wlnu[region]), axis=0)
                            qcdpred_dn = np.sum(np.array(qcd_from_data_full_wlnu[region]) - np.array(qcd_from_data_full_wlnu_err[region]), axis=0)
                            qcdpred_up = np.sum(np.array(qcd_from_data_full_wlnu[region]) + np.array(qcd_from_data_full_wlnu_err[region]), axis=0)
                        templ = (np.array([np.sum(qcdpred)]), mttone.binning, mttone.name)
                    else:
                        templ = (qcdpred, mtt.binning, mtt.name)
                stype = rl.Sample.SIGNAL if sName.name == 'htt125' else rl.Sample.BACKGROUND
                sample = rl.TemplateSample(ch.name + '_' + sName.name, stype, templ)
    
                #sample.setParamEffect(jec, 1.05)
                if sName.name in syst_dict:
                    nom_base = intRegion(thehist[sName],region)[1:-1]
                    if singleBinCR:
                        if includeLowMass:
                            nom_base = intRegion(thehist[sName],region,mrebin=None)
                    for syst in syst_dict[sName.name]:
                        syst_params = syst_dict[sName.name][syst]
                        dnmod = intRegion(thehist[sName],region,systematic=syst_params[1])[1:-1]
                        upmod = intRegion(thehist[sName],region,systematic=syst_params[2])[1:-1]
                        if singleBinCR:
                            if includeLowMass:
                                dnmod = intRegion(thehist[sName],region,mrebin=None,systematic=syst_params[1])
                                upmod = intRegion(thehist[sName],region,mrebin=None,systematic=syst_params[2])
                            dnmod = np.array([np.sum(dnmod)])
                            nom = np.array([np.sum(nom_base)])
                            upmod = np.array([np.sum(upmod)])
                        else:
                            nom = nom_base
                        sample.setParamEffect(syst_params[0], np.divide(dnmod, nom, out=np.ones_like(nom), where=nom>0.) if (dnmod!=nom).all() else np.ones_like(nom)*0.999, np.divide(upmod, nom, out=np.ones_like(nom), where=nom>0.) if (upmod!=nom).all() else np.ones_like(nom)*1.001)

                #if sName.name=='top':
                #    sample.setParamEffect(top_norm, 1.05)
                #if sName.name=='wlnu':
                #    sample.setParamEffect(wlnu_norm, 1.10)
                if sName.name=='multijet':
                    sample.setParamEffect(qcd_norm, 1.20)
                    if singleBinCR:
                        qcdpred_dn = np.array([np.sum(qcdpred_dn)])
                        qcdpred = np.array([np.sum(qcdpred)])
                        qcdpred_up = np.array([np.sum(qcdpred_up)])
                    sample.setParamEffect(qcdwlnu_pass if isPass else qcdwlnu_loosepass if isLoosePass else qcdwlnu_fail, np.divide(qcdpred_dn, qcdpred, out=np.ones_like(qcdpred), where=qcdpred>0.), np.divide(qcdpred_up, qcdpred, out=np.ones_like(qcdpred), where=qcdpred>0.))
                if sName.name=='vqq':
                    sample.setParamEffect(vqq_norm, 1.10)
                #if sName.name=='ztt':
                #    sample.setParamEffect(ztt_norm, 1.10)
                if sName.name=='htt125' or sName.name=='ztt':
                    nom_full = intRegion(thehist[sName],region)
                    shift_dn = nom_full[2:]
                    shift_up = nom_full[:-2]
                    #sample.setParamEffect(m_scale, np.divide(shift_dn, nom_full[1:-1], out=np.ones_like(nom_full[1:-1]), where=nom_full[1:-1]>0.), np.divide(shift_up, nom_full[1:-1], out=np.ones_like(nom_full[1:-1]), where=nom_full[1:-1]>0.))
                sample.setParamEffect(trig, 1.02)
                sample.setParamEffect(lumi, 1.024)
        
                ch.addSample(sample)

    qcdpass = model['wlnuCRpass']['multijet']
    qcdloosepass = model['wlnuCRloosepass']['multijet']
    qcdfail = model['wlnuCRfail']['multijet']
    qcdpass.setParamEffect(qcdnormSF_wlnu, 1.20)
    qcdloosepass.setParamEffect(qcdnormSF_wlnu, 1.20)
    qcdfail.setParamEffect(qcdnormSF_wlnu, 1.20)
    #qcdpass.setParamEffect(qcdnormSF_wlnu, 1*qcdnormSF_wlnu)
    #qcdloosepass.setParamEffect(qcdnormSF_wlnu, 1*qcdnormSF_wlnu)
    #qcdfail.setParamEffect(qcdnormSF_wlnu, 1*qcdnormSF_wlnu)

    toppass = model['wlnuCRpass']['top']
    toploosepass = model['wlnuCRloosepass']['top']
    topfail = model['wlnuCRfail']['top']
    topPF = (toploosepass.getExpectation(nominal=True).sum() + toppass.getExpectation(nominal=True).sum()) / topfail.getExpectation(nominal=True).sum()
    topLPF = toppass.getExpectation(nominal=True).sum() / toploosepass.getExpectation(nominal=True).sum()
    topRLPF = 1./(1.+(1./topLPF)) # P/(L+P)
    toppass.setParamEffect(topLeffSF, 1*topLeffSF)
    toploosepass.setParamEffect(topLeffSF, (1 - topLeffSF) * topLPF + 1)
    toploosepass.setParamEffect(topeffSF, 1*topeffSF)
    topfail.setParamEffect(topeffSF, (1 - topeffSF) * (1-topRLPF) * topPF + 1)
    toppass.setParamEffect(topnormSF, 1*topnormSF)
    toploosepass.setParamEffect(topnormSF, 1*topnormSF)
    topfail.setParamEffect(topnormSF, 1*topnormSF)

    wlnupass = model['wlnuCRpass']['wlnu']
    wlnuloosepass = model['wlnuCRloosepass']['wlnu']
    wlnufail = model['wlnuCRfail']['wlnu']
    wlnuPF = (wlnuloosepass.getExpectation(nominal=True).sum() + wlnupass.getExpectation(nominal=True).sum()) / wlnufail.getExpectation(nominal=True).sum()
    wlnuLPF = wlnupass.getExpectation(nominal=True).sum() / wlnuloosepass.getExpectation(nominal=True).sum()
    wlnuRLPF = 1./(1.+(1./wlnuLPF)) # P/(L+P)
    wlnupass.setParamEffect(wlnuLeffSF, 1*wlnuLeffSF)
    wlnuloosepass.setParamEffect(wlnuLeffSF, (1 - wlnuLeffSF) * wlnuLPF + 1)
    wlnuloosepass.setParamEffect(wlnueffSF, 1*wlnueffSF)
    wlnufail.setParamEffect(wlnueffSF, (1 - wlnueffSF) * (1-wlnuRLPF) * wlnuPF + 1)
    wlnupass.setParamEffect(wlnunormSF, 1*wlnunormSF)
    wlnuloosepass.setParamEffect(wlnunormSF, 1*wlnunormSF)
    wlnufail.setParamEffect(wlnunormSF, 1*wlnunormSF)

    zttpass = model['wlnuCRpass']['ztt']
    zttloosepass = model['wlnuCRloosepass']['ztt']
    zttfail = model['wlnuCRfail']['ztt']
    zttLP = zttpass.getExpectation(nominal=True).sum() / zttloosepass.getExpectation(nominal=True).sum()
    httLP = model['wlnuCRpass']['htt125'].getExpectation(nominal=True).sum() / model['wlnuCRloosepass']['htt125'].getExpectation(nominal=True).sum()
    zttpass.setParamEffect(ztt_eff, 1*ztt_eff)
    zttloosepass.setParamEffect(ztt_eff, (1 - ztt_eff) * zttLP + 1)
    model['wlnuCRpass']['htt125'].setParamEffect(ztt_eff, 1*ztt_eff)
    model['wlnuCRloosepass']['htt125'].setParamEffect(ztt_eff, (1 - ztt_eff) * httLP + 1)
    zttpass.setParamEffect(rztt, 1*rztt)
    zttloosepass.setParamEffect(rztt, 1*rztt)
    zttfail.setParamEffect(rztt, 1*rztt)
    #zttpass.setParamEffect(rztt, 1.05)
    #zttloosepass.setParamEffect(rztt, 1.05)
    #zttfail.setParamEffect(rztt, 1.05)

    with open(os.path.join("%s/%s"%(str(tmpdir),label), 'hadhadModel.pkl'), "wb") as fout:
        pickle.dump(model, fout, protocol=2)

    model.renderCombine(os.path.join("%s/%s"%(str(tmpdir),label), 'hadhadModel'))

def getHist(h,var_name,lumifb,vars_cut,regionsel,blind,sigscale,rebin,debug=False):
    if debug: print(h)
    exceptions = ['process', var_name]
    for var in vars_cut:
        exceptions.append(var)
    if (regionsel is not ''):
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

    hists_mapped = {}
    for h in mcSamples + dataSamples:
        #print(h)
        # open hists
        hists_unmapped = load('%s%s.coffea'%(args.hist,h))
        # map to hists
        for key, val in hists_unmapped.items():
            if isinstance(val, hist.Hist):
                if key in hists_mapped:
                    hists_mapped[key] = hists_mapped[key].add(processmap.apply(val))
                else:
                    hists_mapped[key] = processmap.apply(val)

    nnCut_hadhad = args.nnCutHH
    nnCut_hadhad_loose = args.nnCutHHL
    nnCut_hadhad_met = args.nnCutHHMET
    nnCut_hadhad_met_loose = args.nnCutHHMETL
    nnCut_hadel = args.nnCutEH
    nnCut_hadel_loose = args.nnCutEHL
    nnCut_hadmu = args.nnCutMH
    nnCut_hadmu_loose = args.nnCutMHL
    metCut = args.metCut
    h_pt_min = args.hPtCut
    antimucut = args.antiMuCut
    antielcut = args.antiElCut
    
    includeData = False
    mttbins = np.linspace(40, 210, 18)
    # properties
    region_dict = {
#        'hadhad_signal':['met_nn_kin', 'massreg', {'h_pt':[],'massreg':[mttbins[0], mttbins[-1]],'nn_disc':[],'jetmet_dphi':[], 'antilep':[0.5,None], 'met_pt':[100.,150.],'systematic':[]}, includeData, [1], 'hadhad_signal'],
#        'hadhad_signal_met':['met_nn_kin', 'massreg', {'h_pt':[],'massreg':[mttbins[0], mttbins[-1]],'nn_disc':[],'jetmet_dphi':[], 'antilep':[0.5,None], 'met_pt':[150.,None],'systematic':[]}, includeData, [1], 'hadhad_signal_met'],
#        'hadhad_cr_mu':['met_nn_kin', 'massreg', {'h_pt':[],'massreg':[mttbins[0], mttbins[-1]],'nn_disc':[],'jetmet_dphi':[],'met_pt':[20.,None],'systematic':[]}, includeData, [1], 'hadhad_cr_mu'],
#        'hadhad_cr_b_mu_iso':['met_nn_kin', 'massreg', {'h_pt':[],'massreg':[mttbins[0], mttbins[-1]],'nn_disc':[],'jetmet_dphi':[],'met_pt':[20.,None],'systematic':[]}, includeData, [1], 'hadhad_cr_b_mu_iso'],
#        'hadhad_cr_mu_iso':['met_nn_kin', 'massreg', {'h_pt':[],'massreg':[mttbins[0], mttbins[-1]],'nn_disc':[],'jetmet_dphi':[],'met_pt':[20.,None],'systematic':[]}, includeData, [1], 'hadhad_cr_mu_iso'],
#        'hadel_signal':['met_nn_kin', 'massreg', {'h_pt':[],'massreg':[mttbins[0], mttbins[-1]],'nn_disc':[],'antilep':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadel_signal'],
#        'hadmu_signal':['met_nn_kin', 'massreg', {'h_pt':[],'massreg':[mttbins[0], mttbins[-1]],'nn_disc':[],'antilep':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadmu_signal'],
#        'hadel_cr_b':['met_nn_kin', 'massreg', {'h_pt':[],'massreg':[mttbins[0], mttbins[-1]],'nn_disc':[],'antilep':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadel_cr_b'],
#        'hadmu_cr_b':['met_nn_kin', 'massreg', {'h_pt':[],'massreg':[mttbins[0], mttbins[-1]],'nn_disc':[],'antilep':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadmu_cr_b'],
#        'hadel_cr_w':['met_nn_kin', 'massreg', {'h_pt':[],'massreg':[mttbins[0], mttbins[-1]],'nn_disc':[],'antilep':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadel_cr_w'],
#        'hadmu_cr_w':['met_nn_kin', 'massreg', {'h_pt':[],'massreg':[mttbins[0], mttbins[-1]],'nn_disc':[],'antilep':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadmu_cr_w'],
#        'hadel_cr_qcd':['met_nn_kin', 'massreg', {'h_pt':[],'massreg':[mttbins[0], mttbins[-1]],'nn_disc':[],'antilep':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadel_cr_qcd'],
#        'hadmu_cr_qcd':['met_nn_kin', 'massreg', {'h_pt':[],'massreg':[mttbins[0], mttbins[-1]],'nn_disc':[],'antilep':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadmu_cr_qcd'],
        'hadhad_signal':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[],'jetmet_dphi':[], 'antilep':[0.5,None], 'met_pt':[100.,150.],'systematic':[]}, includeData, [1], 'hadhad_signal'],
        'hadhad_signal_met':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[],'jetmet_dphi':[], 'antilep':[0.5,None], 'met_pt':[150.,None],'systematic':[]}, includeData, [1], 'hadhad_signal_met'],
        'hadhad_cr_mu':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[],'jetmet_dphi':[], 'antilep':[0.5,None],'met_pt':[20.,None],'systematic':[]}, includeData, [1], 'hadhad_cr_mu'],
        'hadhad_cr_b_mu_iso':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[],'jetmet_dphi':[], 'antilep':[0.5,None],'met_pt':[20.,None],'systematic':[]}, includeData, [1], 'hadhad_cr_b_mu_iso'],
        'hadhad_cr_mu_iso':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[],'jetmet_dphi':[], 'antilep':[0.5,None],'met_pt':[20.,None],'systematic':[]}, includeData, [1], 'hadhad_cr_mu_iso'],
        'hadel_signal':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[],'antilep':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadel_signal'],
        'hadmu_signal':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[],'antilep':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadmu_signal'],
        'hadel_cr_b':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[],'antilep':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadel_cr_b'],
        'hadmu_cr_b':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[],'antilep':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadmu_cr_b'],
        'hadel_cr_w':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[],'antilep':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadel_cr_w'],
        'hadmu_cr_w':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[],'antilep':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadmu_cr_w'],
        'hadel_cr_qcd':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[],'antilep':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadel_cr_qcd'],
        'hadmu_cr_qcd':['met_nn_kin', 'massreg', {'h_pt':[],'nn_disc':[],'antilep':[], 'met_pt':[],'systematic':[]}, includeData, [1], 'hadmu_cr_qcd'],
    }

    full_dict = {}
    for reg in region_dict:
        reg_opt = region_dict[reg]
        h = hists_mapped[reg_opt[0]]
        #print(h)
        xhist = getHist(h,reg_opt[1],args.lumi,reg_opt[2],[reg_opt[5]],reg_opt[3],args.sigscale,reg_opt[4])
        #print(xhist)
        full_dict[reg] = xhist

    hadptbins = np.array([h_pt_min]+[float(b) if b!="None" else None for b in args.hPtBinsHad])
    lepptbins = np.array([h_pt_min]+[float(b) if b!="None" else None for b in args.hPtBinsLep])

    createHadHad(full_dict['hadhad_signal'], full_dict['hadhad_signal_met'], full_dict['hadhad_cr_b_mu_iso'], full_dict['hadhad_cr_mu_iso'], full_dict['hadhad_cr_mu'], 'massreg', mttbins, hadptbins, odir, args.label, args.unblind, nnCut_hadhad, nnCut_hadhad_loose, nnCut_hadhad_met, nnCut_hadhad_met_loose, metCut, h_pt_min)
    createLepHad(full_dict['hadel_signal'], full_dict['hadel_cr_b'], full_dict['hadel_cr_w'], full_dict['hadel_cr_qcd'], 'massreg', mttbins, lepptbins, "el", odir, args.label, args.unblind, nnCut_hadel, nnCut_hadel_loose, metCut, h_pt_min, antielcut)
    createLepHad(full_dict['hadmu_signal'], full_dict['hadmu_cr_b'], full_dict['hadmu_cr_w'], full_dict['hadmu_cr_qcd'], 'massreg', mttbins, lepptbins, "mu", odir, args.label, args.unblind, nnCut_hadmu, nnCut_hadmu_loose, metCut, h_pt_min, antimucut)

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
    f.write("\nMET cut: " + str(args.metCut))
    f.write("\nhPt cut: " + str(args.hPtCut))
    f.write("\nLepBins: " + ", ".join(args.hPtBinsLep))
    f.write("\nHadBins: " + ", ".join(args.hPtBinsHad))
    f.write("\nAntiEl cut: " + str(args.antiElCut))
    f.write("\nAntiMu cut: " + str(args.antiMuCut))
    f.close()

    if args.plots:
    
        pwd = os.getcwd()
        os.chdir('%s/%s'%(odir,args.label))
    
        for rn in ['hadhad_signal','hadhad_signal_met','hadhad_cr_b_mu_iso','hadhad_cr_mu_iso','hadhad_cr_mu']:
          if rn=='hadhad_signal':
            nncutl = nnCut_hadhad_loose
            nncutt = nnCut_hadhad
          elif rn=='hadhad_signal_met':
            nncutl = nnCut_hadhad_met_loose
            nncutt = nnCut_hadhad_met
          else:
            nncutl = nnCut_hadhad_met_loose
            nncutt = nnCut_hadhad_met
          plot_stack.drawStack(full_dict[rn],'HadhadModel','massreg',r"$m_{NN}$",'%s fail'%rn,'all',1.,{'nn_disc':[None,nncutl],'jetmet_dphi':[None,1.6]},{},[],rn+'_fail')
          plot_stack.drawStack(full_dict[rn],'HadhadModel','massreg',r"$m_{NN}$",'%s loosepass'%rn,'all',1.,{'nn_disc':[nncutl,nncutt],'jetmet_dphi':[None,1.6]},{},[],rn+'_loosepass')
          plot_stack.drawStack(full_dict[rn],'HadhadModel','massreg',r"$m_{NN}$",'%s pass'%rn,'all',1.,{'nn_disc':[nncutt,None],'jetmet_dphi':[None,1.6]},{},[],rn+'_pass', blind='' if args.unblind or 'signal' not in rn else ['None',110.])
        for lep in ['el','mu']:
          if lep=='el':
            nncutl = nnCut_hadel_loose
            nncutt = nnCut_hadel
          else:
            nncutl = nnCut_hadmu_loose
            nncutt = nnCut_hadmu
          for rn in ['had%s_signal'%lep,'had%s_cr_b'%lep,'had%s_cr_w'%lep,'had%s_cr_qcd'%lep]:
            plot_stack.drawStack(full_dict[rn],'Had%sModel'%lep,'massreg',r"$m_{NN}$",'%s fail'%rn,'all',1.,{'nn_disc':[None,nncutl],'met_pt':[metCut,None],'antilep':[0.5,None]},{},[],rn+'_fail')
            plot_stack.drawStack(full_dict[rn],'Had%sModel'%lep,'massreg',r"$m_{NN}$",'%s loosepass'%rn,'all',1.,{'nn_disc':[nncutl,nncutt],'met_pt':[metCut,None],'antilep':[0.5,None]},{},[],rn+'_loosepass')
            plot_stack.drawStack(full_dict[rn],'Had%sModel'%lep,'massreg',r"$m_{NN}$",'%s pass'%rn,'all',1.,{'nn_disc':[nncutt,None],'met_pt':[metCut,None],'antilep':[0.5,None]},{},[],rn+'_pass', blind='' if args.unblind or 'signal' not in rn else ['None',110.])
    
        os.chdir(pwd)

if __name__ == "__main__":

    if not sys.warnoptions:
        import warnings
        warnings.simplefilter("ignore")
        os.environ["PYTHONWARNINGS"] = "ignore" # Also affect subprocesses

    #ex. python makeCards.py --hist ../../boostedhiggs/condor/Mar08_NN/hists_sum_ --year 2017 --lumi 36.7 --tag Mar08 --label 01
    parser = argparse.ArgumentParser()
    parser.add_argument('--hist',        dest='hist',        default="hists_sum_",   help="hists pickle prefix")
    parser.add_argument('--year',        dest='year',        default="2017",         help="year")
    parser.add_argument('--lumi',        dest='lumi',        default=50.,            help="lumi",       type=float)
    parser.add_argument('--tag',         dest='tag',         default="",             help="tag")
    parser.add_argument('--label',       dest='label',       default="",             help="label")
    parser.add_argument('--sigscale',    dest='sigscale',    default=1.,             help="sigscale",   type=float)
    parser.add_argument('--unblind',     dest='unblind',     action="store_true",    help="unblind")
    parser.add_argument('--plots',       dest='plots',       action="store_true",    help="plots")
    parser.add_argument('--nnCutHH',     dest='nnCutHH',     default=0.9999,         help="nncut hadhad",   type=float)
    parser.add_argument('--nnCutHHMET',  dest='nnCutHHMET',  default=0.9999,         help="nncut hadhad met",   type=float)
    parser.add_argument('--nnCutHHL',    dest='nnCutHHL',    default=0.995,          help="nncut hadhad loose",   type=float)
    parser.add_argument('--nnCutHHMETL', dest='nnCutHHMETL', default=0.995,          help="nncut hadhad met loose",   type=float)
    parser.add_argument('--nnCutEH',     dest='nnCutEH',     default=0.95,           help="nncut elhad",   type=float)
    parser.add_argument('--nnCutMH',     dest='nnCutMH',     default=0.95,           help="nncut muhad",   type=float)
    parser.add_argument('--nnCutEHL',    dest='nnCutEHL',    default=0.5,            help="nncut elhad loose",   type=float)
    parser.add_argument('--nnCutMHL',    dest='nnCutMHL',    default=0.5,            help="nncut muhad loose",   type=float)
    parser.add_argument('--metCut',      dest='metCut',      default=50.,            help="met cut",   type=float)
    parser.add_argument('--hPtBinsHad',  dest='hPtBinsHad',  default=['400','None'], help="h pt bins (had)",   type=str,  nargs='+')
    parser.add_argument('--hPtBinsLep',  dest='hPtBinsLep',  default=['400','None'], help="h pt bins (lep)",   type=str,  nargs='+')
    parser.add_argument('--hPtCut',      dest='hPtCut',      default=None,           help="h pt cut",   type=float)
    parser.add_argument('--antiElCut',   dest='antiElCut',   default=0.5,            help="anti el cut",   type=float)
    parser.add_argument('--antiMuCut',   dest='antiMuCut',   default=0.5,            help="anti mu cut",   type=float)

    args = parser.parse_args()

    makeCards(args)
