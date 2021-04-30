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
import pickle
import ROOT
rl.util.install_roofit_helpers()
rl.ParametericSample.PreferRooParametricHist = False

import plot_stack

overflow_sum = 'allnan'

def expo_sample(norm, scale, obs):
    cdf = scipy.stats.expon.cdf(scale=scale, x=obs.binning) * norm
    return (np.diff(cdf), obs.binning, obs.name)


def gaus_sample(norm, loc, scale, obs):
    cdf = scipy.stats.norm.cdf(loc=loc, scale=scale, x=obs.binning) * norm
    return (np.diff(cdf), obs.binning, obs.name)

nnCut_hadhad = 0.995
nnCut_hadhad_loose = 0.9
nnCut_hadhad_met = 0.99
nnCut_hadhad_met_loose = 0.9
nnCut_hadel = 0.95
nnCut_hadel_loose = 0.1
nnCut_hadmu = 0.95
nnCut_hadmu_loose = 0.1

def getHist(h,var_name,lumifb,vars_cut,regionsel,blind,sigscale,rebin):
    print(h)
    exceptions = ['process', var_name]
    for var in vars_cut:
        exceptions.append(var)
    if (regionsel is not ''):
        exceptions.append('region')
    print([ax.name for ax in h.axes()])
    x = h.sum(*[ax for ax in h.axes() if ax.name not in exceptions],overflow=overflow_sum)
    for reg in regionsel:
        print('integrating ',reg)
        x = x.integrate('region',reg)
    for var,val in vars_cut.items():
        if len(val)==0:
            continue
        if var!=var_name:
            print('integrating ',var,val[0],val[1])
            x = x.integrate(var,slice(val[0],val[1]))
            #x = x.integrate(var,slice(val[0],val[1]),overflow=overflow)
    if var_name in vars_cut.keys():
        x = x[:, vars_cut[var_name][0]:vars_cut[var_name][1]]

    xaxis = var_name
    for ih,hkey in enumerate(x.identifiers('process')):
        x.identifiers('process')[ih].label = process_latex[hkey.name]

    if len(rebin)==1:
        if (rebin[0]>1): 
            x = x.rebin(xaxis,int(rebin[0]))
    else:
        x = x.rebin(xaxis,hist.Bin(xaxis, var_label, rebin))
 

    x_nobkg = x[nobkg]
    x_nosig = x[nosig]

    # normalize to lumi
    x_nosig.scale({p: lumifb for p in x_nosig.identifiers('process')}, axis="process")
    x_nobkg.scale({p: lumifb*float(sigscale) for p in x_nobkg.identifiers('process')}, axis="process")

    x_data = x['data']

    return x_nobkg+x_nosig+x_data

def createLepHad(sig_hist, top_cr_hist, wlnu_cr_hist, qcd_cr_hist, var_name, mttbins, ptbins, leptype, tmpdir, label, usingData):

    #could be made region-dependent

    samp_combinations = {
        'data_obs'   : ['data'],
        'top'        : ['tt-had', 'tt-semilep', 'tt-dilep','st'],
        'htt125'     : ['h125'],
        'multijet'   : ['qcd'],
        'ztt'        : ['zll'],
        'wlnu'       : ['wlnu'],
        'vqq'        : ['vv','vqq'],
        'ignore'     : [],
    }

    sig_faillep = sig_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('antilep',slice(-0.5,0.5))
    top_cr_faillep = top_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('antilep',slice(-0.5,0.5))
    wlnu_cr_faillep = wlnu_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('antilep',slice(-0.5,0.5))
    qcd_cr_faillep = qcd_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('antilep',slice(-0.5,0.5))

    sig_hist = sig_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('antilep',slice(0.5,1.5))
    top_cr_hist = top_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('antilep',slice(0.5,1.5))
    wlnu_cr_hist = wlnu_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('antilep',slice(0.5,1.5))
    qcd_cr_hist = qcd_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations).integrate('antilep',slice(0.5,1.5))

    jec = rl.NuisanceParameter('CMS_jec', 'lnN')
    top_norm = rl.NuisanceParameter('CMS_top_norm', 'lnN')
    wlnu_norm = rl.NuisanceParameter('CMS_wlnu_norm', 'lnN')
    multijet_norm = rl.NuisanceParameter('CMS_multijet_norm', 'lnN')
    ztt_norm = rl.NuisanceParameter('CMS_ztt_norm', 'lnN')
    vqq_norm = rl.NuisanceParameter('CMS_vqq_norm', 'lnN')
    trig = rl.NuisanceParameter('CMS_trig_had%s'%leptype, 'lnN')
    lumi = rl.NuisanceParameter('CMS_lumi', 'lnN')

    qcd_fail = rl.NuisanceParameter('multijet_Rfail', 'shape')
    qcd_loosepass = rl.NuisanceParameter('multijet_Rloosepass', 'shape')
    qcd_pass = rl.NuisanceParameter('multijet_Rpass', 'shape')
    qcdnormSF = rl.IndependentParameter('qcdnormSF', 1., 0, 10)

    topeffSF = rl.IndependentParameter('topeffSF', 1., 0, 10)
    topnormSF = rl.IndependentParameter('topnormSF', 1., 0, 10)
    wlnueffSF = rl.IndependentParameter('wlnueffSF', 1., 0, 10)
    wlnunormSF = rl.IndependentParameter('wlnunormSF', 1., 0, 10)

    topLeffSF = rl.IndependentParameter('topLeffSF', 1., 0, 10)
    wlnuLeffSF = rl.IndependentParameter('wlnuLeffSF', 1., 0, 10)

    rztt = rl.IndependentParameter('r_ztt', 1., 0, 10)
    ztt_eff = rl.IndependentParameter('ztt_eff', 1., 0, 10)

    npt = len(ptbins) - 1
    mtt = rl.Observable(var_name, mttbins)

    model = rl.Model("had%sModel"%leptype)

    if leptype=="el":
        nnCut_lep = nnCut_hadel
        nnCut_lep_loose = nnCut_hadel_loose
    elif leptype=="mu":
        nnCut_lep = nnCut_hadmu
        nnCut_lep_loose = nnCut_hadmu_loose
    else:
        print("Unknown leptype",leptype)
        return

    def intRegion(inhist,theregion):
        if theregion=='pass':
            theslice = slice(nnCut_lep,None)
        elif theregion=='loosepass':
            theslice = slice(nnCut_lep_loose,nnCut_lep)
        elif theregion=='fail':
            theslice = slice(None,nnCut_lep_loose)
        else:
            print("Unknown region",theregion)
            return

        #return inhist.integrate('nn_disc',theslice).integrate('jet_pt',slice(ptbins[0],ptbins[-1]))
        return inhist.integrate('nn_disc',theslice)

    qcd_data_faillep = intRegion(qcd_cr_faillep,"fail").integrate('sample','data_obs')
    qcd_data_faillep = qcd_data_faillep.values()[()]
    qcd_mc_faillep = intRegion(qcd_cr_faillep,"fail").integrate('sample',[s.name for s in qcd_cr_faillep.identifiers('sample') if s.name!='data_obs' and s.name!='multijet'])
    qcd_mc_faillep = qcd_mc_faillep.values()[()]
    qcd_data_faillep = qcd_data_faillep - qcd_mc_faillep
    qcd_data_faillep = np.array([v if v>0. else 1. for v in qcd_data_faillep])

    top_data_faillep = intRegion(top_cr_faillep,"fail").integrate('sample','data_obs')
    top_data_faillep = top_data_faillep.values()[()]
    top_mc_faillep = intRegion(top_cr_faillep,"fail").integrate('sample',[s.name for s in top_cr_faillep.identifiers('sample') if s.name!='data_obs' and s.name!='multijet'])
    top_mc_faillep = top_mc_faillep.values()[()]
    top_data_faillep = top_data_faillep - top_mc_faillep
    top_data_faillep = np.array([v if v>0. else 1. for v in top_data_faillep])

    wlnu_data_faillep = intRegion(wlnu_cr_faillep,"fail").integrate('sample','data_obs')
    wlnu_data_faillep = wlnu_data_faillep.values()[()]
    wlnu_mc_faillep = intRegion(wlnu_cr_faillep,"fail").integrate('sample',[s.name for s in wlnu_cr_faillep.identifiers('sample') if s.name!='data_obs' and s.name!='multijet'])
    wlnu_mc_faillep = wlnu_mc_faillep.values()[()]
    wlnu_data_faillep = wlnu_data_faillep - wlnu_mc_faillep
    wlnu_data_faillep = np.array([v if v>0. else 1. for v in wlnu_data_faillep])

    sig_data_faillep = intRegion(sig_faillep,"fail").integrate('sample','data_obs')
    sig_data_faillep = sig_data_faillep.values()[()]
    sig_mc_faillep = intRegion(sig_faillep,"fail").integrate('sample',[s.name for s in sig_faillep.identifiers('sample') if s.name!='data_obs' and s.name!='multijet'])
    sig_mc_faillep = sig_mc_faillep.values()[()]
    sig_data_faillep = sig_data_faillep - sig_mc_faillep
    sig_data_faillep = np.array([v if v>0. else 1. for v in sig_data_faillep])

    print(leptype)
    print(sig_data_faillep)
    print(top_data_faillep)
    print(wlnu_data_faillep)
    print(qcd_data_faillep)
    print('----------')

    qcdratio_sig = np.divide(sig_data_faillep, qcd_data_faillep)
    qcdratio_top = np.divide(top_data_faillep, qcd_data_faillep)
    qcdratio_wlnu = np.divide(wlnu_data_faillep, qcd_data_faillep)
    print('ratio sig ',qcdratio_sig)
    print('ratio top ',qcdratio_top)
    print('ratio wlnu',qcdratio_wlnu)

    qcd_from_data_tmp = intRegion(qcd_cr_hist['data_obs'],'fail').sum('sample').values()[()]
    mc_other = np.zeros_like(qcd_from_data_tmp)
    for s in qcd_cr_hist.identifiers('sample'):
        if (s.name!='data_obs' and s.name!='multijet'):
            mc_other += intRegion(qcd_cr_hist[s.name],'fail').sum('sample').values()[()]
    qcd_from_data = np.array([qcd_from_data_tmp[bi] - mc_other[bi] if qcd_from_data_tmp[bi] - mc_other[bi] > 0. else 0. for bi in range(len(qcd_from_data_tmp))])
    qcd_from_data_int = hist.poisson_interval(qcd_from_data,qcd_from_data)
    qcd_from_data_dn = np.array([qcd_from_data_int[1][bi] if qcd_from_data_int[1][bi] > 0. and not np.isnan(qcd_from_data_int[1][bi]) else 0. for bi in range(len(qcd_from_data_tmp))])
    qcd_from_data_up = np.array([qcd_from_data_int[0][bi] if qcd_from_data_int[0][bi] > 0. and not np.isnan(qcd_from_data_int[0][bi]) else 0. for bi in range(len(qcd_from_data_tmp))])
    qcd_from_data_dn = qcd_from_data_dn * np.sum(qcd_from_data)/np.sum(qcd_from_data_dn)
    qcd_from_data_up = qcd_from_data_up * np.sum(qcd_from_data)/np.sum(qcd_from_data_up)

    print('nom',qcd_from_data)
    print('dn ',qcd_from_data_dn)
    print('up ',qcd_from_data_up)

    qcdfail_temp = {}
    qcdfail_temp_dn = {}
    qcdfail_temp_up = {}

    for region in ['pass','loosepass','fail']:
        qcd_from_data_tmp = intRegion(qcd_cr_faillep['data_obs'],region).sum('sample').values()[()]
        mc_other = np.zeros_like(qcd_from_data_tmp)
        for s in qcd_cr_faillep.identifiers('sample'):
            if (s.name!='data_obs' and s.name!='multijet'):
                mc_other += intRegion(qcd_cr_faillep[s.name],region).sum('sample').values()[()]
        qcd_from_data_tmp = np.array([qcd_from_data_tmp[bi] - mc_other[bi] if qcd_from_data_tmp[bi] - mc_other[bi] > 1. else 1. for bi in range(len(qcd_from_data_tmp))])
        qcd_from_data_int = hist.poisson_interval(qcd_from_data_tmp,qcd_from_data_tmp)
        qcdfail_temp[region] = qcd_from_data_tmp
        qcdfail_temp_dn[region] = np.array([qcd_from_data_int[1][bi] if qcd_from_data_int[1][bi] > 0. and not np.isnan(qcd_from_data_int[1][bi]) else 0. for bi in range(len(qcd_from_data_tmp))])
        qcdfail_temp_up[region] = np.array([qcd_from_data_int[0][bi] if qcd_from_data_int[0][bi] > 0. and not np.isnan(qcd_from_data_int[0][bi]) else 0. for bi in range(len(qcd_from_data_tmp))])

    print('qcdfail',qcdfail_temp)
    print('qcdfail_dn',qcdfail_temp_dn)
    print('qcdfail_up',qcdfail_temp_up)

    qcdratio_F = {
        'pass':qcdfail_temp['pass']/qcdfail_temp['fail'],
        'loosepass':qcdfail_temp['loosepass']/qcdfail_temp['fail'],
        'fail':np.ones_like(qcd_from_data),
    }
    qcdratio_F_dn = {
        'pass':qcdfail_temp_dn['pass']/qcdfail_temp_dn['fail'],
        'loosepass':qcdfail_temp_dn['loosepass']/qcdfail_temp_dn['fail'],
        'fail':np.ones_like(qcd_from_data),
    }
    qcdratio_F_up = {
        'pass':qcdfail_temp_up['pass']/qcdfail_temp_up['fail'],
        'loosepass':qcdfail_temp_up['loosepass']/qcdfail_temp_up['fail'],
        'fail':np.ones_like(qcd_from_data),
    }
    for region in ['pass', 'loosepass', 'fail']:
        qcdratio_F[region] = np.sum(qcdratio_F[region]*qcd_from_data)/np.sum(qcd_from_data)
        qcdratio_F_dn[region] = np.sum(qcdratio_F_dn[region]*qcd_from_data)/np.sum(qcd_from_data)
        qcdratio_F_up[region] = np.sum(qcdratio_F_up[region]*qcd_from_data)/np.sum(qcd_from_data)
    print('qcdratio_F',qcdratio_F)
    print('qcdratio_F_dn',qcdratio_F_dn)
    print('qcdratio_F_up',qcdratio_F_up)

    #SR               - [F] [L] [P]
    #SR (noAnti)      - [F]
    #mIso CR          - [F] [L] [P]
    #mIso CR (noAnti) - [F] [L] [P]
    #derive ratio from {SR (noAnti) [F]} / {mIso CR (noAnti) [F]} ==> qcdratio

    #CR - SR (R_miso)
    #F - L   (RQCD_FL)
    #F - P   (RQCD_FP)

    #QCD[F] * R_miso = SR[F]
    #QCD[F] * R_miso * RQCD_FL = SR[L]
    #QCD[F] * R_miso * RQCD_FP = SR[L]

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
                if sName.name=='data_obs' and not usingData:
                    continue
                templ = (intRegion(thehist[sName],region).sum('sample').values()[()], mtt.binning, mtt.name)
                for iv,val in enumerate(templ[0]):
                    if val<0.:
                        templ[0][iv] = 0.
                if sName.name=='data_obs':
                    ch.setObservation(templ)
                else:
                    if sName.name=='multijet':
                        print('qcdpred',region,qcd_from_data * qcdratio_sig * np.mean(qcdratio_F[region]))
                        print('qcdpred_dn',region,qcd_from_data_dn * qcdratio_sig * np.mean(qcdratio_F_dn[region]))
                        print('qcdpred_up',region,qcd_from_data_up * qcdratio_sig * np.mean(qcdratio_F_up[region]))
                        templ = (qcd_from_data * qcdratio_sig * np.mean(qcdratio_F[region]), mtt.binning, mtt.name)
                    stype = rl.Sample.SIGNAL if sName.name == 'htt125' else rl.Sample.BACKGROUND
                    sample = rl.TemplateSample(ch.name + '_' + sName.name, stype, templ)
        
                    sample.setParamEffect(jec, 1.05)
                    if sName.name=='top':
                        sample.setParamEffect(top_norm, 1.05)
                    if sName.name=='wlnu':
                        sample.setParamEffect(wlnu_norm, 1.10)
                    if sName.name=='multijet':
                        sample.setParamEffect(multijet_norm, 1.20)
                        sample.setParamEffect(qcd_pass if isPass else qcd_loosepass if isLoosePass else qcd_fail, qcd_from_data_dn * qcdratio_sig * np.mean(qcdratio_F_dn[region]), qcd_from_data_up * qcdratio_sig * np.mean(qcdratio_F_up[region]))
                    if sName.name=='vqq':
                        sample.setParamEffect(vqq_norm, 1.10)
                    if sName.name=='ztt':
                        sample.setParamEffect(ztt_norm, 1.10)
                    sample.setParamEffect(trig, 1.02)
                    sample.setParamEffect(lumi, 1.024)
            
                    ch.addSample(sample)
            if not usingData:
                ch.setObservation((np.zeros(len(mtt.binning)-1),mtt.binning, mtt.name))
    

    for ptbin in range(npt):
        failCh = model['ptbin%dfail' % ptbin]
        passCh = model['ptbin%dpass' % ptbin]
        loosePassCh = model['ptbin%dloosepass' % ptbin]

        qcdpass = passCh['multijet']
        qcdloosepass = loosePassCh['multijet']
        qcdfail = failCh['multijet']
        qcdpass.setParamEffect(qcdnormSF, 1*qcdnormSF)
        qcdloosepass.setParamEffect(qcdnormSF, 1*qcdnormSF)
        qcdfail.setParamEffect(qcdnormSF, 1*qcdnormSF)

        toppass = passCh['top']
        toploosepass = loosePassCh['top']
        topfail = failCh['top']
        topPF = (toploosepass.getExpectation(nominal=True).sum() + toppass.getExpectation(nominal=True).sum()) / topfail.getExpectation(nominal=True).sum()
        topLPF = toppass.getExpectation(nominal=True).sum() / toploosepass.getExpectation(nominal=True).sum()
        topRLPF = 1./(1.+(1./topLPF)) # P/(L+P)
        toppass.setParamEffect(topLeffSF, 1*topLeffSF)
        toploosepass.setParamEffect(topLeffSF, (1 - topLeffSF) * topLPF + 1)
        toppass.setParamEffect(topeffSF, 1*topeffSF*topRLPF)
        toploosepass.setParamEffect(topeffSF, 1*topeffSF * (1 - topRLPF))
        topfail.setParamEffect(topeffSF, (1 - topeffSF) * topPF + 1)
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
        wlnupass.setParamEffect(wlnueffSF, 1*wlnueffSF*wlnuRLPF)
        wlnuloosepass.setParamEffect(wlnueffSF, 1*wlnueffSF * (1 - wlnuRLPF))
        wlnufail.setParamEffect(wlnueffSF, (1 - wlnueffSF) * wlnuPF + 1)
        wlnupass.setParamEffect(wlnunormSF, 1*wlnunormSF)
        wlnuloosepass.setParamEffect(wlnunormSF, 1*wlnunormSF)
        wlnufail.setParamEffect(wlnunormSF, 1*wlnunormSF)

        zttpass = passCh['ztt']
        zttloosepass = loosePassCh['ztt']
        zttfail = failCh['ztt']
        zttLP = zttpass.getExpectation(nominal=True).sum() / zttloosepass.getExpectation(nominal=True).sum()
        zttpass.setParamEffect(ztt_eff, 1*ztt_eff)
        zttloosepass.setParamEffect(ztt_eff, (1 - ztt_eff) * zttLP + 1)
        passCh['htt125'].setParamEffect(ztt_eff, 1*ztt_eff)
        loosePassCh['htt125'].setParamEffect(ztt_eff, (1 - ztt_eff) * zttLP + 1)
        zttpass.setParamEffect(rztt, 1*rztt)
        zttloosepass.setParamEffect(rztt, 1*rztt)
        zttfail.setParamEffect(rztt, 1*rztt)


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
            templ = (intRegion(thehist[sName],region).sum('sample').values()[()], mtt.binning, mtt.name)
            for iv,val in enumerate(templ[0]):
                if val<0.:
                    templ[0][iv] = 0.
            if sName.name=='data_obs':
                ch.setObservation(templ)
            else:
                if sName.name=='multijet':
                    templ = (qcd_from_data * qcdratio_top * np.mean(qcdratio_F[region]), mtt.binning, mtt.name)
                stype = rl.Sample.SIGNAL if sName.name == 'htt125' else rl.Sample.BACKGROUND
                sample = rl.TemplateSample(ch.name + '_' + sName.name, stype, templ)
    
                sample.setParamEffect(jec, 1.05)
                if sName.name=='top':
                    sample.setParamEffect(top_norm, 1.05)
                if sName.name=='wlnu':
                    sample.setParamEffect(wlnu_norm, 1.10)
                if sName.name=='multijet':
                    sample.setParamEffect(multijet_norm, 1.20)
                    sample.setParamEffect(qcd_pass if isPass else qcd_loosepass if isLoosePass else qcd_fail, qcd_from_data_dn * qcdratio_top * np.mean(qcdratio_F_dn[region]), qcd_from_data_up * qcdratio_top * np.mean(qcdratio_F_up[region]))
                if sName.name=='vqq':
                    sample.setParamEffect(vqq_norm, 1.10)
                if sName.name=='ztt':
                    sample.setParamEffect(ztt_norm, 1.10)
                sample.setParamEffect(trig, 1.02)
                sample.setParamEffect(lumi, 1.024)
        
                ch.addSample(sample)

    qcdpass = model['topCRpass']['multijet']
    qcdloosepass = model['topCRloosepass']['multijet']
    qcdfail = model['topCRfail']['multijet']
    qcdpass.setParamEffect(qcdnormSF, 1*qcdnormSF)
    qcdloosepass.setParamEffect(qcdnormSF, 1*qcdnormSF)
    qcdfail.setParamEffect(qcdnormSF, 1*qcdnormSF)

    toppass = model['topCRpass']['top']
    toploosepass = model['topCRloosepass']['top']
    topfail = model['topCRfail']['top']
    topPF = (toploosepass.getExpectation(nominal=True).sum() + toppass.getExpectation(nominal=True).sum()) / topfail.getExpectation(nominal=True).sum()
    topLPF = toppass.getExpectation(nominal=True).sum() / toploosepass.getExpectation(nominal=True).sum()
    topRLPF = 1./(1.+(1./topLPF)) # P/(L+P)
    toppass.setParamEffect(topLeffSF, 1*topLeffSF)
    toploosepass.setParamEffect(topLeffSF, (1 - topLeffSF) * topLPF + 1)
    toppass.setParamEffect(topeffSF, 1*topeffSF*topRLPF)
    toploosepass.setParamEffect(topeffSF, 1*topeffSF * (1 - topRLPF))
    topfail.setParamEffect(topeffSF, (1 - topeffSF) * topPF + 1)
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
    wlnupass.setParamEffect(wlnueffSF, 1*wlnueffSF*wlnuRLPF)
    wlnuloosepass.setParamEffect(wlnueffSF, 1*wlnueffSF * (1 - wlnuRLPF))
    wlnufail.setParamEffect(wlnueffSF, (1 - wlnueffSF) * wlnuPF + 1)
    wlnupass.setParamEffect(wlnunormSF, 1*wlnunormSF)
    wlnuloosepass.setParamEffect(wlnunormSF, 1*wlnunormSF)
    wlnufail.setParamEffect(wlnunormSF, 1*wlnunormSF)

    zttpass = model['topCRpass']['ztt']
    zttloosepass = model['topCRloosepass']['ztt']
    zttfail = model['topCRfail']['ztt']
    zttLP = zttpass.getExpectation(nominal=True).sum() / zttloosepass.getExpectation(nominal=True).sum()
    zttpass.setParamEffect(ztt_eff, 1*ztt_eff)
    zttloosepass.setParamEffect(ztt_eff, (1 - ztt_eff)* zttLP + 1)
    model['topCRpass']['htt125'].setParamEffect(ztt_eff, 1*ztt_eff)
    model['topCRloosepass']['htt125'].setParamEffect(ztt_eff, (1 - ztt_eff) * zttLP + 1)
    zttpass.setParamEffect(rztt, 1*rztt)
    zttloosepass.setParamEffect(rztt, 1*rztt)
    zttfail.setParamEffect(rztt, 1*rztt)

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
            templ = (intRegion(thehist[sName],region).sum('sample').values()[()], mtt.binning, mtt.name)
            for iv,val in enumerate(templ[0]):
                if val<0.:
                    templ[0][iv] = 0.
            if sName.name=='data_obs':
                ch.setObservation(templ)
            else:
                if sName.name=='multijet':
                    templ = (qcd_from_data * qcdratio_wlnu * np.mean(qcdratio_F[region]), mtt.binning, mtt.name)
                stype = rl.Sample.SIGNAL if sName.name == 'htt125' else rl.Sample.BACKGROUND
                sample = rl.TemplateSample(ch.name + '_' + sName.name, stype, templ)
    
                sample.setParamEffect(jec, 1.05)
                if sName.name=='top':
                    sample.setParamEffect(top_norm, 1.05)
                if sName.name=='wlnu':
                    sample.setParamEffect(wlnu_norm, 1.10)
                if sName.name=='multijet':
                    sample.setParamEffect(multijet_norm, 1.20)
                    sample.setParamEffect(qcd_pass if isPass else qcd_loosepass if isLoosePass else qcd_fail, qcd_from_data_dn * qcdratio_wlnu * np.mean(qcdratio_F_dn[region]), qcd_from_data_up * qcdratio_wlnu * np.mean(qcdratio_F_up[region]))
                if sName.name=='vqq':
                    sample.setParamEffect(vqq_norm, 1.10)
                if sName.name=='ztt':
                    sample.setParamEffect(ztt_norm, 1.10)
                sample.setParamEffect(trig, 1.02)
                sample.setParamEffect(lumi, 1.024)
        
                ch.addSample(sample)

    qcdpass = model['wlnuCRpass']['multijet']
    qcdloosepass = model['wlnuCRloosepass']['multijet']
    qcdfail = model['wlnuCRfail']['multijet']
    qcdpass.setParamEffect(qcdnormSF, 1*qcdnormSF)
    qcdloosepass.setParamEffect(qcdnormSF, 1*qcdnormSF)
    qcdfail.setParamEffect(qcdnormSF, 1*qcdnormSF)

    toppass = model['wlnuCRpass']['top']
    toploosepass = model['wlnuCRloosepass']['top']
    topfail = model['wlnuCRfail']['top']
    topPF = (toploosepass.getExpectation(nominal=True).sum() + toppass.getExpectation(nominal=True).sum()) / topfail.getExpectation(nominal=True).sum()
    topLPF = toppass.getExpectation(nominal=True).sum() / toploosepass.getExpectation(nominal=True).sum()
    topRLPF = 1./(1.+(1./topLPF)) # P/(L+P)
    toppass.setParamEffect(topLeffSF, 1*topLeffSF)
    toploosepass.setParamEffect(topLeffSF, (1 - topLeffSF) * topLPF + 1)
    toppass.setParamEffect(topeffSF, 1*topeffSF*topRLPF)
    toploosepass.setParamEffect(topeffSF, 1*topeffSF * (1 - topRLPF))
    topfail.setParamEffect(topeffSF, (1 - topeffSF) * topPF + 1)
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
    wlnupass.setParamEffect(wlnueffSF, 1*wlnueffSF*wlnuRLPF)
    wlnuloosepass.setParamEffect(wlnueffSF, 1*wlnueffSF * (1 - wlnuRLPF))
    wlnufail.setParamEffect(wlnueffSF, (1 - wlnueffSF) * wlnuPF + 1)
    wlnupass.setParamEffect(wlnunormSF, 1*wlnunormSF)
    wlnuloosepass.setParamEffect(wlnunormSF, 1*wlnunormSF)
    wlnufail.setParamEffect(wlnunormSF, 1*wlnunormSF)

    zttpass = model['wlnuCRpass']['ztt']
    zttloosepass = model['wlnuCRloosepass']['ztt']
    zttfail = model['wlnuCRfail']['ztt']
    zttLP = zttpass.getExpectation(nominal=True).sum() / zttloosepass.getExpectation(nominal=True).sum()
    zttpass.setParamEffect(ztt_eff, 1*ztt_eff)
    zttloosepass.setParamEffect(ztt_eff, (1 - ztt_eff) * zttLP + 1)
    model['wlnuCRpass']['htt125'].setParamEffect(ztt_eff, 1*ztt_eff)
    model['wlnuCRloosepass']['htt125'].setParamEffect(ztt_eff, (1 - ztt_eff) * zttLP + 1)
    zttpass.setParamEffect(rztt, 1*rztt)
    zttloosepass.setParamEffect(rztt, 1*rztt)
    zttfail.setParamEffect(rztt, 1*rztt)

    with open(os.path.join("%s/%s"%(str(tmpdir),label), 'had%sModel.pkl'%leptype), "wb") as fout:
        pickle.dump(model, fout, protocol=2)

    model.renderCombine(os.path.join("%s/%s"%(str(tmpdir),label), 'had%sModel'%leptype))

def createHadHad(sig_hist, sig_met_hist, top_cr_hist, var_name, mttbins, ptbins, tmpdir, label, usingData):

    #could be made region-dependent
    samp_combinations = {
        'data_obs'   : ['data'],
        'top'        : ['tt-had', 'tt-semilep', 'tt-dilep'],#,'st'],
        'htt125'     : ['h125'],
        #'multijet'  : ['qcd'],
        'ztt'        : ['zll'],
        'wlnu'       : ['wlnu'],
        'vqq'        : ['vv'],#'vqq'],
        'ignore'     : ['qcd','vqq','st'],
    }

    sig_hist = sig_hist.group("process", hist.Cat("sample", "sample"), samp_combinations)
    sig_met_hist = sig_met_hist.group("process", hist.Cat("sample", "sample"), samp_combinations)
    top_cr_hist = top_cr_hist.group("process", hist.Cat("sample", "sample"), samp_combinations)

    jec = rl.NuisanceParameter('CMS_jec', 'lnN')
    top_norm = rl.NuisanceParameter('CMS_top_norm', 'lnN')
    wlnu_norm = rl.NuisanceParameter('CMS_wlnu_norm', 'lnN')
    multijet_norm = rl.NuisanceParameter('CMS_multijet_norm', 'lnN')
    ztt_norm = rl.NuisanceParameter('CMS_ztt_norm', 'lnN')
    vqq_norm = rl.NuisanceParameter('CMS_vqq_norm', 'lnN')
    trig = rl.NuisanceParameter('CMS_trig_hadhad', 'lnN')
    lumi = rl.NuisanceParameter('CMS_lumi', 'lnN')
    topeffSF = rl.IndependentParameter('topeffSF', 1., 0, 10)
    topnormSF = rl.IndependentParameter('topnormSF', 1., 0, 10)

    rztt = rl.IndependentParameter('r_ztt', 1., 0, 10)

    npt = len(ptbins) - 1
    mtt = rl.Observable(var_name, mttbins)

    model = rl.Model("hadhadModel")

    def intRegion(inhist,theregion):
        if theregion=='pass':
            theslice = slice(nnCut_hadhad,None)
        elif theregion=='loosepass':
            theslice = slice(nnCut_hadhad_loose,nnCut_hadhad)
        elif theregion=='fail':
            theslice = slice(None,nnCut_hadhad_loose)
        else:
            print("Unknown region",theregion)
            return

        return inhist.integrate('nn_disc',theslice)

    for ptbin in range(npt):
        for region in ['pass','fail']:
            ch = rl.Channel("ptbin%d%s" % (ptbin, region))
            model.addChannel(ch)

            isPass = region=='pass'

            thehist = sig_hist if ptbin>0 else sig_met_hist
    
            for sName in thehist.identifiers('sample'):
                if sName.name=='ignore':
                    continue
                if sName.name=='data_obs' and not usingData:
                    continue
                templ = (intRegion(thehist[sName],region).sum('sample').values()[()], mtt.binning, mtt.name)
                for iv,val in enumerate(templ[0]):
                    if val<0.:
                        templ[0][iv] = 0.
                if sName.name=='data_obs':
                    ch.setObservation(templ)
                else:
                    stype = rl.Sample.SIGNAL if sName.name == 'htt125' else rl.Sample.BACKGROUND
                    sample = rl.TemplateSample(ch.name + '_' + sName.name, stype, templ)
        
                    sample.setParamEffect(jec, 1.05)
                    if sName.name=='top':
                        sample.setParamEffect(top_norm, 1.05)
                    if sName.name=='wlnu':
                        sample.setParamEffect(wlnu_norm, 1.10)
                    if sName.name=='multijet':
                        sample.setParamEffect(multijet_norm, 1.20)
                    if sName.name=='vqq':
                        sample.setParamEffect(vqq_norm, 1.10)
                    if sName.name=='ztt':
                        sample.setParamEffect(ztt_norm, 1.10)
                    sample.setParamEffect(trig, 1.02)
                    sample.setParamEffect(lumi, 1.024)
            
                    ch.addSample(sample)
            if not usingData:
                ch.setObservation((np.zeros(len(mtt.binning)-1),mtt.binning, mtt.name))
    

    for ptbin in range(npt):
        failCh = model['ptbin%dfail' % ptbin]
        passCh = model['ptbin%dpass' % ptbin]

        toppass = passCh['top']
        topfail = failCh['top']
        topPF = toppass.getExpectation(nominal=True).sum() / topfail.getExpectation(nominal=True).sum()
        toppass.setParamEffect(topeffSF, 1*topeffSF)
        topfail.setParamEffect(topeffSF, (1 - topeffSF) * topPF + 1)
        toppass.setParamEffect(topnormSF, 1*topnormSF)
        topfail.setParamEffect(topnormSF, 1*topnormSF)

        zttpass = passCh['ztt']
        zttfail = failCh['ztt']
        zttpass.setParamEffect(rztt, 1*rztt)
        zttfail.setParamEffect(rztt, 1*rztt)


    # Fill in muon CR
    for region in ['pass', 'fail']:
        ch = rl.Channel("muonCR%s" % (region, ))
        model.addChannel(ch)

        isPass = region=='pass'

        for sName in top_cr_hist.identifiers('sample'):
            if sName.name=='ignore':
                continue
            if sName.name=='data_obs' and not usingData:
                continue
            templ = (intRegion(top_cr_hist[sName],region).sum('sample').values()[()], mtt.binning, mtt.name)
            for iv,val in enumerate(templ[0]):
                if val<0.:
                    templ[0][iv] = 0.
            if sName.name=='data_obs':
                ch.setObservation(templ)
            else:
                stype = rl.Sample.SIGNAL if sName.name == 'htt125' else rl.Sample.BACKGROUND
                sample = rl.TemplateSample(ch.name + '_' + sName.name, stype, templ)
    
                sample.setParamEffect(jec, 1.05)
                if sName.name=='top':
                    sample.setParamEffect(top_norm, 1.05)
                if sName.name=='wlnu':
                    sample.setParamEffect(wlnu_norm, 1.10)
                if sName.name=='multijet':
                    sample.setParamEffect(multijet_norm, 1.20)
                if sName.name=='vqq':
                    sample.setParamEffect(vqq_norm, 1.10)
                if sName.name=='ztt':
                    sample.setParamEffect(ztt_norm, 1.10)
                sample.setParamEffect(trig, 1.02)
                sample.setParamEffect(lumi, 1.024)
        
                ch.addSample(sample)
        if not usingData:
            ch.setObservation((np.zeros(len(mtt.binning)-1),mtt.binning, mtt.name))

    toppass = model['muonCRpass']['top']
    topfail = model['muonCRfail']['top']
    topPF = toppass.getExpectation(nominal=True).sum() / topfail.getExpectation(nominal=True).sum()
    toppass.setParamEffect(topeffSF, 1*topeffSF)
    topfail.setParamEffect(topeffSF, (1 - topeffSF) * topPF + 1)
    toppass.setParamEffect(topnormSF, 1*topnormSF)
    topfail.setParamEffect(topnormSF, 1*topnormSF)

    zttpass = model['muonCRpass']['ztt']
    zttfail = model['muonCRfail']['ztt']
    zttpass.setParamEffect(rztt, 1*rztt)
    zttfail.setParamEffect(rztt, 1*rztt)

    with open(os.path.join("%s/%s"%(str(tmpdir),label), 'hadhadModel.pkl'), "wb") as fout:
        pickle.dump(model, fout, protocol=2)

    model.renderCombine(os.path.join("%s/%s"%(str(tmpdir),label), 'hadhadModel'))

def getHist(h,var_name,lumifb,vars_cut,regionsel,blind,sigscale,rebin):
    print(h)
    exceptions = ['process', var_name]
    for var in vars_cut:
        exceptions.append(var)
    if (regionsel is not ''):
        exceptions.append('region')
    print([ax.name for ax in h.axes()])
    x = h.sum(*[ax for ax in h.axes() if ax.name not in exceptions],overflow=overflow_sum)
    for reg in regionsel:
        print('integrating ',reg)
        x = x.integrate('region',reg)
    for var,val in vars_cut.items():
        if len(val)==0:
            continue
        if var!=var_name:
            print('integrating ',var,val[0],val[1])
            x = x.integrate(var,slice(val[0],val[1]))
            #x = x.integrate(var,slice(val[0],val[1]),overflow=overflow)
    if var_name in vars_cut.keys():
        print('integrating ',var_name,vars_cut[var_name][0],vars_cut[var_name][1])
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

    mcSamples = [
      #"QCD", 
      #"TT", 
      #"ST", 
      #"WJetsToLNu", 
      #"DYJetsToLL", 
      #"VJetsToQQ", 
      #"HTauTau",
      "DYJetsToLL_Pt-250To400_TuneCP5_13TeV-amcatnloFXFX-pythia8",
      "DYJetsToLL_Pt-400To650_TuneCP5_13TeV-amcatnloFXFX-pythia8",
      "DYJetsToLL_Pt-650ToInf_TuneCP5_13TeV-amcatnloFXFX-pythia8",
      "QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8",
      "QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8",
      "QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8",
      "QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8",
      "QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8",
      "QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8",
      "ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8",
      "ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8",
      "ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8",
      "ST_t-channel_top_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8",
      "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8",
      "ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8",
      "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8",
      "TTToHadronic_TuneCP5_13TeV-powheg-pythia8",
      "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8",
      "WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8",
      "WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8",
      "WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8",
      "WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8",
      "WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8",
      "WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8",
      "WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8",
      "WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8",
      "WJetsToQQ_HT400to600_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8",
      "WJetsToQQ_HT600to800_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8",
      "WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8",
      "ZJetsToQQ_HT400to600_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8",
      "ZJetsToQQ_HT600to800_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8",
      "ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8",
      "GluGluHToTauTau_M125_13TeV_powheg_pythia8",
      "VBFHToTauTau_M125_13TeV_powheg_pythia8",
      "WminusHToTauTau_M125_13TeV_powheg_pythia8",
      "WplusHToTauTau_M125_13TeV_powheg_pythia8",
      "ZHToTauTau_M125_13TeV_powheg_pythia8",
      "ggZH_HToTauTau_ZToLL_M125_13TeV_powheg_pythia8",
      "ggZH_HToTauTau_ZToNuNu_M125_13TeV_powheg_pythia8",
      "ggZH_HToTauTau_ZToQQ_M125_13TeV_powheg_pythia8",
      "ttHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8",
    ]

    dataSamples = [
      #"SingleMuon", 
      #"JetHT", 
      #"SingleElectron",
      #"MET",
      "SingleElectron_pancakes-02-withPF_Run2017C-09Aug2019_UL2017-v1",
      "SingleElectron_pancakes-02-withPF_Run2017D-09Aug2019_UL2017-v1",
      "SingleElectron_pancakes-02-withPF_Run2017E-09Aug2019_UL2017-v1",
      "SingleElectron_pancakes-02-withPF_Run2017F-09Aug2019_UL2017_rsb-v2",
      "JetHT_pancakes-02_Run2017C-09Aug2019_UL2017-v1",
      "JetHT_pancakes-02_Run2017D-09Aug2019_UL2017-v1",
      "JetHT_pancakes-02_Run2017E-09Aug2019_UL2017-v1",
      "JetHT_pancakes-02_Run2017F-09Aug2019_UL2017-v1",
      "SingleMuon_pancakes-02-withPF_Run2017C-09Aug2019_UL2017-v1",
      "SingleMuon_pancakes-02-withPF_Run2017D-09Aug2019_UL2017-v1",
      "SingleMuon_pancakes-02-withPF_Run2017E-09Aug2019_UL2017-v1",
      "SingleMuon_pancakes-02-withPF_Run2017F-09Aug2019_UL2017-v1",
      "MET_pancakes-02-withPF_Run2017C-09Aug2019_UL2017_rsb-v1",
      "MET_pancakes-02-withPF_Run2017D-09Aug2019_UL2017_rsb-v1",
      "MET_pancakes-02-withPF_Run2017E-09Aug2019_UL2017_rsb-v1",
      "MET_pancakes-02-withPF_Run2017F-09Aug2019_UL2017_rsb-v1",
    ]
    dataSamples = []

    tag = args.tag

    odir = 'cards/%s/'%tag
    os.system('mkdir -p %s/%s'%(odir,args.label))
    pwd = os.getcwd()

    hists_mapped = {}
    for h in mcSamples + dataSamples:
        print(h)
        # open hists
        hists_unmapped = load('%s%s.coffea'%(args.hist,h))
        # map to hists
        for key, val in hists_unmapped.items():
            if isinstance(val, hist.Hist):
                if key in hists_mapped:
                    hists_mapped[key] = hists_mapped[key] + processmap.apply(val)
                else:
                    hists_mapped[key] = processmap.apply(val)

    includeData = False
    # properties
    region_dict = {
        'hadhad_signal':['met_nn_kin', 'massreg', {'massreg':[40.,210.],'nn_disc':[],'jetmet_dphi':[None,1.6]}, includeData, [1], 'hadhad_signal'],
        'hadhad_signal_met':['met_nn_kin', 'massreg', {'massreg':[40.,210.],'nn_disc':[],'jetmet_dphi':[None,1.6], 'met_pt':[200.,None]}, includeData, [1], 'hadhad_signal_met'],
        'hadhad_cr_mu':['met_nn_kin', 'massreg', {'massreg':[40.,210.],'nn_disc':[],'jetmet_dphi':[None,1.6]}, includeData, [1], 'hadhad_cr_mu'],
        'hadel_signal':['met_nn_kin', 'massreg', {'massreg':[40.,210.],'nn_disc':[],'antilep':[]}, includeData, [1], 'hadel_signal'],
        'hadmu_signal':['met_nn_kin', 'massreg', {'massreg':[40.,210.],'nn_disc':[],'antilep':[]}, includeData, [1], 'hadmu_signal'],
        'hadel_cr_b':['met_nn_kin', 'massreg', {'massreg':[40.,210.],'nn_disc':[],'antilep':[]}, includeData, [1], 'hadel_cr_b'],
        'hadmu_cr_b':['met_nn_kin', 'massreg', {'massreg':[40.,210.],'nn_disc':[],'antilep':[]}, includeData, [1], 'hadmu_cr_b'],
        'hadel_cr_w':['met_nn_kin', 'massreg', {'massreg':[40.,210.],'nn_disc':[],'antilep':[]}, includeData, [1], 'hadel_cr_w'],
        'hadmu_cr_w':['met_nn_kin', 'massreg', {'massreg':[40.,210.],'nn_disc':[],'antilep':[]}, includeData, [1], 'hadmu_cr_w'],
        'hadel_cr_qcd':['met_nn_kin', 'massreg', {'massreg':[40.,210.],'nn_disc':[],'antilep':[]}, includeData, [1], 'hadel_cr_qcd'],
        'hadmu_cr_qcd':['met_nn_kin', 'massreg', {'massreg':[40.,210.],'nn_disc':[],'antilep':[]}, includeData, [1], 'hadmu_cr_qcd'],
    }

    full_dict = {}
    for reg in region_dict:
        reg_opt = region_dict[reg]
        h = hists_mapped[reg_opt[0]]
        print(h)
        xhist = getHist(h,reg_opt[1],args.lumi,reg_opt[2],[reg_opt[5]],reg_opt[3],args.sigscale,reg_opt[4])
        print(xhist)
        full_dict[reg] = xhist

    mttbins = np.linspace(40, 210, 18)
    hadptbins = np.array([300, 450, 1200])
    lepptbins = np.array([300, 1200])

    createHadHad(full_dict['hadhad_signal'], full_dict['hadhad_signal_met'], full_dict['hadhad_cr_mu'], 'massreg', mttbins, hadptbins, odir, args.label, args.unblind)
    createLepHad(full_dict['hadel_signal'], full_dict['hadel_cr_b'], full_dict['hadel_cr_w'], full_dict['hadel_cr_qcd'], 'massreg', mttbins, lepptbins, "el", odir, args.label, args.unblind)
    createLepHad(full_dict['hadmu_signal'], full_dict['hadmu_cr_b'], full_dict['hadmu_cr_w'], full_dict['hadmu_cr_qcd'], 'massreg', mttbins, lepptbins, "mu", odir, args.label, args.unblind)

    if args.plots:
    
        pwd = os.getcwd()
        os.chdir('%s/%s'%(odir,args.label))
    
        for rn in ['hadhad_signal','hadhad_signal_met','hadhad_cr_mu']:
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
          plot_stack.drawStack(full_dict[rn],'HadhadModel','massreg',r"$m_{NN}$",'%s pass'%rn,'all',1.,{'nn_disc':[nncutt,None]},{},[],rn+'_pass', blind='' if args.unblind or 'signal' not in rn else [40.,110.])
        for lep in ['el','mu']:
          if lep=='el':
            nncutl = nnCut_hadel_loose
            nncutt = nnCut_hadel
          else:
            nncutl = nnCut_hadmu_loose
            nncutt = nnCut_hadmu
          for rn in ['had%s_signal'%lep,'had%s_cr_b'%lep,'had%s_cr_w'%lep,'had%s_cr_qcd'%lep]:
            plot_stack.drawStack(full_dict[rn],'Had%sModel'%lep,'massreg',r"$m_{NN}$",'%s fail'%rn,'all',1.,{'nn_disc':[None,nncutl]},{},[],rn+'_fail')
            plot_stack.drawStack(full_dict[rn],'Had%sModel'%lep,'massreg',r"$m_{NN}$",'%s loosepass'%rn,'all',1.,{'nn_disc':[nncutl,nncutt]},{},[],rn+'_loosepass')
            plot_stack.drawStack(full_dict[rn],'Had%sModel'%lep,'massreg',r"$m_{NN}$",'%s pass'%rn,'all',1.,{'nn_disc':[nncutt,None]},{},[],rn+'_pass', blind='' if args.unblind or 'signal' not in rn else [40.,110.])
    
        os.chdir(pwd)

if __name__ == "__main__":

    if not sys.warnoptions:
        import warnings
        warnings.simplefilter("ignore")
        os.environ["PYTHONWARNINGS"] = "ignore" # Also affect subprocesses

    #ex. python makeCards.py --hist ../../boostedhiggs/condor/Mar08_NN/hists_sum_ --year 2017 --lumi 36.7 --tag Mar08 --label 01
    parser = argparse.ArgumentParser()
    parser.add_argument('--hist',       dest='hist',       default="hists_sum_", help="hists pickle prefix")
    parser.add_argument('--year',       dest='year',       default="2017",       help="year")
    parser.add_argument('--lumi',       dest='lumi',       default=50.,          help="lumi",       type=float)
    parser.add_argument('--tag',        dest='tag',        default="",           help="tag")
    parser.add_argument('--label',      dest='label',      default="",           help="label")
    parser.add_argument('--sigscale',   dest='sigscale',   default=1.,           help="sigscale",   type=float)
    parser.add_argument('--unblind',    dest='unblind',    action="store_true",  help="unblind")
    parser.add_argument('--plots',      dest='plots',      action="store_true",  help="plots")
    args = parser.parse_args()

    makeCards(args)
