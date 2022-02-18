from __future__ import print_function, division
import gzip
import json
import os
import sys
import glob

import uproot
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.styles.ROOT)
import numpy as np

from coffea import hist

import pickle
import gzip
import math
from cycler import cycler

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
err_opts_data = {
    'marker': '.',
    'markersize': 10.,
    'color':'k',
    'elinewidth': 1,
    #'emarker': '-'
}
err_opts_denom = {
    'facecolor':'gray', 
    'alpha':0.6,
}
line_opts = {
    'linewidth':2,
    #'marker':'None',
    'linestyle':'dashed',
    #'drawstyle':'steps'
}

overflow_sum = 'all'
#overflow_sum = 'allnan'

def drawStack(h,sel,var_name,var_label,plottitle,sample,lumifb,vars_cut,sig_cut,regionsel,savename,xlimits='',blind='',solo=False,sigscale=50.,sigstack=False,noratio=False,dosigrat=False,rebin=[1],overflow='none',xlog=False,xexp=False,norm=None,density=False,docomp=False,comp_var=[],comp_cut=[],complabels=[],colormap=True,verbose=False,signame=[],systematic='nominal',qcd_hists=[],qcd_cut_shape={},qcd_cut_num={},qcd_cut_den={}):
    exceptions = ['process', var_name]
    qcd_shape_exceptions = ['process', var_name]
    qcd_num_exceptions = ['process']
    qcd_den_exceptions = ['process']
    if systematic != '':
        exceptions.append('systematic')
        qcd_shape_exceptions.append('systematic')
        qcd_num_exceptions.append('systematic')
        qcd_den_exceptions.append('systematic')
    for var in vars_cut:
        exceptions.append(var)
    for var in qcd_cut_shape:
        qcd_shape_exceptions.append(var)
    for var in qcd_cut_num:
        qcd_num_exceptions.append(var)
    for var in qcd_cut_den:
        qcd_den_exceptions.append(var)
    for var in sig_cut:
        exceptions.append(var)
    if (regionsel != ''):
        exceptions.append('region')
    if (docomp):
        exceptions.extend(comp_var)
    if verbose: print([ax.name for ax in h.axes()])
    x = h.sum(*[ax for ax in h.axes() if ax.name not in exceptions],overflow=overflow_sum)
    if len(qcd_hists)==3:
        qcd_shape = qcd_hists[0].sum(*[ax for ax in qcd_hists[0].axes() if ax.name not in qcd_shape_exceptions],overflow=overflow_sum)
        qcd_num = qcd_hists[1].sum(*[ax for ax in qcd_hists[1].axes() if ax.name not in qcd_num_exceptions],overflow=overflow_sum)
        qcd_den = qcd_hists[2].sum(*[ax for ax in qcd_hists[2].axes() if ax.name not in qcd_den_exceptions],overflow=overflow_sum)
    if systematic != '':
        x = x.integrate('systematic',systematic)
        if len(qcd_hists)==3:
            qcd_shape = qcd_shape.integrate('systematic',systematic)
            qcd_num = qcd_num.integrate('systematic',systematic)
            qcd_den = qcd_den.integrate('systematic',systematic)
    for reg in regionsel:
        if verbose: print('integrating ',reg)
        x = x.integrate('region',reg)
    for var,val in vars_cut.items():
        if var!=var_name:
            if verbose: print('integrating ',var,val[0],val[1])
            x = x.integrate(var,slice(float(val[0]) if val[0] is not None else None,float(val[1]) if val[1] is not None else None),overflow='none' if val[0] is not None and val[1] is not None else 'under' if val[0] is None and val[1] is not None else 'over' if val[0] is not None and val[1] is None else 'allnan')
            #x = x.integrate(var,slice(val[0],val[1]),overflow=overflow)
    if var_name in vars_cut.keys():
        x = x[:, float(vars_cut[var_name][0]) if vars_cut[var_name][0] is not None else None:float(vars_cut[var_name][1]) if vars_cut[var_name][1] is not None else None]

    if len(qcd_hists)==3:
        for var,val in qcd_cut_shape.items():
            if var!=var_name:
                if verbose: print('integrating ',var,val[0],val[1])
                qcd_shape = qcd_shape.integrate(var,slice(float(val[0]) if val[0] is not None else None,float(val[1]) if val[1] is not None else None),overflow='none' if val[0] is not None and val[1] is not None else 'under' if val[0] is None and val[1] is not None else 'over' if val[0] is not None and val[1] is None else 'allnan')
        if var_name in qcd_cut_shape.keys():
            qcd_shape = qcd_shape[:, float(qcd_cut_shape[var_name][0]) if qcd_cut_shape[var_name][0] is not None else None:float(qcd_cut_shape[var_name][1]) if qcd_cut_shape[var_name][1] is not None else None]
    
        for var,val in qcd_cut_num.items():
            if verbose: print('integrating ',var,val[0],val[1])
            qcd_num = qcd_num.integrate(var,slice(float(val[0]) if val[0] is not None else None,float(val[1]) if val[1] is not None else None),overflow='none' if val[0] is not None and val[1] is not None else 'under' if val[0] is None and val[1] is not None else 'over' if val[0] is not None and val[1] is None else 'allnan')
    
        for var,val in qcd_cut_den.items():
            if verbose: print('integrating ',var,val[0],val[1])
            qcd_den = qcd_den.integrate(var,slice(float(val[0]) if val[0] is not None else None,float(val[1]) if val[1] is not None else None),overflow='none' if val[0] is not None and val[1] is not None else 'under' if val[0] is None and val[1] is not None else 'over' if val[0] is not None and val[1] is None else 'allnan')

    xaxis = var_name
    x.axis(xaxis).label = var_label
    for ih,hkey in enumerate(x.identifiers('process')):
        x.identifiers('process')[ih].label = process_latex[hkey.name]

    x.axis('process').sorting = 'integral'
    if (noratio and not dosigrat): fig,ax = plt.subplots()
    elif (not noratio and dosigrat): fig, (ax,axr,axs) = plt.subplots(3, 1, sharex='col', gridspec_kw={'height_ratios': [4, 1, 1],'hspace': 0.1})
    elif (noratio and dosigrat): fig, (ax,axs) = plt.subplots(2, 1, sharex='col', gridspec_kw={'height_ratios': [4, 1],'hspace': 0.1})
    else: fig, (ax,axr) = plt.subplots(2, 1, sharex='col', gridspec_kw={'height_ratios': [4, 1],'hspace': 0.1})

    if len(rebin)==1:
        if (rebin[0]>1): 
            x = x.rebin(xaxis,int(rebin[0]))
    else:
        x = x.rebin(xaxis,hist.Bin(xaxis, var_label, rebin))
 
    if (solo):
        if (sample=='sig'):
            x = x[nobkg]
        elif (sample=='bkg'):
            x = x[nosig]
        elif (sample!='all'):
            x = x[sample]

        hist.plot1d(x,
                overlay='process',ax=ax,
                clear=False,
                stack=False,
                fill_opts=fill_opts,
                error_opts=err_opts,
                density=density,
                binwnorm=norm,
                overflow=overflow
                )

    elif (docomp):
        if (sample=='sig'):
            x = x[nobkg]
        elif (sample=='bkg'):
            x = x[nosig]
        elif (sample!='all'):
            x = x[sample]

        x = x.sum(*[ax for ax in h.axes() if ax.name=='process'])

        comp_axis = hist.Cat('compsel', 'Selections')
        comp_map = {}
        re_ind = []
        for a in x.axes():
            if a.name in comp_var:
                re_ind.append(comp_var.index(a.name))
        comp_var = [comp_var[r] for r in re_ind]
        comp_cut = [comp_cut[r] for r in re_ind]
        complabels = [complabels[r] for r in re_ind]
        clis = [0 for cv in comp_var]
        cmax = [len(c) for c in comp_cut]
        while clis[-1] < cmax[-1]:
            tmpx = x
            for ic in range(len(comp_cut)):
                tmpx = tmpx.integrate(comp_var[ic],slice(comp_cut[ic][clis[ic]][0],comp_cut[ic][clis[ic]][1]))
            slicesum = False
            for key,val in tmpx.values().items():
                if val.sum() > 0.:
                    slicesum = True
            if slicesum:
                comp_map[', '.join([complabels[ic][clis[ic]] for ic in range(len(comp_cut))])] = tuple([slice(comp_cut[ic][clis[ic]][0],comp_cut[ic][clis[ic]][1]) for ic in range(len(comp_cut))])
            clis[0] = clis[0] + 1
            for ic in range(len(clis)-1):
                if clis[ic] >= cmax[ic]:
                    clis[ic+1] = clis[ic+1] + 1
                    clis[ic] = 0
        #{'new_bin': (slice, ...), ...}
        if verbose: print(comp_map)
        x = x.group(old_axes=tuple(comp_var), new_axis=comp_axis, mapping=comp_map)
        if verbose: print(x.values())
        hist.plot1d(x,
                overlay='compsel',ax=ax,
                clear=False,
                stack=False,
                line_opts={
                    'linewidth':2,
                    'drawstyle':'steps'
                },
                density=density,
                binwnorm=norm,
                overflow=overflow
                )

    else:

        x_nobkg = x[nobkg]
        if len(signame)>0:
            x_nobkg = x[nobkg][signame[0]]
            for sign in range(1,len(signame)):
                x_nobkg.add(x[nobkg][signame[sign]])

        for var,val in sig_cut.items():
            if var!=var_name:
                if verbose: print('integrating ',var,val[0],val[1])
                x_nobkg = x_nobkg.integrate(var,slice(val[0],val[1]),overflow='none' if val[0] is not None and val[1] is not None else 'under' if val[0] is None and val[1] is not None else 'over' if val[0] is not None and val[1] is None else 'allnan')
                #x = x.integrate(var,slice(val[0],val[1]),overflow=overflow)
        if var_name in sig_cut.keys():
            x_nobkg = x_nobkg[:, float(sig_cut[var_name][0]) if sig_cut[var_name][0] is not None else None:float(sig_cut[var_name][1]) if sig_cut[var_name][1] is not None else None]
        x = x.sum(*[ax for ax in x.axes() if ax.name in sig_cut],overflow='allnan')

        x_nosig = x[nosig]

        # normalize to lumi
        x_nosig.scale({p: lumifb for p in x_nosig.identifiers('process')}, axis="process")
        all_bkg = 0
        sampcount = 0
        samp_order = []
        for key,val in x_nosig.values().items():
            sampcount = sampcount + 1
            samp_sum = val.sum()
            all_bkg+=samp_sum
            samp_order.append([str(key)[2:-3],samp_sum])
    
        x_nobkg.scale({p: lumifb*float(sigscale) for p in x_nobkg.identifiers('process')}, axis="process")

        if len(qcd_hists)==3:

            qcd_shape_data = qcd_shape['data'].sum('process')
            qcd_shape = qcd_shape[noqcd].sum('process')
            qcd_shape.scale(lumifb)
            qcd_shape = qcd_shape_data.values()[()] - qcd_shape.values()[()]
            qcd_num_data = qcd_num['data'].sum('process')
            qcd_num = qcd_num[noqcd].sum('process')
            qcd_num.scale(lumifb)
            qcd_num = qcd_num_data.values()[()] - qcd_num.values()[()]
            qcd_den_data = qcd_den['data'].sum('process')
            qcd_den = qcd_den[noqcd].sum('process')
            qcd_den.scale(lumifb)
            qcd_den = qcd_den_data.values()[()] - qcd_den.values()[()]
            qcd_ratio = qcd_num/qcd_den if qcd_den!=0. else qcd_num
            if qcd_ratio<0.: qcd_ratio = -qcd_ratio

            x_qcd = x_nosig['qcd'].copy()
            x_qcd.clear()
            x_nosig = x_nosig.remove([iq for iq in x_nosig.axis('process').identifiers() if iq.name=='qcd'],'process')
            qcdfill = np.array(sum([[qcd_shape_data.axis(var_name).centers()[ibin]]*int(qcd_shape_data.values()[()][ibin]) for ibin in range(len(qcd_shape))],[]))
            qcdw = np.array(sum([[qcd_shape[ibin]/(qcd_shape_data.values()[()][ibin]) if qcd_shape[ibin]>0. else 0.]*int(qcd_shape_data.values()[()][ibin]) for ibin in range(len(qcd_shape))],[]))
            print(qcd_shape)
            print(np.mean(qcdw))
            print(qcd_ratio)
            x_qcd.fill(**{'process':'qcd',var_name:qcdfill,'weight':qcdw*qcd_ratio})
            x_nosig = x_nosig.add(x_qcd)
    
        if verbose: print(x_nobkg.identifiers('process'))
    
        sig_order = []
        all_sig = 0
        for key,val in x_nobkg.values().items():
            samp_sum = val.sum()
            all_sig+=samp_sum
            sig_order.append([str(key)[2:-3],samp_sum])
            if verbose: print(key,val)

        samp_order.sort(key=lambda x : x[1],reverse=True)
        sig_order.sort(key=lambda x : x[1],reverse=False)

        if verbose: print(samp_order)
        if verbose: print(sig_order)

        the_order = []
        for lab in samp_order:
            for ident in x_nosig.identifiers('process'):
                if (lab[0] == ident.name):
                    the_order.insert(0,ident)
                    break
        
        sig_color_order = []

        if colormap:
            color_order = [color_map[s[0]] for s in samp_order]
            if (sigstack):
                color_order = [color_map[s[0]] for s in sig_order] + color_order
            else:
                sig_color_order = [sig_color_map[s[0]] for s in sig_order]
                
            custom_cycler = (cycler(color=color_order))
            sig_cycler = (cycler(color=sig_color_order))
            ax.set_prop_cycle(custom_cycler)

        if verbose: print(sig_color_order)

        if (all_bkg>0.): hist.plot1d(x_nosig,
                    overlay='process',ax=ax,
                    clear=False,
                    stack=(sampcount>1),
                    order=the_order,
                    fill_opts=fill_opts,
                    error_opts=err_opts,
                    overflow=overflow
                    )
        if (all_sig>0.): 
            if colormap:
                ax.set_prop_cycle(sig_cycler)
            hist.plot1d(x_nobkg,ax=ax,
                    overlay='process',
                    clear=False,
                    line_opts=line_opts,
                    overflow=overflow)
        if verbose: print('-------> ',ax.get_legend_handles_labels()[1])

        if (sigstack):
            old_handles, old_labels = ax.get_legend_handles_labels()
            ax.cla()
            x_comb = x_nobkg+x_nosig
            #x.axis('process').sorting = 'integral'
            the_order = []
            if colormap:
                ax.set_prop_cycle(custom_cycler)
            for lab in old_labels:
                for ident in x_comb.identifiers('process'):
                    if (lab == ident.label):
                        the_order.insert(0,ident)
                        break
            if (all_bkg+all_sig>0.): hist.plot1d(x_comb,
                    overlay='process',ax=ax,
                    clear=False,
                    stack=True,
                    order=the_order,
                    fill_opts=fill_opts,
                    error_opts=err_opts,
                    overflow=overflow
                    )

        x_data = x['data']
        all_data = 0
        for key,val in x_data.values().items():
            all_data +=val.sum()

        x_allbkg = x_nosig.sum("process",overflow=overflow)
        x_nobkg.scale({p: 1./(lumifb*float(sigscale)) for p in x_nobkg.identifiers('process')}, axis="process")
        x_allsig = x_nobkg.sum("process",overflow=overflow)

        if not blind:
            blindi = 0
        else:
            blindi = len(blind)

        if (all_data>0 and blindi!=1):
            if blindi>1:
                bi = blindi - 2
                x_data_tmp = x_data[:,float(blind[bi]) if blind[bi]!='None' else None:float(blind[bi+1]) if blind[bi+1]!='None' else None]
                x_allbkg_blind_tmp = x_allbkg[float(blind[bi]) if blind[bi]!='None' else None:float(blind[bi+1]) if blind[bi+1]!='None' else None]
                while bi > 1:
                    bi = bi - 2
                    x_data_tmp += x_data[:,float(blind[bi]) if blind[bi]!='None' else None:float(blind[bi+1]) if blind[bi+1]!='None' else None]
                    x_allbkg_blind_tmp += x_allbkg[float(blind[bi]) if blind[bi]!='None' else None:float(blind[bi+1]) if blind[bi+1]!='None' else None]
                x_data = x_data_tmp
                x_allbkg_blind = x_allbkg_blind_tmp
            else:
                x_allbkg_blind = x_allbkg
            hist.plot1d(x_data,ax=ax,
                    overlay='process',
                    clear=False,
                    #line_opts=line_opts,
                    error_opts=err_opts_data,
                    overflow=overflow)

            x_data = x_data.sum("process",overflow=overflow)
            x_data.label = 'Data/MC'
            if (not noratio):
                hist.plotratio(x_data,x_allbkg_blind,ax=axr,
                    unc='num',
                    clear=False,
                    error_opts=err_opts_data,
                    denom_fill_opts=err_opts_denom,
                    guide_opts={'linestyle':'dashed','linewidth':1.5},
                    overflow=overflow)
                axr.set_ylim(0. if axr.get_ylim()[0] < 0. else None, 2. if axr.get_ylim()[1] > 2. else None)
        if (dosigrat):
            x_intsosqrtb = x_allsig.copy(content=False)
            sosqrtb = {}
            sosqrtb[var_name] = np.array(x_allsig.axis(var_name).centers())
            sosqrtb["weight"] = np.clip(np.nan_to_num(np.divide(np.array([sum(x_allsig.values()[()][i:]) for i in range(len(x_allsig.values()[()]))]),np.sqrt([sum(x_allbkg.values()[()][i:]) for i in range(len(x_allbkg.values()[()]))]))),-100.,100.)
            x_intsosqrtb.fill(**sosqrtb)
            x_intsosqrtb.label = r'$S/\sqrt{B}$'
            hist.plot1d(x_intsosqrtb,ax=axs,
                clear=False,
                line_opts=line_opts,
                overflow=overflow)
            axs.set_ylim(0. if axs.get_ylim()[0] < 0. else None, None)
            axs.get_legend().set_visible(False)
    
    
        if verbose: print('MC: %.4f Sig: %.4f - Data: %.4f'%(all_bkg,all_sig,all_data if blindi==0 else 0.))
        if verbose: print('MC: %.4f Sig: %.4f S/sqrt(B): %.4f - Data: %.4f'%(all_bkg,all_sig,all_sig/math.sqrt(all_bkg),all_data if blindi==0 else 0.))

    ax.autoscale(axis='x', tight=True)
    if len(xlimits)==2:
        try:
            ax.set_xlim(float(xlimits[0]), None)
        except:
            pass
        try:
            ax.set_xlim(None, float(xlimits[1]))
        except:
            pass
    if xlog:
        ax.set_xscale('function', functions=(lambda x: np.log10(x), lambda x: 10.**(x)))
    if xexp:
        coeff = 10.
        if (ax.get_xlim()[1]-ax.get_xlim()[0])<0.05:
            coeff = 100.
        ax.set_xscale('function', functions=(lambda x: 10.**(coeff*x), lambda x: np.log10(x/coeff)))
    ax.set_ylim(0, None)
    if norm is not None:
        ax.yaxis.set_label_text('arb.')
    if (not noratio and dosigrat): 
        ax.xaxis.set_label_text('')
        axr.xaxis.set_label_text('')
        axs.xaxis.set_label_text(var_label)
    elif (not noratio and not dosigrat):
        ax.xaxis.set_label_text('')
        axr.xaxis.set_label_text(var_label)
    elif (noratio and dosigrat):
        ax.xaxis.set_label_text('')
        axs.xaxis.set_label_text(var_label)
    else:
        ax.xaxis.set_label_text(var_label)
    ax.ticklabel_format(axis='x', style='sci')
    old_handles, old_labels = ax.get_legend_handles_labels()
    new_labels = []
    for xl in old_labels:
        #if (('H(125)' in xl or '\phi' in xl) and sigscale!=1): xl = xl + " (x " + str(sigscale) + ")"
        if (('\phi' in xl) and sigscale!=1): xl = xl + " (x " + str(sigscale) + ")"
        new_labels.append(xl)
    if not docomp:
        leg = ax.legend(handles=old_handles,labels=new_labels,title=r'%s'%plottitle,frameon=True,framealpha=1.0,facecolor='white',ncol=(2 if len(x.identifiers('process')) > 4 else 1))
    else:
        leg = ax.legend(handles=old_handles,labels=new_labels,title=r'%s'%plottitle,frameon=True,framealpha=1.0,facecolor='white',ncol=(2 if len(x.identifiers('compsel')) > 4 else 1))
    ax.set_ylim(ax.get_ylim()[0],ax.get_ylim()[1]*1.3)
    lumi = ax.text(1., 1., r"%.1f fb$^{-1}$ (13 TeV)"%lumifb,fontsize=16,horizontalalignment='right',verticalalignment='bottom',transform=ax.transAxes)
    cmstext = ax.text(0., 1., "CMS",fontsize=19,horizontalalignment='left',verticalalignment='bottom',transform=ax.transAxes, fontweight='bold')
    issim = True
    if not (solo or docomp):
        if (all_data > 0 and not blind): 
            issim = False
    elif docomp and sample=='data':
        issim = False
    addtext = ax.text(0.085, 1., "%sPreliminary"%("Simulation " if issim else ""),fontsize=16,horizontalalignment='left',verticalalignment='bottom',transform=ax.transAxes, style='italic')
    #hep.cms.cmslabel(ax, data=False, paper=False, year='2017')
    fig.savefig("%s_%s_%s_%s_lumi%i%s.pdf"%(('solo' if solo else 'comp' if docomp else 'stack'),sel,var_name,savename,lumifb,('_xexp' if xexp else '_xlog' if xlog else '')))
    if verbose: print("%s_%s_%s_%s_lumi%i%s.pdf"%(('solo' if solo else 'comp' if docomp else 'stack'),sel,var_name,savename,lumifb,('_xexp' if xexp else '_xlog' if xlog else '')))
    ax.semilogy()
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
        logmin = 0.1**float(len(decsplit[1])-len(decsplit[1].lstrip('0'))+1)
    else:
        logmin = 10.**float(len(decsplit[0])+0)
    if (docomp or solo) and density:
        ax.set_ylim(ax.get_ylim()[1]/10.,ax.get_ylim()[1]*80.)
    else:
        ax.set_ylim(logmin/10. if logmin>1. else 0.1, ax.get_ylim()[1]*80.)
    fig.savefig("%s_%s_%s_%s_lumi%i%s_logy.pdf"%(('solo' if solo else 'comp' if docomp else 'stack'),sel,var_name,savename,lumifb,('_xexp' if xexp else '_xlog' if xlog else '')))
    if verbose: print("%s_%s_%s_%s_lumi%i%s_logy.pdf"%(('solo' if solo else 'comp' if docomp else 'stack'),sel,var_name,savename,lumifb,('_xexp' if xexp else '_xlog' if xlog else '')))

    plt.close('all') 

def getPlots(args,returnHist=False):
    if args.verbose: print(args.lumi)
    lumifb = float(args.lumi)
    tag = args.tag
    savename = args.savetag
    hist_name = args.hist

    odir = 'plots/%s/'%tag
    os.system('mkdir -p %s'%odir)
    pwd = os.getcwd()

    canDoMultijet = len(args.qcd)==3
    qcdRegions = []
    if canDoMultijet:
      qcdRegions = args.qcd
    print(qcdRegions)

    hists_mapped = {}
    if args.dirname is not None:
        args.hists = [h[:-7] for h in glob.glob('%s/*.coffea'%args.dirname)]
    for h in args.hists:
        # open hists
        
        with open("%s.hist"%h, 'rb') as f:
            hists_unmapped = pickle.load(f)
        # map to hists
        for key, val in hists_unmapped.items():
            if isinstance(val, hist.Hist) and key==hist_name:
                #if val.integrate('region',[r]).sum(*[ax for ax in val.axes() if ax.name!='region'],overflow=overflow_sum).values()[()]==0.:
                #    continue
                if key in hists_mapped:
                    hists_mapped[key] = hists_mapped[key].add(processmap.apply(val).integrate('region',args.regions))
                    if canDoMultijet:
                        for r in qcdRegions:
                            hists_mapped[r] = hists_mapped[r].add(processmap.apply(val).integrate('region',[r]))
                else:
                    hists_mapped[key] = processmap.apply(val).integrate('region',args.regions)
                    if canDoMultijet:
                        for r in qcdRegions:
                            hists_mapped[r] = processmap.apply(val).integrate('region',[r])

    hist_mapped = hists_mapped[hist_name]
    qcd_hists = [hists_mapped[r] for r in qcdRegions]
    del hists_mapped

    os.chdir(odir)
    # properties
    var_name = args.var
    var_label = r"%s"%args.varlabel
    #print(args.sel)
    def makeSel(sellist, base_dict = {}):
      if (len(sellist)%3==0):
        for vi in range(int(len(sellist)/3)):
          if (sellist[vi*3+1]=='neginf'):
            base_dict[sellist[vi*3]] = [None, float(sellist[vi*3+2])]
          elif (sellist[vi*3+2]=='inf'):
            base_dict[sellist[vi*3]] = [float(sellist[vi*3+1]), None]
          else:
            base_dict[sellist[vi*3]] = [float(sellist[vi*3+1]), float(sellist[vi*3+2])]
      return base_dict
    vars_cut = makeSel(args.sel, {})
#    if (len(args.sel)%3==0):
#      for vi in range(int(len(args.sel)/3)):
#        if (args.sel[vi*3+1]=='neginf'):
#          vars_cut[args.sel[vi*3]] = [None, float(args.sel[vi*3+2])]
#        elif (args.sel[vi*3+2]=='inf'):
#          vars_cut[args.sel[vi*3]] = [float(args.sel[vi*3+1]), None]
#        else:
#          vars_cut[args.sel[vi*3]] = [float(args.sel[vi*3+1]), float(args.sel[vi*3+2])]
    if args.verbose: print(vars_cut)
    sig_cut = makeSel(args.sigsel, {})
#    if (len(args.sigsel)%3==0):
#      for vi in range(int(len(args.sigsel)/3)):
#        if (args.sigsel[vi*3+1]=='neginf'):
#          sig_cut[args.sigsel[vi*3]] = [None, float(args.sigsel[vi*3+2])]
#        elif (args.sigsel[vi*3+2]=='inf'):
#          sig_cut[args.sigsel[vi*3]] = [float(args.sigsel[vi*3+1]), None]
#        else:
#          sig_cut[args.sigsel[vi*3]] = [float(args.sigsel[vi*3+1]), float(args.sigsel[vi*3+2])]
#    if args.verbose: print(sig_cut)
    qcd_cut_shape = makeSel(args.qcdshapesel, vars_cut.copy())
    qcd_cut_num = makeSel(args.qcdnumsel, vars_cut.copy())
    qcd_cut_den = makeSel(args.qcddensel, vars_cut.copy())
    comp_var = []
    comp_i0 = []
    for iv,v in enumerate(args.compsel):
        try:
            tmp = float(v)
        except:
            if (v!='inf' and v!='neginf'):
                 comp_var.append(v)
                 comp_i0.append(iv)
    comp_i1 = comp_i0[1:]
    comp_i1.append(len(args.compsel))
    comp_cut = []
    comp_labels = []
    label_track = 0
    for ii in range(len(comp_var)):
        comp_cut.append([[float(args.compsel[vi]) if args.compsel[vi]!='neginf' else None, float(args.compsel[vi+1]) if args.compsel[vi+1]!='inf' else None] for vi in range(comp_i0[ii]+1,comp_i1[ii]-1,2)])
        nlabels = int((comp_i1[ii]-comp_i0[ii]-1)/2)
        comp_labels.append(args.complabels[label_track:label_track+nlabels])
        label_track = label_track + nlabels
    if args.verbose: print(comp_labels)
    if args.verbose: print(comp_var,comp_cut)
    if args.verbose: print(hist_mapped)

    normed = args.norm
    if normed is not None:
        normed = float(normed)
    if args.verbose: print('norm')
    if args.verbose: print(normed)

    if args.verbose: print('rebin',args.rebin)

    if not returnHist:
        drawStack(hist_mapped,args.hist,var_name,var_label,args.title,args.sample,lumifb,vars_cut,sig_cut,'',savename,args.xlimits,args.blind,args.solo,args.sigscale,args.sigstack,args.noratio,args.dosigrat,args.rebin,args.overflow,args.xlog,args.xexp,normed,args.density,args.comp,comp_var,comp_cut,comp_labels,args.defcolors,args.verbose,args.signame,systematic=args.systematic,qcd_hists=qcd_hists,qcd_cut_shape=qcd_cut_shape,qcd_cut_num=qcd_cut_num,qcd_cut_den=qcd_cut_den)
        os.chdir(pwd)
        return None

    else:
        os.chdir(pwd)
        if canDoMultijet:
            return hist_mapped,qcd_hists
        else:
            return hist_mapped,None

if __name__ == "__main__":

    if not sys.warnoptions:
        import warnings
        warnings.simplefilter("ignore")
        os.environ["PYTHONWARNINGS"] = "ignore" # Also affect subprocesses

    #ex. python plot_stack.py --hists ../condor/Oct01/hists_sum --tag Oct01 --var jet_msd --varlabel '$m_{SD}(jet)$' --title '$\mu\tau_{h},~300<p_{T}(j)<350$' --lumi 36.7 --regions hadmu_signal --hist mass_kin --savetag hadmu_jet_pt_300_350 --sel jet_pt neginf 350. --sigsel genhtt 2.5 3.5 --dosigrat
    parser = argparse.ArgumentParser()
    parser.add_argument('--hists',       dest='hists',       default="hists",      help="hists pickle name", nargs='+')
    parser.add_argument('--dirname',     dest='dirname',     default=None,         help="hists dir")
    parser.add_argument('--tag',         dest='tag',         default="",           help="tag")
    parser.add_argument('--savetag',     dest='savetag',     default="",           help="savetag")
    parser.add_argument('--var',         dest='var',         default="",           help="var")
    parser.add_argument('--varlabel',    dest='varlabel',    default="",           help="varlabel")
    parser.add_argument('--title',       dest='title',       default="",           help="title")
    parser.add_argument('--lumi',        dest='lumi',        default=50.,          help="lumi",       type=float)
    parser.add_argument('--sel',         dest='sel',         default='',           help='selection',  nargs='+')
    parser.add_argument('--sigsel',      dest='sigsel',      default='',           help='signal selection',  nargs='+')
    parser.add_argument('--signame',     dest='signame',     default=[],           help='signal names',  nargs='+')
    parser.add_argument('--regions',     dest='regions',     default='',           help='regionsel',  nargs='+')
    parser.add_argument('--hist',        dest='hist',        default='',           help='histname')
    parser.add_argument('--xlimits',     dest='xlimits',     default='',           help='xlimits',    nargs='+')
    parser.add_argument('--xlog',        dest='xlog',        action='store_true',  help='xlog')
    parser.add_argument('--xexp',        dest='xexp',        action='store_true',  help='xexp')
    parser.add_argument('--blind',       dest='blind',       default='',           help='blind',      nargs='+')
    parser.add_argument('--solo',        dest='solo',        action='store_true',  help='solo')
    parser.add_argument('--sigscale',    dest='sigscale',    default=50,           help='sigscale',   type=int)
    parser.add_argument('--sigstack',    dest='sigstack',    action='store_true',  help='sigstack')
    parser.add_argument('--noratio',     dest='noratio',     action='store_true',  help='noratio')
    parser.add_argument('--dosigrat',    dest='dosigrat',    action='store_true',  help='dosigrat')
    parser.add_argument('--overflow',    dest='overflow',    default='none',       help='overflow')
    parser.add_argument('--rebin',       dest='rebin',       default=[1],          help='rebin',   type=float,   nargs='+')
    parser.add_argument('--norm',        dest='norm',        default=None,         help='norm')
    parser.add_argument('--density',     dest='density',     action='store_true',  help='density')
    parser.add_argument('--sample',      dest='sample',      default='all',        help='sample')
    parser.add_argument('--systematic',  dest='systematic',  default='nominal',    help='systematic')
    parser.add_argument('--comp',        dest='comp',        action='store_true',  help='comp')
    parser.add_argument('--compsel',     dest='compsel',     default='',           help='compsel',     nargs='+')
    parser.add_argument('--complabels',  dest='complabels',  default='',           help='complabels',  nargs='+')
    parser.add_argument('--defcolors',   dest='defcolors',   action='store_false', help='defcolors')
    parser.add_argument('--verbose',     dest='verbose',     action='store_true',  help='verbose')
    parser.add_argument('--qcd',         dest='qcd',         default=[],           help='qcd (shape,num,den)',    nargs='+')
    parser.add_argument('--qcdshapesel', dest='qcdshapesel', default='',           help='qcdshapesel',            nargs='+')
    parser.add_argument('--qcdnumsel',   dest='qcdnumsel',   default='',           help='qcdnumsel',              nargs='+')
    parser.add_argument('--qcddensel',   dest='qcddensel',   default='',           help='qcddensel',              nargs='+')
    args = parser.parse_args()

    getPlots(args)
