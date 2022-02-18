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

datasets = {
    '2016':['SingleElectron','SingleMuon'],
    '2017':['SingleElectron','SingleMuon'],
    '2018':['EGamma','SingleMuon'],
}
data = {}
for year in datasets:
    data[year] = {}
    for dataset in datasets[year]:
        with open("./hists_sum_%s_%s.hist"%(dataset,year), 'rb') as f:
            data[year][dataset] = pickle.load(f)

print(data)

err_opts = {
    'had%s_jet':{
        #'linestyle':'-',
        'marker': '.',
        'markersize': 15.,
        'color':'k',
        'elinewidth': 1,
        'emarker': '-'
    },
    'had%s_lep':{
        #'linestyle':'-',
        'marker': '.',
        'markersize': 12.,
        'color':'b',
        'elinewidth': 0.8,
        'emarker': '-'
    },
    'had%s_signal':{
        #'linestyle':'-',
        'marker': '.',
        'markersize': 9.,
        'color':'r',
        'elinewidth': 0.6,
        'emarker': '-'
    },
}


for year in data:
    for dataset in data[year]:
        if dataset=='SingleMuon':
            lep = 'mu'
        else:
            lep = 'el'
        fig,ax = plt.subplots(1,1, figsize=(8,8))
        ax.set_title('%s (%s)'%(dataset,year))
        for region in err_opts:
            base = data[year][dataset]['met_nn_kin'].sum('systematic').integrate('dataset',dataset).integrate('region',region%lep)
            hist.plotratio(
                base.integrate('met_trigger',slice(0.5,1.5)),
                base.sum('met_trigger'),
                ax=ax,
                clear=False,
                error_opts=err_opts[region],
                unc='clopper-pearson',
            )
        fig.savefig("plots/%s/met_trigger_%s_%i_logy.pdf"%(args.tag,dataset,year))
