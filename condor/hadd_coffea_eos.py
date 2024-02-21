#!/usr/bin/python
import json
import os
import subprocess

import argparse
import pickle

from os import listdir
from os.path import isfile, join

from coffea import hist

def run_hadd(indir, eosdir, samples, invert, ignore, outname, noscale, noresub, onlyresub, chunk_size, f_resub):

    if not samples:
        onlyfiles = [f[:-4] for f in os.listdir("condor/"+indir+"/") if os.path.isfile(os.path.join("condor/"+indir+"/", f)) and f.endswith(".jdl")]
        samples = ['MET','SingleElectron','SingleMuon', 'JetHT','Tau','EGamma']
    else:
        if (invert): 
            onlyfiles = [f[:-4] for f in os.listdir("condor/"+indir+"/") if os.path.isfile(os.path.join("condor/"+indir+"/", f)) and f.endswith(".jdl") and not any([s in f for s in samples])]
        else: 
            onlyfiles = [f[:-4] for f in os.listdir("condor/"+indir+"/") if os.path.isfile(os.path.join("condor/"+indir+"/", f)) and f.endswith(".jdl") and any([s in f for s in samples])]
            samples = ['MET','SingleElectron','SingleMuon', 'JetHT','Tau','EGamma']
    
    print(len(onlyfiles))
    os.system('xrdfs root://cmseos.fnal.gov/ ls %s%s/outfiles/ | wc -l'%(eosdir,indir))
    
    missing_files = []
    for f in onlyfiles:
        try: 
            output = subprocess.check_output('xrdfs root://cmseos.fnal.gov/ ls %s%s/outfiles/%s.hist'%(eosdir,indir,f), shell=True)
        except:
            missing_files.append(f)
    for f in missing_files:
        onlyfiles.remove(f)
 
    if not missing_files or ignore:
        if missing_files:
            print('Missing files (ignoring them):')
            print(missing_files)
            if (not noresub): 
                for mf in missing_files:
                    f_resub.write('condor_submit condor/%s/%s.jdl\n'%(indir,mf))
        if onlyresub:
            return
        chunk_names = []
        for i in range(0,len(onlyfiles),chunk_size):
          print('Chunk',i)
          flist = []
          if (i+chunk_size<len(onlyfiles)):
            for fi in onlyfiles[i:i+chunk_size]:
              x = "%s/outfiles/%s.hist"%(indir,fi)
              os.system("xrdcp -f root://cmseos.fnal.gov/%s%s condor/%s/"%(eosdir,x,indir))
              with open("condor/%s/%s.hist"%(indir,fi), 'rb') as f:
                try:
                  flist.append(pickle.load(f))
                except:
                  print('Problem loading %s, skipping...'%fi)
          else:
            for fi in onlyfiles[i:]:
              x = "%s/outfiles/%s.hist"%(indir,fi)
              os.system("xrdcp -f root://cmseos.fnal.gov/%s%s condor/%s/"%(eosdir,x,indir))
              with open("condor/%s/%s.hist"%(indir,fi), 'rb') as f:
                try:
                  flist.append(pickle.load(f))
                except:
                  print('Problem loading %s, skipping...'%fi)
        
          for key in flist[0]:
            if isinstance(flist[0][key], hist.Hist):
              for fi in range(1,len(flist)):
                flist[0][key].add(flist[fi][key])
            else:
              for fi in range(1,len(flist)):
                flist[0][key] = flist[0][key] + flist[fi][key]
          
          filename = open('condor/%s/%s_%i.hist' % (indir,outname,i), 'wb')
          pickle.dump(flist[0], filename)
          filename.close()
        
          for f in flist:
            del f
          if (i+chunk_size<len(onlyfiles)):
            for fi in onlyfiles[i:i+chunk_size]:
              x = "condor/%s/%s.hist"%(indir,fi)
              os.system("rm %s"%x)
          else:
            for fi in onlyfiles[i:]:
              x = "condor/%s/%s.hist"%(indir,fi)
              os.system("rm %s"%x)
        
          chunk_names.append('condor/%s/%s_%i.hist' % (indir,outname,i))
        
        print(chunk_names)
        if len(chunk_names)==0:
            return
        
        flist = []
        for x in chunk_names:
          with open(x, 'rb') as f:
            flist.append(pickle.load(f))

        for key in flist[0]:
          if isinstance(key, hist.Hist):
            for fi in range(1,len(flist)):
              flist[0][key].add(flist[fi][key])
          else:
            for fi in range(1,len(flist)):
              flist[0][key] = flist[0][key] + flist[fi][key]
          
        xs = {}
        with open('fileset/xsecs.json', 'r') as f:
            xs = json.load(f)
            
        scale1fb = {k: xs[k] * 1000. / w for k, w in flist[0]['sumw'].items()}

        print('sum weights  ',{k: w for k, w in flist[0]['sumw'].items()})
        print('scaling using',scale1fb)

        dylist = []
        philist = []
        for s in samples:
            if s not in scale1fb: scale1fb[s] = 1.
        for da in scale1fb:
            if "DYJets" in da:
                dylist.append(da)
            elif "Spin0" in da:
                philist.append(da)
        for dy in dylist:
            scale1fb[dy+"_Zee"] = scale1fb[dy]
            scale1fb[dy+"_Zem"] = scale1fb[dy]
            scale1fb[dy+"_Zmm"] = scale1fb[dy]
            scale1fb[dy+"_Ztt"] = scale1fb[dy]
        for phi in philist:
            scale1fb[phi+"_nomatch"] = scale1fb[phi]
        
        print('noscale =',noscale)
        if not noscale:
            for key in flist[0]:
                if isinstance(flist[0][key], hist.Hist):
                    flist[0][key].scale(scale1fb, 'dataset')
                else:
                    if key=='sumw':
                        continue
                    for samp in flist[0][key]:
                        for x in flist[0][key][samp]:
                            flist[0][key][samp][x] = flist[0][key][samp][x]*scale1fb[samp]
            
          
        filename = open('condor/%s/%s.hist' % (indir,outname), 'wb')
        pickle.dump(flist[0], filename)
        filename.close()
        for i,x in enumerate(chunk_names):
          os.system("rm condor/%s/%s_%i.hist" % (indir,outname,i*chunk_size))
    
    else:
        for mf in missing_files:
            if (not noresub): 
                #os.system('condor_submit condor/%s/%s.jdl;'%(indir,mf))
                f_resub.write('condor_submit condor/%s/%s.jdl\n'%(indir,mf))
            print('File: %s.jdl flagged succesfully!'%mf)

#python hadd_coffea_eos.py Apr30
#python hadd_coffea_eos.py Sep17 --samples DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8 DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8 DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8 DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8 DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8 DYJetsToLL_M-50_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8 DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8 --outname hists_sum_zll
samp_dict = {
    '2016APV':{
        "ST":[
            'ST_s-channel_4f_leptonDecays',
            'ST_s-channel_4f_hadronicDecays',
            'ST_t-channel_antitop_4f_InclusiveDecays',
            'ST_t-channel_top_4f_InclusiveDecays',
            'ST_tW_antitop_5f_inclusiveDecays',
            'ST_tW_top_5f_inclusiveDecays',
        ],
        "VJetsToQQ":[
            'WJetsToQQ_HT-400to600',
            'WJetsToQQ_HT-600to800',
            'WJetsToQQ_HT-800toInf',
            'ZJetsToQQ_HT-400to600',
            'ZJetsToQQ_HT-800toInf',
            'ZJetsToQQ_HT-600to800',
        ],
        "VV":[
            'WW',
            'WZ',
            'ZZ',
        ],
        "WJetsToLNu":[
            'WJetsToLNu_HT-100To200',
            'WJetsToLNu_HT-1200To2500',
            'WJetsToLNu_HT-800To1200',
            'WJetsToLNu_HT-600To800',
            'WJetsToLNu_HT-2500ToInf',
            'WJetsToLNu_HT-200To400',
            'WJetsToLNu_HT-70To100',
            'WJetsToLNu_HT-400To600',
        ],
        "QCD":[
            'QCD_Pt_3200toInf',
            'QCD_Pt_1400to1800',
            'QCD_Pt_170to300',
            'QCD_Pt_300to470',
            'QCD_Pt_470to600',
            'QCD_Pt_1000to1400',
            'QCD_Pt_2400to3200',
            'QCD_Pt_600to800',
            'QCD_Pt_800to1000',
            'QCD_Pt_1800to2400',
        ],
        "TT":[
            'TTToSemiLeptonic',
            'TTToHadronic',
            'TTTo2L2Nu',
        ],
        "DYJetsToLL":[
            'DYJetsToLL_Pt-50To100',
            'DYJetsToLL_Pt-100To250',
            'DYJetsToLL_Pt-250To400',
            'DYJetsToLL_Pt-400To650',
            'DYJetsToLL_Pt-650ToInf',
        ],
        "HTauTau":[
            'GluGluHToTauTau',
            'VBFHToTauTau',
            'WminusHToTauTau',
            'WplusHToTauTau',
            'ZHToTauTau',
            'ttHToTauTau',
        ],
        "MET":[
            'MET_Run2016B_ver2_HIPM',
            'MET_Run2016C_HIPM',
            'MET_Run2016D_HIPM',
            'MET_Run2016E_HIPM',
            'MET_Run2016F_HIPM',
        ],
        "SingleElectron":[
            'SingleElectron_Run2016B_ver2_HIPM',
            'SingleElectron_Run2016C_HIPM',
            'SingleElectron_Run2016D_HIPM',
            'SingleElectron_Run2016E_HIPM',
            'SingleElectron_Run2016F_HIPM',
        ],
        "SingleMuon":[
            'SingleMuon_Run2016B_ver2_HIPM',
            'SingleMuon_Run2016C_HIPM',
            'SingleMuon_Run2016D_HIPM',
            'SingleMuon_Run2016E_HIPM',
            'SingleMuon_Run2016F_HIPM',
        ],
        "Spin0ToTauTau":[
            'Spin0ToTauTau_2j_scalar_g1_HT300_M10_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M20_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M30_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M40_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M50_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M75_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M100_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M125_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M150_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M200_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M250_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M300_nodmx_v0_TuneCP5_MLM',
        ],
    },
    '2016':{
        "ST":[
            'ST_s-channel_4f_leptonDecays',
            'ST_s-channel_4f_hadronicDecays',
            'ST_t-channel_antitop_4f_InclusiveDecays',
            'ST_t-channel_top_4f_InclusiveDecays',
            'ST_tW_antitop_5f_inclusiveDecays',
            'ST_tW_top_5f_inclusiveDecays',
        ],
        "VJetsToQQ":[
            'WJetsToQQ_HT-400to600',
            'WJetsToQQ_HT-600to800',
            'WJetsToQQ_HT-800toInf',
            'ZJetsToQQ_HT-400to600',
            'ZJetsToQQ_HT-800toInf',
            'ZJetsToQQ_HT-600to800',
        ],
        "VV":[
            'WW',
            'WZ',
            'ZZ',
        ],
        "WJetsToLNu":[
            'WJetsToLNu_HT-100To200',
            'WJetsToLNu_HT-1200To2500',
            'WJetsToLNu_HT-800To1200',
            'WJetsToLNu_HT-600To800',
            'WJetsToLNu_HT-2500ToInf',
            'WJetsToLNu_HT-200To400',
            'WJetsToLNu_HT-70To100',
            'WJetsToLNu_HT-400To600',
        ],
        "QCD":[
            'QCD_Pt_3200toInf',
            'QCD_Pt_1400to1800',
            'QCD_Pt_170to300',
            'QCD_Pt_300to470',
            'QCD_Pt_470to600',
            'QCD_Pt_1000to1400',
            'QCD_Pt_2400to3200',
            'QCD_Pt_600to800',
            'QCD_Pt_800to1000',
            'QCD_Pt_1800to2400',
        ],
        "TT":[
            'TTToSemiLeptonic',
            'TTToHadronic',
            'TTTo2L2Nu',
        ],
        "DYJetsToLL":[
            'DYJetsToLL_Pt-50To100',
            'DYJetsToLL_Pt-100To250',
            'DYJetsToLL_Pt-250To400',
            'DYJetsToLL_Pt-400To650',
            'DYJetsToLL_Pt-650ToInf',
        ],
        "HTauTau":[
            'GluGluHToTauTau',
            'VBFHToTauTau',
            'WminusHToTauTau',
            'WplusHToTauTau',
            'ZHToTauTau',
            'ttHToTauTau',
        ],
        "MET":[
            'MET_Run2016F',
            'MET_Run2016G',
            'MET_Run2016H',
        ],
        "SingleElectron":[
            'SingleElectron_Run2016F',
            'SingleElectron_Run2016G',
            'SingleElectron_Run2016H',
        ],
        "SingleMuon":[
            'SingleMuon_Run2016F',
            'SingleMuon_Run2016G',
            'SingleMuon_Run2016H',
        ],
        "Spin0ToTauTau":[
            'Spin0ToTauTau_2j_scalar_g1_HT300_M10_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M20_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M30_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M40_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M50_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M75_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M100_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M125_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M150_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M200_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M250_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M300_nodmx_v0_TuneCP5_MLM',
        ],
    },
    '2017':{
        "DYJetsToLL":[
            'DYJetsToLL_Pt-50To100',
            'DYJetsToLL_Pt-100To250',
            'DYJetsToLL_Pt-250To400',
            'DYJetsToLL_Pt-400To650',
            'DYJetsToLL_Pt-650ToInf',
        ],
        #"GluGluHToTauTau":[
        #    'GluGluHToTauTau',
        #],
        "QCD":[
            'QCD_Pt_3200toInf',
            'QCD_Pt_1400to1800',
            'QCD_Pt_170to300',
            'QCD_Pt_300to470',
            'QCD_Pt_470to600',
            'QCD_Pt_1000to1400',
            'QCD_Pt_2400to3200',
            'QCD_Pt_600to800',
            'QCD_Pt_800to1000',
            'QCD_Pt_1800to2400',
        ],
        "ST":[
            'ST_s-channel_4f_leptonDecays',
            'ST_s-channel_4f_hadronicDecays',
            'ST_t-channel_antitop_4f_InclusiveDecays',
            'ST_t-channel_top_4f_InclusiveDecays',
            'ST_tW_antitop_5f_inclusiveDecays',
            'ST_tW_top_5f_inclusiveDecays',
        ],
        "TT":[
            'TTToSemiLeptonic',
            'TTToHadronic',
            'TTTo2L2Nu',
        ],
        "WJetsToLNu":[
            'WJetsToLNu_HT-100To200',
            'WJetsToLNu_HT-1200To2500',
            'WJetsToLNu_HT-800To1200',
            'WJetsToLNu_HT-600To800',
            'WJetsToLNu_HT-2500ToInf',
            'WJetsToLNu_HT-200To400',
            'WJetsToLNu_HT-70To100',
            'WJetsToLNu_HT-400To600',
        ],
        "VJetsToQQ":[
            'ZJetsToQQ_HT-400to600',
            'ZJetsToQQ_HT-800toInf',
            'ZJetsToQQ_HT-600to800',
            'WJetsToQQ_HT-600to800',
            'WJetsToQQ_HT-800toInf',
            'WJetsToQQ_HT-400to600',
        ],
        "VV":[
            'WW',
            'WZ',
            'ZZ',
        ],
        "HTauTau":[
            "GluGluHToTauTau",
            "VBFHToTauTau",
            "WminusHToTauTau",
            "WplusHToTauTau",
            "ZHToTauTau",
            "ttHToTauTau",
        ],
        "Spin0ToTauTau":[
            'Spin0ToTauTau_2j_scalar_g1_HT300_M10_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M20_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M30_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M40_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M50_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M75_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M100_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M125_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M150_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M200_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M250_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M300_nodmx_v0_TuneCP5_MLM',
        ],
        "SingleElectron":[
            'SingleElectron_Run2017B',
            'SingleElectron_Run2017C',
            'SingleElectron_Run2017D',
            'SingleElectron_Run2017E',
            'SingleElectron_Run2017F',
        ],
        "SingleMuon":[
            'SingleMuon_Run2017B',
            'SingleMuon_Run2017C',
            'SingleMuon_Run2017D',
            'SingleMuon_Run2017E',
            'SingleMuon_Run2017F',
        ],
        "MET":[
            'MET_Run2017B',
            'MET_Run2017C',
            'MET_Run2017D',
            'MET_Run2017E',
            'MET_Run2017F',
        ],
    },
    '2018':{
        "DYJetsToLL":[
            'DYJetsToLL_Pt-50To100',
            'DYJetsToLL_Pt-100To250',
            'DYJetsToLL_Pt-250To400',
            'DYJetsToLL_Pt-400To650',
            'DYJetsToLL_Pt-650ToInf',
        ],
        "QCD":[
            'QCD_Pt_3200toInf',
            'QCD_Pt_1400to1800',
            'QCD_Pt_170to300',
            'QCD_Pt_300to470',
            'QCD_Pt_470to600',
            'QCD_Pt_1000to1400',
            'QCD_Pt_2400to3200',
            'QCD_Pt_600to800',
            'QCD_Pt_800to1000',
            'QCD_Pt_1800to2400',
        ],
        "ST":[
            'ST_s-channel_4f_leptonDecays',
            'ST_s-channel_4f_hadronicDecays',
            'ST_t-channel_antitop_4f_InclusiveDecays',
            'ST_t-channel_top_4f_InclusiveDecays',
            'ST_tW_antitop_5f_inclusiveDecays',
            'ST_tW_top_5f_inclusiveDecays',
        ],
        "TT":[
            'TTToSemiLeptonic',
            'TTToHadronic',
            'TTTo2L2Nu',
        ],
        "HTauTau":[
            "GluGluHToTauTau",
            "VBFHToTauTau",
            "WminusHToTauTau",
            "WplusHToTauTau",
            "ZHToTauTau",
            "ttHToTauTau",
        ],
        "WJetsToLNu":[
            'WJetsToLNu_HT-100To200',
            'WJetsToLNu_HT-1200To2500',
            'WJetsToLNu_HT-800To1200',
            'WJetsToLNu_HT-600To800',
            'WJetsToLNu_HT-2500ToInf',
            'WJetsToLNu_HT-200To400',
            'WJetsToLNu_HT-70To100',
            'WJetsToLNu_HT-400To600',
        ],
        "VJetsToQQ":[
            'ZJetsToQQ_HT-400to600',
            'ZJetsToQQ_HT-800toInf',
            'ZJetsToQQ_HT-600to800',
            'WJetsToQQ_HT-600to800',
            'WJetsToQQ_HT-800toInf',
            'WJetsToQQ_HT-400to600',
        ],
        "VV":[
            'WW',
            'WZ',
            'ZZ',
        ],
        "Spin0ToTauTau":[
            'Spin0ToTauTau_2j_scalar_g1_HT300_M10_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M20_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M30_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M40_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M50_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M75_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M100_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M125_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M150_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M200_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M250_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M300_nodmx_v0_TuneCP5_MLM',
        ],
        "SingleMuon":[
            'SingleMuon_Run2018A',
            'SingleMuon_Run2018B',
            'SingleMuon_Run2018C',
            'SingleMuon_Run2018D',
        ],
        "MET":[
            'MET_Run2018A',
            'MET_Run2018B',
            'MET_Run2018C',
            'MET_Run2018D',
        ],
        "EGamma":[
            'EGamma_Run2018A',
            'EGamma_Run2018B',
            'EGamma_Run2018C',
            'EGamma_Run2018D',
        ],
    },
}

datalist = ["SingleMuon", "SingleElectron", "JetHT", "Tau", "MET", "EGamma"]

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('indir', metavar='indir', type=str, help='indir')
parser.add_argument('--eosdir', metavar='eosdir', default='/store/user/eamoreno/boostedhiggs/', type=str, help='eosdir')
parser.add_argument('--samples', metavar='samples', help='samples', nargs='+')
parser.add_argument('-i', '--invert', action='store_true')
parser.add_argument('--ignore', action='store_true')
parser.add_argument('--outname', metavar='outname', default='hists_sum', help='outname', type=str)
parser.add_argument('-n', '--noscale', action='store_true')
parser.add_argument('--noresub', action='store_true')
parser.add_argument('--onlyresub', action='store_true')
parser.add_argument('--sampsplit', action='store_true')
parser.add_argument('--year', metavar='year', default='2017', help='year')
parser.add_argument('--fullsplit', action='store_true')
parser.add_argument('--chunk', metavar='chunk', help='chunk size', type=int, default=10)
args = parser.parse_args()

#python condor/hadd_coffea_eos.py Jun06_2017_UL/ --samples Data MC --sampsplit --fullsplit --year 2017 --ignore

f_resub = open("condor/%s/condor_resub"%(args.indir), "w")

if not args.sampsplit:
    run_hadd(args.indir, args.eosdir, args.samples, args.invert, args.ignore, args.outname, args.noscale, args.noresub, args.onlyresub, args.chunk, f_resub)

else:
    if not args.samples:
        theblocks = [k for k in samp_dict[args.year]]
    
    else:
        samplist = {
            "Data":[s for s in samp_dict[args.year] if s in datalist],
            "MC":[s for s in samp_dict[args.year] if s not in datalist],
        }
        exp_samps = []
        for s in args.samples:
            if (s!="Data" and s!="MC"):
                exp_samps.append(s)
            else:
                exp_samps = exp_samps + samplist[s]
        print(exp_samps)
        theblocks = exp_samps
    
    for block in theblocks:
        if not args.fullsplit:
            run_hadd(args.indir, args.eosdir, samp_dict[args.year][block], args.invert, args.ignore, "%s_%s"%(args.outname,block), True if block in datalist else args.noscale, args.noresub, args.onlyresub, args.chunk, f_resub)
        else:
            for bs in samp_dict[args.year][block]:
                run_hadd(args.indir, args.eosdir, [bs], args.invert, args.ignore, "%s_%s"%(args.outname,bs), True if block in datalist else args.noscale, args.noresub, args.onlyresub, args.chunk, f_resub)

f_resub.close()
