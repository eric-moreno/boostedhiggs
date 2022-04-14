#!/usr/bin/python
from coffea import hist
import json
import os
import subprocess

import argparse
import pickle

from os import listdir
from os.path import isfile, join

def run_hadd(indir, eosdir, samples, invert, ignore, outname, noscale, noresub, chunk_size):

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
        if not not missing_files:
            print('Missing files (ignoring them):')
            print(missing_files)
        chunk_names = []
        for i in range(0,len(onlyfiles),chunk_size):
          print('Chunk',i)
          flist = []
          if (i+chunk_size<len(onlyfiles)):
            for fi in onlyfiles[i:i+chunk_size]:
              x = "%s/outfiles/%s.hist"%(indir,fi)
              os.system("xrdcp -f root://cmseos.fnal.gov/%s%s condor/%s/"%(eosdir,x,indir))
              with open("condor/%s/%s.hist"%(indir,fi), 'rb') as f:
                flist.append(pickle.load(f))
          else:
            for fi in onlyfiles[i:]:
              x = "%s/outfiles/%s.hist"%(indir,fi)
              os.system("xrdcp -f root://cmseos.fnal.gov/%s%s condor/%s/"%(eosdir,x,indir))
              with open("condor/%s/%s.hist"%(indir,fi), 'rb') as f:
                flist.append(pickle.load(f))
        
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
        for s in samples:
            if s not in scale1fb: scale1fb[s] = 1.
        for da in scale1fb:
            if "DYJets" in da:
                dylist.append(da)
        for dy in dylist:
            scale1fb[dy+"_Zee"] = scale1fb[dy]
            scale1fb[dy+"_Zem"] = scale1fb[dy]
            scale1fb[dy+"_Zmm"] = scale1fb[dy]
            scale1fb[dy+"_Ztt"] = scale1fb[dy]
        
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
            #if (not noresub): os.system('condor_submit condor/%s/%s.jdl;'%(indir,mf))
            print('File: %s.jdl resubmitted succesfully!'%mf)

#python hadd_coffea_eos.py Apr30
#python hadd_coffea_eos.py Sep17 --samples DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8 DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8 DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8 DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8 DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8 DYJetsToLL_M-50_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8 DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8 --outname hists_sum_zll
samp_dict = {
    '2016':{
        "SingleElectron":[
            "SingleElectron",
        ],
        "SingleMuon":[
            "SingleMuon",
        ],
        "WJetsToLNu":[
            "WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8",
            "WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8",
            "WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8",
            "WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8",
            "WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8",
            "WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8",
            "WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8",
            "WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8",
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
        "GluGluHToTauTau":[
            'GluGluHToTauTau',
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
            'ST_tW_top_5f_inclusiveDecays',
            'ST_s-channel_4f_leptonDecays',
            'ST_tW_antitop_5f_inclusiveDecays',
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
        "HTauTau":[
            #"GluGluHToTauTau_M125_13TeV_powheg_pythia8",
            #"VBFHToTauTau_M125_13TeV_powheg_pythia8",
            #"WminusHToTauTau_M125_13TeV_powheg_pythia8",
            #"WplusHToTauTau_M125_13TeV_powheg_pythia8",
            #"ZHToTauTau_M125_13TeV_powheg_pythia8",
            #"ggZH_HToTauTau_ZToLL_M125_13TeV_powheg_pythia8",
            #"ggZH_HToTauTau_ZToNuNu_M125_13TeV_powheg_pythia8",
            #"ggZH_HToTauTau_ZToQQ_M125_13TeV_powheg_pythia8",
            #"ttHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8",
            "GluGluHToTauTau",
            "VBFHToTauTau",
            "WminusHToTauTau",
            "WplusHToTauTau",
            "ZHToTauTau",
            "ttHToTauTau",
        ],
        "Spin0ToTauTau":[
            'Spin0ToTauTau_2j_scalar_g1_HT300_M10_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M30_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M50_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M125_nodmx_v0_TuneCP5_MLM',
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
            'DYJetsToLL_Pt-100To250_TuneCP5_13TeV-amcatnloFXFX-pythia8',
            'DYJetsToLL_Pt-250To400_TuneCP5_13TeV-amcatnloFXFX-pythia8',
            'DYJetsToLL_Pt-400To650_TuneCP5_13TeV-amcatnloFXFX-pythia8',
            'DYJetsToLL_Pt-650ToInf_TuneCP5_13TeV-amcatnloFXFX-pythia8',
        ],
        #"GluGluHToTauTau":[
        #    'GluGluHTauTau_13TeV',
        #],
        "QCD":[
            'QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8',
            'QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8',
            'QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8',
            'QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8',
            'QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8',
        ],
        "ST":[
            'ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-madgraph-pythia8',
            'ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-madgraph-pythia8',
            'ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8',
            'ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8',
            'ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8',
            'ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8',
        ],
        "TT":[
            'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8',
            'TTToHadronic_TuneCP5_13TeV-powheg-pythia8',
            'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8',
        ],
        "HTauTau":[
            'GluGluHToTauTau_M125_13TeV_powheg_pythia8',
            'VBFHToTauTau_M125_13TeV_powheg_pythia8',
            'WminusHToTauTau_M125_13TeV_powheg_pythia8',
            'WplusHToTauTau_M125_13TeV_powheg_pythia8',
            'ZHToTauTau_M125_13TeV_powheg_pythia8',
            'ggZH_HToTauTau_ZToLL_M125_13TeV_powheg_pythia8',
            'ggZH_HToTauTau_ZToNuNu_M125_13TeV_powheg_pythia8',
            'ggZH_HToTauTau_ZToQQ_M125_13TeV_powheg_pythia8',
            'ttHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8',
        ],
        "WJetsToLNu":[
            'WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8',
        ],
        "VJetsToQQ":[
            'WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8',
            'WJetsToQQ_HT400to600_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8',
            'WJetsToQQ_HT600to800_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8',
            'ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8',
            'ZJetsToQQ_HT400to600_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8',
            'ZJetsToQQ_HT600to800_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8',
        ],
        "SingleMuon":[
            #'SingleMuon_pancakes-02_Run2018A-12Nov2019_UL2018_rsb-v1',
            #'SingleMuon_pancakes-02_Run2018B-12Nov2019_UL2018-v2',
            #'SingleMuon_pancakes-02_Run2018C-12Nov2019_UL2018-v2',
            #'SingleMuon_pancakes-02_Run2018D-12Nov2019_UL2018-v4',
            "SingleMuon",
        ],
        "MET":[
            'MET_pancakes-02_Run2018A-12Nov2019_UL2018-v3',
            'MET_pancakes-02_Run2018B-12Nov2019_UL2018_rsb-v1',
            'MET_pancakes-02_Run2018C-12Nov2019_UL2018_rsb-v1',
            'MET_pancakes-02_Run2018D-12Nov2019_UL2018_rsb-v2',
        ],
        "EGamma":[
            #'EGamma_pancakes-02_Run2018A-12Nov2019_UL2018-v2',
            #'EGamma_pancakes-02_Run2018B-12Nov2019_UL2018-v2',
            #'EGamma_pancakes-02_Run2018C-12Nov2019_UL2018-v2',
            #'EGamma_pancakes-02_Run2018D-12Nov2019_UL2018-v4',
            "EGamma",
        ],
    },
}

datalist = ["SingleMuon", "SingleElectron", "JetHT", "Tau", "MET", "EGamma"]

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('indir', metavar='indir', type=str, help='indir')
parser.add_argument('--eosdir', metavar='eosdir', default='/store/user/drankin/boostedhiggs/', type=str, help='eosdir')
parser.add_argument('--samples', metavar='samples', help='samples', nargs='+')
parser.add_argument('-i', '--invert', action='store_true')
parser.add_argument('--ignore', action='store_true')
parser.add_argument('--outname', metavar='outname', default='hists_sum', help='outname', type=str)
parser.add_argument('-n', '--noscale', action='store_true')
parser.add_argument('--noresub', action='store_true')
parser.add_argument('--sampsplit', action='store_true')
parser.add_argument('--year', metavar='year', default='2017', help='year')
parser.add_argument('--fullsplit', action='store_true')
parser.add_argument('--chunk', metavar='chunk', help='chunk size', type=int, default=10)
args = parser.parse_args()

#python condor/hadd_coffea_eos.py Apr10_2017_UL/ --samples Data MC --sampsplit --fullsplit --year 2017 --ignore

if not args.sampsplit:
    run_hadd(args.indir, args.eosdir, args.samples, args.invert, args.ignore, args.outname, args.noscale, args.noresub, args.chunk)

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
            run_hadd(args.indir, args.eosdir, samp_dict[args.year][block], args.invert, args.ignore, "%s_%s"%(args.outname,block), True if block in datalist else args.noscale, args.noresub, args.chunk)
        else:
            for bs in samp_dict[args.year][block]:
                run_hadd(args.indir, args.eosdir, [bs], args.invert, args.ignore, "%s_%s"%(args.outname,bs), True if block in datalist else args.noscale, args.noresub, args.chunk)

