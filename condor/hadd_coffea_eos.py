#!/usr/bin/python
from coffea import util,hist
import json
import os
import subprocess

import argparse

from os import listdir
from os.path import isfile, join

def run_hadd(indir, eosdir, samples, invert, ignore, outname, noscale, noresub, chunk_size):

    if not samples:
        onlyfiles = [f[:-7] for f in os.listdir("../condor/"+indir+"/") if os.path.isfile(os.path.join("../condor/"+indir+"/", f)) and f.endswith(".condor")]
        samples = ['MET','SingleElectron','SingleMuon', 'JetHT','Tau','EGamma']
    else:
        if (invert): 
            onlyfiles = [f[:-7] for f in os.listdir("../condor/"+indir+"/") if os.path.isfile(os.path.join("../condor/"+indir+"/", f)) and f.endswith(".condor") and not any([s in f for s in samples])]
        else: 
            onlyfiles = [f[:-7] for f in os.listdir("../condor/"+indir+"/") if os.path.isfile(os.path.join("../condor/"+indir+"/", f)) and f.endswith(".condor") and any([s in f for s in samples])]
            samples = ['MET','SingleElectron','SingleMuon', 'JetHT','Tau','EGamma']
    
    print(len(onlyfiles))
    os.system('eos root://cmseos.fnal.gov/ ls %s%s/ | wc -l'%(eosdir,indir))
    
    missing_files = []
    for f in onlyfiles:
        try: 
            output = subprocess.check_output('eos root://cmseos.fnal.gov/ ls %s%s/%s.coffea'%(eosdir,indir,f), shell=True)
        except:
            missing_files.append(f)
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
              x = "%s/%s.coffea"%(indir,fi)
              os.system("xrdcp -f root://cmseos.fnal.gov/%s%s %s/"%(eosdir,x,indir))
              flist.append(util.load(x))
          else:
            for fi in onlyfiles[i:]:
              x = "%s/%s.coffea"%(indir,fi)
              os.system("xrdcp -f root://cmseos.fnal.gov/%s%s %s/"%(eosdir,x,indir))
              flist.append(util.load(x))
        
          for key in flist[0]:
            if isinstance(flist[0][key], hist.Hist):
              for fi in range(1,len(flist)):
                flist[0][key].add(flist[fi][key])
            else:
              for fi in range(1,len(flist)):
                flist[0][key] = flist[0][key] + flist[fi][key]
          
          util.save(flist[0],'%s/%s_%i.coffea' % (indir,outname,i))
        
          for f in flist:
            del f
          if (i+chunk_size<len(onlyfiles)):
            for fi in onlyfiles[i:i+chunk_size]:
              x = "%s/%s.coffea"%(indir,fi)
              os.system("rm %s"%x)
          else:
            for fi in onlyfiles[i:]:
              x = "%s/%s.coffea"%(indir,fi)
              os.system("rm %s"%x)
        
          chunk_names.append('%s/%s_%i.coffea' % (indir,outname,i))
        
        print(chunk_names)
        
        flist = [ util.load(x) for x in chunk_names ]

        for key in flist[0]:
          if isinstance(key, hist.Hist):
            for fi in range(1,len(flist)):
              flist[0][key].add(flist[fi][key])
          else:
            for fi in range(1,len(flist)):
              flist[0][key] = flist[0][key] + flist[fi][key]
          
        xs = {}
        with open('../data/xsec.json', 'r') as f:
            xs = json.load(f)
            
        scale1fb = {k: xs[k] * 1000. / w for k, w in flist[0]['sumw'].items()}

        print('sum weights  ',{k: w for k, w in flist[0]['sumw'].items()})
        print('scaling using',scale1fb)

        for s in samples:
            if s not in scale1fb: scale1fb[s] = 1.
        
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
            
          
        util.save(flist[0],'%s/%s.coffea' % (indir,outname))
        for i,x in enumerate(chunk_names):
          os.system("rm %s/%s_%i.coffea" % (indir,outname,i*chunk_size))
    
    else:
        for mf in missing_files:
            if (not noresub):
                os.system('rm %s/%s.condor.log;'%(indir,mf))
                os.system('condor_submit %s/%s.condor;'%(indir,mf))
            print('File: %s.condor resubmitted succesfully!'%mf)

#python hadd_coffea_eos.py Apr30
#python hadd_coffea_eos.py Sep17 --samples DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8 DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8 DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8 DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8 DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8 DYJetsToLL_M-50_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8 DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8 --outname hists_sum_zll
samp_dict = {
    '2017':{
        "DYJetsToLL":[
            "DYJetsToLL_Pt-100To250_TuneCP5_13TeV-amcatnloFXFX-pythia8",
            "DYJetsToLL_Pt-250To400_TuneCP5_13TeV-amcatnloFXFX-pythia8",
            "DYJetsToLL_Pt-400To650_TuneCP5_13TeV-amcatnloFXFX-pythia8",
            "DYJetsToLL_Pt-650ToInf_TuneCP5_13TeV-amcatnloFXFX-pythia8",
        ],
        "GluGluHToTauTau":[
            "GluGluHToTauTau",
        ],
        "QCD":[
            "QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8",
            "QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8",
            "QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8",
            "QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8",
            "QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8",
            "QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8",
        ],
        "ST":[
            "ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8",
            "ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8",
            "ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8",
            "ST_t-channel_top_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8",
            "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8",
            "ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8",
        ],
        "TT":[
            "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8",
            "TTToHadronic_TuneCP5_13TeV-powheg-pythia8",
            "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8",
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
        "VJetsToQQ":[
            "WJetsToQQ_HT400to600_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8",
            "WJetsToQQ_HT600to800_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8",
            "WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8",
            "ZJetsToQQ_HT400to600_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8",
            "ZJetsToQQ_HT600to800_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8",
            "ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8",
        ],
        "HTauTau":[
            "GluGluHToTauTau_M125_13TeV_powheg_pythia8",
            "VBFHToTauTau_M125_13TeV_powheg_pythia8",
            "WminusHToTauTau_M125_13TeV_powheg_pythia8",
            "WplusHToTauTau_M125_13TeV_powheg_pythia8",
            "ZHToTauTau_M125_13TeV_powheg_pythia8",
            "ggZH_HToTauTau_ZToLL_M125_13TeV_powheg_pythia8",
            "ggZH_HToTauTau_ZToNuNu_M125_13TeV_powheg_pythia8",
            "ggZH_HToTauTau_ZToQQ_M125_13TeV_powheg_pythia8",
            "ttHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8",
            #"GluGluHToTauTau",
        ],
        "SingleElectron":[
            "SingleElectron_pancakes-02-withPF_Run2017B-09Aug2019_UL2017-v1",
            "SingleElectron_pancakes-02-withPF_Run2017C-09Aug2019_UL2017-v1",
            "SingleElectron_pancakes-02-withPF_Run2017D-09Aug2019_UL2017-v1",
            "SingleElectron_pancakes-02-withPF_Run2017E-09Aug2019_UL2017-v1",
            "SingleElectron_pancakes-02-withPF_Run2017F-09Aug2019_UL2017_rsb-v2",
        ],
        "JetHT":[
            "JetHT_pancakes-02_Run2017B-09Aug2019_UL2017-v1",
            "JetHT_pancakes-02_Run2017C-09Aug2019_UL2017-v1",
            "JetHT_pancakes-02_Run2017D-09Aug2019_UL2017-v1",
            "JetHT_pancakes-02_Run2017E-09Aug2019_UL2017-v1",
            "JetHT_pancakes-02_Run2017F-09Aug2019_UL2017-v1",
        ],
        "SingleMuon":[
            "SingleMuon_pancakes-02-withPF_Run2017B-09Aug2019_UL2017-v1",
            "SingleMuon_pancakes-02-withPF_Run2017C-09Aug2019_UL2017-v1",
            "SingleMuon_pancakes-02-withPF_Run2017D-09Aug2019_UL2017-v1",
            "SingleMuon_pancakes-02-withPF_Run2017E-09Aug2019_UL2017-v1",
            "SingleMuon_pancakes-02-withPF_Run2017F-09Aug2019_UL2017-v1",
        ],
        "MET":[
            "MET_pancakes-02-withPF_Run2017B-09Aug2019_UL2017_rsb-v1",
            "MET_pancakes-02-withPF_Run2017C-09Aug2019_UL2017_rsb-v1",
            "MET_pancakes-02-withPF_Run2017D-09Aug2019_UL2017_rsb-v1",
            "MET_pancakes-02-withPF_Run2017E-09Aug2019_UL2017_rsb-v1",
            "MET_pancakes-02-withPF_Run2017F-09Aug2019_UL2017_rsb-v1",
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
            'SingleMuon_pancakes-02_Run2018A-12Nov2019_UL2018_rsb-v1',
            'SingleMuon_pancakes-02_Run2018B-12Nov2019_UL2018-v2',
            'SingleMuon_pancakes-02_Run2018C-12Nov2019_UL2018-v2',
            'SingleMuon_pancakes-02_Run2018D-12Nov2019_UL2018-v4',
        ],
        "MET":[
            'MET_pancakes-02_Run2018A-12Nov2019_UL2018-v3',
            'MET_pancakes-02_Run2018B-12Nov2019_UL2018_rsb-v1',
            'MET_pancakes-02_Run2018C-12Nov2019_UL2018_rsb-v1',
            'MET_pancakes-02_Run2018D-12Nov2019_UL2018_rsb-v2',
        ],
        "EGamma":[
            'EGamma_pancakes-02_Run2018A-12Nov2019_UL2018-v2',
            'EGamma_pancakes-02_Run2018B-12Nov2019_UL2018-v2',
            'EGamma_pancakes-02_Run2018C-12Nov2019_UL2018-v2',
            'EGamma_pancakes-02_Run2018D-12Nov2019_UL2018-v4',
        ],
    },
}

datalist = ["SingleMuon", "SingleElectron", "JetHT", "Tau", "MET", "EGamma"]

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('indir', metavar='indir', type=str, help='indir')
parser.add_argument('--eosdir', metavar='eosdir', default='/store/user/drankin/coffea/', type=str, help='eosdir')
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

