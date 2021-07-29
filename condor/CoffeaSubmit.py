#!/usr/bin/python

import argparse
import os
import sys
import re
import fileinput

import json
import glob

#python CoffeaSubmit.py Apr23 run_htt.py 50 1
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('settings', metavar='tag script nfiles re-tar', type=str, nargs='+',
                   help='label scriptname (re-tar)')
args = parser.parse_args()

if (not ((len(args.settings) is 3) or (len(args.settings) is 4))):
    print("Wrong number of arguments (must be 3 or 4, found",len(args.settings),")")
    parser.print_help()
    sys.exit()

label = args.settings[0]
if '_2016' in label or '_2017' in label or '_2018' in label:
    print('year listed in label arg (%s), exiting'%label)
    exit()

script = args.settings[1]
files_per_job = int(args.settings[2])

loc_base = os.environ['PWD']

samplelist = {
    '2016': [
     ],

    '2017': [
        #'WW_TuneCP5_13TeV-pythia8',
        #'WZ_TuneCP5_13TeV-pythia8',
        #'ZZ_TuneCP5_13TeV-pythia8',
        'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8',
        'TTToHadronic_TuneCP5_13TeV-powheg-pythia8',
        'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8',
        'QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8',
        #'QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8',
        'QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8',
        'QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8',
        'QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8',
        'QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8',
        #'QCD_HT50to100_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8',
        'QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8',
        'ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8',
        'ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8',
        'ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8',
        'ST_t-channel_top_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8',
        'ST_tW_antitop_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8',
        'ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8',
        #'DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8',
        #'DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8',
        #'DYJetsToLL_M-50_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8',
        #'DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8',
        #'DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8',
        #'DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8',
        #'DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8',
        #'DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8',
        'DYJetsToLL_Pt-100To250_TuneCP5_13TeV-amcatnloFXFX-pythia8',
        'DYJetsToLL_Pt-250To400_TuneCP5_13TeV-amcatnloFXFX-pythia8',
        'DYJetsToLL_Pt-400To650_TuneCP5_13TeV-amcatnloFXFX-pythia8',
        'DYJetsToLL_Pt-650ToInf_TuneCP5_13TeV-amcatnloFXFX-pythia8',
        #'WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8',
        'WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8',
        'WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8',
        'WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8',
        'WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8',
        'WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8',
        'WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8',
        'WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8',
        'WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8',
        'WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8',
        'WJetsToQQ_HT400to600_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8',
        'WJetsToQQ_HT600to800_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8',
        'ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8',
        'ZJetsToQQ_HT400to600_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8',
        'ZJetsToQQ_HT600to800_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8',
        'GluGluHToTauTau_M125_13TeV_powheg_pythia8',
        'VBFHToTauTau_M125_13TeV_powheg_pythia8',
        'WminusHToTauTau_M125_13TeV_powheg_pythia8',
        'WplusHToTauTau_M125_13TeV_powheg_pythia8',
        'ZHToTauTau_M125_13TeV_powheg_pythia8',
        'ggZH_HToTauTau_ZToLL_M125_13TeV_powheg_pythia8',
        'ggZH_HToTauTau_ZToNuNu_M125_13TeV_powheg_pythia8',
        'ggZH_HToTauTau_ZToQQ_M125_13TeV_powheg_pythia8',
        'ttHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8',

        'JetHT_pancakes-02_Run2017B-09Aug2019_UL2017-v1',
        'JetHT_pancakes-02_Run2017C-09Aug2019_UL2017-v1',
        'JetHT_pancakes-02_Run2017D-09Aug2019_UL2017-v1',
        'JetHT_pancakes-02_Run2017E-09Aug2019_UL2017-v1',
        'JetHT_pancakes-02_Run2017F-09Aug2019_UL2017-v1',
        'SingleMuon_pancakes-02-withPF_Run2017B-09Aug2019_UL2017-v1',
        'SingleMuon_pancakes-02-withPF_Run2017C-09Aug2019_UL2017-v1',
        'SingleMuon_pancakes-02-withPF_Run2017D-09Aug2019_UL2017-v1',
        'SingleMuon_pancakes-02-withPF_Run2017E-09Aug2019_UL2017-v1',
        'SingleMuon_pancakes-02-withPF_Run2017F-09Aug2019_UL2017-v1',
        'SingleElectron_pancakes-02-withPF_Run2017B-09Aug2019_UL2017-v1',
        'SingleElectron_pancakes-02-withPF_Run2017C-09Aug2019_UL2017-v1',
        'SingleElectron_pancakes-02-withPF_Run2017D-09Aug2019_UL2017-v1',
        'SingleElectron_pancakes-02-withPF_Run2017E-09Aug2019_UL2017-v1',
        'SingleElectron_pancakes-02-withPF_Run2017F-09Aug2019_UL2017_rsb-v2',
        'MET_pancakes-02-withPF_Run2017B-09Aug2019_UL2017_rsb-v1',
        'MET_pancakes-02-withPF_Run2017C-09Aug2019_UL2017_rsb-v1',
        'MET_pancakes-02-withPF_Run2017D-09Aug2019_UL2017_rsb-v1',
        'MET_pancakes-02-withPF_Run2017E-09Aug2019_UL2017_rsb-v1',
        'MET_pancakes-02-withPF_Run2017F-09Aug2019_UL2017_rsb-v1',
    ],

    '2018': [
        'DYJetsToLL_Pt-100To250_TuneCP5_13TeV-amcatnloFXFX-pythia8',
        'DYJetsToLL_Pt-250To400_TuneCP5_13TeV-amcatnloFXFX-pythia8',
        'DYJetsToLL_Pt-400To650_TuneCP5_13TeV-amcatnloFXFX-pythia8',
        'DYJetsToLL_Pt-650ToInf_TuneCP5_13TeV-amcatnloFXFX-pythia8',
        #'GluGluHTauTau_13TeV',
        'GluGluHToTauTau_M125_13TeV_powheg_pythia8',
        'QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8',
        'QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8',
        'QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8',
        'QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8',
        'QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8',
        'ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-madgraph-pythia8',
        'ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-madgraph-pythia8',
        'ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8',
        'ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8',
        'ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8',
        'ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8',
        'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8',
        'TTToHadronic_TuneCP5_13TeV-powheg-pythia8',
        'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8',
        'VBFHToTauTau_M125_13TeV_powheg_pythia8',
        'WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8',
        'WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8',
        'WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8',
        'WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8',
        'WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8',
        'WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8',
        'WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8',
        'WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8',
        'WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8',
        'WJetsToQQ_HT400to600_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8',
        'WJetsToQQ_HT600to800_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8',
        'WminusHToTauTau_M125_13TeV_powheg_pythia8',
        'WplusHToTauTau_M125_13TeV_powheg_pythia8',
        'ZHToTauTau_M125_13TeV_powheg_pythia8',
        'ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8',
        'ZJetsToQQ_HT400to600_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8',
        'ZJetsToQQ_HT600to800_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8',
        'ggZH_HToTauTau_ZToLL_M125_13TeV_powheg_pythia8',
        'ggZH_HToTauTau_ZToNuNu_M125_13TeV_powheg_pythia8',
        'ggZH_HToTauTau_ZToQQ_M125_13TeV_powheg_pythia8',
        'ttHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8',
    
        'SingleMuon_pancakes-02_Run2018A-12Nov2019_UL2018_rsb-v1',
        'SingleMuon_pancakes-02_Run2018B-12Nov2019_UL2018-v2',
        'SingleMuon_pancakes-02_Run2018C-12Nov2019_UL2018-v2',
        'SingleMuon_pancakes-02_Run2018D-12Nov2019_UL2018-v4',
        'MET_pancakes-02_Run2018A-12Nov2019_UL2018-v3',
        'MET_pancakes-02_Run2018B-12Nov2019_UL2018_rsb-v1',
        'MET_pancakes-02_Run2018C-12Nov2019_UL2018_rsb-v1',
        'MET_pancakes-02_Run2018D-12Nov2019_UL2018_rsb-v2',
        'EGamma_pancakes-02_Run2018A-12Nov2019_UL2018-v2',
        'EGamma_pancakes-02_Run2018B-12Nov2019_UL2018-v2',
        'EGamma_pancakes-02_Run2018C-12Nov2019_UL2018-v2',
        'EGamma_pancakes-02_Run2018D-12Nov2019_UL2018-v4',
    ]
}

samplelist['2017'] = []
#samplelist['2018'] = ['GluGluHToTauTau_M125_13TeV_powheg_pythia8']

#
#################################################

os.chdir('..')
os.system('xrdcp -f test/%s root://cmseos.fnal.gov//store/user/drankin/'%script)
if (len(args.settings) is 4):
    os.system('tar -zcf dazsle_coffea.tgz . --exclude="*nano*.root" --exclude="*.pdf" --exclude="*.pyc" --exclude=tmp --exclude="*.tgz" --exclude-vcs --exclude-caches-all')
    os.system('xrdcp -f dazsle_coffea.tgz root://cmseos.fnal.gov//store/user/drankin/dazsle_coffea.tgz')
os.chdir(loc_base)

for year in ['2016','2017','2018']:
    if (len(samplelist[year])==0):
        continue

    #################################################
    ### Names to give to your output root files
    #################################################
    
    prefix = label + "_" + year

    logdir = label + "_" + year
    outdir = '/store/user/drankin/coffea/'+logdir+'/'
    
    ################################################
    ### Where is your list of root files to run over
    ################################################
    print(label,' ',outdir)
    
    print(str(files_per_job)+' files per job...')
    
    #make local directory
    locdir = logdir
    os.system('mkdir -p  %s' %locdir)
    
    print('CONDOR work dir: '+outdir)
    #os.system('rm -rf '+outdir+label)
    os.system('mkdir -p /eos/uscms'+outdir)
    
    totfiles = {}
    
    if year=='2016':
        print('no 2016 files yet...')

    elif year=='2017':
        with open('../data/fileset2017.json', 'r') as f:
            newfiles = json.load(f)
            totfiles.update(newfiles)
        
        with open('../data/fileset2017UL.json', 'r') as f:
            newfiles = json.load(f)
            totfiles.update(newfiles)
    elif year=='2018':    
        with open('../data/fileset2018UL.json', 'r') as f:
            newfiles = json.load(f)
            totfiles.update(newfiles)
    else:
        print('invalid year specified: %s'%year)
    
    for sample in samplelist[year]:
        totfiles[sample] = len(totfiles[sample])
    
    nsubmit = 0
    
    for sample in samplelist[year]:
    
        prefix = sample
        print('Submitting '+prefix)
    
        njobs = int(totfiles[sample]/files_per_job)+1
        remainder = totfiles[sample]-int(files_per_job*(njobs-1))
    
        for j in range(njobs):
    
            condor_templ_file = open(loc_base+"/Coffea.templ.condor")
            sh_templ_file    = open(loc_base+"/Coffea.templ.sh")
        
            localcondor = locdir+'/'+prefix+"_"+str(j)+".condor"
            condor_file = open(localcondor,"w")
            for line in condor_templ_file:
                line=line.replace('DIRECTORY',locdir)
                line=line.replace('PREFIX',prefix)
                line=line.replace('JOBID',str(j))
                condor_file.write(line)
            condor_file.close()
        
            #copy local to eos
            #os.system('xrdcp -f %s %s' % (localcondor,eoscondor))
            #remove local copy
            #os.system('rm %s' % localcondor)
        
            localsh=locdir+'/'+prefix+"_"+str(j)+".sh"
            eosoutput="root://cmseos.fnal.gov/"+outdir+"/"+prefix+'_'+str(j)+'.coffea'
            sh_file = open(localsh,"w")
            for line in sh_templ_file:
                line=line.replace('SCRIPTNAME',script)
                line=line.replace('FILENUM',str(j))
                line=line.replace('YEAR',year)
                line=line.replace('SAMPLE',sample)
                line=line.replace('STARTNUM',str(j*files_per_job))
                line=line.replace('ENDNUM',str((j+1)*files_per_job))
                line=line.replace('EOSOUT',eosoutput)
                sh_file.write(line)
            sh_file.close()
    
            os.system('chmod u+x '+locdir+'/'+prefix+'_'+str(j)+'.sh')
            #print('condor file is: '+localcondor)
            if (os.path.exists('%s.log'  % localcondor)):
                os.system('rm %s.log' % localcondor)
            os.system('condor_submit %s' % localcondor)
    
            condor_templ_file.close()
            sh_templ_file.close()
    
            nsubmit = nsubmit + 1
    
    print(nsubmit,"jobs submitted.")
