import argparse
import os
import re
import fileinput

import json
import glob
import sys

'''
 Submit condor jobs of processor
 Run as e.g.: python submit.py Feb27 run.py 10 2016
 Arguments:
  = [0]: tag of jobs and output dir in eos e.g. Jul1
  - [1]: script to run e.g. run.py (needs to be included in transfer_files in templ.jdl)
  - [2]: number of files per job e.g. 20
  - [3]: year
  - [4]: additional options for processor
'''
# Note: change username in `drankin` in this script

homedir = "/store/user/drankin/boostedhiggs/"

parser = argparse.ArgumentParser(description='Submit coffea condor jobs.')
parser.add_argument('settings', metavar='S', type=str, nargs='+', help='label scriptname files_per_job year options')
args = parser.parse_args()

if (not ((len(args.settings) == 4) or (len(args.settings) == 5))):
    print("Wrong number of arguments (must be 3 or 4, found", len(args.settings), ")")
    sys.exit()

label = args.settings[0]
script = args.settings[1]  # should be run.py
files_per_job = int(args.settings[2])
year = args.settings[3]
addoptions = ""
if len(args.settings) == 5:
    addoptions = args.settings[4]

loc_base = os.environ['PWD']

# list of samples to run and recos to use
recos = ["UL"] 

# if empty run over all the datasets in both filesets
samples = {
    "UL": {
        '2016APV':[
            'SingleMuon', 
            'SingleElectron', 
            #'WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8', 
            'WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8', 
            'WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8', 
            'WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8', 
            'WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8', 
            'WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8', 
            'WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8', 
        ],
        '2016':[
            #'SingleMuon', 
            #'SingleElectron', 
            #'WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8',
            #'WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8', 
            #'WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8', 
            #'WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8', 
            #'WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8', 
            #'WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8', 
            #'WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8', 
            #'WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8', 
            'Spin0ToTauTau_2j_scalar_g1_HT300_M50_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M100_nodmx_v0_TuneCP5_MLM',
            'Spin0ToTauTau_2j_scalar_g1_HT300_M150_nodmx_v0_TuneCP5_MLM',
        ],
        '2017':[
#            'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8',
#            'TTToHadronic_TuneCP5_13TeV-powheg-pythia8',
#            'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8',
#            'QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8',
#            'QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8',
#            'QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8',
#            'QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8',
#            'QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8',
#            'ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8',
#            'ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8',
#            'ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8',
#            'ST_t-channel_top_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8',
#            'ST_tW_antitop_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8',
#            'ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8',
#            'DYJetsToLL_Pt-100To250_TuneCP5_13TeV-amcatnloFXFX-pythia8',
#            'DYJetsToLL_Pt-250To400_TuneCP5_13TeV-amcatnloFXFX-pythia8',
#            'DYJetsToLL_Pt-400To650_TuneCP5_13TeV-amcatnloFXFX-pythia8',
#            'DYJetsToLL_Pt-650ToInf_TuneCP5_13TeV-amcatnloFXFX-pythia8',
#            'WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8',
#            'WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8',
#            'WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8',
#            'WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8',
#            'WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8',
#            'WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8',
#            'WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8',
#            'WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8',
#            'WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8',
#            'WJetsToQQ_HT400to600_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8',
#            'WJetsToQQ_HT600to800_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8',
#            'ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8',
#            'ZJetsToQQ_HT400to600_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8',
#            'ZJetsToQQ_HT600to800_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8',
#            'GluGluHToTauTau_M125_13TeV_powheg_pythia8',
#            'VBFHToTauTau_M125_13TeV_powheg_pythia8',
#            'WminusHToTauTau_M125_13TeV_powheg_pythia8',
#            'WplusHToTauTau_M125_13TeV_powheg_pythia8',
#            'ZHToTauTau_M125_13TeV_powheg_pythia8',
#            'ggZH_HToTauTau_ZToLL_M125_13TeV_powheg_pythia8',
#            'ggZH_HToTauTau_ZToNuNu_M125_13TeV_powheg_pythia8',
#            'ggZH_HToTauTau_ZToQQ_M125_13TeV_powheg_pythia8',
#            'ttHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8',
#            'Spin0ToTauTau_2j_scalar_g1_HT300_M10_nodmx_v0_TuneCP5_MLM',
#            'Spin0ToTauTau_2j_scalar_g1_HT300_M30_nodmx_v0_TuneCP5_MLM',
#            'Spin0ToTauTau_2j_scalar_g1_HT300_M50_nodmx_v0_TuneCP5_MLM',
#            'Spin0ToTauTau_2j_scalar_g1_HT300_M125_nodmx_v0_TuneCP5_MLM',
#            'Spin0ToTauTau_2j_scalar_g1_HT300_M300_nodmx_v0_TuneCP5_MLM',
#    
#            'SingleMuon_pancakes-02-withPF_Run2017B-09Aug2019_UL2017-v1',
#            'SingleMuon_pancakes-02-withPF_Run2017C-09Aug2019_UL2017-v1',
#            'SingleMuon_pancakes-02-withPF_Run2017D-09Aug2019_UL2017-v1',
#            'SingleMuon_pancakes-02-withPF_Run2017E-09Aug2019_UL2017-v1',
#            'SingleMuon_pancakes-02-withPF_Run2017F-09Aug2019_UL2017-v1',
#            'SingleElectron_pancakes-02-withPF_Run2017B-09Aug2019_UL2017-v1',
#            'SingleElectron_pancakes-02-withPF_Run2017C-09Aug2019_UL2017-v1',
#            'SingleElectron_pancakes-02-withPF_Run2017D-09Aug2019_UL2017-v1',
#            'SingleElectron_pancakes-02-withPF_Run2017E-09Aug2019_UL2017-v1',
#            'SingleElectron_pancakes-02-withPF_Run2017F-09Aug2019_UL2017_rsb-v2',
#            'MET_pancakes-02-withPF_Run2017B-09Aug2019_UL2017_rsb-v1',
#            'MET_pancakes-02-withPF_Run2017C-09Aug2019_UL2017_rsb-v1',
#            'MET_pancakes-02-withPF_Run2017D-09Aug2019_UL2017_rsb-v1',
#            'MET_pancakes-02-withPF_Run2017E-09Aug2019_UL2017_rsb-v1',
#            'MET_pancakes-02-withPF_Run2017F-09Aug2019_UL2017_rsb-v1',
            'SingleMuon', 
            'SingleElectron', 
            #'WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8', 
            'WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8', 
            'WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8', 
            'WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8', 
            'WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8', 
            'WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8', 
            'WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8', 
        ],
        '2018':[
            'SingleMuon', 
            'EGamma', 
            #'WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8',
            'WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8', 
            'WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8', 
            'WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8', 
            'WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8', 
            'WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8', 
            'WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8', 
            'WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8', 
        ],
    },
    "preUL": {
        '2016':[
        ],
        '2017':[
        ],
        '2018':[
        ],
    },
}

# load files from datasets
for reco in recos:
    locdir = 'condor/' + label + '_' + year + '_' + reco
    logdir = locdir + '/logs/'
    os.system('mkdir -p  %s' % logdir)

    outdir = homedir + label + '_' + year + '_' + reco + '/outfiles/'
    os.system('mkdir -p  /eos/uscms/%s' % outdir)

    totfiles = {}
    #with open('fileset/fileset_%s_%s.json'%(year,reco), 'r') as f:
    with open('fileset/v2_2/%s.json'%(year), 'r') as f:
        newfiles = json.load(f)
        totfiles.update(newfiles)

    samplelist = []
    if reco in samples.keys():
        samplelist = samples[reco][year]
    else:
        samplelist = totfiles.keys()

    # write json to metadata.json
    with open("%s/metadata.json"%(locdir), 'w') as json_file:
         json.dump(totfiles, json_file, indent=4, sort_keys=True)

    # copy script to locdir
    os.system('cp %s %s'%(script,locdir))


    # submit jobs
    nsubmit = 0
    for sample in samplelist:    
        prefix = sample + '_%s_%s'%(year,reco)
        print('Submitting '+ prefix)
        
        njobs = int(len(totfiles[sample]) / files_per_job) + 1

        print(njobs)

        yearmod = ""
        if year[:4]=='2016':
            if year=='2016':
                yearmod=' --yearmod postVFP'
            else:
                yearmod=' --yearmod preVFP'
        
        for j in range(njobs):
            condor_templ_file = open(loc_base + "/condor/submit.templ.jdl")
            sh_templ_file     = open(loc_base + "/condor/submit.templ.sh")

            localcondor = '%s/%s_%i.jdl'%(locdir,prefix,j)
            condor_file = open(localcondor, "w")
            for line in condor_templ_file:
                line = line.replace('DIRECTORY', locdir)
                line = line.replace('PREFIX', prefix)
                line = line.replace('JOBID', str(j))
                line = line.replace('JSON', "%s/metadata.json"%(locdir))
                condor_file.write(line)
            condor_file.close()

            localsh = '%s/%s_%i.sh'%(locdir,prefix,j)
            eosoutput = 'root://cmseos.fnal.gov/%s/%s_%i.hist'%(outdir,prefix,j)
            sh_file = open(localsh, "w")
            for line in sh_templ_file:
                line = line.replace('SCRIPTNAME', script)
                line = line.replace('FILENUM', str(j))
                line = line.replace('YEAR', year[:4])
                line = line.replace('SAMPLE', sample)
                line = line.replace('PROCESSOR', 'htt')
                line = line.replace('STARTNUM', str(j * files_per_job))
                line = line.replace('ENDNUM', str((j + 1) * files_per_job))
                line = line.replace('ADDOPTS', addoptions+yearmod)
                line = line.replace('EOSOUT', eosoutput)
                line = line.replace('OUTDIR', outdir)
                sh_file.write(line)
            sh_file.close()

            os.system('chmod u+x %s'%localsh)
            if (os.path.exists('%s.log' % localcondor)):
                os.system('rm %s.log' % localcondor)
            condor_templ_file.close()
            sh_templ_file.close()
            
            print('To submit ', localcondor)
            os.system('condor_submit %s' % localcondor)

            nsubmit = nsubmit + 1

    print(nsubmit,"jobs submitted.")
