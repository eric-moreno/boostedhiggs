import argparse
from audioop import add
import os
import re
import fileinput

import json
import glob
import sys

'''
 Submit condor jobs of processor
 Run as e.g.: python condor/submit.py Jun06 run.py 4 2017
              python condor/submit.py Jun06_Opt2 run.py 4 2017 10000 "--plotopt 2"
              python condor/submit.py Jun21 run.py 4 2016
              python condor/submit.py Jun21 run.py 4 2016APV
              python condor/submit.py Nov04 run.py 4 2018
 Arguments:
  = [0]: tag of jobs and output dir in eos e.g. Jul1
  - [1]: script to run e.g. run.py (needs to be included in transfer_files in templ.jdl)
  - [2]: number of files per job e.g. 20
  - [3]: year
  - [4]: additional options for processor
'''
# Note: change username in `drankin` in this script

homedir = "/store/user/eamoreno/boostedhiggs/"

parser = argparse.ArgumentParser(description='Submit coffea condor jobs.')
parser.add_argument('settings', metavar='S', type=str, nargs='+', help='label scriptname files_per_job year chunksize options')
args = parser.parse_args()

if (not ((len(args.settings) == 4) or (len(args.settings) == 5) or (len(args.settings) == 6))):
    print("Wrong number of arguments (must be 3, 4, or 5, found", len(args.settings), ")")
    sys.exit()

label = args.settings[0]
script = args.settings[1]  # should be run.py
files_per_job = int(args.settings[2])
year = args.settings[3]
nworkers = 1
chunksize = 10000
if len(args.settings) == 5:
    chunksize = args.settings[4]
addoptions = ""
if len(args.settings) == 6:
    addoptions = args.settings[5]

loc_base = os.environ['PWD']

# list of samples to run and recos to use
recos = ["UL"] 

# if empty run over all the datasets in both filesets
samples = {
    "UL": {
        '2016APV':[
            'ZJetsToQQ_HT-400to600',
            'ZJetsToQQ_HT-800toInf',
            'ZJetsToQQ_HT-600to800',
            'ST_s-channel_4f_leptonDecays',
            'ST_s-channel_4f_hadronicDecays',
            'ST_t-channel_antitop_4f_InclusiveDecays',
            'ST_t-channel_top_4f_InclusiveDecays',
            'ST_tW_antitop_5f_inclusiveDecays',
            'ST_tW_top_5f_inclusiveDecays',
            'WJetsToQQ_HT-600to800',
            'WJetsToQQ_HT-800toInf',
            'WJetsToQQ_HT-400to600',
            'WJetsToLNu_HT-100To200',
            'WJetsToLNu_HT-1200To2500',
            'WJetsToLNu_HT-800To1200',
            'WJetsToLNu_HT-600To800',
            'WJetsToLNu_HT-2500ToInf',
            'WJetsToLNu_HT-200To400',
            #'WJetsToLNu_HT-70To100', #deleted in v2_3 apparently we dont need this though 
            'WJetsToLNu_HT-400To600',
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
            'TTToSemiLeptonic',
            'TTToHadronic',
            'TTTo2L2Nu',
            'DYJetsToLL_Pt-50To100',
            'DYJetsToLL_Pt-100To250',
            'DYJetsToLL_Pt-250To400',
            'DYJetsToLL_Pt-400To650',
            'DYJetsToLL_Pt-650ToInf',
            'WW',
            'WZ',
            'ZZ',
            'GluGluHToTauTau',
            'VBFHToTauTau',
            'WminusHToTauTau',
            'WplusHToTauTau',
            'ZHToTauTau',
            'ttHToTauTau',
            'MET_Run2016B_ver2_HIPM',
            'MET_Run2016C_HIPM',
            'MET_Run2016D_HIPM',
            'MET_Run2016E_HIPM',
            'MET_Run2016F_HIPM',
            'SingleElectron_Run2016B_ver2_HIPM',
            'SingleElectron_Run2016C_HIPM',
            'SingleElectron_Run2016D_HIPM',
            'SingleElectron_Run2016E_HIPM',
            'SingleElectron_Run2016F_HIPM',
            'SingleMuon_Run2016B_ver2_HIPM',
            'SingleMuon_Run2016C_HIPM',
            'SingleMuon_Run2016D_HIPM',
            'SingleMuon_Run2016E_HIPM',
            'SingleMuon_Run2016F_HIPM',
            #'JetHT_Run2016B_ver2_HIPM',
            #'JetHT_Run2016C_HIPM',
            #'JetHT_Run2016D_HIPM',
            #'JetHT_Run2016E_HIPM',
            #'JetHT_Run2016F_HIPM',
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
        '2016':[
            'ZJetsToQQ_HT-400to600',
            'ZJetsToQQ_HT-800toInf',
            'ZJetsToQQ_HT-600to800',
            'ST_s-channel_4f_leptonDecays',
            'ST_s-channel_4f_hadronicDecays',
            'ST_t-channel_antitop_4f_InclusiveDecays',
            'ST_t-channel_top_4f_InclusiveDecays',
            'ST_tW_antitop_5f_inclusiveDecays',
            'ST_tW_top_5f_inclusiveDecays',
            'WJetsToQQ_HT-600to800',
            'WJetsToQQ_HT-800toInf',
            'WJetsToQQ_HT-400to600',
            'WJetsToLNu_HT-100To200',
            'WJetsToLNu_HT-1200To2500',
            'WJetsToLNu_HT-800To1200',
            'WJetsToLNu_HT-600To800',
            'WJetsToLNu_HT-2500ToInf',
            'WJetsToLNu_HT-200To400',
            #'WJetsToLNu_HT-70To100', #deleted in v2_3 apparently we dont need this though 
            'WJetsToLNu_HT-400To600',
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
            'TTToSemiLeptonic',
            'TTToHadronic',
            'TTTo2L2Nu',
            'DYJetsToLL_Pt-50To100',
            'DYJetsToLL_Pt-100To250',
            'DYJetsToLL_Pt-250To400',
            'DYJetsToLL_Pt-400To650',
            'DYJetsToLL_Pt-650ToInf',
            'WW',
            'WZ',
            'ZZ',
            'GluGluHToTauTau',
            'VBFHToTauTau',
            'WminusHToTauTau',
            'WplusHToTauTau',
            'ZHToTauTau',
            'ttHToTauTau',
            'MET_Run2016F',
            'MET_Run2016G',
            'MET_Run2016H',
            'SingleElectron_Run2016F',
            'SingleElectron_Run2016G',
            'SingleElectron_Run2016H',
            'SingleMuon_Run2016F',
            'SingleMuon_Run2016G',
            'SingleMuon_Run2016H',
            #'JetHT_Run2016F',
            #'JetHT_Run2016G',
            #'JetHT_Run2016H',
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
        '2017':[
            'ZJetsToQQ_HT-400to600',
            'ZJetsToQQ_HT-800toInf',
            'ZJetsToQQ_HT-600to800',
            'ST_s-channel_4f_leptonDecays',
            'ST_s-channel_4f_hadronicDecays',
            'ST_t-channel_antitop_4f_InclusiveDecays',
            'ST_t-channel_top_4f_InclusiveDecays',
            'ST_tW_antitop_5f_inclusiveDecays',
            'ST_tW_top_5f_inclusiveDecays',
            'WJetsToQQ_HT-600to800',
            'WJetsToQQ_HT-800toInf',
            'WJetsToQQ_HT-400to600',
            'WJetsToLNu_HT-100To200',
            'WJetsToLNu_HT-1200To2500',
            'WJetsToLNu_HT-800To1200',
            'WJetsToLNu_HT-600To800',
            'WJetsToLNu_HT-2500ToInf',
            'WJetsToLNu_HT-200To400',
            'WJetsToLNu_HT-70To100',
            'WJetsToLNu_HT-400To600',
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
            'TTToSemiLeptonic',
            'TTToHadronic',
            'TTTo2L2Nu',
            'DYJetsToLL_Pt-50To100',
            'DYJetsToLL_Pt-100To250',
            'DYJetsToLL_Pt-250To400',
            'DYJetsToLL_Pt-400To650',
            'DYJetsToLL_Pt-650ToInf',
            'WW',
            'WZ',
            'ZZ',
            'SingleElectron_Run2017B',
            'SingleElectron_Run2017C',
            'SingleElectron_Run2017D',
            'SingleElectron_Run2017E',
            'SingleElectron_Run2017F',
            'SingleMuon_Run2017B',
            'SingleMuon_Run2017C',
            'SingleMuon_Run2017D',
            'SingleMuon_Run2017E',
            'SingleMuon_Run2017F',
            'MET_Run2017B',
            'MET_Run2017C',
            'MET_Run2017D',
            'MET_Run2017E',
            'MET_Run2017F',
            'GluGluHToTauTau',
            'VBFHToTauTau',
            'WminusHToTauTau',
            'WplusHToTauTau',
            'ZHToTauTau',
            'ttHToTauTau',
#            #'JetHT_Run2017B',
#            #'JetHT_Run2017C',
#            #'JetHT_Run2017D',
#            #'JetHT_Run2017E',
#            #'JetHT_Run2017F',
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
        '2018':[
            'ZJetsToQQ_HT-400to600',
            'ZJetsToQQ_HT-800toInf',
            'ZJetsToQQ_HT-600to800',
            'ST_s-channel_4f_leptonDecays',
            'ST_s-channel_4f_hadronicDecays',
            'ST_t-channel_antitop_4f_InclusiveDecays',
            'ST_t-channel_top_4f_InclusiveDecays',
            'ST_tW_antitop_5f_inclusiveDecays',
            'ST_tW_top_5f_inclusiveDecays',
            'WJetsToQQ_HT-600to800',
            'WJetsToQQ_HT-800toInf',
            'WJetsToQQ_HT-400to600',
            'WJetsToLNu_HT-100To200',
            'WJetsToLNu_HT-1200To2500',
            'WJetsToLNu_HT-800To1200',
            'WJetsToLNu_HT-600To800',
            'WJetsToLNu_HT-2500ToInf',
            'WJetsToLNu_HT-200To400',
            'WJetsToLNu_HT-70To100',
            'WJetsToLNu_HT-400To600',
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
            'TTToSemiLeptonic',
            'TTToHadronic',
            'TTTo2L2Nu',
            'DYJetsToLL_Pt-50To100',
            'DYJetsToLL_Pt-100To250',
            'DYJetsToLL_Pt-250To400',
            'DYJetsToLL_Pt-400To650',
            'DYJetsToLL_Pt-650ToInf',
            'WW',
            'WZ',
            'ZZ',
            'GluGluHToTauTau',
            'VBFHToTauTau',
            'WminusHToTauTau',
            'WplusHToTauTau',
            'ZHToTauTau',
            'ttHToTauTau',
            'SingleMuon_Run2018A',
            'SingleMuon_Run2018B',
            'SingleMuon_Run2018C',
            'SingleMuon_Run2018D',
            'EGamma_Run2018A',
            'EGamma_Run2018B',
            'EGamma_Run2018C',
            'EGamma_Run2018D',
            'MET_Run2018A',
            'MET_Run2018B',
            'MET_Run2018C',
            'MET_Run2018D',
            #'JetHT_Run2018A',
            #'JetHT_Run2018B',
            #'JetHT_Run2018C',
            #'JetHT_Run2018D',
            #'Spin0ToTauTau_2j_scalar_g1_HT300_M10_nodmx_v0_TuneCP5_MLM',
            #'Spin0ToTauTau_2j_scalar_g1_HT300_M20_nodmx_v0_TuneCP5_MLM',
            #'Spin0ToTauTau_2j_scalar_g1_HT300_M30_nodmx_v0_TuneCP5_MLM',
            #'Spin0ToTauTau_2j_scalar_g1_HT300_M40_nodmx_v0_TuneCP5_MLM',
            #'Spin0ToTauTau_2j_scalar_g1_HT300_M50_nodmx_v0_TuneCP5_MLM',
            #'Spin0ToTauTau_2j_scalar_g1_HT300_M75_nodmx_v0_TuneCP5_MLM',
            #'Spin0ToTauTau_2j_scalar_g1_HT300_M100_nodmx_v0_TuneCP5_MLM',
            #'Spin0ToTauTau_2j_scalar_g1_HT300_M125_nodmx_v0_TuneCP5_MLM',
            #'Spin0ToTauTau_2j_scalar_g1_HT300_M150_nodmx_v0_TuneCP5_MLM',
            #'Spin0ToTauTau_2j_scalar_g1_HT300_M200_nodmx_v0_TuneCP5_MLM',
            #'Spin0ToTauTau_2j_scalar_g1_HT300_M250_nodmx_v0_TuneCP5_MLM',
            #'Spin0ToTauTau_2j_scalar_g1_HT300_M300_nodmx_v0_TuneCP5_MLM',
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
    with open('fileset/v2_2.5/%s.json'%(year), 'r') as f:
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
        
        njobs = int((len(totfiles[sample])-1) / files_per_job) + 1

        print(njobs)

        yearmod = ""
        #if year[:4]=='2016':
        #    if year=='2016':
        #        yearmod=' --yearmod postVFP'
        #    else:
        #        yearmod=' --yearmod preVFP'
        
        for j in range(njobs):
            condor_templ_file = open(loc_base + "/condor/submit.templ.jdl")
            sh_templ_file     = open(loc_base + "/condor/submit.templ.sh")

            localcondor = '%s/%s_%i.jdl'%(locdir,prefix,j)
            condor_file = open(localcondor, "w")
            for line in condor_templ_file:
                line = line.replace('DIRECTORY', locdir)
                line = line.replace('PREFIX', prefix)
                line = line.replace('JOBID', str(j))
                line = line.replace('MEMREQUEST', "7500" if "DYJets" in sample or "HTauTau" in sample else "6000")
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
                line = line.replace('ADDOPTS', addoptions+yearmod+" --chunksize %i"%(chunksize/(2 if "Spin0" in sample else 1))+" --nworkers %i"%nworkers)
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
