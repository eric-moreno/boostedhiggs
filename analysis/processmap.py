from collections import OrderedDict
from coffea import hist

process = hist.Cat("process", "Process", sorting='placement')
process_cat = "dataset"
process_map = OrderedDict()

process_map["wqq"] = [
    'WJetsToQQ_HT400to600_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8',
    'WJetsToQQ_HT600to800_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8',
    'WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8',
]
process_map["wlnu"] = [
    'WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8',
    'WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8',
    'WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8',
    'WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8',
    'WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8',
    'WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8',
    'WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8',
]
process_map["zqq"] = [
    'ZJetsToQQ_HT400to600_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8',
    'ZJetsToQQ_HT600to800_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8',
    'ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8',
]
process_map["qcd"] = [
    'QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8',
    'QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8',
    'QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8',
    'QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8',
    'QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8',
]
process_map["tt"] = [
    'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8',
    'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8',
    'TTToHadronic_TuneCP5_13TeV-powheg-pythia8',
]
process_map["st"] = [
    'ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8',
    'ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8',
    'ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8',
    'ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-madgraph-pythia8',
]
process_map["vv"] = [
    'WW_TuneCP5_13TeV-pythia8',
    'WZ_TuneCP5_13TeV-pythia8',
    'ZZ_TuneCP5_13TeV-pythia8',
]
process_map["h125"] = [
    'hww_mc',
    'hwwmc',
]

def apply(h):
    return h.group(process_cat, process, process_map)
