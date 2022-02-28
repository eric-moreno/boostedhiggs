from collections import OrderedDict
from coffea import hist

process = hist.Cat("process", "Process", sorting='placement')
process_cat = "dataset"
process_map = OrderedDict()

process_map["zll"] = [
    'DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8',
    'DYJetsToLL_M-50_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8',
    'DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8',
    'DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8',
    'DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8',
    'DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8',
    'DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8',
    'DYJetsToLL_Pt-100To250_TuneCP5_13TeV-amcatnloFXFX-pythia8',
    'DYJetsToLL_Pt-250To400_TuneCP5_13TeV-amcatnloFXFX-pythia8',
    'DYJetsToLL_Pt-400To650_TuneCP5_13TeV-amcatnloFXFX-pythia8',
    'DYJetsToLL_Pt-650ToInf_TuneCP5_13TeV-amcatnloFXFX-pythia8',
]
process_map["zee"] = [
    'DYJetsToLL_Pt-100To250_TuneCP5_13TeV-amcatnloFXFX-pythia8_Zee',
    'DYJetsToLL_Pt-250To400_TuneCP5_13TeV-amcatnloFXFX-pythia8_Zee',
    'DYJetsToLL_Pt-400To650_TuneCP5_13TeV-amcatnloFXFX-pythia8_Zee',
    'DYJetsToLL_Pt-650ToInf_TuneCP5_13TeV-amcatnloFXFX-pythia8_Zee',
]
process_map["zem"] = [
    'DYJetsToLL_Pt-100To250_TuneCP5_13TeV-amcatnloFXFX-pythia8_Zem',
    'DYJetsToLL_Pt-250To400_TuneCP5_13TeV-amcatnloFXFX-pythia8_Zem',
    'DYJetsToLL_Pt-400To650_TuneCP5_13TeV-amcatnloFXFX-pythia8_Zem',
    'DYJetsToLL_Pt-650ToInf_TuneCP5_13TeV-amcatnloFXFX-pythia8_Zem',
]
process_map["zmm"] = [
    'DYJetsToLL_Pt-100To250_TuneCP5_13TeV-amcatnloFXFX-pythia8_Zmm',
    'DYJetsToLL_Pt-250To400_TuneCP5_13TeV-amcatnloFXFX-pythia8_Zmm',
    'DYJetsToLL_Pt-400To650_TuneCP5_13TeV-amcatnloFXFX-pythia8_Zmm',
    'DYJetsToLL_Pt-650ToInf_TuneCP5_13TeV-amcatnloFXFX-pythia8_Zmm',
]
process_map["ztt"] = [
    'DYJetsToLL_Pt-100To250_TuneCP5_13TeV-amcatnloFXFX-pythia8_Ztt',
    'DYJetsToLL_Pt-250To400_TuneCP5_13TeV-amcatnloFXFX-pythia8_Ztt',
    'DYJetsToLL_Pt-400To650_TuneCP5_13TeV-amcatnloFXFX-pythia8_Ztt',
    'DYJetsToLL_Pt-650ToInf_TuneCP5_13TeV-amcatnloFXFX-pythia8_Ztt',
]
#process_map["zll-ht100to200"] = [
#    'DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8',
#]
#process_map["zll-ht1200to2500"] = [
#    'DYJetsToLL_M-50_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8',
#]
#process_map["zll-ht200to400"] = [
#    'DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8',
#]
#process_map["zll-ht2500toinf"] = [
#    'DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8',
#]
#process_map["zll-ht400to600"] = [
#    'DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8',
#]
#process_map["zll-ht600to800"] = [
#    'DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8',
#]
#process_map["zll-ht800to1200"] = [
#    'DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8',
#]
process_map["wjets"] = [
    'WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8',
    'WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8',
    'WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8',
    'WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8',
    'WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8',
    'WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8',
    'WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8',
    'boostedTau_WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8',
]
#process_map["wqq"] = [
#    'WJetsToQQ_HT400to600_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8',
#    'WJetsToQQ_HT600to800_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8',
#    'WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8',
#]
#process_map["zqq"] = [
#    'ZJetsToQQ_HT400to600_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8',
#    'ZJetsToQQ_HT600to800_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8',
#    'ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8',
#]
process_map["vqq"] = [
    'WJetsToQQ_HT400to600_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8',
    'WJetsToQQ_HT600to800_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8',
    'WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8',
    'ZJetsToQQ_HT400to600_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8',
    'ZJetsToQQ_HT600to800_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8',
    'ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8',
]
process_map["qcd"] = [
    'QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8',
    #'QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8',
    'QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8',
    'QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8',
    'QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8',
    'QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8',
    #'QCD_HT50to100_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8',
    'QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8',
    'QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8',
    'QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8',
    'QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8',
    'QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8',
    'QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8',
    'boostedTau_QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8',
]
process_map["tt-dilep"] = [
    'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8',
]
process_map["tt-semilep"] = [
    'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8',
]
process_map["tt-had"] = [
    'TTToHadronic_TuneCP5_13TeV-powheg-pythia8',
]
#process_map["tt"] = [
#    'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8',
#    'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8',
#    'TTToHadronic_TuneCP5_13TeV-powheg-pythia8',
#]
process_map["st"] = [
    'ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8',
    'ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8',
    'ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8',
    'ST_t-channel_top_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8',
    'ST_tW_antitop_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8',
    'ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8',
    'ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8',
    'ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8',
    'ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8',
    'ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8',
]
process_map["vv"] = [
    'WW_TuneCP5_13TeV-pythia8',
    'WZ_TuneCP5_13TeV-pythia8',
    'ZZ_TuneCP5_13TeV-pythia8',
]
process_map["h125"] = [
    'GluGluHToWWToLNuQQ_M125_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen710_pythia8',
    'GluGluHToWWToLNuQQ_M125_TuneCP5_PSweight_13TeV-powheg2-jhugen727-pythia8',
    'GluGluHToTauTau',
    'GluGluHToTauTau_M125_13TeV_powheg_pythia8',
    'VBFHToTauTau_M125_13TeV_powheg_pythia8',
    'WminusHToTauTau_M125_13TeV_powheg_pythia8',
    'WplusHToTauTau_M125_13TeV_powheg_pythia8',
    'ZHToTauTau_M125_13TeV_powheg_pythia8',
    'ggZH_HToTauTau_ZToLL_M125_13TeV_powheg_pythia8',
    'ggZH_HToTauTau_ZToNuNu_M125_13TeV_powheg_pythia8',
    'ggZH_HToTauTau_ZToQQ_M125_13TeV_powheg_pythia8',
    'ttHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8',
    'boostedTau_GluGluHTauTau_boostedTaua_13TeV_user',
]
process_map["phi10"] = [
    'Spin0ToTauTau_2j_scalar_g1_HT300_M10_nodmx_v0_TuneCP5_MLM',
]
process_map["phi20"] = [
    'Spin0ToTauTau_2j_scalar_g1_HT300_M20_nodmx_v0_TuneCP5_MLM',
]
process_map["phi30"] = [
    'Spin0ToTauTau_2j_scalar_g1_HT300_M30_nodmx_v0_TuneCP5_MLM',
]
process_map["phi40"] = [
    'Spin0ToTauTau_2j_scalar_g1_HT300_M40_nodmx_v0_TuneCP5_MLM',
]
process_map["phi50"] = [
    'Spin0ToTauTau_2j_scalar_g1_HT300_M50_nodmx_v0_TuneCP5_MLM',
]
process_map["phi75"] = [
    'Spin0ToTauTau_2j_scalar_g1_HT300_M75_nodmx_v0_TuneCP5_MLM',
]
process_map["phi100"] = [
    'Spin0ToTauTau_2j_scalar_g1_HT300_M100_nodmx_v0_TuneCP5_MLM',
]
process_map["phi125"] = [
    'Spin0ToTauTau_2j_scalar_g1_HT300_M125_nodmx_v0_TuneCP5_MLM',
]
process_map["phi150"] = [
    'Spin0ToTauTau_2j_scalar_g1_HT300_M150_nodmx_v0_TuneCP5_MLM',
]
process_map["phi200"] = [
    'Spin0ToTauTau_2j_scalar_g1_HT300_M200_nodmx_v0_TuneCP5_MLM',
]
process_map["phi250"] = [
    'Spin0ToTauTau_2j_scalar_g1_HT300_M250_nodmx_v0_TuneCP5_MLM',
]
process_map["phi300"] = [
    'Spin0ToTauTau_2j_scalar_g1_HT300_M300_nodmx_v0_TuneCP5_MLM',
]
#process_map["ggF-Htt"] = [
#    'GluGluHToTauTau_M125_13TeV_powheg_pythia8',
#]
#process_map["VBF-Htt"] = [
#    'VBFHToTauTau_M125_13TeV_powheg_pythia8',
#]
#process_map["Wm-Htt"] = [
#    'WminusHToTauTau_M125_13TeV_powheg_pythia8',
#]
#process_map["Wp-Htt"] = [
#    'WplusHToTauTau_M125_13TeV_powheg_pythia8',
#]
#process_map["ZH-Htt"] = [
#    'ZHToTauTau_M125_13TeV_powheg_pythia8',
#]
#process_map["ggZll-Htt"] = [
#    'ggZH_HToTauTau_ZToLL_M125_13TeV_powheg_pythia8',
#]
#process_map["ggZvv-Htt"] = [
#    'ggZH_HToTauTau_ZToNuNu_M125_13TeV_powheg_pythia8',
#]
#process_map["ggZqq-Htt"] = [
#    'ggZH_HToTauTau_ZToQQ_M125_13TeV_powheg_pythia8',
#]
#process_map["tt-Htt"] = [
#    'ttHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8',
#]
process_map["data"] = [
    'SingleMuon',
    'SingleElectron',
    'EGamma',
    'MET',
#2017
    #'JetHT',
    'JetHT_pancakes-02_Run2017B-09Aug2019_UL2017-v1',
    'JetHT_pancakes-02_Run2017C-09Aug2019_UL2017-v1',
    'JetHT_pancakes-02_Run2017D-09Aug2019_UL2017-v1',
    'JetHT_pancakes-02_Run2017E-09Aug2019_UL2017-v1',
    'JetHT_pancakes-02_Run2017F-09Aug2019_UL2017-v1',
    #'SingleElectron',
    'SingleElectron_pancakes-02-withPF_Run2017B-09Aug2019_UL2017-v1',
    'SingleElectron_pancakes-02-withPF_Run2017C-09Aug2019_UL2017-v1',
    'SingleElectron_pancakes-02-withPF_Run2017D-09Aug2019_UL2017-v1',
    'SingleElectron_pancakes-02-withPF_Run2017E-09Aug2019_UL2017-v1',
    'SingleElectron_pancakes-02-withPF_Run2017F-09Aug2019_UL2017_rsb-v2',
    #'SingleMuon',
    'SingleMuon_pancakes-02-withPF_Run2017B-09Aug2019_UL2017-v1',
    'SingleMuon_pancakes-02-withPF_Run2017C-09Aug2019_UL2017-v1',
    'SingleMuon_pancakes-02-withPF_Run2017D-09Aug2019_UL2017-v1',
    'SingleMuon_pancakes-02-withPF_Run2017E-09Aug2019_UL2017-v1',
    'SingleMuon_pancakes-02-withPF_Run2017F-09Aug2019_UL2017-v1',
    'Tau',
    #'MET_pancakes-02-withPF_Run2017B-09Aug2019_UL2017_rsb-v1',
    'MET_pancakes-02-withPF_Run2017C-09Aug2019_UL2017_rsb-v1',
    'MET_pancakes-02-withPF_Run2017D-09Aug2019_UL2017_rsb-v1',
    'MET_pancakes-02-withPF_Run2017E-09Aug2019_UL2017_rsb-v1',
    'MET_pancakes-02-withPF_Run2017F-09Aug2019_UL2017_rsb-v1',
#2018
    #'SingleMuon',
    'SingleMuon_pancakes-02_Run2018A-12Nov2019_UL2018_rsb-v1',
    'SingleMuon_pancakes-02_Run2018B-12Nov2019_UL2018-v2',
    'SingleMuon_pancakes-02_Run2018C-12Nov2019_UL2018-v2',
    'SingleMuon_pancakes-02_Run2018D-12Nov2019_UL2018-v4',
    #'MET',
    'MET_pancakes-02_Run2018A-12Nov2019_UL2018-v3',
    'MET_pancakes-02_Run2018B-12Nov2019_UL2018_rsb-v1',
    'MET_pancakes-02_Run2018C-12Nov2019_UL2018_rsb-v1',
    'MET_pancakes-02_Run2018D-12Nov2019_UL2018_rsb-v2',
    #'EGamma',
    'EGamma_pancakes-02_Run2018A-12Nov2019_UL2018-v2',
    'EGamma_pancakes-02_Run2018B-12Nov2019_UL2018-v2',
    'EGamma_pancakes-02_Run2018C-12Nov2019_UL2018-v2',
    'EGamma_pancakes-02_Run2018D-12Nov2019_UL2018-v4',
]

def apply(h):
    return h.group(process_cat, process, process_map)
