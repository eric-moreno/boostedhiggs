import json

xs = {
    "WW_TuneCP5_13TeV-pythia8": 75.83,
    "WZ_TuneCP5_13TeV-pythia8": 27.56,
    "ZZ_TuneCP5_13TeV-pythia8": 12.14,
    "QCD_HT50to100_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8": 248600000.,
    "QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8": 27990000.,
    "QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8": 9999999.,#FIXME
    "QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8": 322600.,# huge weights
    "QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8": 29980.,
    "QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8": 29980.,
    "QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8": 6334.,
    "QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8": 6334.,
    "QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8": 1088.,
    "QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8": 1088.,
    "QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8": 99.11,
    "QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8": 99.11,
    "QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8": 20.23,
    "QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8": 20.23,
    "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8": 6225.42,
    "DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8": 140.1,
    "DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8": 38.35,
    "DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8": 5.217,
    "DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8": 1.267,
    "DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8": 0.5682,
    "DYJetsToLL_M-50_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8": 0.1332,
    "DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8": 0.002978,
    "DYJetsToLL_Pt-100To250_TuneCP5_13TeV-amcatnloFXFX-pythia8": 89.39,
    "DYJetsToLL_Pt-250To400_TuneCP5_13TeV-amcatnloFXFX-pythia8": 3.43,
    "DYJetsToLL_Pt-400To650_TuneCP5_13TeV-amcatnloFXFX-pythia8": 0.464,
    "DYJetsToLL_Pt-650ToInf_TuneCP5_13TeV-amcatnloFXFX-pythia8": 0.0436,
    "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8": 52940.,
    "WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8": 1270.,#FIXME
    "WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8": 1252.,
    "WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8": 336.5,
    "WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8": 45.12,
    "WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8": 10.99,
    "WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8": 4.938,
    "WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8": 1.155,
    "WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8": 0.02625,
    "WJetsToQQ_HT400to600_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8": 277.,
    "WJetsToQQ_HT600to800_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8": 59.06,
    "WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8": 28.75,
    "ZJetsToQQ_HT400to600_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8": 114.5,
    "ZJetsToQQ_HT600to800_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8": 25.41,
    "ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8": 12.91,
    "ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8": 11.24,
    "ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-madgraph-pythia8": 11.24,
    "ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8": 3.74,
    "ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-madgraph-pythia8": 3.74,
    "ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8": 26.2278,
    "ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8": 26.2278,
    "ST_t-channel_top_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8": 44.07048,
    "ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8": 44.07048,
    "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8": 35.6,
    "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8": 35.6,
    "ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8": 35.6,
    "ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8": 35.6,
    "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8": 77.10979749999998,
    "TTToHadronic_TuneCP5_13TeV-powheg-pythia8": 303.8527975,
    "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8": 306.137405,
    "GluGluHToTauTau_M125_13TeV_powheg_pythia8": 48.58*0.0627,
    "VBFHToTauTau_M125_13TeV_powheg_pythia8":3.782*0.0627,
    "WminusHToTauTau_M125_13TeV_powheg_pythia8":0.5328*0.0627,
    "WplusHToTauTau_M125_13TeV_powheg_pythia8":0.8400*0.0627,
    "ZHToTauTau_M125_13TeV_powheg_pythia8":0.7612*0.0627,
    "ggZH_HToTauTau_ZToLL_M125_13TeV_powheg_pythia8":0.1227*0.0627*3*0.033658,
    "ggZH_HToTauTau_ZToNuNu_M125_13TeV_powheg_pythia8":0.1227*0.0627*0.2000,
    "ggZH_HToTauTau_ZToQQ_M125_13TeV_powheg_pythia8":0.1227*0.0627*0.6991,
    "ttHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8":0.5269*0.0627,
    "GluGluHToTauTau": 0.8045*0.0627, #from lhe (x BR)
    "GluGluHToWWToLNuQQ_M125_NNPDF31_TuneCP5_PSweights_13TeV_powheg_JHUGen710_pythia8":48.58*0.2137*2.*0.676*0.108,
    "GluGluHToWWToLNuQQ_M125_TuneCP5_PSweight_13TeV-powheg2-jhugen727-pythia8":48.58*0.2137*2.*0.676*0.108,
    "Spin0ToTauTau_2j_scalar_g1_HT300_M10_nodmx_v0_TuneCP5_MLM": 9.85e-02,
    "Spin0ToTauTau_2j_scalar_g1_HT300_M30_nodmx_v0_TuneCP5_MLM": 1.78e-02,
    "Spin0ToTauTau_2j_scalar_g1_HT300_M50_nodmx_v0_TuneCP5_MLM": 1.61e-02,
    "Spin0ToTauTau_2j_scalar_g1_HT300_M125_nodmx_v0_TuneCP5_MLM": 1.28e-02,
    "Spin0ToTauTau_2j_scalar_g1_HT300_M300_nodmx_v0_TuneCP5_MLM": 8.34e-03,

}
zkeys = []
for d in xs.keys():
    if "DYJetsToLL_Pt" in d:
       zkeys.append(d)
for z in zkeys:
    xs[z+"_Zee"] = xs[z]
    xs[z+"_Zmm"] = xs[z]
    xs[z+"_Zem"] = xs[z]
    xs[z+"_Ztt"] = xs[z]

with open('fileset/xsecs.json', 'w') as outfile:
    json.dump(xs, outfile, indent=4)
