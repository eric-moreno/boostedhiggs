import os
import numpy as np
import awkward as ak
from coffea.util import load

compiled = load(os.path.join(os.path.dirname(__file__), 'data', 'corrections.coffea'))
compiled_trigger = load(os.path.join(os.path.dirname(__file__), 'data', 'trig_sf_corr.coffea'))

puW2017_nonUL = np.array([0.183454, 3.93313, 3.47111, 2.4924, 1.62495, 1.5151, 1.28791, 1.27825, 0.615845, 1.45208, 1.49768, 1.48747, 1.33116, 1.1645, 1.07901, 1.05437, 1.08066, 1.12907, 1.16584, 1.18936, 1.21284, 1.23849, 1.25967, 1.27099, 1.2727, 1.2713, 1.27087, 1.26652, 1.27449, 1.25163, 1.22131, 1.16954, 1.10903, 1.03816, 0.969432, 0.91188, 0.86681, 0.834505, 0.788379, 0.750796, 0.759152, 0.793175, 0.858347, 0.958795, 1.09432, 1.25564, 1.41965, 1.49593, 1.53054, 1.4622, 1.33676, 1.15439, 0.950556, 0.749702, 0.569789, 0.410747, 0.290007, 0.198655, 0.137644, 0.096639, 0.0692314, 0.0508587, 0.038443, 0.0299595, 0.0240949, 0.01712, 0.0124798, 0.0107683, 0.0095972, 0.00881241, 0.00830987, 0.00801759, 0.00788476, 0.0078747, 0.00633731, 0.00533369, 0.00544371, 0.00558365, 0.00574411, 0.00591599, 0.00609007, 0.00625711, 0.00640816, 0.00509462, 0.00422803, 0.00425915, 0.00426558, 0.00424666, 0.00420325, 0.00413738, 0.00405193, 0.00395019, 0.00294303, 0.00229603, 0.00220003, 0.00210138, 0.00200166, 0.00190214, 0.0018038]) #for non-UL
notUL2017 = [
    "WJetsToQQ_HT400to600_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8",
    "WJetsToQQ_HT600to800_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8",
    "WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8",
    "ZJetsToQQ_HT600to800_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8",
    "ZJetsToQQ_HT400to600_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8",
    "ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8",
    "ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8",
    "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8",
    "ST_t-channel_top_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8",
    "ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8",
    "ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8",
    "ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8",
    "WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8",
    "WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8",
    "WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8",
    "WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8",
    "WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8",
    "WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8",
    "WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8",
    "WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8",
    "DYJetsToLL_Pt-250To400_TuneCP5_13TeV-amcatnloFXFX-pythia8",
    "DYJetsToLL_Pt-400To650_TuneCP5_13TeV-amcatnloFXFX-pythia8",
    "DYJetsToLL_Pt-650ToInf_TuneCP5_13TeV-amcatnloFXFX-pythia8",
    "GluGluHToTauTau_M125_13TeV_powheg_pythia8",
    "VBFHToTauTau_M125_13TeV_powheg_pythia8",
    "WminusHToTauTau_M125_13TeV_powheg_pythia8",
    "WplusHToTauTau_M125_13TeV_powheg_pythia8",
    "ZHToTauTau_M125_13TeV_powheg_pythia8",
    "ggZH_HToTauTau_ZToLL_M125_13TeV_powheg_pythia8",
    "ggZH_HToTauTau_ZToNuNu_M125_13TeV_powheg_pythia8",
    "ggZH_HToTauTau_ZToQQ_M125_13TeV_powheg_pythia8",
    "ttHToTauTau_M125_TuneCP5_13TeV-powheg",
]

def _msoftdrop_weight(pt, eta):
    gpar = np.array([1.00626, -1.06161, 0.0799900, 1.20454])
    cpar = np.array([1.09302, -0.000150068, 3.44866e-07, -2.68100e-10, 8.67440e-14, -1.00114e-17])
    fpar = np.array([1.27212, -0.000571640, 8.37289e-07, -5.20433e-10, 1.45375e-13, -1.50389e-17])
    genw = gpar[0] + gpar[1]*np.power(pt*gpar[2], -gpar[3])
    ptpow = np.power.outer(pt, np.arange(cpar.size))
    cenweight = np.dot(ptpow, cpar)
    forweight = np.dot(ptpow, fpar)
    weight = np.where(np.abs(eta) < 1.3, cenweight, forweight)
    return genw*weight


def corrected_msoftdrop(fatjets):
    sf_flat = _msoftdrop_weight(fatjets.pt.flatten(), fatjets.eta.flatten())
    sf_flat = np.maximum(1e-5, sf_flat)
    try:
        # old pancakes
        dazsle_msd = fatjets.msoftdrop_raw
    except AttributeError:
        dazsle_msd = (fatjets.subjets * (1 - fatjets.subjets.rawFactor)).sum().mass
    return dazsle_msd * ak.JaggedArray.fromoffsets(fatjets.array.offsets, sf_flat)


def n2ddt_shift(fatjets, year='2017'):
    return compiled[f'{year}_n2ddt_rho_pt'](fatjets.rho, fatjets.pt)


def add_pileup_weight(weights, nPU, year='2017', dataset=None):
    if year == '2017' and dataset in compiled['2017_pileupweight_dataset']:
        weights.add(
            'pileup_weight',
            compiled['2017_pileupweight_dataset'][dataset](nPU),
            compiled['2017_pileupweight_dataset_puUp'][dataset](nPU),
            compiled['2017_pileupweight_dataset_puDown'][dataset](nPU),
        )
    elif year == '2017' and dataset in notUL2017:
        weights.add('pileup_weight',
            puW2017_nonUL[np.clip(nPU,0,len(puW2017_nonUL)-1)]
        )
    else:
        weights.add(
            'pileup_weight',
            compiled[f'{year}_pileupweight'](nPU),
            compiled[f'{year}_pileupweight_puUp'](nPU),
            compiled[f'{year}_pileupweight_puDown'](nPU),
        )


def add_VJets_NLOkFactor(weights, genBosonPt, year, dataset):
    if year == '2017' and 'ZJetsToQQ_HT' in dataset:
        nlo_over_lo_qcd = compiled['2017_Z_nlo_qcd'](genBosonPt)
        nlo_over_lo_ewk = compiled['Z_nlo_over_lo_ewk'](genBosonPt)
    elif year == '2017' and 'WJetsToQQ_HT' in dataset:
        nlo_over_lo_qcd = compiled['2017_W_nlo_qcd'](genBosonPt)
        nlo_over_lo_ewk = compiled['W_nlo_over_lo_ewk'](genBosonPt)
    elif year == '2016' and 'DYJetsToQQ' in dataset:
        nlo_over_lo_qcd = compiled['2016_Z_nlo_qcd'](genBosonPt)
        nlo_over_lo_ewk = compiled['Z_nlo_over_lo_ewk'](genBosonPt)
    elif year == '2016' and 'WJetsToQQ' in dataset:
        nlo_over_lo_qcd = compiled['2016_W_nlo_qcd'](genBosonPt)
        nlo_over_lo_ewk = compiled['W_nlo_over_lo_ewk'](genBosonPt)
    else:
        return
    weights.add('VJets_NLOkFactor', nlo_over_lo_qcd * nlo_over_lo_ewk)


def add_jetTriggerWeight(weights, jet_msd, jet_pt, year):
    jet_msd = jet_msd.pad(1, clip=True).fillna(0).flatten()
    jet_pt = jet_pt.pad(1, clip=True).fillna(0).flatten()
    nom = compiled[f'{year}_trigweight_msd_pt'](jet_msd, jet_pt)
    up = compiled[f'{year}_trigweight_msd_pt_trigweightUp'](jet_msd, jet_pt)
    down = compiled[f'{year}_trigweight_msd_pt_trigweightDown'](jet_msd, jet_pt)
    weights.add('jet_trigger', nom, up, down)

def add_TriggerWeight(weights, jet_msd, jet_pt, lep_pt, year, channel):

    jet_msd = jet_msd.pad(1, clip=True).fillna(0).flatten()
    jet_pt = jet_pt.pad(1, clip=True).fillna(0).flatten()
    lep_pt = lep_pt.pad(1, clip=True).fillna(0).flatten()
    if (channel=="hadhad"):
        nom = compiled_trigger[f'{year}_trigsf_hadhad_nom/ratio_value'](jet_pt,jet_msd)
        up = compiled_trigger[f'{year}_trigsf_hadhad_up/ratio_value'](jet_pt,jet_msd)
        down = compiled_trigger[f'{year}_trigsf_hadhad_down/ratio_value'](jet_pt,jet_msd)
    if (channel=="hadel"):
        nom = compiled_trigger[f'{year}_trigsf_hadel_nom/ratio_value'](jet_pt,lep_pt)
        up = compiled_trigger[f'{year}_trigsf_hadel_up/ratio_value'](jet_pt,lep_pt)
        down = compiled_trigger[f'{year}_trigsf_hadel_down/ratio_value'](jet_pt,lep_pt)
    if (channel=="hadmu"):
        nom = compiled_trigger[f'{year}_trigsf_hadmu_nom/ratio_value'](jet_pt,lep_pt)
        up = compiled_trigger[f'{year}_trigsf_hadmu_up/ratio_value'](jet_pt,lep_pt)
        down = compiled_trigger[f'{year}_trigsf_hadmu_down/ratio_value'](jet_pt,lep_pt)
    #up = compiled[f'{year}_trigweight_msd_pt_trigweightUp'](jet_msd, jet_pt)
    #down = compiled[f'{year}_trigweight_msd_pt_trigweightDown'](jet_msd, jet_pt)
    weights.add('%s_trigger'%channel, nom, up, down)

def add_TopPtReweighting(weights, topPt, year, dataset):
#$SF(p_T)=e^{0.0615-0.0005\cdot p_T}$ for data/POWHEG+Pythia8
    if 'TT' in dataset:
        toppt_weight1 = np.exp(0.0615-0.0005*np.clip(topPt[:,0],0.,500.))
        toppt_weight2 = np.exp(0.0615-0.0005*np.clip(topPt[:,1],0.,500.))
        weights.add('TopPtReweight', toppt_weight1 * toppt_weight2)
