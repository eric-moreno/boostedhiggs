import os
import numpy as np
import awkward as ak
from coffea.util import load
from coffea.lookup_tools import extractor
import json

compiled = load(os.path.join(os.path.dirname(__file__), 'data', 'corrections.coffea'))
compiled_trigger = load(os.path.join(os.path.dirname(__file__), 'data', 'trig_sf_corr.coffea'))

#from https://www.research-collection.ethz.ch/bitstream/handle/20.500.11850/240924/1/s10052-017-5389-1.pdf
Vpt_corr_bins = np.array([
    30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 950.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1600.0, 1800.0, 2000.0, 2200.0, 2400.0, 2600.0, 2800.0, 3000.0, 6500.0
])
Vpt_corr_value = np.array([
    1.003846859, 1.005303426, 0.9947071550000002, 0.98928453, 0.9850576999999998, 0.9813048799999999, 0.9778652800000001, 0.97445919, 0.9710953100000002, 0.9677317199999999, 0.9642441799999999, 0.9608947800000001, 0.9525394399999999, 0.9365518699999998, 0.9219756100000001, 0.9087814499999999, 0.8968980000000001, 0.8859843000000001, 0.8760160999999999, 0.8671511999999999, 0.8583182999999999, 0.8507631, 0.8431063999999999, 0.8362730999999999, 0.8303883000000001, 0.8242617999999999, 0.8194368000000001, 0.8138635, 0.8077816999999999, 0.8017448, 0.7931413999999999, 0.7852697, 0.7775183999999999, 0.7692830000000002, 0.7558978, 0.7443137, 0.7334527, 0.7233202, 0.7140055999999999, 0.7045699000000001, 0.691076, 0.6890300000000001
])

lepton_exts = {}
# several histograms can be imported at once using wildcards (*)
#available evaluator keys (elec_RECO):
#	 EGamma_SF2D
#	 EGamma_SF2D_error
#            eta - pt
#available evaluator keys (elec_ID):
#	 EGamma_SF2D
#	 EGamma_SF2D_error
#            eta - pt
#available evaluator keys (elec_TRIG32):
#	 EGamma_SF2D
#	 EGamma_SF2D_error
#            eta - pt
#available evaluator keys (elec_TRIG115):
#	 Ele115_PtEtaBins/abseta_pt_SF_value
#	 Ele115_PtEtaBins/abseta_pt_SF_error
#available evaluator keys (muon_ISO):
#	 NUM_LooseRelIso_DEN_MediumID/abseta_pt_value
#	 NUM_LooseRelIso_DEN_MediumID/abseta_pt_error
#available evaluator keys (muon_ID):
#	 NUM_MediumID_DEN_genTracks/abseta_pt_value
#	 NUM_MediumID_DEN_genTracks/abseta_pt_error
#available evaluator keys (muon_TRIG):
#	 IsoMu27_PtEtaBins/pt_abseta_ratio_value
#	 IsoMu27_PtEtaBins/pt_abseta_ratio_error
#	 Mu50_PtEtaBins/pt_abseta_ratio_value
#	 Mu50_PtEtaBins/pt_abseta_ratio_error

lepton_sf_dict = {"elec_RECO":["egammaEffi_txt_EGM2D_runBCDEF_passingRECO.root","EGamma_SF2D","EGamma_SF2D_error",0],#0 for eta-pt, 1 for abseta-pt, 2 for pt-abseta
                  "elec_ID":["egammaEffi_txt_EGM2D_runBCDEF_passingID.root","EGamma_SF2D","EGamma_SF2D_error",0],
                  "elec_TRIG32":["egammaEffi_txt_runBCDEF_passingEle32.root","EGamma_SF2D","EGamma_SF2D_error",0],
                  "elec_TRIG115":["egammaEffi_txt_runBCDEF_passingEle115.json","Ele115_PtEtaBins/abseta_pt_SF_value","Ele115_PtEtaBins/abseta_pt_SF_error",1],
                  "muon_ISO":["muonEff_RunBCDEF_SF_ISO.json","NUM_LooseRelIso_DEN_MediumID/abseta_pt_value","NUM_LooseRelIso_DEN_MediumID/abseta_pt_error",1],
                  "muon_ID":["muonEff_RunBCDEF_SF_ID.json","NUM_MediumID_DEN_genTracks/abseta_pt_value","NUM_MediumID_DEN_genTracks/abseta_pt_error",1],
                  "muon_TRIG27":["muonEff_RunBCDEF_SF_Trig_Nov17Nov2017.json","IsoMu27_PtEtaBins/pt_abseta_ratio_value","IsoMu27_PtEtaBins/pt_abseta_ratio_error",2],
                  "muon_TRIG50":["muonEff_RunBCDEF_SF_Trig_Nov17Nov2017.json","Mu50_PtEtaBins/pt_abseta_ratio_value","Mu50_PtEtaBins/pt_abseta_ratio_error",2],
}
ext = extractor()

for sfname, sfopts in lepton_sf_dict.items():
    ext.add_weight_sets(["%s_value %s %s"%(sfname,sfopts[1],os.path.join(os.path.dirname(__file__), 'data' , sfopts[0]))])
    ext.add_weight_sets(["%s_error %s %s"%(sfname,sfopts[2],os.path.join(os.path.dirname(__file__), 'data' , sfopts[0]))])
ext.finalize()
lepsf_evaluator = ext.make_evaluator()

puW2017_nonUL = np.array([0.183454, 3.93313, 3.47111, 2.4924, 1.62495, 1.5151, 1.28791, 1.27825, 0.615845, 1.45208, 1.49768, 1.48747, 1.33116, 1.1645, 1.07901, 1.05437, 1.08066, 1.12907, 1.16584, 1.18936, 1.21284, 1.23849, 1.25967, 1.27099, 1.2727, 1.2713, 1.27087, 1.26652, 1.27449, 1.25163, 1.22131, 1.16954, 1.10903, 1.03816, 0.969432, 0.91188, 0.86681, 0.834505, 0.788379, 0.750796, 0.759152, 0.793175, 0.858347, 0.958795, 1.09432, 1.25564, 1.41965, 1.49593, 1.53054, 1.4622, 1.33676, 1.15439, 0.950556, 0.749702, 0.569789, 0.410747, 0.290007, 0.198655, 0.137644, 0.096639, 0.0692314, 0.0508587, 0.038443, 0.0299595, 0.0240949, 0.01712, 0.0124798, 0.0107683, 0.0095972, 0.00881241, 0.00830987, 0.00801759, 0.00788476, 0.0078747, 0.00633731, 0.00533369, 0.00544371, 0.00558365, 0.00574411, 0.00591599, 0.00609007, 0.00625711, 0.00640816, 0.00509462, 0.00422803, 0.00425915, 0.00426558, 0.00424666, 0.00420325, 0.00413738, 0.00405193, 0.00395019, 0.00294303, 0.00229603, 0.00220003, 0.00210138, 0.00200166, 0.00190214, 0.0018038]) #for non-UL
notUL2017 = [
    "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8",
    "TTToHadronic_TuneCP5_13TeV-powheg-pythia8",
    "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8",
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
    if (year == '2017' or year == '2018') and dataset in notUL2017:
        weights.add('pileup_weight',
            puW2017_nonUL[np.clip(nPU,0,len(puW2017_nonUL)-1)]
        )
    elif year == '2017' and dataset in compiled['2017_pileupweight_dataset']:
        weights.add(
            'pileup_weight',
            compiled['2017_pileupweight_dataset'][dataset](nPU),
            compiled['2017_pileupweight_dataset_puUp'][dataset](nPU),
            compiled['2017_pileupweight_dataset_puDown'][dataset](nPU),
        )
    else:
        weights.add(
            'pileup_weight',
            compiled[f'{year}_pileupweight'](nPU),
            compiled[f'{year}_pileupweight_puUp'](nPU),
            compiled[f'{year}_pileupweight_puDown'](nPU),
        )


def add_VJets_NLOkFactor(weights, genBosonPt, year, dataset):
    if (year == '2017' or year == '2018') and 'ZJetsToQQ_HT' in dataset:
        nlo_over_lo_qcd = compiled['2017_Z_nlo_qcd'](genBosonPt)
        nlo_over_lo_ewk = compiled['Z_nlo_over_lo_ewk'](genBosonPt)
    elif (year == '2017' or year == '2018') and 'WJetsToQQ_HT' in dataset:
        nlo_over_lo_qcd = compiled['2017_W_nlo_qcd'](genBosonPt)
        nlo_over_lo_ewk = compiled['W_nlo_over_lo_ewk'](genBosonPt)
    elif year == '2016' and 'DYJetsToQQ' in dataset:
        nlo_over_lo_qcd = compiled['2016_Z_nlo_qcd'](genBosonPt)
        nlo_over_lo_ewk = compiled['Z_nlo_over_lo_ewk'](genBosonPt)
    elif year == '2016' and 'WJetsToQQ' in dataset:
        nlo_over_lo_qcd = compiled['2016_W_nlo_qcd'](genBosonPt)
        nlo_over_lo_ewk = compiled['W_nlo_over_lo_ewk'](genBosonPt)
    elif 'DYJetsToLL_Pt' in dataset:
        nlo_over_lo_qcd = np.ones_like(genBosonPt.flatten())
        nlo_over_lo_ewk = Vpt_corr_value[np.digitize(np.clip(genBosonPt.flatten(),Vpt_corr_bins[0], Vpt_corr_bins[-1]), Vpt_corr_bins)-1]
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

#    0        eta - pt
#    1        abseta - pt
#    2        pt - abseta
def add_LeptonSFs(weights, lep_pt, lep_eta, year, match):
    lep_pt  = lep_pt.pad(1, clip=True).fillna(0).flatten()
    lep_eta = lep_eta.pad(1, clip=True).fillna(0).flatten()
    for sf in lepton_sf_dict:
        if match in sf:
            if lepton_sf_dict[sf][3]==0:
                nom = lepsf_evaluator['%s_value'%sf](lep_eta,lep_pt)
                err = lepsf_evaluator['%s_value'%sf](lep_eta,lep_pt)
            elif lepton_sf_dict[sf][3]==1:
                nom = lepsf_evaluator['%s_value'%sf](np.abs(lep_eta),lep_pt)
                err = lepsf_evaluator['%s_value'%sf](np.abs(lep_eta),lep_pt)
            elif lepton_sf_dict[sf][3]==2:
                nom = lepsf_evaluator['%s_value'%sf](lep_pt,np.abs(lep_eta))
                err = lepsf_evaluator['%s_value'%sf](lep_pt,np.abs(lep_eta))
            else: 
                print('invalid type ordering for lepton SF %s'%sf)
                return
            wname = sf
            if "TRIG27" in sf:
                nom[lep_pt>55.] = 1.
                err[lep_pt>55.] = 0.
            if "TRIG50" in sf:
                nom[lep_pt<55.] = 1.
                err[lep_pt<55.] = 0.
            if "TRIG32" in sf:
                nom[lep_pt>120.] = 1.
                err[lep_pt>120.] = 0.
            if "TRIG115" in sf:
                nom[lep_pt<120.] = 1.
                err[lep_pt<120.] = 0.
            weights.add(sf, nom, nom+err, nom-err)

def add_TopPtReweighting(weights, topPt, year, dataset):
#$SF(p_T)=e^{0.0615-0.0005\cdot p_T}$ for data/POWHEG+Pythia8
    if 'TT' in dataset:
        toppt_weight1 = np.exp(0.0615-0.0005*np.clip(topPt[:,0],0.,500.))
        toppt_weight2 = np.exp(0.0615-0.0005*np.clip(topPt[:,1],0.,500.))
    else:
        toppt_weight1 = np.ones_like(topPt[:,0])
        toppt_weight2 = np.ones_like(topPt[:,1])
    weights.add('TopPtReweight', np.sqrt(toppt_weight1 * toppt_weight2), np.ones_like(toppt_weight1), np.sqrt(toppt_weight1 * toppt_weight2))
