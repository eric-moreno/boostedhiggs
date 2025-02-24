import os
import numpy as np
import awkward as ak
import gzip
import pickle
import importlib.resources
import correctionlib
import json
from coffea.lookup_tools.lookup_base import lookup_base
from coffea import lookup_tools
from coffea import util

with importlib.resources.path("boostedhiggs.data", "corrections.pkl.gz") as path:
    with gzip.open(path) as fin:
        compiled = pickle.load(fin)

compiled['2016preVFP_pileupweight']._values = compiled['2016preVFP_pileupweight']._values
compiled['2016preVFP_pileupweight_puUp']._values = compiled['2016preVFP_pileupweight_puUp']._values
compiled['2016preVFP_pileupweight_puDown']._values = compiled['2016preVFP_pileupweight_puDown']._values
compiled['2016postVFP_pileupweight']._values = compiled['2016postVFP_pileupweight']._values
compiled['2016postVFP_pileupweight_puUp']._values = compiled['2016postVFP_pileupweight_puUp']._values
compiled['2016postVFP_pileupweight_puDown']._values = compiled['2016postVFP_pileupweight_puDown']._values
compiled['2017_pileupweight']._values = compiled['2017_pileupweight']._values
compiled['2017_pileupweight_puUp']._values = compiled['2017_pileupweight_puUp']._values
compiled['2017_pileupweight_puDown']._values = compiled['2017_pileupweight_puDown']._values
compiled['2018_pileupweight']._values = compiled['2018_pileupweight']._values
compiled['2018_pileupweight_puUp']._values = compiled['2018_pileupweight_puUp']._values
compiled['2018_pileupweight_puDown']._values = compiled['2018_pileupweight_puDown']._values

class SoftDropWeight(lookup_base):
    def _evaluate(self, pt, eta):
        gpar = np.array([1.00626, -1.06161, 0.0799900, 1.20454])
        cpar = np.array([1.09302, -0.000150068, 3.44866e-07, -2.68100e-10, 8.67440e-14, -1.00114e-17])
        fpar = np.array([1.27212, -0.000571640, 8.37289e-07, -5.20433e-10, 1.45375e-13, -1.50389e-17])
        genw = gpar[0] + gpar[1]*np.power(pt*gpar[2], -gpar[3])
        ptpow = np.power.outer(pt, np.arange(cpar.size))
        cenweight = np.dot(ptpow, cpar)
        forweight = np.dot(ptpow, fpar)
        weight = np.where(np.abs(eta) < 1.3, cenweight, forweight)
        return genw*weight

_softdrop_weight = SoftDropWeight()

def corrected_msoftdrop(fatjets):
    sf = _softdrop_weight(fatjets.pt, fatjets.eta)
    sf = np.maximum(1e-5, sf)
    dazsle_msd = (fatjets.subjets * (1 - fatjets.subjets.rawFactor)).sum()
    return dazsle_msd.mass * sf

def add_pileup_weight(weights, nPU, year='2017'):
    weights.add(
        'pileup_weight',
        compiled[f'{year}_pileupweight'](nPU),
        compiled[f'{year}_pileupweight_puUp'](nPU),
        compiled[f'{year}_pileupweight_puDown'](nPU),
    )

def add_ps_weight(weights, ps_weights):
    nom = np.ones(len(weights.weight()))
    up_isr = np.ones(len(weights.weight()))
    down_isr = np.ones(len(weights.weight()))
    up_fsr = np.ones(len(weights.weight()))
    down_fsr = np.ones(len(weights.weight()))

    if ps_weights is not None:
        if len(ps_weights[0]) == 4:
            up_isr = ps_weights[:, 0]
            down_isr = ps_weights[:, 2]
            up_fsr = ps_weights[:, 1]
            down_fsr = ps_weights[:, 3]
        else:
            warnings.warn(f"PS weight vector has length {len(ps_weights[0])}")
    weights.add('UEPS_ISR', nom, up_isr, down_isr)
    weights.add('UEPS_FSR', nom, up_fsr, down_fsr)

def build_lumimask(filename):
    from coffea.lumi_tools import LumiMask
    with importlib.resources.path("boostedhiggs.data", filename) as path:
        return LumiMask(path)

lumiMasks = {
    "2016": build_lumimask("Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"),
    "2017": build_lumimask("Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"),
    "2018": build_lumimask("Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"),
}

basedir = 'boostedhiggs/data/'
mutriglist = {
    '2016preVFP':{
        'TRIGNOISO':'NUM_Mu50_or_TkMu50_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose_abseta_pt',
        'TRIGISO':'NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt',
    },
    '2016postVFP':{
        'TRIGNOISO':'NUM_Mu50_or_TkMu50_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose_abseta_pt',
        'TRIGISO':'NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt',
    },
    '2017':{
        'TRIGNOISO':'NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose_abseta_pt',
        'TRIGISO':'NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt',
    },
    '2018':{
        'TRIGNOISO':'NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose_abseta_pt',
        'TRIGISO':'NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt',
    },
}
eltriglist = {
    '2016preVFP':{
        'TRIGNOISO':'EGamma_SF2D',
        'TRIGISO':'EGamma_SF2D',
    },
    '2016postVFP':{
        'TRIGNOISO':'EGamma_SF2D',
        'TRIGISO':'EGamma_SF2D',
    },
    '2017':{
        'TRIGNOISO':'EGamma_SF2D',
        'TRIGISO':'EGamma_SF2D',
    },
    '2018':{
        'TRIGNOISO':'EGamma_SF2D',
        'TRIGISO':'EGamma_SF2D',
    },
}

ext = lookup_tools.extractor()
for year in ['2016preVFP','2016postVFP','2017','2018']:
    ext.add_weight_sets([f'muon_ID_{year} NUM_MediumPromptID_DEN_TrackerMuons_abseta_pt {basedir}Efficiencies_muon_generalTracks_Z_Run{year}_UL_ID.root'])
    ext.add_weight_sets([f'muon_ISO_{year} NUM_LooseRelIso_DEN_MediumPromptID_abseta_pt {basedir}Efficiencies_muon_generalTracks_Z_Run{year}_UL_ISO.root'])
    for trigopt in mutriglist[year]:
        trigname = mutriglist[year][trigopt]
        ext.add_weight_sets([f'muon_{trigopt}_{year} {trigname} {basedir}Efficiencies_muon_generalTracks_Z_Run{year}_UL_SingleMuonTriggers.root'])
    ext.add_weight_sets([f'elec_RECO_{year} EGamma_SF2D {basedir}egammaEffi_ptAbove20_txt_EGM2D_UL{year}.root'])
    ext.add_weight_sets([f'elec_ID_{year} EGamma_SF2D {basedir}egammaEffi_txt_Ele_wp90noiso_UL{year}_EGM2D.root'])
    for trigopt in eltriglist[year]:
        trigname = eltriglist[year][trigopt]
        ext.add_weight_sets([f'elec_{trigopt}_{year} {trigname} {basedir}egammaEffi_txt_trigger_EGM2D_UL{year}.root'])
ext.finalize()
lepsf_evaluator = ext.make_evaluator()
lepsf_keys = lepsf_evaluator.keys()

def add_LeptonSFs(weights, lepton, year, match):
    for sf in lepsf_keys:
        if year not in sf:
            continue
        if match not in sf:
            continue
        lep_pt = np.array(ak.fill_none(lepton.pt, 0.))
        lep_eta = np.array(ak.fill_none(lepton.eta, 0.))
        if match in sf:
            if 'muon' in match:
                nom = lepsf_evaluator[sf](np.abs(lep_eta),lep_pt)
            elif 'elec' in match:
                nom = lepsf_evaluator[sf](lep_eta,lep_pt)
            else: 
                print('Error: Invalid type ordering for lepton SF %s'%sf)
                return
            if "_TRIGISO_" in sf:
                if 'muon' in match:
                    nom[lep_pt>55.] = 1.
                else:
                    nom[lep_pt>120.] = 1.
            if "_TRIGNOISO_" in sf:
                if 'muon' in match:
                    nom[lep_pt<55.] = 1.
                else:
                    nom[lep_pt<120.] = 1.
            if "_ISO_" in sf:
                if 'muon' in match:
                    nom[lep_pt>55.] = 1.
                else:
                    nom[lep_pt>120.] = 1.
            weights.add(sf, nom)

def is_overlap(events,dataset,triggers,year):
    dataset_ordering = {
        '2016':['SingleMuon','SingleElectron','MET','JetHT'],
        '2017':['SingleMuon','SingleElectron','MET','JetHT'],
        '2018':['SingleMuon','EGamma','MET','JetHT']
    }
    pd_to_trig = {
        'SingleMuon': ['Mu50',
                       'Mu55',
                       'Mu15_IsoVVVL_PFHT600',
                       'Mu15_IsoVVVL_PFHT450_PFMET50',
                       ],
        'SingleElectron': ['Ele50_CaloIdVT_GsfTrkIdT_PFJet165',
                           'Ele115_CaloIdVT_GsfTrkIdT',
                           'Ele15_IsoVVVL_PFHT600',
                           'Ele35_WPTight_Gsf',
                           'Ele15_IsoVVVL_PFHT450_PFMET50',
                       ],
        'JetHT': ['PFHT800',
                  'PFHT900',
                  'AK8PFJet360_TrimMass30',
                  'AK8PFHT700_TrimR0p1PT0p03Mass50',
                  'PFHT650_WideJetMJJ950DEtaJJ1p5',
                  'PFHT650_WideJetMJJ900DEtaJJ1p5',
                  'PFJet450',
                  'PFHT1050',
                  'PFJet500',
                  'AK8PFJet400_TrimMass30',
                  'AK8PFJet420_TrimMass30',
                  'AK8PFHT800_TrimMass50'
              ],
        'MET': ['PFMETNoMu120_PFMHTNoMu120_IDTight',
                'PFMETNoMu110_PFMHTNoMu110_IDTight',
            ],
    }
    
    overlap = np.ones(len(events), dtype='bool')
    for p in dataset_ordering[year]:
        if dataset.startswith(p):
            pass_pd = np.zeros(len(events), dtype='bool')
            for t in pd_to_trig[p]:
                if t in events.HLT.fields:
                    pass_pd = pass_pd | events.HLT[t]
            overlap = overlap & pass_pd
            break
        else:
            for t in pd_to_trig[p]:
                if t in events.HLT.fields:
                    overlap = overlap & np.logical_not(events.HLT[t])
    return overlap

met_trigger_file = open('boostedhiggs/data/mettrig_sfs.json','r')
met_trigger_sfs = json.load(met_trigger_file)
for key in met_trigger_sfs:
    if key=='binning':
        met_trigger_sfs[key] = np.array(met_trigger_sfs[key])
    else:
        met_trigger_sfs[key]['value'] = np.array(met_trigger_sfs[key]['value'])
        met_trigger_sfs[key]['error'] = np.array(met_trigger_sfs[key]['error'])

def add_METSFs(weights, met, year):
    _met_nom = met_trigger_sfs[year]['value'][np.digitize(np.clip(ak.to_numpy(met).flatten(), met_trigger_sfs['binning'][0], met_trigger_sfs['binning'][-1]-1.), met_trigger_sfs['binning'])-1]
    _met_err = met_trigger_sfs[year]['error'][np.digitize(np.clip(ak.to_numpy(met).flatten(), met_trigger_sfs['binning'][0], met_trigger_sfs['binning'][-1]-1.), met_trigger_sfs['binning'])-1]
    weights.add('MET_TRIG', _met_nom, _met_nom+_met_err, _met_nom-_met_err)
    
def add_JetTriggerSF(weights, jets, year):
    cset = correctionlib.CorrectionSet.from_file(f"boostedhiggs/data/fatjet_triggerSF{year}.json")[f'fatjet_triggerSF{year}']
    _msd = ak.to_numpy(jets.msoftdrop, allow_missing=True)
    _pt = ak.to_numpy(jets.pt, allow_missing=True)
    _sf = cset.evaluate("nominal", _pt, _msd)
    _sfup = cset.evaluate("stat_up",  _pt, _msd)
    _sfdn = cset.evaluate("stat_dn", _pt, _msd)
    weights.add('HAD_TRIG', _sf, _sfup, _sfdn)

def add_pdf_weight(weights, pdf_weights):
    nweights = len(weights.weight())
    nom = np.ones(nweights)
    up = np.ones(nweights)
    down = np.ones(nweights)
    docstring = pdf_weights.__doc__

    # NNPDF31_nnlo_hessian_pdfas
    # https://lhapdfsets.web.cern.ch/current/NNPDF31_nnlo_hessian_pdfas/NNPDF31_nnlo_hessian_pdfas.info
    if "306000 - 306102" in docstring:
        # Hessian PDF weights
        # Eq. 21 of https://arxiv.org/pdf/1510.03865v1.pdf
        arg = pdf_weights[:, 1:-2] - np.ones((nweights, 100))
        summed = ak.sum(np.square(arg), axis=1)
        pdf_unc = np.sqrt((1. / 99.) * summed)
        weights.add('PDF_weight', nom, pdf_unc + nom)

        # alpha_S weights
        # Eq. 27 of same ref
        as_unc = 0.5 * (pdf_weights[:, 102] - pdf_weights[:, 101])
        weights.add('aS_weight', nom, as_unc + nom)

        # PDF + alpha_S weights
        # Eq. 28 of same ref
        pdfas_unc = np.sqrt(np.square(pdf_unc) + np.square(as_unc))
        weights.add('PDFaS_weight', nom, pdfas_unc + nom)

    else:
        weights.add('aS_weight', nom, up, down)
        weights.add('PDF_weight', nom, up, down)
        weights.add('PDFaS_weight', nom, up, down)

# 7-point scale variations
def add_scalevar_7pt(weights,var_weights):
    docstring = var_weights.__doc__
    nweights = len(weights.weight())

    nom   = np.ones(nweights)
    up    = np.ones(nweights)
    down  = np.ones(nweights)
 
    if len(var_weights) > 0:
        if len(var_weights[0]) == 9: 
            up = np.maximum.reduce([var_weights[:,0],var_weights[:,1],var_weights[:,3],var_weights[:,5],var_weights[:,7],var_weights[:,8]])
            down = np.minimum.reduce([var_weights[:,0],var_weights[:,1],var_weights[:,3],var_weights[:,5],var_weights[:,7],var_weights[:,8]])
        elif len(var_weights[0]) == 8:
            up = np.maximum.reduce([var_weights[:,0],var_weights[:,1],var_weights[:,3],var_weights[:,4],var_weights[:,6],var_weights[:,7]])
            down = np.minimum.reduce([var_weights[:,0],var_weights[:,1],var_weights[:,3],var_weights[:,4],var_weights[:,6],var_weights[:,7]])
        elif len(var_weights[0]) > 1:
            print("Scale variation vector has length ", len(var_weights[0]))
    weights.add('scalevar_7pt', nom, up, down)

# 3-point scale variations
def add_scalevar_3pt(weights,var_weights):
    docstring = var_weights.__doc__
    
    nweights = len(weights.weight())

    nom   = np.ones(nweights)
    up    = np.ones(nweights)
    down  = np.ones(nweights)

    if len(var_weights) > 0:
        if len(var_weights[0]) == 9:
            up = np.maximum(var_weights[:,0], var_weights[:,8])
            down = np.minimum(var_weights[:,0], var_weights[:,8])
        elif len(var_weights[0]) == 8:
            up = np.maximum(var_weights[:,0], var_weights[:,7])
            down = np.minimum(var_weights[:,0], var_weights[:,7])
        elif len(var_weights[0]) > 1:
            print("Scale variation vector has length ", len(var_weights[0]))

    weights.add('scalevar_3pt', nom, up, down)

import cloudpickle
with importlib.resources.path("boostedhiggs.data", "jec_compiled.pkl.gz") as path:
    with gzip.open(path) as fin:
        jmestuff = cloudpickle.load(fin)

jet_factory = jmestuff["jet_factory"]
fatjet_factory = jmestuff["fatjet_factory"]
met_factory = jmestuff["met_factory"]


def add_jec_variables(jets, event_rho):
    jets["pt_raw"] = (1 - jets.rawFactor)*jets.pt
    jets["mass_raw"] = (1 - jets.rawFactor)*jets.mass
    jets["pt_gen"] = ak.values_astype(ak.fill_none(jets.matched_gen.pt, 0), np.float32)
    jets["event_rho"] = ak.broadcast_arrays(event_rho, jets.pt)[0]
    return jets

#from https://www.research-collection.ethz.ch/bitstream/handle/20.500.11850/240924/1/s10052-017-5389-1.pdf
Vpt_corr_bins = np.array([
    30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 950.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1600.0, 1800.0, 2000.0, 2200.0, 2400.0, 2600.0, 2800.0, 3000.0, 6500.0
])
Vpt_corr_value = np.array([
    1.003846859, 1.005303426, 0.9947071550000002, 0.98928453, 0.9850576999999998, 0.9813048799999999, 0.9778652800000001, 0.97445919, 0.9710953100000002, 0.9677317199999999, 0.9642441799999999, 0.9608947800000001, 0.9525394399999999, 0.9365518699999998, 0.9219756100000001, 0.9087814499999999, 0.8968980000000001, 0.8859843000000001, 0.8760160999999999, 0.8671511999999999, 0.8583182999999999, 0.8507631, 0.8431063999999999, 0.8362730999999999, 0.8303883000000001, 0.8242617999999999, 0.8194368000000001, 0.8138635, 0.8077816999999999, 0.8017448, 0.7931413999999999, 0.7852697, 0.7775183999999999, 0.7692830000000002, 0.7558978, 0.7443137, 0.7334527, 0.7233202, 0.7140055999999999, 0.7045699000000001, 0.691076, 0.6890300000000001
])

def add_VJets_NLOkFactor(weights, genBosonPt, year, dataset):
    #if (year == '2017' or year == '2018') and 'ZJetsToQQ_HT' in dataset:
    #    nlo_over_lo_qcd = compiled['2017_Z_nlo_qcd'](genBosonPt)
    #    nlo_over_lo_ewk = compiled['Z_nlo_over_lo_ewk'](genBosonPt)
    #elif (year == '2017' or year == '2018') and 'WJetsToQQ_HT' in dataset:
    #    nlo_over_lo_qcd = compiled['2017_W_nlo_qcd'](genBosonPt)
    #    nlo_over_lo_ewk = compiled['W_nlo_over_lo_ewk'](genBosonPt)
    #elif year == '2016' and 'DYJetsToQQ' in dataset:
    #    nlo_over_lo_qcd = compiled['2016_Z_nlo_qcd'](genBosonPt)
    #    nlo_over_lo_ewk = compiled['Z_nlo_over_lo_ewk'](genBosonPt)
    #elif year == '2016' and 'WJetsToQQ' in dataset:
    #    nlo_over_lo_qcd = compiled['2016_W_nlo_qcd'](genBosonPt)
    #    nlo_over_lo_ewk = compiled['W_nlo_over_lo_ewk'](genBosonPt)
    if 'DYJetsToLL_Pt' in dataset:
        nlo_over_lo_qcd = np.ones_like(ak.to_numpy(genBosonPt).flatten())
        nlo_over_lo_ewk = Vpt_corr_value[np.digitize(np.clip(ak.to_numpy(genBosonPt).flatten(),Vpt_corr_bins[0], Vpt_corr_bins[-1]), Vpt_corr_bins)-1]
    else:
        return
    weights.add('VJets_NLOkFactor', nlo_over_lo_qcd * nlo_over_lo_ewk)

with importlib.resources.path("boostedhiggs.data", "ULvjets_corrections.json") as filename:
    vjets_kfactors = correctionlib.CorrectionSet.from_file(str(filename))

#taking from Jennet
def add_VJets_kFactors(weights, genpart, dataset):
    """Revised version of add_VJets_NLOkFactor, for both NLO EW and ~NNLO QCD"""
    def get_vpt(check_offshell=False):
        """Only the leptonic samples have no resonance in the decay tree, and only
        when M is beyond the configured Breit-Wigner cutoff (usually 15*width)
        """
        boson = ak.firsts(genpart[
            ((genpart.pdgId == 23)|(abs(genpart.pdgId) == 24))
            & genpart.hasFlags(["fromHardProcess", "isLastCopy"])
        ])
        if check_offshell:
            offshell = genpart[
                genpart.hasFlags(["fromHardProcess", "isLastCopy"])
                & ak.is_none(boson)
                & (abs(genpart.pdgId) >= 11) & (abs(genpart.pdgId) <= 16)
            ].sum()
            return ak.where(ak.is_none(boson.pt), offshell.pt, boson.pt)
        return np.array(ak.fill_none(boson.pt, 0.))

    common_systs = [
        "d1K_NLO",
        "d2K_NLO",
        "d3K_NLO",
        "d1kappa_EW",
    ]
    zsysts = common_systs + [
        "Z_d2kappa_EW",
        "Z_d3kappa_EW",
    ]
    znlosysts = [
        "d1kappa_EW",
        "Z_d2kappa_EW",
        "Z_d3kappa_EW",
    ]
    wsysts = common_systs + [
        "W_d2kappa_EW",
        "W_d3kappa_EW",
    ]

    def add_systs(systlist, qcdcorr, ewkcorr, vpt):
        ewknom = ewkcorr.evaluate("nominal", vpt)
        weights.add("vjets_nominal", qcdcorr * ewknom if qcdcorr is not None else ewknom)
        ones = np.ones_like(vpt)
        for syst in systlist:
            weights.add(syst, ones, ewkcorr.evaluate(syst + "_up", vpt) / ewknom, ewkcorr.evaluate(syst + "_down", vpt) / ewknom)

    if "ZJetsToQQ_HT" in dataset or "DYJetsToLL_M-50" in dataset:
        vpt = get_vpt()
        #qcdcorr = vjets_kfactors["ULZ_MLMtoFXFX"].evaluate(vpt)
        qcdcorr = vjets_kfactors["ULZ_MLMtoNNLOQCD"].evaluate(vpt)
        ewkcorr = vjets_kfactors["Z_FixedOrderComponent"]
        add_systs(zsysts, qcdcorr, ewkcorr, vpt)
    elif "WJetsToQQ_HT" in dataset or "WJetsToLNu_HT" in dataset:
        vpt = get_vpt()
        #qcdcorr = vjets_kfactors["ULW_MLMtoFXFX"].evaluate(vpt)
        qcdcorr = vjets_kfactors["ULW_MLMtoNNLOQCD"].evaluate(vpt)
        ewkcorr = vjets_kfactors["W_FixedOrderComponent"]
        add_systs(wsysts, qcdcorr, ewkcorr, vpt)
    elif "DYJetsToLL_Pt" in dataset:
        vpt = get_vpt()
        nnloqcd = vjets_kfactors["ULZ_MLMtoNNLOQCD"].evaluate(vpt)/vjets_kfactors["ULZ_MLMtoFXFX"].evaluate(vpt)
        ewkcorr = vjets_kfactors["Z_FixedOrderComponent"]
        add_systs(znlosysts, nnloqcd, ewkcorr, vpt)

def add_TopPtReweighting(weights, topPt, year, dataset):
#$SF(p_T)=e^{0.0615-0.0005\cdot p_T}$ for data/POWHEG+Pythia8
    if 'TT' in dataset:
        toppt_weight1 = np.exp(0.0615-0.0005*np.clip(topPt[:,0],0.,500.))
        toppt_weight2 = np.exp(0.0615-0.0005*np.clip(topPt[:,1],0.,500.))
        weights.add('TopPtReweight', np.sqrt(toppt_weight1 * toppt_weight2), np.ones_like(toppt_weight1), np.sqrt(toppt_weight1 * toppt_weight2))
