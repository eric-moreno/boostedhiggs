import os
import numpy as np
import awkward as ak
import gzip
import pickle
import importlib.resources
import correctionlib
from coffea.lookup_tools.lookup_base import lookup_base
from coffea import lookup_tools
from coffea import util

with importlib.resources.path("boostedhiggs.data", "corrections.pkl.gz") as path:
    with gzip.open(path) as fin:
        compiled = pickle.load(fin)

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

def add_pdf_weight(weights, pdf_weights):
    nom = np.ones(len(weights.weight()))
    up = np.ones(len(weights.weight()))
    down = np.ones(len(weights.weight()))

    # NNPDF31_nnlo_hessian_pdfas
    # https://lhapdfsets.web.cern.ch/current/NNPDF31_nnlo_hessian_pdfas/NNPDF31_nnlo_hessian_pdfas.info
    if pdf_weights is not None and "306000 - 306102" in pdf_weights.__doc__:
        # Hessian PDF weights
        # Eq. 21 of https://arxiv.org/pdf/1510.03865v1.pdf
        arg = pdf_weights[:, 1:-2] - np.ones((len(weights.weight()), 100))
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

# TODO: Update to UL
# Option: 0 for eta-pt, 1 for abseta-pt, 2 for pt-abseta
lepton_sf_dict = {
    "elec_RECO": 0,
    "elec_ID": 0, 
    "elec_TRIG32": 0,
    "elec_TRIG115": 1,
    "muon_ISO": 1,
    "muon_ID": 1,
    "muon_TRIG27": 2,
    "muon_TRIG50": 2,
}

def add_leptonSFs(weights, lepton, year, match):
    for sf in lepton_sf_dict:
        sfoption = lepton_sf_dict[sf]
        lep_pt = np.array(ak.fill_none(lepton.pt, 0.))
        lep_eta = np.array(ak.fill_none(lepton.eta, 0.))
        lep_abseta = np.array(ak.fill_none(abs(lepton.eta), 0.))
        if match in sf:
            if sfoption==0:
                nom = compiled['%s_value'%sf](lep_eta,lep_pt)
                err = compiled['%s_error'%sf](lep_eta,lep_pt)
            elif sfoption==1:
                nom = compiled['%s_value'%sf](np.abs(lep_eta),lep_pt)
                err = compiled['%s_error'%sf](np.abs(lep_eta),lep_pt)
            elif sfoption==2:
                nom = compiled['%s_value'%sf](lep_pt,np.abs(lep_eta))
                err = compiled['%s_error'%sf](lep_pt,np.abs(lep_eta))
            else: 
                print('Error: Invalid type ordering for lepton SF %s'%sf)
                return
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

def add_TopPtReweighting(weights, topPt, year, dataset):
#$SF(p_T)=e^{0.0615-0.0005\cdot p_T}$ for data/POWHEG+Pythia8
    if 'TT' in dataset:
        toppt_weight1 = np.exp(0.0615-0.0005*np.clip(topPt[:,0],0.,500.))
        toppt_weight2 = np.exp(0.0615-0.0005*np.clip(topPt[:,1],0.,500.))
    else:
        toppt_weight1 = np.ones_like(topPt[:,0])
        toppt_weight2 = np.ones_like(topPt[:,1])
    weights.add('TopPtReweight', np.sqrt(toppt_weight1 * toppt_weight2), np.ones_like(toppt_weight1), np.sqrt(toppt_weight1 * toppt_weight2))
