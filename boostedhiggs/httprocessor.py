from functools import partial
import numpy as np
import awkward as ak
import hist as hist2
import json
from copy import deepcopy

from coffea import processor, hist
from coffea.nanoevents.methods import candidate, vector
from coffea.analysis_tools import Weights, PackedSelection

from .corrections import (
    corrected_msoftdrop,
#    n2ddt_shift,
    add_pileup_weight,
#    add_VJets_NLOkFactor,
#    add_TopPtReweighting,
#    add_jetTriggerWeight,
#    add_TriggerWeight,
#    add_LeptonSFs,
)

from .common import (
    getBosons,
    getParticles,
    matchedBosonFlavor,
    matchedBosonFlavorLep,
    getHTauTauDecayInfo,
    isOverlap,
)

import logging
logger = logging.getLogger(__name__)

# function to normalize arrays after a cut or selection
def normalize(val, cut=None):
    if cut is None:
        ar = ak.to_numpy(ak.fill_none(val, np.nan))
        return ar
    else:
        ar = ak.to_numpy(ak.fill_none(val[cut], np.nan))
        return ar

class HttProcessor(processor.ProcessorABC):
    def __init__(self, year="2017", jet_arbitration='met'):
        self._year = year
        self._jet_arbitration = jet_arbitration
        
        self._triggers = {
            '2016': {
                'e': [
                    "Ele50_CaloIdVT_GsfTrkIdT_PFJet165",
                    "Ele115_CaloIdVT_GsfTrkIdT",
                    "Ele15_IsoVVVL_PFHT600",
                    'PFHT800',
                    'PFHT900',
                    'AK8PFJet360_TrimMass30',
                    'AK8PFHT700_TrimR0p1PT0p03Mass50',
                    'PFHT650_WideJetMJJ950DEtaJJ1p5',
                    'PFHT650_WideJetMJJ900DEtaJJ1p5',
                    'PFJet450',

                ],
                'mu': [
                    "Mu50",
                    "Mu55",
                    "Mu15_IsoVVVL_PFHT600",
                    'PFHT800',
                    'PFHT900',
                    'AK8PFJet360_TrimMass30',
                    'AK8PFHT700_TrimR0p1PT0p03Mass50',
                    'PFHT650_WideJetMJJ950DEtaJJ1p5',
                    'PFHT650_WideJetMJJ900DEtaJJ1p5',
                    'PFJet450',
                ],
                'had': [
                    'PFHT800',
                    'PFHT900',
                    'AK8PFJet360_TrimMass30',
                    'AK8PFHT700_TrimR0p1PT0p03Mass50',
                    'PFHT650_WideJetMJJ950DEtaJJ1p5',
                    'PFHT650_WideJetMJJ900DEtaJJ1p5',
                    'PFJet450',
                    'DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg',
                    'DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg',
                    'DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg',
                    'DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg',
                ],
            },
            '2017': {
                'e': [
                    'Ele35_WPTight_Gsf',
                    'Ele115_CaloIdVT_GsfTrkIdT',
                ],
                'mu': [
                    "Mu50",
                    "IsoMu27",
                ],
                'had': [
                    'PFHT1050',
                    'AK8PFJet400_TrimMass30',
                    'AK8PFJet420_TrimMass30',
                    'AK8PFHT800_TrimMass50',
                    'PFJet500',
                    'AK8PFJet500',
                ],
            },
            "2018": {
                'e': [
                    'Ele35_WPTight_Gsf',
                    'Ele115_CaloIdVT_GsfTrkIdT',
                ],
                'mu': [
                    "Mu50",
                    "IsoMu27",
                ],
                'had': [
                    'PFHT1050',
                    'AK8PFJet400_TrimMass30',
                    'AK8PFJet420_TrimMass30',
                    'AK8PFHT800_TrimMass50',
                    'PFJet500',
                    'AK8PFJet500',
                ],
            }
        }[year]

        self._metFilters = {
            '2016': [
                "goodVertices",
                "globalSuperTightHalo2016Filter",
                "HBHENoiseFilter",
                "HBHENoiseIsoFilter",
                "EcalDeadCellTriggerPrimitiveFilter",
                "BadPFMuonFilter",
            ],
            '2017': [
                "goodVertices",
                "globalSuperTightHalo2016Filter",
                "HBHENoiseFilter",
                "HBHENoiseIsoFilter",
                "EcalDeadCellTriggerPrimitiveFilter",
                "BadPFMuonFilter",
                "BadChargedCandidateFilter",
                "eeBadScFilter",
                "ecalBadCalibFilter",
            ],
            '2018': [
                "goodVertices",
                "globalSuperTightHalo2016Filter",
                "HBHENoiseFilter",
                "HBHENoiseIsoFilter",
                "EcalDeadCellTriggerPrimitiveFilter",
                "BadPFMuonFilter",
                "BadChargedCandidateFilter",
                "eeBadScFilter",
                "ecalBadCalibFilter",
            ],
        }[year]

        self._met_triggers = {
            '2016': [
                "PFMETNoMu110_PFMHTNoMu110_IDTight",
            ],
            '2017': [
                "PFMETNoMu120_PFMHTNoMu120_IDTight",
            ],
            '2018': [
                "PFMETNoMu120_PFMHTNoMu120_IDTight",
            ],
        }[year]

        # WPs for btagDeepFlavB (UL)
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation
        self._btagWPs = {
            '2016': {
                'medium': 0.6321,
            },
            '2017': {
                'medium': 0.4941,
            },
            '2018': {
                'medium': 0.4184,
            },
        }[year]

        jet_eta_bin = hist.Bin('jet_eta', r'Jet $\eta$', 20, -3., 3.)
        jet_msd_bin = hist.Bin('jet_msd', r'Jet $m_{sd}$ [GeV]', 8, 50, 210.)
        nn_hadhad_bin = hist.Bin('nn_hadhad',r'$NN_{\tau_{h}\tau_{h}}$', [0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.96,0.97,0.98,0.99,0.995,1.])
        nn_hadel_bin = hist.Bin('nn_hadel',r'$NN_{e\tau_{h}}$', [0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.96,0.97,0.98,0.99,0.995,0.996,0.997,0.998,0.999,1.])
        nn_hadmu_bin = hist.Bin('nn_hadmu',r'$NN_{\mu\tau_{h}}$', [0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.96,0.97,0.98,0.99,0.995,1.])
        nn_disc_bin = hist.Bin('nn_disc',r'$NN$', [0.,0.1,0.5,0.8,0.85,0.9,0.95,0.98,0.99,0.995,0.999,0.9995,0.9999,0.99999,0.999999,1.00001])
        massreg_bin = hist.Bin('massreg',r'$m_{NN}$', 21, 0., 210.)
        mt_lepmet_bin = hist.Bin('mt_lepmet', r'$m_{T}(\ell, MET)$', [0., 60., 500.])
        mt_jetmet_bin = hist.Bin('mt_jetmet', r'$m_{T}(j, MET)$', 10, 0., 1000.)
        oppbjet_pt_bin = hist.Bin('oppbjet_pt', r'Max opp. deepCSV-bjet $p_{T}$ [GeV]', 20, 0., 500)
        oppbtag_bin = hist.Bin('oppbtag', r'Max opp. deepCSV-b', 10, 0., 1)
        lep_pt_bin = hist.Bin('lep_pt', r'Lepton $p_{T}$ [GeV]', 8, 40, 200)
        lep_eta_bin = hist.Bin('lep_eta', r'Lepton $\eta$', 20, -3., 3.)
        jet_lsf3_bin = hist.Bin('lsf3', r'Jet LSF$_3$', 20, 0., 1.)
        lep_jet_dr_bin = hist.Bin('lep_jet_dr', r'$\Delta R(jet,lepton)$', 40, 0., 4.)
        lep_miso_bin = hist.Bin('miso', r'Lepton miniIso', [0.,0.025,0.05,0.075,0.1,0.2,0.4,0.6,0.8,1.])
        jet_jetlep_m_bin = hist.Bin('jetlep_m', r'Jet+lepton $m$ [GeV]', 20, 0, 600.)
        jet_jetmet_m_bin = hist.Bin('jetmet_m', r'Jet+MET $m$ [GeV]', 20, 0, 600.)
        jet_jetlepmet_m_bin = hist.Bin('jetlepmet_m', r'Jet+lepton+MET $m$ [GeV]', 20, 0, 600.)
        jetmet_dphi_bin = hist.Bin('jetmet_dphi', r'$\Delta\phi(jet,MET)$', 2, 0., 3.2)
        met_pt_bin = hist.Bin('met_pt', r'PuppiMET [GeV]', [20.,50.,100.,150.,200.,1000.])
        met_nopup_pt_bin = hist.Bin('met_nopup_pt', r'MET [GeV]', 10, 0, 800)
        h_pt_bin = hist.Bin('h_pt', r'h $p_{T}$ [GeV]', [250,280,300,350,400,500,600,1200])
        ntau_bin = hist.Bin('ntau',r'Number of taus',64,-0.5,63.5)
        antilep_bin = hist.Bin('antilep',r'Anti lepton veto',3,-1.5,1.5)
        genhtt_bin = hist.Bin('genhtt',r'hh,eh,mh,em,ee,mm (- for dr > 0.8)',4,-0.5,3.5)
        gentau1had_bin = hist.Bin('gentau1had',r'1pr,1pr+pi0,3pr',4,-0.5,3.5)
        gentau2had_bin = hist.Bin('gentau2had',r'1pr,1pr+pi0,3pr',4,-0.5,3.5)

        self._accumulator = processor.dict_accumulator({
            # dataset -> sumw
            'sumw': processor.defaultdict_accumulator(float),
            # dataset -> cut -> count
            'cutflow_hadhad': processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            'cutflow_hadhad_met': processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            'cutflow_hadhad_cr_b': processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            'cutflow_hadhad_cr_b_met': processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            'cutflow_hadhad_cr_b_mu': processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            'cutflow_hadhad_cr_b_mu_iso': processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            'cutflow_hadhad_cr_mu': processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            'cutflow_hadhad_cr_mu_iso': processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            'cutflow_hadel': processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            'cutflow_hadmu': processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            'cutflow_hadel_cr_b': processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            'cutflow_hadmu_cr_b': processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            'cutflow_hadel_cr_w': processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            'cutflow_hadmu_cr_w': processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            'cutflow_hadel_cr_qcd': processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            'cutflow_hadmu_cr_qcd': processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            'met_nn_kin': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('systematic', 'Systematic'),
                hist.Cat('region', 'Region'),
                met_pt_bin, massreg_bin, nn_disc_bin, jetmet_dphi_bin, h_pt_bin, antilep_bin,
            ),
        })


    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        dataset = events.metadata['dataset']
        isRealData = not hasattr(events, "genWeight")
        selection = PackedSelection(dtype="uint64")
        nevents = len(events)
        weights = Weights(nevents, storeIndividual=True)
        output = self.accumulator.identity()
        
        if not isRealData:
            output['sumw'][dataset] = ak.sum(events.genWeight)
            
        # trigger
        triggermasks = {}
        for channel in ["e","mu",'had']:
            good_trigs = 0
            #if isRealData:
            trigger = np.zeros(nevents, dtype='bool')
            for t in self._triggers[channel]:
                if t in events.HLT.fields:
                    trigger = trigger | events.HLT[t]
                    good_trigs += 1
            #selection.add('trigger'+channel, trigger)
            #del trigger
            if good_trigs < 1:
                raise ValueError("none of the following triggers found in dataset:", self._triggers[channel])
            #else:
                #selection.add('trigger'+channel, np.ones(nevents, dtype='bool'))
                #trigger = np.ones(nevents, dtype=np.bool)

            triggermasks[channel] = trigger
    
        met_trigger = np.zeros(nevents, dtype='bool')
        good_trigs = 0
        for t in self._met_triggers:
            if t in events.HLT.fields:
                met_trigger = met_trigger | events.HLT[t]
                good_trigs += 1
        if good_trigs < 1:
            raise ValueError("none of the following triggers found in dataset:", self._met_triggers)

        if isRealData:
            overlap_removal = isOverlap(events, dataset, self.triggers['e'] + self.triggers['mu'] + self.triggers['had'] + self.met_triggers, self._year)
        else:
            overlap_removal = np.ones(nevents, dtype=np.bool)

        met_filters = np.ones(nevents, dtype=np.bool)
        for t in self._metFilters:
            met_filters = met_filters & events.Flag[t]

        selection.add('met_filters', met_filters)

        selection.add('met_trigger',    met_trigger         & overlap_removal & met_filters)
        selection.add('hadel_trigger',  triggermasks['e']   & overlap_removal & met_filters)
        selection.add('hadmu_trigger',  triggermasks['mu']  & overlap_removal & met_filters)
        selection.add('hadhad_trigger', triggermasks['had'] & overlap_removal & met_filters)

        #jets
        if hasattr(events, 'FatJet'):
            fatjets = events.FatJet
        else:
            fatjets = events.CustomAK8Puppi

        fatjets['msdcorr'] = corrected_msoftdrop(fatjets)
        fatjets['rho'] = 2*np.log(fatjets.msdcorr/fatjets.pt)

        ak8jets = fatjets[
            (fatjets.pt > 200)
            & (np.abs(fatjets.eta) < 2.5)
        ]

        ak8jets_p4 = ak.zip(
            {
                'pt' : ak8jets.pt,
                'eta' : ak8jets.eta,
                'phi' : ak8jets.phi,
                'mass' : ak8jets.mass
            },
            behavior = vector.behavior,
            with_name='PtEtaPhiMLorentzVector'
        )

        met = events.PuppiMET 
        met_p4 = ak.zip(
            {
                'pt' : met.pt,
                'eta': 0,
                'phi': met.phi,
                'mass' : 0
            },
            behavior = vector.behavior,
            with_name='PtEtaPhiMLorentzVector'
        )

        ak8_met_pair = ak.cartesian( (ak8jets_p4, met_p4) )

        ak8s = ak8_met_pair[:,:,'0']
        mets = ak8_met_pair[:,:,'1']
        ak8_met_dphi = np.abs(ak8s.delta_phi(mets))

        best_ak8_idx = ak.argmin(ak8_met_dphi, axis=1, keepdims=True)
        #indexing ak8jets like ak8s works because there is only one met per event
        best_ak8 = ak.firsts(ak.to_regular(ak8jets[best_ak8_idx]))
        best_ak8_p4 = ak.firsts(ak8jets_p4[best_ak8_idx])
        best_ak8_met_dphi = ak.firsts(ak.to_regular(ak8_met_dphi[best_ak8_idx]))
        mt_jetmet = np.sqrt(2. * ak8s.pt * mets.pt * (1 - np.cos(ak8_met_dphi)))

        nn_disc_hadel = events.IN.hadel_v5p1
        nn_disc_hadmu = events.IN.hadmu_v5p1
        nn_disc_hadhad = events.IN.hadhad_v5p1_multi_Higgs

        massreg_hadel = events.MassReg.hadel_mass
        massreg_hadmu = events.MassReg.hadmu_mass
        massreg_hadhad = events.MassReg.hadhad_mass

        ptreg_hadel = events.MassReg.hadel_pt
        ptreg_hadmu = events.MassReg.hadmu_pt
        ptreg_hadhad = events.MassReg.hadhad_pt

        selection.add('ptreg_hadel', ptreg_hadel>280.0)
        selection.add('ptreg_hadmu', ptreg_hadmu>280.0)
        selection.add('ptreg_hadhad', ptreg_hadhad>280.0)

        selection.add('nn_disc_hadel', nn_disc_hadel>0.95)
        selection.add('nn_disc_hadmu', nn_disc_hadmu>0.95)
        selection.add('nn_disc_hadhad', nn_disc_hadhad>0.9999)
        
        selection.add('jetacceptance', (
            (best_ak8.pt > 200.)
            & (np.abs(best_ak8.eta) < 2.5)
            & (best_ak8.rho > -6.0)
            & (best_ak8.rho < -1.40)
        ))

        selection.add('jetacceptance450Inv', (
            (best_ak8.pt <= 450)
            & (np.abs(best_ak8.eta) < 2.5)
            & (best_ak8.rho > -6.0)
            & (best_ak8.rho < -1.40)
        ))

        selection.add('jetacceptance400', (
            (best_ak8.pt > 400)
            & (np.abs(best_ak8.eta) < 2.5)
            & (best_ak8.rho > -6.0)
            & (best_ak8.rho < -2.0)
        ))

        selection.add('jetacceptance450', (
            (best_ak8.pt > 450)
            & (np.abs(best_ak8.eta) < 2.5)
            & (best_ak8.rho > -6.0)
            & (best_ak8.rho < -2.1)
        ))

        selection.add('jet_msd', (best_ak8.msdcorr > 50.))
        selection.add("jetid", best_ak8.isTight)

        alljets = events.Jet[
            (events.Jet.pt>30.)
        ]

        hem_cleaning =(
            (
                (alljets.eta > -3.)
                & (alljets.eta < -1.3)
                & (alljets.phi > -1.57)
                & (alljets.phi < -0.87)
            )
            | (
                (met.phi > -1.62)
                & (met.pt < 470.)
                & (met.phi < -0.62)
            )
        )
        hem_cleaning = ak.all(hem_cleaning,1)

        selection.add('hem_cleaning', hem_cleaning)

        ak4jets = events.Jet[
            (events.Jet.pt > 30.)
            & np.abs((events.Jet.eta) < 2.5)
            & events.Jet.isTight
        ]

        #only first 5 jets for some reason
        ak4jets = ak4jets[:, :5]

        ak4jets_p4 = ak.zip(
            {
                'pt' : ak4jets.pt,
                'eta' : ak4jets.eta,
                'phi' : ak4jets.phi,
                'mass' : ak4jets.mass
            },
            behavior = vector.behavior,
            with_name = 'PtEtaPhiMLorentzVector'
        )
        
        ak4_ak8_pair = ak.cartesian( (ak4jets_p4, best_ak8_p4) )
        ak4s = ak4_ak8_pair[:,:,'0']
        ak8s = ak4_ak8_pair[:,:,'1']

        ak4_ak8_dr = ak4s.delta_r(ak8s)
        #again, indexing like this works because there is at most 1 best_ak8
        ak4_away = ak4jets[ak4_ak8_dr > 0.8]

        selection.add('antiak4btagMediumOppHem', 
                ak.max(ak4_away.btagDeepB, 1) <= self._btagWPs['medium'])
        selection.add('ak4btagMedium08', 
                ak.max(ak4_away.btagDeepB, 1) > self._btagWPs['medium'])

        ak4_met_pair = ak.cartesian( (ak4jets, met) )
        ak4s = ak4_met_pair[:,:,'0']
        mets = ak4_met_pair[:,:,'1']

        ak4_met_dphi = np.abs(ak4s.delta_phi(mets))
        minidx = ak.argmin(ak4_met_dphi, axis=1, keepdims=True)
        ak4_met_mindphi = ak4_met_dphi[minidx]
        selection.add('jetmet_dphi', ak.firsts(ak4_met_mindphi < np.pi/2))

        selection.add('met', met.pt > 20.)
        selection.add('methard', met.pt>150.)

        #muons
        goodmuon = (
            (events.Muon.pt > 25)
            & (np.abs(events.Muon.eta) < 2.4)
            & events.Muon.mediumId
        )
        ngoodmuons = ak.sum(goodmuon, 1)
        leadingmuon = ak.firsts(events.Muon[goodmuon])

        muons = (
            (events.Muon.pt > 10)
            & (np.abs(events.Muon.eta) < 2.4)
            & events.Muon.looseId
        )
        nmuons = ak.sum(muons ,1)

        #electrons
        goodelectrons = (
            (events.Electron.pt > 25)
            & (np.abs(events.Electron.eta) < 2.4)
            & (events.Electron.mvaFall17V2noIso_WP90)
        )
        ngoodelectrons = ak.sum(goodelectrons, 1)
        leadingelectron = ak.firsts(events.Electron[goodelectrons])

        electrons = (
            (events.Electron.pt > 10)
            & (np.abs(events.Electron.eta) < 2.4)
            & (events.Electron.cutBased >= events.Electron.LOOSE)
        )
        nelectrons = ak.sum(electrons, 1)

        #taus
        if self._year == '2018':
            tauAntiEleId = events.Tau.idAntiEle2018
        else:
            tauAntiEleId = events.Tau.idAntiEle

        loosetaus_el = (
            (events.Tau.pt > 20)
            & (np.abs(events.Tau.eta) < 2.3)
            & (tauAntiEleId >= 2)
        )
        goodtaus_el = (
            loosetaus_el
            & (tauAntiEleId >= 4)
        )
        goodtaus_mu = (
            (events.Tau.pt > 20)
            & (np.abs(events.Tau.eta) < 2.3)
            & (events.Tau.idAntiMu >= 1)
        )
        
        etaus_p4 = ak.zip(
            {
                'pt' : events.Tau[goodtaus_el].pt,
                'eta' : events.Tau[goodtaus_el].eta,
                'phi' : events.Tau[goodtaus_el].phi,
                'mass' : events.Tau[goodtaus_el].mass,
            },
            behavior = vector.behavior,
            with_name='PtEtaPhiMLorentzVector'
        )
        etausloose_p4 = ak.zip(
            {
                'pt' : events.Tau[loosetaus_el].pt,
                'eta' : events.Tau[loosetaus_el].eta,
                'phi' : events.Tau[loosetaus_el].phi,
                'mass' : events.Tau[loosetaus_el].mass,
            },
            behavior = vector.behavior,
            with_name='PtEtaPhiMLorentzVector'
        )
        mtaus_p4 = ak.zip(
            {
                'pt' : events.Tau[goodtaus_mu].pt,
                'eta' : events.Tau[goodtaus_mu].eta,
                'phi' : events.Tau[goodtaus_mu].phi,
                'mass' : events.Tau[goodtaus_mu].mass,
            },
            behavior = vector.behavior,
            with_name='PtEtaPhiMLorentzVector'
        )
        
        etau_ak8_pair = ak.cartesian((etaus_p4, best_ak8_p4))
        etau_ak8_dr = etau_ak8_pair[:,:,'0'].delta_r(etau_ak8_pair[:,:,'1'])

        etauloose_ak8_pair = ak.cartesian((etausloose_p4, best_ak8_p4))
        etauloose_ak8_dr = etauloose_ak8_pair[:,:,'0'].delta_r(etauloose_ak8_pair[:,:,'1'])

        mtau_ak8_pair = ak.cartesian((mtaus_p4, best_ak8_p4))
        mtau_ak8_dr = mtau_ak8_pair[:,:,'0'].delta_r(mtau_ak8_pair[:,:,'1'])

        selection.add('antiElId', ak.any(etau_ak8_dr < 0.8, -1))
        selection.add('antiMuId', ak.any(mtau_ak8_dr < 0.8, -1))

        one_el = (
            ~ak.any(electrons & ~goodelectrons, 1)
            & (nmuons == 0)
            & (ngoodmuons == 0)
            & (ngoodelectrons == 1)
        )

        one_mu = (
            ~ak.any(muons & ~goodmuon, 1)
            & (nelectrons ==0)
            & (ngoodelectrons==0)
            & (ngoodmuons == 1)
        )

        mu_p4 = ak.zip(
            {
                'pt' : leadingmuon.pt,
                'eta' : leadingmuon.eta,
                'phi' : leadingmuon.phi,
                'mass' : leadingmuon.mass
            },
            behavior = vector.behavior,
            with_name='PtEtaPhiMLorentzVector'
        )
        ele_p4 = ak.zip(
            {
                'pt' : leadingelectron.pt,
                'eta' : leadingelectron.eta,
                'phi' : leadingelectron.phi,
                'mass' : leadingelectron.mass
            },
            behavior = vector.behavior,
            with_name='PtEtaPhiMLorentzVector'
        )
        nolep_p4 = ak.zip(
            {
                'pt' : ak.zeros_like(leadingelectron.pt),
                'eta' : ak.zeros_like(leadingelectron.pt),
                'phi' : ak.zeros_like(leadingelectron.pt),
                'mass' : ak.zeros_like(leadingelectron.pt),
            },
            behavior = vector.behavior,
            with_name = 'PtEtaPhiMLorentzVector'
        )

        leadinglep = ak.where(one_el, ele_p4,
                ak.where(one_mu, mu_p4, nolep_p4))

        mu_miso = leadingmuon.miniPFRelIso_all
        el_miso = leadingelectron.miniPFRelIso_all
        leadinglepmiso = ak.where(one_el, el_miso,
                ak.where(one_mu, mu_miso, 0))

        mt_lepmet = np.sqrt(2.*leadinglep.pt*met_p4.pt*(1-np.cos(met_p4.delta_phi(leadinglep))))

        selection.add('mt_lepmet', (mt_lepmet < 60.))
        selection.add('mt_lepmetInv', (mt_lepmet>=60.))

        selection.add('noleptons', ( 
            (nmuons==0)
            & (nelectrons==0)
            & (ngoodmuons ==0)
            & (ngoodelectrons==0)
        ))

        selection.add('onemuon', one_mu)
        selection.add('oneelec', one_el)

        selection.add('muonkin', (
            (leadingmuon.pt > 30.)
            & (np.abs(leadingmuon.eta) < 2.4)
        ))
        selection.add('eleckin', (
            (leadingelectron.pt > 40.)
            & (np.abs(leadingelectron.eta) < 2.4)
        ))

        muon_ak8_dphi = np.abs(ak.fill_none(best_ak8_p4.delta_phi(mu_p4), 0))
        ele_ak8_dphi = np.abs(ak.fill_none(best_ak8_p4.delta_phi(ele_p4), 0)) 
        selection.add('muonDphiAK8', muon_ak8_dphi > 2*np.pi/3)
        selection.add('elecDphiAK8', ele_ak8_dphi > 2*np.pi/3)

        lep_ak8_dr = best_ak8_p4.delta_r(leadinglep)
        selection.add('lepDrAk8', lep_ak8_dr < 0.8)
        selection.add('lepDrAk8Inv', lep_ak8_dr >= 0.8)

        selection.add('muonIso', (
            ((leadingmuon.pt > 30)
            & (leadingmuon.pt < 55)
            & (leadingmuon.pfRelIso04_all < 0.25)
            ) 
            | ((leadingmuon.pt >= 55 )
            & (leadingmuon.miniPFRelIso_all < 0.1))
        ))

        selection.add('muonIsoInv', (
            ((leadingmuon.pt > 30)
            & (leadingmuon.pt < 55)
            & (leadingmuon.pfRelIso04_all >= 0.25)
            ) 
            | ((leadingmuon.pt >= 55) 
            & (leadingmuon.miniPFRelIso_all >= 0.1))
        ))

        selection.add('elecIso', (
            ((leadingelectron.pt > 30 )
            & (leadingelectron.pt < 120)
            & (leadingelectron.pfRelIso03_all < 0.1))
            | ((leadingelectron.pt >= 120)
            & (leadingelectron.miniPFRelIso_all < 0.1))
        ))

        selection.add('elecIsoInv', (
            ((leadingelectron.pt > 30 )
            & (leadingelectron.pt < 120)
            & (leadingelectron.pfRelIso03_all >= 0.1))
            | ((leadingelectron.pt >= 120)
            & (leadingelectron.miniPFRelIso_all >= 0.1))
        ))

        '''
        apply scale factors/weights/something
        I left this as untouched as I could without getting errors
        I make no claims that this does what it's suppossed to do
        This is because I don't really know what it's supposed to do
        See also corrections.py and common.py
        '''
        if isRealData:
            genflavor = ak.zeros_like(best_ak8.pt)
            w_hadhad = deepcopy(weights)
            w_hadhadmet = deepcopy(weights)
            w_hadel = deepcopy(weights)
            w_hadmu = deepcopy(weights)
            genHTauTauDecay = ak.zeros_like(best_ak8.pt)
            genHadTau1Decay = ak.zeros_like(best_ak8.pt)
            genHadTau2Decay = ak.zeros_like(best_ak8.pt)
            genHadTau2Decay = ak.zeros_like(best_ak8.pt)
            genTauTaudecay = ak.zeros_like(best_ak8.pt)
        else:
            weights.add('genweight', events.genWeight)
            add_pileup_weight(weights, events.Pileup.nPU, self._year)
            bosons = getBosons(events)
            genBosonPt = ak.fill_none(ak.pad_none(bosons.pt, 1, clip=True), 0)
            #I don't have an implementation of these
            #add_TopPtReweighting(weights, getParticles(events,6,6,['isLastCopy']).pt.pad(2, clip=True), self._year, dataset) #123 gives a weight of 1
            #add_VJets_NLOkFactor(weights, genBosonPt, self._year, dataset)
            genflavor = matchedBosonFlavor(best_ak8, bosons)
            genlepflavor = matchedBosonFlavorLep(best_ak8, bosons)
            genHTauTauDecay, genHadTau1Decay, genHadTau2Decay = getHTauTauDecayInfo(events,True)
            w_hadhad = deepcopy(weights)
            w_hadhadmet = deepcopy(weights)
            w_hadel = deepcopy(weights)
            w_hadmu = deepcopy(weights)
            #also need implementation here
            #add_LeptonSFs(w_hadel, leadinglep.pt, leadinglep.eta, self._year, "elec")
            #add_LeptonSFs(w_hadmu, leadinglep.pt, leadinglep.eta, self._year, "muon")
            #add_TriggerWeight(w_hadhad, best_ak8.msdcorr, best_ak8.pt, leadinglep.pt, self._year, "hadhad")

        regions = {
            'hadhad_signal': ['hadhad_trigger', 'noleptons', 'jetacceptance450', 'jet_msd', 'jetid', 'antiak4btagMediumOppHem', 'met'],# 'ptreg_hadhad'],
            'hadhad_signal_met': ['met_trigger', 'methard', 'noleptons', 'jetacceptance', 'jet_msd', 'jetid', 'antiak4btagMediumOppHem'],# 'ptreg_hadhad'],
            'hadhad_cr_b': ['hadhad_trigger', 'noleptons', 'jetacceptance450', 'jet_msd', 'jetid', 'ak4btagMedium08', 'met'],# 'ptreg_hadhad'],
            'hadhad_cr_b_met': ['met_trigger', 'methard', 'noleptons', 'jetacceptance', 'jet_msd', 'jetid', 'ak4btagMedium08'],# 'ptreg_hadhad'],
            'hadhad_cr_b_mu': ['hadmu_trigger', 'onemuon', 'muonkin', 'jetacceptance', 'jet_msd', 'jetid', 'ak4btagMedium08', 'met', 'lepDrAk8Inv', 'muonIsoInv'],# 'ptreg_hadhad'],#,'jetlsf'],
            'hadhad_cr_b_mu_iso': ['hadmu_trigger', 'onemuon', 'muonkin', 'jetacceptance', 'jet_msd', 'jetid', 'ak4btagMedium08', 'met', 'lepDrAk8Inv', 'muonIso'],# 'ptreg_hadhad'],#,'jetlsf'],
            'hadhad_cr_mu': ['hadmu_trigger', 'onemuon', 'muonkin', 'jetacceptance', 'jet_msd', 'jetid', 'antiak4btagMediumOppHem', 'met', 'lepDrAk8Inv', 'muonIsoInv'],# 'ptreg_hadhad'],#,'jetlsf'],
            'hadhad_cr_mu_iso': ['hadmu_trigger', 'onemuon', 'muonkin', 'jetacceptance', 'jet_msd', 'jetid', 'antiak4btagMediumOppHem', 'met', 'lepDrAk8Inv', 'muonIso'],# 'ptreg_hadhad'],#,'jetlsf'],
            'hadmu_signal': ['hadmu_trigger', 'onemuon', 'muonkin', 'jetacceptance', 'jet_msd', 'jetid', 'antiak4btagMediumOppHem', 'met', 'lepDrAk8', 'mt_lepmet', 'muonIso'],# 'ptreg_hadmu'],#, 'antiMuId', 'jetlsf'],
            'hadel_signal': ['hadel_trigger', 'oneelec', 'eleckin', 'jetacceptance', 'jet_msd', 'jetid', 'antiak4btagMediumOppHem', 'met', 'lepDrAk8', 'mt_lepmet', 'elecIso'],# 'ptreg_hadel'],#, 'antiElId', 'jetlsf'],
            'hadmu_cr_qcd': ['hadmu_trigger', 'onemuon', 'muonkin', 'jetacceptance', 'jet_msd', 'jetid', 'lepDrAk8', 'muonIsoInv'],# 'ptreg_hadmu'],#, 'antiMuId','jetlsf'],
            'hadel_cr_qcd': ['hadel_trigger', 'oneelec', 'eleckin', 'jetacceptance', 'jet_msd', 'jetid', 'lepDrAk8', 'elecIsoInv'],# 'ptreg_hadel'],#, 'antiElId','jetlsf'],
            'hadmu_cr_b': ['hadmu_trigger', 'onemuon', 'muonkin', 'jetacceptance', 'jet_msd', 'jetid', 'ak4btagMedium08', 'met', 'lepDrAk8', 'muonIso'],# 'ptreg_hadmu'],#, 'antiMuId','jetlsf'],
            'hadel_cr_b': ['hadel_trigger', 'oneelec', 'eleckin', 'jetacceptance', 'jet_msd', 'jetid', 'ak4btagMedium08', 'met', 'lepDrAk8', 'elecIso'],# 'ptreg_hadel'],#, 'antiElId','jetlsf'],
            'hadmu_cr_w': ['hadmu_trigger', 'onemuon', 'muonkin', 'jetacceptance', 'jet_msd', 'jetid', 'antiak4btagMediumOppHem', 'met', 'lepDrAk8', 'mt_lepmetInv', 'muonIso'],# 'ptreg_hadmu'],#, 'antiMuId','jetlsf'],
            'hadel_cr_w': ['hadel_trigger', 'oneelec', 'eleckin', 'jetacceptance', 'jet_msd', 'jetid', 'antiak4btagMediumOppHem', 'met', 'lepDrAk8', 'mt_lepmetInv', 'elecIso'],# 'ptreg_hadel'],#, 'antiElId','jetlsf'],
            #'noselection': [],
        }
        if (self._year == '2018'):
            for r in regions:
                regions[r].append('hem_cleaning')

        w_dict = {
            'hadhad_signal': w_hadhad,
            'hadhad_signal_met': w_hadhadmet,
            'hadhad_cr_b': w_hadhad,
            'hadhad_cr_b_met': w_hadhadmet,
            'hadhad_cr_b_mu': w_hadmu,
            'hadhad_cr_b_mu_iso': w_hadmu,
            'hadhad_cr_mu': w_hadmu,
            'hadhad_cr_mu_iso': w_hadmu,
            'hadmu_signal': w_hadmu,
            'hadel_signal': w_hadel,
            'hadmu_cr_qcd': w_hadmu,
            'hadel_cr_qcd': w_hadel,
            'hadmu_cr_b': w_hadmu,
            'hadel_cr_b': w_hadel,
            'hadmu_cr_w': w_hadmu,
            'hadel_cr_w': w_hadel,
        }

        allcuts_hadel = set()
        allcuts_hadmu = set()
        allcuts_hadel_cr_b = set()
        allcuts_hadmu_cr_b = set()
        allcuts_hadel_cr_w = set()
        allcuts_hadmu_cr_w = set()
        allcuts_hadel_cr_qcd = set()
        allcuts_hadmu_cr_qcd = set()
        allcuts_hadhad = set()
        allcuts_hadhad_met = set()
        allcuts_hadhad_cr_b = set()
        allcuts_hadhad_cr_b_met = set()
        allcuts_hadhad_cr_b_mu = set()
        allcuts_hadhad_cr_b_mu_iso = set()
        allcuts_hadhad_cr_mu = set()
        allcuts_hadhad_cr_mu_iso = set()
        output['cutflow_hadel'][dataset]['none'] += float(w_dict['hadel_signal'].weight().sum())
        output['cutflow_hadmu'][dataset]['none'] += float(w_dict['hadmu_signal'].weight().sum())
        output['cutflow_hadel_cr_b'][dataset]['none'] += float(w_dict['hadel_cr_b'].weight().sum())
        output['cutflow_hadmu_cr_b'][dataset]['none'] += float(w_dict['hadmu_cr_b'].weight().sum())
        output['cutflow_hadel_cr_w'][dataset]['none'] += float(w_dict['hadel_cr_w'].weight().sum())
        output['cutflow_hadmu_cr_w'][dataset]['none'] += float(w_dict['hadmu_cr_w'].weight().sum())
        output['cutflow_hadel_cr_qcd'][dataset]['none'] += float(w_dict['hadel_cr_qcd'].weight().sum())
        output['cutflow_hadmu_cr_qcd'][dataset]['none'] += float(w_dict['hadmu_cr_qcd'].weight().sum())
        output['cutflow_hadhad'][dataset]['none'] += float(w_dict['hadhad_signal'].weight().sum())
        output['cutflow_hadhad_met'][dataset]['none'] += float(w_dict['hadhad_signal_met'].weight().sum())
        output['cutflow_hadhad_cr_b'][dataset]['none'] += float(w_dict['hadhad_cr_b'].weight().sum())
        output['cutflow_hadhad_cr_b_met'][dataset]['none'] += float(w_dict['hadhad_cr_b_met'].weight().sum())
        output['cutflow_hadhad_cr_b_mu'][dataset]['none'] += float(w_dict['hadhad_cr_b_mu'].weight().sum())
        output['cutflow_hadhad_cr_b_mu_iso'][dataset]['none'] += float(w_dict['hadhad_cr_b_mu_iso'].weight().sum())
        output['cutflow_hadhad_cr_mu'][dataset]['none'] += float(w_dict['hadhad_cr_mu'].weight().sum())
        output['cutflow_hadhad_cr_mu_iso'][dataset]['none'] += float(w_dict['hadhad_cr_mu_iso'].weight().sum())

        for cut in regions['hadel_signal']+['ptreg_hadel','nn_disc_hadel']:
            allcuts_hadel.add(cut)
            output['cutflow_hadel'][dataset][cut] += float(w_dict['hadel_signal'].weight()[selection.all(*allcuts_hadel)].sum())
        for cut in regions['hadmu_signal']+['ptreg_hadmu','nn_disc_hadmu']:
            allcuts_hadmu.add(cut)
            output['cutflow_hadmu'][dataset][cut] += float(w_dict['hadmu_signal'].weight()[selection.all(*allcuts_hadmu)].sum())
        for cut in regions['hadel_cr_b']+['ptreg_hadel','nn_disc_hadel']:
            allcuts_hadel_cr_b.add(cut)
            output['cutflow_hadel_cr_b'][dataset][cut] += float(w_dict['hadel_cr_b'].weight()[selection.all(*allcuts_hadel_cr_b)].sum())
        for cut in regions['hadmu_cr_b']+['ptreg_hadmu','nn_disc_hadmu']:
            allcuts_hadmu_cr_b.add(cut)
            output['cutflow_hadmu_cr_b'][dataset][cut] += float(w_dict['hadmu_cr_b'].weight()[selection.all(*allcuts_hadmu_cr_b)].sum())
        for cut in regions['hadel_cr_w']+['ptreg_hadel','nn_disc_hadel']:
            allcuts_hadel_cr_w.add(cut)
            output['cutflow_hadel_cr_w'][dataset][cut] += float(w_dict['hadel_cr_w'].weight()[selection.all(*allcuts_hadel_cr_w)].sum())
        for cut in regions['hadmu_cr_w']+['ptreg_hadmu','nn_disc_hadmu']:
            allcuts_hadmu_cr_w.add(cut)
            output['cutflow_hadmu_cr_w'][dataset][cut] += float(w_dict['hadmu_cr_w'].weight()[selection.all(*allcuts_hadmu_cr_w)].sum())
        for cut in regions['hadel_cr_qcd']+['ptreg_hadel','nn_disc_hadel']:
            allcuts_hadel_cr_qcd.add(cut)
            output['cutflow_hadel_cr_qcd'][dataset][cut] += float(w_dict['hadel_cr_qcd'].weight()[selection.all(*allcuts_hadel_cr_qcd)].sum())
        for cut in regions['hadmu_cr_qcd']+['ptreg_hadmu','nn_disc_hadmu']:
            allcuts_hadmu_cr_qcd.add(cut)
            output['cutflow_hadmu_cr_qcd'][dataset][cut] += float(w_dict['hadmu_cr_qcd'].weight()[selection.all(*allcuts_hadmu_cr_qcd)].sum())
        for cut in regions['hadhad_signal']+['ptreg_hadhad','nn_disc_hadhad']:
            allcuts_hadhad.add(cut)
            output['cutflow_hadhad'][dataset][cut] += float(w_dict['hadhad_signal'].weight()[selection.all(*allcuts_hadhad)].sum())
        for cut in regions['hadhad_signal_met']+['ptreg_hadhad','nn_disc_hadhad']:
            allcuts_hadhad_met.add(cut)
            output['cutflow_hadhad_met'][dataset][cut] += float(w_dict['hadhad_signal_met'].weight()[selection.all(*allcuts_hadhad_met)].sum())
        for cut in regions['hadhad_cr_b']+['ptreg_hadhad','nn_disc_hadhad']:
            allcuts_hadhad_cr_b.add(cut)
            output['cutflow_hadhad_cr_b'][dataset][cut] += float(w_dict['hadhad_cr_b'].weight()[selection.all(*allcuts_hadhad_cr_b)].sum())
        for cut in regions['hadhad_cr_b_met']+['ptreg_hadhad','nn_disc_hadhad']:
            allcuts_hadhad_cr_b_met.add(cut)
            output['cutflow_hadhad_cr_b_met'][dataset][cut] += float(w_dict['hadhad_cr_b_met'].weight()[selection.all(*allcuts_hadhad_cr_b_met)].sum())
        for cut in regions['hadhad_cr_b_mu']+['ptreg_hadhad','nn_disc_hadhad']:
            allcuts_hadhad_cr_b_mu.add(cut)
            output['cutflow_hadhad_cr_b_mu'][dataset][cut] += float(w_dict['hadhad_cr_b_mu'].weight()[selection.all(*allcuts_hadhad_cr_b_mu)].sum())
        for cut in regions['hadhad_cr_b_mu_iso']+['ptreg_hadhad','nn_disc_hadhad']:
            allcuts_hadhad_cr_b_mu_iso.add(cut)
            output['cutflow_hadhad_cr_b_mu_iso'][dataset][cut] += float(w_dict['hadhad_cr_b_mu_iso'].weight()[selection.all(*allcuts_hadhad_cr_b_mu_iso)].sum())
        for cut in regions['hadhad_cr_mu']+['ptreg_hadhad','nn_disc_hadhad']:
            allcuts_hadhad_cr_mu.add(cut)
            output['cutflow_hadhad_cr_mu'][dataset][cut] += float(w_dict['hadhad_cr_mu'].weight()[selection.all(*allcuts_hadhad_cr_mu)].sum())
        for cut in regions['hadhad_cr_mu_iso']+['ptreg_hadhad','nn_disc_hadhad']:
            allcuts_hadhad_cr_mu_iso.add(cut)
            output['cutflow_hadhad_cr_mu_iso'][dataset][cut] += float(w_dict['hadhad_cr_mu_iso'].weight()[selection.all(*allcuts_hadhad_cr_mu_iso)].sum()) 

        systematics = [
            None,
            #'TopPtReweightUp',
        ]
        if isRealData:
            systematics = [None]

        def fill(region, systematic, wmod=None, realData=False):
            selections = regions[region]
            cut = selection.all(*selections)
            sname = 'nominal' if systematic is None else systematic
            if wmod is None:
                if not realData:
                    weight = w_dict[region].weight(modifier=systematic)
                else:
                    weight = w_dict[region].weight(modifier=None)
            else:
                weight = w_dict[region].weight() * wmod

            #weight = zero if cut fails
            #seems reasonable enough, right?
            weight = ak.where(cut, weight, 0)

            '''
            I /think/ this is equivalent to what was here before...
            '''
            if 'hadhad' in region:
                nn_disc = nn_disc_hadhad
                massreg = massreg_hadhad
                ptreg = ptreg_hadhad
                antilep = ak.any(mtau_ak8_dr < 0.8, -1) & ak.any(etauloose_ak8_dr < 0.8, -1)
            if 'hadel' in region:
                nn_disc = nn_disc_hadel
                massreg = massreg_hadel
                ptreg = ptreg_hadel
                antilep = ak.any(etau_ak8_dr < 0.8, -1)
                antilep = ak.any(etauloose_ak8_dr < 0.8, -1) + antilep - 1
            if 'hadmu' in region:
                nn_disc = nn_disc_hadmu
                massreg = massreg_hadmu
                ptreg = ptreg_hadmu
                antilep = ak.any(mtau_ak8_dr < 0.8, -1)

            bmaxind = ak.argmax(ak4_away.btagDeepB, -1)

            output['met_nn_kin'].fill(
                dataset=dataset,
                region=region,
                systematic=sname,
                met_pt=normalize(met_p4.pt),
                massreg=normalize(massreg),
                nn_disc=normalize(nn_disc),
                jetmet_dphi=normalize(best_ak8_met_dphi),
                h_pt=normalize(ptreg),
                antilep=antilep,
                weight=weight,
            )


        for region in regions:
            for systematic in systematics:
                fill(region, systematic, realData=isRealData)

        return output

    def postprocess(self, accumulator):
        return accumulator
