from functools import partial
import numpy as np
from coffea import processor, hist
from uproot_methods import TLorentzVectorArray
import awkward
from copy import deepcopy
from .common import (
    getBosons,
    getParticles,
    matchedBosonFlavor,
    matchedBosonFlavorLep,
    getHTauTauDecayInfo,
    isOverlap,
)
from .corrections import (
    corrected_msoftdrop,
    n2ddt_shift,
    add_pileup_weight,
    add_VJets_NLOkFactor,
    add_TopPtReweighting,
    add_jetTriggerWeight,
    add_TriggerWeight,
    add_LeptonSFs,
)
#from .btag import BTagEfficiency, BTagCorrector

# for old pancakes
from coffea.nanoaod.methods import collection_methods, FatJet
collection_methods['CustomAK8Puppi'] = FatJet
collection_methods['CustomAK8PuppiSubjet'] = FatJet
FatJet.subjetmap['CustomAK8Puppi'] = 'CustomAK8PuppiSubjet'

class HtautauProcessor_Z(processor.ProcessorABC):
    def __init__(self, year='2017'):
        self._year = year

        #self._btagSF = BTagCorrector(year, 'medium')
        self._btagWPs = {
            'medium': {
                '2016': 0.6321,
                '2017': 0.4941,
                '2018': 0.4184,
            },
        }

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
                "ecalBadCalibFilterV2",
            ],
        }

        self._hadel_triggers = {
            '2016': [
                #'Ele35_WPTight_Gsf',
'Ele50_CaloIdVT_GsfTrkIdT_PFJet165','Ele115_CaloIdVT_GsfTrkIdT',
#"Ele15_IsoVVVL_PFHT450_PFMET50",
"Ele15_IsoVVVL_PFHT600",
                'PFHT800',
                'PFHT900',
                'AK8PFJet360_TrimMass30',
                'AK8PFHT700_TrimR0p1PT0p03Mass50',
                'PFHT650_WideJetMJJ950DEtaJJ1p5',
                'PFHT650_WideJetMJJ900DEtaJJ1p5',
                #'AK8DiPFJet280_200_TrimMass30_BTagCSV_p20',
                'PFJet450',
            ],
            '2017': [
                'Ele35_WPTight_Gsf',
#'Ele50_CaloIdVT_GsfTrkIdT_PFJet165',
'Ele115_CaloIdVT_GsfTrkIdT',
#"Ele15_IsoVVVL_PFHT450_PFMET50",
#"Ele15_IsoVVVL_PFHT600",
                #'AK8PFJet330_PFAK8BTagCSV_p17',
                #'PFHT1050',
                #'AK8PFJet400_TrimMass30',
                #'AK8PFJet420_TrimMass30',
                #'AK8PFHT800_TrimMass50',
                #'PFJet500',
                #'AK8PFJet500',
            ],
            '2018': [
                #'Ele35_WPTight_Gsf',
'Ele50_CaloIdVT_GsfTrkIdT_PFJet165','Ele115_CaloIdVT_GsfTrkIdT',
#"Ele15_IsoVVVL_PFHT450_PFMET50",
"Ele15_IsoVVVL_PFHT600",
                'AK8PFJet400_TrimMass30',
                'AK8PFJet420_TrimMass30',
                'AK8PFHT800_TrimMass50',
                'PFHT1050',
                'PFJet500',
                'AK8PFJet500',
                # 'AK8PFJet330_PFAK8BTagCSV_p17', not present in 2018D?
                #'AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4',
                #'AK4PFJet30',
            ],
        }

        self._hadmu_triggers = {
            '2016': [
                'Mu50','Mu55',
#"Mu15_IsoVVVL_PFHT450_PFMET50",
"Mu15_IsoVVVL_PFHT600",
                'PFHT800',
                'PFHT900',
                'AK8PFJet360_TrimMass30',
                'AK8PFHT700_TrimR0p1PT0p03Mass50',
                'PFHT650_WideJetMJJ950DEtaJJ1p5',
                'PFHT650_WideJetMJJ900DEtaJJ1p5',
                #'AK8DiPFJet280_200_TrimMass30_BTagCSV_p20',
                'PFJet450',
            ],
            '2017': [
                'Mu50',#'Mu55',
                #'IsoMu27',
#"Mu15_IsoVVVL_PFHT450_PFMET50",
#"Mu15_IsoVVVL_PFHT600",
                #'AK8PFJet330_PFAK8BTagCSV_p17',
                #'PFHT1050',
                #'AK8PFJet400_TrimMass30',
                #'AK8PFJet420_TrimMass30',
                #'AK8PFHT800_TrimMass50',
                #'PFJet500',
                #'AK8PFJet500',
            ],
            '2018': [
                'Mu50',#'Mu55',
#"Mu15_IsoVVVL_PFHT450_PFMET50",
"Mu15_IsoVVVL_PFHT600",
                'AK8PFJet400_TrimMass30',
                'AK8PFJet420_TrimMass30',
                'AK8PFHT800_TrimMass50',
                'PFHT1050',
                'PFJet500',
                'AK8PFJet500',
                # 'AK8PFJet330_PFAK8BTagCSV_p17', not present in 2018D?
                #'AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4',
                #'AK4PFJet30',
            ],
        }

        self._hadhad_triggers = {
            '2016': [
                'PFHT800',
                'PFHT900',
                'AK8PFJet360_TrimMass30',
                'AK8PFHT700_TrimR0p1PT0p03Mass50',
                'PFHT650_WideJetMJJ950DEtaJJ1p5',
                'PFHT650_WideJetMJJ900DEtaJJ1p5',
                #'AK8DiPFJet280_200_TrimMass30_BTagCSV_p20',
                'PFJet450',
                'DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg',
                'DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg',
                'DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg',
                'DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg',
            ],
            '2017': [
                #'AK8PFJet330_PFAK8BTagCSV_p17',
                'PFHT1050',
                'AK8PFJet400_TrimMass30',
                'AK8PFJet420_TrimMass30',
                'AK8PFHT800_TrimMass50',
                'PFJet500',
                'AK8PFJet500',
                #'DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg',
                #'DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg',
                #'DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg',
                #'DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg',
                #'MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1',
                #'MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr',
            ],
            '2018': [
                'AK8PFJet400_TrimMass30',
                'AK8PFJet420_TrimMass30',
                'AK8PFHT800_TrimMass50',
                'PFHT1050',
                'PFJet500',
                'AK8PFJet500',
                # 'AK8PFJet330_PFAK8BTagCSV_p17', not present in 2018D?
                #'AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4',
                #'AK4PFJet30',
                'DoubleMediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg',
                'DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg',
                'DoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg',
                'DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1_Reg',
                'MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1',
                'MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr',
                'MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1',
                'MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1',
            ],
        }

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
        }

        #jet_pt_bin = hist.Bin('jet_pt', r'Jet $p_{T}$ [GeV]', [200.,250.,300.,350.,400.,450.,500.,550.,600.,1200.])
        jet_pt_bin = hist.Bin('jet_pt', r'Jet $p_{T}$ [GeV]', [200.,450.,1200.])
        jet_eta_bin = hist.Bin('jet_eta', r'Jet $\eta$', 20, -3., 3.)
        jet_msd_bin = hist.Bin('jet_msd', r'Jet $m_{sd}$ [GeV]', 8, 50, 210.)
        #nn_hadhad_bin = hist.Bin('nn_hadhad',r'$NN_{\tau_{h}\tau_{h}}$',20,0.,1.)
        #nn_hadel_bin = hist.Bin('nn_hadel',r'$NN_{e\tau_{h}}$',20,0.,1.)
        #nn_hadmu_bin = hist.Bin('nn_hadmu',r'$NN_{\mu\tau_{h}}$',20,0.,1.)
        nn_hadhad_bin = hist.Bin('nn_hadhad',r'$NN_{\tau_{h}\tau_{h}}$', [0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.96,0.97,0.98,0.99,0.995,1.])
        nn_hadel_bin = hist.Bin('nn_hadel',r'$NN_{e\tau_{h}}$', [0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.96,0.97,0.98,0.99,0.995,0.996,0.997,0.998,0.999,1.])
        nn_hadmu_bin = hist.Bin('nn_hadmu',r'$NN_{\mu\tau_{h}}$', [0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.96,0.97,0.98,0.99,0.995,1.])
        nn_disc_bin = hist.Bin('nn_disc',r'$NN$', [0.,0.1,0.5,0.8,0.85,0.9,0.95,0.96,0.97,0.98,0.99,0.995,0.999,0.9995,0.9999,1.])
        massreg_bin = hist.Bin('massreg',r'$m_{NN}$', 21, 0., 210.)
        mt_lepmet_bin = hist.Bin('mt_lepmet', r'$m_{T}(\ell, MET)$', [0., 80., 500.])
        mt_jetmet_bin = hist.Bin('mt_jetmet', r'$m_{T}(j, MET)$', 10, 0., 1000.)
        oppbjet_pt_bin = hist.Bin('oppbjet_pt', r'Max opp. deepCSV-bjet $p_{T}$ [GeV]', 20, 0., 500)
        oppbtag_bin = hist.Bin('oppbtag', r'Max opp. deepCSV-b', 10, 0., 1)
        lep_pt_bin = hist.Bin('lep_pt', r'Lepton $p_{T}$ [GeV]', 8, 40, 200)
        lep_eta_bin = hist.Bin('lep_eta', r'Lepton $\eta$', 20, -3., 3.)
        jet_lsf3_bin = hist.Bin('lsf3', r'Jet LSF$_3$', 20, 0., 1.)
        lep_jet_dr_bin = hist.Bin('lep_jet_dr', r'$\Delta R(jet,lepton)$', 40, 0., 4.)
        #lep_miso_bin = hist.Bin('miso', r'Lepton miniIso', 20, 0., 0.1)
        lep_miso_bin = hist.Bin('miso', r'Lepton miniIso', [0.,0.025,0.05,0.075,0.1,0.2,0.4,0.6,0.8,1.])
        jet_jetlep_m_bin = hist.Bin('jetlep_m', r'Jet+lepton $m$ [GeV]', 20, 0, 600.)
        jet_jetmet_m_bin = hist.Bin('jetmet_m', r'Jet+MET $m$ [GeV]', 20, 0, 600.)
        jet_jetlepmet_m_bin = hist.Bin('jetlepmet_m', r'Jet+lepton+MET $m$ [GeV]', 20, 0, 600.)
        jetmet_dphi_bin = hist.Bin('jetmet_dphi', r'$\Delta\phi(jet,MET)$', 2, 0., 3.2)
        met_pt_bin = hist.Bin('met_pt', r'PuppiMET [GeV]', [20.,50.,200.,1000.])
        met_nopup_pt_bin = hist.Bin('met_nopup_pt', r'MET [GeV]', 10, 0, 800)
        n2ddt_bin = hist.Bin('n2ddt', r'N_{2}^{DDT}', 2, -1.,1.)
        h_pt_bin = hist.Bin('h_pt', r'h $p_{T}$ [GeV]', [250,300,400,1200])
        ntau_bin = hist.Bin('ntau',r'Number of taus',64,-0.5,63.5)
        antilep_bin = hist.Bin('antilep',r'Anti lepton veto',2,-0.5,1.5)
        genhtt_bin = hist.Bin('genhtt',r'hh,eh,mh,em,ee,mm (- for dr > 0.8)',4,-0.5,3.5)
        gentau1had_bin = hist.Bin('gentau1had',r'1pr,1pr+pi0,3pr',4,-0.5,3.5)
        gentau2had_bin = hist.Bin('gentau2had',r'1pr,1pr+pi0,3pr',4,-0.5,3.5)
        bos_pt_bin = hist.Bin('bos_pt', r'Boson $p_{T}$ [GeV]', [0,200,250,300,350,400,450,500,1500])
        bos_m_bin = hist.Bin('bos_m', r'Boson mass [GeV]', 12, 60., 120.)

        self._accumulator = processor.dict_accumulator({
            # dataset -> sumw
            'sumw': processor.defaultdict_accumulator(float),
            # dataset -> cut -> count
            'cutflow_mumu': processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            'cutflow_mumu_nodr': processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            #'btagWeight': hist.Hist('Events', hist.Cat('dataset', 'Dataset'), hist.Bin('val', 'BTag correction', 50, 0, 2)), #FIXME
            #'jet_nn_kin': hist.Hist(
            #    'Events',
            #    hist.Cat('dataset', 'Dataset'),
            #    hist.Cat('region', 'Region'),
            #    #jet_pt_bin, jet_msd_bin, nn_disc_bin,
            #    jet_pt_bin, massreg_bin, nn_disc_bin, jet_msd_bin, antilep_bin,
            #),
            #'lep_nn_kin': hist.Hist(
            #    'Events',
            #    hist.Cat('dataset', 'Dataset'),
            #    hist.Cat('region', 'Region'),
            #    #lep_miso_bin, jet_msd_bin, nn_disc_bin,
            #    #lep_pt_bin, jet_msd_bin, nn_disc_bin,
            #    lep_pt_bin, massreg_bin, nn_disc_bin, mt_lepmet_bin, antilep_bin,
            #),
            'gen_nn_kin': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                #lep_miso_bin, jet_msd_bin, nn_disc_bin,
                #lep_pt_bin, jet_msd_bin, nn_disc_bin,
                bos_pt_bin, bos_m_bin,
            ),
            
        })

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        dataset = events.metadata['dataset']
        isRealData = 'genWeight' not in events.columns
        selection = processor.PackedSelection()
        weights = processor.Weights(len(events))
        output = self.accumulator.identity()
        if not isRealData:
            output['sumw'][dataset] += events.genWeight.sum()

        trigger_hadmu = np.zeros(events.size, dtype='bool')
        for t in self._hadmu_triggers[self._year]:
            trigger_hadmu = trigger_hadmu | events.HLT[t]

        if (isRealData): overlap_removal = isOverlap(events,dataset,self._hadmu_triggers[self._year])
        else: overlap_removal = np.ones(events.size, dtype='bool')

        met_filters = np.ones(events.size, dtype='bool')
        for t in self._metFilters[self._year]:
            met_filters = met_filters & events.Flag[t]

        selection.add('hadmu_trigger',  trigger_hadmu  & overlap_removal & met_filters)

        try:
            fatjets = events.FatJet
        except AttributeError:
            # early pancakes
            fatjets = events.CustomAK8Puppi
        fatjets['msdcorr'] = corrected_msoftdrop(fatjets)
        fatjets['rho'] = 2 * np.log(fatjets.msdcorr / fatjets.pt)
        fatjets['n2ddt'] = fatjets.n2b1 - n2ddt_shift(fatjets, year=self._year)

        candidatejets = fatjets[
            # https://github.com/DAZSLE/BaconAnalyzer/blob/master/Analyzer/src/VJetLoader.cc#L269
            (fatjets.pt > 200)
            #& (abs(fatjets.eta) < 2.5)
            #& (fatjets.isTight)
        ]#[:, :2]
        met_p4 = TLorentzVectorArray.from_ptetaphim(awkward.JaggedArray.fromiter([[v] for v in events.PuppiMET.pt]), awkward.JaggedArray.fromiter([[v] for v in np.zeros(events.size)]), awkward.JaggedArray.fromiter([[v] for v in events.PuppiMET.phi]), awkward.JaggedArray.fromiter([[v] for v in np.zeros(events.size)]))
        met_nopup_p4 = TLorentzVectorArray.from_ptetaphim(awkward.JaggedArray.fromiter([[v] for v in events.MET.pt]), awkward.JaggedArray.fromiter([[v] for v in np.zeros(events.size)]), awkward.JaggedArray.fromiter([[v] for v in events.MET.phi]), awkward.JaggedArray.fromiter([[v] for v in np.zeros(events.size)]))
        ak8_met_pair = candidatejets.cross(met_p4)
        ak8_met_dphi = abs(ak8_met_pair.i0.delta_phi(ak8_met_pair.i1))
        #aligned_jet = ak8_met_dphi == ak8_met_dphi.min()
        #best_jet_idx = (ak8_met_pair.i0 + aligned_jet * ak8_met_pair.i1).pt.argmax()
        best_jet_idx = ak8_met_dphi.argmin()
        #best_jet_idx = candidatejets.pt.argmax()
        candidatejet = candidatejets[best_jet_idx]
        jetmet_dphi_ak8 = ak8_met_dphi[best_jet_idx]

        mt_jetmet = np.sqrt(2.*ak8_met_pair.i0.pt.fillna(0)*ak8_met_pair.i1.pt*(ak8_met_pair.i1.pt.ones_like()-np.cos(ak8_met_dphi)))

        nn_disc_hadhad = awkward.JaggedArray.fromiter([[v] for v in events.IN.hadhad_v4p1])[candidatejet.pt.pad(1, clip=True).fillna(0.)>200.]
        nn_disc_hadel  = awkward.JaggedArray.fromiter([[v] for v in events.GRU.hadel_v6p1])[candidatejet.pt.pad(1, clip=True).fillna(0.)>200.]
        nn_disc_hadmu  = awkward.JaggedArray.fromiter([[v] for v in events.GRU.hadmu_v6p1])[candidatejet.pt.pad(1, clip=True).fillna(0.)>200.]

        massreg_hadhad = awkward.JaggedArray.fromiter([[v] for v in events.MassReg.hadhad_mass])[candidatejet.pt.pad(1, clip=True).fillna(0.)>200.]
        massreg_hadel  = awkward.JaggedArray.fromiter([[v] for v in events.MassReg.hadel_mass])[candidatejet.pt.pad(1, clip=True).fillna(0.)>200.]
        massreg_hadmu  = awkward.JaggedArray.fromiter([[v] for v in events.MassReg.hadmu_mass])[candidatejet.pt.pad(1, clip=True).fillna(0.)>200.]

        ptreg_hadhad = awkward.JaggedArray.fromiter([[v] for v in events.MassReg.hadhad_pt])[candidatejet.pt.pad(1, clip=True).fillna(0.)>200.]
        ptreg_hadel  = awkward.JaggedArray.fromiter([[v] for v in events.MassReg.hadel_pt])[candidatejet.pt.pad(1, clip=True).fillna(0.)>200.]
        ptreg_hadmu  = awkward.JaggedArray.fromiter([[v] for v in events.MassReg.hadmu_pt])[candidatejet.pt.pad(1, clip=True).fillna(0.)>200.]

        candidatejet_rho = 2 * np.log(candidatejet.msdcorr / candidatejet.pt)
        selection.add('jetacceptance', (
            (candidatejet.pt > 200)
            #& (candidatejet.msdcorr > 40.)
            & (abs(candidatejet.eta) < 2.4)
            & (candidatejet_rho > -6.)
            & (candidatejet_rho < -1.40)
        ).any())
        selection.add('jetacceptance450Inv', (
            (candidatejet.pt <= 450)
            #& (candidatejet.msdcorr > 40.)
            & (abs(candidatejet.eta) < 2.4)
            & (candidatejet_rho > -6.)
            & (candidatejet_rho < -1.40)
        ).any())
        selection.add('jetacceptance400', (
            (candidatejet.pt > 400)
            #& (candidatejet.msdcorr > 40.)
            & (abs(candidatejet.eta) < 2.4)
            & (candidatejet_rho > -6.)
            & (candidatejet_rho < -2.)
        ).any())
        selection.add('jetacceptance450', (
            (candidatejet.pt > 450)
            #& (candidatejet.msdcorr > 40.)
            & (abs(candidatejet.eta) < 2.4)
            & (candidatejet_rho > -6.)
            & (candidatejet_rho < -2.1)
        ).any())
        selection.add('jet_msd', (candidatejet.msdcorr > 40.).any())
        selection.add('jetid', (candidatejet.isTight).any())
        selection.add('n2ddt', (candidatejet.n2ddt < 0.).any())
        #print(np.histogram(candidatejet.pt.fillna(0).flatten()))

        jets = events.Jet[
            (events.Jet.pt > 30.)
            & (abs(events.Jet.eta) < 2.5)
            & (events.Jet.isTight)
        ]
        # only consider first 4 jets to be consistent with old framework
        jets = jets[:, :5]
        ak4_ak8_pair = jets.cross(candidatejet, nested=True)
        ak4_ak8_dphi = abs(ak4_ak8_pair.i0.delta_phi(ak4_ak8_pair.i1))
        ak4_ak8_dr = ak4_ak8_pair.i0.delta_r(ak4_ak8_pair.i1)
        #ak4_opposite = jets[(ak4_ak8_dphi > np.pi / 2).all()]
        ak4_away = jets[(ak4_ak8_dr > 0.8).all()]
        #selection.add('antiak4btagMediumOppHem', ak4_away.btagDeepB.max() < BTagEfficiency.btagWPs[self._year]['medium'])
        #selection.add('ak4btagMedium08', ak4_away.btagDeepB.max() > BTagEfficiency.btagWPs[self._year]['medium'])
        selection.add('antiak4btagMediumOppHem', ak4_away.btagDeepB.max() < self._btagWPs['medium'][self._year])
        selection.add('ak4btagMedium08', ak4_away.btagDeepB.max() > self._btagWPs['medium'][self._year])

        selection.add('met', events.PuppiMET.pt > 20.)
        selection.add('methard', events.PuppiMET.pt > 200.)

        muons = (
            (events.Muon.pt > 10)
            & (abs(events.Muon.eta) < 2.4)
            #& (events.Muon.pfRelIso04_all < 0.25)
            #& (np.abs(events.Muon.dz) < 0.1)
            #& (np.abs(events.Muon.dxy) < 0.05)
            & (events.Muon.looseId).astype(bool)
            #& (events.Muon.highPtId).astype(bool)
        )
        nmuons = muons.sum()
        themuons = events.Muon[muons].pad(2, clip=True)

        #selection.add('twomuons', (nmuons == 2) & (themuons.pt.pad(1, clip=True) > 55).any())
        selection.add('twomuons', (nmuons == 2) & (themuons[:,0].pt > 55).any() & (themuons[:,0].charge!=themuons[:,1].charge).any())

        #selection.add('jetlsf', (
        #    (candidatejet.lsf3 > 0.7).any()
        #))

        mu0_p4 = TLorentzVectorArray.from_ptetaphim(awkward.JaggedArray.fromiter([[v] for v in themuons[:,0].pt]), awkward.JaggedArray.fromiter([[v] for v in themuons[:,0].eta]), awkward.JaggedArray.fromiter([[v] for v in themuons[:,0].phi]), awkward.JaggedArray.fromiter([[v] for v in themuons[:,0].mass]))
        mu1_p4 = TLorentzVectorArray.from_ptetaphim(awkward.JaggedArray.fromiter([[v] for v in themuons[:,1].pt]), awkward.JaggedArray.fromiter([[v] for v in themuons[:,1].eta]), awkward.JaggedArray.fromiter([[v] for v in themuons[:,1].phi]), awkward.JaggedArray.fromiter([[v] for v in themuons[:,1].mass]))
        Zcand = mu0_p4 + mu1_p4

        mu_pair = mu0_p4.cross(mu1_p4)
        mu_dr = mu_pair.i0.delta_r(mu_pair.i1)
        selection.add('mu_dr', ((mu_dr < 0.8) & (mu_dr > 0.1)).any())

        if not isRealData:
            weights.add('genweight', events.genWeight)
            add_pileup_weight(weights, events.Pileup.nPU, self._year, dataset)
            bosons = getBosons(events)
            genBosonPt = bosons.pt.pad(1, clip=True).fillna(0)
            add_TopPtReweighting(weights, getParticles(events,6,6,['isLastCopy']).pt.pad(2, clip=True), self._year, dataset) #123 gives a weight of 1
            add_VJets_NLOkFactor(weights, genBosonPt, self._year, dataset)
            add_LeptonSFs(weights, mu0_p4.pt, mu0_p4.eta, self._year, "muon_TRIG")

        regions = {
            'mumu_signal': ['hadmu_trigger', 'twomuons', 'mu_dr'],
            'mumu_nodr': ['hadmu_trigger', 'twomuons'],
            #'noselection': [],
        }
        w_dict = {
            'mumu_signal': weights,
            'mumu_nodr': weights,
        }

        allcuts_mumu = set()
        allcuts_mumu_nodr = set()
        output['cutflow_mumu'][dataset]['none'] += float(w_dict['mumu_signal'].weight().sum())
        output['cutflow_mumu_nodr'][dataset]['none'] += float(w_dict['mumu_nodr'].weight().sum())
        for cut in regions['mumu_signal']:
            allcuts_mumu.add(cut)
            output['cutflow_mumu'][dataset][cut] += float(w_dict['mumu_signal'].weight()[selection.all(*allcuts_mumu)].sum())
        for cut in regions['mumu_nodr']:
            allcuts_mumu_nodr.add(cut)
            output['cutflow_mumu_nodr'][dataset][cut] += float(w_dict['mumu_nodr'].weight()[selection.all(*allcuts_mumu_nodr)].sum())

        systematics = [
            None,
            #'jet_triggerUp',
            #'jet_triggerDown',
            #'btagWeightUp',
            #'btagWeightDown',
            #'btagEffStatUp',
            #'btagEffStatDown',
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
                    weight = w_dict[region].weight(modifier=systematic)[cut]
                else:
                    weight = w_dict[region].weight(modifier=None)[cut]
            else:
                weight = w_dict[region].weight()[cut] * wmod[cut]

            def normalize(val):
                return val[cut].pad(1, clip=True).fillna(0).flatten()

            output['gen_nn_kin'].fill(
                dataset=dataset,
                region=region,
                bos_pt=normalize(Zcand.pt),
                bos_m=normalize(Zcand.mass),
                weight=weight,
            )


        for region in regions:
            for systematic in systematics:
                fill(region, systematic, realData=isRealData)
        #    if 'GluGluHToTauTau' in dataset:
        #        for i in range(9):
        #            fill(region, 'LHEScale_%d' % i, events.LHEScaleWeight[:, i])
        #        for c in events.LHEWeight.columns[1:]:
        #            fill(region, 'LHEWeight_%s' % c, events.LHEWeight[c])

        return output

    def postprocess(self, accumulator):
        return accumulator
