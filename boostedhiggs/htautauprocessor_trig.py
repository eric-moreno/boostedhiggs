from functools import partial
import numpy as np
from coffea import processor, hist
from uproot_methods import TLorentzVectorArray
import awkward
from .common import (
    getBosons,
    matchedBosonFlavor,
)
from .corrections import (
    corrected_msoftdrop,
    n2ddt_shift,
    add_pileup_weight,
    add_VJets_NLOkFactor,
    add_jetTriggerWeight,
    add_TriggerWeight,
)
#from .btag import BTagEfficiency, BTagCorrector

# for old pancakes
from coffea.nanoaod.methods import collection_methods, FatJet
collection_methods['CustomAK8Puppi'] = FatJet
collection_methods['CustomAK8PuppiSubjet'] = FatJet
FatJet.subjetmap['CustomAK8Puppi'] = 'CustomAK8PuppiSubjet'


class HtautauProcessor_Trigger(processor.ProcessorABC):
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
                #'Ele35_WPTight_Gsf',
                'Ele50_CaloIdVT_GsfTrkIdT_PFJet165','Ele115_CaloIdVT_GsfTrkIdT',
#"Ele15_IsoVVVL_PFHT450_PFMET50",
"Ele15_IsoVVVL_PFHT600",
                #'AK8PFJet330_PFAK8BTagCSV_p17',
                'PFHT1050',
                'AK8PFJet400_TrimMass30',
                'AK8PFJet420_TrimMass30',
                'AK8PFHT800_TrimMass50',
                'PFJet500',
                'AK8PFJet500',
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
                'Mu50',#'Mu55',
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
#"Mu15_IsoVVVL_PFHT450_PFMET50",
"Mu15_IsoVVVL_PFHT600",
                #'AK8PFJet330_PFAK8BTagCSV_p17',
                'PFHT1050',
                'AK8PFJet400_TrimMass30',
                'AK8PFJet420_TrimMass30',
                'AK8PFHT800_TrimMass50',
                'PFJet500',
                'AK8PFJet500',
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
                'DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg',
                'DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg',
                'DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg',
                'DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg',
                'MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1',
                'MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr',
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

        self._ref_triggers = {
            '2016': {
                "hadhad":[
                    #"PFMETNoMu110_PFMHTNoMu110_IDTight",
                    'Mu50',
                ],
                "hadel":[
                    #"PFMETNoMu110_PFMHTNoMu110_IDTight",
                    'Mu50',
                ],
                "hadmu":[
                    #"PFMETNoMu110_PFMHTNoMu110_IDTight",
                    'Ele35_WPTight_Gsf',
                    'Ele50_CaloIdVT_GsfTrkIdT_PFJet165',
                    'Ele115_CaloIdVT_GsfTrkIdT',
                ],
            },
            '2017': {
                "hadhad":[
                    #"PFMETNoMu110_PFMHTNoMu110_IDTight",
                    'Mu50',
                ],
                "hadel":[
                    #"PFMETNoMu110_PFMHTNoMu110_IDTight",
                    'Mu50',
                ],
                "hadmu":[
                    #"PFMETNoMu110_PFMHTNoMu110_IDTight",
                    'Ele35_WPTight_Gsf',
                    'Ele50_CaloIdVT_GsfTrkIdT_PFJet165',
                    'Ele115_CaloIdVT_GsfTrkIdT',
                ],
            },
            '2018': {
                "hadhad":[
                    #"PFMETNoMu110_PFMHTNoMu110_IDTight",
                    'Mu50',
                ],
                "hadel":[
                    #"PFMETNoMu110_PFMHTNoMu110_IDTight",
                    'Mu50',
                ],
                "hadmu":[
                    #"PFMETNoMu110_PFMHTNoMu110_IDTight",
                    #'Ele35_WPTight_Gsf',
                    'Ele50_CaloIdVT_GsfTrkIdT_PFJet165',
                    'Ele115_CaloIdVT_GsfTrkIdT',
                ],
            },
        }

        h_pt_bin = hist.Bin('h_pt', r'Higgs $p_{T}$ [GeV]', 10, 200, 700)
        jet_pt_bin = hist.Bin('jet_pt', r'Jet $p_{T}$ [GeV]', 10, 300, 800)
        lep_pt_bin = hist.Bin('lep_pt', r'Lepton $p_{T}$ [GeV]', 10, 20., 140)
        jet_msd_bin = hist.Bin('jet_msd', r'Jet $m_{SD}$ [GeV]', 10, 0., 100.)
        jet_lep_dr_bin = hist.Bin('jet_lep_dr', r'Jet-lepton $\Delta R$', 10, 0., 1.)
        lep_miso_bin = hist.Bin('lep_miso', r'Lepton miniIso', 10, 0., 1.)
        n2ddt_bin = hist.Bin('n2ddt',r'$N_{2}^{DDT}$',10,-0.25,0.25)

        self._accumulator = processor.dict_accumulator({
            # dataset -> sumw
            'sumw': processor.defaultdict_accumulator(float),
            # dataset -> cut -> count
            'cutflow_hadhad': processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            'cutflow_hadel': processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            'cutflow_hadmu': processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            #'btagWeight': hist.Hist('Events', hist.Cat('dataset', 'Dataset'), hist.Bin('val', 'BTag correction', 50, 0, 2)), #FIXME
            'trigeff_h': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('trig_pass_hadhad', r'Trigger Pass Bit (HadHad)', 2, -0.5, 1.5),
                hist.Bin('trig_pass_hadel', r'Trigger Pass Bit (HadEl)', 2, -0.5, 1.5),
                hist.Bin('trig_pass_hadmu', r'Trigger Pass Bit (HadMu)', 2, -0.5, 1.5),
                hist.Bin('trig_pass_ref', r'Trigger Pass Bit (Reference)', 2, -0.5, 1.5),
                h_pt_bin, jet_pt_bin, lep_pt_bin,
            ),
            'trigeff_m': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('trig_pass_hadhad', r'Trigger Pass Bit (HadHad)', 2, -0.5, 1.5),
                hist.Bin('trig_pass_hadel', r'Trigger Pass Bit (HadEl)', 2, -0.5, 1.5),
                hist.Bin('trig_pass_hadmu', r'Trigger Pass Bit (HadMu)', 2, -0.5, 1.5),
                hist.Bin('trig_pass_ref', r'Trigger Pass Bit (Reference)', 2, -0.5, 1.5),
                jet_msd_bin, jet_pt_bin, lep_pt_bin,
            ),
            'trigeff_dr': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('trig_pass_hadhad', r'Trigger Pass Bit (HadHad)', 2, -0.5, 1.5),
                hist.Bin('trig_pass_hadel', r'Trigger Pass Bit (HadEl)', 2, -0.5, 1.5),
                hist.Bin('trig_pass_hadmu', r'Trigger Pass Bit (HadMu)', 2, -0.5, 1.5),
                hist.Bin('trig_pass_ref', r'Trigger Pass Bit (Reference)', 2, -0.5, 1.5),
                jet_lep_dr_bin, jet_pt_bin, lep_pt_bin,
            ),
            'trigeff_miso': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                hist.Bin('trig_pass_hadhad', r'Trigger Pass Bit (HadHad)', 2, -0.5, 1.5),
                hist.Bin('trig_pass_hadel', r'Trigger Pass Bit (HadEl)', 2, -0.5, 1.5),
                hist.Bin('trig_pass_hadmu', r'Trigger Pass Bit (HadMu)', 2, -0.5, 1.5),
                hist.Bin('trig_pass_ref', r'Trigger Pass Bit (Reference)', 2, -0.5, 1.5),
                lep_miso_bin, jet_pt_bin, lep_pt_bin,
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

        met_filters = np.ones(events.size, dtype='bool')
        for t in self._metFilters[self._year]:
            met_filters = met_filters & events.Flag[t]

        selection.add('met_filters', met_filters)

        trigger_ref = {}
        trigger_ref["hadhad"] = np.zeros(events.size, dtype='bool')
        trigger_ref["hadel"] = np.zeros(events.size, dtype='bool')
        trigger_ref["hadmu"] = np.zeros(events.size, dtype='bool')
        for t in self._ref_triggers[self._year]["hadhad"]:
            trigger_ref["hadhad"] = trigger_ref["hadhad"] | events.HLT[t]
        for t in self._ref_triggers[self._year]["hadel"]:
            trigger_ref["hadel"] = trigger_ref["hadel"] | events.HLT[t]
        for t in self._ref_triggers[self._year]["hadmu"]:
            trigger_ref["hadmu"] = trigger_ref["hadmu"] | events.HLT[t]

        trigger_hadhad = np.zeros(events.size, dtype='bool')
        for t in self._hadhad_triggers[self._year]:
            trigger_hadhad = trigger_hadhad | events.HLT[t]
        if isRealData:
            selection.add('hadhad_trigger', trigger_ref["hadhad"])
        else:
            selection.add('hadhad_trigger', np.ones(events.size, dtype='bool'))

        trigger_hadmu = np.zeros(events.size, dtype='bool')
        for t in self._hadmu_triggers[self._year]:
            trigger_hadmu = trigger_hadmu | events.HLT[t]
        if isRealData:
            selection.add('hadmu_trigger', trigger_ref["hadmu"])
        else:
            selection.add('hadmu_trigger', np.ones(events.size, dtype='bool'))

        trigger_hadel = np.zeros(events.size, dtype='bool')
        for t in self._hadel_triggers[self._year]:
            trigger_hadel = trigger_hadel | events.HLT[t]
        if isRealData:
            selection.add('hadel_trigger', trigger_ref["hadel"])
        else:
            selection.add('hadel_trigger', np.ones(events.size, dtype='bool'))
        #print(np.histogram(trigger))

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
            & (abs(fatjets.eta) < 2.5)
            & (fatjets.isTight)
        ]#[:, :2]
        met_p4 = TLorentzVectorArray.from_ptetaphim(awkward.JaggedArray.fromiter([[v] for v in events.MET.pt]), awkward.JaggedArray.fromiter([[v] for v in np.zeros(events.size)]), awkward.JaggedArray.fromiter([[v] for v in events.MET.phi]), awkward.JaggedArray.fromiter([[v] for v in np.zeros(events.size)]))
        ak8_met_pair = candidatejets.cross(met_p4)
        ak8_met_dphi = abs(ak8_met_pair.i0.delta_phi(ak8_met_pair.i1))
        #aligned_jet = ak8_met_dphi == ak8_met_dphi.min()
        #best_jet_idx = (ak8_met_pair.i0 + aligned_jet * ak8_met_pair.i1).pt.argmax()
        best_jet_idx = ak8_met_dphi.argmin()
        #best_jet_idx = candidatejets.pt.argmax()
        candidatejet = candidatejets[best_jet_idx]
        candidatejet_rho = 2 * np.log(candidatejet.msdcorr / candidatejet.pt)
        selection.add('jetacceptance', (
            (candidatejet.pt > 300)
            #& (candidatejet.msdcorr > 40.)
            & (abs(candidatejet.eta) < 2.4)
            & (candidatejet_rho > -6.)
            & (candidatejet_rho < -2.1)
        ).any())
        selection.add('jetid', (candidatejet.isTight).any())
        #selection.add('n2ddt', (candidatejet.n2ddt < 0.).any())
        #print(np.histogram(candidatejet.pt.fillna(0).flatten()))

        jets = events.Jet[
            (events.Jet.pt > 30.)
            & (abs(events.Jet.eta) < 2.5)
            & (events.Jet.isTight)
        ]
        # only consider first 4 jets to be consistent with old framework
        jets = jets[:, :4]
        ak4_ak8_pair = jets.cross(candidatejet, nested=True)
        dphi = abs(ak4_ak8_pair.i0.delta_phi(ak4_ak8_pair.i1))
        ak4_opposite = jets[(dphi > np.pi / 2).all()]
        #selection.add('antiak4btagMediumOppHem', ak4_opposite.btagDeepB.max() < BTagEfficiency.btagWPs[self._year]['medium'])
        selection.add('antiak4btagMediumOppHem', ak4_opposite.btagDeepB.max() < self._btagWPs['medium'][self._year])
        ak4_away = jets[(dphi > 0.8).all()]
        #selection.add('ak4btagMedium08', ak4_away.btagDeepB.max() > BTagEfficiency.btagWPs[self._year]['medium'])
        selection.add('ak4btagMedium08', ak4_away.btagDeepB.max() > self._btagWPs['medium'][self._year])

        selection.add('met0', events.MET.pt > 0.)
        selection.add('met50', events.MET.pt > 50.)
        selection.add('met100', events.MET.pt > 100.)
        selection.add('met150', events.MET.pt > 150.)
        selection.add('met200', events.MET.pt > 200.)

        el_loose_cuts = [(np.bitwise_and(np.right_shift(events.Electron.vidNestedWPBitmap,events.Electron.vidNestedWPBitmap.ones_like()*(3*k)),events.Electron.vidNestedWPBitmap.ones_like()*7) >= events.Electron.LOOSE) for k in range(10) if k != 7]
        el_tight_cuts = [(np.bitwise_and(np.right_shift(events.Electron.vidNestedWPBitmap,events.Electron.vidNestedWPBitmap.ones_like()*(3*k)),events.Electron.vidNestedWPBitmap.ones_like()*7) >= events.Electron.TIGHT) for k in range(10) if k != 7]
        #                  (MinPtCut,GsfEleSCEtaMultiRangeCut,GsfEleDEtaInSeedCut,GsfEleDPhiInCut,GsfEleFull5x5SigmaIEtaIEtaCut,GsfEleHadronicOverEMEnergyScaledCut,GsfEleEInverseMinusPInverseCut,GsfEleRelPFIsoScaledCut,GsfEleConversionVetoCut,GsfEleMissingHitsCut)
        #                   0       ,1                       ,2                  ,3              ,4                            ,5                                  ,6                             ,7                      ,8                      ,9

        elmask_loose = el_loose_cuts[0].ones_like().astype(bool)
        for m in el_loose_cuts: elmask_loose = elmask_loose & m

        elmask_tight = el_tight_cuts[0].ones_like().astype(bool)
        for m in el_tight_cuts: elmask_tight = elmask_tight & m

        goodmuon = (
            (events.Muon.pt > 20)
            & (np.abs(events.Muon.eta) < 2.4)
            #& (events.Muon.sip3d < 4)
            #& (np.abs(events.Muon.dz) < 0.1)
            #& (np.abs(events.Muon.dxy) < 0.05)
            & (events.Muon.mediumId).astype(bool)
            #& (events.Muon.highPtId).astype(bool)
            #& (events.Muon.miniPFRelIso_all < 0.1)
        )
        ngoodmuons = goodmuon.sum()
        leadingmuon = events.Muon[goodmuon].pad(1, clip=True)

        goodelec = (
            (events.Electron.pt > 20)
            & (abs(events.Electron.eta) < 2.5)
            #& (events.Electron.cutBased >= events.Electron.TIGHT)
            #& (events.Electron.cutBased_HEEP).astype(bool)
            & elmask_tight
            #& (events.Electron.miniPFRelIso_all < 0.1)
        )
        ngoodelecs = goodelec.sum()
        leadingelec = events.Electron[goodelec].pad(1, clip=True)

        nmuons = (
            (events.Muon.pt > 15)
            & (abs(events.Muon.eta) < 2.4)
            #& (events.Muon.pfRelIso04_all < 0.25)
            #& (np.abs(events.Muon.dz) < 0.1)
            #& (np.abs(events.Muon.dxy) < 0.05)
            & (events.Muon.looseId).astype(bool)
            #& (events.Muon.highPtId).astype(bool)
        ).sum()

        nelectrons = (
            (events.Electron.pt > 15)
            & (abs(events.Electron.eta) < 2.5)
            & (events.Electron.cutBased >= events.Electron.VETO)
            #& (events.Electron.cutBased_HEEP).astype(bool)
            #& elmask_loose
        ).sum()

        trigmuon = (
            (events.Muon.pt > 55)
            & (np.abs(events.Muon.eta) < 2.4)
            #& (events.Muon.sip3d < 4)
            #& (np.abs(events.Muon.dz) < 0.1)
            #& (np.abs(events.Muon.dxy) < 0.05)
            & (events.Muon.tightId).astype(bool)
            #& (events.Muon.highPtId).astype(bool)
            #& (events.Muon.miniPFRelIso_all < 0.1)
        )
        ntrigmuons = trigmuon.sum()

        trigelec = (
            (events.Electron.pt > 55)
            & (abs(events.Electron.eta) < 2.5)
            #& (events.Electron.cutBased >= events.Electron.TIGHT)
            #& (events.Electron.cutBased_HEEP).astype(bool)
            & elmask_tight
            #& (events.Electron.miniPFRelIso_all < 0.1)
        )
        ntrigelecs = trigelec.sum()

        ntaus = np.zeros(events.size, dtype='bool')

        #lepsel = ((nmuons <= 1) & (nelectrons == 0) & (ntaus == 0) & (ngoodmuons == 1) & (ngoodelecs == 0)) | ((nmuons == 0) & (nelectrons <= 1) & (ntaus == 0) & (ngoodmuons == 0) & (ngoodelecs == 1))
        mu_p4 = TLorentzVectorArray.from_ptetaphim(leadingmuon.pt.fillna(0),leadingmuon.eta.fillna(0),leadingmuon.phi.fillna(0),leadingmuon.mass.fillna(0))
#[(goodmuon & ((nmuons == 1) & (nelectrons == 0) & (ntaus == 0) & (ngoodmuons == 1)))]
        muon_ak8_pair = mu_p4.cross(candidatejet, nested=True)
        el_p4 = TLorentzVectorArray.from_ptetaphim(leadingelec.pt.fillna(0),leadingelec.eta.fillna(0),leadingelec.phi.fillna(0),leadingelec.mass.fillna(0))
#[(goodelec & ((nmuons == 0) & (nelectrons == 1) & (ntaus == 0) & (ngoodelecs == 1)))]
        elec_ak8_pair = el_p4.cross(candidatejet, nested=True)
        #leadinglep = awkward.concatenate([mu_p4, el_p4], axis=1).pad(1, clip=True)
        leadinglep = {}
        leadinglep["hadhad"] = mu_p4
        leadinglep["hadel"] = el_p4
        leadinglep["hadmu"] = mu_p4

        mu_miso = leadingmuon.miniPFRelIso_all.fillna(0)
        el_miso = leadingelec.miniPFRelIso_all.fillna(0)
        leadinglep_miso = {}
        leadinglep_miso["hadhad"] = mu_miso.pad(1, clip=True)
        leadinglep_miso["hadel"] = el_miso.pad(1, clip=True)
        leadinglep_miso["hadmu"] = mu_miso.pad(1, clip=True)

        selection.add('noleptons', (nmuons == 0) & (nelectrons == 0) & (ntaus == 0) & (ngoodmuons == 0) & (ngoodelecs == 0))
        selection.add('trigmuons', (ntrigmuons >= 1) & (ntaus == 0))
        selection.add('trigelecs', (ntrigelecs >= 1) & (ntaus == 0))
        selection.add('onemuon', (ngoodmuons == 1) & (ntaus == 0))
        selection.add('oneelec', (ngoodelecs == 1) & (ntaus == 0))
        selection.add('muonkin', (
            (leadingmuon.pt > 20.)
            & (abs(leadingmuon.eta) < 2.1)
        ).all())
        selection.add('muonDphiAK8', (
            abs(muon_ak8_pair.i0.delta_phi(muon_ak8_pair.i1)) > 2*np.pi/3
        ).all().all())
        selection.add('eleckin', (
            (leadingelec.pt > 20.)
            & (abs(leadingelec.eta) < 2.4)
        ).all())
        selection.add('elecDphiAK8', (
            abs(elec_ak8_pair.i0.delta_phi(elec_ak8_pair.i1)) > 2*np.pi/3
        ).all().all())

        lep_ak8_pair = {c:leadinglep[c].cross(candidatejet) for c in ["hadhad","hadel","hadmu"]}
        #selection.add('lepDrAK8', (
        #    (lep_ak8_pair.i0.delta_r(lep_ak8_pair.i1) < 0.8).all()
        #))

        #jet_lep_p4 = lep_ak8_pair.i0 + lep_ak8_pair.i1
        #met_jl_pair = met_p4.cross(jet_lep_p4)#, nested=True)
        #jet_lep_met_p4 = met_jl_pair.i0 + met_jl_pair.i1
        #jet_met_p4 = ak8_met_pair.i0[best_jet_idx] + ak8_met_pair.i1[best_jet_idx]

        if isRealData:
            genflavor = candidatejet.pt.zeros_like()
            genBosonPt = candidatejet.pt.zeros_like()
        else:
            weights.add('genweight', events.genWeight)
            add_pileup_weight(weights, events.Pileup.nPU, self._year, dataset)
            bosons = getBosons(events)
            genBosonPt = bosons.pt.pad(1, clip=True).fillna(0)
            add_VJets_NLOkFactor(weights, genBosonPt, self._year, dataset)
            genflavor = matchedBosonFlavor(candidatejet, bosons)
            #add_TriggerWeight(weights, candidatejet.msdcorr, candidatejet.pt, self._year)
            #output['btagWeight'].fill(dataset=dataset, val=self._btagSF.addBtagWeight(weights, ak4_away)) #FIXME

        regions = {
            #'hadhad_signal': ['met_filters','jetacceptance', 'hadhad_trigger', 'jetid', 'n2ddt', 'met', 'noleptons'],#, 'antiak4btagMediumOppHem'
            **{'hadhad_signal_%s'%metcut: ['met_filters','jetacceptance', 'hadhad_trigger', 'jetid', 'met%s'%metcut, 'trigmuons'] for metcut in ['0','50','100']},#,'150','200']},#, 'antiak4btagMediumOppHem'
            **{'hadhad_signal_withb_%s'%metcut: ['met_filters','jetacceptance', 'hadhad_trigger', 'jetid', 'met%s'%metcut, 'ak4btagMedium08', 'trigmuons'] for metcut in ['0','50','100']},#,'150','200']},#, 'antiak4btagMediumOppHem'
            **{'hadmu_signal_%s'%metcut: ['met_filters','jetacceptance', 'hadmu_trigger', 'jetid', 'met%s'%metcut, 'onemuon', 'trigelecs'] for metcut in ['0','50','100']},#,'150','200']},#, 'antiak4btagMediumOppHem'
            **{'hadel_signal_%s'%metcut: ['met_filters','jetacceptance', 'hadel_trigger', 'jetid', 'met%s'%metcut, 'oneelec', 'trigmuons'] for metcut in ['0','50','100']},#,'150','200']},#, 'antiak4btagMediumOppHem'
            #'hadmu_control': ['met_filters','jetacceptance', 'hadmu_trigger', 'jetid', 'ak4btagMedium08', 'met', 'onemuon', 'muonkin', 'muonDphiAK8'],
            #'hadel_control': ['met_filters','jetacceptance', 'hadel_trigger', 'jetid', 'ak4btagMedium08', 'met', 'oneelec', 'eleckin', 'elecDphiAK8'],
            #'noselection': [],
        }

        allcuts_hadel = set()
        allcuts_hadmu = set()
        allcuts_hadhad = set()
        output['cutflow_hadel'][dataset]['none'] += float(weights.weight().sum())
        output['cutflow_hadmu'][dataset]['none'] += float(weights.weight().sum())
        output['cutflow_hadhad'][dataset]['none'] += float(weights.weight().sum())
        for cut in regions['hadel_signal_50']:
            allcuts_hadel.add(cut)
            output['cutflow_hadel'][dataset][cut] += float(weights.weight()[selection.all(*allcuts_hadel)].sum())
        for cut in regions['hadmu_signal_50']:
            allcuts_hadmu.add(cut)
            output['cutflow_hadmu'][dataset][cut] += float(weights.weight()[selection.all(*allcuts_hadmu)].sum())
        for cut in regions['hadhad_signal_withb_50']:
            allcuts_hadhad.add(cut)
            output['cutflow_hadhad'][dataset][cut] += float(weights.weight()[selection.all(*allcuts_hadhad)].sum())

        systematics = [
            None,
            #'jet_triggerUp',
            #'jet_triggerDown',
            #'btagWeightUp',
            #'btagWeightDown',
            #'btagEffStatUp',
            #'btagEffStatDown',
        ]

        def fill(region, systematic, wmod=None):
            if (region.startswith("hadhad")): chan = "hadhad"
            elif (region.startswith("hadel")): chan = "hadel"
            elif (region.startswith("hadmu")): chan = "hadmu"
            else: chan = ""
            selections = regions[region]
            cut = selection.all(*selections)
            sname = 'nominal' if systematic is None else systematic
            if wmod is None:
                weight = weights.weight(modifier=systematic)[cut]
            else:
                weight = weights.weight()[cut] * wmod[cut]

            def normalize(val):
                return val[cut].pad(1, clip=True).fillna(0).flatten()

            output['trigeff_h'].fill(
                dataset=dataset,
                region=region,
                trig_pass_hadhad=trigger_hadhad[cut],
                trig_pass_hadel=trigger_hadel[cut],
                trig_pass_hadmu=trigger_hadmu[cut],
                trig_pass_ref=trigger_ref[chan][cut],
                h_pt=normalize(genBosonPt),
                jet_pt=normalize(candidatejet.pt),
                lep_pt=normalize(leadinglep[chan].pt),
                #weight=weight,
            )

            output['trigeff_m'].fill(
                dataset=dataset,
                region=region,
                trig_pass_hadhad=trigger_hadhad[cut],
                trig_pass_hadel=trigger_hadel[cut],
                trig_pass_hadmu=trigger_hadmu[cut],
                trig_pass_ref=trigger_ref[chan][cut],
                jet_msd=normalize(candidatejet.msdcorr),
                jet_pt=normalize(candidatejet.pt),
                lep_pt=normalize(leadinglep[chan].pt),
                #n2ddt=normalize(candidatejet.n2ddt),
                #weight=weight,
            )

            output['trigeff_dr'].fill(
                dataset=dataset,
                region=region,
                trig_pass_hadhad=trigger_hadhad[cut],
                trig_pass_hadel=trigger_hadel[cut],
                trig_pass_hadmu=trigger_hadmu[cut],
                trig_pass_ref=trigger_ref[chan][cut],
                jet_lep_dr=normalize(lep_ak8_pair[chan].i0.delta_r(lep_ak8_pair[chan].i1)),
                jet_pt=normalize(candidatejet.pt),
                lep_pt=normalize(leadinglep[chan].pt),
#                #weight=weight,
            )

            output['trigeff_miso'].fill(
                dataset=dataset,
                region=region,
                trig_pass_hadhad=trigger_hadhad[cut],
                trig_pass_hadel=trigger_hadel[cut],
                trig_pass_hadmu=trigger_hadmu[cut],
                trig_pass_ref=trigger_ref[chan][cut],
                lep_miso=normalize(leadinglep_miso[chan]),
                jet_pt=normalize(candidatejet.pt),
                lep_pt=normalize(leadinglep[chan].pt),
#                #weight=weight,
            )

        for region in regions:
            for systematic in systematics:
                fill(region, systematic)
        #    if 'GluGluHToTauTau' in dataset:
        #        for i in range(9):
        #            fill(region, 'LHEScale_%d' % i, events.LHEScaleWeight[:, i])
        #        for c in events.LHEWeight.columns[1:]:
        #            fill(region, 'LHEWeight_%s' % c, events.LHEWeight[c])

        return output

    def postprocess(self, accumulator):
        return accumulator
