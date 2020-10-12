from functools import partial
import numpy as np
from coffea import processor, hist
from uproot_methods import TLorentzVectorArray
import awkward
from .common import (
    getBosons,
    matchedBosonFlavor,
    getHTauTauDecayInfo,
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


class HtautauProcessor_Gen(processor.ProcessorABC):
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

        self._hadel_triggers = {
            '2016': [
                #'Ele35_WPTight_Gsf',
'Ele50_CaloIdVT_GsfTrkIdT_PFJet165','Ele115_CaloIdVT_GsfTrkIdT',
"Ele15_IsoVVVL_PFHT450_PFMET50",
"Ele15_IsoVVVL_PFHT600",
                'PFHT800',
                'PFHT900',
                'AK8PFJet360_TrimMass30',
                'AK8PFHT700_TrimR0p1PT0p03Mass50',
                'PFHT650_WideJetMJJ950DEtaJJ1p5',
                'PFHT650_WideJetMJJ900DEtaJJ1p5',
                'AK8DiPFJet280_200_TrimMass30_BTagCSV_p20',
                'PFJet450',
            ],
            '2017': [
                #'Ele35_WPTight_Gsf',
'Ele50_CaloIdVT_GsfTrkIdT_PFJet165','Ele115_CaloIdVT_GsfTrkIdT',
"Ele15_IsoVVVL_PFHT450_PFMET50",
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
"Ele15_IsoVVVL_PFHT450_PFMET50",
"Ele15_IsoVVVL_PFHT600",
                'AK8PFJet400_TrimMass30',
                'AK8PFJet420_TrimMass30',
                'AK8PFHT800_TrimMass50',
                'PFHT1050',
                'PFJet500',
                'AK8PFJet500',
                # 'AK8PFJet330_PFAK8BTagCSV_p17', not present in 2018D?
                'AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4',
                #'AK4PFJet30',
            ],
        }

        self._hadmu_triggers = {
            '2016': [
                'Mu50','Mu55',
"Mu15_IsoVVVL_PFHT450_PFMET50",
"Mu15_IsoVVVL_PFHT600",
                'PFHT800',
                'PFHT900',
                'AK8PFJet360_TrimMass30',
                'AK8PFHT700_TrimR0p1PT0p03Mass50',
                'PFHT650_WideJetMJJ950DEtaJJ1p5',
                'PFHT650_WideJetMJJ900DEtaJJ1p5',
                'AK8DiPFJet280_200_TrimMass30_BTagCSV_p20',
                'PFJet450',
            ],
            '2017': [
                'Mu50',#'Mu55',
"Mu15_IsoVVVL_PFHT450_PFMET50",
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
"Mu15_IsoVVVL_PFHT450_PFMET50",
"Mu15_IsoVVVL_PFHT600",
                'AK8PFJet400_TrimMass30',
                'AK8PFJet420_TrimMass30',
                'AK8PFHT800_TrimMass50',
                'PFHT1050',
                'PFJet500',
                'AK8PFJet500',
                # 'AK8PFJet330_PFAK8BTagCSV_p17', not present in 2018D?
                'AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4',
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
                'AK8DiPFJet280_200_TrimMass30_BTagCSV_p20',
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
                'AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4',
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

        jet_pt_bin = hist.Bin('jet_pt', r'Jet $p_{T}$ [GeV]', 20, 200, 1200)
        jet_eta_bin = hist.Bin('jet_eta', r'Jet $\eta$', 20, -3., 3.)
        jet_msd_bin = hist.Bin('jet_msd', r'Jet $m_{sd}$ [GeV]', 34, 40, 210.)
        oppbjet_pt_bin = hist.Bin('oppbjet_pt', r'Max opp. deepCSV-bjet $p_{T}$ [GeV]', 20, 0., 500)
        oppbtag_bin = hist.Bin('oppbtag', r'Max opp. deepCSV-b ', 20, 0., 1)
        lep_pt_bin = hist.Bin('lep_pt', r'Lepton $p_{T}$ [GeV]', 40, 0, 800)
        lep_eta_bin = hist.Bin('lep_eta', r'Lepton $\eta$', 20, -3., 3.)
        jet_lsf3_bin = hist.Bin('lsf3', r'Jet LSF$_3$', 20, 0., 1.)
        lep_jet_dr_bin = hist.Bin('lep_jet_dr', r'$\Delta R(jet,lepton)$', 40, 0., 4.)
        bos_jet_dr_bin = hist.Bin('bos_jet_dr', r'$\Delta R(boson,jet)$', 40, 0., 4.)
        met_bos_dphi_bin = hist.Bin('met_bos_dphi', r'$\Delta\phi(MET,boson)$', 40, 0., 4.)
        #lep_miso_bin = hist.Bin('miso', r'Lepton miniIso', 20, 0., 0.1)
        lep_miso_bin = hist.Bin('miso', r'Lepton miniIso', [0.,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.])
        jet_jetlep_m_bin = hist.Bin('jetlep_m', r'Jet+lepton $m$ [GeV]', 20, 0, 600.)
        jet_jetmet_m_bin = hist.Bin('jetmet_m', r'Jet+MET $m$ [GeV]', 20, 0, 600.)
        jet_jetlepmet_m_bin = hist.Bin('jetlepmet_m', r'Jet+lepton+MET $m$ [GeV]', 20, 0, 600.)
        met_pt_bin = hist.Bin('met_pt', r'MET [GeV]', 20, 0, 800)
        h_pt_bin = hist.Bin('h_pt', r'h $p_{T}$ [GeV]', 20, 200, 1200)
        bos_pt_bin = hist.Bin('bos_pt', r'boson $p_{T}$ [GeV]', 20, 200, 1200)
        genhtt_bin = hist.Bin('genhtt',r'hh,eh,mh,em,ee,mm (- for dr > 0.8)',13,-6.5,6.5)
        gentau1had_bin = hist.Bin('gentau1had',r'1pr,1pr+pi0,3pr',4,-0.5,3.5)
        gentau2had_bin = hist.Bin('gentau2had',r'1pr,1pr+pi0,3pr',4,-0.5,3.5)

        self._accumulator = processor.dict_accumulator({
            # dataset -> sumw
            'sumw': processor.defaultdict_accumulator(float),
            # dataset -> cut -> count
            'cutflow_hadhad': processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            'cutflow_hadhad_cr_mu': processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            'cutflow_hadel': processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            'cutflow_hadmu': processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            'cutflow_hadel_cr_b': processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            'cutflow_hadmu_cr_b': processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            'cutflow_hadel_cr_qcd': processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            'cutflow_hadmu_cr_qcd': processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)),
            #'btagWeight': hist.Hist('Events', hist.Cat('dataset', 'Dataset'), hist.Bin('val', 'BTag correction', 50, 0, 2)), #FIXME
            'jet_kin': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                bos_pt_bin, jet_pt_bin, jet_msd_bin, genhtt_bin,
            ),
            #'b_kin': hist.Hist(
            #    'Events',
            #    hist.Cat('dataset', 'Dataset'),
            #    hist.Cat('region', 'Region'),
            #    jet_pt_bin, oppbjet_pt_bin, oppbtag_bin, 
            #),
            'lep_kin': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                bos_pt_bin, lep_pt_bin, lep_miso_bin, genhtt_bin,
            ),
            'dr_kin': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                bos_pt_bin, lep_jet_dr_bin, bos_jet_dr_bin, genhtt_bin,
            ),
            'mass_kin': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                bos_pt_bin, jet_jetmet_m_bin, jet_jetlepmet_m_bin, genhtt_bin,
            ),
            'evt_kin': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                bos_pt_bin, met_pt_bin, met_bos_dphi_bin, genhtt_bin,
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

        trigger_hadhad = np.zeros(events.size, dtype='bool')
        for t in self._hadhad_triggers[self._year]:
            trigger_hadhad = trigger_hadhad | events.HLT[t]
        selection.add('hadhad_trigger', trigger_hadhad)

        trigger_hadmu = np.zeros(events.size, dtype='bool')
        for t in self._hadmu_triggers[self._year]:
            trigger_hadmu = trigger_hadmu | events.HLT[t]
        selection.add('hadmu_trigger', trigger_hadmu)

        trigger_hadel = np.zeros(events.size, dtype='bool')
        for t in self._hadel_triggers[self._year]:
            trigger_hadel = trigger_hadel | events.HLT[t]
        selection.add('hadel_trigger', trigger_hadel)
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
        ]#[:, :4]
        met_p4 = TLorentzVectorArray.from_ptetaphim(awkward.JaggedArray.fromiter([[v] for v in events.MET.pt]), awkward.JaggedArray.fromiter([[v] for v in np.zeros(events.size)]), awkward.JaggedArray.fromiter([[v] for v in events.MET.phi]), awkward.JaggedArray.fromiter([[v] for v in np.zeros(events.size)]))
        ak8_met_pair = candidatejets.cross(met_p4)
        ak8_met_dphi = abs(ak8_met_pair.i0.delta_phi(ak8_met_pair.i1))
        #print(np.histogram(ak8_met_dphi.fillna(0).flatten()))
        #aligned_jet = ak8_met_dphi == ak8_met_dphi.min()
        #best_jet_idx = (ak8_met_pair.i0 + aligned_jet * ak8_met_pair.i1).pt.argmax()
        best_jet_idx = ak8_met_dphi.argmin()
        #best_jet_idx = candidatejets.pt.argmax()
        candidatejet = candidatejets[best_jet_idx]
        candidatejet_rho = 2 * np.log(candidatejet.msdcorr / candidatejet.pt)
        selection.add('jetacceptance', (
            (candidatejet.pt > 300)
            & (candidatejet.msdcorr > 40.)
            & (abs(candidatejet.eta) < 2.4)
            & (candidatejet_rho > -6.)
            & (candidatejet_rho < -2.1)
        ).any())
        selection.add('jetid', (candidatejet.isTight).any())
        selection.add('n2ddt', (candidatejet.n2ddt < 0.).any())
        #print(np.histogram(candidatejet.pt.fillna(0).flatten()))

        if isRealData:
            #genflavor = candidatejet.pt.zeros_like()
            genBosonPt = candidatejet.pt.zeros_like()
            ak8_boson_dr = candidatejet.pt.zeros_like()
            met_boson_dphi = candidatejet.pt.zeros_like()
            genHTauTauDecay = candidatejet.pt.zeros_like()
            genHadTau1Decay = candidatejet.pt.zeros_like()
            genHadTau2Decay = candidatejet.pt.zeros_like()

        else:
            weights.add('genweight', events.genWeight)
            add_pileup_weight(weights, events.Pileup.nPU, self._year, dataset)
            bosons = getBosons(events).pad(1, clip=True)
            genHTauTauDecay, genHadTau1Decay, genHadTau2Decay = getHTauTauDecayInfo(events)
            gentautaudecay = awkward.JaggedArray.fromiter([[v] for v in genHTauTauDecay])
            genBosonPt = bosons.pt.fillna(0)
            add_VJets_NLOkFactor(weights, genBosonPt, self._year, dataset)
            #genflavor = matchedBosonFlavor(candidatejet, bosons)
            #add_TriggerWeight(weights, candidatejet.msdcorr, candidatejet.pt, self._year)
            #output['btagWeight'].fill(dataset=dataset, val=self._btagSF.addBtagWeight(weights, ak4_away)) #FIXME

            ak8_boson_pair = candidatejet.cross(bosons)
            met_boson_pair = met_p4.cross(bosons)
            ak8_boson_dr = ak8_boson_pair.i0.delta_r(ak8_boson_pair.i1)
            met_boson_dphi = abs(met_boson_pair.i0.delta_phi(met_boson_pair.i1))
            
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

        selection.add('met', events.MET.pt > 50.)

        goodmuon = (
            (events.Muon.pt > 25)
            & (np.abs(events.Muon.eta) < 2.4)
            #& (events.Muon.sip3d < 4)
            #& (np.abs(events.Muon.dz) < 0.1)
            #& (np.abs(events.Muon.dxy) < 0.05)
            #& (events.Muon.mediumId).astype(bool)
            & (events.Muon.highPtId).astype(bool)
        )
        ngoodmuons = goodmuon.sum()
        leadingmuon = events.Muon[goodmuon].pad(1, clip=True)

        el_loose_cuts = [(np.bitwise_and(np.right_shift(events.Electron.vidNestedWPBitmap,events.Electron.vidNestedWPBitmap.ones_like()*(3*k)),events.Electron.vidNestedWPBitmap.ones_like()*7) >= events.Electron.LOOSE) for k in range(10) if k != 7]
        #                  (MinPtCut,GsfEleSCEtaMultiRangeCut,GsfEleDEtaInSeedCut,GsfEleDPhiInCut,GsfEleFull5x5SigmaIEtaIEtaCut,GsfEleHadronicOverEMEnergyScaledCut,GsfEleEInverseMinusPInverseCut,GsfEleRelPFIsoScaledCut,GsfEleConversionVetoCut,GsfEleMissingHitsCut)
        #                   0       ,1                       ,2                  ,3              ,4                            ,5                                  ,6                             ,7                      ,8                      ,9

        elmask_loose = el_loose_cuts[0].ones_like().astype(bool)
        for m in el_loose_cuts: elmask_loose = elmask_loose & m

        goodelec = (
            (events.Electron.pt > 25)
            & (abs(events.Electron.eta) < 2.5)
            #& (events.Electron.cutBased >= events.Electron.TIGHT)
            #& (events.Electron.cutBased_HEEP).astype(bool)
            & elmask_loose
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
            #& (events.Electron.cutBased >= events.Electron.LOOSE)
            #& (events.Electron.cutBased_HEEP).astype(bool)
            & elmask_loose
        ).sum()

        #ntaus = (
        #    (events.Tau.pt > 20)
        #    & (events.Tau.idDecayMode).astype(bool)
        #    # bacon iso looser than Nano selection
        #).sum()
        ntaus = np.zeros(events.size, dtype='bool')

        lepsel = ((nmuons == 1) & (nelectrons == 0) & (ntaus == 0) & (ngoodmuons == 1)) | ((nmuons == 0) & (nelectrons == 1) & (ntaus == 0) & (ngoodelecs == 1))
        mu_p4 = TLorentzVectorArray.from_ptetaphim(leadingmuon.pt.fillna(0)*lepsel,leadingmuon.eta.fillna(0)*lepsel,leadingmuon.phi.fillna(0)*lepsel,leadingmuon.mass.fillna(0)*lepsel)
#[(goodmuon & ((nmuons == 1) & (nelectrons == 0) & (ntaus == 0) & (ngoodmuons == 1)))]
        muon_ak8_pair = mu_p4.cross(candidatejet, nested=True)
        el_p4 = TLorentzVectorArray.from_ptetaphim(leadingelec.pt.fillna(0)*lepsel,leadingelec.eta.fillna(0)*lepsel,leadingelec.phi.fillna(0)*lepsel,leadingelec.mass.fillna(0)*lepsel)
#[(goodelec & ((nmuons == 0) & (nelectrons == 1) & (ntaus == 0) & (ngoodelecs == 1)))]
        elec_ak8_pair = el_p4.cross(candidatejet, nested=True)
        #leadinglep = awkward.concatenate([mu_p4, el_p4], axis=1).pad(1, clip=True)
        leadinglep = mu_p4 + el_p4

        mu_miso = leadingmuon.miniPFRelIso_all.fillna(0)*lepsel
        el_miso = leadingelec.miniPFRelIso_all.fillna(0)*lepsel
        leadinglep_miso = mu_miso + el_miso
        leadinglep_miso = leadinglep_miso.pad(1, clip=True)

        selection.add('noleptons', (nmuons == 0) & (nelectrons == 0) & (ntaus == 0))
        selection.add('onemuon', (nmuons == 1) & (nelectrons == 0) & (ntaus == 0) & (ngoodmuons == 1))
        selection.add('oneelec', (nmuons == 0) & (nelectrons == 1) & (ntaus == 0) & (ngoodelecs == 1))
        selection.add('muonkin', (
            (leadingmuon.pt > 25.)
            & (abs(leadingmuon.eta) < 2.1)
        ).all())
        selection.add('muonDphiAK8', (
            abs(muon_ak8_pair.i0.delta_phi(muon_ak8_pair.i1)) > 2*np.pi/3
        ).all().all())
        selection.add('eleckin', (
            (leadingelec.pt > 25.)
            & (abs(leadingelec.eta) < 2.4)
        ).all())
        selection.add('elecDphiAK8', (
            abs(elec_ak8_pair.i0.delta_phi(elec_ak8_pair.i1)) > 2*np.pi/3
        ).all().all())

        lep_ak8_pair = leadinglep.cross(candidatejet)#, nested=True)
        selection.add('lepDrAK8', (
            (lep_ak8_pair.i0.delta_r(lep_ak8_pair.i1) < 0.8).all()
            #(lep_ak8_pair.i0.delta_r(lep_ak8_pair.i1) < 99.0).all()
        ))

        #selection.add('jetlsf', (
        #    (candidatejet.lsf3 > 0.7).any()
        #))

        selection.add('miniIso', (
            (leadinglep_miso < 0.1).any()
        ))
        selection.add('miniIsoInv', (
            (leadinglep_miso >= 0.1).any()
        ))

        jet_lep_p4 = lep_ak8_pair.i0 + lep_ak8_pair.i1
        met_jl_pair = met_p4.cross(jet_lep_p4)#, nested=True)
        jet_lep_met_p4 = met_jl_pair.i0 + met_jl_pair.i1
        jet_met_p4 = ak8_met_pair.i0[best_jet_idx] + ak8_met_pair.i1[best_jet_idx]

        regions = {
            'hadhad_signal': ['jetacceptance', 'hadhad_trigger', 'jetid', 'n2ddt', 'antiak4btagMediumOppHem', 'met', 'noleptons'],
            'hadhad_cr_mu': ['jetacceptance', 'hadhad_trigger', 'jetid', 'n2ddt', 'met', 'ak4btagMedium08', 'onemuon', 'muonkin', 'lepDrAK8'],#,'jetlsf'],
            'hadmu_signal': ['jetacceptance', 'hadmu_trigger', 'jetid', 'antiak4btagMediumOppHem', 'met', 'onemuon', 'muonkin', 'lepDrAK8', 'miniIso'],#, 'jetlsf'],
            'hadel_signal': ['jetacceptance', 'hadel_trigger', 'jetid', 'antiak4btagMediumOppHem', 'met', 'oneelec', 'eleckin', 'lepDrAK8', 'miniIso'],#, 'jetlsf'],
            'hadmu_cr_qcd': ['jetacceptance', 'hadmu_trigger', 'jetid', 'antiak4btagMediumOppHem', 'met', 'onemuon', 'muonkin', 'lepDrAK8', 'miniIsoInv'],#,'jetlsf'],
            'hadel_cr_qcd': ['jetacceptance', 'hadel_trigger', 'jetid', 'antiak4btagMediumOppHem', 'met', 'oneelec', 'eleckin', 'lepDrAK8', 'miniIsoInv'],#,'jetlsf'],
            'hadmu_cr_b': ['jetacceptance', 'hadmu_trigger', 'jetid', 'met', 'onemuon', 'muonkin', 'lepDrAK8', 'miniIso'],#,'jetlsf'],
            'hadel_cr_b': ['jetacceptance', 'hadel_trigger', 'jetid', 'met', 'oneelec', 'eleckin', 'lepDrAK8', 'miniIso'],#,'jetlsf'],
            #'noselection': [],
        }

        allcuts_hadel = set()
        allcuts_hadmu = set()
        allcuts_hadel_cr_b = set()
        allcuts_hadmu_cr_b = set()
        allcuts_hadel_cr_qcd = set()
        allcuts_hadmu_cr_qcd = set()
        allcuts_hadhad = set()
        allcuts_hadhad_cr_mu = set()
        output['cutflow_hadel'][dataset]['none'] += float(weights.weight().sum())
        output['cutflow_hadmu'][dataset]['none'] += float(weights.weight().sum())
        output['cutflow_hadel_cr_b'][dataset]['none'] += float(weights.weight().sum())
        output['cutflow_hadmu_cr_b'][dataset]['none'] += float(weights.weight().sum())
        output['cutflow_hadel_cr_qcd'][dataset]['none'] += float(weights.weight().sum())
        output['cutflow_hadmu_cr_qcd'][dataset]['none'] += float(weights.weight().sum())
        output['cutflow_hadhad'][dataset]['none'] += float(weights.weight().sum())
        output['cutflow_hadhad_cr_mu'][dataset]['none'] += float(weights.weight().sum())
        for cut in regions['hadel_signal']:
            allcuts_hadel.add(cut)
            output['cutflow_hadel'][dataset][cut] += float(weights.weight()[selection.all(*allcuts_hadel)].sum())
        for cut in regions['hadmu_signal']:
            allcuts_hadmu.add(cut)
            output['cutflow_hadmu'][dataset][cut] += float(weights.weight()[selection.all(*allcuts_hadmu)].sum())
        for cut in regions['hadel_cr_b']:
            allcuts_hadel_cr_b.add(cut)
            output['cutflow_hadel_cr_b'][dataset][cut] += float(weights.weight()[selection.all(*allcuts_hadel_cr_b)].sum())
        for cut in regions['hadmu_cr_b']:
            allcuts_hadmu_cr_b.add(cut)
            output['cutflow_hadmu_cr_b'][dataset][cut] += float(weights.weight()[selection.all(*allcuts_hadmu_cr_b)].sum())
        for cut in regions['hadel_cr_qcd']:
            allcuts_hadel_cr_qcd.add(cut)
            output['cutflow_hadel_cr_qcd'][dataset][cut] += float(weights.weight()[selection.all(*allcuts_hadel_cr_qcd)].sum())
        for cut in regions['hadmu_cr_qcd']:
            allcuts_hadmu_cr_qcd.add(cut)
            output['cutflow_hadmu_cr_qcd'][dataset][cut] += float(weights.weight()[selection.all(*allcuts_hadmu_cr_qcd)].sum())
        for cut in regions['hadhad_signal']:
            allcuts_hadhad.add(cut)
            output['cutflow_hadhad'][dataset][cut] += float(weights.weight()[selection.all(*allcuts_hadhad)].sum())
        for cut in regions['hadhad_cr_mu']:
            allcuts_hadhad_cr_mu.add(cut)
            output['cutflow_hadhad_cr_mu'][dataset][cut] += float(weights.weight()[selection.all(*allcuts_hadhad_cr_mu)].sum())

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
            selections = regions[region]
            cut = selection.all(*selections)
            sname = 'nominal' if systematic is None else systematic
            if wmod is None:
                weight = weights.weight(modifier=systematic)[cut]
            else:
                weight = weights.weight()[cut] * wmod[cut]

            def normalize(val):
                return val[cut].pad(1, clip=True).fillna(0).flatten()


            #print(dataset)

            output['jet_kin'].fill(
                dataset=dataset,
                region=region,
                bos_pt=normalize(genBosonPt),
                jet_pt=normalize(candidatejet.pt),
                jet_msd=normalize(candidatejet.msdcorr),
                genhtt=normalize(gentautaudecay),
                weight=weight,
            )
            output['lep_kin'].fill(
                dataset=dataset,
                region=region,
                bos_pt=normalize(genBosonPt),
                lep_pt=normalize(leadinglep.pt),
                #lep_eta=normalize(leadinglep.eta),
                #lsf3=normalize(candidatejet.lsf3),
                #lep_jet_dr=normalize(lep_ak8_pair.i0.delta_r(lep_ak8_pair.i1)),
                miso=normalize(leadinglep_miso),
                genhtt=normalize(gentautaudecay),
                weight=weight,
            )
            output['dr_kin'].fill(
                dataset=dataset,
                region=region,
                bos_pt=normalize(genBosonPt),
                lep_jet_dr=normalize(lep_ak8_pair.i0.delta_r(lep_ak8_pair.i1)),
                bos_jet_dr=normalize(ak8_boson_dr),
                genhtt=normalize(gentautaudecay),
                weight=weight,
            )
            output['mass_kin'].fill(
                dataset=dataset,
                region=region,
                bos_pt=normalize(genBosonPt),
                #jetlep_m=normalize(jet_lep_p4.mass),
                jetmet_m=normalize(jet_met_p4.mass),
                jetlepmet_m=normalize(jet_lep_met_p4.mass),
                genhtt=normalize(gentautaudecay),
                weight=weight,
            )
            output['evt_kin'].fill(
                dataset=dataset,
                region=region,
                bos_pt=normalize(genBosonPt),
                met_pt=normalize(met_p4.pt),
                met_bos_dphi=normalize(met_boson_dphi),
                #h_pt=normalize(bosons[events.GenPart.pdgId==25].pt),
                genhtt=normalize(gentautaudecay),
                weight=weight,
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
