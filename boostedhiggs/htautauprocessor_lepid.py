from functools import partial
import numpy as np
from coffea import processor, hist
from uproot_methods import TLorentzVectorArray
import awkward
from .common import (
    getBosons,
    matchedBosonFlavor,
    isOverlap,
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

class HtautauProcessor_LepID(processor.ProcessorABC):
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

        jet_pt_bin = hist.Bin('jet_pt', r'Jet $p_{T}$ [GeV]', 20, 200, 1200)
        jet_eta_bin = hist.Bin('jet_eta', r'Jet $\eta$', 20, -3., 3.)
        jet_msd_bin = hist.Bin('jet_msd', r'Jet $m_{sd}$ [GeV]', 34, 40, 210.)
        oppbjet_pt_bin = hist.Bin('oppbjet_pt', r'Max opp. deepCSV-bjet $p_{T}$ [GeV]', 20, 0., 500)
        oppbtag_bin = hist.Bin('oppbtag', r'Max opp. deepCSV-b ', 20, 0., 1)
        lep_pt_bin = hist.Bin('lep_pt', r'Lepton $p_{T}$ [GeV]', 20, 0, 400)
        lep_eta_bin = hist.Bin('lep_eta', r'Lepton $\eta$', 20, -3., 3.)
        jet_lsf3_bin = hist.Bin('lsf3', r'Jet LSF$_3$', 20, 0., 1.)
        lep_jet_dr_bin = hist.Bin('lep_jet_dr', r'$\Delta R(jet,lepton)$', 10, 0., 1.)
        #lep_miso_bin = hist.Bin('miso', r'Lepton miniIso', 20, 0., 0.1)
        lep_miso_bin = hist.Bin('miso', r'Lepton miniIso', [0.,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.])
        jet_jetlep_m_bin = hist.Bin('jetlep_m', r'Jet+lepton $m$ [GeV]', 20, 0, 600.)
        jet_jetmet_m_bin = hist.Bin('jetmet_m', r'Jet+MET $m$ [GeV]', 20, 0, 600.)
        jet_jetlepmet_m_bin = hist.Bin('jetlepmet_m', r'Jet+lepton+MET $m$ [GeV]', 20, 0, 600.)
        met_pt_bin = hist.Bin('met_pt', r'MET [GeV]', 20, 0, 800)
        h_pt_bin = hist.Bin('h_pt', r'h $p_{T}$ [GeV]', 20, 200, 1200)
        genhtt_bin = hist.Bin('genhtt',r'hh,eh,mh,em,ee,mm (- for dr > 0.8)',13,-6.5,6.5)


        self._accumulator = processor.dict_accumulator({
            # dataset -> sumw
            'sumw': processor.defaultdict_accumulator(float),
            # dataset -> cut -> count
            **{'cutflow_hadhad_%s%s_%s'%(muid,elid,vepre): processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)) for vepre in ['hv','ll','mm'] for muid in ['h','l','m','t'] for elid in ['v','l','m','t']},
            **{'cutflow_hadel_%s%s_%s'%(muid,elid,vepre): processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)) for vepre in ['hv','ll','mm'] for muid in ['h','l','m','t'] for elid in ['v','l','m','t']},
            **{'cutflow_hadmu_%s%s_%s'%(muid,elid,vepre): processor.defaultdict_accumulator(partial(processor.defaultdict_accumulator, float)) for vepre in ['hv','ll','mm'] for muid in ['h','l','m','t'] for elid in ['v','l','m','t']},
            #'btagWeight': hist.Hist('Events', hist.Cat('dataset', 'Dataset'), hist.Bin('val', 'BTag correction', 50, 0, 2)), #FIXME
            'jet_kin': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                jet_pt_bin, jet_msd_bin, genhtt_bin,
            ),
            'lep_kin': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                lep_pt_bin, lep_miso_bin, genhtt_bin,
            ),
            'pt_kin': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                lep_pt_bin, jet_pt_bin, genhtt_bin,
            ),
            'eta_kin': hist.Hist(
                'Events',
                hist.Cat('dataset', 'Dataset'),
                hist.Cat('region', 'Region'),
                lep_eta_bin, jet_eta_bin, genhtt_bin,
            ),
        })

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        dataset = events.metadata['dataset']
        isRealData = 'genWeight' not in events.columns
        selection = {
            '%s_%s'%(k,pre):processor.PackedSelection() for k in ['hadhad','hadel','hadmu'] for pre in ['hv','ll','mm','base']
        }
        weights = processor.Weights(len(events))
        output = self.accumulator.identity()
        if not isRealData:
            output['sumw'][dataset] += events.genWeight.sum()

        trigger_hadhad = np.zeros(events.size, dtype='bool')
        for t in self._hadhad_triggers[self._year]:
            trigger_hadhad = trigger_hadhad | events.HLT[t]

        trigger_hadmu = np.zeros(events.size, dtype='bool')
        for t in self._hadmu_triggers[self._year]:
            trigger_hadmu = trigger_hadmu | events.HLT[t]

        trigger_hadel = np.zeros(events.size, dtype='bool')
        for t in self._hadel_triggers[self._year]:
            trigger_hadel = trigger_hadel | events.HLT[t]
        #print(np.histogram(trigger))

        if (isRealData): overlap_removal = isOverlap(events,dataset,self._hadhad_triggers[self._year]+self._hadmu_triggers[self._year]+self._hadel_triggers[self._year])
        else: overlap_removal = np.ones(events.size, dtype='bool')

        met_filters = np.ones(events.size, dtype='bool')
        for t in self._metFilters[self._year]:
            met_filters = met_filters & events.Flag[t]

        for pre in ['hv','ll','mm','base']: selection['hadhad_%s'%pre].add('hadhad_trigger', trigger_hadhad & overlap_removal & met_filters)
        for pre in ['hv','ll','mm','base']: selection['hadmu_%s'%pre].add('hadmu_trigger', trigger_hadmu & overlap_removal & met_filters)
        for pre in ['hv','ll','mm','base']: selection['hadel_%s'%pre].add('hadel_trigger', trigger_hadel & overlap_removal & met_filters)

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
        for k in selection:
            selection[k].add('jetacceptance', (
                (candidatejet.pt > 300)
                & (candidatejet.msdcorr > 40.)
                & (abs(candidatejet.eta) < 2.4)
                & (candidatejet_rho > -6.)
                & (candidatejet_rho < -2.1)
            ).any())
            selection[k].add('jetid', (candidatejet.isTight).any())
            selection[k].add('n2ddt', (candidatejet.n2ddt < 0.).any())
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
        ak4_away = jets[(dphi > 0.8).all()]
        for k in selection:
            selection[k].add('antiak4btagMediumOppHem', ak4_opposite.btagDeepB.max() < self._btagWPs['medium'][self._year])
            #selection[k].add('ak4btagMedium08', ak4_away.btagDeepB.max() > self._btagWPs['medium'][self._year])

            selection[k].add('met', events.MET.pt > 50.)

        el_tight_cuts = [(np.bitwise_and(np.right_shift(events.Electron.vidNestedWPBitmap,events.Electron.vidNestedWPBitmap.ones_like()*(3*k)),events.Electron.vidNestedWPBitmap.ones_like()*7) >= events.Electron.TIGHT) for k in range(10) if k != 7]
        el_medium_cuts = [(np.bitwise_and(np.right_shift(events.Electron.vidNestedWPBitmap,events.Electron.vidNestedWPBitmap.ones_like()*(3*k)),events.Electron.vidNestedWPBitmap.ones_like()*7) >= events.Electron.MEDIUM) for k in range(10) if k != 7]
        el_loose_cuts = [(np.bitwise_and(np.right_shift(events.Electron.vidNestedWPBitmap,events.Electron.vidNestedWPBitmap.ones_like()*(3*k)),events.Electron.vidNestedWPBitmap.ones_like()*7) >= events.Electron.LOOSE) for k in range(10) if k != 7]
        el_veto_cuts = [(np.bitwise_and(np.right_shift(events.Electron.vidNestedWPBitmap,events.Electron.vidNestedWPBitmap.ones_like()*(3*k)),events.Electron.vidNestedWPBitmap.ones_like()*7) >= events.Electron.VETO) for k in range(10) if k != 7]
        #                  (MinPtCut,GsfEleSCEtaMultiRangeCut,GsfEleDEtaInSeedCut,GsfEleDPhiInCut,GsfEleFull5x5SigmaIEtaIEtaCut,GsfEleHadronicOverEMEnergyScaledCut,GsfEleEInverseMinusPInverseCut,GsfEleRelPFIsoScaledCut,GsfEleConversionVetoCut,GsfEleMissingHitsCut)
        #                   0       ,1                       ,2                  ,3              ,4                            ,5                                  ,6                             ,7                      ,8                      ,9
        elmask_wp90 = events.Electron.mvaFall17V2noIso_WP90
        elmask_wp80 = events.Electron.mvaFall17V2noIso_WP80

        elmask_tight = el_tight_cuts[0].ones_like().astype(bool)
        for m in el_tight_cuts: elmask_tight = elmask_tight & m
        elmask_medium = el_medium_cuts[0].ones_like().astype(bool)
        for m in el_medium_cuts: elmask_medium = elmask_medium & m
        elmask_loose = el_loose_cuts[0].ones_like().astype(bool)
        for m in el_loose_cuts: elmask_loose = elmask_loose & m
        elmask_veto = el_veto_cuts[0].ones_like().astype(bool)
        for m in el_veto_cuts: elmask_veto = elmask_veto & m

        goodmuon_highpt = (
            (events.Muon.pt > 25)
            & (np.abs(events.Muon.eta) < 2.4)
            #& (events.Muon.sip3d < 4)
            #& (np.abs(events.Muon.dz) < 0.1)
            #& (np.abs(events.Muon.dxy) < 0.05)
            #& (events.Muon.mediumId).astype(bool)
            & (events.Muon.highPtId).astype(bool)
        )
        ngoodmuons_highpt = goodmuon_highpt.sum()
        leadingmuon_highpt = events.Muon[goodmuon_highpt].pad(1, clip=True)

        goodmuon_loose = (
            (events.Muon.pt > 25)
            & (np.abs(events.Muon.eta) < 2.4)
            #& (events.Muon.sip3d < 4)
            #& (np.abs(events.Muon.dz) < 0.1)
            #& (np.abs(events.Muon.dxy) < 0.05)
            & (events.Muon.looseId).astype(bool)
            #& (events.Muon.highPtId).astype(bool)
        )
        ngoodmuons_loose = goodmuon_loose.sum()
        leadingmuon_loose = events.Muon[goodmuon_loose].pad(1, clip=True)

        goodmuon_medium = (
            (events.Muon.pt > 25)
            & (np.abs(events.Muon.eta) < 2.4)
            #& (events.Muon.sip3d < 4)
            #& (np.abs(events.Muon.dz) < 0.1)
            #& (np.abs(events.Muon.dxy) < 0.05)
            & (events.Muon.mediumId).astype(bool)
            #& (events.Muon.highPtId).astype(bool)
        )
        ngoodmuons_medium = goodmuon_medium.sum()
        leadingmuon_medium = events.Muon[goodmuon_medium].pad(1, clip=True)

        goodmuon_tight = (
            (events.Muon.pt > 25)
            & (np.abs(events.Muon.eta) < 2.4)
            #& (events.Muon.sip3d < 4)
            #& (np.abs(events.Muon.dz) < 0.1)
            #& (np.abs(events.Muon.dxy) < 0.05)
            & (events.Muon.tightId).astype(bool)
            #& (events.Muon.highPtId).astype(bool)
        )
        ngoodmuons_tight = goodmuon_tight.sum()
        leadingmuon_tight = events.Muon[goodmuon_tight].pad(1, clip=True)

        goodelec_veto = (
            (events.Electron.pt > 25)
            & (abs(events.Electron.eta) < 2.5)
            #& (events.Electron.cutBased >= events.Electron.TIGHT)
            #& (events.Electron.cutBased_HEEP).astype(bool)
            #& elmask_veto ##FIXME HACKY
            & elmask_wp80
        )
        ngoodelecs_veto = goodelec_veto.sum()
        leadingelec_veto = events.Electron[goodelec_veto].pad(1, clip=True)

        goodelec_loose = (
            (events.Electron.pt > 25)
            & (abs(events.Electron.eta) < 2.5)
            #& (events.Electron.cutBased >= events.Electron.TIGHT)
            #& (events.Electron.cutBased_HEEP).astype(bool)
            #& elmask_loose ##FIXME HACKY
            & elmask_wp90
        )
        ngoodelecs_loose = goodelec_loose.sum()
        leadingelec_loose = events.Electron[goodelec_loose].pad(1, clip=True)

        goodelec_medium = (
            (events.Electron.pt > 25)
            & (abs(events.Electron.eta) < 2.5)
            #& (events.Electron.cutBased >= events.Electron.TIGHT)
            #& (events.Electron.cutBased_HEEP).astype(bool)
            & elmask_medium
        )
        ngoodelecs_medium = goodelec_medium.sum()
        leadingelec_medium = events.Electron[goodelec_medium].pad(1, clip=True)

        goodelec_tight = (
            (events.Electron.pt > 25)
            & (abs(events.Electron.eta) < 2.5)
            #& (events.Electron.cutBased >= events.Electron.TIGHT)
            #& (events.Electron.cutBased_HEEP).astype(bool)
            & elmask_tight
        )
        ngoodelecs_tight = goodelec_tight.sum()
        leadingelec_tight = events.Electron[goodelec_tight].pad(1, clip=True)

        nmuons_loose = (
            (events.Muon.pt > 15)
            & (abs(events.Muon.eta) < 2.4)
            #& (events.Muon.pfRelIso04_all < 0.25)
            #& (np.abs(events.Muon.dz) < 0.1)
            #& (np.abs(events.Muon.dxy) < 0.05)
            & (events.Muon.looseId).astype(bool)
            #& (events.Muon.highPtId).astype(bool)
        ).sum()

        nmuons_medium = (
            (events.Muon.pt > 15)
            & (abs(events.Muon.eta) < 2.4)
            #& (events.Muon.pfRelIso04_all < 0.25)
            #& (np.abs(events.Muon.dz) < 0.1)
            #& (np.abs(events.Muon.dxy) < 0.05)
            & (events.Muon.mediumId).astype(bool)
            #& (events.Muon.highPtId).astype(bool)
        ).sum()

        nmuons_highpt = (
            (events.Muon.pt > 15)
            & (abs(events.Muon.eta) < 2.4)
            #& (events.Muon.pfRelIso04_all < 0.25)
            #& (np.abs(events.Muon.dz) < 0.1)
            #& (np.abs(events.Muon.dxy) < 0.05)
            #& (events.Muon.looseId).astype(bool)
            & (events.Muon.highPtId).astype(bool)
        ).sum()

        nelectrons_veto = (
            (events.Electron.pt > 15)
            & (abs(events.Electron.eta) < 2.5)
            & (events.Electron.cutBased >= events.Electron.VETO)
            #& (events.Electron.cutBased_HEEP).astype(bool)
            #& elmask_loose
        ).sum()

        nelectrons_loose = (
            (events.Electron.pt > 15)
            & (abs(events.Electron.eta) < 2.5)
            & (events.Electron.cutBased >= events.Electron.LOOSE)
            #& (events.Electron.cutBased_HEEP).astype(bool)
            #& elmask_loose
        ).sum()

        nelectrons_medium = (
            (events.Electron.pt > 15)
            & (abs(events.Electron.eta) < 2.5)
            & (events.Electron.cutBased >= events.Electron.MEDIUM)
            #& (events.Electron.cutBased_HEEP).astype(bool)
            #& elmask_loose
        ).sum()

        #ntaus = (
        #    (events.Tau.pt > 20)
        #    & (events.Tau.idDecayMode).astype(bool)
        #    # bacon iso looser than Nano selection
        #).sum()
        ntaus = np.zeros(events.size, dtype='bool')

        mu_p4_highpt = TLorentzVectorArray.from_ptetaphim(leadingmuon_highpt.pt.fillna(0),leadingmuon_highpt.eta.fillna(0),leadingmuon_highpt.phi.fillna(0),leadingmuon_highpt.mass.fillna(0))
        muon_ak8_pair_highpt = mu_p4_highpt.cross(candidatejet, nested=True)

        mu_p4_loose = TLorentzVectorArray.from_ptetaphim(leadingmuon_loose.pt.fillna(0),leadingmuon_loose.eta.fillna(0),leadingmuon_loose.phi.fillna(0),leadingmuon_loose.mass.fillna(0))
        muon_ak8_pair_loose = mu_p4_loose.cross(candidatejet, nested=True)

        mu_p4_medium = TLorentzVectorArray.from_ptetaphim(leadingmuon_medium.pt.fillna(0),leadingmuon_medium.eta.fillna(0),leadingmuon_medium.phi.fillna(0),leadingmuon_medium.mass.fillna(0))
        muon_ak8_pair_medium = mu_p4_medium.cross(candidatejet, nested=True)

        mu_p4_tight = TLorentzVectorArray.from_ptetaphim(leadingmuon_tight.pt.fillna(0),leadingmuon_tight.eta.fillna(0),leadingmuon_tight.phi.fillna(0),leadingmuon_tight.mass.fillna(0))
        muon_ak8_pair_tight = mu_p4_tight.cross(candidatejet, nested=True)

        el_p4_veto = TLorentzVectorArray.from_ptetaphim(leadingelec_veto.pt.fillna(0),leadingelec_veto.eta.fillna(0),leadingelec_veto.phi.fillna(0),leadingelec_veto.mass.fillna(0))
        elec_ak8_pair_veto = el_p4_veto.cross(candidatejet, nested=True)

        el_p4_loose = TLorentzVectorArray.from_ptetaphim(leadingelec_loose.pt.fillna(0),leadingelec_loose.eta.fillna(0),leadingelec_loose.phi.fillna(0),leadingelec_loose.mass.fillna(0))
        elec_ak8_pair_loose = el_p4_loose.cross(candidatejet, nested=True)

        el_p4_medium = TLorentzVectorArray.from_ptetaphim(leadingelec_medium.pt.fillna(0),leadingelec_medium.eta.fillna(0),leadingelec_medium.phi.fillna(0),leadingelec_medium.mass.fillna(0))
        elec_ak8_pair_medium = el_p4_medium.cross(candidatejet, nested=True)

        el_p4_tight = TLorentzVectorArray.from_ptetaphim(leadingelec_tight.pt.fillna(0),leadingelec_tight.eta.fillna(0),leadingelec_tight.phi.fillna(0),leadingelec_tight.mass.fillna(0))
        elec_ak8_pair_tight = el_p4_tight.cross(candidatejet, nested=True)

        leadinglep = {}

        leadinglep['hv'] = mu_p4_highpt + el_p4_veto
        leadinglep['lv'] = mu_p4_loose + el_p4_veto
        leadinglep['mv'] = mu_p4_medium + el_p4_veto
        leadinglep['tv'] = mu_p4_tight + el_p4_veto

        leadinglep['hl'] = mu_p4_highpt + el_p4_loose
        leadinglep['ll'] = mu_p4_loose + el_p4_loose
        leadinglep['ml'] = mu_p4_medium + el_p4_loose
        leadinglep['tl'] = mu_p4_tight + el_p4_loose

        leadinglep['hm'] = mu_p4_highpt + el_p4_medium
        leadinglep['lm'] = mu_p4_loose + el_p4_medium
        leadinglep['mm'] = mu_p4_medium + el_p4_medium
        leadinglep['tm'] = mu_p4_tight + el_p4_medium

        leadinglep['ht'] = mu_p4_highpt + el_p4_tight
        leadinglep['lt'] = mu_p4_loose + el_p4_tight
        leadinglep['mt'] = mu_p4_medium + el_p4_tight
        leadinglep['tt'] = mu_p4_tight + el_p4_tight

        leadinglep['base'] = leadinglep['hv'].zeros_like()

        mu_miso_highpt = leadingmuon_highpt.miniPFRelIso_all.fillna(0)
        mu_miso_loose = leadingmuon_loose.miniPFRelIso_all.fillna(0)
        mu_miso_medium = leadingmuon_medium.miniPFRelIso_all.fillna(0)
        mu_miso_tight = leadingmuon_tight.miniPFRelIso_all.fillna(0)

        el_miso_veto = leadingelec_veto.miniPFRelIso_all.fillna(0)
        el_miso_loose = leadingelec_loose.miniPFRelIso_all.fillna(0)
        el_miso_medium = leadingelec_medium.miniPFRelIso_all.fillna(0)
        el_miso_tight = leadingelec_tight.miniPFRelIso_all.fillna(0)

        leadinglep_miso = {}

        leadinglep_miso['hv'] = (mu_miso_highpt + el_miso_veto).pad(1, clip=True)
        leadinglep_miso['lv'] = (mu_miso_loose + el_miso_veto).pad(1, clip=True)
        leadinglep_miso['mv'] = (mu_miso_medium + el_miso_veto).pad(1, clip=True)
        leadinglep_miso['tv'] = (mu_miso_tight + el_miso_veto).pad(1, clip=True)

        leadinglep_miso['hl'] = (mu_miso_highpt + el_miso_loose).pad(1, clip=True)
        leadinglep_miso['ll'] = (mu_miso_loose + el_miso_loose).pad(1, clip=True)
        leadinglep_miso['ml'] = (mu_miso_medium + el_miso_loose).pad(1, clip=True)
        leadinglep_miso['tl'] = (mu_miso_tight + el_miso_loose).pad(1, clip=True)

        leadinglep_miso['hm'] = (mu_miso_highpt + el_miso_medium).pad(1, clip=True)
        leadinglep_miso['lm'] = (mu_miso_loose + el_miso_medium).pad(1, clip=True)
        leadinglep_miso['mm'] = (mu_miso_medium + el_miso_medium).pad(1, clip=True)
        leadinglep_miso['tm'] = (mu_miso_tight + el_miso_medium).pad(1, clip=True)

        leadinglep_miso['ht'] = (mu_miso_highpt + el_miso_tight).pad(1, clip=True)
        leadinglep_miso['lt'] = (mu_miso_loose + el_miso_tight).pad(1, clip=True)
        leadinglep_miso['mt'] = (mu_miso_medium + el_miso_tight).pad(1, clip=True)
        leadinglep_miso['tt'] = (mu_miso_tight + el_miso_tight).pad(1, clip=True)

        leadinglep_miso['base'] = leadinglep_miso['hv'].ones_like()

        selection['hadhad_hv'].add('noleptons_hv_hv', (nmuons_highpt == 0) & (nelectrons_veto == 0) & (ntaus == 0))
        selection['hadmu_hv'].add('onemuon_hv_hv', (nmuons_highpt <= 1) & (nelectrons_veto == 0) & (ntaus == 0) & (ngoodelecs_veto == 0) & (ngoodmuons_highpt == 1))
        selection['hadel_hv'].add('oneelec_hv_hv', (nmuons_highpt == 0) & (nelectrons_veto <= 1) & (ntaus == 0) & (ngoodmuons_highpt == 0) & (ngoodelecs_veto == 1))

        selection['hadhad_hv'].add('noleptons_lv_hv', (nmuons_highpt == 0) & (nelectrons_veto == 0) & (ntaus == 0))
        selection['hadmu_hv'].add('onemuon_lv_hv', (nmuons_highpt <= 1) & (nelectrons_veto == 0) & (ntaus == 0) & (ngoodelecs_veto == 0) & (ngoodmuons_loose == 1))
        selection['hadel_hv'].add('oneelec_lv_hv', (nmuons_highpt == 0) & (nelectrons_veto <= 1) & (ntaus == 0) & (ngoodmuons_loose == 0) & (ngoodelecs_veto == 1))

        selection['hadhad_hv'].add('noleptons_mv_hv', (nmuons_highpt == 0) & (nelectrons_veto == 0) & (ntaus == 0))
        selection['hadmu_hv'].add('onemuon_mv_hv', (nmuons_highpt <= 1) & (nelectrons_veto == 0) & (ntaus == 0) & (ngoodelecs_veto == 0) & (ngoodmuons_medium == 1))
        selection['hadel_hv'].add('oneelec_mv_hv', (nmuons_highpt == 0) & (nelectrons_veto <= 1) & (ntaus == 0) & (ngoodmuons_medium == 0) & (ngoodelecs_veto == 1))

        selection['hadhad_hv'].add('noleptons_tv_hv', (nmuons_highpt == 0) & (nelectrons_veto == 0) & (ntaus == 0))
        selection['hadmu_hv'].add('onemuon_tv_hv', (nmuons_highpt <= 1) & (nelectrons_veto == 0) & (ntaus == 0) & (ngoodelecs_veto == 0) & (ngoodmuons_tight == 1))
        selection['hadel_hv'].add('oneelec_tv_hv', (nmuons_highpt == 0) & (nelectrons_veto <= 1) & (ntaus == 0) & (ngoodmuons_tight == 0) & (ngoodelecs_veto == 1))

        selection['hadhad_ll'].add('noleptons_hv_ll', (nmuons_loose == 0) & (nelectrons_loose == 0) & (ntaus == 0))
        selection['hadmu_ll'].add('onemuon_hv_ll', (nmuons_loose <= 1) & (nelectrons_loose == 0) & (ntaus == 0) & (ngoodelecs_veto == 0) & (ngoodmuons_highpt == 1))
        selection['hadel_ll'].add('oneelec_hv_ll', (nmuons_loose == 0) & (nelectrons_loose <= 1) & (ntaus == 0) & (ngoodmuons_highpt == 0) & (ngoodelecs_veto == 1))

        selection['hadhad_ll'].add('noleptons_lv_ll', (nmuons_loose == 0) & (nelectrons_loose == 0) & (ntaus == 0))
        selection['hadmu_ll'].add('onemuon_lv_ll', (nmuons_loose <= 1) & (nelectrons_loose == 0) & (ntaus == 0) & (ngoodelecs_veto == 0) & (ngoodmuons_loose == 1))
        selection['hadel_ll'].add('oneelec_lv_ll', (nmuons_loose == 0) & (nelectrons_loose <= 1) & (ntaus == 0) & (ngoodmuons_loose == 0) & (ngoodelecs_veto == 1))

        selection['hadhad_ll'].add('noleptons_mv_ll', (nmuons_loose == 0) & (nelectrons_loose == 0) & (ntaus == 0))
        selection['hadmu_ll'].add('onemuon_mv_ll', (nmuons_loose <= 1) & (nelectrons_loose == 0) & (ntaus == 0) & (ngoodelecs_veto == 0) & (ngoodmuons_medium == 1))
        selection['hadel_ll'].add('oneelec_mv_ll', (nmuons_loose == 0) & (nelectrons_loose <= 1) & (ntaus == 0) & (ngoodmuons_medium == 0) & (ngoodelecs_veto == 1))

        selection['hadhad_ll'].add('noleptons_tv_ll', (nmuons_loose == 0) & (nelectrons_loose == 0) & (ntaus == 0))
        selection['hadmu_ll'].add('onemuon_tv_ll', (nmuons_loose <= 1) & (nelectrons_loose == 0) & (ntaus == 0) & (ngoodelecs_veto == 0) & (ngoodmuons_tight == 1))
        selection['hadel_ll'].add('oneelec_tv_ll', (nmuons_loose == 0) & (nelectrons_loose <= 1) & (ntaus == 0) & (ngoodmuons_tight == 0) & (ngoodelecs_veto == 1))

        selection['hadhad_mm'].add('noleptons_hv_mm', (nmuons_medium == 0) & (nelectrons_medium == 0) & (ntaus == 0))
        selection['hadmu_mm'].add('onemuon_hv_mm', (nmuons_medium <= 1) & (nelectrons_medium == 0) & (ntaus == 0) & (ngoodelecs_veto == 0) & (ngoodmuons_highpt == 1))
        selection['hadel_mm'].add('oneelec_hv_mm', (nmuons_medium == 0) & (nelectrons_medium <= 1) & (ntaus == 0) & (ngoodmuons_highpt == 0) & (ngoodelecs_veto == 1))

        selection['hadhad_mm'].add('noleptons_lv_mm', (nmuons_medium == 0) & (nelectrons_medium == 0) & (ntaus == 0))
        selection['hadmu_mm'].add('onemuon_lv_mm', (nmuons_medium <= 1) & (nelectrons_medium == 0) & (ntaus == 0) & (ngoodelecs_veto == 0) & (ngoodmuons_loose == 1))
        selection['hadel_mm'].add('oneelec_lv_mm', (nmuons_medium == 0) & (nelectrons_medium <= 1) & (ntaus == 0) & (ngoodmuons_loose == 0) & (ngoodelecs_veto == 1))

        selection['hadhad_mm'].add('noleptons_mv_mm', (nmuons_medium == 0) & (nelectrons_medium == 0) & (ntaus == 0))
        selection['hadmu_mm'].add('onemuon_mv_mm', (nmuons_medium <= 1) & (nelectrons_medium == 0) & (ntaus == 0) & (ngoodelecs_veto == 0) & (ngoodmuons_medium == 1))
        selection['hadel_mm'].add('oneelec_mv_mm', (nmuons_medium == 0) & (nelectrons_medium <= 1) & (ntaus == 0) & (ngoodmuons_medium == 0) & (ngoodelecs_veto == 1))

        selection['hadhad_mm'].add('noleptons_tv_mm', (nmuons_medium == 0) & (nelectrons_medium == 0) & (ntaus == 0))
        selection['hadmu_mm'].add('onemuon_tv_mm', (nmuons_medium <= 1) & (nelectrons_medium == 0) & (ntaus == 0) & (ngoodelecs_veto == 0) & (ngoodmuons_tight == 1))
        selection['hadel_mm'].add('oneelec_tv_mm', (nmuons_medium == 0) & (nelectrons_medium <= 1) & (ntaus == 0) & (ngoodmuons_tight == 0) & (ngoodelecs_veto == 1))

        selection['hadhad_hv'].add('noleptons_hl_hv', (nmuons_highpt == 0) & (nelectrons_veto == 0) & (ntaus == 0))
        selection['hadmu_hv'].add('onemuon_hl_hv', (nmuons_highpt <= 1) & (nelectrons_veto == 0) & (ntaus == 0) & (ngoodelecs_loose == 0) & (ngoodmuons_highpt == 1))
        selection['hadel_hv'].add('oneelec_hl_hv', (nmuons_highpt == 0) & (nelectrons_veto <= 1) & (ntaus == 0) & (ngoodmuons_highpt == 0) & (ngoodelecs_loose == 1))

        selection['hadhad_hv'].add('noleptons_ll_hv', (nmuons_highpt == 0) & (nelectrons_veto == 0) & (ntaus == 0))
        selection['hadmu_hv'].add('onemuon_ll_hv', (nmuons_highpt <= 1) & (nelectrons_veto == 0) & (ntaus == 0) & (ngoodelecs_loose == 0) & (ngoodmuons_loose == 1))
        selection['hadel_hv'].add('oneelec_ll_hv', (nmuons_highpt == 0) & (nelectrons_veto <= 1) & (ntaus == 0) & (ngoodmuons_loose == 0) & (ngoodelecs_loose == 1))

        selection['hadhad_hv'].add('noleptons_ml_hv', (nmuons_highpt == 0) & (nelectrons_veto == 0) & (ntaus == 0))
        selection['hadmu_hv'].add('onemuon_ml_hv', (nmuons_highpt <= 1) & (nelectrons_veto == 0) & (ntaus == 0) & (ngoodelecs_loose == 0) & (ngoodmuons_medium == 1))
        selection['hadel_hv'].add('oneelec_ml_hv', (nmuons_highpt == 0) & (nelectrons_veto <= 1) & (ntaus == 0) & (ngoodmuons_medium == 0) & (ngoodelecs_loose == 1))

        selection['hadhad_hv'].add('noleptons_tl_hv', (nmuons_highpt == 0) & (nelectrons_veto == 0) & (ntaus == 0))
        selection['hadmu_hv'].add('onemuon_tl_hv', (nmuons_highpt <= 1) & (nelectrons_veto == 0) & (ntaus == 0) & (ngoodelecs_loose == 0) & (ngoodmuons_tight == 1))
        selection['hadel_hv'].add('oneelec_tl_hv', (nmuons_highpt == 0) & (nelectrons_veto <= 1) & (ntaus == 0) & (ngoodmuons_tight == 0) & (ngoodelecs_loose == 1))

        selection['hadhad_ll'].add('noleptons_hl_ll', (nmuons_loose == 0) & (nelectrons_loose == 0) & (ntaus == 0))
        selection['hadmu_ll'].add('onemuon_hl_ll', (nmuons_loose <= 1) & (nelectrons_loose == 0) & (ntaus == 0) & (ngoodelecs_loose == 0) & (ngoodmuons_highpt == 1))
        selection['hadel_ll'].add('oneelec_hl_ll', (nmuons_loose == 0) & (nelectrons_loose <= 1) & (ntaus == 0) & (ngoodmuons_highpt == 0) & (ngoodelecs_loose == 1))

        selection['hadhad_ll'].add('noleptons_ll_ll', (nmuons_loose == 0) & (nelectrons_loose == 0) & (ntaus == 0))
        selection['hadmu_ll'].add('onemuon_ll_ll', (nmuons_loose <= 1) & (nelectrons_loose == 0) & (ntaus == 0) & (ngoodelecs_loose == 0) & (ngoodmuons_loose == 1))
        selection['hadel_ll'].add('oneelec_ll_ll', (nmuons_loose == 0) & (nelectrons_loose <= 1) & (ntaus == 0) & (ngoodmuons_loose == 0) & (ngoodelecs_loose == 1))

        selection['hadhad_ll'].add('noleptons_ml_ll', (nmuons_loose == 0) & (nelectrons_loose == 0) & (ntaus == 0))
        selection['hadmu_ll'].add('onemuon_ml_ll', (nmuons_loose <= 1) & (nelectrons_loose == 0) & (ntaus == 0) & (ngoodelecs_loose == 0) & (ngoodmuons_medium == 1))
        selection['hadel_ll'].add('oneelec_ml_ll', (nmuons_loose == 0) & (nelectrons_loose <= 1) & (ntaus == 0) & (ngoodmuons_medium == 0) & (ngoodelecs_loose == 1))

        selection['hadhad_ll'].add('noleptons_tl_ll', (nmuons_loose == 0) & (nelectrons_loose == 0) & (ntaus == 0))
        selection['hadmu_ll'].add('onemuon_tl_ll', (nmuons_loose <= 1) & (nelectrons_loose == 0) & (ntaus == 0) & (ngoodelecs_loose == 0) & (ngoodmuons_tight == 1))
        selection['hadel_ll'].add('oneelec_tl_ll', (nmuons_loose == 0) & (nelectrons_loose <= 1) & (ntaus == 0) & (ngoodmuons_tight == 0) & (ngoodelecs_loose == 1))

        selection['hadhad_mm'].add('noleptons_hl_mm', (nmuons_medium == 0) & (nelectrons_medium == 0) & (ntaus == 0))
        selection['hadmu_mm'].add('onemuon_hl_mm', (nmuons_medium <= 1) & (nelectrons_medium == 0) & (ntaus == 0) & (ngoodelecs_loose == 0) & (ngoodmuons_highpt == 1))
        selection['hadel_mm'].add('oneelec_hl_mm', (nmuons_medium == 0) & (nelectrons_medium <= 1) & (ntaus == 0) & (ngoodmuons_highpt == 0) & (ngoodelecs_loose == 1))

        selection['hadhad_mm'].add('noleptons_ll_mm', (nmuons_medium == 0) & (nelectrons_medium == 0) & (ntaus == 0))
        selection['hadmu_mm'].add('onemuon_ll_mm', (nmuons_medium <= 1) & (nelectrons_medium == 0) & (ntaus == 0) & (ngoodelecs_loose == 0) & (ngoodmuons_loose == 1))
        selection['hadel_mm'].add('oneelec_ll_mm', (nmuons_medium == 0) & (nelectrons_medium <= 1) & (ntaus == 0) & (ngoodmuons_loose == 0) & (ngoodelecs_loose == 1))

        selection['hadhad_mm'].add('noleptons_ml_mm', (nmuons_medium == 0) & (nelectrons_medium == 0) & (ntaus == 0))
        selection['hadmu_mm'].add('onemuon_ml_mm', (nmuons_medium <= 1) & (nelectrons_medium == 0) & (ntaus == 0) & (ngoodelecs_loose == 0) & (ngoodmuons_medium == 1))
        selection['hadel_mm'].add('oneelec_ml_mm', (nmuons_medium == 0) & (nelectrons_medium <= 1) & (ntaus == 0) & (ngoodmuons_medium == 0) & (ngoodelecs_loose == 1))

        selection['hadhad_mm'].add('noleptons_tl_mm', (nmuons_medium == 0) & (nelectrons_medium == 0) & (ntaus == 0))
        selection['hadmu_mm'].add('onemuon_tl_mm', (nmuons_medium <= 1) & (nelectrons_medium == 0) & (ntaus == 0) & (ngoodelecs_loose == 0) & (ngoodmuons_tight == 1))
        selection['hadel_mm'].add('oneelec_tl_mm', (nmuons_medium == 0) & (nelectrons_medium <= 1) & (ntaus == 0) & (ngoodmuons_tight == 0) & (ngoodelecs_loose == 1))

        selection['hadhad_hv'].add('noleptons_hm_hv', (nmuons_highpt == 0) & (nelectrons_veto == 0) & (ntaus == 0))
        selection['hadmu_hv'].add('onemuon_hm_hv', (nmuons_highpt <= 1) & (nelectrons_veto == 0) & (ntaus == 0) & (ngoodelecs_medium == 0) & (ngoodmuons_highpt == 1))
        selection['hadel_hv'].add('oneelec_hm_hv', (nmuons_highpt == 0) & (nelectrons_veto <= 1) & (ntaus == 0) & (ngoodmuons_highpt == 0) & (ngoodelecs_medium == 1))

        selection['hadhad_hv'].add('noleptons_lm_hv', (nmuons_highpt == 0) & (nelectrons_veto == 0) & (ntaus == 0))
        selection['hadmu_hv'].add('onemuon_lm_hv', (nmuons_highpt <= 1) & (nelectrons_veto == 0) & (ntaus == 0) & (ngoodelecs_medium == 0) & (ngoodmuons_loose == 1))
        selection['hadel_hv'].add('oneelec_lm_hv', (nmuons_highpt == 0) & (nelectrons_veto <= 1) & (ntaus == 0) & (ngoodmuons_loose == 0) & (ngoodelecs_medium == 1))

        selection['hadhad_hv'].add('noleptons_mm_hv', (nmuons_highpt == 0) & (nelectrons_veto == 0) & (ntaus == 0))
        selection['hadmu_hv'].add('onemuon_mm_hv', (nmuons_highpt <= 1) & (nelectrons_veto == 0) & (ntaus == 0) & (ngoodelecs_medium == 0) & (ngoodmuons_medium == 1))
        selection['hadel_hv'].add('oneelec_mm_hv', (nmuons_highpt == 0) & (nelectrons_veto <= 1) & (ntaus == 0) & (ngoodmuons_medium == 0) & (ngoodelecs_medium == 1))

        selection['hadhad_hv'].add('noleptons_tm_hv', (nmuons_highpt == 0) & (nelectrons_veto == 0) & (ntaus == 0))
        selection['hadmu_hv'].add('onemuon_tm_hv', (nmuons_highpt <= 1) & (nelectrons_veto == 0) & (ntaus == 0) & (ngoodelecs_medium == 0) & (ngoodmuons_tight == 1))
        selection['hadel_hv'].add('oneelec_tm_hv', (nmuons_highpt == 0) & (nelectrons_veto <= 1) & (ntaus == 0) & (ngoodmuons_tight == 0) & (ngoodelecs_medium == 1))

        selection['hadhad_ll'].add('noleptons_hm_ll', (nmuons_loose == 0) & (nelectrons_loose == 0) & (ntaus == 0))
        selection['hadmu_ll'].add('onemuon_hm_ll', (nmuons_loose <= 1) & (nelectrons_loose == 0) & (ntaus == 0) & (ngoodelecs_medium == 0) & (ngoodmuons_highpt == 1))
        selection['hadel_ll'].add('oneelec_hm_ll', (nmuons_loose == 0) & (nelectrons_loose <= 1) & (ntaus == 0) & (ngoodmuons_highpt == 0) & (ngoodelecs_medium == 1))

        selection['hadhad_ll'].add('noleptons_lm_ll', (nmuons_loose == 0) & (nelectrons_loose == 0) & (ntaus == 0))
        selection['hadmu_ll'].add('onemuon_lm_ll', (nmuons_loose <= 1) & (nelectrons_loose == 0) & (ntaus == 0) & (ngoodelecs_medium == 0) & (ngoodmuons_loose == 1))
        selection['hadel_ll'].add('oneelec_lm_ll', (nmuons_loose == 0) & (nelectrons_loose <= 1) & (ntaus == 0) & (ngoodmuons_loose == 0) & (ngoodelecs_medium == 1))

        selection['hadhad_ll'].add('noleptons_mm_ll', (nmuons_loose == 0) & (nelectrons_loose == 0) & (ntaus == 0))
        selection['hadmu_ll'].add('onemuon_mm_ll', (nmuons_loose <= 1) & (nelectrons_loose == 0) & (ntaus == 0) & (ngoodelecs_medium == 0) & (ngoodmuons_medium == 1))
        selection['hadel_ll'].add('oneelec_mm_ll', (nmuons_loose == 0) & (nelectrons_loose <= 1) & (ntaus == 0) & (ngoodmuons_medium == 0) & (ngoodelecs_medium == 1))

        selection['hadhad_ll'].add('noleptons_tm_ll', (nmuons_loose == 0) & (nelectrons_loose == 0) & (ntaus == 0))
        selection['hadmu_ll'].add('onemuon_tm_ll', (nmuons_loose <= 1) & (nelectrons_loose == 0) & (ntaus == 0) & (ngoodelecs_medium == 0) & (ngoodmuons_tight == 1))
        selection['hadel_ll'].add('oneelec_tm_ll', (nmuons_loose == 0) & (nelectrons_loose <= 1) & (ntaus == 0) & (ngoodmuons_tight == 0) & (ngoodelecs_medium == 1))

        selection['hadhad_mm'].add('noleptons_hm_mm', (nmuons_medium == 0) & (nelectrons_medium == 0) & (ntaus == 0))
        selection['hadmu_mm'].add('onemuon_hm_mm', (nmuons_medium <= 1) & (nelectrons_medium == 0) & (ntaus == 0) & (ngoodelecs_medium == 0) & (ngoodmuons_highpt == 1))
        selection['hadel_mm'].add('oneelec_hm_mm', (nmuons_medium == 0) & (nelectrons_medium <= 1) & (ntaus == 0) & (ngoodmuons_highpt == 0) & (ngoodelecs_medium == 1))

        selection['hadhad_mm'].add('noleptons_lm_mm', (nmuons_medium == 0) & (nelectrons_medium == 0) & (ntaus == 0))
        selection['hadmu_mm'].add('onemuon_lm_mm', (nmuons_medium <= 1) & (nelectrons_medium == 0) & (ntaus == 0) & (ngoodelecs_medium == 0) & (ngoodmuons_loose == 1))
        selection['hadel_mm'].add('oneelec_lm_mm', (nmuons_medium == 0) & (nelectrons_medium <= 1) & (ntaus == 0) & (ngoodmuons_loose == 0) & (ngoodelecs_medium == 1))

        selection['hadhad_mm'].add('noleptons_mm_mm', (nmuons_medium == 0) & (nelectrons_medium == 0) & (ntaus == 0))
        selection['hadmu_mm'].add('onemuon_mm_mm', (nmuons_medium <= 1) & (nelectrons_medium == 0) & (ntaus == 0) & (ngoodelecs_medium == 0) & (ngoodmuons_medium == 1))
        selection['hadel_mm'].add('oneelec_mm_mm', (nmuons_medium == 0) & (nelectrons_medium <= 1) & (ntaus == 0) & (ngoodmuons_medium == 0) & (ngoodelecs_medium == 1))

        selection['hadhad_mm'].add('noleptons_tm_mm', (nmuons_medium == 0) & (nelectrons_medium == 0) & (ntaus == 0))
        selection['hadmu_mm'].add('onemuon_tm_mm', (nmuons_medium <= 1) & (nelectrons_medium == 0) & (ntaus == 0) & (ngoodelecs_medium == 0) & (ngoodmuons_tight == 1))
        selection['hadel_mm'].add('oneelec_tm_mm', (nmuons_medium == 0) & (nelectrons_medium <= 1) & (ntaus == 0) & (ngoodmuons_tight == 0) & (ngoodelecs_medium == 1))

        selection['hadhad_hv'].add('noleptons_ht_hv', (nmuons_highpt == 0) & (nelectrons_veto == 0) & (ntaus == 0))
        selection['hadmu_hv'].add('onemuon_ht_hv', (nmuons_highpt <= 1) & (nelectrons_veto == 0) & (ntaus == 0) & (ngoodelecs_tight == 0) & (ngoodmuons_highpt == 1))
        selection['hadel_hv'].add('oneelec_ht_hv', (nmuons_highpt == 0) & (nelectrons_veto <= 1) & (ntaus == 0) & (ngoodmuons_highpt == 0) & (ngoodelecs_tight == 1))

        selection['hadhad_hv'].add('noleptons_lt_hv', (nmuons_highpt == 0) & (nelectrons_veto == 0) & (ntaus == 0))
        selection['hadmu_hv'].add('onemuon_lt_hv', (nmuons_highpt <= 1) & (nelectrons_veto == 0) & (ntaus == 0) & (ngoodelecs_tight == 0) & (ngoodmuons_loose == 1))
        selection['hadel_hv'].add('oneelec_lt_hv', (nmuons_highpt == 0) & (nelectrons_veto <= 1) & (ntaus == 0) & (ngoodmuons_loose == 0) & (ngoodelecs_tight == 1))

        selection['hadhad_hv'].add('noleptons_mt_hv', (nmuons_highpt == 0) & (nelectrons_veto == 0) & (ntaus == 0))
        selection['hadmu_hv'].add('onemuon_mt_hv', (nmuons_highpt <= 1) & (nelectrons_veto == 0) & (ntaus == 0) & (ngoodelecs_tight == 0) & (ngoodmuons_medium == 1))
        selection['hadel_hv'].add('oneelec_mt_hv', (nmuons_highpt == 0) & (nelectrons_veto <= 1) & (ntaus == 0) & (ngoodmuons_medium == 0) & (ngoodelecs_tight == 1))

        selection['hadhad_hv'].add('noleptons_tt_hv', (nmuons_highpt == 0) & (nelectrons_veto == 0) & (ntaus == 0))
        selection['hadmu_hv'].add('onemuon_tt_hv', (nmuons_highpt <= 1) & (nelectrons_veto == 0) & (ntaus == 0) & (ngoodelecs_tight == 0) & (ngoodmuons_tight == 1))
        selection['hadel_hv'].add('oneelec_tt_hv', (nmuons_highpt == 0) & (nelectrons_veto <= 1) & (ntaus == 0) & (ngoodmuons_tight == 0) & (ngoodelecs_tight == 1))

        selection['hadhad_ll'].add('noleptons_ht_ll', (nmuons_loose == 0) & (nelectrons_loose == 0) & (ntaus == 0))
        selection['hadmu_ll'].add('onemuon_ht_ll', (nmuons_loose <= 1) & (nelectrons_loose == 0) & (ntaus == 0) & (ngoodelecs_tight == 0) & (ngoodmuons_highpt == 1))
        selection['hadel_ll'].add('oneelec_ht_ll', (nmuons_loose == 0) & (nelectrons_loose <= 1) & (ntaus == 0) & (ngoodmuons_highpt == 0) & (ngoodelecs_tight == 1))

        selection['hadhad_ll'].add('noleptons_lt_ll', (nmuons_loose == 0) & (nelectrons_loose == 0) & (ntaus == 0))
        selection['hadmu_ll'].add('onemuon_lt_ll', (nmuons_loose <= 1) & (nelectrons_loose == 0) & (ntaus == 0) & (ngoodelecs_tight == 0) & (ngoodmuons_loose == 1))
        selection['hadel_ll'].add('oneelec_lt_ll', (nmuons_loose == 0) & (nelectrons_loose <= 1) & (ntaus == 0) & (ngoodmuons_loose == 0) & (ngoodelecs_tight == 1))

        selection['hadhad_ll'].add('noleptons_mt_ll', (nmuons_loose == 0) & (nelectrons_loose == 0) & (ntaus == 0))
        selection['hadmu_ll'].add('onemuon_mt_ll', (nmuons_loose <= 1) & (nelectrons_loose == 0) & (ntaus == 0) & (ngoodelecs_tight == 0) & (ngoodmuons_medium == 1))
        selection['hadel_ll'].add('oneelec_mt_ll', (nmuons_loose == 0) & (nelectrons_loose <= 1) & (ntaus == 0) & (ngoodmuons_medium == 0) & (ngoodelecs_tight == 1))

        selection['hadhad_ll'].add('noleptons_tt_ll', (nmuons_loose == 0) & (nelectrons_loose == 0) & (ntaus == 0))
        selection['hadmu_ll'].add('onemuon_tt_ll', (nmuons_loose <= 1) & (nelectrons_loose == 0) & (ntaus == 0) & (ngoodelecs_tight == 0) & (ngoodmuons_tight == 1))
        selection['hadel_ll'].add('oneelec_tt_ll', (nmuons_loose == 0) & (nelectrons_loose <= 1) & (ntaus == 0) & (ngoodmuons_tight == 0) & (ngoodelecs_tight == 1))

        selection['hadhad_mm'].add('noleptons_ht_mm', (nmuons_medium == 0) & (nelectrons_medium == 0) & (ntaus == 0))
        selection['hadmu_mm'].add('onemuon_ht_mm', (nmuons_medium <= 1) & (nelectrons_medium == 0) & (ntaus == 0) & (ngoodelecs_tight == 0) & (ngoodmuons_highpt == 1))
        selection['hadel_mm'].add('oneelec_ht_mm', (nmuons_medium == 0) & (nelectrons_medium <= 1) & (ntaus == 0) & (ngoodmuons_highpt == 0) & (ngoodelecs_tight == 1))

        selection['hadhad_mm'].add('noleptons_lt_mm', (nmuons_medium == 0) & (nelectrons_medium == 0) & (ntaus == 0))
        selection['hadmu_mm'].add('onemuon_lt_mm', (nmuons_medium <= 1) & (nelectrons_medium == 0) & (ntaus == 0) & (ngoodelecs_tight == 0) & (ngoodmuons_loose == 1))
        selection['hadel_mm'].add('oneelec_lt_mm', (nmuons_medium == 0) & (nelectrons_medium <= 1) & (ntaus == 0) & (ngoodmuons_loose == 0) & (ngoodelecs_tight == 1))

        selection['hadhad_mm'].add('noleptons_mt_mm', (nmuons_medium == 0) & (nelectrons_medium == 0) & (ntaus == 0))
        selection['hadmu_mm'].add('onemuon_mt_mm', (nmuons_medium <= 1) & (nelectrons_medium == 0) & (ntaus == 0) & (ngoodelecs_tight == 0) & (ngoodmuons_medium == 1))
        selection['hadel_mm'].add('oneelec_mt_mm', (nmuons_medium == 0) & (nelectrons_medium <= 1) & (ntaus == 0) & (ngoodmuons_medium == 0) & (ngoodelecs_tight == 1))

        selection['hadhad_mm'].add('noleptons_tt_mm', (nmuons_medium == 0) & (nelectrons_medium == 0) & (ntaus == 0))
        selection['hadmu_mm'].add('onemuon_tt_mm', (nmuons_medium <= 1) & (nelectrons_medium == 0) & (ntaus == 0) & (ngoodelecs_tight == 0) & (ngoodmuons_tight == 1))
        selection['hadel_mm'].add('oneelec_tt_mm', (nmuons_medium == 0) & (nelectrons_medium <= 1) & (ntaus == 0) & (ngoodmuons_tight == 0) & (ngoodelecs_tight == 1))

        lep_ak8_pair = {}

        for k in selection:
            lep_ak8_pair['hv'] = leadinglep['hv'].cross(candidatejet)#, nested=True)
            selection[k].add('lepDrAK8_hv', (
                (lep_ak8_pair['hv'].i0.delta_r(lep_ak8_pair['hv'].i1) < 0.8).all()
            ))
            #selection.add('miniIso_hv', (
                #(leadinglep_miso['hv'] < 0.1).any()
            #))
            #selection.add('miniIsoInv_hv', (
                #(leadinglep_miso['hv'] >= 0.1).any()
            #))
    
            lep_ak8_pair['lv'] = leadinglep['lv'].cross(candidatejet)#, nested=True)
            selection[k].add('lepDrAK8_lv', (
                (lep_ak8_pair['lv'].i0.delta_r(lep_ak8_pair['lv'].i1) < 0.8).all()
            ))
            #selection.add('miniIso_lv', (
                #(leadinglep_miso['lv'] < 0.1).any()
            #))
            #selection.add('miniIsoInv_lv', (
                #(leadinglep_miso['lv'] >= 0.1).any()
            #))
    
            lep_ak8_pair['mv'] = leadinglep['mv'].cross(candidatejet)#, nested=True)
            selection[k].add('lepDrAK8_mv', (
                (lep_ak8_pair['mv'].i0.delta_r(lep_ak8_pair['mv'].i1) < 0.8).all()
            ))
            #selection.add('miniIso_mv', (
                #(leadinglep_miso['mv'] < 0.1).any()
            #))
            #selection.add('miniIsoInv_mv', (
                #(leadinglep_miso['mv'] >= 0.1).any()
            #))
    
            lep_ak8_pair['tv'] = leadinglep['tv'].cross(candidatejet)#, nested=True)
            selection[k].add('lepDrAK8_tv', (
                (lep_ak8_pair['tv'].i0.delta_r(lep_ak8_pair['tv'].i1) < 0.8).all()
            ))
            #selection.add('miniIso_tv', (
                #(leadinglep_miso['tv'] < 0.1).any()
            #))
            #selection.add('miniIsoInv_tv', (
                #(leadinglep_miso['tv'] >= 0.1).any()
            #))
    
            lep_ak8_pair['hl'] = leadinglep['hl'].cross(candidatejet)#, nested=True)
            selection[k].add('lepDrAK8_hl', (
                (lep_ak8_pair['hl'].i0.delta_r(lep_ak8_pair['hl'].i1) < 0.8).all()
            ))
            #selection.add('miniIso_hl', (
                #(leadinglep_miso['hl'] < 0.1).any()
            #))
            #selection.add('miniIsoInv_hl', (
                #(leadinglep_miso['hl'] >= 0.1).any()
            #))
    
            lep_ak8_pair['ll'] = leadinglep['ll'].cross(candidatejet)#, nested=True)
            selection[k].add('lepDrAK8_ll', (
                (lep_ak8_pair['ll'].i0.delta_r(lep_ak8_pair['ll'].i1) < 0.8).all()
            ))
            #selection.add('miniIso_ll', (
                #(leadinglep_miso['ll'] < 0.1).any()
            #))
            #selection.add('miniIsoInv_ll', (
                #(leadinglep_miso['ll'] >= 0.1).any()
            #))
    
            lep_ak8_pair['ml'] = leadinglep['ml'].cross(candidatejet)#, nested=True)
            selection[k].add('lepDrAK8_ml', (
                (lep_ak8_pair['ml'].i0.delta_r(lep_ak8_pair['ml'].i1) < 0.8).all()
            ))
            #selection.add('miniIso_ml', (
                #(leadinglep_miso['ml'] < 0.1).any()
            #))
            #selection.add('miniIsoInv_ml', (
                #(leadinglep_miso['ml'] >= 0.1).any()
            #))
    
            lep_ak8_pair['tl'] = leadinglep['tl'].cross(candidatejet)#, nested=True)
            selection[k].add('lepDrAK8_tl', (
                (lep_ak8_pair['tl'].i0.delta_r(lep_ak8_pair['tl'].i1) < 0.8).all()
            ))
            #selection.add('miniIso_tl', (
                #(leadinglep_miso['tl'] < 0.1).any()
            #))
            #selection.add('miniIsoInv_tl', (
                #(leadinglep_miso['tl'] >= 0.1).any()
            #))
    
            lep_ak8_pair['hm'] = leadinglep['hm'].cross(candidatejet)#, nested=True)
            selection[k].add('lepDrAK8_hm', (
                (lep_ak8_pair['hm'].i0.delta_r(lep_ak8_pair['hm'].i1) < 0.8).all()
            ))
            #selection.add('miniIso_hm', (
                #(leadinglep_miso['hm'] < 0.1).any()
            #))
            #selection.add('miniIsoInv_hm', (
                #(leadinglep_miso['hm'] >= 0.1).any()
            #))
    
            lep_ak8_pair['lm'] = leadinglep['lm'].cross(candidatejet)#, nested=True)
            selection[k].add('lepDrAK8_lm', (
                (lep_ak8_pair['lm'].i0.delta_r(lep_ak8_pair['lm'].i1) < 0.8).all()
            ))
            #selection.add('miniIso_lm', (
                #(leadinglep_miso['lm'] < 0.1).any()
            #))
            #selection.add('miniIsoInv_lm', (
                #(leadinglep_miso['lm'] >= 0.1).any()
            #))
    
            lep_ak8_pair['mm'] = leadinglep['mm'].cross(candidatejet)#, nested=True)
            selection[k].add('lepDrAK8_mm', (
                (lep_ak8_pair['mm'].i0.delta_r(lep_ak8_pair['mm'].i1) < 0.8).all()
            ))
            #selection.add('miniIso_mm', (
                #(leadinglep_miso['mm'] < 0.1).any()
            #))
            #selection.add('miniIsoInv_mm', (
                #(leadinglep_miso['mm'] >= 0.1).any()
            #))
    
            lep_ak8_pair['tm'] = leadinglep['tm'].cross(candidatejet)#, nested=True)
            selection[k].add('lepDrAK8_tm', (
                (lep_ak8_pair['tm'].i0.delta_r(lep_ak8_pair['tm'].i1) < 0.8).all()
            ))
            #selection.add('miniIso_tm', (
                #(leadinglep_miso['tm'] < 0.1).any()
            #))
            #selection.add('miniIsoInv_tm', (
                #(leadinglep_miso['tm'] >= 0.1).any()
            #))
    
            lep_ak8_pair['ht'] = leadinglep['ht'].cross(candidatejet)#, nested=True)
            selection[k].add('lepDrAK8_ht', (
                (lep_ak8_pair['ht'].i0.delta_r(lep_ak8_pair['ht'].i1) < 0.8).all()
            ))
            #selection.add('miniIso_ht', (
                #(leadinglep_miso['ht'] < 0.1).any()
            #))
            #selection.add('miniIsoInv_ht', (
                #(leadinglep_miso['ht'] >= 0.1).any()
            #))
    
            lep_ak8_pair['lt'] = leadinglep['lt'].cross(candidatejet)#, nested=True)
            selection[k].add('lepDrAK8_lt', (
                (lep_ak8_pair['lt'].i0.delta_r(lep_ak8_pair['lt'].i1) < 0.8).all()
            ))
            #selection.add('miniIso_lt', (
                #(leadinglep_miso['lt'] < 0.1).any()
            #))
            #selection.add('miniIsoInv_lt', (
                #(leadinglep_miso['lt'] >= 0.1).any()
            #))
    
            lep_ak8_pair['mt'] = leadinglep['mt'].cross(candidatejet)#, nested=True)
            selection[k].add('lepDrAK8_mt', (
                (lep_ak8_pair['mt'].i0.delta_r(lep_ak8_pair['mt'].i1) < 0.8).all()
            ))
            #selection.add('miniIso_mt', (
                #(leadinglep_miso['mt'] < 0.1).any()
            #))
            #selection.add('miniIsoInv_mt', (
                #(leadinglep_miso['mt'] >= 0.1).any()
            #))
    
            lep_ak8_pair['tt'] = leadinglep['tt'].cross(candidatejet)#, nested=True)
            selection[k].add('lepDrAK8_tt', (
                (lep_ak8_pair['tt'].i0.delta_r(lep_ak8_pair['tt'].i1) < 0.8).all()
            ))
            #selection.add('miniIso_tt', (
                #(leadinglep_miso['tt'] < 0.1).any()
            #))
            #selection.add('miniIsoInv_tt', (
                #(leadinglep_miso['tt'] >= 0.1).any()
            #))

        if isRealData:
            #genflavor = candidatejet.pt.zeros_like()
            w_hadhad = weights
            w_hadel = weights
            w_hadmu = weights
            genHTauTauDecay = candidatejet.pt.zeros_like()
            genHadTau1Decay = candidatejet.pt.zeros_like()
            genHadTau2Decay = candidatejet.pt.zeros_like()
            gentautaudecay = candidatejet.pt.zeros_like()
        else:
            weights.add('genweight', events.genWeight)
            add_pileup_weight(weights, events.Pileup.nPU, self._year, dataset)
            bosons = getBosons(events).pad(1, clip=True)
            genBosonPt = bosons.pt.fillna(0)
            genHTauTauDecay, genHadTau1Decay, genHadTau2Decay = getHTauTauDecayInfo(events)
            gentautaudecay = awkward.JaggedArray.fromiter([[v] for v in genHTauTauDecay])
            add_VJets_NLOkFactor(weights, genBosonPt, self._year, dataset)
            #genflavor = matchedBosonFlavor(candidatejet, bosons)
            w_hadhad = weights
            w_hadel = weights
            w_hadmu = weights
            #add_TriggerWeight(w_hadhad, candidatejet.msdcorr, candidatejet.pt, leadinglep['hl'].pt, self._year, "hadhad")
            #add_TriggerWeight(w_hadel, candidatejet.msdcorr, candidatejet.pt, leadinglep['hl'].pt, self._year, "hadel")
            #add_TriggerWeight(w_hadmu, candidatejet.msdcorr, candidatejet.pt, leadinglep['hl'].pt, self._year, "hadmu")
            #output['btagWeight'].fill(dataset=dataset, val=self._btagSF.addBtagWeight(weights, ak4_away)) #FIXME

        regions = {
            **{'hadhad_%s%s_%s'%(muid,elid,vepre): ['jetacceptance', 'hadhad_trigger', 'jetid', 'n2ddt', 'antiak4btagMediumOppHem', 'met', 'noleptons_%s%s_%s'%(muid,elid,vepre)] for vepre in ['hv','ll','mm'] for muid in ['h','l','m','t'] for elid in ['v','l','m','t']},
            **{'hadmu_%s%s_%s'%(muid,elid,vepre): ['jetacceptance', 'hadmu_trigger', 'jetid', 'antiak4btagMediumOppHem', 'met', 'onemuon_%s%s_%s'%(muid,elid,vepre), 'lepDrAK8_%s%s'%(muid,elid)] for vepre in ['hv','ll','mm'] for muid in ['h','l','m','t'] for elid in ['v','l','m','t']},
            **{'hadel_%s%s_%s'%(muid,elid,vepre): ['jetacceptance', 'hadel_trigger', 'jetid', 'antiak4btagMediumOppHem', 'met', 'oneelec_%s%s_%s'%(muid,elid,vepre), 'lepDrAK8_%s%s'%(muid,elid)] for vepre in ['hv','ll','mm'] for muid in ['h','l','m','t'] for elid in ['v','l','m','t']},
            'hadhad_base': ['jetacceptance', 'hadhad_trigger', 'jetid', 'antiak4btagMediumOppHem', 'met'],
            'hadel_base': ['jetacceptance', 'hadel_trigger', 'jetid', 'antiak4btagMediumOppHem', 'met'],
            'hadmu_base': ['jetacceptance', 'hadmu_trigger', 'jetid', 'antiak4btagMediumOppHem', 'met'],
        }
        w_dict = {
            **{'hadhad_%s%s_%s'%(muid,elid,vepre): w_hadhad for vepre in ['hv','ll','mm'] for muid in ['h','l','m','t'] for elid in ['v','l','m','t']},
            **{'hadmu_%s%s_%s'%(muid,elid,vepre): w_hadmu for vepre in ['hv','ll','mm'] for muid in ['h','l','m','t'] for elid in ['v','l','m','t']},
            **{'hadel_%s%s_%s'%(muid,elid,vepre): w_hadel for vepre in ['hv','ll','mm'] for muid in ['h','l','m','t'] for elid in ['v','l','m','t']},
            'hadhad_base': w_hadhad,
            'hadel_base': w_hadel,
            'hadmu_base': w_hadmu,
        }

        for vepre in ['hv','ll','mm']: 
            for muid in ['h','l','m','t']:
                for elid in ['v','l','m','t']:
                    output['cutflow_hadel_%s%s_%s'%(muid,elid,vepre)][dataset]['none'] += float(w_dict['hadel_%s%s_%s'%(muid,elid,vepre)].weight().sum())
                    output['cutflow_hadmu_%s%s_%s'%(muid,elid,vepre)][dataset]['none'] += float(w_dict['hadmu_%s%s_%s'%(muid,elid,vepre)].weight().sum())
                    output['cutflow_hadhad_%s%s_%s'%(muid,elid,vepre)][dataset]['none'] += float(w_dict['hadhad_%s%s_%s'%(muid,elid,vepre)].weight().sum())
                    for cut in regions['hadel_%s%s_%s'%(muid,elid,vepre)]:
                        allcuts = set()
                        allcuts.add(cut)
                        output['cutflow_hadel_%s%s_%s'%(muid,elid,vepre)][dataset][cut] += float(w_dict['hadel_%s%s_%s'%(muid,elid,vepre)].weight()[selection['hadel_%s'%vepre].all(*allcuts)].sum())
                    for cut in regions['hadmu_%s%s_%s'%(muid,elid,vepre)]:
                        allcuts = set()
                        allcuts.add(cut)
                        output['cutflow_hadmu_%s%s_%s'%(muid,elid,vepre)][dataset][cut] += float(w_dict['hadmu_%s%s_%s'%(muid,elid,vepre)].weight()[selection['hadmu_%s'%vepre].all(*allcuts)].sum())
                    for cut in regions['hadhad_%s%s_%s'%(muid,elid,vepre)]:
                        allcuts = set()
                        allcuts.add(cut)
                        output['cutflow_hadhad_%s%s_%s'%(muid,elid,vepre)][dataset][cut] += float(w_dict['hadhad_%s%s_%s'%(muid,elid,vepre)].weight()[selection['hadhad_%s'%vepre].all(*allcuts)].sum())

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
            if region.startswith('hadhad'):
                chan = 'hadhad'
            elif region.startswith('hadel'):
                chan = 'hadel'
            elif region.startswith('hadmu'):
                chan = 'hadmu'
            if region.endswith('hv'):
                prefix = 'hv'
            elif region.endswith('ll'):
                prefix = 'll'
            elif region.endswith('mm'):
                prefix = 'mm'
            else:
                prefix = 'base'
            if (prefix != 'base'): 
                lepid = region.split('_')[1]
            else:
                lepid = 'base'
            selections = regions[region]
            cut = selection['%s_%s'%(chan,prefix)].all(*selections)
            sname = 'nominal' if systematic is None else systematic
            if wmod is None:
                weight = w_dict[region].weight(modifier=systematic)[cut]
            else:
                weight = w_dict[region].weight()[cut] * wmod[cut]

            def normalize(val):
                return val[cut].pad(1, clip=True).fillna(0).flatten()


            #print(dataset)
            #print(region)
            #print("TRIG")
            #print(np.histogram(trigger_hadel[cut]))
            #print("BOSPT")
            #print(np.histogram(normalize(genBosonPt)))
            #print("JETPT")
            #print(np.histogram(normalize(candidatejet.pt)))
            #print("LEPPT")
            #print(np.histogram(normalize(leadinglep.pt)))
            #print("JLDR")
            #print(np.histogram(normalize(lep_ak8_pair.i0.delta_r(lep_ak8_pair.i1))))
            #print("LSF3")
            #print(np.histogram(normalize(candidatejet.lsf3)))
            #print("WEIGHT")
            #print(np.histogram(weight))
            #print("CUTFLOW")
            #print(output['cutflow_hadhad'][dataset])

            output['jet_kin'].fill(
                dataset=dataset,
                region=region,
                jet_pt=normalize(candidatejet.pt),
                #jet_eta=normalize(candidatejet.eta),
                jet_msd=normalize(candidatejet.msdcorr),
                genhtt=normalize(gentautaudecay),
                #weight=weight,
            )

            output['lep_kin'].fill(
                dataset=dataset,
                region=region,
                lep_pt=normalize(leadinglep[lepid].pt),
                #lep_eta=normalize(leadinglep[lepid].eta),
                #lsf3=normalize(candidatejet.lsf3),
                #lep_jet_dr=normalize(lep_ak8_pair[lepid].i0.delta_r(lep_ak8_pair[lepid].i1)),
                miso=normalize(leadinglep_miso[lepid]),
                genhtt=normalize(gentautaudecay),
                #weight=weight,
            )

            output['pt_kin'].fill(
                dataset=dataset,
                region=region,
                #met_pt=normalize(met_p4.pt),
                lep_pt=normalize(leadinglep[lepid].pt),
                jet_pt=normalize(candidatejet.pt),
                #h_pt=normalize(bosons[events.GenPart.pdgId==25].pt),
                genhtt=normalize(gentautaudecay),
                #weight=weight,
            )

            output['eta_kin'].fill(
                dataset=dataset,
                region=region,
                #met_pt=normalize(met_p4.pt),
                lep_eta=normalize(leadinglep[lepid].eta),
                jet_eta=normalize(candidatejet.eta),
                #h_pt=normalize(bosons[events.GenPart.pdgId==25].pt),
                genhtt=normalize(gentautaudecay),
                #weight=weight,
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
