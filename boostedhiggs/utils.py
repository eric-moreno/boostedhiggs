import awkward as ak
import numpy as np

from coffea.nanoevents.methods import candidate, vector

import tritonclient.grpc as triton_grpc
import tritonclient.http as triton_http
import onnxruntime as rt
import tritongrpcclient


def getParticles(genparticles,lowid=22,highid=25,flags=['fromHardProcess', 'isLastCopy']):
    """
    returns the particle objects that satisfy a low id, 
    high id condition and have certain flags
    """
    absid = abs(genparticles.pdgId)
    return genparticles[
        ((absid >= lowid) & (absid <= highid))
        & genparticles.hasFlags(flags)
    ]

def match_HWWlepqq(genparticles,candidatefj):
    """
    return the number of matched objects (hWW*),daughters, 
    and gen flavor (enuqq, munuqq, taunuqq) 
    """
    higgs = getParticles(genparticles,25)
    is_hWW = ak.all(abs(higgs.children.pdgId)==24,axis=2)

    higgs = higgs[is_hWW]
    higgs_wstar = higgs.children[ak.argmin(higgs.children.mass,axis=2,keepdims=True)]
    higgs_w = higgs.children[ak.argmax(higgs.children.mass,axis=2,keepdims=True)]
    
    prompt_electron = getParticles(genparticles,11,11,['isPrompt','isLastCopy'])
    prompt_muon = getParticles(genparticles,13,13,['isPrompt', 'isLastCopy'])
    prompt_tau = getParticles(genparticles,15,15,['isPrompt', 'isLastCopy'])
    prompt_q = getParticles(genparticles,0,5,['fromHardProcess', 'isLastCopy'])
    prompt_q = prompt_q[abs(prompt_q.distinctParent.pdgId) == 24]
    
    dr_fj_quarks = candidatefj.delta_r(prompt_q)
    dr_fj_electrons = candidatefj.delta_r(prompt_electron)
    dr_fj_muons = candidatefj.delta_r(prompt_muon)
    dr_fj_taus = candidatefj.delta_r(prompt_tau)
    dr_daughters = ak.concatenate([dr_fj_quarks,dr_fj_electrons,dr_fj_muons,dr_fj_taus],axis=1)
    hWWlepqq_nprongs = ak.sum(dr_daughters<0.8,axis=1)
    
    n_electrons = ak.sum(prompt_electron.pt>0,axis=1)
    n_muons = ak.sum(prompt_muon.pt>0,axis=1)
    n_taus = ak.sum(prompt_tau.pt>0,axis=1)
    n_quarks = ak.sum(prompt_q.pt>0,axis=1)

    # 4(elenuqq),6(munuqq),8(taunuqq)
    hWWlepqq_flavor = (n_quarks==2)*1 + (n_electrons==1)*3 + (n_muons==1)*5 + (n_taus==1)*7
    
    matchedH = candidatefj.nearest(higgs, axis=1, threshold=0.8)
    matchedW = candidatefj.nearest(higgs_w, axis=1, threshold=0.8)
    matchedWstar = candidatefj.nearest(higgs_wstar, axis=1, threshold=0.8) 

    # 1 (H only), 4(W), 6(W star), 9(H, W and Wstar)
    hWWlepqq_matched = (
        (ak.sum(matchedH.pt > 0, axis=1)==1) * 1 
        + (ak.sum(ak.flatten(matchedW.pt > 0, axis=2), axis=1)==1) * 3 
        + (ak.sum(ak.flatten(matchedWstar.pt > 0, axis=2), axis=1)==1) * 5
    )
    
    # leptons matched
    dr_leptons = ak.concatenate([dr_fj_electrons,dr_fj_muons], axis=1)
    matched_leptons = dr_leptons < 0.8
    
    leptons = ak.concatenate([prompt_electron, prompt_muon], axis=1)
    leptons = leptons[matched_leptons]
    
    # leptons coming from W or W*
    leptons_mass = ak.firsts(leptons.distinctParent.mass)
    higgs_w_mass = ak.firsts(ak.flatten(higgs_w.mass))[ak.firsts(leptons.pt > 0)]
    higgs_wstar_mass = ak.firsts(ak.flatten(higgs_wstar.mass))[ak.firsts(leptons.pt > 0)]

    iswlepton = leptons_mass == higgs_w_mass
    iswstarlepton = leptons_mass == higgs_wstar_mass
    
    return hWWlepqq_flavor,hWWlepqq_matched,hWWlepqq_nprongs,matchedH,higgs,iswlepton,iswstarlepton

def match_Htt(genparticles,candidatefj):
    higgs = getParticles(genparticles,25)
    is_htt = ak.all(abs(higgs.children.pdgId)==15,axis=2)

    higgs = higgs[is_htt]
    
    fromtau_electron = getParticles(events.GenPart,11,11,['isDirectTauDecayProduct'])
    fromtau_muon = getParticles(events.GenPart,13,13,['isDirectTauDecayProduct'])
    tau_visible = events.GenVisTau
    
    n_visibletaus = ak.sum(tau_visible.pt>0,axis=1)
    n_electrons_fromtaus = ak.sum(fromtau_electron.pt>0,axis=1)
    n_muons_fromtaus = ak.sum(fromtau_muon.pt>0,axis=1)
    # 3(elenuqq),6(munuqq),8(taunuqq)
    htt_flavor = (n_quarks==2)*1 + (n_electrons==1)*3 + (n_muons==1)*5 + (n_taus==1)*7

    matchedH = candidatefj.nearest(higgs, axis=1, threshold=0.8)
    dr_fj_visibletaus = candidatefj.delta_r(tau_visible)
    dr_fj_electrons = candidatefj.delta_r(fromtau_electron)
    dr_fj_muons = candidatefj.delta_r(fromtau_muon)
    dr_daughters = ak.concatenate([dr_fj_visibletaus,dr_fj_electrons,dr_fj_muons],axis=1)
    # 1 (H only), 4 (H and one tau/electron or muon from tau), 5 (H and 2 taus/ele)
    htt_matched = (ak.sum(matchedH.pt>0,axis=1)==1)*1 + (ak.sum(dr_daughters<0.8,axis=1)==1)*3 + (ak.sum(dr_daughters<0.8,axis=1)==2)*5 
    
    return htt_flavor,htt_matched,matchedH,higgs

_Nparts = 30
_Nsvs = 5
_Ntaus = 3
_Nelecs = 2
_Nmuons = 2

def get_pfcands_evt_features(events, fatjet, jet_idx):
    """
    Extracts the pf_candidate features from the
    ``events`` and returns them as numpy arrays
    """
    feature_dict = {}
    evt_feature_dict = {}

    jet_pfcands = events.PFCands[
        events.FatJetPFCands.pFCandsIdx[
            (events.FatJetPFCands.pFCandsIdx != -1)
            * (events.FatJetPFCands.jetIdx == ak.flatten(jet_idx, axis=1))
        ]
    ]

    #print('jet_pfcands.pt',jet_pfcands.pt)
    #print('fatjet.pt',fatjet.pt)

    feature_dict["pf_pt"] = jet_pfcands.pt / fatjet.pt
    feature_dict["pf_eta"] = jet_pfcands.eta - fatjet.eta
    feature_dict["pf_phi"] = -fatjet.delta_phi(jet_pfcands)
    feature_dict["pf_charge"] = jet_pfcands.charge
    feature_dict["pf_pdgId"] = jet_pfcands.pdgId
    feature_dict["pf_dz"] = jet_pfcands.dz
    feature_dict["pf_d0"] = jet_pfcands.d0
    feature_dict["pf_d0Err"] = jet_pfcands.d0Err
    feature_dict["pf_puppiWeight"] = jet_pfcands.puppiWeight
    feature_dict["pf_puppiWeightNoLep"] = jet_pfcands.puppiWeightNoLep
    feature_dict["pf_trkChi2"] = jet_pfcands.trkChi2
    feature_dict["pf_vtxChi2"] = jet_pfcands.vtxChi2

    for iid,pid in enumerate([0., 211., 13., 22., 11., 130., 1., 2., 3., 4., 5.]):
        feature_dict['pf_id%i'%iid] = np.abs(feature_dict['pf_pdgId'])==pid
    feature_dict['pf_idreg'] = feature_dict['pf_id0']*0.
    for iid in range(1,11):
        feature_dict['pf_idreg'] = feature_dict['pf_idreg'] + feature_dict['pf_id0']*float(iid)

    feature_dict['met_covXX'] = ak.unflatten(events.MET.covXX, 1)
    feature_dict['met_covXY'] = ak.unflatten(events.MET.covXY, 1)
    feature_dict['met_covYY'] = ak.unflatten(events.MET.covYY, 1)
    feature_dict['met_dphi'] = ak.unflatten(fatjet.delta_phi(events.MET), 1)
    feature_dict['met_pt'] = ak.unflatten(events.MET.pt, 1)
    feature_dict['met_significance'] = ak.unflatten(events.MET.significance, 1)
    feature_dict['pupmet_pt'] = ak.unflatten(events.PuppiMET.pt, 1)
    feature_dict['pupmet_dphi'] = ak.unflatten(fatjet.delta_phi(events.PuppiMET), 1)
    feature_dict['jet_pt'] = ak.unflatten(fatjet.pt, 1)
    feature_dict['jet_eta'] = ak.unflatten(fatjet.eta, 1)
    feature_dict['jet_phi'] = ak.unflatten(fatjet.phi, 1)
    feature_dict['jet_msd'] = ak.unflatten(fatjet.msoftdrop, 1)
    feature_dict['jet_muonenergy'] = ak.unflatten(ak.sum(ak.mask(feature_dict['pf_pt'], np.abs(feature_dict['pf_pdgId'])==13.), 1), 1)
    feature_dict['jet_elecenergy'] = ak.unflatten(ak.sum(ak.mask(feature_dict['pf_pt'], np.abs(feature_dict['pf_pdgId'])==11.), 1), 1)
    feature_dict['jet_photonenergy'] = ak.unflatten(ak.sum(ak.mask(feature_dict['pf_pt'], feature_dict['pf_pdgId']==22.), 1), 1)
    feature_dict['jet_chhadronenergy'] = ak.unflatten(ak.sum(ak.mask(feature_dict['pf_pt'], np.abs(feature_dict['pf_pdgId'])==211.), 1), 1)
    feature_dict['jet_nehadronenergy'] = ak.unflatten(ak.sum(ak.mask(feature_dict['pf_pt'], feature_dict['pf_pdgId']==130.), 1), 1)
    feature_dict['jet_muonnum'] = ak.unflatten(ak.sum( (np.abs(feature_dict['pf_pdgId'])==13.), 1), 1)
    feature_dict['jet_elecnum'] = ak.unflatten(ak.sum( (np.abs(feature_dict['pf_pdgId'])==11.), 1), 1)
    feature_dict['jet_photonnum'] = ak.unflatten(ak.sum( (feature_dict['pf_pdgId']==22.), 1), 1)
    feature_dict['jet_chhadronnum'] = ak.unflatten(ak.sum( (np.abs(feature_dict['pf_pdgId'])==211.), 1), 1)
    feature_dict['jet_nehadronnum'] = ak.unflatten(ak.sum( (feature_dict['pf_pdgId']==130.), 1), 1)
    feature_dict['jet_unity'] = ak.unflatten(fatjet.pt/fatjet.pt, 1)

    # convert to numpy arrays and normalize features
    for var in feature_dict:
        a = (
            ak.pad_none(
                feature_dict[var], _Nparts if var.startswith('pf') else 1, axis=1, clip=True
            )
            .to_numpy()
            .filled(fill_value=0)
        ).astype(np.float32)

        #info = tagger_vars["sv_features"]["var_infos"][var]
        #a = (a - info["median"]) * info["norm_factor"]
        #a = np.clip(a, info.get("lower_bound", -5), info.get("upper_bound", 5))

        feature_dict[var] = a

    return feature_dict


def get_svs_features(events, fatjet, jet_idx):
    """
    Extracts the sv features specified from the
    ``events`` and returns them as numpy arrays
    """
    feature_dict = {}

    sv_p4 = ak.zip(
        {
            'pt' : events.SV.pt,
            'eta': events.SV.eta,
            'phi': events.SV.phi,
            'mass' : events.SV.mass
        },
        behavior = vector.behavior,
        with_name='PtEtaPhiMLorentzVector'
    )

    sv_ak8_pair = ak.cartesian( (sv_p4, fatjet) )
    jet_svs = events.SV[sv_ak8_pair[:,:,'0'].delta_r(sv_ak8_pair[:,:,'1'])<0.8]

    # get features

    feature_dict["sv_pt"] = jet_svs.pt / fatjet.pt
    feature_dict["sv_eta"] = jet_svs.eta - fatjet.eta
    feature_dict["sv_phi"] = jet_svs.delta_phi(fatjet)
    feature_dict["sv_mass"] = jet_svs.mass
    feature_dict["sv_dlen"] = jet_svs.dlen
    feature_dict["sv_dlenSig"] = jet_svs.dlenSig
    feature_dict["sv_dxy"] = jet_svs.dxy
    feature_dict["sv_dxySig"] = jet_svs.dxySig
    feature_dict["sv_chi2"] = jet_svs.chi2
    feature_dict["sv_pAngle"] = jet_svs.pAngle
    feature_dict["sv_x"] = jet_svs.x
    feature_dict["sv_y"] = jet_svs.y
    feature_dict["sv_z"] = jet_svs.z

    # convert to numpy arrays and normalize features
    for var in feature_dict:
        a = (
            ak.pad_none(
                feature_dict[var], _Nsvs, axis=1, clip=True
            )
            .to_numpy()
            .filled(fill_value=0)
        ).astype(np.float32)

        #info = tagger_vars["sv_features"]["var_infos"][var]
        #a = (a - info["median"]) * info["norm_factor"]
        #a = np.clip(a, info.get("lower_bound", -5), info.get("upper_bound", 5))

        feature_dict[var] = a

    return feature_dict


def get_elecs_features(events, fatjet, jet_idx):
    """
    Extracts the electrons features specified from the
    ``events`` and returns them as numpy arrays
    """
    feature_dict = {}

    elec_p4 = ak.zip(
        {
            'pt' : events.Electron.pt,
            'eta': events.Electron.eta,
            'phi': events.Electron.phi,
            'mass' : events.Electron.mass
        },
        behavior = vector.behavior,
        with_name='PtEtaPhiMLorentzVector'
    )

    elec_ak8_pair = ak.cartesian( (elec_p4, fatjet) )
    jet_elecs = events.Electron[elec_ak8_pair[:,:,'0'].delta_r(elec_ak8_pair[:,:,'1'])<0.8]

    # get features

    feature_dict["elec_pt"] = jet_elecs.pt / fatjet.pt
    feature_dict["elec_eta"] = jet_elecs.eta - fatjet.eta
    feature_dict["elec_phi"] = jet_elecs.delta_phi(fatjet)
    feature_dict["elec_mass"] = jet_elecs.mass
    feature_dict["elec_charge"] = jet_elecs.charge
    feature_dict["elec_convVeto"] = jet_elecs.convVeto
    feature_dict["elec_deltaEtaSC"] = jet_elecs.deltaEtaSC
    feature_dict["elec_dr03EcalRecHitSumEt"] = jet_elecs.dr03EcalRecHitSumEt
    feature_dict["elec_dr03HcalDepth1TowerSumEt"] = jet_elecs.dr03HcalDepth1TowerSumEt
    feature_dict["elec_dr03TkSumPt"] = jet_elecs.dr03TkSumPt
    feature_dict["elec_dxy"] = jet_elecs.dxy
    feature_dict["elec_dxyErr"] = jet_elecs.dxyErr
    feature_dict["elec_dz"] = jet_elecs.dz
    feature_dict["elec_dzErr"] = jet_elecs.dzErr
    feature_dict["elec_eInvMinusPInv"] = jet_elecs.eInvMinusPInv
    feature_dict["elec_hoe"] = jet_elecs.hoe
    feature_dict["elec_ip3d"] = jet_elecs.ip3d
    feature_dict["elec_lostHits"] = jet_elecs.lostHits
    feature_dict["elec_r9"] = jet_elecs.r9
    feature_dict["elec_sieie"] = jet_elecs.sieie
    feature_dict["elec_sip3d"] = jet_elecs.sip3d

    # convert to numpy arrays and normalize features
    for var in feature_dict:
        a = (
            ak.pad_none(
                feature_dict[var], _Nelecs, axis=1, clip=True
            )
            .to_numpy()
            .filled(fill_value=0)
        ).astype(np.float32)

        #info = tagger_vars["sv_features"]["var_infos"][var]
        #a = (a - info["median"]) * info["norm_factor"]
        #a = np.clip(a, info.get("lower_bound", -5), info.get("upper_bound", 5))

        feature_dict[var] = a

    return feature_dict


def get_muons_features(events, fatjet, jet_idx):
    """
    Extracts the muon features specified from the
    ``events`` and returns them as numpy arrays
    """
    feature_dict = {}

    muon_p4 = ak.zip(
        {
            'pt' : events.Muon.pt,
            'eta': events.Muon.eta,
            'phi': events.Muon.phi,
            'mass' : events.Muon.mass
        },
        behavior = vector.behavior,
        with_name='PtEtaPhiMLorentzVector'
    )

    muon_ak8_pair = ak.cartesian( (muon_p4, fatjet) )
    jet_muons = events.Muon[muon_ak8_pair[:,:,'0'].delta_r(muon_ak8_pair[:,:,'1'])<0.8]

    # get features

    feature_dict["muon_pt"] = jet_muons.pt / fatjet.pt
    feature_dict["muon_eta"] = jet_muons.eta - fatjet.eta
    feature_dict["muon_phi"] = jet_muons.delta_phi(fatjet)
    feature_dict["muon_mass"] = jet_muons.mass
    feature_dict["muon_charge"] = jet_muons.charge
    feature_dict["muon_dxy"] = jet_muons.dxy
    feature_dict["muon_dxyErr"] = jet_muons.dxyErr
    feature_dict["muon_dz"] = jet_muons.dz
    feature_dict["muon_dzErr"] = jet_muons.dzErr
    feature_dict["muon_ip3d"] = jet_muons.ip3d
    feature_dict["muon_nStations"] = jet_muons.nStations
    feature_dict["muon_nTrackerLayers"] = jet_muons.nTrackerLayers
    feature_dict["muon_pfRelIso03_all"] = jet_muons.pfRelIso03_all
    feature_dict["muon_pfRelIso03_chg"] = jet_muons.pfRelIso03_chg
    feature_dict["muon_segmentComp"] = jet_muons.segmentComp
    feature_dict["muon_sip3d"] = jet_muons.sip3d
    feature_dict["muon_tkRelIso"] = jet_muons.tkRelIso

    # convert to numpy arrays and normalize features
    for var in feature_dict:
        a = (
            ak.pad_none(
                feature_dict[var], _Nmuons, axis=1, clip=True
            )
            .to_numpy()
            .filled(fill_value=0)
        ).astype(np.float32)

        #info = tagger_vars["muon_features"]["var_infos"][var]
        #a = (a - info["median"]) * info["norm_factor"]
        #a = np.clip(a, info.get("lower_bound", -5), info.get("upper_bound", 5))

        feature_dict[var] = a

    return feature_dict


def get_taus_features(events, fatjet, jet_idx):
    """
    Extracts the tau features specified from the
    ``events`` and returns them as numpy arrays
    """
    feature_dict = {}

    tau_p4 = ak.zip(
        {
            'pt' : events.Tau.pt,
            'eta': events.Tau.eta,
            'phi': events.Tau.phi,
            'mass' : events.Tau.mass
        },
        behavior = vector.behavior,
        with_name='PtEtaPhiMLorentzVector'
    )

    tau_ak8_pair = ak.cartesian( (tau_p4, fatjet) )
    jet_taus = events.Tau[tau_ak8_pair[:,:,'0'].delta_r(tau_ak8_pair[:,:,'1'])<0.8]

    # get features

    feature_dict["tau_pt"] = jet_taus.pt / fatjet.pt
    feature_dict["tau_eta"] = jet_taus.eta - fatjet.eta
    feature_dict["tau_phi"] = fatjet.delta_phi(jet_taus)
    feature_dict["tau_mass"] = jet_taus.mass
    feature_dict["tau_charge"] = jet_taus.charge
    feature_dict["tau_chargedIso"] = jet_taus.chargedIso / fatjet.pt
    feature_dict["tau_dxy"] = jet_taus.dxy
    feature_dict["tau_dz"] = jet_taus.dz
    feature_dict["tau_leadTkDeltaEta"] = jet_taus.leadTkDeltaEta
    feature_dict["tau_leadTkDeltaPhi"] = jet_taus.leadTkDeltaPhi
    feature_dict["tau_leadTkPtOverTauPt"] = jet_taus.leadTkPtOverTauPt
    feature_dict["tau_neutralIso"] = jet_taus.neutralIso / fatjet.pt
    feature_dict["tau_photonsOutsideSignalCone"] = jet_taus.photonsOutsideSignalCone
    #feature_dict["tau_rawAntiEle"] = jet_taus.rawAntiEle
    feature_dict["tau_rawAntiEle"] = jet_taus.rawIso
    feature_dict["tau_rawIso"] = jet_taus.rawIso / fatjet.pt
    feature_dict["tau_rawIsodR03"] = jet_taus.rawIsodR03
    #feature_dict["tau_rawMVAoldDM2017v2"] = jet_taus.rawMVAoldDM2017v2
    feature_dict["tau_rawMVAoldDM2017v2"] = jet_taus.rawIso
    #feature_dict["tau_rawMVAoldDMdR032017v2"] = jet_taus.rawMVAoldDMdR032017v2
    feature_dict["tau_rawMVAoldDMdR032017v2"] = jet_taus.rawIso

    # convert to numpy arrays and normalize features
    for var in feature_dict:
        a = (
            ak.pad_none(
                feature_dict[var], _Ntaus, axis=1, clip=True
            )
            .to_numpy()
            .filled(fill_value=0)
        ).astype(np.float32)

        #info = tagger_vars["tau_features"]["var_infos"][var]
        #a = (a - info["median"]) * info["norm_factor"]
        #a = np.clip(a, info.get("lower_bound", -5), info.get("upper_bound", 5))

        feature_dict[var] = a

    return feature_dict

def runInferenceOnnx(events, fatjet, jet_idx, sessions, presel=None):

    # prepare inputs for both fat jets
    sel_events = events
    sel_fatjet = fatjet
    sel_jet_idx = jet_idx
    if presel is not None:
        sel_events = sel_events[presel]
        sel_fatjet = sel_fatjet[presel]
        sel_jet_idx = sel_jet_idx[presel]

    feature_dict = {
        **get_pfcands_evt_features(sel_events, sel_fatjet, sel_jet_idx),
        **get_svs_features(sel_events, sel_fatjet, sel_jet_idx),
        **get_elecs_features(sel_events, sel_fatjet, sel_jet_idx),
        **get_muons_features(sel_events, sel_fatjet, sel_jet_idx),
        **get_taus_features(sel_events, sel_fatjet, sel_jet_idx),
    }

    inputs_lists = {
        'pf':['pf_pt', 'pf_eta', 'pf_phi', 'pf_charge', 'pf_dz', 'pf_d0', 'pf_d0Err', 'pf_puppiWeight', 'pf_puppiWeightNoLep', 'pf_trkChi2', 'pf_vtxChi2'],
        'pf_reg':['pf_pt', 'pf_eta', 'pf_phi', 'pf_charge', 'pf_dz', 'pf_d0', 'pf_d0Err', 'pf_puppiWeight', 'pf_puppiWeightNoLep','pf_idreg'],
        'sv':['sv_dlen', 'sv_dlenSig', 'sv_dxy', 'sv_dxySig', 'sv_chi2', 'sv_pAngle', 'sv_x', 'sv_y', 'sv_z', 'sv_pt', 'sv_mass', 'sv_eta', 'sv_phi'],
        'muon':['muon_charge', 'muon_dxy', 'muon_dxyErr', 'muon_dz', 'muon_dzErr', 'muon_eta', 'muon_ip3d', 'muon_nStations', 'muon_nTrackerLayers', 'muon_pfRelIso03_all', 'muon_pfRelIso03_chg', 'muon_phi', 'muon_pt', 'muon_segmentComp', 'muon_sip3d', 'muon_tkRelIso'],
        'elec':['elec_charge', 'elec_convVeto', 'elec_deltaEtaSC', 'elec_dr03EcalRecHitSumEt', 'elec_dr03HcalDepth1TowerSumEt', 'elec_dr03TkSumPt', 'elec_dxy', 'elec_dxyErr', 'elec_dz', 'elec_dzErr', 'elec_eInvMinusPInv', 'elec_eta', 'elec_hoe', 'elec_ip3d', 'elec_lostHits', 'elec_phi', 'elec_pt', 'elec_r9', 'elec_sieie', 'elec_sip3d'],
        'tau':['tau_charge', 'tau_chargedIso', 'tau_dxy', 'tau_dz', 'tau_eta', 'tau_leadTkDeltaEta', 'tau_leadTkDeltaPhi', 'tau_leadTkPtOverTauPt', 'tau_mass', 'tau_neutralIso', 'tau_phi', 'tau_photonsOutsideSignalCone', 'tau_pt', 'tau_rawAntiEle', 'tau_rawIso', 'tau_rawIsodR03', 'tau_rawMVAoldDM2017v2', 'tau_rawMVAoldDMdR032017v2'],
        'evt':['jet_muonenergy','jet_elecenergy','jet_photonenergy','jet_chhadronenergy','jet_nehadronenergy','jet_muonnum','jet_elecnum','jet_photonnum','jet_chhadronnum','jet_nehadronnum','jet_unity','jet_unity'],
        'evt_z':['jet_muonenergy','jet_elecenergy','jet_photonenergy','jet_chhadronenergy','jet_nehadronenergy','jet_muonnum','jet_elecnum','jet_photonnum','jet_chhadronnum','jet_nehadronnum','jet_unity'],
        'evt_reg':['met_covXX','met_covXY','met_covYY','met_dphi','met_pt','met_significance','pupmet_pt','pupmet_dphi','jet_msd','jet_pt','jet_eta','jet_phi'],
    }
    inputs_lists['pf'][-3:-3] = ['pf_id%i'%iid for iid in range(11)]
    #inputs_lists['pf_reg'][-1:-1] = ['pf_id%i'%iid for iid in range(11)]
       
    tagger_inputs = {
        input_name: np.concatenate(
            [
                np.expand_dims(feature_dict[key], 1)
                for key in inputs_lists[input_name]
            ],
            axis=1,
        ).transpose(0,2,1)
        for input_name in inputs_lists
    }
    for input_name in tagger_inputs:
        if input_name.startswith('evt'):
            tagger_inputs[input_name] = np.reshape(tagger_inputs[input_name],(-1,len(inputs_lists[input_name])))
        print(len(events))
        print(np.sum(presel))
        print(input_name, tagger_inputs[input_name].shape)

    inference_model_dict = {
        'IN_hadel_v6':['elec','evt_z','pf','sv','tau'],
        'IN_hadmu_v6':['evt_z','muon','pf','sv','tau'],
        'model6_hadhad_multi':['evt_z','pf','sv','tau'],
        'MassReg_hadhad':['evt_reg','pf_reg','sv'],
        'MassReg_hadel':['evt_reg','pf_reg','sv'],
        'MassReg_hadmu':['evt_reg','pf_reg','sv'],
        'Ztagger_Zee_Zhe_v6':['elec','evt_z','pf','sv','tau'],
        'Ztagger_Zmm_Zhm_v6_multi':['evt_z','muon','pf','sv','tau'],
    }

    '''
    nparticles = 200
    pf (200, 30, 22)
    f_reg (200, 30, 10)
    sv (200, 5, 13)
    muon (200, 2, 16)
    elec (200, 2, 20)
    tau (200, 3, 18)
    evt (200, 12)
    evt_z (200, 11)
    evt_reg (200, 12) 

    hadel: 1
    hadmu: 1
    hadhad: 3

    mass: 2

    zee: 1
    zmm: 3

| MassReg_hadel  | 1       | UNAVAILABLE: Invalid argument: unexpected inference input 'pf_reg', allowed inputs are: inputEvent, inputParticle, inputSV |
| MassReg_hadhad | 1       | UNAVAILABLE: Invalid argument: unexpected inference input 'pf_reg', allowed inputs are: inputEvent, inputParticle, inputSV |
| MassReg_hadmu  | 1       | UNAVAILABLE: Invalid argument: unexpected inference input 'pf_reg', allowed inputs are: inputEvent, inputParticle, inputSV |
| Ztag_Zee_Zhe   | 1       | READY                                                                                                                      |
| Ztag_Zmm_Zhm   | 1       | READY                                                                                                                      |
| hadel          | 1       | UNAVAILABLE: Invalid argument: unable to load model 'hadel', configuration expects 6 inputs, model provides 5              |
| hadhad         | 1       | READY                                                                                                                      |
| hadmu          | 1       | READY                                                                                                                      |

    '''

    # run inference for both fat jets
    tagger_outputs = {
        m:sessions[m].run(
            [sessions[m].get_outputs()[0].name],
            {sessions[m].get_inputs()[iin].name: tagger_inputs[nin] for iin,nin in enumerate(inference_model_dict[m])}
        )
        for m in sessions
    }

    return tagger_outputs

URL = "0.0.0.0:8071"
verbose = False
triton_client = tritongrpcclient.InferenceServerClient(url = URL, verbose = verbose)
import multiprocessing

def runInferenceTriton(events, fatjet, jet_idx, sessions, presel=None):

    # prepare inputs for both fat jets
    sel_events = events
    sel_fatjet = fatjet
    sel_jet_idx = jet_idx
    if presel is not None:
        sel_events = sel_events[presel]
        sel_fatjet = sel_fatjet[presel]
        sel_jet_idx = sel_jet_idx[presel]

    feature_dict = {
        **get_pfcands_evt_features(sel_events, sel_fatjet, sel_jet_idx),
        **get_svs_features(sel_events, sel_fatjet, sel_jet_idx),
        **get_elecs_features(sel_events, sel_fatjet, sel_jet_idx),
        **get_muons_features(sel_events, sel_fatjet, sel_jet_idx),
        **get_taus_features(sel_events, sel_fatjet, sel_jet_idx),
    }

    inputs_lists = {
        'pf':['pf_pt', 'pf_eta', 'pf_phi', 'pf_charge', 'pf_dz', 'pf_d0', 'pf_d0Err', 'pf_puppiWeight', 'pf_puppiWeightNoLep', 'pf_trkChi2', 'pf_vtxChi2'],
        'pf_reg':['pf_pt', 'pf_eta', 'pf_phi', 'pf_charge', 'pf_dz', 'pf_d0', 'pf_d0Err', 'pf_puppiWeight', 'pf_puppiWeightNoLep','pf_idreg'],
        'sv':['sv_dlen', 'sv_dlenSig', 'sv_dxy', 'sv_dxySig', 'sv_chi2', 'sv_pAngle', 'sv_x', 'sv_y', 'sv_z', 'sv_pt', 'sv_mass', 'sv_eta', 'sv_phi'],
        'muon':['muon_charge', 'muon_dxy', 'muon_dxyErr', 'muon_dz', 'muon_dzErr', 'muon_eta', 'muon_ip3d', 'muon_nStations', 'muon_nTrackerLayers', 'muon_pfRelIso03_all', 'muon_pfRelIso03_chg', 'muon_phi', 'muon_pt', 'muon_segmentComp', 'muon_sip3d', 'muon_tkRelIso'],
        'elec':['elec_charge', 'elec_convVeto', 'elec_deltaEtaSC', 'elec_dr03EcalRecHitSumEt', 'elec_dr03HcalDepth1TowerSumEt', 'elec_dr03TkSumPt', 'elec_dxy', 'elec_dxyErr', 'elec_dz', 'elec_dzErr', 'elec_eInvMinusPInv', 'elec_eta', 'elec_hoe', 'elec_ip3d', 'elec_lostHits', 'elec_phi', 'elec_pt', 'elec_r9', 'elec_sieie', 'elec_sip3d'],
        'tau':['tau_charge', 'tau_chargedIso', 'tau_dxy', 'tau_dz', 'tau_eta', 'tau_leadTkDeltaEta', 'tau_leadTkDeltaPhi', 'tau_leadTkPtOverTauPt', 'tau_mass', 'tau_neutralIso', 'tau_phi', 'tau_photonsOutsideSignalCone', 'tau_pt', 'tau_rawAntiEle', 'tau_rawIso', 'tau_rawIsodR03', 'tau_rawMVAoldDM2017v2', 'tau_rawMVAoldDMdR032017v2'],
        #'evt':['jet_muonenergy','jet_elecenergy','jet_photonenergy','jet_chhadronenergy','jet_nehadronenergy','jet_muonnum','jet_elecnum','jet_photonnum','jet_chhadronnum','jet_nehadronnum','jet_unity','jet_unity'],
        'evt_z':['jet_muonenergy','jet_elecenergy','jet_photonenergy','jet_chhadronenergy','jet_nehadronenergy','jet_muonnum','jet_elecnum','jet_photonnum','jet_chhadronnum','jet_nehadronnum','jet_unity'],
        'evt_reg':['met_covXX','met_covXY','met_covYY','met_dphi','met_pt','met_significance','pupmet_pt','pupmet_dphi','jet_msd','jet_pt','jet_eta','jet_phi'],
    }
    inputs_lists['pf'][-3:-3] = ['pf_id%i'%iid for iid in range(11)]
    #inputs_lists['pf_reg'][-1:-1] = ['pf_id%i'%iid for iid in range(11)]
       
    tagger_inputs = {
        input_name: np.concatenate(
            [
                np.expand_dims(feature_dict[key], 1)
                for key in inputs_lists[input_name]
            ],
            axis=1,
        ).transpose(0,2,1)
        for input_name in inputs_lists
    }
    for input_name in tagger_inputs:
        if input_name.startswith('evt'):
            tagger_inputs[input_name] = np.reshape(tagger_inputs[input_name],(-1,len(inputs_lists[input_name])))
        #print(len(events))
        #print(np.sum(presel))
        #print(input_name, tagger_inputs[input_name].shape)

    inference_model_dict = {
        'hadel':['elec','evt_z','pf','sv','tau'],
        'hadmu':['evt_z','muon','pf','sv','tau'],
        'hadhad':['evt_z','pf','sv','tau'],
        'MassReg_hadhad':['evt_reg','pf_reg','sv'],
        'MassReg_hadel':['evt_reg','pf_reg','sv'],
        'MassReg_hadmu':['evt_reg','pf_reg','sv'],
        'Ztag_Zee_Zhe':['elec','evt_z','pf','sv','tau'],
        'Ztag_Zmm_Zhm':['evt_z','muon','pf','sv','tau'],
    }

    # run inference for both fat jets

    '''
    nparticles = 200
    pf (200, 30, 22)
    pf_reg (200, 30, 10)
    sv (200, 5, 13)
    muon (200, 2, 16)
    elec (200, 2, 20)
    tau (200, 3, 18)
    evt (200, 12)
    evt_z (200, 11)
    evt_reg (200, 12) 

    hadel: 1
    hadmu: 1
    hadhad: 3

    mass: 2

    zee: 1
    zmm: 3

    '''

    print("PASS","about to make tritoninputs dict", multiprocessing.current_process())
    nevt = np.sum(presel)
    tritoninputs = {
        'elec' : tritongrpcclient.InferInput("inputEl", [nevt, 2, 20], 'FP32'),
        'muon' : tritongrpcclient.InferInput("inputMu", [nevt, 2, 16], 'FP32'),
        'tau' : tritongrpcclient.InferInput('inputTau', [nevt, 3, 18], 'FP32'),
        #'evt' : tritongrpcclient.InferInput('inputEvent', [nevt, 12], 'FP32'),
        'evt_z' : tritongrpcclient.InferInput('inputEvent', [nevt, 11], 'FP32'),
        'evt_reg' : tritongrpcclient.InferInput('inputEvent', [nevt, 12], 'FP32'),
        'sv' : tritongrpcclient.InferInput('inputSV', [nevt, 5, 13], 'FP32'),
        'pf_reg' : tritongrpcclient.InferInput('inputParticle', [nevt, 30, 10], 'FP32'),
        'pf' : tritongrpcclient.InferInput('inputParticle', [nevt, 30, 22], 'FP32'),
    }

    print("PASS", 'filling tritoninputs', multiprocessing.current_process())
    size = 0
    for key in tagger_inputs.keys():
        tritoninputs[key].set_data_from_numpy(tagger_inputs[key])
        size += len(tritoninputs[key]._get_content())

    print("PASS", 'setting up output', multiprocessing.current_process())
    outputs = []
    #for key in inference_model_dict.keys():
        #outputs.append(tritongrpcclient.InferRequestedOutput(key))
    outputs.append(tritongrpcclient.InferRequestedOutput("output"))

    print("PASS", 'getting outputs', multiprocessing.current_process())
    resultdict = {}
    for key in inference_model_dict.keys():
        print("PASS", "\t",key, multiprocessing.current_process())
        results = triton_client.infer(
                model_name = key,
                inputs = [tritoninputs[name] for name in inference_model_dict[key]],
                outputs = outputs)
        result_data = results.as_numpy('output')
        resultdict[key] = result_data
    '''
    results = triton_client.infer(
            model_name = 'ensemble',
            inputs = tritoninputs.values(),
            outputs = outputs)
    for key in inference_model_dict.keys():
        resultdict[key] = results.as_numpy(key)
        print(key, resultdict[key].shape)
    '''

    print("PASS", "end of triton function", multiprocessing.current_process())
    return resultdict
