import awkward as ak
import numpy as np

from coffea.nanoevents.methods import candidate, vector

import tritonclient.grpc as triton_grpc
import tritonclient.http as triton_http
import onnxruntime as rt

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
    ptsorting = ak.argsort(jet_pfcands.pt,axis=-1,ascending=False)
    jet_pfcands = jet_pfcands[ptsorting]

    #print('jet_pfcands.pt',jet_pfcands.pt)
    #print('fatjet.pt',fatjet.pt)

    feature_dict["pf_pt"] = jet_pfcands.pt / fatjet.pt
    feature_dict["pf_eta"] = jet_pfcands.eta - fatjet.eta
    feature_dict["pf_phi"] = -fatjet.delta_phi(jet_pfcands)
    feature_dict["pf_charge"] = jet_pfcands.charge
    feature_dict["pf_pdgId"] = jet_pfcands.pdgId
    feature_dict["pf_dz"] = jet_pfcands.dz
    feature_dict["pf_dzErr"] = jet_pfcands.dzErr
    feature_dict["pf_d0"] = jet_pfcands.d0
    feature_dict["pf_d0Err"] = jet_pfcands.d0Err
    feature_dict["pf_dz"] = ak.fill_none(ak.mask(feature_dict["pf_dz"],feature_dict["pf_charge"]!=0.),0.)
    feature_dict["pf_dzErr"] = ak.fill_none(ak.mask(feature_dict["pf_dzErr"],feature_dict["pf_charge"]!=0.),0.)
    feature_dict["pf_d0"] = ak.fill_none(ak.mask(feature_dict["pf_d0"],feature_dict["pf_charge"]!=0.),0.)
    feature_dict["pf_d0Err"] = ak.fill_none(ak.mask(feature_dict["pf_d0Err"],feature_dict["pf_charge"]!=0.),0.)
    feature_dict["pf_puppiWeight"] = jet_pfcands.puppiWeight
    feature_dict["pf_puppiWeightNoLep"] = jet_pfcands.puppiWeightNoLep
    feature_dict["pf_trkChi2"] = jet_pfcands.trkChi2
    feature_dict["pf_vtxChi2"] = jet_pfcands.vtxChi2

    for iid,pid in enumerate([0., 211., 13., 22., 11., 130., 1., 2., 3., 4., 5.]):
        feature_dict['pf_id%i'%iid] = np.abs(feature_dict['pf_pdgId'])==pid
    feature_dict['pf_idreg'] = feature_dict['pf_id0']*0.
    for iid in range(1,11):
        feature_dict['pf_idreg'] = feature_dict['pf_idreg'] + feature_dict['pf_id%i'%iid]*float(iid)

    feature_dict['met_covXX'] = events.MET.covXX
    feature_dict['met_covXY'] = events.MET.covXY
    feature_dict['met_covYY'] = events.MET.covYY
    #feature_dict['met_dphi'] = fatjet.delta_phi(events.MET)
    feature_dict['met_dphi'] = events.MET.delta_phi(fatjet)
    feature_dict['met_pt'] = events.MET.pt
    feature_dict['met_significance'] = events.MET.significance
    feature_dict['pupmet_pt'] = events.PuppiMET.pt
    #feature_dict['pupmet_dphi'] = fatjet.delta_phi(events.PuppiMET)
    feature_dict['pupmet_dphi'] = events.PuppiMET.delta_phi(fatjet)
    feature_dict['jet_pt'] = fatjet.pt
    feature_dict['jet_eta'] = fatjet.eta
    feature_dict['jet_phi'] = fatjet.phi
    feature_dict['jet_msd'] = fatjet.msoftdrop
    feature_dict['jet_muonenergy'] = ak.sum(ak.mask(feature_dict['pf_pt'], np.abs(feature_dict['pf_pdgId'])==13.), 1)
    feature_dict['jet_elecenergy'] = ak.sum(ak.mask(feature_dict['pf_pt'], np.abs(feature_dict['pf_pdgId'])==11.), 1)
    feature_dict['jet_photonenergy'] = ak.sum(ak.mask(feature_dict['pf_pt'], feature_dict['pf_pdgId']==22.), 1)
    feature_dict['jet_chhadronenergy'] = ak.sum(ak.mask(feature_dict['pf_pt'], np.abs(feature_dict['pf_pdgId'])==211.), 1)
    feature_dict['jet_nehadronenergy'] = ak.sum(ak.mask(feature_dict['pf_pt'], feature_dict['pf_pdgId']==130.), 1)
    feature_dict['jet_muonnum'] = ak.sum( (np.abs(feature_dict['pf_pdgId'])==13.), 1)
    feature_dict['jet_elecnum'] = ak.sum( (np.abs(feature_dict['pf_pdgId'])==11.), 1)
    feature_dict['jet_photonnum'] = ak.sum( (feature_dict['pf_pdgId']==22.), 1)
    feature_dict['jet_chhadronnum'] = ak.sum( (np.abs(feature_dict['pf_pdgId'])==211.), 1)
    feature_dict['jet_nehadronnum'] = ak.sum( (feature_dict['pf_pdgId']==130.), 1)
    feature_dict['jet_unity'] = fatjet.pt/fatjet.pt
    if len(events.MET.covXX)>1:
    # convert to numpy arrays and normalize features
        for var in feature_dict:
            a = (
                ak.pad_none(
                    feature_dict[var] if var.startswith('pf') else ak.unflatten(feature_dict[var], 1), _Nparts if var.startswith('pf') else 1, axis=1, clip=True
                )
                .to_numpy()
                .filled(fill_value=0)
            ).astype(np.float32)
            feature_dict[var] = a
    else:
        for var in feature_dict:
            if var.startswith('pf'):
                a = (
                    ak.pad_none(
                        feature_dict[var], _Nparts, axis=1, clip=True
                    )
                    .to_numpy()
                    .filled(fill_value=0)
                ).astype(np.float32)
            else:
                a = np.array(feature_dict[var]).astype(np.float32)
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
    jet_svs = events.SV[sv_ak8_pair[:,:,'1'].delta_r(sv_ak8_pair[:,:,'0'])<0.8]

    # get features

    feature_dict["sv_pt"] = jet_svs.pt / fatjet.pt
    feature_dict["sv_eta"] = jet_svs.eta - fatjet.eta
    feature_dict["sv_phi"] = jet_svs.p4.delta_phi(fatjet)
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

   #del sv_ak8_pair
    del jet_svs

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

    del elec_ak8_pair
    del jet_elecs

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

    del muon_ak8_pair
    del jet_muons

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
            'pt' : events.boostedTau.pt,
            'eta': events.boostedTau.eta,
            'phi': events.boostedTau.phi,
            'mass' : events.boostedTau.mass
        },
        behavior = vector.behavior,
        with_name='PtEtaPhiMLorentzVector'
    )

    tau_ak8_pair = ak.cartesian( (tau_p4, fatjet) )
    #jet_taus = events.Tau[tau_ak8_pair[:,:,'0'].delta_r(tau_ak8_pair[:,:,'1'])<0.8] #preUL
    jet_taus = events.boostedTau[tau_ak8_pair[:,:,'0'].delta_r(tau_ak8_pair[:,:,'1'])<0.8]

    # get features

    taupt_scale = 1.0

    feature_dict["tau_pt"] = taupt_scale * jet_taus.pt / fatjet.pt
    feature_dict["tau_eta"] = jet_taus.eta - fatjet.eta
    feature_dict["tau_phi"] = fatjet.delta_phi(jet_taus)
    feature_dict["tau_mass"] = jet_taus.mass
    feature_dict["tau_charge"] = jet_taus.charge
    feature_dict["tau_chargedIso"] = jet_taus.chargedIso
    #feature_dict["tau_dxy"] = jet_taus.dxy
    #feature_dict["tau_dz"] = jet_taus.dz
    feature_dict["tau_leadTkDeltaEta"] = jet_taus.leadTkDeltaEta
    feature_dict["tau_leadTkDeltaPhi"] = jet_taus.leadTkDeltaPhi
    feature_dict["tau_leadTkPtOverTauPt"] = jet_taus.leadTkPtOverTauPt
    feature_dict["tau_neutralIso"] = jet_taus.neutralIso
    feature_dict["tau_photonsOutsideSignalCone"] = jet_taus.photonsOutsideSignalCone
    #feature_dict["tau_rawAntiEle"] = jet_taus.rawAntiEle
    feature_dict["tau_rawAntiEle2018"] = jet_taus.rawAntiEle2018
    feature_dict["tau_rawAntiEle"] = jet_taus.rawIso
    feature_dict["tau_rawIso"] = jet_taus.rawIso
    feature_dict["tau_rawIsodR03"] = jet_taus.rawIsodR03
    #feature_dict["tau_rawMVAoldDM2017v2"] = jet_taus.rawMVAoldDM2017v2
    feature_dict["tau_rawMVAoldDM2017v2"] = jet_taus.rawIso
    #feature_dict["tau_rawMVAoldDMdR032017v2"] = jet_taus.rawMVAoldDMdR032017v2
    feature_dict["tau_rawMVAoldDMdR032017v2"] = jet_taus.rawIso

    del tau_ak8_pair
    del jet_taus

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
        'pf':['pf_pt', 'pf_eta', 'pf_phi', 'pf_charge',
            'pf_dz', 'pf_dzErr', 'pf_d0', 'pf_d0Err',
            'pf_puppiWeight', 'pf_puppiWeightNoLep', 'pf_trkChi2', 'pf_vtxChi2'],
        'pf_reg':['pf_pt', 'pf_eta', 'pf_phi', 'pf_charge',
            'pf_dz', 'pf_dzErr', 'pf_d0', 'pf_d0Err', 
            'pf_puppiWeight','pf_puppiWeightNoLep','pf_idreg', 'pf_trkChi2',
            'pf_vtxChi2'],
        'sv':['sv_dlen', 'sv_dlenSig', 'sv_dxy', 'sv_dxySig',
            'sv_chi2', 'sv_pAngle', 'sv_x', 'sv_y',
            'sv_z', 'sv_pt', 'sv_mass', 'sv_eta',
            'sv_phi'],
        'muon':['muon_charge', 'muon_dxy', 'muon_dxyErr', 'muon_dz',
            'muon_dzErr', 'muon_eta', 'muon_ip3d', 'muon_nStations',
            'muon_nTrackerLayers', 'muon_pfRelIso03_all', 'muon_pfRelIso03_chg', 'muon_phi',
            'muon_pt', 'muon_segmentComp', 'muon_sip3d', 'muon_tkRelIso'],
        'elec':['elec_charge', 'elec_convVeto', 'elec_deltaEtaSC', 'elec_dr03EcalRecHitSumEt',
            'elec_dr03HcalDepth1TowerSumEt', 'elec_dr03TkSumPt', 'elec_dxy', 'elec_dxyErr',
            'elec_dz', 'elec_dzErr', 'elec_eInvMinusPInv', 'elec_eta',
            'elec_hoe', 'elec_ip3d', 'elec_lostHits', 'elec_phi',
            'elec_pt', 'elec_r9', 'elec_sieie', 'elec_sip3d'],
        #'tau':['tau_charge', 'tau_chargedIso', 'tau_dxy', 'tau_dz',
        #     'tau_eta', 'tau_leadTkDeltaEta', 'tau_leadTkDeltaPhi', 'tau_leadTkPtOverTauPt',
        #     'tau_mass', 'tau_neutralIso', 'tau_phi', 'tau_photonsOutsideSignalCone',
        #     'tau_pt', 'tau_rawAntiEle', 'tau_rawIso', 'tau_rawIsodR03'], #preUL
        'tau':['tau_charge', 'tau_chargedIso', 'tau_eta', 'tau_leadTkDeltaEta', 
            'tau_leadTkDeltaPhi', 'tau_leadTkPtOverTauPt', 'tau_mass', 'tau_neutralIso', 
            'tau_phi', 'tau_photonsOutsideSignalCone', 'tau_pt', 'tau_rawAntiEle2018', 
            'tau_rawIso', 'tau_rawIsodR03'],
        'evt':['jet_muonenergy','jet_elecenergy','jet_photonenergy','jet_chhadronenergy',
            'jet_nehadronenergy','jet_muonnum','jet_elecnum','jet_photonnum',
            'jet_chhadronnum','jet_nehadronnum'],
        'evt_z':['jet_muonenergy','jet_elecenergy','jet_photonenergy','jet_chhadronenergy',
            'jet_nehadronenergy','jet_muonnum','jet_elecnum','jet_photonnum',
            'jet_chhadronnum','jet_nehadronnum'],
        #'evt':['jet_muonenergy','jet_elecenergy','jet_photonenergy','jet_chhadronenergy',
        #     'jet_nehadronenergy','jet_muonnum','jet_elecnum','jet_photonnum',
        #     'jet_chhadronnum','jet_nehadronnum','jet_unity','jet_unity'],
        #'evt_z':['jet_muonenergy','jet_elecenergy','jet_photonenergy','jet_chhadronenergy',
        #     'jet_nehadronenergy','jet_muonnum','jet_elecnum','jet_photonnum',
        #     'jet_chhadronnum','jet_nehadronnum','jet_unity'],
        'evt_reg':['met_covXX','met_covXY','met_covYY','met_dphi',
            'met_pt','met_significance','pupmet_pt','pupmet_dphi',
            'jet_msd','jet_pt','jet_eta','jet_phi'],
    }
    inputs_lists['pf'][-2:-2] = ['pf_id%i'%iid for iid in range(11)]
    #inputs_lists['pf_reg'][-1:] = ['pf_id%i'%iid for iid in range(11)]
       
    tagger_inputs = {
        input_name: np.concatenate(
            [
                np.expand_dims(feature_dict[key], 1)
                for key in inputs_lists[input_name]
            ],
            axis=1,
        ).transpose((0,2,1) if 'evt' not in input_name or len(sel_events.MET.covXX)>1 else (0,1))
        for input_name in inputs_lists
    }
    
    for feat in list(feature_dict):
        del feature_dict[feat]
    del feature_dict

    for input_name in tagger_inputs:
        if input_name.startswith('evt'):
            tagger_inputs[input_name] = np.reshape(tagger_inputs[input_name],(-1,len(inputs_lists[input_name])))

    inference_model_dict = {
        'IN_hadel_v6':['pf','sv','elec','tau','evt_z'],
        'IN_hadmu_v6':['pf','sv','muon','tau','evt_z'],
        'IN_hadhad_multi_v6':['pf','sv','tau','evt_z'],
        'MassReg_hadhad':['sv','pf_reg','evt_reg'],
        'MassReg_hadel':['sv','pf_reg','evt_reg'],
        'MassReg_hadmu':['sv','pf_reg','evt_reg'],
        'Ztagger_Zee_Zhe_v6':['pf','sv','elec','tau','evt_z'],
        'Ztagger_Zmm_Zhm_v6':['pf','sv','muon','tau','evt_z'],
    }
    
    particleTestData = [[[ 1.86279297e-01  , -7.81250000e-02  , -5.63964844e-02  ,  1.00000000e+00
  ,  2.13623047e-03  ,  1.70326233e-03  ,  8.27789307e-04  ,  9.55581665e-04
  ,  1.00000000e+00  ,  1.00000000e+00  ,  2.00000000e+00  ,  2.00000000e+00
  ,  0.00000000e+00]
  ,[ 6.35528564e-03  ,  1.18408203e-01  ,  1.26464844e-01  ,  1.00000000e+00
  , -1.04751587e-02  ,  1.41716003e-03  ,  1.76525116e-03  ,  1.03092194e-03
  ,  1.00000000e+00  ,  1.00000000e+00  ,  1.00000000e+00  ,  1.00000000e+00
  ,  0.00000000e+00]
  ,[ 5.60379028e-03  ,  1.14746094e-01  ,  1.20117188e-01  , -1.00000000e+00
  ,  7.11059570e-03  ,  1.58882141e-03  ,  6.83593750e-03  ,  1.50966644e-03
  ,  1.00000000e+00  ,  1.00000000e+00  ,  1.00000000e+00  ,  1.00000000e+00
  ,  0.00000000e+00]
  ,[ 2.88391113e-03  ,  1.12792969e-01  ,  1.24023438e-01  , -1.00000000e+00
  ,  1.06353760e-02  ,  1.68228149e-03  ,  1.55687332e-04  ,  1.64604187e-03
  ,  1.00000000e+00  ,  1.00000000e+00  ,  1.00000000e+00  ,  1.00000000e+00
  ,  0.00000000e+00]
  ,[ 7.59124756e-03  , -1.25488281e-01  , -1.95312500e-02  ,  0.00000000e+00
  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00
  ,  1.00000000e+00  ,  1.00000000e+00  ,  5.00000000e+00  , -1.00000000e+00
  , -1.00000000e+00]
  ,[ 2.23350525e-03  , -9.93652344e-02  ,  2.16308594e-01  ,  0.00000000e+00
  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00
  ,  1.00000000e+00  ,  1.00000000e+00  ,  5.00000000e+00  , -1.00000000e+00
  , -1.00000000e+00]
  ,[ 2.65693665e-03  ,  9.17968750e-02  ,  7.20214844e-02  ,  0.00000000e+00
  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00
  ,  1.00000000e+00  ,  1.00000000e+00  ,  3.00000000e+00  , -1.00000000e+00
  , -1.00000000e+00]
  ,[ 1.69944763e-03  ,  3.68652344e-01  , -3.37646484e-01  , -1.00000000e+00
  ,  2.18750000e+00  ,  4.71115112e-03  , -2.67982483e-03  ,  7.06863403e-03
  ,  6.98242188e-01  ,  7.13867188e-01  ,  1.00000000e+00  ,  1.00000000e+00
  ,  0.00000000e+00]
  ,[ 5.27954102e-02  , -1.32324219e-01  , -1.25000000e-01  ,  1.00000000e+00
  , -7.43103027e-03  ,  4.16183472e-03  ,  4.03976440e-03  ,  4.07028198e-03
  ,  1.00000000e+00  ,  1.00000000e+00  ,  1.00000000e+00  ,  1.00000000e+00
  ,  0.00000000e+00]
  ,[ 2.17247009e-03  ,  7.50976562e-01  ,  2.30957031e-01  ,  1.00000000e+00
  , -2.65502930e-02  ,  5.38635254e-03  , -2.25257874e-03  ,  5.45501709e-03
  ,  1.00000000e+00  ,  1.00000000e+00  ,  1.00000000e+00  ,  1.00000000e+00
  ,  0.00000000e+00]
  ,[ 1.91974640e-03  ,  5.65429688e-01  , -3.75000000e-01  ,  1.00000000e+00
  , -3.68690491e-03  ,  6.96182251e-03  , -1.60827637e-02  ,  6.92367554e-03
  ,  1.00000000e+00  ,  1.00000000e+00  ,  1.00000000e+00  ,  0.00000000e+00
  ,  0.00000000e+00]
  ,[ 1.05762482e-03  ,  2.69775391e-01  , -3.95751953e-01  , -1.00000000e+00
  ,  7.84375000e+00  ,  1.62597656e-01  , -9.46777344e-01  ,  1.25732422e-01
  ,  4.70703125e-01  ,  4.94140625e-01  ,  1.00000000e+00  ,  0.00000000e+00
  ,  0.00000000e+00]
  ,[ 1.68991089e-03  ,  2.47314453e-01  ,  4.77294922e-01  ,  0.00000000e+00
  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00
  ,  1.00000000e+00  ,  1.00000000e+00  ,  5.00000000e+00  , -1.00000000e+00
  , -1.00000000e+00]
  ,[ 3.56750488e-02  , -1.00830078e-01  , -2.66845703e-01  ,  0.00000000e+00
  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00
  ,  1.00000000e+00  ,  9.52636719e-01  ,  5.00000000e+00  , -1.00000000e+00
  , -1.00000000e+00]
  ,[ 1.27983093e-03  , -3.73046875e-01  , -7.37304688e-02  , -1.00000000e+00
  ,  7.33566284e-03  ,  5.61904907e-03  , -4.44030762e-03  ,  5.43212891e-03
  ,  1.00000000e+00  ,  1.00000000e+00  ,  1.00000000e+00  ,  2.00000000e+00
  ,  0.00000000e+00]
  ,[ 7.75337219e-04  ,  1.19140625e-01  , -7.13867188e-01  , -1.00000000e+00
  ,  2.14687500e+01  ,  5.58593750e-01  ,  8.07812500e+00  ,  5.25878906e-01
  ,  5.29296875e-01  ,  5.53222656e-01  ,  1.00000000e+00  ,  0.00000000e+00
  ,  0.00000000e+00]
  ,[ 1.24645233e-03  ,  2.17773438e-01  , -6.31835938e-01  , -1.00000000e+00
  ,  4.39147949e-02  ,  1.04522705e-02  , -6.57272339e-03  ,  1.09939575e-02
  ,  1.00000000e+00  ,  1.00000000e+00  ,  1.00000000e+00  ,  1.00000000e+00
  ,  0.00000000e+00]
  ,[ 1.15013123e-03  , -3.74023438e-01  , -5.99609375e-01  , -1.00000000e+00
  ,  1.24588013e-02  ,  9.19342041e-03  , -2.65693665e-03  ,  9.18579102e-03
  ,  1.00000000e+00  ,  1.00000000e+00  ,  1.00000000e+00  ,  0.00000000e+00
  ,  0.00000000e+00]
  ,[ 1.04522705e-03  , -4.39453125e-01  , -6.06445312e-01  ,  1.00000000e+00
  , -5.27832031e-01  ,  8.94927979e-03  , -2.72460938e-01  ,  8.10241699e-03
  ,  9.01855469e-01  ,  9.09667969e-01  ,  1.00000000e+00  ,  0.00000000e+00
  ,  0.00000000e+00]
  ,[ 2.84385681e-03  ,  1.73339844e-01  ,  1.19140625e-01  ,  0.00000000e+00
  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00
  ,  1.00000000e+00  ,  1.00000000e+00  ,  3.00000000e+00  , -1.00000000e+00
  , -1.00000000e+00]
  ,[ 5.60546875e-01  , -1.08398438e-01  , -2.80517578e-01  , -1.00000000e+00
  , -9.16290283e-03  ,  1.47323608e-02  ,  2.47192383e-02  ,  1.72576904e-02
  ,  1.00000000e+00  ,  9.37011719e-01  ,  1.00000000e+00  ,  1.00000000e+00
  ,  0.00000000e+00]
  ,[ 6.51836395e-04  ,  3.33740234e-01  , -1.16455078e-01  ,  0.00000000e+00
  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00
  ,  1.00000000e+00  ,  1.00000000e+00  ,  5.00000000e+00  , -1.00000000e+00
  , -1.00000000e+00]
  ,[ 6.51359558e-04  , -3.38623047e-01  , -6.35742188e-01  , -1.00000000e+00
  , -1.00000000e+00  , -1.00000000e+00  , -1.00000000e+00  , -1.00000000e+00
  ,  9.60937500e-01  ,  9.60937500e-01  ,  1.00000000e+00  , -1.00000000e+00
  , -1.00000000e+00]
  ,[ 2.02775002e-04  ,  6.44531250e-01  , -1.72363281e-01  ,  1.00000000e+00
  , -1.00000000e+00  , -1.00000000e+00  , -1.00000000e+00  , -1.00000000e+00
  ,  1.00000000e+00  ,  1.00000000e+00  ,  1.00000000e+00  , -1.00000000e+00
  , -1.00000000e+00]
  ,[ 0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00
  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00
  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00
  ,  0.00000000e+00]
  ,[ 0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00
  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00
  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00
  ,  0.00000000e+00]
  ,[ 0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00
  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00
  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00
  ,  0.00000000e+00]
  ,[ 0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00
  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00
  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00
  ,  0.00000000e+00]
  ,[ 0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00
  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00
  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00
  ,  0.00000000e+00]
  ,[ 0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00
  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00
  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00  ,  0.00000000e+00
  ,  0.00000000e+00]]]
    svTestData = [[[0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0.]
  ,[0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0.]
  ,[0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0.]
  ,[0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0.]
  ,[0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0. , 0.]]]
    eventTestData = [[ 7.14000000e+02  , -5.18000000e+02  ,  3.20800000e+03  ,  3.30322266e-01
   ,1.16008018e+02  ,  5.56250000e+00  ,  1.47752029e+02  ,  2.81494141e-01
   ,8.56875000e+01  ,  6.72500000e+02  ,  1.36254883e+00  ,  1.73071289e+00]]
#    #evi = np.where((events.event[presel]).to_numpy()==3743783)
#    evi = np.where((events.event[presel]).to_numpy()==2360724)
#    print('EVI: ',evi, sum(evi))
#    if sum(evi)>0:
#        evi = evi[0]
#        print('event ',((events.event[presel])[evi],', lumi ',(events.luminosityBlock[presel])[evi]))
#        print('particleTestData',{v:tagger_inputs['pf'][:,:,iv][evi] for iv,v in enumerate(inputs_lists['pf'])})
#        print('svTestData',{v:tagger_inputs['sv'][:,:,iv][evi] for iv,v in enumerate(inputs_lists['sv'])})
#        #print('elecTestData',{v:tagger_inputs['elec'][:,:,iv][evi] for iv,v in enumerate(inputs_lists['elec'])})
#        print('muonTestData',{v:tagger_inputs['muon'][:,:,iv][evi] for iv,v in enumerate(inputs_lists['muon'])})
#        print('tauTestData',{v:tagger_inputs['tau'][:,:,iv][evi] for iv,v in enumerate(inputs_lists['tau'])})
#        print('eventTestData',{v:tagger_inputs['evt_z'][:,iv][evi] for iv,v in enumerate(inputs_lists['evt_z'])})
#        print('eventRegTestData',{v:tagger_inputs['evt_reg'][:,iv][evi] for iv,v in enumerate(inputs_lists['evt_reg'])})
    #print('particleTestData',{v:np.histogram(tagger_inputs['pf'][:,:,iv]) for iv,v in enumerate(inputs_lists['pf'])})
    #print('svTestData',{v:np.histogram(tagger_inputs['sv'][:,:,iv]) for iv,v in enumerate(inputs_lists['sv'])})
    #print('elecTestData',{v:np.histogram(tagger_inputs['elec'][:,:,iv]) for iv,v in enumerate(inputs_lists['elec'])})
    #print('muonTestData',{v:np.histogram(tagger_inputs['muon'][:,:,iv]) for iv,v in enumerate(inputs_lists['muon'])})
    #print('tauTestData',{v:np.histogram(tagger_inputs['tau'][:,:,iv]) for iv,v in enumerate(inputs_lists['tau'])})
    #print('eventTestData',{v:np.histogram(tagger_inputs['evt_z'][:,iv]) for iv,v in enumerate(inputs_lists['evt_z'])})

    #print(sessions["MassReg_hadmu"].run(
    #        [sessions["MassReg_hadmu"].get_outputs()[0].name],
    #        {
    #            sessions["MassReg_hadmu"].get_inputs()[0].name:svTestData,
    #            sessions["MassReg_hadmu"].get_inputs()[1].name:particleTestData,
    #            sessions["MassReg_hadmu"].get_inputs()[2].name:eventTestData,
    #        }
    #    ))

    # run inference for selected fat jet
    tagger_outputs = {}
    for m in sessions:
        #print('Running',m,'(',[sessions[m].get_inputs()[iin].name for iin in range(len(inference_model_dict[m]))],')')
        tagger_outputs[m] = sessions[m].run(
            [sessions[m].get_outputs()[0].name],
            {sessions[m].get_inputs()[iin].name: tagger_inputs[nin] for iin,nin in enumerate(inference_model_dict[m])}
        )
        #print(m,type(tagger_outputs[m]))
    for input_name in list(tagger_inputs):
        del tagger_inputs[input_name]
    del tagger_inputs
    #for m in sessions:
    #    sessions[m]

    return tagger_outputs
