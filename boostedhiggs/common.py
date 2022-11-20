import numpy as np
import awkward as ak
from coffea.nanoevents.methods import vector

#dataset_ordering = ['JetHT','SingleElectron','SingleMuon','MET','Tau']
dataset_ordering = {
  '2016':['SingleMuon','SingleElectron','MET','JetHT'],
  '2017':['SingleMuon','SingleElectron','MET','JetHT'],
  '2018':['SingleMuon','EGamma','MET','JetHT']
}

pd_to_trig = {
  '2016':{
    "Ele27_WPTight_Gsf":'SingleElectron',
    "Ele115_CaloIdVT_GsfTrkIdT":'SingleElectron',
    "Photon175":'SingleElectron',
    "Mu50":'SingleMuon',
    "TkMu50":'SingleMuon',
    "IsoMu24":'SingleMuon',
    "IsoTkMu24":'SingleMuon',
    'PFHT800':'JetHT',
    'PFHT900':'JetHT',
    'AK8PFJet360_TrimMass30':'JetHT',
    'AK8PFHT700_TrimR0p1PT0p03Mass50':'JetHT',
    'PFHT650_WideJetMJJ950DEtaJJ1p5':'JetHT',
    'PFHT650_WideJetMJJ900DEtaJJ1p5':'JetHT',
    'PFJet450':'JetHT',
    "PFMETNoMu120_PFMHTNoMu120_IDTight":"MET",
    "PFMET120_PFMHT120_IDTight":"MET",
  },
  '2017':{
    'PFHT800':"JetHT",
    'PFHT900':"JetHT",
    'AK8PFJet360_TrimMass30':"JetHT",
    'AK8PFHT700_TrimR0p1PT0p03Mass50':"JetHT",
    'PFHT650_WideJetMJJ950DEtaJJ1p5':"JetHT",
    'PFHT650_WideJetMJJ900DEtaJJ1p5':"JetHT",
    'PFJet450':"JetHT",
    'PFHT1050':"JetHT",
    'AK8PFJet400_TrimMass30':"JetHT",
    'AK8PFJet420_TrimMass30':"JetHT",
    'AK8PFHT800_TrimMass50':"JetHT",
    'PFJet500':"JetHT",
    'AK8PFJet500':"JetHT",
    'Ele50_CaloIdVT_GsfTrkIdT_PFJet165':"SingleElectron",
    'Ele115_CaloIdVT_GsfTrkIdT':"SingleElectron",
    'Photon200':"SingleElectron",
    "Ele15_IsoVVVL_PFHT600":"SingleElectron",
    "Ele35_WPTight_Gsf":"SingleElectron",
    "Ele15_IsoVVVL_PFHT450_PFMET50":"SingleElectron",
    'IsoMu27':"SingleMuon",
    'Mu50':"SingleMuon",
    "OldMu100":"SingleMuon",
    "TkMu100":"SingleMuon",
    'Mu55':"SingleMuon",
    "Mu15_IsoVVVL_PFHT600":"SingleMuon",
    "Mu15_IsoVVVL_PFHT450_PFMET50":"SingleMuon",
    "PFMETNoMu120_PFMHTNoMu120_IDTight":"MET",
    "PFMET120_PFMHT120_IDTight":"MET",
    "PFMET120_PFMHT120_IDTight_PFHT60":"MET",
  },
  '2018':{
    'PFHT800':"JetHT",
    'PFHT900':"JetHT",
    'AK8PFJet360_TrimMass30':"JetHT",
    'AK8PFHT700_TrimR0p1PT0p03Mass50':"JetHT",
    'PFHT650_WideJetMJJ950DEtaJJ1p5':"JetHT",
    'PFHT650_WideJetMJJ900DEtaJJ1p5':"JetHT",
    'PFJet450':"JetHT",
    'PFHT1050':"JetHT",
    'AK8PFJet400_TrimMass30':"JetHT",
    'AK8PFJet420_TrimMass30':"JetHT",
    'AK8PFHT800_TrimMass50':"JetHT",
    'PFJet500':"JetHT",
    'AK8PFJet500':"JetHT",
    'Ele50_CaloIdVT_GsfTrkIdT_PFJet165':"EGamma",
    'Ele115_CaloIdVT_GsfTrkIdT':"EGamma",
    "Ele15_IsoVVVL_PFHT600":"EGamma",
    "Ele35_WPTight_Gsf":"EGamma",
    "Ele32_WPTight_Gsf":"EGamma",
    "Photon200":"EGamma",
    "Ele15_IsoVVVL_PFHT450_PFMET50":"EGamma",
    'IsoMu27':"SingleMuon",
    'IsoMu24':"SingleMuon",
    'Mu50':"SingleMuon",
    'Mu55':"SingleMuon",
    "OldMu100":"SingleMuon",
    "TkMu100":"SingleMuon",
    "Mu15_IsoVVVL_PFHT600":"SingleMuon",
    "Mu15_IsoVVVL_PFHT450_PFMET50":"SingleMuon",
    "PFMETNoMu120_PFMHTNoMu120_IDTight":"MET",
    "PFMET120_PFMHT120_IDTight":"MET",
  },
}

def isOverlap(events,dataset,triggers,year):
    trig_to_pd = {}
    for p in dataset_ordering[year]:
        trig_to_pd[p] = []
    for t in triggers:
        if t not in trig_to_pd[pd_to_trig[year][t]]:
            trig_to_pd[pd_to_trig[year][t]].append(t)
    overlap = np.ones(len(events), dtype='bool')
    for p in dataset_ordering[year]:
        if dataset.startswith(p):
            pass_pd = np.zeros(len(events), dtype='bool')
            for t in trig_to_pd[p]:
                try:
                    pass_pd = pass_pd | events.HLT[t]
                except:
                    pass
            overlap = overlap & pass_pd
            break
        else:
            for t in trig_to_pd[p]:
                try:
                    overlap = overlap & np.logical_not(events.HLT[t])
                except:
                    pass
    return overlap


def getParticles(events,lo_id=22,hi_id=25,flags=['fromHardProcess', 'isLastCopy']):
    absid = np.abs(events.GenPart.pdgId)
    return events.GenPart[
        # no gluons
        (absid >= lo_id)
        & (absid <= hi_id)
        & events.GenPart.hasFlags(flags)
    ]

def getBosons(events):
    return getParticles(events)

def match(left, right, metric, maximum=np.inf):
    '''Matching utility
    For each item in ``left``, find closest item in ``right``, using function ``metric``.
    The function must accept two broadcast-compatible arrays and return a numeric array.
    If maximum is specified, mask matched elements where metric was greater than it.
    '''
    lr = ak.cartesian((left, right))
    mval = metric(lr[:,:,'0'], lr[:,:,'1'])
    idx = ak.argmin(mval, -1, keepdims=True)
    if maximum < np.inf:
        return lr['1'][mval[idx] < maximum]
    else:
        return lr['1'][idx]


def matchedBosonFlavor(candidates, bosons, maxdR=0.8):
    matched = match(candidates, bosons, lambda a, b: a.delta_r(b), maxdR)
    childid = ak.pad_none(abs(matched.children.pdgId), 1, axis=1)
    genflavor = (ak.any(childid == 5, -1)) * 3 \
              + (ak.any(childid == 4, -1)) * 2 \
              + (ak.all(childid < 4, -1)) * 1
    return ak.fill_none(genflavor, 0)

def matchedBosonFlavorLep(candidates, bosons, maxdR=0.8):
    matched = match(candidates, bosons, lambda a, b: a.delta_r(b), maxdR)
    childid = ak.pad_none(abs(matched.children.pdgId), 1, axis=1)
    genflavor = (ak.any(childid == 13, -1)) * 3 \
              + (ak.any(childid == 11, -1)) * 2 \
              + (ak.any(childid == 15, -1)) * 1 \
              + (ak.all((childid != 15) & (childid != 13) & (childid != 11), -1)) * 0
    return ak.fill_none(genflavor, 0)

def getTauTauDecayInfo(events,mod=False):
    genvistau_higgs = events.GenVisTau
    ngenvistau_higgs = ak.sum(ak.fill_none(genvistau_higgs.pt,0)>0., -1)

    el_taus = getParticles(events,11,11,['isDirectTauDecayProduct'])
    mu_taus = getParticles(events,13,13,['isDirectTauDecayProduct'])

    #could be done with ak.count() or ak.num() depending on what is desired here
    nel_taus = ak.sum(ak.fill_none(el_taus.pt, 0) > 0., 1)
    nmu_taus = ak.sum(ak.fill_none(mu_taus.pt, 0) > 0., 1)

    '''
    I'm 99% sure this doesn't do the right thing

    This SHOULD be doing some gen matching
    '''
    tau_pt = ak.pad_none(ak.concatenate((genvistau_higgs.pt, el_taus.pt, mu_taus.pt), 
            axis=1), 2, clip=True)
    tau_eta = ak.pad_none(ak.concatenate((genvistau_higgs.eta, el_taus.eta, mu_taus.eta), 
            axis=1), 2, clip=True)
    tau_phi = ak.pad_none(ak.concatenate((genvistau_higgs.phi, el_taus.phi, mu_taus.phi), 
            axis=1), 2, clip=True)
    tau_mass = ak.pad_none(ak.concatenate((genvistau_higgs.mass, el_taus.mass, mu_taus.mass), 
            axis=1), 2, clip=True)
    tau_p4 = ak.zip(
        {
            'pt' : tau_pt,
            'eta' : tau_eta,
            'phi' : tau_phi,
            'mass' : tau_mass
        },
        behavior = vector.behavior,
        with_name = 'PtEtaPhiMLorentzVector'
    )

    tau_pair = ak.cartesian((tau_p4, tau_p4))
    tau_pair_dr = tau_pair['0'].delta_r(tau_pair['1'])

    genvistau1_decay = ak.flatten(
            ak.fill_none(
                ak.pad_none(
                    genvistau_higgs[:,0:1].status, 1,clip=True), 15))
    genvistau2_decay = ak.flatten(
            ak.fill_none(
                ak.pad_none(
                    genvistau_higgs[:,1:2].status, 1,clip=True), 15))

    genTauTauDecay = np.zeros_like(ngenvistau_higgs)
    genHadTau1Decay = np.zeros_like(ngenvistau_higgs)
    genHadTau2Decay = np.zeros_like(ngenvistau_higgs)
    '''
    Meant to encode what really happened. ie tt->?? or tt->ee or tt->mumu, etc
    '''
    genTauTauDecay = np.zeros_like(ngenvistau_higgs) + 1*np.array((ngenvistau_higgs==2) & (nel_taus==0) & (nmu_taus==0)).astype(int) + 2*np.array((ngenvistau_higgs==1) & (nel_taus==1) & (nmu_taus==0)).astype(int) + 3*np.array((ngenvistau_higgs==1) & (nel_taus==0) & (nmu_taus==1)).astype(int) + 4*np.array((ngenvistau_higgs==0) & (nel_taus==1) & (nmu_taus==1)).astype(int) + 5*np.array((ngenvistau_higgs==0) & (nel_taus==2) & (nmu_taus==0)).astype(int) + 6*np.array((ngenvistau_higgs==0) & (nel_taus==0) & (nmu_taus==2)).astype(int)
    if mod:
        '''
        Flip the sign if it was because of some failed gen matching
        '''
        genTauTauDecay = genTauTauDecay * \
                ak.fill_none(ak.to_numpy(ak.any((tau_pair_dr<0.8) & (tau_pair_dr>0.), -1) 
                            & ak.all(tau_pt[:,:2]>25., -1)), 1)
    else:
        genTauTauDecay = genTauTauDecay * (2*ak.fill_none(ak.to_numpy( 
            ak.any((tau_pair_dr<0.8) & (tau_pair_dr>0.), 1) 
            & ak.all(tau_pt[:,:2]>25., 1)),0) - 1)
    genHadTau1Decay = 0 \
                    + 1*np.array((genvistau1_decay==0)).astype(int) \
                    + 2*np.array((genvistau1_decay==1) | (genvistau1_decay==2)).astype(int) \
                    + 3*np.array((genvistau1_decay==10) | (genvistau1_decay==11)).astype(int)
    genHadTau2Decay = 0 \
                    + 1*np.array((genvistau2_decay==0)).astype(int) \
                    + 2*np.array((genvistau2_decay==1)  | (genvistau2_decay==2)).astype(int) \
                    + 3*np.array((genvistau2_decay==10) | (genvistau2_decay==11)).astype(int)
    return genTauTauDecay, genHadTau1Decay, genHadTau2Decay
