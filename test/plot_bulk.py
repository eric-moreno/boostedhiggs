import os, sys

import plot_stack
import plot_cutflow
import plot_2d
import samplelists

lumi_dict = {
  '2016':36.3,
  '2017':41.5,
  '2018':59.7,
}

doAlt = 0
doAltCutflow = True
altplots = ({ 
    'jet_pt':r"$p_{T}(j)$", 
    'jet_eta':r"$\eta(j)$", 
    'jet_msd':r"$m_{SD}(j)$", 
    'mt_lepmet':r"$m_{T}(LEP,MET)$", 
    'mt_jetmet':r"$m_{T}(j,MET)$", 
    'lep_pt':r"$p_{T}(LEP)$", 
    'lep_eta':r"$\eta(LEP)$", 
    'lep_jet_dr':r"$\Delta R(j,LEP)$", 
    'n2ddt':r"$N_{2}^{DDT}(j)$", 
    'jetlep_m':r"$m(j-LEP)$",
    'met_nopup_pt':r"$p_{T}^{miss}$",
    'met_pup_pt':r"PUPPI $p_{T}^{miss}$",
    'jetmet_dphi':r"$\Delta\phi(j,MET)$",
    'massreg':r"$m_{NN}$",
    'ztagger':r"$Z^{\ell\tau}(NN)$",
  },{
    'lep_miso':r"miniIso$(LEP)$",
    'nn_hadhad_v6':r"$D_{\phi}^{\tau_{h}\tau_{h}}$(v6)",
    'nn_hadhad_v6_qcd':r"$D_{QCD}^{\tau_{h}\tau_{h}}$(v6)",
    'nn_hadhad_v6_wjets':r"$D_{W+jets}^{\tau_{h}\tau_{h}}$(v6)",
    'nn_hadhad_v5p1':r"$D_{\phi}^{\tau_{h}\tau_{h}}$(v5.1)",
    'nn_hadhad_v5p1_qcd':r"$D_{QCD}^{\tau_{h}\tau_{h}}$(v5.1)",
    'nn_hadhad_v5p1_wjets':r"$D_{W+jets}^{\tau_{h}\tau_{h}}$(v5.1)",
    'ztagger_mu_qcd':r"$Z^{\mu\tau}_{QCD}$",
    'ztagger_mu_mm':r"$Z^{\mu\tau}_{\mu\mu}$",
    'ztagger_mu_hm':r"$Z^{\mu\tau}_{\mu\tau}$",
    'nn_hadel':r"$D_{\phi}^{e\tau_{h}}$",
    'nn_hadmu':r"$D_{\phi}^{\mu\tau_{h}}$",
    'antilep':r"Anti $LEP$ ID", 
  })[doAlt-1]
lepExceptions = ['antilep']
useData = True
year = '2017'
lumi = lumi_dict[year]
infilePre = '../condor/Dec01_NN_%s/hists_sum_'%year
tag = "NNTest_Dec01_%s"%year
massVar = "massreg"
massRange = [0.,500.]
massLabel = r"$m_{NN}$"

hadHadCut = 0.999999
hadLepCut = 0.98
lepMETCut = 50.
hadMETCut = 150.
hadCRMETCut = 100.
lowMETCut = 50.

histListBase = [
  "met_nn_kin",
]

lepHists = [
]

if doAlt > 0:
  histListBase = ['%s_kin'%var for var in altplots if 'lep' not in var or var in lepExceptions]
  lepHists = ['%s_kin'%var for var in altplots if 'lep' in var and var not in lepExceptions]

srHist = "met_nn_kin" if doAlt==0 else "massreg_kin"

regionList = [
  "hadhad_signal",
  "hadhad_signal_met",
  "hadhad_cr_anti_inv",
  "hadhad_cr_dphi_inv",
  #"hadhad_cr_b",
  "hadhad_cr_b_met",
  #"hadhad_cr_b_mu",
  "hadhad_cr_mu",
  "hadhad_cr_mu_iso",
  "hadhad_cr_b_mu_iso",
  "hadhad_cr_mu_anti_inv",
  "hadhad_cr_mu_iso_anti_inv",
  "hadhad_cr_b_mu_iso_anti_inv",
  "hadhad_cr_mu_dphi_inv",
  "hadhad_cr_mu_iso_dphi_inv",
  "hadhad_cr_b_mu_iso_dphi_inv",
  "hadmu_signal",
  "hadel_signal",
  "hadmu_cr_ztag_inv",
  "hadel_cr_ztag_inv",
  "hadmu_cr_dphi_inv",
  "hadel_cr_dphi_inv",
  "hadmu_cr_qcd",
  "hadel_cr_qcd",
  "hadmu_cr_qcd_ztag_inv",
  "hadel_cr_qcd_ztag_inv",
  "hadmu_cr_qcd_dphi_inv",
  "hadel_cr_qcd_dphi_inv",
  "hadmu_cr_b",
  "hadel_cr_b",
  "hadmu_cr_b_ztag_inv",
  "hadel_cr_b_ztag_inv",
  "hadmu_cr_b_dphi_inv",
  "hadel_cr_b_dphi_inv",
  "hadmu_cr_w",
  "hadel_cr_w",
  "hadmu_cr_w_ztag_inv",
  "hadel_cr_w_ztag_inv",
  "hadmu_cr_w_dphi_inv",
  "hadel_cr_w_dphi_inv",
]

cutList = [
  "cutflow_hadhad_signal",
  "cutflow_hadhad_signal_met",
  "cutflow_hadhad_cr_anti_inv",
  "cutflow_hadhad_cr_dphi_inv",
  #"cutflow_hadhad_cr_b",
  "cutflow_hadhad_cr_b_met",
  "cutflow_hadhad_cr_b_met_anti_inv",
  "cutflow_hadhad_cr_b_met_dphi_inv",
  #"cutflow_hadhad_cr_b_mu",
  "cutflow_hadhad_cr_b_mu_iso",
  "cutflow_hadhad_cr_mu",
  "cutflow_hadhad_cr_mu_iso",
  "cutflow_hadhad_cr_b_mu_iso_anti_inv",
  "cutflow_hadhad_cr_mu_anti_inv",
  "cutflow_hadhad_cr_mu_iso_anti_inv",
  "cutflow_hadhad_cr_b_mu_iso_dphi_inv",
  "cutflow_hadhad_cr_mu_dphi_inv",
  "cutflow_hadhad_cr_mu_iso_dphi_inv",
  "cutflow_hadel_signal",
  "cutflow_hadmu_signal",
  "cutflow_hadel_cr_ztag_inv",
  "cutflow_hadmu_cr_ztag_inv",
  "cutflow_hadel_cr_dphi_inv",
  "cutflow_hadmu_cr_dphi_inv",
  "cutflow_hadel_cr_b",
  "cutflow_hadmu_cr_b",
  "cutflow_hadmu_cr_b_ztag_inv",
  "cutflow_hadel_cr_b_ztag_inv",
  "cutflow_hadmu_cr_b_dphi_inv",
  "cutflow_hadel_cr_b_dphi_inv",
  "cutflow_hadel_cr_w",
  "cutflow_hadmu_cr_w",
  "cutflow_hadmu_cr_w_ztag_inv",
  "cutflow_hadel_cr_w_ztag_inv",
  "cutflow_hadmu_cr_w_dphi_inv",
  "cutflow_hadel_cr_w_dphi_inv",
  "cutflow_hadel_cr_qcd",
  "cutflow_hadmu_cr_qcd",
  "cutflow_hadel_cr_qcd_ztag_inv",
  "cutflow_hadmu_cr_qcd_ztag_inv",
  "cutflow_hadel_cr_qcd_dphi_inv",
  "cutflow_hadmu_cr_qcd_dphi_inv",
]

plotList = {
  "jet_nn_kin":{
    "varName":["jet_pt", "nn_disc", "nn_disc"],
    "varLabel":[r"$p_{T}(j)$", "NN", "NN"],
    "titleAdd":["", "", ""],
    "selStr":["", "", "_zoom"],
    "sels":[{}, {}, {"nn_disc":[0.9,None]}],
    #"addOpts":[{}, {"rebin":[0.,0.1,0.5,0.8,0.9,1.000001],"dosigrat":True}, {"xlimits":[0.9,1.],"dosigrat":True,"xexp":True}],
    "addOpts":[{}, {"rebin":[0.,0.1,0.5,0.8,0.9,1.000001]}, {"xlimits":[0.9,1.],"xexp":True}],
    "blindSel":["", ["None", hadLepCut], ["None", hadLepCut]],
    "addOptsBase":{"signame":["phi10","phi30","phi50","phi125","phi300"]},
  },
  "lep_nn_kin":{
    "varName":["lep_pt", "mt_lepmet", "mt_lepmet"],
    "varLabel":[r"$p_{T}(LEP)$", r"$m_{T}(LEP,MET)$", r"$m_{T}(LEP,MET)$"],
    "titleAdd":["", "", r", $NN>NNCUT$"],
    "selStr":["", "", "_nnpass"],
    "sels":[{}, {}, {"nn_disc":["NNCUT",None]}],
    "addOpts":[{}, {}, {}],
    "blindSel":["", "", ""],
    "addOptsBase":{"signame":["phi10","phi30","phi50","phi125","phi300"]},
  },
  "met_nn_kin":{
    "varName":["met_pt", "met_pt", "nn_disc", "nn_disc", "h_pt", "h_pt"],
    "varLabel":[r"PUPPI $p_{T}^{miss}$", r"PUPPI $p_{T}^{miss}$", "NN", "NN", r"$p_{T}^{NN}$", r"$p_{T}^{NN}$"],
    "titleAdd":["", r", $NN>NNCUT$", "", "", "", r", $NN>NNCUT$"],
    "selStr":["", "_nnpass", "", "_zoom", "", "_nnpass"],
    "sels":[{}, {"nn_disc":["NNCUT",None]}, {}, {"nn_disc":[0.9,None]}, {}, {"nn_disc":["NNCUT",None]}],
    #"addOpts":[{}, {}, {"rebin":[0.,0.1,0.5,0.8,0.9,1.000001],"dosigrat":True}, {"xlimits":[0.9,1.],"dosigrat":True,"xexp":True,"rebin":[0.9,0.95,0.99,0.995,0.999,0.9995,0.9999,1.000001]}, {"dosigrat":True}, {"dosigrat":True}],
    "addOpts":[{}, {}, {"rebin":[0.,0.1,0.5,0.8,0.9,1.000001]}, {"xlimits":[0.9,1.],"xexp":True,"rebin":[0.9,0.95,0.99,0.995,0.999,0.9995,0.9999,1.000001]}, {}, {}],
    "blindSel":["", "", ["None", hadLepCut], ["None", hadLepCut], "", ""],
    "addOptsBase":{"signame":["phi10","phi30","phi50","phi125","phi300"]},
  },
}

if doAlt > 0:
  plotList.update({
    '%s_kin'%var:{
      "varName":[var, ],
      "varLabel":[altplots[var], ],
      "titleAdd":["", ],
      "selStr":["", ],
      "sels":[{}, ],
      "addOpts":[{},],
      "blindSel":["", ],
      "addOptsBase":{"signame":["phi10","phi30","phi50","phi125","phi300"]},
    } for var in altplots
  })

def StringRemap(inputObj, lepType, nnCut):
  if isinstance(inputObj, dict):
    return {d:StringRemap(inputObj[d], lepType, nnCut) for d in inputObj}
  elif isinstance(inputObj, list):
    return [StringRemap(o, lepType, nnCut) for o in inputObj]
  else:
    newObj = inputObj
    if isinstance(newObj,str):
      if "LEP" in newObj:
        newObj = newObj.replace("LEP",str(lepType))
      if "NNCUT" in newObj:
        newObj = newObj.replace("NNCUT",str(nnCut))
    return newObj
  
#looseList = ["hadhad_signal", "hadhad_signal_met", "hadhad_cr_b", "hadhad_cr_b_met", "hadhad_cr_b_mu", "hadhad_cr_b_mu_iso", "hadhad_cr_mu", "hadhad_cr_mu_iso", "hadel_signal", "hadmu_signal", "hadel_cr_b", "hadel_cr_w", "hadel_cr_qcd", "hadmu_cr_b", "hadmu_cr_w", "hadmu_cr_qcd"]
looseList = regionList
lepLooseNN = "0.9"
hadLooseNN = "0.995"

def addPassFail(origList, region, ptList, keepBlind=True):
  nnCats = ["pass", "fail"]
  nnTitle = ["$NN>NNCUT$", "$NN<NNCUT$"]
  nnSel = [["NNCUT", None], [None, "NNCUT"]]
  nnBlind = [[1] if keepBlind else "", ""]
  if region in looseList:
    nnTitle[nnCats.index("fail")] = "$NN<%s$"%(hadLooseNN if "hadhad" in region else lepLooseNN)
    nnSel[nnCats.index("fail")] = [None, hadLooseNN if "hadhad" in region else lepLooseNN]

    nnCats = nnCats + ["loosepass"]
    nnTitle = nnTitle + [r"$%s<NN<NNCUT$"%(hadLooseNN if "hadhad" in region else lepLooseNN)]
    nnSel = nnSel + [["%s"%(hadLooseNN if "hadhad" in region else lepLooseNN), "NNCUT"]]
    nnBlind = nnBlind + [""]
  print(nnSel)
  if len(ptList)>1:
    passFailList = {
      "varName":[massVar for ic in range(len(nnCats)) for ip in range(len(ptList)-1)],
      "varLabel":[massLabel for ic in range(len(nnCats)) for ip in range(len(ptList)-1)],
      "titleAdd":[r", %s, $%s<p_{T}(j)<%s$"%(nnTitle[ic],str(ptList[ip]),str(ptList[ip+1])) for ic in range(len(nnCats)) for ip in range(len(ptList)-1)],
      "selStr":["_jet%sto%s_nn%s"%(str(ptList[ip]),str(ptList[ip+1]),nnCats[ic]) for ic in range(len(nnCats)) for ip in range(len(ptList)-1)],
      "sels":[{"h_pt":[ptList[ip],ptList[ip+1]], "nn_disc":nnSel[ic], massVar:massRange} for ic in range(len(nnCats)) for ip in range(len(ptList)-1)],
      "sels":[{"h_pt":[ptList[ip],ptList[ip+1]], "nn_disc":nnSel[ic]} for ic in range(len(nnCats)) for ip in range(len(ptList)-1)],
      "addOpts":[{} for ic in range(len(nnCats)) for ip in range(len(ptList)-1)],
      "blindSel":[nnBlind[ic] for ic in range(len(nnCats)) for ip in range(len(ptList)-1)],
    }
  else:
    passFailList = {
      "varName":[massVar for ic in range(len(nnCats))],
      "varLabel":[massLabel for ic in range(len(nnCats))],
      "titleAdd":[r", %s"%(nnTitle[ic]) for ic in range(len(nnCats))],
      "selStr":["_nn%s"%(nnCats[ic]) for ic in range(len(nnCats))],
      "sels":[{"nn_disc":nnSel[ic], massVar:massRange} for ic in range(len(nnCats))],
      "sels":[{"nn_disc":nnSel[ic]} for ic in range(len(nnCats))],
      "addOpts":[{} for ic in range(len(nnCats))],
      "blindSel":[nnBlind[ic] for ic in range(len(nnCats))],
    }
  newList = {obj:origList[obj]+passFailList[obj] if obj in passFailList else origList[obj] for obj in origList}
  return newList

mcSamples = samplelists.getSamplesMC(year)
dataSamples = samplelists.getSamplesData(year)

cutTitleMap = {
  "cutflow_hadhad_signal":r'\tau_{h}\tau_{h},~low~MET',
  "cutflow_hadhad_signal_met":r'\tau_{h}\tau_{h}',
  "cutflow_hadhad_cr_anti_inv":r'\tau_{h}\tau_{h}~CR~(AntiLep~inv)',
  "cutflow_hadhad_cr_dphi_inv":r'\tau_{h}\tau_{h}~CR~(\Delta\phi~inv)',
  "cutflow_hadhad_cr_b":r'\tau_{h}\tau_{h}~CR~b',
  "cutflow_hadhad_cr_b_met":r'\tau_{h}\tau_{h},~high~MET,~CR~b',
  "cutflow_hadhad_cr_b_met_anti_inv":r'\tau_{h}\tau_{h},~high~MET,~CR~b~(AntiLep~inv)',
  "cutflow_hadhad_cr_b_met_dphi_inv":r'\tau_{h}\tau_{h},~high~MET,~CR~b~(\Delta\phi~inv)',
  "cutflow_hadhad_cr_b_mu":r'\tau_{h}\tau_{h}~CR~\mu~(b)',
  "cutflow_hadhad_cr_b_mu_iso":r'\tau_{h}\tau_{h}~CR~\mu~(Iso,~b)',
  "cutflow_hadhad_cr_mu":r'\tau_{h}\tau_{h}~CR~\mu',
  "cutflow_hadhad_cr_mu_iso":r'\tau_{h}\tau_{h}~CR~\mu~(Iso)',
  "cutflow_hadhad_cr_b_mu_iso_anti_inv":r'\tau_{h}\tau_{h}~CR~\mu~(Iso,~b,~AntiLep~inv)',
  "cutflow_hadhad_cr_mu_anti_inv":r'\tau_{h}\tau_{h}~CR~\mu~(AntiLep~inv)',
  "cutflow_hadhad_cr_mu_iso_anti_inv":r'\tau_{h}\tau_{h}~CR~\mu~(Iso,~AntiLep~inv)',
  "cutflow_hadhad_cr_b_mu_iso_dphi_inv":r'\tau_{h}\tau_{h}~CR~\mu~(Iso,~b,~\Delta\phi~inv)',
  "cutflow_hadhad_cr_mu_dphi_inv":r'\tau_{h}\tau_{h}~CR~\mu~(\Delta\phi~inv)',
  "cutflow_hadhad_cr_mu_iso_dphi_inv":r'\tau_{h}\tau_{h}~CR~\mu~(Iso,~\Delta\phi~inv)',
  "cutflow_hadmu_signal":r'\mu\tau_{h}',
  "cutflow_hadel_signal":r'e\tau_{h}',
  "cutflow_hadmu_cr_ztag_inv":r'\mu\tau_{h}~CR~(Z^{\mu\tau}_{NN}~inv)',
  "cutflow_hadel_cr_ztag_inv":r'e\tau_{h}~CR~(Z^{e\tau}_{NN}~inv)',
  "cutflow_hadmu_cr_dphi_inv":r'\mu\tau_{h}~CR~(\Delta\phi~inv)',
  "cutflow_hadel_cr_dphi_inv":r'e\tau_{h}~CR~(\Delta\phi~inv)',
  "cutflow_hadmu_cr_b":r'\mu\tau_{h}~CR~b',
  "cutflow_hadel_cr_b":r'e\tau_{h}~CR~b',
  "cutflow_hadmu_cr_b_ztag_inv":r'\mu\tau_{h}~CR~b~(Z^{\mu\tau}_{NN}~inv)',
  "cutflow_hadel_cr_b_ztag_inv":r'e\tau_{h}~CR~b~(Z^{e\tau}_{NN}~inv)',
  "cutflow_hadmu_cr_b_dphi_inv":r'\mu\tau_{h}~CR~b~(\Delta\phi~inv)',
  "cutflow_hadel_cr_b_dphi_inv":r'e\tau_{h}~CR~b~(\Delta\phi~inv)',
  "cutflow_hadmu_cr_qcd":r'\mu\tau_{h}~CR~QCD',
  "cutflow_hadel_cr_qcd":r'e\tau_{h}~CR~QCD',
  "cutflow_hadmu_cr_qcd_ztag_inv":r'\mu\tau_{h}~CR~QCD~(Z^{\mu\tau}_{NN}~inv)',
  "cutflow_hadel_cr_qcd_ztag_inv":r'e\tau_{h}~CR~QCD~(Z^{e\tau}_{NN}~inv)',
  "cutflow_hadmu_cr_qcd_dphi_inv":r'\mu\tau_{h}~CR~QCD~(\Delta\phi~inv)',
  "cutflow_hadel_cr_qcd_dphi_inv":r'e\tau_{h}~CR~QCD~(\Delta\phi~inv)',
  "cutflow_hadmu_cr_w":r'\mu\tau_{h}~CR~W',
  "cutflow_hadel_cr_w":r'e\tau_{h}~CR~W',
  "cutflow_hadmu_cr_w_ztag_inv":r'\mu\tau_{h}~CR~W~(Z^{\mu\tau}_{NN}~inv)',
  "cutflow_hadel_cr_w_ztag_inv":r'e\tau_{h}~CR~W~(Z^{e\tau}_{NN}~inv)',
  "cutflow_hadmu_cr_w_dphi_inv":r'\mu\tau_{h}~CR~W~(\Delta\phi~inv)',
  "cutflow_hadel_cr_w_dphi_inv":r'e\tau_{h}~CR~W~(\Delta\phi~inv)',
}

infileList=[]
for samp in mcSamples:
  infileList.append("%s%s"%(infilePre,samp))
if useData:
  for samp in dataSamples:
    infileList.append("%s%s"%(infilePre,samp))

class CutflowArgs:
  def __init__(self, hists, tag, title, lumi, regions, defcolors):
    self.hists = hists
    self.tag = tag
    self.title = title
    self.lumi = lumi
    self.regions = regions
    self.defcolors = defcolors

cutflowArgs = CutflowArgs(infileList, tag, [cutTitleMap[c] for c in cutList], lumi, cutList, True)

if doAlt==0 or doAltCutflow:
  plot_cutflow.getPlots(cutflowArgs)

class StackArgs:
  def __init__(self, hists, tag, savetag, var, varlabel, title, regions, hist, lumi=50., sel='', sigsel='', xlimits='', xlog=False, xexp=False, blind='', solo=False, sigscale=50, sigstack=False, noratio=False, dosigrat=False, overflow='allnan', rebin=[1], norm=None, density=False, sample='all', comp=False, compsel='', complabels='', defcolors=True, verbose=False, dirname=None, qcd=[], qcdshapesel='', qcdnumsel='', qcddensel=''):
    self.hists = hists
    self.dirname = dirname
    self.tag = tag
    self.savetag = savetag
    self.var = var
    self.varlabel = varlabel
    self.title = title
    self.lumi = lumi
    self.sel = sel
    self.sigsel = sigsel
    self.regions = regions
    self.hist = hist
    self.xlimits = xlimits
    self.xlog = xlog
    self.xexp = xexp
    self.blind = blind
    self.solo = solo
    self.sigscale = sigscale
    self.sigstack = sigstack
    self.noratio = noratio
    self.dosigrat = dosigrat
    self.overflow = overflow
    self.rebin = rebin
    self.norm = norm
    self.density = density
    self.sample = sample
    self.comp = comp
    self.compsel = compsel
    self.complabels = complabels
    self.defcolors = defcolors
    self.verbose = verbose
    self.qcd = qcd
    self.qcdshapesel = qcdshapesel
    self.qcdnumsel = qcdnumsel
    self.qcddensel = qcddensel

region_opts = {
  "hadhad_signal":{
    "titleBase":r'$\tau_{h}\tau_{h}$, low MET',
    "lep":"",
    "histList":histListBase,
    "blind":True,
    "ptList":[0],#[450,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[None,hadMETCut]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadhad_signal_met":{
    "titleBase":r'$\tau_{h}\tau_{h}$',
    "lep":"",
    "histList":histListBase,
    "blind":True,
    "ptList":[0],#[300,450],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[hadMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadhad_cr_anti_inv":{
    "titleBase":r'$\tau_{h}\tau_{h}$, AntiLep Inv CR',
    "lep":"",
    "histList":histListBase,
    "blind":False,
    "ptList":[0],#[300,450],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[hadMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadhad_cr_dphi_inv":{
    "titleBase":r'$\tau_{h}\tau_{h}$, $\Delta\phi$ Inv CR',
    "lep":"",
    "histList":histListBase,
    "blind":False,
    "ptList":[0],#[300,450],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[hadMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadhad_cr_b":{
    "titleBase":r'$\tau_{h}\tau_{h}$, low MET CR b',
    "lep":"",
    "histList":histListBase,
    "blind":False,
    "ptList":[0],#[450,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[None,hadMETCut]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadhad_cr_b_met":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR b',
    "lep":"",
    "histList":histListBase,
    "blind":False,
    "ptList":[0],#[300,450],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[hadMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadhad_cr_b_met_anti_inv":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR b, AntiLep Inv CR',
    "lep":"",
    "histList":histListBase,
    "blind":False,
    "ptList":[0],#[300,450],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[hadMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadhad_cr_b_met_dphi_inv":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR b, $\Delta\phi$ Inv CR',
    "lep":"",
    "histList":histListBase,
    "blind":False,
    "ptList":[0],#[300,450],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[hadMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadhad_cr_b_mu":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (b)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[hadCRMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadhad_cr_b_mu_iso":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso, b)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[hadCRMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadhad_cr_mu":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "histList":histListBase,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[hadCRMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadhad_cr_mu_iso":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[hadCRMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadhad_cr_b_mu_iso_anti_inv":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso, b, AntiLep inv)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[hadCRMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadhad_cr_mu_anti_inv":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (AntiLep inv)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "histList":histListBase,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[hadCRMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadhad_cr_mu_iso_anti_inv":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso, AntiLep inv)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[hadCRMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadhad_cr_b_mu_iso_dphi_inv":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso, b, $\Delta\phi$ inv)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[hadCRMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadhad_cr_mu_dphi_inv":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ ($\Delta\phi$ inv)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "histList":histListBase,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[hadCRMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadhad_cr_mu_iso_dphi_inv":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso, $\Delta\phi$ inv)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[hadCRMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadmu_signal":{
    "titleBase":r'$\mu\tau_{h}$',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":True,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadel_signal":{
    "titleBase":r'$e\tau_{h}$',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":True,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadmu_cr_ztag_inv":{
    "titleBase":r'$\mu\tau_{h}$ CR ($Z^{\mu\tau}_{NN}$ Inv)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":True,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadel_cr_ztag_inv":{
    "titleBase":r'$e\tau_{h}$ CR ($Z^{e\tau}_{NN}$ Inv)',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":True,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadmu_cr_dphi_inv":{
    "titleBase":r'$\mu\tau_{h}$ CR ($\Delta\phi$ Inv)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":True,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadel_cr_dphi_inv":{
    "titleBase":r'$e\tau_{h}$ CR ($\Delta\phi$ Inv)',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":True,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadmu_cr_b":{
    "titleBase":r'$\mu\tau_{h}$ CR b',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadel_cr_b":{
    "titleBase":r'$e\tau_{h}$ CR b',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadmu_cr_b_ztag_inv":{
    "titleBase":r'$\mu\tau_{h}$ CR b ($Z^{\mu\tau}_{NN}$ Inv)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadel_cr_b_ztag_inv":{
    "titleBase":r'$e\tau_{h}$ CR b ($Z^{e\tau}_{NN}$ Inv)',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadmu_cr_b_dphi_inv":{
    "titleBase":r'$\mu\tau_{h}$ CR b ($\Delta\phi$ Inv)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadel_cr_b_dphi_inv":{
    "titleBase":r'$e\tau_{h}$ CR b ($\Delta\phi$ Inv)',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadmu_cr_qcd":{
    "titleBase":r'$\mu\tau_{h}$ CR QCD',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadel_cr_qcd":{
    "titleBase":r'$e\tau_{h}$ CR QCD',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadmu_cr_qcd_ztag_inv":{
    "titleBase":r'$\mu\tau_{h}$ CR QCD ($Z^{\mu\tau}_{NN}$ Inv)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadel_cr_qcd_ztag_inv":{
    "titleBase":r'$e\tau_{h}$ CR QCD ($Z^{e\tau}_{NN}$ Inv)',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadmu_cr_qcd_dphi_inv":{
    "titleBase":r'$\mu\tau_{h}$ CR QCD ($\Delta\phi$ Inv)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadel_cr_qcd_dphi_inv":{
    "titleBase":r'$e\tau_{h}$ CR QCD ($\Delta\phi$ Inv)',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadmu_cr_w":{
    "titleBase":r'$\mu\tau_{h}$ CR W',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadel_cr_w":{
    "titleBase":r'$e\tau_{h}$ CR W',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadmu_cr_w_ztag_inv":{
    "titleBase":r'$\mu\tau_{h}$ CR W ($Z^{\mu\tau}_{NN}$ Inv)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadel_cr_w_ztag_inv":{
    "titleBase":r'$e\tau_{h}$ CR W ($Z^{e\tau}_{NN}$ Inv)',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadmu_cr_w_dphi_inv":{
    "titleBase":r'$\mu\tau_{h}$ CR W ($\Delta\phi$ Inv)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadel_cr_w_dphi_inv":{
    "titleBase":r'$e\tau_{h}$ CR W ($\Delta\phi$ Inv)',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"",
  },
}


def doPlotting(infileList,tag,lumi,region,hist,plotList,regopts,qcd_args={},add_opts={},sample='all', keepBlind=True):
  if (hist==srHist): 
    thePlotList = StringRemap(addPassFail(plotList[hist], region, regopts["ptList"], keepBlind=keepBlind),regopts["lep"],regopts["cutNN"])
    #thePlotList = StringRemap(plotList[hist],regopts["lep"],regopts["cutNN"])
  else:
    thePlotList = StringRemap(plotList[hist],regopts["lep"],regopts["cutNN"])
  add_opts = StringRemap(add_opts,regopts["lep"],regopts["cutNN"])
  minStackArgs = StackArgs(infileList, tag, "", "", "", "", region, hist, lumi, **qcd_args)
  cofhist,qcdhists = plot_stack.getPlots(minStackArgs,True)
  pwd = os.getcwd()
  os.chdir('plots/%s/'%tag)
  print('###########################',thePlotList)
  for ivar in range(len(thePlotList["varName"])):
    print('------- Running',region,hist,thePlotList["varName"][ivar])
    if regopts["blind"]:
      blindStr = thePlotList["blindSel"][ivar]
    else:
      blindStr=""

    if thePlotList["varName"][ivar] in regopts["sigSel"]:
      print('Variable %s is in sigSel:'%thePlotList["varName"][ivar],regopts["sigSel"],' , skipping...')
      continue
    comb_opt_dict = add_opts.copy()
    comb_opt_dict.update({**thePlotList["addOpts"][ivar],**thePlotList["addOptsBase"],**{"blind":blindStr}})
    if 'comp_cut' in comb_opt_dict:
      comb_opt_dict['comp_cut'] = [[[float(cs[0]) if cs[0] is not None else None, float(cs[1]) if cs[1] is not None else None] for cs in v] for v in comb_opt_dict['comp_cut']]
    thePlotList["sels"][ivar] = {k:[float(thePlotList["sels"][ivar][k][0]) if thePlotList["sels"][ivar][k][0] is not None else None, float(thePlotList["sels"][ivar][k][1]) if thePlotList["sels"][ivar][k][1] is not None else None] for k in thePlotList["sels"][ivar]}
    ndict = regopts["specSel"].copy()
    if doAlt > 0:
      ndict = {}
    ndict.update(thePlotList["sels"][ivar])
    comb_opt_dict['qcd_cut_shape'] = ndict.copy()
    comb_opt_dict['qcd_cut_num'] = ndict.copy()
    comb_opt_dict['qcd_cut_den'] = ndict.copy()
    #comb_opt_dict['qcd_cut_shape'].pop('nn_disc',None)
    #comb_opt_dict['qcd_cut_num'].pop('nn_disc',None)
    #comb_opt_dict['qcd_cut_den'].pop('nn_disc',None)
    print(ndict)
    print(comb_opt_dict)
    print(regopts["sigSel"])
    #plot_stack.drawStack(cofhist,hist,thePlotList["varName"][ivar],thePlotList["varLabel"][ivar],regopts["titleBase"]+thePlotList["titleAdd"][ivar],sample,lumi,ndict,regopts["sigSel"],[region],region+thePlotList["selStr"][ivar]+regopts["saveLabel"],**comb_opt_dict)
    plot_stack.drawStack(cofhist,hist,thePlotList["varName"][ivar],thePlotList["varLabel"][ivar],regopts["titleBase"]+thePlotList["titleAdd"][ivar],sample,lumi,ndict,regopts["sigSel"],'',region+thePlotList["selStr"][ivar]+regopts["saveLabel"],qcd_hists=qcdhists,**comb_opt_dict)

  os.chdir(pwd)
  del cofhist

print('Nominal')
for region in regionList:
  for hist in region_opts[region]["histList"]:
    print(region, hist)
    #doPlotting(infileList,tag,lumi,region,hist,plotList,region_opts[region])
print("Done")

if doAlt > 0:
  exit()

regionList = [
  "hadhad_signal",
  #"hadhad_signal_met",
  #"hadhad_cr_anti_inv",
  #"hadhad_cr_dphi_inv",
  #"hadhad_cr_b",
  #"hadhad_cr_b_met",
  #"hadhad_cr_b_met_anti_inv",
  #"hadhad_cr_b_met_dphi_inv",
  #"hadhad_cr_b_mu",
  "hadhad_cr_mu",
  "hadhad_cr_mu_iso",
  "hadhad_cr_b_mu_iso",
  "hadhad_cr_mu_anti_inv",
  "hadhad_cr_mu_iso_anti_inv",
  "hadhad_cr_b_mu_iso_anti_inv",
  "hadhad_cr_mu_dphi_inv",
  "hadhad_cr_mu_iso_dphi_inv",
  "hadhad_cr_b_mu_iso_dphi_inv",
  "hadmu_signal",
  "hadel_signal",
  "hadmu_cr_ztag_inv",
  "hadel_cr_ztag_inv",
  "hadmu_cr_dphi_inv",
  "hadel_cr_dphi_inv",
  "hadmu_cr_qcd",
  "hadel_cr_qcd",
  "hadmu_cr_qcd_ztag_inv",
  "hadel_cr_qcd_ztag_inv",
  "hadmu_cr_qcd_dphi_inv",
  "hadel_cr_qcd_dphi_inv",
  "hadmu_cr_b",
  "hadel_cr_b",
  "hadmu_cr_w",
  "hadel_cr_w",
  "hadmu_cr_b_ztag_inv",
  "hadel_cr_b_ztag_inv",
  "hadmu_cr_b_dphi_inv",
  "hadel_cr_b_dphi_inv",
  "hadmu_cr_w_ztag_inv",
  "hadel_cr_w_ztag_inv",
  "hadmu_cr_w_dphi_inv",
  "hadel_cr_w_dphi_inv",
]

region_opts = {
  "hadhad_signal":{
    "titleBase":r'$\tau_{h}\tau_{h}$, low MET',
    "lep":"",
    "histList":histListBase,
    "blind":True,
    "ptList":[0],#[450,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[lowMETCut,hadMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadhad_signal_met":{
    "titleBase":r'$\tau_{h}\tau_{h}$, low MET',
    "lep":"",
    "histList":histListBase,
    "blind":True,
    "ptList":[0],#[300,450],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[lowMETCut,hadMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadhad_cr_anti_inv":{
    "titleBase":r'$\tau_{h}\tau_{h}$, AntiLep Inv CR, low MET',
    "lep":"",
    "histList":histListBase,
    "blind":False,
    "ptList":[0],#[300,450],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[lowMETCut,hadMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadhad_cr_dphi_inv":{
    "titleBase":r'$\tau_{h}\tau_{h}$, $\Delta\phi$ Inv CR, low MET',
    "lep":"",
    "histList":histListBase,
    "blind":False,
    "ptList":[0],#[300,450],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[lowMETCut,hadMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadhad_cr_b":{
    "titleBase":r'$\tau_{h}\tau_{h}$, low MET CR b, low MET',
    "lep":"",
    "histList":histListBase,
    "blind":False,
    "ptList":[0],#[450,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[lowMETCut,hadMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadhad_cr_b_met":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR b, low MET',
    "lep":"",
    "histList":histListBase,
    "blind":False,
    "ptList":[0],#[300,450],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[lowMETCut,hadMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadhad_cr_b_met_anti_inv":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR b, (AntiLep inv), low MET',
    "lep":"",
    "histList":histListBase,
    "blind":False,
    "ptList":[0],#[300,450],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[lowMETCut,hadMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadhad_cr_b_met_dphi_inv":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR b, ($\Delta\phi$ inv), low MET',
    "lep":"",
    "histList":histListBase,
    "blind":False,
    "ptList":[0],#[300,450],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[lowMETCut,hadMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadhad_cr_b_mu":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (b), low MET',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[lowMETCut,hadCRMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadhad_cr_b_mu_iso":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso, b), low MET',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[lowMETCut,hadCRMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadhad_cr_mu":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$, low MET',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "histList":histListBase,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[lowMETCut,hadCRMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadhad_cr_mu_iso":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso), low MET',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[lowMETCut,hadCRMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadhad_cr_b_mu_iso_anti_inv":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso, b, AntiLep inv), low MET',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[lowMETCut,hadCRMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadhad_cr_mu_anti_inv":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (AntiLep inv), low MET',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "histList":histListBase,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[lowMETCut,hadCRMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadhad_cr_mu_iso_anti_inv":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso, AntiLep inv), low MET',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[lowMETCut,hadCRMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadhad_cr_b_mu_iso_dphi_inv":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso, b, $\Delta\phi$ inv), low MET',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[lowMETCut,hadCRMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadhad_cr_mu_dphi_inv":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ ($\Delta\phi$ inv), low MET',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "histList":histListBase,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[lowMETCut,hadCRMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadhad_cr_mu_iso_dphi_inv":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso, $\Delta\phi$ inv), low MET',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[lowMETCut,hadCRMETCut]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadmu_signal":{
    "titleBase":r'$\mu\tau_{h}$, low MET',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":True,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lowMETCut,lepMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadel_signal":{
    "titleBase":r'$e\tau_{h}$, low MET',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":True,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lowMETCut,lepMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadmu_cr_ztag_inv":{
    "titleBase":r'$\mu\tau_{h}$ CR ($Z^{\mu\tau}_{NN}$ Inv), low MET',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":True,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lowMETCut,lepMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadel_cr_ztag_inv":{
    "titleBase":r'$e\tau_{h}$ CR ($Z^{e\tau}_{NN}$ Inv), low MET',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":True,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lowMETCut,lepMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadmu_cr_dphi_inv":{
    "titleBase":r'$\mu\tau_{h}$ CR ($\Delta\phi$ Inv), low MET',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":True,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lowMETCut,lepMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadel_cr_dphi_inv":{
    "titleBase":r'$e\tau_{h}$ CR ($\Delta\phi$ Inv), low MET',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":True,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lowMETCut,lepMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadmu_cr_b":{
    "titleBase":r'$\mu\tau_{h}$ CR b, low MET',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lowMETCut,lepMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadel_cr_b":{
    "titleBase":r'$e\tau_{h}$ CR b, low MET',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lowMETCut,lepMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadmu_cr_b_ztag_inv":{
    "titleBase":r'$\mu\tau_{h}$ CR b ($Z^{\mu\tau}_{NN}$ Inv), low MET',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lowMETCut,lepMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadel_cr_b_ztag_inv":{
    "titleBase":r'$e\tau_{h}$ CR b ($Z^{e\tau}_{NN}$ Inv), low MET',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lowMETCut,lepMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadmu_cr_b_dphi_inv":{
    "titleBase":r'$\mu\tau_{h}$ CR b ($\Delta\phi$ Inv), low MET',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lowMETCut,lepMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadel_cr_b_dphi_inv":{
    "titleBase":r'$e\tau_{h}$ CR b ($\Delta\phi$ Inv), low MET',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lowMETCut,lepMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadmu_cr_qcd":{
    "titleBase":r'$\mu\tau_{h}$ CR QCD, low MET',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lowMETCut,lepMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadel_cr_qcd":{
    "titleBase":r'$e\tau_{h}$ CR QCD, low MET',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lowMETCut,lepMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadmu_cr_qcd_ztag_inv":{
    "titleBase":r'$\mu\tau_{h}$ CR QCD ($Z^{\mu\tau}_{NN}$ Inv), low MET',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lowMETCut,lepMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadel_cr_qcd_ztag_inv":{
    "titleBase":r'$e\tau_{h}$ CR QCD ($Z^{e\tau}_{NN}$ Inv), low MET',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lowMETCut,lepMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadmu_cr_qcd_dphi_inv":{
    "titleBase":r'$\mu\tau_{h}$ CR QCD ($\Delta\phi$ Inv), low MET',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lowMETCut,lepMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadel_cr_qcd_dphi_inv":{
    "titleBase":r'$e\tau_{h}$ CR QCD ($\Delta\phi$ Inv), low MET',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lowMETCut,lepMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadmu_cr_w":{
    "titleBase":r'$\mu\tau_{h}$ CR W, low MET',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lowMETCut,lepMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadel_cr_w":{
    "titleBase":r'$e\tau_{h}$ CR W, low MET',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lowMETCut,lepMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadmu_cr_w_ztag_inv":{
    "titleBase":r'$\mu\tau_{h}$ CR W ($Z^{\mu\tau}_{NN}$ Inv), low MET',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lowMETCut,lepMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadel_cr_w_ztag_inv":{
    "titleBase":r'$e\tau_{h}$ CR W ($Z^{e\tau}_{NN}$ Inv), low MET',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lowMETCut,lepMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadmu_cr_w_dphi_inv":{
    "titleBase":r'$\mu\tau_{h}$ CR W ($\Delta\phi$ Inv), low MET',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lowMETCut,lepMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadel_cr_w_dphi_inv":{
    "titleBase":r'$e\tau_{h}$ CR W ($\Delta\phi$ Inv), low MET',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lowMETCut,lepMETCut]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
}

print('LowMET')
for region in regionList:
  for hist in region_opts[region]["histList"]:
    print(region, hist)
    #doPlotting(infileList,tag,lumi,region,hist,plotList,region_opts[region],keepBlind=False)
print("Done")

regionList = [
  "hadhad_signal_met",
  "hadhad_cr_mu_iso",
  "hadhad_cr_b_mu_iso",
  "hadmu_signal",
  "hadel_signal",
  "hadmu_cr_b",
  "hadel_cr_b",
  "hadmu_cr_w",
  "hadel_cr_w",
]

region_opts = {
  "hadhad_signal_met":{
    "titleBase":r'$\tau_{h}\tau_{h}$, QCD est. shape $\Delta\phi$',
    "lep":"",
    "histList":histListBase,
    "blind":True,
    "ptList":[0],#[300,450],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[hadMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstShapeDPhi",
  },
  "hadhad_cr_mu_iso":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso), QCD est. shape $\Delta\phi$',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[hadCRMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstShapeDPhi",
  },
  "hadhad_cr_b_mu_iso":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso, b), QCD est. shape $\Delta\phi$',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[hadCRMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstShapeDPhi",
  },
  "hadmu_signal":{
    "titleBase":r'$\mu\tau_{h}$, QCD est. shape $\Delta\phi$',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":True,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstShapeDPhi",
  },
  "hadel_signal":{
    "titleBase":r'$e\tau_{h}$, QCD est. shape $\Delta\phi$',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":True,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstShapeDPhi",
  },
  "hadmu_cr_b":{
    "titleBase":r'$\mu\tau_{h}$ CR b, QCD est. shape $\Delta\phi$',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstShapeDPhi",
  },
  "hadel_cr_b":{
    "titleBase":r'$e\tau_{h}$ CR b, QCD est. shape $\Delta\phi$',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstShapeDPhi",
  },
  "hadmu_cr_w":{
    "titleBase":r'$\mu\tau_{h}$ CR W, QCD est. shape $\Delta\phi$',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstShapeDPhi",
  },
  "hadel_cr_w":{
    "titleBase":r'$e\tau_{h}$ CR W, QCD est. shape $\Delta\phi$',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstShapeDPhi",
  },
}

qcdRegions = {
  "hadhad_signal_met":["hadhad_cr_mu", "hadhad_cr_dphi_inv", "hadhad_cr_mu_dphi_inv"],
  "hadhad_cr_mu_iso":["hadhad_cr_mu", "hadhad_cr_mu_iso_dphi_inv", "hadhad_cr_mu_dphi_inv"],
  "hadhad_cr_b_mu_iso":["hadhad_cr_mu", "hadhad_cr_b_mu_iso_dphi_inv", "hadhad_cr_mu_dphi_inv"],
  "hadmu_signal":["hadmu_cr_qcd","hadmu_cr_dphi_inv","hadmu_cr_qcd_dphi_inv"],
  "hadel_signal":["hadel_cr_qcd","hadel_cr_dphi_inv","hadel_cr_qcd_dphi_inv"],
  "hadmu_cr_b":["hadmu_cr_qcd","hadmu_cr_b_dphi_inv","hadmu_cr_qcd_dphi_inv"],
  "hadel_cr_b":["hadel_cr_qcd","hadel_cr_b_dphi_inv","hadel_cr_qcd_dphi_inv"],
  "hadmu_cr_w":["hadmu_cr_qcd","hadmu_cr_w_dphi_inv","hadmu_cr_qcd_dphi_inv"],
  "hadel_cr_w":["hadel_cr_qcd","hadel_cr_w_dphi_inv","hadel_cr_qcd_dphi_inv"],
}


print('Shape dPhi QCD estimation')
for region in regionList:
  for hist in region_opts[region]["histList"]:
    print(region, hist)
    #plotList[hist]["addOptsBase"]["qcd"] = 
    qcd_args = {'qcd':qcdRegions[region]}
    doPlotting(infileList,tag,lumi,region,hist,plotList,region_opts[region],keepBlind=True,qcd_args=qcd_args)
print("Done")

region_opts = {
  "hadhad_signal_met":{
    "titleBase":r'$\tau_{h}\tau_{h}$, QCD est. $\Delta\phi$ shape Iso',
    "lep":"",
    "histList":histListBase,
    "blind":True,
    "ptList":[0],#[300,450],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[hadMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstDPhiShapeIso",
  },
  "hadhad_cr_mu_iso":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso), QCD est. $\Delta\phi$ shape Iso',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[hadCRMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstDPhiShapeIso",
  },
  "hadhad_cr_b_mu_iso":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso, b), QCD est. $\Delta\phi$ shape Iso',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[hadCRMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstDPhiShapeIso",
  },
  "hadmu_signal":{
    "titleBase":r'$\mu\tau_{h}$, QCD est. $\Delta\phi$ shape Iso',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":True,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstDPhiShapeIso",
  },
  "hadel_signal":{
    "titleBase":r'$e\tau_{h}$, QCD est. $\Delta\phi$ shape Iso',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":True,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstDPhiShapeIso",
  },
  "hadmu_cr_b":{
    "titleBase":r'$\mu\tau_{h}$ CR b, QCD est. $\Delta\phi$ shape Iso',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstDPhiShapeIso",
  },
  "hadel_cr_b":{
    "titleBase":r'$e\tau_{h}$ CR b, QCD est. $\Delta\phi$ shape Iso',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstDPhiShapeIso",
  },
  "hadmu_cr_w":{
    "titleBase":r'$\mu\tau_{h}$ CR W, QCD est. $\Delta\phi$ shape Iso',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstDPhiShapeIso",
  },
  "hadel_cr_w":{
    "titleBase":r'$e\tau_{h}$ CR W, QCD est. $\Delta\phi$ shape Iso',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstDPhiShapeIso",
  },
}

qcdRegions = {
  "hadhad_signal_met":["hadhad_cr_dphi_inv", "hadhad_cr_mu", "hadhad_cr_mu_dphi_inv"],
  "hadhad_cr_mu_iso":["hadhad_cr_mu", "hadhad_cr_mu_iso_dphi_inv", "hadhad_cr_mu_dphi_inv"],
  "hadhad_cr_b_mu_iso":["hadhad_cr_b_mu_iso_dphi_inv", "hadhad_cr_mu", "hadhad_cr_mu_dphi_inv"],
  "hadmu_signal":["hadmu_cr_dphi_inv", "hadmu_cr_qcd", "hadmu_cr_qcd_dphi_inv"],
  "hadel_signal":["hadel_cr_dphi_inv", "hadel_cr_qcd", "hadel_cr_qcd_dphi_inv"],
  "hadmu_cr_b":["hadmu_cr_b_dphi_inv", "hadmu_cr_qcd", "hadmu_cr_qcd_dphi_inv"],
  "hadel_cr_b":["hadel_cr_b_dphi_inv", "hadel_cr_qcd", "hadel_cr_qcd_dphi_inv"],
  "hadmu_cr_w":["hadmu_cr_w_dphi_inv", "hadmu_cr_qcd", "hadmu_cr_qcd_dphi_inv"],
  "hadel_cr_w":["hadel_cr_w_dphi_inv", "hadel_cr_qcd", "hadel_cr_qcd_dphi_inv"],
}

print('dPhi Shape Iso QCD estimation')
for region in regionList:
  for hist in region_opts[region]["histList"]:
    print(region, hist)
    #plotList[hist]["addOptsBase"]["qcd"] = 
    qcd_args = {'qcd':qcdRegions[region]}
    doPlotting(infileList,tag,lumi,region,hist,plotList,region_opts[region],keepBlind=True,qcd_args=qcd_args)
print("Done")

region_opts = {
  "hadhad_signal_met":{
    "titleBase":r'$\tau_{h}\tau_{h}$, QCD est. shape AntiLep',
    "lep":"",
    "histList":histListBase,
    "blind":True,
    "ptList":[0],#[300,450],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[hadMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstShapeAntiLep",
  },
  "hadhad_cr_mu_iso":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso), QCD est. shape AntiLep',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[hadCRMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstShapeAntiLep",
  },
  "hadhad_cr_b_mu_iso":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso, b), QCD est. shape AntiLep',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[hadCRMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstShapeAntiLep",
  },
  "hadmu_signal":{
    "titleBase":r'$\mu\tau_{h}$, QCD est. shape ZTagger',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":True,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstShapeZTag",
  },
  "hadel_signal":{
    "titleBase":r'$e\tau_{h}$, QCD est. shape ZTagger',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":True,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstShapeZTag",
  },
  "hadmu_cr_b":{
    "titleBase":r'$\mu\tau_{h}$ CR b, QCD est. shape ZTagger',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstShapeZTag",
  },
  "hadel_cr_b":{
    "titleBase":r'$e\tau_{h}$ CR b, QCD est. shape ZTagger',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstShapeZTag",
  },
  "hadmu_cr_w":{
    "titleBase":r'$\mu\tau_{h}$ CR W, QCD est. shape ZTagger',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstShapeZTag",
  },
  "hadel_cr_w":{
    "titleBase":r'$e\tau_{h}$ CR W, QCD est. shape ZTagger',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstShapeZTag",
  },
}

qcdRegions = {
  "hadhad_signal_met":["hadhad_cr_mu", "hadhad_cr_anti_inv", "hadhad_cr_mu_anti_inv"],
  "hadhad_cr_mu_iso":["hadhad_cr_mu", "hadhad_cr_mu_iso_dphi_inv", "hadhad_cr_mu_dphi_inv"],
  "hadhad_cr_b_mu_iso":["hadhad_cr_mu", "hadhad_cr_b_mu_iso_anti_inv", "hadhad_cr_mu_anti_inv"],
  "hadmu_signal":["hadmu_cr_qcd","hadmu_cr_ztag_inv","hadmu_cr_qcd_ztag_inv"],
  "hadel_signal":["hadel_cr_qcd","hadel_cr_ztag_inv","hadel_cr_qcd_ztag_inv"],
  "hadmu_cr_b":["hadmu_cr_qcd","hadmu_cr_b_ztag_inv","hadmu_cr_qcd_ztag_inv"],
  "hadel_cr_b":["hadel_cr_qcd","hadel_cr_b_ztag_inv","hadel_cr_qcd_ztag_inv"],
  "hadmu_cr_w":["hadmu_cr_qcd","hadmu_cr_w_ztag_inv","hadmu_cr_qcd_ztag_inv"],
  "hadel_cr_w":["hadel_cr_qcd","hadel_cr_w_ztag_inv","hadel_cr_qcd_ztag_inv"],
}

print('Shape ZTagger QCD estimation')
for region in regionList:
  for hist in region_opts[region]["histList"]:
    print(region, hist)
    #plotList[hist]["addOptsBase"]["qcd"] = 
    qcd_args = {'qcd':qcdRegions[region]}
    doPlotting(infileList,tag,lumi,region,hist,plotList,region_opts[region],keepBlind=True,qcd_args=qcd_args)
print("Done")

region_opts = {
  "hadhad_signal_met":{
    "titleBase":r'$\tau_{h}\tau_{h}$, QCD est. AntiLep shape Iso',
    "lep":"",
    "histList":histListBase,
    "blind":True,
    "ptList":[0],#[300,450],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[hadMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstAntiLepShapeIso",
  },
  "hadhad_cr_mu_iso":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso), QCD est. AntiLep shape Iso',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[hadCRMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstAntiLepShapeIso",
  },
  "hadhad_cr_b_mu_iso":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso, b), QCD est. AntiLep shape Iso',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'met_pt':[hadCRMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstAntiLepShapeIso",
  },
  "hadmu_signal":{
    "titleBase":r'$\mu\tau_{h}$, QCD est. ZTagger shape Iso',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":True,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstZTagShapeIso",
  },
  "hadel_signal":{
    "titleBase":r'$e\tau_{h}$, QCD est. ZTagger shape Iso',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":True,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstZTagShapeIso",
  },
  "hadmu_cr_b":{
    "titleBase":r'$\mu\tau_{h}$ CR b, QCD est. ZTagger shape Iso',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstZTagShapeIso",
  },
  "hadel_cr_b":{
    "titleBase":r'$e\tau_{h}$ CR b, QCD est. ZTagger shape Iso',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstZTagShapeIso",
  },
  "hadmu_cr_w":{
    "titleBase":r'$\mu\tau_{h}$ CR W, QCD est. ZTagger shape Iso',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstZTagShapeIso",
  },
  "hadel_cr_w":{
    "titleBase":r'$e\tau_{h}$ CR W, QCD est. ZTagger shape Iso',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{'met_pt':[lepMETCut,None]},
    "sigSel":{},
    "saveLabel":"_qcdEstZTagShapeIso",
  },
}

qcdRegions = {
  "hadhad_signal_met":["hadhad_cr_anti_inv", "hadhad_cr_mu", "hadhad_cr_mu_anti_inv"],
  "hadhad_cr_mu_iso":["hadhad_cr_mu", "hadhad_cr_mu_iso_dphi_inv", "hadhad_cr_mu_dphi_inv"],
  "hadhad_cr_b_mu_iso":["hadhad_cr_b_mu_iso_anti_inv", "hadhad_cr_mu", "hadhad_cr_mu_anti_inv"],
  "hadmu_signal":["hadmu_cr_ztag_inv", "hadmu_cr_qcd", "hadmu_cr_qcd_ztag_inv"],
  "hadel_signal":["hadel_cr_ztag_inv", "hadel_cr_qcd", "hadel_cr_qcd_ztag_inv"],
  "hadmu_cr_b":["hadmu_cr_b_ztag_inv", "hadmu_cr_qcd", "hadmu_cr_qcd_ztag_inv"],
  "hadel_cr_b":["hadel_cr_b_ztag_inv", "hadel_cr_qcd", "hadel_cr_qcd_ztag_inv"],
  "hadmu_cr_w":["hadmu_cr_w_ztag_inv", "hadmu_cr_qcd", "hadmu_cr_qcd_ztag_inv"],
  "hadel_cr_w":["hadel_cr_w_ztag_inv", "hadel_cr_qcd", "hadel_cr_qcd_ztag_inv"],
}

print('ZTagger Shape Iso QCD estimation')
for region in regionList:
  for hist in region_opts[region]["histList"]:
    print(region, hist)
    #plotList[hist]["addOptsBase"]["qcdshapesel"] = 
    qcd_args = {'qcd':qcdRegions[region]}
    doPlotting(infileList,tag,lumi,region,hist,plotList,region_opts[region],keepBlind=True,qcd_args=qcd_args)
print("Done")

#regionList = [
#  "hadmu_signal",
#  "hadel_signal",
#  "hadmu_cr_qcd",
#  "hadel_cr_qcd",
#  "hadmu_cr_b",
#  "hadel_cr_b",
#  "hadmu_cr_w",
#  "hadel_cr_w",
#]
#
#region_opts = {
#  "hadmu_signal":{
#    "titleBase":r'$\mu\tau_{h}$ (AntiLep Comp)',
#    "lep":'\mu',
#    "histList":["jet_nn_kin"],
#    "blind":True,
#    "ptList":[0],#[300, 400, 1200],
#    "cutNN":hadLepCut,
#    "specSel":{},
#    "sigSel":{},
#    "saveLabel":"_antilepcomp",
#  },
#  "hadel_signal":{
#    "titleBase":r'$e\tau_{h}$ (AntiLep Comp)',
#    "lep":'e',
#    "histList":["jet_nn_kin"],
#    "blind":True,
#    "ptList":[0],#[300, 400, 1200],
#    "cutNN":hadLepCut,
#    "specSel":{},
#    "sigSel":{},
#    "saveLabel":"_antilepcomp",
#  },
#  "hadmu_cr_b":{
#    "titleBase":r'$\mu\tau_{h}$ CR b (AntiLep Comp)',
#    "lep":'\mu',
#    "histList":["jet_nn_kin"],
#    "blind":False,
#    "ptList":[0],#[300, 400, 1200],
#    "cutNN":hadLepCut,
#    "specSel":{},
#    "sigSel":{},
#    "saveLabel":"_antilepcomp",
#  },
#  "hadel_cr_b":{
#    "titleBase":r'$e\tau_{h}$ CR b (AntiLep Comp)',
#    "lep":'e',
#    "histList":["jet_nn_kin"],
#    "blind":False,
#    "ptList":[0],#[300, 400, 1200],
#    "cutNN":hadLepCut,
#    "specSel":{},
#    "sigSel":{},
#    "saveLabel":"_antilepcomp",
#  },
#  "hadmu_cr_qcd":{
#    "titleBase":r'$\mu\tau_{h}$ CR QCD (AntiLep Comp)',
#    "lep":'\mu',
#    "histList":["jet_nn_kin"],
#    "blind":False,
#    "ptList":[0],#[300, 400, 1200],
#    "cutNN":hadLepCut,
#    "specSel":{},
#    "sigSel":{},
#    "saveLabel":"_antilepcomp",
#  },
#  "hadel_cr_qcd":{
#    "titleBase":r'$e\tau_{h}$ CR QCD (AntiLep Comp)',
#    "lep":'e',
#    "histList":["jet_nn_kin"],
#    "blind":False,
#    "ptList":[0],#[300, 400, 1200],
#    "cutNN":hadLepCut,
#    "specSel":{},
#    "sigSel":{},
#    "saveLabel":"_antilepcomp",
#  },
#  "hadmu_cr_w":{
#    "titleBase":r'$\mu\tau_{h}$ CR W (AntiLep Comp)',
#    "lep":'\mu',
#    "histList":["jet_nn_kin"],
#    "blind":False,
#    "ptList":[0],#[300, 400, 1200],
#    "cutNN":hadLepCut,
#    "specSel":{},
#    "sigSel":{},
#    "saveLabel":"_antilepcomp",
#  },
#  "hadel_cr_w":{
#    "titleBase":r'$e\tau_{h}$ CR W (AntiLep Comp)',
#    "lep":'e',
#    "histList":["jet_nn_kin"],
#    "blind":False,
#    "ptList":[0],#[300, 400, 1200],
#    "cutNN":hadLepCut,
#    "specSel":{},
#    "sigSel":{},
#    "saveLabel":"_antilepcomp",
#  },
#}
#plotList = {
#  "jet_nn_kin":{
#    "varName":["massreg"],
#    "varLabel":[r"$m_{NN}$"],
#    "titleAdd":[""],
#    "selStr":[""],
#    "sels":[{}],
#    "addOpts":[{'noratio':True}],
#    "blindSel":[""],
#    "addOptsBase":"",
#  },
#}
#
#print('AntiLepComp')
#for region in regionList:
#  for hist in region_opts[region]["histList"]:
#    for s in ['bkg','sig','data']:
#      print(region, hist)
#      tmpopts = region_opts[region].copy()
#      tmpopts['saveLabel'] = region_opts[region]['saveLabel']+"_"+s
#      print(tmpopts)
#      newPlotList = plotList.copy()
#      print(newPlotList[hist]["titleAdd"])
#      newPlotList[hist]["titleAdd"] = ["%s, %s"%(t,s) for t in plotList[hist]["titleAdd"]]
#      print(newPlotList)
#      doPlotting(infileList,tag,lumi,region,hist,newPlotList,tmpopts,add_opts={'docomp':True,'density':True,'comp_var':['antilep','nn_disc'],'comp_cut':[[[-0.5,1.5],[0.5,1.5]],[[None,looseNN],[looseNN,"NNCUT"],["NNCUT",None]]],'complabels':[['Inclusive','Pass'],["NN<%s"%looseNN,"%s<NNCUT"%looseNN,"NN>NNCUT"]],'verbose':True},sample=s)
#print("Done")
#
