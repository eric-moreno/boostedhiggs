import os, sys

import plot_stack
import plot_cutflow
import plot_2d
import samplelists

infilePre = '../condor/Jun10_NN/hists_sum_'
useData = True
#lumi = 36.7
lumi = 41.5
year = '2017'
#lumi = 59.7
#tag = "NNTest_Jun30_2018"
tag = "NNTest_Jun10_HWW"
massVar = "massreg"
massRange = [50.,210.]
massLabel = r"$m_{NN}$"

hadHadCut = 0.99999
hadLepCut = 0.99

histListBase = [
  #"jet_nn_kin", 
  "met_nn_kin",
]

lepHists = [
]

srHist = "met_nn_kin"

regionList = [
  "hadhad_signal",
  "hadhad_signal_met",
  #"hadhad_cr_b",
  #"hadhad_cr_b_met",
  #"hadhad_cr_mu",
  #"hadhad_cr_mu_iso",
  #"hadhad_cr_b_mu",
  #"hadhad_cr_b_mu_iso",
  "hadmu_signal",
  "hadel_signal",
  "hadmu_cr_qcd",
  "hadel_cr_qcd",
  "hadmu_cr_b",
  "hadel_cr_b",
  "hadmu_cr_w",
  "hadel_cr_w",
]

cutList = [
  "cutflow_hadhad",
  "cutflow_hadhad_met",
  "cutflow_hadhad_cr_b",
  "cutflow_hadhad_cr_b_met",
  "cutflow_hadhad_cr_b_mu",
  "cutflow_hadhad_cr_b_mu_iso",
  "cutflow_hadhad_cr_mu",
  "cutflow_hadhad_cr_mu_iso",
  "cutflow_hadel",
  "cutflow_hadmu",
  "cutflow_hadel_cr_b",
  "cutflow_hadmu_cr_b",
  "cutflow_hadel_cr_w",
  "cutflow_hadmu_cr_w",
  "cutflow_hadel_cr_qcd",
  "cutflow_hadmu_cr_qcd",
]

plotList = {
  "jet_nn_kin":{
    "varName":["jet_pt", "nn_disc", "nn_disc"],
    "varLabel":[r"$p_{T}(j)$", "NN", "NN"],
    "titleAdd":["", "", ""],
    "selStr":["", "", "_zoom"],
    "sels":[{}, {}, {"nn_disc":[0.9,None]}],
    "addOpts":[{}, {"rebin":[0.,0.1,0.5,0.8,0.9,1.00001],"dosigrat":True}, {"xlimits":[0.9,1.],"dosigrat":True,"xexp":True}],
    "blindSel":["", ["None", hadLepCut], ["None", hadLepCut]],
    "addOptsBase":"",
  },
  "lep_nn_kin":{
    "varName":["lep_pt", "mt_lepmet", "mt_lepmet"],
    "varLabel":[r"$p_{T}(LEP)$", r"$m_{T}(LEP,MET)$", r"$m_{T}(LEP,MET)$"],
    "titleAdd":["", "", r", $NN>NNCUT$"],
    "selStr":["", "", "_nnpass"],
    "sels":[{}, {}, {"nn_disc":["NNCUT",None]}],
    "addOpts":[{}, {}, {}],
    "blindSel":["", "", ""],
    "addOptsBase":"",
  },
  "met_nn_kin":{
    "varName":["met_pt", "met_pt", "nn_disc", "nn_disc", "h_pt", "h_pt"],
    "varLabel":[r"PUPPI $p_{T}^{miss}$", r"PUPPI $p_{T}^{miss}$", "NN", "NN", r"$p_{T}(h)$", r"$p_{T}(h)$"],
    "titleAdd":["", r", $NN>NNCUT$", "", "", "", r", $NN>NNCUT$"],
    "selStr":["", "_nnpass", "", "_zoom", "", "_nnpass"],
    "sels":[{}, {"nn_disc":["NNCUT",None]}, {}, {"nn_disc":[0.9,None]}, {}, {"nn_disc":["NNCUT",None]}],
    "addOpts":[{}, {}, {"rebin":[0.,0.1,0.5,0.8,0.9,1.00001],"dosigrat":True}, {"xlimits":[0.9,1.],"dosigrat":True,"xexp":True,"rebin":[0.9,0.95,0.99,0.995,0.999,0.9995,0.9999,1.00001]}, {"dosigrat":True}, {"dosigrat":True}],
    "blindSel":["", "", ["None", hadLepCut], ["None", hadLepCut], "", ""],
    "addOptsBase":"",
  },
}

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
  
looseList = ["hadhad_signal", "hadhad_signal_met", "hadhad_cr_b", "hadhad_cr_b_met", "hadhad_cr_b_mu", "hadhad_cr_b_mu_iso", "hadhad_cr_mu", "hadhad_cr_mu_iso", "hadel_signal", "hadmu_signal", "hadel_cr_b", "hadel_cr_w", "hadel_cr_qcd", "hadmu_cr_b", "hadmu_cr_w", "hadmu_cr_qcd"]
lepLooseNN = "0.9"
hadLooseNN = "0.995"

def addPassFail(origList, region, ptList, keepBlind=True):
  nnCats = ["pass", "fail"]
  nnTitle = ["$NN>NNCUT$", "$NN<NNCUT$"]
  nnSel = [["NNCUT", None], [None, "NNCUT"]]
  nnBlind = [["None", 110] if keepBlind else "", ""]
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
      "addOpts":[{} for ic in range(len(nnCats))],
      "blindSel":[nnBlind[ic] for ic in range(len(nnCats))],
    }
  newList = {obj:origList[obj]+passFailList[obj] if obj in passFailList else origList[obj] for obj in origList}
  return newList

mcSamples = samplelists.getSamplesMC(year)
dataSamples = samplelists.getSamplesData(year)

cutTitleMap = {
  "cutflow_hadhad":r'\tau_{h}\tau_{h}',
  "cutflow_hadhad_met":r'\tau_{h}\tau_{h},~high~MET',
  "cutflow_hadhad_cr_b":r'\tau_{h}\tau_{h}~CR~b',
  "cutflow_hadhad_cr_b_met":r'\tau_{h}\tau_{h},~high~MET,~CR~b',
  "cutflow_hadhad_cr_b_mu":r'\tau_{h}\tau_{h}~CR~\mu~(b)',
  "cutflow_hadhad_cr_b_mu_iso":r'\tau_{h}\tau_{h}~CR~\mu~(Iso,~b)',
  "cutflow_hadhad_cr_mu":r'\tau_{h}\tau_{h}~CR~\mu',
  "cutflow_hadhad_cr_mu_iso":r'\tau_{h}\tau_{h}~CR~\mu~(Iso)',
  "cutflow_hadmu":r'\mu\tau_{h}',
  "cutflow_hadel":r'e\tau_{h}',
  "cutflow_hadmu_cr_b":r'\mu\tau_{h}~CR~b',
  "cutflow_hadel_cr_b":r'e\tau_{h}~CR~b',
  "cutflow_hadmu_cr_qcd":r'\mu\tau_{h}~CR~QCD',
  "cutflow_hadel_cr_qcd":r'e\tau_{h}~CR~QCD',
  "cutflow_hadmu_cr_w":r'\mu\tau_{h}~CR~W',
  "cutflow_hadel_cr_w":r'e\tau_{h}~CR~W',
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

plot_cutflow.getPlots(cutflowArgs)

class StackArgs:
  def __init__(self, hists, tag, savetag, var, varlabel, title, regions, hist, lumi=50., sel='', sigsel='', xlimits='', xlog=False, xexp=False, blind='', solo=False, sigscale=50, sigstack=False, noratio=False, dosigrat=False, overflow='allnan', rebin=[1], norm=None, density=False, sample='all', comp=False, compsel='', complabels='', defcolors=True, verbose=False, dirname=None):
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

region_opts = {
  "hadhad_signal":{
    "titleBase":r'$\tau_{h}\tau_{h}$',
    "lep":"",
    "histList":histListBase,
    "blind":True,
    "ptList":[0],#[450,1200],
    "cutNN":hadHadCut,
    "specSel":{'jetmet_dphi':[None,1.6],'met_pt':[100.,150.]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadhad_signal_met":{
    "titleBase":r'$\tau_{h}\tau_{h}$, high MET',
    "lep":"",
    "histList":histListBase,
    "blind":True,
    "ptList":[0],#[300,450],
    "cutNN":hadHadCut,
    "specSel":{'jetmet_dphi':[None,1.6],'met_pt':[150.,None]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadhad_cr_b":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR b',
    "lep":"",
    "histList":histListBase,
    "blind":True,
    "ptList":[0],#[450,1200],
    "cutNN":hadHadCut,
    "specSel":{'jetmet_dphi':[None,1.6],'met_pt':[50.,150.]},
    "sigSel":{},
    "saveLabel":"",
  },
  "hadhad_cr_b_met":{
    "titleBase":r'$\tau_{h}\tau_{h}$, high MET CR b',
    "lep":"",
    "histList":histListBase,
    "blind":True,
    "ptList":[0],#[300,450],
    "cutNN":hadHadCut,
    "specSel":{'jetmet_dphi':[None,1.6],'met_pt':[50.,None]},
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
    "specSel":{'jetmet_dphi':[None,1.6],'met_pt':[50.,None]},
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
    "specSel":{'jetmet_dphi':[None,1.6],'met_pt':[50.,None]},
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
    "specSel":{'jetmet_dphi':[None,1.6],'met_pt':[50.,None]},
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
    "specSel":{'jetmet_dphi':[None,1.6],'met_pt':[50.,None]},
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
    "specSel":{"antilep":[0.5,None],'met_pt':[50.,None]},
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
    "specSel":{"antilep":[0.5,None],'met_pt':[50.,None]},
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
    "specSel":{"antilep":[0.5,None],'met_pt':[50.,None]},
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
    "specSel":{"antilep":[0.5,None],'met_pt':[50.,None]},
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
    "specSel":{"antilep":[0.5,None],'met_pt':[50.,None]},
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
    "specSel":{"antilep":[0.5,None],'met_pt':[50.,None]},
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
    "specSel":{"antilep":[0.5,None],'met_pt':[50.,None]},
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
    "specSel":{"antilep":[0.5,None],'met_pt':[50.,None]},
    "sigSel":{},
    "saveLabel":"",
  },
}


def doPlotting(infileList,tag,lumi,region,hist,plotList,regopts,add_opts={},sample='all', keepBlind=True):
  if (hist==srHist): 
    thePlotList = StringRemap(addPassFail(plotList[hist], region, regopts["ptList"], keepBlind=True),regopts["lep"],regopts["cutNN"])
    #thePlotList = StringRemap(plotList[hist],regopts["lep"],regopts["cutNN"])
  else:
    thePlotList = StringRemap(plotList[hist],regopts["lep"],regopts["cutNN"])
  add_opts = StringRemap(add_opts,regopts["lep"],regopts["cutNN"])
  minStackArgs = StackArgs(infileList, tag, "", "", "", "", region, hist, lumi)
  cofhist = plot_stack.getPlots(minStackArgs,True)
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
    comb_opt_dict.update({**thePlotList["addOpts"][ivar],**{"blind":blindStr}})
    if 'comp_cut' in comb_opt_dict:
      comb_opt_dict['comp_cut'] = [[[float(cs[0]) if cs[0] is not None else None, float(cs[1]) if cs[1] is not None else None] for cs in v] for v in comb_opt_dict['comp_cut']]
    thePlotList["sels"][ivar] = {k:[float(thePlotList["sels"][ivar][k][0]) if thePlotList["sels"][ivar][k][0] is not None else None, float(thePlotList["sels"][ivar][k][1]) if thePlotList["sels"][ivar][k][1] is not None else None] for k in thePlotList["sels"][ivar]}
    ndict = regopts["specSel"].copy()
    ndict.update(thePlotList["sels"][ivar])
    print(ndict)
    print(comb_opt_dict)
    print(regopts["sigSel"])
    plot_stack.drawStack(cofhist,hist,thePlotList["varName"][ivar],thePlotList["varLabel"][ivar],regopts["titleBase"]+thePlotList["titleAdd"][ivar],sample,lumi,ndict,regopts["sigSel"],[region],region+thePlotList["selStr"][ivar]+regopts["saveLabel"],**comb_opt_dict)

  os.chdir(pwd)

print('Nominal')
for region in regionList:
  for hist in region_opts[region]["histList"]:
    print(region, hist)
    doPlotting(infileList,tag,lumi,region,hist,plotList,region_opts[region])
print("Done")

#region_opts["hadhad_signal"]["sigSel"]       = {"h_pt":[400,1200]}
#region_opts["hadhad_signal_met"]["sigSel"]   = {"h_pt":[400,1200]}
#region_opts["hadhad_cr_mu"]["sigSel"]        = {"h_pt":[400,1200]}
#region_opts["hadhad_cr_mu_iso"]["sigSel"]   = {"h_pt":[400,1200]}
#region_opts["hadhad_cr_b_mu"]["sigSel"]      = {"h_pt":[400,1200]}
#region_opts["hadhad_cr_b_mu_iso"]["sigSel"] = {"h_pt":[400,1200]}
#region_opts["hadmu_signal"]["sigSel"]        = {"h_pt":[400,1200]}
#region_opts["hadel_signal"]["sigSel"]        = {"h_pt":[400,1200]}
#region_opts["hadmu_cr_qcd"]["sigSel"]        = {"h_pt":[400,1200]}
#region_opts["hadel_cr_qcd"]["sigSel"]        = {"h_pt":[400,1200]}
#region_opts["hadmu_cr_b"]["sigSel"]          = {"h_pt":[400,1200]}
#region_opts["hadel_cr_b"]["sigSel"]          = {"h_pt":[400,1200]}
#region_opts["hadmu_cr_w"]["sigSel"]          = {"h_pt":[400,1200]}
#region_opts["hadel_cr_w"]["sigSel"]          = {"h_pt":[400,1200]}

region_opts["hadhad_signal"]["specSel"]["h_pt"] = [400,1200]
region_opts["hadhad_signal_met"]["specSel"]["h_pt"] = [400,1200]
region_opts["hadhad_cr_mu"]["specSel"]["h_pt"] = [400,1200]
region_opts["hadhad_cr_mu_iso"]["specSel"]["h_pt"] = [400,1200]
region_opts["hadhad_cr_b_mu"]["specSel"]["h_pt"] = [400,1200]
region_opts["hadhad_cr_b_mu_iso"]["specSel"]["h_pt"] = [400,1200]
region_opts["hadmu_signal"]["specSel"]["h_pt"] = [400,1200]
region_opts["hadel_signal"]["specSel"]["h_pt"] = [400,1200]
region_opts["hadmu_cr_qcd"]["specSel"]["h_pt"] = [400,1200]
region_opts["hadel_cr_qcd"]["specSel"]["h_pt"] = [400,1200]
region_opts["hadmu_cr_b"]["specSel"]["h_pt"] = [400,1200]
region_opts["hadel_cr_b"]["specSel"]["h_pt"] = [400,1200]
region_opts["hadmu_cr_w"]["specSel"]["h_pt"] = [400,1200]
region_opts["hadel_cr_w"]["specSel"]["h_pt"] = [400,1200]

region_opts["hadhad_signal"]["titleBase"]       = region_opts["hadhad_signal"]["titleBase"]       + ", $p_{T}(h)>400$"
region_opts["hadhad_signal_met"]["titleBase"]   = region_opts["hadhad_signal_met"]["titleBase"]   + ", $p_{T}(h)>400$"
region_opts["hadhad_cr_mu"]["titleBase"]        = region_opts["hadhad_cr_mu"]["titleBase"]        + ", $p_{T}(h)>400$"
region_opts["hadhad_cr_mu_iso"]["titleBase"]   = region_opts["hadhad_cr_mu_iso"]["titleBase"]   + ", $p_{T}(h)>400$"
region_opts["hadhad_cr_b_mu"]["titleBase"]      = region_opts["hadhad_cr_b_mu"]["titleBase"]      + ", $p_{T}(h)>400$"
region_opts["hadhad_cr_b_mu_iso"]["titleBase"] = region_opts["hadhad_cr_b_mu_iso"]["titleBase"] + ", $p_{T}(h)>400$"
region_opts["hadmu_signal"]["titleBase"]        = region_opts["hadmu_signal"]["titleBase"]        + ", $p_{T}(h)>400$"
region_opts["hadel_signal"]["titleBase"]        = region_opts["hadel_signal"]["titleBase"]        + ", $p_{T}(h)>400$"
region_opts["hadmu_cr_qcd"]["titleBase"]        = region_opts["hadmu_cr_qcd"]["titleBase"]        + ", $p_{T}(h)>400$"
region_opts["hadel_cr_qcd"]["titleBase"]        = region_opts["hadel_cr_qcd"]["titleBase"]        + ", $p_{T}(h)>400$"
region_opts["hadmu_cr_b"]["titleBase"]          = region_opts["hadmu_cr_b"]["titleBase"]          + ", $p_{T}(h)>400$"
region_opts["hadel_cr_b"]["titleBase"]          = region_opts["hadel_cr_b"]["titleBase"]          + ", $p_{T}(h)>400$"
region_opts["hadmu_cr_w"]["titleBase"]          = region_opts["hadmu_cr_w"]["titleBase"]          + ", $p_{T}(h)>400$"
region_opts["hadel_cr_w"]["titleBase"]          = region_opts["hadel_cr_w"]["titleBase"]          + ", $p_{T}(h)>400$"

region_opts["hadhad_signal"]["saveLabel"]       = "_sigsel"
region_opts["hadhad_signal_met"]["saveLabel"]   = "_sigsel"
region_opts["hadhad_cr_mu"]["saveLabel"]        = "_sigsel"
region_opts["hadhad_cr_mu_iso"]["saveLabel"]   = "_sigsel"
region_opts["hadhad_cr_b_mu"]["saveLabel"]      = "_sigsel"
region_opts["hadhad_cr_b_mu_iso"]["saveLabel"] = "_sigsel"
region_opts["hadmu_signal"]["saveLabel"]        = "_sigsel"
region_opts["hadel_signal"]["saveLabel"]        = "_sigsel"
region_opts["hadmu_cr_qcd"]["saveLabel"]        = "_sigsel"
region_opts["hadel_cr_qcd"]["saveLabel"]        = "_sigsel"
region_opts["hadmu_cr_b"]["saveLabel"]          = "_sigsel"
region_opts["hadel_cr_b"]["saveLabel"]          = "_sigsel"
region_opts["hadmu_cr_w"]["saveLabel"]          = "_sigsel"
region_opts["hadel_cr_w"]["saveLabel"]          = "_sigsel"

print('Nominal w/ signal selection')
for region in regionList:
  for hist in region_opts[region]["histList"]:
    print(region, hist)
    doPlotting(infileList,tag,lumi,region,hist,plotList,region_opts[region])
print("Done")

region_opts = {
  "hadhad_signal":{
    "titleBase":r'$\tau_{h}\tau_{h}$ (invert $\Delta\phi(j,MET)$)',
    "lep":"",
    "histList":histListBase,
    "blind":True,
    "ptList":[0],#[450,1200],
    "cutNN":hadHadCut,
    "specSel":{'jetmet_dphi':[1.6,None],'met_pt':[100.,150.]},
    "sigSel":{},
    "saveLabel":"_invdphisel",
  },
  "hadhad_signal_met":{
    "titleBase":r'$\tau_{h}\tau_{h}$, high MET (invert $\Delta\phi(j,MET)$)',
    "lep":"",
    "histList":histListBase,
    "blind":True,
    "ptList":[0],#[300,450],
    "cutNN":hadHadCut,
    "specSel":{'jetmet_dphi':[1.6,None],'met_pt':[50.,None]},
    "sigSel":{},
    "saveLabel":"_invdphisel",
  },
  "hadhad_cr_b":{
    "titleBase":r'$\tau_{h}\tau_{h}$ (invert $\Delta\phi(j,MET)$) CR b',
    "lep":"",
    "histList":histListBase,
    "blind":True,
    "ptList":[0],#[450,1200],
    "cutNN":hadHadCut,
    "specSel":{'jetmet_dphi':[1.6,None],'met_pt':[50.,150.]},
    "sigSel":{},
    "saveLabel":"_invdphisel",
  },
  "hadhad_cr_b_met":{
    "titleBase":r'$\tau_{h}\tau_{h}$, high MET (invert $\Delta\phi(j,MET)$) CR b',
    "lep":"",
    "histList":histListBase,
    "blind":True,
    "ptList":[0],#[300,450],
    "cutNN":hadHadCut,
    "specSel":{'jetmet_dphi':[1.6,None],'met_pt':[50.,None]},
    "sigSel":{},
    "saveLabel":"_invdphisel",
  },
  "hadhad_cr_b_mu":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (b) (invert $\Delta\phi(j,MET)$)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'jetmet_dphi':[1.6,None],'met_pt':[50.,None]},
    "sigSel":{},
    "saveLabel":"_invdphisel",
  },
  "hadhad_cr_b_mu_iso":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso, b) (invert $\Delta\phi(j,MET)$)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'jetmet_dphi':[1.6,None],'met_pt':[50.,None]},
    "sigSel":{},
    "saveLabel":"_invdphisel",
  },
  "hadhad_cr_mu":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (invert $\Delta\phi(j,MET)$)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "histList":histListBase,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'jetmet_dphi':[1.6,None],'met_pt':[50.,None]},
    "sigSel":{},
    "saveLabel":"_invdphisel",
  },
  "hadhad_cr_mu_iso":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso) (invert $\Delta\phi(j,MET)$)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'jetmet_dphi':[1.6,None],'met_pt':[50.,None]},
    "sigSel":{},
    "saveLabel":"_invdphisel",
  },
  "hadmu_signal":{
    "titleBase":r'$\mu\tau_{h}$ (no AntiLep)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":True,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{"antilep":[None,0.5],'met_pt':[50.,None]},
    "sigSel":{},
    "saveLabel":"_noantilep",
  },
  "hadel_signal":{
    "titleBase":r'$e\tau_{h}$ (no AntiLep)',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":True,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{"antilep":[None,0.5],'met_pt':[50.,None]},
    "sigSel":{},
    "saveLabel":"_noantilep",
  },
  "hadmu_cr_b":{
    "titleBase":r'$\mu\tau_{h}$ CR b (no AntiLep)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{"antilep":[None,0.5],'met_pt':[50.,None]},
    "sigSel":{},
    "saveLabel":"_noantilep",
  },
  "hadel_cr_b":{
    "titleBase":r'$e\tau_{h}$ CR b (no AntiLep)',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{"antilep":[None,0.5],'met_pt':[50.,None]},
    "sigSel":{},
    "saveLabel":"_noantilep",
  },
  "hadmu_cr_qcd":{
    "titleBase":r'$\mu\tau_{h}$ CR QCD (no AntiLep)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{"antilep":[None,0.5],'met_pt':[50.,None]},
    "sigSel":{},
    "saveLabel":"_noantilep",
  },
  "hadel_cr_qcd":{
    "titleBase":r'$e\tau_{h}$ CR QCD (no AntiLep)',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{"antilep":[None,0.5],'met_pt':[50.,None]},
    "sigSel":{},
    "saveLabel":"_noantilep",
  },
  "hadmu_cr_w":{
    "titleBase":r'$\mu\tau_{h}$ CR W (no AntiLep)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{"antilep":[None,0.5],'met_pt':[50.,None]},
    "sigSel":{},
    "saveLabel":"_noantilep",
  },
  "hadel_cr_w":{
    "titleBase":r'$e\tau_{h}$ CR W (no AntiLep)',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{"antilep":[None,0.5],'met_pt':[50.,None]},
    "sigSel":{},
    "saveLabel":"_noantilep",
  },
}

print('NoAntiLep')
for region in regionList:
  for hist in region_opts[region]["histList"]:
    print(region, hist)
    doPlotting(infileList,tag,lumi,region,hist,plotList,region_opts[region],keepBlind=False)
print("Done")

regionList = [
  "hadhad_signal",
  "hadhad_signal_met",
  "hadhad_cr_b",
  "hadhad_cr_b_met",
  "hadhad_cr_mu",
  "hadhad_cr_mu_iso",
  "hadhad_cr_b_mu",
  "hadhad_cr_b_mu_iso",
]

region_opts = {
  "hadhad_signal":{
    "titleBase":r'$\tau_{h}\tau_{h}$ (AntiLep)',
    "lep":"",
    "histList":histListBase,
    "blind":True,
    "ptList":[0],#[450,1200],
    "cutNN":hadHadCut,
    "specSel":{'antilep':[0.5,None],'met_pt':[100.,150.],'jetmet_dphi':[None,1.6]},
    "sigSel":{},
    "saveLabel":"_antilep",
  },
  "hadhad_signal_met":{
    "titleBase":r'$\tau_{h}\tau_{h}$, high MET (AntiLep)',
    "lep":"",
    "histList":histListBase,
    "blind":True,
    "ptList":[0],#[300,450],
    "cutNN":hadHadCut,
    "specSel":{'antilep':[0.5,None],'met_pt':[50.,None],'jetmet_dphi':[None,1.6]},
    "sigSel":{},
    "saveLabel":"_antilep",
  },
  "hadhad_cr_b":{
    "titleBase":r'$\tau_{h}\tau_{h}$ (AntiLep) CR b',
    "lep":"",
    "histList":histListBase,
    "blind":True,
    "ptList":[0],#[450,1200],
    "cutNN":hadHadCut,
    "specSel":{'antilep':[0.5,None],'met_pt':[50.,150.],'jetmet_dphi':[None,1.6]},
    "sigSel":{},
    "saveLabel":"_antilep",
  },
  "hadhad_cr_b_met":{
    "titleBase":r'$\tau_{h}\tau_{h}$, high MET (AntiLep) CR b',
    "lep":"",
    "histList":histListBase,
    "blind":True,
    "ptList":[0],#[300,450],
    "cutNN":hadHadCut,
    "specSel":{'antilep':[0.5,None],'met_pt':[50.,None],'jetmet_dphi':[None,1.6]},
    "sigSel":{},
    "saveLabel":"_antilep",
  },
  "hadhad_cr_b_mu":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (b) (AntiLep)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'antilep':[0.5,None],'met_pt':[50.,None],'jetmet_dphi':[None,1.6]},
    "sigSel":{},
    "saveLabel":"_antilep",
  },
  "hadhad_cr_b_mu_iso":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso, b) (AntiLep)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'antilep':[0.5,None],'met_pt':[50.,None],'jetmet_dphi':[None,1.6]},
    "sigSel":{},
    "saveLabel":"_antilep",
  },
  "hadhad_cr_mu":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (AntiLep)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "histList":histListBase,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'antilep':[0.5,None],'met_pt':[50.,None],'jetmet_dphi':[None,1.6]},
    "sigSel":{},
    "saveLabel":"_antilep",
  },
  "hadhad_cr_mu_iso":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso) (AntiLep)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'antilep':[0.5,None],'met_pt':[50.,None],'jetmet_dphi':[None,1.6]},
    "sigSel":{},
    "saveLabel":"_antilep",
  },
}

print('AntiLep')
for region in regionList:
  for hist in region_opts[region]["histList"]:
    print(region, hist)
    #doPlotting(infileList,tag,lumi,region,hist,plotList,region_opts[region],keepBlind=False)
print("Done")

regionList = [
  "hadhad_signal",
  "hadhad_signal_met",
  "hadhad_cr_b",
  "hadhad_cr_b_met",
  "hadhad_cr_b_mu",
  "hadhad_cr_b_mu_iso",
  "hadhad_cr_mu",
  "hadhad_cr_mu_iso",
]

region_opts = {
  "hadhad_signal":{
    "titleBase":r'$\tau_{h}\tau_{h}$ ($N_{2}^{DDT}<0$)',
    "lep":"",
    "histList":histListBase,
    "blind":True,
    "ptList":[0],#[450,1200],
    "cutNN":hadHadCut,
    "specSel":{'n2ddt':[None,0.],'jetmet_dphi':[None,1.6],'met_pt':[100.,150.]},
    "sigSel":{},
    "saveLabel":"_n2sel",
  },
  "hadhad_signal_met":{
    "titleBase":r'$\tau_{h}\tau_{h}$, high MET ($N_{2}^{DDT}<0$)',
    "lep":"",
    "histList":histListBase,
    "blind":True,
    "ptList":[0],#[300,450],
    "cutNN":hadHadCut,
    "specSel":{'n2ddt':[None,0.],'jetmet_dphi':[None,1.6],'met_pt':[50.,None]},
    "sigSel":{},
    "saveLabel":"_n2sel",
  },
  "hadhad_cr_b":{
    "titleBase":r'$\tau_{h}\tau_{h}$ ($N_{2}^{DDT}<0$) CR b',
    "lep":"",
    "histList":histListBase,
    "blind":True,
    "ptList":[0],#[450,1200],
    "cutNN":hadHadCut,
    "specSel":{'n2ddt':[None,0.],'jetmet_dphi':[None,1.6],'met_pt':[50.,150.]},
    "sigSel":{},
    "saveLabel":"_n2sel",
  },
  "hadhad_cr_b_met":{
    "titleBase":r'$\tau_{h}\tau_{h}$, high MET ($N_{2}^{DDT}<0$) CR b',
    "lep":"",
    "histList":histListBase,
    "blind":True,
    "ptList":[0],#[300,450],
    "cutNN":hadHadCut,
    "specSel":{'n2ddt':[None,0.],'jetmet_dphi':[None,1.6],'met_pt':[50.,None]},
    "sigSel":{},
    "saveLabel":"_n2sel",
  },
  "hadhad_cr_b_mu":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (b) ($N_{2}^{DDT}<0$)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'n2ddt':[None,0.],'jetmet_dphi':[None,1.6],'met_pt':[50.,None]},
    "sigSel":{},
    "saveLabel":"_n2sel",
  },
  "hadhad_cr_b_mu_iso":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso, b) ($N_{2}^{DDT}<0$)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'n2ddt':[None,0.],'jetmet_dphi':[None,1.6],'met_pt':[50.,None]},
    "sigSel":{},
    "saveLabel":"_n2sel",
  },
  "hadhad_cr_mu":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ ($N_{2}^{DDT}<0$)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "histList":histListBase,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'n2ddt':[None,0.],'jetmet_dphi':[None,1.6],'met_pt':[50.,None]},
    "sigSel":{},
    "saveLabel":"_n2sel",
  },
  "hadhad_cr_mu_iso":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso) ($N_{2}^{DDT}<0$)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'n2ddt':[None,0.],'jetmet_dphi':[None,1.6],'met_pt':[50.,None]},
    "sigSel":{},
    "saveLabel":"_n2sel",
  },
}

print('N2')
for region in regionList:
  for hist in region_opts[region]["histList"]:
    print(region, hist)
    #doPlotting(infileList,tag,lumi,region,hist,plotList,region_opts[region])
print("Done")

region_opts = {
  "hadhad_signal":{
    "titleBase":r'$\tau_{h}\tau_{h}$ ($N_{2}^{DDT}>0$)',
    "lep":"",
    "histList":histListBase,
    "blind":True,
    "ptList":[0],#[450,1200],
    "cutNN":hadHadCut,
    "specSel":{'n2ddt':[0.,None],'jetmet_dphi':[None,1.6],'met_pt':[100.,150.]},
    "sigSel":{},
    "saveLabel":"_n2cut",
  },
  "hadhad_signal_met":{
    "titleBase":r'$\tau_{h}\tau_{h}$, high MET ($N_{2}^{DDT}>0$)',
    "lep":"",
    "histList":histListBase,
    "blind":True,
    "ptList":[0],#[300,450],
    "cutNN":hadHadCut,
    "specSel":{'n2ddt':[0.,None],'jetmet_dphi':[None,1.6],'met_pt':[50.,None]},
    "sigSel":{},
    "saveLabel":"_n2cut",
  },
  "hadhad_cr_b":{
    "titleBase":r'$\tau_{h}\tau_{h}$ ($N_{2}^{DDT}>0$) CR b',
    "lep":"",
    "histList":histListBase,
    "blind":True,
    "ptList":[0],#[450,1200],
    "cutNN":hadHadCut,
    "specSel":{'n2ddt':[0.,None],'jetmet_dphi':[None,1.6],'met_pt':[50.,150.]},
    "sigSel":{},
    "saveLabel":"_n2cut",
  },
  "hadhad_cr_b_met":{
    "titleBase":r'$\tau_{h}\tau_{h}$, high MET ($N_{2}^{DDT}>0$) CR b',
    "lep":"",
    "histList":histListBase,
    "blind":True,
    "ptList":[0],#[300,450],
    "cutNN":hadHadCut,
    "specSel":{'n2ddt':[0.,None],'jetmet_dphi':[None,1.6],'met_pt':[50.,None]},
    "sigSel":{},
    "saveLabel":"_n2cut",
  },
  "hadhad_cr_b_mu":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (b) ($N_{2}^{DDT}>0$)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'n2ddt':[0.,None],'jetmet_dphi':[None,1.6],'met_pt':[50.,None]},
    "sigSel":{},
    "saveLabel":"_n2cut",
  },
  "hadhad_cr_b_mu_iso":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso, b) ($N_{2}^{DDT}>0$)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'n2ddt':[0.,None],'jetmet_dphi':[None,1.6],'met_pt':[50.,None]},
    "sigSel":{},
    "saveLabel":"_n2cut",
  },
  "hadhad_cr_mu":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ ($N_{2}^{DDT}>0$)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "histList":histListBase,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'n2ddt':[0.,None],'jetmet_dphi':[None,1.6],'met_pt':[50.,None]},
    "sigSel":{},
    "saveLabel":"_n2cut",
  },
  "hadhad_cr_mu_iso":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso) ($N_{2}^{DDT}>0$)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'n2ddt':[0.,None],'jetmet_dphi':[None,1.6],'met_pt':[50.,None]},
    "sigSel":{},
    "saveLabel":"_n2cut",
  },
}

print('N2 Inv')
for region in regionList:
  for hist in region_opts[region]["histList"]:
    print(region, hist)
    #doPlotting(infileList,tag,lumi,region,hist,plotList,region_opts[region])
print("Done")

region_opts = {
  "hadhad_signal":{
    "titleBase":r'$\tau_{h}\tau_{h}$ (low MET)',
    "lep":"",
    "histList":histListBase,
    "blind":True,
    "ptList":[0],#[450,1200],
    "cutNN":hadHadCut,
    "specSel":{'jetmet_dphi':[None,1.6],'met_pt':[None,100.]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadhad_cr_b":{
    "titleBase":r'$\tau_{h}\tau_{h}$ (low MET) CR b',
    "lep":"",
    "histList":histListBase,
    "blind":True,
    "ptList":[0],#[450,1200],
    "cutNN":hadHadCut,
    "specSel":{'jetmet_dphi':[None,1.6],'met_pt':[None,50.]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadhad_cr_b_mu":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (b) (low MET)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'jetmet_dphi':[None,1.6],'met_pt':[None,50.]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadhad_cr_b_mu_iso":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso, b) (low MET)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'jetmet_dphi':[None,1.6],'met_pt':[None,50.]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadhad_cr_mu":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (low MET)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "histList":histListBase,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'jetmet_dphi':[None,1.6],'met_pt':[None,50.]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadhad_cr_mu_iso":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso) (low MET)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'jetmet_dphi':[None,1.6],'met_pt':[None,50.]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadmu_signal":{
    "titleBase":r'$\mu\tau_{h}$ (low MET)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":True,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{"antilep":[0.5,None],'met_pt':[None,50.]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadel_signal":{
    "titleBase":r'$e\tau_{h}$ (low MET)',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":True,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{"antilep":[0.5,None],'met_pt':[None,50.]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadmu_cr_b":{
    "titleBase":r'$\mu\tau_{h}$ CR b (low MET)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{"antilep":[0.5,None],'met_pt':[None,50.]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadel_cr_b":{
    "titleBase":r'$e\tau_{h}$ CR b (low MET)',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{"antilep":[0.5,None],'met_pt':[None,50.]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadmu_cr_qcd":{
    "titleBase":r'$\mu\tau_{h}$ CR QCD (low MET)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{"antilep":[0.5,None],'met_pt':[None,50.]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadel_cr_qcd":{
    "titleBase":r'$e\tau_{h}$ CR QCD (low MET)',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{"antilep":[0.5,None],'met_pt':[None,50.]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadmu_cr_w":{
    "titleBase":r'$\mu\tau_{h}$ CR W (low MET)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{"antilep":[0.5,None],'met_pt':[None,50.]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
  "hadel_cr_w":{
    "titleBase":r'$e\tau_{h}$ CR W (low MET)',
    "lep":'e',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[300, 400, 1200],
    "cutNN":hadLepCut,
    "specSel":{"antilep":[0.5,None],'met_pt':[None,50.]},
    "sigSel":{},
    "saveLabel":"_lowMET",
  },
}

regionList = [
  "hadhad_signal",
  "hadhad_cr_b",
  "hadhad_cr_mu",
  "hadhad_cr_mu_iso",
  "hadhad_cr_b_mu",
  "hadhad_cr_b_mu_iso",
  "hadmu_signal",
  "hadel_signal",
  "hadmu_cr_qcd",
  "hadel_cr_qcd",
  "hadmu_cr_b",
  "hadel_cr_b",
  "hadmu_cr_w",
  "hadel_cr_w",
]

print('Low MET')
for region in regionList:
  for hist in region_opts[region]["histList"]:
    print(region, hist)
    #doPlotting(infileList,tag,lumi,region,hist,plotList,region_opts[region],keepBlind=False)
print("Done")

region_opts = {
  "hadhad_cr_b_mu":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (b) (high MET)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'jetmet_dphi':[None,1.6],'met_pt':[150.,None]},
    "sigSel":{},
    "saveLabel":"_highMET",
  },
  "hadhad_cr_b_mu_iso":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso, b) (high MET)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'jetmet_dphi':[None,1.6],'met_pt':[150.,None]},
    "sigSel":{},
    "saveLabel":"_highMET",
  },
  "hadhad_cr_mu":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (high MET)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "histList":histListBase,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'jetmet_dphi':[None,1.6],'met_pt':[150.,None]},
    "sigSel":{},
    "saveLabel":"_highMET",
  },
  "hadhad_cr_mu_iso":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso) (high MET)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'jetmet_dphi':[None,1.6],'met_pt':[150.,None]},
    "sigSel":{},
    "saveLabel":"_highMET",
  },
}

regionList = [
  "hadhad_cr_mu",
  "hadhad_cr_mu_iso",
  "hadhad_cr_b_mu",
  "hadhad_cr_b_mu_iso",
]

print('High MET')
for region in regionList:
  for hist in region_opts[region]["histList"]:
    print(region, hist)
    #doPlotting(infileList,tag,lumi,region,hist,plotList,region_opts[region],keepBlind=False)
print("Done")

region_opts = {
  "hadhad_cr_b_mu":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (b) (high MET, invert $\Delta\phi(j,MET)$)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'jetmet_dphi':[1.6,None],'met_pt':[150.,None]},
    "sigSel":{},
    "saveLabel":"_highMET_invdphisel",
  },
  "hadhad_cr_b_mu_iso":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso, b) (high MET, invert $\Delta\phi(j,MET)$)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'jetmet_dphi':[1.6,None],'met_pt':[150.,None]},
    "sigSel":{},
    "saveLabel":"_highMET_invdphisel",
  },
  "hadhad_cr_mu":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (high MET, invert $\Delta\phi(j,MET)$)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "histList":histListBase,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'jetmet_dphi':[1.6,None],'met_pt':[150.,None]},
    "sigSel":{},
    "saveLabel":"_highMET_invdphisel",
  },
  "hadhad_cr_mu_iso":{
    "titleBase":r'$\tau_{h}\tau_{h}$ CR $\mu$ (Iso) (high MET, invert $\Delta\phi(j,MET)$)',
    "lep":'\mu',
    "histList":histListBase+lepHists,
    "blind":False,
    "ptList":[0],#[400,1200],
    "cutNN":hadHadCut,
    "specSel":{'jetmet_dphi':[1.6,None],'met_pt':[150.,None]},
    "sigSel":{},
    "saveLabel":"_highMET_invdphisel",
  },
}

regionList = [
  "hadhad_cr_mu",
  "hadhad_cr_mu_iso",
  "hadhad_cr_b_mu",
  "hadhad_cr_b_mu_iso",
]

print('High MET, inv dPhi')
for region in regionList:
  for hist in region_opts[region]["histList"]:
    print(region, hist)
    #doPlotting(infileList,tag,lumi,region,hist,plotList,region_opts[region],keepBlind=False)
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
