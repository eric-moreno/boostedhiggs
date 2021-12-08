#!/bin/bash

INFILE=../condor/Feb05_LepID/hists_sum 
TAG=LepID_Feb05

for VETO in hv ll mm
do
  case ${VETO} in
    hv)
      TITLE='High $p_{T}$ $\mu$ / Veto $e$ ID'
      ;;
    ll)
      TITLE='Loose $\mu$ / Loose $e$ ID'
      ;;
    mm)
      TITLE='Medium $\mu$ / Medium $e$ ID'
      ;;
  esac
  #for VAR in lep_pt lep_eta jet_pt jet_eta
  for VAR in lep_pt
  do
    case ${VAR} in
      lep_pt)
        VAR_LABEL='$p_{T}(lepton)$'
        HISTNAME=lep_kin
        REBIN=2
        SEL_VAR=miso
        SEL_RANGE="None 0.1"
        SEL_STR=_lt_0p1
        SEL_TEXT="miniIso < 0.1"
        ;;
      lep_eta)
        VAR_LABEL='$\eta(lepton)$'
        HISTNAME=lep_kin
        REBIN=1
        SEL_VAR=miso
        SEL_RANGE="None 0.1"
        SEL_STR=_lt_0p1
        SEL_TEXT="miniIso < 0.1"
        ;;
      jet_pt)
        VAR_LABEL='$p_{T}(jet)$'
        HISTNAME=jet_kin
        REBIN=2
        SEL_VAR=miso
        SEL_RANGE="None 0.1"
        SEL_STR=_lt_0p1
        SEL_TEXT="miniIso < 0.1"
        ;;
      jet_eta)
        VAR_LABEL='$\eta(jet)$'
        HISTNAME=jet_kin
        REBIN=1
        SEL_VAR=miso
        SEL_RANGE="None 0.1"
        SEL_STR=_lt_0p1
        SEL_TEXT="miniIso < 0.1"
        ;;
      miso)
        VAR_LABEL='miniIso'
        HISTNAME=lep_kin
        REBIN=1
        SEL_VAR=lep_pt
        SEL_RANGE="50. None"
        SEL_STR=_gt_50
        SEL_TEXT='$p_{T}$ > 50'
        ;;
    esac
    #echo python plot_lepid.py --hists ${INFILE} --histname ${HISTNAME} --varname ${VAR} --varlabel "${VAR_LABEL}" --varcuts region hadel --numsel ${SEL_VAR} ${SEL_RANGE} --title "Electron ID, ${TITLE} (Veto), ${SEL_TEXT}" --label hadel_${SEL_VAR}${SEL_STR}_${VETO} --tag LepID_Feb05 --veto ${VETO} --sigsel genhtt 1.5 2.5 --rebin ${REBIN}
    python plot_lepid.py --hists ${INFILE} --histname ${HISTNAME} --varname ${VAR} --varlabel "${VAR_LABEL}" --varcuts region hadel --numsel ${SEL_VAR} ${SEL_RANGE} --title "Electron ID, ${TITLE} (Veto), ${SEL_TEXT}" --label hadel_${SEL_VAR}${SEL_STR}_${VETO} --tag LepID_Feb05 --veto ${VETO} --sigsel genhtt 1.5 2.5 --rebin ${REBIN}
    python plot_lepid.py --hists ${INFILE} --histname ${HISTNAME} --varname ${VAR} --varlabel "${VAR_LABEL}" --varcuts region hadmu --numsel ${SEL_VAR} ${SEL_RANGE} --title "Muon ID, ${TITLE} (Veto), ${SEL_TEXT}" --label hadmu_${SEL_VAR}${SELSTR}_${VETO} --tag LepID_Feb05 --veto ${VETO} --sigsel genhtt 2.5 3.5 --rebin ${REBIN}
    for SAMP in qcd wjets tt-semilep zll
    do
      python plot_lepid.py --hists ${INFILE} --histname ${HISTNAME} --varname ${VAR} --varlabel "${VAR_LABEL}" --varcuts region hadel --numsel ${SEL_VAR} ${SEL_RANGE} --title "Electron ID, ${TITLE} (Veto), ${SEL_TEXT}" --label hadel_${SEL_VAR}${SEL_STR}_${VETO} --tag LepID_Feb05 --veto ${VETO} --sigsel genhtt 1.5 2.5 --rebin ${REBIN} --sample ${SAMP}
      python plot_lepid.py --hists ${INFILE} --histname ${HISTNAME} --varname ${VAR} --varlabel "${VAR_LABEL}" --varcuts region hadmu --numsel ${SEL_VAR} ${SEL_RANGE} --title "Muon ID, ${TITLE} (Veto), ${SEL_TEXT}" --label hadmu_${SEL_VAR}${SELSTR}_${VETO} --tag LepID_Feb05 --veto ${VETO} --sigsel genhtt 2.5 3.5 --rebin ${REBIN} --sample ${SAMP}
    done
  done
  python plot_lepid_cat.py --hists ${INFILE} --histname lep_kin --varname lep_pt --numsel region hadel --title "Electron ID, ${TITLE} (Veto)" --label hadel_${VETO} --tag LepID_Feb05 --veto ${VETO} --sigsel genhtt 1.5 2.5 --densel region hadel_base
  python plot_lepid_cat.py --hists ${INFILE} --histname lep_kin --varname lep_pt --numsel region hadmu --title "Muon ID, ${TITLE} (Veto)" --label hadmu_${VETO} --tag LepID_Feb05 --veto ${VETO} --sigsel genhtt 2.5 3.5 --densel region hadmu_base
  python plot_lepid_cat.py --hists ${INFILE} --histname lep_kin --varname lep_pt --numsel region hadhad --title "${TITLE} (Veto)" --label hadhad_${VETO} --tag LepID_Feb05 --veto ${VETO} --sigsel genhtt 0.5 1.5 --densel region hadhad_base
  python plot_lepid_cat.py --hists ${INFILE} --histname lep_kin --varname lep_pt --numsel miso None 0.1 region hadel --title "Electron ID + miniIso < 0.1, ${TITLE} (Veto)" --label hadel_${VETO}_miso_lt_0p1 --tag LepID_Feb05 --veto ${VETO} --sigsel genhtt 1.5 2.5 --densel region hadel_base
  python plot_lepid_cat.py --hists ${INFILE} --histname lep_kin --varname lep_pt --numsel miso None 0.1 region hadmu --title "Muon ID + miniIso < 0.1, ${TITLE} (Veto)" --label hadmu_${VETO}_miso_lt_0p1 --tag LepID_Feb05 --veto ${VETO} --sigsel genhtt 2.5 3.5 --densel region hadmu_base
  for SAMP in qcd wjets tt-semilep zll
  do
    echo ${SAMP}
    python plot_lepid_cat.py --hists ${INFILE} --histname lep_kin --varname lep_pt --numsel region hadel --title "Electron ID, ${TITLE} (Veto)" --label hadel_${VETO} --tag LepID_Feb05 --veto ${VETO} --sigsel genhtt 1.5 2.5 --densel region hadel_base --sample ${SAMP}
    python plot_lepid_cat.py --hists ${INFILE} --histname lep_kin --varname lep_pt --numsel region hadmu --title "Muon ID, ${TITLE} (Veto)" --label hadmu_${VETO} --tag LepID_Feb05 --veto ${VETO} --sigsel genhtt 2.5 3.5 --densel region hadmu_base --sample ${SAMP}
    python plot_lepid_cat.py --hists ${INFILE} --histname lep_kin --varname lep_pt --numsel region hadhad --title "${TITLE} (Veto)" --label hadhad_${VETO} --tag LepID_Feb05 --veto ${VETO} --sigsel genhtt 0.5 1.5 --densel region hadhad_base --sample ${SAMP}
    python plot_lepid_cat.py --hists ${INFILE} --histname lep_kin --varname lep_pt --numsel miso None 0.1 region hadel --title "Electron ID + miniIso < 0.1, ${TITLE} (Veto)" --label hadel_${VETO}_miso_lt_0p1 --tag LepID_Feb05 --veto ${VETO} --sigsel genhtt 1.5 2.5 --densel region hadel_base --sample ${SAMP}
    python plot_lepid_cat.py --hists ${INFILE} --histname lep_kin --varname lep_pt --numsel miso None 0.1 region hadmu --title "Muon ID + miniIso < 0.1, ${TITLE} (Veto)" --label hadmu_${VETO}_miso_lt_0p1 --tag LepID_Feb05 --veto ${VETO} --sigsel genhtt 2.5 3.5 --densel region hadmu_base --sample ${SAMP}
  done
  python plot_lepid_roc.py --hists ${INFILE} --histname lep_kin --varname miso --varlabel miniIso --varcuts region hadel --title "Electron ID, ${TITLE} (Veto)" --label hadel_miso_0p05_${VETO} --tag LepID_Feb05 --veto ${VETO} --sigsel genhtt 1.5 2.5 --point 0.05
  python plot_lepid_roc.py --hists ${INFILE} --histname lep_kin --varname miso --varlabel miniIso --varcuts region hadmu --title "Muon ID, ${TITLE} (Veto)" --label hadmu_miso_0p05_${VETO} --tag LepID_Feb05 --veto ${VETO} --sigsel genhtt 2.5 3.5 --point 0.05
done
