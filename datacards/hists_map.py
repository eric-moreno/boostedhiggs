from coffea import hist

hists_properties = {
    "jet0_dHhadhad": {
        "hist_name": "presel_hadseljet0",
        "var_label": r"Jet dak8",
        "var_cut_dir": 1,
    },
    "jet0_dHhadhadMD": {
        "hist_name": "presel_hadseljet0",
        "var_label": r"Jet dak8 MD",
        "var_cut_dir": 1,
    },
    "jet0_msd": {
        "hist_name": "presel_hadseljet0",
        "var_label": r"Jet m$_{SD}$",
        "var_cut_dir": 1,
    },
    "jet0_pt": {
        "hist_name": "presel_hadseljet0",
        "var_label": r"Jet $p_T$",
        "var_cut_dir": 1,
    },
}

process_latex = {
    #'tt': r'$t\bar{t}$',
    "tt-dilep": r"$t\bar{t}~(dileptonic)$",
    "tt-semilep": r"$t\bar{t}~(semileptonic)$",
    "tt-had": r"$t\bar{t}~(hadronic)$",
    "st": "Single-t",
    "vqq": "V(qq)",
    #'zqq': 'Z(qq)',
    #'wqq': 'W(qq)',
    "vv": r"VV",
    "qcd": "Multijet",
    "zll": r"Z($\ell\ell$)",
    "zee": r"Z($ee$)",
    "zem": r"Z($ee$/$\mu\mu$)",
    "zmm": r"Z($\mu\mu$)",
    "ztt": r"Z($\tau\tau$)",
    #'zll-ht100to200': r'Z($\ell\ell$) $100<HT<200$',
    #'zll-ht1200to2500': r'Z($\ell\ell$) $1200<HT<2500$',
    #'zll-ht200to400': r'Z($\ell\ell$) $200<HT<400$',
    #'zll-ht2500toinf': r'Z($\ell\ell$) $2500<HT$',
    #'zll-ht400to600': r'Z($\ell\ell$) $400<HT<600$',
    #'zll-ht600to800': r'Z($\ell\ell$) $600<HT<800$',
    #'zll-ht800to1200': r'Z($\ell\ell$) $800<HT<1200$',
    "wjets": r"W($\ell\nu$)",
    "vv": "Diboson",
    "h125": r"H(125)",
    "phi10": r"$\phi(\tau\tau)$, $m=10$",
    "phi20": r"$\phi(\tau\tau)$, $m=20$",
    "phi30": r"$\phi(\tau\tau)$, $m=30$",
    "phi40": r"$\phi(\tau\tau)$, $m=40$",
    "phi50": r"$\phi(\tau\tau)$, $m=50$",
    "phi75": r"$\phi(\tau\tau)$, $m=75$",
    "phi100": r"$\phi(\tau\tau)$, $m=100$",
    "phi125": r"$\phi(\tau\tau)$, $m=125$",
    "phi150": r"$\phi(\tau\tau)$, $m=150$",
    "phi200": r"$\phi(\tau\tau)$, $m=200$",
    "phi250": r"$\phi(\tau\tau)$, $m=250$",
    "phi300": r"$\phi(\tau\tau)$, $m=300$",
    "phi10_unmatch": r"$\phi(\tau\tau)$, $m=10$ (unmatched)",
    "phi20_unmatch": r"$\phi(\tau\tau)$, $m=20$ (unmatched)",
    "phi30_unmatch": r"$\phi(\tau\tau)$, $m=30$ (unmatched)",
    "phi40_unmatch": r"$\phi(\tau\tau)$, $m=40$ (unmatched)",
    "phi50_unmatch": r"$\phi(\tau\tau)$, $m=50$ (unmatched)",
    "phi75_unmatch": r"$\phi(\tau\tau)$, $m=75$ (unmatched)",
    "phi100_unmatch": r"$\phi(\tau\tau)$, $m=100$ (unmatched)",
    "phi125_unmatch": r"$\phi(\tau\tau)$, $m=125$ (unmatched)",
    "phi150_unmatch": r"$\phi(\tau\tau)$, $m=150$ (unmatched)",
    "phi200_unmatch": r"$\phi(\tau\tau)$, $m=200$ (unmatched)",
    "phi250_unmatch": r"$\phi(\tau\tau)$, $m=250$ (unmatched)",
    "phi300_unmatch": r"$\phi(\tau\tau)$, $m=300$ (unmatched)",
    "ggF-Htt": "ggF-Htt",
    "VBF-Htt": "VBF-Htt",
    "Wm-Htt": "Wm-Htt",
    "Wp-Htt": "Wp-Htt",
    "ZH-Htt": "ZH-Htt",
    "ggZll-Htt": "ggZll-Htt",
    "ggZvv-Htt": "ggZvv-Htt",
    "ggZqq-Htt": "ggZqq-Htt",
    "tt-Htt": "tt-Htt",
    "data": "Data",
    "Stat. Unc.": "Stat. Unc.",
}

color_map = {
    "tt-dilep": "salmon",
    "tt-semilep": "firebrick",
    "tt-had": "rosybrown",
    "st": "sienna",
    "vqq": "darkorange",
    "qcd": "royalblue",
    "zll": "darkorchid",
    "zee": "mediumorchid",
    "zem": "mediumorchid",
    "zmm": "purple",
    "ztt": "indigo",
    "wjets": "forestgreen",
    "vv": "lightpink",
    "h125": "deepskyblue",
    "phi10": "snow",
    "phi20": "snow",
    "phi30": "snow",
    "phi40": "snow",
    "phi50": "snow",
    "phi75": "snow",
    "phi100": "snow",
    "phi125": "snow",
    "phi150": "snow",
    "phi200": "snow",
    "phi250": "snow",
    "phi300": "snow",
    "phi10_unmatch": "snow",
    "phi20_unmatch": "snow",
    "phi30_unmatch": "snow",
    "phi40_unmatch": "snow",
    "phi50_unmatch": "snow",
    "phi75_unmatch": "snow",
    "phi100_unmatch": "snow",
    "phi125_unmatch": "snow",
    "phi150_unmatch": "snow",
    "phi200_unmatch": "snow",
    "phi250_unmatch": "snow",
    "phi300_unmatch": "snow",
}
sig_color_map = {
    "phi10": "g",
    "phi20": "c",
    "phi30": "y",
    "phi40": "m",
    "phi50": "k",
    "phi75": "gray",
    "phi100": "b",
    "phi125": "r",
    "phi150": "tan",
    "phi200": "chartreuse",
    "phi250": "magenta",
    "phi300": "hotpink",
    "phi10_unmatch": "silver",
    "phi20_unmatch": "silver",
    "phi30_unmatch": "silver",
    "phi40_unmatch": "silver",
    "phi50_unmatch": "silver",
    "phi75_unmatch": "silver",
    "phi100_unmatch": "silver",
    "phi125_unmatch": "silver",
    "phi150_unmatch": "silver",
    "phi200_unmatch": "silver",
    "phi250_unmatch": "silver",
    "phi300_unmatch": "silver",
}

import re

# nosig = re.compile("(?!phi125)(?!h125)(?!data)")
nosig = re.compile(
    "(?!phi10)(?!phi20)(?!phi30)(?!phi40)(?!phi50)(?!phi75)(?!phi100)(?!phi125)(?!phi150)(?!phi200)(?!phi250)(?!phi300)(?!data)"
)
noqcd = re.compile(
    "(?!qcd)(?!phi10)(?!phi20)(?!phi30)(?!phi40)(?!phi50)(?!phi75)(?!phi100)(?!phi125)(?!phi150)(?!phi200)(?!phi250)(?!phi300)(?!data)"
)
# nosig = re.compile("(?!ggF-Htt)(?!VBF-Htt)(?!Wm-Htt)(?!Wp-Htt)(?!ZH-Htt)(?!ggZll-Htt)(?!ggZvv-Htt)(?!ggZqq-Htt)(?!tt-Htt)(?!h125)(?!data)")

# nobkg = re.compile("(?!qcd)(?!tt)(?!st)(?!zqq)(?!wjets)(?!vv)(?!wqq)(?!zll)")
# nobkg = re.compile("(?!qcd)(?!tt)(?!st)(?!vqq)(?!wjets)(?!vv)(?!zll)")
# nobkg = re.compile("(?!qcd)(?!tt-dilep)(?!tt-semilep)(?!tt-had)(?!st)(?!zqq)(?!wjets)(?!vv)(?!wqq)(?!zll)")

# nobkg = re.compile("(?!qcd)(?!tt-dilep)(?!tt-semilep)(?!tt-had)(?!st)(?!vqq)(?!wjets)(?!vv)(?!zll)(?!zee)(?!zmm)(?!ztt)(?!data)")
nobkg = re.compile(
    "(?!qcd)(?!tt-dilep)(?!tt-semilep)(?!tt-had)(?!st)(?!vqq)(?!wjets)(?!vv)(?!zll)(?!zem)(?!zee)(?!zmm)(?!ztt)(?!zll)(?!data)(?!h125)"
)

# nobkg = re.compile("(?!qcd)(?!tt-dilep)(?!tt-semilep)(?!tt-had)(?!st)(?!vqq)(?!wjets)(?!vv)(?!zll-ht100to200)(?!zll-ht1200to2500)(?!zll-ht200to400)(?!zll-ht2500toinf)(?!zll-ht400to600)(?!zll-ht600to800)(?!zll-ht800to1200)(?!data)")
