from coffea import hist

hists_properties = {'jet0_dHhadhad':{'hist_name': 'presel_hadseljet0',
                                   'var_label': r'Jet dak8',
                                   'var_cut_dir': 1,
                                 },
                    'jet0_dHhadhadMD':{'hist_name':'presel_hadseljet0',
                                     'var_label': r'Jet dak8 MD',
                                     'var_cut_dir': 1,
                                     },
                    'jet0_msd':{'hist_name':'presel_hadseljet0',
                                'var_label': r'Jet m$_{SD}$',
                                'var_cut_dir': 1,
                                },
                    'jet0_pt':{'hist_name':'presel_hadseljet0',
                               'var_label': r'Jet $p_T$',
                               'var_cut_dir': 1,
                               },
                    }

process_latex = {
                 #'tt': r'$t\bar{t}$',
                 'tt-dilep': r'$t\bar{t}~(dileptonic)$',
                 'tt-semilep': r'$t\bar{t}~(semileptonic)$',
                 'tt-had': r'$t\bar{t}~(hadronic)$',
                 'st': 'Single-t',
                 'vqq': 'V(qq)',
                 #'zqq': 'Z(qq)',
                 #'wqq': 'W(qq)',
                 'qcd': 'Multijet',
                 'zll': r'Z($\ell\ell$)',
                 'zee': r'Z($ee$)',
                 'zem': r'Z($ee$/$\mu\mu$)',
                 'zmm': r'Z($\mu\mu$)',
                 'ztt': r'Z($\tau\tau$)',
                 #'zll-ht100to200': r'Z($\ell\ell$) $100<HT<200$',
                 #'zll-ht1200to2500': r'Z($\ell\ell$) $1200<HT<2500$',
                 #'zll-ht200to400': r'Z($\ell\ell$) $200<HT<400$',
                 #'zll-ht2500toinf': r'Z($\ell\ell$) $2500<HT$',
                 #'zll-ht400to600': r'Z($\ell\ell$) $400<HT<600$',
                 #'zll-ht600to800': r'Z($\ell\ell$) $600<HT<800$',
                 #'zll-ht800to1200': r'Z($\ell\ell$) $800<HT<1200$',
                 'wjets': r'W($\ell\nu$)',
                 'vv': 'Diboson',
                 'h125': r'H(125)',
                 'phi10': r'$\phi(\tau\tau)$, $m=10$',
                 'phi30': r'$\phi(\tau\tau)$, $m=30$',
                 'phi50': r'$\phi(\tau\tau)$, $m=50$',
                 'phi125': r'$\phi(\tau\tau)$, $m=125$',
                 'phi300': r'$\phi(\tau\tau)$, $m=300$',
                 "ggF-Htt":"ggF-Htt",
                 "VBF-Htt":"VBF-Htt",
                 "Wm-Htt":"Wm-Htt",
                 "Wp-Htt":"Wp-Htt",
                 "ZH-Htt":"ZH-Htt",
                 "ggZll-Htt":"ggZll-Htt",
                 "ggZvv-Htt":"ggZvv-Htt",
                 "ggZqq-Htt":"ggZqq-Htt",
                 "tt-Htt":"tt-Htt",
                 'data': 'Data',
                 'Stat. Unc.': 'Stat. Unc.',
                 }

color_map = {
                 'tt-dilep':'salmon',
                 'tt-semilep':'firebrick',
                 'tt-had':'rosybrown',
                 'st':'sienna',
                 'vqq':'darkorange',
                 'qcd':'royalblue',
                 'zll':'darkorchid',
                 'zee':'mediumorchid',
                 'zem':'mediumorchid',
                 'zmm':'purple',
                 'ztt':'indigo',
                 'wjets':'forestgreen',
                 'vv':'lightpink',
                 'h125':'deepskyblue',
                 'phi10':'snow',
                 'phi30':'snow',
                 'phi50':'snow',
                 'phi125':'snow',
                 'phi300':'snow',
}
sig_color_map={
                 'phi10':'g',
                 'phi30':'y',
                 'phi50':'k',
                 'phi125':'r',
                 'phi300':'m',
}

import re
#nosig = re.compile("(?!phi125)(?!h125)(?!data)")
nosig = re.compile("(?!phi10)(?!phi30)(?!phi50)(?!phi125)(?!phi300)(?!data)")
noqcd = re.compile("(?!qcd)(?!phi10)(?!phi30)(?!phi50)(?!phi125)(?!phi300)(?!data)")
#nosig = re.compile("(?!ggF-Htt)(?!VBF-Htt)(?!Wm-Htt)(?!Wp-Htt)(?!ZH-Htt)(?!ggZll-Htt)(?!ggZvv-Htt)(?!ggZqq-Htt)(?!tt-Htt)(?!h125)(?!data)")

#nobkg = re.compile("(?!qcd)(?!tt)(?!st)(?!zqq)(?!wjets)(?!vv)(?!wqq)(?!zll)")
#nobkg = re.compile("(?!qcd)(?!tt)(?!st)(?!vqq)(?!wjets)(?!vv)(?!zll)")
#nobkg = re.compile("(?!qcd)(?!tt-dilep)(?!tt-semilep)(?!tt-had)(?!st)(?!zqq)(?!wjets)(?!vv)(?!wqq)(?!zll)")

#nobkg = re.compile("(?!qcd)(?!tt-dilep)(?!tt-semilep)(?!tt-had)(?!st)(?!vqq)(?!wjets)(?!vv)(?!zll)(?!zee)(?!zmm)(?!ztt)(?!data)") 
nobkg = re.compile("(?!qcd)(?!tt-dilep)(?!tt-semilep)(?!tt-had)(?!st)(?!vqq)(?!wjets)(?!vv)(?!zll)(?!zem)(?!zee)(?!zmm)(?!ztt)(?!zll)(?!data)(?!h125)") 

#nobkg = re.compile("(?!qcd)(?!tt-dilep)(?!tt-semilep)(?!tt-had)(?!st)(?!vqq)(?!wjets)(?!vv)(?!zll-ht100to200)(?!zll-ht1200to2500)(?!zll-ht200to400)(?!zll-ht2500toinf)(?!zll-ht400to600)(?!zll-ht600to800)(?!zll-ht800to1200)(?!data)") 
