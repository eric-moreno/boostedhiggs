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
                 #'zll-ht100to200': r'Z($\ell\ell$) $100<HT<200$',
                 #'zll-ht1200to2500': r'Z($\ell\ell$) $1200<HT<2500$',
                 #'zll-ht200to400': r'Z($\ell\ell$) $200<HT<400$',
                 #'zll-ht2500toinf': r'Z($\ell\ell$) $2500<HT$',
                 #'zll-ht400to600': r'Z($\ell\ell$) $400<HT<600$',
                 #'zll-ht600to800': r'Z($\ell\ell$) $600<HT<800$',
                 #'zll-ht800to1200': r'Z($\ell\ell$) $800<HT<1200$',
                 'wlnu': r'W($\ell\nu$)',
                 'vv': 'Diboson',
                 'h125': 'ggH(125)',
                 'data': 'Data',
                 'Stat. Unc.': 'Stat. Unc.',
                 }

import re
nosig = re.compile("(?!h125)(?!data)")
#nobkg = re.compile("(?!qcd)(?!tt)(?!st)(?!zqq)(?!wlnu)(?!vv)(?!wqq)(?!zll)")
#nobkg = re.compile("(?!qcd)(?!tt)(?!st)(?!vqq)(?!wlnu)(?!vv)(?!zll)")
#nobkg = re.compile("(?!qcd)(?!tt-dilep)(?!tt-semilep)(?!tt-had)(?!st)(?!zqq)(?!wlnu)(?!vv)(?!wqq)(?!zll)")

nobkg = re.compile("(?!qcd)(?!tt-dilep)(?!tt-semilep)(?!tt-had)(?!st)(?!vqq)(?!wlnu)(?!vv)(?!zll)(?!data)") 
#nobkg = re.compile("(?!qcd)(?!tt-dilep)(?!tt-semilep)(?!tt-had)(?!st)(?!vqq)(?!wlnu)(?!vv)(?!zll-ht100to200)(?!zll-ht1200to2500)(?!zll-ht200to400)(?!zll-ht2500toinf)(?!zll-ht400to600)(?!zll-ht600to800)(?!zll-ht800to1200)(?!data)") 
