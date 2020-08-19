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
                 'wlnu': r'W($\ell\nu$)',
                 'vv': 'Diboson',
                 'h125': 'ggH(125)',
                 'Stat. Unc.': 'Stat. Unc.',
                 }

import re
nosig = re.compile("(?!h125)")
#nobkg = re.compile("(?!qcd)(?!tt)(?!st)(?!zqq)(?!wlnu)(?!vv)(?!wqq)(?!zll)")
#nobkg = re.compile("(?!qcd)(?!tt)(?!st)(?!vqq)(?!wlnu)(?!vv)(?!zll)")
#nobkg = re.compile("(?!qcd)(?!tt-dilep)(?!tt-semilep)(?!tt-had)(?!st)(?!zqq)(?!wlnu)(?!vv)(?!wqq)(?!zll)")
nobkg = re.compile("(?!qcd)(?!tt-dilep)(?!tt-semilep)(?!tt-had)(?!st)(?!vqq)(?!wlnu)(?!vv)(?!zll)")

