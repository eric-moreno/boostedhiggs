import sys, os
import json
import pickle
import numpy as np
import math
import collections
import logging

import rhalphalib as rl

import matplotlib.pyplot as plt
import mplhep as hep

import processmap
import samplelists
from coffea import hist
from utils import getHist
from cards import Cards

import click

rl.ParametericSample.PreferRooParametricHist = True
plt.style.use(hep.styles.ROOT)

log_level_dict = {
    "DEBUG": 10,
    "INFO": 20,
    "WARNING": 30,
    "ERROR": 40,
    "CRITICAL": 50,
}


def getQCDShape(c, hists):
    """
    gets data-MC (and dn/up variations from poisson unc.)
    for each "sig", "top_cr", "wlnu_cr", "qcd_cr" region
    """
    hists_qcd = collections.defaultdict(dict)
    for region in c.regions:
        hists_qcd[region] = collections.defaultdict(dict)
        for invregion in c.invregions:
            for category in c.nnregions:
                data_hist = c.get_qcd(
                    hists[region][invregion],
                    category,
                )
                hists_qcd[region][invregion][category] = data_hist
    return hists_qcd


def getRatioFail(c, hists_qcd):
    """
    get ratio to NN fail region
    pass/fail or loosepass/fail or fail/fail

    for each "top_cr", "wlnu_cr", "qcd_cr" region
    """
    ratio_F = collections.defaultdict(dict)
    num_to_denom = {"loosepass": "fail", "pass": "fail", "fail": "fail"}

    def get_ratio(numerator, denominator):
        return np.nan_to_num(
            np.ones_like(numerator) * (np.sum(numerator) / np.sum(denominator))
        )

    regions_to_loop = {
        "noqcdcr_fail": ["noqcdcr", "fail"],
        "qcdcr_nom": ["qcd_cr", "nom"],
        "qcdcr_fail": ["qcd_cr", "fail"],
    }

    for pred_region in ["sig", "top_cr", "wlnu_cr"]:
        for rtag, ritem in regions_to_loop.items():
            ratio_F[pred_region][rtag] = collections.defaultdict(dict)
            region, rinv = ritem
            region = pred_region if region == "noqcdcr" else region

            for num, denom in num_to_denom.items():
                ratio_F[pred_region][rtag][num] = collections.defaultdict(dict)
                for iv, var in enumerate(c.variations):
                    ratio_F[pred_region][rtag][num][var] = get_ratio(
                        hists_qcd[region][rinv][num][iv],
                        hists_qcd[region][rinv][denom][iv],
                    )

    return ratio_F


def getQCDRatio(c, hists_qcd):
    """
    Get ratios:
    sig_to_qcd: sig_fail(p,lp,f)/qcd_cr_fail
    or
    nom_to_fail: qcd_cr_nom(p,lp,f)/qcd_cr_fail
    """
    logger = logging.getLogger("get-qcd-ratio")

    qcd_ratio = collections.defaultdict(dict)
    for region in c.regions:
        if region == "qcd_cr":
            continue

        qcd_ratio[region] = collections.defaultdict(dict)
        for qcdregion in ["sig_to_qcd", "nom_to_fail"]:
            qcd_ratio[region][qcdregion] = collections.defaultdict(dict)
            # denominator always from the QCD CR (fail)
            denom = hists_qcd["qcd_cr"]["fail"]["fail"][0]
            denom_dn = hists_qcd["qcd_cr"]["fail"]["fail"][1]
            denom_up = hists_qcd["qcd_cr"]["fail"]["fail"][2]

            # for each pass/loosepass/fail
            for category in c.nnregions:
                qcd_ratio[region][qcdregion][category] = collections.defaultdict(dict)

                num = {
                    "sig_to_qcd": {
                        "sig": hists_qcd["sig"]["fail"][category],
                        "top_cr": hists_qcd["top_cr"]["fail"][category],
                        "wlnu_cr": hists_qcd["wlnu_cr"]["fail"][category],
                    },
                    "nom_to_fail": {
                        "sig": hists_qcd["qcd_cr"]["nom"][category],
                        "top_cr": hists_qcd["qcd_cr"]["nom"][category],
                        "wlnu_cr": hists_qcd["qcd_cr"]["nom"][category],
                    },
                }[qcdregion][region]

                norm = {
                    "sig_to_qcd": hists_qcd["qcd_cr"]["nom"]["fail"],  # H
                    "nom_to_fail": hists_qcd[region]["fail"]["fail"],  # A
                }[qcdregion]

                for iv, var in enumerate(c.variations):
                    qcd_ratio[region][qcdregion][category][var] = []

                    # take nom/dn/up variations of that template
                    tvar = num[iv]
                    # this will be num/denom
                    ratio = []

                    for ix in range(len(denom)):
                        if denom_dn[ix] > 0.0:
                            # the denominator is always the nominal variation
                            ratio.append(tvar[ix] / denom[ix])
                        else:
                            if ix == 0:
                                ratio.append(1.0)
                            else:
                                ratio.append(ratio[-1])
                    ratio = np.array(ratio)

                    # normalize to norm
                    qcd_ratio[region][qcdregion][category][var] = np.sum(
                        ratio * norm[iv]
                    ) / (np.sum(norm[iv]) if np.sum(norm[iv]) > 0.0 else 1.0)

                    # this should be one single number
                    logger.info(
                        f"qcd_ratio in {region}{qcdregion}{category}{var} is {qcd_ratio[region][qcdregion][category][var]}"
                    )
    return qcd_ratio


def createCards(hist_dict, cat, year, odir, label, unblind=False, no_syst=False):
    """
    Create cards.

    :param: hist_dict
    :type: dict
    :param cat: category (hadhad, hadel, hadmu)
    :type str
    :param year
    :type int
    """

    logger = logging.getLogger("create-cards")

    # create cards class (holds default settings)
    c = Cards(cat, year, no_syst)

    # extract histograms
    logger.info(f"Extracting histograms for each region {hist_dict.keys()}")
    hists = collections.defaultdict(dict)
    sample_cat = hist.Cat("sample", "sample")
    for key in c.regions:
        for invregion, metslices in c.invregions.items():
            metslice = metslices[key]
            hists[key][invregion] = (
                hist_dict[key][invregion]
                .group("process", sample_cat, c.sample_groups)
                .integrate("met_pt", metslice["slice"], metslice["overflow"])
            )

    # data driven QCD
    logger.info(f"Extracting data-driven QCD histograms for each region")
    # get (data-MC)
    hists_qcd = getQCDShape(c, hists)
    # get ratio to fail region (["noqcdcr_fail","qcdcr_nom","qcdcr_fail"])
    ratio_F = getRatioFail(c, hists_qcd)
    # get QCD ratio (["sig_to_qcd", "nom_to_fail"])
    qcd_ratio = getQCDRatio(c, hists_qcd)

    def build_qcdpred(region, islephad):
        """
        Build QCD prediction

        :param region:
        :type region: str
        :param islephad:
        :type islephad: bool
        :param var: Variation, Can be "nom","dn","up"
        :type var: str

        return: qcd predictions for that variation of the QCD estimate (nom/up/dn)
        """

        """
                 sig/top_cr/wlnu_cr             qcd_cr
             fail loosepass pass      fail loosepass pass
        fail  A      C       D         G       I      J
        nom   B      E       F         H       K      L

        Method 1:
        B = H(shape) * A/G*(norm to H) * ratio_to_fail_AA
        E = H(shape) * C/G*(norm to H) * ratio_to_fail_AC
        F = H(shape) * D/G*(norm to H) * ratio_to_fail_AD

        Method 2:
        B = A(shape) * H/G*(norm to H) * ratio_to_fail_HH
        E = A(shape) * K/G*(norm to H) * ratio_to_fail_HK
        F = A(shape) * L/G*(norm to H) * ratio_to_fail_HL

        noqcdcr_fail = ratio_to_fail_A
        qcdcr_nom = ratio_to_fail_H
        """
        # qcdpred: [shape,ratio,ratio_to_fail]
        qcdpred = {
            1: {
                "sig": [
                    hists_qcd["qcd_cr"]["nom"]["fail"],  # H
                    qcd_ratio["sig"]["sig_to_qcd"][region],  # A/G, C/G, D/G
                    ratio_F["sig"]["noqcdcr_fail"][region],  # ratio_to_fail
                ],
                "top_cr": [
                    hists_qcd["qcd_cr"]["nom"]["fail"],
                    qcd_ratio["top_cr"]["sig_to_qcd"][region],
                    ratio_F["top_cr"]["noqcdcr_fail"][region],
                ],
                "wlnu_cr": [
                    hists_qcd["qcd_cr"]["nom"]["fail"],
                    qcd_ratio["wlnu_cr"]["sig_to_qcd"][region],
                    ratio_F["wlnu_cr"]["noqcdcr_fail"][region],
                ],
            },
            2: {
                "sig": [
                    hists_qcd["sig"]["fail"]["fail"],  # A
                    qcd_ratio["sig"]["nom_to_fail"][region],  # H/G,K/G,L/G
                    ratio_F["sig"]["qcdcr_fail"][region],  # ratio_to_fail
                ],
                "top_cr": [
                    hists_qcd["top_cr"]["fail"]["fail"],  # A
                    qcd_ratio["top_cr"]["nom_to_fail"][region],
                    ratio_F["top_cr"]["qcdcr_fail"][region],
                ],
                "wlnu_cr": [
                    hists_qcd["wlnu_cr"]["fail"]["fail"],
                    qcd_ratio["wlnu_cr"]["nom_to_fail"][region],
                    ratio_F["wlnu_cr"]["qcdcr_fail"][region],
                ],
            },
        }

        qcd_index = {"nom": 0, "dn": 1, "up": 2}
        qcd_pred = collections.defaultdict(dict)
        for key in qcdpred[1].keys():
            for var, qid in qcd_index.items():
                # TODO: add prediction from 2
                shape_, ratio_, ratio_F_ = qcdpred[1][key]
                qcd_pred[key][var] = shape_[qid] * ratio_[var] * ratio_F_[var]

        return qcd_pred

    # build regions
    for region in c.nnregions:
        # when to integrate prediction over a single mass bin
        if c.islephad:
            singlebinCR = True
            singlebinFail = False
        else:
            singlebinCR = False
            singlebinFail = True
        singlebin = singlebinFail and region == "fail"
        singlebinCR = singlebinCR or singlebin

        logger.info(f"Building {region}-signal region")

        # build qcd prediction
        qcdpred = build_qcdpred(region, c.islephad)

        # add signal region for that pt bin
        c.build_channel(
            "",
            region,
            hists["sig"]["nom"],
            qcdpred["sig"],
            unblind,
            singlebin=singlebin,
            debug=False,
        )

        # add top control region
        logging.info(f"Building {region}-topCR region")
        c.build_channel(
            "topCR",
            region,
            hists["top_cr"]["nom"],
            qcdpred["top_cr"],
            singlebin=singlebinCR,
        )

        # add wlnu control region
        logging.info(f"Building {region}-wlnuCR region")
        c.build_channel(
            "wlnuCR",
            region,
            hists["wlnu_cr"]["nom"],
            qcdpred["wlnu_cr"],
            singlebin=singlebinCR,
        )

    for region in ["SR", "topCR", "wlnuCR"]:
        region_str = "" if region == "SR" else region
        str_fail = f"{region_str}fail{cat}{year}"
        str_loose = f"{region_str}loosepass{cat}{year}"
        str_pass = f"{region_str}pass{cat}{year}"
        c.set_expectation(str_fail, str_loose, str_pass)

    # save model
    with open(os.path.join(f"{odir}/{label}", f"{cat}Model.pkl"), "wb") as fout:
        pickle.dump(c.model, fout, protocol=2)


def loadHists(histogram, year, sigscale, sigxs, var="met_nn_kin"):
    logger = logging.getLogger("load-hists")

    lumi_dict = {
        "2016": 36.5,
        "2017": 41.5,
        "2018": 59.8,
    }
    lumi = lumi_dict[year]

    logger.info(f"Loading histograms from {year}")
    mcSamples = samplelists.getSamplesMC(year)
    dataSamples = samplelists.getSamplesData(year)

    regions_had = [
        "signal_met",
        "cr_mu",
        # "cr_b_met",
        "cr_b_mu_iso",
        "cr_mu_iso",
        "cr_dphi_inv",
        # "cr_b_met_dphi_inv",
        "cr_mu_dphi_inv",
        "cr_b_mu_iso_dphi_inv",
        "cr_mu_iso_dphi_inv",
    ]
    regions_lep = [
        "signal",
        "cr_dphi_inv",
        "cr_b",
        "cr_w",
        "cr_qcd",
        "cr_b_dphi_inv",
        "cr_w_dphi_inv",
        "cr_qcd_dphi_inv",
    ]
    regions = [f"hadhad_{r}" for r in regions_had]
    regions.extend([f"hadel_{r}" for r in regions_lep])
    regions.extend([f"hadmu_{r}" for r in regions_lep])

    includeData = False

    with open("fileset/xsecs.json", "r") as f:
        xs = json.load(f)
    scale1 = {k: sigxs / xs[k] for k in xs if "Spin0" in k and sigxs > 0.0}

    h_dict = {reg: None for reg in regions}
    for hs in mcSamples + dataSamples:
        try:
            with open(f"{histogram}{hs}.hist", "rb") as f:
                hists_unmapped = pickle.load(f)
        except:
            logger.error(f"Unable to open file {histogram}{hs}.hist")

        try:
            hist_unmapped = hists_unmapped["met_nn_kin"]
        except:
            logger.error(f"No {var} histogram in {histogram}{hs}.hist")

        for reg in regions:
            logger.debug(f"Mapping histograms from histogram {hs} and region {reg}")
            # map to hists
            toadd = hist_unmapped
            toadd.scale(
                {
                    p: scale1[p.name] if p.name in scale1 else 1.0
                    for p in toadd.identifiers("dataset")
                },
                axis="dataset",
            )
            if h_dict[reg] is not None:
                h_dict[reg] = h_dict[reg].add(
                    processmap.apply(toadd).integrate("region", reg)
                )
            else:
                h_dict[reg] = processmap.apply(toadd).integrate("region", reg)
        logger.debug(f"Mapped histograms from histogram {hs} to each region")
    logger.info("Mapped all histograms")

    hist_dict = {}
    empty_dict = {"nn_disc": [], "met_pt": [], "systematic": []}
    for reg in regions:
        logging.debug(
            f"Loading histograms from region {reg}, and scaling MC by {lumi}*{sigscale}"
        )
        xhist = getHist(
            h_dict[reg], "massreg", lumi, empty_dict, "", includeData, sigscale, [1]
        )
        hist_dict[reg] = xhist

    return hist_dict


@click.command()
@click.option("--hist", required=True, help="hists pickle prefix")
@click.option(
    "--year",
    required=True,
    type=click.Choice(["2016", "2016APV", "2017", "2018"]),
    help="year",
)
@click.option("--odir", default="test", help="output directory")
@click.option("--tag", default="", help="tag")
@click.option("--label", default="", help="label")
@click.option("--sigxs", type=float, default=0.0, help="signal xs to scale to")
@click.option("--sigscale", type=float, default=1.0, help="scale signal by")
@click.option(
    "--cat",
    "categories",
    required=True,
    multiple=True,
    type=click.Choice(["hadhad", "hadel", "hadmu"]),
)
@click.option(
    "--loglevel",
    default="INFO",
    type=click.Choice(
        ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], case_sensitive=False
    ),
    help="Sets the verbosity of the logging",
)
def makeCards(hist, year, odir, tag, label, sigscale, sigxs, categories, loglevel):
    logging.basicConfig(
        level=log_level_dict[loglevel],
        format="[%(levelname)-10s:" "%(message)s",
    )

    logger = logging.getLogger("make-cards")

    odir = f"cards/{tag}/"
    os.system(f"mkdir -p {odir}/{label}")
    pwd = os.getcwd()

    # load histograms
    hist_dict = loadHists(
        hist,
        year,
        sigscale,
        sigxs,
        var="met_nn_kin",
    )

    # histograms to be used
    h_dict = {}
    cat = "hadhad"
    h_dict[cat] = {
        "sig": {
            "nom": hist_dict[f"{cat}_signal_met"],
            "fail": hist_dict[f"{cat}_cr_dphi_inv"],
        },
        "top_cr": {
            "nom": hist_dict[f"{cat}_cr_b_mu_iso"],
            "fail": hist_dict[f"{cat}_cr_b_mu_iso_dphi_inv"],
        },
        "wlnu_cr": {
            "nom": hist_dict[f"{cat}_cr_mu_iso"],
            "fail": hist_dict[f"{cat}_cr_mu_iso_dphi_inv"],
        },
        "qcd_cr": {
            "nom": hist_dict[f"{cat}_cr_mu"],
            "fail": hist_dict[f"{cat}_cr_mu_dphi_inv"],
        },
    }
    for cat in ["hadel", "hadmu"]:
        h_dict[cat] = {
            "sig": {
                "nom": hist_dict[f"{cat}_signal"],
                "fail": hist_dict[f"{cat}_cr_dphi_inv"],
            },
            "top_cr": {
                "nom": hist_dict[f"{cat}_cr_b"],
                "fail": hist_dict[f"{cat}_cr_b_dphi_inv"],
            },
            "wlnu_cr": {
                "nom": hist_dict[f"{cat}_cr_w"],
                "fail": hist_dict[f"{cat}_cr_w_dphi_inv"],
            },
            "qcd_cr": {
                "nom": hist_dict[f"{cat}_cr_qcd"],
                "fail": hist_dict[f"{cat}_cr_qcd_dphi_inv"],
            },
        }

    for cat in categories:
        createCards(
            h_dict[cat],
            cat=cat,
            year=year,
            odir=odir,
            label=label,
        )


if __name__ == "__main__":
    if not sys.warnoptions:
        import warnings

        warnings.simplefilter("ignore")
        os.environ["PYTHONWARNINGS"] = "ignore"
        sys.exit(makeCards())
