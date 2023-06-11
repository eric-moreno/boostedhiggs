from hists_map import *
from coffea import hist
import numpy as np
import logging


def intRegion(
    inhist,
    theregion,
    nnCut,
    nnCut_loose,
    systematic="nominal",
    samplelist=None,
    mslice=None,
    mrebin=None,
    debug=False,
):
    """
    Return result of integrating over nnCuts.
    If there is a mass slice (mslice) also integrate over mass.
    """
    # slice on nn
    if theregion == "pass":
        theslice = slice(nnCut, None)
        overflow_str = "over"
    elif theregion == "loosepass":
        theslice = slice(nnCut_loose, nnCut)
        overflow_str = "none"
    elif theregion == "fail":
        theslice = slice(None, nnCut_loose)
        overflow_str = "under"
    else:
        print("Unknown region", theregion)
        return
    the_int = inhist.integrate("nn_disc", theslice, overflow_str).integrate(
        "systematic", systematic
    )
    if debug:
        print("nn slice ", theslice, "nn ", nnCut, " l ", nnCut_loose)
        print("overflow ", overflow_str, the_int.values())

    logging.debug(
        f"Values: {the_int.values()} with nncuts ({nnCut_loose} and {nnCut}) and overflow {overflow_str}"
    )

    # rebin in mass
    if mrebin is not None:
        the_int = the_int.rebin("massreg", hist.Bin("massreg", "massreg", mrebin))

    # slice on mass
    if mslice is not None:
        if mslice.start is not None and mslice.stop is not None:
            overflow_str = "none"
        elif mslice.start is None and mslice.stop is not None:
            overflow_str = "under"
        elif mslice.start is not None and mslice.stop is None:
            overflow_str = "over"
        else:
            overflow_str = "allnan"
        the_int = the_int.integrate("massreg", mslice, overflow_str)
        if debug:
            print("mass ", overflow_str, the_int.values())
        logging.debug(
            f"Values: {the_int.values()} with mass slices and overflow {overflow_str}"
        )

    # slice or sum samples
    if samplelist is not None:
        the_int = the_int.integrate("sample", samplelist).values(sumw2=True)[()]
    else:
        the_int = the_int.sum("sample").values(sumw2=True)
        if () in the_int:
            the_int = the_int[()]
        else:
            the_int = np.zeros(len(mttbins) - 1, dtype=np.float32)

    if debug:
        print("debug", the_int)

    logging.debug(f"After integrating: {the_int}")
    return the_int


def getQCDFromData(
    inhist,
    region,
    nnCut,
    nnCut_loose,
    default=1.0,
    systematic="nominal",
    mslice=None,
    unc_scale=0.0,
    lowmassbin=2,
    highmassbin=-1,
):
    """
    Get Data - other(MC) = QCD.
    :param default: sets event content for bins with 0 events
    :type default: float
    :param unc_scale: can force the errors to be at least some percent (potentially useful for bins with a lot of events)
    :type unc_scale: float
    """

    def clipZeros(data):
        data = np.clip(data, default, None)
        data[data == default] = 0.0
        return data

    qcd_data_full = intRegion(
        inhist["data_obs"],
        region,
        nnCut,
        nnCut_loose,
        systematic=systematic,
        mslice=mslice,
    )[0][lowmassbin:highmassbin]

    other_mc = intRegion(
        inhist,
        region,
        nnCut,
        nnCut_loose,
        samplelist=[
            s.name
            for s in inhist.identifiers("sample")
            if s.name != "data_obs" and s.name != "multijet"
        ],
        systematic=systematic,
        mslice=mslice,
    )

    qcd_data = clipZeros(qcd_data_full - other_mc[0][lowmassbin:highmassbin])
    # logging.debug(f"    Data: %s"%qcd_data_full)
    # logging.debug(f"    Other MC: %s"%other_mc[0][lowmassbin:highmassbin])
    logging.debug(f"    Data - MC: %s" % qcd_data)

    # QCD from data (data - other mc)
    qcd_temp = qcd_data

    # get variations ( from the poisson uncertainty on the number of events in each bin )
    qcd_data_altdn = clipZeros(
        qcd_data_full - other_mc[0][lowmassbin:highmassbin] * (1.0 + unc_scale)
    )
    qcd_data_altup = clipZeros(
        qcd_data_full - other_mc[0][lowmassbin:highmassbin] * (1.0 - unc_scale)
    )

    # QCD from data (data - other mc)
    qcd_data_w2 = np.clip(
        qcd_data_full - other_mc[1][lowmassbin:highmassbin], default, None
    )
    qcd_data_w2[qcd_data_w2 == default] = 0.0

    qcd_data_int = hist.poisson_interval(qcd_data, qcd_data_w2)
    qcd_temp_dn = np.minimum(
        np.array(
            [
                qcd_data_int[0][bi]
                if qcd_data_int[0][bi] >= 0.0 and not np.isnan(qcd_data_int[0][bi])
                else 0.0
                for bi in range(len(qcd_data))
            ]
        ),
        qcd_data_altdn,
    )
    qcd_temp_up = np.maximum(
        np.array(
            [
                qcd_data_int[1][bi]
                if qcd_data_int[1][bi] >= 0.0 and not np.isnan(qcd_data_int[1][bi])
                else 0.0
                for bi in range(len(qcd_data))
            ]
        ),
        qcd_data_altup,
    )
    return qcd_temp, qcd_temp_dn, qcd_temp_up


def getHist(
    h, var_name, lumifb, vars_cut, regionsel, blind, sigscale, rebin, debug=False
):
    """
    Get a histogram after integration and selections

    :param h: histogram
    :type: hist.Hist
    :param var_name: variable to keep
    :type str
    :param lumifb: luminosity
    :type float
    :param vars_cut
    :param regionsel
    :param blind
    :param sigscale
    :param rebin
    """
    if debug:
        print(h)
    exceptions = ["process", var_name]
    exceptions.extend([var for var in vars_cut])
    if regionsel != "":
        exceptions.append("region")
    overflow_sum = "allnan"
    x = h.sum(
        *[ax for ax in h.axes() if ax.name not in exceptions], overflow=overflow_sum
    )
    for reg in regionsel:
        x = x.integrate("region", reg)

    for var, val in vars_cut.items():
        if len(val) == 0:
            continue
        if var != var_name:
            if debug:
                print("integrating ", var, val[0], val[1])
            if val[0] is not None and val[1] is not None:
                overflow_str = "none"
            elif val[0] is None and val[1] is not None:
                overflow_str = "under"
            elif val[0] is not None and val[1] is None:
                overflow_str = "over"
            else:
                overflow_str = "allnan"
            x = x.integrate(var, slice(val[0], val[1]), overflow=overflow_str)

    if var_name in vars_cut.keys():
        if debug:
            print(
                "integrating ", var_name, vars_cut[var_name][0], vars_cut[var_name][1]
            )
        x = x.rebin(
            var_name,
            hist.Bin(
                var_name,
                var_name,
                [
                    xe
                    for xe in x.axis(var_name).edges()
                    if (xe >= vars_cut[var_name][0] and xe <= vars_cut[var_name][1])
                ],
            ),
        )

    xaxis = var_name
    for ih, hkey in enumerate(x.identifiers("process")):
        x.identifiers("process")[ih].label = process_latex[hkey.name]

    if len(rebin) == 1:
        if rebin[0] > 1:
            x = x.rebin(xaxis, int(rebin[0]))
    else:
        x = x.rebin(xaxis, hist.Bin(xaxis, var_name, rebin))

    x_nobkg = x[nobkg]
    x_nosig = x[nosig]
    if debug:
        print("lumifb ", lumifb)
        print("sigscale ", sigscale)
    x_nosig.scale({p: lumifb for p in x_nosig.identifiers("process")}, axis="process")
    x_nobkg.scale(
        {p: lumifb * float(sigscale) for p in x_nobkg.identifiers("process")},
        axis="process",
    )
    x_data = x["data"]

    return x_nobkg + x_nosig + x_data
