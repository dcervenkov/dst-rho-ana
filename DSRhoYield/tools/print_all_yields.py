#!/usr/bin/env python

from __future__ import print_function
import re
from math import sqrt
import ROOT


def recover_yield(filename):
    file = ROOT.TFile(filename)
    sig_plus_cf = float(
        str(file.Get("all;1").GetLineWith("n_signal_plus_cross_model")).split()[1]
    )
    sig_plus_cf_error = float(
        str(file.Get("all;1").GetLineWith("n_signal_plus_cross_model")).split()[3]
    )
    f_sig_cf = float(
        str(
            file.Get("signal_plus_crossfeed;1").GetLineWith("signal_plus_cross_f")
        ).split()[1]
    )

    sig = int(round(sig_plus_cf * f_sig_cf))
    sig_error = int(round(sig_plus_cf_error * f_sig_cf))
    scf = int(round(sig_plus_cf * (1 - f_sig_cf)))
    scf_error = int(round(sig_plus_cf_error * (1 - f_sig_cf)))
    bkg = float(str(file.Get("all;1").GetLineWith("n_bkg")).split()[1])
    bkg_error = float(str(file.Get("all;1").GetLineWith("n_bkg")).split()[3])

    return sig, sig_error, scf, scf_error, bkg, bkg_error


def recover_count(dir):
    chain = ROOT.TChain("h2000")
    chain.Add(dir + "/*.root")
    count = chain.Draw(
        "de", common_cuts() + "&&(evmcflag==1||evmcflag==7||evmcflag==8)", "goff"
    )
    return count


def get_table_midline(header):
    midline = ""
    separators = [m.start() for m in re.finditer("\|", header)]
    for i in range(len(header)):
        if i in separators:
            midline += "|"
        else:
            midline += "-"
    return midline


def common_cuts():
    return "vrusable==1&&vtusable==1&&((vrchi2/vrndf)<50||vrntrk==1)&&((vtchi2/vtndf)<50||vtntrk==1)&&((sqrt(vrerr6)<0.02&&vrntrk>1)||(sqrt(vrerr6)<0.05&&vrntrk==1))&&((sqrt(vterr6)<0.02&&vtntrk>1)||(sqrt(vterr6)<0.05&&vtntrk==1))&&csbdtg>-0.6&&(de>-0.14&&de<0.068)&&(dt>-10&&dt<10)&&thetab>0.5"


def weighted_mean(number_list, weights):
    total = 0
    assert len(number_list) == len(weights)
    for i in range(len(number_list)):
        total += number_list[i] * weights[i]
    return total / sum(weights)


def error_of_sum(errors):
    total = 0
    for error in errors:
        total += error ** 2
    return sqrt(total)


def calculate_fractions(sig_yield, scf_yield, bkg_yield):
    cr_crscf = sig_yield / (sig_yield + scf_yield)
    cr = sig_yield / (sig_yield + scf_yield + bkg_yield)
    scf = scf_yield / (sig_yield + scf_yield + bkg_yield)
    bkg = bkg_yield / (sig_yield + scf_yield + bkg_yield)
    return (cr_crscf, cr, scf, bkg)


def process_MC(channels, streams):
    print("Processing MC results\n")
    for channel in channels:
        print(channel)
        header = "{} +- {} | {} | {} | {} | {} | {} | {} | {}".format(
            "Yield",
            "Err",
            "Count",
            "Bias",
            "Purity",
            "CR/CR+SCF",
            "CR/all",
            "SCF/all",
            "BKG/all",
        )
        print(header)
        print(get_table_midline(header))
        yields = []
        errors = []
        counts = []
        biases = []
        purities = []

        cr_crscf_fractions = []
        cr_fractions = []
        scf_fractions = []
        bkg_fractions = []

        for stream in streams:
            (
                sig_yield,
                sig_yield_error,
                scf_yield,
                scf_yield_error,
                bkg_yield,
                bkg_yield_error,
            ) = recover_yield(
                "plots/" + channel + "_stream" + str(stream) + "/fit_results.root"
            )
            sig_count = recover_count(
                "../data/" + channel + "/realistic_mc/stream" + str(stream)
            )

            cr_crscf, cr, scf, bkg = calculate_fractions(
                sig_yield, scf_yield, bkg_yield
            )

            purity = sig_yield / (sig_yield + scf_yield + bkg_yield)
            print(
                "{:.0f} +- {:.0f} | {:.0f} | {:4.0f} | {:6.2f} | {:9.3f} | {:6.3f} | {:7.3f} | {:7.3f}".format(
                    sig_yield,
                    sig_yield_error,
                    sig_count,
                    sig_yield - sig_count,
                    purity,
                    cr_crscf,
                    cr,
                    scf,
                    bkg,
                )
            )
            yields.append(sig_yield)
            errors.append(sig_yield_error)
            counts.append(sig_count)
            biases.append(sig_yield - sig_count)
            purities.append(purity)

            cr_crscf_fractions.append(cr_crscf)
            cr_fractions.append(cr)
            scf_fractions.append(scf)
            bkg_fractions.append(bkg)

        print(get_table_midline(header))
        print(
            "{:.0f} +- {:.0f} | {:.0f} | {:4.0f} | {:6.2f} | {:9.3f} | {:6.3f} | {:7.3f} | {:7.3f}".format(
                sum(yields) / len(yields),
                sum(errors) / len(errors),
                sum(counts) / len(counts),
                sum(biases) / len(biases),
                sum(purities) / len(purities),
                sum(cr_crscf_fractions) / len(cr_crscf_fractions),
                sum(cr_fractions) / len(cr_fractions),
                sum(scf_fractions) / len(scf_fractions),
                sum(bkg_fractions) / len(bkg_fractions),
            )
        )
        print()


def process_data(channels):
    print("Processing data results\n")
    header = "{} | {} +- {} | {} | {} | {} | {} | {}".format(
        "Channel", "Yield", "Err", "Purity", "CR/CR+SCF", "CR/all", "SCF/all", "BKG/all"
    )
    print(header)
    print(get_table_midline(header))
    yields = []
    errors = []
    purities = []

    cr_crscf_fractions = []
    cr_fractions = []
    scf_fractions = []
    bkg_fractions = []

    for channel in channels:
        (
            sig_yield,
            sig_yield_error,
            scf_yield,
            scf_yield_error,
            bkg_yield,
            bkg_yield_error,
        ) = recover_yield("plots/" + channel + "/fit_results.root")

        purity = sig_yield / (sig_yield + scf_yield + bkg_yield)
        cr_crscf, cr, scf, bkg = calculate_fractions(sig_yield, scf_yield, bkg_yield)

        print(
            "{:7} | {:.0f} +- {:.0f} | {:6.2f} | {:9.3f} | {:6.3f} | {:7.3f} | {:7.3f}".format(
                channel, sig_yield, sig_yield_error, purity, cr_crscf, cr, scf, bkg
            )
        )
        yields.append(sig_yield)
        errors.append(sig_yield_error)
        purities.append(purity)

        cr_crscf_fractions.append(cr_crscf)
        cr_fractions.append(cr)
        scf_fractions.append(scf)
        bkg_fractions.append(bkg)

    print(get_table_midline(header))
    print(
        "{:7} | {:.0f} +- {:.0f} | {:6.2f} | {:9.3f} | {:6.3f} | {:7.3f} | {:7.3f}".format(
            "Total",
            sum(yields),
            error_of_sum(errors),
            weighted_mean(purities, yields),
            weighted_mean(cr_crscf_fractions, yields),
            weighted_mean(cr_fractions, yields),
            weighted_mean(scf_fractions, yields),
            weighted_mean(bkg_fractions, yields),
        )
    )


# Just to init RooFit to make it print its loading line now
r = ROOT.RooRealVar()

channels = ["Kpi", "Kpipi0", "K3pi"]
streams = range(0, 6)

process_MC(channels, streams)
process_data(channels)
