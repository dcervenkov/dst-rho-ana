#!/usr/bin/env python3
"""Module to extract, average, and pretty-print DSRhoYield results."""

from math import sqrt
import json
import os
import sys
import ROOT
from contextlib import contextmanager
from markdowntable import MarkdownTable


def print_line_bold(string):
    """Print a line in bold using terminal ANSI codes."""
    print(u"\u001b[1m", end="")
    print(string)
    print(u"\u001b[0m", end="")


def recover_yield(filename):
    """Extract yield results from a file.

    A total of 6 numbers are returned the fit value and uncertainty of CR,
    SCF, and BKG.
    """
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
    """Extract the MC truth based number of signal events from a directory.

    All ROOT files in the directory are read and summed.
    """
    chain = ROOT.TChain("h2000")
    chain.Add(dir + "/*.root")
    count = chain.Draw(
        "de", common_cuts() + "&&(evmcflag==1||evmcflag==7||evmcflag==8)", "goff"
    )
    return count


def common_cuts():
    """Get common cuts string (should be same as in the C++ programs)."""
    return (
        "vrusable == 1 &&"
        "vtusable==1&&"
        "((vrchi2/vrndf) < 50 | |vrntrk == 1) &&"
        "((vtchi2/vtndf)<50||vtntrk==1)&&"
        "((sqrt(vrerr6)<0.02&&vrntrk>1)||(sqrt(vrerr6)<0.05&&vrntrk==1))&&"
        "((sqrt(vterr6)<0.02&&vtntrk>1)||(sqrt(vterr6)<0.05&&vtntrk==1))&&"
        "csbdtg > -0.6 &&"
        "(de>-0.14&&de<0.068)&&"
        "(dt>-10&&dt<10)&&"
        "(thetab>0.65 && thetab<2.95)&&"
        "nocand<=2"
    )


def weighted_mean(number_list, weights):
    """Calculate a weighted mean of a list of values."""
    total = 0
    assert len(number_list) == len(weights)
    for i in range(len(number_list)):
        total += number_list[i] * weights[i]
    return total / sum(weights)


def error_of_sum(errors):
    """Calculate the uncertainty of a sum, given uncertainty of each term."""
    total = 0
    for error in errors:
        total += error ** 2
    return sqrt(total)


def calculate_fractions(sig_yield, scf_yield, bkg_yield):
    """Translate yield numbers into fractions.

    Returns:
        A list of the following values:
            CR/(CR + SCF),
            CR/(CR + SCF + BKG),
            SCF/(CR + SCF + BKG),
            BKG/(CR + SCF + BKG)

    """
    cr_crscf = float(sig_yield) / (sig_yield + scf_yield)
    cr = float(sig_yield) / (sig_yield + scf_yield + bkg_yield)
    scf = float(scf_yield) / (sig_yield + scf_yield + bkg_yield)
    bkg = float(bkg_yield) / (sig_yield + scf_yield + bkg_yield)
    return (cr_crscf, cr, scf, bkg)


def process_MC(channels, streams):
    """Read, process, and pretty-print MC."""
    for channel in channels:
        print_line_bold(channel)
        header = [
            "Stream",
            "Yield",
            "Err",
            "Count",
            "Bias",
            "CR/CR+SCF",
            "CR/all",
            "SCF/all",
            "BKG/all",
        ]
        precision = [0, 0, 0, 0, 0, 3, 3, 3, 3]
        table = MarkdownTable(header)
        table.set_precision(precision)

        yields = []
        errors = []
        counts = []
        biases = []

        cr_crscf_fractions = []
        cr_fractions = []
        scf_fractions = []
        bkg_fractions = []

        for stream_no, stream in enumerate(streams):
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

            table.add_row(
                [
                    stream_no,
                    sig_yield,
                    sig_yield_error,
                    sig_count,
                    sig_yield - sig_count,
                    cr_crscf,
                    cr,
                    scf,
                    bkg,
                ]
            )
            yields.append(sig_yield)
            errors.append(sig_yield_error)
            counts.append(sig_count)
            biases.append(sig_yield - sig_count)

            cr_crscf_fractions.append(cr_crscf)
            cr_fractions.append(cr)
            scf_fractions.append(scf)
            bkg_fractions.append(bkg)

        table.add_row(
            [
                "Average",
                sum(yields) / len(yields),
                sum(errors) / len(errors),
                sum(counts) / len(counts),
                sum(biases) / len(biases),
                sum(cr_crscf_fractions) / len(cr_crscf_fractions),
                sum(cr_fractions) / len(cr_fractions),
                sum(scf_fractions) / len(scf_fractions),
                sum(bkg_fractions) / len(bkg_fractions),
            ]
        )

        table.print()
        print()

        data = {
            "cr_scf_f": round(sum(cr_crscf_fractions) / len(cr_crscf_fractions), 4),
            "cr_f": round(sum(cr_fractions) / len(cr_fractions), 4),
            "scf_f": round(sum(scf_fractions) / len(scf_fractions), 4),
        }

        save_to_json(data, os.path.join("results", channel + "_mc_fractions.json"))


def process_data(channels):
    """Read, process, and pretty-print data."""
    print_line_bold("Data")
    header = ["Channel", "Yield", "Err", "CR/CR+SCF", "CR/all", "SCF/all", "BKG/all"]
    precision = [None, 0, 0, 3, 3, 3, 3]

    table = MarkdownTable(header)
    table.set_precision(precision)

    yields = []
    errors = []

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

        (
            cr_crscf_fraction,
            cr_fraction,
            scf_fraction,
            bkg_fraction,
        ) = calculate_fractions(sig_yield, scf_yield, bkg_yield)

        table.add_row(
            [
                channel,
                sig_yield,
                sig_yield_error,
                cr_crscf_fraction,
                cr_fraction,
                scf_fraction,
                bkg_fraction,
            ]
        )

        yields.append(sig_yield)
        errors.append(sig_yield_error)

        cr_crscf_fractions.append(cr_crscf_fraction)
        cr_fractions.append(cr_fraction)
        scf_fractions.append(scf_fraction)
        bkg_fractions.append(bkg_fraction)

        data = {
            "cr_scf_f": round(cr_crscf_fraction, 4),
            "cr_f": round(cr_fraction, 4),
            "scf_f": round(scf_fraction, 4),
        }

        save_to_json(data, os.path.join("results", channel + "_data_fractions.json"))

    table.add_row(
        [
            "Tot/Avg",
            sum(yields),
            error_of_sum(errors),
            weighted_mean(cr_crscf_fractions, yields),
            weighted_mean(cr_fractions, yields),
            weighted_mean(scf_fractions, yields),
            weighted_mean(bkg_fractions, yields),
        ]
    )

    table.print()

    data = {
        "cr_scf_f": round(weighted_mean(cr_crscf_fractions, yields), 4),
        "cr_f": round(weighted_mean(cr_fractions, yields), 4),
        "scf_f": round(weighted_mean(scf_fractions, yields), 4)
    }

    save_to_json(data, os.path.join("results", "avg_data_fractions.json"))


def save_to_json(data, filename):
    """Save dictionary as a prettified JSON."""
    with open(filename, "w") as f:
        json.dump(data, f, indent=2, sort_keys=True)


@contextmanager
def stdout_redirected(to=os.devnull):
    """Redirect stdout.

    This can be used to silence unwanted print statements from polluting the
    output, or for redirecting the output into a log file.

    Args:
        to: An optional argument for specifying the redirection destination.
            If none is supplied, the default is os.devnull.

    Usage:
        import os

        with stdout_redirected(to=filename):
            print("from Python")
            os.system("echo non-Python applications are also supported")
    """
    fd = sys.stdout.fileno()

    def _redirect_stdout(to):
        sys.stdout.close()  # + implicit flush()
        os.dup2(to.fileno(), fd)  # fd writes to 'to' file
        sys.stdout = os.fdopen(fd, "w")  # Python writes to fd

    with os.fdopen(os.dup(fd), "w") as old_stdout:
        with open(to, "w") as file:
            _redirect_stdout(to=file)
        try:
            yield  # allow code to be run with the redirected stdout
        finally:
            _redirect_stdout(to=old_stdout)  # restore stdout.
            # buffering and flags such as
            # CLOEXEC may be different


with stdout_redirected():
    # Just to init RooFit to make it print its loading line now
    r = ROOT.RooRealVar()

channels = ["Kpi", "Kpipi0", "K3pi"]
streams = range(0, 6)

process_MC(channels, streams)
process_data(channels)
