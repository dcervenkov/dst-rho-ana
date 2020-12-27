#!/usr/bin/env python3
"""Module to extract, average, and pretty-print DSRhoYield results."""

import argparse
from contextlib import contextmanager
import json
from markdowntable import MarkdownTable
from math import sqrt
import numpy as np
import os
import pathlib
import ROOT
import sys


def decode_arguments():
    """Decode CLI arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d",
        "--directory",
        default="results",
        help="basename of directory where to save results",
    )
    parser.add_argument(
        "-r",
        "--randomize",
        type=int,
        default=0,
        help="number of randomized results to be generated",
    )
    parser.add_argument(
        "-c",
        "--correlations",
        action="store_true",
        help="account for the correlations when randomizing",
    )
    args = parser.parse_args()
    return args


def print_line_bold(string):
    """Print a line in bold using terminal ANSI codes."""
    print("\u001b[1m", end="")
    print(string)
    print("\u001b[0m", end="")


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


def recover_covariance_matrix(filename):
    """Extract the covariance matrix from a file.

    The covariance matrix for the final "all" fit without the width parameter
    is returned.
    """
    file = ROOT.TFile(filename)
    result = file.Get("fitresult_model_data_set")
    root_cov = result.covarianceMatrix()
    cov = np.zeros((2, 2))
    cov[0][0] = root_cov[0][0]
    cov[0][1] = root_cov[0][1]
    cov[1][0] = root_cov[1][0]
    cov[1][1] = root_cov[1][1]
    return cov


def recover_cr_scf_fraction(filename):
    """Extract the fraction of CR/SCF including uncertainty."""
    file = ROOT.TFile(filename)
    f_sig_cf = float(
        str(
            file.Get("signal_plus_crossfeed;1").GetLineWith("signal_plus_cross_f")
        ).split()[1]
    )
    f_sig_cf_error = float(
        str(
            file.Get("signal_plus_crossfeed;1").GetLineWith("signal_plus_cross_f")
        ).split()[3]
    )
    return f_sig_cf, f_sig_cf_error


def recover_count(dir, rbin):
    """Extract the MC truth based number of signal events from a directory.

    All ROOT files in the directory are read and summed.
    """
    chain = ROOT.TChain("h2000")
    chain.Add(dir + "/*.root")
    cuts = common_cuts()
    cuts += "&&(evmcflag==1||evmcflag==7||evmcflag==8)"
    if rbin is not None:
        cuts += f"&&{get_rbin_cut_string(rbin)}"
    count = chain.Draw("de", cuts, "goff")
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


def get_rbin_cut_string(rbin):
    if rbin == 0:
        return "(tagwtag>0.45&&tagwtag<=0.5)"

    elif rbin == 1:
        return "(tagwtag>0.375&&tagwtag<=0.45)"

    elif rbin == 2:
        return "(tagwtag>0.25&&tagwtag<=0.375)"

    elif rbin == 3:
        return "(tagwtag>0.1875&&tagwtag<=0.25)"

    elif rbin == 4:
        return "(tagwtag>0.125&&tagwtag<=0.1875)"

    elif rbin == 5:
        return "(tagwtag>0.0625&&tagwtag<=0.125)"

    elif rbin == 6:
        return "(tagwtag>=0.0&&tagwtag<=0.0625)"


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


def process_MC(channels, directory, streams):
    """Read, process, and pretty-print MC."""
    f_sig = {}
    f_scf = {}
    f_bkg = {}

    f_sig["Kpi"] = [0.216, 0.151, 0.174, 0.116, 0.100, 0.090, 0.153]
    f_scf["Kpi"] = [0.222, 0.171, 0.189, 0.125, 0.099, 0.081, 0.113]
    f_bkg["Kpi"] = [0.254, 0.170, 0.181, 0.121, 0.102, 0.083, 0.090]

    f_sig["Kpipi0"] = [0.216, 0.145, 0.174, 0.114, 0.101, 0.092, 0.158]
    f_scf["Kpipi0"] = [0.228, 0.170, 0.187, 0.117, 0.103, 0.084, 0.112]
    f_bkg["Kpipi0"] = [0.262, 0.175, 0.182, 0.124, 0.103, 0.079, 0.074]

    f_sig["K3pi"] = [0.216, 0.148, 0.171, 0.114, 0.102, 0.092, 0.156]
    f_scf["K3pi"] = [0.233, 0.164, 0.187, 0.122, 0.101, 0.084, 0.110]
    f_bkg["K3pi"] = [0.264, 0.171, 0.187, 0.125, 0.101, 0.075, 0.078]

    for channel in channels:

        cr_crscf_fractions = [[] for _ in range(7)]
        cr_fractions = [[] for _ in range(7)]
        scf_fractions = [[] for _ in range(7)]
        bkg_fractions = [[] for _ in range(7)]

        for stream_no, stream in enumerate(streams):
            print_line_bold(channel)
            header = [
                "Stream",
                "Rbin",
                "Yield",
                "Err",
                "Count",
                "Bias",
                "CR/CR+SCF",
                "CR/all",
                "SCF/all",
                "BKG/all",
            ]
            precision = [0, 0, 0, 0, 0, 0, 3, 3, 3, 3]
            table = MarkdownTable(header)
            table.set_precision(precision)

            data = {}
            results_file = f"plots/{channel}_stream{stream}/fit_results.root"
            (
                sig_yield,
                sig_yield_error,
                scf_yield,
                scf_yield_error,
                bkg_yield,
                bkg_yield_error,
            ) = recover_yield(results_file)

            for rbin in range(7):
                cr_crscf, cr, scf, bkg = calculate_fractions(
                    sig_yield * f_sig[channel][rbin],
                    scf_yield * f_scf[channel][rbin],
                    bkg_yield * f_bkg[channel][rbin],
                )
                sig_count = recover_count(
                    f"../data/{channel}/realistic_mc/stream{stream}", rbin
                )
                table.add_row(
                    [
                        stream_no,
                        rbin,
                        sig_yield * f_sig[channel][rbin],
                        sig_yield_error * f_sig[channel][rbin],
                        sig_count,
                        sig_yield * f_sig[channel][rbin] - sig_count,
                        cr_crscf,
                        cr,
                        scf,
                        bkg,
                    ]
                )

                data[f"cr_scf_rbin{rbin}_f"] = round(cr_crscf, 4)
                data[f"cr_rbin{rbin}_f"] = round(cr, 4)
                data[f"scf_rbin{rbin}_f"] = round(scf, 4)

                cr_crscf_fractions[rbin].append(cr_crscf)
                cr_fractions[rbin].append(cr)
                scf_fractions[rbin].append(scf)
                bkg_fractions[rbin].append(bkg)

            table.print()
            print()
            save_to_json(
                data,
                os.path.join(
                    directory, f"stream{stream_no}/{channel}_mc_fractions.json"
                ),
            )

        print_line_bold(f"{channel} Stream Averaged")
        header = [
            "Rbin",
            # "Yield",
            # "Err",
            # "Count",
            # "Bias",
            "CR/CR+SCF",
            "CR/all",
            "SCF/all",
            "BKG/all",
        ]
        precision = [0, 3, 3, 3, 3]
        table = MarkdownTable(header)
        table.set_precision(precision)

        data = {}
        for rbin in range(7):
            table.add_row(
                [
                    rbin,
                    # sig_yield * f_sig[rbin],
                    # sig_yield_error * f_sig[rbin],
                    # sig_count,
                    # sig_yield * f_sig[rbin] - sig_count,
                    sum(cr_crscf_fractions[rbin]) / len(cr_crscf_fractions[rbin]),
                    sum(cr_fractions[rbin]) / len(cr_fractions[rbin]),
                    sum(scf_fractions[rbin]) / len(scf_fractions[rbin]),
                    sum(bkg_fractions[rbin]) / len(bkg_fractions[rbin]),
                ]
            )
            data[f"cr_scf_rbin{rbin}_f"] = round(
                sum(cr_crscf_fractions[rbin]) / len(cr_crscf_fractions[rbin]), 4
            )
            data[f"cr_rbin{rbin}_f"] = round(
                sum(cr_fractions[rbin]) / len(cr_fractions[rbin]), 4
            )
            data[f"scf_rbin{rbin}_f"] = round(
                sum(scf_fractions[rbin]) / len(scf_fractions[rbin]), 4
            )

        save_to_json(data, os.path.join(directory, channel + "_mc_fractions.json"))

        table.print()

        save_to_json(data, os.path.join(directory, channel + "_mc_fractions.json"))


def process_data(channels, directory, randomize=False, correlations=False):
    """Read, process, and pretty-print data."""
    f_sig = {}
    f_scf = {}
    f_bkg = {}

    # BKG fractions are from data sidebands
    f_sig["Kpi"] = [0.216, 0.151, 0.174, 0.116, 0.100, 0.090, 0.153]
    f_scf["Kpi"] = [0.222, 0.171, 0.189, 0.125, 0.099, 0.081, 0.113]
    f_bkg["Kpi"] = [0.268, 0.171, 0.180, 0.108, 0.105, 0.077, 0.091]

    f_sig["Kpipi0"] = [0.216, 0.145, 0.174, 0.114, 0.101, 0.092, 0.158]
    f_scf["Kpipi0"] = [0.228, 0.170, 0.187, 0.117, 0.103, 0.084, 0.112]
    f_bkg["Kpipi0"] = [0.265, 0.172, 0.189, 0.123, 0.098, 0.071, 0.082]

    f_sig["K3pi"] = [0.216, 0.148, 0.171, 0.114, 0.102, 0.092, 0.156]
    f_scf["K3pi"] = [0.233, 0.164, 0.187, 0.122, 0.101, 0.084, 0.110]
    f_bkg["K3pi"] = [0.267, 0.188, 0.162, 0.118, 0.093, 0.080, 0.092]

    print_line_bold("Data")
    for channel in channels:
        print()
        print_line_bold(channel)
        header = [
            "Rbin",
            "Yield",
            "Err",
            "CR/CR+SCF",
            "CR/all",
            "SCF/all",
            "BKG/all",
        ]
        precision = [0, 0, 0, 3, 3, 3, 3]

        table = MarkdownTable(header)
        table.set_precision(precision)

        yields = []
        errors = []

        cr_crscf_fractions = [[] for _ in range(7)]
        cr_fractions = [[] for _ in range(7)]
        scf_fractions = [[] for _ in range(7)]
        bkg_fractions = [[] for _ in range(7)]

        data = {}
        results_file = f"plots/{channel}"
        results_file += "/fit_results.root"
        (
            sig_yield,
            sig_yield_error,
            scf_yield,
            scf_yield_error,
            bkg_yield,
            bkg_yield_error,
        ) = recover_yield(results_file)

        if randomize:
            rng = np.random.default_rng()
            if correlations:
                print("Randomizing with correlations")
                cov = recover_covariance_matrix(results_file)
                f_cr_scf, f_cr_scf_error = recover_cr_scf_fraction(results_file)
                sig_scf_yield, bkg_yield = rng.multivariate_normal(
                    [sig_yield + scf_yield, bkg_yield], cov
                )
                f_cr_scf = rng.normal(f_cr_scf, f_cr_scf_error)

                sig_yield = sig_scf_yield * f_cr_scf
                scf_yield = sig_scf_yield * (1 - f_cr_scf)
            else:
                print("Randomizing without correlations")
                sig_yield = rng.normal(sig_yield, sig_yield_error)
                scf_yield = rng.normal(scf_yield, scf_yield_error)
                bkg_yield = rng.normal(bkg_yield, bkg_yield_error)

        for rbin in range(7):
            cr_crscf, cr, scf, bkg = calculate_fractions(
                sig_yield * f_sig[channel][rbin],
                scf_yield * f_scf[channel][rbin],
                bkg_yield * f_bkg[channel][rbin],
            )

            table.add_row(
                [
                    rbin,
                    sig_yield * f_sig[channel][rbin],
                    sig_yield_error * f_sig[channel][rbin],
                    cr_crscf,
                    cr,
                    scf,
                    bkg,
                ]
            )

            data[f"cr_scf_rbin{rbin}_f"] = round(cr_crscf, 4)
            data[f"cr_rbin{rbin}_f"] = round(cr, 4)
            data[f"scf_rbin{rbin}_f"] = round(scf, 4)

            yields.append(sig_yield)
            errors.append(sig_yield_error)

            cr_crscf_fractions[rbin].append(cr_crscf)
            cr_fractions[rbin].append(cr)
            scf_fractions[rbin].append(scf)
            bkg_fractions[rbin].append(bkg)

        table.print()
        save_to_json(data, os.path.join(directory, channel + "_data_fractions.json"))

    # table.add_row(
    #     [
    #         "Tot/Avg",
    #         sum(yields),
    #         error_of_sum(errors),
    #         weighted_mean(cr_crscf_fractions, yields),
    #         weighted_mean(cr_fractions, yields),
    #         weighted_mean(scf_fractions, yields),
    #         weighted_mean(bkg_fractions, yields),
    #     ]
    # )

    # table.print()

    # data = {
    #     "cr_scf_f": round(weighted_mean(cr_crscf_fractions, yields), 4),
    #     "cr_f": round(weighted_mean(cr_fractions, yields), 4),
    #     "scf_f": round(weighted_mean(scf_fractions, yields), 4),
    # }

    # save_to_json(data, os.path.join(directory, "avg_data_fractions.json"))


def save_to_json(data, filename):
    """Save dictionary as a prettified JSON."""
    # Create the dir if it doesn't exist
    pathlib.Path(os.path.dirname(filename)).mkdir(parents=True, exist_ok=True)
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


def main():
    args = decode_arguments()
    with stdout_redirected():
        # Just to init RooFit to make it print its loading line now
        ROOT.RooRealVar()

    channels = ["Kpi", "Kpipi0", "K3pi"]
    streams = range(0, 6)

    if not os.path.exists(args.directory):
        os.mkdir(args.directory)

    process_MC(channels, args.directory, streams)
    process_data(channels, args.directory)

    for i in range(args.randomize):
        print(f"Generating random config {i}")
        rnd_directory = os.path.join(args.directory, "rnd_" + str(i))
        if not os.path.exists(rnd_directory):
            os.mkdir(rnd_directory)
        process_data(channels, rnd_directory, True, args.correlations)


if __name__ == "__main__":
    main()
