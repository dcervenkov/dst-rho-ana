#!/usr/bin/env python

from __future__ import print_function
from math import sqrt
import json
import os
import sys
import ROOT
from contextlib import contextmanager


class MarkdownTable:
    def __init__(self, header):
        assert type(header) == list or type(header) == tuple
        self.header = header
        self.rows = []
        self.precision = [3 for _ in header]

    def set_precision(self, precisions):
        assert len(precisions) == len(self.header)
        assert type(precisions) == list or type(precisions) == tuple
        self.precision = precisions

    def add_row(self, row):
        assert len(self.header) == len(row)
        self.rows.append(row)

    # def add_midline(self):
    #     self.rows.append(["" for _ in range(len(header))])

    def get_length_of_number(self, number, precision):
        number_format = "{:." + str(precision) + "f}"
        return len(number_format.format(number))

    def get_column_widths(self):
        columns = []
        for col in range(len(self.header)):
            cells = []
            for row in self.rows:
                if type(row[col]) is str:
                    cells.append(len(row[col]))
                else:
                    cells.append(
                        self.get_length_of_number(row[col], self.precision[col])
                    )
            cells.append(len(self.header[col]))
            columns.append(max(cells))
        return columns

    def print(self):
        self.columns_widths = self.get_column_widths()
        self.print_row(self.header)
        self.print_separator()
        for row in self.rows:
            self.print_row(row)

    def get_format_string(self, width, precision=None):
        format_string = "{:" + str(width)
        if precision is not None:
            format_string += "." + str(precision) + "f} | "
        else:
            format_string += "} | "
        return format_string

    def print_row(self, row):
        row_text = ""
        for col, cell in enumerate(row):
            format_string = self.get_format_string(
                self.columns_widths[col],
                (None if type(row[col]) is str else self.precision[col]),
            )
            row_text += format_string.format(row[col])

        print(row_text.rstrip(" |"))

    def print_separator(self):
        print("|".rjust(self.columns_widths[0] + 2, "-"), end="")
        for col, _ in enumerate(self.header[1:-1]):
            print("|".rjust(self.columns_widths[col + 1] + 3, "-"), end="")
        print("".ljust(self.columns_widths[-1] + 1, "-"))


def print_line_bold(string):
    print(u"\u001b[1m", end="")
    print(string)
    print(u"\u001b[0m", end="")


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


def common_cuts():
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
        "thetab>0.5"
    )


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
    cr_crscf = float(sig_yield) / (sig_yield + scf_yield)
    cr = float(sig_yield) / (sig_yield + scf_yield + bkg_yield)
    scf = float(scf_yield) / (sig_yield + scf_yield + bkg_yield)
    bkg = float(bkg_yield) / (sig_yield + scf_yield + bkg_yield)
    return (cr_crscf, cr, scf, bkg)


def process_MC(channels, streams):
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
            "cr_scf_f": sum(cr_crscf_fractions) / len(cr_crscf_fractions),
            "cr_f": sum(cr_fractions) / len(cr_fractions),
            "scf_f": sum(scf_fractions) / len(scf_fractions),
        }

        save_to_json(data, os.path.join("results", channel + "_mc_fractions.json"))


def process_data(channels):
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
            "cr_scf_f": cr_crscf_fraction,
            "cr_f": cr_fraction,
            "scf_f": scf_fraction,
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
        "cr_scf_f": weighted_mean(cr_crscf_fractions, yields),
        "cr_f": weighted_mean(cr_fractions, yields),
        "scf_f": weighted_mean(scf_fractions, yields)
    }

    save_to_json(data, os.path.join("results", "avg_data_fractions.json"))


def save_to_json(data, filename):
    with open(filename, "w") as f:
        json.dump(data, f, indent=2, sort_keys=True)


@contextmanager
def stdout_redirected(to=os.devnull):
    """
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
