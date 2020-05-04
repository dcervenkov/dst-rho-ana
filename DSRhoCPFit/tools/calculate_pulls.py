#!/usr/bin/env python3

"""Calculate pulls and their means and stddevs from supplied results files."""

import pandas as pd
import itertools
import argparse


VAR_NAMES = (
    "ap",
    "apa",
    "a0",
    "a0a",
    "at",
    "ata",
    "xp",
    "x0",
    "xt",
    "yp",
    "y0",
    "yt",
    "xpb",
    "x0b",
    "xtb",
    "ypb",
    "y0b",
    "ytb",
)

VAR_NAMES_LATEX = (
    r"|a_{\parallel}|",
    r"\arg(a_{\parallel})",
    r"|a_{0}|",
    r"\arg(a_{0})",
    r"|a_{\perp}|",
    r"\arg(a_{\perp})",
    r"x_{\parallel}",
    r"x_{0}",
    r"x_{\perp}",
    r"y_{\parallel}",
    r"y_{0}",
    r"y_{\perp}",
    r"\bar x_{\parallel}",
    r"\bar x_{0}",
    r"\bar x_{\perp}",
    r"\bar y_{\parallel}",
    r"\bar y_{0}",
    r"\bar y_{\perp}",
)

LATEX_HEADER = r"""\begin{table}
  \footnotesize
  \begin{tabular}{lSSr}
    \toprule
    {Var} & {Mean} & {Std. Dev.} \\
    \midrule"""

LATEX_FOOTER = r"""    \bottomrule
  \end{tabular}
\end{table}"""


def decode_arguments():
    """Decode CLI arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+", help="input files")
    parser.add_argument(
        "-l", "--latex", action="store_true", help="output a LaTeX table"
    )
    args = parser.parse_args()
    return args.files, args.latex


def read_in_files(files):
    print("Reading in " + str(len(files)) + " file(s).\n")
    columns = create_column_list()
    dfs = (
        pd.read_csv(file, sep=r"\s+\|+\s+|\s+", names=columns, engine="python")
        for file in files
    )
    df = pd.concat(dfs, ignore_index=True)
    return df


def create_column_list():
    types = ("true", "fit", "error")
    columns = ("_".join(pair) for pair in itertools.product(VAR_NAMES, types))
    return list(columns)


def add_pulls_to_dataframe(df):
    for var in VAR_NAMES:
        df[var + "_pull"] = (df[var + "_fit"] - df[var + "_true"]) / df[var + "_error"]


def print_table(df, latex):
    if latex:
        print(LATEX_HEADER)
    else:
        print("Var | Mean  | Std. Dev.")
        print(":--:|:-----:|:--------:")

    for i, var in enumerate(VAR_NAMES):
        if pd.isna(df[var + "_pull"].mean()):
            continue
        if latex:
            print(
                "    ${:19s}$ & {:+5.2f} & {:+5.2f} \\\\".format(
                    VAR_NAMES_LATEX[i],
                    df[var + "_pull"].mean(),
                    df[var + "_pull"].std(),
                )
            )
        else:
            print(
                "{:3s} | {:+5.2f} | {:+5.2f}".format(
                    var, df[var + "_pull"].mean(), df[var + "_pull"].std()
                )
            )
    if latex:
        print(LATEX_FOOTER)


def main():
    files, latex = decode_arguments()
    df = read_in_files(files)
    add_pulls_to_dataframe(df)
    print_table(df, latex)


main()
