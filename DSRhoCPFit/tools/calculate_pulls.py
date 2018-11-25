#!/usr/bin/env python3

"""Calculate pulls and their means and stddevs from supplied results files."""

import sys
import pandas as pd
import itertools


VAR_NAMES = ("ap", "apa", "a0", "a0a", "at", "ata", "xp", "x0", "xt", "yp",
             "y0", "yt", "xpb", "x0b", "xtb", "ypb", "y0b", "ytb")


def read_in_files(files):
    print("Reading in " + str(len(files)) + " file(s).\n")
    columns = create_column_list()
    dfs = (pd.read_csv(file, sep=r"\s+\|+\s+|\s+",
                       names=columns, engine='python') for file in files)
    df = pd.concat(dfs, ignore_index=True)
    return df


def create_column_list():
    types = ("true", "fit", "error")
    columns = ("_".join(pair) for pair in itertools.product(VAR_NAMES, types))
    return list(columns)


def add_pulls_to_dataframe(df):
    for var in VAR_NAMES:
        df[var + "_pull"] = (df[var + "_fit"] -
                             df[var + "_true"]) / df[var + "_error"]


def main():
    df = read_in_files(sys.argv[1:])
    add_pulls_to_dataframe(df)
    print("Var | Mean  | Std. Dev.")
    print(":--:|:-----:|:--------:")
    for var in VAR_NAMES:
        print("{:3s} | {:+5.2f} | {:+5.2f}".format(var,
                                                   df[var + "_pull"].mean(),
                                                   df[var + "_pull"].std()))


main()
