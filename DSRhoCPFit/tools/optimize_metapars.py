#!/usr/bin/env python3

import argparse
import itertools
import json
import logging
import os
import subprocess

import numpy as np
import pandas as pd
from pandas.compat import StringIO
from scipy.optimize import minimize

VAR_NAMES = ("ap", "apa", "a0", "a0a", "at", "ata", "xp", "x0", "xt", "yp",
             "y0", "yt", "xpb", "x0b", "xtb", "ypb", "y0b", "ytb")

STEPS = []


def decode_arguments():
    """Decode CLI arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+")
    parser.add_argument("-c", "--config", default="config_opt_temp.json",
                        help="filename to use for temporary JSON configs")
    parser.add_argument("-r", "--result", default="results/opt_temp",
                        help="filename to use for temporary results")
    args = parser.parse_args()

    return args.config, args.result, args.files


def run_fit(config):
    subprocess.check_call(config)


def get_pull_metric(filename):
    df = read_in_file(filename)
    add_pulls_to_dataframe(df)
    pull_keys = [key for key in df.keys() if 'pull' in key]
    pull2 = 0
    for pull_key in pull_keys:
        if not pd.isna(df[pull_key][0]):
            pull2 += df[pull_key][0]**2

    if pull2 == 0:
        pull2 = 100
    return pull2


def read_in_file(file):
    columns = create_column_list()
    df = pd.read_csv(file, sep=r"\s+\|+\s+|\s+",
                     names=columns, engine='python')
    return df


def create_column_list():
    types = ("true", "fit", "error")
    columns = ("_".join(pair) for pair in itertools.product(VAR_NAMES, types))
    return list(columns)


def add_pulls_to_dataframe(df):
    for var in VAR_NAMES:
        df[var + "_pull"] = (df[var + "_fit"] -
                             df[var + "_true"]) / df[var + "_error"]


def create_json(pars, filename):
    config = {
        "fitRanges": {
            "dt": {
                "min": -10,
                "max": 10
            },
            "thetat": {
                "min": 0,
                "max": 3.1415
            },
            "thetab": {
                "min": 0.5,
                "max": 2.95
            }
        },
        "modelParameters": {
            "scf_phit_poly_p2": pars[0],
            "scf_phit_f": pars[1],
            "scf_phit_offset": pars[2],
            "scf_thetat_f": pars[3],
            "scf_thetab_gaus_mu": pars[4],
            "scf_thetab_gaus_sigma_l": pars[5],
            "scf_thetab_gaus_sigma_r": pars[6],
            "scf_thetab_exp_alpha": pars[7],
            "scf_thetab_f": pars[8],
            "bkg_phit_poly_p2": 0.040,
            "bkg_phit_f": 0.660,
            "bkg_phit_offset": 0.103,
            "bkg_thetat_f": -0.207,
            "bkg_thetab_gaus_mu": 2.895,
            "bkg_thetab_gaus_sigma_l": 0.902,
            "bkg_thetab_gaus_sigma_r": 0.089,
            "bkg_thetab_exp_alpha": -2.182,
            "bkg_thetab_f": 0.661
        }
    }

    with open(filename, 'w') as f:
        json.dump(config, f, sort_keys=True, indent=4, separators=(',', ': '))


def objective_function(pars, config_filename, result_filename, data_filenames):
    create_json(pars, config_filename)

    config = ["./DSRhoCPFit",
              "--fit=CRSCF",
              "--efficiency-file=efficiency_new",
              "--config=" + config_filename,
              "--efficiency-model=5",
              "--time-independent",
              "--cpus=2",
              result_filename,
              *data_filenames]

    run_fit(config)
    metric = get_pull_metric(result_filename)
    print("metric = " + str(metric))
    print("pars = [" + ", ".join(map(str, pars)) + "]")
    STEPS.append((metric, pars))
    return metric


def main():
    os.chdir("../.")
    config_file, result_file, data_files = decode_arguments()

    x0 = np.array([0.856, 0.147, 0.056, -0.051,
                   2.885, 0.411, 0.094, -4.63, 0.625])
    result = minimize(objective_function, x0, args=(config_file, result_file, data_files), method='BFGS',
                      options={'disp': True, 'eps': 0.1})
    print(result.x)

    print("Steps:")
    for tuple in STEPS:
        print(tuple[0], tuple[1])


main()
