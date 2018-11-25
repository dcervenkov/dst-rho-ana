#!/usr/bin/env python

"""This script prints pull tables from all supplied ROOT result files."""

import argparse
from rootpy.io import root_open
import sys

def decode_arguments():
    """Decode CLI arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+")
    parser.add_argument("-l", "--latex", action="store_true",
                        help="output in LaTeX format")
    args = parser.parse_args()

    return args.latex, args.files

def main():
    latex, files = decode_arguments()

    for path in files:
        with root_open(path, 'r') as f:
            print(path)
            print
            
            if latex:
                print(f.latex_pull_table.GetTitle())
            else:
                print(f.pull_table.GetTitle())

            print

main()
