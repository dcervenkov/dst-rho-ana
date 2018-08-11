#!/usr/bin/env python

"""This script prints pull tables from all supplied ROOT result files."""

from rootpy.io import root_open
import sys

if len(sys.argv) == 1:
    print("ERROR: No arguments supplied!")
    print("USAGE: " + sys.argv[0] + " FILE(S)...")
    sys.exit(1)

for path in sys.argv[1:]:
    with root_open(path, 'r') as f:
        print(path)
        print()
        print(f.pull_table.GetTitle())
        print()
