#!/usr/bin/env python3

import sys

import ROOT


def processLine(line):
    _, chi2_reduced, rest = line.rsplit(" ", 2)
    rest = rest[1:-2]
    chi2, ndof = rest.split("/")
    return chi2_reduced, chi2, ndof


channel = sys.argv[2]
file = ROOT.TFile(sys.argv[1], "r")
canvas = file.Get(f"{channel}_all_histpdf_combined_dataset_thetab")
prim = canvas.GetListOfPrimitives()[0]
for prim in prim.GetListOfPrimitives():
    if prim.GetTitle() == "PaveText. A Pave with several lines of text.":
        line = prim.GetLineWith("chi")

chi2_reduced, chi2, ndof = processLine(str(line))
print(chi2_reduced, chi2, ndof)
