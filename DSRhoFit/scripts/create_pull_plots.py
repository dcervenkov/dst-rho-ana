#!/usr/bin/env python

from ROOT import *
import sys

vars = [["ap", "pull |a_{p}|"], ["apa", "pull arg(a_{p})"], ["a0", "pull |a_{0}|"], ["a0a", "pull arg(a_{0})"], ["at", "pull |a_{t}|"], ["ata", "pull arg(a_{t})"]]

if len(sys.argv) != 3:
    print "ERROR: Wrong argument(s)"
    print "USAGE: " + sys.argv[0] + " RESULTS-FILE" + " OUTPUT-DIR"
    quit(0)


results_file = sys.argv[1]
output_dir = sys.argv[2]

gROOT.ProcessLine(".L ../src/tools.cc")
gROOT.ProcessLine("tools::SetupPlotStyle()")
gStyle.SetOptStat(0)
gStyle.SetOptFit(1)

results = []

with open(results_file, 'r') as f:
    for line in f:
        vars_list = [ float(element) for element in line.split() if element != "|"]
        result = []
        for i in range(len(vars)):
            result.append(vars_list[i*3 : i*3+3])
        results.append(result)
    
c1 = TCanvas("c1", "c1", 500, 500)
for varno, var in enumerate(vars):
    histo = TH1D("histo", "histo", 30, -5, +5)
    for result in results:
        if result[varno][2] == 0:
            continue
        histo.Fill((result[varno][1] - result[varno][0])/result[varno][2])
    histo.SetTitle("")
    histo.GetXaxis().SetTitle(var[1])
    histo.Fit("gaus")
    histo.Draw()
    c1.SaveAs(output_dir + "/" + var[0] + ".pdf")
