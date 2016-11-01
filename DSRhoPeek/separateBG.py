#!/usr/bin/env /usr/bin/python

from ROOT import *
import sys

filenames = sys.argv[1:]

for filename in filenames:
	output_filename = filename.replace(".root","_bg.root")
	file = TFile(filename,"READ")
	tree = file.Get("h2000")
	output_file = TFile(output_filename, "RECREATE")
	tree_bg = tree.CopyTree("mcflag!=1")
	tree_bg.Write()
	output_file.Close()
	print "Separated background from " + filename
