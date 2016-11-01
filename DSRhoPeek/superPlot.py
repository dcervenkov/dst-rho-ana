#!/usr/bin/env python

from ROOT import *
import math
import os
import subprocess
import sys
import glob

if len(sys.argv) < 3:
    print "ERROR: Not enough arguments"
    print "USAGE: " + sys.argv[0] + " PREVIEW DIR(s)..."
    sys.exit(1)

if sys.argv[1] not in ["0","1"]:
    print "ERROR: Unrecognized PREVIEW value; use 0 or 1"
    sys.exit(2)

preview = bool(int(sys.argv[1]))

log_scale = False 

if preview:
    image_format = ".png"
    rowsOnPage = 3
    canvas_width = 1500
    canvas_height = 900
else:
    image_format = ".pdf"
    rowsOnPage = 1
    canvas_width = 300
    canvas_height = 300

gROOT.Reset()
gROOT.ProcessLine(".L include/Colors.h")
gROOT.ProcessLine("setColors()")

def savePage():
    c.Clear()
    c.Divide(len(files),rowsOnPage,0.000001,0.000001)
    for column in range(len(files)):
        file = files[column]

        histogramList = file.GetListOfKeys()
        for key in histogramList:
            if key.GetClassName() != "THStack":
                histogramList.Remove(key)
        histogramList.Sort()
        
        names = []

        for row in range(rowsOnPage):
            maximum = 0.
            if (page*rowsOnPage+row >= len(histogramList)):
                continue
            c.cd(1 + row*len(files) + column)
            if log_scale:
				gPad.SetLogy()
            name = histogramList[page*rowsOnPage+row].GetName()
            names.append(name)
            histo = file.Get(name)
            #histo.SetMaximum(150)
            hists = histo.GetHists()
            for hist in hists:
                if image_format == ".pdf":
                    hist.SetLineWidth(1)
                else:
                    hist.SetLineWidth(2)
            histo.Draw("nostack")

    page_name = "_".join(names)
    if not preview:
        data_types = ["tc1","tc0","charged","charm","uds"]
        plot_types = ["signal","mixed","charged","charm","uds"]
        for i in range(len(data_types)):
            if file.GetName().find(data_types[i]) != -1:
                page_name += "_" + plot_types[i]
                break

    c.SaveAs(os.path.join(svd,channel)+"/" + page_name + image_format)
    c.SetName("c" + str(page+1))
    c.Write()


dirs = [arg.rstrip('/') for arg in sys.argv[2:]]

print "Processing the following dirs:"
print "\n".join(dirs)

for directory in dirs:
    if not os.path.exists(directory):
        print "ERROR: Directory '" + directory + "' doesn't exist! Exiting..."
        sys.exit(1)
    channel = os.path.basename(directory)

    if not glob.glob(directory + "/*mixed*result*.root"):
        subprocess.call("for FILE in " + directory +"/*mixed*sall.root; do Debug/DSRhoPeek $FILE 1; Debug/DSRhoPeek $FILE 0; done", shell=True)

    if not glob.glob(directory + "/*charm*result*.root"):
        subprocess.call("for FILE in " + directory +"/*char*sall.root; do Debug/DSRhoPeek $FILE; done", shell=True)

    if not glob.glob(directory + "/*uds*result*.root"):
        subprocess.call("for FILE in " + directory +"/*uds*sall.root; do Debug/DSRhoPeek $FILE; done", shell=True)


    types = ["evtgen-mixed","evtgen-charged","evtgen-charm","evtgen-uds"]
#    types = ["evtgen-mixed"]#,"evtgen-charged","evtgen-charm","evtgen-uds"]

    c = TCanvas("c", "c", canvas_width, canvas_height)

    for svd in ["svd1", "svd2"]:
        if not os.path.exists(os.path.join(svd,channel)):
            os.makedirs(os.path.join(svd,channel))

        files = []

        for type in types:
            if type == "evtgen-mixed": 
                files.append(TFile(directory+"/DSRhoSkim_"+svd+"_on_resonance_"+type+"_sall_results_tc1.root"))
                files.append(TFile(directory+"/DSRhoSkim_"+svd+"_on_resonance_"+type+"_sall_results_tc0.root"))
                #files.append(TFile(directory+"/DSRhoSkim_"+svd+"_on_resonance_"+type+"_sall_results.root"))
            else:
                files.append(TFile(directory+"/DSRhoSkim_"+svd+"_on_resonance_"+type+"_sall_results.root"))
        
        rootFile = TFile(os.path.join(svd,channel)+"/c.root","RECREATE")

        file = files[0]
        histogramList = file.GetListOfKeys()

        for key in histogramList:
            if key.GetClassName() != "THStack":
                histogramList.Remove(key)

        histogramList.Sort()
        numPages = int(math.ceil(len(histogramList)/float(rowsOnPage)))

        for page in range(numPages):
            if preview:
                savePage()
            else:
                original_files = files
                for file in original_files:
                    files = []
                    files.append(file)
                    savePage()
                files = original_files
                    

        rootFile.Close()

