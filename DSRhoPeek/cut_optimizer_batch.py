#!/usr/bin/env python

from ROOT import *
from math import sqrt
import glob
import sys
import os
import multiprocessing

# How many results to print
print_num_best = 10
print_num_worst = 3

# Define signal region
de_low = -0.14
de_high = 0.068

if len(sys.argv) != 3:
    print "ERROR: Incorrect number of parameters."
    print "USAGE: %s DIR LOG-FILE" % (sys.argv[0])
    sys.exit(1)

dir = sys.argv[1]
log_file = sys.argv[2]

if not os.path.exists(dir):
    print "ERROR: Directory '%s' doesn't exist." % dir
    sys.exit(2)

file_names = glob.glob(dir.rstrip('/')+"/DSRhoSkim_svd*_on_resonance_evtgen-*_sall.root")

signal_region_cut = "(de>"+str(de_low)+")&&(de<"+str(de_high)+")"+""
signal_cut = "(candsel==4)"
background_cut = "(candsel!=4)"

def tree_count(tree, varname, cut_string):
	events = tree.Draw(var_name, cut_string, "goff")
	return events

def print_row(list, file):
	string = ""
	for i in range(len(list)):
		string += "{0[" + str(i) + "]:{padding"
		if len(list[i]) <= 6:
			string += "_small}}"
		else:
			string += "_big}}"
	string += "\n"
	file.write(string.format(list, padding_small = 8, padding_big = 16))

def prepare_row(result, max_sig, max_bkg):
	fom = result[0]
	sig = result[1]
	bkg = result[2]
	cut_low = result[3]
	cut_high = result[4]

	list = [
		str(cut_low), 
		str(cut_high), 
		str(round(fom,2)), 
		str(sig), 
		str(max_sig - sig) + " (" + str(round(100*(1 - sig/float(max_sig)),2)) + "%)",
		str(bkg), 
		str(max_bkg - bkg) + " (" + str(round(100*(1 - bkg/float(max_bkg)),2)) + "%)"]

	return list

def create_steps(mn,mx,step):
    step_list = []
    val = mn
    while val <= mx:
	step_list.append(val)
	# Round because of the floating number 
	# binary representation problem
	val = round(val + step,6)
    return step_list

with open('cuts.cfg', 'r') as f:
    with open(log_file, 'w') as log:
        for line in f:
            var_name, var_low_min, var_low_max, var_high_min, var_high_max, var_step = line.split()
            var_low_min = float(var_low_min)
            var_low_max = float(var_low_max)
            var_high_min = float(var_high_min)
            var_high_max = float(var_high_max)
            var_step = float(var_step)

            log.write("\n\n********* Optimizing var: " + var_name + " *********")

            low_steps  = create_steps(var_low_min, var_low_max, var_step)
            high_steps = create_steps(var_high_min, var_high_max, var_step)

            tree = TChain("h2000")
            for file_name in file_names:
                tree.Add(file_name)

            results = []

            total_steps = len(low_steps)*len(high_steps)
            cur_step = 0

            pool = multiprocessing.Pool(processes=2)

            for cut_low in low_steps:
                for cut_high in high_steps:
                    cur_step += 1
                    sys.stdout.write("\rVariable: {} Evaluated cuts: {}/{}".format(var_name.ljust(8),cur_step,total_steps))
                    sys.stdout.flush()
                    cut_string = "("+var_name+">"+str(cut_low)+")&&("+var_name+"<"+str(cut_high)+")"
                    psig = pool.apply_async(tree_count, args=(tree, var_name, signal_cut+"&&"+signal_region_cut+"&&"+cut_string))
                    pbkg = pool.apply_async(tree_count, args=(tree, var_name, background_cut+"&&"+signal_region_cut+"&&"+cut_string))
                    sig = psig.get()
                    bkg = pbkg.get()
                    if bkg != 0 and sig != 0:
                        fom = abs(sig/sqrt(sig+bkg))
                    else:
                        fom = 0
                    result = [fom, sig, bkg, cut_low, cut_high]
                    results.append(result)

            # Just to print a new line, so the progress counter is not over-written
            # by the next variable
            print ""

            results.sort(key = lambda x: -x[0])

            max_sig = tree.Draw(var_name, signal_cut+"&&"+signal_region_cut,"goff")
            max_bkg = tree.Draw(var_name, background_cut+"&&"+signal_region_cut,"goff")

            log.write("\n\nSignal region: (" + str(de_low) + "," + str(de_high) + ")\n")
            log.write("FOM without cuts: " + str(round(max_sig/sqrt(max_sig+max_bkg),2)) + \
                ", SIG: " + str(max_sig) + ", BKG: " + str(max_bkg) + "\n\n")

            header = ["LOW","HIGH","FOM","SIG","SIGLOST","BKG","BKGLOST"]

            log.write("Best " + str(print_num_best) + " fom(s)\n")
            print_row(header, log)
            for result in results[:print_num_best]:
                print_row(prepare_row(result, max_sig, max_bkg), log)

            log.write("\nWorst " + str(print_num_worst) + " fom(s)\n")
            print_row(header, log)
            for result in results[len(results)-print_num_worst:]:
                print_row(prepare_row(result, max_sig, max_bkg), log)
