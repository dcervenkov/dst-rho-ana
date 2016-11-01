#!/usr/bin/env python

from ROOT import *
from math import sqrt
import glob
import sys
import multiprocessing

file_names = glob.glob("data/D0Kpipi0/DSRhoSkim_svd*_on_resonance_evtgen-*_sall.root")
var_name = "d0pi0chi"
var_low_min = 0
var_low_max = 0
var_high_min = 10
var_high_max = 30
var_step = 1
print_num_best = 10
print_num_worst = 3

# Define signal region
de_low = -0.14
de_high = 0.068

signal_region_cut = "(de>"+str(de_low)+")&&(de<"+str(de_high)+")"+"&&(d0mass>1.834&&d0mass<1.890&&dsd0diff>0.144&&dsd0diff<0.147&&csbdtg>-0.6&&pi0chi2<18&&nocand<6)"
signal_cut = "(candsel==3)"
background_cut = "(candsel!=3)"

def tree_count(tree, varname, cut_string):
	events = tree.Draw(var_name, cut_string, "goff")
	return events

def print_row(list):
	string = ""
	for i in range(len(list)):
		string += "{0[" + str(i) + "]:{padding"
		if len(list[i]) <= 6:
			string += "_small}}"
		else:
			string += "_big}}"
	print string.format(list, padding_small = 8, padding_big = 16)

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

print "Optimizing var: " + var_name

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
		sys.stdout.write("\rEvaluated cuts: {}/{}".format(cur_step,total_steps))
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

results.sort(key = lambda x: -x[0])

max_sig = tree.Draw(var_name, signal_cut+"&&"+signal_region_cut,"goff")
max_bkg = tree.Draw(var_name, background_cut+"&&"+signal_region_cut,"goff")

print "\n\nSignal region: (" + str(de_low) + "," + str(de_high) + ")"
print "FOM without cuts: " + str(round(max_sig/sqrt(max_sig+max_bkg),2)) + \
	", SIG: " + str(max_sig) + ", BKG: " + str(max_bkg) + "\n"

header = ["LOW","HIGH","FOM","SIG","SIGLOST","BKG","BKGLOST"]

print "Best " + str(print_num_best) + " fom(s)"
print_row(header)
for result in results[:print_num_best]:
	print_row(prepare_row(result, max_sig, max_bkg))

print "\nWorst " + str(print_num_worst) + " fom(s)"
print_row(header)
for result in results[len(results)-print_num_worst:]:
	print_row(prepare_row(result, max_sig, max_bkg))
