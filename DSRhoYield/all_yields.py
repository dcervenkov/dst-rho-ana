#!/usr/bin/env python

from ROOT import *


def recover_yield(filename):
    file = TFile(filename)
    sig_plus_cf = float(str(file.Get("all;1").GetLineWith(
        "n_signal_plus_cross_model")).split()[1])
    sig_plus_cf_error = float(
        str(file.Get("all;1").GetLineWith("n_signal_plus_cross_model")).split()[3])
    f_sig_cf = float(str(file.Get("signal_plus_crossfeed;1").GetLineWith(
        "signal_plus_cross_f")).split()[1])

    sig = int(round(sig_plus_cf * f_sig_cf))
    sig_error = int(round(sig_plus_cf_error * f_sig_cf))
    return sig, sig_error


def recover_count(dir):
    chain = TChain("h2000")
    chain.Add(dir + "/*.root")
    count = chain.Draw(
        "de", common_cuts() + "&&(evmcflag==1||evmcflag==7||evmcflag==8)", "goff")
    return count

def common_cuts():
    return "vrusable==1&&vtusable==1&&((vrchi2/vrndf)<50||vrntrk==1)&&((vtchi2/vtndf)<50||vtntrk==1)&&((sqrt(vrerr6)<0.02&&vrntrk>1)||(sqrt(vrerr6)<0.05&&vrntrk==1))&&((sqrt(vterr6)<0.02&&vtntrk>1)||(sqrt(vterr6)<0.05&&vtntrk==1))&&csbdtg>-0.6&&(de>-0.14&&de<0.068)&&(dt>-10&&dt<10)&&thetab>0.5"

# Just to init RooFit to make it print its loading line now
r = RooRealVar()

channels = ["Kpi", "Kpipi0", "K3pi"]
streams = range(0, 6)

for channel in channels:
    print channel
    print "Yield", "+-", "Err", "Count", "Bias"
    yields = []
    errors = []
    counts = []
    biases = []
    for stream in streams:
        sig_yield, sig_yield_error = recover_yield(
            "plots/" + channel + "_stream" + str(stream) + "/fit_results.root")
        sig_count = recover_count("../data/" + channel + "/realistic_mc/stream" + str(stream))
        # sig_yield, sig_yield_error = recover_yield(
        #     "plots/" + channel + "_stream" + str(stream) + "_genMC/fit_results.root")
        # sig_count = recover_count("../data/" + channel + "/mc/stream" + str(stream))
        print sig_yield, "+-", sig_yield_error, sig_count, sig_yield - sig_count
        yields.append(sig_yield)
        errors.append(sig_yield_error)
        counts.append(sig_count)
        biases.append(sig_yield - sig_count)

    print "Averages: "
    print str(sum(yields)/len(yields)), "+-",  str(sum(errors)/len(errors)), str(sum(counts)/len(counts)), str(sum(biases)/len(biases))
    print
