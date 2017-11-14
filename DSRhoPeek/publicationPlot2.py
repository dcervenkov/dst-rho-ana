#!/usr/bin/env python

import ROOT
import os

histograms = []
#histograms.append(["d0massbf","M(D^{0}) [GeV]"])
#histograms.append(["de","#DeltaE [GeV]"])
#histograms.append(["dsd0diff","M(D*)-M(D^{0}) [GeV]"])
#histograms.append(["pi0mass","M(#pi^{0}) [GeV]"])
#histograms.append(["rhomass","M(#rho) [GeV]"])
#histograms.append(["vrvtxz", "vrvtxz", "", "z_{sig} [cm]", 0, "", ""])
#histograms.append(["vrvtxz_log", "vrvtxz", "", "z_{sig} [cm]", 1, "", ""])
#histograms.append(["vtvtxz", "vtvtxz", "vtgexist==1", "z_{tag} [cm]", 0, "", ""])
#histograms.append(["vtvtxz_log", "vtvtxz", "", "z_{tag} [cm]", 1, "", ""])
#histograms.append(["vrgvtxz", "vrgvtxz/10", "vrgexist==1", "z^{gen}_{sig} [cm]", 0, "", ""])
#histograms.append(["vrgvtxz_log", "vrgvtxz/10", "vrgexist==1", "z^{gen}_{sig} [cm]", 1, "", ""])
#histograms.append(["vtgvtxz", "vtgvtxz/10", "vtgexist==1", "z^{gen}_{tag} [cm]", 0, "", ""])
#histograms.append(["vtgvtxz_log", "vtgvtxz/10", "vtgexist==1", "z^{gen}_{tag} [cm]", 1, "", ""])
histograms.append(["vrchi2", "vrchi2", "vrgexist==1&&evmcflag==1&&vrusable==1&&vtusable==1&&vrchi2/vrndf<50", "#chi_{sig}^{2}", 0, 0, 250])
histograms.append(["vtchi2", "vtchi2", "vrgexist==1&&evmcflag==1&&vrusable==1&&vtusable==1&&vtchi2/vtndf<50", "#chi_{tag}^{2}", 0, 0, 250])
histograms.append(["vrh", "vrchi2/vrndf", "vrgexist==1&&evmcflag==1&&vrusable==1&&vtusable==1&&vrchi2/vrndf<50", "h_{sig}", 0, 0, 50])
histograms.append(["vrh_log", "vrchi2/vrndf", "vrgexist==1&&evmcflag==1&&vrusable==1&&vtusable==1&&vrchi2/vrndf<50", "h_{sig}", 1, 0, 50])
histograms.append(["vth", "vtchi2/vtndf", "vrgexist==1&&evmcflag==1&&vrusable==1&&vtusable==1&&vtchi2/vtndf<50", "h_{tag}", 0, 0, 50])
histograms.append(["vth_log", "vtchi2/vtndf", "vrgexist==1&&evmcflag==1&&vrusable==1&&vtusable==1&&vtchi2/vtndf<50", "h_{tag}", 1, 0, 50])
histograms.append(["vrvtxz_res", "(vrvtxz-vrgvtxz/10)*10000", "vrgexist==1&&evmcflag==1&&vrusable==1&&vtusable==1&&vrchi2/vrndf<50", "z_{sig}-z^{gen}_{sig} [#mum]", 0, -200, 200])
#histograms.append(["vrvtxz_res_log", "vrvtxz-vrgvtxz/10", "vrgexist==1", "z_{sig}-z^{gen}_{sig} [cm]", 1, "", ""])
histograms.append(["vrvtxz_error", "sqrt(vrerr6)*10000", "evmcflag==1&&vrusable==1&&vtusable==1&&vrchi2/vrndf<50", "z_{sig}-z^{gen}_{sig} [#mum]", 0, 0, 200])
histograms.append(["vrvtxz_pull", "(vrvtxz-vrgvtxz/10)/sqrt(vrerr6)", "vrgexist==1&&evmcflag==1&&vrusable==1&&vtusable==1&&vrchi2/vrndf<50", "z^{pull}_{sig}", 0, -10, 10])
histograms.append(["vtvtxz_res", "(vtvtxz-vtgvtxz/10)*10000", "vtgexist==1&&evmcflag==1&&vrusable==1&&vtusable==1&&vrchi2/vrndf<50", "z_{tag}-z^{gen}_{tag} [#mum]", 0, -200, 200])
histograms.append(["vtvtxz_res_log", "(vtvtxz-vtgvtxz/10)*10000", "vtgexist==1&&evmcflag==1&&vrusable==1&&vtusable==1&&vrchi2/vrndf<50", "z_{tag}-z^{gen}_{tag} [#mum]", 1, -2000, 2000])
histograms.append(["vtvtxz_error", "sqrt(vterr6)*10000", "evmcflag==1&&vrusable==1&&vtusable==1&&vtchi2/vtndf<50", "z_{sig}-z^{gen}_{sig} [#mum]", 0, 0, 200])
histograms.append(["vtvtxz_pull", "(vtvtxz-vtgvtxz/10)/sqrt(vterr6)", "vtgexist==1&&evmcflag==1&&vrusable==1&&vtusable==1&&vtchi2/vtndf<50", "z^{pull}_{tag}", 0, -10, 10])
histograms.append(["vtvtxz_pull_chi2cut", "(vtvtxz-vtgvtxz/10)/sqrt(vterr6)", "vtgexist==1&&evmcflag==1&&vrusable==1&&vtusable==1&&vrchi2/vrndf<50&&vtchi2<5", "z^{pull}_{tag}", 0, -10, 10])
#histograms.append(["dz", "vrvtxz-vtvtxz", "", "#Deltaz [cm]", 0, "-0.1", "0.1"])
#histograms.append(["dz_log", "vrvtxz-vtvtxz", "", "#Deltaz [cm]", 1, "", ""])
histograms.append(["dt_a", "(vrvtxz-vtvtxz)/(0.425*0.0299792458)", "(vrvtxz-vtvtxz)/(0.425*0.0299792458)<10&&evmcflag==1&&vrusable==1&&vtusable==1&&vrchi2/vrndf<50&&vtchi2/vtndf<50&&brecflav==1&&btagmcli<0", "#Deltat [ps]", 0, -10, 10])
histograms.append(["dt_b", "(vrvtxz-vtvtxz)/(0.425*0.0299792458)", "(vrvtxz-vtvtxz)/(0.425*0.0299792458)<10&&evmcflag==1&&vrusable==1&&vtusable==1&&vrchi2/vrndf<50&&vtchi2/vtndf<50&&brecflav==1&&btagmcli>0", "#Deltat [ps]", 0, -10, 10])
histograms.append(["dt_ab", "(vrvtxz-vtvtxz)/(0.425*0.0299792458)", "(vrvtxz-vtvtxz)/(0.425*0.0299792458)<10&&evmcflag==1&&vrusable==1&&vtusable==1&&vrchi2/vrndf<50&&vtchi2/vtndf<50&&brecflav==-1&&btagmcli>0", "#Deltat [ps]", 0, -10, 10])
histograms.append(["dt_bb", "(vrvtxz-vtvtxz)/(0.425*0.0299792458)", "(vrvtxz-vtvtxz)/(0.425*0.0299792458)<10&&evmcflag==1&&vrusable==1&&vtusable==1&&vrchi2/vrndf<50&&vtchi2/vtndf<50&&brecflav==-1&&btagmcli<0", "#Deltat [ps]", 0, -10, 10])
histograms.append(["dt_a_wtag", "(vrvtxz-vtvtxz)/(0.425*0.0299792458)", "(vrvtxz-vtvtxz)/(0.425*0.0299792458)<10&&evmcflag==1&&vrusable==1&&vtusable==1&&vrchi2/vrndf<50&&vtchi2/vtndf<50&&brecflav==1&&tagflavo<0", "#Deltat [ps]", 0, -10, 10])
histograms.append(["dt_b_wtag", "(vrvtxz-vtvtxz)/(0.425*0.0299792458)", "(vrvtxz-vtvtxz)/(0.425*0.0299792458)<10&&evmcflag==1&&vrusable==1&&vtusable==1&&vrchi2/vrndf<50&&vtchi2/vtndf<50&&brecflav==1&&tagflavo>0", "#Deltat [ps]", 0, -10, 10])
histograms.append(["dt_ab_wtag", "(vrvtxz-vtvtxz)/(0.425*0.0299792458)", "(vrvtxz-vtvtxz)/(0.425*0.0299792458)<10&&evmcflag==1&&vrusable==1&&vtusable==1&&vrchi2/vrndf<50&&vtchi2/vtndf<50&&brecflav==-1&&tagflavo>0", "#Deltat [ps]", 0, -10, 10])
histograms.append(["dt_bb_wtag", "(vrvtxz-vtvtxz)/(0.425*0.0299792458)", "(vrvtxz-vtvtxz)/(0.425*0.0299792458)<10&&evmcflag==1&&vrusable==1&&vtusable==1&&vrchi2/vrndf<50&&vtchi2/vtndf<50&&brecflav==-1&&tagflavo<0", "#Deltat [ps]", 0, -10, 10])
histograms.append(["dt_a_goodtag", "(vrvtxz-vtvtxz)/(0.425*0.0299792458)", "(vrvtxz-vtvtxz)/(0.425*0.0299792458)<10&&evmcflag==1&&vrusable==1&&vtusable==1&&vrchi2/vrndf<50&&vtchi2/vtndf<50&&brecflav==1&&tagflavo<0&&abs(tagqr)>0.5", "#Deltat [ps]", 0, -10, 10])
histograms.append(["dt_b_goodtag", "(vrvtxz-vtvtxz)/(0.425*0.0299792458)", "(vrvtxz-vtvtxz)/(0.425*0.0299792458)<10&&evmcflag==1&&vrusable==1&&vtusable==1&&vrchi2/vrndf<50&&vtchi2/vtndf<50&&brecflav==1&&tagflavo>0&&abs(tagqr)>0.5", "#Deltat [ps]", 0, -10, 10])
histograms.append(["dt_ab_goodtag", "(vrvtxz-vtvtxz)/(0.425*0.0299792458)", "(vrvtxz-vtvtxz)/(0.425*0.0299792458)<10&&evmcflag==1&&vrusable==1&&vtusable==1&&vrchi2/vrndf<50&&vtchi2/vtndf<50&&brecflav==-1&&tagflavo>0&&abs(tagqr)>0.5", "#Deltat [ps]", 0, -10, 10])
histograms.append(["dt_bb_goodtag", "(vrvtxz-vtvtxz)/(0.425*0.0299792458)", "(vrvtxz-vtvtxz)/(0.425*0.0299792458)<10&&evmcflag==1&&vrusable==1&&vtusable==1&&vrchi2/vrndf<50&&vtchi2/vtndf<50&&brecflav==-1&&tagflavo<0&&abs(tagqr)>0.5", "#Deltat [ps]", 0, -10, 10])
histograms.append(["tagqr", "tagqr", "evmcflag==1&&vrusable==1&&vtusable==1&&vrchi2/vrndf<50&&vtchi2/vtndf<50", "qr", 0, "", ""])
histograms.append(["tagwtag", "tagwtag", "evmcflag==1&&vrusable==1&&vtusable==1&&vrchi2/vrndf<50&&vtchi2/vtndf<50", "w_{tag}", 0, "", ""])
histograms.append(["tageff", "(1-2*tagwtag)^2", "evmcflag==1&&vrusable==1&&vtusable==1&&vrchi2/vrndf<50&&vtchi2/vtndf<50", "#epsilon_{eff}", 0, "", ""])

data = []
#data.append(["data/signal/DSRho-mdst_D0Kpi_svd2_baseline.root",1,"baseline"])
#data.append(["data/signal/DSRho-mdst_D0Kpi_svd2_kfitter_baseline.root",1,"kfitter_baseline"])
#data.append(["data/signal/DSRho-mdst_D0Kpi_svd2_wo_iptube.root",1,"wo_iptube"])
#data.append(["data/signal/DSRho-mdst_D0Kpi_svd2.root",1,"svd2"])
#data.append(["data/signal/DSRho-mdst_D0Kpi_exa_svd2.root",1,"exa"])
#data.append(["data/signal/DSRho-mdst_D0Kpi_exa_new_svd2.root",1,"exa_new"])

#data.append(["data/signal/DSRho-mdst_D0Kpi_exa_small_svd1.root",1,"small"])
#data.append(["data/signal/DSRho-mdst_D0Kpi_exa_woip_small_svd1.root",1,"woip_small"])
#data.append(["data/signal/DSRho-mdst_D0Kpi_exa_vnew_small_svd1.root",1,"vnew_small"])

data.append(["../data/DSRho-mdst_basf2_mod_real_unmod_1.root",1,"unmod"])

#extra_text = "_signalMC"
extra_text = ""

plot_dir = "plots"

image_format = ".pdf"
canvas_width = 500
canvas_height = 500

legend_text = ["other", "CR", "bad #rho", "bad D^{0}", "bad D*", "bad D*, D^{0}, #rho", "miss. p. in #rho", "non-res. #pi#pi^{0}", "D**"]

min_num_to_print = 200

ROOT.gROOT.Reset()
ROOT.gROOT.ProcessLine(".L src/colors.cc")
ROOT.gROOT.ProcessLine("colors::setColors()")
ROOT.gROOT.ProcessLine(".L src/tools.cc")
ROOT.gROOT.ProcessLine("tools::SetupPlotStyle()")
ROOT.gStyle.SetOptStat(0);

output_file = ROOT.TFile(os.path.join(plot_dir, "plots.root"),"RECREATE")

c = ROOT.TCanvas("c", "c", canvas_width, canvas_height)
c.cd()

for histogram_name, histogram_formula, histogram_cuts, axis_title, logaritmic, min, max in histograms:
    for file_name, hist_num, channel in data:
        file = ROOT.TFile(file_name)
        tree = file.Get("h2000")
        c.SetLogy(logaritmic)

        if min != "":
            if histogram_cuts != "":
                histogram_cuts += "&&"
            histogram_cuts += histogram_formula + ">" + str(min)
        if max != "":
            if histogram_cuts != "":
                histogram_cuts += "&&"
            histogram_cuts += histogram_formula + "<" + str(max)

        if channel == "exa_new":
            to_same = "same"
        else:
            to_same = ""

        if not tree.Draw(histogram_formula, histogram_cuts, to_same):
            print("WARNING: Nothing to print for {} with {} cuts!".format(histogram_formula, histogram_cuts))
            continue
        histo = ROOT.gPad.GetPrimitive("htemp")
        if channel == "exa_new":
            histo.SetLineColor(ROOT.kRed)
        else:
            histo.SetLineColor(ROOT.kBlack)
        histo.SetStats(True)
        histo.GetXaxis().SetTitle(axis_title)
        histo.SetTitle("")
        histo.SetName(histogram_name)
        histo.Rebin(2)

        output_file.cd()
        histo.Write()

        #c.SetGridx()
        #c.SetGridy()
        c.SetTickx()
        c.SetTicky()
        c.SaveAs(os.path.join(plot_dir, histogram_name + "_" + channel + extra_text + image_format))

output_file.Close()

