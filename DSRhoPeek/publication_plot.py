#!/usr/bin/env python

import os
import ROOT


def create_histogram_list():
    histograms = []
    #histograms.append(["d0massbf","M(D^{0}) [GeV]"])
    # histograms.append(["de","#DeltaE [GeV]"])
    #histograms.append(["dsd0diff","M(D*)-M(D^{0}) [GeV]"])
    # histograms.append(["pi0mass","M(#pi^{0}) [GeV]"])
    # histograms.append(["rhomass","M(#rho) [GeV]"])
    histograms.append(["thetab", "#theta_{b} [rad]"])

    return histograms


def create_data_list():
    data = []
    # data.append(["data/signal/DSRho-mdst_D0Kpi_svd2_results.root",1,"Kpi"])
    # data.append(["data/signal/DSRho-mdst_D0Kpipi0_svd2_results.root",0,"Kpipi0"])
    # data.append(["data/signal/DSRho-mdst_D0K3pi_svd2_results.root",1,"K3pi"])

    # data.append(["data/D0Kpi/mc/DSRhoSkim_svd2_on_resonance_evtgen-mixed_sall_results_tc1.root", 1, "Kpi_signal"])
    # data.append(["data/D0Kpi/mc/DSRhoSkim_svd2_on_resonance_evtgen-mixed_sall_results_tc0.root", 1, "Kpi_mixed"])
    # data.append(["data/D0Kpi/mc/DSRhoSkim_svd2_on_resonance_evtgen-charged_sall_results.root", 1, "Kpi_charged"])
    # data.append(["data/D0Kpi/mc/DSRhoSkim_svd2_on_resonance_evtgen-charm_sall_results.root", 1, "Kpi_charm"])
    # data.append(["data/D0Kpi/mc/DSRhoSkim_svd2_on_resonance_evtgen-uds_sall_results.root", 1, "Kpi_uds"])

    # data.append(["data/D0Kpipi0/mc/DSRhoSkim_svd2_on_resonance_evtgen-mixed_sall_results_tc1.root",1,"Kpipi0_signal"])
    # data.append(["data/D0Kpipi0/mc/DSRhoSkim_svd2_on_resonance_evtgen-mixed_sall_results_tc0.root",1,"Kpipi0_mixed"])
    # data.append(["data/D0Kpipi0/mc/DSRhoSkim_svd2_on_resonance_evtgen-charged_sall_results.root",1,"Kpipi0_charged"])
    # data.append(["data/D0Kpipi0/mc/DSRhoSkim_svd2_on_resonance_evtgen-charm_sall_results.root",1,"Kpipi0_charm"])
    # data.append(["data/D0Kpipi0/mc/DSRhoSkim_svd2_on_resonance_evtgen-uds_sall_results.root",1,"Kpipi0_uds"])

    # data.append(["data/D0K3pi/mc/DSRhoSkim_svd2_on_resonance_evtgen-mixed_sall_results_tc1.root",1,"K3pi_signal"])
    # data.append(["data/D0K3pi/mc/DSRhoSkim_svd2_on_resonance_evtgen-mixed_sall_results_tc0.root",1,"K3pi_mixed"])
    # data.append(["data/D0K3pi/mc/DSRhoSkim_svd2_on_resonance_evtgen-charged_sall_results.root",1,"K3pi_charged"])
    # data.append(["data/D0K3pi/mc/DSRhoSkim_svd2_on_resonance_evtgen-charm_sall_results.root",1,"K3pi_charm"])
    # data.append(["data/D0K3pi/mc/DSRhoSkim_svd2_on_resonance_evtgen-uds_sall_results.root",1,"K3pi_uds"])

    data.append(["../data/DSRho-mdst_basf2_mod_real_unmod_results.root", 1, "Kpi"])

    return data


def setup_plot_style():
    ROOT.gROOT.Reset()
    ROOT.gROOT.ProcessLine(".L src/colors.cc")
    ROOT.gROOT.ProcessLine("colors::setColors()")
    ROOT.gROOT.ProcessLine(".L src/tools.cc")
    ROOT.gROOT.ProcessLine("tools::SetupPlotStyle()")
    ROOT.gStyle.SetOptStat(0)


def make_publication_plots():
    setup_plot_style()
    #extra_text = "_signalMC"
    extra_text = ""

    plot_dir = "plots"

    image_format = ".pdf"
    canvas_width = 500
    canvas_height = 500

    legend_text = ["other", "CR", "bad #rho", "bad D^{0}", "bad D*", "bad D*, D^{0}, #rho",
                   "miss. p. in #rho", "non-res. #pi#pi^{0}", "D**"]

    min_num_to_print = 200

    canvas = ROOT.TCanvas("canvas", "canvas", canvas_width, canvas_height)
    canvas.cd()

    histograms = create_histogram_list()
    data = create_data_list()

    for histogram_name, axis_title in histograms:
        for file_name, hist_num, channel in data:
            file = ROOT.TFile(file_name)
            histo = file.Get(histogram_name + ";1")
            # histo.SetMaximum(150)
            hists = histo.GetHists()
            existing_hists = []
            num_hists_to_draw = 0
            for hist in hists:
                existing_hists.append(int(str(hist.GetName()).split("_")[3]))
                if hist.GetEntries() > min_num_to_print:
                    num_hists_to_draw += 1

            histo.Draw("nostack")
            histo.GetXaxis().SetTitle(axis_title)
            histo.SetTitle("")

            mylegend = ROOT.TLegend(
                0.75, 0.9 - 0.04 * num_hists_to_draw, 0.85, 0.9)
            #mylegend = TLegend(0.15,0.9 - 0.04*num_hists_to_draw,0.25,0.9)
            mylegend.SetBorderSize(0)
            mylegend.SetFillColor(0)
            mylegend.SetTextFont(43)
            mylegend.SetTextSize(16)

            plot_other = False
            for pos, num in enumerate(existing_hists):
                if hists[pos].GetEntries() > min_num_to_print:
                    # A Q&D hack to print 'other' last
                    if num == 0:
                        plot_other = True
                        continue

                    mylegend.AddEntry(hists[pos], legend_text[num], "L")

            if plot_other:
                mylegend.AddEntry(hists[0], legend_text[0], "L")

            mylegend.Draw()

            # c.SetGridx()
            # c.SetGridy()
            canvas.SetTickx()
            canvas.SetTicky()
            canvas.SaveAs(os.path.join(plot_dir, histogram_name +
                                       "_" + channel + extra_text + image_format))


def main():
    """Main function"""
    make_publication_plots()


main()
