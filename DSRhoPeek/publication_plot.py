#!/usr/bin/env python
"""This script generates publication quality plots from ROOT files with
TTrees.
"""

import argparse
import json
import os
import ROOT


def setup_plot_style():
    """Setup a sane ROOT plot style"""
    ROOT.gROOT.Reset()
    ROOT.gROOT.ProcessLine(".L src/colors.cc")
    ROOT.gROOT.ProcessLine("colors::setColors()")
    ROOT.gROOT.ProcessLine(".L src/tools.cc")
    ROOT.gROOT.ProcessLine("tools::SetupPlotStyle()")
    ROOT.gStyle.SetOptStat(0)


def create_legend(plot_data):
    """Create a legend based on labels in plot_data"""
    legend = ROOT.TLegend(
        0.75, 0.9 - 0.04 * len(plot_data['elements']), 0.85, 0.9)
    #mylegend = TLegend(0.15,0.9 - 0.04*num_hists_to_draw,0.25,0.9)
    legend.SetBorderSize(0)
    legend.SetFillColor(0)
    legend.SetTextFont(43)
    legend.SetTextSize(16)

    # Color changes must be done after Draw() because of ROOT
    for i in range(len(plot_data['elements'])):
        legend.AddEntry(ROOT.gPad.GetListOfPrimitives().At(i),
                        plot_data['elements'][i]['label'], "L")

    return legend


def make_publication_plots(json_file, image_format, plot_dir):
    """Make publication quality plots"""
    with open(json_file, 'r') as f:
        plots_data = json.loads(f.read())['plots']

        setup_plot_style()
        canvas_width = 500
        canvas_height = 500

        output_file = ROOT.TFile(os.path.join(
            plot_dir, "plots.root"), "RECREATE")
        canvas = ROOT.TCanvas("canvas", "canvas", canvas_width, canvas_height)
        canvas.cd()

        for plot_data in plots_data:
            root_file = ROOT.TFile(plot_data['file'])
            tree = root_file.Get(plot_data['treeName'])

            for i, element in enumerate(plot_data['elements']):
                # First element in a plot must be drawn without "same" or all
                # the plots would get piled together
                same_opt = "same"
                if i == 0:
                    same_opt = ""

                if 'cut' in element:
                    cut = element['cut']
                else:
                    # A dummy cut so we can always use && when adding something
                    cut = "1==1"

                if 'min' in plot_data:
                    cut += '&&' + plot_data['formula'] + \
                        '>' + str(plot_data['min'])
                if 'max' in plot_data:
                    cut += '&&' + plot_data['formula'] + \
                        '<' + str(plot_data['max'])

                # If this element has it's own 'file' key defined, it
                # superseeds the plot's 'file'
                if 'file' in element:
                    root_file_temp = ROOT.TFile(element['file'])
                    tree_temp = root_file_temp.Get(plot_data['treeName'])
                    tree_temp.Draw(plot_data['formula'], cut, same_opt)
                else:
                    tree.Draw(plot_data['formula'], cut, same_opt)

            # Color changes must be done after Draw() because of ROOT
            for i in range(len(plot_data['elements'])):
                ROOT.gPad.GetListOfPrimitives().At(i).SetLineColor(
                    plot_data['elements'][i]['color'])

            histo = ROOT.gPad.GetPrimitive("htemp")
            histo.SetName(plot_data['fileName'])
            histo.GetXaxis().SetTitle(plot_data['xAxisTitle'])
            histo.SetTitle("")

            legend = create_legend(plot_data)
            legend.Draw()

            output_file.cd()
            histo.Write()
            canvas.SetTickx()
            canvas.SetTicky()
            canvas.SaveAs(
                os.path.join(plot_dir, plot_data['fileName'] + image_format))

        output_file.Close()


def decode_arguments():
    """Decode CLI arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("json_file")
    parser.add_argument("-f", "--format", type=str,
                        help="image format to generate (.pdf, .png, ..., defaults to '.png')")
    parser.add_argument("-d", "--dir", type=str,
                        help="directory for generated plots (defaults to 'plots')")
    args = parser.parse_args()

    if not args.format:
        args.format = '.png'
    if not args.dir:
        args.dir = 'plots'

    return args.json_file, args.format, args.dir


def main():
    """Main function"""
    json_file, image_format, plot_dir = decode_arguments()
    print(json_file, image_format, plot_dir)
    make_publication_plots(json_file, image_format, plot_dir)


main()
