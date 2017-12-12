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

        output_file = ROOT.TFile(os.path.join(
            plot_dir, "plots.root"), "RECREATE")

        for plot_data in plots_data:
            make_plot(plot_data, image_format, plot_dir, output_file)

        output_file.Close()


def make_plot(plot_data, image_format, plot_dir, output_file):
    """Make a single plot, save it as a file and into the ROOT output file"""
    root_file = ROOT.TFile(plot_data['file'])
    tree = root_file.Get(plot_data['treeName'])

    canvas_width = 500
    canvas_height = 500
    canvas = ROOT.TCanvas("canvas", "canvas", canvas_width, canvas_height)
    canvas.cd()

    # 2D plots must have a larger margin to accomodate z legend
    if ':' in plot_data['formula']:
        canvas.SetRightMargin(0.14)

    for i, element in enumerate(plot_data['elements']):
        # First element in a plot must be drawn without "same" or all
        # the plots would get piled together
        same_opt = "same colz"
        if i == 0:
            same_opt = "colz"

        draw_element(tree, element, plot_data, same_opt)

    # 2D plots don't need legends as there is only a single element and it
    # would look bad
    if ':' not in plot_data['formula']:
        change_line_colors(plot_data)
        legend = create_legend(plot_data)
        legend.Draw()

    histo = ROOT.gPad.GetPrimitive("htemp")
    set_histogram_titles(histo, plot_data)

    output_file.cd()
    histo.Write()
    canvas.SetTickx()
    canvas.SetTicky()
    canvas.SaveAs(os.path.join(plot_dir, plot_data['fileName'] + image_format))
    canvas.IsA().Destructor(canvas)


def change_line_colors(plot_data):
    """Change line colors in a plot according to plot_data info"""
    # Color changes must be done after Draw() because of ROOT
    num_histos_found = 0
    for i in range(ROOT.gPad.GetListOfPrimitives().GetEntries()):
        # There are other things than histos in the list
        if ROOT.gPad.GetListOfPrimitives().At(i).GetName() == 'htemp':
            ROOT.gPad.GetListOfPrimitives().At(i).SetLineColor(
                plot_data['elements'][num_histos_found]['color'])
            num_histos_found += 1


def set_histogram_titles(histo, plot_data):
    """Set various titles in the histogram according to plot_data"""
    histo.SetName(plot_data['fileName'])
    histo.SetTitle("")
    if 'xAxisTitle' in plot_data:
        histo.GetXaxis().SetTitle(plot_data['xAxisTitle'])
    if 'yAxisTitle' in plot_data:
        histo.GetYaxis().SetTitle(plot_data['yAxisTitle'])


def draw_element(tree, element, plot_data, same_opt):
    """Draw a single element of a plot, e.g., a histogram"""
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
