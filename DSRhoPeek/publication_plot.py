#!/usr/bin/env python3
"""This script generates publication quality plots from ROOT files with
TTrees.
"""

import argparse
import glob
import json
import os
import ROOT


def setup_plot_style():
    """Setup a sane ROOT plot style"""
    ROOT.gROOT.ProcessLine(".L src/colors.cc")
    ROOT.gROOT.ProcessLine("colors::setColors()")

    ROOT.gStyle.SetLabelFont(43, "xyz")
    ROOT.gStyle.SetLabelSize(18, "xyz")
    ROOT.gStyle.SetLabelOffset(0.01, "xyz")
    ROOT.gStyle.SetTitleFont(43, "xyz")
    ROOT.gStyle.SetTitleSize(18, "xyz")
    ROOT.gStyle.SetTitleOffset(1.2)
    ROOT.gStyle.SetMarkerSize(0.5)
    # Disable perpendicular lines at the end of error bars
    ROOT.gStyle.SetEndErrorSize(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetPadTopMargin(0.05)
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadLeftMargin(0.105)
    ROOT.gStyle.SetPadBottomMargin(0.1)
    ROOT.gStyle.SetOptStat(0)


def create_legend(plot_data, canvas):
    """Create a legend based on labels in plot_data"""
    # Default legend placement is top right
    x_min = 0.75
    x_max = 0.85
    y_min = 0.9 - 0.04 * len(plot_data['elements'])
    y_max = 0.9

    if 'legendPosition' in plot_data:
        if plot_data['legendPosition'] == "top left":
            x_min = 0.15
            x_max = 0.25
            y_min = 0.9 - 0.04 * len(plot_data['elements'])
            y_max = 0.9

        elif plot_data['legendPosition'] == "top right":
            x_min = 0.75
            x_max = 0.85
            y_min = 0.9 - 0.04 * len(plot_data['elements'])
            y_max = 0.9

        elif plot_data['legendPosition'] == "bottom right":
            x_min = 0.75
            x_max = 0.85
            y_min = 0.12
            y_max = 0.12 + 0.04 * len(plot_data['elements'])

        elif plot_data['legendPosition'] == "bottom left":
            x_min = 0.15
            x_max = 0.25
            y_min = 0.12
            y_max = 0.12 + 0.04 * len(plot_data['elements'])

    legend = ROOT.TLegend(x_min, y_min, x_max, y_max)

    legend.SetBorderSize(0)
    legend.SetFillColor(0)
    legend.SetTextFont(43)
    legend.SetTextSize(16)

    # Color changes must be done after Draw() because of ROOT
    num_histos_found = 0
    for i in range(len(canvas.GetListOfPrimitives())):
        # Not all objects in the list are the ones we are looking for
        if canvas.GetListOfPrimitives().At(i).GetName() == 'htemp':
            # This allows plot elements without labels. Useful when plotting a
            # single element.
            if 'label' in plot_data['elements'][num_histos_found]:
                legend.AddEntry(canvas.GetListOfPrimitives().At(i),
                                plot_data['elements'][num_histos_found]['label'], "L")
            num_histos_found += 1

    return legend


def make_publication_plots(json_file, image_format, plot_dir, plots):
    """Make publication quality plots"""
    with open(json_file, 'r') as f:
        plots_data = json.loads(f.read())['plots']

        setup_plot_style()

        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)

        output_file = ROOT.TFile(os.path.join(
            plot_dir, "plots.root"), "RECREATE")

        for plot_data in plots_data:
            if plots and plot_data["fileName"] not in plots:
                continue
            make_plot(plot_data, image_format, plot_dir, output_file)

        output_file.Close()


def make_plot(plot_data, image_format, plot_dir, output_file):
    """Make a single plot, save it as a file and into the ROOT output file"""
    tree = ROOT.TChain(plot_data['treeName'])
    if 'files' in plot_data:
        tree = read_in_files(plot_data['files'], plot_data['treeName'])

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

        if 'files' not in element and 'files' not in plot_data:
            raise Exception("No 'files' array in '{}' plot or element data".format(
                plot_data['fileName']))

        draw_element(tree, element, plot_data, same_opt)

    histograms = get_histograms(canvas)

    # 2D plots don't need legends as there is only a single element and it
    # would look bad
    if ':' not in plot_data['formula']:
        change_line_colors(plot_data, histograms)
        set_optimal_histogram_yrange(histograms)
        legend = create_legend(plot_data, canvas)
        legend.Draw()

    histo = canvas.GetPrimitive("htemp")
    set_histogram_titles(histo, plot_data)

    output_file.cd()
    histo.Write()
    canvas.SetTickx()
    canvas.SetTicky()
    canvas.SaveAs(os.path.join(plot_dir, plot_data['fileName'] + image_format))
    canvas.IsA().Destructor(canvas)


def change_line_colors(plot_data, histograms):
    """Change line colors in a plot according to plot_data info"""
    # Color changes must be done after Draw() because of ROOT
    for i, histo in enumerate(histograms):
        histo.SetLineColor(plot_data['elements'][i]['color'])


def set_optimal_histogram_yrange(histograms):
    """Set an optimal histogram Y axis range"""
    max_value = max([histo.GetMaximum() for histo in histograms])
    for histo in histograms:
        histo.GetYaxis().SetRangeUser(0, max_value * 1.05)


def get_histograms(canvas):
    """Return a list of all current histograms (htemp(s))"""
    histograms = []
    for i in range(canvas.GetListOfPrimitives().GetEntries()):
        # There are other things than histos in the list
        if canvas.GetListOfPrimitives().At(i).GetName() == 'htemp':
            histograms.append(canvas.GetListOfPrimitives().At(i))

    return histograms


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
    # If this element has it's own 'file' key defined, it
    # superseeds the plot's 'file'
    if 'files' in element:
        tree = read_in_files(element['files'], plot_data['treeName'])

    if tree is None:
        print("WARNING: Skipping element '{}' in plot '{}'".format(
            element['label'], plot_data['fileName']))
        return

    cuts = get_cut_string(plot_data, element)

    formula = plot_data['formula']
    if 'formula' in element:
        formula = element['formula']

    num_passing = tree.Draw(formula, cuts, same_opt)
    if num_passing == 0:
        print("WARNING: In plot '" + plot_data['fileName'] +
              "' No events are passing the following cut: " + cuts)


def get_cut_string(plot_data, element):
    """Return a string with all the relevant cuts"""
    cuts = []
    if 'cut' in plot_data:
        cuts.append(plot_data['cut'])
    if 'cut' in element:
        cuts.append(element['cut'])
    if 'min' in plot_data:
        assert ':' not in plot_data['formula'], "'min' can't be used in 2D plots; use 'xmin' or 'ymin'"
        cuts.append(plot_data['formula'] + '>' + str(plot_data['min']))
    if 'max' in plot_data:
        assert ':' not in plot_data['formula'], "'max' can't be used in 2D plots; use 'xmax' or 'ymax'"
        cuts.append(plot_data['formula'] + '<' + str(plot_data['max']))
    if 'xmin' in plot_data:
        assert ':' in plot_data['formula'], "'xmin' can't be used in 1D plots; use 'min'"
        cuts.append(plot_data['formula'].split(':')[1] +
                    '>' + str(plot_data['xmin']))
    if 'xmax' in plot_data:
        assert ':' in plot_data['formula'], "'xmax' can't be used in 1D plots; use 'max'"
        cuts.append(plot_data['formula'].split(':')[1] +
                    '<' + str(plot_data['xmax']))
    if 'ymin' in plot_data:
        assert ':' in plot_data['formula'], "'ymin' can't be used in 1D plots"
        cuts.append(plot_data['formula'].split(':')[0] +
                    '>' + str(plot_data['ymin']))
    if 'ymax' in plot_data:
        assert ':' in plot_data['formula'], "'ymax' can't be used in 1D plots"
        cuts.append(plot_data['formula'].split(':')[0] +
                    '<' + str(plot_data['ymax']))

    return "&&".join(cuts)


def read_in_files(file_entries, tree_name):
    """Read trees from files, while unglobbing the filenames."""
    files = []
    for entry in file_entries:
        files += glob.glob(entry)

    if files:
        tree = ROOT.TChain(tree_name)
        for root_file in files:
            tree.Add(root_file)
        return tree

    print("WARNING: None of the following entries expand to existing files")
    for entry in file_entries:
        print(entry)
    return None


def decode_arguments():
    """Decode CLI arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("json_file")
    parser.add_argument("-f", "--format", type=str,
                        help="image format to generate (.pdf, .png, ..., defaults to '.png')")
    parser.add_argument("-d", "--dir", type=str,
                        help="directory for generated plots (defaults to 'plots')")
    parser.add_argument("-p", "--plots", type=str, nargs='+')
    args = parser.parse_args()

    if not args.format:
        args.format = '.png'
    if not args.dir:
        args.dir = 'plots'

    return args.json_file, args.format, args.dir, args.plots


def main():
    """Main function"""
    json_file, image_format, plot_dir, plots = decode_arguments()
    print(json_file, image_format, plot_dir, plots)
    make_publication_plots(json_file, image_format, plot_dir, plots)


main()
