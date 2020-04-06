/**
 *  @file    tools.cc
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2015-10-07
 *
 *  @brief A collections of various tools and utilities
 *
 */
#include "tools.h"

// Standard includes
#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>

// Boost includes
#include <boost/filesystem.hpp>

// ROOT includes
#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooHist.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TEnv.h"
#include "TFile.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TList.h"
#include "TMath.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"

// Local includes
#include "cksum.h"
#include "constants.h"
#include "gitversion.h"
#include "log.h"
#include "nlohmann/json.hpp"

namespace tools {

/**
 * Return a list of files with a requested extension in a directory.
 *
 * @param dir Directory where to look for files
 * @param ext Files with which extension to list
 */
std::vector<TString> GetListOfFiles(const char* dir, const char* ext) {
    std::vector<TString> file_names;
    TSystemDirectory system_dir(dir, dir);
    TList* files = system_dir.GetListOfFiles();
    if (files) {
        TSystemFile* file;
        TString file_name;
        TIter next(files);
        while ((file = static_cast<TSystemFile*>(next()))) {
            file_name = file->GetName();
            if (!file->IsDirectory() && file_name.EndsWith(ext)) {
                file_names.push_back(file_name);
            }
        }
    }
    delete files;
    return file_names;
}

/**
 * Extract certain trees from data files passing certain conditions,
 * collect them in a TChain and then return a pointer to the TChain.
 *
 * @param dir Directory from which to read in data files
 */
TChain* ReadDataFromDir(const char* dir) {
    TChain* chain = new TChain("h2000");
    std::vector<TString> list_of_files = GetListOfFiles(dir, ".root");
    TString directory = dir;
    // Make the function agnostic to possible '/' at the end of the dir argument
    if (!directory.EndsWith("/")) {
        directory += "/";
    }
    for (auto file_name : list_of_files) {
        if (file_name.Contains("results")) continue;  // Skip results files created by DSRhoPeek
        chain->Add(directory + file_name);
    }
    return chain;
}

/**
 * Setup a sane ROOT plot style.
 */
void SetupPlotStyle() {
    // Fonts ending in 3 are not proportional in size to the pad; their size is in pixels
    // The "xyz" parameter dictates that these settings should be used for all axes
    gStyle->SetLabelFont(43, "xyz");
    gStyle->SetLabelSize(18, "xyz");
    gStyle->SetLabelOffset(0.01, "xyz");
    gStyle->SetTitleFont(43, "xyz");
    gStyle->SetTitleSize(18, "xyz");
    gStyle->SetTitleOffset(1.2);
    gStyle->SetMarkerSize(0.5);
    gStyle->SetEndErrorSize(0);  // Disable perpendicular lines at the end of error bars
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadLeftMargin(0.105);
    gStyle->SetPadBottomMargin(0.1);
    gStyle->SetOptStat(0);
}

/**
 * Set the directory to which to ouput ROOT plots
 *
 * @param plot_dir Directory to which to save plots
 */
void SetPlotDir(const char* plot_dir) {
    if (!boost::filesystem::is_directory(plot_dir)) {
        boost::filesystem::create_directories(plot_dir);
    }
    gEnv->SetValue("Canvas.PrintDirectory", plot_dir);
    Log::print(Log::info, "Setting plot directory: %s\n",
               gEnv->GetValue("Canvas.PrintDirectory", "not found"));
}

/**
 * Create a stat box with supplied fit results.
 *
 * @param chi2 A chi2 value to be included.
 * @param results A list of fit results.
 * @param position_top Place the box at the top (or bottom).
 * @param position_left Place the box to the left (or right).
 */
TPaveText* CreateStatBox(double chi2, int ndof, const RooArgList& results, bool position_top,
                         bool position_left) {
    double x_left, x_right, y_bottom, y_top;
    const double line_height = 0.06;

    if (position_top) {
        y_top = 0.9;
        y_bottom = y_top - (results.getSize() + 2) * line_height;
    } else {
        y_bottom = 0.023;
        y_top = y_bottom + (results.getSize() + 2) * line_height;
    }

    if (position_left) {
        x_left = 0.3;
        x_right = 0.3;
    } else {
        x_left = 0.7;
        x_right = 0.7;
    }

    TPaveText* stat_box = new TPaveText(x_left, y_bottom, x_right, y_top, "NDC");
    stat_box->SetShadowColor(kWhite);
    stat_box->SetBorderSize(0);
    stat_box->SetFillColor(kWhite);
    stat_box->SetTextFont(43);
    stat_box->SetTextSize(14);
    stat_box->SetY1NDC(0.1);

    char line[1000];
    for (int i = 0; i < results.getSize(); i++) {
        snprintf(line, 1000, "%s = %.3f +- %.3f", results[i].GetTitle(),
                 static_cast<RooRealVar*>(results.at(i))->getVal(),
                 static_cast<RooRealVar*>(results.at(i))->getError());
        stat_box->AddText(line);
    }
    snprintf(line, 1000, "#chi^{2}/ndof = %.2f (%.1f/%i)\n", chi2 / ndof, chi2, ndof);
    stat_box->AddText(line);
    snprintf(line, 1000, "p = %.2f\n", TMath::Prob(chi2, ndof));
    stat_box->AddText(line);
    return stat_box;
}

/**
 * Make a simple plot of a variable and save it to a file.
 *
 * NOTE: @var can't be const (without const_cast) because of dumb RooFit design
 */
void PlotVar(const RooRealVar& var, const RooAbsData& data) {
    TCanvas* canvas = new TCanvas(var.GetName(), var.GetTitle(), 500, 500);

    RooPlot* plot = var.frame();

    data.plotOn(plot);
    plot->Draw();

    plot->SetTitle("");
    plot->GetXaxis()->SetTitle(TString(var.GetTitle()));
    plot->GetYaxis()->SetTitle("");

    canvas->Write();
    canvas->SaveAs(constants::format);
}

/**
 * Make a simple plot of a variable and save it to a file.
 *
 * NOTE: @var can't be const (without const_cast) because of dumb RooFit design
 *
 * @param var Variable to be plotted
 * @param pdf PDF to be plotted
 * @param norm_vars Variables which are to be projected out, not used as slices
 */
void PlotVar(const RooRealVar& var, const RooAbsPdf& pdf, const RooArgSet& norm_vars) {
    TString name = pdf.GetName();
    name += "_";
    name += var.GetName();
    TCanvas canvas(name, name, 500, 500);

    RooPlot* plot = var.frame();

    // From RooFit manual: "No variables are projected by default when PDF is
    // plotted on an empty frame" One has to be careful to explicitly project
    // over intended variables in this case, otherwise they are only slices.
    plot->updateNormVars(norm_vars);

    pdf.plotOn(plot);
    plot->Draw();

    plot->SetTitle("");
    plot->GetXaxis()->SetTitle(TString(var.GetTitle()));
    plot->GetYaxis()->SetTitle("");

    canvas.Write();
    canvas.SaveAs(constants::format);
}

/**
 * Create a plot of two RooDataHists with an optional pull or residual subplot
 *
 * @param var Variable to be plotted
 * @param data1 First histogram to be plotted
 * @param data2 Second histogram to be plotted
 * @param draw_pull Whether a pull subplot should be created
 * @param draw_residual Whether a residual subplot should be created
 */
void PlotVar(const RooRealVar& var, const RooDataHist& data1, const RooDataHist& data2,
             bool draw_pull, bool draw_residual) {
    if (draw_pull && draw_residual) {
        Log::print(Log::error,
                   "Requested to plot var with both pull and residual. Select only one!\n");
        return;
    }

    TCanvas canvas(
        TString(var.GetName()) + "_" + TString(data1.GetName()) + "_" + TString(data2.GetName()),
        TString(var.GetTitle()) + " canvas", 500, 500);

    RooPlot* plot = var.frame();

    TPad* pad_var;
    TPad* pad_pull;
    if (draw_pull || draw_residual) {
        pad_var = new TPad("pad_var", "pad_var", 0, 0.25, 1, 1);
        pad_var->SetBottomMargin(0.0);
        pad_pull = new TPad("pad_pull", "pad_pull", 0, 0, 1, 0.25);
        pad_pull->Draw();
        plot->GetXaxis()->SetTitle("");
        plot->GetXaxis()->SetLabelSize(0);
        // This line makes sure the 0 is not drawn as it would overlap with the lower
        // pad
        plot->GetYaxis()->SetRangeUser(0.001, plot->GetMaximum());
    } else {
        pad_var = new TPad("pad_var", "pad_var", 0, 0, 1, 1);
    }
    pad_var->Draw();
    pad_var->cd();
    pad_var->SetLeftMargin(0.12);

    data1.plotOn(plot, RooFit::DataError(RooAbsData::None), RooFit::Name("data1"));
    data2.plotOn(plot, RooFit::MarkerColor(2), RooFit::LineColor(2),
                 RooFit::DataError(RooAbsData::None), RooFit::Name("data2"));

    plot->SetTitle("");
    plot->GetYaxis()->SetTitle("");
    plot->GetYaxis()->SetTitleOffset(1.60);
    plot->Draw();

    TLegend* leg1 = new TLegend(0.65, 0.73, 0.86, 0.87);
    leg1->SetFillColor(kWhite);
    leg1->SetLineColor(kWhite);
    leg1->AddEntry(plot->findObject("data1"), data1.GetTitle(), "LP");
    leg1->AddEntry(plot->findObject("data2"), data2.GetTitle(), "LP");
    leg1->Draw();

    if (draw_pull || draw_residual) {
        pad_pull->cd();
        pad_pull->SetTopMargin(0.0);
        pad_pull->SetBottomMargin(0.35);
        pad_pull->SetLeftMargin(0.12);

        // Create a new frame to draw the pull distribution and add the distribution
        // to the frame
        RooPlot* plot_pull = var.frame(RooFit::Title("Pull Distribution"));
        plot_pull->SetTitle("");

        TH1* data1_histo = data1.createHistogram("data1", var);
        TH1* data2_histo = data2.createHistogram("data2", var);
        TH1* pull_histo = dynamic_cast<TH1*>(data1_histo->Clone("pull_histo"));

        double p1;
        double p2;
        const int questionable_limit = 10;
        int questionable_bins = 0;
        for (int i = 1; i <= pull_histo->GetNbinsX(); i++) {
            p1 = data1_histo->GetBinContent(i);
            p2 = data2_histo->GetBinContent(i);
            if (draw_residual) {
                pull_histo->SetBinContent(i, p1 - p2);
            } else if (draw_pull) {
                pull_histo->SetBinContent(i, (p1 - p2) / std::sqrt(p2));
                if (p2 < questionable_limit) {
                    questionable_bins++;
                }
            }
        }

        if (questionable_bins) {
            const int total_bins = pull_histo->GetNbinsX();
            Log::print(
                Log::warning,
                "There were %i/%i (%.0f%%) bins with less than %i events - Poisson approximation "
                "questionable! Consider using residuals.\n",
                questionable_bins, total_bins, 100 * (double)questionable_bins / total_bins,
                questionable_limit);
        }

        RooHist* hpull = new RooHist(*pull_histo, 0, 1, RooAbsData::ErrorType::None);

        hpull->SetFillColor(kGray);
        // The only working way to get rid of error bars; HIST draw option doesn't
        // work with RooPlot
        for (int i = 0; i < hpull->GetN(); i++) {
            hpull->SetPointError(i, 0, 0, 0, 0);
        }
        plot_pull->addPlotable(hpull, "B");
        // We plot again without bars, so the points are not half covered by
        // bars as in case of "BP" draw option. We need to create and plot a
        // clone, because the ROOT object ownership is transfered to the RooPlot
        // by addPlotable().
        // If we just added the hpull twice, we would get a segfault.
        RooHist* hpull_clone = dynamic_cast<RooHist*>(hpull->Clone());
        // We plot again without bars, so the points are not half covered by bars
        plot_pull->addPlotable(hpull_clone, "P");

        plot_pull->GetXaxis()->SetTickLength(0.03 * pad_var->GetAbsHNDC() / pad_pull->GetAbsHNDC());
        plot_pull->GetXaxis()->SetTitle(TString(var.GetTitle()));
        plot_pull->GetXaxis()->SetTitleOffset(4.0);
        plot_pull->GetXaxis()->SetLabelOffset(0.01 * pad_var->GetAbsHNDC() /
                                              pad_pull->GetAbsHNDC());
        plot_pull->GetYaxis()->SetTitle("");
        // plot_pull->GetYaxis()->SetRangeUser(0.8, 1.19);
        double max = std::max(-pull_histo->GetMinimum(), pull_histo->GetMaximum());
        plot_pull->GetYaxis()->SetRangeUser(-max, +max);
        // plot_pull->GetYaxis()->SetRangeUser(-5, 5);
        plot_pull->GetYaxis()->SetNdivisions(505);
        plot_pull->Draw();

        delete data1_histo;
        delete data2_histo;
        delete pull_histo;
    }

    canvas.Write();
    canvas.SaveAs(constants::format);

    delete plot;

    delete pad_var;
    if (draw_pull) {
        delete pad_pull;
    }
}

/**
 * Make a plot of a variable with both data points and pdf projection and display a pull plot
 * beneath it. Then save to a file.
 *
 * @param var Variable to be plotted
 * @param data Dataset against which to plot
 * @param pdf PDF to use for plotting and pull calculation
 * @param components [optional] Component PDFs to plot alongside the full PDF
 * @param title [optional] Title for the y-axis
 */
void PlotWithPull(const RooRealVar& var, const RooArgSet& projection_vars, const RooDataSet& data,
                  const RooAbsPdf& pdf, const RooFitResult* const result,
                  const std::vector<RooAbsPdf*> components, int numCPUs, std::string prefix,
                  std::string title, std::vector<RooCmdArg> options) {
    TString name;
    if (prefix.size()) {
        name += prefix.c_str();
        name += "_";
    }
    name += pdf.GetName();
    name += "_";
    name += data.GetName();
    name += "_";
    name += var.GetName();
    TCanvas canvas(name, name, 500, 500);

    RooPlot* plot = var.frame();
    TPad* pad_var;
    TPad* pad_pull;
    pad_var = new TPad("pad_var", "pad_var", 0, 0.25, 1, 1);
    pad_pull = new TPad("pad_pull", "pad_pull", 0, 0, 1, 0.25);
    pad_var->Draw();
    pad_pull->Draw();

    pad_var->cd();
    pad_var->SetBottomMargin(0.0);
    pad_var->SetLeftMargin(0.12);

    data.plotOn(plot);

    options.push_back(RooFit::ProjWData(projection_vars, data, kFALSE));
    options.push_back(RooFit::NumCPU(numCPUs));

    if (components.size() > 1) {
        // Plot components before the total PDF as the pull plots are made from the
        // last plotted PDF
        const int colors[] = {4, 7, 5};
        const int styles[] = {3345, 3354, 3395};
        int i = 0;
        for (auto component : components) {
            std::vector<RooCmdArg> component_options;
            component_options.push_back(RooFit::LineColor(colors[i]));
            component_options.push_back(RooFit::FillColor(colors[i]));
            component_options.push_back(RooFit::FillStyle(styles[i]));
            component_options.push_back(RooFit::DrawOption("FL"));
            component_options.push_back(RooFit::VLines());
            RooArgSet set(*component);
            component_options.push_back(RooFit::Components(set));

            std::vector<RooCmdArg> all_options = AddVectors(options, component_options);
            RooLinkedList roo_options = VecToCmdList(all_options);

            pdf.plotOn(plot, roo_options);
            i++;
        }
    }

    RooLinkedList roo_options = VecToCmdList(options);
    pdf.plotOn(plot, roo_options);

    plot->GetXaxis()->SetTitle("");
    plot->GetXaxis()->SetLabelSize(0);

    // This line makes sure the 0 is not drawn as it would overlap with the lower pad
    plot->GetYaxis()->SetRangeUser(0.001, plot->GetMaximum());
    plot->SetTitle("");
    plot->GetYaxis()->SetTitle(title.c_str());
    plot->GetYaxis()->SetTitleOffset(1.60);
    plot->Draw();

    // This is not an OK way of calculating ndof - e.g., in a case where the PDF
    // is defined only on a subset of the var range, or there are empty bins,
    // etc. it will be wrong. However, RooFit doesn't give access to anything
    // like ndof itself, so this will have to do for now.
    const int num_floating_pars = result ? result->floatParsFinal().getSize() : 0;
    const int ndof = var.getBinning().numBins() - num_floating_pars;
    const double chi2 = plot->chiSquare(num_floating_pars) * ndof;
    TPaveText* stat_box = CreateStatBox(chi2, ndof, result->floatParsFinal(), true, true);
    if (stat_box) {
        stat_box->Draw();
    }

    pad_pull->cd();
    pad_pull->SetTopMargin(0.0);
    pad_pull->SetBottomMargin(0.35);
    pad_pull->SetLeftMargin(0.12);

    // Create a new frame to draw the pull distribution and add the distribution to the frame
    RooPlot* plot_pull_ = var.frame(RooFit::Title("Pull Distribution"));
    plot_pull_->SetTitle("");
    RooHist* hpull = plot->pullHist();
    hpull->SetFillColor(kGray);
    // The only working way to get rid of error bars; HIST draw option doesn't work with RooPlot
    for (int i = 0; i < hpull->GetN(); i++) {
        hpull->SetPointError(i, 0, 0, 0, 0);
    }
    plot_pull_->addPlotable(hpull, "B");
    // We plot again without bars, so the points are not half covered by bars as in case of "BP"
    // draw option.
    // We need to create and plot a clone, because the ROOT object ownership is transfered to the
    // RooPlot by addPlotable().
    // If we just added the hpull twice, we would get a segfault.
    RooHist* hpull_clone = dynamic_cast<RooHist*>(hpull->Clone());

    // We plot again without bars, so the points are not half covered by bars
    plot_pull_->addPlotable(hpull_clone, "P");

    plot_pull_->GetXaxis()->SetTickLength(0.03 * pad_var->GetAbsHNDC() / pad_pull->GetAbsHNDC());
    plot_pull_->GetXaxis()->SetTitle(TString(var.GetTitle()));
    plot_pull_->GetXaxis()->SetTitleOffset(4.0);
    plot_pull_->GetXaxis()->SetLabelOffset(0.01 * pad_var->GetAbsHNDC() / pad_pull->GetAbsHNDC());
    plot_pull_->GetYaxis()->SetRangeUser(-5, 5);
    plot_pull_->GetYaxis()->SetNdivisions(505);
    plot_pull_->Draw();

    canvas.Write();
    canvas.SaveAs(constants::format);
}

/**
 * Convert a vector of RooFit options into a list that plotOn understands
 *
 * @param commands A vector of commands such as RooFit::NumCPU(4)
 * @return RooLinkedList The resulting command list that RooFit understands
 */
RooLinkedList VecToCmdList(std::vector<RooCmdArg>& commands) {
    RooLinkedList cmd_list;
    for (auto& cmd : commands) {
        cmd_list.Add(&cmd);
    }
    return cmd_list;
}

/*
 * Create and save 2D plots of supplied vars and data
 *
 * @param var1 First variable.
 * @param var2 Second variable.
 * @param data Data to be plotted.
 * @param prefix [optional] Prefix to be added to filenames.
 * @param format [optional] Format in which to save the images.
 * @param max [optional] Maximum of the range
 */
void PlotVars2D(const RooRealVar& var1, const RooRealVar& var2, const RooAbsData& data,
                const std::string prefix, const char* format, const double max) {
    TString name;
    if (prefix.length()) {
        name = prefix + "_";
    }
    name += data.GetName();
    name += "_";
    name += var1.GetName();
    name += "_";
    name += var2.GetName();

    TCanvas canvas(name, name, 500, 500);
    TH2D* histo = static_cast<TH2D*>(data.createHistogram("histo", var1, RooFit::YVar(var2)));

    canvas.SetRightMargin(0.14);

    histo->SetMinimum(0);
    if (max != 0) histo->SetMaximum(max);
    histo->SetTitle("");
    histo->Draw("colz");
    histo->GetZaxis()->SetTitle("");

    canvas.Write();
    canvas.SaveAs(format);
    delete histo;
}

/*
 * Create and save 2D pulls of supplied vars and two datahists
 *
 * @param var1 First variable.
 * @param var2 Second variable.
 * @param data Data from which to calculate pull
 * @param pdf Histogrammed PDF from which to calculate pull
 * @param format [optional] Format in which to save the images.
 * @param residual [optional] Plot a residual instead of a pull
 */
void PlotPull2D(const RooRealVar& var1, const RooRealVar& var2, const RooAbsData& data,
                const RooAbsData& pdf, const std::string prefix, const char* format,
                const bool residual) {
    TString type = residual ? "residual" : "pull";
    TString name;
    if (prefix.length()) {
        name = prefix + "_";
    }
    name += data.GetName();
    name += "_";
    name += var1.GetName();
    name += "_";
    name += var2.GetName();
    name += "_";
    name += type;

    TCanvas canvas(name, name, 500, 500);
    gStyle->SetPalette(kLightTemperature);

    TH2D* histo1 = static_cast<TH2D*>(data.createHistogram("histo1", var1, RooFit::YVar(var2)));
    TH2D* histo2 = static_cast<TH2D*>(pdf.createHistogram("histo2", var1, RooFit::YVar(var2)));
    TH2D* pull_histo =
        static_cast<TH2D*>(pdf.createHistogram("pull_histo", var1, RooFit::YVar(var2)));

    double p1;
    double p2;
    const int questionable_limit = 10;
    int questionable_bins = 0;
    for (int i = 1; i <= pull_histo->GetNbinsX(); i++) {
        for (int j = 1; j <= pull_histo->GetNbinsY(); j++) {
            p1 = histo1->GetBinContent(i, j);
            p2 = histo2->GetBinContent(i, j);
            // pull_histo->SetBinContent(i, j, p1/p2);
            if (residual) {
                pull_histo->SetBinContent(i, j, p1 - p2);
            } else {
                pull_histo->SetBinContent(i, j, (p1 - p2) / std::sqrt(p2));
                if (p2 < questionable_limit) {
                    questionable_bins++;
                }
            }
        }
    }

    if (questionable_bins) {
        const int total_bins = pull_histo->GetNbinsX() * pull_histo->GetNbinsY();
        Log::print(
            Log::warning,
            "There were %i/%i (%.0f%%) bins with less than %i events - Poisson approximation "
            "questionable! Consider using residuals.\n",
            questionable_bins, total_bins, 100 * (double)questionable_bins / total_bins,
            questionable_limit);
    }

    const double pull_max = std::max(abs(pull_histo->GetMinimum()), abs(pull_histo->GetMaximum()));
    pull_histo->SetMinimum(-pull_max);
    pull_histo->SetMaximum(pull_max);
    canvas.SetRightMargin(0.14);
    pull_histo->SetTitle("");
    pull_histo->Draw("colz");
    pull_histo->GetZaxis()->SetTitle("");

    canvas.Write();
    canvas.SaveAs(format);

    delete histo1;
    delete histo2;
    delete pull_histo;

    gStyle->SetPalette(kViridis);
}

/*
 * Create and save two 2D plots with the SAME range
 *
 * @param var1 First variable.
 * @param var2 Second variable.
 * @param data1 Data to be plotted.
 * @param data2 Data to be plotted.
 * @param format [optional] Format in which to save the images.
 */
void PlotVars2D(const RooRealVar& var1, const RooRealVar& var2, const RooAbsData& data1,
                const RooAbsData& data2, const std::string prefix, const char* format) {
    double max = 0;
    TH2D* histo1 = static_cast<TH2D*>(data1.createHistogram("histo1", var1, RooFit::YVar(var2)));
    TH2D* histo2 = static_cast<TH2D*>(data2.createHistogram("histo2", var1, RooFit::YVar(var2)));
    max = std::max(histo1->GetMaximum(), histo2->GetMaximum());
    delete histo1;
    delete histo2;
    PlotVars2D(var1, var2, data1, prefix, format, max);
    PlotVars2D(var1, var2, data2, prefix, format, max);
}

/**
 * Construct and return a cut string common for all for categories
 */
TString GetCommonCutsString() {
    TString common_cuts("vrusable==1&&vtusable==1&&");
    common_cuts += "((vrchi2/vrndf)<";
    common_cuts += constants::cuts::sig_vtx_h;
    common_cuts += "||vrntrk==1)&&";
    common_cuts += "((vtchi2/vtndf)<";
    common_cuts += constants::cuts::tag_vtx_h;
    common_cuts += "||vtntrk==1)&&";
    common_cuts += "((sqrt(vrerr6)<";
    common_cuts += constants::cuts::sig_vtx_multitrack_sigma_z;
    common_cuts += "&&vrntrk>1)||(sqrt(vrerr6)<";
    common_cuts += constants::cuts::sig_vtx_singletrack_sigma_z;
    common_cuts += "&&vrntrk==1))";
    common_cuts += "&&";
    common_cuts += "((sqrt(vterr6)<";
    common_cuts += constants::cuts::tag_vtx_multitrack_sigma_z;
    common_cuts += "&&vtntrk>1)||(sqrt(vterr6)<";
    common_cuts += constants::cuts::tag_vtx_singletrack_sigma_z;
    common_cuts += "&&vtntrk==1))&&";
    common_cuts += "csbdtg>";
    common_cuts += constants::cuts::cs_bdtg;
    common_cuts += "&&(de>";
    common_cuts += constants::cuts::de_low;
    common_cuts += "&&de<";
    common_cuts += constants::cuts::de_high;
    common_cuts += ")&&(dt>";
    common_cuts += constants::cuts::dt_low;
    common_cuts += "&&dt<";
    common_cuts += constants::cuts::dt_high;
    common_cuts += ")&&(phit>";
    common_cuts += constants::cuts::phit_low;
    common_cuts += "&&phit<";
    common_cuts += constants::cuts::phit_high;
    common_cuts += ")&&(thetab>";
    common_cuts += constants::cuts::thetab_low;
    common_cuts += "&&thetab<";
    common_cuts += constants::cuts::thetab_high;
    common_cuts += ")&&(thetat>";
    common_cuts += constants::cuts::thetat_low;
    common_cuts += "&&thetat<";
    common_cuts += constants::cuts::thetat_high;
    common_cuts += ")";
    return common_cuts;
}

/**
 * Save text to a file
 *
 * @param filename Path to the file to be written
 * @param text Text to be written to the file
 */
void SaveTextToFile(const std::string filename, const std::string text) {
    CreateDirsIfNecessary(filename);
    std::ofstream file(filename);
    file << text;
    file.close();
}

void LogTextFromFile(TFile* file, const char* field_name, const char* filename) {
    file->cd();
    std::ifstream input_file;
    input_file.open(filename);
    std::stringstream buffer;
    buffer << input_file.rdbuf();
    TNamed text(field_name, buffer.str());
    text.Write();
}

void LogFileCRC(TFile* file, const char* field_name, const char* filename) {
    file->cd();
    char buffer[100];
    snprintf(buffer, 100, "%lu", cksum(filename));
    TNamed crc(field_name, buffer);
    crc.Write();
}

void LogText(TFile* file, const char* field_name, const char* text) {
    file->cd();
    TNamed text_field(field_name, text);
    text_field.Write();
}

void LogText(TFile* file, const char* field_name, const std::string text) {
    LogText(file, field_name, text.c_str());
}

void LogCLIArguments(TFile* file, int argc, char* argv[]) {
    file->cd();
    std::string str;
    for (int i = 0; i < argc; i++) {
        str += argv[i];
        str += " ";
    }
    // Remove the final space
    str.pop_back();
    TNamed cli_arguments("cli_arguments", str.c_str());
    cli_arguments.Write();
}

void LogEnvironmentMetadata(TFile* file) {
    file->cd();
    TNamed root_version("root_version", ROOT_RELEASE);
    root_version.Write();

    char buffer[100];
    gethostname(buffer, 100);
    TNamed hostname("hostname", buffer);
    hostname.Write();

    const time_t now = time(0);
    const char* local_time_string = ctime(&now);
    TNamed local_date("local_date", local_time_string);

    tm* gmtm = gmtime(&now);
    const char* utc_time_string = asctime(gmtm);
    TNamed utc_date("utc_date", utc_time_string);

    TNamed git_version("git_version", gitversion);
    git_version.Write();
}

/**
 * Split string into multiple strings based on a delimiter
 */
std::vector<std::string> SplitString(const std::string& input_string, char delimiter) {
    std::vector<std::string> substrings;
    std::stringstream ss(input_string);
    std::string substring;
    while (std::getline(ss, substring, delimiter)) {
        substrings.push_back(substring);
    }
    return substrings;
}

/**
 * Set model parameters (SCF, BKG, etc.) according to a JSON config
 *
 * @param pdf The model whose parameters are to be changed
 * @param prefix String that precedes all model parameter names (e.g. channel name)
 * @param model_parameters JSON config with model parameter values
 */
void ChangeModelParameters(RooAbsPdf* pdf, const nlohmann::json& model_parameters,
                           const std::string prefix) {
    Log::LogLine(Log::info) << "Updating parameter values for channel " << prefix;
    for (auto& parameter : model_parameters.items()) {
        const char* name = parameter.key().c_str();
        const double value = parameter.value().get<double>();

        RooArgSet* pdf_parameters = pdf->getVariables();

        std::string full_name = prefix;
        full_name += name;
        RooRealVar* rooparameter =
            dynamic_cast<RooRealVar*>(pdf_parameters->find(full_name.c_str()));

        if (rooparameter) {
            Log::print(Log::debug, "Changing parameter %s to %f\n", name, value);
            rooparameter->setVal(value);
        } else {
            Log::print(Log::warning, "Parameter %s not found in model, skipping...\n", name);
        }
    }
}

/**
 * Format parameters of the supplied models into a JSON formatted string
 *
 * @param name Name of the parameter section; to be printed in the header
 * @param models Models whose parameters are to be printed
 * @param observables Observables that should be removed from the list of parameters
 *
 * @return std::string JSON formatted string
 */
std::string FormatResultsJSON(std::vector<const RooAbsPdf*> models,
                              const RooArgSet& observables) {
    std::string formatted_result = "{\n";
    const std::size_t buf_size = 1024;
    char buffer[buf_size];
    for (auto& model : models) {
        RooArgSet* vars = model->getVariables();
        vars->remove(observables);
        std::vector<RooRealVar*> vector_of_vars = ToVector<RooRealVar*>(*vars);
        for (auto& var : vector_of_vars) {
            int ret = snprintf(buffer, buf_size, "  \"%s\": %.4f,\n", var->GetName(), var->getVal());
            assert(ret >= 0);
            formatted_result.append(buffer);
        }
    }
    formatted_result.pop_back();
    formatted_result.pop_back();
    formatted_result += "\n}\n";
    return formatted_result;
}

/**
 * Format parameters of the supplied model into a JSON formatted string
 *
 * @param name Name of the parameter section; to be printed in the header
 * @param model Model whose parameters are to be printed
 * @param observables Observables that should be removed from the list of parameters
 *
 * @return std::string JSON formatted string
 */
std::string FormatResultsJSON(const RooAbsPdf* model, const RooArgSet& observables) {
    std::vector<const RooAbsPdf*> models = {model};
    return FormatResultsJSON(models, observables);
}

/**
 * Create a directory structure if it doesn't exist
 *
 * Normally a file couldn't be created if a directory where it's supposed to
 * reside doesn't exist. This function creates the necessary directory
 * structure if it doesn't exist. It does nothing if it exists.
 *
 * @param file Path to a file that might need a directory structure
 */
void CreateDirsIfNecessary(const std::string file) {
    boost::filesystem::path path(file);
    path.remove_filename();
    if (!boost::filesystem::is_directory(path)) {
        Log::LogLine(Log::debug) << "Directory " << path << " doesn't exist; creating";
        boost::filesystem::create_directories(path);
    }
}

/**
 * Rearrange a correlation histogram based on a vector of labels
 *
 * This function is quite finnicky with the indices, because of ROOT's peculiar
 * histogram indexing convention.
 *
 * @param matrix Original correlation histogram
 * @param ordered_labels The requested ordering
 *
 * @return TH2* An ordered copy of the original
 */
TH2* ArrangeCorrelationMatrix(const TH2* matrix, std::vector<std::string> ordered_labels) {
    TH2* ordered_corr_matrix = dynamic_cast<TH2*>(matrix->Clone("ordered_corr_matrix"));
    ordered_corr_matrix->Reset();
    const uint num_bins = matrix->GetNbinsX();
    assert(ordered_labels.size() == num_bins);

    // First we get the "migration" address book (since X axis is labeled from
    // bins 1 to N, while Y is labeled from N to 1, we have to do it separately)
    int new_addr_x[ordered_labels.size()];
    int new_addr_y[ordered_labels.size()];
    for (uint old_bin = 0; old_bin < num_bins; old_bin++) {
        std::string label_x = matrix->GetXaxis()->GetBinLabel(old_bin + 1);
        std::string label_y = matrix->GetXaxis()->GetBinLabel(num_bins - old_bin);
        for (uint new_bin = 0; new_bin < ordered_labels.size(); new_bin++) {
            // Switch the labels themselves
            ordered_corr_matrix->GetXaxis()->SetBinLabel(new_bin + 1, ordered_labels[new_bin].c_str());
            ordered_corr_matrix->GetYaxis()->SetBinLabel(num_bins - new_bin, ordered_labels[new_bin].c_str());

            // Fill in the address book
            if (ordered_labels[new_bin] == label_x) {
                new_addr_x[old_bin] = new_bin + 1;
            }
            if (ordered_labels[new_bin] == label_y) {
                new_addr_y[old_bin] = new_bin + 1;
            }
        }
    }

    // Migrate the contents, keeping in mind the different counting for Y axis
    for(uint x = 1; x <= num_bins; x++) {
        for(uint y = 1; y <= num_bins; y++) {
            int old_bin = matrix->GetBin(x, y);
            double content = matrix->GetBinContent(old_bin);
            int new_bin = ordered_corr_matrix->GetBin(new_addr_x[x - 1], num_bins - new_addr_y[y - 1] + 1);
            ordered_corr_matrix->SetBinContent(new_bin, content);
        }
    }
    return ordered_corr_matrix;
}

/**
 * Create and save a plot of a correlation matrix with ordered labels
 *
 * @param matrix Original correlation histogram
 * @param ordered_labels The requested ordering
 */
void PlotCorrelationMatrix(const RooFitResult& result, std::vector<std::string> ordered_labels) {
    TCanvas correlation_plot("correlation_plot", "correlation_plot", 500, 500);
    gStyle->SetPalette(kLightTemperature);
    gStyle->SetPaintTextFormat("1.1f");

    TH2* corr_histo = tools::ArrangeCorrelationMatrix(result.correlationHist(), ordered_labels);
    corr_histo->SetTitle("");
    corr_histo->SetMaximum(1);
    corr_histo->SetMinimum(-1);
    corr_histo->SetMarkerSize(1.1);
    corr_histo->LabelsOption("v");
    corr_histo->Draw("col text");
    correlation_plot.SaveAs(constants::format);
    gStyle->SetPalette(kViridis);
}

}  // namespace tools
