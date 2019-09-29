/**
 *  @file    tools.cc
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2015-10-07
 *
 *  @brief A collections of various tools and utilities
 *
 */

// Standard includes
#include <algorithm>
#include <fstream>
#include <sstream>

// Boost includes
#include <boost/filesystem.hpp>

// ROOT includes
#include "RooAbsData.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TEnv.h"
#include "TFile.h"
#include "TH2D.h"
#include "TList.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"

// Local includes
#include "cksum.h"
#include "constants.h"
#include "gitversion.h"
#include "log.h"

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
TPaveText* CreateStatBox(double chi2, RooArgList* results, bool position_top, bool position_left) {
    double x_left, x_right, y_bottom, y_top;
    const double line_height = 0.06;

    if (!results) {
        results = new RooArgList();
    }

    if (position_top) {
        y_top = 0.9;
        y_bottom = y_top - results->getSize() * line_height;
    } else {
        y_bottom = 0.05;
        y_top = y_bottom + results->getSize() * line_height;
    }

    if (position_left) {
        x_left = 0.3;
        x_right = 0.301;
    } else {
        x_left = 0.7;
        x_right = 0.701;
    }

    TPaveText* stat_box = new TPaveText(x_left, y_bottom, x_right, y_top, "NDC");
    stat_box->SetShadowColor(kWhite);
    stat_box->SetBorderSize(0);
    stat_box->SetFillColor(kWhite);
    stat_box->SetTextFont(43);
    stat_box->SetTextSize(14);
    stat_box->SetY1NDC(0.1);

    char line[1000];
    for (int i = 0; i < results->getSize(); i++) {
        snprintf(line, 1000, "%s = %.3f +- %.3f", results[i].GetTitle(),
                 static_cast<RooRealVar*>(results->at(i))->getVal(),
                 static_cast<RooRealVar*>(results->at(i))->getError());
        stat_box->AddText(line);
    }
    snprintf(line, 1000, "#chi^{2} = %.2f\n", chi2);
    stat_box->AddText(line);
    return stat_box;
}

/*
 * Create and save 2D plots of supplied vars and data
 *
 * @param var1 First variable.
 * @param var2 Second variable.
 * @param data Data to be plotted.
 * @param format [optional] Format in which to save the images.
 * @param max [optional] Maximum of the range
 */
void PlotVars2D(const RooRealVar& var1, const RooRealVar& var2, const RooAbsData& data,
                const char* format, const double max) {
    TCanvas canvas(
        TString(data.GetName()) + "_" + TString(var1.GetName()) + "_" + TString(var2.GetName()),
        TString(data.GetName()) + "_" + TString(var1.GetName()) + "_" + TString(var2.GetName()),
        500, 500);

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
                const RooAbsData& pdf, const char* format, const bool residual) {
    TString type = residual ? "_residual" : "_pull";
    TCanvas canvas(TString(data.GetName()) + "_" + TString(var1.GetName()) + "_" +
                       TString(var2.GetName()) + type,
                   TString(data.GetName()) + "_" + TString(var1.GetName()) + "_" +
                       TString(var2.GetName()) + type,
                   500, 500);
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
                const RooAbsData& data2, const char* format) {
    double max = 0;
    TH2D* histo1 = static_cast<TH2D*>(data1.createHistogram("histo1", var1, RooFit::YVar(var2)));
    TH2D* histo2 = static_cast<TH2D*>(data2.createHistogram("histo2", var1, RooFit::YVar(var2)));
    max = std::max(histo1->GetMaximum(), histo2->GetMaximum());
    delete histo1;
    delete histo2;
    PlotVars2D(var1, var2, data1, format, max);
    PlotVars2D(var1, var2, data2, format, max);
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
    std::ofstream file(filename);
    file << text;
    file.close();
}

/**
 * Transform RooArgSet of RooRealVars into a std::vector
 */
std::vector<RooRealVar*> ToVector(const RooArgSet& set) {
    std::vector<RooRealVar*> vector;
    TIterator* iterator = set.createIterator();
    RooRealVar* arg;
    while ((arg = dynamic_cast<RooRealVar*>(iterator->Next()))) {
        vector.push_back(arg);
    }
    delete iterator;
    return vector;
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

}  // namespace tools
