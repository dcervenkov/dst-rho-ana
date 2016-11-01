/**
 *  @file    tools.cc
 *  @author  Daniel Cervenkov, cervenkov@ipnp.mff.cuni.cz
 *  @date    2015-10-07
 *
 *  @brief A collections of various tools and utilities
 *
 */

#include "tools.h"

// ROOT includes
#include "TList.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "RooRealVar.h"

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
	TList *files = system_dir.GetListOfFiles();
	if (files) {
		TSystemFile *file;
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
	std::vector<TString> list_of_files = GetListOfFiles(dir,".root");
	TString directory = dir;
	// Make the function agnostic to possible '/' at the end of the dir argument
	if (!directory.EndsWith("/")) {
		directory += "/";
	}
	for (auto file_name : list_of_files) {
		if (file_name.Contains("results")) continue; // Skip results files created by DSRhoPeek
		chain->Add(directory + file_name);
	}
	return chain;
}

/**
 * Setup sane ROOT plot style
 */
void SetupPlotStyle() {
	// Fonts ending in 3 are not proportional in size to the pad; their size is in pixels
	// The "xyz" parameter dictates that these settings should be used for all axes
	gStyle->SetLabelFont(43,"xyz");
	gStyle->SetLabelSize(18,"xyz");
	gStyle->SetLabelOffset(0.01,"xyz");
	gStyle->SetTitleFont(43,"xyz");
	gStyle->SetTitleSize(18,"xyz");
	gStyle->SetTitleOffset(1.2);
	gStyle->SetMarkerSize(0.5);
	gStyle->SetEndErrorSize(0); // Disable perpendicular lines at the end of error bars
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetPadBottomMargin(0.1);
}

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

    TPaveText *stat_box = new TPaveText(x_left, y_bottom, x_right, y_top, "NDC");
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

} // namespace tools

