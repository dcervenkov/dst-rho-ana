/**
 *  @file    tools.cc
 *  @author  Daniel Cervenkov, cervenkov@ipnp.mff.cuni.cz
 *  @date    2015-10-07
 *
 *  @brief A collections of various tools and utilities
 *
 */

// Standard includes
#include <stdio.h>

// ROOT includes
#include "TChain.h"
#include "TList.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TStyle.h"

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
	gStyle->SetPadLeftMargin(0.105);
	gStyle->SetPadBottomMargin(0.1);
	gStyle->SetOptStat(0);

}

} // namespace tools

