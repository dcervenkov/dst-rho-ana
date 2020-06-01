/**
 *  @file    fitterdelta.cc
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2020-05-31
 *
 *  @brief Class that performs the fit itself as well as plotting
 *
 */

#include "fitterdelta.h"

// Standard includes
#include <array>
#include <boost/filesystem.hpp>
#include <sstream>
#include <string>

// ROOT includes
#include "RooArgSet.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "RooHist.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooTreeDataStore.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TEnv.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TTree.h"

// Local includes
#include "constants.h"
#include "log.h"
#include "tools.h"

FitterDelta::FitterDelta(nlohmann::json config) {
    TTree* signal_data = ReadInFile(config["inputFiles"]["signal"]);
    TTree* sidebands_data = ReadInFile(config["inputFiles"]["sidebands"]);

    const double min = config["fitRanges"]["thetab"]["min"].get<double>();
    const double max = config["fitRanges"]["thetab"]["max"].get<double>();
    TH1D* histo = new TH1D("histo", "histo", 20, min, max);

    FillSubtractionHisto(histo, "thetab", signal_data, sidebands_data);

    histo->Fit("chebyshev3");
    fit_function_ = histo->GetFunction("chebyshev3");

    TCanvas canvas("canvas", "canvas", 500, 500);
    histo->Draw("e1");
    canvas.SaveAs("test.png");
}

FitterDelta::~FitterDelta() {
    if (output_file_) {
        output_file_->Close();
    }
}

/**
 * Reads in data from ROOT file(s).
 *
 * @param file_names Vector of paths to the ROOT files
 * @param num_events [optional] Maximum number of events to use (0 to read all)
 */
TTree* FitterDelta::ReadInFile(const nlohmann::json data_files) const {
    TChain* chain = new TChain("h2000");
    int num_files = 0;
    for (auto& file_name : data_files.items()) {
        if (!chain->Add(file_name.value().get<std::string>().c_str(), 0)) {
            Log::LogLine(Log::error) << "File " << file_name.value() << " not found!";
            exit(9);
        }
        num_files++;
    }

    Log::print(Log::info, "Reading %i input files\n", num_files);

    TString common_cuts = tools::GetCommonCutsString();
    common_cuts += "&&evmcflag!=1";

    TTree* tree = nullptr;
    tree = chain->CopyTree(common_cuts);

    // Set the current directory back to the one for writing (ugly ROOT stuff)
    if (output_file_) {
        output_file_->cd();
    }

    Log::print(Log::info, "Input dataset ready\n");
    return tree;
}

/**
 * @brief Fill the supplied histo with tree1 - tree2
 *
 * Stupid ROOT doesn't permit the TTrees to be const as Draw() is not a const function.
 *
 * @param histo Histo to be filled
 * @param branch Name of the tree branch to be used
 * @param tree1 First tree
 * @param tree2 Second tree
 */
void FitterDelta::FillSubtractionHisto(TH1D* histo, TString branch, TTree* tree1, TTree* tree2) const {
    TString name = histo->GetName();
    tree1->Draw(branch + ">>" + name, "", "goff");
    TH1D* temp_histo = static_cast<TH1D*>(histo->Clone("temp_histo"));
    temp_histo->Reset();
    tree2->Draw(branch + ">>temp_histo", "", "goff");
    TCanvas canvas("canvas", "canvas", 500, 500);
    temp_histo->SetLineColor(kRed);
    temp_histo->Draw();
    histo->Draw("same");
    canvas.SaveAs("sig_vs_side.png");
    histo->Add(temp_histo, -1);
}

nlohmann::json FitterDelta::GetJSONResults() const {
    nlohmann::json results = tools::GetResultsJSON(*fit_function_, "thetab_corr_");
    // RooFit's Chebyshev doesn't use the constant term
    results.erase("thetab_corr_p0");
    return results;
};
