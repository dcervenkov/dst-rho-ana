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
#include "RooHistPdf.h"
#include "RooArgSet.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "RooHist.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooTreeDataStore.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TEnv.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
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
    signal_data->Print();
    sidebands_data->Print();
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
