/**
 *  @file    fitterdelta.h
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2020-05-31
 *
 *  @brief This class performs the fitting itself as well as plotting
 *
 */

#pragma once

// Standard includes
#include <array>

// ROOT includes
#include "TCanvas.h"
#include "TChain.h"
#include "TH1D.h"
#include "TString.h"

// Local includes
#include "constants.h"
#include "nlohmann/json.hpp"
#include "tools.h"

class FitterDelta {
   public:
    FitterDelta(nlohmann::json config);
    virtual ~FitterDelta();

    TTree* ReadInFile(const nlohmann::json data_files) const;
    void SetPlotDir(const char* output_dir);
    void PrintResultsJSON() const;
    void SetOutputFile(TFile* file) { output_file_ = file; };
    nlohmann::json GetJSONResults() const;

   private:
    void FillSubtractionHisto(TH1D* histo, TString branch, TTree* tree1, TTree* tree2) const;
    void FillRatioHisto(TH1D* histo, TString branch, TTree* tree1, TTree* tree2) const;

    TFile* output_file_ = nullptr;
    TF1* fit_function_ = nullptr;
};
