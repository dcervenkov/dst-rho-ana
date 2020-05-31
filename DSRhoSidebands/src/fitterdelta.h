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
#include "RooAddPdf.h"
#include "RooBreitWigner.h"
#include "RooBifurGauss.h"
#include "RooChebychev.h"
#include "RooExponential.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooPolynomial.h"
#include "RooProdPdf.h"
#include "RooVoigtian.h"
#include "RooRealVar.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TPaveText.h"

// Local includes
#include "constants.h"
#include "nlohmann/json.hpp"

class FitterDelta {
   public:
    FitterDelta(nlohmann::json config);
    virtual ~FitterDelta();

    TTree* ReadInFile(const nlohmann::json data_files) const;
    void SetPlotDir(const char* output_dir);
    void Fit(RooAbsPdf* pdf, RooDataSet* data);
    void PrintResultsJSON() const;

   private:
    TFile* output_file_ = nullptr;
};
