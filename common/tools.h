/**
 *  @file    tools.h
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2015-10-07
 *
 *  @brief A collections of various tools and utilities
 *
 */

#pragma once

// ROOT includes
#include "RooAbsData.h"
#include "RooArgList.h"
#include "RooRealVar.h"
#include "TChain.h"
#include "TPaveText.h"

// Local includes
#include "nlohmann/json.hpp"

namespace tools {

std::vector<TString> GetListOfFiles(const char* dir, const char* ext);
TChain* ReadDataFromDir(const char* dir);
void SetupPlotStyle();
void SetPlotDir(const char* plot_dir);
TPaveText* CreateStatBox(double chi2, RooArgList* results = nullptr, bool position_top = true,
                         bool position_left = true);
void PlotVars2D(const RooRealVar& var1, const RooRealVar& var2, const RooAbsData& data1,
                const RooAbsData& data2, const std::string prefix = "",
                const char* format = ".pdf");
void PlotVars2D(const RooRealVar& var1, const RooRealVar& var2, const RooAbsData& data,
                const std::string prefix = "", const char* format = ".pdf", const double max = 0);
void PlotPull2D(const RooRealVar& var1, const RooRealVar& var2, const RooAbsData& data,
                const RooAbsData& pdf, const std::string prefix = "", const char* format = ".pdf",
                const bool residual = false);
TString GetCommonCutsString();
void SaveTextToFile(const std::string filename, const std::string text);


void LogEnvironmentMetadata(TFile* file);
void LogCLIArguments(TFile* file, int argc, char* argv[]);
void LogTextFromFile(TFile* file, const char* field_name, const char* filename);
void LogFileCRC(TFile* file, const char* field_name, const char* filename);
void LogText(TFile* file, const char* field_name, const char* text);
void LogText(TFile* file, const char* field_name, const std::string text);
std::vector<std::string> SplitString(const std::string& input_string, char delimiter);
void ChangeModelParameters(RooAbsPdf* pdf, const nlohmann::json& config,
                           const std::string prefix = "");

/**
 * Transform RooArgSet of RooRealVars into a std::vector
 * 
 * Objects other than RooRealVars (e.g., RooCategory are not added to the final
 * vector).
 */
template <class T>
std::vector<T> ToVector(const RooAbsCollection& set) {
    std::vector<T> vector;
    TIterator* iterator = set.createIterator();
    T arg;
    TObject* obj;
    while ((obj = iterator->Next())) {
        if ((arg = dynamic_cast<T>(obj))) {
            vector.push_back(arg);
        }
    }
    delete iterator;
    return vector;
}

}  // namespace tools
