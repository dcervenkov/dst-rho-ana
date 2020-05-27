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
#include "RooArgSet.h"
#include "RooHistPdf.h"
#include "RooRealVar.h"
#include "TChain.h"
#include "TH2.h"
#include "TPaveText.h"

// Local includes
#include "nlohmann/json.hpp"

namespace tools {

std::vector<TString> GetListOfFiles(const char* dir, const char* ext);
TChain* ReadDataFromDir(const char* dir);
void SetupPlotStyle();
void SetPlotDir(const char* plot_dir);
TPaveText* CreateStatBox(double chi2, int ndof, const RooArgList& results,
                         bool position_top = true, bool position_left = true);
void PlotVar(const RooRealVar& var, const RooAbsData& data);
void PlotVar(const RooRealVar& var, const RooAbsPdf&);
void PlotVar(const RooRealVar& var, const RooDataHist& data1, const RooDataHist& data2,
             bool draw_pull, bool draw_residual);
void PlotWithPull(const RooRealVar& var, const RooArgSet& projection_vars, const RooDataSet& data,
                  const RooAbsPdf& pdf, const RooArgList& results,
                  const std::vector<RooAbsPdf*> components = std::vector<RooAbsPdf*>(),
                  int numCPUs = 1, std::string prefix = "", std::string title = "",
                  std::vector<RooCmdArg> options = std::vector<RooCmdArg>());
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
RooArgList BlindResults(const RooArgList& orig_results, double range);


void LogEnvironmentMetadata(TFile* file);
void LogCLIArguments(TFile* file, int argc, char* argv[]);
void LogTextFromFile(TFile* file, const char* field_name, const char* filename);
void LogFileCRC(TFile* file, const char* field_name, const char* filename);
void LogText(TFile* file, const char* field_name, const char* text);
void LogText(TFile* file, const char* field_name, const std::string text);
std::vector<std::string> SplitString(const std::string& input_string, char delimiter);
void ChangeModelParameters(RooAbsPdf* pdf, const nlohmann::json& config,
                           const std::string prefix = "");
std::string FormatResultsJSON(std::vector<const RooAbsPdf*> models, const RooArgSet& observables);
std::string FormatResultsJSON(const RooAbsPdf* model, const RooArgSet& observables);
nlohmann::json GetResultsJSON(const RooAbsPdf* model, const RooArgSet& observables,
                              std::string prefix);
nlohmann::json GetResultsJSON(std::vector<const RooAbsPdf*> models, const RooArgSet& observables,
                              std::string prefix);
void CreateDirsIfNecessary(const std::string file);
RooLinkedList VecToCmdList(std::vector<RooCmdArg>& commands);
TH2* ArrangeCorrelationMatrix(const TH2* matrix, std::vector<std::string> ordered_labels);
void PlotCorrelationMatrix(const RooFitResult& result, std::vector<std::string> ordered_labels);

nlohmann::json MergeJSON (const nlohmann::json& json1, const nlohmann::json& json2);
double RoundToDecimals(double number, int decimals);
RooHistPdf* CreatePdfFromHistos(const char* name, const char* title,
                                std::vector<RooDataHist*> histos, RooArgSet observables);

/**
 * Return a concatenation of two std::vectors
 */
template <class T>
std::vector<T> AddVectors(const std::vector<T>& vector1, const std::vector<T>& vector2) {
    std::vector<T> new_vector(vector1);
    new_vector.reserve(vector1.size() + vector2.size());
    new_vector.insert(new_vector.end(), vector2.begin(), vector2.end());
    return new_vector;
}

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
