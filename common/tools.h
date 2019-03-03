/**
 *  @file    tools.h
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2015-10-07
 *
 *  @brief A collections of various tools and utilities
 *
 */

#ifndef TOOLS_H_
#define TOOLS_H_

// ROOT includes
#include "RooArgList.h"
#include "RooAbsData.h"
#include "RooRealVar.h"
#include "TChain.h"
#include "TPaveText.h"

namespace tools {

std::vector<TString> GetListOfFiles(const char* dir, const char* ext);
TChain* ReadDataFromDir(const char* dir);
void SetupPlotStyle();
void SetPlotDir(const char* plot_dir);
TPaveText* CreateStatBox(double chi2, RooArgList* results = NULL, bool position_top = true,
                         bool position_left = true);
void PlotVars2D(const RooRealVar& var1, const RooRealVar& var2, const RooAbsData& data1,
                const RooAbsData& data2, const char* format = ".pdf");
void PlotVars2D(const RooRealVar& var1, const RooRealVar& var2, const RooAbsData& data,
                const char* format = ".pdf", const double max = 0);
void PlotPull2D(const RooRealVar& var1, const RooRealVar& var2, const RooAbsData& data,
                const RooAbsData& pdf, const char* format = ".pdf", const bool residual = false);
}  // namespace tools

#endif /* TOOLS_H_ */
