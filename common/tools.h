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
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "TChain.h"
#include "TPaveText.h"

namespace tools {

std::vector<TString> GetListOfFiles(const char* dir, const char* ext);
TChain* ReadDataFromDir(const char* dir);
void SetupPlotStyle();
TPaveText* CreateStatBox(double chi2, RooArgList* results = NULL, bool position_top = true,
                         bool position_left = true);
void PlotVars2D(const RooRealVar& var1, const RooRealVar& var2, const RooDataHist& data,
                const char* format = ".pdf");
void PlotPull2D(const RooRealVar& var1, const RooRealVar& var2, const RooDataHist& data,
                const RooDataHist& pdf, const char* format = ".pdf");
}  // namespace tools

#endif /* TOOLS_H_ */
