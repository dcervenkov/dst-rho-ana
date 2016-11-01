/**
 *  @file    tools.h
 *  @author  Daniel Cervenkov, cervenkov@ipnp.mff.cuni.cz
 *  @date    2015-10-07
 *
 *  @brief A collections of various tools and utilities
 *
 */

#ifndef TOOLS_H_
#define TOOLS_H_

// ROOT includes
#include "TChain.h"
#include "TPaveText.h"
#include "RooArgList.h"

namespace tools {

std::vector<TString> GetListOfFiles(const char* dir, const char* ext);
TChain* ReadDataFromDir(const char* dir);
void SetupPlotStyle();
TPaveText* CreateStatBox(double chi2, RooArgList* results = NULL, bool position_top = true, bool position_left = true);

}



#endif /* TOOLS_H_ */
