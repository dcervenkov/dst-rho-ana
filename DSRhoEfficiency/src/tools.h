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
#include "TChain.h"

namespace tools {

std::vector<TString> GetListOfFiles(const char* dir, const char* ext);
TChain* ReadDataFromDir(const char* dir);
void SetupPlotStyle();

}



#endif /* TOOLS_H_ */
