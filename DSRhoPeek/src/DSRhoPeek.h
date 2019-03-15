/*
 * DSRhoPeek.h
 *
 *  Created on: Aug 25, 2014
 *      Author: cervenkov
 */

#ifndef DSRHOPEEK_H_
#define DSRHOPEEK_H_

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

TTree* getBASFTree(TFile* file, const int treeNumber);
unsigned int getHistogramFromTree(TH1F* histogram, TTree* tree, const char* branchName, const int flag, const int candidateSelection = -1);
void getOptimalHistogramScope(TTree* tree, const char* branchName, double& min, double& max);

#endif /* DSRHOPEEK_H_ */
